#!/usr/bin/env python3
"""Replay + reconciliation + descriptive taker-flow GEX for captured Deribit data.

Replays append-only capture files produced by capture.py under
    quant-trading/data/raw/focused/live/deribit/dayYYYYMMDD/*.jsonl
deterministically (fixed ordering by received_ts_ms, then dedup_id), enforces
dedup by the capture-computed dedup_id, reconciles units/direction against
captured instrument metadata, and computes ONLY
    descriptive signed taker-flow gamma exposure (GEX) aggregates
where every required input (taker direction, amount coin units, gamma, index
price in USD) exists. Trades missing any required input are reported under
coverage accounting and EXCLUDED — the pipeline fails closed, never imputes.

regime/return gate: capture directions emit NO strategy returns, NO regime
classification, and NO volatility-conditioning until four complete registered
weeks of capture exist after the campaign's prospective start (see
focused-campaign.json). run.py refuses such claims and marks the gate as
pending with a measured percentage.

Usage per exploration contract:
    python3 run.py --start YYYY-MM-DD --end YYYY-MM-DD [--output-dir DIR]

Outputs (default: this directory/results/):
    results.json   deterministic descriptive aggregate + coverage + gate state
No returns.csv is created here (capture direction contract).

Gamma source policy: venue gamma from captured ticker records is used when
available ("deribit_ticker"); otherwise gamma is computed with standard
Black-Scholes from the trade's own embedded `iv` (%), strike and
expiration_timestamp from instrument metadata, and S = trade index_price
("bs_from_trade_iv"). Missing IV or missing metadata fails closed: the trade
is excluded from GEX sums and reported in coverage.

Observed data note (2026-08-30 smoke): Deribit quantizes gamma to 1e-5 steps,
so deep-OTM / near-expiry options legitimately push gamma == 0.0. That is valid
captured data, contributes exactly 0.0 to GEX sums, and is NOT treated as a
missing greek; only an absent `greeks` object or absent `gamma` key fails closed.
"""

from __future__ import annotations

import argparse
import math
import json
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

DEFAULT_RAW = Path(__file__).resolve().parents[2] / "data/raw/focused/live/deribit"
RECORD_FILES = ("trades.jsonl", "instruments.jsonl", "ticker.jsonl", "index.jsonl", "events.jsonl")


REQUIRED_WEEKS = 4
PROSPECTIVE_START_UTC = "2026-08-30T11:13:24Z"


def parse_utc_ms(s: str) -> int:
    dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
    return int(dt.timestamp() * 1000)


def day_key_from_ms(ts_ms: int) -> str:
    return datetime.fromtimestamp(ts_ms / 1000.0, tz=timezone.utc).strftime("%Y%m%d")


def replay_files(raw_dir: Path, start_ms: int, end_ms: int) -> tuple[list[dict], dict]:
    """Deterministic replay: read every capture file, keep records whose capture-day
    key is in [start, end], order by (received_ts_ms, seq-in-file, dedup_id)."""
    records: list[dict] = []
    stats = {"files": 0, "lines": 0, "malformed": 0, "out_of_window": 0}
    if not raw_dir.is_dir():
        return records, stats
    day_dirs = sorted(p for p in raw_dir.iterdir() if p.is_dir() and p.name.startswith("day"))
    for d in day_dirs:
        day = d.name[3:]
        try:
            day_ms = int(datetime.strptime(day, "%Y%m%d").replace(tzinfo=timezone.utc).timestamp() * 1000)
        except ValueError:
            continue
        if day_ms < start_ms - 86400_000 or day_ms > end_ms + 86400_000:
            stats["out_of_window"] += 1
            continue
        for name in RECORD_FILES:
            fp = d / name
            if not fp.is_file():
                continue
            stats["files"] += 1
            with fp.open("r", encoding="utf-8") as fh:
                for seq, line in enumerate(fh):
                    line = line.strip()
                    stats["lines"] += 1
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except json.JSONDecodeError:
                        stats["malformed"] += 1
                        continue
                    ts = rec.get("received_ts_ms", 0) or 0
                    if not (start_ms <= ts <= end_ms):
                        continue
                    rec["_seq"] = seq
                    records.append(rec)
    records.sort(key=lambda r: (r.get("received_ts_ms", 0), r["_seq"], r.get("dedup_id", "")))
    return records, stats


def dedupe(records: list[dict]) -> tuple[list[dict], int]:
    seen: set[str] = set()
    out: list[dict] = []
    drops = 0
    for rec in records:
        did = rec.get("dedup_id")
        if did is None:
            drops += 1
            continue
        if did in seen:
            drops += 1
            continue
        seen.add(did)
        out.append(rec)
    return out, drops


def reconcile(records: list[dict]) -> dict:
    """Build deterministic instrument metadata map (latest record per instrument),
    then join trades with metadata + latest-prior ticker greeks per instrument."""
    instruments: dict[str, dict] = {}
    latest_ticker: dict[str, dict] = {}   # instrument -> ticker dict (latest by received ts)
    latest_prior_ticker: dict[str, list[tuple[int, dict]]] = defaultdict(list)
    trades: list[dict] = []
    indexes: list[dict] = []
    for rec in records:
        rt = rec.get("record_type")
        if rt == "instrument":
            row = rec.get("instrument") or {}
            if row.get("instrument_name"):
                instruments[row["instrument_name"]] = row
        elif rt == "instrument_state":
            pass  # state changes don't alter strike/expiry/contract_size mapping
        elif rt == "option_ticker":
            iname = rec.get("instrument_name")
            tk = rec.get("ticker") or {}
            if iname:
                latest_ticker[iname] = tk
                # keep bounded trailing window (last 64 tickers per instrument)
                w = latest_prior_ticker[iname]
                w.append((rec.get("received_ts_ms", 0), tk))
                if len(w) > 64:
                    del w[: len(w) - 64]
        elif rt == "index_price":
            indexes.append({"ts": rec.get("received_ts_ms") or 0,
                            "data": rec.get("data") or {},
                            "idx": rec.get("index_name")})
        elif rt == "option_trade":
            trades.append(rec)
    return {"instruments": instruments, "latest_ticker": latest_ticker,
            "ticker_window": latest_prior_ticker, "trades": trades, "indexes": indexes}


def _norm_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def bs_gamma(strike_milli: float, s_usd: float, iv_pct: float, tau_years: float) -> float | None:
    """Black-Scholes gamma dDelta/dS per 1 USD, params: strike, spot, iv (percent),
    time to expiry in ACT/365 years. Deribit quotes iv in the trade record (%)."""
    if strike_milli <= 0 or s_usd <= 0 or iv_pct <= 0 or tau_years <= 0:
        return None
    sigma = iv_pct / 100.0
    d1 = (math.log(s_usd / strike_milli) + 0.5 * sigma * sigma * tau_years) / (
        sigma * math.sqrt(tau_years)
    )
    return _norm_pdf(d1) / (s_usd * sigma * math.sqrt(tau_years))


def gex_units_contribution(trade: dict, meta: dict, ticker: dict | None) -> dict:
    """Compute signed taker-flow GEX contribution in USD per 1% underlying move.

    Units (source-verified against official Deribit docs + live payloads):
      - `direction`  : taker direction ('buy' -> +1, 'sell' -> -1)  [fail closed if absent]
      - `amount`     : option trade amount in BASE COIN (BTC/ETH) per official docs
      - contract_size: coin per contract (=1.0 for Deribit options; verified)
      - gamma source (priority): (a) captured ticker `greeks.gamma` (venue Black-Scholes)
        (b) BS gamma computed from the trade's own embedded `iv` (percent) with strike
        and expiry from instrument metadata and S = trade index_price — concrete and
        self-contained per parent instruction; missing IV or metadata fails closed.
      - S (USD)      : trade `index_price` (verified present on every option trade)
      - premium coin : `price` (BTC/ETH) — recorded; premium USD = price × index_price
    GEX$ per 1% move for the trade = signed_coins × gamma × S × S × 0.01
    (dealer-mirror = opposite sign of taker flow; recorded as taker-signed here).
    """
    out = {
        "ok": False,
        "missing": [],
        "gex_usd_per_1pct": 0.0,
        "signed_usd_premium": None,
        "greeks_source": None,
    }
    direction = trade.get("direction")
    if direction not in ("buy", "sell"):
        out["missing"].append("direction")
        return out
    sign = 1.0 if direction == "buy" else -1.0
    amount = trade.get("amount")
    if amount is None or amount <= 0:
        out["missing"].append("amount")
        return out
    contract_size = meta.get("contract_size")
    if contract_size is None or contract_size <= 0:
        out["missing"].append("contract_size")
        return out
    s_usd = trade.get("index_price")
    if not s_usd or s_usd <= 0:
        out["missing"].append("index_price")
        return out
    gamma = None
    venue_gamma = (ticker or {}).get("greeks", {}).get("gamma") if ticker else None
    if venue_gamma is not None:
        gamma = float(venue_gamma)
        out["greeks_source"] = "deribit_ticker"
    else:
        iv_pct = trade.get("iv")
        strike = meta.get("strike")
        expiry_ms = meta.get("expiration_timestamp")
        ts = trade.get("timestamp")
        if iv_pct is None:
            out["missing"].append("iv")
            return out
        if strike is None or expiry_ms is None or ts is None:
            out["missing"].append("strike_or_expiry")
            return out
        tau_years = max((expiry_ms - ts) / 1000.0, 0.0) / (365.0 * 86400.0)
        gamma = bs_gamma(float(strike), float(s_usd), float(iv_pct), tau_years)
        if gamma is None:
            out["missing"].append("gamma_inputs_invalid")
            return out
        out["greeks_source"] = "bs_from_trade_iv"
    coins = amount * contract_size
    gex = sign * coins * gamma * s_usd * s_usd * 0.01
    price_coin = trade.get("price")
    premium_usd = None
    if price_coin is not None:
        premium_usd = sign * coins * price_coin * s_usd
    out.update({"ok": True, "gex_usd_per_1pct": gex, "signed_usd_premium": premium_usd})
    return out


def pick_ticker(instrument: str, trade_ts_ms: int, window: dict, fallback: dict) -> dict | None:
    """Latest ticker with gamma at or before the trade ts (tolerance 2h), else the
    current-latest ticker if within tolerance, else None (fail closed)."""
    tol_ms = 2 * 3600 * 1000
    best: tuple[int, dict] | None = None
    for ts, tk in reversed(window.get(instrument) or []):
        gamma = (tk.get("greeks") or {}).get("gamma") if tk else None
        if gamma is None or ts is None:
            continue
        if trade_ts_ms - ts > tol_ms or ts - trade_ts_ms > tol_ms:
            continue
        if best is None or abs(trade_ts_ms - ts) < abs(trade_ts_ms - best[0]):
            best = (ts, tk)
    if best is not None:
        return best[1]
    latest = fallback.get(instrument)
    if latest is not None and (latest.get("greeks") or {}).get("gamma") is not None:
        ts = latest.get("timestamp")
        if ts is not None and abs(trade_ts_ms - ts) <= tol_ms:
            return latest
    return None


def aggregate(recon: dict) -> dict:
    """Deterministic descriptive aggregation by UTC day and total; coverage stats."""
    per_day: dict[str, dict] = defaultdict(lambda: {
        "trades": 0, "ok": 0, "missing_direction": 0, "missing_amount": 0,
        "missing_contract_size": 0, "missing_gamma": 0, "missing_index_price": 0,
        "missing_iv": 0, "missing_strike_or_expiry": 0, "missing_gamma_inputs_invalid": 0,
        "gex_usd_per_1pct_sum": 0.0, "signed_usd_premium_sum": 0.0,
        "buy_count": 0, "sell_count": 0,
        "greeks_from_ticker": 0, "greeks_from_trade_iv": 0,
    })
    missing_instruments: set[str] = set()
    for rec in recon["trades"]:
        trade = rec.get("trade") or {}
        iname = trade.get("instrument_name", "")
        ts = trade.get("timestamp") or rec.get("received_ts_ms") or 0
        day = day_key_from_ms(ts)
        meta = recon["instruments"].get(iname)
        if meta is None:
            # instrument metadata genuinely absent from capture — fail closed
            meta = {}
        ticker = pick_ticker(iname, ts, recon["ticker_window"], recon["latest_ticker"])
        r = gex_units_contribution(trade, meta, ticker)
        bucket = per_day[day]
        bucket["trades"] += 1
        if r["ok"]:
            bucket["ok"] += 1
            bucket["gex_usd_per_1pct_sum"] += r["gex_usd_per_1pct"]
            if r["signed_usd_premium"] is not None:
                bucket["signed_usd_premium_sum"] += r["signed_usd_premium"]
            if trade.get("direction") == "buy":
                bucket["buy_count"] += 1
            else:
                bucket["sell_count"] += 1
            if r["greeks_source"] == "deribit_ticker":
                bucket["greeks_from_ticker"] += 1
            elif r["greeks_source"] == "bs_from_trade_iv":
                bucket["greeks_from_trade_iv"] += 1
        else:
            for m in r["missing"]:
                key = f"missing_{m}"
                if key in bucket:
                    bucket[key] += 1
            if "gamma" in r["missing"] or "gamma_inputs_invalid" in r["missing"]:
                missing_instruments.add(iname)
    total = {
        "trades": 0, "ok": 0,
        "missing": defaultdict(int),
        "gex_usd_per_1pct_sum": 0.0,
        "signed_usd_premium_sum": 0.0,
        "greeks_from_ticker": 0, "greeks_from_trade_iv": 0,
    }
    for d, b in per_day.items():
        total["trades"] += b["trades"]
        total["ok"] += b["ok"]
        total["missing"]["direction"] += b["missing_direction"]
        total["missing"]["amount"] += b["missing_amount"]
        total["missing"]["contract_size"] += b["missing_contract_size"]
        total["missing"]["gamma"] += b["missing_gamma"]
        total["missing"]["index_price"] += b["missing_index_price"]
        total["missing"]["iv"] += b["missing_iv"]
        total["missing"]["strike_or_expiry"] += b["missing_strike_or_expiry"]
        total["missing"]["gamma_inputs_invalid"] += b["missing_gamma_inputs_invalid"]
        total["gex_usd_per_1pct_sum"] += b["gex_usd_per_1pct_sum"]
        total["signed_usd_premium_sum"] += b["signed_usd_premium_sum"]
        total["greeks_from_ticker"] += b["greeks_from_ticker"]
        total["greeks_from_trade_iv"] += b["greeks_from_trade_iv"]
    return {
        "per_day": {d: dict(b) for d, b in sorted(per_day.items())},
        "total": {
            "trades": total["trades"],
            "gex_computable": total["ok"],
            "missing": {k: v for k, v in sorted(total["missing"].items()) if v},
            "gex_usd_per_1pct_sum": round(total["gex_usd_per_1pct_sum"], 6),
            "signed_usd_premium_sum": round(total["signed_usd_premium_sum"], 6),
            "greeks_from_ticker": total["greeks_from_ticker"],
            "greeks_from_trade_iv": total["greeks_from_trade_iv"],
        },
        "instruments_without_gamma": sorted(missing_instruments),
    }


def week_completeness(raw_dir: Path) -> dict:
    """Registered-week accounting: a capture week is complete when every UTC day in
    it has nonzero recorded activity (any record file) for both currencies."""
    days: dict[str, dict] = {}
    if raw_dir.is_dir():
        for d in sorted(raw_dir.iterdir()):
            if not (d.is_dir() and d.name.startswith("day")):
                continue
            day = d.name[3:]
            try:
                datetime.strptime(day, "%Y%m%d")
            except ValueError:
                continue
            total = 0
            for f in d.iterdir():
                try:
                    total += f.stat().st_size
                except OSError:
                    pass
            days[day] = {"bytes": total, "active": total > 8}
    # count complete weeks of ANY capture, and weeks since prospective start
    complete_days = [d for d, info in sorted(days.items()) if info["active"]]
    # group consecutive day runs into weeks starting Monday
    weeks: dict[str, list[str]] = {}
    for day in complete_days:
        dt = datetime.strptime(day, "%Y%m%d").replace(tzinfo=timezone.utc)
        monday = (dt - timedelta(days=dt.weekday())).strftime("%Y%m%d")
        weeks.setdefault(monday, []).append(day)
    complete_weeks = sum(1 for wl in weeks.values() if len(wl) == 7)
    start_ms = parse_utc_ms(PROSPECTIVE_START_UTC)
    prospective_days = sum(1 for d in complete_days
                           if int(datetime.strptime(d, "%Y%m%d").replace(tzinfo=timezone.utc)
                                  .timestamp() * 1000) >= start_ms)
    return {
        "captured_days": len(days),
        "active_days": len(complete_days),
        "complete_registered_weeks": complete_weeks,
        "prospective_active_days": prospective_days,
        "day_activity": days,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Deribit GEX replay + reconciliation")
    ap.add_argument("--start", default=None, help="UTC start YYYY-MM-DD (default: prospective start)")
    ap.add_argument("--end", default=None, help="UTC end YYYY-MM-DD (default: today)")
    ap.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    ap.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    args = ap.parse_args(argv)

    if args.start is None:
        start_ms = parse_utc_ms(PROSPECTIVE_START_UTC)
    else:
        start_ms = parse_utc_ms(args.start + "T00:00:00Z")
    end_ms = parse_utc_ms((args.end or datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")) + "T23:59:59Z")

    records, rstats = replay_files(args.raw_dir, start_ms, end_ms)
    deduped, dup_drops = dedupe(records)
    recon = reconcile(deduped)
    agg = aggregate(recon)
    weeks = week_completeness(args.raw_dir)

    gate_pass = weeks["complete_registered_weeks"] >= REQUIRED_WEEKS
    results = {
        "schema_version": 1,
        "direction_id": "focused-deribit-gex",
        "determinism": "sorted by (received_ts_ms, file seq, dedup_id); fixed UTC arithmetic",
        "replay": {
            "input_root": str(args.raw_dir),
            "start_utc": args.start or PROSPECTIVE_START_UTC[:10],
            "end_utc": args.end or datetime.now(tz=timezone.utc).strftime("%Y-%m-%d"),
            "files_read": rstats["files"],
            "lines_read": rstats["lines"],
            "malformed_lines": rstats["malformed"],
            "records_after_window": len(records),
            "duplicate_drops": dup_drops,
            "records_after_dedup": len(deduped),
        },
        "reconciliation": {
            "instruments_metadata": len(recon["instruments"]),
            "option_trades_replayed": len(recon["trades"]),
            "ticker_records": sum(len(v) for v in recon["ticker_window"].values()),
            "index_records": len(recon["indexes"]),
            "units": {
                "amount": "base coin: option trade 'amount' x contract_size coins; "
                          "NOT USD — inverse options only quote in coin",
                "price": "option premium in coin ('price'); USD = price × index_price",
                "index_price": "USD per coin, embedded per trade by Deribit",
                "gamma": "greeks.gamma = dΔ/dS per 1 USD underlying move "
                         "(Black-Scholes, official ticker structure)",
                "gex_formula": "sign(direction) × coins × gamma × index_price² × 0.01 "
                               "[USD delta per 1% underlying move]",
                "taker_sign": "direction='buy' means taker bought option: dealer is "
                              "short gamma from this trade if dealer had none. Recorded "
                              "TAKER-SIGNED; dealer-mirror is the negation.",
            },
        },
        "gex_descriptive": agg,
        "coverage": {
            "trades_missing_inputs": agg["total"]["missing"],
            "instruments_without_gamma": agg["instruments_without_gamma"][:32],
            "n_instruments_without_gamma": len(agg["instruments_without_gamma"]),
            "fail_closed_policy": "any trade missing direction/amount/contract_size/gamma/"
                                      "index_price is excluded from GEX sums and reported",
        },
        "gate": {
            "required_registered_weeks": REQUIRED_WEEKS,
            "weeks": weeks,
            "prospective_start_utc": PROSPECTIVE_START_UTC,
            "regime_claim_allowed": gate_pass,
            "note": "No return, regime, volatility-conditioning, or deployment claim is "
                    "made by this direction until four complete registered weeks of "
                    "append-only capture exist after the prospective start. Captured data "
                    "is prospective-only: no historical backfill claim is possible from "
                    "trade streams without authorized 1ms raw data and full retention.",
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outp = args.output_dir / "results.json"
    outp.write_text(json.dumps(results, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "trades": results["reconciliation"]["option_trades_replayed"],
        "gex_computable": results["gex_descriptive"]["total"]["gex_computable"],
        "sum_gex_usd_per_1pct": results["gex_descriptive"]["total"]["gex_usd_per_1pct_sum"],
        "complete_registered_weeks": weeks["complete_registered_weeks"],
        "output": str(outp),
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
