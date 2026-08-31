#!/usr/bin/env python3
"""Frozen fast screen: deribit_intraday_taker_gex (accelerated-intraday-contract.json).

Implements the registered intraday dealer-mirror taker-GEX vs next-15-minute
realized-volatility screen with every frozen semantic from the contract:
  - immutable_input_prefix: stat byte size first, hash+parse exactly that prefix,
    ignore one trailing partial line, record path/bytes/SHA-256/first/latest
    exchange timestamps/snapshot wall clock
  - clock_semantics: exchange timestamps only (trade.timestamp, index
    data.timestamp); receipt timestamps are provenance only
  - dedup by the existing deterministic dedup_id, replay order
    (received_ts_ms, in-file seq, dedup_id) identical to run.py
  - GEX exactly per final_semantic_amendment.deribit.gex_formula: Black-Scholes
    gamma from the trade's own IV + captured instrument metadata; ticker gamma
    NEVER used (reuses run.py's bs_gamma implementation)
  - dealer-mirror aggregation by underlying x non-overlapping UTC 15-min block
  - next-block RV with anchor/final/>=720-obs/<=5000ms-gap completeness rules
  - earliest outcome-availability stop ordered by outcome-block end, bounded by
    the frozen 90-minute deadline, verified by a truncated-input rerun
  - registered sample conditions + gates with UNKNOWN truth totalization
  - JSON-infinity rule for the sign ratio (null numeric + "+infinity" string)

Research-only fast screen: reports fast-screen-pass / falsified / inconclusive.
No returns file, no durable-alpha claim, no deployment claim, no imputation.

Usage:
    python3 fast_intraday.py [--raw-dir DIR] [--contract PATH] [--output PATH]
Focused invocation only; exits 0 for every valid truth resolution.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import importlib.util
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

OWNER_DIR = Path(__file__).resolve().parent
BLOCK_MS = 900_000                 # non-overlapping UTC 15-minute blocks
ALIGN_TOL_MS = 5_000               # anchor age / final-latency / max-gap tolerance
MIN_IN_BLOCK_OBS = 720
MS_PER_YEAR = 31536000000.0        # 365 * 86400 * 1000 per final amendment
REQUIRED_INDEX_SERIES = ("btc_usd", "eth_usd")

DEFAULT_RAW = OWNER_DIR.parents[1] / "data/raw/focused/live/deribit"
DEFAULT_CONTRACT = OWNER_DIR.parent / "accelerated-intraday-contract.json"
DEFAULT_OUTPUT = OWNER_DIR / "fast_intraday_results.json"


def iso_ms(ts_ms: int) -> str:
    dt = datetime.fromtimestamp(ts_ms / 1000.0, tz=timezone.utc)
    return dt.strftime("%Y-%m-%dT%H:%M:%S.") + f"{int(ts_ms) % 1000:03d}Z"


def now_iso() -> str:
    return datetime.now(tz=timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def r12(x):
    if x is None:
        return None
    v = round(float(x), 12)
    return 0.0 if v == 0 else v


def load_run_module():
    """Reuse run.py's reconciliation-era BS gamma + parsing helpers verbatim."""
    spec = importlib.util.spec_from_file_location("focused_deribit_run", OWNER_DIR / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- snapshots

def snapshot_file(path: Path) -> dict:
    """Immutable-input-prefix rule: stat size BEFORE reading, then read exactly
    that many bytes; append-only file => prefix is stable. Hash parses match."""
    size = path.stat().st_size
    wall = datetime.now(tz=timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")
    with path.open("rb") as fh:
        prefix = fh.read(size)
    actual = len(prefix)
    return {
        "path": path,
        "prefix_bytes": size,
        "read_bytes": actual,
        "truncated_vs_stat": actual < size,
        "sha256": hashlib.sha256(prefix).hexdigest(),
        "snapshot_utc": wall,
        "prefix": prefix,
    }


def parse_prefix(prefix: bytes) -> tuple[list[dict], dict]:
    """Parse exactly the byte prefix; ignore one trailing partial line."""
    text = prefix.decode("utf-8", errors="replace")
    lines = text.split("\n")
    stats = {"complete_lines": 0, "blank_lines": 0, "partial_lines_ignored": 0,
             "malformed_json": 0, "records": 0}
    if lines and lines[-1].strip() != "":
        stats["partial_lines_ignored"] = 1
        lines = lines[:-1]
    recs: list[dict] = []
    for ln in lines:
        s = ln.strip()
        if not s:
            stats["blank_lines"] += 1
            continue
        stats["complete_lines"] += 1
        try:
            recs.append(json.loads(s))
        except json.JSONDecodeError:
            stats["malformed_json"] += 1
            continue
    stats["records"] = len(recs)
    return recs, stats


def exchange_ts_of(kind: str, rec: dict):
    """Per contract clock_semantics: the exchange-authored timestamp field."""
    if kind == "trades":
        return (rec.get("trade") or {}).get("timestamp")
    if kind == "index":
        return (rec.get("data") or {}).get("timestamp")
    if kind == "instruments":
        return (rec.get("instrument") or {}).get("creation_timestamp")
    return rec.get("received_ts_ms")


# ------------------------------------------------------------ replay + recon

def replay_order(recs: list[dict]) -> list[dict]:
    """Deterministic total order identical to run.py replay:
    (received_ts_ms, in-file seq, dedup_id)."""
    for r in recs:
        r["_seq"] = r.get("_seq", 0)
    recs.sort(key=lambda r: (r.get("received_ts_ms", 0), r["_seq"], r.get("dedup_id", "")))
    return recs


def dedupe(recs: list[dict]) -> tuple[list[dict], int]:
    seen: set[str] = set()
    out: list[dict] = []
    drops = 0
    for rec in recs:
        did = rec.get("dedup_id")
        if did is None or did in seen:
            drops += 1
            continue
        seen.add(did)
        out.append(rec)
    return out, drops


def build_instrument_map(recs: dict[str, list[dict]]) -> tuple[dict[str, dict], int]:
    """Latest-wins deterministic instrument metadata map (run.py reconcile rule)."""
    instruments: dict[str, dict] = {}
    n = 0
    for rec in recs["instruments"]:
        if rec.get("record_type") != "instrument":
            continue
        row = rec.get("instrument") or {}
        if row.get("instrument_name"):
            instruments[row["instrument_name"]] = row
            n += 1
    return instruments, n


# ----------------------------------------------------------------- GEX units

def compute_trade_gex(t: dict, meta: dict | None) -> dict:
    """final_semantic_amendment.deribit.gex_formula, exactly:
    S=trade.index_price, sigma=trade.iv/100, tau=max((expiry_ms-trade.timestamp)/
    31536000000,0), d1=(ln(S/K)+0.5*sigma^2*tau)/(sigma*sqrt(tau)),
    unit_gamma=pdf(d1)/(S*sigma*sqrt(tau)), coins=amount*contract_size,
    taker_gex=sign*coins*unit_gamma*S^2*0.01; dealer-mirror = negative sum.
    Finite/positivity requirements; missing or nonfinite inputs excluded+counted.
    Zero amount or gamma IS GEX-computable. Ticker gamma never used."""
    name = t.get("instrument_name") or ""
    u = name.split("-")[0] if name else ""
    out = {"underlying": u, "instrument": name, "ts": t.get("timestamp"),
           "ok": False, "reason": None, "taker_gex": None}
    ts = t.get("timestamp")
    if not isinstance(ts, (int, float)) or not math.isfinite(ts):
        out["reason"] = "missing_or_nonfinite_timestamp"
        return out
    if u not in ("BTC", "ETH"):
        out["reason"] = "unknown_underlying_prefix"
        return out
    iv = t.get("iv")
    if iv is None:
        out["reason"] = "missing_iv"
        return out
    amount = t.get("amount")
    if amount is None:
        out["reason"] = "missing_amount"
        return out
    s_usd = t.get("index_price")
    if s_usd is None:
        out["reason"] = "missing_index_price"
        return out
    direction = t.get("direction")
    if direction == "buy":
        sign = 1.0
    elif direction == "sell":
        sign = -1.0
    else:
        out["reason"] = "missing_or_unknown_direction"
        return out
    vals = {"iv": iv, "amount": amount, "index_price": s_usd}
    if any((not isinstance(v, (int, float))) or not math.isfinite(v) for v in vals.values()):
        out["reason"] = "nonfinite_trade_input"
        return out
    if amount < 0:
        out["reason"] = "negative_amount"
        return out
    if s_usd <= 0:
        out["reason"] = "nonpositive_index_price"
        return out
    if iv <= 0:
        out["reason"] = "nonpositive_iv"
        return out
    if meta is None:
        out["reason"] = "missing_instrument_metadata"
        return out
    strike = meta.get("strike")
    contract_size = meta.get("contract_size")
    expiry_ms = meta.get("expiration_timestamp")
    if strike is None or contract_size is None or expiry_ms is None:
        out["reason"] = "missing_metadata_strike_contract_or_expiry"
        return out
    if not isinstance(strike, (int, float)) or not math.isfinite(strike) or strike <= 0:
        out["reason"] = "invalid_strike"
        return out
    if (not isinstance(contract_size, (int, float)) or not math.isfinite(contract_size)
            or contract_size <= 0):
        out["reason"] = "invalid_contract_size"
        return out
    if not isinstance(expiry_ms, (int, float)) or not math.isfinite(expiry_ms):
        out["reason"] = "invalid_expiry"
        return out
    tau = max((expiry_ms - ts) / MS_PER_YEAR, 0.0)
    if tau <= 0:
        out["reason"] = "nonpositive_tau"
        return out
    gamma = RUN.bs_gamma(float(strike), float(s_usd), float(iv), tau)
    if gamma is None or not math.isfinite(gamma):
        out["reason"] = "invalid_gamma_inputs"
        return out
    coins = float(amount) * float(contract_size)
    gex = sign * coins * gamma * float(s_usd) * float(s_usd) * 0.01
    if not math.isfinite(gex):
        out["reason"] = "nonfinite_gex"
        return out
    out.update({"ok": True, "taker_gex": gex, "tau_years": tau, "gamma": gamma})
    return out


def build_index_series(recs: list[dict]) -> tuple[dict[str, tuple[list, list]], dict]:
    """Per series: dedupe equal exchange timestamps keeping the last in stable
    replay order, then sort by exchange ts. Nonpositive prices are retained so
    the affected block fails closed as unavailable."""
    counts = defaultdict(int)
    last = {}
    for rec in recs:
        if rec.get("record_type") != "index_price":
            continue
        d = rec.get("data") or {}
        name = d.get("index_name")
        ts = d.get("timestamp")
        counts["index_records"] += 1
        if name not in REQUIRED_INDEX_SERIES:
            counts["other_or_missing_series"] += 1
            continue
        if ts is None or not isinstance(ts, (int, float)) or not math.isfinite(ts):
            counts["missing_or_nonfinite_timestamp"] += 1
            continue
        counts[f"records_{name}"] += 1
        last[(name, ts)] = d.get("price")  # equal-ts keep-last in stable order
    series: dict[str, tuple[list, list]] = {}
    for name in REQUIRED_INDEX_SERIES:
        items = sorted(((ts, px) for (nm, ts), px in last.items() if nm == name),
                       key=lambda v: v[0])
        merged = counts[f"records_{name}"] - len(items)
        counts[f"duplicate_exchange_timestamps_merged_{name}"] = max(merged, 0)
        series[name] = ([ts for ts, _ in items], [px for _, px in items])
    return series, dict(counts)


# ---------------------------------------------------------------- RV blocks

def rv_for_block(series: dict[str, tuple[list, list]], name: str, k: int) -> dict:
    """rv_completeness (semantic_amendment.deribit): anchor <= block start within
    5000ms, final in-block obs >= end-5000 and < end, >=720 in-block obs, no
    consecutive gap > 5000ms; RV = sqrt(sum log-return^2) anchor->in-block."""
    ts_arr, px_arr = series[name]
    start = k * BLOCK_MS
    end = start + BLOCK_MS
    base = {"underlying_index": name, "block": k, "start_utc": iso_ms(start),
            "end_utc": iso_ms(end)}
    i = bisect.bisect_right(ts_arr, start) - 1
    if i < 0 or start - ts_arr[i] > ALIGN_TOL_MS:
        return {**base, "ok": False, "reason": "no_anchor_within_5000ms_before_block_start"}
    j0 = bisect.bisect_left(ts_arr, start)
    j1 = bisect.bisect_left(ts_arr, end)
    n_in = j1 - j0
    if n_in == 0:
        return {**base, "ok": False, "reason": "no_in_block_observations"}
    if n_in < MIN_IN_BLOCK_OBS:
        return {**base, "ok": False, "n_in_block": n_in,
                "reason": f"too_few_in_block_observations_lt_{MIN_IN_BLOCK_OBS}"}
    last_ts = ts_arr[j1 - 1]
    if last_ts < end - ALIGN_TOL_MS:
        return {**base, "ok": False, "n_in_block": n_in,
                "reason": "final_observation_more_than_5000ms_before_block_end"}
    seq = list(range(i, j1))
    max_gap = 0
    for a, b in zip(seq, seq[1:]):
        max_gap = max(max_gap, ts_arr[b] - ts_arr[a])
    if max_gap > ALIGN_TOL_MS:
        return {**base, "ok": False, "n_in_block": n_in,
                "reason": f"consecutive_gap_gt_5000ms_max_{max_gap}"}
    for j in seq:
        p = px_arr[j]
        if p is None or not isinstance(p, (int, float)) or not math.isfinite(p) or p <= 0:
            return {**base, "ok": False, "n_in_block": n_in,
                    "reason": "nonpositive_or_invalid_price"}
    logs = 0.0
    for a, b in zip(seq, seq[1:]):
        ratio = px_arr[b] / px_arr[a]
        logs += math.log(ratio) ** 2
    return {**base, "ok": True, "rv": math.sqrt(logs), "n_in_block": n_in,
            "anchor_ts": ts_arr[i], "anchor_utc": iso_ms(ts_arr[i]),
            "final_ts": last_ts, "final_utc": iso_ms(last_ts),
            "max_gap_ms": max_gap, "reason": None}


# ----------------------------------------------------------- block pipeline

def compute_blocks(gex_rows: list[dict], ts_max: int | None) -> dict:
    """Aggregate taker-signed GEX by (underlying, signal block); dealer = -sum."""
    agg: dict[tuple[str, int], dict] = {}
    excluded_block_counts = defaultdict(int)
    for row in gex_rows:
        if ts_max is not None and row["ts"] > ts_max:
            continue
        k = int(row["ts"]) // BLOCK_MS
        key = (row["underlying"], k)
        b = agg.setdefault(key, {"trades": 0, "computable": 0, "taker_sum": 0.0,
                                 "dealer_sum": 0.0,
                                 "excluded_reasons": defaultdict(int)})
        b["trades"] += 1
        if row["ok"]:
            b["computable"] += 1
            b["taker_sum"] += row["taker_gex"]
            b["dealer_sum"] = -b["taker_sum"]
        else:
            b["excluded_reasons"][row["reason"]] += 1
            excluded_block_counts[row["reason"]] += 1
    return agg


def build_transitions(agg: dict, series: dict) -> tuple[list[dict], dict]:
    """Eligible signal blocks (>=1 computable trade, nonzero dealer GEX) paired
    with the immediately following block's RV as the outcome. Zero-GEX blocks are
    reported and excluded; missing/unavailable RV means no transition (reported)."""
    eligible, zero_gex = [], []
    for (u, k), b in sorted(agg.items()):
        if b["computable"] >= 1 and b["dealer_sum"] != 0.0:
            eligible.append((u, k))
        elif b["computable"] >= 1:
            zero_gex.append({"underlying": u, "block": k, "block_start_utc": iso_ms(k * BLOCK_MS),
                             "dealer_gex": r12(b["dealer_sum"]), "computable_trades": b["computable"]})
    transitions, rv_journal = [], {}
    for (u, k) in eligible:
        idx_series = "btc_usd" if u == "BTC" else "eth_usd"
        r = rv_for_block(series, idx_series, k + 1)
        rv_journal[f"{u}:{k + 1}"] = r
        if r["ok"]:
            b = agg[(u, k)]
            transitions.append({
                "underlying": u, "signal_block": k,
                "signal_start_utc": iso_ms(k * BLOCK_MS),
                "outcome_block": k + 1,
                "outcome_start_utc": iso_ms((k + 1) * BLOCK_MS),
                "outcome_end_utc": iso_ms((k + 2) * BLOCK_MS),
                "outcome_end_ms": (k + 2) * BLOCK_MS,
                "computable_trades_in_signal_block": b["computable"],
                "trades_in_signal_block": b["trades"],
                "dealer_gex": b["dealer_sum"], "dealer_gex_rounded": r12(b["dealer_sum"]),
                "next_rv": r["rv"], "next_rv_rounded": r12(r["rv"]),
                "rv_n_in_block": r["n_in_block"], "rv_max_gap_ms": r["max_gap_ms"],
                "rv_final_utc": r["final_utc"],
            })
    return transitions, {"rv_journal": rv_journal, "zero_gex_blocks": zero_gex,
                         "eligible_signal_blocks": len(eligible)}


def conditions_at(transitions: list[dict], ss: dict) -> dict:
    n = len(transitions)
    per_asset = {u: sum(1 for t in transitions if t["underlying"] == u) for u in ("BTC", "ETH")}
    pos = sum(1 for t in transitions if t["dealer_gex"] > 0)
    neg = sum(1 for t in transitions if t["dealer_gex"] < 0)
    trades = sum(t["computable_trades_in_signal_block"] for t in transitions)
    return {
        "gex_computable_trades": trades, "total_transitions": n,
        "transitions_btc": per_asset["BTC"], "transitions_eth": per_asset["ETH"],
        "positive_dealer_gex_blocks": pos, "negative_dealer_gex_blocks": neg,
        "met": (trades >= ss["minimum_gex_computable_trades"]
                and n >= ss["minimum_complete_asset_block_transitions"]
                and per_asset["BTC"] >= ss["minimum_transitions_per_asset"]
                and per_asset["ETH"] >= ss["minimum_transitions_per_asset"]
                and pos >= ss["minimum_positive_dealer_gex_blocks"]
                and neg >= ss["minimum_negative_dealer_gex_blocks"]),
        "conditions_detail": {
            "minimum_gex_computable_trades": {"observed": trades, "minimum": ss["minimum_gex_computable_trades"], "pass": trades >= ss["minimum_gex_computable_trades"]},
            "minimum_complete_asset_block_transitions": {"observed": n, "minimum": ss["minimum_complete_asset_block_transitions"], "pass": n >= ss["minimum_complete_asset_block_transitions"]},
            "minimum_transitions_per_asset": {"observed": per_asset, "minimum": ss["minimum_transitions_per_asset"], "pass": per_asset["BTC"] >= ss["minimum_transitions_per_asset"] and per_asset["ETH"] >= ss["minimum_transitions_per_asset"]},
            "minimum_positive_dealer_gex_blocks": {"observed": pos, "minimum": ss["minimum_positive_dealer_gex_blocks"], "pass": pos >= ss["minimum_positive_dealer_gex_blocks"]},
            "minimum_negative_dealer_gex_blocks": {"observed": neg, "minimum": ss["minimum_negative_dealer_gex_blocks"], "pass": neg >= ss["minimum_negative_dealer_gex_blocks"]},
        },
    }


def first_failure(c: dict, ss: dict) -> str | None:
    d = c["conditions_detail"]
    if d["minimum_gex_computable_trades"]["pass"] is False:
        return "minimum_gex_computable_trades"
    if d["minimum_complete_asset_block_transitions"]["pass"] is False:
        return "minimum_complete_asset_block_transitions"
    if d["minimum_transitions_per_asset"]["pass"] is False:
        return "minimum_transitions_per_asset"
    if d["minimum_positive_dealer_gex_blocks"]["pass"] is False:
        return "minimum_positive_dealer_gex_blocks"
    if d["minimum_negative_dealer_gex_blocks"]["pass"] is False:
        return "minimum_negative_dealer_gex_blocks"
    return None


def find_earliest_stop(transitions: list[dict], ss: dict, deadline_ms: int) -> tuple[int | None, list]:
    """earliest_stop: order by outcome-block end; earliest end meeting every
    registered minimum, searched only at or before the frozen deadline."""
    ends = sorted({t["outcome_end_ms"] for t in transitions if t["outcome_end_ms"] <= deadline_ms})
    trail = []
    stop = None
    for e in ends:
        subset = [t for t in transitions if t["outcome_end_ms"] <= e]
        c = conditions_at(subset, ss)
        fails = first_failure(c, ss)
        trail.append({"outcome_end_ms": e, "outcome_end_utc": iso_ms(e),
                      "cumulative": {k: v for k, v in c.items() if k not in ("met", "conditions_detail")},
                      "all_conditions_met": c["met"], "first_unmet": fails})
        if c["met"]:
            stop = e
            break
    return stop, trail


# ----------------------------------------------------------- rank statistics

def avg_ranks(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: (vals[i], i))
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for kk in range(i, j + 1):
            ranks[order[kk]] = avg
        i = j + 1
    return ranks


def pearson(x: list[float], y: list[float]):
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0.0 or syy == 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def spearman(x: list[float], y: list[float]):
    """Average ranks for ties + Pearson of ranks; UNKNOWN (<2 obs or constant
    rank vector) returns None per semantic_amendment.deribit.spearman."""
    if len(x) < 2 or len(x) != len(y):
        return None
    return pearson(avg_ranks(x), avg_ranks(y))


def evaluate_gates(transitions: list[dict], gate_cfg: dict) -> tuple[dict, dict]:
    xs = [t["dealer_gex"] for t in transitions]
    ys = [t["next_rv"] for t in transitions]
    n = len(transitions)

    rho = spearman(xs, ys)
    g1 = {
        "gate": "spearman_dealer_gex_vs_next_rv_lte",
        "threshold": gate_cfg["spearman_dealer_gex_vs_next_rv_lte"],
        "n": n,
        "value": r12(rho), "extended_value": None,
        "pass": (None if rho is None else rho <= gate_cfg["spearman_dealer_gex_vs_next_rv_lte"]),
    }

    neg = [t["next_rv"] for t in transitions if t["dealer_gex"] < 0]
    pos = [t["next_rv"] for t in transitions if t["dealer_gex"] > 0]
    thr2 = gate_cfg["negative_vs_positive_dealer_gex_next_rv_ratio_gte"]
    if not neg or not pos:
        g2 = {"gate": "negative_vs_positive_dealer_gex_next_rv_ratio_gte", "threshold": thr2,
              "n_negative_cohort": len(neg), "n_positive_cohort": len(pos),
              "value": None, "extended_value": None, "pass": None}
    else:
        med_neg = statistics.median(neg)
        med_pos = statistics.median(pos)
        if med_pos == 0.0:
            if med_neg > 0.0:
                # json_infinity: null numeric, exact string "+infinity", comparison true
                g2 = {"gate": "negative_vs_positive_dealer_gex_next_rv_ratio_gte", "threshold": thr2,
                      "n_negative_cohort": len(neg), "n_positive_cohort": len(pos),
                      "median_negative_rv": r12(med_neg), "median_positive_rv": 0.0,
                      "value": None, "extended_value": "+infinity", "pass": True}
            else:
                g2 = {"gate": "negative_vs_positive_dealer_gex_next_rv_ratio_gte", "threshold": thr2,
                      "n_negative_cohort": len(neg), "n_positive_cohort": len(pos),
                      "median_negative_rv": 0.0, "median_positive_rv": 0.0,
                      "value": None, "extended_value": None, "pass": None}
        else:
            ratio = med_neg / med_pos
            g2 = {"gate": "negative_vs_positive_dealer_gex_next_rv_ratio_gte", "threshold": thr2,
                  "n_negative_cohort": len(neg), "n_positive_cohort": len(pos),
                  "median_negative_rv": r12(med_neg), "median_positive_rv": r12(med_pos),
                  "value": r12(ratio), "extended_value": None, "pass": ratio >= thr2}

    per_asset = {}
    for u, min_t in (("BTC", None), ("ETH", None)):
        sub = [(t["dealer_gex"], t["next_rv"]) for t in transitions if t["underlying"] == u]
        rho_u = spearman([a for a, _ in sub], [b for _, b in sub]) if len(sub) >= 2 else None
        per_asset[u] = {"n": len(sub), "spearman": r12(rho_u), "defined": rho_u is not None}
    all_def = per_asset["BTC"]["defined"] and per_asset["ETH"]["defined"]
    g3 = {"gate": "both_asset_spearman_lt_zero",
          "btc": {"n": per_asset["BTC"]["n"], "spearman": per_asset["BTC"]["spearman"]},
          "eth": {"n": per_asset["ETH"]["n"], "spearman": per_asset["ETH"]["spearman"]},
          "value": None if not all_def else {"btc": per_asset["BTC"]["spearman"],
                                             "eth": per_asset["ETH"]["spearman"]},
          "extended_value": None,
          "pass": (None if not all_def else
                   (per_asset["BTC"]["spearman"] < 0 and per_asset["ETH"]["spearman"] < 0))}

    defined = 0
    nonpos = 0
    any_undefined = False
    for i in range(n):
        xs2 = xs[:i] + xs[i + 1:]
        ys2 = ys[:i] + ys[i + 1:]
        rho_i = spearman(xs2, ys2)
        if rho_i is None:
            any_undefined = True
        else:
            defined += 1
            if rho_i <= 0:
                nonpos += 1
    thr4 = gate_cfg["leave_one_block_out_nonpositive_fraction_gte"]
    if any_undefined or defined == 0:
        g4 = {"gate": "leave_one_block_out_nonpositive_fraction_gte", "threshold": thr4,
              "n": n, "defined_folds": defined, "undefined_folds": n - defined,
              "nonpositive_folds": nonpos, "value": None, "extended_value": None, "pass": None}
    else:
        frac = nonpos / defined
        g4 = {"gate": "leave_one_block_out_nonpositive_fraction_gte", "threshold": thr4,
              "n": n, "defined_folds": defined, "undefined_folds": 0,
              "nonpositive_folds": nonpos, "value": r12(frac),
              "extended_value": None, "pass": frac >= thr4}

    gates = {"spearman_dealer_gex_vs_next_rv_lte": g1,
             "negative_vs_positive_dealer_gex_next_rv_ratio_gte": g2,
             "both_asset_spearman_lt_zero": g3,
             "leave_one_block_out_nonpositive_fraction_gte": g4}
    extras = {"per_asset_spearman_detail": per_asset, "pooled_n": n}
    return gates, extras


def resolve_verdict(conditions_detail: dict, gates: dict) -> tuple[str, list[str]]:
    reasons: list[str] = []
    for name, c in conditions_detail.items():
        if c["pass"] is not True:
            reasons.append(f"sample condition not met: {name} "
                           f"(observed {json.dumps(c['observed'])} vs {json.dumps(c['minimum'])})")
    if reasons:
        return "inconclusive", reasons
    for name, g in gates.items():
        if g["pass"] is False:
            reasons.append(f"gate observed false: {name}")
    if reasons:
        return "falsified", reasons
    for name, g in gates.items():
        if g["pass"] is None:
            reasons.append(f"gate UNKNOWN (unavailable value): {name}")
    if reasons:
        return "inconclusive", reasons
    return "fast-screen-pass", ["every registered sample condition and gate observed true"]


# --------------------------------------------------------------------- main

RUN = None  # set in main(); holds run.py module for bs_gamma reuse


def main(argv=None) -> int:
    global RUN
    RUN = load_run_module()
    ap = argparse.ArgumentParser(description="Frozen fast screen: deribit_intraday_taker_gex")
    ap.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    ap.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args(argv)

    invocation_utc = now_iso()
    contract_bytes = args.contract.read_bytes()
    contract = json.loads(contract_bytes)
    screen = contract["screens"]["deribit_intraday_taker_gex"]
    amend = contract["semantic_amendment"]
    final = contract["final_semantic_amendment"]
    ss = screen["sample_stop"]
    gate_cfg = screen["gate"]
    deadline_ms = RUN.parse_utc_ms(final["deadline_utc"])

    # --- 1. Snapshot every live JSONL input (day dirs, required record files) ----
    day_dirs = sorted(p for p in args.raw_dir.iterdir() if p.is_dir() and p.name.startswith("day"))
    snap_infos = []
    parsed: dict[str, list[dict]] = {"trades": [], "instruments": [], "index": []}
    parse_stats = {}
    exchanges = {}
    for d in day_dirs:
        for fname, kind in (("trades.jsonl", "trades"), ("instruments.jsonl", "instruments"),
                            ("index.jsonl", "index")):
            fp = d / fname
            if not fp.is_file():
                continue
            snap = snapshot_file(fp)
            recs, stats = parse_prefix(snap.pop("prefix"))
            # replay-order marker: file-line order within same received timestamp
            for seq, r in enumerate(recs):
                r["_seq"] = seq
            parsed[kind].extend(recs)
            first_ts = None
            last_ts = None
            per_series_latest: dict[str, int] = {}
            for r in recs:
                if kind == "index":
                    d = r.get("data") or {}
                    nm = d.get("index_name")
                    ets = d.get("timestamp")
                    ts = ets
                else:
                    ts = exchange_ts_of(kind, r)
                if isinstance(ts, (int, float)) and math.isfinite(ts):
                    first_ts = ts if first_ts is None else min(first_ts, ts)
                    last_ts = ts if last_ts is None else max(last_ts, ts)
                    if kind == "index" and nm in ("btc_usd", "eth_usd"):
                        prev = per_series_latest.get(nm)
                        per_series_latest[nm] = ts if prev is None else max(prev, ts)
            exchanges[kind] = {"file": str(fp), "first_exchange_ts_ms": first_ts,
                               "latest_exchange_ts_ms": last_ts,
                               "per_series_latest_exchange_ts_ms": per_series_latest}
            parse_stats[str(fp)] = stats
            snap_infos.append({
                "path": str(fp), "prefix_bytes": snap["prefix_bytes"],
                "read_bytes": snap["read_bytes"],
                "truncated_vs_stat": bool(snap["truncated_vs_stat"]),
                "sha256": snap["sha256"], "snapshot_utc": snap["snapshot_utc"],
                "exchange_clock_field": {"trades": "trade.timestamp",
                                         "instruments": "instrument.creation_timestamp",
                                         "index": "data.timestamp"}[kind],
                "first_parsed_exchange_ts_ms": first_ts,
                "latest_parsed_exchange_ts_ms": last_ts,
                "first_parsed_exchange_ts_utc": None if first_ts is None else iso_ms(first_ts),
                "latest_parsed_exchange_ts_utc": None if last_ts is None else iso_ms(last_ts),
                "parse": stats,
            })
    if not snap_infos:
        raise SystemExit(f"no capture files found under {args.raw_dir}")

    # --- 2. Common cutoff: min of (trade latest, BTC latest, ETH latest), fail closed
    trade_latest = exchanges.get("trades", {}).get("latest_exchange_ts_ms")
    idx = exchanges.get("index", {}).get("per_series_latest_exchange_ts_ms") or {}
    btc_latest = idx.get("btc_usd")
    eth_latest = idx.get("eth_usd")
    if not trade_latest or not btc_latest or not eth_latest:
        raise SystemExit("common cutoff unavailable: missing trade/BTC/ETH latest exchange "
                         f"timestamps (trades={trade_latest}, btc_usd={btc_latest}, "
                         f"eth_usd={eth_latest}) — fail closed per contract")
    cutoff = min(trade_latest, btc_latest, eth_latest)
    per_series_cutoffs = {"trades": trade_latest, "btc_usd": btc_latest, "eth_usd": eth_latest}

    # --- 3. Deterministic replay order + dedup by dedup_id (run.py semantics)
    all_recs = parsed["trades"] + parsed["instruments"] + parsed["index"]
    replay_order(all_recs)
    deduped, dup_drops = dedupe(all_recs)
    recs_by_kind: dict[str, list[dict]] = {"trades": [], "instruments": [], "index": []}
    for r in deduped:
        rt = r.get("record_type")
        if rt == "option_trade":
            recs_by_kind["trades"].append(r)
        elif rt == "instrument":
            recs_by_kind["instruments"].append(r)
        elif rt == "index_price":
            recs_by_kind["index"].append(r)

    # deadline_rule: earliest-stop search may only use observations available at or
    # before the frozen deadline; never extend sampling past it.
    analysis_cutoff = min(cutoff, deadline_ms)
    pre_cutoff = {k: len(v) for k, v in recs_by_kind.items()}
    recs_by_kind["trades"] = [r for r in recs_by_kind["trades"]
                              if (r.get("trade") or {}).get("timestamp", 0) <= analysis_cutoff]
    recs_by_kind["index"] = [r for r in recs_by_kind["index"]
                             if ((r.get("data") or {}).get("timestamp") or 0) <= analysis_cutoff]
    cutoff_dropped = {k: pre_cutoff[k] - len(recs_by_kind[k]) for k in pre_cutoff}

    instruments, n_meta = build_instrument_map(recs_by_kind)

    # --- 4. Per-trade GEX exactly per frozen formula
    reason_counts = defaultdict(int)
    gex_rows = []
    for r in recs_by_kind["trades"]:
        t = r.get("trade") or {}
        row = compute_trade_gex(t, instruments.get(t.get("instrument_name")))
        row_dropped_after_cutoff = row["ts"] is not None and row["ts"] > analysis_cutoff
        gex_rows.append(row)
        if not row["ok"] and not row_dropped_after_cutoff:
            reason_counts[row["reason"]] += 1
    n_trades_total = len(gex_rows)
    n_trades_computable = sum(1 for row in gex_rows if row["ok"])

    series, index_counts = build_index_series(recs_by_kind["index"])

    # --- 5. Blocks + transitions; earliest stop scan bounded by frozen deadline
    agg_full = compute_blocks(gex_rows, ts_max=None)
    series = truncate_series(series, analysis_cutoff)
    transitions_full, aux_full = build_transitions(agg_full, series)
    stop_ms, trail = find_earliest_stop(transitions_full, ss, deadline_ms)

    ts_max = stop_ms if stop_ms is not None else None
    agg = compute_blocks(gex_rows, ts_max=ts_max)
    series_stop = truncate_series(series, stop_ms) if stop_ms is not None else series
    transitions, aux = build_transitions(agg, series_stop)
    if stop_ms is not None:
        assert [(t["underlying"], t["signal_block"], t["outcome_block"]) for t in transitions] == \
               [(t["underlying"], t["signal_block"], t["outcome_block"]) for t in transitions_full
                if t["outcome_end_ms"] <= stop_ms], "truncation rerun mismatch"

    cond = conditions_at(transitions, ss)
    gates, gate_extras = evaluate_gates(transitions, gate_cfg)
    verdict, verdict_reasons = resolve_verdict(cond["conditions_detail"], gates)

    # --- 6. Frozen report -----------------------------------------------------
    fm = final["deribit"]
    results = {
        "schema_version": 1,
        "screen_id": "deribit_intraday_taker_gex",
        "direction_id": "focused-deribit-gex",
        "owner_dir": screen["owner_dir"],
        "script": "fast_intraday.py",
        "output": str(args.output.name),
        "fast_screen_verdict": verdict,
        "verdict_reasons": verdict_reasons,
        "created_utc": invocation_utc,
        "claim_boundary": contract["claim_boundary"],
        "screen_question": screen["question"],
        "frozen_semantics_applied": {
            "contract": {
                "path": str(args.contract),
                "sha256": hashlib.sha256(contract_bytes).hexdigest(),
                "created_at_utc": contract["created_at_utc"],
                "status": contract["status"],
                "semantic_amendment_created_at_utc": amend["created_at_utc"],
                "final_semantic_amendment_created_at_utc": final["created_at_utc"],
                "final_schema_clarified_at_utc": final["schema_clarified_at_utc"],
            },
            "deribit_amendment": amend["deribit"],
            "deribit_final_amendment": fm,
            "deadline_utc": final["deadline_utc"],
            "deadline_rule": final["deadline_rule"],
            "deadline_respected": (stop_ms is None) or (stop_ms <= deadline_ms),
            "json_infinity_rule": final["deribit"]["json_infinity"],
            "gex_units_note": "reuses run.py units (coin amounts, index_price USD, USD-per-1%-GEX) with the amendment-mandated BS-from-trade-IV gamma source",
        },
        "input_prefixes": {
            "rule": amend["immutable_input_prefix"]["rule"],
            "common_cutoff_rule": amend["immutable_input_prefix"]["common_cutoff"],
            "earliest_stop_rule": amend["immutable_input_prefix"]["earliest_stop"],
            "snapshot_wall_clock_utc": min(s["snapshot_utc"] for s in snap_infos),
            "files": snap_infos,
            "common_cutoff_ms": cutoff,
            "common_cutoff_utc": iso_ms(cutoff),
            "common_cutoff_rule": "min(trades latest, BTC index latest, ETH index latest); fails closed if any is missing",
            "per_series_latest_exchange_ts_ms": per_series_cutoffs,
            "records_after_cutoff_drop": cutoff_dropped,
        },
        "reconciliation": {
            "replay_order": "(received_ts_ms, in-file seq, dedup_id) — identical to run.py",
            "dedup_rule": "existing deterministic dedup_id; duplicates and missing-id records dropped",
            "lines_read_complete": {k: v for k, v in parse_stats.items()},
            "records_by_kind_after_dedup_and_cutoff": {k: len(v) for k, v in recs_by_kind.items()},
            "duplicate_drops": dup_drops,
            "instruments_metadata_map": n_meta,
            "units": {
                "amount": "base coin: option trade 'amount' x contract_size coins (run.py units; NOT USD)",
                "index_price": "USD per coin, embedded per trade by Deribit",
                "gex_formula": "taker_gex = sign(direction) x coins x BS unit gamma(trade IV, metadata strike/expiry) x index_price^2 x 0.01 [USD per 1% move]; dealer-mirror block GEX = -sum(taker_gex)",
                "gamma_policy": "final_semantic_amendment.deribit.gex_formula only — Black-Scholes from the trade's own IV; ticker/venue gamma NEVER used",
                "bs_implementation": "reused verbatim from run.py:bs_gamma (r=0, q=0, ACT/365, normal pdf)",
            },
        },
        "gex_trade_coverage": {
            "option_trades_parsed": n_trades_total,
            "gex_computable_total": n_trades_computable,
            "excluded_by_reason": {k: v for k, v in sorted(reason_counts.items())},
            "fail_closed_policy": "any trade missing or carrying nonfinite/invalid inputs is excluded from GEX sums and counted above; zero amount or zero gamma is GEX-computable per final amendment; nothing imputed",
        },
        "index_series": {
            "counts": index_counts,
            "rv_completeness_rule": amend["deribit"]["rv_completeness"],
            "required_in_block_observations": MIN_IN_BLOCK_OBS,
        },
        "sample_stop": {
            "registered": ss,
            "earliest_stop_rule": amend["deribit"]["earliest_stop"],
            "stop_outcome_end_ms": stop_ms,
            "stop_outcome_end_utc": None if stop_ms is None else iso_ms(stop_ms),
            "stop_basis": ("earliest outcome-block end at which all registered sample conditions "
                           "first became true, searched only at or before the frozen deadline"
                           if stop_ms is not None else
                           "no registered earliest stop found: no outcome-block end satisfied "
                           "all registered sample conditions at or before the frozen deadline; "
                           "analysis uses the full frozen prefix through "
                           "min(common cutoff, deadline)"),
            "candidate_trail": trail,
        },
        "block_counts": build_block_report(agg, aux, transitions, bool(stop_ms is None)),
        "transitions": {
            "unit_of_observation": amend["deribit"]["unit_of_observation"],
            "signal_eligibility": amend["deribit"]["signal_eligibility"],
            "zero_gex_blocks_reported": aux["zero_gex_blocks"],
            "rv_journal": {k: {kk: vv for kk, vv in v.items()} for k, v in aux["rv_journal"].items()},
            "included": [{k: (r12(v) if k in ("dealer_gex", "next_rv") else v)
                          for k, v in t.items()} for t in transitions],
            "included_count": len(transitions),
        },
        "sample_conditions": {
            "note": ("evaluated on the included cohort at the frozen earliest stop"
                     if stop_ms is not None else
                     "evaluated on the included cohort over the full frozen prefix through "
                     "min(common cutoff, deadline); no registered earliest stop exists"),
            "conditions": cond["conditions_detail"],
            "all_met": cond["met"],
        },
        "gates": {
            "registered": gate_cfg,
            "evaluated": gates,
            "pooled_observations": gate_extras["pooled_n"],
            "per_asset_detail": gate_extras["per_asset_spearman_detail"],
        },
        "decision": {
            "deadline_ms": deadline_ms,
            "deadline_utc": final["deadline_utc"],
            "common_cutoff_utc": iso_ms(cutoff),
            "analysis_cutoff_ms": analysis_cutoff,
            "analysis_cutoff_utc": iso_ms(analysis_cutoff),
            "analysis_cutoff_rule": "min(common_cutoff, frozen deadline_utc) — earliest-stop search may only use observations available at or before the deadline",
            "records_after_cutoff_drop": cutoff_dropped,
            "truth_table": screen["truth_table"],
            "unavailable_rule": amend["deribit"]["unavailable"],
            "verdict": verdict,
            "verdict_reasons": verdict_reasons,
            "no_returns": "no returns.csv and no return stream is produced by this screen; cumulative 15-minute realized volatility of a public index is a risk descriptive, not a tradable strategy return",
            "no_deployment_claim": True,
            "no_alpha_claim": "a fast screen may falsify or prioritize a mechanism; it cannot establish durable alpha or replace prospective confirmation",
            "no_imputation": "missing-block/no-anchor/incomplete blocks are reported with reasons and never imputed",
            "existing_week_gate_untouched": "run.py four-complete-registered-weeks gate and registry remain the durable-alpha path; this fast screen does not modify them",
        },
    }
    args.output.write_text(json.dumps(results, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "screen": "deribit_intraday_taker_gex",
        "verdict": verdict,
        "trades_parsed": n_trades_total,
        "gex_computable_total": n_trades_computable,
        "transitions": {"total": len(transitions),
                        "btc": sum(1 for t in transitions if t["underlying"] == "BTC"),
                        "eth": sum(1 for t in transitions if t["underlying"] == "ETH")},
        "stop_utc": None if stop_ms is None else iso_ms(stop_ms),
        "deadline_utc": final["deadline_utc"],
        "output": str(args.output),
    }, indent=1))
    return 0


def truncate_series(series: dict, ts_max: int | None) -> dict:
    if ts_max is None:
        return series
    out = {}
    for name, (ts_arr, px_arr) in series.items():
        keep = [i for i, t in enumerate(ts_arr) if t <= ts_max]
        out[name] = ([ts_arr[i] for i in keep], [px_arr[i] for i in keep])
    return out


def build_block_report(agg: dict, aux: dict, transitions: list[dict], full_prefix: bool) -> dict:
    per_u: dict[str, dict] = {}
    for (u, k), b in sorted(agg.items()):
        entry = per_u.setdefault(u, {"blocks_with_trades": 0, "blocks_eligible": 0,
                                     "trades_total": 0, "computable_total": 0,
                                     "zero_gex_blocks": 0, "block_detail": {}})
        entry["blocks_with_trades"] += 1
        entry["trades_total"] += b["trades"]
        entry["computable_total"] += b["computable"]
        eligible = b["computable"] >= 1 and b["dealer_sum"] != 0.0
        zero = b["computable"] >= 1 and b["dealer_sum"] == 0.0
        if eligible:
            entry["blocks_eligible"] += 1
        if zero:
            entry["zero_gex_blocks"] += 1
        entry["block_detail"][str(k)] = {
            "block_start_utc": iso_ms(k * BLOCK_MS),
            "trades": b["trades"], "computable": b["computable"],
            "dealer_gex": r12(b["dealer_sum"]),
            "eligible_signal_block": eligible, "zero_gex_block": zero,
            "excluded_reasons": dict(b["excluded_reasons"]),
        }
    for u in per_u:
        per_u[u]["transitions_included"] = sum(1 for t in transitions if t["underlying"] == u)
        per_u[u]["positive_dealer_gex_transitions"] = sum(
            1 for t in transitions if t["underlying"] == u and t["dealer_gex"] > 0)
        per_u[u]["negative_dealer_gex_transitions"] = sum(
            1 for t in transitions if t["underlying"] == u and t["dealer_gex"] < 0)
        detail = per_u[u].pop("block_detail")
        per_u[u]["block_detail"] = {k: detail[k] for k in sorted(detail, key=int)}
    return {"per_underlying": per_u,
            "note": "block_detail keyed by absolute UTC 15-minute block index floor(exchange_ms/900000); signal blocks without a completed passing outcome block have no transition (reported, not imputed)"}


if __name__ == "__main__":
    raise SystemExit(main())
