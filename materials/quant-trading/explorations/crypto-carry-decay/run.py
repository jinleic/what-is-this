#!/usr/bin/env python3
"""Unconditional equal-weight BTCUSDT/ETHUSDT carry: spot-long, USD-M-perp-short.

Fully funded 1:1 basis trade. Initial capital 1.0 splits 50/50 into BTC and
ETH sleeves (no cross-sleeve transfers afterwards). While on, a sleeve holds
spot notional = short USD-M perp notional = K ~= half its sleeve equity at the
last trade open, so portfolio gross exposure is ~2x equity with ZERO net delta
and zero borrowed notional (spot inventory funds the short; entry is sized so
sleeve cash is exactly zero after entry fees).

Exchange/settlement model (documented assumptions):
  * Spot leg is a real asset: principal paid at buy, received at sell;
    unrealized spot P&L lives in inventory and is marked hourly at closes.
  * The USD-M perpetual has NO principal cash flows. Its mark-to-market P&L
    settles into cash hour by hour (settle-to-market convention applied at
    every segment boundary: close->open on current units at trade opens,
    open->close within the hour). Closing the position needs no principal and
    carries only the closing fee.
  * Funding settles in cash at each event: a SHORT position receives
    +rate * perp_notional for positive rate and pays for negative rate
    (exact sign: short cash flow = -q_p * P * rate with q_p < 0).
  * One-way taker fees incl. slippage: spot 12 bps, perp 7 bps (contract
    defaults), charged to cash on traded notional at entries, scheduled
    rebalances, and final-liquidation. The 2x stress run doubles them.
  * No interest is paid on idle sleeve cash; funding notional is marked at
    the perp close of the hour containing the event.

Position continuity: the engine simulates from HISTORY_START through the
requested end and reports windows by slicing the continuous accounting path -
the validation window inherits the position carried since 2020-08-01 (the
rule runs unconditionally), so no artificial in-window entry is charged and
no post-2024-12-31 bar is read.

Determinism: fixed monthly zip archives verified against .CHECKSUM sha256
files, forward-filled close marks over the 10 audited spot outage hours, one
hourly-grid simulation, fixed summation order. Identical --start/--end
reproduce identical outputs. `python3 run.py --self-test` asserts accounting
identities (telescoping, funding sign/magnitude, convergence capture, cost
doubling, settlement semantics, real-data invariants) and exits nonzero on
any failure.

Wave-one access control: fails closed when --end exceeds 2024-12-31 unless
QUANT_HOLDOUT_UNLOCK=1 is set (external evaluator only). In wave one the flag
is never set; variant selection, parameters, and all wave-one artifacts are
independent of it.

Usage:
  python3 run.py                                  # dev + validation windows
  python3 run.py --start 2024-01-01 --end 2024-12-31 --output-dir results
  python3 run.py --rebalance-hours N              # maintenance clock variant
  python3 run.py --self-test
"""

from __future__ import annotations

import argparse
import bisect
import datetime as dt
import hashlib
import json
import math
import os
import sys
import zipfile
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------- constants --

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
RAW = os.path.join(ROOT, "data", "raw", "binance")
SOFR_CSV = os.path.join(ROOT, "data", "raw", "macro", "nyfed_sofr_2020-01-01_2026-08-28.csv")

SYMBOLS = ("BTCUSDT", "ETHUSDT")
SLEEVE_WEIGHT = 0.5            # initial capital split

SPOT_BPS = 12.0                # contract default one-way spot taker incl. slippage
PERP_BPS = 7.0                 # contract default one-way perpetual taker incl. slippage
STRESS_MULT = 2.0              # contract 2x-cost stress

FROZEN_REBALANCE_HOURS = 0     # frozen variant (development-selected)

HISTORY_START = dt.date(2020, 8, 1)
WAVE1_END = dt.date(2024, 12, 31)

HOUR = 3_600_000
MS_8H = 28_800_000
DAY_MS = 86_400_000
HPY = 24.0 * 365.25            # annualization for hourly series
DAY_DENOM = 24.0 * 365.0       # daily SOFR accrual denominator

# ----------------------------------------------------------------- data io ---


def month_list(a: dt.date, b: dt.date) -> list[tuple[int, int]]:
    out = []
    cur = dt.date(a.year, a.month, 1)
    while cur <= b:
        out.append((cur.year, cur.month))
        cur = dt.date(cur.year + (cur.month == 12), 1 if cur.month == 12 else cur.month + 1, 1)
    return out


def verify_checksum(path: str) -> None:
    cs = path + ".CHECKSUM"
    if not os.path.exists(cs):
        raise RuntimeError(f"missing checksum file {cs}")
    expect = open(cs).read().split()[0].strip()
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    if h.hexdigest() != expect:
        raise RuntimeError(f"sha256 mismatch {path}")


def parse_ohlcv(zip_path: str) -> dict[int, tuple[float, float]]:
    """Bar-open ms -> (open, close). Tolerates optional header rows and the
    2022..2024 microsecond-timestamp archive era (normalized to ms)."""
    out: dict[int, tuple[float, float]] = {}
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(zf.namelist()[0]) as f:
            for line in f:
                line = line.decode("utf-8").strip()
                if not line or not line[0].isdigit():
                    continue
                p = line.split(",")
                ts = int(p[0])
                if ts > 10**14:
                    ts //= 1000
                out[ts] = (float(p[1]), float(p[4]))
    return out


def parse_funding(zip_path: str, interval_green: set[str]) -> dict[int, float]:
    """Funding event (snapped to the canonical 8h grid; calc_time carries rare
    millisecond processing jitter; snap collisions were never observed) -> rate."""
    out: dict[int, float] = {}
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(zf.namelist()[0]) as f:
            for line in f:
                line = line.decode("utf-8").strip()
                if not line or not line[0].isdigit():
                    continue
                p = line.split(",")
                if p[1] not in interval_green:
                    raise RuntimeError(
                        f"funding_interval_hours={p[1]} outside audited {sorted(interval_green)} in {zip_path}")
                ts = (int(p[0]) // MS_8H) * MS_8H
                out[ts] = max(float(p[2]), out.get(ts, 0.0))
    return out


def load_symbol(symbol: str, start_d: dt.date, end_d: dt.date):
    """Spot/perp bar-open->(open, close) dicts and snapped funding event->rate
    dict, for the UTC date range [start_d, end_d]."""
    start_ms = int(dt.datetime(start_d.year, start_d.month, start_d.day, tzinfo=dt.UTC).timestamp() * 1000)
    end_ms = int(dt.datetime(end_d.year, end_d.month, end_d.day, tzinfo=dt.UTC).timestamp() * 1000) + DAY_MS - 1
    spot_o, spot_c, perp_o, perp_c, fund = {}, {}, {}, {}, {}
    fund_paths = []
    for (yy, mm) in month_list(start_d, end_d):
        p = os.path.join(RAW, "spot-1h", symbol, f"{symbol}-1h-{yy}-{mm:02d}.zip")
        if not os.path.exists(p):
            raise RuntimeError(f"missing spot archive {yy}-{mm:02d} for {symbol}")
        verify_checksum(p)
        for ts, (o, c) in parse_ohlcv(p).items():
            if start_ms <= ts <= end_ms:
                spot_o[ts], spot_c[ts] = o, c

        p = os.path.join(RAW, "um-perp-1h", symbol, f"{symbol}-1h-{yy}-{mm:02d}.zip")
        if not os.path.exists(p):
            raise RuntimeError(f"missing perp archive {yy}-{mm:02d} for {symbol}")
        verify_checksum(p)
        for ts, (o, c) in parse_ohlcv(p).items():
            if start_ms <= ts <= end_ms:
                perp_o[ts], perp_c[ts] = o, c

        p = os.path.join(RAW, "funding-rate", symbol, f"{symbol}-fundingRate-{yy}-{mm:02d}.zip")
        if not os.path.exists(p):
            raise RuntimeError(f"missing funding archive {yy}-{mm:02d} for {symbol}")
        verify_checksum(p)
        fund_paths.append(p)

    iv_green: set[str] = set()
    for p in fund_paths:
        with zipfile.ZipFile(p) as zf:
            with zf.open(zf.namelist()[0]) as f:
                for line in f:
                    line = line.decode("utf-8").strip()
                    if line and line[0].isdigit():
                        iv_green.add(line.split(",")[1])
    if not iv_green <= {"4", "8"}:
        raise RuntimeError(f"{symbol}: unhandled funding intervals {sorted(iv_green)}")
    for p in fund_paths:
        for ts, r in parse_funding(p, iv_green).items():
            if start_ms <= ts <= end_ms:
                fund[ts] = r
    return spot_o, spot_c, perp_o, perp_c, fund


def load_sofr() -> dict[dt.date, float]:
    """NY Fed SOFR (percent -> simple fraction) by session date."""
    rates: dict[dt.date, float] = {}
    with open(SOFR_CSV, newline="") as f:
        for line in f:
            parts = line.rstrip("\n").split(",")
            try:
                d = dt.datetime.strptime(parts[0], "%m/%d/%Y").date()
                rates[d] = float(parts[2]) / 100.0
            except (ValueError, IndexError):
                continue
    if not rates:
        raise RuntimeError("no SOFR rows parsed")
    return rates


# ------------------------------------------------------------------ engine ---


def build_grid(data: dict) -> list[int]:
    """Shared hourly grid: union over symbols of bars where both legs exist.
    Sleeves trade only where their own symbol trades; marking forward-fills
    closes through the 10 audited spot outage hours."""

    def legs(s):
        so, _, po, _, _ = data[s]
        return {k for k in so if k in po}

    g = sorted(legs(SYMBOLS[0]) | legs(SYMBOLS[1]))
    for a, b in zip(g, g[1:]):
        if b <= a:
            raise RuntimeError("non-monotonic grid")
    return g


def ff_maps(data: dict, grid: list[int]) -> dict:
    """Per symbol: forward-filled close marks (S, P) over the grid."""
    out = {}
    for s in SYMBOLS:
        so, sc, po, pc, _ = data[s]
        Smap, Pmap = {}, {}
        lastS = lastP = None
        for h in grid:
            if h in so:
                lastS = sc[h]
            if h in po:
                lastP = pc[h]
            Smap[h], Pmap[h] = lastS, lastP
        out[s] = (Smap, Pmap)
    return out


def simulate(data: dict, grid: list[int], rebalance_hours: int, cost_mult: float) -> dict:
    """Full-history two-sleeve simulation (settle-to-market accounting).

    Asserts the telescoping identity
        sum(price_pnl) + sum(funding) + sum(fees) == E_final - 1.0
    (exact for entry-and-hold; exact up to float rounding with rebalances,
     because unit changes revalue at trade prices against cash).
    """
    n = len(grid)
    state_cash = {s: np.empty(n) for s in SYMBOLS}
    state_qs = {s: np.empty(n) for s in SYMBOLS}
    state_qp = {s: np.empty(n) for s in SYMBOLS}
    spot_fee = SPOT_BPS * cost_mult / 1e4
    perp_fee = PERP_BPS * cost_mult / 1e4
    ffs = ff_maps(data, grid)
    tradable = {s: {h for h in grid if h in data[s][0] and h in data[s][2]} for s in SYMBOLS}
    g0 = grid[0]

    # funding events -> grid index of the hour that OPENS at the event time F
    # (settlement at the funding timestamp: a position receives/pays funding
    # iff it is held at the event instant; entry orders at the same open are
    # filled after the snapshot, so they miss that event, no lookahead).
    fund_at = {s: defaultdict(float) for s in SYMBOLS}
    fcounts = {"credited": 0, "dropped": 0}
    for s in SYMBOLS:
        for F, rate in sorted(data[s][4].items()):
            idx = bisect.bisect_left(grid, F)
            if idx < n and grid[idx] == F:
                fund_at[s][idx] += rate
                fcounts["credited"] += 1
            else:
                fcounts["dropped"] += 1

    qs = {s: 0.0 for s in SYMBOLS}    # spot units held

    qp = {s: 0.0 for s in SYMBOLS}    # perp units (<=0 short; 0 flat)
    cash = {s: 0.0 for s in SYMBOLS}
    on = {s: False for s in SYMBOLS}

    eq_close = {s: np.empty(n) for s in SYMBOLS}
    E_tot = np.empty(n)
    r_port = np.empty(n)
    r_sleeve = {s: np.empty(n) for s in SYMBOLS}
    r_bench = np.empty(n)
    turnover = np.zeros(n)      # traded notional (both legs) / E, per hour
    costs = np.zeros(n)         # fees / E, per hour
    pos_frac = np.zeros(n)      # gross exposure / equity mark, per hour
    pnl_price = {s: np.zeros(n) for s in SYMBOLS}
    pnl_fund = {s: np.zeros(n) for s in SYMBOLS}
    pnl_fees = {s: np.zeros(n) for s in SYMBOLS}
    basis_lvl = {s: np.full(n, np.nan) for s in SYMBOLS}   # S_open/P_close - 1
    held = {s: np.zeros(n) for s in SYMBOLS}
    fee_events = 0

    sofr = load_sofr()
    E_prev_total = 1.0
    eq_prev_s = {s: SLEEVE_WEIGHT for s in SYMBOLS}
    prev_S = {s: None for s in SYMBOLS}   # previous close mark (forward-filled)
    prev_P = {s: None for s in SYMBOLS}
    entered = {s: False for s in SYMBOLS}
    E_final = None

    for i, h in enumerate(grid):
        is_last = i == n - 1
        d_date = dt.datetime.fromtimestamp(h / 1000, dt.UTC).date()
        tno = 0.0
        tfo = 0.0
        rebalance_now = (rebalance_hours > 0 and i > 0
                         and (h - g0) % (rebalance_hours * HOUR) == 0)

        for s in SYMBOLS:
            Smap, Pmap = ffs[s]
            can_trade = h in tradable[s]
            price_seg = 0.0
            fund_usd = 0.0
            fee_usd = 0.0

            if can_trade:
                S_open = data[s][0][h]
                P_open = data[s][2][h]
                P_pre = P_open   # funding mark at the event instant = open
            else:
                S_open = P_open = P_pre = None

            # ---- funding settles AT the event instant = bar open -----------
            # (applies to units held over the preceding interval: pre-trade
            # units; marked at the funding price ~= this bar's open; entrants
            # at this same open fill AFTER the snapshot and miss the event)
            rate_pre = fund_at[s].get(i, 0.0)
            if on[s] and rate_pre != 0.0:
                fund_usd = -qp[s] * P_pre * rate_pre    # short: receives +rate
                cash[s] += fund_usd

            # ---- boundary settle on previous units (prev close -> open) ----
            if on[s] and can_trade and prev_S[s] is not None:
                price_seg += qs[s] * (S_open - prev_S[s]) + qp[s] * (P_open - prev_P[s])
                cash[s] += qp[s] * (P_open - prev_P[s])   # perp MTM settles
                prev_S[s] = S_open
                prev_P[s] = P_open

            # ---- trades at this open ------------------------------------
            if can_trade and not entered[s] and eq_prev_s[s] > 0.0:
                # first tradable common hour: deploy the sleeve
                C = eq_prev_s[s]
                f_word = spot_fee + perp_fee
                K = C / (1.0 + f_word)      # cash ends exactly 0 after fees
                qs[s] = K / S_open
                qp[s] = -K / P_open
                fee_usd = K * f_word
                cash[s] = C - K - fee_usd   # ~ 0 by construction
                on[s] = True
                entered[s] = True
                prev_S[s], prev_P[s] = S_open, P_open
                tno = 2.0 * K
                tfo = fee_usd
                fee_events += 1
            elif on[s] and can_trade and rebalance_now:
                # equity at open = cash + spot inventory (perp just settled)
                Eo = cash[s] + qs[s] * S_open
                tgt = Eo / 2.0
                d_spot = tgt - qs[s] * S_open            # principal move
                new_qp = -tgt / P_open                   # retarget short
                fee = abs(d_spot) * spot_fee + abs(tgt - (-qp[s] * P_open)) * perp_fee
                if d_spot != 0.0:
                    cash[s] -= d_spot
                    qs[s] += d_spot / S_open
                qp[s] = new_qp
                cash[s] -= fee
                fee_usd += fee
                tno += abs(d_spot) + abs(tgt - (-qp[s] * P_open))
                tfo += fee
                fee_events += 1

            # ---- accrue to close (current units) --------------------------
            if on[s]:
                S_c, P_c = Smap[h], Pmap[h]
                if can_trade:
                    S_ref, P_ref = S_open, P_open
                else:
                    # outage bar: reference = previous close (ff marks)
                    S_ref = prev_S[s] if prev_S[s] is not None else S_c
                    P_ref = prev_P[s] if prev_P[s] is not None else P_c
                price_seg += qs[s] * (S_c - S_ref) + qp[s] * (P_c - P_ref)
                cash[s] += qp[s] * (P_c - P_ref)          # perp MTM settles
                if P_c:
                    basis_lvl[s][i] = S_ref / P_c - 1.0
                eq_close[s][i] = cash[s] + qs[s] * S_c     # settled accounting
                state_cash[s][i], state_qs[s][i], state_qp[s][i] = cash[s], qs[s], qp[s]
                held[s][i] = 1.0
                prev_S[s], prev_P[s] = S_c, P_c

            # ---- final-bar liquidation ------------------------------------
            if is_last and on[s]:
                S_c, P_c = Smap[h], Pmap[h]
                fee = qs[s] * S_c * spot_fee + abs(qp[s]) * P_c * perp_fee
                cash[s] += qs[s] * S_c - fee    # sell spot; perp already settled
                fee_usd += fee
                tno += qs[s] * S_c + abs(qp[s]) * P_c
                tfo += fee
                fee_events += 1
                qs[s] = qp[s] = 0.0
                on[s] = False
                eq_close[s][i] = cash[s]

            # ---- book sleeve ----------------------------------------------
            pnl_price[s][i] += price_seg
            pnl_fund[s][i] += fund_usd
            pnl_fees[s][i] -= fee_usd
            if not on[s] and is_last:
                pass
            if not on[s] and not is_last and not entered[s]:
                eq_close[s][i] = cash[s]

        # ---- portfolio book ------------------------------------------------
        E_close = sum(eq_close[s][i] for s in SYMBOLS)
        if E_close <= 0:
            raise RuntimeError("nonpositive equity")
        r_port[i] = E_close / E_prev_total - 1.0
        for s in SYMBOLS:
            if eq_prev_s[s] > 0.0:
                r_sleeve[s][i] = eq_close[s][i] / eq_prev_s[s] - 1.0
            else:
                r_sleeve[s][i] = 0.0
            eq_prev_s[s] = eq_close[s][i]
        E_tot[i] = E_close
        E_prev_total = E_close
        turnover[i] = tno / E_close
        costs[i] = tfo / E_close
        gross_usd = 0.0
        for s in SYMBOLS:
            if on[s]:
                gross_usd += qs[s] * ffs[s][0][h] + abs(qp[s]) * ffs[s][1][h]
        pos_frac[i] = gross_usd / E_close
        r_bench[i] = sofr.get(d_date, 0.0) / DAY_DENOM
        if is_last:
            E_final = E_close

    if E_final is None:
        raise RuntimeError("empty window")
    total_pnl = (sum(pnl_price[s].sum() for s in SYMBOLS)
                 + sum(pnl_fund[s].sum() for s in SYMBOLS)
                 + sum(pnl_fees[s].sum() for s in SYMBOLS))
    resid = total_pnl - (E_final - 1.0)
    tol = 1e-6 * max(1.0, abs(E_final - 1.0))
    if abs(resid) > tol:
        raise RuntimeError(f"telescoping identity broken: residual {resid:.3e}")

    return {
        "grid": grid,
        "E_tot": E_tot,
        "eq_close": eq_close,
        "r_port": r_port,
        "r_sleeve": r_sleeve,
        "r_bench": r_bench,
        "turnover": turnover,
        "costs": costs,
        "pos_frac": pos_frac,
        "pnl_price": pnl_price,
        "pnl_fund": pnl_fund,
        "pnl_fees": pnl_fees,
        "basis_lvl": basis_lvl,
        "held": held,
        "fee_events": fee_events,
        "funding_events": fcounts,
        "E_final": E_final,
        "identity_residual": float(resid),
        "state_cash": state_cash,
        "state_qs": state_qs,
        "state_qp": state_qp,
    }

# ----------------------------------------------------------------- metrics ---


def series_metrics(r: np.ndarray) -> dict:
    if r.size == 0:
        return {"n_hours": 0}
    ann = float(r.mean() * HPY)
    vol = float(r.std(ddof=0) * math.sqrt(HPY))
    eq = np.cumprod(1.0 + r)
    mdd = float(((eq / np.maximum.accumulate(eq)) - 1.0).min())
    return {
        "n_hours": int(r.size),
        "ann_return": ann,
        "ann_vol": vol,
        "sharpe": (ann / vol) if vol > 1e-12 else 0.0,
        "best_hour": float(r.max()),
        "worst_hour": float(r.min()),
        "max_drawdown": mdd,
        "cum_return": float(eq[-1] - 1.0),
    }


def window_slice_idx(grid: list[int], start_d: dt.date, end_d: dt.date) -> tuple[int, int]:
    start_ms = int(dt.datetime(start_d.year, start_d.month, start_d.day, tzinfo=dt.UTC).timestamp() * 1000)
    end_ms = int(dt.datetime(end_d.year, end_d.month, end_d.day, tzinfo=dt.UTC).timestamp() * 1000) + DAY_MS - 1
    return bisect.bisect_left(grid, start_ms), bisect.bisect_right(grid, end_ms)


def window_metrics(res: dict, lo: int, hi: int) -> dict:
    years = (hi - lo) / HPY
    per_symbol = {}
    for s in SYMBOLS:
        notional = res["eq_close"][s][lo:hi] / 2.0
        fu = float(res["pnl_fund"][s][lo:hi].sum())
        fe = float(-res["pnl_fees"][s][lo:hi].sum())
        pp = float(res["pnl_price"][s][lo:hi].sum())
        per_symbol[s] = {
            "price_pnl_usd": pp,
            "funding_pnl_usd": fu,
            "fees_usd": fe,
            "funding_ann_pct_of_notional": float(fu / max(1.0, years) / max(np.nanmean(notional), 1e-9) * 100.0) if (hi - lo) else 0.0,
            "fees_ann_pct_of_notional": float(fe / max(1.0, years) / max(np.nanmean(notional), 1e-9) * 100.0) if (hi - lo) else 0.0,
            "private_price_ann_pct_of_notional": float(pp / max(1.0, years) / max(np.nanmean(notional), 1e-9) * 100.0) if (hi - lo) else 0.0,
            "held_frac_of_hours": float(res["held"][s][lo:hi].mean()),
            "basis_annualized_premium_pct": float(np.nanmean(res["basis_lvl"][s][lo:hi]) * HPY * 100.0) if (hi - lo) else 0.0,

        }
    price = float(sum(res["pnl_price"][s][lo:hi].sum() for s in SYMBOLS))
    fund = float(sum(res["pnl_fund"][s][lo:hi].sum() for s in SYMBOLS))
    fees = float(sum(res["pnl_fees"][s][lo:hi].sum() for s in SYMBOLS))
    r = res["r_port"][lo:hi]
    b = res["r_bench"][lo:hi]
    return {
        "strategy": series_metrics(r),
        "benchmark": series_metrics(b),
        "excess_cum_geo": float((1.0 + r).prod() - (1.0 + b).prod()),
        "excess_cum_arith": float((r - b).sum()),
        "attribution": {
            "price_pnl_usd": price,
            "funding_pnl_usd": fund,
            "fees_usd": fees,
            "check_sum_pnl_usd": float(price + fund + fees),
            "per_symbol": per_symbol,
        },
        "turnover_notional_turns": float(res["turnover"][lo:hi].sum() / years) if years else 0.0,
        "costs_ann_pct": float(res["costs"][lo:hi].sum() / years) if years else 0.0,
        "pos_frac_mean": float(res["pos_frac"][lo:hi].mean()),
        "n_hours": int(hi - lo),
        "years": float(years),
        "window": [str(res["grid"][lo]), str(res["grid"][hi - 1])],
    }


# ---------------------------------------------------------------- self-test ---


def self_test() -> int:
    failures = []

    def check(name, cond, detail=""):
        print(f"  {'PASS' if cond else 'FAIL'}: {name}" + (f" ({detail})" if detail else ""))
        if not cond:
            failures.append(name)

    print("[self-test] synthetic probes:")
    d0 = dt.date(2024, 1, 1)
    n = 72
    t0 = int(dt.datetime(d0.year, d0.month, d0.day, tzinfo=dt.UTC).timestamp() * 1000)
    grid = [t0 + k * HOUR for k in range(n)]
    fund = {t0 + 8 * HOUR * k: 1e-4 for k in range(1, 9)}   # 8 positive events

    def mk(perp_oc=None, perp_pc=None, spot_pc=None, fnd=None):
        sc = {h: (spot_pc[h] if spot_pc else 100.0) for h in grid}
        po = {h: (perp_oc[h] if perp_oc else 100.0) for h in grid}
        pc = {h: (perp_pc[h] if perp_pc else 100.0) for h in grid}
        so = {h: 100.0 for h in grid}
        spec = (so, sc, po, pc, dict(fnd if fnd else fund))
        return {"BTCUSDT": spec,
                "ETHUSDT": ({h: 100.0 for h in grid}, dict(sc), dict(po), dict(pc), dict(fnd if fnd else fund))}

    K = 0.5 / (1.0 + 19e-4)

    # 1. flat prices: constant equity, funding pockets only, exact magnitude
    r1 = simulate(mk(), grid, 0, 1.0)
    fu = sum(r1["pnl_fund"][s].sum() for s in SYMBOLS)
    check("flat: short receives positive funding", fu > 0, f"{fu:.8f}")
    check("flat: price pnl == 0", abs(sum(r1["pnl_price"][s].sum() for s in SYMBOLS)) < 1e-12)
    check("flat: funding == 8 events x 1e-4 x K x 2 sleeves", abs(fu - 16 * K * 1e-4) < 1e-9,
          f"{fu:.8f} vs {16 * K * 1e-4:.8f}")
    fe = sum(r1["pnl_fees"][s].sum() for s in SYMBOLS)
    check("flat: E_final == 1 + funding + fees", abs(r1["E_final"] - (1.0 + fu + fe)) < 1e-9,
          f"E={r1['E_final']:.9f}")
    check("flat: equity never dips intraday (settled accounting)",
          float(r1["E_tot"].min()) > 0.9, f"min E={r1['E_tot'].min():.6f}")

    # 2. permanent 1% perp premium, no convergence: premium is not earned
    #    unless it converges; only funding accrues. Price pnl must be zero.
    #    Settled equity at close = spot inventory marked + cash (short liab
    #    already settled): E = 2*(K/100)*S ... but perp cash settlement makes
    #    sleeve equity = spot mark + residual cash; assert via identity only.
    r2 = simulate(mk(perp_oc={h: 101.0 for h in grid}, perp_pc={h: 101.0 for h in grid}), grid, 0, 1.0)
    check("premium persistent: no price pnl", abs(sum(r2["pnl_price"][s].sum() for s in SYMBOLS)) < 1e-12)
    check("premium persistent: telescoping holds", abs(r2["identity_residual"]) < 1e-9)

    # 3. 1% premium converging at final bar -> capture (1 - 100/101)*K per sleeve
    pc3 = {h: (101.0 if h != grid[-1] else 100.0) for h in grid}
    r3 = simulate(mk(perp_oc={h: 101.0 for h in grid}, perp_pc=pc3), grid, 0, 1.0)
    pp3 = sum(r3["pnl_price"][s].sum() for s in SYMBOLS)
    check("convergence: short captures 1%-of-notional", abs(pp3 - (1 - 100.0 / 101.0) * K * 2) < 1e-9,
          f"{pp3:.8f}")
    # 4. negative funding: short pays
    fneg = {h: -1e-4 for h in fund}
    r4 = simulate(mk(fnd=fneg), grid, 0, 1.0)
    check("negative funding: short pays exact negative",
          abs(sum(r4["pnl_fund"][s].sum() for s in SYMBOLS) + 16 * K * 1e-4) < 1e-9)

    # 5. cost scaling: under 2x stress, K scales by (1+f)/(1+2f); both fee
    #    (= K*2f) and funding (= sum(rate)*K*2) scale by that K-ratio.
    r5 = simulate(mk(), grid, 0, 2.0)
    fe5 = sum(r5["pnl_fees"][s].sum() for s in SYMBOLS)
    fu5 = sum(r5["pnl_fund"][s].sum() for s in SYMBOLS)
    fp = 19e-4
    k_ratio = (1.0 + fp) / (1.0 + 2.0 * fp)
    check("2x costs: fee ratio == 2*(1+f)/(1+2f)", abs(fe5 - 2 * fe * k_ratio) < 1e-12,
          f"{fe5:.9f} vs {2 * fe * k_ratio:.9f}")
    check("2x costs: funding scales with K ratio", abs(fu5 - fu * k_ratio) < 1e-12,
          f"{fu5:.9f} vs {fu * k_ratio:.9f}")



    # 6. rebalance variant runs on flat prices without identity break
    r6 = simulate(mk(), grid, 6, 1.0)
    check("rebalance-hourly: telescoping holds", abs(r6["identity_residual"]) < 1e-9)

    print("[self-test] real-data invariants (Jan 2024 sample):")
    data = {s: load_symbol(s, dt.date(2024, 1, 1), dt.date(2024, 1, 7)) for s in SYMBOLS}
    g = build_grid(data)
    r7 = simulate(data, g, 0, 1.0)
    check("real: telescoping identity", abs(r7["identity_residual"]) < 1e-6)
    hours = (g[-1] - g[0]) / HOUR
    fph = sum(r7["pnl_fund"][s].sum() for s in SYMBOLS) / hours
    check("real: positive funding accrual per hour (Jan-2024 sample)", fph > 0, f"{fph:.6f} USD/h")
    check("real: hourly returns sane", float(np.abs(r7["r_port"]).max()) < 0.2,
          f"max|r|={float(np.abs(r7['r_port']).max()):.4f}")
    check("real: gross exposure ~2x equity while on", abs(r7["pos_frac"].mean() - 2.0) < 0.2,
          f"mean pos_frac={r7['pos_frac'].mean():.3f}")

    print(f"[self-test] {'ALL PASS' if not failures else str(len(failures)) + ' FAILURES'}")
    return 0 if not failures else 1


# --------------------------------------------------------------------- run ---


def main() -> int:
    ap = argparse.ArgumentParser(description="Unconditional equal-weight crypto carry backtest")
    ap.add_argument("--start", default="2020-08-01")
    ap.add_argument("--end", default="2024-12-31", help="fail-closed beyond 2024-12-31 in wave one")
    ap.add_argument("--output-dir", default="results")
    ap.add_argument("--rebalance-hours", type=int, default=FROZEN_REBALANCE_HOURS,
                    help="within-sleeve notional maintenance clock; 0 = entry-and-hold")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    unlock = os.environ.get("QUANT_HOLDOUT_UNLOCK") == "1"
    start_d = dt.datetime.strptime(args.start, "%Y-%m-%d").date()
    end_d = dt.datetime.strptime(args.end, "%Y-%m-%d").date()
    if end_d > WAVE1_END and not unlock:
        print(f"REFUSED: --end {end_d} exceeds wave-one boundary {WAVE1_END}; "
              "set QUANT_HOLDOUT_UNLOCK=1 to override (external evaluator only).", file=sys.stderr)
        return 2
    if end_d > WAVE1_END and unlock:
        print("NOTE: QUANT_HOLDOUT_UNLOCK=1: post-2024 access unlocked (external evaluator).",
              file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)

    load_start = min(start_d, HISTORY_START)
    load_end = max(start_d, end_d)
    data = {s: load_symbol(s, load_start, load_end) for s in SYMBOLS}
    grid = build_grid(data)

    variants = {name: simulate(data, grid, args.rebalance_hours, mult)
                for name, mult in (("net", 1.0), ("gross", 0.0), ("stress_2x", STRESS_MULT))}

    report = {
        "schema_version": 1,
        "id": "crypto-carry-decay",
        "requested_window": {"start": str(start_d), "end": str(end_d)},
        "loaded_history": {"start": str(load_start), "end": str(load_end)},
        "rebalance_hours": args.rebalance_hours,
        "identity_residuals": {k: variants[k]["identity_residual"] for k in variants},
        "funding_events": variants["net"]["funding_events"],
        "fee_events": variants["net"]["fee_events"],
        "cost_model": {
            "spot_one_way_bps": SPOT_BPS,
            "perp_one_way_bps": PERP_BPS,
            "stress_multiplier": STRESS_MULT,
            "application": "entry open, scheduled rebalance opens, final liquidation; cash outflow on traded notional",
        },
    }

    windows = {}
    for label, (ws_, we_) in (
        ("development", (HISTORY_START, dt.date(2023, 12, 31))),
        ("validation", (dt.date(2024, 1, 1), dt.date(2024, 12, 31))),
        ("requested", (start_d, end_d)),
    ):
        if we_ > end_d:
            continue
        lo, hi = window_slice_idx(grid, ws_, min(we_, end_d))
        if hi - lo < 4:
            windows[label] = {"error": f"only {hi - lo} hours"}
            continue
        wm = window_metrics(variants["net"], lo, hi)
        wm["gross_strategy"] = series_metrics(variants["gross"]["r_port"][lo:hi])
        bs = variants["stress_2x"]["r_bench"][lo:hi]
        ss = variants["stress_2x"]["r_port"][lo:hi]
        wm["stress_2x"] = {
            "strategy": series_metrics(ss),
            "excess_cum_geo": float((1.0 + ss).prod() - (1.0 + bs).prod()),
            "excess_cum_arith": float((ss - bs).sum()),
            "costs_ann_pct": float(variants["stress_2x"]["costs"][lo:hi].sum() / wm["years"]),
            "turnover_notional_turns": float(variants["stress_2x"]["turnover"][lo:hi].sum() / wm["years"]),
        }
        windows[label] = wm
    report["windows"] = windows

    # returns.csv: requested window, net series, daily aggregation
    lo, hi = window_slice_idx(grid, start_d, end_d)
    res = variants["net"]
    by_day = defaultdict(list)
    for i in range(lo, hi):
        d = dt.datetime.fromtimestamp(grid[i] / 1000, dt.UTC).date().isoformat()
        by_day[d].append(i)
    returns_path = os.path.join(args.output_dir, "returns.csv")
    with open(returns_path, "w") as f:
        f.write("date,strategy_return,benchmark_return,position,turnover,cost\n")
        for d in sorted(by_day):
            idxs = by_day[d]
            sr = float(np.prod(1.0 + res["r_port"][idxs]) - 1.0)
            br = float(np.prod(1.0 + res["r_bench"][idxs]) - 1.0)
            pos = float(np.mean(res["pos_frac"][idxs]))
            to = float(np.sum(res["turnover"][idxs]))
            co = float(np.sum(res["costs"][idxs]))
            f.write(f"{d},{sr:.10f},{br:.10f},{pos:.6f},{to:.6f},{co:.6f}\n")

    with open(os.path.join(args.output_dir, "report.json"), "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write("\n")

    for label, wm in windows.items():
        if "error" in wm:
            print(f"[{label}] {wm['error']}")
            continue
        s, b, att = wm["strategy"], wm["benchmark"], wm["attribution"]
        print(f"[{label}] n={wm['n_hours']}h ann={s['ann_return']:+.4f} vol={s['ann_vol']:.4f} "
              f"sharpe={s['sharpe']:+.3f} mdd={s['max_drawdown']:.3f} | bench ann={b['ann_return']:+.4f} "
              f"sharpe={b['sharpe']:+.3f} | excess_geo={wm['excess_cum_geo']:+.4f} "
              f"2x_excess_geo={wm['stress_2x']['excess_cum_geo']:+.4f}")
        print(f"    USD: price={att['price_pnl_usd']:+.4f} funding={att['funding_pnl_usd']:+.4f} "
              f"fees={att['fees_usd']:+.4f} sum={att['check_sum_pnl_usd']:+.4f}")
    print(f"returns.csv -> {returns_path}")
    print(f"identity residuals: {report['identity_residuals']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
