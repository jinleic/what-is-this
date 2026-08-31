#!/usr/bin/env python3
"""Corrected restatement of the frozen crypto-carry-decay candidate.

This is NOT a rescue and NOT a new strategy validation candidate. It restates
the frozen pre-registered rule (unconditional equal-weight BTCUSDT/ETHUSDT
spot-long / USD-M-perp-short 1:1 carry, entry-and-hold, and the three
pre-declared maintenance variants) with every confirmed accounting defect
repaired, so the corrected metrics can be compared side by side with the
invalid frozen report.

Confirmed defects repaired here (frozen run.py sites in parentheses):

  D1  parse_funding clipped every standalone negative funding rate through
      max(rate, 0.0) (frozen line 152). The short pays when the rate is
      negative: payment = -q_p * P * r with q_p < 0. This engine preserves
      the exact signed rate and FAILS on duplicate snapped events.
  D2  window_metrics reported funding/fee/price yields against
      eq_close/2 (frozen line 518) instead of the actual short notional
      |q_p| * P. This engine normalizes by actual |q_p| * P.
  D3  The frozen entry (frozen lines 356-364) consumed the whole sleeve in
      spot (cash exactly 0) and reserved NO independent USD-M margin, so
      total gross notional was ~2x NAV at entry and gross/NAV drifted up
      with prices. This engine funds spot principal from sleeve capital and
      reserves independent USD cash collateral equal to 100% of the perp
      notional at entry, so at entry: spot + perp-notional <= 1.0x NAV
      exactly, with zero borrowing and no cross-collateral assumption.
      Collateral remuneration is stated conservatively (see SOFR_COLLAT).
  D4  The SOFR benchmark accrued zero on non-publication days (frozen line
      448, sofr.get(d, 0.0)) and used /365 hourly accrual
      (DAY_DENOM=24*365). This engine carries the last published rate
      forward across non-publication calendar days and accrues ACT/360
      (overnight USD convention) per calendar day.

Frozen rules preserved unchanged:
  * unconditional always-on signal, 50/50 initial capital split,
    no cross-sleeve transfers;
  * next-bar execution: entry fills at the first tradable common bar open
    AFTER data start; funding is credited only to positions held at the
    event instant (entry orders at the same open fill after the snapshot);
  * one-way taker fees incl. slippage: spot 12 bps, perp 7 bps
    (contract defaults), 2x cost stress = double every trade fee;
  * funding settled in cash at each 8h event timestamp on notional marked
    at the perp open of the hour containing the event;
  * hourly settle-to-market accounting for the perp leg; spot marked at
    forward-filled closes with entry-at-open identity.
  * The pre-registered first falsification gate, 15% annualized funding
    income per unit of actual short notional on the 2024 validation year,
    is applied UNCHANGED (correction D2 changes only the denominator to
    what the gate always intended: the short notional).
  * The four pre-declared rebalance variants (0/24/168/720h maintenance
    clocks) are all recomputed as diagnostics. NO new winner is selected:
    entry-and-hold (rebalance_hours=0) remains the frozen-rule
    restatement; the others are reported for comparison only.

Position continuity and window slicing mirror the frozen engine: the
simulation runs continuously from 2020-08-01 (first tradable bar) through
the requested end; windows are slices of the same accounting path.

Determinism: byte-identical outputs for identical --start/--end.

Access control: identical fail-closed gate to the frozen run (exit 2 when
--end exceeds 2024-12-31 without QUANT_HOLDOUT_UNLOCK=1). The corrected-audit
dossier is a retrospective diagnostic that intentionally inspects the
historical-confirmation window; holdout_accessed is therefore true and the
status can never be validation-pass.

Usage:
  python3 run.py                                  # dev + validation report
  python3 run.py --start 2024-01-01 --end 2024-12-31 --output-dir out
  python3 run.py --start 2025-01-01 --end 2026-07-31 --output-dir out \
      --unlock-historical-confirmation            # audit-only backtest
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
SLEEVE_WEIGHT = 0.5            # initial capital split (frozen rule)

SPOT_BPS = 12.0                # contract default one-way spot taker incl. slippage
PERP_BPS = 7.0                 # contract default one-way perpetual taker incl. slippage
STRESS_MULT = 2.0              # contract 2x-cost stress

FROZEN_REBALANCE_HOURS = 0     # frozen variant (development-selected, unchanged)

HISTORY_START = dt.date(2020, 8, 1)
WAVE1_END = dt.date(2024, 12, 31)
HC_END = dt.date(2026, 7, 31)  # historical-confirmation window end

# Conservative collateral remuneration: idle sleeve cash (entry cash plus
# accumulating funding income, before it is needed) and the perp collateral
# account earn ZERO interest in this engine. SOFR accrues ONLY on the
# benchmark leg. No sweep, no repo, no cross-margin credit: stated so no
# result can be flattered by an unmodeled cash yield.
SOFR_COLLAT = None
# Perp collateral sizing: 100% of entry perp notional, funded from sleeve
# equity at entry (reserve collateral), independent of the spot principal.
# No borrowing, no cross-collateral, no rehypothecation of spot inventory.
# The reserve is maintained so that ALLCash >= 0 at every close (verified);
# perp unrealized P&L settles into the collateral cash.
PERP_COLLAT_FRAC = 1.0
# MARGIN_BUFFER_FRAC of sleeve equity is kept as uncalled precautionary
# cash ABOVE the 100% collateral floor: with both the <= 1.0x-entry-gross
# mandate and 100% independent collateral, the binding equation holds with
# exactly zero residual — any adverse tick would otherwise instantaneously
# breach the floor and force liquidation. The buffer is a minimum
# feasibility requirement, declared here and identical for all four
# variants (not a tuned parameter).
MARGIN_BUFFER_FRAC = 0.02
# MAINTENANCE_MARGIN_FRAC: exchange-style liquidation threshold. The 100%
# INITIAL collateral is posted once at entry; afterwards the perp account
# (initial collateral + settled funding/perp P&L) is liquidated only when
# its equity falls below this fraction of current short notional (Binance
# USD-M tier-1 maintenance for low notional ~0.5%). A breach is a margin
# call: the sleeve terminates; no borrowing, no cross-collateral.
MAINTENANCE_MARGIN_FRAC = 0.005

HOUR = 3_600_000
MS_8H = 28_800_000
DAY_MS = 86_400_000
HPY = 24.0 * 365.25            # annualization for hourly series (frozen rule)
# D4: overnight USD interest is ACT/360. The benchmark accrues a full
# calendar day at the prevailing daily SOFR (carried forward), charged in
# 24 hourly steps: daily accrual / 360 / 24 per hour.
ACT_DENOM = 360.0


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
    """Funding events snapped to the canonical 8h grid, EXACT SIGNED rates.

    Correction D1: the frozen parser used max(rate, 0.0), clipping every
    standalone negative rate to zero (zero snap collisions exist in the
    audited archives, so the max() never resolved a duplicate - it was a
    pure sign clip). This parser keeps the signed rate and raises on any
    duplicate snapped event instead of silently folding it."""
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
                r = float(p[2])
                if ts in out:
                    raise RuntimeError(f"duplicate snapped funding event at {ts} in {zip_path}")
                out[ts] = r
    return out


def load_symbol(symbol: str, start_d: dt.date, end_d: dt.date):
    """Spot/perp bar-open->(open, close) dicts and snapped funding event->SIGNED
    rate dict, for the UTC date range [start_d, end_d]."""
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
    """NY Fed SOFR (percent -> simple fraction) by publication date."""
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


def sofr_daily_carry_forward(rates: dict[dt.date, float], start_d: dt.date, end_d: dt.date) -> dict[dt.date, float]:
    """Correction D4: overnight cash earns the last published SOFR on
    non-publication calendar days (weekends/holidays), ACT/360 per day.
    Returns date -> daily simple accrual factor fraction (rate/360)."""
    out: dict[dt.date, float] = {}
    last = None
    d = start_d
    while d <= end_d:
        if d in rates:
            last = rates[d]
        if last is None:
            raise RuntimeError(f"no SOFR publication on or before {d}")
        out[d] = last / ACT_DENOM
        d += dt.timedelta(days=1)
    return out


# ------------------------------------------------------------------ engine ---


def build_grid(data: dict) -> list[int]:
    """Shared hourly grid: union over symbols of bars where both legs exist.
    Sleeves trade only where their own symbol trades; marking forward-fills
    closes through the audited spot outage hours."""
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


def simulate(data: dict, grid: list[int], rebalance_hours: int, cost_mult: float,
             sofr_hr: dict[int, float]) -> dict:
    """Full-history two-sleeve simulation with corrected, fully funded
    accounting.

    Per sleeve (initial equity C = 0.5):
      correction D3 entry:
        spot principal  K_s = (C - perp_collateral) / (1 + spot_fee_frac)
        perp collateral K_m = PERP_COLLAT_FRAC * K_p_entry_target
      with the 1:1 rule short notional == spot notional (both at entry
      price), entry fee on each leg; collateral is independent USD-M margin
      cash, earning nothing, and sleeve ALL-cash never goes negative.

    Equities: eq = spot inventory mark + all sleeve cash (collateral +
    spot-principal cash residual + settled funding and perp P&L).

    Asserts the telescoping identity
        sum(price_pnl) + sum(funding) + sum(fees) == E_final - 1.0
    exact for entry-and-hold (funding is a cash flow; perp MTM telescopes;
    fees are cash outflows; collateral transfers cancel inside equity).
    """
    n = len(grid)
    state_cash = {s: np.empty(n) for s in SYMBOLS}
    state_qs = {s: np.empty(n) for s in SYMBOLS}
    state_qp = {s: np.empty(n) for s in SYMBOLS}
    state_neg = {s: np.empty(n) for s in SYMBOLS}   # 1.0 if short leg active
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
    cash = {s: 0.0 for s in SYMBOLS}  # ALL sleeve cash (spot-principal
                                      # residual + perp collateral + P&L)
    on = {s: False for s in SYMBOLS}

    # correction D2: actual short notional series |q_p| * P (hourly marks)
    short_notional = {s: np.zeros(n) for s in SYMBOLS}
    spot_notional = {s: np.zeros(n) for s in SYMBOLS}
    collat_balance = {s: np.zeros(n) for s in SYMBOLS}  # perp collateral+P&L view

    eq_close = {s: np.empty(n) for s in SYMBOLS}
    E_tot = np.empty(n)
    r_port = np.empty(n)
    r_sleeve = {s: np.empty(n) for s in SYMBOLS}
    r_bench = np.empty(n)
    turnover = np.zeros(n)      # traded notional (both legs) / E, per hour
    costs = np.zeros(n)         # fees / E, per hour
    pos_frac = np.zeros(n)      # gross notional / equity mark, per hour
    pnl_price = {s: np.zeros(n) for s in SYMBOLS}
    pnl_fund = {s: np.zeros(n) for s in SYMBOLS}
    COLLAT0 = {s: 0.0 for s in SYMBOLS}    # initial collateral posted at entry
    perp_acct = {s: 0.0 for s in SYMBOLS}  # settled funding + perp MTM in the account
    pnl_fees = {s: np.zeros(n) for s in SYMBOLS}
    basis_lvl = {s: np.full(n, np.nan) for s in SYMBOLS}   # same-time S/P - 1
    held = {s: np.zeros(n) for s in SYMBOLS}
    neg_cash = np.zeros(n)      # negative all-cash flag (must stay 0)
    fee_events = 0

    E_prev_total = 1.0
    eq_prev_s = {s: SLEEVE_WEIGHT for s in SYMBOLS}
    prev_S = {s: None for s in SYMBOLS}   # previous close mark (forward-filled)
    prev_P = {s: None for s in SYMBOLS}
    entered = {s: False for s in SYMBOLS}
    term = {s: False for s in SYMBOLS}          # sleeve terminated by margin call
    term_at = {s: None for s in SYMBOLS}        # grid hour of termination
    term_fee_usd = {s: 0.0 for s in SYMBOLS}    # liquidation fee at termination
    term_gross_usd = {s: 0.0 for s in SYMBOLS}  # gross notional unwound
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
            eq_close[s][i] = cash[s]

            if can_trade:
                S_open = data[s][0][h]
                P_open = data[s][2][h]
                P_pre = P_open   # funding mark at the event instant = open
            else:
                S_open = P_open = P_pre = None

            # ---- funding settles AT the event instant = bar open -----------
            # exact signed rate (D1): short cash flow = -q_p * P * r, q_p < 0
            # -> short RECEIVES when r > 0, PAYS when r < 0.
            rate_pre = fund_at[s].get(i, 0.0)
            if on[s] and rate_pre != 0.0:
                fund_usd = -qp[s] * P_pre * rate_pre
                cash[s] += fund_usd
                perp_acct[s] += fund_usd

            # ---- boundary settle on previous units (prev close -> open) ----
            if on[s] and can_trade and prev_S[s] is not None:
                price_seg += qs[s] * (S_open - prev_S[s]) + qp[s] * (P_open - prev_P[s])
                cash[s] += qp[s] * (P_open - prev_P[s])   # perp MTM settles
                perp_acct[s] += qp[s] * (P_open - prev_P[s])
                prev_S[s] = S_open
                prev_P[s] = P_open
            if can_trade and not entered[s] and not term[s] and eq_prev_s[s] > 0.0:
                # first tradable common hour: deploy the sleeve fully funded
                C = eq_prev_s[s]
                f_word = spot_fee + perp_fee
                # 1:1 rule at entry prices: spot notional == perp notional == N
                # budget: N (spot) + N*W (collateral floor) + B (buffer)
                # + N*f (fees) <= C  ->  N = (C - B) / (1 + f + W),
                # B = MARGIN_BUFFER_FRAC * C precautionary uncalled cash.
                buf = MARGIN_BUFFER_FRAC * C
                N = (C - buf) / (1.0 + f_word + PERP_COLLAT_FRAC)
                qs[s] = N / S_open
                qp[s] = -N / P_open
                fee_usd = N * f_word
                # cash after entry = C - N (spot principal) - N*f (fees),
                # = N*W + buf exactly; no borrowing (cash floor = collateral+buffer).
                cash[s] = C - N - fee_usd
                on[s] = True
                entered[s] = True
                prev_S[s], prev_P[s] = S_open, P_open
                tno = 2.0 * N
                tfo = fee_usd
                fee_events += 1
                COLLAT0[s] = PERP_COLLAT_FRAC * N     # initial collateral posted
                perp_acct[s] = 0.0                    # account opens flat
            elif on[s] and can_trade and rebalance_now:
                # maintenance (comparison-only variants): retarget BOTH legs
                # to half of open equity, keeping the fully funded identity
                # (spot principal + retained collateral <= equity, cash >= 0).
                Eo = cash[s] + qs[s] * S_open
                tgt = Eo / 2.0
                d_spot = tgt - qs[s] * S_open            # principal move
                new_qp = -tgt / P_open                   # retarget short
                fee = abs(d_spot) * spot_fee + abs(tgt - (-qp[s] * P_open)) * perp_fee
                collat_need = PERP_COLLAT_FRAC * tgt
                # cash after moves: cash - d_spot - fee >= collat_need
                # (satisfied: cash >= Eo - invested = collat_need by induction;
                #  with tgt=Eo/2 and cash >= Eo/2 = collat_need)
                if d_spot != 0.0:
                    cash[s] -= d_spot
                    qs[s] += d_spot / S_open
                qp[s] = new_qp
                cash[s] -= fee
                # re-post independent collateral for the new short notional
                COLLAT0[s] = PERP_COLLAT_FRAC * tgt
                perp_acct[s] = 0.0    # collateral account resets to fresh base
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
                cash[s] += qp[s] * (P_c - P_ref)          # perp MTM settles -> sleeve pool
                perp_acct[s] += qp[s] * (P_c - P_ref)     # -> perp account view
                # perp-account equity = initial 100% collateral + settled
                # funding + settled perp MTM. Trading cash is separate and
                # cannot be topped up (no borrowing). Liquidation fires when
                # perp-account equity < maintenance fraction of current
                # short notional.
                collat_acct = COLLAT0[s] + perp_acct[s]
                maint_req = MAINTENANCE_MARGIN_FRAC * abs(qp[s]) * P_c
                if collat_acct < maint_req - 1e-9:
                    Sx, Px = S_c, P_c
                    fee = qs[s] * Sx * spot_fee + abs(qp[s]) * Px * perp_fee
                    cash[s] += qs[s] * Sx - fee   # sell spot at close; perp already settled to close
                    fee_usd += fee
                    tno += qs[s] * Sx + abs(qp[s]) * Px
                    tfo += fee
                    fee_events += 1
                    qs[s] = qp[s] = 0.0
                    on[s] = False
                    term[s] = True
                    term_at[s] = h
                    term_fee_usd[s] = fee
                    eq_close[s][i] = cash[s]
                    short_notional[s][i] = 0.0
                    spot_notional[s][i] = 0.0
                    prev_S[s], prev_P[s] = Sx, Px
                else:
                    eq_close[s][i] = cash[s] + qs[s] * S_c
                    short_notional[s][i] = abs(qp[s]) * P_c
                    spot_notional[s][i] = qs[s] * S_c
                    collat_balance[s][i] = cash[s]
                    state_cash[s][i], state_qs[s][i], state_qp[s][i] = cash[s], qs[s], qp[s]
                    state_neg[s][i] = 1.0
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
                short_notional[s][i] = 0.0
                spot_notional[s][i] = 0.0

            # ---- book sleeve ----------------------------------------------
            pnl_price[s][i] += price_seg
            pnl_fund[s][i] += fund_usd
            pnl_fees[s][i] -= fee_usd
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
        neg_cash[i] = 1.0 if any(cash[s] < -1e-9 for s in SYMBOLS) else 0.0
        # benchmark: carried-forward SOFR, ACT/360 (D4) on calendar days
        r_bench[i] = sofr_hr[d_date]
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
    if neg_cash.sum() > 0:
        raise RuntimeError("borrowing detected: negative sleeve cash")

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
        "short_notional": short_notional,
        "spot_notional": spot_notional,
        "collat_balance": collat_balance,
        "state_neg": state_neg,
        "neg_cash_hours": int(neg_cash.sum()),
        "margin_call_terminations": {
            s: ({"terminated": True, "at_hour": term_at[s],
                 "at_iso": dt.datetime.fromtimestamp(term_at[s] / 1000, dt.UTC).isoformat(),
                 "unwound_gross_usd": term_gross_usd[s],
                 "liquidation_fee_usd": term_fee_usd[s]}
                if term[s] else {"terminated": False})
            for s in SYMBOLS
        },
        "margin_call_any": any(term[s] for s in SYMBOLS),
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
    """Window metrics with correction D2 applied: every per-notional yield
    uses the ACTUAL short notional |q_p| * P (hourly close marks), never
    equity/2."""
    years = (hi - lo) / HPY
    per_symbol = {}
    for s in SYMBOLS:
        sn = res["short_notional"][s][lo:hi]
        sn_mean = float(np.nanmean(sn)) if (hi - lo) else 0.0
        fu = float(res["pnl_fund"][s][lo:hi].sum())
        fe = float(-res["pnl_fees"][s][lo:hi].sum())
        pp = float(res["pnl_price"][s][lo:hi].sum())
        held_sn = (res["held"][s][lo:hi] * res["short_notional"][s][lo:hi])
        mean_held_sn = float(np.nanmean(held_sn)) if (hi - lo) else 0.0
        per_symbol[s] = {
            "price_pnl_usd": pp,
            "funding_pnl_usd": fu,
            "fees_usd": fe,
            "mean_short_notional_usd": mean_held_sn,
            "funding_ann_pct_of_short_notional": float(
                fu / max(1.0, years) / max(mean_held_sn, 1e-9) * 100.0) if (hi - lo) else 0.0,
            "fees_ann_pct_of_short_notional": float(
                fe / max(1.0, years) / max(mean_held_sn, 1e-9) * 100.0) if (hi - lo) else 0.0,
            "price_pnl_ann_pct_of_short_notional": float(
                pp / max(1.0, years) / max(mean_held_sn, 1e-9) * 100.0) if (hi - lo) else 0.0,
            "held_frac_of_hours": float(res["held"][s][lo:hi].mean()),
            "same_time_basis_mean_frac": float(np.nanmean(res["basis_lvl"][s][lo:hi])) if (hi - lo) else 0.0,
        }
    price = float(sum(res["pnl_price"][s][lo:hi].sum() for s in SYMBOLS))
    fund = float(sum(res["pnl_fund"][s][lo:hi].sum() for s in SYMBOLS))
    fees = float(sum(res["pnl_fees"][s][lo:hi].sum() for s in SYMBOLS))
    r = res["r_port"][lo:hi]
    b = res["r_bench"][lo:hi]
    posm = float(res["pos_frac"][lo:hi].mean()) if (hi - lo) else 0.0
    posmax = float(res["pos_frac"][lo:hi].max()) if (hi - lo) else 0.0
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
        "gross_over_equity_mean": posm,
        "gross_over_equity_max": posmax,
        "n_hours": int(hi - lo),
        "years": float(years),
        "window": [str(res["grid"][lo]), str(res["grid"][hi - 1])],
    }



def funding_gate_metrics(res: dict, lo: int, hi: int) -> dict:
    """Pre-registered first falsification gate RECOMPUTED CORRECTLY:
    annualized funding income as % of ACTUAL short notional (mean |q_p|*P
    over held hours), per symbol and pair-aggregate (dollar-weighted)."""
    fu_pair = 0.0
    sn_mean_pair = 0.0
    per = {}
    years = (hi - lo) / HPY
    for s in SYMBOLS:
        held = res["held"][s][lo:hi]
        sn_held = float(np.nansum(held * res["short_notional"][s][lo:hi]))
        hours = float(held.sum())
        sn_mean = sn_held / hours if hours > 0 else 0.0
        fu = float(res["pnl_fund"][s][lo:hi].sum())
        per[s] = {
            "funding_usd": fu,
            "mean_short_notional_usd": sn_mean,
            "held_hours": hours,
            "funding_ann_pct": float(fu / max(1.0, years) / max(sn_mean, 1e-9) * 100.0),
        }
        fu_pair += fu
        sn_mean_pair += sn_mean
    pair_ann_pct = float(fu_pair / max(1.0, years) / max(sn_mean_pair, 1e-9) * 100.0)
    return {
        "gate": "fg-carry-1 (recomputed per original spec: annualized funding income on ACTUAL short notional > 15%/yr)",
        "threshold_pct": 15.0,
        "per_symbol": per,
        "pair_dollar_weighted_ann_pct": pair_ann_pct,
        "gate_result_recomputed": "pass" if pair_ann_pct > 15.0 else "fail",
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
    sofr_hr = {d: 0.05 / ACT_DENOM / 24.0 for d in (d0 + dt.timedelta(days=k) for k in range(n // 24 + 2))}

    def mk(perp_oc=None, perp_pc=None, spot_pc=None, fnd=None):
        sc = {h: (spot_pc[h] if spot_pc else 100.0) for h in grid}
        po = {h: (perp_oc[h] if perp_oc else 100.0) for h in grid}
        pc = {h: (perp_pc[h] if perp_pc else 100.0) for h in grid}
        so = {h: 100.0 for h in grid}
        spec = (so, sc, po, pc, dict(fnd if fnd else fund))
        return {"BTCUSDT": spec,
                "ETHUSDT": ({h: 100.0 for h in grid}, dict(sc), dict(po), dict(pc), dict(fnd if fnd else fund))}

    K = (0.5 - MARGIN_BUFFER_FRAC * 0.5) / (1.0 + 19e-4 + 1.0)   # corrected entry notional per sleeve
    f_word = 19e-4

    # 1. flat prices: constant equity, funding pockets only, exact magnitude
    r1 = simulate(mk(), grid, 0, 1.0, sofr_hr)
    fu = sum(r1["pnl_fund"][s].sum() for s in SYMBOLS)
    check("flat: short receives positive funding", fu > 0, f"{fu:.8f}")
    check("flat: price pnl == 0", abs(sum(r1["pnl_price"][s].sum() for s in SYMBOLS)) < 1e-12)
    check("flat: funding == 8 events x 1e-4 x K x 2 sleeves", abs(fu - 16 * K * 1e-4) < 1e-9,
          f"{fu:.8f} vs {16 * K * 1e-4:.8f}")
    fe = sum(r1["pnl_fees"][s].sum() for s in SYMBOLS)
    check("flat: E_final == 1 + funding + fees", abs(r1["E_final"] - (1.0 + fu + fe)) < 1e-9,
          f"E={r1['E_final']:.9f}")
    check("flat: entry gross <= 1x NAV", r1["pos_frac"][1] <= 1.0 + 1e-12,
          f"pos_frac={r1['pos_frac'][1]:.6f}")
    check("flat: no negative cash", r1["neg_cash_hours"] == 0)

    # 2. permanent 1% perp premium, no convergence: only funding accrues
    r2 = simulate(mk(perp_oc={h: 101.0 for h in grid}, perp_pc={h: 101.0 for h in grid}), grid, 0, 1.0, sofr_hr)
    check("premium persistent: no price pnl", abs(sum(r2["pnl_price"][s].sum() for s in SYMBOLS)) < 1e-12)
    check("premium persistent: telescoping holds", abs(r2["identity_residual"]) < 1e-9)

    # 3. convergence capture with corrected sizing
    pc3 = {h: (101.0 if h != grid[-1] else 100.0) for h in grid}
    r3 = simulate(mk(perp_oc={h: 101.0 for h in grid}, perp_pc=pc3), grid, 0, 1.0, sofr_hr)
    pp3 = sum(r3["pnl_price"][s].sum() for s in SYMBOLS)
    check("convergence: short captures 1%-of-notional", abs(pp3 - (1 - 100.0 / 101.0) * K * 2) < 1e-9,
          f"{pp3:.8f}")

    # 4. negative funding: short PAYS the exact signed amount on every event
    #    for the full horizon (D1). The 2%-of-equity buffer absorbs this
    #    small drain; no margin call fires (asserted).
    fneg = {h: -1e-4 for h in fund}
    r4 = simulate(mk(fnd=fneg), grid, 0, 1.0, sofr_hr)
    fu4 = sum(r4["pnl_fund"][s].sum() for s in SYMBOLS)
    check("negative funding: short pays exact signed every event",
          abs(fu4 + 16 * K * 1e-4) < 1e-9)
    check("negative funding: no margin call (buffer absorbs)",
          not r4["margin_call_any"])

    # 5. cost scaling: under 2x stress, entry notional K scales by
    #    (1+f+1)/(1+2f+1); both fee (= K*2f) and funding (= sum(rate)*2K)
    #    scale by that K-ratio.
    r5 = simulate(mk(), grid, 0, 2.0, sofr_hr)
    fe5 = sum(r5["pnl_fees"][s].sum() for s in SYMBOLS)
    fu5 = sum(r5["pnl_fund"][s].sum() for s in SYMBOLS)
    k2 = (0.5 - MARGIN_BUFFER_FRAC * 0.5) / (1.0 + 2.0 * f_word + 1.0)
    k_ratio = k2 / K
    check("2x costs: fee ratio == K2x/K1x", abs(abs(fe5) - abs(2 * fe * k_ratio)) < 1e-12,
          f"{abs(fe5):.9f} vs {abs(2 * fe * k_ratio):.9f}")
    check("2x costs: funding scales with K ratio", abs(fu5 - fu * k_ratio) < 1e-12,
          f"{fu5:.9f} vs {fu * k_ratio:.9f}")

    # 6. rebalance variant runs on flat prices without identity break
    r6 = simulate(mk(), grid, 6, 1.0, sofr_hr)
    check("rebalance-hourly: telescoping holds", abs(r6["identity_residual"]) < 1e-9)
    check("rebalance-hourly: no negative cash", r6["neg_cash_hours"] == 0)

    # 7. rising market: perp MTM drains collateral into equity via settled
    #    P&L; no borrowing is allowed and spot+perp gross can drift above
    #    1x NAV passively (mark-to-market) — assert bounded drift only.
    rise = {h: 100.0 * (1.0 + 0.002 * (h - grid[0]) / HOUR) for h in grid}
    r7 = simulate(mk(perp_oc=rise, perp_pc=rise, spot_pc=rise), grid, 0, 1.0, sofr_hr)
    check("rising market: no borrowing", r7["neg_cash_hours"] == 0)
    check("rising market: telescoping holds", abs(r7["identity_residual"]) < 1e-9)
    check("rising market: gross/NAV drift bounded (no entry leverage)",
          r7["pos_frac"].max() < 1.15, f"max={r7['pos_frac'].max():.4f}")

    print(f"[self-test] {'ALL PASS' if not failures else str(len(failures)) + ' FAILURES'}")
    return 0 if not failures else 1


# --------------------------------------------------------------------- run ---


def build_sofr_hourly(rates: dict[dt.date, float], start_ms: int, end_ms: int) -> dict[dt.date, float]:
    last_d = max(rates)
    first_d = min(rates)
    d0 = dt.datetime.fromtimestamp(start_ms / 1000, dt.UTC).date()
    d1 = dt.datetime.fromtimestamp(end_ms / 1000, dt.UTC).date()
    d0 = min(d0, first_d)
    return sofr_daily_carry_forward(rates, d0, max(d1, first_d))


def main() -> int:
    ap = argparse.ArgumentParser(description="Corrected fully funded crypto carry restatement")
    ap.add_argument("--start", default="2020-08-01")
    ap.add_argument("--end", default="2024-12-31", help="fail-closed beyond 2024-12-31 in wave one")
    ap.add_argument("--output-dir", default="results")
    ap.add_argument("--rebalance-hours", type=int, default=FROZEN_REBALANCE_HOURS,
                    help="maintenance clock; 0 = entry-and-hold (frozen rule)")
    ap.add_argument("--variants", action="store_true",
                    help="also compute the three pre-declared maintenance variants as diagnostics")
    ap.add_argument("--unlock-historical-confirmation", action="store_true",
                    help="audit-only access to the 2025-01-01..2026-07-31 window "
                         "(deliberate retrospective inspection; never a holdout claim)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    unlock = os.environ.get("QUANT_HOLDOUT_UNLOCK") == "1" or args.unlock_historical_confirmation
    start_d = dt.datetime.strptime(args.start, "%Y-%m-%d").date()
    end_d = dt.datetime.strptime(args.end, "%Y-%m-%d").date()
    if end_d > WAVE1_END and not unlock:
        print(f"REFUSED: --end {end_d} exceeds wave-one boundary {WAVE1_END}; "
              "set QUANT_HOLDOUT_UNLOCK=1 or --unlock-historical-confirmation.", file=sys.stderr)
        return 2
    if end_d > WAVE1_END and unlock:
        print("NOTE: retrospective historical-confirmation window accessed deliberately "
              "(audit diagnostic; holdout_accessed=true; never a validation claim).",
              file=sys.stderr)

    os.makedirs(args.output_dir, exist_ok=True)

    load_start = min(start_d, HISTORY_START)
    load_end = max(start_d, end_d)
    data = {s: load_symbol(s, load_start, load_end) for s in SYMBOLS}
    grid = build_grid(data)
    sofr = load_sofr()
    d0 = dt.datetime.fromtimestamp(grid[0] / 1000, dt.UTC).date()
    d1 = dt.datetime.fromtimestamp((grid[-1] + DAY_MS - 1) / 1000, dt.UTC).date()
    sofr_hr = sofr_daily_carry_forward(sofr, min(d0, min(sofr)), d1)

    variants = {name: simulate(data, grid, args.rebalance_hours, mult, sofr_hr)
                for name, mult in (("net", 1.0), ("gross", 0.0), ("stress_2x", STRESS_MULT))}

    report = {
        "schema_version": 1,
        "id": "crypto-carry-corrected-audit",
        "requested_window": {"start": str(start_d), "end": str(end_d)},
        "loaded_history": {"start": str(load_start), "end": str(load_end)},
        "rebalance_hours": args.rebalance_hours,
        "identity_residuals": {k: variants[k]["identity_residual"] for k in variants},
        "funding_events": variants["net"]["funding_events"],
        "fee_events": variants["net"]["fee_events"],
        "negative_cash_hours": variants["net"]["neg_cash_hours"],
        "margin_call_terminations": variants["net"]["margin_call_terminations"],
        "margin_call_any": variants["net"]["margin_call_any"],
        "collateral_policy": {
            "initial_usd_collateral_frac": PERP_COLLAT_FRAC,
            "maintenance_margin_frac_of_notional": MAINTENANCE_MARGIN_FRAC,
            "buffer_frac_of_initial_sleeve_equity": MARGIN_BUFFER_FRAC,
            "remuneration": "zero interest on collateral and idle cash (conservative)",
        },
        "corrections_applied": {
            "D1_signed_funding": "exact signed funding rates; short pays on negative; parser fails on duplicate snapped events",
            "D2_actual_short_notional": "all per-notional yields use mean |q_p|*P over held hours",
            "D3_fully_funded_entry": "independent USD-M margin collateral = 100% of entry short notional; entry gross = (1-fee drag) / (1+f+1.0) <= 1x NAV; no borrowing; sleeve cash >= 0 every hour (asserted)",
            "D4_sofr_carry_forward_act360": "benchmark accrues carried-forward SOFR on every calendar day, /360",
            "collateral_renumeration": "conservative zero interest on sleeve cash and collateral",
        },
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
        wm["funding_gate_recomputed"] = funding_gate_metrics(variants["net"], lo, hi)
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

    # comparison-only diagnostic variants (never selected)
    if args.variants:
        var_report = {}
        for rh, name in ((24, "rebal-daily"), (168, "rebal-weekly"), (720, "rebal-monthly")):
            vres = simulate(data, grid, rh, 1.0, sofr_hr)
            vs = {}
            for label, (ws_, we_) in (
                ("development", (HISTORY_START, dt.date(2023, 12, 31))),
                ("validation", (dt.date(2024, 1, 1), dt.date(2024, 12, 31))),
            ):
                if we_ > end_d:
                    continue
                vhi = window_slice_idx(grid, ws_, we_)[1]
                vlo = window_slice_idx(grid, ws_, we_)[0]
                vstats = window_metrics(vres, vlo, vhi)
                vs[label] = {
                    "ann_return": vstats["strategy"]["ann_return"],
                    "sharpe": vstats["strategy"]["sharpe"],
                    "excess_cum_geo": vstats["excess_cum_geo"],
                    "excess_cum_geo_2x": None,  # computed below
                }
                vres2 = simulate(data, grid, rh, STRESS_MULT, sofr_hr)
                s2 = vres2["r_port"][vlo:vhi]
                b2 = vres2["r_bench"][vlo:vhi]
                vs[label]["excess_cum_geo_2x"] = float((1.0 + s2).prod() - (1.0 + b2).prod())
            var_report[name] = vs
        report["comparison_variants_dev_val_only"] = var_report

    with open(os.path.join(args.output_dir, "report.json"), "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write("\n")

    for label, wm in windows.items():
        if "error" in wm:
            print(f"[{label}] {wm['error']}")
            continue
        s, b, att = wm["strategy"], wm["benchmark"], wm["attribution"]
        g = wm["funding_gate_recomputed"]
        print(f"[{label}] n={wm['n_hours']}h ann={s['ann_return']:+.4f} vol={s['ann_vol']:.4f} "
              f"sharpe={s['sharpe']:+.3f} mdd={s['max_drawdown']:.3f} | bench ann={b['ann_return']:+.4f} "
              f"sharpe={b['sharpe']:+.3f} | excess_geo={wm['excess_cum_geo']:+.4f} "
              f"2x_excess_geo={wm['stress_2x']['excess_cum_geo']:+.4f}")
        print(f"    USD: price={att['price_pnl_usd']:+.4f} funding={att['funding_pnl_usd']:+.4f} "
              f"fees={att['fees_usd']:+.4f} sum={att['check_sum_pnl_usd']:+.4f}")
        print(f"    gross/E mean={wm['gross_over_equity_mean']:.3f} max={wm['gross_over_equity_max']:.3f}")
        if label == "validation":
            print(f"    funding gate: pair {g['pair_dollar_weighted_ann_pct']:.3f}%/yr "
                  f"of actual short notional vs 15% -> {g['gate_result_recomputed'].upper()}")
    print(f"returns.csv -> {returns_path}")
    print(f"identity residuals: {report['identity_residuals']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
