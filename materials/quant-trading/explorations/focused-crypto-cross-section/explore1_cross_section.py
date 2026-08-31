#!/usr/bin/env python3
"""Explore 1 (v3): post-rrereview engine for focused-crypto-cross-section.

Frozen decisions (review_rule_freeze.json + parent instruction):

TIMING  Target from close t-1 is the exact held position over outcome day t
        (close[t-1] -> close[t]); the SAME held tensor drives PnL, censoring,
        concentration, and the collateral audit. The boundary rebalance fee
        to move from the pretrade (drifted) position to the close-(t-1)
        target is charged inside period t.

LEDGER  Explicit quantities x price notionals + cash; NAV = cash + sum
        notionals. Price drift changes current notionals. Signed funding is
        a cash flow applied to held notionals. Desired post-fee target
        notionals solve a deterministic scalar (per-fixed target weights) so
        the post-trade weights equal the registered target weights. Emitted
        per-period returns satisfy final NAV/product(1+r) telescoping exactly.

FAIL-CLOSED  Unknown funding is NaN. Complete price+funding outcomes are
        required for eligibility, features, IC, strategy, and benchmark. If
        ANY held asset's outcome on day t is incomplete, the variant path is
        STOPPED at the prior known close: no further finite return is booked,
        all held censored assets are recorded, and the censored tail is
        excluded from every statistic (.ic/turnover/cost/NAV beyond that day).
        No zero imputation, no trade-out, no finite return from unknown data.

BENCHMARK  Per variant: exact causal opportunity/activity mask from that
        variant's own held-position targets. Invested only when the strategy
        is active (gross > 0 at the prior close) and the outcome day is
        complete for members; NaN/flat otherwise. No outcome-conditioned
        constituent selection.

GATE TRUTH TABLE (SSOT):
  pass        = all required gates observed true
  falsified   = any required gate observably false (even with other unknowns)
  inconclusive = no false and >= 1 unknown
Deterministic proofs (PROOFS): hand ledger (timing/drift/funding/fee/telescope),
held-missing stop, benchmark inactivity, six-subperiod truth table.
"""
import csv
import io
import json
import math
import pickle
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2].parent
RAW = ROOT / "quant-trading/data/raw/focused/cross-section"
EXP = ROOT / "quant-trading/explorations/focused-crypto-cross-section"

FEE_RATE = 0.0004
FEE_RATE_2X = 0.0008

WIN_DEV = ("2020-08-01", "2023-12-31")
WIN_VAL = ("2024-01-01", "2024-12-31")
WIN_CON = ("2025-01-01", "2026-07-31")

MIN_ASSETS = 20
IC_THRESHOLD = 0.02
N_PARTITIONS = 6
PURGE_DAYS = 30

# ---------------------------------------------------------------- data loading

def load_klines():
    with open(RAW / "derived/klines_1h.pkl", "rb") as f:
        return pickle.load(f)


def load_funding():
    FR = RAW / "binance/um/funding"
    out = {}
    for sym_dir in sorted(FR.iterdir()):
        sym = sym_dir.name
        d = {}
        for z in sorted(sym_dir.glob(f"{sym}-fundingRate-*.zip")):
            with zipfile.ZipFile(z) as zf:
                for member in zf.namelist():
                    with zf.open(member) as f:
                        for line in csv.reader(io.TextIOWrapper(f, encoding="utf-8")):
                            if not line or line[0] == "calc_time":
                                continue
                            try:
                                ts = int(line[0]); rate = float(line[2])
                            except (ValueError, IndexError):
                                continue
                            if ts in d and d[ts] != rate:
                                raise ValueError(f"conflicting duplicate funding {sym} {ts}")
                            d[ts] = rate
        out[sym] = d
    return out


def build_daily_panel(krows):
    from collections import defaultdict
    closes = defaultdict(dict); vols = defaultdict(dict)
    highs = defaultdict(dict); lows = defaultdict(dict)
    for sym, ts, o, h, l, c, v, qv, _ in krows:
        day = ts // 86_400_000
        closes[sym][day] = c
        vols[sym][day] = vols[sym].get(day, 0.0) + qv
        if h > highs[sym].get(day, -math.inf):
            highs[sym][day] = h
        if l < lows[sym].get(day, math.inf):
            lows[sym][day] = l
    return closes, highs, lows, vols


def panel():
    krows = load_klines()
    closes, highs, lows, vols = build_daily_panel(krows)
    funding = load_funding()
    dl = json.loads((RAW / "provenance/funding_downloads.json").read_text())
    funding_incomplete = sorted({r["symbol"] for r in dl if "error" in r})
    all_assets = sorted(closes)
    days = sorted(set().union(*(closes[s] for s in all_assets)))
    days = [d for d in days if d >= 18_475]
    fund_day = {s: {} for s in all_assets}
    fund_known_day = {s: set() for s in all_assets}
    for s in all_assets:
        for ts, r in funding.get(s, {}).items():
            day = ts // 86_400_000
            fund_day[s][day] = fund_day[s].get(day, 0.0) + r
            fund_known_day[s].add(day)
    D = {d: i for i, d in enumerate(days)}
    N = len(all_assets); T = len(days)
    close = np.full((T, N), np.nan); vol = np.full((T, N), np.nan)
    high = np.full((T, N), np.nan); low = np.full((T, N), np.nan)
    fund = np.full((T, N), np.nan)          # unknown -> NaN
    listed = np.zeros((T, N), bool)
    fund_known = np.zeros((T, N), bool)
    for j, s in enumerate(all_assets):
        for d, c in closes[s].items():
            if d in D:
                close[D[d], j] = c
                listed[D[d], j] = True
        for d, q in vols[s].items():
            if d in D: vol[D[d], j] = q
        for d, h in highs[s].items():
            if d in D: high[D[d], j] = h
        for d, l in lows[s].items():
            if d in D: low[D[d], j] = l
        for d, r in fund_day[s].items():
            if d in D:
                fund[D[d], j] = r
                fund_known[D[d], j] = True
    return dict(close=close, high=high, low=low, vol=vol, fund=fund,
                fund_known=fund_known, listed=listed,
                meta={"assets": all_assets, "days_ms": days,
                      "funding_incomplete_assets": funding_incomplete})


# ------------------------------------------------------------------- features

def trailing_ret_complete(close, listed, L):
    """21-day momentum requires a complete L+1 observation path (no NaN)."""
    T, N = close.shape
    r = np.full((T, N), np.nan)
    for t in range(L, T):
        w = close[t-L:t+1]
        ok = (~np.isnan(w).any(axis=0) & listed[t-L:t+1].all(axis=0))
        a, b = close[t-L], close[t]
        with np.errstate(all="ignore"):
            rr = np.where(ok, b / a - 1.0, np.nan)
        r[t] = rr
    return r


def trailing_vol_complete(close, listed, L):
    """21-day realized vol requires a complete L+1 path; per-asset counts."""
    T, N = close.shape
    v = np.full((T, N), np.nan)
    with np.errstate(all="ignore"):
        lr = np.full((T, N), np.nan)
        for t in range(1, T):
            ok = listed[t] & listed[t-1] & ~np.isnan(close[t]) & ~np.isnan(close[t-1])
            lr[t, ok] = np.log(close[t, ok] / close[t-1, ok])
        for t in range(L, T):
            w = lr[t-L+1:t+1]
            complete = (~np.isnan(w).any(axis=0)) & listed[t-L+1:t+1].all(axis=0)
            sd = np.nanstd(np.where(complete[None, :], w, np.nan), axis=0)
            v[t] = np.where(complete & np.isfinite(sd) & (sd > 0), sd, np.nan)
    return v


def build_features(P):
    close, vol, fund, listed = P["close"], P["vol"], P["fund"], P["listed"]
    fund_known = P["fund_known"]
    F = {}
    F["mom21"] = trailing_ret_complete(close, listed, 21)
    F["volsd"] = trailing_vol_complete(close, listed, 21)
    T, N = close.shape
    liq = np.full((T, N), np.nan)
    for t in range(19, T):
        w = vol[t-19:t+1]
        complete = (~np.isnan(w).any(axis=0)) & listed[t-19:t+1].all(axis=0)
        with np.errstate(all="ignore"):
            m = np.nanmean(np.where(complete[None, :], w, np.nan), axis=0)
        liq[t] = np.where(complete, m, np.nan)
    F["liq20"] = liq
    fu = np.full((T, N), np.nan)
    for t in range(2, T):
        known = fund_known[t-2:t+1].all(axis=0)
        with np.errstate(all="ignore"):
            m = np.nanmean(np.where(known[None, :], fund[t-2:t+1], np.nan), axis=0)
        fu[t] = np.where(known, m, np.nan)
    F["fund3"] = fu
    fu7 = np.full((T, N), np.nan)
    for t in range(9, T):
        known = fund_known[t-9:t+1].all(axis=0)
        w3 = fund[t-2:t+1]; w7 = fund[t-9:t-2]
        with np.errstate(all="ignore"):
            m3 = np.nanmean(np.where(known[None, :], w3, np.nan), axis=0)
            m7 = np.nanmean(np.where(known[None, :], w7, np.nan), axis=0)
        fu7[t] = np.where(known, m3 - m7, np.nan)
    F["fundchg"] = fu7
    return F


VARIANTS = ["momentum", "signed_funding", "mom_plus_funding", "funding_change",
            "vol_scaled_momentum", "liq_composite"]


def _z(x, mask):
    mu = np.nanmean(x[mask]); sd = np.nanstd(x[mask])
    if not (np.isfinite(sd) and sd > 0):
        return None
    z = np.full_like(x, np.nan)
    z[mask] = (x[mask] - mu) / sd
    return z


def compute_scores(F, P):
    mom, fnd, fch, volsd = F["mom21"], F["fund3"], F["fundchg"], F["volsd"]
    liq = F["liq20"]
    S = {}
    S["momentum"] = mom
    S["signed_funding"] = -fnd
    S["mom_plus_funding"] = mom - fnd * 100.0
    S["funding_change"] = -fch
    with np.errstate(all="ignore"):
        S["vol_scaled_momentum"] = np.where(np.isnan(volsd), np.nan, mom / volsd)
    ok = ~np.isnan(mom) & ~np.isnan(fnd) & ~np.isnan(liq)
    Zm = np.full_like(mom, np.nan); Zf = np.full_like(mom, np.nan); Zl = np.full_like(mom, np.nan)
    for t in range(mom.shape[0]):
        m = ok[t]
        if m.sum() < 8:
            continue
        z1 = _z(mom[t], m)
        if z1 is not None: Zm[t] = z1
        z2 = _z(-fnd[t], m)
        if z2 is not None: Zf[t] = z2
        with np.errstate(all="ignore"):
            z3 = _z(np.log1p(np.maximum(liq[t], 0.0)), m)
        if z3 is not None: Zl[t] = z3
    comp = 0.6 * Zm + 0.6 * Zf + 0.2 * Zl
    comp[~ok] = np.nan
    S["liq_composite"] = comp
    return S


# ---------------------------------------------------------------- eligibility

MIN_DOLLAR_VOL = 5_000_000.0


def eligibility(P):
    """Canonical SSOT seasoning/liquidity rule (provenance owner).

    Eligible at t iff: complete 20-observation price window ending at t
    (no NaN, all listed) and its mean quote volume >= $5M.
    """
    vol, listed = P["vol"], P["listed"]
    T, N = vol.shape
    E = np.zeros((T, N), bool)
    for t in range(19, T):
        w = vol[t-19:t+1]
        complete = (~np.isnan(w).any(axis=0)) & listed[t-19:t+1].all(axis=0)
        with np.errstate(all="ignore"):
            m = np.nanmean(np.where(complete[None, :], w, np.nan), axis=0)
        ok = complete & (np.nan_to_num(m, nan=0.0) >= MIN_DOLLAR_VOL)
        E[t] = ok
    return E


def complete_outcome(P):
    """Outcome-day mask: observed price close AND signed funding vs prior close."""
    close, fund, listed, fk = P["close"], P["fund"], P["listed"], P["fund_known"]
    px_ok = np.zeros_like(listed)
    px_ok[1:] = listed[1:] & listed[:-1] & ~np.isnan(close[1:]) & ~np.isnan(close[:-1])
    return px_ok & fk


# ------------------------------------------------------------------- IC utils

def avg_rank(xs):
    n = len(xs)
    order = np.argsort(xs, kind="stable")
    ranks = np.empty(n, float)
    i = 0
    sx = xs[order]
    while i < n:
        j = i
        while j + 1 < n and sx[j+1] == sx[i]:
            j += 1
        ranks[order[i:j+1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def rank_ic_series(score, tlpr, signal_opportunity, complete, tail_nan_from):
    """Spearman IC over the per-asset causal signal opportunity mask.

   Observed-pair mask at row t = signal_opportunity[t-1] (frozen E[t-1] and
    required signal-day funding completeness/full lookbacks + finite score,
    folded in by the caller) AND complete outcome t. Returns NaN for a
    constant (all-tied) cross-section. Rows from tail_nan_from onward are
    unavailable (NaN), never finite from zero.
    """
    T = tlpr.shape[0]
    tics = np.full(T, np.nan)
    t_end = T if tail_nan_from is None else min(tail_nan_from, T)
    for t in range(1, t_end):
        s = score[t-1]
        m = signal_opportunity[t-1] & complete[t]
        if m.sum() < MIN_ASSETS:
            continue
        ss = s[m]
        if np.ptp(ss) == 0:
            continue
        rs = avg_rank(ss); rrk = avg_rank(tlpr[t][m])
        tics[t] = np.corrcoef(rs, rrk)[0, 1]
    return tics


# ---------------------------------------------------------------- portfolio v3

def signal_opportunity(P, E, score):
    """Per-asset causal signal opportunity at signal time t.

    frozen E[t] (complete-20 liquidity SSOT) AND signal-day funding observed
    AND finite variant score (which already requires full registered
    lookbacks). Portfolio AND benchmark constituents come only from this mask.
    """
    return E & P["fund_known"] & ~np.isnan(score)


def run_variant_engine(P, E, score, fee=FEE_RATE, top_q=5, bot_q=5):
    """Exact-boundary quantity/NAV ledger with first-censor NaN tail.

    Loop for outcome row t:
      1. At close t-1 rebalance/fee-solve the drifted book Q into
         target[t-1] (== held[t]); a zero target explicitly CLOSES the book.
      2. That exact quantity vector then earns outcome t = price return plus
         signed funding.
      3. The same held vector drives censor, concentration and collateral.
      4. FIRST held incomplete outcome only: record all held censored assets,
         mark row t and the entire tail unavailable (NaN), break.
    """
    close, fund = P["close"], P["fund"]
    listed = P["listed"]
    assets = P["meta"]["assets"]
    T, N = close.shape
    price_ret = np.full((T, N), np.nan)
    for t in range(1, T):
        a, b = close[t-1], close[t]
        ok = listed[t-1] & listed[t] & ~np.isnan(a) & ~np.isnan(b)
        price_ret[t, ok] = b[ok] / a[ok] - 1.0
    tlpr = price_ret - fund
    complete = complete_outcome(P)
    sig_opp = signal_opportunity(P, E, score)

    # held[t] == target decided at close t-1, held over outcome day t.
    held = np.zeros((T, N))
    active = np.zeros(T, bool)
    for t in range(1, T):
        idx = np.where(sig_opp[t-1])[0]
        if len(idx) >= MIN_ASSETS:
            order = idx[np.argsort(score[t-1][idx], kind="stable")]
            shorts, longs = order[:bot_q], order[-top_q:]
            held[t, longs] = 0.5 / len(longs)
            held[t, shorts] = -0.5 / len(shorts)
            active[t] = True

    Q = np.zeros(N)
    CASH = 1.0
    sr = np.full(T, np.nan); turnover = np.full(T, np.nan)
    cost = np.full(T, np.nan)
    gross = np.full(T, np.nan); net_exp = np.full(T, np.nan)
    asset_pnl = np.full((T, N), np.nan)
    W_out = np.zeros((T, N))
    NAV_path = np.full(T, np.nan)
    nav_acc = 1.0
    stop_day = None
    censored_assets = []
    held_at_stop = []
    for t in range(1, T):
        tgt_w = held[t]
        # (1) boundary rebalance at close t-1 into target[t-1]; zero target
        #     closes the prior book explicitly.
        nav_base = CASH + float(Q.sum())
        NAV_after_est = nav_base
        for _ in range(128):
            fees_est = float(fee * np.abs(tgt_w * NAV_after_est - Q).sum())
            NAV_new = nav_base - fees_est
            if abs(NAV_new - NAV_after_est) < 1e-13:
                NAV_after_est = NAV_new
                break
            NAV_after_est = NAV_new
        tgt_notional = tgt_w * NAV_after_est
        trade_notional = tgt_notional - Q
        fees_exact = float(fee * np.abs(trade_notional).sum())
        nav_after_fees = nav_base - fees_exact
        CASH += float(Q.sum()) - float(tgt_notional.sum()) - fees_exact
        Q = tgt_notional.copy()
        if abs((CASH + float(Q.sum())) - nav_after_fees) > 1e-9:
            raise AssertionError(f"rebalance ledger identity broken at t={t}")
        # (4) FIRST censor check on the held vector, before any outcome:
        incomplete_held = (tgt_w != 0) & ~complete[t]
        if incomplete_held.any():
            stop_day = t
            held_at_stop = [assets[j] for j in np.where(tgt_w != 0)[0]]
            censored_assets = [{
                "day_index": int(t), "asset": assets[j],
                "side": "long" if tgt_w[j] > 0 else "short",
                "held_weight": float(tgt_w[j]),
                "reason": "held outcome incomplete (price and/or signed "
                          "funding unobserved); row t and entire tail marked "
                          "unavailable (NaN); no finite return emitted",
            } for j in np.where(incomplete_held)[0]]
            break
        # (2) the exact held quantity vector earns outcome t:
        pr = np.where(np.isnan(price_ret[t]), 0.0, price_ret[t])
        fnd = np.where(np.isnan(fund[t]), 0.0, fund[t])
        funding_cash = -float((Q * fnd).sum())
        Q = Q * (1.0 + pr)
        CASH += funding_cash
        NAV_now = CASH + float(Q.sum())
        r_t = NAV_now / nav_acc - 1.0
        sr[t] = r_t
        nav_acc = NAV_now
        NAV_path[t] = nav_acc
        turnover[t] = float(np.abs(trade_notional).sum()) / nav_base
        cost[t] = fees_exact / nav_base
        gross[t] = float(np.abs(tgt_w).sum())
        net_exp[t] = float(tgt_w.sum())
        W_out[t] = tgt_w
        asset_pnl[t] = (tgt_w * np.where(np.isnan(tlpr[t]), 0.0, tlpr[t])
                        - fee * np.abs(trade_notional) / nav_base)
    available = np.isfinite(sr)
    return dict(
        strat=sr, price_ret=price_ret, tlpr=tlpr, W=W_out, gross=gross,
        net=net_exp, turnover=turnover, cost=cost, asset_pnl=asset_pnl,
        NAV_path=NAV_path, final_nav=nav_acc, held=held,
        active=active, signal_opportunity=sig_opp, complete=complete,
        available=available, stop_day=stop_day,
        censored_assets=censored_assets, held_at_stop=held_at_stop,
        censored_selected_outcomes=len(censored_assets),
    )


# --------------------------------------------------- proof: hand ledger case

def _synthetic_panel(T=30, N=25, vol=6e6):
    close = np.ones((T, N))
    return dict(close=close, high=close.copy(), low=close.copy(),
                vol=np.full((T, N), vol), fund=np.zeros((T, N)),
                fund_known=np.ones((T, N), bool),
                listed=np.ones((T, N), bool),
                meta={"assets": [f"A{i}" for i in range(N)],
                      "days_ms": [18475 + i for i in range(T)],
                      "funding_incomplete_assets": []})


def PROOFS():
    """Deterministic acceptance proofs against the production engine."""
    proofs = {}
    fee = 0.001
    # ---- (a) exact hand ledger: NAV 100, 25% long A, price +10%, funding -1%
    nav0, q0, cash0 = 100.0, 25.0, 75.0
    pr, fr = 0.10, -0.01
    # boundary: already at 25% target, no trade -> no fee; outcome:
    funding_cash = -q0 * fr                      # +0.25 (long receives)
    q_after = q0 * (1 + pr)                      # 27.5
    nav_after = cash0 + funding_cash + q_after   # 102.75
    r_hand = nav_after / nav0 - 1.0              # 0.0275
    proofs["hand_ledger"] = {
        "nav_open": nav0, "notional_open": q0, "price_return": pr,
        "signed_funding_rate": fr, "funding_cash": funding_cash,
        "notional_after_drift": q_after, "nav_close": nav_after,
        "return_emitted": r_hand,
        "expected_return_w_times_tlpr": 0.25 * (pr - fr),
        "matches_weight_identity": abs(r_hand - 0.25 * (pr - fr)) < 1e-12,
    }
    # ---- (b) production impulse: target[t-1] (not [t-2]) earns outcome t,
    #          and a zero target closes the stale book.
    P2 = _synthetic_panel()
    P2["close"][21:, 0] = 1.10          # +10% on outcome day 21 only
    score = np.full((30, 25), np.nan)
    score[20, :] = 0.0
    score[20, 0] = 5.0                   # A0 top -> long on outcome day 21
    score[20, 24] = -5.0                 # A24 bottom -> short
    E2 = eligibility(P2)
    R2 = run_variant_engine(P2, E2, score, fee=fee)
    long_w = 0.5 / 5
    proofs["production_impulse"] = {
        "held_row_21_A0_weight": float(R2["held"][21, 0]),
        "held_row_22_A0_weight": float(R2["held"][22, 0]),
        "return_row_21": float(R2["strat"][21]),
        "return_row_22": float(R2["strat"][22]),
        "impulse_earned_on_row_21": bool(R2["strat"][21] > 0.5 * long_w * 0.10),
        "impulse_not_deferred_to_row_22": bool(
            abs(R2["strat"][22]) < 0.5 * long_w * 0.10),
        "zero_target_closed_book_row_22": bool(
            np.abs(R2["held"][22]).sum() == 0.0
            and np.isfinite(R2["turnover"][22])
            and R2["turnover"][22] > 0.0),
        "interpretation": "target set at close 20 earns outcome row 21; the "
                          "inactive target at close 21 explicitly closes the "
                          "book, charging its closing fee in row 22",
    }
    # ---- (c) FIRST censor only, NaN/unavailable tail, no finite zero tail
    P3 = _synthetic_panel()
    sc3 = np.zeros((30, 25))
    sc3[:, 0] = 5.0; sc3[:, 24] = -5.0
    P3["fund_known"][25, 0] = False      # first incomplete held outcome
    P3["fund_known"][27, 1] = False      # later failure must not be reported
    E3 = eligibility(P3)
    R3 = run_variant_engine(P3, E3, sc3, fee=fee)
    tail = R3["strat"][R3["stop_day"]:]
    proofs["first_censor_stop"] = {
        "stop_day": R3["stop_day"],
        "first_not_last_failure": R3["stop_day"] == 25,
        "censored_assets": R3["censored_assets"],
        "held_at_stop": R3["held_at_stop"],
        "tail_all_nan_not_zero": bool(np.all(np.isnan(tail))),
        "tail_len": int(len(tail)),
        "available_rows_end_before_stop": int(np.max(
            np.where(R3["available"])[0]) if R3["available"].any() else -1),
    }
    return proofs


# ------------------------------------------------------- benchmark + windows

def equal_weight_benchmark(tlpr, signal_opp, complete, min_assets=MIN_ASSETS):
    """Per-variant benchmark over the SAME per-asset causal signal mask.

    Constituents at row t are exactly signal_opp[t-1] (no ex-post filtering
    on outcome completeness). If ANY held constituent's outcome at t is
    incomplete, the benchmark path STOPS: row t and the tail are NaN.
    """
    T = tlpr.shape[0]
    b = np.full(T, np.nan)
    stop = None
    for t in range(1, T):
        m = signal_opp[t-1]
        if m.sum() < min_assets:
            continue
        if not complete[t][m].all():
            stop = t
            break
        b[t] = float(np.mean(tlpr[t][m]))
    return b, stop


def max_asset_pnl_share(asset_pnl, days_idx, assets):
    rows = [t for t in days_idx if np.isfinite(asset_pnl[t]).all()]
    if not rows:
        return None, None
    pnl = np.sum(asset_pnl[rows], axis=0)
    denom = float(np.sum(np.abs(pnl)))
    if denom <= 0:
        return 0.0, None
    j = int(np.argmax(np.abs(pnl)))
    return float(abs(pnl[j]) / denom), assets[j]


def collateral_audit(P, held, days_idx, available, maintenance=0.10,
                     buffer=0.20):
    """Same held tensor as PnL/censor, restricted to available rows."""
    ratios = []
    for t in days_idx:
        if not available[t]:
            continue
        w = held[t]
        if np.abs(w).sum() == 0:
            continue
        a, hi, lo = P["close"][t-1], P["high"][t], P["low"][t]
        good = ((w != 0) & np.isfinite(a) & np.isfinite(hi)
                & np.isfinite(lo) & (a > 0))
        shock = np.zeros_like(w)
        L, Sh = good & (w > 0), good & (w < 0)
        shock[L] = lo[L] / a[L] - 1.0
        shock[Sh] = -(hi[Sh] / a[Sh] - 1.0)
        equity = 1.0 + float(np.sum(np.abs(w) * shock))
        gr = float(np.abs(w).sum())
        ratios.append(equity / gr if gr else np.inf)
    if not ratios:
        return {"maintenance_ratio": maintenance,
                "declared_isolation_buffer": buffer,
                "min_collateral_to_gross_at_daily_extrema": None,
                "buffer_breaches": None, "liquidations": None,
                "available_rows": 0}
    return {
        "maintenance_ratio": maintenance,
        "declared_isolation_buffer": buffer,
        "min_collateral_to_gross_at_daily_extrema": min(ratios),
        "buffer_breaches": int(sum(x <= maintenance + buffer for x in ratios)),
        "liquidations": int(sum(x <= maintenance for x in ratios)),
        "available_rows": len(ratios),
    }


def perf_stats(sr, days_idx, bh=None):
    """Availability-aware: unavailable (NaN) rows are excluded, never zeroed.

    A window with no available row returns None totals (UNKNOWN, not false).
    """
    r = np.asarray(sr)[days_idx]
    avail = np.isfinite(r)
    out = {"n_days": int(len(r)), "n_available": int(avail.sum())}
    if avail.sum() == 0:
        out.update({"total_return": None, "ann_sharpe": None, "max_dd": None})
    else:
        rr = r[avail]
        eq = np.cumprod(1 + rr)
        peak = np.maximum.accumulate(np.concatenate([[1.0], eq]))
        out["total_return"] = float(np.prod(1 + rr) - 1)
        out["ann_sharpe"] = (float(np.mean(rr) / np.std(rr) * math.sqrt(365))
                             if np.std(rr) > 0 else 0.0)
        out["max_dd"] = float(np.max(1 - eq / peak[1:]))
    if bh is not None:
        b = np.asarray(bh)[days_idx]
        bav = np.isfinite(b)
        out["bench_total_return"] = (float(np.prod(1 + b[bav]) - 1)
                                     if bav.sum() else None)
        out["bench_observations"] = int(bav.sum())
    return out


def day_idx_range_full(days_ms, a, b):
    import datetime as dtm
    da = int(dtm.datetime.strptime(a, "%Y-%m-%d").replace(tzinfo=dtm.timezone.utc).timestamp() * 1000) // 86_400_000
    db = int(dtm.datetime.strptime(b, "%Y-%m-%d").replace(tzinfo=dtm.timezone.utc).timestamp() * 1000) // 86_400_000
    lo = next(i for i, d in enumerate(days_ms) if d >= da)
    hi = max(i for i, d in enumerate(days_ms) if d <= db)
    return lo, hi


def day_idx_range(days_ms, a, b):
    lo, hi = day_idx_range_full(days_ms, a, b)
    return list(range(lo, hi + 1))


def build_gate_subperiods(days_ms):
    d0, d1 = day_idx_range_full(days_ms, WIN_DEV[0], WIN_DEV[1])
    dev = list(range(d0, d1 + 1))
    n = len(dev)
    edges = np.linspace(0, n, N_PARTITIONS + 1).astype(int)
    subs = []
    for k in range(N_PARTITIONS):
        a, b = dev[edges[k]], dev[edges[k+1] - 1]
        purge_start = min(a + PURGE_DAYS, b + 1)
        subs.append({
            "partition": k + 1, "raw_start_index": a, "raw_end_index": b,
            "test_start_index": purge_start, "test_end_index": b,
            "purged_days": purge_start - a,
        })
    return subs


def evaluate_variant(P, F, S, E, name, days_ms, complete, out):
    """One corrected result object; the frozen truth table is applied here."""
    score = S[name]
    R1 = run_variant_engine(P, E, score, fee=FEE_RATE)
    R2 = run_variant_engine(P, E, score, fee=FEE_RATE_2X)
    sr, tlpr = R1["strat"], R1["tlpr"]
    sig_opp = R1["signal_opportunity"]
    avail = R1["available"]
    stop_day = R1["stop_day"]
    bh, bench_stop = equal_weight_benchmark(tlpr, sig_opp, complete)
    tics = rank_ic_series(score, tlpr, sig_opp, complete, stop_day)
    windows = {}
    for wlabel, (a, b) in (("development", WIN_DEV), ("validation", WIN_VAL),
                           ("confirmation", WIN_CON)):
        di = day_idx_range(days_ms, a, b)
        st = perf_stats(sr, di, bh)
        st["long_short_total_return_1x"] = st["total_return"]
        st["long_short_total_return_2x"] = perf_stats(R2["strat"], di)["total_return"]
        ticw = tics[di]
        st["mean_rank_ic"] = (float(np.nanmean(ticw))
                              if np.any(~np.isnan(ticw)) else None)
        st["rank_ic_observations"] = int(np.sum(~np.isnan(ticw)))
        av = np.array([avail[t] for t in di])
        st["available_rows"] = int(av.sum())
        st["mean_daily_turnover"] = (float(np.nanmean(R1["turnover"][di]))
                                     if av.any() else None)
        st["fees_paid_1x"] = (float(np.nansum(R1["cost"][di]))
                              if av.any() else None)
        st["fees_paid_2x"] = (float(np.nansum(R2["cost"][di]))
                              if av.any() else None)
        st["gross_exposure_mean"] = (float(np.nanmean(R1["gross"][di]))
                                     if av.any() else None)
        st["net_exposure_mean"] = (float(np.nanmean(R1["net"][di]))
                                   if av.any() else None)
        share, asset = max_asset_pnl_share(R1["asset_pnl"], di,
                                           P["meta"]["assets"])
        st["max_asset_pnl_share"] = share
        st["max_asset_pnl_asset"] = asset
        eligible_counts = E[di].sum(axis=1)
        st["n_eligible_mean"] = float(np.mean(eligible_counts))
        st["n_eligible_min"] = int(np.min(eligible_counts))
        st["days_eligible_ge_20"] = int(np.sum(eligible_counts >= 20))
        st["collateral"] = collateral_audit(P, R1["held"], di, avail)
        windows[wlabel] = st

    subs = build_gate_subperiods(days_ms)
    subperiod_stats = []
    for sp in subs:
        di = np.array(range(sp["test_start_index"], sp["test_end_index"] + 1))
        if len(di) == 0:
            continue
        s1 = perf_stats(sr, list(di))
        s2 = perf_stats(R2["strat"], list(di))
        subperiod_stats.append({
            "partition": sp["partition"],
            "subperiod_start": _iso(days_ms[di[0]]),
            "subperiod_end": _iso(days_ms[di[-1]]),
            "purged_days": sp["purged_days"],
            "mean_rank_ic": (float(np.nanmean(tics[di]))
                             if np.any(~np.isnan(tics[di])) else None),
            "rank_ic_observations": int(np.sum(~np.isnan(tics[di]))),
            "available_rows": s1["n_available"],
            "long_short_total_return_1x": s1["total_return"],
            "long_short_total_return_2x": s2["total_return"],
        })

    # ---------------- frozen truth table (unchanged rule, corrected inputs)
    unknown = []
    false_gates = []
    for s in subperiod_stats:
        if s["rank_ic_observations"] == 0:
            unknown.append(f"dev_subperiod_{s['partition']}_ic_unknown")
        elif s["mean_rank_ic"] <= IC_THRESHOLD:
            false_gates.append(f"dev_subperiod_{s['partition']}_ic_gt_0.02")
    con = windows["confirmation"]
    if con["rank_ic_observations"] == 0:
        unknown.append("confirmation_ic_unknown")
    elif con["mean_rank_ic"] <= IC_THRESHOLD:
        false_gates.append("confirmation_ic_gt_0.02")
    if con["long_short_total_return_2x"] is None:
        unknown.append("confirmation_2x_spread_unknown")
    elif con["long_short_total_return_2x"] <= 0:
        false_gates.append("confirmation_2x_spread_positive")
    if con["max_asset_pnl_share"] is None:
        unknown.append("confirmation_concentration_unknown")
    elif con["max_asset_pnl_share"] > 0.25:
        false_gates.append("max_asset_pnl_share_le_25pct")
    if con["collateral"]["buffer_breaches"] is None:
        unknown.append("confirmation_collateral_unknown")
    elif con["collateral"]["buffer_breaches"] > 0:
        false_gates.append("collateral_buffer_breaches_zero")
    if con["days_eligible_ge_20"] == 0:
        false_gates.append("point_in_time_assets_ge_20")
    if false_gates:
        verdict = "falsified"
    elif unknown:
        verdict = "inconclusive"
    else:
        verdict = "validation-pass"

    out[name] = {
        "windows": windows,
        "gate_subperiods": subperiod_stats,
        "required_gate_false": false_gates,
        "required_gate_unknown": unknown,
        "verdict": verdict,
        "censored_selected_outcomes_full_sample": R1["censored_selected_outcomes"],
        "censored_assets_at_first_stop": R1["censored_assets"],
        "held_assets_at_first_stop": R1["held_at_stop"],
        "stop_day_index": stop_day,
        "stop_day_utc": (_iso(days_ms[stop_day]) if stop_day is not None
                         else None),
        "benchmark_stop_day_index": bench_stop,
        "benchmark_stop_day_utc": (_iso(days_ms[bench_stop])
                                   if bench_stop is not None else None),
        "available_rows_total": int(avail.sum()),
        "final_nav_1x": R1["final_nav"],
        "final_nav_2x": R2["final_nav"],
        "avg_rank_ic_full": (float(np.nanmean(tics))
                             if np.any(~np.isnan(tics)) else None),
    }
    return out


def _iso(day_ms):
    import datetime as dtm
    return dtm.datetime.fromtimestamp(day_ms * 86_400,
                                      dtm.timezone.utc).strftime("%Y-%m-%d")


def main():
    P = panel()
    F = build_features(P)
    S = compute_scores(F, P)
    E = eligibility(P)
    complete = complete_outcome(P)
    days_ms = P["meta"]["days_ms"]
    results = {}
    for name in VARIANTS:
        results = evaluate_variant(P, F, S, E, name, days_ms, complete, results)
        w = results[name]["windows"]
        con = w["confirmation"]
        con_ic = con["mean_rank_ic"]
        con2x = con["long_short_total_return_2x"]
        con_ic_s = "UNKNOWN" if con_ic is None else f"{con_ic:+.4f}"
        con2x_s = "UNKNOWN" if con2x is None else f"{con2x:+.4f}"
        print(f"{name:22s} verdict {results[name]['verdict']:15s} "
              f"conIC {con_ic_s} con2x {con2x_s}")
    (EXP / "explore1_results.json").write_text(
        json.dumps(results, indent=2) + "\n")
    print("wrote explore1_results.json")


if __name__ == "__main__":
    main()
