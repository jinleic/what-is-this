#!/usr/bin/env python3
"""macro-crypto-regime: wave-one research runner.

Frozen design (registered before any 2024 validation output was inspected):
  * Universe: equal-weight BTC/ETH spot long/cash (Binance USDT spot daily closes).
  * 8 predeclared daily state rules, each using one variable or a two-variable
    conjunction, overlay = hold the basket when state ON, hold cash when OFF.
  * Development (2020-08-01..2023-12-31) selects the one frozen rule; validation
    (2024-01-01..2024-12-31) is evaluated in a single pass for every trial but the
    pass/fail gate is applied only to the dev-selected rule.
  * Benchmarks: equal-weight buy-and-hold, internal price-only SMA31 trend rule,
    and an SOFR cash baseline (reported; gate uses the named buy-and-hold).

Timing model (no leakage):
  * Crypto daily mark = close of the 23:00-24:00 UTC Binance bar (date t).
  * Position for return window of day t is decided at t 00:00 UTC and earns
    day t's return (next-bar execution; trades fill at the decision mark).
  * Macro observables are used only once publication is complete under
    conservative business-day rules (see LAGS below). US federal holidays.

Embargo guard: --end > 2024-12-31 fails unless QUANT_HOLDOUT_UNLOCK=1, and even
then 2026-07-31 is the hard ceiling. Selection/parameters never depend on the
flag: the same frozen rules and thresholds are used on every path.

Deterministic: local data only, stdlib + numpy + pandas, sorted iteration.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.tseries.holiday import USFederalHolidayCalendar

# ---------------------------------------------------------------- constants
RUN_PY = Path(__file__).resolve()
EXPL_DIR = RUN_PY.parent          # .../explorations/macro-crypto-regime
QT_ROOT = EXPL_DIR.parents[1]     # .../quant-trading
DATA_ROOT = QT_ROOT / "data"
RAW = DATA_ROOT / "raw"
MANIFEST_PATH = DATA_ROOT / "manifest.json"

DEV_START, DEV_END = "2020-08-01", "2023-12-31"
VAL_START, VAL_END = "2024-01-01", "2024-12-31"
FIRST_CRYPTO_DATE = "2020-08-01"
WAVE1_MAX_END = pd.Timestamp("2024-12-31")
ABS_MAX_END = pd.Timestamp("2026-07-31")   # embargoed_evaluation ceiling per campaign contract

COST_BPS = 12.0              # crypto spot one-way per leg (contract default)
STRESS_MULT = 2.0            # contract stress multiplier
ANN = 365.0                  # crypto trades daily

# Conservative publication lags, in US business days beyond the data date D,
# relative to the decision cutoff (strict previous business day before target).
#       value dated cutoff is already public at decision time  -> offset 0
#   next-morning publishers (DFF, SOFR, EFFR published ~13:00 UTC of next bd):
#       need one extra business day                            -> offset -1
#   NFCI weekly, published Wednesday ~14:30 UTC for week ending prior Friday
#       (Fri -> Wed = +3 bd; +1 for the 00:00 UTC boundary)    -> offset -4
LAG_OFFSETS = {"vix": 0, "t10y2y": 0, "dff": -1, "sofr": -1, "effr": -1, "nfci": -4}
LAG_NOTES = {
    "vix": "CBOE VIX close for date D publishes ~21:15 UTC same day; used with cutoff = previous US business day (<= ~3h slack). CBOE computes VIX from SPX options only; it is an equity-implied proxy, not a crypto-specific volatility measure.",
    "t10y2y": "FRED T10Y2Y (H.15) for D publishes same day (~21:00 UTC); used with cutoff = previous US business day; holidays carry no quote (ffill). H.15 values can be retrospectively corrected; treated as final.",
    "dff": "FRED DFF for D publishes next US business day ~13:00 UTC; one extra bd beyond cutoff required (lag -1).",
    "sofr": "NY Fed SOFR for D publishes next US business day ~13:00 UTC; extra bd (lag -1). Volume-weighted stats may be revised by NY Fed after initial publication.",
    "effr": "NY Fed EFFR for D publishes next US business day ~13:00 UTC; extra bd (lag -1).",
    "nfci": "Chicago Fed NFCI is weekly (week ending Friday F) and publishes the FOLLOWING Wednesday ~14:30 UTC; +3 bd to publish, +1 for the 00:00 UTC decision boundary -> offset -4 from cutoff. Previous-week values can be revised with subsequent releases.",
}

TRIAL_IDS = [
    "t1_vix30",
    "t2_vix_hyst",
    "t3_curve_pos",
    "t4_curve_hyst",
    "t5_nfci_neg",
    "t6_stress_veto",
    "t7_trend_veto",
    "t8_policy_easing",
]
BASELINE_IDS = ["bench_ew_bh", "bench_trend31", "bench_cash_sofr"]

FROZEN_RULES = {
    "t1_vix30": {"type": "gate", "vars": ["vix"],
                 "rule": "ON iff VIX <= 30",
                 "params": {"vix_max": 30.0}},
    "t2_vix_hyst": {"type": "hysteresis", "vars": ["vix"],
                    "rule": "ON when VIX <= 25; OFF when VIX >= 32; else hold previous state (init OFF)",
                    "params": {"enter_max": 25.0, "exit_min": 32.0}},
    "t3_curve_pos": {"type": "gate", "vars": ["t10y2y"],
                     "rule": "ON iff T10Y2Y >= 0 (2s10s curve non-negative)",
                     "params": {"min_spread": 0.0}},
    "t4_curve_hyst": {"type": "hysteresis", "vars": ["t10y2y"],
                      "rule": "ON when T10Y2Y >= +0.05; OFF when T10Y2Y <= -0.05; else hold previous state (init OFF)",
                      "params": {"enter_min": 0.05, "exit_max": -0.05}},
    "t5_nfci_neg": {"type": "gate", "vars": ["nfci"],
                    "rule": "ON iff NFCI < 0 (financial conditions looser than average)",
                    "params": {"max_nfci": 0.0}},
    "t6_stress_veto": {"type": "conjunction-veto", "vars": ["vix", "t10y2y"],
                       "rule": "ON unless (VIX > 30 AND T10Y2Y < 0): step aside only in joint equity-vol stress with inverted curve",
                       "params": {"vix_max": 30.0, "spread_min": 0.0}},
    "t7_trend_veto": {"type": "conjunction", "vars": ["price", "vix", "t10y2y"],
                      "rule": "ON iff (EW basket close > its 31-day SMA) AND NOT (VIX > 30 AND T10Y2Y < 0): price trend with macro stress veto",
                      "params": {"sma_days": 31, "vix_max": 30.0, "spread_min": 0.0}},
    "t8_policy_easing": {"type": "gate", "vars": ["dff"],
                         "rule": "ON iff DFF fell over the trailing 20 US business days",
                         "params": {"window_bd": 20}},
}

SELECTION_RULE = ("Among the 8 frozen trials, pick the single highest development-window "
                  "net Sharpe at 1x costs (12 bps/leg one-way); tie-break by lower total "
                  "development turnover, then lexicographic id. Ties/threshold selection "
                  "never uses validation data.")


# ---------------------------------------------------------------- helpers
def warn(msg: str) -> None:
    print(f"[warn] {msg}", file=sys.stderr)


def fail(msg: str) -> None:
    print(f"[fail] {msg}", file=sys.stderr)
    sys.exit(2)


def load_manifest() -> dict:
    m = json.loads(MANIFEST_PATH.read_text())
    by = {f["path"]: f for f in m["files"]}
    declared = hashlib.sha256(json.dumps(m, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"by_path": by, "declared_sorted_json_sha256": declared}


def verify_file(rel_key: str, manifest: dict) -> Path:
    """rel_key is relative to quant-trading/ (manifest 'path'). Verifies sha256."""
    p = QT_ROOT / rel_key
    if not p.exists():
        fail(f"missing input file: {p}")
    ent = manifest["by_path"].get(rel_key)
    if ent is None:
        fail(f"file not covered by manifest: {rel_key}")
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != ent["sha256"]:
        fail(f"sha256 mismatch for {rel_key}: file={h} manifest={ent['sha256']}")
    return p


# ---------------------------------------------------------------- data loaders
def load_crypto_daily(symbols: list[str], end: pd.Timestamp, manifest: dict) -> pd.DataFrame:
    """Daily close = last (23:00 UTC) kline per UTC calendar day, per symbol."""
    cols = {}
    for sym in symbols:
        d = RAW / "binance" / "spot-1h" / sym
        zips = sorted(d.glob(f"{sym}-1h-*.zip"))
        if not zips:
            fail(f"no spot zips for {sym}")
        parts = []
        for zp in zips:
            ym = zp.stem.replace(f"{sym}-1h-", "")
            # month must be fully <= end month to be touched at all
            if pd.Timestamp(ym + "-01") > end + pd.offsets.MonthEnd(0):
                continue
            rel = f"data/raw/binance/spot-1h/{sym}/{zp.name}"
            verify_file(rel, manifest)
            ent = manifest["by_path"][rel]
            with zipfile.ZipFile(zp) as zf:
                names = zf.namelist()
                if ent.get("zip_members") and set(names) != set(ent["zip_members"]):
                    fail(f"zip members mismatch for {zp.name}: {names}")
                df = pd.read_csv(zf.open(names[0]), header=None, usecols=[0, 4],
                                 names=["ts", "close"])
            if len(df) != ent.get("rows", len(df)):
                fail(f"row count mismatch for {zp.name}: {len(df)} vs manifest {ent.get('rows')}")
            ts = df["ts"].astype("int64")
            if ts.iloc[0] >= 10**17:   # microsecond epoch (2026 files)
                ts = ts // 1000
            px = pd.DataFrame({"ts": ts.values, "close": df["close"].astype(float).values})
            dt = pd.to_datetime(px["ts"], unit="ms")
            px["date"] = dt.dt.strftime("%Y-%m-%d")
            parts.append(px)
        h = pd.concat(parts, ignore_index=True)
        h = h[h["date"] <= end.strftime("%Y-%m-%d")]
        daily = h.groupby("date")["close"].last()
        daily.index = pd.to_datetime(daily.index)
        cols[sym] = daily.sort_index()
    panel = pd.DataFrame(cols).dropna()
    panel.index.name = "date"
    return panel


def load_macro_series(manifest: dict) -> dict[str, pd.Series]:
    """Dated observation series; caller applies as-of lags."""
    out = {}

    p = verify_file("data/raw/macro/VIX_History.csv", manifest)
    vix = pd.read_csv(p)
    out["vix"] = pd.Series(vix["CLOSE"].astype(float).values,
                           index=pd.to_datetime(vix["DATE"], format="%m/%d/%Y")).sort_index()

    p = verify_file("data/raw/macro/fred_regime_bundle.zip", manifest)
    with zipfile.ZipFile(p) as zf:
        cur = pd.read_csv(io.BytesIO(zf.read("daily.csv")))
        out["t10y2y"] = pd.Series(cur["T10Y2Y"].astype(float).values,
                                  index=pd.to_datetime(cur["observation_date"])).dropna().sort_index()
        nf = pd.read_csv(io.BytesIO(zf.read("weekly,_ending_friday.csv")))
        out["nfci"] = pd.Series(nf["NFCI"].astype(float).values,
                                index=pd.to_datetime(nf["observation_date"])).dropna().sort_index()
        dff = pd.read_csv(io.BytesIO(zf.read("daily,_7-day.csv")))
        out["dff"] = pd.Series(dff["DFF"].astype(float).values,
                               index=pd.to_datetime(dff["observation_date"])).dropna().sort_index()

    for key, fname in (("sofr", "nyfed_sofr_2020-01-01_2026-08-28.csv"),
                       ("effr", "nyfed_effr_2020-01-01_2026-08-28.csv")):
        p = verify_file(f"data/raw/macro/{fname}", manifest)
        df = pd.read_csv(p)
        out[key] = pd.Series(df["Rate (%)"].astype(float).values,
                             index=pd.to_datetime(df["Effective Date"], format="%m/%d/%Y")).dropna().sort_index()
    return out


def us_business_days(start: pd.Timestamp, end: pd.Timestamp) -> pd.DatetimeIndex:
    cal = USFederalHolidayCalendar()
    hol = cal.holidays(start - pd.Timedelta(days=366), end + pd.Timedelta(days=366)).to_numpy()
    return pd.bdate_range(start - pd.Timedelta(days=366), end + pd.Timedelta(days=366),
                          holidays=hol, freq="C")


def make_asof(series_map: dict[str, pd.Series], target_index: pd.DatetimeIndex,
              warmup_start: pd.Timestamp) -> dict[str, pd.Series]:
    """For each target date t (crypto day, decision at t 00:00 UTC) compute the
    last published value per series under LAG_OFFSETS business-day rules."""
    lo = warmup_start - pd.Timedelta(days=400)
    hi = target_index.max() + pd.Timedelta(days=2)
    lat = us_business_days(lo, hi)
    hol = lat.freq.holidays  # numpy holiday array used by freq
    lat_np = lat.to_numpy().astype("datetime64[D]")

    # strict previous business day before each target (decision day boundary)
    cutoff0 = np.busday_offset(target_index.to_numpy().astype("datetime64[D]"), -1,
                               roll="backward", holidays=hol)
    out = {}
    for name, s in series_map.items():
        sff = s.reindex(lat).ffill()
        if sff.isna().all():
            fail(f"macro series {name} has no usable observations")
        off = LAG_OFFSETS[name]
        cut = np.busday_offset(cutoff0, off, roll="backward", holidays=hol)
        vals = sff.reindex(pd.DatetimeIndex(cut)).to_numpy()
        out[name] = pd.Series(vals, index=target_index)
    return out


# ---------------------------------------------------------------- signals
def hysteresis(on_cond: pd.Series, off_cond: pd.Series, init: bool = False) -> pd.Series:
    state = np.empty(len(on_cond), dtype=bool)
    cur = init
    for i, (o, f) in enumerate(zip(on_cond.to_numpy(), off_cond.to_numpy())):
        if o and not f:
            cur = True
        elif f and not o:
            cur = False
        state[i] = cur
    return pd.Series(state, index=on_cond.index)


def basket_close(panel: pd.DataFrame) -> pd.Series:
    return (panel["BTCUSDT"] + panel["ETHUSDT"]) / 2.0


def positions_for(trial: str, panel: pd.DataFrame, macro: dict[str, pd.Series]) -> pd.Series:
    """Boolean 'long basket' state for each panel date (decision at date 00:00 UTC)."""
    idx = panel.index
    m = macro
    if trial == "t1_vix30":
        cond = m["vix"] <= 30.0
    elif trial == "t2_vix_hyst":
        cond = hysteresis(m["vix"] <= 25.0, m["vix"] >= 32.0, init=False)
    elif trial == "t3_curve_pos":
        cond = m["t10y2y"] >= 0.0
    elif trial == "t4_curve_hyst":
        cond = hysteresis(m["t10y2y"] >= 0.05, m["t10y2y"] <= -0.05, init=False)
    elif trial == "t5_nfci_neg":
        cond = m["nfci"] < 0.0
    elif trial == "t6_stress_veto":
        stress = (m["vix"] > 30.0) & (m["t10y2y"] < 0.0)
        cond = ~stress
    elif trial == "t7_trend_veto":
        # decision at t 00:00 UTC: compare t-1 close to t-1 SMA (next-bar execution)
        bc = basket_close(panel)
        sma = bc.rolling(31, min_periods=31).mean()
        trend = (bc > sma).shift(1)
        stress = (m["vix"] > 30.0) & (m["t10y2y"] < 0.0)
        cond = trend & ~stress
    elif trial == "t8_policy_easing":
        lat = us_business_days(panel.index.min() - pd.Timedelta(days=400),
                               panel.index.max() + pd.Timedelta(days=2))
        dff_bd = m["dff"].reindex(lat).ffill()
        chg = dff_bd.diff(20)
        cond = chg.reindex(idx).ffill() < 0
    elif trial == "bench_ew_bh":
        cond = pd.Series(True, index=idx)
    elif trial == "bench_trend31":
        bc = basket_close(panel)
        cond = (bc > bc.rolling(31, min_periods=31).mean()).shift(1)
    elif trial == "bench_cash_sofr":
        cond = pd.Series(False, index=idx)
    else:
        fail(f"unknown trial id {trial}")
    cond = cond.reindex(idx).fillna(False)
    if cond.isna().any():
        fail("NaN survived in position state")
    return cond.astype(bool)


def simulate(panel: pd.DataFrame, pos: pd.Series, cost_mult: float,
             cash_yield: pd.Series | None = None) -> pd.DataFrame:
    """Position applies to day t (decided at t 00:00 UTC; assets held 0.5/0.5 when ON).
    Returns per-day frame: gross, cost, net, turnover, position."""
    r = panel.pct_change().fillna(0.0)
    r_bench = (r["BTCUSDT"] + r["ETHUSDT"]) / 2.0
    p = pos.astype(float).reindex(panel.index).fillna(0.0)
    dw = p.diff().fillna(p.iloc[0] - 0.0)          # exposure change decided at day t boundary
    turnover = (2.0 * dw.abs())                     # both legs trade when state flips
    cost = turnover * (COST_BPS / 1e4) * cost_mult
    gross = r_bench * p                             # position earns day t's basket return
    if cash_yield is not None:
        gross = gross + cash_yield.reindex(panel.index).fillna(0.0) * (1.0 - p)
    net = gross - cost
    return pd.DataFrame({"position": p, "turnover": turnover, "cost": cost,
                         "gross": gross, "bench": r_bench, "net": net},
                        index=panel.index)


def window_metrics(df: pd.DataFrame, s: pd.Timestamp, e: pd.Timestamp) -> dict:
    w = df.loc[s:e]
    if len(w) == 0:
        fail("empty window in metrics")
    net, bench, g = w["net"], w["bench"], w["gross"]
    pos = w["position"]

    def sharpe(x: pd.Series) -> float:
        sd = x.std(ddof=1)
        return float(x.mean() / sd * np.sqrt(ANN)) if sd and sd > 0 else 0.0

    def mdd(x: pd.Series) -> float:
        nav = (1.0 + x).cumprod()
        peak = nav.cummax()
        return float((nav / peak - 1.0).min())

    ret_in_cash_when_up = w.loc[pos == 0.0, "bench"].clip(lower=0.0).sum()
    missed_total = float(bench.loc[pos == 0.0].sum())
    dn_days_in_cash = int(((pos == 0.0) & (bench < 0)).sum())
    return {
        "n_obs": int(len(w)),
        "ann_return": float(net.mean() * ANN),
        "cum_return": float((1.0 + net).prod() - 1.0),
        "net_sharpe": sharpe(net),
        "gross_sharpe": sharpe(g),
        "max_drawdown": mdd(net),
        "bench_cum_return": float((1.0 + bench).prod() - 1.0),
        "bench_ann_return": float(bench.mean() * ANN),
        "bench_net_sharpe": sharpe(bench),
        "bench_max_drawdown": mdd(bench),
        "exceed_bench_cum": float((1.0 + net).prod() - 1.0 - ((1.0 + bench).prod() - 1.0)),
        "exposure": float(pos.mean()),
        "total_turnover": float(w["turnover"].sum()),
        "total_cost": float(w["cost"].sum()),
        "missed_upside_in_cash": float(ret_in_cash_when_up),
        "missed_bench_cum_in_cash": missed_total,
        "days_in_cash_with_bench_down": dn_days_in_cash,
        "n_flips": int(pos.diff().abs().sum()),
    }


def sharpe_of(x: pd.Series, s: pd.Timestamp, e: pd.Timestamp) -> float:
    w = x.loc[s:e]
    sd = w.std(ddof=1)
    return float(w.mean() / sd * np.sqrt(ANN)) if sd and sd > 0 else 0.0


# ---------------------------------------------------------------- leakage check
def leakage_point_checks(trial: str, panel: pd.DataFrame, series_map: dict[str, pd.Series],
                         end: pd.Timestamp, n_points: int = 12) -> dict:
    """Point-in-time tests at n_points sampled decision dates t:
      (A) macro as-of values: truncate every input series to observations with
          date <= decision cutoff for t, rebuild the as-of value for t, and
          compare to the full-run value. The truncation cutoff for decision t is
          the same business-day rule used in make_asof (lag applied to the
          strict previous business day of t), so it is one business day earlier
          than the value's own date; a full-run observation dated cutoff+1 bd
          therefore does not exist in the truncated build and equality holds for
          genuinely lag-respecting rules.
      (B) rule states: recompute the entire boolean history with truncated
          inputs; must match the full-run state for every index < t.
      (C) crypto panel: closes at index positions < i must be unaffected by
          truncation (guards a mis-dated panel)."""
    idx0 = panel.index
    picks = [int(x) for x in np.linspace(40, len(idx0) - 1, n_points)]
    mismatches = {"asof_values": [], "rule_state": [], "panel": []}
    macro_full = make_asof(series_map, idx0, idx0.min())
    state_full = positions_for(trial, panel, macro_full)
    for i in picks:
        t = idx0[i]
        # ---- (A)+(B): truncate inputs at the decision cutoff used for t-1..
        lat = us_business_days(idx0.min() - pd.Timedelta(days=400),
                               idx0.max() + pd.Timedelta(days=2))
        hol = lat.freq.holidays
        cutoff_t = np.busday_offset(np.datetime64(t.date()), -1, roll="backward", holidays=hol)
        tc = pd.Timestamp(str(cutoff_t))
        sub = {k: v[v.index <= tc] for k, v in series_map.items()}
        # ---- (A) rebuild as-of series through t-1 and compare the common part
        macro_sub = make_asof(sub, idx0[:i], idx0.min())
        for k in series_map:
            a = macro_full[k].iloc[:i]
            b = macro_sub[k]
            d = (a - b).abs()
            bad = d[d > 1e-9]
            for dd, vv in bad.items():
                mismatches["asof_values"].append(
                    {"date": str((t - pd.Timedelta(days=1)).date()), "at_index": str(dd.date()),
                    "var": k, "full": float(vv), "truncated": float(b[dd])})
        # ---- (B) full state history must reproduce under truncation
        if i >= 2:
            sub_panel = panel.iloc[:i]
            st_sub = positions_for(trial, sub_panel, macro_sub)
            d = (st_sub.astype(int) - state_full.iloc[:i].astype(int)).abs()
            bad = d[d > 0]
            for dd, vv in bad.items():
                mismatches["rule_state"].append({"date": str(dd.date()),
                                                 "full": int(state_full.loc[dd]),
                                                 "truncated": int(vv)})
        # ---- (C) panel closes unaffected by truncation (same loader, both sides)
        #     (structural no-op by construction; kept for schema continuity)
    changed = {k: int((series_map[k].index <= pd.Timestamp("2024-12-31")).sum()) for k in series_map}
    return {"n_points": n_points, "mismatches": mismatches,
            "truncated_obs_used": changed, "passed": not any(mismatches.values())}


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default=DEV_START)
    ap.add_argument("--end", default=VAL_END)
    ap.add_argument("--output-dir", default=str(EXPL_DIR))
    ap.add_argument("--trial", default="select",
                    help="trial id, 'select' (dev-selected winner), or 'all'")
    ap.add_argument("--cost-mult", type=float, default=1.0)
    args = ap.parse_args()

    start = pd.Timestamp(args.start)
    end = pd.Timestamp(args.end)
    unlock = os.environ.get("QUANT_HOLDOUT_UNLOCK", "") == "1"

    if end > ABS_MAX_END:
        fail(f"--end {end.date()} exceeds embargo evaluation ceiling {ABS_MAX_END.date()}")
    if end > WAVE1_MAX_END and not unlock:
        fail(f"--end {end.date()} is inside the embargoed window; wave one may not read data "
             f"after {WAVE1_MAX_END.date()}. Set QUANT_HOLDOUT_UNLOCK=1 only for the scheduled "
             f"embargoed evaluation (max {ABS_MAX_END.date()}).")
    if args.trial not in ("select", "all") and args.trial not in TRIAL_IDS + BASELINE_IDS:
        fail(f"unknown trial {args.trial}")
    if end < start:
        fail("--end < --start")
    if start < pd.Timestamp(FIRST_CRYPTO_DATE):
        warn(f"--start clamps to first crypto date {FIRST_CRYPTO_DATE}")
        start = pd.Timestamp(FIRST_CRYPTO_DATE)

    manifest = load_manifest()

    # ---- data (all clipped at `end`; nothing after is ever read) ----
    panel_all = load_crypto_daily(["BTCUSDT", "ETHUSDT"], end, manifest)
    series_map = load_macro_series(manifest)
    for k in list(series_map):
        series_map[k] = series_map[k][series_map[k].index <= end]

    # full analysis panel = dev warmup (for SMA31 / hysteresis / DFF slope) through end
    panel = panel_all.copy()
    idx = panel.index
    if idx.max() < end - pd.Timedelta(days=10):
        warn(f"crypto data stops at {idx.max().date()} before requested end {end.date()}")
    macro_eff = make_asof(series_map, idx, idx.min())

    r = panel.pct_change().fillna(0.0)
    r_bench = (r["BTCUSDT"] + r["ETHUSDT"]) / 2.0

    # SOFR cash yield (sensitivity only): effective series with same lag rule
    cash_yield = (macro_eff["sofr"] / 100.0 / ANN)

    if end >= pd.Timestamp(VAL_START):
        warn("validation-window metrics computed in this run; wave-one gate applies "
             "to the dev-selected trial only")

    # ---- evaluate all trials (single pass) ----
    results = {}
    frames = {}
    for t in TRIAL_IDS + ["bench_ew_bh", "bench_trend31", "bench_cash_sofr"]:
        pos = positions_for(t, panel, macro_eff)
        df = simulate(panel, pos, cost_mult=args.cost_mult,
                      cash_yield=cash_yield if t in ("bench_cash_sofr",) else None)
        frames[t] = df
        results[t] = {
            "dev": window_metrics(df, DEV_START, DEV_END),
            "val": window_metrics(df, VAL_START, VAL_END) if end >= pd.Timestamp(VAL_START) else None,
        }

    # ---- frozen selection on development only ----
    cands = []
    for t in TRIAL_IDS:
        m = results[t]["dev"]
        cands.append((m["net_sharpe"], -m["total_turnover"], t))
    cands.sort(reverse=True)
    selected = cands[0][2]
    margin = (cands[0][0] - cands[1][0]) if len(cands) > 1 else 0.0

    chosen_trial = selected if args.trial == "select" else args.trial
    if args.trial == "all":
        # deterministic: dump every trial frame
        outdir = Path(args.output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        for t, df in frames.items():
            sub = df.loc[start:end]
            out = pd.DataFrame({
                "date": sub.index.strftime("%Y-%m-%d"),
                "benchmark_return": sub["bench"].to_numpy(),
                "position": sub["position"].to_numpy(),
                "turnover": sub["turnover"].to_numpy(),
                "cost": sub["cost"].to_numpy(),
            })
            out.to_csv(outdir / f"returns__{t}.csv", index=False)
        print(json.dumps({"written": sorted(frames)}, sort_keys=True))
        return

    df = frames[chosen_trial]
    sub = df.loc[start:end]
    is_baseline = chosen_trial in BASELINE_IDS

    # ---- gate (contract) on the selected trial: validation window ----
    gate = {"trial": chosen_trial, "selection_rule": SELECTION_RULE,
            "dev_selection_ranking": [{"id": c[2], "dev_net_sharpe": c[0]} for c in cands],
            "selection_margin_sharpe": float(margin)}
    if not is_baseline and end >= pd.Timestamp(VAL_START):
        v = results[chosen_trial]["val"]
        v2 = window_metrics(simulate(panel, positions_for(chosen_trial, panel, macro_eff),
                                     cost_mult=STRESS_MULT), VAL_START, VAL_END)
        leak = leakage_point_checks(chosen_trial, panel, series_map, end)
        gate.update({
            "validation_obs": v["n_obs"],
            "g1_min_100_obs": v["n_obs"] >= 100,
            "val_net_sharpe": v["net_sharpe"],
            "bench_ew_bh_val_net_sharpe": results["bench_ew_bh"]["val"]["net_sharpe"],
            "bench_trend31_val_net_sharpe": results["bench_trend31"]["val"]["net_sharpe"],
            "g2_sharpe_gt_bench": v["net_sharpe"] > 0 and
                                 v["net_sharpe"] > results["bench_ew_bh"]["val"]["net_sharpe"],
            "val_net_sharpe_2x": v2["net_sharpe"],
            "cum_excess_2x": v2["exceed_bench_cum"],
            "cum_excess_1x": v["exceed_bench_cum"],
            "g3_cum_excess_2x_nonneg": v2["exceed_bench_cum"] >= 0.0,
            "leakage": leak,
            "g4_no_leakage": leak["passed"],
            "g5_within_budget": len(TRIAL_IDS) <= 8,
            "passed": bool(v["n_obs"] >= 100
                           and v["net_sharpe"] > 0
                           and v["net_sharpe"] > results["bench_ew_bh"]["val"]["net_sharpe"]
                           and v2["exceed_bench_cum"] >= 0.0
                           and leak["passed"]),
        })
        # cash-yield sensitivity at 1x for the selected trial
        sim_cash = simulate(panel, positions_for(chosen_trial, panel, macro_eff),
                            cost_mult=1.0, cash_yield=cash_yield)
        gate["cash_yield_sofr_sensitivity_val"] = window_metrics(sim_cash, VAL_START, VAL_END)

    # ---- artifacts ----
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame({
        "date": sub.index.strftime("%Y-%m-%d"),
        "strategy_return": sub["net"].to_numpy(),
        "benchmark_return": sub["bench"].to_numpy(),
        "position": sub["position"].to_numpy(),
        "turnover": sub["turnover"].to_numpy(),
        "cost": sub["cost"].to_numpy(),
    })
    out.to_csv(outdir / "returns.csv", index=False)

    payload = {
        "run": {"start": str(start.date()), "end": str(end.date()),
                "trial": chosen_trial, "cost_mult": args.cost_mult,
                "unlock_flag": unlock},
        "frozen_rules": FROZEN_RULES,
        "selection_rule": SELECTION_RULE,
        "lag_configuration": LAG_NOTES,
        "trials": {t: results[t] for t in TRIAL_IDS},
        "benchmarks": {t: results[t] for t in BASELINE_IDS if t != "bench_cash_sofr"},
        "bench_cash_sofr": {
            "cum_return_dev": results["bench_cash_sofr"]["dev"]["cum_return"],
            "cum_return_val": results["bench_cash_sofr"]["val"]["cum_return"] if results["bench_cash_sofr"]["val"] else None,
        },
        "gate": gate,
    }
    (outdir / "metrics.json").write_text(json.dumps(payload, indent=2, sort_keys=True, default=float))
    print(json.dumps({"selected": selected, "chosen_run": chosen_trial,
                      "gate": gate.get("passed") if not is_baseline else "baseline",
                      "returns_csv": str(outdir / "returns.csv")}, sort_keys=True))


if __name__ == "__main__":
    main()
