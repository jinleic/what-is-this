#!/usr/bin/env python3
"""crypto-trend-vol: deterministic backtest runner (wave one).

Strategy: slow time-series trend on daily bars aggregated from 1h Binance spot
closes (BTCUSDT, ETHUSDT). Six predeclared variants: lookback L in {20, 60, 120}
crossed with realized-vol window V in {20, 60}; long/cash equal-weight active
assets; scale-down-only volatility control to a 50% annualized target using only
lagged data; next-day execution at the daily settlement price; 12 bps one-way
costs on turnover (contract default).

Embargo policy: fails closed if the requested end date exceeds 2024-12-31,
unless QUANT_HOLDOUT_UNLOCK=1 is set in the environment. Wave one never sets
that flag; the blind evaluator may. Selection and parameters are independent of
the flag: it only widens the --end ceiling.

Usage:
  python3 run.py --start YYYY-MM-DD --end YYYY-MM-DD --output-dir DIR
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]          # .../quant-trading
SCRATCH = Path(__file__).resolve().parent / "scratch"
MANIFEST_SHA = "e1b760e96ba7d71a70188afdfa171db6ff33621d64f769fe182b37645f70a72d"
WAVE_END = pd.Timestamp("2024-12-31", tz="UTC")
COST_BPS = 12.0
TARGET_VOL = 0.50

VARIANTS = [  # (name, lookback, vol window) - the six predeclared trials
    ("trend_L20_V20", 20, 20),
    ("trend_L20_V60", 20, 60),
    ("trend_L60_V20", 60, 20),
    ("trend_L60_V60", 60, 60),
    ("trend_L120_V20", 120, 20),
    ("trend_L120_V60", 120, 60),
]
SELECTED = "trend_L120_V20"  # chosen on development net Sharpe (frozen before validation)


def fail(msg: str) -> "None":
    print(json.dumps({"error": msg}))
    sys.exit(2)


def verify_manifest() -> None:
    import hashlib
    p = REPO / "data" / "manifest.json"
    if not p.exists():
        fail(f"data manifest not found at {p}")
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != MANIFEST_SHA:
        fail(f"data manifest sha256 mismatch: {h}")


def load_daily_settlement() -> pd.DataFrame:
    """Daily settlement price per asset = price at 00:00 UTC of each date,
    aggregated locally from the collected 1h spot zips (bar start open_time,
    bar close = price one hour later). Raises on structural data problems."""
    frames = {}
    for sym in ("BTCUSDT", "ETHUSDT"):
        f = SCRATCH / f"{sym}_hourly_px.parquet"
        if not f.exists():
            fail(f"hourly cache missing: {f} (see scratch build notes in trials.json)")
        hp = pd.read_parquet(f)["close"]
        days = hp.index[hp.index.hour == 0]
        s = hp.reindex(days).dropna()
        if s.empty:
            fail(f"no settleable midnights for {sym}")
        frames[sym] = s
    px = pd.DataFrame(frames)
    if not px.index.is_monotonic_increasing:
        fail("settlement index not monotonic")
    return px


def simulate(px: pd.DataFrame, lookback: int, volwin: int, cost_bps: float):
    """pos(D) uses only Close(<=D) (both momentum endpoints and vol through the
    return that completed at 00:00 D). Earns r_next(D) = C(D+1)/C(D)-1 over day D."""
    C = px
    rets = C / C.shift(1) - 1
    mom = C / C.shift(lookback) - 1
    sig = (mom > 0).astype(float)
    n_active = sig.sum(axis=1)
    w = sig.div(n_active.where(n_active > 0), axis=0).fillna(0.0)
    vols = rets.rolling(volwin, min_periods=volwin).std() * np.sqrt(365)
    portvol = (w * vols).sum(axis=1)
    scale = (TARGET_VOL / portvol).clip(upper=1.0)          # scale-down only
    scale = scale.where(vols.notna().all(axis=1), 1.0)
    w = w.mul(scale, axis=0)
    r_next = C.shift(-1) / C - 1
    turn = (w - w.shift(1)).abs().sum(axis=1)
    cost = turn * (cost_bps / 1e4)
    sr = (w * r_next).sum(axis=1) - cost
    valid = r_next.notna().all(axis=1) & mom.notna().all(axis=1)
    bench = r_next.mean(axis=1)                             # EW buy-and-hold gross
    return pd.DataFrame({
        "strategy_return": sr.where(valid),
        "benchmark_return": bench.where(valid),
        "position": w.sum(axis=1).where(valid),
        "turnover": turn.where(valid),
        "cost": cost.where(valid),
        "scale": scale.where(valid),
    })


def metrics(r: pd.Series) -> dict:
    r = r.dropna()
    if r.empty or r.std() == 0:
        return {"n": int(len(r)), "sharpe": 0.0, "annret": 0.0, "annvol": 0.0,
                "cum": 0.0, "maxdd": 0.0}
    eq = (1 + r).cumprod()
    return {"n": int(len(r)), "sharpe": float(r.mean() / r.std() * np.sqrt(365)),
            "annret": float(r.mean() * 365), "annvol": float(r.std() * np.sqrt(365)),
            "cum": float((1 + r).prod() - 1),
            "maxdd": float((eq / eq.cummax() - 1).min())}


def main() -> None:
    ap = argparse.ArgumentParser(description="crypto-trend-vol backtest")
    ap.add_argument("--start", required=True, help="YYYY-MM-DD (inclusive)")
    ap.add_argument("--end", required=True, help="YYYY-MM-DD (inclusive)")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    verify_manifest()

    start = pd.Timestamp(args.start, tz="UTC").normalize()
    end = pd.Timestamp(args.end, tz="UTC").normalize()
    if end > WAVE_END and os.environ.get("QUANT_HOLDOUT_UNLOCK") != "1":
        fail(f"refusing --end {args.end}: exceeds wave-one ceiling 2024-12-31; "
             "set QUANT_HOLDOUT_UNLOCK=1 only if you are the frozen evaluator")

    px = load_daily_settlement()
    # Returns on date D span 00:00(D) -> 00:00(D+1); allow D+1 price strictly past --end.
    end_vir = end + pd.Timedelta(days=2)
    mask = (px.index >= start) & (px.index <= end)
    usable = px[(px.index >= start) & (px.index <= end_vir)]
    if usable.empty:
        fail("empty window after filter")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = {"start": args.start, "end": args.end, "cost_bps_one_way": COST_BPS,
              "stress_cost_bps": 2 * COST_BPS, "selected": SELECTED,
              "target_vol": TARGET_VOL, "holdout_unlocked": os.environ.get("QUANT_HOLDOUT_UNLOCK") == "1",
              "variants": {}}

    returns_frame = None
    for name, L, V in VARIANTS:
        sim = simulate(px, L, V, COST_BPS)
        sim2 = simulate(px, L, V, 2 * COST_BPS)
        sel = sim[mask].dropna(subset=["strategy_return"])
        sel2 = sim2[mask].dropna(subset=["strategy_return"])
        bh = sel["benchmark_return"].dropna()
        bh_net = bh.copy()
        if len(bh_net):
            bh_net.iloc[0] -= (COST_BPS / 1e4)  # one-time EW entry cost
        ex = sel["strategy_return"] - bh_net.reindex(sel.index).fillna(0.0)
        ex2 = sel2["strategy_return"] - bh_net.reindex(sel2.index).fillna(0.0)
        m = {
            "strategy": metrics(sel["strategy_return"]),
            "strategy_2x": metrics(sel2["strategy_return"]),
            "benchmark_ew_bh": metrics(bh_net),
            "excess_cum_1x": float((1 + ex).prod() - 1) if len(ex) else 0.0,
            "excess_cum_2x": float((1 + ex2).prod() - 1) if len(ex2) else 0.0,
            "avg_position": float(sel["position"].mean()) if len(sel) else 0.0,
            "avg_turnover": float(sel["turnover"].mean()) if len(sel) else 0.0,
        }
        report["variants"][name] = m

    # returns.csv for the selected variant (contract schema)
    sim = simulate(px, 120, 20, COST_BPS)
    rf = sim.loc[mask, ["strategy_return", "benchmark_return", "position", "turnover", "cost"]].dropna(
        subset=["strategy_return"]).copy()
    rf.index.name = "date"
    rf = rf.reset_index()
    rf["date"] = rf["date"].dt.strftime("%Y-%m-%d")
    rf["benchmark_return"] = (rf["benchmark_return"] - (COST_BPS / 1e4) * (rf.index == 0).astype(float))
    rf.to_csv(out_dir / "returns.csv", index=False)

    (out_dir / "run_report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps({k: report["variants"][SELECTED] for k in [SELECTED]}))


if __name__ == "__main__":
    main()
