#!/usr/bin/env python3
"""Deterministic engine for factor-momentum-timing (Ken French daily FF3 + MOM).

Wave-one discipline (campaign quant-frontier-2026-08-29):
  - development window ends 2009-12-31 (long history, earliest usable 1930-01-02)
  - validation window 2010-01-01..2024-12-31 (one frozen comparison, chosen variant only)
  - embargo: 2025-01-01 onward. run.py FAILS CLOSED if --end > 2024-12-31 unless the
    environment variable QUANT_HOLDOUT_UNLOCK=1 is set. This script must never be run
    with that variable during wave one; no selection/parameter logic reads the flag.

Units / risk-free handling (verified against raw files):
  - All Ken French daily series are PERCENT per day; converted to decimals (/100).
  - RF is the simple daily rate that compounds to the 1-month T-bill (Ibbotson
    through 2024-05, ICE BofA 1-Month T-Bill Index from 2024-06).
  - MOM is the daily 6-portfolio momentum factor (Umd), an excess-return long-short
    research spread, not a tradable instrument price.

Portfolio convention (no leverage, cash-and-overlay):
  NAV return day t = RF_t + w_{t-1} * spread_t - cost_t, with w in [0,1]
  (scale-down-only exposure to the momentum spread; collateral earns T-bills).
  Signal for day t uses data through day t-1 ONLY (strict next-bar execution).
Costs: cost_t = |w_t - w_{t-1}| * one_way_bps (contract daily ETF 5 bps; 2x stress
  = 10 bps). Benchmark = always-on momentum: RF_t + MOM_t, zero trading cost.
Secondary reference: market = RF_t + (Mkt-RF)_t.
"""
import argparse
import os
import statistics
import sys
import zipfile
from pathlib import Path

PROJ = Path(__file__).resolve().parent
DATA_DIR = PROJ.parents[1] / "data" / "raw" / "factors"
FF3_ZIP = DATA_DIR / "F-F_Research_Data_Factors_daily_CSV.zip"
MOM_ZIP = DATA_DIR / "F-F_Momentum_Factor_daily_CSV.zip"

DEV_END = 20091231
VAL_END = 20241231
DATA_CEILING = 20241231          # wave-one embargo ceiling when unlocked flag absent
HARD_DATA_MAX = 20260630         # last date present in the 2026-06 vintage files
DAYS = 252.0
COST_BPS = 5.0                   # contract: daily ETF one-way bps
LOOKBACK = 126                   # volatility management + state lookback
VOL_MEDIAN_WIN = 756             # slow state window for T05

# Mutation of the 20x-long high-volatility U.S. history is minor in money terms;
# factors here are treated as the documented Ken French research series.

# ---------------------------------------------------------------- frozen variants
# Predeclared (8 total, <= task cap of 8 and contract cap of 12):
#   T01 always-on momentum           : benchmark incumbent
#   T02 vol-managed, scale-down-only : w = min(1, 0.09 / sigma126_ann(MOM))
#   T03 vol-managed, scale-down-only : w = min(1, 0.06 / sigma126_ann(MOM))
#   T04 market trend gate            : w = 1 if cum126(Mkt-RF) > 0 else 0
#   T05 market vol gate              : w = 1 if vol126(Mkt) < median756(vol126(Mkt)) else 0
#   T06 T02 AND T04                  : vol-managed w, forced to 0 in bad market trend
#   T07 T03 AND T05                  : vol-managed w, forced to 0 in high market vol
#   T08 factor TS-momentum gate      : w = 1 if cum126(MOM) > 0 else 0
VARIANTS = [
    {"id": "T01", "name": "always_on_momentum", "kind": "always", "target": None},
    {"id": "T02", "name": "vol_managed_126d_target_0.09", "kind": "vm", "target": 0.09},
    {"id": "T03", "name": "vol_managed_126d_target_0.06", "kind": "vm", "target": 0.06},
    {"id": "T04", "name": "market_trend_gate_126d", "kind": "gate_mkt_trend", "target": None},
    {"id": "T05", "name": "market_vol_gate_126_756", "kind": "gate_mkt_vol", "target": None},
    {"id": "T06", "name": "vm_0.09_and_mkt_trend", "kind": "combo_trend", "target": 0.09},
    {"id": "T07", "name": "vm_0.06_and_mkt_vol", "kind": "combo_vol", "target": 0.06},
    {"id": "T08", "name": "factor_tsmom_gate_126d", "kind": "gate_mom_trend", "target": None},
]

# FROZEN SELECTION (set after development analysis of <=2009 data only; see
# trials.json selection_note). Must remain independent of QUANT_HOLDOUT_UNLOCK.
SELECTED = "T07"  # frozen from dev-window (<2010) analysis only


# ---------------------------------------------------------------- data loading
def _parse_pct(tok):
    v = float(tok)
    if v <= -99.0:
        raise ValueError("missing marker -99.99/-999 in parsed range")
    return v / 100.0


def load_factors(end_int):
    """Parse Ken French zips; return dicts date(int yyyymmdd) -> decimal returns."""
    dates, mom, mkt_rf, rf, smb, hml = [], {}, {}, {}, {}, {}
    with zipfile.ZipFile(FF3_ZIP) as z:
        name = [n for n in z.namelist() if n.endswith(".csv")][0]
        lines = z.read(name).decode("latin-1").splitlines()
    idx = next(i for i, ln in enumerate(lines) if ln.strip().startswith(",Mkt-RF"))
    for ln in lines[idx + 1:]:
        parts = [p for p in ln.split(",")]
        if len(parts) < 5 or not parts[0].strip().isdigit():
            break
        d = int(parts[0].strip())
        if d > end_int:
            break
        dates.append(d)
        mkt_rf[d] = _parse_pct(parts[1])
        smb[d] = _parse_pct(parts[2])
        hml[d] = _parse_pct(parts[3])
        rf[d] = _parse_pct(parts[4])
    with zipfile.ZipFile(MOM_ZIP) as z:
        name = [n for n in z.namelist() if n.endswith(".csv")][0]
        lines = z.read(name).decode("latin-1").splitlines()
    idx = next(i for i, ln in enumerate(lines) if ln.strip().startswith(",Mom"))
    for ln in lines[idx + 1:]:
        parts = ln.split(",")
        if len(parts) < 2 or not parts[0].strip().isdigit():
            break
        d = int(parts[0].strip())
        if d > end_int:
            break
        mom[d] = _parse_pct(parts[1])
    common = sorted(set(dates) & set(mom))
    if len(common) < LOOKBACK + VOL_MEDIAN_WIN + 10:
        raise SystemExit("insufficient aligned history after date cap")
    return common, mom, mkt_rf, rf, smb, hml


# ---------------------------------------------------------------- signal engine
def weights_for(variant, dates, mom, mkt):
    """w[i] is the position held on day dates[i]; uses data through i-1 only."""
    n = len(dates)
    w = [0.0] * n
    # Precompute 126d stats at each index (window [i-LOOKBACK, i) excludes day i).
    if variant["kind"] in ("vm", "combo_trend", "combo_vol"):
        sigma = [0.0] * n
        for i in range(LOOKBACK, n):
            seg = mom[i - LOOKBACK:i]
            mu = sum(seg) / LOOKBACK
            sigma[i] = (sum((x - mu) ** 2 for x in seg) / (LOOKBACK - 1)) ** 0.5 * DAYS ** 0.5
    if variant["kind"] in ("gate_mkt_vol", "combo_vol"):
        mvol = [0.0] * n
        for i in range(LOOKBACK, n):
            seg = mkt[i - LOOKBACK:i]
            mu = sum(seg) / LOOKBACK
            mvol[i] = (sum((x - mu) ** 2 for x in seg) / (LOOKBACK - 1)) ** 0.5 * DAYS ** 0.5
        med_vol = [None] * n
        for i in range(LOOKBACK + VOL_MEDIAN_WIN, n):
            med_vol[i] = statistics.median(mvol[i - VOL_MEDIAN_WIN:i])
    kind = variant["kind"]
    tgt = variant["target"]
    for i in range(n):
        if kind == "always":
            w[i] = 1.0
            continue
        wv = 1.0
        if kind in ("vm", "combo_trend", "combo_vol"):
            if i < LOOKBACK:
                continue
            wv = min(1.0, tgt / sigma[i])
        if kind == "vm":
            w[i] = wv
        elif kind == "gate_mkt_trend":
            if i < LOOKBACK:
                continue
            w[i] = 1.0 if sum(mkt[i - LOOKBACK:i]) > 0 else 0.0
        elif kind == "gate_mkt_vol":
            if med_vol[i] is None:
                continue
            w[i] = 1.0 if mvol[i] < med_vol[i] else 0.0
        elif kind == "gate_mom_trend":
            if i < LOOKBACK:
                continue
            w[i] = 1.0 if sum(mom[i - LOOKBACK:i]) > 0 else 0.0
        elif kind == "combo_trend":
            if i < LOOKBACK:
                continue
            w[i] = wv if sum(mkt[i - LOOKBACK:i]) > 0 else 0.0
        elif kind == "combo_vol":
            if med_vol[i] is None:
                continue
            w[i] = wv if mvol[i] < med_vol[i] else 0.0
    return w


def simulate(dates, mom, rf, w, cost_bps):
    """Gross/net returns on [start_idx..end_idx] plus diagnostics (net of costs)."""
    n = len(dates)
    gross, net, turnover, cost = [], [], [], []
    prev_w = 0.0
    for i in range(n):
        g = rf[i] + w[i] * mom[i]
        t = abs(w[i] - prev_w)
        c = t * cost_bps / 10000.0
        gross.append(g)
        net.append(g - c)
        turnover.append(t)
        cost.append(c)
        prev_w = w[i]
    return gross, net, turnover, cost


# ---------------------------------------------------------------- metrics
def ema(xs):
    out = 1.0
    for x in xs:
        out *= 1.0 + x
    return out


def ann_return(rs):
    n = len(rs)
    if n == 0:
        return float("nan")
    return ema(rs) ** (DAYS / n) - 1.0


def ann_vol(rs):
    n = len(rs)
    if n < 2:
        return float("nan")
    return statistics.stdev(rs) * DAYS ** 0.5


def sharpe(rs):
    n = len(rs)
    if n < 2:
        return float("nan")
    mu = sum(rs) / n
    sd = statistics.stdev(rs)
    return mu / sd * DAYS ** 0.5 if sd > 0 else float("nan")


def max_drawdown(rs):
    peak, mdd = 1.0, 0.0
    nav = 1.0
    for x in rs:
        nav *= 1.0 + x
        peak = max(peak, nav)
        mdd = min(mdd, nav / peak - 1.0)
    return mdd


def metrics(rs, turnover=None, cost=None):
    m = {
        "n_days": len(rs),
        "ann_return": round(ann_return(rs), 6),
        "ann_vol": round(ann_vol(rs), 6),
        "sharpe": round(sharpe(rs), 4),
        "max_drawdown": round(max_drawdown(rs), 6),
        "cum_return": round(ema(rs) - 1.0, 6),
        "skew": None,
        "kurt": None,
    }
    if len(rs) > 3:
        mu = sum(rs) / len(rs)
        sd = statistics.stdev(rs)
        z = [(x - mu) / sd for x in rs]
        n = len(z)
        sk = sum(v ** 3 for v in z) / n
        ku = sum(v ** 4 for v in z) / n - 3.0
        m["skew"] = round(sk, 4)
        m["kurt"] = round(ku, 4)
    if turnover is not None:
        m["ann_turnover_1w"] = round(sum(turnover) / len(turnover) * DAYS, 4)
        m["total_cost_pct"] = round(sum(cost) * 100.0, 4)
        m["ann_cost_drag_pct"] = round(sum(cost) / len(cost) * DAYS * 100.0, 4)
        m["state_flips_per_year"] = round(
            sum(1 for x in turnover if x > 0.5) / len(turnover) * DAYS, 2)
    return m

def corr(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = (sum((x - ma) ** 2 for x in a)) ** 0.5
    db = (sum((y - mb) ** 2 for y in b)) ** 0.5
    return num / (da * db) if da > 0 and db > 0 else float("nan")


# ---------------------------------------------------------------- phase runners
def run_window(dates, mom, rf, mkt, start_int, end_int, cost_bps):
    lo = next(i for i, d in enumerate(dates) if d >= start_int)
    hi = next((i for i, d in enumerate(dates) if d > end_int), len(dates))
    out = {}
    bmarks = {}
    for v in VARIANTS:
        w = weights_for(v, dates, mom, mkt)
        gross, net, turnover, cost = simulate(dates, mom, rf, w, cost_bps)
        out[v["id"]] = {
            "variant": v,
            "weights": w,
            "gross": gross[lo:hi],
            "net": net[lo:hi],
            "turnover": turnover[lo:hi],
            "cost": cost[lo:hi],
        }
    wone = [1.0] * len(dates)
    bm = [rf[i] + mom[i] for i in range(len(dates))]
    mk = [rf[i] + mkt[i] for i in range(len(dates))]
    return out, bm[lo:hi], mk[lo:hi], dates[lo:hi]


def print_dev(dates, mom, rf, mkt, start_int, end_int):
    sims, bm, mk, wdates = run_window(dates, mom, rf, mkt, start_int, end_int, COST_BPS)
    print(f"dev window {wdates[0]}..{wdates[-1]} ({len(wdates)} days) cost={COST_BPS}bps")
    print(f"{'id':4} {'sharpe_net':>10} {'sharpe_gross':>12} {'ann_ret':>8} "
          f"{'ann_vol':>8} {'maxDD':>8} {'to_1w':>7} {'flips_y':>7}")
    for tid in sorted(sims):
        s = sims[tid]
        m = metrics(s["net"], s["turnover"], s["cost"])
        mg = metrics(s["gross"])
        print(f"{tid:4} {m['sharpe']:>10} {mg['sharpe']:>12} {m['ann_return']*100:>8.2f} "
              f"{m['ann_vol']*100:>8.2f} {m['max_drawdown']*100:>8.2f} "
              f"{m['ann_turnover_1w']:>7.2f} {m['state_flips_per_year']:>7.1f}")
    mb = metrics([b for b in bm])  # always-on benchmark has zero turnover, no cost
    mmk = metrics([k for k in mk])
    print(f"benchmark always-on MOM:   sharpe={mb['sharpe']} ann_ret={mb['ann_return']*100:.2f}% "
          f"vol={mb['ann_vol']*100:.2f}% maxDD={mb['max_drawdown']*100:.2f}%")
    print(f"reference RF+Mkt-RF:       sharpe={mmk['sharpe']} ann_ret={mmk['ann_return']*100:.2f}% "
          f"vol={mmk['ann_vol']*100:.2f}% maxDD={mmk['max_drawdown']*100:.2f}%")
    # subperiod stability for every variant
    for tag, a, b in (("1930-1954", 19300101, 19541231), ("1955-1979", 19550101, 19791231),
                      ("1980-2009", 19800101, 20091231)):
        cells = []
        for tid in sorted(sims):
            s = sims[tid]
            seg = [s["net"][i] for i, d in enumerate(wdates) if a <= d <= b]
            mm = metrics(seg)
            cells.append(f"{tid}:{mm['sharpe']:.2f}/{mm['max_drawdown']*100:.1f}")
        print("  " + "  ".join(cells))
    print("\nRF annualized mean per decade (unit sanity):")
    dec = {}
    for i, d in enumerate(wdates):
        dec.setdefault(d // 10000 * 10, []).append(rf[i])
    for k in sorted(dec):
        print(f"  {k}s: {sum(dec[k])/len(dec[k])*DAYS*100:.2f}%")


def run_validation(dates, mom, rf, mkt, start_int, end_int, outdir):
    if SELECTED is None:
        raise SystemExit("validation mode requires frozen SELECTED id in run.py; "
                         "selection must come from development analysis only")
    sel = next(v for v in VARIANTS if v["id"] == SELECTED)
    sims, bm, mk, wdates = run_window(dates, mom, rf, mkt, start_int, end_int, COST_BPS)
    s = sims[SELECTED]
    mnet = metrics(s["net"], s["turnover"], s["cost"])
    mgross = metrics(s["gross"])
    mbm = metrics([b for b in bm])
    mmk = metrics(k for k in mk) if False else metrics([k for k in mk])
    nav_s = ema(s["net"])
    nav_b = ema(bm)
    excess1 = nav_s - nav_b
    # 2x cost stress
    _, net2, _, _ = (None, None, None, None)
    gross2, net2, to2, c2 = simulate(dates, mom, rf, s["weights"], COST_BPS * 2.0)
    lo = next(i for i, d in enumerate(dates) if d >= start_int)
    hi = next((i for i, d in enumerate(dates) if d > end_int), len(dates))
    net2w = net2[lo:hi]
    excess2 = ema(net2w) - nav_b
    print(f"validation {wdates[0]}..{wdates[-1]} ({len(wdates)} days) SELECTED={SELECTED}")
    print(f"selected  net : {mnet}")
    print(f"selected gross: {mgross}")
    print(f"benchmark MOM : {mbm}")
    print(f"reference Mkt : {mmk}")
    print(f"NAV final: strat(net)={nav_s:.6f} benchmark={nav_b:.6f} excess1x={excess1:.6f} "
          f"excess@2x-costs={excess2:.6f}")
    for tag, a, b in (("2010-2017", 20100101, 20171231), ("2018-2024", 20180101, 20241231)):
        segs = [s["net"][i] for i, d in enumerate(wdates) if a <= d <= b]
        segb = [bm[i] for i, d in enumerate(wdates) if a <= d <= b]
        msv, mbv = metrics(segs), metrics(segb)
        print(f"  {tag}: strat sharpe={msv['sharpe']} ann_ret={msv['ann_return']*100:.2f}% "
              f"maxDD={msv['max_drawdown']*100:.2f}% | bm sharpe={mbv['sharpe']} "
              f"ann_ret={mbv['ann_return']*100:.2f}% maxDD={mbv['max_drawdown']*100:.2f}%")
    print(f"\nnet corr vs bm={corr(s['net'], bm):+.3f} vs mkt={corr(s['net'], mk):+.3f}")
    # weight diagnostics for the selected variant
    ws = s["weights"][lo:hi]
    at_full = sum(1 for x in ws if x >= 0.999)
    mean_w = sum(ws) / len(ws)
    m2 = metrics(net2w, to2[lo:hi], c2[lo:hi])
    print(f"weight diag: mean_w={mean_w:.3f} frac_days_at_full={at_full/len(ws):.3f} "
          f"min_w={min(ws):.3f} max_w={max(ws):.3f}")
    print(f"2x-cost stress: sharpe={m2['sharpe']} cum_net={ema(net2w):.6f} "
          f"cum_bench={nav_b:.6f} excess={excess2:.6f} maxDD={m2['max_drawdown']*100:.2f}%")
    # first falsification gate verdicts
    g1 = mnet["sharpe"] > mbm["sharpe"]
    g2 = excess2 > 0.0
    g3 = mnet["max_drawdown"] > mbm["max_drawdown"] and mnet["ann_vol"] <= 0.08
    print(f"gate: sharpe_exceeds_benchmark={str(g1).lower()} "
          f"excess_positive_at_2x_costs={str(g2).lower()} "
          f"dd_below_benchmark_and_vol_within_target_tol={str(g3).lower()}")
    print(f"overall_first_falsification_gate="
          f"{'pass' if (g1 and g2 and g3) else 'fail'}")
    # returns.csv
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "returns.csv"
    with path.open("w") as f:
        f.write("date,strategy_return,benchmark_return,position,turnover,cost\n")
        for i, d in enumerate(wdates):
            f.write(f"{d//10000:04d}-{(d//100)%100:02d}-{d%100:02d},"
                    f"{s['net'][i]:.8f},{bm[i]:.8f},{s['weights'][lo+i]:.4f},"
                    f"{s['turnover'][i]:.6f},{s['cost'][i]:.8f}\n")
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="window start YYYYMMDD")
    ap.add_argument("--end", required=True, help="window end YYYYMMDD (inclusive)")
    ap.add_argument("--output-dir", default=".")
    a = ap.parse_args()
    start_int = int(a.start.replace("-", ""))
    end_int = int(a.end.replace("-", ""))
    unlock = os.environ.get("QUANT_HOLDOUT_UNLOCK") == "1"
    if end_int > VAL_END and not unlock:
        print(f"FAIL: --end {end_int} exceeds embargo ceiling {VAL_END}. "
              "Wave one must not read observations after 2024-12-31. "
              "Set QUANT_HOLDOUT_UNLOCK=1 only in the frozen evaluator.", file=sys.stderr)
        sys.exit(2)
    cap = HARD_DATA_MAX if unlock else min(VAL_END, HARD_DATA_MAX)
    dates, mom, rf, smb, hml, mkt = None, None, None, None, None, None
    common, mom, mkt_rf, rf, smb, hml = load_factors(cap)
    mkt = [mkt_rf[d] for d in common]
    rf_s = [rf[d] for d in common]
    mom_s = [mom[d] for d in common]
    # hard-reassert embargo inside the parsed data
    rd = [d for d in common if d > cap]
    if rd:
        raise SystemExit("internal embargo violation: parsed dates beyond cap")
    if end_int <= DEV_END:
        print_dev(common, mom_s, rf_s, mkt, start_int, end_int)
    elif start_int >= VAL_START_INT if False else start_int >= 20100101:
        run_validation(common, mom_s, rf_s, mkt, start_int, min(end_int, cap), a.output_dir)
    else:
        print("window spans the dev/validation boundary; pass either end<=20091231 "
              "or start>=20100101", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
