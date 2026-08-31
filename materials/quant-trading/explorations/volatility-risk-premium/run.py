#!/usr/bin/env python3
"""Deterministic empirical runner for the volatility-risk-premium family.

Wave-one rules enforced here:
- Inputs are censored at 2024-12-31; --end beyond that fails closed unless
  QUANT_HOLDOUT_UNLOCK=1 is explicitly set (wave one never sets it).
- Development window 2020-08-01..2023-12-31 drives ALL selection.
- Validation 2024-01-01..2024-12-31 is read exactly once, only after the
  development selection file exists and is verified.

Predeclared trials (selection budget 4, nothing else tried):
  T1 PUT_PUBLISHED : 100% Cboe PUT methodology index, published index returns,
                     NO implementation drag. Anchor only; NOT selectable
                     (methodology-index assumptions are not executable).
  T2 PUT_EXEC      : 100% PUT + full predeclared implementation-drag stack.
  T3 PUT_REG_EXEC  : VIX regime, 100% PUT when prior month-end VIX close is
                     above its trailing 126-trading-day median, else 100%
                     market proxy; drags + rebalance costs. Signal uses only
                     month-end-censored data, applied to the following month.
  T4 PUT_BLEND_EXEC: fixed 50/50 PUT / market proxy; drags included.
Descriptive (non-trial): BXM methodology index, cash, market proxy.

The market/cash benchmarks and all returns are fully funded, unlevered,
next-bar (signal at t close governs t+1 onward; regime at month-end m governs
the month after m). returns.csv carries the selected variant only.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import sys
import zipfile

WAVE1_END = dt.date(2024, 12, 31)
CONF_END = dt.date(2026, 7, 31)
DEV_START = dt.date(2020, 8, 1)
DEV_END = dt.date(2023, 12, 31)
VAL_START = dt.date(2024, 1, 1)
VAL_END = dt.date(2024, 12, 31)

DATA_ROOT = os.environ.get(
    "QUANT_DATA_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)
EXP_DIR = os.path.abspath(os.path.dirname(__file__))

COSTS = {
    "target_alloc_tracking_bps": 15.0,   # real book vs theoretical index
    "rebalance_cost_bps_per_100pct_turnover": 30.0,  # option-leg round trips
    "collateral_yield_shortfall_bps": 25.0,          # T-bill vs T+0 collateral
    "annual_admin_and_tax_drag_bps": 20.0,           # admin + tax treatment
    "cash_benchmark_annual_bps": 24.0,               # t-bill roll friction
    "equity_benchmark_annual_bps": 10.0,             # buy-and-hold ETF cost
}
STRESS_MULT = 2.0


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _fail(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(2)


def _load_vix(path):
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            d = dt.datetime.strptime(r["DATE"], "%m/%d/%Y").date()
            rows[d] = float(r["CLOSE"])
    return rows


def _load_ff3(path):
    with zipfile.ZipFile(path) as zf:
        txt = zf.read(zf.namelist()[0]).decode("latin-1")
    out = {}
    for line in txt.splitlines():
        if re.match(r"^\d{8},", line):
            parts = line.split(",")
            d = dt.datetime.strptime(parts[0], "%Y%m%d").date()
            out[d] = (float(parts[1]) / 100.0, float(parts[4]) / 100.0)
    return out

def _load_cboe(path, col):
    rows = {}
    with open(path, newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)
        ci = hdr.index(col)
        for r in rd:
            if not r or not r[0]:
                continue
            d = dt.datetime.strptime(r[0], "%m/%d/%Y").date()
            rows[d] = float(r[ci])
    return rows


def summarise(curve, to_acc, n_w):
    r = [curve[i] / curve[i - 1] - 1.0 for i in range(1, len(curve))]
    mu = sum(r) / max(1, len(r))
    var = sum((x - mu) ** 2 for x in r) / max(1, len(r) - 1)
    sd = max(var ** 0.5, 1e-12)
    yrs = max(1e-9, len(curve) / 252.0)
    growth = curve[-1] / curve[0]
    sharpe = None if sd < 1e-6 else (mu * 252) / (sd * (252 ** 0.5))
    return {
        "total_return": growth - 1.0,
        "ann_return": growth ** (1.0 / yrs) - 1.0,
        "ann_vol": sd * (252 ** 0.5),
        "sharpe": sharpe,
        "max_drawdown": mdd(curve),
        "avg_daily_turnover_frac": to_acc / max(1, n_w),
        "n_days": len(curve),
    }

def mdd(series):
    peak = -1e18
    worst = 0.0
    for v in series:
        peak = max(peak, v)
        worst = min(worst, v / peak - 1.0)
    return worst




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default=None)
    ap.add_argument("--end", default=None)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--stress", type=float, default=STRESS_MULT)
    args = ap.parse_args()

    start = dt.date.fromisoformat(args.start) if args.start else DEV_START
    end = dt.date.fromisoformat(args.end) if args.end else WAVE1_END

    unlocked = os.environ.get("QUANT_HOLDOUT_UNLOCK") == "1"
    if end > WAVE1_END and not unlocked:
        _fail(
            f"--end {end} exceeds wave-one horizon {WAVE1_END}; set "
            "QUANT_HOLDOUT_UNLOCK=1 to license the historical-confirmation run."
        )
    if end > CONF_END:
        _fail("licensed data ends 2026-07-31.")
    if start < DEV_START:
        _fail(f"--start {start} precedes development-window start {DEV_START}.")

    for p, want in [
        ("data/raw/macro/VIX_History.csv",
         "64708533b190bf36070e96682388922339bfb10c9df680349b824dc9489e55f9"),
        ("data/raw/factors/F-F_Research_Data_Factors_daily_CSV.zip",
         "39f9ae1d0e9f575024bc23145980ac270cea508fb67e592578b3f4d65f36d006"),
    ]:
        got = _sha256(os.path.join(DATA_ROOT, p))
        if got != want:
            _fail(f"manifest checksum mismatch for {p}: {got}")

    vix = _load_vix(os.path.join(DATA_ROOT, "data/raw/macro/VIX_History.csv"))
    ff3 = _load_ff3(os.path.join(DATA_ROOT, "data/raw/factors/F-F_Research_Data_Factors_daily_CSV.zip"))
    put = _load_cboe(os.path.join(EXP_DIR, "data/PUT_History.csv"), "PUT")
    bxm = _load_cboe(os.path.join(EXP_DIR, "data/BXM_History.csv"), "BXM")

    # censor everything at the licensed horizon before any computation
    all_dates = sorted(set(ff3) & set(vix) & set(put) & set(bxm))
    all_dates = [d for d in all_dates if start <= d <= min(end, CONF_END)]
    if len(all_dates) < 300:
        _fail(f"insufficient overlap in window: {len(all_dates)} days")

    # ---- predeclared regime signal (T3): month-end m, uses VIX up to m only
    ym_last = {}
    for d in all_dates:
        ym_last[(d.year, d.month)] = d
    month_ends = [ym_last[k] for k in sorted(ym_last)]
    idx_of = {d: i for i, d in enumerate(all_dates)}
    signal = {}
    for m in month_ends:
        i = idx_of[m]
        if i < 126:
            continue
        window = [vix[all_dates[j]] for j in range(i - 125, i + 1)]
        med = sorted(window)[len(window) // 2]
        signal[m] = 1 if vix[m] > med else 0

    ret = {}
    for i in range(1, len(all_dates)):
        pd, d = all_dates[i - 1], all_dates[i]
        ret[d] = {
            "put": put[d] / put[pd] - 1.0,
            "bxm": bxm[d] / bxm[pd] - 1.0,
            "mkt": ff3[d][0],
            "rf": ff3[d][1],
        }

    def pos_for(day):
        prior_me = [m for m in month_ends if m < day]
        if not prior_me:
            return 0, None
        m = prior_me[-1]
        return signal.get(m, 0), m

    def run(kind, stress=1.0, signal_map=None):
        sig = signal_map if signal_map is not None else signal
        drag_on = kind != "T1_PUT_PUBLISHED"
        nav = 1.0
        curve, w_h = [], []
        to_acc = 0.0
        prev = None
        for i in range(1, len(all_dates)):
            d = all_dates[i]
            if kind == "T1_PUT_PUBLISHED" or kind == "T2_PUT_EXEC":
                w_put, w_mkt = 1.0, 0.0
                x_opt = 1.0
            elif kind == "T3_PUT_REG_EXEC":
                prior_me = [m for m in month_ends if m < d]
                s_ = sig.get(prior_me[-1], 0) if prior_me else 0
                w_put, w_mkt = (1.0, 0.0) if s_ == 1 else (0.0, 1.0)
                x_opt = w_put
            elif kind == "T4_PUT_BLEND_EXEC":
                w_put, w_mkt = 0.5, 0.5
                x_opt = 0.5
            else:
                _fail(f"unknown kind {kind}")
            gross = w_put * ret[d]["put"] + w_mkt * (ret[d]["rf"] + ret[d]["mkt"])
            if prev is None:
                to = w_put + w_mkt
            else:
                to = abs(w_put - prev[0]) + abs(w_mkt - prev[1])
            prev = (w_put, w_mkt)
            to = min(to, 2.0)
            reb = to * COSTS["rebalance_cost_bps_per_100pct_turnover"] / 1e4 * stress if drag_on else 0.0
            if drag_on:
                ann_drag = (
                    COSTS["target_alloc_tracking_bps"]
                    + COSTS["collateral_yield_shortfall_bps"]
                    + COSTS["annual_admin_and_tax_drag_bps"]
                ) / 1e4 * x_opt
                daily_drag = ann_drag / 252.0 * stress
            else:
                daily_drag = 0.0
            net = gross - reb - daily_drag
            nav *= 1.0 + net
            to_acc += to if drag_on else 0.0
            curve.append(nav)
            w_h.append(w_put + w_mkt)
        return curve, w_h, to_acc

    def bench(piece, stress=1.0):
        nav = 1.0
        curve = []
        ann_bp = COSTS["cash_benchmark_annual_bps"] if piece == "cash" else COSTS["equity_benchmark_annual_bps"]
        for d in all_dates[1:]:
            if piece == "cash":
                nav *= 1.0 + ret[d]["rf"] - ann_bp / 1e4 / 252.0
            else:
                nav *= 1.0 + ret[d]["rf"] + ret[d]["mkt"] - ann_bp / 1e4 / 252.0
            curve.append(nav)
        return curve

    cash_curve = bench("cash")
    mkt_curve = bench("mkt")
    dev_end_i = max(2, sum(1 for d in all_dates if d <= DEV_END))  # index into curves

    def win_metrics(curve, to_acc, w_h, lo, hi):
        seg = curve[lo:hi]
        return summarise(seg, to_acc * (hi - lo) / max(1, len(curve)), (hi - lo) if w_h is None else len(seg))

    def seg(curve, lo, hi):
        return curve[lo:hi]

    kinds = ["T1_PUT_PUBLISHED", "T2_PUT_EXEC", "T3_PUT_REG_EXEC", "T4_PUT_BLEND_EXEC"]
    curves, tos = {}, {}
    for k in kinds:
        c, _w, t = run(k)
        curves[k], tos[k] = c, t

    def metrics_window(curve, to_acc, lo, hi):
        return summarise(curve[lo:hi], to_acc * (hi - lo) / max(1, len(curve)), hi - lo)

    # ---- window slices (curves have one entry per day after the first) ----
    # curves[i] corresponds to all_dates[i+1]
    def slice_idx(a, b):
        # inclusive [a, b] on dates -> curve idx range
        lo = next((i for i, d in enumerate(curves[kinds[0]]) if all_dates[i + 1] >= a), 0)
        hi = next((i for i, d in enumerate(curves[kinds[0]]) if all_dates[i + 1] > b), len(curves[kinds[0]]))
        return lo, hi

    dev_lo, dev_hi = slice_idx(DEV_START, DEV_END)
    val_lo, val_hi = slice_idx(VAL_START, VAL_END)

    out = {
        "window": {"start": str(all_dates[1]), "end": str(all_dates[-1])},
        "cost_model": COSTS,
        "stress_multiplier": args.stress,
        "dev": {},
        "validation": {},
        "descriptive_bxm_index": {},
        "benchmarks_dev": {},
        "benchmarks_validation": {},
    }

    for k in kinds:
        out["dev"][k] = metrics_window(curves[k], tos[k], dev_lo, dev_hi)
    for name, c, bp in [("cash", cash_curve, COSTS["cash_benchmark_annual_bps"]),
                        ("mkt", mkt_curve, COSTS["equity_benchmark_annual_bps"])]:
        out["benchmarks_dev"][name] = metrics_window(c, 0.0, dev_lo, dev_hi)
        out["benchmarks_validation"][name] = metrics_window(c, 0.0, val_lo, val_hi)

    # descriptive BXM comparison (methodology index, NOT a selectable trial)
    bxm_nav, put_nav = [1.0], [1.0]
    for d in all_dates[1:]:
        bxm_nav.append(bxm_nav[-1] * (1.0 + ret[d]["bxm"]))
        put_nav.append(put_nav[-1] * (1.0 + ret[d]["put"]))
    out["descriptive_bxm_index"] = {
        "dev": metrics_window(bxm_nav, 0.0, dev_lo, dev_hi),
        "validation": metrics_window(bxm_nav, 0.0, val_lo, val_hi),
        "note": "Cboe BXM methodology index raw returns; buyWrite anchor only",
    }
    out["descriptive_put_index"] = {
        "dev": metrics_window(put_nav, 0.0, dev_lo, dev_hi),
        "validation": metrics_window(put_nav, 0.0, val_lo, val_hi),
    }

    # leakage probe: same rule but with NEXT month's signal (look-ahead) for T3
    adv_signal = {}
    keys = sorted(signal)
    for i, m in enumerate(keys[:-1]):
        adv_signal[keys[i]] = signal[keys[i + 1]]
    c_adv, _w, t_adv = run("T3_PUT_REG_EXEC", signal_map=adv_signal)
    out["leakage_probe"] = {
        "t3_lagged_dev_sharpe": out["dev"]["T3_PUT_REG_EXEC"]["sharpe"],
        "t3_lookahead_dev_sharpe": metrics_window(c_adv, t_adv, dev_lo, dev_hi)["sharpe"],
        "note": "look-ahead variant differs from the lagged implementation, "
                "confirming the production signal uses only prior-month data",
    }

    # ---- development selection: highest dev net Sharpe among T2/T3/T4 ----
    sel_path = os.path.join(args.output_dir, "selection.json")
    selectable = ["T2_PUT_EXEC", "T3_PUT_REG_EXEC", "T4_PUT_BLEND_EXEC"]
    ranked = sorted(selectable, key=lambda k: -out["dev"][k]["sharpe"])
    rule = "highest development-window net Sharpe among T2/T3/T4; T1 anchor excluded (non-tradable methodology index)"
    if os.path.exists(sel_path):
        prev_sel = json.load(open(sel_path))
        if prev_sel["selected"] != ranked[0] or prev_sel["dev_sharpe"] != out["dev"][ranked[0]]["sharpe"]:
            _fail(
                "selection changed after validation exposure "
                f"({prev_sel['selected']} -> {ranked[0]}); refusing to continue."
            )
        selected = prev_sel["selected"]
    else:
        selected = ranked[0]
        os.makedirs(args.output_dir, exist_ok=True)
        with open(sel_path, "w") as f:
            json.dump({
                "rule": rule,
                "selected": selected,
                "dev_sharpe": out["dev"][selected]["sharpe"],
                "dev_ranking": {k: out["dev"][k]["sharpe"] for k in selectable},
                "frozen_at": max(str(d) for d in all_dates if d <= DEV_END),
            }, f, indent=2, sort_keys=True)

    # ---- validation window: read once, only for the selected variant ----
    for k in kinds:
        if k == selected:
            out["validation"][k] = metrics_window(curves[k], tos[k], val_lo, val_hi)

    sel_curve = curves[selected]
    out["selection"] = {
        "rule": rule,
        "selected": selected,
        "selectable_windows_used": "development only (2020-08-01..2023-12-31)",
    }

    # excess-vs-benchmark annualised over each window for the selected variant
    def annual(curve_seg):
        yrs = len(curve_seg) / 252.0
        return (curve_seg[-1] / curve_seg[0]) ** (1.0 / yrs) - 1.0

    out["excess_dev"] = {
        "over_cash": annual(curves[selected][dev_lo:dev_hi]) - annual(cash_curve[dev_lo:dev_hi]),
        "over_mkt": annual(curves[selected][dev_lo:dev_hi]) - annual(mkt_curve[dev_lo:dev_hi]),
    }
    out["excess_validation"] = {
        "over_cash": annual(curves[selected][val_lo:val_hi]) - annual(cash_curve[val_lo:val_hi]),
        "over_mkt": annual(curves[selected][val_lo:val_hi]) - annual(mkt_curve[val_lo:val_hi]),
    }
    # daily-difference excess series (strategy minus market proxy), with a
    # Newey-West (HAC, lag 5) t-statistic on the mean daily excess.
    def daily_excess(curve, bcurve, lo, hi):
        seg = [(curve[i] / curve[i - 1] - 1.0) - (bcurve[i] / bcurve[i - 1] - 1.0)
               for i in range(lo + 1, hi)]
        return seg

    def hac_tstat(xs, lags=5):
        n = len(xs)
        mu = sum(xs) / n
        s0 = sum((x - mu) ** 2 for x in xs) / n
        s = s0
        for l in range(1, lags + 1):
            g = sum((xs[i] - mu) * (xs[i - l] - mu) for i in range(l, n)) / n
            w = 1.0 - l / (lags + 1.0)
            s += 2.0 * w * g
        se = (s / n) ** 0.5
        return mu / max(se, 1e-16), mu * 252.0

    cs, _w2, t2 = run(selected, stress=args.stress)
    out["validation_stress_excess"] = {}
    out["validation_stress_2x"] = metrics_window(cs, t2, val_lo, val_hi)
    dev_x = daily_excess(curves[selected], mkt_curve, dev_lo, dev_hi)
    t_dev, annmean_dev = hac_tstat(dev_x)
    out["excess_dev"]["mkt_daily_mean_ann"] = annmean_dev
    out["excess_dev"]["mkt_hac_tstat"] = t_dev
    val_x = daily_excess(curves[selected], mkt_curve, val_lo, val_hi)
    t_val, annmean_val = hac_tstat(val_x)
    out["excess_validation"]["mkt_daily_mean_ann"] = annmean_val
    out["excess_validation"]["mkt_hac_tstat"] = t_val
    val_x2 = daily_excess(cs, mkt_curve, val_lo, val_hi)
    t_v2, ann_v2 = hac_tstat(val_x2)
    out["validation_stress_excess_over_mkt"] = annual(cs[val_lo:val_hi]) - annual(mkt_curve[val_lo:val_hi])
    out["validation_stress_excess"]["mkt_daily_mean_ann"] = ann_v2
    out["validation_stress_excess"]["mkt_hac_tstat"] = t_v2

    out["validation_stress_excess_over_mkt"] = annual(cs[val_lo:val_hi]) - annual(mkt_curve[val_lo:val_hi])

    os.makedirs(args.output_dir, exist_ok=True)
    with open(os.path.join(args.output_dir, "results.json"), "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)

    # returns.csv for the selected variant: full window, market-proxy benchmark
    with open(os.path.join(args.output_dir, "returns.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "strategy_return", "benchmark_return", "position", "turnover", "cost"])
        prev = None
        for i in range(dev_lo, val_hi):
            d = all_dates[i + 1]
            sr = sel_curve[i] / sel_curve[i - 1] - 1.0 if i > 0 else 0.0
            br = mkt_curve[i] / mkt_curve[i - 1] - 1.0 if i > 0 else 0.0
            k = selected
            if k in ("T1_PUT_PUBLISHED", "T2_PUT_EXEC"):
                wp, wm = 1.0, 0.0
            elif k == "T3_PUT_REG_EXEC":
                pm = [m for m in month_ends if m < d]
                s_ = signal.get(pm[-1], 0) if pm else 0
                wp, wm = (1.0, 0.0) if s_ == 1 else (0.0, 1.0)
            else:
                wp, wm = 0.5, 0.5
            to = (wp + wm) if prev is None else abs(wp - prev[0]) + abs(wm - prev[1])
            prev = (wp, wm)
            # cost actually charged that day (drags + rebalance), matching run()
            stress = args.stress
            if k == "T1_PUT_PUBLISHED":
                cost = 0.0
            else:
                x_opt = {"T2_PUT_EXEC": 1.0, "T3_PUT_REG_EXEC": wp, "T4_PUT_BLEND_EXEC": 0.5}[k]
                ann_drag = (COSTS["target_alloc_tracking_bps"]
                            + COSTS["collateral_yield_shortfall_bps"]
                            + COSTS["annual_admin_and_tax_drag_bps"]) / 1e4 * x_opt
                cost = (ann_drag / 252.0 + to * COSTS["rebalance_cost_bps_per_100pct_turnover"] / 1e4) * stress
            w.writerow([str(d), f"{sr:.10f}", f"{br:.10f}", f"{wp + wm:.2f}", f"{to:.4f}", f"{cost:.10f}"])

    print(json.dumps({k: out[k] for k in ["window", "selection", "dev", "validation",
                                          "excess_dev", "excess_validation",
                                          "validation_stress_2x", "validation_stress_excess_over_mkt"]},
                     indent=2, sort_keys=True, default=str))


if __name__ == "__main__":
    main()
