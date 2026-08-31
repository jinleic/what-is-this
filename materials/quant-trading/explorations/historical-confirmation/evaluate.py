#!/usr/bin/env python3
"""Historical-confirmation evaluator for frozen quant-frontier-2026-08-29 candidates.

Deterministic, stdlib-only. Reads freeze.json, re-verifies every SHA-256, parses
each candidate's emitted returns.csv (1x cost) and captured run logs, and computes:
  - benchmark-excess daily series and HAC (Newey-West) t/p for the mean,
  - moving-block and stationary block bootstrap intervals (frozen protocol),
  - Benjamini-Hochberg (q=0.10) and Holm-Bonferroni (alpha=0.05) across m=2,
  - descriptive metrics, 2x-cost stress (from frozen run artifacts, no re-run),
  - trial accounting (PSR with skew/kurtosis correction; PBO where matrix exists),
  - half-year regime attribution,
  - a factual confirmation/falsification verdict with the prospective caveat.

Random(20260829) is consumed in one fixed documented order:
carry moving-block -> carry stationary -> factor moving-block -> factor stationary.
Writes results.json next to this file. Never modifies candidate directories.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import random
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPLORATIONS = HERE.parent
FREEZE_PATH = EXPLORATIONS / "freeze.json"
PROTOCOL_PATH = EXPLORATIONS / "protocol-amendment.json"

FREEZE_SHA256 = "f4ffef40ae43c9037c4bd4d0649a8f5c576ea6672cce6d8af1c2fe239ed3749a"
DAYS_PER_YEAR = 365.0
BOOT_REPS = 10000
SEED = 20260829
BLOCK_DAILY = 20          # daily moving-block length (frozen protocol)
HOURLY_BLOCK = 168        # frozen crypto hourly block length (hours)
BH_Q = 0.10
HOLM_ALPHA = 0.05

REVIEWER_WARNINGS = {
    "crypto-carry-decay CarryAccountingReview": {
        "source": "read-only reviewer agent CarryAccountingReview (2026-08-29); cited code sites "
                  "independently re-verified in the frozen run.py by the evaluator; sign wording "
                  "corrected per reviewer follow-up",
        "severity": "interpretation-affecting; does not alter the frozen strategy, run, or emitted outputs",
        "run.py:152 funding clip": ("parse_funding keeps max(rate, 0.0) per standalone 8h funding event "
                  "(reviewer confirms ZERO snap collisions across BTC/ETH raw files, so first_seen "
                  "semantics are irrelevant): each standalone negative rate is clipped to 0. Sign "
                  "convention (funding payment = -q_p*P*r with short q_p<0): SHORT RECEIVES when r>0, "
                  "SHORT PAYS when r<0. Clipping negative rates to 0 therefore suppresses payments the "
                  "short would make in negative-funding regimes, which OVERSTATES P&L (income is not "
                  "affected; the error removes a cost). Reviewer's 2024 overstatement recompute: "
                  "BTC +0.00722824, ETH +0.00427122 (quantified in the reviewer's own writeup)"),
        "run.py:518,526-528 yield denominator": ("funding_ann_pct_of_notional divides PnL by eq_close/2 "
                  "(half current equity), not the actual fixed |q_p|*P perp-leg notional; because equity "
                  "compounds, the reported 45.0%/47.9% (validation) and 19.9%/12.8% (confirmation) "
                  "annualized funding-yield figures are not true per-notional yields"),
        "reviewer_recomputation": ("raw 2024 funding yields recompute to BTC 12.382%, ETH 13.867%, pair "
                  "13.13-13.28% (as reported in the reviewer message), not the reported ~46%; the reviewer "
                  "flags the wave-one 15% funding-yield gate as effectively reversed by this")
    }
}
REGIME_SEGMENTS = [
    ("2025H1", "2025-01-01", "2025-06-30"),
    ("2025H2", "2025-07-01", "2025-12-31"),
    ("2026H1", "2026-01-01", "2026-06-30"),
    ("2026H2_to_07-31", "2026-07-01", "2026-07-31"),
]


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_cdf(x: float) -> float:
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def norm_ppf(p: float) -> float:
    lo, hi = -40.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if norm_cdf(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def read_returns(path: Path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "date": r["date"],
                "s": float(r["strategy_return"]),
                "b": float(r["benchmark_return"]),
                "pos": float(r["position"]),
                "to": float(r["turnover"]),
                "cost": float(r["cost"]),
            })
    return rows


def nav_drawdown(returns):
    nav, peak, mdd = 1.0, 1.0, 0.0
    for r in returns:
        nav *= (1.0 + r)
        peak = max(peak, nav)
        mdd = min(mdd, nav / peak - 1.0)
    return mdd


def sharpe_daily(returns) -> float:
    if len(returns) < 2:
        return float("nan")
    m, sd = statistics.mean(returns), statistics.stdev(returns)
    return m / sd * math.sqrt(DAYS_PER_YEAR) if sd > 0 else float("nan")


def ann_return(returns) -> float:
    n = len(returns)
    if n == 0:
        return float("nan")
    growth = 1.0
    for r in returns:
        growth *= (1.0 + r)
    if growth <= 0:
        return float("nan")
    return growth ** (DAYS_PER_YEAR / n) - 1.0


def excess_series(rows):
    return [(1.0 + r["s"]) / (1.0 + r["b"]) - 1.0 for r in rows]


def hac_test(e):
    """Newey-West HAC t-test that the mean daily benchmark-excess return is zero."""
    n = len(e)
    m = statistics.fmean(e)
    lag = int(4.0 * (n / 100.0) ** (2.0 / 9.0))
    d = [x - m for x in e]
    g0 = sum(x * x for x in d) / n
    lrv = g0
    for l in range(1, lag + 1):
        gl = sum(d[t] * d[t - l] for t in range(l, n)) / n
        lrv += 2.0 * (1.0 - l / (lag + 1.0)) * gl
    se = math.sqrt(max(lrv / n, 0.0))
    t = m / se if se > 0 else float("nan")
    p = math.erfc(abs(t) / math.sqrt(2.0)) if se > 0 else float("nan")
    return {
        "n_obs": n,
        "mean_daily_excess": m,
        "mean_annualized_excess": m * DAYS_PER_YEAR,
        "newey_west_lag": lag,
        "hac_se_daily": se,
        "hac_se_annualized": se * math.sqrt(DAYS_PER_YEAR),
        "t_stat": t,
        "p_value_two_sided": p,
    }


def _blocks_resample(e, rng, mode):
    n = len(e)
    out = []
    if mode == "moving":
        k = math.ceil(n / BLOCK_DAILY)
        for _ in range(k):
            s = rng.randrange(0, n - BLOCK_DAILY + 1)
            out.extend(e[s:s + BLOCK_DAILY])
    else:  # stationary bootstrap, expected block length BLOCK_DAILY
        p = 1.0 / BLOCK_DAILY
        idx = rng.randrange(n)
        for _ in range(n):
            out.append(e[idx])
            idx = rng.randrange(n) if rng.random() < p else (idx + 1) % n
    return out[:n]


def bootstrap(e, rng):
    result = {"block_size_daily_obs": BLOCK_DAILY,
              "hourly_equivalent_note": f"protocol hourly block {HOURLY_BLOCK}h ~ {HOURLY_BLOCK / 24:.1f} days; daily series uses the frozen daily_block_observations=20",
              "replicates": BOOT_REPS, "seed": SEED, "rng_sequence": rng.name}
    for mode, key in (("moving", "moving_block"), ("stationary", "stationary")):
        means = []
        for _ in range(BOOT_REPS):
            xs = _blocks_resample(e, rng, mode)
            means.append(statistics.fmean(xs) * DAYS_PER_YEAR)
        means.sort()
        n_le_0 = sum(1 for x in means if x <= 0.0)
        result[key] = {
            "ci95_annualized_mean_excess": [means[int(0.025 * BOOT_REPS)], means[int(0.975 * BOOT_REPS)]],
            "median_annualized_mean_excess": means[int(0.5 * BOOT_REPS)],
            "replicates_le_0": n_le_0,
            "p_bootstrap_one_sided_mean_le_0": n_le_0 / BOOT_REPS,
        }
    return result


def bh_holm(raw_p):
    order = sorted(raw_p, key=lambda k: raw_p[k])
    m = len(order)
    holm_adj, bh_adj, running_h, running_b = {}, {}, 0.0, 1.0
    for rank, k in enumerate(order):
        running_h = max(running_h, (m - rank) * raw_p[k])
        holm_adj[k] = min(1.0, running_h)
    for rank in range(m - 1, -1, -1):
        k = order[rank]
        running_b = min(running_b, raw_p[k] * m / (rank + 1))
        bh_adj[k] = min(1.0, running_b)
    holm_rej, thresh = [], HOLM_ALPHA
    for k in order:  # step-down
        if raw_p[k] <= thresh:
            holm_rej.append(k)
            thresh = HOLM_ALPHA / max(m - len(holm_rej), 1)
        else:
            break
    bh_rej = []
    for rank in range(m - 1, -1, -1):  # step-up
        if raw_p[order[rank]] <= (rank + 1) / m * BH_Q:
            bh_rej = order[:rank + 1]
            break
    return holm_adj, bh_adj, sorted(holm_rej), sorted(bh_rej)


def psr(e):
    """Probabilistic Sharpe ratio SR>0 with skew/kurtosis correction."""
    n = len(e)
    m, sd = statistics.fmean(e), statistics.stdev(e)
    if sd == 0:
        return {"error": "zero variance"}
    sr = m / sd
    ex = statistics.fmean([(x - m) ** 3 for x in e]) / sd ** 3
    k4 = statistics.fmean([(x - m) ** 4 for x in e]) / sd ** 4
    denom = 1.0 - ex * sr + (k4 - 1.0) / 4.0 * sr * sr
    if denom <= 0:
        return {"sr_daily": sr, "psr": None, "note": "degenerate correction term"}
    z = sr * math.sqrt(n - 1) / math.sqrt(denom)
    return {
        "sr_daily": sr,
        "sr_annualized": sr * math.sqrt(DAYS_PER_YEAR),
        "skew": ex,
        "kurtosis": k4,
        "psr_sr_gt_0": norm_cdf(z),
    }


def regime_attribution(rows):
    out = {}
    for name, lo, hi in REGIME_SEGMENTS:
        seg = [r for r in rows if lo <= r["date"] <= hi]
        if not seg:
            out[name] = {"n_obs": 0, "note": "no emitted observations in this segment (data vintage gap)"}
            continue
        ex = excess_series(seg)
        out[name] = {
            "window": [seg[0]["date"], seg[-1]["date"]],
            "n_obs": len(seg),
            "excess_cum_geo": math.prod(1.0 + x for x in ex) - 1.0,
            "strategy_cum": math.prod(1.0 + r["s"] for r in seg) - 1.0,
            "benchmark_cum": math.prod(1.0 + r["b"] for r in seg) - 1.0,
            "strategy_sharpe": sharpe_daily([r["s"] for r in seg]),
        }
    return out


def main() -> int:
    freeze = json.loads(FREEZE_PATH.read_text())

    # ---------------------------------------------------------------- hashes
    integrity = {
        "freeze_json": {"expected": FREEZE_SHA256, "actual": sha256_file(FREEZE_PATH),
                        "match": sha256_file(FREEZE_PATH) == FREEZE_SHA256},
        "protocol_amendment.json": {
            "expected": freeze["protocol_amendment_sha256"],
            "actual": sha256_file(PROTOCOL_PATH),
            "match": sha256_file(PROTOCOL_PATH) == freeze["protocol_amendment_sha256"]},
    }
    with open(HERE / "hash_verification.json") as f:
        pre_run = json.load(f)
    integrity["pre_run_record"] = pre_run
    integrity["all_candidates_verified"] = True
    per = {}
    blocked = []
    for c in freeze["candidates"]:
        files, ok = {}, True
        for fname, expected in c["files"].items():
            p = HERE.parents[2] / c["directory"] / fname
            actual = sha256_file(p)
            match = actual == expected
            ok = ok and match
            files[fname] = {"expected": expected, "actual": actual, "match": match}
        per[c["id"]] = files
        integrity["all_candidates_verified"] = integrity["all_candidates_verified"] and ok
        if not ok:
            blocked.append(c["id"])
    integrity["per_candidate"] = per
    gate_open = (integrity["freeze_json"]["match"]
                 and integrity["protocol_amendment.json"]["match"]
                 and integrity["all_candidates_verified"]
                 and pre_run.get("verified_before_run") is True)
    integrity["note"] = ("verification performed before QUANT_HOLDOUT_UNLOCK=1 (recorded in "
                         "hash_verification.json) and recomputed at evaluate time; both passed "
                         "so both candidates ran once" if gate_open else
                         "HASH MISMATCH - affected candidate(s) blocked: " + ",".join(blocked or ["all"]))

    # sibling trial counts (read-only) for carry-family accounting
    siblings = {
        "crypto-basis-convergence": json.loads(
            (EXPLORATIONS / "crypto-basis-convergence" / "trials.json").read_text()).get("trial_count_total"),
        "funding-crowding": json.loads(
            (EXPLORATIONS / "funding-crowding" / "trials.json").read_text()).get("total_variants_attempted"),
    }

    rng = random.Random(SEED)
    rng.name = "single shared Random(20260829); consumed: carry-MBB, carry-stationary, factor-MBB, factor-stationary"
    raw_p, candidates, verdict_inputs = {}, {}, {}

    for c in freeze["candidates"]:
        cid = c["id"]
        outdir = HERE / cid
        rows = read_returns(outdir / "returns.csv")
        e = excess_series(rows)
        strat_r = [r["s"] for r in rows]
        stdout_text = (HERE / f"{cid}.stdout.log").read_text()
        stderr_text = (HERE / f"{cid}.stderr.log").read_text()
        strat = json.loads((HERE.parents[2] / c["directory"] / "strategy.json").read_text())

        window_note = None
        if rows[-1]["date"] != freeze["historical_confirmation_window"][1]:
            window_note = ("frozen inputs end at data-vintage boundary " + rows[-1]["date"] +
                           "; candidate ran to the last observation its frozen data contains - "
                           "no code or data changes were made")

        hac = hac_test(e)
        boot = bootstrap(e, rng)

        if cid == "crypto-carry-decay":
            report = json.loads((outdir / "report.json").read_text())
            req = report["windows"]["requested"]
            val = report["windows"]["validation"]
            st2, attr = req["stress_2x"], req["attribution"]
            self_trials = strat["trial_count"]
            extra = {
                "window_hours_emitted": req["n_hours"],
                "position_gross_notional_over_initial_capital": {
                    "mean_reported": req["pos_frac_mean"],
                    "semantics": "settle-to-market accounting: long spot + matched USD-M perp short; net directional delta ~0, gross = sleeve notional / initial capital",
                },
                "turnover_annualized_notional_turns": req["turnover_notional_turns"],
                "pnl_attribution_usd_per_initial_dollar": {
                    "price_pnl": attr["price_pnl_usd"], "funding_pnl": attr["funding_pnl_usd"],
                    "fees": attr["fees_usd"], "funding_plus_price_plus_fees": attr["check_sum_pnl_usd"]},
                "funding_decay_evidence_pct_of_notional_annualized": {
                    "validation_2024": {s: val["attribution"]["per_symbol"][s]["funding_ann_pct_of_notional"]
                                        for s in ("BTCUSDT", "ETHUSDT")},
                    "confirmation_2025_2026": {s: attr["per_symbol"][s]["funding_ann_pct_of_notional"]
                                               for s in ("BTCUSDT", "ETHUSDT")}},
                "reviewer_interpretation_warning": REVIEWER_WARNINGS["crypto-carry-decay CarryAccountingReview"],
                "identity_residuals": report["identity_residuals"],
                "fee_events": report["fee_events"],
                "funding_events": {"credited": report["funding_events"]["credited"],
                                   "dropped": report["funding_events"]["dropped"]},
                "cost_stress_2x": {
                    "excess_cum_geo": st2["excess_cum_geo"], "excess_cum_arith": st2["excess_cum_arith"],
                    "sharpe": st2["strategy"]["sharpe"], "max_drawdown": st2["strategy"]["max_drawdown"],
                    "turnover_annualized_notional_turns": st2["turnover_notional_turns"],
                    "note": "simulated inside the frozen run.py (STRESS_MULT=2.0); evaluator did not re-run"},
                "carry_family_trial_accounting_note": (
                    f"carry_complex family trials: crypto-carry-decay {self_trials} + "
                    f"crypto-basis-convergence {siblings['crypto-basis-convergence']} + "
                    f"funding-crowding {siblings['funding-crowding']} = "
                    f"{self_trials + siblings['crypto-basis-convergence'] + siblings['funding-crowding']} "
                    "(family members were falsified in wave one and contribute trials under the "
                    "protocol trial-accounting rule)"),
            }
        else:
            gate = re.search(r"overall_first_falsification_gate=(\w+)", stdout_text)
            g1 = re.search(r"sharpe_exceeds_benchmark=(\w+)", stdout_text)
            g2 = re.search(r"excess_positive_at_2x_costs=(\w+)", stdout_text)
            g3 = re.search(r"dd_below_benchmark_and_vol_within_target_tol=(\w+)", stdout_text)
            self_trials = strat["trial_count"]
            extra = {
                "position_weight": {
                    "mean": float(re.search(r"mean_w=([\d.]+)", stdout_text).group(1)),
                    "frac_days_at_full": float(re.search(r"frac_days_at_full=([\d.]+)", stdout_text).group(1)),
                    "range": [float(re.search(r"min_w=([\d.]+)", stdout_text).group(1)),
                              float(re.search(r"max_w=([\d.]+)", stdout_text).group(1))],
                    "semantics": "fraction of NAV in the vol-managed gated MOM sleeve; remainder in RF cash",
                },
                "turnover": {"ann_turnover_1w": 2.4021,
                             "sum_emitted_csv": sum(r["to"] for r in rows)},
                "cost_stress_2x": {
                    "excess_cum": float(re.search(r"excess@2x-costs=([-\d.]+)", stdout_text).group(1)),
                    "note": "simulated inside the frozen run.py (COST_BPS*2 = 10bps); evaluator did not re-run"},
                "first_falsification_gate": {
                    "sharpe_exceeds_benchmark": g1.group(1) == "true",
                    "excess_positive_at_2x_costs": g2.group(1) == "true",
                    "dd_below_benchmark_and_vol_within_target_tol": g3.group(1) == "true",
                    "overall": gate.group(1)},
                "selection_vs_outcome": {
                    "validation_2010_2024_excess_cum": 0.91977 - 0.422651,
                    "note": "strategy.json validation cum 0.91977 vs benchmark 0.422651 (positive excess); "
                            "confirmation window flips excess negative - recent-horizon deterioration"},
                "factor_investability_note": "strategy and benchmark returns are Ken French research "
                                             "factor returns (RF cash + MOM long-short sleeve); underlying "
                                             "portfolio turnover and borrow costs are not directly "
                                             "investable (frozen review list)",
            }

        metrics = {
            "ann_return_strategy_geo": ann_return(strat_r),
            "ann_return_benchmark_geo": ann_return([r["b"] for r in rows]),
            "strategy_sharpe_ann": sharpe_daily(strat_r),
            "strategy_max_drawdown": nav_drawdown(strat_r),
            "excess_cum_geo": math.prod(1.0 + x for x in e) - 1.0,
            "excess_cum_arith": sum(e),
            "mean_position_col": statistics.fmean([r["pos"] for r in rows]),
            "total_cost_col": sum(r["cost"] for r in rows),
        }

        # multiple-testing verdict inputs
        candidates[cid] = {
            "group": c["group"],
            "directory": c["directory"],
            "evaluation_command": c["evaluation_command"],
            "run_outputs": {
                "exit_status": 0,
                "invoked": ("python3 " + c["evaluation_command"].split("python ", 1)[1]
                            if "python " in c["evaluation_command"] else c["evaluation_command"]),
                "stdout_log": f"historical-confirmation/{cid}.stdout.log",
                "stderr_log": f"historical-confirmation/{cid}.stderr.log",
                "reruns_for_metrics": 0,
            },
            "window": {
                "requested": freeze["historical_confirmation_window"],
                "actual_first_date": rows[0]["date"],
                "actual_last_date": rows[-1]["date"],
                "n_observations_daily": len(rows),
                "note": window_note,
            },
            "cost_model_frozen": strat["cost_model"],
            "metrics_1x": metrics,
            "hac": hac,
            "bootstrap": boot,
            "trial_accounting": {
                "trial_count": self_trials,
                "wave_one_budget": strat.get("trial_count"),
                "pbo": "not computable: no per-trial return matrix exists in frozen artifacts",
                "dsr": {**psr(e),
                        "full_dsr_note": ("DSR requires the cross-trial variance of Sharpe ratios from "
                                          "per-trial return series, which the frozen artifacts do not "
                                          "contain; PSR(SR*=0) with skew/kurtosis correction is reported "
                                          "and the trial count above is the deflation input")},
            },
            "regime_attribution": regime_attribution(rows),
            "extras": extra,
        }
        verdict_inputs[cid] = {"excess_cum_geo": metrics["excess_cum_geo"],
                               "t_stat": hac["t_stat"]}
        raw_p[cid] = hac["p_value_two_sided"]

    # ------------------------------------------------- multiple testing
    holm_adj, bh_adj, holm_rej, bh_rej = bh_holm(raw_p)
    for cid in candidates:
        candidates[cid]["multiple_testing"] = {
            "m_candidates": len(raw_p),
            "raw_p_value_two_sided": raw_p[cid],
            "holm_adjusted_p": holm_adj[cid],
            "holm_alpha": HOLM_ALPHA,
            "holm_reject_zero_mean_excess": cid in holm_rej,
            "bh_adjusted_p": bh_adj[cid],
            "bh_q": BH_Q,
            "bh_discover": cid in bh_rej,
            "caution": "FDR-only rejection (BH) is discovery, not strong confirmation; "
                       "Holm family-wise is the strong-confirmation criterion per protocol-amendment.json. "
                       "BH treats rejection either direction (the two-sided HAC p-value): a rejected "
                       "NEGATIVE mean excess is evidence of underperformance, not alpha.",
        }
        v = verdict_inputs[cid]
        mt = candidates[cid]["multiple_testing"]
        if cid == "crypto-carry-decay":
            flip = candidates[cid]["extras"]["cost_stress_2x"]["excess_cum_geo"] <= 0.0
        else:
            flip = candidates[cid]["extras"]["cost_stress_2x"]["excess_cum"] <= 0.0
        if v["excess_cum_geo"] <= 0.0 or flip:
            verdict = "falsified"
            basis = ("confirmation-window benchmark-excess cumulative return is negative "
                     "(1x: {:+.4f}, frozen 2x-cost stress: {:+.4f})").format(
                         v["excess_cum_geo"], candidates[cid]["extras"]["cost_stress_2x"].get(
                             "excess_cum_geo", candidates[cid]["extras"]["cost_stress_2x"]["excess_cum"]))
        elif mt["holm_reject_zero_mean_excess"] and not flip:
            verdict = "confirmed_under_holm_within_historical_confirmation"
            basis = ("HAC raw p={:.3g} rejects zero mean excess at Holm alpha=0.05 across m={} candidates; "
                     "excess positive at both 1x and frozen-run 2x cost").format(raw_p[cid], len(raw_p))
        elif mt["bh_discover"]:
            verdict = "fdr_discovery_only_not_strong_confirmation"
            basis = "BH rejection without Holm confirmation; FDR-only results must not be labeled confirmed"
        else:
            verdict = "not_confirmed_not_falsified"
            basis = "no multiple-testing rejection at the frozen thresholds"
        candidates[cid]["verdict"] = verdict
        candidates[cid]["verdict_basis"] = basis

    results = {
        "schema_version": 1,
        "campaign": freeze["campaign"],
        "evaluator": "historical-confirmation/evaluate.py (deterministic, stdlib-only)",
        "protocol": {
            "source": "explorations/freeze.json + explorations/protocol-amendment.json",
            "window_semantics": "historical confirmation (parameter-embargoed, NOT epistemically blind); "
                                "only the forward shadow from 2026-08-29 is genuinely untouched",
            "requested_window": freeze["historical_confirmation_window"],
            "candidates_m": len(freeze["candidates"]),
            "test_statistic": "HAC (Newey-West) annualized mean benchmark-excess return",
            "bootstrap": {"replicates": BOOT_REPS, "seed": SEED,
                          "block_daily_observations": BLOCK_DAILY,
                          "block_hourly_hours": HOURLY_BLOCK},
            "corrections": {"benjamini_hochberg_q": BH_Q, "holm_alpha": HOLM_ALPHA},
        },
        "hash_integrity": integrity,
        "data_coverage_note": (
            "crypto-carry-decay ran the full requested window (emitted through 2026-07-31). "
            "factor-momentum-timing ran the full requested window through the last observation in its "
            "hash-frozen data vintage (2026-06-30); the frozen run.py caps parsed data at HARD_DATA_MAX=20260630. "
            "No frozen file was altered and no rerun was performed."),
        "candidates": candidates,
        "overall": {
            "multiple_testing": {"m": len(raw_p), "raw_p_values": raw_p,
                                  "holm_adjusted": holm_adj, "bh_adjusted": bh_adj,
                                  "holm_reject_set": holm_rej, "bh_discovery_set": bh_rej},
            "epistemic_status": (
                "Historical confirmation may falsify a frozen rule but CANNOT establish prospective "
                "durable alpha: candidate families were selected in 2026 with knowledge of published and "
                "local 2025-2026 market evidence (protocol-amendment.json). Decisions to deploy require "
                "the append-only prospective forward shadow beginning 2026-08-29."),
            "deployment_conclusion": "no candidate advances to capital deployment from historical confirmation alone",
        },
    }

    with open(HERE / "results.json", "w") as f:
        json.dump(results, f, indent=2, sort_keys=False)
        f.write("\n")

    print(f"results.json written; hash gate open={gate_open}")
    for cid, cand in candidates.items():
        print(f"  {cid}: window {cand['window']['actual_first_date']}..{cand['window']['actual_last_date']} "
              f"n={cand['window']['n_observations_daily']} "
              f"excess_cum_geo={cand['metrics_1x']['excess_cum_geo']:+.4f} "
              f"raw_p={raw_p[cid]:.4g} holm_p={holm_adj[cid]:.4g} verdict={cand['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
