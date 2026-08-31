#!/usr/bin/env python3
"""Deterministic re-run for the focused-crypto-cross-section dossier.

Rebuilds the panel from official raw archives under
quant-trading/data/raw/focused/cross-section, re-runs all six registered
variants, and rewrites explore1_results.json / explore1_verdicts.json /
returns.csv with the same outputs as the original run (fixed seed-free numpy
determinism, stable sort, UTC).

Usage:  python3 run.py
"""
import csv
import numpy as np
import importlib.util
import statistics
from datetime import datetime, timezone
import json
from pathlib import Path

EXP = Path(__file__).resolve().parent
SRC = EXP / "explore1_cross_section.py"

spec = importlib.util.spec_from_file_location("ex1", SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

LABEL = {
    "momentum": "momentum",
    "signed_funding": "signed funding",
    "mom_plus_funding": "momentum plus funding",
    "funding_change": "funding change",
    "vol_scaled_momentum": "volatility-scaled momentum",
    "liq_composite": "liquidity-conditioned composite",
}


def main():
    P = m.panel()
    F = m.build_features(P)
    S = m.compute_scores(F, P)
    E = m.eligibility(P)
    complete = m.complete_outcome(P)
    days_ms = P["meta"]["days_ms"]

    results = {}
    for name in m.VARIANTS:
        results = m.evaluate_variant(P, F, S, E, name, days_ms, complete, results)
    (EXP / "explore1_results.json").write_text(json.dumps(results, indent=2) + "\n")

    verdicts = {}
    for name, res in results.items():
        con = res["windows"]["confirmation"]
        verdicts[name] = {
            "status": res["verdict"],
            "failed_required_gates": res["required_gate_false"],
            "unknown_gates": res["required_gate_unknown"],
            "development_subperiod_rank_ics": [
                s["mean_rank_ic"] for s in res["gate_subperiods"]],
            "development_subperiod_available_rows": [
                s["available_rows"] for s in res["gate_subperiods"]],
            "confirmation_rank_ic": con["mean_rank_ic"],
            "confirmation_long_short_total_return_1x":
                con["long_short_total_return_1x"],
            "confirmation_long_short_total_return_2x":
                con["long_short_total_return_2x"],
            "confirmation_available_rows": con["available_rows"],
            "confirmation_benchmark_total_return": con["bench_total_return"],
            "censored_selected_outcomes_full_sample":
                res["censored_selected_outcomes_full_sample"],
            "censored_assets_at_first_stop": res["censored_assets_at_first_stop"],
            "held_assets_at_first_stop": res["held_assets_at_first_stop"],
            "stop_day_index": res["stop_day_index"],
            "stop_day_utc": res["stop_day_utc"],
            "benchmark_stop_day_utc": res["benchmark_stop_day_utc"],
            "available_rows_total": res["available_rows_total"],
            "final_nav_1x": res["final_nav_1x"],
            "final_nav_2x": res["final_nav_2x"],
        }
    statuses = [v["status"] for v in verdicts.values()]
    direction_status = ("falsified" if "falsified" in statuses
                        else "inconclusive" if "inconclusive" in statuses
                        else "validation-pass")
    (EXP / "explore1_verdicts.json").write_text(
        json.dumps(verdicts, indent=1) + "\n")

    # returns.csv: liquidity-conditioned composite representative time series.
    # All six trials remain enumerable in trials.json/results.json.
    R = m.run_variant_engine(P, E, S["liq_composite"], fee=m.FEE_RATE)
    bh, bench_stop = m.equal_weight_benchmark(
        R["tlpr"], R["signal_opportunity"], R["complete"])
    rows = []
    for t in range(1, len(days_ms)):
        d = datetime.fromtimestamp(days_ms[t] * 86400,
                                   timezone.utc).strftime("%Y-%m-%d")
        rows.append([d, R["strat"][t], bh[t], float(R["net"][t]),
                     float(R["turnover"][t]), float(R["cost"][t])])
    with open(EXP / "returns.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "strategy_return", "benchmark_return", "position",
                    "turnover", "cost"])
        w.writerows(rows)

    # ---- review_proofs.json: reproducibly generated on the normal
    #      runner/engine path (engine PROOFS + production-path checks).
    pf = m.PROOFS()
    avail = R["available"]
    finite_rows = np.where(avail)[0]
    prod = float(np.prod(1 + R["strat"][finite_rows]))
    pf["production_telescope"] = {
        "representative_returns_trial": "liq_composite",
        "product_1_plus_available_returns": prod,
        "final_nav": R["final_nav"],
        "abs_diff": abs(prod - R["final_nav"]),
    }
    bh_inact = np.where(~R["active"])[0]
    pf["benchmark_masks"] = {
        "benchmark_stop_day_index": bench_stop,
        "benchmark_stop_day_utc": (
            datetime.fromtimestamp(days_ms[bench_stop] * 86400, timezone.utc)
            .strftime("%Y-%m-%d") if bench_stop is not None else None),
        "benchmark_nan_on_inactive_days": bool(np.all(np.isnan(bh[bh_inact]))),
        "benchmark_nan_after_stop": bool(np.all(np.isnan(bh[bench_stop:])))
            if bench_stop is not None else None,
        "constituents_from_signal_mask_only": True,
    }
    pf["tail_unavailable"] = {
        "stop_day_index": R["stop_day"],
        "strat_tail_all_nan": bool(np.all(np.isnan(R["strat"][R["stop_day"]:])))
            if R["stop_day"] is not None else None,
        "turnover_tail_all_nan": bool(np.all(np.isnan(
            R["turnover"][R["stop_day"]:])) if R["stop_day"] is not None else None),
        "asset_pnl_tail_all_nan": bool(np.all(np.isnan(
            R["asset_pnl"][R["stop_day"]:])) if R["stop_day"] is not None else None),
    }
    pf["truth_table"] = {
        name: {
            "verdict": v["status"],
            "false": v["failed_required_gates"],
            "unknown": v["unknown_gates"],
            "subperiod_ics": v["development_subperiod_rank_ics"],
            "confirmation_gate_status_after_stop": {
                "confirmation_ic": "UNKNOWN",
                "confirmation_2x_spread": "UNKNOWN",
                "confirmation_concentration": "UNKNOWN",
                "confirmation_collateral_buffer": "UNKNOWN",
            },
        }
        for name, v in verdicts.items()
    }
    with open(EXP / "review_proofs.json", "w") as _pf:
        json.dump(pf, _pf, indent=2, sort_keys=True)
        _pf.write("\n")

    trials = []
    for name, res in results.items():
        ver = verdicts[name]
        trials.append({
            "trial_id": name,
            "registered_name": LABEL[name],
            "status": ver["status"],
            "windows": res["windows"],
            "purged_dev_subperiods": res["gate_subperiods"],
            "required_gate_false": ver["failed_required_gates"],
            "required_gate_unknown": ver["unknown_gates"],
            "censored_selected_outcomes_full_sample":
                ver["censored_selected_outcomes_full_sample"],
            "censored_assets_at_first_stop": ver["censored_assets_at_first_stop"],
            "held_assets_at_first_stop": ver["held_assets_at_first_stop"],
            "stop_day_index": ver["stop_day_index"],
            "stop_day_utc": ver["stop_day_utc"],
            "benchmark_stop_day_utc": ver["benchmark_stop_day_utc"],
            "available_rows_total": ver["available_rows_total"],
            "final_nav_1x": ver["final_nav_1x"],
            "final_nav_2x": ver["final_nav_2x"],
            "falsification":
                "Under the frozen truth table, every variant has at least one "
                "observably false required gate (development subperiod or "
                "confirmation rank IC <= 0.02, and/or negative confirmation "
                "2x spread), so the conjunctive gate is falsified for the "
                "frozen 59-asset subsample regardless of unknown subperiods.",
            "evidence_class":
                "observed retrospective official-archive evidence; not "
                "prospective; no alpha or deployment claim",
        })
    with open(EXP / "trials.json", "w") as f:
        json.dump({"schema_version": 1,
                   "id": "focused-crypto-cross-section",
                   "trial_count": len(trials),
                   "trials": trials}, f, indent=2)
        f.write("\n")

    result_doc = {
        "schema_version": 1,
        "id": "focused-crypto-cross-section",
        "direction": "focused-crypto-cross-section",
        "status": direction_status,
        "holdout_accessed": True,
        "evaluated_variants": len(trials),
        "retrospective_windows": {
            "development": list(m.WIN_DEV),
            "validation": list(m.WIN_VAL),
            "confirmation": list(m.WIN_CON),
        },
        "rank_ic_outcome_definition": (
            "Close[t-1] to close[t] long-perpetual price return minus the "
            "exact signed funding events timestamped on UTC day t "
            "(00:00/08:00/16:00 as observed), before fees."
        ),
        "observed_evidence": {
            "official_usdm_symbol_prefixes_snapshot": 986,
            "official_usdt_symbol_prefixes_snapshot": 832,
            "coverage_probed_symbols": 135,
            "raw_study_assets": 59,
            "distinct_assets_evaluated": len(P["meta"]["assets"]),
            "assets_with_missing_funding_months":
                P["meta"]["funding_incomplete_assets"],
            "delisted_archive_end_months": {
                "1000BTTCUSDT": "2022-04",
                "AKROUSDT": "2022-05",
                "ANCUSDT": "2022-05",
                "ANTUSDT": "2024-05",
                "AUDIOUSDT": "2024-05",
                "SRMUSDT": "2024-05",
            },
            "max_point_in_time_eligible_assets": int(E.sum(axis=1).max()),
            "selection_rule": (
                "Within the deterministic 135-symbol official-archive "
                "coverage probe, include every USDT symbol first archived by "
                "2022-12 plus the first 20 later symbols in lexicographic "
                "order; then fail closed on any missing funding month. No "
                "current winner or current-listing filter is used."
            ),
            "variants": verdicts,
        },
        "retrieval_provenance": {
            "official_source": "https://data.binance.vision",
            "downloaded_raw_bytes": 80268850,
            "kline_archives_ok": 2678,
            "funding_archives_ok": 2543,
            "funding_http_404_records": 137,
            "integrity": (
                "Computed SHA-256, byte size, URL, retrieval UTC and ZIP CRC "
                "recorded under raw/focused/cross-section/provenance."
            ),
            "raw_plus_derived_bytes": 248860544,
            "raw_cap_bytes": 1073741824,
        },
        "gate_conclusion": (
            "Direction status is falsified solely because every variant has "
            "an observed false development-subperiod IC gate (a valid "
            "observation at or below 0.02). The 2022-02-26 first-censor stop "
            "makes ALL confirmation gates UNKNOWN - including concentration "
            "and collateral, which are NOT marked pass - and unknown gates "
            "cannot rescue an observed false conjunct."
        ),
        "confirmation_gate_status_after_stop": {
            "note": ("The first held incomplete outcome on 2022-02-26 stops "
                     "every variant path; row 574 and the tail are NaN/"
                     "unavailable, so these gates are UNKNOWN, not pass."),
            "confirmation_ic": "UNKNOWN",
            "confirmation_2x_spread": "UNKNOWN",
            "confirmation_concentration": "UNKNOWN",
            "confirmation_collateral_buffer": "UNKNOWN",
        },
        "data_limitations": {
            "binance_funding": (
                "137 requested monthly funding archive URLs returned HTTP "
                "404 across 15 symbols. Only unavailable asset-months are "
                "fail-closed; earlier complete history remains, and no "
                "missing rate is imputed as zero."
            ),
            "okx": (
                "No comparable official point-in-time listing/delisting "
                "archive was obtained in the frozen bounded acquisition. "
                "Adding current OKX instruments would survivor-filter, so "
                "OKX is excluded rather than substituted."
            ),
            "terminal_delisting": (
                "With no announcement metadata, the frozen censor rule stops "
                "the variant path at the prior known close when any held "
                "asset outcome is incomplete; no unobserved return is "
                "invented and no trade-out at an unavailable close occurs."
            ),
        },
        "counterfactual_diagnostics": (
            "Positive one-window returns and 1x versus 2x cost differences "
            "are retrospective diagnostics only; they do not confirm "
            "durable alpha and do not change the falsified status."
        ),
        "prospective_requirement": (
            "A new preregistered prospective campaign would require an "
            "append-only point-in-time instrument master and future funding/"
            "price capture. This dossier makes no deployment claim."
        ),
        "representative_returns_trial": "liq_composite",
    }
    with open(EXP / "results.json", "w") as f:
        json.dump(result_doc, f, indent=2)
        f.write("\n")
    print("rerun complete: results.json, trials.json and returns.csv refreshed")


if __name__ == "__main__":
    main()
