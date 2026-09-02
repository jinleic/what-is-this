#!/usr/bin/env python3
"""Independent postcompute validator for the frozen MM3 census artifacts."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CAMPAIGN = "2026-08-31T083323Z_55447c11-5887-4fdd-aa2f-30e31f2332e8_b3046e6a6c50"
CODE_SHA256 = "b3046e6a6c506ffd5c78871a9f057770a552968314ac6fd225ad60e2b1d4e77c"
DECOMPOSITIONS = ("mws59", "stapleton60")
MONOMIAL_NODE = 125
EXPECTED_SURVIVORS = {"mws59": 5796, "stapleton60": 4800}
EXPECTED_ALL_MONOMIAL = {"mws59": 56, "stapleton60": 58}
EXPECTED_PUBLISHED = {
    "paper55": 55,
    "perminov58": 58,
    "sun56": 56,
    "mws59": 59,
    "stapleton60": 60,
}
EXPECTED_CERTIFIED = {
    "paper55": [([12, 13, 13], [False, False, False], 55),
                ([13, 13, 12], [False, False, False], 55),
                ([13, 12, 13], [False, False, False], 55)],
    "sun56": [([12, 11, 16], [False, False, True], 55),
              ([11, 16, 12], [False, True, False], 55),
              ([16, 12, 11], [True, False, False], 55)],
    "mws59": [([13, 12, 14], [False, False, False], 56),
              ([12, 14, 13], [False, False, False], 56),
              ([14, 13, 12], [False, False, False], 56)],
    "stapleton60": [([14, 14, 13], [False, False, False], 58),
                    ([14, 13, 14], [False, False, False], 58),
                    ([13, 14, 14], [False, False, False], 58)],
}
CORE_OUTPUTS = (
    "run.log",
    "phaseA_controls.json",
    "decision_mws59.npz",
    "decision_mws59_summary.json",
    "decision_stapleton60.npz",
    "decision_stapleton60_summary.json",
    "pair_tables.npz",
    "campaign_result.json",
)

checks: list[str] = []


def check(condition: bool, name: str) -> None:
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(name: str) -> dict[str, object]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def verify_checksum_ledger() -> None:
    ledger = HERE / "checksums_precompute.sha256"
    count = 0
    for line in ledger.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        expected, relative = line.split("  ", 1)
        check(sha256(HERE / relative) == expected, f"precompute hash {relative}")
        count += 1
    check(count == 17, "17 precompute checksum rows")


def validate_phase_a(phase: dict[str, object]) -> None:
    check(
        phase["ternary_matrix_census"] == {
            "ternary_matrices": 19683,
            "invertible": 11808,
            "unimodular": 6960,
            "monomial": 48,
            "nonmonomial_unimodular": 6912,
            "inverse_ternary": 4656,
            "inverse_ternary_nonmonomial": 4608,
        },
        "ternary matrix census",
    )
    factor = phase["data_factorization"]
    check(factor["data_nodes"] == 145, "145 data nodes")
    check(factor["right_orbit_size"] == 48, "right orbit size 48")
    check(factor["partition_size"] == 6960, "factorization partition 6960")
    check(factor["unique_factorization"] is True, "unique factorization")
    check(factor["monomial_data_node"] == MONOMIAL_NODE, "monomial node 125")

    anchors = phase["anchors"]
    for name, total in EXPECTED_PUBLISHED.items():
        row = anchors["published"][name]
        check(row["published_total"] == total, f"{name} published total")
        check(row["ternary"] is True, f"{name} ternary")
        check(
            row["brent_int"] == row["brent_fmpz"] == row["brent_frozen"] == 729,
            f"{name} Brent 729",
        )
    perminov = anchors["published"]["perminov58"]["count_provenance"]
    check(perminov["decoded_factor_mismatches"] == 0, "Perminov decoded factors")
    check(perminov["complexity"]["reduced"] == 58, "Perminov reduced complexity")

    for name, expected_rows in EXPECTED_CERTIFIED.items():
        rows = anchors["certified"][name]["sigma_classes"]
        actual = [(row["d"], row["floor_exists"], row["total_lb"]) for row in rows]
        check(actual == expected_rows, f"{name} certified sigma table")
    check(anchors["mandatory_55_anchor_asserted_before_new_counts"] is True,
          "mandatory 55 anchor order")
    check(phase["canonical_data_node_anchors"] == {
        "paper55": 55, "sun56": 55, "mws59": 56, "stapleton60": 58
    }, "canonical data-node anchors")
    check(anchors["frozen_vs_standard_letter_controls"]["checks"] == 25,
          "25 convention checks")

    controls = phase["semantic_controls"]
    check(controls["honest_final_family"] == {
        "brent_int": 729, "brent_fmpz": 729, "ternary": True
    }, "honest semantic plant")
    for key in ("wrong_W_plant", "transpose_as_inverse_plant",
                "direct_brent_corruption_plant"):
        check(controls[key]["brent_failures"] > 0, f"{key} detected")
    check(controls["nonunimodular_det2_rejected"] is True, "det-2 plant rejected")
    expected_side = {"nonzero_rows": 23, "rank_Q": 9, "nonzero_columns": 9}
    for name, rows in controls["activity"].items():
        check(len(rows) == 3, f"{name} three activity sigma rows")
        for row in rows:
            check(row["ok"] is True, f"{name} activity ok sigma {row['sigma']}")
            check(all(side == expected_side for side in row["sides"]),
                  f"{name} activity side invariants sigma {row['sigma']}")


def expected_survivor_rows(U: np.ndarray, V: np.ndarray,
                           W: np.ndarray) -> np.ndarray:
    rows: list[tuple[int, int, int]] = []
    for p in range(145):
        for q in np.flatnonzero(U[p]):
            for r in np.flatnonzero(V[q] & W[p]):
                rows.append((p, int(q), int(r)))
    return np.asarray(rows, dtype=np.uint16)


def validate_decomposition(name: str, summary: dict[str, object],
                           census: dict[str, object], pairs: np.lib.npyio.NpzFile
                           ) -> dict[str, int]:
    expected_names = {
        f"{side}_{suffix}" for side in "UVW"
        for suffix in ("d", "floor", "lb", "states")
    } | {"survivors", "side_lbs", "totals"}
    with np.load(HERE / f"decision_{name}.npz", allow_pickle=False) as decision:
        check(set(decision.files) == expected_names, f"{name} 15 decision arrays")
        for side in "UVW":
            mask = pairs[f"{name}_{side}_ternary"]
            check(mask.shape == (145, 145) and mask.dtype == np.bool_,
                  f"{name} {side} pair mask schema")
            for suffix, dtype in (("d", np.int16), ("floor", np.int8),
                                  ("lb", np.int16), ("states", np.int64)):
                array = decision[f"{side}_{suffix}"]
                check(array.shape == (145, 145) and array.dtype == dtype,
                      f"{name} {side}_{suffix} schema")
                check(np.all(array[~mask] == -1), f"{name} {side}_{suffix} sentinels")
            floor = decision[f"{side}_floor"]
            d = decision[f"{side}_d"]
            lb = decision[f"{side}_lb"]
            check(np.all((floor[mask] == 0) | (floor[mask] == 1)),
                  f"{name} {side} floor range")
            check(np.array_equal(lb[mask], d[mask] + (floor[mask] == 0)),
                  f"{name} {side} lower-bound formula")
            check(np.all(decision[f"{side}_states"][mask] > 0),
                  f"{name} {side} DFS states positive")

        survivors = decision["survivors"]
        side_lbs = decision["side_lbs"]
        totals = decision["totals"]
        count = EXPECTED_SURVIVORS[name]
        check(survivors.shape == (count, 3) and survivors.dtype == np.uint16,
              f"{name} survivor schema")
        check(side_lbs.shape == (count, 3) and side_lbs.dtype == np.int16,
              f"{name} side-lb schema")
        check(totals.shape == (count,) and totals.dtype == np.int16,
              f"{name} totals schema")
        check(np.array_equal(survivors, pairs[f"{name}_survivors"]),
              f"{name} survivor files agree")
        U = pairs[f"{name}_U_ternary"]
        V = pairs[f"{name}_V_ternary"]
        W = pairs[f"{name}_W_ternary"]
        check(np.array_equal(survivors, expected_survivor_rows(U, V, W)),
              f"{name} complete lexicographic survivor enumeration")
        indexed = np.column_stack((
            decision["U_lb"][survivors[:, 0], survivors[:, 1]],
            decision["V_lb"][survivors[:, 1], survivors[:, 2]],
            decision["W_lb"][survivors[:, 0], survivors[:, 2]],
        ))
        check(np.array_equal(side_lbs, indexed), f"{name} side-lb indexing")
        check(np.array_equal(totals, side_lbs.sum(axis=1) + 14),
              f"{name} +14 total formula")
        all_monomial_rows = np.flatnonzero(np.all(survivors == MONOMIAL_NODE, axis=1))
        check(len(all_monomial_rows) == 1, f"{name} unique all-monomial row")
        mono_index = int(all_monomial_rows[0])
        check(int(totals[mono_index]) == EXPECTED_ALL_MONOMIAL[name],
              f"{name} all-monomial total")

        histogram = {str(k): v for k, v in sorted(Counter(map(int, totals)).items())}
        minimizers = np.flatnonzero(totals == totals.min()).astype(int).tolist()
        nonmonomial = np.any(survivors != MONOMIAL_NODE, axis=1)
        check(summary["decomposition"] == name, f"{name} summary name")
        check(summary["decision_formula"] ==
              "U_lb + V_lb + Wfac_lb + 14, side_lb=d+[floor impossible]",
              f"{name} decision formula")
        check(summary["survivors_decided"] == count, f"{name} summary survivors")
        check(summary["histogram"] == histogram, f"{name} histogram")
        check(summary["min_total_lb"] == int(totals.min()), f"{name} minimum")
        check(summary["max_total_lb"] == int(totals.max()), f"{name} maximum")
        check(summary["upper_median_total_lb"] == int(np.sort(totals)[count // 2]),
              f"{name} upper median")
        check(summary["minimizer_count_data_triples"] == len(minimizers),
              f"{name} minimizer count")
        check(summary["minimizer_survivor_indices"] == minimizers,
              f"{name} minimizer indices")
        check(summary["nonmonomial_involving_min_total_lb"] ==
              int(totals[nonmonomial].min()), f"{name} nonmonomial minimum")
        check(summary["flags_le_54"] == int(np.count_nonzero(totals <= 54)) == 0,
              f"{name} zero <=54")
        check(summary["all_monomial_data_triple"] == {
            "index": mono_index,
            "total_lb": EXPECTED_ALL_MONOMIAL[name],
            "asserted_before_histogram": True,
        }, f"{name} all-monomial summary")

        for side in "UVW":
            mask = pairs[f"{name}_{side}_ternary"]
            row = summary["pair_decisions"][side]
            check(row["ternary_pairs"] == int(mask.sum()),
                  f"{name} {side} pair count")
            check(row["floor_exists"] + row["floor_impossible_exact_census"] ==
                  row["ternary_pairs"], f"{name} {side} floor census partition")
            check(row["floor_exists"] ==
                  int(np.count_nonzero(decision[f"{side}_floor"][mask] == 1)),
                  f"{name} {side} floor exists count")
            d_hist = {str(k): v for k, v in
                      sorted(Counter(map(int, decision[f"{side}_d"][mask])).items())}
            check(row["d_histogram"] == d_hist, f"{name} {side} d histogram")

        check(len(summary["independent_pair_crosschecks"]) == 18,
              f"{name} 18 pair crosschecks")
        for row in summary["independent_pair_crosschecks"]:
            side = row["side"]
            a, b = row["a"], row["b"]
            check(row["d"] == int(decision[f"{side}_d"][a, b]),
                  f"{name} pair crosscheck d")
            check(row["floor_exists"] == bool(decision[f"{side}_floor"][a, b]),
                  f"{name} pair crosscheck floor")
        check(len(summary["fresh_total_crosschecks"]) == 24,
              f"{name} 24 total crosschecks")
        for index, row in enumerate(summary["fresh_total_crosschecks"]):
            survivor_index = row["survivor_index"]
            total = int(totals[survivor_index])
            check(row["total_lb"] == total,
                  f"{name} total crosscheck {index}")
            if index < 12:
                check(row["fiber_total_lb"] == total,
                      f"{name} fiber crosscheck {index}")
                check(row["sigma_totals"] == [total, total, total],
                      f"{name} sigma crosscheck {index}")
            else:
                check(set(row) == {"survivor_index", "total_lb"},
                      f"{name} direct-only crosscheck schema {index}")
        check(len(summary["survivor_brent_crosschecks"]) == 12,
              f"{name} 12 Brent crosschecks")
        check(all(row["brent_int"] == row["brent_fmpz"] == 729
                  for row in summary["survivor_brent_crosschecks"]),
              f"{name} Brent crosschecks 729")

    check(census["survivor_data_triples"] == count, f"{name} census survivors")
    check(census["independent_dense_count"] == count, f"{name} dense count")
    check(sum(census["by_monomial_node_count"].values()) == count,
          f"{name} monomial bins sum")
    check(census["by_monomial_node_count"]["3"] == 1,
          f"{name} one triple-monomial survivor")
    check(census["fiber_multiplier_per_sigma"] == 110592,
          f"{name} fiber multiplier")
    check(census["sigma_multiplier"] == 3, f"{name} sigma multiplier")
    check(census["sandwich_instances_per_sigma"] == count * 110592,
          f"{name} per-sigma instances")
    check(census["sandwich_instances_all_sigma"] == count * 331776,
          f"{name} all-sigma instances")
    check(census["full_raw_domain_per_decomposition"] == 1011460608000,
          f"{name} raw domain")
    pair_info = census["pair_tables"]
    check(pair_info["int64_hard_bounds"] == {"U": 36, "V": 18, "W": 9},
          f"{name} int64 hard bounds")
    check(pair_info["int64_overflow_impossible"] is True,
          f"{name} overflow impossible")
    check(all(pair_info["max_abs_observed"][side] <= bound
              for side, bound in pair_info["int64_hard_bounds"].items()),
          f"{name} observed values within bounds")
    check(census["controls"]["anchored_subcube"]["triples"] == 21025,
          f"{name} anchored subcube size")
    check(census["controls"]["anchored_subcube"]["mismatches"] == 0,
          f"{name} anchored subcube agreement")
    check(census["controls"]["anchored_subcube"]["direct_positives"] ==
          census["controls"]["anchored_subcube"]["table_positives"],
          f"{name} anchored positive count")
    check(census["controls"]["random_concordance"]["trials"] == 200 and
          census["controls"]["random_concordance"]["mismatches"] == 0,
          f"{name} random concordance")
    check(census["controls"]["numpy_vs_python_entry_checks"] == 24,
          f"{name} entry checks")
    check(census["controls"]["right_monomial_fiber_predicate"] ==
          {"trials": 25, "mismatches": 0}, f"{name} fiber checks")
    return {"survivors": count, "valid": count * 331776,
            "minimum": summary["min_total_lb"],
            "nonmonomial_minimum": summary["nonmonomial_involving_min_total_lb"]}


def validate_result(result: dict[str, object], summaries: dict[str, dict[str, object]],
                    headlines: dict[str, dict[str, int]]) -> None:
    check(result["campaign"] == CAMPAIGN, "campaign name")
    check(result["code_sha256"] == CODE_SHA256, "campaign code hash")
    check(result["domain"] == {
        "decompositions": ["mws59", "stapleton60"],
        "T_size": 6960,
        "data_nodes": 145,
        "full_raw_scheme_instances": 2022921216000,
        "data_triples_enumerated": 6097250,
        "fiber_multiplier_per_sigma": 110592,
        "sigma_multiplier": 3,
        "alphabet": [-1, 0, 1],
        "action_standard": "U'=P^-1 U Q^-T; V'=Q^T V R^-T; W'=P^T W R",
    }, "campaign domain")
    for name in DECOMPOSITIONS:
        check(result["outputs"][name]["decision"] == summaries[name],
              f"{name} embedded decision")
    survivors = sum(row["survivors"] for row in headlines.values())
    valid = sum(row["valid"] for row in headlines.values())
    check(result["aggregate"] == {
        "survivor_data_triples": survivors,
        "valid_scheme_instances_all_sigma": valid,
        "min_total_lb": min(row["minimum"] for row in headlines.values()),
        "nonmonomial_involving_min_total_lb": min(
            row["nonmonomial_minimum"] for row in headlines.values()),
        "flags_le_54": 0,
    }, "campaign aggregate")
    check(result["evidence"] ==
          "MACHINE-VERIFIED exact integer census + complete finite subset DFS; no float",
          "campaign evidence label")


def validate_log() -> None:
    lines = (HERE / "run.log").read_text(encoding="utf-8").splitlines()
    needles = [
        "MANDATORY STARTUP ANCHOR PASS: paper55=55, sun56=55, mws59=56, stapleton60=58",
        "PHASE A PASS",
        "mws59: census controls pass",
        "mws59: decision complete, min=56 max=74 flags<=54=0",
        "stapleton60: census controls pass",
        "stapleton60: decision complete, min=58 max=76 flags<=54=0",
        "RUN PASS: valid orientations=3515498496, aggregate min_lb=56, zero <=54",
    ]
    positions = []
    for needle in needles:
        matches = [i for i, line in enumerate(lines) if needle in line]
        check(len(matches) == 1, f"run log line: {needle}")
        positions.append(matches[0])
    check(positions == sorted(positions), "run log ordering")
    check(not any("Traceback" in line or "RecordFlag" in line for line in lines),
          "run log has no escalation trace")


def main() -> None:
    verify_checksum_ledger()
    check(all((HERE / name).is_file() for name in CORE_OUTPUTS), "eight core outputs exist")
    validate_log()
    phase = load_json("phaseA_controls.json")
    validate_phase_a(phase)
    summaries = {name: load_json(f"decision_{name}_summary.json")
                 for name in DECOMPOSITIONS}
    result = load_json("campaign_result.json")
    with np.load(HERE / "pair_tables.npz", allow_pickle=False) as pairs:
        check(pairs["data_reps"].shape == pairs["data_inverses"].shape == (145, 3, 3),
              "pair-table representative schemas")
        check(pairs["data_reps"].dtype == pairs["data_inverses"].dtype == np.int8,
              "pair-table representative dtypes")
        headlines = {
            name: validate_decomposition(
                name, summaries[name], result["outputs"][name]["census"], pairs
            ) for name in DECOMPOSITIONS
        }
    validate_result(result, summaries, headlines)
    output_hashes = {name: sha256(HERE / name) for name in CORE_OUTPUTS}
    report = {
        "schema": "mm3-postcompute-validation-v1",
        "status": "PASS",
        "checks_passed": len(checks),
        "core_output_hashes": output_hashes,
        "headline": headlines,
        "aggregate": result["aggregate"],
        "escalation_triggered": False,
        "timing_fields_not_used_for_adjudication": True,
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
