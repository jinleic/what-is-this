#!/usr/bin/env python3
"""Produce the exact branch-11111 section-orbit audit and complete F_5-unit census.

The symmetry result is an honest negative quotient theorem: after the existing
gauge and diagonal normalizations, only the identity acts on e156's seven-held
section fibration.  The complete finite census is nevertheless feasible through
e170's exact character evaluator.  A process-time cap is never promoted to a
mathematical conclusion.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import platform
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

import e169_kw_section_group as group_engine
import e170_kw_unit_census as census_engine

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/kac_ward/orbit_quotient.json"
E156_ARTIFACT = ROOT / "results/kac_ward/full_thin_census.json"
SCRIPT_PATHS = [
    ROOT / "experiments/e169_kw_section_group.py",
    ROOT / "experiments/e170_kw_unit_census.py",
    ROOT / "experiments/e171_kw_orbit_quotient.py",
]
DEFAULT_CAP_SECONDS = 120.0


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def prior_preflight() -> dict[str, Any]:
    payload = json.loads(E156_ARTIFACT.read_text(encoding="utf-8"))
    return {
        "artifact": str(E156_ARTIFACT.relative_to(ROOT)),
        "sha256": sha256_file(E156_ARTIFACT),
        "sections_launched": int(payload["section_coverage"]["sections_launched"]),
        "total_sections": int(payload["section_coverage"]["total_sections"]),
        "first_section_process_seconds": payload["preflight"]["first_section_seconds"],
        "projected_full_process_seconds": payload["preflight"]["estimated_full_section_seconds"],
        "status": payload["status"],
    }


def volatile_free(value: dict[str, Any]) -> dict[str, Any]:
    stable = copy.deepcopy(value)
    meta = stable.get("meta", {})
    meta.pop("generated_utc", None)
    meta.pop("decision_process_seconds", None)
    stable["meta"] = meta
    method = stable.get("data", {}).get("method_comparison", {})
    method.pop("new_character_enumeration_process_seconds", None)
    method.pop("full_producer_decision_process_seconds", None)
    method.pop("projected_to_character_enumeration_ratio", None)
    return stable


def evaluate_catalog_action(
    point: tuple[int, ...],
    coordinate_action: dict[str, dict[str, Any]],
    coordinate_names: list[str],
) -> tuple[int, ...]:
    image: list[int] = []
    for name in coordinate_names:
        record = coordinate_action[name]
        value = int(record["coefficient"]) % 5
        for coordinate, exponent in enumerate(record["exponents"]):
            if exponent:
                value = value * pow(point[coordinate], int(exponent) % 4, 5) % 5
        image.append(value)
    return tuple(image)


def construction_fiber_counterexamples(
    e139,
    symmetry: dict[str, Any],
    solution_indices: list[int],
) -> list[dict[str, Any]]:
    """Disprove section factorization even after restriction to the construction zero set."""
    points = sorted(census_engine.decode_point(index) for index in solution_indices)
    point_set = set(points)
    counterexamples: list[dict[str, Any]] = []
    for element in symmetry["cubic_point_group"]["elements"]:
        action = element["coordinate_action"]
        seen: dict[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]] = {}
        witness = None
        for point in points:
            image = evaluate_catalog_action(point, action, list(e139.NONDIAGONAL_NAMES))
            if image not in point_set:
                raise AssertionError("catalog action did not preserve the exact F_5 construction zero set")
            if element["induces_section_map"]:
                continue
            section = point[6:]
            target = image[6:]
            if section in seen and seen[section][1] != target:
                first_point, first_target = seen[section]
                first_image = evaluate_catalog_action(
                    first_point, action, list(e139.NONDIAGONAL_NAMES)
                )
                witness = {
                    "axis_permutation": element["axis_permutation"],
                    "axis_signs": element["axis_signs"],
                    "input_section": list(section),
                    "first_point13_mod5": list(first_point),
                    "second_point13_mod5": list(point),
                    "first_image13_mod5": list(first_image),
                    "second_image13_mod5": list(image),
                    "distinct_target_sections": [list(first_target), list(target)],
                }
                break
            seen[section] = (point, target)
        if not element["induces_section_map"]:
            if witness is None:
                raise AssertionError("missing construction-variety fiber counterexample")
            counterexamples.append(witness)
    return counterexamples




def make_payload(cap_seconds: float, chunk_size: int) -> dict[str, Any]:
    started = time.process_time()
    deadline = started + cap_seconds
    e139 = group_engine.load_e139()
    construction = group_engine.build_construction(e139)
    symmetry = group_engine.analyze_group(e139, construction)
    census = census_engine.run_complete_census(e139, construction, deadline, chunk_size)
    solution_indices = census["enumeration"]["solution_indices"]
    variety_counterexamples = (
        construction_fiber_counterexamples(e139, symmetry, solution_indices)
        if census["enumeration"]["complete"]
        else []
    )
    symmetry["section_action"]["construction_variety_counterexamples_mod5"] = variety_counterexamples
    symmetry["checks"].append({
        "name": "construction_variety_section_factorization_counterexamples",
        "passed": len(variety_counterexamples) == 15,
        "detail": len(variety_counterexamples),
    })
    decision_seconds = time.process_time() - started
    prior = prior_preflight()

    point_records = [] if census["holdout"] is None else census["holdout"]["point_records"]
    sections = census["sections"]
    disposition_counts = Counter(record["disposition"] for record in point_records)
    status_counts = (
        {"EMPTY_OVER_Q": 0, "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT": 0, "UNDECIDED": 16_384}
        if sections is None
        else sections["status_counts"]
    )
    remaining_indices = list(range(16_384)) if sections is None else sections["undecided_lex_indices"]
    classification_rows = [] if sections is None else sections["classification_rows"]
    anchor_row = next(
        (row for row in classification_rows if row["section"] == [1, 1, 1, 1, 1, 4, 4]),
        None,
    )

    checks: list[dict[str, Any]] = []
    checks.extend(symmetry["checks"])
    checks.extend([
        {
            "name": "hard_process_time_cap_honored",
            "passed": decision_seconds <= cap_seconds,
            "detail": {"decision_seconds": decision_seconds, "cap_seconds": cap_seconds},
        },
        {
            "name": "complete_4_to_13_unit_grid",
            "passed": bool(census["enumeration"]["complete"]),
            "detail": {
                "next_global_index": census["enumeration"]["next_global_index"],
                "total_points": census_engine.TOTAL_POINTS,
            },
        },
        {
            "name": "exact_solution_count",
            "passed": len(census["enumeration"]["solution_indices"]) == 2_960,
            "detail": len(census["enumeration"]["solution_indices"]),
        },
        {
            "name": "all_local_points_have_exact_disposition",
            "passed": len(point_records) == 2_960 and all(
                record["disposition"].startswith("KILLED_") for record in point_records
            ),
            "detail": dict(sorted(disposition_counts.items())),
        },
        {
            "name": "three_holdout_certificate_mechanisms_present",
            "passed": all(
                disposition_counts[name] > 0
                for name in (
                    "KILLED_H8_MOD5",
                    "KILLED_H12_MOD5",
                    "KILLED_H8_MOD625_NONSINGULAR",
                    "KILLED_BY_SINGULAR_LIFT_OBSTRUCTION",
                )
            ),
            "detail": dict(sorted(disposition_counts.items())),
        },
        {
            "name": "complete_section_classification",
            "passed": len(classification_rows) == 16_384,
            "detail": len(classification_rows),
        },
        {
            "name": "section_status_partition",
            "passed": status_counts == {
                "EMPTY_OVER_Q": 0,
                "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT": 2_437,
                "UNDECIDED": 13_947,
            },
            "detail": status_counts,
        },
        {
            "name": "anchor_section_control",
            "passed": anchor_row is not None
            and anchor_row["solution_count_mod5"] == 3
            and anchor_row["status"] == "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT",
            "detail": anchor_row,
        },
        {
            "name": "no_empty_over_Q_overclaim",
            "passed": status_counts["EMPTY_OVER_Q"] == 0,
            "detail": (
                "No missing F_5-unit point and no process-time event is labelled EMPTY_OVER_Q."
            ),
        },
        {
            "name": "prior_preflight_control",
            "passed": prior["sections_launched"] == 1
            and prior["total_sections"] == 16_384
            and prior["status"] == "PREFLIGHT",
            "detail": prior,
        },
    ])

    complete = census["complete"] and all(check["passed"] for check in checks)
    projected = prior["projected_full_process_seconds"]
    character_seconds = census["enumeration"]["process_seconds"]
    payload: dict[str, Any] = {
        "schema_version": 1,
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e171_kw_orbit_quotient.py",
            "support_modules": [
                "experiments/e169_kw_section_group.py",
                "experiments/e170_kw_unit_census.py",
            ],
            "interpreter_contract": ".venv/bin/python",
            "sys_executable": sys.executable,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "sympy_version": sp.__version__,
            "numpy_version": np.__version__,
            "arithmetic": "exact integer/Laurent and exact residues modulo 5,25,125,625",
            "floating_point_decides_no_statement": True,
            "process_time_cap_seconds": cap_seconds,
            "decision_process_seconds": decision_seconds,
            "chunk_size": chunk_size,
            "peak_rss_design_bound": "<4 GB; largest census chunk arrays are <100 MB",
            "source_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path) for path in SCRIPT_PATHS
            }
            | {str(group_engine.E139_PATH.relative_to(ROOT)): sha256_file(group_engine.E139_PATH)},
        },
        "data": {
            "status": "COMPLETE_F5_UNIT_CENSUS_BRANCH_UNRESOLVED" if complete else "PARTIAL_CAP_OR_CHECK_FAILURE",
            "headline": (
                "[THEOREM] On the diagonal-normalized 13-coordinate thin chart, the exact "
                "group acting on the 4^7 e156 sections is the identity, so there are 16384 "
                "singleton orbits and no quotient reduction. [COMPUTATION] The complete F_5-unit "
                "census has 2960 construction points in 2437 sections, and every recorded local "
                "point fails an independent order-8 or order-12 holdout by an exact residue or "
                "finite lift obstruction. [UNRESOLVED] This is not branch-11111 emptiness over Q."
            ),
            "symmetry": symmetry,
            "census": {
                "field": "F_5^x",
                "unit_values": list(census_engine.UNIT_VALUES),
                "coordinate_count": 13,
                "solve_coordinate_count": 6,
                "held_coordinate_count": 7,
                "total_grid_points": census_engine.TOTAL_POINTS,
                "total_sections": census_engine.TOTAL_SECTIONS,
                "solve_grid_size_per_section": census_engine.SOLVE_GRID_SIZE,
                "enumeration": census["enumeration"],
                "holdout_orders": [8, 10, 12],
                "holdout_log_P_coefficients": {str(key): value for key, value in census_engine.HOLDOUT_LOG_P.items()},
                "holdout_classification": census["holdout"],
                "section_classification": None
                if sections is None
                else {key: value for key, value in sections.items() if key != "classification_rows"},
            },
            "classification_table": classification_rows,
            "classification_summary": {
                "status_counts": status_counts,
                "point_disposition_counts": dict(sorted(disposition_counts.items())),
                "remaining_section_count": len(remaining_indices),
                "remaining_undecided_lex_indices": remaining_indices,
                "remaining_definition": (
                    "lex_index is base-4 on the seven-tuple in lexicographic value order "
                    "1,2,3,4; exactly the listed indices have status UNDECIDED"
                ),
                "empty_over_Q_certificates": [],
            },
            "scope": {
                "proved": (
                    "[THEOREM] For each of the 16384 slices whose seven held coordinates are "
                    "fixed exactly at a tuple in {1,2,3,4}^7, no Q_5 point with all six solve "
                    "coordinates 5-adic units satisfies both the 42 construction equations and "
                    "the recorded order-8/10/12 holdouts: 13947 slices have no F_5-unit "
                    "construction reduction, and every reduction in the other 2437 is killed."
                ),
                "unresolved": (
                    "[UNRESOLVED] Other higher 5-adic digits in the seven held coordinates, "
                    "non-unit 5-adic valuations, coordinates with zero reductions, chart-degenerate "
                    "points, Q_p for p other than 5, and characteristic-zero points are outside "
                    "this census; branch 11111 therefore remains open."
                ),
                "nullstellensatz_route": (
                    "[THEOREM] The construction ideal is proper (the prior anchor Q_5 construction "
                    "point is an exact witness), so no rational Nullstellensatz certificate can "
                    "prove that construction ideal empty at any multiplier degree."
                ),
                "timeout_is_not_theorem": True,
            },
            "method_comparison": {
                "prior": prior,
                "new_character_enumeration_process_seconds": character_seconds,
                "full_producer_decision_process_seconds": decision_seconds,
                "projected_to_character_enumeration_ratio": None
                if not projected or not character_seconds
                else projected / character_seconds,
                "earlier_complete_run_before_variety_counterexample_extension": {
                    "character_enumeration_process_seconds": 8.650371,
                    "full_producer_decision_process_seconds": 24.405328,
                    "projected_to_character_enumeration_ratio": 11560.50157663758,
                    "scope": (
                        "[COMPUTATION] Grounded prior PASS run of the same complete census, before "
                        "the final artifact added the 15 variety-fibration counterexample records."
                    ),
                },
                "statement": (
                    "[COMPUTATION] Exact finite-group character evaluation replaces e156's "
                    "per-section SymPy substitution loop; process times are measured host-specific "
                    "observations, not mathematical claims."
                ),
            },
            "checks": checks,
        },
    }
    stable = volatile_free(payload)
    payload["meta"]["canonical_body_sha256"] = hashlib.sha256(canonical_json(stable)).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap-seconds", type=float, default=DEFAULT_CAP_SECONDS)
    parser.add_argument("--chunk-size", type=int, default=census_engine.DEFAULT_CHUNK_SIZE)
    args = parser.parse_args()
    if args.cap_seconds <= 0 or args.chunk_size <= 0:
        raise ValueError("cap and chunk size must be positive")
    payload = make_payload(args.cap_seconds, args.chunk_size)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    checks = payload["data"]["checks"]
    if all(check["passed"] for check in checks) and payload["data"]["status"].startswith("COMPLETE"):
        print(
            f"orbit_count={payload['data']['symmetry']['section_action']['orbit_count']} "
            f"solutions={len(payload['data']['census']['holdout_classification']['point_records'])} "
            f"statuses={payload['data']['classification_summary']['status_counts']}"
        )
        print(f"artifact={OUTPUT}")
        print("PASS")
    else:
        print(f"PARTIAL artifact={OUTPUT}")


if __name__ == "__main__":
    main()
