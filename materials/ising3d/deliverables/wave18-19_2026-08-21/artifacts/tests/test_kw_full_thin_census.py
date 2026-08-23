"""Standalone verifier for the bounded Kac--Ward thin-section census artifact.

This test intentionally imports neither the producer nor e139.  It verifies the
artifact's exact residue/schema controls, exhaustive-section accounting, anchor
reproduction, deterministic manifest, and the distinction between a completed
census and an honest process-capped preflight.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/full_thin_census.json"
PRIOR = ROOT / "results/kac_ward/components.json"
E139 = ROOT / "experiments/e139_kw_components.py"
TOTAL_SECTIONS = 4**7
SOLVE_GRID = 4**6
UNIT_VALUES = {1, 2, 3, 4}


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifest(payload: dict) -> None:
    manifest = payload["manifest"]
    stable = copy.deepcopy(payload)
    stable.pop("manifest")
    provenance = stable["provenance"]
    provenance.pop("generated_utc", None)
    provenance.pop("process_seconds", None)
    expected = hashlib.sha256(canonical_json(stable)).hexdigest()
    assert manifest["canonical_body_sha256"] == expected
    assert manifest["prior_components_sha256"] == sha256(PRIOR)
    assert manifest["e139_source_sha256"] == sha256(E139)


def check_sections(payload: dict) -> None:
    arithmetic = payload["arithmetic"]
    assert arithmetic == {
        "exact_integer_residues_only": True,
        "field_for_section_census": "F_5^x",
        "modulus": 625,
        "prime": 5,
        "solve_grid_size_per_section": SOLVE_GRID,
        "unit_values": [1, 2, 3, 4],
    }
    coverage = payload["section_coverage"]
    assert coverage["method"] == "exhaustive_4^7_unit_sections"
    assert coverage["orbit_quotient_used"] is False
    assert coverage["total_sections"] == TOTAL_SECTIONS
    assert coverage["sections_launched"] + coverage["sections_not_launched"] == TOTAL_SECTIONS
    representatives = coverage["representatives"]
    assert len(representatives) == coverage["sections_launched"]
    held_seen = set()
    for section in representatives:
        held = tuple(section["held_values_mod5"])
        assert len(held) == 7 and set(held) <= UNIT_VALUES
        assert held not in held_seen
        held_seen.add(held)
        assert section["solve_grid_size"] == SOLVE_GRID
        solutions = section["solutions"]
        assert section["solution_count_mod5"] == len(solutions)
        assert section["nonsingular_solution_count_mod5"] == sum(
            int(record["jacobian_determinant_mod5"] % 5 != 0) for record in solutions
        )
        for record in solutions:
            point = record["point13_mod5"]
            assert len(point) == 13 and set(point) <= UNIT_VALUES
            assert tuple(point[6:]) == held
            determinant = record["jacobian_determinant_mod5"]
            assert 0 <= determinant < 5
            assert len(record["jacobian_mod5"]) == 6
            assert all(len(row) == 6 and all(0 <= value < 5 for value in row)
                       for row in record["jacobian_mod5"])
            if record["lift_success_mod625"]:
                lift = record["lift_solve6_mod625"]
                assert lift is not None and len(lift) == 6
                assert all(0 <= value < 625 for value in lift)
                assert all(value % 5 == point[index] for index, value in enumerate(lift))
            else:
                assert record["lift_solve6_mod625"] is None
            assert record["holdout_status"] == "NOT_LAUNCHED"
    if coverage["all_sections_completed"]:
        assert coverage["sections_launched"] == TOTAL_SECTIONS
        assert coverage["sections_not_launched"] == 0
        assert len(held_seen) == TOTAL_SECTIONS
    else:
        assert coverage["sections_not_launched"] > 0


def check_anchor(payload: dict, prior: dict) -> None:
    anchor = payload["anchor_cross_check"]
    assert anchor["status"] == "PASS"
    assert anchor["independent_recomputed"] is True
    assert anchor["held_values_mod5"] == [1, 1, 1, 1, 1, 4, 4]
    assert anchor["solve_grid_size"] == SOLVE_GRID
    assert anchor["solution_count"] == 3
    assert anchor["solutions_mod5"] == prior["thinvariants"]["column"]["solutions"]
    assert anchor["jacobian_determinants_mod5"] == [2, 4, 1]
    assert all(anchor["checks"].values())
    holdout = anchor["holdout_residues_mod625_from_components"]
    assert holdout["orders"] == [4, 6, 8, 10, 12]
    for name in ("k8_residues", "k10_residues", "k12_residues"):
        residues = holdout[name]
        assert len(residues) == 3
        assert all(isinstance(value, int) and 0 <= value < 625 for value in residues)
    assert holdout["k8_residues"] == [250, 456, 430]
    assert holdout["k10_residues"] == [55, 350, 480]
    assert holdout["k12_residues"] == [266, 132, 286]


def check_status(payload: dict) -> None:
    assert payload["status"] in {"PREFLIGHT", "COMPLETE"}
    assert payload["claim_tag"] == "[UNRESOLVED]"
    assert "[UNRESOLVED]" in payload["headline"]
    assert payload["theorem_scope"]["timeout_is_not_theorem"] is True
    assert "[UNRESOLVED]" in payload["theorem_scope"]["full_branch"]
    holdout = payload["holdout"]
    if payload["status"] == "PREFLIGHT":
        assert payload["preflight"]["reason"]
        assert holdout["status"] == "NOT_LAUNCHED"
        assert holdout["complete_for_all_lifts"] is False
        assert payload["preflight"]["what_was_launched"]["holdout_routes"] is False
        assert payload["preflight"]["what_was_not_launched"]["full_branch_proof"] is True
    else:
        coverage = payload["section_coverage"]
        assert coverage["all_sections_completed"] is True
        assert holdout["complete_for_all_lifts"] is True


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    check_manifest(payload)
    check_sections(payload)
    check_anchor(payload, prior)
    check_status(payload)
    if payload["status"] == "COMPLETE":
        print("PASS complete thin-section census artifact")
    else:
        coverage = payload["section_coverage"]
        print(
            f"PREFLIGHT artifact verified honestly: {coverage['sections_launched']}/"
            f"{coverage['total_sections']} sections; unresolved scope preserved"
        )


if __name__ == "__main__":
    main()
