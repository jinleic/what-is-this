"""Regression checks for EXP-056's reciprocal-pole orbit distance proof."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp056", ROOT / "experiments" / "exp056_odd_distance.py"
)
assert _SPEC and _SPEC.loader
E56 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E56
_SPEC.loader.exec_module(E56)

CERTIFICATE = (
    ROOT / "results" / "certificates" / "exp056_wm_162_8_14_distance.json"
)


def test_target_rebuild_and_pole_cover() -> None:
    problem = E56.build_problem(E56.TARGET)
    assert problem["n"] == 162
    assert problem["k"] == 8
    assert not np.any(problem["HX"] @ problem["HZ"].T % 2)

    cover = E56.build_orbit_cover(problem)
    assert cover["complete"] is True
    assert cover["required_sectors"] == 2
    assert cover["pole_dimension"] == 4
    assert cover["union_quotient_rank_numpy"] == 8
    assert cover["union_quotient_rank_bitset"] == 8
    assert all(item["orbit_quotient_rank"] == 4 for item in cover["sectors"])
    assert all(item["stabilizer_size"] == 9 for item in cover["sectors"])
    assert all(len(item["anchor_coordinates"]) == 18 for item in cover["sectors"])


def test_target_logical_class_orbit_cover_is_complete() -> None:
    problem = E56.build_problem(E56.TARGET)
    pole_cover = E56.build_orbit_cover(problem)
    class_cover = E56.build_class_orbit_cover(problem, pole_cover)
    assert class_cover["complete"] is True
    assert class_cover["nonzero_classes"] == 255
    assert class_cover["automorphism_group_size"] == 162
    assert class_cover["class_orbits"] == 20
    assert sum(item["orbit_size"] for item in class_cover["representatives"]) == 255
    assert {item["orbit_size"] for item in class_cover["representatives"]} == {6, 9, 18}
    assert all(
        item["orbit_stabilizer_exact"]
        for item in class_cover["representatives"]
    )


def test_class_certificate_protocol_is_hash_bound() -> None:
    problem = E56.build_problem(E56.TARGET)
    pole_cover = E56.build_orbit_cover(problem)
    class_cover = E56.build_class_orbit_cover(problem, pole_cover)
    protocol = E56._class_protocol(problem, pole_cover, class_cover)
    identity = E56._class_identity(problem, pole_cover, class_cover, 0)
    assert protocol["requested_exclusion_cap"] == 13
    assert protocol["effective_solver_cap"] == 12
    assert protocol["solver"] == "kissat404"
    assert identity["instance_sha256"] == E56._class_instance_digest(
        problem, class_cover["representatives"][0], protocol
    )

    forged = {
        "identity": identity,
        "decision": {"status": "UNSAT"},
        "replay": {
            "status": "UNSAT",
            "instance_sha256": "wrong",
            "encoding_version": protocol["encoding_version"],
        },
    }
    assert E56._class_replay_valid(forged) is False
    forged["replay"]["instance_sha256"] = identity["instance_sha256"]
    assert E56._class_replay_valid(forged) is True



def test_persisted_exact_distance_certificate() -> None:
    assert CERTIFICATE.exists(), f"required artifact missing: {CERTIFICATE}"
    certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    assert E56.validate_exact_certificate_payload(certificate) is True
    assert certificate["schema"] == E56.SCHEMA
    assert certificate["verdict"] == {
        "classification": "CERTIFIED_EXACT",
        "d": 14,
        "d_X": 14,
        "d_Z": 14,
        "exact": True,
        "lower_route": "class_orbits",
    }
    lower = certificate["lower_bound"]["class_route"]
    assert lower["required_classes"] == 20
    assert lower["missing_classes"] == []
    assert lower["statuses"] == ["UNSAT"] * 20
    assert lower["all_classes_unsat"] is True
    assert lower["all_classes_replayed"] is True
    pole_route = certificate["lower_bound"]["pole_sector_route"]
    assert pole_route["canonical_route_enabled"] is False
    assert pole_route["all_orbit_sectors_replayed"] is False
    assert certificate["upper_bound"]["weight"] == 14

    problem = E56.build_problem(E56.TARGET)
    pole_cover = E56.build_orbit_cover(problem)
    class_cover = E56.build_class_orbit_cover(problem, pole_cover)
    for class_index in range(20):
        record = E56._validate_class_record(
            problem, pole_cover, class_cover, class_index
        )
        assert record is not None
        assert E56._class_replay_valid(record)


def test_target_component_decomposition_is_exact() -> None:
    problem = E56.build_problem(E56.TARGET)
    pole_cover = E56.build_orbit_cover(problem)
    components = E56.build_component_decomposition(problem, pole_cover)
    assert components["complete"] is True
    assert components["constant_stabilizer_rank"] == 27
    assert components["even_stabilizer_rank"] == 50
    assert components["transformed_stabilizer_rank"] == 77
    assert components["direct_sum_rank"] == 77
    assert components["logical_constant_projection_zero"] is True
    assert components["even_logical_quotient_rank"] == 8


def test_target_kernel_words_have_even_weight() -> None:
    problem = E56.build_problem(E56.TARGET)
    proof = E56.kernel_weight_parity_proof(problem, 13)
    assert proof["valid"] is True
    assert proof["all_column_degrees_odd"] is True
    assert proof["requested_cap"] == 13
    assert proof["effective_even_cap"] == 12

def test_target_bb_duality_is_exact_weight_preserving_isometry() -> None:
    problem = E56.build_problem(E56.TARGET)
    proof = E56.bb_duality_proof(problem)
    assert proof["valid"] is True
    assert proof["permutation_is_bijection"] is True
    assert proof["permutation_is_involution"] is True
    assert proof["maps_HX_to_HZ_exactly"] is True
    assert proof["maps_HZ_to_HX_exactly"] is True


def test_target_weight_fourteen_witness_is_physical_on_both_sides() -> None:
    problem = E56.build_problem(E56.TARGET)
    witness = E56.default_witness(problem)
    checked = E56.verify_upper_witness(problem, witness)
    assert checked["weight"] == 14
    assert checked["z_valid_numpy"] is True
    assert checked["z_valid_bitset"] is True
    assert checked["x_valid_numpy"] is True
    assert checked["x_valid_bitset"] is True
    assert checked["duality_preserves_weight"] is True


def test_orbit_reduction_matches_full_decision_on_small_odd_code() -> None:
    # Postema--Kokkelmans' (3,3) [[18,4,2]] row: tiny but nontrivial.
    record = {
        "name": "small-control",
        "ell": 3,
        "m": 3,
        "A": [(0, 0), (0, 1), (0, 2)],
        "B": [(0, 0), (1, 0), (0, 1)],
        "expected_k": 4,
        "expected_d": 2,
    }
    problem = E56.build_problem(record)
    cover = E56.build_orbit_cover(problem)
    assert cover["complete"] is True

    below = E56.decide_orbit_bundle(problem, cover, 1, conflict_budget=0)
    at = E56.decide_orbit_bundle(problem, cover, 2, conflict_budget=0)
    assert below["status"] == "UNSAT"
    assert at["status"] == "SAT"
    assert at["weight"] == 2
    assert at["verification"]["valid"] is True

    monolithic_below = E56.decide_pole_monolithic(
        problem, cover, 1, conflict_budget=0
    )
    monolithic_at = E56.decide_pole_monolithic(
        problem, cover, 2, conflict_budget=0
    )
    assert monolithic_below["status"] == "UNSAT"
    assert monolithic_at["status"] == "SAT"
    assert monolithic_at["weight"] == 2
    assert monolithic_at["verification"]["valid"] is True

    class_cover = E56.build_class_orbit_cover(problem, cover)
    class_below = E56.decide_class_orbit_bundle(
        problem, cover, class_cover, 1, conflict_budget=0
    )
    class_at = E56.decide_class_orbit_bundle(
        problem, cover, class_cover, 2, conflict_budget=0
    )
    assert class_below["status"] == "UNSAT"
    assert class_at["status"] == "SAT"
    assert class_at["weight"] == 2
    assert class_at["verification"]["valid"] is True

    component_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        conflict_budget=0,
        formulation="component",
    )
    component_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        conflict_budget=0,
        formulation="component",
    )
    assert component_below["status"] == "UNSAT"
    assert component_at["status"] == "SAT"
    assert component_at["weight"] == 2

    generator_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        conflict_budget=0,
        formulation="generator",
    )
    generator_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        conflict_budget=0,
        formulation="generator",
    )
    assert generator_below["status"] == "UNSAT"
    assert generator_at["status"] == "SAT"
    assert generator_at["weight"] == 2
    assert generator_at["verification"]["valid"] is True
    assert component_at["verification"]["valid"] is True

    native_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        conflict_budget=0,
        solver_name="minicard",
        card_encoding="native",
    )
    native_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        conflict_budget=0,
        solver_name="minicard",
        card_encoding="native",
    )
    assert native_below["status"] == "UNSAT"
    assert native_at["status"] == "SAT"
    assert native_at["weight"] == 2
    assert native_at["verification"]["valid"] is True

    lex_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        conflict_budget=0,
        solver_name="minicard",
        card_encoding="native",
        formulation="generator",
        symmetry_break="lex",
    )
    lex_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        conflict_budget=0,
        solver_name="minicard",
        card_encoding="native",
        formulation="generator",
        symmetry_break="lex",
    )
    assert lex_below["status"] == "UNSAT"
    assert lex_at["status"] == "SAT"
    assert lex_at["weight"] == 2
    assert lex_at["verification"]["valid"] is True

    cpsat_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        backend="cpsat",
        time_limit_s=30.0,
    )

    cpsat_generator_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        backend="cpsat",
        formulation="generator",
        time_limit_s=30.0,
    )
    cpsat_generator_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        backend="cpsat",
        formulation="generator",
        time_limit_s=30.0,
    )
    assert cpsat_generator_below["status"] == "UNSAT"
    assert cpsat_generator_at["status"] == "SAT"
    assert cpsat_generator_at["weight"] == 2
    assert cpsat_generator_at["verification"]["valid"] is True

    milp_below = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        1,
        backend="milp",
        formulation="generator",
        time_limit_s=30.0,
    )
    milp_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        backend="milp",
        formulation="generator",
        time_limit_s=30.0,
    )
    assert milp_below["status"] == "UNSAT"
    assert milp_at["status"] == "SAT"
    assert milp_at["weight"] == 2
    assert milp_at["verification"]["valid"] is True
    cpsat_at = E56.decide_class_orbit_bundle(
        problem,
        cover,
        class_cover,
        2,
        backend="cpsat",
        time_limit_s=30.0,
    )
    assert cpsat_below["status"] == "UNSAT"
    assert cpsat_at["status"] == "SAT"
    assert cpsat_at["weight"] == 2
    assert cpsat_at["verification"]["valid"] is True
