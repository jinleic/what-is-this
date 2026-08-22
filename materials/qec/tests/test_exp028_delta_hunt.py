"""Focused contracts for EXP-028's search predicates and artifact routing."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from qec_research.codes.bicycle import poly_matrix


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp028", ROOT / "experiments" / "exp028_delta_hunt.py"
)
assert SPEC and SPEC.loader
EXP028 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP028)


def test_translation_keys_preserve_expected_orbits() -> None:
    ell, m = 6, 4
    A = [(0, 0), (1, 2), (2, 1)]
    B = [(0, 0), (1, 0), (2, 3)]
    C = [(2, 2), (3, 0)]
    D = [(1, 3)]
    shifted_A = [((a + 2) % ell, (b + 1) % m) for a, b in A]
    shifted_C = [((a + 2) % ell, (b + 1) % m) for a, b in C]
    shifted_B = [((a + 4) % ell, (b + 3) % m) for a, b in B]
    shifted_D = [((a + 4) % ell, (b + 3) % m) for a, b in D]

    assert EXP028.parent_translation_key(ell, m, A, B) == EXP028.parent_translation_key(
        ell, m, shifted_A, shifted_B
    )
    assert EXP028.pbb_translation_key(ell, m, A, B, C, D) == EXP028.pbb_translation_key(
        ell, m, shifted_A, shifted_B, shifted_C, shifted_D
    )


def test_selected_perturbations_commute_and_are_mixed() -> None:
    _, catalogue_orbits = EXP028.catalogue_summary()
    parents, _ = EXP028.sample_parents()
    parent = next(row for row in parents if row["lattice"] == "6x4" and row["k_parent"] == 8)
    seed = int(
        np.random.SeedSequence(
            [
                EXP028.SEED,
                parent["ell"],
                parent["m"],
                parent["sample_index_within_lattice"],
            ]
        ).generate_state(1)[0]
    )
    candidates, accounting = EXP028.enumerate_perturbations(
        parent, catalogue_orbits, seed
    )
    assert 0 < len(candidates) <= EXP028.PERTURBATIONS_PER_PARENT_CAP
    assert accounting["commutation_valid_supports"] >= len(candidates)

    A = poly_matrix(parent["ell"], parent["m"], parent["A"])
    B = poly_matrix(parent["ell"], parent["m"], parent["B"])
    for candidate in candidates:
        C = poly_matrix(parent["ell"], parent["m"], candidate["C"])
        D = poly_matrix(parent["ell"], parent["m"], candidate["D"])
        assert not EXP028.commutation_defect(A, B, C, D).any()
        code = EXP028.build_pbb(
            EXP028.PBBSpec(
                ell=parent["ell"],
                m=parent["m"],
                A=parent["A"],
                B=parent["B"],
                C=candidate["C"],
                D=candidate["D"],
            )
        )
        mixed = code.H[:, : code.n].any(axis=1) & code.H[:, code.n :].any(axis=1)
        assert mixed.any()


def test_partial_run_routes_away_from_canonical() -> None:
    route, path = EXP028.routed_output_path(
        EXP028.OUTPUT,
        clean=True,
        full_declared_scope=False,
        stop_reason="test_partial",
        timestamp=20260812,
    )
    assert route == "partial_runs"
    assert path != EXP028.OUTPUT
    assert path.parent == ROOT / "results" / "partial_runs"

    route, path = EXP028.routed_output_path(
        EXP028.OUTPUT,
        clean=False,
        full_declared_scope=True,
        stop_reason="test_failure",
        timestamp=20260812,
    )
    assert route == "quarantine"
    assert path != EXP028.OUTPUT
    assert path.parent == ROOT / "results" / "quarantine"


def test_coverage_statement_does_not_multiply_or_imply_empty_lattices() -> None:
    per_lattice = {
        "6x6": {"perturbations_solver_tested": 247},
        "6x4": {"perturbations_solver_tested": 0},
        "4x6": {"perturbations_solver_tested": 0},
    }
    text = EXP028.format_coverage_statement(
        verified_hit=False,
        parents_reached=3,
        perturbations_tested=247,
        per_lattice=per_lattice,
    )
    assert "247 perturbations across 3 parents" in text
    assert "solver-tested lattice(s) (6,6)" in text
    assert "zero solver-tested perturbations on (6,4), (4,6)" in text
    assert "3 parents x 247" not in text
