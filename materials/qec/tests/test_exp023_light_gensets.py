"""Focused soundness tests for EXP-023's physical-space closure certificate."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "exp023_light_gensets", ROOT / "experiments" / "exp023_light_gensets.py"
)
assert SPEC is not None and SPEC.loader is not None
EXP023 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXP023)


def test_physical_functionals_are_complete_and_pulled_back() -> None:
    # Independent 3-row physical basis in a 5-bit ambient space.
    h_basis = np.array(
        [[1, 0, 1, 0, 0], [0, 1, 1, 0, 0], [0, 0, 0, 1, 1]],
        dtype=np.uint8,
    )
    found_v = np.array([h_basis[0] ^ h_basis[1]], dtype=np.uint8)

    functionals, physical_h, checks = EXP023.physical_annihilator_functionals(
        found_v, h_basis
    )

    assert functionals.shape == (2, 3)
    assert EXP023.numpy_gf2_rank(functionals) == 2
    assert not EXP023.numpy_gf2_matmul(found_v, physical_h.T).any()
    assert np.array_equal(
        functionals, EXP023.numpy_gf2_matmul(physical_h, h_basis.T)
    )
    assert checks["physical_functionals_verified"] is True


def test_query_excludes_exactly_the_found_physical_span() -> None:
    # Two qubits, symplectic coordinates (x0,x1|z0,z1).  The three basis
    # elements all have weight 1, so V_1 is the full 3-dimensional rowspace.
    h_basis = np.array(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]], dtype=np.uint8
    )
    found_v = h_basis[:1]
    functionals, _, _ = EXP023.physical_annihilator_functionals(found_v, h_basis)
    model, _ = EXP023.build_query_model(h_basis, 1, functionals)
    solver = cp_model.CpSolver()
    assert solver.solve(model) in (cp_model.OPTIMAL, cp_model.FEASIBLE)

    full_functionals, _, checks = EXP023.physical_annihilator_functionals(
        h_basis, h_basis
    )
    assert full_functionals.shape == (0, 3)
    assert checks["physical_functionals_verified"] is True


def test_dependent_published_rows_map_through_bijective_basis() -> None:
    h_basis = np.array(
        [[1, 0, 1, 0], [0, 1, 0, 1], [1, 1, 0, 0]], dtype=np.uint8
    )
    original_h = np.vstack([h_basis, h_basis[0] ^ h_basis[2]])
    coordinates = EXP023.coefficient_map(h_basis, original_h)

    assert np.array_equal(
        EXP023.numpy_gf2_matmul(coordinates, h_basis), original_h
    )
    assert EXP023.numpy_gf2_rank(coordinates) == 3
