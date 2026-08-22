"""Standalone regression checks for the full Kac--Ward branch artifact."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    exact_kac_ward_polynomial,
    full_weight_finite_system,
    polynomial_square,
)
from ising.lattices import square

RESULT = Path("results/kac_ward/full_family.json")


def check_2d_control() -> None:
    for shape in ((3, 3), (3, 4), (4, 4)):
        lattice = square(*shape, periodic=False)
        observed, certificate = exact_kac_ward_polynomial(lattice)
        target = polynomial_square(even_subgraph_polynomial(lattice))
        assert tuple(target) == observed, shape
        assert certificate.max_nonrational_component == 0, shape
        print(f"2D {shape[0]}x{shape[1]} free: exact coefficient discrepancy 0")


def check_exact_minimal_system() -> None:
    system = full_weight_finite_system(
        ((3, 3, 2), (2, 2, 3)), (4, 6, 8), gauge_fix=True
    )
    assert len(system.variables) == 25
    assert len(system.equations) == 6
    assert system.labels == (
        ((3, 3, 2), 4), ((3, 3, 2), 6), ((3, 3, 2), 8),
        ((2, 2, 3), 4), ((2, 2, 3), 6), ((2, 2, 3), 8),
    )
    assert [len(equation.as_ordered_terms()) for equation in system.equations] == [
        7, 41, 350, 7, 37, 179
    ]


def check_result_scope() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["data"]["headline"].startswith("UNRESOLVED")
    witnesses = payload["data"]["exact_modular_solution_witnesses"]
    assert {row["prime"] for row in witnesses} == {3, 5, 7, 11}
    assert all(row["construction_solution_verified"] for row in witnesses)
    assert all(not any(row["construction_residues"]) for row in witnesses)
    assert any(not row["independent_prediction_passed"] for row in witnesses)
    branches = payload["data"]["zero_pattern_branches"]
    assert len(branches) == 32
    assert {row["branch"] for row in branches} == {
        "".join(bits) for bits in itertools.product("01", repeat=5)
    }
    assert all(row["decision"] == "unresolved" for row in branches)
    assert all(row["ambient_dimension_before_equations"] == 25 for row in branches)
    for row in branches:
        ranks = row["random_point_jacobian_ranks"]
        assert set(ranks) == {"101", "1009", "10007"}
        assert len(set(ranks.values())) == 1
    attempts = payload["data"]["modular_elimination_resource_wall"]["attempts"]
    assert len(attempts) == 2
    assert all(attempt["status"] == "timeout" for attempt in attempts)
    numerical = payload["data"]["numerical_surviving_point"]
    assert numerical["claim_class"].startswith("NUMERICAL")
    assert numerical["first_failed_prediction_order"] == 8
    assert float(numerical["independent_holdout"][2]["absolute_residual"]) > 100
    assert all(check["passed"] for check in payload["checks"])


def main() -> None:
    check_2d_control()
    check_exact_minimal_system()
    check_result_scope()
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise
