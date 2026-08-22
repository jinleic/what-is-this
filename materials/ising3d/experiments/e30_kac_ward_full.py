"""Exact diagnostics for the unresolved full scalar Kac--Ward family.

This computation deliberately does not turn a resource limit into a no-go claim.
It writes ``results/kac_ward/full_family.json`` and prints PASS when every
recorded exact invariant and numerical falsification reproduces.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import random
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    exact_kac_ward_polynomial,
    full_weight_finite_system,
    polynomial_square,
)
from ising.lattices import square

SQUARE_CASES = ((3, 3), (3, 4), (4, 4))
CONSTRUCTION_SHAPES = ((3, 3, 2), (2, 2, 3))
FIT_SHAPES = (
    (2, 2, 2), (3, 2, 2), (2, 3, 2), (2, 2, 3), (4, 2, 2),
    (2, 4, 2), (2, 2, 4), (3, 3, 2), (3, 2, 3), (2, 3, 3),
)
ORDERS = (4, 6, 8)
PRIMES = (101, 1009, 10007)
SEED = 314159
MODULAR_WITNESSES = {
    3: (2, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 0, 0, 1, 1, 2, 2, 2, 1, 0, 1),
    5: (2, 4, 0, 2, 4, 3, 2, 2, 4, 0, 4, 3, 0, 4, 4, 3, 2, 2, 2, 0, 0, 4, 4, 3, 3),
    7: (2, 3, 3, 4, 3, 1, 3, 5, 3, 5, 3, 1, 3, 5, 1, 5, 1, 0, 3, 0, 2, 6, 2, 6, 4),
    11: (7, 10, 8, 5, 1, 7, 7, 8, 0, 3, 3, 5, 9, 10, 6, 1, 2, 5, 2, 6, 4, 3, 0, 2, 4),
}


def _pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def _rank_mod(rows: list[list[int]], prime: int) -> int:
    matrix = [[int(value) % prime for value in row] for row in rows]
    rank = 0
    for column in range(len(matrix[0]) if matrix else 0):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]), None
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[rank], strict=True)
            ]
        rank += 1
    return rank


def two_dimensional_controls() -> list[dict]:
    records = []
    for shape in SQUARE_CASES:
        lattice = square(*shape, periodic=False)
        target = polynomial_square(even_subgraph_polynomial(lattice))
        observed, certificate = exact_kac_ward_polynomial(lattice)
        discrepancy = max(
            (abs(left - right) for left, right in zip(observed, target, strict=True)),
            default=0,
        )
        records.append(
            {
                "shape": list(shape),
                "boundary": "free",
                "directed_edges": certificate.directed_edge_count,
                "maximum_exact_coefficient_discrepancy": discrepancy,
                "max_nonrational_component": certificate.max_nonrational_component,
                "crt_modulus": str(certificate.modulus),
                "coefficient_bound": str(certificate.coefficient_bound),
            }
        )
    return records


def minimal_system() -> tuple[object, dict]:
    system = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=True
    )
    canonical = "\n".join(
        f"{shape}:{order}:{sp.srepr(equation)}"
        for (shape, order), equation in zip(
            system.labels, system.equations, strict=True
        )
    )
    record = {
        "field": "Q",
        "normalization": "Tr_shape(U^k)+2*k*[v^k]log(P_shape(v))=0",
        "construction_shapes": [list(shape) for shape in CONSTRUCTION_SHAPES],
        "orders": list(ORDERS),
        "variable_count": len(system.variables),
        "variables": [str(variable) for variable in system.variables],
        "equation_count": len(system.equations),
        "labels": [
            {"shape": list(shape), "order": order} for shape, order in system.labels
        ],
        "equations": [str(sp.expand(equation)) for equation in system.equations],
        "equation_term_counts": [
            len(sp.expand(equation).as_ordered_terms()) for equation in system.equations
        ],
        "sha256_sympy_srepr": hashlib.sha256(canonical.encode()).hexdigest(),
    }
    return system, record


def modular_solution_witnesses(system) -> list[dict]:
    """Verify stored exact F_p points on the minimal system and a holdout."""
    holdout = full_weight_finite_system(((3, 3, 3),), ORDERS, gauge_fix=True)
    records = []
    for prime, values in MODULAR_WITNESSES.items():
        substitutions = dict(zip(system.variables, values, strict=True))
        construction = [
            int(equation.subs(substitutions)) % prime for equation in system.equations
        ]
        prediction = [
            int(equation.subs(substitutions)) % prime for equation in holdout.equations
        ]
        records.append(
            {
                "prime": prime,
                "values_in_variable_order": list(values),
                "construction_residues": construction,
                "independent_3x3x3_residues_k4_k6_k8": prediction,
                "construction_solution_verified": not any(construction),
                "independent_prediction_passed": not any(prediction),
                "decision": "nonempty_reduction_mod_prime",
                "scope": (
                    "An exact F_p solution proves this reduced finite-field variety is nonempty. "
                    "It neither lifts to characteristic zero automatically nor decides the Q-variety."
                ),
            }
        )
    return records


def zero_pattern_strata() -> list[dict]:
    ungauged = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
    )
    symbol_map = dict(ungauged.all_symbols)
    remaining = [
        symbol for pair, symbol in ungauged.all_symbols if pair not in DIRECTION_GAUGE_TREE
    ]
    records = []
    for bits in itertools.product((0, 1), repeat=len(DIRECTION_GAUGE_TREE)):
        substitutions = {
            symbol_map[pair]: value
            for pair, value in zip(DIRECTION_GAUGE_TREE, bits, strict=True)
        }
        equations = [sp.expand(equation.subs(substitutions)) for equation in ungauged.equations]
        active = sorted(
            set().union(*(equation.free_symbols for equation in equations)), key=str
        )
        globally_linear = [
            variable
            for variable in remaining
            if max(sp.degree(equation, variable) for equation in equations) <= 1
            and any(sp.degree(equation, variable) == 1 for equation in equations)
        ]
        jacobian = sp.Matrix(equations).jacobian(active)
        ranks = []
        for prime in PRIMES:
            rng = random.Random(
                SEED + prime + sum((index + 1) * bit for index, bit in enumerate(bits))
            )
            point = {variable: rng.randrange(prime) for variable in active}
            rows = [
                [int(value.subs(point)) for value in row] for row in jacobian.tolist()
            ]
            ranks.append(_rank_mod(rows, prime))
        records.append(
            {
                "branch": "".join(map(str, bits)),
                "bit_convention": "0=fixed tree weight vanishes; 1=normalized to one",
                "zero_pattern": [
                    _pair_name(pair)
                    for pair, bit in zip(DIRECTION_GAUGE_TREE, bits, strict=True)
                    if not bit
                ],
                "normalized_nonzero": [
                    _pair_name(pair)
                    for pair, bit in zip(DIRECTION_GAUGE_TREE, bits, strict=True)
                    if bit
                ],
                "ambient_dimension_before_equations": 25,
                "active_variable_count_in_six_equations": len(active),
                "globally_linear_variable_count": len(globally_linear),
                "random_point_jacobian_ranks": dict(zip(map(str, PRIMES), ranks, strict=True)),
                "decision": "unresolved",
                "method": (
                    "exact substitution and degree scan; exact modular Jacobian rank at one "
                    "deterministic random point per prime"
                ),
                "scope": (
                    "Jacobian rank is local differential information only and certifies neither "
                    "existence nor emptiness; no branch is discarded by this diagnostic"
                ),
            }
        )
    return records


def numerical_survivor_and_holdout() -> dict:
    construction = full_weight_finite_system(FIT_SHAPES, ORDERS, gauge_fix=True)
    evaluator = sp.lambdify(
        construction.variables, construction.equations, "numpy", cse=True
    )
    scales = np.array(
        [1000 if order == 4 else 10000 if order == 6 else 100000
         for _, order in construction.labels]
    )

    def residual(real_vector: np.ndarray) -> np.ndarray:
        values = real_vector[:25] + 1j * real_vector[25:]
        raw = np.asarray(evaluator(*values), dtype=complex) / scales
        return np.r_[raw.real, raw.imag]

    initial = np.random.default_rng(12345).normal(size=50) * 0.7
    fit = least_squares(
        residual, initial, max_nfev=5000,
        ftol=1e-13, xtol=1e-13, gtol=1e-13,
    )
    values = fit.x[:25] + 1j * fit.x[25:]
    raw = np.asarray(evaluator(*values), dtype=complex)

    holdout = full_weight_finite_system(((3, 3, 3),), ORDERS, gauge_fix=True)
    holdout_evaluator = sp.lambdify(
        holdout.variables, holdout.equations, "numpy", cse=True
    )
    predicted = np.asarray(holdout_evaluator(*values), dtype=complex)
    rows = [
        {
            "order": order,
            "trace_equation_residual_real": format(value.real, ".17g"),
            "trace_equation_residual_imag": format(value.imag, ".17g"),
            "absolute_residual": format(abs(value), ".17g"),
        }
        for (_, order), value in zip(holdout.labels, predicted, strict=True)
    ]
    return {
        "status": "numerically_falsified_on_independent_prediction",
        "claim_class": "NUMERICAL; not an algebraic solution or no-go certificate",
        "fit_shapes": [list(shape) for shape in FIT_SHAPES],
        "fit_orders": list(ORDERS),
        "seed": 12345,
        "least_squares_function_evaluations": fit.nfev,
        "maximum_scaled_construction_residual": format(
            max(abs(residual(fit.x))), ".17g"
        ),
        "maximum_raw_construction_residual": format(max(abs(raw)), ".17g"),
        "weights": {
            str(variable): {
                "real": format(value.real, ".17g"),
                "imag": format(value.imag, ".17g"),
            }
            for variable, value in zip(construction.variables, values, strict=True)
        },
        "independent_holdout_shape": [3, 3, 3],
        "independent_holdout": rows,
        "first_failed_prediction_order": 8,
        "detail": (
            "The fitted point satisfies the construction equations numerically but misses the "
            "exact even-subgraph-derived 3x3x3 k=8 trace target. It falsifies this point only."
        ),
    }


def main() -> None:
    controls = two_dimensional_controls()
    system, system_record = minimal_system()
    modular_witnesses = modular_solution_witnesses(system)
    branches = zero_pattern_strata()
    numerical = numerical_survivor_and_holdout()

    checks = [
        {
            "name": "2d_exact_control",
            "passed": all(
                row["maximum_exact_coefficient_discrepancy"] == 0
                and row["max_nonrational_component"] == 0
                for row in controls
            ),
            "detail": "3x3, 3x4, 4x4 free square lattices",
        },
        {
            "name": "minimal_two_box_system",
            "passed": system_record["variable_count"] == 25
            and system_record["equation_count"] == 6,
            "detail": "two different boxes, k=4,6,8, exact Q polynomials",
        },
        {
            "name": "exact_modular_solution_witnesses",
            "passed": all(
                row["construction_solution_verified"] for row in modular_witnesses
            ),
            "detail": "explicit solutions of the six construction equations over F_3,F_5,F_7,F_11",
        },
        {
            "name": "all_fixed_tree_zero_patterns_enumerated",
            "passed": len(branches) == 32
            and {row["branch"] for row in branches}
            == {"".join(bits) for bits in itertools.product("01", repeat=5)},
            "detail": "all 2^5 zero/nonzero patterns of the previous gauge tree",
        },
        {
            "name": "no_modular_rank_overclaim",
            "passed": all(row["decision"] == "unresolved" for row in branches),
            "detail": "random-point Jacobian ranks are explicitly non-decisive",
        },
        {
            "name": "independent_prediction_exercised",
            "passed": numerical["first_failed_prediction_order"] == 8
            and float(numerical["independent_holdout"][2]["absolute_residual"]) > 100,
            "detail": "numerical construction survivor fails exact 3x3x3 k=8 target",
        },
    ]
    payload = {
        "provenance": {
            "script": "experiments/e30_kac_ward_full.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": (
                "exact Python integers/sympy QQ; modular ranks over F_101,F_1009,F_10007; "
                "IEEE float64 only for explicitly numerical least-squares survivor"
            ),
            "random_seeds": {"jacobian": SEED, "least_squares": 12345},
        },
        "data": {
            "headline": "UNRESOLVED full translation-invariant scalar family at k=8",
            "two_dimensional_controls": controls,
            "minimal_unresolved_system": system_record,
            "exact_modular_solution_witnesses": modular_witnesses,
            "zero_pattern_branches": branches,
            "linear_elimination": {
                "status": "partial_nonzero_pivot_chart_only",
                "independent_equations_before": 20,
                "eliminated_variables": [
                    "u_mx_my", "u_my_mz", "u_mx_mx", "u_mx_mz", "u_my_my", "u_mz_mz"
                ],
                "remaining_equations": 9,
                "remaining_variables": 19,
                "maximum_terms_in_remaining_equation": 311,
                "scope": (
                    "Each division adds a nonzero pivot guard and therefore creates a zero "
                    "branch. Continuing the exact substitutions caused expression growth; this "
                    "does not decide any guarded or zero branch."
                ),
            },
            "modular_elimination_resource_wall": {
                "attempts": [
                    {
                        "system": "30 equations, 25 variables, ten free boxes, k=4,6,8",
                        "field": "F_101",
                        "method": "sympy groebner grevlex",
                        "limit_seconds": 3600,
                        "status": "timeout",
                    },
                    {
                        "system": "6 equations, 25 variables, 3x3x2 and 2x2x3, k=4,6,8",
                        "field": "F_101",
                        "method": "sympy groebner grevlex",
                        "limit_seconds": 3600,
                        "status": "timeout",
                    },
                ],
                "conclusion": "timeouts certify nothing about consistency or emptiness",
            },
            "numerical_surviving_point": numerical,
            "scope": {
                "proved_or_exactly_reproduced": (
                    "the listed 2D determinant identities; the exact construction/branch "
                    "polynomials and modular Jacobian ranks; and explicit solutions of the "
                    "minimal construction system over F_3,F_5,F_7,F_11"
                ),
                "not_proved": (
                    "existence or nonexistence of a full translation-invariant scalar solution; "
                    "coverage of arbitrary zero patterns among all 30 weights; any k>8 identity"
                ),
                "unchanged_prior_fact": (
                    "the strict cubic-covariant two-weight slice is exactly impossible at k=8"
                ),
            },
        },
        "checks": checks,
    }
    output = Path("results/kac_ward/full_family.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for check in checks:
        print(f"{check['name']}: {'PASS' if check['passed'] else 'FAIL'} -- {check['detail']}")
    if not all(check["passed"] for check in checks):
        print("FAIL")
        raise SystemExit(1)
    print(f"wrote {output}")
    print("PASS")


if __name__ == "__main__":
    main()
