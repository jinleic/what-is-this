"""Clean-room exact verifier for the wave-10 branch-11111 cubic study.

This file deliberately imports no producer experiment.  It reconstructs the
finite Kac--Ward system from stable ``ising`` APIs, independently derives the
Macaulay and thin-torus records, and verifies the F_5 Hensel certificate.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    full_weight_finite_system,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/branch11111_cubic.json"
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(
        set(itertools.permutations((2, 2, 1)))
        | set(itertools.permutations((3, 2, 1)))
    )
)
FULL_SHAPES = ((3, 3, 2), (2, 2, 3), (2, 2, 2), (3, 2, 2), (2, 3, 2))
CONSTRUCTION_SHAPES = PRIMARY_SHAPES + FULL_SHAPES
HOLDOUT_SHAPE = (3, 3, 3)

EXPECTED_BASIS_INDICES = [
    0,
    2,
    4,
    5,
    7,
    9,
    11,
    13,
    14,
    15,
    17,
    19,
    20,
    22,
    25,
    28,
    29,
    32,
    35,
    38,
    41,
]
EXPECTED_DENSE = {
    "multiplier_count": 3276,
    "column_count": 68796,
    "basis_term_count": 1001,
    "nonzero_incidence_count": 3279276,
    "distinct_monomial_row_count": 1612152,
}
EXPECTED_SYZYGY = {
    "base_relation_count": 6,
    "lift_multiplier_count": 26,
    "lifted_relation_count": 156,
    "lifted_relation_rank": 156,
    "relation_sha256": "a101ddfc652e6b1e811ef90889528e8c44a48c2089f2d460e3478f05de620011",
    "pivot_sha256": "d7d8bc25087a8e35f9be150d6914d5ddbb469b1928f2f125dcbe9acdb195a196",
    "pruned_column_count": 68640,
}
EXPECTED_SPARSE = {
    "all_degree_le_2_plus_slice6_degree_3": {
        "extra_degree_three_multiplier_count": 56,
        "multiplier_count": 407,
        "column_count": 8547,
        "column_rank": 8541,
        "dependent_column_count": 6,
        "pivot_monomial_sha256": "23668ce7e67d1f7b5d78e2a9b67c6308b1b49d9095336b29cc7fec9a9347a4c8",
    },
    "all_degree_le_2_plus_slice8_degree_3": {
        "extra_degree_three_multiplier_count": 120,
        "multiplier_count": 471,
        "column_count": 9891,
        "column_rank": 9885,
        "dependent_column_count": 6,
        "pivot_monomial_sha256": "a73059dc9f6d25fc9a2e17707a9a80d3980559aaafbee3c7e6c87fd4e2d5df0e",
    },
}
SELECTED_NUMERATOR_INDICES = [28, 29, 32, 35, 38, 41]
DIAGONAL_NAMES = ["u_px_px", "u_py_py", "u_pz_pz"]
NONDIAGONAL_NAMES = [
    "u_py_px",
    "u_py_pz",
    "u_py_mz",
    "u_my_px",
    "u_my_mx",
    "u_pz_px",
    "u_pz_mx",
    "u_pz_py",
    "u_pz_my",
    "u_mz_px",
    "u_mz_mx",
    "u_mz_py",
    "u_mz_my",
]
SLICE6_NAMES = NONDIAGONAL_NAMES[:6]
SLICE8_NAMES = SLICE6_NAMES + ["u_pz_mx", "u_mz_mx"]
HENSEL_POINT_MOD_5 = [1, 2, 3, 3, 1, 1, 1, 1, 1, 1, 1, 4, 4]
HENSEL_LIFT_MOD_625 = [226, 197, 98, 63, 161, 216]

Monomial = tuple[int, ...]
SparseInteger = dict[Monomial, int]
SparseRational = dict[Monomial, Fraction]


def monomial_key(monomial: Monomial) -> tuple[int, Monomial]:
    return (sum(monomial), monomial)


def sparse_poly(expression: sp.Expr, variables: list[sp.Symbol]) -> SparseInteger:
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return {
        tuple(exponents): int(coefficient)
        for exponents, coefficient in polynomial.terms()
        if coefficient
    }


def sparse_hash(columns: list[SparseInteger]) -> str:
    digest = hashlib.sha256()
    for index, column in enumerate(columns):
        digest.update(f"C{index}:".encode())
        for monomial, coefficient in sorted(column.items()):
            digest.update((f"{coefficient}@{','.join(map(str, monomial))};").encode())
    return digest.hexdigest()


def branch_11111_equations(system) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    fixed = set(DIRECTION_GAUGE_TREE)
    substitutions = {symbol: 1 for pair, symbol in system.all_symbols if pair in fixed}
    variables = [symbol for pair, symbol in system.all_symbols if pair not in fixed]
    equations = [sp.expand(expression.subs(substitutions)) for expression in system.equations]
    return variables, equations


@lru_cache(maxsize=1)
def construction_data() -> tuple[list[sp.Symbol], list[sp.Expr], list[dict], object]:
    system = full_weight_finite_system(
        CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False
    )
    variables, equations = branch_11111_equations(system)
    labels = [{"shape": list(shape), "order": order} for shape, order in system.labels]
    return variables, equations, labels, system


def independent_indices(expressions: list[sp.Expr], variables: list[sp.Symbol]) -> list[int]:
    basis: dict[Monomial, SparseRational] = {}
    selected = []
    for index, expression in enumerate(expressions):
        vector: SparseRational = {
            monomial: Fraction(coefficient)
            for monomial, coefficient in sparse_poly(expression, variables).items()
        }
        while vector:
            pivot = max(vector, key=monomial_key)
            coefficient = vector[pivot]
            existing = basis.get(pivot)
            if existing is None:
                inverse = Fraction(1, 1) / coefficient
                basis[pivot] = {
                    monomial: value * inverse
                    for monomial, value in vector.items()
                    if value
                }
                selected.append(index)
                break
            for monomial, value in existing.items():
                reduced = vector.get(monomial, Fraction(0)) - coefficient * value
                if reduced:
                    vector[monomial] = reduced
                else:
                    vector.pop(monomial, None)
    return selected


def multiplier_exponents(
    variable_count: int,
    degree_bound: int,
    *,
    exact_degree: int | None = None,
    allowed_positions: list[int] | None = None,
) -> list[Monomial]:
    positions = list(range(variable_count)) if allowed_positions is None else allowed_positions
    output = []
    for degree in range(degree_bound + 1):
        if exact_degree is not None and degree != exact_degree:
            continue
        for factors in itertools.combinations_with_replacement(positions, degree):
            exponents = [0] * variable_count
            for position in factors:
                exponents[position] += 1
            output.append(tuple(exponents))
    return output


def macaulay_columns(
    polynomials: list[SparseInteger], multipliers: list[Monomial]
) -> list[SparseInteger]:
    return [
        {
            tuple(left + right for left, right in zip(monomial, multiplier, strict=True)): coefficient
            for monomial, coefficient in polynomial.items()
        }
        for polynomial in polynomials
        for multiplier in multipliers
    ]


def fraction_reduction(
    columns: list[SparseInteger], variable_count: int, *, capture_relations: bool
) -> dict:
    basis: dict[Monomial, SparseRational] = {}
    basis_relations: dict[Monomial, dict[int, Fraction]] = {}
    dependent = []
    relations = []
    pivot_order = []
    for index, source in enumerate(columns):
        vector: SparseRational = {
            monomial: Fraction(coefficient) for monomial, coefficient in source.items()
        }
        relation = {index: Fraction(1)}
        while vector:
            pivot = max(vector, key=monomial_key)
            coefficient = vector[pivot]
            existing = basis.get(pivot)
            if existing is None:
                inverse = Fraction(1, 1) / coefficient
                basis[pivot] = {
                    monomial: value * inverse
                    for monomial, value in vector.items()
                    if value
                }
                if capture_relations:
                    basis_relations[pivot] = {
                        column: value * inverse
                        for column, value in relation.items()
                        if value
                    }
                pivot_order.append(pivot)
                break
            for monomial, value in existing.items():
                reduced = vector.get(monomial, Fraction(0)) - coefficient * value
                if reduced:
                    vector[monomial] = reduced
                else:
                    vector.pop(monomial, None)
            if capture_relations:
                for column, value in basis_relations[pivot].items():
                    reduced = relation.get(column, Fraction(0)) - coefficient * value
                    if reduced:
                        relation[column] = reduced
                    else:
                        relation.pop(column, None)
        if not vector:
            dependent.append(index)
            if capture_relations:
                relations.append({column: value for column, value in relation.items() if value})

    target: SparseRational = {(0,) * variable_count: Fraction(1)}
    while target:
        pivot = max(target, key=monomial_key)
        coefficient = target[pivot]
        existing = basis.get(pivot)
        if existing is None:
            break
        for monomial, value in existing.items():
            reduced = target.get(monomial, Fraction(0)) - coefficient * value
            if reduced:
                target[monomial] = reduced
            else:
                target.pop(monomial, None)
    return {
        "rank": len(basis),
        "dependent": dependent,
        "relations": relations,
        "target": target,
        "pivot_monomial_sha256": hashlib.sha256(
            "".join(",".join(map(str, monomial)) + ";" for monomial in pivot_order).encode()
        ).hexdigest(),
    }


def sparse_vector_rank(vectors: list[dict[int, Fraction]]) -> tuple[int, list[int]]:
    basis: dict[int, dict[int, Fraction]] = {}
    pivots = []
    for source in vectors:
        vector = dict(source)
        while vector:
            pivot = max(vector)
            coefficient = vector[pivot]
            existing = basis.get(pivot)
            if existing is None:
                inverse = Fraction(1, 1) / coefficient
                basis[pivot] = {key: value * inverse for key, value in vector.items() if value}
                pivots.append(pivot)
                break
            for key, value in existing.items():
                reduced = vector.get(key, Fraction(0)) - coefficient * value
                if reduced:
                    vector[key] = reduced
                else:
                    vector.pop(key, None)
    return len(basis), pivots


def relation_hash(vectors: list[dict[int, Fraction]]) -> str:
    digest = hashlib.sha256()
    for index, vector in enumerate(vectors):
        for column, coefficient in sorted(vector.items()):
            digest.update(
                f"{index}:{column}:{coefficient.numerator}/{coefficient.denominator};".encode()
            )
    return digest.hexdigest()


def lifted_degree_two_syzygies(
    degree_two_relations: list[dict[int, Fraction]],
    degree_two_multipliers: list[Monomial],
    degree_three_multipliers: list[Monomial],
    basis_count: int,
) -> list[dict[int, Fraction]]:
    del basis_count  # The relation indices already encode the basis block.
    multiplier_position = {multiplier: index for index, multiplier in enumerate(degree_three_multipliers)}
    variable_count = len(degree_two_multipliers[0])
    degree_one = multiplier_exponents(variable_count, 1)
    output = []
    for relation in degree_two_relations:
        for factor in degree_one:
            lifted: dict[int, Fraction] = {}
            for column, coefficient in relation.items():
                basis_position, multiplier_index = divmod(column, len(degree_two_multipliers))
                multiplier = tuple(
                    left + right
                    for left, right in zip(degree_two_multipliers[multiplier_index], factor, strict=True)
                )
                new_column = (
                    basis_position * len(degree_three_multipliers)
                    + multiplier_position[multiplier]
                )
                reduced = lifted.get(new_column, Fraction(0)) + coefficient
                if reduced:
                    lifted[new_column] = reduced
                else:
                    lifted.pop(new_column, None)
            output.append(lifted)
    return output


def thin_torus_substitution(
    variables: list[sp.Symbol],
) -> tuple[list[sp.Symbol], dict[sp.Symbol, sp.Expr]]:
    symbol = {str(value): value for value in variables}
    substitutions = {
        symbol["u_mx_my"]: -1 / symbol["u_my_px"],
        symbol["u_mx_py"]: -1 / (symbol["u_my_mx"] * symbol["u_py_px"]),
        symbol["u_mx_mx"]: 1 / symbol["u_px_px"],
        symbol["u_mx_mz"]: -1 / (symbol["u_mz_px"] * symbol["u_pz_mx"]),
        symbol["u_mx_pz"]: -1 / (symbol["u_mz_mx"] * symbol["u_pz_px"]),
        symbol["u_mz_mz"]: 1 / symbol["u_pz_pz"],
        symbol["u_my_mz"]: -1
        / (symbol["u_mz_py"] * symbol["u_py_pz"] * symbol["u_pz_my"]),
        symbol["u_my_pz"]: -1
        / (symbol["u_mz_my"] * symbol["u_py_mz"] * symbol["u_pz_py"]),
        symbol["u_my_my"]: 1 / symbol["u_py_py"],
    }
    free = [variable for variable in variables if variable not in substitutions]
    assert len(free) == 16
    return free, substitutions


def primitive_torus_numerators(
    variables: list[sp.Symbol], equations: list[sp.Expr]
) -> tuple[list[sp.Symbol], dict[sp.Symbol, sp.Expr], list[sp.Expr], list[sp.Poly]]:
    free, substitutions = thin_torus_substitution(variables)
    diagonal = {free[index]: 1 for index in (0, 2, 11)}
    selected_numerators = []
    primitive = []
    non_diagonal = [variable for variable in free if str(variable) not in DIAGONAL_NAMES]
    for index in SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, _denominator = sp.fraction(reduced)
        selected_numerators.append(sp.expand(numerator))
        polynomial = sp.Poly(
            sp.expand(numerator.subs(diagonal)), *non_diagonal, domain=sp.ZZ
        )
        _content, primitive_polynomial = polynomial.primitive()
        primitive.append(primitive_polynomial)
    return free, substitutions, selected_numerators, primitive


def primitive_integer_vector(vector: sp.Matrix) -> list[int]:
    denominator = sp.ilcm(*[value.q for value in vector])
    values = [int(value * denominator) for value in vector]
    divisor = 0
    for value in values:
        divisor = gcd(divisor, abs(value))
    values = [value // divisor for value in values]
    first = next(value for value in values if value)
    if first < 0:
        values = [-value for value in values]
    return values


def eval_sparse_mod(
    polynomial: SparseInteger, values: list[int], modulus: int) -> int:
    total = 0
    for exponents, coefficient in polynomial.items():
        term = coefficient % modulus
        for exponent, value in zip(exponents, values, strict=True):
            if exponent:
                term = term * pow(value, exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def eval_expr_mod(
    expression: sp.Expr,
    variables: list[sp.Symbol],
    values: dict[sp.Symbol, int],
    modulus: int,
) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    numerator_value = eval_sparse_mod(sparse_poly(numerator, variables), [values[v] for v in variables], modulus)
    denominator_value = eval_sparse_mod(sparse_poly(denominator, variables), [values[v] for v in variables], modulus)
    assert denominator_value % 5
    return numerator_value * pow(denominator_value, -1, modulus) % modulus


def solve_linear_system_mod(matrix: list[list[int]], rhs: list[int], prime: int) -> list[int]:
    size = len(matrix)
    augmented = [
        [entry % prime for entry in row] + [value % prime]
        for row, value in zip(matrix, rhs, strict=True)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column] % prime)
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, prime)
        augmented[column] = [value * inverse % prime for value in augmented[column]]
        for row in range(size):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (left - factor * right) % prime
                for left, right in zip(augmented[row], augmented[column], strict=True)
            ]
    return [row[-1] for row in augmented]


def eval_poly_integer(polynomial: sp.Poly, values: list[int]) -> int:
    total = 0
    for exponents, coefficient in polynomial.terms():
        term = int(coefficient)
        for exponent, value in zip(exponents, values, strict=True):
            if exponent:
                term *= value**exponent
        total += term
    return total


def hensel_lift_slice(
    polynomials: list[sp.Poly], solve_variables: list[sp.Symbol], initial: list[int], prime: int, exponent: int
) -> tuple[list[int], int, list[list[int]]]:
    jacobian = [
        [sp.Poly(sp.diff(polynomial.as_expr(), variable), *solve_variables, domain=sp.ZZ) for variable in solve_variables]
        for polynomial in polynomials
    ]
    jacobian_mod_prime = [
        [eval_poly_integer(entry, initial) % prime for entry in row]
        for row in jacobian
    ]
    assert sp.Matrix(jacobian_mod_prime).det() % prime == 2
    values = list(initial)
    modulus = prime
    for _ in range(1, exponent):
        residuals = [eval_poly_integer(polynomial, values) for polynomial in polynomials]
        assert all(residual % modulus == 0 for residual in residuals)
        correction = solve_linear_system_mod(
            jacobian_mod_prime,
            [-(residual // modulus) for residual in residuals],
            prime,
        )
        values = [
            value + modulus * digit
            for value, digit in zip(values, correction, strict=True)
        ]
        modulus *= prime
    return values, modulus, jacobian_mod_prime


def check_envelope(payload: dict) -> None:
    assert payload["provenance"]["script"] == "experiments/e90_kw_cubic.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert payload["provenance"]["no_float_decisions"] is True
    assert payload["data"]["branch"] == "11111"
    assert payload["data"]["claim_tag"] == "[UNRESOLVED]"
    assert payload["data"]["status"] == "UNRESOLVED"


def check_degree_three_data(payload: dict) -> None:
    data = payload["data"]
    variables, equations, _labels, _system = construction_data()
    basis_indices = independent_indices(equations, variables)
    assert basis_indices == EXPECTED_BASIS_INDICES
    basis_polynomials = [sparse_poly(equations[index], variables) for index in basis_indices]
    degree_two = multiplier_exponents(len(variables), 2)
    degree_three = multiplier_exponents(len(variables), 3)
    dense = data["dense_degree_three"]
    assert dense["basis_constraint_indices"] == basis_indices
    assert len(degree_three) == EXPECTED_DENSE["multiplier_count"]
    assert dense["multiplier_count"] == len(degree_three)
    assert dense["column_count"] == len(basis_polynomials) * len(degree_three)
    assert dense["basis_term_count"] == sum(map(len, basis_polynomials)) == EXPECTED_DENSE["basis_term_count"]
    assert dense["nonzero_incidence_count"] == (
        len(degree_three) * sum(map(len, basis_polynomials))
    ) == EXPECTED_DENSE["nonzero_incidence_count"]
    rows = {
        tuple(left + right for left, right in zip(monomial, multiplier, strict=True))
        for polynomial in basis_polynomials
        for monomial in polynomial
        for multiplier in degree_three
    }
    assert len(rows) == dense["distinct_monomial_row_count"] == EXPECTED_DENSE["distinct_monomial_row_count"]

    degree_two_columns = macaulay_columns(basis_polynomials, degree_two)
    reduction = fraction_reduction(
        degree_two_columns, len(variables), capture_relations=True
    )
    assert reduction["rank"] == 7365
    assert len(reduction["relations"]) == EXPECTED_SYZYGY["base_relation_count"]
    lifted = lifted_degree_two_syzygies(
        reduction["relations"], degree_two, degree_three, len(basis_polynomials)
    )
    rank, pivots = sparse_vector_rank(lifted)
    pruning = dense["syzygy_pruning"]
    assert len(lifted) == EXPECTED_SYZYGY["lifted_relation_count"]
    assert rank == EXPECTED_SYZYGY["lifted_relation_rank"]
    assert relation_hash(lifted) == EXPECTED_SYZYGY["relation_sha256"]
    pivot_hash = hashlib.sha256(";".join(map(str, pivots)).encode()).hexdigest()
    assert pivot_hash == EXPECTED_SYZYGY["pivot_sha256"]
    for key, expected in EXPECTED_SYZYGY.items():
        assert pruning[key] == expected
    assert pruning["pruned_column_count"] == dense["column_count"] - rank

    variable_positions = {str(variable): position for position, variable in enumerate(variables)}
    sparse_by_name = {record["name"]: record for record in data["sparse_degree_three"]}
    assert set(sparse_by_name) == set(EXPECTED_SPARSE)
    for name, expected in EXPECTED_SPARSE.items():
        record = sparse_by_name[name]
        selected_names = record["degree_three_variable_names"]
        selected_positions = [variable_positions[value] for value in selected_names]
        extra = multiplier_exponents(
            len(variables), 3, exact_degree=3, allowed_positions=selected_positions
        )
        multipliers = degree_two + extra
        columns = macaulay_columns(basis_polynomials, multipliers)
        result = fraction_reduction(columns, len(variables), capture_relations=False)
        assert record["degree_two_multiplier_count"] == len(degree_two) == 351
        assert record["extra_degree_three_multiplier_count"] == len(extra) == expected["extra_degree_three_multiplier_count"]
        assert record["multiplier_count"] == len(multipliers) == expected["multiplier_count"]
        assert record["column_count"] == len(columns) == expected["column_count"]
        assert result["rank"] == record["column_rank"] == expected["column_rank"]
        assert len(result["dependent"]) == record["dependent_column_count"] == expected["dependent_column_count"]
        assert result["target"] == {(0,) * len(variables): Fraction(1)}
        assert record["constant_target_in_span"] is False
        assert result["pivot_monomial_sha256"] == record["pivot_monomial_sha256"] == expected["pivot_monomial_sha256"]


def check_torus_hensel_certificate(payload: dict) -> None:
    data = payload["data"]
    variables, equations, labels, _system = construction_data()
    free, substitutions, selected_numerators, primitive = primitive_torus_numerators(variables, equations)
    assert [str(value) for value in free] == data["torus_slice"]["free_variable_order"]
    assert [str(value) for value in free if str(value) not in DIAGONAL_NAMES] == NONDIAGONAL_NAMES
    assert [str(value) for value in free if str(value) in DIAGONAL_NAMES] == DIAGONAL_NAMES
    diagonal = {free[index]: 1 for index in (0, 2, 11)}
    non_diagonal = [value for value in free if str(value) not in DIAGONAL_NAMES]
    primitive_contents = [
        int(
            sp.Poly(
                sp.expand(numerator.subs(diagonal)),
                *non_diagonal,
                domain=sp.ZZ,
            ).primitive()[0]
        )
        for numerator in selected_numerators
    ]
    assert primitive_contents == [24, 16, 8, 16, 8, 8]
    assert [len(polynomial.terms()) for polynomial in primitive] == [33, 277, 155, 7, 155, 155]

    full_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in FULL_SHAPES
    ]
    nonzero_numerators = []
    for index in full_indices:
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, _denominator = sp.fraction(reduced)
        if numerator:
            nonzero_numerators.append(sp.expand(numerator))
    assert len(nonzero_numerators) == 10
    assert len(independent_indices(nonzero_numerators, free)) == 6
    assert len(independent_indices(selected_numerators, free)) == 6
    primary_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in PRIMARY_SHAPES
    ]
    assert all(sp.cancel(equations[index].subs(substitutions)) == 0 for index in primary_indices)

    support_rows = []
    for numerator in selected_numerators:
        terms = sp.Poly(numerator, *free, domain=sp.ZZ).terms()
        base = terms[0][0]
        support_rows.extend(
            [exponent - origin for exponent, origin in zip(monomial, base, strict=True)]
            for monomial, _coefficient in terms[1:]
        )
    support_matrix = sp.Matrix(support_rows)
    assert support_matrix.rank() == data["torus_slice"]["support_difference_rank"] == 13
    kernel = [primitive_integer_vector(vector) for vector in support_matrix.nullspace()]
    assert kernel == data["torus_slice"]["torus_action_kernel_vectors"]
    diagonal_positions = [0, 2, 11]
    diagonal_minor = sp.Matrix(kernel)[:, diagonal_positions].det()
    assert abs(int(diagonal_minor)) == data["torus_slice"]["diagonal_normalization_minor_abs"] == 1

    section = data["hensel_five_adic_section"]
    assert section["prime"] == 5
    assert section["modulus"] == 625
    assert section["numerator_constraint_indices"] == SELECTED_NUMERATOR_INDICES
    assert section["solve_variable_names"] == SLICE6_NAMES
    assert section["held_variable_names"] == NONDIAGONAL_NAMES[6:]
    assert section["initial_non_diagonal_values_mod_5"] == HENSEL_POINT_MOD_5
    assert section["primitive_contents"] == [24, 16, 8, 16, 8, 8]

    non_diagonal = [value for value in free if str(value) not in DIAGONAL_NAMES]
    solve_variables = non_diagonal[:6]
    held = dict(zip(non_diagonal[6:], HENSEL_POINT_MOD_5[6:], strict=True))
    slice_polynomials = [
        sp.Poly(polynomial.as_expr().subs(held), *solve_variables, domain=sp.ZZ)
        for polynomial in primitive
    ]
    initial = HENSEL_POINT_MOD_5[:6]
    assert all(eval_poly_integer(polynomial, initial) % 5 == 0 for polynomial in slice_polynomials)
    lift, modulus, jacobian = hensel_lift_slice(
        slice_polynomials, solve_variables, initial, 5, 4
    )
    assert modulus == 625
    assert lift == HENSEL_LIFT_MOD_625 == section["lifted_solve_values_mod_625"]
    assert jacobian == section["jacobian_mod_5"]
    assert int(sp.Matrix(jacobian).det()) % 5 == section["jacobian_determinant_mod_5"] == 2
    assert all(eval_poly_integer(polynomial, lift) % modulus == 0 for polynomial in slice_polynomials)

    free_values: dict[sp.Symbol, int] = {value: 1 for value in free if str(value) in DIAGONAL_NAMES}
    free_values.update(held)
    free_values.update(dict(zip(solve_variables, lift, strict=True)))
    branch_values = [
        eval_expr_mod(substitutions.get(variable, variable), free, free_values, modulus)
        for variable in variables
    ]
    construction_residues = [
        eval_sparse_mod(sparse_poly(equation, variables), branch_values, modulus)
        for equation in equations
    ]
    assert construction_residues == [0] * 42
    assert section["construction_equation_count"] == len(construction_residues) == 42
    assert section["all_construction_residues_zero_mod_625"] is True

    holdout_system = full_weight_finite_system((HOLDOUT_SHAPE,), BASE_ORDERS, gauge_fix=False)
    holdout_variables, holdout_equations = branch_11111_equations(holdout_system)
    assert [str(value) for value in holdout_variables] == [str(value) for value in variables]
    holdout_residues = [
        eval_sparse_mod(sparse_poly(equation, holdout_variables), branch_values, modulus)
        for equation in holdout_equations
    ]
    assert holdout_residues == [0, 0, 250]
    assert section["holdout_shape"] == [3, 3, 3]
    assert section["holdout_orders"] == [4, 6, 8]
    assert section["holdout_residues_mod_625"] == holdout_residues
    assert section["holdout_k8_p_adic_valuation"] == 3
    assert section["holdout_k8_unit_mod_5"] == 2


def main() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    check_envelope(payload)
    check_degree_three_data(payload)
    check_torus_hensel_certificate(payload)
    print("PASS dense degree-three dimensions and lifted syzygies")
    print("PASS exact sparse degree-three selections")
    print("PASS F5 Hensel construction certificate and 3x3x3 k=8 falsification")
    print("PASS")


if __name__ == "__main__":
    main()
