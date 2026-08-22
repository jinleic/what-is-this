"""Exact cubic follow-up for the generic Kac--Ward branch ``11111``.

The experiment has three connected parts:

* it sizes the full degree-three Macaulay block and lifts the six stored
  degree-two syzygies through every degree-at-most-one monomial;
* it exhausts two exact sparse graded degree-three selections over ``Q``;
* it gives a finite ``F_5``/Hensel certificate for a thin-torus construction
  point and an exact ``3x3x3`` order-eight falsification of that local slice.

No floating-point calculation decides any statement in this file.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import DIRECTION_GAUGE_TREE, full_weight_finite_system


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/kac_ward/branch11111_cubic.json"
PRIOR = ROOT / "results/kac_ward/branch11111.json"
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
    started = time.monotonic()
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
        "wall_seconds": time.monotonic() - started,
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


def verify_relation(columns: list[SparseInteger], relation: dict[int, Fraction]) -> None:
    total: SparseRational = {}
    for column, coefficient in relation.items():
        for monomial, value in columns[column].items():
            reduced = total.get(monomial, Fraction(0)) + coefficient * value
            if reduced:
                total[monomial] = reduced
            else:
                total.pop(monomial, None)
    if total:
        raise AssertionError("stored degree-two relation does not vanish")


def prior_syzygy_vectors(prior: dict, degree_two_columns: list[SparseInteger]) -> list[dict[int, Fraction]]:
    records = prior["data"]["degree_two_linear_algebra"]["dependency_relations"]
    vectors = []
    for record in records:
        vector: dict[int, Fraction] = {}
        for term in record["terms"]:
            coefficient = Fraction(
                int(term["coefficient"]["numerator"]),
                int(term["coefficient"]["denominator"]),
            )
            column = int(term["column_index"])
            vector[column] = vector.get(column, Fraction(0)) + coefficient
        vector = {column: coefficient for column, coefficient in vector.items() if coefficient}
        verify_relation(degree_two_columns, vector)
        vectors.append(vector)
    return vectors


def lifted_degree_two_syzygies(
    degree_two_relations: list[dict[int, Fraction]],
    degree_two_multipliers: list[Monomial],
    degree_three_multipliers: list[Monomial],
) -> list[dict[int, Fraction]]:
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
    if len(free) != 16:
        raise AssertionError("thin torus has an unexpected free-coordinate count")
    return free, substitutions


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
    polynomial: SparseInteger, values: list[int], modulus: int
) -> int:
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
    numerator_value = eval_sparse_mod(
        sparse_poly(numerator, variables), [values[variable] for variable in variables], modulus
    )
    denominator_value = eval_sparse_mod(
        sparse_poly(denominator, variables), [values[variable] for variable in variables], modulus
    )
    if not denominator_value % 5:
        raise AssertionError("thin-torus denominator is not a five-adic unit")
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
    polynomials: list[sp.Poly],
    solve_variables: list[sp.Symbol],
    initial: list[int],
    prime: int,
    exponent: int,
) -> tuple[list[int], int, list[list[int]]]:
    jacobian = [
        [sp.Poly(sp.diff(polynomial.as_expr(), variable), *solve_variables, domain=sp.ZZ) for variable in solve_variables]
        for polynomial in polynomials
    ]
    jacobian_mod_prime = [
        [eval_poly_integer(entry, initial) % prime for entry in row]
        for row in jacobian
    ]
    if sp.Matrix(jacobian_mod_prime).det() % prime != 2:
        raise AssertionError("selected five-adic slice is not nonsingular mod five")
    values = list(initial)
    modulus = prime
    for _ in range(1, exponent):
        residuals = [eval_poly_integer(polynomial, values) for polynomial in polynomials]
        if any(residual % modulus for residual in residuals):
            raise AssertionError("Hensel invariant failed")
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


def q_record(value: Fraction | int) -> dict[str, int]:
    rational = Fraction(value)
    return {"numerator": rational.numerator, "denominator": rational.denominator}


def dense_degree_three_record(
    variables: list[sp.Symbol], basis_indices: list[int], basis_polynomials: list[SparseInteger], prior: dict
) -> tuple[dict, list[Monomial], list[Monomial], list[dict[int, Fraction]]]:
    degree_two = multiplier_exponents(len(variables), 2)
    degree_three = multiplier_exponents(len(variables), 3)
    degree_two_columns = macaulay_columns(basis_polynomials, degree_two)
    source_relations = prior_syzygy_vectors(prior, degree_two_columns)
    local_reduction = fraction_reduction(
        degree_two_columns, len(variables), capture_relations=True
    )
    if local_reduction["rank"] != 7365 or len(source_relations) != 6:
        raise AssertionError("prior degree-two structure no longer matches the exact system")
    lifted = lifted_degree_two_syzygies(source_relations, degree_two, degree_three)
    lifted_rank, pivots = sparse_vector_rank(lifted)
    support_started = time.monotonic()
    rows = {
        tuple(left + right for left, right in zip(monomial, multiplier, strict=True))
        for polynomial in basis_polynomials
        for monomial in polynomial
        for multiplier in degree_three
    }
    support_wall = time.monotonic() - support_started
    record = {
        "claim_tag": "[COMPUTATION]",
        "scope": (
            "Exact dimensions for the full multiplier-degree-at-most-three block. "
            "The all-degree p-adic certificate below makes a dense Fraction solve "
            "logically unnecessary; no dense rank equality is claimed."
        ),
        "basis_constraint_indices": basis_indices,
        "variable_count": len(variables),
        "multiplier_degree_bound": 3,
        "multiplier_count": len(degree_three),
        "column_count": len(basis_polynomials) * len(degree_three),
        "basis_term_count": sum(map(len, basis_polynomials)),
        "nonzero_incidence_count": len(degree_three) * sum(map(len, basis_polynomials)),
        "distinct_monomial_row_count": len(rows),
        "support_count_wall_seconds": support_wall,
        "syzygy_pruning": {
            "base_relation_count": len(source_relations),
            "lift_multiplier_count": len(multiplier_exponents(len(variables), 1)),
            "lifted_relation_count": len(lifted),
            "lifted_relation_rank": lifted_rank,
            "relation_sha256": relation_hash(lifted),
            "pivot_sha256": hashlib.sha256(";".join(map(str, pivots)).encode()).hexdigest(),
            "pruned_column_count": len(basis_polynomials) * len(degree_three) - lifted_rank,
            "claim": (
                "[LEMMA] The displayed lifted degree-two syzygies are 156 independent "
                "relations among full degree-three Macaulay columns, so eliminating one "
                "pivot column per relation preserves the dense column span."
            ),
        },
    }
    expected = {
        "multiplier_count": 3276,
        "column_count": 68796,
        "basis_term_count": 1001,
        "nonzero_incidence_count": 3279276,
        "distinct_monomial_row_count": 1612152,
    }
    for key, value in expected.items():
        if record[key] != value:
            raise AssertionError(f"unexpected dense degree-three {key}: {record[key]}")
    if record["syzygy_pruning"]["lifted_relation_rank"] != 156:
        raise AssertionError("unexpected lifted-syzygy rank")
    return record, degree_two, degree_three, source_relations


def sparse_degree_three_records(
    variables: list[sp.Symbol], basis_polynomials: list[SparseInteger], degree_two: list[Monomial]
) -> list[dict]:
    positions = {str(variable): index for index, variable in enumerate(variables)}
    output = []
    for name, selected_names in (
        ("all_degree_le_2_plus_slice6_degree_3", SLICE6_NAMES),
        ("all_degree_le_2_plus_slice8_degree_3", SLICE8_NAMES),
    ):
        extra = multiplier_exponents(
            len(variables),
            3,
            exact_degree=3,
            allowed_positions=[positions[value] for value in selected_names],
        )
        multipliers = degree_two + extra
        columns = macaulay_columns(basis_polynomials, multipliers)
        reduction = fraction_reduction(columns, len(variables), capture_relations=False)
        if reduction["target"]:
            target_record = [
                {
                    "monomial_exponents": list(monomial),
                    "coefficient": q_record(coefficient),
                }
                for monomial, coefficient in sorted(reduction["target"].items())
            ]
        else:
            target_record = []
        record = {
            "claim_tag": "[COMPUTATION]",
            "name": name,
            "degree_three_variable_names": selected_names,
            "degree_two_multiplier_count": len(degree_two),
            "extra_degree_three_multiplier_count": len(extra),
            "multiplier_count": len(multipliers),
            "column_count": len(columns),
            "column_rank": reduction["rank"],
            "dependent_column_count": len(reduction["dependent"]),
            "constant_target_in_span": not reduction["target"],
            "constant_residual": target_record,
            "pivot_monomial_sha256": reduction["pivot_monomial_sha256"],
            "wall_seconds": reduction["wall_seconds"],
            "scope": (
                "[COMPUTATION] This is a selected subspace of degree-three multipliers, "
                "not a rank calculation for the full dense block."
            ),
        }
        if reduction["target"] != {(0,) * len(variables): Fraction(1)}:
            raise AssertionError("unexpected sparse target residual")
        output.append(record)
    return output


def torus_and_hensel_record(
    variables: list[sp.Symbol], equations: list[sp.Expr], labels: list[dict]
) -> tuple[dict, dict]:
    free, substitutions = thin_torus_substitution(variables)
    diagonal = {free[index]: 1 for index in (0, 2, 11)}
    non_diagonal = [variable for variable in free if str(variable) not in DIAGONAL_NAMES]
    if [str(value) for value in non_diagonal] != NONDIAGONAL_NAMES:
        raise AssertionError("unexpected thin-torus coordinate order")

    selected_numerators = []
    primitive_polynomials = []
    contents = []
    support_rows = []
    for index in SELECTED_NUMERATOR_INDICES:
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, _denominator = sp.fraction(reduced)
        numerator = sp.expand(numerator)
        selected_numerators.append(numerator)
        terms = sp.Poly(numerator, *free, domain=sp.ZZ).terms()
        origin = terms[0][0]
        support_rows.extend(
            [exponent - base for exponent, base in zip(monomial, origin, strict=True)]
            for monomial, _coefficient in terms[1:]
        )
        normalized = sp.Poly(
            sp.expand(numerator.subs(diagonal)), *non_diagonal, domain=sp.ZZ
        )
        content, primitive = normalized.primitive()
        contents.append(int(content))
        primitive_polynomials.append(primitive)

    support_matrix = sp.Matrix(support_rows)
    support_rank = support_matrix.rank()
    kernel = [primitive_integer_vector(vector) for vector in support_matrix.nullspace()]
    diagonal_minor = sp.Matrix(kernel)[:, [0, 2, 11]].det()
    if support_rank != 13 or abs(int(diagonal_minor)) != 1:
        raise AssertionError("thin-torus action computation changed")

    full_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in FULL_SHAPES
    ]
    nonzero_numerators = []
    for index in full_indices:
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, _denominator = sp.fraction(reduced)
        if numerator:
            nonzero_numerators.append(sp.expand(numerator))
    if len(nonzero_numerators) != 10:
        raise AssertionError("unexpected number of nonzero three-dimensional rows")
    if len(independent_indices(nonzero_numerators, free)) != 6:
        raise AssertionError("three-dimensional Laurent span rank changed")
    if len(independent_indices(selected_numerators, free)) != 6:
        raise AssertionError("selected numerator basis no longer spans")
    primary_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in PRIMARY_SHAPES
    ]
    if not all(sp.cancel(equations[index].subs(substitutions)) == 0 for index in primary_indices):
        raise AssertionError("primary equations no longer vanish on the thin torus")

    solve_variables = non_diagonal[:6]
    held_values = dict(zip(non_diagonal[6:], HENSEL_POINT_MOD_5[6:], strict=True))
    slice_polynomials = [
        sp.Poly(polynomial.as_expr().subs(held_values), *solve_variables, domain=sp.ZZ)
        for polynomial in primitive_polynomials
    ]
    initial = HENSEL_POINT_MOD_5[:6]
    if any(eval_poly_integer(polynomial, initial) % 5 for polynomial in slice_polynomials):
        raise AssertionError("declared F5 point does not solve primitive slice equations")
    lifted, modulus, jacobian_mod_5 = hensel_lift_slice(
        slice_polynomials, solve_variables, initial, 5, 4
    )
    if lifted != [226, 197, 98, 63, 161, 216] or modulus != 625:
        raise AssertionError("finite Hensel lift changed")
    if any(eval_poly_integer(polynomial, lifted) % modulus for polynomial in slice_polynomials):
        raise AssertionError("finite Hensel lift does not solve slice equations")

    free_values: dict[sp.Symbol, int] = {value: 1 for value in free if str(value) in DIAGONAL_NAMES}
    free_values.update(held_values)
    free_values.update(dict(zip(solve_variables, lifted, strict=True)))
    branch_values = [
        eval_expr_mod(substitutions.get(variable, variable), free, free_values, modulus)
        for variable in variables
    ]
    construction_residues = [
        eval_sparse_mod(sparse_poly(equation, variables), branch_values, modulus)
        for equation in equations
    ]
    if construction_residues != [0] * len(equations):
        raise AssertionError("finite lift misses a construction equation")
    holdout_system = full_weight_finite_system((HOLDOUT_SHAPE,), BASE_ORDERS, gauge_fix=False)
    holdout_variables, holdout_equations = branch_11111_equations(holdout_system)
    if [str(value) for value in holdout_variables] != [str(value) for value in variables]:
        raise AssertionError("holdout variable order changed")
    holdout_residues = [
        eval_sparse_mod(sparse_poly(equation, holdout_variables), branch_values, modulus)
        for equation in holdout_equations
    ]
    if holdout_residues != [0, 0, 250]:
        raise AssertionError("expected order-eight holdout residue changed")

    torus_record = {
        "claim_tag": "[COMPUTATION]",
        "free_variable_order": [str(value) for value in free],
        "non_diagonal_variable_order": [str(value) for value in non_diagonal],
        "selected_numerator_indices": SELECTED_NUMERATOR_INDICES,
        "support_difference_rank": support_rank,
        "torus_action_kernel_vectors": kernel,
        "diagonal_normalization_variables": DIAGONAL_NAMES,
        "diagonal_normalization_minor_abs": abs(int(diagonal_minor)),
        "claim": (
            "[LEMMA] The six selected thin-torus numerators have a three-dimensional "
            "integer torus action.  Its diagonal-coordinate minor is unimodular, so the "
            "three diagonal free coordinates may be normalized to one on this action chart."
        ),
        "scope": (
            "This is a coordinate reduction for the displayed finite thin-torus system; "
            "it is not an all-box identity."
        ),
    }
    hensel_record = {
        "claim_tag": "[THEOREM]",
        "prime": 5,
        "modulus": modulus,
        "numerator_constraint_indices": SELECTED_NUMERATOR_INDICES,
        "primitive_contents": contents,
        "solve_variable_names": [str(value) for value in solve_variables],
        "held_variable_names": [str(value) for value in non_diagonal[6:]],
        "initial_non_diagonal_values_mod_5": HENSEL_POINT_MOD_5,
        "lifted_solve_values_mod_625": lifted,
        "jacobian_mod_5": jacobian_mod_5,
        "jacobian_determinant_mod_5": int(sp.Matrix(jacobian_mod_5).det()) % 5,
        "construction_equation_count": len(construction_residues),
        "all_construction_residues_zero_mod_625": not any(construction_residues),
        "holdout_shape": list(HOLDOUT_SHAPE),
        "holdout_orders": list(BASE_ORDERS),
        "holdout_residues_mod_625": holdout_residues,
        "holdout_k8_p_adic_valuation": 3,
        "holdout_k8_unit_mod_5": 2,
        "claim": (
            "[THEOREM] The primitive six-equation thin-torus slice has a nonsingular "
            "F_5 point.  Multivariate Hensel lifting therefore gives a Q_5 point of all "
            "42 construction equations.  At that unique local lift, the 3x3x3 order-eight "
            "residual is congruent to 2*5^3 modulo 5^4 and is therefore nonzero."
        ),
        "consequences": [
            "[THEOREM] The recorded 42-equation ideal is proper over Q, so no rational Nullstellensatz certificate exists at any multiplier degree for this catalogue.",
            "[FALSIFIED] This nonsingular five-adic construction specialization cannot satisfy the 3x3x3 order-eight holdout.",
        ],
        "scope": (
            "The Hensel theorem concerns one zero-dimensional section after seven thin-torus "
            "coordinates are held.  It neither supplies an explicit quartic number-field point "
            "nor decides every construction point or the full scalar family."
        ),
    }
    return torus_record, hensel_record


def main() -> None:
    prior = json.loads(PRIOR.read_text(encoding="utf-8"))
    if prior["data"]["branch"] != "11111":
        raise AssertionError("wrong prior branch artifact")
    variables, equations, labels, _system = construction_data()
    basis_indices = independent_indices(equations, variables)
    if len(basis_indices) != 21:
        raise AssertionError("unexpected exact constraint-span rank")
    basis_polynomials = [sparse_poly(equations[index], variables) for index in basis_indices]

    dense_record, degree_two, _degree_three, _relations = dense_degree_three_record(
        variables, basis_indices, basis_polynomials, prior
    )
    sparse_records = sparse_degree_three_records(
        variables, basis_polynomials, degree_two
    )
    torus_record, hensel_record = torus_and_hensel_record(variables, equations, labels)

    payload = {
        "provenance": {
            "script": "experiments/e90_kw_cubic.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "no_float_decisions": True,
            "method": (
                "exact Fraction sparse linear algebra, integer Laurent support analysis, "
                "exact F_5 Gaussian elimination, and finite multivariate Hensel lifting"
            ),
            "prior_branch_artifact": "results/kac_ward/branch11111.json",
            "prior_branch_artifact_sha256": hashlib.sha256(PRIOR.read_bytes()).hexdigest(),
        },
        "data": {
            "claim_tag": "[UNRESOLVED]",
            "branch": "11111",
            "status": "UNRESOLVED",
            "headline": (
                "[THEOREM] A nonsingular F_5 thin-torus section Hensel-lifts to a Q_5 point "
                "of the exact 42-constraint catalogue, excluding rational Nullstellensatz "
                "certificates at every degree; [FALSIFIED] that local section fails the exact "
                "3x3x3 order-eight holdout; [UNRESOLVED] the full branch remains undecided."
            ),
            "construction_shapes": [list(shape) for shape in CONSTRUCTION_SHAPES],
            "orders": list(BASE_ORDERS),
            "variable_order": [str(value) for value in variables],
            "exact_Q_constraint_span_rank": len(basis_indices),
            "dense_degree_three": dense_record,
            "sparse_degree_three": sparse_records,
            "torus_slice": torus_record,
            "hensel_five_adic_section": hensel_record,
            "quartic_extension_search": {
                "claim_tag": "[UNRESOLVED]",
                "statement": (
                    "No exact Q(alpha) coordinate presentation of degree four was constructed. "
                    "The Hensel certificate is a p-adic nonemptiness and holdout-falsification "
                    "certificate, not an explicit quartic point."
                ),
                "historical_modular_groebner_wall": {
                    "claim_tag": "[COMPUTATION]",
                    "field": "F_5",
                    "order": "grevlex",
                    "slice_non_diagonal_values_mod_5": [3, 4, 1, 1, 1, 1, 2, 2, 3, 2, 4, 1, 1],
                    "solve_variable_names": [
                        "u_py_px",
                        "u_py_mz",
                        "u_my_px",
                        "u_my_mx",
                        "u_pz_mx",
                        "u_mz_mx",
                    ],
                    "wall_seconds": 600,
                    "status": "timed_out_without_basis",
                    "scope": "[UNRESOLVED] A timeout is not an emptiness or existence theorem.",
                },
            },
        },
        "checks": [
            {
                "name": "dense_degree_three_sized_and_pruned",
                "passed": (
                    dense_record["column_count"] == 68796
                    and dense_record["syzygy_pruning"]["lifted_relation_rank"] == 156
                ),
                "detail": {
                    "rows": dense_record["distinct_monomial_row_count"],
                    "columns": dense_record["column_count"],
                    "pruned_columns": dense_record["syzygy_pruning"]["pruned_column_count"],
                },
            },
            {
                "name": "sparse_degree_three_fraction_checks",
                "passed": all(not record["constant_target_in_span"] for record in sparse_records),
                "detail": {
                    record["name"]: {
                        "rank": record["column_rank"],
                        "columns": record["column_count"],
                    }
                    for record in sparse_records
                },
            },
            {
                "name": "five_adic_hensel_holdout_falsification",
                "passed": (
                    hensel_record["jacobian_determinant_mod_5"] != 0
                    and hensel_record["all_construction_residues_zero_mod_625"]
                    and hensel_record["holdout_residues_mod_625"] == [0, 0, 250]
                ),
                "detail": {
                    "jacobian_determinant_mod_5": hensel_record["jacobian_determinant_mod_5"],
                    "holdout_k8_p_adic_valuation": hensel_record["holdout_k8_p_adic_valuation"],
                    "holdout_k8_unit_mod_5": hensel_record["holdout_k8_unit_mod_5"],
                },
            },
        ],
    }
    if not all(check["passed"] for check in payload["checks"]):
        raise AssertionError("one cubic branch check failed")
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "PASS dense degree-three "
        f"columns={dense_record['column_count']} rows={dense_record['distinct_monomial_row_count']} "
        f"syzygy_pruned={dense_record['syzygy_pruning']['pruned_column_count']}"
    )
    for record in sparse_records:
        print(
            "PASS sparse degree-three "
            f"{record['name']} rank={record['column_rank']} columns={record['column_count']}"
        )
    print(
        "PASS F5 Hensel "
        f"det={hensel_record['jacobian_determinant_mod_5']} "
        f"holdout_k8={hensel_record['holdout_residues_mod_625'][2]} mod 625"
    )
    print("PASS")


if __name__ == "__main__":
    main()
