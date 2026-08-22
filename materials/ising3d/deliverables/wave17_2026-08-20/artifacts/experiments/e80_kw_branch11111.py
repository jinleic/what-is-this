"""Exact degree-two Kac--Ward investigation on the generic branch ``11111``.

The calculation deliberately works with the same finite 42-equation catalogue
as the preceding degree-one certificate, then constructs every multiplier of
total degree at most two.  All decisive characteristic-zero reduction uses
``fractions.Fraction`` sparse Gaussian elimination.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    OPPOSITE_DIRECTION,
    full_weight_finite_system,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/kac_ward/branch11111.json"
BASE_ORDERS = (4, 6, 8)
PRIMARY_SHAPES = tuple(
    sorted(
        set(itertools.permutations((2, 2, 1)))
        | set(itertools.permutations((3, 2, 1)))
    )
)
FULL_SHAPES = ((3, 3, 2), (2, 2, 3), (2, 2, 2), (3, 2, 2), (2, 3, 2))
CONSTRUCTION_SHAPES = PRIMARY_SHAPES + FULL_SHAPES

Monomial = tuple[int, ...]
SparseVector = dict[Monomial, Fraction]
SparseIntegerVector = dict[Monomial, int]


def branch_11111_equations(system) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    """Take the generic nonzero gauge-tree chart directly from an ungauged system."""
    fixed = set(DIRECTION_GAUGE_TREE)
    substitutions = {symbol: 1 for pair, symbol in system.all_symbols if pair in fixed}
    variables = [symbol for pair, symbol in system.all_symbols if pair not in fixed]
    return variables, [sp.expand(expression.subs(substitutions)) for expression in system.equations]


def monomial_key(monomial: Monomial) -> tuple[int, Monomial]:
    return (sum(monomial), monomial)


def sparse_poly(expression: sp.Expr, variables: list[sp.Symbol]) -> SparseIntegerVector:
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return {
        tuple(monomial): int(coefficient)
        for monomial, coefficient in polynomial.terms()
        if coefficient
    }


def sparse_matrix_hash(columns: list[SparseIntegerVector]) -> str:
    digest = hashlib.sha256()
    for index, column in enumerate(columns):
        digest.update(f"C{index}:".encode())
        for monomial, coefficient in sorted(column.items()):
            digest.update((f"{coefficient}@{','.join(map(str, monomial))};").encode())
    return digest.hexdigest()


def reduce_fraction_columns(
    columns: list[SparseIntegerVector], variable_count: int
) -> dict:
    """Exactly reduce columns and retain every zero-column relation over Q."""
    basis: dict[Monomial, SparseVector] = {}
    basis_relations: dict[Monomial, dict[int, Fraction]] = {}
    dependent_indices: list[int] = []
    dependency_relations: list[dict] = []
    pivot_order: list[Monomial] = []
    started = time.monotonic()
    for index, source in enumerate(columns):
        vector: SparseVector = {
            monomial: Fraction(coefficient)
            for monomial, coefficient in source.items()
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
            for column, value in basis_relations[pivot].items():
                reduced = relation.get(column, Fraction(0)) - coefficient * value
                if reduced:
                    relation[column] = reduced
                else:
                    relation.pop(column, None)
        if not vector:
            dependent_indices.append(index)
            dependency_relations.append(
                {
                    "dependent_column_index": index,
                    "terms": [
                        {
                            "column_index": column,
                            "coefficient": {
                                "numerator": value.numerator,
                                "denominator": value.denominator,
                            },
                        }
                        for column, value in sorted(relation.items())
                        if value
                    ],
                }
            )
    target: SparseVector = {(0,) * variable_count: Fraction(1)}
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
        "column_rank": len(basis),
        "dependent_column_indices": dependent_indices,
        "dependency_relations": dependency_relations,
        "constant_target_in_span": not target,
        "constant_residual": [
            {
                "monomial_exponents": list(monomial),
                "coefficient": {
                    "numerator": value.numerator,
                    "denominator": value.denominator,
                },
            }
            for monomial, value in sorted(target.items())
        ],
        "pivot_monomial_sha256": hashlib.sha256(
            "".join(",".join(map(str, monomial)) + ";" for monomial in pivot_order).encode()
        ).hexdigest(),
        "wall_seconds": time.monotonic() - started,
    }


def independent_constraint_indices(
    equations: list[sp.Expr], variables: list[sp.Symbol]
) -> list[int]:
    columns = [sparse_poly(expression, variables) for expression in equations]
    basis: dict[Monomial, SparseVector] = {}
    indices: list[int] = []
    for index, source in enumerate(columns):
        vector: SparseVector = {monomial: Fraction(coefficient) for monomial, coefficient in source.items()}
        while vector:
            pivot = max(vector, key=monomial_key)
            coefficient = vector[pivot]
            existing = basis.get(pivot)
            if existing is None:
                inverse = Fraction(1, 1) / coefficient
                basis[pivot] = {monomial: value * inverse for monomial, value in vector.items() if value}
                indices.append(index)
                break
            for monomial, value in existing.items():
                reduced = vector.get(monomial, Fraction(0)) - coefficient * value
                if reduced:
                    vector[monomial] = reduced
                else:
                    vector.pop(monomial, None)
    return indices


def multiplier_exponents(variable_count: int, degree_bound: int = 2) -> list[Monomial]:
    if degree_bound != 2:
        raise ValueError("this recorded catalogue is fixed at multiplier degree two")
    zero = (0,) * variable_count
    output = [zero]
    output.extend(
        tuple(1 if coordinate == index else 0 for coordinate in range(variable_count))
        for index in range(variable_count)
    )
    output.extend(
        tuple(int(first == index) + int(second == index) for index in range(variable_count))
        for first in range(variable_count)
        for second in range(first, variable_count)
    )
    return output


def macaulay_columns(
    polynomials: list[SparseIntegerVector], multipliers: list[Monomial]
) -> list[SparseIntegerVector]:
    return [
        {
            tuple(exponent + multiplier for exponent, multiplier in zip(monomial, factor, strict=True)): coefficient
            for monomial, coefficient in polynomial.items()
        }
        for polynomial in polynomials
        for factor in multipliers
    ]


def column_metadata(
    column_index: int,
    basis_indices: list[int],
    multipliers: list[Monomial],
    labels: list[dict],
) -> dict:
    basis_position, multiplier_position = divmod(column_index, len(multipliers))
    basis_index = basis_indices[basis_position]
    return {
        "basis_constraint_index": basis_index,
        "basis_constraint_label": labels[basis_index],
        "multiplier_exponents": list(multipliers[multiplier_position]),
    }


def q_record(value: sp.Rational | Fraction | int) -> dict[str, int]:
    rational = sp.Rational(value)
    return {"numerator": int(rational.p), "denominator": int(rational.q)}


def thin_radical_relations(variables: list[sp.Symbol]) -> list[sp.Expr]:
    """Rebuild the nine saturated thin-plane binomials, without reading e62."""
    symbol = {str(value): value for value in variables}
    return [
        symbol["u_mx_my"] * symbol["u_my_px"] + 1,
        symbol["u_mx_py"] * symbol["u_my_mx"] * symbol["u_py_px"] + 1,
        symbol["u_mx_mx"] * symbol["u_px_px"] - 1,
        symbol["u_mx_mz"] * symbol["u_mz_px"] * symbol["u_pz_mx"] + 1,
        symbol["u_mx_pz"] * symbol["u_mz_mx"] * symbol["u_pz_px"] + 1,
        symbol["u_mz_mz"] * symbol["u_pz_pz"] - 1,
        symbol["u_my_mz"]
        * symbol["u_mz_py"]
        * symbol["u_py_pz"]
        * symbol["u_pz_my"]
        + 1,
        symbol["u_my_pz"]
        * symbol["u_mz_my"]
        * symbol["u_py_mz"]
        * symbol["u_pz_py"]
        + 1,
        symbol["u_my_my"] * symbol["u_py_py"] - 1,
    ]


def thin_torus_substitution(
    variables: list[sp.Symbol],
) -> tuple[list[sp.Symbol], dict[sp.Symbol, sp.Expr]]:
    """Solve each radical binomial for one distinct active coordinate."""
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
        raise AssertionError(f"thin torus should have 16 free coordinates, got {len(free)}")
    return free, substitutions


def denominator_record(denominator: sp.Expr, variables: list[sp.Symbol]) -> dict:
    polynomial = sp.Poly(denominator, *variables, domain=sp.QQ)
    terms = polynomial.terms()
    if len(terms) != 1:
        raise AssertionError(f"torus denominator is not a monomial: {denominator}")
    exponents, coefficient = terms[0]
    return {
        "coefficient": q_record(coefficient),
        "exponents": list(exponents),
    }


def torus_elimination_record(
    variables: list[sp.Symbol],
    equations: list[sp.Expr],
    labels: list[dict],
) -> dict:
    """Push all finite 3D equations to the explicit 16-dimensional thin torus."""
    free, substitutions = thin_torus_substitution(variables)
    radical = thin_radical_relations(variables)
    if not all(sp.cancel(relation.subs(substitutions)) == 0 for relation in radical):
        raise AssertionError("declared torus substitution does not solve the radical")

    primary_indices = [
        index
        for index, label in enumerate(labels)
        if tuple(label["shape"]) in PRIMARY_SHAPES
    ]
    full_indices = [
        index
        for index, label in enumerate(labels)
        if tuple(label["shape"]) in FULL_SHAPES
    ]
    primary_zero = [
        sp.cancel(equations[index].subs(substitutions)) == 0
        for index in primary_indices
    ]
    if not all(primary_zero):
        raise AssertionError("a primary thin constraint survives the radical torus")

    rows = []
    nonzero_numerators: list[sp.Expr] = []
    nonzero_row_indices: list[int] = []
    for index in full_indices:
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, denominator = sp.fraction(reduced)
        numerator = sp.expand(numerator)
        terms = sparse_poly(numerator, free)
        row = {
            "constraint_index": index,
            "constraint_label": labels[index],
            "is_zero_on_torus": not terms,
            "numerator_term_count": len(terms),
            "numerator_total_degree_max": max(map(sum, terms), default=0),
            "numerator_sha256": sparse_matrix_hash([terms]),
            "denominator": denominator_record(denominator, free),
        }
        rows.append(row)
        if terms:
            nonzero_numerators.append(numerator)
            nonzero_row_indices.append(index)

    independent = independent_constraint_indices(nonzero_numerators, free)
    return {
        "claim_tag": "[COMPUTATION]",
        "claim": (
            "[COMPUTATION] The nine exact thin-radical binomials were solved as "
            "Laurent-coordinate substitutions, and every 3D construction equation "
            "was reduced on that 16-dimensional torus."
        ),
        "input_coordinate_count": len(variables),
        "radical_relation_count": len(radical),
        "free_torus_coordinate_count": len(free),
        "free_variable_order": [str(variable) for variable in free],
        "pivot_substitutions": [
            {"variable": str(variable), "value": str(value)}
            for variable, value in substitutions.items()
        ],
        "radical_binomials": [str(relation) for relation in radical],
        "radical_substitution_verified": True,
        "primary_constraint_count": len(primary_indices),
        "primary_constraints_vanish_on_torus": all(primary_zero),
        "three_dimensional_constraint_count": len(full_indices),
        "three_dimensional_laurent_rows": rows,
        "nonzero_three_dimensional_constraint_count": len(nonzero_numerators),
        "cleared_numerator_Q_span_rank": len(independent),
        "cleared_numerator_independent_constraint_indices": [
            nonzero_row_indices[index] for index in independent
        ],
        "scope": (
            "[UNRESOLVED] Clearing a Laurent monomial denominator is valid only on "
            "this torus chart. The recorded exact support and Q-linear ranks do not "
            "by themselves prove an algebraic zero or emptiness."
        ),
    }


def pair_text(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def reversal_pair(pair: tuple[int, int]) -> tuple[int, int]:
    return (OPPOSITE_DIRECTION[pair[1]], OPPOSITE_DIRECTION[pair[0]])


def evaluate_mu4(
    polynomial: SparseIntegerVector, exponents: list[int]
) -> tuple[int, int]:
    """Evaluate an integer polynomial at coordinates i**exponents[j]."""
    power_sums = [0, 0, 0, 0]
    for monomial, coefficient in polynomial.items():
        phase = sum(
            exponent * coordinate
            for exponent, coordinate in zip(monomial, exponents, strict=True)
        ) % 4
        power_sums[phase] += coefficient
    return (power_sums[0] - power_sums[2], power_sums[1] - power_sums[3])


def mu4_value_record(value: tuple[int, int]) -> dict[str, int]:
    return {"real": value[0], "imaginary": value[1]}


def root_of_unity_chart(
    system, variables: list[sp.Symbol]
) -> tuple[
    list[tuple[int, int]],
    list[tuple[int, int] | None],
    list[int],
    list[dict],
]:
    """Build the branch-normalized reversal-conjugate Q(i) coordinate chart."""
    symbol_by_pair = dict(system.all_symbols)
    pair_by_symbol = {symbol: pair for pair, symbol in system.all_symbols}
    tree = set(DIRECTION_GAUGE_TREE)
    representatives: list[tuple[int, int]] = []
    coordinate_rule: dict[tuple[int, int], tuple[int, int] | None] = {}
    fixed_orbits: list[dict] = []
    visited: set[tuple[int, int]] = set()
    for pair in sorted(symbol_by_pair):
        if pair in visited:
            continue
        reverse = reversal_pair(pair)
        if reverse not in symbol_by_pair:
            raise AssertionError(f"missing reversal partner for {pair}")
        visited.update((pair, reverse))
        representative, partner = sorted((pair, reverse))
        if representative in tree or partner in tree:
            coordinate_rule[representative] = None
            coordinate_rule[partner] = None
            fixed_orbits.append(
                {
                    "representative": pair_text(representative),
                    "partner": pair_text(partner),
                    "fixed_value": "1",
                }
            )
            continue
        coordinate = len(representatives)
        representatives.append(representative)
        coordinate_rule[representative] = (coordinate, 1)
        coordinate_rule[partner] = (coordinate, -1)

    active_rules = [coordinate_rule[pair_by_symbol[variable]] for variable in variables]
    if len(representatives) != 11:
        raise AssertionError(f"expected 11 free reversal orbits, got {len(representatives)}")
    return representatives, active_rules, [0] * len(variables), fixed_orbits


def radical_mu4_congruences(
    radical: list[sp.Expr],
    variables: list[sp.Symbol],
    active_rules: list[tuple[int, int] | None],
    free_count: int,
) -> list[dict]:
    """Translate every monomial-equals-sign radical relation to Z/4Z."""
    zero = (0,) * len(variables)
    output = []
    for relation in radical:
        polynomial = sparse_poly(relation, variables)
        nonconstant = [
            (monomial, coefficient)
            for monomial, coefficient in polynomial.items()
            if monomial != zero
        ]
        constant = polynomial.get(zero, 0)
        if len(nonconstant) != 1 or abs(nonconstant[0][1]) != 1 or abs(constant) != 1:
            raise AssertionError(f"not a signed monomial radical relation: {relation}")
        monomial, coefficient = nonconstant[0]
        target_value = -constant // coefficient
        if target_value not in (-1, 1):
            raise AssertionError("unexpected root-of-unity target")
        coefficients = [0] * free_count
        for variable_index, power in enumerate(monomial):
            rule = active_rules[variable_index]
            if rule is not None:
                coordinate, sign = rule
                coefficients[coordinate] = (
                    coefficients[coordinate] + sign * power
                ) % 4
        output.append(
            {
                "relation": str(relation),
                "coefficients_mod_4": coefficients,
                "target_mod_4": 2 if target_value == -1 else 0,
            }
        )
    return output


def mod4_unit_rref(
    congruences: list[dict], variable_count: int
) -> tuple[list[list[int]], list[int]]:
    """Row-reduce using only odd (hence invertible) pivots in Z/4Z."""
    rows = [
        [int(value) % 4 for value in relation["coefficients_mod_4"]]
        + [int(relation["target_mod_4"]) % 4]
        for relation in congruences
    ]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(variable_count):
        candidate = next(
            (
                row
                for row in range(pivot_row, len(rows))
                if rows[row][column] % 2
            ),
            None,
        )
        if candidate is None:
            continue
        rows[pivot_row], rows[candidate] = rows[candidate], rows[pivot_row]
        inverse = 1 if rows[pivot_row][column] == 1 else 3
        rows[pivot_row] = [
            value * inverse % 4 for value in rows[pivot_row]
        ]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % 4
                    for value, pivot_value in zip(
                        rows[row], rows[pivot_row], strict=True
                    )
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivots


def mod4_presolved_assignments(
    congruences: list[dict], variable_count: int
) -> tuple[itertools.product, list[int], list[int], list[list[int]]]:
    """Enumerate only the free coordinates after exact unit-pivot reduction."""
    rows, pivots = mod4_unit_rref(congruences, variable_count)
    free_columns = [column for column in range(variable_count) if column not in pivots]
    return itertools.product(range(4), repeat=len(free_columns)), pivots, free_columns, rows


def reversal_mu4_scan_record(
    system,
    variables: list[sp.Symbol],
    equations: list[sp.Expr],
    labels: list[dict],
) -> dict:
    """Exhaust the branch-normalized reversal-conjugate fourth-root chart."""
    representatives, active_rules, fixed_exponents, fixed_orbits = root_of_unity_chart(
        system, variables
    )
    radical = thin_radical_relations(variables)
    congruences = radical_mu4_congruences(
        radical, variables, active_rules, len(representatives)
    )
    radical_sparse = [sparse_poly(relation, variables) for relation in radical]
    construction_sparse = [sparse_poly(equation, variables) for equation in equations]
    scan_values, pivots, free_columns, reduced_rows = mod4_presolved_assignments(
        congruences, len(representatives)
    )

    solution_samples: list[dict] = []
    passing_construction: list[dict] = []
    digest = hashlib.sha256()
    raw_assignment_count = 4 ** len(representatives)
    solution_count = 0
    for free_values in scan_values:
        free_exponents = [0] * len(representatives)
        for column, value in zip(free_columns, free_values, strict=True):
            free_exponents[column] = value
        for row, pivot in zip(reduced_rows, pivots, strict=False):
            free_exponents[pivot] = (
                row[-1]
                - sum(
                    row[column] * free_exponents[column]
                    for column in free_columns
                )
            ) % 4
        if any(
            sum(
                coefficient * exponent
                for coefficient, exponent in zip(
                    relation["coefficients_mod_4"], free_exponents, strict=True
                )
            )
            % 4
            != relation["target_mod_4"]
            for relation in congruences
        ):
            continue
        active_exponents = [
            fixed
            if rule is None
            else (rule[1] * free_exponents[rule[0]]) % 4
            for fixed, rule in zip(fixed_exponents, active_rules, strict=True)
        ]
        radical_residuals = [
            evaluate_mu4(polynomial, active_exponents) for polynomial in radical_sparse
        ]
        if any(residual != (0, 0) for residual in radical_residuals):
            raise AssertionError("mod-four radical solve disagrees with exact Q(i) evaluation")
        first_failure = None
        for index, polynomial in enumerate(construction_sparse):
            residual = evaluate_mu4(polynomial, active_exponents)
            if residual != (0, 0):
                first_failure = (index, residual)
                break
        candidate = {
            "free_exponents_mod_4": free_exponents,
            "active_exponents_in_variable_order_mod_4": active_exponents,
            "first_nonzero_construction_residual": (
                None
                if first_failure is None
                else {
                    "constraint_index": first_failure[0],
                    "constraint_label": labels[first_failure[0]],
                    "value": mu4_value_record(first_failure[1]),
                }
            ),
        }
        solution_count += 1
        if len(solution_samples) < 16:
            solution_samples.append(candidate)
        digest.update(
            (
                ",".join(map(str, free_exponents))
                + ":"
                + (
                    "PASS"
                    if first_failure is None
                    else f"{first_failure[0]}:{first_failure[1][0]},{first_failure[1][1]}"
                )
                + "\n"
            ).encode()
        )
        if first_failure is None:
            passing_construction.append(candidate)

    holdout_records = []
    if passing_construction:
        holdout_system = full_weight_finite_system(
            ((3, 3, 3),), BASE_ORDERS, gauge_fix=False
        )
        holdout_variables, holdout_equations = branch_11111_equations(holdout_system)
        if [str(value) for value in holdout_variables] != [str(value) for value in variables]:
            raise AssertionError("holdout variable order differs from construction order")
        holdout_sparse = [
            sparse_poly(equation, variables) for equation in holdout_equations
        ]
        for candidate in passing_construction:
            residuals = [
                evaluate_mu4(
                    polynomial, candidate["active_exponents_in_variable_order_mod_4"]
                )
                for polynomial in holdout_sparse
            ]
            holdout_records.append(
                {
                    "free_exponents_mod_4": candidate["free_exponents_mod_4"],
                    "shape": [3, 3, 3],
                    "orders": list(BASE_ORDERS),
                    "residuals": [mu4_value_record(value) for value in residuals],
                    "passed": not any(value != (0, 0) for value in residuals),
                }
            )

    rref_hash = hashlib.sha256(
        "\n".join(",".join(map(str, row)) for row in reduced_rows).encode()
    ).hexdigest()
    return {
        "claim_tag": "[COMPUTATION]",
        "claim": (
            "[COMPUTATION] The finite Q(i) fourth-root slice of the "
            "reversal-conjugate branch-11111 thin-torus chart was exhaustively "
            "enumerated with exact Gaussian-integer residuals."
        ),
        "additional_chart_assumptions": [
            "U(-d',-d)=conjugate(U(d,d'))",
            "every weight lies in {1, i, -1, -i}",
        ],
        "reversal_free_orbit_count": len(representatives),
        "reversal_free_orbit_representatives": [
            pair_text(pair) for pair in representatives
        ],
        "tree_forced_reversal_orbits": fixed_orbits,
        "raw_assignment_count": raw_assignment_count,
        "thin_radical_congruences_mod_4": congruences,
        "mod_4_unit_pivot_columns": pivots,
        "mod_4_free_columns": free_columns,
        "mod_4_presolved_assignment_count": 4 ** len(free_columns),
        "mod_4_unit_rref_sha256": rref_hash,
        "forced_holonomy_consequences": {
            "oriented_plaquette_products": [-1] * 6,
            "sum_of_oriented_plaquette_products": -6,
            "reversal_conjugation_plane_real_parts": [-1, -1, -1],
        },
        "radical_solution_count": solution_count,
        "radical_solution_samples": solution_samples,
        "first_failure_sha256": digest.hexdigest(),
        "construction_passing_candidate_count": len(passing_construction),
        "holdout_3x3x3_records": holdout_records,
        "surviving_candidate_count": sum(
            int(record["passed"]) for record in holdout_records
        ),
        "scope": (
            "[UNRESOLVED] This is an exact exhaustion only of the displayed finite "
            "mu_4 subset of the reversal-conjugate complex chart. It neither proves "
            "nor disproves a general complex solution."
        ),
    }


def main() -> None:
    system = full_weight_finite_system(CONSTRUCTION_SHAPES, BASE_ORDERS, gauge_fix=False)
    variables, equations = branch_11111_equations(system)
    labels = [{"shape": list(shape), "order": order} for shape, order in system.labels]
    basis_indices = independent_constraint_indices(equations, variables)
    basis_polynomials = [sparse_poly(equations[index], variables) for index in basis_indices]
    multipliers = multiplier_exponents(len(variables))
    columns = macaulay_columns(basis_polynomials, multipliers)
    exact = reduce_fraction_columns(columns, len(variables))
    for relation in exact["dependency_relations"]:
        for term in relation["terms"]:
            term.update(
                column_metadata(
                    int(term["column_index"]), basis_indices, multipliers, labels
                )
            )
    torus = torus_elimination_record(variables, equations, labels)
    mu4 = reversal_mu4_scan_record(system, variables, equations, labels)
    degree_two_passed = (
        len(basis_indices) == 21
        and len(multipliers) == 351
        and len(columns) == 7371
        and len(exact["dependent_column_indices"]) == 6
        and exact["column_rank"] == len(columns) - 6
        and not exact["constant_target_in_span"]
    )
    torus_passed = (
        torus["radical_substitution_verified"]
        and torus["primary_constraints_vanish_on_torus"]
    )
    mu4_passed = (
        mu4["construction_passing_candidate_count"] == 0
        or len(mu4["holdout_3x3x3_records"])
        == mu4["construction_passing_candidate_count"]
    )
    if not (degree_two_passed and torus_passed and mu4_passed):
        raise AssertionError("one exact branch-11111 check failed")
    payload = {
        "provenance": {
            "script": "experiments/e80_kw_branch11111.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "method": (
                "exact sparse Fraction Gaussian elimination over Q on degree-two "
                "Macaulay columns; exact Laurent reduction through the stored thin "
                "torus; exhaustive Q(i) fourth-root reversal-chart scan"
            ),
            "no_float_decisions": True,
        },
        "data": {
            "claim_tag": "[UNRESOLVED]",
            "branch": "11111",
            "status": "UNRESOLVED",
            "headline": (
                "[THEOREM] No multiplier-degree<=2 rational Nullstellensatz "
                "certificate exists for the recorded 42-constraint system; "
                "[UNRESOLVED] branch 11111 remains undecided over Q/C."
            ),
            "construction_shapes": [list(shape) for shape in CONSTRUCTION_SHAPES],
            "orders": list(BASE_ORDERS),
            "constraint_count": len(equations),
            "variable_order": [str(variable) for variable in variables],
            "exact_Q_constraint_span_rank": len(basis_indices),
            "basis_constraint_indices": basis_indices,
            "basis_constraint_labels": [labels[index] for index in basis_indices],
            "multiplier_degree_bound": 2,
            "multiplier_count": len(multipliers),
            "column_count": len(columns),
            "integer_macaulay_matrix_sha256": sparse_matrix_hash(columns),
            "degree_two_linear_algebra": {
                **exact,
                "observed_syzygy_count": len(exact["dependent_column_indices"]),
                "claim_tag": "[THEOREM]",
                "claim": (
                    "[THEOREM] For this exact 42-constraint branch-11111 catalogue, "
                    "the constant is absent from the Q-span of all multiplier-degree "
                    "at-most-two Macaulay columns."
                ),
                "rank_consequence": (
                    "The six stored exact zero-reduction relations give rank <= "
                    "7371-6=7365, while Fraction elimination constructs 7365 pivots "
                    "and leaves the constant target nonzero. Thus the Q column rank is "
                    "exactly 7365 and no degree<=2 rational Nullstellensatz certificate "
                    "exists for this recorded finite system."
                ),
            },
            "torus_elimination": torus,
            "reversal_mu4_scan": mu4,
        },
        "checks": [
            {
                "name": "degree_two_fraction_no_certificate",
                "passed": degree_two_passed,
                "detail": {
                    "constraint_span_rank": len(basis_indices),
                    "column_count": len(columns),
                    "column_rank": exact["column_rank"],
                    "observed_syzygy_count": len(exact["dependent_column_indices"]),
                    "constant_target_in_span": exact["constant_target_in_span"],
                    "wall_seconds": exact["wall_seconds"],
                },
            },
            {
                "name": "thin_torus_elimination_completed",
                "passed": torus_passed,
                "detail": {
                    "free_torus_coordinate_count": torus["free_torus_coordinate_count"],
                    "cleared_numerator_Q_span_rank": torus[
                        "cleared_numerator_Q_span_rank"
                    ],
                },
            },
            {
                "name": "reversal_mu4_chart_exhausted",
                "passed": mu4_passed,
                "detail": {
                    "raw_assignment_count": mu4["raw_assignment_count"],
                    "radical_solution_count": mu4["radical_solution_count"],
                    "construction_passing_candidate_count": mu4[
                        "construction_passing_candidate_count"
                    ],
                    "surviving_candidate_count": mu4["surviving_candidate_count"],
                },
            },
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "PASS degree-two Fraction elimination "
        f"rank={exact['column_rank']} columns={len(columns)} "
        f"constant_in_span={exact['constant_target_in_span']}"
    )
    print(
        "PASS thin torus "
        f"free={torus['free_torus_coordinate_count']} "
        f"numerator_rank={torus['cleared_numerator_Q_span_rank']}"
    )
    print(
        "PASS reversal Q(i) mu4 "
        f"radical_solutions={mu4['radical_solution_count']} "
        f"construction_passes={mu4['construction_passing_candidate_count']}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
