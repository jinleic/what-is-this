"""Clean-room exact verifier for the branch-11111 Kac--Ward investigation.

This file deliberately imports no producer module.  It independently rebuilds
all finite-box polynomials, every degree-two Macaulay column, the thin torus,
and the finite reversal-conjugate Q(i) slice.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    OPPOSITE_DIRECTION,
    full_weight_finite_system,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/kac_ward/branch11111.json"
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
SparseInteger = dict[Monomial, int]
SparseRational = dict[Monomial, Fraction]


def q_record(value: Fraction | sp.Rational | int) -> dict[str, int]:
    value = sp.Rational(value)
    return {"numerator": int(value.p), "denominator": int(value.q)}


def q_value(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def monomial_key(monomial: Monomial) -> tuple[int, Monomial]:
    return (sum(monomial), monomial)


def sparse_poly(expression: sp.Expr, variables: list[sp.Symbol]) -> SparseInteger:
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return {
        tuple(monomial): int(coefficient)
        for monomial, coefficient in polynomial.terms()
        if coefficient
    }


def matrix_hash(columns: list[SparseInteger]) -> str:
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


def independent_constraint_indices(
    equations: list[sp.Expr], variables: list[sp.Symbol]
) -> list[int]:
    basis: dict[Monomial, SparseRational] = {}
    selected = []
    for index, expression in enumerate(equations):
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


def multiplier_exponents(variable_count: int) -> list[Monomial]:
    zero = (0,) * variable_count
    output = [zero]
    output.extend(
        tuple(1 if index == coordinate else 0 for index in range(variable_count))
        for coordinate in range(variable_count)
    )
    output.extend(
        tuple(int(first == index) + int(second == index) for index in range(variable_count))
        for first in range(variable_count)
        for second in range(first, variable_count)
    )
    return output


def macaulay_columns(
    basis_polynomials: list[SparseInteger], multipliers: list[Monomial]
) -> list[SparseInteger]:
    return [
        {
            tuple(a + b for a, b in zip(monomial, multiplier, strict=True)): coefficient
            for monomial, coefficient in polynomial.items()
        }
        for polynomial in basis_polynomials
        for multiplier in multipliers
    ]


def fraction_reduction(columns: list[SparseInteger], variable_count: int) -> dict:
    """Independent exact sparse Q reduction, including the constant target."""
    basis: dict[Monomial, SparseRational] = {}
    dependent = []
    pivots = []
    for index, source in enumerate(columns):
        vector: SparseRational = {
            monomial: Fraction(coefficient) for monomial, coefficient in source.items()
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
                pivots.append(pivot)
                break
            for monomial, value in existing.items():
                reduced = vector.get(monomial, Fraction(0)) - coefficient * value
                if reduced:
                    vector[monomial] = reduced
                else:
                    vector.pop(monomial, None)
        if not vector:
            dependent.append(index)

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
    pivot_hash = hashlib.sha256(
        "".join(",".join(map(str, monomial)) + ";" for monomial in pivots).encode()
    ).hexdigest()
    return {
        "rank": len(basis),
        "dependent": dependent,
        "target": target,
        "pivot_hash": pivot_hash,
    }


def column_metadata(
    index: int, basis_indices: list[int], multipliers: list[Monomial], labels: list[dict]
) -> dict:
    basis_position, multiplier_position = divmod(index, len(multipliers))
    basis_index = basis_indices[basis_position]
    return {
        "basis_constraint_index": basis_index,
        "basis_constraint_label": labels[basis_index],
        "multiplier_exponents": list(multipliers[multiplier_position]),
    }


def verify_sparse_relation(columns: list[SparseInteger], terms: list[dict]) -> None:
    total: SparseRational = {}
    for term in terms:
        coefficient = q_value(term["coefficient"])
        for monomial, value in columns[int(term["column_index"])].items():
            reduced = total.get(monomial, Fraction(0)) + coefficient * value
            if reduced:
                total[monomial] = reduced
            else:
                total.pop(monomial, None)
    assert not total


def thin_radical_relations(variables: list[sp.Symbol]) -> list[sp.Expr]:
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


def denominator_record(denominator: sp.Expr, variables: list[sp.Symbol]) -> dict:
    terms = sp.Poly(denominator, *variables, domain=sp.QQ).terms()
    assert len(terms) == 1
    exponents, coefficient = terms[0]
    return {"coefficient": q_record(coefficient), "exponents": list(exponents)}


def pair_text(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def reversal_pair(pair: tuple[int, int]) -> tuple[int, int]:
    return (OPPOSITE_DIRECTION[pair[1]], OPPOSITE_DIRECTION[pair[0]])


def root_of_unity_chart(system, variables: list[sp.Symbol]):
    symbol_by_pair = dict(system.all_symbols)
    pair_by_symbol = {symbol: pair for pair, symbol in system.all_symbols}
    tree = set(DIRECTION_GAUGE_TREE)
    representatives = []
    coordinate_rule = {}
    fixed_orbits = []
    visited = set()
    for pair in sorted(symbol_by_pair):
        if pair in visited:
            continue
        reverse = reversal_pair(pair)
        assert reverse in symbol_by_pair
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
        else:
            coordinate = len(representatives)
            representatives.append(representative)
            coordinate_rule[representative] = (coordinate, 1)
            coordinate_rule[partner] = (coordinate, -1)
    active_rules = [coordinate_rule[pair_by_symbol[variable]] for variable in variables]
    assert len(representatives) == 11
    return representatives, active_rules, [0] * len(variables), fixed_orbits


def radical_mu4_congruences(
    radical: list[sp.Expr],
    variables: list[sp.Symbol],
    active_rules,
    free_count: int,
) -> list[dict]:
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
        assert len(nonconstant) == 1
        monomial, coefficient = nonconstant[0]
        assert abs(coefficient) == abs(constant) == 1
        target_value = -constant // coefficient
        coefficients = [0] * free_count
        for variable_index, power in enumerate(monomial):
            rule = active_rules[variable_index]
            if rule is not None:
                coordinate, sign = rule
                coefficients[coordinate] = (coefficients[coordinate] + sign * power) % 4
        output.append(
            {
                "relation": str(relation),
                "coefficients_mod_4": coefficients,
                "target_mod_4": 2 if target_value == -1 else 0,
            }
        )
    return output


def mod4_unit_rref(congruences: list[dict], variable_count: int):
    rows = [
        [int(value) % 4 for value in relation["coefficients_mod_4"]]
        + [int(relation["target_mod_4"]) % 4]
        for relation in congruences
    ]
    pivots = []
    pivot_row = 0
    for column in range(variable_count):
        candidate = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column] % 2),
            None,
        )
        if candidate is None:
            continue
        rows[pivot_row], rows[candidate] = rows[candidate], rows[pivot_row]
        inverse = 1 if rows[pivot_row][column] == 1 else 3
        rows[pivot_row] = [value * inverse % 4 for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row:
                continue
            factor = rows[row][column]
            if factor:
                rows[row] = [
                    (value - factor * pivot_value) % 4
                    for value, pivot_value in zip(rows[row], rows[pivot_row], strict=True)
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivots


def evaluate_mu4(polynomial: SparseInteger, exponents: list[int]) -> tuple[int, int]:
    powers = [0, 0, 0, 0]
    for monomial, coefficient in polynomial.items():
        phase = sum(
            exponent * coordinate
            for exponent, coordinate in zip(monomial, exponents, strict=True)
        ) % 4
        powers[phase] += coefficient
    return powers[0] - powers[2], powers[1] - powers[3]


def mu4_value_record(value: tuple[int, int]) -> dict[str, int]:
    return {"real": value[0], "imaginary": value[1]}


def check_envelope() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["provenance"]["script"] == "experiments/e80_kw_branch11111.py"
    assert payload["provenance"]["no_float_decisions"] is True
    assert payload["data"]["branch"] == "11111"
    assert payload["data"]["status"] == "UNRESOLVED"
    assert payload["data"]["claim_tag"] == "[UNRESOLVED]"


def check_degree_two_fraction_certificate() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    data = payload["data"]
    record = data["degree_two_linear_algebra"]
    variables, equations, labels, _system = construction_data()
    assert data["variable_order"] == [str(variable) for variable in variables]
    basis_indices = independent_constraint_indices(equations, variables)
    assert basis_indices == data["basis_constraint_indices"]
    assert len(basis_indices) == data["exact_Q_constraint_span_rank"] == 21
    basis_polynomials = [sparse_poly(equations[index], variables) for index in basis_indices]
    multipliers = multiplier_exponents(len(variables))
    columns = macaulay_columns(basis_polynomials, multipliers)
    assert len(multipliers) == data["multiplier_count"] == 351
    assert len(columns) == record["column_rank"] + record["observed_syzygy_count"] == 7371
    assert matrix_hash(columns) == data["integer_macaulay_matrix_sha256"]

    reduction = fraction_reduction(columns, len(variables))
    assert reduction["rank"] == record["column_rank"] == 7365
    assert reduction["dependent"] == record["dependent_column_indices"]
    assert reduction["pivot_hash"] == record["pivot_monomial_sha256"]
    assert reduction["target"]
    assert record["constant_target_in_span"] is False
    stored_target = {
        tuple(term["monomial_exponents"]): q_value(term["coefficient"])
        for term in record["constant_residual"]
    }
    assert reduction["target"] == stored_target == {(0,) * len(variables): Fraction(1)}

    assert len(record["dependency_relations"]) == record["observed_syzygy_count"] == 6
    for relation in record["dependency_relations"]:
        assert relation["dependent_column_index"] in reduction["dependent"]
        verify_sparse_relation(columns, relation["terms"])
        for term in relation["terms"]:
            assert {
                key: term[key]
                for key in ("basis_constraint_index", "basis_constraint_label", "multiplier_exponents")
            } == column_metadata(
                int(term["column_index"]), basis_indices, multipliers, labels
            )


def check_torus_elimination() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    stored = payload["data"]["torus_elimination"]
    variables, equations, labels, _system = construction_data()
    free, substitutions = thin_torus_substitution(variables)
    radical = thin_radical_relations(variables)
    assert all(sp.cancel(relation.subs(substitutions)) == 0 for relation in radical)
    assert stored["input_coordinate_count"] == len(variables) == 25
    assert stored["radical_relation_count"] == len(radical) == 9
    assert stored["free_torus_coordinate_count"] == len(free) == 16
    assert stored["free_variable_order"] == [str(variable) for variable in free]
    assert stored["radical_binomials"] == [str(relation) for relation in radical]
    assert stored["radical_substitution_verified"] is True

    primary_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in PRIMARY_SHAPES
    ]
    full_indices = [
        index for index, label in enumerate(labels) if tuple(label["shape"]) in FULL_SHAPES
    ]
    assert len(primary_indices) == stored["primary_constraint_count"] == 27
    assert all(sp.cancel(equations[index].subs(substitutions)) == 0 for index in primary_indices)
    assert stored["primary_constraints_vanish_on_torus"] is True
    assert len(full_indices) == stored["three_dimensional_constraint_count"] == 15

    nonzero = []
    nonzero_indices = []
    for stored_row, index in zip(stored["three_dimensional_laurent_rows"], full_indices, strict=True):
        reduced = sp.cancel(equations[index].subs(substitutions))
        numerator, denominator = sp.fraction(reduced)
        numerator = sp.expand(numerator)
        terms = sparse_poly(numerator, free)
        assert stored_row["constraint_index"] == index
        assert stored_row["constraint_label"] == labels[index]
        assert stored_row["is_zero_on_torus"] == (not terms)
        assert stored_row["numerator_term_count"] == len(terms)
        assert stored_row["numerator_total_degree_max"] == max(map(sum, terms), default=0)
        assert stored_row["numerator_sha256"] == matrix_hash([terms])
        assert stored_row["denominator"] == denominator_record(denominator, free)
        if terms:
            nonzero.append(numerator)
            nonzero_indices.append(index)
    independent = independent_constraint_indices(nonzero, free)
    assert len(nonzero) == stored["nonzero_three_dimensional_constraint_count"] == 10
    assert len(independent) == stored["cleared_numerator_Q_span_rank"] == 6
    assert [nonzero_indices[index] for index in independent] == stored[
        "cleared_numerator_independent_constraint_indices"
    ]


def check_reversal_mu4_exhaustion() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    stored = payload["data"]["reversal_mu4_scan"]
    variables, equations, labels, system = construction_data()
    representatives, active_rules, fixed_exponents, fixed_orbits = root_of_unity_chart(
        system, variables
    )
    radical = thin_radical_relations(variables)
    congruences = radical_mu4_congruences(
        radical, variables, active_rules, len(representatives)
    )
    radical_sparse = [sparse_poly(relation, variables) for relation in radical]
    construction_sparse = [sparse_poly(equation, variables) for equation in equations]
    rows, pivots = mod4_unit_rref(congruences, len(representatives))
    free_columns = [column for column in range(len(representatives)) if column not in pivots]
    assert stored["reversal_free_orbit_count"] == len(representatives) == 11
    assert stored["reversal_free_orbit_representatives"] == [
        pair_text(pair) for pair in representatives
    ]
    assert stored["tree_forced_reversal_orbits"] == fixed_orbits
    assert stored["raw_assignment_count"] == 4 ** len(representatives) == 4_194_304
    assert stored["thin_radical_congruences_mod_4"] == congruences
    assert stored["mod_4_unit_pivot_columns"] == pivots == [1, 3, 6]
    assert stored["mod_4_free_columns"] == free_columns
    assert stored["mod_4_presolved_assignment_count"] == 4 ** len(free_columns) == 65_536
    assert stored["mod_4_unit_rref_sha256"] == hashlib.sha256(
        "\n".join(",".join(map(str, row)) for row in rows).encode()
    ).hexdigest()
    assert stored["forced_holonomy_consequences"] == {
        "oriented_plaquette_products": [-1] * 6,
        "sum_of_oriented_plaquette_products": -6,
        "reversal_conjugation_plane_real_parts": [-1, -1, -1],
    }

    samples = []
    digest = hashlib.sha256()
    solution_count = 0
    passing = []
    for free_values in itertools.product(range(4), repeat=len(free_columns)):
        free_exponents = [0] * len(representatives)
        for column, value in zip(free_columns, free_values, strict=True):
            free_exponents[column] = value
        for row, pivot in zip(rows, pivots, strict=False):
            free_exponents[pivot] = (
                row[-1]
                - sum(row[column] * free_exponents[column] for column in free_columns)
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
        assert all(
            evaluate_mu4(polynomial, active_exponents) == (0, 0)
            for polynomial in radical_sparse
        )
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
        if len(samples) < 16:
            samples.append(candidate)
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
            passing.append(candidate)

    assert solution_count == stored["radical_solution_count"] == 65_536
    assert samples == stored["radical_solution_samples"]
    assert digest.hexdigest() == stored["first_failure_sha256"]
    assert not passing
    assert stored["construction_passing_candidate_count"] == 0
    assert stored["holdout_3x3x3_records"] == []
    assert stored["surviving_candidate_count"] == 0


def run_check(name: str, function) -> None:
    function()
    print(f"PASS {name}")


def main() -> None:
    run_check("envelope", check_envelope)
    run_check("degree-two Fraction certificate", check_degree_two_fraction_certificate)
    run_check("thin-torus elimination", check_torus_elimination)
    run_check("reversal Q(i) mu4 exhaustion", check_reversal_mu4_exhaustion)
    print("PASS")


if __name__ == "__main__":
    main()
