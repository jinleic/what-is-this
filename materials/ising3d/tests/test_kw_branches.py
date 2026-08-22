"""Independent exact checks for the gauge-tree Kac--Ward branch ledger."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    DIRECTION_LABELS,
    full_weight_finite_system,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/kac_ward/full_family.json"
RESULT = ROOT / "results/kac_ward/branches.json"
CONSTRUCTION_SHAPES = ((3, 3, 2), (2, 2, 3))
CERTIFICATE_SHAPES = (
    (3, 3, 2),
    (2, 2, 3),
    (3, 2, 3),
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
    (4, 2, 2),
    (2, 4, 2),
    (2, 2, 4),
    (2, 3, 3),
)
ORDERS = (4, 6, 8)

SparsePolynomial = tuple[tuple[int, tuple[int, ...]], ...]


def pair_name(pair: tuple[int, int]) -> str:
    return f"{DIRECTION_LABELS[pair[0]]}->{DIRECTION_LABELS[pair[1]]}"


def branch_equations(system, branch: str) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    symbol_map = dict(system.all_symbols)
    substitutions = {
        symbol_map[pair]: int(bit)
        for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
    }
    equations = [sp.expand(equation.subs(substitutions)) for equation in system.equations]
    active = sorted(set().union(*(equation.free_symbols for equation in equations)), key=str)
    return active, equations


def to_sparse(expressions, variables: list[sp.Symbol]) -> list[SparsePolynomial]:
    output = []
    for expression in expressions:
        polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
        output.append(
            tuple(
                (int(coefficient), tuple(int(exponent) for exponent in monomial))
                for monomial, coefficient in polynomial.terms()
            )
        )
    return output


def eval_sparse(poly: SparsePolynomial, point: list[int], modulus: int) -> int:
    total = 0
    for coefficient, exponents in poly:
        term = coefficient % modulus
        for value, exponent in zip(point, exponents, strict=True):
            if exponent:
                term = term * pow(value, exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def eval_system(
    polynomials: list[SparsePolynomial], point: list[int], modulus: int
) -> list[int]:
    return [eval_sparse(poly, point, modulus) for poly in polynomials]


def derivative(poly: SparsePolynomial, coordinate: int) -> SparsePolynomial:
    terms = []
    for coefficient, exponents in poly:
        power = exponents[coordinate]
        if not power:
            continue
        reduced = list(exponents)
        reduced[coordinate] -= 1
        terms.append((coefficient * power, tuple(reduced)))
    return tuple(terms)


def jacobian_mod(
    polynomials: list[SparsePolynomial], point: list[int], prime: int
) -> list[list[int]]:
    return [
        [eval_sparse(derivative(poly, coordinate), point, prime) for coordinate in range(len(point))]
        for poly in polynomials
    ]


def rational_reconstruct_independent(
    residue: int, modulus: int, bound: int
) -> Fraction | None:
    residue %= modulus
    if residue == 0:
        return Fraction(0, 1)
    old_r, current_r = modulus, residue
    old_t, current_t = 0, 1
    while current_r > bound:
        quotient = old_r // current_r
        old_r, current_r = current_r, old_r - quotient * current_r
        old_t, current_t = current_t, old_t - quotient * current_t
    numerator, denominator = current_r, current_t
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        abs(numerator) > bound
        or not 1 <= denominator <= bound
        or math.gcd(numerator, denominator) != 1
        or math.gcd(denominator, modulus) != 1
        or (numerator - residue * denominator) % modulus
    ):
        return None
    return Fraction(numerator, denominator)


def valuation_mod(residue: int, prime: int, exponent: int) -> int:
    residue %= prime**exponent
    if residue == 0:
        return exponent
    valuation = 0
    while residue % prime == 0:
        residue //= prime
        valuation += 1
    return valuation


def load_payloads() -> tuple[dict, dict]:
    return (
        json.loads(SOURCE.read_text(encoding="utf-8")),
        json.loads(RESULT.read_text(encoding="utf-8")),
    )


def check_envelope_and_branch_decomposition() -> None:
    source, payload = load_payloads()
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e50_kw_branches.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert payload["provenance"]["groebner_timeout_seconds_per_branch"] <= 900
    assert all(check["passed"] for check in payload["checks"])

    expected_labels = [format(index, "05b") for index in range(32)]
    source_labels = [row["branch"] for row in source["data"]["zero_pattern_branches"]]
    records = payload["data"]["branches"]
    assert payload["data"]["total_branch_count"] == 32
    assert expected_labels == source_labels == [row["branch"] for row in records]
    assert payload["data"]["status_counts"] == {
        "EMPTY_OVER_Q": 17,
        "RATIONAL_POINT_FOUND": 0,
        "UNRESOLVED": 15,
    }
    assert payload["data"]["global_claim"].startswith("[UNRESOLVED]")

    construction = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=False
    )
    remaining = [
        symbol
        for pair, symbol in construction.all_symbols
        if pair not in DIRECTION_GAUGE_TREE
    ]
    gauge_names = [pair_name(pair) for pair in DIRECTION_GAUGE_TREE]
    source_by_branch = {
        row["branch"]: row for row in source["data"]["zero_pattern_branches"]
    }
    for record in records:
        branch = record["branch"]
        active, equations = branch_equations(construction, branch)
        linear = [
            variable
            for variable in remaining
            if max(sp.degree(equation, variable) for equation in equations) <= 1
            and any(sp.degree(equation, variable) == 1 for equation in equations)
        ]
        zeros = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "0"
        ]
        ones = [
            name for name, bit in zip(gauge_names, branch, strict=True) if bit == "1"
        ]
        assert record["zero_pattern"] == zeros
        assert record["normalized_nonzero"] == ones
        assert record["active_variable_count_in_six_equations"] == len(active)
        assert record["globally_linear_variable_count"] == len(linear)
        assert source_by_branch[branch]["zero_pattern"] == zeros
        assert source_by_branch[branch]["normalized_nonzero"] == ones


def check_all_empty_certificates_fresh() -> None:
    _source, payload = load_payloads()
    certificate_system = full_weight_finite_system(
        CERTIFICATE_SHAPES, ORDERS, gauge_fix=False
    )
    expected_labels = [
        {"shape": list(shape), "order": order}
        for shape, order in certificate_system.labels
    ]
    empty_records = [
        row for row in payload["data"]["branches"] if row["status"] == "EMPTY_OVER_Q"
    ]
    assert len(empty_records) == 17
    for record in empty_records:
        _active, equations = branch_equations(certificate_system, record["branch"])
        certificate = record["empty_certificate"]
        assert certificate["type"] in {
            "two_equation_constant_identity",
            "single_equation_nonzero_constant",
            "full_Q_linear_span_constant_identity",
        }
        coefficients = certificate["coefficient_vector"]
        assert len(coefficients) == len(equations) == len(expected_labels)
        combination = sp.expand(
            sum(
                coefficient * equation
                for coefficient, equation in zip(coefficients, equations, strict=True)
                if coefficient
            )
        )
        assert not combination.free_symbols
        assert combination == certificate["constant"]
        assert certificate["constant"] != 0
        stored_terms = [
            {
                "coefficient": coefficient,
                "shape": expected_labels[index]["shape"],
                "order": expected_labels[index]["order"],
            }
            for index, coefficient in enumerate(coefficients)
            if coefficient
        ]
        assert certificate["terms"] == stored_terms

    unresolved = [
        row for row in payload["data"]["branches"] if row["status"] == "UNRESOLVED"
    ]
    assert len(unresolved) == 15
    assert all(row["groebner"]["status"] == "timeout" for row in unresolved)
    assert all(row["groebner"]["field"] == "Q" for row in unresolved)
    assert all(row["groebner"]["limit_seconds"] == 60 for row in unresolved)


def check_padic_lifts_reconstruction_and_falsifications() -> None:
    source, payload = load_payloads()
    construction = full_weight_finite_system(
        CONSTRUCTION_SHAPES, ORDERS, gauge_fix=True
    )
    variables = list(construction.variables)
    construction_polynomials = to_sparse(construction.equations, variables)
    holdout = full_weight_finite_system(((3, 3, 3),), ORDERS, gauge_fix=True)
    assert list(map(str, holdout.variables)) == list(map(str, variables))
    holdout_polynomials = to_sparse(holdout.equations, variables)
    witnesses = {
        int(row["prime"]): [int(value) for value in row["values_in_variable_order"]]
        for row in source["data"]["exact_modular_solution_witnesses"]
    }
    records = {
        int(row["prime"]): row
        for row in payload["data"]["p_adic_records_for_branch_11111"]
    }
    assert set(records) == set(witnesses) == {3, 5, 7, 11}
    for prime, point in witnesses.items():
        assert eval_system(construction_polynomials, point, prime) == [0] * 6
        assert records[prime]["stored_solution_verified"]

    assert {prime for prime, row in records.items() if row["status"] == "lifted"} == {
        3,
        5,
        11,
    }
    for prime in (3, 5, 11):
        record = records[prime]
        exponent = int(record["precision_exponent"])
        modulus = int(record["modulus"])
        residues = [int(value) for value in record["residues"]]
        assert exponent >= 60
        assert modulus == prime**exponent
        assert eval_system(construction_polynomials, residues, modulus) == [0] * 6

        reconstruction = record["rational_reconstruction"]
        bound = int(reconstruction["height_bound"])
        assert 2 * bound * bound < modulus
        assert 2 * (bound + 1) * (bound + 1) >= modulus
        independently_reconstructed = []
        for coordinate, residue in zip(
            reconstruction["coordinates"], residues, strict=True
        ):
            observed = rational_reconstruct_independent(residue, modulus, bound)
            independently_reconstructed.append(observed)
            if coordinate["status"] == "failure":
                assert observed is None
            else:
                assert observed == Fraction(
                    int(coordinate["numerator"]), int(coordinate["denominator"])
                )
        assert len(independently_reconstructed) == 25
        assert sum(value is not None for value in independently_reconstructed) == reconstruction[
            "successful_coordinate_count"
        ]
        assert not reconstruction["all_coordinates_reconstructed"]

        observed_holdout = eval_system(holdout_polynomials, residues, modulus)
        stored_holdout = record["holdout_3x3x3"]
        assert observed_holdout == [int(value) for value in stored_holdout["residues"]]
        assert any(observed_holdout)
        assert not stored_holdout["passed_to_full_precision"]
        assert [
            valuation_mod(value, prime, exponent) for value in observed_holdout
        ] == stored_holdout["p_adic_valuations_capped_at_precision"]

    seven = records[7]
    assert seven["status"] == "no_lift_mod_p_squared"
    point7 = witnesses[7]
    residuals49 = eval_system(construction_polynomials, point7, 49)
    right_hand_side7 = [-(value // 7) % 7 for value in residuals49]
    jacobian7 = jacobian_mod(construction_polynomials, point7, 7)
    certificate7 = seven["certificate"]
    left7 = certificate7["left_kernel_vector"]
    assert all(
        sum(left7[row] * jacobian7[row][column] for row in range(6)) % 7 == 0
        for column in range(25)
    )
    assert sum(left7[row] * right_hand_side7[row] for row in range(6)) % 7 == certificate7[
        "pairing_with_rhs_mod_p"
    ] != 0

    combined = construction_polynomials + holdout_polynomials
    point3 = witnesses[3]
    assert eval_system(combined, point3, 3) == [0] * 9
    residuals9 = eval_system(combined, point3, 9)
    right_hand_side3 = [-(value // 3) % 3 for value in residuals9]
    jacobian3 = jacobian_mod(combined, point3, 3)
    combined_record = payload["data"]["p3_combined_3x3x3_nonlift"]
    assert combined_record["status"] == "no_lift_mod_p_squared"
    certificate3 = combined_record["certificate"]
    left3 = certificate3["left_kernel_vector"]
    assert all(
        sum(left3[row] * jacobian3[row][column] for row in range(9)) % 3 == 0
        for column in range(25)
    )
    assert sum(left3[row] * right_hand_side3[row] for row in range(9)) % 3 == certificate3[
        "pairing_with_rhs_mod_p"
    ] != 0
    assert payload["data"]["exact_rational_reconstruction_candidates"] == []


def run_check(name: str, function) -> None:
    try:
        function()
    except Exception as error:
        print(f"FAIL {name}: {type(error).__name__}: {error}")
        raise
    print(f"PASS {name}")


def main() -> None:
    run_check("result envelope and all 32 branch labels", check_envelope_and_branch_decomposition)
    run_check("all EMPTY_OVER_Q constant certificates", check_all_empty_certificates_fresh)
    run_check(
        "p-adic lifts, reconstructions, and exact holdout falsifications",
        check_padic_lifts_reconstruction_and_falsifications,
    )
    print("PASS")


if __name__ == "__main__":
    main()
