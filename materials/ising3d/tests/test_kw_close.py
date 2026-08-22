"""Independent exact checks for the wave-8 Kac--Ward branch closure artifact."""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from ising.exact_enumeration import even_subgraph_polynomial
from ising.fermions.kac_ward import (
    DIRECTION_GAUGE_TREE,
    formal_log_coefficients,
    full_weight_finite_system,
    translation_invariant_trace_expression,
)
from ising.lattices import cubic

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "results/kac_ward/branches.json"
RESULT = ROOT / "results/kac_ward/close.json"
BASE_ORDERS = (4, 6, 8)
SCAN_ORDERS = (4, 6, 8, 10)
SCAN_CANONICAL = (
    (2, 2, 1),
    (3, 2, 1),
    (4, 2, 1),
    (3, 3, 1),
    (4, 3, 1),
    (4, 4, 1),
)
SCAN_SHAPES = tuple(
    sorted(set().union(*(set(itertools.permutations(shape)) for shape in SCAN_CANONICAL)))
)
PRIMARY_SHAPES = tuple(
    sorted(
        set(itertools.permutations((2, 2, 1)))
        | set(itertools.permutations((3, 2, 1)))
    )
)
MACAULAY_SHAPES = PRIMARY_SHAPES + (
    (3, 3, 2),
    (2, 2, 3),
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
)

SparsePolynomial = tuple[tuple[int, tuple[int, ...]], ...]


def load_payloads() -> tuple[dict, dict]:
    return (
        json.loads(PRIOR.read_text(encoding="utf-8")),
        json.loads(RESULT.read_text(encoding="utf-8")),
    )


def branch_equations(system, branch: str) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    symbol_map = dict(system.all_symbols)
    substitutions = {
        symbol_map[pair]: int(bit)
        for pair, bit in zip(DIRECTION_GAUGE_TREE, branch, strict=True)
    }
    equations = [sp.expand(equation.subs(substitutions)) for equation in system.equations]
    variables = sorted(set().union(*(equation.free_symbols for equation in equations)), key=str)
    return variables, equations


def q_value(record: dict) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def constraint_id(shape, order: int) -> str:
    return "x".join(map(str, shape)) + f":k{order}"


def to_sparse(expressions: list[sp.Expr], variables: list[sp.Symbol]) -> list[SparsePolynomial]:
    return [
        tuple(
            (int(coefficient), tuple(int(exponent) for exponent in monomial))
            for monomial, coefficient in sp.Poly(expression, *variables, domain=sp.ZZ).terms()
        )
        for expression in expressions
    ]


def eval_sparse(poly: SparsePolynomial, point: list[int], modulus: int) -> int:
    total = 0
    for coefficient, exponents in poly:
        term = coefficient % modulus
        for value, exponent in zip(point, exponents, strict=True):
            if exponent:
                term = term * pow(value, exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def eval_system(polynomials: list[SparsePolynomial], point: list[int], modulus: int) -> list[int]:
    return [eval_sparse(polynomial, point, modulus) for polynomial in polynomials]


def derivative(poly: SparsePolynomial, coordinate: int) -> SparsePolynomial:
    output = []
    for coefficient, exponents in poly:
        exponent = exponents[coordinate]
        if not exponent:
            continue
        reduced = list(exponents)
        reduced[coordinate] -= 1
        output.append((coefficient * exponent, tuple(reduced)))
    return tuple(output)


def jacobian(polynomials: list[SparsePolynomial], point: list[int], prime: int) -> list[list[int]]:
    return [
        [eval_sparse(derivative(polynomial, coordinate), point, prime) for coordinate in range(len(point))]
        for polynomial in polynomials
    ]


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    rows = [[value % prime for value in row] for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (left - factor * right) % prime
                for left, right in zip(rows[row], rows[rank], strict=True)
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


@lru_cache(maxsize=1)
def scan_system():
    return full_weight_finite_system(SCAN_SHAPES, SCAN_ORDERS, gauge_fix=False)


@lru_cache(maxsize=1)
def primary_system():
    return full_weight_finite_system(PRIMARY_SHAPES, BASE_ORDERS, gauge_fix=False)


@lru_cache(maxsize=1)
def macaulay_system():
    return full_weight_finite_system(MACAULAY_SHAPES, BASE_ORDERS, gauge_fix=False)


@lru_cache(maxsize=1)
def holdout_system():
    return full_weight_finite_system(((3, 3, 3),), BASE_ORDERS, gauge_fix=False)


def independent_indices(expressions: list[sp.Expr], variables: list[sp.Symbol]) -> list[int]:
    polynomials = [sp.Poly(expression, *variables, domain=sp.QQ) for expression in expressions]
    monomials = sorted(set().union(*(set(polynomial.monoms()) for polynomial in polynomials)))
    matrix = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for monomial in monomials]
            for polynomial in polynomials
        ]
    )
    return list(matrix.T.rref()[1])


def build_macaulay_columns(
    polynomials: list[sp.Poly], variable_count: int
) -> list[dict[tuple[int, ...], int]]:
    zero = (0,) * variable_count
    multipliers = [zero] + [
        tuple(1 if index == coordinate else 0 for index in range(variable_count))
        for coordinate in range(variable_count)
    ]
    columns = []
    for polynomial in polynomials:
        terms = {monomial: int(coefficient) for monomial, coefficient in polynomial.terms()}
        for multiplier in multipliers:
            columns.append(
                {
                    tuple(a + b for a, b in zip(monomial, multiplier, strict=True)): coefficient
                    for monomial, coefficient in terms.items()
                }
            )
    return columns


def matrix_hash(columns: list[dict[tuple[int, ...], int]]) -> str:
    digest = hashlib.sha256()
    for index, column in enumerate(columns):
        digest.update(f"C{index}:".encode())
        for monomial, coefficient in sorted(column.items()):
            digest.update((f"{coefficient}@{','.join(map(str, monomial))};").encode())
    return digest.hexdigest()


def column_rank_and_constant(columns, variable_count: int, prime: int) -> tuple[int, int]:
    basis = {}
    for source in columns:
        vector = {m: c % prime for m, c in source.items() if c % prime}
        while vector:
            pivot = max(vector, key=lambda monomial: (sum(monomial), monomial))
            value = vector[pivot]
            if pivot not in basis:
                inverse = pow(value, -1, prime)
                vector = {m: c * inverse % prime for m, c in vector.items() if c % prime}
                basis[pivot] = vector
                break
            for monomial, coefficient in basis[pivot].items():
                reduced = (vector.get(monomial, 0) - value * coefficient) % prime
                if reduced:
                    vector[monomial] = reduced
                elif monomial in vector:
                    del vector[monomial]
    target = {(0,) * variable_count: 1}
    while target:
        pivot = max(target, key=lambda monomial: (sum(monomial), monomial))
        value = target[pivot]
        if pivot not in basis:
            return len(basis), len(basis) + 1
        for monomial, coefficient in basis[pivot].items():
            reduced = (target.get(monomial, 0) - value * coefficient) % prime
            if reduced:
                target[monomial] = reduced
            elif monomial in target:
                del target[monomial]
    return len(basis), len(basis)


def check_envelope_and_status_map() -> None:
    prior, payload = load_payloads()
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e62_kw_close.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert payload["provenance"]["no_float_decisions"] is True
    assert all(check["passed"] for check in payload["checks"])
    records = payload["data"]["branches"]
    assert [record["branch"] for record in records] == [format(index, "05b") for index in range(32)]
    assert payload["data"]["status_counts"] == {
        "EMPTY_OVER_Q": 31,
        "RATIONAL_POINT_FOUND": 0,
        "UNRESOLVED": 1,
    }
    assert [record["branch"] for record in records if record["status"] == "UNRESOLVED"] == [
        "11111"
    ]
    assert payload["data"]["global_claim"].startswith("[UNRESOLVED]")

    prior_empty = {
        record["branch"]: record["empty_certificate"]
        for record in prior["data"]["branches"]
        if record["status"] == "EMPTY_OVER_Q"
    }
    imported = {
        record["branch"]: record["empty_certificate"]
        for record in records
        if record.get("certificate_origin") == "e50_imported_unchanged"
    }
    assert imported == prior_empty
    assert len(imported) == 17


def check_new_constant_certificates_from_raw_lattice_data() -> None:
    _prior, payload = load_payloads()
    raw_records = {
        record["constraint_id"]: record
        for record in payload["data"]["raw_lattice_constants_for_new_certificates"]
    }
    assert raw_records
    raw_equations = {}
    symbol_system = full_weight_finite_system(((2, 2, 1),), (4,), gauge_fix=False)
    symbol_map = dict(symbol_system.all_symbols)
    for identifier, record in raw_records.items():
        shape = tuple(record["shape"])
        order = int(record["order"])
        polynomial = even_subgraph_polynomial(cubic(*shape, periodic=False))
        assert polynomial[: order + 1] == record["even_subgraph_coefficients_through_order"]
        observed_log = formal_log_coefficients(polynomial, order)[order]
        assert observed_log == q_value(record["log_even_subgraph_coefficient"])
        required = -2 * order * observed_log
        assert required == q_value(record["required_trace"])
        assert identifier == constraint_id(shape, order)
        raw_equations[identifier] = sp.expand(
            translation_invariant_trace_expression(shape, order, gauge_fix=False)
            - sp.Rational(required.numerator, required.denominator)
        )

    new_records = [
        record
        for record in payload["data"]["branches"]
        if record.get("certificate_origin") == "e62_new_thin_box_identity"
    ]
    assert len(new_records) == 14
    expected = {
        "01010",
        "01011",
        "01100",
        "01101",
        "01110",
        "01111",
        "10011",
        "10101",
        "10111",
        "11010",
        "11011",
        "11100",
        "11101",
        "11110",
    }
    assert {record["branch"] for record in new_records} == expected
    for record in new_records:
        substitutions = {
            symbol_map[pair]: int(bit)
            for pair, bit in zip(DIRECTION_GAUGE_TREE, record["branch"], strict=True)
        }
        certificate = record["empty_certificate"]
        combination = sp.expand(
            sum(
                int(term["coefficient"])
                * raw_equations[term["constraint_id"]].subs(substitutions)
                for term in certificate["terms"]
            )
        )
        assert not combination.free_symbols
        assert combination == certificate["constant"] != 0
        assert certificate["identity_verified"] is True


def check_systematic_span_and_primary_radical() -> None:
    _prior, payload = load_payloads()
    generic = next(record for record in payload["data"]["branches"] if record["branch"] == "11111")
    variables, equations = branch_equations(scan_system(), "11111")
    stored = generic["degree_zero_thin_span"]
    polynomials = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
    zero = (0,) * len(variables)
    nonconstant = sorted(
        set().union(*(set(polynomial.monoms()) - {zero} for polynomial in polynomials))
    )
    all_monomials = sorted(set(nonconstant) | {zero})
    nonconstant_matrix = sp.Matrix(
        [[polynomial.coeff_monomial(m) for polynomial in polynomials] for m in nonconstant]
    )
    full_matrix = sp.Matrix(
        [[polynomial.coeff_monomial(m) for polynomial in polynomials] for m in all_monomials]
    )
    assert nonconstant_matrix.rank() == stored["nonconstant_coefficient_rank"]
    assert full_matrix.rank() == stored["full_coefficient_rank"]
    assert stored["full_coefficient_rank"] == stored["nonconstant_coefficient_rank"]
    assert stored["constant_in_equation_span"] is False

    pvariables, pequations = branch_equations(primary_system(), "11111")
    namespace = {str(variable): variable for variable in pvariables}
    labels = [
        {"shape": list(shape), "order": order} for shape, order in primary_system().labels
    ]
    by_id = {
        constraint_id(label["shape"], label["order"]): equation
        for label, equation in zip(labels, pequations, strict=True)
    }
    primary = generic["thin_primary_radical"]
    for plane in primary["planes"]:
        for relation_name in ("sum_relation", "product_relation", "straight_relation"):
            relation = plane[relation_name]
            target = sp.expand(sp.sympify(relation["polynomial"], locals=namespace))
            observed = sp.expand(
                sum(
                    sp.Rational(term["coefficient"]) * by_id[term["constraint_id"]]
                    for term in relation["linear_combination_of_raw_constraints"]
                )
            )
            assert observed == target
            assert relation["identity_verified"] is True
    radical = [sp.Poly(sp.sympify(text, locals=namespace), *pvariables) for text in primary["radical_binomials"]]
    radical_basis = sp.groebner(
        [polynomial.as_expr() for polynomial in radical],
        *pvariables,
        order="grevlex",
        domain=sp.QQ,
    )
    assert all(
        sp.expand(radical_basis.reduce(equation)[1]) == 0 for equation in pequations
    )
    assert primary["all_raw_constraints_reduce_to_zero_mod_radical"] is True
    exponent_rows = []
    for polynomial in radical:
        terms = [monomial for monomial, _ in polynomial.terms() if any(monomial)]
        assert len(terms) == 1
        exponent_rows.append(terms[0])
    diagonal = [
        abs(int(value))
        for value in smith_normal_form(sp.Matrix(exponent_rows), domain=ZZ).diagonal()
        if value
    ]
    assert diagonal == primary["exponent_lattice_smith_diagonal"] == [1] * 9
    assert primary["prime_radical_dimension"] == len(pvariables) - 9 == 16


def check_bounded_nullstellensatz_matrix() -> None:
    _prior, payload = load_payloads()
    generic = next(record for record in payload["data"]["branches"] if record["branch"] == "11111")
    record = generic["bounded_nullstellensatz"]
    variables, equations = branch_equations(macaulay_system(), "11111")
    indices = independent_indices(equations, variables)
    assert indices == record["basis_constraint_indices"]
    assert len(indices) == record["exact_Q_constraint_span_rank"] == 21
    polynomials = [sp.Poly(equations[index], *variables, domain=sp.ZZ) for index in indices]
    columns = build_macaulay_columns(polynomials, len(variables))
    assert len(columns) == record["column_count"] == 546
    assert matrix_hash(columns) == record["integer_macaulay_matrix_sha256"]
    for stored in record["modular_same_matrix_rank_certificates"]:
        rank, augmented = column_rank_and_constant(columns, len(variables), int(stored["prime"]))
        assert rank == stored["column_rank"] == 546
        assert augmented == stored["augmented_rank"] == 547
        assert stored["constant_target_independent"] is True
    assert record["certificate_absent_at_budget"] is True


def check_modular_points_padic_certificates_and_holdout() -> None:
    _prior, payload = load_payloads()
    generic = next(record for record in payload["data"]["branches"] if record["branch"] == "11111")
    variables, equations = branch_equations(macaulay_system(), "11111")
    polynomials = to_sparse(equations, variables)
    hvariables, hequations = branch_equations(holdout_system(), "11111")
    assert [str(variable) for variable in variables] == [str(variable) for variable in hvariables]
    holdout_polynomials = to_sparse(hequations, hvariables)
    for record in generic["modular_and_padic_diagnostics"]:
        prime = int(record["prime"])
        initial = [int(value) for value in record["values_in_variable_order"]]
        assert eval_system(polynomials, initial, prime) == record["construction_residues_mod_p"]
        assert not any(record["construction_residues_mod_p"])
        matrix = jacobian(polynomials, initial, prime)
        assert rank_mod(matrix, prime) == record["jacobian_rank_mod_p"]

        lift = record["selected_zero_free_padic_path"]
        point = [int(value) for value in lift["residues"]]
        exponent = int(lift["last_solved_exponent"])
        modulus = prime**exponent
        assert modulus == int(lift["modulus"])
        assert not any(eval_system(polynomials, point, modulus))
        residuals = eval_system(polynomials, point, modulus * prime)
        right = [-(residual // modulus) % prime for residual in residuals]
        certificate = lift["certificate"]
        left = [int(value) for value in certificate["left_kernel_vector"]]
        assert all(
            sum(left[row] * matrix[row][column] for row in range(len(matrix))) % prime == 0
            for column in range(len(variables))
        )
        assert sum(left[row] * right[row] for row in range(len(matrix))) % prime == certificate[
            "pairing_with_rhs_mod_p"
        ] != 0
        assert right == certificate["linearized_rhs_mod_p"]
        assert record["holdout_3x3x3_residues_at_last_precision"] == eval_system(
            holdout_polynomials, point, modulus
        )

        substitution = dict(zip(variables, initial, strict=True))
        exact_construction = [
            int(sp.Poly(equation, *variables).eval(substitution)) for equation in equations
        ]
        holdout_substitution = dict(zip(hvariables, initial, strict=True))
        exact_holdout = [
            int(sp.Poly(equation, *hvariables).eval(holdout_substitution))
            for equation in hequations
        ]
        assert exact_construction == record["least_residue_integer_candidate"][
            "construction_residuals"
        ]
        assert exact_holdout == record["least_residue_integer_candidate"]["holdout_residuals"]
        assert not record["least_residue_integer_candidate"]["passed_all_exact_equations"]
    assert generic["exact_rational_candidate_count"] == 0


def run_check(name: str, function) -> None:
    try:
        function()
    except Exception as error:
        print(f"FAIL {name}: {type(error).__name__}: {error}")
        raise
    print(f"PASS {name}")


def main() -> None:
    run_check("envelope, imported map, and honest survivor", check_envelope_and_status_map)
    run_check(
        "fourteen constant certificates from raw lattice constants",
        check_new_constant_certificates_from_raw_lattice_data,
    )
    run_check("systematic span and thin primary radical", check_systematic_span_and_primary_radical)
    run_check("bounded rational Nullstellensatz matrix", check_bounded_nullstellensatz_matrix)
    run_check(
        "modular points, selected p-adic certificates, and 3x3x3 holdout",
        check_modular_points_padic_certificates_and_holdout,
    )
    print("PASS")


if __name__ == "__main__":
    main()
