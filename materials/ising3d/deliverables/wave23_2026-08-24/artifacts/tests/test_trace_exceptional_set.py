#!/usr/bin/env python3
"""Independent exact verifier for the symbolic trace exceptional set.

This file deliberately does not import the producer.  Its elimination path is a
degree-filtered affine Macaulay reduction: at each total degree it solves the
constant leading-form coefficient matrix directly.  The producer instead uses
lifted grevlex H-basis relations.  Its independent root path rebuilds the
integer binomial shift E(1+x), applies Descartes' rule, and checks three
disjoint rational sign brackets.  This replaces an exact Sturm call whose
resource cost was observed to exceed the verifier wall.

Set ``TRACE_EXCEPTIONAL_STAGE=elimination`` or ``roots`` to run either bounded
half independently; the default runs both.
"""

from __future__ import annotations

import hashlib
import math
import os
import json
import sys
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.polys.domains import QQ, ZZ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
PRODUCER = ROOT / "experiments" / "e238_trace_exceptional_set.py"
N_SITES = 6
BONDS_2X3 = ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5))
BONDS_CHAIN = tuple((index, index + 1) for index in range(N_SITES - 1))

QPAR, U = sp.symbols("q u")
A5, A4, A3 = sp.symbols("a5 a4 a3")
VARS = (A5, A4, A3)
QQQ = QQ.frac_field(QPAR)
QQ_POLY_Q = QQ.poly_ring(QPAR)
ZZ_POLY_Q = ZZ.poly_ring(QPAR)
ZERO = QQQ.zero
ONE = QQQ.one
BASIS_MONOMIALS = (
    (0, 0, 0),
    (0, 0, 1),
    (0, 0, 2),
    (0, 0, 3),
    (1, 0, 0),
    (1, 0, 1),
    (0, 1, 0),
    (0, 1, 1),
)
BASIS_INDEX = {monomial: index for index, monomial in enumerate(BASIS_MONOMIALS)}
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_rational(text: str) -> sp.Rational:
    numerator, denominator = text.split("/")
    return sp.Rational(int(numerator), int(denominator))


def coefficients_ascending(poly: sp.Poly) -> list[str]:
    return [str(coefficient) for coefficient in reversed(poly.all_coeffs())]


def monomial_expression(monomial: tuple[int, int, int]) -> sp.Expr:
    return A5 ** monomial[0] * A4 ** monomial[1] * A3 ** monomial[2]


def add_monomials(
    left: tuple[int, int, int], right: tuple[int, int, int]
) -> tuple[int, int, int]:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def monomials_of_degree(degree: int) -> list[tuple[int, int, int]]:
    return [
        (first, second, degree - first - second)
        for first in range(degree, -1, -1)
        for second in range(degree - first, -1, -1)
    ]


def spin_disagreements(state: int, bonds: tuple[tuple[int, int], ...]) -> int:
    def bit(site: int) -> int:
        return (state >> (N_SITES - 1 - site)) & 1

    return sum(bit(left) != bit(right) for left, right in bonds)


def traces_over_zzq(
    bonds: tuple[tuple[int, int], ...], maximum: int = 7
) -> tuple[int, int, list[sp.Expr]]:
    epsilon = len(bonds) % 2
    d_values: list[int] = []
    for state in range(1 << N_SITES):
        bond_value = len(bonds) - 2 * spin_disagreements(state, bonds)
        d_values.append((bond_value - epsilon) // 2)
    shift = -min(d_values)
    entries = [
        [
            QPAR
            ** (
                d_values[row]
                + N_SITES
                - (row ^ column).bit_count()
                + shift
            )
            for column in range(1 << N_SITES)
        ]
        for row in range(1 << N_SITES)
    ]
    matrix = DomainMatrix.from_list_sympy(64, 64, entries).convert_to(ZZ_POLY_Q)
    power = DomainMatrix.eye(64, ZZ_POLY_Q).to_dense()
    traces: list[sp.Expr] = []
    for _ in range(maximum):
        power = power.matmul(matrix)
        trace = ZZ_POLY_Q.zero
        for index in range(64):
            trace += power[index, index].element
        traces.append(trace.as_expr())
    return epsilon, shift, traces


def normalized_targets(
    traces: list[sp.Expr], epsilon: int, shift: int
) -> dict[str, sp.Expr]:
    delta = QPAR**2 - 1
    one = QPAR ** (epsilon - 2 * shift) / delta**6
    two = QPAR ** (2 * epsilon - 4 * shift) / delta**12
    three = QPAR ** (3 * epsilon - 6 * shift) / delta**18
    return {
        "r1_squared": sp.cancel(traces[0] ** 2 * one),
        "r2": sp.cancel(traces[1] * one),
        "r3_over_r1": sp.cancel(traces[2] * one / traces[0]),
        "r4": sp.cancel(traces[3] * two),
        "r5_over_r1": sp.cancel(traces[4] * two / traces[0]),
        "r6_over_r2": sp.cancel(traces[5] * two / traces[1]),
        "r7_over_r1": sp.cancel(traces[6] * three / traces[0]),
    }


def norm_equations(
    targets: dict[str, sp.Expr],
) -> tuple[list[sp.Poly], sp.Expr, sp.Expr]:
    a0 = targets["r1_squared"]
    c2 = targets["r2"] - a0 - 64 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = targets["r3_over_r1"] - a0 - 729 - 243 * A5 - 81 * A4 - 27 * A3
    a2 = sp.cancel(c3 / 3 - c2 / 2)
    a1 = sp.cancel(3 * c2 / 2 - 2 * c3 / 3)
    sextic = sp.expand(
        U**6 + A5 * U**5 + A4 * U**4 + A3 * U**3 + a2 * U**2 + a1 * U + a0
    )
    minimals = (
        ("r4", U**2 - 4 * U + 2),
        ("r5_over_r1", U**2 - 5 * U + 5),
        ("r6_over_r2", U**2 - 4 * U + 1),
        ("r7_over_r1", U**3 - 7 * U**2 + 14 * U - 7),
    )
    equations = [
        sp.Poly(
            sp.cancel(sp.resultant(minimal, sextic, U) - targets[name]),
            *VARS,
            domain=QQQ,
        )
        for name, minimal in minimals
    ]
    return equations, a2, a1


def leading_quadrics(equations: list[sp.Poly]) -> list[sp.Poly]:
    result: list[sp.Poly] = []
    for equation in equations[:3]:
        expression = sp.Integer(0)
        for monomial, coefficient in equation.rep.to_dict().items():
            if sum(monomial) == 2:
                expression += coefficient.as_expr() * monomial_expression(monomial)
        result.append(sp.Poly(expression, *VARS, domain=QQ))
    return result


def macaulay_plans(
    homogeneous: list[sp.Poly], maximum_degree: int
) -> dict[int, dict[str, object]]:
    """Independent graded coefficient solves, with no Groebner lifting."""

    plans: dict[int, dict[str, object]] = {}
    for degree in range(maximum_degree + 1):
        monomials = monomials_of_degree(degree)
        columns: list[sp.Matrix] = []
        labels: list[tuple[int, tuple[int, int, int]]] = []
        if degree >= 2:
            for equation_index, equation in enumerate(homogeneous):
                for multiplier in monomials_of_degree(degree - 2):
                    product = sp.Poly(
                        equation.as_expr() * monomial_expression(multiplier),
                        *VARS,
                        domain=QQ,
                    ).rep.to_dict()
                    columns.append(
                        sp.Matrix(
                            [product.get(monomial, sp.Rational(0)) for monomial in monomials]
                        )
                    )
                    labels.append((equation_index, multiplier))
        ideal_matrix = (
            sp.Matrix.hstack(*columns) if columns else sp.zeros(len(monomials), 0)
        )
        _, pivots = ideal_matrix.rref()
        selected_labels = [labels[index] for index in pivots]
        selected_columns = [columns[index] for index in pivots]
        complement = [
            monomial for monomial in BASIS_MONOMIALS if sum(monomial) == degree
        ]
        complement_columns = []
        for complement_monomial in complement:
            complement_columns.append(
                sp.Matrix(
                    [
                        sp.Rational(1) if monomial == complement_monomial else sp.Rational(0)
                        for monomial in monomials
                    ]
                )
            )
        square = sp.Matrix.hstack(*(selected_columns + complement_columns))
        assert square.rows == square.cols == len(monomials)
        plans[degree] = {
            "monomials": monomials,
            "ideal_labels": selected_labels,
            "complement": complement,
            "inverse": square.inv(),
        }
    return plans


def macaulay_reduce(
    polynomial: dict[tuple[int, int, int], Any],
    equations: list[sp.Poly],
    homogeneous: list[sp.Poly],
    plans: dict[int, dict[str, object]],
) -> dict[tuple[int, int, int], Any]:
    pending = {monomial: coefficient for monomial, coefficient in polynomial.items() if coefficient}
    answer: dict[tuple[int, int, int], Any] = {}
    lower_parts: list[dict[tuple[int, int, int], Any]] = []
    for equation in equations[:3]:
        lower_parts.append(
            {
                monomial: coefficient
                for monomial, coefficient in equation.rep.to_dict().items()
                if sum(monomial) < 2
            }
        )

    maximum = max((sum(monomial) for monomial in pending), default=0)
    for degree in range(maximum, -1, -1):
        plan = plans[degree]
        monomials = plan["monomials"]
        assert isinstance(monomials, list)
        top = [pending.pop(monomial, ZERO) for monomial in monomials]
        inverse = plan["inverse"]
        assert isinstance(inverse, sp.MatrixBase)
        solution: list[Any] = []
        for row in range(inverse.rows):
            value = ZERO
            for column in range(inverse.cols):
                if inverse[row, column]:
                    value += QQQ.convert(inverse[row, column]) * top[column]
            solution.append(value)
        labels = plan["ideal_labels"]
        complement = plan["complement"]
        assert isinstance(labels, list) and isinstance(complement, list)
        ideal_count = len(labels)
        for offset, monomial in enumerate(complement):
            coefficient = solution[ideal_count + offset]
            if coefficient:
                answer[monomial] = answer.get(monomial, ZERO) + coefficient
        for index, (equation_index, multiplier) in enumerate(labels):
            coefficient = solution[index]
            if not coefficient:
                continue
            for monomial, lower_coefficient in lower_parts[equation_index].items():
                shifted = add_monomials(monomial, multiplier)
                updated = pending.get(shifted, ZERO) - coefficient * lower_coefficient
                if updated:
                    pending[shifted] = updated
                elif shifted in pending:
                    del pending[shifted]
    assert not pending
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def independent_exceptional_polynomial(equations: list[sp.Poly]) -> sp.Poly:
    homogeneous = leading_quadrics(equations)
    plans = macaulay_plans(homogeneous, maximum_degree=6)
    columns: list[list[Any]] = []
    seventh = equations[3].rep.to_dict()
    for basis_monomial in BASIS_MONOMIALS:
        product = {
            add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in seventh.items()
        }
        reduced = macaulay_reduce(product, equations, homogeneous, plans)
        columns.append([reduced.get(monomial, ZERO) for monomial in BASIS_MONOMIALS])
    rows = [[columns[column][row] for column in range(8)] for row in range(8)]

    denominators: list[sp.Poly] = []
    for column in range(8):
        denominator = sp.Integer(1)
        for row in range(8):
            denominator = sp.lcm(denominator, rows[row][column].denom.as_expr())
        denominators.append(sp.Poly(denominator, QPAR, domain=QQ).monic())
    cleared_rows: list[list[Any]] = [[None] * 8 for _ in range(8)]
    for column, denominator in enumerate(denominators):
        multiplier = QQ_POLY_Q.from_sympy(denominator.as_expr())
        for row in range(8):
            value = rows[row][column] * QQQ.new(multiplier)
            assert value.denom.degree() == 0
            constant = next(iter(value.denom.to_dict().values()))
            cleared_rows[row][column] = value.numer / constant
    determinant = DomainMatrix(cleared_rows, (8, 8), QQ_POLY_Q).det()
    determinant_poly = sp.Poly(determinant.as_expr(), QPAR, domain=QQ)
    total_denominator = QQ_POLY_Q.one
    for denominator in denominators:
        total_denominator *= QQ_POLY_Q.from_sympy(denominator.as_expr())
    denominator_poly = sp.Poly(total_denominator.as_expr(), QPAR, domain=QQ)
    numerator = determinant_poly.exquo(sp.gcd(determinant_poly, denominator_poly))
    _, numerator = numerator.primitive()
    _, integer = numerator.clear_denoms(convert=True)
    _, integer = integer.primitive()
    if integer.LC() < 0:
        integer = -integer
    return integer


def verify_chain_branch(artifact: dict[str, object]) -> bool:
    epsilon, shift, traces = traces_over_zzq(BONDS_CHAIN)
    targets = normalized_targets(traces, epsilon, shift)
    equations, _, _ = norm_equations(targets)
    control = artifact["data"]["open_chain_control"]["explicit_branch"]
    substitution = {
        A5: sp.sympify(control["a5"], locals={"q": QPAR}),
        A4: sp.sympify(control["a4"], locals={"q": QPAR}),
        A3: sp.sympify(control["a3"], locals={"q": QPAR}),
    }
    return all(sp.cancel(equation.as_expr().subs(substitution)) == 0 for equation in equations)
def integer_shift_at_one(coefficients: list[int]) -> list[int]:
    shifted = [0] * len(coefficients)
    for degree, coefficient in enumerate(coefficients):
        binomial = 1
        for new_degree in range(degree + 1):
            shifted[new_degree] += coefficient * binomial
            if new_degree < degree:
                binomial = (
                    binomial * (degree - new_degree) // (new_degree + 1)
                )
    return shifted


def sign_variations(coefficients: list[int]) -> int:
    signs = [1 if value > 0 else -1 for value in coefficients if value]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def rational_polynomial_sign(coefficients: list[int], value: sp.Rational) -> int:
    numerator = int(value.p)
    denominator = int(value.q)
    total = coefficients[-1]
    denominator_power = 1
    for coefficient in reversed(coefficients[:-1]):
        denominator_power *= denominator
        total = total * numerator + coefficient * denominator_power
    return 1 if total > 0 else -1 if total < 0 else 0


def independent_t_pullback(coefficients: list[int]) -> list[int]:
    degree = len(coefficients) - 1
    pulled = [0] * (2 * degree + 1)
    powers_of_two = [1] * (degree + 1)
    for exponent in range(1, degree + 1):
        powers_of_two[exponent] = 2 * powers_of_two[exponent - 1]
    for q_degree, coefficient in enumerate(coefficients):
        binomial = 1
        base_degree = degree - q_degree
        scale = coefficient * powers_of_two[base_degree]
        for choice in range(q_degree + 1):
            pulled[base_degree + 2 * choice] += scale * binomial
            if choice < q_degree:
                binomial = (
                    binomial * (q_degree - choice) // (choice + 1)
                )
    content = 0
    for coefficient in pulled:
        content = math.gcd(content, abs(coefficient))
    if content == 0:
        raise AssertionError("zero pullback polynomial")
    primitive = [coefficient // content for coefficient in pulled]
    if primitive[-1] < 0:
        primitive = [-coefficient for coefficient in primitive]
    return primitive




def main() -> int:
    stage = os.environ.get("TRACE_EXCEPTIONAL_STAGE", "all")
    if stage not in {"all", "elimination", "roots"}:
        raise ValueError("TRACE_EXCEPTIONAL_STAGE must be all, elimination, or roots")
    artifact = json.loads(ARTIFACT.read_text())
    source_hashes = artifact["meta"]["source_sha256"]
    check(
        "producer_source_hash",
        source_hashes["experiments/e238_trace_exceptional_set.py"] == file_sha256(PRODUCER),
    )
    check(
        "verifier_source_hash",
        source_hashes["tests/test_trace_exceptional_set.py"] == file_sha256(Path(__file__).resolve()),
    )
    check("benchmark_absent", artifact["meta"]["benchmark_Kc_used"] is False)

    stored_q = artifact["data"]["elimination"]["exceptional_q_polynomial"]
    q_coefficients = [int(value) for value in stored_q["coefficients"]]
    exceptional_q = sp.Poly.from_list(list(reversed(q_coefficients)), QPAR, domain=ZZ)
    check("stored_q_degree", exceptional_q.degree() == stored_q["degree"] == 971)
    check(
        "stored_q_digest",
        canonical_sha256([str(value) for value in q_coefficients]) == stored_q["sha256"],
    )
    check("stored_q_primitive", sp.gcd_list(exceptional_q.all_coeffs()) == 1)
    check("stored_q_squarefree", sp.gcd(exceptional_q, exceptional_q.diff()).degree() == 0)

    if stage != "roots":
        epsilon, shift, traces = traces_over_zzq(BONDS_2X3)
        targets = normalized_targets(traces, epsilon, shift)
        equations, _, _ = norm_equations(targets)
        independent = independent_exceptional_polynomial(equations)
        check(
            "independent_macaulay_elimination",
            independent == exceptional_q,
            "degree-filtered affine Macaulay reduction reproduces the primitive quotient norm",
        )
        if stage == "elimination":
            if FAILURES:
                print("FAILURES:", ", ".join(FAILURES))
                return 1
            print("PASS: independent trace exceptional-set elimination verified")
            return 0

    root_data = artifact["data"]["physical_root_certificate"]
    shifted_coefficients = integer_shift_at_one(q_coefficients)
    variations = sign_variations(shifted_coefficients)
    check(
        "independent_descartes_upper_bound",
        variations == root_data["sign_variations_E_of_one_plus_x"] == 3,
        "pure-integer binomial shift gives at most three q>1 roots",
    )
    bracket_signs: list[list[int]] = []
    q_intervals: list[tuple[sp.Rational, sp.Rational]] = []
    t_intervals: list[tuple[sp.Rational, sp.Rational]] = []
    mapping_ok = True
    for row in root_data["roots"]:
        q_lower, q_upper = map(parse_rational, row["q_interval"])
        q_intervals.append((q_lower, q_upper))
        t_lower, t_upper = map(parse_rational, row["t_interval"])
        t_intervals.append((t_lower, t_upper))
        signs = [
            rational_polynomial_sign(q_coefficients, q_lower),
            rational_polynomial_sign(q_coefficients, q_upper),
        ]
        bracket_signs.append(signs)
        mapping_ok &= signs == row["q_endpoint_signs"]
        mapping_ok &= signs[0] * signs[1] < 0
        q_of_t_lower = (1 + t_lower**2) / (2 * t_lower)
        q_of_t_upper = (1 + t_upper**2) / (2 * t_upper)
        mapping_ok &= q_of_t_lower > q_upper and q_of_t_upper < q_lower
    check(
        "independent_three_sign_brackets",
        len(bracket_signs) == 3
        and all(left * right < 0 for left, right in bracket_signs),
        str(bracket_signs),
    )
    physical_disjoint = all(
        sp.Rational(1) < lower < upper for lower, upper in q_intervals
    ) and all(
        left_upper < right_lower
        for (_, left_upper), (right_lower, _) in zip(
            q_intervals, q_intervals[1:]
        )
    )
    check(
        "independent_brackets_physical_and_disjoint",
        len(q_intervals) == 3 and physical_disjoint,
        str(q_intervals),
    )
    sorted_t_intervals = sorted(t_intervals)
    physical_t_disjoint = all(
        sp.Rational(0) < lower < upper < sp.Rational(1)
        for lower, upper in t_intervals
    ) and all(
        left_upper < right_lower
        for (_, left_upper), (right_lower, _) in zip(
            sorted_t_intervals, sorted_t_intervals[1:]
        )
    )
    check(
        "independent_t_brackets_physical_and_disjoint",
        len(t_intervals) == 3 and physical_t_disjoint,
        str(t_intervals),
    )
    check(
        "exactly_three_physical_q_roots",
        variations == len(bracket_signs) == len(q_intervals) == 3
        and physical_disjoint,
        "Descartes upper bound plus three disjoint brackets",
    )
    check(
        "exact_q_to_t_interval_map",
        mapping_ok and physical_t_disjoint,
    )
    check(
        "surviving_branch_dispositions",
        all(
            "survives" in row["disposition"]
            and "unresolved" in row["disposition"]
            for row in root_data["roots"]
        ),
    )

    stored_t = artifact["data"]["elimination"]["exceptional_t_polynomial"]
    t_coefficients = [int(value) for value in stored_t["coefficients"]]
    reconstructed_t = independent_t_pullback(q_coefficients)
    check(
        "independent_q_to_t_polynomial_pullback",
        reconstructed_t == t_coefficients
        and len(t_coefficients) == 1943
        and stored_t["degree"] == 1942,
    )
    check(
        "stored_t_reciprocity",
        t_coefficients == list(reversed(t_coefficients)),
    )
    check(
        "stored_t_digest",
        canonical_sha256([str(value) for value in t_coefficients])
        == stored_t["sha256"],
    )
    check("symbolic_chain_branch", verify_chain_branch(artifact))

    theorem = artifact["data"]["theorem"]
    check(
        "scope_preserved",
        "every physical coupling" in theorem["scope"]
        and "no all-size" in theorem["scope"]
        and "does not prove an empty" in theorem["method_boundary"],
    )
    if FAILURES:
        print("FAILURES:", ", ".join(FAILURES))
        return 1
    print("PASS: independent trace exceptional-set certificate verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
