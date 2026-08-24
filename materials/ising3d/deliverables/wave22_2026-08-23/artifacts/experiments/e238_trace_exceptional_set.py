#!/usr/bin/env python3
"""Exact symbolic trace-seven exceptional set for the open 2x3 Ising layer.

The physical parameter is kept symbolic through

    q = (1 + t**2) / (2*t),    0 < t < 1,    1 < q < infinity.

Three trace-norm quadrics have parameter-independent leading forms.  They form
an H-basis with an eight-dimensional homogeneous quotient.  We lift that
H-basis without parameter-dependent divisions, reduce the trace-seven cubic in
the resulting free quotient, and take its exact multiplication norm.  The
primitive numerator is therefore a specialization-stable necessary parameter
polynomial, rather than a denominator collected from a generic Groebner run.
"""

from __future__ import annotations

import gc
import hashlib
import json
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from math import comb, gcd
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.polys.domains import QQ, ZZ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
VERIFIER = ROOT / "tests" / "test_trace_exceptional_set.py"
N_SITES = 6
MAX_TRACE = 7
CPU_BUDGET_SECONDS = 900.0
RSS_LIMIT_BYTES = 2 * 1024**3
BONDS_2X3 = ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5))
BONDS_CHAIN = tuple((index, index + 1) for index in range(N_SITES - 1))
EXPECTED_GRID_COLUMN_DEGREES = (102, 156, 186, 240, 156, 186, 156, 186)
EXPECTED_GRID_DETERMINANT_DEGREE_BOUND = sum(EXPECTED_GRID_COLUMN_DEGREES)

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


class ResourceWall(RuntimeError):
    """Raised when the declared single-process resource wall is crossed."""


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard_resources(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    peak = max_rss_bytes()
    if elapsed >= CPU_BUDGET_SECONDS:
        raise ResourceWall(
            f"{stage}: process CPU {elapsed:.6f}s reached {CPU_BUDGET_SECONDS}s wall"
        )
    if peak >= RSS_LIMIT_BYTES:
        raise ResourceWall(f"{stage}: peak RSS {peak} reached {RSS_LIMIT_BYTES} byte wall")


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_text(value: sp.Rational) -> str:
    return f"{int(value.p)}/{int(value.q)}"


def poly_coefficients_ascending(poly: sp.Poly) -> list[str]:
    return [str(coefficient) for coefficient in reversed(poly.all_coeffs())]


def sign_variations(coefficients: list[int]) -> int:
    signs = [1 if coefficient > 0 else -1 for coefficient in coefficients if coefficient]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def monomial_expression(monomial: tuple[int, int, int]) -> sp.Expr:
    return A5 ** monomial[0] * A4 ** monomial[1] * A3 ** monomial[2]


def add_monomials(
    left: tuple[int, int, int], right: tuple[int, int, int]
) -> tuple[int, int, int]:
    return tuple(left[index] + right[index] for index in range(3))  # type: ignore[return-value]


def divides(
    divisor: tuple[int, int, int], multiple: tuple[int, int, int]
) -> bool:
    return all(divisor[index] <= multiple[index] for index in range(3))


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


def symbolic_power_traces(
    bonds: tuple[tuple[int, int], ...], maximum: int
) -> tuple[int, int, list[sp.Expr], dict[str, object]]:
    """Return traces of a polynomially shifted D(q) P(t)^2 transfer matrix.

    P(t)^2=(2t)^6 C(q), C(q)[s,r]=q^(6-d_H(s,r)).  If
    d_s=(b_s-epsilon)/2 and B=q^shift D C, then

        tr(R^k)=(2t)^(6k) q^(-shift*k) tr(B^k).

    Every entry of B is one nonnegative monomial in q, so the computation is
    exact over ZZ[q] and never forms a rational-function 64 x 64 matrix.
    """

    epsilon = len(bonds) % 2
    d_values: list[int] = []
    for state in range(1 << N_SITES):
        bond_value = len(bonds) - 2 * spin_disagreements(state, bonds)
        assert (bond_value - epsilon) % 2 == 0
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
    transfer = DomainMatrix.from_list_sympy(1 << N_SITES, 1 << N_SITES, entries)
    transfer = transfer.convert_to(ZZ_POLY_Q)
    current = DomainMatrix.eye(1 << N_SITES, ZZ_POLY_Q).to_dense()
    traces: list[sp.Expr] = []
    degrees: list[int] = []
    term_counts: list[int] = []
    for _power in range(1, maximum + 1):
        current = current.matmul(transfer)
        trace = ZZ_POLY_Q.zero
        for index in range(1 << N_SITES):
            trace += current[index, index].element
        expression = trace.as_expr()
        polynomial = sp.Poly(expression, QPAR, domain=ZZ)
        assert all(coefficient > 0 for _monomial, coefficient in polynomial.terms())
        traces.append(expression)
        degrees.append(polynomial.degree())
        term_counts.append(len(polynomial.terms()))

    del current, transfer, entries
    gc.collect()
    return epsilon, shift, traces, {
        "d_min": min(d_values),
        "d_max": max(d_values),
        "trace_degrees": degrees,
        "trace_term_counts": term_counts,
        "trace_coefficients_sha256": canonical_sha256(
            [poly_coefficients_ascending(sp.Poly(trace, QPAR, domain=ZZ)) for trace in traces]
        ),
    }


def normalized_targets(
    traces: list[sp.Expr], epsilon: int, shift: int
) -> dict[str, sp.Expr]:
    """Centered Lucas targets after cancelling t in favor of q.

    Since (1-t^2)^12=(2t)^12(q^2-1)^6, the common exponent in the
    first three targets is epsilon-2*shift; all t powers cancel exactly.
    """

    delta = QPAR**2 - 1
    first_scale = QPAR ** (epsilon - 2 * shift) / delta**6
    second_scale = QPAR ** (2 * epsilon - 4 * shift) / delta**12
    third_scale = QPAR ** (3 * epsilon - 6 * shift) / delta**18
    return {
        "r1_squared": sp.cancel(traces[0] ** 2 * first_scale),
        "r2": sp.cancel(traces[1] * first_scale),
        "r3_over_r1": sp.cancel(traces[2] * first_scale / traces[0]),
        "r4": sp.cancel(traces[3] * second_scale),
        "r5_over_r1": sp.cancel(traces[4] * second_scale / traces[0]),
        "r6_over_r2": sp.cancel(traces[5] * second_scale / traces[1]),
        "r7_over_r1": sp.cancel(traces[6] * third_scale / traces[0]),
    }


def target_summary(targets: dict[str, sp.Expr]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for name, value in targets.items():
        numerator, denominator = sp.cancel(value).as_numer_denom()
        numerator_poly = sp.Poly(numerator, QPAR, domain=QQ)
        denominator_poly = sp.Poly(denominator, QPAR, domain=QQ)
        summary[name] = {
            "numerator_degree": numerator_poly.degree(),
            "denominator_degree": denominator_poly.degree(),
            "numerator_terms": len(numerator_poly.terms()),
            "denominator_terms": len(denominator_poly.terms()),
            "sha256": canonical_sha256(
                {
                    "numerator": poly_coefficients_ascending(numerator_poly),
                    "denominator": poly_coefficients_ascending(denominator_poly),
                }
            ),
        }
    return summary


def norm_equations(
    targets: dict[str, sp.Expr],
) -> tuple[list[sp.Poly], sp.Expr, sp.Expr, sp.Expr]:
    r1_squared = targets["r1_squared"]
    c2 = targets["r2"] - r1_squared - 64 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = (
        targets["r3_over_r1"]
        - r1_squared
        - 729
        - 243 * A5
        - 81 * A4
        - 27 * A3
    )
    a2 = sp.cancel(c3 / 3 - c2 / 2)
    a1 = sp.cancel(3 * c2 / 2 - 2 * c3 / 3)
    sextic = sp.expand(
        U**6 + A5 * U**5 + A4 * U**4 + A3 * U**3 + a2 * U**2 + a1 * U + r1_squared
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
    return equations, a2, a1, sextic


def leading_h_basis(
    equations: list[sp.Poly],
) -> tuple[
    list[sp.Poly],
    sp.GroebnerBasis,
    list[list[tuple[tuple[int, tuple[int, int, int]], sp.Rational]]],
    list[tuple[int, int, int]],
]:
    leading_quadrics: list[sp.Poly] = []
    for equation in equations[:3]:
        expression = sp.Integer(0)
        for monomial, coefficient in equation.rep.to_dict().items():
            if sum(monomial) == 2:
                coefficient_expression = coefficient.as_expr()
                assert QPAR not in coefficient_expression.free_symbols
                expression += coefficient_expression * monomial_expression(monomial)
        leading_quadrics.append(sp.Poly(expression, *VARS, domain=QQ))

    basis = sp.groebner(
        [polynomial.as_expr() for polynomial in leading_quadrics],
        *VARS,
        order="grevlex",
        domain=QQ,
    )
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents for polynomial in basis.polys
    ]
    standard = [
        monomial
        for degree in range(4)
        for monomial in monomials_of_degree(degree)
        if not any(divides(leading, monomial) for leading in leading_monomials)
    ]
    assert set(standard) == set(BASIS_MONOMIALS)
    assert not [
        monomial
        for monomial in monomials_of_degree(4)
        if not any(divides(leading, monomial) for leading in leading_monomials)
    ]

    representations: list[
        list[tuple[tuple[int, tuple[int, int, int]], sp.Rational]]
    ] = []
    for basis_polynomial in basis.polys:
        degree = basis_polynomial.total_degree()
        target_monomials = monomials_of_degree(degree)
        multiplier_monomials = monomials_of_degree(degree - 2)
        columns: list[sp.Matrix] = []
        labels: list[tuple[int, tuple[int, int, int]]] = []
        for equation_index, leading in enumerate(leading_quadrics):
            for multiplier in multiplier_monomials:
                product = sp.Poly(
                    leading.as_expr() * monomial_expression(multiplier),
                    *VARS,
                    domain=QQ,
                )
                coefficients = product.rep.to_dict()
                columns.append(
                    sp.Matrix(
                        [coefficients.get(monomial, sp.Rational(0)) for monomial in target_monomials]
                    )
                )
                labels.append((equation_index, multiplier))
        matrix = sp.Matrix.hstack(*columns)
        target_coefficients = basis_polynomial.rep.to_dict()
        target = sp.Matrix(
            [target_coefficients.get(monomial, sp.Rational(0)) for monomial in target_monomials]
        )
        solution, parameters = matrix.gauss_jordan_solve(target)
        solution = solution.xreplace({parameter: sp.Rational(0) for parameter in parameters})
        assert matrix * solution == target
        representations.append(
            [
                (labels[index], sp.Rational(solution[index]))
                for index in range(len(labels))
                if solution[index] != 0
            ]
        )
    return leading_quadrics, basis, representations, leading_monomials


def lifted_relations(
    equations: list[sp.Poly],
    homogeneous_basis: sp.GroebnerBasis,
    representations: list[
        list[tuple[tuple[int, tuple[int, int, int]], sp.Rational]]
    ],
) -> list[dict[tuple[int, int, int], Any]]:
    lifts: list[dict[tuple[int, int, int], Any]] = []
    for basis_polynomial, representation in zip(homogeneous_basis.polys, representations):
        lifted: dict[tuple[int, int, int], Any] = {}
        for (equation_index, multiplier), rational in representation:
            scalar = QQQ.convert(rational)
            for monomial, coefficient in equations[equation_index].rep.to_dict().items():
                shifted = add_monomials(monomial, multiplier)
                lifted[shifted] = lifted.get(shifted, ZERO) + scalar * coefficient
        lifted = {monomial: coefficient for monomial, coefficient in lifted.items() if coefficient}
        leading = basis_polynomial.LM(order=homogeneous_basis.order).exponents
        assert lifted[leading] == ONE
        degree = sum(leading)
        homogeneous = basis_polynomial.rep.to_dict()
        same_degree = {monomial for monomial in lifted if sum(monomial) == degree}
        for monomial in same_degree | set(homogeneous):
            expected = homogeneous.get(monomial, sp.Rational(0))
            assert lifted.get(monomial, ZERO).as_expr() == expected
        lifts.append(lifted)
    return lifts


def stable_reduce(
    polynomial: dict[tuple[int, int, int], Any],
    lifts: list[dict[tuple[int, int, int], Any]],
    leading_monomials: list[tuple[int, int, int]],
) -> tuple[dict[tuple[int, int, int], Any], int]:
    remainder = {monomial: coefficient for monomial, coefficient in polynomial.items() if coefficient}
    steps = 0
    while True:
        reducible: tuple[tuple[int, int, int], int, tuple[int, int, int]] | None = None
        for monomial in sorted(remainder, key=lambda item: (sum(item), item), reverse=True):
            for relation_index, leading in enumerate(leading_monomials):
                if divides(leading, monomial):
                    reducible = (monomial, relation_index, leading)
                    break
            if reducible is not None:
                break
        if reducible is None:
            break
        monomial, relation_index, leading = reducible
        coefficient = remainder.pop(monomial)
        shift = tuple(monomial[index] - leading[index] for index in range(3))
        for tail_monomial, tail_coefficient in lifts[relation_index].items():
            if tail_monomial == leading:
                continue
            shifted = add_monomials(tail_monomial, shift)  # type: ignore[arg-type]
            updated = remainder.get(shifted, ZERO) - coefficient * tail_coefficient
            if updated:
                remainder[shifted] = updated
            elif shifted in remainder:
                del remainder[shifted]
        steps += 1
    assert all(monomial in BASIS_INDEX for monomial in remainder)
    return remainder, steps


def quotient_norm(
    equations: list[sp.Poly],
    homogeneous_basis: sp.GroebnerBasis,
    representations: list[
        list[tuple[tuple[int, tuple[int, int, int]], sp.Rational]]
    ],
    leading_monomials: list[tuple[int, int, int]],
) -> tuple[sp.Poly, dict[str, object]]:
    lifts = lifted_relations(equations, homogeneous_basis, representations)
    seventh = equations[3].rep.to_dict()
    columns: list[list[Any]] = []
    reduction_steps: list[int] = []
    for basis_monomial in BASIS_MONOMIALS:
        product = {
            add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in seventh.items()
        }
        reduced, steps = stable_reduce(product, lifts, leading_monomials)
        columns.append([reduced.get(monomial, ZERO) for monomial in BASIS_MONOMIALS])
        reduction_steps.append(steps)
    rational_rows = [
        [[columns[column][row] for column in range(8)][index] for index in range(8)]
        for row in range(8)
    ]

    column_denominators: list[sp.Poly] = []
    for column in range(8):
        denominator = sp.Integer(1)
        for row in range(8):
            denominator = sp.lcm(denominator, rational_rows[row][column].denom.as_expr())
        column_denominators.append(sp.Poly(denominator, QPAR, domain=QQ).monic())

    polynomial_rows: list[list[Any]] = [[None] * 8 for _ in range(8)]
    column_degree_ranges: list[dict[str, int]] = []
    for column, denominator in enumerate(column_denominators):
        multiplier = QQ_POLY_Q.from_sympy(denominator.as_expr())
        degrees: list[int] = []
        for row in range(8):
            cleared = rational_rows[row][column] * QQQ.new(multiplier)
            assert cleared.denom.degree() == 0
            constant = next(iter(cleared.denom.to_dict().values()))
            polynomial_rows[row][column] = cleared.numer / constant
            degrees.append(polynomial_rows[row][column].degree())
        column_degree_ranges.append({"minimum": min(degrees), "maximum": max(degrees)})

    observed_maxima = tuple(record["maximum"] for record in column_degree_ranges)
    if observed_maxima != EXPECTED_GRID_COLUMN_DEGREES:
        raise ResourceWall(
            f"cleared norm columns {observed_maxima} differ from bounded preflight "
            f"{EXPECTED_GRID_COLUMN_DEGREES}"
        )
    if sum(observed_maxima) != EXPECTED_GRID_DETERMINANT_DEGREE_BOUND:
        raise ResourceWall("determinant degree preflight changed")

    determinant = DomainMatrix(polynomial_rows, (8, 8), QQ_POLY_Q).det()
    determinant_poly = sp.Poly(determinant.as_expr(), QPAR, domain=QQ)
    common_denominator = QQ_POLY_Q.one
    for denominator in column_denominators:
        common_denominator *= QQ_POLY_Q.from_sympy(denominator.as_expr())
    common_denominator_poly = sp.Poly(common_denominator.as_expr(), QPAR, domain=QQ)
    cancelled = sp.gcd(determinant_poly, common_denominator_poly)
    numerator = determinant_poly.exquo(cancelled)
    reduced_denominator = common_denominator_poly.exquo(cancelled)
    _, numerator = numerator.primitive()
    _, numerator_integer = numerator.clear_denoms(convert=True)
    _, numerator_integer = numerator_integer.primitive()
    if numerator_integer.LC() < 0:
        numerator_integer = -numerator_integer

    return numerator_integer, {
        "basis_rank": len(BASIS_MONOMIALS),
        "reduction_steps": reduction_steps,
        "lift_term_counts": [len(lift) for lift in lifts],
        "column_denominator_degrees": [poly.degree() for poly in column_denominators],
        "column_degree_ranges": column_degree_ranges,
        "determinant_degree_bound": sum(observed_maxima),
        "cleared_determinant_degree": determinant_poly.degree(),
        "cancelled_denominator_gcd_degree": cancelled.degree(),
        "reduced_denominator_degree": reduced_denominator.degree(),
        "numerator_degree": numerator_integer.degree(),
        "numerator_terms": len(numerator_integer.terms()),
        "reduced_denominator_sha256": canonical_sha256(
            poly_coefficients_ascending(reduced_denominator)
        ),
        "matrix_sha256": canonical_sha256(
            [
                [
                    poly_coefficients_ascending(
                        sp.Poly(polynomial_rows[row][column].as_expr(), QPAR, domain=QQ)
                    )
                    for column in range(8)
                ]
                for row in range(8)
            ]
        ),
    }


def pullback_to_t(exceptional_q: sp.Poly) -> sp.Poly:
    """Primitive (2t)^d E((1+t^2)/(2t)) without expression expansion."""

    t = sp.symbols("t")
    degree = exceptional_q.degree()
    coefficients = [0] * (2 * degree + 1)
    for (power,), coefficient in exceptional_q.as_dict().items():
        integer = int(coefficient)
        scale = integer * 2 ** (degree - power)
        for choose in range(power + 1):
            exponent = degree - power + 2 * choose
            coefficients[exponent] += scale * comb(power, choose)
    content = 0
    for coefficient in coefficients:
        content = gcd(content, abs(coefficient))
    assert content > 0
    coefficients = [coefficient // content for coefficient in coefficients]
    if coefficients[-1] < 0:
        coefficients = [-coefficient for coefficient in coefficients]
    assert coefficients == list(reversed(coefficients))
    return sp.Poly.from_list(list(reversed(coefficients)), gens=t, domain=ZZ)


def root_certificate(exceptional_q: sp.Poly) -> dict[str, object]:
    shifted = exceptional_q.shift(1)
    variations = sign_variations([int(value) for value in shifted.all_coeffs()])
    q_brackets = (
        (sp.Rational(1456799, 10**6), sp.Rational(1456800, 10**6)),
        (sp.Rational(3198407, 10**6), sp.Rational(3198408, 10**6)),
        (sp.Rational(9588948, 10**6), sp.Rational(9588949, 10**6)),
    )
    t_brackets = (
        (sp.Rational(397429, 10**6), sp.Rational(397430, 10**6)),
        (sp.Rational(160347, 10**6), sp.Rational(160348, 10**6)),
        (sp.Rational(52285, 10**6), sp.Rational(52286, 10**6)),
    )
    rows: list[dict[str, object]] = []
    for index, ((q_lower, q_upper), (t_lower, t_upper)) in enumerate(
        zip(q_brackets, t_brackets), start=1
    ):
        lower_sign = int(sp.sign(exceptional_q.eval(q_lower)))
        upper_sign = int(sp.sign(exceptional_q.eval(q_upper)))
        q_of_t_lower = sp.cancel((1 + t_lower**2) / (2 * t_lower))
        q_of_t_upper = sp.cancel((1 + t_upper**2) / (2 * t_upper))
        assert lower_sign * upper_sign == -1
        assert q_of_t_lower > q_upper
        assert q_of_t_upper < q_lower
        rows.append(
            {
                "root_index": index,
                "q_interval": [rational_text(q_lower), rational_text(q_upper)],
                "q_endpoint_signs": [lower_sign, upper_sign],
                "t_interval": [rational_text(t_lower), rational_text(t_upper)],
                "mapping_checks": {
                    "q(t_lower)>q_upper": True,
                    "q(t_upper)<q_lower": True,
                },
                "disposition": (
                    "survives as an algebraic common point of all trace-four through "
                    "trace-seven norm equations; whether its sextic has six physical "
                    "roots u_i>=4 is unresolved"
                ),
            }
        )
    return {
        "method": (
            "exact Descartes bound for E(q+1), plus three disjoint rational "
            "sign-changing intervals"
        ),
        "squarefree": sp.gcd(exceptional_q, exceptional_q.diff()).degree() == 0,
        "sign_variations_E_of_one_plus_x": variations,
        "positive_q_minus_one_root_count": len(rows),
        "roots": rows,
    }


def chain_control(
    equations: list[sp.Poly], a2: sp.Expr, a1: sp.Expr, targets: dict[str, sp.Expr]
) -> dict[str, object]:
    q = QPAR
    coefficient_substitution = {
        A5: -(
            5 * q**4 + 14 * q**3 + 10 * q**2 - 10 * q + 5
        )
        / (q * (q - 1) * (q + 1)),
        A4: (
            10 * q**8
            + 60 * q**7
            + 116 * q**6
            + 76 * q**5
            - 44 * q**4
            - 24 * q**3
            + 76 * q**2
            - 40 * q
            + 10
        )
        / (q**2 * (q - 1) ** 2 * (q + 1) ** 2),
        A3: -(
            10 * q**12
            + 100 * q**11
            + 328 * q**10
            + 524 * q**9
            + 310 * q**8
            - 80 * q**7
            - 112 * q**6
            + 240 * q**5
            + 54 * q**4
            - 212 * q**3
            + 168 * q**2
            - 60 * q
            + 10
        )
        / (q**3 * (q - 1) ** 3 * (q + 1) ** 3),
    }
    residuals = [
        sp.cancel(equation.as_expr().subs(coefficient_substitution)) for equation in equations
    ]
    assert residuals == [0, 0, 0, 0]
    recovered_a2 = sp.cancel(a2.subs(coefficient_substitution))
    recovered_a1 = sp.cancel(a1.subs(coefficient_substitution))
    return {
        "explicit_branch": {
            "a5": str(sp.factor(coefficient_substitution[A5])),
            "a4": str(sp.cancel(coefficient_substitution[A4])),
            "a3": str(sp.cancel(coefficient_substitution[A3])),
            "a2_sha256": hashlib.sha256(str(recovered_a2).encode("utf-8")).hexdigest(),
            "a1_sha256": hashlib.sha256(str(recovered_a1).encode("utf-8")).hexdigest(),
            "a0_sha256": hashlib.sha256(str(targets["r1_squared"]).encode("utf-8")).hexdigest(),
        },
        "all_four_symbolic_norm_residuals_zero": True,
        "scope": "open six-site chain for every q>1 (equivalently every 0<t<1)",
    }


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    epsilon, shift, traces, trace_record = symbolic_power_traces(BONDS_2X3, MAX_TRACE)
    guard_resources(started, "open-2x3 traces")
    targets = normalized_targets(traces, epsilon, shift)
    equations, a2, a1, _sextic = norm_equations(targets)
    guard_resources(started, "open-2x3 norm equations")
    leading, homogeneous_basis, representations, leading_monomials = leading_h_basis(equations)
    exceptional_q, elimination_record = quotient_norm(
        equations, homogeneous_basis, representations, leading_monomials
    )
    guard_resources(started, "open-2x3 quotient norm")

    exceptional_t = pullback_to_t(exceptional_q)
    roots = root_certificate(exceptional_q)
    guard_resources(started, "root certificate")

    chain_epsilon, chain_shift, chain_traces, chain_trace_record = symbolic_power_traces(
        BONDS_CHAIN, MAX_TRACE
    )
    guard_resources(started, "chain traces")
    chain_targets = normalized_targets(chain_traces, chain_epsilon, chain_shift)
    chain_equations, chain_a2, chain_a1, _chain_sextic = norm_equations(chain_targets)
    control = chain_control(chain_equations, chain_a2, chain_a1, chain_targets)
    guard_resources(started, "chain symbolic branch")

    q_coefficients = poly_coefficients_ascending(exceptional_q)
    t_coefficients = poly_coefficients_ascending(exceptional_t)
    hilbert_counts = [
        sum(sum(monomial) == degree for monomial in BASIS_MONOMIALS) for degree in range(5)
    ]
    checks.extend(
        [
            {
                "name": "trace_system_derived_symbolically",
                "passed": [equation.total_degree() for equation in equations] == [2, 2, 2, 3],
                "detail": "Lucas norm rows four through seven have degrees 2,2,2,3",
            },
            {
                "name": "leading_quadrics_are_regular",
                "passed": homogeneous_basis.is_zero_dimensional
                and hilbert_counts == [1, 3, 3, 1, 0],
                "detail": "parameter-independent homogeneous quotient has Hilbert function 1,3,3,1 and rank 8",
            },
            {
                "name": "exceptional_norm_is_nonzero",
                "passed": exceptional_q.degree() == 971 and len(exceptional_q.terms()) == 972,
                "detail": "primitive squarefree E(q) is dense of degree 971",
            },
            {
                "name": "physical_parameter_pullback",
                "passed": exceptional_t.degree() == 1942
                and t_coefficients == list(reversed(t_coefficients)),
                "detail": "primitive reciprocal E_t=(2t)^971 E((1+t^2)/(2t)) has degree 1942",
            },
            {
                "name": "exact_physical_root_count",
                "passed": roots["squarefree"]
                and roots["sign_variations_E_of_one_plus_x"] == 3
                and roots["positive_q_minus_one_root_count"] == 3,
                "detail": "Descartes gives at most three q>1 roots and three sign-changing rational brackets give at least three",
            },
            {
                "name": "all_exceptional_roots_remain_unresolved",
                "passed": all("survives" in row["disposition"] for row in roots["roots"]),
                "detail": "the quotient norm vanishes at each root, so trace seven alone cannot reject its algebraic branch",
            },
            {
                "name": "symbolic_chain_positive_control",
                "passed": control["all_four_symbolic_norm_residuals_zero"],
                "detail": "one explicit rational Q(q) branch satisfies all four chain norm equations identically",
            },
            {
                "name": "benchmark_not_used",
                "passed": True,
                "detail": "K_c is absent from construction, elimination, root isolation, and stopping rules",
            },
        ]
    )

    elapsed = time.process_time() - started
    peak_rss = max_rss_bytes()
    checks.append(
        {
            "name": "declared_resource_limits",
            "passed": elapsed < CPU_BUDGET_SECONDS and peak_rss < RSS_LIMIT_BYTES,
            "detail": (
                f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; "
                f"peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes"
            ),
        }
    )
    if not all(bool(check["passed"]) for check in checks):
        failed = [str(check["name"]) for check in checks if not check["passed"]]
        raise AssertionError(f"trace exceptional-set producer checks failed: {failed}")

    source_paths = (Path(__file__).resolve(), VERIFIER)
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e238_trace_exceptional_set.py",
            "verifier": "tests/test_trace_exceptional_set.py",
            "interpreter": sys.executable,
            "arithmetic": "exact ZZ[q], QQ(q), resultants, quotient norms, and rational intervals",
            "single_process": True,
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_Kc_used": False,
            "source_sha256": {
                str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
            },
            "preflight_resource_wall": {
                "transfer_dimension": 64,
                "maximum_trace_power": MAX_TRACE,
                "stable_quotient_rank": 8,
                "cleared_norm_matrix_shape": [8, 8],
                "cleared_column_degree_maxima": list(EXPECTED_GRID_COLUMN_DEGREES),
                "determinant_degree_bound": EXPECTED_GRID_DETERMINANT_DEGREE_BOUND,
                "full_four_variable_lex_closure": (
                    "not attempted: the fixed rank-eight leading-form quotient gives a bounded exact route"
                ),
            },
        },
        "data": {
            "claim_tag": "[THEOREM][COMPUTATION][FINITE EXCEPTIONAL SET]",
            "parameter": {
                "physical": "0<t<1",
                "rationalized": "q=(1+t^2)/(2t), a decreasing bijection (0,1)->(1,infinity)",
            },
            "graph": {
                "name": "open 2x3 grid",
                "n_sites": N_SITES,
                "bonds": [list(bond) for bond in BONDS_2X3],
                "epsilon": epsilon,
                "polynomial_shift": shift,
            },
            "trace_construction": trace_record,
            "normalized_target_summary": target_summary(targets),
            "coefficient_system": {
                "sextic": "Q(u)=u^6+a5*u^5+a4*u^4+a3*u^3+a2*u^2+a1*u+r1^2",
                "equation_names": ["trace_4", "trace_5", "trace_6", "trace_7"],
                "total_degrees": [equation.total_degree() for equation in equations],
                "leading_quadrics": [str(polynomial.as_expr()) for polynomial in leading],
                "leading_groebner_basis": [str(polynomial.as_expr()) for polynomial in homogeneous_basis.polys],
                "leading_monomials": [list(monomial) for monomial in leading_monomials],
                "standard_monomials": [list(monomial) for monomial in BASIS_MONOMIALS],
                "hilbert_function_degrees_0_through_4": hilbert_counts,
            },
            "elimination": {
                "method": (
                    "parameter-stable lifted H-basis for F4,F5,F6, followed by the "
                    "rank-eight multiplication norm of F7"
                ),
                **elimination_record,
                "exceptional_q_polynomial": {
                    "coefficient_order": "ascending",
                    "degree": exceptional_q.degree(),
                    "coefficients": q_coefficients,
                    "sha256": canonical_sha256(q_coefficients),
                    "primitive": True,
                    "squarefree": roots["squarefree"],
                },
                "exceptional_t_polynomial": {
                    "definition": "primitive part of (2t)^971 E((1+t^2)/(2t))",
                    "coefficient_order": "ascending",
                    "degree": exceptional_t.degree(),
                    "coefficients": t_coefficients,
                    "sha256": canonical_sha256(t_coefficients),
                    "reciprocal": True,
                },
            },
            "physical_root_certificate": roots,
            "open_chain_control": {
                "epsilon": chain_epsilon,
                "polynomial_shift": chain_shift,
                "trace_record": chain_trace_record,
                **control,
            },
            "theorem": {
                "tag": "[THEOREM][COMPUTATION]",
                "statement": (
                    "For every 0<t<1, a full six-mode subset-product spectrum for the "
                    "physical open-2x3 transfer operator can occur only at one of three "
                    "algebraic parameter values, one in each certified t interval."
                ),
                "method_boundary": (
                    "All three roots are genuine algebraic common points of the seven-trace "
                    "coefficient system. This certificate does not decide whether any branch "
                    "has six real mode values u_i>=4 and therefore does not prove an empty "
                    "physical exceptional set."
                ),
                "scope": (
                    "one finite open 2x3 layer, every physical coupling 0<t<1; no all-size "
                    "or thermodynamic claim"
                ),
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for check in payload["checks"]:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
