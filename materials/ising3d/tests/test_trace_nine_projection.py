#!/usr/bin/env python3
"""Clean-room verifier for the e251 trace-nine projection theorem.

This file imports no producer helper.  It uses direct Sylvester resultants,
Buchberger, a separately written normal-form loop, and scalar finite-field
elimination to reconstruct every stored modular certificate.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
import math
import platform
import resource
import time
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.polys.domains import QQ, ZZ

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
PRODUCER = ROOT / "experiments" / "e251_trace_nine_projection.py"
SELF = ROOT / "tests" / "test_trace_nine_projection.py"
PROOF = ROOT / "proofs" / "trace_nine_projection.md"
PLAN = ROOT / "checkpoints" / "wave26_research_plan.md"
BASE_PRODUCER = ROOT / "experiments" / "e248_replica_trace_nine.py"
BASE_VERIFIER = ROOT / "tests" / "test_replica_trace_nine.py"
BASE_PROOF = ROOT / "proofs" / "replica_trace_nine.md"
BASE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
LOCKFILE = ROOT / "uv.lock"
CPU_LIMIT_SECONDS = 840.0
RSS_LIMIT_BYTES = 7 * 1024**3 // 4
Q, U, X = sp.symbols("q u x")
A7, A6, A5, A4, A3 = sp.symbols("a7 a6 a5 a4 a3")
VARIABLES = (A7, A6, A5, A4, A3)
S0, S2, S3, S4, S5, S6, S7, S8, S9 = sp.symbols("s0 s2 s3 s4 s5 s6 s7 s8 s9")
SYMBOLS = (S0, S2, S3, S4, S5, S6, S7, S8, S9)
TARGET_NAMES = (
    "r1_squared",
    "r2",
    "r3_over_r1",
    "r4",
    "r5_over_r1",
    "r6_over_r2",
    "r7_over_r1",
    "r8",
    "r9_over_r1",
)
Monomial = tuple[int, int, int, int, int]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        library = ctypes.CDLL("/usr/lib/libSystem.B.dylib")

        class Info(ctypes.Structure):
            _fields_ = [
                ("virtual", ctypes.c_uint64),
                ("resident", ctypes.c_uint64),
                ("resident_max", ctypes.c_uint64),
                ("u_s", ctypes.c_int32),
                ("u_us", ctypes.c_int32),
                ("s_s", ctypes.c_int32),
                ("s_us", ctypes.c_int32),
                ("policy", ctypes.c_int32),
                ("suspend", ctypes.c_int32),
            ]

        library.mach_task_self.restype = ctypes.c_uint32
        library.task_info.argtypes = [
            ctypes.c_uint32,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        library.task_info.restype = ctypes.c_int
        info = Info()
        count = ctypes.c_uint32(ctypes.sizeof(info) // 4)
        result = library.task_info(library.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count))
        if result:
            raise RuntimeError(f"task_info failed: {result}")
        return int(info.resident_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def polynomial(coefficients: Iterable[int | str]) -> sp.Expr:
    values = [int(value) for value in coefficients]
    return sp.Poly.from_list(list(reversed(values)), gens=Q, domain=ZZ).as_expr()


def target_expressions(rows: list[list[int]]) -> dict[str, sp.Expr]:
    traces = [polynomial(row) for row in rows]
    delta = Q**2 - 1
    scales = {power: Q ** (-10 * power) / delta ** (8 * power) for power in range(1, 5)}
    return {
        "r1_squared": sp.cancel(traces[0] ** 2 * scales[1]),
        "r2": sp.cancel(traces[1] * scales[1]),
        "r3_over_r1": sp.cancel(traces[2] * scales[1] / traces[0]),
        "r4": sp.cancel(traces[3] * scales[2]),
        "r5_over_r1": sp.cancel(traces[4] * scales[2] / traces[0]),
        "r6_over_r2": sp.cancel(traces[5] * scales[2] / traces[1]),
        "r7_over_r1": sp.cancel(traces[6] * scales[3] / traces[0]),
        "r8": sp.cancel(traces[7] * scales[4]),
        "r9_over_r1": sp.cancel(traces[8] * scales[4] / traces[0]),
    }


def target_record(expression: sp.Expr) -> dict[str, object]:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    pn = sp.Poly(numerator, Q, domain=QQ)
    pd = sp.Poly(denominator, Q, domain=QQ)
    ascending_n = [str(value) for value in reversed(pn.all_coeffs())]
    ascending_d = [str(value) for value in reversed(pd.all_coeffs())]
    return {
        "numerator_degree": pn.degree(),
        "denominator_degree": pd.degree(),
        "numerator_term_count": len(pn.terms()),
        "denominator_term_count": len(pd.terms()),
        "sha256": canonical_digest({"numerator": ascending_n, "denominator": ascending_d}),
    }


def independent_lucas_factors() -> tuple[tuple[str, sp.Expr], ...]:
    tau = sp.symbols("tau")
    sequence = [sp.Integer(2), tau]
    while len(sequence) <= 9:
        sequence.append(sp.expand(tau * sequence[-1] - sequence[-2]))

    def even_quotient(index: int, divisor: sp.Expr) -> sp.Expr:
        quotient, remainder = sp.div(sequence[index], divisor, tau, domain=ZZ)
        if remainder:
            raise AssertionError(f"L{index} divisor failed")
        output = 0
        for (power,), coefficient in sp.Poly(quotient, tau, domain=ZZ).terms():
            if power & 1:
                raise AssertionError(f"L{index} quotient is odd")
            output += coefficient * U ** (power // 2)
        return sp.expand(output)

    return (
        ("r4", even_quotient(4, 1)),
        ("r5_over_r1", even_quotient(5, tau)),
        ("r6_over_r2", even_quotient(6, sequence[2])),
        ("r7_over_r1", even_quotient(7, tau)),
        ("r8", even_quotient(8, 1)),
        ("r9_over_r1", even_quotient(9, tau)),
    )


def independent_incidence() -> tuple[list[sp.Poly], tuple[tuple[str, sp.Expr], ...]]:
    c2 = S2 - S0 - 256 - 128 * A7 - 64 * A6 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = S3 - S0 - 6561 - 2187 * A7 - 729 * A6 - 243 * A5 - 81 * A4 - 27 * A3
    a2 = c3 / 3 - c2 / 2
    a1 = 3 * c2 / 2 - 2 * c3 / 3
    mode = (
        U**8 + A7 * U**7 + A6 * U**6 + A5 * U**5 + A4 * U**4 + A3 * U**3
        + a2 * U**2 + a1 * U + S0
    )
    factors = independent_lucas_factors()
    targets = dict(zip(TARGET_NAMES[3:], SYMBOLS[3:]))
    field = QQ.frac_field(*SYMBOLS)
    equations = [
        sp.Poly(sp.resultant(factor, mode, U) - targets[name], *VARIABLES, domain=field)
        for name, factor in factors
    ]
    return equations, factors


def leading_quotient(equations: list[sp.Poly]) -> dict[str, Any]:
    forms = []
    for equation in equations[:5]:
        degree = equation.total_degree()
        expression = sum(
            coefficient.as_expr()
            * math.prod(variable**power for variable, power in zip(VARIABLES, monomial))
            for monomial, coefficient in equation.terms()
            if sum(monomial) == degree
        )
        if set(expression.free_symbols) & set(SYMBOLS):
            raise AssertionError("leading form retained a target")
        forms.append(sp.Poly(expression, *VARIABLES, domain=QQ))
    basis = sp.groebner([form.as_expr() for form in forms], *VARIABLES, order="grevlex", domain=QQ)
    leading = [polynomial.LM(order=basis.order).exponents for polynomial in basis.polys]
    bounds = []
    for index in range(5):
        pure = [m[index] for m in leading if m[index] and sum(m) == m[index]]
        bounds.append(min(pure))
    standard = [
        exponents
        for exponents in product(*(range(bound) for bound in bounds))
        if not any(all(left <= right for left, right in zip(lm, exponents)) for lm in leading)
    ]
    return {"forms": forms, "basis": basis, "leading": leading, "bounds": bounds, "standard": standard}


def compile_rows(equations: list[sp.Poly]) -> list[list[tuple[Monomial, list[tuple[tuple[int, ...], sp.Rational]]]]]:
    output = []
    for equation in equations:
        row = []
        for a_monomial, coefficient in equation.terms():
            coefficient_poly = sp.Poly(coefficient.as_expr(), *SYMBOLS, domain=QQ)
            row.append((a_monomial, coefficient_poly.terms()))
        output.append(row)
    return output


def horner_mod(coefficients: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def physical_targets_mod(
    rows: list[list[int]], q_value: Fraction, prime: int
) -> tuple[tuple[int, ...], dict[str, int]]:
    q_mod = q_value.numerator * pow(q_value.denominator, -1, prime) % prime
    if q_mod == 0:
        raise AssertionError(f"q={q_value} vanishes modulo p={prime}")
    delta = (q_mod * q_mod - 1) % prime
    traces = [horner_mod(row, q_mod, prime) for row in rows]
    denominators = {
        "q_denominator": q_value.denominator % prime,
        "delta": delta,
        "H1": traces[0],
        "H2": traces[1],
    }
    if any(value == 0 for value in denominators.values()):
        raise AssertionError(f"bad independent specialization {q_value}, {prime}")

    def scale(power: int) -> int:
        numerator = pow(q_mod, -10 * power, prime)
        denominator = pow(delta, 8 * power, prime)
        return numerator * pow(denominator, -1, prime) % prime

    values = (
        traces[0] * traces[0] % prime * scale(1) % prime,
        traces[1] * scale(1) % prime,
        traces[2] * scale(1) % prime * pow(traces[0], -1, prime) % prime,
        traces[3] * scale(2) % prime,
        traces[4] * scale(2) % prime * pow(traces[0], -1, prime) % prime,
        traces[5] * scale(2) % prime * pow(traces[1], -1, prime) % prime,
        traces[6] * scale(3) % prime * pow(traces[0], -1, prime) % prime,
        traces[7] * scale(4) % prime,
        traces[8] * scale(4) % prime * pow(traces[0], -1, prime) % prime,
    )
    return values, {"q_residue": q_mod, **denominators}


def evaluate_target_polynomial(
    terms: list[tuple[tuple[int, ...], sp.Rational]], values: tuple[int, ...], prime: int
) -> int:
    result = 0
    for powers, coefficient in terms:
        term = int(coefficient.p) % prime * pow(int(coefficient.q) % prime, -1, prime) % prime
        for value, power in zip(values, powers):
            term = term * pow(value, power, prime) % prime
        result = (result + term) % prime
    return result


def specialize(
    rows: list[list[tuple[Monomial, list[tuple[tuple[int, ...], sp.Rational]]]]],
    values: tuple[int, ...],
    prime: int,
) -> list[sp.Poly]:
    domain = sp.GF(prime)
    output = []
    for row in rows:
        coefficients = {}
        for monomial, terms in row:
            value = evaluate_target_polynomial(terms, values, prime)
            if value:
                coefficients[monomial] = value
        output.append(sp.Poly.from_dict(coefficients, VARIABLES, domain=domain))
    return output


def divides(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    return all(a <= b for a, b in zip(left, right))


def add(left: Monomial, right: Monomial) -> Monomial:
    return tuple(left[index] + right[index] for index in range(5))  # type: ignore[return-value]


def exact_monomials(count: int, degree: int) -> Iterable[tuple[int, ...]]:
    if count == 1:
        yield (degree,)
    else:
        for first in range(degree + 1):
            for tail in exact_monomials(count - 1, degree - first):
                yield (first, *tail)


def rule_lookup(leading: list[Monomial], maximum: int = 11) -> dict[Monomial, int]:
    lookup = {}
    for degree in range(maximum + 1):
        for raw in exact_monomials(5, degree):
            monomial: Monomial = raw  # type: ignore[assignment]
            for index, divisor in enumerate(leading):
                if divides(divisor, monomial):
                    lookup[monomial] = index
                    break
    return lookup


def grevlex_key(monomial: Monomial) -> tuple[int, tuple[int, ...]]:
    return sum(monomial), tuple(-value for value in reversed(monomial))


def independent_reduce(
    polynomial: dict[Monomial, int],
    rules: list[tuple[Monomial, dict[Monomial, int]]],
    lookup: dict[Monomial, int],
    prime: int,
) -> tuple[dict[Monomial, int], int]:
    remainder = {m: c % prime for m, c in polynomial.items() if c % prime}
    steps = 0
    while True:
        candidates = [monomial for monomial in remainder if monomial in lookup]
        if not candidates:
            return remainder, steps
        monomial = max(candidates, key=grevlex_key)
        rule_index = lookup[monomial]
        leading, tail = rules[rule_index]
        coefficient = remainder.pop(monomial)
        shift: Monomial = tuple(monomial[index] - leading[index] for index in range(5))  # type: ignore[assignment]
        for tail_monomial, tail_coefficient in tail.items():
            shifted = add(tail_monomial, shift)
            value = (remainder.get(shifted, 0) - coefficient * tail_coefficient) % prime
            if value:
                remainder[shifted] = value
            elif shifted in remainder:
                del remainder[shifted]
        steps += 1


def scalar_det(matrix: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    determinant = 1
    swaps = 0
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            swaps ^= 1
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        work[column][column:] = [value * inverse % prime for value in work[column][column:]]
        for row in range(column + 1, len(work)):
            factor = work[row][column]
            if factor:
                work[row][column:] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row][column:], work[column][column:])
                ]
    return (-determinant) % prime if swaps else determinant


def scalar_rank(matrix: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row][column:] = [value * inverse % prime for value in work[pivot_row][column:]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row][column:] = [
                (left - factor * right) % prime
                for left, right in zip(work[row][column:], work[pivot_row][column:])
            ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def rebuild_certificate(
    equations: list[sp.Poly],
    prime: int,
    generic_leading: list[Monomial],
    standard: list[Monomial],
    lookup: dict[Monomial, int],
) -> dict[str, object]:
    basis = sp.groebner(
        [equation.as_expr() for equation in equations[:5]],
        *VARIABLES,
        order="grevlex",
        domain=sp.GF(prime),
        method="buchberger",
    )
    leading = [polynomial.LM(order=basis.order).exponents for polynomial in basis.polys]
    if leading != generic_leading:
        raise AssertionError("independent specialized leading ideal changed")
    rules = []
    groebner_record = []
    for polynomial in basis.polys:
        terms = polynomial.terms(order="grevlex")
        leading_monomial, leading_coefficient = terms[0]
        inverse = pow(int(leading_coefficient) % prime, -1, prime)
        normalized = [
            [list(monomial), int(coefficient) % prime * inverse % prime]
            for monomial, coefficient in terms
        ]
        groebner_record.append(normalized)
        rules.append(
            (
                leading_monomial,
                {
                    monomial: int(coefficient) % prime * inverse % prime
                    for monomial, coefficient in terms[1:]
                },
            )
        )
    ninth = {monomial: int(coefficient) % prime for monomial, coefficient in equations[5].terms()}
    basis_index = {monomial: index for index, monomial in enumerate(standard)}
    columns = []
    step_counts = []
    for basis_monomial in standard:
        source = {add(monomial, basis_monomial): coefficient for monomial, coefficient in ninth.items()}
        reduced, steps = independent_reduce(source, rules, lookup, prime)
        if any(monomial not in basis_index for monomial in reduced):
            raise AssertionError("independent reduction left the standard basis")
        columns.append([reduced.get(monomial, 0) for monomial in standard])
        step_counts.append(steps)
    matrix = [[columns[column][row] for column in range(96)] for row in range(96)]
    return {
        "prime": prime,
        "groebner_basis_size": len(basis.polys),
        "groebner_sha256": canonical_digest(groebner_record),
        "leading_monomials_sha256": canonical_digest([list(value) for value in leading]),
        "quotient_basis_rank": len(standard),
        "matrix_shape": [96, 96],
        "matrix_sha256": canonical_digest(matrix),
        "matrix_rank": scalar_rank(matrix, prime),
        "determinant_residue": scalar_det(matrix, prime),
        "reduction_steps": {
            "minimum": min(step_counts),
            "maximum": max(step_counts),
            "sum": sum(step_counts),
        },
    }


def control_data(factors: tuple[tuple[str, sp.Expr], ...]) -> tuple[sp.Poly, tuple[int, ...]]:
    mode = sp.Poly(math.prod(U - root for root in range(4, 12)), U, domain=ZZ)
    values = [int(mode.eval(point)) for point in (0, 2, 3)]
    values.extend(int(sp.resultant(factor, mode.as_expr(), U)) for _name, factor in factors)
    return mode, tuple(values)


def primitive(expression: sp.Expr) -> sp.Poly:
    polynomial = sp.Poly(expression, Q, domain=QQ)
    _scale, integral = polynomial.clear_denoms(convert=True)
    _content, result = integral.primitive()
    return result


def shift_certificate(expression: sp.Expr) -> dict[str, object]:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    if sp.Poly(denominator, Q, domain=QQ).eval(2) < 0:
        numerator, denominator = -numerator, -denominator
    pn, pd = primitive(numerator), primitive(denominator)
    shifted_n = sp.Poly(sp.expand(pn.as_expr().subs(Q, X + 1)), X, domain=ZZ)
    shifted_d = sp.Poly(sp.expand(pd.as_expr().subs(Q, X + 1)), X, domain=ZZ)
    cn = [int(value) for value in reversed(shifted_n.all_coeffs())]
    cd = [int(value) for value in reversed(shifted_d.all_coeffs())]
    if any(value <= 0 for value in cn) or any(value < 0 for value in cd) or not any(cd):
        raise AssertionError("independent shifted positivity failed")
    return {
        "numerator_degree": pn.degree(),
        "denominator_degree": pd.degree(),
        "shifted_numerator_term_count": len(shifted_n.terms()),
        "shifted_denominator_term_count": len(shifted_d.terms()),
        "shifted_numerator_minimum_coefficient": str(min(cn)),
        "shifted_denominator_zero_coefficients": sum(value == 0 for value in cd),
        "shifted_numerator_sha256": canonical_digest([str(value) for value in cn]),
        "shifted_denominator_sha256": canonical_digest([str(value) for value in cd]),
    }


def gate_data(
    targets: dict[str, sp.Expr], factors: tuple[tuple[str, sp.Expr], ...]
) -> dict[str, object]:
    shifted = {
        name: sp.Poly(sp.expand(factor.subs(U, X + 4)), X, domain=ZZ)
        for name, factor in factors
    }
    f4, f5, f6, f7, f8, f9 = [shifted[name].as_expr() for name in TARGET_NAMES[3:]]
    differences = {
        "f4_minus_f6": f4 - f6,
        "f7_minus_f5": f7 - f5,
        "f8_minus_f9": f8 - f9,
        "f4_squared_minus_f8": f4**2 - f8,
        "f4_f6_minus_f5_squared": f4 * f6 - f5**2,
        "f5_f7_minus_f6_squared": f5 * f7 - f6**2,
        "f6_f8_minus_f7_squared": f6 * f8 - f7**2,
    }
    pointwise = {}
    for name, expression in differences.items():
        coefficients = [int(value) for value in reversed(sp.Poly(sp.expand(expression), X, domain=ZZ).all_coeffs())]
        if any(value < 0 for value in coefficients) or not any(coefficients):
            raise AssertionError(f"independent pointwise gate failed: {name}")
        pointwise[name] = coefficients
    r4, r5, r6, r7, r8, r9 = [targets[name] for name in TARGET_NAMES[3:]]
    physical = {
        "r4_gt_r6": r4 - r6,
        "r7_gt_r5": r7 - r5,
        "r8_gt_r9": r8 - r9,
        "r4_squared_gt_r8": r4**2 - r8,
        "r4_r6_gt_r5_squared": r4 * r6 - r5**2,
        "r5_r7_ge_r6_squared": r5 * r7 - r6**2,
        "r6_r8_gt_r7_squared": r6 * r8 - r7**2,
    }
    return {
        "shifted_factor_coefficients_ascending": {
            name: [int(value) for value in reversed(poly.all_coeffs())]
            for name, poly in shifted.items()
        },
        "pointwise_difference_coefficients_ascending": pointwise,
        "physical_target_certificates": {
            name: shift_certificate(sp.cancel(expression)) for name, expression in physical.items()
        },
        "all_strict_for_literal_open_2x4_targets": True,
    }


def projection_record(
    targets: dict[str, sp.Expr], equations: list[sp.Poly]
) -> dict[str, object]:
    summaries = {name: target_record(expression) for name, expression in targets.items()}
    first = TARGET_NAMES[:3]
    d_degrees = [int(summaries[name]["denominator_degree"]) for name in first]
    d_total = sum(d_degrees)
    dq_degree = max(
        int(summaries[name]["numerator_degree"]) + d_total - int(summaries[name]["denominator_degree"])
        for name in first
    )
    clearing_degrees = [2, 2, 2, 3, 4, 4]
    row_summaries = [summaries[name] for name in TARGET_NAMES[3:]]
    row_bounds = [
        max(
            dq_degree * degree + int(summary["denominator_degree"]),
            d_total * degree + int(summary["numerator_degree"]),
        )
        for degree, summary in zip(clearing_degrees, row_summaries)
    ]
    equation_degrees = [equation.total_degree() for equation in equations]
    weights = [
        math.prod(equation_degrees[index] for index in range(6) if index != excluded)
        for excluded in range(6)
    ]
    return {
        "first_target_numerator_degrees": [int(summaries[name]["numerator_degree"]) for name in first],
        "first_target_denominator_degrees": d_degrees,
        "common_D_degree": d_total,
        "DQ_coefficient_degree_bound": dq_degree,
        "row_target_numerator_denominator_degrees": [
            [int(summary["numerator_degree"]), int(summary["denominator_degree"])]
            for summary in row_summaries
        ],
        "lucas_resultant_degrees_for_clearing": clearing_degrees,
        "cleared_row_q_degree_bounds": row_bounds,
        "affine_equation_degrees": equation_degrees,
        "projective_resultant_multidegrees": weights,
        "univariate_degree_bound": sum(a * b for a, b in zip(row_bounds, weights)),
        "denominator_scope": (
            "the cleared resultant may contain denominator factors, but all inherited target "
            "denominators are strictly positive for physical q>1"
        ),
    }


def semantic_certificate(record: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in record.items() if key != "process_cpu_seconds"}


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    base = json.loads(BASE_ARTIFACT.read_text())
    failures: list[str] = []

    expected_sources = {
        str(path.relative_to(ROOT.parent.parent)): digest(path)
        for path in (
            PRODUCER,
            SELF,
            PROOF,
            PLAN,
            BASE_PRODUCER,
            BASE_VERIFIER,
            BASE_PROOF,
            BASE_ARTIFACT,
            LOCKFILE,
        )
    }
    if artifact["meta"]["source_sha256"] != expected_sources:
        failures.append("complete source hash map mismatch")
    if artifact["meta"]["data_sha256"] != canonical_digest(artifact["data"]):
        failures.append("complete data hash mismatch")
    if not all(check["passed"] for check in base["checks"]):
        failures.append("inherited e248 artifact contains a failed check")
    for relative, expected in base["meta"]["source_sha256"].items():
        if digest(ROOT / relative) != expected:
            failures.append(f"inherited source drift: {relative}")

    rows = [[int(value) for value in row] for row in base["data"]["trace_coefficients_ascending"]]
    trace_hash = canonical_digest([[str(value) for value in row] for row in rows])
    targets = target_expressions(rows)
    summaries = {name: target_record(value) for name, value in targets.items()}
    if summaries != base["data"]["normalized_targets"]:
        failures.append("independent normalized target summary mismatch")
    expected_inherited = {
        "artifact_sha256": digest(BASE_ARTIFACT),
        "trace_coefficients_sha256": trace_hash,
        "trace_count": 9,
        "base_check_names": [check["name"] for check in base["checks"]],
        "normalized_target_sha256": canonical_digest(summaries),
    }
    if artifact["data"]["inherited_e248"] != expected_inherited:
        failures.append("inherited e248 semantic record mismatch")

    equations, factors = independent_incidence()
    quotient = leading_quotient(equations)
    leading = quotient["leading"]
    standard = quotient["standard"]
    hilbert = [sum(sum(monomial) == degree for monomial in standard) for degree in range(9)]
    expected_system = {
        "coefficient_variables": [str(value) for value in VARIABLES],
        "equation_names": [f"F{index}" for index in range(4, 10)],
        "equation_total_degrees": [2, 2, 2, 3, 4, 3],
        "leading_form_term_counts": [len(form.terms()) for form in quotient["forms"]],
        "groebner_basis_size": len(quotient["basis"].polys),
        "leading_monomials_sha256": canonical_digest([list(value) for value in leading]),
        "pure_power_bounds": quotient["bounds"],
        "quotient_basis_rank": len(standard),
        "standard_monomials_sha256": canonical_digest([list(value) for value in standard]),
        "maximum_standard_degree": max(map(sum, standard)),
        "hilbert_vector_degrees_0_through_8": hilbert,
        "reduction_monomials_through_degree_11": sum(
            1 for degree in range(12) for _ in exact_monomials(5, degree)
        ),
    }
    if expected_system != artifact["data"]["coefficient_system"]:
        failures.append("coefficient quotient/Hilbert record mismatch")
    if hilbert != [1, 5, 12, 19, 22, 19, 12, 5, 1] or len(standard) != 96:
        failures.append("independent quotient dimension invariant failed")

    compiled = compile_rows(equations)
    lookup = rule_lookup(leading)
    stored_physical = artifact["data"]["modular_certificates"]["physical_witnesses"]
    specifications = (("2", Fraction(2), 2_147_483_647), ("5/3", Fraction(5, 3), 2_147_483_629))
    if [(entry["q"], entry["prime"]) for entry in stored_physical] != [
        (label, prime) for label, _value, prime in specifications
    ]:
        failures.append("physical witness coverage/order mismatch")
    rebuilt_physical = []
    for label, q_value, prime in specifications:
        if not sp.isprime(prime) or prime * prime >= 1 << 63:
            failures.append(f"unsafe claimed prime {prime}")
            continue
        values, denominators = physical_targets_mod(rows, q_value, prime)
        specialized = specialize(compiled, values, prime)
        certificate = rebuild_certificate(specialized, prime, leading, standard, lookup)
        certificate.update(
            {
                "q": label,
                "q_residue": denominators.pop("q_residue"),
                "target_residues": list(values),
                "target_residues_sha256": canonical_digest(list(values)),
                "nonzero_denominator_residues": denominators,
            }
        )
        rebuilt_physical.append(certificate)
    for rebuilt, stored in zip(rebuilt_physical, stored_physical):
        if semantic_certificate(stored) != rebuilt:
            failures.append(f"physical witness mismatch at q={stored['q']}")
        if not 0 < float(stored["process_cpu_seconds"]) < CPU_LIMIT_SECONDS:
            failures.append(f"invalid stored physical stage CPU at q={stored['q']}")
        if not rebuilt["determinant_residue"] or rebuilt["matrix_rank"] != 96:
            failures.append(f"physical witness is not invertible at q={rebuilt['q']}")
    if [entry["determinant_residue"] for entry in rebuilt_physical] != [1_660_951_362, 1_100_728_375]:
        failures.append("frozen determinant residues changed")

    mode, integer_targets = control_data(factors)
    control_prime = 2_147_483_647
    control_values = tuple(value % control_prime for value in integer_targets)
    control_equations = specialize(compiled, control_values, control_prime)
    point = [int(value) % control_prime for value in mode.all_coeffs()[1:6]]
    residuals = [int(equation.eval(dict(zip(VARIABLES, point)))) % control_prime for equation in control_equations]
    rebuilt_control = rebuild_certificate(control_equations, control_prime, leading, standard, lookup)
    rebuilt_control.update(
        {
            "mode_roots": list(range(4, 12)),
            "mode_polynomial_coefficients_descending": [str(value) for value in mode.all_coeffs()],
            "target_integers": [str(value) for value in integer_targets],
            "target_residues": list(control_values),
            "coefficient_point_residues": point,
            "point_residuals": residuals,
        }
    )
    stored_control = artifact["data"]["modular_certificates"]["synthetic_control"]
    if semantic_certificate(stored_control) != rebuilt_control:
        failures.append("synthetic eight-mode control mismatch")
    if any(residuals) or rebuilt_control["determinant_residue"] or rebuilt_control["matrix_rank"] >= 96:
        failures.append("synthetic eight-mode singularity control failed")
    if not 0 < float(stored_control["process_cpu_seconds"]) < CPU_LIMIT_SECONDS:
        failures.append("invalid stored control stage CPU")

    degree = projection_record(targets, equations)
    degree["nonzero_specialization"] = {
        "q": "2",
        "prime": 2_147_483_647,
        "determinant_residue": 1_660_951_362,
    }
    if degree != artifact["data"]["projection_degree"]:
        failures.append("projection degree derivation mismatch")
    if degree["cleared_row_q_degree_bounds"] != [192, 200, 224, 296, 384, 392]:
        failures.append("cleared row degree vector mismatch")
    if degree["projective_resultant_multidegrees"] != [144, 144, 144, 96, 72, 96]:
        failures.append("resultant multidegree vector mismatch")
    if degree["univariate_degree_bound"] != 182_400:
        failures.append("univariate projection degree mismatch")

    gates = gate_data(targets, factors)
    if gates != artifact["data"]["lucas_product_gates"]:
        failures.append("Lucas product gate record mismatch")
    if set(gates["physical_target_certificates"]) != {
        "r4_gt_r6",
        "r7_gt_r5",
        "r8_gt_r9",
        "r4_squared_gt_r8",
        "r4_r6_gt_r5_squared",
        "r5_r7_ge_r6_squared",
        "r6_r8_gt_r7_squared",
    }:
        failures.append("Lucas product gate coverage mismatch")

    expected_theorem = {
        "statement": (
            "For the open 2x4 layer, the trace-four-through-nine full-eight-mode incidence is "
            "supported over a nonzero univariate projection polynomial in q of degree at most "
            "182400. Hence it can occur at no more than 182400 complex q values at which the "
            "normalized targets are defined, and at only finitely many physical q>1 couplings."
        ),
        "exact_witness": (
            "At q=2, reduction modulo p=2147483647 gives a 96-dimensional F4..F8 quotient "
            "and an invertible multiplication-by-F9 matrix. The resulting absence of affine "
            "and projective roots proves the cleared characteristic-zero resultant nonzero."
        ),
        "emptiness_status": "UNRESOLVED",
        "unresolved": (
            "The primitive F9 norm is not reconstructed; its q>1 roots are not isolated; no "
            "coefficient branch has been subjected to complete shifted-Stieltjes disposition."
        ),
        "scope": (
            "exact trace-four-through-nine necessary incidence for a positive full eight-mode "
            "subset-product spectrum of the finite open 2x4 transfer matrix; this is not an "
            "exceptional-root list, an emptiness theorem, an all-size result, or a thermodynamic solution"
        ),
    }
    if artifact["data"]["theorem"] != expected_theorem:
        failures.append("theorem/scope semantic record mismatch")
    expected_resource = {
        "full_norm_status": "not materialized",
        "coarse_interpolation_nodes": 182_401,
        "direct_generic_interpolation": (
            "forbidden by the frozen resource gate; every node would require a specialized "
            "five-variable quotient without a compiled lifted template"
        ),
        "next_exact_action": (
            "construct a fixed-leading-form lift or a fraction-free border template, then "
            "derive a tighter denominator/degree/height envelope before interpolation"
        ),
    }
    if artifact["data"]["resource_boundary"] != expected_resource:
        failures.append("resource boundary semantic record mismatch")
    if artifact["data"]["graph"] != {
        "name": "open 2x4 grid",
        "sites": 8,
        "in_layer_edges": 10,
        "transfer_dimension": 256,
        "physical_parameter": "q=(1+t^2)/(2t)>1",
    }:
        failures.append("graph/parameter record mismatch")
    if artifact["data"]["claim_tag"] != "[THEOREM][FINITE TRACE-NINE PROJECTION][EXACT MODULAR WITNESS]":
        failures.append("claim tag mismatch")

    expected_check_names = {
        "C0_inherited_e248_provenance",
        "C1_abstract_incidence_and_hilbert",
        "C2_two_physical_modular_witnesses",
        "C3_synthetic_eight_mode_control",
        "C4_nonzero_projection_degree_bound",
        "C5_pointwise_lucas_product_gates",
        "C6_literal_target_gate_positivity",
        "C7_theorem_scope_gate",
        "C8_cpu_budget",
        "C9_rss_budget",
    }
    if {check["name"] for check in artifact["checks"]} != expected_check_names:
        failures.append("producer check-name coverage mismatch")
    if not all(check["passed"] for check in artifact["checks"]):
        failures.append("producer artifact contains a failed check")

    elapsed = time.process_time() - started
    rss = peak_rss_bytes()
    if elapsed >= CPU_LIMIT_SECONDS:
        failures.append(f"verifier CPU {elapsed:.6f}s reached {CPU_LIMIT_SECONDS}s")
    if rss >= RSS_LIMIT_BYTES:
        failures.append(f"verifier RSS {rss} reached {RSS_LIMIT_BYTES}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print(
        "PASS: e251 finite trace-nine projection; "
        f"two nonzero physical witnesses, degree<=182400; CPU={elapsed:.6f}s RSS={rss} bytes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
