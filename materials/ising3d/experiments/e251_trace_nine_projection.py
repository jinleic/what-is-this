#!/usr/bin/env python3
"""Exact finite trace-nine projection theorem for the open 2x4 Ising layer.

The script consumes e248's landed trace coefficients, rebuilds the abstract
Lucas-resultant incidence, and specializes the 96-dimensional F4..F8 quotient
at exact finite-field points.  A nonzero determinant of multiplication by F9
is a characteristic-zero nonvanishing witness, not a numerical approximation.
"""
from __future__ import annotations

import ctypes
import hashlib
import heapq
import json
import math
import os
import platform
import resource
import time
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ, ZZ

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e251_trace_nine_projection.py"
OUTPUT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
VERIFIER = ROOT / "tests" / "test_trace_nine_projection.py"
PROOF = ROOT / "proofs" / "trace_nine_projection.md"
PLAN = ROOT / "checkpoints" / "wave26_research_plan.md"
BASE_PRODUCER = ROOT / "experiments" / "e248_replica_trace_nine.py"
BASE_VERIFIER = ROOT / "tests" / "test_replica_trace_nine.py"
BASE_PROOF = ROOT / "proofs" / "replica_trace_nine.md"
BASE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
LOCKFILE = ROOT / "uv.lock"
CPU_LIMIT_SECONDS = 840.0
RSS_LIMIT_BYTES = 7 * 1024**3 // 4
PRIMES = (2_147_483_647, 2_147_483_629)
QPAR, U = sp.symbols("q u")
X = sp.symbols("x")
A7, A6, A5, A4, A3 = sp.symbols("a7 a6 a5 a4 a3")
COEFFICIENT_VARIABLES = (A7, A6, A5, A4, A3)
S0, S2, S3, S4, S5, S6, S7, S8, S9 = sp.symbols("s0 s2 s3 s4 s5 s6 s7 s8 s9")
TARGET_SYMBOLS = (S0, S2, S3, S4, S5, S6, S7, S8, S9)
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
ROW_TARGET_NAMES = TARGET_NAMES[3:]

Monomial = tuple[int, int, int, int, int]
CompiledCoefficient = tuple[tuple[tuple[int, ...], sp.Rational], ...]
CompiledRow = tuple[tuple[Monomial, CompiledCoefficient], ...]


class ResourceWall(RuntimeError):
    """Raised before the declared process resource ceiling is crossed."""


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


def peak_rss_measurement() -> str:
    if platform.system() == "Darwin":
        return "mach_task_basic_info.resident_size_max (task_info flavor 20)"
    return "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"


def guard_resources(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    peak = peak_rss_bytes()
    if elapsed >= CPU_LIMIT_SECONDS:
        raise ResourceWall(f"{stage}: process CPU {elapsed:.6f}s reached {CPU_LIMIT_SECONDS}s")
    if peak >= RSS_LIMIT_BYTES:
        raise ResourceWall(f"{stage}: peak RSS {peak} reached {RSS_LIMIT_BYTES} bytes")


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def polynomial_from_ascending(coefficients: Iterable[int | str], variable: sp.Symbol) -> sp.Expr:
    values = [int(value) for value in coefficients]
    return sp.Poly.from_list(list(reversed(values)), gens=variable, domain=ZZ).as_expr()


def normalized_targets(trace_rows: list[list[int]]) -> dict[str, sp.Expr]:
    traces = [polynomial_from_ascending(row, QPAR) for row in trace_rows]
    delta = QPAR**2 - 1

    def scale(center_power: int) -> sp.Expr:
        return QPAR ** (-10 * center_power) / delta ** (8 * center_power)

    return {
        "r1_squared": sp.cancel(traces[0] ** 2 * scale(1)),
        "r2": sp.cancel(traces[1] * scale(1)),
        "r3_over_r1": sp.cancel(traces[2] * scale(1) / traces[0]),
        "r4": sp.cancel(traces[3] * scale(2)),
        "r5_over_r1": sp.cancel(traces[4] * scale(2) / traces[0]),
        "r6_over_r2": sp.cancel(traces[5] * scale(2) / traces[1]),
        "r7_over_r1": sp.cancel(traces[6] * scale(3) / traces[0]),
        "r8": sp.cancel(traces[7] * scale(4)),
        "r9_over_r1": sp.cancel(traces[8] * scale(4) / traces[0]),
    }


def rational_summary(value: sp.Expr) -> dict[str, object]:
    numerator, denominator = sp.cancel(value).as_numer_denom()
    numerator_poly = sp.Poly(numerator, QPAR, domain=QQ)
    denominator_poly = sp.Poly(denominator, QPAR, domain=QQ)
    numerator_coefficients = [str(value) for value in reversed(numerator_poly.all_coeffs())]
    denominator_coefficients = [str(value) for value in reversed(denominator_poly.all_coeffs())]
    return {
        "numerator_degree": int(numerator_poly.degree()),
        "denominator_degree": int(denominator_poly.degree()),
        "numerator_term_count": len(numerator_poly.terms()),
        "denominator_term_count": len(denominator_poly.terms()),
        "sha256": canonical_sha256(
            {"numerator": numerator_coefficients, "denominator": denominator_coefficients}
        ),
    }


def lucas_factors() -> tuple[tuple[str, sp.Expr], ...]:
    tau = sp.symbols("tau")
    lucas = [sp.Integer(2), tau]
    for _index in range(2, 10):
        lucas.append(sp.expand(tau * lucas[-1] - lucas[-2]))

    def quotient_in_u(index: int, divisor: sp.Expr) -> sp.Expr:
        quotient, remainder = sp.div(lucas[index], divisor, tau, domain=ZZ)
        if remainder:
            raise AssertionError(f"Lucas L_{index} has a nonzero requested remainder")
        output = sp.Integer(0)
        for (power,), coefficient in sp.Poly(quotient, tau, domain=ZZ).terms():
            if power % 2:
                raise AssertionError(f"Lucas L_{index} quotient is not even")
            output += coefficient * U ** (power // 2)
        return sp.expand(output)

    return (
        ("r4", quotient_in_u(4, sp.Integer(1))),
        ("r5_over_r1", quotient_in_u(5, tau)),
        ("r6_over_r2", quotient_in_u(6, lucas[2])),
        ("r7_over_r1", quotient_in_u(7, tau)),
        ("r8", quotient_in_u(8, sp.Integer(1))),
        ("r9_over_r1", quotient_in_u(9, tau)),
    )


def small_quotient_norm(factor: sp.Expr, value: sp.Expr) -> sp.Expr:
    """Return Res(factor,value) as multiplication by value mod monic factor."""
    factor_poly = sp.Poly(factor, U, domain=ZZ)
    degree = factor_poly.degree()
    if factor_poly.LC() != 1 or not 1 <= degree <= 4:
        raise ValueError("small quotient norm requires a monic factor of degree one through four")
    columns: list[list[sp.Expr]] = []
    for power in range(degree):
        remainder = sp.rem(sp.expand(value * U**power), factor, U)
        remainder_poly = sp.Poly(remainder, U, domain="EX")
        columns.append([remainder_poly.nth(row) for row in range(degree)])
    matrix = sp.Matrix(degree, degree, lambda row, column: columns[column][row])
    return sp.expand(matrix.det())


def abstract_mode_equations() -> tuple[list[sp.Poly], tuple[tuple[str, sp.Expr], ...]]:
    c2 = S2 - S0 - 256 - 128 * A7 - 64 * A6 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = S3 - S0 - 6561 - 2187 * A7 - 729 * A6 - 243 * A5 - 81 * A4 - 27 * A3
    a2 = c3 / 3 - c2 / 2
    a1 = 3 * c2 / 2 - 2 * c3 / 3
    mode = (
        U**8
        + A7 * U**7
        + A6 * U**6
        + A5 * U**5
        + A4 * U**4
        + A3 * U**3
        + a2 * U**2
        + a1 * U
        + S0
    )
    factors = lucas_factors()
    target_by_name = dict(zip(ROW_TARGET_NAMES, TARGET_SYMBOLS[3:]))
    field = QQ.frac_field(*TARGET_SYMBOLS)
    equations = [
        sp.Poly(
            small_quotient_norm(factor, mode) - target_by_name[name],
            *COEFFICIENT_VARIABLES,
            domain=field,
        )
        for name, factor in factors
    ]
    return equations, factors


def compile_equations(equations: list[sp.Poly]) -> tuple[CompiledRow, ...]:
    rows: list[CompiledRow] = []
    for equation in equations:
        compiled: list[tuple[Monomial, CompiledCoefficient]] = []
        for monomial, coefficient in equation.terms():
            polynomial = sp.Poly(coefficient.as_expr(), *TARGET_SYMBOLS, domain=QQ)
            compiled.append((monomial, tuple(polynomial.terms())))
        rows.append(tuple(compiled))
    return tuple(rows)


def generic_quotient(equations: list[sp.Poly]) -> dict[str, Any]:
    leading_forms: list[sp.Poly] = []
    for equation in equations[:5]:
        degree = equation.total_degree()
        expression = sum(
            coefficient.as_expr()
            * math.prod(variable**power for variable, power in zip(COEFFICIENT_VARIABLES, monomial))
            for monomial, coefficient in equation.terms()
            if sum(monomial) == degree
        )
        if set(expression.free_symbols) & set(TARGET_SYMBOLS):
            raise AssertionError("leading form depends on trace targets")
        leading_forms.append(sp.Poly(expression, *COEFFICIENT_VARIABLES, domain=QQ))
    basis = sp.groebner(
        [form.as_expr() for form in leading_forms],
        *COEFFICIENT_VARIABLES,
        order="grevlex",
        domain=QQ,
    )
    leading_monomials = [polynomial.LM(order=basis.order).exponents for polynomial in basis.polys]
    pure_power_bounds: list[int] = []
    for variable_index in range(5):
        powers = [
            monomial[variable_index]
            for monomial in leading_monomials
            if monomial[variable_index] and sum(monomial) == monomial[variable_index]
        ]
        if not powers:
            raise AssertionError(f"missing pure leading power for variable {variable_index}")
        pure_power_bounds.append(min(powers))
    standard_monomials = [
        monomial
        for monomial in product(*(range(bound) for bound in pure_power_bounds))
        if not any(divides(leading, monomial) for leading in leading_monomials)
    ]
    hilbert = [sum(sum(monomial) == degree for monomial in standard_monomials) for degree in range(9)]
    return {
        "leading_forms": leading_forms,
        "basis": basis,
        "leading_monomials": leading_monomials,
        "pure_power_bounds": pure_power_bounds,
        "standard_monomials": standard_monomials,
        "hilbert_vector": hilbert,
    }


def divides(divisor: tuple[int, ...], multiple: tuple[int, ...]) -> bool:
    return all(left <= right for left, right in zip(divisor, multiple))


def add_monomials(left: Monomial, right: Monomial) -> Monomial:
    return tuple(left[index] + right[index] for index in range(5))  # type: ignore[return-value]


def exact_monomials(variable_count: int, degree: int) -> Iterable[tuple[int, ...]]:
    if variable_count == 1:
        yield (degree,)
        return
    for first in range(degree + 1):
        for tail in exact_monomials(variable_count - 1, degree - first):
            yield (first, *tail)


def monomials_through_degree(variable_count: int, maximum: int) -> Iterable[tuple[int, ...]]:
    for degree in range(maximum + 1):
        yield from exact_monomials(variable_count, degree)


def reduction_priority(monomial: Monomial) -> tuple[int, tuple[int, ...]]:
    return -sum(monomial), tuple(reversed(monomial))


def reduction_rule_map(
    leading_monomials: list[Monomial], maximum_degree: int
) -> dict[Monomial, int]:
    result: dict[Monomial, int] = {}
    for raw in monomials_through_degree(5, maximum_degree):
        monomial: Monomial = raw  # type: ignore[assignment]
        for index, leading in enumerate(leading_monomials):
            if divides(leading, monomial):
                result[monomial] = index
                break
    return result


def reduce_mod_prime(
    polynomial: dict[Monomial, int],
    rules: list[tuple[Monomial, dict[Monomial, int]]],
    rule_map: dict[Monomial, int],
    prime: int,
) -> tuple[dict[Monomial, int], int]:
    remainder = {monomial: coefficient % prime for monomial, coefficient in polynomial.items() if coefficient % prime}
    heap = [(reduction_priority(monomial), monomial) for monomial in remainder if monomial in rule_map]
    heapq.heapify(heap)
    steps = 0
    while heap:
        _priority, monomial = heapq.heappop(heap)
        if monomial not in remainder or monomial not in rule_map:
            continue
        rule_index = rule_map[monomial]
        leading, tail = rules[rule_index]
        coefficient = remainder.pop(monomial)
        shift: Monomial = tuple(monomial[index] - leading[index] for index in range(5))  # type: ignore[assignment]
        for tail_monomial, tail_coefficient in tail.items():
            shifted = add_monomials(tail_monomial, shift)
            previous = remainder.get(shifted, 0)
            updated = (previous - coefficient * tail_coefficient) % prime
            if updated:
                remainder[shifted] = updated
                if not previous and shifted in rule_map:
                    heapq.heappush(heap, (reduction_priority(shifted), shifted))
            elif previous:
                del remainder[shifted]
        steps += 1
    return remainder, steps


def modular_polynomial_value(coefficients: list[int], value: int, prime: int) -> int:
    accumulator = 0
    for coefficient in reversed(coefficients):
        accumulator = (accumulator * value + coefficient) % prime
    return accumulator


def physical_target_residues(
    trace_rows: list[list[int]], q_value: Fraction, prime: int
) -> tuple[tuple[int, ...], dict[str, int]]:
    q_residue = q_value.numerator * pow(q_value.denominator, -1, prime) % prime
    if not q_residue:
        raise ZeroDivisionError(f"q={q_value} vanishes modulo p={prime}")
    delta = (q_residue * q_residue - 1) % prime
    traces = [modular_polynomial_value(row, q_residue, prime) for row in trace_rows]
    denominators = {"q_denominator": q_value.denominator % prime, "delta": delta, "H1": traces[0], "H2": traces[1]}
    if any(not value for value in denominators.values()):
        raise ZeroDivisionError(f"bad specialization at q={q_value}, p={prime}: {denominators}")

    def scale(center_power: int) -> int:
        return (
            pow(q_residue, -10 * center_power, prime)
            * pow(pow(delta, 8 * center_power, prime), -1, prime)
            % prime
        )

    values = (
        traces[0] ** 2 * scale(1) % prime,
        traces[1] * scale(1) % prime,
        traces[2] * scale(1) * pow(traces[0], -1, prime) % prime,
        traces[3] * scale(2) % prime,
        traces[4] * scale(2) * pow(traces[0], -1, prime) % prime,
        traces[5] * scale(2) * pow(traces[1], -1, prime) % prime,
        traces[6] * scale(3) * pow(traces[0], -1, prime) % prime,
        traces[7] * scale(4) % prime,
        traces[8] * scale(4) * pow(traces[0], -1, prime) % prime,
    )
    return values, {"q_residue": q_residue, **denominators}


def evaluate_compiled_coefficient(
    terms: CompiledCoefficient, targets: tuple[int, ...], prime: int
) -> int:
    output = 0
    for monomial, coefficient in terms:
        value = int(coefficient.p) % prime * pow(int(coefficient.q) % prime, -1, prime) % prime
        for target, power in zip(targets, monomial):
            value = value * pow(target, power, prime) % prime
        output = (output + value) % prime
    return output


def specialize_equations(
    compiled: tuple[CompiledRow, ...], targets: tuple[int, ...], prime: int
) -> list[sp.Poly]:
    domain = sp.GF(prime)
    equations = []
    for row in compiled:
        coefficients = {
            monomial: value
            for monomial, terms in row
            if (value := evaluate_compiled_coefficient(terms, targets, prime))
        }
        equations.append(sp.Poly.from_dict(coefficients, COEFFICIENT_VARIABLES, domain=domain))
    return equations


def determinant_mod_prime(matrix: np.ndarray, prime: int) -> int:
    if prime * prime >= 1 << 63:
        raise ValueError("finite-field product does not fit signed int64")
    work = matrix.astype(np.uint64, copy=True)
    modulus = np.uint64(prime)
    determinant = 1
    odd_swaps = False
    for pivot_index in range(work.shape[0]):
        candidates = np.nonzero(work[pivot_index:, pivot_index])[0]
        if not len(candidates):
            return 0
        pivot_row = pivot_index + int(candidates[0])
        if pivot_row != pivot_index:
            work[[pivot_index, pivot_row]] = work[[pivot_row, pivot_index]]
            odd_swaps = not odd_swaps
        pivot = int(work[pivot_index, pivot_index])
        determinant = determinant * pivot % prime
        inverse = pow(pivot, -1, prime)
        work[pivot_index, pivot_index:] = (
            work[pivot_index, pivot_index:] * np.uint64(inverse)
        ) % modulus
        for row in range(pivot_index + 1, work.shape[0]):
            factor = int(work[row, pivot_index])
            if factor:
                product_row = (
                    work[pivot_index, pivot_index:] * np.uint64(factor)
                ) % modulus
                work[row, pivot_index:] = (
                    work[row, pivot_index:] + modulus - product_row
                ) % modulus
    return (-determinant) % prime if odd_swaps else determinant


def rank_mod_prime(matrix: np.ndarray, prime: int) -> int:
    work = matrix.astype(np.uint64, copy=True)
    modulus = np.uint64(prime)
    pivot_row = 0
    for column in range(work.shape[1]):
        candidates = np.nonzero(work[pivot_row:, column])[0]
        if not len(candidates):
            continue
        selected = pivot_row + int(candidates[0])
        if selected != pivot_row:
            work[[pivot_row, selected]] = work[[selected, pivot_row]]
        inverse = pow(int(work[pivot_row, column]), -1, prime)
        work[pivot_row, column:] = (
            work[pivot_row, column:] * np.uint64(inverse)
        ) % modulus
        for row in range(work.shape[0]):
            if row == pivot_row:
                continue
            factor = int(work[row, column])
            if factor:
                product_row = (work[pivot_row, column:] * np.uint64(factor)) % modulus
                work[row, column:] = (work[row, column:] + modulus - product_row) % modulus
        pivot_row += 1
        if pivot_row == work.shape[0]:
            break
    return pivot_row


def quotient_certificate(
    equations: list[sp.Poly],
    prime: int,
    expected_leading: list[Monomial],
    standard_monomials: list[Monomial],
    rule_map: dict[Monomial, int],
    started: float,
    label: str,
) -> dict[str, object]:
    stage_started = time.process_time()
    basis = sp.groebner(
        [equation.as_expr() for equation in equations[:5]],
        *COEFFICIENT_VARIABLES,
        order="grevlex",
        domain=sp.GF(prime),
        method="buchberger",
    )
    guard_resources(started, f"{label} Groebner basis")
    leading = [polynomial.LM(order=basis.order).exponents for polynomial in basis.polys]
    if leading != expected_leading:
        raise AssertionError(f"{label}: unlucky leading-monomial signature")
    rules: list[tuple[Monomial, dict[Monomial, int]]] = []
    groebner_record: list[list[list[object]]] = []
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
    basis_index = {monomial: index for index, monomial in enumerate(standard_monomials)}
    ninth = {monomial: int(coefficient) % prime for monomial, coefficient in equations[5].terms()}
    columns: list[list[int]] = []
    reduction_steps: list[int] = []
    for basis_monomial in standard_monomials:
        product_polynomial = {
            add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in ninth.items()
        }
        reduced, steps = reduce_mod_prime(product_polynomial, rules, rule_map, prime)
        if any(monomial not in basis_index for monomial in reduced):
            raise AssertionError(f"{label}: reduction escaped the 96 standard monomials")
        columns.append([reduced.get(monomial, 0) for monomial in standard_monomials])
        reduction_steps.append(steps)
    matrix = np.asarray(columns, dtype=np.uint64).T
    determinant = determinant_mod_prime(matrix, prime)
    rank = rank_mod_prime(matrix, prime)
    guard_resources(started, f"{label} determinant")
    return {
        "prime": prime,
        "groebner_basis_size": len(basis.polys),
        "groebner_sha256": canonical_sha256(groebner_record),
        "leading_monomials_sha256": canonical_sha256([list(value) for value in leading]),
        "quotient_basis_rank": len(standard_monomials),
        "matrix_shape": [int(value) for value in matrix.shape],
        "matrix_sha256": canonical_sha256([[int(value) for value in row] for row in matrix]),
        "matrix_rank": rank,
        "determinant_residue": determinant,
        "reduction_steps": {
            "minimum": min(reduction_steps),
            "maximum": max(reduction_steps),
            "sum": sum(reduction_steps),
        },
        "process_cpu_seconds": time.process_time() - stage_started,
    }


def synthetic_mode_targets(factors: tuple[tuple[str, sp.Expr], ...]) -> tuple[sp.Poly, tuple[int, ...]]:
    mode = sp.Poly(math.prod(U - root for root in range(4, 12)), U, domain=ZZ)
    targets = [int(mode.eval(value)) for value in (0, 2, 3)]
    targets.extend(int(small_quotient_norm(factor, mode.as_expr())) for _name, factor in factors)
    return mode, tuple(targets)


def primitive_integer_poly(expression: sp.Expr, variable: sp.Symbol) -> sp.Poly:
    polynomial = sp.Poly(expression, variable, domain=QQ)
    _denominator, integer = polynomial.clear_denoms(convert=True)
    _content, primitive = integer.primitive()
    return primitive


def positive_shift_record(expression: sp.Expr) -> dict[str, object]:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    if sp.Poly(denominator, QPAR, domain=QQ).eval(2) < 0:
        numerator, denominator = -numerator, -denominator
    numerator_poly = primitive_integer_poly(numerator, QPAR)
    denominator_poly = primitive_integer_poly(denominator, QPAR)
    shifted_numerator = sp.Poly(sp.expand(numerator_poly.as_expr().subs(QPAR, 1 + X)), X, domain=ZZ)
    shifted_denominator = sp.Poly(
        sp.expand(denominator_poly.as_expr().subs(QPAR, 1 + X)), X, domain=ZZ
    )
    numerator_coefficients = [int(value) for value in reversed(shifted_numerator.all_coeffs())]
    denominator_coefficients = [int(value) for value in reversed(shifted_denominator.all_coeffs())]
    if not numerator_coefficients or any(value <= 0 for value in numerator_coefficients):
        raise AssertionError("shifted numerator is not coefficientwise positive")
    if not denominator_coefficients or any(value < 0 for value in denominator_coefficients):
        raise AssertionError("shifted denominator is not coefficientwise nonnegative")
    if not any(denominator_coefficients):
        raise AssertionError("shifted denominator is zero")
    return {
        "numerator_degree": numerator_poly.degree(),
        "denominator_degree": denominator_poly.degree(),
        "shifted_numerator_term_count": len(shifted_numerator.terms()),
        "shifted_denominator_term_count": len(shifted_denominator.terms()),
        "shifted_numerator_minimum_coefficient": str(min(numerator_coefficients)),
        "shifted_denominator_zero_coefficients": sum(not value for value in denominator_coefficients),
        "shifted_numerator_sha256": canonical_sha256([str(value) for value in numerator_coefficients]),
        "shifted_denominator_sha256": canonical_sha256([str(value) for value in denominator_coefficients]),
    }


def lucas_gate_record(
    targets: dict[str, sp.Expr], factors: tuple[tuple[str, sp.Expr], ...]
) -> dict[str, object]:
    shifted_factors = {
        name: sp.Poly(sp.expand(factor.subs(U, X + 4)), X, domain=ZZ)
        for name, factor in factors
    }
    f4 = shifted_factors["r4"].as_expr()
    f5 = shifted_factors["r5_over_r1"].as_expr()
    f6 = shifted_factors["r6_over_r2"].as_expr()
    f7 = shifted_factors["r7_over_r1"].as_expr()
    f8 = shifted_factors["r8"].as_expr()
    f9 = shifted_factors["r9_over_r1"].as_expr()
    pointwise = {
        "f4_minus_f6": sp.expand(f4 - f6),
        "f7_minus_f5": sp.expand(f7 - f5),
        "f8_minus_f9": sp.expand(f8 - f9),
        "f4_squared_minus_f8": sp.expand(f4**2 - f8),
        "f4_f6_minus_f5_squared": sp.expand(f4 * f6 - f5**2),
        "f5_f7_minus_f6_squared": sp.expand(f5 * f7 - f6**2),
        "f6_f8_minus_f7_squared": sp.expand(f6 * f8 - f7**2),
    }
    pointwise_coefficients: dict[str, list[int]] = {}
    for name, expression in pointwise.items():
        coefficients = [int(value) for value in reversed(sp.Poly(expression, X, domain=ZZ).all_coeffs())]
        if any(value < 0 for value in coefficients) or not any(coefficients):
            raise AssertionError(f"pointwise Lucas gate {name} is not nonnegative")
        pointwise_coefficients[name] = coefficients
    r4 = targets["r4"]
    r5 = targets["r5_over_r1"]
    r6 = targets["r6_over_r2"]
    r7 = targets["r7_over_r1"]
    r8 = targets["r8"]
    r9 = targets["r9_over_r1"]
    physical = {
        "r4_gt_r6": sp.cancel(r4 - r6),
        "r7_gt_r5": sp.cancel(r7 - r5),
        "r8_gt_r9": sp.cancel(r8 - r9),
        "r4_squared_gt_r8": sp.cancel(r4**2 - r8),
        "r4_r6_gt_r5_squared": sp.cancel(r4 * r6 - r5**2),
        "r5_r7_ge_r6_squared": sp.cancel(r5 * r7 - r6**2),
        "r6_r8_gt_r7_squared": sp.cancel(r6 * r8 - r7**2),
    }
    return {
        "shifted_factor_coefficients_ascending": {
            name: [int(value) for value in reversed(polynomial.all_coeffs())]
            for name, polynomial in shifted_factors.items()
        },
        "pointwise_difference_coefficients_ascending": pointwise_coefficients,
        "physical_target_certificates": {
            name: positive_shift_record(expression) for name, expression in physical.items()
        },
        "all_strict_for_literal_open_2x4_targets": True,
    }


def projection_degree_record(
    targets: dict[str, sp.Expr], equations: list[sp.Poly]
) -> dict[str, object]:
    summaries = {name: rational_summary(value) for name, value in targets.items()}
    first_names = TARGET_NAMES[:3]
    first_denominator_degrees = [int(summaries[name]["denominator_degree"]) for name in first_names]
    common_denominator_degree = sum(first_denominator_degrees)
    dq_degree_bound = max(
        int(summaries[name]["numerator_degree"])
        + common_denominator_degree
        - int(summaries[name]["denominator_degree"])
        for name in first_names
    )
    clearing_resultant_degrees = [2, 2, 2, 3, 4, 4]
    row_summaries = [summaries[name] for name in ROW_TARGET_NAMES]
    cleared_row_bounds = [
        max(
            dq_degree_bound * degree + int(summary["denominator_degree"]),
            common_denominator_degree * degree + int(summary["numerator_degree"]),
        )
        for degree, summary in zip(clearing_resultant_degrees, row_summaries)
    ]
    equation_degrees = [int(equation.total_degree()) for equation in equations]
    resultant_multidegrees = [
        math.prod(equation_degrees[index] for index in range(6) if index != excluded)
        for excluded in range(6)
    ]
    univariate_degree_bound = sum(
        degree * weight for degree, weight in zip(cleared_row_bounds, resultant_multidegrees)
    )
    return {
        "first_target_numerator_degrees": [
            int(summaries[name]["numerator_degree"]) for name in first_names
        ],
        "first_target_denominator_degrees": first_denominator_degrees,
        "common_D_degree": common_denominator_degree,
        "DQ_coefficient_degree_bound": dq_degree_bound,
        "row_target_numerator_denominator_degrees": [
            [int(summary["numerator_degree"]), int(summary["denominator_degree"])]
            for summary in row_summaries
        ],
        "lucas_resultant_degrees_for_clearing": clearing_resultant_degrees,
        "cleared_row_q_degree_bounds": cleared_row_bounds,
        "affine_equation_degrees": equation_degrees,
        "projective_resultant_multidegrees": resultant_multidegrees,
        "univariate_degree_bound": univariate_degree_bound,
        "denominator_scope": (
            "the cleared resultant may contain denominator factors, but all inherited target "
            "denominators are strictly positive for physical q>1"
        ),
    }


def source_hashes() -> dict[str, str]:
    paths = (
        ROOT / SCRIPT,
        VERIFIER,
        PROOF,
        PLAN,
        BASE_PRODUCER,
        BASE_VERIFIER,
        BASE_PROOF,
        BASE_ARTIFACT,
        LOCKFILE,
    )
    return {str(path.relative_to(ROOT.parent.parent)): file_sha256(path) for path in paths}


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    base = json.loads(BASE_ARTIFACT.read_text())
    if not all(check["passed"] for check in base["checks"]):
        raise AssertionError("inherited e248 artifact contains a failed check")
    for relative, expected in base["meta"]["source_sha256"].items():
        if file_sha256(ROOT / relative) != expected:
            raise AssertionError(f"inherited e248 source drift: {relative}")
    trace_rows = [[int(value) for value in row] for row in base["data"]["trace_coefficients_ascending"]]
    if len(trace_rows) != 9:
        raise AssertionError("e248 trace artifact does not contain nine rows")
    trace_digest = canonical_sha256([[str(value) for value in row] for row in trace_rows])
    if trace_digest != base["data"]["trace_coefficients_sha256"]:
        raise AssertionError("e248 trace coefficient digest mismatch")
    targets = normalized_targets(trace_rows)
    target_summaries = {name: rational_summary(value) for name, value in targets.items()}
    if target_summaries != base["data"]["normalized_targets"]:
        raise AssertionError("independent normalized target records differ from e248")
    add_check(checks, "C0_inherited_e248_provenance", True, trace_digest)
    guard_resources(started, "inherited input")

    equations, factors = abstract_mode_equations()
    compiled = compile_equations(equations)
    quotient = generic_quotient(equations)
    leading_monomials: list[Monomial] = quotient["leading_monomials"]
    standard_monomials: list[Monomial] = quotient["standard_monomials"]
    expected_equation_degrees = [2, 2, 2, 3, 4, 3]
    expected_hilbert = [1, 5, 12, 19, 22, 19, 12, 5, 1]
    if [equation.total_degree() for equation in equations] != expected_equation_degrees:
        raise AssertionError("abstract equation degree vector changed")
    if len(quotient["basis"].polys) != 35 or len(standard_monomials) != 96:
        raise AssertionError("generic quotient dimensions changed")
    if quotient["hilbert_vector"] != expected_hilbert:
        raise AssertionError("generic quotient Hilbert vector changed")
    rule_map = reduction_rule_map(leading_monomials, maximum_degree=11)
    coefficient_system = {
        "coefficient_variables": [str(value) for value in COEFFICIENT_VARIABLES],
        "equation_names": [f"F{index}" for index in range(4, 10)],
        "equation_total_degrees": expected_equation_degrees,
        "leading_form_term_counts": [len(value.terms()) for value in quotient["leading_forms"]],
        "groebner_basis_size": len(quotient["basis"].polys),
        "leading_monomials_sha256": canonical_sha256([list(value) for value in leading_monomials]),
        "pure_power_bounds": quotient["pure_power_bounds"],
        "quotient_basis_rank": len(standard_monomials),
        "standard_monomials_sha256": canonical_sha256([list(value) for value in standard_monomials]),
        "maximum_standard_degree": max(map(sum, standard_monomials)),
        "hilbert_vector_degrees_0_through_8": expected_hilbert,
        "reduction_monomials_through_degree_11": sum(1 for _ in monomials_through_degree(5, 11)),
    }
    add_check(checks, "C1_abstract_incidence_and_hilbert", True, "35 relations; rank 96")
    guard_resources(started, "abstract incidence")

    physical_witnesses = []
    witness_specs = (("2", Fraction(2), PRIMES[0]), ("5/3", Fraction(5, 3), PRIMES[1]))
    for label, q_value, prime in witness_specs:
        if not sp.isprime(prime) or prime * prime >= 1 << 63:
            raise AssertionError(f"unsafe modular prime {prime}")
        target_residues, denominator_record = physical_target_residues(trace_rows, q_value, prime)
        specialized = specialize_equations(compiled, target_residues, prime)
        certificate = quotient_certificate(
            specialized,
            prime,
            leading_monomials,
            standard_monomials,
            rule_map,
            started,
            f"physical q={label}",
        )
        certificate.update(
            {
                "q": label,
                "q_residue": denominator_record.pop("q_residue"),
                "target_residues": list(target_residues),
                "target_residues_sha256": canonical_sha256(list(target_residues)),
                "nonzero_denominator_residues": denominator_record,
            }
        )
        if not certificate["determinant_residue"] or certificate["matrix_rank"] != 96:
            raise AssertionError(f"physical q={label} did not give a nonzero norm witness")
        physical_witnesses.append(certificate)
    add_check(
        checks,
        "C2_two_physical_modular_witnesses",
        True,
        "; ".join(
            f"q={entry['q']}, p={entry['prime']}, det={entry['determinant_residue']}"
            for entry in physical_witnesses
        ),
    )
    guard_resources(started, "physical witnesses")

    control_mode, control_targets_integer = synthetic_mode_targets(factors)
    control_prime = PRIMES[0]
    control_targets = tuple(value % control_prime for value in control_targets_integer)
    control_equations = specialize_equations(compiled, control_targets, control_prime)
    control_coefficients = [int(value) % control_prime for value in control_mode.all_coeffs()[1:6]]
    substitution = dict(zip(COEFFICIENT_VARIABLES, control_coefficients))
    point_residuals = [int(equation.eval(substitution)) % control_prime for equation in control_equations]
    if any(point_residuals):
        raise AssertionError("synthetic mode coefficient point does not annihilate all six rows")
    control_certificate = quotient_certificate(
        control_equations,
        control_prime,
        leading_monomials,
        standard_monomials,
        rule_map,
        started,
        "synthetic eight-mode control",
    )
    control_certificate.update(
        {
            "mode_roots": list(range(4, 12)),
            "mode_polynomial_coefficients_descending": [str(value) for value in control_mode.all_coeffs()],
            "target_integers": [str(value) for value in control_targets_integer],
            "target_residues": list(control_targets),
            "coefficient_point_residues": control_coefficients,
            "point_residuals": point_residuals,
        }
    )
    if control_certificate["determinant_residue"] != 0 or control_certificate["matrix_rank"] >= 96:
        raise AssertionError("synthetic physical mode control did not force a singular F9 matrix")
    add_check(
        checks,
        "C3_synthetic_eight_mode_control",
        True,
        f"rank {control_certificate['matrix_rank']}; determinant zero",
    )
    guard_resources(started, "synthetic control")

    degree_record = projection_degree_record(targets, equations)
    if degree_record["univariate_degree_bound"] != 182_400:
        raise AssertionError("projection degree bound changed")
    degree_record["nonzero_specialization"] = {
        "q": physical_witnesses[0]["q"],
        "prime": physical_witnesses[0]["prime"],
        "determinant_residue": physical_witnesses[0]["determinant_residue"],
    }
    add_check(checks, "C4_nonzero_projection_degree_bound", True, "degree <= 182400")

    gate_record = lucas_gate_record(targets, factors)
    add_check(checks, "C5_pointwise_lucas_product_gates", True, "seven shifted-factor gates")
    add_check(checks, "C6_literal_target_gate_positivity", True, "seven q=1+x positive numerators")
    guard_resources(started, "Lucas product gates")

    theorem = {
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
    add_check(checks, "C7_theorem_scope_gate", theorem["emptiness_status"] == "UNRESOLVED", theorem["scope"])

    data = {
        "claim_tag": "[THEOREM][FINITE TRACE-NINE PROJECTION][EXACT MODULAR WITNESS]",
        "graph": {
            "name": "open 2x4 grid",
            "sites": 8,
            "in_layer_edges": 10,
            "transfer_dimension": 256,
            "physical_parameter": "q=(1+t^2)/(2t)>1",
        },
        "inherited_e248": {
            "artifact_sha256": file_sha256(BASE_ARTIFACT),
            "trace_coefficients_sha256": trace_digest,
            "trace_count": len(trace_rows),
            "base_check_names": [check["name"] for check in base["checks"]],
            "normalized_target_sha256": canonical_sha256(target_summaries),
        },
        "coefficient_system": coefficient_system,
        "modular_certificates": {
            "physical_witnesses": physical_witnesses,
            "synthetic_control": control_certificate,
        },
        "projection_degree": degree_record,
        "lucas_product_gates": gate_record,
        "resource_boundary": {
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
        },
        "theorem": theorem,
    }

    elapsed = time.process_time() - started
    peak = peak_rss_bytes()
    add_check(checks, "C8_cpu_budget", elapsed < CPU_LIMIT_SECONDS, f"{elapsed:.6f} < {CPU_LIMIT_SECONDS}")
    add_check(checks, "C9_rss_budget", peak < RSS_LIMIT_BYTES, f"{peak} < {RSS_LIMIT_BYTES}")
    if not all(check["passed"] for check in checks):
        failed = [check["name"] for check in checks if not check["passed"]]
        raise AssertionError(f"producer checks failed: {failed}")
    payload = {
        "checks": checks,
        "data": data,
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": SCRIPT,
            "verifier": str(VERIFIER.relative_to(ROOT)),
            "proof": str(PROOF.relative_to(ROOT)),
            "plan": str(PLAN.relative_to(ROOT)),
            "interpreter": os.path.realpath(os.sys.executable),
            "python_version": platform.python_version(),
            "sympy_version": sp.__version__,
            "numpy_version": np.__version__,
            "arithmetic": "exact ZZ[q], QQ targets, GF(p) Groebner bases, and GF(p) determinants",
            "single_process": True,
            "process_cpu_seconds": elapsed,
            "process_cpu_limit_seconds": CPU_LIMIT_SECONDS,
            "peak_rss_bytes": peak,
            "peak_rss_measurement": peak_rss_measurement(),
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "source_sha256": source_hashes(),
            "data_sha256": canonical_sha256(data),
            "benchmark_Kc_used": False,
        },
    }
    guard_resources(started, "final write")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> int:
    payload = run()
    for check in payload["checks"]:
        print(f"{'PASS' if check['passed'] else 'FAIL'} {check['name']}: {check['detail']}")
    meta = payload["meta"]
    print(
        "PASS e251 trace-nine projection; "
        f"CPU={meta['process_cpu_seconds']:.6f}s RSS={meta['peak_rss_bytes']} bytes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
