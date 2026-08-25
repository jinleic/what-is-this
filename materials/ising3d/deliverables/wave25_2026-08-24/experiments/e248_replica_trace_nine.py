#!/usr/bin/env python3
"""Exact replica-column traces for the open 2x4 Ising layer.

The high-risk second stage builds the eight-mode Lucas/resultant incidence through
trace nine.  It deliberately refuses the underdetermined trace-eight shortcut.
"""
from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import platform
import resource
import sys
import time
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ, ZZ

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e248_replica_trace_nine.py"
OUTPUT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
VERIFIER = ROOT / "tests" / "test_replica_trace_nine.py"
PROOF = ROOT / "proofs" / "replica_trace_nine.md"
E238_PRODUCER = ROOT / "experiments" / "e238_trace_exceptional_set.py"
E238_ARTIFACT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
CPU_LIMIT_SECONDS = 900.0
RSS_LIMIT_BYTES = 2 * 1024**3
CONTROL_DIGEST_2X3 = "0b3a1f464f8776ba2f0ede1214538e998874e021ce74531fe6f519d12be714ff"
Q = sp.symbols("q")
U = sp.symbols("u")
A7, A6, A5, A4, A3 = sp.symbols("a7 a6 a5 a4 a3")
COEFFICIENT_VARIABLES = (A7, A6, A5, A4, A3)


class ResourceWall(RuntimeError):
    """Raised before a declared process resource ceiling is crossed."""


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        mach = ctypes.CDLL("/usr/lib/libSystem.B.dylib", use_errno=True)

        class MachTaskBasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_time_seconds", ctypes.c_int32),
                ("user_time_microseconds", ctypes.c_int32),
                ("system_time_seconds", ctypes.c_int32),
                ("system_time_microseconds", ctypes.c_int32),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        mach.mach_task_self.restype = ctypes.c_uint32
        mach.task_info.argtypes = [
            ctypes.c_uint32,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        mach.task_info.restype = ctypes.c_int
        info = MachTaskBasicInfo()
        count = ctypes.c_uint32(ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_uint32))
        result = mach.task_info(
            mach.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count)
        )
        if result != 0:
            error = ctypes.get_errno()
            raise RuntimeError(f"task_info failed with kern_return={result}, errno={error}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def peak_rss_measurement() -> str:
    if platform.system() == "Darwin":
        return "mach_task_basic_info.resident_size_max (task_info flavor 20)"
    return "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"


def guard_resources(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    if elapsed >= CPU_LIMIT_SECONDS:
        raise ResourceWall(f"{stage}: process CPU {elapsed:.6f}s reached {CPU_LIMIT_SECONDS}s")
    peak = peak_rss_bytes()
    if peak >= RSS_LIMIT_BYTES:
        raise ResourceWall(f"{stage}: peak RSS {peak} reached {RSS_LIMIT_BYTES} bytes")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def integer_polynomial_sha256(rows: Iterable[Iterable[int]]) -> str:
    return canonical_sha256([[str(int(value)) for value in row] for row in rows])


def column_local_degrees(replica_count: int) -> np.ndarray:
    """Satisfied rung plus directed replica-cycle edges for one 2-spin column."""
    state_count = 1 << (2 * replica_count)
    states = np.arange(state_count, dtype=np.uint32)
    degrees = np.zeros(state_count, dtype=np.uint16)
    for replica in range(replica_count):
        top = (states >> (2 * replica)) & 1
        bottom = (states >> (2 * replica + 1)) & 1
        successor = (replica + 1) % replica_count
        degrees += top == bottom
        degrees += top == ((states >> (2 * successor)) & 1)
        degrees += bottom == ((states >> (2 * successor + 1)) & 1)
    return degrees


def replica_trace_coefficients(length: int, replica_count: int) -> list[int]:
    """Return H_k(q)=tr(B(q)^k) for the open 2xL layer.

    The vector for a fixed final column counts at most 2^(2*k*(L-1)) partial
    configurations.  At L=4,k=9 this is 2^54, so packed uint64 additions are
    exact.  Only the final state sum can exceed uint64; it is accumulated in
    128-row blocks and then in Python integers.
    """
    if length < 1 or replica_count < 1:
        raise ValueError("length and replica_count must be positive")
    state_count = 1 << (2 * replica_count)
    total_degree = replica_count * (5 * length - 2)
    partial_bound_bits = 2 * replica_count * (length - 1)
    if partial_bound_bits >= 64:
        raise ResourceWall(
            f"uint64 partial-count proof fails: 2^{partial_bound_bits} configurations"
        )
    slot_count = state_count * (total_degree + 1)
    packed_bytes = 2 * slot_count * np.dtype(np.uint64).itemsize
    if packed_bytes >= 1024**3:
        raise ResourceWall(
            f"two packed trace arrays require {packed_bytes} bytes, at least 1 GiB"
        )

    first = np.zeros((state_count, total_degree + 1), dtype=np.uint64)
    second = np.zeros_like(first)
    degrees = column_local_degrees(replica_count)
    first[np.arange(state_count), degrees] = 1
    current, spare = first, second
    width = 3 * replica_count + 1

    for _column in range(1, length):
        for bit in range(2 * replica_count):
            new_width = width + 1
            spare[:, :new_width] = 0
            block = 1 << bit
            source = current.reshape((-1, 2, block, total_degree + 1))
            target = spare.reshape((-1, 2, block, total_degree + 1))
            target[:, 0, :, :width] = source[:, 1, :, :width]
            target[:, 1, :, :width] = source[:, 0, :, :width]
            np.add(
                target[:, 0, :, 1:new_width],
                source[:, 0, :, :width],
                out=target[:, 0, :, 1:new_width],
            )
            np.add(
                target[:, 1, :, 1:new_width],
                source[:, 1, :, :width],
                out=target[:, 1, :, 1:new_width],
            )
            current, spare = spare, current
            width = new_width

        new_width = width + 3 * replica_count
        spare[:, :new_width] = 0
        for degree in np.unique(degrees):
            indices = np.flatnonzero(degrees == degree)
            offset = int(degree)
            spare[indices, offset : offset + width] = current[indices, :width]
        current, spare = spare, current
        width = new_width

    if width != total_degree + 1:
        raise AssertionError(f"degree width {width} != {total_degree + 1}")
    block_rows = min(128, state_count)
    if state_count % block_rows:
        raise AssertionError("state count must be divisible by final summation block")
    final_block_bound_bits = partial_bound_bits + (block_rows.bit_length() - 1)
    if final_block_bound_bits >= 64:
        raise ResourceWall(
            f"uint64 final block-sum proof fails at 2^{final_block_bound_bits}"
        )
    partials = current.reshape(
        (state_count // block_rows, block_rows, width)
    ).sum(axis=1, dtype=np.uint64)
    coefficients = [sum(int(value) for value in partials[:, degree]) for degree in range(width)]
    if sum(coefficients) != 1 << (2 * replica_count * length):
        raise AssertionError("H_k(1) does not count all replica configurations")
    return coefficients


def dense_trace_coefficients(length: int, maximum: int) -> list[list[int]]:
    """Independent small dense ZZ[q] construction used only at L<=2."""
    if length > 2:
        raise ValueError("dense control is intentionally limited to L<=2")
    site_count = 2 * length
    bonds = [(row * length + column, row * length + column + 1)
             for row in range(2) for column in range(length - 1)]
    bonds += [(column, length + column) for column in range(length)]
    state_count = 1 << site_count

    def disagreements(state: int) -> int:
        return sum(((state >> left) ^ (state >> right)) & 1 for left, right in bonds)

    matrix: list[list[list[int]]] = []
    for row in range(state_count):
        satisfied = len(bonds) - disagreements(row)
        matrix.append([])
        for column in range(state_count):
            equal = site_count - (row ^ column).bit_count()
            polynomial = [0] * (len(bonds) + site_count + 1)
            polynomial[satisfied + equal] = 1
            matrix[-1].append(polynomial)

    def add_product(left: list[int], right: list[int], out: list[int]) -> None:
        for i, left_value in enumerate(left):
            if not left_value:
                continue
            for j, right_value in enumerate(right):
                if right_value:
                    out[i + j] += left_value * right_value

    identity = [[[1] if row == column else [0] for column in range(state_count)]
                for row in range(state_count)]
    current = identity
    traces: list[list[int]] = []
    for power in range(1, maximum + 1):
        degree = power * (len(bonds) + site_count)
        product = [[[0] * (degree + 1) for _ in range(state_count)]
                   for _ in range(state_count)]
        for row in range(state_count):
            for middle in range(state_count):
                left = current[row][middle]
                if not any(left):
                    continue
                for column in range(state_count):
                    add_product(left, matrix[middle][column], product[row][column])
        current = product
        trace = [0] * (degree + 1)
        for index in range(state_count):
            for degree_index, value in enumerate(current[index][index]):
                trace[degree_index] += value
        traces.append(trace)
    return traces


def expression_from_coefficients(coefficients: list[int]) -> sp.Expr:
    return sp.Poly.from_list(list(reversed(coefficients)), gens=Q, domain=ZZ).as_expr()


def normalized_targets(traces: list[sp.Expr], site_count: int, epsilon: int, shift: int) -> dict[str, sp.Expr]:
    delta = Q**2 - 1

    def scale(center_power: int) -> sp.Expr:
        return Q ** (center_power * epsilon - 2 * center_power * shift) / delta ** (
            center_power * site_count
        )

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


def rational_record(value: sp.Expr) -> dict[str, object]:
    numerator, denominator = sp.cancel(value).as_numer_denom()
    numerator_poly = sp.Poly(numerator, Q, domain=QQ)
    denominator_poly = sp.Poly(denominator, Q, domain=QQ)
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


def mode_polynomial_data(
    targets: dict[str, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, tuple[sp.Expr, ...]]:
    a0 = targets["r1_squared"]
    c2 = targets["r2"] - a0 - 256 - 128 * A7 - 64 * A6 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = targets["r3_over_r1"] - a0 - 6561 - 2187 * A7 - 729 * A6 - 243 * A5 - 81 * A4 - 27 * A3
    a2 = sp.cancel(c3 / 3 - c2 / 2)
    a1 = sp.cancel(3 * c2 / 2 - 2 * c3 / 3)
    coefficients_descending = (
        sp.Integer(1),
        A7,
        A6,
        A5,
        A4,
        A3,
        a2,
        a1,
        a0,
    )
    return a2, a1, coefficients_descending


def lucas_factor_rows() -> tuple[tuple[str, sp.Expr], ...]:
    tau = sp.symbols("tau")
    lucas = [sp.Integer(2), tau]
    for _index in range(2, 10):
        lucas.append(sp.expand(tau * lucas[-1] - lucas[-2]))

    def quotient_in_u(index: int, divisor: sp.Expr) -> sp.Expr:
        quotient, remainder = sp.div(lucas[index], divisor, tau, domain=ZZ)
        if remainder:
            raise AssertionError(f"Lucas L_{index} has nonzero requested remainder")
        polynomial = sp.Poly(sp.expand(quotient), tau, domain=ZZ)
        expression = sp.Integer(0)
        for (power,), coefficient in polynomial.terms():
            if power % 2:
                raise AssertionError(f"Lucas L_{index} quotient is not even in tau")
            expression += coefficient * U ** (power // 2)
        return sp.expand(expression)

    return (
        ("r4", quotient_in_u(4, sp.Integer(1))),
        ("r5_over_r1", quotient_in_u(5, tau)),
        ("r6_over_r2", quotient_in_u(6, lucas[2])),
        ("r7_over_r1", quotient_in_u(7, tau)),
        ("r8", quotient_in_u(8, sp.Integer(1))),
        ("r9_over_r1", quotient_in_u(9, tau)),
    )


def mode_equations(
    targets: dict[str, sp.Expr],
) -> tuple[list[sp.Poly], sp.Expr, sp.Expr, tuple[sp.Expr, ...]]:
    target_symbols = sp.symbols("s0 s2 s3 s4 s5 s6 s7 s8 s9")
    s0, s2, s3, s4, s5, s6, s7, s8, s9 = target_symbols
    c2 = s2 - s0 - 256 - 128 * A7 - 64 * A6 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = s3 - s0 - 6561 - 2187 * A7 - 729 * A6 - 243 * A5 - 81 * A4 - 27 * A3
    a2_abstract = c3 / 3 - c2 / 2
    a1_abstract = 3 * c2 / 2 - 2 * c3 / 3
    mode_abstract = (
        U**8 + A7 * U**7 + A6 * U**6 + A5 * U**5 + A4 * U**4 + A3 * U**3
        + a2_abstract * U**2 + a1_abstract * U + s0
    )
    target_symbols_by_name = {
        "r4": s4,
        "r5_over_r1": s5,
        "r6_over_r2": s6,
        "r7_over_r1": s7,
        "r8": s8,
        "r9_over_r1": s9,
    }
    factors = [
        (factor, target_symbols_by_name[name])
        for name, factor in lucas_factor_rows()
    ]
    abstract_field = QQ.frac_field(*target_symbols)
    abstract_equations = [
        sp.Poly(
            sp.resultant(factor, mode_abstract, U) - target,
            *COEFFICIENT_VARIABLES,
            domain=abstract_field,
        )
        for factor, target in factors
    ]
    substitutions = {
        s0: targets["r1_squared"],
        s2: targets["r2"],
        s3: targets["r3_over_r1"],
        s4: targets["r4"],
        s5: targets["r5_over_r1"],
        s6: targets["r6_over_r2"],
        s7: targets["r7_over_r1"],
        s8: targets["r8"],
        s9: targets["r9_over_r1"],
    }
    field = QQ.frac_field(Q)
    equations = [
        sp.Poly.from_dict(
            {
                monomial: sp.cancel(coefficient.as_expr().subs(substitutions))
                for monomial, coefficient in abstract.terms()
            },
            COEFFICIENT_VARIABLES,
            domain=field,
        )
        for abstract in abstract_equations
    ]
    a2, a1, mode_coefficients = mode_polynomial_data(targets)
    return equations, a2, a1, mode_coefficients


def leading_form(poly: sp.Poly) -> sp.Poly:
    degree = poly.total_degree()
    expression = sum(
        coefficient * sp.prod(variable**power for variable, power in zip(COEFFICIENT_VARIABLES, monomial))
        for monomial, coefficient in poly.terms()
        if sum(monomial) == degree
    )
    simplified = sp.Poly(expression, *COEFFICIENT_VARIABLES, domain=QQ.frac_field(Q))
    for _monomial, coefficient in simplified.terms():
        if Q in coefficient.as_expr().free_symbols:
            raise AssertionError("leading coefficient unexpectedly depends on q")
    return sp.Poly(simplified.as_expr(), *COEFFICIENT_VARIABLES, domain=QQ)


def leading_ideal_record(equations: list[sp.Poly]) -> dict[str, object]:
    first_five = equations[:5]
    forms = [leading_form(equation) for equation in first_five]
    basis = sp.groebner(
        [form.as_expr() for form in forms], *COEFFICIENT_VARIABLES, order="grevlex", domain=QQ
    )
    return {
        "equation_total_degrees": [int(equation.total_degree()) for equation in equations],
        "leading_form_term_counts": [len(form.terms()) for form in forms],
        "groebner_basis_size": len(basis.polys),
        "is_zero_dimensional": bool(basis.is_zero_dimensional),
        "groebner_sha256": canonical_sha256([str(poly.as_expr()) for poly in basis.polys]),
    }


def abstract_leading_ideal_record() -> dict[str, object]:
    target_symbols = sp.symbols("s0 s2 s3 s4 s5 s6 s7 s8 s9")
    s0, s2, s3, s4, s5, s6, s7, s8, s9 = target_symbols
    c2 = s2 - s0 - 256 - 128 * A7 - 64 * A6 - 32 * A5 - 16 * A4 - 8 * A3
    c3 = s3 - s0 - 6561 - 2187 * A7 - 729 * A6 - 243 * A5 - 81 * A4 - 27 * A3
    a2 = c3 / 3 - c2 / 2
    a1 = 3 * c2 / 2 - 2 * c3 / 3
    mode = (
        U**8 + A7 * U**7 + A6 * U**6 + A5 * U**5 + A4 * U**4 + A3 * U**3
        + a2 * U**2 + a1 * U + s0
    )
    target_symbols_by_name = {
        "r4": s4,
        "r5_over_r1": s5,
        "r6_over_r2": s6,
        "r7_over_r1": s7,
        "r8": s8,
        "r9_over_r1": s9,
    }
    factors = [
        (factor, target_symbols_by_name[name])
        for name, factor in lucas_factor_rows()
    ]
    field = QQ.frac_field(*target_symbols)
    equations = [
        sp.Poly(
            sp.resultant(factor, mode, U) - target,
            *COEFFICIENT_VARIABLES,
            domain=field,
        )
        for factor, target in factors
    ]
    first_five = equations[:5]
    forms: list[sp.Poly] = []
    for equation in first_five:
        degree = equation.total_degree()
        expression = sum(
            coefficient.as_expr()
            * sp.prod(
                variable**power
                for variable, power in zip(COEFFICIENT_VARIABLES, monomial)
            )
            for monomial, coefficient in equation.terms()
            if sum(monomial) == degree
        )
        if set(expression.free_symbols) & set(target_symbols):
            raise AssertionError("abstract leading form depends on trace targets")
        forms.append(sp.Poly(expression, *COEFFICIENT_VARIABLES, domain=QQ))
    basis = sp.groebner(
        [form.as_expr() for form in forms],
        *COEFFICIENT_VARIABLES,
        order="grevlex",
        domain=QQ,
    )
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents for polynomial in basis.polys
    ]
    pure_power_bounds: list[int] = []
    for variable_index in range(len(COEFFICIENT_VARIABLES)):
        candidates = [
            monomial[variable_index]
            for monomial in leading_monomials
            if monomial[variable_index] > 0
            and sum(monomial) == monomial[variable_index]
        ]
        if not candidates:
            raise AssertionError(f"no pure leading power for variable {variable_index}")
        pure_power_bounds.append(min(candidates))
    standard_monomials = [
        exponents
        for exponents in product(*(range(bound) for bound in pure_power_bounds))
        if not any(
            all(leading[index] <= exponents[index] for index in range(5))
            for leading in leading_monomials
        )
    ]
    return {
        "equation_total_degrees": [int(equation.total_degree()) for equation in equations],
        "leading_form_term_counts": [len(form.terms()) for form in forms],
        "groebner_basis_size": len(basis.polys),
        "is_zero_dimensional": bool(basis.is_zero_dimensional),
        "pure_power_bounds": pure_power_bounds,
        "quotient_basis_rank": len(standard_monomials),
        "maximum_standard_degree": max(map(sum, standard_monomials)),
        "standard_monomials_sha256": canonical_sha256(standard_monomials),
        "leading_monomials_sha256": canonical_sha256(leading_monomials),
        "groebner_sha256": canonical_sha256(
            [
                [[list(monomial), str(coefficient)] for monomial, coefficient in polynomial.terms()]
                for polynomial in basis.polys
            ]
        ),
    }


def actual_quotient_probe(
    equations: list[sp.Poly], abstract_record: dict[str, object]
) -> dict[str, object]:
    field = QQ.frac_field(Q)
    basis = sp.groebner(
        [equation.as_expr() for equation in equations[:5]],
        *COEFFICIENT_VARIABLES,
        order="grevlex",
        domain=field,
        method="f5b",
    )
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents for polynomial in basis.polys
    ]
    bounds = [int(value) for value in abstract_record["pure_power_bounds"]]
    standard_monomials = [
        exponents
        for exponents in product(*(range(bound) for bound in bounds))
        if not any(
            all(leading[index] <= exponents[index] for index in range(5))
            for leading in leading_monomials
        )
    ]
    numerator_degrees: list[int] = []
    denominator_degrees: list[int] = []
    coefficient_count = 0
    for polynomial in basis.polys:
        for _monomial, coefficient in polynomial.terms():
            coefficient_count += 1
            numerator_degrees.append(coefficient.numer.degree())
            denominator_degrees.append(coefficient.denom.degree())
    leading_digest = canonical_sha256(leading_monomials)
    return {
        "basis_size": len(basis.polys),
        "basis_term_count": coefficient_count,
        "is_zero_dimensional": bool(basis.is_zero_dimensional),
        "leading_monomials_sha256": leading_digest,
        "leading_monomials_match_abstract": (
            leading_digest == abstract_record["leading_monomials_sha256"]
        ),
        "quotient_basis_rank": len(standard_monomials),
        "standard_monomials_sha256": canonical_sha256(standard_monomials),
        "maximum_numerator_degree": max(numerator_degrees),
        "maximum_denominator_degree": max(denominator_degrees),
    }


def add_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def run(probe_only: bool = False, quotient_probe: bool = False) -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    dense_control = dense_trace_coefficients(2, 4)
    recurrence_control = [replica_trace_coefficients(2, power) for power in range(1, 5)]
    add_check(checks, "C1_dense_recurrence_control", dense_control == recurrence_control, "open 2x2 traces k=1..4")

    traces_2x3 = [replica_trace_coefficients(3, power) for power in range(1, 8)]
    digest_2x3 = integer_polynomial_sha256(traces_2x3)
    add_check(checks, "C2_landed_2x3_digest", digest_2x3 == CONTROL_DIGEST_2X3, digest_2x3)
    guard_resources(started, "controls")

    trace_timings: list[float] = []
    traces_2x4: list[list[int]] = []
    for power in range(1, 10):
        stage_started = time.process_time()
        traces_2x4.append(replica_trace_coefficients(4, power))
        trace_timings.append(time.process_time() - stage_started)
        guard_resources(started, f"2x4 trace {power}")
    expected_degrees = [18 * power for power in range(1, 10)]
    observed_degrees = [
        max(index for index, coefficient in enumerate(row) if coefficient)
        for row in traces_2x4
    ]
    add_check(
        checks,
        "C3_trace_degrees_2x4",
        observed_degrees == expected_degrees and all(row[-1] == 2 for row in traces_2x4),
        f"observed {observed_degrees}; top coefficients {[row[-1] for row in traces_2x4]}",
    )
    add_check(
        checks,
        "C4_trace_configuration_counts",
        all(sum(row) == 1 << (8 * power) for power, row in enumerate(traces_2x4, 1)),
        "H_k(1)=2^(8k), k=1..9",
    )

    trace_expressions = [expression_from_coefficients(row) for row in traces_2x4]
    targets = normalized_targets(trace_expressions, site_count=8, epsilon=0, shift=5)
    target_records = {name: rational_record(value) for name, value in targets.items()}
    a2, a1, mode_coefficients = mode_polynomial_data(targets)
    lucas_factor_records = {
        name: [int(coefficient) for coefficient in reversed(sp.Poly(factor, U, domain=ZZ).all_coeffs())]
        for name, factor in lucas_factor_rows()
    }
    expected_lucas_factor_records = {
        "r4": [2, -4, 1],
        "r5_over_r1": [5, -5, 1],
        "r6_over_r2": [1, -4, 1],
        "r7_over_r1": [-7, 14, -7, 1],
        "r8": [2, -16, 20, -8, 1],
        "r9_over_r1": [9, -30, 27, -9, 1],
    }
    add_check(
        checks,
        "C9_lucas_recurrence_factors",
        lucas_factor_records == expected_lucas_factor_records,
        "L0=2,L1=tau recurrence through L9 gives all six exact u-factors",
    )
    leading_record = abstract_leading_ideal_record()
    add_check(
        checks,
        "C5_trace_nine_dimension_count",
        len(COEFFICIENT_VARIABLES) == 5
        and len(leading_record["equation_total_degrees"]) == 6,
        "five coefficient variables; F4..F8 define the coefficient ideal and F9 is the first additional row that can cut q",
    )
    add_check(
        checks,
        "C6_expected_equation_degrees",
        leading_record["equation_total_degrees"] == [2, 2, 2, 3, 4, 3],
        str(leading_record["equation_total_degrees"]),
    )
    add_check(
        checks,
        "C7_leading_complete_intersection_rank",
        bool(leading_record["is_zero_dimensional"])
        and leading_record["quotient_basis_rank"] == 96,
        (
            f"basis size {leading_record['groebner_basis_size']}; "
            f"quotient rank {leading_record['quotient_basis_rank']}"
        ),
    )
    coefficient_quotient_theorem = {
        "field": "Q(q), with every recorded target denominator treated as a nonzero field element",
        "statement": "F4 through F8 form a zero-dimensional complete intersection of scheme length 96 in the five mode-coefficient variables",
        "proof": (
            "their q-independent leading forms have no nonzero common affine zero, hence no "
            "projective zero; homogenization adds no point at infinity, and projective Bezout "
            "gives length 2*2*2*3*4=96"
        ),
        "basis_boundary": (
            "length 96 is proved; an explicit lifted Q(q) Groebner basis and the F9 "
            "multiplication matrix are not materialized"
        ),
    }
    add_check(
        checks,
        "C8_filtered_complete_intersection_length",
        bool(leading_record["is_zero_dimensional"])
        and leading_record["equation_total_degrees"][:5] == [2, 2, 2, 3, 4]
        and leading_record["quotient_basis_rank"] == 96,
        "no projective leading-form zero; filtered/projective Bezout length 96 over Q(q)",
    )
    guard_resources(started, "abstract leading ideal")
    quotient_record: dict[str, object] | None = None
    if quotient_probe:
        equations, _a2, _a1, _mode = mode_equations(targets)
        quotient_record = actual_quotient_probe(equations, leading_record)
        add_check(
            checks,
            "C10_actual_quotient_matches_abstract",
            bool(quotient_record["leading_monomials_match_abstract"])
            and quotient_record["quotient_basis_rank"] == 96,
            (
                f"actual basis {quotient_record['basis_size']}; "
                f"rank {quotient_record['quotient_basis_rank']}"
            ),
        )
        guard_resources(started, "actual quotient")

    if probe_only:
        return {
            "checks": checks,
            "probe": {
                "trace_timings_process_seconds": trace_timings,
                "trace_term_counts": [sum(value != 0 for value in row) for row in traces_2x4],
                "target_records": target_records,
                "leading_ideal": leading_record,
                "actual_quotient": quotient_record,
                "process_cpu_seconds": time.process_time() - started,
                "peak_rss_bytes": peak_rss_bytes(),
            },
        }

    all_passed = all(bool(check["passed"]) for check in checks)
    if not all_passed:
        failed = [str(check["name"]) for check in checks if not check["passed"]]
        raise AssertionError(f"refusing to write failed artifact: {failed}")
    payload = {
        "meta": {
            "experiment": "e248",
            "script": SCRIPT,
            "source_sha256": {
                SCRIPT: file_sha256(ROOT / SCRIPT),
                "tests/test_replica_trace_nine.py": file_sha256(VERIFIER),
                "proofs/replica_trace_nine.md": file_sha256(PROOF),
                "experiments/e238_trace_exceptional_set.py": file_sha256(E238_PRODUCER),
                "results/spectral/trace_exceptional_set.json": file_sha256(E238_ARTIFACT),
            },
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": peak_rss_bytes(),
            "peak_rss_measurement": peak_rss_measurement(),
            "cpu_limit_seconds": CPU_LIMIT_SECONDS,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
        },
        "checks": checks,
        "data": {
            "scope": "exact open 2x4 trace-nine incidence preflight; no exceptional-set emptiness claim",
            "replica_column_identity": {
                "statement": "tr(B(q)^k) is the equality-enumerator of the 2 x L x directed-C_k replica graph and is evaluated by a 4^k-state column recurrence",
                "maximum_state_count": 4**9,
                "maximum_degree": 162,
                "uint64_partial_count_bound_bits": 54,
                "final_block_rows": 128,
                "uint64_final_block_sum_bound_bits": 61,
            },
            "trace_coefficients_ascending": [[str(value) for value in row] for row in traces_2x4],
            "trace_coefficients_sha256": integer_polynomial_sha256(traces_2x4),
            "trace_term_counts": [sum(value != 0 for value in row) for row in traces_2x4],
            "trace_timings_process_seconds": trace_timings,
            "normalized_targets": target_records,
            "mode_incidence": {
                "coefficient_variables": [str(variable) for variable in COEFFICIENT_VARIABLES],
                "first_rows": ["Q(0)=r1^2", "Q(2)=r2", "Q(3)=r3/r1"],
                "remaining_rows": ["F4", "F5", "F6", "F7", "F8", "F9"],
                "a2_sha256": canonical_sha256(str(a2)),
                "a1_sha256": canonical_sha256(str(a1)),
                "lucas_factor_coefficients_ascending": lucas_factor_records,
                "mode_polynomial_coefficients_sha256": canonical_sha256(
                    [str(coefficient) for coefficient in mode_coefficients]
                ),
                "leading_ideal": leading_record,
                "coefficient_quotient_theorem": coefficient_quotient_theorem,
                "status": "F4..F8 have generic Q(q) scheme length 96 by filtered/projective Bezout; an explicit lifted basis and the F9 multiplication norm are not materialized",
            },
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--probe",
        action="store_true",
        help="run through the abstract leading-ideal gate without writing",
    )
    parser.add_argument(
        "--quotient-probe",
        action="store_true",
        help="also build the actual Q(q) coefficient quotient without writing",
    )
    arguments = parser.parse_args()
    payload = run(
        probe_only=arguments.probe or arguments.quotient_probe,
        quotient_probe=arguments.quotient_probe,
    )
    print(
        json.dumps(
            payload["probe"] if arguments.probe or arguments.quotient_probe else payload["meta"],
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
