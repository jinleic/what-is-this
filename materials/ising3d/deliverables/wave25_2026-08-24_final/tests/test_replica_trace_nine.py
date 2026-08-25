#!/usr/bin/env python3
"""Clean-room verifier for the e248 replica-column trace-nine artifact.

This file does not import the producer.  It uses a different column-state bit
ordering and independently rebuilds every stored integer coefficient.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
import platform
import resource
import time
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ, ZZ

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
PRODUCER = ROOT / "experiments" / "e248_replica_trace_nine.py"
SELF = ROOT / "tests" / "test_replica_trace_nine.py"
PROOF = ROOT / "proofs" / "replica_trace_nine.md"
E238_PRODUCER = ROOT / "experiments" / "e238_trace_exceptional_set.py"
E238_ARTIFACT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
CPU_LIMIT_SECONDS = 180.0
RSS_LIMIT_BYTES = 2 * 1024**3
CONTROL_DIGEST_2X3 = "0b3a1f464f8776ba2f0ede1214538e998874e021ce74531fe6f519d12be714ff"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        system = ctypes.CDLL("/usr/lib/libSystem.B.dylib", use_errno=True)

        class BasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_seconds", ctypes.c_int32),
                ("user_microseconds", ctypes.c_int32),
                ("system_seconds", ctypes.c_int32),
                ("system_microseconds", ctypes.c_int32),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        system.mach_task_self.restype = ctypes.c_uint32
        system.task_info.argtypes = [
            ctypes.c_uint32,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_uint32),
        ]
        system.task_info.restype = ctypes.c_int
        info = BasicInfo()
        count = ctypes.c_uint32(ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_uint32))
        result = system.task_info(
            system.mach_task_self(), 20, ctypes.byref(info), ctypes.byref(count)
        )
        if result:
            raise RuntimeError(f"task_info failed with kern_return={result}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def local_weight_degrees(replica_count: int) -> np.ndarray:
    """Top replicas occupy low bits; bottom replicas occupy high bits."""
    size = 1 << (2 * replica_count)
    configurations = np.arange(size, dtype=np.uint32)
    weights = np.zeros(size, dtype=np.uint16)
    for replica in range(replica_count):
        successor = (replica + 1) % replica_count
        top = (configurations >> replica) & 1
        bottom = (configurations >> (replica_count + replica)) & 1
        top_next = (configurations >> successor) & 1
        bottom_next = (configurations >> (replica_count + successor)) & 1
        weights += top == bottom
        weights += top == top_next
        weights += bottom == bottom_next
    return weights


def independent_trace(length: int, replica_count: int) -> list[int]:
    size = 1 << (2 * replica_count)
    degree = replica_count * (5 * length - 2)
    weights = local_weight_degrees(replica_count)
    left = np.zeros((size, degree + 1), dtype=np.uint64)
    right = np.zeros_like(left)
    left[np.arange(size), weights] = 1
    active = 3 * replica_count + 1

    for _ in range(length - 1):
        for coordinate in range(2 * replica_count):
            stride = 1 << coordinate
            groups = size // (2 * stride)
            right[:, : active + 1] = 0
            source = left.reshape((groups, 2, stride, degree + 1))
            target = right.reshape((groups, 2, stride, degree + 1))
            target[:, 0, :, :active] = source[:, 1, :, :active]
            target[:, 1, :, :active] = source[:, 0, :, :active]
            target[:, 0, :, 1 : active + 1] += source[:, 0, :, :active]
            target[:, 1, :, 1 : active + 1] += source[:, 1, :, :active]
            left, right = right, left
            active += 1
        right[:, : active + 3 * replica_count] = 0
        for weight in sorted({int(value) for value in weights}):
            rows = np.flatnonzero(weights == weight)
            right[rows, weight : weight + active] = left[rows, :active]
        left, right = right, left
        active += 3 * replica_count

    if active != degree + 1:
        raise AssertionError("independent recurrence degree mismatch")
    block = min(64, size)
    partials = left.reshape((size // block, block, active)).sum(axis=1, dtype=np.uint64)
    result = [sum(map(int, partials[:, index])) for index in range(active)]
    if sum(result) != 1 << (2 * length * replica_count):
        raise AssertionError("independent recurrence configuration count mismatch")
    return result


def coefficient_digest(rows: list[list[int]]) -> str:
    return canonical_sha256([[str(value) for value in row] for row in rows])


def polynomial(coefficients: list[int], variable: sp.Symbol) -> sp.Expr:
    return sp.Poly.from_list(list(reversed(coefficients)), gens=variable, domain=ZZ).as_expr()


def normalized_targets(rows: list[list[int]]) -> dict[str, sp.Expr]:
    q = sp.symbols("q")
    traces = [polynomial(row, q) for row in rows]
    delta = q**2 - 1

    def scale(power: int) -> sp.Expr:
        return q ** (-10 * power) / delta ** (8 * power)

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


def rational_digest(value: sp.Expr) -> tuple[int, int, int, int, str]:
    q = sp.symbols("q")
    numerator, denominator = sp.cancel(value).as_numer_denom()
    numerator_poly = sp.Poly(numerator, q, domain=QQ)
    denominator_poly = sp.Poly(denominator, q, domain=QQ)
    numerator_coefficients = [str(entry) for entry in reversed(numerator_poly.all_coeffs())]
    denominator_coefficients = [str(entry) for entry in reversed(denominator_poly.all_coeffs())]
    digest = canonical_sha256(
        {"numerator": numerator_coefficients, "denominator": denominator_coefficients}
    )
    return (
        int(numerator_poly.degree()),
        int(denominator_poly.degree()),
        len(numerator_poly.terms()),
        len(denominator_poly.terms()),
        digest,
    )


def independent_lucas_factors() -> tuple[tuple[tuple[str, sp.Expr], ...], dict[str, list[int]]]:
    tau, u = sp.symbols("tau u")
    lucas = [sp.Integer(2), tau]
    for _index in range(2, 10):
        lucas.append(sp.expand(tau * lucas[-1] - lucas[-2]))

    def quotient_in_u(index: int, divisor: sp.Expr) -> sp.Expr:
        quotient, remainder = sp.div(lucas[index], divisor, tau, domain=ZZ)
        if remainder:
            raise AssertionError(f"independent Lucas L_{index} division failed")
        expression = sp.Integer(0)
        for (power,), coefficient in sp.Poly(quotient, tau, domain=ZZ).terms():
            if power % 2:
                raise AssertionError(f"independent Lucas L_{index} quotient is not even")
            expression += coefficient * u ** (power // 2)
        return sp.expand(expression)

    factors = (
        ("r4", quotient_in_u(4, sp.Integer(1))),
        ("r5_over_r1", quotient_in_u(5, tau)),
        ("r6_over_r2", quotient_in_u(6, lucas[2])),
        ("r7_over_r1", quotient_in_u(7, tau)),
        ("r8", quotient_in_u(8, sp.Integer(1))),
        ("r9_over_r1", quotient_in_u(9, tau)),
    )
    records = {
        name: [int(coefficient) for coefficient in reversed(sp.Poly(factor, u, domain=ZZ).all_coeffs())]
        for name, factor in factors
    }
    return factors, records


def abstract_incidence_control() -> tuple[list[int], bool, int, str, str, dict[str, list[int]]]:
    u = sp.symbols("u")
    variables = sp.symbols("b7 b6 b5 b4 b3")
    b7, b6, b5, b4, b3 = variables
    targets = sp.symbols("s0 s2 s3 s4 s5 s6 s7 s8 s9")
    s0, s2, s3, s4, s5, s6, s7, s8, s9 = targets
    c2 = s2 - s0 - 256 - 128 * b7 - 64 * b6 - 32 * b5 - 16 * b4 - 8 * b3
    c3 = s3 - s0 - 6561 - 2187 * b7 - 729 * b6 - 243 * b5 - 81 * b4 - 27 * b3
    b2 = c3 / 3 - c2 / 2
    b1 = 3 * c2 / 2 - 2 * c3 / 3
    mode = (
        u**8 + b7 * u**7 + b6 * u**6 + b5 * u**5 + b4 * u**4 + b3 * u**3
        + b2 * u**2 + b1 * u + s0
    )
    independent_factors, factor_records = independent_lucas_factors()
    target_by_name = {
        "r4": s4,
        "r5_over_r1": s5,
        "r6_over_r2": s6,
        "r7_over_r1": s7,
        "r8": s8,
        "r9_over_r1": s9,
    }
    factors = [(factor, target_by_name[name]) for name, factor in independent_factors]
    field = QQ.frac_field(*targets)
    equations = [
        sp.Poly(sp.resultant(factor, mode, u) - target, *variables, domain=field)
        for factor, target in factors
    ]
    degrees = [int(equation.total_degree()) for equation in equations]
    leading = []
    for equation in equations[:5]:
        degree = equation.total_degree()
        expression = sum(
            coefficient.as_expr()
            * sp.prod(variable**power for variable, power in zip(variables, monomial))
            for monomial, coefficient in equation.terms()
            if sum(monomial) == degree
        )
        if set(expression.free_symbols) & set(targets):
            raise AssertionError("abstract leading form depends on trace targets")
        leading.append(sp.Poly(expression, *variables, domain=QQ))
    basis = sp.groebner(
        [entry.as_expr() for entry in leading],
        *variables,
        order="grevlex",
        domain=QQ,
    )
    leading_monomials = [
        entry.LM(order=basis.order).exponents for entry in basis.polys
    ]
    bounds = []
    for index in range(len(variables)):
        pure = [
            monomial[index]
            for monomial in leading_monomials
            if monomial[index] and sum(monomial) == monomial[index]
        ]
        if not pure:
            raise AssertionError(f"no pure leading power for variable {index}")
        bounds.append(min(pure))
    standard = [
        exponents
        for exponents in product(*(range(bound) for bound in bounds))
        if all(
            any(leading[index] > exponents[index] for index in range(5))
            for leading in leading_monomials
        )
    ]
    digest = canonical_sha256(
        [
            [[list(monomial), str(coefficient)] for monomial, coefficient in entry.terms()]
            for entry in basis.polys
        ]
    )
    return (
        degrees,
        bool(basis.is_zero_dimensional),
        len(standard),
        canonical_sha256(standard),
        digest,
        factor_records,
    )


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    failures: list[str] = []

    expected_hashes = {
        "experiments/e248_replica_trace_nine.py": file_sha256(PRODUCER),
        "tests/test_replica_trace_nine.py": file_sha256(SELF),
        "proofs/replica_trace_nine.md": file_sha256(PROOF),
        "experiments/e238_trace_exceptional_set.py": file_sha256(E238_PRODUCER),
        "results/spectral/trace_exceptional_set.json": file_sha256(E238_ARTIFACT),
    }
    if artifact["meta"]["source_sha256"] != expected_hashes:
        failures.append("source hash map mismatch")

    control = [independent_trace(3, power) for power in range(1, 8)]
    if coefficient_digest(control) != CONTROL_DIGEST_2X3:
        failures.append("independent 2x3 digest mismatch")

    rows = [independent_trace(4, power) for power in range(1, 10)]
    stored_rows = [[int(value) for value in row]
                   for row in artifact["data"]["trace_coefficients_ascending"]]
    if rows != stored_rows:
        failures.append("independent 2x4 coefficient reconstruction mismatch")
    if coefficient_digest(rows) != artifact["data"]["trace_coefficients_sha256"]:
        failures.append("2x4 coefficient digest mismatch")
    if [sum(value != 0 for value in row) for row in rows] != artifact["data"]["trace_term_counts"]:
        failures.append("2x4 term-count mismatch")
    observed_degrees = [
        max(index for index, coefficient in enumerate(row) if coefficient)
        for row in rows
    ]
    if observed_degrees != [18 * power for power in range(1, 10)]:
        failures.append(f"independent actual trace degrees mismatch: {observed_degrees}")
    if any(row[-1] != 2 for row in rows):
        failures.append("independent top trace coefficient is not two")

    targets = normalized_targets(rows)
    for name, value in targets.items():
        observed = rational_digest(value)
        record = artifact["data"]["normalized_targets"][name]
        expected = (
            record["numerator_degree"],
            record["denominator_degree"],
            record["numerator_term_count"],
            record["denominator_term_count"],
            record["sha256"],
        )
        if observed != expected:
            failures.append(f"normalized target mismatch: {name}")
    a7, a6, a5, a4, a3 = sp.symbols("a7 a6 a5 a4 a3")
    a0 = targets["r1_squared"]
    top_at_two = 256 + 128 * a7 + 64 * a6 + 32 * a5 + 16 * a4 + 8 * a3 + a0
    top_at_three = 6561 + 2187 * a7 + 729 * a6 + 243 * a5 + 81 * a4 + 27 * a3 + a0
    rhs_two = targets["r2"] - top_at_two
    rhs_three = targets["r3_over_r1"] - top_at_three
    solved_a2, solved_a1 = sp.Matrix([[4, 2], [9, 3]]).inv() * sp.Matrix(
        [rhs_two, rhs_three]
    )
    solved_a2 = sp.cancel(solved_a2)
    solved_a1 = sp.cancel(solved_a1)
    if (
        sp.cancel(4 * solved_a2 + 2 * solved_a1 - rhs_two) != 0
        or sp.cancel(9 * solved_a2 + 3 * solved_a1 - rhs_three) != 0
    ):
        failures.append("independent Q(0),Q(2),Q(3) reconstruction failed")
    mode_coefficients = (
        sp.Integer(1),
        a7,
        a6,
        a5,
        a4,
        a3,
        solved_a2,
        solved_a1,
        a0,
    )
    mode_incidence = artifact["data"]["mode_incidence"]
    if mode_incidence["a2_sha256"] != canonical_sha256(str(solved_a2)):
        failures.append("physical a2 digest mismatch")
    if mode_incidence["a1_sha256"] != canonical_sha256(str(solved_a1)):
        failures.append("physical a1 digest mismatch")
    if mode_incidence["mode_polynomial_coefficients_sha256"] != canonical_sha256(
        [str(coefficient) for coefficient in mode_coefficients]
    ):
        failures.append("physical mode-polynomial coefficient digest mismatch")
    if mode_incidence["coefficient_variables"] != ["a7", "a6", "a5", "a4", "a3"]:
        failures.append("mode coefficient-variable schema mismatch")
    if mode_incidence["first_rows"] != ["Q(0)=r1^2", "Q(2)=r2", "Q(3)=r3/r1"]:
        failures.append("mode first-row schema mismatch")
    if mode_incidence["remaining_rows"] != ["F4", "F5", "F6", "F7", "F8", "F9"]:
        failures.append("mode remaining-row schema mismatch")

    (
        degrees,
        zero_dimensional,
        quotient_rank,
        standard_digest,
        leading_digest,
        independent_factor_records,
    ) = abstract_incidence_control()
    leading_record = artifact["data"]["mode_incidence"]["leading_ideal"]
    if degrees != [2, 2, 2, 3, 4, 3]:
        failures.append(f"abstract equation degrees are {degrees}")
    if degrees != leading_record["equation_total_degrees"]:
        failures.append("artifact equation-degree mismatch")
    if zero_dimensional != leading_record["is_zero_dimensional"]:
        failures.append("leading-ideal dimension mismatch")
    if quotient_rank != 96 or quotient_rank != leading_record["quotient_basis_rank"]:
        failures.append(f"leading quotient rank mismatch: {quotient_rank}")
    if standard_digest != leading_record["standard_monomials_sha256"]:
        failures.append("standard-monomial digest mismatch")
    if leading_digest != leading_record["groebner_sha256"]:
        failures.append("leading Groebner digest mismatch")
    if mode_incidence["lucas_factor_coefficients_ascending"] != independent_factor_records:
        failures.append("independent Lucas factor recurrence mismatch")
    expected_quotient_theorem = {
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
    if mode_incidence["coefficient_quotient_theorem"] != expected_quotient_theorem:
        failures.append("coefficient-quotient theorem semantics mismatch")
    if mode_incidence["status"] != (
        "F4..F8 have generic Q(q) scheme length 96 by filtered/projective Bezout; "
        "an explicit lifted basis and the F9 multiplication norm are not materialized"
    ):
        failures.append("coefficient-quotient status mismatch")
    if not all(row and all(value >= 0 for value in row) and any(value > 0 for value in row) for row in rows):
        failures.append("physical denominator positivity prerequisite failed")
    expected_replica_identity = {
        "statement": "tr(B(q)^k) is the equality-enumerator of the 2 x L x directed-C_k replica graph and is evaluated by a 4^k-state column recurrence",
        "maximum_state_count": 4**9,
        "maximum_degree": 162,
        "uint64_partial_count_bound_bits": 54,
        "final_block_rows": 128,
        "uint64_final_block_sum_bound_bits": 61,
    }
    if artifact["data"]["replica_column_identity"] != expected_replica_identity:
        failures.append("replica-column identity/resource record mismatch")
    if 54 + ((128).bit_length() - 1) != 61 or 61 >= 64:
        failures.append("independent uint64 block-sum proof failed")

    expected_check_names = {
        "C1_dense_recurrence_control",
        "C2_landed_2x3_digest",
        "C3_trace_degrees_2x4",
        "C4_trace_configuration_counts",
        "C5_trace_nine_dimension_count",
        "C6_expected_equation_degrees",
        "C7_leading_complete_intersection_rank",
        "C8_filtered_complete_intersection_length",
        "C9_lucas_recurrence_factors",
    }
    if {check["name"] for check in artifact["checks"]} != expected_check_names:
        failures.append("producer check-name set mismatch")
    checks_by_name = {check["name"]: check for check in artifact["checks"]}
    if checks_by_name["C5_trace_nine_dimension_count"]["detail"] != (
        "five coefficient variables; F4..F8 define the coefficient ideal and "
        "F9 is the first additional row that can cut q"
    ):
        failures.append("trace-nine dimension scope detail mismatch")
    if artifact["data"]["scope"] != (
        "exact open 2x4 trace-nine incidence preflight; no exceptional-set emptiness claim"
    ):
        failures.append("artifact scope mismatch")

    if not all(check["passed"] for check in artifact["checks"]):
        failures.append("artifact contains a failed producer gate")
    elapsed = time.process_time() - started
    peak = peak_rss_bytes()
    if elapsed >= CPU_LIMIT_SECONDS:
        failures.append(f"verifier CPU {elapsed:.6f}s reached {CPU_LIMIT_SECONDS}s")
    if peak >= RSS_LIMIT_BYTES:
        failures.append(f"verifier RSS {peak} reached {RSS_LIMIT_BYTES}")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print(
        "PASS: e248 clean-room trace-nine reconstruction; "
        f"CPU={elapsed:.6f}s RSS={peak} bytes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
