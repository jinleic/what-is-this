#!/usr/bin/env python3
"""Clean-room verifier for the H676 trace-nine modular polynomial determinant."""
from __future__ import annotations

import copy
import ctypes
import hashlib
import importlib.util
import json
import math
import platform
import resource
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import flint
import sympy as sp
from flint import nmod_poly

ROOT = Path(__file__).resolve().parents[1]
SELF = ROOT / "tests" / "test_trace_nine_structured_determinant.py"
PRODUCER = ROOT / "experiments" / "e253_trace_nine_structured_determinant.py"
WAVE27_VERIFIER = ROOT / "tests" / "test_trace_nine_norm_envelope.py"
WAVE27_PRODUCER = ROOT / "experiments" / "e252_trace_nine_norm_envelope.py"
ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_structured_determinant.json"
WAVE27_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_norm_envelope.json"
E251_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
TRACE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
TEMPLATE = ROOT / "results" / "spectral" / "trace_nine_lift_template.json.xz"
LOCKFILE = ROOT / "uv.lock"
RUNTIME_FREEZE = ROOT.parent / "requirements-freeze.txt"
PRIME = 2_147_483_647
CPU_LIMIT_SECONDS = 1_800.0
RSS_LIMIT_BYTES = 2 * 1024**3
EVALUATION_CAP = 10_000_000
Monomial = tuple[int, int, int, int, int]

spec = importlib.util.spec_from_file_location("wave28_independent_wave27", WAVE27_VERIFIER)
wave27 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = wave27
spec.loader.exec_module(wave27)
e251 = wave27.base


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        class TimeValue(ctypes.Structure):
            _fields_ = [("seconds", ctypes.c_int32), ("microseconds", ctypes.c_int32)]

        class TaskBasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_time", TimeValue),
                ("system_time", TimeValue),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        library = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        task = library.mach_task_self()
        info = TaskBasicInfo()
        count = ctypes.c_uint32(ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_int32))
        if library.task_info(task, 20, ctypes.byref(info), ctypes.byref(count)):
            raise RuntimeError("mach task_info failed")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def guard_resources(started: float, stage: str) -> None:
    cpu_seconds = time.process_time() - started
    resident = peak_rss_bytes()
    if cpu_seconds >= CPU_LIMIT_SECONDS:
        raise RuntimeError(
            f"{stage}: process CPU {cpu_seconds:.3f}s reached {CPU_LIMIT_SECONDS:.0f}s"
        )
    if resident >= RSS_LIMIT_BYTES:
        raise RuntimeError(
            f"{stage}: peak RSS {resident} reached {RSS_LIMIT_BYTES} bytes"
        )


def zero_poly() -> nmod_poly:
    return nmod_poly([], PRIME)


def one_poly() -> nmod_poly:
    return nmod_poly([1], PRIME)


def monic(poly: nmod_poly) -> nmod_poly:
    if not poly:
        return poly
    return poly * pow(int(poly.leading_coefficient()), -1, PRIME)


@dataclass(frozen=True)
class ModRat:
    numerator: nmod_poly
    denominator: nmod_poly

    @staticmethod
    def make(numerator: nmod_poly, denominator: nmod_poly | None = None) -> "ModRat":
        denominator = one_poly() if denominator is None else denominator
        if not denominator:
            raise ZeroDivisionError
        if not numerator:
            return ModRat(zero_poly(), one_poly())
        common = numerator.gcd(denominator)
        numerator = numerator // common
        denominator = denominator // common
        inverse = pow(int(denominator.leading_coefficient()), -1, PRIME)
        return ModRat(numerator * inverse, denominator * inverse)

    @staticmethod
    def scalar(value: Fraction | int | sp.Rational) -> "ModRat":
        rational = Fraction(int(sp.numer(value)), int(sp.denom(value)))
        denominator = rational.denominator % PRIME
        if not denominator:
            raise ZeroDivisionError("scalar denominator vanishes modulo the prime")
        residue = rational.numerator % PRIME * pow(denominator, -1, PRIME) % PRIME
        return ModRat.make(nmod_poly([residue], PRIME))

    @staticmethod
    def from_sympy(expression: sp.Expr) -> "ModRat":
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        left, left_scale = wave27.zpoly_from_sympy(numerator)
        right, right_scale = wave27.zpoly_from_sympy(denominator)
        if left_scale % PRIME == 0 or right_scale % PRIME == 0:
            raise ZeroDivisionError("polynomial coefficient scale vanishes modulo the prime")
        return ModRat.make(
            nmod_poly([int(value) % PRIME for value in left], PRIME)
            * (right_scale % PRIME),
            nmod_poly([int(value) % PRIME for value in right], PRIME)
            * (left_scale % PRIME),
        )

    def __bool__(self) -> bool:
        return bool(self.numerator)

    def __neg__(self) -> "ModRat":
        return ModRat(-self.numerator, self.denominator)

    def __add__(self, other: "ModRat" | Fraction | int | sp.Rational) -> "ModRat":
        if not isinstance(other, ModRat):
            other = ModRat.scalar(other)
        common = self.denominator.gcd(other.denominator)
        left = self.denominator // common
        right = other.denominator // common
        return ModRat.make(
            self.numerator * right + other.numerator * left,
            left * other.denominator,
        )

    __radd__ = __add__

    def __sub__(self, other: "ModRat" | Fraction | int | sp.Rational) -> "ModRat":
        return self + (-other if isinstance(other, ModRat) else -ModRat.scalar(other))

    def __mul__(self, other: "ModRat" | Fraction | int | sp.Rational) -> "ModRat":
        if not isinstance(other, ModRat):
            other = ModRat.scalar(other)
        first = self.numerator.gcd(other.denominator)
        second = other.numerator.gcd(self.denominator)
        return ModRat.make(
            (self.numerator // first) * (other.numerator // second),
            (self.denominator // second) * (other.denominator // first),
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "ModRat":
        result = ModRat.scalar(1)
        base_value = self
        while exponent:
            if exponent & 1:
                result = result * base_value
            exponent //= 2
            if exponent:
                base_value = base_value * base_value
        return result

    def evaluate(self, value: int) -> int:
        numerator = int(self.numerator(value))
        denominator = int(self.denominator(value))
        if not denominator:
            raise ZeroDivisionError("rational matrix denominator vanishes at check point")
        return numerator * pow(denominator, -1, PRIME) % PRIME

    def canonical(self) -> list[list[int]]:
        return [list(map(int, self.numerator)), list(map(int, self.denominator))]


def evaluate_target_polynomial(terms, target_values: tuple[ModRat, ...]) -> ModRat:
    result = ModRat.scalar(0)
    for powers, coefficient in terms:
        term = ModRat.scalar(coefficient)
        for value, exponent in zip(target_values, powers):
            if exponent:
                term = term * (value**exponent)
        result = result + term
    return result


def reduce_scan(
    polynomial: dict[Monomial, ModRat],
    rules: list[tuple[Monomial, dict[Monomial, ModRat]]],
) -> tuple[dict[Monomial, ModRat], int]:
    remainder = {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if coefficient
    }
    steps = 0
    while True:
        candidates = [
            monomial
            for monomial in remainder
            if any(wave27.divides(leading, monomial) for leading, _ in rules)
        ]
        if not candidates:
            return remainder, steps
        monomial = max(candidates, key=wave27.grevlex_order)
        coefficient = remainder.pop(monomial)
        leading, tail = next(
            (leading, tail)
            for leading, tail in rules
            if wave27.divides(leading, monomial)
        )
        shift = tuple(a - b for a, b in zip(monomial, leading))
        for tail_monomial, tail_coefficient in tail.items():
            output = wave27.add_monomials(shift, tail_monomial)
            updated = remainder.get(output, ModRat.scalar(0)) - coefficient * tail_coefficient
            if updated:
                remainder[output] = updated
            else:
                remainder.pop(output, None)
        steps += 1


def polynomial_lcm(left: nmod_poly, right: nmod_poly) -> nmod_poly:
    if not left or not right:
        return zero_poly()
    return monic((left // left.gcd(right)) * right)


class NMod(ctypes.Structure):
    _fields_ = [
        ("n", ctypes.c_ulong),
        ("ninv", ctypes.c_ulong),
        ("norm", ctypes.c_ulong),
    ]


class NmodPolyStruct(ctypes.Structure):
    _fields_ = [
        ("coeffs", ctypes.POINTER(ctypes.c_ulong)),
        ("alloc", ctypes.c_long),
        ("length", ctypes.c_long),
        ("mod", NMod),
    ]


class NmodPolyMatStruct(ctypes.Structure):
    _fields_ = [
        ("entries", ctypes.POINTER(NmodPolyStruct)),
        ("r", ctypes.c_long),
        ("c", ctypes.c_long),
        ("stride", ctypes.c_long),
        ("modulus", ctypes.c_ulong),
    ]


def flint_library_path() -> Path:
    candidates = sorted((Path(flint.__file__).resolve().parent / ".dylibs").glob("libflint.*.dylib"))
    if len(candidates) != 1:
        raise AssertionError(f"expected one bundled FLINT library, found {candidates}")
    if flint.__FLINT_VERSION__ != "3.6.0" or flint.__FLINT_RELEASE__ != 30_600:
        raise AssertionError("native polynomial-matrix ABI requires FLINT 3.6.0")
    if ctypes.sizeof(NmodPolyStruct) != 48 or ctypes.sizeof(NmodPolyMatStruct) != 40:
        raise AssertionError("native FLINT structure layout changed")
    return candidates[0]


def native_entry(
    matrix: NmodPolyMatStruct, row: int, column: int
) -> ctypes.POINTER(NmodPolyStruct):
    index = row * matrix.stride + column
    return ctypes.cast(
        ctypes.byref(matrix.entries.contents, index * ctypes.sizeof(NmodPolyStruct)),
        ctypes.POINTER(NmodPolyStruct),
    )


def native_interpolation_determinant(rows: list[list[nmod_poly]]) -> nmod_poly:
    path = flint_library_path()
    library = ctypes.CDLL(str(path))
    library.nmod_poly_mat_init.argtypes = [
        ctypes.POINTER(NmodPolyMatStruct),
        ctypes.c_long,
        ctypes.c_long,
        ctypes.c_ulong,
    ]
    library.nmod_poly_mat_clear.argtypes = [ctypes.POINTER(NmodPolyMatStruct)]
    library.nmod_poly_init.argtypes = [ctypes.POINTER(NmodPolyStruct), ctypes.c_ulong]
    library.nmod_poly_clear.argtypes = [ctypes.POINTER(NmodPolyStruct)]
    library.nmod_poly_set_coeff_ui.argtypes = [
        ctypes.POINTER(NmodPolyStruct),
        ctypes.c_long,
        ctypes.c_ulong,
    ]
    library.nmod_poly_mat_det_interpolate.argtypes = [
        ctypes.POINTER(NmodPolyStruct),
        ctypes.POINTER(NmodPolyMatStruct),
    ]
    size = len(rows)
    matrix = NmodPolyMatStruct()
    determinant = NmodPolyStruct()
    library.nmod_poly_mat_init(ctypes.byref(matrix), size, size, PRIME)
    library.nmod_poly_init(ctypes.byref(determinant), PRIME)
    try:
        for row_index, row in enumerate(rows):
            for column_index, polynomial in enumerate(row):
                if not polynomial:
                    continue
                entry = native_entry(matrix, row_index, column_index)
                for degree, coefficient in enumerate(polynomial):
                    if coefficient:
                        library.nmod_poly_set_coeff_ui(entry, degree, int(coefficient))
        library.nmod_poly_mat_det_interpolate(
            ctypes.byref(determinant), ctypes.byref(matrix)
        )
        return nmod_poly(
            [int(determinant.coeffs[index]) for index in range(determinant.length)],
            PRIME,
        )
    finally:
        library.nmod_poly_clear(ctypes.byref(determinant))
        library.nmod_poly_mat_clear(ctypes.byref(matrix))


def source_hashes() -> dict[str, str]:
    workspace = ROOT.parent.parent
    paths = (
        PRODUCER,
        SELF,
        WAVE27_PRODUCER,
        WAVE27_VERIFIER,
        WAVE27_ARTIFACT,
        E251_ARTIFACT,
        TRACE_ARTIFACT,
        TEMPLATE,
        LOCKFILE,
        RUNTIME_FREEZE,
    )
    return {str(path.relative_to(workspace)): file_sha256(path) for path in paths}


def recompute() -> dict[str, object]:
    started = time.process_time()
    template, template_payload_sha256 = wave27.load_template()
    equations, _ = e251.independent_incidence()
    trace_artifact = json.loads(TRACE_ARTIFACT.read_text())
    trace_rows = [
        [int(value) for value in row]
        for row in trace_artifact["data"]["trace_coefficients_ascending"]
    ]
    target_values = tuple(
        ModRat.from_sympy(value) for value in e251.target_expressions(trace_rows).values()
    )
    compiled = e251.compile_rows(equations)
    physical = [
        {
            monomial: evaluate_target_polynomial(terms, target_values)
            for monomial, terms in row
        }
        for row in compiled
    ]
    leading_monomials = [tuple(record["leading_monomial"]) for record in template["records"]]
    pure_power_bounds, standard = wave27.standard_monomials(leading_monomials)
    if pure_power_bounds != [2, 2, 3, 5, 9] or len(standard) != 96:
        raise AssertionError("fixed leading quotient changed")
    remainders = []
    for index, equation in enumerate(equations[:5]):
        degree = equation.total_degree()
        remainders.append(
            {
                monomial: value
                for monomial, value in physical[index].items()
                if sum(monomial) < degree
            }
        )

    lifted_initial: list[tuple[Monomial, dict[Monomial, ModRat]]] = []
    for record in template["records"]:
        polynomial = {
            monomial: ModRat.scalar(coefficient)
            for monomial, coefficient in wave27.parse_q_polynomial(record["basis"]).items()
        }
        multipliers = [
            wave27.parse_q_polynomial(entries) for entries in record["transform"]
        ]
        for source_index, multiplier in enumerate(multipliers):
            for left_monomial, scalar in multiplier.items():
                for right_monomial, coefficient in remainders[source_index].items():
                    monomial = wave27.add_monomials(left_monomial, right_monomial)
                    polynomial[monomial] = (
                        polynomial.get(monomial, ModRat.scalar(0)) + coefficient * scalar
                    )
                    if not polynomial[monomial]:
                        del polynomial[monomial]
        leading = tuple(record["leading_monomial"])
        if polynomial.get(leading) != ModRat.scalar(1):
            raise AssertionError("independent modular lift lost a monic leader")
        lifted_initial.append(
            (
                leading,
                {
                    monomial: coefficient
                    for monomial, coefficient in polynomial.items()
                    if monomial != leading
                },
            )
        )

    reduced_by_index = {}
    processed_rules: list[tuple[Monomial, dict[Monomial, ModRat]]] = []
    for record_index in sorted(
        range(len(lifted_initial)),
        key=lambda index: (sum(lifted_initial[index][0]), lifted_initial[index][0]),
    ):
        leading, tail = lifted_initial[record_index]
        reduced_tail, _ = reduce_scan(tail, processed_rules)
        reduced_by_index[record_index] = (leading, reduced_tail)
        processed_rules.append((leading, reduced_tail))
    lifted = [reduced_by_index[index] for index in range(len(lifted_initial))]
    for equation_index, equation in enumerate(physical[:5]):
        remainder, _ = reduce_scan(equation, lifted)
        if remainder:
            raise AssertionError(
                f"physical equation {equation_index} has nonzero independent normal form"
            )
    guard_resources(started, "independent modular lift")

    row_index = {monomial: index for index, monomial in enumerate(standard)}
    rational_rows = [[ModRat.scalar(0) for _ in standard] for _ in standard]
    step_count = 0
    for column_index, basis_monomial in enumerate(standard):
        initial = {
            wave27.add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in physical[5].items()
        }
        reduced, steps = reduce_scan(initial, lifted)
        step_count += steps
        for monomial, coefficient in reduced.items():
            if monomial not in row_index:
                raise AssertionError(f"nonstandard independent matrix monomial {monomial}")
            rational_rows[row_index[monomial]][column_index] = coefficient
        guard_resources(started, f"independent matrix column {column_index}")

    cleared_rows = []
    numerator_scale = one_poly()
    denominator_scale = one_poly()
    row_denominator_degrees = []
    row_gcd_degrees = []
    row_max_degrees = []
    for row in rational_rows:
        denominator = one_poly()
        for value in row:
            if value:
                denominator = polynomial_lcm(denominator, value.denominator)
        cleared = [
            zero_poly()
            if not value
            else value.numerator * (denominator // value.denominator)
            for value in row
        ]
        row_gcd = zero_poly()
        for polynomial in cleared:
            if polynomial:
                row_gcd = polynomial if not row_gcd else row_gcd.gcd(polynomial)
        row_gcd = monic(row_gcd)
        primitive_row = [
            zero_poly() if not polynomial else polynomial // row_gcd
            for polynomial in cleared
        ]
        numerator_scale *= row_gcd
        denominator_scale *= denominator
        row_denominator_degrees.append(denominator.degree())
        row_gcd_degrees.append(row_gcd.degree())
        row_max_degrees.append(
            max(polynomial.degree() for polynomial in primitive_row if polynomial)
        )
        cleared_rows.append(primitive_row)
    matrix_digest = canonical_sha256(
        [
            canonical_sha256([entry.canonical() for entry in row])
            for row in rational_rows
        ]
    )
    matrix_nonzero = sum(bool(polynomial) for row in cleared_rows for polynomial in row)
    determinant = native_interpolation_determinant(cleared_rows)
    reduced_norm = ModRat.make(determinant * numerator_scale, denominator_scale)
    inherited = json.loads(E251_ARTIFACT.read_text())
    inherited_witness = next(
        witness
        for witness in inherited["data"]["modular_certificates"]["physical_witnesses"]
        if witness["q"] == "2"
    )
    direct_point_checks = []
    for q_value in (Fraction(2), Fraction(3), Fraction(5, 3)):
        q_residue = q_value.numerator * pow(q_value.denominator, -1, PRIME) % PRIME
        direct = e251.scalar_det(
            [[entry.evaluate(q_residue) for entry in row] for row in rational_rows],
            PRIME,
        )
        if direct != reduced_norm.evaluate(q_residue):
            raise AssertionError(f"independent point check fails at q={q_value}")
        if q_value == 2 and direct != inherited_witness["determinant_residue"]:
            raise AssertionError("independent q=2 check misses the inherited witness")
        direct_point_checks.append(
            {"q": str(q_value), "q_residue": q_residue, "determinant_residue": direct}
        )
    guard_resources(started, "independent native determinant")
    return {
        "template_payload_sha256": template_payload_sha256,
        "inherited_data_sha256": json.loads(WAVE27_ARTIFACT.read_text())["meta"]["data_sha256"],
        "standard_monomials": [list(monomial) for monomial in standard],
        "matrix_nonzero": matrix_nonzero,
        "step_count": step_count,
        "row_denominator_degrees": row_denominator_degrees,
        "row_gcd_degrees": row_gcd_degrees,
        "row_max_degrees": row_max_degrees,
        "matrix_digest": matrix_digest,
        "determinant": list(map(int, determinant)),
        "reduced_numerator": list(map(int, reduced_norm.numerator)),
        "reduced_denominator": list(map(int, reduced_norm.denominator)),
        "direct_point_checks": direct_point_checks,
        "cpu_seconds": time.process_time() - started,
    }


def validate_artifact(artifact: dict[str, object], recomputed: dict[str, object]) -> None:
    if artifact.get("schema") != "ising3d.trace-nine-structured-determinant/v1":
        raise AssertionError("unexpected H676 artifact schema")
    expected_meta = {
        "arithmetic": "exact F_2147483647(q) reduction and F_2147483647[q] determinant",
        "runtime": wave27.runtime_manifest(),
        "producer": str(PRODUCER.relative_to(ROOT)),
        "verifier": str(SELF.relative_to(ROOT)),
        "source_sha256": source_hashes(),
        "data_sha256": canonical_sha256(artifact["data"]),
    }
    if artifact["meta"] != expected_meta:
        raise AssertionError("H676 artifact metadata mismatch")
    data = artifact["data"]
    if data["field"] != {"prime": PRIME, "parameter": "q"}:
        raise AssertionError("H676 field changed")
    if data["template_payload_sha256"] != recomputed["template_payload_sha256"]:
        raise AssertionError("H676 template digest mismatch")
    matrix = data["matrix"]
    expected_matrix = {
        "shape": [96, 96],
        "standard_monomials": recomputed["standard_monomials"],
        "exact_nonzero_entries": recomputed["matrix_nonzero"],
        "ordered_reduction_steps": recomputed["step_count"],
        "row_denominator_degrees": recomputed["row_denominator_degrees"],
        "row_gcd_degrees": recomputed["row_gcd_degrees"],
        "row_max_entry_degrees": recomputed["row_max_degrees"],
        "row_denominator_degree_sum": sum(recomputed["row_denominator_degrees"]),
        "row_gcd_degree_sum": sum(recomputed["row_gcd_degrees"]),
        "row_degree_sum_bound": sum(recomputed["row_max_degrees"]),
        "exact_rational_matrix_sha256": recomputed["matrix_digest"],
    }
    if matrix != expected_matrix:
        raise AssertionError("H676 exact matrix record mismatch")
    method = data["structured_determinant"]
    determinant = recomputed["determinant"]
    maximum_entry_degree = max(recomputed["row_max_degrees"])
    expected_evaluations = 96 * maximum_entry_degree + 1
    if method["algorithm"] != (
        "FLINT 3.6.0 nmod_poly_mat_det_interpolate via pinned ctypes ABI"
    ):
        raise AssertionError("H676 determinant algorithm changed")
    if method["native_library_sha256"] != file_sha256(flint_library_path()):
        raise AssertionError("H676 native FLINT library digest mismatch")
    if method["maximum_entry_degree"] != maximum_entry_degree:
        raise AssertionError("H676 maximum matrix-entry degree mismatch")
    if method["scalar_determinant_evaluations"] != expected_evaluations:
        raise AssertionError("H676 scalar determinant count mismatch")
    if expected_evaluations > EVALUATION_CAP:
        raise AssertionError("H676 native interpolation exceeds the cap")
    if method["determinant_coefficients_ascending"] != determinant:
        raise AssertionError("H676 determinant coefficients mismatch")
    if method["determinant_degree"] != len(determinant) - 1:
        raise AssertionError("H676 determinant degree mismatch")
    if method["determinant_nonzero_coefficients"] != sum(bool(value) for value in determinant):
        raise AssertionError("H676 determinant support count mismatch")
    if method["determinant_coefficients_sha256"] != canonical_sha256(determinant):
        raise AssertionError("H676 determinant digest mismatch")
    norm = data["reduced_modular_norm"]
    numerator = recomputed["reduced_numerator"]
    denominator = recomputed["reduced_denominator"]
    expected_norm = {
        "numerator_degree": len(numerator) - 1,
        "denominator_degree": len(denominator) - 1,
        "numerator_coefficients_sha256": canonical_sha256(numerator),
        "denominator_coefficients_sha256": canonical_sha256(denominator),
        "numerator_coefficients_ascending": numerator,
        "denominator_coefficients_ascending": denominator,
        "status": "MATERIALIZED_MODULO_ONE_CERTIFIED_PRIME_ONLY",
    }
    if norm != expected_norm:
        raise AssertionError("H676 reduced modular norm mismatch")
    expected_mutations = {
        "rejected_count": 6,
        "records": [
            {
                "id": "M1_matrix_digest_last_hex_flip",
                "rejected": True,
                "falsified_claim": "exact modular matrix provenance",
            },
            {
                "id": "M2_determinant_constant_plus_one",
                "rejected": True,
                "falsified_claim": "full modular determinant coefficients",
            },
            {
                "id": "M3_determinant_degree_minus_one",
                "rejected": True,
                "falsified_claim": "materialized modular determinant degree",
            },
            {
                "id": "M4_reduced_denominator_constant_plus_one",
                "rejected": True,
                "falsified_claim": "reduced modular norm denominator",
            },
            {
                "id": "M5_scalar_evaluation_cap_minus_count",
                "rejected": True,
                "falsified_claim": "native interpolation scalar determinant count",
            },
            {
                "id": "M6_runtime_python_flint_0_6",
                "rejected": True,
                "falsified_claim": "installed runtime matches the frozen environment",
            },
        ],
    }
    if data["mutations"] != expected_mutations:
        raise AssertionError("H676 mutation transcript changed")
    expected_theorem = {
        "statement": (
            "The row-cleared rank-96 trace-nine multiplication determinant is materialized "
            "exactly over F_2147483647[q]; after restoring row scales and cancelling "
            "the denominator, the rational norm has a nonzero reduced modular numerator "
            "of degree 19846."
        ),
        "characteristic_zero_status": "NOT_RECONSTRUCTED",
        "physical_root_status": "UNRESOLVED",
        "scope": (
            "one exact good-prime reduction of the finite open-2x4 trace-nine norm; "
            "no characteristic-zero coefficient, real-root, endpoint, all-size, "
            "sampling, or thermodynamic claim"
        ),
    }
    if data["theorem"] != expected_theorem:
        raise AssertionError("H676 theorem scope changed")
    expected_checks = [
        {
            "name": "C0_inherited_wave27_provenance",
            "passed": True,
            "detail": recomputed["inherited_data_sha256"],
        },
        {
            "name": "C1_fixed_leading_modular_lift",
            "passed": True,
            "detail": "35 monic rules; F4-F8 normal forms are zero; rank 96",
        },
        {
            "name": "C2_exact_modular_multiplication_matrix",
            "passed": True,
            "detail": (
                f"{recomputed['matrix_nonzero']} entries; "
                f"{recomputed['step_count']} reductions; {recomputed['matrix_digest']}"
            ),
        },
        {
            "name": "C3_row_primitive_polynomial_matrix",
            "passed": True,
            "detail": (
                f"denominator degree {sum(recomputed['row_denominator_degrees'])}; "
                f"row gcd degree {sum(recomputed['row_gcd_degrees'])}"
            ),
        },
        {
            "name": "C4_native_flint_polynomial_determinant",
            "passed": True,
            "detail": (
                f"degree {len(determinant) - 1}; {expected_evaluations} exact scalar "
                f"determinants; FLINT {flint.__FLINT_VERSION__} dylib "
                f"{file_sha256(flint_library_path())}"
            ),
        },
        {
            "name": "C5_reduced_modular_norm",
            "passed": True,
            "detail": (
                f"numerator degree {len(numerator) - 1}; "
                f"denominator degree {len(denominator) - 1}"
            ),
        },
        {
            "name": "C6_direct_point_determinants",
            "passed": True,
            "detail": "; ".join(
                f"q={record['q']}: det={record['determinant_residue']}"
                for record in recomputed["direct_point_checks"]
            ),
        },
        {
            "name": "C7_claim_falsifying_mutations",
            "passed": True,
            "detail": "6/6 rejected",
        },
    ]
    if artifact["checks"] != expected_checks:
        raise AssertionError("H676 check transcript changed")


def reseal_data(artifact: dict[str, object]) -> None:
    artifact["meta"]["data_sha256"] = canonical_sha256(artifact["data"])


def mutation_rejections(
    artifact: dict[str, object], recomputed: dict[str, object]
) -> list[str]:
    mutations = []

    mutated = copy.deepcopy(artifact)
    digest = mutated["data"]["matrix"]["exact_rational_matrix_sha256"]
    mutated["data"]["matrix"]["exact_rational_matrix_sha256"] = digest[:-1] + (
        "0" if digest[-1] != "0" else "1"
    )
    reseal_data(mutated)
    mutations.append("M1" if rejects(mutated, recomputed) else "")

    mutated = copy.deepcopy(artifact)
    coefficients = mutated["data"]["structured_determinant"][
        "determinant_coefficients_ascending"
    ]
    coefficients[0] = (coefficients[0] + 1) % PRIME
    mutated["data"]["structured_determinant"][
        "determinant_coefficients_sha256"
    ] = canonical_sha256(coefficients)
    reseal_data(mutated)
    mutations.append("M2" if rejects(mutated, recomputed) else "")

    mutated = copy.deepcopy(artifact)
    mutated["data"]["structured_determinant"]["determinant_degree"] -= 1
    reseal_data(mutated)
    mutations.append("M3" if rejects(mutated, recomputed) else "")

    mutated = copy.deepcopy(artifact)
    coefficients = mutated["data"]["reduced_modular_norm"][
        "denominator_coefficients_ascending"
    ]
    coefficients[0] = (coefficients[0] + 1) % PRIME
    mutated["data"]["reduced_modular_norm"][
        "denominator_coefficients_sha256"
    ] = canonical_sha256(coefficients)
    reseal_data(mutated)
    mutations.append("M4" if rejects(mutated, recomputed) else "")

    mutated = copy.deepcopy(artifact)
    mutated["data"]["structured_determinant"]["scalar_determinant_evaluations"] -= 1
    reseal_data(mutated)
    mutations.append("M5" if rejects(mutated, recomputed) else "")

    mutated = copy.deepcopy(artifact)
    mutated["meta"]["runtime"]["distributions"]["python-flint"] = "0.6.0"
    mutations.append("M6" if rejects(mutated, recomputed) else "")
    return mutations


def rejects(artifact: dict[str, object], recomputed: dict[str, object]) -> bool:
    try:
        validate_artifact(artifact, recomputed)
    except AssertionError:
        return True
    return False


def main() -> int:
    artifact_bytes = ARTIFACT.read_bytes()
    artifact = json.loads(artifact_bytes)
    if artifact_bytes != canonical_bytes(artifact):
        raise AssertionError("H676 artifact is not canonical JSON")
    recomputed = recompute()
    validate_artifact(artifact, recomputed)

    if artifact["data"]["direct_point_checks"] != recomputed["direct_point_checks"]:
        raise AssertionError("H676 direct point checks mismatch")
    rejected = mutation_rejections(artifact, recomputed)
    if rejected != ["M1", "M2", "M3", "M4", "M5", "M6"]:
        raise AssertionError(f"H676 mutation rejection mismatch: {rejected}")
    print(
        "H676_TRACE_NINE_STRUCTURED_DETERMINANT_VERIFIER PASS "
        f"checks=8 mutations=6 numerator_degree={len(recomputed['reduced_numerator']) - 1} "
        f"cpu={recomputed['cpu_seconds']:.3f}s peak_rss={peak_rss_bytes()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
