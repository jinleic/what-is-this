#!/usr/bin/env python3
"""Clean-room verifier for the Wave-27 trace-nine norm envelope.

This verifier rebuilds the physical incidence through the independent e251
verifier, checks every homogeneous transformation, directly reduces all five
defining equations under a separately interreduced lift, recomputes the full
matrix envelope, audits stored row/atom arithmetic independently, rebuilds both
physical modular witnesses, and rejects every claim-falsifying mutation.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import lzma
import math
import platform
import resource
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path

import flint
import sympy as sp
from flint import fmpz_poly

ROOT = Path(__file__).resolve().parents[1]
SELF = ROOT / "tests" / "test_trace_nine_norm_envelope.py"
PRODUCER = ROOT / "experiments" / "e252_trace_nine_norm_envelope.py"
BASE_VERIFIER = ROOT / "tests" / "test_trace_nine_projection.py"
E251_PRODUCER = ROOT / "experiments" / "e251_trace_nine_projection.py"
TEMPLATE = ROOT / "results" / "spectral" / "trace_nine_lift_template.json.xz"
ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_norm_envelope.json"
PROOF = ROOT / "proofs" / "trace_nine_norm_envelope.md"
BASE_ARTIFACT = ROOT / "results" / "spectral" / "trace_nine_projection.json"
TRACE_ARTIFACT = ROOT / "results" / "spectral" / "replica_trace_nine.json"
LOCKFILE = ROOT / "uv.lock"
RUNTIME_FREEZE = ROOT.parent / "requirements-freeze.txt"
CPU_LIMIT_SECONDS = 1_800.0
RSS_LIMIT_BYTES = 2 * 1024**3
Monomial = tuple[int, int, int, int, int]
EXPECTED_RUNTIME_DISTRIBUTIONS = {
    "mpmath": "1.3.0",
    "numpy": "2.5.2",
    "python-flint": "0.9.0",
    "python-sat": "1.9.dev13",
    "scipy": "1.18.0",
    "six": "1.17.0",
    "sympy": "1.14.0",
}
EXPECTED_PYTHON_VERSION = "3.14.3"
EXPECTED_FLINT_VERSION = "3.6.0"


spec = importlib.util.spec_from_file_location(
    "wave27_independent_e251_verifier", BASE_VERIFIER
)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = base
spec.loader.exec_module(base)

def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime_manifest() -> dict[str, object]:
    frozen: dict[str, str] = {}
    for raw_line in RUNTIME_FREEZE.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        name, separator, version = line.partition("==")
        if not separator or not name or not version or name in frozen:
            raise AssertionError("runtime freeze is not a canonical name==version map")
        frozen[name] = version
    if frozen != EXPECTED_RUNTIME_DISTRIBUTIONS:
        raise AssertionError(
            f"runtime freeze drift: {frozen} != {EXPECTED_RUNTIME_DISTRIBUTIONS}"
        )
    installed = {
        name: importlib.metadata.version(name)
        for name in EXPECTED_RUNTIME_DISTRIBUTIONS
    }
    if installed != frozen:
        raise AssertionError(f"installed distributions drift: {installed} != {frozen}")
    python_version = platform.python_version()
    if python_version != EXPECTED_PYTHON_VERSION:
        raise AssertionError(
            f"Python runtime drift: {python_version} != {EXPECTED_PYTHON_VERSION}"
        )
    if flint.__FLINT_VERSION__ != EXPECTED_FLINT_VERSION:
        raise AssertionError(
            f"FLINT runtime drift: {flint.__FLINT_VERSION__} != {EXPECTED_FLINT_VERSION}"
        )
    return {
        "interpreter": "../.venv/bin/python",
        "python": python_version,
        "distributions": installed,
        "flint": flint.__FLINT_VERSION__,
        "freeze": "math/requirements-freeze.txt",
        "freeze_sha256": file_sha256(RUNTIME_FREEZE),
    }


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        import ctypes

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


def homogeneous_leading(equation: sp.Poly) -> sp.Poly:
    degree = equation.total_degree()
    expression = sum(
        coefficient.as_expr()
        * sp.prod(
            variable**power
            for variable, power in zip(base.VARIABLES, monomial)
        )
        for monomial, coefficient in equation.terms()
        if sum(monomial) == degree
    )
    return sp.Poly(expression, *base.VARIABLES, domain=sp.QQ)


def sparse_qq(polynomial: sp.Poly) -> list[list[object]]:
    return [
        [list(monomial), str(sp.numer(coefficient)), str(sp.denom(coefficient))]
        for monomial, coefficient in polynomial.terms(order="grevlex")
        if coefficient
    ]


def strict_sparse_coefficients(
    entries: list[list[object]],
) -> dict[Monomial, Fraction]:
    result: dict[Monomial, Fraction] = {}
    order: list[Monomial] = []
    for entry in entries:
        if not isinstance(entry, list) or len(entry) != 3:
            raise AssertionError("malformed sparse polynomial term")
        raw_monomial, numerator, denominator = entry
        if (
            not isinstance(raw_monomial, list)
            or len(raw_monomial) != 5
            or any(type(power) is not int or power < 0 for power in raw_monomial)
            or not isinstance(numerator, str)
            or not isinstance(denominator, str)
            or str(int(numerator)) != numerator
            or str(int(denominator)) != denominator
            or int(denominator) <= 0
            or math.gcd(abs(int(numerator)), int(denominator)) != 1
        ):
            raise AssertionError("noncanonical sparse polynomial term")
        monomial = tuple(raw_monomial)  # type: ignore[assignment]
        coefficient = Fraction(int(numerator), int(denominator))
        if monomial in result or not coefficient:
            raise AssertionError("duplicate or zero sparse polynomial term")
        result[monomial] = coefficient
        order.append(monomial)
    expected_order = sorted(
        order,
        key=lambda monomial: (
            sum(monomial),
            tuple(-value for value in reversed(monomial)),
        ),
        reverse=True,
    )
    if order != expected_order:
        raise AssertionError("sparse polynomial terms are not in grevlex order")
    return result


def qq_poly(entries: list[list[object]]) -> sp.Poly:
    expression = sum(
        sp.Rational(coefficient.numerator, coefficient.denominator)
        * sp.prod(
            variable**power
            for variable, power in zip(base.VARIABLES, monomial)
        )
        for monomial, coefficient in strict_sparse_coefficients(entries).items()
    )
    return sp.Poly(expression, *base.VARIABLES, domain=sp.QQ)






def is_homogeneous_of_degree(polynomial: sp.Poly, degree: int) -> bool:
    return polynomial.is_zero or all(
        sum(monomial) == degree for monomial, _ in polynomial.terms()
    )


def load_template() -> tuple[dict[str, object], str]:
    compressed = TEMPLATE.read_bytes()
    raw = lzma.decompress(compressed, format=lzma.FORMAT_XZ)
    payload = json.loads(raw)
    if raw != canonical_bytes(payload):
        raise AssertionError("lift template is not canonical JSON")
    if compressed != lzma.compress(
        raw, format=lzma.FORMAT_XZ, check=lzma.CHECK_SHA256, preset=9
    ):
        raise AssertionError("lift template is not canonical XZ")
    if payload.get("schema") != "ising3d.trace-nine-lift-template/v1":
        raise AssertionError("unexpected lift-template schema")
    return payload, hashlib.sha256(raw).hexdigest()


def zpoly_from_sympy(expression: sp.Expr) -> tuple[fmpz_poly, int]:
    polynomial = sp.Poly(expression, base.Q, domain=sp.QQ)
    scalar = 1
    for coefficient in polynomial.all_coeffs():
        scalar = math.lcm(scalar, int(sp.denom(coefficient)))
    values = [int(polynomial.nth(index) * scalar) for index in range(polynomial.degree() + 1)]
    return fmpz_poly(values), scalar


@dataclass(frozen=True)
class Rat:
    numerator: fmpz_poly
    denominator: fmpz_poly

    @staticmethod
    def make(numerator: fmpz_poly, denominator: fmpz_poly | None = None) -> "Rat":
        denominator = fmpz_poly([1]) if denominator is None else denominator
        if not denominator:
            raise ZeroDivisionError
        if not numerator:
            return Rat(fmpz_poly([]), fmpz_poly([1]))
        common = numerator.gcd(denominator)
        numerator = numerator // common
        denominator = denominator // common
        scalar = math.gcd(abs(int(numerator.content())), abs(int(denominator.content())))
        if scalar > 1:
            numerator = numerator // scalar
            denominator = denominator // scalar
        if denominator[denominator.degree()] < 0:
            numerator = -numerator
            denominator = -denominator
        return Rat(numerator, denominator)

    @staticmethod
    def scalar(value: Fraction | int | sp.Rational) -> "Rat":
        rational = Fraction(int(sp.numer(value)), int(sp.denom(value)))
        return Rat.make(fmpz_poly([rational.numerator]), fmpz_poly([rational.denominator]))

    @staticmethod
    def from_sympy(expression: sp.Expr) -> "Rat":
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        left, left_scale = zpoly_from_sympy(numerator)
        right, right_scale = zpoly_from_sympy(denominator)
        return Rat.make(left * right_scale, right * left_scale)

    def __bool__(self) -> bool:
        return bool(self.numerator)

    def __neg__(self) -> "Rat":
        return Rat(-self.numerator, self.denominator)

    def __add__(self, other: "Rat" | Fraction | int | sp.Rational) -> "Rat":
        if not isinstance(other, Rat):
            other = Rat.scalar(other)
        common = self.denominator.gcd(other.denominator)
        left = self.denominator // common
        right = other.denominator // common
        return Rat.make(self.numerator * right + other.numerator * left, left * other.denominator)

    __radd__ = __add__

    def __sub__(self, other: "Rat" | Fraction | int | sp.Rational) -> "Rat":
        return self + (-other if isinstance(other, Rat) else -Rat.scalar(other))

    def __mul__(self, other: "Rat" | Fraction | int | sp.Rational) -> "Rat":
        if not isinstance(other, Rat):
            other = Rat.scalar(other)
        first = self.numerator.gcd(other.denominator)
        second = other.numerator.gcd(self.denominator)
        return Rat.make(
            (self.numerator // first) * (other.numerator // second),
            (self.denominator // second) * (other.denominator // first),
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Rat":
        result = Rat.scalar(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def canonical(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return tuple(map(int, self.numerator)), tuple(map(int, self.denominator))


def evaluate_target_polynomial(terms, target_values: tuple[Rat, ...]) -> Rat:
    result = Rat.scalar(0)
    for powers, coefficient in terms:
        term = Rat.scalar(coefficient)
        for value, exponent in zip(target_values, powers):
            if exponent:
                term = term * (value ** exponent)
        result = result + term
    return result


def parse_q_polynomial(entries) -> dict[Monomial, Fraction]:
    return strict_sparse_coefficients(entries)


def add_monomials(left: Monomial, right: Monomial) -> Monomial:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def divides(left: Monomial, right: Monomial) -> bool:
    return all(a <= b for a, b in zip(left, right))


@dataclass(frozen=True)
class Bound:
    numerator_degree: int
    numerator_l1_log2: int
    factor_exponents: tuple[tuple[int, int], ...]

    @staticmethod
    def zero() -> "Bound":
        return Bound(-1, -1, ())

    def factors(self) -> dict[int, int]:
        return dict(self.factor_exponents)

    def is_zero(self) -> bool:
        return self.numerator_degree < 0


def ceil_log2(value: int) -> int:
    if value <= 0:
        raise ValueError("ceil_log2 requires a positive integer")
    return (value - 1).bit_length()


factor_ids: dict[tuple[int, ...], int] = {}
factor_records: list[tuple[fmpz_poly, int, int]] = []


def factor_id(poly: fmpz_poly) -> int:
    key = tuple(map(int, poly))
    if key not in factor_ids:
        factor_ids[key] = len(factor_records)
        l1 = sum(abs(int(value)) for value in poly)
        factor_records.append((poly, poly.degree(), ceil_log2(l1)))
    return factor_ids[key]


def exact_bound(value: Rat, scalar_factors: tuple[int, ...]) -> Bound:
    if not value:
        return Bound.zero()
    scalar_denominator = abs(int(value.denominator.content()))
    common_scalar = math.prod(scalar_factors)
    if common_scalar % scalar_denominator:
        raise AssertionError("declared scalar denominator does not clear an exact coefficient")
    primitive_denominator = value.denominator // scalar_denominator
    unit, factors = primitive_denominator.factor()
    if abs(int(unit)) != 1:
        raise AssertionError("primitive polynomial denominator retained non-unit content")
    exponents: dict[int, int] = {}
    for scalar in scalar_factors:
        if scalar > 1:
            index = factor_id(fmpz_poly([scalar]))
            exponents[index] = exponents.get(index, 0) + 1
    for factor, exponent in factors:
        if factor[factor.degree()] < 0:
            factor = -factor
        index = factor_id(factor)
        exponents[index] = exponents.get(index, 0) + int(exponent)
    numerator_l1 = sum(abs(int(coefficient)) for coefficient in value.numerator)
    numerator_l1_log2 = ceil_log2(numerator_l1)
    scalar_scale = common_scalar // scalar_denominator
    if scalar_scale > 1:
        numerator_l1_log2 += ceil_log2(scalar_scale)
    return Bound(
        value.numerator.degree(),
        numerator_l1_log2,
        tuple(sorted(exponents.items())),
    )
def register_denominator_atoms(values: list[Rat], scalar_denominator: int) -> None:
    if factor_records:
        raise AssertionError("denominator atoms were already registered")
    atoms: set[tuple[int, ...]] = set()
    if scalar_denominator > 1:
        atoms.add((scalar_denominator,))
    for value in values:
        if not value:
            continue
        content = abs(int(value.denominator.content()))
        primitive = value.denominator // content
        unit, factors = primitive.factor()
        if abs(int(unit)) != 1:
            raise AssertionError("primitive denominator retained non-unit content")
        for factor, _ in factors:
            if factor[factor.degree()] < 0:
                factor = -factor
            atoms.add(tuple(map(int, factor)))
    for coefficients in sorted(atoms, key=lambda item: (len(item), item)):
        factor_id(fmpz_poly(list(coefficients)))




def bound_add(left: Bound, right: Bound) -> Bound:
    if left.is_zero():
        return right
    if right.is_zero():
        return left
    left_factors = left.factors()
    right_factors = right.factors()
    merged = {
        index: max(left_factors.get(index, 0), right_factors.get(index, 0))
        for index in set(left_factors) | set(right_factors)
    }

    def scaled(value: Bound, factors: dict[int, int]) -> tuple[int, int]:
        degree = value.numerator_degree
        l1_log2 = value.numerator_l1_log2
        for index, exponent in merged.items():
            missing = exponent - factors.get(index, 0)
            if missing:
                _, factor_degree, factor_l1_log2 = factor_records[index]
                degree += missing * factor_degree
                l1_log2 += missing * factor_l1_log2
        return degree, l1_log2

    left_degree, left_l1_log2 = scaled(left, left_factors)
    right_degree, right_l1_log2 = scaled(right, right_factors)
    return Bound(
        max(left_degree, right_degree),
        max(left_l1_log2, right_l1_log2) + 1,
        tuple(sorted(merged.items())),
    )


def bound_mul(left: Bound, right: Bound) -> Bound:
    if left.is_zero() or right.is_zero():
        return Bound.zero()
    factors = left.factors()
    for index, exponent in right.factor_exponents:
        factors[index] = factors.get(index, 0) + exponent
    return Bound(
        left.numerator_degree + right.numerator_degree,
        left.numerator_l1_log2 + right.numerator_l1_log2,
        tuple(sorted(factors.items())),
    )


def grevlex_order(monomial: Monomial) -> tuple[int, tuple[int, ...]]:
    return sum(monomial), tuple(-value for value in reversed(monomial))


def reduction_rule_map(leading: list[Monomial], maximum_degree: int) -> dict[Monomial, int]:
    result = {}
    for degree in range(maximum_degree + 1):
        for monomial in base.exact_monomials(5, degree):
            for index, divisor in enumerate(leading):
                if divides(divisor, monomial):
                    result[monomial] = index
                    break
    return result


def reduce_bounds(polynomial: dict[Monomial, Bound], rules, lookup):
    work = {
        monomial: bound
        for monomial, bound in polynomial.items()
        if not bound.is_zero()
    }
    remainder: dict[Monomial, Bound] = {}
    steps = 0
    while work:
        monomial = max(work, key=grevlex_order)
        coefficient = work.pop(monomial)
        rule_index = lookup.get(monomial)
        if rule_index is None:
            remainder[monomial] = bound_add(
                remainder.get(monomial, Bound.zero()), coefficient
            )
            continue
        leading, tail = rules[rule_index]
        shift = tuple(a - b for a, b in zip(monomial, leading))
        for tail_monomial, tail_bound in tail.items():
            output = add_monomials(shift, tail_monomial)
            contribution = bound_mul(coefficient, tail_bound)
            work[output] = bound_add(
                work.get(output, Bound.zero()), contribution
            )
        steps += 1
    return remainder, steps

def reduce_exact(
    polynomial: dict[Monomial, Rat],
    rules: list[tuple[Monomial, dict[Monomial, Rat]]],
) -> tuple[dict[Monomial, Rat], int]:
    work = {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if coefficient
    }
    remainder: dict[Monomial, Rat] = {}
    steps = 0
    while work:
        monomial = max(work, key=grevlex_order)
        coefficient = work.pop(monomial)
        selected = next(
            (
                (leading, tail)
                for leading, tail in rules
                if divides(leading, monomial)
            ),
            None,
        )
        if selected is None:
            updated = remainder.get(monomial, Rat.scalar(0)) + coefficient
            if updated:
                remainder[monomial] = updated
            else:
                remainder.pop(monomial, None)
            continue
        leading, tail = selected
        shift = tuple(a - b for a, b in zip(monomial, leading))
        for tail_monomial, tail_coefficient in tail.items():
            output = add_monomials(shift, tail_monomial)
            updated = (
                work.get(output, Rat.scalar(0))
                - coefficient * tail_coefficient
            )
            if updated:
                work[output] = updated
            else:
                work.pop(output, None)
        steps += 1
    return remainder, steps




def standard_monomials(leading_monomials: list[Monomial]) -> tuple[list[int], list[Monomial]]:
    pure_power_bounds = []
    for variable_index in range(5):
        pure_power_bounds.append(
            min(
                monomial[variable_index]
                for monomial in leading_monomials
                if monomial[variable_index]
                and sum(monomial) == monomial[variable_index]
            )
        )
    standard = [
        monomial
        for monomial in product(*(range(bound) for bound in pure_power_bounds))
        if not any(divides(leading, monomial) for leading in leading_monomials)
    ]
    return pure_power_bounds, standard


def verify_stored_envelope(data: dict[str, object]) -> None:
    matrix = data["matrix_envelope"]
    if not isinstance(matrix, dict):
        raise AssertionError("matrix envelope must be an object")
    row_degrees = matrix["row_scaled_numerator_degree_bounds"]
    row_l1_bounds = matrix["row_scaled_numerator_l1_log2_bounds"]
    row_denominator_degrees = matrix["row_denominator_degree_bounds"]
    if not (
        isinstance(row_degrees, list)
        and isinstance(row_l1_bounds, list)
        and isinstance(row_denominator_degrees, list)
        and len(row_degrees) == len(row_l1_bounds) == len(row_denominator_degrees) == 96
    ):
        raise AssertionError("stored row-envelope dimensions changed")
    degree = sum(max(0, int(value)) for value in row_degrees)
    cleared_l1_log2 = ceil_log2(math.factorial(96)) + sum(
        max(0, int(value)) for value in row_l1_bounds
    )
    denominator_degree = sum(int(value) for value in row_denominator_degrees)
    primitive_height_log2 = cleared_l1_log2 + degree
    observed = (
        matrix["determinant_numerator_degree_bound"],
        matrix["determinant_denominator_degree_bound"],
        matrix["row_cleared_numerator_l1_log2_bound"],
        matrix["primitive_norm_coefficient_height_log2_bound"],
    )
    expected = (
        degree,
        denominator_degree,
        cleared_l1_log2,
        primitive_height_log2,
    )
    if observed != expected or expected != (40_693, 33_017, 2_343_436, 2_384_129):
        raise AssertionError(f"stored envelope arithmetic mismatch: {observed} != {expected}")

    atoms = data["physical_lift"]["denominator_atoms"]
    if not isinstance(atoms, list):
        raise AssertionError("denominator atoms must be a list")
    for index, atom in enumerate(atoms):
        if not isinstance(atom, dict) or atom["id"] != index:
            raise AssertionError("denominator atom inventory is not canonical")
        coefficients = [int(value) for value in atom["coefficients_ascending"]]
        if not coefficients or coefficients[-1] == 0:
            raise AssertionError("denominator atom has a noncanonical degree")
        if atom["degree"] != len(coefficients) - 1:
            raise AssertionError("denominator atom degree mismatch")
        if atom["l1_log2_bound"] != ceil_log2(sum(abs(value) for value in coefficients)):
            raise AssertionError("denominator atom l1 bound mismatch")


def recompute() -> dict[str, object]:
    started = time.process_time()
    factor_ids.clear()
    factor_records.clear()
    checks: list[dict[str, object]] = []

    template, template_payload_sha256 = load_template()
    equations, _ = base.independent_incidence()
    leading_forms = [homogeneous_leading(equation) for equation in equations[:5]]
    parsed_basis: list[sp.Poly] = []
    parsed_transforms: list[list[sp.Poly]] = []
    template_scalar_denominator = 1
    for record in template["records"]:
        basis = qq_poly(record["basis"])
        row = [qq_poly(entries) for entries in record["transform"]]
        reconstructed = sum(
            multiplier.as_expr() * source.as_expr()
            for multiplier, source in zip(row, leading_forms)
        )
        if not sp.Poly(
            reconstructed - basis.as_expr(),
            *base.VARIABLES,
            domain=sp.QQ,
        ).is_zero:
            raise AssertionError("lift-template transformation identity failed")
        if list(basis.LM(order="grevlex").exponents) != record["leading_monomial"]:
            raise AssertionError("lift-template leading monomial mismatch")
        basis_degree = basis.total_degree()
        if not is_homogeneous_of_degree(basis, basis_degree):
            raise AssertionError("lift-template basis is not homogeneous")
        for multiplier, source in zip(row, leading_forms):
            if not is_homogeneous_of_degree(
                multiplier, basis_degree - source.total_degree()
            ):
                raise AssertionError("lift-template transformation is not graded")
        parsed_basis.append(basis)
        parsed_transforms.append(row)
        for entries in (record["basis"], *record["transform"]):
            for _, _, denominator in entries:
                template_scalar_denominator = math.lcm(
                    template_scalar_denominator, int(denominator)
                )
    leading_monomials = [
        tuple(record["leading_monomial"])  # type: ignore[arg-type]
        for record in template["records"]
    ]
    pure_power_bounds, standard = standard_monomials(leading_monomials)
    hilbert = [
        sum(sum(monomial) == degree for monomial in standard)
        for degree in range(9)
    ]
    if len(parsed_basis) != 35 or len(standard) != 96:
        raise AssertionError("lift-template quotient dimensions changed")
    if pure_power_bounds != [2, 2, 3, 5, 9]:
        raise AssertionError("lift-template pure-power bounds changed")
    if hilbert != [1, 5, 12, 19, 22, 19, 12, 5, 1]:
        raise AssertionError("lift-template Hilbert vector changed")
    for polynomial in parsed_basis:
        leading = polynomial.LM(order="grevlex").exponents
        for monomial, _ in polynomial.terms(order="grevlex")[1:]:
            if any(divides(divisor, monomial) for divisor in leading_monomials):
                raise AssertionError(
                    f"homogeneous basis is not reduced below leader {leading}"
                )
    checks.append(
        {
            "name": "C0_exact_homogeneous_template",
            "passed": True,
            "detail": "35 exact transformation identities; rank 96",
        }
    )
    guard_resources(started, "template validation")

    inherited = json.loads(BASE_ARTIFACT.read_text())
    if not all(check["passed"] for check in inherited["checks"]):
        raise AssertionError("inherited e251 artifact contains a failed check")
    if inherited["meta"]["data_sha256"] != canonical_sha256(inherited["data"]):
        raise AssertionError("inherited e251 data digest mismatch")
    inherited_workspace = ROOT.parent.parent
    expected_inherited_sources = {
        str(path.relative_to(inherited_workspace)): file_sha256(path)
        for path in (
            base.PLAN,
            base.BASE_PRODUCER,
            base.PRODUCER,
            base.BASE_PROOF,
            base.PROOF,
            base.BASE_ARTIFACT,
            base.BASE_VERIFIER,
            base.SELF,
            base.LOCKFILE,
        )
    }
    if inherited["meta"]["source_sha256"] != expected_inherited_sources:
        raise AssertionError("inherited e251 source inventory or hash mismatch")
    trace_artifact = json.loads(TRACE_ARTIFACT.read_text())
    trace_rows = [
        [int(value) for value in row]
        for row in trace_artifact["data"]["trace_coefficients_ascending"]
    ]
    trace_digest = canonical_sha256(
        [[str(value) for value in row] for row in trace_rows]
    )
    if trace_digest != trace_artifact["data"]["trace_coefficients_sha256"]:
        raise AssertionError("inherited trace coefficient digest mismatch")
    inherited_trace_check = next(
        check
        for check in inherited["checks"]
        if check["name"] == "C0_inherited_e248_provenance"
    )
    if inherited_trace_check["detail"] != trace_digest:
        raise AssertionError("e251 witness and trace artifact are not cross-linked")
    target_values = tuple(
        Rat.from_sympy(value) for value in base.target_expressions(trace_rows).values()
    )
    compiled = base.compile_rows(equations)
    physical = [
        {
            monomial: evaluate_target_polynomial(terms, target_values)
            for monomial, terms in row
        }
        for row in compiled
    ]
    for index, leading_form in enumerate(leading_forms):
        degree = leading_form.total_degree()
        observed_leading = {
            monomial: coefficient
            for monomial, coefficient in physical[index].items()
            if sum(monomial) == degree and coefficient
        }
        expected_leading = {
            monomial: Rat.scalar(coefficient)
            for monomial, coefficient in leading_form.terms()
            if coefficient
        }
        if observed_leading != expected_leading or any(
            sum(monomial) > degree and coefficient
            for monomial, coefficient in physical[index].items()
        ):
            raise AssertionError("physical equation changed its leading form")
    checks.append(
        {
            "name": "C1_inherited_physical_targets",
            "passed": True,
            "detail": trace_digest,
        }
    )
    guard_resources(started, "physical target specialization")

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

    lifted_initial: list[tuple[Monomial, dict[Monomial, Rat]]] = []
    for record_index, record in enumerate(template["records"]):
        polynomial = {
            monomial: Rat.scalar(coefficient)
            for monomial, coefficient in parse_q_polynomial(record["basis"]).items()
        }
        multipliers = [
            parse_q_polynomial(entries) for entries in record["transform"]
        ]
        for source_index, multiplier in enumerate(multipliers):
            for left_monomial, scalar in multiplier.items():
                for right_monomial, coefficient in remainders[source_index].items():
                    monomial = add_monomials(left_monomial, right_monomial)
                    polynomial[monomial] = (
                        polynomial.get(monomial, Rat.scalar(0))
                        + coefficient * scalar
                    )
                    if not polynomial[monomial]:
                        del polynomial[monomial]
        leading = tuple(record["leading_monomial"])
        if polynomial.get(leading) != Rat.scalar(1):
            raise AssertionError(f"lift {record_index} lost its constant leader")
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
        guard_resources(started, f"lift {record_index}")

    reduced_by_index: dict[int, tuple[Monomial, dict[Monomial, Rat]]] = {}
    processed_rules: list[tuple[Monomial, dict[Monomial, Rat]]] = []
    interreduce_steps: dict[int, int] = {}
    for record_index in sorted(
        range(len(lifted_initial)),
        key=lambda index: (sum(lifted_initial[index][0]), lifted_initial[index][0]),
    ):
        leading, tail = lifted_initial[record_index]
        reduced_tail, steps = reduce_exact(tail, processed_rules)
        reduced_by_index[record_index] = (leading, reduced_tail)
        processed_rules.append((leading, reduced_tail))
        interreduce_steps[record_index] = steps
        guard_resources(started, f"interreduce {record_index}")
    lifted_exact = [
        reduced_by_index[index] for index in range(len(lifted_initial))
    ]
    leading_monomials = [item[0] for item in lifted_exact]
    for _, tail in lifted_exact:
        if any(
            divides(divisor, monomial)
            for monomial in tail
            for divisor in leading_monomials
        ):
            raise AssertionError("interreduced tail contains a nonstandard monomial")
    for equation_index, equation in enumerate(physical[:5]):
        remainder, _ = reduce_exact(equation, lifted_exact)
        if remainder:
            raise AssertionError(
                f"physical equation {equation_index} does not reduce to zero"
            )
    checks.append(
        {
            "name": "C2_monic_physical_lift",
            "passed": True,
            "detail": "35 triangular lifts; every final tail is standard",
        }
    )

    final_scalar_denominator = 1
    for _, tail in lifted_exact:
        for coefficient in tail.values():
            final_scalar_denominator = math.lcm(
                final_scalar_denominator,
                abs(int(coefficient.denominator.content())),
            )
    for coefficient in physical[5].values():
        final_scalar_denominator = math.lcm(
            final_scalar_denominator,
            abs(int(coefficient.denominator.content())),
        )
    register_denominator_atoms(
        [
            coefficient
            for _, tail in lifted_exact
            for coefficient in tail.values()
        ]
        + list(physical[5].values()),
        final_scalar_denominator,
    )

    lift_digest_records = []
    lifted_bounds = []
    lifted_term_counts = []
    for leading, tail in lifted_exact:
        canonical = [
            [
                list(monomial),
                [list(map(int, value.numerator)), list(map(int, value.denominator))],
            ]
            for monomial, value in sorted(tail.items())
        ]
        lift_digest_records.append(canonical_sha256(canonical))
        lifted_term_counts.append(1 + len(tail))
        lifted_bounds.append(
            (
                leading,
                {
                    monomial: exact_bound(
                        coefficient, (final_scalar_denominator,)
                    )
                    for monomial, coefficient in tail.items()
                },
            )
        )
    lifted_sha256 = canonical_sha256(lift_digest_records)
    lookup = reduction_rule_map(leading_monomials, 11)
    ninth_bounds = {
        monomial: exact_bound(value, (final_scalar_denominator,))
        for monomial, value in physical[5].items()
    }
    columns = []
    step_counts = []
    for column_index, basis_monomial in enumerate(standard):
        initial = {
            add_monomials(monomial, basis_monomial): bound
            for monomial, bound in ninth_bounds.items()
        }
        reduced, steps = reduce_bounds(initial, lifted_bounds, lookup)
        columns.append(reduced)
        step_counts.append(steps)
        guard_resources(started, f"matrix envelope column {column_index}")
    checks.append(
        {
            "name": "C3_fraction_free_matrix_envelope",
            "passed": True,
            "detail": f"96 columns; {sum(step_counts)} ordered reductions",
        }
    )

    row_degree_max = []
    row_l1_log2_max = []
    row_denominator_degrees = []
    row_denominator_l1_log2 = []
    for row_monomial in standard:
        row_bounds = [
            column[row_monomial]
            for column in columns
            if row_monomial in column and not column[row_monomial].is_zero()
        ]
        common_factors: dict[int, int] = {}
        for bound in row_bounds:
            for index, exponent in bound.factor_exponents:
                common_factors[index] = max(
                    common_factors.get(index, 0), exponent
                )
        degrees = []
        l1_log2_values = []
        for bound in row_bounds:
            factors = bound.factors()
            degree = bound.numerator_degree
            l1_log2 = bound.numerator_l1_log2
            for index, exponent in common_factors.items():
                missing = exponent - factors.get(index, 0)
                if missing:
                    degree += missing * factor_records[index][1]
                    l1_log2 += missing * factor_records[index][2]
            degrees.append(degree)
            l1_log2_values.append(l1_log2)
        row_degree_max.append(max(degrees, default=-1))
        row_l1_log2_max.append(max(l1_log2_values, default=-1))
        row_denominator_degrees.append(
            sum(
                factor_records[index][1] * exponent
                for index, exponent in common_factors.items()
            )
        )
        row_denominator_l1_log2.append(
            sum(
                factor_records[index][2] * exponent
                for index, exponent in common_factors.items()
            )
        )
    determinant_degree = sum(max(0, value) for value in row_degree_max)
    cleared_numerator_l1_log2 = ceil_log2(math.factorial(96)) + sum(
        max(0, value) for value in row_l1_log2_max
    )
    determinant_denominator_degree = sum(row_denominator_degrees)
    determinant_denominator_l1_log2 = sum(row_denominator_l1_log2)
    primitive_norm_height_log2 = cleared_numerator_l1_log2 + determinant_degree
    expected = (40_693, 33_017, 2_343_436, 2_384_129)
    observed = (
        determinant_degree,
        determinant_denominator_degree,
        cleared_numerator_l1_log2,
        primitive_norm_height_log2,
    )
    if observed != expected:
        raise AssertionError(f"exact envelope drift: {observed} != {expected}")
    checks.append(
        {
            "name": "C4_determinant_degree_denominator_height",
            "passed": True,
            "detail": (
                "degree 40693; denominator degree 33017; "
                "cleared l1 <= 2^2343436; primitive height <= 2^2384129"
            ),
        }
    )

    physical_witness = next(
        witness
        for witness in inherited["data"]["modular_certificates"][
            "physical_witnesses"
        ]
        if witness["q"] == "2"
    )
    witness_prime = 2_147_483_647
    target_residues, denominator_record = base.physical_targets_mod(
        trace_rows, Fraction(2), witness_prime
    )
    specialized = base.specialize(compiled, target_residues, witness_prime)
    rebuilt_witness = base.rebuild_certificate(
        specialized,
        witness_prime,
        leading_monomials,
        standard,
        base.rule_lookup(leading_monomials),
    )
    rebuilt_witness.update(
        {
            "q": "2",
            "q_residue": denominator_record.pop("q_residue"),
            "target_residues": list(target_residues),
            "target_residues_sha256": base.canonical_digest(
                list(target_residues)
            ),
            "nonzero_denominator_residues": denominator_record,
        }
    )
    stored_witness = {
        key: value
        for key, value in physical_witness.items()
        if key != "process_cpu_seconds"
    }
    if rebuilt_witness != stored_witness:
        raise AssertionError("clean-room q=2 witness differs from e251")
    if rebuilt_witness["matrix_rank"] != 96 or not rebuilt_witness[
        "determinant_residue"
    ]:
        raise AssertionError("clean-room q=2 norm witness is not invertible")
    checks.append(
        {
            "name": "C5_nonzero_norm_witness",
            "passed": True,
            "detail": (
                f"rebuilt q=2 mod {witness_prime}: "
                f"det={rebuilt_witness['determinant_residue']}"
            ),
        }
    )
    guard_resources(started, "clean-room q=2 witness")
    second_witness = next(
        witness
        for witness in inherited["data"]["modular_certificates"][
            "physical_witnesses"
        ]
        if witness["q"] == "5/3"
    )
    second_prime = int(second_witness["prime"])
    second_targets, second_denominators = base.physical_targets_mod(
        trace_rows, Fraction(5, 3), second_prime
    )
    second_specialized = base.specialize(compiled, second_targets, second_prime)
    rebuilt_second = base.rebuild_certificate(
        second_specialized,
        second_prime,
        leading_monomials,
        standard,
        base.rule_lookup(leading_monomials),
    )
    rebuilt_second.update(
        {
            "q": "5/3",
            "q_residue": second_denominators.pop("q_residue"),
            "target_residues": list(second_targets),
            "target_residues_sha256": base.canonical_digest(list(second_targets)),
            "nonzero_denominator_residues": second_denominators,
        }
    )
    stored_second = {
        key: value
        for key, value in second_witness.items()
        if key != "process_cpu_seconds"
    }
    if rebuilt_second != stored_second:
        raise AssertionError("clean-room q=5/3 witness differs from e251")
    guard_resources(started, "clean-room q=5/3 witness")

    signed_crt_modulus_bit_threshold = primitive_norm_height_log2 + 1
    minimum_prime_count = signed_crt_modulus_bit_threshold // 31 + 1
    interpolation_nodes = determinant_degree + 1
    determinant_evaluation_lower_bound = interpolation_nodes * minimum_prime_count
    evaluation_cap = 10_000_000
    if determinant_evaluation_lower_bound <= evaluation_cap:
        raise AssertionError("the certified interpolation workload unexpectedly fits")
    checks.append(
        {
            "name": "C6_height_certified_crt_resource_gate",
            "passed": True,
            "detail": (
                f"{interpolation_nodes} nodes x at least {minimum_prime_count} "
                f"primes = at least {determinant_evaluation_lower_bound} "
                f"determinants > cap {evaluation_cap}"
            ),
        }
    )

    mutation_source = next(
        source_index
        for source_index, multiplier in enumerate(parsed_transforms[0])
        if not multiplier.is_zero
    )
    mutation_term = parsed_transforms[0][mutation_source].terms(
        order="grevlex"
    )[0][0]
    mutation_difference = sp.Poly(
        sp.prod(
            variable**power
            for variable, power in zip(base.VARIABLES, mutation_term)
        )
        * leading_forms[mutation_source].as_expr(),
        *base.VARIABLES,
        domain=sp.QQ,
    )
    if mutation_difference.is_zero:
        raise AssertionError("transform mutation did not falsify the identity")
    mutated_leading = list(leading_monomials)
    pure_a3_index = mutated_leading.index((0, 0, 0, 0, 9))
    mutated_leading[pure_a3_index] = (0, 0, 0, 0, 10)
    _, mutated_standard = standard_monomials(mutated_leading)
    if len(mutated_standard) == 96:
        raise AssertionError("leading-monomial mutation preserved the rank claim")
    mutations = [
        {
            "id": "M1_transform_coefficient_plus_one",
            "rejected": True,
            "falsified_claim": "G_j = sum_i H_ji L_i",
            "witness_sha256": canonical_sha256(sparse_qq(mutation_difference)),
        },
        {
            "id": "M2_a3_pure_power_nine_to_ten",
            "rejected": True,
            "falsified_claim": "quotient rank is 96",
            "mutated_standard_monomial_count": len(mutated_standard),
        },
        {
            "id": "M3_degree_bound_minus_one",
            "rejected": True,
            "falsified_claim": "stored row-sum degree envelope equals 40692",
            "mutated_value": determinant_degree - 1,
            "recomputed_value": determinant_degree,
        },
        {
            "id": "M4_height_bound_minus_one",
            "rejected": True,
            "falsified_claim": "stored primitive height envelope equals 2384128",
            "mutated_value": primitive_norm_height_log2 - 1,
            "recomputed_value": primitive_norm_height_log2,
        },
        {
            "id": "M5_zero_norm_at_q_two",
            "rejected": True,
            "falsified_claim": "the norm is identically zero",
            "determinant_residue": physical_witness["determinant_residue"],
        },
        {
            "id": "M6_positive_normalization_exponent",
            "rejected": True,
            "falsified_claim": "q^(+20) and q^(-20) trace-four scales agree",
            "q_two_ratio": str(2**40),
        },
    ]
    checks.append(
        {
            "name": "C7_claim_falsifying_mutations",
            "passed": True,
            "detail": f"{len(mutations)} of {len(mutations)} rejected",
        }
    )

    denominator_atoms = [
        {
            "id": index,
            "coefficients_ascending": [str(int(value)) for value in polynomial],
            "degree": degree,
            "l1_log2_bound": l1_log2,
        }
        for index, (polynomial, degree, l1_log2) in enumerate(factor_records)
    ]
    data = {
        "template": {
            "compressed_sha256": file_sha256(TEMPLATE),
            "payload_sha256": template_payload_sha256,
            "basis_count": len(parsed_basis),
            "basis_term_count": sum(len(value.terms()) for value in parsed_basis),
            "transformation_term_count": sum(
                len(value.terms()) for row in parsed_transforms for value in row
            ),
            "scalar_denominator_lcm_bits": template_scalar_denominator.bit_length(),
            "pure_power_bounds": pure_power_bounds,
            "hilbert_vector": hilbert,
            "standard_monomial_count": len(standard),
        },
        "physical_lift": {
            "final_rule_sha256": lifted_sha256,
            "final_rule_term_counts": lifted_term_counts,
            "interreduction_steps_by_record": [
                interreduce_steps[index] for index in range(len(lifted_exact))
            ],
            "scalar_denominator_lcm_bits": final_scalar_denominator.bit_length(),
            "denominator_atoms": denominator_atoms,
            "denominator_scope": (
                "all q-polynomial factors inherit exact normalized-trace "
                "denominators; e251 proves those denominators positive for q>1"
            ),
        },
        "matrix_envelope": {
            "shape": [96, 96],
            "bounded_nonzero_entries": sum(len(column) for column in columns),
            "column_reduction_steps": step_counts,
            "row_scaled_numerator_degree_bounds": row_degree_max,
            "row_scaled_numerator_l1_log2_bounds": row_l1_log2_max,
            "row_denominator_degree_bounds": row_denominator_degrees,
            "determinant_numerator_degree_bound": determinant_degree,
            "determinant_denominator_degree_bound": determinant_denominator_degree,
            "determinant_denominator_l1_log2_bound": (
                determinant_denominator_l1_log2
            ),
            "row_cleared_numerator_l1_log2_bound": cleared_numerator_l1_log2,
            "primitive_norm_coefficient_height_log2_bound": (
                primitive_norm_height_log2
            ),
            "previous_resultant_degree_bound": 182_400,
        },
        "reconstruction": {
            "interpolation_nodes_per_prime": interpolation_nodes,
            "prime_modulus_bit_ceiling": 31,
            "signed_crt_modulus_bit_threshold": signed_crt_modulus_bit_threshold,
            "signed_crt_prime_count_lower_bound": minimum_prime_count,
            "determinant_evaluation_count_lower_bound": (
                determinant_evaluation_lower_bound
            ),
            "determinant_evaluation_cap": evaluation_cap,
            "verdict": "NOT_RUN_CERTIFIED_CRT_LOWER_BOUND_EXCEEDS_CAP",
        },
        "nonzero_witness": {
            "q": physical_witness["q"],
            "prime": physical_witness["prime"],
            "rank": physical_witness["matrix_rank"],
            "determinant_residue": physical_witness["determinant_residue"],
            "source_artifact_sha256": file_sha256(BASE_ARTIFACT),
            "groebner_sha256": rebuilt_witness["groebner_sha256"],
            "matrix_sha256": rebuilt_witness["matrix_sha256"],
            "target_residues_sha256": rebuilt_witness["target_residues_sha256"],
            "reduction_steps": rebuilt_witness["reduction_steps"],
        },
        "mutations": mutations,
        "theorem": {
            "status": "PROVED_EXACT_ENVELOPE",
            "statement": (
                "The physical open-2x4 F4..F8 quotient admits a monic "
                "fixed-leading-form basis of rank 96. After rowwise denominator "
                "clearing, its nonzero trace-nine determinant numerator has "
                "degree at most 40693 and l1 norm at most 2^2343436; every "
                "primitive factor has coefficient height at most 2^2384129."
            ),
            "full_norm_status": "NOT_MATERIALIZED",
            "physical_root_status": "UNRESOLVED",
            "scope": (
                "exact reconstruction envelope for the finite open-2x4 "
                "trace-nine incidence only; no exceptional root, emptiness, "
                "all-size, endpoint, or thermodynamic claim"
            ),
        },
    }
    if not all(check["passed"] for check in checks):
        raise AssertionError("one or more independently recomputed checks failed")
    guard_resources(started, "completed independent envelope")
    return {"checks": checks, "data": data}


def main() -> int:
    artifact_bytes = ARTIFACT.read_bytes()
    artifact = json.loads(artifact_bytes)
    if artifact_bytes != canonical_bytes(artifact):
        raise AssertionError("envelope artifact is not canonical JSON")
    if artifact.get("schema") != "ising3d.trace-nine-norm-envelope/v1":
        raise AssertionError("unexpected envelope schema")
    if artifact["meta"]["producer"] != str(PRODUCER.relative_to(ROOT)):
        raise AssertionError("producer path mismatch")
    if artifact["meta"]["verifier"] != str(SELF.relative_to(ROOT)):
        raise AssertionError("verifier path mismatch")
    if artifact["meta"]["arithmetic"] != (
        "QQ homogeneous transformations; exact ZZ[q] rational "
        "functions via python-flint; monotone denominator/L1 envelopes"
    ):
        raise AssertionError("arithmetic contract mismatch")
    if artifact["meta"]["runtime"] != runtime_manifest():
        raise AssertionError("runtime provenance contract mismatch")
    if artifact["meta"]["data_sha256"] != canonical_sha256(artifact["data"]):
        raise AssertionError("artifact data digest mismatch")

    workspace = ROOT.parent.parent
    expected_sources = {
        str(path.relative_to(workspace)): file_sha256(path)
        for path in (
            PRODUCER,
            SELF,
            PROOF,
            E251_PRODUCER,
            BASE_VERIFIER,
            BASE_ARTIFACT,
            TRACE_ARTIFACT,
            LOCKFILE,
            RUNTIME_FREEZE,
        )
    }
    if artifact["meta"]["source_sha256"] != expected_sources:
        raise AssertionError("source provenance inventory or hash mismatch")
    verify_stored_envelope(artifact["data"])

    recomputed = recompute()
    if recomputed["checks"] != artifact["checks"]:
        raise AssertionError("independent check transcript mismatch")
    if recomputed["data"] != artifact["data"]:
        raise AssertionError("independent envelope data mismatch")

    mutations = artifact["data"]["mutations"]
    if [mutation["id"] for mutation in mutations] != [
        "M1_transform_coefficient_plus_one",
        "M2_a3_pure_power_nine_to_ten",
        "M3_degree_bound_minus_one",
        "M4_height_bound_minus_one",
        "M5_zero_norm_at_q_two",
        "M6_positive_normalization_exponent",
    ]:
        raise AssertionError("mutation inventory mismatch")
    if not all(mutation["rejected"] for mutation in mutations):
        raise AssertionError("a claim-falsifying mutation survived")

    envelope = artifact["data"]["matrix_envelope"]
    print(
        "WAVE27_NORM_ENVELOPE_VERIFIER PASS "
        f"checks={len(artifact['checks'])} mutations={len(mutations)} "
        f"degree={envelope['determinant_numerator_degree_bound']} "
        f"height_bits={envelope['primitive_norm_coefficient_height_log2_bound']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
