#!/usr/bin/env python3
"""Independent verifier for the e242 trace-branch disposition theorem.

The quotient stage uses e238's independent degree-filtered Macaulay reducer,
not the producer's lifted H-basis reducer.  The disposition stage then rebuilds
all dyadic enclosures, Hermite signs, and the shifted Descartes bound from the
verified cofactor polynomials.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import platform
import resource
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "trace_branch_disposition.json"
BASE_ARTIFACT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
BASE_VERIFIER = ROOT / "tests" / "test_trace_exceptional_set.py"
PRODUCER = ROOT / "experiments" / "e242_trace_branch_disposition.py"
RSS_LIMIT_BYTES = 2 * 1024**3
QUOTIENT_CPU_LIMIT = 1800.0
DISPOSITION_CPU_LIMIT = 600.0
FAILURES: list[str] = []

spec = importlib.util.spec_from_file_location("e238_independent_verifier", BASE_VERIFIER)
if spec is None or spec.loader is None:
    raise RuntimeError("could not load the e238 independent verifier")
independent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(independent)


class Interval:
    __slots__ = ("lo", "hi")

    def __init__(self, lo: Fraction | int, hi: Fraction | int | None = None):
        self.lo = Fraction(lo)
        self.hi = self.lo if hi is None else Fraction(hi)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    def __add__(self, other: Interval | Fraction | int) -> Interval:
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: Interval | Fraction | int) -> Interval:
        other = other if isinstance(other, Interval) else Interval(other)
        return self + (-other)

    def __rsub__(self, other: Interval | Fraction | int) -> Interval:
        return Interval(other) - self

    def __mul__(self, other: Interval | Fraction | int) -> Interval:
        other = other if isinstance(other, Interval) else Interval(other)
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other: Interval | Fraction | int) -> Interval:
        other = other if isinstance(other, Interval) else Interval(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        return self * Interval(
            min(Fraction(1, other.lo), Fraction(1, other.hi)),
            max(Fraction(1, other.lo), Fraction(1, other.hi)),
        )

    def __rtruediv__(self, other: Interval | Fraction | int) -> Interval:
        return Interval(other) / self

    def sign(self) -> int:
        return 1 if self.lo > 0 else -1 if self.hi < 0 else 0


PolyData = tuple[int, tuple[int, ...]]


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def parse_interval(values: list[str]) -> Interval:
    return Interval(parse_fraction(values[0]), parse_fraction(values[1]))


def integer_poly_data(value: Any) -> PolyData:
    expression = value.as_expr() if hasattr(value, "as_expr") else value
    polynomial = sp.Poly(expression, independent.QPAR, domain=QQ)
    scale, integer = polynomial.clear_denoms(convert=True)
    return int(scale), tuple(int(coefficient) for coefficient in integer.all_coeffs())

def qpoly_coefficients_ascending(value: Any) -> list[str]:
    expression = value.as_expr() if hasattr(value, "as_expr") else value
    polynomial = sp.Poly(expression, independent.QPAR, domain=QQ)
    return [str(coefficient) for coefficient in reversed(polynomial.all_coeffs())]


def cofactor_record(value: Any) -> dict[str, object]:
    scale, descending = integer_poly_data(value)
    ascending = [str(coefficient) for coefficient in reversed(descending)]
    return {
        "degree": len(descending) - 1,
        "term_count": sum(coefficient != 0 for coefficient in descending),
        "integer_scale_denominator": str(scale),
        "integer_coefficients_ascending": ascending,
        "sha256": canonical_sha256(ascending),
    }


def artifact_poly_data(record: dict[str, object]) -> PolyData:
    ascending = [int(value) for value in record["integer_coefficients_ascending"]]
    return int(record["integer_scale_denominator"]), tuple(reversed(ascending))


def sign_at_dyadic(coefficients_ascending: list[int], numerator: int, bits: int) -> int:
    accumulator = coefficients_ascending[-1]
    denominator_power = 1 << bits
    power = denominator_power
    for coefficient in reversed(coefficients_ascending[:-1]):
        accumulator = accumulator * numerator + coefficient * power
        power *= denominator_power
    return (accumulator > 0) - (accumulator < 0)


def polynomial_interval(data: PolyData, numerator: int, bits: int) -> Interval:
    scale, coefficients = data
    degree = len(coefficients) - 1
    midpoint_numerator = 2 * numerator + 1
    midpoint_denominator = 1 << (bits + 1)
    accumulator = coefficients[0]
    denominator_power = midpoint_denominator
    for coefficient in coefficients[1:]:
        accumulator = accumulator * midpoint_numerator + coefficient * denominator_power
        denominator_power *= midpoint_denominator
    midpoint = Fraction(accumulator, denominator_power // midpoint_denominator) / scale
    if degree == 0:
        return Interval(midpoint)
    derivative_weight = sum(
        (degree - index) * abs(coefficient)
        for index, coefficient in enumerate(coefficients[:-1])
    )
    error = Fraction(
        derivative_weight * pow(numerator + 1, degree - 1),
        1 << (bits * degree + 1),
    ) / scale
    return Interval(midpoint - error, midpoint + error)


def outward(interval: Interval, bits: int = 220) -> Interval:
    magnitude = max(abs(interval.lo), abs(interval.hi))
    if magnitude == 0:
        return Interval(0)
    exponent = magnitude.numerator.bit_length() - magnitude.denominator.bit_length()
    shift = bits - exponent
    if shift >= 0:
        denominator = 1 << shift
        lower = (interval.lo.numerator * denominator) // interval.lo.denominator
        upper = -((-interval.hi.numerator * denominator) // interval.hi.denominator)
        return Interval(Fraction(lower, denominator), Fraction(upper, denominator))
    unit = 1 << (-shift)
    lower = (interval.lo.numerator // (interval.lo.denominator * unit)) * unit
    upper = -((-interval.hi.numerator // (interval.hi.denominator * unit))) * unit
    return Interval(lower, upper)


def rational_data(expression: sp.Expr) -> tuple[PolyData, PolyData]:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    return integer_poly_data(numerator), integer_poly_data(denominator)


def rational_interval(
    data: tuple[PolyData, PolyData], numerator: int, bits: int
) -> Interval:
    return outward(
        polynomial_interval(data[0], numerator, bits)
        / polynomial_interval(data[1], numerator, bits)
    )


def independent_cleared_matrix(
    equations: list[sp.Poly],
) -> tuple[list[list[Any]], list[sp.Poly]]:
    homogeneous = independent.leading_quadrics(equations)
    plans = independent.macaulay_plans(homogeneous, maximum_degree=6)
    seventh = equations[3].rep.to_dict()
    columns: list[list[Any]] = []
    for basis_monomial in independent.BASIS_MONOMIALS:
        product = {
            independent.add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in seventh.items()
        }
        reduced = independent.macaulay_reduce(product, equations, homogeneous, plans)
        columns.append(
            [
                reduced.get(monomial, independent.ZERO)
                for monomial in independent.BASIS_MONOMIALS
            ]
        )
    rows = [[columns[column][row] for column in range(8)] for row in range(8)]
    denominators: list[sp.Poly] = []
    for column in range(8):
        denominator = sp.Integer(1)
        for row in range(8):
            denominator = sp.lcm(denominator, rows[row][column].denom.as_expr())
        denominators.append(
            sp.Poly(denominator, independent.QPAR, domain=QQ).monic()
        )
    cleared: list[list[Any]] = [[None] * 8 for _ in range(8)]
    for column, denominator in enumerate(denominators):
        multiplier = independent.QQ_POLY_Q.from_sympy(denominator.as_expr())
        for row in range(8):
            value = rows[row][column] * independent.QQQ.new(multiplier)
            if value.denom.degree() != 0:
                raise AssertionError("independent column clearing failed")
            constant = next(iter(value.denom.to_dict().values()))
            cleared[row][column] = value.numer / constant
    return cleared, denominators


def independent_cofactor(rows: list[list[Any]], basis_index: int) -> Any:
    submatrix = [
        [rows[row][column] for column in range(1, 8)]
        for row in range(8)
        if row != basis_index
    ]
    determinant = DomainMatrix(
        submatrix, (7, 7), independent.QQ_POLY_Q
    ).det()
    return -determinant if basis_index % 2 else determinant


def interval_determinant(matrix: list[list[Interval]]) -> Interval:
    work = [[Interval(entry.lo, entry.hi) for entry in row] for row in matrix]
    dimension = len(work)
    parity = 1
    for column in range(dimension):
        pivots = [row for row in range(column, dimension) if work[row][column].sign()]
        if not pivots:
            raise ArithmeticError(f"no certified pivot at {column}")
        pivot = max(
            pivots,
            key=lambda row: min(
                abs(work[row][column].lo), abs(work[row][column].hi)
            ),
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            parity = -parity
        for row in range(column + 1, dimension):
            multiplier = work[row][column] / work[column][column]
            for offset in range(column + 1, dimension):
                work[row][offset] = (
                    work[row][offset] - multiplier * work[column][offset]
                )
            work[row][column] = Interval(0)
    result = Interval(parity)
    for index in range(dimension):
        result = result * work[index][index]
    return result


def power_sums(coefficients: list[Interval]) -> list[Interval]:
    degree = len(coefficients) - 1
    result = [Interval(degree)]
    for power in range(1, 2 * degree - 1):
        value = Interval(0)
        if power <= degree:
            for index in range(1, power):
                value = value + coefficients[index] * result[power - index]
            value = value + power * coefficients[power]
        else:
            for index in range(1, degree + 1):
                value = value + coefficients[index] * result[power - index]
        result.append(-value)
    return result


def hermite_signs(coefficients: list[Interval]) -> tuple[list[int], int]:
    sums = power_sums(coefficients)
    matrix = [[sums[row + column] for column in range(6)] for row in range(6)]
    signs: list[int] = []
    for dimension in range(1, 7):
        minor = outward(
            interval_determinant([row[:dimension] for row in matrix[:dimension]]),
            bits=180,
        )
        signs.append(minor.sign())
    changes = sum(
        left != right for left, right in zip([1] + signs, signs)
    )
    return signs, 6 - 2 * changes


def shifted_signs(coefficients: list[Interval]) -> tuple[list[int], int]:
    ascending = list(reversed(coefficients))
    shifted_ascending: list[Interval] = []
    for target in range(7):
        value = Interval(0)
        for source in range(target, 7):
            value = value + (
                ascending[source]
                * math.comb(source, target)
                * 4 ** (source - target)
            )
        shifted_ascending.append(outward(value, bits=180))
    signs = [entry.sign() for entry in reversed(shifted_ascending)]
    nonzero = [sign for sign in signs if sign]
    variations = sum(left != right for left, right in zip(nonzero, nonzero[1:]))
    return signs, variations


def build_trace_data() -> tuple[dict[str, sp.Expr], list[sp.Poly]]:
    epsilon, shift, traces = independent.traces_over_zzq(independent.BONDS_2X3)
    targets = independent.normalized_targets(traces, epsilon, shift)
    equations, _, _ = independent.norm_equations(targets)
    return targets, equations


def main() -> int:
    stage = os.environ.get("TRACE_BRANCH_STAGE", "all")
    if stage not in {"all", "quotient", "disposition"}:
        raise ValueError("TRACE_BRANCH_STAGE must be all, quotient, or disposition")
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    base_artifact = json.loads(BASE_ARTIFACT.read_text())
    quotient = artifact["data"]["quotient_certificate"]
    targets: dict[str, sp.Expr] | None = None
    equations: list[sp.Poly] | None = None

    check(
        "current producer hash",
        artifact["meta"]["source_sha256"]["experiments/e242_trace_branch_disposition.py"]
        == file_sha256(PRODUCER),
    )
    check(
        "current verifier hash",
        artifact["meta"]["source_sha256"]["tests/test_trace_branch_disposition.py"]
        == file_sha256(Path(__file__).resolve()),
    )
    check(
        "current e238 dependency hash",
        artifact["meta"]["dependency_sha256"][
            "results/spectral/trace_exceptional_set.json"
        ]
        == file_sha256(BASE_ARTIFACT),
    )

    if stage in {"all", "quotient"}:
        targets, equations = build_trace_data()
        rows, denominators = independent_cleared_matrix(equations)
        denominator_degrees = [denominator.degree() for denominator in denominators]
        check(
            "independent column denominator degrees",
            denominator_degrees == quotient["column_denominator_degrees"],
            str(denominator_degrees),
        )
        check(
            "independent column denominator polynomials",
            all(
                all(
                    cofactor_record(denominators[index])[key]
                    == quotient["column_denominators"][index][key]
                    for key in (
                        "degree",
                        "term_count",
                        "integer_scale_denominator",
                        "integer_coefficients_ascending",
                        "sha256",
                    )
                )
                for index in range(8)
            ),
        )
        matrix_digest = canonical_sha256(
            [
                [
                    qpoly_coefficients_ascending(rows[row][column])
                    for column in range(8)
                ]
                for row in range(8)
            ]
        )
        check(
            "independent cleared matrix digest",
            matrix_digest
            == quotient["matrix_sha256"]
            == quotient["inherited_matrix_sha256"],
            matrix_digest,
        )
        for basis_index in (0, 1, 4, 6):
            observed = cofactor_record(independent_cofactor(rows, basis_index))
            stored = quotient["cofactors"][str(basis_index)]
            check(
                f"independent cofactor {basis_index}",
                all(
                    observed[key] == stored[key]
                    for key in (
                        "degree",
                        "term_count",
                        "integer_scale_denominator",
                        "integer_coefficients_ascending",
                        "sha256",
                    )
                ),
                f"degree={observed['degree']} sha256={observed['sha256']}",
            )
        quotient_cpu = time.process_time() - started
        check(
            "quotient-stage CPU wall",
            quotient_cpu < QUOTIENT_CPU_LIMIT,
            f"{quotient_cpu:.6f}s/{QUOTIENT_CPU_LIMIT}s",
        )

    if stage in {"all", "disposition"}:
        if targets is None or equations is None:
            targets, equations = build_trace_data()
        exceptional = base_artifact["data"]["elimination"][
            "exceptional_q_polynomial"
        ]
        exceptional_coefficients = [int(value) for value in exceptional["coefficients"]]
        check(
            "inherited E digest",
            exceptional["sha256"]
            == canonical_sha256(exceptional["coefficients"])
            == artifact["data"]["inherited_exceptional_set"][
                "coefficient_sha256"
            ],
        )
        target_data = {
            name: rational_data(targets[name])
            for name in ("r1_squared", "r2", "r3_over_r1")
        }
        denominator_data = [
            artifact_poly_data(record) for record in quotient["column_denominators"]
        ]
        branches = artifact["data"]["branches"]
        inherited_branches = base_artifact["data"]["physical_root_certificate"][
            "roots"
        ]
        check(
            "three aligned inherited branches",
            len(branches) == len(inherited_branches) == 3
            and [int(branch["root_index"]) for branch in branches] == [1, 2, 3],
        )
        signatures: list[int] = []
        third_shifted: tuple[list[int], int] | None = None
        for branch, inherited_branch in zip(branches, inherited_branches):
            bits = int(branch["dyadic_bits"])
            numerator = int(branch["lower_numerator"])
            endpoint_signs = [
                sign_at_dyadic(exceptional_coefficients, numerator + offset, bits)
                for offset in (0, 1)
            ]
            check(
                f"branch {branch['root_index']} exact E sign change",
                endpoint_signs == branch["E_endpoint_signs"]
                and endpoint_signs[0] * endpoint_signs[1] == -1,
                str(endpoint_signs),
            )
            inherited_lower, inherited_upper = map(
                parse_fraction, inherited_branch["q_interval"]
            )
            dyadic_lower = Fraction(numerator, 1 << bits)
            dyadic_upper = Fraction(numerator + 1, 1 << bits)
            check(
                f"branch {branch['root_index']} inherited bracket nesting",
                inherited_lower < dyadic_lower < dyadic_upper < inherited_upper
                and bool(branch["inside_inherited_interval"]),
            )
            denominator_signs = [
                polynomial_interval(data, numerator, bits).sign()
                for data in denominator_data
            ]
            check(
                f"branch {branch['root_index']} clearing denominators",
                all(denominator_signs)
                and denominator_signs
                == branch["rank_certificate"]["clearing_denominator_signs"],
                str(denominator_signs),
            )
            cofactor_intervals = {
                index: outward(
                    polynomial_interval(
                        artifact_poly_data(quotient["cofactors"][str(index)]),
                        numerator,
                        bits,
                    )
                )
                for index in (0, 1, 4, 6)
            }
            c00 = cofactor_intervals[0]
            check(
                f"branch {branch['root_index']} nonzero C00",
                c00.sign() == branch["rank_certificate"]["C00_sign"] != 0,
            )
            coefficient_intervals = {
                "a5": outward(cofactor_intervals[4] / c00),
                "a4": outward(cofactor_intervals[6] / c00),
                "a3": outward(cofactor_intervals[1] / c00),
            }
            a0 = rational_interval(target_data["r1_squared"], numerator, bits)
            c2 = (
                rational_interval(target_data["r2"], numerator, bits)
                - a0
                - 64
                - 32 * coefficient_intervals["a5"]
                - 16 * coefficient_intervals["a4"]
                - 8 * coefficient_intervals["a3"]
            )
            c3 = (
                rational_interval(target_data["r3_over_r1"], numerator, bits)
                - a0
                - 729
                - 243 * coefficient_intervals["a5"]
                - 81 * coefficient_intervals["a4"]
                - 27 * coefficient_intervals["a3"]
            )
            coefficient_intervals["a2"] = outward(c3 / 3 - c2 / 2)
            coefficient_intervals["a1"] = outward(3 * c2 / 2 - 2 * c3 / 3)
            coefficient_intervals["a0"] = a0
            ordered = [
                Interval(1),
                coefficient_intervals["a5"],
                coefficient_intervals["a4"],
                coefficient_intervals["a3"],
                coefficient_intervals["a2"],
                coefficient_intervals["a1"],
                coefficient_intervals["a0"],
            ]
            minor_signs, signature = hermite_signs(ordered)
            signatures.append(signature)
            check(
                f"branch {branch['root_index']} Hermite signature",
                minor_signs
                == branch["hermite_certificate"][
                    "leading_principal_minor_signs"
                ]
                and signature == branch["hermite_certificate"]["signature"],
                f"minor signs={minor_signs}; signature={signature}",
            )
            if int(branch["root_index"]) == 3:
                third_shifted = shifted_signs(ordered)
                check(
                    "branch 3 shifted Descartes signs",
                    third_shifted[0]
                    == branch["shifted_descartes_certificate"][
                        "coefficient_signs_descending"
                    ]
                    and third_shifted[1]
                    == branch["shifted_descartes_certificate"]["sign_variations"],
                    f"signs={third_shifted[0]}; variations={third_shifted[1]}",
                )
        check(
            "all branch dispositions",
            signatures == [2, 2, 6]
            and third_shifted is not None
            and third_shifted[1] == 2
            and all(
                branch["physical_branch_rejected"]
                for branch in artifact["data"]["branches"]
            ),
            f"Hermite signatures={signatures}",
        )
        disposition_cpu = time.process_time() - started
        limit = (
            QUOTIENT_CPU_LIMIT + DISPOSITION_CPU_LIMIT
            if stage == "all"
            else DISPOSITION_CPU_LIMIT
        )
        check(
            "disposition-stage CPU wall",
            disposition_cpu < limit,
            f"{disposition_cpu:.6f}s/{limit}s",
        )

    check(
        "RSS wall",
        peak_rss_bytes() < RSS_LIMIT_BYTES,
        f"{peak_rss_bytes()}/{RSS_LIMIT_BYTES} bytes",
    )
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print(f"PASS: trace branch disposition ({stage})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
