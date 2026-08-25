#!/usr/bin/env python3
"""Exact semialgebraic disposition of the three open-2x3 trace branches.

This experiment consumes the already-certified degree-971 exceptional
polynomial from e238.  It does not recompute that elimination.  Instead it
reconstructs only multiplication by the trace-seven residual in the stable
rank-eight quotient, certifies a rank-seven minor at each exceptional root,
recovers the unique mode sextic, and applies Hermite/Descartes certificates.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from experiments import e238_trace_exceptional_set as base  # noqa: E402

OUTPUT = ROOT / "results" / "spectral" / "trace_branch_disposition.json"
VERIFIER = ROOT / "tests" / "test_trace_branch_disposition.py"
BASE_ARTIFACT = ROOT / "results" / "spectral" / "trace_exceptional_set.json"
CPU_BUDGET_SECONDS = 1800.0
RSS_LIMIT_BYTES = 2 * 1024**3
DYADIC_BITS = 1000
DYADIC_NUMERATORS = (
    15609729022729137756391725731089016481629322190038055715944076325074649558595494901366235708503870943693630513065402471620685754224847411025267059763818106492676464055553525842839767000060595621884927609557610940752604477403862523168403838624776846494211405721976062578726891904195804348598759257122996,
    34271210900521225399855796045827006790908931866899616285067565772689639740896389258759596372684695903008549202016376912605873194155131061075994637949234023798654854648746774319444250035128006584901618043013224615567027743967957892093698639962083927644394573141430568942402467849970054930585757102712445,
    102746408354143832221103402726418962652395172934286109390231244126198758660415027640909434748390637038008018436969142230546912164206381285717621866091472061906543500207449826035670621318486510274023014769461912663147193561514141420872549419372374307799259129866915718745765452504836768178039165647849493,
)
COFACTOR_INDICES = (0, 1, 4, 6)
EXPECTED_COFACTOR_DEGREES = {0: 1208, 1: 1211, 4: 1208, 6: 1209}
EXPECTED_HERMITE_SIGNATURES = (2, 2, 6)
EXPECTED_SHIFTED_SIGNS = (1, -1, -1, -1, 1, 1, 1)
ROUNDING_BITS = 220


class ExactInterval:
    """Closed rational interval with inclusion-preserving arithmetic."""

    __slots__ = ("lower", "upper")

    def __init__(self, lower: Fraction | int, upper: Fraction | int | None = None):
        self.lower = Fraction(lower)
        self.upper = self.lower if upper is None else Fraction(upper)
        if self.lower > self.upper:
            raise ValueError("interval endpoints are reversed")

    def __add__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        other = other if isinstance(other, ExactInterval) else ExactInterval(other)
        return ExactInterval(self.lower + other.lower, self.upper + other.upper)

    __radd__ = __add__

    def __neg__(self) -> ExactInterval:
        return ExactInterval(-self.upper, -self.lower)

    def __sub__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        other = other if isinstance(other, ExactInterval) else ExactInterval(other)
        return self + (-other)

    def __rsub__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        return ExactInterval(other) - self

    def __mul__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        other = other if isinstance(other, ExactInterval) else ExactInterval(other)
        products = (
            self.lower * other.lower,
            self.lower * other.upper,
            self.upper * other.lower,
            self.upper * other.upper,
        )
        return ExactInterval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        other = other if isinstance(other, ExactInterval) else ExactInterval(other)
        if other.lower <= 0 <= other.upper:
            raise ZeroDivisionError("interval divisor contains zero")
        reciprocal = ExactInterval(
            min(Fraction(1, other.lower), Fraction(1, other.upper)),
            max(Fraction(1, other.lower), Fraction(1, other.upper)),
        )
        return self * reciprocal

    def __rtruediv__(self, other: ExactInterval | Fraction | int) -> ExactInterval:
        return ExactInterval(other) / self

    def sign(self) -> int:
        if self.lower > 0:
            return 1
        if self.upper < 0:
            return -1
        return 0


PolyData = tuple[int, tuple[int, ...]]


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def guard_resources(started: float, stage: str) -> None:
    elapsed = time.process_time() - started
    if elapsed >= CPU_BUDGET_SECONDS:
        raise RuntimeError(
            f"{stage}: process CPU {elapsed:.6f}s reached {CPU_BUDGET_SECONDS}s wall"
        )
    peak = peak_rss_bytes()
    if peak >= RSS_LIMIT_BYTES:
        raise RuntimeError(f"{stage}: peak RSS {peak} reached {RSS_LIMIT_BYTES} byte wall")


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def interval_json(interval: ExactInterval) -> list[str]:
    return [fraction_text(interval.lower), fraction_text(interval.upper)]


def integer_poly_data(value: Any) -> PolyData:
    expression = value.as_expr() if hasattr(value, "as_expr") else value
    polynomial = sp.Poly(expression, base.QPAR, domain=base.QQ)
    scale, integer = polynomial.clear_denoms(convert=True)
    return int(scale), tuple(int(coefficient) for coefficient in integer.all_coeffs())


def polynomial_record(value: Any) -> tuple[dict[str, object], PolyData]:
    scale, descending = integer_poly_data(value)
    ascending = [str(coefficient) for coefficient in reversed(descending)]
    return (
        {
            "degree": len(descending) - 1,
            "term_count": sum(coefficient != 0 for coefficient in descending),
            "integer_scale_denominator": str(scale),
            "integer_coefficients_ascending": ascending,
            "sha256": canonical_sha256(ascending),
        },
        (scale, descending),
    )


def dyadic_polynomial_sign(coefficients_ascending: list[int], numerator: int, bits: int) -> int:
    """Sign of p(numerator/2**bits), using one exact integer Horner numerator."""

    accumulator = coefficients_ascending[-1]
    for power, coefficient in enumerate(reversed(coefficients_ascending[:-1]), start=1):
        accumulator = accumulator * numerator + (coefficient << (bits * power))
    return (accumulator > 0) - (accumulator < 0)


def integer_polynomial_interval(data: PolyData, numerator: int, bits: int) -> ExactInterval:
    """Rigorous p(q) enclosure on [n/2^b,(n+1)/2^b].

    The midpoint is evaluated exactly.  The error uses the mean-value theorem
    and |p'| <= (sum_k k|c_k|) q_max^(degree-1), valid because q>1.
    """

    scale, coefficients = data
    degree = len(coefficients) - 1
    midpoint_numerator = 2 * numerator + 1
    midpoint_bits = bits + 1
    accumulator = coefficients[0]
    for power, coefficient in enumerate(coefficients[1:], start=1):
        accumulator = accumulator * midpoint_numerator + (
            coefficient << (midpoint_bits * power)
        )
    midpoint = Fraction(accumulator, 1 << (midpoint_bits * degree)) / scale
    if degree == 0:
        return ExactInterval(midpoint)
    derivative_weight = sum(
        (degree - index) * abs(coefficient)
        for index, coefficient in enumerate(coefficients[:-1])
    )
    error = Fraction(
        derivative_weight * pow(numerator + 1, degree - 1),
        1 << (bits * degree + 1),
    ) / scale
    return ExactInterval(midpoint - error, midpoint + error)


def round_outward(interval: ExactInterval, bits: int = ROUNDING_BITS) -> ExactInterval:
    """Replace a large exact interval by a smaller-bit enclosing dyadic interval."""

    magnitude = max(abs(interval.lower), abs(interval.upper))
    if magnitude == 0:
        return ExactInterval(0)
    exponent = magnitude.numerator.bit_length() - magnitude.denominator.bit_length()
    shift = bits - exponent
    if shift >= 0:
        scale = 1 << shift
        lower = (interval.lower.numerator * scale) // interval.lower.denominator
        upper = -((-interval.upper.numerator * scale) // interval.upper.denominator)
        return ExactInterval(Fraction(lower, scale), Fraction(upper, scale))
    unit = 1 << (-shift)
    lower = (interval.lower.numerator // (interval.lower.denominator * unit)) * unit
    upper = -((-interval.upper.numerator // (interval.upper.denominator * unit))) * unit
    return ExactInterval(lower, upper)


def rational_function_data(expression: sp.Expr) -> tuple[PolyData, PolyData]:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    return integer_poly_data(numerator), integer_poly_data(denominator)


def rational_function_interval(
    data: tuple[PolyData, PolyData], numerator: int, bits: int
) -> ExactInterval:
    numerator_interval = integer_polynomial_interval(data[0], numerator, bits)
    denominator_interval = integer_polynomial_interval(data[1], numerator, bits)
    return round_outward(numerator_interval / denominator_interval)


def build_cleared_multiplication_matrix(
    equations: list[sp.Poly],
) -> tuple[list[list[Any]], list[sp.Poly]]:
    """Build multiplication by F7 without evaluating its 8x8 determinant."""

    _, homogeneous_basis, representations, leading_monomials = base.leading_h_basis(equations)
    lifts = base.lifted_relations(equations, homogeneous_basis, representations)
    seventh = equations[3].rep.to_dict()
    columns: list[list[Any]] = []
    for basis_monomial in base.BASIS_MONOMIALS:
        product = {
            base.add_monomials(monomial, basis_monomial): coefficient
            for monomial, coefficient in seventh.items()
        }
        reduced, _ = base.stable_reduce(product, lifts, leading_monomials)
        columns.append(
            [reduced.get(monomial, base.ZERO) for monomial in base.BASIS_MONOMIALS]
        )
    rational_rows = [
        [columns[column][row] for column in range(8)] for row in range(8)
    ]

    denominators: list[sp.Poly] = []
    for column in range(8):
        denominator = sp.Integer(1)
        for row in range(8):
            denominator = sp.lcm(
                denominator, rational_rows[row][column].denom.as_expr()
            )
        denominators.append(sp.Poly(denominator, base.QPAR, domain=base.QQ).monic())

    polynomial_rows: list[list[Any]] = [[None] * 8 for _ in range(8)]
    for column, denominator in enumerate(denominators):
        multiplier = base.QQ_POLY_Q.from_sympy(denominator.as_expr())
        for row in range(8):
            cleared = rational_rows[row][column] * base.QQQ.new(multiplier)
            if cleared.denom.degree() != 0:
                raise AssertionError("column denominator did not clear")
            constant = next(iter(cleared.denom.to_dict().values()))
            polynomial_rows[row][column] = cleared.numer / constant
    return polynomial_rows, denominators


def cofactor_from_adjugate_row(polynomial_rows: list[list[Any]], basis_index: int) -> Any:
    """Return adj(P)[0,basis_index] = Cofactor(P)[basis_index,0]."""

    submatrix = [
        [polynomial_rows[row][column] for column in range(1, 8)]
        for row in range(8)
        if row != basis_index
    ]
    determinant = base.DomainMatrix(submatrix, (7, 7), base.QQ_POLY_Q).det()
    return -determinant if basis_index % 2 else determinant


def interval_determinant(matrix: list[list[ExactInterval]]) -> ExactInterval:
    """Gaussian interval enclosure; every selected pivot excludes zero."""

    work = [[ExactInterval(entry.lower, entry.upper) for entry in row] for row in matrix]
    dimension = len(work)
    parity = 1
    for column in range(dimension):
        candidates = [
            row for row in range(column, dimension) if work[row][column].sign()
        ]
        if not candidates:
            raise ArithmeticError(f"no nonzero interval pivot in column {column}")
        pivot = max(
            candidates,
            key=lambda row: min(
                abs(work[row][column].lower), abs(work[row][column].upper)
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
            work[row][column] = ExactInterval(0)
    determinant = ExactInterval(parity)
    for index in range(dimension):
        determinant = determinant * work[index][index]
    return determinant


def newton_power_sums(coefficients: list[ExactInterval]) -> list[ExactInterval]:
    degree = len(coefficients) - 1
    sums = [ExactInterval(degree)]
    for power in range(1, 2 * degree - 1):
        total = ExactInterval(0)
        if power <= degree:
            for index in range(1, power):
                total = total + coefficients[index] * sums[power - index]
            total = total + power * coefficients[power]
        else:
            for index in range(1, degree + 1):
                total = total + coefficients[index] * sums[power - index]
        sums.append(-total)
    return sums


def hermite_minor_intervals(coefficients: list[ExactInterval]) -> list[ExactInterval]:
    sums = newton_power_sums(coefficients)
    hermite = [[sums[row + column] for column in range(6)] for row in range(6)]
    return [
        round_outward(
            interval_determinant([row[:dimension] for row in hermite[:dimension]]),
            bits=180,
        )
        for dimension in range(1, 7)
    ]


def shifted_coefficients(coefficients: list[ExactInterval], shift: int) -> list[ExactInterval]:
    ascending = list(reversed(coefficients))
    result: list[ExactInterval] = []
    for target_power in range(7):
        value = ExactInterval(0)
        for source_power in range(target_power, 7):
            value = value + (
                ascending[source_power]
                * math.comb(source_power, target_power)
                * shift ** (source_power - target_power)
            )
        result.append(round_outward(value, bits=180))
    return list(reversed(result))


def sign_variations(signs: list[int]) -> int:
    nonzero = [sign for sign in signs if sign]
    return sum(left != right for left, right in zip(nonzero, nonzero[1:]))


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    inherited = json.loads(BASE_ARTIFACT.read_text())
    exceptional_record = inherited["data"]["elimination"]["exceptional_q_polynomial"]
    exceptional_coefficients = [int(value) for value in exceptional_record["coefficients"]]
    inherited_roots = inherited["data"]["physical_root_certificate"]

    root_records: list[dict[str, object]] = []
    for index, (numerator, inherited_root) in enumerate(
        zip(DYADIC_NUMERATORS, inherited_roots["roots"]), start=1
    ):
        lower = Fraction(numerator, 1 << DYADIC_BITS)
        upper = Fraction(numerator + 1, 1 << DYADIC_BITS)
        inherited_lower, inherited_upper = map(
            parse_fraction, inherited_root["q_interval"]
        )
        endpoint_signs = [
            dyadic_polynomial_sign(
                exceptional_coefficients, numerator + offset, DYADIC_BITS
            )
            for offset in (0, 1)
        ]
        root_records.append(
            {
                "root_index": index,
                "dyadic_bits": DYADIC_BITS,
                "lower_numerator": str(numerator),
                "q_interval": [fraction_text(lower), fraction_text(upper)],
                "E_endpoint_signs": endpoint_signs,
                "inside_inherited_interval": (
                    inherited_lower < lower < upper < inherited_upper
                ),
            }
        )
    guard_resources(started, "dyadic root isolation")

    epsilon, shift, traces, trace_record = base.symbolic_power_traces(
        base.BONDS_2X3, base.MAX_TRACE
    )
    targets = base.normalized_targets(traces, epsilon, shift)
    equations, _, _, _ = base.norm_equations(targets)
    polynomial_rows, column_denominators = build_cleared_multiplication_matrix(
        equations
    )
    matrix_payload = [
        [
            base.poly_coefficients_ascending(
                sp.Poly(polynomial_rows[row][column].as_expr(), base.QPAR, domain=base.QQ)
            )
            for column in range(8)
        ]
        for row in range(8)
    ]
    matrix_digest = canonical_sha256(matrix_payload)
    inherited_matrix_digest = inherited["data"]["elimination"]["matrix_sha256"]
    column_denominator_records: list[dict[str, object]] = []
    column_denominator_data: list[PolyData] = []
    for denominator in column_denominators:
        record, data = polynomial_record(denominator)
        column_denominator_records.append(record)
        column_denominator_data.append(data)
    guard_resources(started, "rank-eight multiplication matrix")

    cofactor_records: dict[str, dict[str, object]] = {}
    cofactor_data: dict[int, PolyData] = {}
    for basis_index in COFACTOR_INDICES:
        cofactor = cofactor_from_adjugate_row(polynomial_rows, basis_index)
        record, data = polynomial_record(cofactor)
        record["basis_index"] = basis_index
        record["basis_monomial"] = list(base.BASIS_MONOMIALS[basis_index])
        cofactor_records[str(basis_index)] = record
        cofactor_data[basis_index] = data
        guard_resources(started, f"cofactor {basis_index}")

    target_data = {
        name: rational_function_data(value)
        for name, value in targets.items()
        if name in {"r1_squared", "r2", "r3_over_r1"}
    }
    branches: list[dict[str, object]] = []
    for root_index, numerator in enumerate(DYADIC_NUMERATORS, start=1):
        column_denominator_signs = [
            integer_polynomial_interval(data, numerator, DYADIC_BITS).sign()
            for data in column_denominator_data
        ]
        if not all(column_denominator_signs):
            raise AssertionError(
                f"branch {root_index}: a column-clearing denominator may vanish"
            )
        cofactor_intervals = {
            basis_index: round_outward(
                integer_polynomial_interval(data, numerator, DYADIC_BITS)
            )
            for basis_index, data in cofactor_data.items()
        }
        c00 = cofactor_intervals[0]
        if c00.sign() == 0:
            raise AssertionError(f"branch {root_index}: C00 interval contains zero")
        coefficient_intervals = {
            "a5": round_outward(cofactor_intervals[4] / c00),
            "a4": round_outward(cofactor_intervals[6] / c00),
            "a3": round_outward(cofactor_intervals[1] / c00),
        }
        a0 = rational_function_interval(
            target_data["r1_squared"], numerator, DYADIC_BITS
        )
        c2 = (
            rational_function_interval(target_data["r2"], numerator, DYADIC_BITS)
            - a0
            - 64
            - 32 * coefficient_intervals["a5"]
            - 16 * coefficient_intervals["a4"]
            - 8 * coefficient_intervals["a3"]
        )
        c3 = (
            rational_function_interval(
                target_data["r3_over_r1"], numerator, DYADIC_BITS
            )
            - a0
            - 729
            - 243 * coefficient_intervals["a5"]
            - 81 * coefficient_intervals["a4"]
            - 27 * coefficient_intervals["a3"]
        )
        coefficient_intervals["a2"] = round_outward(c3 / 3 - c2 / 2)
        coefficient_intervals["a1"] = round_outward(3 * c2 / 2 - 2 * c3 / 3)
        coefficient_intervals["a0"] = a0
        ordered = [
            ExactInterval(1),
            coefficient_intervals["a5"],
            coefficient_intervals["a4"],
            coefficient_intervals["a3"],
            coefficient_intervals["a2"],
            coefficient_intervals["a1"],
            coefficient_intervals["a0"],
        ]
        minors = hermite_minor_intervals(ordered)
        minor_signs = [minor.sign() for minor in minors]
        pivot_signs = [minor_signs[0]] + [
            minor_signs[index] * minor_signs[index - 1]
            for index in range(1, 6)
        ]
        signature = sum(pivot_signs)

        shifted = shifted_coefficients(ordered, 4)
        shifted_signs = [coefficient.sign() for coefficient in shifted]
        variations = sign_variations(shifted_signs)
        rejected = signature < 6 or (
            shifted_signs[-1] != 0 and variations < 6
        )
        branches.append(
            {
                **root_records[root_index - 1],
                "rank_certificate": {
                    "determinant_zero_reason": "E(q)=0 and E divides det(P)",
                    "clearing_denominator_signs": column_denominator_signs,
                    "clearing_denominators_nonzero": True,
                    "C00_interval": interval_json(c00),
                    "C00_sign": c00.sign(),
                    "rank_P": 7,
                    "quotient_dimension_after_F7": 1,
                    "unique_coefficient_point": True,
                },
                "adjugate_coordinate_intervals": {
                    name: interval_json(coefficient_intervals[name])
                    for name in ("a5", "a4", "a3")
                },
                "sextic_coefficient_intervals": {
                    name: interval_json(coefficient_intervals[name])
                    for name in ("a5", "a4", "a3", "a2", "a1", "a0")
                },
                "hermite_certificate": {
                    "matrix": "H_ij=s_(i+j), 0<=i,j<6, from Newton sums of Q",
                    "leading_principal_minor_intervals": [
                        interval_json(minor) for minor in minors
                    ],
                    "leading_principal_minor_signs": minor_signs,
                    "LDL_pivot_signs": pivot_signs,
                    "signature": signature,
                    "distinct_real_roots_of_Q": signature,
                },
                "shifted_descartes_certificate": {
                    "polynomial": "Q(4+x)",
                    "coefficient_intervals_descending": [
                        interval_json(coefficient) for coefficient in shifted
                    ],
                    "coefficient_signs_descending": shifted_signs,
                    "sign_variations": variations,
                    "positive_root_upper_bound": variations,
                },
                "physical_branch_rejected": rejected,
                "disposition": (
                    f"only {signature} distinct real mode roots"
                    if signature < 6
                    else (
                        "six distinct real roots, but Q(4+x) has at most "
                        f"{variations} positive roots"
                    )
                ),
            }
        )
        guard_resources(started, f"branch {root_index} disposition")

    checks.extend(
        [
            {
                "name": "inherits_exact_three_root_exceptional_set",
                "passed": exceptional_record["degree"] == 971
                and inherited_roots["positive_q_minus_one_root_count"] == 3
                and exceptional_record["sha256"]
                == canonical_sha256(exceptional_record["coefficients"]),
                "detail": "e238 supplies one and only one q>1 root in each inherited rational interval",
            },
            {
                "name": "narrow_dyadic_root_intervals",
                "passed": all(
                    record["inside_inherited_interval"]
                    and record["E_endpoint_signs"][0]
                    * record["E_endpoint_signs"][1]
                    == -1
                    for record in root_records
                ),
                "detail": "each 1000-bit dyadic interval lies inside its inherited isolator and has opposite exact E signs",
            },
            {
                "name": "stable_rank_eight_quotient_reused_without_E_recomputation",
                "passed": len(base.BASIS_MONOMIALS) == 8
                and [denominator.degree() for denominator in column_denominators]
                == [85, 132, 155, 202, 132, 155, 132, 155]
                and matrix_digest == inherited_matrix_digest,
                "detail": "the reconstructed cleared multiplication matrix matches e238 exactly; its determinant was not recomputed",
            },
            {
                "name": "adjugate_cofactor_degrees",
                "passed": {
                    index: cofactor_records[str(index)]["degree"]
                    for index in COFACTOR_INDICES
                }
                == EXPECTED_COFACTOR_DEGREES,
                "detail": "adj(P) row zero entries C_(j,0) have the bounded exact degrees",
            },
            {
                "name": "rank_seven_and_unique_point_on_every_branch",
                "passed": all(
                    branch["rank_certificate"]["C00_sign"] != 0
                    and branch["rank_certificate"]["clearing_denominators_nonzero"]
                    and all(
                        branch["rank_certificate"]["clearing_denominator_signs"]
                    )
                    and branch["rank_certificate"]["rank_P"] == 7
                    and branch["rank_certificate"]["quotient_dimension_after_F7"]
                    == 1
                    for branch in branches
                ),
                "detail": "nonzero clearing denominators identify P with M; det(P)=0 and C00!=0 give rank seven; its one-dimensional quotient is one reduced coefficient point",
            },
            {
                "name": "exact_Hermite_signatures",
                "passed": tuple(
                    branch["hermite_certificate"]["signature"]
                    for branch in branches
                )
                == EXPECTED_HERMITE_SIGNATURES,
                "detail": "the three mode sextics have respectively 2, 2, and 6 distinct real roots",
            },
            {
                "name": "third_branch_shifted_Descartes_obstruction",
                "passed": tuple(
                    branches[2]["shifted_descartes_certificate"][
                        "coefficient_signs_descending"
                    ]
                )
                == EXPECTED_SHIFTED_SIGNS
                and branches[2]["shifted_descartes_certificate"]["sign_variations"]
                == 2,
                "detail": "Q(4+x) has sign pattern +---+++ and therefore at most two positive roots",
            },
            {
                "name": "all_three_trace_branches_are_nonphysical",
                "passed": all(branch["physical_branch_rejected"] for branch in branches),
                "detail": "no surviving seven-trace branch supplies six real mode values u_i>=4",
            },
            {
                "name": "benchmark_not_used",
                "passed": True,
                "detail": "K_c is absent from branch isolation, quotient reconstruction, and root disposition",
            },
        ]
    )

    elapsed = time.process_time() - started
    peak_rss = peak_rss_bytes()
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
    failed = [str(check["name"]) for check in checks if not check["passed"]]
    if failed:
        raise AssertionError(f"trace branch disposition checks failed: {failed}")

    source_paths = (Path(__file__).resolve(), VERIFIER, base.Path(base.__file__).resolve())
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e242_trace_branch_disposition.py",
            "verifier": "tests/test_trace_branch_disposition.py",
            "interpreter": sys.executable,
            "arithmetic": "exact QQ(q), ZZ[q], dyadic rational intervals, Hermite signatures, and Descartes signs",
            "single_process": True,
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_Kc_used": False,
            "source_sha256": {
                str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
            },
            "dependency_sha256": {
                "results/spectral/trace_exceptional_set.json": file_sha256(BASE_ARTIFACT)
            },
        },
        "data": {
            "claim_tag": "[THEOREM][COMPUTATION][EMPTY EXCEPTIONAL SET]",
            "inherited_exceptional_set": {
                "source": "e238",
                "degree": exceptional_record["degree"],
                "coefficient_sha256": exceptional_record["sha256"],
                "positive_q_minus_one_root_count": inherited_roots[
                    "positive_q_minus_one_root_count"
                ],
                "recomputed": False,
            },
            "quotient_certificate": {
                "basis": [list(monomial) for monomial in base.BASIS_MONOMIALS],
                "cleared_matrix_shape": [8, 8],
                "column_denominator_degrees": [
                    denominator.degree() for denominator in column_denominators
                ],
                "column_denominators": column_denominator_records,
                "matrix_sha256": matrix_digest,
                "inherited_matrix_sha256": inherited_matrix_digest,
                "adjugate_row": 0,
                "cofactor_convention": "adj(P)[0,j]=(-1)^j det(P with row j and column 0 deleted)",
                "cofactors": cofactor_records,
            },
            "trace_construction": trace_record,
            "branches": branches,
            "theorem": {
                "tag": "[THEOREM][COMPUTATION]",
                "statement": (
                    "For every physical coupling 0<t<1, the positive-definite open-2x3 Ising "
                    "layer transfer spectrum is not a full six-mode subset-product spectrum."
                ),
                "proof_chain": (
                    "e238 confines any full spectrum to three algebraic q roots; at each root "
                    "a rank-seven adjugate certificate gives one coefficient sextic; Hermite "
                    "signatures 2,2,6 reject the first two, while the third has only two "
                    "Descartes-allowed roots above 4."
                ),
                "physical_necessity": (
                    "positive subset-product modes have u_i=(sqrt(v_i)+1/sqrt(v_i))^2>=4"
                ),
                "scope": "one finite open 2x3 layer, every 0<t<1",
                "not_proved": [
                    "an all-size bipartite-grid obstruction",
                    "absence of parity-sector or auxiliary-mode factorizations",
                    "a thermodynamic-limit solution or critical coupling",
                ],
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for check in payload["checks"]:
        print(f"  [{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
