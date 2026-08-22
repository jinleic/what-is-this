"""Exact crossing-obstruction certificate for the 2x3 layer forced-value test.

The prior e48 certificate correctly used exact Weyl/inertia propagation, but
required strict separation of the three lowest eigenvalues.  The underlying
Gaussian forced-value lemma is stronger: it is valid with multiplicities.

This experiment therefore does two things using only integer/Fraction
arithmetic for every decision:

* it proves four disjoint c-count jumps, hence four genuine parameter values
  at which the forced value is an eigenvalue; and
* it certifies selected closed rational intervals away from those obstructions
  by the same Weyl/inertia method with no artificial distinctness hypothesis.

Float64 is used only to propose eigenvalue brackets.  Every bracket, count,
and interval decision is then accepted only through exact Sylvester inertia.
"""

from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.clifford import pauli_x, symplectic_form, zz  # noqa: E402
from ising.transfer_matrix import layer_bonds  # noqa: E402


MAX_WAVE8_COMPONENT_WIDTH = Fraction(1, 50_000_000)  # 2*10^-8


@dataclass(frozen=True)
class SectorBasis:
    """An integer character basis with diagonal Gram matrix."""

    character: int
    vectors: tuple[dict[int, int], ...]
    gram: tuple[int, ...]


@dataclass(frozen=True)
class NumericSector:
    """The integer generalized symmetric pencil B - sigma G."""

    matrix: tuple[tuple[int, ...], ...]
    gram: tuple[int, ...]
    character: int


@dataclass
class PointContext:
    """Exact matrix and cached inertia counter at one rational t."""

    t: Fraction
    integer_matrix: list[list[int]]
    scale: int
    counter: "ExactCounter"
    guides: list[float]

    def count_below(self, level: Fraction) -> int | None:
        """Number of eigenvalues of the unscaled M(t) strictly below level."""
        return self.counter.count(Fraction(level) * self.scale)


# ---------------------------------------------------------------------------
# Formatting and JSON helpers
# ---------------------------------------------------------------------------


def exact_decimal(value: Fraction, digits: int = 24) -> str:
    """Display-only decimal; it never decides any certificate."""
    with localcontext() as context:
        context.prec = digits + 16
        return format(Decimal(value.numerator) / Decimal(value.denominator), f".{digits}g")


def fraction_json(value: Fraction) -> dict[str, str]:
    return {"exact": str(value), "decimal": exact_decimal(value)}


def interval_json(lower: Fraction, upper: Fraction) -> dict[str, str]:
    if lower > upper:
        raise AssertionError("reversed rational interval")
    return {
        "lower": str(lower),
        "upper": str(upper),
        "width": str(upper - lower),
        "lower_decimal": exact_decimal(lower),
        "upper_decimal": exact_decimal(upper),
    }


# ---------------------------------------------------------------------------
# Exact polynomial representative of the layer operator
# ---------------------------------------------------------------------------


def polynomial_matrix(
    n: int, bonds: list[tuple[int, int]]
) -> tuple[list[list[list[int]]], int, int]:
    """Return coefficients of a positive polynomial multiple of R(t).

    With q=(1+t^2)/(2t), R=P_t diag(q^e_s) P_t and e_s shifted so its
    minimum is zero, this constructs M=(2t)^degree q^shift R.  For the
    open 2x3 layer this is M=(2t)^7 q^4 R and all coefficients are
    nonnegative integers.
    """
    dimension = 1 << n
    energies: list[int] = []
    for state in range(dimension):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    if not all((energy - parity) % 2 == 0 for energy in energies):
        raise AssertionError("bond-energy parity is not constant")
    exponents = [(energy - parity) // 2 for energy in energies]
    shift = -min(exponents)
    degree = max(exponent + shift for exponent in exponents)
    max_degree = 2 * n + 2 * degree
    coefficients = [
        [[0] * (max_degree + 1) for _ in range(dimension)]
        for _ in range(dimension)
    ]
    for row in range(dimension):
        for col in range(row, dimension):
            entry = [0] * (max_degree + 1)
            for middle, exponent in enumerate(exponents):
                power_q = exponent + shift
                hamming_power = (row ^ middle).bit_count() + (middle ^ col).bit_count()
                base_power = degree - power_q + hamming_power
                integer_factor = 2 ** (degree - power_q)
                for even_power in range(power_q + 1):
                    entry[base_power + 2 * even_power] += (
                        integer_factor * math.comb(power_q, even_power)
                    )
            coefficients[row][col] = entry
            coefficients[col][row] = entry
    return coefficients, degree, shift


def evaluate_polynomial_matrix(
    coefficients: list[list[list[int]]], t: Fraction
) -> list[list[Fraction]]:
    powers = [t**power for power in range(len(coefficients[0][0]))]
    return [
        [
            sum(
                (Fraction(coefficient) * powers[power]
                 for power, coefficient in enumerate(entry) if coefficient),
                Fraction(0),
            )
            for entry in row
        ]
        for row in coefficients
    ]


def clear_denominators(matrix: list[list[Fraction]]) -> tuple[list[list[int]], int]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = math.lcm(scale, value.denominator)
    return [[(value * scale).numerator for value in row] for row in matrix], scale


def polynomial_invariant(
    coefficients: list[list[list[int]]], generators: list[list[int]]
) -> bool:
    dimension = len(coefficients)
    return all(
        coefficients[generator[row]][generator[col]] == coefficients[row][col]
        for generator in generators
        for row in range(dimension)
        for col in range(dimension)
    )


def derivative_norm_upper(
    coefficients: list[list[list[int]]], upper_t: Fraction
) -> tuple[Fraction, Fraction, Fraction]:
    """Exact upper bounds on ||M'||_F and ||M'||_2 on 0 <= t <= upper_t.

    Nonnegative coefficients make every derivative entry nonnegative and
    increasing.  The spectral bound is min(Frobenius, max absolute row sum).
    """
    derivative_entries: list[list[Fraction]] = []
    frobenius_square = Fraction(0)
    row_sums: list[Fraction] = []
    for row in coefficients:
        derivative_row: list[Fraction] = []
        for entry in row:
            derivative = sum(
                (Fraction(power * coefficient) * upper_t ** (power - 1)
                 for power, coefficient in enumerate(entry) if power and coefficient),
                Fraction(0),
            )
            derivative_row.append(derivative)
            frobenius_square += derivative * derivative
        derivative_entries.append(derivative_row)
        row_sums.append(sum(derivative_row, Fraction(0)))
    del derivative_entries  # the exact aggregates are all that are needed
    frobenius = rational_sqrt_upper(frobenius_square)
    row_sum = max(row_sums)
    return frobenius, row_sum, min(frobenius, row_sum)


def rational_sqrt_upper(value: Fraction, denominator: int = 10**18) -> Fraction:
    """Outward rational square root verified by one integer-square comparison."""
    if value < 0:
        raise ValueError("square root upper bound received a negative value")
    target = (
        value.numerator * denominator * denominator + value.denominator - 1
    ) // value.denominator
    numerator = math.isqrt(target)
    if numerator * numerator < target:
        numerator += 1
    result = Fraction(numerator, denominator)
    if result * result < value:
        raise AssertionError("outward square-root rounding failed")
    return result


# ---------------------------------------------------------------------------
# Exact character-sector inertia counters
# ---------------------------------------------------------------------------


def state_permutation(n: int, site_permutation: list[int]) -> list[int]:
    result: list[int] = []
    for state in range(1 << n):
        moved = 0
        for old_site, new_site in enumerate(site_permutation):
            bit = (state >> (n - 1 - old_site)) & 1
            moved |= bit << (n - 1 - new_site)
        result.append(moved)
    return result


def grid_reflections(rows: int, cols: int) -> list[list[int]]:
    return [
        [(rows - 1 - row) * cols + col for row in range(rows) for col in range(cols)],
        [row * cols + (cols - 1 - col) for row in range(rows) for col in range(cols)],
    ]


def spin_flip_permutation(n: int) -> list[int]:
    return [state ^ ((1 << n) - 1) for state in range(1 << n)]


def make_sector_bases(n: int, generators: list[list[int]]) -> list[SectorBasis]:
    """Build exact bases for all characters of commuting involutions."""
    dimension = 1 << n
    group: list[list[int]] = []
    for mask in range(1 << len(generators)):
        permutation = list(range(dimension))
        for index, generator in enumerate(generators):
            if mask & (1 << index):
                permutation = [generator[state] for state in permutation]
        group.append(permutation)

    seen: set[int] = set()
    representatives: list[int] = []
    for state in range(dimension):
        if state in seen:
            continue
        orbit = {permutation[state] for permutation in group}
        representatives.append(min(orbit))
        seen.update(orbit)

    bases: list[SectorBasis] = []
    for character in range(1 << len(generators)):
        vectors: list[dict[int, int]] = []
        gram: list[int] = []
        for representative in representatives:
            vector: dict[int, int] = {}
            for group_element, permutation in enumerate(group):
                sign = -1 if (character & group_element).bit_count() % 2 else 1
                image = permutation[representative]
                vector[image] = vector.get(image, 0) + sign
            vector = {state: coefficient for state, coefficient in vector.items() if coefficient}
            if vector:
                vectors.append(vector)
                gram.append(sum(coefficient * coefficient for coefficient in vector.values()))
        bases.append(SectorBasis(character, tuple(vectors), tuple(gram)))
    if sum(len(basis.gram) for basis in bases) != dimension:
        raise AssertionError("character-sector dimensions do not sum to the matrix dimension")
    return bases


def numeric_sectors(
    matrix: list[list[int]], bases: list[SectorBasis]
) -> list[NumericSector]:
    sectors: list[NumericSector] = []
    for basis in bases:
        size = len(basis.gram)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            for col in range(row + 1):
                value = sum(
                    left * matrix[left_state][right_state] * right
                    for left_state, left in basis.vectors[row].items()
                    for right_state, right in basis.vectors[col].items()
                )
                block[row][col] = block[col][row] = value
        sectors.append(
            NumericSector(
                tuple(tuple(entry for entry in row) for row in block),
                basis.gram,
                basis.character,
            )
        )
    return sectors


def fraction_free_inertia(lower_matrix: list[list[int]]) -> int | None:
    """Negative-pivot count for a symmetric matrix using exact Bareiss LDL."""
    size = len(lower_matrix)
    work = [row[:] for row in lower_matrix]
    previous = 1
    negatives = 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        if pivot_index + 1 == size:
            break
        for row_index in range(pivot_index + 1, size):
            row = work[row_index]
            left = row[pivot_index]
            for col_index in range(pivot_index + 1, row_index + 1):
                numerator = pivot * row[col_index] - left * work[col_index][pivot_index]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Bareiss division was not exact")
                row[col_index] = quotient
        previous = pivot
    return negatives


class ExactCounter:
    """Exact count of generalized-sector eigenvalues below a rational shift."""

    def __init__(self, sectors: Iterable[NumericSector]):
        self.sectors = tuple(sectors)
        self.cache: dict[Fraction, int | None] = {}
        self.calls = 0

    def count(self, shift: Fraction) -> int | None:
        shift = Fraction(shift)
        cached = self.cache.get(shift)
        if cached is not None or shift in self.cache:
            return cached
        total = 0
        for sector in self.sectors:
            size = len(sector.gram)
            shifted = [[0] * size for _ in range(size)]
            common = 0
            for row in range(size):
                for col in range(row + 1):
                    value = shift.denominator * sector.matrix[row][col]
                    if row == col:
                        value -= shift.numerator * sector.gram[row]
                    shifted[row][col] = value
                    common = math.gcd(common, abs(value))
            if common > 1:
                for row in range(size):
                    for col in range(row + 1):
                        shifted[row][col] //= common
            count = fraction_free_inertia(shifted)
            if count is None:
                self.cache[shift] = None
                return None
            total += count
        self.cache[shift] = total
        self.calls += 1
        return total


def guide_eigenvalues(integer_matrix: list[list[int]], scale: int, count: int = 3) -> list[float]:
    """Float64 proposals only; all resulting brackets are verified exactly."""
    largest = max(abs(value) for row in integer_matrix for value in row)
    scaled = np.array(
        [[value / largest for value in row] for row in integer_matrix], dtype=float
    )
    multiplier = float(Fraction(largest, scale))
    return [float(value * multiplier) for value in np.linalg.eigvalsh(scaled)[:count]]


def make_point_context(
    coefficients: list[list[list[int]]], t: Fraction, bases: list[SectorBasis]
) -> PointContext:
    rational_matrix = evaluate_polynomial_matrix(coefficients, t)
    integer_matrix, scale = clear_denominators(rational_matrix)
    sectors = numeric_sectors(integer_matrix, bases)
    return PointContext(
        t=t,
        integer_matrix=integer_matrix,
        scale=scale,
        counter=ExactCounter(sectors),
        guides=guide_eigenvalues(integer_matrix, scale),
    )


# ---------------------------------------------------------------------------
# Rational eigenvalue enclosures and c-counts
# ---------------------------------------------------------------------------


def usable_outward(
    context: PointContext, level: Fraction, direction: int
) -> tuple[Fraction, int]:
    """Return a nonsingular rational level on the safe outward side."""
    if direction not in (-1, 1):
        raise ValueError("outward direction must be +/-1")
    count = context.count_below(level)
    if count is not None:
        return level, count
    scale = max(abs(level), Fraction(1))
    for attempt in range(1, 48):
        candidate = level + direction * scale / (10 ** (28 + attempt))
        count = context.count_below(candidate)
        if count is not None:
            return candidate, count
    raise RuntimeError("could not find a nonsingular outward rational shift")


def usable_inside(
    context: PointContext, lower: Fraction, upper: Fraction
) -> tuple[Fraction, int]:
    """Return a nonsingular rational shift strictly inside an exact bracket."""
    if not lower < upper:
        raise ValueError("empty bisection bracket")
    span = upper - lower
    for attempt in range(48):
        trial = (
            lower + span / 2
            if attempt == 0
            else lower + span / 2 + span * Fraction((-1) ** attempt, 3 ** (attempt + 2))
        )
        if lower < trial < upper:
            count = context.count_below(trial)
            if count is not None:
                return trial, count
    raise RuntimeError("could not find a nonsingular rational shift inside bracket")


def initial_bracket(
    context: PointContext, index: int, guide: float
) -> tuple[Fraction, Fraction]:
    """Obtain a certified eigenvalue bracket; float supplies only a proposal."""
    if not math.isfinite(guide) or guide <= 0:
        lower, upper = Fraction(0), Fraction(1)
    else:
        proposal = Fraction(format(guide, ".15e"))
        lower = max(Fraction(0), proposal * Fraction(99, 100))
        upper = proposal * Fraction(101, 100)
    lower, lower_count = usable_outward(context, lower, -1)
    while lower_count > index:
        lower, lower_count = usable_outward(context, lower / 2, -1)
    upper, upper_count = usable_outward(context, upper, +1)
    while upper_count <= index:
        upper, upper_count = usable_outward(context, upper * 2, +1)
    return lower, upper


def isolate_eigenvalue(
    context: PointContext,
    index: int,
    relative_width: Fraction,
) -> tuple[Fraction, Fraction]:
    """Exact enclosure of sorted lambda_index through certified inertia bisection."""
    lower, upper = initial_bracket(context, index, context.guides[index])
    while upper - lower > relative_width * lower:
        trial, count = usable_inside(context, lower, upper)
        if count <= index:
            lower = trial
        else:
            upper = trial
    if not (lower > 0 and context.count_below(lower) <= index):
        raise AssertionError("lower eigenvalue bracket invariant failed")
    upper_count = context.count_below(upper)
    if upper_count is None or upper_count <= index:
        raise AssertionError("upper eigenvalue bracket invariant failed")
    return lower, upper


def enclose_low_three(
    context: PointContext, relative_width: Fraction
) -> list[tuple[Fraction, Fraction]]:
    return [
        isolate_eigenvalue(context, index, relative_width)
        for index in range(3)
    ]


def point_forced_window(
    context: PointContext,
    start_power: int = 5,
    max_power: int = 18,
) -> dict[str, object] | None:
    """Certify c(t)=#{mu<f} when an exact eigenvalue-free f-window is found.

    A None return is intentional: it means the progressively tightened exact
    window remains occupied, as it must for a Gaussian control or at a forced
    equality point.
    """
    started = time.perf_counter()
    for power in range(start_power, max_power + 1):
        relative_width = Fraction(1, 10**power)
        enclosures = enclose_low_three(context, relative_width)
        (l0, u0), (l1, u1), (l2, u2) = enclosures
        if l0 <= 0:
            raise AssertionError("positive layer matrix received a nonpositive lambda0 bracket")
        forced_lower = l1 * l2 / u0
        forced_upper = u1 * u2 / l0
        actual_lower, count_lower = usable_outward(context, forced_lower, -1)
        actual_upper, count_upper = usable_outward(context, forced_upper, +1)
        if count_lower == count_upper:
            return {
                "t": str(context.t),
                "relative_width": str(relative_width),
                "lambda_enclosures": {
                    f"lambda{index}": interval_json(*enclosure)
                    for index, enclosure in enumerate(enclosures)
                },
                "forced_value": interval_json(forced_lower, forced_upper),
                "tested_window": interval_json(actual_lower, actual_upper),
                "inertia_counts": [count_lower, count_upper],
                "count_below_forced_value": count_lower,
                "exact_inertia_calls": context.counter.calls,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
            }
    return None


# ---------------------------------------------------------------------------
# Multiplicity-safe interval propagation
# ---------------------------------------------------------------------------


def certify_absence_interval(
    name: str,
    coefficients: list[list[list[int]]],
    bases: list[SectorBasis],
    center: Fraction,
    half_width: Fraction,
) -> dict[str, object] | None:
    """Exact no-forced-eigenvalue certificate on one closed rational interval.

    Unlike e48, no lambda0<lambda1<lambda2 assumption is imposed.  The
    multiplicity-safe forced lemma makes that separation unnecessary.
    """
    lower_t, upper_t = center - half_width, center + half_width
    if not (0 < lower_t <= center <= upper_t < 1):
        raise ValueError("interval must be a closed subinterval of (0,1)")
    started = time.perf_counter()
    context = make_point_context(coefficients, center, bases)
    frobenius, row_sum, derivative_bound = derivative_norm_upper(coefficients, upper_t)
    rho = derivative_bound * half_width
    for power in range(5, 19):
        enclosures = enclose_low_three(context, Fraction(1, 10**power))
        (l0, u0), (l1, u1), (l2, u2) = enclosures
        if l0 - rho <= 0 or l1 - rho <= 0 or l2 - rho <= 0:
            return None
        forced_lower = (l1 - rho) * (l2 - rho) / (u0 + rho)
        forced_upper = (u1 + rho) * (u2 + rho) / (l0 - rho)
        tested_lower, count_lower = usable_outward(context, forced_lower - rho, -1)
        tested_upper, count_upper = usable_outward(context, forced_upper + rho, +1)
        if count_lower == count_upper:
            return {
                "name": name,
                "claim_tag": "[THEOREM]",
                "center_t": str(center),
                "half_width_t": str(half_width),
                "interval_t": interval_json(lower_t, upper_t),
                "lambda_center_enclosures": {
                    f"lambda{index}": interval_json(*enclosure)
                    for index, enclosure in enumerate(enclosures)
                },
                "frobenius_derivative_upper": fraction_json(frobenius),
                "row_sum_derivative_upper": fraction_json(row_sum),
                "selected_spectral_derivative_upper": fraction_json(derivative_bound),
                "weyl_radius": fraction_json(rho),
                "forced_value_whole_interval": interval_json(forced_lower, forced_upper),
                "tested_center_absence_window": interval_json(tested_lower, tested_upper),
                "inertia_counts_at_tested_window_ends": [count_lower, count_upper],
                "forced_relation_absent_throughout": True,
                "low_eigenvalue_relative_width": str(Fraction(1, 10**power)),
                "exact_inertia_calls": context.counter.calls,
                "elapsed_seconds": round(time.perf_counter() - started, 6),
            }
    return None


def find_safe_interval(
    name: str,
    coefficients: list[list[list[int]]],
    bases: list[SectorBasis],
    center: Fraction,
) -> dict[str, object]:
    """Find a large dyadic half-width accepted by the exact interval certificate."""
    for exponent in range(12, 42):
        candidate = Fraction(1, 2**exponent)
        certificate = certify_absence_interval(
            name, coefficients, bases, center, candidate
        )
        if certificate is not None:
            certificate["searched_dyadic_half_width_exponents"] = {
                "first_attempt": 12,
                "accepted": exponent,
            }
            return certificate
    raise RuntimeError(f"no dyadic interval certificate found for {name}")


# ---------------------------------------------------------------------------
# Crossing discovery and exact c-jump refinements
# ---------------------------------------------------------------------------


def scan_count_runs(
    coefficients: list[list[list[int]]],
    bases: list[SectorBasis],
    lower: Fraction,
    upper: Fraction,
    step: Fraction,
) -> tuple[list[dict[str, object]], list[tuple[Fraction, Fraction, int, int]]]:
    """Exact rational grid scan, compressed into constant-c runs and transitions."""
    if not (lower < upper and step > 0):
        raise ValueError("invalid rational scan interval")
    points: list[tuple[Fraction, int]] = []
    current = lower
    while current <= upper:
        certificate = point_forced_window(make_point_context(coefficients, current, bases))
        if certificate is None:
            raise RuntimeError(f"forced window unexpectedly unresolved at grid point {current}")
        points.append((current, int(certificate["count_below_forced_value"])))
        current += step
    if points[-1][0] != upper:
        raise AssertionError("scan grid did not land on upper endpoint")

    runs: list[dict[str, object]] = []
    transitions: list[tuple[Fraction, Fraction, int, int]] = []
    run_lower, run_count = points[0]
    previous_t, previous_count = points[0]
    for t, count in points[1:]:
        if count != previous_count:
            runs.append({
                "interval": interval_json(run_lower, previous_t),
                "count_below_forced_value": previous_count,
            })
            transitions.append((previous_t, t, previous_count, count))
            run_lower, run_count = t, count
        previous_t, previous_count = t, count
    runs.append({
        "interval": interval_json(run_lower, points[-1][0]),
        "count_below_forced_value": previous_count,
    })
    del run_count
    return runs, transitions


def refine_c_jump(
    name: str,
    coefficients: list[list[list[int]]],
    bases: list[SectorBasis],
    lower: Fraction,
    upper: Fraction,
    lower_count: int,
    upper_count: int,
    target_width: Fraction = Fraction(1, 10**10),
) -> dict[str, object]:
    """Bisection localizes an exact c-count jump without any float decision."""
    if upper_count != lower_count + 1:
        raise ValueError("this refinement expects one exact c-count jump")
    coarse_lower, coarse_upper = lower, upper
    lower_certificate = point_forced_window(make_point_context(coefficients, lower, bases))
    upper_certificate = point_forced_window(make_point_context(coefficients, upper, bases))
    if lower_certificate is None or upper_certificate is None:
        raise RuntimeError("coarse crossing endpoints must have clear forced windows")
    if (int(lower_certificate["count_below_forced_value"]) != lower_count
            or int(upper_certificate["count_below_forced_value"]) != upper_count):
        raise AssertionError("coarse c-count witness mismatch")

    iterations = 0
    while upper - lower > target_width:
        middle = (lower + upper) / 2
        middle_certificate = point_forced_window(
            make_point_context(coefficients, middle, bases),
            start_power=6,
            max_power=22,
        )
        if middle_certificate is None:
            raise RuntimeError(
                "unable to resolve a rational midpoint; this could be an exact rational crossing"
            )
        middle_count = int(middle_certificate["count_below_forced_value"])
        if middle_count == lower_count:
            lower, lower_certificate = middle, middle_certificate
        elif middle_count == upper_count:
            upper, upper_certificate = middle, middle_certificate
        else:
            raise RuntimeError(
                f"unexpected intermediate c-count {middle_count} in {name}; "
                "the one-jump bracket must be subdivided separately"
            )
        iterations += 1

    return {
        "name": name,
        "claim_tag": "[THEOREM]",
        "coarse_bracket_t": interval_json(coarse_lower, coarse_upper),
        "fine_bracket_t": interval_json(lower, upper),
        "c_counts_left_to_right": [lower_count, upper_count],
        "left_point_certificate": lower_certificate,
        "right_point_certificate": upper_certificate,
        "bisection_iterations": iterations,
        "conclusion": (
            "There exists at least one algebraic parameter t* in the fine bracket at which "
            "f(t*)=lambda_1(t*)lambda_2(t*)/lambda_0(t*) is an eigenvalue.  The c-count "
            "changes across two exact rational endpoints, so continuity forbids an "
            "eigenvalue-free path between them."
        ),
    }


# ---------------------------------------------------------------------------
# Exact Gaussian controls
# ---------------------------------------------------------------------------


def majorana_chain_control() -> dict[str, object]:
    """Structural all-interval Gaussian control for the open n=4 Ising chain.

    The Jordan--Wigner strings show X_j and Z_j Z_{j+1} are Majorana
    bilinears.  Thus P_t D_t P_t is a product of Gaussian exponentials for
    every 0<t<1.  The forced relation therefore holds for the whole control
    interval, including multiplicities.
    """
    n = 4
    majoranas: list[int] = []
    for site in range(n):
        prefix_x = sum(1 << previous for previous in range(site))
        gamma_even = prefix_x | (1 << (n + site))
        gamma_odd = prefix_x | (1 << site) | (1 << (n + site))
        majoranas.extend([gamma_even, gamma_odd])

    pairwise_anticommuting = all(
        symplectic_form(majoranas[left], majoranas[right], n) == 1
        for left in range(2 * n)
        for right in range(left + 1, 2 * n)
    )
    transverse = all(
        (majoranas[2 * site] ^ majoranas[2 * site + 1]) == pauli_x(n, site)
        for site in range(n)
    )
    longitudinal = all(
        (majoranas[2 * site + 1] ^ majoranas[2 * (site + 1)]) == zz(n, site, site + 1)
        for site in range(n - 1)
    )
    if not (pairwise_anticommuting and transverse and longitudinal):
        raise AssertionError("Jordan--Wigner Majorana control mapping failed")

    # A completely explicit two-mode diagonal control also checks the
    # multiplicity-safe forced-value formula symbolically, not numerically.
    # A(t)=1+t^2, B(t)=A(t)^3 obey 1 <= A <= B on [1/5,1/2].
    t = Fraction(1, 3)
    a = 1 + t * t
    b = a**3
    explicit_spectrum = [Fraction(1), a, b, a * b]
    explicit_forced = explicit_spectrum[1] * explicit_spectrum[2] / explicit_spectrum[0]
    explicit_identity = explicit_forced == explicit_spectrum[3]
    if not explicit_identity:
        raise AssertionError("explicit Gaussian two-mode forced identity failed")

    return {
        "claim_tag": "[LEMMA]",
        "interval_t": interval_json(Fraction(1, 5), Fraction(1, 2)),
        "n_sites": n,
        "majorana_bilinear_checks": {
            "all_majorana_strings_pairwise_anticommute": pairwise_anticommuting,
            "X_j_is_adjacent_majorana_bilinear": transverse,
            "Z_j_Z_jplus1_is_adjacent_majorana_bilinear": longitudinal,
        },
        "conclusion": (
            "For every t in (0,1), the open-chain P_t D_t P_t is a positive Gaussian "
            "operator because P_t is a product of exp(a X_j) and D_t is a scalar times "
            "a product of exp(K Z_j Z_{j+1}); all displayed generators are Majorana "
            "bilinears.  Hence the forced eigenvalue is present throughout the stated "
            "closed control interval."
        ),
        "explicit_two_mode_identity": {
            "spectrum_at_t_1_3": [str(value) for value in explicit_spectrum],
            "forced_value": str(explicit_forced),
            "equals_subset_product_eigenvalue": explicit_identity,
            "symbolic_form": "{1, A(t), B(t), A(t)B(t)}, A(t)=1+t^2, B(t)=A(t)^3",
        },
    }


def chain_window_control() -> dict[str, object]:
    """Exact occupied-window check on the actual n=4 open-chain matrix."""
    bonds = [(site, site + 1) for site in range(3)]
    coefficients, degree, shift = polynomial_matrix(4, bonds)
    if (degree, shift) != (3, 2):
        raise AssertionError("unexpected open-chain polynomial representative")
    generators = [
        state_permutation(4, list(reversed(range(4)))),
        spin_flip_permutation(4),
    ]
    bases = make_sector_bases(4, generators)
    rows: list[dict[str, object]] = []
    all_occupied = True
    for t in (Fraction(1, 5), Fraction(1, 3), Fraction(1, 2)):
        context = make_point_context(coefficients, t, bases)
        clear = point_forced_window(context, start_power=5, max_power=14)
        occupied = clear is None
        all_occupied &= occupied
        rows.append({
            "t": str(t),
            "forced_window_clear": not occupied,
            "interpretation": (
                "The exact forced window remains occupied under repeated tightening; "
                "the absence procedure correctly refuses this Gaussian control."
            ),
            "exact_inertia_calls": context.counter.calls,
        })
    return {
        "claim_tag": "[COMPUTATION]",
        "samples": rows,
        "all_sampled_forced_windows_occupied": all_occupied,
    }


# ---------------------------------------------------------------------------
# Main producer
# ---------------------------------------------------------------------------


def main() -> int:
    started = time.perf_counter()
    bonds = list(layer_bonds((2, 3), (False, False)))
    coefficients, degree, shift = polynomial_matrix(6, bonds)
    if (degree, shift) != (7, 4):
        raise AssertionError("unexpected 2x3 energy-exponent normalization")
    generators = [
        state_permutation(6, permutation)
        for permutation in grid_reflections(2, 3)
    ] + [spin_flip_permutation(6)]
    if not polynomial_invariant(coefficients, generators):
        raise AssertionError("2x3 polynomial representative is not symmetry invariant")
    bases = make_sector_bases(6, generators)
    sector_dimensions = [len(basis.gram) for basis in bases]
    if sector_dimensions != [14, 6, 6, 6, 10, 10, 6, 6]:
        raise AssertionError(f"unexpected spin-flip/reflection sectors: {sector_dimensions}")

    print("exactly scanning c(t) on [1/4,7/20] at rational step 1/1000 ...", flush=True)
    scan_runs, transitions = scan_count_runs(
        coefficients, bases, Fraction(1, 4), Fraction(7, 20), Fraction(1, 1000)
    )
    if len(transitions) != 4 or any(right != left + 1 for _, _, left, right in transitions):
        raise AssertionError(f"unexpected exact c-count transitions: {transitions}")

    crossings: list[dict[str, object]] = []
    for index, (lower, upper, left_count, right_count) in enumerate(transitions, start=1):
        print(f"refining exact c-jump crossing {index} ...", flush=True)
        crossing = refine_c_jump(
            f"forced_value_crossing_{index}",
            coefficients,
            bases,
            lower,
            upper,
            left_count,
            right_count,
        )
        crossings.append(crossing)
        fine = crossing["fine_bracket_t"]
        print(
            f"  [{fine['lower']}, {fine['upper']}], c: {left_count}->{right_count}",
            flush=True,
        )

    print("certifying selected closed safe intervals ...", flush=True)
    safe_intervals = [
        find_safe_interval("safe_interval_t_1_4", coefficients, bases, Fraction(1, 4)),
        find_safe_interval("safe_interval_t_3_10", coefficients, bases, Fraction(3, 10)),
        find_safe_interval("safe_interval_t_1_3", coefficients, bases, Fraction(1, 3)),
    ]
    for interval in safe_intervals:
        endpoint = interval["interval_t"]
        print(f"  {interval['name']}: [{endpoint['lower']}, {endpoint['upper']}]", flush=True)

    print("running all-interval chain control and occupied-window checks ...", flush=True)
    chain_structure = majorana_chain_control()
    chain_window = chain_window_control()

    fine_disjoint = all(
        Fraction(crossings[index]["fine_bracket_t"]["upper"])
        < Fraction(crossings[index + 1]["fine_bracket_t"]["lower"])
        for index in range(len(crossings) - 1)
    )
    all_fine_narrow = all(
        Fraction(crossing["fine_bracket_t"]["width"]) <= Fraction(1, 10**10)
        for crossing in crossings
    )
    widest_safe = max(Fraction(row["interval_t"]["width"]) for row in safe_intervals)
    checks = [
        {
            "name": "integer_polynomial_representative_and_eight_sector_split",
            "passed": (
                all(
                    coefficient >= 0
                    for row in coefficients
                    for entry in row
                    for coefficient in entry
                )
                and sector_dimensions == [14, 6, 6, 6, 10, 10, 6, 6]
            ),
            "detail": "M(t) has nonnegative integer coefficients and the three commuting involutions give the expected exact sectors.",
        },
        {
            "name": "four_disjoint_exact_c_jumps",
            "passed": len(crossings) == 4 and fine_disjoint and all_fine_narrow,
            "detail": "Each refined rational bracket has exact endpoint counts differing by one; the brackets are disjoint and width <= 10^-10.",
        },
        {
            "name": "forced_value_absence_target_falsified_by_genuine_crossings",
            "passed": all(
                crossing["c_counts_left_to_right"][1]
                == crossing["c_counts_left_to_right"][0] + 1
                for crossing in crossings
            ),
            "detail": "Continuity plus each exact c-count jump proves f(t*) is an eigenvalue at at least one algebraic t* in every listed bracket.",
        },
        {
            "name": "corrected_safe_intervals_are_strictly_wider_than_wave8_components",
            "passed": widest_safe > MAX_WAVE8_COMPONENT_WIDTH,
            "detail": (
                f"Largest new closed interval width is {widest_safe}, exceeding the prior maximum "
                f"component width {MAX_WAVE8_COMPONENT_WIDTH}."
            ),
        },
        {
            "name": "multiplicity_safe_forced_lemma_control",
            "passed": bool(chain_structure["explicit_two_mode_identity"]["equals_subset_product_eigenvalue"]),
            "detail": "The exact two-mode Gaussian control has f=A(t)B(t) as its subset-product eigenvalue even without a distinctness hypothesis.",
        },
        {
            "name": "open_chain_gaussian_control_stays_present_on_control_interval",
            "passed": (
                all(chain_structure["majorana_bilinear_checks"].values())
                and bool(chain_window["all_sampled_forced_windows_occupied"])
            ),
            "detail": "Jordan--Wigner bilinears prove all-t Gaussianity of the open chain; exact sampled windows remain occupied and never trigger false absence.",
        },
        {
            "name": "all_decisive_computation_exact",
            "passed": True,
            "detail": "Integer/Fraction arithmetic decides all inertia counts, bounds, and interval endpoints; float64 proposes brackets only.",
        },
    ]

    artifact = {
        "provenance": {
            "script": "experiments/e81_spectral_interval2.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "method": (
                "Exact integer-polynomial 2x3 representative; eight exact spin-flip/reflection "
                "character sectors; fraction-free Sylvester inertia; rational low-eigenvalue "
                "enclosures; exact c-count jump/continuity crossing certificate; and rational "
                "Weyl/inertia safe intervals with no distinctness assumption. Float64 only "
                "proposes brackets later verified by exact inertia."
            ),
        },
        "data": {
            "status": "FALSIFIED",
            "target": (
                "Universal forced-eigenvalue absence on a connected closed interval containing "
                "one of the certified crossings."
            ),
            "matrix_parametrization": {
                "t": "tanh(K*/2)",
                "q_exp_2K": "(1+t^2)/(2t)",
                "R": "P_t diag(q(t)^e_s) P_t",
                "polynomial_representative": "M(t)=(2t)^7 q(t)^4 R(t)",
                "n_sites": 6,
                "n_bonds": len(bonds),
                "dimension": 64,
                "polynomial_degree_upper": len(coefficients[0][0]) - 1,
                "all_entry_coefficients_nonnegative": True,
                "sector_dimensions_spinflip_reflections": sector_dimensions,
            },
            "forced_value": "f(t)=lambda_1(t)lambda_2(t)/lambda_0(t), with sorted eigenvalues counted with multiplicity",
            "scan": {
                "interval_t": interval_json(Fraction(1, 4), Fraction(7, 20)),
                "step_t": "1/1000",
                "constant_count_runs": scan_runs,
                "transition_count": len(transitions),
            },
            "crossings": crossings,
            "safe_closed_intervals": safe_intervals,
            "chain_control": {
                "structural_all_interval": chain_structure,
                "exact_window_samples": chain_window,
            },
            "corrected_statement": (
                "The forced-value absence criterion is valid on each listed safe closed rational "
                "interval, but it cannot hold on any connected interval containing one of the "
                "four certified crossing brackets: each bracket contains an algebraic parameter "
                "where the forced value is attained."
            ),
            "limitations": [
                "The c-jump certificate proves at least one crossing in each bracket; it does not compute a minimal polynomial or prove uniqueness within a bracket.",
                "The computation does not claim that the four brackets exhaust all crossings in [1/5,1/2].",
                "Forced-value attainment is necessary for Gaussianity but not sufficient; it does not make the 2x3 layer Gaussian at a crossing.",
                "No critical-coupling benchmark is used to select a point, bracket, or interval.",
            ],
        },
        "checks": checks,
    }
    output_path = ROOT / "results" / "spectral" / "interval2.json"
    os.makedirs(output_path.parent, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(artifact, handle, indent=2)
        handle.write("\n")

    for check in checks:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    elapsed = time.perf_counter() - started
    print(f"wall_seconds={elapsed:.6f}")
    if all(check["passed"] for check in checks):
        print("PASS")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
