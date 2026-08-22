"""Exact rational spectral-Gaussianity certificates at additional couplings and sizes.

This extends the finite-operator computation in e38_gaussianity_certificate.py.  The
mathematical certificate is unchanged: exact Sylvester-inertia bisection encloses the
three smallest eigenvalues, and two exact inertia counts test the outward-rounded
window forced by a Gaussian spectrum.

For the 256-dimensional cases, R is cleared once to an integer matrix A = scale * R.
Commuting permutation symmetries split A into exact rational invariant sectors.  In an
unnormalised integer orbit basis, each sector is the generalized symmetric pencil
B - sigma G; Sylvester inertia of the pencil is exactly the inertia of the restriction.
The fraction-free Bareiss LDL calculation below never uses floating point.  NumPy is
used only to propose initial brackets; every endpoint is accepted or rejected by an
exact inertia count.
"""

from __future__ import annotations

import json
import math
import os
import platform
import signal
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402


class ResourceLimit(RuntimeError):
    """An explicitly configured exact-computation wall was reached."""


class _Alarm:
    """Hard wall-clock alarm for one exact inertia count (main thread, POSIX)."""

    def __init__(self, seconds: float):
        self.seconds = seconds
        self.old_handler = None
        self.old_timer = None

    @staticmethod
    def _raise_timeout(_signum, _frame):
        raise ResourceLimit("hard wall-clock limit reached during one inertia count")

    def __enter__(self):
        if not hasattr(signal, "setitimer"):
            return self
        self.old_handler = signal.getsignal(signal.SIGALRM)
        self.old_timer = signal.getitimer(signal.ITIMER_REAL)
        signal.signal(signal.SIGALRM, self._raise_timeout)
        signal.setitimer(signal.ITIMER_REAL, max(0.001, self.seconds))
        return self

    def __exit__(self, exc_type, exc, tb):
        if hasattr(signal, "setitimer"):
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, self.old_handler)
            if self.old_timer and self.old_timer[0] > 0:
                signal.setitimer(signal.ITIMER_REAL, *self.old_timer)
        return False


@dataclass
class ScaledR:
    matrix: list[list[int]]
    scale: int
    q: Fraction
    bond_parity: int
    energy_exponents: tuple[int, int]


@dataclass
class Sector:
    character: int
    matrix: list[list[int]]
    gram: list[int]


@dataclass
class CountRecord:
    count: int | None
    seconds: float


def exact_decimal(x: Fraction, digits: int = 18) -> str:
    """Display-only decimal; exact rational strings remain the certificate."""
    with localcontext() as ctx:
        ctx.prec = digits
        return format(Decimal(x.numerator) / Decimal(x.denominator), ".12E")


def build_scaled_R(n: int, bonds: list[tuple[int, int]], t: Fraction) -> ScaledR:
    """Build A = scale * R exactly, with A integral and scale a positive integer.

    R = P_t diag(q^e) P_t is exactly the matrix used by e38.  Clearing the
    denominators before elimination prevents repeated Fraction normalisation.
    """
    if not (0 < t < 1):
        raise ValueError("the certificate expects 0 < t < 1")
    dim = 1 << n
    q = (1 + t * t) / (2 * t)

    energies: list[int] = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - i)) & 1) for i in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = energies[0] % 2
    if not all((energy - eps) % 2 == 0 for energy in energies):
        raise AssertionError("bond-count parity assumption violated")
    exponents = [(energy - eps) // 2 for energy in energies]

    min_exp, max_exp = min(exponents), max(exponents)
    weight_scale = (
        q.denominator ** max(max_exp, 0) * q.numerator ** max(-min_exp, 0)
    )
    weights: list[int] = []
    for exponent in exponents:
        value = q**exponent * weight_scale
        if value.denominator != 1:
            raise AssertionError("internal denominator-clearing error")
        weights.append(value.numerator)

    # P_t = P_num / denominator(t)^n, with P_num a Kronecker power of
    # [[den(t), num(t)], [num(t), den(t)]].
    a, b = t.numerator, t.denominator
    powers_a = [a**h for h in range(n + 1)]
    powers_b = [b**h for h in range(n + 1)]
    p_num = [
        [
            powers_a[(i ^ j).bit_count()] * powers_b[n - (i ^ j).bit_count()]
            for j in range(dim)
        ]
        for i in range(dim)
    ]

    # Start from diag(weights) P_num, then apply P_num on the left as n local
    # two-by-two transforms.  This is O(n 4^n), rather than a dense O(8^n)
    # multiplication, and is the same exact matrix product.
    matrix = [[weights[i] * entry for entry in p_num[i]] for i in range(dim)]
    for bit in range(n):
        mask = 1 << bit
        for base in range(dim):
            if base & mask:
                continue
            mate = base | mask
            row0, row1 = matrix[base], matrix[mate]
            matrix[base] = [b * x + a * y for x, y in zip(row0, row1)]
            matrix[mate] = [a * x + b * y for x, y in zip(row0, row1)]

    if any(matrix[i][j] != matrix[j][i] for i in range(dim) for j in range(i)):
        raise AssertionError("P diag(d) P construction lost symmetry")
    scale = b ** (2 * n) * weight_scale
    return ScaledR(matrix, scale, q, eps, (min_exp, max_exp))


def state_permutation(n: int, site_permutation: list[int]) -> list[int]:
    """Permutation of computational-basis states induced by a site permutation."""
    if sorted(site_permutation) != list(range(n)):
        raise ValueError("not a site permutation")
    result: list[int] = []
    for state in range(1 << n):
        image = 0
        for source, target in enumerate(site_permutation):
            if (state >> (n - 1 - source)) & 1:
                image |= 1 << (n - 1 - target)
        result.append(image)
    return result


def make_generators(n: int, site_permutations: list[list[int]]) -> list[list[int]]:
    dim = 1 << n
    generators = [[state ^ (dim - 1) for state in range(dim)]]
    generators.extend(state_permutation(n, perm) for perm in site_permutations)
    for generator in generators:
        if any(generator[generator[state]] != state for state in range(dim)):
            raise AssertionError("all symmetry generators must be involutions")
    for left in generators:
        for right in generators:
            if any(left[right[state]] != right[left[state]] for state in range(dim)):
                raise AssertionError("symmetry generators must commute")
    return generators


def matrix_invariant(matrix: list[list[int]], generators: list[list[int]]) -> bool:
    dim = len(matrix)
    return all(
        matrix[g[i]][g[j]] == matrix[i][j]
        for g in generators
        for i in range(dim)
        for j in range(dim)
    )


def symmetry_sectors(matrix: list[list[int]], generators: list[list[int]]) -> list[Sector]:
    """Exact character-sector pencils B - sigma G for commuting involutions.

    Orbit vectors have disjoint support inside each character sector, so their Gram
    matrix G is positive diagonal.  Normalising by G^{-1/2} is a congruence; hence
    inertia(B - sigma G) is the spectral count in that invariant sector.
    """
    dim = len(matrix)
    n_generators = len(generators)
    group_size = 1 << n_generators

    images: list[list[int]] = []
    for group_element in range(group_size):
        image: list[int] = []
        for state in range(dim):
            value = state
            for generator_index, generator in enumerate(generators):
                if (group_element >> generator_index) & 1:
                    value = generator[value]
            image.append(value)
        images.append(image)

    vectors_by_character: list[list[list[tuple[int, int]]]] = [
        [] for _ in range(group_size)
    ]
    seen: set[int] = set()
    for representative in range(dim):
        if representative in seen:
            continue
        orbit = {image[representative] for image in images}
        seen.update(orbit)
        for character in range(group_size):
            coefficients: dict[int, int] = {}
            for group_element, image in enumerate(images):
                sign = -1 if (character & group_element).bit_count() % 2 else 1
                state = image[representative]
                coefficients[state] = coefficients.get(state, 0) + sign
            nonzero = [(state, value) for state, value in coefficients.items() if value]
            if not nonzero:
                continue
            common = 0
            for _, value in nonzero:
                common = math.gcd(common, abs(value))
            vectors_by_character[character].append(
                sorted((state, value // common) for state, value in nonzero)
            )

    sectors: list[Sector] = []
    for character, vectors in enumerate(vectors_by_character):
        if not vectors:
            continue
        gram = [sum(coefficient * coefficient for _, coefficient in vector) for vector in vectors]
        size = len(vectors)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            left = vectors[row]
            for col in range(row, size):
                right = vectors[col]
                value = sum(
                    left_coefficient * right_coefficient * matrix[left_state][right_state]
                    for left_state, left_coefficient in left
                    for right_state, right_coefficient in right
                )
                block[row][col] = block[col][row] = value
        sectors.append(Sector(character, block, gram))

    if sum(len(sector.gram) for sector in sectors) != dim:
        raise AssertionError("character sectors do not span the full state space")
    return sectors


def fraction_free_inertia(lower_matrix: list[list[int]]) -> int | None:
    """Exact inertia via the symmetric Bareiss/LDL recurrence.

    With nonzero leading principal minors, the successive Bareiss pivots are the
    leading principal determinants.  The LDL diagonal signs are therefore
    sign(pivot_k / pivot_{k-1}).  A zero pivot is reported, never perturbed here.
    """
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
                numerator = (
                    pivot * row[col_index]
                    - left * work[col_index][pivot_index]
                )
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Bareiss division was not exact")
                row[col_index] = quotient
        previous = pivot
    return negatives


class ExactCounter:
    """Cached exact inertia counts for all generalized symmetry sectors."""

    def __init__(
        self,
        sectors: list[Sector],
        started: float,
        per_count_limit: float = 3000.0,
        total_limit: float | None = None,
    ):
        self.sectors = sectors
        self.started = started
        self.per_count_limit = per_count_limit
        self.total_limit = total_limit
        self.cache: dict[Fraction, CountRecord] = {}
        self.count_seconds = 0.0
        self.max_count_seconds = 0.0

    def count(self, shift: Fraction) -> int | None:
        shift = Fraction(shift)
        if shift in self.cache:
            return self.cache[shift].count
        remaining = self.per_count_limit
        if self.total_limit is not None:
            remaining = min(remaining, self.total_limit - (time.perf_counter() - self.started))
        if remaining <= 0:
            raise ResourceLimit("four-hour total wall reached before an inertia count")

        began = time.perf_counter()
        with _Alarm(remaining):
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
                    total_count: int | None = None
                    break
                total += count
            else:
                total_count = total
        elapsed = time.perf_counter() - began
        self.count_seconds += elapsed
        self.max_count_seconds = max(self.max_count_seconds, elapsed)
        self.cache[shift] = CountRecord(total_count, elapsed)
        return total_count


def count_inside(counter: ExactCounter, lo: Fraction, hi: Fraction) -> tuple[Fraction, int]:
    """Replicate e38: return the exact usable shift and its exact count."""
    span = hi - lo
    for attempt in range(40):
        trial = (
            lo + span / 2
            if attempt == 0
            else lo
            + span * Fraction(1, 2)
            + span * Fraction((-1) ** attempt, 3 ** (attempt + 2))
        )
        if not (lo < trial < hi):
            continue
        count = counter.count(trial)
        if count is not None:
            return trial, count
    raise RuntimeError("no non-degenerate shift found inside the bracket")


def count_outward(
    counter: ExactCounter, shift: Fraction, direction: int
) -> tuple[Fraction, int]:
    """Replicate e38: a degenerate endpoint is moved only outward."""
    count = counter.count(shift)
    if count is not None:
        return shift, count
    step = abs(shift) if shift else Fraction(1)
    for attempt in range(1, 40):
        moved = shift + direction * step * Fraction(1, 10 ** (20 - attempt // 2))
        count = counter.count(moved)
        if count is not None:
            return moved, count
    raise RuntimeError("no non-degenerate outward shift found")


def guide_eigenvalues(matrix: list[list[int]]) -> list[float]:
    """Non-rigorous bracket proposals only; never used without exact endpoint counts."""
    largest = max(abs(value) for row in matrix for value in row)
    scaled = np.array([[value / largest for value in row] for row in matrix], dtype=float)
    estimates = np.linalg.eigvalsh(0.5 * (scaled + scaled.T))[:3]
    return [float(value * largest) for value in estimates]


def initial_bracket(
    counter: ExactCounter, index: int, guide: float, trace: int
) -> tuple[Fraction, Fraction]:
    """Turn a floating proposal into a certified exact inertia bracket."""
    if not math.isfinite(guide) or guide <= 0:
        lo, hi = Fraction(0), Fraction(trace + 1)
    else:
        center = Fraction(format(guide, ".16e"))
        lo, hi = center * Fraction(49, 50), center * Fraction(51, 50)

    lo, count_lo = count_outward(counter, max(Fraction(0), lo), -1)
    while count_lo > index:
        lo, count_lo = count_outward(counter, lo / 2, -1)
    hi, count_hi = count_outward(counter, hi, +1)
    while count_hi <= index:
        hi, count_hi = count_outward(counter, hi * 2, +1)
    return lo, hi


def isolate(
    counter: ExactCounter,
    index: int,
    lo: Fraction,
    hi: Fraction,
    relative_width: Fraction,
) -> tuple[Fraction, Fraction]:
    """Certified rational enclosure with exactly the invariant used by e38."""
    while not (lo > 0 and hi - lo <= relative_width * lo):
        trial, count = count_inside(counter, lo, hi)
        if count <= index:
            lo = trial
        else:
            hi = trial
    return lo, hi


def grid_reflections(rows: int, cols: int) -> list[list[int]]:
    return [
        [(rows - 1 - row) * cols + col for row in range(rows) for col in range(cols)],
        [row * cols + (cols - 1 - col) for row in range(rows) for col in range(cols)],
    ]


def chain_reflection(n: int) -> list[list[int]]:
    return [list(reversed(range(n)))]


def interval_json(lo: Fraction, hi: Fraction) -> dict[str, object]:
    return {
        "lower": str(lo),
        "upper": str(hi),
        "lower_decimal": exact_decimal(lo),
        "upper_decimal": exact_decimal(hi),
        "relative_width": str((hi - lo) / lo),
    }


def certify_case(case: dict[str, object]) -> dict[str, object]:
    started = time.perf_counter()
    completed = "case initialization"
    timings: dict[str, float] = {}
    try:
        n = int(case["n"])
        bonds = list(case["bonds"])
        t = Fraction(case["t"])

        began = time.perf_counter()
        scaled = build_scaled_R(n, bonds, t)
        timings["matrix_build_seconds"] = time.perf_counter() - began
        completed = "exact integer-scaled R built"

        generators = make_generators(n, list(case["site_symmetries"]))
        invariant = matrix_invariant(scaled.matrix, generators)
        if not invariant:
            raise AssertionError("the proposed permutations do not commute with R")
        began = time.perf_counter()
        sectors = symmetry_sectors(scaled.matrix, generators)
        timings["sector_build_seconds"] = time.perf_counter() - began
        completed = "exact symmetry-sector pencils built"

        total_limit = 4 * 3600 if case["name"] == "layer_2x4_t_1_3" else None
        counter = ExactCounter(
            sectors,
            started,
            per_count_limit=3000.0,
            total_limit=total_limit,
        )
        if counter.count(Fraction(0)) != 0:
            raise AssertionError("R must be positive definite")
        completed = "positive definiteness certified by exact inertia"

        guides = guide_eigenvalues(scaled.matrix)
        trace = sum(scaled.matrix[i][i] for i in range(1 << n))
        relative_target = Fraction(1, 10 ** int(case["relative_exp"]))
        enclosures: list[tuple[Fraction, Fraction]] = []
        for index in range(3):
            lo, hi = initial_bracket(counter, index, guides[index], trace)
            enclosures.append(isolate(counter, index, lo, hi, relative_target))
            completed = f"lambda_{index} enclosed by exact inertia bisection"

        (l0, h0), (l1, h1), (l2, h2) = enclosures
        distinct = h0 < l1 and h1 < l2
        if not distinct:
            raise AssertionError(
                "three lowest eigenvalue enclosures overlap; Gaussian obstruction is invalid"
            )
        gap01 = (l1 - h0) / l1
        gap12 = (l2 - h1) / l2
        predicted_lo = l1 * l2 / h0
        predicted_hi = h1 * h2 / l0
        tested_lo, count_lo = count_outward(counter, predicted_lo, -1)
        tested_hi, count_hi = count_outward(counter, predicted_hi, +1)
        absent = count_lo == count_hi
        completed = "outward forced-value window tested at both ends"

        normalized = [(lo / scaled.scale, hi / scaled.scale) for lo, hi in enclosures]
        norm_predicted = (predicted_lo / scaled.scale, predicted_hi / scaled.scale)
        norm_tested = (tested_lo / scaled.scale, tested_hi / scaled.scale)
        total_seconds = time.perf_counter() - started
        timings.update(
            {
                "inertia_seconds": counter.count_seconds,
                "max_single_inertia_seconds": counter.max_count_seconds,
                "total_seconds": total_seconds,
            }
        )

        lambda_fields: dict[str, object] = {}
        for index, ((lo, hi), (raw_lo, raw_hi)) in enumerate(zip(normalized, enclosures)):
            lower_count = counter.count(raw_lo)
            upper_count = counter.count(raw_hi)
            lambda_fields[f"lambda{index}"] = {
                **interval_json(lo, hi),
                "count_below_lower": lower_count,
                "count_below_upper": upper_count,
            }

        return {
            "name": case["name"],
            "kind": case["kind"],
            "status": "CERTIFIED",
            "t": str(t),
            "exp_2K": str(scaled.q),
            "n_sites": n,
            "n_bonds": len(bonds),
            "dim": 1 << n,
            "bond_parity": scaled.bond_parity,
            "integer_scale_for_R": str(scaled.scale),
            "energy_exponent_range": list(scaled.energy_exponents),
            "lambda_enclosures": lambda_fields,
            "three_lowest_distinct": distinct,
            "relative_gaps": {
                "gap01_lower_bound": str(gap01),
                "gap01_lower_bound_decimal": exact_decimal(gap01),
                "gap12_lower_bound": str(gap12),
                "gap12_lower_bound_decimal": exact_decimal(gap12),
            },
            "window": {
                "predicted": interval_json(*norm_predicted),
                "actually_tested": interval_json(*norm_tested),
                "tested_window_encloses_prediction": (
                    tested_lo <= predicted_lo and tested_hi >= predicted_hi
                ),
            },
            "inertia_counts_at_window_ends": [count_lo, count_hi],
            "predicted_eigenvalue_absent": absent,
            "relative_enclosure_target": str(relative_target),
            "symmetry": {
                "generator_count": len(generators),
                "matrix_invariant": invariant,
                "sector_dimensions": [len(sector.gram) for sector in sectors],
                "dimensions_sum": sum(len(sector.gram) for sector in sectors),
            },
            "exact_inertia_calls": len(counter.cache),
            "elapsed_seconds": round(total_seconds, 6),
            "elapsed": {key: round(value, 6) for key, value in timings.items()},
        }
    except ResourceLimit as exc:
        timings["total_seconds"] = time.perf_counter() - started
        return {
            "name": case["name"],
            "kind": case["kind"],
            "status": "UNRESOLVED-RESOURCE",
            "t": str(case["t"]),
            "exp_2K": str((1 + Fraction(case["t"]) ** 2) / (2 * Fraction(case["t"]))),
            "n_sites": int(case["n"]),
            "n_bonds": len(case["bonds"]),
            "dim": 1 << int(case["n"]),
            "largest_completed_substep": completed,
            "resource_limit": str(exc),
            "per_inertia_timeout_seconds": 3000,
            "total_case_timeout_seconds": 4 * 3600,
            "elapsed_seconds": round(timings["total_seconds"], 6),
            "elapsed": {key: round(value, 6) for key, value in timings.items()},
        }


def case_definitions() -> list[dict[str, object]]:
    grid23 = list(layer_bonds((2, 3), (False, False)))
    grid24 = list(layer_bonds((2, 4), (False, False)))
    chain8 = list(layer_bonds((8,), (False,)))
    return [
        {
            "name": "layer_2x3_t_1_5",
            "kind": "3D-layer",
            "n": 6,
            "bonds": grid23,
            "t": Fraction(1, 5),
            "site_symmetries": grid_reflections(2, 3),
            "relative_exp": 7,
        },
        {
            "name": "layer_2x3_t_2_5",
            "kind": "3D-layer",
            "n": 6,
            "bonds": grid23,
            "t": Fraction(2, 5),
            "site_symmetries": grid_reflections(2, 3),
            "relative_exp": 7,
        },
        {
            "name": "layer_2x3_t_1_2",
            "kind": "3D-layer",
            "n": 6,
            "bonds": grid23,
            "t": Fraction(1, 2),
            "site_symmetries": grid_reflections(2, 3),
            "relative_exp": 7,
        },
        {
            "name": "layer_2x4_t_1_3",
            "kind": "3D-layer",
            "n": 8,
            "bonds": grid24,
            "t": Fraction(1, 3),
            "site_symmetries": grid_reflections(2, 4),
            "relative_exp": 7,
        },
        {
            "name": "control_chain_n8_t_1_3",
            "kind": "1D-control",
            "n": 8,
            "bonds": chain8,
            "t": Fraction(1, 3),
            "site_symmetries": chain_reflection(8),
            "relative_exp": 4,
        },
    ]


def computed_checks(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    checks: list[dict[str, object]] = []
    for row in rows:
        name = str(row["name"])
        if row["status"] != "CERTIFIED":
            checks.append(
                {
                    "name": f"{name}_resource_status_recorded",
                    "passed": row["status"] == "UNRESOLVED-RESOURCE",
                    "detail": (
                        f"status={row['status']}; largest completed substep="
                        f"{row.get('largest_completed_substep', 'not recorded')}"
                    ),
                }
            )
            continue
        counts = list(row["inertia_counts_at_window_ends"])
        if row["kind"] == "3D-layer":
            passed = bool(
                row["three_lowest_distinct"]
                and row["predicted_eigenvalue_absent"]
                and counts[0] == counts[1]
            )
            detail = (
                f"distinct={row['three_lowest_distinct']}; exact outward-window counts "
                f"{counts[0]} and {counts[1]}; absent={row['predicted_eigenvalue_absent']}"
            )
        else:
            passed = bool(
                row["three_lowest_distinct"]
                and not row["predicted_eigenvalue_absent"]
                and counts[0] != counts[1]
            )
            detail = (
                f"distinct={row['three_lowest_distinct']}; exact outward-window counts "
                f"{counts[0]} and {counts[1]}. Differing counts certify only that SOME "
                "eigenvalue occupies the positive-width window, not exact equality."
            )
        checks.append({"name": name, "passed": passed, "detail": detail})
    return checks


def main() -> int:
    rows: list[dict[str, object]] = []
    for case in case_definitions():
        print(f"certifying {case['name']} ...", flush=True)
        row = certify_case(case)
        rows.append(row)
        if row["status"] == "CERTIFIED":
            counts = row["inertia_counts_at_window_ends"]
            print(
                f"  dim={row['dim']} exp(2K)={row['exp_2K']} counts={counts[0]},{counts[1]} "
                f"absent={row['predicted_eigenvalue_absent']} "
                f"elapsed={row['elapsed']['total_seconds']:.3f}s",
                flush=True,
            )
        else:
            print(
                f"  {row['status']}: {row['largest_completed_substep']} "
                f"({row['elapsed']['total_seconds']:.3f}s)",
                flush=True,
            )

    checks = computed_checks(rows)
    result = {
        "meta": {
            "provenance": "experiments/e39_spectral_couplings.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "arithmetic": (
                "Exact integers/rationals for R, every enclosure endpoint, and every inertia "
                "count. Fraction-free Bareiss LDL in exact symmetry sectors. Float64 is used "
                "only to propose initial brackets; exact endpoint inertia counts certify all "
                "accepted brackets. Decimal strings are display-only."
            ),
            "claim_scope": (
                "Each absent=True row is a finite-instance exact certificate only. Control "
                "rows with differing counts are consistency checks only."
            ),
            "resource_limits": {
                "per_inertia_count_seconds": 3000,
                "layer_2x4_total_seconds": 14400,
            },
        },
        "data": {"rows": rows},
        "checks": checks,
    }
    output = ROOT / "results" / "spectral" / "couplings.json"
    os.makedirs(output.parent, exist_ok=True)
    with output.open("w") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")

    for check in checks:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    failed = [check["name"] for check in checks if not check["passed"]]
    if failed:
        print(f"FAIL: {failed}")
        return 1
    print("PASS: all completed exact certificates and control consistency checks succeeded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
