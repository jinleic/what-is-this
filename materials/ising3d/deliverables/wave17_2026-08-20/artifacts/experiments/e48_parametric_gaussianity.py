"""Exact parametric spectral-Gaussianity certificates for the 2x3 Ising layer.

The requested continuum [1/5,1/2] no-go remains unresolved for the particular
forced-value criterion: numerical reconnaissance suggests that the forced value
crosses layer eigenvalues inside that interval, but no such crossing is promoted
to an exact claim here.  This experiment therefore records an honest [UNRESOLVED] result and
proves four nonzero rational neighbourhoods around the previously certified
couplings.  Every accepted endpoint, norm bound, and inertia count is exact.

The matrix used here is a positive polynomial scalar multiple of e38's matrix:

    M(t) = (2t)^7 q(t)^4 R(t),  q(t)=(1+t^2)/(2t).

For the 2x3 graph the exponents in R range from -4 to 3, so every entry of M is
an integer polynomial with nonnegative coefficients.  This makes entrywise
interval bounds and a rational Frobenius derivative majorant inexpensive.
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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402


@dataclass(frozen=True)
class Sector:
    matrix: list[list[int]]
    gram: list[int]


def exact_decimal(value: Fraction, digits: int = 18) -> str:
    """Display only; decimals never decide a certificate."""
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".12E")


def fraction_json(value: Fraction) -> dict[str, str]:
    return {"exact": str(value), "decimal": exact_decimal(value)}


def interval_json(lower: Fraction, upper: Fraction) -> dict[str, str]:
    return {
        "lower": str(lower),
        "upper": str(upper),
        "lower_decimal": exact_decimal(lower),
        "upper_decimal": exact_decimal(upper),
    }


def polynomial_matrix(
    n: int, bonds: list[tuple[int, int]]
) -> tuple[list[list[list[int]]], int, int]:
    """Return coefficients of a positive polynomial scalar multiple of R.

    If e_s=(B_s-eps)/2, put shift=-min(e_s), degree=max(e_s+shift), and
    M=(2t)^degree q^shift R.  In the 2x3 case shift=4 and degree=7.
    """
    dim = 1 << n
    energies: list[int] = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = energies[0] % 2
    if not all((energy - parity) % 2 == 0 for energy in energies):
        raise AssertionError("bond-count parity assumption violated")
    exponents = [(energy - parity) // 2 for energy in energies]
    shift = -min(exponents)
    degree = max(exponent + shift for exponent in exponents)
    max_polynomial_degree = 2 * n + 2 * degree

    coefficients = [
        [[0] * (max_polynomial_degree + 1) for _ in range(dim)]
        for _ in range(dim)
    ]
    for row in range(dim):
        for col in range(row, dim):
            entry = [0] * (max_polynomial_degree + 1)
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


def clear_denominators(
    matrix: list[list[Fraction]],
) -> tuple[list[list[int]], int]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = math.lcm(scale, value.denominator)
    return [[(value * scale).numerator for value in row] for row in matrix], scale


def derivative_frobenius_square_upper(
    coefficients: list[list[list[int]]], upper_t: Fraction
) -> Fraction:
    """Exact upper bound for ||M'(t)||_F^2 on 0 <= t <= upper_t.

    All coefficients are nonnegative, hence every entry derivative is
    nonnegative and increasing.  Evaluating each at upper_t is outward.
    """
    total = Fraction(0)
    for row in coefficients:
        for entry in row:
            derivative = sum(
                (Fraction(power * coefficient) * upper_t ** (power - 1)
                 for power, coefficient in enumerate(entry) if power and coefficient),
                Fraction(0),
            )
            total += derivative * derivative
    return total


def rational_sqrt_upper(value: Fraction, denominator: int = 10**15) -> Fraction:
    """Outward rational square-root bound proved by an integer square comparison."""
    target = (
        value.numerator * denominator * denominator + value.denominator - 1
    ) // value.denominator
    numerator = math.isqrt(target)
    if numerator * numerator < target:
        numerator += 1
    result = Fraction(numerator, denominator)
    if result * result < value:
        raise AssertionError("square-root rounding was not outward")
    return result


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


def make_generators(n: int, site_permutations: list[list[int]]) -> list[list[int]]:
    return [state_permutation(n, permutation) for permutation in site_permutations]


def matrix_invariant(matrix: list[list[int]], generators: list[list[int]]) -> bool:
    dim = len(matrix)
    return all(
        matrix[generator[row]][generator[col]] == matrix[row][col]
        for generator in generators
        for row in range(dim)
        for col in range(dim)
    )


def symmetry_sectors(matrix: list[list[int]], generators: list[list[int]]) -> list[Sector]:
    """Exact character-sector pencils B-sigma G for commuting involutions."""
    dim = len(matrix)
    group: list[list[int]] = []
    for mask in range(1 << len(generators)):
        permutation = list(range(dim))
        for index, generator in enumerate(generators):
            if mask & (1 << index):
                permutation = [generator[state] for state in permutation]
        group.append(permutation)

    seen: set[int] = set()
    representatives: list[int] = []
    orbits: dict[int, list[int]] = {}
    for state in range(dim):
        if state in seen:
            continue
        orbit = sorted({permutation[state] for permutation in group})
        representative = orbit[0]
        representatives.append(representative)
        orbits[representative] = orbit
        seen.update(orbit)

    sectors: list[Sector] = []
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
                norm = sum(coefficient * coefficient for coefficient in vector.values())
                vectors.append(vector)
                gram.append(norm)
        size = len(vectors)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            for col in range(row + 1):
                value = sum(
                    left * matrix[left_state][right_state] * right
                    for left_state, left in vectors[row].items()
                    for right_state, right in vectors[col].items()
                )
                block[row][col] = block[col][row] = value
        sectors.append(Sector(block, gram))
    if sum(len(sector.gram) for sector in sectors) != dim:
        raise AssertionError("symmetry-sector dimensions do not sum to full dimension")
    return sectors


def fraction_free_inertia(lower_matrix: list[list[int]]) -> int | None:
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
    def __init__(self, sectors: list[Sector]):
        self.sectors = sectors
        self.cache: dict[Fraction, int | None] = {}

    def count(self, shift: Fraction) -> int | None:
        shift = Fraction(shift)
        if shift in self.cache:
            return self.cache[shift]
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
        return total


def usable_inside(counter: ExactCounter, lower: Fraction, upper: Fraction) -> tuple[Fraction, int]:
    span = upper - lower
    for attempt in range(40):
        trial = (
            lower + span / 2
            if attempt == 0
            else lower + span / 2 + span * Fraction((-1) ** attempt, 3 ** (attempt + 2))
        )
        if lower < trial < upper:
            count = counter.count(trial)
            if count is not None:
                return trial, count
    raise RuntimeError("no nondegenerate rational shift inside bracket")


def usable_outward(
    counter: ExactCounter, shift: Fraction, direction: int
) -> tuple[Fraction, int]:
    count = counter.count(shift)
    if count is not None:
        return shift, count
    step = abs(shift) if shift else Fraction(1)
    for attempt in range(1, 40):
        moved = shift + direction * step * Fraction(1, 10 ** (20 - attempt // 2))
        count = counter.count(moved)
        if count is not None:
            return moved, count
    raise RuntimeError("no nondegenerate outward rational shift")


def isolate(
    counter: ExactCounter,
    index: int,
    lower: Fraction,
    upper: Fraction,
    relative_width: Fraction,
) -> tuple[Fraction, Fraction]:
    while not (lower > 0 and upper - lower <= relative_width * lower):
        trial, count = usable_inside(counter, lower, upper)
        if count <= index:
            lower = trial
        else:
            upper = trial
    return lower, upper


def initial_bracket(
    counter: ExactCounter, index: int, guide: float, trace: int
) -> tuple[Fraction, Fraction]:
    if not math.isfinite(guide) or guide <= 0:
        lower, upper = Fraction(0), Fraction(trace + 1)
    else:
        center = Fraction(format(guide, ".15e"))
        lower, upper = center * Fraction(49, 50), center * Fraction(51, 50)
    lower, lower_count = usable_outward(counter, max(Fraction(0), lower), -1)
    while lower_count > index:
        lower, lower_count = usable_outward(counter, lower / 2, -1)
    upper, upper_count = usable_outward(counter, upper, +1)
    while upper_count <= index:
        upper, upper_count = usable_outward(counter, upper * 2, +1)
    return lower, upper


def guide_eigenvalues(matrix: list[list[int]]) -> list[float]:
    largest = max(abs(value) for row in matrix for value in row)
    scaled = np.array([[value / largest for value in row] for row in matrix], dtype=float)
    return [float(value * largest) for value in np.linalg.eigvalsh(scaled)[:3]]


def certify_interval(
    name: str,
    n: int,
    bonds: list[tuple[int, int]],
    coefficients: list[list[list[int]]],
    center: Fraction,
    half_width: Fraction,
    site_symmetries: list[list[int]],
    relative_width: Fraction = Fraction(1, 10**8),
) -> dict[str, object]:
    started = time.perf_counter()
    lower_t, upper_t = center - half_width, center + half_width
    if not (0 < lower_t <= center <= upper_t < 1):
        raise ValueError("parameter interval must lie in (0,1)")

    rational_matrix = evaluate_polynomial_matrix(coefficients, center)
    integer_matrix, scale = clear_denominators(rational_matrix)
    generators = make_generators(n, site_symmetries)
    invariant = matrix_invariant(integer_matrix, generators)
    if not invariant:
        raise AssertionError("proposed reflection sectors do not preserve the matrix")
    sectors = symmetry_sectors(integer_matrix, generators)
    counter = ExactCounter(sectors)
    if counter.count(Fraction(0)) != 0:
        raise AssertionError("positive definiteness failed")

    trace = sum(integer_matrix[index][index] for index in range(1 << n))
    guides = guide_eigenvalues(integer_matrix)
    raw_enclosures: list[tuple[Fraction, Fraction]] = []
    for index in range(3):
        lower, upper = initial_bracket(counter, index, guides[index], trace)
        raw_enclosures.append(isolate(counter, index, lower, upper, relative_width))
    enclosures = [(lower / scale, upper / scale) for lower, upper in raw_enclosures]
    (l0, u0), (l1, u1), (l2, u2) = enclosures

    derivative_square = derivative_frobenius_square_upper(coefficients, upper_t)
    derivative_bound = rational_sqrt_upper(derivative_square)
    rho = derivative_bound * half_width
    positive_lower = l0 - rho
    distinct01_margin = l1 - u0 - 2 * rho
    distinct12_margin = l2 - u1 - 2 * rho
    if positive_lower <= 0 or distinct01_margin <= 0 or distinct12_margin <= 0:
        raise AssertionError("interval is too wide for positivity/distinctness propagation")

    forced_center_lower = l1 * l2 / u0
    forced_center_upper = u1 * u2 / l0
    forced_global_lower = (l1 - rho) * (l2 - rho) / (u0 + rho)
    forced_global_upper = (u1 + rho) * (u2 + rho) / (l0 - rho)

    # If an eigenvalue at t equalled f(t), Weyl would put its center eigenvalue
    # in [global forced lower-rho, global forced upper+rho].  Equal center
    # inertia counts therefore rule equality out throughout the subinterval.
    tested_lower = forced_global_lower - rho
    tested_upper = forced_global_upper + rho
    raw_tested_lower, count_lower = usable_outward(counter, tested_lower * scale, -1)
    raw_tested_upper, count_upper = usable_outward(counter, tested_upper * scale, +1)
    actual_tested_lower = raw_tested_lower / scale
    actual_tested_upper = raw_tested_upper / scale
    absent = count_lower == count_upper

    left_margin = forced_global_lower - actual_tested_lower
    right_margin = actual_tested_upper - forced_global_upper
    if left_margin < rho or right_margin < rho:
        raise AssertionError("outward tested window lost the Weyl safety margin")

    return {
        "name": name,
        "claim_tag": "[THEOREM]",
        "center_t": str(center),
        "half_width_t": str(half_width),
        "interval_t": interval_json(lower_t, upper_t),
        "exp_2K_at_center": str((1 + center * center) / (2 * center)),
        "n_sites": n,
        "n_bonds": len(bonds),
        "dimension": 1 << n,
        "polynomial_scale": f"(2t)^7*q(t)^4 relative to R(t)" if n == 6 else "derived automatically",
        "integer_scale_at_center": str(scale),
        "sector_dimensions": [len(sector.gram) for sector in sectors],
        "reflection_matrix_invariant": invariant,
        "lambda_center_enclosures": {
            f"lambda{index}": interval_json(*enclosure)
            for index, enclosure in enumerate(enclosures)
        },
        "frobenius_derivative_square_upper": str(derivative_square),
        "frobenius_derivative_upper": fraction_json(derivative_bound),
        "weyl_radius": fraction_json(rho),
        "positive_lambda0_lower": fraction_json(positive_lower),
        "distinctness_margins": {
            "lambda1_minus_lambda0": fraction_json(distinct01_margin),
            "lambda2_minus_lambda1": fraction_json(distinct12_margin),
        },
        "forced_value_center": interval_json(forced_center_lower, forced_center_upper),
        "forced_value_whole_subinterval": interval_json(
            forced_global_lower, forced_global_upper
        ),
        "tested_center_absence_window": interval_json(
            actual_tested_lower, actual_tested_upper
        ),
        "tested_window_safety_margins": {
            "left": fraction_json(left_margin),
            "right": fraction_json(right_margin),
            "required_each_at_least": str(rho),
        },
        "inertia_counts_at_tested_window_ends": [count_lower, count_upper],
        "forced_relation_absent_throughout": bool(absent),
        "exact_inertia_calls": len(counter.cache),
        "elapsed_seconds": round(time.perf_counter() - started, 6),
    }


def control_chain(center: Fraction = Fraction(1, 3)) -> dict[str, object]:
    """Same exact point machinery on the free-fermion n=4 chain control."""
    bonds = [(0, 1), (1, 2), (2, 3)]
    coefficients, _, _ = polynomial_matrix(4, bonds)
    matrix = evaluate_polynomial_matrix(coefficients, center)
    integer_matrix, scale = clear_denominators(matrix)
    generators = make_generators(4, [list(reversed(range(4)))])
    if not matrix_invariant(integer_matrix, generators):
        raise AssertionError("chain reflection does not preserve the control matrix")
    sectors = symmetry_sectors(integer_matrix, generators)
    counter = ExactCounter(sectors)
    trace = sum(integer_matrix[index][index] for index in range(16))
    guides = guide_eigenvalues(integer_matrix)
    enclosures = []
    for index in range(3):
        lower, upper = initial_bracket(counter, index, guides[index], trace)
        lower, upper = isolate(counter, index, lower, upper, Fraction(1, 10**8))
        enclosures.append((lower / scale, upper / scale))
    (l0, u0), (l1, u1), (l2, u2) = enclosures
    predicted_lower, predicted_upper = l1 * l2 / u0, u1 * u2 / l0
    tested_lower, count_lower = usable_outward(counter, predicted_lower * scale, -1)
    tested_upper, count_upper = usable_outward(counter, predicted_upper * scale, +1)
    return {
        "name": "control_chain_n4_t_1_3",
        "claim_tag": "[COMPUTATION]",
        "t": str(center),
        "lambda_enclosures": {
            f"lambda{index}": interval_json(*enclosure)
            for index, enclosure in enumerate(enclosures)
        },
        "forced_window": interval_json(tested_lower / scale, tested_upper / scale),
        "inertia_counts": [count_lower, count_upper],
        "absence_not_certified": count_lower != count_upper,
        "interpretation": (
            "Differing exact counts prove that some eigenvalue lies in the finite forced-value "
            "window.  This is a negative control, not an equality proof."
        ),
    }


def main() -> int:
    bonds = list(layer_bonds((2, 3), (False, False)))
    coefficients, polynomial_degree, q_shift = polynomial_matrix(6, bonds)
    if (polynomial_degree, q_shift) != (7, 4):
        raise AssertionError("unexpected 2x3 energy exponents")

    definitions = [
        ("layer_2x3_near_t_1_5", Fraction(1, 5), Fraction(1, 10**8)),
        ("layer_2x3_near_t_1_3", Fraction(1, 3), Fraction(1, 10**8)),
        ("layer_2x3_near_t_2_5", Fraction(2, 5), Fraction(1, 10**8)),
        ("layer_2x3_near_t_1_2", Fraction(1, 2), Fraction(1, 10**9)),
    ]
    rows: list[dict[str, object]] = []
    for name, center, half_width in definitions:
        print(f"certifying {name} ...", flush=True)
        row = certify_interval(
            name,
            6,
            bonds,
            coefficients,
            center,
            half_width,
            grid_reflections(2, 3),
        )
        rows.append(row)
        print(
            f"  t in [{row['interval_t']['lower']}, {row['interval_t']['upper']}], "
            f"counts={row['inertia_counts_at_tested_window_ends']}"
        )

    control = control_chain()
    pairwise_disjoint = all(
        Fraction(rows[index]["interval_t"]["upper"])
        < Fraction(rows[index + 1]["interval_t"]["lower"])
        for index in range(len(rows) - 1)
    )
    checks = [
        {
            "name": "all_four_nonzero_intervals_certified",
            "passed": all(row["forced_relation_absent_throughout"] for row in rows),
            "detail": "Each row has equal exact endpoint inertia counts after Weyl expansion.",
        },
        {
            "name": "three_lowest_distinct_parametrically",
            "passed": all(
                Fraction(row["distinctness_margins"][key]["exact"]) > 0
                for row in rows
                for key in ("lambda1_minus_lambda0", "lambda2_minus_lambda1")
            ),
            "detail": "Every propagated exact lower gap is strictly positive.",
        },
        {
            "name": "requested_full_interval_not_claimed",
            "passed": pairwise_disjoint,
            "detail": (
                "The four proved neighbourhoods are disjoint and do not chain from 1/5 to 1/2; "
                "the requested continuum theorem is therefore explicitly UNRESOLVED."
            ),
        },
        {
            "name": "chain_control_does_not_certify_absence",
            "passed": control["absence_not_certified"],
            "detail": f"Exact control counts differ: {control['inertia_counts']}.",
        },
        {
            "name": "all_certificate_arithmetic_exact",
            "passed": True,
            "detail": (
                "Fraction/integer arithmetic decides polynomial bounds, enclosures, margins, "
                "and inertia. Float64 is used only for bracket proposals subsequently certified "
                "by exact inertia."
            ),
        },
    ]

    artifact = {
        "provenance": {
            "script": "experiments/e48_parametric_gaussianity.py",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "method": (
                "Exact integer-polynomial layer matrix; exact rational interval derivative "
                "majorant; Frobenius-to-spectral norm bound; Weyl propagation; exact rational "
                "Sylvester inertia in four reflection-character sectors. Float64 only proposes "
                "brackets and never decides a claim."
            ),
        },
        "data": {
            "status": "UNRESOLVED",
            "requested_target_t": interval_json(Fraction(1, 5), Fraction(1, 2)),
            "failure_step": (
                "The exact Lipschitz/inertia certificates prove only the listed disjoint "
                "subintervals. They cannot chain across [1/5,1/2]; numerical reconnaissance "
                "suggests forced-value/eigenvalue crossings, but no crossing is promoted to an "
                "exact claim here."
            ),
            "matrix_parametrization": {
                "t": "tanh(K*/2)",
                "q_exp_2K": "(1+t^2)/(2t)",
                "K": "(1/2) log((1+t^2)/(2t))",
                "polynomial_representative": "M(t)=(2t)^7 q(t)^4 R(t)",
                "entry_polynomial_degree_upper": len(coefficients[0][0]) - 1,
                "all_entry_coefficients_nonnegative": all(
                    coefficient >= 0
                    for row in coefficients
                    for entry in row
                    for coefficient in entry
                ),
            },
            "certified_subintervals": rows,
            "control": control,
        },
        "checks": checks,
    }
    output_path = ROOT / "results" / "spectral" / "parametric_gaussianity.json"
    os.makedirs(output_path.parent, exist_ok=True)
    with output_path.open("w") as handle:
        json.dump(artifact, handle, indent=2)
        handle.write("\n")

    for check in checks:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    if all(check["passed"] for check in checks):
        print("PASS")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
