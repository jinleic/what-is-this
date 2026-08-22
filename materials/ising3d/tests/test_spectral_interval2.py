"""Clean-room standalone verifier for e81's spectral-interval certificate.

This file intentionally does not import experiments.e81_spectral_interval2.
It independently builds the 2x3 matrix from R=P_t diag(q^e) P_t, reconstructs
character-sector inertia counts for the c-jump witnesses, and uses a full
64x64 exact inertia calculation to recheck the widest safe interval.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "interval2.json"

import sys

sys.path.insert(0, str(ROOT))
from ising.clifford import pauli_x, symplectic_form, zz  # noqa: E402
from ising.transfer_matrix import layer_bonds  # noqa: E402


# ---------------------------------------------------------------------------
# Independent direct construction of R(t) and its polynomial multiple M(t).
# ---------------------------------------------------------------------------


def layer_exponents(n: int, bonds: list[tuple[int, int]]) -> tuple[list[int], int, int]:
    energies: list[int] = []
    for state in range(1 << n):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    parity = energies[0] % 2
    exponents = [(energy - parity) // 2 for energy in energies]
    shift = -min(exponents)
    degree = max(exponent + shift for exponent in exponents)
    return exponents, degree, shift


def direct_m_matrix(t: Fraction) -> list[list[Fraction]]:
    """Build M=(2t)^7 q^4 R directly, rather than through coefficient lists."""
    n = 6
    bonds = list(layer_bonds((2, 3), (False, False)))
    exponents, degree, shift = layer_exponents(n, bonds)
    assert (degree, shift) == (7, 4)
    q = (1 + t * t) / (2 * t)
    scalar = (2 * t) ** degree * q**shift
    dimension = 1 << n
    return [
        [
            scalar * sum(
                (
                    t ** ((row ^ middle).bit_count() + (middle ^ col).bit_count())
                    * q**exponents[middle]
                    for middle in range(dimension)
                ),
                Fraction(0),
            )
            for col in range(dimension)
        ]
        for row in range(dimension)
    ]


def independent_polynomial_matrix() -> list[list[list[int]]]:
    """A second construction, used only for the safe-interval derivative bound."""
    n = 6
    bonds = list(layer_bonds((2, 3), (False, False)))
    exponents, degree, shift = layer_exponents(n, bonds)
    dimension = 1 << n
    maximum_degree = 2 * n + 2 * degree
    coefficients = [
        [[0] * (maximum_degree + 1) for _ in range(dimension)]
        for _ in range(dimension)
    ]
    for row in range(dimension):
        for col in range(row, dimension):
            poly = [0] * (maximum_degree + 1)
            for middle, exponent in enumerate(exponents):
                power = exponent + shift
                distance = (row ^ middle).bit_count() + (middle ^ col).bit_count()
                for even in range(power + 1):
                    poly[degree - power + distance + 2 * even] += (
                        2 ** (degree - power) * math.comb(power, even)
                    )
            coefficients[row][col] = poly
            coefficients[col][row] = poly
    return coefficients


def evaluate_coefficients(coefficients: list[list[list[int]]], t: Fraction) -> list[list[Fraction]]:
    powers = [t**index for index in range(len(coefficients[0][0]))]
    return [
        [
            sum((Fraction(value) * powers[index] for index, value in enumerate(poly) if value), Fraction(0))
            for poly in row
        ]
        for row in coefficients
    ]


def clear_denominators(matrix: list[list[Fraction]]) -> tuple[list[list[int]], int]:
    scale = 1
    for row in matrix:
        for value in row:
            scale = math.lcm(scale, value.denominator)
    return [[(value * scale).numerator for value in row] for row in matrix], scale


# ---------------------------------------------------------------------------
# Independent character-sector inertia implementation.
# ---------------------------------------------------------------------------


def state_permutation(n: int, site_permutation: list[int]) -> list[int]:
    answer: list[int] = []
    for state in range(1 << n):
        moved = 0
        for old_site, new_site in enumerate(site_permutation):
            bit = (state >> (n - 1 - old_site)) & 1
            moved |= bit << (n - 1 - new_site)
        answer.append(moved)
    return answer


def independent_bases() -> list[tuple[list[list[tuple[int, int]]], list[int]]]:
    """Return character bases as sparse rows, independently of the producer."""
    n = 6
    row_reflection = state_permutation(
        n, [(1 - row) * 3 + col for row in range(2) for col in range(3)]
    )
    col_reflection = state_permutation(
        n, [row * 3 + (2 - col) for row in range(2) for col in range(3)]
    )
    spin_flip = [state ^ 63 for state in range(64)]
    generators = [row_reflection, col_reflection, spin_flip]
    group: list[list[int]] = []
    for mask in range(8):
        permutation = list(range(64))
        for generator_index, generator in enumerate(generators):
            if mask & (1 << generator_index):
                permutation = [generator[state] for state in permutation]
        group.append(permutation)

    seen: set[int] = set()
    representatives: list[int] = []
    for state in range(64):
        if state not in seen:
            orbit = {permutation[state] for permutation in group}
            representatives.append(min(orbit))
            seen.update(orbit)

    bases: list[tuple[list[list[tuple[int, int]]], list[int]]] = []
    for character in range(8):
        vectors: list[list[tuple[int, int]]] = []
        gram: list[int] = []
        for representative in representatives:
            coefficients: dict[int, int] = {}
            for group_element, permutation in enumerate(group):
                sign = -1 if (character & group_element).bit_count() % 2 else 1
                image = permutation[representative]
                coefficients[image] = coefficients.get(image, 0) + sign
            sparse = sorted((state, value) for state, value in coefficients.items() if value)
            if sparse:
                vectors.append(sparse)
                gram.append(sum(value * value for _, value in sparse))
        bases.append((vectors, gram))
    assert [len(gram) for _, gram in bases] == [14, 6, 6, 6, 10, 10, 6, 6]
    return bases


def sector_blocks(
    matrix: list[list[int]], bases: list[tuple[list[list[tuple[int, int]]], list[int]]]
) -> list[tuple[list[list[int]], list[int]]]:
    blocks: list[tuple[list[list[int]], list[int]]] = []
    for vectors, gram in bases:
        size = len(gram)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            for col in range(row + 1):
                value = sum(
                    left_coefficient * matrix[left_state][right_state] * right_coefficient
                    for left_state, left_coefficient in vectors[row]
                    for right_state, right_coefficient in vectors[col]
                )
                block[row][col] = block[col][row] = value
        blocks.append((block, gram))
    return blocks


def inertia_count_lower_triangular(lower: list[list[int]]) -> int | None:
    """Independent fraction-free LDL negative-pivot count."""
    work = [row[:] for row in lower]
    previous = 1
    negatives = 0
    for pivot_index in range(len(work)):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        for row_index in range(pivot_index + 1, len(work)):
            left = work[row_index][pivot_index]
            for col_index in range(pivot_index + 1, row_index + 1):
                numerator = pivot * work[row_index][col_index] - left * work[col_index][pivot_index]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                work[row_index][col_index] = quotient
        previous = pivot
    return negatives


class SectorCounter:
    def __init__(self, blocks: list[tuple[list[list[int]], list[int]]], scale: int):
        self.blocks = blocks
        self.scale = scale
        self.cache: dict[Fraction, int | None] = {}

    def count_below(self, level: Fraction) -> int | None:
        level = Fraction(level)
        if level in self.cache:
            return self.cache[level]
        shift = level * self.scale
        result = 0
        for block, gram in self.blocks:
            size = len(gram)
            shifted = [[0] * size for _ in range(size)]
            common = 0
            for row in range(size):
                for col in range(row + 1):
                    value = shift.denominator * block[row][col]
                    if row == col:
                        value -= shift.numerator * gram[row]
                    shifted[row][col] = value
                    common = math.gcd(common, abs(value))
            if common > 1:
                for row in range(size):
                    for col in range(row + 1):
                        shifted[row][col] //= common
            count = inertia_count_lower_triangular(shifted)
            if count is None:
                self.cache[level] = None
                return None
            result += count
        self.cache[level] = result
        return result


def full_count_below(matrix: list[list[int]], scale: int, level: Fraction) -> int | None:
    """No-sector full 64x64 exact check used for the widest safe interval."""
    shift = Fraction(level) * scale
    dimension = len(matrix)
    shifted = [
        [
            shift.denominator * matrix[row][col] - (shift.numerator if row == col else 0)
            for col in range(dimension)
        ]
        for row in range(dimension)
    ]
    return inertia_count_lower_triangular(shifted)


def float_guides(matrix: list[list[int]], scale: int) -> list[float]:
    largest = max(abs(value) for row in matrix for value in row)
    values = np.linalg.eigvalsh(
        np.array([[value / largest for value in row] for row in matrix], dtype=float)
    )
    multiplier = float(Fraction(largest, scale))
    return [float(value * multiplier) for value in values[:3]]


def outward(counter: SectorCounter, level: Fraction, direction: int) -> tuple[Fraction, int]:
    count = counter.count_below(level)
    if count is not None:
        return level, count
    magnitude = max(abs(level), Fraction(1))
    for exponent in range(30, 76):
        candidate = level + direction * magnitude / (10**exponent)
        count = counter.count_below(candidate)
        if count is not None:
            return candidate, count
    raise AssertionError("could not move an exact endpoint away from a spectral root")


def inside(counter: SectorCounter, lower: Fraction, upper: Fraction) -> tuple[Fraction, int]:
    width = upper - lower
    for attempt in range(48):
        candidate = (
            lower + width / 2
            if attempt == 0
            else lower + width / 2 + width * Fraction((-1) ** attempt, 3 ** (attempt + 2))
        )
        count = counter.count_below(candidate)
        if lower < candidate < upper and count is not None:
            return candidate, count
    raise AssertionError("could not find a nonsingular interior shift")


def isolate_low_three(counter: SectorCounter, guides: list[float], relative: Fraction) -> list[tuple[Fraction, Fraction]]:
    enclosures: list[tuple[Fraction, Fraction]] = []
    for index, guide in enumerate(guides):
        proposal = Fraction(format(guide, ".15e"))
        lower = max(Fraction(0), proposal * Fraction(99, 100))
        upper = proposal * Fraction(101, 100)
        lower, lower_count = outward(counter, lower, -1)
        while lower_count > index:
            lower, lower_count = outward(counter, lower / 2, -1)
        upper, upper_count = outward(counter, upper, +1)
        while upper_count <= index:
            upper, upper_count = outward(counter, upper * 2, +1)
        while upper - lower > relative * lower:
            trial, count = inside(counter, lower, upper)
            if count <= index:
                lower = trial
            else:
                upper = trial
        assert counter.count_below(lower) is not None and counter.count_below(lower) <= index
        assert counter.count_below(upper) is not None and counter.count_below(upper) > index
        enclosures.append((lower, upper))
    return enclosures


def independently_certify_c(t: Fraction, expected: int) -> int:
    matrix = direct_m_matrix(t)
    integer_matrix, scale = clear_denominators(matrix)
    counter = SectorCounter(sector_blocks(integer_matrix, independent_bases()), scale)
    enclosures = isolate_low_three(counter, float_guides(integer_matrix, scale), Fraction(1, 1_000_000))
    (l0, u0), (l1, u1), (l2, u2) = enclosures
    forced_lower = l1 * l2 / u0
    forced_upper = u1 * u2 / l0
    _, lower_count = outward(counter, forced_lower, -1)
    _, upper_count = outward(counter, forced_upper, +1)
    assert lower_count == upper_count, (t, forced_lower, forced_upper, lower_count, upper_count)
    assert lower_count == expected, (t, lower_count, expected)
    return lower_count


# ---------------------------------------------------------------------------
# Safe interval and chain-control checks
# ---------------------------------------------------------------------------


def sqrt_upper(value: Fraction, denominator: int = 10**18) -> Fraction:
    target = (value.numerator * denominator * denominator + value.denominator - 1) // value.denominator
    numerator = math.isqrt(target)
    if numerator * numerator < target:
        numerator += 1
    answer = Fraction(numerator, denominator)
    assert answer * answer >= value
    return answer


def independent_derivative_bound(coefficients: list[list[list[int]],], upper_t: Fraction) -> Fraction:
    frobenius_square = Fraction(0)
    row_sums: list[Fraction] = []
    for row in coefficients:
        total = Fraction(0)
        for poly in row:
            derivative = sum(
                (Fraction(power * coefficient) * upper_t ** (power - 1)
                 for power, coefficient in enumerate(poly) if power and coefficient),
                Fraction(0),
            )
            frobenius_square += derivative * derivative
            total += derivative
        row_sums.append(total)
    return min(sqrt_upper(frobenius_square), max(row_sums))


def verify_widest_safe_interval(data: dict[str, object]) -> None:
    rows = data["safe_closed_intervals"]
    widest = max(rows, key=lambda row: Fraction(row["interval_t"]["width"]))
    assert widest["name"] == "safe_interval_t_1_4"
    center = Fraction(widest["center_t"])
    half_width = Fraction(widest["half_width_t"])
    lower_t = Fraction(widest["interval_t"]["lower"])
    upper_t = Fraction(widest["interval_t"]["upper"])
    assert (lower_t + upper_t) / 2 == center
    assert (upper_t - lower_t) / 2 == half_width

    direct = direct_m_matrix(center)
    coefficients = independent_polynomial_matrix()
    assert evaluate_coefficients(coefficients, center) == direct
    integer_matrix, scale = clear_denominators(direct)
    counter = SectorCounter(sector_blocks(integer_matrix, independent_bases()), scale)

    stored_enclosures = widest["lambda_center_enclosures"]
    enclosures: list[tuple[Fraction, Fraction]] = []
    for index in range(3):
        stored = stored_enclosures[f"lambda{index}"]
        lower, upper = Fraction(stored["lower"]), Fraction(stored["upper"])
        lower_count = counter.count_below(lower)
        upper_count = counter.count_below(upper)
        assert lower_count is not None and lower_count <= index
        assert upper_count is not None and upper_count > index
        enclosures.append((lower, upper))

    derivative = independent_derivative_bound(coefficients, upper_t)
    assert derivative == Fraction(widest["selected_spectral_derivative_upper"]["exact"])
    rho = derivative * half_width
    assert rho == Fraction(widest["weyl_radius"]["exact"])
    (l0, u0), (l1, u1), (l2, u2) = enclosures
    assert l0 > rho
    forced_lower = (l1 - rho) * (l2 - rho) / (u0 + rho)
    forced_upper = (u1 + rho) * (u2 + rho) / (l0 - rho)
    stored_forced = widest["forced_value_whole_interval"]
    assert forced_lower == Fraction(stored_forced["lower"])
    assert forced_upper == Fraction(stored_forced["upper"])

    tested = widest["tested_center_absence_window"]
    low, high = Fraction(tested["lower"]), Fraction(tested["upper"])
    # This deliberately bypasses all symmetry sectors.
    low_count = full_count_below(integer_matrix, scale, low)
    high_count = full_count_below(integer_matrix, scale, high)
    assert [low_count, high_count] == widest["inertia_counts_at_tested_window_ends"] == [4, 4]
    assert widest["forced_relation_absent_throughout"] is True


def verify_majorana_control(data: dict[str, object]) -> None:
    control = data["chain_control"]["structural_all_interval"]
    n = 4
    majoranas: list[int] = []
    for site in range(n):
        prefix = sum(1 << previous for previous in range(site))
        majoranas.extend([
            prefix | (1 << (n + site)),
            prefix | (1 << site) | (1 << (n + site)),
        ])
    assert all(
        symplectic_form(majoranas[left], majoranas[right], n) == 1
        for left in range(2 * n)
        for right in range(left + 1, 2 * n)
    )
    assert all(
        majoranas[2 * site] ^ majoranas[2 * site + 1] == pauli_x(n, site)
        for site in range(n)
    )
    assert all(
        majoranas[2 * site + 1] ^ majoranas[2 * site + 2] == zz(n, site, site + 1)
        for site in range(n - 1)
    )
    assert all(control["majorana_bilinear_checks"].values())
    explicit = control["explicit_two_mode_identity"]
    assert explicit["equals_subset_product_eigenvalue"] is True
    assert explicit["forced_value"] == explicit["spectrum_at_t_1_3"][3]


# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------


def main() -> None:
    assert ARTIFACT.is_file(), "e81 producer has not created its certificate artifact"
    artifact = json.loads(ARTIFACT.read_text())
    assert set(artifact) == {"provenance", "data", "checks"}
    assert artifact["provenance"]["script"] == "experiments/e81_spectral_interval2.py"
    data = artifact["data"]
    assert data["status"] == "FALSIFIED"
    assert all(check["passed"] for check in artifact["checks"])
    assert data["matrix_parametrization"]["dimension"] == 64
    assert data["matrix_parametrization"]["n_bonds"] == 7
    assert data["matrix_parametrization"]["sector_dimensions_spinflip_reflections"] == [14, 6, 6, 6, 10, 10, 6, 6]

    # Four independent exact coarse c-jumps are sufficient to prove four
    # distinct forced-value crossings.  Fine brackets are producer refinements.
    crossings = data["crossings"]
    assert len(crossings) == 4
    expected = [(4, 5), (5, 6), (6, 7), (7, 8)]
    previous_upper = Fraction(0)
    for crossing, expected_counts in zip(crossings, expected, strict=True):
        coarse = crossing["coarse_bracket_t"]
        lower, upper = Fraction(coarse["lower"]), Fraction(coarse["upper"])
        assert previous_upper < lower
        assert independently_certify_c(lower, expected_counts[0]) == expected_counts[0]
        assert independently_certify_c(upper, expected_counts[1]) == expected_counts[1]
        assert crossing["c_counts_left_to_right"] == list(expected_counts)
        fine = crossing["fine_bracket_t"]
        assert Fraction(fine["width"]) <= Fraction(1, 10**10)
        previous_upper = upper

    verify_widest_safe_interval(data)
    verify_majorana_control(data)
    maximum_wave8_width = Fraction(1, 50_000_000)
    assert max(Fraction(row["interval_t"]["width"]) for row in data["safe_closed_intervals"]) > maximum_wave8_width
    print("OK")


if __name__ == "__main__":
    main()
