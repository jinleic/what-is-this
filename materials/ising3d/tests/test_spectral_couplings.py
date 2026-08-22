"""Standalone exact checks for the extended spectral-Gaussianity certificates.

The 2x3, t=1/5 certificate is rebuilt here from scratch: this file does not import the
experiment.  It constructs the dense rational R directly, clears denominators only
after construction, encloses the three lowest eigenvalues by exact full-matrix inertia
bisection, and tests the forced-value window.  Two slower Fraction-LDL counts provide an
independent check of the fraction-free final counts.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


# ---------------------------------------------------------------- independent dense construction

def build_dense_R(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Direct definition R = P_t diag(q^e) P_t, independent of e39's fast transform."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    energies: list[int] = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = energies[0] % 2
    assert all((energy - parity) % 2 == 0 for energy in energies)
    diagonal = [q ** ((energy - parity) // 2) for energy in energies]
    powers = [t**h for h in range(n + 1)]
    p_matrix = [
        [powers[(row ^ col).bit_count()] for col in range(dim)] for row in range(dim)
    ]
    result = [[Fraction(0)] * dim for _ in range(dim)]
    for row in range(dim):
        for col in range(row, dim):
            value = sum(
                (p_matrix[row][middle] * diagonal[middle] * p_matrix[middle][col]
                 for middle in range(dim)),
                Fraction(0),
            )
            result[row][col] = result[col][row] = value
    return result, q


def clear_matrix_denominators(matrix: list[list[Fraction]]):
    scale = 1
    for row in matrix:
        for value in row:
            scale = math.lcm(scale, value.denominator)
    integer_matrix = [
        [(value * scale).numerator for value in row] for row in matrix
    ]
    return integer_matrix, scale


# --------------------------------------------------------------- independent exact inertia

def bareiss_inertia(integer_matrix: list[list[int]], shift: Fraction) -> int | None:
    """Full-matrix exact fraction-free LDL; no symmetry reduction is used in this test."""
    size = len(integer_matrix)
    lower = [[0] * size for _ in range(size)]
    common = 0
    for row in range(size):
        for col in range(row + 1):
            value = shift.denominator * integer_matrix[row][col]
            if row == col:
                value -= shift.numerator
            lower[row][col] = value
            common = math.gcd(common, abs(value))
    if common > 1:
        for row in range(size):
            for col in range(row + 1):
                lower[row][col] //= common

    previous = 1
    negatives = 0
    for pivot_index in range(size):
        pivot = lower[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        if pivot_index + 1 == size:
            break
        for row_index in range(pivot_index + 1, size):
            row = lower[row_index]
            left = row[pivot_index]
            for col_index in range(pivot_index + 1, row_index + 1):
                numerator = (
                    pivot * row[col_index]
                    - left * lower[col_index][pivot_index]
                )
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                row[col_index] = quotient
        previous = pivot
    return negatives


def fraction_ldl_inertia(matrix: list[list[Fraction]], shift: Fraction) -> int | None:
    """Slow exact Fraction implementation matching e38, used at the two decisive endpoints."""
    size = len(matrix)
    work = [
        [matrix[row][col] - (shift if row == col else 0) for col in range(size)]
        for row in range(size)
    ]
    negatives = 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot < 0:
            negatives += 1
        pivot_inverse = 1 / pivot
        pivot_row = work[pivot_index]
        for row_index in range(pivot_index + 1, size):
            factor = work[row_index][pivot_index] * pivot_inverse
            if factor:
                row = work[row_index]
                for col_index in range(pivot_index, size):
                    row[col_index] -= factor * pivot_row[col_index]
    return negatives


def usable_inside(
    integer_matrix: list[list[int]], lo: Fraction, hi: Fraction
) -> tuple[Fraction, int]:
    span = hi - lo
    for attempt in range(40):
        trial = (
            lo + span / 2
            if attempt == 0
            else lo + span / 2 + span * Fraction((-1) ** attempt, 3 ** (attempt + 2))
        )
        if lo < trial < hi:
            count = bareiss_inertia(integer_matrix, trial)
            if count is not None:
                return trial, count
    raise RuntimeError("no usable exact shift")


def usable_outward(
    integer_matrix: list[list[int]], shift: Fraction, direction: int
) -> tuple[Fraction, int]:
    count = bareiss_inertia(integer_matrix, shift)
    if count is not None:
        return shift, count
    step = abs(shift) if shift else Fraction(1)
    for attempt in range(1, 40):
        moved = shift + direction * step * Fraction(1, 10 ** (20 - attempt // 2))
        count = bareiss_inertia(integer_matrix, moved)
        if count is not None:
            return moved, count
    raise RuntimeError("no usable outward shift")


def isolate(
    integer_matrix: list[list[int]],
    index: int,
    lo: Fraction,
    hi: Fraction,
    relative_width: Fraction,
) -> tuple[Fraction, Fraction]:
    while not (lo > 0 and hi - lo <= relative_width * lo):
        trial, count = usable_inside(integer_matrix, lo, hi)
        if count <= index:
            lo = trial
        else:
            hi = trial
    return lo, hi


def exact_bracket_from_guide(
    integer_matrix: list[list[int]], index: int, guide: float
) -> tuple[Fraction, Fraction]:
    """Float chooses a proposal only; exact counts establish both bracket invariants."""
    center = Fraction(format(guide, ".15e"))
    lo, hi = center * Fraction(49, 50), center * Fraction(51, 50)
    lo, count_lo = usable_outward(integer_matrix, lo, -1)
    while count_lo > index:
        lo, count_lo = usable_outward(integer_matrix, lo / 2, -1)
    hi, count_hi = usable_outward(integer_matrix, hi, +1)
    while count_hi <= index:
        hi, count_hi = usable_outward(integer_matrix, hi * 2, +1)
    return lo, hi


def rederive_2x3_t_one_fifth():
    bonds = list(layer_bonds((2, 3), (False, False)))
    rational_matrix, q = build_dense_R(6, bonds, Fraction(1, 5))
    integer_matrix, scale = clear_matrix_denominators(rational_matrix)

    largest = max(abs(value) for row in integer_matrix for value in row)
    floating = np.array(
        [[value / largest for value in row] for row in integer_matrix], dtype=float
    )
    guides = [float(value * largest) for value in np.linalg.eigvalsh(floating)[:3]]

    enclosures = []
    for index in range(3):
        lo, hi = exact_bracket_from_guide(integer_matrix, index, guides[index])
        enclosures.append(isolate(integer_matrix, index, lo, hi, Fraction(1, 100000)))

    (l0, h0), (l1, h1), (l2, h2) = enclosures
    predicted_lo = l1 * l2 / h0
    predicted_hi = h1 * h2 / l0
    tested_lo, count_lo = usable_outward(integer_matrix, predicted_lo, -1)
    tested_hi, count_hi = usable_outward(integer_matrix, predicted_hi, +1)

    # Independently repeat the two decisive counts on the original rational R.
    fraction_count_lo = fraction_ldl_inertia(rational_matrix, tested_lo / scale)
    fraction_count_hi = fraction_ldl_inertia(rational_matrix, tested_hi / scale)
    return {
        "q": q,
        "positive_count": bareiss_inertia(integer_matrix, Fraction(0)),
        "enclosures": enclosures,
        "distinct": h0 < l1 and h1 < l2,
        "tested_window_contains_prediction": tested_lo <= predicted_lo <= predicted_hi <= tested_hi,
        "counts": (count_lo, count_hi),
        "fraction_counts": (fraction_count_lo, fraction_count_hi),
    }


print("1. independent end-to-end exact certificate: 2x3 at t=1/5")
independent = rederive_2x3_t_one_fifth()
check("exp(2K) is 13/5", independent["q"] == Fraction(13, 5))
check("R is positive definite", independent["positive_count"] == 0)
check("three lowest eigenvalues have disjoint exact enclosures", independent["distinct"])
check(
    "tested window encloses the forced-value window",
    independent["tested_window_contains_prediction"],
)
check(
    "forced-value window is exactly empty",
    independent["counts"][0] == independent["counts"][1] == 4,
    f"fraction-free exact counts={independent['counts']}",
)
check(
    "independent Fraction LDL confirms decisive counts",
    independent["fraction_counts"] == independent["counts"],
    f"Fraction counts={independent['fraction_counts']}",
)


# -------------------------------------------------------------- stored-result consistency
print("\n2. exact internal consistency of results/spectral/couplings.json")
result_path = ROOT / "results" / "spectral" / "couplings.json"
with result_path.open() as handle:
    stored = json.load(handle)

check("top-level schema", set(stored) == {"meta", "data", "checks"})
rows = stored["data"]["rows"]
check("all five requested cases are present", len(rows) == 5)
expected_names = {
    "layer_2x3_t_1_5",
    "layer_2x3_t_2_5",
    "layer_2x3_t_1_2",
    "layer_2x4_t_1_3",
    "control_chain_n8_t_1_3",
}
check("case names match the requested family", {row["name"] for row in rows} == expected_names)

for row in rows:
    name = row["name"]
    required = {
        "name",
        "t",
        "exp_2K",
        "dim",
        "lambda_enclosures",
        "three_lowest_distinct",
        "relative_gaps",
        "window",
        "inertia_counts_at_window_ends",
        "predicted_eigenvalue_absent",
        "elapsed",
    }
    check(f"{name}: required fields", required <= set(row))
    if row.get("status") != "CERTIFIED":
        check(f"{name}: unresolved resource record is explicit", row.get("status") == "UNRESOLVED-RESOURCE")
        continue

    t = Fraction(row["t"])
    check(f"{name}: exp(2K) identity", Fraction(row["exp_2K"]) == (1 + t * t) / (2 * t))
    check(f"{name}: dimension identity", row["dim"] == 1 << row["n_sites"])

    intervals: list[tuple[Fraction, Fraction]] = []
    endpoint_counts_ok = True
    widths_ok = True
    for index in range(3):
        enclosure = row["lambda_enclosures"][f"lambda{index}"]
        lo, hi = Fraction(enclosure["lower"]), Fraction(enclosure["upper"])
        intervals.append((lo, hi))
        widths_ok &= 0 < lo < hi and Fraction(enclosure["relative_width"]) == (hi - lo) / lo
        endpoint_counts_ok &= (
            enclosure["count_below_lower"] <= index
            and enclosure["count_below_upper"] > index
        )
    check(f"{name}: exact lambda interval widths", widths_ok)
    check(f"{name}: stored endpoint-count invariants", endpoint_counts_ok)

    (l0, h0), (l1, h1), (l2, h2) = intervals
    distinct = h0 < l1 and h1 < l2
    check(f"{name}: distinctness boolean is exact", row["three_lowest_distinct"] == distinct)
    gap01 = (l1 - h0) / l1
    gap12 = (l2 - h1) / l2
    gaps = row["relative_gaps"]
    check(
        f"{name}: exact relative gap lower bounds",
        Fraction(gaps["gap01_lower_bound"]) == gap01
        and Fraction(gaps["gap12_lower_bound"]) == gap12,
    )

    predicted = row["window"]["predicted"]
    predicted_lo, predicted_hi = Fraction(predicted["lower"]), Fraction(predicted["upper"])
    check(
        f"{name}: forced-value interval arithmetic",
        predicted_lo == l1 * l2 / h0 and predicted_hi == h1 * h2 / l0,
    )
    tested = row["window"]["actually_tested"]
    tested_lo, tested_hi = Fraction(tested["lower"]), Fraction(tested["upper"])
    encloses = tested_lo <= predicted_lo <= predicted_hi <= tested_hi
    check(
        f"{name}: outward tested window",
        encloses == row["window"]["tested_window_encloses_prediction"],
    )

    count_lo, count_hi = row["inertia_counts_at_window_ends"]
    absent = count_lo == count_hi
    check(f"{name}: absence boolean matches exact counts", absent == row["predicted_eigenvalue_absent"])
    if row["kind"] == "3D-layer":
        check(f"{name}: certified non-Gaussian obstruction", distinct and absent)
    else:
        check(f"{name}: occupied-window consistency only", distinct and not absent and count_lo != count_hi)

    sector_dimensions = row["symmetry"]["sector_dimensions"]
    check(
        f"{name}: exact sector dimensions cover the matrix",
        row["symmetry"]["matrix_invariant"]
        and sum(sector_dimensions) == row["symmetry"]["dimensions_sum"] == row["dim"],
    )
    check(
        f"{name}: elapsed timing recorded",
        row["elapsed"]["total_seconds"] >= 0,
    )

check("every stored top-level check passed", all(item["passed"] for item in stored["checks"]))
check("one top-level check per row", {item["name"] for item in stored["checks"]} == expected_names)

print()
if FAILURES:
    print(f"FAIL: {len(FAILURES)} checks failed: {FAILURES}")
    raise SystemExit(1)
print("PASS: independent t=1/5 obstruction and all stored exact-certificate invariants verified")
