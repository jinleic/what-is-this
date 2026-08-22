"""Standalone exact regression for the physical-parity Gaussian obstruction.

This deliberately reimplements the matrix, sector split, rational bisection,
and inertia certificate rather than importing the experiment or reading JSON.
"""

from __future__ import annotations

import signal
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ising.transfer_matrix import layer_bonds  # noqa: E402

FAILURES: list[str] = []

def hard_timeout(_signum, _frame):
    raise TimeoutError("NON-DECISIVE: standalone parity-sector test exceeded 120 seconds")


def check(name: str, passed: bool, detail: str) -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {detail}")
    if not passed:
        FAILURES.append(name)


def build_rational_matrix(n: int, bonds: list[tuple[int, int]], t: Fraction):
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    energies = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energies.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = energies[0] % 2
    weights = [q ** ((energy - parity) // 2) for energy in energies]
    powers = [t**h for h in range(n + 1)]
    kernel = [[powers[bin(i ^ j).count("1")] for j in range(dim)] for i in range(dim)]
    result = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            value = sum((kernel[i][k] * weights[k] * kernel[j][k] for k in range(dim)), Fraction(0))
            result[i][j] = result[j][i] = value
    return result


def sector_blocks(R: list[list[Fraction]]):
    dim = len(R)
    half = dim // 2
    mask = dim - 1
    blocks = {}
    for sign in (1, -1):
        block = [[Fraction(0)] * half for _ in range(half)]
        for a in range(half):
            pa = mask ^ a
            for b in range(half):
                pb = mask ^ b
                block[a][b] = (R[a][b] + sign * R[a][pb] + sign * R[pa][b] + R[pa][pb]) / 2
        blocks[sign] = block
    return blocks


def inertia_below(R, shift: Fraction):
    n = len(R)
    work = [[R[i][j] - (shift if i == j else 0) for j in range(n)] for i in range(n)]
    negative = 0
    for k in range(n):
        pivot = work[k][k]
        if pivot == 0:
            return None
        negative += int(pivot < 0)
        row = work[k]
        for i in range(k + 1, n):
            factor = work[i][k] / pivot
            if factor:
                for j in range(k, n):
                    work[i][j] -= factor * row[j]
    return negative


def isolate(R, index: int, relative: Fraction):
    lo = Fraction(0)
    hi = sum((R[i][i] for i in range(len(R))), Fraction(0)) + 1
    while not (lo > 0 and hi - lo <= relative * lo):
        span = hi - lo
        accepted = None
        for trial_number in range(40):
            trial = (
                lo + span / 2
                if trial_number == 0
                else lo + span / 2 + span * Fraction((-1) ** trial_number, 3 ** (trial_number + 2))
            )
            count = inertia_below(R, trial)
            if count is not None:
                accepted = (trial, count)
                break
        if accepted is None:
            raise RuntimeError("could not find a nondegenerate bisection shift")
        trial, count = accepted
        if count <= index:
            lo = trial
        else:
            hi = trial
    return lo, hi


def assignment_certificate(blocks, even_sign: int, odd_sign: int):
    even_min = isolate(blocks[even_sign], 0, Fraction(1, 10**6))
    odd_min = isolate(blocks[odd_sign], 0, Fraction(1, 10**6))
    odd_second = isolate(blocks[odd_sign], 1, Fraction(1, 10**6))
    relevant = (even_min, odd_min, odd_second)
    pairwise_distinct = all(
        left[1] < right[0] or right[1] < left[0]
        for i, left in enumerate(relevant)
        for right in relevant[i + 1 :]
    )
    gaussian_ordered = even_min[1] < odd_min[0] and odd_min[1] < odd_second[0]
    lo = odd_min[0] * odd_second[0] / even_min[1]
    hi = odd_min[1] * odd_second[1] / even_min[0]
    count_lo = inertia_below(blocks[even_sign], lo)
    count_hi = inertia_below(blocks[even_sign], hi)
    if count_lo is None or count_hi is None:
        raise RuntimeError("unexpected endpoint eigenvalue; test requires explicit outward handling")
    return pairwise_distinct, gaussian_ordered, count_lo, count_hi


signal.signal(signal.SIGALRM, hard_timeout)
signal.alarm(120)

t = Fraction(1, 3)
bonds = list(layer_bonds((2, 3), (False, False)))
R = build_rational_matrix(6, bonds, t)
flip = [63 ^ state for state in range(64)]
signed_permutation = (
    sorted(flip) == list(range(64))
    and all(flip[flip[state]] == state for state in range(64))
    and all(flip[state] != state for state in range(64))
)
commutes = all(R[flip[i]][flip[j]] == R[i][j] for i in range(64) for j in range(64))
blocks = sector_blocks(R)
check(
    "physical P is an exact commuting permutation and sectors are 32x32 rational matrices",
    signed_permutation and commutes and all(len(blocks[s]) == 32 for s in (1, -1)),
    f"signed permutation: {signed_permutation}; PR=RP: {commutes}; dimensions {[len(blocks[s]) for s in (1, -1)]}",
)

# Required end-to-end assignment: P=+ as fermionic even and P=- as fermionic odd.
distinct, ordered, count_lo, count_hi = assignment_certificate(blocks, 1, -1)
check(
    "P+=even/P-=odd factor hypotheses follow from ordered disjoint enclosures",
    distinct and ordered,
    f"three relevant values pairwise disjoint={distinct}; even_min<odd_min<odd_second={ordered}",
)
check(
    "P+=even/P-=odd forced cross-sector value is absent",
    count_lo == count_hi,
    f"exact target-sector inertia counts {count_lo} and {count_hi}",
)

# Assignment symmetry is also exercised, although only one assignment was required end-to-end.
swapped_distinct, swapped_ordered, swapped_lo, swapped_hi = assignment_certificate(blocks, -1, 1)
check(
    "P-=even/P+=odd assignment is rejected and its formal forced window is empty",
    swapped_distinct and not swapped_ordered and swapped_lo == swapped_hi,
    f"pairwise distinct={swapped_distinct}; required order={swapped_ordered}; exact target-sector inertia counts {swapped_lo} and {swapped_hi}",
)

if FAILURES:
    print(f"FAIL: {len(FAILURES)} checks failed: {FAILURES}")
    raise SystemExit(1)
print("PASS: independently rederived exact 2x3 physical-parity certificate for both assignments.")
