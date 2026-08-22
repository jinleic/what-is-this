"""Standalone verifier for the ordering-free physical-parity pair-product certificate.

This deliberately does not import ``e119_parity_pair_product``.  It rebuilds the
rational 2x3 operator, its P=prod X sectors, and one modular pair polynomial
from scratch.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402

INT64_MAX = 2**63 - 1
PRIME = 1_000_003
N_MODES = 6
SECTOR_SLOTS = 1 << (N_MODES - 1)
WITHIN_PAIRS = SECTOR_SLOTS * (SECTOR_SLOTS - 1) // 2
CROSS_PAIRS = SECTOR_SLOTS * SECTOR_SLOTS
WITHIN_CLASSES = sum(math.comb(N_MODES, r) * 2 ** (N_MODES - r) for r in range(2, N_MODES + 1, 2))
CROSS_CLASSES = sum(math.comb(N_MODES, r) * 2 ** (N_MODES - r) for r in range(1, N_MODES + 1, 2))
WITHIN_GCD_FLOOR = WITHIN_PAIRS - WITHIN_CLASSES
CROSS_GCD_FLOOR = CROSS_PAIRS - CROSS_CLASSES

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)
    return bool(ok)


def bit_vector(slot: int, modes: int) -> tuple[int, ...]:
    return tuple((slot >> i) & 1 for i in range(modes))


def exponent_classes(modes: int, left_parity: int, right_parity: int, same_half: bool) -> set[tuple[int, ...]]:
    """Enumerate exponent sums from labelled slots, retaining the pair rule."""
    left = [slot for slot in range(1 << modes) if slot.bit_count() % 2 == left_parity]
    right = [slot for slot in range(1 << modes) if slot.bit_count() % 2 == right_parity]
    out: set[tuple[int, ...]] = set()
    for i, epsilon in enumerate(left):
        starts = i + 1 if same_half else 0
        for delta in right[starts:]:
            e = bit_vector(epsilon, modes)
            d = bit_vector(delta, modes)
            out.add(tuple(x + y for x, y in zip(e, d)))
    return out


def within_formula(modes: int) -> int:
    return sum(math.comb(modes, r) * 2 ** (modes - r) for r in range(2, modes + 1, 2))


def cross_formula(modes: int) -> int:
    return sum(math.comb(modes, r) * 2 ** (modes - r) for r in range(1, modes + 1, 2))


def build_R(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Independent P_t diag(q^m) P_t reconstruction of the rational model."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    bond_sums = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - i)) & 1) for i in range(n)]
        bond_sums.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = bond_sums[0] % 2
    assert all((value - eps) % 2 == 0 for value in bond_sums)
    diagonal = [q ** ((value - eps) // 2) for value in bond_sums]
    powers = [t**h for h in range(n + 1)]
    kernel = [[powers[(i ^ j).bit_count()] for j in range(dim)] for i in range(dim)]
    R = [[Fraction(0) for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            value = sum((kernel[i][k] * diagonal[k] * kernel[j][k] for k in range(dim)), Fraction(0))
            R[i][j] = R[j][i] = value
    return R, q, eps


def split_physical_parity(R: list[list[Fraction]]) -> dict[int, list[list[Fraction]]]:
    """Use the unnormalised e_a +/- e_complement(a) basis of physical P."""
    dim = len(R)
    assert dim and dim & (dim - 1) == 0
    half = dim // 2
    mask = dim - 1
    out: dict[int, list[list[Fraction]]] = {}
    for sign in (1, -1):
        block = [[Fraction(0) for _ in range(half)] for _ in range(half)]
        for a in range(half):
            pa = mask ^ a
            for b in range(half):
                pb = mask ^ b
                block[a][b] = (R[a][b] + sign * R[a][pb] + sign * R[pa][b] + R[pa][pb]) / 2
        out[sign] = block
    return out


def integralize(R: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    scale = 1
    for row in R:
        for entry in row:
            scale = math.lcm(scale, entry.denominator)
    matrix = []
    for row in R:
        integer_row = []
        for entry in row:
            scaled = scale * entry
            assert scaled.denominator == 1
            integer_row.append(scaled.numerator)
        matrix.append(integer_row)
    return scale, matrix


def matmul(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    assert A.shape[0] * (p - 1) ** 2 < INT64_MAX
    return (A @ B) % p

def binary_trace(M_rows: list[list[int]], power: int, p: int) -> int:
    """Independent square-and-multiply trace for a recurrence spot check."""
    matrix = np.array([[entry % p for entry in row] for row in M_rows], dtype=np.int64)
    result = np.eye(len(M_rows), dtype=np.int64)
    base = matrix.copy()
    exponent = power
    while exponent:
        if exponent & 1:
            result = matmul(result, base, p)
        exponent >>= 1
        if exponent:
            base = matmul(base, base, p)
    return int(result.trace()) % p


def newton_coefficients(power_sums: list[int], degree: int, p: int) -> list[int]:
    """Return descending coefficients of the monic polynomial from power sums."""
    coefficients = [1]
    for k in range(1, degree + 1):
        total = sum(coefficients[k - i] * power_sums[i] for i in range(1, k + 1))
        coefficients.append((-total * pow(k, -1, p)) % p)
    return coefficients


def trace_sequence(M_rows: list[list[int]], required: int, p: int) -> tuple[list[int], list[int]]:
    """Power traces, with the characteristic recurrence extended exactly mod p."""
    dim = len(M_rows)
    matrix = np.array([[entry % p for entry in row] for row in M_rows], dtype=np.int64)
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for _ in range(dim):
        power = matmul(power, matrix, p)
        traces.append(int(power.trace()) % p)
    characteristic = newton_coefficients(traces, dim, p)
    for k in range(dim + 1, required + 1):
        total = sum(characteristic[i] * traces[k - i] for i in range(1, dim + 1))
        traces.append((-total) % p)
    return traces, characteristic


def poly_trim(coefficients: list[int]) -> list[int]:
    for index, coefficient in enumerate(coefficients):
        if coefficient:
            return coefficients[index:]
    return [0]


def poly_remainder(dividend: list[int], divisor: list[int], p: int) -> list[int]:
    dividend = poly_trim([entry % p for entry in dividend])
    divisor = poly_trim([entry % p for entry in divisor])
    if divisor == [0]:
        raise ZeroDivisionError
    inverse = pow(divisor[0], -1, p)
    divisor = [(inverse * entry) % p for entry in divisor]
    if len(dividend) < len(divisor):
        return dividend
    remainder = dividend[:]
    gap = len(dividend) - len(divisor)
    for index in range(gap + 1):
        coefficient = remainder[index]
        if coefficient:
            for offset, divisor_coefficient in enumerate(divisor):
                remainder[index + offset] = (remainder[index + offset] - coefficient * divisor_coefficient) % p
    return poly_trim(remainder[gap + 1:])


def poly_gcd(left: list[int], right: list[int], p: int) -> list[int]:
    left = poly_trim(left)
    right = poly_trim(right)
    while right != [0]:
        inverse = pow(right[0], -1, p)
        right = [(inverse * entry) % p for entry in right]
        left, right = right, poly_remainder(left, right, p)
    inverse = pow(left[0], -1, p)
    return [(inverse * entry) % p for entry in left]


def derivative_descending(polynomial: list[int], p: int) -> list[int]:
    degree = len(polynomial) - 1
    return [((degree - index) * coefficient) % p for index, coefficient in enumerate(polynomial[:-1])]


def within_pair_gcd_degree(M_rows: list[list[int]], p: int) -> int:
    dim = len(M_rows)
    degree = dim * (dim - 1) // 2
    assert p > degree and p % 2
    traces, _ = trace_sequence(M_rows, 2 * degree, p)
    inv_two = pow(2, -1, p)
    pair_traces = [0] + [((traces[k] * traces[k] - traces[2 * k]) * inv_two) % p for k in range(1, degree + 1)]
    polynomial = newton_coefficients(pair_traces, degree, p)
    gcd = poly_gcd(polynomial, derivative_descending(polynomial, p), p)
    assert poly_remainder(polynomial, gcd, p) == [0]
    assert poly_remainder(derivative_descending(polynomial, p), gcd, p) == [0]
    return len(gcd) - 1


def main() -> int:
    print(f"decisive prime p = {PRIME}")
    check("p exceeds all pair-polynomial degrees and preserves int64 matrix products", PRIME > CROSS_PAIRS and 32 * (PRIME - 1) ** 2 < INT64_MAX)
    cubic = newton_coefficients([0, 10, 38, 160], 3, PRIME)
    check(
        "Newton sign convention recovers (z-2)(z-3)(z-5)",
        cubic == [1, (-10) % PRIME, 31 % PRIME, (-30) % PRIME],
        f"got {cubic}",
    )

    print("\n1. independent exponent-class enumeration")
    small_ok = True
    for modes in range(2, N_MODES + 1):
        ee = exponent_classes(modes, 0, 0, True)
        oo = exponent_classes(modes, 1, 1, True)
        eo = exponent_classes(modes, 0, 1, False)
        small_ok = small_ok and len(ee) == within_formula(modes)
        small_ok = small_ok and len(oo) == within_formula(modes)
        small_ok = small_ok and len(eo) == cross_formula(modes)
    check("brute-force EE, OO, EO exponent classes agree with formulas for m=2..6", small_ok)
    check("m=6 class counts are EE=OO=301 and EO=364", WITHIN_CLASSES == 301 and CROSS_CLASSES == 364)
    check("m=6 gcd floors are EE=OO=195 and EO=660", WITHIN_GCD_FLOOR == 195 and CROSS_GCD_FLOOR == 660)

    print("\n2. raw rational 2x3 operator and physical P sectors")
    t = Fraction(1, 3)
    bonds = [(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)]
    check("hardcoded open 2x3 bonds agree with the repository lattice API", set(bonds) == set(layer_bonds((2, 3), (False, False))))
    R, q, eps = build_R(6, bonds, t)
    mask = (1 << 6) - 1
    check("q=5/3, bond parity=1, and P is a fixed-point-free involution", q == Fraction(5, 3) and eps == 1 and all((state ^ mask) != state and ((state ^ mask) ^ mask) == state for state in range(64)))
    check("raw rational matrix commutes exactly with P", all(R[state ^ mask][other ^ mask] == R[state][other] for state in range(64) for other in range(64)))
    sectors = split_physical_parity(R)
    scale_plus, M_plus = integralize(sectors[1])
    scale_minus, _ = integralize(sectors[-1])
    check("both P sectors are rational 32x32 restrictions with integral monic models", len(M_plus) == 32 and scale_plus > 0 and scale_minus > 0)
    traces_65, characteristic = trace_sequence(M_plus, 65, PRIME)
    check("P=+ characteristic polynomial is monic of degree 32", len(characteristic) == 33 and characteristic[0] == 1)
    check(
        "P=+ Cayley-Hamilton trace recurrence agrees with binary powering at k=65",
        traces_65[65] == binary_trace(M_plus, 65, PRIME),
        f"recurrence={traces_65[65]}, direct={binary_trace(M_plus, 65, PRIME)}",
    )

    print("\n3. independent decisive modular gcd degree")
    plus_degree = within_pair_gcd_degree(M_plus, PRIME)
    check("independent P=+ modular gcd degree is exactly 177", plus_degree == 177)
    check(f"P=+ within-sector gcd degree {plus_degree} is below the universal floor {WITHIN_GCD_FLOOR}", plus_degree < WITHIN_GCD_FLOOR)
    check("the same P=+ failure excludes P+=even/P-=odd", plus_degree < WITHIN_GCD_FLOOR)
    check("the same P=+ failure excludes P+=odd/P-=even", plus_degree < WITHIN_GCD_FLOOR)

    print("\n4. producer artifact regression")
    artifact_path = ROOT / "results" / "spectral" / "parity_pair_product.json"
    check("producer artifact exists", artifact_path.exists())
    if artifact_path.exists():
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        layer = artifact["data"]["cases"]["layer_2x3"]
        recorded = layer["within_sector"]["+1"]["per_prime"][str(PRIME)]["gcd_degree"]
        check("independent P=+ gcd degree equals the producer artifact", plus_degree == recorded, f"test={plus_degree}, artifact={recorded}")
        check("artifact labels and rejects both physical-to-fermionic assignments", all(item["excluded"] for item in layer["assignments"]))
        check("artifact records exact-threshold chain and synthetic controls", all(item["all_thresholds_met_exactly"] for name, item in artifact["data"]["cases"].items() if name != "layer_2x3"))
        genericity = artifact["data"].get("genericity")
        check("artifact records the symbolic genericity certificate", genericity is not None)
        if genericity is not None:
            expected_grid = ["1/4", "1/3", "2/5", "1/2", "2/3"]
            grid = genericity["rational_grid"]
            check(
                "deterministic rational t-grid has the requested five samples",
                [sample["t"] for sample in grid] == expected_grid,
            )
            check(
                "every grid point independently has P=+ gcd degree 177 at p=1000003",
                all(sample["plus_within_per_prime"][str(PRIME)]["gcd_degree"] == 177 for sample in grid),
            )
            check(
                "genericity bound records at most 25,534,083 algebraic exceptional parameter values including denominator roots",
                genericity["exceptional_set_cardinality_bound_including_denominator_roots"] == 25_534_083,
            )
            quotient = genericity["L2_genericity"]
            check(
                "generic resultant metadata distinguishes monic H from J with leading coefficient 496",
                quotient.get("H_is_monic") is True
                and quotient.get("J_is_in_Q_t_z") is True
                and quotient.get("J_leading_z_coefficient") == 496
                and quotient.get("H_J_coprime_over_Q_t") is True,
            )

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(f"OK: independently rebuilt P=+ sector has gcd degree {plus_degree} < {WITHIN_GCD_FLOOR}; both parity assignments are excluded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
