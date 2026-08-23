"""Standalone verifier for the layer-group spectral obstruction.

This test imports none of e188--e190.  It independently checks the exact
varying-conjugator witness and rebuilds the 2x3/2x4 t=1/3 modular spectral
controls directly in F_p from P_t D_q P_t.  The producer instead constructs a
scaled integer matrix by Kronecker sweeps.
"""

from __future__ import annotations

import itertools
import json
import platform
import resource
import time
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "layer_group_spectral.json"
PRIME = 1_000_003
CPU_LIMIT = 180.0
RSS_LIMIT = 2_000_000_000
FLOAT_EXACT = 1 << 53
INT64_MAX = (1 << 63) - 1
STARTED = time.process_time()
FAILURES: list[str] = []


def rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget(stage: str) -> None:
    used = time.process_time() - STARTED
    if used >= CPU_LIMIT:
        raise RuntimeError(f"CPU budget exceeded at {stage}: {used:.3f}s")
    if rss_bytes() >= RSS_LIMIT:
        raise RuntimeError(f"RSS budget exceeded at {stage}: {rss_bytes()}")


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not passed:
        FAILURES.append(name)


# ---------------------------------------------------------------- formalization, no producer import

def verify_spin_weights() -> None:
    for modes in range(1, 9):
        signs = list(itertools.product((-1, 1), repeat=modes))
        normalized = {tuple((entry + 1) // 2 for entry in row) for row in signs}
        boolean_cube = set(itertools.product((0, 1), repeat=modes))
        check(
            f"n={modes}: spin weights normalize to all subset exponents",
            len(signs) == 1 << modes and normalized == boolean_cube,
        )


def physical_coordinates(t: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    x = (1 + t) / (1 - t)
    y = (1 + t * t) / (2 * t)
    curve = (y - 1) * x * x - (y + 1)
    return x, y, curve


def matmul2(left: list[list[Fraction]], right: list[list[Fraction]]) -> list[list[Fraction]]:
    return [
        [sum((left[i][k] * right[k][j] for k in range(2)), Fraction(0)) for j in range(2)]
        for i in range(2)
    ]


def matsub2(left: list[list[Fraction]], right: list[list[Fraction]]) -> list[list[Fraction]]:
    return [[left[i][j] - right[i][j] for j in range(2)] for i in range(2)]


def verify_varying_conjugator() -> None:
    H = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(-1)]]
    X = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    commutator = matsub2(matmul2(H, X), matmul2(X, H))
    flattened = [
        [H[i][j], X[i][j], commutator[i][j]] for i in range(2) for j in range(2)
    ]
    # The 3x3 minor from matrix slots (00,01,10) is nonzero.
    determinant_minor = (
        flattened[0][0]
        * (flattened[1][1] * flattened[2][2] - flattened[1][2] * flattened[2][1])
        - flattened[0][1]
        * (flattened[1][0] * flattened[2][2] - flattened[1][2] * flattened[2][0])
        + flattened[0][2]
        * (flattened[1][0] * flattened[2][1] - flattened[1][1] * flattened[2][0])
    )
    check("H,X,[H,X] are linearly independent", determinant_minor != 0, str(determinant_minor))

    for r in (Fraction(2), Fraction(3, 2), Fraction(5, 3)):
        D = [[r, Fraction(0)], [Fraction(0), 1 / r]]
        Y = [
            [Fraction(0), -(r * r) / (r * r - 1)],
            [1 / (r * r - 1), Fraction(0)],
        ]
        orbit_velocity = matsub2(matmul2(Y, D), matmul2(D, Y))
        transverse_velocity = matmul2(D, X)
        check(
            f"r={r}: conjugator velocity absorbs D X",
            orbit_velocity == transverse_velocity,
        )

    r, u = Fraction(2), Fraction(2)
    c, s = (u + 1 / u) / 2, (u - 1 / u) / 2
    sample = [[r * c, r * s], [s / r, c / r]]
    determinant = sample[0][0] * sample[1][1] - sample[0][1] * sample[1][0]
    trace = sample[0][0] + sample[1][1]
    discriminant = trace * trace - 4 * determinant
    check(
        "regular-semisimple non-diagonal witness sample",
        determinant == 1 and discriminant == Fraction(369, 64) and sample[0][1] * sample[1][0] != 0,
        f"det={determinant}, Delta={discriminant}",
    )


# ---------------------------------------------------------- independent F_p layer builder

def grid_edges(cols: int) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for row in range(2):
        for col in range(cols):
            site = row * cols + col
            if col + 1 < cols:
                edges.append((site, site + 1))
            if row == 0:
                edges.append((site, site + cols))
    return edges


def direct_layer_mod(cols: int, prime: int) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Build rational R=P_t D_q P_t directly modulo p, with no integer sweep."""
    sites = 2 * cols
    dimension = 1 << sites
    edges = grid_edges(cols)
    t = pow(3, -1, prime)
    q = 5 * pow(3, -1, prime) % prime
    t_powers = [pow(t, degree, prime) for degree in range(sites + 1)]
    P = np.empty((dimension, dimension), dtype=np.int64)
    for row in range(dimension):
        for col in range(dimension):
            P[row, col] = t_powers[(row ^ col).bit_count()]
    diagonal = np.empty(dimension, dtype=np.int64)
    for state in range(dimension):
        aligned = 0
        for left, right in edges:
            left_bit = (state >> (sites - 1 - left)) & 1
            right_bit = (state >> (sites - 1 - right)) & 1
            aligned += left_bit == right_bit
        diagonal[state] = pow(q, aligned, prime)
    budget(f"direct {2}x{cols} matrix construction")
    return mod_product((P * diagonal[np.newaxis, :]) % prime, P, prime), edges


def mod_product(left: np.ndarray, right: np.ndarray, prime: int) -> np.ndarray:
    assert left.shape[1] * (prime - 1) ** 2 < FLOAT_EXACT
    product = np.asarray(left, dtype=np.float64) @ np.asarray(right, dtype=np.float64)
    return np.mod(product, float(prime)).astype(np.int64)


def dot_mod(left: np.ndarray, right: np.ndarray, prime: int) -> int:
    assert left.size * (prime - 1) ** 2 < INT64_MAX
    return int(np.dot(left, right)) % prime


def strip_leading_zeros(poly: np.ndarray) -> np.ndarray:
    nonzero = np.flatnonzero(poly)
    return np.zeros(1, dtype=np.int64) if nonzero.size == 0 else poly[nonzero[0] :]


def make_monic(poly: np.ndarray, prime: int) -> np.ndarray:
    poly = strip_leading_zeros(np.asarray(poly, dtype=np.int64) % prime)
    if poly.size == 1 and int(poly[0]) == 0:
        return poly
    return poly if int(poly[0]) == 1 else poly * pow(int(poly[0]), -1, prime) % prime


def polynomial_remainder(dividend: np.ndarray, divisor: np.ndarray, prime: int) -> np.ndarray:
    work = np.asarray(dividend, dtype=np.int64).copy() % prime
    divisor_degree = divisor.size - 1
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return strip_leading_zeros(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset % 256 == 0:
            budget("independent polynomial division")
        leading = int(work[offset])
        if leading:
            block = work[offset : offset + divisor_degree + 1]
            work[offset : offset + divisor_degree + 1] = (block - leading * divisor) % prime
    return strip_leading_zeros(work[dividend_degree - divisor_degree + 1 :])


def gcd_polynomials(left: np.ndarray, right: np.ndarray, prime: int) -> np.ndarray:
    left, right = make_monic(left, prime), make_monic(right, prime)
    while not (right.size == 1 and int(right[0]) == 0):
        left, right = right, make_monic(polynomial_remainder(left, right, prime), prime)
    return left


def polynomial_derivative(poly: np.ndarray, prime: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    return np.ascontiguousarray(
        ((np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % prime)[::-1]
    )


def coefficients_from_traces(
    traces: np.ndarray, degree: int, prime: int, inverses: list[int]
) -> np.ndarray:
    chunk = (FLOAT_EXACT - 1) // ((prime - 1) ** 2)
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    sums = np.asarray(traces[: degree + 1], dtype=np.float64)
    reversed_coefficients = np.zeros(degree + 1, dtype=np.float64)
    reversed_coefficients[degree] = 1.0
    for k in range(1, degree + 1):
        if k % 1024 == 0:
            budget(f"independent Newton degree {degree}: {k}")
        total = 0
        start = 1
        base = degree - k
        while start <= k:
            stop = min(k, start + chunk - 1)
            total += int(
                np.dot(
                    sums[start : stop + 1],
                    reversed_coefficients[base + start : base + stop + 1],
                )
            ) % prime
            start = stop + 1
        coefficient = (-(total % prime) * inverses[k]) % prime
        coefficients[k] = coefficient
        reversed_coefficients[degree - k] = float(coefficient)
    return coefficients


def independent_pair_gcd_degree(matrix: np.ndarray, prime: int) -> tuple[int, int, bool]:
    dimension = matrix.shape[0]
    slots = comb(dimension, 2)
    needed = 2 * slots
    power = np.eye(dimension, dtype=np.int64)
    traces = [0]
    for exponent in range(1, dimension + 1):
        power = mod_product(power, matrix, prime)
        traces.append(int(power.trace()) % prime)
        if exponent % 32 == 0:
            budget(f"independent traces {exponent}/{dimension}")
    inverses = [0] + [pow(index, -1, prime) for index in range(1, needed + 1)]
    characteristic = coefficients_from_traces(
        np.array(traces, dtype=np.int64), dimension, prime, inverses
    )
    tail = characteristic[1:]
    for exponent in range(dimension + 1, needed + 1):
        window = np.array(
            traces[exponent - 1 : exponent - dimension - 1 : -1], dtype=np.int64
        )
        traces.append((-dot_mod(tail, window, prime)) % prime)
        if exponent % 2048 == 0:
            budget(f"independent trace recurrence {exponent}/{needed}")
    half = pow(2, -1, prime)
    pair_traces = np.zeros(slots + 1, dtype=np.int64)
    for exponent in range(1, slots + 1):
        pair_traces[exponent] = (
            (traces[exponent] * traces[exponent] - traces[2 * exponent]) * half
        ) % prime
    pair_polynomial = coefficients_from_traces(pair_traces, slots, prime, inverses)
    gcd_polynomial = gcd_polynomials(
        pair_polynomial, polynomial_derivative(pair_polynomial, prime), prime
    )
    divides = polynomial_remainder(pair_polynomial, gcd_polynomial, prime)
    divides_ok = divides.size == 1 and int(divides[0]) == 0
    return gcd_polynomial.size - 1, slots - (gcd_polynomial.size - 1), divides_ok


def rederive_controls() -> dict[str, dict[str, int | bool]]:
    expected = {
        "2x3": (3, 6, 7, 385, 1631, 665),
        "2x4": (4, 8, 10, 9329, 23311, 6305),
    }
    rows: dict[str, dict[str, int | bool]] = {}
    for label, (cols, sites, edge_count, expected_gcd, expected_r, ceiling) in expected.items():
        matrix, edges = direct_layer_mod(cols, PRIME)
        check(f"{label}: independently built matrix is symmetric", bool(np.array_equal(matrix, matrix.T)))
        check(f"{label}: raw edge count", len(edges) == edge_count, str(len(edges)))
        gcd_degree, lower_bound, divides = independent_pair_gcd_degree(matrix, PRIME)
        check(
            f"{label}: independent exact modular gcd degree",
            gcd_degree == expected_gcd,
            f"{gcd_degree} at p={PRIME}",
        )
        check(f"{label}: independent gcd division audit", divides)
        check(
            f"{label}: independent one-sided spectral certificate",
            lower_bound == expected_r and lower_bound > ceiling and ceiling == 3**sites - 2**sites,
            f"r_Q>={lower_bound}>{ceiling}",
        )
        rows[label] = {
            "sites": sites,
            "edge_count": len(edges),
            "pair_gcd_degree": gcd_degree,
            "distinct_pair_products_lower_bound": lower_bound,
            "gaussian_ceiling": ceiling,
        }
    return rows


print("1. formalization and physical curve")
verify_spin_weights()
for value in (Fraction(1, 5), Fraction(1, 3), Fraction(2, 5), Fraction(3, 4)):
    x_value, y_value, residue = physical_coordinates(value)
    check(
        f"t={value}: exact physical curve equation",
        residue == 0 and (x_value - 1) / (x_value + 1) == value,
        f"x={x_value}, y={y_value}",
    )

print("\n2. exact varying-conjugator witness")
verify_varying_conjugator()

print("\n3. independent 2x3 and 2x4 controls")
independent_rows = rederive_controls()

print("\n4. artifact consistency and scope discipline")
with ARTIFACT.open() as handle:
    artifact = json.load(handle)
check("artifact top-level schema", set(artifact) == {"meta", "data", "checks"})
check("every producer check passed", all(bool(item["passed"]) for item in artifact["checks"]))
check(
    "artifact does not claim the isotropic all-size theorem",
    artifact["data"]["conclusion"]["isotropic_all_size_spectral_nogo_proved"] is False
    and artifact["data"]["conclusion"]["exceptional_set"] is None,
)
stored_rows = {row["layer"]: row for row in artifact["data"]["control_computations"]}
check("artifact has exactly the two declared controls", set(stored_rows) == set(independent_rows))
for label, independent in independent_rows.items():
    stored = stored_rows[label]
    check(
        f"{label}: stored control matches independent reconstruction",
        all(int(stored[key]) == int(value) for key, value in independent.items()),
    )
check(
    "artifact records the exact obstruction rather than a spectral converse",
    artifact["data"]["varying_conjugator_obstruction"]["tag"] == "[THEOREM]"
    and artifact["data"]["conclusion"]["tag"] == "[UNRESOLVED]",
)

elapsed = time.process_time() - STARTED
check(
    "standalone verifier resource budget",
    elapsed < CPU_LIMIT and rss_bytes() < RSS_LIMIT,
    f"CPU={elapsed:.3f}s<{CPU_LIMIT}; RSS={rss_bytes()}<{RSS_LIMIT}",
)

if FAILURES:
    print(f"FAIL test_layer_group_spectral ({len(FAILURES)} failures: {', '.join(FAILURES)})")
    raise SystemExit(1)
print("PASS test_layer_group_spectral")
