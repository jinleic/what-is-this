"""Standalone independent verifier for the open-2x5 pair-product certificate.

This test deliberately imports neither experiments/e136_pair_product_2x5.py,
e122_pair_product_3x3.py, nor e38_gaussianity_certificate.py.  It independently:

  1. re-derives the canonical rational R = P_t diag(q^e) P_t for every control
     row from raw inputs -- as actual Fraction objects through dim 512 (2x3 and
     the 3x3 digest-replay row), and through an exact pure-Python integer
     numerator reconstruction W = P~ diag(d~) P~ at dim 1024 (the decisive 2x5
     layer and the n=10 open-chain floor control), where every one of the
     1024x1024 numerators is an exact math.sumprod of Python integers with no
     numpy matmul, no split, and no producer code;
  2. re-derives N0 = den(t)^{2n} num(q)^E den(q)^F from its own exponent scan,
     computes the content G by a FULL scan (no early break), the actual LCM D
     both as N0/G and directly as the LCM of every entry's reduced denominator,
     and anchors the integer route to true Fraction arithmetic on 23
     deterministic entries;
  3. reduces M = W/G modulo each declared prime, rebuilds charpoly(M), the
     pair polynomial charpoly(wedge^2 M), and its gcd with the derivative,
     using strictly reduced long division at every update and naive
     reversed-slice Newton loops (different code paths from the producer);
  4. and only then reads the stored artifacts and compares values and
     coefficient digests row by row.

Compute gates read time.process_time(), never wall clock and never
signal.alarm; the strict-gcd stage runs on its own local clock.

Run: PYTHONPATH=src .venv/bin/python tests/test_pair_product_2x5.py
"""

from __future__ import annotations

import gc
import hashlib
import json
import platform
import resource
import sys
import time
from fractions import Fraction
from math import gcd as math_gcd
from math import isqrt, sumprod
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402

INT64_MAX = 2**63 - 1
FLOAT_EXACT_MAX = 2**53
PRIMES = (1_000_003, 2_000_003)
RSS_LIMIT_BYTES = 6_000_000_000
TEST_WALL_CPU_SECONDS = 10_800
GCD_WALL_CPU_SECONDS = 5_400
ARTIFACT = ROOT / "results" / "spectral" / "pair_product_2x5.json"
WAVE13_ARTIFACT = ROOT / "results" / "spectral" / "pair_product_3x3.json"

FAILS: list[str] = []
STARTED_CPU = time.process_time()
STARTED_WALL = time.perf_counter()


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)
    return bool(ok)


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def enforce(stage: str, wall_seconds: int = TEST_WALL_CPU_SECONDS, *, started: float | None = None) -> None:
    """Process-time gate: CPU seconds, never wall clock and never signals."""
    anchor = STARTED_CPU if started is None else started
    elapsed = time.process_time() - anchor
    assert elapsed < wall_seconds, f"CPU wall {wall_seconds}s fired at {stage} after {elapsed:.3f} CPU s"
    assert max_rss_bytes() < RSS_LIMIT_BYTES, f"RSS wall at {stage}: {max_rss_bytes()} >= {RSS_LIMIT_BYTES}"


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def assert_int64_guard(length: int, p: int, where: str) -> None:
    assert length * (p - 1) ** 2 < INT64_MAX, f"int64 guard failed at {where}"


def assert_float_guard(length: int, p: int, where: str) -> None:
    assert length * (p - 1) ** 2 < FLOAT_EXACT_MAX, f"binary64 exactness guard failed at {where}"


# ------------------------------------------------- shared raw rational ingredients

def raw_geometry(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Independent spin/bond-sum scan: q, parity eps, e exponents, energy range."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    energy = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energy.append(sum(spins[left] * spins[right] for left, right in bonds))
    epsilon = energy[0] % 2
    assert all((value - epsilon) % 2 == 0 for value in energy)
    exponents = [(value - epsilon) // 2 for value in energy]
    return q, epsilon, exponents, min(energy), max(energy)


def raw_build_r_fractions(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Full Fraction-object R = P_t diag(q^e) P_t (feasible through dim 512)."""
    q, epsilon, exponents, energy_min, energy_max = raw_geometry(n, bonds, t)
    dim = 1 << n
    weights = [q**e for e in exponents]
    t_powers = [t**hamming for hamming in range(n + 1)]
    kernel = [[t_powers[(left ^ right).bit_count()] for right in range(dim)] for left in range(dim)]
    matrix = [[Fraction(0) for _ in range(dim)] for _ in range(dim)]
    for left in range(dim):
        if left & 31 == 0:
            enforce("raw Fraction build")
        for right in range(left, dim):
            entry = Fraction(0)
            for middle in range(dim):
                entry += kernel[left][middle] * weights[middle] * kernel[middle][right]
            matrix[left][right] = matrix[right][left] = entry
    return matrix, q, epsilon, energy_min, energy_max


def actual_integral_lcm_fractions(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    """Actual LCM of all entry denominators, followed by exact M = D*R."""
    scale = 1
    for row in matrix:
        for entry in row:
            scale = scale * entry.denominator // math_gcd(scale, entry.denominator)
    integer_rows = []
    for row in matrix:
        output = []
        for entry in row:
            value = scale * entry
            assert value.denominator == 1
            output.append(value.numerator)
        integer_rows.append(output)
    return scale, integer_rows


def raw_build_w_exact_integers(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Exact pure-Python integer numerators W = P~ diag(d~) P~ at dim 1024.

    P~ = den(t)^n P_t has integer entries num(t)^h den(t)^(n-h); the diagonal
    d~_k = q^{e_k} num(q)^E den(q)^F = num^{e+E} den^{F-e} is integral for the
    measured exponent range.  Every W_ij below is a math.sumprod of Python
    integers -- no numpy matmul, no split, no producer code.  Returns W, the
    independently derived N0, the full-scan content G, and both derivations of
    the actual LCM D.
    """
    q, epsilon, exponents, energy_min, energy_max = raw_geometry(n, bonds, t)
    e_min, e_max = min(exponents), max(exponents)
    num, den = q.numerator, q.denominator
    a, bden = t.numerator, t.denominator
    E, F = max(0, -e_min), max(0, e_max)
    dim = 1 << n
    enforce("integer kernel tables")
    kernel_table = [a**h * bden ** (n - h) for h in range(n + 1)]
    kernel = [[kernel_table[(left ^ right).bit_count()] for right in range(dim)] for left in range(dim)]
    diagonal = [num ** (e + E) * den ** (F - e) for e in exponents]
    enforce("integer numerator rows")
    rows = []
    for left in range(dim):
        if left & 63 == 0:
            enforce("integer numerator rows")
        scaled = [kernel[left][k] * diagonal[k] for k in range(dim)]
        rows.append([sumprod(scaled, kernel[right]) for right in range(dim)])
    del kernel
    N0 = bden ** (2 * n) * num**E * den**F
    # FULL content scan: no early break, every entry participates.
    content = 0
    for row in rows:
        for value in row:
            content = math_gcd(content, value)
    G = math_gcd(content, N0)
    # Independent actual-LCM derivation: LCM over entries of N0/gcd(W_ij, N0).
    actual_lcm = 1
    for row in rows:
        for value in row:
            reduced_denominator = N0 // math_gcd(value, N0)
            actual_lcm = actual_lcm * reduced_denominator // math_gcd(actual_lcm, reduced_denominator)
        enforce("actual LCM scan")
    assert actual_lcm == N0 // G, "the two actual-LCM derivations disagree"
    return rows, N0, G, N0 // G, q, epsilon, e_min, e_max, E, F, energy_min, energy_max


def fraction_anchor_entries(
    n: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    rows: list[list[int]],
    N0: int,
    random_entries: int = 16,
) -> list[tuple[int, int, bool]]:
    """Anchor the integer route to true Fraction arithmetic on a fixed stride and seeded entries."""
    q, epsilon, exponents, _, _ = raw_geometry(n, bonds, t)
    dim = 1 << n
    weights = [q**e for e in exponents]
    t_powers = [t**hamming for hamming in range(n + 1)]
    rng = np.random.default_rng(20260818)
    entries = [(0, j) for j in range(0, dim, 149)] + [
        (int(i), int(j))
        for i, j in zip(rng.integers(0, dim, random_entries), rng.integers(0, dim, random_entries))
    ]
    results = []
    for (i, j) in entries:
        left = [t_powers[(i ^ m).bit_count()] for m in range(dim)]
        right = [t_powers[(j ^ m).bit_count()] for m in range(dim)]
        acc = Fraction(0)
        for m in range(dim):
            acc += left[m] * weights[m] * right[m]
        results.append((i, j, acc * N0 == rows[i][j]))
    return results


# ------------------------------------------------ exact modular arithmetic / traces

def float_matmul_mod(left: np.ndarray, right: np.ndarray, p: int) -> np.ndarray:
    """Exact modular product through binary64 with a runtime integral round trip."""
    assert_float_guard(left.shape[1], p, "float matrix product")
    raw = left.astype(np.float64) @ right.astype(np.float64)
    integer = raw.astype(np.int64)
    assert np.array_equal(raw, integer.astype(np.float64)), "non-integral binary64 matrix product"
    return integer % p


def int_matmul_mod(left: np.ndarray, right: np.ndarray, p: int) -> np.ndarray:
    assert_int64_guard(left.shape[1], p, "int64 matrix product")
    return (left @ right) % p


def power_traces(matrix: np.ndarray, p: int) -> np.ndarray:
    dim = matrix.shape[0]
    assert_float_guard(dim, p, "power traces")
    power = np.eye(dim, dtype=np.float64)
    base = matrix.astype(np.float64)
    traces = np.zeros(dim + 1, dtype=np.int64)
    for exponent in range(1, dim + 1):
        if exponent & 63 == 0:
            enforce("float64 power traces")
        raw = power @ base
        reduced = raw.astype(np.int64)
        assert np.array_equal(raw, reduced.astype(np.float64)), "non-integral binary64 trace power"
        reduced %= p
        traces[exponent] = int(reduced.trace()) % p
        power = reduced.astype(np.float64)
    return traces


def binary_trace_int64(matrix: np.ndarray, exponent: int, p: int) -> int:
    """Separate all-int64 binary route for recurrence checks."""
    result = np.eye(matrix.shape[0], dtype=np.int64)
    base = matrix.copy()
    e = exponent
    while e:
        if e & 1:
            result = int_matmul_mod(result, base, p)
        e >>= 1
        if e:
            base = int_matmul_mod(base, base, p)
    return int(result.trace()) % p


# ---------------------------------------------------- independent Newton/recurrent path

def newton(power_sums: np.ndarray, degree: int, p: int) -> np.ndarray:
    """[1,a1,...,a_degree] via naive reversed-slice Newton loops."""
    assert_int64_guard(degree, p, "Newton")
    coefficients = np.zeros(degree + 1, dtype=np.int64)
    coefficients[0] = 1
    for k in range(1, degree + 1):
        if k & 1023 == 0:
            enforce("Newton reconstruction")
        total = int(np.dot(power_sums[1 : k + 1], np.ascontiguousarray(coefficients[k - 1 :: -1]))) % p
        coefficients[k] = (-total * pow(k, -1, p)) % p
    return coefficients


def extend_traces(initial: np.ndarray, characteristic: np.ndarray, upto: int, p: int) -> np.ndarray:
    dim = characteristic.size - 1
    assert_int64_guard(dim, p, "trace recurrence")
    output = np.zeros(upto + 1, dtype=np.int64)
    output[: initial.size] = initial
    tail = characteristic[1:]
    for exponent in range(dim + 1, upto + 1):
        if exponent & 4095 == 0:
            enforce("trace recurrence")
        output[exponent] = (-int(np.dot(tail, output[exponent - 1 : exponent - dim - 1 : -1]))) % p
    return output


def horner_zero(characteristic: np.ndarray, matrix: np.ndarray, p: int) -> bool:
    """Independent Cayley-Hamilton check through exact binary64 products."""
    dim = matrix.shape[0]
    value = np.eye(dim, dtype=np.int64)
    identity = np.eye(dim, dtype=np.int64)
    for coefficient in characteristic[1:]:
        enforce("Horner Cayley-Hamilton", 10_800)
        value = (float_matmul_mod(value, matrix, p) + int(coefficient) * identity) % p
    return bool(np.all(value == 0))


# ------------------------------------------------ strict (non-lazy) F_p[z] Euclid

def trim(poly: np.ndarray) -> np.ndarray:
    poly = np.asarray(poly, dtype=np.int64)
    if poly.size == 1:
        return poly
    if poly[0] != 0:
        return np.ascontiguousarray(poly)
    indices = np.flatnonzero(poly)
    return np.ascontiguousarray(poly[indices[0] :]) if indices.size else np.zeros(1, dtype=np.int64)


def zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = trim(poly % p)
    return poly if zero(poly) else (poly * pow(int(poly[0]), -1, p)) % p


def strict_remainder(dividend: np.ndarray, divisor_monic: np.ndarray, p: int, gcd_started: float) -> np.ndarray:
    """Fully reduced long division; unlike the producer, reduces every update mod p.

    The divisor must already be monic with values in [0, p), so the lead
    coefficient itself is the quotient factor.
    """
    divisor_degree = divisor_monic.size - 1
    work = np.asarray(dividend, dtype=np.int64).copy() % p
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset & 1023 == 0:
            enforce("strict Euclidean division", GCD_WALL_CPU_SECONDS, started=gcd_started)
        lead = int(work[offset])
        if lead:
            work[offset : offset + divisor_degree + 1] = (
                work[offset : offset + divisor_degree + 1] - lead * divisor_monic
            ) % p
    return trim(work[dividend_degree - divisor_degree + 1 :])


def strict_gcd(first: np.ndarray, second: np.ndarray, p: int) -> tuple[np.ndarray, int]:
    gcd_started = time.process_time()
    first, second = monic(first, p), monic(second, p)
    steps = 0
    while not zero(second):
        if steps & 127 == 0:
            enforce("strict Euclidean gcd", GCD_WALL_CPU_SECONDS, started=gcd_started)
        first, second = second, monic(strict_remainder(first, second, p, gcd_started), p)
        steps += 1
    return first, steps


def local_gcd_wall_probe() -> bool:
    """A strict-gcd stage must use its own process-time clock, not prior work."""
    try:
        strict_gcd(
            np.array([1, PRIMES[0] - 2, 1], dtype=np.int64),
            np.array([1, PRIMES[0] - 1], dtype=np.int64),
            PRIMES[0],
        )
    except AssertionError:
        return False
    # A genuinely exhausted local clock must fire on a fresh division.
    fake_started = time.process_time() - GCD_WALL_CPU_SECONDS - 1
    try:
        strict_remainder(
            np.zeros(4096, dtype=np.int64), np.ones(2048, dtype=np.int64), PRIMES[0], fake_started
        )
        return False
    except AssertionError:
        return True


def naive_polynomial_remainder(dividend: list[int], divisor: list[int], p: int) -> list[int]:
    """Pure-Python ascending-coefficient reference division for self-validation."""
    divisor = list(divisor)
    while divisor and divisor[-1] % p == 0:
        divisor.pop()
    assert divisor and divisor[-1] % p != 0
    inv_lead = pow(divisor[-1], -1, p)
    work = [c % p for c in dividend]
    while len(work) >= len(divisor):
        if len(work) & 1023 == 0:
            enforce("naive reference division")
        factor = work[-1] * inv_lead % p
        shift = len(work) - len(divisor)
        for index, coefficient in enumerate(divisor):
            work[shift + index] = (work[shift + index] - factor * coefficient) % p
        work.pop()
    while len(work) > 1 and work[-1] == 0:
        work.pop()
    return work


def derivative_descending(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    return np.ascontiguousarray((np.arange(degree, 0, -1, dtype=np.int64) * poly[:-1]) % p)


def digest(ascending) -> str:
    return hashlib.sha256(",".join(str(int(value)) for value in ascending).encode()).hexdigest()


def strict_ops_self_validation(p: int) -> bool:
    """strict_remainder must agree with the pure-Python reference on random input."""
    rng = np.random.default_rng(20260819)
    started = time.process_time()
    for trial in range(25):
        dd = int(rng.integers(2, 300))
        dv = int(rng.integers(1, dd))
        dividend = [int(v) for v in rng.integers(0, p, dd + 1)]
        divisor = [int(v) for v in rng.integers(1, p, dv + 1)]
        reference = naive_polynomial_remainder(dividend, divisor, p)
        got = strict_remainder(
            np.array(dividend[::-1], dtype=np.int64), monic(np.array(divisor[::-1], dtype=np.int64), p), p, started
        )
        if [int(v) for v in got[::-1] % p] != reference:
            return False
    return True


# ---------------------------------------------------------------- independent full pipeline

def independent_pipeline(integer_rows: list[list[int]], p: int) -> dict[str, object]:
    started = time.perf_counter()
    dim = len(integer_rows)
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > pairs and p % 2 == 1
    matrix = np.array([[entry % p for entry in row] for row in integer_rows], dtype=np.int64)

    first_product_equal = np.array_equal(float_matmul_mod(matrix, matrix, p), int_matmul_mod(matrix, matrix, p))
    traces = power_traces(matrix, p)
    characteristic = newton(traces, dim, p)
    horner_ok = horner_zero(characteristic, matrix, p)
    extended = extend_traces(traces, characteristic, 2 * pairs, p)

    # Direct int64 powers hit an immediate post-charpoly index and the full
    # 2*pairs endpoint: this catches a recurrence sign, a dropped factor of
    # two, and an off-by-one range defect.
    spot_indices = tuple(sorted({min(x, 2 * pairs) for x in (dim + 1, 2 * pairs)}))
    spots = {k: binary_trace_int64(matrix, k, p) == int(extended[k]) for k in spot_indices}

    inverse_two = pow(2, -1, p)
    pair_sums = np.zeros(pairs + 1, dtype=np.int64)
    pair_sums[1:] = (
        (extended[1 : pairs + 1] * extended[1 : pairs + 1] - extended[2 : 2 * pairs + 1 : 2]) % p
    ) * inverse_two % p
    pair_poly = newton(pair_sums, pairs, p)
    derivative = derivative_descending(pair_poly, p)
    common, steps = strict_gcd(pair_poly, derivative, p)
    division_check_started = time.process_time()
    divides = zero(strict_remainder(pair_poly, common, p, division_check_started)) and zero(
        strict_remainder(derivative, common, p, division_check_started)
    )
    return {
        "gcd_degree": int(common.size) - 1,
        "charpoly_digest": digest(characteristic[::-1] % p),
        "pair_poly_digest": digest(pair_poly[::-1] % p),
        "gcd_digest": digest(common[::-1] % p),
        "charpoly_monic": int(characteristic[0]) == 1,
        "pair_poly_monic": int(pair_poly[0]) == 1,
        "horner_zero": horner_ok,
        "first_product_equal": first_product_equal,
        "recurrence_spots": spots,
        "gcd_divides_both": divides,
        "euclid_steps": steps,
        "seconds_wall_reference_only": round(time.perf_counter() - started, 3),
        "cpu_seconds": round(time.process_time(), 3),
        "rss_bytes": max_rss_bytes(),
    }


def run_fraction_case(name: str, n: int, bonds: list[tuple[int, int]], primes: tuple[int, ...]) -> dict[str, object]:
    print(f"raw Fraction reconstruction {name} ...", flush=True)
    started = time.perf_counter()
    rational, q, epsilon, energy_min, energy_max = raw_build_r_fractions(n, bonds, Fraction(1, 3))
    D, rows = actual_integral_lcm_fractions(rational)
    del rational
    gc.collect()
    build_seconds = round(time.perf_counter() - started, 3)
    check(f"{name}: q=5/3 and canonical parity", q == Fraction(5, 3) and epsilon in {0, 1})
    results = {}
    for p in primes:
        print(f"  independent strict pipeline p={p} ...", flush=True)
        results[str(p)] = independent_pipeline(rows, p)
        print(
            f"    gcd={results[str(p)]['gcd_degree']} steps={results[str(p)]['euclid_steps']} "
            f"in {results[str(p)]['seconds_wall_reference_only']} s",
            flush=True,
        )
    del rows
    gc.collect()
    return {
        "D": D,
        "q": q,
        "epsilon": epsilon,
        "energy_min": energy_min,
        "energy_max": energy_max,
        "build_seconds": build_seconds,
        "per_prime": results,
    }


def run_integer_case(name: str, n: int, bonds: list[tuple[int, int]], prime: int) -> dict[str, object]:
    print(f"raw exact-integer reconstruction {name} ...", flush=True)
    started = time.perf_counter()
    rows, N0, G, D, q, epsilon, e_min, e_max, E, F, energy_min, energy_max = raw_build_w_exact_integers(
        n, bonds, Fraction(1, 3)
    )
    anchors = fraction_anchor_entries(n, bonds, Fraction(1, 3), rows, N0)
    build_seconds = round(time.perf_counter() - started, 3)
    check(f"{name}: q=5/3, canonical parity, and measured exponent range", q == Fraction(5, 3))
    check(
        f"{name}: all {len(anchors)} Fraction anchor entries equal W/N0 exactly",
        all(ok for _, _, ok in anchors),
        f"failed: {[e for e, _, ok in anchors if not ok][:4]}",
    )
    print(f"  independent strict pipeline p={prime} ...", flush=True)
    # M = D*R = W/G exactly (G divides every W_ij); the pipeline reduces M.
    monic_rows = [[value // G for value in row] for row in rows]
    result = independent_pipeline(monic_rows, prime)
    print(
        f"    gcd={result['gcd_degree']} steps={result['euclid_steps']} "
        f"in {result['seconds_wall_reference_only']} s",
        flush=True,
    )
    del rows
    gc.collect()
    return {
        "N0": N0,
        "G": G,
        "D": D,
        "epsilon": epsilon,
        "e_min": e_min,
        "e_max": e_max,
        "energy_min": energy_min,
        "energy_max": energy_max,
        "E": E,
        "F": F,
        "build_seconds": build_seconds,
        "per_prime": {str(prime): result},
    }


def main() -> int:
    check("producer is not imported by this standalone verifier", "e136_pair_product_2x5" not in sys.modules)
    check("both declared primes are prime", all(is_prime(p) for p in PRIMES))
    check(
        "primes exceed the degree 523776 and avoid 3 and 5",
        all(p > 523776 and p % 3 and p % 5 for p in PRIMES),
    )
    check(
        "strict-gcd stage wall is a local process-time clock",
        local_gcd_wall_probe(),
    )
    check("strict long division matches the pure-Python reference", strict_ops_self_validation(PRIMES[0]))

    # Independent sign convention control before any expensive matrix work.
    sign = newton(np.array([0, 10, 38, 160], dtype=np.int64), 3, PRIMES[0])
    check(
        "Newton signs recover (z-2)(z-3)(z-5)",
        [int(value) for value in sign] == [1, -10 % PRIMES[0], 31 % PRIMES[0], -30 % PRIMES[0]],
    )

    bonds23 = list(layer_bonds((2, 3), (False, False)))
    bonds33 = list(layer_bonds((3, 3), (False, False)))
    bonds25 = list(layer_bonds((2, 5), (False, False)))
    bonds_n10 = [(site, site + 1) for site in range(9)]
    degree25 = sorted(sum(site in bond for bond in bonds25) for site in range(10))
    check(
        "2x5 target is the open 13-bond layer with degree profile 2^4 3^6",
        len(bonds25) == 13 and degree25 == [2, 2, 2, 2, 3, 3, 3, 3, 3, 3],
        f"bonds={len(bonds25)}, degrees={degree25}",
    )
    check("n=10 chain control has exactly nine consecutive bonds", bonds_n10 == [(i, i + 1) for i in range(9)])

    # Every claimed finite row is independently rebuilt from raw inputs.
    cases = {
        "layer_2x3_t_1_3": run_fraction_case("layer_2x3_t_1_3", 6, bonds23, PRIMES),
        "layer_3x3_t_1_3": run_fraction_case("layer_3x3_t_1_3", 9, bonds33, (PRIMES[0],)),
        "chain_n10_t_1_3": run_integer_case("chain_n10_t_1_3", 10, bonds_n10, PRIMES[0]),
        "layer_2x5_t_1_3": run_integer_case("layer_2x5_t_1_3", 10, bonds25, PRIMES[0]),
    }
    expected = {
        "layer_2x3_t_1_3": {"n": 6, "gcd": 385, "pairs": 2016, "max": 665, "energy": (-7, 7)},
        "layer_3x3_t_1_3": {"n": 9, "gcd": 59641, "pairs": 130816, "max": 19171, "energy": (-12, 12)},
        "chain_n10_t_1_3": {"n": 10, "gcd": 465751, "pairs": 523776, "max": 58025, "energy": (-9, 9)},
        "layer_2x5_t_1_3": {"n": 10, "gcd": None, "pairs": 523776, "max": 58025, "energy": (-13, 13)},
    }
    for name, case in cases.items():
        want = expected[name]
        energy = (case["energy_min"], case["energy_max"])
        check(
            f"{name}: raw energy range catches geometry/direction defects",
            energy == want["energy"],
            f"got {energy}",
        )
        for p_str, result in case["per_prime"].items():
            check(
                f"{name}@{p_str}: independent monicity, float/int, Horner, endpoint recurrence, strict gcd division",
                bool(result["charpoly_monic"])
                and bool(result["pair_poly_monic"])
                and bool(result["first_product_equal"])
                and bool(result["horner_zero"])
                and all(result["recurrence_spots"].values())
                and bool(result["gcd_divides_both"]),
            )
            if want["gcd"] is not None:
                check(
                    f"{name}@{p_str}: independently reconstructed gcd degree",
                    int(result["gcd_degree"]) == want["gcd"],
                    f"got {result['gcd_degree']}",
                )

    # The decisive inference and the same-size Gaussian control are deliberately
    # separate assertions: an off-by-one in C(1024,2) or a missing distinct-slot
    # convention makes at least one fail.
    pair_slots = 1024 * 1023 // 2
    gaussian_max = 3**10 - 2**10
    gaussian_floor = pair_slots - gaussian_max
    layer_gcd = int(cases["layer_2x5_t_1_3"]["per_prime"][str(PRIMES[0])]["gcd_degree"])
    chain_gcd = int(cases["chain_n10_t_1_3"]["per_prime"][str(PRIMES[0])]["gcd_degree"])
    check(
        "ten-site arithmetic is 523776 slots, 58025 Gaussian maximum, floor 465751",
        (pair_slots, gaussian_max, gaussian_floor) == (523776, 58025, 465751),
    )
    check(
        "DECISIVE: raw 2x5 reconstruction gives gcd < 465751",
        layer_gcd < gaussian_floor,
        f"gcd={layer_gcd}",
    )
    check(
        "DECISIVE: raw 2x5 reconstruction gives > 58025 distinct pair products",
        pair_slots - layer_gcd > gaussian_max,
        f"distinct >= {pair_slots - layer_gcd}",
    )
    check(
        "same-size open-chain control reaches exactly the Gaussian floor",
        chain_gcd == gaussian_floor,
        f"gcd={chain_gcd}",
    )

    # Independent check of the generic-t bookkeeping arithmetic for 2x5:
    # c = max(-e_min, e_max) = 7 over the measured range [-7, 6].
    e_min = int(cases["layer_2x5_t_1_3"]["e_min"])
    e_max = int(cases["layer_2x5_t_1_3"]["e_max"])
    c = max(-e_min, e_max)
    scalar_degree = 3 * c + e_max
    entry_degree = 2 * 10 + scalar_degree
    check(
        "generic-t clearing exponent and degree bounds: c=7, scalar 27, entry 47, wedge 94",
        c == 7 and scalar_degree == 27 and entry_degree == 47 and 2 * entry_degree == 94,
        f"c={c}, scalar={scalar_degree}, entry={entry_degree}",
    )
    resultant_bound = 2 * (2 * entry_degree) * pair_slots * (pair_slots - 1)
    check(
        "crude exceptional resultant bound equals 2*94*523776*523775",
        resultant_bound == 2 * 94 * 523776 * 523775,
        f"{resultant_bound}",
    )

    # Read the producer artifact only after independent raw reconstruction.
    artifact = json.loads(ARTIFACT.read_text())
    check(
        "artifact producer check ledger is nonempty and all recorded checks passed",
        bool(artifact["checks"]) and all(record.get("ok") is True for record in artifact["checks"]),
    )
    stored_cases = artifact["data"]["cases"]
    for name, case in cases.items():
        want = expected[name]
        stored = stored_cases[name]
        check(
            f"artifact {name}: no stage was blocked and both primes completed",
            stored.get("blocked_stage") is None and len(stored.get("per_prime", {})) == 2,
        )
        check(
            f"artifact {name}: actual LCM agrees with the raw rebuild",
            int(stored["D"]) == case["D"],
        )
        n = want["n"]
        dim = 1 << n
        pair_count = dim * (dim - 1) // 2
        gaussian_maximum = 3**n - 2**n
        gaussian_floor_case = pair_count - gaussian_maximum
        raw_gcds = [int(result["gcd_degree"]) for result in case["per_prime"].values()]
        check(
            f"artifact {name}: every stored threshold and decisive-inference field matches the raw rebuild",
            int(stored["n_sites"]) == n
            and int(stored["dim"]) == dim
            and int(stored["pair_count"]) == pair_count
            and int(stored["gaussian_max_distinct_pair_products"]) == gaussian_maximum
            and int(stored["gaussian_gcd_floor"]) == gaussian_floor_case
            and int(stored["modular_gcd_degree_upper_bound"]) == min(
                int(value) for value in stored["modular_gcd_degrees"]
            )
            and int(stored["distinct_pair_products_lower_bound"])
            == pair_count - int(stored["modular_gcd_degree_upper_bound"])
            and all(
                int(stored["modular_gcd_degrees"][index]) == raw_gcds[index]
                for index in range(len(raw_gcds))
            ),
        )
        build_resource = stored["build_resource"]
        check(
            f"artifact {name}: observed build CPU/RSS recorded under its predeclared wall",
            build_resource["predeclared_wall_cpu_seconds"] == 900
            and 0 <= build_resource["observed_cpu_seconds"] <= build_resource["predeclared_wall_cpu_seconds"]
            and build_resource["ru_maxrss_bytes_after_stage"] < RSS_LIMIT_BYTES,
        )
        for p_str, result in case["per_prime"].items():
            saved = stored["per_prime"][p_str]
            check(
                f"artifact {name}@{p_str}: all independent coefficient digests agree",
                result["charpoly_digest"] == saved["charpoly_sha256"]
                and result["pair_poly_digest"] == saved["pair_poly_sha256"]
                and result["gcd_digest"] == saved["gcd_sha256"],
            )
            stage_keys = {"traces", "charpoly", "recurrence", "pair_newton", "gcd"}
            stages = saved["stage_resources"]
            check(
                f"artifact {name}@{p_str}: every predeclared pipeline stage records observed CPU time and RSS",
                set(stages) == stage_keys
                and saved["matrix_mod_p_seconds"] >= 0
                and saved["matrix_mod_p_ru_maxrss_bytes_after_construction"] < RSS_LIMIT_BYTES
                and all(
                    stages[stage]["predeclared_wall_cpu_seconds"]
                    == artifact["data"]["resource_policy"]["predeclared_stage_wall_cpu_seconds"][stage]
                    and 0
                    <= stages[stage]["observed_cpu_seconds"]
                    <= stages[stage]["predeclared_wall_cpu_seconds"]
                    and stages[stage]["ru_maxrss_bytes_after_stage"] < RSS_LIMIT_BYTES
                    for stage in stage_keys
                ),
            )

    # The artifact's decisive row and generic-t bookkeeping against the raw values.
    stored25 = stored_cases["layer_2x5_t_1_3"]
    check(
        "artifact 2x5 row is decisive against the raw rebuild",
        bool(stored25["decisive_non_gaussian"])
        and not bool(stored25["at_gaussian_floor"])
        and int(stored25["modular_gcd_degree_upper_bound"]) == layer_gcd,
        f"stored gcd bound {stored25['modular_gcd_degree_upper_bound']}, raw {layer_gcd}",
    )
    generic = artifact["data"]["generic_coupling_theorem_2x5"]
    check(
        "artifact generic-t block: c=7, entry 47, wedge 94, distinct lower bound above 58025",
        generic["clearing_exponent_c"] == 7
        and generic["matrix_entry_t_degree_upper_bound"] == 47
        and generic["wedge2_entry_t_degree_upper_bound"] == 94
        and generic["resultant_exceptional_root_bound"] == resultant_bound
        and generic["generic_gcd_degree_upper_bound"] == int(stored25["modular_gcd_degree_upper_bound"])
        and generic["generic_distinct_pair_products_lower_bound"] == pair_slots - layer_gcd
        and generic["generic_distinct_pair_products_lower_bound"] > generic["gaussian_max_distinct_pair_products"],
    )
    check(
        "artifact resource policy gates on process_time with per-case abort semantics",
        artifact["data"]["resource_policy"]["gate"].startswith("every compute gate reads time.process_time()")
        and artifact["data"]["resource_policy"]["rss_limit_bytes"] == 6_000_000_000,
    )

    # The wave-13 digest replay claim: the artifact's 3x3 row equals the stored
    # wave-13 digests at both primes, independently re-read here.
    if WAVE13_ARTIFACT.exists():
        wave13 = json.loads(WAVE13_ARTIFACT.read_text())["data"]["cases"]["layer_3x3_t_1_3"]
        replay_ok = True
        for p in PRIMES:
            saved = stored_cases["layer_3x3_t_1_3"]["per_prime"][str(p)]
            replay_ok = replay_ok and (
                saved["charpoly_sha256"] == wave13["per_prime"][str(p)]["charpoly_sha256"]
                and saved["pair_poly_sha256"] == wave13["per_prime"][str(p)]["pair_poly_sha256"]
                and saved["gcd_sha256"] == wave13["per_prime"][str(p)]["gcd_sha256"]
            )
        check("artifact 3x3 row digests replay the stored wave-13 digests at both primes", replay_ok)
        check(
            "raw 3x3 rebuild at p1 matches the stored wave-13 digest chain",
            cases["layer_3x3_t_1_3"]["per_prime"][str(PRIMES[0])]["charpoly_digest"]
            == wave13["per_prime"][str(PRIMES[0])]["charpoly_sha256"]
            and cases["layer_3x3_t_1_3"]["per_prime"][str(PRIMES[0])]["gcd_digest"]
            == wave13["per_prime"][str(PRIMES[0])]["gcd_sha256"],
        )

    elapsed_wall = time.perf_counter() - STARTED_WALL
    elapsed_cpu = time.process_time() - STARTED_CPU
    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(
        "PASS: independently rebuilt raw 2x5 certificate, chain floor control, 3x3 digest replay, "
        f"and 2x3 regression in {elapsed_wall:.2f} s wall / {elapsed_cpu:.2f} CPU s; "
        f"peak RSS {max_rss_bytes()} bytes."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
