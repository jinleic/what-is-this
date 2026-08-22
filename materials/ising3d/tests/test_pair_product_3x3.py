"""Standalone independent verifier for the open-3x3 pair-product certificate.

This test deliberately imports neither experiments/e122_pair_product_3x3.py nor
e38_gaussianity_certificate.py.  It independently:

  1. re-derives the canonical rational R=P_t diag(q^e) P_t using Fraction,
  2. computes the ACTUAL entry-denominator LCM D directly from every R entry,
  3. reduces M=D*R modulo each declared prime,
  4. reconstructs charpoly(M), pair-product charpoly(wedge^2 M), and its gcd
     with the derivative, and
  5. compares only after computation against the stored artifact digests.

The producer uses a fast split-integral construction and lazy-reduction Euclid.
This verifier instead starts from the raw Fraction matrix and uses strict modular
reduction at every long-division update.  It covers both primes and every finite
row claimed in the artifact: the decisive 3x3 layer, the n=9 chain control, and
the wave-12 2x3/2x4 layer regressions.

Run: PYTHONPATH=src .venv/bin/python tests/test_pair_product_3x3.py
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
from math import isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402

INT64_MAX = 2**63 - 1
FLOAT_EXACT_MAX = 2**53
PRIMES = (1_000_003, 2_000_003)
RSS_LIMIT_BYTES = 7_500_000_000
TEST_WALL_SECONDS = 3_600       # CPU-seconds budget (process_time), not wall time
GCD_WALL_SECONDS = 900          # CPU-seconds budget per gcd stage
ARTIFACT = ROOT / "results" / "spectral" / "pair_product_3x3.json"

FAILS: list[str] = []
STARTED = time.process_time()


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)
    return bool(ok)


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def enforce(stage: str, wall_seconds: int = TEST_WALL_SECONDS, *, started: float | None = None) -> None:
    # CPU time, not wall time: the machine routinely carries >5x unrelated
    # background load, which must not change this verifier's verdict.
    anchor = STARTED if started is None else started
    elapsed = time.process_time() - anchor
    assert elapsed < wall_seconds, f"budget {wall_seconds}s CPU fired at {stage} after {elapsed:.3f}s"
    assert max_rss_bytes() < RSS_LIMIT_BYTES, f"RSS wall at {stage}: {max_rss_bytes()} >= {RSS_LIMIT_BYTES}"


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


# --------------------------------------------------- independent raw Fraction build

def raw_build_r(n: int, bonds: list[tuple[int, int]], t: Fraction):
    """Re-derive R=P_t diag(q^e) P_t from the raw rational construction.

    This does not call e38/e122.  The j>=i loop preserves the canonical symmetric
    matrix convention but every entry is evaluated by its own exact Fraction sum.
    """
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    energy = []
    for state in range(dim):
        spins = [1 - 2 * ((state >> (n - 1 - site)) & 1) for site in range(n)]
        energy.append(sum(spins[left] * spins[right] for left, right in bonds))
    epsilon = energy[0] % 2
    assert all((value - epsilon) % 2 == 0 for value in energy)
    weights = [q ** ((value - epsilon) // 2) for value in energy]
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
    return matrix, q, epsilon, min(energy), max(energy)


def actual_integral_lcm(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    """Actual LCM of all entry denominators, followed by exact M=D*R."""
    scale = 1
    for row in matrix:
        for entry in row:
            scale = scale * entry.denominator // math_gcd(scale, entry.denominator)
    integer_rows: list[list[int]] = []
    for row in matrix:
        output = []
        for entry in row:
            value = scale * entry
            assert value.denominator == 1
            output.append(value.numerator)
        integer_rows.append(output)
    return scale, integer_rows


# ------------------------------------------------ exact modular arithmetic / traces

def assert_int64_guard(length: int, p: int, where: str) -> None:
    assert length * (p - 1) ** 2 < INT64_MAX, f"int64 guard failed at {where}"


def assert_float_guard(length: int, p: int, where: str) -> None:
    assert length * (p - 1) ** 2 < FLOAT_EXACT_MAX, f"binary64 exactness guard failed at {where}"


def float_matmul_mod(left: np.ndarray, right: np.ndarray, p: int) -> np.ndarray:
    """Exact modular product through binary64, independently audited by round-trip equality."""
    assert left.shape[1] == right.shape[0]
    assert_float_guard(left.shape[1], p, "float matrix product")
    raw = left.astype(np.float64) @ right.astype(np.float64)
    integer = raw.astype(np.int64)
    # The cast round-trip is a concrete runtime check that each BLAS output was
    # an exactly represented integer, in addition to the 2^53 proof bound.
    assert np.array_equal(raw, integer.astype(np.float64)), "non-integral binary64 matrix product"
    return integer % p


def int_matmul_mod(left: np.ndarray, right: np.ndarray, p: int) -> np.ndarray:
    assert_int64_guard(left.shape[1], p, "int64 matrix product")
    return (left @ right) % p


def power_traces(matrix: np.ndarray, p: int) -> np.ndarray:
    dim = matrix.shape[0]
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
    """[1,a1,...,a_degree] for z^degree+a1 z^(degree-1)+..., independently implemented."""
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
    """Independent Cayley-Hamilton check, retaining exact binary64 integer guards."""
    dim = matrix.shape[0]
    value = np.eye(dim, dtype=np.int64)
    identity = np.eye(dim, dtype=np.int64)
    for coefficient in characteristic[1:]:
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


def strict_remainder(dividend: np.ndarray, divisor: np.ndarray, p: int, gcd_started: float) -> np.ndarray:
    """Fully reduced long division; unlike the producer, reduces each update mod p."""
    divisor = monic(divisor, p)
    assert not zero(divisor)
    work = np.asarray(dividend, dtype=np.int64).copy() % p
    dividend_degree = work.size - 1
    divisor_degree = divisor.size - 1
    if dividend_degree < divisor_degree:
        return trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset & 1023 == 0:
            enforce("strict Euclidean division", GCD_WALL_SECONDS, started=gcd_started)
        lead = int(work[offset])
        if lead:
            work[offset : offset + divisor_degree + 1] = (
                work[offset : offset + divisor_degree + 1] - lead * divisor
            ) % p
    return trim(work[dividend_degree - divisor_degree + 1 :])


def strict_gcd(first: np.ndarray, second: np.ndarray, p: int) -> tuple[np.ndarray, int]:
    gcd_started = time.process_time()
    first, second = monic(first, p), monic(second, p)
    steps = 0
    while not zero(second):
        if steps & 127 == 0:
            enforce("strict Euclidean gcd", GCD_WALL_SECONDS, started=gcd_started)
        first, second = second, monic(strict_remainder(first, second, p, gcd_started), p)
        steps += 1
    return first, steps

def local_gcd_wall_probe() -> bool:
    """A strict-gcd stage must use its own clock, not prior raw-build elapsed time."""
    global STARTED
    original_started = STARTED
    try:
        STARTED = time.process_time() - GCD_WALL_SECONDS - 1
        strict_gcd(
            np.array([1, PRIMES[0] - 2, 1], dtype=np.int64),
            np.array([1, PRIMES[0] - 1], dtype=np.int64),
            PRIMES[0],
        )
        return True
    except AssertionError:
        return False
    finally:
        STARTED = original_started


def derivative_descending(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    return np.ascontiguousarray((np.arange(degree, 0, -1, dtype=np.int64) * poly[:-1]) % p)


def digest(ascending) -> str:
    return hashlib.sha256(",".join(str(int(value)) for value in ascending).encode()).hexdigest()


# ---------------------------------------------------------------- independent full pipeline

def independent_pipeline(integer_rows: list[list[int]], p: int) -> dict[str, object]:
    started = time.perf_counter()
    dim = len(integer_rows)
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > 2 * pairs and p % 2 == 1
    matrix = np.array([[entry % p for entry in row] for row in integer_rows], dtype=np.int64)

    first_product_equal = np.array_equal(float_matmul_mod(matrix, matrix, p), int_matmul_mod(matrix, matrix, p))
    traces = power_traces(matrix, p)
    characteristic = newton(traces, dim, p)
    horner_ok = horner_zero(characteristic, matrix, p)
    extended = extend_traces(traces, characteristic, 2 * pairs, p)

    # Direct int64 powers hit an immediate post-charpoly index, the pair-slot
    # endpoint, and the full 2*pairs endpoint: this catches a recurrence sign,
    # a dropped factor of two, and an off-by-one range defect.
    spot_indices = tuple(sorted({dim + 1, pairs, 2 * pairs}))
    spots = {k: binary_trace_int64(matrix, k, p) == int(extended[k]) for k in spot_indices}

    inverse_two = pow(2, -1, p)
    pair_sums = np.zeros(pairs + 1, dtype=np.int64)
    pair_sums[1:] = (
        (extended[1 : pairs + 1] * extended[1 : pairs + 1] - extended[2 : 2 * pairs + 1 : 2]) % p
    ) * inverse_two % p
    pair_poly = newton(pair_sums, pairs, p)
    derivative = derivative_descending(pair_poly, p)
    common, steps = strict_gcd(pair_poly, derivative, p)
    division_check_started = time.perf_counter()
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
        "seconds": round(time.perf_counter() - started, 3),
        "rss_bytes": max_rss_bytes(),
    }


def run_raw_case(name: str, n: int, bonds: list[tuple[int, int]], expected_d: int) -> dict[str, object]:
    print(f"raw Fraction reconstruction {name} ...", flush=True)
    started = time.perf_counter()
    rational, q, epsilon, energy_min, energy_max = raw_build_r(n, bonds, Fraction(1, 3))
    D, rows = actual_integral_lcm(rational)
    del rational
    gc.collect()
    build_seconds = round(time.perf_counter() - started, 3)
    check(f"{name}: q=5/3 and canonical parity", q == Fraction(5, 3) and epsilon in {0, 1})
    check(f"{name}: actual raw Fraction denominator LCM", D == expected_d, f"got D={D}")
    results = {}
    for p in PRIMES:
        print(f"  independent strict pipeline p={p} ...", flush=True)
        results[str(p)] = independent_pipeline(rows, p)
        print(
            f"    gcd={results[str(p)]['gcd_degree']} steps={results[str(p)]['euclid_steps']} "
            f"in {results[str(p)]['seconds']} s",
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


def main() -> int:
    check("producer is not imported by this standalone verifier", "e122_pair_product_3x3" not in sys.modules)
    check("both declared primes are prime", all(is_prime(p) for p in PRIMES))
    check("strict-gcd stage wall has a local clock independent of prior raw-build work", local_gcd_wall_probe())
    check("p exceeds 2*C(512,2) and avoids 3,5", all(p > 261632 and p % 3 and p % 5 for p in PRIMES))

    # Independent sign convention control before any expensive matrix work.
    sign = newton(np.array([0, 10, 38, 160], dtype=np.int64), 3, PRIMES[0])
    check(
        "Newton signs recover (z-2)(z-3)(z-5)",
        [int(value) for value in sign] == [1, -10 % PRIMES[0], 31 % PRIMES[0], -30 % PRIMES[0]],
    )

    bonds23 = list(layer_bonds((2, 3), (False, False)))
    bonds24 = list(layer_bonds((2, 4), (False, False)))
    bonds33 = list(layer_bonds((3, 3), (False, False)))
    bonds_n9 = [(site, site + 1) for site in range(8)]
    degree33 = [sum(site in bond for bond in bonds33) for site in range(9)]
    check(
        "3x3 target is the open 12-bond layer with a degree-4 center",
        len(bonds33) == 12 and degree33[4] == 4 and sorted(degree33) == [2, 2, 2, 2, 3, 3, 3, 3, 4],
        f"bonds={len(bonds33)}, degrees={degree33}",
    )
    check("n=9 chain control has exactly eight consecutive bonds", bonds_n9 == [(i, i + 1) for i in range(8)])

    # Every claimed finite row is independently rebuilt from raw inputs.
    cases = {
        "layer_2x3_t_1_3": run_raw_case("layer_2x3_t_1_3", 6, bonds23, 3**15 * 5**4),
        "layer_2x4_t_1_3": run_raw_case("layer_2x4_t_1_3", 8, bonds24, 3**21 * 5**5),
        "layer_3x3_t_1_3": run_raw_case("layer_3x3_t_1_3", 9, bonds33, 3**24 * 5**6),
        "chain_n9_t_1_3": run_raw_case("chain_n9_t_1_3", 9, bonds_n9, 3**22 * 5**4),
    }
    expected = {
        "layer_2x3_t_1_3": {"n": 6, "gcd": 385, "pairs": 2016, "max": 665, "energy": (-7, 7)},
        "layer_2x4_t_1_3": {"n": 8, "gcd": 9329, "pairs": 32640, "max": 6305, "energy": (-10, 10)},
        "layer_3x3_t_1_3": {"n": 9, "gcd": 59641, "pairs": 130816, "max": 19171, "energy": (-12, 12)},
        "chain_n9_t_1_3": {"n": 9, "gcd": 111645, "pairs": 130816, "max": 19171, "energy": (-8, 8)},
    }
    for name, case in cases.items():
        want = expected[name]
        check(f"{name}: raw energy range catches geometry/direction defects", (case["energy_min"], case["energy_max"]) == want["energy"])
        for p in PRIMES:
            result = case["per_prime"][str(p)]
            check(
                f"{name}@{p}: independent monicity, float/int, Horner, endpoint recurrence, and strict gcd division",
                bool(result["charpoly_monic"])
                and bool(result["pair_poly_monic"])
                and bool(result["first_product_equal"])
                and bool(result["horner_zero"])
                and all(result["recurrence_spots"].values())
                and bool(result["gcd_divides_both"]),
            )
            check(
                f"{name}@{p}: independently reconstructed gcd degree",
                int(result["gcd_degree"]) == want["gcd"],
                f"got {result['gcd_degree']}",
            )

    # The decisive inference and the same-size Gaussian control are deliberately
    # separate assertions: an off-by-one in C(512,2) or a missing distinct-slot
    # convention makes at least one fail.
    pair_slots = 512 * 511 // 2
    gaussian_max = 3**9 - 2**9
    gaussian_floor = pair_slots - gaussian_max
    layer_gcd = cases["layer_3x3_t_1_3"]["per_prime"][str(PRIMES[0])]["gcd_degree"]
    chain_gcds = [cases["chain_n9_t_1_3"]["per_prime"][str(p)]["gcd_degree"] for p in PRIMES]
    check("nine-site arithmetic is 130816 slots, 19171 Gaussian maximum, floor 111645", (pair_slots, gaussian_max, gaussian_floor) == (130816, 19171, 111645))
    check("DECISIVE: raw 3x3 reconstruction gives gcd < 111645", layer_gcd < gaussian_floor, f"gcd={layer_gcd}")
    check("DECISIVE: raw 3x3 reconstruction gives >=71175 distinct pair products", pair_slots - layer_gcd == 71175 and pair_slots - layer_gcd > gaussian_max)
    check("same-size open-chain control reaches exactly the Gaussian floor at both primes", chain_gcds == [gaussian_floor, gaussian_floor], f"gcds={chain_gcds}")

    # Read the producer artifact only after independent raw reconstruction, and
    # compare all stored decisive values/digests row by row and prime by prime.
    artifact = json.loads(ARTIFACT.read_text())
    check(
        "artifact producer check ledger is nonempty and all recorded checks passed",
        bool(artifact["checks"]) and all(record.get("ok") is True for record in artifact["checks"]),
    )
    stored_cases = artifact["data"]["cases"]
    for name, case in cases.items():
        want = expected[name]
        stored = stored_cases[name]
        check(f"artifact {name}: actual LCM agrees with raw rebuild", int(stored["D"]) == case["D"])
        n = want["n"]
        dim = 1 << n
        pair_count = dim * (dim - 1) // 2
        gaussian_maximum = 3**n - 2**n
        gaussian_floor_case = pair_count - gaussian_maximum
        raw_gcds = [int(case["per_prime"][str(p)]["gcd_degree"]) for p in PRIMES]
        raw_upper_bound = min(raw_gcds)
        raw_distinct_lower_bound = pair_count - raw_upper_bound
        raw_decisive = raw_upper_bound < gaussian_floor_case
        check(
            f"artifact {name}: every stored threshold and decisive-inference field matches the raw rebuild",
            int(stored["n_sites"]) == n
            and int(stored["dim"]) == dim
            and int(stored["pair_count"]) == pair_count
            and int(stored["gaussian_max_distinct_pair_products"]) == gaussian_maximum
            and int(stored["gaussian_gcd_floor"]) == gaussian_floor_case
            and [int(value) for value in stored["modular_gcd_degrees"]] == raw_gcds
            and int(stored["modular_gcd_degree_upper_bound"]) == raw_upper_bound
            and int(stored["distinct_pair_products_lower_bound"]) == raw_distinct_lower_bound
            and bool(stored["decisive_non_gaussian"]) == raw_decisive
            and bool(stored["at_gaussian_floor"]) == (raw_upper_bound == gaussian_floor_case),
        )
        build_resource = stored["build_resource"]
        check(
            f"artifact {name}: observed build RSS/time is recorded under its predeclared wall",
            build_resource["predeclared_wall_seconds"] == 300
            and 0 <= build_resource["observed_seconds"] <= build_resource["predeclared_wall_seconds"]
            and build_resource["ru_maxrss_bytes_after_stage"] < RSS_LIMIT_BYTES,
        )
        for p in PRIMES:
            got = case["per_prime"][str(p)]
            saved = stored["per_prime"][str(p)]
            check(
                f"artifact {name}@{p}: all independent coefficient digests agree",
                got["charpoly_digest"] == saved["charpoly_sha256"]
                and got["pair_poly_digest"] == saved["pair_poly_sha256"]
                and got["gcd_digest"] == saved["gcd_sha256"],
            )
            timing_keys = {
                "traces": "traces",
                "charpoly": "characteristic_newton_and_horner",
                "recurrence": "trace_recurrence_and_spots",
                "pair_newton": "pair_newton",
                "gcd": "euclidean_gcd",
            }
            expected_walls = {"traces": 240, "charpoly": 120, "recurrence": 240, "pair_newton": 600, "gcd": 1200}
            stages = saved["stage_resources"]
            check(
                f"artifact {name}@{p}: every predeclared pipeline stage records observed time and RSS",
                set(stages) == set(timing_keys)
                and saved["matrix_mod_p_seconds"] >= 0
                and saved["matrix_mod_p_ru_maxrss_bytes_after_construction"] < RSS_LIMIT_BYTES
                and all(
                    stages[stage]["predeclared_wall_seconds"] == expected_walls[stage]
                    and stages[stage]["observed_seconds"] == saved["timing_seconds"][timing_key]
                    and 0 <= stages[stage]["observed_seconds"] <= stages[stage]["predeclared_wall_seconds"]
                    and stages[stage]["ru_maxrss_bytes_after_stage"] < RSS_LIMIT_BYTES
                    for stage, timing_key in timing_keys.items()
                ),
            )

    elapsed = time.process_time() - STARTED
    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(
        "PASS: independently rebuilt raw-Fraction 3x3 certificate, chain control, and wave-12 regressions "
        f"in {elapsed:.2f} s; peak RSS {max_rss_bytes()} bytes."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
