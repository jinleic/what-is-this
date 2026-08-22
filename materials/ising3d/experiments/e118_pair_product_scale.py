"""Scale the exact ordering-free pair-product obstruction.

For an exact rational transfer representative R, this script clears *all*
entry denominators to M = D R, then reconstructs
C_2(z) = charpoly(wedge^2 M) modulo two primes from power traces.  Because
C_2 is monic integral, a modular gcd degree is an upper bound on the
characteristic-zero gcd degree; it is never treated as an equality over Q.

The implementation deliberately follows the wave-11 route:
tr(M^k) -> Newton charpoly -> Cayley-Hamilton trace recurrence -> pair traces
-> Newton C_2 -> explicit Euclidean gcd(C_2, C_2').  There is no eigenvalue
ordering or floating-point calculation in a certificate.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import resource
import signal
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd as math_gcd
from math import isqrt
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from ising.transfer_matrix import layer_bonds  # noqa: E402
from e38_gaussianity_certificate import build_R  # noqa: E402 -- canonical rational construction

INT64_MAX = 2**63 - 1
PRIMES = (1_000_003, 2_000_003)
RSS_LIMIT_BYTES = 7_500_000_000
CASE_WALL_SECONDS = 1_800
THREE_BY_THREE_WALL_SECONDS = int(os.environ.get("PAIR_PRODUCT_3X3_WALL_SECONDS", "90"))

# Predeclared before any case is computed.  The additional regular grid is
# k/20, k=4,...,10.  Unioning it with the prior-coupling set yields these eight
# distinct exact values in increasing order.
REPOSITORY_COUPLINGS = (Fraction(1, 5), Fraction(1, 3), Fraction(2, 5), Fraction(1, 2))
GRID_DENOMINATOR = 20
GRID_NUMERATORS = tuple(range(4, 11))
GRID_COUPLINGS = tuple(Fraction(k, GRID_DENOMINATOR) for k in GRID_NUMERATORS)
COUPLINGS_2X3 = tuple(sorted(set(REPOSITORY_COUPLINGS + GRID_COUPLINGS)))

CHECKS: list[dict[str, object]] = []
FAILS: list[str] = []


class ResourceWall(RuntimeError):
    """A predeclared time or RSS wall stopped an intentionally bounded route."""


@dataclass
class Progress:
    stage: str = "initialization"
    started: float = 0.0
    wall_seconds: int = CASE_WALL_SECONDS


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    CHECKS.append({"name": name, "passed": bool(ok), "detail": detail})
    if not ok:
        FAILS.append(name)
    return bool(ok)


def max_rss_bytes() -> int:
    """Return ru_maxrss in bytes on Darwin and Linux."""
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def assert_resources(progress: Progress) -> None:
    elapsed = time.perf_counter() - progress.started
    if elapsed > progress.wall_seconds:
        raise ResourceWall(f"time wall {progress.wall_seconds}s at {progress.stage} (elapsed {elapsed:.2f}s)")
    rss = max_rss_bytes()
    if rss > RSS_LIMIT_BYTES:
        raise ResourceWall(f"RSS wall {RSS_LIMIT_BYTES} bytes at {progress.stage} (ru_maxrss {rss})")


@contextmanager
def alarm_wall(seconds: int):
    """Interrupt Python-level exact work at a real-time wall on POSIX platforms."""
    if seconds <= 0:
        yield
        return

    def handler(_signum, _frame):
        raise ResourceWall(f"signal time wall {seconds}s")

    previous = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


# --------------------------------------------------------------------------- arithmetic

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for divisor in range(2, isqrt(n) + 1):
        if n % divisor == 0:
            return False
    return True


def assert_int64_guard(length: int, p: int, kind: str) -> None:
    assert length * (p - 1) ** 2 < INT64_MAX, f"int64 {kind} guard violated"


def sequential_power_traces(matrix: np.ndarray, max_power: int, p: int, progress: Progress) -> list[int]:
    """[0,tr(M),...,tr(M^max_power)] by direct modular matrix products."""
    dim = matrix.shape[0]
    assert_int64_guard(dim, p, "matrix multiplication")
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for exponent in range(1, max_power + 1):
        if exponent & 15 == 0:
            assert_resources(progress)
        power = (power @ matrix) % p
        traces.append(int(power.trace()) % p)
    return traces


def trace_power_binary(matrix: np.ndarray, exponent: int, p: int) -> int:
    """Independent binary-power trace used only to audit the recurrence."""
    dim = matrix.shape[0]
    assert_int64_guard(dim, p, "matrix multiplication")
    result = np.eye(dim, dtype=np.int64)
    base = matrix.copy()
    while exponent:
        if exponent & 1:
            result = (result @ base) % p
        exponent >>= 1
        if exponent:
            base = (base @ base) % p
    return int(result.trace()) % p


def newton_from_power_sums(
    power_sums: np.ndarray, degree: int, p: int, inverses: list[int], progress: Progress
) -> np.ndarray:
    """Return DESCENDING [1,a1,...,a_degree] from power sums via Newton identities."""
    assert_int64_guard(degree, p, "Newton dot product")
    coeffs = np.zeros(degree + 1, dtype=np.int64)
    coeffs[0] = 1
    for k in range(1, degree + 1):
        if k & 255 == 0:
            assert_resources(progress)
        total = int(np.dot(power_sums[1 : k + 1], np.ascontiguousarray(coeffs[k - 1 :: -1]))) % p
        coeffs[k] = (-total * inverses[k]) % p
    return coeffs


def extend_power_sums(
    traces: list[int], charpoly_desc: np.ndarray, upto: int, p: int, progress: Progress
) -> list[int]:
    """Apply the exact Cayley-Hamilton trace recurrence beyond matrix dimension."""
    dim = charpoly_desc.size - 1
    assert_int64_guard(dim, p, "trace recurrence")
    coefficients = charpoly_desc[1:]
    for k in range(dim + 1, upto + 1):
        if k & 255 == 0:
            assert_resources(progress)
        window = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
        total = int(np.dot(coefficients, window)) % p
        traces.append((-total) % p)
    return traces


def horner_matrix_zero(charpoly_desc: np.ndarray, matrix: np.ndarray, p: int) -> bool:
    """Cayley-Hamilton validation in F_p^{dim x dim}."""
    dim = matrix.shape[0]
    assert_int64_guard(dim, p, "matrix multiplication")
    value = int(charpoly_desc[0]) * np.eye(dim, dtype=np.int64) % p
    identity = np.eye(dim, dtype=np.int64)
    for coefficient in charpoly_desc[1:]:
        value = ((value @ matrix) + int(coefficient) * identity) % p
    return bool(np.all(value == 0))


# --------------------------------------------------------------- F_p[z], descending

def poly_trim(poly: np.ndarray) -> np.ndarray:
    nonzero = np.flatnonzero(poly)
    if nonzero.size == 0:
        return np.zeros(1, dtype=np.int64)
    return poly[nonzero[0] :]


def poly_is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def poly_monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = poly_trim(np.asarray(poly, dtype=np.int64) % p)
    if poly_is_zero(poly):
        return poly
    return (poly * pow(int(poly[0]), -1, p)) % p


def poly_remainder(dividend: np.ndarray, divisor: np.ndarray, p: int, progress: Progress) -> np.ndarray:
    """Remainder by monic long division, all arrays in descending order."""
    assert int(divisor[0]) == 1
    work = (np.asarray(dividend, dtype=np.int64) % p).copy()
    divisor_degree = divisor.size - 1
    dividend_degree = work.size - 1
    if dividend_degree < divisor_degree:
        return poly_trim(work)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset & 255 == 0:
            assert_resources(progress)
        leading = int(work[offset])
        if leading:
            work[offset : offset + divisor_degree + 1] = (
                work[offset : offset + divisor_degree + 1] - leading * divisor
            ) % p
    return poly_trim(work[dividend_degree - divisor_degree + 1 :])


def poly_gcd(a: np.ndarray, b: np.ndarray, p: int, progress: Progress) -> np.ndarray:
    """Explicit monic Euclidean gcd in F_p[z]."""
    a, b = poly_monic(a, p), poly_monic(b, p)
    while not poly_is_zero(b):
        assert_resources(progress)
        a, b = b, poly_monic(poly_remainder(a, b, p, progress), p)
    return a


def derivative_desc(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    derivative_ascending = (np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p
    return np.ascontiguousarray(derivative_ascending[::-1])


def coefficient_digest(ascending_coefficients: list[int]) -> str:
    return hashlib.sha256(",".join(str(int(c)) for c in ascending_coefficients).encode()).hexdigest()


# ---------------------------------------------------------- exact rational integral model

def integralize(matrix: list[list[Fraction]]) -> tuple[int, list[list[int]]]:
    """Use the actual LCM of every Fraction denominator; never infer D by shortcut."""
    scale = 1
    for row in matrix:
        for entry in row:
            scale = scale * entry.denominator // math_gcd(scale, entry.denominator)
    output: list[list[int]] = []
    for row in matrix:
        integer_row: list[int] = []
        for entry in row:
            scaled = scale * entry
            assert scaled.denominator == 1, "D failed to clear a canonical R entry"
            integer_row.append(scaled.numerator)
        output.append(integer_row)
    return scale, output


def factorize(n: int) -> dict[str, int]:
    result: dict[str, int] = {}
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            result[str(divisor)] = result.get(str(divisor), 0) + 1
            n //= divisor
        divisor += 1
    if n > 1:
        result[str(n)] = result.get(str(n), 0) + 1
    return result


def t_key(t: Fraction) -> str:
    return f"{t.numerator}_{t.denominator}"


def case_threshold(n_sites: int) -> dict[str, int]:
    dimension = 1 << n_sites
    pairs = dimension * (dimension - 1) // 2
    gaussian_maximum = 3**n_sites - 2**n_sites
    return {
        "n_sites": n_sites,
        "dim": dimension,
        "pair_count": pairs,
        "gaussian_max_distinct_pair_products": gaussian_maximum,
        "gaussian_gcd_floor": pairs - gaussian_maximum,
    }


def generic_2x3_bound(generic_gcd_upper_bound: int) -> dict[str, object]:
    """Effective degree bookkeeping for the Q(t) specialization theorem.

    For the 2x3 canonical family, e=(b-epsilon)/2 lies in [-4,3].
    M(t)=(2t)^7 q(t)^4 R(t)=(2t)^3(1+t^2)^4 R(t) is therefore an
    integer-polynomial matrix.  Its entries have t-degree at most 26; an
    exterior-square entry has degree at most 52.  The remaining values are
    deliberately coarse upper bounds, sufficient to bound the roots of the
    generic resultant without constructing that enormous resultant.
    """
    threshold = case_threshold(6)
    pair_degree = int(threshold["pair_count"])
    entry_degree = 26
    wedge_entry_degree = 2 * entry_degree
    charpoly_coefficient_degree_bound = wedge_entry_degree * pair_degree
    # If A/G has z-degree n0 and (partial_z A)/G has z-degree n0-1, their
    # Sylvester determinant has t-degree at most
    # 2*wedge_entry_degree*n0*(n0-1), hence the displayed n0<=2016 bound.
    resultant_degree_bound = 2 * wedge_entry_degree * pair_degree * (pair_degree - 1)
    return {
        "status": "THEOREM",
        "domain": "all positive real t outside a finite exceptional set E",
        "polynomial_scalar": "(2t)^7 q(t)^4 = (2t)^3(1+t^2)^4",
        "matrix_entry_t_degree_upper_bound": entry_degree,
        "wedge2_entry_t_degree_upper_bound": wedge_entry_degree,
        "A_coefficient_degree_bound": (
            "for A(t,z)=charpoly(wedge^2 M(t))=z^N+a_1(t)z^(N-1)+... "
            f"with N={pair_degree}, deg_t a_k <= {wedge_entry_degree} k"
        ),
        "A_max_coefficient_t_degree_upper_bound": charpoly_coefficient_degree_bound,
        "generic_gcd_degree_upper_bound": generic_gcd_upper_bound,
        "generic_distinct_pair_products_lower_bound": pair_degree - generic_gcd_upper_bound,
        "gaussian_max_distinct_pair_products": threshold["gaussian_max_distinct_pair_products"],
        "resultant_exceptional_root_bound": resultant_degree_bound,
        "resultant_bound_derivation": (
            "Writing G=gcd_Q(t)[z](A,partial_z A), put A0=A/G and B0=(partial_z A)/G. "
            "Monicity puts all three factors in Q[t][z].  The degree-at-infinity "
            "bound deg_t(coefficient index k)<=52k passes to monic factors and "
            "monic long division gives the same bound for B0.  If n0=deg_z A0, "
            "the Sylvester determinant Res_z(A0,B0) has degree_t at most "
            "2*52*n0*(n0-1) <= 2*52*2016*2015."
        ),
        "specialization_argument": (
            "At t=1/3, the actual-LCM integral model has modular gcd degree 385; "
            "monic integral reduction gives deg gcd_Q <=385.  Every specialization "
            "of the monic generic gcd remains a common monic divisor, so "
            "deg_z G<=385.  For rho=Res_z(A0,B0) nonzero, rho(t)!=0 forces the "
            "specialized gcd to equal G(t,z)."
        ),
        "scope": (
            "At every t>0 with rho(t)!=0, the canonical R(t) differs from the "
            "polynomial M(t) by a nonzero scalar, so it has at least "
            f"{pair_degree - generic_gcd_upper_bound} distinct pair products and "
            "is not a full six-mode Gaussian subset-product spectrum.  No value, "
            "location, or minimal polynomial of an exceptional parameter is claimed."
        ),
    }

# ----------------------------------------------------------------------- reconstruction

def run_pipeline(matrix_rows: list[list[int]], p: int, progress: Progress) -> dict[str, object]:
    """The wave-11 trace/Newton/recurrence/pair-Newton/gcd pipeline at one prime."""
    started = time.perf_counter()
    dim = len(matrix_rows)
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > pairs and p % 2 == 1
    assert_int64_guard(dim, p, "matrix multiplication")
    assert_int64_guard(pairs, p, "Newton dot product")
    matrix = np.array([[entry % p for entry in row] for row in matrix_rows], dtype=np.int64)
    inverses = [0] + [pow(k, -1, p) for k in range(1, pairs + 1)]

    progress.stage = "direct power traces"
    traces = sequential_power_traces(matrix, dim, p, progress)
    after_traces = time.perf_counter()

    progress.stage = "characteristic Newton reconstruction"
    characteristic = newton_from_power_sums(np.array(traces, dtype=np.int64), dim, p, inverses, progress)
    horner_ok = horner_matrix_zero(characteristic, matrix, p)
    after_characteristic = time.perf_counter()

    progress.stage = "Cayley-Hamilton trace recurrence"
    traces = extend_power_sums(traces, characteristic, 2 * pairs, p, progress)
    spot_indices = tuple(sorted({dim + 1, dim + 2, 2 * dim, pairs, 2 * pairs}))
    spots = []
    for exponent in spot_indices:
        direct = trace_power_binary(matrix, exponent, p)
        spots.append({"k": exponent, "direct": direct, "recurrence": int(traces[exponent]), "match": direct == int(traces[exponent])})
    after_recurrence = time.perf_counter()

    progress.stage = "pair-product power sums and Newton reconstruction"
    inverse_two = pow(2, -1, p)
    pair_sums = np.zeros(pairs + 1, dtype=np.int64)
    for exponent in range(1, pairs + 1):
        pair_sums[exponent] = (
            (int(traces[exponent]) * int(traces[exponent]) - int(traces[2 * exponent])) % p
        ) * inverse_two % p
    pair_poly = newton_from_power_sums(pair_sums, pairs, p, inverses, progress)
    after_pair_poly = time.perf_counter()

    progress.stage = "Euclidean gcd(C_2,C_2')"
    derivative = derivative_desc(pair_poly, p)
    common = poly_gcd(pair_poly, derivative, p, progress)
    after_gcd = time.perf_counter()

    divides = poly_is_zero(poly_remainder(pair_poly, common, p, progress)) and poly_is_zero(
        poly_remainder(derivative, common, p, progress)
    )
    characteristic_ascending = [int(value) for value in characteristic[::-1] % p]
    pair_ascending = [int(value) for value in pair_poly[::-1] % p]
    gcd_ascending = [int(value) for value in common[::-1] % p]
    return {
        "prime": p,
        "charpoly_degree": dim,
        "pair_poly_degree": pairs,
        "charpoly_monic": int(characteristic[0]) == 1,
        "pair_poly_monic": int(pair_poly[0]) == 1,
        "charpoly_sha256": coefficient_digest(characteristic_ascending),
        "pair_poly_sha256": coefficient_digest(pair_ascending),
        "gcd_sha256": coefficient_digest(gcd_ascending),
        "gcd_degree": int(common.size) - 1,
        "gcd_monic": int(common[0]) == 1,
        "gcd_divides_both": bool(divides),
        "horner_charpoly_zero_matrix": bool(horner_ok),
        "spot_trace_checks": spots,
        "spot_trace_checks_all_match": all(spot["match"] for spot in spots),
        "timing_seconds": {
            "direct_traces": round(after_traces - started, 3),
            "characteristic_newton": round(after_characteristic - after_traces, 3),
            "trace_recurrence_and_spots": round(after_recurrence - after_characteristic, 3),
            "pair_newton": round(after_pair_poly - after_recurrence, 3),
            "euclidean_gcd": round(after_gcd - after_pair_poly, 3),
            "total": round(after_gcd - started, 3),
        },
        "ru_maxrss_bytes": max_rss_bytes(),
    }


def build_and_run(
    *,
    name: str,
    geometry: str,
    n_sites: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    progress: Progress,
) -> dict[str, object]:
    """Build canonical R, exact-integralize, and run both-prime reconstruction."""
    threshold = case_threshold(n_sites)
    progress.stage = f"canonical build_R for {name}"
    before_build = time.perf_counter()
    rational_matrix, q, epsilon = build_R(n_sites, bonds, t)
    progress.stage = f"LCM integralization for {name}"
    scale, integer_matrix = integralize(rational_matrix)
    build_elapsed = time.perf_counter() - before_build
    del rational_matrix
    per_prime: dict[str, object] = {}
    for p in PRIMES:
        progress.stage = f"{name}, p={p}"
        print(f"pipeline {name} at p={p} ...", flush=True)
        result = run_pipeline(integer_matrix, p, progress)
        per_prime[str(p)] = result
        print(f"  gcd degree {result['gcd_degree']} in {result['timing_seconds']['total']} s", flush=True)
    gcd_degrees = [int(per_prime[str(p)]["gcd_degree"]) for p in PRIMES]
    modular_upper_bound = min(gcd_degrees)
    distinct_lower_bound = int(threshold["pair_count"]) - modular_upper_bound
    record: dict[str, object] = {
        "name": name,
        "geometry": geometry,
        "t": str(t),
        "exp_2K": str(q),
        "bond_count": len(bonds),
        "bond_parity": epsilon,
        "D": str(scale),
        "D_factorization": factorize(scale),
        "build_and_integralize_seconds": round(build_elapsed, 3),
        **threshold,
        "per_prime": per_prime,
        "modular_gcd_degrees": gcd_degrees,
        "modular_gcd_degree_upper_bound": modular_upper_bound,
        "distinct_pair_products_lower_bound": distinct_lower_bound,
        "decisive_non_gaussian": bool(modular_upper_bound < int(threshold["gaussian_gcd_floor"])),
        "ru_maxrss_bytes_after_case": max_rss_bytes(),
    }
    return record


def test_newton_signs() -> bool:
    p = PRIMES[0]
    progress = Progress(stage="Newton sign unit check", started=time.perf_counter())
    # Power sums of roots 2,3,5 are 10,38,160.
    coeffs = newton_from_power_sums(np.array([0, 10, 38, 160], dtype=np.int64), 3, p, [0, 1, pow(2, -1, p), pow(3, -1, p)], progress)
    return [int(x) for x in coeffs] == [1, -10 % p, 31 % p, -30 % p]


def attempt_three_by_three() -> dict[str, object]:
    """Bounded exact attempt; any wall is an honest non-result, not a no-go claim."""
    threshold = case_threshold(9)
    started = time.perf_counter()
    progress = Progress(stage="canonical build_R", started=started, wall_seconds=THREE_BY_THREE_WALL_SECONDS)
    before_rss = max_rss_bytes()
    partial: dict[str, object] = {}
    try:
        with alarm_wall(THREE_BY_THREE_WALL_SECONDS):
            record = build_and_run(
                name="layer_3x3_t_1_3",
                geometry="open 3x3 layer",
                n_sites=9,
                bonds=list(layer_bonds((3, 3), (False, False))),
                t=Fraction(1, 3),
                progress=progress,
            )
        return {
            "status": "COMPUTED",
            "wall_seconds": THREE_BY_THREE_WALL_SECONDS,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "ru_maxrss_bytes_before": before_rss,
            "ru_maxrss_bytes_after": max_rss_bytes(),
            **record,
        }
    except ResourceWall as error:
        return {
            "status": "UNRESOLVED_RESOURCE",
            "reason": str(error),
            "stage": progress.stage,
            "wall_seconds": THREE_BY_THREE_WALL_SECONDS,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
            "ru_maxrss_bytes_before": before_rss,
            "ru_maxrss_bytes_after": max_rss_bytes(),
            "n_sites": 9,
            "dim": threshold["dim"],
            "pair_count": threshold["pair_count"],
            "gaussian_max_distinct_pair_products": threshold["gaussian_max_distinct_pair_products"],
            "gaussian_gcd_floor": threshold["gaussian_gcd_floor"],
            "partial": partial,
        }


# ------------------------------------------------------------------------------- driver

def main() -> int:
    started = time.perf_counter()
    progress = Progress(stage="main", started=started, wall_seconds=CASE_WALL_SECONDS)

    check("both declared modular primes are prime", all(is_prime(p) for p in PRIMES))
    check("2x3 threshold is C(64,2)=2016, 3^6-2^6=665, floor=1351", case_threshold(6) == {"n_sites": 6, "dim": 64, "pair_count": 2016, "gaussian_max_distinct_pair_products": 665, "gaussian_gcd_floor": 1351})
    check("2x4 threshold is C(256,2)=32640, 3^8-2^8=6305, floor=26335", case_threshold(8) == {"n_sites": 8, "dim": 256, "pair_count": 32640, "gaussian_max_distinct_pair_products": 6305, "gaussian_gcd_floor": 26335})
    check("3x3 planned threshold is C(512,2)=130816, 3^9-2^9=19171, floor=111645", case_threshold(9) == {"n_sites": 9, "dim": 512, "pair_count": 130816, "gaussian_max_distinct_pair_products": 19171, "gaussian_gcd_floor": 111645})
    check("Newton sign convention recovers (z-2)(z-3)(z-5)", test_newton_signs())

    guards: dict[str, object] = {}
    for n_sites in (6, 8, 9):
        threshold = case_threshold(n_sites)
        guards[str(n_sites)] = {
            str(p): {
                "dim_times_pm1_sq": threshold["dim"] * (p - 1) ** 2,
                "pairs_times_pm1_sq": threshold["pair_count"] * (p - 1) ** 2,
                "int64_max": INT64_MAX,
                "safe": threshold["dim"] * (p - 1) ** 2 < INT64_MAX and threshold["pair_count"] * (p - 1) ** 2 < INT64_MAX,
            }
            for p in PRIMES
        }
    check("all matrix-product and Newton int64 overflow guards hold", all(value["safe"] for per_prime in guards.values() for value in per_prime.values()))

    cases: dict[str, dict[str, object]] = {}
    bonds23 = list(layer_bonds((2, 3), (False, False)))
    bonds_n6 = [(site, site + 1) for site in range(5)]
    for t in COUPLINGS_2X3:
        layer_name = f"layer_2x3_t_{t_key(t)}"
        chain_name = f"chain_n6_t_{t_key(t)}"
        cases[layer_name] = build_and_run(
            name=layer_name,
            geometry="open 2x3 layer",
            n_sites=6,
            bonds=bonds23,
            t=t,
            progress=progress,
        )
        cases[chain_name] = build_and_run(
            name=chain_name,
            geometry="open chain n=6",
            n_sites=6,
            bonds=bonds_n6,
            t=t,
            progress=progress,
        )

    # The required larger exact layer and its same-site-count free-fermion control.
    cases["layer_2x4_t_1_3"] = build_and_run(
        name="layer_2x4_t_1_3",
        geometry="open 2x4 layer",
        n_sites=8,
        bonds=list(layer_bonds((2, 4), (False, False))),
        t=Fraction(1, 3),
        progress=progress,
    )
    cases["chain_n8_t_1_3"] = build_and_run(
        name="chain_n8_t_1_3",
        geometry="open chain n=8",
        n_sites=8,
        bonds=[(site, site + 1) for site in range(7)],
        t=Fraction(1, 3),
        progress=progress,
    )

    for key, record in cases.items():
        floor = int(record["gaussian_gcd_floor"])
        gcds = [int(value) for value in record["modular_gcd_degrees"]]
        for p in PRIMES:
            pipeline = record["per_prime"][str(p)]
            check(
                f"{key}@{p}: monic degrees, Cayley-Hamilton, recurrence spots, and gcd division",
                bool(pipeline["charpoly_monic"])
                and bool(pipeline["pair_poly_monic"])
                and int(pipeline["charpoly_degree"]) == int(record["dim"])
                and int(pipeline["pair_poly_degree"]) == int(record["pair_count"])
                and bool(pipeline["horner_charpoly_zero_matrix"])
                and bool(pipeline["spot_trace_checks_all_match"])
                and bool(pipeline["gcd_monic"])
                and bool(pipeline["gcd_divides_both"]),
            )
        if key.startswith("chain_"):
            check(f"{key}: corresponding open-chain control equals Gaussian gcd floor {floor}", all(value == floor for value in gcds), f"gcds={gcds}")
        else:
            check(f"{key}: layer modular upper bound is below Gaussian floor {floor}", min(gcds) < floor, f"gcds={gcds}")

    expected_names = {
        "layer_2x3_t_1_5",
        "layer_2x3_t_1_3",
        "layer_2x3_t_2_5",
        "layer_2x3_t_1_2",
    }
    check("all repository-used rational 2x3 couplings were covered", expected_names.issubset(cases))
    grid_names = {f"layer_2x3_t_{t_key(t)}" for t in GRID_COUPLINGS}
    check("all predeclared additional k/20 grid cases were covered", grid_names.issubset(cases))
    check("wave-11 2x3 t=1/3 gcd regression is 385 at both primes", cases["layer_2x3_t_1_3"]["modular_gcd_degrees"] == [385, 385], f"got {cases['layer_2x3_t_1_3']['modular_gcd_degrees']}")
    generic_gcd_upper_bound = int(cases["layer_2x3_t_1_3"]["modular_gcd_degree_upper_bound"])
    generic_theorem = generic_2x3_bound(generic_gcd_upper_bound)
    grid_gcds = {
        f"layer_2x3_t_{t_key(t)}": cases[f"layer_2x3_t_{t_key(t)}"]["modular_gcd_degrees"]
        for t in COUPLINGS_2X3
    }
    check(
        "specialization at t=1/3 gives the generic gcd upper bound 385",
        generic_gcd_upper_bound == 385,
        f"got {generic_gcd_upper_bound}",
    )
    check(
        "all predeclared rational layer points reproduce modular gcd degree 385 at both primes",
        all(gcds == [385, 385] for gcds in grid_gcds.values()),
        f"gcds={grid_gcds}",
    )
    check(
        "generic theorem's distinct-pair lower bound 1631 exceeds the Gaussian maximum 665",
        generic_theorem["generic_distinct_pair_products_lower_bound"]
        > generic_theorem["gaussian_max_distinct_pair_products"],
    )

    three_by_three = attempt_three_by_three()
    check(
        "3x3 route either completed or hit its explicit resource wall",
        three_by_three["status"] in {"COMPUTED", "UNRESOLVED_RESOURCE"},
        f"status={three_by_three['status']}",
    )

    elapsed = time.perf_counter() - started
    artifact = {
        "provenance": {
            "script": "experiments/e118_pair_product_scale.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "precision": "exact Fraction construction and exact integer integralization; exact finite-field arithmetic only after R_int=D*R. No floating point participates in any certificate.",
            "construction": "canonical e38_gaussianity_certificate.build_R; D is the actual LCM of every resulting Fraction denominator, not a guessed scale.",
        },
        "data": {
            "status_tags": {
                "theorem": "[THEOREM] the Q(t) specialization/resultant argument proves a full six-mode no-go for every positive t outside a finite set of at most 422472960 parameters.",
                "computation": "[COMPUTATION] modular values are exact values in F_p; by the monic-integral lemma each is an upper bound on the characteristic-zero gcd degree at that named rational point.",
                "unresolved": "[UNRESOLVED] neither the generic gcd polynomial nor the exceptional resultant is expanded here; the exact exceptional parameters and the 3x3 case remain undecided.",
            },
            "coupling_plan": {
                "repository_rational_points": [str(t) for t in REPOSITORY_COUPLINGS],
                "predeclared_grid": {"formula": "t=k/20 for k=4,5,6,7,8,9,10, reduced as Fractions", "points": [str(t) for t in GRID_COUPLINGS]},
                "computed_2x3_union": [str(t) for t in COUPLINGS_2X3],
            },
            "threshold_derivation": "For a full n-mode subset-product spectrum, distinct eigenvalue slots yield pair-product exponents in {0,1,2}^n with at least one coordinate 1. Thus at most 3^n-2^n pair values occur; C(2^n,2)-(3^n-2^n) is the necessary characteristic-zero gcd floor.",
            "integral_model": "R_int=D R is integer for the actual LCM D of all canonical Fraction denominators. C_2=charpoly(wedge^2 R_int) is monic in Z[z]. Therefore its monic characteristic-zero gcd with C_2' remains monic integral and reduces at every prime, so degree gcd_Fp(C_2,C_2') >= degree gcd_Q(C_2,C_2'). Nonzero scaling by D preserves pair-product equality multiplicities.",
            "pipeline": [
                "s_k=tr(R_int^k), k=1,...,dim, by modular matrix products",
                "Newton identities reconstruct monic charpoly(R_int); matrix Horner checks Cayley-Hamilton",
                "the charpoly recurrence extends s_k through 2*C(dim,2), audited by binary matrix powers",
                "q_m=(s_m^2-s_{2m})/2 reconstructs pair-product traces",
                "Newton identities reconstruct monic C_2=charpoly(wedge^2 R_int)",
                "explicit finite-field Euclidean gcd(C_2,C_2') returns its degree and is checked to divide both polynomials",
            ],
            "primes": {
                "values": list(PRIMES),
                "primality": "trial division through integer sqrt",
                "constraints": "p is odd and exceeds every pair-polynomial degree; all Newton divisors 1,...,degree and 2 are invertible.",
                "int64_guards": guards,
            },
            "cases": cases,
            "generic_coupling_theorem_2x3": generic_theorem,
            "grid_cross_check": {
                "status": "COMPUTATION",
                "per_layer_case_modular_gcd_degrees": grid_gcds,
                "interpretation": "Each rational point is an exact finite certificate. Their uniform values are a deterministic implementation cross-check; the all-but-finite theorem comes instead from monic specialization and the generic resultant.",
            },
            "three_by_three": three_by_three,
            "symbolic_t_remainder": {
                "status": "UNRESOLVED",
                "statement": "The theorem gives no explicit coefficients, locations, or minimal polynomials for the roots of rho(t), and it does not decide any exceptional positive coupling. A direct Z[t][z] subresultant computation could refine that finite set but was not needed for the all-but-finite statement.",
            },
            "resource_policy": {
                "rss_limit_bytes": RSS_LIMIT_BYTES,
                "ordinary_case_wall_seconds": CASE_WALL_SECONDS,
                "three_by_three_wall_seconds": THREE_BY_THREE_WALL_SECONDS,
                "wall_seconds_total": round(elapsed, 3),
                "ru_maxrss_bytes": max_rss_bytes(),
            },
            "inference": "For each layer case, min_p deg gcd_Fp is an upper bound on deg gcd_Q. Hence the characteristic-zero number of distinct pair products is at least pair_count-min_p(deg gcd_Fp). If that lower bound exceeds 3^n-2^n, the named finite layer and rational coupling cannot have a full n-mode Gaussian subset-product spectrum. Two primes are a cross-check only; their agreement is not an equality claim over Q.",
        },
        "checks": CHECKS,
    }
    os.makedirs(ROOT / "results" / "spectral", exist_ok=True)
    with (ROOT / "results" / "spectral" / "pair_product_scale.json").open("w") as handle:
        json.dump(artifact, handle, indent=1)

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(f"PASS: scaled pair-product certificate wrote exact modular evidence in {elapsed:.2f} s.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
