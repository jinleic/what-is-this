"""Ordering-free pair-product certificate for the open 3x3 layer at t=1/3.

Wave 12 (e118) proved the pair-product obstruction for the open 2x3 layer
generically in t and for the open 2x4 layer at t=1/3, but its bounded 3x3
attempt hit the predeclared 90-second wall inside the direct modular
power-trace stage (dim 512, int64 matrix products, ~69 s per prime on this
hardware class).

This experiment removes each measured bottleneck while keeping every step
exact, then lands the 3x3 certificate:

  * exact integral construction without Fraction triple loops:
    P~ = den(t)^n * P_t and d~ = diag(q^e * num(q)^E * den(q)^F) are integer
    matrices, and W = P~ diag(d~) P~ is computed exactly in int64 after a
    2^26 split of the left factor.  R = W / N0 with
    N0 = den(t)^{2n} * num(q)^E * den(q)^F, and the ACTUAL entry-denominator
    LCM is D = N0 / G with G = gcd(gcd_all W_ij, N0), so M = D R = W / G is
    the wave-12 monic integral model (Lemma 2 of proofs/pair_product_scale.md).
  * power traces s_k = tr(M^k), k = 1..512, by float64 BLAS matrix products.
    Every partial sum is an integer below 2^53 (asserted numerically), so the
    double arithmetic is exact; int64 binary powers audit the whole chain.
  * Newton identities give the degree-512 characteristic polynomial; its
    Cayley-Hamilton trace recurrence extends s_k to k = 2*C(512,2) = 261632.
  * q_m = (s_m^2 - s_{2m})/2 and Newton identities reconstruct the monic
    degree-130816 pair polynomial C_2 = charpoly(wedge^2 M) in F_p[z].
  * Euclidean gcd(C_2, C_2') with lazy (drifted) int64 reduction: each
    eliminated coefficient subtracts lead*divisor without reducing the window
    mod p; the per-entry drift is bounded by p + N*(p-1)^2 < 2^63 (asserted).

Decisive quantity: a full nine-mode Gaussian subset-product spectrum forces
deg gcd >= C(512,2) - (3^9 - 2^9) = 111645.  The monic integral lemma turns
each modular gcd degree into an upper bound on the characteristic-zero value,
so a modular degree below 111645 at t=1/3 certifies more than 19171 distinct
pair products and kills the Gaussian form for the first layer with a
degree-4 vertex.

Controls: the open 9-site chain must return exactly the Gaussian floor
111645 at both primes, and the wave-12 2x3/2x4 layer values (385 and 9329,
with coefficient digests) must replay through this new pipeline.

Run: PYTHONPATH=src .venv/bin/python experiments/e122_pair_product_3x3.py
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
from dataclasses import dataclass, field
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd as math_gcd
from math import isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from ising.transfer_matrix import layer_bonds  # noqa: E402
from e38_gaussianity_certificate import build_R  # noqa: E402 -- canonical rational construction

INT64_MAX = 2**63 - 1
FLOAT_EXACT_MAX = 2**53  # integers below this are exactly representable in float64
PRIMES = (1_000_003, 2_000_003)
RSS_LIMIT_BYTES = 7_500_000_000
SPLIT = 1 << 26

# Predeclared per-stage walls (seconds), fixed before any case is computed.
STAGE_WALL_SECONDS = {
    "build": 300,
    "traces": 240,
    "charpoly": 120,
    "recurrence": 240,
    "pair_newton": 600,
    "gcd": 1200,
}
CASE_WALL_SECONDS = 2400

WAVE12_ARTIFACT = ROOT / "results" / "spectral" / "pair_product_scale.json"
WAVE12_LAYER_D = {
    "layer_2x3_t_1_3": 8_968_066_875,
    "chain_n6_t_1_3": 597_871_125,
    "layer_2x4_t_1_3": 32_688_603_759_375,
    "chain_n8_t_1_3": 726_413_416_875,
}

CHECKS: list[dict[str, object]] = []
FAILS: list[str] = []


class ResourceWall(RuntimeError):
    """A predeclared time or RSS wall stopped a stage."""


@dataclass
class Progress:
    stage: str = "initialization"
    started: float = field(default_factory=time.perf_counter)
    wall_seconds: int = CASE_WALL_SECONDS


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    CHECKS.append({"name": name, "ok": bool(ok), "detail": detail})
    if not ok:
        FAILS.append(name)
    return bool(ok)


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


def assert_resources(progress: Progress) -> None:
    elapsed = time.perf_counter() - progress.started
    assert elapsed < progress.wall_seconds, (
        f"case time wall {progress.wall_seconds}s at {progress.stage} after {elapsed:.3f}s"
    )
    rss = max_rss_bytes()
    assert rss < RSS_LIMIT_BYTES, f"RSS wall {RSS_LIMIT_BYTES} bytes at {progress.stage} (ru_maxrss {rss})"


@contextmanager
def alarm_wall(seconds: int, stage: str):
    """Interrupt Python-level exact work at a real-time wall on POSIX."""

    def handler(signum, frame):
        raise ResourceWall(f"stage wall {seconds}s fired at {stage}")

    previous = signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


# ------------------------------------------------------------------ arithmetic guards

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for divisor in range(2, isqrt(n) + 1):
        if n % divisor == 0:
            return False
    return True


def assert_int64_guard(length: int, p: int, kind: str) -> None:
    assert length * (p - 1) ** 2 < INT64_MAX, f"int64 {kind} guard violated"


def assert_float_guard(length: int, p: int, kind: str) -> None:
    """Every partial sum of `length` products of [0,p) integers stays below 2^53."""
    assert length * (p - 1) ** 2 < FLOAT_EXACT_MAX, f"float64 exactness {kind} guard violated"


# ------------------------------------------------- exact integral model (fast route)

@dataclass
class IntegralModel:
    n_sites: int
    bonds: list[tuple[int, int]]
    t: Fraction
    num: int
    den: int
    eps: int
    e_values: list[int]
    E: int
    F: int
    N0: int
    G: int
    D: int
    W_hi: np.ndarray
    W_lo: np.ndarray
    max_hamming_sum: int


def geometry_exponents(n: int, bonds, t: Fraction):
    """q = exp(2K), bond parity eps, and exponents e = (b - eps)/2 of the canonical model."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    b: list[int] = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = b[0] % 2
    assert all((x - eps) % 2 == 0 for x in b), "bond-count parity assumption violated"
    return q, eps, [(x - eps) // 2 for x in b]


def fast_integral_model(n: int, bonds, t: Fraction, progress: Progress) -> IntegralModel:
    """Exact integer W with R = W/N0 and the actual LCM D = N0/G; no Fraction matmul.

    P~ = den(t)^n * P_t has entries num(t)^h * den(t)^(n-h); the scaled diagonal
    d~_k = q^{e_k} * num(q)^E * den(q)^F = num^{e+E} * den^{F-e} is integral for
    e in [e_min, e_max] with E = max(0, -e_min), F = max(0, e_max).  Then

        W = P~ diag(d~) P~,   N0 = den(t)^{2n} * num(q)^E * den(q)^F,

    R = W/N0 exactly, and the LCM of the denominators of R equals N0/G with
    G = gcd(gcd_{ij} W_ij, N0), so M = D*R = W/G is the monic integral model.
    """
    progress.stage = "build"
    q, eps, e_values = geometry_exponents(n, bonds, t)
    num, den = q.numerator, q.denominator
    a, bden = t.numerator, t.denominator
    E, F = max(0, -min(e_values)), max(0, max(e_values))
    dim = 1 << n

    # P~ entries num(t)^h * den(t)^(n-h); max_hamming_sum bounds h(i,k)+h(j,k).
    max_h_sum = 0
    powers = np.empty((dim, dim), dtype=np.int64)
    for k in range(dim):
        for l in range(dim):
            h = bin(k ^ l).count("1")
            if 2 * h > max_h_sum:
                max_h_sum = 2 * h
            powers[k, l] = a**h * bden ** (n - h)

    diagonal = np.array([num ** (x + E) * den ** (F - x) for x in e_values], dtype=np.int64)
    left = powers * diagonal[None, :]
    # int64 split matmul: left = left_hi * SPLIT + left_lo.  Both products and
    # both 512-term partial sums stay far below 2^63 (asserted), so the plain
    # int64 matmuls below are exact integer arithmetic.
    left_hi = left >> 26
    left_lo = left & (SPLIT - 1)
    max_power = int(powers.max())
    assert int(np.abs(left_hi).max()) * max_power * dim < INT64_MAX, "split hi int64 guard violated"
    assert (SPLIT - 1) * max_power * dim < INT64_MAX, "split lo int64 guard violated"
    W_hi = left_hi @ powers
    W_lo = left_lo @ powers
    assert int(np.abs(W_hi).max()) < INT64_MAX and int(np.abs(W_lo).max()) < INT64_MAX

    N0 = bden ** (2 * n) * num**E * den**F
    content = 0
    for hi, lo in zip(W_hi.ravel(), W_lo.ravel()):
        content = math_gcd(content, (int(hi) << 26) + int(lo))
        if content == 1:
            break
    G = math_gcd(content, N0)
    return IntegralModel(
        n_sites=n,
        bonds=list(bonds),
        t=t,
        num=num,
        den=den,
        eps=eps,
        e_values=e_values,
        E=E,
        F=F,
        N0=N0,
        G=G,
        D=N0 // G,
        W_hi=W_hi,
        W_lo=W_lo,
        max_hamming_sum=max_h_sum,
    )


def W_mod_p(model: IntegralModel, p: int) -> np.ndarray:
    """W mod p from the split halves; the combined sum stays below 2^63."""
    combined = (model.W_hi % p) * (SPLIT % p) + model.W_lo
    assert int(np.abs(combined).max()) < INT64_MAX
    return combined % p


def matrix_mod_p(model: IntegralModel, p: int) -> np.ndarray:
    """M = W/G mod p; G divides N0, whose prime factors exclude p."""
    G = model.G % p
    assert G != 0
    return (W_mod_p(model, p) * pow(int(G), -1, p)) % p


def spot_fraction_check(model: IntegralModel, entries: list[tuple[int, int]]) -> list[bool]:
    """Verify R_ij = W_ij/N0 on named entries by direct Fraction row sums."""
    n, t, q = model.n_sites, model.t, Fraction(model.num, model.den)
    dim = 1 << n
    bonds = model.bonds
    b = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(spins[i] * spins[j] for i, j in bonds))
    parity = b[0] % 2
    d = [q ** ((x - parity) // 2) for x in b]
    results = []
    for (i, j) in entries:
        acc = Fraction(0)
        for k in range(dim):
            acc += (t ** bin(i ^ k).count("1")) * d[k] * (t ** bin(j ^ k).count("1"))
        w_entry = (int(model.W_hi[i, j]) << 26) + int(model.W_lo[i, j])
        results.append(acc * model.N0 == w_entry)
    return results


# ------------------------------------------------------- modular matrix arithmetic

def float_matmul_mod(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """Exact (a@b) mod p through float64 BLAS; all partial sums are integers < 2^53."""
    assert_float_guard(a.shape[0], p, "float matmul")
    product = a.astype(np.float64) @ b.astype(np.float64)
    return product.astype(np.int64) % p


def int_matmul_mod(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    assert_int64_guard(a.shape[0], p, "int64 matmul")
    return (a @ b) % p


def sequential_power_traces(matrix: np.ndarray, max_power: int, p: int, progress: Progress) -> np.ndarray:
    """[0, tr(M), ..., tr(M^max_power)] by float64-BLAS modular matrix products."""
    dim = matrix.shape[0]
    assert_float_guard(dim, p, "power traces")
    base = matrix.astype(np.float64)
    power = np.eye(dim, dtype=np.float64)
    traces = np.zeros(max_power + 1, dtype=np.int64)
    for exponent in range(1, max_power + 1):
        if exponent & 63 == 0:
            assert_resources(progress)
        power = power @ base
        reduced = power.astype(np.int64) % p
        traces[exponent] = int(reduced.trace()) % p
        power = reduced.astype(np.float64)
    return traces


def trace_power_binary_int64(matrix: np.ndarray, exponent: int, p: int) -> int:
    """Independent int64 binary-power trace; audits the float64 chain exactly."""
    dim = matrix.shape[0]
    assert_int64_guard(dim, p, "binary power")
    result = np.eye(dim, dtype=np.int64)
    base = matrix.copy()
    e = exponent
    while e:
        if e & 1:
            result = (result @ base) % p
        e >>= 1
        if e:
            base = (base @ base) % p
    return int(result.trace()) % p


def horner_matrix_zero(charpoly_desc: np.ndarray, matrix: np.ndarray, p: int) -> bool:
    """Cayley-Hamilton validation chi(M) = 0 in F_p^{dim x dim} via float64 BLAS."""
    dim = matrix.shape[0]
    assert_float_guard(dim, p, "Horner")
    identity_f = np.eye(dim, dtype=np.float64)
    value = np.eye(dim, dtype=np.float64)
    for coefficient in charpoly_desc[1:]:
        value = value @ matrix
        value = value + (int(coefficient) * identity_f)
        value = (value.astype(np.int64) % p).astype(np.float64)
    return bool(np.all(value == 0))


# ---------------------------------------------------------------- Newton machinery

def newton_from_power_sums(
    power_sums: np.ndarray, degree: int, p: int, inverses: list[int], progress: Progress
) -> np.ndarray:
    """DESCENDING [1,a1,...,a_degree] from power sums via Newton identities."""
    assert_int64_guard(degree, p, "Newton dot product")
    coeffs = np.zeros(degree + 1, dtype=np.int64)
    coeffs[0] = 1
    for k in range(1, degree + 1):
        if k & 1023 == 0:
            assert_resources(progress)
        total = int(np.dot(power_sums[1 : k + 1], np.ascontiguousarray(coeffs[k - 1 :: -1]))) % p
        coeffs[k] = (-total * inverses[k]) % p
    return coeffs


def extend_power_sums(
    traces: np.ndarray, charpoly_desc: np.ndarray, upto: int, p: int, progress: Progress
) -> np.ndarray:
    """Cayley-Hamilton trace recurrence s_k = -sum a_i s_{k-i} beyond matrix dimension."""
    dim = charpoly_desc.size - 1
    assert_int64_guard(dim, p, "trace recurrence")
    coefficients = charpoly_desc[1:]
    out = np.zeros(upto + 1, dtype=np.int64)
    out[: traces.size] = traces
    for k in range(dim + 1, upto + 1):
        if k & 4095 == 0:
            assert_resources(progress)
        total = int(np.dot(coefficients, out[k - 1 : k - dim - 1 : -1])) % p
        out[k] = (-total) % p
    return out


# ------------------------------------------------- F_p[z], descending, lazy reduction

def poly_trim(poly: np.ndarray) -> np.ndarray:
    if poly.size == 1:
        return poly
    index = int(np.argmax(poly != 0))
    if poly[index] == 0:
        return np.zeros(1, dtype=np.int64)
    return np.ascontiguousarray(poly[index:])


def poly_is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def poly_monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = poly_trim(np.asarray(poly, dtype=np.int64) % p)
    if poly_is_zero(poly):
        return poly
    return (poly * pow(int(poly[0]), -1, p)) % p


def remainder_lazy(dividend: np.ndarray, divisor: np.ndarray, p: int, progress: Progress) -> np.ndarray:
    """Remainder of division in F_p[z] with lazy int64 reduction.

    Both inputs are descending int64 arrays with values in [0, p); the divisor
    need not be monic.  Each eliminated coefficient subtracts lead*divisor from
    its window WITHOUT reducing the window mod p.  Every intermediate entry
    stays within p + (#touches)*(p-1)^2 <= p + N*(p-1)^2 < 2^63 (asserted), so
    the lazy arithmetic is exact and the tail is reduced once at the end.
    """
    dividend_degree = dividend.size - 1
    divisor_degree = divisor.size - 1
    touches = dividend_degree - divisor_degree + 1
    assert p + max(touches, 0) * (p - 1) ** 2 < INT64_MAX, "lazy-division drift guard violated"
    if dividend_degree < divisor_degree:
        return poly_trim(dividend)
    work = dividend.copy()
    lead_inv = pow(int(divisor[0]), -1, p)
    scratch = np.empty(divisor_degree + 1, dtype=np.int64)
    for offset in range(touches):
        if offset & 1023 == 0:
            assert_resources(progress)
        top = int(work[offset]) % p
        if top:
            lead = top * lead_inv % p
            np.multiply(divisor, lead, out=scratch)
            work[offset : offset + divisor_degree + 1] -= scratch
    return poly_trim(work[dividend_degree - divisor_degree + 1 :] % p)


def poly_gcd(a: np.ndarray, b: np.ndarray, p: int, progress: Progress) -> tuple[np.ndarray, int]:
    """Euclidean gcd in F_p[z] with lazy remainders; the result is monic."""
    a = poly_trim(np.asarray(a, dtype=np.int64) % p)
    b = poly_trim(np.asarray(b, dtype=np.int64) % p)
    steps = 0
    while not poly_is_zero(b):
        assert_resources(progress)
        a, b = b, remainder_lazy(a, b, p, progress)
        steps += 1
    return poly_monic(a, p), steps


def derivative_desc(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    derivative_ascending = (np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p
    return np.ascontiguousarray(derivative_ascending[::-1])


def coefficient_digest(ascending_coefficients) -> str:
    return hashlib.sha256(",".join(str(int(c)) for c in ascending_coefficients).encode()).hexdigest()


# ----------------------------------------------------------------------- thresholds

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


# ------------------------------------------------------------------------ pipeline

def run_pipeline(matrix: np.ndarray, p: int, progress: Progress) -> dict[str, object]:
    """Trace/Newton/recurrence/pair-Newton/gcd pipeline at one prime."""
    started = time.perf_counter()
    dim = matrix.shape[0]
    pairs = dim * (dim - 1) // 2
    assert is_prime(p) and p > 2 * pairs and p % 2 == 1
    assert_int64_guard(dim, p, "int64 binary powers")
    assert_int64_guard(pairs, p, "Newton dot product")
    assert_float_guard(dim, p, "float64 matmul")
    inverses = [0] + [pow(k, -1, p) for k in range(1, pairs + 1)]

    progress.stage = "traces"
    # One full int64 product cross-checks the float64 BLAS path on this matrix.
    float_check = np.array_equal(float_matmul_mod(matrix, matrix, p), int_matmul_mod(matrix, matrix, p))
    with alarm_wall(STAGE_WALL_SECONDS["traces"], "traces"):
        traces = sequential_power_traces(matrix, dim, p, progress)
    after_traces = time.perf_counter()
    rss_after_traces = max_rss_bytes()

    progress.stage = "charpoly"
    with alarm_wall(STAGE_WALL_SECONDS["charpoly"], "charpoly"):
        characteristic = newton_from_power_sums(traces, dim, p, inverses, progress)
        horner_ok = horner_matrix_zero(characteristic, matrix, p)
    after_characteristic = time.perf_counter()
    rss_after_characteristic = max_rss_bytes()

    progress.stage = "recurrence"
    with alarm_wall(STAGE_WALL_SECONDS["recurrence"], "recurrence"):
        full_traces = extend_power_sums(traces, characteristic, 2 * pairs, p, progress)
        spot_indices = tuple(
            sorted({min(x, 2 * pairs) for x in (dim + 1, dim + 2, 2 * dim, 4096, 65536, pairs, 2 * pairs)})
        )
        spots = []
        for exponent in spot_indices:
            direct = trace_power_binary_int64(matrix, exponent, p)
            spots.append(
                {
                    "k": exponent,
                    "direct_int64_binary": direct,
                    "recurrence": int(full_traces[exponent]),
                    "match": direct == int(full_traces[exponent]),
                }
            )
    after_recurrence = time.perf_counter()
    rss_after_recurrence = max_rss_bytes()

    progress.stage = "pair_newton"
    with alarm_wall(STAGE_WALL_SECONDS["pair_newton"], "pair_newton"):
        inverse_two = pow(2, -1, p)
        first = full_traces[1 : pairs + 1]
        second = full_traces[2 : 2 * pairs + 1 : 2]
        assert first.size == pairs and second.size == pairs
        body = ((first * first - second) % p) * inverse_two % p  # int64: (p-1)^2 < 2^42
        pair_sums = np.concatenate(([0], body))  # index 0 is the Newton placeholder
        pair_poly = newton_from_power_sums(pair_sums, pairs, p, inverses, progress)
    after_pair_poly = time.perf_counter()
    rss_after_pair_poly = max_rss_bytes()

    progress.stage = "gcd"
    with alarm_wall(STAGE_WALL_SECONDS["gcd"], "gcd"):
        derivative = derivative_desc(pair_poly, p)
        common, gcd_steps = poly_gcd(pair_poly, derivative, p, progress)
        after_gcd_algorithm = time.perf_counter()
        # Belt-and-braces: the reported monic gcd must divide both inputs exactly.
        divides = poly_is_zero(remainder_lazy(pair_poly, common, p, progress)) and poly_is_zero(
            remainder_lazy(derivative, common, p, progress)
        )
    after_gcd = time.perf_counter()
    rss_after_gcd = max_rss_bytes()
    characteristic_ascending = [int(v) for v in characteristic[::-1] % p]
    pair_ascending = [int(v) for v in pair_poly[::-1] % p]
    gcd_ascending = [int(v) for v in common[::-1] % p]
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
        "gcd_euclid_steps": gcd_steps,
        "float64_vs_int64_first_product_equal": bool(float_check),
        "horner_charpoly_zero_matrix": bool(horner_ok),
        "spot_trace_checks": spots,
        "spot_trace_checks_all_match": all(spot["match"] for spot in spots),
        "timing_seconds": {
            "traces": round(after_traces - started, 3),
            "characteristic_newton_and_horner": round(after_characteristic - after_traces, 3),
            "trace_recurrence_and_spots": round(after_recurrence - after_characteristic, 3),
            "pair_newton": round(after_pair_poly - after_recurrence, 3),
            "euclidean_gcd": round(after_gcd - after_pair_poly, 3),
            "euclidean_gcd_algorithm": round(after_gcd_algorithm - after_pair_poly, 3),
            "gcd_division_checks": round(after_gcd - after_gcd_algorithm, 3),
            "total": round(after_gcd - started, 3),
        },
        "stage_resources": {
            "traces": {
                "predeclared_wall_seconds": STAGE_WALL_SECONDS["traces"],
                "observed_seconds": round(after_traces - started, 3),
                "ru_maxrss_bytes_after_stage": rss_after_traces,
            },
            "charpoly": {
                "predeclared_wall_seconds": STAGE_WALL_SECONDS["charpoly"],
                "observed_seconds": round(after_characteristic - after_traces, 3),
                "ru_maxrss_bytes_after_stage": rss_after_characteristic,
            },
            "recurrence": {
                "predeclared_wall_seconds": STAGE_WALL_SECONDS["recurrence"],
                "observed_seconds": round(after_recurrence - after_characteristic, 3),
                "ru_maxrss_bytes_after_stage": rss_after_recurrence,
            },
            "pair_newton": {
                "predeclared_wall_seconds": STAGE_WALL_SECONDS["pair_newton"],
                "observed_seconds": round(after_pair_poly - after_recurrence, 3),
                "ru_maxrss_bytes_after_stage": rss_after_pair_poly,
            },
            "gcd": {
                "predeclared_wall_seconds": STAGE_WALL_SECONDS["gcd"],
                "observed_seconds": round(after_gcd - after_pair_poly, 3),
                "ru_maxrss_bytes_after_stage": rss_after_gcd,
            },
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
    validation: dict[str, object] | None = None,
) -> dict[str, object]:
    """Fast exact integral model plus both-prime reconstruction."""
    threshold = case_threshold(n_sites)
    progress.stage = "build"
    before_build = time.perf_counter()
    with alarm_wall(STAGE_WALL_SECONDS["build"], "build"):
        model = fast_integral_model(n_sites, bonds, t, progress)
    build_elapsed = time.perf_counter() - before_build
    build_rss = max_rss_bytes()
    per_prime: dict[str, object] = {}
    for p in PRIMES:
        progress.stage = f"{name}, p={p}"
        print(f"pipeline {name} at p={p} ...", flush=True)
        matrix_started = time.perf_counter()
        matrix = matrix_mod_p(model, p)
        matrix_mod_elapsed = time.perf_counter() - matrix_started
        matrix_mod_rss = max_rss_bytes()
        result = run_pipeline(matrix, p, progress)
        result["matrix_mod_p_seconds"] = round(matrix_mod_elapsed, 3)
        result["matrix_mod_p_ru_maxrss_bytes_after_construction"] = matrix_mod_rss
        per_prime[str(p)] = result
        print(
            f"  gcd degree {result['gcd_degree']} (Euclid steps {result['gcd_euclid_steps']})"
            f" in {result['timing_seconds']['total']} s",
            flush=True,
        )
        del matrix
    gcd_degrees = [int(per_prime[str(p)]["gcd_degree"]) for p in PRIMES]
    modular_upper_bound = min(gcd_degrees)
    distinct_lower_bound = int(threshold["pair_count"]) - modular_upper_bound
    record: dict[str, object] = {
        "name": name,
        "geometry": geometry,
        "t": str(t),
        "exp_2K": f"{model.num}/{model.den}",
        "bond_count": len(bonds),
        "bond_parity": model.eps,
        "N0": str(model.N0),
        "N0_factorization": factorize(model.N0),
        "content_gcd_G": str(model.G),
        "D": str(model.D),
        "D_factorization": factorize(model.D),
        "e_min": min(model.e_values),
        "e_max": max(model.e_values),
        "E": model.E,
        "F": model.F,
        "max_hamming_sum_bound": model.max_hamming_sum,
        "build_seconds": round(build_elapsed, 3),
        "build_resource": {
            "predeclared_wall_seconds": STAGE_WALL_SECONDS["build"],
            "observed_seconds": round(build_elapsed, 3),
            "ru_maxrss_bytes_after_stage": build_rss,
        },
        **threshold,
        "per_prime": per_prime,
        "modular_gcd_degrees": gcd_degrees,
        "modular_gcd_degree_upper_bound": modular_upper_bound,
        "distinct_pair_products_lower_bound": distinct_lower_bound,
        "decisive_non_gaussian": bool(modular_upper_bound < int(threshold["gaussian_gcd_floor"])),
        "at_gaussian_floor": bool(modular_upper_bound == int(threshold["gaussian_gcd_floor"])),
        "ru_maxrss_bytes_after_case": max_rss_bytes(),
    }
    if validation is not None:
        record["validation"] = validation
    return record


# ------------------------------------------------------------- generic-t bookkeeping

def generic_3x3_bound(generic_gcd_upper_bound: int, model: IntegralModel) -> dict[str, object]:
    """Effective degree bookkeeping for the Q(t) specialization theorem, 3x3.

    M(t) = (2t)^{2F} q(t)^F R(t) is an integer-polynomial matrix because
    (2t)^{2F} q(t)^{F+e} = (2t)^{F-e} (1+t^2)^{F+e} for e in [e_min, e_max].
    """
    threshold = case_threshold(9)
    pair_degree = int(threshold["pair_count"])
    n = 9
    two_F = 2 * model.F
    entry_degree = 2 * n + two_F + model.F + max(model.e_values)  # h1+h2 + (F-e) + 2(F+e), maximized
    wedge_entry_degree = 2 * entry_degree
    resultant_degree_bound = 2 * wedge_entry_degree * pair_degree * (pair_degree - 1)
    return {
        "status": "THEOREM",
        "domain": "all positive real t outside a finite exceptional set E",
        "polynomial_scalar": f"(2t)^{two_F} q(t)^{model.F}",
        "e_range": [min(model.e_values), max(model.e_values)],
        "matrix_entry_t_degree_upper_bound": entry_degree,
        "wedge2_entry_t_degree_upper_bound": wedge_entry_degree,
        "A_coefficient_degree_bound": (
            "for A(t,z)=charpoly(wedge^2 M(t))=z^N+a_1(t)z^(N-1)+... "
            f"with N={pair_degree}, deg_t a_k <= {wedge_entry_degree} k"
        ),
        "generic_gcd_degree_upper_bound": generic_gcd_upper_bound,
        "generic_distinct_pair_products_lower_bound": pair_degree - generic_gcd_upper_bound,
        "gaussian_max_distinct_pair_products": threshold["gaussian_max_distinct_pair_products"],
        "resultant_exceptional_root_bound": resultant_degree_bound,
        "resultant_bound_derivation": (
            "Writing G=gcd_Q(t)[z](A,partial_z A), put A0=A/G and B0=(partial_z A)/G. "
            "Monicity puts all three factors in Q[t][z]. The degree-at-infinity bound "
            f"deg_t(coefficient index k)<={wedge_entry_degree} k passes to monic factors, and "
            "monic long division gives the same bound for B0. If n0=deg_z A0, "
            "the Sylvester determinant Res_z(A0,B0) has degree_t at most "
            f"2*{wedge_entry_degree}*n0*(n0-1) <= 2*{wedge_entry_degree}*{pair_degree}*{pair_degree - 1}."
        ),
        "specialization_argument": (
            "At t=1/3 the actual-LCM integral model has modular gcd degree "
            f"{generic_gcd_upper_bound}; monic integral reduction gives deg gcd_Q <= that value. "
            "Every specialization of the monic generic gcd remains a common monic divisor, so "
            "deg_z G <= that value. For rho=Res_z(A0,B0) nonzero, rho(t)!=0 forces the "
            "specialized gcd to equal G(t,z)."
        ),
        "scope": (
            "At every t>0 with rho(t)!=0, the canonical R(t) differs from the polynomial M(t) "
            "by a nonzero scalar, so it has at least "
            f"{pair_degree - generic_gcd_upper_bound} distinct pair products and is not a full "
            "nine-mode Gaussian subset-product spectrum. No value, location, or minimal "
            "polynomial of an exceptional parameter is claimed."
        ),
    }


# ------------------------------------------------------------------------------ driver

def validate_against_canonical(name: str, n: int, bonds, t: Fraction, model: IntegralModel) -> dict[str, object]:
    """Full-entry comparison for dim 64; actual-LCM D regression for dim 256."""
    if n <= 6:
        R, _, _ = build_R(n, bonds, t)
        dim = 1 << n
        ok = True
        for i in range(dim):
            for j in range(dim):
                value = R[i][j] * model.N0
                if value.denominator != 1 or value.numerator != (
                    (int(model.W_hi[i, j]) << 26) + int(model.W_lo[i, j])
                ):
                    ok = False
                    break
        return {"mode": "full entrywise build_R comparison", "full_match": ok}
    expected = WAVE12_LAYER_D.get(name)
    return {
        "mode": "actual-LCM D regression against the wave-12 artifact",
        "wave12_D": expected,
        "D_match": expected == model.D,
    }


def main() -> int:
    started = time.perf_counter()
    progress = Progress(stage="main", started=started, wall_seconds=CASE_WALL_SECONDS)

    check("both declared modular primes are prime", all(is_prime(p) for p in PRIMES))
    check(
        "3x3 threshold is C(512,2)=130816, 3^9-2^9=19171, floor=111645",
        case_threshold(9)
        == {
            "n_sites": 9,
            "dim": 512,
            "pair_count": 130816,
            "gaussian_max_distinct_pair_products": 19171,
            "gaussian_gcd_floor": 111645,
        },
    )
    check("case_threshold(6) matches the wave-12 2x3 values", case_threshold(6)["gaussian_gcd_floor"] == 1351)
    check("case_threshold(8) matches the wave-12 2x4 values", case_threshold(8)["gaussian_gcd_floor"] == 26335)

    guards: dict[str, object] = {}
    for n_sites in (6, 8, 9):
        threshold = case_threshold(n_sites)
        guards[str(n_sites)] = {
            str(p): {
                "dim_times_pm1_sq": threshold["dim"] * (p - 1) ** 2,
                "pairs_times_pm1_sq": threshold["pair_count"] * (p - 1) ** 2,
                "lazy_drift_bound": p + threshold["pair_count"] * (p - 1) ** 2,
                "int64_max": INT64_MAX,
                "float_exact_max": FLOAT_EXACT_MAX,
                "safe_int64": threshold["dim"] * (p - 1) ** 2 < INT64_MAX
                and threshold["pair_count"] * (p - 1) ** 2 < INT64_MAX
                and p + threshold["pair_count"] * (p - 1) ** 2 < INT64_MAX,
                "safe_float64_matmul": threshold["dim"] * (p - 1) ** 2 < FLOAT_EXACT_MAX,
            }
            for p in PRIMES
        }
    check(
        "all int64 matmul/Newton/lazy-drift and float64-BLAS exactness guards hold",
        all(v["safe_int64"] and v["safe_float64_matmul"] for per_n in guards.values() for v in per_n.values()),
    )
    check("primes exceed 2*130816 and avoid 3 and 5", all(p > 2 * 130816 and p % 3 and p % 5 for p in PRIMES))

    wave12_cases = json.loads(WAVE12_ARTIFACT.read_text())["data"]["cases"] if WAVE12_ARTIFACT.exists() else {}

    cases: dict[str, dict[str, object]] = {}
    bonds23 = list(layer_bonds((2, 3), (False, False)))
    bonds_n6 = [(site, site + 1) for site in range(5)]
    bonds24 = list(layer_bonds((2, 4), (False, False)))
    bonds_n8 = [(site, site + 1) for site in range(7)]
    bonds33 = list(layer_bonds((3, 3), (False, False)))
    bonds_n9 = [(site, site + 1) for site in range(8)]
    t = Fraction(1, 3)

    for name, geometry, n, bonds in (
        ("layer_2x3_t_1_3", "open 2x3 layer", 6, bonds23),
        ("chain_n6_t_1_3", "open chain n=6", 6, bonds_n6),
        ("layer_2x4_t_1_3", "open 2x4 layer", 8, bonds24),
        ("chain_n8_t_1_3", "open chain n=8", 8, bonds_n8),
    ):
        model = fast_integral_model(n, bonds, t, progress)
        validation = validate_against_canonical(name, n, bonds, t, model)
        del model
        cases[name] = build_and_run(
            name=name, geometry=geometry, n_sites=n, bonds=bonds, t=t, progress=progress, validation=validation
        )

    # The decisive 3x3 layer and its same-site-count free-fermion control.
    model33 = fast_integral_model(9, bonds33, t, progress)
    rng = np.random.default_rng(20260816)
    spot_entries = [(0, j) for j in range(0, 512, 37)] + [
        (int(i), int(j)) for i, j in zip(rng.integers(0, 512, 16), rng.integers(0, 512, 16))
    ]
    spot_ok = spot_fraction_check(model33, spot_entries)
    validation33 = {
        "mode": "first-row stride plus 16 random entries against direct Fraction row sums",
        "checked_entries": len(spot_entries),
        "all_match": all(spot_ok),
        "failed_entries": [e for e, ok in zip(spot_entries, spot_ok) if not ok],
    }
    cases["layer_3x3_t_1_3"] = build_and_run(
        name="layer_3x3_t_1_3",
        geometry="open 3x3 layer (first layer with a degree-4 vertex)",
        n_sites=9,
        bonds=bonds33,
        t=t,
        progress=progress,
        validation=validation33,
    )

    model_n9 = fast_integral_model(9, bonds_n9, t, progress)
    spot_ok_n9 = spot_fraction_check(model_n9, spot_entries)
    validation_n9 = {
        "mode": "first-row stride plus 16 random entries against direct Fraction row sums",
        "checked_entries": len(spot_entries),
        "all_match": all(spot_ok_n9),
        "failed_entries": [e for e, ok in zip(spot_entries, spot_ok_n9) if not ok],
    }
    cases["chain_n9_t_1_3"] = build_and_run(
        name="chain_n9_t_1_3",
        geometry="open chain n=9",
        n_sites=9,
        bonds=bonds_n9,
        t=t,
        progress=progress,
        validation=validation_n9,
    )

    # ------------------------------------------------------------------- checks
    for key, record in cases.items():
        floor = int(record["gaussian_gcd_floor"])
        gcds = [int(v) for v in record["modular_gcd_degrees"]]
        for p in PRIMES:
            pipeline = record["per_prime"][str(p)]
            check(
                f"{key}@{p}: monic degrees, float/int product equality, Cayley-Hamilton, "
                "recurrence spots, and gcd division",
                bool(pipeline["charpoly_monic"])
                and bool(pipeline["pair_poly_monic"])
                and int(pipeline["charpoly_degree"]) == int(record["dim"])
                and int(pipeline["pair_poly_degree"]) == int(record["pair_count"])
                and bool(pipeline["float64_vs_int64_first_product_equal"])
                and bool(pipeline["horner_charpoly_zero_matrix"])
                and bool(pipeline["spot_trace_checks_all_match"])
                and bool(pipeline["gcd_monic"])
                and bool(pipeline["gcd_divides_both"]),
            )
        if key.startswith("chain_"):
            check(
                f"{key}: open-chain control equals the Gaussian gcd floor {floor}",
                all(v == floor for v in gcds),
                f"gcds={gcds}",
            )
        else:
            check(
                f"{key}: layer modular upper bound is below the Gaussian floor {floor}",
                min(gcds) < floor,
                f"gcds={gcds}",
            )

    check("2x3 fast construction matches canonical build_R entrywise", bool(cases["layer_2x3_t_1_3"]["validation"]["full_match"]))
    check("chain n=6 fast construction matches canonical build_R entrywise", bool(cases["chain_n6_t_1_3"]["validation"]["full_match"]))
    check(
        "2x4 actual-LCM D regression against the wave-12 artifact",
        bool(cases["layer_2x4_t_1_3"]["validation"]["D_match"]),
        f"D={cases['layer_2x4_t_1_3']['D']}",
    )
    check(
        "chain n=8 actual-LCM D regression against the wave-12 artifact",
        bool(cases["chain_n8_t_1_3"]["validation"]["D_match"]),
        f"D={cases['chain_n8_t_1_3']['D']}",
    )
    check("3x3 integral model anchored to direct Fraction entries", bool(cases["layer_3x3_t_1_3"]["validation"]["all_match"]))
    check("chain n=9 integral model anchored to direct Fraction entries", bool(cases["chain_n9_t_1_3"]["validation"]["all_match"]))
    check(
        "wave-12 2x3 t=1/3 gcd regression is 385 at both primes",
        cases["layer_2x3_t_1_3"]["modular_gcd_degrees"] == [385, 385],
        f"got {cases['layer_2x3_t_1_3']['modular_gcd_degrees']}",
    )
    check(
        "wave-12 2x4 t=1/3 gcd regression is 9329 at both primes",
        cases["layer_2x4_t_1_3"]["modular_gcd_degrees"] == [9329, 9329],
        f"got {cases['layer_2x4_t_1_3']['modular_gcd_degrees']}",
    )
    for name in ("layer_2x3_t_1_3", "layer_2x4_t_1_3"):
        if name in wave12_cases:
            for p in PRIMES:
                reference = wave12_cases[name]["per_prime"][str(p)]
                check(
                    f"{name}@{p}: charpoly/pair/gcd digests replay the wave-12 values",
                    cases[name]["per_prime"][str(p)]["charpoly_sha256"] == reference["charpoly_sha256"]
                    and cases[name]["per_prime"][str(p)]["pair_poly_sha256"] == reference["pair_poly_sha256"]
                    and cases[name]["per_prime"][str(p)]["gcd_sha256"] == reference["gcd_sha256"],
                )

    layer33 = cases["layer_3x3_t_1_3"]
    floor9 = int(layer33["gaussian_gcd_floor"])
    bound9 = int(layer33["modular_gcd_degree_upper_bound"])
    check(
        "DECISIVE: 3x3 layer modular gcd upper bound is strictly below the Gaussian floor 111645",
        bound9 < floor9,
        f"min gcd = {bound9}, distinct >= {int(layer33['distinct_pair_products_lower_bound'])} > 19171",
    )
    check(
        "3x3 distinct-pair lower bound exceeds the nine-mode Gaussian maximum 19171",
        int(layer33["distinct_pair_products_lower_bound"]) > int(layer33["gaussian_max_distinct_pair_products"]),
    )
    generic_theorem = generic_3x3_bound(bound9, model33)
    check(
        "generic 3x3 bookkeeping: distinct-pair lower bound exceeds 19171",
        generic_theorem["generic_distinct_pair_products_lower_bound"]
        > generic_theorem["gaussian_max_distinct_pair_products"],
    )
    check(
        "generic 3x3 entry-degree bound uses the measured e range [-6,6], F=6, and h-sum bound 18",
        generic_theorem["e_range"] == [-6, 6] and model33.F == 6 and model33.max_hamming_sum == 18,
    )

    elapsed = time.perf_counter() - started
    artifact = {
        "provenance": {
            "script": "experiments/e122_pair_product_3x3.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "blas": "float64 matrix products route through BLAS dgemm; exactness is guaranteed by the asserted all-partial-sums-below-2^53 guard and audited against int64 products and int64 binary powers",
            "precision": "exact Fraction anchors and exact integer integralization; exact finite-field arithmetic only after M = W/G. Floating point enters only through exactly-representable integer partial sums below 2^53.",
            "construction": "canonical e38_gaussianity_certificate.build_R definition computed through the split-integral route W = P~ diag(d~) P~; D is the actual LCM N0/G of every resulting Fraction denominator.",
        },
        "data": {
            "status_tags": {
                "theorem": "[THEOREM] the t=1/3 modular certificate plus the Q(t) specialization/resultant argument proves a full nine-mode no-go for the open 3x3 layer at t=1/3 and for every positive t outside a finite exceptional set.",
                "computation": "[COMPUTATION] modular values are exact values in F_p; by the monic-integral lemma each is an upper bound on the characteristic-zero gcd degree at that named rational point.",
                "unresolved": "[UNRESOLVED] neither the generic gcd polynomial nor the exceptional resultant is expanded; exact exceptional parameters, parity-projected sectors, and thermodynamic limits remain undecided.",
            },
            "threshold_derivation": "For a full n-mode subset-product spectrum, distinct eigenvalue slots yield pair-product exponents in {0,1,2}^n with at least one coordinate 1. Thus at most 3^n-2^n pair values occur; C(2^n,2)-(3^n-2^n) is the necessary characteristic-zero gcd floor.",
            "integral_model": "R = W/N0 with N0 = den(t)^{2n} num(q)^E den(q)^F; the actual LCM of the entry denominators is D = N0/G with G = gcd(gcd_ij W_ij, N0), so M = D R = W/G is integral and charpoly(wedge^2 M) in Z[z] is monic. Nonzero scaling preserves pair-product equality multiplicities, and reduction at every prime keeps the monic char-zero gcd degree below the modular one.",
            "fast_construction_guards": {
                "split": "W = W_hi * 2^26 + W_lo with the left factor split before the int64 matmul; asserted |left_hi|*max(P~)*dim and (2^26-1)*max(P~)*dim below 2^63, so both int64 matmuls are exact.",
                "float64": "all BLAS partial sums are integers below 2^53 (asserted per use); int64 binary powers and one full int64 product audit the route.",
                "lazy_division": "each window entry drifts by at most one (p-1)^2 per eliminated coefficient, giving p + N*(p-1)^2 < 2^63 (asserted); the tail is reduced once.",
            },
            "pipeline": [
                "s_k = tr(M^k), k = 1..dim, by float64-BLAS modular matrix products with per-step reduction",
                "Newton identities reconstruct the monic degree-dim characteristic polynomial; float64 Horner checks Cayley-Hamilton",
                "the charpoly recurrence extends s_k through 2*C(dim,2), audited by int64 binary matrix powers",
                "q_m = (s_m^2 - s_{2m})/2 reconstructs pair-product traces",
                "Newton identities reconstruct the monic pair polynomial C_2 = charpoly(wedge^2 M)",
                "lazy-reduction Euclidean gcd(C_2, C_2') returns its degree and is checked to divide both inputs",
            ],
            "primes": {
                "values": list(PRIMES),
                "primality": "trial division through integer sqrt",
                "constraints": "p is odd, exceeds 2*130816, and avoids 3 and 5 so that every Newton divisor, 2, and G invert.",
                "guards": guards,
            },
            "cases": cases,
            "generic_coupling_theorem_3x3": generic_theorem,
            "wave12_3x3_wall_reference": {
                "status": "context for this experiment",
                "stage": "direct power traces",
                "wall_seconds": 90,
                "elapsed_seconds": 90.058,
                "note": "wave-12 e118 hit its predeclared wall inside the int64 direct power-trace stage; this experiment replaces that stage and the Fraction build with exact faster routes and completes the certificate.",
            },
            "resource_policy": {
                "rss_limit_bytes": RSS_LIMIT_BYTES,
                "predeclared_stage_wall_seconds": STAGE_WALL_SECONDS,
                "ordinary_case_wall_seconds": CASE_WALL_SECONDS,
                "wall_seconds_total": round(elapsed, 3),
                "ru_maxrss_bytes": max_rss_bytes(),
                "observed_vs_preflight": "every timing_seconds and ru_maxrss_bytes field below is a measured observation on this machine; the stage walls are predeclared budgets, not measurements.",
            },
            "inference": "For each layer case, min_p deg gcd_Fp is an upper bound on deg gcd_Q. Hence the characteristic-zero number of distinct pair products is at least pair_count - min_p(deg gcd_Fp). If that lower bound exceeds 3^9-2^9 = 19171, the named finite layer and rational coupling cannot have a full nine-mode Gaussian subset-product spectrum. Two primes are a cross-check only; their agreement is not an equality claim over Q.",
        },
        "checks": CHECKS,
    }
    os.makedirs(ROOT / "results" / "spectral", exist_ok=True)
    with (ROOT / "results" / "spectral" / "pair_product_3x3.json").open("w") as handle:
        json.dump(artifact, handle, indent=1)

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(
        f"PASS: 3x3 pair-product certificate landed in {elapsed:.2f} s "
        f"(gcd degrees {layer33['modular_gcd_degrees']}, distinct >= {layer33['distinct_pair_products_lower_bound']})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
