"""Ordering-free pair-product certificate for the open 2x5 layer at t=1/3.

Wave 12 (e118) proved the pair-product obstruction for the open 2x3 layer
generically in t and the open 2x4 layer at t=1/3; wave 13 (e122) landed the
open 3x3 layer at t=1/3 (dim 512, 130,816 pair slots, modular gcd 59,641 <
floor 111,645).  This experiment is the next size: the open 2x5 layer, the
first open layer with FIVE columns, n = 10 sites, dim = 1024, and
C(1024,2) = 523,776 pair slots, whose nine-mode-style Gaussian floor is
523,776 - (3^10 - 2^10) = 465,751.

Pipeline (identical in structure to e122, rescaled and re-guarded):

  * exact integral construction without Fraction triple loops: the vectorized
    Hamming kernel gives P~ = den(t)^n P_t and d~ = q^e num(q)^E den(q)^F as
    integer matrices, W = P~ diag(d~) P~ is exact after a 2^26 split of the
    left factor, and M = D R with D = N0/G, G = gcd(content W, N0),
    N0 = den(t)^{2n} num(q)^E den(q)^F — exactly the e107 actual-LCM monic
    integral model R_int = D R, carried through the later wave-12/13 pipelines.
  * power traces s_k = tr(M^k), k = 1..1024, by float64 BLAS matrix products
    with a runtime integral round-trip audit; int64 binary powers audit spots.
  * Newton identities give the degree-1024 characteristic polynomial; its
    Cayley-Hamilton trace recurrence extends s_k to 2*523,776 = 1,047,552.
  * q_m = (s_m^2 - s_{2m})/2 and Newton identities reconstruct the monic
    degree-523,776 pair polynomial C_2 = charpoly(wedge^2 M) in F_p[z].  Both
    Newton loops run on int64 dots against a fill-from-the-end descending
    coefficient buffer, so every dot is contiguous (guard: N(p-1)^2 < 2^63).
  * Euclidean gcd(C_2, C_2') with lazy int64 reduction, drift bounded by
    p + (deg+1)(p-1)^2 < 2^63 (asserted).

Controls: the open 10-site chain must return EXACTLY the Gaussian floor
465,751 at both primes; the wave-13 3x3 layer must replay gcd 59,641 with its
stored coefficient digests; the wave-12 2x3 layer must replay 385.

Resource policy (project convention): every compute gate reads
time.process_time(), never wall clock and never signal.alarm.  Each stage has
a PREDECLARED process-time budget; a budget overrun raises and aborts only its
own case, recording the observed numbers, so a stalled stage can never stall
the whole run.  RSS is capped at 6,000,000,000 bytes.

Run: PYTHONPATH=src .venv/bin/python experiments/e136_pair_product_2x5.py
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import resource
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
RSS_LIMIT_BYTES = 6_000_000_000
SPLIT = 1 << 26

# Predeclared per-stage process-time budgets (seconds), fixed before any case.
STAGE_WALL_CPU_SECONDS = {
    "build": 900,
    "traces": 900,
    "charpoly": 600,
    "recurrence": 900,
    "pair_newton": 1200,
    "gcd": 4500,
}
CASE_WALL_CPU_SECONDS = 9000

WAVE12_ARTIFACT = ROOT / "results" / "spectral" / "pair_product_scale.json"
WAVE13_ARTIFACT = ROOT / "results" / "spectral" / "pair_product_3x3.json"

CHECKS: list[dict[str, object]] = []
FAILS: list[str] = []


class ResourceWall(RuntimeError):
    """A predeclared process-time or RSS budget stopped a stage."""


@dataclass
class Budget:
    """Process-time compute gate: predeclared per-stage and per-case budgets.

    tick() is called at polling intervals inside every long loop.  A budget
    overrun raises ResourceWall with the measured elapsed numbers; the case
    driver catches it, so a stalled stage aborts only its own case.
    """

    stage_walls: dict[str, int] = field(default_factory=lambda: dict(STAGE_WALL_CPU_SECONDS))
    case_wall: int = CASE_WALL_CPU_SECONDS
    case_cpu0: float = field(default_factory=time.process_time)
    stage: str = "initialization"
    stage_cpu0: float = field(default_factory=time.process_time)
    stage_wall: int = 10**9
    polls: int = 0

    def start_stage(self, name: str) -> None:
        assert name in self.stage_walls, f"unknown stage {name}"
        self.stage = name
        self.stage_cpu0 = time.process_time()
        self.stage_wall = self.stage_walls[name]

    def stage_cpu_elapsed(self) -> float:
        return time.process_time() - self.stage_cpu0

    def tick(self) -> None:
        self.polls += 1
        cpu = time.process_time()
        stage_elapsed = cpu - self.stage_cpu0
        if stage_elapsed > self.stage_wall:
            raise ResourceWall(
                f"stage '{self.stage}' exceeded its predeclared {self.stage_wall}s "
                f"process_time budget after {stage_elapsed:.1f} CPU s"
            )
        case_elapsed = cpu - self.case_cpu0
        if case_elapsed > self.case_wall:
            raise ResourceWall(
                f"case exceeded its predeclared {self.case_wall}s process_time budget "
                f"after {case_elapsed:.1f} CPU s at stage '{self.stage}'"
            )
        if max_rss_bytes() >= RSS_LIMIT_BYTES:
            raise ResourceWall(
                f"RSS cap {RSS_LIMIT_BYTES} bytes reached at stage '{self.stage}' "
                f"(ru_maxrss {max_rss_bytes()})"
            )


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    CHECKS.append({"name": name, "ok": bool(ok), "detail": detail})
    if not ok:
        FAILS.append(name)
    return bool(ok)


def max_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if platform.system() == "Darwin" else value * 1024


@contextmanager
def measured_stage(budget: Budget, name: str):
    """Enter a named stage; exit yields {cpu, wall, rss, predeclared} record fields."""
    budget.start_stage(name)
    cpu0, wall0 = time.process_time(), time.perf_counter()
    record: dict[str, object] = {}
    try:
        yield record
    finally:
        record["predeclared_wall_cpu_seconds"] = STAGE_WALL_CPU_SECONDS[name]
        record["observed_cpu_seconds"] = round(time.process_time() - cpu0, 3)
        record["observed_wall_seconds_reference_only"] = round(time.perf_counter() - wall0, 3)
        record["ru_maxrss_bytes_after_stage"] = max_rss_bytes()


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


def hamming_matrix(n: int) -> np.ndarray:
    """Vectorized d_H(k,l) for k,l < 2^n as int8."""
    dim = 1 << n
    index = np.arange(dim, dtype=np.int64)
    xor = index[:, None] ^ index[None, :]
    return np.unpackbits(
        xor.astype(np.uint32).view(np.uint8).reshape(dim, dim, 4), axis=-1, bitorder="little"
    ).sum(-1, dtype=np.int8)


def fast_integral_model(n: int, bonds, t: Fraction, budget: Budget) -> IntegralModel:
    """Exact integer W with R = W/N0 and the actual LCM D = N0/G; no Fraction matmul.

    P~ = den(t)^n * P_t has entries num(t)^h * den(t)^(n-h); the scaled diagonal
    d~_k = q^{e_k} * num(q)^E * den(q)^F = num^{e+E} * den^{F-e} is integral for
    e in [e_min, e_max] with E = max(0, -e_min), F = max(0, e_max).  Then

        W = P~ diag(d~) P~,   N0 = den(t)^{2n} * num(q)^E * den(q)^F,

    R = W/N0 exactly, and the LCM of the denominators of R equals N0/G with
    G = gcd(gcd_{ij} W_ij, N0), so M = D*R = W/G is the monic integral model.
    """
    budget.start_stage("build")
    q, eps, e_values = geometry_exponents(n, bonds, t)
    num, den = q.numerator, q.denominator
    a, bden = t.numerator, t.denominator
    E, F = max(0, -min(e_values)), max(0, max(e_values))
    dim = 1 << n
    budget.tick()

    hamming = hamming_matrix(n)
    # P~ entries num(t)^h * den(t)^(n-h) gathered from an (n+1)-table; the
    # kernel-degree bound max(h(i,k)+h(k,j)) = 2*max h follows from the table.
    h_table = np.array([a**h * bden ** (n - h) for h in range(n + 1)], dtype=np.int64)
    powers = h_table[hamming]
    max_h_sum = 2 * int(hamming.max())
    del hamming
    budget.tick()

    diagonal = np.array([num ** (x + E) * den ** (F - x) for x in e_values], dtype=np.int64)
    left = powers * diagonal[None, :]
    # int64 split matmul: left = left_hi * SPLIT + left_lo.  Both products and
    # both dim-term partial sums stay far below 2^63 (asserted), so the plain
    # int64 matmuls below are exact integer arithmetic.
    left_hi = left >> 26
    left_lo = left & (SPLIT - 1)
    max_power = int(powers.max())
    assert int(np.abs(left_hi).max()) * max_power * dim < INT64_MAX, "split hi int64 guard violated"
    assert (SPLIT - 1) * max_power * dim < INT64_MAX, "split lo int64 guard violated"
    budget.tick()
    W_hi = left_hi @ powers
    budget.tick()
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

def float_matmul_mod(a: np.ndarray, b: np.ndarray, p: int, budget: Budget) -> np.ndarray:
    """Exact (a@b) mod p through float64 BLAS; all partial sums are integers < 2^53.

    The int64 round-trip comparison is a runtime audit that every BLAS output
    was an exactly represented integer, on top of the asserted 2^53 bound.
    """
    assert_float_guard(a.shape[0], p, "float matmul")
    product = a.astype(np.float64) @ b.astype(np.float64)
    integer = product.astype(np.int64)
    assert np.array_equal(product, integer.astype(np.float64)), "non-integral binary64 matrix product"
    budget.tick()
    return integer % p


def int_matmul_mod(a: np.ndarray, b: np.ndarray, p: int, budget: Budget) -> np.ndarray:
    assert_int64_guard(a.shape[0], p, "int64 matmul")
    return (a @ b) % p


def sequential_power_traces(matrix: np.ndarray, max_power: int, p: int, budget: Budget) -> np.ndarray:
    """[0, tr(M), ..., tr(M^max_power)] by float64-BLAS modular matrix products."""
    dim = matrix.shape[0]
    assert_float_guard(dim, p, "power traces")
    base = matrix.astype(np.float64)
    power = np.eye(dim, dtype=np.float64)
    traces = np.zeros(max_power + 1, dtype=np.int64)
    for exponent in range(1, max_power + 1):
        if exponent & 63 == 0:
            budget.tick()
        power = power @ base
        reduced = power.astype(np.int64)
        assert np.array_equal(power, reduced.astype(np.float64)), "non-integral binary64 trace power"
        reduced %= p
        traces[exponent] = int(reduced.trace()) % p
        power = reduced.astype(np.float64)
    return traces


def binary_power_trace(matrix: np.ndarray, exponent: int, p: int, budget: Budget, *, route: str) -> int:
    """Independent binary-powering trace audit ('int64' exact, 'float64' BLAS)."""
    dim = matrix.shape[0]
    if route == "int64":
        assert_int64_guard(dim, p, "int64 binary power")
        result = np.eye(dim, dtype=np.int64)
        base = matrix.copy()
        e = exponent
        while e:
            if e & 1:
                result = int_matmul_mod(result, base, p, budget)
            e >>= 1
            if e:
                base = int_matmul_mod(base, base, p, budget)
        return int(result.trace()) % p
    assert route == "float64"
    assert_float_guard(dim, p, "float64 binary power")
    result = np.eye(dim, dtype=np.float64)
    base = matrix.astype(np.float64)
    e = exponent
    while e:
        if e & 1:
            result = float_matmul_mod(result, base, p, budget).astype(np.float64)
        e >>= 1
        if e:
            base = float_matmul_mod(base, base, p, budget).astype(np.float64)
    return int(result.astype(np.int64).trace()) % p


def horner_matrix_zero(charpoly_desc: np.ndarray, matrix: np.ndarray, p: int, budget: Budget) -> bool:
    """Cayley-Hamilton validation chi(M) = 0 in F_p^{dim x dim} via float64 BLAS."""
    dim = matrix.shape[0]
    assert_float_guard(dim, p, "Horner")
    identity_f = np.eye(dim, dtype=np.float64)
    value = np.eye(dim, dtype=np.float64)
    for coefficient in charpoly_desc[1:]:
        value = value @ matrix
        value = value + (int(coefficient) * identity_f)
        integer = value.astype(np.int64)
        assert np.array_equal(value, integer.astype(np.float64)), "non-integral binary64 Horner step"
        value = (integer % p).astype(np.float64)
        budget.tick()
    return bool(np.all(value == 0))


# ---------------------------------------------------------------- Newton machinery

def newton_from_power_sums(
    power_sums: np.ndarray, degree: int, p: int, inverses: list[int], budget: Budget
) -> np.ndarray:
    """DESCENDING [1,a1,...,a_degree] from power sums via Newton identities.

    Internally the coefficients live in a fill-from-the-end ascending buffer
    desc[degree-j] = a_j, so that every int64 dot touches two contiguous
    slices (guard degree*(p-1)^2 < 2^63 asserted by the caller); the buffer is
    reversed once on return to the descending convention.
    """
    assert_int64_guard(degree, p, "Newton dot product")
    desc = np.zeros(degree + 1, dtype=np.int64)
    desc[degree] = 1  # a_0 = 1 lives at index degree
    sums = np.ascontiguousarray(power_sums[1 : degree + 1])
    for k in range(1, degree + 1):
        if k & 2047 == 0:
            budget.tick()
        # sum_{i=1..k} a_{k-i} s_i pairs desc[degree-k+i] = a_{k-i} with s_i.
        total = int(np.dot(desc[degree - k + 1 :], sums[:k])) % p
        desc[degree - k] = (-total * inverses[k]) % p
    return np.ascontiguousarray(desc[::-1])


def extend_power_sums(
    traces: np.ndarray, charpoly_desc: np.ndarray, upto: int, p: int, budget: Budget
) -> np.ndarray:
    """Cayley-Hamilton trace recurrence s_k = -sum a_i s_{k-i} beyond matrix dimension.

    The history is stored reversed (index capacity-1-m holds s_m) so that each
    int64 dot pairs the contiguous a_1..a_d window with s_{k-1}..s_{k-d}.
    """
    dim = charpoly_desc.size - 1
    assert_int64_guard(dim, p, "trace recurrence")
    coefficients = np.ascontiguousarray(charpoly_desc[1:])
    out = np.zeros(upto + 1, dtype=np.int64)
    out[: traces.size] = traces
    capacity = upto + 1
    history = np.zeros(capacity, dtype=np.int64)  # history[capacity-1-m] = s_m
    history[capacity - traces.size :] = traces[::-1]
    for k in range(dim + 1, upto + 1):
        if k & 8191 == 0:
            budget.tick()
        total = int(np.dot(coefficients, history[capacity - k : capacity - k + dim])) % p
        value = (-total) % p
        out[k] = value
        history[capacity - 1 - k] = value
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


def remainder_lazy(dividend: np.ndarray, divisor: np.ndarray, p: int, budget: Budget) -> np.ndarray:
    """Remainder of division in F_p[z] with lazy int64 reduction.

    Both inputs are descending int64 arrays with values in [0, p); the divisor
    need not be monic.  Each eliminated coefficient subtracts lead*divisor from
    its window WITHOUT reducing the window mod p.  One entry can be touched by
    up to deg(divisor)+1 windows, so the drift bound is
    p + (deg(divisor)+1)*(p-1)^2 < 2^63 (asserted), and the tail is reduced
    once at the end.
    """
    dividend_degree = dividend.size - 1
    divisor_degree = divisor.size - 1
    touches_per_entry = min(divisor_degree + 1, dividend_degree - divisor_degree + 1)
    assert p + max(touches_per_entry, 0) * (p - 1) ** 2 < INT64_MAX, "lazy-division drift guard violated"
    if dividend_degree < divisor_degree:
        return poly_trim(dividend)
    work = dividend.copy()
    lead_inv = pow(int(divisor[0]), -1, p)
    scratch = np.empty(divisor_degree + 1, dtype=np.int64)
    for offset in range(dividend_degree - divisor_degree + 1):
        if offset & 1023 == 0:
            budget.tick()
        top = int(work[offset]) % p
        if top:
            lead = top * lead_inv % p
            np.multiply(divisor, lead, out=scratch)
            work[offset : offset + divisor_degree + 1] -= scratch
    return poly_trim(work[dividend_degree - divisor_degree + 1 :] % p)


def poly_gcd(a: np.ndarray, b: np.ndarray, p: int, budget: Budget) -> tuple[np.ndarray, int]:
    """Euclidean gcd in F_p[z] with lazy remainders; the result is monic."""
    a = poly_trim(np.asarray(a, dtype=np.int64) % p)
    b = poly_trim(np.asarray(b, dtype=np.int64) % p)
    steps = 0
    while not poly_is_zero(b):
        budget.tick()
        a, b = b, remainder_lazy(a, b, p, budget)
        steps += 1
    return poly_monic(a, p), steps


def derivative_desc(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    derivative_ascending = (np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p
    return np.ascontiguousarray(derivative_ascending[::-1])


def coefficient_digest(ascending_coefficients) -> str:
    return hashlib.sha256(",".join(str(int(c)) for c in ascending_coefficients).encode()).hexdigest()


# ------------------------------------------------------------------------ thresholds

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

def run_pipeline(matrix: np.ndarray, p: int, budget: Budget) -> dict[str, object]:
    """Trace/Newton/recurrence/pair-Newton/gcd pipeline at one prime."""
    started = time.perf_counter()
    dim = matrix.shape[0]
    pairs = dim * (dim - 1) // 2
    # p must exceed the polynomial degree N = C(dim,2): every Newton divisor k
    # and every derivative exponent 1..N must be a unit mod p.  The extended
    # trace index 2*N is an ARRAY INDEX, never a value mod p, so p > N is the
    # exact requirement (wave-12/13 wrote p > 2N, which was sufficient but not
    # necessary at their sizes).
    assert is_prime(p) and p > pairs and p % 2 == 1
    assert_int64_guard(dim, p, "int64 binary powers")
    assert_int64_guard(pairs, p, "Newton dot product")
    assert_float_guard(dim, p, "float64 matmul")
    inverses = [0] + [pow(k, -1, p) for k in range(1, pairs + 1)]

    stage_records: dict[str, dict[str, object]] = {}
    with measured_stage(budget, "traces") as record:
        # One full int64 product cross-checks the float64 BLAS path on this matrix.
        float_check = np.array_equal(
            float_matmul_mod(matrix, matrix, p, budget), int_matmul_mod(matrix, matrix, p, budget)
        )
        traces = sequential_power_traces(matrix, dim, p, budget)
    stage_records["traces"] = record

    with measured_stage(budget, "charpoly") as record:
        characteristic = newton_from_power_sums(traces, dim, p, inverses, budget)
        horner_ok = horner_matrix_zero(characteristic, matrix.astype(np.float64), p, budget)
    stage_records["charpoly"] = record

    with measured_stage(budget, "recurrence") as record:
        full_traces = extend_power_sums(traces, characteristic, 2 * pairs, p, budget)
        # Audited recurrence spots, clamped to the extended range 2*pairs.  Two
        # int64 binary powers anchor the whole float64 route; the remaining
        # spots ride the audited float64 powers.
        int64_spots = tuple(sorted({min(x, 2 * pairs) for x in (dim + 1, 65536)}))
        float_spots = tuple(sorted({min(x, 2 * pairs) for x in (2 * dim, 65537, pairs, 2 * pairs)}))
        spots = []
        for exponent in int64_spots:
            direct = binary_power_trace(matrix, exponent, p, budget, route="int64")
            spots.append(
                {
                    "k": exponent,
                    "route": "int64",
                    "direct": direct,
                    "recurrence": int(full_traces[exponent]),
                    "match": direct == int(full_traces[exponent]),
                }
            )
        for exponent in float_spots:
            direct = binary_power_trace(matrix, exponent, p, budget, route="float64")
            spots.append(
                {
                    "k": exponent,
                    "route": "float64",
                    "direct": direct,
                    "recurrence": int(full_traces[exponent]),
                    "match": direct == int(full_traces[exponent]),
                }
            )
    stage_records["recurrence"] = record

    with measured_stage(budget, "pair_newton") as record:
        inverse_two = pow(2, -1, p)
        first = full_traces[1 : pairs + 1]
        second = full_traces[2 : 2 * pairs + 1 : 2]
        assert first.size == pairs and second.size == pairs
        body = ((first * first - second) % p) * inverse_two % p  # int64: (p-1)^2 < 2^42
        pair_sums = np.concatenate(([0], body))  # index 0 is the Newton placeholder
        pair_poly = newton_from_power_sums(pair_sums, pairs, p, inverses, budget)
    stage_records["pair_newton"] = record

    with measured_stage(budget, "gcd") as record:
        derivative = derivative_desc(pair_poly, p)
        common, gcd_steps = poly_gcd(pair_poly, derivative, p, budget)
        after_gcd_algorithm = time.perf_counter()
        # Belt-and-braces: the reported monic gcd must divide both inputs exactly.
        divides = poly_is_zero(remainder_lazy(pair_poly, common, p, budget)) and poly_is_zero(
            remainder_lazy(derivative, common, p, budget)
        )
        gcd_check_seconds = round(time.perf_counter() - after_gcd_algorithm, 3)
    stage_records["gcd"] = record

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
        "gcd_division_checks_wall_seconds": gcd_check_seconds,
        "gcd_euclid_steps": gcd_steps,
        "float64_vs_int64_first_product_equal": bool(float_check),
        "horner_charpoly_zero_matrix": bool(horner_ok),
        "spot_trace_checks": spots,
        "spot_trace_checks_all_match": all(spot["match"] for spot in spots),
        "timing_seconds_wall_reference_only": {
            "traces": stage_records["traces"]["observed_wall_seconds_reference_only"],
            "characteristic_newton_and_horner": stage_records["charpoly"]["observed_wall_seconds_reference_only"],
            "trace_recurrence_and_spots": stage_records["recurrence"]["observed_wall_seconds_reference_only"],
            "pair_newton": stage_records["pair_newton"]["observed_wall_seconds_reference_only"],
            "euclidean_gcd": stage_records["gcd"]["observed_wall_seconds_reference_only"],
            "total": round(time.perf_counter() - started, 3),
        },
        "stage_resources": stage_records,
        "ru_maxrss_bytes": max_rss_bytes(),
    }


def build_and_run(
    *,
    name: str,
    geometry: str,
    n_sites: int,
    bonds: list[tuple[int, int]],
    t: Fraction,
    validation: dict[str, object] | None = None,
) -> dict[str, object]:
    """Fast exact integral model plus both-prime reconstruction.

    A predeclared-stage budget overrun aborts THIS case only, recording the
    blocking stage and its measured numbers; the run continues.
    """
    threshold = case_threshold(n_sites)
    budget = Budget()
    record: dict[str, object] = {
        "name": name,
        "geometry": geometry,
        "t": str(t),
        **threshold,
    }
    print(f"case {name} ...", flush=True)
    try:
        with measured_stage(budget, "build") as build_record:
            model = fast_integral_model(n_sites, bonds, t, budget)
        build_record["observed_N0"] = str(model.N0)
        record.update(
            {
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
                "build_resource": dict(build_record),
            }
        )
        per_prime: dict[str, object] = {}
        for p in PRIMES:
            budget.start_stage("build")  # matrix reduction rides the build budget
            matrix_started = time.perf_counter()
            matrix = matrix_mod_p(model, p)
            matrix_mod_elapsed = time.perf_counter() - matrix_started
            print(f"  pipeline p={p} ...", flush=True)
            result = run_pipeline(matrix, p, budget)
            result["matrix_mod_p_seconds"] = round(matrix_mod_elapsed, 3)
            result["matrix_mod_p_ru_maxrss_bytes_after_construction"] = max_rss_bytes()
            per_prime[str(p)] = result
            print(
                f"    gcd degree {result['gcd_degree']} (Euclid steps {result['gcd_euclid_steps']})",
                flush=True,
            )
            del matrix
        gcd_degrees = [int(per_prime[str(p)]["gcd_degree"]) for p in PRIMES]
        modular_upper_bound = min(gcd_degrees)
        record.update(
            {
                "per_prime": per_prime,
                "modular_gcd_degrees": gcd_degrees,
                "modular_gcd_degree_upper_bound": modular_upper_bound,
                "distinct_pair_products_lower_bound": int(threshold["pair_count"]) - modular_upper_bound,
                "decisive_non_gaussian": bool(
                    modular_upper_bound < int(threshold["gaussian_gcd_floor"])
                ),
                "at_gaussian_floor": bool(
                    modular_upper_bound == int(threshold["gaussian_gcd_floor"])
                ),
                "ru_maxrss_bytes_after_case": max_rss_bytes(),
                "case_cpu_seconds": round(time.process_time() - budget.case_cpu0, 3),
                "blocked_stage": None,
            }
        )
    except (ResourceWall, AssertionError) as stop:
        record.update(
            {
                "blocked_stage": budget.stage,
                "blocked_detail": str(stop),
                "observed_cpu_seconds_at_block": round(budget.stage_cpu_elapsed(), 3),
                "predeclared_wall_cpu_seconds": STAGE_WALL_CPU_SECONDS.get(budget.stage),
                "ru_maxrss_bytes_at_block": max_rss_bytes(),
                "case_cpu_seconds": round(time.process_time() - budget.case_cpu0, 3),
            }
        )
        print(f"  BLOCKED at stage '{budget.stage}': {stop}", flush=True)
    if validation is not None:
        record["validation"] = validation
    return record


# ------------------------------------------------------------- generic-t bookkeeping

def generic_layer_bound(generic_gcd_upper_bound: int, model: IntegralModel, n_sites: int) -> dict[str, object]:
    """Effective degree bookkeeping for the Q(t) specialization theorem.

    M(t) = (2t)^{2c} q(t)^c R(t) is an integer-polynomial matrix when the
    clearing exponent satisfies c >= max(-e_min, e_max), because
    (2t)^{2c} q^{c+e} = (2t)^{c-e} (1+t^2)^{c+e} and both exponents are then
    nonnegative for every e in the measured range.  For the open 2x5 layer the
    measured range is e in [-7,6], so c = 7 (= E), NOT F = 6: the wave-13 3x3
    value c = F was valid there only because its range was symmetric.
    """
    threshold = case_threshold(n_sites)
    pair_degree = int(threshold["pair_count"])
    e_min, e_max = min(model.e_values), max(model.e_values)
    c = max(-e_min, e_max)
    assert c == max(model.E, model.F)
    assert all(c - e >= 0 and c + e >= 0 for e in (e_min, e_max)), (
        "clearing exponent does not cover the measured e range"
    )
    # scalar degree s(e) = (c-e) + 2(c+e) = 3c+e is maximized at e_max;
    # entries also carry kernel degrees h1+h2 <= 2n.
    scalar_degree = 3 * c + e_max
    entry_degree = 2 * model.n_sites + scalar_degree
    wedge_entry_degree = 2 * entry_degree
    resultant_degree_bound = 2 * wedge_entry_degree * pair_degree * (pair_degree - 1)
    return {
        "status": "THEOREM",
        "domain": "all positive real t outside a finite exceptional set E",
        "clearing_exponent_c": c,
        "polynomial_scalar": f"(2t)^{2 * c} q(t)^{c}",
        "e_range": [e_min, e_max],
        "scalar_t_degree_upper_bound": scalar_degree,
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
            "ten-mode Gaussian subset-product spectrum. No value, location, or minimal "
            "polynomial of an exceptional parameter is claimed."
        ),
    }


# ------------------------------------------------------------------------------ driver

def validate_against_canonical(name: str, n: int, bonds, t: Fraction, model: IntegralModel) -> dict[str, object]:
    """Full-entry comparison for dim <= 64; Fraction spot anchors above that."""
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
    dim = 1 << n
    rng = np.random.default_rng(20260817)
    spot_entries = [(0, j) for j in range(0, dim, 37)] + [
        (int(i), int(j)) for i, j in zip(rng.integers(0, dim, 16), rng.integers(0, dim, 16))
    ]
    spot_ok = spot_fraction_check(model, spot_entries)
    return {
        "mode": "first-row stride plus 16 seeded random entries against direct Fraction row sums",
        "checked_entries": len(spot_entries),
        "all_match": all(spot_ok),
        "failed_entries": [e for e, ok in zip(spot_entries, spot_ok) if not ok],
    }


def main() -> int:
    started = time.perf_counter()
    budget = Budget()

    check("both declared modular primes are prime", all(is_prime(p) for p in PRIMES))
    check(
        "2x5 threshold is C(1024,2)=523776, 3^10-2^10=58025, floor=465751",
        case_threshold(10)
        == {
            "n_sites": 10,
            "dim": 1024,
            "pair_count": 523776,
            "gaussian_max_distinct_pair_products": 58025,
            "gaussian_gcd_floor": 465751,
        },
    )
    check("case_threshold(6) matches the wave-12 2x3 values", case_threshold(6)["gaussian_gcd_floor"] == 1351)
    check("case_threshold(9) matches the wave-13 3x3 values", case_threshold(9)["gaussian_gcd_floor"] == 111645)

    guards: dict[str, object] = {}
    for n_sites in (6, 9, 10):
        threshold = case_threshold(n_sites)
        guards[str(n_sites)] = {
            str(p): {
                "dim_times_pm1_sq": threshold["dim"] * (p - 1) ** 2,
                "pairs_times_pm1_sq": threshold["pair_count"] * (p - 1) ** 2,
                "two_pairs_times_pm1_sq": 2 * threshold["pair_count"] * (p - 1) ** 2,
                "lazy_drift_bound": p + (threshold["pair_count"] + 1) * (p - 1) ** 2,
                "int64_max": INT64_MAX,
                "float_exact_max": FLOAT_EXACT_MAX,
                "safe_int64": threshold["dim"] * (p - 1) ** 2 < INT64_MAX
                and threshold["pair_count"] * (p - 1) ** 2 < INT64_MAX
                and 2 * threshold["pair_count"] * (p - 1) ** 2 < INT64_MAX
                and p + (threshold["pair_count"] + 1) * (p - 1) ** 2 < INT64_MAX,
                "safe_float64_matmul": threshold["dim"] * (p - 1) ** 2 < FLOAT_EXACT_MAX,
            }
            for p in PRIMES
        }
    check(
        "all int64 matmul/Newton/lazy-drift and float64-BLAS exactness guards hold",
        all(v["safe_int64"] and v["safe_float64_matmul"] for per_n in guards.values() for v in per_n.values()),
    )
    check(
        "primes exceed the degree 523776 (so Newton divisors and derivative exponents are units) and avoid 3 and 5",
        all(p > 523776 and p % 3 and p % 5 for p in PRIMES),
    )

    wave12_cases = json.loads(WAVE12_ARTIFACT.read_text())["data"]["cases"] if WAVE12_ARTIFACT.exists() else {}
    wave13_cases = json.loads(WAVE13_ARTIFACT.read_text())["data"]["cases"] if WAVE13_ARTIFACT.exists() else {}

    cases: dict[str, dict[str, object]] = {}
    t = Fraction(1, 3)

    # Cheap warm-up control with a stored wave-12 digest chain.
    bonds23 = list(layer_bonds((2, 3), (False, False)))
    model23 = fast_integral_model(6, bonds23, t, budget)
    validation23 = validate_against_canonical("layer_2x3_t_1_3", 6, bonds23, t, model23)
    del model23
    cases["layer_2x3_t_1_3"] = build_and_run(
        name="layer_2x3_t_1_3", geometry="open 2x3 layer", n_sites=6, bonds=bonds23, t=t,
        validation=validation23,
    )

    # Wave-13 digest replay control.
    bonds33 = list(layer_bonds((3, 3), (False, False)))
    model33 = fast_integral_model(9, bonds33, t, budget)
    validation33 = validate_against_canonical("layer_3x3_t_1_3", 9, bonds33, t, model33)
    del model33
    cases["layer_3x3_t_1_3"] = build_and_run(
        name="layer_3x3_t_1_3", geometry="open 3x3 layer (wave-13 digest replay)", n_sites=9,
        bonds=bonds33, t=t, validation=validation33,
    )

    # The Gaussian-floor control at the decisive size: the open 10-site chain.
    bonds_n10 = [(site, site + 1) for site in range(9)]
    model_n10 = fast_integral_model(10, bonds_n10, t, budget)
    validation_n10 = validate_against_canonical("chain_n10_t_1_3", 10, bonds_n10, t, model_n10)
    del model_n10
    cases["chain_n10_t_1_3"] = build_and_run(
        name="chain_n10_t_1_3", geometry="open chain n=10 (Gaussian floor control)", n_sites=10,
        bonds=bonds_n10, t=t, validation=validation_n10,
    )

    # The decisive case: the open 2x5 layer.
    bonds25 = list(layer_bonds((2, 5), (False, False)))
    degree_25 = sorted(sum(site in bond for bond in bonds25) for site in range(10))
    check(
        "2x5 target is the open 13-bond layer with degree profile 2^4 3^6",
        len(bonds25) == 13 and degree_25 == [2, 2, 2, 2, 3, 3, 3, 3, 3, 3],
        f"bonds={len(bonds25)}, degrees={degree_25}",
    )
    model25 = fast_integral_model(10, bonds25, t, budget)
    validation25 = validate_against_canonical("layer_2x5_t_1_3", 10, bonds25, t, model25)
    del model25
    cases["layer_2x5_t_1_3"] = build_and_run(
        name="layer_2x5_t_1_3",
        geometry="open 2x5 layer (first open layer with five columns)",
        n_sites=10, bonds=bonds25, t=t, validation=validation25,
    )

    # ------------------------------------------------------------------- checks
    for key, record in cases.items():
        if record.get("blocked_stage") is not None:
            check(
                f"{key}: blocked stage reports measured process_time and RSS",
                record["observed_cpu_seconds_at_block"] > 0 and record["ru_maxrss_bytes_at_block"] > 0,
                f"stage {record['blocked_stage']} after {record['observed_cpu_seconds_at_block']} CPU s",
            )
            continue
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

    check(
        "2x3 fast construction matches canonical build_R entrywise",
        bool(cases["layer_2x3_t_1_3"]["validation"]["full_match"]),
    )
    check("3x3 integral model anchored to direct Fraction entries", bool(cases["layer_3x3_t_1_3"]["validation"]["all_match"]))
    check("chain n=10 integral model anchored to direct Fraction entries", bool(cases["chain_n10_t_1_3"]["validation"]["all_match"]))
    check("2x5 integral model anchored to direct Fraction entries", bool(cases["layer_2x5_t_1_3"]["validation"]["all_match"]))
    check(
        "wave-12 2x3 t=1/3 gcd regression is 385 at both primes",
        cases["layer_2x3_t_1_3"].get("modular_gcd_degrees") == [385, 385],
        f"got {cases['layer_2x3_t_1_3'].get('modular_gcd_degrees')}",
    )
    check(
        "wave-13 3x3 t=1/3 gcd regression is 59641 at both primes",
        cases["layer_3x3_t_1_3"].get("modular_gcd_degrees") == [59641, 59641],
        f"got {cases['layer_3x3_t_1_3'].get('modular_gcd_degrees')}",
    )
    for name, reference_cases, label in (
        ("layer_2x3_t_1_3", wave12_cases, "wave-12"),
        ("layer_3x3_t_1_3", wave13_cases, "wave-13"),
    ):
        if name in reference_cases:
            for p in PRIMES:
                reference = reference_cases[name]["per_prime"][str(p)]
                check(
                    f"{name}@{p}: charpoly/pair/gcd digests replay the {label} values",
                    cases[name]["per_prime"][str(p)]["charpoly_sha256"] == reference["charpoly_sha256"]
                    and cases[name]["per_prime"][str(p)]["pair_poly_sha256"] == reference["pair_poly_sha256"]
                    and cases[name]["per_prime"][str(p)]["gcd_sha256"] == reference["gcd_sha256"],
                )

    layer25 = cases["layer_2x5_t_1_3"]
    decisive_value: int | None = None
    if layer25.get("blocked_stage") is None:
        floor10 = int(layer25["gaussian_gcd_floor"])
        bound10 = int(layer25["modular_gcd_degree_upper_bound"])
        decisive_value = bound10
        check(
            "DECISIVE: 2x5 layer modular gcd upper bound is strictly below the Gaussian floor 465751",
            bound10 < floor10,
            f"min gcd = {bound10}, distinct >= {int(layer25['distinct_pair_products_lower_bound'])} > 58025",
        )
        check(
            "2x5 distinct-pair lower bound exceeds the ten-mode Gaussian maximum 58025",
            int(layer25["distinct_pair_products_lower_bound"])
            > int(layer25["gaussian_max_distinct_pair_products"]),
        )
        # Generic-t bookkeeping on the same measured model.
        model_for_generic = fast_integral_model(10, list(layer_bonds((2, 5), (False, False))), t, budget)
        generic_theorem = generic_layer_bound(bound10, model_for_generic, 10)
        check(
            "generic 2x5 bookkeeping: distinct-pair lower bound exceeds 58025",
            generic_theorem["generic_distinct_pair_products_lower_bound"]
            > generic_theorem["gaussian_max_distinct_pair_products"],
        )
        check(
            "generic 2x5 clearing exponent is c = max(E,F) = 7 over the measured e range [-7,6]",
            generic_theorem["clearing_exponent_c"] == 7
            and generic_theorem["e_range"] == [-7, 6]
            and model_for_generic.E == 7
            and model_for_generic.F == 6
            and model_for_generic.max_hamming_sum == 20,
        )
        del model_for_generic
    else:
        generic_theorem = None

    elapsed = time.perf_counter() - started
    artifact = {
        "provenance": {
            "script": "experiments/e136_pair_product_2x5.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "blas": "float64 matrix products route through BLAS dgemm; exactness is guaranteed by the asserted all-partial-sums-below-2^53 guard and audited by int64 round-trip equality, one full int64 product, and int64 binary powers",
            "precision": "exact Fraction anchors and exact integer integralization; exact finite-field arithmetic only after M = W/G. Floating point enters only through exactly-representable integer partial sums below 2^53.",
            "construction": "canonical e38_gaussianity_certificate.build_R definition computed through the split-integral route W = P~ diag(d~) P~; D is the actual LCM N0/G of every resulting Fraction denominator.",
        },
        "data": {
            "status_tags": {
                "theorem": "[THEOREM] the t=1/3 modular certificate plus the Q(t) specialization/resultant argument proves a full ten-mode no-go for the open 2x5 layer at t=1/3 and for every positive t outside a finite exceptional set.",
                "computation": "[COMPUTATION] modular values are exact values in F_p; by the monic-integral lemma each is an upper bound on the characteristic-zero gcd degree at that named rational point.",
                "unresolved": "[UNRESOLVED] neither the generic gcd polynomial nor the exceptional resultant is expanded; exact exceptional parameters, parity-projected sectors, and thermodynamic limits remain undecided.",
            },
            "threshold_derivation": "For a full n-mode subset-product spectrum, distinct eigenvalue slots yield pair-product exponents in {0,1,2}^n with at least one coordinate 1. Thus at most 3^n-2^n pair values occur; C(2^n,2)-(3^n-2^n) is the necessary characteristic-zero gcd floor.",
            "integral_model": "R = W/N0 with N0 = den(t)^{2n} num(q)^E den(q)^F; the actual LCM of the entry denominators is D = N0/G with G = gcd(gcd_ij W_ij, N0), so M = D R = W/G is integral and charpoly(wedge^2 M) in Z[z] is monic. Nonzero scaling preserves pair-product equality multiplicities, and reduction at every prime keeps the monic char-zero gcd degree below the modular one.",
            "fast_construction_guards": {
                "split": "W = W_hi * 2^26 + W_lo with the left factor split before the int64 matmul; asserted |left_hi|*max(P~)*dim and (2^26-1)*max(P~)*dim below 2^63, so both int64 matmuls are exact.",
                "float64": "all BLAS partial sums are integers below 2^53 (asserted per use) with int64 round-trip equality checks; int64 binary powers and one full int64 product audit the route.",
                "newton_buffers": "both Newton loops dot contiguous int64 slices against fill-from-the-end descending buffers; guard N*(p-1)^2 < 2^63 asserted.",
                "lazy_division": "one entry is touched by at most deg(divisor)+1 windows, so its drift is bounded by p + (deg+1)*(p-1)^2 < 2^63 (asserted); the tail is reduced once.",
            },
            "pipeline": [
                "s_k = tr(M^k), k = 1..1024, by float64-BLAS modular matrix products with per-step int64 round-trip audits",
                "Newton identities reconstruct the monic degree-1024 characteristic polynomial; float64 Horner checks Cayley-Hamilton",
                "the charpoly recurrence extends s_k through 2*C(1024,2) = 1047552, audited by int64 and float64 binary powers",
                "q_m = (s_m^2 - s_{2m})/2 reconstructs pair-product traces",
                "Newton identities reconstruct the monic pair polynomial C_2 = charpoly(wedge^2 M) of degree 523776",
                "lazy-reduction Euclidean gcd(C_2, C_2') returns its degree and is checked to divide both inputs",
            ],
            "primes": {
                "values": list(PRIMES),
                "primality": "trial division through integer sqrt",
                "constraints": "p is odd, exceeds the degree N=523776 so every Newton divisor and derivative exponent 1..N is a unit mod p, and avoids 3 and 5 so 2 and G invert; the extended trace index 2N=1047552 is an array index, not a residue.",
                "guards": guards,
            },
            "cases": cases,
            "generic_coupling_theorem_2x5": generic_theorem,
            "wave13_3x3_reference": {
                "status": "control replay",
                "gcd_degrees": [59641, 59641],
                "digests": "charpoly/pair/gcd SHA-256 digests replay the wave-13 artifact at both primes when that artifact is present",
            },
            "resource_policy": {
                "gate": "every compute gate reads time.process_time() (CPU seconds), never wall clock and never signal.alarm; wall-clock seconds are recorded as reference only",
                "abort_policy": "a predeclared per-stage budget overrun raises and aborts only its own case, recording the blocking stage and measured numbers; the run always completes",
                "rss_limit_bytes": RSS_LIMIT_BYTES,
                "predeclared_stage_wall_cpu_seconds": STAGE_WALL_CPU_SECONDS,
                "ordinary_case_wall_cpu_seconds": CASE_WALL_CPU_SECONDS,
                "wall_seconds_total_reference_only": round(elapsed, 3),
                "cpu_seconds_total": round(time.process_time(), 3),
                "ru_maxrss_bytes": max_rss_bytes(),
            },
            "inference": "For each layer case, min_p deg gcd_Fp is an upper bound on deg gcd_Q. Hence the characteristic-zero number of distinct pair products is at least pair_count - min_p(deg gcd_Fp). If that lower bound exceeds 3^10-2^10 = 58025, the named finite layer and rational coupling cannot have a full ten-mode Gaussian subset-product spectrum. Two primes are a cross-check only; their agreement is not an equality claim over Q.",
        },
        "checks": CHECKS,
    }
    os.makedirs(ROOT / "results" / "spectral", exist_ok=True)
    with (ROOT / "results" / "spectral" / "pair_product_2x5.json").open("w") as handle:
        json.dump(artifact, handle, indent=1)

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    decisive_message = (
        f"PASS: 2x5 pair-product certificate landed in {elapsed:.2f} s wall "
        f"(gcd degrees {layer25.get('modular_gcd_degrees')}, "
        f"distinct >= {layer25.get('distinct_pair_products_lower_bound')})."
        if decisive_value is not None
        else f"PASS (controls only): the decisive 2x5 case was blocked at stage "
        f"'{layer25.get('blocked_stage')}' after {layer25.get('observed_cpu_seconds_at_block')} CPU s."
    )
    print(decisive_message)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
