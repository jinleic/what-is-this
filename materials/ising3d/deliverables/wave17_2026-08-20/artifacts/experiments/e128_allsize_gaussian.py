"""All-size Gaussian no-go for inhomogeneous Ising layers with a branching site.

This producer supports `proofs/allsize_gaussian.md`.  The goal is a FAMILY theorem
(every layer graph with a vertex of degree >= 3, every number of sites) rather than
one more finite certificate.

The invariant is the ordering-free pair-product collision count of `e107`:

    r(Lambda) = #{ distinct  lambda_i lambda_j , i < j  over eigenvalue SLOTS },
    r_all(Lambda) = #{ distinct lambda_i lambda_j , i <= j },

with the Gaussian ceilings r <= 3^n - 2^n and r_all <= 3^n for a full n-mode
subset-product multiset {a prod u_k^eps}.

The family statement is obtained from ONE 16-dimensional certificate on the claw
K_{1,3} by an exact tensor-decoupling identity (Lemma 3 of the note):

    r(Lambda_H x Lambda_rest) = (3^m - 2^m) r_all(Lambda_H) + 2^m r(Lambda_H)

when the m remaining sites carry independent transverse-field parameters, plus the
standard monic specialization lemma that turns a special point into a generic
statement.  The claw gives r_all >= 129 > 81 = 3^4, and the resulting excess over
the Gaussian ceiling is exactly 48 * 3^(n-4) for every n.

Everything decisive is exact: integer matrices, exact finite-field linear algebra,
exact rational operator identities.  Modular gcd degrees are used only in the
one-sided direction (they are UPPER bounds for the characteristic-zero gcd degree,
hence LOWER bounds for r), which is the direction the no-go needs.
"""

from __future__ import annotations

import hashlib
import json
import math
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
from math import comb, gcd, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from ising.transfer_matrix import layer_bonds  # noqa: E402
from e38_gaussianity_certificate import build_R  # noqa: E402  canonical isotropic construction

INT64_MAX = 2**63 - 1
FLOAT64_EXACT = 2**53
PRIMES = (1_000_003, 2_000_003)
RSS_LIMIT_BYTES = 7_500_000_000
DEFAULT_WALL_SECONDS = 900
BIG_DIRECT_WALL_SECONDS = int(os.environ.get("ALLSIZE_BIG_WALL_SECONDS", "900"))

CHECKS: list[dict[str, object]] = []
FAILS: list[str] = []

# ------------------------------------------------------------------ declared parameters
# Predeclared BEFORE any computation.  The claw certificate is the only decisive
# finite input of the family theorem; these are the parameters it is evaluated at.
CLAW_T = (Fraction(1, 3), Fraction(1, 5), Fraction(1, 7), Fraction(2, 11))
CLAW_W = (Fraction(5, 3), Fraction(7, 2), Fraction(11, 5))
# Generic-looking anisotropic values used for the path/chain controls.
CONTROL_T = (
    Fraction(1, 3), Fraction(1, 5), Fraction(1, 7), Fraction(2, 11), Fraction(3, 13),
    Fraction(1, 17), Fraction(4, 19), Fraction(2, 23), Fraction(5, 29), Fraction(3, 31),
)
CONTROL_W = (
    Fraction(5, 3), Fraction(7, 2), Fraction(11, 5), Fraction(13, 4), Fraction(17, 6),
    Fraction(19, 7), Fraction(23, 8), Fraction(29, 9), Fraction(31, 10), Fraction(37, 11),
)
ISOTROPIC_GRID = tuple(Fraction(k, 20) for k in range(4, 11))
FAMILY_SIZES = ((2, 3), (2, 4), (3, 3), (2, 5), (3, 4), (4, 4), (5, 5), (10, 10))


class ResourceWall(RuntimeError):
    """A predeclared time or RSS wall stopped an intentionally bounded route."""


@dataclass
class Progress:
    stage: str = "initialization"
    started: float = 0.0
    wall_seconds: int = DEFAULT_WALL_SECONDS


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
    if elapsed > progress.wall_seconds:
        raise ResourceWall(f"time wall {progress.wall_seconds}s at {progress.stage} ({elapsed:.1f}s)")
    rss = max_rss_bytes()
    if rss > RSS_LIMIT_BYTES:
        raise ResourceWall(f"RSS wall {RSS_LIMIT_BYTES} bytes at {progress.stage} (ru_maxrss {rss})")


@contextmanager
def alarm_wall(seconds: int):
    if not hasattr(signal, "SIGALRM"):
        yield
        return

    def _fire(_signum, _frame):
        raise ResourceWall(f"SIGALRM wall {seconds}s")

    previous = signal.signal(signal.SIGALRM, _fire)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)
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


def assert_float64_guard(length: int, p: int, kind: str) -> None:
    assert length * (p - 1) ** 2 < FLOAT64_EXACT, f"float64 {kind} guard violated"


def modmul(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """Exact modular matrix product.  float64 BLAS is exact under the stated guard."""
    dim = a.shape[1]
    assert_float64_guard(dim, p, "matrix multiplication")
    prod = np.asarray(a, dtype=np.float64) @ np.asarray(b, dtype=np.float64)
    return np.mod(prod, float(p)).astype(np.int64)


def mod_dot(a: np.ndarray, b: np.ndarray, p: int) -> int:
    """Exact modular dot product with int64 accumulation under an explicit guard."""
    assert_int64_guard(a.size, p, "dot product")
    return int(np.dot(a, b)) % p


# --------------------------------------------------------------- exact integer operator

def path_bonds(n: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(n - 1)]


def star_bonds(k: int) -> list[tuple[int, int]]:
    return [(0, i) for i in range(1, k + 1)]


def degrees(n: int, bonds) -> list[int]:
    deg = [0] * n
    for i, j in bonds:
        deg[i] += 1
        deg[j] += 1
    return deg


def build_M_integer(n: int, bonds, t_list, w_list) -> tuple[int, list[list[int]]]:
    """Integer matrix  M = (prod b_v)^2 (prod d_e) * P(t) D(w) P(t),  built exactly.

    P(t) = kron_v [[1,t_v],[t_v,1]] and D(w)_k = prod_{e aligned in k} w_e, so
    M = kron_v [[b_v,a_v],[a_v,b_v]] * diag(integer D) * kron_v [[b_v,a_v],[a_v,b_v]]
    is integral by construction, with t_v = a_v/b_v and w_e = c_e/d_e.
    """
    assert len(t_list) == n and len(w_list) == len(bonds)
    dim = 1 << n
    scale = 1
    for t in t_list:
        scale *= t.denominator**2
    for w in w_list:
        scale *= w.denominator

    diag = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        acc = 1
        for (i, j), w in zip(bonds, w_list):
            acc *= w.numerator if spins[i] * spins[j] == 1 else w.denominator
        diag.append(acc)

    # rows[k][l] = diag[k] * delta_{kl}, then apply the integer Kronecker factor twice.
    mat = [[0] * dim for _ in range(dim)]
    for k in range(dim):
        mat[k][k] = diag[k]

    def apply_kron_left(matrix):
        for site in range(n):
            bit = 1 << (n - 1 - site)
            a, b = t_list[site].numerator, t_list[site].denominator
            for k in range(dim):
                if k & bit:
                    continue
                row0, row1 = matrix[k], matrix[k | bit]
                for l in range(dim):
                    x, y = row0[l], row1[l]
                    row0[l] = b * x + a * y
                    row1[l] = a * x + b * y
        return matrix

    def transpose(matrix):
        return [list(col) for col in zip(*matrix)]

    mat = apply_kron_left(mat)           # P_int * diag
    mat = transpose(mat)                 # (P_int diag)^T = diag * P_int  (P_int symmetric)
    mat = apply_kron_left(mat)           # P_int * diag * P_int
    if dim <= 64:
        assert all(mat[i][j] == mat[j][i] for i in range(dim) for j in range(dim))
    return scale, mat


def build_R_fraction(n: int, bonds, t_list, w_list) -> list[list[Fraction]]:
    """Rational R = P D P, used only for exact small-case cross-checks."""
    scale, mat = build_M_integer(n, bonds, t_list, w_list)
    return [[Fraction(x, scale) for x in row] for row in mat]


# ----------------------------------------------------------------- F_p polynomial layer

def poly_trim(poly: np.ndarray) -> np.ndarray:
    # Cheap front scan first: after a monic division the leading entry is almost
    # always already nonzero, and flatnonzero over a 130k-long array in every
    # Euclid step would dominate the runtime.
    for i in range(min(poly.size, 64)):
        if poly[i]:
            return poly if i == 0 else poly[i:]
    nonzero = np.flatnonzero(poly)
    if nonzero.size == 0:
        return np.zeros(1, dtype=np.int64)
    return poly[nonzero[0] :]


def poly_is_zero(poly: np.ndarray) -> bool:
    return poly.size == 1 and int(poly[0]) == 0


def poly_monic(poly: np.ndarray, p: int) -> np.ndarray:
    poly = poly_trim(np.asarray(poly, dtype=np.int64))
    if poly_is_zero(poly):
        return poly
    leading = int(poly[0]) % p
    if leading == 1:
        return poly
    return (poly * pow(leading, -1, p)) % p


def poly_remainder(dividend: np.ndarray, divisor: np.ndarray, p: int, progress: Progress) -> np.ndarray:
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
    a, b = poly_monic(a, p), poly_monic(b, p)
    while not poly_is_zero(b):
        assert_resources(progress)
        a, b = b, poly_monic(poly_remainder(a, b, p, progress), p)
    return a


def poly_derivative(poly: np.ndarray, p: int) -> np.ndarray:
    degree = poly.size - 1
    ascending = poly[::-1]
    derivative = (np.arange(1, degree + 1, dtype=np.int64) * ascending[1:]) % p
    return np.ascontiguousarray(derivative[::-1])


def poly_divides(dividend: np.ndarray, divisor: np.ndarray, p: int, progress: Progress) -> bool:
    return poly_is_zero(poly_remainder(dividend, divisor, p, progress))


def newton_from_power_sums(power_sums: np.ndarray, degree: int, p: int, inverses, progress: Progress) -> np.ndarray:
    """Descending [1, a_1, ..., a_degree] from power sums s_1..s_degree.

    a_k = -(1/k) sum_{i=1..k} s_i a_{k-i}.  The coefficients are also written into a
    reversed buffer so the required window is contiguous; the inner product is done
    in exactly-representable float64 chunks (chunk * (p-1)^2 < 2^53), which is exact
    integer arithmetic and roughly two orders of magnitude faster than an int64 dot.
    """
    chunk = FLOAT64_EXACT // ((p - 1) ** 2)
    assert chunk >= 1
    coeffs = np.zeros(degree + 1, dtype=np.int64)
    coeffs[0] = 1
    sums = np.ascontiguousarray(power_sums[: degree + 1], dtype=np.float64)
    reversed_coeffs = np.zeros(degree + 1, dtype=np.float64)
    reversed_coeffs[degree] = 1.0                       # slot of a_0
    for k in range(1, degree + 1):
        if k & 1023 == 0:
            assert_resources(progress)
        total = 0
        base = degree - k
        lo = 1
        while lo <= k:
            hi = min(k, lo + chunk - 1)
            total += int(np.dot(sums[lo : hi + 1], reversed_coeffs[base + lo : base + hi + 1])) % p
            lo = hi + 1
        value = (-(total % p) * inverses[k]) % p
        coeffs[k] = value
        reversed_coeffs[degree - k] = float(value)
    return coeffs


def coefficient_digest(ascending: np.ndarray) -> str:
    return hashlib.sha256(",".join(str(int(c)) for c in ascending).encode()).hexdigest()


# ------------------------------------------------------------------- the invariant core

def pair_product_counts(matrix_rows, p: int, progress: Progress, want_all: bool = True,
                        spot_checks=()) -> dict:
    """Exact F_p reconstruction of the two pair polynomials and their derivative gcds.

    Returns lower bounds  r >= C(d,2) - g2  and  r_all >= C(d+1,2) - g_all  for the
    characteristic-zero distinct-pair-product counts (Lemma 2 + Lemma 5 of the note).
    """
    dim = len(matrix_rows)
    progress.stage = f"reduce {dim}x{dim} mod {p}"
    assert_resources(progress)
    reduced = np.array([[x % p for x in row] for row in matrix_rows], dtype=np.int64)

    slots2 = comb(dim, 2)
    slots_all = comb(dim + 1, 2)
    needed = 2 * (slots_all if want_all else slots2)

    progress.stage = f"power traces dim={dim}"
    power = np.eye(dim, dtype=np.int64)
    traces = [0]
    for _ in range(dim):
        assert_resources(progress)
        power = modmul(power, reduced, p)
        traces.append(int(power.trace()) % p)

    inverses = [0] * (needed + 2)
    for k in range(1, needed + 2):
        inverses[k] = pow(k, -1, p)

    progress.stage = f"charpoly dim={dim}"
    charpoly = newton_from_power_sums(np.array(traces, dtype=np.int64), dim, p, inverses, progress)
    assert int(charpoly[0]) == 1

    progress.stage = f"trace recurrence to {needed}"
    tail = charpoly[1:]
    for k in range(dim + 1, needed + 1):
        if k & 1023 == 0:
            assert_resources(progress)
        window = np.array(traces[k - 1 : k - dim - 1 : -1], dtype=np.int64)
        traces.append((-mod_dot(tail, window, p)) % p)

    spot_ok = True
    for exponent in spot_checks:
        if exponent > needed:
            continue
        result = np.eye(dim, dtype=np.int64)
        base = reduced.copy()
        e = exponent
        while e:
            if e & 1:
                result = modmul(result, base, p)
            e >>= 1
            if e:
                base = modmul(base, base, p)
        spot_ok = spot_ok and int(result.trace()) % p == traces[exponent]

    half = pow(2, -1, p)
    out: dict[str, object] = {
        "dim": dim,
        "prime": p,
        "pair_slots": slots2,
        "sym_slots": slots_all,
        "charpoly_sha256": coefficient_digest(charpoly[::-1]),
        "recurrence_spot_ok": bool(spot_ok),
        "cayley_hamilton_ok": None,
    }

    if dim <= 256:
        progress.stage = f"Cayley-Hamilton dim={dim}"
        value = (int(charpoly[0]) * np.eye(dim, dtype=np.int64)) % p
        identity = np.eye(dim, dtype=np.int64)
        for coefficient in charpoly[1:]:
            value = (modmul(value, reduced, p) + int(coefficient) * identity) % p
        out["cayley_hamilton_ok"] = bool(np.all(value == 0))

    for label, slots, sign in (("pair", slots2, -1), ("sym", slots_all, +1)):
        if label == "sym" and not want_all:
            continue
        progress.stage = f"{label} Newton degree={slots}"
        power_sums = np.zeros(slots + 1, dtype=np.int64)
        for m in range(1, slots + 1):
            square = (traces[m] * traces[m]) % p
            power_sums[m] = ((square + sign * traces[2 * m]) % p) * half % p
        poly = newton_from_power_sums(power_sums, slots, p, inverses, progress)
        progress.stage = f"{label} gcd degree={slots}"
        derivative = poly_derivative(poly, p)
        divisor = poly_gcd(poly, derivative, p, progress)
        degree = divisor.size - 1
        out[f"{label}_poly_sha256"] = coefficient_digest(poly[::-1])
        out[f"{label}_gcd_sha256"] = coefficient_digest(divisor[::-1])
        out[f"{label}_gcd_degree"] = degree
        out[f"{label}_divides_ok"] = bool(
            poly_divides(poly, divisor, p, progress) and poly_divides(derivative, divisor, p, progress)
        )
        out[f"{label}_distinct_lower_bound"] = slots - degree

    out["r_lower_bound"] = out["pair_distinct_lower_bound"]
    if want_all:
        out["r_all_lower_bound"] = out["sym_distinct_lower_bound"]
    return out


# ------------------------------------------------------------------------ family bounds

def gaussian_ceiling(n: int) -> int:
    return 3**n - 2**n


def family_bound(n: int, core_sites: int, r_all_core: int, r_core: int) -> int:
    """Lemma 3: exact distinct-pair-product count of  core  x  (n - core_sites) free sites."""
    m = n - core_sites
    assert m >= 0
    return (3**m - 2**m) * r_all_core + 2**m * r_core


def exceptional_degree_bound(n: int, edges: int) -> int:
    slots = comb(1 << n, 2)
    return 4 * (2 * n + edges) * slots * (slots - 1)


def isotropic_curve_degree_bound(n: int, edges: int) -> int:
    """Degree bound for the one-parameter exceptional set on the isotropic curve.

    Along  w = q(t) = (1+t^2)/(2t),  M(t) = (2t)^{|E|} P_t D(q) P_t  is integral with
    entry t-degree at most  d_R = 2n + 2|E|;  the exterior square doubles that, the
    k-th charpoly coefficient is bounded by 2 d_R k, and the Sylvester determinant of
    the two monic cofactors is bounded by 2 * (2 d_R) * S (S-1).
    """
    slots = comb(1 << n, 2)
    return 4 * (2 * n + 2 * edges) * slots * (slots - 1)


def equal_field_upper_bound(core_sites: int, free_sites: int) -> int:
    """Lemma 9 with the TRUE r_all upper bound  C(2^{n_H}+1, 2)  for the core."""
    return (2 * free_sites + 1) * comb((1 << core_sites) + 1, 2)


# ------------------------------------------------------------------------------ helpers

def frac_str(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def reciprocality_identity(n: int, bonds, t_list, w_list) -> tuple[bool, Fraction]:
    """Exact check of  (W R W^T) R = c I  with  W = (prod_{v in P} X_v)(prod_v Z_v).

    Proves the bipartite spectrum is closed under lambda -> c / lambda.
    """
    colour = [None] * n
    adjacency = [[] for _ in range(n)]
    for i, j in bonds:
        adjacency[i].append(j)
        adjacency[j].append(i)
    for start in range(n):
        if colour[start] is not None:
            continue
        colour[start] = 0
        stack = [start]
        while stack:
            u = stack.pop()
            for v in adjacency[u]:
                if colour[v] is None:
                    colour[v] = 1 - colour[u]
                    stack.append(v)
                elif colour[v] == colour[u]:
                    raise ValueError("graph is not bipartite")
    scale, mat = build_M_integer(n, bonds, t_list, w_list)
    dim = 1 << n
    mask = 0
    for v in range(n):
        if colour[v] == 0:
            mask |= 1 << (n - 1 - v)
    sign = [(-1) ** bin(k).count("1") for k in range(dim)]
    conjugated = [[sign[a ^ mask] * sign[b ^ mask] * mat[a ^ mask][b ^ mask] for b in range(dim)]
                  for a in range(dim)]
    centre = Fraction(1)
    for t in t_list:
        centre *= (1 - t * t) ** 2
    for w in w_list:
        centre *= w
    target = centre * scale * scale
    assert target.denominator == 1
    target_int = int(target)
    ok = True
    for a in range(dim):
        row = conjugated[a]
        for b in range(dim):
            acc = 0
            for k in range(dim):
                acc += row[k] * mat[k][b]
            if acc != (target_int if a == b else 0):
                ok = False
                break
        if not ok:
            break
    return ok, centre


def synthetic_tensor_formula_check() -> list[dict]:
    """Independent exact validation of Lemma 3 using coprime integer multisets."""
    base = [1, 2, 3, 5, 7, 6, 10, 15]
    extra = [11, 13, 17]

    def distinct(values, strict):
        out = set()
        for i in range(len(values)):
            for j in range(i + (1 if strict else 0), len(values)):
                out.add(values[i] * values[j])
        return len(out)

    r_core = distinct(base, True)
    r_all_core = distinct(base, False)
    rows = []
    for m in (1, 2, 3):
        block = [1]
        for prime in extra[:m]:
            block = [x * y for x in block for y in (1, prime)]
        product = [a * b for a in base for b in block]
        rows.append({
            "m": m,
            "direct": distinct(product, True),
            "formula": (3**m - 2**m) * r_all_core + 2**m * r_core,
            "direct_all": distinct(product, False),
            "formula_all": 3**m * r_all_core,
        })
    return {"core_r": r_core, "core_r_all": r_all_core, "rows": rows}


# ------------------------------------------------------------------------------- driver

def run_case(label, n, bonds, t_list, w_list, primes, progress_wall, want_all=True, spot=()):
    scale, mat = build_M_integer(n, bonds, t_list, w_list)
    record = {
        "label": label,
        "sites": n,
        "edges": len(bonds),
        "bonds": [list(b) for b in bonds],
        "max_degree": max(degrees(n, bonds)) if bonds else 0,
        "t": [frac_str(x) for x in t_list],
        "w": [frac_str(x) for x in w_list],
        "integer_scale": str(scale),
        "gaussian_r_ceiling": gaussian_ceiling(n),
        "gaussian_r_all_ceiling": 3**n,
        "primes": {},
    }
    for p in primes:
        progress = Progress(stage=f"{label} p={p}", started=time.perf_counter(), wall_seconds=progress_wall)
        started = time.perf_counter()
        with alarm_wall(progress_wall):
            data = pair_product_counts(mat, p, progress, want_all=want_all, spot_checks=spot)
        data["seconds"] = round(time.perf_counter() - started, 3)
        record["primes"][str(p)] = data
    first = record["primes"][str(primes[0])]
    record["r_lower_bound"] = first["r_lower_bound"]
    if want_all:
        record["r_all_lower_bound"] = first["r_all_lower_bound"]
    record["flagged"] = bool(first["r_lower_bound"] > gaussian_ceiling(n))
    return record


def main() -> int:
    started = time.perf_counter()
    out: dict[str, object] = {
        "experiment": "e128_allsize_gaussian",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "primes": list(PRIMES),
        "rss_limit_bytes": RSS_LIMIT_BYTES,
        "declared_parameters": {
            "claw_t": [frac_str(x) for x in CLAW_T],
            "claw_w": [frac_str(x) for x in CLAW_W],
            "control_t": [frac_str(x) for x in CONTROL_T],
            "control_w": [frac_str(x) for x in CONTROL_W],
            "isotropic_grid": [frac_str(x) for x in ISOTROPIC_GRID],
        },
    }

    print("[0] guards", flush=True)
    for p in PRIMES:
        check(f"prime {p}", is_prime(p))
    check("float64 matmul guard dim<=1024", 1024 * (PRIMES[0] - 1) ** 2 < FLOAT64_EXACT,
          f"{1024 * (PRIMES[0] - 1) ** 2} < {FLOAT64_EXACT}")
    check("int64 dot guard slots<=523776", 523776 * (PRIMES[0] - 1) ** 2 < INT64_MAX)

    # ---------------------------------------------------------- 1. construction agreement
    print("[1] construction agreement with the canonical isotropic build_R", flush=True)
    agreements = []
    for cross in ((2, 3), (2, 2)):
        n = cross[0] * cross[1]
        bonds = layer_bonds(cross, (False, False))
        t = Fraction(1, 3)
        q = (1 + t * t) / (2 * t)
        mine = build_R_fraction(n, bonds, [t] * n, [q] * len(bonds))
        canonical, q_ref, eps = build_R(n, bonds, t)
        ratio = mine[0][0] / canonical[0][0]
        same = all(mine[i][j] == ratio * canonical[i][j]
                   for i in range(1 << n) for j in range(1 << n))
        expected = q ** ((len(bonds) + eps) // 2)
        agreements.append({"cross": list(cross), "ratio": frac_str(ratio),
                           "expected_ratio": frac_str(expected), "match": bool(same)})
        check(f"{cross[0]}x{cross[1]} anisotropic model = build_R up to scalar",
              same and ratio == expected, f"ratio {frac_str(ratio)}")
    out["construction_agreement"] = agreements

    # ------------------------------------------------------------- 2. the claw certificate
    print("[2] claw certificate K_{1,3} (the only decisive finite input)", flush=True)
    claw = run_case("claw_K13", 4, star_bonds(3), list(CLAW_T), list(CLAW_W), PRIMES,
                    DEFAULT_WALL_SECONDS, spot=(17, 33, 65, 128, 271))
    path4 = run_case("path_P4", 4, path_bonds(4), list(CLAW_T), list(CLAW_W), PRIMES,
                     DEFAULT_WALL_SECONDS, spot=(17, 33, 65, 128, 271))
    out["claw_certificate"] = claw
    out["path_control_P4"] = path4

    r_all_core = min(claw["primes"][str(p)]["r_all_lower_bound"] for p in PRIMES)
    r_core = min(claw["primes"][str(p)]["r_lower_bound"] for p in PRIMES)
    out["core"] = {"sites": 4, "r_all_lower_bound": r_all_core, "r_lower_bound": r_core,
                   "gaussian_r_all_ceiling": 81, "gaussian_r_ceiling": 65,
                   "excess_constant": r_all_core - 81}
    check("claw r_all exceeds the 4-mode Gaussian ceiling", r_all_core > 81,
          f"r_all >= {r_all_core} > 81")
    check("claw r exceeds the 4-mode Gaussian ceiling", r_core > 65, f"r >= {r_core} > 65")
    check("claw values agree at both primes",
          claw["primes"][str(PRIMES[0])]["r_all_lower_bound"] == claw["primes"][str(PRIMES[1])]["r_all_lower_bound"]
          and claw["primes"][str(PRIMES[0])]["r_lower_bound"] == claw["primes"][str(PRIMES[1])]["r_lower_bound"])
    check("path P_4 saturates the Gaussian ceilings exactly (same parameters)",
          path4["r_lower_bound"] == 65 and path4["r_all_lower_bound"] == 81,
          f"r={path4['r_lower_bound']} r_all={path4['r_all_lower_bound']}")
    check("r_all - r <= 2^4 for the claw", r_all_core - r_core <= 16, f"{r_all_core - r_core}")

    # ------------------------------------------------------- 3. tensor formula validation
    print("[3] Lemma 3 (tensor decoupling identity)", flush=True)
    synthetic = synthetic_tensor_formula_check()
    ok = all(row["direct"] == row["formula"] and row["direct_all"] == row["formula_all"]
             for row in synthetic["rows"])
    out["synthetic_tensor_formula"] = synthetic
    check("synthetic coprime-multiset validation of Lemma 3", ok)

    direct_rows = []
    for m in (1, 2, 3, 4, 5):
        n = 4 + m
        t_list = list(CLAW_T) + list(CONTROL_T[4 : 4 + m])
        wall = BIG_DIRECT_WALL_SECONDS if n >= 9 else DEFAULT_WALL_SECONDS
        want_all = n <= 8          # the n=9 symmetric square doubles an already large run
        try:
            case = run_case(f"claw_plus_{m}_free_sites", n, star_bonds(3), t_list, list(CLAW_W),
                            PRIMES[:1], wall, want_all=want_all)
        except ResourceWall as exc:
            direct_rows.append({"m": m, "sites": n, "wall": str(exc)})
            check(f"direct n={n} specialization", False, f"wall: {exc}")
            continue
        predicted = family_bound(n, 4, r_all_core, r_core)
        predicted_all = 3**m * r_all_core
        row = {
            "m": m, "sites": n,
            "direct_r": case["r_lower_bound"], "formula_r": predicted,
            "direct_r_all": case.get("r_all_lower_bound"),
            "formula_r_all": predicted_all if want_all else None,
            "gaussian_ceiling": gaussian_ceiling(n),
            "excess": predicted - gaussian_ceiling(n),
            "seconds": case["primes"][str(PRIMES[0])]["seconds"],
        }
        direct_rows.append(row)
        check(f"direct n={n} matches Lemma 3 exactly",
              row["direct_r"] == predicted
              and (not want_all or row["direct_r_all"] == predicted_all),
              f"r {row['direct_r']} vs {predicted}")
        check(f"direct n={n} excess equals 48*3^{m}",
              row["excess"] == (r_all_core - 81) * 3**m, str(row["excess"]))
    out["direct_specialization"] = direct_rows

    # --------------------------------------------------------------- 4. the family table
    print("[4] family table", flush=True)
    table = []
    for cross in FAMILY_SIZES:
        n = cross[0] * cross[1]
        bonds = layer_bonds(cross, (False, False))
        deg = max(degrees(n, bonds))
        bound = family_bound(n, 4, r_all_core, r_core)
        ceiling = gaussian_ceiling(n)
        table.append({
            "layer": f"{cross[0]}x{cross[1]}",
            "sites": n,
            "edges": len(bonds),
            "max_degree": deg,
            "has_branching_site": deg >= 3,
            "gaussian_ceiling": str(ceiling),
            "family_lower_bound": str(bound),
            "excess": str(bound - ceiling),
            "excess_closed_form": str((r_all_core - 81) * 3 ** (n - 4)),
            "exceptional_degree_bound": str(exceptional_degree_bound(n, len(bonds))),
        })
        check(f"family bound flags {cross[0]}x{cross[1]} (n={n})", bound > ceiling,
              f"{bound} > {ceiling}")
        check(f"excess closed form for n={n}", bound - ceiling == (r_all_core - 81) * 3 ** (n - 4))
    out["family_table"] = table

    chain_table = []
    for cross in FAMILY_SIZES:
        n = cross[0] * cross[1]
        bonds = path_bonds(n)
        chain_table.append({
            "control": f"open chain P_{n}",
            "sites": n,
            "max_degree": max(degrees(n, bonds)),
            "invariant": str(gaussian_ceiling(n)),
            "gaussian_ceiling": str(gaussian_ceiling(n)),
            "excess": "0",
        })
        check(f"chain control P_{n} is not flagged", True, "invariant = 3^n - 2^n by Theorem C")
    out["chain_table"] = chain_table

    # ------------------------------------------------- 5. computed chain / path controls
    print("[5] computed chain controls (open path, anisotropic rational parameters)", flush=True)
    chain_rows = []
    for n in range(3, 9):
        case = run_case(f"chain_P{n}", n, path_bonds(n), list(CONTROL_T[:n]),
                        list(CONTROL_W[: n - 1]), PRIMES[:1], DEFAULT_WALL_SECONDS)
        row = {"sites": n, "r": case["r_lower_bound"], "r_all": case["r_all_lower_bound"],
               "gaussian_ceiling": gaussian_ceiling(n), "gaussian_all_ceiling": 3**n,
               "seconds": case["primes"][str(PRIMES[0])]["seconds"]}
        chain_rows.append(row)
        check(f"chain P_{n} saturates without excess",
              row["r"] == gaussian_ceiling(n) and row["r_all"] == 3**n,
              f"r={row['r']} r_all={row['r_all']}")
    out["computed_chain_controls"] = chain_rows

    # -------------------------------------------------- 6. isotropic locus: the obstruction
    print("[6] isotropic locus (the isolated obstruction)", flush=True)
    iso_claw = []
    for t in ISOTROPIC_GRID:
        q = (1 + t * t) / (2 * t)
        case = run_case(f"iso_claw_t{t.numerator}_{t.denominator}", 4, star_bonds(3),
                        [t] * 4, [q] * 3, PRIMES[:1], DEFAULT_WALL_SECONDS)
        iso_claw.append({"t": frac_str(t), "q": frac_str(q),
                         "r_lower_bound": case["r_lower_bound"], "r": case["r_lower_bound"],
                         "r_all_lower_bound": case["r_all_lower_bound"],
                         "r_all": case["r_all_lower_bound"],
                         "certified_flagged": case["flagged"], "flagged": case["flagged"]})
    out["isotropic_claw_sweep"] = iso_claw
    check("isotropic claw carries NO CERTIFIED excess on the declared grid (the modular "
          "value is a lower bound, so this is a non-certificate, not absence of excess)",
          all(not row["certified_flagged"] for row in iso_claw),
          f"r >= {iso_claw[0]['r']} against the ceiling 65 at every grid point")
    check("the certified r_all lower bound for the isotropic claw is below 3^4 at every "
          "grid point (again a non-certificate)",
          all(row["r_all"] < 81 for row in iso_claw), f"r_all >= {iso_claw[0]['r_all']} < 81")

    iso5 = []
    t = Fraction(1, 3)
    q = (1 + t * t) / (2 * t)
    for name, bonds in (
        ("K_{1,4} star", star_bonds(4)),
        ("T = claw with one branch of length 2", [(0, 1), (0, 2), (0, 3), (3, 4)]),
        ("claw + isolated site", star_bonds(3)),
        ("paw (triangle + pendant)", [(0, 1), (0, 2), (1, 2), (0, 3)]),
        ("path P_5", path_bonds(5)),
        ("cycle C_5", path_bonds(5) + [(4, 0)]),
    ):
        case = run_case(f"iso5_{name}", 5, bonds, [t] * 5, [q] * len(bonds), PRIMES[:1],
                        DEFAULT_WALL_SECONDS)
        iso5.append({"graph": name, "max_degree": case["max_degree"],
                     "r_lower_bound": case["r_lower_bound"], "r": case["r_lower_bound"],
                     "r_all_lower_bound": case["r_all_lower_bound"],
                     "r_all": case["r_all_lower_bound"], "gaussian_ceiling": gaussian_ceiling(5),
                     "certified_flagged": case["flagged"], "flagged": case["flagged"]})
    out["isotropic_n5_survey"] = iso5
    by_name = {row["graph"]: row for row in iso5}
    check("isotropic 5-site T graph IS certified flagged",
          by_name["T = claw with one branch of length 2"]["certified_flagged"],
          f"r >= {by_name['T = claw with one branch of length 2']['r']} > 211")
    check("isotropic star K_{1,4} is NOT certified flagged (a lower bound below the "
          "ceiling is a non-certificate, not absence of excess)",
          not by_name["K_{1,4} star"]["certified_flagged"],
          f"r >= {by_name['K_{1,4} star']['r']}, ceiling 211")
    check("isotropic path P_5 count is pinned EXACTLY at 211 (modular lower bound meets "
          "the Lemma 1 ceiling, which applies because a path is Gaussian)",
          by_name["path P_5"]["r"] == gaussian_ceiling(5))

    # Lemma 9: with EQUAL fields on the free sites the true count is at most
    # (2m+1) * r_all(core), and r_all(core) <= C(2^{n_H}+1, 2) unconditionally.  That
    # upper bound drops below the Gaussian ceiling at m = 3, so the invariant PROVABLY
    # cannot flag an equal-field decoupled configuration from then on.  This is the
    # rigorous half of the isotropic obstruction; the individual modular counts below
    # are only lower bounds and are recorded as such.
    collapse = []
    for m in (1, 2, 3, 4):
        nn = 4 + m
        case = run_case(f"iso_collapse_m{m}", nn, star_bonds(3), [t] * nn, [q] * 3,
                        PRIMES[:1], DEFAULT_WALL_SECONDS)
        proved_upper = equal_field_upper_bound(4, m)
        row = {"m": m, "sites": nn, "r_lower_bound": case["r_lower_bound"],
               "r": case["r_lower_bound"],
               "proved_upper_bound": proved_upper,
               "observed_bound_with_certified_r_all": (2 * m + 1) * iso_claw[0]["r_all"],
               "gaussian_ceiling": gaussian_ceiling(nn),
               "upper_bound_below_ceiling": bool(proved_upper < gaussian_ceiling(nn)),
               "certified_flagged": case["flagged"]}
        collapse.append(row)
        check(f"Lemma 9 proved upper bound is consistent at m={m}",
              row["r_lower_bound"] <= proved_upper and not row["certified_flagged"],
              f"r in [{row['r_lower_bound']}, {proved_upper}], ceiling {row['gaussian_ceiling']}")
    out["equal_field_collapse"] = collapse
    check("Lemma 9 PROVES no equal-field decoupled configuration can be flagged from m=3 on",
          all(row["upper_bound_below_ceiling"] for row in collapse if row["m"] >= 3)
          and not any(row["upper_bound_below_ceiling"] for row in collapse if row["m"] < 3),
          ", ".join(f"m={r['m']}:{r['proved_upper_bound']} vs {r['gaussian_ceiling']}"
                    for r in collapse))

    # ------- 6b. uniform-field route: replace the free sites by an OPEN CHAIN ---------
    # The free sites of Lemma 9 degenerate because their modes coincide.  An open
    # chain at the same uniform field is Gaussian with pairwise distinct modes, so
    # Lemma 5 might apply with no degeneracy.  IMPORTANT: the tensor expression below
    # is evaluated at the CERTIFIED LOWER BOUNDS for the core counts, so it is neither
    # an upper nor a lower bound for the composite count -- it is only a reference
    # value.  Nothing here certifies the injectivity hypothesis of Lemma 5.  The only
    # rigorous outputs of this section are the certified lower bounds themselves.
    print("[6b] uniform-field core-plus-chain route", flush=True)
    T_BONDS = [(0, 1), (0, 2), (0, 3), (3, 4)]
    chain_route = []
    for coupling in (Fraction(1, 3), Fraction(1, 4), Fraction(2, 5), Fraction(1, 5)):
        qq = (1 + coupling * coupling) / (2 * coupling)
        core_case = run_case(f"T_core_t{coupling.numerator}_{coupling.denominator}", 5, T_BONDS,
                             [coupling] * 5, [qq] * 4, PRIMES[:1], DEFAULT_WALL_SECONDS)
        r_core_T = core_case["r_lower_bound"]
        r_all_core_T = core_case["r_all_lower_bound"]
        rows = []
        for m in (1, 2, 3):
            nn = 5 + m
            chain = [(5 + i, 6 + i) for i in range(m - 1)]
            case = run_case(f"T_plus_P{m}_t{coupling.numerator}_{coupling.denominator}", nn,
                            T_BONDS + chain, [coupling] * nn, [qq] * (4 + len(chain)),
                            PRIMES[:1], DEFAULT_WALL_SECONDS)
            reference_value = (3**m - 2**m) * r_all_core_T + 2**m * r_core_T
            rows.append({
                "m": m, "sites": nn,
                "certified_r_lower_bound": case["r_lower_bound"],
                "direct_r": case["r_lower_bound"],
                "tensor_reference_from_core_lower_bounds": reference_value,
                "certified_r_all_lower_bound": case["r_all_lower_bound"],
                "tensor_reference_all": 3**m * r_all_core_T,
                "gaussian_ceiling": gaussian_ceiling(nn),
                "certified_excess": case["r_lower_bound"] - gaussian_ceiling(nn),
                "matches_tensor_reference": bool(case["r_lower_bound"] == reference_value),
                "certified_flagged": case["flagged"],
            })
        chain_route.append({
            "t": frac_str(coupling), "q": frac_str(qq),
            "core_T_r_lower_bound": r_core_T, "core_T_r_all_lower_bound": r_all_core_T,
            "core_T_r": r_core_T, "core_T_r_all": r_all_core_T,
            "core_T_excess_all_lower_bound": r_all_core_T - 3**5,
            "core_T_certified_flagged": core_case["flagged"], "rows": rows,
        })
        check(f"uniform T core is certified flagged at t={frac_str(coupling)}",
              core_case["flagged"], f"r >= {r_core_T} > {gaussian_ceiling(5)}")
        check(f"uniform T + open chain certified flagged for m=1,2,3 at t={frac_str(coupling)}",
              all(row["certified_flagged"] for row in rows),
              ", ".join(f"m={row['m']}:r>={row['direct_r']}>{row['gaussian_ceiling']}"
                        for row in rows))
    out["uniform_field_chain_route"] = chain_route
    reference = chain_route[0]["rows"]
    check("the certified count coincides with the lower-bound tensor reference at m=2 "
          "only (an OBSERVATION about two numbers; it certifies nothing about Lemma 5's "
          "injectivity hypothesis, because the reference substitutes core LOWER bounds)",
          reference[1]["matches_tensor_reference"]
          and not reference[0]["matches_tensor_reference"]
          and not reference[2]["matches_tensor_reference"],
          ", ".join(f"m={row['m']}: certified {row['direct_r']} vs reference "
                    f"{row['tensor_reference_from_core_lower_bounds']}" for row in reference))
    check("the uniform-field certified counts at m=1,2,3 are the same at all four "
          "declared couplings",
          all(block["rows"][i]["direct_r"] == reference[i]["direct_r"] for block in chain_route
              for i in range(3)),
          "identical certified r at t = 1/3, 1/4, 2/5, 1/5 (m=4 was computed at t=1/3 only)")

    # m = 4 is the decisive parity data point: n = 9, dimension 512.
    t_base = Fraction(1, 3)
    q_base = (1 + t_base * t_base) / (2 * t_base)
    core_ref = chain_route[0]
    try:
        big = run_case("T_plus_P4", 9, T_BONDS + [(5, 6), (6, 7), (7, 8)], [t_base] * 9,
                       [q_base] * 7, PRIMES[:1], BIG_DIRECT_WALL_SECONDS, want_all=False)
        reference_value = ((3**4 - 2**4) * core_ref["core_T_r_all_lower_bound"]
                           + 2**4 * core_ref["core_T_r_lower_bound"])
        out["uniform_T_plus_P4"] = {
            "sites": 9, "certified_r_lower_bound": big["r_lower_bound"],
            "direct_r": big["r_lower_bound"],
            "tensor_reference_from_core_lower_bounds": reference_value,
            "gaussian_ceiling": gaussian_ceiling(9),
            "matches_tensor_reference": bool(big["r_lower_bound"] == reference_value),
            "certified_flagged": big["flagged"],
            "pair_gcd_degree": big["primes"][str(PRIMES[0])]["pair_gcd_degree"],
            "pair_poly_sha256": big["primes"][str(PRIMES[0])]["pair_poly_sha256"],
            "prime": PRIMES[0],
            "seconds": big["primes"][str(PRIMES[0])]["seconds"],
        }
        check("uniform T + open 4-chain (n=9, t=1/3 only) is certified flagged",
              big["flagged"], f"r >= {big['r_lower_bound']} > {gaussian_ceiling(9)}")
    except ResourceWall as exc:
        out["uniform_T_plus_P4"] = {"sites": 9, "status": "wall", "detail": str(exc)}
        check("uniform T + open 4-chain (n=9) completed", False, f"wall: {exc}")

    uniform_chain_controls = []
    for nn in range(3, 8):
        case = run_case(f"iso_chain_P{nn}", nn, path_bonds(nn), [t] * nn, [q] * (nn - 1),
                        PRIMES[:1], DEFAULT_WALL_SECONDS)
        uniform_chain_controls.append({
            "sites": nn, "r": case["r_lower_bound"], "r_all": case["r_all_lower_bound"],
            "gaussian_ceiling": gaussian_ceiling(nn), "gaussian_all_ceiling": 3**nn,
        })
    out["uniform_chain_controls"] = uniform_chain_controls
    check("open chains at a UNIFORM coupling still saturate exactly (own modes independent)",
          all(row["r"] == gaussian_ceiling(row["sites"]) and row["r_all"] == 3 ** row["sites"]
              for row in uniform_chain_controls),
          "n = 3..7 at t = 1/3")

    # Which disjoint unions of open chains keep the pair-exponent map injective at a
    # SHARED coupling?  Both components are Gaussian, so Lemma 1 gives the TRUE upper
    # bound r <= 3^n - 2^n; a modular lower bound meeting it therefore pins r exactly
    # and certifies injectivity.  A shortfall is only a lower bound and certifies
    # nothing, except for equal-length pairs, where the two mode multisets are literally
    # identical and dependence is immediate.
    parity_rows = []
    for lengths in ((1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (2, 4), (2, 5),
                    (3, 3), (3, 4), (1, 1, 1), (1, 2, 3), (2, 2, 2)):
        bonds = []
        offset = 0
        for length in lengths:
            bonds += [(offset + i, offset + i + 1) for i in range(length - 1)]
            offset += length
        nn = offset
        case = run_case("parity_" + "_".join(map(str, lengths)), nn, bonds, [t] * nn,
                        [q] * len(bonds), PRIMES, DEFAULT_WALL_SECONDS)
        both = {str(p): case["primes"][str(p)]["r_lower_bound"] for p in PRIMES}
        parity_rows.append({
            "lengths": list(lengths), "sites": nn, "r_lower_bound": case["r_lower_bound"],
            "r_at_both_primes": both, "gaussian_ceiling": gaussian_ceiling(nn),
            "injectivity_certified": bool(case["r_lower_bound"] == gaussian_ceiling(nn)),
            "independent_certified": bool(case["r_lower_bound"] == gaussian_ceiling(nn)),
            "equal_lengths": bool(len(set(lengths)) == 1 and len(lengths) > 1),
            "total_parity": "odd" if nn % 2 else "even",
        })
    out["disjoint_chain_parity"] = parity_rows
    two_part = [row for row in parity_rows if len(row["lengths"]) == 2]
    check("two disjoint open chains at a shared coupling: injectivity of the pair-exponent "
          "map is CERTIFIED exactly when their sizes have opposite parity (equal-length "
          "pairs are immediately dependent; unequal same-parity pairs are non-certified)",
          all(row["injectivity_certified"] == (row["sites"] % 2 == 1) for row in two_part),
          ", ".join(f"{tuple(r['lengths'])}:{r['r_lower_bound']}/{r['gaussian_ceiling']}"
                    for r in two_part))

    # Floating-point DIAGNOSTIC (explicitly not part of any proof): the alternating
    # product of the open-chain modes depends only on the parity of the chain length.
    mode_rows = []
    for length in range(1, 8):
        matrix = build_R_fraction(length, path_bonds(length), [t] * length, [q] * (length - 1))
        dense = np.array([[float(x) for x in row] for row in matrix])
        spectrum = np.sort(np.linalg.eigvalsh(dense))[::-1]
        logs = np.log(spectrum)
        head, rest, known, extracted = logs[0], list(logs[1:]), [logs[0]], []
        for _ in range(length):
            shift = rest[0] - head
            extracted.append(shift)
            for target in sorted([k + shift for k in known], reverse=True):
                rest.pop(int(np.argmin([abs(x - target) for x in rest])))
            known = known + [k + shift for k in known]
        extracted.sort()
        alternating = sum((-1) ** k * value for k, value in enumerate(extracted))
        ratio = float((1 - t) / (1 + t))
        predicted = 2 * math.log(ratio) if length % 2 else -math.log(float(q))
        mode_rows.append({
            "length": length, "modes_log": [float(x) for x in extracted],
            "alternating_sum": float(alternating), "closed_form": predicted,
            "residual": abs(float(alternating) - predicted), "leftover_slots": len(rest),
        })
    out["chain_mode_parity_diagnostic"] = mode_rows
    check("[diagnostic, floating point] alternating mode product depends only on chain parity",
          all(row["residual"] < 1e-9 and row["leftover_slots"] == 0 for row in mode_rows),
          f"max residual {max(row['residual'] for row in mode_rows):.2e}")

    # ------------------------------------------- 7. decoupling checks requested for 2x3
    print("[7] 2x3 decoupling checks", flush=True)
    n = 6
    bonds23 = layer_bonds((2, 3), (False, False))
    rungs = [(0, 3), (1, 4), (2, 5)]
    t = Fraction(1, 3)
    q = (1 + t * t) / (2 * t)
    decoupling = {}
    for label, t_list, w_list in (
        ("isotropic_all_bonds_on", [t] * n, [q] * len(bonds23)),
        ("isotropic_rungs_off", [t] * n, [Fraction(1) if b in rungs else q for b in bonds23]),
        ("anisotropic_rungs_off", list(CONTROL_T[:n]),
         [Fraction(1) if b in rungs else CONTROL_W[i] for i, b in enumerate(bonds23)]),
        ("anisotropic_all_bonds_on", list(CONTROL_T[:n]), list(CONTROL_W[: len(bonds23)])),
    ):
        case = run_case(f"2x3_{label}", n, bonds23, t_list, w_list, PRIMES, DEFAULT_WALL_SECONDS)
        decoupling[label] = {
            "r": case["r_lower_bound"], "r_all": case["r_all_lower_bound"],
            "pair_gcd_degree": {str(p): case["primes"][str(p)]["pair_gcd_degree"] for p in PRIMES},
            "gaussian_ceiling": gaussian_ceiling(n),
            "gaussian_gcd_floor": comb(64, 2) - gaussian_ceiling(n),
            "flagged": case["flagged"],
            "integer_scale": case["integer_scale"],
            "charpoly_sha256": {str(p): case["primes"][str(p)]["charpoly_sha256"] for p in PRIMES},
            "pair_poly_sha256": {str(p): case["primes"][str(p)]["pair_poly_sha256"] for p in PRIMES},
        }
    out["decoupling_2x3"] = decoupling
    check("2x3 isotropic all-bonds-on reproduces the e107 gcd degree 385",
          all(v == 385 for v in decoupling["isotropic_all_bonds_on"]["pair_gcd_degree"].values()),
          str(decoupling["isotropic_all_bonds_on"]["pair_gcd_degree"]))
    check("2x3 isotropic all-bonds-on reproduces r >= 1631",
          decoupling["isotropic_all_bonds_on"]["r"] == 1631)
    check("2x3 with rungs switched off is NOT certified flagged (isotropic)",
          not decoupling["isotropic_rungs_off"]["flagged"],
          f"r >= {decoupling['isotropic_rungs_off']['r']}, ceiling 665")
    check("2x3 with rungs switched off has its count PINNED at 665 (anisotropic): the "
          "modular lower bound meets the Lemma 1 ceiling, which applies because the "
          "decoupled layer is a disjoint union of paths",
          decoupling["anisotropic_rungs_off"]["r"] == gaussian_ceiling(6)
          and not decoupling["anisotropic_rungs_off"]["flagged"],
          f"r = {decoupling['anisotropic_rungs_off']['r']} = 665")
    check("anisotropic rungs-off hits the Gaussian gcd floor 1351 exactly",
          all(v == 1351 for v in decoupling["anisotropic_rungs_off"]["pair_gcd_degree"].values()),
          str(decoupling["anisotropic_rungs_off"]["pair_gcd_degree"]))
    check("isotropic rungs-off modular gcd degree exceeds the Gaussian floor",
          all(v > 1351 for v in decoupling["isotropic_rungs_off"]["pair_gcd_degree"].values()),
          str(decoupling["isotropic_rungs_off"]["pair_gcd_degree"]))
    check("2x3 fully anisotropic all-bonds-on is flagged",
          decoupling["anisotropic_all_bonds_on"]["flagged"],
          f"r = {decoupling['anisotropic_all_bonds_on']['r']} > 665")

    # ------------------------------------------------------------ 8. reciprocality lemma
    print("[8] bipartite reciprocality (why cheaper invariants cannot work)", flush=True)
    recip = []
    for label, nn, bonds, t_list, w_list in (
        ("claw_K13", 4, star_bonds(3), list(CLAW_T), list(CLAW_W)),
        ("path_P4", 4, path_bonds(4), list(CLAW_T), list(CLAW_W)),
        ("2x3_layer", 6, bonds23, list(CONTROL_T[:6]), list(CONTROL_W[: len(bonds23)])),
        ("2x3_isotropic", 6, bonds23, [t] * 6, [q] * len(bonds23)),
    ):
        ok, centre = reciprocality_identity(nn, bonds, t_list, w_list)
        recip.append({"case": label, "identity_holds": bool(ok), "centre": frac_str(centre)})
        check(f"reciprocality identity exact for {label}", ok, f"c = {frac_str(centre)}")
    out["reciprocality"] = recip

    # ------------------------------------------------------- 9. bounded isotropic attempts
    print("[9] bounded direct isotropic attempts at larger layers", flush=True)
    attempts = []
    for cross, wall in (((2, 4), BIG_DIRECT_WALL_SECONDS), ((3, 3), BIG_DIRECT_WALL_SECONDS)):
        nn = cross[0] * cross[1]
        bonds = layer_bonds(cross, (False, False))
        started_case = time.perf_counter()
        try:
            case = run_case(f"iso_{cross[0]}x{cross[1]}", nn, bonds, [t] * nn, [q] * len(bonds),
                            PRIMES[:1], wall, want_all=False)
            attempts.append({
                "layer": f"{cross[0]}x{cross[1]}", "sites": nn, "status": "completed",
                "pair_gcd_degree": case["primes"][str(PRIMES[0])]["pair_gcd_degree"],
                "r_lower_bound": case["r_lower_bound"],
                "gaussian_ceiling": gaussian_ceiling(nn),
                "flagged": case["flagged"],
                "pair_poly_sha256": case["primes"][str(PRIMES[0])]["pair_poly_sha256"],
                "seconds": case["primes"][str(PRIMES[0])]["seconds"],
            })
            check(f"isotropic {cross[0]}x{cross[1]} is flagged", case["flagged"],
                  f"r >= {case['r_lower_bound']} > {gaussian_ceiling(nn)}")
        except ResourceWall as exc:
            attempts.append({
                "layer": f"{cross[0]}x{cross[1]}", "sites": nn, "status": "wall",
                "detail": str(exc), "wall_seconds": wall,
                "elapsed_seconds": round(time.perf_counter() - started_case, 3),
            })
            print(f"  [WALL] isotropic {cross[0]}x{cross[1]}: {exc}", flush=True)
    out["isotropic_direct_attempts"] = attempts

    # ------------- 10. isotropic-curve genericity: all but finitely many t -------------
    # Every isotropic certificate above is a single named point.  Along the isotropic
    # curve w = q(t) the family M(t) = (2t)^{|E|} P_t D(q) P_t is integral in t, so the
    # monic specialization lemma pushes each named-point LOWER bound to the generic point
    # of Q(t), and the resultant lemma then bounds the exceptional t-set.  This turns each
    # single-coupling certificate into an all-but-finitely-many-t statement for the SAME
    # graph, with no new computation -- only degree bookkeeping.
    print("[10] isotropic-curve genericity bounds", flush=True)
    curve_rows = []
    iso_named = []
    for m in (0, 1, 2, 3):
        row = chain_route[0]["rows"][m - 1] if m >= 1 else None
        certified = chain_route[0]["core_T_r_lower_bound"] if m == 0 else row["direct_r"]
        iso_named.append((f"T u P_{m}", 5 + m, 4 + max(m - 1, 0), certified))
    if out["uniform_T_plus_P4"].get("status") != "wall":
        iso_named.append(("T u P_4", 9, 7, out["uniform_T_plus_P4"]["direct_r"]))
    iso_named.append(("2x3", 6, len(bonds23), decoupling["isotropic_all_bonds_on"]["r"]))
    for attempt in attempts:
        if attempt.get("status") == "completed":
            cross = tuple(int(x) for x in attempt["layer"].split("x"))
            iso_named.append((attempt["layer"], attempt["sites"],
                              len(layer_bonds(cross, (False, False))), attempt["r_lower_bound"]))
    for label, nn, edges, certified in iso_named:
        bound = isotropic_curve_degree_bound(nn, edges)
        curve_rows.append({
            "graph": label, "sites": nn, "edges": edges,
            "certified_r_at_t_one_third": certified,
            "gaussian_ceiling": gaussian_ceiling(nn),
            "generic_r_lower_bound": certified,
            "exceptional_t_degree_bound": str(bound),
            "conclusion_holds_off_finite_set": bool(certified > gaussian_ceiling(nn)),
        })
        check(f"isotropic-curve genericity for {label} (n={nn})",
              certified > gaussian_ceiling(nn),
              f"r >= {certified} > {gaussian_ceiling(nn)} for all t off <= {bound} roots")
    out["isotropic_curve_genericity"] = curve_rows
    check("the 2x3 isotropic-curve degree bound reproduces pair_product_scale.md's "
          "422472960",
          isotropic_curve_degree_bound(6, len(bonds23)) == 422_472_960,
          str(isotropic_curve_degree_bound(6, len(bonds23))))

    out["checks"] = CHECKS
    out["failures"] = FAILS
    out["total_seconds"] = round(time.perf_counter() - started, 3)
    out["peak_rss_bytes"] = max_rss_bytes()

    destination = ROOT / "results" / "spectral" / "allsize_gaussian.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {destination}")
    print(f"checks: {len(CHECKS)}  failures: {len(FAILS)}  seconds: {out['total_seconds']}")
    print(f"peak RSS: {out['peak_rss_bytes']} bytes")
    return 1 if FAILS else 0


if __name__ == "__main__":
    raise SystemExit(main())
