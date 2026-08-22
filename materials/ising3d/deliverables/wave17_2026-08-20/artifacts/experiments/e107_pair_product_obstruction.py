"""EXACT multiplicative pair-product obstruction to six-mode spectral Gaussianity.

This is a second, ORDERING-FREE certificate of Theorem S
(`proofs/spectral_gaussianity.md`): the 2x3 3D-Ising layer transfer spectrum at
t = tanh(K*/2) = 1/3 is not a full six-mode fermionic-Gaussian subset-product
spectrum.  It uses no eigenvalue sorting, no inertia windows, no enclosures and
no subset-sum reconstruction -- only the multiset of eigenvalues enters, through
MULTIPLICATIVE COLLISIONS among pair products.

NECESSARY CONDITION (proved in proofs/pair_product_obstruction.md).
A full Gaussian spectrum on n = 6 modes is the multiset

    { a * prod_{i in S} u_i :  S subset of {1..6} }        (64 values, a, u_i arbitrary).

Every pair product over two DISTINCT eigenvalue slots has exponent vector in
{0,1,2}^6 containing at least one 1 (the all-{0,2} vectors arise only from S = T,
which is excluded), so among the C(64,2) = 2016 pair products at most
3^6 - 2^6 = 665 DISTINCT values can occur -- whatever the u_i (zero, coincident,
negative, complex) and whatever the degeneracies.  Over characteristic zero,

    deg gcd(C_2, C_2') = 2016 - #distinct(pair products),
    C_2(z) = prod_{i<j} (z - mu_i mu_j),

so ANY six-mode Gaussian satisfies deg gcd(C_2, C_2') >= 2016 - 665 = 1351.
Contrapositive: a smaller gcd degree is an ordering-free certificate of
non-Gaussianity.  (Contrast `proofs/spectral_polynomial.md`: that candidate
polynomial condition was [FALSIFIED]; the collision bound above is proved, and
it survives the all-degenerate control u_i = 1, where all pair products collide
to a single value.)

INTEGRAL MODEL -- reduction can only ADD common factors, at EVERY prime.
R is the exact rational 2x3 layer matrix of e38 (imported, not duplicated).  Let
D > 0 be the exact least common multiple of all entry denominators of R and
R_int = D R in Z^{64x64}, asserted integral entry by entry.  (Do NOT guess the
shape of D: at t = 1/3 the exponents m = (b_k - eps)/2 are negative as well, so
denominators contain powers of 5, not just of 3; the observed factorization is
stored in the artifact.)  Then

    C_2^int(z) := prod_{i<j} (z - D^2 mu_i mu_j) = charpoly( wedge^2 R_int )

is MONIC with INTEGER coefficients.  Its monic characteristic-zero gcd g with the
formal derivative lies in Z[z] (Gauss) and divides both factors with integral
quotients, hence deg g and both divisibilities survive reduction modulo ANY
prime p:  deg gcd_{F_p} >= deg gcd_Q for every p, with no exceptional
denominator primes.  (Reducing the rational coefficients of C_2 directly, while
merely avoiding the entry denominators, would NOT give this: the gcd's own
factor denominators could be hit by p.  This is why the integer scaling is used
and stored.)  Collisions are invariant under the scaling mu -> D^2 mu.

PIPELINE (all exact, modular, at a prime p > 2016):
  1. Rb = R_int mod p; power traces s_k = tr(Rb^k) for k = 1..64 (int64 matmul
     with the proved overflow guard 64 (p-1)^2 < 2^63).
  2. Newton identities -> charpoly of Rb (monic, degree 64); validated by
     Horner evaluation: charpoly(Rb) == 0 as a matrix.
  3. Extend s_k to k = 4032 = 2*2016 by the exact Cayley--Hamilton trace
     recurrence; validated against direct binary-powering traces at 13 indices
     k > 64.
  4. Pair power sums q_m = (s_m^2 - s_{2m}) / 2, m = 1..2016.
  5. Newton identities -> C_2 mod p (monic, degree 2016).
  6. Explicit Euclidean gcd of C_2 with its formal derivative in F_p[z]; the
     result is verified to divide both inputs with zero remainder.

One decisive prime with deg gcd < 1351 proves deg gcd_Q < 1351, i.e. more than
665 distinct pair products, i.e. NOT a full six-mode Gaussian spectrum.  A second
prime is an independent implementation cross-check; mathematically (integral
model) a single prime already suffices.

CONTROLS.
* synthetic six-mode Gaussian: the diagonal family with u = 2,3,5,7,11,13
  (a = 1) conjugated by a unimodular integer matrix.  Its pipeline charpoly and
  pair polynomial must equal the DIRECT polynomial products over the known
  eigenvalues, and its gcd degree must be >= 1351 (exactly 1351 over Q).
* 1D open chain n = 6 (2D Ising strip): the physics free-fermion control,
  expected >= 1351 (consistency only, same caveat as Theorem S's control).
* small-size ground truth: the pair reconstruction equals the determinant-
  interpolated characteristic polynomial of the exterior square wedge^2(M) for a
  random integer matrix and for the chain n = 3 -- an independent route with no
  Newton identities and no power sums.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import random
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from math import gcd as math_gcd, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from ising.transfer_matrix import layer_bonds  # noqa: E402
from e38_gaussianity_certificate import build_R  # noqa: E402  the frozen rational construction

INT64_MAX = 2**63 - 1
DIM = 64                                   # 2^6
PAIRS = DIM * (DIM - 1) // 2               # C(64,2) = 2016
EXPONENT_VECTORS = 3**6 - 2**6             # 665
GAUSSIAN_GCD_FLOOR = PAIRS - EXPONENT_VECTORS   # 1351
MAX_POWER = 2 * PAIRS                      # 4032
SPOT_KS = (65, 66, 67, 100, 128, 256, 512, 1000, 1024, 2016, 2017, 3000, 4032)

CHECK_RECORDS: list[dict] = []
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    CHECK_RECORDS.append({"name": name, "passed": bool(ok), "detail": detail})
    if not ok:
        FAILS.append(name)
    return bool(ok)


# --------------------------------------------------------------------------- primes
def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            return False
    return True


def first_prime_ge(n: int) -> int:
    while not is_prime(n):
        n += 1
    return n


# ------------------------------------------------------------------- modular algebra
def assert_matmul_guard(dim: int, p: int) -> None:
    assert dim * (p - 1) ** 2 < INT64_MAX, "int64 matmul overflow guard violated"


def sequential_power_traces(Rp: np.ndarray, kmax: int, p: int) -> list[int]:
    """[0, tr(R), tr(R^2), ..., tr(R^kmax)] mod p by explicit products."""
    dim = Rp.shape[0]
    assert_matmul_guard(dim, p)
    out = [0]
    M = np.eye(dim, dtype=np.int64)
    for _ in range(kmax):
        M = (M @ Rp) % p
        out.append(int(M.trace()) % p)
    return out


def trace_power_binary(Rp: np.ndarray, k: int, p: int) -> int:
    """tr(R^k) mod p by square-and-multiply (independent of the sequential route)."""
    dim = Rp.shape[0]
    assert_matmul_guard(dim, p)
    R = np.eye(dim, dtype=np.int64)
    B = Rp.copy()
    e = k
    while e:
        if e & 1:
            R = (R @ B) % p
        e >>= 1
        if e:
            B = (B @ B) % p
    return int(R.trace()) % p


def newton_from_power_sums(q: np.ndarray, N: int, p: int, inv: list[int]) -> np.ndarray:
    """Monic P(z) = z^N + a_1 z^{N-1} + ... + a_N from power sums q[1..N]
    (index 0 unused) via Newton identities:
        k a_k = -( sum_{i=1..k} a_{k-i} q_i ),  a_0 = 1.
    Returns the DESCENDING coefficient array [1, a_1, ..., a_N] (int64 mod p)."""
    assert N * (p - 1) ** 2 < INT64_MAX, "int64 Newton dot overflow guard violated"
    a = np.zeros(N + 1, dtype=np.int64)
    a[0] = 1
    for k in range(1, N + 1):
        acc = int(np.dot(q[1:k + 1], np.ascontiguousarray(a[k - 1::-1]))) % p
        a[k] = (-acc * inv[k]) % p
    return a


def extend_power_sums(s: list[int], a: np.ndarray, upto: int, p: int) -> list[int]:
    """Extend s[k] = tr(R^k) beyond the charpoly degree N = len(a)-1 using the
    exact Cayley--Hamilton recurrence s_k = -sum_{i=1..N} a_i s_{k-i}."""
    N = a.size - 1
    assert N * (p - 1) ** 2 < INT64_MAX
    coeffs = a[1:].astype(np.int64)
    for k in range(N + 1, upto + 1):
        # sum_{i=1..N} a_i * s_{k-i} pairs a_1 with s_{k-1}, ..., a_N with s_{k-N}
        window = np.array(s[k - 1:k - N - 1:-1], dtype=np.int64)
        acc = int(np.dot(coeffs, window)) % p
        s.append((-acc) % p)
    return s


def horner_matrix_zero(chi_desc: np.ndarray, Rp: np.ndarray, p: int) -> bool:
    """Cayley--Hamilton validation: chi(Rp) == 0 mod p."""
    dim = Rp.shape[0]
    assert_matmul_guard(dim, p)
    H = int(chi_desc[0]) * np.eye(dim, dtype=np.int64) % p
    for c in chi_desc[1:]:
        H = ((H @ Rp) + int(c) * np.eye(dim, dtype=np.int64)) % p
    return bool(np.all(H == 0))


# ------------------------------------------------------------------ F_p[z] utilities
def poly_trim(a: np.ndarray) -> np.ndarray:
    nz = np.flatnonzero(a)
    if nz.size == 0:
        return np.zeros(1, dtype=np.int64)
    return a[nz[0]:]


def poly_is_zero(a: np.ndarray) -> bool:
    return a.size == 1 and int(a[0]) == 0


def poly_monic(a: np.ndarray, p: int) -> np.ndarray:
    a = poly_trim(np.asarray(a, dtype=np.int64) % p)
    if poly_is_zero(a):
        return a
    return (a * pow(int(a[0]), -1, p)) % p


def poly_rem(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """Remainder of a modulo the MONIC polynomial b.  Both DESCENDING arrays.
    Synthetic division: each pivot costs one vectorized multiply-subtract."""
    assert int(b[0]) == 1, "divisor must be monic"
    r = (np.asarray(a, dtype=np.int64) % p).copy()
    db, da = b.size - 1, r.size - 1
    if da < db:
        return poly_trim(r)
    for t in range(da - db + 1):
        c = int(r[t])
        if c:
            r[t:t + db + 1] = (r[t:t + db + 1] - c * b) % p
    return poly_trim(r[da - db + 1:])


def poly_gcd(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """Monic Euclidean gcd in F_p[z], explicit."""
    a = poly_monic(a, p)
    b = poly_monic(b, p)
    while not poly_is_zero(b):
        a, b = b, poly_monic(poly_rem(a, b, p), p)
    return a


def poly_deriv_desc(P_desc: np.ndarray, p: int) -> np.ndarray:
    M = P_desc.size - 1
    asc = P_desc[::-1]
    d_asc = (np.arange(1, M + 1, dtype=np.int64) * asc[1:]) % p
    return np.ascontiguousarray(d_asc[::-1])


def coeff_digest(coeffs) -> str:
    """SHA256 of the ASCENDING coefficient sequence c_0, c_1, ..., c_N (c_N = 1)."""
    return hashlib.sha256(",".join(str(int(c)) for c in coeffs).encode()).hexdigest()


def poly_from_roots(roots, p: int) -> np.ndarray:
    """prod (z - r) mod p, ascending coefficients, vectorized update."""
    a = np.array([1], dtype=np.int64)
    for r in roots:
        b = np.zeros(a.size + 1, dtype=np.int64)
        b[1:] = a
        b[:a.size] = (b[:a.size] - int(r) * a) % p
        a = b % p
    return a


# --------------------------------------------- exterior square + determinant ground truth
def wedge2(M: np.ndarray, p: int) -> np.ndarray:
    """Matrix of the exterior square in the basis e_i ^ e_j (i < j):
    W[(i,j),(k,l)] = M[k,i] M[l,j] - M[l,i] M[k,j].  Its eigenvalues are exactly
    the pair products of the eigenvalues of M (with multiplicity)."""
    n = M.shape[0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    W = np.zeros((m, m), dtype=np.int64)
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            W[a, b] = (int(M[k, i]) * int(M[l, j]) - int(M[l, i]) * int(M[k, j])) % p
    return W


def det_modp(A: np.ndarray, p: int) -> int:
    """Exact determinant mod p by Gaussian elimination with pivoting (Python ints)."""
    n = A.shape[0]
    a = [[int(x) % p for x in row] for row in A]
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            det = -det
        det = det * a[c][c] % p
        inv = pow(a[c][c] % p, -1, p)
        for r in range(c + 1, n):
            f = a[r][c] * inv % p
            if f:
                rowc, ar = a[c], a[r]
                for j in range(c, n):
                    ar[j] = (ar[j] - f * rowc[j]) % p
    return det % p


def charpoly_by_interpolation(M: np.ndarray, p: int) -> list[int]:
    """det(z I - M) in F_p[z], ASCENDING coefficients, by evaluating at
    deg+1 integer points and Lagrange interpolating.  No Newton identities, no
    power sums -- an independent ground-truth route."""
    m = M.shape[0]
    xs = list(range(m + 1))
    ys = []
    for x in xs:
        A = ((x % p) * np.eye(m, dtype=np.int64) - (M % p)) % p
        ys.append(det_modp(A, p))
    total = [0] * (m + 1)
    for i, x in enumerate(xs):
        denom = 1
        for j, xj in enumerate(xs):
            if j != i:
                denom = denom * ((x - xj) % p) % p
        scale = ys[i] * pow(denom, -1, p) % p
        basis = [1]
        for j, xj in enumerate(xs):
            if j == i:
                continue
            nb = [0] * (len(basis) + 1)
            for d, c in enumerate(basis):
                nb[d + 1] = (nb[d + 1] + c) % p
                nb[d] = (nb[d] - c * xj) % p
            basis = nb
        for d, c in enumerate(basis):
            total[d] = (total[d] + scale * c) % p
    return total


# ------------------------------------------------------------------------ constructions
def integralize(R) -> tuple[int, list[list[int]]]:
    """D = exact lcm of ALL entry denominators of R (computed, never guessed:
    at t = 1/3 the exponents m = (b_k - eps)/2 are negative for some states, so
    5-powers enter alongside 3-powers); R_int = D R, asserted integral entrywise."""
    D = 1
    for row in R:
        for f in row:
            D = D * f.denominator // math_gcd(D, f.denominator)
    assert D > 0
    R_int = []
    for row in R:
        r = []
        for f in row:
            g = D * f
            assert isinstance(g, Fraction) and g.denominator == 1, (
                f"D*R not integral at entry with denominator {f.denominator} (D = {D})")
            r.append(int(g))
        R_int.append(r)
    return D, R_int


def factorize(n: int) -> dict[int, int]:
    f: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def synthetic_gaussian_matrix(seed: int = 20260816) -> tuple[list[list[int]], list[int]]:
    """Integral matrix similar (over Z) to diag( prod_{i in S} u_i ), u = first six
    primes, a = 1: a bona fide full six-mode Gaussian spectrum with exactly
    3^6 - 2^6 = 665 distinct pair products over Q (unique factorization).  A
    product of unimodular shears makes it a NON-diagonal, NON-symmetric integer
    matrix, so the pipeline exercises a generic integer input."""
    u = [2, 3, 5, 7, 11, 13]
    dim = 1 << 6
    vals = []
    for e in range(dim):
        v = 1
        for i in range(6):
            if (e >> i) & 1:
                v *= u[i]
        vals.append(v)
    rng = random.Random(seed)
    U = [[1 if r == c else 0 for c in range(dim)] for r in range(dim)]
    Uinv = [row[:] for row in U]
    for _ in range(8):
        i = rng.randrange(dim)
        j = rng.randrange(dim)
        while j == i:
            j = rng.randrange(dim)
        s = rng.choice((-2, -1, 1, 2))
        # U <- U (I + s E_ij):  column j of U += s * column i   (det stays 1)
        for r in range(dim):
            U[r][j] += s * U[r][i]
        # Uinv <- (I - s E_ij) Uinv:  row i of Uinv -= s * row j   (exact inverse)
        for c in range(dim):
            Uinv[i][c] -= s * Uinv[j][c]
    G = [[sum(U[r][k] * vals[k] * Uinv[k][c] for k in range(dim)) for c in range(dim)]
         for r in range(dim)]
    return G, vals


# ---------------------------------------------------------------------------- pipeline
def run_pipeline(M_rows: list[list[int]], p: int) -> dict:
    """Full trace -> Newton -> recurrence -> pair power sums -> Newton -> gcd
    pipeline at prime p for an integral matrix.  Returns coefficient lists for
    further direct comparison, but only digests should be stored in artifacts."""
    t0 = time.time()
    dim = len(M_rows)
    pairs = dim * (dim - 1) // 2
    assert p > pairs and p > dim and p % 2 == 1, "prime must exceed both degrees and be odd"
    Rp = np.array([[int(x) % p for x in row] for row in M_rows], dtype=np.int64)

    # 1. power traces through dim
    s = sequential_power_traces(Rp, dim, p)

    # 2. charpoly by Newton; validate by Cayley--Hamilton (Horner)
    inv = [0] + [pow(k, -1, p) for k in range(1, pairs + 1)]
    chi_desc = newton_from_power_sums(np.array(s, dtype=np.int64), dim, p, inv)
    horner_ok = horner_matrix_zero(chi_desc, Rp, p)

    # 3. extend traces to 2*pairs by the exact recurrence; validate by binary powering
    s = extend_power_sums(s, chi_desc, 2 * pairs, p)
    spots = []
    for k in SPOT_KS:
        if k <= 2 * pairs:
            d = trace_power_binary(Rp, k, p)
            spots.append({"k": k, "direct": d, "recurrence": int(s[k]), "match": d == int(s[k])})

    # 4. pair power sums  q_m = (s_m^2 - s_{2m}) / 2
    inv2 = pow(2, -1, p)
    q = np.zeros(pairs + 1, dtype=np.int64)
    for m in range(1, pairs + 1):
        q[m] = ((int(s[m]) * int(s[m]) - int(s[2 * m])) % p) * inv2 % p

    # 5. pair polynomial by Newton
    P_desc = newton_from_power_sums(q, pairs, p, inv)

    # 6. gcd with the formal derivative
    dP_desc = poly_deriv_desc(P_desc, p)
    g = poly_gcd(P_desc, dP_desc, p)
    deg_gcd = 0 if poly_is_zero(g) else int(g.size) - 1
    divides = poly_is_zero(poly_rem(P_desc, g, p)) and poly_is_zero(poly_rem(dP_desc, g, p))
    monic = poly_is_zero(g) or int(g[0]) == 1

    chi_asc = [int(c) for c in (chi_desc[::-1] % p)]
    P_asc = [int(c) for c in (P_desc[::-1] % p)]
    assert int(chi_desc[0]) == 1 and int(P_desc[0]) == 1
    return {
        "prime": p,
        "charpoly_degree": dim,
        "pair_poly_degree": pairs,
        "charpoly_monic": True,
        "pair_poly_monic": True,
        "charpoly_sha256": coeff_digest(chi_asc),
        "pair_poly_sha256": coeff_digest(P_asc),
        "gcd_sha256": coeff_digest([int(c) for c in (g[::-1] % p)]),
        "gcd_degree": deg_gcd,
        "gcd_divides_both": bool(divides),
        "gcd_monic": bool(monic),
        "horner_charpoly_zero_matrix": bool(horner_ok),
        "spot_trace_checks": spots,
        "spot_trace_checks_all_match": all(x["match"] for x in spots),
        "elapsed_seconds": round(time.time() - t0, 2),
        "_charpoly_asc": chi_asc,
        "_pair_asc": P_asc,
    }


def wedge_ground_truth(M_rows: list[list[int]], p: int) -> dict:
    """Small-size validation: the pair polynomial reconstructed from traces of M
    must equal the determinant-interpolated charpoly of wedge^2(M), an
    independent route with no Newton identities and no power sums."""
    dim = len(M_rows)
    pairs = dim * (dim - 1) // 2
    M = np.array([[int(x) % p for x in row] for row in M_rows], dtype=np.int64)
    W = wedge2(M, p)
    gt_pair_asc = charpoly_by_interpolation(W, p)      # ground truth
    gt_chi_asc = charpoly_by_interpolation(M, p)       # extra: charpoly ground truth

    inv = [0] + [pow(k, -1, p) for k in range(1, pairs + 1)]
    inv2 = pow(2, -1, p)
    s_direct = sequential_power_traces(M, 2 * pairs, p)     # direct traces, no recurrence
    chi_desc = newton_from_power_sums(np.array(s_direct[:dim + 1]), dim, p, inv)
    s = extend_power_sums(s_direct[:dim + 1], chi_desc, 2 * pairs, p)
    recurrence_ok = all(s[k] == s_direct[k] for k in range(dim + 1, 2 * pairs + 1))
    q = np.zeros(pairs + 1, dtype=np.int64)
    for m in range(1, pairs + 1):
        q[m] = ((int(s[m]) * int(s[m]) - int(s[2 * m])) % p) * inv2 % p
    P_desc = newton_from_power_sums(q, pairs, p, inv)
    P_asc = [int(c) for c in (P_desc[::-1] % p)]

    # power sums of the pair products must also equal traces of powers of W
    wtrace = sequential_power_traces(W % p, pairs, p)
    q_match_wedge = all(int(q[m]) == wtrace[m] for m in range(1, pairs + 1))

    return {
        "dimension": dim,
        "pair_degree": pairs,
        "charpoly_matches_det_interpolation":
            [int(c) for c in (chi_desc[::-1] % p)] == gt_chi_asc,
        "pair_poly_matches_wedge2_det_interpolation": P_asc == gt_pair_asc,
        "pair_power_sums_match_wedge2_traces": bool(q_match_wedge),
        "trace_recurrence_matches_direct": bool(recurrence_ok),
    }


# ------------------------------------------------------------------------------- driver
def main() -> int:
    t_start = time.time()
    t = Fraction(1, 3)

    bonds_layer = list(layer_bonds((2, 3), (False, False)))
    bonds_chain = [(i, i + 1) for i in range(5)]

    print("building exact rational matrices (e38 construction, t = 1/3) ...", flush=True)
    R_layer, q_exp, eps_l = build_R(6, bonds_layer, t)
    R_chain, _, eps_c = build_R(6, bonds_chain, t)
    R_small, _, eps_s = build_R(3, [(0, 1), (1, 2)], t)

    D_layer, Rl_int = integralize(R_layer)
    D_chain, Rc_int = integralize(R_chain)
    D_small, Rs_int = integralize(R_small)
    G_syn, syn_vals = synthetic_gaussian_matrix()

    print(f"   bond parities eps: layer {eps_l}, chain {eps_c}, chain n=3 {eps_s}")
    print(f"   D(2x3 layer) = {D_layer} = {factorize(D_layer)}")
    print(f"   D(chain n=6) = {D_chain} = {factorize(D_chain)}")
    print(f"   D(chain n=3) = {D_small} = {factorize(D_small)}")

    p1 = first_prime_ge(1_000_000)
    p2 = first_prime_ge(2_000_000)
    primes = [p1, p2]
    guards = {
        "matmul_guard": "dim*(p-1)^2 < 2^63",
        "newton_guard": "N*(p-1)^2 < 2^63 for dot lengths N <= 2016",
        "values": {str(p): {"dim_times_pm1_sq": DIM * (p - 1) ** 2,
                            "pairs_times_pm1_sq": PAIRS * (p - 1) ** 2,
                            "int64_max": INT64_MAX} for p in primes},
    }
    for p in primes:
        assert is_prime(p) and p > PAIRS and p % 2 == 1
        assert DIM * (p - 1) ** 2 < INT64_MAX and PAIRS * (p - 1) ** 2 < INT64_MAX

    # ---------------------------------------------------------------- small ground truth
    print("small-size exterior-square ground truth (no Newton, no power sums) ...", flush=True)
    rng = random.Random(1234)
    M_rand = [[rng.randrange(-50, 51) for _ in range(6)] for _ in range(6)]
    gt_results = {}
    for p in primes:
        gt_results[str(p)] = {
            "random_6x6": wedge_ground_truth(M_rand, p),
            "chain_n3_layer": wedge_ground_truth(Rs_int, p),
        }
        for case, r in gt_results[str(p)].items():
            check(f"ground_truth[{case}] pair poly == det-interpolated charpoly(wedge^2) @ p={p}",
                  r["pair_poly_matches_wedge2_det_interpolation"])
            check(f"ground_truth[{case}] charpoly == det-interpolated charpoly @ p={p}",
                  r["charpoly_matches_det_interpolation"])
            check(f"ground_truth[{case}] pair power sums == tr(wedge^2 powers) @ p={p}",
                  r["pair_power_sums_match_wedge2_traces"])
            check(f"ground_truth[{case}] trace recurrence == direct traces @ p={p}",
                  r["trace_recurrence_matches_direct"])

    # ---------------------------------------------------------------------- main runs
    cases = {
        "layer_2x3": {"matrix": Rl_int, "D": D_layer, "eps": eps_l,
                      "desc": "2x3 open grid (3D Ising layer), t=1/3, THE target case"},
        "chain_n6_control": {"matrix": Rc_int, "D": D_chain, "eps": eps_c,
                             "desc": "1D open chain n=6 (2D Ising strip) -- free-fermion control"},
        "synthetic_gaussian_6mode": {"matrix": G_syn, "D": 1, "eps": None,
                                     "desc": "unimodular conjugate of diag(prod of primes), a=1"},
    }
    results = {}
    for name, spec in cases.items():
        results[name] = {"description": spec["desc"], "D": spec["D"],
                         "exp_2K": str(q_exp) if spec["eps"] is not None else "n/a",
                         "bond_parity": spec["eps"],
                         "per_prime": {}}
        for p in primes:
            print(f"pipeline {name} @ p = {p} ...", flush=True)
            r = run_pipeline(spec["matrix"], p)
            results[name]["per_prime"][str(p)] = r
            print(f"   gcd(C_2, C_2') degree = {r['gcd_degree']}   "
                  f"({r['elapsed_seconds']} s)")
        gcds = [results[name]["per_prime"][str(p)]["gcd_degree"] for p in primes]
        gmin = min(gcds)
        results[name]["gcd_degree_min"] = gmin
        results[name]["gcd_degrees"] = gcds
        results[name]["distinct_pair_products_lower_bound"] = PAIRS - gmin
        results[name]["gaussian_gcd_floor"] = GAUSSIAN_GCD_FLOOR
        results[name]["decisive_non_gaussian"] = bool(gmin < GAUSSIAN_GCD_FLOOR)

    # ------------------------------------------------- synthetic full-size direct products
    syn_direct = {}
    for p in primes:
        r = results["synthetic_gaussian_6mode"]["per_prime"][str(p)]
        vals = [v % p for v in syn_vals]
        chi_direct = [int(c) for c in poly_from_roots(vals, p)]
        pair_roots = [(vals[i] * vals[j]) % p for i in range(DIM) for j in range(i + 1, DIM)]
        pair_direct = [int(c) for c in poly_from_roots(pair_roots, p)]
        syn_direct[str(p)] = {
            "charpoly_equals_direct_product": r["_charpoly_asc"] == chi_direct,
            "pair_poly_equals_direct_product": r["_pair_asc"] == pair_direct,
        }
    results["synthetic_gaussian_6mode"]["direct_product_checks"] = syn_direct

    # ---------------------------------------------------------------------- final checks
    for p in primes:
        for name in cases:
            r = results[name]["per_prime"][str(p)]
            check(f"{name}@p={p}: charpoly monic degree 64, pair poly monic degree 2016",
                  r["charpoly_monic"] and r["charpoly_degree"] == DIM
                  and r["pair_poly_monic"] and r["pair_poly_degree"] == PAIRS)
            check(f"{name}@p={p}: Cayley-Hamilton Horner charpoly(Rb) == 0",
                  r["horner_charpoly_zero_matrix"])
            check(f"{name}@p={p}: all {len(r['spot_trace_checks'])} trace-recurrence spot "
                  f"checks match binary powering", r["spot_trace_checks_all_match"])
            check(f"{name}@p={p}: gcd divides both inputs, monic",
                  r["gcd_divides_both"] and r["gcd_monic"])
        check(f"synthetic@p={p}: charpoly == direct product over known eigenvalues",
              syn_direct[str(p)]["charpoly_equals_direct_product"])
        check(f"synthetic@p={p}: pair polynomial == direct pair-product product",
              syn_direct[str(p)]["pair_poly_equals_direct_product"])

    syn_gcds = results["synthetic_gaussian_6mode"]["gcd_degrees"]
    check(f"synthetic six-mode Gaussian control: gcd degree >= {GAUSSIAN_GCD_FLOOR} "
          f"(got {syn_gcds}; exactly {GAUSSIAN_GCD_FLOOR} over Q by unique factorization)",
          all(g >= GAUSSIAN_GCD_FLOOR for g in syn_gcds))
    ch_gcds = results["chain_n6_control"]["gcd_degrees"]
    check(f"1D chain n=6 control: gcd degree >= {GAUSSIAN_GCD_FLOOR} (got {ch_gcds}) "
          f"[consistency only]", all(g >= GAUSSIAN_GCD_FLOOR for g in ch_gcds))
    lay = results["layer_2x3"]
    check(f"2x3 layer DECISIVE: gcd degree < {GAUSSIAN_GCD_FLOOR} at BOTH primes "
          f"(got {lay['gcd_degrees']})",
          all(g < GAUSSIAN_GCD_FLOOR for g in lay["gcd_degrees"]))
    check(f"2x3 layer: at least {lay['distinct_pair_products_lower_bound']} > {EXPONENT_VECTORS} "
          f"distinct pair products over Q",
          lay["distinct_pair_products_lower_bound"] > EXPONENT_VECTORS)

    # strip the private coefficient lists before storing (digests + algorithm suffice)
    for name in cases:
        for p in primes:
            results[name]["per_prime"][str(p)].pop("_charpoly_asc", None)
            results[name]["per_prime"][str(p)].pop("_pair_asc", None)

    lay_gcds = lay["gcd_degrees"]
    out = {
        "provenance": {
            "script": "experiments/e107_pair_product_obstruction.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": "exact modular arithmetic (int64 matmul with proved overflow "
                         "guards dim*(p-1)^2 < 2^63; Python ints elsewhere); no floating "
                         "point anywhere in the certificate. Coefficient sequences stored "
                         "as SHA256 digests over the ASCENDING order "
                         "c_0 + c_1 z + ... + c_N z^N with c_N = 1.",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "construction": "experiments/e38_gaussianity_certificate.build_R imported "
                            "(frozen rational layer construction, t = 1/3); no second "
                            "transfer convention.",
        },
        "data": {
            "threshold": {
                "pair_products": PAIRS,
                "exponent_vectors_with_a_one": EXPONENT_VECTORS,
                "gaussian_gcd_floor": GAUSSIAN_GCD_FLOOR,
                "derivation": "pair products of a full 6-mode Gaussian spectrum have "
                              "exponent vectors in {0,1,2}^6 with at least one 1 "
                              "(all-{0,2} vectors need S = T, excluded); 3^6 - 2^6 = 665 "
                              "classes; deg gcd(C_2, C_2') = 2016 - #distinct(pair "
                              "products) over characteristic 0, multiplicities included.",
            },
            "integral_model": {
                "note": "R_int = D * R with D the exact lcm of ALL entry denominators "
                        "(computed from the Fractions, never guessed); C_2^int = "
                        "charpoly(wedge^2 R_int) is MONIC in Z[z], so the monic char-0 "
                        "gcd with the derivative lies in Z[z] (Gauss) with integral "
                        "quotients, and deg gcd mod p >= deg gcd over Q at EVERY prime; "
                        "reduction can add common factors, never remove them. No "
                        "exceptional denominator primes exist for this model.",
                "D_layer_2x3": str(D_layer),
                "D_layer_2x3_factorization": {str(k): v for k, v in factorize(D_layer).items()},
                "D_chain_n6": str(D_chain),
                "D_chain_n6_factorization": {str(k): v for k, v in factorize(D_chain).items()},
                "D_chain_n3": str(D_small),
                "D_chain_n3_factorization": {str(k): v for k, v in factorize(D_small).items()},
                "scaling_invariance": "eigenvalues scale by D, pair products by D^2; "
                                      "collision multiplicities and both gcd degrees are "
                                      "invariant under nonzero scaling.",
            },
            "primes": {"p1": p1, "p2": p2, "primality": "exact trial division to sqrt",
                       "constraints": "p > 2016 (Newton divisions by k <= 2016 and by 2), "
                                      "p odd", "guards": guards},
            "pipeline": [
                "s_k = tr(Rb^k), k=1..64, by int64 modular matrix products",
                "charpoly by Newton identities; validated by Horner (Cayley-Hamilton)",
                "s_k extended to 4032 by the exact charpoly trace recurrence; validated "
                f"against binary-powering traces at k in {list(SPOT_KS)}",
                "q_m = (s_m^2 - s_{2m})/2, m=1..2016",
                "C_2 by Newton identities (monic degree 2016)",
                "explicit Euclidean gcd with the formal derivative in F_p[z], verified "
                "to divide both inputs",
            ],
            "small_ground_truth": gt_results,
            "cases": results,
            "verdict": {
                "layer_2x3_non_gaussian_by_pair_products": bool(
                    lay["gcd_degree_min"] < GAUSSIAN_GCD_FLOOR),
                "gcd_degrees_layer_2x3": lay_gcds,
                "distinct_pair_products_lower_bound": lay["distinct_pair_products_lower_bound"],
                "gaussian_maximum_distinct": EXPONENT_VECTORS,
                "scope": "full six-mode Gaussian subset-product spectra of the 2x3 open "
                         "layer at the single exact rational coupling t=1/3. NOT parity "
                         "restrictions, NOT restrictions to submultisets of larger "
                         "Gaussians, NOT other couplings or layers.",
                "relation_to_theorem_S": "independent second certificate; ordering-free "
                                         "(no sorted eigenvalues, no enclosures); one "
                                         "decisive prime suffices by the integral-model "
                                         "lemma, the second prime is an implementation "
                                         "cross-check.",
            },
            "wall_seconds": round(time.time() - t_start, 2),
        },
        "checks": CHECK_RECORDS,
    }
    os.makedirs("results/spectral", exist_ok=True)
    with open("results/spectral/pair_product_obstruction.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print()
    print(f"threshold: pairs = {PAIRS}, Gaussian max distinct = {EXPONENT_VECTORS}, "
          f"gcd floor = {GAUSSIAN_GCD_FLOOR}")
    print(f"2x3 layer gcd degrees: {lay_gcds}")
    print(f"chain n=6 control gcd degrees: {ch_gcds}")
    print(f"synthetic control gcd degrees: {syn_gcds}")
    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(f"PASS: ordering-free pair-product certificate -- the 2x3 layer spectrum has "
          f"more than {EXPONENT_VECTORS} distinct pair products, exceeding the six-mode "
          f"Gaussian maximum; total wall {out['data']['wall_seconds']} s.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
