"""Standalone verifier for the pair-product Gaussianity obstruction (e107).

Independently rebuilds the decisive prime certificate for the 2x3 3D-Ising layer
from the RAW rational transfer construction: nothing is imported from the
producer (experiments/e107_pair_product_obstruction.py) or from e38.  The exact
rational matrix R = P_t diag(q^m) P_t at t = 1/3 is rebuilt here from scratch,
the integer model R_int = D R is formed from the actual Fraction denominators,
and the whole modular pipeline (traces -> Newton -> recurrence -> pair power
sums -> Newton -> Euclidean gcd) is re-implemented with independent code (pure
Python polynomial arithmetic; numpy is used only for int64 matrix products under
the proved overflow guard dim*(p-1)^2 < 2^63).

The results JSON is read ONLY for regression comparison of independently
recomputed quantities (gcd degrees, digests), never to derive them.

Verified here:
  * threshold combinatorics: C(64,2) = 2016, 3^6 - 2^6 = 665, floor 1351;
  * Newton identity sign convention on a known cubic;
  * pair power-sum formula q_m = (s_m^2 - s_{2m})/2 against direct enumeration;
  * trace recurrence against direct binary-powering traces beyond degree 64;
  * Cayley--Hamilton Horner check charpoly(Rb) == 0;
  * pair polynomial == determinant-interpolated charpoly of wedge^2(M) (an
    independent route with no Newton identities and no power sums);
  * full decisive prime: monic degrees, gcd degree 385 < 1351 for the 2x3 layer
    (=> at least 1631 > 665 distinct pair products over Q => NOT a full
    six-mode Gaussian spectrum), regression-pinned to the producer artifact;
  * 1D open chain n=6 control: gcd degree 1351 >= 1351 (consistency only);
  * exact four-mode threshold: a synthetic Gaussian family with u = 2,3,5,7 has
    exactly 3^4 - 2^4 = 65 distinct pair products, so its gcd degree is exactly
    120 - 65 = 55 -- the threshold combinatorics, end to end.

Run:  PYTHONPATH=src .venv/bin/python tests/test_pair_product_obstruction.py
"""

from __future__ import annotations

import hashlib
import json
import random
import sys
from fractions import Fraction
from math import gcd as math_gcd, isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ising.transfer_matrix import layer_bonds  # noqa: E402  repo API only, no producer code

INT64_MAX = 2**63 - 1
DIM = 64
PAIRS = 2016
EXPONENT_VECTORS = 665
GAUSSIAN_GCD_FLOOR = 1351

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)
    return bool(ok)


# ------------------------------------------------------------------ tiny exact helpers
def is_prime(n: int) -> bool:
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def first_prime_ge(n: int) -> int:
    while not is_prime(n):
        n += 1
    return n


def mat_mul_mod(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    dim = A.shape[0]
    assert dim * (p - 1) ** 2 < INT64_MAX, "int64 matmul overflow guard violated"
    return (A @ B) % p


def mat_pow_trace(M: np.ndarray, k: int, p: int) -> int:
    """tr(M^k) mod p by square-and-multiply."""
    dim = M.shape[0]
    R = np.eye(dim, dtype=np.int64)
    B = M.copy()
    e = k
    while e:
        if e & 1:
            R = mat_mul_mod(R, B, p)
        e >>= 1
        if e:
            B = mat_mul_mod(B, B, p)
    return int(R.trace()) % p


# ------------------------------------------- polynomial arithmetic in F_p[z] (pure Python)
def newton_coeffs(power_sums: list[int], N: int, p: int) -> list[int]:
    """Monic P(z) = z^N + a_1 z^{N-1} + ... + a_N from s_1..s_N (1-based list
    with a dummy at index 0) via k a_k = -(sum_{i=1..k} a_{k-i} s_i).
    Returns DESCENDING [1, a_1, ..., a_N].  Pure Python ints (no overflow)."""
    a = [1]
    for k in range(1, N + 1):
        acc = 0
        for i in range(1, k + 1):
            acc += a[k - i] * power_sums[i]
        a.append(-acc % p * pow(k, -1, p) % p)
    return a


def extend_traces(s: list[int], a_desc: list[int], upto: int, p: int) -> list[int]:
    """s_k = -(a_1 s_{k-1} + ... + a_N s_{k-N}) for k > N (Cayley--Hamilton)."""
    N = len(a_desc) - 1
    tail = a_desc[1:]
    for k in range(N + 1, upto + 1):
        acc = sum(tail[i - 1] * s[k - i] for i in range(1, N + 1))
        s.append(-acc % p)
    return s


def horner_matrix(chi_desc: list[int], M: np.ndarray, p: int) -> np.ndarray:
    H = np.eye(M.shape[0], dtype=np.int64) * chi_desc[0] % p
    for c in chi_desc[1:]:
        H = (mat_mul_mod(H, M, p) + c * np.eye(M.shape[0], dtype=np.int64)) % p
    return H


def poly_trim(a: list[int]) -> list[int]:
    i = 0
    while i < len(a) - 1 and a[i] == 0:
        i += 1
    return a[i:] if any(a[i:]) else [0]


def poly_div_rem(a: list[int], b: list[int], p: int) -> list[int]:
    """Remainder of a (descending) modulo b; b is made monic first."""
    ilb = pow(b[0] % p, -1, p)
    b = [x * ilb % p for x in b]
    r = [x % p for x in a]
    db = len(b) - 1
    while len(r) - 1 >= db and any(r):
        c = r[0]
        if c:
            for j in range(db + 1):
                r[j] = (r[j] - c * b[j]) % p
        r = poly_trim(r[1:])  # the monic divisor zeroes the processed leading slot
    return r


def poly_gcd(a: list[int], b: list[int], p: int) -> list[int]:
    a = poly_trim(a)
    a = [x * pow(a[0], -1, p) % p for x in a] if any(a) else [0]
    b = poly_trim(b)
    b = [x * pow(b[0], -1, p) % p for x in b] if any(b) else [0]
    while any(b):
        r = poly_div_rem(a, b, p)
        a, b = b, ([x * pow(r[0], -1, p) % p for x in r] if any(r) else [0])
    return a if any(a) else [1]


def digest(coeffs: list[int]) -> str:
    return hashlib.sha256(",".join(str(int(c)) for c in coeffs).encode()).hexdigest()


# ------------------------------------------- independent rational transfer construction
def build_R(n: int, bonds, t: Fraction):
    """The e38 rational layer construction, re-derived here (P_t diag(q^m) P_t)."""
    q = (1 + t * t) / (2 * t)
    dim = 1 << n
    b = []
    for k in range(dim):
        spins = [1 - 2 * ((k >> (n - 1 - i)) & 1) for i in range(n)]
        b.append(sum(spins[i] * spins[j] for i, j in bonds))
    eps = b[0] % 2
    assert all((x - eps) % 2 == 0 for x in b)
    d = [q ** ((x - eps) // 2) for x in b]
    tp = [t ** h for h in range(n + 1)]
    P = [[tp[bin(k ^ l).count("1")] for l in range(dim)] for k in range(dim)]
    R = [[Fraction(0)] * dim for _ in range(dim)]
    for i in range(dim):
        for j in range(i, dim):
            acc = Fraction(0)
            for k in range(dim):
                acc += P[i][k] * d[k] * P[k][j]
            R[i][j] = R[j][i] = acc
    return R, q, eps


def integralize(R):
    D = 1
    for row in R:
        for f in row:
            D = D * f.denominator // math_gcd(D, f.denominator)
    R_int = []
    for row in R:
        out = []
        for f in row:
            g = D * f
            assert g.denominator == 1
            out.append(int(g))
        R_int.append(out)
    return D, R_int


# --------------------------------------------------------------------- pipeline (test)
def pipeline(M_rows: list[list[int]], p: int) -> dict:
    dim = len(M_rows)
    pairs = dim * (dim - 1) // 2
    M = np.array([[x % p for x in row] for row in M_rows], dtype=np.int64)
    # direct sequential traces to dim
    T = np.eye(dim, dtype=np.int64)
    s = [0]
    for _ in range(dim):
        T = mat_mul_mod(T, M, p)
        s.append(int(T.trace()) % p)
    chi_desc = newton_coeffs(s, dim, p)
    horner_zero = bool(np.all(horner_matrix(chi_desc, M, p) == 0))
    s = extend_traces(s, chi_desc, 2 * pairs, p)
    inv2 = pow(2, -1, p)
    q = [0] + [((s[m] * s[m] - s[2 * m]) % p) * inv2 % p for m in range(1, pairs + 1)]
    P_desc = newton_coeffs(q, pairs, p)
    # derivative, descending: d/dz of z^M + a1 z^{M-1} + ... has coefficients
    # (M, (M-1) a_1, ..., a_{M-1})
    Mdeg = pairs
    dP_desc = [(Mdeg - i) % p * P_desc[i] % p for i in range(Mdeg)]
    g = poly_gcd(P_desc, dP_desc, p)
    deg_gcd = len(poly_trim(g)) - 1 if any(g) else 0
    return {
        "chi_desc": chi_desc,
        "P_desc": P_desc,
        "dP_desc": dP_desc,
        "gcd": g,
        "gcd_degree": deg_gcd,
        "horner_zero": horner_zero,
        "s": s,
        "M": M,
    }


def det_modp(A: np.ndarray, p: int) -> int:
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
                for j in range(c, n):
                    a[r][j] = (a[r][j] - f * a[c][j]) % p
    return det % p


def charpoly_interp(M: np.ndarray, p: int) -> list[int]:
    m = M.shape[0]
    xs = list(range(m + 1))
    ys = [det_modp((x * np.eye(m, dtype=np.int64) - M) % p, p) for x in xs]
    total = [0] * (m + 1)
    for i, x in enumerate(xs):
        denom = 1
        for j, xj in enumerate(xs):
            if j != i:
                denom = denom * (x - xj) % p
        scale = ys[i] * pow(denom, -1, p) % p
        basis = [1]
        for j, xj in enumerate(xs):
            if j == i:
                continue
            nb = [0] * (len(basis) + 1)
            for ddx, c in enumerate(basis):
                nb[ddx + 1] = (nb[ddx + 1] + c) % p
                nb[ddx] = (nb[ddx] - c * xj) % p
            basis = nb
        for ddx, c in enumerate(basis):
            total[ddx] = (total[ddx] + scale * c) % p
    return total


def wedge2(M: np.ndarray, p: int) -> np.ndarray:
    n = M.shape[0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    W = np.zeros((len(pairs), len(pairs)), dtype=np.int64)
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            W[a, b] = (int(M[k, i]) * int(M[l, j]) - int(M[l, i]) * int(M[k, j])) % p
    return W


# ====================================================================================
def main() -> int:
    p = first_prime_ge(1_000_000)
    print(f"decisive prime p = {p} (first prime >= 10^6, recomputed here)")
    check("prime constraints and int64 guard: p > 2016, 64*(p-1)^2 < 2^63",
          p > PAIRS and DIM * (p - 1) ** 2 < INT64_MAX)

    # ---------------------------------------------------------------- 1. combinatorics
    print("\n1. threshold combinatorics (exact integers)")
    check("C(64,2) = 2016", DIM * (DIM - 1) // 2 == 2016)
    check("3^6 - 2^6 = 665", 3**6 - 2**6 == 665)
    check("Gaussian gcd floor 2016 - 665 = 1351", PAIRS - EXPONENT_VECTORS == GAUSSIAN_GCD_FLOOR)

    # ------------------------------------------------------- 2. Newton sign convention
    print("\n2. Newton identities recover a known cubic exactly")
    roots = [2, 3, 5]
    s = [0] + [sum(r ** k for r in roots) for k in (1, 2, 3)]
    chi = newton_coeffs(s, 3, p)
    # (z-2)(z-3)(z-5) = z^3 - 10 z^2 + 31 z - 30
    check("Newton signs: z^3 - 10 z^2 + 31 z - 30 from power sums of {2,3,5}",
          chi == [1, -10 % p, 31 % p, -30 % p], f"got {chi}")

    # --------------------------------------------------- 3. pair power-sum formula
    print("\n3. pair power-sum formula q_m = (s_m^2 - s_{2m})/2 vs direct enumeration")
    vals = [1, 2, 3, 4]
    Dm = np.diag([v % p for v in vals]).astype(np.int64)
    s_direct = [0] + [mat_pow_trace(Dm, k, p) for k in range(1, 21)]
    ok = True
    for m in range(1, 11):
        direct = sum((vals[i] * vals[j]) ** m for i in range(4) for j in range(i + 1, 4)) % p
        formula = (s_direct[m] * s_direct[m] - s_direct[2 * m]) % p * pow(2, -1, p) % p
        ok = ok and direct == formula
    check("q_m formula matches direct pair enumeration for m = 1..10", ok)

    # ------------------------------------------------------------- 4. wedge^2 ground truth
    print("\n4. pair polynomial == det-interpolated charpoly of wedge^2(M) (dim 4)")
    rng = random.Random(99)
    Mrows = [[rng.randrange(-50, 51) for _ in range(4)] for _ in range(4)]
    M4 = np.array([[x % p for x in row] for row in Mrows], dtype=np.int64)
    gt_pair = charpoly_interp(wedge2(M4, p), p)          # ascending
    r4 = pipeline(Mrows, p)
    got_pair = [c % p for c in r4["P_desc"][::-1]]
    check("pair poly from traces == charpoly(wedge^2 M) by det interpolation",
          got_pair == gt_pair, f"degrees {len(got_pair) - 1} vs {len(gt_pair) - 1}")
    chi_gt = charpoly_interp(M4, p)
    check("charpoly from traces == det-interpolated charpoly",
          [c % p for c in r4["chi_desc"][::-1]] == chi_gt)

    # ------------------------------------------------- 5. independent 2x3 rebuild + prime
    print("\n5. independent rational rebuild of the 2x3 layer at t = 1/3")
    t = Fraction(1, 3)
    bonds23 = [(0, 1), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (4, 5)]
    check("hardcoded 2x3 bond set == layer_bonds((2,3), open)",
          set(bonds23) == set(layer_bonds((2, 3), (False, False))))
    R, q, eps = build_R(6, bonds23, t)
    D, R_int = integralize(R)
    check("exp(2K) = 5/3 and bond parity 1", q == Fraction(5, 3) and eps == 1,
          f"q = {q}, eps = {eps}")
    check(f"D = {D} clears every denominator (3^15 * 5^4)",
          D == 3**15 * 5**4, f"D = {D}")
    r = pipeline(R_int, p)
    check("charpoly monic degree 64 and Cayley--Hamilton Horner zero",
          r["chi_desc"][0] == 1 and len(r["chi_desc"]) == 65 and r["horner_zero"])
    spot_ok = True
    for k in (65, 100, 2016, 4032):
        spot_ok = spot_ok and mat_pow_trace(r["M"], k, p) == r["s"][k]
    check("trace recurrence matches binary-powering traces at k in "
          "{65, 100, 2016, 4032}", spot_ok)
    check("pair polynomial monic degree 2016",
          r["P_desc"][0] == 1 and len(r["P_desc"]) == 2017)
    check("explicit finite-field gcd divides both C_2 and C_2'",
          poly_div_rem(r["P_desc"], r["gcd"], p) == [0]
          and poly_div_rem(r["dP_desc"], r["gcd"], p) == [0])
    gdeg = r["gcd_degree"]
    check(f"2x3 layer gcd degree = {gdeg} < {GAUSSIAN_GCD_FLOOR} (DECISIVE)",
          gdeg < GAUSSIAN_GCD_FLOOR)
    distinct = PAIRS - gdeg
    check(f"inference chain: deg gcd_Q <= {gdeg} => distinct pair products "
          f">= {distinct} > {EXPONENT_VECTORS} => not a full 6-mode Gaussian spectrum",
          distinct > EXPONENT_VECTORS)

    # ------------------------------------------------------------- 6. 1D chain control
    print("\n6. 1D open chain n=6 (2D Ising strip) control at the same prime")
    Rc, _, _ = build_R(6, [(i, i + 1) for i in range(5)], t)
    Dc, Rc_int = integralize(Rc)
    rc = pipeline(Rc_int, p)
    check(f"chain control gcd degree = {rc['gcd_degree']} >= {GAUSSIAN_GCD_FLOOR} "
          f"[consistency only]", rc["gcd_degree"] >= GAUSSIAN_GCD_FLOOR,
          f"D_chain = {Dc}")

    # ------------------------------------------- 7. exact four-mode threshold (synthetic)
    print("\n7. synthetic four-mode Gaussian family: gcd degree exactly 120 - 65 = 55")
    u4 = [2, 3, 5, 7]
    vals4 = []
    for e in range(16):
        v = 1
        for i in range(4):
            if (e >> i) & 1:
                v *= u4[i]
        vals4.append(v)
    n_distinct = len({vals4[i] * vals4[j] for i in range(16) for j in range(i + 1, 16)})
    check("four-mode family has exactly 3^4 - 2^4 = 65 distinct pair products",
          n_distinct == 3**4 - 2**4, f"counted {n_distinct}")
    rng2 = random.Random(424242)
    dim4 = 16
    U = [[1 if a == b else 0 for b in range(dim4)] for a in range(dim4)]
    Uinv = [row[:] for row in U]
    for _ in range(6):
        i = rng2.randrange(dim4)
        j = rng2.randrange(dim4)
        while j == i:
            j = rng2.randrange(dim4)
        sc = rng2.choice((-2, -1, 1, 2))
        for rr in range(dim4):
            U[rr][j] += sc * U[rr][i]
        for cc in range(dim4):
            Uinv[i][cc] -= sc * Uinv[j][cc]
    G4 = [[sum(U[a][k] * vals4[k] * Uinv[k][b] for k in range(dim4)) for b in range(dim4)]
          for a in range(dim4)]
    rg = pipeline(G4, p)
    check("four-mode synthetic gcd degree is exactly 120 - 65 = 55",
          rg["gcd_degree"] == 120 - (3**4 - 2**4), f"got {rg['gcd_degree']}")

    # --------------------------------------------------------- 8. regression vs artifact
    print("\n8. regression against the producer artifact (content keys, not file bytes)")
    art = json.load(open(ROOT / "results" / "spectral" / "pair_product_obstruction.json"))
    data = art["data"]
    check("artifact prime p1 equals the prime recomputed here",
          data["primes"]["p1"] == p)
    lay = data["cases"]["layer_2x3"]["per_prime"][str(p)]
    check(f"2x3 gcd degree matches artifact ({lay['gcd_degree']})",
          gdeg == lay["gcd_degree"])
    check("2x3 charpoly digest matches artifact",
          digest([c % p for c in r["chi_desc"][::-1]]) == lay["charpoly_sha256"])
    check("2x3 pair-poly digest matches artifact",
          digest([c % p for c in r["P_desc"][::-1]]) == lay["pair_poly_sha256"])
    chc = data["cases"]["chain_n6_control"]["per_prime"][str(p)]
    check(f"chain control gcd degree matches artifact ({chc['gcd_degree']})",
          rc["gcd_degree"] == chc["gcd_degree"])
    check("artifact records the same D for the layer",
          int(data["integral_model"]["D_layer_2x3"]) == D)
    check("all producer checks passed in the stored artifact",
          all(c["passed"] for c in art["checks"]))

    print()
    if FAILS:
        print(f"FAIL: {FAILS}")
        return 1
    print(f"OK: independent rebuild confirms the ordering-free pair-product "
          f"certificate: 2x3 gcd degree {gdeg} < {GAUSSIAN_GCD_FLOOR} at p = {p}, "
          f"chain control {rc['gcd_degree']}, four-mode exact threshold 55.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
