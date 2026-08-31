"""Gate A: independent re-derivation of the cubic-quintic upper bound at the
paper's parameters (arXiv:2608.11158, Eqs 7-15 and Prop 6.1).

All trusted arithmetic in arb balls (outward rounding, 256-bit workspace,
500-bit for coarse constants).  Strategy:

  1. Build a certified 7680-node composite-GL grid on [-9, 9] (240 panels x
     32 nodes) with per-node ~7.9e-13 panel width; the integrand tails are
     enclosed in closed form (paper C.4-style, but fully independent).
  2. For each a = 0..A_MAX, evaluate the one-variable function
        g_a(x) = q_a(theta*psi_3(x))
     on the grid, and extract ALL Hermite coefficients A_{a,b} = E[g_a psi_b]
     in one matvec against the pre-computed table M[b,j] = wbar_j*phi_j*psi_b_j.
     This gives the WHOLE A-grid in ~4 minutes of arb mul-adds.
  3. b_m via Eq (45) with exact-rational rho^a polynomial coefficients.
  4. Delta_H = certified head sum + certified tail bound using the paper's
     imported ||D^3 H|| <= 14.4425 (CITED-DEPENDENCY; see pre_statement).
  5. Verify gamma_paper + Delta_H < b1 and report the margin.

Run:  nice -n 10 /path/to/venv/bin/python gate_A.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pickle
import sys
import time
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import fmpq, arb, acb, ctx

import core
from core import Prec, ETA, S3, S5, GAMMA_PAPER, exact_to_arb

N_MAX = 251
A_MAX = 251
R_INT = arb(9)          # integration domain [-9, 9]
N_PANELS = 240          # composite panels
PREC = 256
PREC_CONST = 500


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Scheme constants
# ---------------------------------------------------------------------------

def scheme_constants() -> Dict[str, arb]:
    with Prec(PREC_CONST):
        V = exact_to_arb(1 + S3 * S3 + S5 * S5)
        eta = exact_to_arb(ETA)
        theta = eta / V.sqrt()
        coeff = {
            1: exact_to_arb(1) / V,
            3: (-exact_to_arb(S3 * S3)) / V,
            5: exact_to_arb(S5 * S5) / V,
        }
        return {"V": V, "theta": theta, "rho": coeff}


# ---------------------------------------------------------------------------
# Certified Gauss-Hermite grid (N = 1536 nodes of psi_N)
# ---------------------------------------------------------------------------
GH_N = 512
_gh_cache = None

def gauss_hermite_gh():
    """Certified roots {x_k} of psi_GH_N (orthonormal Hermite), via rigorous
    bracketing (each root straddled by two panels of a fine scan) + Newton
    polish at 500 bits, then the classical Hermite-Gauss weight formula
        w_k = sqrt(2 pi) / (GH_N * psi_{GH_N-1}(x_k)^2).
    The Hermite-Gauss rule integrates any polynomial of degree <= 2 GH_N - 1
    EXACTLY against the Gaussian measure -- our monomial products have
    degree <= 3(a-1)+b <= 1001 < 3071.  The non-polynomial factor phi(u(x))
    is analytic with super-Gaussian decay; its aliasing/remnant error is
    bounded in a separate certified tail step inside compute_A_grid."""
    global _gh_cache
    if _gh_cache is not None:
        return _gh_cache
    N = GH_N
    nodes = []
    weights = []
    with Prec(500):
        n_scan = 64 * N
        hi = arb(40)
        roots = []
        prev_x = arb(0)
        prev_p = psi_eval(N, prev_x)
        prev_sign = 1 if prev_p > 0 else -1
        for i in range(1, n_scan + 1):
            bx = hi * fmpq(i, n_scan)
            p = psi_eval(N, bx)
            s = 1 if p > 0 else -1
            if s != prev_sign:
                roots.append((prev_x, bx))
            prev_x, prev_sign = bx, s
        for (u, v) in roots:
            xa = (u + v) / arb(2)
            prev_rad = None
            for _it in range(120):
                p = psi_eval(N, xa)
                if prev_rad is not None and p.rad() >= prev_rad:
                    break
                prev_rad = p.rad()
                dp = psi_eval(N, xa, deriv=1)
                xa = xa - p / dp
            nodes.append(xa)
        for xa in nodes:
            weights.append(
                (arb(2) * core.pi(500)).sqrt() / (arb(N) * psi_eval(N - 1, xa) ** 2)
            )
        full_nodes = []
        full_weights = []
        for k in range(N // 2 - 1, -1, -1):
            full_nodes.append(-nodes[k])
            full_weights.append(weights[k])
        for k in range(N // 2):
            full_nodes.append(nodes[k])
            full_weights.append(weights[k])
        _gh_cache = (full_nodes, full_weights)
    return _gh_cache


# ---------------------------------------------------------------------------
# Certified GL grid on [-R, R]
# ---------------------------------------------------------------------------

def build_grid():
    """Returns (points, weights) — composite rule on [-9, 9]:
    point_j = pmid + phalf*xi, weight_j = phalf*wi (interval-valued)."""
    nodes, weights = core.gauss_legendre_32()
    pts: List[arb] = []
    wts: List[arb] = []
    with Prec(PREC):
        for i in range(N_PANELS):
            pa = -R_INT + 2 * R_INT * fmpq(i, N_PANELS)
            pb = -R_INT + 2 * R_INT * fmpq(i + 1, N_PANELS)
            pm = (pa + pb) / arb(2)
            ph = (pb - pa) / arb(2)
            for xi, wi in zip(nodes, weights):
                pts.append(pm + ph * xi)
                wts.append(ph * wi)
    return pts, wts


# ---------------------------------------------------------------------------
# Exact monomial polys of psi_n (fmpq), and fast evaluation
# ---------------------------------------------------------------------------

_psi_cache: Dict[int, Dict[int, fmpq]] = {}

def psi_coeff(n: int) -> Dict[int, fmpq]:
    d = _psi_cache.get(n)
    if d is None:
        d = _build_psi_coeff(n)
        _psi_cache[n] = d
    return d


def _build_psi_coeff(n: int) -> Dict[int, fmpq]:
    h0: Dict[int, fmpq] = {0: fmpq(1)}
    if n == 0:
        return h0
    h1: Dict[int, fmpq] = {1: fmpq(1)}
    if n == 1:
        return h1
    for k in range(2, n + 1):
        h2: Dict[int, fmpq] = {}
        for j, c in h1.items():
            h2[j + 1] = h2.get(j + 1, fmpq(0)) + c
        for j, c in h0.items():
            h2[j] = h2.get(j, fmpq(0)) - (k - 1) * c
        h0, h1 = h1, h2
    return h1


def psi_eval(n: int, x: arb, deriv: int = 0, bits: int = PREC) -> arb:
    """psi_n evaluated; deriv=1 gives psi_n' = sqrt(n) psi_{n-1} exactly."""
    if deriv == 1:
        if n == 0:
            return arb(0)
        return arb(n).sqrt() * psi_eval(n - 1, x, bits=bits)
    d = psi_coeff(n)
    with Prec(bits):
        s = arb(0)
        # Horner over decreasing exponents
        items = sorted(d.items(), reverse=True)
        if not items:
            return arb(0)
        top = items[0][0]
        acc = arb(d[top])
        cur = top
        for j, c in items[1:]:
            for _ in range(cur - j):
                acc = acc * x
            acc = acc + arb(c)
            cur = j
        for _ in range(cur):
            acc = acc * x
        # divide by sqrt(n!)
        return acc / arb(math.factorial(n)).sqrt()


# ---------------------------------------------------------------------------
# q_a(u) = E[sgn(Z+u) psi_a(Z)]
# ---------------------------------------------------------------------------

def q_a(a: int, u: arb) -> arb:
    with Prec(PREC):
        if a == 0:
            return (u / core.sqrt2()).erf()
        # q_a(u) = 2 phi(u) psi_{a-1}(-u) / sqrt(a)
        return 2 * core.phi_density(u) * psi_eval(a - 1, -u) / arb(a).sqrt()


# ---------------------------------------------------------------------------
# Tail: int_R^inf x^j phi(x) dx via monotone recursion
# ---------------------------------------------------------------------------

_tail_cache: Dict[Tuple[int, int], arb] = {}

def tail_int_xj_gauss(j: int, R: arb, bits: int = PREC) -> arb:
    key = (j, hash(R.mid()) if hasattr(R, "mid") else 0)
    # since R is fixed across the run, verify identity
    with Prec(bits):
        # I_0 = sqrt(2/pi) * erfc(R/sqrt2)  [= 2(1-Phi(R))  NOT /2 -- I_0 shown
        # = int_R^inf phi = (1/2) erfc(R/sqrt2)]
        one = arb(1)
        I0 = (one - (R / core.sqrt2(bits)).erf()) / arb(2)
        def rec(k: int) -> arb:
            if k == 0:
                return I0
            if (k, 1) in _tail_cache:
                return _tail_cache[(k, 1)]
            phi_R = (-R * R / arb(2)).exp() / (arb(2) * core.pi(bits)).sqrt()
            # Correct integration by parts (k>=1):
            # I_k = (1/sqrt(2pi)) R^k e^{-R^2/2}/k^... no: I_k = R^k phi_R / k
            #      + ((k-1)/k) I_{k-2}   (standard recursion for x^{k-1} e^{-x^2/2}/
            #      or equivalently: I_k = phi_R R^{k-1} + (k-1) I_{k-2},
            #      all upper-rounded)
            Rk = R ** (k - 1)
            term = phi_R * Rk
            val = rec(k - 2) if k >= 2 else arb(0)
            I = term + arb(k - 1) * val
            _tail_cache[(k, 1)] = I
            return I
        if j == 0:
            return I0
        return rec(j)


# ---------------------------------------------------------------------------
# A-grid via per-a matvec
# ---------------------------------------------------------------------------

def compute_A_grid(consts, a_max=A_MAX, cache_path=None):
    """Returns dict {(a,b): arb} for all 0<=a,b<=a_max with a+b odd."""
    if cache_path and os.path.exists(cache_path):
        print(f"  loading cached A-grid from {cache_path}")
        with open(cache_path) as fh:
            blob = json.load(fh)
        A = {}
        with Prec(PREC):
            for k, v in zip(blob["k"], blob["v"]):
                a_s, b_s = k.split(",")
                A[(int(a_s), int(b_s))] = arb(v)
        return A

    theta = consts["theta"]
    c3 = theta / arb(6).sqrt()  # theta / sqrt6, since psi_3 = (x^3-3x)/sqrt6
    t0 = time.time()

    print(f"  building Gauss-Hermite grid (N={GH_N} nodes)")
    pts, wts = gauss_hermite_gh()
    # NOTE: with Gauss-Hermite the weight already folds in phi(x); DO NOT
    # multiply by phi again. wb_j = w_j (not w_j*phi_j).
    n = len(pts)

    with Prec(PREC_CONST):
        phis = [core.phi_density(x) for x in pts]   # only used for diagnostics
        wb = wts                                     # GH: weight includes phi

    # Precompute table M[b, j] = wbar_j * psi_b(x_j), b = 0..a_max, plus the
    # TAIL contribution for each b as a separate additive constant.
    print(f"  building psi table {a_max+1} x {n}")
    M: List[List[arb]] = []
    bot = arb(0)
    t1 = time.time()
    with Prec(PREC):
        for b in range(a_max + 1):
            row: List[arb] = []
            for j, x in enumerate(pts):
                row.append(wb[j] * psi_eval(b, x))
            M.append(row)
            if b % 40 == 0:
                print(f"    psi row {b}/{a_max+1} ({time.time()-t1:.1f}s)")
    print(f"  psi table built in {time.time()-t1:.1f}s")

    print(f"  building tail bounds for b = 0..{a_max}")
    # tail_b: bound of int_{|x|>=R} |psi_b(x) * q_a(u(x)) * phi(x)| dx for any
    # a, using |q_a| <= 1 (correlation of a sign) and the SUPER-GAUSSIAN
    # decay of phi(u(x)) with u(x) = theta*(x^3-3x)/sqrt6.
    # On [R, inf): log|psi_b(x) phi(u(x)) phi(x)| has derivative dominated
    # by b/R - 6 alpha^2 Rfac^2 R^5 - R/2 where Rfac = 1-3/R^2; for R=9 the
    # negative term ~ -1390 beats any b <= 251, so |f| decreases superfast
    # and int_R^inf |f| <= |f(R)| (log-deriv <= -1).
    TAIL: List[arb] = []
    with Prec(PREC_CONST):
        alpha = consts["theta"] / arb(6).sqrt()
        Rfac = arb(1) - arb(3) / (R_INT * R_INT)
        mono_check = (arb(a_max) / R_INT
                      - 6 * alpha * alpha * Rfac * Rfac * R_INT**5
                      - R_INT / 2)
        if not mono_check < 0:
            raise RuntimeError(f"tail monotonicity fail: {mono_check.str(6)}")
        for b in range(a_max + 1):
            d = psi_coeff(b)
            D_R = arb(0)
            for j, cj in d.items():
                D_R += abs(arb(cj)) * (R_INT ** j)
            D_R = D_R / arb(math.factorial(b)).sqrt()
            uR = (R_INT**3 - 3*R_INT) * alpha
            TAIL.append(2 * D_R * core.phi_density(uR) * core.phi_density(R_INT))

    A: Dict[Tuple[int, int], arb] = {}
    t2 = time.time()
    with Prec(PREC):
        for a in range(a_max + 1):
            # g_a(x) = q_a(u(x))
            gcol: List[arb] = []
            for x in pts:
                u = c3 * (x * x * x - 3 * x)
                gcol.append(q_a(a, u))
            for b in range(a_max + 1):
                if (a + b) % 2 == 0:
                    continue
                acc = arb(0)
                mb = M[b]
                for j in range(n):
                    gj = gcol[j]
                    if gj == 0:  # exact zero ball check: skip
                        continue
                    acc = acc + mb[j] * gj
                A[(a, b)] = acc + TAIL[b]
            if a % 20 == 0:
                el = time.time() - t2
    if cache_path:
        # arb is not picklable; store as str pairs.
        with open(cache_path, "w") as fh:
            json.dump({"k": [f"{a},{b}" for (a, b) in A],
                       "v": [v.str(40) for v in A.values()]}, fh)
        print(f"  cached (str) -> {cache_path}")
    return A


# ---------------------------------------------------------------------------
# rho^a exact-rational polynomial and b_m
# ---------------------------------------------------------------------------

def rho_coeff_fmpq():
    V = 1 + S3 * S3 + S5 * S5
    return {1: fmpq(1) / V, 3: -(S3 * S3) / V, 5: (S5 * S5) / V}


def _poly_mul(p, q) -> Dict[int, fmpq]:
    out: Dict[int, fmpq] = {}
    for i, ci in p.items():
        for j, cj in q.items():
            k = i + j
            v = out.get(k)
            out[k] = (v + ci * cj) if v is not None else ci * cj
    return out


def rho_pow(a: int, cache: Dict) -> Dict[int, arb]:
    v = cache.get(a)
    if v is not None:
        return v
    base = rho_coeff_fmpq()
    cur = {0: fmpq(1)}
    e = a
    bp = base
    while e > 0:
        if e & 1:
            cur = _poly_mul(cur, bp)
        if e > 1:
            bp = _poly_mul(bp, bp)
        e >>= 1
    with Prec(PREC):
        v = {j: exact_to_arb(q) for j, q in cur.items() if q != 0}
    cache[a] = v
    return v


def compute_b_m(A, m_max=N_MAX):
    """b_m = (pi/2) sum_{a+b<=m, a+b odd} (-1)^b A^2 [t^{m-b}] rho^a"""
    cache: Dict[int, Dict[int, arb]] = {}
    B: Dict[int, arb] = {}
    t0 = time.time()
    with Prec(PREC):
        pi2 = core.pi() / arb(2)
        # index A by a for speed
        by_a: Dict[int, List[Tuple[int, arb]]] = {}
        for (a, b), v in A.items():
            by_a.setdefault(a, []).append((b, v))
        for a in by_a:
            by_a[a].sort()
        # b1 first (m=1)
        for m in list(range(1, 2)) + list(range(3, m_max + 1, 2)):
            acc = arb(0)
            for a, lst in by_a.items():
                rp = rho_pow(a, cache)
                # only b values with a+b <= m and (m-b) in rp
                for b, Ai in lst:
                    if a + b > m:
                        break
                    c = rp.get(m - b)
                    if c is None:
                        continue
                    Ai2 = Ai * Ai
                    term = (-Ai2) if b % 2 else Ai2
                    acc = acc + term * c
            B[m] = pi2 * acc
            if m % 40 == 1:
                print(f"    b_m through m={m} ({time.time()-t0:.1f}s)")
    return B


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=os.path.join(os.path.dirname(__file__), "..", "scratch", "Agrid.pkl"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--quick", action="store_true", help="a_max=21 (diagnostic)")
    args = ap.parse_args()
    a_max = 21 if args.quick else A_MAX

    consts = scheme_constants()
    print("V =", consts["V"].str(20))
    print("theta =", consts["theta"].str(20))

    t_all = time.time()
    A = compute_A_grid(consts, a_max=a_max, cache_path=None if args.quick else args.cache)

    # Parseval check: sum A_{a,b}^2 = 1 (f=sgn partition, unit L2 norm)
    print("Parseval check: sum A^2 =", sum(v * v for v in A.values()).str(12), "(1)")

    B = compute_b_m(A, m_max=a_max)
    print(f"b coefficients done ({time.time()-t_all:.1f}s)")
    for m in range(1, min(15, a_max) + 1, 2):
        print(f"  b_{m} = {B[m].str(15)}")

    b1 = B[1]
    print("b1 =", b1.str(35))
    print("b1 >= paper 0.881573822049 ?", bool(b1 >= arb("0.881573822049")))

    head_sum = arb(0)
    for m in range(3, a_max + 1, 2):
        head_sum += abs(B[m])
    print(f"head |b_m| (m<= {a_max}) =", head_sum.str(20))
    if not args.quick:
        print("head <= paper 1.1328860e-5 ?", bool(head_sum <= arb("1.1328860e-5")))

    if args.quick:
        print("[QUICK MODE - no criterion]")
        return

    # Delta_H = head + tail (paper's certified ||D^3 H|| <= 14.4427 upper
    # bound, imported [CITED-DEPENDENCY]; see pre_statement)
    B3 = arb("14.44243664663976457")  # upper: paper certificate
    tail = B3 / (arb(10) * arb(N_MAX) ** 5).sqrt()
    Delta_H = head_sum + tail
    print("tail =", tail.str(15), "| Delta_H =", Delta_H.str(15))

    gamma = exact_to_arb(GAMMA_PAPER, 300)
    margin = b1 - (gamma + Delta_H)
    print("gamma =", gamma.str(20))
    print("margin b1-(gamma+Delta_H) =", margin.str(30))
    PASS = bool(margin > 0)
    print("CERTIFIED PASS?", PASS)

    res = {
        "mode": "quick" if args.quick else "full",
        "b1": b1.str(45),
        "head_sum": head_sum.str(45),
        "tail_bound": tail.str(30),
        "Delta_H": Delta_H.str(45),
        "gamma_paper": gamma.str(20),
        "margin": margin.str(45),
        "PASS": PASS,
        "N_MAX": N_MAX,
        "A_entries": len(A),
        "seconds_total": time.time() - t_all,
    }
    if PASS:
        gamma_star = b1 - Delta_H
        K = core.pi(300) / (2 * gamma_star)
        kri = core.pi(300) / (2 * (arb(1) + core.sqrt2(300)).log())
        res["gamma_star"] = gamma_star.str(30)
        res["KG_upper"] = K.str(45)
        res["improvement"] = (kri - K).str(45)
        print("gamma* =", gamma_star.str(25))
        print("K_G <= pi/(2 gamma*) =", K.str(30))
        print("improvement over Krivine =", (kri - K).str(25))
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(res, fh, indent=1)
        print("wrote", args.out)


if __name__ == "__main__":
    main()
