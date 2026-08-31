"""gate_c_head.py — Gate C(a) head-side certified machinery (septic).

Certified (Arb, outward-rounded) head pipeline for the septic family:
  rho(t) = (t - s3^2 t^3 + s5^2 t^5 - s7^2 t^7)/V, V = 1+s3^2+s5^2+s7^2.
b_m chain (paper Eqs 44/45, lines 1723-1727):
  b_m = (pi/2) sum_{a+b<=m} (-1)^b A_{a,b}^2 [t^{m-b}] rho(t)^a.
The certified A-grid (grid_M251.json, a+b<=251, authors' frozen certificate) depends on
(eta, s3, s5) only through vartheta = eta/sqrt(V); the septic family changes V but NOT
the A_{a,b} themselves — the partitions f,g are built from vartheta*psi_3 alone
(paper lines 1717-1719, A_{a,b} = E[psi_b(X) q_a(vartheta psi_3(X))]). So the septic head
recomputes [t^{m-b}] rho^a EXACTLY (fmpq, septic coefficients) while importing the
certified A_{a,b} balls (CITED-DEPENDENCY), then accumulates b_m as a ball product.

Ball bookkeeping mirrors gate_A_final.compute_b_m (README-locked validated path):
  b_m ball = (pi/2) * (acc_mid ± acc_half) with
    acc_mid = sum term*c  (term = +-upper2 by (-1)^b sign convention, |.| taken per-term)
    acc_half = sum upper2 * |c| * 2*half_radius_propagation (crude outward)
Actually, to stay conservative AND match gate A: per-term upper bound of |A^2 c| used for
|b_m|: b_m_head_term = (pi/2) * upper2 * |c|. b_m's own ball (for the head midpoint) is
(pi/2)*(term*c ± slack). Head = sum |b_m| <= sum (|mid|+half)_outward  (ball abs).

b1 identity (README-locked, paper-validated at 4.86e-17):
  b1 = (pi/2) (A_{1,0}^2 / V - A_{0,1}^2)  with V the septic V (c1 = 1/V unchanged by
  the septic family!). Note: at septic parameters the identity holds for the same reason
  ([t^1]rho = c1 = 1/V exactly).
"""
import os, sys, math, time, json
from typing import Dict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import fmpq, arb
from core import _dec_fmpq, Prec, pi, exact_to_arb

ETA = _dec_fmpq("0.136419125")
PREC = 256
N_MAX = 251
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scratch", "repo")


def _dec(s: str) -> fmpq:
    return _dec_fmpq(s)


def rho_coeff_fmpq_septic(s3q: fmpq, s5q: fmpq, s7q: fmpq):
    """Exact septic rho coefficients {1,3,5,7} -> (V_q, coeff dict of FMPQ)."""
    V = fmpq(1) + s3q * s3q + s5q * s5q + s7q * s7q
    return V, {1: fmpq(1) / V, 3: -(s3q * s3q) / V, 5: (s5q * s5q) / V, 7: -(s7q * s7q) / V}


def load_grid():
    p = os.path.join(REPO, "grid_M251.json")
    with open(p) as f:
        raw = json.load(f)
    return raw["entries"]


def A_ball(entry):
    with Prec(PREC + 64):
        mid = arb(int(entry["m"])) * arb(10) ** entry["e"]
        rad = arb(int(entry["r"])) * arb(10) ** entry["e"]
    return mid, rad


def _poly_mul_q(p: Dict[int, fmpq], q: Dict[int, fmpq]) -> Dict[int, fmpq]:
    out: Dict[int, fmpq] = {}
    for i, ci in p.items():
        for j, cj in q.items():
            k = i + j
            v = out.get(k)
            out[k] = (v + ci * cj) if v is not None else ci * cj
    return out


def rho_pow_fmpq_septic(a: int, base: Dict[int, fmpq], cache: Dict):
    """rho^a as exact FMPQ dict — EXACT coefficients, converted to arb at the end."""
    if a == 0:
        return {0: fmpq(1)}
    if a in cache:
        return cache[a]
    cur = {0: fmpq(1)}
    bp = dict(base)
    e = a
    while e > 0:
        if e & 1:
            cur = _poly_mul_q(cur, bp)
        if e > 1:
            bp = _poly_mul_q(bp, bp)
        e >>= 1
    res = {j: q for j, q in cur.items() if q != 0}
    with Prec(PREC):
        cache[a] = {j: exact_to_arb(q) for j, q in res.items()}
    return cache[a]


def compute_head(s3q: fmpq, s5q: fmpq, s7q: fmpq, grid_entries, verbose=False):
    """Certified (head_up = sum_{3<=m<=251 odd} (|b_m|+width)_outward, plus b1 ball.
    Returns (b1_mid, b1_half, head_upper, meta)."""
    t0 = time.time()
    V_q, coeff_q = rho_coeff_fmpq_septic(s3q, s5q, s7q)
    cache: Dict[int, Dict[int, arb]] = {}
    by_a: Dict[int, list] = {}
    for entry in grid_entries:
        a, b = entry["a"], entry["b"]
        if a + b > 251:
            continue
        by_a.setdefault(a, []).append(entry)

    found10 = found01 = None
    for entry in by_a.get(1, []):
        if entry["b"] == 0:
            found10 = entry
    for entry in by_a.get(0, []):
        if entry["b"] == 1:
            found01 = entry
    if found10 is None or found01 is None:
        raise RuntimeError("grid missing (1,0) or (0,1) cell")

    pi2 = pi(PREC) / arb(2)
    head_upper = arb(0)
    with Prec(PREC):
        m10, r10 = A_ball(found10)
        m01, r01 = A_ball(found01)
        c1 = arb(fmpq(1) / V_q)
        a10_lo2 = (abs(m10) - r10).nonnegative_part() ** 2
        a10_hi2 = (abs(m10) + r10) ** 2
        a01_lo2 = (abs(m01) - r01).nonnegative_part() ** 2
        a01_hi2 = (abs(m01) + r01) ** 2
        # b1 = (pi/2)(A10^2 · c1 − A01^2): interval endpoints outward
        t_lo = a10_lo2 * c1.lower() - a01_hi2
        t_hi = a10_hi2 * c1.upper() - a01_lo2
        lo = pi2 * t_lo
        hi = pi2 * t_hi
        b1_mid = (lo + hi) / 2
        b1_half = (hi - lo) / 2

        # ---- head ----
        with Prec(PREC):
            for m in range(3, N_MAX + 1, 2):
                # per-m ball of b_m: mimic gate_A (signed mid + crude half)
                acc_mid = arb(0)
                acc_half = arb(0)
                for a, lst in by_a.items():
                    if a > m:
                        continue
                    rp = rho_pow_fmpq_septic(a, coeff_q, cache)
                    for entry in lst:
                        b = entry["b"]
                        if a + b > m:
                            continue
                        c = rp.get(m - b)
                        if c is None:
                            continue
                        mm, rr = A_ball(entry)
                        upper2 = (abs(mm) + rr) ** 2
                        term = (-upper2) if b % 2 else upper2
                        acc_mid += term * c
                        acc_half += upper2 * abs(c) * abs(rr) * 2  # gate A's crude outward
                bm = pi2 * acc_mid
                bh = pi2 * abs(acc_half)
                head_upper += abs(bm) + bh
                if verbose and m % 40 == 1:
                    print(f"  m={m}: |b_m|+w up to {head_upper.str(12)} ({time.time()-t0:.0f}s)")
    meta = dict(seconds=time.time() - t0, V_str=str(V_q),
                c_str={str(k): str(v) for k, v in coeff_q.items()})
    return b1_mid, b1_half, head_upper, meta


def gamma_star_head(b1_mid, b1_half, head_upper):
    """gamma*_head = b1_lower - head_upper (head-only, no tail term). CONDITIONAL objective."""
    return b1_mid - b1_half - head_upper


def gamma_from_s(s3: float, s5: float, s7: float, grid_entries):
    s3q, s5q, s7q = _dec(repr(s3)), _dec(repr(s5)), _dec(repr(s7))
    b1_mid, b1_half, head_up, meta = compute_head(s3q, s5q, s7q, grid_entries)
    g = gamma_star_head(b1_mid, b1_half, head_up)
    return g, dict(b1_mid=b1_mid.str(20), b1_half=b1_half.str(6),
                   head_up=head_up.str(16), **meta)


if __name__ == "__main__":
    ents = load_grid()
    print(f"grid entries: {len(ents)}  (a+b<=251 filter applied inside)")
    t0 = time.time()
    # Anchor: s7=0 → should match Gate A exactly: gamma*_head ≈ b1_lo - 1.1328860e-5
    g, info = gamma_from_s(0.34101124, 0.05276111, 0.0, ents)
    print(f"anchor s7=0: gamma*_head = {g.str(30)}  [{time.time()-t0:.0f}s]")
    print("  b1 =", info["b1_mid"], "±", info["b1_half"], "  head_up =", info["head_up"])
