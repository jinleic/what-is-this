"""gate_c_fast.py — Gate C(a): FAST certified head machinery via arb_poly convolution,
plus the certified box bound used by the branch-and-bound.

Identity used (paper Eq (44)/(45), lines 1723-1727):
    H(t) = (pi/2) sum_{a,b>=0} A_{a,b}^2 rho(t)^a (-t)^b
  ⇒ sum_m b_m t^m = (pi/2) sum_a [ Q_a(t) · rho(t)^a ],   Q_a(t) := sum_b (-1)^b A_{a,b}^2 t^b.
So the whole coefficient head is ONE truncated polynomial convolution per a — done in
arb_poly (C level, ball coefficients), not Python loops. Truncation at degree N=251
matches the paper's head cutoff (line 1895).

Certified per-BOX evaluation: rho's coefficients are intervals over the box
  c1 = 1/V,  c3 = -s3^2/V,  c5 = +s5^2/V,  c7 = -s7^2/V,   V = 1+s3^2+s5^2+s7^2
with s_d ranging over the box; each coefficient is enclosed as an arb ball built from
exact rational endpoints (mid/rad form — NEVER arb(lo,hi), which is mid+radius).

Objective and box bound:
  gamma*(s) = b1(s) - head(s) - tail(s)  with tail(s) >= 0 and head(s) >= 0,
  b1(s) = (pi/2)(A_{1,0}^2 / V(s) - A_{0,1}^2)  is DECREASING in V,
  hence for every box B:   gamma*(s) <= b1(V_min(B)) - head_lo(B)   for all s in B.
That upper bound is what the branch-and-bound prunes on; it needs no tail input, so a
certified CAP is reachable without any S(theta) integral. (An IMPROVEMENT claim does
need the tail, which is the float sweep — labelled COMPUTATIONAL-EVIDENCE.)
"""
import os, sys, json, time
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import fmpq, arb, arb_poly, ctx
from core import _dec_fmpq, Prec, pi

PREC = 256
N_MAX = 251
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scratch", "repo")
ETA = _dec_fmpq("0.136419125")


# --------------------------------------------------------------------------- grid
def load_grid_entries():
    with open(os.path.join(REPO, "grid_M251.json")) as f:
        return json.load(f)["entries"]


def build_Q_polys(entries, prec=PREC):
    """Q_a(t) = sum_b (-1)^b A_{a,b}^2 t^b as arb_poly with certified ball coefficients.
    A^2 ball: [max(0,|mid|-rad)^2, (|mid|+rad)^2] -> (midpoint, radius) form."""
    ctx.prec = prec + 64
    by_a: Dict[int, Dict[int, arb]] = {}
    a10 = a01 = None
    for e in entries:
        a, b = e["a"], e["b"]
        if a + b > N_MAX:
            continue
        mid = arb(int(e["m"])) * arb(10) ** e["e"]
        rad = arb(int(e["r"])) * arb(10) ** e["e"]
        lo = (abs(mid) - rad).nonnegative_part() ** 2
        hi = (abs(mid) + rad) ** 2
        c_mid = (lo + hi) / 2
        c_rad = (hi - lo) / 2
        A2 = arb(c_mid, c_rad)          # midpoint/radius — correct ball construction
        if a == 1 and b == 0:
            a10 = A2
        if a == 0 and b == 1:
            a01 = A2
        by_a.setdefault(a, {})[b] = (-A2) if (b % 2) else A2
    Qs = {}
    for a, d in by_a.items():
        bmax = max(d)
        coeffs = [arb(0)] * (bmax + 1)
        for b, v in d.items():
            coeffs[b] = v
        Qs[a] = arb_poly(coeffs)
    ctx.prec = prec
    return Qs, a10, a01


# --------------------------------------------------------- septic rho over a box
def _exact_q(v) -> fmpq:
    """Exact rational from a decimal string OR a Python float (binary value, exact).
    Never approximates: Fraction(float) is the exact binary rational."""
    if isinstance(v, fmpq):
        return v
    if isinstance(v, int):
        return fmpq(v)
    if isinstance(v, float):
        from fractions import Fraction
        fr = Fraction(v)
        return fmpq(fr.numerator, fr.denominator)
    s = str(v)
    if "e" in s or "E" in s:
        from fractions import Fraction
        fr = Fraction(float(s))
        return fmpq(fr.numerator, fr.denominator)
    return _dec_fmpq(s)


def rho_poly_box(s3lo, s3hi, s5lo, s5hi, s7lo, s7hi, prec=PREC):
    """rho(t) coefficients as certified balls over the box [s3lo,s3hi]x[s5lo,s5hi]x[s7lo,s7hi].
    Endpoints exact rationals; balls built mid/rad at working precision."""
    ctx.prec = prec + 64
    q3lo, q3hi = _exact_q(s3lo), _exact_q(s3hi)
    q5lo, q5hi = _exact_q(s5lo), _exact_q(s5hi)
    q7lo, q7hi = _exact_q(s7lo), _exact_q(s7hi)
    u3lo, u3hi = q3lo * q3lo, q3hi * q3hi           # s3^2 range (s>=0 so monotone)
    u5lo, u5hi = q5lo * q5lo, q5hi * q5hi
    u7lo, u7hi = q7lo * q7lo, q7hi * q7hi
    Vlo = fmpq(1) + u3lo + u5lo + u7lo
    Vhi = fmpq(1) + u3hi + u5hi + u7hi

    def ball_from_q(lo: fmpq, hi: fmpq) -> arb:
        mid = arb(lo + hi) / 2
        rad = arb(hi - lo) / 2
        return arb(mid, rad)

    # c1 = 1/V  in [1/Vhi, 1/Vlo]
    c1 = ball_from_q(fmpq(1) / Vhi, fmpq(1) / Vlo)
    # c3 = -u3/V in [-u3hi/Vlo, -u3lo/Vhi]
    c3 = ball_from_q(-(u3hi / Vlo), -(u3lo / Vhi))
    # c5 = +u5/V in [u5lo/Vhi, u5hi/Vlo]
    c5 = ball_from_q(u5lo / Vhi, u5hi / Vlo)
    # c7 = -u7/V in [-u7hi/Vlo, -u7lo/Vhi]
    c7 = ball_from_q(-(u7hi / Vlo), -(u7lo / Vhi))
    ctx.prec = prec
    coeffs = [arb(0)] * 8
    coeffs[1], coeffs[3], coeffs[5], coeffs[7] = c1, c3, c5, c7
    return arb_poly(coeffs), Vlo, Vhi


# ------------------------------------------------------------------ head & gamma
def head_and_b1(rho_p: arb_poly, Vlo: fmpq, Vhi: fmpq, Qs, a10, a01,
                prec=PREC, nmax=N_MAX):
    """Certified (b1_lo, head_lo, head_hi) for the given rho-interval polynomial.
    head = sum_{3<=m<=nmax odd} |b_m|; b1 from the two-cell identity at V = Vhi (lower bound
    on b1 over the box, since b1 decreases in V) and V = Vlo (upper bound)."""
    ctx.prec = prec
    pi2 = pi(prec) / arb(2)
    # b1 bounds: b1(V) = (pi/2)(A10^2/V - A01^2)
    b1_lo = pi2 * (a10 / arb(Vhi) - a01)
    b1_hi = pi2 * (a10 / arb(Vlo) - a01)
    # coefficient accumulation:   T(t) = sum_a Q_a(t) * rho(t)^a   truncated at nmax
    T = [arb(0)] * (nmax + 1)
    amax = max(Qs)
    Pa = arb_poly([arb(1)])              # rho^0
    for a in range(0, amax + 1):
        Q = Qs.get(a)
        if Q is not None:
            prod = Q * Pa
            for k in range(min(len(prod), nmax + 1)):
                T[k] += prod[k]
        if a < amax:
            Pa = Pa * rho_p
            if len(Pa) > nmax + 1:       # truncate
                Pa = arb_poly([Pa[k] for k in range(nmax + 1)])
    head_lo = arb(0)
    head_hi = arb(0)
    for m in range(3, nmax + 1, 2):
        bm = pi2 * T[m]
        amid = abs(bm)
        head_hi += amid.upper()
        lo = amid.lower()
        if lo > 0:
            head_lo += lo
    return b1_lo, b1_hi, head_lo, head_hi


def gamma_box_upper(box, Qs, a10, a01, prec=PREC, nmax=N_MAX):
    """Certified UPPER bound on gamma*(s) = b1 - head - tail over the box (tail >= 0)."""
    (s3lo, s3hi), (s5lo, s5hi), (s7lo, s7hi) = box
    rho_p, Vlo, Vhi = rho_poly_box(s3lo, s3hi, s5lo, s5hi, s7lo, s7hi, prec)
    b1_lo, b1_hi, head_lo, head_hi = head_and_b1(rho_p, Vlo, Vhi, Qs, a10, a01, prec, nmax)
    ub = b1_hi.upper() - head_lo          # gamma* <= b1_max - head_min - 0
    return dict(ub=ub, b1_lo=b1_lo, b1_hi=b1_hi, head_lo=head_lo, head_hi=head_hi,
                Vlo=Vlo, Vhi=Vhi)


def gamma_point(s3: str, s5: str, s7: str, Qs, a10, a01, prec=PREC, nmax=N_MAX):
    """Certified head-only gamma at a single parameter point (degenerate box)."""
    return gamma_box_upper(((s3, s3), (s5, s5), (s7, s7)), Qs, a10, a01, prec, nmax)


if __name__ == "__main__":
    t0 = time.time()
    ents = load_grid_entries()
    Qs, a10, a01 = build_Q_polys(ents)
    print(f"Q polys: {len(Qs)} slices, build {time.time()-t0:.1f}s")
    print(f"  A_(1,0)^2 = {a10.str(15)}   A_(0,1)^2 = {a01.str(15)}")

    # ---- ANCHOR at the paper point (s7 = 0) ----
    t1 = time.time()
    res = gamma_point("0.34101124", "0.05276111", "0", Qs, a10, a01)
    dt = time.time() - t1
    print(f"\nANCHOR (s7=0): [{dt:.2f}s]")
    print(f"  b1        = {res['b1_lo'].str(25)}   (paper >= 0.881573822049)")
    print(f"  head      = [{res['head_lo'].str(16)}, {res['head_hi'].str(16)}]"
          f"   (paper <= 1.1328860e-5)")
    print(f"  gamma_head UB = {res['ub'].str(20)}")
    tail_paper = arb("14.2459830173") / (arb(10) * arb(251) ** 5).sqrt()
    g_full = res['b1_lo'] - res['head_hi'] - tail_paper
    print(f"  with our corrected B3 = 14.2459830173: tail = {tail_paper.str(12)}")
    print(f"  gamma* (head+tail)   = {g_full.str(20)}   (gate A: 0.881557917504162)")
    K = pi(300) / (2 * g_full)
    kri = pi(300) / (2 * (arb(1) + arb(2).sqrt()).log())
    print(f"  K_G <= {K.str(22)}   improvement = {(kri - K).str(12)}")
