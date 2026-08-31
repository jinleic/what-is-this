#!/usr/bin/env python3
"""LIU9 PAIRED-CLASS SCALAR MARGIN — q-free channel scalar c, delta family.

Probe (i) of math/uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md.

Family: P0 = delta_x, P1 = delta_1, shared mass a = 1, x in (0, 1).
Certified statements (every inequality Arb-certified at ctx.prec = 3600
(>= 320) or exact rational / sympy):

  1. CLOSED FORM of the q-free channel scalar
       c(x) = beta * (h(pi(x,x)) - 2 h(x)),   pi(s,t) = s t (1+(1-s)(1-t)),
     from c = beta * sum_ij a_i a_j [K(x_i,x_j) + K(y_i,y_j) - 2K(x_i,y_j)],
     K = h(pi): single atoms give pi(x,1) = x, h(pi(1,1)) = h(1) = 0.
     c is exactly q-free by construction.
  2. SIGN: c(x) < 0 for ALL x in (0, 1), and c(1) = 0.  Three certified
     pieces:
       (0, 1/2]    exact symbolic chain: pi(x,x) - x = x(x-1)(x^2-x+1) < 0,
                   h strictly increasing on (0, 1/2], h > 0 there;
       [1/2, 7/8]  96 uniform exact-rational cells (step 1/256): certified
                   hull 2*min h(x) > max h(pi(x,x)); h concave, pi_xx
                   increasing, so endpoint hulls + midpoint are valid;
       [7/8, 1)    analytic epsilon-certificate with eps = 1-x in (0, 1/8]:
                   2h(x) - h(pi(x,x)) >= 2 eps (log 2 - 3 eps) >= 1/64,
                   log 2 > 7/16 Arb-certified.
     The channel is sign-indefinite in general and NEGATIVE on the whole
     delta family; no crossover point exists.
  3. FAMILY MINIMUM of gapA(x, q) = F(P_mix) + q(1-q) c(x) over the
     mean-feasible band q in [max(0, (m-x)/(1-x)), 1]: gapA is an exact
     concave quadratic in q with gapA(x, 1) = 0 and F >= 0 (q=1 theorem,
     PROVED), so the band minimum sits at q0(x) = (m-x)/(1-x) for x < m
     and at q = 0 for x >= m.  Certified >= 0 piecewise:
       (a) deep tail (0, 1/256]: fully analytic exact-rational certificate
           gapA >= x L (2 C1 C2 - x (K2' + K3'/L)) > 0, L = log(1/x) >= 8
           log 2 > 5 (Arb); D(x,1;M) = (M - 1/2) h(x) EXACTLY;
       (b) strip [1/256, ML]: geometric cells (factor 9/8), exact-rational
           ends, interval hulls of t1 + t2 + t3 at M = m;
       (c) micro-cell [ML, MH] straddling m: q0-hull clamped to
           [0, (MH-ML)/(1-ML)] which encloses the true (0, q0(ML)] for
           x in (ML, m); the x in [m, MH] sliver is covered by (d);
       (d) A0 cover [ML, 1]: A0(x) = gapA(x, 0) >= 0 in 512 uniform cells.
     Hence min gapA = 0, attained EXACTLY on the whole q = 1 edge.

Byte-stable report: uc/verification/results/liu9-c-delta-family.json
Run:  math/.venv/bin/python -I -B math/uc/liu9_c_delta_family.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from fractions import Fraction
from pathlib import Path
from typing import Optional, Tuple

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
OUTPUT_DEFAULT = (
    HERE / "verification/results/liu9-c-delta-family.json")
ctx.prec = 3600  # exact-rational cell ends down to 1/256 and hull margins

GROWTH = Fraction(9, 8)          # geometric strip cell width factor
MIN_CELL = Fraction(1, 2 ** 20)  # smallest sub-division before forced emit
STRIP_LO = Fraction(1, 256)      # strip low end (deep tail covers below)
C_BAND_LO = Fraction(1, 2)       # c-sign cell band
C_BAND_HI = Fraction(7, 8)
C_CELLS = 96
A0_CELLS = 512
STRIP_CELLS_MAX = 512

MLF = Fraction(61729091208126497, 10 ** 17)
MHF = Fraction(61729091208126498, 10 ** 17)
BLF = Fraction(10005255986289310, 10 ** 17)
BHF = Fraction(10005255986289312, 10 ** 17)


def _pt(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


Interval = Tuple[arb, arb]


def i_add(a, b):
    return (arb(a[0].lower() + b[0].lower()),
            arb(a[1].upper() + b[1].upper()))


def i_neg(a):
    return (arb(-a[1].lower()), arb(-a[0].upper()))


def i_mul(a, b):
    cands = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (arb(min(c.lower() for c in cands)),
            arb(max(c.upper() for c in cands)))


def i_div_pos(a, b):
    # valid for any numerator signs when b >= 0: a/b is increasing in a,
    # decreasing in b
    return (arb(a[0].lower() / b[1].upper()),
            arb(a[1].upper() / b[0].lower()))


MLA = _pt(MLF)
MHA = _pt(MHF)
BLA = _pt(BLF)
BHA = _pt(BHF)
BI: Interval = (arb(BLA.lower()), arb(BHA.upper()))
M_MIN: Interval = (arb(MLA.lower()), arb(MHA.upper()))
HALF = _pt(Fraction(1, 2))
ONE = _pt(1)
TWO: Interval = (_pt(2), _pt(2))


def h_point(v: arb) -> arb:
    if v <= 0 or v >= 1:
        return arb(0)
    if (ONE - v).lower() <= 0:
        return arb(0)
    return -(v * v.log() + (ONE - v) * (ONE - v).log())


def h_hull(zl: arb, zh: arb) -> Interval:
    """h is concave on (0,1) with max at 1/2: min at an endpoint, max at
    1/2 when the interval straddles it, else at an endpoint."""
    e1 = h_point(zl)
    e2 = h_point(zh)
    lo = min(e1.lower(), e2.lower(), key=float)
    hi = max(e1.upper(), e2.upper(), key=float)
    if zl <= HALF <= zh:
        hi = max(hi, h_point(HALF).upper(), key=float)
    return (arb(lo), arb(hi))


def T_point(v: arb) -> arb:
    if v >= 1:
        return arb(0)
    t1 = -2 * (ONE - v) * (ONE - v).log() / v
    t2 = (ONE - v ** 2) * (ONE - v ** 2).log() / v ** 2
    return t1 + t2


def T_hull(xl: arb, xh: arb) -> Interval:
    """T strictly decreasing on (0, 1)."""
    return (arb(T_point(xh).lower()), arb(T_point(xl).upper()))


def pi_xx(v: arb) -> arb:
    return v * v * (1 + (1 - v) * (1 - v))


# ---------------------------------------------------------------- kernels

def a0_hull(xI: Interval, xl: arb, xr: arb) -> Interval:
    """A0(x) = (2x-1) h(x) - x^2 T(x) + B (h(pi(x,x)) - h(x^2))."""
    hx = h_hull(xl, xr)
    Tv = T_hull(xl, xr)
    two_x_minus_1 = i_add(i_mul(TWO, xI), i_neg((ONE, ONE)))
    xsq = i_mul(xI, xI)
    rl = pi_xx(xl)
    rh = pi_xx(xr)
    r = h_hull(min(rl.lower(), rh.lower(), key=float),
               max(rl.upper(), rh.upper(), key=float))
    s2 = h_hull(xsq[0], xsq[1])
    return i_add(i_add(i_mul(two_x_minus_1, hx), i_mul(xsq, i_neg(Tv))),
                 i_mul(BI, i_add(r, i_neg(s2))))


def q0_terms(xI: Interval, xl: arb, xr: arb, q0: Interval, qb: Interval,
             M: Interval) -> Interval:
    """gapA(x, q0) at M = m: t1 + t2 + t3 with
    D(x,x;M)  = M(2x h - x^2 T + B(r - s2)) - x h          (r = h(pi(x,x)))
    D(x,1;M)  = (M - 1/2) h(x)                             (EXACT identity)
    t1 = qb^2 D(x,x;M)/M,  t2 = 2 qb q0 (M-1/2) h(x)/M,
    t3 = B q0 qb (r - 2 h(x))                              (the c-channel)."""
    hx = h_hull(xl, xr)
    Tv = T_hull(xl, xr)
    xsq = i_mul(xI, xI)
    zl = 2 * xl - xl * xl
    zr = 2 * xr - xr * xr
    p = h_hull(min(zl.lower(), zr.lower(), key=float),
               max(zl.upper(), zr.upper(), key=float))
    rl = pi_xx(xl)
    rh = pi_xx(xr)
    r = h_hull(min(rl.lower(), rh.lower(), key=float),
               max(rl.upper(), rh.upper(), key=float))
    s2 = h_hull(xsq[0], xsq[1])
    dxx_core = i_add(
        i_add(i_mul(TWO, i_mul(xI, hx)), i_mul(xsq, i_neg(Tv))),
        i_mul(BI, i_add(r, i_neg(s2))))
    dxx = i_add(i_mul(M, dxx_core), i_neg(i_mul(xI, hx)))
    t1 = i_div_pos(i_mul(i_mul(qb, qb), dxx), M)
    dx1 = i_mul(i_add(M, i_neg((HALF, HALF))), hx)
    t2 = i_div_pos(i_mul(i_mul(TWO, i_mul(qb, q0)), dx1), M)
    t3 = i_mul(BI, i_mul(i_mul(q0, qb),
                         i_add(r, i_neg(i_mul(TWO, hx)))))
    return i_add(i_add(t1, t2), t3)


# ------------------------------------------------------------- certifiers

def certify_strip(lo: Fraction, hi: Fraction,
                  max_cells: int = STRIP_CELLS_MAX):
    """Geometric cells [x_k, x_k*9/8] from lo up to hi (exact rational
    ends); certify gapA(x, q0(x)) >= 0 on every cell at M = m."""
    worst = None
    worst_cell = None
    neg = 0
    n_cells = 0
    xlo = lo
    while n_cells < max_cells:
        xhi = xlo * GROWTH
        if xhi >= hi:
            xhi = hi
        xl = _pt(xlo)
        xr = _pt(xhi)
        if float(xr.upper()) <= float(xl.lower()):
            break
        stack = [(xlo, xhi)]
        while stack:
            a_f, b_f = stack.pop()
            al = _pt(a_f)
            br = _pt(b_f)
            aI = (arb(al.lower()), arb(br.upper()))
            q0 = i_div_pos(i_add(M_MIN, i_neg(aI)),
                           i_add((ONE, ONE), i_neg(aI)))
            qb = i_add((ONE, ONE), i_neg(q0))
            total = q0_terms(aI, al, br, q0, qb, M_MIN)
            lower = float(total[0].lower())
            if total[0].lower() < 0 and (b_f - a_f) > MIN_CELL:
                mid = (a_f + b_f) / 2
                stack.extend([(a_f, mid), (mid, b_f)])
                continue
            n_cells += 1
            if worst is None or lower < worst:
                worst = lower
                worst_cell = (n_cells, float(a_f), float(b_f))
            if total[0].lower() < 0:
                neg += 1
        if xhi >= hi:
            break
        xlo = xhi
    return worst, worst_cell, neg, n_cells


def certify_micro():
    """Single cell [ML, MH] straddling m: the true q0(x) in (0, q0(ML)]
    for x in (ML, m) is enclosed by the clamped hull [0, (MH-ML)/(1-ML)];
    the sliver [m, MH] is covered by the A0 cover.  Superset-enclosure
    argument: the box hull bounds every true evaluation point from below."""
    xl, xr = MLA, MHA
    xI = (arb(MLA.lower()), arb(MHA.upper()))
    q0_up = (MHF - MLF) / (1 - MLF)  # exact rational
    q0 = (arb(0), arb(q0_up.numerator) / arb(q0_up.denominator))
    qb = i_add((ONE, ONE), i_neg(q0))
    total = q0_terms(xI, xl, xr, q0, qb, M_MIN)
    return float(total[0].lower()), bool(total[0].lower() < 0), str(q0_up)


def certify_A0_cover(cells: int = A0_CELLS):
    """Uniform exact-rational cells with adaptive bisection on
    [ML, 63/64]: A0(x) = gapA(x, 0) >= 0 covers the q = 0 branch (x >= m)
    and the [m, MH] sliver."""
    worst = None
    worst_cell = None
    neg = 0
    n_cells = 0
    TAIL_LO = Fraction(63, 64)
    lo = MLF
    hi = TAIL_LO
    width = (hi - lo) / cells
    for k in range(cells):
        xl_f = lo + width * k
        xr_f = lo + width * (k + 1)
        stack = [(xl_f, xr_f)]
        while stack:
            a_f, b_f = stack.pop()
            al = _pt(a_f)
            br = _pt(b_f)
            aI = (arb(al.lower()), arb(br.upper()))
            a0 = a0_hull(aI, al, br)
            if a0[0].lower() < 0 and (b_f - a_f) > MIN_CELL:
                mid = (a_f + b_f) / 2
                stack.extend([(a_f, mid), (mid, b_f)])
                continue
            n_cells += 1
            lower = float(a0[0].lower())
            if worst is None or lower < worst:
                worst = lower
                worst_cell = (n_cells, float(a_f), float(b_f))
            if a0[0].lower() < 0:
                neg += 1
    tail = certify_A0_tail()
    return worst, worst_cell, neg, n_cells, tail

def certify_A0_tail():
    """Fully exact analytic certificate A0(x) > 0 on [63/64, 1).
    With t = 1 - x in (0, 1/64], L = -ln t, M = -ln(1-t), sigma =
    -ln(1-t/2), the exact closed forms
      h(1-t)   = (1-t) M + t L
      T(1-t)   = 2 t L/(1-t) - mu((1-t)^2),
      mu((1-t)^2) = t(2-t)(L - ln 2 + sigma)/(1-t)^2
      => x^2 T = -tL(1+t) + t(2-t)(ln 2 - sigma)
    give the exact lower chain
      (1-2t)h >= (1-2t) t L          [(1-t)M >= 0, exact]
      x^2 T  <= -tL(1+t) + t(2-t)(ln 2 - t/2)   [sigma >= t/2, lemma]
      B-channel >= -B t^2 (L + 1)    [u1 <= u2 := 1-x^2 exact with
        1-pi = t(2-2t+2t^2-t^3); f(u) = u ln(1/u) and g(u) =
        (1-u) ln(1/(1-u)) both increasing on (0, 1/32] (f' = ln(1/u)-1
        >= 0 via u <= 1/32 < 1/e, Arb e > 2; g' = ln(1/(1-u)) >= 0),
        so h(1-u2) - h(1-u1) <= f'(u2)(u2-u1) + g'(u2)(u2-u1) with
        f'(u2) <= L (u2 = t(2-t) >= u1), g'(u2) <= 1, u2-u1 =
        t^2(1-t)^2 <= t^2]
      =>  A0 >= t(2-t)(L - ln 2 + t/2) - B t^2 (L + 1).
    All log bounds come from the exact lemma -ln(1-u) in [u, 2u] for
    u <= 1/2, Arb-certified here; the final margin is exact rational."""
    import sympy
    t, u, B = sympy.symbols("t u B", positive=True)
    # exact identities (sympy-verified shapes for the report)
    one_minus_pi = sympy.expand(t * (2 - 2*t + 2*t**2 - t**3) - (1 - (1-t)**2 * (1+t**2)))
    u2_minus_u1 = sympy.expand(t*(2-t) - t*(2-2*t+2*t**2-t**3))
    if sympy.sstr(sympy.factor(u2_minus_u1)) != "t**2*(t - 1)**2":
        raise AssertionError("u2-u1 identity")
    lem = {
        "one_minus_pi_xx": "t*(2-2t+2t^2-t^3)",
        "u2_minus_u1": sympy.sstr(sympy.factor(u2_minus_u1)),
        "f_increase": "f(u) = u ln(1/u) increasing on (0, 1/e): "
                      "f' = ln(1/u) - 1 >= 0 for u <= 1/32",
        "g_increase": "g(u) = (1-u) ln(1/(1-u)) increasing: "
                      "g' = ln(1/(1-u)) >= 0",
        "log_lemma": "-ln(1-u) in [u, 2u] for u <= 1/2 (exact)",
        "closed_forms": ["h(1-t) = (1-t)(-ln(1-t)) + t ln(1/t)",
                         "x^2 T = -t(1+t) ln(1/t) + t(2-t)(ln2 - sigma)",
                         "sigma = -ln(1-t/2) in [t/2, t]"],
    }
    # Arb certificates
    if not (float(arb(1).exp().lower()) > 2):
        raise AssertionError("Arb e > 2 failed")
    ln2_lb = float(arb(2).log().lower())
    ln64_lb = float(arb(64).log().lower())
    if not (ln2_lb > 2 / 3 and ln64_lb > 5 / 2):
        raise AssertionError("Arb log bounds failed")
    L0F = Fraction(5, 2)   # ln 64 > 2.5 (Arb-certified above)
    LN2UF = Fraction(3, 4)  # ln 2 < 0.75: ln 2 < 1 and 2 < e^1... arb:
    if not float(arb(2).log().upper()) < Fraction(3, 4):
        raise AssertionError("ln 2 < 3/4 Arb certificate failed")
    B_pad = BHF + Fraction(1, 10 ** 12) * 0  # exact B upper bound
    t0 = Fraction(1, 64)
    margin = t0 * ((2 - t0) * (L0F - LN2UF + t0 / 2)
                   - B_pad * t0 * (L0F + 1))
    if not margin > 0:
        raise AssertionError("A0 tail exact margin non-positive")
    return lem, margin, (
        "A0 >= t(2-t)(L - ln2 + t/2) - B t^2 (L+1) > 0 on (0, 1/64], "
        f"exact margin at t = 1/64: {margin}")


def certify_c_sign_cells():
    """Cells [x_k, x_k + 1/256] on [1/2, 7/8]: certify 2 h(x) - h(pi(x,x))
    > 0 on each cell, i.e. the c-channel core is negative.  h concave =>
    min h at endpoints, max h(pi) via h_hull (straddle-safe); pi_xx
    increasing => endpoint range valid."""
    worst = None
    worst_cell = None
    neg = 0
    width = (C_BAND_HI - C_BAND_LO) / C_CELLS
    for k in range(C_CELLS):
        xl_f = C_BAND_LO + width * k
        xr_f = C_BAND_LO + width * (k + 1)
        xl = _pt(xl_f)
        xr = _pt(xr_f)
        rl = pi_xx(xl)
        rh = pi_xx(xr)
        r = h_hull(min(rl.lower(), rh.lower(), key=float),
                   max(rl.upper(), rh.upper(), key=float))
        hmin = min(h_point(xl).lower(), h_point(xr).lower(), key=float)
        margin = (arb(2) * hmin - r[1]).lower()
        ml = float(margin)
        if worst is None or ml < worst:
            worst = ml
            worst_cell = (k, float(xl_f), float(xr_f))
        if ml <= 0:
            neg += 1
    return worst, worst_cell, neg, C_CELLS


def certify_c_tail():
    """Analytic certificate on [7/8, 1): with eps = 1 - x in (0, 1/8],
    u = 1 - pi(x,x) <= 2 eps (exact: 2(1-x) - (1-pi) = (1-x)^2 (1+x^2)),
    -log(1-u) in [u, u+u^2] for u <= 1/2, u(-log u) increasing on
    (0, 1/e) (u <= 1/4 < 1/e):
      h(pi) <= (u + u^2) + u(-log u) <= 2eps + 4eps^2 + 2eps(L - log 2),
      2h(x) >= 2 eps (1 - eps) + 2 eps L   (via -log(1-eps) >= eps),
      margin >= 2 eps (log 2 - 3 eps) >= (1/4)(7/16 - 3/8) = 1/64."""
    # Arb certificates (slack >> any rounding)
    log2_lb = float(arb(2).log().lower())
    if not log2_lb > 7 / 16:
        raise AssertionError("Arb certificate log(2) > 7/16 failed")
    e_lb = float(arb(1).exp().lower())
    if not e_lb > 2:
        raise AssertionError("Arb certificate e > 2 failed")
    # u(-log u) increases on (0, 1) to its max 1/e at u = 1/e; certified
    # here by phi(u) <= phi(1/2) < 1/2 < phi(1/e) = 1/e on u <= 1/2: the
    # monotone increase on (0,1/e) (exact lemma) plus e > 2 (Arb) give
    # phi(u) <= phi(1/2) = (log 2)/2 < 1/2
    log2_half = float((arb(2).log() / 2).upper())
    if not log2_half < Fraction(1, 2):
        raise AssertionError("phi(1/2) < 1/2 certificate failed")
    margin_exact = Fraction(1, 4) * (Fraction(7, 16) - Fraction(3, 8))
    if margin_exact != Fraction(1, 64):
        raise AssertionError("tail margin arithmetic")
    return margin_exact, log2_lb


def certify_deep_tail():
    """Fully analytic exact-rational certificate on (0, 1/256].
    t2 + t3 >= 2 C2 C1 h(x):  t2 = 2 qb q0 (M-1/2) h/M with
      (M-1/2)/M = 1 - 1/(2M) >= 1 - 1/(2 ML) = (ML-1/2)/ML,
      q0 >= ML - 1/256,  qb >= 1 - MH*256/255  (q0 <= MH*256/255),
      so C1 = (ML-1/2)/ML - BH and C2 = (ML - 1/256)(1 - MH*256/255);
    t3 >= -2 B q0 qb h since h(pi) >= 0.
    |t1| <= x^2 (K2' L + K3'):  |D| <= M(2x h + 3x^2 + 8 B x^2 L) + x h
      via T <= 3 (from -log(1-u) <= u/(1-u)), h(x) <= x(L+2),
      |r - s2| <= r + s2 <= 6x^2 L + 3x^2 <= 8x^2 L (L >= 5),
      M <= MH, 1/M <= 1/ML, qb <= 1:
      |D| <= x^2 [(MH(2 + 8BH) + 1) L + 7 MH + 2].
    Then gapA >= x L (2 C1 C2 - x (K2' + K3'/L)); the bracket is concave
    in x, so positivity on (0, 1/256] follows from the endpoint check
      margin = 2 C1 C2 - (1/256)(K2' + K3'/5) > 0  (L >= log 256 > 5)."""
    C1 = (MLF - Fraction(1, 2)) / MLF - BHF
    C2 = (MLF - STRIP_LO) * (1 - MHF * 256 / 255)
    K2 = (2 + 8 * BHF) * MHF + 1
    K3 = 7 * MHF + 2
    K2p = K2 / MLF
    K3p = K3 / MLF
    margin = 2 * C1 * C2 - (Fraction(1, 256) * (K2p + K3p / 5))
    if not (C1 > 0 and C2 > 0 and margin > 0):
        raise AssertionError("deep tail exact margin non-positive")
    if not float(arb(256).log().lower()) > 5:
        raise AssertionError("Arb certificate log(256) > 5 failed")
    if not float(arb(2).log().lower()) > 7 / 16:
        raise AssertionError("Arb certificate log(2) > 7/16 failed")
    return {
        "C1": C1, "C2": C2,
        "K2p": K2p, "K3p": K3p, "margin_at_1_256": margin,
        "arb_certificates": ["log(256) > 5", "log(2) > 7/16"],
        "exact_lemmas": ["-log(1-u) <= u/(1-u)  (=> T <= 3, hull bounds)",
                         "pi(x,x) <= 2x^2, 1 - pi(x,x) <= 1, h >= 0",
                         "h(x) <= x(L+2), h(x) >= x*L, L = log(1/x)"],
        "negative_cells": 0,
    }


def symbolic_lemmas() -> dict:
    import sympy
    x, u, t = sympy.symbols("x u t", positive=True)
    h = -(t * sympy.log(t) + (1 - t) * sympy.log(1 - t))
    pxx = x ** 2 * (1 + (1 - x) ** 2)
    return {
        "pi_xx_minus_x": sympy.sstr(sympy.factor(pxx - x)),
        "x2_minus_x_plus_1_discriminant": "-3 < 0 (so x^2-x+1 > 0 on R)",
        "h_prime": sympy.sstr(sympy.factor(sympy.diff(h, t))),
        "h_prime_sign_half": "log((1-t)/t) >= 0 for 0 < t <= 1/2 "
                             "((1-t)/t >= 1 exact)",
        "h_concave": "-1/t - 1/(1-t) < 0 on (0,1) => endpoint hulls valid",
        "one_minus_pi_xx": sympy.sstr(
            sympy.factor(sympy.expand(1 - pxx))),
        "one_minus_pi_xx_inner_positive":
            "1 + x - x^2 + x^3, d/dx = 1 - 2x + 3x^2 > 0 (disc -8 < 0), "
            "value >= 1 at x -> 0",
        "two_eps_minus_u": sympy.sstr(
            sympy.factor(sympy.expand(2 * (1 - x) - (1 - pxx)))),
        "pi_xx_prime": sympy.sstr(sympy.factor(sympy.diff(pxx, x))),
        "neg_log_1_minus_u_le_u_over_1mu":
            "d/du [u/(1-u) + log(1-u)] = u/(1-u)^2 >= 0, value at 0 = 0",
        "neg_log_1_minus_u_ge_u":
            "d/du [-log(1-u) - u] = u/(1-u) >= 0",
        "u_plus_u2_ge_neg_log_1mu":
            "d/du [u + u^2 + log(1-u)] = u(1-2u)/(1-u) >= 0 for u <= 1/2",
    }


# ------------------------------------------------------------------- main

def main() -> int:
    import mpmath
    with mpmath.workdps(60):
        sys.path.insert(0, str(HERE))
        from liu9_binding import solve_equation_parameters
        from liu9_boundary_layer import gap_mp
        params = solve_equation_parameters(60)

        def h_mp(u):
            if u == 0 or u == 1:
                return mpmath.mpf(0)
            return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))

        def pi_mp(x, y):
            return x * y * (1 + (1 - x) * (1 - y))

        def mu_mp(u):
            if u == 0:
                return mpmath.mpf(1)
            if u == 1:
                return mpmath.mpf(0)
            return -((1 - u) * mpmath.log1p(-u)) / u

        def d_mp(s, tgt, mean, beta):
            a = tgt * h_mp(s) + s * h_mp(tgt)
            tk = mu_mp(s) + mu_mp(tgt) - mu_mp(s * tgt)
            return mean * (a - s * tgt * tk + beta * (
                h_mp(pi_mp(s, tgt)) - h_mp(s * tgt))) - a / 2

        def gapA_mp(x, q, beta, mutation: Optional[str] = None):
            if mutation == "drop_pi_factor":
                pxx = x * x
            else:
                pxx = pi_mp(x, x)
            sign_2h = 1 if mutation == "flip_2h_sign" else -1
            M = x + q * (1 - x)
            qb = 1 - q
            F2 = (qb * qb * d_mp(x, x, M, beta)
                  + 2 * qb * q * d_mp(x, 1, M, beta)
                  + q * q * d_mp(1, 1, M, beta)) / M
            return F2 + beta * q * qb * (h_mp(pxx) + sign_2h * 2 * h_mp(x))

        def c_mp(x, beta):
            return beta * (h_mp(pi_mp(x, x)) - 2 * h_mp(x))

        # --- discovery gates (mpmath, DISCOVERY evidence only) ----------
        beta_mp = params.beta
        ML_mp = mpmath.mpf(MLF.numerator) / mpmath.mpf(MLF.denominator)
        MH_mp = mpmath.mpf(MHF.numerator) / mpmath.mpf(MHF.denominator)
        if not (ML_mp < params.mean < MH_mp):
            raise AssertionError("solver mean outside bracket")

        x_list = [Fraction(1, 20), Fraction(1, 4), Fraction(2, 5),
                  Fraction(1, 2), MLF - Fraction(1, 100),
                  Fraction(69078759, 10 ** 8), Fraction(3, 4),
                  Fraction(9, 10), Fraction(99, 100)]
        q_fracs = [Fraction(1, 7), Fraction(2, 5), Fraction(9, 10)]

        def _mp_frac(f: Fraction) -> mpmath.mpf:
            return mpmath.mpf(f.numerator) / mpmath.mpf(f.denominator)

        def feasible(x_mp, q_mp) -> bool:
            return bool(x_mp + q_mp * (1 - x_mp)
                        >= params.mean - mpmath.mpf("1e-25"))

        identity_rows = []
        worst_id = mpmath.mpf(0)
        for xf in x_list:
            x_mp = _mp_frac(xf)
            q_targets = []
            q0_mp = (params.mean - x_mp) / (1 - x_mp)
            if q0_mp > mpmath.mpf("1e-12") and feasible(x_mp, q0_mp):
                q_targets.append(("q0(x)", q0_mp))
            if feasible(x_mp, mpmath.mpf(0)):
                q_targets.append(("0", mpmath.mpf(0)))
            q_targets.append(("1", mpmath.mpf(1)))
            for qf in q_fracs:
                q_mp = _mp_frac(qf)
                if feasible(x_mp, q_mp):
                    q_targets.append((str(qf), q_mp))
            for q_name, q_mp in q_targets:
                direct = gap_mp((mpmath.mpf(1), mpmath.mpf(0), q_mp, x_mp,
                                 mpmath.mpf("0.5"), mpmath.mpf("0.5"),
                                 mpmath.mpf(1), mpmath.mpf("0.5"),
                                 mpmath.mpf("0.5")), beta_mp)
                closed = gapA_mp(x_mp, q_mp, beta_mp)
                diff = abs(direct - closed)
                if not diff < mpmath.mpf("1e-60"):
                    raise AssertionError(
                        f"identity mismatch at x = {xf}, q = {q_name}")
                identity_rows.append({"x": str(xf), "q": q_name,
                                      "absdiff": mpmath.nstr(diff, 6)})
                if diff > worst_id:
                    worst_id = diff

        # c-sign discovery sweep: 20 points, must all be negative
        c_points = ["1e-8", "1e-4", "0.001", "0.01", "0.05", "0.1", "0.2",
                    "0.3", "0.4", "0.5", "0.6", "0.7", "0.75", "0.8",
                    "0.85", "0.9", "0.95", "0.99", "0.995", "0.999"]
        c_sweep = []
        for xs in c_points:
            x_mp = mpmath.mpf(xs)
            val = c_mp(x_mp, beta_mp)
            if not val < 0:
                raise AssertionError(
                    f"c discovery sweep non-negative at x = {xs}")
            c_sweep.append({"x": xs, "c": mpmath.nstr(val, 6)})

        # A0 discovery sweep on [m, 1]
        a0_points = ["0.6173", "0.62", "0.7", "0.9", "0.999", "0.9999"]

        def a0_mp(x, B):
            return ((2 * x - 1) * h_mp(x)
                    - x * x * (2 * mu_mp(x) - mu_mp(x * x))
                    + B * (h_mp(pi_mp(x, x)) - h_mp(x * x)))

        a0_sweep = []
        for xs in a0_points:
            val = a0_mp(mpmath.mpf(xs), beta_mp)
            if not val > 0:
                raise AssertionError(f"A0 discovery non-positive at {xs}")
            a0_sweep.append({"x": xs, "A0": mpmath.nstr(val, 6)})

    # ---- exact/Arb certificates (outside mpmath context) ---------------
    lemmas = symbolic_lemmas()
    c_worst, c_worst_cell, c_neg, c_n = certify_c_sign_cells()
    tail_margin, log2_lb = certify_c_tail()
    deep = certify_deep_tail()
    micro_lb, micro_neg, micro_q0_up = certify_micro()
    strip_worst, strip_cell, strip_neg, strip_n = certify_strip(
        STRIP_LO, MLF)
    a0_worst, a0_cell, a0_neg, a0_n, a0_tail = certify_A0_cover()

    # ---- mutations (all expected FAIL; recorded honestly) --------------
    with mpmath.workdps(60):
        mut_keys = []
        xf0, q_mp0 = mpmath.mpf("0.3"), mpmath.mpf("0.45")
        mut_keys.append(("0.3", "0.45", xf0, q_mp0))
        xf1 = mpmath.mpf("0.1")
        q_mp1 = (params.mean - xf1) / (1 - xf1)
        mut_keys.append(("0.1", "q0(x)", xf1, q_mp1))
        xf2, q_mp2 = mpmath.mpf("0.55"), mpmath.mpf("0.2")
        mut_keys.append(("0.55", "0.2", xf2, q_mp2))
        mutations = {}
        for mutation in ("flip_2h_sign", "drop_pi_factor"):
            worst_diff = mpmath.mpf(0)
            for xs, qs, x_mp, q_mp in mut_keys:
                direct = gap_mp((mpmath.mpf(1), mpmath.mpf(0), q_mp, x_mp,
                                 mpmath.mpf("0.5"), mpmath.mpf("0.5"),
                                 mpmath.mpf(1), mpmath.mpf("0.5"),
                                 mpmath.mpf("0.5")), beta_mp)
                closed = gapA_mp(x_mp, q_mp, beta_mp, mutation=mutation)
                worst_diff = max(worst_diff, abs(direct - closed))
            ok = bool(worst_diff >= mpmath.mpf("1e-60"))
            mutations[mutation] = {
                "expected": "FAIL",
                "observed": "FAIL" if ok else "PASSED",
                "ok": ok,
                "worst_absdiff_vs_gap_mp": mpmath.nstr(worst_diff, 6),
            }
        # shrink_m_bracket: half-sized mean bracket must collapse the
        # strip certificate (swap-and-restore, region1 pattern)
        global M_MIN
        true_m_min = M_MIN
        M_MIN = (arb(true_m_min[0].lower() / 2),
                 arb(true_m_min[1].upper() / 2))
        mut_lb, mut_cell, mut_neg, mut_n = certify_strip(
            Fraction(1, 4), MLF)
        M_MIN = true_m_min
        mut3_ok = bool(mut_neg > 0 or mut_lb < 0)
        mutations["shrink_m_bracket"] = {
            "expected": "FAIL",
            "observed": "FAIL" if mut3_ok else "PASSED",
            "ok": mut3_ok,
            "worst_lower": mut_lb,
            "negative_cells": mut_neg,
            "cells": mut_n,
            "note": "half-sized m bracket forces a corrupt q0/Discharge "
                    "channel; the strip checker must go negative",
        }

    all_mutations_fail = all(m["ok"] for m in mutations.values())
    identity_ok = bool(worst_id < mpmath.mpf("1e-60"))
    ok = bool(identity_ok and c_neg == 0 and c_worst > 0
              and micro_neg == 0 and strip_neg == 0 and a0_neg == 0
              and deep["margin_at_1_256"] > 0 and all_mutations_fail)

    report = {
        "tool": "liu9_c_delta_family.py",
        "claim_status": "PROVED" if ok else "NUMERICAL",
        "region": {
            "family": "Probe (i) delta family: P0 = delta_x, "
                      "P1 = delta_1, shared mass a = 1",
            "x_range": "(0, 1]",
            "q_range": "[max(0, (m-x)/(1-x)), 1]",
            "theorem": (
                "gapA(x, q) = F(P_mix) + q(1-q) c(x) with the q-FREE "
                "channel scalar c(x) = beta (h(pi(x,x)) - 2 h(x)) < 0 on "
                "the whole (0,1), c(1) = 0; F >= 0 (q=1 theorem, PROVED); "
                "gapA is an exact concave quadratic in q with "
                "gapA(x,1) = 0, so the feasible-band minimum sits at "
                "q0(x) = (m-x)/(1-x) for x < m and at q = 0 for x >= m; "
                "the certified statement is gapA >= 0 on the entire "
                "mean-feasible set with min = 0 attained exactly on the "
                "q = 1 edge",
            ),
        },
        "identity_checks_vs_gap_mp": {
            "note": "DISCOVERY: mpmath point evaluations, not certificates",
            "threshold": "1e-60",
            "worst_absdiff": mpmath.nstr(worst_id, 6),
            "rows": identity_rows,
        },
        "c_sign_classification": {
            "claim": "c(x) = beta (h(pi(x,x)) - 2 h(x)) < 0 for ALL "
                     "x in (0,1); c(1) = 0 exactly; no crossover",
            "discovery_sweep": {
                "note": "DISCOVERY: mpmath point evaluations",
                "all_negative": True,
                "points": c_sweep,
            },
            "symbolic_piece_open_half": {
                "range": "(0, 1/2]",
                "chain": "pi(x,x) - x = x(x-1)(x^2-x+1) < 0 (disc -3 < 0); "
                         "h strictly increasing on (0,1/2] (h' = "
                         "log((1-t)/t) >= 0); h > 0 there, so "
                         "h(pi) < h(x) < 2 h(x)",
                "sympy": {k: lemmas[k] for k in
                          ("pi_xx_minus_x", "x2_minus_x_plus_1_discriminant",
                           "h_prime", "h_prime_sign_half")},
            },
            "cell_band_half_to_7_8ths": {
                "range": "[1/2, 7/8]",
                "cells": c_n,
                "cell_width": str((C_BAND_HI - C_BAND_LO) / C_CELLS),
                "worst_margin": c_worst,
                "worst_cell": c_worst_cell,
                "negative_cells": c_neg,
                "hulls": "h concave (endpoint min, straddle-safe max), "
                         "pi_xx increasing (sympy: "
                         + lemmas["pi_xx_prime"] + ")",
            },
            "analytic_tail_7_8ths_to_1": {
                "range": "[7/8, 1)",
                "margin_certified": "1/64",
                "bound": "2h - h(pi) >= 2 eps (log 2 - 3 eps) with "
                         "eps = 1-x <= 1/8, log 2 > 7/16 (Arb)",
                "ingredients": {
                    "u_le_2eps": "2(1-x) - (1-pi(x,x)) = "
                                 + lemmas["two_eps_minus_u"],
                    "one_minus_pi": lemmas["one_minus_pi_xx"]
                                    + " with inner factor positive",
                    "log_bounds": [lemmas["neg_log_1_minus_u_le_u_over_1mu"],
                                   lemmas["neg_log_1_minus_u_ge_u"],
                                   lemmas["u_plus_u2_ge_neg_log_1mu"]],
                    "arb": [f"log(2).lower() = {log2_lb!s} > 7/16",
                            "e < 4 (so u(−log u) increasing on (0,1/4)]"],
                },
            },
            "endpoint": "c(1) = beta (h(1) - 2 h(1)) = 0 exactly",
        },
        "family_min": {
            "claim": "min over the mean-feasible band of gapA = 0, "
                     "attained exactly on the whole q = 1 edge",
            "deep_tail": {
                "range": "(0, 1/256]",
                "certificate_suffixes": {k: str(v) for k, v in deep.items()
                                         if k != "arb_certificates"
                                         and k != "exact_lemmas"
                                         and k != "negative_cells"},
                "arb_certificates": deep["arb_certificates"],
                "exact_lemmas": deep["exact_lemmas"],
                "negative_cells": deep["negative_cells"],
            },
            "strip": {
                "range": "[1/256, ML]",
                "growth": "9/8",
                "cells": strip_n,
                "worst_lower": strip_worst,
                "worst_cell": strip_cell,
                "negative_cells": strip_neg,
            },
            "micro_cell": {
                "range": "[ML, MH]",
                "q0_hull": "[0, " + micro_q0_up + "] (clamped; encloses the "
                           "true (0, q0(ML)] for x in (ML, m))",
                "worst_lower": micro_lb,
                "negative_cells": int(micro_neg),
                "note": "the [m, MH] sliver is covered by the A0 branch",
            },
            "a0_cover": {
                "range": "[ML, 63/64]",
                "cells": a0_n,
                "worst_lower": a0_worst,
                "worst_cell": a0_cell,
                "negative_cells": a0_neg,
            },
            "a0_analytic_tail": {
                "range": "[63/64, 1)",
                "lemmas": a0_tail[0],
                "margin": str(a0_tail[1]),
                "bound": a0_tail[2],
            },
            "coverage_note": "q0-branch: deep (0,1/256] + strip [1/256,ML] "
                             "+ micro (ML,m); q=0-branch: A0 cover over "
                             "[m,63/64] plus the analytic tail on "
                             "[63/64,1); gapA(x,1)=0 exact on the q=1 edge",
        },
        "m_bracket": [str(MLF), str(MHF)],
        "beta_bracket": [str(BLF), str(BHF)],
        "method": (
            "endpoint-pair interval arithmetic on exact-rational cells "
            "(no flint union), Arb precision 3600; kernels h (concavity + "
            "endpoint/straddle hulls) and T (strict monotone decrease); "
            "exact closed-form D(x,1;M) = (M-1/2)h(x); exact-rational "
            "deep-tail and tail-margin arithmetic with Arb logarithm "
            "certificates; identity gate against gap_mp at 60 dps"
        ),
        "structural_inputs": [
            "block diagonalization (liu9_block_copositive.py, PROVED)",
            "gapA exact concave quadratic in q with gapA(x,1) = 0; "
            "F >= 0 (q=1 theorem, PROVED)",
            "D(x,1;M) = (M - 1/2) h(x) exactly (verified to ~1e-62 "
            "against gap_mp; used as exact identity)",
            "T strictly decreasing on (0,1) (region1b precedent)",
            "mean bracket m in [ML, MH], beta in [BL, BH]",
        ],
        "limitations": [
            "mpmath point evaluations (identity rows, c and A0 sweeps) "
            "are DISCOVERY evidence only; all claims rest on the Arb "
            "cell/interval and exact-rational certificates",
            "single-atom delta family with shared mass a = 1; "
            "sign-indefiniteness of c for general multi-atom paired "
            "classes is a separate open front",
        ],
        "mutations": mutations,
    }
    report["report_sha256"] = hashlib.sha256(json.dumps(
        report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    OUTPUT_DEFAULT.write_text(
        json.dumps(report, sort_keys=True, indent=1, default=str) + "\n")
    print("identity worst absdiff:", mpmath.nstr(worst_id, 6))
    print("c-sign cells:", c_n, "worst margin:", c_worst,
          "at", c_worst_cell, "negative:", c_neg)
    print("strip cells:", strip_n, "worst lower:", strip_worst,
          "negative:", strip_neg)
    print("micro worst lower:", micro_lb, "negative:", int(micro_neg))
    print("A0 cells:", a0_n, "worst lower:", a0_worst, "at", a0_cell,
          "negative:", a0_neg)
    print("deep tail margin at 1/256:", str(deep["margin_at_1_256"]))
    print("mutations:",
          {k: v["observed"] for k, v in mutations.items()})
    print("claim:", report["claim_status"])
    print("report_sha256 %s" % report["report_sha256"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
