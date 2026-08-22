"""Certified lower-bound probes for the five-parameter Phi functional.

This file is self-contained on purpose: importing cert.py executes its long validation
suite.  All interval arithmetic used for lower bounds remains in Arb at 160-bit
precision; mpmath is used only for numerical reference checks.
"""

import os
import sys
import time

import numpy as np
import sympy as sp
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy import PSI, h as h_mp, mp, mpf


ctx.prec = 160

ALPHA = arb("0.0356069")
ONE = arb(1)
TWO = arb(2)
HALF = arb(1) / 2
LOG2 = TWO.log()
ZERO = arb(0)
B_OBS = mpf("0.32945473850303697239")
LAMBDAS = ("0", "0.5", "1.0", "1.626", "2.2", "3.0")
RH_GMAX = None


def _arb(x):
    if isinstance(x, arb):
        return arb(x)
    return arb(str(x))


# Copied from cert.py rather than importing that executable module.
def hull(*balls):
    out = balls[0]
    for ball in balls[1:]:
        out = out.union(ball)
    return out


def h_pt(x):
    """h at a point, as an Arb ball; h(0) = h(1) = 0 exactly."""
    x = _arb(x)
    if x <= 0 or x >= 1:
        return arb(0)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def h_encl(u, v):
    """Enclosure of the range of h over [u,v]."""
    u, v = _arb(u), _arb(v)
    enclosure = hull(h_pt(u), h_pt(v))
    if v <= HALF or u >= HALF:
        return enclosure
    return hull(enclosure, ONE)


def amin(a, b):
    a, b = _arb(a), _arb(b)
    if a <= b:
        return a
    if b <= a:
        return b
    return hull(a, b)


def amax(a, b):
    a, b = _arb(a), _arb(b)
    if a >= b:
        return a
    if b >= a:
        return b
    return hull(a, b)


def sstar_pt(p, r):
    p, r = _arb(p), _arb(r)
    return amin(amax(HALF, amax(p, r)), amin(p + r, ONE))


def sstar_encl(pl, ph, rl, rh):
    """sstar is coordinatewise nondecreasing, so two corners give its range."""
    return hull(sstar_pt(pl, rl), sstar_pt(ph, rh))


def rh_encl(u, v, gmax=None):
    """Enclose rh(z)=2(1-z)h(z)-h((1-z)^2) over [u,v]."""
    u, v = _arb(u), _arb(v)
    z = hull(u, v)
    hz = h_encl(u, v)
    x_lo, x_hi = (ONE - v) * (ONE - v), (ONE - u) * (ONE - u)
    h2 = h_encl(x_lo, x_hi)
    enclosure = TWO * (ONE - z) * hz - h2
    lo = amax(enclosure.lower(), ZERO)
    hi = enclosure.upper() if gmax is None else amin(enclosure.upper(), gmax)
    return hull(lo, hi)


def sqrt_rh_enc(u, v, gmax=None):
    """NaN-free two-sided enclosure of sqrt(rh) over [u,v]."""
    hi = rh_encl(u, v, gmax).upper()
    # Never sqrt a ball that may straddle zero: hi is a positive endpoint.
    return hull(ZERO, hi.sqrt())


def _rh_global(n=20000):
    """Rigorous upper bound for max rh from an Arb cover of [0,1]."""
    best = arb(0)
    for i in range(n):
        enclosure = rh_encl(arb(i) / n, arb(i + 1) / n)
        upper = enclosure.upper()
        if not (upper <= best):
            best = amax(best, upper)
    return best.upper()


def _prepare_box(box):
    (p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh) = box
    p1l, p1h = _arb(p1l), _arb(p1h)
    q1l, q1h = _arb(q1l), _arb(q1h)
    p2l, p2h = _arb(p2l), _arb(p2h)
    q2l, q2h = _arb(q2l), _arb(q2h)
    wl, wh = _arb(wl), _arb(wh)

    # The orbit parametrization uses p_i <= q_i.  As in cert.py, this rectangle
    # is a safe superset of the ordered part of the input box.
    q1l, q2l = amax(q1l, p1l), amax(q2l, p2l)
    if q1l > q1h or q2l > q2h:
        return None
    return p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh


def phi_lam_encl(box, t, lam):
    """Enclose Phi + lam*(M-t) on the full box, with lam >= 0.

    On feasible points M <= t, the shift is nonpositive, so the lower endpoint
    is also a lower bound for Phi on the feasible part of the box.
    """
    prepared = _prepare_box(box)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    t, lam = _arb(t), _arb(lam)
    if not (lam >= ZERO):
        raise ValueError("lambda must be certified nonnegative")

    w = hull(wl, wh)
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    b1, b2 = ONE - m1, ONE - m2
    B = w * b1 + (ONE - w) * b2
    # Deliberately no B >= 1-t intersection here.  This enclosure is over the
    # full box, which includes M > t; imposing the feasibility inequality would
    # cease to enclose Phi_lam there and invalidate the full-box infimum.
    L = (w * (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
         + (ONE - w) * (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO)
    L = hull(amax(L.lower(), ZERO), L.upper())
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    C = (w * h_encl(s1.lower(), s1.upper())
         + (ONE - w) * h_encl(s2.lower(), s2.upper()))
    C = hull(amax(C.lower(), ZERO), C.upper())
    r1 = (sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    r2 = (sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    sig = w * r1 + (ONE - w) * r2
    phi = ((TWO * (ONE - ALPHA) * B - ONE) * L + ALPHA * C
           - (ONE - ALPHA) * sig * sig)
    return phi + lam * (mean - t)


def phi_feasible_encl(box, t):
    """Original feasibility-using Phi enclosure, or None if infeasible."""
    prepared = _prepare_box(box)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    t = _arb(t)
    w = hull(wl, wh)
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    if mean.lower() > t:
        return None
    b1, b2 = ONE - m1, ONE - m2
    B = w * b1 + (ONE - w) * b2
    B = hull(amax(B.lower(), ONE - t), B.upper())
    L = (w * (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
         + (ONE - w) * (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO)
    L = hull(amax(L.lower(), ZERO), L.upper())
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    C = (w * h_encl(s1.lower(), s1.upper())
         + (ONE - w) * h_encl(s2.lower(), s2.upper()))
    C = hull(amax(C.lower(), ZERO), C.upper())
    r1 = (sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    r2 = (sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    sig = w * r1 + (ONE - w) * r2
    return ((TWO * (ONE - ALPHA) * B - ONE) * L + ALPHA * C
            - (ONE - ALPHA) * sig * sig)


def phi_true(p1, q1, p2, q2, w, alpha="0.0356069"):
    """Fifty-digit mpmath reference implementation of Phi."""
    old_dps = mp.dps
    mp.dps = 50
    try:
        al = mpf(alpha)
        P = [mpf(p1), mpf(q1), mpf(p2), mpf(q2)]
        wv = mpf(w)

        def rr(z):
            x = 1 - z
            return max(2 * x * h_mp(z) - h_mp(x * x), mpf(0))

        def st(a, b):
            return min(max(mpf("0.5"), max(a, b)), min(a + b, mpf(1)))

        B = (wv * (1 - (P[0] + P[1]) / 2)
             + (1 - wv) * (1 - (P[2] + P[3]) / 2))
        L = (wv * (h_mp(P[0]) + h_mp(P[1])) / 2
             + (1 - wv) * (h_mp(P[2]) + h_mp(P[3])) / 2)
        C = wv * h_mp(st(P[0], P[1])) + (1 - wv) * h_mp(st(P[2], P[3]))
        S = (wv * (mp.sqrt(rr(P[0])) + mp.sqrt(rr(P[1]))) / 2
             + (1 - wv) * (mp.sqrt(rr(P[2])) + mp.sqrt(rr(P[3]))) / 2)
        return ((2 * (1 - al) * B - 1) * L + al * C
                - (1 - al) * S * S)
    finally:
        mp.dps = old_dps


def hprime_encl(u, v):
    """Enclose h'(z)=log2((1-z)/z), or return None at an endpoint."""
    u, v = _arb(u), _arb(v)
    if not (u > ZERO and v < ONE):
        return None

    def hp(z):
        return ((ONE - z).log() - z.log()) / LOG2

    # h' is strictly decreasing.
    return hull(hp(v), hp(u))


def rhprime_encl(u, v):
    """Enclose rh' on a strict interior interval."""
    u, v = _arb(u), _arb(v)
    hp = hprime_encl(u, v)
    if hp is None:
        return None
    z = hull(u, v)
    x_lo, x_hi = (ONE - v) * (ONE - v), (ONE - u) * (ONE - u)
    hp2 = hprime_encl(x_lo, x_hi)
    if hp2 is None:
        return None
    # rh' = -2h(z) + 2(1-z)(h'(z) + h'((1-z)^2)).
    return -TWO * h_encl(u, v) + TWO * (ONE - z) * (hp + hp2)


def _positive_rh_encl(u, v):
    """Find a certified-positive rh enclosure by subdividing, if practical."""
    u, v = _arb(u), _arb(v)
    if not (u > ZERO and v < ONE):
        return None
    direct = rh_encl(u, v, RH_GMAX)
    if direct.lower() > ZERO:
        return direct
    for n in (8, 32, 128):
        step = (v - u) / n
        pieces = []
        positive = True
        for i in range(n):
            left = u + i * step
            right = u + (i + 1) * step
            piece = rh_encl(left, right, RH_GMAX)
            if not (piece.lower() > ZERO):
                positive = False
                break
            pieces.append(piece)
        if positive:
            return hull(*pieces)
    return None


def sqrt_rh_prime_encl(u, v):
    """Enclose d/dz sqrt(rh(z)), returning None if division may hit zero."""
    r = _positive_rh_encl(u, v)
    if r is None:
        return None
    rp = rhprime_encl(u, v)
    if rp is None:
        return None
    # r has certified-positive lower endpoint, so this sqrt cannot be NaN and
    # the denominator excludes zero.
    root = r.sqrt()
    if not root.is_finite() or not (root.lower() > ZERO):
        return None
    out = rp / (TWO * root)
    return out if out.is_finite() else None


def _sstar_branch(pl, ph, ql, qh):
    """Return one certified affine branch (lo,hi,dp,dq,name), else None."""
    if ph <= HALF and qh <= HALF:
        A = (HALF, HALF, ZERO, ZERO, "half")
    elif pl >= HALF and pl >= qh:
        A = (pl, ph, ONE, ZERO, "p")
    elif ql >= HALF and ql >= ph:
        A = (ql, qh, ZERO, ONE, "q")
    else:
        return None

    if ph + qh <= ONE:
        D = (pl + ql, ph + qh, ONE, ONE, "sum")
    elif pl + ql >= ONE:
        D = (ONE, ONE, ZERO, ZERO, "one")
    else:
        return None

    if A[1] <= D[0]:
        return A
    if D[1] <= A[0]:
        return D
    return None


def interval_gradient(box):
    """Interval enclosure of all five partial derivatives, or None."""
    prepared = _prepare_box(box)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    endpoints = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h))

    hps = [hprime_encl(lo, hi) for lo, hi in endpoints]
    gps = [sqrt_rh_prime_encl(lo, hi) for lo, hi in endpoints]
    if any(value is None for value in hps + gps):
        return None

    branch1 = _sstar_branch(p1l, p1h, q1l, q1h)
    branch2 = _sstar_branch(p2l, p2h, q2l, q2h)
    if branch1 is None or branch2 is None:
        return None

    def c_derivatives(branch):
        sl, sh, dp, dq, name = branch
        if name in ("half", "one"):
            return ZERO, ZERO
        hp = hprime_encl(sl, sh)
        if hp is None:
            return None
        return hp * dp, hp * dq

    cd1 = c_derivatives(branch1)
    cd2 = c_derivatives(branch2)
    if cd1 is None or cd2 is None:
        return None

    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    B = w * (ONE - m1) + v * (ONE - m2)
    l1 = (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
    l2 = (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO
    L = w * l1 + v * l2
    c1 = h_encl(branch1[0], branch1[1])
    c2 = h_encl(branch2[0], branch2[1])
    u1 = (sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    u2 = (sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    S = w * u1 + v * u2
    kappa = TWO * (ONE - ALPHA) * B - ONE

    g1 = (-(ONE - ALPHA) * w * L + kappa * w * hps[0] / TWO
          + ALPHA * w * cd1[0] - (ONE - ALPHA) * w * S * gps[0])
    g2 = (-(ONE - ALPHA) * w * L + kappa * w * hps[1] / TWO
          + ALPHA * w * cd1[1] - (ONE - ALPHA) * w * S * gps[1])
    g3 = (-(ONE - ALPHA) * v * L + kappa * v * hps[2] / TWO
          + ALPHA * v * cd2[0] - (ONE - ALPHA) * v * S * gps[2])
    g4 = (-(ONE - ALPHA) * v * L + kappa * v * hps[3] / TWO
          + ALPHA * v * cd2[1] - (ONE - ALPHA) * v * S * gps[3])
    gw = (TWO * (ONE - ALPHA) * (m2 - m1) * L
          + kappa * (l1 - l2) + ALPHA * (c1 - c2)
          - TWO * (ONE - ALPHA) * S * (u1 - u2))
    gradient = (g1, g2, g3, g4, gw)
    if not all(component.is_finite() for component in gradient):
        return None
    return gradient


def _phi_lam_point(point, t, lam=ZERO):
    """Direct Arb point evaluation without p/q canonicalization."""
    p1, q1, p2, q2, w = map(_arb, point)
    t, lam = _arb(t), _arb(lam)
    v = ONE - w
    m1, m2 = (p1 + q1) / TWO, (p2 + q2) / TWO
    l1 = (h_pt(p1) + h_pt(q1)) / TWO
    l2 = (h_pt(p2) + h_pt(q2)) / TWO
    c1, c2 = h_pt(sstar_pt(p1, q1)), h_pt(sstar_pt(p2, q2))
    u1 = (sqrt_rh_enc(p1, p1, RH_GMAX)
          + sqrt_rh_enc(q1, q1, RH_GMAX)) / TWO
    u2 = (sqrt_rh_enc(p2, p2, RH_GMAX)
          + sqrt_rh_enc(q2, q2, RH_GMAX)) / TWO
    mean = w * m1 + v * m2
    entropy = w * l1 + v * l2
    correction = w * c1 + v * c2
    sigma = w * u1 + v * u2
    coeff = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    value = (coeff * entropy + ALPHA * correction
             - (ONE - ALPHA) * sigma * sigma + lam * (mean - t))
    return value if value.is_finite() else None


def centered_lower(box, t):
    """Mean-value lower bound Phi(mid)-sum_i sup|d_i Phi|*radius_i."""
    prepared = _prepare_box(box)
    if prepared is None:
        return None
    gradient = interval_gradient(box)
    if gradient is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    intervals = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh))
    mids = tuple((lo + hi) / TWO for lo, hi in intervals)
    radii = tuple((hi - lo) / TWO for lo, hi in intervals)
    # Keep the true geometric midpoint.  Sorting reversed p/q midpoints can
    # move an anchor outside one coordinate interval while retaining the old
    # radius, invalidating the MVT penalty.
    point_phi = _phi_lam_point(mids, t, ZERO)
    if point_phi is None:
        return None
    penalty = arb(0)
    for derivative, radius in zip(gradient, radii):
        abs_upper = amax(-derivative.lower(), derivative.upper())
        penalty += abs_upper * radius
    lower = _arb(point_phi.lower()) - penalty
    return lower if lower.is_finite() else None


def phi_corner_lower(box, t):
    """Feasible-restricted corner lower bound from component extrema."""
    prepared = _prepare_box(box)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    t = _arb(t)
    w = hull(wl, wh)
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    if mean.lower() > t:
        return None
    m_star = amin(mean.upper(), t)
    coeff_star = TWO * (ONE - ALPHA) * (ONE - m_star) - ONE
    if not (coeff_star > ZERO):
        return None
    L = (w * (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
         + (ONE - w) * (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO)
    L_lo = amax(L.lower(), ZERO)
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    C = (w * h_encl(s1.lower(), s1.upper())
         + (ONE - w) * h_encl(s2.lower(), s2.upper()))
    C_lo = amax(C.lower(), ZERO)
    r1 = (sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    r2 = (sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    sig = w * r1 + (ONE - w) * r2
    S_hi = sig.upper()
    return (coeff_star * L_lo + ALPHA * C_lo
            - (ONE - ALPHA) * S_hi * S_hi)


def _lower(enclosure):
    if enclosure is None:
        return None
    return _arb(enclosure.lower())


def _best(values):
    assert values
    best = values[0]
    for value in values[1:]:
        if value > best:
            best = value
    assert all(best >= value for value in values)
    return best


def collect_bounds(box, t):
    """Evaluate individual bounds and both composites on one box."""
    values = {}
    feasible = _lower(phi_feasible_encl(box, t))
    if feasible is not None:
        values["feas"] = feasible
    lambda_names = []
    for lam_text in LAMBDAS:
        name = "l" + lam_text
        value = _lower(phi_lam_encl(box, t, lam_text))
        if value is not None:
            values[name] = value
            lambda_names.append(name)
    if lambda_names:
        values["lammax"] = _best([values[name] for name in lambda_names])
    shift_sources = (["feas"] if "feas" in values else []) \
        + (["lammax"] if "lammax" in values else [])
    if shift_sources:
        values["shift"] = _best([values[name] for name in shift_sources])
    centered = _lower(centered_lower(box, t))
    if centered is not None:
        values["center"] = centered
    corner = _lower(phi_corner_lower(box, t))
    if corner is not None:
        values["corner"] = corner
    all_sources = [name for name in ("shift", "center", "corner") if name in values]
    if all_sources:
        values["composite"] = _best([values[name] for name in all_sources])
    return values


def _derive_symbolic_gradient():
    """Use SymPy to differentiate Phi for arbitrary affine sstar branches."""
    p1, q1, p2, q2, w, alpha = sp.symbols("p1 q1 p2 q2 w alpha", real=True)
    c1, d1p, d1q, c2, d2p, d2q = sp.symbols(
        "c1 d1p d1q c2 d2p d2q", real=True)

    def hs(z):
        return -(z * sp.log(z) + (1 - z) * sp.log(1 - z)) / sp.log(2)

    def rs(z):
        return 2 * (1 - z) * hs(z) - hs((1 - z) ** 2)

    s1 = c1 + d1p * p1 + d1q * q1
    s2 = c2 + d2p * p2 + d2q * q2
    m1, m2 = (p1 + q1) / 2, (p2 + q2) / 2
    M = w * m1 + (1 - w) * m2
    L = (w * (hs(p1) + hs(q1)) / 2
         + (1 - w) * (hs(p2) + hs(q2)) / 2)
    C = w * hs(s1) + (1 - w) * hs(s2)
    S = (w * (sp.sqrt(rs(p1)) + sp.sqrt(rs(q1))) / 2
         + (1 - w) * (sp.sqrt(rs(p2)) + sp.sqrt(rs(q2))) / 2)
    phi = (2 * (1 - alpha) * (1 - M) - 1) * L + alpha * C - (1 - alpha) * S ** 2
    variables = (p1, q1, p2, q2, w)
    partials = tuple(sp.diff(phi, variable) for variable in variables)
    arguments = variables + (alpha, c1, d1p, d1q, c2, d2p, d2q)
    function = sp.lambdify(arguments, partials, modules="mpmath", cse=True)
    assert len(partials) == 5
    assert all(expression.has(sp.Derivative) is False for expression in partials)
    return function, partials


def _point_sstar_branch(p, q):
    half = mpf("0.5")
    if half >= p and half >= q:
        A = (half, (mpf("0.5"), mpf(0), mpf(0), "half"))
    elif p >= q:
        A = (p, (mpf(0), mpf(1), mpf(0), "p"))
    else:
        A = (q, (mpf(0), mpf(0), mpf(1), "q"))
    if p + q <= 1:
        D = (p + q, (mpf(0), mpf(1), mpf(1), "sum"))
    else:
        D = (mpf(1), (mpf(1), mpf(0), mpf(0), "one"))
    return A[1] if A[0] <= D[0] else D[1]


def _branch_signature(point):
    b1 = _point_sstar_branch(point[0], point[1])
    b2 = _point_sstar_branch(point[2], point[3])
    return b1[3], b2[3]


def verify_symbolic_gradient(symbolic_gradient):
    """Numerically compare the SymPy partials with central differences."""
    old_dps = mp.dps
    mp.dps = 60
    try:
        rng = np.random.default_rng(20260812)
        eps = mpf("1e-7")
        max_relative_error = mpf(0)
        accepted = 0
        checked = 0
        attempts = 0
        while accepted < 100:
            attempts += 1
            assert attempts < 10000
            point = [mpf(str(x)) for x in rng.uniform(0.03, 0.97, 4)]
            point.append(mpf(str(rng.uniform(0.05, 0.95))))
            signature = _branch_signature(point)
            stable = True
            perturbed = []
            for index in range(5):
                plus, minus = list(point), list(point)
                plus[index] += eps
                minus[index] -= eps
                if (_branch_signature(plus) != signature
                        or _branch_signature(minus) != signature):
                    stable = False
                    break
                perturbed.append((plus, minus))
            if not stable:
                continue
            b1 = _point_sstar_branch(point[0], point[1])
            b2 = _point_sstar_branch(point[2], point[3])
            args = tuple(point) + (mpf("0.0356069"),) + b1[:3] + b2[:3]
            symbolic = tuple(mpf(value) for value in symbolic_gradient(*args))
            local_errors = []
            for index, (plus, minus) in enumerate(perturbed):
                finite_difference = (phi_true(*plus) - phi_true(*minus)) / (2 * eps)
                scale = max(abs(symbolic[index]), abs(finite_difference), mpf("1e-12"))
                relative_error = abs(symbolic[index] - finite_difference) / scale
                local_errors.append(relative_error)
            if max(local_errors) >= mpf("1e-6"):
                raise AssertionError((point, symbolic, local_errors))
            max_relative_error = max(max_relative_error, max(local_errors))
            accepted += 1
            checked += 5
        assert accepted == 100
        assert checked == 500
        assert max_relative_error < mpf("1e-6")
        return accepted, checked, max_relative_error
    finally:
        mp.dps = old_dps


def _mean_float(point):
    p1, q1, p2, q2, w = point
    return w * (p1 + q1) / 2 + (1 - w) * (p2 + q2) / 2


def _random_box_and_samples(rng, t):
    for _ in range(10000):
        pair1 = sorted(rng.uniform(0.01, 0.72, 2))
        pair2 = sorted(rng.uniform(0.01, 0.72, 2))
        if pair1[1] - pair1[0] < 0.04 or pair2[1] - pair2[0] < 0.04:
            continue
        center = (pair1[0], pair1[1], pair2[0], pair2[1], rng.uniform(0.05, 0.95))
        if _mean_float(center) >= t - 0.02:
            continue
        gap1, gap2 = pair1[1] - pair1[0], pair2[1] - pair2[0]
        radius_caps = (
            min(center[0] - 0.002, 0.998 - center[0], gap1 / 4, 0.05),
            min(center[1] - 0.002, 0.998 - center[1], gap1 / 4, 0.05),
            min(center[2] - 0.002, 0.998 - center[2], gap2 / 4, 0.05),
            min(center[3] - 0.002, 0.998 - center[3], gap2 / 4, 0.05),
            min(center[4] - 0.002, 0.998 - center[4], 0.06),
        )
        if min(radius_caps) <= 1e-4:
            continue
        radii = [cap * rng.uniform(0.08, 0.95) for cap in radius_caps]
        box = tuple((value - radius, value + radius)
                    for value, radius in zip(center, radii))
        samples = []
        for _ in range(10000):
            point = tuple(rng.uniform(lo + (hi - lo) * 1e-6,
                                      hi - (hi - lo) * 1e-6)
                          for lo, hi in box)
            if (point[0] < point[1] and point[2] < point[3]
                    and _mean_float(point) < t):
                samples.append(point)
                if len(samples) == 20:
                    return box, samples
    raise AssertionError("could not generate a sampled-soundness box")


def sampled_soundness(t):
    """Check every returned bound on 150 boxes and 20 feasible samples each."""
    rng = np.random.default_rng(31415926)
    violations = []
    boxes = 0
    samples_checked = 0
    bound_checks = 0
    centered_boxes = 0
    for _ in range(150):
        box, samples = _random_box_and_samples(rng, t)
        values = collect_bounds(box, t)
        assert "feas" in values and "shift" in values and "corner" in values
        assert "composite" in values
        if "center" in values:
            centered_boxes += 1
        boxes += 1
        assert len(samples) == 20
        for point in samples:
            assert all(lo < x < hi for x, (lo, hi) in zip(point, box))
            assert point[0] < point[1] and point[2] < point[3]
            mean = _mean_float(point)
            assert mean < t
            for lam_text in LAMBDAS:
                assert float(lam_text) * (mean - t) <= 0
            truth = phi_true(*point)
            rhs = _arb(mp.nstr(truth, 45)) + arb("1e-20")
            for name, lower in values.items():
                bound_checks += 1
                if not (lower <= rhs):
                    violations.append((boxes, samples_checked, name, lower, truth))
            samples_checked += 1
    assert boxes == 150
    assert samples_checked == 3000
    assert bound_checks > 0
    assert 0 <= centered_boxes <= boxes
    assert len(violations) == 0, violations[:3]
    return boxes, samples_checked, bound_checks, centered_boxes, len(violations)


def full_width_box(center, width):
    """A box in [0,1]^5 whose five actual full widths equal width."""
    assert 0 < width <= 1
    box = []
    for value in center:
        lo = value - width / 2
        if lo < 0:
            lo = 0.0
        if lo + width > 1:
            lo = 1.0 - width
        hi = lo + width
        assert 0 <= lo < hi <= 1
        assert abs((hi - lo) - width) <= 4e-16
        assert lo <= value <= hi
        box.append((lo, hi))
    return tuple(box)


def _format_lower(value):
    if value is None:
        return "       n/a"
    assert value.is_finite()
    return "%+10.3e" % float(value.lower())


def probe_table(label, center, t, widths):
    boxes = [full_width_box(center, width) for width in widths]
    assert all(all(abs((hi - lo) - width) <= 4e-16 for lo, hi in box)
               for width, box in zip(widths, boxes))
    print("\n%s; d is the actual FULL width in every coordinate" % label)
    columns = ("d", "feas", "l0", "l0.5", "l1.0", "l1.626", "l2.2", "l3.0",
               "shift", "center", "corner", "composite")
    assert len(columns) == 12
    print(" ".join("%10s" % column for column in columns))
    rows = []
    for width, box in zip(widths, boxes):
        values = collect_bounds(box, t)
        assert "composite" in values and "shift" in values and "corner" in values
        source_names = [name for name in ("shift", "center", "corner") if name in values]
        assert values["composite"] >= _best([values[name] for name in source_names])
        assert all(value.is_finite() for value in values.values())
        fields = ["%10.1e" % width]
        for name in columns[1:]:
            fields.append(_format_lower(values.get(name)))
        assert len(fields) == len(columns)
        print(" ".join(fields))
        rows.append((width, box, values))
    assert len(rows) == len(widths)
    return rows


def _bound_at_width(center, t, width, name):
    values = collect_bounds(full_width_box(center, width), t)
    assert name in values
    return float(values[name].lower())


def numerical_flip_bracket(center, t, name):
    """Bisect the observed positive-to-negative lower-bound sign transition."""
    positive_width = 1e-8
    negative_width = 1e-2
    positive_value = _bound_at_width(center, t, positive_width, name)
    negative_value = _bound_at_width(center, t, negative_width, name)
    assert positive_value >= 0
    assert negative_value < 0
    for _ in range(42):
        midpoint = (positive_width + negative_width) / 2
        value = _bound_at_width(center, t, midpoint, name)
        if value >= 0:
            positive_width, positive_value = midpoint, value
        else:
            negative_width, negative_value = midpoint, value
    assert positive_value >= 0
    assert negative_value < 0
    assert positive_width < negative_width
    assert negative_width - positive_width < 3e-15
    return positive_width, negative_width


def main():
    global RH_GMAX
    started = time.monotonic()

    RH_GMAX = _rh_global()
    assert _arb("0.2342294") <= RH_GMAX <= _arb("0.2350")
    assert _arb(RH_GMAX).is_finite()
    print("Arb helper PASS: certified global rh upper = %.12f" % float(RH_GMAX))

    symbolic_gradient, partials = _derive_symbolic_gradient()
    assert len(partials) == 5
    point_count, derivative_count, max_error = verify_symbolic_gradient(symbolic_gradient)
    assert point_count == 100 and derivative_count == 500
    assert max_error < mpf("1e-6")
    print("Gradient check PASS: 100 random interior points, 500 partials, max relative error %.3e"
          % float(max_error))

    t = float(PSI)
    boxes, samples, checks, centered_boxes, violations = sampled_soundness(t)
    assert boxes == 150 and samples == 3000 and checks > 0
    assert 0 <= centered_boxes <= boxes and violations == 0
    print("SAMPLED EVIDENCE ONLY (not a universal proof): 150 random boxes x 20 feasible interior samples, %d bound checks, 0 violations; centered usable on %d boxes"
          % (checks, centered_boxes))

    widths = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
    assert widths == (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
    b = float(B_OBS)
    a = (t - b) / (1 - b)
    w_star = 1 - 2 * a
    obstruction = (b, b, b, 1.0, w_star)
    assert abs(_mean_float(obstruction) - t) < 2e-16
    obstruction_rows = probe_table("OBSTRUCTION (b,b,b,1,w*) at t=psi", obstruction, t, widths)
    obstruction_centered = sum("center" in values for _, _, values in obstruction_rows)
    assert obstruction_centered == 0
    assert all(box[3][1] == 1.0 for _, box, _ in obstruction_rows)
    print("Obstruction centered availability: 0/7 widths (q2 interval touches 1)")

    q_interior = 0.999
    m2 = (b + q_interior) / 2
    w_feasible = (m2 - t) / (m2 - b)
    interior_center = (b, b, b, q_interior, w_feasible)
    assert abs(_mean_float(interior_center) - t) < 2e-16
    interior_rows = probe_table("INTERIOR q2=0.999 with mean-binding w_feasible", interior_center,
                                t, widths)
    usable_widths = [width for width, _, values in interior_rows if "center" in values]
    positive_widths = [width for width, _, values in interior_rows
                       if "center" in values and values["center"] >= ZERO]
    winning_widths = [width for width, _, values in interior_rows
                      if "center" in values
                      and values["center"] >= values["shift"]
                      and values["center"] >= values["corner"]]
    assert all(width in widths for width in usable_widths)
    assert all(width in usable_widths for width in positive_widths)
    assert all(width in usable_widths for width in winning_widths)
    print("Interior centered measurements: usable d=%s; nonnegative d=%s; composite-winning d=%s"
          % (usable_widths, positive_widths, winning_widths))

    feasible_bracket = numerical_flip_bracket(obstruction, t, "feas")
    lambda_bracket = numerical_flip_bracket(obstruction, t, "lammax")
    shift_bracket = numerical_flip_bracket(obstruction, t, "shift")
    corner_bracket = numerical_flip_bracket(obstruction, t, "corner")
    composite_bracket = numerical_flip_bracket(obstruction, t, "composite")
    brackets = (feasible_bracket, lambda_bracket, shift_bracket,
                corner_bracket, composite_bracket)
    for lo, hi in brackets:
        assert lo < hi and hi - lo < 3e-15
    old_flip = 2e-4
    feasible_mid = sum(feasible_bracket) / 2
    lambda_mid = sum(lambda_bracket) / 2
    shift_mid = sum(shift_bracket) / 2
    corner_mid = sum(corner_bracket) / 2
    composite_mid = sum(composite_bracket) / 2
    ratio = composite_mid / old_flip
    assert abs(ratio * old_flip - composite_mid) < 1e-18
    print("NUMERICAL flip brackets in actual FULL width: feasible=[%.12g, %.12g], lambda-max=[%.12g, %.12g], feasible+lambda=[%.12g, %.12g], corner=[%.12g, %.12g], composite=[%.12g, %.12g]"
          % (feasible_bracket + lambda_bracket + shift_bracket
             + corner_bracket + composite_bracket))
    print("NUMERICAL composite/old-2e-4 flip-width ratio = %.6f" % ratio)

    # This isolates the incremental effect of adding the lambda family to the
    # pre-existing feasibility enclosure, rather than crediting that enclosure
    # (or the corner bound) to lambda shifting.
    lambda_gain = shift_mid - feasible_mid
    assert abs(lambda_gain) < 3e-15
    assert lambda_mid < feasible_mid
    print("NUMERICAL lambda-shift answer: no obstruction-width gain over the feasibility enclosure; delta %.3e, while lambda-only max flips at %.12g"
          % (lambda_gain, lambda_mid))
    assert composite_mid > old_flip
    print("NUMERICAL overall composite versus old 2e-4: factor %.6f (delta %.6e)"
          % (ratio, composite_mid - old_flip))
    if shift_mid > corner_mid:
        assert shift_mid / corner_mid > 1
        print("NUMERICAL feasibility+lambda/corner flip-width factor = %.6f"
              % (shift_mid / corner_mid))
    else:
        assert shift_mid / corner_mid <= 1
        print("NUMERICAL feasibility+lambda does not beat the reimplemented corner flip; factor %.6f"
              % (shift_mid / corner_mid))

    elapsed = time.monotonic() - started
    assert elapsed < 600
    print("Runtime PASS: %.2f seconds < 600 seconds" % elapsed)


if __name__ == "__main__":
    main()
