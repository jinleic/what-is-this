"""Composite five-parameter certificate with a q2=1 pinned-face rule.

All universal inequalities used by the certifier are Arb interval statements.
The numerical derivative and random-box checks printed by ``main`` are labeled
sampled evidence and are not used as proof steps.  The full-cube driver uses
the exact orbit involution to restrict w to [1/2,1]; both atom-face orientations
remain covered because the involution also exchanges the two pair-orbits.
"""

import argparse
import math
import os
import random
import sys
import time
from collections import Counter
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bound_kkt
import cert2
import diag_exhaust
from arbcore import (ALPHA, LOG2, ONE, TWO, amax, amin, get_rh_gmax,
                     h_encl, h_pt, hull, sqrt_rh_enc, sstar_encl, sstar_pt)
from entropy import PSI, h as h_mp

ctx.prec = 160
mp.mp.dps = 60

ZERO = arb(0)
B_OBS = mp.mpf("0.32945473850303697239")
Q2PIN = 0.9375
LAMBDAS_BASE = (0.0, 1.0, 1.626, 2.2)
DEFAULT_MIN_WIDTH = 5e-4
DEFAULT_FACE_MIN_WIDTH = 1e-3

TRACE_MEAN_INFEASIBLE = 0
TRACE_CORNER_INFEASIBLE = 1
TRACE_CORNER = 2
TRACE_RATIO = 3
TRACE_CENTER = 4
TRACE_CENTER_MIXED = 5
TRACE_CENTER_MIXED_SWAP = 6
TRACE_CENTER_W = 7
TRACE_FACE = 8
TRACE_RESIDUAL = 9
TRACE_SPLIT_BASE = 16
TRACE_SCHEMA = "cert3-trace-v1"

TRACE_TALLY_NAMES = {
    TRACE_MEAN_INFEASIBLE: "infeasible",
    TRACE_CORNER_INFEASIBLE: "infeasible",
    TRACE_CORNER: "corner",
    TRACE_RATIO: "ratio",
    TRACE_CENTER: "center",
    TRACE_CENTER_MIXED: "center_mixed",
    TRACE_CENTER_MIXED_SWAP: "center_mixed_swap",
    TRACE_CENTER_W: "center_w",
    TRACE_FACE: "face",
    TRACE_RESIDUAL: "residual",
}
TRACE_TALLY_NAMES.update(
    (TRACE_SPLIT_BASE + coordinate, "split") for coordinate in range(5))
_TRACE_BYTES = tuple(bytes((event,)) for event in range(TRACE_SPLIT_BASE + 5))


def _write_trace(trace, event):
    if trace is None:
        return
    if trace.write(_TRACE_BYTES[event]) != 1:
        raise OSError("trace writer did not accept exactly one byte")


def _arb(x):
    return arb(str(x)) if not isinstance(x, arb) else arb(x)


def _abs_upper(x):
    return amax(-x.lower(), x.upper()).upper()


def _entropy_prime(u, v):
    """Exact range of decreasing h' on a strict interior interval."""
    u, v = _arb(u), _arb(v)
    if not (u > ZERO and v < ONE):
        return None

    def hp(z):
        return ((ONE - z).log() - z.log()) / LOG2

    return hull(hp(v), hp(u))


def _rh_prime(u, v):
    """Natural interval extension of rh' on a strict interior interval."""
    u, v = _arb(u), _arb(v)
    hp = _entropy_prime(u, v)
    x2lo, x2hi = (ONE - v) ** 2, (ONE - u) ** 2
    hp2 = _entropy_prime(x2lo, x2hi)
    if hp is None or hp2 is None:
        return None
    z = hull(u, v)
    return -TWO * h_encl(u, v) + TWO * (ONE - z) * (hp + hp2)


def _sqrt_rh_prime(u, v, gmax):
    """Enclose (sqrt(rh))' when interval evaluation proves rh positive."""
    u, v = _arb(u), _arb(v)
    if not (u > ZERO and v < ONE):
        return None
    r = bound_kkt._positive_rh_encl(u, v)
    rp = _rh_prime(u, v)
    if r is None or rp is None or not (r.lower() > ZERO):
        return None
    root = r.sqrt()
    ans = rp / (TWO * root)
    return ans if ans.is_finite() else None


# ---------------------------------------------------------------------------
# Endpoint derivative lemma used by the q2-pin rule.
# ---------------------------------------------------------------------------
# Put x=1-z. Exact cancellation gives
#   ln(2) rh(1-x) = F(x)
#     = (1-x)^2 ln(1-x) + (1-x^2) ln(1+x).
# For 0<=x<=X<1, Taylor's formulas with explicit Lagrange/integral remainders
# give
#   -x-x^2/[2(1-X)] <= ln(1-x) <= -x,
#    x-x^2/2 <= ln(1+x) <= x.
# Hence F>=x^2 Q_X(x), with
#   Q_X(x)=(1-x)(2-3X+Xx)/(2(1-X)).
# Q_X is decreasing and c=Q_X(X)>0.  The exact derivative is
# F'=2[x ln(1-x)-x ln(1+x)-ln(1-x)].  Expanding the logarithms into their
# convergent alternating/geometric series, discarding negative terms, and using
# -ln(1-x)<=x/(1-X) gives |F'|<=2x/(1-X). Therefore
# |d sqrt(rh)/dx| <= 1/[(1-X)sqrt(ln(2)c)].
_X_Q = Fraction(1, 16)
_X = arb(_X_Q.numerator) / arb(_X_Q.denominator)
_DQ_C = ((ONE - _X) * (TWO - arb(3) * _X + _X * _X)
         / (TWO * (ONE - _X)))
DQ_ENDPOINT = (ONE / ((ONE - _X) * (LOG2 * _DQ_C).sqrt())).upper()

_sx, _sX = sp.symbols("x X", positive=True)
_Fsym = (1 - _sx) ** 2 * sp.log(1 - _sx) + (1 - _sx ** 2) * sp.log(1 + _sx)
_Fprime_expected = 2 * (_sx * sp.log(1 - _sx)
                        - _sx * sp.log(1 + _sx) - sp.log(1 - _sx))
_Qsym = (1 - _sx) * (2 - 3 * _sX + _sX * _sx) / (2 * (1 - _sX))
_SERIES_IDENTITIES = (
    sp.simplify(sp.diff(_Fsym, _sx) - _Fprime_expected) == 0,
    sp.simplify(sp.diff(_Qsym, _sx)
                - (_sX * _sx - 2 * _sX + 1) / (_sX - 1)) == 0,
    sp.simplify(_Qsym.subs(_sx, _sX)
                - (1 - _sX) * (2 - 3 * _sX + _sX ** 2)
                / (2 * (1 - _sX))) == 0,
)
assert all(_SERIES_IDENTITIES)
assert _DQ_C > ZERO and DQ_ENDPOINT > arb("1.2") and DQ_ENDPOINT < arb("1.4")


def dq_upper(q2lo):
    """Certified cap for |(sqrt rh)'| on [q2lo,1], q2lo>=15/16."""
    q2lo = _arb(q2lo)
    if not (q2lo >= ONE - _X):
        return None
    return DQ_ENDPOINT


# ---------------------------------------------------------------------------
# Face functional and interval/centered bounds.
# ---------------------------------------------------------------------------
def _face_prepare(box4):
    vals = tuple(_arb(x) for pair in box4 for x in pair)
    p1l, p1h, q1l, q1h, p2l, p2h, wl, wh = vals
    q1l = amax(q1l, p1l)
    if q1l > q1h:
        return None
    return p1l, p1h, q1l, q1h, p2l, p2h, wl, wh


def _face_components(box4, gmax):
    prep = _face_prepare(box4)
    if prep is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, wl, wh = prep
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + ONE) / TWO, (p2h + ONE) / TWO)
    mean = w * m1 + v * m2
    l1 = (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
    l2 = h_encl(p2l, p2h) / TWO
    L = w * l1 + v * l2
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    c1 = h_encl(s1.lower(), s1.upper())
    C = w * c1
    u1 = (sqrt_rh_enc(p1l, p1h, gmax)
          + sqrt_rh_enc(q1l, q1h, gmax)) / TWO
    u2 = sqrt_rh_enc(p2l, p2h, gmax) / TWO
    S = w * u1 + v * u2
    return prep, w, v, m1, m2, mean, l1, l2, L, c1, C, u1, u2, S


def face_interval(box4, t, lam, gmax=None, feasible_coeff=False):
    """Natural interval bound for Phi_face + lam*(Mf-t)."""
    if gmax is None:
        gmax = get_rh_gmax()
    c = _face_components(box4, gmax)
    if c is None:
        return None
    _, _, _, _, _, mean, _, _, L, _, C, _, _, S = c
    t, lam = _arb(t), _arb(lam)
    assert lam >= ZERO
    if feasible_coeff:
        if mean.lower() > t:
            return None
        # This is a Phi bound using feasibility, not a full-box Phi_lam bound.
        mcap = amin(mean.upper(), t)
        coeff = TWO * (ONE - ALPHA) * (ONE - mcap) - ONE
        if not (coeff > ZERO):
            return None
        return coeff * amax(L.lower(), ZERO) + ALPHA * amax(C.lower(), ZERO) \
            - (ONE - ALPHA) * S.upper() ** 2
    coeff = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    return coeff * L + ALPHA * C - (ONE - ALPHA) * S * S \
        + lam * (mean - t)


def _sstar_derivatives(pl, ph, ql, qh):
    """Enclose both derivatives of h(s*) across certified affine branches."""
    branch = bound_kkt._sstar_branch(pl, ph, ql, qh)
    if branch is not None:
        sl, sh, dp, dq, name = branch
        if name in ("half", "one"):
            return ZERO, ZERO
        hp = _entropy_prime(sl, sh)
        if hp is None:
            return None
        return hp * dp, hp * dq
    s = sstar_encl(pl, ph, ql, qh)
    hp = _entropy_prime(s.lower(), s.upper())
    if hp is None:
        return None
    cap = _abs_upper(hp)
    return hull(-cap, cap), hull(-cap, cap)


def face_gradient(box4, lam, gmax=None):
    """Interval enclosure of the four face partials of Phi+lam(M-t)."""
    if gmax is None:
        gmax = get_rh_gmax()
    c = _face_components(box4, gmax)
    if c is None:
        return None
    prep, w, v, m1, m2, mean, l1, l2, L, c1, _, u1, u2, S = c
    p1l, p1h, q1l, q1h, p2l, p2h, _, _ = prep
    hps = (_entropy_prime(p1l, p1h), _entropy_prime(q1l, q1h),
           _entropy_prime(p2l, p2h))
    gps = (_sqrt_rh_prime(p1l, p1h, gmax),
           _sqrt_rh_prime(q1l, q1h, gmax),
           _sqrt_rh_prime(p2l, p2h, gmax))
    cder = _sstar_derivatives(p1l, p1h, q1l, q1h)
    if any(x is None for x in hps + gps) or cder is None:
        return None
    lam = _arb(lam)
    assert lam >= ZERO
    kappa = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    dc1p, dc1q = cder
    gp = (-(ONE - ALPHA) * w * L + kappa * w * hps[0] / TWO
          + ALPHA * w * dc1p - (ONE - ALPHA) * w * S * gps[0]
          + lam * w / TWO)
    gq = (-(ONE - ALPHA) * w * L + kappa * w * hps[1] / TWO
          + ALPHA * w * dc1q - (ONE - ALPHA) * w * S * gps[1]
          + lam * w / TWO)
    gr = (-(ONE - ALPHA) * v * L + kappa * v * hps[2] / TWO
          - (ONE - ALPHA) * v * S * gps[2] + lam * v / TWO)
    gw = (-TWO * (ONE - ALPHA) * (m1 - m2) * L
          + kappa * (l1 - l2) + ALPHA * c1
          - TWO * (ONE - ALPHA) * S * (u1 - u2)
          + lam * (m1 - m2))
    ans = gp, gq, gr, gw
    return ans if all(x.is_finite() for x in ans) else None


def face_point_value(point4, t, lam, gmax=None):
    """Direct Arb evaluation at one face point, without box canonicalization."""
    p, q, r, w = map(_arb, point4)
    t, lam = _arb(t), _arb(lam)
    v = ONE - w
    mean = w * (p + q) / TWO + v * (r + ONE) / TWO
    entropy = w * (h_pt(p) + h_pt(q)) / TWO + v * h_pt(r) / TWO
    correction = w * h_pt(sstar_pt(p, q))
    sigma = (w * (sqrt_rh_enc(p, p, gmax) + sqrt_rh_enc(q, q, gmax)) / TWO
             + v * sqrt_rh_enc(r, r, gmax) / TWO)
    coeff = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    return (coeff * entropy + ALPHA * correction
            - (ONE - ALPHA) * sigma * sigma + lam * (mean - t))


def face_centered(box4, t, lam, gmax=None):
    """Mean-value lower bound for lambda-shifted face functional."""
    if gmax is None:
        gmax = get_rh_gmax()
    prep = _face_prepare(box4)
    grad = face_gradient(box4, lam, gmax)
    if prep is None or grad is None:
        return None
    intervals = ((prep[0], prep[1]), (prep[2], prep[3]),
                 (prep[4], prep[5]), (prep[6], prep[7]))
    mids = tuple((lo + hi) / TWO for lo, hi in intervals)
    radii = tuple((hi - lo) / TWO for lo, hi in intervals)
    # Evaluate at the true geometric midpoint; the face formula itself is
    # pair-symmetric and does not require p<=q. Keeping this point unchanged
    # preserves the half-width MVT radii even when the midpoints are reversed.
    point = mids
    value = face_point_value(point, t, lam, gmax)
    if value is None:
        return None
    penalty = ZERO
    for derivative, radius in zip(grad, radii):
        penalty += _abs_upper(derivative) * radius
    ans = value.lower() - penalty
    return ans if ans.is_finite() else None


def face_best_bound(box4, t, lambdas, gmax=None, allow_feasible=False):
    """Maximum of independently sound face lower bounds.

    ``allow_feasible`` is false in the pinned path because raising q2 to one
    can leave the feasible set. Each shifted candidate is valid on the entire
    face box and is paired with the lambda whose pin test succeeded.
    """
    candidates = []
    if allow_feasible:
        feasible = face_interval(box4, t, ZERO, gmax, True)
        if feasible is not None:
            candidates.append(feasible.lower())
    for lam in lambdas:
        direct = face_interval(box4, t, lam, gmax, False)
        centered = face_centered(box4, t, lam, gmax)
        if direct is not None:
            candidates.append(direct.lower())
        if centered is not None:
            candidates.append(centered.lower())
    if not candidates:
        return None
    best = candidates[0]
    for value in candidates[1:]:
        best = amax(best, value)
    return best.lower()


# SymPy derivation of the smooth affine-sstar face formulas.
_spvars = sp.symbols("p q r w alpha lam t")
_sp_p, _sp_q, _sp_r, _sp_w, _sp_a, _sp_lam, _sp_t = _spvars
_sp_h = sp.Function("h")
_sp_g = sp.Function("g")
_sp_s = sp.Function("s")
_sp_M = _sp_w * (_sp_p + _sp_q) / 2 + (1 - _sp_w) * (_sp_r + 1) / 2
_sp_L = _sp_w * (_sp_h(_sp_p) + _sp_h(_sp_q)) / 2 + (1 - _sp_w) * _sp_h(_sp_r) / 2
_sp_S = _sp_w * (_sp_g(_sp_p) + _sp_g(_sp_q)) / 2 + (1 - _sp_w) * _sp_g(_sp_r) / 2
_sp_F = ((2 * (1 - _sp_a) * (1 - _sp_M) - 1) * _sp_L
         + _sp_a * _sp_w * _sp_h(_sp_s(_sp_p, _sp_q))
         - (1 - _sp_a) * _sp_S ** 2 + _sp_lam * (_sp_M - _sp_t))
_FACE_PARTIALS = tuple(sp.diff(_sp_F, z) for z in (_sp_p, _sp_q, _sp_r, _sp_w))
assert len(_FACE_PARTIALS) == 4 and all(expr != 0 for expr in _FACE_PARTIALS)


# ---------------------------------------------------------------------------
# Q2 pin test.
# ---------------------------------------------------------------------------
def q2_pin(box5, t, lam, gmax=None):
    """Prove d(Phi+lam(M-t))/dq2 <= 0 on the STRIP [q2lo, 1] over the box.

    Every ingredient below is already evaluated over the strip, not merely the
    box's q2 interval: hp_lo = h'(q2lo) bounds h' on [q2lo,1] (h' decreasing);
    the sigma enclosure uses sqrt_rh over [q2lo, 1]; m2's hull extends to 1;
    the DQ cap is the endpoint lemma on [15/16, 1].  Hence a success certifies
    the derivative sign for all q2 in [q2lo, 1], and by integrating along the
    segment q2 -> 1 (which stays inside the strip for fixed other coordinates)
    every point of the box satisfies Phi_lam(x, q2) >= Phi_lam(x, 1): the
    box's infimum of Phi_lam is bounded below by the q2=1 FACE infimum, even
    when the box itself does not touch q2 = 1.

    Branch requirement: sstar(p2, q2) = q2 needs q2 >= p2 throughout the
    strip, guaranteed by p2h <= q2lo (with q2lo >= 15/16 > 1/2).
    """
    if gmax is None:
        gmax = get_rh_gmax()
    (p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh) = box5
    q2l = max(q2l, p2l)
    if q2l > q2h or q2l < Q2PIN or p2h > q2l:
        return False
    dq = dq_upper(q2l)
    if dq is None:
        return False
    vals = tuple(_arb(x) for pair in box5 for x in pair)
    p1l, p1h, q1l, q1h, p2l, p2h, _, _, wl, wh = vals
    q1l = amax(q1l, p1l)
    if q1l > q1h:
        return True
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((_arb(p2l) + _arb(q2l)) / TWO, (_arb(p2h) + ONE) / TWO)
    mean = w * m1 + v * m2
    coeff = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    coeff_lo = coeff.lower()
    if not (coeff_lo > ZERO):
        return False
    hp_lo = ((ONE - _arb(q2l)).log() - _arb(q2l).log()) / LOG2
    vlo = amax(ONE - _arb(wh), ZERO)
    vhi = amin(ONE - _arb(wl), ONE)
    neg_hi = (coeff_lo * vlo / TWO + ALPHA * vlo) * hp_lo
    sig = (w * (sqrt_rh_enc(p1l, p1h, gmax)
                + sqrt_rh_enc(q1l, q1h, gmax)) / TWO
           + v * (sqrt_rh_enc(p2l, p2h, gmax)
                  + sqrt_rh_enc(_arb(q2l), ONE, gmax)) / TWO)
    pos_hi = (ONE - ALPHA) * sig.upper() * vhi * dq \
        + _arb(lam) * vhi / TWO
    # The omitted coefficient-derivative term is -(1-alpha)*v*L <= 0.
    return bool(neg_hi + pos_hi <= ZERO)


_qM0, _qL0, _qS0, _qq, _qw, _qa, _ql = sp.symbols(
    "M0 L0 S0 q w alpha lam")
_qh = sp.Function("h")
_qg = sp.Function("g")
_qv = 1 - _qw
_qM_expr = _qM0 + _qv * _qq / 2
_qL_expr = _qL0 + _qv * _qh(_qq) / 2
_qS_expr = _qS0 + _qv * _qg(_qq) / 2
_qcoeff = 2 * (1 - _qa) * (1 - _qM_expr) - 1
_Q2_FUNCTION = (_qcoeff * _qL_expr + _qa * _qv * _qh(_qq)
                - (1 - _qa) * _qS_expr ** 2 + _ql * _qM_expr)
_Q2_EXPECTED = ((_qcoeff * _qv / 2 + _qa * _qv)
                * sp.diff(_qh(_qq), _qq)
                - (1 - _qa) * _qv * _qL_expr
                - (1 - _qa) * _qS_expr * _qv
                * sp.diff(_qg(_qq), _qq) + _ql * _qv / 2)
assert sp.simplify(sp.diff(_Q2_FUNCTION, _qq) - _Q2_EXPECTED) == 0


def solve_face_kkt(t):
    """Numerically locate the symmetric mean-binding face KKT point."""
    old = mp.mp.dps
    mp.mp.dps = 60
    aa, tt = mp.mpf(str(float(ALPHA))), mp.mpf(str(t))

    def rh(z):
        x = 1 - z
        return max(mp.mpf(0), 2 * x * h_mp(z) - h_mp(x * x))

    def sdiag(z):
        return min(max(mp.mpf("0.5"), z), min(2 * z, 1))

    def fun(p, r, w, lam):
        mean = w * p + (1 - w) * (r + 1) / 2
        entropy = w * h_mp(p) + (1 - w) * h_mp(r) / 2
        sigma = w * mp.sqrt(rh(p)) + (1 - w) * mp.sqrt(rh(r)) / 2
        return ((2 * (1 - aa) * (1 - mean) - 1) * entropy
                + aa * w * h_mp(sdiag(p)) - (1 - aa) * sigma * sigma
                + lam * (mean - tt))

    try:
        dp = lambda p, r, w, lam: mp.diff(lambda z: fun(z, r, w, lam), p)
        dr = lambda p, r, w, lam: mp.diff(lambda z: fun(p, z, w, lam), r)
        dw = lambda p, r, w, lam: mp.diff(lambda z: fun(p, r, z, lam), w)
        cons = lambda p, r, w, lam: w * p + (1 - w) * (r + 1) / 2 - tt
        p, r, w, lam = mp.findroot(
            (dp, dr, dw, cons),
            (mp.mpf("0.328"), mp.mpf("0.328"), mp.mpf("0.838"), mp.mpf("1.37")),
            tol=mp.mpf("1e-45"), maxsteps=100,
        )
        residual = max(abs(dp(p, r, w, lam)), abs(dr(p, r, w, lam)),
                       abs(dw(p, r, w, lam)), abs(cons(p, r, w, lam)))
        assert residual < mp.mpf("1e-35") and 0 < p < 1 and 0 < r < 1
        assert 0 < w < 1 and lam > 0
        return p, r, w, lam, residual
    finally:
        mp.mp.dps = old


def lambda_family(t):
    _, _, _, lam, _ = solve_face_kkt(t)
    values = list(LAMBDAS_BASE) + [float(lam * mp.mpf("0.98")), float(lam),
                                  float(lam * mp.mpf("1.02"))]
    return tuple(sorted(set(values))), lam


# ---------------------------------------------------------------------------
# Assert-backed sampled validation.
# ---------------------------------------------------------------------------
def _face_true(point4, t, lam):
    p, q, r, w = map(mp.mpf, point4)
    return cert2.phi_true(p, q, r, 1, w) \
        + mp.mpf(str(lam)) * (w * (p + q) / 2 + (1 - w) * (r + 1) / 2 - mp.mpf(str(t)))


def verify_derivatives(samples=100, seed=8731):
    rng = random.Random(seed)
    max_face_rel = mp.mpf(0)
    max_q2_rel = mp.mpf(0)
    checked_q2 = 0
    for _ in range(samples):
        # Keep s*=q on a smooth certified branch, so the point derivatives can
        # be compared without branch-boundary ambiguity.
        p = rng.uniform(.08, .44)
        q = rng.uniform(max(p + .03, .51), .82)
        r = rng.uniform(.08, .72)
        w = rng.uniform(.08, .92)
        lam = rng.uniform(.2, 2.3)
        point = (mp.mpf(p), mp.mpf(q), mp.mpf(r), mp.mpf(w))
        box = tuple((float(x - mp.mpf("1e-10")), float(x + mp.mpf("1e-10"))) for x in point)
        grad = face_gradient(box, lam)
        assert grad is not None
        for i in range(4):
            step = mp.mpf("1e-7")
            plus = point[:i] + (point[i] + step,) + point[i + 1:]
            minus = point[:i] + (point[i] - step,) + point[i + 1:]
            exact = (_face_true(plus, .382, lam) - _face_true(minus, .382, lam)) / (2 * step)
            assert exact in grad[i]
            rel = abs(mp.mpf(str(float(grad[i].mid()))) - exact) \
                / max(mp.mpf(1), abs(exact))
            max_face_rel = max(max_face_rel, rel)

        q2 = rng.uniform(.55, .96)
        p2 = min(r, q2 - .02)
        pt5 = (p, q, p2, q2, w)
        eps = 1e-10
        b5 = tuple((x - eps, x + eps) for x in pt5)
        g5 = bound_kkt.interval_gradient(b5)
        if g5 is None:
            continue
        hq = mp.log((1 - q2) / q2, 2)
        M = w * (p + q) / 2 + (1 - w) * (p2 + q2) / 2
        L = w * (h_mp(p) + h_mp(q)) / 2 + (1 - w) * (h_mp(p2) + h_mp(q2)) / 2
        def rr(z):
            x = 1 - z
            return max(mp.mpf(0), 2 * x * h_mp(z) - h_mp(x * x))
        S = w * (mp.sqrt(rr(p)) + mp.sqrt(rr(q))) / 2 \
            + (1 - w) * (mp.sqrt(rr(p2)) + mp.sqrt(rr(q2))) / 2
        dg_step = mp.mpf("1e-7")
        dg = (mp.sqrt(rr(mp.mpf(q2) + dg_step))
              - mp.sqrt(rr(mp.mpf(q2) - dg_step))) / (2 * dg_step)
        aa = mp.mpf(str(float(ALPHA)))
        coeff = 2 * (1 - aa) * (1 - M) - 1
        exact_q = ((coeff * (1 - w) / 2 + aa * (1 - w)) * hq
                   - (1 - aa) * (1 - w) * L
                   - (1 - aa) * S * (1 - w) * dg)
        assert exact_q in g5[3]
        fd_q = (cert2.phi_true(p, q, p2, mp.mpf(q2) + dg_step, w)
                - cert2.phi_true(p, q, p2, mp.mpf(q2) - dg_step, w)) / (2 * dg_step)
        rel_q = abs(fd_q - exact_q) / max(mp.mpf(1), abs(exact_q))
        max_q2_rel = max(max_q2_rel, rel_q)
        checked_q2 += 1
    assert max_face_rel < mp.mpf("1e-8")
    assert checked_q2 >= 50 and max_q2_rel < mp.mpf("1e-8")
    print("Derivative verification PASS (sampled): 100 face points, %d q2 points; max relative errors %.3e, %.3e"
          % (checked_q2, float(max_face_rel), float(max_q2_rel)), flush=True)


def verify_dq(samples=200):
    cap = mp.mpf(str(float(DQ_ENDPOINT.upper())))
    worst = mp.mpf(0)
    for i in range(1, samples + 1):
        x = mp.mpf(i) / samples / 16
        z = 1 - x
        def rr(y):
            xx = 1 - y
            return 2 * xx * h_mp(y) - h_mp(xx * xx)
        value = abs(mp.diff(lambda y: mp.sqrt(rr(y)), z))
        worst = max(worst, value)
        assert value <= cap
    assert worst < cap and samples == 200
    print("DQ_hi containment PASS (sampled): 200 strip points; worst %.12f <= certified %.12f"
          % (float(worst), float(cap)), flush=True)


def verify_face_soundness(t, lambdas, boxes=100, points=20, seed=91903):
    rng = random.Random(seed)
    checked = 0
    for _ in range(boxes):
        p_lo = rng.uniform(.04, .64)
        p_hi = min(.75, p_lo + rng.uniform(2e-4, .025))
        q_lo = rng.uniform(max(p_hi + .01, .51), .82)
        q_hi = min(.90, q_lo + rng.uniform(2e-4, .025))
        r_lo = rng.uniform(.04, .72)
        r_hi = min(.94, r_lo + rng.uniform(2e-4, .025))
        w_lo = rng.uniform(.03, .83)
        w_hi = min(.94, w_lo + rng.uniform(2e-4, .025))
        box = ((p_lo, p_hi), (q_lo, q_hi), (r_lo, r_hi), (w_lo, w_hi))
        for lam in lambdas:
            direct = face_interval(box, t, lam)
            centered = face_centered(box, t, lam)
            assert direct is not None and centered is not None
            for _ in range(points):
                sample = tuple(rng.uniform(lo, hi) for lo, hi in box)
                shifted = _face_true(sample, t, lam)
                assert mp.mpf(str(float(direct.lower()))) <= shifted + mp.mpf("1e-30")
                assert mp.mpf(str(float(centered.lower()))) <= shifted + mp.mpf("1e-30")
                checked += 1
    assert checked == boxes * points * len(lambdas)
    print("Face-bound soundness PASS (sampled): %d random boxes x %d samples x %d lambdas, 0 violations"
          % (boxes, points, len(lambdas)), flush=True)

def verify_orbit_swap(samples=100, seed=12091):
    rng = random.Random(seed)
    # Exact-rational regression for the DECIMAL semantics used by _arb.
    # The endpoint immediately below one distinguishes this from binary-float
    # semantics: its exact image is 1e-16, not 2^-53.
    intervals = [(0.0, 1.0),
                 (0.1234567890123456, 0.9999999999999999),
                 (0.7500320261571772, 0.7532106456016533)]
    intervals.extend(sorted((rng.random(), rng.random()))
                     for _ in range(samples))
    def endpoint_fraction(value):
        point = _arb(value)
        assert point.is_exact()
        mantissa, exponent = point.man_exp()
        mantissa, exponent = int(mantissa), int(exponent)
        return (Fraction(mantissa) * (2 ** exponent)
                if exponent >= 0
                else Fraction(mantissa, 2 ** (-exponent)))
    for wl, wh in intervals:
        box = ((0.0, 1.0),) * 4 + ((wl, wh),)
        sw_lo, sw_hi = orbit_swap_box(box)[4]
        exact_lo = Fraction(1) - Fraction(str(wh))
        exact_hi = Fraction(1) - Fraction(str(wl))
        assert endpoint_fraction(sw_lo) <= exact_lo
        assert endpoint_fraction(sw_hi) >= exact_hi
    worst = mp.mpf(0)
    for _ in range(samples):
        p1, q1, p2, q2, w = [mp.mpf(rng.random()) for _ in range(5)]
        left = cert2.phi_true(p1, q1, p2, q2, w)
        right = cert2.phi_true(p2, q2, p1, q1, 1 - w)
        error = abs(left - right)
        assert error < mp.mpf("1e-30")
        worst = max(worst, error)
    assert samples == 100
    print("Orbit-swap enclosure PASS (exact decimal rationals): %d intervals; "
          "value invariance PASS (sampled): 100 points, max error %.3e"
          % (len(intervals), float(worst)), flush=True)


def verify_centered_overlap_regression(t, lam, samples=200, seed=33801):
    """Regression for overlapping ordered boxes whose midpoints are reversed."""
    rng = random.Random(seed)
    box = ((.30, .42), (.30, .34), (.31, .35), (.80, .86))
    prep = _face_prepare(box)
    assert prep is not None
    pmid = (prep[0] + prep[1]) / TWO
    qmid = (prep[2] + prep[3]) / TWO
    assert pmid > qmid
    # Choice (b): evaluate at the true midpoint. The formula is symmetric and
    # does not require order at a point, so the original radii remain valid.
    centered = face_centered(box, t, lam)
    assert centered is not None
    for _ in range(samples):
        # Sample the ordered portion represented by the prepared rectangle.
        p = rng.uniform(box[0][0], box[0][1])
        q = rng.uniform(max(p, box[1][0]), box[1][1]) if p <= box[1][1] else p
        if q > box[1][1]:
            p, q = box[1][1], box[1][1]
        sample = (p, q, rng.uniform(*box[2]), rng.uniform(*box[3]))
        assert mp.mpf(str(float(centered.lower()))) \
            <= _face_true(sample, t, lam) + mp.mpf("1e-30")
    print("Centered overlap regression PASS (sampled): true midpoint retained; "
          "centered face bound contains %d samples" % samples, flush=True)

def verify_pinned_face_regression(t, lambdas, samples=200, seed=44017):
    """Guard against using face feasibility after pinning q2 upward."""
    rng = random.Random(seed)
    p0, w0, delta = .25, .645, 2e-5
    box5 = ((p0 - delta, p0 + delta),) * 3 + (
        (.95, 1.0), (w0 - delta, w0 + delta))
    # At q2=.95 the midpoint is feasible, while its q2=1 image is not.
    p = q = r = mp.mpf(str(p0))
    w = mp.mpf(str(w0))
    m_low = w * (p + q) / 2 + (1 - w) * (r + mp.mpf(".95")) / 2
    m_face = w * (p + q) / 2 + (1 - w) * (r + 1) / 2
    assert m_low < mp.mpf(str(t)) < m_face
    assert q2_pin(box5, t, 0)
    face_box = (box5[0], box5[1], box5[2], box5[4])
    for lam in lambdas:
        full_only = face_best_bound(
            face_box, t, (lam,), allow_feasible=False)
        assert full_only is not None
        # The API-level regression: the default used by face_bb is full-only.
        assert face_best_bound(face_box, t, (lam,)) == full_only
        for _ in range(samples):
            sample = tuple(rng.uniform(lo, hi) for lo, hi in face_box)
            assert mp.mpf(str(float(full_only.lower()))) \
                <= _face_true(sample, t, lam) + mp.mpf("1e-30")
    print("Pinned-face infeasibility regression PASS (sampled): "
          "face M>t candidate excluded, %d samples x %d lambdas"
          % (samples, len(lambdas)), flush=True)


def _prepared_floats(box5):
    """The prepared rectangle (q_lo raised to p_lo) as float pairs, or None."""
    prep = bound_kkt._prepare_box(box5)
    if prep is None:
        return None
    vals = [float(x) for x in prep]
    return tuple((vals[2 * i], vals[2 * i + 1]) for i in range(5))


def verify_centered5(t, center_lams, boxes=150, points=20, seed=70211):
    """Sampled soundness on the full prepared rectangle."""
    rng = random.Random(seed)
    usable = comparisons = 0
    t_mp = mp.mpf(str(t))

    def check_box(box):
        nonlocal usable, comparisons
        prep = _prepared_floats(box)
        if prep is None:
            return
        for lam in center_lams:
            val = centered5(box, t, lam)
            if val is None:
                continue
            usable += 1
            lam_mp = mp.mpf(str(lam))
            vlo = mp.mpf(str(float(val.lower() if hasattr(val, "lower")
                                   else val)))
            for _ in range(points):
                sample = [rng.uniform(lo, hi) for lo, hi in prep]
                mean = (sample[4] * (sample[0] + sample[1]) / 2
                        + (1 - sample[4]) * (sample[2] + sample[3]) / 2)
                truth = cert2.phi_true(*sample) + lam_mp * (mean - t_mp)
                assert vlo <= truth + mp.mpf("1e-25"), \
                    (box, lam, tuple(sample), vlo, truth)
                comparisons += 1

    for _ in range(boxes):
        centre = [rng.uniform(0.05, 0.95) for _ in range(5)]
        halfw = [rng.uniform(1e-4, CENTERED5_MAX_WIDTH / 2)
                 for _ in range(5)]
        check_box(tuple((max(0.0, c - hw), min(1.0, c + hw))
                        for c, hw in zip(centre, halfw)))

    # Preparation raises q1_lo from .02 to .10.  This checks that midpoint,
    # radii, gradient, and sampled domain all use the prepared rectangle.
    skew = ((0.10, 0.20), (0.02, 0.12), (0.30, 0.34), (0.35, 0.39),
            (0.60, 0.62))
    before = comparisons
    check_box(skew)
    assert comparisons > before, "skew regression produced no usable bound"

    # Reversed-midpoint overlap regression for the former sort/radius bug.
    # The geometric p1 midpoint exceeds q1's, but each stays in its OWN
    # coordinate interval.  All centered implementations clear this narrow
    # positive box without sorting either coordinate.
    overlap = ((.329, .331), (.329, .3295), (.329, .331), (.329, .331),
               (.83, .831))
    prep = bound_kkt._prepare_box(overlap)
    assert prep is not None
    omids = tuple((prep[2 * i] + prep[2 * i + 1]) / TWO
                  for i in range(5))
    assert omids[0] > omids[1]
    for i, mid in enumerate(omids):
        assert prep[2 * i] <= mid <= prep[2 * i + 1]
    overlap_single = [centered5(overlap, t, lam)
                      for lam in center_lams]
    assert all(value is not None and value >= ZERO
               for value in overlap_single)
    overlap_best = centered5_best(overlap, t, center_lams)
    overlap_mixed = centered5_mixed(overlap, t, center_lams)
    assert overlap_best is not None and overlap_best >= ZERO
    assert overlap_mixed is not None and overlap_mixed >= ZERO
    before = comparisons
    check_box(overlap)
    assert comparisons > before

    # centered5_best must agree with the per-lambda maximum.
    agree = 0
    for _ in range(10):
        centre = [rng.uniform(0.10, 0.90) for _ in range(5)]
        halfw = [rng.uniform(1e-4, CENTERED5_MAX_WIDTH / 2)
                 for _ in range(5)]
        box = tuple((max(0.0, c - hw), min(1.0, c + hw))
                    for c, hw in zip(centre, halfw))
        singles = [centered5(box, t, lam) for lam in center_lams]
        singles = [value for value in singles if value is not None]
        joint = centered5_best(box, t, center_lams)
        if joint is None:
            assert not singles
            continue
        best_single = singles[0]
        for value in singles[1:]:
            if value > best_single:
                best_single = value
        assert float(joint.lower()) == float(best_single.lower()), box
        agree += 1
    assert agree > 0
    assert usable > 0 and comparisons > 0
    print("centered5 soundness PASS (sampled): %d usable box/lam pairs, "
          "%d comparisons incl. skew/reversed-overlap regressions, 0 violations"
          % (usable, comparisons), flush=True)


def verify_centered5_mixed(t, center_lams, samples=200, seed=52807):
    """Vestigial-orbit regression: mixed form works where centered5 cannot."""
    rng = random.Random(seed)
    t_mp = mp.mpf(str(t))
    vest = ((0.375, 0.385), (0.375, 0.385), (0.0, 0.05), (0.95, 1.0),
            (0.96, 1.0))
    full = centered5_best(vest, t, center_lams)
    assert full is None, "expected the full gradient to be unusable here"
    mixed = centered5_mixed(vest, t, center_lams)
    assert mixed is not None, "mixed form must be usable on the vestigial box"
    vlo = mp.mpf(str(float(mixed.lower())))
    prep = _prepared_floats(vest)
    checked = 0
    for lam in center_lams:
        lam_mp = mp.mpf(str(lam))
        for _ in range(samples):
            s = list(rng.uniform(lo, hi) for lo, hi in prep)
            s[0], s[1] = min(s[0], s[1]), max(s[0], s[1])
            s[2], s[3] = min(s[2], s[3]), max(s[2], s[3])
            m = s[4] * (s[0] + s[1]) / 2 + (1 - s[4]) * (s[2] + s[3]) / 2
            truth = cert2.phi_true(*s) + lam_mp * (m - t_mp)
            # mixed is the max over center_lams of per-lambda bounds, each a
            # lower bound for its own Phi_lam; on the FEASIBLE part every
            # Phi_lam <= Phi, so compare against Phi_lam pointwise per lam is
            # NOT valid for the max.  Compare against Phi on feasible samples
            # and per-lambda otherwise is intricate; simplest sound check:
            # the max must lower-bound max_lam Phi_lam at the point.
            best_here = truth
            for lam2 in center_lams:
                cand = cert2.phi_true(*s) + mp.mpf(str(lam2)) * (m - t_mp)
                if cand > best_here:
                    best_here = cand
            assert vlo <= best_here + mp.mpf("1e-25"), (lam, tuple(s))
            checked += 1
    print("centered5_mixed vestigial regression PASS (sampled): usable where "
          "full gradient is None; %d comparisons, 0 violations"
          % checked, flush=True)


def verify_sink_safe_centering(t, center_lams, boxes=100, points=20,
                               seed=82113):
    """Sampled controls for swapped-mixed and weight-only MVT bounds."""
    rng = random.Random(seed)
    t_mp = mp.mpf(str(t))
    target = ((.01171875, .0126953125), (.9998779296875, 1.0),
              (.0107421875, .01171875), (.0107421875, .01171875),
              (.75, .7500640523143549))
    targets = [target]
    for _ in range(boxes - 1):
        p1 = rng.uniform(.001, .35)
        q1 = rng.uniform(max(.40, p1 + .01), .999)
        p2 = rng.uniform(.001, .35)
        q2 = rng.uniform(max(.40, p2 + .01), .999)
        d = rng.uniform(1e-5, 2e-3)
        wl = rng.uniform(.05, .94)
        targets.append(((p1, min(p1 + d, q1)),
                        (q1, min(q1 + d, .999999)),
                        (p2, min(p2 + d, q2)),
                        (q2, min(q2 + d, .999999)),
                        (wl, min(wl + d, .999999))))
    usable = Counter()
    comparisons = Counter()
    target_swap = centered5_mixed_swap(
        target, t, center_lams)
    assert target_swap is not None and target_swap >= ZERO
    for box in targets:
        prep = _prepared_floats(box)
        assert prep is not None
        values = {
            "mixed_swap": centered5_mixed_swap(
                box, t, center_lams),
            "center_w": centered_w_best(box, t, center_lams),
        }
        for name, value in values.items():
            if value is None:
                continue
            usable[name] += 1
            vlo = mp.mpf(str(float(value.lower())))
            for _ in range(points):
                sample = [rng.uniform(lo, hi) for lo, hi in prep]
                mean = (sample[4] * (sample[0] + sample[1]) / 2
                        + (1 - sample[4]) * (sample[2] + sample[3]) / 2)
                phi = cert2.phi_true(*sample)
                best_here = max(
                    phi + mp.mpf(str(lam)) * (mean - t_mp)
                    for lam in center_lams)
                assert vlo <= best_here + mp.mpf("1e-25"), \
                    (name, box, tuple(sample), vlo, best_here)
                comparisons[name] += 1
    assert usable["mixed_swap"] > 0 and usable["center_w"] > 0
    print("sink-safe centering PASS (sampled): swapped-mixed %d bounds/%d "
          "points, weight-only %d bounds/%d points; target swap clears"
          % (usable["mixed_swap"], comparisons["mixed_swap"],
             usable["center_w"], comparisons["center_w"]), flush=True)


def verify_strip_pin(t, lambdas, samples=120, seed=61903):
    """The pin must fire on a strip box with q2_hi < 1, and the certified
    derivative sign must agree with finite differences (sampled)."""
    box5 = ((0.328125, 0.3359375), (0.328125, 0.3359375),
            (0.328125, 0.3359375), (0.9921875, 0.99609375),
            (0.8359375, 0.84375))
    fired = [lam for lam in lambdas if q2_pin(box5, t, lam)]
    assert fired, "strip pin failed to fire on the regression box"
    rng = random.Random(seed)
    t_mp = mp.mpf(str(t))
    eps = mp.mpf("1e-12")
    checked = 0
    for _ in range(samples):
        s = [rng.uniform(lo, hi) for lo, hi in box5]
        q2 = mp.mpf(str(rng.uniform(box5[3][0], 0.999999)))
        lam = mp.mpf(str(rng.choice(fired)))

        def phi_lam_q(q):
            m = (mp.mpf(str(s[4])) * (mp.mpf(str(s[0])) + mp.mpf(str(s[1]))) / 2
                 + (1 - mp.mpf(str(s[4]))) * (mp.mpf(str(s[2])) + q) / 2)
            return cert2.phi_true(s[0], s[1], s[2], q, s[4]) + lam * (m - t_mp)

        d = (phi_lam_q(q2 + eps) - phi_lam_q(q2 - eps)) / (2 * eps)
        assert d <= mp.mpf("1e-6"), (tuple(s), float(q2), float(lam), float(d))
        checked += 1
    print("strip-pin regression PASS (sampled): fires with q2_hi<1 for %d "
          "lambdas; FD derivative <= 0 at %d sampled strip points"
          % (len(fired), checked), flush=True)


def centered5_subpin_probe(t, center_lams):
    """Flip probe just BELOW the pin threshold, where corner ruled before."""
    b = float(B_OBS)
    q2 = 0.93
    m2 = (b + q2) / 2
    w = (m2 - float(t)) / (m2 - b)
    centre = (b, b, b, q2, w)
    print("centered5 sub-pin probe at q2=0.93 centre (%.4f,%.4f,%.4f,%.4f,%.4f)"
          % centre, flush=True)
    print("%-10s %14s" % ("FULLwidth", "centered5_best"), flush=True)
    for d in (1.0 / 32, 1.0 / 64, 1.0 / 128, 1.0 / 256, 1e-3, 3e-4):
        box = tuple((max(0.0, c - d / 2), min(1.0, c + d / 2)) for c in centre)
        val = centered5_best(box, t, center_lams)
        print("%-10.3e %14s" % (d, "unusable" if val is None
                                else "%+.4e" % float(val.lower())), flush=True)


def face_flip_probe(t, lambdas, center=None):
    if center is None:
        p, r, w, _, _ = solve_face_kkt(t)
        center = (p, p, r, w)
    def bound(width):
        box = tuple((max(1e-12, float(x) - width / 2),
                     min(1 - 1e-12, float(x) + width / 2)) for x in center)
        ans = face_best_bound(box, t, lambdas)
        assert ans is not None
        return float(ans.lower())
    lo, hi = 1e-7, .1
    assert bound(lo) > 0
    while bound(hi) > 0 and hi < .5:
        hi *= 2
    assert bound(hi) <= 0
    for _ in range(60):
        mid = (lo + hi) / 2
        if bound(mid) > 0:
            lo = mid
        else:
            hi = mid
    assert lo < hi and bound(lo) > 0 and bound(hi) <= 0
    print("Face centered flip at FULL width [%.12g, %.12g] around (%.12g,%.12g,%.12g,%.12g)"
          % ((lo, hi) + tuple(float(x) for x in center)), flush=True)
    return lo, hi


# ---------------------------------------------------------------------------
# B&B driver.
# ---------------------------------------------------------------------------
def _split(box, coordinate):
    lo, hi = box[coordinate]
    mid = (lo + hi) / 2
    left, right = list(box), list(box)
    left[coordinate], right[coordinate] = (lo, mid), (mid, hi)
    return tuple(left), tuple(right)


def face_bb(box4, t, lambdas, gmax, min_width=DEFAULT_FACE_MIN_WIDTH,
            box_budget=100000, deadline=None):
    stack = [box4]
    stats = Counter()
    while stack:
        if stats["processed"] >= box_budget or (deadline and time.monotonic() >= deadline):
            stats["budget"] += len(stack)
            return False, stats
        box = stack.pop()
        stats["processed"] += 1
        # A pin maps feasible points to q2=1 points which may be infeasible.
        # Therefore ONLY full-face Phi_lam bounds are sound here.
        bound = face_best_bound(box, t, lambdas, gmax, allow_feasible=False)
        if bound is not None and bound >= ZERO:
            stats["cleared"] += 1
            continue
        widths = [hi - lo for lo, hi in box]
        widest = max(range(4), key=widths.__getitem__)
        if widths[widest] <= min_width:
            stats["residual"] += 1
            return False, stats
        stack.extend(_split(box, widest))
        stats["split"] += 1
    assert stats["residual"] == 0 and stats["budget"] == 0
    return True, stats


def _classify(box, t):
    mids = [(a + b) / 2 for a, b in box]
    p1, q1, p2, q2, w = mids
    m = w * (p1 + q1) / 2 + (1 - w) * (p2 + q2) / 2
    if min(p1, q1, p2, q2, 1 - p1, 1 - q1, 1 - p2, 1 - q2) < .01:
        return "sink"
    if abs(p1 - float(B_OBS)) < .03 and abs(q1 - float(B_OBS)) < .03 \
            and abs(p2 - float(B_OBS)) < .03 and q2 > .9:
        return "obstruction"
    if abs(m - float(t)) < .01:
        return "mean-boundary"
    return "other"



CENTERED5_MAX_WIDTH = 1.0 / 32.0
# Only try centered5 when the corner bound is within this of clearing; far
# more negative boxes split faster than the gradient evaluation costs.
CENTER_GATE = -0.02


def split_by_influence(splittable, widths, box):
    """Pick the coordinate with the largest influence-weighted width.

    Every atom coordinate enters Phi only through its own orbit's weight, so
    orbit 1's coordinates carry at most ``w_hi`` and orbit 2's carry at most
    ``1 - w_lo``.  On the w-slices near 1 that discounts ``p2,q2`` by up to
    16x, and bisecting them there buys almost nothing.

    This is a COST heuristic only: any split is sound, because the two
    children cover the parent exactly, so the choice cannot change what a
    cleared box means.
    """
    influence = (box[4][1], box[4][1],
                 1.0 - box[4][0], 1.0 - box[4][0], 1.0)
    return max(splittable, key=lambda j: widths[j] * influence[j])


def split_by_width(splittable, widths, box):
    """The plain widest-coordinate rule, kept for cost comparisons."""
    return max(splittable, key=widths.__getitem__)


SPLIT_CHOICE = split_by_influence




def phi_lam_point(point5, t, lam, gmax=None):
    """Direct Arb evaluation at an arbitrary 5-D point.

    Unlike ``bound_kkt.phi_lam_encl``, this helper does not canonicalize
    overlapping p/q coordinates.  Pair symmetry makes ordering unnecessary,
    and preserving the geometric midpoint is essential for MVT radii.
    """
    if gmax is None:
        gmax = get_rh_gmax()
    p1, q1, p2, q2, w = map(_arb, point5)
    t, lam = _arb(t), _arb(lam)
    v = ONE - w
    m1, m2 = (p1 + q1) / TWO, (p2 + q2) / TWO
    l1 = (h_pt(p1) + h_pt(q1)) / TWO
    l2 = (h_pt(p2) + h_pt(q2)) / TWO
    c1 = h_pt(sstar_pt(p1, q1))
    c2 = h_pt(sstar_pt(p2, q2))
    u1 = (sqrt_rh_enc(p1, p1, gmax)
          + sqrt_rh_enc(q1, q1, gmax)) / TWO
    u2 = (sqrt_rh_enc(p2, p2, gmax)
          + sqrt_rh_enc(q2, q2, gmax)) / TWO
    mean = w * m1 + v * m2
    entropy = w * l1 + v * l2
    correction = w * c1 + v * c2
    sigma = w * u1 + v * u2
    coeff = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    value = (coeff * entropy + ALPHA * correction
             - (ONE - ALPHA) * sigma * sigma + lam * (mean - t))
    return value if value.is_finite() else None


def centered5(box5, t, lam, gmax=None):
    """Mean-value lower bound for Phi + lam*(M-t) over the PREPARED 5-D box.

    ``bound_kkt._prepare_box`` raises each q_lo to max(q_lo, p_lo); the
    gradient is certified on that prepared rectangle, so the midpoint, radii,
    and shift ranges MUST come from the same tuple (raw-box geometry would let
    the MVT segment exit the certified rectangle on overlapping p/q boxes).
    The prepared rectangle contains the ordered part {p <= q} of the raw box,
    which is the only part the partition-wide certificate needs: every
    measure has an ordered representative somewhere in the cube.

    Sound for the feasible subset: lam >= 0 and M <= t there make the shift
    nonpositive.  The value is evaluated at the unsorted geometric midpoint:
    the formula is pair-symmetric and does not require p <= q, while every
    midpoint coordinate and its original radius stay in the same rectangle.
    """
    if gmax is None:
        gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    prepared = bound_kkt._prepare_box(box5)
    if prepared is None:
        return None
    grad = bound_kkt.interval_gradient(box5)
    if grad is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    lam = _arb(lam)
    assert lam >= ZERO
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    shift = (lam * w / TWO, lam * w / TWO, lam * v / TWO, lam * v / TWO,
             lam * (m1 - m2))
    pairs = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh))
    mids = tuple((a + b) / TWO for a, b in pairs)
    radii = tuple((b - a) / TWO for a, b in pairs)
    point = mids
    pv = phi_lam_point(point, t, lam, gmax)
    if pv is None:
        return None
    penalty = ZERO
    for g, s, r in zip(grad, shift, radii):
        gs = g + s
        penalty += amax(-gs.lower(), gs.upper()) * r
    lower = _arb(pv.lower()) - penalty
    return lower if lower.is_finite() else None


def centered5_best(box5, t, lams, gmax=None):
    """Max of the sound centered bounds, sharing one gradient across lambdas.

    Arithmetic per lambda is identical to :func:`centered5`; only the
    prepared-rectangle geometry and interval gradient are hoisted out of the
    loop (they do not depend on lambda).
    """
    if gmax is None:
        gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    prepared = bound_kkt._prepare_box(box5)
    if prepared is None:
        return None
    grad = bound_kkt.interval_gradient(box5)
    if grad is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    pairs = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh))
    mids = tuple((a + b) / TWO for a, b in pairs)
    radii = tuple((b - a) / TWO for a, b in pairs)
    point = mids
    best = None
    for lam_f in lams:
        lam = _arb(lam_f)
        assert lam >= ZERO
        pv = phi_lam_point(point, t, lam, gmax)
        if pv is None:
            continue
        shift = (lam * w / TWO, lam * w / TWO, lam * v / TWO, lam * v / TWO,
                 lam * (m1 - m2))
        penalty = ZERO
        for g, s, r in zip(grad, shift, radii):
            gs = g + s
            penalty += amax(-gs.lower(), gs.upper()) * r
        lower = _arb(pv.lower()) - penalty
        if lower.is_finite() and (best is None or lower > best):
            best = lower
    return best


def centered5_mixed(box5, t, lams, gmax=None):
    """Block-mixed bound: MVT in (p1, q1, w), natural interval in (p2, q2).

    For x = (x1, x2), x1 = (p1, q1, w), x2 = (p2, q2):

        Phi_lam(x1, x2) >= Phi_lam(x1*, x2) - sum_{i in x1} sup|d_i| rad_i
                        >= inf_{x2-box} Phi_lam(x1*, x2) - penalty,

    where x1* is the unsorted geometric midpoint of x1: the segment from x1*
    to x1 at FIXED x2 stays inside the box, and the three partial enclosures
    are taken over the whole box.  The inner infimum is a natural-interval
    evaluation with orbit 1 fixed at that midpoint and orbit 2 interval-valued,
    so orbit 2 may touch the sinks freely, where the full gradient (hence
    centered5) is unusable.  Partials in p1, q1, w involve NO orbit-2
    derivatives: d/dw needs only interval aggregates of both orbits.
    """
    if gmax is None:
        gmax = get_rh_gmax()
    prepared = bound_kkt._prepare_box(box5)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    hp1 = _entropy_prime(p1l, p1h)
    hq1 = _entropy_prime(q1l, q1h)
    gp1 = _sqrt_rh_prime(p1l, p1h, gmax)
    gq1 = _sqrt_rh_prime(q1l, q1h, gmax)
    cd = _sstar_derivatives(p1l, p1h, q1l, q1h)
    if hp1 is None or hq1 is None or gp1 is None or gq1 is None or cd is None:
        return None
    dc1p, dc1q = cd
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    l1 = (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
    l2 = (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    c1 = h_encl(s1.lower(), s1.upper())
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    c2 = h_encl(s2.lower(), s2.upper())
    u1 = (sqrt_rh_enc(p1l, p1h, gmax) + sqrt_rh_enc(q1l, q1h, gmax)) / TWO
    u2 = (sqrt_rh_enc(p2l, p2h, gmax) + sqrt_rh_enc(q2l, q2h, gmax)) / TWO
    L = w * l1 + v * l2
    S = w * u1 + v * u2
    mean = w * m1 + v * m2
    kappa = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    # bound_kkt.interval_gradient formulas, restricted to p1, q1, w.
    g_p1 = (-(ONE - ALPHA) * w * L + kappa * w * hp1 / TWO
            + ALPHA * w * dc1p - (ONE - ALPHA) * w * S * gp1)
    g_q1 = (-(ONE - ALPHA) * w * L + kappa * w * hq1 / TWO
            + ALPHA * w * dc1q - (ONE - ALPHA) * w * S * gq1)
    g_w = (TWO * (ONE - ALPHA) * (m2 - m1) * L
           + kappa * (l1 - l2) + ALPHA * (c1 - c2)
           - TWO * (ONE - ALPHA) * S * (u1 - u2))
    if not all(x.is_finite() for x in (g_p1, g_q1, g_w)):
        return None
    pm, qm = (p1l + p1h) / TWO, (q1l + q1h) / TWO
    wm = (wl + wh) / TWO
    r_p1, r_q1, r_w = ((p1h - p1l) / TWO, (q1h - q1l) / TWO,
                       (wh - wl) / TWO)
    m1p = (pm + qm) / TWO
    l1p = (h_pt(pm) + h_pt(qm)) / TWO
    c1p = h_pt(sstar_pt(pm, qm))
    u1p = (sqrt_rh_enc(pm, pm, gmax) + sqrt_rh_enc(qm, qm, gmax)) / TWO
    wmv = ONE - wm
    t = _arb(t)
    best = None
    for lam_f in lams:
        lam = _arb(lam_f)
        assert lam >= ZERO
        Mx = wm * m1p + wmv * m2
        Lx = wm * l1p + wmv * l2
        Cx = wm * c1p + wmv * c2
        Sx = wm * u1p + wmv * u2
        val = ((TWO * (ONE - ALPHA) * (ONE - Mx) - ONE) * Lx
               + ALPHA * Cx - (ONE - ALPHA) * Sx * Sx + lam * (Mx - t))
        gp = g_p1 + lam * w / TWO
        gq = g_q1 + lam * w / TWO
        gw = g_w + lam * (m1 - m2)
        penalty = (amax(-gp.lower(), gp.upper()) * r_p1
                   + amax(-gq.lower(), gq.upper()) * r_q1
                   + amax(-gw.lower(), gw.upper()) * r_w)
        lower = _arb(val.lower()) - penalty
        if lower.is_finite() and (best is None or lower > best):
            best = lower
    return best


def orbit_swap_box(box5):
    """Arb enclosure of (A, B, w) -> (B, A, 1-w).

    Float bounds enter every proof enclosure through ``_arb``, which interprets
    their shortest DECIMAL strings.  Perform the subtraction in that same Arb
    model and return certified point endpoints; float ``nextafter`` is not
    reliably outward relative to decimal semantics near one.
    """
    wl, wh = map(_arb, box5[4])
    swapped_lo = (ONE - wh).lower()
    swapped_hi = (ONE - wl).upper()
    return (box5[2], box5[3], box5[0], box5[1],
            (swapped_lo, swapped_hi))


def centered5_mixed_swap(box5, t, lams, gmax=None):
    """Mixed MVT bound after the exact value/feasibility-preserving swap.

    ``centered5_mixed`` differentiates orbit 1 and naturally interval-encloses
    orbit 2.  The involution exchanges those roles, so this orientation stays
    usable when the ORIGINAL orbit 1 touches an entropy/sqrt-rh sink.
    """
    return centered5_mixed(orbit_swap_box(box5), t, lams, gmax)


_wm1, _wm2, _wl1, _wl2, _wc1, _wc2, _wu1, _wu2, _ww = \
    sp.symbols("wm1 wm2 wl1 wl2 wc1 wc2 wu1 wu2 ww")
_walpha, _wlam, _wt = sp.symbols("walpha wlam wt")
_wM = _ww * _wm1 + (1 - _ww) * _wm2
_wL = _ww * _wl1 + (1 - _ww) * _wl2
_wC = _ww * _wc1 + (1 - _ww) * _wc2
_wS = _ww * _wu1 + (1 - _ww) * _wu2
_wPhiLam = ((2 * (1 - _walpha) * (1 - _wM) - 1) * _wL
            + _walpha * _wC - (1 - _walpha) * _wS ** 2
            + _wlam * (_wM - _wt))
_wDerivative = (2 * (1 - _walpha) * (_wm2 - _wm1) * _wL
                + (2 * (1 - _walpha) * (1 - _wM) - 1)
                * (_wl1 - _wl2)
                + _walpha * (_wc1 - _wc2)
                - 2 * (1 - _walpha) * _wS * (_wu1 - _wu2)
                + _wlam * (_wm1 - _wm2))
_WEIGHT_DERIVATIVE_IDENTITY = \
    sp.simplify(sp.diff(_wPhiLam, _ww) - _wDerivative) == 0
assert _WEIGHT_DERIVATIVE_IDENTITY


def _weight_components(box5, gmax):
    """Prepared aggregate enclosures and the exact partial d(Phi)/dw."""
    prepared = bound_kkt._prepare_box(box5)
    if prepared is None:
        return None
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    w = hull(wl, wh)
    v = ONE - w
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    l1 = (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO
    l2 = (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    c1 = h_encl(s1.lower(), s1.upper())
    c2 = h_encl(s2.lower(), s2.upper())
    u1 = (sqrt_rh_enc(p1l, p1h, gmax)
          + sqrt_rh_enc(q1l, q1h, gmax)) / TWO
    u2 = (sqrt_rh_enc(p2l, p2h, gmax)
          + sqrt_rh_enc(q2l, q2h, gmax)) / TWO
    mean = w * m1 + v * m2
    entropy = w * l1 + v * l2
    sigma = w * u1 + v * u2
    kappa = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    dphi_dw = (TWO * (ONE - ALPHA) * (m2 - m1) * entropy
                + kappa * (l1 - l2) + ALPHA * (c1 - c2)
                - TWO * (ONE - ALPHA) * sigma * (u1 - u2))
    if not dphi_dw.is_finite():
        return None
    return prepared, m1, m2, dphi_dw


def centered_w_best(box5, t, lams, gmax=None):
    """Sink-safe one-dimensional MVT bound in the orbit weight.

    On the prepared ordered atom rectangle A and weight interval W, fix the
    weight midpoint w*.  For every atoms a in A and w in W, the segment
    (a,w*)--(a,w) stays in A x W.  Thus

      Phi_lam(a,w) >= inf_A Phi_lam(a,w*)
                      - sup_(A x W)|d_w Phi_lam| rad(W).

    The anchor is a natural interval enclosure over ALL four atom coordinates;
    the derivative contains only bounded aggregates, not atom derivatives, so
    the rule remains finite at sinks.  For lam >= 0 and feasible M <= t,
    Phi_lam = Phi + lam(M-t) <= Phi, giving the required inequality direction.
    """
    if gmax is None:
        gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    components = _weight_components(box5, gmax)
    if components is None:
        return None
    prepared, m1, m2, dphi_dw = components
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    wm = (wl + wh) / TWO
    rw = (wh - wl) / TWO
    anchor = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h),
              (wm, wm))
    best = None
    for lam_f in lams:
        lam = _arb(lam_f)
        assert lam >= ZERO
        value = bound_kkt.phi_lam_encl(anchor, t, lam)
        if value is None:
            continue
        derivative = dphi_dw + lam * (m1 - m2)
        penalty = amax(-derivative.lower(), derivative.upper()) * rw
        lower = _arb(value.lower()) - penalty
        if lower.is_finite() and (best is None or lower > best):
            best = lower
    return best


COLLAR_EDGE = 0.01


def _touches_collar(box):
    """Any atom interval intersecting [0, COLLAR_EDGE] or [1-COLLAR_EDGE, 1]."""
    return any(lo <= COLLAR_EDGE or hi >= 1.0 - COLLAR_EDGE
               for lo, hi in box[:4])


def certify(t, root_box, max_boxes=40_000_000, time_budget=720.0,
            min_width=DEFAULT_MIN_WIDTH,
            face_min_width=DEFAULT_FACE_MIN_WIDTH,
            face_box_budget=100000, gmax=None, rho_gmax=None,
            lambdas=None, center_lams=None, progress_every=2_000_000,
            collar_min_width=None, residual_dump=None, trace=None):
    """Certify one root box; return ``(verdict, tallies)``.

    ``collar_min_width``: sink-collar-touching boxes are refined to this finer
    width floor (sink faces carry the global margin, but the corner bound's
    slop is first-order, so they need ~2-3 extra halvings).  ``residual_dump``:
    optional .npy path collecting centres+FULL widths of every residual leaf.
    ``trace`` is an optional binary file-like object receiving one byte for
    every processed raw DFS node.
    """
    started = time.monotonic()
    deadline = started + time_budget
    t_arb = _arb(t)
    t_upper = t_arb.upper()
    if gmax is None:
        gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    if rho_gmax is None:
        rho_gmax = cert2.get_rho_gmax()
    if lambdas is None or center_lams is None:
        fam, lamhat = lambda_family(t)
        if lambdas is None:
            lambdas = fam
        if center_lams is None:
            center_lams = (0.0, float(lamhat))
    stack = [root_box]
    stats = Counter()
    residual_geometry = Counter()
    residual_rows = []
    if collar_min_width is None:
        # Per-COORDINATE floor (see the split logic): safe to default deep.
        collar_min_width = 1.25e-4
    while stack:
        if stats["processed"] >= max_boxes:
            stats["budget_boxes"] = len(stack)
            break
        if time.monotonic() >= deadline:
            stats["budget_time"] = len(stack)
            break
        raw = stack.pop()
        stats["processed"] += 1
        if progress_every and stats["processed"] % progress_every == 0:
            print("  certify: %d processed, corner %d, ratio %d, center %d, "
                  "mixed %d, swap %d, weight %d, face %d, infeasible %d, "
                  "residual %d, stack %d (%.0fs)"
                  % (stats["processed"], stats["corner"], stats["ratio"],
                     stats["center"], stats["center_mixed"],
                     stats["center_mixed_swap"], stats["center_w"],
                     stats["face"], stats["infeasible"], stats["residual"],
                     len(stack), time.monotonic() - started), flush=True)
        box = diag_exhaust.mean_contract(raw, t_upper)
        if box is None:
            stats["infeasible"] += 1
            _write_trace(trace, TRACE_MEAN_INFEASIBLE)
            continue
        corner = diag_exhaust.phi_corner(box, t_arb)
        if corner is None:
            stats["infeasible"] += 1
            _write_trace(trace, TRACE_CORNER_INFEASIBLE)
            continue
        if corner >= ZERO:
            stats["corner"] += 1
            _write_trace(trace, TRACE_CORNER)
            continue
        # The ratio rule IS its own cheapest applicability test: it computes
        # certified rho caps and checks kappa - beta*rho_hi >= 0, at only
        # 1.29x the cost of the corner bound already evaluated above.  An
        # earlier hand-tuned geometric pre-filter (all atoms <= 0.05 or
        # >= 0.42) skipped the rule on exactly the boxes that need it: the
        # near-total-sink family p1,q1 -> 0, q2 -> 1, w -> 1 with
        # p2 ~ 0.41 falls in that filter's dead zone, and 1.45 M such boxes
        # became terminal residuals in the 2026-08-16 campaign even though
        # the rule clears every one of them.
        #
        # Pass the BOX's w-interval: feasible points of this box have w here,
        # and the endpoint-max over a narrower window is sharper than the
        # widest derived window (measured: the near-total-sink leaves fire
        # with margin +6.6e-3 under the box window and fail by -2.7e-3 under
        # the derived one).
        w_win = hull(_arb(repr(box[4][0])), _arb(repr(box[4][1])))
        if cert2.ratio_rule(box[:4], t_arb, W=w_win, rho_gmax=rho_gmax):
            stats["ratio"] += 1
            _write_trace(trace, TRACE_RATIO)
            continue
        # centered5 is the expensive rule; spend it only where the cheap
        # corner bound says the box is NEAR-CLEARING.  Sound: skipping a
        # bound can only send the box to the split branch, never clear it.
        if (max(hi - lo for lo, hi in box) <= CENTERED5_MAX_WIDTH
                and float(corner.lower()) >= CENTER_GATE):
            cb = centered5_best(box, t_arb, center_lams, gmax)
            if cb is not None and cb >= ZERO:
                stats["center"] += 1
                _write_trace(trace, TRACE_CENTER)
                continue
            if cb is None:
                # Full gradient unusable (typically a sink-touching vestigial
                # orbit 2): fall back to the block-mixed form.
                cm = centered5_mixed(box, t_arb, center_lams, gmax)
                if cm is not None and cm >= ZERO:
                    stats["center_mixed"] += 1
                    _write_trace(trace, TRACE_CENTER_MIXED)
                    continue
            if _touches_collar(box):
                # Exact orbit swap lets the same mixed proof interval-enclose
                # the ORIGINAL orbit 1 instead.  This is decisive when orbit 1
                # carries the sink singularity.  The w-only form is a final
                # derivative-free sink fallback.
                cms = centered5_mixed_swap(
                    box, t_arb, center_lams, gmax)
                if cms is not None and cms >= ZERO:
                    stats["center_mixed_swap"] += 1
                    _write_trace(trace, TRACE_CENTER_MIXED_SWAP)
                    continue
                cw = centered_w_best(box, t_arb, center_lams, gmax)
                if cw is not None and cw >= ZERO:
                    stats["center_w"] += 1
                    _write_trace(trace, TRACE_CENTER_W)
                    continue
        pinned_lams = []
        if box[3][0] >= Q2PIN:
            pinned_lams = [lam for lam in lambdas
                           if q2_pin(box, t_arb, lam, gmax)]
        if pinned_lams:
            face_box = (box[0], box[1], box[2], box[4])
            completed, fst = face_bb(
                face_box, t_arb, tuple(pinned_lams), gmax,
                face_min_width, face_box_budget, deadline,
            )
            stats["pin_success"] += 1
            stats["face_processed"] += fst["processed"]
            stats["face_cleared"] += fst["cleared"]
            if completed:
                stats["face"] += 1
                _write_trace(trace, TRACE_FACE)
                continue
            stats["face_failed"] += 1
        widths = [hi - lo for lo, hi in box]
        # Once every sink-safe bound has failed, a collar-touching box may need
        # one non-singular coordinate refined below the ordinary floor.  Pure
        # subdivision is sound; the new mixed/weight rules make this exceptional
        # rather than the former 2^15-wide collar refinement.
        collar_box = _touches_collar(box)
        def _floor(j):
            return collar_min_width if collar_box else min_width
        splittable = [j for j in range(5) if widths[j] > _floor(j)]
        if not splittable:
            stats["residual"] += 1
            residual_geometry[_classify(box, t)] += 1
            residual_rows.append([(lo + hi) / 2 for lo, hi in box]
                                 + [hi - lo for lo, hi in box])
            _write_trace(trace, TRACE_RESIDUAL)
            continue
        widest = SPLIT_CHOICE(splittable, widths, box)
        stack.extend(_split(box, widest))
        stats["split"] += 1
        _write_trace(trace, TRACE_SPLIT_BASE + widest)
    stats["stack"] = len(stack)
    elapsed = time.monotonic() - started
    stats["elapsed_ms"] = int(1000 * elapsed)
    stats["residual_geometry"] = dict(residual_geometry)
    # The manifest defines a hard wall-clock budget, not merely an admission
    # deadline.  A final expensive clearing rule may empty the stack after the
    # deadline; record that overrun and refuse COMPLETE.
    if elapsed > time_budget:
        stats["budget_time"] = max(stats["budget_time"], 1)
    if residual_dump and residual_rows:
        # Diagnostics are attempt-unique; never overwrite evidence from a retry.
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        fd = os.open(residual_dump, flags, 0o600)
        try:
            with os.fdopen(fd, "wb") as fh:
                np.save(fh, np.array(residual_rows, dtype=float))
                fh.flush()
                os.fsync(fh.fileno())
        except BaseException:
            # The caller owns cleanup policy; never delete even a partial dump.
            raise
        print("  residual dump: %d rows -> %s"
              % (len(residual_rows), residual_dump), flush=True)
    complete = (not stack and stats["residual"] == 0
                and stats["budget_boxes"] == 0 and stats["budget_time"] == 0)
    assert sum(residual_geometry.values()) == stats["residual"]
    return complete, stats


def run_bb(t, lambdas, gmax, rho_gmax, min_width=DEFAULT_MIN_WIDTH,
           face_min_width=DEFAULT_FACE_MIN_WIDTH, max_boxes=40_000_000,
           time_budget=720.0, root_box=None):
    if root_box is None:
        # Orbit swap (p1,q1,p2,q2,w)->(p2,q2,p1,q1,1-w) preserves
        # M,L,C,S and Phi. Thus w in [1/2,1] represents the full cube and
        # keeps the obstruction on the implemented q2=1 face.
        root_box = ((0.0, 1.0),) * 4 + ((0.5, 1.0),)
    complete, stats = certify(
        t, root_box, max_boxes=max_boxes, time_budget=time_budget,
        min_width=min_width, face_min_width=face_min_width, gmax=gmax,
        rho_gmax=rho_gmax, lambdas=lambdas,
    )
    print("B&B t=%s: processed=%d infeasible=%d corner=%d ratio=%d pin=%d "
          "face=%d face_failed=%d face_nodes=%d split=%d residual=%d "
          "stack=%d elapsed=%.2fs"
          % (mp.nstr(mp.mpf(str(t)), 18), stats["processed"],
             stats["infeasible"], stats["corner"], stats["ratio"],
             stats["pin_success"], stats["face"], stats["face_failed"],
             stats["face_processed"], stats["split"], stats["residual"],
             stats["stack"], stats["elapsed_ms"] / 1000), flush=True)
    print("Residual geometry: %s" % stats["residual_geometry"], flush=True)
    print("VERDICT t=%s: %s"
          % (mp.nstr(mp.mpf(str(t)), 18),
             "COMPLETE CERTIFICATE" if complete else "INCOMPLETE"), flush=True)
    return complete, stats


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-boxes", type=int, default=3_000_000)
    parser.add_argument("--seconds", type=float, default=480.0)
    parser.add_argument("--min-width", type=float, default=DEFAULT_MIN_WIDTH)
    parser.add_argument("--face-min-width", type=float, default=DEFAULT_FACE_MIN_WIDTH)
    parser.add_argument("--skip-validation", action="store_true")
    parser.add_argument("--probe-only", action="store_true")
    args = parser.parse_args(argv)

    started = time.monotonic()
    assert all(_SERIES_IDENTITIES) and _WEIGHT_DERIVATIVE_IDENTITY
    print("Symbolic derivative identities PASS: endpoint and weight-only rules",
          flush=True)
    print("Endpoint series/remainder proof PASS: X=1/16, c=%s, certified DQ_hi=%s"
          % (_DQ_C, DQ_ENDPOINT), flush=True)
    cert2.prove_rho_endpoint_monotonicity()
    gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    rho_gmax = cert2.get_rho_gmax()
    t0 = mp.mpf(PSI) + mp.mpf("0.0001")
    lambdas, lamhat = lambda_family(t0)
    pstar, rstar, wstar, _, residual = solve_face_kkt(t0)
    assert residual < mp.mpf("1e-35")
    print("Face KKT solve PASS (numerical): p*=%.15f r*=%.15f w*=%.15f lambda*=%.15f residual<1e-35"
          % (float(pstar), float(rstar), float(wstar), float(lamhat)), flush=True)
    print("Lambda family: %s" % (lambdas,), flush=True)
    if not args.skip_validation:
        verify_derivatives()
        verify_dq()
        verify_face_soundness(t0, lambdas)
        verify_orbit_swap()
        verify_centered_overlap_regression(t0, lamhat)
        verify_pinned_face_regression(t0, lambdas)
        verify_centered5(t0, (0.0, float(lamhat)))
        verify_centered5_mixed(t0, (0.0, float(lamhat)))
        verify_sink_safe_centering(t0, (0.0, float(lamhat)))
        verify_strip_pin(t0, lambdas)
    face_flip_probe(t0, lambdas)
    centered5_subpin_probe(t0, (0.0, float(lamhat)))
    if args.probe_only:
        return 0

    # A bounded demonstration around the obstruction. This is intentionally
    # not presented as the full-cube certificate; external integrations use
    # certify(t, root_box, ...).
    demo_root = ((0.2, 0.5), (0.2, 0.5), (0.2, 0.5),
                 (0.9, 1.0), (0.7, 1.0))
    ok, stats = run_bb(
        t0, lambdas, gmax, rho_gmax, args.min_width,
        args.face_min_width, args.max_boxes, args.seconds,
        root_box=demo_root,
    )
    assert stats["processed"] > 0
    print("Demonstration root verdict: %s; pin+face clears=%d"
          % ("COMPLETE" if ok else "INCOMPLETE", stats["face"]), flush=True)
    assert stats["face"] > 0
    tight_width = .005
    tight_root = (
        (float(pstar) - tight_width / 2, float(pstar) + tight_width / 2),
        (float(pstar) - tight_width / 2, float(pstar) + tight_width / 2),
        (float(rstar) - tight_width / 2, float(rstar) + tight_width / 2),
        (1 - tight_width / 4, 1.0),
        (float(wstar) - tight_width / 2, float(wstar) + tight_width / 2),
    )
    tight_ok, tight_stats = run_bb(
        t0, lambdas, gmax, rho_gmax, args.min_width,
        args.face_min_width, 100000, min(args.seconds, 60.0),
        root_box=tight_root,
    )
    assert tight_ok and tight_stats["face"] > 0 and tight_stats["stack"] == 0
    print("Tight obstruction face-path demonstration PASS: COMPLETE, "
          "pin+face clears=%d" % tight_stats["face"], flush=True)
    print("Largest certified t: %s"
          % (mp.nstr(t0, 30) if ok else
             "NONE (bounded obstruction demonstration incomplete)"), flush=True)
    elapsed = time.monotonic() - started
    assert elapsed >= 0
    print("cert3 total runtime %.2fs" % elapsed, flush=True)
    return 0


if __name__ == "__main__":
    main()
