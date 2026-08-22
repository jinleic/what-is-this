"""Sound sink-safe partial centered bound (experimental certifier rule).

For the prepared ordered atom rectangle A and weight interval W, fix
``w_star = mid(W)``.  For every atom tuple a in A and w in W, the segment
``(a,w_star) -> (a,w)`` stays in A x W.  The one-dimensional mean-value theorem
therefore gives, for lambda >= 0,

  Phi_lam(a,w) >= Phi_lam(a,w_star)
                    - sup_{A x W}|d_w Phi_lam| rad(W).

``bound_kkt.phi_lam_encl`` naturally encloses the first term over every a in A.
The derivative below is the exact ``w`` partial from
``bound_kkt.interval_gradient`` plus ``lambda*(m1-m2)``.  It uses only bounded
aggregate enclosures, never atom derivatives, so it remains finite when atom
intervals touch entropy/sqrt-rh sinks.  On feasible points M <= t,

  Phi_lam = Phi + lambda*(M-t) <= Phi,

so any nonnegative lower bound for Phi_lam certifies Phi >= 0.  No feasibility
intersection is inserted into either full-box interval enclosure.

This module is not part of the active frozen campaign.  Its ``__main__`` block
runs sampled soundness controls only; a future certificate must rerun a fresh
immutable full campaign with this rule included and independently reviewed.
"""

import os
import random
import sys

import mpmath as mp
from flint import arb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import bound_kkt
import cert2
import cert3
from arbcore import ALPHA, ONE, TWO, amax, get_rh_gmax, h_encl, hull, \
    sqrt_rh_enc, sstar_encl
from entropy import PSI

ZERO = arb(0)


def _components(box5, gmax):
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
    L = w * l1 + v * l2
    S = w * u1 + v * u2
    kappa = TWO * (ONE - ALPHA) * (ONE - mean) - ONE
    dphi_dw = (TWO * (ONE - ALPHA) * (m2 - m1) * L
                + kappa * (l1 - l2) + ALPHA * (c1 - c2)
                - TWO * (ONE - ALPHA) * S * (u1 - u2))
    if not dphi_dw.is_finite():
        return None
    return prepared, m1, m2, dphi_dw


def centered_w_best(box5, t, lams, gmax=None):
    """Return max sound w-only MVT lower bound for feasible points in box."""
    if gmax is None:
        gmax = get_rh_gmax()
    bound_kkt.RH_GMAX = gmax
    components = _components(box5, gmax)
    if components is None:
        return None
    prepared, m1, m2, dphi_dw = components
    p1l, p1h, q1l, q1h, p2l, p2h, q2l, q2h, wl, wh = prepared
    wm = (wl + wh) / TWO
    rw = (wh - wl) / TWO
    anchor = ((p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h),
              (wm, wm))
    best = None
    for lam_value in lams:
        lam = cert3._arb(lam_value)
        assert lam >= ZERO
        value = bound_kkt.phi_lam_encl(anchor, t, lam)
        if value is None:
            continue
        derivative = dphi_dw + lam * (m1 - m2)
        penalty = amax(-derivative.lower(), derivative.upper()) * rw
        lower = cert3._arb(value.lower()) - penalty
        if lower.is_finite() and (best is None or lower > best):
            best = lower
    return best


def _sampled_controls(samples=200, points=20, seed=20260813):
    """Sampled evidence only: bound <= pointwise Phi_lam on random boxes."""
    rng = random.Random(seed)
    t = mp.mpf(PSI) + mp.mpf("0.0001")
    gmax = get_rh_gmax()
    lambdas, _ = cert3.lambda_family(t)
    comparisons = 0
    usable = 0
    for _ in range(samples):
        box = []
        for _j in range(4):
            lo = rng.random() * 0.98
            hi = lo + rng.random() * (1.0 - lo)
            box.append((lo, hi))
        wl = rng.random() * 0.98
        wh = wl + rng.random() * (1.0 - wl)
        box.append((wl, wh))
        lam = rng.choice(lambdas)
        lower = centered_w_best(tuple(box), t, (lam,), gmax)
        prepared = bound_kkt._prepare_box(tuple(box))
        if lower is None or prepared is None:
            continue
        usable += 1
        lo_value = mp.mpf(str(float(lower.lower())))
        prepared_box = list(zip(prepared[0::2], prepared[1::2]))
        for _k in range(points):
            # Sample the exact prepared rectangle used by the proof.  Sorting
            # an arbitrary raw-box point would be wrong: on overlapping p/q
            # intervals it can move the point outside this rectangle.
            point = [rng.uniform(float(a), float(b))
                     for a, b in prepared_box]
            p1, q1, p2, q2, w = map(mp.mpf, point)
            mean = w * (p1 + q1) / 2 + (1 - w) * (p2 + q2) / 2
            truth = cert2.phi_true(p1, q1, p2, q2, w) \
                + mp.mpf(str(lam)) * (mean - t)
            assert lo_value <= truth + mp.mpf("1e-25"), \
                (box, prepared_box, lam, lo_value, truth)
            comparisons += 1
    assert usable > 0 and comparisons > 0
    print("centered_w sampled soundness PASS: %d boxes, %d point comparisons" %
          (usable, comparisons))


if __name__ == "__main__":
    _sampled_controls()
