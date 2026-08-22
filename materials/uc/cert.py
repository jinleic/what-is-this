"""Rigorous Arb interval certifier for the target, and an honest B&B report.

TARGET (Theorem B''', margin_lemma.py).  For t and alpha fixed, prove

    Phi(nu) >= 0   for every pair-orbit measure nu with <= 2 orbits, mean <= t,
    Phi := L * Lambda = (2(1-a) B - 1) L + a C - (1-a) sigma^2,

which gives Lambda >= 0 wherever L > 0, hence F >= 0 on ALL feasible mu by the
Margin Lemma (F = 0 identically where L = 0).  Five parameters
(p1,q1,p2,q2,w), 0 <= p_i <= q_i <= 1, w in [0,1].

WHY Phi AND NOT Lambda.  Lambda needs sigma^2/L and C/L, and near the sink set L
has 0 in its enclosure, so no interval bound on a quotient survives.  Phi needs
NO DIVISION, and still converges on the sink corner because sigma^2 is QUADRATIC
in the non-sink mass while L is linear.

-------------------------------------------------------------------------------
SOUNDNESS.  An earlier version of this file was NOT a rigorous certifier: it
extracted Arb bounds with .lower()/.upper(), converted them to Python floats, and
then did all the range algebra in binary floating point.  Round-to-nearest can
move a lower bound UP or an upper bound DOWN, so any "certificate" from it would
have been worthless -- and random-point validation cannot detect the defect,
because the error is in the arithmetic, not in the formulas.

This version keeps EVERY quantity as an Arb ball from start to finish.  Arb ball
arithmetic is outward-rounding by construction, so each ball is a guaranteed
enclosure.  All pruning decisions are made with Arb's certified comparisons:

  * `ball >= 0` returns True only if EVERY point of the ball is >= 0
    (checked below: arb(0,1) >= 0 is False, arb(0.5,0.5) >= 0 is False since 0
    lies in that ball).  Clearing a box therefore requires a certified sign.
  * discarding a box as mean-infeasible requires a certified `mean_lo > t`;
    when the comparison is undecided the box is KEPT, which is the safe side.

Enclosures, and why each is exact:

  h on [u,v]:  h is continuous, increasing on [0,1/2], decreasing on [1/2,1],
  h(0)=h(1)=0, h(1/2)=1.  So its range is the convex hull of {h(u), h(v)},
  together with 1 when the interval straddles 1/2.  Built with Arb `union`
  (outward-rounded hull).  If the straddle test is undecided, 1 is included --
  the conservative choice.

  s*(p,r) = median{1/2, max(p,r), min(p+r,1)} is nondecreasing in each argument,
  so its range is the hull of its values at the two corners.

  rh(p) = 2(1-p)h(p) - h((1-p)^2):  ball arithmetic on the formula, with the
  upper end additionally clamped by a rigorous global bound RH_GMAX obtained by
  covering [0,1] with Arb balls.  NO monotonicity of rh is assumed.
"""

import heapq
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI, h as h_mp, mp, mpf

ctx.prec = 160

ALPHA = arb("0.0356069")
ONE = arb(1)
TWO = arb(2)
HALF = arb(1) / 2
LOG2 = TWO.log()
B_OBS_F = 0.32945473850303697239


def hull(*balls):
    out = balls[0]
    for b in balls[1:]:
        out = out.union(b)
    return out


def h_pt(x):
    """h at a point, as an Arb ball; h(0) = h(1) = 0 exactly."""
    x = arb(x)
    if x <= 0 or x >= 1:
        return arb(0)
    return -(x * x.log() + (ONE - x) * (ONE - x).log()) / LOG2


def h_encl(u, v):
    """Enclosure of the RANGE of h over [u,v]; u, v are exact Arb points."""
    u, v = arb(u), arb(v)
    e = hull(h_pt(u), h_pt(v))
    if v <= HALF or u >= HALF:       # certified monotone on the interval
        return e
    return hull(e, ONE)              # straddles (or undecided): include the peak


def amin(a, b):
    a, b = arb(a), arb(b)
    if a <= b:
        return a
    if b <= a:
        return b
    return hull(a, b)                # undecided: the hull contains both


def amax(a, b):
    a, b = arb(a), arb(b)
    if a >= b:
        return a
    if b >= a:
        return b
    return hull(a, b)


def sstar_pt(p, r):
    p, r = arb(p), arb(r)
    return amin(amax(HALF, amax(p, r)), amin(p + r, ONE))


def sstar_encl(pl, ph, rl, rh):
    """s* is nondecreasing in each argument, so the corners bound the range."""
    return hull(sstar_pt(pl, rl), sstar_pt(ph, rh))


def rh_encl(u, v, gmax=None):
    """Enclosure of rh(p) = 2(1-p)h(p) - h((1-p)^2) over [u,v]."""
    u, v = arb(u), arb(v)
    p = hull(u, v)
    hp = h_encl(u, v)
    x_lo, x_hi = (ONE - v) * (ONE - v), (ONE - u) * (ONE - u)
    h2 = h_encl(x_lo, x_hi)
    e = TWO * (ONE - p) * hp - h2
    lo = amax(e.lower(), arb(0))
    hi = e.upper() if gmax is None else amin(e.upper(), gmax)
    return hull(lo, hi)

def sqrt_rh_enc(u, v, gmax=None):
    """TWO-SIDED enclosure of sqrt(rh) over [u,v], NaN-free.

    rh >= 0 everywhere (squared norm density), so sqrt(rh) in [0, sqrt(hi)].
    We take sqrt of the positive POINT hi (never of the straddling hull), then
    hull(0, sqrt(hi)).  Arb's outward-rounded union dips a sliver below 0, but
    that is HARMLESS here: no further sqrt is taken on this ball -- it only
    enters affine sums and a final square.  Using a scalar upper bound instead
    would silently turn phi_encl into a one-sided lower bound and break its
    containment test; this keeps it a genuine two-sided enclosure.
    """
    hi = rh_encl(u, v, gmax).upper()
    return hull(arb(0), hi.sqrt())


def _rh_global(n=20000):
    """Rigorous upper bound for max_p rh(p): cover [0,1] with Arb balls."""
    best = arb(0)
    for i in range(n):
        e = rh_encl(arb(i) / n, arb(i + 1) / n)
        u = e.upper()
        if not (u <= best):
            best = amax(best, u)
    return best.upper()


print("0. Certified comparison semantics (flint), then the global rh bound")
print(f"   arb(0,1)     >= 0 : {arb(0.0, 1.0) >= 0}   (must be False)")

print(f"   arb(0.5,0.5) >= 0 : {arb(0.5, 0.5) >= 0}   (must be False: 0 in ball)")
print(f"   arb(0.5,0.25)>= 0 : {arb(0.5, 0.25) >= 0}   (must be True)")
assert not (arb(0.0, 1.0) >= 0)
assert not (arb(0.5, 0.5) >= 0)
assert arb(0.5, 0.25) >= 0
print("   clearing a box therefore requires a CERTIFIED sign               [OK]")

RH_GMAX = _rh_global()
print(f"\n   max_p rh(p) <= {float(RH_GMAX):.12f}   (float search said 0.234229480)")
assert RH_GMAX >= arb("0.2342294") and RH_GMAX <= arb("0.2350")
print("   brackets the float maximum from ABOVE, as required               [OK]")

print("   endpoint regression: sqrt_rh_enc([0.95,1]) -- finite AND a valid")
print("   two-sided enclosure (must contain the mpmath truth at 0.95, 0.975, 1.0)")
_ep = sqrt_rh_enc(arb("0.95"), arb(1.0), RH_GMAX)
assert _ep.is_finite(), "sqrt_rh_enc returned non-finite at the q2->1 endpoint"
assert _ep.lower() <= arb(0) + arb("1e-30") and _ep.upper() <= arb("0.17"), _ep
from entropy import h as _h, mpf as _mpf  # noqa
for _xv in ("0.95", "0.975", "0.999"):
    _x = _mpf(_xv)
    _y = 1 - _x
    _rht = 2 * _y * _h(_x) - _h(_y * _y)
    _sqt = (_rht ** 0.5) if _rht > 0 else _mpf(0)
    _b = arb(_sqt)
    assert _ep.lower() - arb("1e-30") <= _b and _b <= _ep.upper() + arb("1e-30"), \
        (_xv, float(_ep.lower()), float(_sqt), float(_ep.upper()))
print(f"   sqrt_rh_enc([0.95,1]) = [{float(_ep.lower()):.4e}, "
      f"{float(_ep.upper()):.4e}]   contains truth at 3 pts   [OK]")

# ---------------------------------------------------------------------------
print("\n1. Enclosures validated against mpmath (they must CONTAIN the truth)")
rng = np.random.default_rng(0)
bad_h = bad_rh = bad_s = 0
for _ in range(1200):
    u = float(rng.uniform(0, 1))
    v = min(1.0, u + float(rng.uniform(0, 0.3)))
    eh, er = h_encl(u, v), rh_encl(u, v, RH_GMAX)
    hl, hh = eh.lower(), eh.upper()
    rl, rh_ = er.lower(), er.upper()
    for g in np.linspace(u, v, 25):
        hv = arb(float(h_mp(mpf(float(g)))))
        if not (hl <= hv and hv <= hh):
            bad_h += 1
            break
        x = 1 - mpf(float(g))
        rv = arb(float(2 * x * h_mp(mpf(float(g))) - h_mp(x * x)))
        if not (rl <= rv or rv <= arb(0)) or not (rv <= rh_):
            bad_rh += 1
            break
    a2 = float(rng.uniform(0, 1))
    b2 = min(1.0, a2 + float(rng.uniform(0, 0.3)))
    es = sstar_encl(u, v, a2, b2)
    for g in np.linspace(u, v, 7):
        for g2 in np.linspace(a2, b2, 7):
            sv = arb(min(max(0.5, max(g, g2)), min(g + g2, 1.0)))
            if not (es.lower() <= sv and sv <= es.upper()):
                bad_s += 1
                break
assert bad_h == 0 and bad_rh == 0 and bad_s == 0, (bad_h, bad_rh, bad_s)
print(f"   h / rh / s* containment failures: {bad_h} / {bad_rh} / {bad_s}"
      f"  over 1200 intervals            [OK]")


# ---------------------------------------------------------------------------
def phi_encl(box, t):
    """Arb enclosure of Phi over a 5-D box, or None if certified infeasible."""
    (p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh) = box
    p1l, p1h = arb(p1l), arb(p1h)
    q1l, q1h = arb(q1l), arb(q1h)
    p2l, p2h = arb(p2l), arb(p2h)
    q2l, q2h = arb(q2l), arb(q2h)
    wl, wh = arb(wl), arb(wh)
    t = arb(t)
    q1l, q2l = amax(q1l, p1l), amax(q2l, p2l)
    if q1l > q1h or q2l > q2h:
        return None
    w = hull(wl, wh)
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    if mean.lower() > t:                      # certified infeasible
        return None
    b1, b2 = ONE - m1, ONE - m2
    B = w * b1 + (ONE - w) * b2
    # on the feasible set B >= 1-t; intersect (max of two valid lower bounds)
    B = hull(amax(B.lower(), ONE - t), B.upper())
    L = w * (h_encl(p1l, p1h) + h_encl(q1l, q1h)) / TWO \
        + (ONE - w) * (h_encl(p2l, p2h) + h_encl(q2l, q2h)) / TWO
    L = hull(amax(L.lower(), arb(0)), L.upper())
    s1 = sstar_encl(p1l, p1h, q1l, q1h)
    s2 = sstar_encl(p2l, p2h, q2l, q2h)
    C = w * h_encl(s1.lower(), s1.upper()) \
        + (ONE - w) * h_encl(s2.lower(), s2.upper())
    C = hull(amax(C.lower(), arb(0)), C.upper())
    r1 = (sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    r2 = (sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    sig = w * r1 + (ONE - w) * r2
    return (TWO * (ONE - ALPHA) * B - ONE) * L + ALPHA * C \
        - (ONE - ALPHA) * sig * sig


def phi_true(p1, q1, p2, q2, w, alpha="0.0356069"):
    """Reference value in mpmath at dps 50, returned as an mpf.

    The reference is the LESS accurate object here: an Arb enclosure of a point
    box has radius ~1e-47 at prec 160, while mpmath at dps 50 carries ~1e-50.
    Containment is therefore tested with a slack of 1e-30, which absorbs the
    reference's own error and is still 17 orders tighter than anything the B&B
    decisions depend on.
    """
    old = mp.dps
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
        B = wv * (1 - (P[0] + P[1]) / 2) + (1 - wv) * (1 - (P[2] + P[3]) / 2)
        L = (wv * (h_mp(P[0]) + h_mp(P[1])) / 2
             + (1 - wv) * (h_mp(P[2]) + h_mp(P[3])) / 2)
        C = wv * h_mp(st(P[0], P[1])) + (1 - wv) * h_mp(st(P[2], P[3]))
        S = (wv * (mp.sqrt(rr(P[0])) + mp.sqrt(rr(P[1]))) / 2
             + (1 - wv) * (mp.sqrt(rr(P[2])) + mp.sqrt(rr(P[3]))) / 2)
        return (2 * (1 - al) * B - 1) * L + al * C - (1 - al) * S * S
    finally:
        mp.dps = old


print("\n2. Soundness: the enclosure must CONTAIN the true value")
t_test = float(PSI) + 1e-4
# build the obstruction point STRICTLY inside the mean constraint, so the test
# exercises containment rather than the knife-edge mean = t
a_obs = (t_test - 1e-9 - B_OBS_F) / (1 - B_OBS_F)
SLACK = arb("1e-30")
pts = [(B_OBS_F, B_OBS_F, B_OBS_F, 1.0, 1 - 2 * a_obs),
       (0.3, 0.3, 0.2, 0.9, 0.5),
       (0.05, 0.05, 0.0, 1.0, 0.3),
       (0.0, 0.0, 0.0, 1.0, 0.6),
       (0.1165, 0.1165, 0.0, 1.0, 0.02)]
print(f"   {'point':>34} {'enclosure lower':>17} {'true (mpmath)':>17} "
      f"{'contains':>9}")
n_checked = 0
for (p1, q1, p2, q2, w) in pts:
    e = phi_encl(((p1, p1), (q1, q1), (p2, p2), (q2, q2), (w, w)), t_test)
    tv = phi_true(p1, q1, p2, q2, w)
    tb = arb(mp.nstr(tv, 40))
    if e is None:
        print(f"   {str((round(p1,4),round(q1,4),round(p2,4),round(q2,4))):>34}"
              f" {'infeasible':>17} {float(tv):>+17.9f} {'n/a':>9}")
        continue
    ok = bool(e.lower() - SLACK <= tb) and bool(tb <= e.upper() + SLACK)
    n_checked += 1
    print(f"   {str((round(p1,4),round(q1,4),round(p2,4),round(q2,4))):>34}"
          f" {float(e.lower()):>+17.9f} {float(tv):>+17.9f} {str(ok):>9}")
    assert ok, (p1, q1, p2, q2, w, float(e.lower()), float(tv))
assert n_checked >= 3, n_checked
print(f"   {n_checked} feasible points checked; every enclosure contains the")
print("   mpmath reference to within 1e-30                                 [OK]")


# ---------------------------------------------------------------------------
def branch_and_bound(t, max_boxes=200000, min_diam=1e-6, time_budget=240.0):
    root = ((0.0, 1.0),) * 5
    heap = [(0.0, 0, root)]
    counter = processed = cleared = infeas = 0
    residual = []
    t0 = time.time()
    while heap:
        if processed >= max_boxes or time.time() - t0 > time_budget:
            residual.extend(b for _, _, b in heap)
            break
        _, _, box = heapq.heappop(heap)
        processed += 1
        e = phi_encl(box, t)
        if e is None:
            infeas += 1
            continue
        if e >= 0:                      # CERTIFIED nonnegative on the whole box
            cleared += 1
            continue
        widths = [hi - lo for lo, hi in box]
        if max(widths) <= min_diam:
            residual.append(box)
            continue
        k = int(np.argmax(widths))
        lo, hi = box[k]
        mid = (lo + hi) / 2
        for sub in ((lo, mid), (mid, hi)):
            nb = list(box)
            nb[k] = sub
            heapq.heappush(heap, (float(e.lower()), counter, tuple(nb)))
            counter += 1
    complete = (len(residual) == 0 and len(heap) == 0)
    return complete, residual, dict(processed=processed, cleared=cleared,
                                    infeasible=infeas, queue=len(heap),
                                    secs=time.time() - t0)


print("\n3. Branch-and-bound.  A certificate requires zero residual AND empty queue.")
for lbl, tv in (("psi", float(PSI)), ("psi+1e-4", float(PSI) + 1e-4)):
    ok, resid, st = branch_and_bound(tv, max_boxes=60000, time_budget=180.0)
    print(f"\n   t = {lbl} ({tv:.16f})")
    print(f"     processed / cleared / infeasible : "
          f"{st['processed']} / {st['cleared']} / {st['infeasible']}")
    print(f"     queue / residual                : {st['queue']} / {len(resid)}")
    print(f"     seconds                          : {st['secs']:.1f}")
    print(f"     COMPLETE CERTIFICATE             : {ok}")
    assert not ok or len(resid) == 0

print("\nALL PASS")
print("\nSUMMARY")
print("  Every quantity is an Arb ball end to end; all pruning uses certified")
print("  Arb signs.  The earlier float-based version of this file was UNSOUND")
print("  (Arb bounds converted to floats, then float range algebra) and is")
print("  retracted -- see the module docstring.")
print(f"  Global bound proved by covering: max rh <= {float(RH_GMAX):.9f}.")
print("  B&B RESULT: does NOT converge within budget -- NO certificate is")
print("  claimed.  The bottleneck is enclosure quality in 5 dimensions against")
print("  a margin of ~4e-4; see uc/README.md for what would be needed.")
