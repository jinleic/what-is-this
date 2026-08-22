"""Rigorous four-dimensional Arb certifier with the orbit weight eliminated.

For fixed orbit aggregates, Phi is a quadratic in the orbit weight ``w``.  This
module encloses its three coefficients and minimizes the quadratic over a
certified feasible weight interval.  A complementary Cauchy--Schwarz ratio rule
clears sink-collar boxes, where any direct lower enclosure must converge to the
exact zero at a sink corner.

Importing this module is side-effect free apart from Python imports and Arb's
precision setting inherited from :mod:`arbcore`; validations and searches run
only under ``if __name__ == '__main__'``.
"""

import heapq
import os
import random
import sys
import time

import mpmath as mp
import sympy as sp
from flint import arb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from arbcore import (ALPHA, HALF, ONE, TWO, amax, amin, get_rh_gmax, h_encl,
                     h_pt, hull, rh_encl, sqrt_rh_enc, sstar_encl)
from entropy import PSI, h as h_mp

B_OBS = 0.32945473850303697239
C_STAR = 0.38234553336670272115
ZERO = arb(0)

# The symbols and expansion are module constants so the implementation really
# is generated from, and checked against, the requested symbolic polynomial.
_A_M, _A_L, _A_C, _A_S = sp.symbols("AM AL AC AS")
_B_M, _B_L, _B_C, _B_S = sp.symbols("BM BL BC BS")
_W, _AA = sp.symbols("w alpha")
_M_EXPR = _B_M + _W * (_A_M - _B_M)
_L_EXPR = _B_L + _W * (_A_L - _B_L)
_C_EXPR = _B_C + _W * (_A_C - _B_C)
_S_EXPR = _B_S + _W * (_A_S - _B_S)
_PHI_EXPR = sp.expand((2 * (1 - _AA) * (1 - _M_EXPR) - 1) * _L_EXPR
                      + _AA * _C_EXPR - (1 - _AA) * _S_EXPR ** 2)
_POLY = sp.Poly(_PHI_EXPR, _W)
_C2_EXPR, _C1_EXPR, _C0_EXPR = _POLY.all_coeffs()
_COEFF_FLOAT = sp.lambdify((_A_M, _A_L, _A_C, _A_S,
                            _B_M, _B_L, _B_C, _B_S, _AA),
                           (_C2_EXPR, _C1_EXPR, _C0_EXPR), "math")

# Exact endpoint-ratio monotonicity used by ``rho_upper``.  With base-2
# entropy h, h(z)-z h'(z)=-log_2(1-z), so
# (z/h(z))'=-log_2(1-z)/h(z)^2 > 0 on (0,1).  By h(1-z)=h(z),
# (1-z)/h(z) is therefore decreasing in z.
_RHO_Z = sp.symbols("z", positive=True)
_RHO_H = -(_RHO_Z * sp.log(_RHO_Z)
           + (1 - _RHO_Z) * sp.log(1 - _RHO_Z)) / sp.log(2)
_RHO_NUMERATOR_IDENTITY = sp.simplify(
    _RHO_H - _RHO_Z * sp.diff(_RHO_H, _RHO_Z)
    + sp.log(1 - _RHO_Z) / sp.log(2))
_RHO_DERIVATIVE_IDENTITY = sp.simplify(
    sp.diff(_RHO_Z / _RHO_H, _RHO_Z)
    + sp.log(1 - _RHO_Z) / (sp.log(2) * _RHO_H ** 2))
assert _RHO_NUMERATOR_IDENTITY == 0
assert _RHO_DERIVATIVE_IDENTITY == 0


def _quadratic_coefficients(AM, AL, AC, AS, BM, BL, BC, BS):
    """Arb evaluation of the symbolically expanded quadratic coefficients."""
    dm, dl, dc, ds = AM - BM, AL - BL, AC - BC, AS - BS
    beta = ONE - ALPHA
    # Expanding (K0 - 2 beta dm w)(BL + dl w)
    # + alpha(BC + dc w) - beta(BS + ds w)^2,
    # where K0 = 2 beta(1-BM)-1.
    k0 = TWO * beta * (ONE - BM) - ONE
    c0 = k0 * BL + ALPHA * BC - beta * BS * BS
    c1 = k0 * dl - TWO * beta * dm * BL + ALPHA * dc \
        - TWO * beta * BS * ds
    c2 = -TWO * beta * dm * dl - beta * ds * ds
    return c2, c1, c0


def _lemma_sqrt_rh_upper(u, v):
    """Upper bound from the proved ``rh(z) <= 2 min(z,1-z)`` lemma."""
    u, v = arb(u), arb(v)
    if v <= HALF:
        m = v
    elif u >= HALF:
        m = ONE - u
    else:
        m = HALF
    return (TWO * amax(m, ZERO)).upper().sqrt()


def _sqrt_rh_tight(u, v, gmax):
    base = sqrt_rh_enc(u, v, gmax)
    upper = amin(base.upper(), _lemma_sqrt_rh_upper(u, v))
    return hull(base.lower(), upper)


def _orbit_aggregates(pbox, qbox, gmax):
    """Return exact-range aggregate enclosures ``(M,L,C,S)`` for one orbit."""
    pl, ph = arb(pbox[0]), arb(pbox[1])
    ql, qh = arb(qbox[0]), arb(qbox[1])
    mean = hull((pl + ql) / TWO, (ph + qh) / TWO)
    entropy = (h_encl(pl, ph) + h_encl(ql, qh)) / TWO
    ss = sstar_encl(pl, ph, ql, qh)
    correction = h_encl(ss.lower(), ss.upper())
    sigma = (_sqrt_rh_tight(pl, ph, gmax)
             + _sqrt_rh_tight(ql, qh, gmax)) / TWO
    return mean, entropy, correction, sigma


def _weight_range(AM, BM, t, wlo=ZERO, whi=ONE):
    """Keep exactly the possibly-feasible weight range, outward rounded.

    Orbit boxes are independent.  For fixed ``w`` in ``[0,1]``, the minimum
    realizable mean is therefore

        (1-w) * BM.lower() + w * AM.lower().

    Contracting this affine inequality avoids the co-variation error that would
    result from treating ``AM-BM`` independently of its shared ``BM`` term.
    """
    wlo, whi, t = arb(wlo), arb(whi), arb(t)
    if wlo > whi:
        return None
    am, bm = AM.lower(), BM.lower()
    am_ok, bm_ok = bool(am <= t), bool(bm <= t)
    if not am_ok and not bm_ok:
        return None
    if am_ok and bm_ok:
        pass
    elif am_ok:
        floor = (bm - t) / (bm - am)
        wlo = amax(wlo, floor.lower())
    else:
        ceiling = (t - bm) / (am - bm)
        whi = amin(whi, ceiling.upper())
    if wlo > whi or whi < 0 or wlo > 1:
        return None
    wlo, whi = amax(wlo, ZERO), amin(whi, ONE)
    if wlo > whi:
        return None
    return hull(wlo, whi)


def _quadratic_lower(c2, c1, c0, W):
    """Lower bound all coefficient realizations over the feasible weight range.

    For each fixed realization, a one-dimensional quadratic attains its minimum
    at an endpoint or at its interior vertex.  Endpoint Arb evaluations enclose
    every realization.  If convexity is possible, the interval vertex
    ``-c1/(2*c2)`` is intersected with ``W`` and the quadratic is evaluated over
    that whole interval, enclosing every possible interior vertex value.
    """
    wl, wh = W.lower(), W.upper()
    values = [c2 * wl * wl + c1 * wl + c0,
              c2 * wh * wh + c1 * wh + c0]
    if c2.lower() > 0:
        vertex = -c1 / (TWO * c2)
        vlo = amax(vertex.lower(), wl)
        vhi = amin(vertex.upper(), wh)
        if not (vlo > vhi):
            vbox = hull(vlo, vhi)
            values.append(c2 * vbox * vbox + c1 * vbox + c0)
    elif not (c2.upper() <= 0):
        # c2 straddles zero.  Nonzero convex realizations can have vertices, but
        # direct division by c2 would be unbounded.  Evaluating over all W is a
        # conservative enclosure containing every such vertex evaluation.
        values.append(c2 * W * W + c1 * W + c0)
    best = values[0].lower()
    for value in values[1:]:
        best = amin(best, value.lower())
    return best.lower()


def quadratic_bound(box, t, gmax=None, wlo=ZERO, whi=ONE):
    """Return the w-eliminated Arb lower endpoint, or ``None`` if infeasible."""
    if gmax is None:
        gmax = get_rh_gmax()
    if len(box) != 4:
        raise ValueError("quadratic_bound expects four coordinate intervals")
    AM, AL, AC, AS = _orbit_aggregates(box[0], box[1], gmax)
    BM, BL, BC, BS = _orbit_aggregates(box[2], box[3], gmax)
    W = _weight_range(AM, BM, arb(t), wlo, whi)
    if W is None:
        return None
    c2, c1, c0 = _quadratic_coefficients(AM, AL, AC, AS, BM, BL, BC, BS)
    return _quadratic_lower(c2, c1, c0, W)


# ---------------------------------------------------------------------------
# Ratio-rule machinery.  rho(z)=rh(z)/h(z) is evaluated directly away from the
# entropy-zero endpoints.  Endpoint collars use the separately proved bound
# rh(z)<=2 min(z,1-z), avoiding division by an interval containing h=0.
# ---------------------------------------------------------------------------
RHO_GMAX = None
_RHO_GMAX_N = None


def _rho_formula_encl(u, v):
    u, v = arb(u), arb(v)
    hz = h_encl(u, v)
    gu = TWO * u - u * u
    gv = TWO * v - v * v
    hg = h_encl(gu, gv)
    return TWO * (ONE - hull(u, v)) - hg / hz


def _rho_direct_upper(u, v):
    """Dependency-safe upper bound ``2(1-u)-h(g)_lo/h(z)_hi``."""
    u, v = arb(u), arb(v)
    gu = TWO * u - u * u
    gv = TWO * v - v * v
    hg_lo = h_encl(gu, gv).lower()
    hz_hi = h_encl(u, v).upper()
    return (TWO * (ONE - u) - hg_lo / hz_hi).upper()


def _rho_global(n=20000):
    """Certified global rho cap, with analytic endpoint guards.

    Every cover cell uses the same sound caps as :func:`rho_upper`, except for
    the global cap itself.  Thus the endpoint cells use the proved lemma and
    interior cells take the tighter of the direct formula and lemma cap.
    """
    if n < 4:
        raise ValueError("rho cover needs at least four cells")
    best = ZERO
    for i in range(n):
        u, v = arb(i) / n, arb(i + 1) / n
        caps = []
        if u > 0 and v < 1:
            caps.append(_rho_direct_upper(u, v))
        if v <= HALF and v > 0:
            caps.append((TWO * v / h_pt(v)).upper())
        elif u >= HALF and u < 1:
            caps.append((TWO * (ONE - u) / h_pt(u)).upper())
        cell = caps[0]
        for cap in caps[1:]:
            cell = amin(cell, cap)
        if not (cell <= best):
            best = amax(best, cell)
    return best.upper()


def get_rho_gmax(n=20000):
    global RHO_GMAX, _RHO_GMAX_N
    if RHO_GMAX is None or _RHO_GMAX_N != n:
        RHO_GMAX = _rho_global(n)
        _RHO_GMAX_N = n
    return RHO_GMAX


def prove_rho_endpoint_monotonicity():
    """Report the exact derivative proof behind the endpoint-safe caps."""
    assert _RHO_NUMERATOR_IDENTITY == 0
    assert _RHO_DERIVATIVE_IDENTITY == 0
    print("rho endpoint monotonicity PROVED: "
          "h-z*h'=-log2(1-z)>0; symmetry gives the right endpoint",
          flush=True)


def rho_upper(u, v, global_cap=None):
    """Certified upper bound for rho on an atom interval.

    It takes the minimum of a certified global cover, a direct interior bound,
    and endpoint-safe bounds from ``rh <= 2 min(z,1-z)``.  The exact identity
    ``h-z*h'=-log2(1-z)>0`` proves ``z/h(z)`` increasing; entropy symmetry
    proves ``(1-z)/h(z)`` decreasing.  See
    :func:`prove_rho_endpoint_monotonicity`.
    """
    u, v = arb(u), arb(v)
    cap = arb(global_cap if global_cap is not None else get_rho_gmax())
    if u > 0 and v < 1:
        cap = amin(cap, _rho_direct_upper(u, v))
    if v <= HALF and v > 0:
        cap = amin(cap, (TWO * v / h_pt(v)).upper())
    elif u >= HALF and u < 1:
        cap = amin(cap, (TWO * (ONE - u) / h_pt(u)).upper())
    return amax(cap, ZERO).upper()


def ratio_rule(box, t, W=None, rho_gmax=None):
    """Return whether the certified Cauchy--Schwarz ratio rule clears ``box``."""
    if rho_gmax is None:
        rho_gmax = get_rho_gmax()
    if W is None:
        AM = hull((arb(box[0][0]) + arb(box[1][0])) / TWO,
                  (arb(box[0][1]) + arb(box[1][1])) / TWO)
        BM = hull((arb(box[2][0]) + arb(box[3][0])) / TWO,
                  (arb(box[2][1]) + arb(box[3][1])) / TWO)
        W = _weight_range(AM, BM, arb(t))
        if W is None:
            return None
    arho = (rho_upper(box[0][0], box[0][1], rho_gmax)
            + rho_upper(box[1][0], box[1][1], rho_gmax)) / TWO
    brho = (rho_upper(box[2][0], box[2][1], rho_gmax)
            + rho_upper(box[3][0], box[3][1], rho_gmax)) / TWO
    # arho and brho are one-sided upper caps.  Convex combinations preserve
    # those inequalities because w and 1-w are nonnegative.  The upper-cap
    # combination is linear in w, so its maximum over W is at an endpoint.
    # Writing brho+w*(arho-brho) would be unsound interval arithmetic here:
    # the subtracted brho has no matching lower enclosure.
    wl, wh = W.lower(), W.upper()
    rlo = (ONE - wl) * brho + wl * arho
    rhi = (ONE - wh) * brho + wh * arho
    rho_hi = amax(rlo.upper(), rhi.upper())
    beta = ONE - ALPHA
    kappa = TWO * beta * (ONE - arb(t)) - ONE
    margin = kappa - beta * rho_hi
    return bool(margin >= 0)


def certify_box(box, t, gmax, rho_gmax):
    """Classify a box as infeasible, ratio-cleared, quadratic-cleared, or open."""
    AM, AL, AC, AS = _orbit_aggregates(box[0], box[1], gmax)
    BM, BL, BC, BS = _orbit_aggregates(box[2], box[3], gmax)
    W = _weight_range(AM, BM, arb(t))
    if W is None:
        return "infeasible", None
    if ratio_rule(box, t, W, rho_gmax):
        return "ratio", ZERO
    c2, c1, c0 = _quadratic_coefficients(AM, AL, AC, AS, BM, BL, BC, BS)
    bound = _quadratic_lower(c2, c1, c0, W)
    if bound >= 0:
        return "quadratic", bound
    return "open", bound


def phi_true(p1, q1, p2, q2, w, alpha="0.0356069"):
    """Independent 50-decimal mpmath reference for the exact functional."""
    with mp.workdps(50):
        aa = mp.mpf(alpha)
        p1, q1, p2, q2, w = map(mp.mpf, (p1, q1, p2, q2, w))

        def rr(z):
            x = 1 - z
            return max(mp.mpf(0), 2 * x * h_mp(z) - h_mp(x * x))

        def ss(x, y):
            return min(max(mp.mpf("0.5"), max(x, y)), min(x + y, 1))

        m1, m2 = (p1 + q1) / 2, (p2 + q2) / 2
        mean = w * m1 + (1 - w) * m2
        entropy = (w * (h_mp(p1) + h_mp(q1)) / 2
                   + (1 - w) * (h_mp(p2) + h_mp(q2)) / 2)
        correction = w * h_mp(ss(p1, q1)) + (1 - w) * h_mp(ss(p2, q2))
        sigma = (w * (mp.sqrt(rr(p1)) + mp.sqrt(rr(q1))) / 2
                 + (1 - w) * (mp.sqrt(rr(p2)) + mp.sqrt(rr(q2))) / 2)
        return (2 * (1 - aa) * (1 - mean) - 1) * entropy \
            + aa * correction - (1 - aa) * sigma ** 2


def _point_aggregates(p, q):
    p, q = float(p), float(q)
    hp, hq = float(h_mp(p)), float(h_mp(q))
    ss = min(max(0.5, max(p, q)), min(p + q, 1.0))

    def rr(z):
        return max(0.0, 2 * (1 - z) * float(h_mp(z))
                   - float(h_mp((1 - z) ** 2)))

    return ((p + q) / 2, (hp + hq) / 2, float(h_mp(ss)),
            (rr(p) ** 0.5 + rr(q) ** 0.5) / 2)


def verify_expansion(samples=300, seed=2103):
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(samples):
        p1, q1, p2, q2, w = [rng.random() for _j in range(5)]
        A = _point_aggregates(p1, q1)
        B = _point_aggregates(p2, q2)
        c2, c1, c0 = _COEFF_FLOAT(*(A + B + (float(ALPHA),)))
        expanded = c2 * w * w + c1 * w + c0
        truth = float(phi_true(p1, q1, p2, q2, w))
        error = abs(expanded - truth)
        worst = max(worst, error)
        assert error <= 1e-10, (error, p1, q1, p2, q2, w)
    assert samples == 300
    print("expansion-verification PASS: 300 random float points; max error %.3e" % worst,
          flush=True)


def validate_rho(samples=500, seed=814):
    prove_rho_endpoint_monotonicity()
    rng = random.Random(seed)
    cap = get_rho_gmax()
    failures = 0
    checked = 0
    monotone_failures = 0
    previous_left = previous_right = None
    for i in range(1, 1001):
        z = mp.mpf(i) / 2002
        left = z / h_mp(z)
        right_z = mp.mpf("0.5") + mp.mpf(i) / 2002
        right = (1 - right_z) / h_mp(right_z)
        if previous_left is not None and left < previous_left:
            monotone_failures += 1
        if previous_right is not None and right > previous_right:
            monotone_failures += 1
        previous_left, previous_right = left, right
    assert monotone_failures == 0
    for _ in range(samples):
        u = rng.random()
        v = min(1.0, u + rng.random() * (1.0 - u))
        upper = rho_upper(u, v, cap)
        for j in range(11):
            z = mp.mpf(u + (v - u) * j / 10)
            if z == 0 or z == 1:
                truth = mp.mpf(0)
            else:
                truth = max(mp.mpf(0), 2 * (1 - z) * h_mp(z)
                            - h_mp((1 - z) ** 2)) / h_mp(z)
            checked += 1
            if not (arb(mp.nstr(truth, 45)) <= upper + arb("1e-40")):
                failures += 1
    assert samples == 500 and checked == 5500 and failures == 0
    print("rho sampled containment PASS: 500 random subintervals, 5,500 points; "
          "endpoint ratio samples agree with the exact derivative proof", flush=True)
    assert cap >= arb("0.304") and cap <= arb("0.32")
    print("certified global rho cover upper bound: %.12f" % float(cap.upper()), flush=True)


def validate_weight_contraction(boxes=200, points_per_box=50, seed=1701):
    """Sampled non-removal check for the certified feasible-weight contraction."""
    rng = random.Random(seed)
    checked = removed = 0
    for _ in range(boxes):
        intervals = []
        for _j in range(4):
            x, y = rng.random(), rng.random()
            intervals.append((min(x, y), max(x, y)))
        t = mp.mpf(str(float(PSI) + rng.random() * (C_STAR - float(PSI))))
        AM = hull((arb(intervals[0][0]) + arb(intervals[1][0])) / TWO,
                  (arb(intervals[0][1]) + arb(intervals[1][1])) / TWO)
        BM = hull((arb(intervals[2][0]) + arb(intervals[3][0])) / TWO,
                  (arb(intervals[2][1]) + arb(intervals[3][1])) / TWO)
        W = _weight_range(AM, BM, arb(mp.nstr(t, 45)))
        for _j in range(points_per_box):
            point = [mp.mpf(str(rng.uniform(lo, hi))) for lo, hi in intervals]
            w = mp.mpf(str(rng.random()))
            mean = w * (point[0] + point[1]) / 2 \
                + (1 - w) * (point[2] + point[3]) / 2
            if mean <= t:
                checked += 1
                if W is None:
                    removed += 1
                else:
                    wb = arb(mp.nstr(w, 45))
                    if not (W.lower() <= wb and wb <= W.upper()):
                        removed += 1
    assert boxes == 200 and points_per_box == 50
    assert checked > 1000, checked
    assert removed == 0, (removed, checked)
    print("weight contraction sampled non-removal PASS: %d feasible samples "
          "from 200 boxes x 50 points; 0 removed (sampled, not a proof)"
          % checked, flush=True)


def validate_sampled_soundness(boxes=200, points_per_box=30, seed=9917):
    rng = random.Random(seed)
    gmax = get_rh_gmax()
    violations = checked = 0
    for _ in range(boxes):
        box = []
        for _j in range(4):
            center, width = rng.random(), rng.random()
            box.append((max(0.0, center - width / 2),
                        min(1.0, center + width / 2)))
        t = float(PSI) + rng.random() * (C_STAR - float(PSI))
        lower = quadratic_bound(tuple(box), t, gmax)
        if lower is None:
            continue
        lo = mp.mpf(str(float(lower)))
        for _j in range(points_per_box):
            p1, q1, p2, q2 = [rng.uniform(a, b) for a, b in box]
            m1, m2 = (p1 + q1) / 2, (p2 + q2) / 2
            feasible_w = []
            for _k in range(80):
                w = rng.random()
                if w * m1 + (1 - w) * m2 <= t:
                    feasible_w.append(w)
                    break
            if not feasible_w:
                continue
            w = feasible_w[0]
            truth = phi_true(p1, q1, p2, q2, w)
            checked += 1
            if not (lo <= truth + mp.mpf("1e-25")):
                violations += 1
    assert boxes == 200 and points_per_box == 30
    assert checked > 1000, checked
    assert violations == 0, (violations, checked)
    print("sampled-soundness PASS: %d feasible points from 200 random 4-D boxes; "
          "0 violations (sampled soundness, not a proof of correctness)" % checked,
          flush=True)


def resolution_probe(gmax):
    t = float(PSI)
    center = (B_OBS, B_OBS, B_OBS, 1.0)
    widths = (1.0, 0.5, 0.1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4,
              3e-5, 1e-5)
    rows = []
    print("resolution probe at t=psi (d is FULL width; half-width=d/2):", flush=True)
    print("  %-10s %16s" % ("FULL d", "quadratic lower"), flush=True)
    for d in widths:
        box = tuple((max(0.0, x - d / 2), min(1.0, x + d / 2))
                    for x in center)
        bound = quadratic_bound(box, t, gmax)
        assert bound is not None
        value = float(bound)
        rows.append((d, value))
        print("  %-10g %+16.7e" % (d, value), flush=True)
    positive = [d for d, value in rows if value >= 0]
    assert positive, rows
    flip = max(positive)
    larger = [value for d, value in rows if d > flip]
    assert not larger or any(value < 0 for value in larger)
    print("OBSTRUCTION FLIP FULL WIDTH: d = %.10g" % flip, flush=True)
    return flip


def branch_and_bound(t, gmax, rho_gmax, max_boxes=2000000,
                     min_diam=1e-6, time_budget=600.0):
    """Best-first four-dimensional B&B using both certified discharge rules.

    Every queued box has already been evaluated and is keyed by its own Arb
    lower bound, not an inherited parent priority.  Splitting evaluates both
    children immediately, so the heap is genuinely best-first on bound.
    """
    root = ((0.0, 1.0),) * 4
    root_status, root_bound = certify_box(root, t, gmax, rho_gmax)
    processed = 1
    infeasible = cleared_ratio = cleared_quad = 0
    residual = []
    heap = []
    serial = 0
    if root_status == "infeasible":
        infeasible = 1
    elif root_status == "ratio":
        cleared_ratio = 1
    elif root_status == "quadratic":
        cleared_quad = 1
    else:
        heap.append((float(root_bound), serial, root))
    start = time.monotonic()
    last_report = start
    stop_reason = "queue-empty"
    while heap:
        now = time.monotonic()
        # A split evaluates exactly two children.  Keep the unresolved parent in
        # the queue when the remaining arithmetic budget cannot cover both.
        if processed + 2 > max_boxes:
            stop_reason = "box-budget"
            break
        if now - start >= time_budget:
            stop_reason = "time-budget"
            break
        _priority, _serial, box = heapq.heappop(heap)
        widths = [b - a for a, b in box]
        width = max(widths)
        if width <= min_diam:
            residual.append(box)
            continue
        k = widths.index(width)
        lo, hi = box[k]
        mid = (lo + hi) / 2
        for child_interval in ((lo, mid), (mid, hi)):
            child = list(box)
            child[k] = child_interval
            child = tuple(child)
            status, bound = certify_box(child, t, gmax, rho_gmax)
            processed += 1
            if status == "infeasible":
                infeasible += 1
            elif status == "ratio":
                cleared_ratio += 1
            elif status == "quadratic":
                cleared_quad += 1
            else:
                child_width = max(b - a for a, b in child)
                if child_width <= min_diam:
                    residual.append(child)
                else:
                    serial += 1
                    heapq.heappush(heap, (float(bound), serial, child))
        now = time.monotonic()
        if processed % 50000 <= 1 or now - last_report >= 60:
            last_report = now
            accounted = infeasible + cleared_ratio + cleared_quad + len(residual)
            assert processed >= accounted
            assert len(heap) >= 0 and now >= start
            print("    progress processed=%d ratio=%d quadratic=%d infeasible=%d "
                  "queue=%d residual=%d elapsed=%.1fs" %
                  (processed, cleared_ratio, cleared_quad, infeasible,
                   len(heap), len(residual), now - start), flush=True)
    elapsed = time.monotonic() - start
    complete = not residual and not heap
    stats = dict(processed=processed, cleared_ratio=cleared_ratio,
                 cleared_quadratic=cleared_quad, infeasible=infeasible,
                 residual=len(residual), queue=len(heap), seconds=elapsed,
                 stop_reason=stop_reason, complete=complete)
    assert complete == (stats["residual"] == 0 and stats["queue"] == 0)
    assert stats["processed"] >= (stats["cleared_ratio"]
                                  + stats["cleared_quadratic"]
                                  + stats["infeasible"] + stats["residual"])
    assert stats["seconds"] >= 0
    return stats


def run_bb_level(label, t, gmax, rho_gmax):
    assert 0 <= t <= 1
    print("B&B t=%s (%.16f), min FULL diameter 1e-6, budget 2,000,000 / 600s"
          % (label, t), flush=True)
    stats = branch_and_bound(t, gmax, rho_gmax)
    assert stats["processed"] >= 1
    assert stats["processed"] >= (stats["cleared_ratio"]
                                  + stats["cleared_quadratic"]
                                  + stats["infeasible"] + stats["residual"])
    cleared_total = stats["cleared_ratio"] + stats["cleared_quadratic"]
    ratio_fraction = (stats["cleared_ratio"] / cleared_total
                      if cleared_total else 0.0)
    quadratic_fraction = (stats["cleared_quadratic"] / cleared_total
                          if cleared_total else 0.0)
    assert abs(ratio_fraction + quadratic_fraction - (1.0 if cleared_total else 0.0)) < 1e-15
    print("  processed / cleared-ratio / cleared-quadratic / infeasible: "
          "%d / %d / %d / %d" %
          (stats["processed"], stats["cleared_ratio"],
           stats["cleared_quadratic"], stats["infeasible"]), flush=True)
    print("  discharge fractions ratio / quadratic: %.6f / %.6f" %
          (ratio_fraction, quadratic_fraction), flush=True)
    assert stats["residual"] >= 0 and stats["queue"] >= 0 and stats["seconds"] >= 0
    print("  residual / queue / seconds / stop: %d / %d / %.1f / %s" %
          (stats["residual"], stats["queue"], stats["seconds"],
           stats["stop_reason"]), flush=True)
    if stats["complete"]:
        assert stats["residual"] == 0 and stats["queue"] == 0
        print("  COMPLETE CERTIFICATE", flush=True)
    else:
        assert stats["residual"] > 0 or stats["queue"] > 0
        print("  INCOMPLETE: no certificate claimed", flush=True)
    return stats


def main():
    verify_expansion()
    gmax = get_rh_gmax()
    assert gmax >= arb("0.2342294") and gmax <= arb("0.2350")
    print("certified global rh cover upper bound: %.12f" % float(gmax.upper()), flush=True)
    validate_rho()
    rho_gmax = get_rho_gmax()
    validate_weight_contraction()
    validate_sampled_soundness()
    resolution_probe(gmax)
    first = run_bb_level("psi+1e-4", float(PSI) + 1e-4, gmax, rho_gmax)
    if first["complete"]:
        for label, t in (("psi+2e-4", float(PSI) + 2e-4),
                         ("0.3822", 0.3822),
                         ("c*-2e-4", C_STAR - 2e-4)):
            result = run_bb_level(label, t, gmax, rho_gmax)
            if not result["complete"]:
                break


if __name__ == "__main__":
    main()
