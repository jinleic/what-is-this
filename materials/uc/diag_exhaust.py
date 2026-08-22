"""Depth-first corner-bound B&B with a sound mean contractor.

The contractor keeps every point satisfying M <= t.  Residual records are saved
as ten float64 columns: five centres followed by five FULL box widths.
"""

import functools
import math
import os
import struct
from fractions import Fraction
import sys
import tempfile
import time

import numpy as np
from flint import arb, ctx
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy import PSI


ctx.prec = 160
mp.dps = 50

ALPHA = arb("0.0356069")
ONE = arb(1)
TWO = arb(2)
HALF = arb(1) / 2
LOG2 = TWO.log()
ZERO = arb(0)
B_OBS_TEXT = "0.32945473850303697239"

MIN_DIAM = 2e-3
BOX_BUDGET = 8_000_000
TIME_BUDGET = 600.0
PROGRESS_EVERY = 500_000
OUTPUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "exhaust_residuals.npy")

T_MP = mp.mpf(PSI) + mp.mpf("0.0001")
T_TEXT = mp.nstr(T_MP, 45)
T_ARB = arb(T_TEXT)
# float(arf) is round-to-nearest.  One step toward +infinity therefore remains
# at least as permissive as the upper endpoint of the Arb enclosure for t.
T_UPPER_FLOAT = math.nextafter(float(T_ARB.upper()), math.inf)
T_SAMPLE_FLOAT = float(T_MP)


# Copied locally from cert.py: importing cert executes its long validation.
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
    """Enclosure of the range of h over [u,v]; u and v are Arb points."""
    u, v = arb(u), arb(v)
    e = hull(h_pt(u), h_pt(v))
    if v <= HALF or u >= HALF:
        return e
    return hull(e, ONE)


def amin(a, b):
    a, b = arb(a), arb(b)
    if a <= b:
        return a
    if b <= a:
        return b
    return hull(a, b)


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
    """s* is nondecreasing in each argument, so corners bound its range."""
    return hull(sstar_pt(pl, rl), sstar_pt(ph, rh))


def rh_encl(u, v, gmax=None):
    """Enclosure of rh(p) over [u,v]."""
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
    """Two-sided, NaN-free enclosure of sqrt(rh) over [u,v]."""
    hi = rh_encl(u, v, gmax).upper()
    # rh >= 0.  Never take sqrt of a ball that may straddle a negative value:
    # take sqrt only at the positive upper endpoint and hull it with zero.
    return hull(arb(0), hi.sqrt())


def _rh_global(n=20000):
    """Rigorous upper bound for max_p rh(p), by an Arb cover of [0,1]."""
    best = arb(0)
    for i in range(n):
        e = rh_encl(arb(i) / n, arb(i + 1) / n)
        u = e.upper()
        if not (u <= best):
            best = amax(best, u)
    return best.upper()


RH_GMAX = _rh_global()
assert RH_GMAX >= arb("0.2342294") and RH_GMAX <= arb("0.2350")


def _arb_float_point(x):
    """Construct an Arb point/enclosure from the round-trip decimal of a float."""
    return arb(repr(float(x)))


@functools.lru_cache(maxsize=None)
def _h_float_range(lo, hi):
    return h_encl(_arb_float_point(lo), _arb_float_point(hi))


@functools.lru_cache(maxsize=None)
def _sqrt_rh_float_range(lo, hi):
    return sqrt_rh_enc(_arb_float_point(lo), _arb_float_point(hi), RH_GMAX)


@functools.lru_cache(maxsize=131072)
def _sstar_h_float_range(pl, ph, ql, qh):
    p_lo, p_hi = _arb_float_point(pl), _arb_float_point(ph)
    q_lo, q_hi = _arb_float_point(ql), _arb_float_point(qh)
    s = sstar_encl(p_lo, p_hi, q_lo, q_hi)
    return h_encl(s.lower(), s.upper())


def phi_corner(box, t):
    """Feasible-restricted corner lower bound, or None if infeasible."""
    (p1lf, p1hf), (q1lf, q1hf), (p2lf, p2hf), (q2lf, q2hf), (wlf, whf) = box

    # Phi is symmetric in p_i,q_i.  This is the same p_i <= q_i reduction used
    # by diag_v2.py's corner bound.  On point endpoints, max is exact.
    q1lf, q2lf = max(q1lf, p1lf), max(q2lf, p2lf)
    if q1lf > q1hf or q2lf > q2hf:
        return None

    p1l, p1h = _arb_float_point(p1lf), _arb_float_point(p1hf)
    q1l, q1h = _arb_float_point(q1lf), _arb_float_point(q1hf)
    p2l, p2h = _arb_float_point(p2lf), _arb_float_point(p2hf)
    q2l, q2h = _arb_float_point(q2lf), _arb_float_point(q2hf)
    wl, wh = _arb_float_point(wlf), _arb_float_point(whf)

    w = hull(wl, wh)
    m1 = hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    if mean.lower() > t:
        return None
    m_star = amin(mean.upper(), t)
    coeff_star = TWO * (ONE - ALPHA) * (ONE - m_star) - ONE
    if not (coeff_star > ZERO):
        return arb(-1)

    L = (w * (_h_float_range(p1lf, p1hf)
              + _h_float_range(q1lf, q1hf)) / TWO
         + (ONE - w) * (_h_float_range(p2lf, p2hf)
                        + _h_float_range(q2lf, q2hf)) / TWO)
    L_lo = amax(L.lower(), ZERO)

    c1 = _sstar_h_float_range(p1lf, p1hf, q1lf, q1hf)
    c2 = _sstar_h_float_range(p2lf, p2hf, q2lf, q2hf)
    Cc = w * c1 + (ONE - w) * c2
    C_lo = amax(Cc.lower(), ZERO)

    r1 = (_sqrt_rh_float_range(p1lf, p1hf)
          + _sqrt_rh_float_range(q1lf, q1hf)) / TWO
    r2 = (_sqrt_rh_float_range(p2lf, p2hf)
          + _sqrt_rh_float_range(q2lf, q2hf)) / TWO
    sig = w * r1 + (ONE - w) * r2
    S_hi = sig.upper()
    return coeff_star * L_lo + ALPHA * C_lo - (ONE - ALPHA) * S_hi * S_hi


def _decimal_fraction(value):
    """Exact rational denoted by this endpoint's shortest decimal string."""
    return Fraction(repr(float(value)))


def _exact_fraction(value):
    """Exact rational for an Arb point or a shortest-decimal float endpoint."""
    if isinstance(value, arb):
        point = arb(value)
        if not point.is_exact():
            raise ValueError("contractor threshold must be an Arb point")
        mantissa, exponent = point.man_exp()
        mantissa, exponent = int(mantissa), int(exponent)
        return (Fraction(mantissa) * (2 ** exponent)
                if exponent >= 0
                else Fraction(mantissa, 2 ** (-exponent)))
    return _decimal_fraction(value)


def _float_decimal_down(value):
    """Largest practical float whose repr-rational is <= ``value``."""
    candidate = float(value)
    while _decimal_fraction(candidate) > value:
        candidate = math.nextafter(candidate, -math.inf)
    return candidate


def _float_decimal_up(value):
    """Smallest practical float whose repr-rational is >= ``value``."""
    candidate = float(value)
    while _decimal_fraction(candidate) < value:
        candidate = math.nextafter(candidate, math.inf)
    return candidate


def mean_contract(box, t_upper):
    """Contract w by an exact-rational necessary mean inequality.

    Every float endpoint used by the proof is interpreted as its shortest
    decimal string.  Convert those strings to :class:`Fraction` and perform the
    complete cutoff calculation exactly; binary float arithmetic padded by one
    ULP is not outward under this decimal semantics after cancellation.

    Let l_i be the exact mean of the two lower atom endpoints.  Every box point
    has m_i >= l_i and hence

        M >= l2 + w(l1-l2) = A(w).

    ``t_upper`` is itself a permissive decimal upper endpoint for t.  Solving
    A(w) <= t_upper exactly gives a necessary interval for every feasible w.
    Returned float cutoffs are then stepped until their DECIMAL meanings are
    outward from those exact rational cutoffs.
    """
    p1, q1, p2, q2 = box[:4]
    l1 = (_decimal_fraction(p1[0]) + _decimal_fraction(q1[0])) / 2
    l2 = (_decimal_fraction(p2[0]) + _decimal_fraction(q2[0])) / 2
    t_cap = _exact_fraction(t_upper)
    wl, wh = box[4]

    if l1 > t_cap and l2 > t_cap:
        return None

    new_wl, new_wh = wl, wh
    if l1 > t_cap and l2 <= t_cap:
        cutoff = (t_cap - l2) / (l1 - l2)
        new_wh = min(new_wh, _float_decimal_up(cutoff))
    elif l2 > t_cap and l1 <= t_cap:
        cutoff = (l2 - t_cap) / (l2 - l1)
        new_wl = max(new_wl, _float_decimal_down(cutoff))

    if new_wl > new_wh:
        return None
    if new_wl == wl and new_wh == wh:
        return box
    out = list(box)
    out[4] = (new_wl, new_wh)
    return tuple(out)


def _random_box(rng, case):
    low_pair = rng.uniform(0.02, 0.18, size=2)
    high_pair = rng.uniform(0.50, 0.72, size=2)
    if case == 0:
        lows = np.concatenate((low_pair, high_pair))
    elif case == 1:
        lows = np.concatenate((high_pair, low_pair))
    else:
        lows = rng.uniform(0.02, 0.18, size=4)
    highs = lows + (1.0 - lows) * rng.uniform(0.15, 0.95, size=4)
    wl = float(rng.uniform(0.0, 0.05))
    wh = float(rng.uniform(0.95, 1.0))
    return tuple((float(lows[i]), float(highs[i])) for i in range(4)) + ((wl, wh),)


def contractor_decimal_regression():
    """Exact checks for decimal semantics and near-cancelling cutoffs."""
    cases = (
        # Upper contraction: t-l2 is 1e-16 in decimal arithmetic.
        (((0.5000000000000002, 0.9), (0.5000000000000002, 0.9),
          (0.4999999999999999, 0.8), (0.4999999999999999, 0.8),
          (0.0, 1.0)), 0.5, "upper"),
        # Lower contraction: l2-t is 1e-16.
        (((0.4999999999999999, 0.8), (0.4999999999999999, 0.8),
          (0.5000000000000002, 0.9), (0.5000000000000002, 0.9),
          (0.0, 1.0)), 0.5, "lower"),
    )
    for box, t_cap, direction in cases:
        l1 = (_decimal_fraction(box[0][0])
              + _decimal_fraction(box[1][0])) / 2
        l2 = (_decimal_fraction(box[2][0])
              + _decimal_fraction(box[3][0])) / 2
        exact_t = _decimal_fraction(t_cap)
        contracted = mean_contract(box, t_cap)
        assert contracted is not None
        if direction == "upper":
            cutoff = (exact_t - l2) / (l1 - l2)
            assert _decimal_fraction(contracted[4][1]) >= cutoff
        else:
            cutoff = (l2 - exact_t) / (l2 - l1)
            assert _decimal_fraction(contracted[4][0]) <= cutoff
    # Direct conversion guards, including the endpoint used by the swap bug.
    values = (Fraction(1, 10**16), Fraction(1, 3),
              Fraction(9999999999999999, 10**16))
    for value in values:
        assert _decimal_fraction(_float_decimal_down(value)) <= value
        assert _decimal_fraction(_float_decimal_up(value)) >= value
    print("contractor decimal containment PASS (exact rationals): "
          "near-cancelling lower/upper cutoffs and 3 conversion guards",
          flush=True)


def contractor_nonremoval_test():
    rng = np.random.default_rng(20260812)
    tested_boxes = 100
    points_per_box = 50
    feasible_points = 0
    removals = 0
    lower_contracted = 0
    upper_contracted = 0
    unchanged = 0

    for i in range(tested_boxes):
        box = _random_box(rng, i % 3)
        contracted = mean_contract(box, T_UPPER_FLOAT)
        assert contracted is not None
        old_wl, old_wh = box[4]
        new_wl, new_wh = contracted[4]
        if new_wl > old_wl:
            lower_contracted += 1
        elif new_wh < old_wh:
            upper_contracted += 1
        else:
            unchanged += 1

        accepted = 0
        attempts = 0
        while accepted < points_per_box:
            attempts += 1
            assert attempts <= 200000
            coords = []
            for lo, hi in box[:4]:
                # Bias toward lower endpoints so all three contractor cases
                # yield feasible samples without changing their randomness.
                u = float(rng.random()) ** 4
                coords.append(lo + (hi - lo) * u)
            w = old_wl + (old_wh - old_wl) * float(rng.random())
            m1 = (coords[0] + coords[1]) / 2.0
            m2 = (coords[2] + coords[3]) / 2.0
            mean = w * m1 + (1.0 - w) * m2
            if mean > T_SAMPLE_FLOAT:
                continue
            coords_mp = [mp.mpf(repr(x)) for x in coords]
            w_mp = mp.mpf(repr(w))
            m1_mp = (coords_mp[0] + coords_mp[1]) / 2
            m2_mp = (coords_mp[2] + coords_mp[3]) / 2
            mean_mp = w_mp * m1_mp + (1 - w_mp) * m2_mp
            if mean_mp > T_MP:
                continue
            assert mean <= T_SAMPLE_FLOAT and mean_mp <= T_MP
            feasible_points += 1
            accepted += 1
            if not (new_wl <= w <= new_wh):
                removals += 1

    assert tested_boxes == 100
    assert points_per_box == 50
    assert feasible_points == tested_boxes * points_per_box == 5000
    assert removals == 0
    assert lower_contracted == 34 and upper_contracted == 33 and unchanged == 33
    print("contractor non-removal PASS on 100 sampled boxes x 50 sampled feasible "
          "points: 0 removals; lower/upper/unchanged boxes = 34/33/33",
          flush=True)


def _assert_accounting(processed, cleared, infeasible, residual, branched, stack_len):
    assert processed == cleared + infeasible + residual + branched
    assert stack_len == 1 + branched - cleared - infeasible - residual
    assert min(processed, cleared, infeasible, residual, branched, stack_len) >= 0


def run_bb():
    root = ((0.0, 1.0),) * 5
    stack = [root]
    processed = 0
    cleared = 0
    infeasible = 0
    contractor_infeasible = 0
    corner_infeasible = 0
    contractor_lower = 0
    contractor_upper = 0
    residual = 0
    branched = 0
    next_progress = PROGRESS_EVERY
    residual_file = tempfile.TemporaryFile(mode="w+b")
    started = time.monotonic()
    deadline = started + TIME_BUDGET
    stop_reason = "empty"

    while stack:
        if processed >= BOX_BUDGET:
            stop_reason = "box-budget"
            break
        if time.monotonic() >= deadline:
            stop_reason = "time-budget"
            break

        box = stack.pop()
        processed += 1
        contracted = mean_contract(box, T_UPPER_FLOAT)
        if contracted is None:
            infeasible += 1
            contractor_infeasible += 1
        else:
            if contracted[4][0] > box[4][0]:
                contractor_lower += 1
            if contracted[4][1] < box[4][1]:
                contractor_upper += 1
            bound = phi_corner(contracted, T_ARB)
            if bound is None:
                infeasible += 1
                corner_infeasible += 1
            elif bound >= ZERO:
                cleared += 1
            else:
                widths = tuple(hi - lo for lo, hi in contracted)
                max_width = max(widths)
                if max_width <= MIN_DIAM:
                    centres = tuple((lo + hi) / 2.0 for lo, hi in contracted)
                    assert max(widths) <= MIN_DIAM
                    residual_file.write(struct.pack("=10d", *(centres + widths)))
                    residual += 1
                else:
                    k = max(range(5), key=lambda j: widths[j])
                    lo, hi = contracted[k]
                    mid = (lo + hi) / 2.0
                    assert lo < mid < hi
                    low_child = list(contracted)
                    high_child = list(contracted)
                    low_child[k] = (lo, mid)
                    high_child[k] = (mid, hi)
                    # Push high first so the low child is visited next: genuine
                    # depth-first traversal with O(depth) live-tree memory.
                    stack.append(tuple(high_child))
                    stack.append(tuple(low_child))
                    branched += 1

        if processed == next_progress:
            elapsed = time.monotonic() - started
            _assert_accounting(processed, cleared, infeasible, residual,
                               branched, len(stack))
            assert processed % PROGRESS_EVERY == 0 and elapsed >= 0.0
            assert contractor_lower + contractor_upper <= processed
            print("[running] processed=%d cleared=%d infeasible=%d "
                  "residual-at-min=%d stack=%d contracted-w=%d elapsed=%.1fs"
                  % (processed, cleared, infeasible, residual, len(stack),
                     contractor_lower + contractor_upper, elapsed),
                  flush=True)
            next_progress += PROGRESS_EVERY

    elapsed = time.monotonic() - started
    if not stack:
        stop_reason = "empty"
    _assert_accounting(processed, cleared, infeasible, residual, branched, len(stack))
    assert infeasible == contractor_infeasible + corner_infeasible
    assert stop_reason in ("empty", "box-budget", "time-budget")
    assert (stop_reason == "empty") == (len(stack) == 0)
    if stop_reason == "box-budget":
        assert processed == BOX_BUDGET
    if stop_reason == "time-budget":
        assert elapsed >= TIME_BUDGET

    residual_file.flush()
    residual_file.seek(0)
    flat = np.fromfile(residual_file, dtype=np.float64, count=residual * 10)
    residual_file.close()
    assert flat.size == residual * 10
    rows = flat.reshape((residual, 10))
    if residual:
        assert np.all(rows[:, 5:] >= 0.0)
        assert np.all(rows[:, 5:] <= MIN_DIAM)
        assert np.all(rows[:, :5] >= 0.0) and np.all(rows[:, :5] <= 1.0)

    return rows, {
        "processed": processed,
        "cleared": cleared,
        "infeasible": infeasible,
        "contractor_infeasible": contractor_infeasible,
        "corner_infeasible": corner_infeasible,
        "contractor_lower": contractor_lower,
        "contractor_upper": contractor_upper,
        "residual": residual,
        "branched": branched,
        "stack": len(stack),
        "elapsed": elapsed,
        "reason": stop_reason,
    }


def _quantiles(values):
    probabilities = np.array([0.0, 0.10, 0.25, 0.50, 0.75, 0.90, 1.0])
    result = np.quantile(values, probabilities)
    assert result.shape == probabilities.shape
    assert np.all(np.diff(result) >= 0.0)
    assert result[0] == np.min(values) and result[-1] == np.max(values)
    return result


def save_and_report(rows, stats):
    complete = stats["reason"] == "empty"
    label = "complete" if complete else "partial"
    assert complete == (stats["stack"] == 0)
    assert stats["processed"] == (stats["cleared"] + stats["infeasible"]
                                  + stats["residual"] + stats["branched"])
    assert stats["residual"] == rows.shape[0]
    assert stats["elapsed"] >= 0.0
    assert stats["contractor_lower"] + stats["contractor_upper"] <= stats["processed"]

    print("[%s] stop=%s; processed=%d branched=%d elapsed=%.1fs"
          % (label, stats["reason"], stats["processed"], stats["branched"],
             stats["elapsed"]), flush=True)
    print("[%s] four-way tally: cleared=%d; infeasible=%d; "
          "residual-at-min-diam=%d; unprocessed-stack-remainder=%d"
          % (label, stats["cleared"], stats["infeasible"], stats["residual"],
             stats["stack"]), flush=True)
    exhausted_text = "YES" if complete else "NO"
    assert (exhausted_text == "YES") == complete
    print("[%s] stack exhausted=%s; geometry scope=%s"
          % (label, exhausted_text,
             "all refined residuals" if complete
             else "refined residual subset only (stack remainder excluded)"),
          flush=True)
    assert stats["infeasible"] == (stats["contractor_infeasible"]
                                    + stats["corner_infeasible"])
    print("[%s] infeasible split: mean-contractor=%d; corner-bound=%d"
          % (label, stats["contractor_infeasible"],
             stats["corner_infeasible"]), flush=True)
    print("[%s] mean-contractor w contractions: lower=%d; upper=%d; total=%d"
          % (label, stats["contractor_lower"], stats["contractor_upper"],
             stats["contractor_lower"] + stats["contractor_upper"]), flush=True)

    np.save(OUTPUT_PATH, rows)
    loaded = np.load(OUTPUT_PATH, mmap_mode="r")
    assert loaded.shape == rows.shape == (stats["residual"], 10)
    if rows.size:
        assert np.array_equal(np.asarray(loaded), rows)
    assert MIN_DIAM == 2e-3
    print("[%s] saved %d refined residual records to %s; columns are five "
          "centres then five FULL widths, each max width <= %.6g"
          % (label, rows.shape[0], OUTPUT_PATH, MIN_DIAM), flush=True)

    if rows.shape[0] == 0:
        assert stats["residual"] == 0
        print("[%s] refined residual subset is empty; obstruction-local, "
              "sink-boundary, and other counts are 0/0/0" % label, flush=True)
        return

    centres = rows[:, :5]
    b_obs = mp.mpf(B_OBS_TEXT)
    a_obs = (T_MP - b_obs) / (1 - b_obs)
    w_obs = 1 - 2 * a_obs
    obstruction = np.array([float(b_obs), float(b_obs), float(b_obs), 1.0,
                            float(w_obs)])
    d_obstruction = np.max(np.abs(centres - obstruction), axis=1)
    d_sink = np.min(np.minimum(centres[:, :4], 1.0 - centres[:, :4]), axis=1)
    assert np.all(d_obstruction >= 0.0)
    assert np.all(d_sink >= 0.0) and np.all(d_sink <= 0.5)

    q_obs = _quantiles(d_obstruction)
    q_sink = _quantiles(d_sink)
    assert len(q_obs) == len(q_sink) == 7
    print("[%s] L_inf centre-distance to obstruction quantiles "
          "min/10%%/25%%/50%%/75%%/90%%/max = %s"
          % (label, np.array2string(q_obs, precision=6, separator=",")),
          flush=True)
    print("[%s] centre-distance to nearest p/q sink boundary quantiles "
          "min/10%%/25%%/50%%/75%%/90%%/max = %s"
          % (label, np.array2string(q_sink, precision=6, separator=",")),
          flush=True)

    obstruction_local = d_obstruction <= MIN_DIAM
    sink_local_raw = d_sink <= MIN_DIAM
    obstruction_count = int(np.count_nonzero(obstruction_local))
    overlap_count = int(np.count_nonzero(obstruction_local & sink_local_raw))
    sink_only = int(np.count_nonzero((~obstruction_local) & sink_local_raw))
    other = int(rows.shape[0] - obstruction_count - sink_only)
    raw_sink_count = int(np.count_nonzero(sink_local_raw))
    assert obstruction_count + sink_only + other == rows.shape[0]
    assert raw_sink_count == sink_only + overlap_count
    assert 0 <= overlap_count <= obstruction_count

    denom = float(rows.shape[0])
    obs_fraction = obstruction_count / denom
    sink_raw_fraction = raw_sink_count / denom
    sink_only_fraction = sink_only / denom
    other_fraction = other / denom
    assert abs((obstruction_count + sink_only + other) / denom - 1.0) < 1e-15
    assert abs(raw_sink_count / denom - sink_raw_fraction) < 1e-15
    print("[%s] exclusive geometry at threshold %.6g: obstruction-local=%d "
          "(%.6f); sink-boundary-only=%d (%.6f); other=%d (%.6f)"
          % (label, MIN_DIAM, obstruction_count, obs_fraction, sink_only,
             sink_only_fraction, other, other_fraction), flush=True)
    print("[%s] raw sink-boundary fraction (including obstruction overlap) "
          "= %d/%d = %.6f; obstruction/sink overlap=%d"
          % (label, raw_sink_count, rows.shape[0], sink_raw_fraction,
             overlap_count), flush=True)


def main():
    total_started = time.monotonic()
    assert T_UPPER_FLOAT >= T_SAMPLE_FLOAT
    assert T_MP == mp.mpf(PSI) + mp.mpf("0.0001")
    assert T_TEXT == mp.nstr(T_MP, 45)
    assert MIN_DIAM == 0.002 and BOX_BUDGET == 8_000_000
    assert TIME_BUDGET == 600.0 and PROGRESS_EVERY == 500_000
    print("B&B configuration: t=psi+1e-4=%s; min_diam=0.002 FULL width; "
          "budgets=8000000 boxes or 600s" % T_TEXT, flush=True)
    contractor_decimal_regression()
    contractor_nonremoval_test()
    rows, stats = run_bb()
    save_and_report(rows, stats)
    total_elapsed = time.monotonic() - total_started
    assert total_elapsed <= 720.0
    print("asserted total script wall time %.1fs <= 720s" % total_elapsed,
          flush=True)


if __name__ == "__main__":
    main()
