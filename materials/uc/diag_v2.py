"""Head-to-head: current phi_encl (v1) vs the feasible-restricted corner bound (v2).

WHY v2.  The t=psi residual diagnostic showed the B&B clears NOTHING: every
processed box has Phi_lo < 0 and gets subdivided, because phi_encl evaluates
Phi over the whole box including the high-mean infeasible part.

THE BOUND (sound; uses inf(f+g+h) >= inf f + inf g + inf h, NO independence).
On the feasible region M <= t, coeff(M)L=(2(1-a)(1-M)-1)L has d/dM=-2(1-a)L<=0,
so it is minimised at M* = min(m_hi, t); coeff* = 2(1-a)(1-M*)-1 > 0 for t<=c*.
  => coeff*L >= coeff* L_lo ;  aC >= a C_lo ;  -(1-a)S^2 >= -(1-a)S_hi^2.
So inf_{box cap M<=t} Phi >= coeff* L_lo + a C_lo - (1-a) S_hi^2, coeff* a point.

This module is import-safe: the benchmark runs only under __main__.
"""

import heapq
import os
import sys
import time

import numpy as np
from flint import arb
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert as C
from entropy import PSI

mp.dps = 40
T = float(PSI)

ONE = arb(1)
TWO = arb(2)
ALPHA = C.ALPHA
RH_GMAX = C.RH_GMAX


def phi_corner(box, t):
    """Feasible-restricted corner lower bound, or None if certified infeasible."""
    (p1l, p1h), (q1l, q1h), (p2l, p2h), (q2l, q2h), (wl, wh) = box
    p1l, p1h = arb(p1l), arb(p1h)
    q1l, q1h = arb(q1l), arb(q1h)
    p2l, p2h = arb(p2l), arb(p2h)
    q2l, q2h = arb(q2l), arb(q2h)
    wl, wh = arb(wl), arb(wh)
    q1l, q2l = C.amax(q1l, p1l), C.amax(q2l, p2l)
    if q1l > q1h or q2l > q2h:
        return None
    w = C.hull(wl, wh)
    m1 = C.hull((p1l + q1l) / TWO, (p1h + q1h) / TWO)
    m2 = C.hull((p2l + q2l) / TWO, (p2h + q2h) / TWO)
    mean = w * m1 + (ONE - w) * m2
    if mean.lower() > arb(t):
        return None
    m_star = C.amin(mean.upper(), arb(t))
    coeff_star = TWO * (ONE - ALPHA) * (ONE - m_star) - ONE
    if not (coeff_star > arb(0)):
        return arb(-1)
    L = w * (C.h_encl(p1l, p1h) + C.h_encl(q1l, q1h)) / TWO \
        + (ONE - w) * (C.h_encl(p2l, p2h) + C.h_encl(q2l, q2h)) / TWO
    L_lo = C.amax(L.lower(), arb(0))
    s1 = C.sstar_encl(p1l, p1h, q1l, q1h)
    s2 = C.sstar_encl(p2l, p2h, q2l, q2h)
    Cc = w * C.h_encl(s1.lower(), s1.upper()) \
        + (ONE - w) * C.h_encl(s2.lower(), s2.upper())
    C_lo = C.amax(Cc.lower(), arb(0))
    r1 = (C.sqrt_rh_enc(p1l, p1h, RH_GMAX)
          + C.sqrt_rh_enc(q1l, q1h, RH_GMAX)) / TWO
    r2 = (C.sqrt_rh_enc(p2l, p2h, RH_GMAX)
          + C.sqrt_rh_enc(q2l, q2h, RH_GMAX)) / TWO
    sig = w * r1 + (ONE - w) * r2
    S_hi = sig.upper()
    return coeff_star * L_lo + ALPHA * C_lo - (ONE - ALPHA) * S_hi * S_hi


def bb(lower, t, max_boxes=120000, min_diam=1e-6, time_budget=200.0):
    root = ((0.0, 1.0),) * 5
    heap = [(0.0, 0, root)]
    counter = processed = cleared = infeas = 0
    resid = []
    t0 = time.time()
    while heap:
        if processed >= max_boxes or time.time() - t0 > time_budget:
            resid.extend(b for _, _, b in heap)
            break
        _, _, box = heapq.heappop(heap)
        processed += 1
        e = lower(box, t)
        if e is None:
            infeas += 1
            continue
        if e >= 0:
            cleared += 1
            continue
        widths = [hi - lo for lo, hi in box]
        if max(widths) <= min_diam:
            resid.append(box)
            continue
        k = int(np.argmax(widths))
        lo, hi = box[k]
        mid = (lo + hi) / 2
        for sub in ((lo, mid), (mid, hi)):
            nb = list(box)
            nb[k] = sub
            heapq.heappush(heap, (float(e.lower()), counter, tuple(nb)))
            counter += 1
    return len(resid), dict(proc=processed, clr=cleared, infeas=infeas,
                            q=len(heap), secs=time.time() - t0)


if __name__ == "__main__":
    print("0. soundness of the corner bound (must be <= mpmath truth)")
    t_test = T + 1e-4
    bad = 0
    for (p1, q1, p2, q2, w) in [(0.3295, 0.3295, 0.3295, 1.0, 0.84),
                                (0.2, 0.3, 0.4, 0.9, 0.6),
                                (0.1, 0.2, 0.3, 0.5, 0.5),
                                (0.05, 0.05, 0.0, 1.0, 0.3)]:
        e = phi_corner(((p1, p1), (q1, q1), (p2, p2), (q2, q2), (w, w)), t_test)
        tv = float(C.phi_true(p1, q1, p2, q2, w))
        if e is None:
            continue
        lo = float(e.lower())
        ok = lo <= tv + 1e-9
        print("   pt (%.3f,%.3f,%.3f,%.3f,%.3f)  corner=%.6e  true=%.6e  %s"
              % (p1, q1, p2, q2, w, lo, tv, "OK" if ok else "FAIL"))
        if not ok:
            bad += 1
    print("   soundness failures:", bad)
    assert bad == 0

    print("\n1. head-to-head B&B at t = psi  (v1 = current phi_encl, v2 = corner)")
    for name, fn in (("v1 current", C.phi_encl), ("v2 corner  ", phi_corner)):
        nres, st = bb(fn, T, max_boxes=120000, time_budget=180.0)
        print("   %s  cleared %6d / infeas %6d / residual %6d  (proc %d, %.1fs)"
              % (name, st["clr"], st["infeas"], nres, st["proc"], st["secs"]))
