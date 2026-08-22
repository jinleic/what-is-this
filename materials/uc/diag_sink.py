"""Which coordinate's width actually blocks the sink-collar leaves?

Phase 1: run cert3.certify (full v3 rule chain) on the p1-collar slab of
slice 0 -- root (0,1/16) x (0,1)^3 x (1/2, 9/16) -- dumping residual leaves.
Phase 2: for a sample of dumped leaves and each coordinate j, refine ONLY
coordinate j to 1/8 of its width and test whether every child then clears
via any rule.  This measures per-coordinate responsibility and decides
between anisotropic refinement, a new analytic rule, or full refinement.
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import ctx
from mpmath import mp

import cert2
import cert3
from arbcore import get_rh_gmax
from entropy import PSI

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, "sinkdiag_resid.npy")
T = mp.mpf(PSI) + mp.mpf("0.0001")
ROOT = ((0.0, 0.0625), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.5, 0.5625))


def rules_clear(box, t_arb, gmax, rho_gmax, lambdas, center_lams):
    """True iff any sound rule certifies the box (mirrors certify's chain)."""
    import diag_exhaust
    import math
    t_upper = math.nextafter(float(t_arb.upper()), math.inf)
    cb = diag_exhaust.mean_contract(box, t_upper)
    if cb is None:
        return True                      # certified infeasible
    box = cb
    corner = diag_exhaust.phi_corner(box, t_arb)
    if corner is None:
        return True
    if corner >= 0:
        return True
    if cert2.ratio_rule(box[:4], t_arb, rho_gmax=rho_gmax):
        return True
    val = cert3.centered5_best(box, t_arb, center_lams, gmax)
    if val is not None and val >= 0:
        return True
    if val is None:
        vm = cert3.centered5_mixed(box, t_arb, center_lams, gmax)
        if vm is not None and vm >= 0:
            return True
    if box[3][0] >= cert3.Q2PIN:
        pinned = [lam for lam in lambdas if cert3.q2_pin(box, t_arb, lam, gmax)]
        if pinned:
            face_box = (box[0], box[1], box[2], box[4])
            done, _ = cert3.face_bb(face_box, t_arb, tuple(pinned), gmax,
                                    1e-3, 100000, None)
            if done:
                return True
    return False


if __name__ == "__main__":
    t0 = time.time()
    gmax = get_rh_gmax()
    rho_gmax = cert2.get_rho_gmax()
    ctx.prec = 80
    lambdas, lamhat = cert3.lambda_family(T)
    center_lams = (0.0, float(lamhat))
    print("phase 1: certify the p1-collar slab, dumping residuals", flush=True)
    ok, stats = cert3.certify(
        T, ROOT, max_boxes=30_000_000, time_budget=2400.0,
        min_width=1e-3, face_min_width=1e-3, gmax=gmax, rho_gmax=rho_gmax,
        lambdas=lambdas, center_lams=center_lams,
        progress_every=2_000_000, residual_dump=DUMP,
    )
    print("phase 1: complete=%s stats=%s (%.0fs)"
          % (ok, dict(stats), time.time() - t0), flush=True)
    if not os.path.exists(DUMP):
        print("no residuals dumped; nothing to diagnose", flush=True)
        sys.exit(0)
    rows = np.load(DUMP)
    print("phase 2: %d residual leaves; sampling for per-coordinate "
          "responsibility" % len(rows), flush=True)
    rng = np.random.default_rng(7)
    sample = rows[rng.choice(len(rows), size=min(200, len(rows)),
                             replace=False)]
    t_arb = cert3._arb(mp.nstr(T, 45))
    names = ("p1", "q1", "p2", "q2", "w")
    solo = np.zeros(5, dtype=int)
    tested = 0
    unresolved = []
    for row in sample:
        cen, wid = row[:5], row[5:]
        leaf = tuple((cen[j] - wid[j] / 2, cen[j] + wid[j] / 2)
                     for j in range(5))
        tested += 1
        resolved_by = []
        for j in range(5):
            lo, hi = leaf[j]
            parts = np.linspace(lo, hi, 9)
            all_clear = True
            for k in range(8):
                child = list(leaf)
                child[j] = (float(parts[k]), float(parts[k + 1]))
                if not rules_clear(tuple(child), t_arb, gmax, rho_gmax,
                                   lambdas, center_lams):
                    all_clear = False
                    break
            if all_clear:
                solo[j] += 1
                resolved_by.append(j)
        if not resolved_by:
            unresolved.append(leaf)
        if tested % 25 == 0:
            print("  %d leaves tested; solo-clear counts %s; unresolved %d"
                  % (tested, dict(zip(names, solo.tolist())),
                     len(unresolved)), flush=True)
    print("phase 2 result over %d sampled leaves:" % tested, flush=True)
    for j in range(5):
        print("  refine ONLY %-2s (x8): clears %5d / %d  (%.2f)"
              % (names[j], solo[j], tested, solo[j] / max(tested, 1)),
              flush=True)
    print("  no single coordinate suffices: %d / %d"
          % (len(unresolved), tested), flush=True)
    for leaf in unresolved[:6]:
        print("   unresolved leaf:", ["(%.6f,%.6f)" % b for b in leaf],
              flush=True)
    print("total %.0fs" % (time.time() - t0), flush=True)
