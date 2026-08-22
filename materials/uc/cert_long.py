"""Two-phase decisive experiment: corner-exhaust then ratio-discharge the leaves.

Phase A: depth-first corner-bound B&B at t = psi+1e-4 over [0,1]^5, min FULL
width 2e-3, generous budget, collecting every residual leaf.
Phase B: apply cert2's certified Cauchy-Schwarz ratio rule to each leaf.
Phase C: classify the survivors geometrically: sink-adjacent / obstruction-
local (both orbit-swap images) / mean-boundary / other.

Sound rules only; tallies are exact counts, geometry is exact per-leaf
arithmetic.  No Phi >= 0 claim is made unless the tree is exhausted AND no
survivor remains (assert-gated).
"""

import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import diag_exhaust as DX          # float-fast corner bound + contractor
import cert2 as C2                 # certified ratio rule + rho caps
from flint import arb
from entropy import PSI
from mpmath import mp

T_MP = mp.mpf(PSI) + mp.mpf("0.0001")
T_ARB = arb(mp.nstr(T_MP, 45))
T_UP = math.nextafter(float(T_ARB.upper()), math.inf)

B0 = 0.32945473850303697239
A0 = (float(T_MP) - B0) / (1.0 - B0)
W0 = 1.0 - 2.0 * A0
OBS_A = np.array([B0, B0, B0, 1.0, W0])       # (orbit1=(b,b), orbit2=(b,1), w)
OBS_B = np.array([B0, 1.0, B0, B0, 1.0 - W0])  # orbit-swapped image

MIN_W = 2e-3          # FULL width floor per coordinate
BOX_BUDGET = 200_000_000
TIME_BUDGET = 4500.0

# Geometric pre-filter for the inline ratio rule: rho <= kappa/(1-alpha)
# ~ 0.19894 roughly requires mass in [0,~0.002] U [~0.48,1]; call the (arb,
# ~1ms) rule only when every atom interval lies in the fire zone's outer
# approximation, so the bulk pays nothing.
RATIO_LO_MAX = 0.05
RATIO_HI_MIN = 0.42


def ratio_plausible(box):
    for lo, hi in box[:4]:
        if not (hi <= RATIO_LO_MAX or lo >= RATIO_HI_MIN):
            return False
    return True


def phase_a(rho_gmax):
    # Phi is invariant under the orbit swap (p1,q1,w) <-> (p2,q2,1-w): every
    # aggregate is a w-convex combination of the two symmetric orbit terms.
    # Certifying w in [0, 1/2] therefore certifies the whole cube.
    root = ((0.0, 1.0),) * 4 + ((0.0, 0.5),)
    stack = [root]
    leaves = []
    processed = cleared = infeasible = ratio_cleared = 0
    t0 = time.time()
    while stack:
        if processed >= BOX_BUDGET or time.time() - t0 > TIME_BUDGET:
            break
        box = stack.pop()
        processed += 1
        if processed % 2_000_000 == 0:
            print("  A: %d processed, %d corner, %d ratio, %d infeasible, "
                  "%d leaves, stack %d (%.0fs)"
                  % (processed, cleared, ratio_cleared, infeasible,
                     len(leaves), len(stack), time.time() - t0), flush=True)
        cb = DX.mean_contract(box, T_UP)
        if cb is None:
            infeasible += 1
            continue
        box = cb
        e = DX.phi_corner(box, T_ARB)
        if e is None:
            infeasible += 1
            continue
        if e >= 0:
            cleared += 1
            continue
        if ratio_plausible(box):
            W = DX.hull(arb(repr(box[4][0])), arb(repr(box[4][1])))
            if C2.ratio_rule(box[:4], float(T_MP), W=W, rho_gmax=rho_gmax):
                ratio_cleared += 1
                continue
        widths = [hi - lo for lo, hi in box]
        k = int(np.argmax(widths))
        if widths[k] <= MIN_W:
            leaves.append(box)
            continue
        lo, hi = box[k]
        mid = (lo + hi) / 2.0
        for part in ((lo, mid), (mid, hi)):
            nb = list(box)
            nb[k] = part
            stack.append(tuple(nb))
    exhausted = not stack and processed < BOX_BUDGET
    return leaves, dict(processed=processed, cleared=cleared,
                        ratio=ratio_cleared, infeasible=infeasible,
                        stack=len(stack), exhausted=exhausted,
                        secs=time.time() - t0)


def phase_b(leaves):
    """Certified ratio rule on each leaf; returns the surviving leaves."""
    rho_gmax = C2.get_rho_gmax()
    print("  B: certified global rho cap = %.9f" % float(arb(rho_gmax)),
          flush=True)
    cache = {}

    def cached_fire(leaf):
        atoms = leaf[:4]
        wl, wh = leaf[4]
        key = (atoms, round(wl, 12), round(wh, 12))
        r = cache.get(key)
        if r is None:
            W = DX.hull(arb(repr(wl)), arb(repr(wh)))
            r = C2.ratio_rule(atoms, float(T_MP), W=W, rho_gmax=rho_gmax)
            cache[key] = bool(r)
        return r

    survivors = []
    t0 = time.time()
    fired = 0
    for i, leaf in enumerate(leaves):
        if i % 200_000 == 0 and i:
            print("  B: %d/%d leaves, %d ratio-cleared (%.0fs)"
                  % (i, len(leaves), fired, time.time() - t0), flush=True)
        if cached_fire(leaf):
            fired += 1
        else:
            survivors.append(leaf)
    return survivors, fired, time.time() - t0


def phase_c(survivors):
    if not survivors:
        print("C: no survivors to classify")
        return
    cen = np.array([[(lo + hi) / 2 for lo, hi in b] for b in survivors])
    wid = np.array([[hi - lo for lo, hi in b] for b in survivors])
    edge = 2 * MIN_W
    sink = np.zeros(len(survivors), dtype=bool)
    for j in range(4):                          # atom coordinates only
        sink |= (cen[:, j] - wid[:, j] / 2 <= edge)
        sink |= (cen[:, j] + wid[:, j] / 2 >= 1 - edge)
    d_obs = np.minimum(np.max(np.abs(cen - OBS_A), axis=1),
                       np.max(np.abs(cen - OBS_B), axis=1))
    obs = d_obs <= 0.05
    # mean interval straddles t?
    m1l = (cen[:, 0] - wid[:, 0] / 2 + cen[:, 1] - wid[:, 1] / 2) / 2
    m1h = (cen[:, 0] + wid[:, 0] / 2 + cen[:, 1] + wid[:, 1] / 2) / 2
    m2l = (cen[:, 2] - wid[:, 2] / 2 + cen[:, 3] - wid[:, 3] / 2) / 2
    m2h = (cen[:, 2] + wid[:, 2] / 2 + cen[:, 3] + wid[:, 3] / 2) / 2
    wl = cen[:, 4] - wid[:, 4] / 2
    wh = cen[:, 4] + wid[:, 4] / 2
    # outer mean range over the box: coefficients of m1, m2 are nonnegative,
    # so extremes use (m1l, m2l) or (m1h, m2h) at a w endpoint (float,
    # classification only -- certification never uses these numbers)
    Ml = np.minimum(wl * m1l + (1 - wl) * m2l, wh * m1l + (1 - wh) * m2l)
    Mh = np.maximum(wl * m1h + (1 - wl) * m2h, wh * m1h + (1 - wh) * m2h)
    tf = float(T_MP)
    meanb = (Ml <= tf) & (Mh >= tf)
    only_sink = sink & ~obs
    only_obs = obs & ~sink
    both = sink & obs
    other = ~sink & ~obs
    print("C: survivor classification (n = %d):" % len(survivors))
    print("   sink-adjacent only : %8d  (%.3f)"
          % (only_sink.sum(), only_sink.mean()))
    print("   obstruction-local  : %8d  (%.3f)  [Chebyshev <= 0.05 to either image]"
          % (only_obs.sum(), only_obs.mean()))
    print("   both               : %8d" % both.sum())
    print("   other              : %8d  (%.3f)" % (other.sum(), other.mean()))
    print("   mean-boundary flag : %8d  (%.3f of survivors)"
          % (meanb.sum(), meanb.mean()))
    if other.sum():
        idx = np.where(other)[0][:8]
        print("   sample 'other' centres:")
        for i in idx:
            print("     (%.4f,%.4f,%.4f,%.4f,%.4f)  d_obs=%.3f"
                  % (*cen[i], d_obs[i]))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "certlong_survivors.npy")
    np.save(out, np.hstack([cen, wid]))
    print("   saved centres+widths ->", out)


if __name__ == "__main__":
    print("cert_long v2: t = psi+1e-4 = %s" % mp.nstr(T_MP, 30), flush=True)
    print("Phase A: corner+ratio DFS on w<=1/2 (orbit-swap symmetry), min FULL"
          " width %g, budget %d boxes / %.0fs"
          % (MIN_W, BOX_BUDGET, TIME_BUDGET), flush=True)
    RG = C2.get_rho_gmax()
    print("certified global rho cap = %.9f" % float(arb(RG)), flush=True)
    leaves, st = phase_a(RG)
    print("A: processed %(processed)d corner %(cleared)d ratio %(ratio)d "
          "infeasible %(infeasible)d stack %(stack)d (%(secs).0fs)" % st,
          flush=True)
    print("A: exhausted = %s, residual leaves = %d"
          % (st["exhausted"], len(leaves)), flush=True)
    survivors, fired, tb = phase_b(leaves)
    print("B: ratio rule cleared %d / %d leaves (%.0fs); survivors %d"
          % (fired, len(leaves), tb, len(survivors)), flush=True)
    phase_c(survivors)
    if st["exhausted"] and not survivors:
        print("COMPLETE CERTIFICATE at t = psi+1e-4 (corner + ratio)")
    else:
        print("NO certificate claimed (exhausted=%s, survivors=%d)"
              % (st["exhausted"], len(survivors)))
