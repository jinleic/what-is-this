"""Are the far-from-obstruction residuals real minimisers or loose-bound artifacts?

Re-runs the corner B&B briefly, then for residual boxes at various distances
reports the TRUE Phi (mpmath) at the centre.  If true Phi is large where the
bound is negative, the bound is just loose there (an enclosure problem, possibly
fixable); if true Phi is small, those are genuine second minimisers.
"""

import heapq
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import diag_v2 as V
import cert as C
from entropy import PSI

T = float(PSI)
B0 = 0.32945473850303697239
a0 = (T - B0) / (1 - B0)
w0 = 1 - 2 * a0
OBS = np.array([B0, B0, B0, 1.0, w0])


def bb(lower, t, max_boxes=400_000, min_diam=1e-5, time_budget=120.0):
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
    return resid, dict(proc=processed, clr=cleared, infeas=infeas,
                       secs=time.time() - t0)


resid, st = bb(V.phi_corner, T)
print("processed %d cleared %d infeas %d residual %d (%.0fs)"
      % (st["proc"], st["clr"], st["infeas"], len(resid), st["secs"]))

cen = np.array([[(lo + hi) / 2 for lo, hi in b] for b in resid])
d_inf = np.max(np.abs(cen - OBS), axis=1)
true_phi = np.array([float(C.phi_true(*c)) for c in cen])
def _elo(c):
    box = tuple((max(0.0, c[j] - 1e-4), min(1.0, c[j] + 1e-4)) for j in range(5))
    e = V.phi_corner(box, T)
    return float("nan") if e is None else float(e.lower())
e_lo = np.array([_elo(c) for c in cen])

print("\nby distance-to-obstruction bucket: median true Phi, median bound gap")
print("%-22s %10s %14s %14s" % ("dist bucket", "count", "med truePhi", "med gap"))
edges = [0.0, 0.05, 0.1, 0.3, 0.6, 1.01]
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (d_inf >= lo) & (d_inf < hi)
    if m.sum():
        print("  [%4.2f,%4.2f)        %8d  %+13.3e  %+13.3e"
              % (lo, hi, m.sum(), np.median(true_phi[m]),
                 np.median((true_phi - e_lo)[m])))

print("\nsmallest-true-Phi residual centres (potential second minimisers):")
idx = np.argsort(true_phi)[:10]
for i in idx:
    feas = "FEAS" if (cen[i, 4] * (cen[i, 0] + cen[i, 1]) / 2
                      + (1 - cen[i, 4]) * (cen[i, 2] + cen[i, 3]) / 2) <= T + 1e-9 else "infeas"
    print("   (%.3f,%.3f,%.3f,%.3f,%.3f) truePhi=%+.2e d=%.3f %s"
          % (cen[i, 0], cen[i, 1], cen[i, 2], cen[i, 3], cen[i, 4],
             true_phi[i], d_inf[i], feas))
