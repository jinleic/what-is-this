"""Does the corner-bound B&B localise residuals to the obstruction?

If yes -> hybrid certificate is viable (clear complement coarsely, handle the
obstruction neighbourhood specially).  If residuals stay diffuse -> dead end.

Runs the SOUND corner bound (diag_v2.phi_corner) at t = psi with a large budget
and reports the residual distribution by distance to the obstruction (b,b,b,1,w*).
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


def bb(lower, t, max_boxes=2_000_000, min_diam=1e-5, time_budget=600.0):
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
                       q=len(heap), secs=time.time() - t0)


resid, st = bb(V.phi_corner, T, max_boxes=2_000_000, time_budget=480.0)
print("t = psi   corner-bound B&B")
print("processed %d / cleared %d / infeasible %d / queue %d / residual %d  (%.0fs)"
      % (st["proc"], st["clr"], st["infeas"], st["q"], len(resid), st["secs"]))

if not resid:
    print("NO RESIDUALS -> would be a complete certificate (check independently)")
else:
    cen = np.array([[(lo + hi) / 2 for lo, hi in b] for b in resid])
    d_inf = np.max(np.abs(cen - OBS), axis=1)          # Chebyshev dist to obs
    widths = np.array([max(hi - lo for lo, hi in b) for b in resid])
    print("\nresidual centres: distance (Chebyshev) to obstruction (b,b,b,1,w*):")
    for p in (1, 10, 25, 50, 75, 90, 99):
        print("   p%2d  %8.4f" % (p, np.percentile(d_inf, p)))
    print("   frac within 0.02 of obstruction : %.3f" % np.mean(d_inf < 0.02))
    print("   frac within 0.10 of obstruction : %.3f" % np.mean(d_inf < 0.10))
    print("\nresidual box widths:  min %.2e  median %.2e  max %.2e"
          % (widths.min(), np.median(widths), widths.max()))
    print("residuals at min_diam (fully refined): %.3f"
          % np.mean(widths <= 1.01e-5))
    # how many residuals are FEASIBLE-straddling vs deep-interior?
    n = len(resid)
    print("\nfirst 8 residual centres (p1,q1,p2,q2,w) + dist + width:")
    idx = np.argsort(d_inf)[:8]
    for i in idx:
        print("   (%.4f,%.4f,%.4f,%.4f,%.4f)  d=%.4f  w=%.2e"
              % (cen[i, 0], cen[i, 1], cen[i, 2], cen[i, 3], cen[i, 4],
                 d_inf[i], widths[i]))
