"""Diagnostic: WHERE do the t=psi residuals sit, and are they real?

Imports the SOUND enclosures from cert.py and runs a focused branch-and-bound
at t = psi, collecting every residual box's centre, width, the TRUE Phi at the
centre (mpmath), the enclosure lower bound, and the gap.  No new arithmetic --
this only reports what the existing bound leaves behind.

Reads:
  * if residual centres have TRUE Phi clearly > 0  -> they are enclosure
    looseness, not real violations; tightening the bound (lemma/centred) clears
    them.
  * if they cluster near p in {0,1}  -> sink zeros stalled it (lemma helps).
  * if they cluster near (b,b,b,1)    -> the obstruction stalled it (different).
"""

import os
import sys
import time

import numpy as np
from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert as C            # sound enclosures: phi_encl, phi_true
from entropy import PSI

mp.dps = 40
T = float(PSI)


def center(box):
    return tuple((lo + hi) / 2 for lo, hi in box)


def widths(box):
    return [hi - lo for lo, hi in box]


def bb_dump(t, max_boxes=60000, min_diam=1e-4, time_budget=200.0):
    import heapq
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
        e = C.phi_encl(box, t)
        if e is None:
            infeas += 1
            continue
        if e >= 0:
            cleared += 1
            continue
        if max(widths(box)) <= min_diam:
            resid.append(box)
            continue
        k = int(np.argmax(widths(box)))
        lo, hi = box[k]
        mid = (lo + hi) / 2
        for sub in ((lo, mid), (mid, hi)):
            nb = list(box)
            nb[k] = sub
            heapq.heappush(heap, (float(e.lower()), counter, tuple(nb)))
            counter += 1
    return resid, dict(proc=processed, clr=cleared, infeas=infeas,
                       q=len(heap), secs=time.time() - t0)


resid, st = bb_dump(T)
print("t = psi =", T)
print("processed/cleared/infeas/queue :", st["proc"], st["clr"], st["infeas"],
      st["q"], " residual boxes:", len(resid), "  (%.1fs)" % st["secs"])

# true Phi at each residual centre
rows = []
for b in resid:
    c = center(b)
    tv = float(C.phi_true(c[0], c[1], c[2], c[3], c[4]))
    e = C.phi_encl(b, T)
    lo = float("nan") if e is None else float(e.lower())
    rows.append((tv, lo, tv - lo, c, max(widths(b))))

rows.sort()
tv_arr = np.array([r[0] for r in rows])
gap = np.array([r[2] for r in rows])
cen = np.array([r[3] for r in rows])

print("\n-- TRUE Phi at residual centres (are they real violations?) --")
print("   min true Phi   : %.3e" % tv_arr.min())
print("   1%%  true Phi   : %.3e" % np.percentile(tv_arr, 1))
print("   median true Phi: %.3e" % np.median(tv_arr))
print("   => if min true Phi >> 0, residuals are ENCLOSURE LOOSENESS, not violations")

print("\n-- looseness: (true Phi - enclosure lower) --")
print("   median gap : %.3e" % np.median(gap))
print("   max gap    : %.3e" % gap.max())

print("\n-- WHERE do they cluster?  centre-coordinate percentiles --")
labs = ["p1", "q1", "p2", "q2", "w"]
for j, nm in enumerate(labs):
    col = cen[:, j]
    print("   %-3s  min %.3f  p10 %.3f  med %.3f  p90 %.3f  max %.3f"
          % (nm, col.min(), np.percentile(col, 10), np.median(col),
             np.percentile(col, 90), col.max()))

# distance to the obstruction (b,b,b,1) and to the sinks {0,1}
B0 = 0.32945473850303697
d_obs = np.max(np.abs(cen[:, :4] - np.array([B0, B0, B0, 1.0])), axis=1)
near_obs = np.mean(d_obs < 0.05)
# 'near a sink' if any of p1,q1,p2,q2 within 0.02 of 0 or 1
nr_edge = cen[:, :4]
near_sink = np.mean(np.any((nr_edge < 0.02) | (nr_edge > 0.98), axis=1))
print("\n-- clustering --")
print("   frac of residuals within 0.05 of obstruction (b,b,b,1): %.3f" % near_obs)
print("   frac with some coord within 0.02 of {0,1} (sink)      : %.3f" % near_sink)

print("\n-- the 12 TIGHTEST residuals (smallest true Phi) --")
print("   %-9s %-11s %-9s %-22s %s" %
      ("truePhi", "encLo", "width", "centre(p1,q1,p2,q2,w)", "tag"))
for tv, lo, g, c, w in rows[:12]:
    tag = ("OBS" if max(abs(c[0]-B0), abs(c[1]-B0), abs(c[2]-B0),
                       abs(c[3]-1)) < 0.05 else
           "SINK" if min(min(c[:4]), 1-max(c[:4])) < 0.02 else "mid")
    print("   %+9.3e %+11.3e %.2e (%.3f,%.3f,%.3f,%.3f,%.3f) %s"
          % (tv, lo, w, c[0], c[1], c[2], c[3], c[4], tag))
