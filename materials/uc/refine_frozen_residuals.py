"""DIAGNOSTIC ONLY: replay new sound rules on dumped residual boxes.

The input .npy stores float64 centres/full widths.  Endpoints are reconstructed
outward, so results are guidance; only a fresh immutable full campaign can
certify the original exact boxes.

Rule additions tested:
  1. exact orbit swap followed by the existing sound centered5_mixed bound;
  2. optional sink-safe w-only MVT bound from partial_center.centered_w_best;
  3. if an old terminal leaf still touches a sink collar, split its widest
     coordinate down to collar_min_width (pure refinement, hence sound).

Usage:
  refine_frozen_residuals.py LAUNCH_JSON RESIDUAL_NPY [MAX_BOXES] [USE_WEIGHT]
where USE_WEIGHT is 0 or 1 (default 1).
"""

import collections
import math
import os
import sys
from fractions import Fraction

import numpy as np

if len(sys.argv) not in (3, 4, 5):
    raise SystemExit(__doc__)
launch_path = os.path.abspath(sys.argv[1])
residual_path = os.path.abspath(sys.argv[2])
max_boxes = int(sys.argv[3]) if len(sys.argv) >= 4 else 2_000_000
use_weight = bool(int(sys.argv[4])) if len(sys.argv) == 5 else True
snapshot = os.path.join(os.path.dirname(launch_path), "snapshot")
sys.path.insert(0, snapshot)

from flint import ctx
from mpmath import mp

import cert2
import cert3
import cert3_par as driver
import diag_exhaust
from arbcore import get_rh_gmax, hull

# All proof functions resolve from the hash-pinned snapshot.
from entropy import PSI


def outward(value, direction, steps=4):
    for _ in range(steps):
        value = math.nextafter(value, direction)
    return value


def row_box(row):
    return tuple((max(0.0, outward(float(row[j] - row[5 + j] / 2),
                                     -math.inf)),
                  min(1.0, outward(float(row[j] + row[5 + j] / 2),
                                     math.inf)))
                 for j in range(5))


def swap_orbits(box):
    return cert3.orbit_swap_box(box)


def split(box, coordinate):
    lo, hi = box[coordinate]
    mid = (lo + hi) / 2
    left, right = list(box), list(box)
    left[coordinate] = (lo, mid)
    right[coordinate] = (mid, hi)
    return tuple(left), tuple(right)


launch = driver.verify_launch(launch_path)
t = mp.mpf(PSI) + mp.mpf(launch["level_offset"])
t_arb = cert3._arb(mp.nstr(t, 45))
t_upper = math.nextafter(float(t_arb.upper()), math.inf)
gmax = get_rh_gmax()
rho_gmax = cert2.get_rho_gmax()
ctx.prec = driver.WORK_PREC
lambdas, lamhat = cert3.lambda_family(t)
center_lams = (0.0, float(lamhat))
rows = np.load(residual_path)
assert rows.ndim == 2 and rows.shape[1] == 10
stack = [row_box(row) for row in rows]
stats = collections.Counter()
remaining = []
print("DIAGNOSTIC ONLY rows=%d max_boxes=%d hash=%s" %
      (len(stack), max_boxes, launch["code_sha256"]), flush=True)
while stack and stats["processed"] < max_boxes:
    box = stack.pop()
    stats["processed"] += 1
    box = diag_exhaust.mean_contract(box, t_upper)
    if box is None:
        stats["infeasible"] += 1
        continue
    corner = diag_exhaust.phi_corner(box, t_arb)
    if corner is None:
        stats["infeasible"] += 1
        continue
    if corner >= 0:
        stats["corner"] += 1
        continue
    w = hull(cert3._arb(repr(box[4][0])), cert3._arb(repr(box[4][1])))
    if cert2.ratio_rule(box[:4], t_arb, W=w, rho_gmax=rho_gmax):
        stats["ratio"] += 1
        continue
    # Existing full and oriented-mixed bounds.
    full = cert3.centered5_best(box, t_arb, center_lams, gmax)
    if full is not None and full >= 0:
        stats["full"] += 1
        continue
    mixed = cert3.centered5_mixed(box, t_arb, center_lams, gmax)
    if mixed is not None and mixed >= 0:
        stats["mixed"] += 1
        continue
    swapped = cert3.centered5_mixed_swap(
        box, t_arb, center_lams, gmax)
    if swapped is not None and swapped >= 0:
        stats["mixed_swap"] += 1
        continue
    if use_weight:
        weight = cert3.centered_w_best(box, t_arb, center_lams, gmax)
        if weight is not None and weight >= 0:
            stats["weight"] += 1
            continue
    widths = [hi - lo for lo, hi in box]
    collar_box = cert3._touches_collar(box)
    floor = driver.COLLAR_MIN_WIDTH if collar_box else driver.MIN_WIDTH
    candidates = [j for j, width in enumerate(widths) if width > floor]
    if not candidates:
        stats["residual"] += 1
        remaining.append(box)
        continue
    widest = max(candidates, key=widths.__getitem__)
    stack.extend(split(box, widest))
    stats["split"] += 1
    if stats["processed"] % 10000 == 0:
        print("processed=%d stack=%d stats=%s" %
              (stats["processed"], len(stack), dict(stats)), flush=True)

stats["stack"] = len(stack)
print("DIAGNOSTIC RESULT stats=%s" % dict(stats), flush=True)
for box in remaining[:20]:
    print("remaining", box)
assert stats["processed"] <= max_boxes
