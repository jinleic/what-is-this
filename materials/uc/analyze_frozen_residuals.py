"""DIAGNOSTIC ONLY: analyze a frozen campaign's dumped residual leaves.

The .npy format stores float64 centres and full widths, not exact endpoints.
Boxes are reconstructed with four outward ``nextafter`` steps per endpoint;
all conclusions are NUMERICAL GUIDANCE, never a proof artifact.  A subsequent
full immutable campaign must certify any rule suggested here.

Usage: analyze_frozen_residuals.py LAUNCH_JSON RESIDUAL_NPY
"""

import collections
import json
import math
import os
import sys

import numpy as np

if len(sys.argv) != 3:
    raise SystemExit(__doc__)

launch_path = os.path.abspath(sys.argv[1])
residual_path = os.path.abspath(sys.argv[2])
snapshot = os.path.join(os.path.dirname(launch_path), "snapshot")
sys.path.insert(0, snapshot)

from flint import ctx
from mpmath import mp

import cert2
import cert3
import cert3_par as driver
import partial_center
from arbcore import get_rh_gmax
from entropy import PSI


def outward(value, direction, steps=4):
    for _ in range(steps):
        value = math.nextafter(value, direction)
    return value


def row_box(row):
    centre, width = row[:5], row[5:]
    return tuple((max(0.0, outward(float(centre[j] - width[j] / 2), -math.inf)),
                  min(1.0, outward(float(centre[j] + width[j] / 2), math.inf)))
                 for j in range(5))


def swap_orbits(box):
    wl, wh = box[4]
    return (box[2], box[3], box[0], box[1], (1.0 - wh, 1.0 - wl))


launch = driver.verify_launch(launch_path)
rows = np.load(residual_path)
assert rows.ndim == 2 and rows.shape[1] == 10
print("DIAGNOSTIC ONLY -- NUMERICAL GUIDANCE, NOT A CERTIFICATE")
print("rows=%d source=%s code_sha256=%s" %
      (len(rows), residual_path, launch["code_sha256"]))
if not len(rows):
    raise SystemExit("empty residual array")

centres = rows[:, :5]
widths = rows[:, 5:]
names = ("p1", "q1", "p2", "q2", "w")
print("centre min:", dict(zip(names, np.min(centres, axis=0))))
print("centre max:", dict(zip(names, np.max(centres, axis=0))))
print("width values/counts:")
for j, name in enumerate(names):
    values, counts = np.unique(widths[:, j], return_counts=True)
    print("  %s: %s" % (name, list(zip(values.tolist(), counts.tolist()))))
edge_distance = np.minimum(centres[:, :4], 1.0 - centres[:, :4])
print("nearest atom-edge distance quantiles:",
      np.quantile(np.min(edge_distance, axis=1), [0, .1, .5, .9, 1]).tolist())

boxes = [row_box(row) for row in rows]
classes = collections.Counter(cert3._classify(box, float(mp.mpf(PSI) +
                                                       mp.mpf(launch["level_offset"])))
                              for box in boxes)
touch_patterns = collections.Counter(
    tuple(names[j] for j in range(4)
          if box[j][0] <= cert3.COLLAR_EDGE or
          box[j][1] >= 1.0 - cert3.COLLAR_EDGE)
    for box in boxes)
print("classes:", dict(classes))
print("collar-coordinate patterns:", dict(touch_patterns))
print("q2 pin-eligible:", sum(box[3][0] >= cert3.Q2PIN for box in boxes))

# Evaluate three sound lower-bound formulas on outward-expanded reconstructed
# boxes.  This is still diagnostic because the source .npy lost exact endpoints.
t = mp.mpf(PSI) + mp.mpf(launch["level_offset"])
gmax = get_rh_gmax()
ctx.prec = driver.WORK_PREC
_, lamhat = cert3.lambda_family(t)
center_lams = (0.0, float(lamhat))
t_arb = cert3._arb(mp.nstr(t, 45))
counts = collections.Counter()
worst = {"full": math.inf, "mixed": math.inf, "mixed_swap": math.inf,
         "weight": math.inf}
examples = []
true_phi = []
for i, (row, box) in enumerate(zip(rows, boxes)):
    values = {
        "full": cert3.centered5_best(box, t_arb, center_lams, gmax),
        "mixed": cert3.centered5_mixed(box, t_arb, center_lams, gmax),
        "mixed_swap": cert3.centered5_mixed(swap_orbits(box), t_arb,
                                             center_lams, gmax),
        "weight": partial_center.centered_w_best(
            box, t_arb, center_lams, gmax),
    }
    fired = []
    for name, value in values.items():
        if value is not None:
            lo = float(value.lower())
            worst[name] = min(worst[name], lo)
            if value >= 0:
                counts[name] += 1
                fired.append(name)
        else:
            counts[name + "_none"] += 1
    if fired:
        counts["any"] += 1
    elif len(examples) < 12:
        examples.append((i, row.tolist(),
                         {k: None if v is None else float(v.lower())
                          for k, v in values.items()}))
    true_phi.append(float(cert2.phi_true(*map(mp.mpf, row[:5]))))
    if (i + 1) % 100 == 0:
        print("bounds %d/%d any=%d full=%d mixed=%d mixed_swap=%d weight=%d" %
              (i + 1, len(rows), counts["any"], counts["full"],
               counts["mixed"], counts["mixed_swap"], counts["weight"]),
              flush=True)
print("bound firing counts:", dict(counts))
print("minimum finite lower bounds:", worst)
print("centre true-Phi quantiles (NUMERICAL):",
      np.quantile(true_phi, [0, .01, .1, .5, .9, 1]).tolist())
print("first unresolved examples:")
for example in examples:
    print(json.dumps(example))
