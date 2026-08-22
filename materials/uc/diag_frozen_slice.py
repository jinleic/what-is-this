"""DIAGNOSTIC ONLY: truncate one frozen cert3 slice and dump residual leaves.

This does not produce or modify proof results.  It imports the immutable
campaign snapshot, replays its exact rule chain, stops at a chosen box budget,
and saves centres plus full widths for terminal residual leaves.

Usage: diag_frozen_slice.py LAUNCH_JSON SLICE MAX_BOXES SECONDS OUTPUT_NPY
"""

import json
import os
import sys
from fractions import Fraction

if len(sys.argv) != 6:
    raise SystemExit(__doc__)

launch_path = os.path.abspath(sys.argv[1])
slice_index = int(sys.argv[2])
max_boxes = int(sys.argv[3])
seconds = float(sys.argv[4])
output = os.path.abspath(sys.argv[5])
snapshot = os.path.join(os.path.dirname(launch_path), "snapshot")
sys.path.insert(0, snapshot)

from flint import ctx
from mpmath import mp

import cert2
import cert3
import cert3_par as driver
from arbcore import get_rh_gmax
from entropy import PSI

launch = driver.verify_launch(launch_path)
run = launch["runs"][slice_index]
assert run["slice"] == slice_index
assert not os.path.exists(output), "refusing to overwrite diagnostic output"

t = mp.mpf(PSI) + mp.mpf(launch["level_offset"])
gmax = get_rh_gmax()
rho_gmax = cert2.get_rho_gmax()
ctx.prec = driver.WORK_PREC
lambdas, _ = cert3.lambda_family(t)
assert [repr(float(value)) for value in lambdas] \
    == launch["parameters"]["lambda_family"]
lo = float(Fraction(run["w_lo"]))
hi = float(Fraction(run["w_hi"]))
root = ((0.0, 1.0),) * 4 + ((lo, hi),)
print("DIAGNOSTIC ONLY slice=%d w=[%s,%s] max_boxes=%d seconds=%g hash=%s" %
      (slice_index, run["w_lo"], run["w_hi"], max_boxes, seconds,
       launch["code_sha256"]), flush=True)
complete, stats = cert3.certify(
    t,
    root,
    max_boxes=max_boxes,
    time_budget=seconds,
    min_width=driver.MIN_WIDTH,
    face_min_width=driver.FACE_MIN_WIDTH,
    collar_min_width=driver.COLLAR_MIN_WIDTH,
    gmax=gmax,
    rho_gmax=rho_gmax,
    lambdas=lambdas,
    progress_every=1_000_000,
    residual_dump=output,
)
print("DIAGNOSTIC VERDICT (NOT A CERTIFICATE): complete=%s stats=%s output=%s" %
      (complete, json.dumps(dict(stats), sort_keys=True), output), flush=True)
