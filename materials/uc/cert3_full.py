"""Full-cube certification driver at increasing levels t.

Runs cert3.certify (via run_bb) on the default symmetry-halved root
((0,1)^4 x [1/2,1]) with the complete rule chain: mean contraction ->
float-fast corner -> prefiltered ratio -> lambda-centered interior ->
q2-pin + centered face.  Prints per-level tallies and an honest verdict;
stops at the first INCOMPLETE level.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import ctx

from mpmath import mp

import cert2
import cert3
from arbcore import get_rh_gmax
from entropy import PSI

LEVELS = ("0.0001", "0.0002", "0.0003")
MAX_BOXES = 600_000_000
TIME_BUDGET = 10800.0
MIN_WIDTH = 1e-3
WORK_PREC = 80
if __name__ == "__main__":
    t0 = time.time()
    print("cert3_full: levels psi+%s, min FULL width %g, budget %d boxes / %.0fs each"
          % (", psi+".join(LEVELS), MIN_WIDTH, MAX_BOXES, TIME_BUDGET),
          flush=True)
    gmax = get_rh_gmax()
    rho_gmax = cert2.get_rho_gmax()
    print("global covers ready (%.0fs)" % (time.time() - t0), flush=True)
    # Covers above were built at 160 bits (their cached values remain valid).
    # The runtime bounds are sound at ANY precision (outward rounding); 80
    # bits only widens enclosures slightly and roughly doubles throughput.
    ctx.prec = WORK_PREC
    print("working precision set to %d bits" % WORK_PREC, flush=True)
    certified = []
    for off in LEVELS:
        t = mp.mpf(PSI) + mp.mpf(off)
        lambdas, lamhat = cert3.lambda_family(t)
        print("\n=== level t = psi + %s = %s ===\nlambda family %s (KKT %.9f)"
              % (off, mp.nstr(t, 30), lambdas, float(lamhat)), flush=True)
        ok, stats = cert3.run_bb(
            t, lambdas, gmax, rho_gmax,
            min_width=MIN_WIDTH, face_min_width=1e-3,
            max_boxes=MAX_BOXES, time_budget=TIME_BUDGET,
        )
        if ok:
            certified.append(off)
        else:
            print("stopping ladder: level psi+%s INCOMPLETE" % off, flush=True)
            break
    print("\ncert3_full RESULT: certified levels: %s"
          % (["psi+" + o for o in certified] or "NONE"), flush=True)
    print("total wall %.0fs" % (time.time() - t0), flush=True)
