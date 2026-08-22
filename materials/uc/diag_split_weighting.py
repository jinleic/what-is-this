"""Does influence-weighted split selection shrink the near-degenerate slices?

Atom coordinates enter Phi only through their own orbit's weight, so on the
w-slices near 1 the second orbit's coordinates are discounted by 1-w.  Choosing
the split coordinate by width*influence instead of raw width is a pure cost
heuristic: any split is sound because the two children cover the parent
exactly.  This probe times the same slab under both rules via the
``cert3.SPLIT_CHOICE`` hook.

Sampled COST measurement only; it proves nothing about Phi.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cert3
from entropy import PSI
from mpmath import mp

mp.dps = 40
T_TEXT = repr(float(mp.mpf(PSI) + mp.mpf("0.0001")))

# Slabs of the hardest slice (w near 1).  The first did not finish in 5400 s
# under raw-width splitting.
SLABS = (
    ("w[0.9375,0.9400]", ((0.0, 1.0),) * 4 + ((0.9375, 0.9400),)),
    ("w[0.9900,1.0000]", ((0.0, 1.0),) * 4 + ((0.9900, 1.0000),)),
)


def main():
    print("t = %s" % T_TEXT, flush=True)
    try:
        for label, slab in SLABS:
            for name, choice in (("weighted", cert3.split_by_influence),
                                 ("raw", cert3.split_by_width)):
                cert3.SPLIT_CHOICE = choice
                assert cert3.SPLIT_CHOICE is choice
                started = time.monotonic()
                complete, stats = cert3.certify(
                    T_TEXT, slab, max_boxes=6_000_000, time_budget=2700.0,
                    progress_every=0,
                )
                elapsed = time.monotonic() - started
                print("%s %-9s complete=%-5s processed=%9d residual=%d "
                      "elapsed=%.0fs"
                      % (label, name, complete, stats["processed"],
                         stats["residual"], elapsed), flush=True)
    finally:
        cert3.SPLIT_CHOICE = cert3.split_by_influence
    print("PROBE DONE", flush=True)


if __name__ == "__main__":
    main()
