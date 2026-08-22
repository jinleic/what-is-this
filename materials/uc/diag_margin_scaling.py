"""How does certification cost scale with the margin c* - t?

The Margin Lemma gives min Lambda = 1.6298 (c* - t) to four digits, so a target
closer to psi has a LARGER margin and should clear at shallower depth.  Any
t > psi is already a new record, so the cheapest record is the smallest
positive offset.  This probe measures the actual tree size on one fixed hard
slab of the w in [15/16, 1] slice at several offsets, to decide which offset a
full campaign should target.

Sampled measurement of COST only; it proves nothing about Phi.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cert3
from entropy import PSI
from mpmath import mp

mp.dps = 40
C_STAR = mp.mpf("0.38234553336670272115")

# A narrow slab of the hardest slice (w near 1), small enough to finish.
SLAB = ((0.0, 1.0),) * 4 + ((0.9375, 0.9400),)
OFFSETS = ("0.0001", "0.00003", "0.00001")


def main():
    print("slab = %s" % (SLAB,), flush=True)
    for offset in OFFSETS:
        t = mp.mpf(PSI) + mp.mpf(offset)
        t_text = repr(float(t))
        margin = 1.6298 * float(C_STAR - t)
        started = time.monotonic()
        complete, stats = cert3.certify(
            t_text, SLAB, max_boxes=8_000_000, time_budget=5400.0,
            progress_every=0,
        )
        elapsed = time.monotonic() - started
        print("offset=%-9s t=%s margin=%.4e complete=%s processed=%d "
              "residual=%d ratio=%d center=%d elapsed=%.0fs"
              % (offset, t_text, margin, complete, stats["processed"],
                 stats["residual"], stats["ratio"], stats["center"], elapsed),
              flush=True)
    print("PROBE DONE", flush=True)


if __name__ == "__main__":
    main()
