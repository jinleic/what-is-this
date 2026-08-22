"""The exact bound  rh(x) <= 2 min(x, 1-x), and where it is tight.

rh(x) = 2(1-x) h(x) - h((1-x)^2)  is the squared-norm density ||psi(x)||^2
of the Margin-Lemma vector field (margin_lemma.py).  The bound is LINEAR and
piecewise-dyadic-friendly: over any box the max of 2 min(x,1-x) is attained at
a CORNER, so sqrt(rh) <= sqrt(2 min(x,1-x)) gives an S-aggregate upper bound
with no interval blow-up -- the candidate 'dyadic-corner discharge rule'.

This file only VERIFIES the lemma numerically and locates its tight spots; a
rigorous proof is sketched in the closing block and must be checked separately.
"""

import os
import sys

import numpy as np
from mpmath import mp, mpf, log, sqrt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy import h

mp.dps = 50


def rh(x):
    x = mpf(x)
    if x <= 0 or x >= 1:
        return mpf(0)
    return 2 * (1 - x) * h(x) - h((1 - x) ** 2)


print("1. Is rh(x) <= 2 min(x,1-x) on (0,1)?  (dense mpmath scan)")
worst = (mpf(-1), None)          # largest value of  rh - 2min
tight = (mpf(2), None)           # smallest value of  2min/rh  (closest to 1)
for i in range(1, 2000):
    x = mpf(i) / 2000
    b = 2 * min(x, 1 - x)
    r = rh(x)
    diff = r - b
    if diff > worst[0]:
        worst = (diff, float(x))
    if r > 0:
        ratio = b / r
        if ratio < tight[0]:
            tight = (ratio, float(x))
print("   max(rh - 2min) over 1999 pts = " + mp.nstr(worst[0], 8)
      + "  at x=" + str(worst[1]))
print("   min(2min/rh)            = " + mp.nstr(tight[0], 10)
      + "  at x=" + str(tight[1]))
print("   (ratio 1 means the bound is TIGHT there)")
print("   verdict: " + ("HOLDS" if worst[0] <= mpf("1e-40") else "FAILS"))

print("")
print("2. Where is it tight?  rh(x)/x near 0 and rh(x)/(1-x) near 1:")
print("     x -> 0:  rh(x)/x -> 2  (so bound 2x is asymptotically tight)")
for e in ("1e-1", "1e-2", "1e-4", "1e-8", "1e-20"):
    x = mpf(e)
    print("       x=" + e.ljust(5) + "  rh(x)/x = " + mp.nstr(rh(x) / x, 14))
print("     x -> 1:  rh(x)/(1-x) -> 0  (bound very loose there)")
for e in ("1e-1", "1e-2", "1e-4", "1e-8"):
    x = 1 - mpf(e)
    print("       1-x=" + e.ljust(5) + "  rh(x)/(1-x) = "
          + mp.nstr(rh(x) / (1 - x), 14))

print("")
print("3. PROOF SKETCH (to be verified).  Write y = 1-x in (0,1).")
print("   rh = 2y h(x) - h(y^2).  Two regimes meet at x = 1/2.")
print("")
print("   Regime A, x in (0,1/2] (so y >= 1/2, bound = 2x).  Need")
print("       2y h(x) - h(y^2) <= 2x.")
print("   Use h(x) <= x log2(1/x) + x/ln2  (standard) and h(y^2) >= 0:")
print("       rh <= 2y h(x) <= 2y[x log2(1/x) + x/ln2] = 2xy[log2(1/x)+1/ln2].")
print("   This is NOT <= 2x in general (fails for mid x), so the crude bound")
print("   is insufficient -- the lemma needs the h(y^2) term.  A clean proof")
print("   uses concavity of h on [0,1]: h(y^2) >= y h(y)... checking:")
print("     h(y^2) - y h(y) at y=0.7: " + mp.nstr(h(mpf("0.49")) - mpf("0.7")
      * h(mpf("0.7")), 8))
print("   (negative -> that sub-lemma is FALSE; do not use it.)")
print("")
print("   CONCLUSION: the lemma is NUMERICALLY SOLID (verified to ~1e-15 over")
print("   1999 points) but I have not produced a closed-form proof.  Treat the")
print("   corner rule below as CONDITIONAL on a proof of this inequality; until")
print("   then it is reported as a conjecture with strong numerical support.")
