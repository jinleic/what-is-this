"""Resolution probe: Phi_lo (corner) vs FULL box width around the obstruction.

box_around(c, delta) spans [v - delta/2, v + delta/2]: `delta` is the FULL box
width (half-width = delta/2).  An earlier ledger entry mislabelled it.
"""

import os
import sys

import cert as C
import diag_v2 as V
from entropy import PSI

B0 = 0.32945473850303697239
T = float(PSI)
a = (T - B0) / (1 - B0)
w0 = 1 - 2 * a
cen = (B0, B0, B0, 1.0, w0)
tv = float(C.phi_true(*cen))
print("centre (b,b,b,1,w*)  true Phi = %.6e" % tv)
print("(delta = FULL box width; half-width = delta/2)")
print("%-9s %13s" % ("delta", "Phi_lo (corner)"))


def box_around(c, delta):
    return tuple((max(0.0, v - delta / 2), min(1.0, v + delta / 2)) for v in c)


for delta in (1.0, 0.5, 0.1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
    box = box_around(cen, delta)
    e = V.phi_corner(box, T)
    if e is None:
        print("%-9g  infeasible-certified" % delta)
        continue
    print("%-9g  %+13.4e" % (delta, float(e.lower())))
