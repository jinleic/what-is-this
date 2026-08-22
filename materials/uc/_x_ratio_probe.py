"""Probe: why does the ratio rule not fire on the unresolved sink leaves?"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb
from mpmath import mp

import cert2
import cert3
from arbcore import ALPHA, ONE, TWO, get_rh_gmax
from entropy import PSI

T = mp.mpf(PSI) + mp.mpf("0.0001")
t_arb = cert3._arb(mp.nstr(T, 45))

LEAF = ((0.000977, 0.001953), (0.999023, 1.0),
        (0.000977, 0.001953), (0.002930, 0.003906),
        (0.514648, 0.515625))

gmax = get_rh_gmax()
rho_gmax = cert2.get_rho_gmax()
print("global rho cap:", float(arb(rho_gmax)))

atoms = LEAF[:4]
for i, (u, v) in enumerate(atoms):
    cap = cert2.rho_upper(u, v, rho_gmax)
    print("atom %d [%.6f, %.6f]: rho cap = %.6f" % (i, u, v, float(arb(cap))))

AM = cert3.hull((arb(repr(atoms[0][0])) + arb(repr(atoms[1][0]))) / TWO,
                (arb(repr(atoms[0][1])) + arb(repr(atoms[1][1]))) / TWO)
BM = cert3.hull((arb(repr(atoms[2][0])) + arb(repr(atoms[3][0]))) / TWO,
                (arb(repr(atoms[2][1])) + arb(repr(atoms[3][1]))) / TWO)
W = cert2._weight_range(AM, BM, t_arb)
print("AM = [%.6f, %.6f]  BM = [%.6f, %.6f]"
      % (float(AM.lower()), float(AM.upper()),
         float(BM.lower()), float(BM.upper())))
print("weight window W =", None if W is None else
      "[%.6f, %.6f]" % (float(W.lower()), float(W.upper())))

fire = cert2.ratio_rule(atoms, t_arb, rho_gmax=rho_gmax)
print("ratio_rule fires:", fire)

arho = (cert2.rho_upper(*atoms[0], rho_gmax)
        + cert2.rho_upper(*atoms[1], rho_gmax)) / TWO
brho = (cert2.rho_upper(*atoms[2], rho_gmax)
        + cert2.rho_upper(*atoms[3], rho_gmax)) / TWO
print("arho = %.6f  brho = %.6f" % (float(arb(arho)), float(arb(brho))))
if W is not None:
    wl, wh = W.lower(), W.upper()
    rlo = (ONE - wl) * brho + wl * arho
    rhi = (ONE - wh) * brho + wh * arho
    rho_hi = cert3.amax(rlo.upper(), rhi.upper())
    beta = ONE - ALPHA
    kappa = TWO * beta * (ONE - t_arb) - ONE
    margin = kappa - beta * rho_hi
    print("RHO_hi = %.6f  kappa = %.6f  margin = %+.6f"
          % (float(arb(rho_hi)), float(kappa), float(arb(margin))))
