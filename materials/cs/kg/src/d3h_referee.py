"""referee — third-route window sweep for leg1/leg2 at theta=0, adjudicating herm140 vs GL480."""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.legendre import leggauss
import d3h_leg2_indep as M

VT = 0.128957369921412
r = 0.7921695401463447


def legnorm_window(leg2, XW, P, m, B=220):
    edges = np.linspace(-XW, XW, P + 1)
    xg, wg = leggauss(m)
    nl, wl = [], []
    for i in range(P):
        a, b = edges[i], edges[i + 1]
        xm, xw = 0.5 * (a + b), 0.5 * (b - a)
        nds = xm + xw * xg
        wts = xw * wg * np.exp(-0.5 * nds ** 2) / math.sqrt(2 * math.pi)
        nl.append(nds); wl.append(wts)
    xs = np.concatenate(nl); ws = np.concatenate(wl)
    F = np.zeros((len(xs), len(xs)))
    nus = range(1, 5) if leg2 else range(0, 4)
    for nu in nus:
        D = (M.Dcoefs if leg2 else M.Dcoefs_leg1)(nu, 1.0)[:B]
        Q = np.array([[M.Qjet(nu, b, x)[nu] for b in range(B)] for x in xs])
        F += (math.pi / 2) * (Q * D) @ Q.T
    Wt = ws[:, None] * ws[None, :]
    return float(np.sum(Wt * F * F))


if __name__ == "__main__":
    for leg2 in (False, True):
        name = "leg2" if leg2 else "leg1"
        print(f"--- {name}(theta=0): GL composite, window sweep")
        for XW in (4.4, 6.0, 8.0, 11.0, 14.0):
            v = legnorm_window(leg2, XW, P=int(3 * XW), m=24)
            print(f"  XW={XW:5.1f}: {name}(0) = {v:.6f}")
        print(f"  herm140/trim14 (fully replicated legs49):",
              128.24802622928533 if not leg2 else 51.52475132200948)
        print(f"  GL480 on [-4.4,4.4] (first run):", 1.80717552 if not leg2 else 34.9601411)
