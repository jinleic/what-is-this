"""kpow — discriminator: legs49 assembly with kernel power 1/2 (frozen) vs 1 (paper K^2 norm)."""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.hermite import hermgauss
import d3h_leg2_indep as M
from d3h_phi3 import G_pa_terms
from d3h_kernels import mult_terms
from d3h_dxdy import weight_eval
from flint import acb

VT = 0.136419125 / math.sqrt(1 + 0.34101124**2 + 0.05276111**2)


def rho_c(t):
    return (t - 0.34101124**2 * t**3 + 0.05276111**2 * t**5) / (1 + 0.34101124**2 + 0.05276111**2)


def cfun(pa, t):
    s3 = 0.34101124; s5 = 0.05276111
    den = 1 + s3*s3 + s5*s5
    r1 = (1 - 3*s3*s3*t*t + 5*s5*s5*t**4) / den
    r2 = (6*s3*s3*t*t + 20*s5*s5*t**4) / den
    r3 = (12*s3*s3*t*t + 80*s5*s5*t**4) / den
    return {(1, 0): r3, (2, 0): 3*r1*r2, (3, 0): r1**3, (0, 1): -t, (0, 2): 3*t*t, (0, 3): -(t**3),
            (1, 1): -3*t*(r1+r2), (2, 1): -3*t*r1*r1, (1, 2): 3*t*t*r1}[pa]


NN = 140
yg, wg = hermgauss(NN)
xs = math.sqrt(2) * yg
w = wg / math.sqrt(math.pi)
msk = np.abs(xs) <= 14
xs, w = xs[msk], w[msk]
ux = VT * (xs**3 - 3*xs) / math.sqrt(6)
Wij = np.outer(w, w)
_wv = {}


def wvec(wt):
    if wt not in _wv:
        _wv[wt] = np.array([float(weight_eval(wt, acb(x)).mid().real) for x in xs])
    return _wv[wt]


def leg2_at(th, kpow):
    t = complex(math.cos(th), math.sin(th))
    r = rho_c(t)
    s = 1 - r * r
    import cmath
    Km = np.exp(-(ux[:, None]**2 - 2*r*ux[:, None]*ux[None, :] + ux[None, :]**2) / (2*s)) / (2*math.pi*np.sqrt(complex(s)))
    # fractional power of a complex array: principal branch, like legs49's sqrt
    Kp = Km ** kpow
    Ph2 = np.zeros((len(xs), len(xs)), dtype=complex)
    for (p, a) in M.IO_LIST:
        cval = cfun((p, a), t)
        for (wx, wy, kind, tp) in G_pa_terms(p, a + 1):
            Mmat = np.zeros((len(xs), len(xs)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r**ar) if ar else 1.0
                sp = (s**(-pw)) if pw else 1.0
                up = (ux**au) if au else np.ones_like(ux)
                vp = (ux**av) if av else np.ones_like(ux)
                Mmat += cc*rp*sp*np.outer(up, vp)
            Ph2 += cval*(2*math.pi)*(VT**tp)*np.outer(wvec(wx), wvec(wy))*Mmat*Kp
    return float(np.real(np.sum(Wij*np.abs(Ph2)**2)))


if __name__ == "__main__":
    for th in (0.0, math.pi):
        for name, kpow in (("K^(1/2) [frozen legs49]", 0.5), ("K^1 [paper |G|^2 K^2]", 1.0)):
            v = leg2_at(th, kpow)
            print(f"theta={th:.4f}  {name}: leg2 = {v:.4f}")
        print()
