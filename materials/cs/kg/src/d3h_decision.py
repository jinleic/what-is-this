"""decision — Main's decider: fix the rho/c-row defects in a CLEAN COPY of the legs49 assembly
(never mutating the frozen campaign) and test whether leg1(0)/leg2(0) land on the
independently-validated 1.807176 / 34.960447.

Variants at theta = 0 and theta = pi:
  A"as-frozen": cfun verbatim from legs49 (r2, r3 sign-flipped), r = rho_c (correct 0.79217).
      -> regression anchor: must reproduce 128.2480 / 51.5248.
  B"sign-fixed": c-rows from the paper's D-derivatives of eq (8) (C_TABLE-equivalent), same r.
      -> THE DECIDER: compare against 1.807176 / 34.960447.
  C"pointwise identity": with B's parameters, compare assembly Ph(x,y) against the independent
      Q-jet Phi3 at the same herm nodes (core region |x|,|y| <= 4): sup-norm difference.
"""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.hermite import hermgauss
from d3h_phi3 import G_pa_terms
from d3h_kernels import mult_terms
from d3h_dxdy import weight_eval
from flint import acb
import d3h_leg2_indep as M

S3 = 0.34101124
S5 = 0.05276111
V = 1 + S3 * S3 + S5 * S5
C1 = 1 / V
C3 = -(S3 * S3) / V
C5 = (S5 * S5) / V
VT = 0.136419125 / math.sqrt(V)

R_TRUE = VT + 0  # placeholder, real value below


def rho_paper(t):
    return C1 * t + C3 * t ** 3 + C5 * t ** 5


def rhoj_paper(j, t):
    if j == 0:
        return rho_paper(t)
    if j == 1:
        return C1 + 3 * C3 * t * t + 5 * C5 * t ** 4
    if j == 2:
        return 6 * C3 * t * t + 20 * C5 * t ** 4
    if j == 3:
        return 12 * C3 * t * t + 80 * C5 * t ** 4
    raise ValueError(j)


def c_paper(t):
    return {(1, 0): rhoj_paper(3, t), (2, 0): 3 * rhoj_paper(1, t) * rhoj_paper(2, t),
            (3, 0): rhoj_paper(1, t) ** 3, (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (rhoj_paper(1, t) + rhoj_paper(2, t)),
            (2, 1): -3 * t * rhoj_paper(1, t) ** 2, (1, 2): 3 * t * t * rhoj_paper(1, t)}


def c_broken(t):
    """legs49.cfun verbatim."""
    den = V
    r1 = (1 - 3 * S3 * S3 * t * t + 5 * S5 * S5 * t ** 4) / den
    r2 = (6 * S3 * S3 * t * t + 20 * S5 * S5 * t ** 4) / den
    r3 = (12 * S3 * S3 * t * t + 80 * S5 * S5 * t ** 4) / den
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3, (0, 1): -t, (0, 2): 3 * t * t,
            (0, 3): -(t ** 3), (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 * r1,
            (1, 2): 3 * t * t * r1}


NN = 140
yg, wg = hermgauss(NN)
xs_all = math.sqrt(2) * yg
w_all = wg / math.sqrt(math.pi)
msk = np.abs(xs_all) <= 14
xs, w = xs_all[msk], w_all[msk]
ux = VT * (xs ** 3 - 3 * xs) / math.sqrt(6)
Wij = np.outer(w, w)
_wv = {}


def wvec(wt):
    if wt not in _wv:
        _wv[wt] = np.array([float(weight_eval(wt, acb(x)).mid().real) for x in xs])
    return _wv[wt]


def assembly_leg(c, r, shift):
    """legs49-style assembly with given c-table rows and kernel parameter r (complex ok)."""
    s = 1 - r * r
    sr = complex(s).real if isinstance(s, complex) else float(np.real(s))
    Km = np.exp(-(ux[:, None] ** 2 - 2 * r * ux[:, None] * ux[None, :] + ux[None, :] ** 2) / (2 * s)) \
        / (2 * math.pi * np.sqrt(complex(s)))
    Ph = np.zeros((len(xs), len(xs)), dtype=complex)
    for (p, aa) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
        cval = c[(p, aa)]
        for (wx, wy, kind, tp) in G_pa_terms(p, aa + shift):
            Mmat = np.zeros((len(xs), len(xs)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (s ** (-pw)) if pw else 1.0
                if isinstance(sp, complex):
                    sp = complex(sp)
                up = (ux ** au) if au else np.ones_like(ux)
                vp = (ux ** av) if av else np.ones_like(ux)
                Mmat += cc * rp * sp * np.outer(up, vp)
            Ph += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
    return float(np.real(np.sum(Wij * np.abs(Ph) ** 2))), Ph


if __name__ == "__main__":
    r_t = rho_paper(1.0)
    print(f"rho_paper(1) = {r_t:.16f}   (paper eq (8): (1 - s3^2 + s5^2)/V)")
    print(f"legs49 rho_c(1) = {(1 - S3*S3 + S5*S5)/V:.16f}  identical: {abs(r_t - (1 - S3*S3 + S5*S5)/V) < 1e-15}")
    print(f"kernels/phi3 rho_x(1) = {1 + C3 + C5:.16f}  <- defective (missing 1/V on t)")
    print()
    for name, th in (("theta=0", 1.0), ("theta=pi", -1.0)):
        r = rho_paper(th)
        cA = c_broken(th)
        cB = c_paper(th)
        l1a, Pha = assembly_leg(cA, r, 0)
        l2a, _ = assembly_leg(cA, r, 1)
        l1b, Phb = assembly_leg(cB, r, 0)
        l2b, _ = assembly_leg(cB, r, 1)
        print(f"--- {name} (t = {th:+.0f}, r = rho({th:+.0f}) = {r:+.9f})")
        print(f"  A as-frozen  : leg1 = {l1a:.4f}   leg2 = {l2a:.4f}   (anchor 128.2480/51.5248 at theta=0)")
        print(f"  B sign-fixed : leg1 = {l1b:.4f}   leg2 = {l2b:.4f}   (target 1.807176/34.960447 at theta=0)")
        # C: pointwise identity check of assembly vs independent Q-jet Phi3 (core nodes)
        F = np.zeros((len(xs), len(xs)))
        B = 220
        core = np.abs(xs) <= 4.0
        idx_core = np.where(core)[0]
        Ds = []
        for nu in range(0, 4):
            D = np.zeros(B)
            for (p, aa) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
                if aa != nu:
                    continue
                cc = c_paper(th)[(p, aa)]
                for b in range(p, B):
                    fall = math.factorial(b) / math.factorial(b - p)
                    D[b] += cc * fall * (r ** (b - p))
            Ds.append(D)
        sup = 0.0
        worst = None
        for i in idx_core[::4]:
            for j in idx_core[::4]:
                val = 0.0
                for nu in range(0, 4):
                    for b in range(B):
                        jx = M.Qjet(nu, b, xs[i])[nu]
                        jy = M.Qjet(nu, b, xs[j])[nu]
                        val += (math.pi / 2) * Ds[nu][b] * jx * jy
                dd = abs(Phb[i, j].real - val)
                sc = max(abs(val), 1.0)
                if dd / sc > sup:
                    sup = dd / sc
                    worst = (xs[i], xs[j], Phb[i, j].real, val)
        print(f"  C pointwise sup|assembly - QjetPhi3| (|x|,|y| <= 4) = {sup:.3e}  at {worst[0]:+.3f},{worst[1]:+.3f}: {worst[2]:+.5f} vs {worst[3]:+.5f}")
        print()
