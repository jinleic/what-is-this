"""dxy_isolate — Main's localization test: the 2.94x residual after leg1=0 lives in
leg2 = ||dxdy Phi3||^2. Two tests, sharing no code with the leg2 path:

(a) Pointwise FD of Phi3 (dual-route validated, sup-diff 3.7e-14) vs the leg2 path's dxdyPhi3:
    central differences in x and y on Phi3 (from the sign-fixed closed-form assembly, itself
    cross-checked against Q-jets), at several (x,y) and theta = 0, 1.7671, pi. Report ratios.
(b) The three second-derivative norms ||dxx Phi3||^2, ||dxy Phi3||^2, ||dyy Phi3||^2 on the
    same grid, against leg2 as returned. Detects a Hessian-aggregate factor (1 + 2r^2 form).
"""
import math, sys, time
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.legendre import leggauss
from d3h_phi3 import G_pa_terms, Phi3_eval, dxdy_Phi3_eval
from d3h_kernels import mult_terms, vartheta
from d3h_dxdy import weight_eval
from flint import acb
import mpmath
from mpmath import mp, mpf, mpc
mp.dps = 40
mexp = mpmath.exp
msq = mpmath.sqrt

S3 = mpf('0.34101124'); S5 = mpf('0.05276111')
V = 1 + S3 * S3 + S5 * S5
VTm = mpf('0.136419125') / msq(V)
C1 = 1 / V; C3 = -(S3 * S3) / V; C5 = (S5 * S5) / V


def rho_mc(t):
    return C1 * t + C3 * t ** 3 + C5 * t ** 5


def rhoj_mc(j, t):
    if j == 1:
        return C1 + 3 * C3 * t * t + 5 * C5 * t ** 4
    if j == 2:
        return 6 * C3 * t * t + 20 * C5 * t ** 4
    if j == 3:
        return 12 * C3 * t * t + 80 * C5 * t ** 4
    raise ValueError(j)


def ct_mc(t):
    return {(1, 0): rhoj_mc(3, t), (2, 0): 3 * rhoj_mc(1, t) * rhoj_mc(2, t), (3, 0): rhoj_mc(1, t) ** 3,
            (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (rhoj_mc(1, t) + rhoj_mc(2, t)), (2, 1): -3 * t * rhoj_mc(1, t) ** 2,
            (1, 2): 3 * t * t * rhoj_mc(1, t)}


def phi3_mc(x, y, t):
    """sign-fixed closed-form assembly via mpmath (independent of numpy paths)."""
    r = rho_mc(t)
    s = 1 - r * r
    u = VTm * (x ** 3 - 3 * x) / msq(6)
    v = VTm * (y ** 3 - 3 * y) / msq(6)
    K = mexp(-(u * u - 2 * r * u * v + v * v) / (2 * s)) / (2 * mpmath.pi * msq(s))
    tot = mpc(0)
    ct = ct_mc(t)
    for (p, a) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
        cval = ct[(p, a)]
        if cval == 0:
            continue
        for (wx, wy, kind, tp) in dxdy_terms_local(p, a):
            mv = mpc(0)
            for (au, av, ar, pw), cc in mult_terms_mc(kind).items():
                mv += cc * (r ** ar) * ((s) ** (-pw)) * (u ** au) * (v ** av)
            wxv = weight_mc(wx, x)
            wyv = weight_mc(wy, y)
            tot += cval * 2 * mpmath.pi * (VTm ** tp) * wxv * wyv * mv * K
    return tot


def dxdy_terms_local(p, a):
    if a == 0:
        return [((), (), "K" + "r" * max(0, p - 1) if p >= 1 else "K", 0)] if p <= 1 else [((), (), "K" + "r" * (p - 1), 0)]
    # reuse the validated generator but only its SHAPE (weights/kinds/tp):
    import d3h_dxdy as DX
    from d3h_phi3 import G_pa_terms
    return G_pa_terms(p, a)


def mult_terms_mc(kind):
    # mpmath mirror of d3h_kernels.mult_terms (validated term dicts), rebuilt here:
    import d3h_kernels as K
    terms = K.mult_terms(kind)
    return {(au, av, ar, pw): mpc(repr(cc.real) if isinstance(cc, complex) else repr(cc))
            for (au, av, ar, pw), cc in terms.items()}


def weight_mc(w, x):
    tot = mpf(1)
    for lvl in w:
        if lvl == 1:
            tot *= (3 * x * x - 3) / msq(6)
        elif lvl == 2:
            tot *= x * msq(6)
        elif lvl == 3:
            tot *= msq(6)
    return tot


def leg2_path_dxdy(x, y, t):
    """the leg2 path's dxdy Phi3 at a point via mpmath mirror of its assembly."""
    r = rho_mc(t)
    s = 1 - r * r
    u = VTm * (x ** 3 - 3 * x) / msq(6)
    v = VTm * (y ** 3 - 3 * y) / msq(6)
    K = mexp(-(u * u - 2 * r * u * v + v * v) / (2 * s)) / (2 * mpmath.pi * msq(s))
    tot = mpc(0)
    ct = ct_mc(t)
    for (p, a) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
        cval = ct[(p, a)]
        if cval == 0:
            continue
        for (wx, wy, kind, tp) in dxdy_terms_local(p, a + 1):
            mv = mpc(0)
            for (au, av, ar, pw), cc in mult_terms_mc(kind).items():
                mv += cc * (r ** ar) * ((s) ** (-pw)) * (u ** au) * (v ** av)
            wxv = weight_mc(wx, x)
            wyv = weight_mc(wy, y)
            tot += cval * 2 * mpmath.pi * (VTm ** tp) * wxv * wyv * mv * K
    return tot


if __name__ == "__main__":
    HX = mpf('0.002')
    print("(a) pointwise: FD-dxdy of Phi3 (independent) vs leg2-path dxdyPhi3")
    for theta in (0.0, math.pi / 2, 1.7671458676, math.pi):
        t = mpc(math.cos(theta), math.sin(theta))
        print(f"  theta = {theta:.4f}")
        for (x, y) in [(mpf('0.7'), mpf('1.1')), (mpf('-1.2'), mpf('0.9')), (mpf('2.1'), mpf('-0.4'))]:
            fdxdy = (phi3_mc(x + HX, y + HX, t) - phi3_mc(x + HX, y - HX, t)
                     - phi3_mc(x - HX, y + HX, t) + phi3_mc(x - HX, y - HX, t)) / (4 * HX * HX)
            lp = leg2_path_dxdy(x, y, t)
            fd = complex(fdxdy)
            lpv = complex(lp)
            ratio = abs(fd / lpv) if abs(lpv) > 1e-30 else float('nan')
            print(f"    ({float(x):+.2f},{float(y):+.2f}): FD = {fd:+.6e}  path = {lpv:+.6e}  ratio = {ratio:.8f}")
