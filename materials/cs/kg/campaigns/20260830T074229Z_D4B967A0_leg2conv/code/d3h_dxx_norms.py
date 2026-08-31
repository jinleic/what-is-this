"""d3h_dxx_norms — part (b): the three second-derivative norms of Phi3 via mpmath FD of the
sign-fixed closed-form assembly, sharing no numpy leg2 machinery:
  ||dxx Phi3||^2, ||dxy Phi3||^2, ||dyy Phi3||^2 at theta = 0, pi/2, pi on the GL grid points,
  computed by applying FD to the pointwise phi3_mc function on a 2-D tensor grid.
"""
import math, sys, time
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.legendre import leggauss
import mpmath
from mpmath import mp, mpf, mpc
mp.dps = 30
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


import d3h_kernels as K


def phi3_mc(x, y, t):
    r = rho_mc(t)
    s = 1 - r * r
    u = VTm * (x ** 3 - 3 * x) / msq(6)
    v = VTm * (y ** 3 - 3 * y) / msq(6)
    Kk = mexp(-(u * u - 2 * r * u * v + v * v) / (2 * s)) / (2 * mpmath.pi * msq(s))
    tot = mpc(0)
    ct = ct_mc(t)
    for (p, a) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
        cval = ct[(p, a)]
        if cval == 0:
            continue
        # G rows: a=0: kind K+r*(p-1); a>=1: dxdy^a entries with r^p folded
        if a == 0:
            rows = [((), (), "K" + "r" * (p - 1), 0)]
        else:
            import d3h_dxdy as DX
            from d3h_phi3 import G_pa_terms
            rows = G_pa_terms(p, a)
        for (wx, wy, kind, tp) in rows:
            mv = mpc(0)
            for (au, av, ar, pw), cc in K.mult_terms(kind).items():
                mv += mpc(repr(cc)) * (r ** ar) * (s ** (-pw)) * (u ** au) * (v ** av)
            tot += cval * 2 * mpmath.pi * (VTm ** tp) * weight_mc(wx, x) * weight_mc(wy, y) * mv * Kk
    return tot


def dphi3(x, y, t, which, h=mpf('0.01')):
    """second derivative of Phi3: which in {'xx','xy','yy'} via centered differences."""
    if which == 'xx':
        return (phi3_mc(x + h, y, t) - 2 * phi3_mc(x, y, t) + phi3_mc(x - h, y, t)) / (h * h)
    if which == 'yy':
        return (phi3_mc(x, y + h, t) - 2 * phi3_mc(x, y, t) + phi3_mc(x, y - h, t)) / (h * h)
    if which == 'xy':
        return (phi3_mc(x + h, y + h, t) - phi3_mc(x + h, y - h, t)
                - phi3_mc(x - h, y + h, t) + phi3_mc(x - h, y - h, t)) / (4 * h * h)
    raise ValueError(which)


def norm_fd(t, which, XW=mpf('3.0'), n=17):
    """trapezoid-in-tensor over [-3,3]^2 with weight dmu — the mass is far inside 3 (validated)."""
    xs = [mpf(-XW) + 2 * XW * mpf(i) / (n - 1) for i in range(n)]
    tot = mpc(0)
    for xi in xs:
        for yj in xs:
            dvals = dphi3(xi, yj, t, which)
            wt = mexp(-(xi * xi) / 2) / msq(2 * mpmath.pi) * mexp(-(yj * yj) / 2) / msq(2 * mpmath.pi)
            tot += wt * dvals * mpmath.conj(dvals)
    dx = (2 * XW) / (n - 1)
    return float((tot * dx * dx).real)


if __name__ == "__main__":
    print("(b) the three Second-derivative norms of Phi3 (independent mpmath FD route), vs leg2.")
    import d3h_sweep49_asm as A
    for tname, tval in (("theta=0", 1.0), ("theta=pi/2", 1j), ("theta=pi", -1.0)):
        t = mpc(math.cos(math.atan2(tval.imag if isinstance(tval, complex) else 0,
                                    tval.real if isinstance(tval, complex) else tval)))
        l1, l2 = A.legs_at(tval)
        nxx = norm_fd(t, 'xx')
        nxy = norm_fd(t, 'xy')
        nyy = norm_fd(t, 'yy')
        frob = nxx + 2 * nxy + nyy
        print(f"{tname}: leg1 = {l1:.6f}  leg2 = {l2:.6f}")
        print(f"   ||dxx||^2 = {nxx:.6f}  ||dxy||^2 = {nxy:.6f}  ||dyy||^2 = {nyy:.6f}")
        print(f"   Frobenius-ish (nxx+2nxy+nyy) = {frob:.6f}   leg2/Frobenius = {l2/frob:.6f}   leg2/nxy = {l2/nxy:.6f}")
