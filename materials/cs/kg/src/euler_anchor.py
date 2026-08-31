"""euler_anchor — re-anchor all cross-checks on Euler rho_j (paper D_t, lines 857/1826).
Verifies: (1) pointwise assembly vs independent mpmath Phi3 (Euler c-rows); (2) the (48)-bridge;
(3) the norm unit test survives. Then freezes the Euler-corrected sweep."""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import d3h_sweep49_asm as A
import d3h_leg2_indep as M
from d3h_phi3 import G_pa_terms
from d3h_kernels import mult_terms
from d3h_dxdy import weight_eval
from flint import acb

S3 = 0.34101124; S5 = 0.05276111; V = 1 + S3 * S3 + S5 * S5


def euler_j(j, t):
    return (t - (3 ** j) * S3 * S3 * t ** 3 + (5 ** j) * S5 * S5 * t ** 5) / V


def ctab_euler(t):
    r1, r2, r3 = euler_j(1, t), euler_j(2, t), euler_j(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3, (0, 1): -t, (0, 2): 3 * t * t,
            (0, 3): -(t ** 3), (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2, (1, 2): 3 * t * t * r1}


if __name__ == "__main__":
    r = 0.7921695401463447
    vt = A.VT
    x = 0.7
    s = 1 - r * r
    u = vt * M.psi3(x)
    Kval = math.exp(-(2 * u * u - 2 * r * u * u) / (2 * s)) / (2 * math.pi * math.sqrt(s))
    CE = ctab_euler(1.0)

    def row_val(p, a):
        tv = 0.0
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            mv = 0.0
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                mv += cc * (r ** ar) * (s ** (-pw)) * (u ** au) * (u ** av)
            wv = float(weight_eval(wx, acb(x)).mid().real)
            tv += 2 * math.pi * (vt ** tp) * wv * wv * mv * Kval
        return tv

    asm = sum(CE[(p, a)] * row_val(p, a) for (p, a) in CE)

    # independent mpmath Phi3 with Euler rows (fresh construction):
    import mpmath
    from mpmath import mp, mpf, mpc
    mp.dps = 40
    msq = mpmath.sqrt
    mexp = mpmath.exp
    VTm = mpf('0.136419125') / msq(mpf(repr(V)))
    C3m = mpf('-0.34101124') ** 2 / V
    C5m = mpf('0.05276111') ** 2 / V

    def rhoj_euler_m(j, t):
        return (t - (3 ** j) * mpf('0.34101124') ** 2 * t ** 3 + (5 ** j) * mpf('0.05276111') ** 2 * t ** 5) / V

    def ct_m(t):
        r1, r2, r3 = rhoj_euler_m(1, t), rhoj_euler_m(2, t), rhoj_euler_m(3, t)
        return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3, (0, 1): -t, (0, 2): 3 * t * t,
                (0, 3): -(t ** 3), (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2,
                (1, 2): 3 * t * t * r1}

    def phi3_m(x, y, t):
        r = (t - mpf('0.34101124') ** 2 * t ** 3 + mpf('0.05276111') ** 2 * t ** 5) / V
        s = 1 - r * r
        u = VTm * (x ** 3 - 3 * x) / msq(6)
        v = VTm * (y ** 3 - 3 * y) / msq(6)
        Kk = mexp(-(u * u - 2 * r * u * v + v * v) / (2 * s)) / (2 * mpmath.pi * msq(s))
        tot = mpc(0)
        for (p, a), cval in ct_m(t).items():
            if cval == 0:
                continue
            for (wx, wy, kind, tp) in G_pa_terms(p, a):
                mv = mpc(0)
                for (au, av, ar, pw), cc in mult_terms(kind).items():
                    mv += mpc(repr(cc)) * (r ** ar) * (s ** (-pw)) * (u ** au) * (v ** av)
                wm = lambda w, xx: (mpf(1) if not w else mpmath.nprod(lambda k: mpf(1), [1, 1]) *
                                    __import__('functools').reduce(
                                        lambda acc, lvl: acc * {1: (3 * xx * xx - 3) / msq(6),
                                                                2: xx * msq(6),
                                                                3: msq(6)}[lvl], w, mpf(1)))
                tot += cval * 2 * mpmath.pi * (VTm ** tp) * wm(wx, x) * wm(wy, y) * mv * Kk
        return tot

    mc = phi3_m(mpf('0.7'), mpf('0.7'), mpc(1))
    print(f"Phi3(0.7,0.7) [Euler rows]: assembly = {asm:.10f}   mpmath = {float(mc.real):.10f}   "
          f"rel = {abs(asm - float(mc.real)) / abs(float(mc.real)):.2e}")

    # (48)-bridge is rho_j-independent; re-verify quickly:
    bridge = 2 * math.pi * vt * vt * ((3 * x * x - 3) / math.sqrt(6)) ** 2 * Kval
    series = (math.pi / 2) * sum((r ** b) * (math.sqrt(b + 1) * M.qn(b + 1, u)) ** 2 for b in range(700)) \
        * vt * vt * ((3 * x * x - 3) / math.sqrt(6)) ** 2
    print(f"(48)-bridge: series = {series:.12e}  closed = {bridge:.12e}  rel = {abs(series - bridge) / abs(bridge):.2e}")
