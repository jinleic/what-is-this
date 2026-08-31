"""euler_fresh_anchor — CLOSURE 2: fresh independent evaluator for pointwise Phi3 with Euler
rows at theta=0, vs the committed assembly, to high order."""
import math, sys
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import d3h_sweep49_asm as A
from d3h_phi3 import G_pa_terms
from d3h_kernels import mult_terms
from flint import acb
import mpmath
from mpmath import mp, mpf, mpc
mp.dps = 40
msq = mpmath.sqrt
mexp = mpmath.exp

VTm = mpf('0.12895736992141200631702124692030')
S3m = mpf('0.34101124')
S5m = mpf('0.05276111')


def rho_e(t):
    return (t - S3m * S3m * t ** 3 + S5m * S5m * t ** 5) / (1 + S3m * S3m + S5m * S5m)


def rhoj_e(j, t):
    return (t - (3 ** j) * S3m * S3m * t ** 3 + (5 ** j) * S5m * S5m * t ** 5) / (1 + S3m * S3m + S5m * S5m)


def ct_e(t):
    r1, r2, r3 = rhoj_e(1, t), rhoj_e(2, t), rhoj_e(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3, (0, 1): -t, (0, 2): 3 * t * t,
            (0, 3): -(t ** 3), (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2,
            (1, 2): 3 * t * t * r1}


def phi3_e(x, y, t):
    r = rho_e(t)
    s = 1 - r * r
    u = VTm * (x ** 3 - 3 * x) / msq(6)
    v = VTm * (y ** 3 - 3 * y) / msq(6)
    KK = mexp(-(u * u - 2 * r * u * v + v * v) / (2 * s)) / (2 * mpmath.pi * msq(s))
    tot = mpc(0)
    for (p, a), cval in ct_e(t).items():
        if cval == 0:
            continue
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            mv = mpc(0)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                mv += mpc(repr(cc)) * (r ** ar) * (s ** (-pw)) * (u ** au) * (v ** av)
            wt = mpf(1)
            for lvl in wx:
                wt *= {1: (3 * x * x - 3) / msq(6), 2: x * msq(6), 3: msq(6)}[lvl]
            wt2 = mpf(1)
            for lvl in wy:
                wt2 *= {1: (3 * y * y - 3) / msq(6), 2: y * msq(6), 3: msq(6)}[lvl]
            tot += cval * 2 * mpmath.pi * (VTm ** tp) * wt * wt2 * mv * KK
    return tot


def phi3_asm_np(x, y):
    r = 0.7921695401463447
    s = 1 - r * r
    ux = A.VT * (x ** 3 - 3 * x) / math.sqrt(6)
    uy = A.VT * (y ** 3 - 3 * y) / math.sqrt(6)
    Km = math.exp(-(ux ** 2 - 2 * r * ux * uy + uy ** 2) / (2 * s)) / (2 * math.pi * math.sqrt(s))
    ct = A.ctab(1.0)
    tot = 0.0
    for (p, a), cval in ct.items():
        if cval == 0:
            continue
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            mv = 0.0
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                mv += cc * (r ** ar) * (s ** (-pw)) * (ux ** au) * (uy ** av)
            wv = 1.0
            for lvl in wx:
                wv *= float(A.wvec((lvl,))[list(A.XS).index(min(A.XS, key=lambda q: abs(q - x)))] * 0 +
                            (1 if lvl == 0 else {1: (3 * x * x - 3) / math.sqrt(6),
                                                 2: x * math.sqrt(6),
                                                 3: math.sqrt(6)}[lvl]))
            wv2 = 1.0
            for lvl in wy:
                wv2 *= (1 if lvl == 0 else {1: (3 * y * y - 3) / math.sqrt(6),
                                            2: y * math.sqrt(6),
                                            3: math.sqrt(6)}[lvl])
            tot += cval * (2 * math.pi) * (A.VT ** tp) * wv * wv2 * mv * Km
    return tot


if __name__ == "__main__":
    worst = 0.0
    for (x, y) in [(0.7, 0.7), (1.2, -0.9), (-2.0, 0.5)]:
        vm = complex(phi3_e(mpf(repr(x)), mpf(repr(y)), mpc(1)))
        va = phi3_asm_np(x, y)
        rel = abs(vm - va) / max(abs(vm), 1e-12)
        worst = max(worst, rel)
        print(f'({x:+.2f},{y:+.2f}): mpmath-fresh = {vm:+.12e}  assembly = {va:+.12e}  rel = {rel:.2e}')
    print(f'worst rel (theta=0, Euler rows): {worst:.2e}')
