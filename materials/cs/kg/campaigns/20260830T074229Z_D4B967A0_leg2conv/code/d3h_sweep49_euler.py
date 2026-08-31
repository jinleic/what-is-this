"""sweep49_euler — DEFINITIVE corrected sweep: paper-printed c-table verbatim + Euler
rho_j = D_t^j rho (paper lines 857/1826: D = t d/dt), i.e. rho_j(t) = (t - 3^j s3^2 t^3 +
5^j s5^2 t^5)/V. This is the 3-line diff from d3h_sweep49_asm.py (ctab definition replaced).
RESULT (49-pt uniform trapezoid, theta-converged): avg(leg1) = 1.939156, avg(leg2) = 30.359613,
avg(S) = 123.377609, ratio to the paper's certified 126.80385221 = 0.972980 (inside).
All values COMPUTATIONAL-EVIDENCE (float64).
"""
import math, sys, time
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import d3h_leg2_indep as M   # noqa: F401 (used indirectly)

S3 = 0.34101124
S5 = 0.05276111
V = 1 + S3 * S3 + S5 * S5
VT = 0.136419125 / math.sqrt(V)


def rho_c(t):
    return (t - S3 * S3 * t ** 3 + S5 * S5 * t ** 5) / V


def rhoj(j, t):
    """Euler derivative D_t^j rho (paper D = t d/dt, line 857; rho_j := D_t^j rho, line 1826)."""
    return (t - (3 ** j) * S3 * S3 * t ** 3 + (5 ** j) * S5 * S5 * t ** 5) / V


def ctab(t):
    r1, r2, r3 = rhoj(1, t), rhoj(2, t), rhoj(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3, (0, 1): -t, (0, 2): 3 * t * t,
            (0, 3): -(t ** 3), (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2,
            (1, 2): 3 * t * t * r1}


IO = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def legs_at(t, XW=6.0, P=18, m=24, B=200):
    from numpy.polynomial.legendre import leggauss
    from d3h_phi3 import G_pa_terms
    from d3h_kernels import mult_terms
    from d3h_dxdy import weight_eval
    from flint import acb
    edges = np.linspace(-XW, XW, P + 1)
    xg, wg = leggauss(m)
    nl, wl = [], []
    for i in range(P):
        a_, b_ = edges[i], edges[i + 1]
        xm, xw = 0.5 * (a_ + b_), 0.5 * (b_ - a_)
        nds = xm + xw * xg
        nl.append(nds)
        wl.append(xw * wg * np.exp(-0.5 * nds ** 2) / math.sqrt(2 * math.pi))
    xs = np.concatenate(nl)
    ws = np.concatenate(wl)
    r = rho_c(t)
    s = 1 - r * r
    ux = VT * (xs ** 3 - 3 * xs) / math.sqrt(6)
    Km = np.exp(-(ux[:, None] ** 2 - 2 * r * ux[:, None] * ux[None, :] + ux[None, :] ** 2) / (2 * s)) \
        / (2 * math.pi * np.sqrt(complex(s)))
    Wt = ws[:, None] * ws[None, :]
    _wv = {}

    def wvec(wt):
        if wt not in _wv:
            _wv[wt] = np.array([float(weight_eval(wt, acb(x)).mid().real) for x in xs])
        return _wv[wt]

    ct = ctab(t.real if isinstance(t, float) else t)
    F1 = np.zeros((len(xs), len(xs)), dtype=complex)
    F2 = np.zeros((len(xs), len(xs)), dtype=complex)
    for (p, a) in IO:
        cval = ct[(p, a)]
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            Mmat = np.zeros((len(xs), len(xs)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (complex(s) ** (-pw)) if pw else 1.0
                up = (ux ** au) if au else np.ones_like(ux)
                vp = (ux ** av) if av else np.ones_like(ux)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F1 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
        for (wx, wy, kind, tp) in G_pa_terms(p, a + 1):
            Mmat = np.zeros((len(xs), len(xs)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (complex(s) ** (-pw)) if pw else 1.0
                up = (ux ** au) if au else np.ones_like(ux)
                vp = (ux ** av) if av else np.ones_like(ux)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F2 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
    l1 = float(np.real(np.sum(Wt * np.abs(F1) ** 2)))
    l2 = float(np.real(np.sum(Wt * np.abs(F2) ** 2)))
    return l1, l2


def sweep(n):
    th = np.linspace(0, math.pi, n + 1)
    L1 = np.zeros(n + 1)
    L2 = np.zeros(n + 1)
    t0 = time.time()
    for i, x in enumerate(th):
        L1[i], L2[i] = legs_at(complex(math.cos(x), math.sin(x)))
    h = th[1] - th[0]
    a1 = h * (L1[0] / 2 + L1[-1] / 2 + L1[1:-1].sum()) / math.pi
    a2 = h * (L2[0] / 2 + L2[-1] / 2 + L2[1:-1].sum()) / math.pi
    print(f"{n + 1}-pt: avg(leg1) = {a1:.6f}  avg(leg2) = {a2:.6f}  avg(S) = {a1 + 4 * a2:.6f}  [{time.time() - t0:.0f}s]")
    return a1, a2, a1 + 4 * a2, th, L1, L2


if __name__ == "__main__":
    print("rho_1(1) =", rhoj(1, 1.0), " rho_2(1) =", rhoj(2, 1.0), " rho_3(1) =", rhoj(3, 1.0))
    print("rho_0(i) =", rho_c(1j), " |rho_0(i)| =", abs(rho_c(1j)))
    for n in (12, 24, 48):
        sweep(n)
    a1, a2, Sv, th, L1, L2 = sweep(48)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2euler_th49.npy', th)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2euler_leg1_49.npy', L1)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2euler_leg2_49.npy', L2)
    print(f"ratio to 126.80385221 = {Sv / 126.80385221:.6f}")
    print(f"S(0) = {L1[0] + 4 * L2[0]:.4f}  S(pi/2) = {L1[24] + 4 * L2[24]:.4f}  S(pi) = {L1[48] + 4 * L2[48]:.4f}")
