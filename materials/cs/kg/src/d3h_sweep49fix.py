"""sweep49fix — corrected 49-point S(theta) sweep (AUTHORISED by Main post-decider).
Fixed parameters: rho(t) = (t - s3^2 t^3 + s5^2 t^5)/V (paper eq (8), line 687 printed form),
c-rows from D-derivatives of that rho (C_TABLE-equivalent, = d3h_decision.c_paper),
independent Q-jet evaluation with F-matrix (sum-then-square) norms.
Grids: 49-pt (pi/48), 25-pt (pi/24), 13-pt (pi/12) on identical machinery; trapezoid averages.
All numbers COMPUTATIONAL-EVIDENCE (float64).
"""
import math, sys, time
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
import d3h_leg2_indep as M

S3 = 0.34101124
S5 = 0.05276111
V = 1 + S3 * S3 + S5 * S5
C1 = 1 / V
C3 = -(S3 * S3) / V
C5 = (S5 * S5) / V
VT = 0.136419125 / math.sqrt(V)


def rho_c(t):
    return C1 * t + C3 * t ** 3 + C5 * t ** 5


def rhoj(j, t):
    if j == 0:
        return rho_c(t)
    if j == 1:
        return C1 + 3 * C3 * t * t + 5 * C5 * t ** 4
    if j == 2:
        return 6 * C3 * t * t + 20 * C5 * t ** 4
    if j == 3:
        return 12 * C3 * t * t + 80 * C5 * t ** 4
    raise ValueError(j)


def ctab(t):
    return {(1, 0): rhoj(3, t), (2, 0): 3 * rhoj(1, t) * rhoj(2, t), (3, 0): rhoj(1, t) ** 3,
            (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (rhoj(1, t) + rhoj(2, t)), (2, 1): -3 * t * rhoj(1, t) ** 2,
            (1, 2): 3 * t * t * rhoj(1, t)}


IO = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def Dnu(nu, t, B):
    """leg1 (nu=a) / leg2 (nu=a+1) coefficient arrays, COMPLEX at |t|=1."""
    r = rho_c(t)
    D = np.zeros(B, dtype=complex)
    for (p, a) in IO:
        nuv = a if nu < 4 and a == nu else (a + 1 if nu >= 1 and a + 1 == nu else -1)
    # do it directly:
    D = np.zeros(B, dtype=complex)
    for (p, a) in IO:
        pass
    return D


def D_arr(nus_wanted, t, B):
    """D^nu_b for the given nu value (nu>=0), aggregated over (p,a) rows."""
    r = rho_c(t)
    D = np.zeros(B, dtype=complex)
    for (p, a) in IO:
        if a == nus_wanted:
            c = ctab(t)[(p, a)]
            for b in range(p, B):
                fall = math.factorial(b) / math.factorial(b - p)
                D[b] += c * fall * (r ** (b - p)) if b >= p else 0
    return D


def D_arr2(nus_wanted, t, B):
    """leg2 nu = a+1 aggregation."""
    r = rho_c(t)
    D = np.zeros(B, dtype=complex)
    for (p, a) in IO:
        if a + 1 == nus_wanted:
            c = ctab(t)[(p, a)]
            for b in range(p, B):
                if b >= p:
                    fall = math.factorial(b) / math.factorial(b - p)
                    D[b] += c * fall * (r ** (b - p))
    return D


def legs_at(t, XW=6.0, P=18, m=24, B=200):
    """leg1, leg2 at complex t on the circle via Q-jet F-matrix on GL grid."""
    edges = np.linspace(-XW, XW, P + 1)
    from numpy.polynomial.legendre import leggauss
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
    Wt = ws[:, None] * ws[None, :]
    F1 = np.zeros((len(xs), len(xs)), dtype=complex)
    F2 = np.zeros((len(xs), len(xs)), dtype=complex)
    for nu in range(0, 4):
        D = D_arr(nu, t, B)
        Q = np.array([[M.Qjet(nu, b, x)[nu] for b in range(B)] for x in xs])
        F1 += (math.pi / 2) * (Q * D) @ Q.T
    for nu in range(1, 5):
        D = D_arr2(nu, t, B)
        Q = np.array([[M.Qjet(nu, b, x)[nu] for b in range(B)] for x in xs])
        F2 += (math.pi / 2) * (Q * D) @ Q.T
    l1 = float(np.real(np.sum(Wt * np.abs(F1) ** 2)))
    l2 = float(np.real(np.sum(Wt * np.abs(F2) ** 2)))
    return l1, l2


def sweep(n_panels_theta):
    th = np.linspace(0, math.pi, n_panels_theta + 1)
    L1 = np.zeros(len(th))
    L2 = np.zeros(len(th))
    t0 = time.time()
    for i, x in enumerate(th):
        L1[i], L2[i] = legs_at(complex(math.cos(x), math.sin(x)))
        print(f"  th[{i}/{len(th)-1}] = {x:.4f}: leg1 = {L1[i]:10.4f}  leg2 = {L2[i]:10.4f}", flush=True)
    h = th[1] - th[0]
    a1 = h * (L1[0] / 2 + L1[-1] / 2 + L1[1:-1].sum()) / math.pi
    a2 = h * (L2[0] / 2 + L2[-1] / 2 + L2[1:-1].sum()) / math.pi
    print(f"grid {n_panels_theta + 1}-pt: avg(leg1) = {a1:.4f}  avg(leg2) = {a2:.4f}  "
          f"avg(S) = {a1 + 4 * a2:.4f}  [{time.time() - t0:.0f}s]")
    return a1, a2, a1 + 4 * a2, th, L1, L2


if __name__ == "__main__":
    print("CONSTANTS: rho(1) =", rho_c(1.0), " rho'1 =", rhoj(1, 1.0))
    print("counter-check vs paper eq (8) at 120 bits: (1 - s3^2 + s5^2)/V =", (1 - S3 * S3 + S5 * S5) / V)
    print()
    print("=== 13-pt sweep ===")
    a1_13, a2_13, S_13, *_ = sweep(12)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2fix_th13.npy',
            np.linspace(0, math.pi, 13))
    print()
    print("=== 25-pt sweep ===")
    a1_25, a2_25, S_25, *_ = sweep(24)
    print()
    print("=== 49-pt sweep ===")
    a1_49, a2_49, S_49, th49, L1, L2 = sweep(48)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2fix_th49.npy', th49)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2fix_leg1_49.npy', L1)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2fix_leg2_49.npy', L2)
    print()
    B3sq = (math.pi ** 2 / 6) * S_49
    print(f"49-pt: avg(S) = {S_49:.4f}  B3^2 = {B3sq:.4f}  B3 = {math.sqrt(B3sq):.6f}")
    print(f"paper certificate comparison: avg(S) vs 126.80385221: excess factor = {S_49 / 126.80385221:.6f}")
