"""sweep49_asm — corrected 49-pt S(theta) sweep via the CLOSED-FORM assembly (no b-series),
which is immune to the |r| = 1 boundary (theta ~ pi/2) where the b-series diverges.
Assembly = sign-fixed legs49 evaluator + C_TABLE rows + rho_true; GL nodes; vectorised per theta.
All values COMPUTATIONAL-EVIDENCE (float64).
"""
import math, sys, time
sys.path.insert(0, '/Users/jinleic/jinleic-workspace/cs/kg/src')
import numpy as np
from numpy.polynomial.legendre import leggauss
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


def rho_c(t):
    return C1 * t + C3 * t ** 3 + C5 * t ** 5




# GL grid (fixed, window 6.0 — the decider showed the integrand dies by |x| ~ 5)
XW = 6.0
P = 18
GM = 24
_edges = np.linspace(-XW, XW, P + 1)
_xg, _wg = leggauss(GM)
nl, wl = [], []
for _i in range(P):
    _a, _b = _edges[_i], _edges[_i + 1]
    _xm, _xw = 0.5 * (_a + _b), 0.5 * (_b - _a)
    _nds = _xm + _xw * _xg
    nl.append(_nds)
    wl.append(_xw * _wg * np.exp(-0.5 * _nds ** 2) / math.sqrt(2 * math.pi))
XS = np.concatenate(nl)
WS = np.concatenate(wl)


def rhoj(j, t):
    """Euler derivative D_t^j rho (paper lines 857/1826): rho_j(t) =
    (t - 3^j s3^2 t^3 + 5^j s5^2 t^5)/V. Plain derivatives were the defect."""
    return (C1 * t + (3 ** j) * C3 * t ** 3 + (5 ** j) * C5 * t ** 5)


def ctab(t):
    r1, r2, r3 = rhoj(1, t), rhoj(2, t), rhoj(3, t)
    return {(1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1 ** 3,
            (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 ** 2, (1, 2): 3 * t * t * r1}


def ctab_PLAIN_DERIVS(t):
    """DEFECTED historical variant (plain d^j/dt^j instead of Euler D_t) — provenance only.

    THIS VARIANT PRODUCED THE RETRACTED avg(S) = 379.718041 / 397.667200 AND THE VOIDED
    B3 = 24.9922. DO NOT USE. See campaign 20260830T074229Z_D4B967A0_leg2conv addenda."""
    return {(1, 0): _plain_rho_j(3, t), (2, 0): 3 * _plain_rho_j(1, t) * _plain_rho_j(2, t),
            (3, 0): _plain_rho_j(1, t) ** 3, (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -(t ** 3),
            (1, 1): -3 * t * (_plain_rho_j(1, t) + _plain_rho_j(2, t)), (2, 1): -3 * t * _plain_rho_j(1, t) ** 2,
            (1, 2): 3 * t * t * _plain_rho_j(1, t)}


def _plain_rho_j(j, t):
    if j == 1:
        return C1 + 3 * C3 * t * t + 5 * C5 * t ** 4
    if j == 2:
        return 6 * C3 * t * t + 20 * C5 * t ** 4
    if j == 3:
        return 12 * C3 * t * t + 80 * C5 * t ** 4
    raise ValueError(j)


UX = VT * (XS ** 3 - 3 * XS) / math.sqrt(6)
WIJ = np.outer(WS, WS)
_wv = {}


def wvec(wt):
    if wt not in _wv:
        _wv[wt] = np.array([float(weight_eval(wt, acb(x)).mid().real) for x in XS])
    return _wv[wt]


def legs_at(t):
    r = rho_c(t)
    s = 1 - r * r
    sq_s = np.sqrt(complex(s))
    Km = np.exp(-(UX[:, None] ** 2 - 2 * r * UX[:, None] * UX[None, :] + UX[None, :] ** 2) / (2 * s)) \
        / (2 * math.pi * sq_s)
    ct = ctab(t)
    F1 = np.zeros((len(XS), len(XS)), dtype=complex)
    F2 = np.zeros((len(XS), len(XS)), dtype=complex)
    for (p, a) in [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]:
        cval = ct[(p, a)]
        for (wx, wy, kind, tp) in G_pa_terms(p, a):
            Mmat = np.zeros((len(XS), len(XS)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (s ** (-pw)) if pw else 1.0
                up = (UX ** au) if au else np.ones_like(UX)
                vp = (UX ** av) if av else np.ones_like(UX)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F1 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
        for (wx, wy, kind, tp) in G_pa_terms(p, a + 1):
            Mmat = np.zeros((len(XS), len(XS)), dtype=complex)
            for (au, av, ar, pw), cc in mult_terms(kind).items():
                rp = (r ** ar) if ar else 1.0
                sp = (s ** (-pw)) if pw else 1.0
                up = (UX ** au) if au else np.ones_like(UX)
                vp = (UX ** av) if av else np.ones_like(UX)
                Mmat += cc * rp * sp * np.outer(up, vp)
            F2 += cval * (2 * math.pi) * (VT ** tp) * np.outer(wvec(wx), wvec(wy)) * Mmat * Km
    l1 = float(np.real(np.sum(WIJ * np.abs(F1) ** 2)))
    l2 = float(np.real(np.sum(WIJ * np.abs(F2) ** 2)))
    return l1, l2


def sweep(n, tag):
    th = np.linspace(0, math.pi, n + 1)
    L1 = np.zeros(n + 1)
    L2 = np.zeros(n + 1)
    t0 = time.time()
    for i, x in enumerate(th):
        L1[i], L2[i] = legs_at(complex(math.cos(x), math.sin(x)))
        if i % 4 == 0 or i == n:
            print(f"  {tag} th[{i}/{n}] = {x:.4f}: leg1 = {L1[i]:12.5f}  leg2 = {L2[i]:12.5f}", flush=True)
    h = th[1] - th[0]
    a1 = h * (L1[0] / 2 + L1[-1] / 2 + L1[1:-1].sum()) / math.pi
    a2 = h * (L2[0] / 2 + L2[-1] / 2 + L2[1:-1].sum()) / math.pi
    print(f"{tag}: avg(leg1) = {a1:.4f}  avg(leg2) = {a2:.4f}  avg(S) = {a1 + 4 * a2:.4f}  [{time.time() - t0:.0f}s]", flush=True)
    return a1, a2, a1 + 4 * a2, th, L1, L2


def assert_symmetries(S, th, tol=1e-10):
    """Standing tripwires (would have caught the rho_j defect on day one):
    (1) pi-periodicity: S(theta + pi) == S(theta)  [checked textually as S(0) == S(pi) on [0,pi]]
    (2) reflection about pi/2: S(theta) == S(pi - theta)  (even + pi-periodic), 24 independent
        equalities on the 49-pt grid. Both follow from rho_j oddness under the Euler operator."""
    n = len(S)
    ok = True
    # (1) endpoints:
    if abs(S[0] - S[-1]) > tol * max(1.0, abs(S[-1])):
        print(f"SYMMETRY VIOLATION: S(0) = {S[0]} vs S(pi) = {S[-1]}")
        ok = False
    # (2) reflection: for odd n, S[k] == S[n-1-k]
    worst = 0.0
    for k in range(n // 2 + 1):
        d = abs(S[k] - S[n - 1 - k]) / max(1.0, abs(S[n - 1 - k]))
        worst = max(worst, d)
    if worst > tol:
        print(f"REFLECTION VIOLATION: max |S(k) - S(n-1-k)|/|S| = {worst:.3e}")
        ok = False
    else:
        print(f"symmetries OK: S(0)=S(pi) and reflection about pi/2, worst rel = {worst:.2e}")
    return ok


if __name__ == "__main__":
    print("rho(1) =", rho_c(1.0), " |imax check: |rho(i)| =", abs(rho_c(1j)))
    print()
    a1_13, a2_13, S_13, *_ = sweep(12, "13-pt")
    print()
    a1_25, a2_25, S_25, *_ = sweep(24, "25-pt")
    print()
    a1_49, a2_49, S_49, th49, L1s, L2s = sweep(48, "49-pt")
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2asm_th49.npy', th49)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2asm_leg1_49.npy', L1s)
    np.save('/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/leg2asm_leg2_49.npy', L2s)
    print()
    B3sq = (math.pi ** 2 / 6) * S_49
    print(f"49-pt: avg(S) = {S_49:.4f}  B3^2 = {B3sq:.4f}  B3 = {math.sqrt(B3sq):.6f}")
    print(f"paper avg(S) = 126.80385221: ratio = {S_49 / 126.80385221:.6f}")
