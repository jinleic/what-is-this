"""d3h_circle.py — S(θ) probes at t = e^{iθ} (unit circle, complex), float64 stage.

⟨Gi,Gj⟩ for complex t: G = 2π ϑ^tp Wx(x) Wy(y) m_[kind](u,v; r) K_r(u,v) with r = ρ(t) complex.
inner = (1/|s|)·Σ_n γ^n/n!·Σ_pairs [C(m1,m2)]·Mn_x·Mn_y with
  C(m1,m2) = Σ c1·conj-term pairing:
    G1·conj(G2) → m1(u,v;r)·conj(m2(u,v;r)) = m1(u,v;r)·m2(u,v;r̄) (real coefficients) and K_r·K_r̄
    = (4π²|s|)^{-1} exp(−α(u²+v²) + 2γuv), α = Re(1/s), γ = Re(r/s).
  m-side expansion: m1·m2(r̄) = Σ c1c2 u^{au1+au2} v^{av1+av2} r^{ar1} r̄^{ar2} s^{−pw1} s̄^{−pw2}
    = Σ c1c2 u^{a} v^{b} Z where Z = r^{ar1}r̄^{ar2}s^{−pw1}s̄^{−pw2} (complex scalar).
  So ⟨Gi,Gj⟩ = (1/|s|)·Σ_entries 2πϑ^{tp1}·2πϑ^{tp2}·conj^* bookkeeping … assemble in code as:
     Σ_entries Σ_n γ^n/n! Mn(au+n) Mn(av+n) × Z_complex × real-Mn-products — final complex value;
  for i == j (self inner) the result must be REAL ≥ 0 — assertion in code.
"""
import math, sys, cmath
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
from numpy.polynomial.hermite import hermgauss
from d3h_kernels import mult_terms
from d3h_phi3 import G_pa_terms, _IO_LIST

V = 1 + 0.34101124**2 + 0.05276111**2
VT0 = 0.136419125 / math.sqrt(V)

yg, wg = hermgauss(170)
XS = math.sqrt(2) * yg
WTS = wg / math.sqrt(math.pi)
mask = np.abs(XS) <= 12.0
XSn, WTSn = XS[mask], WTS[mask]
US = VT0 * (XSn**3 - 3 * XSn) / math.sqrt(6)
LOGUS = np.log(np.abs(US) + 1e-300)
SGN = np.sign(US) + (US == 0)

_wv = {}
def wpv(w):
    if w not in _wv:
        def wpoly(x, w=w):
            tot = 1.0
            for lvl in w:
                if lvl == 1:
                    tot *= (3 * x * x - 3) / math.sqrt(6)
                elif lvl == 2:
                    tot *= x * math.sqrt(6)
                else:
                    tot *= math.sqrt(6)
            return tot
        _wv[w] = np.array([wpoly(x) for x in XSn])
    return _wv[w]


_PM = {}
def Pmat(alpha, KMAX=120):
    key = round(alpha, 10)
    if key in _PM:
        return _PM[key]
    P = np.empty((KMAX + 1, len(XSn)))
    for k in range(KMAX + 1):
        e = -alpha * US * US + k * LOGUS
        m = float(np.max(e))
        if not np.isfinite(m):
            P[k] = 0.0
            continue
        row = np.exp(e - m)
        if k % 2 == 1:
            row = row * SGN
        row = np.nan_to_num(row, nan=0.0, posinf=0.0, neginf=0.0)
        P[k] = row * math.exp(m)
    _PM[key] = P
    return P


_MnC = {}
def Mn(k, alpha, W1, W2):
    P = Pmat(alpha)
    key = (k, W1 if isinstance(W1, tuple) else tuple(W1), W2 if isinstance(W2, tuple) else tuple(W2), round(alpha, 9))
    if key in _MnC:
        return _MnC[key]
    v = float(np.dot(WTSn * wpv(W1) * wpv(W2), P[min(k, P.shape[0]-1)])) if k < P.shape[0] else 0.0
    _MnC[key] = v
    return v


def cfun(pa, t):
    s3 = 0.34101124; s5 = 0.05276111
    r1 = (1 - 3*s3*s3*t*t + 5*s5*s5*t**4) / V
    r2 = (6*s3*s3*t*t + 20*s5*s5*t**4) / V
    r3 = (12*s3*s3*t*t + 80*s5*s5*t**4) / V
    return {(1,0): r3, (2,0): 3*r1*r2, (3,0): r1**3,
            (0,1): -t, (0,2): 3*t*t, (0,3): -(t**3),
            (1,1): -3*t*(r1+r2), (2,1): -3*t*r1*r1, (1,2): 3*t*t*r1}[pa]


_MC = {}
def inner_c(e1, e2, t, NMAX=60, tol=1e-20):
    """Complex inner product ⟨G_e1, conj(G_e2)⟩ at t on the circle."""
    r = (t - (0.34101124**2)*t**3 + (0.05276111**2)*t**5) / V
    s = 1 - r * r
    alpha = (1/s).real
    gamma = (r/s).real
    key = (tuple(map(tuple, e1)), tuple(map(tuple, e2)), round(t.real, 9), round(t.imag, 9))
    if key in _MC:
        return _MC[key]
    tot = 0.0
    for (wx1, wy1, k1, tp1) in e1:
        m1 = mult_terms(k1)
        for (wx2, wy2, k2, tp2) in e2:
            m2 = mult_terms(k2)
            wth = VT0 ** (tp1 + tp2)
            for (au1, av1, ar1, pw1), c1 in m1.items():
                for (au2, av2, ar2, pw2), c2 in m2.items():
                    Z = (r**ar1) * ((r.conjugate())**ar2) * ((s**(-pw1)) * ((s.conjugate())**(-pw2)))
                    cbase = c1 * c2 * Z * wth
                    if abs(cbase) < 1e-30:
                        continue
                    terms = []
                    gpow = 1.0
                    for n in range(NMAX):
                        if n > 0:
                            gpow *= gamma / n
                        mk_x = au1 + au2 + n
                        mk_y = av1 + av2 + n
                        if mk_x >= 120 or mk_y >= 120:
                            terms.append(0.0)
                            continue
                        terms.append(gpow * Mn(mk_x, alpha, wx1, wx2) * Mn(mk_y, alpha, wy1, wy2))
                    S = math.fsum(terms)
                    tot += cbase * S

    # (2pi)^2 x (4pi)^(-2) cancel; final |s|^{-1}:
    res = tot / abs(s)
    # self-inner must be real ≥ 0:
    if e1 == e2 and abs(res.imag) > 1e-9 * max(1e-12, abs(res.real)):
        print(f"WARN inner_c not real for self pair at t={t}: {res}")
    if e1 == e2 and res.real < 0:
        print(f"WARN inner_c NEGATIVE for self pair at t={t}: {res}")
    _MC[key] = res
    return res


def ctable_c(t):
    return {pa: cfun(pa, t) for pa in _IO_LIST}


def S_c(t):
    c = ctable_c(t)
    s1 = 0.0
    for (p1, a1) in _IO_LIST:
        e1 = G_pa_terms(p1, a1)
        for (p2, a2) in _IO_LIST:
            e2 = G_pa_terms(p2, a2)
            s1 += c[(p1, a1)] * (c[(p2, a2)].conjugate()) * inner_c(e1, e2, t)
    s2 = 0.0
    for (p1, a1) in _IO_LIST:
        e1 = G_pa_terms(p1, a1 + 1)
        for (p2, a2) in _IO_LIST:
            e2 = G_pa_terms(p2, a2 + 1)
            s2 += c[(p1, a1)] * (c[(p2, a2)].conjugate()) * inner_c(e1, e2, t)
    return s1, 4 * s2, s1 + 4 * s2


if __name__ == "__main__":
    import time
    print("θ-sweep on the unit circle t = e^{iθ}:", flush=True)
    for th in [0.0, 0.4, 0.8, 1.2, 1.57, 2.0, 2.4, 2.8, 3.14159]:
        t = complex(math.cos(th), math.sin(th))
        t0 = time.time()
        a, b, s = S_c(t)
        r = (t - (0.34101124**2)*t**3 + (0.05276111**2)*t**5)/V
        print(f"θ={th:.3f}: S = {a.real:.4f} + 4·{b.real/4:.4f}·→ {s.real:.4f} (|s(1−r²)|={abs(1-r*r):.3f}) "
              f"[{time.time()-t0:.0f}s]", flush=True)
    print("gate: avg(S) must land in [71.553, 126.804]")
