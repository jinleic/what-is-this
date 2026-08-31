"""d3h_inner.py — float64 stage for ⟨G_i, G_j⟩ and S(θ) = ‖Φ3‖² + 4‖∂x∂yΦ3‖² (structural probe).

⟨G_i, G_j⟩ factorization (real t first; complex both-sides handled by the caller):
  G entries: (wx, wy, kind, tp): value = 2π ϑ^tp Wx(x) Wy(y) m_[kind](u,v; r) K_r(u,v).
  inner = ∫∫ Gi·Gj μ×μ = Σ_(e1,e2) 4π² ϑ^{tp1+tp2} · Pn(x-combo) · Pn(y-combo) · CORE(m1,m2)
  where CORE absorbs m1·m2·|K|²: with α = Re(1/s), ρh = Re(r) (REAL case: 1/s, r):
   |K|² = (4π²|s|)^{-1} e^{−α(u²+v²−2ρh uv)} — write β = 2αρh, expand e^{βuv} = Σ β^n u^n v^n/n!:
   CORE(m1,m2) = |s|^{-1} Σ_n (β^n/n!) · Σ_(mon pairs) c1c2 s^{−(pw1+pw2)} r^{ar1+ar2}
                 · Mn(au1+au2+n) Mn(av1+av2+n).
  Mn(k) = ∫ W-combo(x)·u(x)^k·e^{−αu(x)²} dμ(x) — hermite nodes, LOG-DOMAIN for u^k (overflow guard),
  nodes |x| > XMAX dropped (weights e^{−x²/2} < e^{−98} ⇒ negligible).
S(θ) structurally: S = Σ_ij c_i c̄_j [⟨Gi,Gj⟩] + 4·Σ_ij c_i c̄_j [⟨Gi',Gj'⟩], Gi' = G_{p,a+1}.
"""
import math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
from numpy.polynomial.hermite import hermgauss
from d3h_kernels import mult_terms
from d3h_phi3 import G_pa_terms, _IO_LIST

V = 1 + 0.34101124**2 + 0.05276111**2
VT0 = 0.136419125 / math.sqrt(V)

_yg, _wg = hermgauss(180)
XS = math.sqrt(2) * _yg
WTS = _wg / math.sqrt(math.pi)
XMAX = 12.0
_mask = np.abs(XS) <= XMAX
XSn = XS[_mask]
WTSn = WTS[_mask]
US = VT0 * (XSn**3 - 3 * XSn) / math.sqrt(6)
LOGUS = np.log(np.abs(US) + 1e-300)
SGNUS = np.sign(US) + (US == 0)



_PMAT = {}
def _pmatrix(alpha, KMAX=110):
    """P[k, node] = sign(u)·|u|^k·e^{−αu²} for k = 0..KMAX (sign tracked separately for odd k)."""
    key = round(alpha, 12)
    if key in _PMAT:
        return _PMAT[key]
    ex = -alpha * US * US
    # log|u|^k − αu² computed via iterated multiply (stable here because the row product is dominated
    # by moderate nodes; large-|u| nodes decay like e^{−αu²} faster than u^k grows for k ≤ 110 since
    # |u| ≤ 11 ⇒ u^110 ≤ 1e114 but e^{−α·121} = e^{−259} kills it — combined ok in float via max-exp
    # LSE per row).
    lgabs = np.log(np.abs(US) + 1e-300)


    with np.errstate(over='ignore', invalid='ignore'):
        P = np.empty((KMAX + 1, len(US)))
        for k in range(KMAX + 1):
            e = ex + k * lgabs
            m = np.max(e)
            if not np.isfinite(m):
                P[k] = 0.0
                continue
            ev = np.exp(e - m)
            row = ev * (SGNUS ** k if k % 2 == 1 else 1.0)
            row[~np.isfinite(row)] = 0.0
            P[k] = row * math.exp(m)
    _PMAT[key] = P
    return P


def Mn_fast(k, alpha, Wx1, Wx2, logu_scale=None):
    """∫ W1(x)W2(x)·u(x)^k·e^{−αu²} dμ (k ≥ 0) via precomputed P-matrix."""
    P = _pmatrix(alpha)
    wprod = wpoly_vec(Wx1) * wpoly_vec(Wx2)
    return float(np.dot(WTSn * wprod, P[k]))


def wpoly(w, x):
    tot = 1.0
    for lvl in w:
        if lvl == 1:
            tot *= (3 * x * x - 3) / math.sqrt(6)
        elif lvl == 2:
            tot *= x * math.sqrt(6)
        else:
            tot *= math.sqrt(6)
    return tot



_wp_cache = {}
def wpoly_vec(w):
    key = w
    if key in _wp_cache:
        return _wp_cache[key]
    vals = np.array([wpoly(w, x) for x in XSn])
    _wp_cache[w] = vals
    return vals



_MC = {}
def inner_terms(e1, e2, r, NMAX=64, tol=1e-22):
    """⟨G_e1, G_e2⟩ for REAL r (t real), entries lists."""
    s = 1 - r * r
    alpha = 1.0 / s
    beta = 2 * alpha * r
    key = (tuple(map(tuple, e1)), tuple(map(tuple, e2)), round(r, 12))
    if key in _MC:
        return _MC[key]
    tot = 0.0
    for (wx1, wy1, k1, tp1) in e1:
        m1 = mult_terms(k1)
        for (wx2, wy2, k2, tp2) in e2:
            m2 = mult_terms(k2)
            wtheta = VT0 ** (tp1 + tp2)
            for (au1, av1, ar1, pw1), c1 in m1.items():
                for (au2, av2, ar2, pw2), c2 in m2.items():
                    cbase = c1 * c2 * (s ** (-(pw1 + pw2))) * (r ** (ar1 + ar2)) * wtheta
                    terms = []
                    bpow = 1.0
                    for n in range(NMAX):
                        if n > 0:
                            bpow *= beta / n
                        t1 = Mn_fast(au1 + au2 + n, alpha, wx1, wx2)
                        t2 = Mn_fast(av1 + av2 + n, alpha, wy1, wy2)
                        terms.append(bpow * t1 * t2)
                    S = math.fsum(terms)
                    tot += cbase * S
    res = tot / s   # (2pi)^2 G-prefactor over (4pi^2) of |K|^2 cancel; s^{-1} remains
    _MC[key] = res
    return res


def ctable_float(t):
    r1 = (1 - 3 * 0.34101124**2 * t * t + 5 * 0.05276111**2 * t**4) / V
    r2 = (6 * 0.34101124**2 * t * t + 20 * 0.05276111**2 * t**4) / V
    r3 = (12 * 0.34101124**2 * t * t + 80 * 0.05276111**2 * t**4) / V
    return {
        (1, 0): r3, (2, 0): 3 * r1 * r2, (3, 0): r1**3,
        (0, 1): -t, (0, 2): 3 * t * t, (0, 3): -t**3,
        (1, 1): -3 * t * (r1 + r2), (2, 1): -3 * t * r1 * r1, (1, 2): 3 * t * t * r1,
    }


def S_float(t):
    c = ctable_float(t)
    s1 = 0.0
    for (p1, a1) in _IO_LIST:
        e1 = G_pa_terms(p1, a1)
        for (p2, a2) in _IO_LIST:
            e2 = G_pa_terms(p2, a2)
            s1 += c[(p1, a1)] * c[(p2, a2)] * inner_terms(e1, e2, r=t)
    s2 = 0.0
    for (p1, a1) in _IO_LIST:
        e1 = G_pa_terms(p1, a1 + 1)
        for (p2, a2) in _IO_LIST:
            e2 = G_pa_terms(p2, a2 + 1)
            s2 += c[(p1, a1)] * c[(p2, a2)] * inner_terms(e1, e2, r=t)
    return s1, 4 * s2, s1 + 4 * s2


if __name__ == "__main__":
    import time
    for t in [0.0, 0.25, 0.5, 0.75, 0.9]:
        t0 = time.time()
        a, b, s = S_float(t)
        print(f"S({t}) = {a:.6f} + 4·{b/4:.6f} = {s:.6f}   ({time.time()-t0:.0f}s)")
    print("budget (π²/6)·126.804 =", math.pi**2/6*126.80385221)
