"""d3h_phi3.py — R3 stage: Φ3 assembly and S(θ) = ||Φ3||² + 4||∂x∂yΦ3||².

Φ3(t;x,y) = Σ_{(p,a)∈ℐ} c_{p,a}(t)·G_{p,a}(t;x,y),  G_{p,a} = [∂_r^p (∂x∂y)^a C_r]_{r=ρ(t)},
ℐ = {(1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1),(2,1),(1,2)},
c-table (paper 1835–1852), ρ_j = D_t^j ρ:
  c_{1,0} = ρ3            c_{0,1} = −t
  c_{2,0} = 3ρ1ρ2         c_{0,2} = 3t²
  c_{3,0} = ρ1³           c_{0,3} = −t³
  c_{1,1} = −3t(ρ1+ρ2)    c_{2,1} = −3tρ1²     c_{1,2} = 3t²ρ1

ρ(t) = c1 t + c3 t³ + c5 t⁵ ⇒
  ρ1 = Dρ = c1 + 3c3t² + 5c5t⁴
  ρ2 = D²ρ = 6c3t² + 20c5t⁴
  ρ3 = D³ρ = 6c3t²·2? — D(t^k) = k t^k: Dρ = c1 + 3c3 t² + 5c5 t⁴;
        D²ρ = 6c3 t² + 20 c5 t⁴; D³ρ = 12 c3 t² + 80 c5 t⁴.  (D keeps t-powers.)

G_{p,a} evaluation = 2π·Σ dxdy entries (weights Wx(x)Wy(y), ϑ^tp) × mult("K"+"r"*p + base-kind) × K_r.
    (r-steps touch ONLY the kind mult; weights/ϑ unchanged.)

‖Φ3‖²_{L²(μ×μ)} = Σ_{i,j} c_i(t)·conj(c_j(t))·⟨G_i, G_j⟩ where
⟨G_i,G_j⟩ = ∫∫ Wi(x)Wj(x)-products... — each G is a SUM of entries; the inner product EXPANDS to
Σ_entry-pairs 2πϑ^{tp_i+tp_j}·(Wx_i Wx_j-products)(x)·(Wy_i Wy_j-products)(y)·[mult_i·mult_j̄]·(K K̄).
K(u,v)·conj(K(u,v)) with r = ρ(t), r̄ = conj: K_r·K_r̄ = 1/(4π²|s|)·exp(−Re[(u²−2ruv+v²)/s·conj-part])
where the exponent −(u²−2ruv+v²)/(2s) − (u²−2r̄uv+v²)/(2s̄) = −Re[(u²+v²) − 2ruv]/s... = real Gaussian:
    −[u²+v²]·Re[1/(2s)]·2·... — work out: E(u,v) = −(u²−2ruv+v²)/(2s) − (u²−2r̄uv+v²)/(2s̄)
    = −(u²+v²)/2·(1/s + 1/s̄) + ruv·(1/s + 1/s̄) = −(u²+v²)·Re(1/s) + uv·2? NO: −(u²+v²)(1/s+1/s̄)/2
    + ruv(1/s+1/s̄) = −(u²+v²)Re(1/s)·2/2... let a = Re(1/s), and (1/s+1/s̄) = 2Re(1/s) = 2a:
    E = −a(u²+v²) + 2a·Re(r)uv·? : ruv(1/s+1/s̄) = 2a·Re(r·uv)? r·uv: uv real ⇒ Re(r·uv) = Re(r)uv:
    E = −a(u²+v²) + 2a·Re(r)·uv = −a[u²+v² − 2Re(r)uv]. Designate α = Re(1/s) > 0 (|r|<1 guarantees),
    ρh = Re(r): K·K̄ = (1/(4π²))·|s|^{−1}·e^{−α(u²+v²−2ρh uv)} — the bivariate normal density with
    correlation ρh up to normalization: exactly α = 1/(1−ρh²)-shaped ⇒ K·K̄ = |K_ρh|²·... indeed
    = K_{ρh}(u,v)² ·( VALUE): K_ρh = (2π)^{−1}(1−ρh²)^{−1/2}e^{−...}: K²: (2π)^{−2}(1−ρh²)^{−1}e^{−2[..]},
    while K·K̄ here has e^{−[u²+v²−2ρh uv]·a} with a = Re(1/s) = Re[(1−r̄²)/|s|²] = (1−Re(r²))/|s|²·...
    Compute exactly in CODE as real arb: α = Re(1/s); then e^{−α(u²+v²)+2αρh uv} and prefactor (4π²|s|)^{−1}.
    SEPARABLE: e^{−αu²}e^{−αv²}e^{2αρh uv}: SERIES e^{2αρhuv} = Σ_n (2αρh)^n u^n v^n/n! ⇒
    ⟨G_i,G_j⟩ = (2π)^{−1}|s|^{−1}·Σ_n β^n/n!·M_n^{|G_i,G_j|-x-weights}·M_n^{y} with β = 2αρh,
    M_n = ∫ W(x)^2-product·u(x)^n·e^{−αu(x)²} μ(dx) — 1-D certified integrals per n.
    |β| bound: |r| ≤ 1 boundary ⇒ |s| ≥ 0.3724 (per-panel verified) ⇒ β ≤ 2·4.6·1 ≈ 9.2 worst-case
    — series n! counter: β^n/n! ≤ 9.2^n/n! — decays after n ≥ 10; N_MAX ~ 60 for 1e-60 tails. FINE.
"""

import math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from core import Prec, PREC, PREC_HI, pi, _dec_fmpq, integrate_certified
from d3h_kernels import (vartheta, rho_eval, psi3f, psi3pf, psi3ppf, psi3pppf,
                         Kf, eval_terms, mult_terms, chain_mult, _T_ONE)
from d3h_dxdy import dxdy_terms, weight_eval, bump_kind, d_weight
from flint import fmpq, arb, acb, ctx

C1 = 1 / (1 + fmpq(34101124, 10**8)**2 + fmpq(5276111, 10**8)**2)
C3 = -(fmpq(34101124, 10**8)**2) * C1
C5 = (fmpq(5276111, 10**8)**2) * C1


def rho_j(j: int, t: acb) -> acb:
    """D^j ρ(t) as acb (exact polynomial in t)."""
    if j == 0:
        return t + C3 * t**3 + C5 * t**5
    if j == 1:
        return C1 + 3 * C3 * t**2 + 5 * C5 * t**4
    if j == 2:
        return 6 * C3 * t**2 + 20 * C5 * t**4
    if j == 3:
        return 12 * C3 * t**2 + 80 * C5 * t**4
    raise ValueError(j)


C_TABLE = {
    (1, 0): lambda t: rho_j(3, t),
    (2, 0): lambda t: 3 * rho_j(1, t) * rho_j(2, t),
    (3, 0): lambda t: rho_j(1, t)**3,
    (0, 1): lambda t: -t,
    (0, 2): lambda t: 3 * t**2,
    (0, 3): lambda t: -(t**3),
    (1, 1): lambda t: -3 * t * (rho_j(1, t) + rho_j(2, t)),
    (2, 1): lambda t: -3 * t * rho_j(1, t)**2,
    (1, 2): lambda t: 3 * t**2 * rho_j(1, t),
}

_IO_LIST = [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2)]


def G_pa_terms(p: int, a: int):
    """Entries for (∂r^p)(∂x∂y)^a C_r. For a=0 (p ≥ 1): ∂r^p C = 2π·∂r^{p−1}K (bridge (47))
    ⇒ single entry, kind "K"+"r"*(p−1), no weights, ϑ⁰. For a ≥ 1: dxdy entries with r^p folded
    into the kind (r-steps don't touch x/y weights)."""
    if a == 0:
        assert p >= 1
        return [((), (), "K" + "r" * (p - 1), 0)]
    out = []
    for (wx, wy, kind, tp) in dxdy_terms(a):
        p2 = kind.count('r'); i2 = kind.count('u'); j2 = kind.count('v')
        assert p2 == 0
        kind_full = "K" + "r" * p + "u" * i2 + "v" * j2
        out.append((wx, wy, kind_full, tp))
    return out



def G_pa_eval(p: int, a: int, x: acb, y: acb, t: acb, vt: arb, bits: int = PREC) -> acb:
    """(∂r^p)(∂x∂y)^a C_r(x,y) at r = ρ(t). NOT including any c-factor."""
    with Prec(bits):
        r = rho_eval(t)
        u = acb(vt) * psi3f(x)
        v = acb(vt) * psi3f(y)
        Kval = Kf(u, v, r, bits)
        tot = acb(0)
        for (wx, wy, kind_full, tp) in G_pa_terms(p, a):
            mval = eval_terms(mult_terms(kind_full), u, v, r, bits)
            wprod = weight_eval(wx, x, bits) * weight_eval(wy, y, bits)
            tot = tot + acb(2) * acb.pi() * (acb(vt) ** tp) * wprod * mval * Kval
        return tot


def Phi3_eval(x: acb, y: acb, t: acb, vt: arb, bits: int = PREC) -> acb:
    """Φ3(t;x,y) = Σ c_{p,a}(t)·G_{p,a}(t;x,y)."""
    with Prec(bits):
        tot = acb(0)
        for (p, a) in _IO_LIST:
            c = C_TABLE[(p, a)](t)
            if c == 0:
                continue
            tot = tot + c * G_pa_eval(p, a, x, y, t, vt, bits)
        return tot


def dxdy_Phi3_eval(x: acb, y: acb, t: acb, vt: arb, bits: int = PREC) -> acb:
    """(∂x∂y)Φ3 = Σ c_{p,a}·G_{p,a+1}."""
    with Prec(bits):
        tot = acb(0)
        for (p, a) in _IO_LIST:
            c = C_TABLE[(p, a)](t)
            if c == 0:
                continue
            tot = tot + c * G_pa_eval(p, a + 1, x, y, t, vt, bits)
        return tot
