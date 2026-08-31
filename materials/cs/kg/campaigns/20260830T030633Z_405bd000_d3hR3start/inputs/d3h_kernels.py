"""d3h_kernels.py — R3: independent kernel-side derivation of ||D^3 H||_{L^2(T)}.

T = t·∂t Euler operator. C_r(x,y) = (π/2)Σ_{a≥0} r^a q_a(u)q_a(v) (49) with q_a orthonormal-form
q_a(u) = 2φ(u)ψ_{a−1}(−u)/√a; ψ_n = He_n/√(n!); ψ3(x) = (x³−3x)/√6, ψ3'(x) = (3x²−3)/√6,
ψ3''(x) = √6 x, ψ3'''(x) = √6. u = ϑψ3(x), v = ϑψ3(y).

Mehler-kernel identities (paper Lemma C.1):  ∂_r C_r = 2π K_r(u,v),
    ∂_x∂_y C_r = 2π ϑ² ψ3'(x)ψ3'(y) K_r(u,v),
K_r(u,v) = (2π)^{−1}(1−r²)^{−1/2} e^{−(u²−2ruv+v²)/(2(1−r²))}.

Recursive differentiation (∂x∂y acts on the product):
  ∂x (w(x) K_r) = w'(x) K + w(x) u' ∂uK,  ∂uK = L·K,  L = (−u+rv)/(1−r²),
  ∂u∂v K = M'·K + L·M·K,  M = (−v+ru)/(1−r²),  ∂vM = r/(1−r²).
So every G_{p,a} is a FINITE sum Σ_k α_k·W_k(x)·W'_k(y)·F_k(K, ∂uK, ∂vK, ∂u∂vK, ∂r^j K...).

Implementation strategy: represent (∂x∂y)^a C_r as a linear combination
  Σ_k c_k(t) · X_k(x) · Y_k(y) · 𝒦_k(u,v; r)
where 𝒦_k ∈ {K, ∂uK, ∂vK, ∂u∂vK} and X,Y are products of ψ3-derivatives (polynomials). Then
∂_r^p hits the 𝒦_k factors via ∂r K = ℓ·K (ℓ = (log K)_r), ∂r∂uK = (∂uℓ + ℓL)K etc:
  ∂rK = ℓK;            ∂r∂uK = (∂uℓ)K + ℓ·∂uK = K·(∂uℓ + ℓL)
  ∂r∂vK = K·(∂vℓ + ℓM)
  ∂r∂u∂vK = K·(∂u∂vℓ + ∂uℓ·M + ∂vℓ·L + ℓ(r/(1−r²) + L·M))
with ℓ = r/s − Q'/2 (Q = n0/s, n0 = u²−2ruv+v², s = 1−r²), and
  ∂uℓ = (−2u+2rv)/s = −2L,          ∂vℓ = −2M,
  ∂u∂vℓ = 2r/s.
(These last follow from ∂uℓ = ∂u[(r)/s] + ∂u[−Q'/2]; Q' = ∂r Q = ∂r[n0/s] = (−2uv)/s + 2rn0/s²;
 ∂uQ' = −2v/s + 2r(2u−2rv)/s² = 2(−v+ru)/s = 2M... check sign: ∂u n0 = 2u−2rv = 2(u−rv) = −2(−u+rv)
 ⇒ ∂uℓ = ∂u[r/s] − (1/2)∂uQ' = 0 − (1/2)(−2v/s + 2r(2u−2rv)/s²). Let me recompute in code directly via
 symbolic-differentiation test rather than trusting algebra here.)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from core import Prec, PREC, PREC_HI, pi, _dec_fmpq
from flint import fmpq, arb, acb

ETA = _dec_fmpq("0.136419125")
S3 = _dec_fmpq("0.34101124")
S5 = _dec_fmpq("0.05276111")
V = 1 + S3 * S3 + S5 * S5
C1 = 1 / V
C3 = -(S3 * S3) / V
C5 = (S5 * S5) / V
THETA_VAR = None  # vartheta arb computed lazily


def vartheta(bits: int = PREC_HI) -> arb:
    global THETA_VAR
    if THETA_VAR is None:
        with Prec(bits):
            THETA_VAR = arb(ETA) / arb(V).sqrt()
    return THETA_VAR


def rho_eval(t: acb) -> acb:
    return t + C3 * t**3 + C5 * t**5


def psi3f(x):
    return (x**3 - 3 * x) / arb(6).sqrt()


def psi3pf(x):
    return (3 * x**2 - 3) / arb(6).sqrt()


def psi3ppf(x):
    return x * arb(6).sqrt()


def psi3pppf(x):
    return arb(6).sqrt()


# ---------------------------------------------------------------------------
# K_r and its partials, all as (multiplier)·K with multipliers built by ball arithmetic.
# ℓ ≡ (∂r K)/K,  L ≡ (∂u K)/K,  M ≡ (∂v K)/K  — computed EXACTLY (as balls) below.
# ---------------------------------------------------------------------------

def _s(r):
    return 1 - r * r


def Lf(u, v, r):
    """L = ∂uK/K."""
    return (-u + r * v) / _s(r)


def Mf(u, v, r):
    """M = ∂vK/K."""
    return (-v + r * u) / _s(r)


def Kf(u, v, r, bits=PREC):
    with Prec(bits):
        s = _s(r)
        A = 1 / (2 * acb.pi() * s.sqrt())
        Q = (u * u - 2 * r * u * v + v * v) / s
        return A * (-Q / 2).exp()


def ellf(u, v, r, bits=PREC):
    """ℓ = (d/dr K)/K = r/s + uv/s − r·n0/s², s = 1−r², n0 = u²−2ruv+v². [FD-validated]"""
    with Prec(bits):
        s = _s(r)
        n0 = u * u - 2 * r * u * v + v * v
        return r / s + (u * v) / s - r * n0 / (s * s)


def L_puf(u, v, r):
    """∂u L = −1/s."""
    return -1 / _s(r)


def M_vf(u, v, r):
    return -1 / _s(r)


def du_ellf(u, v, r, bits=PREC):
    """∂u ℓ = −Q'_u/2 with Q'_u = (−2v)/s + 2r(2u−2rv)/s². [FD-validated]"""
    with Prec(bits):
        s = _s(r)
        dQpu = (-2 * v) / s + 2 * r * (2 * u - 2 * r * v) / (s * s)
        return -dQpu / 2


def dv_ellf(u, v, r, bits=PREC):
    with Prec(bits):
        s = _s(r)
        dQpv = (-2 * u) / s + 2 * r * (2 * v - 2 * r * u) / (s * s)
        return -dQpv / 2


def duv_ellf(u, v, r, bits=PREC):
    """∂u∂v ℓ [FD-validated]."""
    with Prec(bits):
        s = _s(r)
        return 1 / s + 2 * r * r / (s * s)


def ell1f(u, v, r, bits=PREC):
    """ℓ' = d/dr ℓ — term-dict {(0,0,0,1):1,(2,0,0,2):2,(1,1,0,2):6,(0,0,1,2):−1,(2,0,1,3):−4,(3,1,0,3):8},
    FD-VALIDATED."""
    with Prec(bits):
        s = _s(r); uv = u * v; qq = u * u + v * v
        return (1.0 * s ** -1
                + 2.0 * r**2 * s ** -2
                + 6.0 * r * uv * s ** -2
                - 1.0 * qq * s ** -2
                - 4.0 * r**2 * qq * s ** -3
                + 8.0 * r**3 * uv * s ** -3)



def ell2f(u, v, r, bits=PREC):
    """ℓ''. [FD-VALIDATED]"""
    with Prec(bits):
        s = _s(r); uv = u * v; qq = u * u + v * v
        return (6.0 * r * s ** -2
                + 8.0 * r**3 * s ** -3
                + 6.0 * uv * s ** -2
                + 48.0 * r**2 * uv * s ** -3
                - 12.0 * r * qq * s ** -3
                - 24.0 * r**3 * qq * s ** -4
                + 48.0 * r**4 * uv * s ** -4)


def ell3f(u, v, r, bits=PREC):
    """ℓ'''. [FD-VALIDATED]"""
    with Prec(bits):
        s = _s(r); uv = u * v; qq = u * u + v * v
        return (6.0 * s ** -2
                + 48.0 * r**2 * s ** -3
                + 48.0 * r**4 * s ** -4
                + 120.0 * r * uv * s ** -3
                + 480.0 * r**3 * uv * s ** -4
                - 12.0 * qq * s ** -3
                - 144.0 * r**2 * qq * s ** -4
                - 192.0 * r**4 * qq * s ** -5
                + 384.0 * r**5 * uv * s ** -5)

# ---------------------------------------------------------------------------
# (∂x∂y)^a C_r and ∂_r^p thereof as FINITE WEIGHT SUMS over K-partials.
# Representation: a "term" = (s_x, s_y, kind) with s_x, s_y := products of ψ3-derivative labels,
# kind ∈ {"K","Ku","Kv","Kuv","Kr","Kru","Krv","Kruv","Kr2","Kr2u","Kr2v","Kr2uv","Kr3",...}.
# We only need p ≤ 3, a ≤ 3. ∂r acts on any 𝒦 by the chain through ell/L/M multipliers.
# ---------------------------------------------------------------------------

_DR = {"K": "Kr", "Kr": "Kr2", "Kr2": "Kr3", "Ku": "Kru", "Kru": "Kr2u",
       "Kr2u": "Kr3u", "Kv": "Krv", "Krv": "Kr2v", "Kr2v": "Kr3v",
       "Kuv": "Kruv", "Kruv": "Kr2uv", "Kr2uv": "Kr3uv"}


def _kind_value(kind, u, v, r, bits=PREC):
    """Return the multiplier (ball) such that the operator = multiplier × K."""
    with Prec(bits):
        s = _s(r)
        L, M, ell = Lf(u, v, r), Mf(u, v, r), ellf(u, v, r, bits)
        du_ell, dv_ell, duv_ell = (du_ellf(u, v, r, bits), dv_ellf(u, v, r, bits),
                                   duv_ellf(u, v, r, bits))
        table = {
            "K":      acb(1),
            "Ku":     L,
            "Kv":     M,
            "Kuv":    M_vf(u, v, r) / 1 * 0,  # placeholder replaced below
        }
        # Kuv = ∂u∂vK / K = r/s + L*M
        table["Kuv"] = r / s + L * M
        ell1 = ell1f(u, v, r, bits)
        ell2 = ell2f(u, v, r, bits)
        ell3 = ell3f(u, v, r, bits)
        # ---- mults by RECURSIVE COMPOSITION (all FD/series-validated) ----
        # NOTATION: mult(∂r^p ∂u^i ∂v^j K)/K built from base mults L, M, g0 = r/s + L·M and
        # base functions ℓ, ℓ', ℓ''; recursion: next_r(m) = (∂r m) + m·ℓ where ∂r of each ingredient:
        #   ∂r L = v/s + 2rL/s ; ∂r²L = 2rv/s³ + 2L/s + 2r(∂rL)/s − 4r²L/s²
        #   ∂r M symmetric ; ∂r ℓ = ℓ' ; ∂r ℓ' = ℓ'' ; ∂r g0 = d/dr(r/s) + (∂rL)M + L(∂rM)
        drL = v / s + 2 * r * L / s
        dr2L = (4 * r * v + 2 * (-u + r * v)) / (s * s) + 8 * r * r * (-u + r * v) / (s**3)
        drM = u / s + 2 * r * M / s
        dr2M = (4 * r * u + 2 * (-v + r * u)) / (s * s) + 8 * r * r * (-v + r * u) / (s**3)
        g0 = r / s + L * M
        dr_g0 = (1 + r * r) / (s * s) + drL * M + L * drM
        mKu  = L                                   # ∂uK/K
        mKru = drL + L * ell                       # ∂r∂uK/K  = next_r(L)
        mKr2u = dr2L + 2 * drL * ell + L * (ell1 + ell * ell)      # next_r(mKru)
        mKv  = M
        mKrv = drM + M * ell
        mKr2v = dr2M + 2 * drM * ell + M * (ell1 + ell * ell)
        mKruv = dr_g0 + g0 * ell                   # next_r(g0)
        dr2_g0 = 2 * r / (s * s) + (1 + r * r) * 4 * r / (s**3) + dr2L * M + 2 * drL * drM + L * dr2M
        mKr2uv = dr2_g0 + dr_g0 * ell + g0 * ell1 + mKruv * ell
        table["Kru"] = mKru
        table["Krv"] = mKrv
        table["Kruv"] = mKruv
        table["Kr2"] = ell1 + ell * ell
        table["Kr2u"] = mKr2u
        table["Kr2v"] = mKr2v
        table["Kr3"] = ell2 + 3 * ell * ell1 + ell**3
        table["Kr"] = ell
        table["Kr2uv"] = mKr2uv


        if kind in table and table[kind] is not None:
            return table[kind]
        raise NotImplementedError(f"kind {kind}")




# ---------------------------------------------------------------------------
# SYMBOLIC G_{p,a} BUILDER (a ≤ 3). Represent (∂x∂y)^a C_r as a list of
#   (scalar, wx_fn, wy_fn, k-chain) where k-chain = tuple like ("K",) / ("K","Ku") etc.
# Denoting 𝒦 = K_r, ∂x𝒦 = u'·L𝒦, ∂y𝒦 = v'·M𝒦, and ∂u(L𝒦) = (∂uL)𝒦 + L(L𝒦) ⇒ the u/v operator
# MULTIPLIERS compose alphabetically: a "u-step" on mult m gives (∂u m) + m·L, a "v-step": (∂v m)+m·M.
# Chain applied in the order the derivatives occurred (commutative — partials commute — so any fixed
# order works as long as applied consistently). For a chain like (u,uv): start mult 1, apply steps.
# Steps:
#   u-step:  m ↦ du(m) + m·L        with du(1)=0, du(L)=∂uL=−1/s, du(M)=r/s, du(g0)= du(g0):
#            g0 = r/s + L·M: du(g0) = duL·M + L·duM = −M/s + L·r/s
#   v-step:  m ↦ dv(m) + m·M        with dv(L) = r/s, dv(M) = −1/s, dv(g0) = −L/s + r·M/s·? — derive:
#            dv(L) = ∂v[(−u+rv)/s] = r/s ✓; dv(M) = −1/s; dv(g0) = 0 + dv(L)·M + L·dv(M) = rM/s − L/s
#   and ∂u(ℓ), ∂v(ℓ), ∂u∂v(ℓ): closed forms (du_ellf etc.), ∂u of ℓ', ∂v of ℓ'' needed for r-chains:
#   r-step:  m ↦ dr(m) + m·ℓ         dr(ℓ) = ℓ', dr(L) = drL, dr(M) = drM, dr(g0) = (1+r²)/s² + drL·M + L·drM
# The mult OBJECT is an acb-evaluable function of (u,v,r) via a small term-DAG; we represent it as a
# PYTHON closure over (u,v,r) — using only ball arithmetic inside. Composition is exact.
# ---------------------------------------------------------------------------
from typing import Callable, List, Tuple

MultFcn = Callable[[acb, acb, acb], acb]


def _div_s(x: acb, s: acb) -> acb:
    return x / s


def make_mult_monomial() -> MultFcn:
    return lambda u, v, r: acb(1)


def step_u(m: MultFcn) -> MultFcn:
    """u-step on mult m: m ↦ ∂u m + m·L, where ∂u hits L,M,g0,ℓ-chain constants."""
    def f(u, v, r):
        bits = ctx.prec
        with Prec(bits):
            s = _s(r)
            L = Lf(u, v, r)
            M = Mf(u, v, r)
            # du of base quantities known analytically:
            duL = acb(-1) / s
            duM = r / s
            # du applied to m — m is a general mult; we need m's structural derivative.
            # Represent m as (kind-specific closed form) — for the mults that occur (built only
            # from 1, L, M, g0, ℓ by u/v/r-steps), we dispatch on m.mult_kind metadata.
            raise NotImplementedError
    return f


# STOP — the closure-composition route needs structural metadata. Cleaner: represent mults as
# LINE COMBINATIONS over a FIXED BASIS of "operator eigenfunctions". Observed: every mult is a
# POLYNOMIAL in (u, v) of degree ≤ (derivatives) times (rational in r) — because each step of the
# chain on (L, M, ℓ, g0) yields terms (1/s^k)·u^a·v^b·(polynomials in r). The exact mult for
# ∂r^p∂u^i∂v^jK = [poly(u,v) degree ≤ 2p+i+j−(i+j)…]/·K. Instead of symbolic diff, EXPLOIT:
#   mult(∂u^i ∂v^j ∂r^p K) — FD-VALIDATED FORMULAS EXIST for exactly the 13 kinds needed (all n ≤ 3
#   totals): K, Ku, Kv, Kuv, Kr, Kru, Krv, Kruv, Kr2, Kr2u, Kr2v, Kr2uv, Kr3. These ARE the needed
#   set for ℐ = {(p,a): p ≤ 3, a ≤ 3} because (∂x∂y)^a C_r's K-side kinds are precisely the u/v-steps
#   of count i+j ≤ 2(a−1)+… — ENUMERATE what's needed:
#     a=1: K ; a=2: K, Ku, Kv, Kuv ; a=3: K, Ku2, Kv2, Ku2v, Kuv2, Kuv (from ∂x∂y of a=2 terms)
#   with p = 0..3 r-steps: TOTAL KINDS = {K, Ku, Kv, Kuv} × r-steps 0..3 PLUS {Ku2, Kv2, Ku2v, Kuv2} × 0..3.
#   EXISTING kinds cover: K/Ku/Kv/Kuv × {0..3 r} = 13 minus Kr3u/Kr3v/Kr3uv = 10 present; need
#   additionally: Ku2 (u²), Kv2, Ku2v, Kuv2 for a=3 at r-steps 0..3 — 12 more kinds (16 total).
# Build them by the SAME validated pattern (tree of letters u/v printed in a canonical order with
# falling-factorial bookkeeping)... — decision: implement u2/v2/u2v/uv2 mults by recursive composition
# with EXPLICIT dict-based term algebra (the approach that worked for ell chain). See mult_algebra.


# ---------------------------------------------------------------------------
# mult_algebra: dict-based term algebra for u/v/r chains (mirrors the validated ell/d… approach).
# A mult is TERMS: dict {(a_u, a_v, a_r, pw): coeff} meaning coeff·u^{a_u} v^{a_v} r^{a_r} s^{−pw}·(1)
# — factor (1/√s) prefactors A(r) are tracked SEPARATELY (they multiply K... careful: mults are
# relative to K, so mult = Σ terms/(s^{pw})). u/v appear ONLY polynomially in L/M/ℓ.
# All derivatives in u/v/r act term-by-term with the rules:
#   ∂u: a_u += 1 → coeff·a_u ;  EXCEPT the (1/s^{pw}) factor is u-free (u,v enter only polynomially
#       in L,M,ℓ,g0) ✓.
#   ∂v: a_v += 1 similarly.
#   ∂r: two effects: r^{a_r} → a_r r^{a_r−1} (if a_r>0) AND s^{−pw} → +2r·pw·s^{−pw−1}.
# Basis-closure: is the algebra closed under these? ∂r of u^a v^b r^c s^{−p}: YES (r-powers change,
# s-power grows) — CLOSED. Closure needs uv-DEGREES bounded: u-derivatives INCREASE a_u — but we only
# ever APPLY at most 3 total steps, so degrees ≤ 3 ✓. REPRESENTATION COMPLETE AND CLOSED.
# Baseline mult values needed to start: 1 (for K) and per the u/v/r chain from the paper:
#   Ku = L: dict {(0,0,0,1): −u-term?} — NO: L = (−u + rv)/s: terms {(1,0,0,1): −1, (0,1,1,1): +1} ✓
#   Kv = M: {(0,1,0,1): −1, (1,1,1,1): +1} ✓
#   Kr = ℓ: ℓ = r s^{−1} + uv s^{−1} − r·n0 s^{−2} with n0 = qq − 2ruv (qq = u²+v²):
#        = r s^{−1} + uv s^{−1} − r qq s^{−2} + 2r² uv s^{−2}
#        terms: {(0,0,1,1): 1, (1,1,0,1): 1, (0,0,1,2): −1·qq — PROBLEM: qq = u²+v² is a SUM of two
#        monomials — represent as TWO terms: {(0,0,1,2): −1 (u²), (2,0,... } — i.e. expand: ℓ =
#        r·s^{−1} + uv·s^{−1} − u²·r·s^{−2} − v²·r·s^{−2} + 2r²uv·s^{−2}
#        terms: {(0,0,1,1):1, (1,1,0,1):1, (2,0,1,2):−1, (0,2,1,2):−1, (1,1,2,2):2} ✓ CLOSED.
#   Kuv = ∂u∂vK/K = ∂u(L)·? = du(L) + L·M = −1/s + L·M:
#        L·M = (−u+rv)(−v+ru) s^{−2} = (uv − r u² − r v² + r²uv) s^{−2}: terms {(1,1,0,2):1,
#        (2,0,1,2):−1, (0,2,1,2):−1, (1,1,2,2):1} + {du(L) = −1/s}: {(0,0,0,1):−1} ✓.
# IMPLEMENTATION: functions u_step(terms), v_step(terms), r_step(terms) returning new term dicts;
# evaluate(terms, u, v, r) via ball arithmetic. Build the 25 needed kinds programmatically, then
# FD-VALIDATE EVERY ONE automatically (single harness). This replaces _kind_value's hand table.


def _norm(terms: dict) -> dict:
    return {k: c for k, c in terms.items() if c != 0}


def u_step(terms: dict) -> dict:
    out = {}
    for (au, av, ar, pw), c in terms.items():
        if au > 0:
            k = (au - 1, av, ar, pw)
            out[k] = out.get(k, 0) + c * au
    return _norm(out)


def v_step(terms: dict) -> dict:
    out = {}
    for (au, av, ar, pw), c in terms.items():
        if av > 0:
            k = (au, av - 1, ar, pw)
            out[k] = out.get(k, 0) + c * av
    return _norm(out)


def r_step(terms: dict) -> dict:
    out = {}
    for (au, av, ar, pw), c in terms.items():
        if ar > 0:
            k = (au, av, ar - 1, pw)
            out[k] = out.get(k, 0) + c * ar
        if pw > 0:
            k = (au, av, ar + 1, pw + 1)
            out[k] = out.get(k, 0) + c * 2 * pw
    return _norm(out)


_T_ONE = {(0, 0, 0, 0): 1}
_T_L = {(1, 0, 0, 1): -1, (0, 1, 1, 1): 1}
_T_M = {(0, 1, 0, 1): -1, (1, 0, 1, 1): 1}
_T_ELL = {(0, 0, 1, 1): 1, (1, 1, 0, 1): 1, (2, 0, 1, 2): -1, (0, 2, 1, 2): -1, (1, 1, 2, 2): 2}


def eval_terms(terms: dict, u: acb, v: acb, r: acb, bits: int = PREC) -> acb:
    with Prec(bits):
        s = _s(r)
        tot = acb(0)
        for (au, av, ar, pw), c in terms.items():
            tot = tot + acb(c) * (u ** au) * (v ** av) * (r ** ar) * (s ** (-pw) if pw else acb(1))
        return tot


def kind_terms(chain: Tuple[str, ...]) -> dict:
    """mult for operator = K·(chain) where chain letters: 'r','u','v' (canonical COMMUTATIVE order:
    since partials commute we canonicalize r's first, then u's, then v's — as computed mults of
    ∂r^p ∂u^i ∂v^jK: start from... CAREFUL: the mult of ∂u^i∂v^jK = u^i-then-v^j composition:
    u^i(K): 1→L→(du(L)+L·L)→... — u-steps on a mult term-dict: du(m) is u_step(m) PLUS handle the
    u-dependence INSIDE L, M, ℓ via the PRODUCT RULE: mult composition m_new = du(m) + m·L —
    requires du(m) = u_step(m) ONLY IF m's term dict already contains ALL u-monomials — TRUE by
    construction (L, M, ℓ, g0, and products thereof are polynomials in (u,v) over s^{−p}) —
    so du(m) = u_step(m) ✓ and the composition is m ↦ u_step(m) + m·L_terms (product of term dicts).
    Product of dicts: convolvable — implement _tmul."""
    raise NotImplementedError


def _tmul(A: dict, B: dict) -> dict:
    out = {}
    for (au, av, ar, pw), c in A.items():
        for (bu, bv, br, bpw), d in B.items():
            k = (au + bu, av + bv, ar + br, pw + bpw)
            out[k] = out.get(k, 0) + c * d
    return _norm(out)


def _tadd(A: dict, B: dict) -> dict:
    out = dict(A)
    for k, c in B.items():
        out[k] = out.get(k, 0) + c
    return _norm(out)


_T_G0 = _tadd({(0, 0, 1, 1): 1}, _tmul(_T_L, _T_M))   # Kuv mult = r/s + L·M


def chain_mult(chain: Tuple[str, ...]) -> dict:
    """Terms for mult(∂_chain K)/K where chain letters ∈ {'u','v','r'} — commutative partials, so
    fix the order: all r's first, then u's, then v's for CANONICALITY... but the recursion for r-steps
    on a mult containing u,v terms needs ∂r of u-polynomials — r_step handles (u,v-free in s) ✓.
    Composition rule for a r-step: m ↦ r_step(m) + m·ℓ (product rule: ∂r(m·K) = (∂r m)K + m·ℓK).
    For a u-step: m ↦ u_step(m) + m·L. For v-step: m ↦ v_step(m) + m·M.
    ANY order of steps gives the same final mult (partials commute) — process the chain AS GIVEN."""
    m = dict(_T_ONE)
    for ch in chain:
        if ch == 'u':
            m = _tadd(u_step(m), _tmul(m, _T_L))
        elif ch == 'v':
            m = _tadd(v_step(m), _tmul(m, _T_M))
        elif ch == 'r':
            m = _tadd(r_step(m), _tmul(m, _T_ELL))
        else:
            raise ValueError(ch)
    return m


# The 25 kinds needed for ℐ (p ≤ 3, a ≤ 3): chains sorted (commutative ⇒ canonical sort rs|us|vs):
_KIND_CHAINS = {}
for _p in range(0, 4):
    for _i in range(0, 4):
        for _j in range(0, 4):
            if _p + _i + _j == 0:
                continue
            _name = "K" + "r" * _p + "u" * _i + "v" * _j
            _KIND_CHAINS[_name] = ("r",) * _p + ("u",) * _i + ("v",) * _j


def mult_terms(kind: str) -> dict:
    if kind == "K":
        return {(0, 0, 0, 0): 1}
    chain = _KIND_CHAINS.get(kind)
    if chain is None:
        raise KeyError(f"unknown kind {kind}")
    return chain_mult(chain)


# ---------------------------------------------------------------------------
# (dx dy)^a C_r ASSEMBLY (paper Lemma C.1 generalized). C_r(x,y) = (π/2)Σ r^b q_b(u)q_b(v) with
# u = ϑψ3(x). Differentiating in x hits u' = ϑψ3'(x); each ∂x gives weight-update W ↦ W'·? via the
# PRODUCT RULE on (Wx(x)𝒦(u(x),v)): ∂x[Wx(x Wy(y))𝒦] = Wx'·𝒦 + Wx·u'·(𝒦-derivatives).
# We track SYMBOLIC term lists: entries (Wx_fn, Wy_fn, kind) where kind = "K" + "r"*p + "u"*i + "v"*j,
# Wx_fn/Wy_fn are ψ3-derivative-weight FUNCTIONS of x (resp. y). ∂x∂y of an entry yields up to 2·2
# sub-entries — build a=1,2,3 programmatically:
#   upgrade "u": kind K → Ku (i += 1) and weight gets u' = ϑψ3'(x); an ∂x of an existing u-weight
#   produces ψ3''-weights WITHOUT kind change... precisely: ∂x[Wx·𝒦] where 𝒦 has u-degree i:
#     = Wx'·𝒦 + Wx·(ϑψ3')·𝒦_{u+1}.
# Raise power ϑ FACTORS: each u' adds a ϑ. We bundle ϑ^{power} into the scalar coefficient (evaluated
# at call time). Term: (Wx_fn, Wy_fn, kind) with GLOBAL ϑ-power tracked separately per (a) — collect
# ϑ^{2a}·(counts of u'/v' factors): actually each kind letter u/v contributes exactly one u'/v' factor;
# ∂x of WxFn contributes NO ϑ. So ϑ-power = (i + j)/1 where i = #u-letters, j = #v-letters? NO:
# 𝒦 = ∂u^i∂v^j K is a mult·K where the u,v in the mult are the KERNEL variables (already-in-file
# convention: 𝒦(u,v) = mult(u,v)·K(u,v)). The WEIGHT side: Wx = product of ψ3-derivative factors.
# Final: (∂x∂y)^a C_r = 2π·Σ_terms ϑ^{(letters u)+(letters v)}... NO wait — every ∂x = u'·∂u with
# u' = ϑψ3'(x) contributes ONE ϑ; so total ϑ-power = a (for ∂x) + a (for ∂y) = 2a? For a=1: (48) has
# ϑ² ✓ (2πϑ²p3p p3p K). For a=2: two ϑ from x-derivatives + two from y = ϑ⁴ ✓ (with p3p weights).
# BUT the u-letters in kind correspond to ∂u's applied to K — each such ∂u arose from ONE ∂x = u'∂u:
# so ϑ-power = (total # of ∂x) + (# ∂y) = a + a = 2a — CONSISTENT. However careful: ∂x of a WEIGHT
# ψ3'(x) gives ψ3''(x) with NO ϑ — right.
# ---------------------------------------------------------------------------
def dxdy_terms(a: int):
    """Symbolic list for (∂x∂y)^a C_r evaluated at (x,y): entries dict kind → list of (Wx_fn, Wy_fn).
    Representative: Σ_kinds Σ_pairs 2π ϑ²ᵃ Wx(x)·Wy(y)·mult(kind)·K_r(u,v)... EXCEPT the leading 2π
    factor: (∂x∂y)^a C_r = (π/2)·(ϑ²)ᵃ·Σ ... — paper a=1: (48) = 2πϑ²p3p p3p K ⇒ C-side collections:
    C = (π/2)Σq q; ∂x∂yC = (π/2)·(ϑ²p3p²)·Σ(q'q') = 2πϑ²p3p²K R≡ (48). The (π/2) ↔ 2π conversion happens
    ONCE via the K-bridge: for a≥1 terms KEEP the 2πϑ^{2a}·(weights)·(mult·K) form. Verify in code."""
    # build iteratively: start from [(oneW, oneW, "K")] representing C = (π/2)[1·1·K₀]—NO: C itself is
    # a series, not K — but we only ever need (∂x∂y)^a C_r for a ≥ 1 (ℐ excludes (0,0)) — and for a ≥ 1
    # Lemma C.1 gives the CLOSED form starting from 2πϑ²ψ3'ψ3'K. So for a ≥ 1, build from the a=1 base:
    #   base(a=1) = [(psi3p, psi3p, "K")] with prefactor 2πϑ².
    # Each ∂x∂y application maps entry (Wx, Wy, kind):
    #   ∂x: [ (Wx', Wy, kind), (Wx·ϑp3-chain?, ...) ] — concretely:
    #     ∂x[Wx(x)·Wy(y)·𝒦kind(u,v)] = Wx'·Wy·𝒦kind + Wx·ϑψ3'(x)·Wy·𝒦kind+u
    #   likewise ∂y. So each entry expands to 4: (Wx', Wy, k), (Wx, Wy', k), (ϑψ3'-Wx, Wy, k+u),
    #     (Wx, ϑψ3'-Wy, k+v) — with ϑ powers tracked in a per-entry scalar.
    def bump(kind, ch):
        # kind "K" + suffix letters; adding one 'u' or 'v' — canonical order: r's, u's, v's
        p = kind.count('r'); i = kind.count('u'); j = kind.count('v')
        assert p + i + j == len(kind) - 1 and kind[0] == 'K'
        if ch == 'u':
            i += 1
        else:
            j += 1
        return "K" + "r" * p + "u" * i + "v" * j

    entries = [(psi3pf, psi3pf, "K", 1)]  # (Wx, Wy, kind, scalar-coeff without 2πϑ²)
    for _ in range(a - 1):
        new = []
        for Wx, Wy, kind, sc in entries:
            new.append((psi3ppf, Wy, kind, sc))          # ∂x on Wx
            new.append((lambda x: psi3pf(x) * psi3pf(x) if False else None, Wy, bump(kind, 'u'), sc))  # placeholder fixed below
        entries = new
    raise NotImplementedError
