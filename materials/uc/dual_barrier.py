"""Dual certificates for the coupled term: all of them cap at psi.  (Corrected.)

RETRACTION.  An earlier version of this file claimed "linear dual certificates
cannot exceed psi" by requiring copositivity of a kernel Q and evaluating it at
delta_0, delta_1 and delta_s.  That CONFLATED copositivity over all probability
measures with copositivity over { E_mu[p] <= t }, and was wrong both ways:

  * Global reading (what a multiplier actually forces): delta_s is admissible for
    EVERY s.  Since gamma = 0 is forced and phi <= 0, this needs
    h(2s-s^2) >= h(s) for every s, FALSE at s = 0.45.  The class is EMPTY, so
    nothing is "capped"; and the table of "largest t with Q(t,t) >= 0" was just
    re-deriving lambda(t) >= 1 at the single point s = t, which is the DEFINITION
    of psi and has no dual content.  That was the source of the spurious
    "alpha = 0 recovers psi exactly" punchline.
  * Constrained reading: delta_1 has mean 1 > t, so Q(1,1) >= 0 is not required,
    phi(1) = 0 does not follow, and the old proof dies at step 2.

The correct results are below.  THEOREM C is new and settles the constrained
case, which the retraction had left open: the route is closed, by a different
and much simpler mechanism -- the SINK PAIR {0,1}.

-------------------------------------------------------------------------------
SETUP.  Q(mu) = int int h(p (+) q) dmu dmu,  L(mu) = int h dmu, and

    C(mu) = min over symmetric couplings M of mu with itself of int h(s*) dM,
    s*(p,r) = median{1/2, max(p,r), min(p+r,1)}.

Monge-Kantorovich duality on the transportation polytope gives C(mu) >= 2 int phi
dmu for every phi obeying

    (A)   phi(p) + phi(r) <= h(s*(p,r))     for all p, r in [0,1].

A "linear dual certificate at t" is a pair (alpha, phi) with alpha in (0,1] and
phi obeying (A) such that

    Fhat(mu) := (1-alpha) Q(mu) + 2 alpha int phi dmu - L(mu)  >=  0
                                        for every mu with E_mu[p] <= t.       (D)

-------------------------------------------------------------------------------
THEOREM C (sink-pair barrier).  If (alpha, phi) is a linear dual certificate at
t, then t <= psi = (3-sqrt5)/2.  This holds even when phi is allowed to depend on
E_mu[p], because the whole argument lives inside the single slice of mean t.

Proof.  Consider the SINK MIXTURE

    nu_t = (1-t) delta_0 + t delta_1,        E[p] = t,  so nu_t is feasible.

Every entropy term vanishes on it: h(0) = h(1) = 0, and 0 (+) 0 = 0 while
0 (+) 1 = 1 (+) 1 = 1, so Q(nu_t) = L(nu_t) = 0.  Hence (D) reads

    Fhat(nu_t) = 2 alpha [ (1-t) phi(0) + t phi(1) ] >= 0.

But s*(0,0) = 0 and s*(1,1) = 1, both with h = 0, so (A) forces phi(0) <= 0 AND
phi(1) <= 0.  A nonnegative combination of two nonpositive numbers is >= 0 only
if both vanish:

    phi(0) = phi(1) = 0.

Since s*(p,1) = 1 for EVERY p, (A) at (p,1) now gives

    phi(p) <= -phi(1) = 0     for all p in [0,1].

Finally delta_t is feasible, and (D) at delta_t gives

    (1-alpha) h(2t - t^2) + 2 alpha phi(t) - h(t) >= 0,

so with phi(t) <= 0 we need (1-alpha) h(2t-t^2) >= h(t), i.e.

    lambda(t) := h(2t-t^2)/h(t) >= 1/(1-alpha) >= 1.

lambda is strictly decreasing with lambda(psi) = 1, so t <= psi.            []

INTERPRETATION.  The sink pair {0,1} is a two-point set on which the ENTIRE
entropy functional vanishes identically, and the mixtures nu_u, u in [0,t], are
feasible for every mean u <= t.  So a linear surrogate for C must be exactly
right on a whole curve of tight measures, which pins it to 0 at both endpoints;
(A) then drags it below 0 everywhere.  Point masses alone then cap the iid term
at psi.  This is why the literature cannot avoid a support reduction: no linear
lower bound on the coupled term -- localised by the mean or not -- can help.

-------------------------------------------------------------------------------
PROPOSITION D (decoupled certificates).  The same conclusion holds, by a
different route, if the iid term is first replaced by Sawin's sharp bound
Q >= lambda(t) L and the coupled term by 2 int phi dmu, then the mean constraint
is discharged with a multiplier: the multiplier then acts on a LINEAR functional,
where Lagrangian duality is tight, and the requirement G(s) + gamma (s-t) >= 0
with G = [(1-alpha)lambda(t) - 1] h + 2 alpha phi forces gamma = 0, phi(1) = 0,
phi <= 0, and finally lambda(t) >= 1/(1-alpha).
"""

import os
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI, h, mp, mpf, sawin_lambda


def maxent_or(p, r):
    p, r = mpf(p), mpf(r)
    return min(max(mpf("0.5"), max(p, r)), min(p + r, mpf(1)))


def union(p, r):
    p, r = mpf(p), mpf(r)
    return p + r - p * r


# ---------------------------------------------------------------------------
print("PART 1 -- why the GLOBAL-multiplier class is EMPTY (the retraction)")
print("  A multiplier absorbing E[p] <= t demands the inequality for ALL")
print("  probability measures, so delta_s is admissible for every s.  Then")
print("  gamma = 0 is forced (delta_0 with phi(0) <= 0) and phi <= 0, so")
print("  Q(s,s) <= (1-a) h(2s-s^2) - h(s) must be >= 0 for EVERY s.")
print(f"\n  {'s':>7} {'h(2s-s^2)':>16} {'h(s)':>16} {'difference':>14}")
witness = None
for s in ("0.30", "0.3819660112501051", "0.40", "0.45", "0.49"):
    sv = mpf(s)
    lo, hi = h(union(sv, sv)), h(sv)
    d = lo - hi
    print(f"  {s:>7} {mp.nstr(lo, 14):>16} {mp.nstr(hi, 14):>16} "
          f"{mp.nstr(d, 10):>14}")
    if d < 0 and witness is None:
        witness = (sv, d)
sv, d = witness
print(f"\n  witness s = {mp.nstr(sv, 6)}: difference {mp.nstr(d, 8)} < 0, and")
print("  (1-a) only shrinks the positive term, so this holds for EVERY a >= 0:")
for alpha in ("0.0", "0.0356069", "0.5", "0.9"):
    al = mpf(alpha)
    q = (1 - al) * h(union(sv, sv)) - h(sv)
    assert q < 0, (alpha, q)
    print(f"    alpha = {alpha:>10}: Q(s,s) <= {mp.nstr(q, 8)} < 0")
print("  => NO certificate of that form exists at ANY t > 0; the class is")
print("     empty and the retracted table measured lambda(t) >= 1, not")
print("     copositivity.                                                [OK]")

# ---------------------------------------------------------------------------
print("\nPART 2 -- THEOREM C: the sink pair closes the CONSTRAINED case too")

print("\n  (i) the sink pair kills every entropy term")
for p, r in ((0, 0), (0, 1), (1, 1)):
    print(f"      h({p}) = {h(p)}   {p} (+) {r} = {union(p, r)}   "
          f"h({p} (+) {r}) = {h(union(p, r))}   s*({p},{r}) = {maxent_or(p, r)}")
assert h(0) == 0 and h(1) == 0
assert union(0, 0) == 0 and union(0, 1) == 1 and union(1, 1) == 1
assert h(union(0, 0)) == 0 and h(union(0, 1)) == 0 and h(union(1, 1)) == 0
print("      so Q(nu_t) = L(nu_t) = 0 for nu_t = (1-t) d_0 + t d_1     [OK]")

print("\n  (ii) nu_t is feasible for every t, with mean exactly t")
for t in ("0.1", "0.3", "0.3819660112501051", "0.45"):
    tv = mpf(t)
    mean = (1 - tv) * 0 + tv * 1
    assert mean == tv
    print(f"      t = {t:>20}: mean(nu_t) = {mp.nstr(mean, 18)}  = t")
print("      so the argument stays inside the slice of mean t, hence it")
print("      applies even to a phi that depends on the mean            [OK]")

print("\n  (iii) (A) forces phi(0) <= 0 and phi(1) <= 0")
s00, s11 = maxent_or(0, 0), maxent_or(1, 1)
print(f"      s*(0,0) = {s00}, h = {h(s00)}  => 2 phi(0) <= 0")
print(f"      s*(1,1) = {s11}, h = {h(s11)}  => 2 phi(1) <= 0")
assert s00 == 0 and s11 == 1 and h(s00) == 0 and h(s11) == 0
print("      a nonnegative combination of two nonpositive numbers is >= 0")
print("      only if both vanish  =>  phi(0) = phi(1) = 0              [OK]")

print("\n  (iv) s*(p,1) = 1 for every p, so (A) at (p,1) gives phi <= 0")
for p in ("0", "0.2", "0.5", "0.8", "1"):
    v = maxent_or(p, 1)
    assert v == 1 and h(v) == 0, (p, v)
    print(f"      s*({p:>4},1) = {v}, h = {h(v)}")
print("      phi(p) <= -phi(1) = 0 everywhere                          [OK]")

print("\n  (v) delta_t then forces lambda(t) >= 1/(1-alpha) >= 1")
print(f"      {'t':>22} {'lambda(t)':>18} {'>= 1?':>7}")
for t in ("0.30", "0.35", "0.3819660112501051", "0.3823455333667027", "0.40"):
    lam = sawin_lambda(t)
    print(f"      {t:>22} {mp.nstr(lam, 14):>18} {str(lam >= 1):>7}")

xs = [mpf(i) / 2000 for i in range(1, 2000)]
lams = [sawin_lambda(x) for x in xs]
assert all(lams[i] > lams[i + 1] for i in range(len(lams) - 1))
print("      lambda strictly decreasing on (0,1) (2000-point check)    [OK]")

p_sym = sp.Symbol("p", positive=True)
psi_sym = (3 - sp.sqrt(5)) / 2
assert sp.simplify((2 * psi_sym - psi_sym**2) - (1 - psi_sym)) == 0
print("      lambda(psi) = 1 exactly, since 2psi-psi^2 = 1-psi (sympy)  [OK]")
assert sawin_lambda(PSI - mpf(10) ** -8) > 1
assert sawin_lambda(PSI + mpf(10) ** -8) < 1
print("      lambda(t) >= 1  <=>  t <= psi                              [OK]")
assert sawin_lambda(mpf("0.3823455333667027")) < 1
print("      lambda(c*) < 1, so no dual certificate reaches c*          [OK]")

# ---------------------------------------------------------------------------
print("\nPART 3 -- the sink family is tight for EVERY mean, not just t")
print("  nu_u = (1-u) d_0 + u d_1 has Q = C = L = 0 for every u, so the")
print("  tight set is a whole CURVE.  A linear surrogate must be exact on")
print("  all of it; that is what pins phi at both endpoints.")
print(f"\n  {'u':>8} {'Q':>10} {'L':>10} {'mean':>10} {'G(nu_u)':>12}")
for u in ("0", "0.1", "0.25", "0.3819660112501051"):
    uv = mpf(u)
    Q = ((1 - uv) ** 2 * h(union(0, 0)) + 2 * uv * (1 - uv) * h(union(0, 1))
         + uv**2 * h(union(1, 1)))
    L = (1 - uv) * h(0) + uv * h(1)
    assert Q == 0 and L == 0
    print(f"  {u:>8} {mp.nstr(Q, 3):>10} {mp.nstr(L, 3):>10} "
          f"{mp.nstr(uv, 8):>10} {mp.nstr(1 - uv, 8):>12}")
print("\n  G(nu_u) = 1-u sweeps [1-t, 1], so no single tangent to the convex")
print("  direction G^2 can be exact on the tight set either -- which is why")
print("  the support reduction must FREEZE G rather than linearise it.  [OK]")

print("\nALL PASS")
print("\nSUMMARY")
print("  * global-multiplier dual certificates: class is EMPTY (vacuous)")
print("  * decoupled certificates (Sawin lambda + linear coupled bound): <= psi")
print("  * THEOREM C: mean-constrained linear dual certificates: <= psi,")
print("    even with mean-dependent phi.  Mechanism: the sink pair {0,1}.")
print("  => every linear treatment of the coupled term is capped at psi;")
print("     a support reduction (Theorem B) is unavoidable.")
