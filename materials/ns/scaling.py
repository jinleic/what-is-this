"""Symbolic checks for Lemma 0 and Lemma 1 in README.md.

Lemma 1 assumes the first-kind condition alpha + beta = 1 and the standard pressure
scaling gamma = 2 alpha.  This script does NOT derive alpha + beta = 1 -- powers of
(T-t) are linearly independent, so the PDE only forces the profile coefficients sharing
each DISTINCT exponent to cancel, and terms may cancel in groups.  What is checked here:

  0. scaling invariance under u -> lambda u(lambda^2 t, lambda x) forces
     alpha = beta = 1/2 outright (Lemma 0);
  1. under the assumed alpha + beta = 1 and gamma = 2 alpha, the three inviscid
     exponents coincide, so dividing by (T-t)^{-(alpha+1)} leaves them t-independent;
  2. the viscous term then carries exactly (T-t)^(1-2beta);
  3. hence beta != 1/2 leaves an explicitly time-dependent equation -- the hypothesis
     Lemma 1 converts into Delta U = 0 by subtracting at two times.

The Liouville step is in README.md; this is a check on the algebra only.
"""

import sympy as sp

t, T, nu = sp.symbols("t T nu", positive=True)
al, be = sp.symbols("alpha beta", real=True)
tau = T - t  # time to blowup

# Formal scaling weights. A term written as (T-t)^e * (profile) is tracked by e alone;
# every profile factor is order 1 by construction.
w_u = -al                  # u ~ tau^-alpha U
w_grad = -be               # d/dx = tau^-beta d/dy
w_p = -2 * al              # p ~ tau^-2alpha P  (forced by the nonlinearity)

e_dt = w_u - 1             # d/dt u : the chain rule adds one power of tau^-1
e_conv = 2 * w_u + w_grad  # (u.grad)u
e_press = w_p + w_grad     # grad p
e_visc = w_u + 2 * w_grad  # nu Laplacian u

print("Exponents of (T-t) in each term of the momentum equation")
for name, e in [("d_t u", e_dt), ("(u.grad)u", e_conv), ("grad p", e_press),
                ("nu*Lap u", e_visc)]:
    print(f"   {name:12s} {sp.simplify(e)}")

# 0. Lemma 0: invariance under u -> lam*u(lam^2 t, lam x) forces alpha = beta = 1/2.
lam = sp.Symbol("lam", positive=True)
# lam * u(lam^2 t, lam x) = lam^(1-2alpha) tau^-alpha U( lam^(1-2beta) x / tau^beta )
amp_exp = sp.simplify(1 - 2 * al)      # power of lam in front
arg_exp = sp.simplify(1 - 2 * be)      # power of lam inside U's argument
assert sp.solve(sp.Eq(amp_exp, 0), al) == [sp.Rational(1, 2)]
assert sp.solve(sp.Eq(arg_exp, 0), be) == [sp.Rational(1, 2)]
print("\n0. scaling invariance: lam^(1-2alpha) and lam^(1-2beta) must both be 1")
print("   => alpha = beta = 1/2 (Lemma 0; Leray's ansatz, nothing else)   [OK]")

# 1. First-kind condition. NOTE: this is a HYPOTHESIS of Lemma 1, not a consequence.
# Powers of (T-t) are linearly independent, so the PDE only forces coefficients sharing
# a distinct exponent to cancel; terms may cancel in groups. Here we merely record that
# alpha + beta = 1 is exactly the condition under which d_t u and (u.grad)u share an
# exponent, and check it is consistent.
sol = sp.solve(sp.Eq(e_dt, e_conv), be)
assert sol == [1 - al], sol
print(f"\n1. d_t u and (u.grad)u share an exponent iff beta = {sol[0]},")
print("   i.e. alpha + beta = 1  -- ASSUMED by Lemma 1, not derived        [OK]")

sub = {be: 1 - al}

# 2. Pressure rides along automatically at the same weight.
assert sp.simplify(e_press.subs(sub) - e_dt.subs(sub)) == 0
print("2. grad p carries the same weight as d_t u             [OK]")

# 3. Residual weight of the viscous term after dividing by (T-t)^(-(alpha+1)).
resid = sp.simplify(e_visc.subs(sub) - e_dt.subs(sub))
print(f"3. viscous term relative weight = {resid}")
assert sp.simplify(resid - (1 - 2 * (1 - al))) == 0
resid_beta = sp.simplify(resid.subs(al, 1 - be))
assert sp.simplify(resid_beta - (1 - 2 * be)) == 0
print(f"   in terms of beta:              {resid_beta}            [OK]")
# 4. The trichotomy. Reparametrize in tau = T - t and take tau -> 0+ directly: sympy
# does not know t < T, so limit(..., t, T, "-") on (T-t)^(-1/2) picks the wrong branch
# and returns -oo. Working in tau removes the ambiguity.
s = sp.symbols("tau_s", positive=True)
print("\n4. ratio nu*Lap/(u.grad)u = nu*tau^(1-2beta) as tau = T-t -> 0+")
lims = {}
for bval, label in [(sp.Rational(1, 4), "beta < 1/2"),
                    (sp.Rational(1, 2), "beta = 1/2"),
                    (sp.Rational(3, 4), "beta > 1/2")]:
    lims[bval] = sp.limit(nu * s ** (1 - 2 * bval), s, 0, "+")
    print(f"   {label}  (beta={bval}):  exponent {1-2*bval},  limit = {lims[bval]}")

assert lims[sp.Rational(1, 4)] == 0
assert lims[sp.Rational(1, 2)] == nu
assert lims[sp.Rational(3, 4)] == sp.oo
print("\n   beta<1/2 -> 0 (viscosity negligible, Re_l -> inf)      [OK]")
print("   beta=1/2 -> nu (critical; Leray ansatz)                [OK]")
print("   beta>1/2 -> oo (viscosity dominates, Re_l -> 0)        [OK]")

# Local Reynolds number at length l = tau^beta, velocity tau^-alpha, alpha = 1-beta.
Re = sp.simplify(s ** (-(1 - be)) * s**be / nu)
assert sp.simplify(Re - s ** (2 * be - 1) / nu) == 0
print(f"\n5. Re_l = U*l/nu = {Re}  = tau^(2beta-1)/nu             [OK]")
for bval in (sp.Rational(1, 4), sp.Rational(3, 4)):
    print(f"   beta={bval}:  Re_l -> {sp.limit(s ** (2 * bval - 1) / nu, s, 0, '+')}")

print("\nLemma 1 algebra verified. beta != 1/2 leaves an explicitly time-dependent")
print("equation; subtracting it at two times gives Delta U = 0, hence U == 0 by")
print("Liouville for functions vanishing at infinity. Only beta = 1/2 survives,")
print("which is exactly Leray's ansatz -- Lemma 2 territory.")
