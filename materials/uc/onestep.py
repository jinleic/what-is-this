"""The independent ("one-step") union-closed functional, verified from scratch.

Everything here is re-derived rather than transcribed.  Three claims:

  (1) EXACT.  psi = (3-sqrt5)/2 is the exact threshold of h(2p-p^2) >= h(p),
      because 1-(2p-p^2) = (1-p)^2 and h(x)=h(1-x), so the condition collapses
      to (1-p)^2 >= p on [0,1/2].  Verified symbolically in sympy.

  (2) EXACT.  Sawin's two branches both equal 1 at psi:
        h(2psi-psi^2)/h(psi) = 1   since 2psi-psi^2 = 1-(1-psi)^2 = 1-psi,
        PHI*(1-psi)         = 1   since 1-psi = 1/PHI.

  (3) NUMERICAL.  The variational problem
        G(mu) = E_{p,q~mu(x)mu} h(p+q-pq) - E_{p~mu} h(p),   E_mu[p] <= t
      has inf_mu G(mu) >= 0 exactly for t <= psi.  We locate the minimising
      mu by multistart optimisation over k-atom laws and report which law wins,
      rather than assuming the extremiser's shape.

Run directly; every assertion is checked.
"""

import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PHI, PSI, h, mp, mpf, onestep_gain, onestep_mean, sawin_lambda

# ---------------------------------------------------------------------------
# (1) psi is EXACTLY the threshold, symbolically
# ---------------------------------------------------------------------------
print("1. psi = (3-sqrt5)/2 is the exact one-step threshold")

p = sp.Symbol("p", positive=True)
psi_sym = (3 - sp.sqrt(5)) / 2

# The reflection identity that drives everything.
assert sp.simplify(1 - (2 * p - p**2) - (1 - p) ** 2) == 0
print("   1-(2p-p^2) = (1-p)^2                                     [OK]")

# h(2p-p^2) = h(p)  <=>  (1-p)^2 = p  (nontrivial branch), whose root is psi.
assert sp.simplify(psi_sym**2 - 3 * psi_sym + 1) == 0
assert sp.simplify((1 - psi_sym) ** 2 - psi_sym) == 0
print("   psi^2-3psi+1 = 0  and  (1-psi)^2 = psi                    [OK]")

roots = sp.solve(sp.Eq((1 - p) ** 2, p), p)
assert psi_sym in [sp.nsimplify(sp.radsimp(r)) for r in roots], roots
print(f"   (1-p)^2 = p has roots {roots}, containing psi            [OK]")

# The nontrivial reflection is the ONLY way h(a)=h(b) with a!=b on [0,1].
sym_p = sp.nsimplify(psi_sym)
print(f"   psi = {sp.N(sym_p, 25)}")

# ---------------------------------------------------------------------------
# (2) Sawin's branches meet at 1, exactly
# ---------------------------------------------------------------------------
print("\n2. Sawin's lambda(u): both branches equal 1 at u = psi")

# 2psi - psi^2 = 1 - psi, so h(2psi-psi^2) = h(1-psi) = h(psi).
assert sp.simplify((2 * psi_sym - psi_sym**2) - (1 - psi_sym)) == 0
print("   2psi-psi^2 = 1-psi, so h(2psi-psi^2) = h(psi)             [OK]")

phi_sym = (1 + sp.sqrt(5)) / 2
assert sp.simplify(phi_sym * (1 - psi_sym) - 1) == 0
print("   PHI*(1-psi) = 1                                           [OK]")

lo = sawin_lambda(PSI - mpf(10) ** -12)
hi = sawin_lambda(PSI + mpf(10) ** -12)
assert abs(lo - 1) < mpf(10) ** -10, lo
assert abs(hi - 1) < mpf(10) ** -10, hi
print(f"   numerically lambda(psi-) = {mp.nstr(lo, 12)}, "
      f"lambda(psi+) = {mp.nstr(hi, 12)}    [OK]")

# lambda < 1 strictly above psi: the method dies exactly there.
assert sawin_lambda(mpf("0.40")) < 1
assert sawin_lambda(mpf("0.35")) > 1
print("   lambda(0.35) > 1 > lambda(0.40)                           [OK]")

# ---------------------------------------------------------------------------
# (3) Which measure is extremal?  Search, do not assume.
# ---------------------------------------------------------------------------
print("\n3. Variational problem: inf over k-atom laws with mean <= t")

LOG2 = np.log(2.0)


def hf(x):
    """float64 binary entropy, safe at the endpoints."""
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def gain_f(atoms, weights):
    a = np.asarray(atoms, dtype=float)
    w = np.asarray(weights, dtype=float)
    orm = a[:, None] + a[None, :] - a[:, None] * a[None, :]
    pair = w @ hf(orm.ravel()).reshape(orm.shape) @ w
    return pair - w @ hf(a)


def min_gain(t, k, restarts=400, seed=0):
    """Minimise G over k-atom laws with mean <= t.  Returns (value, atoms, weights)."""
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)

    def unpack(z):
        atoms = 1.0 / (1.0 + np.exp(-z[:k]))
        e = np.exp(z[k:] - z[k:].max())
        return atoms, e / e.sum()

    def obj(z):
        return gain_f(*unpack(z))

    def mean_con(z):
        atoms, w = unpack(z)
        return t - float(w @ atoms)          # >= 0

    for _ in range(restarts):
        z0 = rng.normal(0, 2.0, 2 * k)
        try:
            r = minimize(obj, z0, method="SLSQP",
                         constraints=[{"type": "ineq", "fun": mean_con}],
                         options={"maxiter": 500, "ftol": 1e-14})
        except Exception:
            continue
        if r.success and r.fun < best[0]:
            atoms, w = unpack(r.x)
            best = (float(r.fun), atoms, w)
    return best


print(f"   {'t':>10} {'k':>3} {'min G':>14}   extremal law (weight@atom)")
results = {}
for t in ["0.30", "0.37", "0.381966011250105", "0.39", "0.45"]:
    tf = float(t)
    row = []
    for k in (1, 2, 3, 4):
        val, atoms, w = min_gain(tf, k, restarts=200, seed=k)
        row.append(val)
        if k == 4:
            keep = w > 1e-4
            desc = "  ".join(
                f"{wi:.4f}@{ai:.6f}" for ai, wi in
                sorted(zip(atoms[keep], w[keep]), key=lambda s: -s[1])
            )
            print(f"   {t:>10} {k:>3} {val:>14.3e}   {desc}")
    results[t] = row

# The threshold: G can be driven negative iff t > psi.
psi_f = float(PSI)
v_below, _, _ = min_gain(psi_f - 1e-4, 4, restarts=300, seed=11)
v_above, _, _ = min_gain(psi_f + 1e-4, 4, restarts=300, seed=12)
print(f"\n   inf G at t = psi - 1e-4 : {v_below:+.3e}")
print(f"   inf G at t = psi + 1e-4 : {v_above:+.3e}")
assert v_below > -1e-9, f"method should survive below psi, got {v_below}"
assert v_above < -1e-9, f"method should fail above psi, got {v_above}"
print("   => one-step reach is exactly psi                          [OK]")

# At t = psi the extremiser is the point mass delta_psi (mean constraint tight).
val, atoms, w = min_gain(psi_f, 1, restarts=50, seed=3)
assert abs(atoms[0] - psi_f) < 1e-6, atoms
assert abs(val) < 1e-9, val
print(f"   delta_psi gives G = {val:+.2e} with atom {atoms[0]:.9f}    [OK]")

# High-precision confirmation of the point-mass value.
g_exact = onestep_gain([PSI], [mpf(1)])
assert abs(g_exact) < mpf(10) ** -30, g_exact
print(f"   mpmath G(delta_psi) = {mp.nstr(g_exact, 6)} (dps={mp.dps})     [OK]")

print("\nALL PASS")
