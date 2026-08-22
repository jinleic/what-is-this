"""Is P_p -> E h(pbar qbar) concave ON THE FEASIBLE SET?

Yu (arXiv:2212.00658, Section 4) reduces the union-closed variational problem
to a 5-parameter family by Krein-Milman, resting on:

    "Note that Ppq -> g(Ppq, alpha) is concave, since by [4, Lemma 5],
     Pp -> E_{(p,q) ~ Pp^{tensor 2}} h(p+q-pq) is concave, and Ppq -> Pp is
     linear."

and Krein-Milman is applied to

    P_B := { symmetric Ppq concentrated on B^2 such that E[p] <= t }.

SCOPE MATTERS.  Concavity is only ever used ON P_B, i.e. for laws obeying the
mean constraint E[p] <= t.  Writing x = 1-p and F(mu) = int int h(xy) dmu dmu,

    F concave on a convex set S
      <=>  Q(mu_1 - mu_2) <= 0 for all mu_1, mu_2 in S,
      Q(nu) := int int h(xy) dnu(x) dnu(y),

since F(lam mu_1 + lam' mu_2) - lam F(mu_1) - lam' F(mu_2) = -lam lam' Q(nu).

WITHOUT the mean constraint F is NOT concave: mu_1 = delta_{0.3},
mu_2 = delta_1 (i.e. x = 0.7 and x = 0) give F = 0.9997, 0 and midpoint
F = 0.2499 < 0.4999.  But that pair has E[p] = 0.65, far outside E[p] <= t,
so it says nothing about Yu's argument.  The question that matters is whether
concavity holds on the FEASIBLE set, and that is what this file tests.
"""

import numpy as np
from scipy.optimize import minimize

LOG2 = np.log(2.0)
T_STAR = 0.382345533366702721
PSI = 0.3819660112501051517954132


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def F_p(atoms, w):
    """E h(p+q-pq) for the law sum w_i delta_{atoms_i}, in p-coordinates."""
    a = np.asarray(atoms, float)
    w = np.asarray(w, float)
    orm = a[:, None] + a[None, :] - a[:, None] * a[None, :]
    return float(w @ hf(orm.ravel()).reshape(len(w), len(w)) @ w)


# ---------------------------------------------------------------------------
print("0. The unrestricted kernel is not conditionally negative semidefinite")
print("   witness (OUT OF SCOPE, mean 0.65 > t): mu1 = d_0.3, mu2 = d_1")
f1, f2 = F_p([0.3], [1.0]), F_p([1.0], [1.0])
fm = F_p([0.3, 1.0], [0.5, 0.5])
print(f"   F(mu1) = {f1:.6f}  F(mu2) = {f2:.6f}  F(mid) = {fm:.6f}"
      f"  defect = {(f1+f2)/2 - fm:+.6f}")
print("   => concavity genuinely FAILS globally; scope is what saves it.")

# ---------------------------------------------------------------------------
print("\n1. Two-point sub-question, exactly")
print("   mu_i = delta_{p_i} with p_i <= t gives, in x = 1-p >= 1-t,")
print("      Q <= 0  <=>  h(x^2) + h(y^2) <= 2 h(xy)   for x,y in [1-t, 1].")
lo = 1 - T_STAR
worst, arg = -np.inf, None
g = np.linspace(lo, 1.0, 1400)
for x in g:
    d = hf(x * x) + hf(g * g) - 2 * hf(x * g)
    j = int(np.argmax(d))
    if d[j] > worst:
        worst, arg = float(d[j]), (x, float(g[j]))
print(f"   max over [1-t,1]^2 of h(x^2)+h(y^2)-2h(xy) = {worst:+.3e}")
print(f"   attained at (x,y) = ({arg[0]:.6f}, {arg[1]:.6f})")
print(f"   (equality along x = y, so 0 is the expected maximum)")
assert worst < 1e-9, "two-point concavity fails on the feasible range"
print("   => two-point case: concavity HOLDS on the feasible range   [OK]")

print("\n   asymptotics as x,y -> 1 (the tight corner):")
for eps in (1e-2, 1e-3, 1e-4, 1e-6):
    x = 1 - eps
    d = float(hf(x * x) + hf(1.0) - 2 * hf(x))
    print(f"     x = 1-{eps:<7g}, y = 1: h(x^2)-2h(x) = {d:+.3e}")

# ---------------------------------------------------------------------------
print("\n2. Full search: maximise Q(mu1 - mu2) over FEASIBLE pairs")
print("   both mu_i are k-atom laws in p-coordinates with mean <= t")


def defect(z, k, t):
    """(F(mu1)+F(mu2))/2 - F(midpoint);  > 0 means concavity fails."""
    a1 = 1 / (1 + np.exp(-np.clip(z[0:k], -30, 30)))
    w1 = np.exp(z[k:2 * k] - z[k:2 * k].max()); w1 /= w1.sum()
    a2 = 1 / (1 + np.exp(-np.clip(z[2 * k:3 * k], -30, 30)))
    w2 = np.exp(z[3 * k:4 * k] - z[3 * k:4 * k].max()); w2 /= w2.sum()
    mid_a = np.concatenate([a1, a2])
    mid_w = np.concatenate([w1, w2]) / 2
    return (F_p(a1, w1) + F_p(a2, w2)) / 2 - F_p(mid_a, mid_w), (a1, w1, a2, w2)


for t in (PSI, T_STAR, 0.45, 0.5):
    rng = np.random.default_rng(4)
    best, barg = -np.inf, None
    for k in (1, 2, 3):
        for _ in range(300):
            z0 = rng.normal(0, 2.5, 4 * k)

            def neg(z):
                return -defect(z, k, t)[0]

            def c1(z):
                a1 = 1 / (1 + np.exp(-np.clip(z[0:k], -30, 30)))
                w1 = np.exp(z[k:2 * k] - z[k:2 * k].max()); w1 /= w1.sum()
                return t - float(w1 @ a1)

            def c2(z):
                a2 = 1 / (1 + np.exp(-np.clip(z[2 * k:3 * k], -30, 30)))
                w2 = np.exp(z[3 * k:4 * k] - z[3 * k:4 * k].max()); w2 /= w2.sum()
                return t - float(w2 @ a2)

            try:
                r = minimize(neg, z0, method="SLSQP",
                             constraints=[{"type": "ineq", "fun": c1},
                                          {"type": "ineq", "fun": c2}],
                             options={"maxiter": 200, "ftol": 1e-14})
            except Exception:
                continue
            if r.success and -r.fun > best:
                best, barg = -float(r.fun), defect(r.x, k, t)[1]
    tag = "CONCAVE" if best < 1e-9 else "*** FAILS ***"
    print(f"   t = {t:.10f}: max defect = {best:+.3e}   {tag}")
    if best > 1e-9:
        a1, w1, a2, w2 = barg
        print(f"      mu1 atoms {np.array2string(a1, precision=6)} "
              f"w {np.array2string(w1, precision=4)} mean {w1@a1:.6f}")
        print(f"      mu2 atoms {np.array2string(a2, precision=6)} "
              f"w {np.array2string(w2, precision=4)} mean {w2@a2:.6f}")

# ---------------------------------------------------------------------------
print("\n3. Where does concavity start to fail?  Sweep the mean bound.")
print(f"   {'t':>8} {'max defect':>13}")
rng = np.random.default_rng(11)
for t in (0.38, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65):
    best = -np.inf
    for k in (1, 2):
        for _ in range(200):
            z0 = rng.normal(0, 2.5, 4 * k)

            def neg(z):
                return -defect(z, k, t)[0]

            def c1(z):
                a1 = 1 / (1 + np.exp(-np.clip(z[0:k], -30, 30)))
                w1 = np.exp(z[k:2 * k] - z[k:2 * k].max()); w1 /= w1.sum()
                return t - float(w1 @ a1)

            def c2(z):
                a2 = 1 / (1 + np.exp(-np.clip(z[2 * k:3 * k], -30, 30)))
                w2 = np.exp(z[3 * k:4 * k] - z[3 * k:4 * k].max()); w2 /= w2.sum()
                return t - float(w2 @ a2)

            try:
                r = minimize(neg, z0, method="SLSQP",
                             constraints=[{"type": "ineq", "fun": c1},
                                          {"type": "ineq", "fun": c2}],
                             options={"maxiter": 150, "ftol": 1e-13})
            except Exception:
                continue
            if r.success and -r.fun > best:
                best = -float(r.fun)
    print(f"   {t:>8.3f} {best:>13.3e}")

print("\nDONE")
