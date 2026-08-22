"""Robustness of c*: does a richer adversary beat the claimed optimum?

The published record c* = 0.3823455333667027 is the value of

    max_{alpha in [0,1]}  inf_{P_S, P_SR}  [ (1-alpha) E h(Sbar Tbar)
                                             + alpha E h(s*(S,R)) - E h(S) ]

subject to E[S] <= t, where T is an independent copy of S and P_SR is any
symmetric coupling of P_S with itself.  Yu reduces the inner infimum to a
5-parameter family via Krein-Milman; Yu and Cambie then locate the optimum
NUMERICALLY.

If that reduction were merely an ansatz, a law outside the family could drive
the slack negative at t = c*, and c* would not be a valid lower bound at all.
This file searches k-atom laws for k = 2..7 with NO structural restriction.

The inner minimisation over couplings is EXACT: for fixed atoms and weights,
E h(s*) is linear in the joint matrix M, and {M >= 0, M1 = w, M^T 1 = w} is the
transportation polytope, so the worst coupling is an LP.  (Restricting to
SYMMETRIC couplings loses nothing: h(s*) is symmetric, so if M is optimal so is
(M + M^T)/2, with the same cost and the same marginals.)
"""

import os
import sys

import numpy as np
from scipy.optimize import linprog, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

C_STAR = 0.382345533366702721
ALPHA_STAR = 0.0356069
PSI = 0.3819660112501051517954132
B_STAR = 0.329454738503036972
A_STAR = 0.0788772927059231734

LOG2 = np.log(2.0)


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def maxent_or(pv, rv):
    """s*(p,r) = median{1/2, max(p,r), min(p+r,1)}."""
    return np.clip(0.5, np.maximum(pv, rv), np.minimum(pv + rv, 1.0))


def _ot_min(costs, w):
    """min_M <costs, M> over couplings of w with itself.  Exact LP."""
    k = len(w)
    Aeq = np.zeros((2 * k, k * k))
    for i in range(k):
        Aeq[i, i * k:(i + 1) * k] = 1.0
        Aeq[k + i, i::k] = 1.0
    beq = np.concatenate([w, w])
    res = linprog(costs.ravel(), A_eq=Aeq, b_eq=beq,
                  bounds=(0, None), method="highs")
    return float(res.fun) if res.success else np.inf


def slack(atoms, w, alpha):
    atoms = np.asarray(atoms, float)
    w = np.asarray(w, float)
    orm = atoms[:, None] + atoms[None, :] - atoms[:, None] * atoms[None, :]
    iid = float(w @ hf(orm.ravel()).reshape(orm.shape) @ w)
    cpl = _ot_min(hf(maxent_or(atoms[:, None], atoms[None, :]).ravel()
                     ).reshape(len(w), len(w)), w)
    return (1 - alpha) * iid + alpha * cpl - float(w @ hf(atoms))


def min_slack(t, alpha, k, restarts, seed):
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)

    def unpack(z):
        a = 1.0 / (1.0 + np.exp(-np.clip(z[:k], -30, 30)))
        e = np.exp(z[k:] - z[k:].max())
        return a, e / e.sum()

    def obj(z):
        return slack(*unpack(z), alpha)

    def con(z):
        a, w = unpack(z)
        return t - float(w @ a)

    for _ in range(restarts):
        z0 = rng.normal(0, 2.5, 2 * k)
        try:
            r = minimize(obj, z0, method="SLSQP",
                         constraints=[{"type": "ineq", "fun": con}],
                         options={"maxiter": 200, "ftol": 1e-13})
        except Exception:
            continue
        if r.success and r.fun < best[0]:
            a, w = unpack(r.x)
            best = (float(r.fun), a, w)
    return best


# ---------------------------------------------------------------------------
print("1. Sanity: the published extremal law reproduces slack 0 at t = c*")
atoms0 = np.array([1.0, B_STAR])
w0 = np.array([A_STAR, 1 - A_STAR])
print(f"   mean = {float(w0 @ atoms0):.16f}   (c* = {C_STAR:.16f})")
print(f"   slack at alpha* = {slack(atoms0, w0, ALPHA_STAR):+.3e}")

# ---------------------------------------------------------------------------
print("\n2. Unrestricted k-atom search at t = c* (alpha = alpha*)")
print(f"   {'k':>3} {'min slack':>14}   extremal law")
worst = np.inf
for k in range(2, 8):
    val, a, w = min_slack(C_STAR, ALPHA_STAR, k, restarts=60, seed=100 + k)
    worst = min(worst, val)
    keep = w > 1e-3
    desc = "  ".join(f"{wi:.4f}@{ai:.5f}" for ai, wi in
                     sorted(zip(a[keep], w[keep]), key=lambda s: -s[1]))
    print(f"   {k:>3} {val:>14.3e}   {desc}")
print(f"   worst over all k: {worst:+.3e}")

# ---------------------------------------------------------------------------
print("\n3. Is c* really the threshold?  Sweep t, maximising over alpha")
print(f"   {'t':>20} {'best over alpha':>16} {'argmax alpha':>13}")
alphas = np.linspace(0.0, 0.12, 13)
for dt, lbl in [(-4e-4, "c*-4e-4 (<psi)"), (-2e-4, "c*-2e-4"),
                (-2e-5, "c*-2e-5"), (0.0, "c*"), (+2e-5, "c*+2e-5"),
                (+2e-4, "c*+2e-4")]:
    t = C_STAR + dt
    best_a, best_v = None, -np.inf
    for al in alphas:
        v, _, _ = min_slack(t, al, 4, restarts=25, seed=7)
        if v > best_v:
            best_v, best_a = v, al
    flag = "" if best_v > -1e-9 else "   <-- BROKEN"
    print(f"   {t:.16f} {best_v:>16.3e} {best_a:>13.4f}{flag}   [{lbl}]")

# ---------------------------------------------------------------------------
print("\n4. Does the 5-parameter Yu family capture the true infimum?")
print("   Yu family: P_pq = (1-beta) Q_{a1,a2} + beta Q_{b1,b2},")
print("   Q_{x,y} = (delta_(x,y) + delta_(y,x))/2.")


def yu_slack(params, alpha):
    a1, a2, b1, b2, beta = params
    atoms = np.array([a1, a2, b1, b2])
    w = np.array([(1 - beta) / 2, (1 - beta) / 2, beta / 2, beta / 2])
    orm = atoms[:, None] + atoms[None, :] - atoms[:, None] * atoms[None, :]
    iid = float(w @ hf(orm.ravel()).reshape(4, 4) @ w)
    # The family's own coupling: pairs (a1,a2) with weight 1-beta, (b1,b2) with beta.
    cpl = (1 - beta) * float(hf(maxent_or(a1, a2))) + beta * float(hf(maxent_or(b1, b2)))
    return (1 - alpha) * iid + alpha * cpl - float(w @ hf(atoms))


def yu_mean(params):
    a1, a2, b1, b2, beta = params
    return (1 - beta) * (a1 + a2) / 2 + beta * (b1 + b2) / 2


rng = np.random.default_rng(0)
best = (np.inf, None)
for _ in range(4000):
    z = rng.uniform(0, 1, 5)
    r = minimize(lambda p: yu_slack(p, ALPHA_STAR), z, method="SLSQP",
                 bounds=[(0, 1)] * 5,
                 constraints=[{"type": "ineq", "fun": lambda p: C_STAR - yu_mean(p)}],
                 options={"maxiter": 200, "ftol": 1e-14})
    if r.success and r.fun < best[0]:
        best = (float(r.fun), r.x)
print(f"   Yu-family min slack at t = c*: {best[0]:+.3e}")
print(f"   at (a1,a2,b1,b2,beta) = {np.array2string(best[1], precision=6)}")
print(f"   mean = {yu_mean(best[1]):.16f}")
print(f"   general k-atom search found  : {worst:+.3e}")
gap = worst - best[0]
print(f"   difference (general - family): {gap:+.3e}")
if gap < -1e-7:
    print("   => richer laws BEAT the family; the reduction would be unsound.")
else:
    print("   => the family attains the infimum, to numerical accuracy.")

print("\nDONE")
