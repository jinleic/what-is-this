"""The Sawin/Yu/Cambie coupling bound, derived and reproduced from scratch.

STRUCTURE OF THE METHOD.  Let p be the law of the conditional inclusion
probability.  Two strategies are combined:

  * IID:  B is an independent copy of A.  Coordinate OR probability p+q-pq.
  * MAX-ENTROPY:  C is coupled to A to maximise the entropy of the OR.
    Given marginals p,r the OR probability s = p+r-Pr[both] ranges over
    [max(p,r), min(p+r,1)] because Pr[both] in [max(0,p+r-1), min(p,r)].
    Entropy is maximised at s closest to 1/2, hence

        s*(p,r) = median{ 1/2, max(p,r), min(p+r,1) }.

    (Equivalently max(p, r, min(p+r, 1/2)); we use the median form since it is
    the one the derivation produces.)

Sawin's Question 2, as formalised by Yu and Cambie: find the largest c such
that for some alpha in [0,1], for EVERY family of identically distributed
[0,1]-valued p,q,r with E[p] <= c, q independent of p, and (p,r) an ARBITRARY
coupling,

    (1-alpha) E[h(p+q-pq)] + alpha E[h(s*(p,r))]  >=  E[h(p)].          (Q2)

CLOSED FORM FOR THE CEILING.  Take p supported on {1, b} with Pr[p=1] = a,
and let (p,r) be maximally negatively correlated, Pr[p=r=1] = 0 (possible iff
a <= 1/2).  Then, using h(1) = 0 and s*(1,x) = 1:

    E[h(p)]            = (1-a) h(b)
    E[h(p+q-pq)]       = (1-a)^2 h(2b-b^2)
    E[h(s*(p,r))]      = (1-2a) h(1/2) = 1-2a          [when 2b > 1/2]

The obstruction is the case where all three coincide, since then (Q2) holds
with equality for EVERY alpha and no choice of alpha helps:

    (1-a) h(2b-b^2) = h(b)          =>  1-a = h(b)/h(2b-b^2)
    1-2a = (1-a) h(b)               =>  1-a = 1/(2-h(b))

Eliminating a gives the defining equation for b, and then c:

    h(b) (2 - h(b)) = h(2b - b^2)                                        (*)
    c* = a + (1-a) b = 1 - (1-a)(1-b) = 1 - (1-b)/(2-h(b))               (**)

Everything below is checked numerically against Cambie's reported
0.382345533366703.
"""

import os
import sys

import numpy as np
from mpmath import findroot
from scipy.optimize import linprog, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI, h, mp, mpf

# ---------------------------------------------------------------------------
# 1.  The ceiling in closed form
# ---------------------------------------------------------------------------
print("1. Cambie's ceiling for the two-strategy bound, from equation (*)")


def cambie_eq(x):
    """h(x)(2-h(x)) - h(2x-x^2); its larger root in (0,1/2) defines b."""
    hx = h(x)
    return hx * (2 - hx) - h(2 * x - x * x)


# Two roots are reported in (0, 1/2); bracket each.
b1 = findroot(cambie_eq, mpf("0.14"))
b2 = findroot(cambie_eq, mpf("0.33"))
print(f"   roots of h(x)(2-h(x)) = h(2x-x^2):")
print(f"     b1 = {mp.nstr(b1, 18)}")
print(f"     b2 = {mp.nstr(b2, 18)}")
assert abs(cambie_eq(b1)) < mpf(10) ** -30
assert abs(cambie_eq(b2)) < mpf(10) ** -30

b = b2
a = 1 - 1 / (2 - h(b))
c_star = 1 - (1 - b) / (2 - h(b))
print(f"   a  = {mp.nstr(a, 18)}   (needs a <= 1/2: {a <= mpf('0.5')})")
print(f"   c* = 1 - (1-b)/(2-h(b)) = {mp.nstr(c_star, 18)}")

# Consistency: all three entropy terms must agree at the obstruction.
term_solo = (1 - a) * h(b)
term_iid = (1 - a) ** 2 * h(2 * b - b * b)
term_cpl = 1 - 2 * a
print(f"   E h(p)          = {mp.nstr(term_solo, 16)}")
print(f"   E h(p+q-pq)     = {mp.nstr(term_iid, 16)}")
print(f"   E h(s*(p,r))    = {mp.nstr(term_cpl, 16)}")
assert abs(term_iid - term_solo) < mpf(10) ** -30, term_iid - term_solo
assert abs(term_cpl - term_solo) < mpf(10) ** -30, term_cpl - term_solo
print("   all three equal => no alpha can help                      [OK]")

CAMBIE = mpf("0.382345533366703")
assert abs(c_star - CAMBIE) < mpf(10) ** -15, c_star - CAMBIE
print(f"   matches Cambie's 0.382345533366703 to 1e-15               [OK]")
assert c_star > PSI
print(f"   c* - psi = {mp.nstr(c_star - PSI, 6)} > 0                        [OK]")

# ---------------------------------------------------------------------------
# 2.  Independent search: is that obstruction really the worst case?
# ---------------------------------------------------------------------------
print("\n2. Search over k-atom laws and ALL couplings (inner LP is exact)")

LOG2 = np.log(2.0)


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def maxent_or(pv, rv):
    """s*(p,r) = median{1/2, max(p,r), min(p+r,1)}, vectorised."""
    lo = np.maximum(pv, rv)
    hi = np.minimum(pv + rv, 1.0)
    return np.clip(0.5, lo, hi)


def worst_coupling(atoms, w, alpha):
    """Exact inner minimisation over couplings of w with itself, by LP.

    The coupled term is linear in the joint matrix M, and the feasible set
    {M >= 0, M 1 = w, M^T 1 = w} is the transportation polytope, so the
    minimum is an LP and needs no heuristic search.
    """
    k = len(atoms)
    S = maxent_or(atoms[:, None], atoms[None, :])
    cost = (alpha * hf(S.ravel())).astype(float)
    # Row sums and column sums both equal w.
    Aeq, beq = [], []
    for i in range(k):
        row = np.zeros((k, k)); row[i, :] = 1.0
        Aeq.append(row.ravel()); beq.append(w[i])
    for j in range(k):
        col = np.zeros((k, k)); col[:, j] = 1.0
        Aeq.append(col.ravel()); beq.append(w[j])
    res = linprog(cost, A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * (k * k), method="highs")
    if not res.success:
        return np.inf
    return float(res.fun)


def slack(atoms, w, alpha):
    """(1-alpha) E h(p+q-pq) + alpha min_M E h(s*) - E h(p).  Negative = broken."""
    atoms = np.asarray(atoms, dtype=float)
    w = np.asarray(w, dtype=float)
    orm = atoms[:, None] + atoms[None, :] - atoms[:, None] * atoms[None, :]
    iid = w @ hf(orm.ravel()).reshape(orm.shape) @ w
    cpl = worst_coupling(atoms, w, 1.0)
    return (1 - alpha) * iid + alpha * cpl - w @ hf(atoms)


def min_slack(t, alpha, k, restarts=120, seed=0):
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)

    def unpack(z):
        atoms = 1.0 / (1.0 + np.exp(-np.clip(z[:k], -30, 30)))
        e = np.exp(z[k:] - z[k:].max())
        return atoms, e / e.sum()

    def obj(z):
        return slack(*unpack(z), alpha)

    def mean_con(z):
        atoms, w = unpack(z)
        return t - float(w @ atoms)

    for _ in range(restarts):
        z0 = rng.normal(0, 2.0, 2 * k)
        try:
            r = minimize(obj, z0, method="SLSQP",
                         constraints=[{"type": "ineq", "fun": mean_con}],
                         options={"maxiter": 300, "ftol": 1e-13})
        except Exception:
            continue
        if r.success and r.fun < best[0]:
            atoms, w = unpack(r.x)
            best = (float(r.fun), atoms, w)
    return best


ALPHA = 0.0356069        # Cambie's reported optimal alpha
cs = float(c_star)
print(f"   using alpha = {ALPHA} (Cambie's optimum)")
for t, lbl in [(cs - 2e-4, "c* - 2e-4"), (cs, "c*"), (cs + 2e-4, "c* + 2e-4")]:
    val, atoms, w = min_slack(t, ALPHA, 3, restarts=120, seed=7)
    keep = w > 1e-4
    desc = "  ".join(f"{wi:.4f}@{ai:.6f}" for ai, wi in
                     sorted(zip(atoms[keep], w[keep]), key=lambda s: -s[1]))
    print(f"   t = {lbl:>10}: min slack = {val:+.3e}   {desc}")

# Reproduce the obstruction law directly and confirm zero slack for every alpha.
print("\n   direct check of the closed-form obstruction law {1, b}:")
oa, ob = float(a), float(b)
atoms_obs = np.array([1.0, ob])
w_obs = np.array([oa, 1 - oa])
print(f"     law: {oa:.6f}@1.0  {1-oa:.6f}@{ob:.6f},  mean = "
      f"{float(w_obs @ atoms_obs):.15f}")
for al in (0.0, 0.0356069, 0.2, 0.5, 1.0):
    s = slack(atoms_obs, w_obs, al)
    print(f"     alpha = {al:<9}: slack = {s:+.3e}")
    assert abs(s) < 1e-9, (al, s)
print("     slack is zero for EVERY alpha                           [OK]")

print("\nALL PASS")
