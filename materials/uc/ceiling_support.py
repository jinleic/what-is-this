"""SUPPORT-RESTRICTED CEILING -- A TEMPTING ROUTE PAST c*, AND WHY IT FAILS.

STATUS: the motivating hypothesis of this script is REFUTED.  The numbers it
prints are correct; the interpretation that would have made them a gain over
c* is not.  Kept as a documented negative result, because the failure mode is
precise and rules out a whole family of "just restrict the support" ideas.

WHAT WAS ASKED.  coupling.py machine-checks the ceiling of the
mean-constrained problem: the law p in {1, b} with
    1 - a = h(b)/h(2b-b^2) = 1/(2-h(b)),   c* = 1 - (1-b)/(2-h(b))
makes ALL THREE entropy terms coincide, so (Q2) holds with equality for EVERY
alpha -- no choice of alpha can push past c* = 0.3823455333667027.  That
obstruction places an atom at p = 1.  A family with an element of frequency 1
satisfies the conjecture outright, so it is tempting to argue that the
functional inequality only ever needs to hold for laws supported in [0,c),
which would exclude the obstruction and raise the ceiling.  Section 3 below
measures exactly that: the sampled support-restricted ceiling is 0.4003611
at alpha = 0.5, i.e. +1.8e-2 above c*.

WHY IT FAILS.  In the Gilmer/Sawin/Cambie argument the random variable in the
functional is NOT the element frequency.  With A, B iid uniform on the family
and the chain rule taken coordinatewise,
    H(A) = sum_i E[h(a_i)],   a_i := Pr[A_i = 1 | A_{<i}],
and the union bound compares this against E[h(a_i + b_i - a_i b_i)] with
(a_i, b_i) conditionally independent.  The contradiction hypothesis "every
element has frequency < c" bounds the UNCONDITIONAL frequency
    Pr[i in A] = E[a_i] < c,
i.e. it bounds the MEAN of the conditional-probability law -- not its
support.  A coordinate can be forced (a_i = 1) on some prefixes while its
frequency stays below c, so the p = 1 atom is genuinely admissible and the
mean constraint E[p] <= c is exactly the right formulation.  c* stands.

CONSEQUENCE FOR THE PROGRAMME.  Nothing about the certified results changes:
our certificates are proved under the (weaker, correct) mean constraint, so
they are unaffected -- a mean-constrained certificate implies the
support-restricted one, never the reverse.  What this rules out is the cheap
route past c*: any argument that hopes to gain by shrinking the feasible set
to laws supported below c is invalid at the source.  Going past c* needs a
genuinely stronger inequality (a third strategy, or Liu-style auxiliary
randomness whose hypotheses we would have to certify), not a tighter domain.

EVERYTHING BELOW IS NUMERICAL: sampled minima over <=3-atom laws, reported to
quantify the (unusable) gap.  No claim in this file is machine-checked.

Run: ./.venv/bin/python uc/ceiling_support.py
"""
import os
import sys

import numpy as np
from scipy.optimize import linprog, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI  # noqa: E402

LOG2 = np.log(2.0)
C_STAR = 0.3823455333667027          # coupling.py, mean-constrained ceiling
PSI_F = float(PSI)


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def maxent_or(pv, rv):
    """s*(p,r) = median{1/2, max(p,r), min(p+r,1)} (Cambie's max-entropy OR)."""
    lo = np.maximum(pv, rv)
    hi = np.minimum(pv + rv, 1.0)
    return np.clip(0.5, lo, hi)


def worst_coupling(atoms, w):
    """min over couplings M of sum M_ij h(s*(a_i,a_j)); exact LP."""
    k = len(atoms)
    S = maxent_or(atoms[:, None], atoms[None, :])
    cost = hf(S.ravel()).astype(float)
    Aeq, beq = [], []
    for i in range(k):
        row = np.zeros((k, k))
        row[i, :] = 1.0
        Aeq.append(row.ravel())
        beq.append(w[i])
    for j in range(k):
        col = np.zeros((k, k))
        col[:, j] = 1.0
        Aeq.append(col.ravel())
        beq.append(w[j])
    res = linprog(cost, A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * (k * k), method="highs")
    return float(res.fun) if res.success else np.inf


def slack(atoms, w, alpha):
    """(1-a) E h(p+q-pq) + a min_M E h(s*) - E h(p).  Negative = broken."""
    atoms = np.asarray(atoms, dtype=float)
    w = np.asarray(w, dtype=float)
    orm = atoms[:, None] + atoms[None, :] - atoms[:, None] * atoms[None, :]
    iid = w @ hf(orm.ravel()).reshape(orm.shape) @ w
    cpl = worst_coupling(atoms, w)
    return (1 - alpha) * iid + alpha * cpl - w @ hf(atoms)


def min_slack_support(c, alpha, k, restarts=40, seed=0):
    """Minimise slack over laws with k atoms, ALL atoms in [0,c].

    The support restriction implies mean <= c, so no separate mean constraint
    is needed: this is a strictly smaller feasible set than the
    mean-constrained problem coupling.py explores.
    """
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)

    def unpack(z):
        atoms = np.clip(z[:k], 0.0, c)
        raw = np.abs(z[k:]) + 1e-12
        return atoms, raw / raw.sum()

    def obj(z):
        return slack(*unpack(z), alpha)

    for _ in range(restarts):
        z0 = np.concatenate([rng.uniform(0.0, c, k), rng.uniform(0.1, 1.0, k)])
        res = minimize(obj, z0, method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-10,
                                "fatol": 1e-12})
        if res.fun < best[0]:
            atoms, w = unpack(res.x)
            best = (float(res.fun), atoms, w)
    return best


def ceiling_for_alpha(alpha, kmax=3, lo=0.37, hi=0.60, iters=26, seed=1):
    """Largest c (bisection) with sampled min slack >= 0 over support [0,c]."""
    def feasible(c):
        for k in (1, 2, kmax):
            val = min_slack_support(c, alpha, k, restarts=18, seed=seed + k)[0]
            if val < -1e-9:
                return False, val
        return True, None

    ok_lo, _ = feasible(lo)
    if not ok_lo:
        return lo, "infeasible already at lo"
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ok, _ = feasible(mid)
        if ok:
            lo = mid
        else:
            hi = mid
    return lo, None


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    print("psi    = %.16f" % PSI_F)
    print("c*     = %.16f   (mean-constrained ceiling, alpha-independent)"
          % C_STAR)
    print()

    # --- 1. the mean-constrained obstruction is NOT support-admissible ------
    print("1. The c* obstruction has an atom at p = 1 (excluded by support)")
    B_OBS = 0.32945473850303697239        # pinned in cert3_par.py / margin_lemma.py
    a_obs = 1 - 1 / (2 - float(hf(np.array([B_OBS]))[0]))
    atoms_obs = np.array([1.0, B_OBS])
    w_obs = np.array([a_obs, 1 - a_obs])
    print("   obstruction law: %.6f @ p=1   +   %.6f @ p=%.17f"
          % (a_obs, 1 - a_obs, B_OBS))
    print("   mean = %.16f  (<= c* : %s)"
          % (w_obs @ atoms_obs, bool(w_obs @ atoms_obs <= C_STAR + 1e-12)))
    for al in (0.0, 0.0356069, 0.25, 0.5, 1.0):
        print("     alpha = %-9s slack = %+.3e" % (al, slack(atoms_obs, w_obs, al)))
    print("   => zero for every alpha (coupling.py's finding), and it needs the")
    print("      p=1 atom, which the support hypothesis forbids.")
    print()
    print("2. Point mass at c (support-admissible): slack as a function of c")
    print("   slack(c,alpha) = (1-alpha) h(2c-c^2) + alpha h(s*(c,c)) - h(c)")
    print("   %-10s %-14s %-14s %-14s" % ("c", "alpha=0.0356", "alpha=0.25",
                                          "alpha=0.50"))
    for c in (0.3823, 0.39, 0.40, 0.42, 0.45, 0.48):
        vals = []
        for al in (0.0356069, 0.25, 0.50):
            vals.append(slack(np.array([c]), np.array([1.0]), al))
        print("   %-10.4f %+14.6e %+14.6e %+14.6e" % (c, *vals))
    print()

    # --- 3. sampled support-restricted ceiling per alpha -------------------
    print("3. Sampled support-restricted ceiling c_sup(alpha)  [NUMERICAL]")
    print("   %-10s %-14s %s" % ("alpha", "c_sup(alpha)", "vs c*"))
    best = (0.0, None)
    for alpha in (0.0356069, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60):
        c_sup, note = ceiling_for_alpha(alpha)
        flag = "ABOVE c*" if c_sup > C_STAR + 1e-6 else "at/below c*"
        print("   %-10.6f %-14.10f %s%s" % (alpha, c_sup, flag,
                                            "" if note is None else " (%s)" % note))
        if c_sup > best[0]:
            best = (c_sup, alpha)
    print()
    print("   best sampled ceiling %.10f at alpha = %s" % (best[0], best[1]))
    print("   gain over c*: %+.3e" % (best[0] - C_STAR))
    print()
    print("NUMERICAL ONLY -- sampled minima over <=3-atom laws; not a proof.")
    print()
    print("VERDICT: this gain is NOT AVAILABLE.  The variable in the")
    print("functional is the CONDITIONAL probability Pr[A_i=1 | A_{<i}], whose")
    print("MEAN is the element frequency; the contradiction hypothesis bounds")
    print("that mean, not the support, so the p=1 atom stays admissible and")
    print("c* = %.16f remains the two-strategy ceiling." % C_STAR)
    print("Certificates proved under the mean constraint (ours) are unaffected;")
    print("they imply the support-restricted statement, never the reverse.")
    print("Going above c* needs a strictly stronger inequality -- e.g. Liu's")
    print("auxiliary-randomness route, whose two numerical hypotheses are")
    print("exactly the kind of object this project's interval machinery can")
    print("discharge (see PROGRESS.md, next-target analysis).")
