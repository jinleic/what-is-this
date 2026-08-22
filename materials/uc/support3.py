"""Theorem B' in action: the two-strategy problem on <= 3 PAIR-ORBITS.

Theorem B (decomposition.py) claimed <= 3 MARGINAL atoms and is RETRACTED: it
needed C to be concave, and C is convex.  The surviving reduction keeps the
symmetric coupling itself as the variable.  Write

    Delta = { (p,q) : 0 <= p <= q <= 1 },      c(p,q) = h(s*(p,q)),
    mu_nu = int (delta_p + delta_q)/2 dnu      (the induced marginal),
    Phi(nu) = (1-alpha) Q(mu_nu) + alpha int c dnu - L(mu_nu).

Then inf F over feasible mu equals inf Phi over feasible nu, and the latter is
attained on nu with at most 3 atoms -- at most 3 unordered pairs, hence at most
6 marginal atoms.  8 parameters, not 5; weaker than the retracted claim, and in
turn SUPERSEDED by Theorem B''' (margin_lemma.py), which reaches 2 orbits and 5
parameters by freezing B = 1-mean.  This file's searches remain valid evidence.

SCALE-FREE FORM.  Phi = 0 identically on the SINK FAMILY nu_u carried by the
pairs {0,0} and {0,1} and {1,1}, whose marginal is (1-u)delta_0 + u delta_1,
because h(0) = h(1) = 0, 0 (+) 1 = 1 (+) 1 = 1, and s*(p,1) = 1.  So min Phi = 0
for every t and carries no information; the real object is the ratio

    R(nu) = [ (1-alpha) Q(mu_nu) + alpha int c dnu ] / L(mu_nu),

and the bound holds at t exactly when inf R >= 1.
"""

import os
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI, h, mp, mpf

LOG2 = np.log(2.0)
A_OBS = mpf("0.078877292705923173412")
B_OBS = mpf("0.32945473850303697239")
C_STAR = mpf("0.38234553336670272115")
ALPHA = 0.0356069


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def sstar(p, r):
    """s*(p,r) = median{1/2, max(p,r), min(p+r,1)}."""
    return np.clip(0.5, np.maximum(p, r), np.minimum(p + r, 1.0))


def coupled_min(atoms, w):
    """Exact min over symmetric couplings of w with itself, by LP."""
    k = len(atoms)
    P, R = np.meshgrid(atoms, atoms, indexing="ij")
    cost = hf(sstar(P, R)).ravel()
    rows, rhs = [], []
    for i in range(k):
        row = np.zeros(k * k)
        row[i * k:(i + 1) * k] = 1.0
        rows.append(row)
        rhs.append(w[i])
    for i in range(k):
        for j in range(i + 1, k):
            row = np.zeros(k * k)
            row[i * k + j] = 1.0
            row[j * k + i] = -1.0
            rows.append(row)
            rhs.append(0.0)
    res = linprog(cost, A_eq=np.array(rows), b_eq=np.array(rhs),
                  bounds=[(0, None)] * (k * k), method="highs")
    assert res.status == 0, res.message
    return float(res.fun)


def marginal(pairs, wts):
    """mu_nu as (atoms, weights), each pair splitting its mass in half."""
    at, wt = [], []
    for (p, q), wi in zip(pairs, wts):
        at += [float(p), float(q)]
        wt += [wi / 2.0, wi / 2.0]
    return np.array(at), np.array(wt)


def orbit_parts(pairs, wts):
    at, wt = marginal(pairs, wts)
    P, R = np.meshgrid(at, at, indexing="ij")
    Q = float(wt @ hf(P + R - P * R) @ wt)
    L = float(wt @ hf(at))
    Cc = float(sum(wi * hf(sstar(np.array([p]), np.array([q])))[0]
                   for (p, q), wi in zip(pairs, wts)))
    mean = float(sum(wi * (p + q) / 2.0 for (p, q), wi in zip(pairs, wts)))
    return Q, Cc, L, mean


def orbit_ratio(pairs, wts, alpha=ALPHA):
    Q, Cc, L, _ = orbit_parts(pairs, wts)
    return np.inf if L <= 1e-14 else ((1 - alpha) * Q + alpha * Cc) / L


def orbit_slack(pairs, wts, alpha=ALPHA):
    Q, Cc, L, _ = orbit_parts(pairs, wts)
    return (1 - alpha) * Q + alpha * Cc - L


# --- the exact 1-D obstruction family (marginal form) ----------------------
def family_ratio(b, t, alpha):
    """R for mu = a d_1 + (1-a) d_b with mean exactly t, in mpmath."""
    b, t, al = mpf(b), mpf(t), mpf(alpha)
    a = (t - b) / (1 - b)
    if a < 0 or a > mpf("0.5"):
        return None
    hb = h(b)
    if hb == 0:
        return None
    sbb = min(max(mpf("0.5"), b), min(2 * b, mpf(1)))
    Q = (1 - a) ** 2 * h(2 * b - b * b)
    L = (1 - a) * hb
    C = (1 - 2 * a) * h(sbb)
    return ((1 - al) * Q + al * C) / L


def family_min(t, alpha, N=3000):
    lo, hi = mpf("0.05"), mpf("0.499")
    best = (mp.inf, None)
    for i in range(N + 1):
        b = lo + (hi - lo) * mpf(i) / N
        v = family_ratio(b, t, alpha)
        if v is not None and v < best[0]:
            best = (v, b)
    b0, step = best[1], (hi - lo) / N
    for _ in range(120):
        for s in (-step, step):
            v = family_ratio(b0 + s, t, alpha)
            if v is not None and v < best[0]:
                best, b0 = (v, b0 + s), b0 + s
        step /= 2
    return best


def min_orbit_ratio(t, korb, alpha=ALPHA, restarts=60, seed=0):
    """Minimise R over nu with korb pair-orbits, mean(mu_nu) <= t."""
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)

    def project(raw, w):
        w = np.clip(w, 0.0, None)
        s = w.sum()
        if s <= 0:
            return None
        w = w / s
        pr = [(min(u, v), max(u, v)) for u, v in raw]
        mean = sum(w[i] * (pr[i][0] + pr[i][1]) / 2.0 for i in range(len(pr)))
        if mean > t:
            sc = t / mean
            pr = [(p * sc, q * sc) for p, q in pr]
        return pr, w

    seeds = []
    fb = family_min(t, alpha)[1]
    if fb is not None:
        bb, av = float(fb), float((mpf(t) - fb) / (1 - fb))
        base = [(bb, bb), (bb, 1.0)] + [(bb, bb)] * (korb - 2)
        bw = [1 - 2 * av, 2 * av] + [0.0] * (korb - 2)
        if korb >= 2:
            seeds.append((base[:korb], np.array(bw[:korb])))
    for _ in range(restarts):
        seeds.append((list(rng.uniform(0, 1, (korb, 2))),
                      rng.dirichlet(np.ones(korb))))

    for raw, w in seeds:
        pr = project(list(raw), np.array(w, dtype=float))
        if pr is None:
            continue
        prs, w = pr
        cur, step = orbit_ratio(prs, w, alpha), 0.25
        for _ in range(500):
            jr = [(min(max(p + rng.normal(0, step), 0.0), 1.0),
                   min(max(q + rng.normal(0, step), 0.0), 1.0)) for p, q in prs]
            pj = project(jr, w + rng.normal(0, step * 0.5, korb))
            if pj is None:
                continue
            jp, jw = pj
            val = orbit_ratio(jp, jw, alpha)
            if val < cur:
                prs, w, cur = jp, jw, val
            else:
                step *= 0.985
            if step < 1e-10:
                break
        if cur < best[0]:
            best = (cur, list(prs), np.array(w))
    return best


# ---------------------------------------------------------------------------
print("1. inf Phi over nu equals inf F over mu -- the reformulation is exact")
print("   for each mu, the OPTIMAL coupling gives a nu with Phi(nu) = F(mu)")
print(f"   {'law':>34} {'F(mu)':>16} {'Phi(nu_opt)':>16}")
for atoms, wts in (([0.3, 0.6], [0.5, 0.5]), ([1.0, float(B_OBS)],
                                              [float(A_OBS), float(1 - A_OBS)]),
                   ([0.1, 0.35, 0.9], [0.2, 0.5, 0.3])):
    at, wv = np.array(atoms), np.array(wts)
    P, R = np.meshgrid(at, at, indexing="ij")
    Q = float(wv @ hf(P + R - P * R) @ wv)
    L = float(wv @ hf(at))
    F = (1 - ALPHA) * Q + ALPHA * coupled_min(at, wv) - L
    # rebuild the optimal coupling as a pair-orbit measure
    k = len(at)
    cost = hf(sstar(P, R)).ravel()
    rows, rhs = [], []
    for i in range(k):
        row = np.zeros(k * k)
        row[i * k:(i + 1) * k] = 1.0
        rows.append(row)
        rhs.append(wv[i])
    for i in range(k):
        for j in range(i + 1, k):
            row = np.zeros(k * k)
            row[i * k + j] = 1.0
            row[j * k + i] = -1.0
            rows.append(row)
            rhs.append(0.0)
    sol = linprog(cost, A_eq=np.array(rows), b_eq=np.array(rhs),
                  bounds=[(0, None)] * (k * k), method="highs").x.reshape(k, k)
    prs, pw = [], []
    for i in range(k):
        for j in range(i, k):
            mass = sol[i, j] + (sol[j, i] if j > i else 0.0)
            if mass > 1e-12:
                prs.append((min(at[i], at[j]), max(at[i], at[j])))
                pw.append(mass)
    pw = np.array(pw)
    ph = orbit_slack(prs, pw)
    print(f"   {str(atoms)[:32]:>34} {F:>16.12f} {ph:>16.12f}")
    assert abs(F - ph) < 1e-9, (atoms, F, ph)
print("   agree, so minimising over pair orbits is lossless             [OK]")

# ---------------------------------------------------------------------------
print("\n2. The sink family: Phi = 0 there exactly, for every u")
for u in (0.0, 0.1, 0.25, 0.38):
    prs = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
    pw = np.array([(1 - u) ** 2, 2 * u * (1 - u), u**2])
    Q, Cc, L, mean = orbit_parts(prs, pw)
    print(f"   u = {u:<5}: Q = {Q:.3e}  int c = {Cc:.3e}  L = {L:.3e}  "
          f"mean = {mean:.4f}  Phi = {orbit_slack(prs, pw):+.3e}")
    assert abs(orbit_slack(prs, pw)) < 1e-12
print("   => min Phi = 0 for every t; the ratio R is the real object    [OK]")

print("\n   R near the sink, mu = (1-u-e) d_0 + u d_1 + e d_0.3:")
for u in (0.0, 0.2, 0.38):
    for e in (1e-3, 1e-5):
        at = np.array([0.0, 1.0, 0.3])
        wv = np.array([1 - u - e, u, e])
        P, R = np.meshgrid(at, at, indexing="ij")
        Q = float(wv @ hf(P + R - P * R) @ wv)
        L = float(wv @ hf(at))
        r = ((1 - ALPHA) * Q + ALPHA * coupled_min(at, wv)) / L
        pred = 2 * (1 - ALPHA) * (1 - u)
        print(f"     u = {u:<5} e = {e:<7}: R = {r:.6f}   predicted {pred:.6f}")
        assert abs(r - pred) < 0.05
print("   margin at u = t is ~0.19, not a near-violation                [OK]")

# ---------------------------------------------------------------------------
print("\n3. The obstruction, as a 2-pair-orbit measure, reaches R = 1")
prs = [(float(B_OBS), float(B_OBS)), (float(B_OBS), 1.0)]
pw = np.array([float(1 - 2 * A_OBS), float(2 * A_OBS)])
Q, Cc, L, mean = orbit_parts(prs, pw)
print(f"   nu*: {{b,b}} @ {pw[0]:.12f}   {{b,1}} @ {pw[1]:.12f}")
print(f"   b = {float(B_OBS):.18f}")
print(f"   mean(mu_nu) = {mean:.18f}")
print(f"   c*          = {float(C_STAR):.18f}")
assert abs(mean - float(C_STAR)) < 1e-15
print(f"   Q = {Q:.15f}   int c = {Cc:.15f}   L = {L:.15f}")
assert abs(Q - L) < 1e-12 and abs(Cc - L) < 1e-12
print("   all three coincide, so no alpha can help:")
for al in (0.0, ALPHA, 0.2, 0.5, 1.0):
    r = orbit_ratio(prs, pw, al)
    print(f"     alpha = {al:<8}: R = {r:.15f}")
    assert abs(r - 1) < 1e-11
print("   c* is the ceiling of the whole two-strategy family            [OK]")
print("   and it uses only 2 pair-orbits, inside B'                     [OK]")

# ---------------------------------------------------------------------------
print("\n4. Search over k pair-orbits (evidence, NOT a proof of optimality)")
print(f"   alpha = {ALPHA},  t = c*")
print(f"   {'orbits':>7} {'min R':>18}")
rows = {}
for korb in (2, 3, 4, 5):
    val, prs_b, pw_b = min_orbit_ratio(float(C_STAR), korb, restarts=40,
                                       seed=5 + korb)
    rows[korb] = val
    print(f"   {korb:>7} {val:>18.12f}")
best23 = min(rows[2], rows[3])
for korb in (4, 5):
    assert rows[korb] > best23 - 1e-9, (korb, rows[korb], best23)
assert abs(best23 - 1) < 1e-8, best23
print("   no k > 3 beats k <= 3, and k <= 3 attains exactly 1 at c*     [OK]")
print("   Theorem B' guarantees k <= 3 suffices; the search only shows")
print("   the reduction is not vacuous at the binding instance.")

# ---------------------------------------------------------------------------
print("\n5. Room below c*: exact 1-D family vs the pair-orbit search")
print(f"   {'t':>22} {'t - psi':>11} {'family min R':>18} {'orbit search':>15}")
grid = [PSI, PSI + mpf(10) ** -5, PSI + mpf(10) ** -4,
        PSI + 2 * mpf(10) ** -4, PSI + 3 * mpf(10) ** -4, C_STAR]
margins = {}
for t in grid:
    fv = family_min(t, ALPHA)[0]
    sv = min(min_orbit_ratio(float(t), 2, restarts=25, seed=17)[0],
             min_orbit_ratio(float(t), 3, restarts=25, seed=19)[0])
    margins[t] = fv - 1
    print(f"   {float(t):>22.16f} {float(t - PSI):>11.3e} {float(fv):>18.12f} "
          f"{sv:>15.12f}")
    assert sv > float(fv) - 1e-7, (t, sv, fv)
print("   the search never beats the 1-D family                         [OK]")

print("\n   margin of the 1-D family is essentially linear in (c* - t):")
print(f"   {'t':>22} {'c* - t':>12} {'margin':>13} {'ratio':>9}")
for t in grid[:-1]:
    gapv = float(C_STAR - t)
    print(f"   {float(t):>22.16f} {gapv:>12.3e} {float(margins[t]):>13.3e} "
          f"{float(margins[t]) / gapv:>9.4f}")

print("\nALL PASS")
print("\nSUMMARY")
print(f"  psi = {float(PSI):.16f}   (record, four independent proofs)")
print(f"  c*  = {float(C_STAR):.16f}   (claimed; Yu's reduction invalid)")
print("  Theorem B' reduces the remaining task to <= 3 pair-orbits:")
print("  Theorem B''' (margin_lemma.py) supersedes this: <= 2 pair-orbits,")
print("  5 parameters, by freezing B = 1-mean rather than G.")
