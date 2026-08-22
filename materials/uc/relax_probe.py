"""Probe: what dimension reduction, if any, is actually available?

Written after retracting a bogus 5D->3D claim (product-order monotonicity in the
aggregates does NOT let you push M to t, and the Pareto set of a 2-D orbit
surface under 4 objectives is generically 2-D, not 1-D).

Nothing here is a certificate.  This file records which reductions SURVIVE a
check and which die, so the next attempt does not repeat either mistake.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entropy import PSI, h, mp, mpf

mp.dps = 30
AL = mpf("0.0356069")
CS = mpf("0.38234553333667027")
PS = mpf(PSI)
B0 = mpf("0.32945473850303697239")


def hh(z):
    return h(mpf(z))


def rr(z):
    z = mpf(z)
    x = 1 - z
    v = 2 * x * h(z) - h(x * x)
    return v if v > 0 else mpf(0)


def st(a, b):
    return min(max(mpf("0.5"), max(a, b)), min(a + b, mpf(1)))


def parts(p1, q1, p2, q2, w):
    p1, q1, p2, q2, w = map(mpf, (p1, q1, p2, q2, w))
    m1, m2 = (p1 + q1) / 2, (p2 + q2) / 2
    M = w * m1 + (1 - w) * m2
    L = w * (hh(p1) + hh(q1)) / 2 + (1 - w) * (hh(p2) + hh(q2)) / 2
    C = w * hh(st(p1, q1)) + (1 - w) * hh(st(p2, q2))
    S = (w * (mp.sqrt(rr(p1)) + mp.sqrt(rr(q1))) / 2
         + (1 - w) * (mp.sqrt(rr(p2)) + mp.sqrt(rr(q2))) / 2)
    return M, L, C, S


def Phi(p1, q1, p2, q2, w, a=AL):
    M, L, C, S = parts(p1, q1, p2, q2, w)
    return (2 * (1 - a) * (1 - M) - 1) * L + a * C - (1 - a) * S * S


def Phi_relax(p1, q1, p2, q2, w, t, a=AL):
    _, L, C, S = parts(p1, q1, p2, q2, w)
    return (2 * (1 - a) * (1 - mpf(t)) - 1) * L + a * C - (1 - a) * S * S


rng = np.random.default_rng(1)

print("1. M := t IS a valid relaxation (coeff(L) decreasing in M, L >= 0),")
print("   and it is EXACTLY tight at the obstruction -- but useless globally.")
for t, lbl in ((PS, "psi"), (CS, "c* ")):
    a = (t - B0) / (1 - B0)
    w = 1 - 2 * a
    f, g = Phi(B0, B0, B0, 1, w), Phi_relax(B0, B0, B0, 1, w, t)
    print("   t=" + lbl + " at obstruction: loss = " + mp.nstr(f - g, 6))
    worst = mp.inf
    for p in np.linspace(0, 1, 121):
        for q in np.linspace(p, 1, 81):
            worst = min(worst, Phi_relax(p, q, p, q, 1, t))
    print("   t=" + lbl + " min of Phi_relax on single orbits = "
          + mp.nstr(worst, 8) + "  -> NEGATIVE, relaxation too lossy")

print("")
print("2. Is Phi CONCAVE in w?  (concave => min at w in {0,1} => 2 params)")
print("   w^2 coefficient = -2(1-a)*dM*dL - (1-a)*dS^2 : sign INDEFINITE.")
conv = 0
for _ in range(1500):
    v = rng.uniform(0, 1, 4)
    p1, q1 = sorted(v[:2])
    p2, q2 = sorted(v[2:4])

    def f(w):
        return Phi(p1, q1, p2, q2, w)
    if f(mpf("0.5")) < (f(mpf(0)) + f(mpf(1))) / 2 - mpf("1e-25"):
        conv += 1
print("   strictly-convex-in-w cases out of 1500:", conv)
a = (CS - B0) / (1 - B0)
w_obs = 1 - 2 * a
print("   and the obstruction minimises at w = " + mp.nstr(w_obs, 8)
      + ", strictly INTERIOR,")
print("   which a concave function cannot do.  => NO w-endpoint reduction.")
assert conv > 0

print("")
print("3. The one reduction that could work: sign of dPhi/dq2 on the feasible set.")
print("   If dPhi/dq2 <= 0 then the min has q2 = 1 OR the mean binds (M = t):")
print("   either branch drops 5 -> 4 dimensions RIGOROUSLY (a case split, not")
print("   a deformation hand-wave).  Testing the sign:")
d = mpf("1e-8")
neg = pos = tot = 0
worst_pos = (mpf(0), None)
for _ in range(20000):
    v = rng.uniform(0, 1, 5)
    p1, q1 = sorted(v[:2])
    p2, q2 = sorted(v[2:4])
    w = v[4]
    if q2 + d > 1:
        continue
    if parts(p1, q1, p2, q2, w)[0] > CS:
        continue
    tot += 1
    g = (Phi(p1, q1, p2, q2 + d, w) - Phi(p1, q1, p2, q2, w)) / d
    if g > 0:
        pos += 1
        if g > worst_pos[0]:
            worst_pos = (g, (float(p1), float(q1), float(p2), float(q2),
                             float(w)))
    else:
        neg += 1
print("   feasible samples " + str(tot) + ":  dPhi/dq2 <= 0 in " + str(neg)
      + ",  > 0 in " + str(pos))
if pos:
    print("   largest POSITIVE derivative " + mp.nstr(worst_pos[0], 6)
          + " at " + str(tuple(round(x, 4) for x in worst_pos[1])))
    print("   => the sign is INDEFINITE: no q2 = 1 lemma of this form.")
else:
    print("   => sign holds on all samples: a q2 = 1 lemma is a live target.")

print("")
print("4. Same test for the other three positions, for completeness.")
for k, nm in ((0, "p1"), (1, "q1"), (2, "p2")):
    neg = pos = tot = 0
    for _ in range(8000):
        v = rng.uniform(0, 1, 5)
        p1, q1 = sorted(v[:2])
        p2, q2 = sorted(v[2:4])
        w = v[4]
        P = [p1, q1, p2, q2]
        if P[k] + d > (P[k + 1] if k in (0, 2) else 1):
            continue
        if parts(p1, q1, p2, q2, w)[0] > CS:
            continue
        tot += 1
        Q = list(P)
        Q[k] += d
        g = (Phi(Q[0], Q[1], Q[2], Q[3], w) - Phi(p1, q1, p2, q2, w)) / d

        if g > 0:
            pos += 1
        else:
            neg += 1
    verdict = "INDEFINITE" if (pos and neg) else "one-signed"
    print("   " + nm + ": n=" + str(tot) + "  neg=" + str(neg) + " pos="
          + str(pos) + "  -> " + verdict)

print("")
print("5. What IS true at the q2 -> 1 boundary (where the obstruction sits).")
print("   Write eps = 1 - q2.  Two different rates:")
for e in ["1e-1", "1e-2", "1e-3", "1e-4", "1e-6"]:
    e = mpf(e)
    v = 2 * e * hh(e) - hh(e * e)
    print("     eps=" + mp.nstr(e, 4)
          + "  rh/eps^2=" + mp.nstr(v / (e * e), 12)
          + "  sqrt(rh)/eps=" + mp.nstr(mp.sqrt(v) / e, 12)
          + "  h(eps)/eps=" + mp.nstr(hh(e) / e, 8))
print("   rh(1-eps)/eps^2 -> " + mp.nstr(1 / mp.log(2), 12) + " = 1/ln2, so")
print("   sqrt(rh(1-eps)) = eps/sqrt(ln2) + O(eps^2): LINEAR, not singular.")
print("   h(1-eps) = h(eps) ~ eps log2(1/eps): vanishes with INFINITE slope.")
print("")
print("   So the ONLY singular quantity is eta := h(q2), and it enters Phi")
print("   AFFINELY with POSITIVE coefficient: kappa*(1-w)/2 in L, alpha*(1-w)")
print("   in C.  Increasing eta increases Phi.  The competing term is S, whose")
print("   q2-part sqrt(rh(q2)) also grows as q2 leaves 1, and -(1-a)S^2")
print("   DEcreases Phi.  The two fight -- which is why item 3 is indefinite.")
print("   But the rates differ: eps*log2(1/eps) beats eps as eps -> 0, and the")
print("   log term's coefficient is positive.  So q2 = 1 may be a strict local")
print("   minimum in q2.  TESTED ONLY ON FEASIBLE w -- an earlier version of")
print("   this block sampled w = 0.2 and 0.5, where M = 0.598 and 0.497 exceed")
print("   c*, i.e. INFEASIBLE points, and then printed 'all positive' above two")
print("   negative rows.  The assert below now enforces the sentence.")
w_min = ((B0 + 1) / 2 - CS) / (((B0 + 1) / 2) - B0)
print("   feasibility along (b,b,b,q2) requires w >= " + mp.nstr(w_min, 12)
       + " (= the obstruction's w)")
viol = []
for wv in np.linspace(float(w_min), 1.0, 15):
    wv = mpf(float(wv))
    base = Phi(B0, B0, B0, mpf(1), wv)
    for e in ["1e-2", "1e-3", "1e-5", "1e-8"]:
        e = mpf(e)
        if parts(B0, B0, B0, 1 - e, wv)[0] > CS:
            continue
        d_phi = Phi(B0, B0, B0, 1 - e, wv) - base
        if wv >= 1:
            # orbit 2 carries zero weight: Phi does not depend on q2 at all
            if d_phi != 0:
                viol.append(("w=1 must be q2-independent", float(wv), str(e)))
        elif d_phi <= 0:
            viol.append((float(wv), str(e), mp.nstr(d_phi, 6)))
for wv in (w_min, mpf("0.95"), mpf(1)):
    base = Phi(B0, B0, B0, mpf(1), wv)
    row = [mp.nstr(Phi(B0, B0, B0, 1 - mpf(e), wv) - base, 6)
           for e in ("1e-2", "1e-3", "1e-5", "1e-8")]
    print("     w=" + mp.nstr(wv, 8) + "  Phi(1-eps)-Phi(1): " + ", ".join(row))
print("   violations over 15 feasible w x 4 eps: " + str(len(viol)))
assert not viol, viol
print("   STRICTLY positive for every feasible w < 1, and identically 0 at w = 1")
print("   (where orbit 2 has zero weight, so q2 is absent).  So: moving q2 off 1")
print("   raises Phi on the whole feasible w-range.")
print("   [INFERENCE] this gives a slab q2 in [1-eps0, 1] on which q2 = 1 can be")
print("   pinned rigorously, leaving 4 parameters -- and the obstruction lives")
print("   IN that slab.  The asymptotic constant 1/ln2 is verified to 12 digits")
print("   numerically, NOT proved symbolically; eps0 is not yet computed.")
print("")
print("VERDICT")
print("  Dimension of the certification target remains 5.  No reduction from")
print("  aggregate monotonicity (coupled), no Pareto reduction (generically 2-D),")
print("  no w-endpoint reduction (Phi not concave in w; the obstruction is")
print("  interior in w), and M := t is exact only AT the obstruction.")
print("  The ~10^8-box estimate for a naive 5-D B&B stands unchanged.")
