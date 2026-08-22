"""High-precision check of the counterexample to Yu's decomposition step.

Yu (arXiv:2212.00658, Sec. 4) asserts, for P_pq written as a Caratheodory
combination of extreme points Q_i of P_B = {symmetric P_pq : E[p] <= t}:

    g(sum_i gamma_i Q_i, alpha)  >=  sum_i gamma_i g(Q_i, alpha),          (*)

justified by "Ppq -> g(Ppq, alpha) is concave".  This file exhibits explicit
extreme points violating (*), at mpmath precision, and separately records that
the RATIO conclusion Yu actually needs survives on this instance.

THE COUNTEREXAMPLE (all parameters exact rationals except t):

    Q1 = Q_{3/20, 3/20} = delta_{(3/20, 3/20)}          (beta = 0, a = 3/20 <= t)
    Q2 = (1-beta) Q_{1/10,1/10} + beta Q_{1,1}
       = (1-beta) delta_{(1/10,1/10)} + beta delta_{(1,1)}
                                       (a = 1/10 <= t < b = 1, beta = (t-1/10)/(9/10))
    gamma = 1/2

Both have exactly Yu's extreme-point shape, and P = Q1/2 + Q2/2 lies in P_B
because its marginal mean is (3/20 + t)/2 <= t.
"""

import os
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from entropy import PSI, h

mp.dps = 40

C_STAR = mpf("0.382345533366702721")
ALPHA = mpf("0.0356069")


def maxent_or(p, r):
    """s*(p,r) = median{1/2, max(p,r), min(p+r,1)} -- the max-entropy OR prob."""
    lo, hi = max(p, r), min(p + r, mpf(1))
    return min(max(mpf("0.5"), lo), hi)


def law_terms(atoms, weights, pairs, alpha):
    """(iid, coupled, Eh, g) for a joint law given as atom/weight marginal + pairs."""
    iid = sum(
        wi * wj * h(atoms[i] + atoms[j] - atoms[i] * atoms[j])
        for i, wi in enumerate(weights)
        for j, wj in enumerate(weights)
    )
    cpl = sum(wt * h(maxent_or(p, r)) for p, r, wt in pairs)
    eh = sum(w * h(a) for a, w in zip(atoms, weights))
    return iid, cpl, eh, (1 - alpha) * iid + alpha * cpl


def check(t, label):
    a = mpf(3) / 20
    c = mpf(1) / 10
    beta = (t - c) / (mpf(9) / 10)
    gamma = mpf(1) / 2
    assert 0 < beta < 1, beta

    # Q1: single atom a, joint delta_(a,a)
    q1 = law_terms([a], [mpf(1)], [(a, a, mpf(1))], ALPHA)
    # Q2: marginal (1-beta) d_c + beta d_1, joint (1-beta) d_(c,c) + beta d_(1,1)
    q2 = law_terms([c, mpf(1)], [1 - beta, beta],
                   [(c, c, 1 - beta), (mpf(1), mpf(1), beta)], ALPHA)
    # P = gamma Q1 + (1-gamma) Q2
    p_atoms = [a, c, mpf(1)]
    p_w = [gamma, (1 - gamma) * (1 - beta), (1 - gamma) * beta]
    p_pairs = [(a, a, gamma), (c, c, (1 - gamma) * (1 - beta)),
               (mpf(1), mpf(1), (1 - gamma) * beta)]
    pp = law_terms(p_atoms, p_w, p_pairs, ALPHA)

    mean_p = sum(w * x for x, w in zip(p_atoms, p_w))
    g_mix = pp[3]
    g_avg = gamma * q1[3] + (1 - gamma) * q2[3]
    defect = g_mix - g_avg

    print(f"\n=== {label}:  t = {mp.nstr(t, 18)} ===")
    print(f"  beta = {mp.nstr(beta, 15)}   marginal mean of P = "
          f"{mp.nstr(mean_p, 15)}  (must be <= t)")
    assert mean_p <= t + mpf(10) ** -30, (mean_p, t)

    print(f"  g(Q1)              = {mp.nstr(q1[3], 15)}")
    print(f"  g(Q2)              = {mp.nstr(q2[3], 15)}")
    print(f"  g(P)               = {mp.nstr(g_mix, 15)}")
    print(f"  gamma g(Q1)+..     = {mp.nstr(g_avg, 15)}")
    print(f"  DEFECT g(P) - avg  = {mp.nstr(defect, 12)}")
    assert defect < -mpf("0.001"), f"expected a clear violation, got {defect}"
    print("  => (*) FAILS: the concavity step in Yu Sec. 4 is invalid   [OK]")

    # The ratio conclusion Yu actually needs, on this same instance.
    r_p = g_mix / pp[2]
    r_1 = q1[3] / q1[2]
    r_2 = q2[3] / q2[2]
    print(f"  ratio(P)  = {mp.nstr(r_p, 12)}")
    print(f"  ratio(Q1) = {mp.nstr(r_1, 12)},  ratio(Q2) = {mp.nstr(r_2, 12)}")
    print(f"  min over extreme points = {mp.nstr(min(r_1, r_2), 12)}")
    if r_p >= min(r_1, r_2):
        print("  => the RATIO conclusion nevertheless HOLDS here:")
        print("     the counterexample refutes the PROOF, not the result.")
    else:
        print("  => the RATIO conclusion ALSO fails: the reduction is unsound.")
    return defect, r_p >= min(r_1, r_2)


print("Counterexample to the concavity step of Yu, arXiv:2212.00658 Sec. 4")
d1, ok1 = check(C_STAR, "at Cambie/Yu's optimum c*")
d2, ok2 = check(PSI, "at the rigorous barrier psi")

print("\n--- summary ---")
print(f"  defect at c*  : {mp.nstr(d1, 10)}")
print(f"  defect at psi : {mp.nstr(d2, 10)}")
print(f"  ratio conclusion survived on both instances: {ok1 and ok2}")
print("\nInterpretation: the inequality (*) that Yu's Krein-Milman reduction")
print("relies on is false on his own extreme points, so the published proof of")
print("Gamma(t) >= Gammahat(t) is invalid.  Independent unrestricted search")
print("(robust.py) finds no violation of the CONCLUSION, so c* is likely the")
print("true threshold -- but it currently has no valid proof.")
print("\nALL PASS")
