"""gate_c_b5.py — Gate C(b): the degree-5 analogue of the affine coefficient inequality.

DERIVATION (from the paper's own chain; line numbers refer to scratch/paper_full.txt).

The paper's barrier (Prop 10.2, lines 1265-1291) uses exactly two inputs:
    0 < b1 <= 1        and       b3 >= 2 b1 - 11/6,
and concludes gamma <= 11/12, hence K_G >= 6 pi/11 (lines 1328-1339).
Its degree-3 input is the AFFINE functional (lines 966-981)
    A_0(f,g) := 2 b1 - b3 <= 11/6,
sharp at the hyperplane pair f = g = sgn(X_1), where (b1,b3) = (1, 1/6) (lines 969-973).

Degree-5 analogue, for lambda >= 0:
    A_lambda(f,g) := 2 b1 - b3 - lambda b5 <= C_lambda,
equivalently b3 + lambda b5 >= 2 b1 - C_lambda. A_lambda is affine in the pair's
coefficient vector, so like the paper's own inequality (line 1257) it survives mixtures
and coefficientwise limits — the property the Naor-Regev transfer (Section 11) needs.

BARRIER ARITHMETIC WITH b5 (derived here; numerically re-checked in __main__).
Reversion of H(t) = b1 t + b3 t^3 + b5 t^5 + ... gives (line 1276 for a1,a3; a5 by the
same coefficient comparison, cross-checked against Heilman 2606.00247 section 15.5):
    a1 = 1/b1,  a3 = -b3/b1^4,  a5 = (3 b3^2 - b1 b5)/b1^7.
An adversary keeping M_H(c) = sum |a_{2j+1}| c^{2j+1} <= 1 sets b3 = max(0, 2b1 - 11/6)
(smallest |a3|) and b5 = max(0, (2 b1 - C_lambda)/lambda). With b3 = b5 = 0,
M_H(c) >= c/b1 exceeds 1 only when b1 < c; so the largest b1 with both coefficients
annihilated is b1 = min(11/12, C_lambda/2) and the barrier becomes
    gamma <= min(11/12, C_lambda/2).
=> the b5 route TIGHTENS the barrier iff C_lambda < 11/6, and then
    gamma <= C_lambda / 2,     K_G >= pi / C_lambda.

DECISIVE QUANTITY. The hyperplane has b5 = 3/40 > 0, so A_lambda(hyperplane)
= 11/6 - (3/40) lambda < 11/6: lambda > 0 lowers the constant there. Since
A_lambda = A_0 - lambda b5, the constant is pushed back to 11/6 exactly by pairs with
b5 <= 0 and A_0 at (or approaching) 11/6. Hence Gate C(b) reduces to
    C_neg := sup { A_0(f,g) = 2 b1 - b3  :  b5(f,g) <= 0 }.
  * C_neg = 11/6  =>  C_lambda >= 11/6 for every lambda > 0: the b5 route is DEAD.
  * C_neg < 11/6  =>  any 0 < lambda <= (11/6 - C_neg)/|b5| at the maximiser gives
    C_lambda < 11/6 and a strictly better barrier.
Every family member below is a genuine odd sign pair, so each value is a valid LOWER
bound on C_neg (and on C_lambda).

Exact moment identities used:
    d/dx [He_{m-1}(x) phi(x)] = -He_m(x) phi(x)
  => int_a^b psi_m phi = [He_{m-1}(a) phi(a) - He_{m-1}(b) phi(b)] / sqrt(m!)
    He_0 = 1, He_1 = x, He_2 = x^2-1, He_3 = x^3-3x, He_4 = x^4-6x^2+3.
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flint import arb, ctx
from core import pi, nu as nu_const, phi_density, gauss_cdf

PREC = 300
ELEVEN_SIXTHS = arb(11) / arb(6)


def phi(x: arb) -> arb:
    return phi_density(x, PREC)


def He(m: int, x: arb) -> arb:
    if m == 0:
        return arb(1)
    if m == 1:
        return x
    h0, h1 = arb(1), x
    for n in range(1, m):
        h0, h1 = h1, x * h1 - arb(n) * h0
    return h1


def fact_sqrt(m: int) -> arb:
    f = 1
    for i in range(2, m + 1):
        f *= i
    return arb(f).sqrt()


def seg_psi(m: int, a: arb, b) -> arb:
    """int_a^b psi_m phi  (b = None means +infinity)."""
    ta = He(m - 1, a) * phi(a)
    tb = arb(0) if b is None else He(m - 1, b) * phi(b)
    return (ta - tb) / fact_sqrt(m)


def odd_coeffs(segs, ms=(1, 3, 5)):
    """<u, psi_m> for ODD u equal to sigma_i on {a_i <= |x| < a_{i+1}}."""
    out = {}
    for m in ms:
        tot = arb(0)
        for (a, b, sg) in segs:
            if sg == 0:
                continue
            tot += arb(sg) * seg_psi(m, a, b)
        out[m] = 2 * tot
    return out


def even_moment(m: int, c: arb) -> arb:
    """E[1_{|S|<c} psi_m(S)] for even m."""
    if m == 0:
        return 2 * gauss_cdf(c, PREC) - 1
    return 2 * (He(m - 1, arb(0)) * phi(arb(0)) - He(m - 1, c) * phi(c)) / fact_sqrt(m)


def pair_from_1d(segs_h, segs_k):
    """1-D pair f = h+k, g = h-k; b_m = (pi/2)(<h,psi_m>^2 - <k,psi_m>^2)."""
    H = odd_coeffs(segs_h)
    K = odd_coeffs(segs_k)
    pi2 = pi(PREC) / 2
    b = {m: pi2 * (H[m] ** 2 - K[m] ** 2) for m in (1, 3, 5)}
    return b[1], b[3], b[5], 2 * b[1] - b[3]


def pair_transverse(c: arb):
    """2-D Davie-Reeds-shaped pair: h = sgn(S)1_{|S|>=c}, k = sgn(T)1_{|S|<c}."""
    nu = nu_const(PREC)
    pi2 = pi(PREC) / 2
    H = odd_coeffs([(c, None, 1)])
    q = 2 * gauss_cdf(c, PREC) - 1
    e2 = even_moment(2, c)
    e4 = even_moment(4, c)
    h1, h3, h5 = H[1], H[3], H[5]
    k1sq = (q * nu) ** 2
    k3sq = (q * nu / arb(6).sqrt()) ** 2 + (e2 * nu) ** 2
    k5sq = (3 * q * nu / arb(120).sqrt()) ** 2 + (e4 * nu) ** 2 + (e2 * nu / arb(6).sqrt()) ** 2
    b1 = pi2 * (h1 ** 2 - k1sq)
    b3 = pi2 * (h3 ** 2 - k3sq)
    b5 = pi2 * (h5 ** 2 - k5sq)
    return b1, b3, b5, 2 * b1 - b3


if __name__ == "__main__":
    ctx.prec = PREC
    t0 = time.time()
    res = {}
    nu = nu_const(PREC)
    print(f"nu = {nu.str(20)}   11/6 = {ELEVEN_SIXTHS.str(20)}")

    print("\n=== hyperplane pair f = g = sgn(X1)  (paper lines 969-973) ===")
    b1, b3, b5, A0 = pair_from_1d([(arb(0), None, 1)], [])
    print(f"  b1 = {b1.str(18)}   (expect 1)")
    print(f"  b3 = {b3.str(18)}   (expect 1/6 = {(arb(1)/6).str(14)})")
    print(f"  b5 = {b5.str(18)}   (expect 3/40 = {(arb(3)/40).str(14)})")
    print(f"  A_0 = {A0.str(18)}   (expect 11/6 — the paper's SHARP constant)")
    res["hyperplane"] = dict(b1=b1.str(20), b3=b3.str(20), b5=b5.str(20), A0=A0.str(20))

    print("\n=== F1: 1-D strip pair (h = sgn 1_{|S|>=c}, k = sgn 1_{|S|<c}) ===")
    rows = []
    for cv in ("0.02", "0.05", "0.1", "0.2", "0.25573", "0.3", "0.4", "0.5", "0.7", "1.0", "1.5"):
        c = arb(cv)
        b1, b3, b5, A0 = pair_from_1d([(c, None, 1)], [(arb(0), c, 1)])
        gap = ELEVEN_SIXTHS - A0
        print(f"  c={cv:<8} b1 = {b1.str(12):>15} b3 = {b3.str(10):>14} b5 = {b5.str(10):>14} "
              f"A0 = {A0.str(12):>15} 11/6-A0 = {gap.str(8):>12} b5<=0: {bool(b5 <= 0)}")
        rows.append(dict(c=cv, b1=b1.str(18), b3=b3.str(18), b5=b5.str(18),
                         A0=A0.str(18), gap=gap.str(12), b5_nonpos=bool(b5 <= 0)))
    res["F1"] = rows

    print("\n=== F3: 2-D transverse pair (Davie-Reeds shape) ===")
    rows = []
    best_neg = None
    for cv in ("0.02", "0.05", "0.1", "0.15", "0.2", "0.25573", "0.3", "0.4", "0.5",
               "0.7", "1.0", "1.5", "2.0"):
        c = arb(cv)
        b1, b3, b5, A0 = pair_transverse(c)
        gap = ELEVEN_SIXTHS - A0
        neg = bool(b5 <= 0)
        print(f"  c={cv:<8} b1 = {b1.str(12):>15} b3 = {b3.str(10):>14} b5 = {b5.str(10):>14} "
              f"A0 = {A0.str(12):>15} 11/6-A0 = {gap.str(8):>12} b5<=0: {neg}")
        rows.append(dict(c=cv, b1=b1.str(18), b3=b3.str(18), b5=b5.str(18),
                         A0=A0.str(18), gap=gap.str(12), b5_nonpos=neg))
        if neg and (best_neg is None or A0 > best_neg[0]):
            best_neg = (A0, cv, b5)
    res["F3"] = rows
    if best_neg:
        A0n, cn, b5n = best_neg
        print(f"\n  best A_0 among b5 <= 0 members: A_0 = {A0n.str(18)} at c = {cn} "
              f"(b5 = {b5n.str(12)})")
        print(f"  11/6 - A_0 = {(ELEVEN_SIXTHS - A0n).str(12)}")
        res["C_neg_lower_bound"] = dict(A0=A0n.str(20), c=cn, b5=b5n.str(20),
                                        gap_to_11_6=(ELEVEN_SIXTHS - A0n).str(20))

    print(f"\n[{time.time()-t0:.1f}s]")
    with open("/Users/jinleic/jinleic-workspace/cs/kg/scratch/logs/gc_b5.json", "w") as f:
        json.dump(res, f, indent=1)
    print("wrote scratch/logs/gc_b5.json")


def lambda_ceiling(KG_upper_str="1.7818413238804372880"):
    """CERTIFIED structural bound on the b5 program.

    Any admissible gamma satisfies gamma <= C_lambda/2 (barrier, derived above), and the
    Naor-Regev transfer (paper Lemma 11.2, lines 1301-1327) supplies admissible gammas
    converging to pi/(2 K_G). Hence for EVERY lambda >= 0:
            C_lambda >= pi / K_G  >=  pi / K_G^upper.
    The hyperplane pair gives A_lambda = 11/6 - (3/40) lambda <= C_lambda. Therefore
            11/6 - (3/40) lambda  <=  C_lambda      is consistent only while
            11/6 - (3/40) lambda  >=  pi / K_G^upper   is NOT violated,
    i.e. a b5-strip whose constant is ATTAINED AT THE HYPERPLANE can exist only for
            lambda <= lambda* := (11/6 - pi/K_G^upper) * 40/3.
    For lambda > lambda*, sign pairs with A_lambda > 11/6 - (3/40) lambda MUST exist:
    the hyperplane is provably not the maximiser there. All arithmetic in Arb."""
    KG = arb(KG_upper_str)
    floor_C = pi(PREC) / KG                      # >= pi/K_G  (certified floor on C_lambda)
    lam_star = (ELEVEN_SIXTHS - floor_C) * arb(40) / arb(3)
    gamma_at_star = ELEVEN_SIXTHS / 2 - arb(3) * lam_star / arb(80)
    KG_from_star = pi(PREC) / (2 * gamma_at_star)
    return dict(KG_upper=KG, C_floor=floor_C, lam_star=lam_star,
                gamma_at_star=gamma_at_star, KG_from_star=KG_from_star)


def family_lambda_max(rows_list):
    """lambda_max over the swept families: min over members with b5 <= 0 of
    (11/6 - A_0)/|b5| — the largest lambda for which NO swept member exceeds 11/6."""
    best = None
    for rows in rows_list:
        for r in rows:
            if not r["b5_nonpos"]:
                continue
            b5 = arb(r["b5"].split(" ")[0].lstrip("["))
            A0 = arb(r["A0"].split(" ")[0].lstrip("["))
            if b5 >= 0:
                continue
            lam = (ELEVEN_SIXTHS - A0) / abs(b5)
            if best is None or lam < best[0]:
                best = (lam, r)
    return best
