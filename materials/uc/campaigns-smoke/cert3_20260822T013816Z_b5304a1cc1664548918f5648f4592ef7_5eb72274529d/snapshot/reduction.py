"""The Gilmer iid term, in closed form:  Q = 2(1-mean) L - (PSD defect).

This is Theorem A of decomposition.py, rearranged into the form that actually
does work.  Everything is exact; the series is used only to prove positivity.

-------------------------------------------------------------------------------
NOTATION.  mu is the law of the inclusion probability p on [0,1]; x = 1-p.

    Q(mu) = int int h(p (+) q) dmu dmu       (Gilmer's iid / OR term)
    L(mu) = int h(p) dmu                     (the entropy being beaten)
    B(mu) = int (1-p) dmu = 1 - mean(mu)
    T(u)  = u + (1-u) ln(1-u) = sum_{n>=2} u^n/(n(n-1))

-------------------------------------------------------------------------------
THEOREM A' (exact reduction).  Define

    W(x,y) = xy - x T(y) - y T(x) + T(xy).

Then W is a POSITIVE-DEFINITE kernel, because

    W(x,y) = sum_{n>=2} (x - x^n)(y - y^n) / (n(n-1)),                     (1)

and for every law mu, writing E(mu) = (1/ln2) int int W dmu dmu >= 0,

    Q(mu) = 2 B(mu) L(mu) - E(mu).                                         (2)

In moment form, with M_n = int x^n dmu (so M_1 = B),

    ln2 * E(mu) = sum_{n>=2} (M_1 - M_n)^2 / (n(n-1)),                     (3)

every term of which is a squared *defect* M_1 - M_n >= 0.

Proof.  Theorem A says ln2*Q = (A+B)^2 - A^2 - int int T(xy), A = int(-x ln x).
Also ln2*L = int H(x) dmu = A + B - int T(x) dmu since H(u) = -u ln u + u - T(u).
Expanding W and using int int x T(y) = B int T(x) dmu,

    int int W dmu dmu = B^2 - 2 B int T(x) dmu + int int T(xy),

so  2 B (ln2 L) - int int W = 2AB + 2B^2 - 2B intT - B^2 + 2B intT - int int T
                            = 2AB + B^2 - int int T(xy)  =  ln2 Q.          []

DIAGONAL.  Putting y = x in (2) with mu = delta_p (so x = 1-p, B = x, L = h(p)
= h(x), Q = h(x^2)):

    W(x,x) / ln2 = 2 x h(x) - h(x^2).                                      (4)

-------------------------------------------------------------------------------
COROLLARY 1 (the whole two-strategy functional).  With C the coupled term,

    F(mu) = (1-alpha) Q + alpha C - L
          = [ 2(1-alpha) B - 1 ] L  +  alpha C  -  (1-alpha) E.            (5)

Every term is explicit, and all four of L, C, E, and F vanish simultaneously on
the sink family (1-u) delta_0 + u delta_1.

COROLLARY 2 (the identity is SHARP: it reproduces the record constant).  Take
alpha = 0 and mu = delta_p.  Then (2) and (4) give

    Q - L = (2x-1) h(x) - [2x h(x) - h(x^2)] = h(x^2) - h(x),   x = 1-p,

so F(delta_p) >= 0  <=>  h((1-p)^2) >= h(1-p)  <=>  (1-p)^2 >= p  <=>  p <= psi,

with equality exactly at the golden point x = 1/phi, where x^2 = 1-x.  So the
barrier psi = (3-sqrt5)/2 is *precisely* the statement that the PSD defect E
first exceeds (2B-1)L, and the golden ratio enters only through x^2 = 1-x.

THE alpha = 0 HYPOTHESIS IS ESSENTIAL.  For alpha > 0 the coupled term does NOT
drop out of F(delta_p): s*(p,p) = 1/2 for p in (1/4,1/2), so h(s*) = 1 near psi
and

    F(delta_p) = (1-alpha) h((1-p)^2) + alpha h(s*(p,p)) - h(p),
    F(delta_psi) = alpha [1 - h(psi)] = 0.0405813 alpha  >  0.

So point masses do not cap at psi once alpha > 0.  But the point-mass threshold
does NOT exceed c* for every alpha > 0: at p = c* the pure-iid part is already
negative, h((1-c*)^2) - h(c*) = -5.8929e-4 < 0 because c* > psi, and
F_point(c*, .) is AFFINE in alpha, so it stays negative for all small alpha.  The
exact crossover is

    alpha_crit = [h(c*) - h((1-c*)^2)] / [1 - h((1-c*)^2)] = 0.0144054585140081,

and
    point-mass threshold > c*   <=>   alpha > alpha_crit.

Cambie's alpha = 0.0356069 exceeds alpha_crit, so at THAT alpha point masses are
not the binding case and the 2-atom obstruction {1,b} is; for alpha < alpha_crit
point masses bind strictly below c*.  Corollary 2 itself is a statement about the
PURE IID problem (alpha = 0), and that is where it is sharp.
Verified numerically in section 6 below, with alpha sampled on BOTH sides of
alpha_crit.

COROLLARY 3 (equivalent moment form of the psi theorem).  For alpha = 0 the
one-step bound at t is exactly the assertion

    sum_{n>=2} (M_1 - M_n)^2 / (n(n-1))  <=  (2 M_1 - 1) * int H(x) dmu     (6)

for every law with M_1 >= 1 - t.  This is a pure moment inequality: no entropy
appears on the left, and the right side is linear in mu.

-------------------------------------------------------------------------------
WHAT THE IDENTITY DOES NOT GIVE.  Bounding E by dropping the variance,

    E <= (1/ln2) int W(x,x) dmu = int [ 2(1-p) h(p) - h(2p - p^2) ] dmu     (7)

(Jensen, tight exactly at point masses) yields only the constant

    min_p (2p + lambda(p) - 1)/2 = 0.3475266...  <  psi,

because the maximiser of the pointwise ratio sits at x = 0.8835, i.e. p = 0.117,
far below the mean bound.  The variance term is therefore load-bearing, which is
why beating psi needs the support reduction (Theorem B', pair-orbit form) and not
a pointwise inequality.  That is checked below rather than asserted.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import PSI, h, mp, mpf, sawin_lambda
from mpmath import diff, findroot, log, sqrt

LN2 = log(2)


def T(u):
    """T(u) = u + (1-u) ln(1-u) = sum_{n>=2} u^n/(n(n-1))."""
    u = mpf(u)
    return mpf(1) if u >= 1 else u + (1 - u) * log(1 - u)


def W(x, y):
    x, y = mpf(x), mpf(y)
    return x * y - x * T(y) - y * T(x) + T(x * y)


def parts(atoms, weights):
    at = [mpf(a) for a in atoms]
    wv = [mpf(v) for v in weights]
    n = len(at)
    Q = sum(wv[i] * wv[j] * h(at[i] + at[j] - at[i] * at[j])
            for i in range(n) for j in range(n))
    L = sum(wv[i] * h(at[i]) for i in range(n))
    B = sum(wv[i] * (1 - at[i]) for i in range(n))
    E = sum(wv[i] * wv[j] * W(1 - at[i], 1 - at[j])
            for i in range(n) for j in range(n)) / LN2
    return Q, L, B, E


# ---------------------------------------------------------------------------
print("1. W is positive definite: W(x,y) = sum_{n>=2} (x-x^n)(y-y^n)/(n(n-1))")
print("   The rearrangement is termwise over four convergent series:")
print("     sum (xy)/(n(n-1)) = xy   (telescoping 1/(n-1) - 1/n)")
print("     sum x y^n/(n(n-1)) = x T(y),  sum y x^n/... = y T(x),")
print("     sum (xy)^n/(n(n-1)) = T(xy).")
tele = mpf(0)
NT = 20000
for nn in range(2, NT + 1):
    tele += mpf(1) / (nn * (nn - 1))
assert abs(tele - (1 - mpf(1) / NT)) < mpf(10) ** -30
print(f"   telescoping check: partial_{NT} = 1 - 1/{NT} exactly           [OK]")

NS = 600
print(f"\n   numeric check with a RIGOROUS tail bracket (N = {NS}):")
print(f"   {'x':>7} {'y':>7} {'W closed':>20} {'partial + xy/N':>20} {'|diff|':>11}")
for xv, yv in (("0.3", "0.7"), ("0.9", "0.2"), ("0.5", "0.5"), ("0.85", "0.9")):
    x, y = mpf(xv), mpf(yv)
    part = sum((x - x**nn) * (y - y**nn) / (nn * (nn - 1))
               for nn in range(2, NS + 1))
    approx = part + x * y / NS
    # the three remaining tails are each bounded by max(x,y)^(N+1)
    bound = 4 * max(x, y) ** (NS + 1)
    d = abs(W(x, y) - approx)
    print(f"   {xv:>7} {yv:>7} {mp.nstr(W(x, y), 16):>20} "
          f"{mp.nstr(approx, 16):>20} {mp.nstr(d, 3):>11}")
    assert d <= bound + mpf(10) ** -30, (xv, yv, d, bound)
print("   agrees within the analytic tail bound                          [OK]")
assert W(1, 1) == 0 and W(0, 0) == 0
print("   W(0,0) = W(1,1) = 0 exactly (the sink endpoints)                [OK]")

n = 160
xg = np.linspace(1e-9, 1 - 1e-9, n)
Wm = np.array([[float(W(mpf(float(a)), mpf(float(b)))) for b in xg] for a in xg])
ev = np.linalg.eigvalsh(Wm)
print(f"   Gram spectrum on a {n}-point grid: min eig = {ev[0]:+.3e}")
assert ev[0] > -1e-11
print("   numerically PSD as well                                        [OK]")

# ---------------------------------------------------------------------------
print("\n2. THEOREM A':  Q = 2 B L - E")
print(f"   {'law':>44} {'Q':>19} {'2BL - E':>19} {'err':>10}")
laws = [
    ([0.2, 0.6], [0.5, 0.5]),
    ([0.0, 1.0], [0.618, 0.382]),
    ([0.1, 0.35, 0.9], [0.2, 0.5, 0.3]),
    ([1.0, 0.32945473850303697239], [0.078877292705923173, 0.921122707294076827]),
    ([0.4, 0.95], [0.5, 0.5]),
    ([0.05, 0.5, 0.5, 0.99], [0.4, 0.2, 0.1, 0.3]),
]
for atoms, wts in laws:
    Q, L, B, E = parts(atoms, wts)
    rhs = 2 * B * L - E
    print(f"   {str(atoms)[:42]:>44} {mp.nstr(Q, 13):>19} {mp.nstr(rhs, 13):>19} "
          f"{mp.nstr(abs(Q - rhs), 3):>10}")
    assert abs(Q - rhs) < mpf(10) ** -30, (atoms, Q - rhs)

rng = np.random.default_rng(2)
worst = mpf(0)
for _ in range(300):
    k = int(rng.integers(1, 6))
    at = [mpf(float(v)) for v in rng.uniform(0, 1, k)]
    wv = [mpf(float(v)) for v in rng.dirichlet(np.ones(k))]
    Q, L, B, E = parts(at, wv)
    worst = max(worst, abs(Q - (2 * B * L - E)))
print(f"   worst error over 300 random laws: {mp.nstr(worst, 4)}          [OK]")
assert worst < mpf(10) ** -28

# ---------------------------------------------------------------------------
print("\n3. moment form:  ln2 E = sum (M_1 - M_n)^2/(n(n-1)), defects >= 0")
print("   d_n = M_1 - M_n increases in n to d_inf = M_1 - mu({p=0}), and")
print("   sum_{n>N} 1/(n(n-1)) = 1/N exactly, so the tail is bracketed by")
print("   d_N^2/N  <=  tail  <=  d_inf^2/N.  That is a rigorous enclosure.")
NM = 4000
print(f"\n   {'law':>38} {'lower':>17} {'ln2*E':>17} {'upper':>17}")
for atoms, wts in laws[:5]:
    at = [mpf(a) for a in atoms]
    wv = [mpf(v) for v in wts]
    xs = [1 - a for a in at]
    M1 = sum(wv[i] * xs[i] for i in range(len(xs)))

    def defect(nn):
        return M1 - sum(wv[i] * xs[i] ** nn for i in range(len(xs)))

    part = sum(defect(nn) ** 2 / (nn * (nn - 1)) for nn in range(2, NM + 1))
    d_inf = M1 - sum(wv[i] for i in range(len(xs)) if xs[i] >= 1)
    lo = part + defect(NM) ** 2 / NM
    hi = part + d_inf**2 / NM
    _, _, _, E = parts(atoms, wts)
    print(f"   {str(atoms)[:36]:>38} {mp.nstr(lo, 14):>17} "
          f"{mp.nstr(E * LN2, 14):>17} {mp.nstr(hi, 14):>17}")
    assert lo - mpf(10) ** -30 <= E * LN2 <= hi + mpf(10) ** -30, (atoms, lo, hi)
    for nn in (2, 3, 10, 100):
        assert defect(nn) >= -mpf(10) ** -30, (atoms, nn)
    for nn in (2, 5, 50):
        assert defect(nn) <= defect(nn + 1) + mpf(10) ** -30, (atoms, nn)
print("   ln2*E lies inside the enclosure for every law tested            [OK]")
print("   every defect M_1 - M_n >= 0 (x^n <= x on [0,1])                [OK]")

# ---------------------------------------------------------------------------
print("\n4. diagonal identity:  W(x,x)/ln2 = 2 x h(x) - h(x^2)")
for xv in ("0.1", "0.3", "0.6180339887498948482", "0.9"):
    x = mpf(xv)
    lhs, rhs = W(x, x) / LN2, 2 * x * h(x) - h(x * x)
    assert abs(lhs - rhs) < mpf(10) ** -30
print("   verified to 1e-30 at four points, incl. the golden point       [OK]")

# ---------------------------------------------------------------------------
print("\n5. COROLLARY 1: F = [2(1-a)B - 1] L + a C - (1-a) E   (a = alpha)")
ALPHA = mpf("0.0356069")


def coupled_two_atom(atoms, wts):
    """Exact min-coupling for a <=2-atom law, in closed form."""
    at = [mpf(a) for a in atoms]
    wv = [mpf(v) for v in wts]
    if len(at) == 1:
        return h(smax(at[0], at[0]))
    c11, c12, c22 = (h(smax(at[0], at[0])), h(smax(at[0], at[1])),
                     h(smax(at[1], at[1])))
    base = wv[0] * c11 + wv[1] * c22
    gain = 2 * c12 - c11 - c22
    return base + min(wv[0], wv[1]) * min(mpf(0), gain)


def smax(p, r):
    p, r = mpf(p), mpf(r)
    return min(max(mpf("0.5"), max(p, r)), min(p + r, mpf(1)))


for atoms, wts in [([1.0, 0.32945473850303697239],
                    [0.078877292705923173, 0.921122707294076827]),
                   ([0.3819660112501051518], [1.0]),
                   ([0.2, 0.6], [0.5, 0.5])]:
    Q, L, B, E = parts(atoms, wts)
    C = coupled_two_atom(atoms, wts)
    direct = (1 - ALPHA) * Q + ALPHA * C - L
    viaid = (2 * (1 - ALPHA) * B - 1) * L + ALPHA * C - (1 - ALPHA) * E
    print(f"   {str(atoms)[:40]:>42}: F = {mp.nstr(direct, 13):>17}  "
          f"err {mp.nstr(abs(direct - viaid), 3)}")
    assert abs(direct - viaid) < mpf(10) ** -28
print("   Corollary 1 holds exactly                                      [OK]")

# ---------------------------------------------------------------------------
print("\n6. COROLLARY 2: at alpha = 0 the identity reproduces psi EXACTLY")
print(f"   {'p':>22} {'F(delta_p)':>18} {'h(x^2)-h(x)':>18} {'sign':>6}")
for pv in ("0.2", "0.35", "0.3819660112501051518", "0.39", "0.45"):
    p = mpf(pv)
    x = 1 - p
    Q, L, B, E = parts([p], [1])
    F0 = Q - L
    alt = h(x * x) - h(x)
    sign = "+" if F0 > mpf(10) ** -30 else ("0" if abs(F0) < mpf(10) ** -30 else "-")
    print(f"   {pv:>22} {mp.nstr(F0, 12):>18} {mp.nstr(alt, 12):>18} {sign:>6}")
    assert abs(F0 - alt) < mpf(10) ** -28
x_gold = (sqrt(5) - 1) / 2
assert abs(x_gold**2 - (1 - x_gold)) < mpf(10) ** -35
assert abs(1 - x_gold - PSI) < mpf(10) ** -35
print(f"   golden point x = 1/phi = {mp.nstr(x_gold, 20)}")
print(f"   x^2 = 1-x exactly, and 1-x = psi                              [OK]")
print("   => psi is exactly where the PSD defect E overtakes (2B-1)L")

print("\n   BUT the alpha = 0 hypothesis is ESSENTIAL.  For alpha > 0 the")
print("   coupled term survives at point masses: s*(p,p) = 1/2 near psi.")


def sstar_diag(p):
    p = mpf(p)
    return min(max(mpf("0.5"), p), min(2 * p, mpf(1)))


def F_point(p, al):
    p, al = mpf(p), mpf(al)
    return (1 - al) * h(2 * p - p * p) + al * h(sstar_diag(p)) - h(p)


print(f"   s*(psi,psi) = {mp.nstr(sstar_diag(PSI), 6)}, "
      f"h(s*) = {mp.nstr(h(sstar_diag(PSI)), 6)}, h(psi) = "
      f"{mp.nstr(h(PSI), 16)}")
C_STAR = mpf("0.38234553336670272115")

# F_point(p, .) is AFFINE in alpha, so the crossover at any fixed p is exact.
h_iid = h(2 * C_STAR - C_STAR * C_STAR)     # h((1-c*)^2)
h_c = h(C_STAR)
ALPHA_CRIT = (h_c - h_iid) / (1 - h_iid)
print(f"\n   h((1-c*)^2) = {mp.nstr(h_iid, 18)}")
print(f"   h(c*)       = {mp.nstr(h_c, 18)}")
print(f"   F_point(c*, 0) = {mp.nstr(h_iid - h_c, 10)}  < 0, since c* > psi")
print("   F_point(c*, .) is AFFINE in alpha, so it stays negative for all")
print("   sufficiently small alpha > 0.  The exact crossover is")
print("     alpha_crit = [h(c*) - h((1-c*)^2)] / [1 - h((1-c*)^2)]")
print(f"                = {mp.nstr(ALPHA_CRIT, 18)}")
assert abs(F_point(C_STAR, ALPHA_CRIT)) < mpf(10) ** -25
print("   verified: F_point(c*, alpha_crit) = 0                          [OK]")

# --- the missing step: UNIQUENESS of the root, so that "the threshold" exists
# and its position is decided by the sign at c*. --------------------------------
P0 = 1 - 1 / sqrt(2)          # 0.29289...: the p where 2p - p^2 = 1/2
assert abs((2 * P0 - P0 * P0) - mpf("0.5")) < mpf(10) ** -30
print(f"\n   p_0 := 1 - 1/sqrt2 = {mp.nstr(P0, 18)}, where 2p - p^2 = 1/2")
print(f"   1/4 < p_0 < c* < 1/2:  {mpf('0.25') < P0 < C_STAR < mpf('0.5')}")


def hprime(u):
    """h'(u) = log2((1-u)/u): positive for u < 1/2, negative for u > 1/2."""
    u = mpf(u)
    return log((1 - u) / u) / LN2


def dF_dp(p, al):
    """On (1/4,1/2), s*(p,p) = 1/2 is CONSTANT, so the alpha term drops out."""
    p, al = mpf(p), mpf(al)
    return (1 - al) * (2 - 2 * p) * hprime(2 * p - p * p) - hprime(p)


print("\n   LEMMA (uniqueness).  On (1/4,1/2), s*(p,p) = 1/2 is constant, so")
print("     dF/dp = (1-a)(2-2p) h'(2p-p^2) - h'(p).")
print("   For p in [p_0, 1/2) and 0 <= a < 1 every factor's sign is elementary:")
print("     2p-p^2 >= 1/2 <=> (1-p)^2 <= 1/2 <=> p >= p_0, so h'(2p-p^2) <= 0;")
print("     2-2p > 0; and h'(p) > 0 since p < 1/2.")
print("   Hence dF/dp < 0 STRICTLY on [p_0,1/2) -- proved, not sampled.")
print("   Also F >= 0 on [0,p_0] for EVERY a in [0,1]: there 2p-p^2 <= 1/2 and")
print("   2p-p^2 >= p give h(2p-p^2) >= h(p); and s*(p,p) is 2p (p<=1/4) or 1/2,")
print("   both in [p,1/2], so h(s*) >= h(p); hence F >= (1-a)h(p)+a h(p)-h(p) = 0.")
print("   And F(1/2, a) = -(1-a)(1 - h(3/4)) < 0 for a < 1.  So for each")
print("   a in [0,1) there is EXACTLY ONE root in (p_0,1/2), F >= 0 to its left.")

print(f"\n   {'alpha':>10} {'min dF/dp':>18} {'F(p_0)':>14} {'F(1/2)':>14}")
for al in ("0.0", "0.001", "0.0144", "0.0356069", "0.5", "0.99"):
    a = mpf(al)
    grid = [P0 + (mpf("0.4999") - P0) * mpf(i) / 400 for i in range(401)]
    mind = min(dF_dp(pp, a) for pp in grid)
    fp0, fhalf = F_point(P0, a), -(1 - a) * (1 - h(mpf(3) / 4))
    print(f"   {al:>10} {mp.nstr(mind, 12):>18} {mp.nstr(fp0, 10):>14} "
          f"{mp.nstr(fhalf, 10):>14}")
    assert mind < 0 and fp0 > 0 and fhalf < 0, (al, mind, fp0, fhalf)
    assert abs(F_point(mpf("0.5"), a) - fhalf) < mpf(10) ** -25, al
print("   corroborates the proof at 401 points per alpha                 [OK]")


def threshold(al):
    """The unique root in (p_0,1/2) by bisection -- licensed by the lemma."""
    a = mpf(al)
    lo, hi = P0, mpf("0.5")
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if F_point(mid, a) > 0 else (lo, mid)
    return (lo + hi) / 2


print(f"\n   {'alpha':>11} {'F(delta_psi)':>18} {'F(delta_c*)':>18} "
      f"{'threshold':>20} {'vs c*':>7}")
for al in ("0.0", "0.001", "0.005", "0.01", "0.0144", "0.0145", "0.02",
           "0.0356069", "0.1", "0.3", "0.99"):
    a = mpf(al)
    fpsi, fcs = F_point(PSI, a), F_point(C_STAR, a)
    thr = threshold(a)
    tag = "ABOVE" if thr > C_STAR else "BELOW"
    print(f"   {al:>11} {mp.nstr(fpsi, 12):>18} {mp.nstr(fcs, 12):>18} "
          f"{mp.nstr(thr, 18):>20} {tag:>7}")
    if a > 0:
        assert fpsi > 0, (al, fpsi)
    # position decided by the SIGN AT c* (licensed by the lemma), and
    # F_point(c*, .) is affine in alpha with root alpha_crit
    assert (thr > C_STAR) == (fcs > 0), (al, thr, fcs)
    assert (fcs > 0) == (a > ALPHA_CRIT), (al, fcs)
assert abs(threshold(mpf(0)) - PSI) < mpf(10) ** -25
print("   threshold(0) = psi to 1e-25, as Corollary 2 requires           [OK]")

print("\n   alpha = 1 SEPARATELY: F(p,1) = h(s*(p,p)) - h(p) >= 0 for all p,")
print("   with equality on [1/2,1] where s*(p,p) = p.  The strict decrease")
print("   fails there, and the threshold is 1/2 > c* -- still consistent")
print("   with the equivalence, since 1 > alpha_crit.")
for pv in ("0.1", "0.3", "0.45", "0.5", "0.7", "0.9"):
    p = mpf(pv)
    f1 = h(sstar_diag(p)) - h(p)
    assert f1 >= -mpf(10) ** -30, (pv, f1)
    tag = "   (= 0, p >= 1/2)" if p >= mpf("0.5") else ""
    print(f"     F({pv:>4}, 1) = {mp.nstr(f1, 10):>16}{tag}")
print("   alpha = 1 handled                                              [OK]")
assert abs(F_point(PSI, mpf("0.0356069"))
           - mpf("0.0356069") * (1 - h(PSI))) < mpf(10) ** -25
print(f"\n   F(delta_psi) = alpha*(1 - h(psi)) = {mp.nstr(1 - h(PSI), 8)} * alpha")
print(f"   c* = {mp.nstr(C_STAR, 18)}")
print("   REGRESSION (an earlier version of this file sampled only")
print("   alpha in {0, 0.0356069, 0.1, 0.3} and then claimed 'for EVERY")
print("   alpha > 0', skipping the whole interval (0, alpha_crit) where the")
print("   claim is FALSE).  Correct statement:")
print(f"     point-mass threshold > c*   <=>   alpha > {mp.nstr(ALPHA_CRIT, 12)}")
print(f"   Cambie's alpha = 0.0356069 > alpha_crit, so at THAT alpha point")
print("   masses are not binding and the 2-atom obstruction is.  For")
print("   alpha < alpha_crit point masses bind BELOW c*.               [OK]")
print("   Corollary 2 itself is about the PURE IID problem (alpha = 0),")
print("   and is sharp there.")

# ---------------------------------------------------------------------------
print("\n7. COROLLARY 3: the psi bound as a pure moment inequality")
print("   sum (M_1-M_n)^2/(n(n-1))  <=  (2 M_1 - 1) int H dmu   for M_1 >= 1/phi")


def moment_ineq(atoms, wts):
    Q, L, B, E = parts(atoms, wts)
    return (2 * B - 1) * L * LN2 - E * LN2


print(f"   {'law':>40} {'mean':>10} {'slack':>16} {'ok?':>5}")
for atoms, wts in [([0.3], [1.0]), ([0.3819660112501051518], [1.0]),
                   ([0.39], [1.0]), ([0.1, 0.5], [0.5, 0.5]),
                   ([0.0, 1.0], [0.618, 0.382]),
                   ([1.0, 0.32945473850303697239],
                    [0.078877292705923173, 0.921122707294076827])]:
    m = sum(mpf(wts[i]) * mpf(atoms[i]) for i in range(len(atoms)))
    s = moment_ineq(atoms, wts)
    ok = "yes" if s >= -mpf(10) ** -25 else "NO"
    print(f"   {str(atoms)[:38]:>40} {mp.nstr(m, 8):>10} {mp.nstr(s, 10):>16} "
          f"{ok:>5}")

viol = 0
for _ in range(4000):
    k = int(rng.integers(1, 4))
    at = rng.uniform(0, 1, k)
    wv = rng.dirichlet(np.ones(k))
    if float(wv @ at) > float(PSI):
        continue
    s = moment_ineq([mpf(float(v)) for v in at], [mpf(float(v)) for v in wv])
    if s < -mpf(10) ** -22:
        viol += 1
print(f"   random search over laws with mean <= psi: {viol} violations of (6)")
assert viol == 0
print("   consistent with the known psi theorem                          [OK]")

# ---------------------------------------------------------------------------
print("\n8. WHAT THE POINTWISE BOUND GIVES (and why it is not enough)")
print("   dropping the variance: E <= int [2(1-p)h(p) - h(2p-p^2)] dmu")


def ratio_pointwise(x):
    x = mpf(x)
    return 2 * x - h(x * x) / h(x)


best = (mpf(-10), None)
N = 4000
for i in range(1, N):
    xx = mpf(i) / N
    v = ratio_pointwise(xx)
    if v > best[0]:
        best = (v, xx)
x0, step = best[1], mpf(1) / N
for _ in range(200):
    for s in (-step, step):
        if 0 < x0 + s < 1 and ratio_pointwise(x0 + s) > best[0]:
            best = (ratio_pointwise(x0 + s), x0 + s)
            x0 = x0 + s
    step /= 2
kap, xstar = best
print(f"   kappa* = max_x [2x - h(x^2)/h(x)] = {mp.nstr(kap, 18)}")
print(f"   attained at x* = {mp.nstr(xstar, 18)}  (p* = {mp.nstr(1 - xstar, 12)})")
print(f"   implied constant (1-kappa*)/2  = {mp.nstr((1 - kap) / 2, 18)}")
print(f"   psi                            = {mp.nstr(PSI, 18)}")
print(f"   shortfall = {mp.nstr((1 - kap) / 2 - PSI, 6)}  -- WEAKER than psi")
assert (1 - kap) / 2 < PSI
print("   the maximiser sits at p = 0.117, far below the mean bound, so")
print("   the pointwise bound wastes the constraint.  The variance term is")
print("   load-bearing => a support reduction is required, not an inequality.")

print("\n   for comparison, at the golden point the pointwise ratio is exact:")
print(f"     2x - h(x^2)/h(x) at x = 1/phi : "
      f"{mp.nstr(ratio_pointwise(x_gold), 18)}")
print(f"     sqrt(5) - 2                   : {mp.nstr(sqrt(5) - 2, 18)}")
assert abs(ratio_pointwise(x_gold) - (sqrt(5) - 2)) < mpf(10) ** -28
print("     equal, and (1-(sqrt5-2))/2 = psi                             [OK]")
print(f"   d/dx[2x - h(x^2)/h(x)] at 1/phi = "
      f"{mp.nstr(diff(ratio_pointwise, x_gold), 12)}  (= psi, so 1/phi is")
print("     NOT the maximiser -- the pointwise route leaks there)")
assert abs(sawin_lambda(PSI) - 1) < mpf(10) ** -20

print("\nALL PASS")
print("\nSUMMARY")
print("  Q = 2(1-mean) L - E,   E = (1/ln2) int int W dmu dmu >= 0")
print("  F = [2(1-alpha)(1-mean) - 1] L + alpha C - (1-alpha) E")
print("  sharp on point masses: F >= 0 <=> p <= psi, via x^2 = 1-x")
print("  pointwise relaxation gives only 0.34753, so the variance matters")
