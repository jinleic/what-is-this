"""The OR-entropy kernel has AT MOST ONE POSITIVE SQUARE, and what follows.

Theorem A is elementary and exact.  Theorem B' is the support reduction; an
earlier THEOREM B claimed a stronger conclusion from a SIGN ERROR and is
retracted below.

-------------------------------------------------------------------------------
NOTATION.  Natural-log entropy H(u) = -u ln u - (1-u) ln(1-u), so H = ln2 * h.
For a law mu of the inclusion probability p, substitute x = 1-p and push mu
forward to nu.  Since p (+) q := p+q-pq = 1-(1-p)(1-q) and H(1-z) = H(z),

    H(p (+) q) = H(xy),        x = 1-p,  y = 1-q,

so the Gilmer/Sawin iid term is the quadratic form of the kernel H(xy).

-------------------------------------------------------------------------------
LEMMA 1 (elementary).  Put  T(u) := u + (1-u) ln(1-u),  T(0) = 0, T(1) = 1.
Then H(u) = -u ln u + u - T(u), and

    T(u) = sum_{n>=2} u^n / (n(n-1))                                       (*)

with every coefficient positive, so T(xy) is a POSITIVE-DEFINITE kernel:
int int T(xy) dnu dnu = sum_{n>=2} M_n(nu)^2/(n(n-1)) >= 0,  M_n = int x^n dnu.

Proof of (*).  -(1-u)ln(1-u) = (1-u) sum_{n>=1} u^n/n
              = sum_{n>=1} u^n/n - sum_{n>=2} u^n/(n-1)
              = u - sum_{n>=2} u^n/(n(n-1)).                                []

-------------------------------------------------------------------------------
THEOREM A (one positive square).  For every finite signed measure nu on [0,1],

    int int H(xy) dnu dnu  =  G(nu)^2 - A(nu)^2 - int int T(xy) dnu dnu,

    A = int (-x ln x) dnu,   B = int x dnu,   G = A + B = int x(1 - ln x) dnu.

Since the last term is positive definite, the quadratic form of H(xy) is
(one positive square) minus (a PSD form): H(xy) has AT MOST ONE POSITIVE SQUARE,
and the unique positive direction is the explicit linear functional G.

Proof.  By Lemma 1 with u = xy,
    H(xy) = a(x)b(y) + b(x)a(y) + b(x)b(y) - T(xy),  a = -x ln x,  b = x,
and integrating gives 2AB + B^2 - int int T = (A+B)^2 - A^2 - int int T.    []

-------------------------------------------------------------------------------
RETRACTED THEOREM B.  An earlier version of this file claimed the infimum of

    F(mu) = (1-alpha) Q(mu) + alpha C(mu) - L(mu),
    Q = int int h(p (+) q) dmu dmu,   L = int h dmu,
    C(mu) = min over symmetric couplings M of mu with itself of int h(s*) dM,

over { E_mu[p] <= t } is attained on <= 3 atoms, via Bauer's minimum principle
after freezing G.  That needs Psi = (1-alpha)Psi_Q + alpha C - L to be CONCAVE.
It is not.  The file's own argument -- if M_i couples mu_i then a convex
combination couples the combination, so

    C(theta mu_1 + (1-theta) mu_2)  <=  theta C(mu_1) + (1-theta) C(mu_2)

-- is the definition of CONVEXITY, and I read it as concavity.  Monge-Kantorovich
duality confirms it independently: C(mu) = sup { 2 int phi dmu : phi obeys (A) }
is a supremum of linear functionals, hence convex.  A witness is computed below:
C(d_0.3) = 1, C(d_0.6) = 0.97095, but C of the half-half mixture is 0.97095,
below the chord 0.98548 by 1.45e-2.  Concavity would need it ABOVE.

So Bauer does not apply, and every consequence is withdrawn: the "<= 3 marginal
atoms" claim, the "same dimension as Yu's 5-parameter family" claim, and the
reading of the k = 2..6 searches as confirmation of a theorem.  Those searches
remain valid as evidence; they were never a proof.

-------------------------------------------------------------------------------
THEOREM B' (pair-orbit support reduction -- the correct version).  Let

    Delta = { (p,q) : 0 <= p <= q <= 1 },      c(p,q) = h(s*(p,q)).

For nu in P(Delta) let mu_nu = int (delta_p + delta_q)/2 dnu be the INDUCED
MARGINAL, and put

    Phi(nu) = (1-alpha) Q(mu_nu) + alpha int c dnu - L(mu_nu).

Then
    inf { F(mu) : mu in P[0,1], E_mu[p] <= t }
  = inf { Phi(nu) : nu in P(Delta), E_{mu_nu}[p] <= t },                    (1)

and the right-hand infimum is attained at some nu with AT MOST 3 ATOMS, i.e. a
symmetric coupling carried by at most 3 unordered pairs -- hence an induced
marginal with at most 6 atoms.

Proof.  (i) Symmetric probability measures M on [0,1]^2 correspond to nu in
P(Delta), with the marginal of M equal to mu_nu.  For any nu, int c dnu >= C(mu_nu)
since M is *some* symmetric coupling of its own marginal, so Phi(nu) >= F(mu_nu);
conversely for each mu an optimal coupling gives nu with Phi(nu) = F(mu).  Hence
(1).

(ii) The point of keeping M as a variable: int c dnu is now LINEAR in nu.  So is
mu_nu itself, as a linear image of nu, hence so are L(mu_nu), E_{mu_nu}[p] and
G(mu_nu).  By Theorem A,

    Q(mu_nu) = G(mu_nu)^2/ln2 + Psi_Q(mu_nu),
    Psi_Q(mu) = -[ A(mu)^2 + int int T dmu dmu ] / ln2,

and Psi_Q is concave on signed measures -- minus a sum of squares of LINEAR
functionals -- hence concave in nu (concave composed with linear).  Therefore

    Phi(nu) = (1-alpha) G(nu)^2 / ln2  +  [ concave in nu ],

with the single convex direction G.  The convex coupling term has become linear,
which is exactly what the reformulation buys.

(iii) Freeze G(nu) = gamma.  On K_gamma = { nu in P(Delta) : E[p] <= t,
G(nu) = gamma } -- convex and weak-* compact -- Phi differs from a concave
continuous function by a constant (Q and L are weak-* continuous, and int c dnu
is too since c is continuous on Delta), so Bauer's minimum principle puts the
minimum at an extreme point.  K_gamma is cut from the positive cone by three
moment conditions (mass, the mean inequality, G), so its extreme points are
atomic with at most 3 atoms.  Every feasible nu lies in some K_gamma.       []

SCOPE, STATED PLAINLY.  This is WEAKER than the retracted claim.  Three pairs are
6 coordinates plus 2 free weights = 8 parameters, against Yu's 5, so B' does not
match the dimension of the literature's family; it is a valid finite reduction
where the literature has an invalid one, and nothing more.
SUPERSEDED: Theorem B''' (margin_lemma.py) reaches 2 orbits = 5 parameters, by
freezing B = 1-mean instead of G.  B' is kept because its proof is the one that
first isolated the pair-orbit reformulation.
"""

import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from entropy import h, mp, mpf
from mpmath import log

LN2 = log(2)
LOG2 = np.log(2.0)
A_OBS = 0.078877292705923173412
B_OBS = 0.32945473850303697239
C_STAR = 0.38234553336670272115
ALPHA = 0.0356069


def H(u):
    u = mpf(u)
    if u <= 0 or u >= 1:
        return mpf(0)
    return -(u * log(u) + (1 - u) * log(1 - u))


def T(u):
    u = mpf(u)
    return mpf(1) if u >= 1 else u + (1 - u) * log(1 - u)


def a_fn(x):
    x = mpf(x)
    return mpf(0) if x <= 0 else -x * log(x)


def g_fn(x):
    x = mpf(x)
    return mpf(0) if x <= 0 else x * (1 - log(x))


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def sstar(p, r):
    return np.clip(0.5, np.maximum(p, r), np.minimum(p + r, 1.0))


def cmin(atoms, w):
    """Exact min over symmetric couplings, by LP."""
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


# ---------------------------------------------------------------------------
print("LEMMA 1 -- T(u) = u + (1-u)ln(1-u) = sum_{n>=2} u^n/(n(n-1))")
for uv in ("0.05", "0.25", "0.5", "0.75", "0.9"):
    u = mpf(uv)
    ser = mp.nsum(lambda n: u ** n / (n * (n - 1)), [2, mp.inf])
    assert abs(T(u) - ser) < mpf(10) ** -30, (uv,)
print("   closed form = series at five points, to 1e-30                  [OK]")
one = mp.nsum(lambda n: mpf(1) / (n * (n - 1)), [2, mp.inf])
assert abs(one - 1) < mpf(10) ** -30
print(f"   sum 1/(n(n-1)) = {mp.nstr(one, 18)} = T(1)                    [OK]")
for uv in ("0.05", "0.3819660112501051", "0.5", "0.7", "0.95"):
    u = mpf(uv)
    assert abs(H(u) - (-u * log(u) + u - T(u))) < mpf(10) ** -32
print("   H(u) = -u ln u + u - T(u)                                      [OK]")

# ---------------------------------------------------------------------------
print("\nTHEOREM A -- symbolic verification of the quadratic-form identity")
k = 3
w = sp.symbols("w1 w2 w3", positive=True)
xs = sp.symbols("x1 x2 x3", positive=True)
Hs = lambda u: -u * sp.log(u) - (1 - u) * sp.log(1 - u)
Ts = lambda u: u + (1 - u) * sp.log(1 - u)
lhs = sum(w[i] * w[j] * Hs(xs[i] * xs[j]) for i in range(k) for j in range(k))
A_s = sum(w[i] * (-xs[i] * sp.log(xs[i])) for i in range(k))
B_s = sum(w[i] * xs[i] for i in range(k))
TT = sum(w[i] * w[j] * Ts(xs[i] * xs[j]) for i in range(k) for j in range(k))
rhs = (A_s + B_s) ** 2 - A_s**2 - TT
diff = sp.simplify(sp.expand(sp.expand_log(sp.expand(lhs - rhs), force=True)))
print(f"   LHS - RHS for a symbolic 3-atom measure: {diff}")
assert diff == 0
print("   identity is an ALGEBRAIC TAUTOLOGY -- exact, no truncation     [OK]")

rng = np.random.default_rng(11)
worst = mpf(0)
for _ in range(200):
    m = int(rng.integers(2, 7))
    xv = [mpf(float(v)) for v in rng.uniform(0.001, 0.999, m)]
    wv = [mpf(float(v)) for v in rng.normal(0, 1, m)]
    lhs_n = sum(wv[i] * wv[j] * H(xv[i] * xv[j]) for i in range(m) for j in range(m))
    A = sum(wv[i] * a_fn(xv[i]) for i in range(m))
    B = sum(wv[i] * xv[i] for i in range(m))
    TTn = sum(wv[i] * wv[j] * T(xv[i] * xv[j]) for i in range(m) for j in range(m))
    worst = max(worst, abs(lhs_n - ((A + B) ** 2 - A**2 - TTn)))
print(f"   worst error over 200 random SIGNED measures: {mp.nstr(worst, 6)} [OK]")
assert worst < mpf(10) ** -28

for n in (40, 120, 300):
    xg = np.linspace(1e-9, 1.0 - 1e-9, n)
    Hm = np.array([[float(H(mpf(float(i)) * mpf(float(j)))) for j in xg] for i in xg])
    Tm = np.array([[float(T(mpf(float(i)) * mpf(float(j)))) for j in xg] for i in xg])
    eH, eT = np.linalg.eigvalsh(Hm), np.linalg.eigvalsh(Tm)
    tol = 1e-9 * n
    assert int((eH > tol).sum()) == 1
    assert int((eT < -tol).sum()) == 0
    print(f"   n = {n:>4}: #positive eig of H(xy) = 1, T(xy) PSD           [OK]")

n = 200
xg = np.linspace(1e-9, 1 - 1e-9, n)
Hm = np.array([[float(H(mpf(float(i)) * mpf(float(j)))) for j in xg] for i in xg])
Tm = np.array([[float(T(mpf(float(i)) * mpf(float(j)))) for j in xg] for i in xg])
gv = np.array([float(g_fn(mpf(float(i)))) for i in xg])
av = np.array([float(a_fn(mpf(float(i)))) for i in xg])
err = np.abs(np.outer(gv, gv) - np.outer(av, av) - Tm - Hm).max()
print(f"   max |H(xy) - (gg^T - aa^T - T)| on a {n}x{n} grid: {err:.3e}   [OK]")
assert err < 1e-12

# ---------------------------------------------------------------------------
print("\nRETRACTION -- C is CONVEX, not concave.  Explicit witness:")
C1 = cmin(np.array([0.3]), np.array([1.0]))
C2 = cmin(np.array([0.6]), np.array([1.0]))
Cm = cmin(np.array([0.3, 0.6]), np.array([0.5, 0.5]))
print(f"   C(delta_0.3)                 = {C1:.15f}")
print(f"   C(delta_0.6)                 = {C2:.15f}")
print(f"   C(half-half mixture)         = {Cm:.15f}")
print(f"   chord                        = {0.5 * C1 + 0.5 * C2:.15f}")
print(f"   C(mix) - chord               = {Cm - (0.5 * C1 + 0.5 * C2):+.4e}")
assert Cm < 0.5 * C1 + 0.5 * C2 - 1e-6
print("   BELOW the chord => convex direction; concavity REFUTED         [OK]")

worst_v = (0.0, None)
for _ in range(1500):
    p, q = rng.uniform(0, 1, 2)
    th = float(rng.uniform(0.1, 0.9))
    d = (cmin(np.array([p, q]), np.array([th, 1 - th]))
         - th * cmin(np.array([p]), np.array([1.0]))
         - (1 - th) * cmin(np.array([q]), np.array([1.0])))
    if d < worst_v[0]:
        worst_v = (d, (p, q, th))
print(f"   worst C(mix) - chord over 1500 random segments: {worst_v[0]:+.4e}")
print(f"     at (p,q,theta) = ({worst_v[1][0]:.5f}, {worst_v[1][1]:.5f}, "
      f"{worst_v[1][2]:.5f})")
assert worst_v[0] < -0.1
print("   MK duality agrees: C = sup_phi 2 int phi dmu is a sup of linear")
print("   functionals, hence convex.  Bauer cannot be applied to alpha*C. [OK]")

# ---------------------------------------------------------------------------
print("\nTHEOREM B' -- what IS concave, in the pair-orbit variable nu")


def psi_q(atoms, weights):
    xv = [1 - mpf(p) for p in atoms]
    wv = [mpf(v) for v in weights]
    A = sum(wv[i] * a_fn(xv[i]) for i in range(len(xv)))
    TTv = sum(wv[i] * wv[j] * T(xv[i] * xv[j])
              for i in range(len(xv)) for j in range(len(xv)))
    return -(A**2 + TTv) / LN2


def q_direct(atoms, weights):
    return sum(mpf(weights[i]) * mpf(weights[j])
               * h(mpf(atoms[i]) + mpf(atoms[j]) - mpf(atoms[i]) * mpf(atoms[j]))
               for i in range(len(atoms)) for j in range(len(atoms)))


def g_of(atoms, weights):
    return sum(mpf(weights[i]) * g_fn(1 - mpf(atoms[i])) for i in range(len(atoms)))


print("\n   (i) Q = G^2/ln2 + Psi_Q in the p variable:")
for atoms, wts in (([0.2, 0.6], [0.5, 0.5]), ([0.0, 1.0], [0.618, 0.382]),
                   ([0.1, 0.35, 0.9], [0.2, 0.5, 0.3])):
    err2 = abs(q_direct(atoms, wts) - (g_of(atoms, wts) ** 2 / LN2
                                       + psi_q(atoms, wts)))
    print(f"     atoms {str(atoms):<22} err {mp.nstr(err2, 3)}")
    assert err2 < mpf(10) ** -28

print("\n   (ii) Psi_Q concave (midpoint test, 400 random segments)")
gap = mpf(0)
for _ in range(400):
    m = int(rng.integers(1, 5))
    atoms = ([mpf(float(v)) for v in rng.uniform(0, 1, m)]
             + [mpf(float(v)) for v in rng.uniform(0, 1, m)])
    wa = [mpf(float(v)) for v in rng.dirichlet(np.ones(m))] + [mpf(0)] * m
    wb = [mpf(0)] * m + [mpf(float(v)) for v in rng.dirichlet(np.ones(m))]
    th = mpf(float(rng.uniform(0.05, 0.95)))
    wm = [th * wa[i] + (1 - th) * wb[i] for i in range(2 * m)]
    gap = min(gap, psi_q(atoms, wm)
              - (th * psi_q(atoms, wa) + (1 - th) * psi_q(atoms, wb)))
print(f"     worst (Psi_Q(mid) - chord) = {mp.nstr(gap, 6)}  (>= 0)       [OK]")
assert gap >= -mpf(10) ** -25

print("\n   (iii) in the nu variable the coupling term is LINEAR:")
print("     int c dnu is linear by construction, so the only convex")
print("     direction left is G, which Bauer-freezing removes.")


def orbit_parts(pairs, wts, alpha=ALPHA):
    at, wt = [], []
    for (p, q), wi in zip(pairs, wts):
        at += [float(p), float(q)]
        wt += [wi / 2, wi / 2]
    at, wt = np.array(at), np.array(wt)
    P, R = np.meshgrid(at, at, indexing="ij")
    Q = float(wt @ hf(P + R - P * R) @ wt)
    L = float(wt @ hf(at))
    Cc = float(sum(wi * hf(sstar(np.array([p]), np.array([q])))[0]
                   for (p, q), wi in zip(pairs, wts)))
    mean = float(sum(wi * (p + q) / 2 for (p, q), wi in zip(pairs, wts)))
    return (1 - alpha) * Q + alpha * Cc - L, mean, L


pairs = [(B_OBS, B_OBS), (B_OBS, 1.0)]
wts = [1 - 2 * A_OBS, 2 * A_OBS]
val, mean, L = orbit_parts(pairs, wts)
print(f"\n   the obstruction law as a pair-orbit measure nu*:")
print(f"     {{b,b}} @ {wts[0]:.12f}   {{b,1}} @ {wts[1]:.12f}")
print(f"     mean(mu_nu) = {mean:.18f}")
print(f"     c*          = {C_STAR:.18f}")
assert abs(mean - C_STAR) < 1e-15
print(f"     Phi(nu*)    = {val:+.3e}   ratio - 1 = {val / L:+.3e}")
assert abs(val) < 1e-12
print("   only 2 pair-orbits, so B' (<= 3) covers the obstruction        [OK]")

print("\nALL PASS")
print("\nSUMMARY")
print("  H(xy) = g(x)g(y) - a(x)a(y) - T(xy),  T PSD, g(x) = x(1 - ln x)")
print("  => the OR-entropy kernel has AT MOST ONE POSITIVE SQUARE")
print("  C is CONVEX (retraction), so the reduction must keep the coupling")
print("  as a variable: <= 3 PAIR-ORBITS, hence <= 6 marginal atoms.")
