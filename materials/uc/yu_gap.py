"""Does Yu's decomposition INEQUALITY fail on genuine extreme points?

Yu (arXiv:2212.00658, Sec. 4) writes P_pq as a Caratheodory combination of
extreme points Q_i of

    P_B = { symmetric P_pq on B^2 : E[p] <= t },

and then uses concavity of P_pq -> g(P_pq, alpha) to assert

    g(sum_i gamma_i Q_i, alpha)  >=  sum_i gamma_i g(Q_i, alpha).        (*)

concavity.py already refutes the CONCAVITY CLAIM itself.  But refuting the
cited justification is weaker than refuting the inference (*): (*) could still
hold for the particular decompositions that arise.  This file tests (*) directly
on Yu's own extreme points.

YU'S EXTREME POINTS, verbatim shape:
    Q = (1-beta) Q_{a1,a2} + beta Q_{b1,b2},   Q_{x,y} = (d_(x,y)+d_(y,x))/2,
    0 <= a := (a1+a2)/2 <= t < b := (b1+b2)/2 <= 1,
    beta = 0   or   beta = (t-a)/(b-a) > 0.
A beta>0 extreme point has marginal mean EXACTLY t, since
(1-beta) a + beta b = a + (t-a) = t.

TWO POINT MASSES CANNOT BREAK IT.  For Q_1 = d_(a,a), Q_2 = d_(c,c) the defect
of (*) is  -gamma(1-gamma)(1-alpha)[h(2a-a^2)+h(2c-c^2)-2h(a+c-ac)], and in
x = 1-a, y = 1-c that bracket is h(x^2)+h(y^2)-2h(xy), whose maximum over
[1-t,1]^2 is 0, attained only at x=y.  So any genuine failure of (*) must use
an extreme point with beta > 0, i.e. one carrying an atom above t.
"""

import numpy as np
from scipy.optimize import minimize

LOG2 = np.log(2.0)
T_STAR = 0.382345533366702721
PSI = 0.3819660112501051517954132
ALPHA_STAR = 0.0356069


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def maxent_or(p, r):
    return float(np.clip(0.5, max(p, r), min(p + r, 1.0)))


class Ext:
    """A Yu extreme point: atoms with joint weights, plus its marginal."""

    def __init__(self, a1, a2, b1, b2, beta):
        self.pairs = []                      # (p, r, weight) on the joint
        if beta < 1.0:
            self.pairs += [(a1, a2, (1 - beta) / 2), (a2, a1, (1 - beta) / 2)]
        if beta > 0.0:
            self.pairs += [(b1, b2, beta / 2), (b2, b1, beta / 2)]
        self.atoms = np.array([a1, a2, b1, b2])
        self.w = np.array([(1 - beta) / 2, (1 - beta) / 2, beta / 2, beta / 2])

    def marg_mean(self):
        return float(self.w @ self.atoms)

    def Eh(self):
        return float(self.w @ hf(self.atoms))

    def coupled(self):
        return sum(wt * float(hf(maxent_or(p, r))) for p, r, wt in self.pairs)


def iid_term(atoms, w):
    a = np.asarray(atoms, float); w = np.asarray(w, float)
    orm = a[:, None] + a[None, :] - a[:, None] * a[None, :]
    return float(w @ hf(orm.ravel()).reshape(len(w), len(w)) @ w)


def g_of(atoms, w, coupled, alpha):
    return (1 - alpha) * iid_term(atoms, w) + alpha * coupled


def defect(Q1, Q2, gamma, alpha):
    """g(mixture) - [gamma g(Q1) + (1-gamma) g(Q2)];  < 0 refutes (*)."""
    atoms = np.concatenate([Q1.atoms, Q2.atoms])
    w = np.concatenate([gamma * Q1.w, (1 - gamma) * Q2.w])
    cpl = gamma * Q1.coupled() + (1 - gamma) * Q2.coupled()   # linear, exact
    g_mix = g_of(atoms, w, cpl, alpha)
    g_avg = gamma * g_of(Q1.atoms, Q1.w, Q1.coupled(), alpha) \
        + (1 - gamma) * g_of(Q2.atoms, Q2.w, Q2.coupled(), alpha)
    return g_mix - g_avg


print("1. Confirm two point masses CANNOT break (*)")
lo = 1 - T_STAR
g = np.linspace(lo, 1.0, 1200)
worst = -np.inf
for x in g:
    d = hf(x * x) + hf(g * g) - 2 * hf(x * g)
    worst = max(worst, float(d.max()))
print(f"   max over [1-t,1]^2 of h(x^2)+h(y^2)-2h(xy) = {worst:+.3e}   [<= 0]")
assert worst < 1e-9

print("\n2. Search for a violation of (*) using genuine extreme points")
print("   Q1: beta=0 (a1,a2 <= t on average).  Q2: beta>0, mean exactly t.")


def build(z, t):
    """z -> (Q1, Q2, gamma) respecting Yu's extreme-point constraints."""
    a1, a2, c1, c2, d1, d2, gam = z
    # Q1: beta = 0, needs (a1+a2)/2 <= t
    s = (a1 + a2) / 2
    if s > t:
        a1 *= t / s * 0.999999
        a2 *= t / s * 0.999999
    Q1 = Ext(a1, a2, 0.0, 0.0, 0.0)
    # Q2: beta > 0, needs a=(c1+c2)/2 <= t < b=(d1+d2)/2
    a = (c1 + c2) / 2
    b = (d1 + d2) / 2
    if not (a <= t < b):
        return None
    beta = (t - a) / (b - a)
    if not (0.0 < beta <= 1.0):
        return None
    Q2 = Ext(c1, c2, d1, d2, beta)
    return Q1, Q2, float(np.clip(gam, 1e-6, 1 - 1e-6))


best = (np.inf, None)
rng = np.random.default_rng(0)
for trial in range(30000):
    z = rng.uniform(0, 1, 7)
    z[4] = rng.uniform(T_STAR, 1.0)          # push d1 above t
    z[5] = rng.uniform(T_STAR, 1.0)
    got = build(z, T_STAR)
    if got is None:
        continue
    Q1, Q2, gam = got
    d = defect(Q1, Q2, gam, ALPHA_STAR)
    if d < best[0]:
        best = (d, (Q1, Q2, gam, z.copy()))

print(f"   best (most negative) defect from random search: {best[0]:+.6e}")


def polish(z0, t, alpha):
    def obj(z):
        got = build(z, t)
        if got is None:
            return 1.0
        Q1, Q2, gam = got
        return defect(Q1, Q2, gam, alpha)

    r = minimize(obj, z0, method="Nelder-Mead",
                 options={"maxiter": 6000, "xatol": 1e-12, "fatol": 1e-15})
    return r


if best[1] is not None:
    r = polish(best[1][3], T_STAR, ALPHA_STAR)
    got = build(r.x, T_STAR)
    if got is not None:
        Q1, Q2, gam = got
        d = defect(Q1, Q2, gam, ALPHA_STAR)
        print(f"   after polishing: defect = {d:+.6e}")
        if d < -1e-10:
            print("\n   *** (*) FAILS on genuine extreme points ***")
            print(f"   Q1 (beta=0): atoms {np.array2string(Q1.atoms[:2], precision=9)}"
                  f"  mean {Q1.marg_mean():.9f}  Eh {Q1.Eh():.6f}")
            b1, b2 = Q2.atoms[2], Q2.atoms[3]
            a1, a2 = Q2.atoms[0], Q2.atoms[1]
            bb = (b1 + b2) / 2
            aa = (a1 + a2) / 2
            beta = (T_STAR - aa) / (bb - aa)
            print(f"   Q2: a=({a1:.9f},{a2:.9f}) b=({b1:.9f},{b2:.9f}) "
                  f"beta={beta:.9f}")
            print(f"       mean {Q2.marg_mean():.9f} (must equal t) "
                  f"Eh {Q2.Eh():.6f}")
            print(f"   gamma = {gam:.9f}")
            print(f"   g(mix)            = "
                  f"{g_of(np.concatenate([Q1.atoms,Q2.atoms]), np.concatenate([gam*Q1.w,(1-gam)*Q2.w]), gam*Q1.coupled()+(1-gam)*Q2.coupled(), ALPHA_STAR):.12f}")
            print(f"   gamma g(Q1)+..    = "
                  f"{gam*g_of(Q1.atoms,Q1.w,Q1.coupled(),ALPHA_STAR)+(1-gam)*g_of(Q2.atoms,Q2.w,Q2.coupled(),ALPHA_STAR):.12f}")
        else:
            print("\n   no violation of (*) found on extreme points")
            print("   => the inference step may hold even though concavity does not")

print("\n3. Simplest structured probe: Q1 = d_(a,a), Q2 = (1-beta)d_(c,c)+beta Q_{d,1}")
print(f"   {'a':>9} {'c':>9} {'d':>9} {'gamma':>7} {'defect':>13}")
rows = 0
worst_struct = (np.inf, None)
for a in (0.05, 0.15, 0.25, 0.32, 0.38):
    for c in (0.10, 0.20, 0.30, 0.3294547385):
        for d in (0.5, 0.7, 0.9, 1.0):
            bb = (d + 1.0) / 2
            if not (c <= T_STAR < bb):
                continue
            beta = (T_STAR - c) / (bb - c)
            if not (0 < beta <= 1):
                continue
            Q1 = Ext(a, a, 0, 0, 0.0)
            Q2 = Ext(c, c, d, 1.0, beta)
            for gam in (0.25, 0.5, 0.75):
                dd = defect(Q1, Q2, gam, ALPHA_STAR)
                if dd < worst_struct[0]:
                    worst_struct = (dd, (a, c, d, gam))
                if dd < -1e-12 and rows < 12:
                    print(f"   {a:>9.4f} {c:>9.4f} {d:>9.4f} {gam:>7.2f} {dd:>13.3e}")
                    rows += 1
print(f"   most negative structured defect: {worst_struct[0]:+.6e} at "
      f"{worst_struct[1]}")

print("\nDONE")
