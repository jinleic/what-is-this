"""Does the concavity failure survive away from the degenerate {0,1} corner?

concavity.py established, with a hand-verified witness, that

    F(P_p) := E_{P_p x P_p} h(p + q - pq)

is NOT concave on { E[p] <= t }, even at t = psi and t = c*.  Exact defect for
mu_1 = (1-t) d_0 + t d_1  and  mu_2 = (1-w) d_0 + w d_{p2}:

    F(mu_1) = 0,   defect = (w/4) [ w h(2 p2 - p2^2) + 2 h(p2) (t - w) ].

But that witness has E_{mu_1}[h(p)] = 0: mu_1 is supported on {0,1}, which is
exactly the degenerate case Liu's Lemma 7 quarantines and which Yu excludes via
the side condition E h(p) > 0.  So the question that decides whether the gap is
fatal or cosmetic is:

    does concavity fail for laws with E[h(p)] bounded AWAY from 0?

We add the constraint E[h(p)] >= eta to BOTH measures and re-run the search.
"""

import numpy as np
from scipy.optimize import minimize

LOG2 = np.log(2.0)
T_STAR = 0.382345533366702721
PSI = 0.3819660112501051517954132


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


def F_p(a, w):
    a = np.asarray(a, float); w = np.asarray(w, float)
    orm = a[:, None] + a[None, :] - a[:, None] * a[None, :]
    return float(w @ hf(orm.ravel()).reshape(len(w), len(w)) @ w)


def Eh(a, w):
    return float(np.asarray(w, float) @ hf(a))


print("1. Closed-form defect for the degenerate witness, as a function of eta")
print("   mu_1 = (1-t) d_0 + t d_1  has  E h(p) = 0  exactly.")
t = T_STAR
for p2, w in [(0.2179, 0.531), (0.25, 0.5), (0.3, 0.45)]:
    d = (w / 4) * (w * float(hf(2 * p2 - p2**2)) + 2 * float(hf(p2)) * (t - w))
    print(f"   p2 = {p2:<7} w = {w:<6}: defect = {d:+.6f}   "
          f"(mu_2 has E h(p) = {w*float(hf(p2)):.4f})")

print("\n2. Search with E[h(p)] >= eta imposed on BOTH measures")
print(f"   {'eta':>8} {'t':>14} {'max defect':>13}  verdict")


def run(t, eta, ks=(1, 2, 3), tries=250, seed=5):
    rng = np.random.default_rng(seed)
    best, barg = -np.inf, None
    for k in ks:
        def split(z):
            a1 = 1 / (1 + np.exp(-np.clip(z[0:k], -30, 30)))
            w1 = np.exp(z[k:2 * k] - z[k:2 * k].max()); w1 /= w1.sum()
            a2 = 1 / (1 + np.exp(-np.clip(z[2 * k:3 * k], -30, 30)))
            w2 = np.exp(z[3 * k:4 * k] - z[3 * k:4 * k].max()); w2 /= w2.sum()
            return a1, w1, a2, w2

        def neg(z):
            a1, w1, a2, w2 = split(z)
            mid_a = np.concatenate([a1, a2]); mid_w = np.concatenate([w1, w2]) / 2
            return -((F_p(a1, w1) + F_p(a2, w2)) / 2 - F_p(mid_a, mid_w))

        cons = [
            {"type": "ineq", "fun": lambda z: t - float(split(z)[1] @ split(z)[0])},
            {"type": "ineq", "fun": lambda z: t - float(split(z)[3] @ split(z)[2])},
            {"type": "ineq", "fun": lambda z: Eh(split(z)[0], split(z)[1]) - eta},
            {"type": "ineq", "fun": lambda z: Eh(split(z)[2], split(z)[3]) - eta},
        ]
        for _ in range(tries):
            z0 = rng.normal(0, 2.5, 4 * k)
            try:
                r = minimize(neg, z0, method="SLSQP", constraints=cons,
                             options={"maxiter": 250, "ftol": 1e-14})
            except Exception:
                continue
            if r.success and -r.fun > best:
                ok = all(c["fun"](r.x) > -1e-8 for c in cons)
                if ok:
                    best, barg = -float(r.fun), split(r.x)
    return best, barg


for eta in (0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8):
    best, barg = run(T_STAR, eta)
    verdict = "CONCAVE" if best < 1e-9 else "FAILS"
    print(f"   {eta:>8.2f} {T_STAR:>14.10f} {best:>13.3e}  {verdict}")
    if best > 1e-9 and eta >= 0.4:
        a1, w1, a2, w2 = barg
        print(f"      mu1 {np.array2string(a1, precision=6)} "
              f"w {np.array2string(w1, precision=4)} "
              f"mean {w1@a1:.6f} Eh {Eh(a1,w1):.4f}")
        print(f"      mu2 {np.array2string(a2, precision=6)} "
              f"w {np.array2string(w2, precision=4)} "
              f"mean {w2@a2:.6f} Eh {Eh(a2,w2):.4f}")

print("\n3. The bottom line: can a feasible law drive the RATIO below 1?")
print("   Concavity is only a proof device.  What matters for soundness is")
print("   whether inf over ALL feasible laws of the two-strategy slack is >= 0")
print("   at t = c*.  That is tested directly in robust.py; here we check the")
print("   specific laws that break concavity.")

print(f"\n   {'law':<46} {'mean':>9} {'E h(p)':>9} {'iid ratio':>11}")
for a, w, lbl in [
    (np.array([0.0, 1.0]), np.array([1 - T_STAR, T_STAR]), "(1-t)d_0 + t d_1  [degenerate]"),
    (np.array([0.0, 0.2179]), np.array([0.469, 0.531]), "0.469 d_0 + 0.531 d_0.2179"),
    (np.array([T_STAR]), np.array([1.0]), "d_{c*}"),
    (np.array([1.0, 0.329454738503037]),
     np.array([0.0788772927059232, 1 - 0.0788772927059232]), "obstruction {1, b}"),
]:
    m = float(w @ a); e = Eh(a, w); f = F_p(a, w)
    ratio = f / e if e > 1e-12 else float("nan")
    print(f"   {lbl:<46} {m:>9.6f} {e:>9.6f} {ratio:>11.6f}")

print("\n   A ratio < 1 for the iid term alone is expected above psi; the")
print("   coupled term is what must rescue it.  See robust.py.")
print("\nDONE")
