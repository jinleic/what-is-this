"""Is the entropy kernel PSD?  The decisive question for rigour.

WHY THIS MATTERS.  The union-closed variational problem has a quadratic term

    Q(mu) = E_{p,q ~ mu x mu} h(p + q - pq).

A quadratic form cannot be bounded below by a linear functional in general, so
every attempt to certify inf_mu >= 0 rigorously wants the linearisation

    Q(mu) >= 2 B(mu_0, mu) - B(mu_0, mu_0),

which is valid EXACTLY WHEN the kernel is positive semidefinite on differences
of measures, since the gap is B(mu - mu_0, mu - mu_0).  Liu's Theorem 13 is
conditional on a positive-semidefiniteness hypothesis; this file asks the
underlying question directly.

CHANGE OF VARIABLE.  Put x = 1-p, y = 1-q.  Then

    p + q - pq = 1 - (1-p)(1-q) = 1 - xy,     h(1-xy) = h(xy),

so the kernel is exactly  K(x,y) = h(xy)  on [0,1]^2.

SERIES DECOMPOSITION.  With natural log H(u) = -u ln u - (1-u) ln(1-u),
expanding ln(1/(1-u)) = sum_{k>=1} u^k / k gives the exact identity

    H(u) = u ln(1/u) + u - sum_{k>=2} u^k / (k(k-1)),

(checked at u=1, where sum_{k>=2} 1/(k(k-1)) = 1).  Substituting u = xy:

    H(xy) = [ x l(x) * y  +  x * y l(y) ]   <- rank 2, signature (1,1)
            + x y                            <- rank 1, PSD
            - sum_{k>=2} x^k y^k / (k(k-1))  <- PSD, entering NEGATIVELY

with l(x) = ln(1/x).  The last block is a positive combination of the rank-one
PSD kernels x^k y^k and enters with a MINUS sign, so the kernel is expected to
be indefinite.  We settle it numerically.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LOG2 = np.log(2.0)


def hf(x):
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    out = np.zeros_like(x)
    m = (x > 0) & (x < 1)
    xm = x[m]
    out[m] = -(xm * np.log(xm) + (1 - xm) * np.log1p(-xm)) / LOG2
    return out


print("1. Series identity  H(u) = u ln(1/u) + u - sum_{k>=2} u^k/(k(k-1))")
for u in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
    direct = -(u * np.log(u) + (1 - u) * np.log1p(-u))
    k = np.arange(2, 4000)
    series = u * np.log(1 / u) + u - np.sum(u**k / (k * (k - 1)))
    err = abs(direct - series)
    print(f"   u = {u:<5}: |direct - series| = {err:.3e}")
    assert err < 1e-12, (u, err)
print("   identity verified                                          [OK]")

print("\n2. Spectrum of the kernel K(x,y) = h(xy) on [0,1]")
print("   Gram matrix on a Chebyshev-like grid, weighted for quadrature.")

for n in (60, 120, 240):
    # Midpoint grid with uniform weights: the Gram matrix W^{1/2} K W^{1/2}
    # has the same signature as the integral operator's discretisation.
    xs = (np.arange(n) + 0.5) / n
    w = np.full(n, 1.0 / n)
    K = hf(np.outer(xs, xs).ravel()).reshape(n, n)
    G = np.sqrt(w)[:, None] * K * np.sqrt(w)[None, :]
    ev = np.linalg.eigvalsh(G)
    neg = ev[ev < -1e-12]
    print(f"   n = {n:>3}: lambda_min = {ev[0]:+.6e}   lambda_max = {ev[-1]:+.6e}"
          f"   #(eig < -1e-12) = {len(neg)}")

print("\n3. Explicit violating measure (certificate of indefiniteness)")
# Find a signed measure nu with sum nu = 0 and nu^T K nu < 0, i.e. a difference
# of two probability measures on which the linearisation FAILS.
n = 200
xs = (np.arange(n) + 0.5) / n
K = hf(np.outer(xs, xs).ravel()).reshape(n, n)
# Project onto the mean-zero subspace, then take the most negative direction.
P = np.eye(n) - np.ones((n, n)) / n
Kp = P @ K @ P
ev, evec = np.linalg.eigh(Kp)
v = evec[:, 0]
val = float(v @ K @ v)
print(f"   most negative mean-zero direction: nu^T K nu = {val:+.6e}")
print(f"   sum(nu) = {float(v.sum()):+.2e}  (mean-zero)")

if val < -1e-10:
    # Turn the eigenvector into an explicit difference of two probability
    # measures, so the failure is exhibited on genuine laws.
    pos = np.clip(v, 0, None); neg = np.clip(-v, 0, None)
    pos /= pos.sum(); neg /= neg.sum()
    d = pos - neg
    q = float(d @ K @ d)
    print(f"   normalised as mu_+ - mu_-: value = {q:+.6e}")
    sup_p = xs[pos > pos.max() * 0.05]
    sup_n = xs[neg > neg.max() * 0.05]
    print(f"   mu_+ concentrated near x in [{sup_p.min():.3f}, {sup_p.max():.3f}]")
    print(f"   mu_- concentrated near x in [{sup_n.min():.3f}, {sup_n.max():.3f}]")
    print("\n   => K(x,y) = h(xy) is NOT positive semidefinite.")
    print("      The naive linearisation of the quadratic term is INVALID;")
    print("      any rigorous certificate must avoid it or restrict the cone.")
else:
    print("\n   => no negative direction found; kernel may be PSD.")

# Minimal-support version: how few atoms suffice to see indefiniteness?
print("\n4. Smallest atom count exhibiting a negative direction")
found = None
for k in range(2, 9):
    best = 0.0
    rng = np.random.default_rng(0)
    for _ in range(20000):
        pts = np.sort(rng.uniform(0, 1, k))
        c = rng.normal(0, 1, k)
        c -= c.mean()
        Kk = hf(np.outer(pts, pts).ravel()).reshape(k, k)
        q = float(c @ Kk @ c) / float(c @ c)
        if q < best:
            best, arg = q, (pts.copy(), c.copy())
    print(f"   k = {k}: min normalised form = {best:+.6e}")
    if best < -1e-9 and found is None:
        found = (k, arg)

if found:
    k, (pts, c) = found
    print(f"\n   minimal witness at k = {k}:")
    print(f"     points  = {np.array2string(pts, precision=6)}")
    print(f"     weights = {np.array2string(c, precision=6)} (sum {c.sum():+.1e})")

print("\nDONE")
