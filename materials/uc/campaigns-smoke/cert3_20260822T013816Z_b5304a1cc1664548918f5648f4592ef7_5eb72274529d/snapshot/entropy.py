"""Binary entropy and the one-step union-closed functional.

Base-2 binary entropy throughout:

    h(x) = -x log2 x - (1-x) log2 (1-x),   h(0) = h(1) = 0.

Papers differ in log base (Gilmer/Sawin/Cambie/Yu: base 2; Wakhare/Boppana:
natural).  Every inequality here is homogeneous of degree 1 in h, so the base
cancels and the constants are unaffected.  We fix base 2 and never mix.

CONSTANTS.  With phi = (1+sqrt5)/2 the golden ratio,

    psi   = (3 - sqrt5)/2 = 0.381966...   the one-step barrier
    ohi   = (sqrt5 - 1)/2 = 1/phi = 1 - psi

psi is the nontrivial root of  x^2 - 3x + 1 = 0, equivalently of (1-x)^2 = x.
"""

from mpmath import mp, mpf, log, sqrt

mp.dps = 40


def h(x):
    """Base-2 binary entropy, with h(0) = h(1) = 0 exactly."""
    x = mpf(x)
    if x <= 0 or x >= 1:
        # Guard the endpoints: x log x -> 0 but mpmath would raise on log(0).
        if x == 0 or x == 1:
            return mpf(0)
        raise ValueError(f"h: argument {x} outside [0,1]")
    return -(x * log(x) + (1 - x) * log(1 - x)) / log(2)


PHI = (1 + sqrt(5)) / 2          # 1.618...  golden ratio
OHI = (sqrt(5) - 1) / 2          # 0.618...  = 1/PHI = 1 - PSI
PSI = (3 - sqrt(5)) / 2          # 0.381966... the one-step barrier


def union_prob(p, q):
    """Pr[(A u B)_i = 1] when Pr[A_i=1]=p, Pr[B_i=1]=q independently."""
    return p + q - p * q


def onestep_gain(atoms, weights):
    """G(mu) = E_{p,q~mu x mu} h(p+q-pq) - E_{p~mu} h(p).

    mu is the law of the conditional inclusion probability
    p_i(a) = Pr[A_i = 1 | A_{<i} = a].  The Gilmer-style argument needs
    G(mu) >= 0 for every mu obeying the mean constraint E_mu[p] <= t; the
    largest such t is the reach of the one-step method.
    """
    assert len(atoms) == len(weights)
    tot = sum(weights)
    assert abs(tot - 1) < mpf(10) ** (-mp.dps + 5), f"weights sum to {tot}"
    pair = sum(
        wi * wj * h(union_prob(atoms[i], atoms[j]))
        for i, wi in enumerate(weights)
        for j, wj in enumerate(weights)
    )
    solo = sum(w * h(a) for a, w in zip(atoms, weights))
    return pair - solo


def onestep_mean(atoms, weights):
    return sum(mpf(a) * w for a, w in zip(atoms, weights))


def sawin_lambda(u):
    """Sawin's sharp factor: E h(p+q-pq) >= lambda(u) E h(p) for E[p] <= u.

    Two branches, meeting at u = PSI where both equal 1 exactly.
    """
    u = mpf(u)
    if u <= PSI:
        return h(union_prob(u, u)) / h(u)
    return PHI * (1 - u)
