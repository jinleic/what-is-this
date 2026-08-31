"""shift-group orbit reduction (valid for Lambda = Id).

The instance with Lambda_i = Id is invariant under independent index shifts
a_i -> a_i + 1 mod s_i on each axis: on the S_i axis, index shift ai -> ai+1
is field multiplication alpha -> sigma*alpha (sigma a generator of the order-s
subgroup), and f(x) -> f(sigma^{-1} x) maps RS(S, t) to itself (deg preserved),
with Lambda = Id absorbing the scaling. Hence the group
G = Z_{s0} x Z_{s1} x Z_{s2} (order N) acts on points, lines, line-support
tuples, preserving wt, ell_i, delta, and W(B)-membership structure.
So delta(M) and r(M) are constant on G-orbits, and w_min(W(B)) is constant on
orbit classes of B-triples. Orbit representatives are enumerated in
lexicographic index space.
"""
from __future__ import annotations

import itertools


def shift(v: tuple, s: tuple, k: tuple) -> tuple:
    return tuple((v[i] + k[i]) % s[i] for i in range(3))


def orbit(v: tuple, s: tuple) -> list[tuple]:
    return [shift(v, s, k) for k in itertools.product(*[range(x) for x in s])]


def orbit_rep(v: tuple, s: tuple) -> tuple:
    return min(orbit(v, s))


def canonical_points(pts: list[tuple], s: tuple) -> tuple:
    """Canonical form of a point-set under G: min over images of sorted tuple."""
    best = None
    for k in itertools.product(*[range(x) for x in s]):
        img = tuple(sorted(shift(p, s, k) for p in pts))
        if best is None or img < best:
            best = img
    return best


def point_orbits(s: tuple) -> dict[tuple, int]:
    """Map each point to its orbit id (canonical rep)."""
    reps = {}
    for a in itertools.product(*[range(x) for x in s]):
        reps[a] = orbit_rep(a, s)
    return reps
