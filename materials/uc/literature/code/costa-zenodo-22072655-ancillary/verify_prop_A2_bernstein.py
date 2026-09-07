#!/usr/bin/env python3
"""Exact rational verification of the Bernstein coefficients in Proposition A.2."""
from fractions import Fraction as F
from math import comb

# Polynomials are dictionaries (i,j) -> coefficient of x^i y^j.
def add(P, Q):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, F(0)) + v
        if R[k] == 0:
            del R[k]
    return R

def scale(P, c):
    return {k: c*v for k, v in P.items() if c*v}

def mul(P, Q):
    R = {}
    for (i,j), a in P.items():
        for (k,l), b in Q.items():
            key = (i+k, j+l)
            R[key] = R.get(key, F(0)) + a*b
    return {k:v for k,v in R.items() if v}

one = {(0,0): F(1)}
x = {(1,0): F(1)}
y = {(0,1): F(1)}
xy = mul(x,y)

# W=2-x-y+xy, U=2-2x-y+2xy, V=2-x-2y+2xy.
W = add(add(add({(0,0):F(2)}, scale(x,F(-1))), scale(y,F(-1))), xy)
U = add(add(add({(0,0):F(2)}, scale(x,F(-2))), scale(y,F(-1))), scale(xy,F(2)))
V = add(add(add({(0,0):F(2)}, scale(x,F(-1))), scale(y,F(-2))), scale(xy,F(2)))
r = mul(xy, W)
Pxy = add(mul(mul(W,U),V), add(scale(one,F(-1)), r))  # WUV-(1-r)

# Substitute x=(43/64)s, y=(43/64)t.
a = F(43,64)
power = {}
for (i,j), c in Pxy.items():
    power[(i,j)] = c * a**(i+j)

max_i = max(i for i,j in power)
max_j = max(j for i,j in power)
assert (max_i, max_j) == (3,3)

# Power-to-Bernstein conversion in degree n=m=3:
# s^p = sum_{i=p}^3 C(i,p)/C(3,p) B_i^3(s), and similarly for t.
bern = [[F(0) for _ in range(4)] for _ in range(4)]
for i in range(4):
    for j in range(4):
        total = F(0)
        for p in range(i+1):
            for q in range(j+1):
                c = power.get((p,q), F(0))
                if c:
                    total += c * F(comb(i,p), comb(3,p)) * F(comb(j,q), comb(3,q))
        bern[i][j] = total

expected = [
    [F(7), F(41,12), F(8221,6144), F(20653,131072)],
    [F(41,12), F(18803,9216), F(402175,393216), F(1894129,6291456)],
    [F(8221,6144), F(402175,393216), F(52130497,75497472), F(587913113,1610612736)],
    [F(20653,131072), F(1894129,6291456), F(587913113,1610612736), F(6389233321,17179869184)],
]

print("bidegree: (3, 3)")
print("Bernstein coefficient matrix:")
for row in bern:
    print("  " + "  ".join(str(v) for v in row))

assert bern == expected, "computed coefficient matrix does not match Proposition A.2"
assert all(v > 0 for row in bern for v in row), "a Bernstein coefficient is not positive"
print("all 16 coefficients are strictly positive")
print("PROPOSITION A.2 CERTIFICATE PASSED")
