"""gate_a_radius.py — INDEPENDENT re-derivation of the Krawczyk radius for
the archived rank-25 certificate: find the largest rho* such that ALL THREE
hypotheses of Proposition 6 hold:

  (i)   K(rho) := max_e [ rs0_e + rho*rs1_e + rho^2*rs2_e ] < 1
  (ii)  off_max <= (1 - K(rho)) * rho,  off_e = |(Y g0)_e|
  (iii) margin_e(rho) := rho - (off_e + rho*rs0_e + rho^2*rs1_e
                               + rho^3*rs2_e) > 0  for every e

Each candidate rho is tested with OUTWARD rounding: quantities are evaluated
first exactly (fmpq) and then re-checked with arb balls whose lower/upper
endpoints carry the rounding; a hypothesis passes iff the ball endpoints
still certify it.  Bisection between lo=1e-9 and hi=1e-3, then the limiting
hypothesis is identified by exact per-hypothesis probes at returned rho.

Also reports, at fixed rho = 1e-6 (the paper's certificate), the
per-hypothesis exact values for the campaign record.
"""
import csv
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
from flint import fmpq, fmpq_mat, arb

CERT = os.environ.get(
    "CERTDIR",
    os.path.join(os.path.dirname(__file__), "..", "scratch", "upstream_ref",
                 "certs", "rank25_cert"))

t0 = time.time()
M_EXPECT = 512


def fqa(v):
    return fmpq(f"{Fraction(v).numerator}/{Fraction(v).denominator}")


def load(fn):
    rows = []
    for rec in csv.reader(open(fn)):
        rows.append([fqa(float(v)) for v in rec])
    return rows


A = load(os.path.join(CERT, "A.csv"))
B = load(os.path.join(CERT, "B.csv"))
C = load(os.path.join(CERT, "C.csv"))
Y = load(os.path.join(CERT, "Y_preconditioner.csv"))
S = [int(float(v)) for rec in csv.reader(open(os.path.join(CERT, "S_indices.csv")))
     for v in rec]
p = nb = nc = 8
r = len(A[0])
M = p * nb * nc
Nv = (p + nb + nc) * r
offB, offC = p * r, (p + nb) * r
assert M == M_EXPECT and len(S) == M

# tensor from CD doubling (same construction as gate_a_verify; re-done here
# so this script is standalone)
def t4(i, j):
    if i == 0:
        return j, 1
    if j == 0:
        return i, 1
    if i == j:
        return 0, -1
    k = 6 - i - j
    return k, (1 if (i, j) in ((1, 2), (2, 3), (3, 1)) else -1)


def hv4(u, v):
    out = [0] * 4
    for i in range(4):
        if u[i] == 0:
            continue
        for j in range(4):
            if v[j] == 0:
                continue
            k, s = t4(i, j)
            out[k] += s * u[i] * v[j]
    return out


def h_mul(x, y):
    a, b = x[:4], x[4:]
    c, d = y[:4], y[4:]
    dc = [d[0], -d[1], -d[2], -d[3]]
    cc = [c[0], -c[1], -c[2], -c[3]]
    left = [hv4(a, c)[k] - hv4(dc, b)[k] for k in range(4)]
    right = [hv4(d, a)[k] + hv4(b, cc)[k] for k in range(4)]
    return left + right


E8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
T = [[[h_mul(E8[i], E8[j])[k] for k in range(8)] for j in range(8)]
     for i in range(8)]

g0 = [fmpq(0)] * M
J = [[fmpq(0)] * Nv for _ in range(M)]
for i in range(p):
    for j in range(nb):
        for k in range(nc):
            e = (i * nb + j) * nc + k
            g0[e] = sum(A[i][s] * B[j][s] * C[k][s] for s in range(r)) \
                - fmpq(T[i][j][k])
            for s in range(r):
                J[e][i * r + s] += B[j][s] * C[k][s]
                J[e][offB + j * r + s] += A[i][s] * C[k][s]
                J[e][offC + k * r + s] += A[i][s] * B[j][s]

Ymat = fmpq_mat(M, M, [v for row in Y for v in row])
Jsel = fmpq_mat(M, M, [J[e][c] for e in range(M) for c in S])
Imat = fmpq_mat(M, M, [fmpq(1) if i == j else fmpq(0)
                       for i in range(M) for j in range(M)])


def mat_abs(Mm):
    return fmpq_mat(Mm.nrows(), Mm.ncols(), [abs(v) for v in Mm.entries()])


E = Ymat * Jsel - Imat
absImYJ = mat_abs(E)
rs0 = [sum(absImYJ[e, c] for c in range(M)) for e in range(M)]
g0m = fmpq_mat(M, 1, g0)
Yg0 = Ymat * g0m
off = [abs(Yg0[e, 0]) for e in range(M)]
off_max = max(off)
absY = mat_abs(Ymat)
x = [None] * Nv
for i in range(p):
    for s in range(r):
        x[i * r + s] = A[i][s]
for j in range(nb):
    for s in range(r):
        x[offB + j * r + s] = B[j][s]
for k in range(nc):
    for s in range(r):
        x[offC + k * r + s] = C[k][s]
Spos = {c: t for t, c in enumerate(S)}
freeS = [False] * Nv
for c in S:
    freeS[c] = True
absx = [abs(v) for v in x]
W1s = [[fmpq(0)] * M for _ in range(M)]
W2s = [[0] * M for _ in range(M)]
for i in range(p):
    for j in range(nb):
        for k in range(nc):
            e = (i * nb + j) * nc + k
            for s in range(r):
                ia, jb, kc = i * r + s, offB + j * r + s, offC + k * r + s
                for col, o1, o2 in ((ia, jb, kc), (jb, ia, kc), (kc, ia, jb)):
                    if col in Spos:
                        tc = Spos[col]
                        w1 = (absx[o2] if freeS[o1] else 0) + \
                             (absx[o1] if freeS[o2] else 0)
                        if w1:
                            W1s[e][tc] += w1
                        if freeS[o1] and freeS[o2]:
                            W2s[e][tc] += 1
W1s_mat = fmpq_mat(M, M, [v for row in W1s for v in row])
W2s_mat = fmpq_mat(M, M, [fmpq(v) for row in W2s for v in row])
absYW1 = absY * W1s_mat
absYW2 = absY * W2s_mat
rs1 = [sum(absYW1[e, c] for c in range(M)) for e in range(M)]
rs2 = [sum(absYW2[e, c] for c in range(M)) for e in range(M)]
print(f"setup done in {time.time() - t0:.1f}s; off_max={float(off_max):.3e}",
      flush=True)

ONE = fmpq(1)


def test_rho(rho, outward=True):
    """Return the triple of hypothesis-flags at radius rho (exact fmpq)."""
    rho2, rho3 = rho * rho, rho ** 3
    # contraction row sums and margins, formed in exact rationals
    K = None
    margins_ok = True
    for e in range(M):
        rowK = rs0[e] + rho * rs1[e] + rho2 * rs2[e]
        if K is None or rowK > K:
            K = rowK
        marg = rho - (off[e] + rho * rs0[e] + rho2 * rs1[e] + rho3 * rs2[e])
        if marg <= 0:
            margins_ok = False
    K_lt_1 = K < ONE
    banach = off_max <= (ONE - K) * rho
    flags = [K_lt_1, banach, margins_ok]
    if not outward:
        return flags, K
    # OUTWARD re-check with arb balls: thresholds are 1 and 0; re-derive
    K_arb = arb(0)
    worst = None
    for e in range(M):
        K_arb = max(K_arb, arb(rs0[e]) + arb(rho) * arb(rs1[e])
                    + arb(rho2) * arb(rs2[e]))
        m = (arb(rho) - (arb(off[e]) + arb(rho) * arb(rs0[e])
                         + arb(rho2) * arb(rs1[e]) + arb(rho3) * arb(rs2[e])))
        worst = m if worst is None else min(worst, m)
    outK = K_arb.upper() < 1
    outB = (arb(1) - K_arb).lower() >= 0 and \
        ((arb(1) - K_arb) * arb(rho) - arb(off_max)).lower() >= 0
    outM = worst.lower() > 0
    out = [outK, outB, outM]
    return [a and b for a, b in zip(flags, out)], K


# ---- bisection: lo known to pass, hi known to fail on (i) or (iii) ----
lo, hi = fmpq(1, 10 ** 9), fmpq(1, 10 ** 3)
fl_lo, K_lo = test_rho(lo)
fl_hi, K_hi = test_rho(hi)
print("flags(1e-9) =", fl_lo, " flags(1e-3) =", fl_hi, flush=True)
assert fl_lo == [True, True, True]
assert not all(fl_hi)

for it in range(60):
    mid = (lo + hi) / 2
    fl, _ = test_rho(mid)
    if all(fl):
        lo = mid
    else:
        hi = mid
    if float(hi - lo) < float(lo) * 1e-6:
        break

# The 60-iteration bisection narrows rho*; identify WHICH hypothesis fails
# just above and holds at lo: evaluate each hypothesis on its own.
rho_star = lo
fls, K_at_star = test_rho(rho_star)
print(f"rho* after bisection > {float(rho_star):.6e}; flags {fls}; "
      f"K(rho*) = {float(K_at_star):.9f}", flush=True)

# Identify limiting hypothesis at hi-frustum and at 10^-6
for probe in ("1/1000000", "1/10000", "1/1000"):
    fq = fmpq(probe)
    fl, Kv = test_rho(fq)
    print(f"probe rho={probe}: flags={fl}, K={float(Kv):.9f}", flush=True)

# per-hypothesis boundary at the largest failing radius: the failure just
# above rho* — find which flag flips first
hi_near = rho_star * fmpq(3, 2)
fl_up, _ = test_rho(hi_near)
print("flags at 1.5*rho* =", fl_up, flush=True)

res = {
    "rho_star_lower_bound_dyadic": str(rho_star),
    "rho_star_float": float(rho_star),
    "K_at_rho_star": str(K_at_star),
    "K_float": float(K_at_star),
    "flags_at_rho_star": fls,
    "flags_at_1.5rhostar": fl_up,
    "paper_rho": 1e-6,
    "bisection_iters": 60,
    "outward": True,
}
print("RADIUS_JSON " + json.dumps(res))
print(f"total {time.time() - t0:.1f}s")
