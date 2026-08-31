"""gate_a_verify.py — INDEPENDENT re-verification of the rank-25 Krawczyk
certificate of arXiv:2608.16649 / commit 816a01e, with EXACT rationals
(python-flint fmpq) and outward-rounded interval arithmetic (arb balls).

Everything is implemented from the PAPER's definitions (Prop 6 + section 3),
not by porting the upstream script: the T_O tensor is rebuilt independently
from the octonion multiplication table by Cayley-Dickson doubling, the
residual, Jacobian, W1/W2 enclosures, Krawczyk radii and the three
hypotheses are formed in exact rational arithmetic with fmpq/fmpz big-int
matrices, and every final comparison is re-evaluated as a strictly-outward
arb ball comparison.

Outputs a JSON verdict to stdout (and a human-readable log).
"""
import csv
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
from flint import fmpz, fmpz_mat, fmpq, fmpq_mat, arb, fmpq_poly

CERT = os.environ.get(
    "CERTDIR",
    os.path.join(os.path.dirname(__file__), "..", "scratch", "upstream_ref",
                 "certs", "rank25_cert"))

t_start = time.time()
log_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    log_lines.append(s)
    print(s, flush=True)


# ----------------------------------------------------------------------
# Independent octonion structure tensor via Cayley-Dickson doubling.
# Basis (1, i, j, k, l, il, jl, kl) i.e. H (+) H l with
# (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c));
# this is the standard normalization, entries in {-1,0,1}.
# ----------------------------------------------------------------------
def qmul(t, u):
    return (t[0] * u[0] - t[1] * u[1] - t[2] * u[2] - t[3] * u[3],
            t[0] * u[1] + t[1] * u[0] + t[2] * u[3] - t[3] * u[2],
            t[0] * u[2] - t[1] * u[3] + t[2] * u[0] + t[3] * u[1],
            t[0] * u[3] + t[1] * u[2] - t[2] * u[1] + t[3] * u[0])


def qconj(t):
    return (t[0], -t[1], -t[2], -t[3])


def t4(i, j):
    """basis-product table for H: e_i * e_j as index into (1,i,j,k)."""
    if i == 0:
        return (j, 1)
    if j == 0:
        return (i, 1)
    if i == j:
        return (0, -1)
    # i,j in {1,2,3}, i != j: e_i e_j = +/- e_k, k = 6 - i - j...
    # with basis (1,i,j,k): e1 e2 = e3; e2 e3 = e1; e3 e1 = e2; anti-symmetric.
    k = 6 - i - j
    sign = 1 if (i, j) in ((1, 2), (2, 3), (3, 1)) else -1
    return (k, sign)


def omul_pair(x, y):
    """octonion multiply given as pairs of quaternions."""
    a, b = x
    c, d = y
    left = [0] * 4
    right = [0] * 4
    qdc = qconj(d)
    qcc = qconj(c)
    for k in range(4):
        # left = a*c - qdc*b  (quaternion mult; sign-table scalar * e_idx)
        s1, i1 = t4(*[0, 0]) if False else (None, None)
    # product a*c via index table
    def qmul_idx(pq, uq, vq):
        """product of two quaternions given as (idx, sign) sparse pairs."""
        # expand: sum over index pairs
        return None
    # do quaternion arithmetic with explicit integer vectors
    av = [a[i] for i in range(4)] if isinstance(a, tuple) else a
    raise RuntimeError("unused path")


def to_vec(t):
    return list(t)


def h_mul(x, y):
    """8-dim integer octonion product via CD on integer vectors."""
    a, b = x[:4], x[4:]
    c, d = y[:4], y[4:]
    # left = a*c - conj(d)*b ; right = d*a + b*conj(c)
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

    dc = [d[0], -d[1], -d[2], -d[3]]
    cc = [c[0], -c[1], -c[2], -c[3]]
    left = [hv4(a, c)[k] - hv4(dc, b)[k] for k in range(4)]
    right = [hv4(d, a)[k] + hv4(b, cc)[k] for k in range(4)]
    return left + right


E8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
T = [[[h_mul(E8[i], E8[j])[k] for k in range(8)] for j in range(8)]
     for i in range(8)]
# consistency: integer entries, unit, norm multiplicativity at random pts
assert all(v in (-1, 0, 1) for sl in T for row in sl for v in row)
assert all(T[0][j][k] == (1 if j == k else 0) for j in range(8)
           for k in range(8))
rng = np.random.default_rng(12345)
for _ in range(20):
    x = rng.integers(-3, 4, 8)
    y = rng.integers(-3, 4, 8)
    xy = [sum(T[xi][yi_idx][k] * int(x[xi]) * int(y[yi_idx])
              for xi in range(8) for yi_idx in range(8)) for k in range(8)]
    nx, ny, nxy = sum(int(v) ** 2 for v in x), sum(int(v) ** 2 for v in y), \
        sum(v * v for v in xy)
    assert nx * ny == nxy
log("[A1-T1] independent Cayley-Dickson T_O built; unit + norm-mult OK")

# Cross-check my CD tensor against the upstream table (must be identical:
# same standard basis ordering (1,i,j,k,l,il,jl,kl)).
sys.path.insert(0, os.path.join(os.environ.get("UPSTREAM", ""), "verify")
                if os.environ.get("UPSTREAM") else
                os.path.join(os.path.dirname(__file__), "..", "scratch",
                             "upstream_ref", "verify"))
from octonion_core import octonion_tensor as upstream_octonion_tensor
Tup = upstream_octonion_tensor().tolist()
assert Tup == T, "my CD tensor differs from upstream table"
log("[A1-T2] my T_O == upstream octonion_tensor() entrywise")

# read certificate as EXACT dyadic rationals
def load_dyadic(fn):
    rows = []
    with open(fn) as f:
        for rec in csv.reader(f):
            rows.append([fmpq(frac_str) for frac_str in rec])
    return rows


def frac_from_float_str(s):
    """parse csv float text exactly as the dyadic float64 value."""
    fr = Fraction(float(s))
    return fmpq(f"{fr.numerator}/{fr.denominator}")


def load_dyadic_floats(fn):
    rows = []
    with open(fn) as f:
        for rec in csv.reader(f):
            rows.append([frac_from_float_str(v) for v in rec])
    return rows


A_rows = load_dyadic_floats(os.path.join(CERT, "A.csv"))
B_rows = load_dyadic_floats(os.path.join(CERT, "B.csv"))
C_rows = load_dyadic_floats(os.path.join(CERT, "C.csv"))
Y_rows = load_dyadic_floats(os.path.join(CERT, "Y_preconditioner.csv"))
S_idx = []
with open(os.path.join(CERT, "S_indices.csv")) as f:
    for rec in csv.reader(f):
        S_idx.extend(int(float(v)) for v in rec)
p, n = 8, 8
r = 25   # rank of the certified decomposition; A,B,C are 8x25
nb, nc = n, n
M = p * nb * nc
Nv = (p + nb + nc) * r
offB, offC = p * r, (p + nb) * r
S = list(S_idx)   # column order of J[:, S] exactly as archived
assert len(S) == M == 512 and len(set(S)) == M, (len(S), M)
log(f"[A2-T0] cert loaded: A,B,C {len(A_rows)}x{len(A_rows[0])}, "
    f"Y {len(Y_rows)}x{len(Y_rows[0])}, |S|={len(S)}, r={r}")

# build flat x = (A,B,C) as fmpq
x = [None] * Nv
for i in range(p):
    for s in range(r):
        x[i * r + s] = A_rows[i][s]
for j in range(nb):
    for s in range(r):
        x[offB + j * r + s] = B_rows[j][s]
for k in range(nc):
    for s in range(r):
        x[offC + k * r + s] = C_rows[k][s]

# exact residual g0 and Jacobian J (M x Nv) in fmpq
g0 = [fmpq(0)] * M
J = [[fmpq(0)] * Nv for _ in range(M)]
for i in range(p):
    for j in range(nb):
        for k in range(nc):
            e = (i * nb + j) * nc + k
            g0[e] = sum(A_rows[i][s] * B_rows[j][s] * C_rows[k][s]
                        for s in range(r)) - fmpq(T[i][j][k])
            for s in range(r):
                J[e][i * r + s] += B_rows[j][s] * C_rows[k][s]
                J[e][offB + j * r + s] += A_rows[i][s] * C_rows[k][s]
                J[e][offC + k * r + s] += A_rows[i][s] * B_rows[j][s]
r = len(A_rows[0])   # 25 rank-one terms; A is 8 x r
log(f"[A2-T0] cert loaded: A {len(A_rows)}x{r}, "
    f"Y {len(Y_rows)}x{len(Y_rows[0])}, |S|={len(S)}, r={r}")
# fmpq_mat versions for exact products
S = list(S_idx)   # keep upstream order; columns of J[:, S] in this order
Jsel_mat = fmpq_mat(M, M, [J[e][c] for e in range(M) for c in S])
Y_mat = fmpq_mat(M, M, [v for row in Y_rows for v in row])
YJ = Y_mat * Jsel_mat
Imat = fmpq_mat(M, M, [fmpq(1) if i == j else fmpq(0)
                       for i in range(M) for j in range(M)])
E = YJ - Imat
# A2: consistency of stored Y, exact
Emax = max(abs(E[i, j]) for i in range(M) for j in range(M))
log(f"[A2] ||Y*J_S - I||_inf exact = {Emax} (= {float(Emax):.3e})")
assert Emax < fmpq(1, 10 ** 6), "stored Y inconsistent (exact)"

# |Y| as fmpq_mat; g0 as fmpq_mat column; Y*g0 exact
def mat_abs(Mm):
    return fmpq_mat(Mm.nrows(), Mm.ncols(),
                    [abs(v) for v in Mm.entries()])


absY = mat_abs(Y_mat)
g0m = fmpq_mat(M, 1, g0)
Yg0 = Y_mat * g0m
off = [abs(Yg0[e, 0]) for e in range(M)]
off_max = max(off)
log(f"[A4-prep] off_max exact = {off_max} (= {float(off_max):.3e})")

absImYJ = mat_abs(Imat - YJ)
rs0 = [sum(absImYJ[e, c] for c in range(M)) for e in range(M)]

# W1, W2 restricted to S: exact fmpq matrices M x M (columns indexed S)
freeS = [False] * Nv
for c in S:
    freeS[c] = True
absx = [abs(v) for v in x]
W1s = [[fmpq(0)] * M for _ in range(M)]
W2s = [[0] * M for _ in range(M)]  # integer counts; exact (0/1 entries summed)
Spos = {c: t for t, c in enumerate(S)}
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

# A3/A4/A5 exact at rho = 1/10^6
rho = fmpq(1, 10 ** 6)
K = max(rs0[e] + rho * rs1[e] + rho ** 2 * rs2[e] for e in range(M))
K_lt_1 = K < fmpq(1)
margin_banach = (fmpq(1) - K) * rho - off_max
banach_ok = margin_banach >= 0
worst_margin = min(rho - (off[e] + rho * rs0[e] + rho ** 2 * rs1[e]
                          + rho ** 3 * rs2[e]) for e in range(M))
margins_ok = worst_margin > 0
log(f"[A3] exact K = {K} = {float(K):.9f};  K<1: {K_lt_1}")
log(f"[A4] exact (1-K)*rho - off_max = {margin_banach} "
    f"(= {float(margin_banach):.3e}); hold: {banach_ok}")
log(f"[A5] exact worst per-equation margin = {worst_margin} "
    f"(= {float(worst_margin):.3e}); hold: {margins_ok}")

# ---- outward-rounded arb re-evaluation of the three FINAL comparisons ----
# every input that is an exact rational goes in as an arb endpoint-ball;
# products/sums are carried at sufficient precision and rounded OUTWARD by
# comparing with explicit slack: we recompute K and margins as arb and
# require K.get_ub() < 1 (strict) etc.

def frac_to_arb_ends(qv):
    a = arb(qv)
    return a


K_arb = arb(0)
for e in range(M):
    row = rs0[e] + rho * rs1[e] + rho ** 2 * rs2[e]
    K_arb = max(K_arb, arb(row))
K_lt_1_out = K_arb.upper() < 1
banach_out = ((arb(1) - K_arb).lower() >= 0) and \
    (((arb(1) - K_arb) * arb(rho) - arb(off_max)).lower() >= 0)
worst_arb = None
for e in range(M):
    m = arb(rho) - (arb(off[e]) + arb(rho) * arb(rs0[e])
                    + arb(rho ** 2) * arb(rs1[e])
                    + arb(rho ** 3) * arb(rs2[e]))
    worst_arb = m if worst_arb is None else min(worst_arb, m)
margins_out = worst_arb.lower() > 0
log(f"[OUT-RND] arb K upper < 1: {K_lt_1_out} (K_arb={str(K_arb)[:40]})")
log(f"[OUT-RND] arb banach lb >= 0: {banach_out}")
log(f"[OUT-RND] arb worst margin lb > 0: {margins_out} "
    f"({str(worst_arb)[:40]})")

verdict = {
    "gate": "A",
    "A1_tensor_and_table": True,
    "A2_Y_consistency_exact": bool(Emax < fmpq(1, 10 ** 6)),
    "A2_Y_inf_norm": str(Emax),
    "A3_K_exact": str(K), "A3_K_float": float(K), "A3_K_lt_1": bool(K_lt_1),
    "A4_margin_banach": str(margin_banach), "A4_ok": bool(banach_ok),
    "A5_worst_margin": str(worst_margin), "A5_ok": bool(margins_ok),
    "outward_K_lt_1": bool(K_lt_1_out),
    "outward_banach": bool(banach_out),
    "outward_margins": bool(margins_out),
    "rho": "1/1000000",
    "M": M, "r": r, "Nv": Nv,
}
print("GATE_A_JSON " + json.dumps(verdict))
ok = all([verdict["A1_tensor_and_table"], verdict["A2_Y_consistency_exact"],
          verdict["A3_K_lt_1"], verdict["A4_ok"], verdict["A5_ok"],
          verdict["outward_K_lt_1"], verdict["outward_banach"],
          verdict["outward_margins"]])
log("GATE_A_VERDICT:", "PASS" if ok else "FAIL",
    f"({time.time() - t_start:.1f}s)")
sys.exit(0 if ok else 1)
