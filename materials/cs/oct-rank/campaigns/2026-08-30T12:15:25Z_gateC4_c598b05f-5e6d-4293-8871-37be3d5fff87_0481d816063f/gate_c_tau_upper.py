"""gate_c_tau_upper.py — INDEPENDENT exact re-verification that
R_R(tau_H) <= 7 via the upstream Krawczyk certificate, reimplemented in
python-flint (fmpq) rather than sympy, with outward-rounded arb
re-evaluation of the three final inequalities.

tau_H is the 3x4x4 tensor: slices are H left-multiplication by 1, i, j
(basis (1,i,j,k) of H). Certificate: certs/tower_cert/tau_r7 (dyadic
float64 CSVs read as exact dyadic rationals).

Verdict GATE_C_TAU_UPPER_JSON with all exact comparisons.
"""
import csv
import json
import os
import sys
import time
from fractions import Fraction

from flint import arb, fmpq, fmpq_mat

t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "scratch", "upstream_ref", "certs",
                    "tower_cert", "tau_r7")
log = lambda *a: print(*a, flush=True)


def qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])


E4 = [[1 if k == m else 0 for k in range(4)] for m in range(4)]
w = (E4[0], E4[1], E4[2])
tauT = [[[qmul(w[p], E4[b])[c] for c in range(4)] for b in range(4)]
        for p in range(3)]
assert all(v in (-1, 0, 1) for sl in tauT for row in sl for v in row)
log("[C-T1] tau_H rebuilt independently (3x4x4, entries in {-1,0,1})")


def frac_from_float_str(s):
    fr = Fraction(float(s))
    return fmpq(f"{fr.numerator}/{fr.denominator}")


def load(fn):
    with open(fn) as f:
        return [[frac_from_float_str(v) for v in rec]
                for rec in csv.reader(f)]


A_rows = load(os.path.join(CERT, "A.csv"))
B_rows = load(os.path.join(CERT, "B.csv"))
C_rows = load(os.path.join(CERT, "C.csv"))
Y_rows = load(os.path.join(CERT, "Y_preconditioner.csv"))
S_idx = []
with open(os.path.join(CERT, "S_indices.csv")) as f:
    for rec in csv.reader(f):
        S_idx.extend(int(float(v)) for v in rec)
S = list(S_idx)
p, nb, nc = 3, 4, 4
r = len(A_rows[0])
assert (len(A_rows), len(B_rows), len(C_rows)) == (p, nb, nc) == (3, 4, 4)
M = p * nb * nc
Nv = (p + nb + nc) * r
offB, offC = p * r, (p + nb) * r
assert len(S) == M and len(set(S)) == M
log(f"[C-T0] cert loaded r={r}, |S|={len(S)}")
assert r == 7

# exact residual and Jacobian
g0 = [fmpq(0)] * M
J = [[fmpq(0)] * Nv for _ in range(M)]
T = [[[fmpq(tauT[i][j][k]) for k in range(nc)] for j in range(nb)]
     for i in range(p)]
for i in range(p):
    for j in range(nb):
        for k in range(nc):
            e = (i * nb + j) * nc + k
            g0[e] = sum(A_rows[i][s] * B_rows[j][s] * C_rows[k][s]
                        for s in range(r)) - T[i][j][k]
            for s in range(r):
                J[e][i * r + s] += B_rows[j][s] * C_rows[k][s]
                J[e][offB + j * r + s] += A_rows[i][s] * C_rows[k][s]
                J[e][offC + k * r + s] += A_rows[i][s] * B_rows[j][s]

Jsel = fmpq_mat(M, M, [J[e][c] for e in range(M) for c in S])
Y = fmpq_mat(M, M, [v for row in Y_rows for v in row])
E = Y * Jsel - fmpq_mat(M, M, [fmpq(1 if i == j else 0)
                               for i in range(M) for j in range(M)])
Emax = max(abs(E[i, j]) for i in range(M) for j in range(M))
log(f"[C-T2] ||Y*J_S - I||_inf exact = {Emax} = {float(Emax):.3e}")
assert Emax < fmpq(1, 10**6)


def mat_abs(Mm):
    return fmpq_mat(Mm.nrows(), Mm.ncols(), [abs(v) for v in Mm.entries()])


g0m = fmpq_mat(M, 1, g0)
Yg0 = Y * g0m
off = [abs(Yg0[e, 0]) for e in range(M)]
off_max = max(off)

absY = mat_abs(Y)
ImYJ = fmpq_mat(M, M, [(fmpq(1) if i == j else fmpq(0)) - (Y * Jsel)[i, j]
                       for i in range(M) for j in range(M)])
absImYJ = mat_abs(ImYJ)
rs0 = [sum(absImYJ[e, c] for c in range(M)) for e in range(M)]

freeS = [False] * Nv
for c in S:
    freeS[c] = True
Spos = {c: t for t, c in enumerate(S)}
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
W1m = fmpq_mat(M, M, [v for row in W1s for v in row])
W2m = fmpq_mat(M, M, [fmpq(v) for row in W2s for v in row])
rs1 = [sum((absY * W1m)[e, c] for c in range(M)) for e in range(M)]
rs2 = [sum((absY * W2m)[e, c] for c in range(M)) for e in range(M)]

rho = fmpq(1, 10**5)  # RAD = 1e-5 (matches tau_r7 cert box)
K = max(rs0[e] + rho * rs1[e] + rho**2 * rs2[e] for e in range(M))
K_lt_1 = K < fmpq(1)
margin_banach = (fmpq(1) - K) * rho - off_max
banach_ok = margin_banach >= 0
worst = min(rho - (off[e] + rho * rs0[e] + rho**2 * rs1[e]
                   + rho**3 * rs2[e]) for e in range(M))
margins_ok = worst > 0
log(f"[C-T3] exact K = {float(K):.9f}  K<1: {K_lt_1}")
log(f"[C-T4] banach margin = {float(margin_banach):.3e}  ok: {banach_ok}")
log(f"[C-T5] worst margin = {float(worst):.3e}  ok: {margins_ok}")

# outward-rounded arb re-evaluation of the three final comparisons
K_arb = max(arb(rs0[e] + rho * rs1[e] + rho**2 * rs2[e]) for e in range(M))
out_K = bool(K_arb.upper() < 1)
out_banach = bool(((arb(1) - K_arb) * arb(rho) - arb(off_max)).lower() >= 0)
worst_arb = min(arb(rho) - (arb(off[e]) + arb(rho) * arb(rs0[e])
                            + arb(rho**2) * arb(rs1[e])
                            + arb(rho**3) * arb(rs2[e])) for e in range(M))
out_margins = bool(worst_arb.lower() > 0)
log(f"[C-OUT] K_lb<1: {out_K} | banach_lb>=0: {out_banach} | "
    f"margins_lb>0: {out_margins}")

verdict = {
    "gate": "C.tau_upper",
    "tensor": "tau_H 3x4x4 (1,i,j left-multiplication slices)",
    "claim": "R_R(tau_H) <= 7 Krawczyk-certified (upstream cert re-checked)",
    "certificate": os.path.relpath(CERT, HERE),
    "r": r, "rho": "1/100000",
    "Y_inf_norm": str(Emax), "Y_consistency": bool(Emax < fmpq(1, 10**6)),
    "K_exact": str(K), "K_float": float(K), "K_lt_1": bool(K_lt_1),
    "banach_margin": str(margin_banach), "banach_ok": bool(banach_ok),
    "worst_margin": str(worst), "margins_ok": bool(margins_ok),
    "outward_rounded": {"K_lt_1": out_K, "banach": out_banach,
                        "margins": out_margins},
    "elapsed_s": time.time() - t0,
}
print("GATE_C_TAU_UPPER_JSON " + json.dumps(verdict))
ok = all([verdict["Y_consistency"], K_lt_1, banach_ok, margins_ok,
          out_K, out_banach, out_margins])
log("GATE_C_TAU_UPPER_VERDICT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
