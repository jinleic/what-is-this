"""rf_krawczyk.py — Route F certification instrument (pre-registered).

Implements the square-slice Krawczyk/interval-Newton existence certificate of
campaigns/2026-09-01T04:30:00Z_routeF_kraw (pre_statement.md commit 739e63f,
REFINEMENT_LOG.md commit <this-freeze>): exact fmpq arithmetic end to end,
outward-arb final comparisons, explicit midpoint/radius control.

Systems:
  MAIN  : T_F = (L_1, L_i, L_j) octonion triple, target rank 13, 192 eqs,
          certified from seeds S_CP / S_EXT (polish pipeline P1..P5).
  C-TAU7: tau quaternion triple, r=7, 64 eqs, frozen tau_r7 certificate
          factors as seed (no polish; containment expected).
  C-TAU6: tau at r=6, 64 eqs, truncation seed + polish (containment MUST be
          absent; certified exclusion expected at small rungs).
  P-POS : synthetic rank-13 control (S-slice perturbation only).
  P-WRNG: corrupted table T_cor (entry [0,0,0] = 2), main seed.

Every rung of the fixed ladder is evaluated with: containment
  |c0_c| + R2_c(rho) < rho for all c (strict, fmpq, arb-outward re-checked)
  exclusion:   |c0_c| > rho + R2_c(rho)  for some c   (certified no-root)
Outputs per-system verdict rows to stdout and krawczyk_results.json;
deterministic; exits 0 iff the behavior contract holds.
"""
import os

_BLAS1 = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
          "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
          "NUMEXPR_NUM_THREADS": "1"}
for _k, _v in _BLAS1.items():
    os.environ.setdefault(_k, _v)

import csv
import json
import sys
import time
from fractions import Fraction

import numpy as np
from flint import arb, fmpq, fmpq_mat
from scipy import linalg as sla
from scipy.optimize import least_squares

T0 = time.process_time()
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))


LADDER = [Fraction(1, 1), Fraction(1, 10), Fraction(1, 100),
          Fraction(1, 10**3), Fraction(1, 10**4), Fraction(1, 10**5),
          Fraction(1, 10**6), Fraction(1, 10**7), Fraction(1, 10**8),
          Fraction(3, 10**9), Fraction(1, 10**9), Fraction(3, 10**10),
          Fraction(1, 10**10), Fraction(1, 10**11), Fraction(1, 10**12)]

def log(*a):
    print(*a, flush=True)

def dy(q_str):
    """Exact dyadic rational from a float or decimal string."""
    if isinstance(q_str, (int,)):
        return fmpq(q_str)
    fr = Fraction(float(q_str)) if not isinstance(q_str, str) else Fraction(q_str)
    return fmpq(f"{fr.numerator}/{fr.denominator}")

def dy_from_float(v):
    fr = Fraction(float(v))
    return fmpq(f"{fr.numerator}/{fr.denominator}")

def f2s(q):
    return str(q)

# ---------------------------------------------------------------- octonions
def cd_mul_list(x, y):
    """S3/RouteAF Cayley-Dickson recursion on lists."""
    n = len(x)
    if n == 1:
        return [x[0] * y[0]]
    m = n // 2
    a, b = x[:m], x[m:]
    c, d = y[:m], y[m:]

    def cj(z):
        if len(z) == 1:
            return z[:]
        m2 = len(z) // 2
        return cj(z[:m2]) + [-v for v in z[m2:]]

    left = [p - r for p, r in zip(cd_mul_list(a, c), cd_mul_list(cj(d), b))]
    right = [p + r for p, r in zip(cd_mul_list(d, a), cd_mul_list(b, cj(c)))]
    return left + right

def gateA_table():
    """gate-A index-table construction (t4/h_mul)."""
    def t4(i, j):
        if i == 0:
            return (j, 1)
        if j == 0:
            return (i, 1)
        if i == j:
            return (0, -1)
        k = 6 - i - j
        sign = 1 if (i, j) in ((1, 2), (2, 3), (3, 1)) else -1
        return (k, sign)

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
    return [[h_mul(E8[i], E8[j]) for j in range(8)] for i in range(8)]

def upstream_table():
    sys.path.insert(0, os.path.join(REPO, "scratch", "upstream_ref", "verify"))
    import octonion_core  # noqa
    return octonion_core.octonion_tensor().tolist()

def anchor_table():
    T_cd = [[None] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            e = [0] * 8
            e[i] = 1
            f = [0] * 8
            f[j] = 1
            T_cd[i][j] = cd_mul_list(e, f)
    T_ga = gateA_table()
    T_up = upstream_table()
    ok = (T_cd == T_ga == T_up)
    log("[ANCHOR-TBL] cd==gateA:", T_cd == T_ga, "| cd==upstream:",
        T_cd == T_up, "| ALL:", ok)
    # unit + norm multiplicativity at fixed points
    rng = np.random.default_rng(777)
    oknorm = True
    for _ in range(20):
        xv = [int(v) for v in rng.integers(-3, 4, 8)]
        yv = [int(v) for v in rng.integers(-3, 4, 8)]
        p = [0] * 8
        for i in range(8):
            for j in range(8):
                if xv[i] and yv[j]:
                    pr = T_cd[i][j]
                    for k in range(8):
                        p[k] += xv[i] * yv[j] * pr[k]
        nx = sum(v * v for v in xv)
        ny = sum(v * v for v in yv)
        nxy = sum(v * v for v in p)
        oknorm &= (nx * ny == nxy)
    ok_unit = all(T_cd[0][j] == [1 if k == j else 0 for k in range(8)]
                  and T_cd[j][0] == [1 if k == j else 0 for k in range(8)]
                  for j in range(8))
    log("[ANCHOR-TBL] e0 unit:", ok_unit, "| norm mult 20pts:", oknorm)
    return T_cd, (ok and oknorm and ok_unit)

def target_tensor(T_tbl, pmax=2):
    """3x8x8 tensor: T[p][b][c] = coeff of e_c in e_p * e_b."""
    return [[[T_tbl[p][b][c] for c in range(8)] for b in range(8)]
            for p in range(pmax + 1)]

# -------------------------------------------------- CP system machinery
class Sys:
    """Square-slice CP system: exact data + float polish utilities."""

    def __init__(self, Tint, r, seed_factors, name):
        self.name = name
        self.p = len(Tint)
        self.b = len(Tint[0])
        self.c = len(Tint[0][0])
        self.r = r
        self.pe, self.be, self.ce = self.p, self.b, self.c
        self.M = self.p * self.b * self.c
        self.NV = (self.p + self.b + self.c) * r
        self.Tint = Tint
        self.nsolved = self.M
        # seed_factors: list of (a,b,c) float arrays (each of the polisher's
        # conventions). We keep ONE canonical layout:
        #   z = [a0,.., a_{r-1}; b0,..; c0,..] with a 3xr etc.
        self.seed = seed_factors

    # ---- exact system evaluation
    # a, b, c are FLAT lists of fmpq in the canonical layout (a: p*r, then
    # b: b*r, then c: c*r); helper accessors index them.
    def resid_exact(self, a, b, c):
        M = self.M
        R = self.r
        g0 = [None] * M
        e = 0
        for pp in range(self.p):
            for bb in range(self.b):
                for cc in range(self.c):
                    s = fmpq(0)
                    for rr in range(R):
                        s += a[pp * R + rr] * b[bb * R + rr] * c[cc * R + rr]
                    g0[e] = s - fmpq(self.Tint[pp][bb][cc])
                    e += 1
        return g0

    def jac_exact(self, a, b, c):
        M, NV = self.M, self.NV
        R = self.r
        J = [[fmpq(0)] * NV for _ in range(M)]
        e = 0
        for pp in range(self.p):
            for bb in range(self.b):
                for cc in range(self.c):
                    for rr in range(R):
                        J[e][pp * R + rr] += b[bb * R + rr] * c[cc * R + rr]
                        J[e][(self.p + bb) * R + rr] += a[pp * R + rr] \
                            * c[cc * R + rr]
                        J[e][(self.p + self.b + cc) * R + rr] \
                            += a[pp * R + rr] * b[bb * R + rr]
                    e += 1
        return J

    # ---- float side (polish)
    def split(self, z):
        r = self.r
        p, b, c = self.p, self.b, self.c
        return (z[:p * r].reshape(p, r), z[p * r:(p + b) * r].reshape(b, r),
                z[(p + b) * r:(p + b + c) * r].reshape(c, r))

    def resid_f(self, z):
        a, b, c = self.split(z)
        T = np.asarray(self.Tint, dtype=float)
        return (np.einsum("pr,br,cr->pbc", a, b, c) - T).ravel()

    def jac_f(self, z):
        a, b, c = self.split(z)
        p, b_, c_, r = self.p, self.b, self.c, self.r
        M, NV = self.M, self.NV
        J = np.zeros((M, NV))
        e = 0
        for pp in range(p):
            for bb in range(b_):
                for cc in range(c_):
                    for rr in range(r):
                        J[e, pp * r + rr] += b[bb, rr] * c[cc, rr]
                        J[e, (p + bb) * r + rr] += a[pp, rr] * c[cc, rr]
                        J[e, (p + b_ + cc) * r + rr] += a[pp, rr] * b[bb, rr]
                    e += 1
        return J

    def balance(self, z):
        """Gauge move: geometric-mean equalization per rank term (float)."""
        a, b, c = self.split(z.copy())
        for s in range(self.r):
            na, nb, nc = (np.linalg.norm(a[:, s]), np.linalg.norm(b[:, s]),
                          np.linalg.norm(c[:, s]))
            if min(na, nb, nc) == 0:
                continue
            g = (na * nb * nc) ** (1.0 / 3.0)
            a[:, s] *= g / na
            b[:, s] *= g / nb
            c[:, s] *= g / nc
        return np.concatenate([a.ravel(), b.ravel(), c.ravel()])

    def rel(self, z):
        return float(np.linalg.norm(self.resid_f(z)) /
                     np.linalg.norm(np.asarray(self.Tint, dtype=float)))

# -------------------------------------------------- certification core
def certify(sysobj, z_dy, slice_rule, tag):
    """Exact square-slice Krawczyk over the fixed ladder.

    z_dy: list of exact fmpq seed coordinates (length NV), canonical layout.
    slice_rule: 'qrcp' (QRCP pivot order) or natural fallback; selection
    diagnostics float-only; the classification itself is recorded.
    Returns dict with per-rung data and verdict.
    """
    p, b, c, r = sysobj.p, sysobj.b, sysobj.c, sysobj.r
    NV, M = sysobj.NV, sysobj.M
    # flat exact vectors (canonical layout) for the exact evaluators
    a0f = z_dy[:p * r]
    b0f = z_dy[p * r:(p + b) * r]
    c0f = z_dy[(p + b) * r:(p + b + c) * r]
    g0 = sysobj.resid_exact(a0f, b0f, c0f)
    Jfull = sysobj.jac_exact(a0f, b0f, c0f)
    # 2D views for the margin algebra
    a0 = [[z_dy[pp * r + rr] for rr in range(r)] for pp in range(p)]
    b0 = [[z_dy[(p + bb) * r + rr] for rr in range(r)] for bb in range(b)]
    c0c = [[z_dy[(p + b + cc) * r + rr] for rr in range(r)]
           for cc in range(c)]

    # ---- slice choice
    # fmpq -> float64 via the string form (exact rational, deliberately
    # rounded to the nearest float for the FLOAT diagnostics only)
    Jf = sysobj.jac_f(np.asarray([float(v) for v in z_dy]))
    sv = np.linalg.svd(Jf, compute_uv=False)
    rankJ = int((sv > 1e-9).sum())
    pivot_note = {}
    if slice_rule == "qrcp" and rankJ == M:
        _, Rq, piv = sla.qr(Jf, mode="economic", pivoting=True)
        diagabs = np.abs(np.diag(Rq))
        if diagabs.min() >= 1e-9:
            S_cols = piv[:M].tolist()
            F_cols = sorted(piv[M:].tolist())
            rule_used = "qrcp"
            pivot_note = {"min_pivot": float(diagabs.min())}
        else:
            S_cols = list(range(M))
            F_cols = list(range(M, NV))
            rule_used = "natural-fallback-minpivot"
    else:
        S_cols = list(range(M))
        F_cols = list(range(M, NV))
        rule_used = "natural-fallback-rank" if slice_rule == "qrcp" \
            else slice_rule
    log(f"[{tag}] rank(J)@1e-9={rankJ}/{NV} rule={rule_used} {pivot_note}")

    Sset = set(S_cols)
    J_S = fmpq_mat(M, M, [Jfull[e][col] for e in range(M) for col in S_cols])
    tinv0 = time.process_time()
    Y = J_S.inv()
    tinv = time.process_time() - tinv0
    # exact check Y * J_S == I
    I_ok = Y * J_S == fmpq_mat(M, M, [fmpq(1) if i == j else fmpq(0)
                                      for i in range(M) for j in range(M)])
    log(f"[{tag}] exact inverse ok={I_ok} ({tinv:.2f}s cpu)")

    g0m = fmpq_mat(M, 1, g0)
    Yg0 = Y * g0m
    cbou = [Yg0[e, 0] for e in range(M)]
    absY = fmpq_mat(M, M, [abs(v) for v in Y.entries()])

    # seed factor magnitudes for W_e
    absa = [[abs(v) for v in row] for row in a0]
    absb = [[abs(v) for v in row] for row in b0]
    absc = [[abs(v) for v in row] for row in c0c]

    # (sv_kind helper folded into W_e below; see REFINEMENT_LOG R1.)

    # J_F transversality diagnostic (exact rank via rref on rational block)
    J_F_cols = F_cols
    JF = fmpq_mat(M, len(J_F_cols), [Jfull[e][col] for e in range(M)
                                     for col in J_F_cols])
    rankJF = JF.rref()[1] if len(J_F_cols) else 0
    log(f"[{tag}] exact rank(J_F)={rankJF}/{len(J_F_cols)}")

    rungs = []
    verdict = "NO-CONTAINMENT-ANY-RUNG"
    excl_any = False
    excl_rung = None
    pass_rung = None
    rho_q = {fr: fmpq(f"{fr.numerator}/{fr.denominator}") for fr in LADDER}
    for rho in LADDER:
        rho = rho_q[rho]          # exact fmpq from here on
        rho_f = float(rho)
        # W_e(rho)
        We = [fmpq(0)] * M
        e = 0
        for pp in range(p):
            for bb in range(b):
                for cc in range(c):
                    wr = fmpq(0)
                    for rr in range(r):
                        Pa = rho if (pp * r + rr) in Sset else fmpq(0)
                        Pb = rho if ((p + bb) * r + rr) in Sset else fmpq(0)
                        Pc = rho if ((p + b + cc) * r + rr) in Sset \
                            else fmpq(0)
                        # quadratic + cubic remainder terms
                        wr += Pa * (Pb * absb[bb][rr] + Pc * absc[cc][rr])
                        wr += Pb * Pc * absa[pp][rr] + (Pa * Pb * Pc
                                                        if (Pa != 0 and Pb != 0
                                                            and Pc != 0)
                                                        else fmpq(0))
                    We[e] = wr
                    e += 1
        # R2_c = sum_e |Y[c][e]| * We
        Wem = fmpq_mat(M, 1, We)
        absYW = absY * Wem
        R2 = [absYW[e, 0] for e in range(M)]
        # containment / exclusion per coordinate
        n_pass = 0
        n_excl = 0
        margins = []
        for cc in range(M):
            lhs = abs(cbou[cc]) + R2[cc]
            if lhs < rho:
                n_pass += 1
            if abs(cbou[cc]) > rho + R2[cc]:
                n_excl += 1
            margins.append((cc, str(abs(cbou[cc])), str(R2[cc])))
        row = {"tag": tag, "rho": f2s(rho), "rho_f": rho_f,
               "n_contain": n_pass, "n_excl": n_excl, "rule": rule_used,
               "rankJF": int(rankJF)}
        if n_pass == M:
            verdict = "CONTAINMENT"
            pass_rung = f2s(rho)
            # arb outward re-verification of the strictness for ALL coords
            ok_arb = True
            worst = None
            for cc in range(M):
                # every input enters arb as an EXACT fmpq midpoint (the
                # exact rational itself; no float endpoint subtraction).
                lhs_arb = arb(abs(cbou[cc])) + arb(R2[cc])
                gap_arb = arb(rho) - lhs_arb
                lb = gap_arb.lower()
                if not (lb > 0):
                    ok_arb = False
                if worst is None or lb < worst:
                    worst = lb
            row["arb_recheck_all_strict"] = ok_arb
            row["worst_gaparb_upper"] = str(worst) if worst is not None else None
            if not ok_arb:
                verdict = "ARB-DISAGREE"
        elif n_excl > 0:
            excl_any = True
            if excl_rung is None:
                excl_rung = f2s(rho)
            row["excl_coords"] = [cc for cc in range(M)
                                  if abs(cbou[cc]) > rho + R2[cc]][:8]
        rungs.append(row)
        log(f"[{tag}] rho={rho_f:.0e} contain={n_pass}/{M} excl={n_excl}/{M}")
        if verdict == "CONTAINMENT":
            break
    return {"tag": tag, "verdict": verdict, "pass_rung": pass_rung,
            "excl_any": excl_any, "excl_rung": excl_rung,
            "rule": rule_used, "rankJ": rankJ, "rankJF": int(rankJF),
            "inv_ok": bool(I_ok), "sv_min": float(sv[-1]),
            "sv_max": float(sv[0]), "rungs": rungs,
            "S_cols": S_cols, "F_cols": F_cols,
            "bounce_max": f2s(max(abs(v) for v in cbou)),
            "cpu_inv_s": tinv}

# -------------------------------------------------- pipeline driver
def polish_pipeline(sysobj, z0, rounds, max_nfev, tag):
    z = z0
    hist = []
    for k in range(rounds):
        res = least_squares(sysobj.resid_f, z, jac=sysobj.jac_f, method="trf",
                            tr_solver="exact", x_scale="jac",
                            max_nfev=max_nfev, ftol=3e-15, xtol=3e-15,
                            gtol=3e-15)
        z = sysobj.balance(res.x)
        rr = sysobj.rel(z)
        hist.append({"round": k + 1, "rel": rr,
                     "nfev": int(res.nfev), "status": int(res.status)})
        log(f"[{tag}] round {k+1}: rel={rr:.3e} nfev={res.nfev}")
        if rr < 1e-11:
            break
    return z, hist

def to_dy(z):
    return [dy_from_float(v) for v in z]

def main():
    Ttbl, anchor_ok = anchor_table()
    log(f"[ANCHOR] table 3-way byte-exact + unit + norm: {anchor_ok}")
    if not anchor_ok:
        log("[ANCHOR] FAIL — campaign HALTS per pre-statement sec 1")
        sys.exit(3)

    out = {"env": sys.version.split()[0], "ladder": [f2s(r) for r in LADDER]}
    results = {}
    # ---------------- MAIN system: the CONJUGATED Route-F tensor
    # tau boxtimes s = blockdiag(tau_p, tau_p) 3x8x8 (Route AF equivalence:
    # the frozen seeds were optimized against THIS target; its npz
    # 'convention' field states it). Anchor the equivalence entrywise:
    # T_sim = diag(1..1,-1) gives T_sim L_{e_n} T_sim == blockdiag(Lq_n,S2
    # Lq_n S2) for n=0..2 (s3_routeF3_global V1).
    # tau_4 = raw (L_1, L_i, L_j) 4x4 coordinate tensors on H
    tau4 = [[[Ttbl[pp][bb][cc] for cc in range(4)] for bb in range(4)]
            for pp in range(3)]
    # blockdiag doubled tensor: slice pp, block (bb,cc) >= 4 maps to
    # tau4[pp][bb-4][cc-4]; top-left block tau4[pp][bb][cc].
    TF = [[[ (tau4[pp][bb][cc] if (bb < 4 and cc < 4) else
              (tau4[pp][bb - 4][cc - 4] if (bb >= 4 and cc >= 4) else 0))
            for cc in range(8)] for bb in range(8)] for pp in range(3)]
    # ENTRYWISE rank-equivalence anchors (the corrected form of the frozen
    # similitude chain; each is an exact diagonal-similitude identity):
    #  (1) L-matrix of (L_1,L_i,L_j), L[i][j] = coeff of e_i in e_pp * e_j,
    #      IS ALREADY blockdiag(Lq_p, S2 Lq_p S2) with S2 = diag(1,1,1,-1);
    #  (2) Sbig = diag(1,1,1,1,1,1,1,-1) conjugates blockdiag(Lq_p, Lq_p)
    #      (the certification target) to blockdiag(Lq_p, S2 Lq_p S2).
    # Both checked entrywise on all 3 slices (fmpq-identical integers).
    # Consequence: rank((L_1,L_i,L_j)) = rank(blockdiag(tau,tau)) — the
    # certified box for the target certifies the Route-F triple.
    S2 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
    S2[3][3] = -1
    Sbig = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
    Sbig[7][7] = -1
    def matmul(A, B):
        n = len(A); m = len(B[0]); k = len(B)
        return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
                for i in range(n)]
    anchor_blk = True
    for pp in range(3):
        Lq4 = [[Ttbl[pp][j][i] for j in range(4)] for i in range(4)]
        L8 = [[Ttbl[pp][j][i] for j in range(8)] for i in range(8)]
        S2LqS2 = matmul(matmul(S2, Lq4), S2)
        blk1 = [[Lq4[i][j] if (i < 4 and j < 4) else
                 (S2LqS2[i - 4][j - 4] if (i >= 4 and j >= 4) else 0)
                 for j in range(8)] for i in range(8)]
        f1 = blk1 == L8
        blk0 = [[Lq4[i][j] if (i < 4 and j < 4) else
                 (Lq4[i - 4][j - 4] if (i >= 4 and j >= 4) else 0)
                 for j in range(8)] for i in range(8)]
        f2 = matmul(matmul(Sbig, blk0), Sbig) == blk1
        anchor_blk &= f1 and f2
        log(f"[ANCHOR-BLK] slice {pp}: L==blk(Lq,S2LqS2)={f1} Sbig={f2}")
    log(f"[ANCHOR-BLK] equivalence chain holds on 3 slices: {anchor_blk}")
    if not anchor_blk:
        log("[ANCHOR-BLK] FAIL — HALT (equivalence anchor)")
        sys.exit(3)
    print("TF_NONZERO", np.count_nonzero(np.asarray(TF)))

    # load the two frozen seeds
    d_cp = np.load(os.path.join(REPO, "campaigns",
                                "2026-08-31T08:02:18Z_routeAF",
                                "rank13_candidate.npz"))
    d_ex = np.load(os.path.join(REPO, "campaigns",
                                "2026-08-31T08:02:18Z_routeAF",
                                "rank13_extension_candidate.npz"))
    sysmain = Sys(TF, 13, None, "MAIN")
    # verify stored rel anchors
    rel_cp = sysmain.rel(d_cp["x"])
    rel_ex = sysmain.rel(d_ex["x"])
    log(f"[MAIN] seed rel: CP={rel_cp:.12e} (stored 1.3252060344282049e-3) "
        f"EXT={rel_ex:.12e} (stored 9.99874493968644e-5)")
    anchor_cp = abs(rel_cp - 1.3252060344282049e-3) < 1e-12
    anchor_ex = abs(rel_ex - 9.99874493968644e-5) < 1e-12
    log(f"[MAIN] seed rel anchors: {anchor_cp} {anchor_ex}")

    # polish both pipelines
    z_cp, h_cp = polish_pipeline(sysmain, d_cp["x"], 5, 5000, "MAIN-CP")
    z_ex, h_ex = polish_pipeline(sysmain, d_ex["x"], 5, 5000, "MAIN-EXT")
    rel_cp_f, rel_ex_f = sysmain.rel(z_cp), sysmain.rel(z_ex)
    log(f"[MAIN] polished rel: CP={rel_cp_f:.6e} EXT={rel_ex_f:.6e}")
    # tie -> the CP-direct pipeline (pre-statement 2.3)
    pick = "CP" if rel_cp_f <= rel_ex_f else "EXT"
    zsel = z_cp if pick == "CP" else z_ex
    hist_sel = h_cp if pick == "CP" else h_ex
    hist_oth = h_ex if pick == "CP" else h_cp
    log(f"[MAIN] selected pipeline: {pick} rel={sysmain.rel(zsel):.6e} "
        f"maxentry={np.abs(zsel).max():.3e}")

    # ---- freeze the selected candidate factors as exact decimal strings
    with open(os.path.join(HERE, "candidate_main_factors.csv"), "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["FACTOR"] + [repr(float(v)) for v in zsel])


    # certification attempt (qrcp rule)
    res_main = certify(sysmain, to_dy(zsel), "qrcp", "MAIN")
    res_main["seed_rel_history_sel"] = hist_sel
    res_main["seed_rel_history_oth"] = hist_oth
    res_main["selected"] = pick
    res_main["rel_final"] = sysmain.rel(zsel)
    results["MAIN"] = {k: v for k, v in res_main.items() if k != "S_cols"}

    # ---------------- C-TAU7: tau r=7, frozen certificate seed, no polish
    E4 = [[1 if k == m else 0 for k in range(4)] for m in range(4)]

    def qmul(a, b):
        return (a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
                a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
                a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
                a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0])

    tauT = [[[qmul(E4[p_], E4[bb])[cc] for cc in range(4)]
             for bb in range(4)] for p_ in range(3)]
    cert = os.path.join(REPO, "scratch", "upstream_ref", "certs",
                        "tower_cert", "tau_r7")

    def load(fn):
        with open(os.path.join(cert, fn)) as f:
            return [[dy(v) for v in rec] for rec in csv.reader(f)]

    A7, B7, C7 = load("A.csv"), load("B.csv"), load("C.csv")
    z_tau7 = ([A7[p_][rr] for p_ in range(3) for rr in range(7)]
              + [B7[bb][rr] for bb in range(4) for rr in range(7)]
              + [C7[cc][rr] for cc in range(4) for rr in range(7)])
    systau7 = Sys(tauT, 7, None, "TAU7")
    log("[TAU7] natural-slice attempt would be singular (rank(J)=48/77); using the registered qrcp rule (identical path)")
    res_tau7 = certify(systau7, z_tau7, "qrcp", "TAU7")

    with open(os.path.join(HERE, "candidate_tau7_factors.csv"), "w",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["FACTOR"] + [str(v) for v in z_tau7])
    results["TAU7"] = res_tau7

    # ---------------- C-TAU6: truncation seed + polish; MUST NOT contain
    z_tau6_f = np.concatenate([
        np.asarray(A7, dtype=object)[:, :6].astype(float).ravel(),
        np.asarray(B7, dtype=object)[:, :6].astype(float).ravel(),
        np.asarray(C7, dtype=object)[:, :6].astype(float).ravel()])
    systau6 = Sys(tauT, 6, None, "TAU6")
    z_tau6, h6 = polish_pipeline(systau6, z_tau6_f, 5, 5000, "TAU6")
    log(f"[TAU6] polished rel={systau6.rel(z_tau6):.6e}")
    res_tau6 = certify(systau6, to_dy(z_tau6), "qrcp", "TAU6")
    res_tau6["rel_final"] = systau6.rel(z_tau6)
    results["TAU6"] = {k: v for k, v in res_tau6.items() if k != "S_cols"}

    # ---------------- P-POS: synthetic rank-13 control
    # Root factors are exact dyadics (float64 values); the target tensor U
    # is the EXACT dyadic CP product of those factors, so (a*,b*,c*) solves
    # the P-pos system exactly in fmpq arithmetic. The plant perturbs ONLY
    # S-coordinates (REFINEMENT_LOG R2), by scale 1e-8 (R3).
    rngA = np.random.default_rng(20260901)
    astar = rngA.normal(scale=0.7, size=(3, 13))
    bstar = rngA.normal(scale=0.7, size=(8, 13))
    cstar = rngA.normal(scale=0.7, size=(8, 13))
    astr_f = astar.ravel()
    bstr_f = bstar.ravel()
    cstr_f = cstar.ravel()
    a_dy = [dy_from_float(v) for v in astr_f]
    b_dy = [dy_from_float(v) for v in bstr_f]
    c_dy = [dy_from_float(v) for v in cstr_f]
    U_dy = [[[sum(a_dy[pp * 13 + rr] * b_dy[bb * 13 + rr]
                  * c_dy[cc * 13 + rr] for rr in range(13))
              for cc in range(8)] for bb in range(8)] for pp in range(3)]
    syspos = Sys(U_dy, 13, None, "PPOS")
    z_rootf = np.concatenate([astr_f, bstr_f, cstr_f])
    # S-slice at the ROOT (the root itself is a rank-13 solution)
    Jf_root = syspos.jac_f(z_rootf)
    _, Rq_r, piv_r = sla.qr(Jf_root, mode="economic", pivoting=True)
    Scols = piv_r[:syspos.M].tolist()
    z_plant = z_rootf.copy()
    rngB = np.random.default_rng(20260902)
    pert = 1e-8 * rngB.normal(size=syspos.NV)
    for col in Scols:
        z_plant[col] += pert[col]
    rel_plant = syspos.rel(z_plant)
    log(f"[PPOS] plant rel={rel_plant:.6e} (S-only perturbation 1e-8)")
    res_pos = certify(syspos, to_dy(z_plant), "qrcp", "PPOS")
    results["PPOS"] = {k: v for k, v in res_pos.items() if k != "S_cols"}

    # ---------------- P-WRNG: corrupted table (entry [0,0,0] = 2)
    Tcor = [[[v for v in row] for row in mat] for mat in TF]
    Tcor[0][0][0] = 2
    syswrong = Sys(Tcor, 13, None, "PWRNG")
    res_wrong = certify(syswrong, to_dy(zsel), "qrcp", "PWRNG")
    results["PWRNG"] = {k: v for k, v in res_wrong.items() if k != "S_cols"}

    # ---------------- behavior contract
    ok = True
    if results["TAU7"]["verdict"] != "CONTAINMENT":
        log("[CONTRACT] TAU7 failed to certify — instrument broken?")
        ok = False
    if results["TAU6"]["verdict"] == "CONTAINMENT":
        log("[CONTRACT] TAU6 CONTAINED — FALSE CERTIFICATE — quarantine")
        ok = False
    if results["PPOS"]["verdict"] != "CONTAINMENT":
        log("[CONTRACT] PPOS failed to certify — instrument broken?")
        ok = False
    if res_main["verdict"] == "CONTAINMENT":
        log(f"[MAIN] CONTAINMENT at rung {res_main['pass_rung']} "
            f"arb={res_main['rungs'][-1].get('arb_recheck_all_strict')}")
    else:
        log("[MAIN] no containment at any rung -> FAILURE TO CERTIFY")

    out["results"] = {k: {kk: vv for kk, vv in v.items()
                          if kk not in ("rungs",)}
                      for k, v in results.items()}
    out["rungs_full"] = {k: v["rungs"] for k, v in results.items()}
    out["behavior_ok"] = ok
    out["cpu_s"] = time.process_time() - T0
    with open(os.path.join(HERE, "krawczyk_results.json"), "w") as f:
        json.dump(out, f, indent=1)
    log(f"[DONE] cpu={out['cpu_s']:.2f}s behavior_ok={ok}")
    sys.exit(0 if ok else 4)

if __name__ == "__main__":
    main()
