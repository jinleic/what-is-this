"""n1_run.py — N1 exact-CP-completion campaign driver.

Pre-registered: pre_statement.md (== cs/oct-rank/n1_prereg.md, commit
d6c7e44d6901eb025835313fe8bf9c7b01ee3a81, sha256 d569057e...f82).

Sequence (strict order, per prereg):
  1. CONTROLS (section 4) through the FROZEN instrument path
     (n1_instrument_base.py == frozen rf_krawczyk.py, sha256-pinned):
       C-TAU7 must CONTAIN; C-TAU6 must NOT contain (certified exclusion);
       P-POS must CONTAIN; P-WRNG must NOT contain at any rung.
     Any break -> INVALID INSTRUMENT, exit 5, main sweep not run.
  2. MAIN sweep: 40 seeds (32 random class R, 7 merge class M, 1 frozen
     class F), 3 TRF rounds x 5000 nfev each, identical solver settings to
     the frozen instrument; per-seed exact residual admission at
     ||g0||_max <= 1e-7 (exact fmpq).
  3. Admitted seeds -> exact Newton (max 3 iters, exact fmpq, fresh slice
     + inverse each iter). Exact-zero iterate (CP13(x)==TF entrywise) ->
     WITNESS: rank 13 EXACTLY [MACHINE-VERIFIED].
  4. Admitted non-exact candidates -> frozen 15-rung Krawczyk ladder through
     certify() (containment PASS also settles rank 13).
  5. Kill otherwise: FAILURE TO CERTIFY (no rank claim either direction).

All in-process accounting: time.process_time(). Global CPU cap 2 h.
"""
import os

_BLAS1 = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
          "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
          "NUMEXPR_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}
for _k, _v in _BLAS1.items():
    os.environ[_k] = _v

import csv
import json
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
ROUTEF = os.path.normpath(os.path.join(HERE, "..",
                        "2026-09-01T04:30:00Z_routeF_kraw"))
ROUTEAF = os.path.normpath(os.path.join(HERE, "..",
                           "2026-08-31T08:02:18Z_routeAF"))
CERTDIR = os.path.normpath(os.path.join(REPO, "scratch", "upstream_ref",
                          "certs", "tower_cert", "tau_r7"))
sys.path.insert(0, HERE)
import n1_instrument_base as base  # frozen certification core (no main())

import n1_instrument_base as rfk  # alias symmetry with frozen naming

T0 = time.process_time()
CPU_CAP = 2 * 3600.0
ADMIT = Fraction(1, 10**7)          # ||g0||_max admission (exact)
NEWTON_ITERS = 3
TRF_ROUNDS = 3
TRF_NFEV = 5000

def log(*a):
    print(*a, flush=True)

def cpu_left():
    return CPU_CAP - (time.process_time() - T0)

def dy_from_float(v):
    fr = Fraction(float(v))
    return rfk.fmpq(f"{fr.numerator}/{fr.denominator}")

def to_dy(z):
    return [dy_from_float(v) for v in z]

def gmax_exact(sysobj, z_dy):
    p, b, c, r = sysobj.p, sysobj.b, sysobj.c, sysobj.r
    g0 = sysobj.resid_exact(z_dy[:p * r], z_dy[p * r:(p + b) * r],
                            z_dy[(p + b) * r:(p + b + c) * r])
    return max((abs(v) for v in g0), default=rfk.fmpq(0)), g0

def exact_witness_check(sysobj, z_dy):
    """CP13(x) == T entrywise in fmpq?"""
    _, g0 = gmax_exact(sysobj, z_dy)
    return all(v == 0 for v in g0)

def exact_newton(sysobj, z_dy, iters=NEWTON_ITERS):
    """Exact fmpq Newton on the square slice; returns (z_dy, history).

    Slice rule and margin algebra identical to certify(): qrcp with natural
    fallback; Y = J_S^{-1} exact each iteration; x <- x - Y g0.
    """
    hist = []
    for it in range(iters):
        p, b, c, r = sysobj.p, sysobj.b, sysobj.c, sysobj.r
        a0f = z_dy[:p * r]
        b0f = z_dy[p * r:(p + b) * r]
        c0f = z_dy[(p + b) * r:(p + b + c) * r]
        g0 = sysobj.resid_exact(a0f, b0f, c0f)
        gm = max((abs(v) for v in g0), default=rfk.fmpq(0))
        hist.append({"iter": it, "gmax": str(gm)})
        log(f"[NEWTON] iter {it}: ||g0||_max = {gm}")
        if gm == 0:
            return z_dy, hist, True
        Jfull = sysobj.jac_exact(a0f, b0f, c0f)
        M, NV = sysobj.M, sysobj.NV
        # float diagnostics for the slice choice (frozen rule)
        Jf = sysobj.jac_f(np.asarray([float(v) for v in z_dy]))
        sv = np.linalg.svd(Jf, compute_uv=False)
        rankJ = int((sv > 1e-9).sum())
        if rankJ == M:
            from scipy import linalg as sla
            _, Rq, piv = sla.qr(Jf, mode="economic", pivoting=True)
            diagabs = np.abs(np.diag(Rq))
            if diagabs.min() >= 1e-9:
                S_cols = piv[:M].tolist()
                rule = "qrcp"
            else:
                S_cols = list(range(M)); rule = "natural"
        else:
            S_cols = list(range(M)); rule = "natural"
        J_S = rfk.fmpq_mat(M, M, [Jfull[e][col] for e in range(M)
                                  for col in S_cols])
        Y = J_S.inv()
        ok = Y * J_S == rfk.fmpq_mat(M, M, [rfk.fmpq(1) if i == j
                                            else rfk.fmpq(0)
                                            for i in range(M)
                                            for j in range(M)])
        if not ok:
            log("[NEWTON] exact inverse check FAILED"); break
        g0m = rfk.fmpq_mat(M, 1, g0)
        step = Y * g0m
        # x <- x - Y g0 on the S coordinates (F coordinates untouched)
        z_new = list(z_dy)
        for e, col in enumerate(S_cols):
            z_new[col] = z_new[col] - step[e, 0]
        z_dy = z_new
        p_, b_, c_, r_ = sysobj.p, sysobj.b, sysobj.c, sysobj.r
        g1 = sysobj.resid_exact(z_dy[:p_ * r_], z_dy[p_ * r_:(p_ + b_) * r_],
                                z_dy[(p_ + b_) * r_:(p_ + b_ + c_) * r_])
        gm1 = max((abs(v) for v in g1), default=rfk.fmpq(0))
        hist[-1]["gmax_after"] = str(gm1)
        log(f"[NEWTON] iter {it}: ||g1||_max = {gm1} rule={rule}")
        if gm1 == 0:
            return z_dy, hist, True
    return z_dy, hist, exact_witness_check(sysobj, z_dy)

def polish(sysobj, z0, tag):
    from scipy.optimize import least_squares
    z = np.asarray(z0, dtype=float)
    hist = []
    for k in range(TRF_ROUNDS):
        res = least_squares(sysobj.resid_f, z, jac=sysobj.jac_f,
                            method="trf", tr_solver="exact", x_scale="jac",
                            max_nfev=TRF_NFEV, ftol=3e-15, xtol=3e-15,
                            gtol=3e-15)
        z = sysobj.balance(res.x)
        rr = sysobj.rel(z)
        hist.append({"round": k + 1, "rel": rr, "nfev": int(res.nfev),
                     "status": int(res.status)})
        log(f"[{tag}] round {k+1}: rel={rr:.3e} nfev={res.nfev}")
    return z, hist

# -------------------------------------------------------------- targets
def build_targets():
    Ttbl, ok = rfk.anchor_table()
    if not ok:
        log("[ANCHOR] table FAIL — HALT"); sys.exit(3)
    tau4 = [[[Ttbl[pp][bb][cc] for cc in range(4)] for bb in range(4)]
            for pp in range(3)]
    TF = [[[ (tau4[pp][bb][cc] if (bb < 4 and cc < 4) else
              (tau4[pp][bb - 4][cc - 4] if (bb >= 4 and cc >= 4) else 0))
            for cc in range(8)] for bb in range(8)] for pp in range(3)]
    tauT = [[[rfk_qmul(rfk_E()[pp], rfk_E()[bb])[cc] for cc in range(4)]
             for bb in range(4)] for pp in range(3)]
    return Ttbl, TF, tauT

def rfk_E():
    return [[1 if k == m else 0 for k in range(4)] for m in range(4)]

def rfk_qmul(a, b):
    return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
            a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
            a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
            a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])

def tau_r7_factors():
    """Exact dyadic parse of the frozen tau_r7 certificate factors,
    layout z = [a:3x7; b:4x7; c:4x7] (canonical Sys layout)."""
    import csv as _csv
    def load(fn):
        with open(os.path.join(CERTDIR, fn)) as f:
            return [[rfk.dy(v) for v in rec] for rec in _csv.reader(f)]
    A7, B7, C7 = load("A.csv"), load("B.csv"), load("C.csv")
    z = ([A7[p][rr] for p in range(3) for rr in range(7)]
         + [B7[bb][rr] for bb in range(4) for rr in range(7)]
         + [C7[cc][rr] for cc in range(4) for rr in range(7)])
    return z, A7, B7, C7

def merge_seed_from_tau(A7, B7, C7, j):
    """Merge tau-column j with tau-column (j+7) mod 14 of the block-
    duplicated 14-term layout; returns a float seed (canonical layout)."""
    R = 13
    # block-duplicated 14-term factors: for term s in 0..13,
    #   a[s] = A7[:, s % 7] (3-vector), b = (B7 rows; same) with block copy,
    #   c likewise.  tau-column j (top j) merges with top (j+7)%14... here:
    # merged term: same b-column and c-column for the two blocks, single a.
    # Layout: ranks 0..12; merged term sits at index 12.
    a = np.zeros((3, R)); b = np.zeros((8, R)); c = np.zeros((8, R))
    s = 0
    A7f = np.array([[float(v) for v in row] for row in A7])
    B7f = np.array([[float(v) for v in row] for row in B7])
    C7f = np.array([[float(v) for v in row] for row in C7])
    merged_top = j            # top-block term j (0..6)
    merged_bot = j            # bottom-block term j (same tau column)
    for t in range(7):
        if t == merged_top:
            continue
        # top block term t -> rank s
        a[:, s] = A7f[:, t]; b[:4, s] = B7f[:, t]; c[:4, s] = C7f[:, t]
        s += 1
    for t in range(7):
        if t == merged_bot:
            continue
        a[:, s] = A7f[:, t]; b[4:, s] = B7f[:, t]; c[4:, s] = C7f[:, t]
        s += 1
    # merged rank 12: single 3-vector a serving both blocks with the SAME
    # b/c columns (the merge hypothesis)
    a[:, 12] = A7f[:, merged_top]
    b[:4, 12] = B7f[:, merged_top]; c[:4, 12] = C7f[:, merged_top]
    b[4:, 12] = B7f[:, merged_bot]; c[4:, 12] = C7f[:, merged_bot]
    return np.concatenate([a.ravel(), b.ravel(), c.ravel()])

def main():
    out = {"env": sys.version.split()[0],
           "prereg_commit": "d6c7e44d6901eb025835313fe8bf9c7b01ee3a81"}
    Ttbl, TF, tauT = build_targets()
    sysmain = rfk.Sys(TF, 13, None, "MAIN")
    log(f"[ANCHOR] table ok; TF nonzero={sum(1 for s in TF for rowv in s for v in rowv if v)}")

    # ---------------- 1. controls (exact same paths as frozen instrument)
    z_tau7, A7, B7, C7 = tau_r7_factors()
    resc7 = rfk.certify(rfk.Sys(tauT, 7, None, "TAU7"), z_tau7, "qrcp", "N1-TAU7")
    log(f"[CONTROL] TAU7 verdict={resc7['verdict']} pass_rung={resc7['pass_rung']}")
    z_tau6_f = np.concatenate([
        np.asarray(A7, dtype=object)[:, :6].astype(float).ravel(),
        np.asarray(B7, dtype=object)[:, :6].astype(float).ravel(),
        np.asarray(C7, dtype=object)[:, :6].astype(float).ravel()])
    systau6 = rfk.Sys(tauT, 6, None, "TAU6")
    z_tau6, _ = polish(systau6, z_tau6_f, "N1-TAU6")
    resc6 = rfk.certify(systau6, to_dy(z_tau6), "qrcp", "N1-TAU6")
    log(f"[CONTROL] TAU6 verdict={resc6['verdict']}")

    rngA = np.random.default_rng(20260901)
    astar = rngA.normal(scale=0.7, size=(3, 13))
    bstar = rngA.normal(scale=0.7, size=(8, 13))
    cstar = rngA.normal(scale=0.7, size=(8, 13))
    astr_f, bstr_f, cstr_f = astar.ravel(), bstar.ravel(), cstar.ravel()
    a_dy = [dy_from_float(v) for v in astr_f]
    b_dy = [dy_from_float(v) for v in bstr_f]
    c_dy = [dy_from_float(v) for v in cstr_f]
    U_dy = [[[sum(a_dy[pp * 13 + rr] * b_dy[bb * 13 + rr]
                  * c_dy[cc * 13 + rr] for rr in range(13))
              for cc in range(8)] for bb in range(8)] for pp in range(3)]
    syspos = rfk.Sys(U_dy, 13, None, "PPOS")
    z_rootf = np.concatenate([astr_f, bstr_f, cstr_f])
    from scipy import linalg as sla
    Jf_root = syspos.jac_f(z_rootf)
    _, _, piv_r = sla.qr(Jf_root, mode="economic", pivoting=True)
    Scols = piv_r[:syspos.M].tolist()
    z_plant = z_rootf.copy()
    rngB = np.random.default_rng(20260902)
    pert = 1e-8 * rngB.normal(size=syspos.NV)
    for col in Scols:
        z_plant[col] += pert[col]
    respos = rfk.certify(syspos, to_dy(z_plant), "qrcp", "N1-PPOS")
    log(f"[CONTROL] PPOS verdict={respos['verdict']} pass_rung={respos['pass_rung']}")

    # P-WRNG probe seed: the frozen polished candidate (same as frozen inst.)
    rows = [r for r in csv.reader(open(os.path.join(ROUTEF,
            "candidate_main_factors.csv"))) if r and r[0] == "FACTOR"]
    z_frozen = np.array([float(v) for v in rows[0][1:]])
    Tcor = [[[v for v in rowv] for rowv in mat] for mat in TF]
    Tcor[0][0][0] = 2
    reswrng = rfk.certify(rfk.Sys(Tcor, 13, None, "PWRNG"),
                          to_dy(z_frozen), "qrcp", "N1-PWRNG")
    log(f"[CONTROL] PWRNG verdict={reswrng['verdict']}")

    behavior_ok = (resc7["verdict"] == "CONTAINMENT"
                   and resc6["verdict"] != "CONTAINMENT"
                   and respos["verdict"] == "CONTAINMENT"
                   and reswrng["verdict"] != "CONTAINMENT")
    out["controls"] = {"tau7": resc7["verdict"], "tau6": resc6["verdict"],
                       "ppos": respos["verdict"], "pwrng": reswrng["verdict"],
                       "behavior_ok": behavior_ok}
    if not behavior_ok:
        log("[CONTRACT] CONTROL BREAK — INVALID INSTRUMENT; main sweep ABORTED")
        out["verdict"] = "INVALID INSTRUMENT"
        with open(os.path.join(HERE, "n1_results.json"), "w") as f:
            json.dump(out, f, indent=1)
        sys.exit(5)

    # ---------------- 2. main sweep: 40 seeds
    seeds = []   # (class, label, float vector)
    for k in range(32):
        rng = np.random.default_rng(20260901 + k)
        a = rng.normal(size=(3, 13)) / np.sqrt(8)
        b = rng.normal(size=(8, 13)) / np.sqrt(8)
        c = rng.normal(size=(8, 13)) / np.sqrt(8)
        scale = 4.898979485566356 ** (1.0 / 3.0)
        seeds.append(("R", f"R{k:02d}", np.concatenate([(a * scale).ravel(),
                     (b * scale).ravel(), (c * scale).ravel()])))
    for j in range(7):
        seeds.append(("M", f"M{j}", merge_seed_from_tau(A7, B7, C7, j)))
    seeds.append(("F", "F_frozen", z_frozen))

    admitted = []
    summary = []
    for cls, label, z0 in seeds:
        if cpu_left() < 300:
            log("[BUDGET] CPU cap nearing — remaining seeds SKIPPED (recorded)")
            summary.append({"cls": cls, "label": label, "skipped": "budget"})
            continue
        z, hist = polish(sysmain, z0, f"N1-{label}")
        z_dy = to_dy(z)
        gm, _ = gmax_exact(sysmain, z_dy)
        exact = exact_witness_check(sysmain, z_dy)
        row = {"cls": cls, "label": label, "rel_final": sysmain.rel(z),
               "gmax_exact": str(gm), "exact": exact, "hist": hist}
        summary.append(row)
        log(f"[N1-{label}] cls={cls} rel_final={row['rel_final']:.3e} "
            f"gmax_exact={gm} admitted={gm <= rfk.fmpq(1, 10**7)}")
        if gm <= rfk.fmpq(1, 10**7):
            admitted.append((cls, label, z_dy, row))
        with open(os.path.join(HERE, "n1_progress.json"), "w") as f:
            json.dump(out | {"summary": summary}, f, indent=1)

    out["summary"] = summary
    out["n_admitted"] = len(admitted)

    # ---------------- 3./4. admitted: exact Newton then ladder
    out["admitted_outcomes"] = []
    main_verdict = "FAILURE TO CERTIFY"
    for cls, label, z_dy, row in admitted:
        zn, nh, is_witness = exact_newton(sysmain, z_dy)
        outcome = {"cls": cls, "label": label, "newton_hist": nh,
                   "witness": is_witness}
        if is_witness:
            gm, _ = gmax_exact(sysmain, zn)
            outcome["gmax_after_newton"] = str(gm)
            out["witness_factors"] = [str(v) for v in zn]
            main_verdict = "RANK 13 EXACTLY (exact witness)"
            log(f"[N1-{label}] EXACT WITNESS: ||g||_max = {gm}")
            out["admitted_outcomes"].append(outcome)
            break
        resn = rfk.certify(sysmain, zn, "qrcp", f"N1-{label}-NEWTON")
        outcome["ladder"] = {"verdict": resn["verdict"],
                             "pass_rung": resn["pass_rung"],
                             "excl_rung": resn["excl_rung"]}
        log(f"[N1-{label}] NEWTON ladder: {resn['verdict']} "
            f"pass={resn['pass_rung']}")
        if resn["verdict"] == "CONTAINMENT":
            main_verdict = "RANK 13 EXACTLY (containment)"
            out["containment_candidate"] = [str(v) for v in zn]
            out["admitted_outcomes"].append(outcome)
            break
        out["admitted_outcomes"].append(outcome)

    out["verdict"] = main_verdict
    out["cpu_s"] = time.process_time() - T0
    with open(os.path.join(HERE, "n1_results.json"), "w") as f:
        json.dump(out, f, indent=1)
    log(f"[DONE] verdict={main_verdict} cpu={out['cpu_s']:.1f}s")
    sys.exit(0)

if __name__ == "__main__":
    main()
