"""gate_b_deflate.py — structured rank-24 attempts starting from the
verified rank-25 certificate x0 (deflation from an exact solution),
NOT random starts. Per Main's redirection (optional, under declared
budget). Attempts, all jointly instrumented (residual + factor norms):

  MODE=merge     : freeze 24 terms of x0 as FIXED, drive the 25th term's
                   scale s -> 0 by 1-D continuation (residual is quadratic
                   in s); at the s* minimizing the model, LM on all 600
                   coords from that neighbourhood.
  MODE=compress  : CD-structured merge — find the pair of columns of the
                   certified A/B/C with minimal "merge distortion"
                   || a_u b_u^T - a_v b_v^T || style score, merge them
                   (u <- (u+v)/2), run LM.
  MODE=alphabet  : snap each factor column scale to a small rational grid,
                   run exact residual check at 600 coords (quick probe);
                   no bound claim — proximity evidence only.

Output: JSON rows appended to stdout; best candidate npz.  No Krawczyk
attempt unless rel <= 1e-13 (escalation rule: message Main first).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gate_b_search3 import T, TF, N1, N2, N3, M, residual, jacobian, \
    factor_stats  # noqa: E402

CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "scratch", "upstream_ref", "certs", "rank25_cert")
MODE = os.environ.get("MODE", "merge")
R = 25


def load_cert_x0():
    import csv
    from fractions import Fraction
    A = []
    with open(os.path.join(CERT, "A.csv")) as f:
        for rec in csv.reader(f):
            A.append([float(v) for v in rec])
    B = []
    with open(os.path.join(CERT, "B.csv")) as f:
        for rec in csv.reader(f):
            B.append([float(v) for v in rec])
    C = []
    with open(os.path.join(CERT, "C.csv")) as f:
        for rec in csv.reader(f):
            C.append([float(v) for v in rec])
    x0 = np.concatenate([np.asarray(A).ravel(order="C"),
                         np.asarray(B).ravel(order="C"),
                         np.asarray(C).ravel(order="C")])
    return x0


def split(x, r=R):
    oB, oC = N1 * r, (N1 + N2) * r
    return (x[:oB].reshape(N1, r), x[oB:oC].reshape(N2, r),
            x[oC:].reshape(N3, r))


def lm_refine(x0, r, max_rounds=400):
    from scipy.optimize import least_squares
    res = least_squares(lambda z: residual(z, r), x0,
                        jac=lambda z: jacobian(z, r), method="trf",
                        tr_solver="lsmr", x_scale="jac",
                        max_nfev=max_rounds,
                        xtol=3e-16, ftol=3e-16, gtol=3e-16)
    x = np.asarray(res.x, dtype=float)
    rel = float(np.linalg.norm(residual(x, r)) / TF)
    mx, per, ratio, prod = factor_stats(x, r, rel)
    return x, rel, int(res.nfev), int(res.status), mx, per, ratio, prod


def main():
    t0 = time.time()
    x0 = load_cert_x0()
    rel0 = float(np.linalg.norm(residual(x0, R)) / TF)
    mx0, per0, ratio0, prod0 = factor_stats(x0, R, rel0)
    print(json.dumps({"phase": "x0_check", "rank": R, "rel": rel0,
                      "factor_norms": per0, "max_factor_norm": mx0}))
    out = {"x0": {"rel": rel0, "factor_norms": per0}}
    best = None

    if MODE == "merge":
        # continuation: scale one term toward 0 while the other 24 stay free
        a, b, c = split(x0)
        for term in range(R):
            # 1-D quadratic model in s: residual(s) = F0 + s D + s^2 E
            # where D = dF/ds at s=1 along (a_t b_t c_t scaled), E fixed.
            # Direct approach: for grid of s, do full LM from perturbed x0.
            s_grid = [0.9, 0.7, 0.5, 0.3, 0.1, 0.05, 0.0]
            for s in s_grid:
                xs = x0.copy()
                oB, oC = N1 * R, (N1 + N2) * R
                sl = slice(term, term + 1)
                for blk, off in ((0, 0), (1, oB), (2, oC)):
                    idx0 = off + term
                    # scale term's column in each factor block
                    xs = xs.copy()
                # scale columns of A, B, C by s for this term:
                aa, bb, cc = split(xs)
                aa[:, term] *= s
                # NOTE: scaling only one factor's column keeps the RANK-1
                # term scaled by s (not s^3): correct deflation geometry.
                xs = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
                rel_start = float(np.linalg.norm(residual(xs, R)) / TF)
                x, rel, nfev, status, mx, per, ratio, prod = lm_refine(xs, R)
                row = {"mode": "merge", "term": term, "s": s,
                       "rel_after_LM": rel, "nfev": nfev, "status": status,
                       "max_factor_norm": mx, "factor_norms": per}
                print(json.dumps(row), flush=True)
                if best is None or best[0] > rel:
                    best = (rel, x, row)
                if rel <= 1e-13:
                    np.save("gate_b_deflate_candidate.npz", x, rel)
                    print("SUB-1E-13 — STOP, escalate before writing")
                    return
        # try full-term annihilation with all three factors scaled by s^(1/3)
    elif MODE == "compress":
        a, b, c = split(x0)
        # score column pairs by how close their outer products are
        scores = []
        for u in range(R):
            for v in range(u + 1, R):
                d = (np.linalg.norm(np.outer(a[:, u], b[:, u])
                                    - np.outer(a[:, v], b[:, v]))
                     + np.linalg.norm(np.outer(a[:, u], c[:, u])
                                      - np.outer(a[:, v], c[:, v])))
                scores.append((float(d), u, v))
        scores.sort()
        for d, u, v in scores[:5]:
            xc = x0.copy()
            aa, bb, cc = split(xc)
            for Mat in (aa, bb, cc):
                Mat[:, u] = 0.5 * (Mat[:, u] + Mat[:, v])
            # zero column v: reduces CP rank to 24 terms (exactly)
            for Mat in (aa, bb, cc):
                Mat[:, v] = 0.0
            xs = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
            rel_start = float(np.linalg.norm(residual(xs, 25)) / TF)
            # now treat as rank-24 by dropping zero column
            aa24 = aa[:, [k for k in range(25) if k != v]]
            bb24 = bb[:, [k for k in range(25) if k != v]]
            cc24 = cc[:, [k for k in range(25) if k != v]]
            x24 = np.concatenate([aa24.ravel(), bb24.ravel(), cc24.ravel()])
            rel0 = float(np.linalg.norm(residual(x24, 24)) / TF)
            x, rel, nfev, status, mx, per, ratio, prod = lm_refine(x24, 24)
            row = {"mode": "compress", "pair": [u, v], "merge_dist": d,
                   "rel_start_r24": rel0, "rel_after_LM": rel,
                   "nfev": nfev, "status": status, "max_factor_norm": mx,
                   "factor_norms": per}
            print(json.dumps(row), flush=True)
            if best is None or best[0] > rel:
                best = (rel, x, row)
            if rel <= 1e-13:
                print("SUB-1E-13 — STOP, escalate before writing")
                np.savez("gate_b_deflate_candidate.npz", x=x, rel=rel)
                return
    elif MODE == "alphabet":
        grid = [-2, -1, -0.5, 0, 0.5, 1, 2]
        a, b, c = split(x0)
        # snap each entry to nearest grid value, check residual damage
        def snap(Mm):
            Mg = np.zeros_like(Mm)
            for i in range(Mm.shape[0]):
                for j in range(Mm.shape[1]):
                    Mg[i, j] = min(grid, key=lambda g: abs(Mm[i, j] - g))
            return Mg
        best_snap = None
        for which in ("A", "B", "C", "ABC"):
            aa, bb, cc = snap(a), b, c
            if which in ("B", "ABC"):
                bb = snap(b)
            if which in ("C", "ABC"):
                cc = snap(c)
            xs = np.concatenate([aa.ravel(), bb.ravel(), cc.ravel()])
            rel = float(np.linalg.norm(residual(xs, R)) / TF)
            row = {"mode": "alphabet", "snapped": which, "rel": rel}
            print(json.dumps(row), flush=True)
            if best_snap is None or rel < best_snap[0]:
                best_snap = (rel, row)
        out["alphabet_best"] = best_snap[1]
        return

    if best is not None:
        rel, x, row = best
        np.savez("gate_b_deflate_best.npz", x=x, rel=rel, mode=MODE)
        print("BEST " + json.dumps({"rel": rel, "row": row}))


if __name__ == "__main__":
    main()
