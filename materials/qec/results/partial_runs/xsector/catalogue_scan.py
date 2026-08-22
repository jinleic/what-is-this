"""Catalogue-wide X-sector scan: dims, demotion, and exact d_X(Q) vs d_X(P).

For every delta>0 catalogue row (155 of 368):
  - dim Xcen(P), dim Xcen(Q), rho_X, dim S_X(P), dim S_X(Q), dim Deltabar
  - identity checks (Theorem J.3)
  - demotion_dim = dim(rowspace H_X cap ker[C D]) - dim S_X(Q)
  - exact d_X(Q) by coset enumeration when dim Xcen(Q) <= 22
  - parent d_Z(P) (=> d_X(P) by reversal symmetry J.0) from EXP-039 certs.
GF(2)/numpy only, one thread, no SAT.
"""
from __future__ import annotations

import importlib.util
import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

_spec39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py"
)
E39 = importlib.util.module_from_spec(_spec39)
sys.modules[_spec39.name] = E39
_spec39.loader.exec_module(E39)

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np  # noqa: E402


def mats(row):
    ell, m = int(row["ell"]), int(row["m"])
    A = poly_matrix(ell, m, [tuple(t) for t in row["A_terms"]])
    B = poly_matrix(ell, m, [tuple(t) for t in row["B_terms"]])
    C = poly_matrix(ell, m, [tuple(t) for t in (row.get("C_terms") or [])])
    D = poly_matrix(ell, m, [tuple(t) for t in (row.get("D_terms") or [])])
    return ell, m, A, B, C, D


def min_wt_coset(basis, exclude, cap_bits=22):
    d = basis.shape[0]
    if d > cap_bits:
        return None
    exn = nullspace_np(exclude) if (exclude is not None and exclude.shape[0]) else None
    best = None
    for start in range(1, 1 << d, 1 << 18):
        idx = np.arange(start, min(start + (1 << 18), 1 << d), dtype=np.int64)
        V = np.zeros((idx.size, basis.shape[1]), dtype=np.uint8)
        for j in range(d):
            V ^= (((idx >> j) & 1)[:, None] & basis[j]).astype(np.uint8)
        if exn is not None and exn.shape[0]:
            V = V[~(((V @ exn.T) % 2).sum(axis=1) == 0)]
        if V.shape[0]:
            b = int(V.sum(axis=1).min())
            best = b if best is None else min(best, b)
    return best


def analyse(row):
    ell, m, A, B, C, D = mats(row)
    dim, N = ell * m, 2 * ell * m
    Z2 = np.zeros((dim, dim), np.uint8)
    HX, HZ, CD = np.hstack([A, B]), np.hstack([B.T, A.T]), np.hstack([C, D])
    HQ = np.vstack([np.hstack([A, B, C, D]), np.hstack([Z2, Z2, B.T, A.T])])
    kP = N - rank_np(HX) - rank_np(HZ)
    kQ = N - rank_np(HQ)
    L = nullspace_np(HX.T)
    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, N), np.uint8)
    dbar = rank_np(np.vstack([HZ, Delta])) - rank_np(HZ)
    XcenP = nullspace_np(HZ)
    XcenQ = nullspace_np(np.vstack([HZ, CD]))
    rho = XcenP.shape[0] - XcenQ.shape[0]
    SZD = np.vstack([HZ, Delta]) if Delta.shape[0] else HZ
    Wperp = nullspace_np(SZD)
    G = (CD @ Wperp.T) % 2
    Lam = nullspace_np(G.T)
    SXQ = (Lam @ HX) % 2 if Lam.shape[0] else np.zeros((0, N), np.uint8)
    dimSXQ = rank_np(SXQ)
    kerCD = nullspace_np(CD)
    inter = rank_np(HX) + rank_np(kerCD) - rank_np(np.vstack([HX, kerCD]))
    dem = inter - dimSXQ
    ok = {
        "kQ": kQ == kP - dbar,
        "XcenQ": XcenQ.shape[0] == kQ + dimSXQ,
        "Xbar": XcenQ.shape[0] - dimSXQ == kQ,
        "SX_shrink": (rank_np(HX) - dimSXQ) == (rho - dbar),
    }
    dXQ = min_wt_coset(XcenQ, SXQ) if XcenQ.shape[0] <= 22 else None
    return {
        "n": N, "k_P": kP, "k_Q": kQ, "dim_Deltabar": dbar,
        "dim_XcenP": XcenP.shape[0], "dim_XcenQ": XcenQ.shape[0], "rho_X": rho,
        "dim_SXP": rank_np(HX), "dim_SXQ": dimSXQ, "demotion_dim": dem,
        "identities": all(ok.values()), "identity_detail": ok,
        "d_X_Q_exact": dXQ, "dim_XcenQ_basis": XcenQ.shape[0],
    }


def main():
    rows = E27.load_catalogue()
    certs = E39.load_certificates()
    fp_cache = {}
    results = []
    n_with_demotion = 0
    n_exact_dXQ = 0
    for row in rows:
        if not ((row.get("C_terms") or []) + (row.get("D_terms") or [])):
            continue
        label = row.get("code_id")
        out = analyse(row)
        key = (int(row["ell"]), int(row["m"]),
               tuple(sorted(map(tuple, row["A_terms"]))),
               tuple(sorted(map(tuple, row["B_terms"]))))
        if key not in fp_cache:
            _, HX, HZ = E27.parent_matrices(row)
            fp_cache[key] = E27.matrix_fingerprint(HX, HZ)
        cert = certs.get(fp_cache[key])
        dzP = cert.get("d_z_parent") if cert else None
        dzP_exact = cert.get("d_z_exact") if cert else None
        out.update({"label": label, "d_Z_P_cert": dzP, "d_Z_P_cert_exact": dzP_exact,
                    "cat_d": row.get("d"), "cat_d_exact": row.get("d_is_exact")})
        if out["demotion_dim"] > 0:
            n_with_demotion += 1
        if out["d_X_Q_exact"] is not None:
            n_exact_dXQ += 1
        results.append(out)
        interesting = out["demotion_dim"] > 0 or (
            out["d_X_Q_exact"] is not None and dzP is not None and out["d_X_Q_exact"] < dzP)
        if interesting:
            print("HIT", json.dumps(out), flush=True)
    print(json.dumps({"rows_scanned": len(results), "with_demotion": n_with_demotion,
                      "exact_dXQ": n_exact_dXQ,
                      "identity_failures": sum(0 if r["identities"] else 1 for r in results)}))
    (ROOT / "results/partial_runs/xsector/catalogue_scan.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
