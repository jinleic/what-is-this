"""Theorem J (X-sector) verification: pure-X centralizer/stabilizer/logical dims.

Scratch probe for notes/theorem_j_xsector.md.  GF(2)/numpy only, one thread.

Conventions (match reports/paper_pbb_nogo.md):
  H_X = [A B], H_Z = [B^T A^T], H_Q = [[A B | C D], [0 0 | B^T A^T]].
  Pure-X centralizer of a code with check matrix H (rows in symplectic (x|z)):
      Xcen = {x : (x|0) commutes with every row of H}.
  Pure-X stabilizer space: S_X = {x : (x|0) in rowspace(H)}.
  Pure-X logicals: Xbar = Xcen / S_X   (as vector spaces via x <-> (x|0)).
"""
from __future__ import annotations

import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

_SPEC = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py"
)
E27 = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = E27
_SPEC.loader.exec_module(E27)

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np  # noqa: E402


def mats(row):
    ell, m = int(row["ell"]), int(row["m"])
    A = poly_matrix(ell, m, [tuple(t) for t in row["A_terms"]])
    B = poly_matrix(ell, m, [tuple(t) for t in row["B_terms"]])
    C = poly_matrix(ell, m, [tuple(t) for t in (row.get("C_terms") or [])])
    D = poly_matrix(ell, m, [tuple(t) for t in (row.get("D_terms") or [])])
    return ell, m, A, B, C, D


def sub_dim(U: np.ndarray, V: np.ndarray) -> int:
    """dim(U cap V) for row-space bases U, V."""
    return rank_np(U) + rank_np(V) - rank_np(np.vstack([U, V]))


def row_filter(M: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Rows of M restricted to the row space of keep-basis combos: here M rows
    whose coefficient vectors lie in rowspace(keep) are formed by projecting
    coefficients.  Simpler: return basis of {M[i] rows} over coeffs in rowspace(keep)."""
    raise NotImplementedError


def sector_dims(ell, m, A, B, C, D):
    dim = ell * m
    N = 2 * dim
    HX = np.hstack([A, B])
    HZ = np.hstack([B.T, A.T])
    CD = np.hstack([C, D])
    Z2 = np.zeros((dim, dim), np.uint8)
    HQ = np.vstack([np.hstack([A, B, C, D]), np.hstack([Z2, Z2, B.T, A.T])])

    out = {}
    out["n"] = N
    out["rank_HX"] = rank_np(HX)
    out["rank_HZ"] = rank_np(HZ)
    out["k_P"] = N - out["rank_HX"] - out["rank_HZ"]
    out["k_Q"] = N - rank_np(HQ)

    # ---- Z-side (Theorem G) ----
    L = nullspace_np(HX.T)             # {lambda : lambda [A B] = 0}  (left kernel)
    out["dim_L"] = L.shape[0]
    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, N), np.uint8)
    out["dim_Delta"] = rank_np(Delta)
    out["dim_Deltabar"] = rank_np(np.vstack([HZ, Delta])) - out["rank_HZ"]
    out["kQ_identity_holds"] = out["k_Q"] == out["k_P"] - out["dim_Deltabar"]

    # ---- X-side centralizers ----
    XcenP = nullspace_np(HZ)           # {x : H_Z x^T = 0}
    XcenQ = nullspace_np(np.vstack([HZ, CD]))
    out["dim_XcenP"] = XcenP.shape[0]
    out["dim_XcenQ"] = XcenQ.shape[0]
    out["rho_X"] = out["dim_XcenP"] - out["dim_XcenQ"]
    out["rho_identity_holds"] = out["rho_X"] == rank_np(np.vstack([HZ, CD])) - out["rank_HZ"]
    # membership cross-check: Xcen(Q) subset of Xcen(P) cap ker(CD)
    assert sub_dim(XcenQ, XcenP) == out["dim_XcenQ"]

    # ---- X-side stabilizers ----
    SX_P = HX
    out["dim_SXP"] = rank_np(SX_P)
    # Pure-X stabilizers of Q: row combos (lambda first block, mu second block)
    # with vanishing total Z-part:  x = lambda[A B], lambda[C D] in S_Z + Delta
    # (Delta = L[C D] accounts for representative shifts by the left kernel L).
    SZD = np.vstack([HZ, Delta]) if Delta.shape[0] else HZ   # S_Z + Delta basis
    Wperp = nullspace_np(SZD)          # rows: basis of (S_Z + Delta)^perp
    if C.any() or D.any():
        G = (CD @ Wperp.T) % 2         # lambda in Lambda' <=> lambda @ G == 0
        Lam = nullspace_np(G.T)        # {lambda : lambda[C D] in S_Z + Delta}
    else:
        Lam = np.eye(dim, dtype=np.uint8)
    SX_Q = (Lam @ HX) % 2
    out["dim_SXQ"] = rank_np(SX_Q)
    # sanity: every S_X(Q) row, as (x|0), commutes with every row of H_Q
    c_parts = HQ[:, N:]  # z-parts of all rows of H_Q (x|z blocks each length N)
    viol = (SX_Q @ c_parts.T) % 2 if SX_Q.shape[0] else np.zeros((1, 1), np.uint8)
    assert not viol.any(), "S_X(Q) not inside Xcen(Q)?!"

    # ---- dimension identities ----
    out["XcenP_eq_kP_plus_SXP"] = out["dim_XcenP"] == out["k_P"] + out["dim_SXP"]
    out["XcenQ_eq_kQ_plus_SXQ"] = out["dim_XcenQ"] == out["k_Q"] + out["dim_SXQ"]
    out["dim_Xbar_Q"] = out["dim_XcenQ"] - out["dim_SXQ"]
    out["dim_Zbar_Q"] = (N - out["rank_HX"]) - rank_np(np.vstack([HZ, Delta]))
    out["Xbar_eq_kQ"] = out["dim_Xbar_Q"] == out["k_Q"]
    out["Zbar_eq_kQ"] = out["dim_Zbar_Q"] == out["k_Q"]
    out["SX_shrink"] = out["dim_SXP"] - out["dim_SXQ"]
    out["SX_shrink_identity"] = out["SX_shrink"] == out["rho_X"] - out["dim_Deltabar"]
    out["rho_ge_deltabar"] = out["rho_X"] >= out["dim_Deltabar"]

    # ---- demotion space: rowspace(H_X) cap ker(CD) vs S_X(Q) ----
    if C.any() or D.any():
        kerCD = nullspace_np(CD)
        inter = sub_dim(HX, kerCD)       # dim rowspace(H_X) cap ker(CD)
    else:
        inter = out["dim_SXP"]
    out["dim_rowcapkerCD"] = inter
    out["demotion_dim"] = inter - out["dim_SXQ"] if inter >= out["dim_SXQ"] else None
    return out, XcenP, XcenQ, SX_Q, HX, HZ, CD


def min_weight_enumerate(basis: np.ndarray, exclude_basis: np.ndarray | None, cap_bits: int = 22):
    """Exact min weight of {v in span(basis) minus span(exclude_basis)}, v != 0."""
    d = basis.shape[0]
    if d > cap_bits:
        return None, d
    ex = exclude_basis if exclude_basis is not None and exclude_basis.shape[0] else None
    exn = nullspace_np(ex) if ex is not None else None
    best = None
    CH = 1 << 18
    for start in range(1, 1 << d, CH):
        idx = np.arange(start, min(start + CH, 1 << d), dtype=np.int64)
        V = np.zeros((idx.size, basis.shape[1]), dtype=np.uint8)
        for j in range(d):
            sel = (idx >> j) & 1
            V ^= (sel[:, None] & basis[j]).astype(np.uint8)
        if exn is not None and exn.shape[0]:
            # v in rowspace(ex) iff v . w = 0 for every w in nullspace(ex)
            inside = ((V @ exn.T) % 2).sum(axis=1) == 0
            V = V[~inside]
        elif ex is not None:
            return 0, d  # ex spans everything: no logicals at all
        if V.shape[0]:
            bmin = int(V.sum(axis=1).min())
            best = bmin if best is None else min(best, bmin)
    return best, d


def dX_by_support(HZ: np.ndarray, SXP_perp: np.ndarray, n: int, cap: int):
    """Min wt of {x : H_Z x^T = 0, x not in rowspace(H_X)} for wt <= cap.
    SXP_perp: basis of ker([A B]) (columns); x in rowspace H_X iff x.v=0 for all."""
    import itertools
    for w in range(1, cap + 1):
        CH = 200_000
        buf = itertools.islice(itertools.combinations(range(n), w), 10**9)
        while True:
            chunk = list(itertools.islice(buf, CH))
            if not chunk:
                break
            X = np.zeros((len(chunk), n), dtype=np.uint8)
            for r, sup in enumerate(chunk):
                X[r, list(sup)] = 1
            synd = (X @ HZ.T) % 2
            ok = ~synd.any(axis=1)
            if ok.any():
                Xin = X[ok]
                if SXP_perp.shape[0]:
                    inside = ((Xin @ SXP_perp.T) % 2).sum(axis=1) == 0
                    Xin = Xin[~inside]
                if Xin.shape[0]:
                    return w, Xin[[0]]
    return None, None


def main():
    rows = E27.load_catalogue()
    targets = ["phase2_58", "phase2_60", "phase2_88", "12_6_0193"]
    table = {}
    for row in rows:
        if row.get("code_id") not in targets and row.get("label") not in targets:
            continue
        label = row.get("code_id") or row.get("label")
        ell, m, A, B, C, D = mats(row)
        out, XcenP, XcenQ, SX_Q, HX, HZ, CD = sector_dims(ell, m, A, B, C, D)
        out["delta"] = len(row.get("C_terms") or []) + len(row.get("D_terms") or [])
        out["cat_k"] = row.get("k")
        out["cat_d"] = row.get("d")
        out["cat_d_exact"] = row.get("d_is_exact")
        # exact pure-X distance of Q when the coset is small enough
        dXQ, dd = min_weight_enumerate(XcenQ, SX_Q)
        out["d_X_Q_exact"] = dXQ
        out["d_X_Q_basisdim"] = dd
        # parent pure-X distance: exhaustive up to cap 4 (n <= 108 rows only)
        if out["n"] <= 108:
            Zcen = nullspace_np(HX)          # ker[A B] columns (rows of basis)
            w, wit = dX_by_support(HZ, Zcen, out["n"], 4)
            out["d_X_P_search4"] = w
            if wit is not None:
                out["d_X_P_witness_wt"] = int(wit.sum())
        table[label] = out
        print(label, json.dumps(out))
    (ROOT / "results/partial_runs/xsector").mkdir(parents=True, exist_ok=True)
    (ROOT / "results/partial_runs/xsector/crosstable.json").write_text(json.dumps(table, indent=1))


if __name__ == "__main__":
    main()
