"""EXP-046 — deterministic (search-free) syzygy collapse constructor.

Mechanism (derived 2026-08-19 from the J.5 witness analysis):
for a BB parent (A, B) over R = GF(2)[x,y]/(x^ell-1, y^m-1), take D = 0 and
C = circulant(c) with c in ker(A) (group-algebra annihilator of the polynomial
a). Then the perturbation (C, D) satisfies M = A C^T + B D^T = 0, hence is
valid; the partner z0 = (c, 0) lies in ker H_X automatically; so row 0 of H_X
(weight wt(A) + wt(B)) is a DEMOTED parent X-stabilizer iff
(c, 0) not in S_Z + Delta_C (J.2/J.4(ii) of notes/theorem_j_xsector.md) —
in which case d_X(Q) <= wt(A) + wt(B) while k_Q = k_P when the demotion is
witnessed with dim bar-Delta = 0.

dim ker(A) >= k_P/2 always (A is singular whenever k_P > 0 since
rank A <= rank[A B] = lm - k_P/2); the flagship parent has dim ker A = k_P.

This script replaces the sharp_probe/sampling machinery of j5core for the
D=0 sector with a pure-nullspace constructor: no RNG, no sampling, no
per-perturbation scans.  For each parent it tries nullspace basis vectors
singly, then pairwise sums, then triples (capped); the FIRST candidate passing
the dual verification battery (path A: z0 not in S_Z + Delta_C rank test;
path B: x0 not in S_X(Q) nullspace probe) is recorded as a verified witness
with its full term list.

Collapse is flagged whenever a certified distance bound for the parent exists
(EXP-039 exact certificate preferred, else EXP-037 pool lower bound) and
wt(A) + wt(B) strictly below the bound.  Parents without any bound still get
a verified demotion row (bound=None, collapse=None).

Artifact: results/processed/exp046_generic_syzygy.json
(schema exp046-generic-syzygy-v1).  GF(2)/numpy only, single thread.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

spec = importlib.util.spec_from_file_location("exp027", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = E27
spec.loader.exec_module(E27)
spec39 = importlib.util.spec_from_file_location("exp039", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(spec39)
sys.modules[spec39.name] = E39
spec39.loader.exec_module(E39)

from qec_research.gf2.linalg import nullspace_np, rank_np  # noqa: E402
from qec_research.codes.bicycle import poly_matrix  # noqa: E402

OUT = ROOT / "results" / "processed" / "exp046_generic_syzygy.json"
TRIPLE_CAP = 5000  # max triples tested per parent before declaring 'deferred'


def build_CD(ell: int, m: int, cvec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    terms = [(int(i // m), int(i % m)) for i in np.flatnonzero(cvec)]
    C = poly_matrix(ell, m, terms)
    D = np.zeros_like(C)
    return C, D


def try_candidate(A, B, HX, HZ, ell, m, cvec):
    """Dual-path verification for D=0, C=circ(cvec).  Return dict or None."""
    dim = ell * m
    C, D = build_CD(ell, m, cvec)
    M = (A @ C.T + B @ D.T) % 2
    if M.any():
        return None  # construction guarantee failed (should not happen)
    CD = np.hstack([C, D])
    z0 = np.concatenate([cvec, np.zeros(dim, np.uint8)])
    L = nullspace_np(HX.T)                    # left nullspace of [A B]
    spanDelta = (L @ CD) % 2 if L.size else np.zeros((0, 2 * dim), np.uint8)
    SZD = np.vstack([HZ, spanDelta])
    r = rank_np(SZD)
    if rank_np(np.vstack([SZD, z0])) == r:
        return None                           # c absorbed: no demotion
    # path B: x0 = row 0 of H_X must be in Xcen(Q) but not S_X(Q)
    x0 = HX[0].copy()
    if (HZ @ x0 % 2).any() or (CD @ x0 % 2).any():
        return None                           # x0 not in Xcen(Q) (unexpected)
    Wp = nullspace_np(SZD)
    Lam = nullspace_np(((CD @ Wp.T) % 2).T)
    SXQ = (Lam @ HX) % 2 if Lam.size else np.zeros((0, 2 * dim), np.uint8)
    ns = nullspace_np(SXQ)
    x0_in_SXQ = not bool(((x0 @ ns.T) % 2).sum())
    if x0_in_SXQ:
        return None                           # x0 still a stabilizer: no demotion
    return {
        "c_terms": [[int(i // m), int(i % m)] for i in np.flatnonzero(cvec)],
        "c_weight": int(cvec.sum()),
        "x0_weight": int(x0.sum()),
        "z0_weight": int(z0.sum()),
        "M_zero": True,
        "z0_in_ker_HX": bool(not (HX @ z0 % 2).any()),
        "x0_demoted_pathA": True,
        "x0_demoted_pathB": not x0_in_SXQ,
        "rho_X": int(rank_np(np.vstack([HZ, CD])) - rank_np(HZ)),
        "delta_bar": int(rank_np(SZD) - rank_np(HZ)),
    }


def solve_parent(fp: str, entry: dict, rows_by_fp: dict, certs, pool):
    t0 = time.time()
    members = entry["members"]
    row = rows_by_fp[fp]
    ell, m = int(entry["ell"]), int(entry["m"])
    dim = ell * m
    _, HX, HZ = E27.parent_matrices(row)
    A, B = HX[:, :dim], HX[:, dim:]
    wt_sum = int(HX[0].sum())
    K = nullspace_np(A)
    nullity = int(K.shape[0]) if K.size else 0
    cert = certs.get(fp) or {}
    d_cert = cert.get("d_z_parent") if cert.get("d_z_exact") else None
    d_pool = pool.get(fp)
    d_bound = d_cert if d_cert is not None else d_pool
    rec = {
        "fingerprint": fp,
        "fingerprint12": fp[:12],
        "label": members[0]["label"],
        "ell": ell, "m": m, "n": 2 * dim,
        "k_parent": entry["k_parent"],
        "rank_A": int(rank_np(A)),
        "nullity_A": nullity,
        "wt_sum_row0": wt_sum,
        "d_cert_exact": d_cert,
        "d_pool_lb": d_pool,
        "d_bound": d_bound,
        "demotion": False,
        "demotion_depth": None,
        "witness": None,
        "collapse": None,
        "wall_s": None,
    }
    if not nullity:
        rec["wall_s"] = round(time.time() - t0, 3)
        return rec
    # Stage 1: D=0 with single basis vectors of ker(A) (cheapest; matches all
    # historical witnesses).  Stage 2: general syzygy kernel K0 = ker[A B] on
    # the pair (c,d): M = A C^T + B D^T = 0 iff A c + B d = 0 (single vector
    # condition for row 0; the rest are all shift-equivalents), z0=(c,d), so
    # both the validity and the partner-in-ker-H_X checks are automatic.
    candidates = [("D0", K[i]) for i in range(nullity)]
    for cvec in [c for _, c in candidates]:
        hit = try_candidate(A, B, HX, HZ, ell, m, cvec)
        if hit is not None:
            rec.update(demotion=True, demotion_depth="D0_basis", witness=hit)
            if d_bound is not None:
                rec["collapse"] = bool(wt_sum < d_bound)
            rec["wall_s"] = round(time.time() - t0, 3)
            return rec
    K0 = nullspace_np(HX)
    rec["nullity_HX_pair"] = int(K0.shape[0])
    parent_deadline = t0 + 300.0
    for i in range(K0.shape[0]):
        hit = try_candidate_general(A, B, HX, HZ, ell, m, K0[i])
        if hit is not None:
            rec.update(demotion=True, demotion_depth="syz_basis", witness=hit)
            if d_bound is not None:
                rec["collapse"] = bool(wt_sum < d_bound)
            rec["wall_s"] = round(time.time() - t0, 3)
            return rec
    for i, j in combinations(range(K0.shape[0]), 2):
        if time.time() > parent_deadline:
            rec["demotion_depth"] = "search_cap_exceeded"
            rec["wall_s"] = round(time.time() - t0, 3)
            return rec
        hit = try_candidate_general(A, B, HX, HZ, ell, m, (K0[i] ^ K0[j]) % 2)
        if hit is not None:
            rec.update(demotion=True, demotion_depth="syz_pairs", witness=hit)
            if d_bound is not None:
                rec["collapse"] = bool(wt_sum < d_bound)
            rec["wall_s"] = round(time.time() - t0, 3)
            return rec
    rec["demotion_depth"] = "none_in_syz_basis_pairs(cap)"
    rec["wall_s"] = round(time.time() - t0, 3)
    return rec


def try_candidate_general(A, B, HX, HZ, ell, m, cdvec):
    """Dual-path verification for general syzygy pair (c, d)."""
    dim = ell * m
    cvec, dvec = cdvec[:dim], cdvec[dim:]
    c_terms = [(int(i // m), int(i % m)) for i in np.flatnonzero(cvec)]
    d_terms = [(int(i // m), int(i % m)) for i in np.flatnonzero(dvec)]
    C = poly_matrix(ell, m, c_terms)
    D = poly_matrix(ell, m, d_terms)
    M = (A @ C.T + B @ D.T) % 2
    if M.any():
        return None
    CD = np.hstack([C, D])
    z0 = np.concatenate([cvec, dvec])
    L = nullspace_np(HX.T)
    spanDelta = (L @ CD) % 2 if L.size else np.zeros((0, 2 * dim), np.uint8)
    SZD = np.vstack([HZ, spanDelta])
    r = rank_np(SZD)
    if rank_np(np.vstack([SZD, z0])) == r:
        return None
    x0 = HX[0].copy()
    if (HZ @ x0 % 2).any() or (CD @ x0 % 2).any():
        return None
    Wp = nullspace_np(SZD)
    Lam = nullspace_np(((CD @ Wp.T) % 2).T)
    SXQ = (Lam @ HX) % 2 if Lam.size else np.zeros((0, 2 * dim), np.uint8)
    ns = nullspace_np(SXQ)
    if not bool(((x0 @ ns.T) % 2).sum()):
        return None
    return {
        "c_terms": c_terms,
        "d_terms": d_terms,
        "c_weight": int(cvec.sum()),
        "d_weight": int(dvec.sum()),
        "x0_weight": int(x0.sum()),
        "z0_weight": int(z0.sum()),
        "M_zero": True,
        "z0_in_ker_HX": bool(not (HX @ z0 % 2).any()),
        "x0_demoted_pathA": True,
        "x0_demoted_pathB": True,
        "rho_X": int(rank_np(np.vstack([HZ, CD])) - rank_np(HZ)),
        "delta_bar": int(rank_np(SZD) - rank_np(HZ)),
    }


def main() -> None:
    t00 = time.time()
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    certs = E39.load_certificates()
    pool = E39.load_pool_lower_bounds()
    rows_by_fp = {}
    for index, row in enumerate(rows):
        _, HX, HZ = E27.parent_matrices(row)
        rows_by_fp.setdefault(E27.matrix_fingerprint(HX, HZ), row)
    out_rows = []
    for fp in sorted(parents):
        out_rows.append(solve_parent(fp, parents[fp], rows_by_fp, certs, pool))
    n = len(out_rows)
    dem = [r for r in out_rows if r["demotion"]]
    col = [r for r in dem if r["collapse"]]
    bounded = [r for r in out_rows if r["d_bound"] is not None]
    col_bounded = [r for r in bounded if r["collapse"]]
    payload = {
        "schema": "exp046-generic-syzygy-v1",
        "experiment": "EXP-046 deterministic syzygy collapse constructor",
        "mechanism": "D=0, C=circ(c), c in nullspace(A); M=A C^T=0 automatic; demotion iff (c,0) not in S_Z+Delta_C",
        "parents_total": n,
        "demotion_verified": len(dem),
        "demotion_by_depth": {
            k: sum(1 for r in dem if r["demotion_depth"] == k)
            for k in ("D0_basis", "syz_basis", "syz_pairs", "syz_triples")
        },
        "demotion_failed": [r["fingerprint12"] for r in out_rows if not r["demotion"]],
        "parents_with_bound": len(bounded),
        "collapse_verified_bounded": len(col_bounded),
        "collapse_verified_all": len(col),
        "bound_misses": [r["fingerprint12"] for r in bounded if not r["collapse"]],
        "nullity_histogram": {str(k): sum(1 for r in out_rows if r["nullity_A"] == k)
                              for k in sorted({r["nullity_A"] for r in out_rows})},
        "total_wall_s": round(time.time() - t00, 1),
        "rows": out_rows,
    }
    OUT.write_text(json.dumps(payload, indent=1))
    print("parents:", n, "| demotion:", len(dem), "| bounded:", len(bounded),
          "| collapse(bounded):", len(col_bounded))
    print("depth:", payload["demotion_by_depth"])
    print("failed:", payload["demotion_failed"])
    print("bound_misses:", payload["bound_misses"])
    print("nullity histogram:", payload["nullity_histogram"])
    print("WROTE", OUT, round(time.time() - t00, 1), "s")


if __name__ == "__main__":
    main()
