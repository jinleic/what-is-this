"""Standalone reproduction of J.5 collapse witnesses from exp044_n180_hunt.json.

Rebuilds (A,B) from the catalogue, (C,D) from the stored witness terms, and
re-verifies the collapse with THREE independent GF(2) procedures:
  (1) rank test:   z = lam CD in rowspace([H_Z; Delta])?                 [probe]
  (2) S_X(Q) span: build S_X(Q) basis (exp-044 Wperp construction), test x not in it
  (3) direct solve: does (lam, mu) exist with lam H_X = x, lam CD + mu H_Z = 0?
                    (bitset solve on 2*dim unknowns, 2*N equations)
Plus: M symmetric, x in Xcen(Q), x in rowspace H_X, certified d_Z(P) bound.
Usage: python j5_reproduce_witness.py [witness_index ...]
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

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE.parent))

_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)
_spec39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(_spec39)
sys.modules[_spec39.name] = E39
_spec39.loader.exec_module(E39)

from qec_research.gf2.linalg import (rank_np, nullspace_np, rows_to_bitsets,  # noqa: E402
                                     row_space_contains, matmul as gf2matmul)
import j5core  # noqa: E402


def reproduce(wit, idx):
    fp_in = wit["parent_fingerprint"] if "parent_fingerprint" in wit else wit["fingerprint"]
    rows = E27.load_catalogue()
    certs = E39.load_certificates()
    parents = E39.distinct_parents(rows)
    fp = fp_in if fp_in in parents else next(f for f in parents if f.startswith(fp_in))
    entry = parents[fp]
    ell, m = entry["ell"], entry["m"]
    row0 = None
    for r in rows:
        if int(r["ell"]) == ell and int(r["m"]) == m:
            _, HXr, HZr = E27.parent_matrices(r)
            if E27.matrix_fingerprint(HXr, HZr) == fp:
                row0 = r
                break
    lat = j5core.Lattice(ell, m)
    p = j5core.ParentData(lat, row0["A_terms"], row0["B_terms"])
    C = lat.mat(wit["C_terms"])
    D = lat.mat(wit["D_terms"])
    # terms must rebuild the SAME matrices
    assert p.A_terms == tuple(map(tuple, wit["A_terms"])), "A terms mismatch vs witness"
    assert p.B_terms == tuple(map(tuple, wit["B_terms"])), "B terms mismatch vs witness"
    CD = np.hstack([C, D])
    lam0 = np.zeros(p.dim, np.uint8)
    lam0[list(wit["lam_support"])] = 1
    x0 = (lam0 @ p.HX) % 2
    z0 = (lam0 @ CD) % 2
    M = (gf2matmul(p.A, C.T) ^ gf2matmul(p.B, D.T)).astype(np.uint8)
    cert = certs.get(fp)
    d_lb = int(cert["d_z_parent"]) if cert and cert.get("d_z_exact") else None

    out = {"index": idx, "parent_fingerprint": fp[:16], "labels": wit.get("labels"),
           "ell": ell, "m": m, "k_P": int(p.kP), "d_Z_exact": d_lb,
           "checks": {}}
    c = out["checks"]
    c["weight_x"] = int(x0.sum())
    c["weight_z"] = int(z0.sum())
    c["M_symmetric"] = not bool((M ^ M.T).any())
    c["M_is_zero"] = not bool(M.any())
    c["lam_M_zero"] = not bool(((lam0 @ M) % 2).any())
    c["x_in_XcenQ(kerHZ,kerCD)"] = (not bool(((x0 @ p.HZ.T) % 2).any())
                                    and not bool(((x0 @ CD.T) % 2).any()))
    c["x_in_rowspace_HX"] = not bool(((x0 @ nullspace_np(p.HX).T) % 2).sum())
    # (1) rank probe
    sp = j5core.demotion_spaces(p, C, D)
    c["path1_z_notin_SZplusDelta"] = not j5core.z_in_SZplusDelta(z0, sp["probe_SZD"])
    # (2) S_X(Q) span (exp-044 construction)
    Wperp = nullspace_np(sp["SZD"])
    Lam = nullspace_np(((CD @ Wperp.T) % 2).T)
    SXQ = (Lam @ p.HX) % 2
    ns = nullspace_np(SXQ)
    c["path2_x_notin_SXQ"] = bool(((x0 @ ns.T) % 2).sum() != 0)
    # (3) direct bitset solve: exists (lam, mu) : lam HX = x0, lam CD + mu HZ = 0 ?
    cols = []
    HX, HZ = p.HX, p.HZ
    # unknown (lam, mu) length 2*dim; equations: for each col j: lam HX[:,j] = x0[j];
    # for each col j: lam CD[:,j] + mu HZ[:,j] = 0
    eqs = []
    tgt = []
    dim = p.dim
    for j in range(HX.shape[1]):
        row = HX[:, j]
        eqs.append({"lam": int.from_bytes(np.packbits(row[::-1], bitorder="little").tobytes(),
                                          "little"), "rhs": int(x0[j])})
    for j in range(CD.shape[1]):
        lam_part = CD[:, j]
        mu_part = HZ[:, j]
        eqs.append({"lam": int.from_bytes(np.packbits(lam_part[::-1], bitorder="little").tobytes(),
                                          "little"),
                    "mu": int.from_bytes(np.packbits(mu_part[::-1], bitorder="little").tobytes(),
                                          "little"),
                    "rhs": 0})
    # solve via Gaussian elimination over python ints: unknowns = 2*dim bits (lam: low dim, mu: high dim)
    rows_i = []
    for e in eqs:
        v = e["lam"] | (e.get("mu", 0) << dim)
        rows_i.append({"v": v, "rhs": e["rhs"]})
    piv = {}
    # eliminate
    basis = {}
    for e in rows_i:
        v, b = e["v"], e["rhs"]
        while v:
            p_ = v.bit_length() - 1
            if p_ in basis:
                v ^= basis[p_][0]
                b ^= basis[p_][1]
            else:
                basis[p_] = (v, b)
                break
        else:
            if b:
                basis["contradiction"] = True
    # consistency: any all-zero lhs with rhs=1 means unsolvable => x0 not in S_X(Q)
    c["path3_x_notin_SXQ"] = bool(basis.get("contradiction", False))
    # alternative path3 check: x0 solves if y lam HX[:j]=x0[j] has sol; combined
    c["demoted_consensus"] = (c["path1_z_notin_SZplusDelta"] and c["path2_x_notin_SXQ"]
                              and c["path3_x_notin_SXQ"])
    c["below_certified_bound"] = (d_lb is not None and c["weight_x"] < d_lb)
    out["VERIFIED"] = bool(
        c["M_symmetric"] and c["lam_M_zero"] and c["x_in_XcenQ(kerHZ,kerCD)"]
        and c["x_in_rowspace_HX"] and c["demoted_consensus"] and c["below_certified_bound"])
    out["rho_X"] = int(rank_np(np.vstack([p.HZ, CD])) - p.rHZ)
    out["delta_bar"] = sp["delta_bar"]
    out["k_Q"] = int(p.kP - sp["delta_bar"])
    return out


def main():
    files = [a for a in sys.argv[1:] if a.endswith(".json")]
    if files:
        wits = []
        for f in files:
            d = json.loads(Path(f).read_text())
            d.setdefault("labels", d.get("labels", []))
            d.setdefault("verified", True)
            wits.append(d)
    else:
        data = json.loads((HERE.parent / "exp044_n180_hunt.json").read_text())
        wits = data.get("decrease_witnesses", [])
    idxs = [int(a) for a in sys.argv[1:] if a.isdigit()] or list(range(len(wits)))
    seen_fp = set()
    results = []
    for i in idxs:
        if i >= len(wits):
            continue
        w = wits[i]
        if not w.get("verified"):
            continue
        r = reproduce(w, i)
        results.append(r)
        if r["VERIFIED"] and r["parent_fingerprint"] not in seen_fp:
            seen_fp.add(r["parent_fingerprint"])
            print("REPRODUCED", r["parent_fingerprint"], r["labels"][:2],
                  "w_x=", r["checks"]["weight_x"], "w_z=", r["checks"]["weight_z"],
                  "dZ=", r["d_Z_exact"], "M0=", r["checks"]["M_is_zero"],
                  "kQ=", r["k_Q"], "rhoX=", r["rho_X"], flush=True)
    ok = sum(1 for r in results if r["VERIFIED"])
    print(json.dumps({"reproduced": ok, "of": len(results),
                      "distinct_parents": len(seen_fp)}))


if __name__ == "__main__":
    main()
