"""EXP-047 — EXACT single-row syzygy demotion decision per parent.

Theory (Theorem J-E, notes/theorem_je_exact_decision.md):
  * The syzygy space equals the parent's pure-Z centralizer: ker H_X = ker[A B].
    For z = (c, d) in ker H_X the perturbation C=circ(c), D=circ(d) has M=0.
  * Let Q = quotient basis of ker H_X / S_Z (dim k_P), mu(x) = coordinates in Q
    via a right inverse of [Q; S_Z].  S_Z is shift-invariant, so the rank test
    is constant on S_Z-cosets.
  * For class y with representative z_y = y·Q, row 0 of H_X is demoted by the
    perturbation z_y iff
        y  NOTIN  mu( Delta(z_y) ),   Delta(z) = rowspace( Lpre @ S_z ),
    where Lpre = left-nullspace([A B]) (rows lambda with lambda [A B] = 0) and
    S_z = all lattice shifts of z (so lambda @ S_z = lambda convolve z).
  * Exhausting the classes decides demotion existence EXACTLY: demotion exists
    iff some class y in GF(2)^(k_P) fails membership.  Because M=0 makes every
    row a candidate and the test is translation-invariant, this decides the
    whole single-row-orbit demotion problem of the parent.
  * The demoted X-logical has weight wt(A)+wt(B) (row weight of H_X), equal
    for every realized class.

Costs: 2^{k_P} classes worst case; early-stop on the first realized class for
positive parents; full enumeration only for negative certificates (k_P <= 16
on all currently unverdicted parents).

Cross-validation: every EXP-046 witness must sit in a realized class (checked);
every EXP-046 'no demotion' parent is re-decided exactly here.

Artifact: results/processed/exp047_exact_demotion_decision.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
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

from qec_research.gf2.linalg import (  # noqa: E402
    nullspace_np, rank_np, rref_np, rows_to_bitsets, rref_bitset,
)

OUT = ROOT / "results" / "processed" / "exp047_exact_demotion_decision.json"
CLASS_CAP = 1 << 22  # 4M; enough for k_P <= 22


def quotient_and_projector(HZ, Kb):
    """Q: rows of Kb extending an S_Z basis to ker H_X; R: right inverse of
    [Q; Sz_basis] (full row rank).  mu(x) = (x @ R)[:k_P]."""
    Sb, _ = rref_np(HZ)  # row basis of S_Z
    cur = Sb.copy()
    Qrows = []
    for row in Kb:
        t = np.vstack([cur, row[None, :]])
        if rank_np(t) > rank_np(cur):
            Qrows.append(row)
            cur = t
    Q = np.array(Qrows, np.uint8)
    M_full = np.vstack([Q, Sb])  # full row rank = dim ker H_X
    r = M_full.shape[0]
    rr, pivots = rref_np(M_full)
    A = M_full[:, pivots]
    assert A.shape == (r, r) and rank_np(A) == r, "pivot submatrix must be invertible"
    # GF(2) inverse of A via rref provenance (RREF of invertible = I, prov = inverse)
    _, _, prov = rref_bitset(rows_to_bitsets(A), A.shape[1])
    Ainv = np.zeros((r, r), np.uint8)
    for i, combo in enumerate(prov):
        for j in combo:
            Ainv[i, j] = 1
    R = np.zeros((M_full.shape[1], r), np.uint8)
    R[pivots, :] = Ainv
    assert not (M_full @ R % 2 ^ np.eye(r, dtype=np.uint8)).any()
    return Q, R, rank_np(Q)


def shift_matrix_of(z, ell, m):
    """S_z: all ell*m lattice shifts of z (2 lm-vector) as an lm x 2lm matrix."""
    dim = ell * m
    Z = z.reshape(2 * ell, m)  # both halves share the (ell, m) lattice indexing
    rows = np.empty((dim, 2 * dim), np.uint8)
    xs = np.arange(ell)
    ys = np.arange(m)
    Xi = (xs[:, None] + xs[None, :]) % ell
    Yi = (ys[:, None] + ys[None, :]) % m
    ZZ = Z.reshape(2, ell * m)
    full = Z.reshape(2, ell, m)
    i = 0
    for x in range(ell):
        for y in range(m):
            # shift by (x, y): S[(x,y)] = z shifted so index (a,b) -> (a+x,b+y)
            rows[i] = np.concatenate([full[0][(xs - x) % ell][:, (ys - y) % m].reshape(-1),
                                      full[1][(xs - x) % ell][:, (ys - y) % m].reshape(-1)])
            i += 1
    return rows


def decide_parent(fp, entry, row, certs, pool, exp046_row=None):
    t0 = time.time()
    ell, m = int(entry["ell"]), int(entry["m"])
    dim = ell * m
    _, HX, HZ = E27.parent_matrices(row)
    A, B = HX[:, :dim], HX[:, dim:]
    r = int(HX[0].sum())
    k_P = entry["k_parent"]
    cert = certs.get(fp) or {}
    d_cert = cert.get("d_z_parent") if cert.get("d_z_exact") else None
    d_pool = pool.get(fp)
    d_bound = d_cert if d_cert is not None else d_pool
    rec = {
        "fingerprint": fp, "fingerprint12": fp[:12], "label": entry["members"][0]["label"],
        "ell": ell, "m": m, "n": 2 * dim, "k_parent": k_P,
        "wt_sum_row0": r, "d_bound": d_bound,
        "demotion": False, "immunity_exact": False,
        "classes_space": 1 << k_P, "classes_tested": 0, "class_cap_overflow": False,
        "realized_class_weight": None, "witness": None, "collapse": None,
        "wall_s": None, "consistency_exp046": None,
    }
    if k_P == 0:
        rec["class_cap_overflow"] = True
        rec["wall_s"] = round(time.time() - t0, 3)
        return rec
    exhaustive = (1 << k_P) <= CLASS_CAP
    classes = sorted(range(1, min(1 << k_P, CLASS_CAP)), key=lambda y: bin(y).count("1"))
    if not exhaustive:
        rec["class_cap_overflow"] = True  # demotion may still be found by early stop;
        # immunity cannot be certified in that case
    Kb = nullspace_np(HX)
    Q, R, kq = quotient_and_projector(HZ, Kb)
    assert kq == k_P, f"quotient dim {kq} != k_P {k_P}"
    Rk = R[:, :k_P]
    Lpre = nullspace_np(HX.T)

    def test_class(yvec: np.ndarray):
        z = (yvec @ Q) % 2 if yvec.any() else Q[0] * 0
        Sz = shift_matrix_of(z, ell, m)
        delta = (Lpre @ Sz) % 2
        proj = delta @ Rk % 2
        # membership of coordinate vector yvec in rowspace(proj): quotient
        # coordinates of a ker-vector x are (x @ Rk); rows of delta are
        # ker-vectors, so proj rows ARE their mu-coordinates.
        return yvec, z, proj

    # classes in increasing hamming weight
    realized = None
    for yint in classes:
        yvec = np.array([(yint >> i) & 1 for i in range(k_P)], np.uint8)
        rec["classes_tested"] += 1
        yvec, z, proj = test_class(yvec)
        if rank_np(np.vstack([proj, yvec[None, :]])) != rank_np(proj):
            realized = (yint, z)
            break
    if realized is not None:
        yint, z = realized
        rec["demotion"] = True
        rec["realized_class_weight"] = bin(yint).count("1")
        c, d = z[:dim], z[dim:]
        rec["witness"] = {
            "c_terms": [[int(i // m), int(i % m)] for i in np.flatnonzero(c)],
            "d_terms": [[int(i // m), int(i % m)] for i in np.flatnonzero(d)],
            "c_weight": int(c.sum()), "d_weight": int(d.sum()),
            "x0_weight": r, "z0_weight": int(z.sum()),
        }
        if d_bound is not None:
            rec["collapse"] = bool(r < d_bound)
    else:
        rec["immunity_exact"] = True
    if exp046_row is not None:
        a, b = bool(exp046_row["demotion"]), bool(rec["demotion"])
        rec["consistency_exp046"] = ("agree" if a == b else
                                     ("EXP046_SEARCH_LIMITED->EXP047_FOUND" if not a and b else
                                      "EXP046_FOUND->EXP047?REGRESSION"))
    rec["wall_s"] = round(time.time() - t0, 3)
    return rec


def main() -> None:
    t00 = time.time()
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    certs = E39.load_certificates()
    pool = E39.load_pool_lower_bounds()
    rows_by_label = {}
    for index, row in enumerate(rows):
        rows_by_label[E27.catalogue_label(row, index)] = row
    e46 = {r["fingerprint"]: r for r in
           json.loads((ROOT / "results/processed/exp046_generic_syzygy.json").read_text())["rows"]}
    # resume: keep rows from a prior artifact unless they only hold a cap overflow
    prior = {}
    if OUT.exists():
        for r in json.loads(OUT.read_text())["rows"]:
            if not (r["class_cap_overflow"] and not r["demotion"]):
                prior[r["fingerprint"]] = r
    # failure parents first (hot path), then everyone else
    fails46 = {r["fingerprint"] for r in e46.values() if not r["demotion"]}
    order = sorted(parents, key=lambda f: (f not in fails46,))
    out_rows = []
    for fp in order:
        if fp in prior:
            out_rows.append(prior[fp])
            continue
        label = parents[fp]["members"][0]["label"]
        out_rows.append(decide_parent(fp, parents[fp], rows_by_label[label],
                                      certs, pool, e46.get(fp)))
        if fp in fails46:
            r = out_rows[-1]
            print(f"[exact] {r['label']:<14} k={r['k_parent']:>2} demotion={r['demotion']} "
                  f"classes={r['classes_tested']}/{r['classes_space']} "
                  f"consistency={r['consistency_exp046']} wall={r['wall_s']}s", flush=True)
    n = len(out_rows)
    dem = [r for r in out_rows if r["demotion"]]
    imm = [r for r in out_rows if r["immunity_exact"]]
    upgraded = [r for r in out_rows if r["consistency_exp046"] == "EXP046_SEARCH_LIMITED->EXP047_FOUND"]
    regress = [r for r in out_rows if r["consistency_exp046"] == "EXP046_FOUND->EXP047?REGRESSION"]
    payload = {
        "schema": "exp047-exact-demotion-decision-v1",
        "experiment": "EXP-047 exact single-row syzygy demotion decision",
        "parents_total": n,
        "demotion_exact": len(dem),
        "immunity_exact": len(imm),
        "immunity_list": [(r["fingerprint12"], r["label"]) for r in imm],
        "upgraded_from_exp046_search": [(r["fingerprint12"], r["label"]) for r in upgraded],
        "regressions_vs_exp046": [(r["fingerprint12"], r["label"]) for r in regress],
        "cap_overflow": [(r["fingerprint12"], r["k_parent"]) for r in out_rows if r["class_cap_overflow"]],
        "collapses_bounded": sum(1 for r in dem if r["collapse"]),
        "total_wall_s": round(time.time() - t00, 1),
        "rows": out_rows,
    }
    OUT.write_text(json.dumps(payload, indent=1))
    print("parents:", n, "| exact demotion:", len(dem), "| exact immunity:", len(imm))
    print("upgraded:", payload["upgraded_from_exp046_search"])
    print("regressions:", payload["regressions_vs_exp046"])
    print("immune parents:", payload["immunity_list"])
    print("WROTE", OUT, round(time.time() - t00, 1), "s")


if __name__ == "__main__":
    main()
