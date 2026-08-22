"""EXP-044: does any valid perturbation of any catalogue parent exhibit d_X(Q) < d_X(P)?

By Theorem J.4, d_X(Q) < d_X(P) iff some demoted parent X-stabilizer has weight
< d_X(P) ("pure-X collapse").  Demotion set = span(kerM)*[A B] minus S_X(Q), with
kerM = left-null(M), M = AC^T + BD^T (symmetric for valid (C,D)).

Valid (C,D) form the linear space V = ker(nu), nu(C,D) = skew(AC^T + BD^T).
Measured across all 202 parents: dim V ranges 29 (n=36) to 282 (n=360), so full
enumeration of V is out of reach everywhere (2^29 at the smallest parent); we
sample uniform elements of V per parent (dense + small basis sums), add the
catalogue's own perturbations, and for parents with ell*m <= 72 additionally
rejection-sample light (C,D) of weight <= 8, then run the exact demotion check
(full 2^dk span enumeration, dk = dim kerM capped at 18) against a certified
d_Z(P) lower bound.  GF(2)/numpy only, one thread, no SAT.  Hard deadline.
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

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

from qec_research.codes.bicycle import monomial_matrix, poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np, matmul as gf2matmul  # noqa: E402

DEADLINE = time.monotonic() + 26 * 60
OUT = ROOT / "results/processed/exp044_xsector_collapse_hunt.json"
KERM_CAP = 18
PROBE_SAMPLES = 4096


class MonoTables:
    def __init__(self, ell, m):
        self.ell, self.m, self.dim = ell, m, ell * m
        self.X = np.zeros((self.dim, self.dim, self.dim), dtype=np.uint8)
        for k, (a, b) in enumerate(itertools.product(range(ell), range(m))):
            self.X[k] = monomial_matrix(ell, m, a, b)

    def mat(self, termset):
        out = np.zeros((self.dim, self.dim), np.uint8)
        for a, b in termset:
            out ^= self.X[(a % self.ell) * self.m + b % self.m]
        return out

    def sparse(self, coeff_vec):
        nz = np.flatnonzero(coeff_vec)
        out = np.zeros((self.dim, self.dim), np.uint8)
        for k in nz:
            out ^= self.X[k]
        return out


def validity_kernel(mt: MonoTables, A, B):
    dim = mt.dim
    ii, jj = np.triu_indices(dim, k=1)
    # stack all products A X_k^T and B X_k^T in two big matmuls
    XT = np.moveaxis(mt.X, -1, 1).reshape(dim, dim * dim)  # same as X.T per block
    colsA = np.empty((2 * dim, ii.size), dtype=np.uint8)
    PA = (A @ mt.X.reshape(dim, dim * dim)) % 2  # wrong axis; do per-k loop instead
    # per-k loop (cheap enough: 2*dim matmuls of dim x dim)
    for k in range(dim):
        P = (A @ mt.X[k].T) % 2
        colsA[k] = (P ^ P.T)[ii, jj]
    for k in range(dim):
        P = (B @ mt.X[k].T) % 2
        colsA[dim + k] = (P ^ P.T)[ii, jj]
    N = colsA.T
    return nullspace_np(N)


def demotion_classify(dim, N, HX, HZ, L, CD, M, D_lb, rng):
    """Returns (cls, best_w, witness_x_or_None)."""
    kerM = nullspace_np(M.T)
    dk = kerM.shape[0]
    if dk == 0:
        return "no_demotion", None, None
    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, N), np.uint8)
    SZD = np.vstack([HZ, Delta]) if Delta.shape[0] else HZ
    Wperp = nullspace_np(SZD)
    Lam = nullspace_np(((CD @ Wperp.T) % 2).T)
    SXQ = (Lam @ HX) % 2 if Lam.shape[0] else np.zeros((0, N), np.uint8)
    ns = nullspace_np(SXQ) if SXQ.shape[0] else None
    best, bestx = None, None
    if dk <= KERM_CAP:
        for start in range(1, 1 << dk, 1 << 16):
            idx = np.arange(start, min(start + (1 << 16), 1 << dk), dtype=np.int64)
            V = np.zeros((idx.size, dim), dtype=np.uint8)
            for j in range(dk):
                V ^= (((idx >> j) & 1)[:, None] & kerM[j]).astype(np.uint8)
            X = (V @ HX) % 2
            w = X.sum(axis=1)
            keep = w > 0
            if ns is not None and ns.shape[0]:
                keep &= ((X @ ns.T) % 2).sum(axis=1) != 0
            if keep.any():
                ww = w[keep]
                r = int(ww.argmin())
                b = int(ww[r])
                if best is None or b < best:
                    best, bestx = b, X[keep][r]
    else:
        for _ in range(PROBE_SAMPLES):
            sel = rng.integers(0, 2, dk).astype(np.uint8)
            lam = (sel @ kerM) % 2
            x = (lam @ HX) % 2
            w = int(x.sum())
            if w == 0 or (best is not None and w >= best):
                continue
            if ns is not None and ns.shape[0] and ((x @ ns.T) % 2).sum() == 0:
                continue
            best, bestx = w, x
    if best is None:
        return "no_demotion", None, None
    if D_lb is not None and best < D_lb:
        return "decrease_witness", best, bestx
    if dk > KERM_CAP:
        return "unresolved_capped", best, None
    if D_lb is None:
        return "safe_exact_no_bound", best, None
    return "safe_exact", best, None


def verify_witness(ell, m, A, B, C, D, x, D_lb):
    HX = np.hstack([A, B])
    HZ = np.hstack([B.T, A.T])
    CD = np.hstack([C, D])
    M = (gf2matmul(A, C.T) ^ gf2matmul(B, D.T)).astype(np.uint8)
    return all([
        not ((M ^ M.T).any()),
        not ((x @ HZ.T) % 2).any(),
        not ((x @ CD.T) % 2).any(),
        (((x @ nullspace_np(HX).T) % 2).sum()) == 0,
        int(x.sum()) < D_lb,
    ])


def main():
    t0 = time.monotonic()
    rng = np.random.default_rng(0xC011A9)
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    certs = E39.load_certificates()
    pool = E39.load_pool_lower_bounds()

    cat_cds: dict = {}
    for row in rows:
        _, HXr, HZr = E27.parent_matrices(row)
        fp = E27.matrix_fingerprint(HXr, HZr)
        cat_cds.setdefault(fp, []).append(
            (tuple(map(tuple, row.get("C_terms") or [])), tuple(map(tuple, row.get("D_terms") or []))))

    # small parents first
    order = sorted(parents.items(), key=lambda kv: kv[1]["n"])
    per_parent = []
    totals = {"instances": 0, "no_demotion": 0, "safe_exact": 0, "safe_exact_no_bound": 0,
              "unresolved_capped": 0, "decrease_witness": 0}
    witnesses = []
    global_budget_max = 9_000
    for fp, entry in order:
        if time.monotonic() > DEADLINE or totals["instances"] >= global_budget_max:
            break
        tp = time.monotonic()
        ell, m = entry["ell"], entry["m"]
        row0 = None
        for r in rows:
            if int(r["ell"]) == ell and int(r["m"]) == m:
                _, HXr, HZr = E27.parent_matrices(r)
                if E27.matrix_fingerprint(HXr, HZr) == fp:
                    row0 = r
                    break
        if row0 is None:
            continue
        mt = MonoTables(ell, m)
        dim, N = ell * m, 2 * ell * m
        A = poly_matrix(ell, m, [tuple(t) for t in row0["A_terms"]])
        B = poly_matrix(ell, m, [tuple(t) for t in row0["B_terms"]])
        HX, HZ = np.hstack([A, B]), np.hstack([B.T, A.T])
        kP = N - rank_np(HX) - rank_np(HZ)
        if kP < 2:
            continue
        L = nullspace_np(HX.T)
        cert = certs.get(fp)
        D_lb = int(cert["d_z_parent"]) if cert and cert.get("d_z_parent") else (int(pool[fp]) if fp in pool else None)
        dZ_src = "exp039-cert(exact)" if cert and cert.get("d_z_exact") else (
            "exp039-cert(bound)" if cert else ("exp037-pool(bound)" if fp in pool else None))
        V = validity_kernel(mt, A, B)
        dimV = V.shape[0]

        inst = []  # list of (tag, C, D)
        for Ct, Dt in cat_cds.get(fp, [])[:6]:
            C = mt.mat(Ct) if Ct else np.zeros((dim, dim), np.uint8)
            D = mt.mat(Dt) if Dt else np.zeros((dim, dim), np.uint8)
            if C.any() or D.any():
                inst.append(("catalogue", C, D))
        # sampled V elements
        d = dimV
        n_dense = 25 if dim <= 90 else 15
        n_small = 10 if dim <= 90 else 5
        for _ in range(n_dense + n_small):
            sel = rng.integers(0, 2, d).astype(np.uint8)
            if not sel.any():
                continue
            cd = (sel @ V) % 2
            inst.append(("V-sample", mt.sparse(cd[:dim]), mt.sparse(cd[dim:])))
        for _ in range(n_small):
            sel = np.zeros(d, dtype=np.uint8)
            sel[rng.choice(d, min(3, d), replace=False)] = 1
            cd = (sel @ V) % 2
            inst.append(("V-small", mt.sparse(cd[:dim]), mt.sparse(cd[dim:])))
        # light rejection sampling only for small lattices
        if dim <= 72:
            got, tried = 0, 0
            while got < 10 and tried < 800 and time.monotonic() < DEADLINE:
                tried += 1
                delta = int(rng.integers(2, 9))
                C = np.zeros((dim, dim), np.uint8)
                D = np.zeros((dim, dim), np.uint8)
                bits = rng.integers(0, dim, delta)
                for k in bits:
                    if rng.random() < 0.5:
                        C ^= mt.X[k]
                    else:
                        D ^= mt.X[k]
                if not (C.any() or D.any()):
                    continue
                M0 = (gf2matmul(A, C.T) ^ gf2matmul(B, D.T)).astype(np.uint8)
                if (M0 ^ M0.T).any():
                    continue
                inst.append(("light-reject", C, D))
                got += 1

        rec = {"fingerprint": fp[:16], "labels": [mm["label"] for mm in entry["members"][:4]],
               "ell": ell, "m": m, "n": N, "k_P": kP, "dim_V": int(dimV),
               "d_Z_lower": D_lb, "d_Z_source": dZ_src,
               "classes": {}, "min_w_dem": None, "instances": 0,
               "truncated": False}
        for tag, C, D in inst:
            if time.monotonic() > DEADLINE:
                rec["truncated"] = True
                break
            CD = np.hstack([C, D])
            M = (gf2matmul(A, C.T) ^ gf2matmul(B, D.T)).astype(np.uint8)
            if (M ^ M.T).any():
                rec["classes"]["invalid_skipped"] = rec["classes"].get("invalid_skipped", 0) + 1
                continue
            cls, w, wx = demotion_classify(dim, N, HX, HZ, L, CD, M, D_lb, rng)
            rec["instances"] += 1
            totals["instances"] += 1
            key = {"decrease_witness": "decrease_witness"}.get(cls, cls)
            rec["classes"][key] = rec["classes"].get(key, 0) + 1
            totals[key] = totals.get(key, 0) + 1
            if w is not None and (rec["min_w_dem"] is None or w < rec["min_w_dem"]):
                rec["min_w_dem"] = w
            if cls == "decrease_witness":
                ok = verify_witness(ell, m, A, B, C, D, wx, D_lb) if wx is not None else False
                wit = {"fingerprint": fp[:16], "labels": rec["labels"], "w_dem": w,
                       "d_Z_lower": D_lb, "verified": bool(ok),
                       "witness_support": sorted(int(i) for i in np.nonzero(wx)[0]) if wx is not None else None}
                witnesses.append(wit)
                print("WITNESS", json.dumps(wit), flush=True)
        rec["seconds"] = round(time.monotonic() - tp, 2)
        per_parent.append(rec)
        if len(per_parent) % 20 == 0:
            print(f"progress {len(per_parent)}/202 elapsed {time.monotonic()-t0:.0f}s "
                  f"instances {totals['instances']}", flush=True)

    payload = {
        "schema": "exp044-xsector-collapse-hunt-v1",
        "generated_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "question": "does d_X(Q) < d_X(P) occur for any valid perturbation of any catalogue parent?",
        "method": {
            "valid_space": "V = ker(nu(C,D)=skew(AC^T+BD^T)); dim V recorded per parent (29..282)",
            "sampling": "catalogue-own (C,D) <=6 + uniform dense V elements + small basis sums + light rejection (dim<=72 only)",
            "check": "d_X(Q) < d_X(P) iff light demotion (Thm J.4); demotion = span(kerM)[A B] minus S_X(Q); exact 2^dk enumeration (dk<=18) else uniform probe",
            "decrease_threshold": "certified d_Z(P) lower bound (exp039 certs, else exp037 pool); d_X(P) = d_Z(P) by Thm J.0",
            "kerm_cap": KERM_CAP,
            "probe_samples": PROBE_SAMPLES,
        },
        "parents_covered": len(per_parent),
        "wall_seconds": round(time.monotonic() - t0, 1),
        "totals": totals,
        "decrease_witnesses": witnesses,
        "per_parent": per_parent,
    }
    if witnesses:
        payload["verdict"] = "X_COLLAPSE_WITNESSES_FOUND"
    elif totals["unresolved_capped"] == 0:
        payload["verdict"] = "X_COLLAPSE_NEVER_OBSERVED"
    else:
        payload["verdict"] = "X_COLLAPSE_NEVER_OBSERVED_WITHIN_COMPUTED_BUDGET"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=1))
    print(json.dumps({"verdict": payload["verdict"], "totals": totals,
                      "parents": len(per_parent), "seconds": payload["wall_seconds"]}))


if __name__ == "__main__":
    main()
