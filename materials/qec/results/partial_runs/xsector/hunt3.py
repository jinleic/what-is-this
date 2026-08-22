"""Round 3: EXHAUSTIVE perturbation space on the smallest lattices.

For lattices with ell*m in {4, 6}: every (C, D) pair of nonempty subsets up to
simultaneous translation, against all parents with wt(A), wt(B) <= 3 (plus
singletons).  For ell*m = 8: full subset space too.  For ell*m = 9: |C|+|D| <= 7.
A complete computed statement: no valid perturbation on these lattices
exhibits d_X(Q) < d_X(P).  GF(2)/numpy only, one thread, no SAT.
"""
from __future__ import annotations

import itertools
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

from qec_research.codes.bicycle import poly_matrix  # noqa: E402
from qec_research.gf2.linalg import rank_np, nullspace_np, matmul as gf2matmul  # noqa: E402
from hunt_enriched import MonoCache, dX_support, span_outside, menu  # noqa: E402


def all_subset_perts(ell, m, max_total):
    grid = [(a, b) for a in range(ell) for b in range(m)]
    seen, out = set(), []
    for tot in range(1, max_total + 1):
        for split in range(tot + 1):
            if split > len(grid) or tot - split > len(grid):
                continue
            for Cs in itertools.combinations(grid, split):
                for Ds in itertools.combinations(grid, tot - split):
                    allm = sorted(Cs + Ds)
                    t = allm[0]
                    def shift(s):
                        return tuple(sorted(((a - t[0]) % ell, (b - t[1]) % m) for a, b in s))
                    key = (shift(Cs), shift(Ds))
                    if key not in seen:
                        seen.add(key)
                        out.append(key)
    return out


def parent_terms(ell, m, max_wt):
    grid = [(a, b) for a in range(ell) for b in range(m)]
    S = []
    for w in range(1, max_wt + 1):
        for c in itertools.combinations(grid, w):
            S.append(tuple(c))
    return S


def main():
    spec = [((2, 2), 99), ((2, 3), 99), ((3, 2), 99), ((2, 4), 99), ((4, 2), 99), ((3, 3), 7)]
    stats = {"parents": 0, "parents_used": 0, "cand_perts": 0, "valid_perturbations": 0,
             "with_demotion": 0, "counterexamples": 0}
    gaps = []
    best = None
    for (ell, m), max_total in spec:
        dim, N = ell * m, 2 * ell * m
        mc = MonoCache(ell, m)
        cap_total = min(max_total, 2 * dim)
        perts = all_subset_perts(ell, m, cap_total)
        print(f"lattice ({ell},{m}): {len(perts)} anchored perturbations", flush=True)
        terms = parent_terms(ell, m, 3) if dim <= 6 else sorted(
            set(menu(ell, m)) | {((a, 0),) for a in range(1, ell)}
            | {((0, b),) for b in range(1, m)} | {((0, 0),)})
        seen = set()
        for At in terms:
            for Bt in terms:
                A, B = mc.mono(At), mc.mono(Bt)
                key = (A.tobytes(), B.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                stats["parents"] += 1
                HX = np.hstack([A, B])
                HZ = np.hstack([B.T, A.T])
                kP = N - rank_np(HX) - rank_np(HZ)
                if kP < 2:
                    continue
                Zcen = nullspace_np(HX)
                dXP, _ = dX_support(HZ, Zcen, N, 8)
                eff = dXP if dXP is not None else 9
                if eff < 2:
                    continue
                stats["parents_used"] += 1
                for (Ct, Dt) in perts:
                    stats["cand_perts"] += 1
                    C, D = mc.mono(Ct), mc.mono(Dt)
                    M = (gf2matmul(A, C.T) ^ gf2matmul(B, D.T)).astype(np.uint8)
                    if (M ^ M.T).any():
                        continue
                    stats["valid_perturbations"] += 1
                    CD = np.hstack([C, D])
                    L = nullspace_np(HX.T)
                    Delta = (L @ CD) % 2 if L.shape[0] else np.zeros((0, N), np.uint8)
                    SZD = np.vstack([HZ, Delta]) if Delta.shape[0] else HZ
                    Wperp = nullspace_np(SZD)
                    G = (CD @ Wperp.T) % 2
                    Lam = nullspace_np(G.T)
                    SXQ = (Lam @ HX) % 2 if Lam.shape[0] else np.zeros((0, N), np.uint8)
                    kerM = nullspace_np(M.T)
                    mind, minx = span_outside(kerM, HX, SXQ)
                    if mind is None:
                        continue
                    stats["with_demotion"] += 1
                    gaps.append(mind - eff)
                    if mind < eff:
                        stats["counterexamples"] += 1
                        rec = {
                            "ell": ell, "m": m,
                            "A": [list(t) for t in At], "B": [list(t) for t in Bt],
                            "C": [list(t) for t in Ct], "D": [list(t) for t in Dt],
                            "k_P": int(kP), "d_X_P": int(dXP) if dXP else ">8",
                            "w_dem": int(mind),
                            "support": sorted(int(i) for i in np.nonzero(minx[0])[0]),
                        }
                        if best is None or mind - eff < best[0]:
                            best = (mind - eff, rec)
                            print("COUNTEREXAMPLE", json.dumps(rec), flush=True)
        print(f"lattice ({ell},{m}) done: {json.dumps(stats)}, min gap: "
              f"{min(gaps) if gaps else None}", flush=True)
    gaps.sort()
    print(json.dumps(stats))
    print("gap min / p1 / median:", gaps[0] if gaps else None,
          gaps[max(0, len(gaps) // 100)] if gaps else None,
          gaps[len(gaps) // 2] if gaps else None)
    out = {"stats": stats, "best": best[1] if best else None,
           "gap_min": gaps[0] if gaps else None, "n_gaps": len(gaps)}
    (ROOT / "results/partial_runs/xsector/hunt3.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
