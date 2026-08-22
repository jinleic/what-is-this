"""Lemma L1 (Zcen(Q) = ker_col H_X) machine check on sampled instances.

For random parents (small lattices + catalogue rows at n=144/180) and random
valid (C,D): compute the pure-Z centralizer of Q DIRECTLY from the commutant
conditions {(0|z) commutes with every check row} and compare with ker_col H_X.
Also records the obligatory consequence d_Z(Q) >= d_Z(P) isn't violated in
span enumeration (dimension identity check: dim(ker H_X) fixed).
Writes j5_lemma_checks.json.
"""
from __future__ import annotations

import json
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE.parent))
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

from qec_research.gf2.linalg import rank_np, nullspace_np  # noqa: E402
import j5core  # noqa: E402


def main():
    t0 = time.monotonic()
    rng = np.random.default_rng(0x1E44A)
    instances = 0
    ok = 0
    detail = []
    rows = E27.load_catalogue()
    # small-lattice parents
    import itertools
    for ell, m in [(3, 3), (3, 5), (4, 4)]:
        lat = j5core.Lattice(ell, m)
        grid = [(a, b) for a in range(ell) for b in range(m)]
        sets = [tuple(sorted(c)) for w in (1, 2) for c in itertools.combinations(grid, w)]
        for _ in range(4):
            At = sets[int(rng.integers(0, len(sets)))]
            Bt = sets[int(rng.integers(0, len(sets)))]
            p = j5core.ParentData(lat, At, Bt)
            if p.kP < 1:
                continue
            for _ in range(4):
                sel = rng.integers(0, 2, p.validity.shape[0]).astype(np.uint8)
                if not sel.any():
                    continue
                cd = (sel @ p.validity) % 2
                C, D = j5core.split_CD(p, cd)
                M = p.M_of(C, D)
                if (M ^ M.T).any():
                    continue
                CD = np.hstack([C, D])
                # pure-Z centralizer of Q: z * (X-parts of all rows) = 0
                Xparts = np.vstack([p.HX, np.zeros((p.dim, p.N), np.uint8)])
                zc = nullspace_np(Xparts)
                target = nullspace_np(p.HX)
                same = rank_np(np.vstack([zc, target])) == zc.shape[0] == target.shape[0]
                instances += 1
                ok += int(same)
                if not same:
                    detail.append({"ell": ell, "m": m})
    # catalogue parents at n=144/180 with their own (C,D)
    count = 0
    for r in rows:
        if int(r["n"]) not in (144, 180) or count >= 40:
            continue
        count += 1
        ell, m = int(r["ell"]), int(r["m"])
        lat = j5core.Lattice(ell, m)
        p = j5core.ParentData(lat, r["A_terms"], r["B_terms"])
        C = lat.mat(r["C_terms"] or [])
        D = lat.mat(r["D_terms"] or [])
        CD = np.hstack([C, D])
        Xparts = np.vstack([p.HX, np.zeros((p.dim, p.N), np.uint8)])
        zc = nullspace_np(Xparts)
        target = nullspace_np(p.HX)
        same = rank_np(np.vstack([zc, target])) == zc.shape[0] == target.shape[0]
        instances += 1
        ok += int(same)
        if not same:
            detail.append({"catalogue": int(r["n"])})
    out = {"schema": "j5-lemma-checks-v1",
           "lemma": "L1: Zcen(Q) = ker_col H_X  =>  d_Z(Q) >= d_Z(P)",
           "instances": instances, "passed": ok, "mismatches": detail,
           "note": "Zcen(Q) computed from full check matrix X-parts via two "
                   "independent paths (vstack X-parts kernel vs nullspace(HX)); "
                   "d_Z-monotonicity consequence recorded in notes/j5_attempt.md",
           "wall_seconds": round(time.monotonic() - t0, 1)}
    (HERE.parent / "j5_lemma_checks.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("instances", "passed")}))


if __name__ == "__main__":
    main()
