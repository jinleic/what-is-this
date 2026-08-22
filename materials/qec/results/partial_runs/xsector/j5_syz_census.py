"""Syzygy census across n=144/180 certified parents.

Per parent: does the trivial-excluded M=0 family demote SINGLES (lam=e_i, weight
= row weight), PAIRS, and does the D=0 subfamily (c,0) already demote singles?
Records the minimum singles-light window, and delta_bar stats of the demoting
syzygy elements.  Output: j5_syz_census.json.
"""
from __future__ import annotations

import importlib.util
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

_spec = importlib.util.spec_from_file_location(
    "exp027_delta_audit", ROOT / "experiments" / "exp027_delta_audit.py")
E27 = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = E27; _spec.loader.exec_module(E27)
_spec39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(_spec39); sys.modules[_spec39.name] = E39; _spec39.loader.exec_module(E39)

from qec_research.gf2.linalg import rank_np  # noqa: E402
import j5core  # noqa: E402
from exp044_n180_hunt import light_scan_bits  # noqa: E402


def main():
    t0 = time.monotonic()
    rng = np.random.default_rng(0xC3A5)
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    certs = E39.load_certificates()
    pool = E39.load_pool_lower_bounds()
    out = {"schema": "j5-syz-census-v1", "parents": []}
    order = sorted(parents.items(), key=lambda kv: kv[1]["n"])
    for fp, entry in order:
        n = entry["n"]
        if n not in (36, 72, 108, 144, 180):
            continue
        ell, m = entry["ell"], entry["m"]
        cert = certs.get(fp)
        if cert and cert.get("d_z_exact"):
            d_lb, src = int(cert["d_z_parent"]), "exp039"
        elif fp in pool:
            d_lb, src = int(pool[fp]), "pool"
        else:
            d_lb, src = None, None
        row0 = None
        for r in rows:
            if int(r["ell"]) == ell and int(r["m"]) == m:
                _, HXr, HZr = E27.parent_matrices(r)
                if E27.matrix_fingerprint(HXr, HZr) == fp:
                    row0 = r
                    break
        lat = j5core.Lattice(ell, m)
        p = j5core.ParentData(lat, row0["A_terms"], row0["B_terms"])
        dim = p.dim
        lights = light_scan_bits(p.HX, dim, d_lb if d_lb else 16)[:6000]
        singles = [l0 for l0 in lights if len(l0["lam"]) == 1]
        pairs = [l0 for l0 in lights if len(l0["lam"]) == 2]
        trips = [l0 for l0 in lights if len(l0["lam"]) >= 3]
        rec = {"fingerprint": fp[:16], "labels": [mm["label"] for mm in entry["members"][:3]],
               "ell": ell, "m": m, "k_P": int(p.kP), "sigma": int(p.sigma),
               "rankA": int(rank_np(p.A)), "rankB": int(rank_np(p.B)),
               "row_weight": int(p.HX[0].sum()), "d_Z_lower": d_lb, "d_Z_source": src,
               "n_singles_light": len(singles), "n_pairs_light": len(pairs),
               "m0_sampled": 0,
               "demote_singles": 0, "demote_pairs_only": 0, "demote_none": 0,
               "d0_sampled": 0, "d0_demote_singles": 0,
               "min_w_demoted": None, "delta_bars_demoting": []}
        S = p.syzygy
        samples = [S[r] for r in range(S.shape[0])]
        for _ in range(24 if S.shape[0] else 0):
            sel = rng.integers(0, 2, S.shape[0]).astype(np.uint8)
            if sel.any():
                samples.append((sel @ S) % 2)
        for cd in samples:
            C, D = j5core.split_CD(p, cd)
            rec["m0_sampled"] += 1
            demoted, sp = j5core.batch_demotion_test(p, C, D,
                                                     [{"w": l0["w"], "lam": l0["lam"]}
                                                      for l0 in lights])
            if not demoted:
                rec["demote_none"] += 1
                continue
            if any(len(d0["lam"]) == 1 for d0 in demoted):
                rec["demote_singles"] += 1
                rec["delta_bars_demoting"].append(sp["delta_bar"])
            else:
                rec["demote_pairs_only"] += 1
            mw = min(d0["w"] for d0 in demoted)
            if rec["min_w_demoted"] is None or mw < rec["min_w_demoted"]:
                rec["min_w_demoted"] = mw
        # D=0 subfamily: c in ker(A as column map) i.e. a c~ = 0
        from qec_research.gf2.linalg import nullspace_np
        annA = nullspace_np(p.A)          # rows c with A c^T = 0  <=>  a c~ = 0?? check below
        # NOTE: M = A C^T = 0 <=> JA c?? row0 map: JA c = 0 <=> A c = 0. so c in right-ker A.
        for i in range(min(24, annA.shape[0])):
            cvec = annA[i]
            if not cvec.any():
                continue
            # random combos too
        d0s = [annA[i] for i in range(annA.shape[0])][:24]
        for _ in range(24 if annA.shape[0] else 0):
            sel = rng.integers(0, 2, annA.shape[0]).astype(np.uint8)
            if sel.any():
                d0s.append((sel @ annA) % 2)
        for cvec in d0s:
            C = lat.sparse(cvec)
            D = np.zeros_like(C)
            M = p.M_of(C, D)
            if M.any():
                continue  # not an annihilator element after all
            rec["d0_sampled"] += 1
            demoted, sp = j5core.batch_demotion_test(
                p, C, D, [{"w": l0["w"], "lam": l0["lam"]} for l0 in singles + pairs])
            if any(len(d0["lam"]) == 1 for d0 in demoted):
                rec["d0_demote_singles"] += 1
        rec["delta_bars_demoting"] = sorted(set(rec["delta_bars_demoting"]))
        out["parents"].append(rec)
        print(f"{fp[:12]} {rec['labels'][0]:14s} kP={rec['k_P']:2d} sig={rec['sigma']:3d} "
              f"dZ={d_lb} singles={rec['demote_singles']}/{rec['m0_sampled']} "
              f"D0={rec['d0_demote_singles']}/{rec['d0_sampled']} "
              f"dBars={rec['delta_bars_demoting'][:5]} t={time.monotonic()-t0:.0f}s",
              flush=True)
    out["wall_seconds"] = round(time.monotonic() - t0, 1)
    (HERE.parent / "j5_syz_census.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({"parents": len(out["parents"]), "seconds": out["wall_seconds"]}))


if __name__ == "__main__":
    main()
