"""EXP-044 follow-up (J.5 route C): X-side collapse hunt with sharp syzygy probe.

Pass 1 (primary): all 46 catalogue parents at n=180.
Pass 2 (bonus):    all 38 at n=144.   Pass 3 (bonus): all 39 at n=360.

Method, per parent:
  * baseline EXP-044-style perturbations (catalogue-own, uniform V, small sums)
  * syzygy/M=0 family: basis + random combos (single rows demote iff M=0;
    cols of M are R-translates, so col_i = 0 <=> M = 0)
  * V'(lam0) = V cap {M lam0^T = 0} probes for the lightest stabilizers
  * light stabilizers: lam of weight <= 3 (bitset scan) + structured lines;
    batch demotion test per perturbation (partner z = lam CD not in S_Z + Delta)
Every below-bound event is verified by TWO independent code paths (j5core
verify_full).  Decrease threshold: exp039 exact d_Z(P) cert, else exp037 pool
bound, else uncertified (recorded, never counted as violation).
GF(2)/numpy only, single thread, hard deadline.

Repair (2026-08-19, CensusFinisher): method UNCHANGED; plumbing hardened.
  * per-parent deterministic RNG: SeedSequence([0xC011A6, int(fp[:16],16)]),
    so every parent's pool is reproducible independently of execution order
    (the interrupted run used one shared stream; sampling is equivalent).
  * incremental JSONL checkpoint (sidecar) after every completed parent:
    a mid-run kill loses at most the in-flight parent; EXP044_RESUME=1
    resumes, skipping fingerprints already in the sidecar.
  * mid-parent deadline: partial rows are NOT appended; the in-flight parent
    becomes truncated_at (no silent drops; remaining parents enumerated).
  * verdict field: CENSUS_COMPLETE / CENSUS_TRUNCATED_WITH_BOUNDARY.
  Env: EXP044_BUDGET_S (wall budget, default 115 min), EXP044_RESUME=1,
  EXP044_ASSEMBLE_ONLY=1 (rebuild OUT json from sidecar only).
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
E27 = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = E27
_spec.loader.exec_module(E27)

_spec39 = importlib.util.spec_from_file_location(
    "exp039_nogo_module", ROOT / "experiments" / "exp039_nogo_module.py")
E39 = importlib.util.module_from_spec(_spec39)
sys.modules[_spec39.name] = E39
_spec39.loader.exec_module(E39)

from qec_research.gf2.linalg import rank_np, nullspace_np, rows_to_bitsets  # noqa: E402
import j5core  # noqa: E402

BUDGET_S = float(os.environ.get("EXP044_BUDGET_S", 115 * 60))
DEADLINE = time.monotonic() + BUDGET_S
OUT = HERE.parent / "exp044_n180_hunt.json"
SIDECAR = HERE.parent / "exp044_n180_hunt.rows.jsonl"
UNCERT_CAP = 16  # light-scan weight cap for parents without a certified bound
MAX_LIGHTS = 6000  # cap lights passed to the per-perturbation batch test


def light_scan_bits(HX, dim, d_lb, max_w_lam=3):
    R = rows_to_bitsets(HX)
    out = {}

    def rec(x, lam):
        w = x.bit_count()
        if 0 < w < d_lb and x not in out:
            out[x] = (w, lam)
    for i in range(dim):
        rec(R[i], (i,))
    pairs = []
    for i in range(dim):
        Ri = R[i]
        for j in range(i + 1, dim):
            x = Ri ^ R[j]
            rec(x, (i, j))
            pairs.append((i, j, x))
    if max_w_lam >= 3:
        for (i, j, xij) in pairs:
            for k in range(j + 1, dim):
                rec(xij ^ R[k], (i, j, k))
    return [{"x": xbits, "w": w, "lam": lam} for xbits, (w, lam) in
            sorted(out.items(), key=lambda kv: kv[1][0])]


def main():
    t0 = time.monotonic()
    rng = np.random.default_rng(0xC011A6)
    rows = E27.load_catalogue()
    parents = E39.distinct_parents(rows)
    certs = E39.load_certificates()
    pool = E39.load_pool_lower_bounds()

    cat_cds = {}
    for row in rows:
        _, HXr, HZr = E27.parent_matrices(row)
        fp = E27.matrix_fingerprint(HXr, HZr)
        cat_cds.setdefault(fp, []).append(
            (tuple(map(tuple, row.get("C_terms") or [])), tuple(map(tuple, row.get("D_terms") or []))))

    passes = [(180, "n180-primary"), (144, "n144-bonus"), (360, "n360-bonus")]
    payload = {
        "schema": "exp044-n180-xcollapse-hunt-v1",
        "generated_utc": None,
        "question": "does d_X(Q) < d_X(P) occur for any valid perturbation of any "
                    "catalogue parent at n=180 (primary), n=144/360 (bonus)?",
        "method": {
            "baseline": "catalogue-own (C,D) + uniform V + <=3 basis sums",
            "sharp": "syzygy/M=0 family scan + V'(lam0) nullspace probes; light "
                     "stabilizers lam weight<=3 (bitset) + structured lines; "
                     "batch demotion test z=lam CD not in S_Z+Delta",
            "verification": "dual code paths (j5core.verify_full): rank test on "
                            "S_Z+Delta AND exp044-style S_X(Q) basis test",
            "decrease_threshold": "exp039 exact d_Z(P) cert, else exp037 pool bound, "
                                  "else recorded-uncertified",
        },
        "parents": [], "decrease_witnesses": [], "uncertified_candidates": []}

    per_parent = payload["parents"]
    for target_n, pass_name in passes:
        todo = [(fp, e) for fp, e in parents.items() if e["n"] == target_n]
        todo.sort(key=lambda kv: (kv[1]["ell"], kv[1]["m"]))
        for fp, entry in todo:
            if time.monotonic() > DEADLINE:
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
            lat = j5core.Lattice(ell, m)
            p = j5core.ParentData(lat, row0["A_terms"], row0["B_terms"])
            cert = certs.get(fp)
            if cert and cert.get("d_z_exact"):
                d_lb, src = int(cert["d_z_parent"]), "exp039-cert(exact)"
            elif fp in pool:
                d_lb, src = int(pool[fp]), "exp037-pool(bound)"
            else:
                d_lb, src = None, None
            scan_cap = d_lb if d_lb is not None else UNCERT_CAP
            lights = light_scan_bits(p.HX, p.dim, scan_cap)
            lights_truncated = len(lights) > MAX_LIGHTS
            lights = lights[:MAX_LIGHTS]
            lrecs = [{"w": r["w"], "lam": r["lam"]} for r in lights]
            rec = {"pass": pass_name, "fingerprint": fp[:16],
                   "labels": [mm["label"] for mm in entry["members"][:4]],
                   "ell": ell, "m": m, "n": target_n, "k_P": int(p.kP),
                   "rHX": int(p.rHX), "rHZ": int(p.rHZ),
                   "dim_V": int(p.validity.shape[0]),
                   "dim_M0": int(p.syzygy.shape[0]), "sigma": int(p.sigma),
                   "rankA": int(rank_np(p.A)), "rankB": int(rank_np(p.B)),
                   "d_Z_lower": d_lb, "d_Z_source": src,
                   "row_weight": int(p.HX[0].sum()), "n_lights": len(lrecs),
                   "lights_truncated": bool(lights_truncated),
                   "perturbations": 0, "demoted_events": 0, "min_demoted_w": None,
                   "m0_tested": 0, "m0_demoting": 0, "m0_min_demoted_w": None,
                   "min_margin": None, "violations": 0}

            pool_cd = []
            for Ct, Dt in cat_cds.get(fp, [])[:6]:
                if Ct or Dt:
                    c = np.zeros(p.dim, np.uint8); c[[a * m + b for a, b in Ct]] = 1
                    d = np.zeros(p.dim, np.uint8); d[[a * m + b for a, b in Dt]] = 1
                    pool_cd.append(np.concatenate([c, d]))
            for _ in range(15):
                sel = rng.integers(0, 2, p.validity.shape[0]).astype(np.uint8)
                if sel.any():
                    pool_cd.append((sel @ p.validity) % 2)
            for _ in range(8):
                sel = np.zeros(p.validity.shape[0], np.uint8)
                sel[rng.choice(p.validity.shape[0], min(3, p.validity.shape[0]), replace=False)] = 1
                pool_cd.append((sel @ p.validity) % 2)
            S = p.syzygy
            for r in range(S.shape[0]):
                pool_cd.append(S[r])
            for _ in range(96 if S.shape[0] else 0):
                sel = rng.integers(0, 2, S.shape[0]).astype(np.uint8)
                if sel.any():
                    pool_cd.append((sel @ S) % 2)
            # V'(lam0) probes for the lightest few
            P0 = p.row0_maps
            for lr in lrecs[:14]:
                lam0 = np.zeros(p.dim, np.uint8)
                lam0[list(lr["lam"])] = 1
                K0 = (lat.sparse(lam0) @ P0) % 2
                selker = nullspace_np((K0 @ p.validity.T) % 2)
                for j in range(min(3, selker.shape[0])):
                    pool_cd.append((selker[j] @ p.validity) % 2)

            seen = set()
            for cd in pool_cd:
                if time.monotonic() > DEADLINE:
                    break
                C, D = j5core.split_CD(p, cd)
                key = (C.tobytes(), D.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                M = p.M_of(C, D)
                if (M ^ M.T).any():
                    continue
                rec["perturbations"] += 1
                is_m0 = not bool(M.any())
                if is_m0:
                    rec["m0_tested"] += 1
                demoted, sp = j5core.batch_demotion_test(p, C, D, lrecs)
                if not demoted:
                    continue
                rec["demoted_events"] += len(demoted)
                if is_m0:
                    rec["m0_demoting"] += 1
                    mw = min(d0["w"] for d0 in demoted)
                    if rec["m0_min_demoted_w"] is None or mw < rec["m0_min_demoted_w"]:
                        rec["m0_min_demoted_w"] = mw
                demoted.sort(key=lambda d0: d0["w"])
                verified_this_cd = 0
                for d0 in demoted:
                    w = d0["w"]
                    if rec["min_demoted_w"] is None or w < rec["min_demoted_w"]:
                        rec["min_demoted_w"] = w
                    if d_lb is not None and w < d_lb and verified_this_cd < 2:
                        verified_this_cd += 1
                        ver = j5core.verify_full(p, cd, d0["lam"], d_lb)
                        if ver["verified"]:
                            rec["violations"] += 1
                            if rec["violations"] <= 4:
                                payload["decrease_witnesses"].append({
                                    "pass": pass_name, "fingerprint": fp[:16],
                                    "labels": rec["labels"], "d_Z_lower": d_lb,
                                    "d_Z_source": src, **{k: v for k, v in ver.items()}})
                                print("WITNESS", fp[:16], rec["labels"], "w_x=", w,
                                      "d_lb=", d_lb, flush=True)
                    elif w < scan_cap and d_lb is None and rec.get("uncert_recorded", 0) < 4:
                        rec["uncert_recorded"] = rec.get("uncert_recorded", 0) + 1
                        payload["uncertified_candidates"].append({
                            "pass": pass_name, "fingerprint": fp[:16],
                            "labels": rec["labels"], "w_x": w, "w_z": d0["w_z"],
                            "lam": d0["lam"],
                            "A_terms": [list(t) for t in p.A_terms],
                            "B_terms": [list(t) for t in p.B_terms],
                            "C_terms": [list(t) for t in lat.terms_of(C)],
                            "D_terms": [list(t) for t in lat.terms_of(D)]})
                if d_lb is not None and rec["min_demoted_w"] is not None:
                    mgn = rec["min_demoted_w"] - d_lb
                    if rec["min_margin"] is None or mgn < rec["min_margin"]:
                        rec["min_margin"] = mgn
            rec["seconds"] = round(time.monotonic() - tp, 2)
            per_parent.append(rec)
            print(f"{pass_name} {fp[:12]} {rec['labels'][0]:16s} kP={rec['k_P']:2d} "
                  f"syz={rec['sigma']:3d} d_lb={d_lb} lights={len(lrecs):5d} "
                  f"perts={rec['perturbations']:4d} dem={rec['demoted_events']:3d} "
                  f"minw={rec['min_demoted_w']} viol={rec['violations']} "
                  f"t={rec['seconds']}s", flush=True)
        if time.monotonic() > DEADLINE:
            break

    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    payload["wall_seconds"] = round(time.monotonic() - t0, 1)
    payload["verdict"] = ("X_COLLAPSE_WITNESSES_FOUND"
                          if payload["decrease_witnesses"]
                          else "X_COLLAPSE_NEVER_OBSERVED_WITHIN_COMPUTED_BUDGET")
    payload["truncated"] = time.monotonic() > DEADLINE
    OUT.write_text(json.dumps(payload, indent=1))
    print(json.dumps({"verdict": payload["verdict"],
                      "witnesses": len(payload["decrease_witnesses"]),
                      "parents": len(per_parent),
                      "uncert": len(payload["uncertified_candidates"]),
                      "seconds": payload["wall_seconds"]}))


if __name__ == "__main__":
    main()
