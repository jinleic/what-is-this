"""EXP-003: does the CSS shadow dominate every published PBB code?

For each catalogue entry Q = PBB(A,B,C,D) with reported distance D
(an upper bound on d(Q) in every trust level), build the CSS shadow

    Q' :  X-checks HX = [A B],   Z-checks S_Z = pure-Z stabilisers of Q.

Theorem 2 gives  n(Q')=n(Q), k(Q')=k(Q), and d(Q) <= d_Z(Q').
The domination question reduces to two *refutation* queries with a tight cap:

    (i)  is there a Z-logical of Q' of weight <= D-1 ?
    (ii) is there an X-logical of Q' of weight <= D-1 ?

If both are INFEASIBLE then d(Q') >= D >= d(Q):  the CSS shadow matches Q in
(n,k) and is no worse in distance.  A SAT answer to (i) additionally refutes
the catalogue distance, because a Z-logical of Q' *is* a logical of Q.
"""

from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import PBBSpec, build_pbb  # noqa: E402
from qec_research.codes.pbb_theory import analyse_pbb, css_shadow, parent_bb_matrices  # noqa: E402
from qec_research.distance.exact import min_weight_with_parity  # noqa: E402
from qec_research.gf2.linalg import matmul, nullspace_np, rank_np, rref_np  # noqa: E402
from qec_research.symplectic.core import StabilizerCode, symplectic_weight  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def quotient_basis(Hc: np.ndarray, Hs: np.ndarray) -> list[np.ndarray]:
    """Basis of ker(Hc) / rowspace(Hs)."""
    ker = nullspace_np(Hc)
    R, _ = rref_np(Hs)
    cur = [r for r in R]
    r0 = len(cur)
    reps: list[np.ndarray] = []
    for row in ker:
        if rank_np(np.array(cur + [row], dtype=np.uint8)) > r0:
            reps.append(row)
            cur.append(row)
            r0 += 1
    return reps


def refute_below(Hc: np.ndarray, detectors: list[np.ndarray], cap: int,
                 time_limit: float, workers: int) -> dict:
    """Is there a nontrivial logical of weight <= cap?  Returns witness or proof."""
    n = Hc.shape[1]
    groups = [[j] for j in range(n)]
    all_infeasible = True
    best = None
    best_sol = None
    for det in detectors:
        if not det.any():
            continue
        val, lb, status, sol = min_weight_with_parity(
            Hc, groups, det, n, time_limit_s=time_limit, workers=workers,
            upper_bound=cap)
        if status == "INFEASIBLE":
            continue
        all_infeasible = False
        if val is not None and (best is None or val < best):
            best, best_sol = val, sol
        if status not in ("OPTIMAL", "FEASIBLE"):
            return {"decided": False, "status": status, "witness_weight": best,
                    "witness": best_sol}
    return {"decided": True, "no_logical_at_or_below_cap": all_infeasible,
            "witness_weight": best, "witness": best_sol}


def analyse_one(r: dict, time_limit: float, workers: int) -> dict:
    t0 = time.time()
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    st = analyse_pbb(spec)
    HX, SZ = css_shadow(spec)
    D = int(r["d"])
    n = HX.shape[1]

    pbb = build_pbb(spec)
    logicals = pbb.logical_basis()

    # (i) Z-logicals of Q'  == pure-Z logicals of Q.  Detector = x-part of logicals.
    zdet = [row for row in logicals[:, :n]]
    zres = refute_below(HX, zdet, D - 1, time_limit, workers)

    # (ii) X-logicals of Q' : v in ker(S_Z), detector = Z-logical reps of Q'.
    xdet = quotient_basis(HX, SZ)          # Z-type reps pair with X-type by dot product
    xres = refute_below(SZ, xdet, D - 1, time_limit, workers)

    dominated = (zres.get("decided") and zres.get("no_logical_at_or_below_cap")
                 and xres.get("decided") and xres.get("no_logical_at_or_below_cap"))
    refuted = bool(zres.get("witness_weight") is not None
                   and zres["witness_weight"] < D)
    rec = {
        "code_id": r.get("code_id", f"{r['ell']}_{r['m']}_?"),
        "n": r["n"], "k": r["k"], "d_reported": D,
        "trust_level": r.get("trust_level"), "d_is_exact": r.get("d_is_exact"),
        "delta": st.delta, "k_parent_bb": st.k_bb,
        "shadow_dZ_ge_D": bool(zres.get("no_logical_at_or_below_cap")),
        "shadow_dX_ge_D": bool(xres.get("no_logical_at_or_below_cap")),
        "css_shadow_dominates": bool(dominated),
        "catalogue_distance_refuted": refuted,
        "refutation_weight": zres.get("witness_weight"),
        "refutation_witness": zres.get("witness"),
        "z_decided": zres.get("decided"), "x_decided": xres.get("decided"),
        "wall_time_s": round(time.time() - t0, 1),
    }
    if refuted:
        v = np.zeros(2 * n, dtype=np.uint8)
        for j in zres["witness"]:
            v[n + j] = 1
        rec["witness_is_nontrivial_logical"] = bool(pbb.is_logical(v))
        rec["witness_symplectic_weight"] = symplectic_weight(v)
    return rec


def _job(args):
    r, tl, wk = args
    try:
        return analyse_one(r, tl, wk)
    except Exception as exc:  # keep the campaign alive, record the failure
        return {"code_id": r.get("code_id"), "n": r.get("n"), "error": repr(exc)}


if __name__ == "__main__":
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 216
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    per = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    rows = [json.loads(l) for l in open(CATALOG)]
    rows = [r for r in rows if r["n"] <= max_n]
    print(f"analysing {len(rows)} catalogue entries with n <= {max_n} "
          f"({procs} processes x {per} CP-SAT workers)", flush=True)
    out = []
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = {ex.submit(_job, (r, tl, per)): r for r in rows}
        for i, f in enumerate(as_completed(futs)):
            rec = f.result()
            out.append(rec)
            if "error" in rec:
                print(f"[{i+1}/{len(rows)}] {rec['code_id']} ERROR {rec['error']}", flush=True)
            else:
                tag = "DOMINATED" if rec["css_shadow_dominates"] else "undecided"
                if rec["catalogue_distance_refuted"]:
                    tag = f"*** REFUTED d<={rec['refutation_weight']} (claim {rec['d_reported']})"
                print(f"[{i+1}/{len(rows)}] {rec['code_id']} n={rec['n']} k={rec['k']} "
                      f"d={rec['d_reported']} delta={rec['delta']} -> {tag} "
                      f"[{rec['wall_time_s']}s]", flush=True)
    (OUT / f"exp003_css_shadow_domination_n{max_n}.json").write_text(json.dumps(out, indent=2))
    dom = sum(1 for x in out if x.get("css_shadow_dominates"))
    ref = sum(1 for x in out if x.get("catalogue_distance_refuted"))
    err = sum(1 for x in out if "error" in x)
    print(f"\nSUMMARY  dominated={dom}/{len(out)}  refuted={ref}  errors={err}")
