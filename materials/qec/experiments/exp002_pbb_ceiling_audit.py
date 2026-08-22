"""EXP-002: audit catalogue PBB distances against Theorem 1.

For every PBB entry with delta = 0 the nontrivial pure-Z logical operators are
identical to those of the parent CSS BB code.  Therefore the true PBB distance
satisfies  d <= d_Z(BB(A,B)).

Test: ask CP-SAT whether the parent BB code has a nontrivial Z-logical of
weight <= d_reported - 1.  A SAT answer is a *certificate* that the catalogue
distance is an over-estimate (the witness is an explicit logical operator of
the PBB code itself, verified independently).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import PBBSpec, build_pbb  # noqa: E402
from qec_research.codes.pbb_theory import analyse_pbb, parent_bb_matrices  # noqa: E402
from qec_research.distance.exact import min_weight_with_parity  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402
from qec_research.symplectic.core import StabilizerCode, symplectic_weight  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def z_logical_reps(HX: np.ndarray, SZ: np.ndarray) -> list[np.ndarray]:
    """Basis of ker(HX)/rowspace(SZ) -- the pure-Z logical classes."""
    ker = nullspace_np(HX)
    R, _ = rref_np(SZ)
    cur = [r for r in R]
    r0 = len(cur)
    reps = []
    for row in ker:
        if rank_np(np.array(cur + [row], dtype=np.uint8)) > r0:
            reps.append(row)
            cur.append(row)
            r0 += 1
    return reps


def audit_row(r: dict, time_limit: float, workers: int) -> dict:
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    st = analyse_pbb(spec)
    HX, HZ, P = parent_bb_matrices(spec)
    d_rep = r["d"]

    # pure-Z stabilizer subgroup of the PBB code
    from qec_research.gf2.linalg import matmul
    U = nullspace_np(HX.T)
    SZ = np.vstack([HZ, matmul(U, P)]) if U.shape[0] else HZ

    reps = z_logical_reps(HX, SZ)
    n = HX.shape[1]
    groups = [[j] for j in range(n)]
    best = None
    best_vec = None
    proved_all = True
    t0 = time.time()
    for rep in reps:
        cap = (best - 1) if best is not None else (d_rep - 1)
        if cap < 1:
            continue
        val, lb, status, sol = min_weight_with_parity(
            HX, groups, rep, n, time_limit_s=time_limit, workers=workers,
            upper_bound=cap)
        if status not in ("OPTIMAL", "INFEASIBLE"):
            proved_all = False
        if val is not None and (best is None or val < best):
            best, best_vec = val, sol
    return {
        "code_id": r["code_id"], "n": r["n"], "k": r["k"], "d_reported": d_rep,
        "trust_level": r.get("trust_level"), "d_is_exact": r.get("d_is_exact"),
        "delta": st.delta, "theorem1_applies": st.theorem1_applies,
        "z_logical_upper_bound_found": best,
        "violation": (best is not None and best < d_rep),
        "witness_support": best_vec,
        "all_subproblems_decided": proved_all,
        "wall_time_s": round(time.time() - t0, 2),
    }


def verify_witness(r: dict, support: list[int]) -> dict:
    """Independently confirm the witness is a genuine nontrivial PBB logical."""
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    code = build_pbb(spec)
    n = code.n
    v = np.zeros(2 * n, dtype=np.uint8)
    for j in support:
        v[n + j] = 1          # pure-Z operator
    return {
        "is_nontrivial_logical": bool(code.is_logical(v)),
        "symplectic_weight": symplectic_weight(v),
    }


if __name__ == "__main__":
    only_n = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
    wk = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    rows = [json.loads(l) for l in open(CATALOG)]
    if only_n:
        rows = [r for r in rows if r["n"] == only_n]
    out = []
    for i, r in enumerate(rows):
        rec = audit_row(r, tl, wk)
        if rec["violation"]:
            rec["witness_check"] = verify_witness(r, rec["witness_support"])
        out.append(rec)
        flag = "  *** VIOLATION" if rec["violation"] else ""
        print(f"[{i+1}/{len(rows)}] {rec['code_id']} n={rec['n']} k={rec['k']} "
              f"d_rep={rec['d_reported']} delta={rec['delta']} "
              f"dZ_ub={rec['z_logical_upper_bound_found']} "
              f"trust={rec['trust_level']}{flag}", flush=True)
    tag = f"n{only_n}" if only_n else "all"
    (OUT / f"exp002_pbb_ceiling_audit_{tag}.json").write_text(json.dumps(out, indent=2))
    nv = sum(1 for x in out if x["violation"])
    print(f"\nSUMMARY: {nv}/{len(out)} catalogue distances refuted by an explicit witness")
