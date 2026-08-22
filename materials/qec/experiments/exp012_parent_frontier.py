"""EXP-012: universal parent-domination test (works for every delta).

For a PBB code Q = PBB(A,B,C,D) with reported distance D (an upper bound on
d(Q) at every trust level), its parent P = BB(A,B) satisfies

    n(P) = n(Q)            (same qubits)
    k(P) = k(Q) + delta    >= k(Q)          (Proposition 3)
    d(P) = d_Z(P)                            (Proposition 4)

so P dominates Q in [[n,k,d]] as soon as d_Z(P) >= D.  That is a single
*refutation* query per logical sector:

    is there a nontrivial Z-logical of P of weight <= D-1 ?

INFEASIBLE on every sector  ->  d(P) >= D  ->  P dominates Q.
A witness instead gives an explicit low-weight parent logical, which is
recorded (it does NOT refute the PBB distance unless delta = 0, in which case
Theorem 1 makes the two logical sets identical).

Note this subsumes exp006: for delta = 0 the answer is guaranteed to be
"dominates" by Corollary 1, so those cases act as a built-in control.
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import BBSpec, PBBSpec, bb_stabilizer, build_bb  # noqa: E402
from qec_research.codes.pbb_theory import analyse_pbb  # noqa: E402
from qec_research.distance.exact import min_weight_with_parity  # noqa: E402
from qec_research.gf2.linalg import nullspace_np, rank_np, rref_np  # noqa: E402

CATALOG = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def quotient_basis(Hc: np.ndarray, Hs: np.ndarray) -> list[np.ndarray]:
    ker = nullspace_np(Hc)
    R, _ = rref_np(Hs)
    cur = [r for r in R]
    r0 = len(cur)
    reps = []
    for row in ker:
        if rank_np(np.array(cur + [row], dtype=np.uint8)) > r0:
            reps.append(row)
            cur.append(row)
            r0 += 1
    return reps


def test_one(r: dict, tl: float, workers: int) -> dict:
    t0 = time.time()
    spec = PBBSpec(r["ell"], r["m"],
                   [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                   [tuple(t) for t in (r["C_terms"] or [])],
                   [tuple(t) for t in (r["D_terms"] or [])])
    st = analyse_pbb(spec)
    parent = BBSpec(spec.ell, spec.m, spec.A, spec.B)
    HX, HZ = build_bb(parent)
    D = int(r["d"])
    n = HX.shape[1]

    # Z-logicals of P live in ker(H_X)/rowspace(H_Z); nontriviality is detected
    # by pairing with an X-type representative (opposite sector).
    detectors = quotient_basis(HZ, HX)
    groups = [[j] for j in range(n)]
    all_inf, witness, wweight, decided = True, None, None, True
    for det in detectors:
        if not det.any():
            continue
        val, lb, status, sol = min_weight_with_parity(
            HX, groups, det, n, time_limit_s=tl, workers=workers, upper_bound=D - 1)
        if status == "INFEASIBLE":
            continue
        if status in ("OPTIMAL", "FEASIBLE"):
            all_inf = False
            if wweight is None or (val is not None and val < wweight):
                wweight, witness = val, sol
        else:
            decided = False
            break
    parent_k = bb_stabilizer(parent).k
    return {
        "code_id": r.get("code_id"), "n": r["n"], "k_pbb": r["k"], "d_reported": D,
        "trust_level": r.get("trust_level"), "delta": st.delta, "parent_k": parent_k,
        "parent_dZ_ge_D": bool(decided and all_inf),
        "parent_low_weight_logical": wweight,
        "parent_dominates": bool(decided and all_inf and parent_k >= r["k"]),
        "decided": decided,
        "wall_time_s": round(time.time() - t0, 1),
    }


def _job(a):
    r, tl, w = a
    try:
        return test_one(r, tl, w)
    except Exception as e:
        return {"code_id": r.get("code_id"), "n": r.get("n"), "error": repr(e)}


def main() -> None:
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 108
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    per = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    delta_pos_only = len(sys.argv) > 5 and sys.argv[5] == "deltapos"

    rows = [json.loads(l) for l in open(CATALOG) if json.loads(l)["n"] <= max_n]
    if delta_pos_only:
        out_rows = []
        for r in rows:
            spec = PBBSpec(r["ell"], r["m"],
                           [tuple(t) for t in r["A_terms"]], [tuple(t) for t in r["B_terms"]],
                           [tuple(t) for t in (r["C_terms"] or [])],
                           [tuple(t) for t in (r["D_terms"] or [])])
            if analyse_pbb(spec).delta > 0:
                out_rows.append(r)
        rows = out_rows
    print(f"testing {len(rows)} PBB codes (n<={max_n}, delta>0 only={delta_pos_only})",
          flush=True)

    out = []
    with ProcessPoolExecutor(max_workers=procs) as ex:
        futs = [ex.submit(_job, (r, tl, per)) for r in rows]
        for i, f in enumerate(as_completed(futs)):
            rec = f.result()
            out.append(rec)
            if "error" in rec:
                print(f"[{i+1}/{len(rows)}] ERROR {rec.get('code_id')} {rec['error'][:80]}", flush=True)
            else:
                tag = ("PARENT DOMINATES" if rec["parent_dominates"]
                       else ("parent weaker" if rec["decided"] else "undecided"))
                print(f"[{i+1}/{len(rows)}] {rec['code_id']} n={rec['n']} "
                      f"k={rec['k_pbb']}->{rec['parent_k']} d={rec['d_reported']} "
                      f"delta={rec['delta']} -> {tag}"
                      + (f" (parent logical wt {rec['parent_low_weight_logical']})"
                         if rec["parent_low_weight_logical"] is not None else ""),
                      flush=True)

    tag = f"n{max_n}" + ("_deltapos" if delta_pos_only else "")
    (OUT / f"exp012_parent_frontier_{tag}.json").write_text(json.dumps(out, indent=2))
    good = [x for x in out if "error" not in x]
    print(f"\nSUMMARY over {len(good)} codes")
    print(f"  parent dominates : {sum(1 for x in good if x['parent_dominates'])}")
    print(f"  parent weaker    : {sum(1 for x in good if x['decided'] and not x['parent_dominates'])}")
    print(f"  undecided        : {sum(1 for x in good if not x['decided'])}")
    print("  by delta:", dict(sorted(Counter(
        (x["delta"], x["parent_dominates"]) for x in good).items())))


if __name__ == "__main__":
    main()
