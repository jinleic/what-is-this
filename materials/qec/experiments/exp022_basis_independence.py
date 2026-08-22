"""EXP-022: is Theorem C4 basis-independent, or only true for the published generators?

Proposition C2 bounds the one-ancilla depth by the weight of the generators the
circuit ACTUALLY MEASURES:

    T  >=  max_a |supp(g_a)| .

Theorem C4 applied this with the catalogue's canonical generators, which have
max weight >= 8, and concluded T >= 8.  But a circuit may measure ANY generating
set of the same stabiliser group, and stabiliser row operations can change
generator weights.  So as stated the argument proves

    T >= 8  for one-ancilla circuits measuring THE PUBLISHED GENERATORS,

not for every one-ancilla circuit for the code.  This script decides which.

THE DECIDING QUANTITY.  Let

    w_mix  =  min { |supp(v)| : v in S,  x-part(v) != 0 } .

Any generating set must contain at least r_X = rank(H_x) elements with nonzero
X-part (they have to span the X-part of the group).  Every such element has
weight >= w_mix.  Hence for EVERY generating set,

    max_a |supp(g_a)|  >=  w_mix        =>        T >= w_mix .

So:
  * w_mix >= 8  ->  Theorem C4 is basis-independent; the claim "any one-ancilla
                    schedule" is justified.
  * w_mix <= 7  ->  a lower-weight generating set may exist and C4 must stay
                    scoped to the published generators.  (w_mix <= 7 does not by
                    itself construct such a basis -- it only removes the proof.)

The pure-Z part needs no separate treatment: the parent's Z-checks already have
weight 6, so the minimum over all of S is at most 6 and says nothing.  The
X-supported elements are where the separation must live.

Encoding: v = c^T H with c in GF(2)^{rows}.  "x-part != 0" is the linear
constraint sum(x-part) >= 1, so the whole thing is one bounded integer program
solved to proven optimality.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import (  # noqa: E402
    BRAVYI_BB, PBBSpec, bb_stabilizer, build_pbb)
from qec_research.symplectic.core import symplectic_weight  # noqa: E402

OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def min_weight_with_x_support(H: np.ndarray, n: int, *, time_limit_s: float,
                              workers: int, upper_bound: int | None = None) -> dict:
    """min |supp(v)| over v in rowspace(H) with nonzero X-part, proven optimal."""
    rows = H.shape[0]
    mdl = cp_model.CpModel()
    c = [mdl.new_bool_var(f"c{i}") for i in range(rows)]
    x = [mdl.new_bool_var(f"x{j}") for j in range(n)]
    z = [mdl.new_bool_var(f"z{j}") for j in range(n)]
    q = [mdl.new_bool_var(f"q{j}") for j in range(n)]

    def parity(col: np.ndarray, target, tag: str) -> None:
        idx = [i for i in range(rows) if col[i]]
        if not idx:
            mdl.add(target == 0)
            return
        t = mdl.new_int_var(0, len(idx) // 2, f"t_{tag}")
        mdl.add(sum(c[i] for i in idx) - 2 * t == target)

    for j in range(n):
        parity(H[:, j], x[j], f"x{j}")
        parity(H[:, n + j], z[j], f"z{j}")
        mdl.add(q[j] >= x[j])
        mdl.add(q[j] >= z[j])
        mdl.add(q[j] <= x[j] + z[j])

    mdl.add(sum(x) >= 1)                      # nonzero X-part
    if upper_bound is not None:
        mdl.add(sum(q) <= upper_bound)
    mdl.minimize(sum(q))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = workers
    st = solver.solve(mdl)
    name = solver.status_name(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        vec = np.zeros(2 * n, dtype=np.uint8)
        for j in range(n):
            vec[j] = solver.value(x[j])
            vec[n + j] = solver.value(z[j])
        return {"status": name, "value": int(solver.objective_value),
                "bound": int(solver.best_objective_bound), "witness": vec}
    return {"status": name, "value": None,
            "bound": int(solver.best_objective_bound) if st != cp_model.UNKNOWN else 0,
            "witness": None}



def xrank_of_light_elements(H: np.ndarray, n: int, wcap: int, *,
                            time_limit_s: float, workers: int,
                            rX: int) -> dict:
    """DECISIVE test for basis-independence of the depth bound.

    Enumerate stabiliser elements v with |supp(v)| <= wcap and x-part != 0, and
    accumulate the GF(2) rank of their x-parts.

      rank_X(V_wcap) <  r_X   =>  the weight-<=wcap elements cannot span the
                                  X-part, so EVERY generating set must contain
                                  an element of weight > wcap.  PROVEN.
      rank_X(V_wcap) == r_X   =>  NECESSARY but NOT SUFFICIENT for refutation.
                                  A generating set must also span the pure-Z
                                  subgroup S_Z, and for a PBB code S_Z is
                                  ENLARGED beyond the parent's H_Z rows by the
                                  {uP} terms (Prop. 2), so the parent's weight-6
                                  Z-checks need not suffice.  Refutation requires
                                  the light vectors -- X-supported AND pure-Z --
                                  to span all of rowspace(H), i.e. rank n-k.

    This is exactly the gap that a single-element w_mix computation cannot close.
    """
    from qec_research.gf2.linalg import rank_np
    rows = H.shape[0]
    mdl = cp_model.CpModel()
    c = [mdl.new_bool_var(f"c{i}") for i in range(rows)]
    x = [mdl.new_bool_var(f"x{j}") for j in range(n)]
    z = [mdl.new_bool_var(f"z{j}") for j in range(n)]
    q = [mdl.new_bool_var(f"q{j}") for j in range(n)]

    def parity(col, target, tag):
        idx = [i for i in range(rows) if col[i]]
        if not idx:
            mdl.add(target == 0); return
        t = mdl.new_int_var(0, len(idx) // 2, f"t_{tag}")
        mdl.add(sum(c[i] for i in idx) - 2 * t == target)

    for j in range(n):
        parity(H[:, j], x[j], f"x{j}")
        parity(H[:, n + j], z[j], f"z{j}")
        mdl.add(q[j] >= x[j]); mdl.add(q[j] >= z[j]); mdl.add(q[j] <= x[j] + z[j])
    mdl.add(sum(x) >= 1)
    mdl.add(sum(q) <= wcap)

    found_x: list[np.ndarray] = []       # x-parts, for the X-rank test
    found_full: list[np.ndarray] = []    # FULL symplectic vectors, for the span test
    state = {"rank": 0, "hit": False, "count": 0}

    class Collect(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            xs = np.array([self.value(v) for v in x], dtype=np.uint8)
            zs = np.array([self.value(v) for v in z], dtype=np.uint8)
            state["count"] += 1
            cand = np.vstack(found_x + [xs]) if found_x else xs[None, :]
            r = rank_np(cand)
            if r > state["rank"]:
                found_x.append(xs)
                found_full.append(np.concatenate([xs, zs]).astype(np.uint8))
                state["rank"] = r
            if state["rank"] >= rX:
                state["hit"] = True
                self.stop_search()

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = 1
    st = solver.solve(mdl, Collect())
    name = solver.status_name(st)
    complete = (name == "OPTIMAL" and not state["hit"])
    return {"weight_cap": wcap, "solutions_seen": state["count"],
            "x_rank_of_light_span": state["rank"], "rank_X_needed": rX,
            "spans_X": state["hit"], "enumeration_complete": complete,
            "status": name,
            "light_vectors": [v.tolist() for v in found_full],
            "NOTE": ("spans_X is NECESSARY but NOT SUFFICIENT for refutation; the "
                     "pure-Z subgroup must also be generated at weight <= cap.")}

def analyse(label: str, code, tl: float, workers: int) -> dict:
    n = code.n
    H = code.H
    supp = ((H[:, :n] | H[:, n:]) != 0)
    published_max_w = int(supp.sum(1).max())
    rX = int(np.linalg.matrix_rank(H[:, :n].astype(float)))  # coarse; refined below
    from qec_research.gf2.linalg import rank_np
    rX = rank_np(H[:, :n])

    t0 = time.time()
    r = min_weight_with_x_support(H, n, time_limit_s=tl, workers=workers)
    w = r["value"]
    ok = r["status"] == "OPTIMAL"

    # verify the witness really is in the group with nonzero X-part
    verified = None
    if r["witness"] is not None:
        v = r["witness"]
        from qec_research.gf2.linalg import row_space_contains, rows_to_bitsets
        in_group = row_space_contains(rows_to_bitsets(H), 2 * n,
                                      rows_to_bitsets(v[None, :])[0])
        verified = bool(in_group and v[:n].any() and symplectic_weight(v) == w)

    # DECISIVE test: can weight-<=7 elements span the X-part at all?
    dec = xrank_of_light_elements(H, n, 7, time_limit_s=tl, workers=workers, rX=rX)

    lb = r["bound"]
    # SOUND criterion: a proven LOWER bound >= 8 makes C4 basis-independent,
    # whether or not the search reached OPTIMAL.  A FEASIBLE value is only an
    # UPPER bound on w_mix and can never establish it.
    # Sufficient certificate A: every X-supported element already weighs >= 8.
    cert_A = bool(lb is not None and lb >= 8)
    # Decisive test B: weight-<=7 elements cannot span the X-part.
    cert_B = bool(dec["enumeration_complete"] and not dec["spans_X"])
    # Refutation requires MORE than full X-span: the light vectors, together with
    # light PURE-Z vectors, must span all of rowspace(H) (rank n-k).  For a PBB
    # code S_Z is enlarged by {uP}, so the parent's weight-6 Z-checks may not do it.
    refuted = False
    span_check = {"attempted": False}
    if dec["spans_X"] and dec.get("light_vectors"):
        from qec_research.gf2.linalg import rank_np as _rk
        LV = np.array(dec["light_vectors"], dtype=np.uint8)
        target_rank = 2 * n - code.k - n      # = n - k  (rank of the stabiliser group)
        target_rank = H.shape[1] // 2 * 0 + (n - code.k)
        # add every published pure-Z row of weight <= cap as a light generator
        pureZ = H[~H[:, :n].any(axis=1)]
        lightZ = pureZ[((pureZ[:, :n] | pureZ[:, n:]) != 0).sum(1) <= dec["weight_cap"]]
        stack = np.vstack([LV] + ([lightZ] if lightZ.size else []))
        got = _rk(stack)
        refuted = bool(got >= target_rank)
        span_check = {"attempted": True, "rank_from_light_vectors": int(got),
                      "rank_needed_n_minus_k": int(target_rank),
                      "light_pure_z_rows_used": int(lightZ.shape[0]),
                      "spans_whole_group": refuted}
    basis_indep = cert_A or cert_B
    if cert_A:
        verdict = "BASIS-INDEPENDENT (all X-supported elements weigh >= 8)"
    elif cert_B:
        verdict = ("BASIS-INDEPENDENT (weight-<=7 elements span only rank "
                   f"{dec['x_rank_of_light_span']} < {rX} of the X-part)")
    elif refuted:
        verdict = ("REFUTED — light vectors span the whole stabiliser group "
                   f"(rank {span_check.get('rank_from_light_vectors')} = n-k), so a "
                   "generating set of max weight <= 7 exists")
    elif dec["spans_X"]:
        verdict = ("UNDECIDED — light elements span the X-part but were NOT shown to "
                   "generate the whole group; the pure-Z subgroup S_Z is enlarged by "
                   "{uP} and may need a heavier generator")
    else:
        verdict = ("UNDECIDED — the sufficient certificate fails and the light-span "
                   "enumeration did not complete; this neither proves nor refutes")
    return {
        "code": label, "n": n, "k": code.k,
        "published_max_check_weight": published_max_w,
        "rank_X_part": rX,
        "w_mix_upper_bound": w, "w_mix_proven_lower_bound": lb,
        "w_mix_exact": w if ok else None,
        "solver_status": r["status"], "proven_optimal": ok,
        "witness_verified": verified,
        "C4_basis_independent": basis_indep,
        "certificate_A_all_X_elements_heavy": cert_A,
        "certificate_B_light_cannot_span": cert_B,
        "refuted_light_basis_exists": refuted,
        "light_span_test": {k2: v2 for k2, v2 in dec.items() if k2 != "light_vectors"},
        "whole_group_span_check": span_check,
        "verdict": verdict,
        "wall_s": round(time.time() - t0, 1),
    }


if __name__ == "__main__":
    tl = float(sys.argv[1]) if len(sys.argv) > 1 else 900.0
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 3

    cat = ROOT / "third_party" / "qcode-discovery" / "results" / "campaign7_publication_merged.jsonl"
    rows = []
    for line in cat.open():
        r = json.loads(line)
        if r.get("n") == 144 and r.get("k") == 12 and r.get("d") == 12 \
                and all(r.get(t) is not None for t in ("A_terms", "B_terms", "C_terms", "D_terms")):
            rows.append(r)

    out = []
    gross = bb_stabilizer(BRAVYI_BB["[[144,12,12]]"])
    print("=== EXP-022: is the depth-8 bound basis-independent? ===\n")
    print("Reference: CSS Gross code")
    g = analyse("CSS Gross [[144,12,12]]", gross, tl, workers)
    out.append(g)
    print(f"  published max check weight = {g['published_max_check_weight']}, "
          f"w_mix <= {g['w_mix_upper_bound']}, proven >= {g['w_mix_proven_lower_bound']} "
          f"({g['solver_status']}) -> {g['verdict']}\n")

    print(f"PBB [[144,12,12]] members (first {limit} of {len(rows)}):")
    for r in rows[:limit]:
        Q = build_pbb(PBBSpec(ell=r["ell"], m=r["m"], A=r["A_terms"], B=r["B_terms"],
                              C=r["C_terms"], D=r["D_terms"]))
        rec = analyse(r.get("code_id", "PBB"), Q, tl, workers)
        out.append(rec)
        print(f"  {rec['code']:12s} published max w = {rec['published_max_check_weight']:2d}, "
              f"w_mix <= {rec['w_mix_upper_bound']}, proven >= {rec['w_mix_proven_lower_bound']} "
              f"[{rec['solver_status']}, witness_ok={rec['witness_verified']}]"
              f"  -> {rec['verdict']}")

    pbb = [o for o in out if o["code"] != "CSS Gross [[144,12,12]]"]
    all_indep = bool(pbb) and all(o["C4_basis_independent"] for o in pbb)
    print(f"\n=== VERDICT ===")
    if all_indep:
        print("  Every probed PBB has a PROVEN lower bound w_mix >= 8, so ANY generating set")
        print("  contains an X-supported generator of weight >= 8 and Theorem C4 holds for")
        print("  any one-ancilla circuit for the code, not just the published generators.")
    else:
        und = [o for o in pbb if o["verdict"].startswith("UNDECIDED")]
        ref = [o for o in pbb if o["verdict"].startswith("REFUTED")]
        print(f"  NOT established: {len(und)} undecided, {len(ref)} refuted.")
        print("  Theorem C4 MUST stay scoped to circuits measuring the PUBLISHED generators.")
        print("  An undecided result means the solver has only an upper bound on w_mix; it")
        print("  neither proves nor refutes basis-independence.")

    (OUT / "exp022_basis_independence.json").write_text(json.dumps(
        {"results": out, "all_pbb_basis_independent": all_indep,
         "interpretation": (
             "w_mix = min support weight of a stabiliser element with nonzero X-part. "
             "Any generating set needs rank(H_x) such elements, so max generator weight "
             ">= w_mix for EVERY basis. w_mix >= 8 makes Theorem C4 basis-independent; "
             "w_mix <= 7 means C4 is proved only for the published generators.")},
        indent=2, default=str))
    print(f"\nwrote results/processed/exp022_basis_independence.json")
