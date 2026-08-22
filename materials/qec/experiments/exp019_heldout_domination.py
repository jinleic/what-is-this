"""EXP-019: does any held-out PBB code escape CSS domination?  (Track B, decisive part)

EXP-018 established the *structural* predicates on (12,12) and (15,12).  It did
not answer Track B's actual question, which is about distance.  This does.

The test used here is the universal one (valid for every delta), not
Corollary 1 (delta = 0 only).  For Q = PBB(A,B,C,D) with parent P = BB(A,B):

    n(P) = n(Q)                                  same qubits
    k(P) = k(Q) + delta  >=  k(Q)                Prop. 3
    d(P) >= d(Q)  ?                              <-- the thing to decide

If the last line holds, P is a CSS code that weakly dominates Q in [[n,k,d]],
so Q is not on the frontier and cannot exceed the CSS envelope.

Both halves are cheap *in the right direction*:

  upper bound on d(Q)   a single pure-Z logical of Q is a witness; by Props 1-2
                        the pure-Z minimum of Q equals d_Z of its shadow, and
                        d(Q) <= that minimum.  Finding a witness is a
                        feasibility query, not an optimisation.
  lower bound on d(P)   ask CP-SAT for a logical of P of weight < d_up(Q).
                        INFEASIBLE in every sector proves d(P) >= d_up(Q).
                        This is cheap precisely because d_up(Q) is small --
                        which is the regime that matters.

Outcomes per code:
  DOMINATED    d(P) >= d_up(Q) >= d(Q) proved, and k(P) >= k(Q) at equal n, so P
               weakly dominates Q in [[n,k,d]].  This direction is SOUND: both
               inequalities are certified, one by INFEASIBLE-in-every-sector and
               one by an explicit witness.
  ESCALATE     a logical of P lighter than d_up(Q) was exhibited.  This does
               NOT show d(P) < d(Q) and must never be read that way: d_up(Q) is
               only an UPPER bound on d(Q), so Q may have an equally light or
               lighter logical -- including a mixed one that the pure-Z query
               cannot see.  Deciding it needs a certified lower bound on d(Q),
               which is the expensive direction; these are queued for it.
  UNDECIDED    solver did not finish; recorded, never counted as success
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "e18", str(ROOT / "experiments" / "exp018_heldout_lattices.py"))
e18 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(e18)

from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec, PBBSpec, build_bb, build_pbb, poly_matrix)
from qec_research.codes.pbb_theory import analyse_pbb  # noqa: E402
from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.distance.exact import min_weight_with_parity  # noqa: E402
from qec_research.symplectic.core import lambda_swap, symplectic_weight  # noqa: E402

OUT = ROOT / "results" / "processed"


def exists_logical_within(code, W: int, tl: float, workers: int) -> tuple:
    """Decide  d(code) <= W  by feasibility, with early exit.

    Returns (found, weight, all_decided).  `found=True` exhibits a nontrivial
    logical of weight <= W, so d(code) <= W -- certified by a witness.
    `found=False` with all_decided=True proves d(code) > W.
    """
    n = code.n
    L = code.logical_basis()
    # Membership in the symplectic centralizer S^perp is (Lambda H) v = 0, NOT
    # H v = 0.  Using the ordinary dot product here silently admits vectors that
    # are not logicals at all on mixed (non-CSS) checks.
    constraints = lambda_swap(code.H)
    groups = [[j, j + n] for j in range(n)]      # symplectic weight groups
    all_decided = True
    for j in range(L.shape[0]):
        prow = lambda_swap(L[j][None, :])[0]
        val, lb, status, sol = min_weight_with_parity(
            constraints, groups, prow, 2 * n,
            time_limit_s=tl, workers=workers, upper_bound=W)
        if val is not None:
            # never trust the solver: verify the witness independently
            v = np.zeros(2 * n, dtype=np.uint8)
            v[sol] = 1
            if not code.is_logical(v):
                raise AssertionError("solver returned a non-logical witness")
            if symplectic_weight(v) != int(val):
                raise AssertionError("witness weight disagrees with objective")
            return True, int(val), True
        if status != "INFEASIBLE":
            all_decided = False
    return False, None, all_decided


def parent_distance(ell: int, m: int, A: list, B: list, tl: float,
                    workers: int) -> dict:
    """Exact distance of the CSS parent, computed ONCE and reused by every
    perturbation of it.  This is where the expensive work now lives."""
    HX, HZ = build_bb(BBSpec(ell=ell, m=m, A=A, B=B))
    t0 = time.time()
    r = exact_distance_css(HX, HZ, time_limit_s=tl, workers=workers)
    return {"d": r["d"], "d_X": r["d_X"], "d_Z": r["d_Z"],
            "exact": bool(r["d_exact"]), "wall_s": round(time.time() - t0, 1)}


def domination_test(ell: int, m: int, A: list, B: list, C: list, D: list,
                    dP: dict, tl: float, workers: int) -> dict:
    """Is the PBB child dominated by its already-measured CSS parent?

    Parent P has n(P)=n(Q) and k(P)=k(Q)+delta >= k(Q) (Prop. 3).  So P weakly
    dominates Q in [[n,k,d]] as soon as d(P) >= d(Q).  We test that by asking
    whether Q has ANY nontrivial logical of weight <= d(P):

      feasible   -> d(Q) <= d(P), witness in hand      -> DOMINATED   (sound)
      infeasible -> d(Q) >  d(P), proved in every sector -> BEATS_PARENT
                    (a real finding: the perturbation raised the distance,
                     though k fell by delta -- still has to beat the whole CSS
                     envelope at its own k to matter for Track B)
    """
    t0 = time.time()
    pspec = PBBSpec(ell=ell, m=m, A=A, B=B, C=C, D=D)
    st = analyse_pbb(pspec)
    Q = build_pbb(pspec)
    rec = {"A": A, "B": B, "C": C, "D": D, "n": st.n, "k_pbb": st.k_pbb,
           "k_parent": st.k_bb, "delta": st.delta,
           "parent_d": dP["d"], "parent_d_exact": dP["exact"]}
    if dP["d"] is None or not dP["exact"]:
        rec.update(verdict="UNDECIDED", reason="parent distance not certified",
                   wall_s=round(time.time() - t0, 1))
        return rec
    found, w, decided = exists_logical_within(Q, dP["d"], tl, workers)
    if found:
        rec.update(verdict="DOMINATED", pbb_logical_weight_found=w)
    elif decided:
        rec.update(verdict="BEATS_PARENT", pbb_d_lower_bound=dP["d"] + 1)
    else:
        rec.update(verdict="UNDECIDED", reason="feasibility query hit the time limit")
    rec["wall_s"] = round(time.time() - t0, 1)
    return rec


if __name__ == "__main__":
    lattice = sys.argv[1]
    kmin = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    n_parents = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    max_w = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    tl = float(sys.argv[6]) if len(sys.argv) > 6 else 120.0
    workers = int(sys.argv[7]) if len(sys.argv) > 7 else 12
    max_codes = int(sys.argv[8]) if len(sys.argv) > 8 else 10**9
    n_sectors = int(sys.argv[9]) if len(sys.argv) > 9 else 4
    ell, m = (int(x) for x in lattice.split("x"))

    src = json.loads((OUT / f"exp018_heldout_{lattice}.json").read_text())
    reps = src["stage1_dimension_envelope"]["representatives_by_k"]
    As, Bs = e18.canonical_A(ell, m), e18.canonical_B(ell, m)

    # QEC-relevant parents only.  Selecting by max k (as EXP-018 stage 2 did)
    # picks degenerate codes like k=128 at n=288, which have tiny distance and
    # are irrelevant to the envelope question.
    parents = []
    for k in sorted((int(x) for x in reps), reverse=True):
        if not (kmin <= k <= kmax):
            continue
        for rep in reps[str(k)][:n_parents]:
            # resolved from TERMS, never from an index into a regenerated family
            parents.append((rep["A"], rep["B"], k))
    print(f"=== EXP-019 ({ell},{m}) n={2*ell*m}: {len(parents)} QEC-relevant parents "
          f"with {kmin} <= k <= {kmax} ===", flush=True)

    rows = []
    for A, B, k in parents:
        if len(rows) >= max_codes:
            print(f"  budget cap {max_codes} reached; remaining parents not probed",
                  flush=True)
            break
        dP = parent_distance(ell, m, A, B, tl, workers)
        basis = e18.commutation_nullspace(ell, m, poly_matrix(ell, m, A),
                                          poly_matrix(ell, m, B))
        cands, complete = e18.all_members_up_to_weight(basis, ell, m, max_w, cap=10000)
        print(f"\n  parent k={k} d={dP['d']} (exact={dP['exact']}, {dP['wall_s']}s) "
              f"A={A} B={B}: {len(cands)} perturbations of weight<={max_w} "
              f"(enumeration complete={complete})", flush=True)
        for C, D in cands:
            if len(rows) >= max_codes:
                break
            rec = domination_test(ell, m, A, B, C, D, dP, tl, workers)
            rec["enumeration_complete"] = complete
            rows.append(rec)
            if rec["verdict"] != "DOMINATED":
                print(f"    {rec['verdict']}: k={rec.get('k_pbb')} delta={rec.get('delta')} "
                      f"d_up={rec.get('d_pbb_upper_bound')} "
                      f"parent_lighter={rec.get('parent_lighter_logical_weight')}",
                      flush=True)
        done = sum(1 for r in rows if r["verdict"] == "DOMINATED")
        print(f"    running total: {done}/{len(rows)} dominated", flush=True)

    from collections import Counter
    tally = Counter(r["verdict"] for r in rows)
    print(f"\n=== TRACK B VERDICT on ({ell},{m}) ===")
    print(f"  PBB codes tested : {len(rows)}")
    for v, c in tally.most_common():
        print(f"  {v:16s} : {c}")
    esc = [r for r in rows if r["verdict"] in ("BEATS_PARENT", "UNDECIDED")]
    dom = tally.get("DOMINATED", 0)
    print(f"\n  SCOPE: {len(parents)} parents (the first {n_parents} stored "
          f"representatives per k in [{kmin},{kmax}]) and perturbations of weight "
          f"<= {max_w}. This is NOT lattice-wide coverage and is never reported "
          f"as such.")
    print(f"  Within that scope, {dom}/{len(rows)} PBB codes are provably weakly "
          f"dominated by their parent CSS BB code.")
    if esc:
        print(f"  {len(esc)} codes are UNDECIDED in the domination direction: the "
              f"parent has a lighter logical than the PBB upper bound, which does "
              f"NOT make the parent weaker. Each needs a certified lower bound on "
              f"d(Q) before any frontier claim.")
        for r in esc[:8]:
            print(f"    A={r['A']} B={r['B']} C={r['C']} D={r['D']} "
                  f"k={r['k_pbb']} delta={r['delta']} d_up(Q)={r['d_pbb_upper_bound']} "
                  f"parent logical weight {r['parent_lighter_logical_weight']}")

    (OUT / f"exp019_domination_{lattice}.json").write_text(json.dumps(
        {"lattice": lattice, "n": 2 * ell * m, "k_range": [kmin, kmax],
         "max_perturbation_weight": max_w, "time_limit_s": tl,
         "tally": dict(tally), "escalated_undecided": esc, "rows": rows,
         "scope": (f"{len(parents)} parents = first {n_parents} stored representatives "
                   f"per k in [{kmin},{kmax}]; perturbations of weight <= {max_w}; "
                   f"budget cap {max_codes} codes. Not lattice-wide coverage."),
         "soundness_note": ("DOMINATED is certified in both directions. ESCALATE means "
                            "undecided, NOT that the parent is weaker: d_pbb_upper_bound "
                            "is an upper bound on d(Q), so a lighter parent logical is "
                            "consistent with the parent still dominating.")},
        indent=2, default=str))
    print(f"\nwrote results/processed/exp019_domination_{lattice}.json")
