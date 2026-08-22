"""EXP-018: Track B -- the previously excluded (ell,m) = (12,12) and (15,12) lattices.

WHAT THE THEORY ACTUALLY PREDICTS HERE (quoting proofs/pbb_structure.md, not a
paraphrase of it):

  Prop. 3   (all delta)  k(Q) = k(BB(A,B)) - delta,  delta >= 0.
  Prop. 1-2 (all delta)  the pure-Z elements of rowspace(H) are exactly
                         S_Z = rowspace(H_Z^parent) + <u[C,D]>, and when
                         delta = 0 the pure-Z logicals of Q are literally the
                         parent's.
  Prop. 4   (BB only)    d_X = d_Z for a bivariate-bicycle code.
  Theorem 2 (all delta)  the CSS code Q' = (X-checks H_X, Z-checks S_Z) has
                         n(Q') = n(Q), k(Q') = k(Q), and d(Q) <= d_Z(Q').
                         This is a pure-Z CEILING.  It does NOT bound d_X(Q'),
                         so it does NOT say Q' dominates Q.
  Cor. 1    (delta = 0)  the *parent* BB(A,B) has the same n and k and
                         d(parent) >= d(Q), so it weakly dominates Q in
                         [[n,k,d]].  This needs Prop. 4 and fails for delta>0:
                         `phase2_58`/`phase2_60` are explicit escapers.

So the Track B question splits cleanly, and only half of it needs distances:

  delta = 0 codes : Corollary 1 settles it with NO distance computation --
                    the parent is a CSS code at the same (n,k) that is at
                    least as good.  A held-out failure would have to show up
                    as a violation of Prop. 3 or Prop. 4, both exact linear
                    algebra.
  delta > 0 codes : NOT settled by theory.  Needs real distances: we compute
                    d(Q) (symplectic), d_Z(Q') and d_X(Q') and ask whether any
                    code escapes its own shadow, as two catalogue codes do.

Stages
------
1   exhaustive exact k over the Bravyi canonical family on the lattice.
    This is the DIMENSION envelope.  It is not a rate-distance envelope and is
    never reported as one.
1b  certified distances (CP-SAT, two-sided with witness) for the best-k
    parents, giving genuine CSS rate-distance points on this lattice.
2   the admissible perturbations (C,D) are an exact GF(2) null space; all
    members of weight <= W are enumerated EXHAUSTIVELY by CP-SAT
    all-solution search, not sampled.
3   held-out structural tests of Prop. 3 / Prop. 4 / Prop. 1-2 on every code
    built in stage 2.
4   for delta > 0 codes, certified sector minima to test the Theorem-2
    ceiling d(Q) <= d_Z(Q') and to look for shadow escapers.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_research.codes.bicycle import (  # noqa: E402
    BBSpec, PBBSpec, bb_stabilizer, build_bb, build_pbb, commutation_defect,
    monomial_matrix, poly_matrix)
from qec_research.codes.pbb_theory import analyse_pbb, css_shadow  # noqa: E402
from qec_research.distance.exact import exact_distance_css  # noqa: E402
from qec_research.distance.sectors import min_pure_z_logical  # noqa: E402
from qec_research.gf2.linalg import (  # noqa: E402
    matmul, nullspace_np, rank_bitset, rank_np, rows_to_bitsets, rref_np)

OUT = ROOT / "results" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- stage 1
def canonical_A(ell: int, m: int) -> list[list[tuple[int, int]]]:
    """A = x^a1 + y^a2 + y^a3 ,  a2 < a3   (Bravyi et al. convention).

    Members with a repeated monomial are EXCLUDED: a1 = a2 = 0 makes (0,0)
    appear twice and XOR-cancel, so the "trinomial" is secretly a monomial and
    the resulting code is not in the intended family at all.
    """
    out = [[(a1, 0), (0, a2), (0, a3)]
           for a1 in range(ell) for a2, a3 in itertools.combinations(range(m), 2)]
    return [t for t in out if len(set(t)) == 3]


def canonical_B(ell: int, m: int) -> list[list[tuple[int, int]]]:
    """B = y^b1 + x^b2 + x^b3 ,  b2 < b3.  Repeated monomials excluded."""
    out = [[(0, b1), (b2, 0), (b3, 0)]
           for b1 in range(m) for b2, b3 in itertools.combinations(range(ell), 2)]
    return [t for t in out if len(set(t)) == 3]


def _k_of_pair(args) -> tuple[int, int, int]:
    ell, m, ia, A, ib, B = args
    HX = np.hstack([poly_matrix(ell, m, A), poly_matrix(ell, m, B)]).astype(np.uint8)
    r = rank_bitset(rows_to_bitsets(HX), HX.shape[1])
    return ia, ib, 2 * (ell * m - r)


def stage1(ell: int, m: int, workers: int) -> dict:
    As, Bs = canonical_A(ell, m), canonical_B(ell, m)
    t0 = time.time()
    total = len(As) * len(Bs)
    jobs = ((ell, m, ia, A, ib, B) for ia, A in enumerate(As) for ib, B in enumerate(Bs))
    hist: dict[int, int] = {}
    best: dict[int, list] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (ia, ib, k) in enumerate(ex.map(_k_of_pair, jobs, chunksize=1024)):
            hist[k] = hist.get(k, 0) + 1
            if k > 0:
                best.setdefault(k, [])
                if len(best[k]) < 40:
                    best[k].append({"A": As[ia], "B": Bs[ib]})
            if (i + 1) % 200000 == 0:
                print(f"    {i+1:,}/{total:,}  {time.time()-t0:.0f}s", flush=True)
    return {"ell": ell, "m": m, "n": 2 * ell * m, "pairs_enumerated": total,
            "A_family_size": len(As), "B_family_size": len(Bs),
            "exhaustive_over": ("Bravyi canonical family A=x^a1+y^a2+y^a3, "
                                "B=y^b1+x^b2+x^b3, members with a repeated "
                                "monomial excluded (they reduce to monomials)"),
            "k_histogram": {str(k): v for k, v in sorted(hist.items())},
            "representatives_by_k": {str(k): v for k, v in sorted(best.items())},
            "wall_s": round(time.time() - t0, 1)}


# ---------------------------------------------------------------- stage 2
def commutation_nullspace(ell: int, m: int, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """EXACT solution space of  A C^T + B D^T symmetric.

    The map (C,D) -> M + M^T is GF(2)-linear, so the admissible perturbations
    form a subspace; we return a basis.  Column j is C's monomial x^(j//m)y^(j%m),
    column dim+j is D's.
    """
    dim = ell * m
    iu = np.triu_indices(dim, k=1)
    cols = []
    for which in ("C", "D"):
        for j in range(dim):
            e = monomial_matrix(ell, m, j // m, j % m)
            M = matmul(A, e.T) if which == "C" else matmul(B, e.T)
            cols.append((M ^ M.T).astype(np.uint8)[iu])
    L = np.array(cols, dtype=np.uint8).T
    R, piv = rref_np(L)                      # <= 2*dim independent constraints
    R = R[:len(piv)]
    return nullspace_np(R) if len(piv) else np.eye(2 * dim, dtype=np.uint8)


def all_members_up_to_weight(basis: np.ndarray, ell: int, m: int, max_weight: int,
                             cap: int, time_limit_s: float = 120.0) -> tuple[list, bool]:
    """EXHAUSTIVE enumeration of every (C,D) in the space with 1 <= wt <= W.

    CP-SAT all-solution search over the basis coefficients, so this is a
    complete list, not a sample.  Returns (members, complete) where `complete`
    is False only if the cap or the time limit truncated the search.
    """
    dim = ell * m
    nb, width = basis.shape
    if nb == 0:
        return [], True
    mdl = cp_model.CpModel()
    c = [mdl.new_bool_var(f"c{i}") for i in range(nb)]
    v = [mdl.new_bool_var(f"v{j}") for j in range(width)]
    for j in range(width):
        idx = [i for i in range(nb) if basis[i, j]]
        if not idx:
            mdl.add(v[j] == 0)
            continue
        t = mdl.new_int_var(0, len(idx) // 2, f"t{j}")
        mdl.add(sum(c[i] for i in idx) - 2 * t == v[j])
    mdl.add(sum(v) >= 1)
    mdl.add(sum(v) <= max_weight)

    found: list[tuple[list, list]] = []

    class Collect(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.truncated = False

        def on_solution_callback(self):
            vec = np.array([self.value(x) for x in v], dtype=np.uint8)
            C = [(int(j) // m, int(j) % m) for j in np.flatnonzero(vec[:dim])]
            D = [(int(j) // m, int(j) % m) for j in np.flatnonzero(vec[dim:])]
            found.append((C, D))
            if len(found) >= cap:
                self.truncated = True
                self.stop_search()

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = 1                 # required for enumeration
    cb = Collect()
    status = solver.solve(mdl, cb)
    complete = (status == cp_model.OPTIMAL and not cb.truncated)
    return found, complete


# ---------------------------------------------------------------- stage 3
def structural_tests(ell: int, m: int, A: list, B: list, C: list, D: list) -> dict:
    """Held-out Prop. 3 / Prop. 4 / Prop. 1-2 checks.  Exact linear algebra."""
    pspec = PBBSpec(ell=ell, m=m, A=A, B=B, C=C, D=D)
    st = analyse_pbb(pspec)
    HXp, HZp = build_bb(BBSpec(ell=ell, m=m, A=A, B=B))
    HXs, HZs = css_shadow(pspec)
    defect = commutation_defect(poly_matrix(ell, m, A), poly_matrix(ell, m, B),
                                poly_matrix(ell, m, C), poly_matrix(ell, m, D))
    shadow_k = int(HXs.shape[1] - rank_np(HXs) - rank_np(HZs))
    return {
        "ell": ell, "m": m, "A": A, "B": B, "C": C, "D": D, "n": st.n,
        "k_pbb": st.k_pbb, "k_parent_bb": st.k_bb, "delta": st.delta,
        "commutes_exactly": bool(not defect.any()),
        "P3_k_identity": st.k_pbb == st.k_bb - st.delta,
        "P3_delta_nonneg": st.delta >= 0,
        # Prop. 4 on the PARENT (a BB code), where it is claimed
        "P4_parent_rank_symmetry": bool(
            rank_np(HXp) == rank_np(HZp)),
        # Theorem 2's shadow must have the SAME k as Q -- not k + delta
        "T2_shadow_k_equals_k_pbb": shadow_k == st.k_pbb,
        "shadow_k": shadow_k,
        "P2_applies_delta_zero": st.delta == 0,
        "P2_z_logicals_identical": bool(st.z_logicals_identical),
        "max_check_weight_pbb": int(((build_pbb(pspec).H[:, :st.n]
                                      | build_pbb(pspec).H[:, st.n:]) != 0).sum(1).max()),
        "max_check_weight_parent": int(len(A) + len(B)),
    }


# ---------------------------------------------------------------- stage 4
def distance_tests(rec: dict, tl: float, workers: int) -> dict:
    """Certified sector minima: the Theorem-2 ceiling and shadow escape.

    Only meaningful for delta > 0, where Corollary 1 does not apply.
    """
    ell, m = rec["ell"], rec["m"]
    pspec = PBBSpec(ell=ell, m=m, A=rec["A"], B=rec["B"], C=rec["C"], D=rec["D"])
    code = build_pbb(pspec)
    HXs, HZs = css_shadow(pspec)
    t0 = time.time()
    dz_q = min_pure_z_logical(code, time_limit_s=tl, workers=workers)
    sh = exact_distance_css(HXs, HZs, time_limit_s=tl, workers=workers)
    out = {
        "pure_z_min_weight_Q": dz_q.weight, "pure_z_exact": dz_q.exact,
        "shadow_d_X": sh["d_X"], "shadow_d_X_exact": sh["d_X_exact"],
        "shadow_d_Z": sh["d_Z"], "shadow_d_Z_exact": sh["d_Z_exact"],
        "wall_s": round(time.time() - t0, 1),
    }
    # Theorem 2 says d(Q) <= d_Z(Q').  d(Q) <= pure-Z minimum of Q, and by
    # Props 1-2 that minimum IS d_Z(Q'); check the identity the theorem rests on.
    if dz_q.exact and sh["d_Z_exact"]:
        out["T2_identity_holds"] = (dz_q.weight == sh["d_Z"])
    if sh["d_X_exact"] and sh["d_Z_exact"]:
        out["shadow_d"] = min(sh["d_X"], sh["d_Z"])
        out["shadow_escapes"] = sh["d_X"] < sh["d_Z"]
    return out


if __name__ == "__main__":
    lattice = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    max_w = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    n_parents = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    ell, m = (int(x) for x in lattice.split("x"))

    print(f"=== Stage 1: EXHAUSTIVE exact k, ({ell},{m}), n={2*ell*m} ===", flush=True)
    s1 = stage1(ell, m, workers)
    print(f"  k histogram (dimension envelope only, NO distances): {s1['k_histogram']}")
    print(f"  [{s1['wall_s']}s]", flush=True)

    As, Bs = canonical_A(ell, m), canonical_B(ell, m)
    ks = sorted((int(k) for k in s1["representatives_by_k"]), reverse=True)
    parents = []
    for k in ks:
        for rep in s1["representatives_by_k"][str(k)][:3]:
            parents.append((rep["A"], rep["B"], k))
        if len(parents) >= n_parents:
            break
    parents = parents[:n_parents]

    print(f"\n=== Stage 2/3: exact perturbation spaces on {len(parents)} parents ===",
          flush=True)
    results, spaces = [], []
    for A, B, k in parents:
        t0 = time.time()
        basis = commutation_nullspace(ell, m, poly_matrix(ell, m, A),
                                      poly_matrix(ell, m, B))
        cands, complete = all_members_up_to_weight(basis, ell, m, max_w, cap=200)
        spaces.append({"A": A, "B": B, "k_parent": k,
                       "solution_space_dim": int(basis.shape[0]),
                       "members_up_to_weight": max_w,
                       "count": len(cands), "enumeration_complete": complete})
        print(f"  k={k} A={A} B={B}: dim={basis.shape[0]}, "
              f"{len(cands)} perturbations of weight<={max_w} "
              f"(complete={complete})  [{time.time()-t0:.1f}s]", flush=True)
        for C, D in cands:
            try:
                results.append(structural_tests(ell, m, A, B, C, D))
            except Exception as e:                          # noqa: BLE001
                results.append({"A": A, "B": B, "C": C, "D": D, "error": str(e)})

    ok = [r for r in results if "error" not in r]
    fails = {name: [r for r in ok if not r[name]] for name in
             ("commutes_exactly", "P3_k_identity", "P3_delta_nonneg",
              "P4_parent_rank_symmetry", "T2_shadow_k_equals_k_pbb")}
    p2 = [r for r in ok if r["P2_applies_delta_zero"]]
    p2bad = [r for r in p2 if not r["P2_z_logicals_identical"]]

    print(f"\n=== HELD-OUT STRUCTURAL VERDICT on ({ell},{m}) ===")
    print(f"  non-CSS PBB codes constructed and tested : {len(ok)}")
    for name, bad in fails.items():
        print(f"  {name:32s}: {len(ok)-len(bad)}/{len(ok)}"
              + (f"   *** {len(bad)} FAILURES ***" if bad else ""))
    print(f"  {'P2 (delta=0 => Z-logicals parent’s)':32s}: {len(p2)-len(p2bad)}/{len(p2)}"
          + (f"   *** {len(p2bad)} FAILURES ***" if p2bad else ""))
    dd = {d: sum(1 for r in ok if r["delta"] == d) for d in sorted({r["delta"] for r in ok})}
    print(f"  delta distribution: {dd}")
    n0 = sum(1 for r in ok if r["delta"] == 0)
    print(f"\n  Corollary 1 settles {n0}/{len(ok)} of these with no distance needed:")
    print(f"    the parent BB code has the same (n,k) and d(parent) >= d(PBB).")
    print(f"  The remaining {len(ok)-n0} have delta>0 and are NOT settled by theory.")

    (OUT / f"exp018_heldout_{lattice}.json").write_text(json.dumps(
        {"stage1_dimension_envelope": s1, "perturbation_spaces": spaces,
         "structural_results": results,
         "verdict": {
             "tested": len(ok),
             "failures": {k2: len(v) for k2, v in fails.items()},
             "p2_tested": len(p2), "p2_failures": len(p2bad),
             "delta_distribution": {str(k2): v for k2, v in dd.items()},
             "settled_by_corollary1_delta_zero": n0,
             "unsettled_delta_positive": len(ok) - n0},
         "caveat": ("stage 1 is a DIMENSION envelope (k only). No rate-distance "
                    "claim is made from it. Distances are computed separately in "
                    "stage 4 for delta>0 codes only.")},
        indent=2, default=str))
    print(f"\nwrote results/processed/exp018_heldout_{lattice}.json")
