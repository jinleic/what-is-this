"""Certified exact distance for stabilizer codes.

Key identity used throughout
----------------------------
Let S be the stabilizer group with check matrix H, and let {L_1..L_{2k}} be any
basis of the logical quotient S^perp / S.  A vector v in S^perp is a
*nontrivial* logical iff <v, L_j>_s = 1 for at least one j.  Hence

    d = min_{j=1..2k}  min { wt_s(v) : v in S^perp,  <v,L_j>_s = 1 }.

Each inner problem is a bounded integer program that CP-SAT solves to proven
optimality, so the overall minimum is *exact* (not an upper bound), provided
every subproblem reaches OPTIMAL.  If some subproblem only reaches a feasible
solution within the time limit we return a certified upper bound plus the best
lower bound CP-SAT proved.

Parity constraints  sum_j a_j v_j = b (mod 2)  are linearised as
    sum_j a_j v_j - 2 t = b ,   t integer in [0, floor(|a|/2)] .
This is exact integer arithmetic; no floating point enters the model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Iterable, Sequence

import numpy as np
from ortools.sat.python import cp_model

from ..gf2.linalg import rows_to_bitsets
from ..symplectic.core import StabilizerCode, lambda_swap, symplectic_weight

__all__ = ["DistanceResult", "exact_distance_symplectic", "exact_distance_css",
           "min_weight_with_parity"]


@dataclass
class DistanceResult:
    value: int | None            # best objective found (None if nothing found)
    lower_bound: int             # proven lower bound on the true distance
    exact: bool                  # True iff value == true distance is proven
    status: str
    witness: list[int] | None    # support (indices) of the minimising operator
    witness_vector: list[int] | None
    per_sector: list[dict]
    wall_time_s: float

    def to_dict(self) -> dict:
        return asdict(self)


def _add_parity(model: cp_model.CpModel, lits: Sequence, coeffs: Sequence[int],
                rhs: int, tag: str) -> None:
    """sum coeffs_j * lits_j == rhs (mod 2)."""
    idx = [i for i, c in enumerate(coeffs) if c & 1]
    if not idx:
        if rhs & 1:
            model.add_bool_or([])  # infeasible
        return
    t = model.new_int_var(0, len(idx) // 2, f"t_{tag}")
    model.add(sum(lits[i] for i in idx) - 2 * t == rhs & 1)


def min_weight_with_parity(
    constraints: np.ndarray,
    objective_groups: list[list[int]],
    parity_row: np.ndarray,
    nvars: int,
    *,
    time_limit_s: float = 300.0,
    workers: int = 8,
    upper_bound: int | None = None,
) -> tuple[int | None, int, str, list[int] | None]:
    """Minimise |{g : some var in g is 1}| subject to
       constraints @ v = 0 (mod 2)  and  parity_row . v = 1 (mod 2).

    ``objective_groups[g]`` lists variable indices whose OR forms weight unit g.
    Returns (objective, proven_lower_bound, status, solution).
    """
    model = cp_model.CpModel()
    v = [model.new_bool_var(f"v{j}") for j in range(nvars)]

    for i, row in enumerate(np.asarray(constraints, dtype=np.uint8)):
        _add_parity(model, v, row.tolist(), 0, f"c{i}")
    _add_parity(model, v, np.asarray(parity_row, dtype=np.uint8).tolist(), 1, "obj")

    w = []
    for g, grp in enumerate(objective_groups):
        if len(grp) == 1:
            w.append(v[grp[0]])
        else:
            wg = model.new_bool_var(f"w{g}")
            for j in grp:
                model.add_implication(v[j], wg)
            model.add_bool_or([v[j] for j in grp] + [wg.negated()])
            w.append(wg)
    model.minimize(sum(w))
    if upper_bound is not None:
        model.add(sum(w) <= upper_bound)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = workers
    st = solver.solve(model)
    name = solver.status_name(st)
    if st == cp_model.OPTIMAL:
        sol = [j for j in range(nvars) if solver.value(v[j])]
        val = int(round(solver.objective_value))
        return val, val, name, sol
    if st == cp_model.FEASIBLE:
        sol = [j for j in range(nvars) if solver.value(v[j])]
        return int(round(solver.objective_value)), int(math.ceil(solver.best_objective_bound)), name, sol
    if st == cp_model.INFEASIBLE:
        return None, 10**9, name, None
    return None, 0, name, None


def exact_distance_symplectic(
    code: StabilizerCode,
    *,
    time_limit_s: float = 300.0,
    workers: int = 8,
    logicals: np.ndarray | None = None,
    upper_bound: int | None = None,
    sector_filter: Iterable[int] | None = None,
) -> DistanceResult:
    """Exact minimum symplectic weight of a nontrivial logical operator.

    Works for arbitrary (CSS or non-CSS) stabilizer codes.
    """
    import time

    t0 = time.time()
    n = code.n
    H = code.H
    if logicals is None:
        logicals = code.logical_basis()
    if logicals.shape[0] == 0:
        return DistanceResult(None, 10**9, True, "K_ZERO", None, None, [], 0.0)

    # v ranges over GF(2)^{2n}; membership in S^perp is  (Lambda H) v = 0.
    constraints = lambda_swap(H)
    # symplectic weight groups: qubit j -> {x_j, z_j}
    groups = [[j, n + j] for j in range(n)]

    best = upper_bound
    best_sol: list[int] | None = None
    best_lb = 0
    per_sector: list[dict] = []
    all_exact = True
    sectors = range(logicals.shape[0]) if sector_filter is None else list(sector_filter)

    cap0 = best            # `upper_bound` is a SEARCH CAP, not an achieved value
    best = None
    for j in sectors:
        # <v, L_j>_s = v . (Lambda L_j)
        prow = lambda_swap(logicals[j][None, :])[0]
        ub = best if best is not None else cap0
        if ub is not None and ub < 1:
            per_sector.append({"sector": int(j), "status": "PRUNED", "value": None})
            continue
        val, lb, status, sol = min_weight_with_parity(
            constraints, groups, prow, 2 * n,
            time_limit_s=time_limit_s, workers=workers, upper_bound=ub,
        )
        per_sector.append({"sector": int(j), "status": status, "value": val, "lb": int(lb)})
        if status not in ("OPTIMAL", "INFEASIBLE"):
            all_exact = False
        if val is not None and (best is None or val < best):
            best, best_sol = val, sol

    witness_vec = None
    support = None
    if best_sol is not None:
        vv = np.zeros(2 * n, dtype=np.uint8)
        vv[best_sol] = 1
        witness_vec = vv.tolist()
        support = sorted({j % n for j in best_sol})
    # exact iff every subproblem decided AND a witness of exactly that weight
    # is in hand; "all sectors infeasible below the cap" proves only a bound.
    exact = bool(all_exact and best is not None and witness_vec is not None
                 and symplectic_weight(np.asarray(witness_vec, dtype=np.uint8)) == best)
    lower = best if exact else ((cap0 + 1) if (all_exact and cap0 and best is None) else 0)
    return DistanceResult(
        value=best, lower_bound=lower, exact=exact,
        status="OPTIMAL" if exact else ("BOUND" if all_exact else "PARTIAL"),
        witness=support, witness_vector=witness_vec,
        per_sector=per_sector, wall_time_s=time.time() - t0,
    )


def exact_distance_css(
    HX: np.ndarray, HZ: np.ndarray, *,
    time_limit_s: float = 300.0, workers: int = 8,
    upper_bound: int | None = None,
) -> dict:
    """Exact d_X and d_Z for a CSS code with H_X H_Z^T = 0.

    d_X = min{|v| : H_Z v = 0, v not in rowspace(H_X)}   (an X-type logical)
    d_Z = min{|v| : H_X v = 0, v not in rowspace(H_Z)}
    d   = min(d_X, d_Z)
    """
    import time
    from ..gf2.linalg import nullspace_np, rank_np, rref_np

    t0 = time.time()
    HX = np.asarray(HX, dtype=np.uint8) & 1
    HZ = np.asarray(HZ, dtype=np.uint8) & 1
    n = HX.shape[1]
    out: dict = {}

    def quotient_basis(Hc: np.ndarray, Hs: np.ndarray) -> list[np.ndarray]:
        """Basis of ker(Hc) / rowspace(Hs)."""
        ker = nullspace_np(Hc)
        Rs, _ = rref_np(Hs)
        reps: list[np.ndarray] = []
        cur = [r for r in Rs]
        r0 = len(cur)
        for row in ker:
            if rank_np(np.array(cur + [row], dtype=np.uint8)) > r0:
                reps.append(row)
                cur.append(row)
                r0 += 1
        return reps

    # X-type logicals live in ker(H_Z)/rowspace(H_X); Z-type in ker(H_X)/rowspace(H_Z).
    LX = quotient_basis(HZ, HX)
    LZ = quotient_basis(HX, HZ)

    # An X-type v is nontrivial iff it pairs to 1 with some *Z*-type logical.
    for tag, Hc, detectors in (("X", HZ, LZ), ("Z", HX, LX)):
        reps = detectors
        # `upper_bound` is a SEARCH CAP, never an achieved value.  `best` stays
        # None until an actual minimising operator is produced, so a run in
        # which every sector is infeasible below the cap can only ever report
        # a lower bound -- it must not silently return the cap as the distance.
        cap0 = upper_bound
        best: int | None = None
        best_sol = None
        all_decided = True
        sect = []
        groups = [[j] for j in range(n)]
        for j, rep in enumerate(reps):
            ub = best if best is not None else cap0
            if ub is not None and ub < 1:
                sect.append({"sector": j, "status": "PRUNED"})
                continue
            val, lb, status, sol = min_weight_with_parity(
                Hc, groups, rep, n, time_limit_s=time_limit_s,
                workers=workers, upper_bound=ub)
            sect.append({"sector": j, "status": status, "value": val})
            if status not in ("OPTIMAL", "INFEASIBLE"):
                all_decided = False
            if val is not None and (best is None or val < best):
                best, best_sol = val, sol
        # exact iff every subproblem was decided AND we hold a witness whose
        # weight equals the reported value
        exact = bool(all_decided and best is not None and best_sol is not None
                     and len(best_sol) == best)
        out[f"d_{tag}"] = best
        out[f"d_{tag}_exact"] = exact
        out[f"d_{tag}_witness"] = best_sol
        out[f"d_{tag}_all_sectors_decided"] = all_decided
        # A found `best` is only an UPPER bound while any sector is undecided:
        # a lighter operator could still live in an undecided sector.
        if all_decided:
            lb_tag = best if best is not None else ((cap0 + 1) if cap0 else 0)
        else:
            lb_tag = 0
        out[f"d_{tag}_lower_bound"] = lb_tag
        out[f"d_{tag}_is_upper_bound_only"] = bool(best is not None and not all_decided)
        out[f"d_{tag}_search_cap"] = cap0
        out[f"d_{tag}_sectors"] = sect
    dv = [out["d_X"], out["d_Z"]]
    dv = [x for x in dv if x is not None]
    out["d"] = min(dv) if dv else None
    out["d_exact"] = bool(out["d_X_exact"] and out["d_Z_exact"])
    out["d_lower_bound"] = min(out["d_X_lower_bound"], out["d_Z_lower_bound"])
    out["wall_time_s"] = time.time() - t0
    return out
