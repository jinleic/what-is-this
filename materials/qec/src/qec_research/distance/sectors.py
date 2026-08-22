"""Minimum-weight logical operators restricted to a Pauli sector.

CORRECTNESS NOTE (this cost us one falsified claim, see notes/failed_routes.md).
A vector v in the centralizer is a *nontrivial* logical iff

        <v, L_j>_s = 1   for some logical basis element L_j,

and  <v, L>_s = v_x . L_z + v_z . L_x.  Hence for a **pure-Z** operator
(v_x = 0) the detector is the **x-part** of L_j, NOT a Z-type operator.
Pairing two Z-type operators with a plain dot product is degenerate and
admits stabilizer elements as false witnesses.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..gf2.linalg import matmul, nullspace_np, rank_np
from ..symplectic.core import StabilizerCode, lambda_swap
from .exact import min_weight_with_parity

__all__ = ["SectorResult", "min_pure_z_logical", "min_pure_x_logical"]


@dataclass
class SectorResult:
    weight: int | None
    exact: bool
    support: list[int] | None
    n_sectors: int
    statuses: list[str]


def _min_pure_sector(code: StabilizerCode, *, pure: str,
                     logicals: np.ndarray | None, time_limit_s: float,
                     workers: int, upper_bound: int | None) -> SectorResult:
    n = code.n
    H = code.H
    if logicals is None:
        logicals = code.logical_basis()
    if pure == "z":
        # v = (0|w): centralizer condition  H_x w^T = 0 ; detector = L_x
        constraints = H[:, :n]
        detectors = logicals[:, :n]
    elif pure == "x":
        # v = (w|0): centralizer condition  H_z w^T = 0 ; detector = L_z
        constraints = H[:, n:]
        detectors = logicals[:, n:]
    else:
        raise ValueError(pure)

    groups = [[j] for j in range(n)]
    # `upper_bound` is a SEARCH CAP.  `best` must stay None until an actual
    # minimising operator is produced, otherwise a run in which every sector is
    # infeasible below the cap would return the cap as though it were the
    # distance -- proving only a lower bound while claiming exactness.
    cap0 = upper_bound
    best: int | None = None
    best_sol = None
    all_decided = True
    statuses: list[str] = []
    for row in detectors:
        if not row.any():
            statuses.append("TRIVIAL_DETECTOR")
            continue
        cap = best if best is not None else cap0
        if cap is not None and cap < 1:
            statuses.append("PRUNED")
            continue
        val, lb, status, sol = min_weight_with_parity(
            constraints, groups, row, n, time_limit_s=time_limit_s,
            workers=workers, upper_bound=cap)
        statuses.append(status)
        if status not in ("OPTIMAL", "INFEASIBLE"):
            all_decided = False
        if val is not None and (best is None or val < best):
            best, best_sol = val, sol
    exact = bool(all_decided and best is not None and best_sol is not None
                 and len(best_sol) == best)
    return SectorResult(best, exact, best_sol, len(detectors), statuses)


def min_pure_z_logical(code: StabilizerCode, *, logicals: np.ndarray | None = None,
                       time_limit_s: float = 120.0, workers: int = 12,
                       upper_bound: int | None = None) -> SectorResult:
    return _min_pure_sector(code, pure="z", logicals=logicals,
                            time_limit_s=time_limit_s, workers=workers,
                            upper_bound=upper_bound)


def min_pure_x_logical(code: StabilizerCode, *, logicals: np.ndarray | None = None,
                       time_limit_s: float = 120.0, workers: int = 12,
                       upper_bound: int | None = None) -> SectorResult:
    return _min_pure_sector(code, pure="x", logicals=logicals,
                            time_limit_s=time_limit_s, workers=workers,
                            upper_bound=upper_bound)
