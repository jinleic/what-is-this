"""Structural domination criterion for perturbed bivariate-bicycle codes.

The whole point: decide ``d_PBB <= d_Z(parent)`` by *linear algebra alone*, with
no weight-bounded search on the PBB, for a large fraction of the catalogue.

Setup.  A parent CSS BB code has ``H_X = [A B]`` and ``H_Z = [B^T A^T]``.  The
PBB built on it is

    H = [ A B | C D ]
        [ 0 0 | B^T A^T ]

Three facts drive everything (each machine-checked in
``tests/test_pbb_survival.py``):

1. **The pure-Z centralizer is unchanged.**  A vector ``(0|z)`` commutes with
   the second row block automatically (Z against Z) and with a first-block row
   ``(a|c)`` iff ``a . z = 0``.  So the pure-Z centralizer of the PBB is
   ``ker[A B]`` -- exactly the parent's.

2. **The pure-Z stabilizer group grows by a dressing space.**  A combination of
   first-block rows is pure Z iff its X part vanishes, i.e. iff its coefficient
   vector lies in the left kernel of ``[A B]``.  Such a combination contributes
   its Z part, so the PBB's pure-Z stabilizers are

       rowspace(H_Z) + Delta,   Delta = { lambda [C D] : lambda [A B] = 0 }.

   Consequently ``k_PBB = k_parent - dim Delta`` (verified on every catalogue
   row, and it gives an independent cross-check on both dimensions).

3. **Survival bounds the PBB distance.**  If a minimum-weight Z-logical ``z`` of
   the parent is *not* in ``rowspace(H_Z) + Delta`` then ``z`` is still a
   nontrivial logical of the PBB, of the same weight.  Hence

       d_PBB <= wt(z) = d_Z(parent).

   Perturbation can only *remove* logical qubits, never lighten a surviving one.

Counting form (the cheap test).  A strict increase ``d_PBB > d_Z(parent)``
requires Delta to kill *every* minimum-weight Z-logical, so Delta must contain
their whole span modulo stabilizers.  If that span has dimension ``t`` then
``dim Delta >= t``, i.e.

    d_PBB > d_Z(parent)  ==>  k_PBB <= k_parent - t.

So exhibiting ``dim Delta + 1`` independent minimum-weight Z-logicals of the
parent certifies ``d_PBB <= d_Z(parent)`` outright.  BB translations are code
automorphisms, so one witness generates its whole orbit for free.

This explains the empirical pattern: reversals only ever occur where the
perturbation halves ``k`` (``dim Delta = k_parent/2``), because only then is
Delta big enough to absorb the minimum-weight span.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..gf2.linalg import nullspace_np, rank_np


def dressing_space(HX: np.ndarray, CD: np.ndarray) -> np.ndarray:
    """``Delta = { lambda [C D] : lambda [A B] = 0 }`` as a row matrix.

    ``HX`` is the parent's ``[A B]``; ``CD`` is the Z-part of the PBB's dressed
    first row block.  Rows of the result span Delta (not necessarily reduced).
    """

    left = nullspace_np(np.ascontiguousarray(HX.T))
    if left.size == 0:
        return np.zeros((0, CD.shape[1]), dtype=np.uint8)
    return (left.astype(np.uint8) @ CD.astype(np.uint8)) % 2


def translation_orbit(
    vector: np.ndarray, ell: int, m: int, *, blocks: int = 2
) -> np.ndarray:
    """All ``ell*m`` translates of a Z-sector vector under ``x^a y^b``.

    Qubits are indexed block-major with position ``i*m + j`` inside each block,
    matching :func:`qec_research.codes.bicycle.poly_matrix`.  Translations act
    the same way on every block, commute with the circulant blocks of ``H_X``
    and ``H_Z``, and therefore map logicals to logicals of equal weight.
    """

    v = np.asarray(vector, dtype=np.uint8).reshape(blocks, ell, m)
    out = np.empty((ell * m, blocks * ell * m), dtype=np.uint8)
    row = 0
    for a in range(ell):
        rolled_a = np.roll(v, a, axis=1)
        for b in range(m):
            out[row] = np.roll(rolled_a, b, axis=2).reshape(-1)
            row += 1
    return out


def quotient_rank(vectors: np.ndarray, stabilizers: np.ndarray) -> int:
    """``dim span(vectors) mod rowspace(stabilizers)``."""

    base = rank_np(stabilizers) if stabilizers.size else 0
    if vectors.size == 0:
        return 0
    stack = np.vstack([stabilizers, vectors]) if stabilizers.size else vectors
    return rank_np(stack) - base


@dataclass
class SurvivalVerdict:
    """Outcome of the structural test on one catalogue row."""

    dim_delta: int
    orbit_rank: int
    k_parent: int
    k_pbb: int
    dimension_identity_holds: bool
    survivor: np.ndarray | None = None
    survivor_weight: int | None = None
    detail: dict[str, Any] = field(default_factory=dict)

    @property
    def dominated(self) -> bool:
        """True when ``d_PBB <= wt(survivor) = d_Z(parent)`` is certified."""

        return self.survivor is not None

    @property
    def counting_margin(self) -> int:
        """``orbit_rank - dim_delta``; positive forces a survivor to exist."""

        return self.orbit_rank - self.dim_delta


def survival_test(
    HX: np.ndarray,
    HZ: np.ndarray,
    CD: np.ndarray,
    witness: np.ndarray,
    ell: int,
    m: int,
    *,
    k_parent: int | None = None,
    k_pbb: int | None = None,
) -> SurvivalVerdict:
    """Decide structural domination from one minimum-weight parent Z-logical.

    ``witness`` must be a minimum-weight Z-logical of the parent: in
    ``ker[A B]``, outside ``rowspace(H_Z)``, of weight ``d_Z(parent)``.  Both
    conditions are re-checked here rather than trusted.

    Returns a verdict carrying an explicit surviving logical whenever one is
    found; that survivor *is* the certificate for ``d_PBB <= d_Z(parent)``.
    """

    HX = np.asarray(HX, dtype=np.uint8)
    HZ = np.asarray(HZ, dtype=np.uint8)
    CD = np.asarray(CD, dtype=np.uint8)
    w = np.asarray(witness, dtype=np.uint8).reshape(-1)

    if ((HX.astype(np.int64) @ w.astype(np.int64)) % 2).any():
        raise ValueError("witness is not in ker[A B]: not a Z-centralizer element")
    if quotient_rank(w.reshape(1, -1), HZ) == 0:
        raise ValueError("witness is a parent stabilizer, not a logical")

    delta = dressing_space(HX, CD)
    pbb_stab = np.vstack([HZ, delta]) if delta.size else HZ
    dim_delta = rank_np(pbb_stab) - rank_np(HZ)

    orbit = translation_orbit(w, ell, m)
    # Translations preserve weight and the logical quotient, so these are
    # guards against a wrong index convention, not filters.
    weights = orbit.sum(axis=1)
    if (weights != weights[0]).any():
        raise ValueError("translation changed the weight: index convention is wrong")
    orbit_rank = quotient_rank(orbit, HZ)

    survivor = None
    survivor_weight = None
    for candidate in orbit:
        if quotient_rank(candidate.reshape(1, -1), pbb_stab) > 0:
            survivor = candidate.copy()
            survivor_weight = int(candidate.sum())
            break

    identity_ok = True
    if k_parent is not None and k_pbb is not None:
        identity_ok = (k_parent - k_pbb) == dim_delta

    return SurvivalVerdict(
        dim_delta=int(dim_delta),
        orbit_rank=int(orbit_rank),
        k_parent=int(k_parent) if k_parent is not None else -1,
        k_pbb=int(k_pbb) if k_pbb is not None else -1,
        dimension_identity_holds=bool(identity_ok),
        survivor=survivor,
        survivor_weight=survivor_weight,
        detail={
            "orbit_size": int(orbit.shape[0]),
            "witness_weight": int(w.sum()),
            "distinct_orbit_weights": sorted({int(x) for x in weights}),
        },
    )
