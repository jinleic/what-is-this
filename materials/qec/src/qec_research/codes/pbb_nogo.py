"""Parent-level no-go bound for perturbed bivariate-bicycle codes (Theorem H).

Theorem G decides one catalogue row at a time: given a minimum-weight
``Z``-logical of the parent that survives the dressing space, that survivor
certifies ``d_PBB <= d_Z(parent)``.  The cost is one exact parent distance per
row, and the conclusion covers exactly one perturbation ``[C D]``.

This module removes the dependence on ``[C D]`` entirely.

The observation.  ``Delta = { lambda [C D] : lambda [A B] = 0 }`` is not merely
a subspace: the coefficient set ``{lambda : lambda[A B] = 0}`` is an *ideal* of
``R = GF(2)[x,y]/(x^l-1, y^m-1)`` (if ``lambda A = lambda B = 0`` then
``(x lambda) A = x (lambda A) = 0``), so ``Delta`` is an ``R``-submodule of the
``Z``-sector -- it is closed under the translations ``x^a y^b``.  The parent's
``Z``-stabilizer space ``S_Z = rowspace[B^T A^T]`` is an ``R``-submodule too.
Therefore *absorbing one element of a translation orbit absorbs the whole
orbit*, and more generally the absorbed set is always a submodule.

Consequently, define

    M(P) = span of the translation orbits of ALL minimum-weight Z-logicals
           of the parent P,       T(P) = dim (M + S_Z)/S_Z .

``T(P)`` depends only on the parent.  A strict distance increase forces
``Delta`` to absorb every minimum-weight Z-logical, hence to contain ``M``
modulo ``S_Z``, hence ``dim Delta >= T``; with ``k_Q = k_P - dim Delta``
(Theorem G(ii)) this gives the perturbation-independent no-go

    d_Q > d_Z(P)   ==>   k_Q <= k_P - T(P)                       (Theorem H)

for **every** perturbation of ``P``.  In particular ``T(P) = k_P`` closes the
entire family: no PBB over ``P`` keeps a logical qubit while exceeding
``d_Z(P)``.

Computing ``T`` exactly.  Grow ``W = S_Z + M`` by repeatedly asking a SAT
solver for a weight-``<= d_Z`` element of ``ker[A B]`` *outside* ``W``.  Each
SAT hit is a new minimum-weight logical (weight cannot be below ``d_Z`` outside
``S_Z``), whose whole orbit joins ``M``.  The terminating UNSAT certifies that
``W`` contains every minimum-weight Z-logical, so ``T`` is exact rather than a
lower bound.  Nontriviality-with-respect-to-``W`` is encoded exactly like the
usual logical nontriviality: pair against a basis of ``W^perp``, since
``(exists f in W^perp with f.z = 1)`` iff ``z not in W``.

The translation symmetry break stays sound throughout: ``ker[A B]``, ``S_Z``
and every ``M`` built from full orbits are translation invariant, so the
solution set of each query is translation invariant and any solution can be
translated to anchor support at index 0 of one block.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..distance.sat_decide import (
    DecisionInstance,
    decide_weight_bounded,
    verify_witness_two_paths,
)
from ..gf2.linalg import nullspace_np, rank_np
from .pbb_survival import translation_orbit


@dataclass
class NoGoCertificate:
    """Parent-level distance-increase obstruction.

    ``T`` is exact when ``closed`` is true (the enumeration ended in UNSAT);
    otherwise it is a certified *lower* bound, which is the sound direction:
    Theorem H only ever needs ``dim Delta >= T``.
    """

    k_parent: int
    d_z_parent: int
    T: int
    closed: bool
    witnesses: list[np.ndarray] = field(default_factory=list)
    calls: list[dict[str, Any]] = field(default_factory=list)

    @property
    def k_ceiling(self) -> int:
        """Any perturbation with ``d_Q > d_Z(P)`` has ``k_Q <= k_ceiling``."""

        return self.k_parent - self.T

    @property
    def d_z_exact(self) -> bool:
        """``d_z_parent`` is pinned exactly iff a witness was actually found."""

        return bool(self.witnesses)

    @property
    def family_closed(self) -> bool:
        """True when no perturbation can both beat ``d_Z(P)`` and encode."""

        return self.k_ceiling <= 0

    def permits(self, k_pbb: int) -> bool:
        """Whether a perturbation with this ``k`` may exceed ``d_Z(parent)``.

        ``False`` is a proof, for that ``k``, that no such perturbation exists
        -- no search over ``[C D]`` is performed or needed.
        """

        return int(k_pbb) <= self.k_ceiling


def orthogonal_complement(W: np.ndarray, n: int) -> np.ndarray:
    """Basis of ``W^perp`` = {f : f . w = 0 for every row w of W}."""

    W = np.asarray(W, dtype=np.uint8) & 1
    if W.size == 0:
        return np.eye(n, dtype=np.uint8)
    return nullspace_np(W)


def outside_instance(
    HX: np.ndarray, W: np.ndarray, *, block_length: int | None = None
) -> DecisionInstance:
    """``exists z: H_X z = 0, wt(z) <= cap, z not in rowspace(W)``.

    ``W`` must be translation invariant for ``symmetry_clause`` to stay sound;
    callers build it as ``S_Z`` plus full translation orbits.
    """

    HX = np.asarray(HX, dtype=np.uint8) & 1
    n = int(HX.shape[1])
    pairing = orthogonal_complement(W, n)
    if pairing.shape[0] == 0:
        raise ValueError("W spans the whole space: 'outside W' is unsatisfiable")
    symmetry = None
    if block_length is not None:
        if 2 * int(block_length) != n:
            raise ValueError(f"block_length {block_length} incompatible with n={n}")
        symmetry = [0, int(block_length)]
    return DecisionInstance(
        parity_rows=HX,
        pairing_rows=pairing,
        groups=[[j] for j in range(n)],
        kind="css_z",
        meta={"n": n, "block_length": block_length, "modulus_rank": int(rank_np(W))},
        symmetry_clause=symmetry,
    )


def minimum_weight_module(
    HX: np.ndarray,
    HZ: np.ndarray,
    d_z_lower: int,
    ell: int,
    m: int,
    *,
    k_parent: int,
    solver_name: str = "cadical195",
    conflict_budget: int = 0,
    max_witnesses: int = 64,
    max_climb: int = 32,
) -> NoGoCertificate:
    """Close the module generated by all minimum-weight parent ``Z``-logicals.

    ``d_z_lower`` must be a *certified* lower bound on ``d_Z(parent)`` -- the
    exact CSS distance ``min(d_X,d_Z)`` qualifies, since the UNSAT-below proof
    covers the ``Z`` side.  The first query doubles as the exactness test: an
    UNSAT at cap ``c`` certifies ``d_Z > c``, so the cap climbs until the first
    SAT pins ``d_Z`` exactly.  Nothing here assumes ``d_X = d_Z``.

    Every witness is re-verified through two independent GF(2) paths before it
    is allowed to grow ``W``.
    """

    HX = np.asarray(HX, dtype=np.uint8) & 1
    HZ = np.asarray(HZ, dtype=np.uint8) & 1
    base_rank = int(rank_np(HZ))
    W = HZ.copy()
    witnesses: list[np.ndarray] = []
    calls: list[dict[str, Any]] = []
    closed = False
    cap = int(d_z_lower)
    d_z = None
    climbs = 0

    while len(witnesses) < max_witnesses:
        instance = outside_instance(HX, W, block_length=ell * m)
        record = decide_weight_bounded(
            instance, cap, solver_name=solver_name, conflict_budget=conflict_budget
        )
        calls.append(
            {
                "status": record["status"],
                "weight_cap": cap,
                "modulus_rank": int(rank_np(W)),
                "cnf_sha256": record.get("cnf_sha256"),
                "encoding_version": record.get("encoding_version"),
                "solver": record.get("solver"),
                "wall_time_s": record.get("wall_time_s"),
            }
        )
        if record["status"] == "UNSAT":
            if d_z is None:
                # No Z-logical at this cap: d_Z > cap.  Climb; the module is
                # still empty, so nothing computed so far is invalidated.
                climbs += 1
                if climbs > max_climb:
                    break
                cap += 1
                continue
            closed = True
            break
        if record["status"] != "SAT":
            break  # budget exhausted: T below is a certified lower bound
        z = np.asarray(record["vector"], dtype=np.uint8).reshape(-1)
        checked = verify_witness_two_paths(instance, z)
        if not checked["valid"]:
            raise RuntimeError("SAT witness failed independent re-verification")
        weight = int(z.sum())
        if d_z is None:
            # First witness outside S_Z: the preceding UNSATs (at cap-1 from the
            # supplied certified bound, plus any climbs) pin d_Z = weight.
            d_z = weight
            cap = weight
        elif weight != d_z:
            raise RuntimeError(
                f"witness weight {weight} != d_Z {d_z}: cap discipline broken"
            )
        orbit = translation_orbit(z, ell, m)
        if (orbit.sum(axis=1) != weight).any():
            raise RuntimeError("translation changed the weight: wrong index convention")
        witnesses.append(z)
        W = np.vstack([W, orbit])

    T = int(rank_np(W)) - base_rank
    if T > k_parent:
        raise RuntimeError(
            f"T={T} exceeds k_parent={k_parent}: the minimum-weight module cannot "
            "be larger than the logical space"
        )
    return NoGoCertificate(
        k_parent=int(k_parent),
        d_z_parent=int(d_z if d_z is not None else d_z_lower),
        T=T,
        closed=closed,
        witnesses=witnesses,
        calls=calls,
    )
