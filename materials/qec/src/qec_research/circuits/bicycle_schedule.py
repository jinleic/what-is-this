"""Translation-symmetry-reduced scheduling for BB / PBB codes.

Both families are invariant under the abelian group G = Z_ell x Z_m acting
simultaneously on checks and on qubits.  We therefore look for a
*translation-invariant* schedule: the time slot of the gate between a check
and a data qubit depends only on

    (check block, qubit half, group displacement)

i.e. on the monomial that generated the incidence, not on which translate of
the check we are looking at.  This collapses the CP-SAT model from ~10^3
integer variables to one per monomial (16 for a weight-6 PBB with |C|=|D|=2),
and it is exactly the structure of the hand-built depth-7 BB schedule.

Also provides pure-Z logical representatives, which are the only observables a
Z-basis memory experiment can read out.
"""

from __future__ import annotations

import numpy as np

from ..codes.bicycle import BBSpec, PBBSpec, bb_stabilizer, build_pbb
from ..codes.pbb_theory import css_shadow, parent_bb_matrices
from ..gf2.linalg import matmul, nullspace_np, rank_np, rref_np
from ..symplectic.core import StabilizerCode
from .mixed_stabilizer import PauliSupport, generator_supports

__all__ = ["translation_orbits", "pure_z_logical_basis", "pbb_supports_and_orbits",
           "bb_supports_and_orbits"]


def _group_of(index: int, m: int) -> tuple[int, int]:
    return (index // m, index % m)


def translation_orbits(H: np.ndarray, ell: int, m: int) -> dict[tuple[int, int], int]:
    """Map each (check, qubit) edge to its translation orbit.

    Checks are indexed block-major: rows [0,dim) are block 1, [dim,2dim) are
    block 2.  Qubits are [0,dim) = left half, [dim,2dim) = right half.
    Orbit key = (check_block, qubit_half, displacement in Z_ell x Z_m).
    """
    dim = ell * m
    sup = generator_supports(H)
    ids: dict[tuple, int] = {}
    out: dict[tuple[int, int], int] = {}
    for s in sup:
        cb, cg = divmod(s.index, dim)
        ci, cj = _group_of(cg, m)
        for q, pauli in s.paulis.items():
            qh, qg = divmod(q, dim)
            qi, qj = _group_of(qg, m)
            key = (cb, qh, (qi - ci) % ell, (qj - cj) % m, pauli)
            if key not in ids:
                ids[key] = len(ids)
            out[(s.index, q)] = ids[key]
    return out


def pure_z_logical_basis(code: StabilizerCode) -> np.ndarray:
    """Basis of the pure-Z logical classes, as full (2n) symplectic vectors.

    Pure-Z centraliser = ker(H_x-part);  pure-Z stabilisers = the pure-Z
    subgroup of rowspace(H).  Returns one representative per quotient
    dimension.
    """
    n = code.n
    H = code.H
    HX = H[:, :n]
    U = nullspace_np(HX.T)                       # u with u @ HX = 0
    SZ = matmul(U, H[:, n:]) if U.shape[0] else np.zeros((0, n), np.uint8)
    ker = nullspace_np(HX)
    R, _ = rref_np(SZ) if SZ.shape[0] else (np.zeros((0, n), np.uint8), [])
    cur = [r for r in R]
    r0 = len(cur)
    reps = []
    for row in ker:
        if rank_np(np.array(cur + [row], dtype=np.uint8)) > r0:
            reps.append(row)
            cur.append(row)
            r0 += 1
    out = np.zeros((len(reps), 2 * n), dtype=np.uint8)
    for i, w in enumerate(reps):
        out[i, n:] = w
    return out


def bb_supports_and_orbits(spec: BBSpec):
    code = bb_stabilizer(spec)
    sup = generator_supports(code.H)
    orb = translation_orbits(code.H, spec.ell, spec.m)
    return code, sup, orb


def pbb_supports_and_orbits(spec: PBBSpec):
    code = build_pbb(spec)
    sup = generator_supports(code.H)
    orb = translation_orbits(code.H, spec.ell, spec.m)
    return code, sup, orb
