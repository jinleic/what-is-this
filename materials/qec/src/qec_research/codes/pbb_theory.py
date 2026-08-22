"""Structural analysis of perturbed bivariate-bicycle (PBB) codes.

Main object: the relationship between a PBB code (A,B,C,D) and its *parent*
CSS BB code (A,B).

Notation
--------
    dim = ell*m,  n = 2*dim
    HX  = [A B]                (dim x n)   -- x-part of block 1
    HZ  = [B^T A^T]            (dim x n)   -- z-part of block 2
    P   = [C D]                (dim x n)   -- z-part of block 1 (perturbation)

    H_PBB = [ HX | P  ]
            [ 0  | HZ ]

Proposition 1 (pure-Z centralizer is unchanged).
    A pure-Z operator (0|w) commutes with every PBB generator iff HX w^T = 0.
    Proof: <(0|w),(HX_i|P_i)> = w . HX_i ; <(0|w),(0|HZ_i)> = 0.  The pure-Z
    centralizer of the PBB code therefore equals ker(HX), i.e. exactly the
    pure-Z centralizer of the parent CSS BB code.  []

Proposition 2 (pure-Z stabilizer subgroup).
    The pure-Z elements of rowspace(H_PBB) are
        S_Z = { u P + w HZ : u HX = 0 }
            = rowspace(HZ) + { u P : u in leftnull(HX) }.
    Write delta = dim( S_Z ) - rank(HZ)  in [0, dim - rank(HX)].  []

Proposition 3 (dimension bookkeeping).
    rank(H_PBB) = rank(HX) + dim(S_Z) = rank(HX) + rank(HZ) + delta,
    hence
        k_PBB = k_BB - delta.                                        []

THEOREM 1 (PBB distance ceiling).
    If k_PBB = k_BB (equivalently delta = 0) then the nontrivial pure-Z
    logical operators of the PBB code are *exactly* those of the parent CSS
    BB code.  Consequently

        d(PBB(A,B,C,D))  <=  d_Z(BB(A,B)) .

    Proof.  delta = 0 means S_Z = rowspace(HZ).  By Prop. 1 the pure-Z
    centralizer is ker(HX) for both codes, and by Prop. 2 the pure-Z
    stabilizers coincide.  So the two codes have literally the same set of
    nontrivial pure-Z logical operators.  A minimum-weight one has weight
    d_Z(BB(A,B)) and is a valid (symplectic-weight equal) nontrivial logical
    of the PBB code, so the PBB distance cannot exceed it.  []

COROLLARY 1.  For any BB code with d_X = d_Z = d (true for every BB instance
    in the Bravyi catalogue by the A<->B transpose symmetry), no k-preserving
    perturbation can raise the distance:  d(PBB) <= d(BB).

COROLLARY 2 (why perturbation cannot buy rate either).
    delta > 0 strictly *lowers* k.  So a PBB code either keeps k and is
    capped in distance by its parent, or trades k away.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from ..gf2.linalg import matmul, nullspace_np, rank_np, rref_np
from .bicycle import PBBSpec, poly_matrix

__all__ = ["PBBStructure", "analyse_pbb", "parent_bb_matrices", "css_shadow",
           "bb_transpose_swap_permutation"]


def parent_bb_matrices(spec: PBBSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (HX, HZ, P) = ([A B], [B^T A^T], [C D])."""
    ell, m = spec.ell, spec.m
    dim = ell * m
    A = poly_matrix(ell, m, spec.A)
    B = poly_matrix(ell, m, spec.B)
    C = poly_matrix(ell, m, spec.C) if spec.C else np.zeros((dim, dim), np.uint8)
    D = poly_matrix(ell, m, spec.D) if spec.D else np.zeros((dim, dim), np.uint8)
    HX = np.hstack([A, B]).astype(np.uint8)
    HZ = np.hstack([B.T, A.T]).astype(np.uint8)
    P = np.hstack([C, D]).astype(np.uint8)
    return HX, HZ, P


def _left_nullspace(M: np.ndarray) -> np.ndarray:
    """Basis of {u : u M = 0}."""
    return nullspace_np(np.asarray(M, dtype=np.uint8).T)


@dataclass
class PBBStructure:
    n: int
    dim: int
    rank_HX: int
    rank_HZ: int
    k_bb: int
    k_pbb: int
    delta: int
    z_logicals_identical: bool
    theorem1_applies: bool

    def to_dict(self) -> dict:
        return asdict(self)


def analyse_pbb(spec: PBBSpec) -> PBBStructure:
    HX, HZ, P = parent_bb_matrices(spec)
    dim = spec.ell * spec.m
    n = 2 * dim
    rX = rank_np(HX)
    rZ = rank_np(HZ)
    k_bb = n - rX - rZ

    U = _left_nullspace(HX)                     # dim(U) = dim - rX
    if U.shape[0]:
        extra = matmul(U, P)                    # candidate new pure-Z stabilizers
        stacked = np.vstack([HZ, extra])
    else:
        stacked = HZ
    dim_SZ = rank_np(stacked)
    delta = dim_SZ - rZ
    k_pbb = n - rX - dim_SZ
    return PBBStructure(
        n=n, dim=dim, rank_HX=rX, rank_HZ=rZ,
        k_bb=k_bb, k_pbb=k_pbb, delta=int(delta),
        z_logicals_identical=(delta == 0),
        theorem1_applies=(delta == 0),
    )


def css_shadow(spec: PBBSpec) -> tuple[np.ndarray, np.ndarray]:
    """The CSS shadow Q' of a PBB code Q.

    THEOREM 2 (CSS shadow).  Let Q = PBB(A,B,C,D) and define the CSS code Q'
    by
            X-checks:  HX  = [A B]
            Z-checks:  S_Z = rowspace(HZ) + { u [C D] : u HX = 0 }
    (S_Z is exactly the pure-Z subgroup of Q's stabilizer group, so
    HX S_Z^T = 0 and Q' is a legitimate CSS code on the same n qubits).
    Then
        (a) n(Q') = n(Q)  and  k(Q') = k(Q);
        (b) the nontrivial pure-Z logical operators of Q are *exactly* the
            Z-logical operators of Q';
        (c) d(Q) <= d_Z(Q').

    Proof.
    (a) k(Q') = n - rank(HX) - rank(S_Z) = k(Q) by Proposition 3.
    (b) By Proposition 1 the pure-Z centralizer of Q is ker(HX), which is the
        Z-centralizer of Q'; by Proposition 2 the pure-Z stabilizers of Q are
        S_Z, the Z-stabilizers of Q'.  The quotients therefore coincide.
    (c) Take a minimum-weight Z-logical w of Q'.  By (b), (0|w) is a
        nontrivial logical of Q with symplectic weight |w| = d_Z(Q').  []

    COROLLARY (domination).  If additionally d_X(Q') >= d(Q) then
        d(Q')  =  min(d_X(Q'), d_Z(Q'))  >=  d(Q),
    so the CSS code Q' matches Q in n and k and is no worse in distance --
    while requiring only pure-X and pure-Z checks.  Any circuit-level
    advantage claimed for Q must then come from somewhere other than its
    code parameters.
    """
    HX, HZ, P = parent_bb_matrices(spec)
    U = _left_nullspace(HX)
    SZ = np.vstack([HZ, matmul(U, P)]).astype(np.uint8) if U.shape[0] else HZ
    return HX, SZ


def bb_transpose_swap_permutation(ell: int, m: int) -> np.ndarray:
    """The involution sigma on the n = 2*ell*m qubits that exchanges the two
    blocks and inverts the lattice:  (L,(i,j)) <-> (R,(-i,-j)).

    PROPOSITION 4.  For any BB code,  sigma maps H_X = [A B] onto H_Z =
    [B^T A^T] and vice versa.  Hence d_X = d_Z for every bivariate-bicycle
    code.  (Swapping the halves sends [A B] -> [B A]; the lattice inversion
    then transposes each circulant, giving [B^T A^T].)
    """
    dim = ell * m
    inv = np.empty(dim, dtype=int)
    for i in range(ell):
        for j in range(m):
            inv[i * m + j] = ((-i) % ell) * m + ((-j) % m)
    return np.concatenate([inv + dim, inv])
