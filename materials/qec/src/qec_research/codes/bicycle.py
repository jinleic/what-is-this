"""Independent construction of bivariate-bicycle (BB) and perturbed
bivariate-bicycle (PBB) codes directly from polynomial shift operators.

Ring:   R = GF(2)[x,y]/(x^ell - 1, y^m - 1)
Qubits: index (i,j) in Z_ell x Z_m, flattened idx = i*m + j.

A monomial x^a y^b acts as the permutation (i,j) -> (i+a, j+b).  We represent
it by the matrix  P[idx(i,j), idx(i+a,j+b)] = 1, so that (row = source).
Transposition corresponds to inversion:  (x^a y^b)^T = x^{-a} y^{-b}.

This module intentionally shares NO code with qLDPC so that it can serve as an
independent validation path.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..gf2.linalg import matmul
from ..symplectic.core import StabilizerCode

Terms = list[tuple[int, int]]

__all__ = ["monomial_matrix", "poly_matrix", "BBSpec", "PBBSpec",
           "build_bb", "build_pbb", "poly_transpose_terms", "commutation_defect"]


def monomial_matrix(ell: int, m: int, a: int, b: int) -> np.ndarray:
    dim = ell * m
    P = np.zeros((dim, dim), dtype=np.uint8)
    for i in range(ell):
        for j in range(m):
            P[i * m + j, ((i + a) % ell) * m + ((j + b) % m)] = 1
    return P


def poly_matrix(ell: int, m: int, terms: Terms) -> np.ndarray:
    dim = ell * m
    M = np.zeros((dim, dim), dtype=np.uint8)
    for a, b in terms:
        M ^= monomial_matrix(ell, m, a % ell, b % m)
    return M


def poly_transpose_terms(ell: int, m: int, terms: Terms) -> Terms:
    """Terms of the polynomial whose matrix is the transpose."""
    return [((-a) % ell, (-b) % m) for a, b in terms]


def commutation_defect(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray) -> np.ndarray:
    """M + M^T where M = A C^T + B D^T.  Zero <=> PBB rows commute."""
    M = (matmul(A, C.T) ^ matmul(B, D.T)).astype(np.uint8)
    return (M ^ M.T).astype(np.uint8)


@dataclass
class BBSpec:
    ell: int
    m: int
    A: Terms
    B: Terms
    name: str = ""

    @property
    def n(self) -> int:
        return 2 * self.ell * self.m


@dataclass
class PBBSpec:
    ell: int
    m: int
    A: Terms
    B: Terms
    C: Terms = field(default_factory=list)
    D: Terms = field(default_factory=list)
    name: str = ""

    @property
    def n(self) -> int:
        return 2 * self.ell * self.m

    @property
    def is_css_by_construction(self) -> bool:
        return not self.C and not self.D


def build_bb(spec: BBSpec) -> tuple[np.ndarray, np.ndarray]:
    """Return (H_X, H_Z) for the CSS BB code:  H_X=[A B], H_Z=[B^T A^T]."""
    A = poly_matrix(spec.ell, spec.m, spec.A)
    B = poly_matrix(spec.ell, spec.m, spec.B)
    HX = np.hstack([A, B])
    HZ = np.hstack([B.T, A.T])
    return HX.astype(np.uint8), HZ.astype(np.uint8)


def bb_stabilizer(spec: BBSpec) -> StabilizerCode:
    HX, HZ = build_bb(spec)
    dim = spec.ell * spec.m
    n = 2 * dim
    Z = np.zeros((dim, n), dtype=np.uint8)
    top = np.hstack([HX, Z])
    bot = np.hstack([Z, HZ])
    return StabilizerCode(np.vstack([top, bot]), name=spec.name or "BB")


def build_pbb(spec: PBBSpec, *, check: bool = True) -> StabilizerCode:
    """Build the PBB symplectic check matrix

        H = [ A  B | C    D   ]
            [ 0  0 | B^T  A^T ]
    """
    ell, m = spec.ell, spec.m
    dim = ell * m
    A = poly_matrix(ell, m, spec.A)
    B = poly_matrix(ell, m, spec.B)
    C = poly_matrix(ell, m, spec.C) if spec.C else np.zeros((dim, dim), np.uint8)
    D = poly_matrix(ell, m, spec.D) if spec.D else np.zeros((dim, dim), np.uint8)
    if check:
        defect = commutation_defect(A, B, C, D)
        if defect.any():
            raise ValueError(
                f"PBB rows do not commute: A C^T + B D^T is not symmetric "
                f"({int(defect.sum())} violated entries)"
            )
    zero = np.zeros((dim, dim), dtype=np.uint8)
    top = np.hstack([A, B, C, D])
    bot = np.hstack([zero, zero, B.T, A.T])
    H = np.vstack([top, bot]).astype(np.uint8)
    return StabilizerCode(H, name=spec.name or "PBB")


# --------------------------------------------------------------------------
# Reference instances from the primary artifacts.
# Source: github.com/sbravyi/BivariateBicycleCodes  decoder_setup.py
#         (commit fa77e3333d3ec44c79d8f914dd24c040d1da471b)
# --------------------------------------------------------------------------
BRAVYI_BB: dict[str, BBSpec] = {
    # [[72,12,6]]  ell=6, m=6,  A = x^3+y+y^2, B = y^3+x+x^2
    "[[72,12,6]]": BBSpec(6, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)],
                          name="[[72,12,6]]"),
    # [[90,8,10]]  ell=15, m=3, A = x^9+y+y^2, B = 1+x^2+x^7
    "[[90,8,10]]": BBSpec(15, 3, [(9, 0), (0, 1), (0, 2)], [(0, 0), (2, 0), (7, 0)],
                          name="[[90,8,10]]"),
    # [[108,8,10]] ell=9, m=6,  A = x^3+y+y^2, B = y^3+x+x^2
    "[[108,8,10]]": BBSpec(9, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)],
                           name="[[108,8,10]]"),
    # [[144,12,12]] ell=12, m=6, A = x^3+y+y^2, B = y^3+x+x^2   (the "Gross code")
    "[[144,12,12]]": BBSpec(12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)],
                            name="[[144,12,12]]"),
    # [[288,12,18]] ell=12, m=12, A = x^3+y^2+y^7, B = y^3+x+x^2
    "[[288,12,18]]": BBSpec(12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)],
                            name="[[288,12,18]]"),
    # [[360,12,<=24]] ell=30, m=6, A = x^9+y+y^2, B = y^3+x^25+x^26
    "[[360,12,<=24]]": BBSpec(30, 6, [(9, 0), (0, 1), (0, 2)], [(0, 3), (25, 0), (26, 0)],
                              name="[[360,12,<=24]]"),
    # [[784,24,24]] ell=28, m=14, A = x^26+y^6+y^8, B = y^7+x^9+x^20
    "[[784,24,24]]": BBSpec(28, 14, [(26, 0), (0, 6), (0, 8)], [(0, 7), (9, 0), (20, 0)],
                            name="[[784,24,24]]"),
}
