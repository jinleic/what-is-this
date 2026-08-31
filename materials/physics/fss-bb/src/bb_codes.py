"""Bivariate bicycle (BB) code parity-check matrices.

Family studied (arXiv:2603.19062v3 Sec 2.1, following Bravyi et al.
arXiv:2308.07915): A = x^3 + y + y^2, B = y^3 + x + x^2, with

    Hx = [A | B],   Hz = [B^T | A^T],   N = 2*L*M.

A and B are L*M x L*M binary circulant matrices over Z_L x Z_M: the monomial
x^a y^b contributes kron(shift_L(a), shift_M(b)), shift_k(s)[i, (i+s) % k] = 1.

All matrices are GF(2). Only matrix structure matters for the erasure channel;
any qubit/check permutation of the canonical code has identical WER.
"""

from __future__ import annotations

import numpy as np

# Monomial shift lists (dx, dy) for the study family.
A_SHIFTS = ((3, 0), (0, 1), (0, 2))  # x^3 + y + y^2
B_SHIFTS = ((0, 3), (1, 0), (2, 0))  # y^3 + x + x^2

# (L, M) -> (N, K) for the five study sizes; K from arXiv:2603.19062v3 Table 1
# (K = 12, 8, 16, 8, 12). Used as a construction self-check.
STUDY_SIZES = {
    (12, 6): (144, 12),
    (18, 9): (324, 8),
    (24, 12): (576, 16),
    (30, 15): (900, 8),
    (36, 18): (1296, 12),
}


def _shift_matrix(size: int, shift: int) -> np.ndarray:
    """size x size binary cyclic shift matrix, shift[s][ (i+shift) % size ] = 1."""
    S = np.zeros((size, size), dtype=np.uint8)
    rows = np.arange(size)
    S[rows, (rows + shift) % size] = 1
    return S


def _circulant_block(shape: tuple[int, int], shifts: tuple) -> np.ndarray:
    """LM x LM circulant over Z_L x Z_M from monomial shifts (dx, dy)."""
    L, M = shape
    block = np.zeros((L * M, L * M), dtype=np.uint8)
    for dx, dy in shifts:
        block += np.kron(_shift_matrix(L, dx), _shift_matrix(M, dy))
    return block % 2


def bb_parity_checks(L: int, M: int, a_shifts=A_SHIFTS, b_shifts=B_SHIFTS):
    """Return (Hx, Hz) as uint8 arrays, shapes (L*M, 2*L*M) each."""
    A = _circulant_block((L, M), a_shifts)
    B = _circulant_block((L, M), b_shifts)
    Hx = np.hstack([A, B]).astype(np.uint8)
    Hz = np.hstack([B.T, A.T]).astype(np.uint8)
    return Hx, Hz


def gf2_rref(mat: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Reduced row echelon form over GF(2). Returns (rref, pivot_columns)."""
    A = mat.copy().astype(np.uint8) % 2
    n_rows, n_cols = A.shape
    pivots = []
    row = 0
    for col in range(n_cols):
        if row >= n_rows:
            break
        pivot_rows = np.nonzero(A[row:, col])[0]
        if pivot_rows.size == 0:
            continue
        pr = row + pivot_rows[0]
        if pr != row:
            A[[row, pr]] = A[[pr, row]]
        # eliminate all other rows
        rows_with_one = np.nonzero(A[:, col])[0]
        rows_with_one = rows_with_one[rows_with_one != row]
        A[rows_with_one] ^= A[row]
        pivots.append(col)
        row += 1
    return A, pivots


def gf2_nullspace_basis(mat: np.ndarray) -> np.ndarray:
    """Basis (rows) of the right kernel of `mat` over GF(2)."""
    A, pivots = gf2_rref(mat)
    n_rows, n_cols = A.shape
    free = [c for c in range(n_cols) if c not in set(pivots)]
    basis = []
    for f in free:
        v = np.zeros(n_cols, dtype=np.uint8)
        v[f] = 1
        for r, pc in enumerate(pivots):
            if A[r, f]:
                v[pc] = 1
        basis.append(v)
    return np.array(basis, dtype=np.uint8) if basis else np.zeros((0, n_cols), dtype=np.uint8)


def gf2_rank(mat: np.ndarray) -> int:
    _, pivots = gf2_rref(mat)
    return len(pivots)


def logical_operator_basis(Hx: np.ndarray, Hz: np.ndarray):
    """K x N bases of Z-type and X-type logical operators (rows).

    Z-logicals span ker(Hx) / rowspace(Hz) (CSS: rowspace(Hz) ⊆ ker(Hx) since
    Hx·Hz^T = 0). Construction: basis of ker(Hx), reduce each vector against
    rref(Hz) (pivot-column elimination = canonical coset representatives),
    then rref the survivors and keep nonzero rows — an independent basis.
    X-logicals symmetrically: ker(Hz) / rowspace(Hx).
    """
    def _quotient_basis(ker_mat, mod_mat):
        ker = gf2_nullspace_basis(ker_mat)
        R, pivots = gf2_rref(mod_mat)
        survivors = []
        for v in ker:
            w = v.copy()
            for r, pc in enumerate(pivots):
                if w[pc]:
                    w ^= R[r]
            if w.any():
                survivors.append(w)
        if not survivors:
            return np.zeros((0, ker_mat.shape[1]), dtype=np.uint8)
        S, spiv = gf2_rref(np.array(survivors, dtype=np.uint8))
        keep = [r for r in range(S.shape[0]) if S[r].any()]
        return S[keep]

    Lz = _quotient_basis(Hx, Hz)
    Lx = _quotient_basis(Hz, Hx)
    return Lz, Lx



def bb_code(L: int, M: int):
    """Build checks + logical bases for one BB code. Returns dict."""
    Hx, Hz = bb_parity_checks(L, M)
    N = Hx.shape[1]
    K_expected = STUDY_SIZES.get((L, M))
    K = N - gf2_rank(Hx) - gf2_rank(Hz)
    if K_expected is not None:
        assert K == K_expected[1], (
            f"K mismatch for (L,M)=({L},{M}): constructed K={K}, expected {K_expected[1]} "
            "(arXiv:2603.19062v3 Table 1)"
        )
    Lz, Lx = logical_operator_basis(Hx, Hz)
    assert Lz.shape[0] == K, f"Z-logical count {Lz.shape[0]} != K {K}"
    assert Lx.shape[0] == K, f"X-logical count {Lx.shape[0]} != K {K}"
    return {
        "L": L,
        "M": M,
        "N": N,
        "K": K,
        "Hx": Hx,
        "Hz": Hz,
        "Lz": Lz,  # rows: Z-logicals; used to pass/fail Z-sector (Hx-syndrome) residuals
        "Lx": Lx,  # rows: X-logicals; pass/fail X-sector residuals
    }
