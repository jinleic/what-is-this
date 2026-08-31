"""GF(2) linear algebra, bit-packed. Scope: Gate-D pilot only.

Conventions: a binary matrix is stored dense (rows x cols) as uint8 0/1.
The packed form packs rows into uint8 bytes along axis 1: bit j of byte j//8
of row r is entry (r, j). Under this layout, XOR-eliminating a pivot row into
another row is a single vectorized `P[r2] ^= P[r1]` -- the low-level trick
every routine here relies on.
"""
from __future__ import annotations

import numpy as np


def _mul2_bits(M: np.ndarray, x: np.ndarray) -> np.ndarray:
    """A @ x over GF(2), exact for any sizes (dot in wider int then mod 2)."""
    return (M.astype(np.uint64) @ x.astype(np.uint64)) % 2


def gf2_rank(M: np.ndarray) -> int:
    """GF(2) rank of a dense 0/1 matrix (rows x cols)."""
    P = np.packbits(M.astype(np.uint8), axis=1, bitorder="little")
    m, n = M.shape
    r = 0
    for c in range(n):
        byte = c >> 3
        mask = np.uint8(1 << (c & 7))
        cand = np.nonzero(P[:, byte] & mask)[0]
        if cand.size == 0:
            continue
        cand = cand[cand >= r]
        if cand.size == 0:
            continue
        pr = cand[0]
        if pr != r:
            P[[r, pr]] = P[[pr, r]]
        rows = np.nonzero(P[:, byte] & mask)[0]
        rows = rows[rows != r]
        if rows.size:
            P[rows] ^= P[r]
        r += 1
        if r == m:
            break
    return r


def solve_gf2(A_mat: np.ndarray, b: np.ndarray) -> np.ndarray | None:
    """Any particular solution of A x = b over GF(2); None if inconsistent.

    A_mat dense (m, n) uint8/bool; b (m,) 0/1. Bit-packed row-major
    elimination. The b column is padded to its own byte so it never
    aliases A bits during packed elimination.
    """
    A_mat = np.ascontiguousarray(A_mat, dtype=np.uint8)
    m, n = A_mat.shape
    npad = 8 * ((n + 7) // 8)
    aug = np.zeros((m, npad + 1), dtype=np.uint8)
    aug[:, :n] = A_mat
    aug[:, npad] = np.asarray(b, dtype=np.uint8).ravel()
    P = np.packbits(aug, axis=1, bitorder="little")
    nb = P.shape[1]
    pivots: list[int] = []
    r = 0
    for c in range(n):
        byte = c >> 3
        mask = np.uint8(1 << (c & 7))
        cand = np.nonzero(P[:, byte] & mask)[0]
        cand = cand[cand >= r]
        if cand.size == 0:
            continue
        pr = cand[0]
        if pr != r:
            P[[r, pr]] = P[[pr, r]]
        rows = np.nonzero(P[:, byte] & mask)[0]
        rows = rows[rows != r]
        if rows.size:
            P[rows] ^= P[r]
        pivots.append(c)
        r += 1
        if r == m:
            break
    # consistency: rows r..m-1 have zero A part; their b bit must be 0
    if r < m and np.any(P[r:, nb - 1]):
        return None
    # back-substitution over pivots (free vars stay 0)
    x = np.zeros(n, dtype=np.uint8)
    Aech = np.unpackbits(P[:, : nb - 1], axis=1, bitorder="little")[:, :n]
    bech = np.unpackbits(P[:, nb - 1][:, None].copy(), axis=1, bitorder="little")[:, 0]
    for i, c in enumerate(reversed(pivots)):
        row = r - 1 - i
        s = int(bech[row])
        sel = np.nonzero(Aech[row, c + 1:])[0]
        for j in sel:
            s ^= int(x[c + 1 + j])
        x[c] = s & 1
    if not np.array_equal((A_mat @ x.astype(np.uint64)) % 2,
                          (np.asarray(b, dtype=np.uint64) % 2)):
        return None
    return x


def kernel_basis(M: np.ndarray) -> np.ndarray:
    """Basis (rows) of the right-kernel of M: {x : M x = 0}, dense uint8 0/1.

    Column reduction with row echelon: bring M to [I_r B; 0 C] form; kernel
    basis rows are [B^T-col rows; C] read off the free columns.
    """
    M = M.astype(np.uint8)
    m, n = M.shape
    A = M.copy()
    r = 0
    pivot_cols: list[int] = []
    for c in range(n):
        sel = np.nonzero(A[r:, c])[0]
        if sel.size == 0:
            continue
        pr = r + int(sel[0])
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
        rows = np.nonzero(A[:, c])[0]
        rows = rows[rows != r]
        if rows.size:
            A[rows] ^= A[r]
        pivot_cols.append(c)
        r += 1
        if r == m:
            break
    piv = set(pivot_cols)
    free_cols = [c for c in range(n) if c not in piv]
    K = np.zeros((len(free_cols), n), dtype=np.uint8)
    # forward-substitute: for each free col f, kernel vector has 1 at f and
    # 1 at every pivot col p whose pivot row has a 1 in col f (looking at the
    # reduced A: rows are ordered by pivot discovery, pivot i at pivot_cols[i])
    for fi, f in enumerate(free_cols):
        K[fi, f] = 1
        for i, p in enumerate(pivot_cols):
            if A[i, f]:
                K[fi, p] = 1
    if not np.any(K):
        return np.zeros((0, n), dtype=np.uint8)
    assert np.array_equal((M @ K.T % 2).astype(np.uint64), np.zeros((m, K.shape[0]), dtype=np.uint64))
    return K
