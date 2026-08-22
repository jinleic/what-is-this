"""Exact GF(2) linear algebra.

Two independent representations are maintained deliberately so that every
algebraic claim in this project can be cross-checked by two code paths that do
NOT share an implementation:

  * ``bitset``  -- rows are Python ``int`` bitmasks (arbitrary precision, exact).
  * ``numpy``   -- rows are ``uint8`` arrays mod 2.

Nothing here depends on floating point.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "rows_to_bitsets",
    "bitsets_to_rows",
    "rref_bitset",
    "rank_bitset",
    "rank_np",
    "rref_np",
    "nullspace_bitset",
    "nullspace_np",
    "solve_bitset",
    "row_space_contains",
    "matmul",
    "popcount",
]


# --------------------------------------------------------------------------
# conversions
# --------------------------------------------------------------------------
def rows_to_bitsets(M: np.ndarray) -> list[int]:
    """Pack each row of a 0/1 matrix into a Python int (bit j == column j)."""
    M = np.asarray(M, dtype=np.uint8) & 1
    if M.ndim != 2:
        raise ValueError(f"expected 2-D matrix, got shape {M.shape}")
    out: list[int] = []
    for row in M:
        # little-endian: column j -> bit j
        out.append(int.from_bytes(np.packbits(row, bitorder="little").tobytes(), "little"))
    return out


def bitsets_to_rows(bs: list[int], ncols: int) -> np.ndarray:
    M = np.zeros((len(bs), ncols), dtype=np.uint8)
    for i, b in enumerate(bs):
        nbytes = (ncols + 7) // 8
        arr = np.frombuffer(int(b).to_bytes(nbytes, "little"), dtype=np.uint8)
        M[i] = np.unpackbits(arr, bitorder="little")[:ncols]
    return M


def popcount(x: int) -> int:
    return int(x).bit_count()


# --------------------------------------------------------------------------
# bitset path
# --------------------------------------------------------------------------
def rref_bitset(rows: list[int], ncols: int) -> tuple[list[int], list[int], list[list[int]]]:
    """Reduced row echelon form over GF(2).

    Returns ``(rref_rows, pivot_cols, transform)`` where ``transform[i]`` lists
    the indices of the *original* rows whose sum equals ``rref_rows[i]``.
    """
    rows = list(rows)
    prov: list[set[int]] = [{i} for i in range(len(rows))]
    pivots: list[int] = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(rows)):
            if (rows[i] >> c) & 1:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        prov[r], prov[piv] = prov[piv], prov[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> c) & 1):
                rows[i] ^= rows[r]
                prov[i] ^= prov[r]
        pivots.append(c)
        r += 1
        if r == len(rows):
            break
    # drop zero rows (they carry provenance of a dependency, kept separately)
    keep = [i for i in range(len(rows)) if rows[i] != 0]
    return ([rows[i] for i in keep], pivots, [sorted(prov[i]) for i in keep])


def rank_bitset(rows: list[int], ncols: int) -> int:
    rr, _, _ = rref_bitset(rows, ncols)
    return len(rr)


def nullspace_bitset(rows: list[int], ncols: int) -> list[int]:
    """Basis of {v : M v^T = 0} as bitsets of length ``ncols``."""
    rr, pivots, _ = rref_bitset(rows, ncols)
    pivset = set(pivots)
    free = [c for c in range(ncols) if c not in pivset]
    basis: list[int] = []
    for f in free:
        v = 1 << f
        for i, pc in enumerate(pivots):
            if (rr[i] >> f) & 1:
                v |= 1 << pc
        basis.append(v)
    return basis


def solve_bitset(rows: list[int], ncols: int, target: int) -> int | None:
    """Find x with sum_{i in x} rows[i] == target; return bitmask over rows."""
    m = len(rows)
    aug = [rows[i] | (1 << (ncols + i)) for i in range(m)]
    rr, pivots, _ = rref_bitset(aug, ncols + m)
    mask = (1 << ncols) - 1
    cur = target
    combo = 0
    for i, pc in enumerate(pivots):
        if pc >= ncols:
            break
        if (cur >> pc) & 1:
            cur ^= rr[i] & mask
            combo ^= rr[i] >> ncols
    return combo if cur == 0 else None


def row_space_contains(rows: list[int], ncols: int, target: int) -> bool:
    return solve_bitset(rows, ncols, target) is not None


# --------------------------------------------------------------------------
# numpy path (independent implementation)
# --------------------------------------------------------------------------
def rref_np(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    A = (np.asarray(M, dtype=np.uint8) & 1).copy()
    nrows, ncols = A.shape
    pivots: list[int] = []
    r = 0
    for c in range(ncols):
        col = A[r:, c]
        nz = np.flatnonzero(col)
        if nz.size == 0:
            continue
        p = r + int(nz[0])
        if p != r:
            A[[r, p]] = A[[p, r]]
        sel = np.flatnonzero(A[:, c])
        sel = sel[sel != r]
        if sel.size:
            A[sel] ^= A[r]
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    return A[:r], pivots


def rank_np(M: np.ndarray) -> int:
    if np.asarray(M).size == 0:
        return 0
    return rref_np(M)[0].shape[0]


def nullspace_np(M: np.ndarray) -> np.ndarray:
    A = np.asarray(M, dtype=np.uint8) & 1
    ncols = A.shape[1]
    R, pivots = rref_np(A)
    pivset = set(pivots)
    free = [c for c in range(ncols) if c not in pivset]
    if not free:
        return np.zeros((0, ncols), dtype=np.uint8)
    N = np.zeros((len(free), ncols), dtype=np.uint8)
    for j, f in enumerate(free):
        N[j, f] = 1
        for i, pc in enumerate(pivots):
            N[j, pc] = R[i, f]
    return N


def matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """GF(2) matrix product, overflow-safe for large inner dimensions."""
    A = np.asarray(A, dtype=np.uint8) & 1
    B = np.asarray(B, dtype=np.uint8) & 1
    if A.shape[1] != B.shape[0]:
        raise ValueError(f"shape mismatch {A.shape} @ {B.shape}")
    # int64 accumulate then reduce; inner dim <= ~9e18 is fine
    return (A.astype(np.int64) @ B.astype(np.int64) % 2).astype(np.uint8)
