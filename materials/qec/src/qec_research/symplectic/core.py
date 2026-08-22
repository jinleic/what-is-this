"""Binary symplectic representation of stabilizer codes.

Convention (fixed for the whole project):
    an n-qubit Pauli is  v = (x | z) in GF(2)^{2n}
    (x_j,z_j) = (0,0)->I  (1,0)->X  (0,1)->Z  (1,1)->Y   (up to phase)

    symplectic product  <(x|z),(x'|z')> = x.z' + z.x'  (mod 2)

A stabilizer check matrix is  H = (H_X | H_Z)  with r rows.
Validity:  H_X H_Z^T + H_Z H_X^T = 0 (mod 2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..gf2.linalg import (
    matmul,
    nullspace_bitset,
    nullspace_np,
    rank_bitset,
    rank_np,
    rows_to_bitsets,
    bitsets_to_rows,
    row_space_contains,
)

__all__ = [
    "symplectic_product_matrix",
    "is_valid_stabilizer",
    "symplectic_weight",
    "pauli_type_counts",
    "StabilizerCode",
    "lambda_swap",
]


def lambda_swap(V: np.ndarray) -> np.ndarray:
    """Apply the symplectic form matrix Lambda: (x|z) -> (z|x)."""
    V = np.asarray(V, dtype=np.uint8)
    n = V.shape[-1] // 2
    return np.concatenate([V[..., n:], V[..., :n]], axis=-1)


def symplectic_product_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Return S with S[i,j] = <A_i, B_j>_s."""
    A = np.asarray(A, dtype=np.uint8) & 1
    B = np.asarray(B, dtype=np.uint8) & 1
    n = A.shape[1] // 2
    Ax, Az = A[:, :n], A[:, n:]
    Bx, Bz = B[:, :n], B[:, n:]
    return (matmul(Ax, Bz.T) ^ matmul(Az, Bx.T)).astype(np.uint8)


def is_valid_stabilizer(H: np.ndarray) -> bool:
    S = symplectic_product_matrix(H, H)
    return not S.any()


def symplectic_weight(v: np.ndarray) -> int:
    """Number of qubits acted on non-trivially (I is free, X/Y/Z each cost 1)."""
    v = np.asarray(v, dtype=np.uint8) & 1
    n = v.shape[-1] // 2
    return int(np.count_nonzero(v[..., :n] | v[..., n:]))


def pauli_type_counts(v: np.ndarray) -> dict[str, int]:
    v = np.asarray(v, dtype=np.uint8) & 1
    n = len(v) // 2
    x, z = v[:n].astype(bool), v[n:].astype(bool)
    return {
        "I": int(np.count_nonzero(~x & ~z)),
        "X": int(np.count_nonzero(x & ~z)),
        "Z": int(np.count_nonzero(~x & z)),
        "Y": int(np.count_nonzero(x & z)),
    }


def pauli_string(v: np.ndarray) -> str:
    v = np.asarray(v, dtype=np.uint8) & 1
    n = len(v) // 2
    tbl = {(0, 0): "_", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}
    return "".join(tbl[(int(v[j]), int(v[n + j]))] for j in range(n))


@dataclass
class StabilizerCode:
    """A stabilizer code given by an arbitrary (possibly redundant) check matrix."""

    H: np.ndarray  # (r, 2n)
    name: str = ""

    def __post_init__(self) -> None:
        self.H = np.asarray(self.H, dtype=np.uint8) & 1
        if self.H.ndim != 2 or self.H.shape[1] % 2:
            raise ValueError(f"bad check matrix shape {self.H.shape}")

    # ---- basic invariants -------------------------------------------------
    @property
    def n(self) -> int:
        return self.H.shape[1] // 2

    @property
    def num_checks(self) -> int:
        return self.H.shape[0]

    def rank(self, path: str = "bitset") -> int:
        if path == "bitset":
            return rank_bitset(rows_to_bitsets(self.H), 2 * self.n)
        return rank_np(self.H)

    @property
    def k(self) -> int:
        return self.n - self.rank()

    def validate(self) -> dict:
        """Two independent code paths for every algebraic invariant."""
        n = self.n
        # path 1: numpy
        S_np = symplectic_product_matrix(self.H, self.H)
        comm_np = not S_np.any()
        rank_a = rank_np(self.H)
        # path 2: pure-python bitsets, independent arithmetic
        bs = rows_to_bitsets(self.H)
        comm_bs = True
        for i, a in enumerate(bs):
            ax, az = a & ((1 << n) - 1), a >> n
            for b in bs[i:]:
                bx, bz = b & ((1 << n) - 1), b >> n
                if ((ax & bz).bit_count() + (az & bx).bit_count()) & 1:
                    comm_bs = False
                    break
            if not comm_bs:
                break
        rank_b = rank_bitset(bs, 2 * n)
        return {
            "n": n,
            "num_checks": self.num_checks,
            "commutes_numpy": comm_np,
            "commutes_bitset": comm_bs,
            "rank_numpy": rank_a,
            "rank_bitset": rank_b,
            "rank_agree": rank_a == rank_b,
            "k": n - rank_a,
            "check_weights": self.check_weights().tolist(),
            "max_check_weight": int(self.check_weights().max()) if self.num_checks else 0,
            "qubit_degrees_max": int(self.qubit_degrees().max()) if n else 0,
            "is_css": self.is_css(),
        }

    def check_weights(self) -> np.ndarray:
        n = self.n
        return np.count_nonzero(self.H[:, :n] | self.H[:, n:], axis=1)

    def qubit_degrees(self) -> np.ndarray:
        n = self.n
        return np.count_nonzero(self.H[:, :n] | self.H[:, n:], axis=0)

    def is_css(self) -> bool:
        """True iff every generator (as given) is pure-X or pure-Z."""
        n = self.n
        px = self.H[:, :n].any(axis=1)
        pz = self.H[:, n:].any(axis=1)
        return not bool((px & pz).any())

    # ---- centralizer / logicals ------------------------------------------
    def centralizer_basis(self) -> np.ndarray:
        """Basis of S^{perp_s}: all v with <v, s>=0 for every stabilizer s.

        <v, h> = v_x . h_z + v_z . h_x = v . (Lambda h).  So the constraint
        matrix is Lambda applied to H, i.e. columns swapped.
        """
        M = lambda_swap(self.H)  # (r, 2n); constraint  M v^T = 0
        return nullspace_np(M)

    def logical_basis(self) -> np.ndarray:
        """A basis of 2k representatives of S^perp / S (rows of length 2n)."""
        n = self.n
        C = self.centralizer_basis()
        Hb = rows_to_bitsets(self.H)
        Cb = rows_to_bitsets(C)
        # greedily pick centralizer elements independent modulo the stabilizer group
        chosen: list[int] = []
        cur = list(Hb)
        base_rank = rank_bitset(cur, 2 * n)
        for c in Cb:
            trial = cur + [c]
            r = rank_bitset(trial, 2 * n)
            if r > base_rank:
                chosen.append(c)
                cur = trial
                base_rank = r
        return bitsets_to_rows(chosen, 2 * n)

    def is_logical(self, v: np.ndarray) -> bool:
        """v is a nontrivial logical: in centralizer, not in stabilizer group."""
        v = np.asarray(v, dtype=np.uint8) & 1
        if symplectic_product_matrix(v[None, :], self.H).any():
            return False
        return not row_space_contains(rows_to_bitsets(self.H), 2 * self.n,
                                      rows_to_bitsets(v[None, :])[0])

    # ---- structure --------------------------------------------------------
    def connected_components(self) -> list[list[int]]:
        """Partition qubits by the check-sharing relation (detects direct sums)."""
        n = self.n
        supp = self.H[:, :n] | self.H[:, n:]
        parent = list(range(n))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for row in supp:
            js = np.flatnonzero(row)
            for j in js[1:]:
                union(int(js[0]), int(j))
        groups: dict[int, list[int]] = {}
        for j in range(n):
            groups.setdefault(find(j), []).append(j)
        return sorted(groups.values(), key=len, reverse=True)

    def is_direct_sum(self) -> bool:
        return len(self.connected_components()) > 1
