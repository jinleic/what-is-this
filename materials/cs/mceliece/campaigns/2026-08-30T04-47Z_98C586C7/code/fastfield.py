"""Vectorized exact linear algebra over E = F_{2^m} for the Gate A census.

Representation: field elements are bitmasks in [0, 2^m) (uint16). A full
multiplication table MUL[q, q] (uint16) makes all vectorized ops gathers.
Everything is exact: table construction is verified by MUL[x,y] == gf.mul(x,y)
against the scalar engine on random samples at import of the census runner.
"""
from __future__ import annotations

import numpy as np

from gfield import GF


class EField:
    """numpy engine over F_{2^m}; tables of size q^2 (<= 16.7M entries, m<=12)."""

    def __init__(self, m: int, prim: int | None = None):
        gf = GF(m) if prim is None else GF(m, prim=prim)
        self.gf = gf
        self.m, self.q, self.prim = m, gf.q, gf.prim
        q = self.q
        self.LOG = np.full(q, q - 1, dtype=np.int64)
        self.EXP = np.zeros(2 * q, dtype=np.int64)
        x = 1
        for i in range(q - 1):
            self.EXP[i] = x
            self.LOG[x] = i
            x <<= 1
            if x & q:
                x ^= self.prim
        assert x == 1
        for i in range(q - 1, 2 * q):
            self.EXP[i] = self.EXP[i - (q - 1)]
        # MUL table (q x q uint16); row/col 0 are zero
        self.MUL = np.zeros((q, q), dtype=np.uint16)
        nz = np.arange(1, q, dtype=np.int64)
        lg = self.LOG[1:q]
        e = self.EXP
        for a in range(1, q):
            self.MUL[a, 1:q] = e[lg + self.LOG[a]].astype(np.uint16)
        self.INV = np.zeros(q, dtype=np.uint16)
        self.INV[1:q] = e[(q - 1) - lg].astype(np.uint16)
        self.EXP16 = self.EXP.astype(np.uint16)

    # ---------- scalar ----------
    def smul(self, a, b):
        return int(self.MUL[a, b])

    def sinv(self, a):
        assert a
        return int(self.INV[a])

    # ---------- vectors / matrices ----------
    def mulvv(self, a, b):
        return self.MUL[a, b] if np.ndim(a) == np.ndim(b) else None

    def mulsv(self, s, v):
        return self.MUL[int(s), v]

    def rref(self, M):
        """Reduced row echelon over E. M: uint16 2D. Returns (M', pivots)."""
        M = M.copy()
        rows, cols = M.shape
        piv = []
        r = 0
        for c in range(cols):
            if r >= rows:
                break
            colv = M[r:, c]
            nz = np.nonzero(colv)[0]
            if nz.size == 0:
                continue
            pr = r + int(nz[0])
            if pr != r:
                M[[r, pr]] = M[[pr, r]]
            iv = self.INV[M[r, c]]
            if iv != 1:
                M[r] = self.MUL[iv, M[r]]
            f = M[:, c].copy()
            f[r] = 0
            R = np.nonzero(f)[0]
            if R.size:
                M[R] ^= self.MUL[f[R][:, None], M[r][None, :]]
            piv.append(c)
            r += 1
        return M, piv

    def rank(self, M):
        return len(self.rref(M)[1])

    def row_nullspace(self, M):
        """basis (rows) of {x : M x = 0}."""
        M, piv = self.rref(M)
        rows, cols = M.shape
        pivset = set(piv)
        free = [c for c in range(cols) if c not in pivset]
        out = []
        for fc in free:
            v = np.zeros(cols, dtype=np.uint16)
            v[fc] = 1
            for ri, pc in enumerate(piv):
                if M[ri, fc]:
                    v[pc] = M[ri, fc]
            out.append(v)
        return out

    def left_complement(self, W):
        """rows spanning {u : u·w = 0 for all rows w of W} (standard dot)."""
        W = np.asarray(W, dtype=np.uint16)
        k = W.shape[1]
        if W.shape[0] == 0:
            return np.eye(k, dtype=np.uint16)
        aug = np.zeros((W.shape[0], W.shape[0] + k), dtype=np.uint16)
        aug[:, : W.shape[0]] = np.eye(W.shape[0], dtype=np.uint16)
        # find u with W u^T = 0 ... equivalently solve (W) x = 0 for x as
        # column: u·w = sum u_i w_i = w·u.  So complement = right nullspace of W.
        K = self.row_nullspace(W)  # rows are vectors in E^k with W v^T = 0
        return np.array(K, dtype=np.uint16) if K else np.zeros((0, k), dtype=np.uint16)

    def powers(self, a, D):
        """[a^0, ..., a^D] uint16"""
        out = np.zeros(D + 1, dtype=np.uint16)
        out[0] = 1
        if a != 0:
            aa = np.full(D + 1, a, dtype=np.uint16)
            out[1:] = np.cumprod(aa[1:]) if False else out[1:]
            # iterative exact pow via MUL
            for d in range(1, D + 1):
                out[d] = self.MUL[out[d - 1], a]
        else:
            out[1:] = 0
        return out

    def jet_row(self, D, j, a, apow):
        """coeffs w_d = binom(d, j) * a^{d-j} for d in 0..D (binom mod 2 via Lucas)."""
        d = np.arange(D + 1)
        mask = ((d & j) == j) & (d >= j)
        w = np.zeros(D + 1, dtype=np.uint16)
        if j == 0:
            w[mask] = self.MUL[np.ones(int(mask.sum()), dtype=np.uint16), apow[: D + 1][mask]] if a != 0 else np.where(mask, np.uint16(1), np.uint16(0))
            return w
        w[mask] = apow[d[mask] - j]
        return w

    def jet_eval_poly(self, coeffs, j, a, apow):
        """Hasse jet value (∂[j]f)(a) = sum_d c_d binom(d,j) a^{d-j}."""
        coeffs = np.asarray(coeffs, dtype=np.uint16)
        w = self.jet_row(len(coeffs) - 1, j, a, apow)
        nz = np.nonzero(coeffs)[0]
        if nz.size == 0 or not w[nz].any():
            # still must multiply
            pass
        terms = self.MUL[coeffs, w]
        return np.bitwise_xor.reduce(terms) if terms.size else 0

    def series_inv(self, s, order):
        """invert a power series (uint16 coeffs, s[0] != 0) to given order.
        Recurrence: (s * out)_e = 0 for e >= 1, out_0 = 1/s_0:
            out[e] = s_0^{-1} * sum_{i=1..e} s_i out_{e-i}
        """
        s = np.asarray(s, dtype=np.uint16)
        out = np.zeros(order + 1, dtype=np.uint16)
        out[0] = self.INV[s[0]]
        for e in range(1, order + 1):
            i_hi = min(e, len(s) - 1)
            i_idx = np.arange(1, i_hi + 1)
            if i_idx.size == 0:
                out[e] = 0
                continue
            terms = self.MUL[s[i_idx], out[e - i_idx]]
            acc = np.bitwise_xor.reduce(terms) if terms.size else 0
            out[e] = self.MUL[acc, self.INV[s[0]]]
        return out

    def series_mul(self, u, v, order):
        w = np.zeros(order + 1, dtype=np.uint16)
        u = u[: order + 1]
        v = v[: order + 1]
        for e in range(order + 1):
            lo = max(0, e - len(v) + 1)
            hi = min(e, len(u) - 1)
            if hi < lo:
                continue
            uu = u[lo : hi + 1]
            vv = v[e - hi : e - lo + 1]
            terms = self.MUL[uu, vv[::-1]]
            if terms.size:
                w[e] = np.bitwise_xor.reduce(terms)
        return w
