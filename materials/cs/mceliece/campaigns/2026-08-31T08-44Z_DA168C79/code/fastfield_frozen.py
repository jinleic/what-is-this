"""Vectorized exact linear algebra over E = F_{2^m} for the Gate A census.

Representation: field elements are bitmasks in [0, 2^m) (uint16). A full
multiplication table MUL[q, q] (uint16) makes all vectorized ops gathers.
Everything is exact: table construction is verified by MUL[x,y] == gf.mul(x,y)
against the scalar engine on random samples at import of the census runner.
"""
from __future__ import annotations

import numpy as np

from gfield_frozen import GF   # frozen byte copy (campaign provenance lock)


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
        """[a^0, ..., a^D] uint16 — vectorized via log-domain (exact)."""
        out = np.zeros(D + 1, dtype=np.uint16)
        out[0] = 1
        if a != 0:
            la = int(self.LOG[a])
            d = np.arange(1, D + 1, dtype=np.int64)
            out[1:] = self.EXP[(la * d) % (self.q - 1)].astype(np.uint16)
        # a == 0: remaining entries already 0
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


def vec_pis_irreducible(ef, G_full):
    """Vectorized Rabin irreducibility over E = F_{2^m} (deg t poly).

    Coefficients live in E = F_q and are FIXED by the q-Frobenius, so
    x -> x^q is pure exponent stretching. We precompute w = Z^q mod G
    (ONE stretch-reduce) and the basis B_j = w^j mod G; then frob(x) =
    sum_j x_j B_j is O(t^2). Rabin: Z^{q^t} == Z (mod G) and
    gcd(Z^{q^{t/f}} - Z, G) constant for prime f | t.

    Maintained in mceliece/src only (ownership rule; NO cross-folder imports).
    Validated against the scalar gfield.GF.pis_irreducible (61/61 over
    F_{2^12}, degree 1-3 random) which agrees with flint fq_default
    is_irreducible on the same inputs.
    """
    import numpy as np
    from gfield_frozen import GF, _factor   # frozen byte copy
    G = np.asarray(G_full, dtype=np.uint16)
    while len(G) and G[-1] == 0:
        G = G[:-1]
    t = len(G) - 1
    if t < 1:
        return False
    if t == 1:
        return True
    ig = ef.INV[G[t]]

    def polymod(f):
        f = f.copy()
        dg = t + 1
        for i in range(len(f) - 1, dg - 2, -1):
            c = f[i]
            if c == 0:
                continue
            fct = ef.MUL[c, ig]
            f[i - dg + 1:i + 1] ^= ef.MUL[fct, G[:dg]]
        return f[:dg - 1]

    def mulmod2(a, b):
        out = np.zeros(2 * t + 1, dtype=np.uint16)
        for i in np.nonzero(a)[0]:
            seg = ef.MUL[a[i], b]
            L = min(len(seg), len(out) - i)
            out[i:i + L] ^= seg[:L]
        return polymod(out)

    # w = Z^q mod G: single exponent stretch of Z, then reduce
    stretch = np.zeros(ef.q + 1, dtype=np.uint16)
    stretch[ef.q] = 1
    w = polymod(stretch)
    B = [np.zeros(t, dtype=np.uint16) for _ in range(t)]
    B[0][0] = 1
    for j in range(1, t):
        B[j] = mulmod2(B[j - 1], w)

    def frob(xr):
        out = np.zeros(t, dtype=np.uint16)
        for j in np.nonzero(xr)[0]:
            out ^= ef.MUL[xr[j], B[j]]
        return out

    xr = np.zeros(t, dtype=np.uint16)
    if t >= 2:
        xr[1] = 1
    else:
        xr[0] = 1
    for _ in range(t):
        xr = frob(xr)
    if not (xr[1] == 1 and xr[0] == 0 and not np.any(xr[2:])):
        return False
    gfc = GF(ef.m, prim=ef.prim)
    Gl = [int(v) for v in G]
    for f in _factor(t):
        e = t // f
        yr = xr.copy()
        for _ in range(e):
            yr = frob(yr)
        diff = [int(v) for v in yr] + [0] * max(0, 2 - len(yr))
        diff[1] ^= 1
        g = gfc.pgcd(gfc.ptrim(diff), Gl)
        if gfc.pdeg(g) != 0:
            return False
    return True
