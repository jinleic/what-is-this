"""finite field F_{2^m} utilities, byte-exact, no dependencies.

Elements are ints in [0, 2^m) — bitmask of the polynomial-basis
coefficients (basis 1, alpha, alpha^2, ..., alpha^{m-1} where alpha is a
root of the primitive modulus). Multiplication via log/exp tables.

Polynomials are lists of field elements, index = degree.
"""
from __future__ import annotations

import random


class GF:
    def __init__(self, m: int, prim: int | None = None, seed: int = 0x5EED):
        self.m = m
        self.q = 1 << m
        if prim is None:
            prim = self._find_primitive(seed)
        self.prim = prim
        # log/exp tables
        self.exp = [0] * (2 * self.q)
        self.log = [0] * self.q
        x = 1
        for i in range(self.q - 1):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & self.q:
                x ^= prim
        assert x == 1, "modulus not primitive"
        for i in range(self.q - 1, 2 * self.q - 1):
            self.exp[i] = self.exp[i - (self.q - 1)]

    def _find_primitive(self, seed: int) -> int:
        """find a primitive (=> irreducible) degree-m polynomial over F2"""
        rng = random.Random(seed)

        def prim_ok(p):
            # p primitive iff order of x in F2[x]/p is 2^m - 1
            order = self.q - 1
            # test: for each prime factor f of order, x^(order/f) != 1
            fs = _factor(order)
            for f in fs:
                if self._pow_mod(x=2, e=order // f, p=p) == 1:
                    return False
            return self._pow_mod(x=2, e=order, p=p) == 1

        # search: must be irreducible of degree m, highest bit set, constant 1
        cands = [p for p in range(self.q + 1, 2 * self.q, 2)]
        rng.shuffle(cands)
        for p in cands:
            if bin(p).count("1") % 2 == 1 and (p & 1):  # odd weight (const 1 needed)
                # irreducibility: x^(2^m) == x mod p, and gcd conditions; the
                # order test below implies irreducibility
                pass
            if prim_ok(p):
                return p
        raise RuntimeError("no primitive poly found")

    def _poly_mul_mod(self, a, b, p):
        # multiply bit-polys a, b mod p (p includes high bit at 2^m)
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= p
        return r

    def _pow_mod(self, x, e, p):
        r = 1
        while e:
            if e & 1:
                r = self._poly_mul_mod(r, x, p)
            x = self._poly_mul_mod(x, x, p)
            e >>= 1
        return r

    # ---------- scalar ops ----------
    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        assert a != 0, "inverse of zero"
        return self.exp[self.q - 1 - self.log[a]]

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def pow(self, a, e):
        if e == 0:
            return 1
        if a == 0:
            return 0
        return self.exp[(self.log[a] * e) % (self.q - 1)]

    def eq(self, a, b):
        return a == b

    # ---------- polynomial ops (lists, index = degree) ----------
    def padd(self, p, r):
        n = max(len(p), len(r))
        out = [0] * n
        for i in range(n):
            a = p[i] if i < len(p) else 0
            b = r[i] if i < len(r) else 0
            out[i] = a ^ b
        return out

    def pscale(self, p, s):
        if s == 0:
            return [0] * len(p)
        ls = self.log[s]
        return [0 if c == 0 else self.exp[self.log[c] + ls] for c in p]

    def pmul(self, p, r):
        if not p or not r or all(c == 0 for c in p) or all(c == 0 for c in r):
            return [0] * (len(p) + len(r) - 1)
        out = [0] * (len(p) + len(r) - 1)
        for i, a in enumerate(p):
            if a == 0:
                continue
            la = self.log[a]
            for j, b in enumerate(r):
                if b:
                    out[i + j] ^= self.exp[la + self.log[b]]
        return out

    def ptrim(self, p):
        d = len(p) - 1
        while d >= 0 and p[d] == 0:
            d -= 1
        return p[: d + 1] if d >= 0 else []

    def pdeg(self, p):
        d = len(p) - 1
        while d >= 0 and p[d] == 0:
            d -= 1
        return d

    def pdivmod(self, p, r):
        """exact poly division with remainder; r must be nonzero"""
        p = self.ptrim(p[:])
        r = self.ptrim(r[:])
        dr = self.pdeg(r)
        assert dr >= 0
        if dr == 0:
            # division by scalar
            s = self.inv(r[0])
            return self.pscale(p, s), [0]
        if self.pdeg(p) < dr:
            return [0], p
        ri = self.inv(r[dr])
        rem = p[:]
        quo = [0] * (len(p) - dr)
        for i in range(len(p) - 1, dr - 1, -1):
            c = rem[i]
            if c == 0:
                continue
            f = self.mul(c, ri)
            quo[i - dr] = f
            lf = self.log[f]
            for j in range(dr + 1):
                if r[j]:
                    rem[i - dr + j] ^= self.exp[lf + self.log[r[j]]]
        return self.ptrim(quo), self.ptrim(rem)

    def pmod(self, p, r):
        return self.pdivmod(p, r)[1]

    def pderiv(self, p):
        # formal derivative over char 2: d/dZ sum c_i Z^i = sum_{i odd} c_i Z^{i-1}
        out = [0] * max(0, len(p) - 1)
        for i in range(1, len(p)):
            if i % 2 == 1:
                # i odd: i ≡ 1 mod 2, so the scalar i is 1 in E
                out[i - 1] = p[i]
            # i even: scalar i ≡ 0, coefficient kills the term
        return self.ptrim(out)

    def peval(self, p, a):
        v = 0
        for c in reversed(p):
            v = self.mul(v, a) ^ c
        return v

    def pcompose(self, p, r):
        # compose r into p, coefficient-recorded
        out = [0]
        for c in reversed(p):
            out = self.padd(self.pmul(out, r), [c])
        return self.ptrim(out)

    def pgcd(self, p, r):
        p = self.ptrim(p[:])
        r = self.ptrim(r[:])
        while True:
            d = self.pdeg(r)
            if d < 0:
                return p
            if d == 0:
                return [1]
            _, rem = self.pdivmod(p, r)
            p, r = r, rem

    def pinv_mod(self, a, mod):
        """a^{-1} mod m for polynomials with gcd(a, mod) = 1 (extended
        Euclid); returns ptrimmed inverse with degree < deg(mod)."""
        a = self.ptrim(a[:])
        mod = self.ptrim(mod[:])
        assert self.pdeg(a) >= 0
        # ext-gcd: s,t with s*a + t*mod = 1
        r0, r1 = mod, a
        s0, s1 = [], [1]
        while True:
            if self.pdeg(r1) < 0:
                break
            q, r = self.pdivmod(r0, r1)
            r0, r1 = r1, r
            # s <- s0 - q*s1
            s = self.ptrim(self.padd(s0, self.pmul(q, s1)))
            s0, s1 = s1, s
        # gcd is a scalar (unit) since inversible; scale s0 by inv(gcd)
        if self.pdeg(r0) == 0 and r0 != []:
            g = r0[0] if r0 else 1
            s0 = self.pscale(s0, self.inv(g))
        return self.ptrim(s0)


    def pis_irreducible(self, p) -> bool:
        """Rabin test over E: Z^{q^d} = Z mod p, gcd's coprime (q = |E|)."""
        p = self.ptrim(p[:])
        d = self.pdeg(p)
        if d < 1:
            return False
        x = [0, 1]
        xq = [self.pow(c, self.q) for c in x]  # = x (F2 coeffs) — for general
        # polynomials coefficients must also be raised. Frobenius map
        # fZ = sum c_j Z^j -> sum c_j^q Z^{q j}, then reduce mod p.
        def frob_apply(u):
            out = [0]
            for j in range(len(u)):
                if u[j] == 0:
                    continue
                cq = self.frob_scalar(u[j])
                # add cq * Z^{q j}
                term = [0] * (self.q * j) + [cq]
                out = self.padd(out, term)
            return self.pmod(self.ptrim(out), p)
        xfrob = x[:]
        for _ in range(d):
            xfrob = frob_apply(xfrob)
        if self.ptrim(self.padd(xfrob, x)) != []:
            return False
        for f in _factor(d):
            e = d // f
            xf = x[:]
            for _ in range(e):
                xf = frob_apply(xf)
            g = self.pgcd(self.ptrim(self.padd(xf, x)), p)
            if self.pdeg(g) != 0:
                return False
        return True

    def frob_scalar(self, a):
        # a -> a^(2^m) in F_{2^m}; = a^(q)
        if a == 0:
            return 0
        return self.exp[(self.log[a] * self.q) % (self.q - 1)]

    def prandom(self, rng, deg, monic=False):
        coef = [rng.randrange(self.q) for _ in range(deg)] + ([1] if monic else [])
        return self.ptrim(coef)

    # ---------- matrices over E, exact Gaussian elimination ----------
    def rref(self, rows, ncols):
        """in-place row-reduce a list of rows (lists of field elements);
        returns (rank, pivot_cols, rows)"""
        rows = [r[:] for r in rows]
        nr = len(rows)
        piv = []
        r = 0
        for c in range(ncols):
            pr = None
            for i in range(r, nr):
                if rows[i][c] != 0:
                    pr = i
                    break
            if pr is None:
                continue
            rows[r], rows[pr] = rows[pr], rows[r]
            iv = self.inv(rows[r][c])
            rows[r] = self.pscale(rows[r], iv) if False else [self.mul(v, iv) for v in rows[r]]
            for i in range(nr):
                if i != r and rows[i][c] != 0:
                    f = rows[i][c]
                    lf = self.log[f]
                    rows[i] = [rows[i][j] ^ (0 if rows[r][j] == 0 else self.exp[lf + self.log[rows[r][j]]]) for j in range(ncols)]
            piv.append(c)
            r += 1
            if r == nr:
                break
        return r, piv, rows

    def nullspace(self, rows, ncols):
        """basis (as lists) of {x : M x = 0} for matrix M given as rows"""
        rank, piv, red = self.rref(rows, ncols)
        pivset = set(piv)
        free = [c for c in range(ncols) if c not in pivset]
        basis = []
        for fc in free:
            v = [0] * ncols
            v[fc] = 1
            for ri, pc in enumerate(piv):
                if red[ri][fc] != 0:
                    v[pc] = red[ri][fc]
            basis.append(v)
        return basis, rank, piv

    def rank(self, rows, ncols):
        return self.rref(rows, ncols)[0]


def _factor(n):
    fs = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            fs.append(d)
            n //= d
        d += 1
    if n > 1:
        fs.append(n)
    return fs
