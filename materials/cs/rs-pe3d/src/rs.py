"""RS tensor-code construction and line-space machinery (exact F_q).

An Inst fixes q, s=(s0,s1,s2), t=(t0,t1,t2), and optional diagonal rescalings
lam[i] (length s[i], nonzero entries) playing Lambda_i. Coordinate order:
    idx = (a0*s1 + a1)*s2 + a2, with a_i an index into S_i.
An i-axis LINE fixes all coordinates except a_i and varies a_i.
"""
from __future__ import annotations

import itertools

from .field import subgroup, primitive_root


class Inst:
    def __init__(self, q: int, s: tuple, t: tuple, lam=None):
        self.q = q
        self.s = tuple(int(x) for x in s)
        self.t = tuple(int(x) for x in t)
        assert self.q <= self.N_bound_guard()  # guard against silly misuse
        self.N = self.s[0] * self.s[1] * self.s[2]
        self.g = primitive_root(q)
        self.S = [subgroup(q, self.s[i], self.g) for i in range(3)]
        if lam is None:
            self.lam = [tuple([1] * self.s[i]) for i in range(3)]
        else:
            self.lam = [tuple(int(x) % q for x in l) for l in lam]
            for i in range(3):
                assert len(self.lam[i]) == self.s[i]
                assert all(self.lam[i])

    def N_bound_guard(self):
        return 10 ** 9

    # ---------------- index helpers ----------------

    def idx(self, a: tuple) -> int:
        return (a[0] * self.s[1] + a[1]) * self.s[2] + a[2]

    def coord_of(self, flat: int) -> tuple:
        a2 = flat % self.s[2]
        rest = flat // self.s[2]
        a1 = rest % self.s[1]
        a0 = rest // self.s[1]
        return (a0, a1, a2)

    def line_indices(self, i: int, x: int) -> list[int]:
        """Line of direction i indexed by other-axis flat index x in
        [0, N/s_i): the OTHER coordinates in row-major order of axes < i, > i.
        Returns list of flat indices (length s[i])."""
        dim_i = self.s[i]
        front = self.s[0] if i > 0 else 1
        # decompose x over the other axes in axis order (skip i)
        others = [self.s[j] for j in range(3) if j != i]
        r0 = x % others[0]
        r1 = x // others[0]
        aa = [0, 0, 0]
        pos = 0
        out = []
        # rebuild the other coords
        oc = []
        oc.append(r0)
        oc.append(r1 % others[1])
        # oc[0] -> axis 0 if 0 in others else axis1; general: axis order
        axes = [j for j in range(3) if j != i]
        jj = 0
        for j in axes:
            aa[j] = oc[jj]
            jj += 1
        for ai in range(dim_i):
            aa[i] = ai
            out.append(self.idx(aa))
        return out

    def lines(self, i: int):
        for x in range(self.N // self.s[i]):
            yield self.line_indices(i, x)

    def line_of_point(self, i: int, flat: int) -> int:
        """Other-axis index of the i-line through flat point."""
        c = self.coord_of(flat)
        axes = [j for j in range(3) if j != i]
        others = [self.s[j] for j in range(3) if j != i]
        oc = [c[axes[0]], c[axes[1]]]
        return oc[0] + oc[1] * others[0]

    def num_lines(self, i: int) -> int:
        return self.N // self.s[i]

    # ---------------- line spaces (Lambda-scaled RS), exact ----------------

    def line_space(self, i: int) -> list[list[int]]:
        """Basis of the direction-i line space P_i in F_q^{s_i}:
        rows = Lambda_i(a_alpha) * alpha^d for d = 0..t_i-1, alpha over S_i
        in index order."""
        q, t, S = self.q, self.t[i], self.S[i]
        rows = []
        for d in range(t):
            rows.append([self.lam[i][ai] * pow(S[ai], d, q) % q
                         for ai in range(self.s[i])])
        return rows

    def line_basis_flat(self, i: int, line_other_idx: int) -> list[list[int]]:
        """Basis vectors (flat, in F_q^N) for one line of direction i."""
        idxs = self.line_indices(i, line_other_idx)
        out = []
        for v in self.line_space(i):
            w = [0] * self.N
            for j, flat in enumerate(idxs):
                w[flat] = v[j] % self.q
            out.append(w)
        return out

    def lift_basis(self, i: int) -> list[list[int]]:
        out = []
        for x in range(self.num_lines(i)):
            out.extend(self.line_basis_flat(i, x))
        return out
