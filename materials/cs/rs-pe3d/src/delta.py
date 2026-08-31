"""Exact delta(M) and exact V-membership algebra (pilot-scale exhaustive core).

All exact F_q arithmetic. delta(M) = min sum_i s_i |B_i| over line-support
tuples B with M in W(B). Computed by enumerating ALL compositions (n_0,n_1,n_2)
in ASCENDING s-weighted cost and, for each, ALL B_i of sizes n_i, testing
membership M in W(B) exactly. First hit = exact delta (all cheaper costs were
fully enumerated). Optional cutoff: returns None if exceeded.
"""
from __future__ import annotations

import itertools

from .field import rref, rank_mod


class DeltaEngine:
    """Caches per-line bases and membership subroutines for one Inst."""

    def __init__(self, inst):
        self.inst = inst
        self.q = inst.q
        # per direction: list over lines of reduced basis (rows in F_q^N)
        self.line_rows: list[list[list[list[int]]]] = []
        for i in range(3):
            per = []
            for l in range(inst.num_lines(i)):
                rows = inst.line_basis_flat(i, l)
                rr, _ = rref(rows, self.q)
                per.append(rr)
            self.line_rows.append(per)
        # compositions sorted by ascending weighted cost
        s = inst.s
        comps = []
        cap = [inst.num_lines(i) for i in range(3)]
        for n0 in range(cap[0] + 1):
            for n1 in range(cap[1] + 1):
                base = s[0] * n0 + s[1] * n1
                if base > 3 * inst.N:
                    continue
                for n2 in range(cap[2] + 1):
                    c = base + s[2] * n2
                    if c <= 3 * inst.N:
                        comps.append((c, n0, n1, n2))
        comps.sort()
        self.comps = comps

    def wt(self, M):
        return sum(1 for x in M if x % self.q)

    def _membership(self, rows_acc: list[list[int]], M) -> bool:
        """Exact membership of M in span(rows_acc)."""
        stacked = [list(r) for r in rows_acc] + [list(M)]
        _, piv = rref(stacked, self.q)
        _, piv0 = rref(rows_acc, self.q)
        return len(piv) == len(piv0)

    def delta_exact(self, M, max_cost: int | None = None):
        """Exact delta by exhaustive ascending-cost composition sweep.

        Returns (delta, (B0,B1,B2)) or (None, None) if max_cost exceeded
        without success.
        """
        q = self.q
        for cost, n0, n1, n2 in self.comps:
            if cost == 0:
                continue
            if max_cost is not None and cost > max_cost:
                return None, None
            for B0 in itertools.combinations(range(self.inst.num_lines(0)), n0):
                rows0 = [r for l in B0 for r in self.line_rows[0][l]]
                r0k = rank_mod(rows0, q) if rows0 else 0
                for B1 in itertools.combinations(range(self.inst.num_lines(1)), n1):
                    rows01 = rows0 + [r for l in B1 for r in self.line_rows[1][l]]
                    r01k = rank_mod(rows01, q) if n1 else r0k
                    for B2 in itertools.combinations(range(self.inst.num_lines(2)), n2):
                        if n2:
                            rows = rows01 + [r for l in B2
                                             for r in self.line_rows[2][l]]
                        else:
                            rows = rows01
                        if self._membership(rows, M):
                            return cost, (B0, B1, B2)
        return None, None

    # ---------- exact V-intersection algebra ----------

    def lift_rref(self, i: int):
        """Reduced basis of L_i(C_i) (cached rref over all its line bases)."""
        if not hasattr(self, "_lift_cache"):
            self._lift_cache = {}
        if i not in self._lift_cache:
            rows = []
            for l in range(self.inst.num_lines(i)):
                rows.extend(self.line_rows[i][l])
            self._lift_cache[i] = rref(rows, self.q)[0]
        return self._lift_cache[i]

    def V_basis(self) -> list[list[int]]:
        if not hasattr(self, "_V_basis"):
            rows = []
            for i in range(3):
                rows.extend(self.lift_rref(i))
            self._V_basis = rref(rows, self.q)[0]
        return self._V_basis

    def V_dim(self) -> int:
        return len(self.V_basis())

    def V_inter_support(self, S: set[int]) -> list[list[int]]:
        """Reduced basis of V ∩ F^S (elements of V supported in point set S).

        Exact: compute kernel of V-basis restricted+projected. Method: a
        vector v ∈ V with supp(v) ⊆ S iff v ∈ V and v has zeros off S:
        solve  V-basis * x = v  with v_off = 0. Parametrize v = B^T x.
        Implement via: stack = rows of V-basis; free-support constraint:
        v = sum x_j b_j; require v[c] = 0 for c ∉ S: linear system
        (B_off) x = 0 where B_off = columns of basis off S; kernel of B_off
        gives coefficient vectors x; then v = x @ basis.
        """
        q = self.q
        basis = self.V_basis()
        off = [c for c in range(self.inst.N) if c not in S]
        B_off = [[basis[j][c] for c in off] for j in range(len(basis))]
        # kernel of B_off (as matrix with rows = basis vectors, cols = off)
        # we need x with B_off @ ... wait B_off[j][c'] = b_j[off[c']]:
        # (x @ B_off)[c'] = 0 -> x in left-kernel of B_off.
        return self._left_kernel(B_off, basis)

    def _left_kernel(self, mat: list[list[int]], basis) -> list[list[int]]:
        """{x : x @ mat = 0} as coefficient vectors; return the V-elements."""
        q = self.q
        # There are no equations when S is the whole point set.
        if not mat or not mat[0]:
            ker = [[1 if j == k else 0 for j in range(len(basis))]
                   for k in range(len(basis))]
        else:
            # transpose problem: x @ mat = 0 <=> mat^T @ x^T = 0
            rows = [[mat[j][c] for j in range(len(mat))]
                    for c in range(len(mat[0]))]
            from .field import kernel_mod
            ker = kernel_mod(rows, q, len(mat))
        out = []
        for x in ker:
            v = [0] * self.inst.N
            for j, coef in enumerate(x):
                if coef:
                    bj = basis[j]
                    for c in range(self.inst.N):
                        v[c] = (v[c] + coef * bj[c]) % q
            if any(v):
                out.append(v)
        if out:
            from .field import rref as _rref
            out = _rref(out, q)[0]
        return out
