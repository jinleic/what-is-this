"""Exact F_q arithmetic and subgroup construction. No floats anywhere."""
from __future__ import annotations

from fractions import Fraction

import sympy
from sympy.ntheory.residue_ntheory import primitive_root as _sympy_primitive_root


def primitive_root(q: int) -> int:
    """Primitive root of F_q^x for prime q (sympy returns the smallest)."""
    return int(_sympy_primitive_root(q))


def divisors(n: int) -> list[int]:
    ds = set()
    d = 1
    while d * d <= n:
        if n % d == 0:
            ds.add(d)
            ds.add(n // d)
        d += 1
    return sorted(ds)


def subgroup(q: int, s: int, g: int) -> list[int]:
    """Elements of the subgroup of F_q^x of order s (s | q-1), sorted ascending.

    Generated as <g^((q-1)/s)>.
    """
    assert s >= 1 and (q - 1) % s == 0, (q, s)
    h = pow(g, (q - 1) // s, q)
    els = [1]
    el = h
    while el != 1:
        els.append(el)
        el = (el * h) % q
    assert len(els) == s, (q, s, len(els))
    return sorted(els)


def coprime_triples(n: int, smin: int = 2) -> list[tuple[int, int, int]]:
    """All triples (s0,s1,s2) of divisors of n, pairwise coprime, each >= smin."""
    out = []
    ds = [d for d in divisors(n) if d >= smin]
    for a in ds:
        for b in ds:
            if b <= a:
                continue
            for c in ds:
                if c <= b:
                    continue
                from math import gcd
                if gcd(a, b) == 1 and gcd(a, c) == 1 and gcd(b, c) == 1:
                    out.append((a, b, c))
    return out


# ---------------- exact F_q linear algebra ----------------

def rref(mat: list[list[int]], q: int) -> tuple[list[list[int]], list[int]]:
    """Row-reduce mod q (row echelon, unit pivots). Returns (rows, pivot_cols)."""
    rows = [[x % q for x in r] for r in mat if any(r)]
    ncol = len(mat[0]) if mat else 0
    res: list[list[int]] = []
    piv: list[int] = []
    r = 0
    for c in range(ncol):
        p = None
        for i in range(r, len(rows)):
            if rows[i][c]:
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = pow(rows[r][c], -1, q)
        rr = rows[r]
        for j in range(c, ncol):
            rr[j] = (rr[j] * inv) % q
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c]
                ri = rows[i]
                for j in range(c, ncol):
                    ri[j] = (ri[j] - f * rr[j]) % q
        res.append(rows[r])
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return res, piv


def rank_mod(rows: list[list[int]], q: int) -> int:
    if not rows:
        return 0
    rr, piv = rref(rows, q)
    return len(piv)


def kernel_mod(rows: list[list[int]], q: int, ncol: int) -> list[list[int]]:
    """Basis of {x : rows @ x = 0} over F_q, in terms of unit pivots/free vars."""
    rr, piv = rref(rows, q) if rows else ([], [])
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * ncol
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-rr[i][fc]) % q
        basis.append(v)
    return basis


def in_span(rows: list[list[int]], vec: list[int], q: int) -> bool:
    """Is vec in the F_q-span of rows? Exact."""
    if not rows:
        return not any(vec)
    stacked = [list(r) for r in rows] + [list(vec)]
    rr, piv = rref(stacked, q)
    # vec in span iff while eliminating, the last row became zero
    # equivalent: rank(stacked) == rank(rows)
    rr2, piv2 = rref(rows, q)
    return len(piv) == len(piv2)


def solve_mod(rows: list[list[int]], rhs: list[int], q: int):
    """Solve rows @ x = rhs over F_q. Returns one solution list or None."""
    ncol = len(rows[0]) if rows else len(rhs)
    aug = [list(r) + [b % q] for r, b in zip(rows, rhs)]
    rr, piv = rref(aug, q)
    # inconsistency: nonzero row with all-zero left part
    for r in rr:
        if all(v == 0 for v in r[:ncol]) and r[ncol] % q:
            return None
    x = [0] * ncol
    for i, pc in enumerate(piv):
        x[pc] = rr[i][ncol] % q
    # verify
    for r, b in zip(rows, rhs):
        s = sum(a * xi for a, xi in zip(r, x)) % q
        assert s == b % q
    return x


def frac_str(fr: Fraction) -> str:
    return f"{fr.numerator}/{fr.denominator}"


def frac_cmp(a: Fraction, b: Fraction) -> int:
    return (a > b) - (a < b)
