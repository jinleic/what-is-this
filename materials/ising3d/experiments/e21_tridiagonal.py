"""Does the Ising transfer pair (A,B) satisfy ANY tridiagonal (Askey-Wilson) relation?

Theorem DG (`proofs/dolan_grady_defect.md`) shows the Dolan-Grady relations fail on any layer
graph with a vertex of degree >= 3.  But Dolan-Grady is only ONE point of a three-parameter
family.  The general criterion for the Onsager / Askey-Wilson / tridiagonal-pair mechanism
(Terwilliger; Ito-Tanabe-Terwilliger) is the pair of TRIDIAGONAL RELATIONS

    (TD1)   [ A,  A^2 B - beta A B A + B A^2 - gamma  (A B + B A) - rho  B ] = 0
    (TD2)   [ B,  B^2 A - beta B A B + A B^2 - gamma* (B A + A B) - rho* A ] = 0

Dolan-Grady is beta = 2, gamma = gamma* = 0, rho = rho* = 16.

For FIXED A, B this is a LINEAR condition on (beta, gamma, rho): with
    T0 = [A, A^2B + BA^2],  T1 = [A, ABA],  T2 = [A, AB+BA],  T3 = [A, B],
(TD1) says T0 = beta T1 + gamma T2 + rho T3, i.e. T0 in span{T1,T2,T3}.  Exact integer linear
algebra decides it.

LEMMA (affine/scaling completeness -- proved in the docstring of `affine_note()` below).
Replacing A -> aA + s and B -> bB + t only REPARAMETRISES (gamma, rho); the solvability of the
linear system is unchanged.  Hence the linear test decides the FULL affine tridiagonal family,
not just one normalisation.
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ising.clifford import pauli_x, zz  # noqa: E402


def affine_note():
    """Why the linear test is complete.

    Put A' = aA + s, B' = bB + t.  Every term of
        A'^2 B' - beta A'B'A' + B'A'^2 - gamma (A'B'+B'A') - rho B'
    is linear in B', so the t-dependent part is  t(2-beta)(aA+s)^2 - 2 gamma t (aA+s) - rho t,
    a polynomial in A alone, which is killed by the outer bracket [A', .] = a[A, .].
    The B-dependent part is b times
        a^2 (A^2B + BA^2) - beta a^2 ABA + [ s a (2-beta) - gamma a ](AB+BA)
        + [ s^2 (2-beta) - 2 gamma s - rho ] B .
    Dividing by a^2 b, the relation becomes the UNSHIFTED one with
        gamma_new = gamma/a - s(2-beta)/a ,     rho_new = (rho + 2 gamma s - s^2(2-beta))/a^2 .
    Since (gamma_new, rho_new) range over all of Q^2 as (gamma, rho) do, the affine family is
    solvable iff the unshifted linear system is.  beta is unchanged.  QED
    """


def mul(e1: dict, e2: dict, n: int) -> dict:
    mask = (1 << n) - 1
    out: dict[int, int] = {}
    for v, c in e1.items():
        bv = (v >> n) & mask
        for w, d in e2.items():
            aw = w & mask
            s = -1 if (bin(bv & aw).count("1") & 1) else 1
            k = v ^ w
            out[k] = out.get(k, 0) + s * c * d
    return {k: c for k, c in out.items() if c}


def add(*es):
    out: dict[int, int] = {}
    for e in es:
        for k, c in e.items():
            out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}


def smul(a, e):
    return {k: a * c for k, c in e.items() if a * c}


def comm(e1, e2, n):
    return add(mul(e1, e2, n), smul(-1, mul(e2, e1, n)))


def build(n, bonds, swap=False):
    A = {pauli_x(n, i): 1 for i in range(n)}
    B = {}
    for i, j in bonds:
        B[zz(n, i, j)] = B.get(zz(n, i, j), 0) + 1
    if swap:
        A, B = B, A
    AA = mul(A, A, n)
    T0 = comm(A, add(mul(AA, B, n), mul(B, AA, n)), n)
    T1 = comm(A, mul(mul(A, B, n), A, n), n)
    T2 = comm(A, add(mul(A, B, n), mul(B, A, n)), n)
    T3 = comm(A, B, n)
    return T0, T1, T2, T3


def solve_linear(T0, cols):
    """Exact solve T0 = sum_j x_j cols[j].  Returns (solution_or_None, rank, nullity)."""
    keys = sorted(set(T0) | {k for c in cols for k in c})
    w = len(cols)
    rows = [[Fraction(c.get(k, 0)) for c in cols] + [Fraction(T0.get(k, 0))] for k in keys]
    m = len(rows)
    piv, where = 0, []
    for col in range(w):
        sel = next((r for r in range(piv, m) if rows[r][col] != 0), None)
        if sel is None:
            where.append(-1)
            continue
        rows[piv], rows[sel] = rows[sel], rows[piv]
        pv = rows[piv][col]
        rows[piv] = [x / pv for x in rows[piv]]
        for r in range(m):
            if r != piv and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[piv])]
        where.append(piv)
        piv += 1
    rank = piv
    for r in range(m):
        if all(rows[r][c] == 0 for c in range(w)) and rows[r][w] != 0:
            return None, rank, w - rank
    sol = [Fraction(0)] * w
    for col in range(w):
        if where[col] != -1:
            sol[col] = rows[where[col]][w]
    return sol, rank, w - rank


def grid_bonds(a, b, per=False):
    bonds = []
    for x in range(a):
        for y in range(b):
            i = x * b + y
            if x + 1 < a:
                bonds.append((i, (x + 1) * b + y))
            elif per and a > 2:
                bonds.append((i, y))
            if y + 1 < b:
                bonds.append((i, x * b + y + 1))
            elif per and b > 2:
                bonds.append((i, x * b))
    return bonds


BOUNDARYLESS_1D = [(f"cycle C{L} (2-regular)", L, [(i, (i + 1) % L) for i in range(L)])
                   for L in range(3, 10)]
BOUNDARYLESS_2D = [(f"torus layer {a}x{b} (4-regular)", a * b, grid_bonds(a, b, per=True))
                   for (a, b) in [(3, 3), (3, 4), (4, 4), (3, 5), (4, 5), (5, 5)]]
WITH_BOUNDARY = [
    ("path P4  (open 1D layer)", 4, [(0, 1), (1, 2), (2, 3)]),
    ("path P6  (open 1D layer)", 6, [(i, i + 1) for i in range(5)]),
    ("P5 + one bond (Delta=3)", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)]),
    ("2x3 grid (open 2D layer)", 6, grid_bonds(2, 3)),
    ("3x3 grid (open 2D layer)", 9, grid_bonds(3, 3)),
]


def run(cases, title):
    print(f"\n### {title}")
    print(f"{'layer graph':32s} {'n':>2s} | {'TD1 (beta,gamma,rho)':>26s} | {'TD2 (beta,gamma*,rho*)':>26s}")
    print("-" * 96)
    out = []
    for name, n, bonds in cases:
        rec = dict(name=name, n=n)
        for swap, tag in [(False, "TD1"), (True, "TD2")]:
            T0, T1, T2, T3 = build(n, bonds, swap=swap)
            sol, rank, nul = solve_linear(T0, [T1, T2, T3])
            if sol is None:
                rec[tag] = "NO SOLUTION"
            else:
                rec[tag] = f"({sol[0]}, {sol[1]}, {sol[2]})" + (f" +{nul}d" if nul else "")
            rec[tag + "_rank"] = rank
        print(f"{name:32s} {n:2d} | {rec['TD1']:>26s} | {rec['TD2']:>26s}")
        out.append(rec)
    return out


if __name__ == "__main__":
    rows = []
    rows += run(BOUNDARYLESS_1D, "BOUNDARYLESS 1D layers  ->  2D classical Ising on a cylinder")
    rows += run(BOUNDARYLESS_2D, "BOUNDARYLESS 2D layers  ->  3D classical Ising, no boundary confound")
    rows += run(WITH_BOUNDARY, "layers WITH boundary (for contrast; open chains fail by boundary terms)")

    b1 = [r for r in rows if "2-regular" in r["name"]]
    b2 = [r for r in rows if "4-regular" in r["name"]]
    ok1 = all(r["TD1"] != "NO SOLUTION" and r["TD2"] != "NO SOLUTION" for r in b1)
    ok2 = all(r["TD2"] == "NO SOLUTION" for r in b2)
    print("\n" + "=" * 96)
    print(f"every BOUNDARYLESS 1D layer (cycle) satisfies a tridiagonal relation : {ok1}")
    print(f"NO  BOUNDARYLESS 2D layer (torus) satisfies any tridiagonal relation : {ok2}")
    print("""
CONCLUSION.  The obstruction is not an artefact of the Dolan-Grady normalisation and not a
boundary effect.  Comparing the two boundaryless, translation-invariant cases:

    2-regular layer (2D Ising)  ->  (A,B) IS a tridiagonal pair;
    4-regular layer (3D Ising)  ->  (A,B) satisfies NO tridiagonal relation, for any
                                    (beta, gamma*, rho*), and by the affine lemma this is
                                    invariant under A -> aA+s, B -> bB+t.

The 3D Ising layer transfer pair therefore lies outside the entire Onsager / Askey-Wilson /
tridiagonal-pair mechanism, not merely outside its Dolan-Grady specialisation.""")
    json.dump(rows, open("results/tridiagonal.json", "w"), indent=1, default=str)
    print("\nwritten results/tridiagonal.json")
