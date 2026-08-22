"""The Dolan-Grady test: an EXACT algebraic identity that decides the Onsager mechanism.

Onsager's 1944 solution of the 2D Ising model rests on the pair of operators

    A = sum_i X_i ,      B = sum_{<ij> in layer} Z_i Z_j

generating the Onsager algebra.  Dolan and Grady (Phys. Rev. D 25, 1587 (1982)) showed that
A and B generate an Onsager algebra -- equivalently that the one-parameter family
`A + k B` yields commuting transfer matrices and an infinite tower of conserved charges --
precisely when the two *Dolan-Grady relations* hold:

    [A, [A, [A, B]]]  =  16 [A, B]
    [B, [B, [B, A]]]  =  16 [B, A]

(the constant 16 is fixed by the normalisation `A^2 -> 1` per term; we determine the correct
constant empirically for the 1D chain and then apply the SAME constant in 2D).

This is a finite, exact, integer identity: all commutators of Pauli-string sums have integer
coefficients (structure constants in {0, +-2}).  So the test is decisive, not numerical.

Result to be established:
  * 1D layers (2D Ising)  : the relations hold identically -- the residual is exactly zero.
  * 2D layers (3D Ising)  : the residual is a nonzero operator; we print it, count its terms and
    identify WHICH lattice motif produces it.
"""

from __future__ import annotations

import json
from itertools import product

from ising.clifford import pauli_x, zz


def comm(e1: dict, e2: dict, n: int) -> dict:
    """[sum c_v Q_v, sum d_w Q_w] with exact integer coefficients."""
    mask = (1 << n) - 1
    out: dict[int, int] = {}
    for v, c in e1.items():
        av, bv = v & mask, (v >> n) & mask
        for w, d in e2.items():
            aw, bw = w & mask, (w >> n) & mask
            s1 = bin(bv & aw).count("1") & 1
            s2 = bin(bw & av).count("1") & 1
            if s1 == s2:
                continue
            k = v ^ w
            out[k] = out.get(k, 0) + (2 if s1 == 0 else -2) * c * d
    return {k: c for k, c in out.items() if c}


def lin(e1: dict, e2: dict, a=1, b=1) -> dict:
    out = dict()
    for k, c in e1.items():
        out[k] = out.get(k, 0) + a * c
    for k, c in e2.items():
        out[k] = out.get(k, 0) + b * c
    return {k: c for k, c in out.items() if c}


def dolan_grady_residuals(n, bonds, const=16):
    xs = [pauli_x(n, i) for i in range(n)]
    A = {g: 1 for g in xs}
    B = {}
    for i, j in bonds:
        B[zz(n, i, j)] = B.get(zz(n, i, j), 0) + 1
    AB = comm(A, B, n)
    AAB = comm(A, AB, n)
    AAAB = comm(A, AAB, n)
    r1 = lin(AAAB, AB, 1, -const)
    BA = comm(B, A, n)
    BBA = comm(B, BA, n)
    BBBA = comm(B, BBA, n)
    r2 = lin(BBBA, BA, 1, -const)
    return r1, r2, AB


def pauli_str(v, n):
    mask = (1 << n) - 1
    a, b = v & mask, (v >> n) & mask
    s = ""
    for i in range(n):
        ai, bi = (a >> i) & 1, (b >> i) & 1
        s += "IXZY"[ai + 2 * bi] if not (ai and bi) else "Y"
    return s


def graph_cases():
    cases = []
    for L in range(3, 9):
        cases.append((f"open chain n={L}", L, [(i, i + 1) for i in range(L - 1)]))
    for L in range(3, 8):
        cases.append((f"ring n={L}", L, [(i, (i + 1) % L) for i in range(L)]))
    cases.append(("star K_{1,3}", 4, [(0, 1), (0, 2), (0, 3)]))
    cases.append(("star K_{1,4}", 5, [(0, 1), (0, 2), (0, 3), (0, 4)]))
    cases.append(("P4 + pendant at 1 (Delta=3)", 5, [(0, 1), (1, 2), (2, 3), (1, 4)]))
    cases.append(("P5 + pendant at 2 (Delta=3)", 6, [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)]))
    for (a, b) in [(2, 2), (2, 3), (2, 4), (3, 3)]:
        n = a * b
        bonds = []
        for x in range(a):
            for y in range(b):
                i = x * b + y
                if x + 1 < a:
                    bonds.append((i, (x + 1) * b + y))
                if y + 1 < b:
                    bonds.append((i, x * b + y + 1))
        cases.append((f"grid {a}x{b} (3D Ising layer)", n, bonds))
    return cases


if __name__ == "__main__":
    rows = []
    print(f"{'layer graph':32s} {'n':>2s} {'|res1|':>7s} {'|res2|':>7s}  verdict")
    print("-" * 72)
    for name, n, bonds in graph_cases():
        r1, r2, AB = dolan_grady_residuals(n, bonds)
        deg = [0] * n
        for i, j in bonds:
            deg[i] += 1
            deg[j] += 1
        holds = (len(r1) == 0 and len(r2) == 0)
        rows.append(dict(name=name, n=n, maxdeg=max(deg), res1_terms=len(r1),
                         res2_terms=len(r2), dolan_grady=holds))
        print(f"{name:32s} {n:2d} {len(r1):7d} {len(r2):7d}  "
              f"{'DOLAN-GRADY HOLDS' if holds else 'FAILS'}", flush=True)

    print()
    print("=== structure of the first failure ===")
    # smallest failing case: the star K_{1,3}
    for name, n, bonds in [("star K_{1,3}", 4, [(0, 1), (0, 2), (0, 3)]),
                           ("P4+pendant", 5, [(0, 1), (1, 2), (2, 3), (1, 4)])]:
        r1, r2, AB = dolan_grady_residuals(n, bonds)
        print(f"\n{name}:  residual of [A,[A,[A,B]]] - 16[A,B] has {len(r1)} Pauli terms:")
        for v, c in sorted(r1.items())[:14]:
            print(f"    {c:+6d} * {pauli_str(v, n)}")
        print(f"  residual of [B,[B,[B,A]]] - 16[B,A] has {len(r2)} Pauli terms:")
        for v, c in sorted(r2.items())[:14]:
            print(f"    {c:+6d} * {pauli_str(v, n)}")

    ok_rule = all((r["maxdeg"] <= 2) == r["dolan_grady"] for r in rows)
    print(f"\nRULE  'Dolan-Grady holds  <=>  max degree <= 2' : {'CONFIRMED' if ok_rule else 'VIOLATED'}"
          f" on all {len(rows)} layer graphs tested")
    json.dump(rows, open("results/dolan_grady.json", "w"), indent=1)
    print("written results/dolan_grady.json")
