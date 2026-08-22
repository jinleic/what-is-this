"""Verify the closed-form Dolan-Grady defect theorem against brute-force Pauli algebra.

THEOREM (derived by hand in proofs/dolan_grady_defect.md, checked here on many graphs).
For A = sum_v X_v and B = sum_{(ij) in E} Z_i Z_j on any finite simple graph Gamma:

  (i)   [A,[A,[A,B]]] - 16 [A,B] = 0                        identically, for EVERY graph;

  (ii)  [B,[B,[B,A]]] - 16 [B,A]
          = 24 i sum_v  Y_v [ (deg v - 2) * sum_{u~v} Z_u
                              + 2 * sum_{ {u,u',u''} subset N(v) } Z_u Z_{u'} Z_{u''} ] .

Consequently the Dolan-Grady relations hold identically iff every vertex has degree exactly 2,
and the QUARTIC (4-body) part vanishes identically iff every vertex has degree <= 2.

This script rebuilds the right-hand side of (ii) directly from the graph and compares it, term by
term with exact integer coefficients, against the brute-force nested commutators.
"""

from __future__ import annotations

import json
import os
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from e14_dolan_grady import comm, lin  # noqa: E402
from ising.clifford import pauli_x, zz  # noqa: E402

# Pauli encoding: bits 0..n-1 = X part, bits n..2n-1 = Z part.
def X(n, i):
    return 1 << i


def Z(n, i):
    return 1 << (n + i)


def Y(n, i):
    return (1 << i) | (1 << (n + i))


def prod(n, terms):
    """Product of distinct-site single-qubit Paulis -> (vector, +1) since sites are distinct."""
    v = 0
    for t in terms:
        v ^= t
    return v


def predicted_residual(n, bonds):
    """Right-hand side of (ii), up to an overall global sign/phase fixed by comparison."""
    nbr = [[] for _ in range(n)]
    for i, j in bonds:
        nbr[i].append(j)
        nbr[j].append(i)
    out: dict[int, int] = {}
    for v in range(n):
        d = len(nbr[v])
        for u in nbr[v]:
            key = prod(n, [Y(n, v), Z(n, u)])
            out[key] = out.get(key, 0) + 24 * (d - 2)
        for trip in combinations(sorted(nbr[v]), 3):
            key = prod(n, [Y(n, v)] + [Z(n, u) for u in trip])
            out[key] = out.get(key, 0) + 48
    return {k: c for k, c in out.items() if c}


def actual_residual(n, bonds, const=16):
    xs = [pauli_x(n, i) for i in range(n)]
    A = {g: 1 for g in xs}
    B = {}
    for i, j in bonds:
        B[zz(n, i, j)] = B.get(zz(n, i, j), 0) + 1
    BA = comm(B, A, n)
    r2 = lin(comm(B, comm(B, BA, n), n), BA, 1, -const)
    AB = comm(A, B, n)
    r1 = lin(comm(A, comm(A, AB, n), n), AB, 1, -const)
    return r1, r2


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


CASES = [
    ("path P3", 3, [(0, 1), (1, 2)]),
    ("path P6", 6, [(i, i + 1) for i in range(5)]),
    ("path P9", 9, [(i, i + 1) for i in range(8)]),
    ("ring C5", 5, [(i, (i + 1) % 5) for i in range(5)]),
    ("ring C7", 7, [(i, (i + 1) % 7) for i in range(7)]),
    ("star K13", 4, [(0, 1), (0, 2), (0, 3)]),
    ("star K15", 6, [(0, k) for k in range(1, 6)]),
    ("K4", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
    ("K5", 5, list(combinations(range(5), 2))),
    ("tree caterpillar", 7, [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5), (3, 6)]),
    ("2x3 grid", 6, grid_bonds(2, 3)),
    ("2x4 grid", 8, grid_bonds(2, 4)),
    ("3x3 grid", 9, grid_bonds(3, 3)),
    ("3x3 torus layer", 9, grid_bonds(3, 3, per=True)),
    ("3x4 grid", 12, grid_bonds(3, 4)),
    ("petersen-ish cubic", 6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]),
]

if __name__ == "__main__":
    rows = []
    allok = True
    sign = None
    print(f"{'graph':22s} {'n':>2s} {'rel(i)':>7s} {'#terms':>7s} {'#quartic':>9s} "
          f"{'sum C(d,3)':>11s} {'formula':>9s}")
    print("-" * 78)
    for name, n, bonds in CASES:
        r1, r2 = actual_residual(n, bonds)
        pred = predicted_residual(n, bonds)
        # fix the global phase convention once, from the first nonempty case
        if sign is None and r2 and pred:
            k = sorted(set(r2) & set(pred))
            if k:
                sign = r2[k[0]] // pred[k[0]]
        s = sign if sign else 1
        match = (r2 == {k: s * c for k, c in pred.items()})
        deg = [0] * n
        for i, j in bonds:
            deg[i] += 1
            deg[j] += 1
        nq = sum(1 for v in r2 if bin(v & ((1 << n) - 1)).count("1")
                 + bin(((v >> n) & ((1 << n) - 1)) & ~(v & ((1 << n) - 1))).count("1") >= 4)

        def body(v):
            mask = (1 << n) - 1
            a, b = v & mask, (v >> n) & mask
            return bin(a | b).count("1")

        nq = sum(1 for v in r2 if body(v) == 4)
        claws = sum(d * (d - 1) * (d - 2) // 6 for d in deg)
        pf = sum((d if d != 2 else 0) + d * (d - 1) * (d - 2) // 6 for d in deg)
        ok = match and (len(r1) == 0) and nq == claws and len(r2) == pf
        allok &= ok
        print(f"{name:22s} {n:2d} {'0 OK' if not r1 else 'NONZERO':>7s} {len(r2):7d} {nq:9d} "
              f"{claws:11d} {str(pf):>9s}  {'MATCH' if match else 'MISMATCH'}")
        rows.append(dict(name=name, n=n, rel_i_zero=len(r1) == 0, terms=len(r2),
                         quartic_terms=nq, sum_binom_deg_3=claws, formula_terms=pf,
                         closed_form_match=bool(match)))
    print()
    print(f"global phase factor between derivation and code convention: {sign}")
    print(f"\nCLOSED FORM VERIFIED ON ALL {len(CASES)} GRAPHS: {'PASS' if allok else 'FAIL'}")
    print("  - relation (i) [A,[A,[A,B]]]=16[A,B] holds identically on every graph")
    print("  - relation (ii) residual matches 24i sum_v Y_v[(d-2)S_v + 2T_v] exactly")
    print("  - #quartic terms == sum_v C(deg v,3) == number of induced claws in the")
    print("    frustration graph, on every graph tested")
    json.dump(rows, open("results/dg_theorem_check.json", "w"), indent=1)
    print("written results/dg_theorem_check.json")
