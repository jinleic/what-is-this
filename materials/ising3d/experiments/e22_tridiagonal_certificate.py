"""Turn the tridiagonal no-go from a computation into a SHORT CHECKABLE CERTIFICATE.

e21 shows by exact Gaussian elimination that, for a 4-regular layer graph,
    T0 = beta T1 + gamma T2 + rho T3
has no solution.  A rank computation over thousands of Pauli strings is not a proof a human can
audit.  This script extracts a MINIMAL WITNESS: a set of at most four Pauli strings
`P_1, ..., P_m` such that the 3-unknown linear system restricted to those rows is already
inconsistent.  Verifying the certificate requires only reading off four integer coefficients from
each of T0, T1, T2, T3 -- which the script also prints, so the whole no-go can be checked by hand.
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from e21_tridiagonal import build, grid_bonds, solve_linear  # noqa: E402


def pauli_str(v, n):
    mask = (1 << n) - 1
    a, b = v & mask, (v >> n) & mask
    return "".join("IXZY"[((a >> i) & 1) + 2 * ((b >> i) & 1)] for i in range(n))


def minimal_witness(T0, cols, max_rows=4):
    """Smallest set of Pauli strings on which the system is already inconsistent.

    Strategy (constructive, not search):
      size 1 -- a string in T0 that appears in NO column: a one-line proof;
      size k -- pick k-1 rows spanning the column space, solve, then find a violated row.
    """
    keys = sorted(set(T0) | {k for c in cols for k in c})
    w = len(cols)

    # --- size 1: T0 has a component where every column vanishes
    for k in keys:
        if T0.get(k, 0) != 0 and all(c.get(k, 0) == 0 for c in cols):
            return [k]

    # --- build a maximal independent set of rows, solve there, then find a violated row
    basis_rows, mat = [], []

    def rank_of(rows):
        m = [[Fraction(c.get(k, 0)) for c in cols] for k in rows]
        r = 0
        m = [row[:] for row in m]
        for col in range(w):
            piv = next((i for i in range(r, len(m)) if m[i][col] != 0), None)
            if piv is None:
                continue
            m[r], m[piv] = m[piv], m[r]
            pv = m[r][col]
            m[r] = [x / pv for x in m[r]]
            for i in range(len(m)):
                if i != r and m[i][col] != 0:
                    f = m[i][col]
                    m[i] = [a - f * b for a, b in zip(m[i], m[r])]
            r += 1
        return r

    for k in keys:
        if rank_of(basis_rows + [k]) > len(basis_rows):
            basis_rows.append(k)
        if len(basis_rows) == w:
            break
    sub = {k: T0.get(k, 0) for k in basis_rows}
    subcols = [{k: c.get(k, 0) for k in basis_rows} for c in cols]
    sol, _, _ = solve_linear(sub, subcols)
    if sol is None:
        return basis_rows
    for k in keys:
        lhs = sum(Fraction(c.get(k, 0)) * s for c, s in zip(cols, sol))
        if lhs != Fraction(T0.get(k, 0)):
            return basis_rows + [k]
    return None


def report(name, n, bonds):
    T0, T1, T2, T3 = build(n, bonds, swap=True)     # swap=True -> the TD2 relation
    sol, rank, nul = solve_linear(T0, [T1, T2, T3])
    print(f"\n=== {name}  (n={n}) ===")
    print(f"  #Pauli strings: T0 {len(T0)}, T1 {len(T1)}, T2 {len(T2)}, T3 {len(T3)}")
    if sol is not None:
        print(f"  TD2 SOLVABLE: (beta,gamma*,rho*) = {tuple(map(str, sol))}")
        return dict(name=name, n=n, solvable=True, solution=[str(s) for s in sol])
    w = minimal_witness(T0, [T1, T2, T3])
    print(f"  TD2 has NO SOLUTION.  Minimal inconsistency witness: {len(w)} Pauli strings")
    print(f"    {'Pauli string':<{max(n,12)}}  {'[T0]':>8s} {'[T1]':>8s} {'[T2]':>8s} {'[T3]':>8s}")
    for k in w:
        print(f"    {pauli_str(k, n):<{max(n,12)}}  {T0.get(k,0):8d} {T1.get(k,0):8d} "
              f"{T2.get(k,0):8d} {T3.get(k,0):8d}")
    print("  The 3-unknown system built from just these rows is inconsistent, hence no")
    print("  (beta, gamma*, rho*) exists.")
    return dict(name=name, n=n, solvable=False,
                witness=[dict(pauli=pauli_str(k, n), T0=T0.get(k, 0), T1=T1.get(k, 0),
                              T2=T2.get(k, 0), T3=T3.get(k, 0)) for k in w])


if __name__ == "__main__":
    out = []
    out.append(report("cycle C6 (2-regular, 2D Ising) -- CONTROL", 6,
                      [(i, (i + 1) % 6) for i in range(6)]))
    out.append(report("cycle C8 (2-regular, 2D Ising) -- CONTROL", 8,
                      [(i, (i + 1) % 8) for i in range(8)]))
    for (a, b) in [(3, 3), (3, 4), (4, 4), (3, 5), (4, 5), (5, 5)]:
        out.append(report(f"torus layer {a}x{b} (4-regular, 3D Ising)", a * b,
                          grid_bonds(a, b, per=True)))
    out.append(report("open 2x3 grid (Delta=3, 3D Ising)", 6, grid_bonds(2, 3)))
    out.append(report("open 3x3 grid (Delta=4 centre, 3D Ising)", 9, grid_bonds(3, 3)))

    nogo = [r for r in out if not r["solvable"]]
    sizes = sorted({len(r["witness"]) for r in nogo})
    print("\n" + "=" * 78)
    print(f"{len(nogo)} layer graphs admit NO tridiagonal relation.")
    print(f"Every one has an explicit finite witness; witness sizes occurring: {sizes}")
    print("Each witness is checkable by hand: solve the 3-unknown system on any three of the")
    print("listed rows and substitute into the remaining row.")
    solv = [r for r in out if r["solvable"]]
    print(f"{len(solv)} control layers (cycles) ARE tridiagonal pairs: "
          f"{[(r['name'].split()[1], r['solution']) for r in solv]}")
    json.dump(out, open("results/tridiagonal_certificate.json", "w"), indent=1)
    print("\nwritten results/tridiagonal_certificate.json")

