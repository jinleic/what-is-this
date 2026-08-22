"""Where EXACTLY does 2D-solvability break?  Isolate the minimal obstruction.

Two independent diagnostics are applied to arbitrary layer graphs `Gamma`:

  D1  the frustration graph of the standard generator set {X_i} u {Z_iZ_j} contains an induced
      claw  <=>  Gamma has a vertex of degree >= 3   (Theorem 5, proofs/algebraic_obstruction.md);
  D2  dim <A,B>_Lie with A = sum_i X_i, B = sum_{<ij>} Z_iZ_j -- the Onsager algebra dimension.

If the two diagnostics agree on the SAME threshold, the obstruction is pinned down exactly.
The decisive test is the star K_{1,3}: four sites, three bonds, one vertex of degree 3.  It is the
smallest connected graph that is not a path or a cycle.  Compare it with the two other connected
4-vertex layer graphs that ARE paths/cycles.
"""

from __future__ import annotations

import json

from ising.clifford import pauli_x, zz, tfim_generators, frustration_graph, find_induced_claw
from ising.clifford.fast_lie import GeneratedAlgebra, P1


def analyse(name, n, bonds, max_dim=400000):
    xs = [pauli_x(n, i) for i in range(n)]
    zs = [zz(n, i, j) for i, j in bonds]
    gens, labels, _ = tfim_generators(list(range(n)), bonds, n)
    adj = frustration_graph(gens, n)
    claw = find_induced_claw(adj)
    deg = [0] * n
    for i, j in bonds:
        deg[i] += 1
        deg[j] += 1
    singles = xs + zs
    A = {g: 1 for g in xs}
    B = {g: 1 for g in zs}
    alg = GeneratedAlgebra([A, B], singles, n, p=P1)
    gA = list(range(len(xs)))
    gB = list(range(len(xs), len(singles)))
    d = alg.closure_dim([A, B], [gA, gB], max_dim=max_dim)
    claw_str = "-"
    if claw:
        c, lv = claw
        claw_str = f"{labels[c]} | " + ", ".join(str(labels[t]) for t in lv)
    return dict(name=name, n=n, bonds=len(bonds), maxdeg=max(deg), claw=bool(claw),
                claw_witness=claw_str, dim=d, support=alg.m)


CASES = [
    # --- all connected graphs on 4 vertices, ordered by max degree ---
    ("path P4        (Delta=2)", 4, [(0, 1), (1, 2), (2, 3)]),
    ("cycle C4       (Delta=2)", 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    ("star K_{1,3}   (Delta=3)", 4, [(0, 1), (0, 2), (0, 3)]),
    ("paw            (Delta=3)", 4, [(0, 1), (1, 2), (2, 0), (0, 3)]),
    ("K4             (Delta=3)", 4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
    # --- 5 vertices ---
    ("path P5        (Delta=2)", 5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    ("cycle C5       (Delta=2)", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    ("star K_{1,4}   (Delta=4)", 5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
    ("T (spider)     (Delta=3)", 5, [(0, 1), (1, 2), (2, 3), (1, 4)]),
    ("bull           (Delta=3)", 5, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4)]),
    # --- 6 vertices: chain vs one extra rung (the very first 3D-like layer) ---
    ("path P6        (Delta=2)", 6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]),
    ("cycle C6       (Delta=2)", 6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]),
    ("P6 + one rung  (Delta=3)", 6, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (1, 4)]),
    ("2x3 ladder     (Delta=3)", 6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
    # --- 7,8 vertices: growth of the exponential branch ---
    ("path P7        (Delta=2)", 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]),
    ("P7 + one rung  (Delta=3)", 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 5)]),
    ("path P8        (Delta=2)", 8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]),
    ("P8 + one rung  (Delta=3)", 8, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (1, 6)]),
]

if __name__ == "__main__":
    rows = []
    hdr = f"{'layer graph':26s} {'n':>2s} {'|E|':>3s} {'Dmax':>4s} {'claw':>5s} {'dim<A,B>':>9s} {'support':>9s}"
    print(hdr)
    print("-" * len(hdr))
    for name, n, bonds in CASES:
        r = analyse(name, n, bonds)
        rows.append(r)
        print(f"{r['name']:26s} {r['n']:2d} {r['bonds']:3d} {r['maxdeg']:4d} "
              f"{'YES' if r['claw'] else 'no':>5s} {r['dim']:9d} {r['support']:9d}", flush=True)
    print()
    print("claw witness for the star K_{1,3}:", [r for r in rows if r["name"].startswith("star K_{1,3}")][0]["claw_witness"])
    ok = all((r["maxdeg"] >= 3) == r["claw"] for r in rows)
    print(f"\nDIAGNOSTIC AGREEMENT (claw <=> Delta>=3): {'PASS' if ok else 'FAIL'}")
    poly = [r for r in rows if not r["claw"]]
    expo = [r for r in rows if r["claw"]]
    print("claw-free layers  : dim<A,B> =", {r["n"]: r["dim"] for r in poly})
    print("clawed layers     : dim<A,B> =", {r["name"].split()[0] + str(r["n"]): r["dim"] for r in expo})
    json.dump(rows, open("results/minimal_obstruction.json", "w"), indent=1)
    print("\nwritten results/minimal_obstruction.json")
