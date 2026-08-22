"""Exact LOCAL structure of the Dolan-Grady defect.

e14 established: for A = sum_v X_v, B = sum_{(ij) in E} Z_iZ_j on a layer graph Gamma,

    R(Gamma) := [B,[B,[B,A]]] - 16 [B,A]

vanishes identically for cycles and is nonzero otherwise.  Here we determine R exactly:

  (H1)  R is a SUM OF LOCAL TERMS, one per vertex: R = sum_v R_v with R_v supported on the
        closed neighbourhood of v.  Test: R is additive over disjoint unions, and R_v depends on
        v only through deg(v) and the identity of its neighbours.
  (H2)  R_v = 0  <=>  deg(v) = 2.
  (H3)  For deg(v) = d >= 3 the residual contains the (d+1)-BODY term  Y_v prod_{u~v} Z_u ,
        which is the minimal operator that cannot be written as a Majorana bilinear.

If (H1)-(H3) hold, the obstruction is *extensive*: for a 1D layer only the two endpoints
contribute (a surface term, O(1)), whereas for a 2D layer a finite FRACTION of all vertices
contributes (a bulk term, Theta(n)).  That is the precise, quantitative sense in which the
2D Ising model is solvable and the 3D Ising model is not, within the Onsager mechanism.
"""

from __future__ import annotations

import json

from ising.clifford import pauli_x, zz
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e14_dolan_grady import comm, lin, pauli_str  # noqa: E402


def residual(n, bonds, const=16):
    xs = [pauli_x(n, i) for i in range(n)]
    A = {g: 1 for g in xs}
    B = {}
    for i, j in bonds:
        B[zz(n, i, j)] = B.get(zz(n, i, j), 0) + 1
    BA = comm(B, A, n)
    r2 = lin(comm(B, comm(B, BA, n), n), BA, 1, -const)
    return r2


def degrees(n, bonds):
    d = [0] * n
    for i, j in bonds:
        d[i] += 1
        d[j] += 1
    return d


def support_of(v, n):
    mask = (1 << n) - 1
    a, b = v & mask, (v >> n) & mask
    return frozenset(i for i in range(n) if ((a >> i) & 1) or ((b >> i) & 1))


def main():
    out = {}
    print("=== (H2) stars K_{1,d}: residual as a function of the single high-degree vertex ===")
    for d in range(1, 7):
        n = d + 1
        bonds = [(0, k) for k in range(1, n)]
        r = residual(n, bonds)
        top = {v: c for v, c in r.items() if len(support_of(v, n)) == n}
        print(f"  K_1,{d}: n={n:2d}  #terms={len(r):3d}   max body size="
              f"{max(len(support_of(v,n)) for v in r) if r else 0}   "
              f"top-body terms: {[(c, pauli_str(v,n)) for v,c in sorted(top.items())]}")
        out[f"star_{d}"] = dict(n=n, terms=len(r),
                                maxbody=max((len(support_of(v, n)) for v in r), default=0))

    print("\n=== (H1) additivity over disjoint unions ===")
    # two disjoint triangles vs one triangle
    r1 = residual(3, [(0, 1), (1, 2), (2, 0)])
    r2 = residual(6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)])
    print(f"  triangle: {len(r1)} terms; two disjoint triangles: {len(r2)} terms -> "
          f"{'ADDITIVE (both zero)' if len(r1) == 0 and len(r2) == 0 else 'check'}")
    rp = residual(4, [(0, 1), (1, 2), (2, 3)])
    rpp = residual(8, [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7)])
    print(f"  P4: {len(rp)} terms; P4 + P4 (disjoint): {len(rpp)} terms -> "
          f"{'ADDITIVE' if len(rpp) == 2 * len(rp) else 'NOT additive'}")

    print("\n=== (H2)/(H3) residual term count vs degree sequence ===")
    cases = [
        ("ring C3", 3, [(0, 1), (1, 2), (2, 0)]),
        ("ring C6", 6, [(i, (i + 1) % 6) for i in range(6)]),
        ("path P6", 6, [(i, i + 1) for i in range(5)]),
        ("path P9", 9, [(i, i + 1) for i in range(8)]),
        ("2x3 grid", 6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
        ("2x4 grid", 8, [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7), (0, 4), (1, 5), (2, 6), (3, 7)]),
        ("3x3 grid", 9, [(0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8),
                         (0, 3), (3, 6), (1, 4), (4, 7), (2, 5), (5, 8)]),
        ("3x3 torus layer", 9, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (6, 7), (7, 8), (8, 6),
                                (0, 3), (3, 6), (6, 0), (1, 4), (4, 7), (7, 1), (2, 5), (5, 8), (8, 2)]),
    ]
    print(f"  {'graph':18s} {'degseq':>26s} {'#deg!=2':>8s} {'#terms':>7s} {'maxbody':>8s}")
    rows = []
    for name, n, bonds in cases:
        r = residual(n, bonds)
        deg = degrees(n, bonds)
        nbad = sum(1 for d in deg if d != 2)
        mb = max((len(support_of(v, n)) for v in r), default=0)
        cnt = {}
        for d in deg:
            cnt[d] = cnt.get(d, 0) + 1
        print(f"  {name:18s} {str(sorted(cnt.items())):>26s} {nbad:8d} {len(r):7d} {mb:8d}")
        rows.append(dict(name=name, n=n, degcount=cnt, nbad=nbad, terms=len(r), maxbody=mb))
    out["degree_table"] = rows

    print("\n=== per-degree residual weight (fit terms = sum_v f(deg v)) ===")
    # solve for f(d) from the table by least squares over integers
    import numpy as np
    degs = sorted({d for r in rows for d in r["degcount"]})
    Amat = np.array([[r["degcount"].get(d, 0) for d in degs] for r in rows], dtype=float)
    y = np.array([r["terms"] for r in rows], dtype=float)
    sol, res, rank, sv = np.linalg.lstsq(Amat, y, rcond=None)
    pred = Amat @ sol
    print(f"  degrees present: {degs}")
    print(f"  fitted f(d)    : {[round(s, 6) for s in sol]}")
    print(f"  exact fit?     : {'YES' if np.allclose(pred, y, atol=1e-9) else 'NO'}  (residual {np.abs(pred-y).max():.2e})")
    out["f_of_degree"] = {int(d): float(s) for d, s in zip(degs, sol)}
    out["fit_exact"] = bool(np.allclose(pred, y, atol=1e-9))

    print("\n=== (H3) explicit minimal 4-body obstruction at a degree-3 vertex ===")
    n, bonds = 7, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (3, 0)]  # C4 with a tail: has deg-3
    n, bonds = 5, [(0, 1), (1, 2), (1, 3), (3, 4)]
    r = residual(n, bonds)
    print(f"  graph: bonds {bonds}, degrees {degrees(n,bonds)}")
    for v, c in sorted(r.items(), key=lambda kv: -len(support_of(kv[0], n))):
        print(f"    {c:+6d} * {pauli_str(v, n)}   (body {len(support_of(v,n))})")

    json.dump(out, open("results/dg_local_structure.json", "w"), indent=1, default=str)
    print("\nwritten results/dg_local_structure.json")


if __name__ == "__main__":
    main()
