"""Regression test for `plus_boundary=True` -- the path the LOW-TEMPERATURE series depends on.

The adversarial audit found that `box_broken_bond_poly(..., plus_boundary=True)` undercounted one
ghost bond per site whenever the transfer length was `c = 1` (the single layer is both the first
and the last layer and must carry its z-ghosts twice).  There was no test covering this path at
all.  This file supplies one: an INDEPENDENT brute-force implementation of the `+`-boundary
partition function, compared exactly (integer for integer) against the transfer matrix.

Independent definition used here.  For a box with FREE internal boundaries, every site `i` has
`ghost(i) = 2*dim - deg_box(i)` extra bonds to spins frozen at `+1`.  For a spin configuration
`sigma`,

    q(sigma) = #{internal bonds with sigma_i != sigma_j}  +  sum_{i : sigma_i = -1} ghost(i)

and the polynomial is `Xi(x) = sum_sigma x^{q(sigma)}`.  This is exactly the position-independent
droplet weighting the low-temperature finite-lattice method requires: `Xi(0) = 1` because the
all-up state is the unique configuration with no broken bond.
"""

from __future__ import annotations

import itertools

import numpy as np

from ising.lattices import hyperrect
from ising.transfer_matrix import box_broken_bond_poly

FAILS = []


def brute_plus_boundary(shape):
    """Independent brute force: enumerate all 2^N states, count internal + ghost broken bonds."""
    dim = len(shape)
    lat = hyperrect(shape, periodic=False)
    N = lat.n_sites
    bonds = lat.bonds
    deg = [0] * N
    for i, j in bonds:
        deg[i] += 1
        deg[j] += 1
    ghost = [2 * dim - d for d in deg]
    assert all(g >= 0 for g in ghost)
    maxq = len(bonds) + sum(ghost)
    out = [0] * (maxq + 1)
    for s in range(1 << N):
        q = 0
        for i, j in bonds:
            if ((s >> i) ^ (s >> j)) & 1:
                q += 1
        for i in range(N):
            if not ((s >> i) & 1):          # bit 0 == spin down
                q += ghost[i]
        out[q] += 1
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def check(name, a, b):
    ok = (a == b)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        n = max(len(a), len(b))
        aa = a + [0] * (n - len(a))
        bb = b + [0] * (n - len(b))
        first = next((k for k in range(n) if aa[k] != bb[k]), None)
        print(f"        first difference at x^{first}: transfer {aa[first]} vs brute {bb[first]}")
        print(f"        transfer {a[:12]}")
        print(f"        brute    {b[:12]}")
        FAILS.append(name)


if __name__ == "__main__":
    print("plus_boundary transfer matrix vs independent brute force (exact integers)")
    print("\n3D boxes (ghost count = 6 - deg):")
    for shape in [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1), (2, 1, 2), (1, 2, 2),
                  (2, 2, 2), (3, 1, 1), (3, 2, 1), (2, 2, 3), (3, 3, 1), (3, 2, 2), (3, 3, 2)]:
        tm = box_broken_bond_poly(shape, plus_boundary=True)
        bf = brute_plus_boundary(shape)
        check(f"3D {shape[0]}x{shape[1]}x{shape[2]}", tm, bf)

    print("\n2D boxes (ghost count = 4 - deg):")
    for shape in [(1, 1), (2, 1), (1, 2), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3)]:
        tm = box_broken_bond_poly(shape, plus_boundary=True)
        bf = brute_plus_boundary(shape)
        check(f"2D {shape[0]}x{shape[1]}", tm, bf)

    print("\nstructural invariants:")
    for shape in [(2, 2, 2), (3, 3, 2), (3, 2, 1)]:
        tm = box_broken_bond_poly(shape, plus_boundary=True)
        N = shape[0] * shape[1] * shape[2]
        ok1 = tm[0] == 1
        ok2 = sum(tm) == 1 << N
        print(f"  [{'PASS' if ok1 else 'FAIL'}] {shape}: Xi(0)=1 (unique all-up ground state)")
        print(f"  [{'PASS' if ok2 else 'FAIL'}] {shape}: sum of coefficients = 2^{N}")
        if not ok1:
            FAILS.append(f"Xi(0) {shape}")
        if not ok2:
            FAILS.append(f"sum {shape}")

    print()
    if FAILS:
        print(f"FAIL: {len(FAILS)} checks failed: {FAILS}")
        raise SystemExit(1)
    print("PASS: plus_boundary agrees with independent brute force on all boxes, including c=1")
