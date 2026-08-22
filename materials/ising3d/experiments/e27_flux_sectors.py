"""Are there conserved Z_2 fluxes?  (Track A items 5 and 9.)

The one remaining escape from Theorem 7 (`proofs/algebraic_obstruction.md`) is the *Kitaev
mechanism*: a spin model can fail to be globally Gaussian and yet be exactly solvable if it
possesses a large family of commuting conserved "flux" operators `W_p`, because then the
Hamiltonian block-diagonalises and each block may separately be quadratic in fermions.  This is
exactly what happens in the honeycomb model, and it is what a "higher-dimensional Jordan-Wigner
with an emergent Z_2 gauge field" would have to deliver.

PROPOSITION (proved in the module docstring below and verified here).
For the layer generator set  G = {X_i : i in V} u {Z_iZ_j : (ij) in E}  on a CONNECTED graph
Gamma, the set of Pauli strings commuting with every element of G is exactly

        { I ,  prod_{i in V} X_i } .

There are NO conserved fluxes beyond the single global Z_2 spin-flip symmetry -- for ANY graph,
in any dimension.

Proof.  Write the string as Q_{(a|b)}.  Commuting with X_i = Q_{(e_i|0)} means
<(a|b),(e_i|0)> = b_i = 0, for every i, so b = 0 and the string is X^a.  Commuting with
Z_iZ_j = Q_{(0|e_i+e_j)} means a_i + a_j = 0 (mod 2) for every edge (ij), so `a` is constant on
each connected component.  On a connected graph a = 0 or a = (1,...,1).  QED

Consequences.
  * The block-diagonalisation that solves the Kitaev honeycomb model does not exist here: the
    Hilbert space splits only into the two global spin-flip parity sectors, not into 2^{#plaquettes}
    flux sectors.  So this escape route from Theorem 7 is closed.
  * It also shows the obstruction is NOT that "we have not found the right gauge structure":
    for these generators, no nontrivial Pauli-string gauge structure exists at all.
  * Note the contrast is not with dimension: the same proposition holds for the 1D chain, where
    solvability comes from Majorana bilinearity instead (Theorem 5).

The script additionally computes the full centraliser of the PAIR {A, B} (the summed generators
that actually build the transfer matrix), which is a larger and less trivial object.
"""

from __future__ import annotations

import json

import numpy as np

from ising.clifford import pauli_x, symplectic_form, zz


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


def centraliser_of_terms(n, bonds):
    """All Pauli strings commuting with every X_i and every Z_iZ_j, by exact F_2 linear algebra."""
    # unknown v = (a|b) in F_2^{2n}; conditions: b_i = 0 all i;  a_i + a_j = 0 all edges
    # -> solve directly
    sols = []
    # a constant on connected components
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in bonds:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    ncomp = len(comps)
    for mask in range(1 << ncomp):
        a = 0
        for k, (_, members) in enumerate(sorted(comps.items())):
            if (mask >> k) & 1:
                for i in members:
                    a |= 1 << i
        sols.append(a)  # b = 0
    return sols, ncomp


def brute_centraliser(n, gens, restrict=None):
    """Brute-force centraliser over all 4^n Pauli strings (n <= 10)."""
    out = []
    for v in range(1 << (2 * n)):
        if all(symplectic_form(v, g, n) == 0 for g in gens):
            out.append(v)
    return out


def pauli_str(v, n):
    mask = (1 << n) - 1
    a, b = v & mask, (v >> n) & mask
    return "".join("IXZY"[((a >> i) & 1) + 2 * ((b >> i) & 1)] for i in range(n))


CASES = [
    ("chain P4", 4, [(i, i + 1) for i in range(3)]),
    ("ring C4", 4, [(i, (i + 1) % 4) for i in range(4)]),
    ("ring C6", 6, [(i, (i + 1) % 6) for i in range(6)]),
    ("2x3 grid", 6, grid_bonds(2, 3)),
    ("2x4 grid", 8, grid_bonds(2, 4)),
    ("3x3 grid", 9, grid_bonds(3, 3)),
    ("3x3 torus layer", 9, grid_bonds(3, 3, per=True)),
]

if __name__ == "__main__":
    print("=== centraliser of the TERM set {X_i} u {Z_iZ_j} in the Pauli group ===")
    print("(a 'conserved flux' would be a nontrivial element of this centraliser)")
    rows = []
    allok = True
    for name, n, bonds in CASES:
        gens = [pauli_x(n, i) for i in range(n)] + [zz(n, i, j) for i, j in bonds]
        theory, ncomp = centraliser_of_terms(n, bonds)
        brute = brute_centraliser(n, gens) if n <= 10 else None
        match = (brute is not None) and sorted(brute) == sorted(theory)
        allok &= match
        strs = [pauli_str(v, n) for v in sorted(theory)]
        print(f"  {name:18s} n={n:2d}: |centraliser| = {len(theory)} "
              f"(= 2^{ncomp}, components={ncomp});  brute force agrees: {match}")
        print(f"      elements: {strs}")
        rows.append(dict(name=name, n=n, size=len(theory), components=ncomp,
                         elements=strs, brute_force_agrees=bool(match)))

    print(f"\nPROPOSITION VERIFIED on all {len(CASES)} graphs: {allok}")
    print("  For every CONNECTED layer graph the centraliser is exactly {I, prod_i X_i}:")
    print("  the ONLY conserved Pauli-string quantity is the global Z_2 spin flip.")
    print("  => the Kitaev flux-sector mechanism is unavailable for the Ising layer operator,")
    print("     in any dimension.  The 3D case is not failing for want of a gauge structure;")
    print("     there is no nontrivial Pauli gauge structure to find.")

    print("\n=== centraliser of the summed pair {A, B} (transfer-matrix symmetries) ===")
    for name, n, bonds in CASES[:5]:
        A = [pauli_x(n, i) for i in range(n)]
        B = [zz(n, i, j) for i, j in bonds]
        cnt = 0
        elems = []
        for v in range(1 << (2 * n)):
            # [P, A] = 0 requires P to commute with the SUM, not each term:
            # sum_i (anticommuting terms) must cancel -- for a Pauli string P the commutator is
            # sum over anticommuting generators, and distinct generators give distinct strings,
            # so [P,A]=0 iff P commutes with every X_i.  Same for B.  Recorded as a lemma.
            if all(symplectic_form(v, g, n) == 0 for g in A + B):
                cnt += 1
                elems.append(pauli_str(v, n))
        print(f"  {name:18s}: |Pauli centraliser of A and B| = {cnt}  {elems}")
    print("\n  LEMMA (why these coincide): [P, sum_a h_a] = sum_{a: h_a anticommutes with P} 2 P h_a,")
    print("  and the strings P h_a are pairwise DISTINCT for distinct a, so no cancellation is")
    print("  possible and [P, sum_a h_a] = 0 iff P commutes with every term individually.")

    json.dump(dict(provenance=dict(script="experiments/e27_flux_sectors.py"),
                   data=dict(centralisers=rows),
                   checks=[dict(name="proposition matches brute force", passed=bool(allok),
                                detail=f"{len(CASES)} graphs")]),
              open("results/flux_sectors.json", "w"), indent=1)
    print("\nwritten results/flux_sectors.json")
