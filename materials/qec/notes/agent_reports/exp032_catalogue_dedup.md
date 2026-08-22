POSITIVE — 368 catalogue rows form 368 classes under G; the 14 headline rows form 14 classes under G.
Scope: these are **distinct under group G**, not proved inequivalent under any broader equivalence relation.

## Result

- Full catalogue: **368 -> 368 classes under G**; 0 multi-member classes, largest class 1.
- Stored `(n,k,d)=(144,12,12)` family: **14 -> 14 classes under G**; mutually distinct under G = `true`.
- Exact GF(2) rank audit: **368/368** rows agree with stored `n,k`; parameter corrections found: **0**.
- Generator preservation checks: **12/12 PASS** on three small instances at lattices `((6, 3), (3, 6), (6, 6))`.  Each comparison built both PBB matrices and checked exact GF(2) row-space equality after the stated qubit permutation.
- Multi-member sanity: 0 classes required the conditional check.  Because there are no multi-member classes, no `n <= 72` exact-distance job was applicable and no larger class was skipped.
- Distance caveat: selection of the headline family uses the catalogue's stored `d=12`.  Its metadata marks 7/14 exact and 7/14 not exact; this experiment does not upgrade the latter to exact(OPTIMAL).

## Group G used (and only this group)

1. Common lattice translation `tau_(i,j)`: shift every term of `A,B,C,D` by `(i,j)` modulo `(ell,m)`.  It is an invertible row relabelling (`T` on the mixed block and `T^T` on the pure-Z block).
2. Variable inversion `iota`: `(A,B,C,D) -> (A*,B*,C*,D*)`, where `F*(x,y)=F(x^-1,y^-1)=F^T`.  The exhibited qubit permutation is `(block,g)->(block,-g)`.
3. Transpose/swap `sigma`: `(A,B,C,D) -> (B*,A*,D*,C*)`.  The exhibited qubit permutation is `(L,g)<->(R,-g)`.  The `D*,C*` action is required: it sends the mixed row to `[B^T A^T | D^T C^T]`, while the pure-Z row becomes `[0 0 | A B]`.

The closure contains the plain block swap as `iota sigma`; it was not introduced as an additional assumption.  Excluded: axis exchange, general torus automorphisms, arbitrary qubit permutations, local Clifford transformations, and any other code-equivalence operation.  Proof basis: `proofs/equations.md` E4, E5(R1-R2), E9 and `proofs/pbb_structure.md` R1-R2, Proposition 4.

## Canonicalization and audits

Each polynomial is reduced as a GF(2) exponent set modulo `(ell,m)`.  The script enumerates the discrete closure of `iota,sigma` and all `ell*m` common translations, then takes the lexicographically least `(A,B,C,D)` tuple.  Rows are tabulated by stored `(n,k,d)` and canonical form; a global canonical table also detects cross-parameter collisions.  No parameter-table correction or term-normalization issue was found.

Distance checking was pre-registered only for multi-member classes with `n <= 72`, using `exact_distance_symplectic`, 120 s per logical sector, 8 workers, and exact results labelled only as `exact(OPTIMAL)`.  No class met the multi-member premise.

## Reproduce

```bash
cd /Users/jinleic/jinleic-workspace/qec-codesign
PYTHONPATH=src .venv/bin/python experiments/exp032_catalogue_dedup.py
```

Observed wall time: **1.318 s** on the shared workstation.  The classification is exhaustive for the 368 input rows and finite group G, not for broader code equivalence.
