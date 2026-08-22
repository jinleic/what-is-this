# FROZEN SPEC — edge-extremal Ramsey(4,5,n) census by vertex gluing

Purpose: gate 2 of the R(5,5) campaign. Reproduce strata of the
Angeltveit–McKay ≤46 proof (arXiv:2409.15709 Sec. 3) from an independent
implementation. This spec is the single source of truth for BOTH independent
implementations (A: `src/`, B: `srcB/`). Neither implementation may read the
other's code. Frozen 2026-08-13; changes require a dated erratum section.

## Problem

For given (n, e): enumerate ALL graphs in R(4,5,n, e edges) up to
isomorphism — every graph on n vertices with **no K4, no independent 5-set,
exactly e edges**. Output may contain isomorphic duplicates; completeness
(every iso class appears at least once) and soundness (every output graph is
in the class) are required. Dedup happens downstream via canonical labeling.

## Decomposition (PROVED frame; verified on the 31 published D3 graphs)

Let G ∈ R(4,5,n). For any vertex v of degree d:
- H := G[N(v)] has no K3 (a triangle in N(v) plus v is a K4) and no I5,
  so H ∈ R(3,5,d), hence d ≤ 13 (R(3,5) = 14).
- K := G[V∖(N(v)∪{v})] has no K4 and no I4 (an I4 avoiding N(v) plus v is
  an I5), so K ∈ R(4,4, n−1−d), hence n−1−d ≤ 17 (R(4,4) = 18).
- The remaining unknown is the bipartite cone: for each k ∈ V(K) a set
  S_k ⊆ V(H) of cone-neighbors.

**Max-degree convention.** Enumerate d over Δ_range = [⌈2e/n⌉, min(13, n−1)]
and require every glued degree ≤ d. Every target graph G is produced at
d = Δ(G) (glue over a maximum-degree vertex); it is produced at no other d
(at d′ > Δ(G) the center vertex would need degree d′ > Δ; at d′ < Δ(G) the
degree cap excludes it). Δ(G) ≥ ⌈2e/n⌉ by the handshake lemma.

**Degree caps.** deg_G(h) = 1 + deg_H(h) + |{k : h ∈ S_k}| ≤ d and
deg_G(k) = deg_K(k) + |S_k| ≤ d.

**Edge count.** |cone| := Σ_k |S_k| = e − d − e(H) − e(K); skip (H,K) pairs
where this is outside [0, d·(n−1−d)].

## Constraint set (complete case analysis)

Write splits of a K4 / I5 across (v, H, K). H is K3-free and I5-free; K is
K4-free and I4-free; v is adjacent to all of H and none of K.

K4 splits (|∩H|, |∩K|), v excluded automatically except noted:
- (4,0), (3,1), any with v + ≥2 H-vertices adjacent: impossible, H K3-free.
- (0,4): impossible, K K4-free. v + 3 K-vertices: v–K non-adjacent.
- **(2,2)**: H-edge {h1,h2} × K-edge {k1,k2} with all four cross edges.
  Constraint: for every K-edge {k1,k2}: S_{k1} ∩ S_{k2} contains no H-edge
  (i.e. is an independent set of H).
- **(1,3)**: h × K-triangle {k1,k2,k3} fully crossed.
  Constraint: for every K-triangle: S_{k1} ∩ S_{k2} ∩ S_{k3} = ∅.

I5 splits (|∩H|, |∩K|); v-containing sets die (v adjacent to all H; K has
no I4); (5,0) dies (H has no I5); (0,5),(1,4) die (K has no I4):
- **(4,1)**: independent 4-set Q ⊆ H × single k, no cross edges.
  Constraint: every S_k intersects every independent 4-set of H.
- **(3,2)**: independent triple T ⊆ H × K-independent pair {k1,k2}.
  Constraint: for every independent pair of K: S_{k1} ∪ S_{k2} intersects
  every independent triple of H.
- **(2,3)**: independent pair P ⊆ H × K-independent triple {k1,k2,k3}.
  Constraint: S_{k1} ∪ S_{k2} ∪ S_{k3} intersects every independent pair
  of H.

These, plus the catalog properties of H and K, are **equivalent** to the
glued graph being in R(4,5,n) (soundness and completeness of the split is
the case analysis above). Implementations MUST also re-verify every output
with an independent full-graph check (A uses `src/check_ramsey.py`; B must
implement or reuse its own full check) — the constraint machinery is for
search, the final check is the safety net.

## Inputs (all gate-1 validated; see data/VALIDATION.json)

- `data/r35_{d}.g6`, d ≤ 13 — complete R(3,5,d) catalogs.
- `data/r44_{q}.g6`, 1 ≤ q ≤ 17 — complete R(4,4,q) catalogs.
  R(4,4,0) = the empty graph on 0 vertices (synthesize; not a file).
- Expected answers for testing: `data/r45extreme/r45{n}.{e}.g6`.

## Output contract (identical for A and B — enables exact cross-diff)

Vertex order of each output graph: index 0 = v; 1..d = H vertices in
catalog line order; d+1..n−1 = K vertices in catalog line order. Emit:
1. `out{A|B}/r45_{n}_{e}.g6` — one graph6 line per solution (raw, dup-laden).
2. `out{A|B}/r45_{n}_{e}.counts.csv` — rows `d,h_idx,k_idx,n_solutions`
   for every (H,K) pair with ≥1 solution (h_idx,k_idx = 0-based line index
   in the catalog file).
Because both implementations enumerate over identically-labeled catalog
graphs, the raw g6 line MULTISETS and the counts tables must be IDENTICAL.
Any mismatch = a bug in one side. This is a far stronger check than
comparing iso classes.

## Test ladder (published counts, gate-1 validated; run in order)

| n | e | expected iso classes |
|---|---|---|
| 12 | 48 | 1 |
| 13 | 53 | 2 |
| 13 | 52 | 10 |
| 14 | 60 | 1 |
| 15 | 66 | 1 |
| 16 | 72 | 5 |
| 16 | 71 | 138 |
| 17 | 79 | 1 |
| 17 | 78 | 86 |
| **21** | **107** | **31 ← gate 2 target (D3 stratum)** |

Iso-class comparison: `shortg -u` canonical count must equal expected, and
the canonical set must equal `labelg`-canonicalized published file. Raw
multiset cross-diff A vs B on every stratum both have run.

## Notes for implementers

- Bitmask adjacency; the I5 constraints reduce to 3 precomputed hitting
  tables over subsets of V(H) (hits-all-I2s / I3s / I4s of H) — union the
  S_k over an independent pair/triple of K and look up. (2,2) reduces to an
  is-independent-in-H table on S_{k1} ∩ S_{k2}.
- Enumerate K-vertices in a fixed order (record it if not catalog order —
  the counts table must still be per catalog labeling), assigning S_k with
  DFS + degree-cap, edge-budget, and constraint propagation. Group pair/
  triple constraints by their largest K-index so each is checked exactly
  once, when it becomes fully assigned.
- Candidate S_k pre-filter: |S_k| ≤ d − deg_K(k), S_k hits all I4(H),
  |S_k| ≥ 1 whenever H has at least one independent 4-set.
- Suggested budget pruning: maintain min/max achievable Σ|S_k| over
  unassigned k from static per-k candidate size ranges.
