# Adversarial audit of `census_gluing_spec.md` (frozen 2026-08-13)

Audit date: 2026-08-13. Auditor: Claude (independent subagent; did not read
`src/` or `srcB/`; all checkers written from scratch under `scratch_audit/`).
Mandate: break the spec. Verdict summary at the end.

## 1. Independent re-derivation of the K4/I5 split case analysis

Setup: G on n vertices, v of degree d, H = G[N(v)], K = G[V∖(N(v)∪{v})],
S_k = cone neighbourhood of k ∈ V(K) in V(H). v is adjacent to all of H,
none of K.

### K4 splits (from scratch)

A K4 is a pairwise-adjacent 4-set. Partition its vertices as
(#{v}, #H, #K):

| split | verdict | reason |
|---|---|---|
| (1,3,0) | impossible | 3 H-vertices pairwise adjacent = triangle in H; H is K3-free |
| (1,2,1), (1,1,2), (1,0,3) | impossible | contains a v–K pair; v is non-adjacent to K |
| (0,4,0) | impossible | K4 in H; H K3-free ⇒ K4-free |
| (0,3,1) | impossible | triangle in H |
| (0,2,2) | **constraint** | H-edge {h1,h2} ⊆ S_{k1} ∩ S_{k2} for a K-edge {k1,k2} |
| (0,1,3) | **constraint** | h ∈ S_{k1} ∩ S_{k2} ∩ S_{k3} for a K-triangle |
| (0,0,4) | impossible | K is K4-free |

Matches the spec exactly: the only surviving splits are (2,2) and (1,3),
with exactly the constraints the spec states. The (2,2) constraint
"S_{k1} ∩ S_{k2} contains no H-edge, i.e. is independent in H" is a correct
iff-rewriting; the (1,3) constraint "triple intersection empty over
K-triangles" likewise. Note (1,3) is NOT implied by (2,2) applied to the
triangle's three edges (a single common vertex violates nothing pairwise),
so both are needed — the spec has both.

Clarity note (no math error): the spec line
"(4,0), (3,1), any with v + ≥2 H-vertices adjacent: impossible, H K3-free"
attributes the v+2H+1K and v+1H+2K cases to K3-freeness, but those actually
die because the K-vertex would need a v–K edge. Every v-containing case does
die (any K4 vertex adjacent to v lies in H by definition of H, so v-cases
reduce to v+3H = triangle in H); the conclusion is right, the justification
is compressed/scattered.

### I5 splits (from scratch)

An I5 is a pairwise-non-adjacent 5-set. v-containing: the other 4 are
non-neighbours of v, hence all in K, forming an I4 — impossible (K I4-free).
Without v, (|∩H|, |∩K|):

| split | verdict | reason |
|---|---|---|
| (5,0) | impossible | H is I5-free |
| (4,1) | **constraint** | I4 Q ⊆ H with S_k ∩ Q = ∅ |
| (3,2) | **constraint** | I3 T ⊆ H, K-nonadjacent pair, (S_{k1} ∪ S_{k2}) ∩ T = ∅ |
| (2,3) | **constraint** | I2 P ⊆ H, K-independent triple, union ∩ P = ∅ |
| (1,4) | impossible | I4 in K |
| (0,5) | impossible | contains I4 in K |

Matches the spec exactly. No missed splits; no split dismissed incorrectly.
The "equivalent" claim (constraints + H ∈ R(3,5,d) + K ∈ R(4,4,q) ⟺ glued
graph ∈ R(4,5,n)) follows since the splits are exhaustive; verified
empirically in both directions (Section 5, check D).

## 2. Max-degree convention and degree-cap arithmetic

- deg_G(h) = 1 (edge to v) + deg_H(h) + #{k : h ∈ S_k}: **correct**, and
  verified as an exact identity on every published decomposition (not just
  ≤ d). deg_G(k) = deg_K(k) + |S_k|: same. No off-by-one.
- Production at exactly d = Δ(G): at level d the centre has degree exactly d
  and the cap forces all degrees ≤ d, so every output at level d has
  Δ = d exactly; conversely G is reachable at d = Δ(G) via any max-degree
  vertex (its H and K are in the complete catalogs — verified empirically,
  Section 5 check C). Sound in both directions.
- Lower endpoint ⌈2e/n⌉: Δ ≥ 2e/n (max ≥ mean) and Δ ∈ ℤ give
  Δ ≥ ⌈2e/n⌉. **The bound is tight in the ladder**: the unique r45(12,48)
  graph is 8-regular with ⌈96/12⌉ = 8, so a stricter bound (e.g.
  ⌊2e/n⌋ + 1 misapplied to the regular case, or 2e/n + 1) would MISS a real
  graph. The spec's endpoint is exactly right.
- Upper endpoint min(13, n−1): d ≤ 13 since H ∈ R(3,5,d) and R(3,5)=14;
  d ≤ n−1 trivially. Correct. If ⌈2e/n⌉ > min(13, n−1) the range is empty
  and the stratum is genuinely empty (Δ would have to exceed 13) — correct
  behaviour.
- d = n−1 ⇒ q = 0: handled by the synthesized R(4,4,0) = empty graph. This
  case is actually reachable in the ladder (n=12: range [8,11], 11 = n−1),
  so the "synthesize" instruction is load-bearing — good that it is there.
- **Gap (minor, robustness):** the range does not clamp d ≥ n−18. For
  d with q = n−1−d ≥ 18, R(4,4,q) = ∅ but there is no catalog file
  (`data/r44_{q}.g6` exists only for q ≤ 17); a literal implementation would
  try to open a nonexistent file. Harmless for the entire test ladder
  (worst case: n=21, e=107 gives q ∈ [7,9]; verified all ladder strata have
  q ≤ 9), but the spec text should say d ≥ max(⌈2e/n⌉, n−18) or "treat
  q > 17 as an empty catalog" if it is ever run at larger n / smaller e.
- Catalog-order-d-has-no-graphs concern: r35_d is nonempty for every
  d ≤ 13 (counts 1,2,3,7,13,32,71,179,290,313,105,12,1) and r44_q nonempty
  for every q ≤ 17, so the concern is vacuous for the given inputs; a d
  whose (H,K) pairs all fail the edge budget simply contributes nothing,
  which is correct (e.g. n=12, e=48 at d=11 forces |cone| < 0 ⇒ skipped).

## 3. Constraint-to-table reductions (implementer notes)

- "S_k hits all I4(H)" ⟺ no (4,1) violation for that k: exact.
- Pair-union hits-all-I3(H) over K-nonadjacent pairs ⟺ no (3,2): exact.
- Triple-union hits-all-I2(H) over K-independent triples ⟺ no (2,3): exact.
- Intersection-is-independent-in-H over K-edges ⟺ no (2,2): exact.
- **Gap (minor):** the notes' reduction list covers the three I5 tables and
  (2,2) but never mentions how (1,3) (K-triangle triple-intersection-empty)
  is to be checked. The normative "Constraint set" section does include it,
  and "group pair/triple constraints by largest K-index" implicitly covers
  K-triangles, but an implementer treating the notes' table list as the
  complete checking machinery would omit (1,3) and emit non-Ramsey graphs
  (caught by the mandated final full-graph check and by the A/B cross-diff,
  so wasted work, not a wrong census). Recommend an erratum sentence.
- Pre-filter "|S_k| ≥ 1 whenever H has at least one independent 4-set" is
  correctly conditional (if I4(H) = ∅, empty S_k is legitimately allowed).

## 4. Expected-count table vs data

Line counts of `data/r45extreme/r45{n}.{e}.g6`:
12.48→1, 13.53→2, 13.52→10, 14.60→1, 15.66→1, 16.72→5, 16.71→138,
17.79→1, 17.78→86, 21.107→31. **All ten match the spec's table exactly.**
(These files are iso-deduped, consistent with using them as the `labelg`
reference set.)

## 5. Empirical checks (all scripts in `scratch_audit/`, own g6 parser,
nothing imported from src/)

A. `audit_necessity.py` — all 276 graphs of all ten ladder strata:
   - published graph really is K4-free, I5-free, with the stated e;
   - at EVERY vertex v (frame claim "for any vertex"): H K3-free and
     I5-free, K K4-free and I4-free, d ≤ 13, n−1−d ≤ 17;
   - at every max-degree vertex: all five cone constraints hold on the real
     cone; degree-cap identities exact; |cone| = e − d − e(H) − e(K) exact
     and within [0, d(n−1−d)]; Δ within the spec range.
   **Result: 0 failures** (276 graphs, every vertex).

B. `audit_drange.py` — per-stratum Δ_range endpoints, implied q range
   (never > 17 or < 0 in the ladder), observed Δ multisets all in range,
   all input catalogs nonempty. **Pass.** Observed Δ distributions, e.g.
   21.107: {11: 30, 12: 1}; 16.71: {9:1, 10:106, 11:26, 12:5}; 12.48: {8:1}.

C. `audit_catalog_membership.py` — own invariant-refined backtracking
   isomorphism test: for every published graph and EVERY max-degree vertex,
   H is isomorphic to a line of `r35_d.g6` and K to a line of `r44_q.g6`
   (q=0 synthesized). 1243 (v,H,K) decompositions across strata
   12.48/13.53/13.52/16.71/17.78/21.107. **0 failures** — the search space
   provably contains a producing (d, h_idx, k_idx, cone) tuple for every
   published graph.

D. `audit_equiv3.py` — the "equivalent" claim tested as a biconditional:
   - Part A: EXHAUSTIVE enumeration of all 2^(d·q) cone assignments for
     every catalog (H,K) pair with d·q ≤ 12 (d ≤ 4), comparing the
     constraint-set verdict against a direct has-K4/has-I5 check of the
     glued graph.
   - Part B: 600 boundary mutants per published graph (flip 1–3 cone bits
     of the real decomposition) for strata 13.52, 16.71, 17.78, 21.107,
     with an incremental exact Ramsey check (a new K4/I5 must contain a
     flipped pair).
   **Result: 240,906 exhaustive assignments (184,245 satisfying both sides)
   + 159,000 mutants — 0 mismatches in either direction.**

## 6. Findings

| # | severity | finding |
|---|---|---|
| 1 | minor | Implementer-notes reduction list omits the (1,3) K-triangle check; normative section has it, but notes read as if the 4 listed tables were the whole machinery. Erratum suggested. |
| 2 | minor | Δ_range lower end not clamped to d ≥ n−18; q = n−1−d > 17 would reference a nonexistent catalog file. Unreachable in the entire test ladder (max q = 9), so robustness only. |
| 3 | note | K4-with-v impossibility justification compressed: v+2H+1K / v+1H+2K die by v–K non-adjacency, not "H K3-free" as the grouped line suggests. Conclusion unaffected. |
| 4 | note | Degenerate inputs (e = 0 ⇒ d = 0, H on 0 vertices) unspecified; irrelevant to all ladder strata. |

No fatal or major finding. Every claimed-impossible split is impossible,
every constraint is exactly the surviving split it claims to encode, the
degree/edge arithmetic is exact on all real data, the Δ-range is tight at
both ends on real data, and the constraint set is empirically a true
biconditional on ~241k exhaustive assignments plus 159k boundary mutants.
The spec as frozen will produce a correct census on the ladder strata.
