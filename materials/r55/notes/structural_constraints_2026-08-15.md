# Structural constraints on hypothetical R(5,5,n) graphs

Session record, 2026-08-15. All computational claims below were executed and
verified in-session; scripts and data in this directory reproduce them.

## Setup

A Ramsey(5,5,n) graph is an n-vertex graph with no clique of size 5 and no
independent set of size 5. `R(5,5)` is the least n for which none exists.
Current published state: `43 <= R(5,5) <= 46` (lower: Exoo 1989; upper:
Angeltveit–McKay 2024, arXiv:2409.15709).

Notation: for a vertex v of graph F, `F_v^+` = subgraph induced on N(v)
(order d(v)), `F_v^-` = subgraph induced on non-neighbors of v (order
n−1−d(v)). `e45(m)/E45(m)` = min/max edge counts over all Ramsey(4,5,m)
graphs.

## Input data and verification status

| Fact | Source | In-session verification |
|---|---|---|
| 328 (+ complements = 656) Ramsey(5,5,42) graphs | McKay's data page | all 656 verified K5-free, I5-free |
| R(4,5,24) census: 352,366 graphs | McKay–Radziszowski / Angeltveit–McKay | count + edge distribution confirmed; every 100th graph + all with e ≥ 130 property-checked by `src/structural_constraints.py` (in-session additionally: 3000 random); full property validation = gate 1 |
| Extremal classes of R(4,5,m), m ≤ 23 | r45extreme.tar.gz | all 2,887 graphs in the outermost min/max classes verified |
| **Completeness** of census and extremal classes | published computations | NOT re-verified (infeasible locally); dependency noted below |
| R(4,5) = 25 | McKay–Radziszowski 1995 | assumed |

Extremal edge tables assembled (orders needed here):

| m | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
|---|---|---|---|---|---|---|---|---|
| e45(m) | 41 | 50 | 57 | 68 | 77 | 88 | 101 | 116 |
| E45(m) | 79 | 85 | 92 | 100 | 107 | 114 | 122 | 132 |

Order 24 computed from the full census; exactly **2 graphs attain e = 132**,
both 11-regular with 176 triangles.

## Theorem 1 (degree window)

In any Ramsey(5,5,n) graph, every vertex degree d satisfies
`n − 25 <= d <= 24`.

*Proof.* N(v) contains no K4 (else K5 with v) and no I5, so it is a
Ramsey(4,5) graph, hence d <= R(4,5) − 1 = 24. The non-neighborhood contains
no I4 (else I5 with v) and no K5, hence n − 1 − d <= R(5,4) − 1 = 24. ∎

**Corollary.** `R(5,5) <= 50` (at n = 50 the window is empty). Verified: the
656 graphs at n = 42 have degrees in {19,…,22} ⊂ [17,24].

## Theorem 2 (excess identity)

For EVERY graph F on n vertices,
`Σ_v [ e(F_v^-) − e(F_v^+) − ½·d(v)·(n − 2d(v)) ] = 0`.

*Proof.* Fix an edge xy. It lies inside F_v^+ for exactly c(x,y) vertices v
(the common neighbors) and inside F_v^- for exactly n − d(x) − d(y) + c(x,y)
vertices. Summing over edges,
`Σ_v [e(F_v^-) − e(F_v^+)] = Σ_{xy∈E}(n − d(x) − d(y)) = n·e − Σ_v d(v)²`,
which equals `½ Σ_v d(v)(n − 2d(v))`. ∎

Verified to hold exactly on all 656 Ramsey(5,5,42) graphs.

**Equivalent form: the neighbour-degree identity (recorded 2026-08-21).** The
proof line `Σ_v [e(F_v^-) − e(F_v^+)] = n·e − Σ_v d(v)²` is exactly the
aggregate of the per-vertex counting identity

`Σ_{u ∈ N(v)} d(u) = e + e(G[N(v)]) − e(G[D(v)])`  for every vertex v,

which follows from `e = d(v) + e_x(v) + e_D(v) + e_xz(v)` together with
`Σ_{u∈N(v)} d(u) = d(v) + 2 e_x(v) + e_xz(v)`. Summing and using
`Σ_v Σ_{u∈N(v)} d(u) = Σ_u d(u)²` gives
`Σ_v d(v)² = n·e + Σ_v e(G[N(v)]) − Σ_v e(G[D(v)])`. Verified on all
2,131,018 labelled graphs up to n = 7 and on random graphs at n = 12, 20, 45;
independently re-verified over all graphs on n = 4, 5, 6 (209,184 per-vertex
checks, zero violations).

**Consequence — this row cannot be added twice.** In the n = 45 cones the state
functional of the identity is `2d² − 45d − 2e_x + 2e_y` with
`e_x = E45[d] − a` and `e_y = C(m,2) − (E45[m] − b)`, and that expression equals
the stored `excess_balance(s)` identically on all 3215 states; the exact rank of
`{1, excess_balance, g_lo, g_hi, h_lo, h_hi, f_lo, f_hi}` is 8 before and after
adjoining it. So `excess_balance` *is* the degree-squared/handshake identity,
and no handshake-, assortativity- or degree-correlation row can strengthen the
relaxation at the aggregate level. Only the per-degree-class form is strictly
finer; the level-2 pair lift prices that refinement at 1/349
(see `mixed_campaign_result_2026-08-21.md` §7 and `src/pair_lift_mixed.py`).

**Trap.** Pairing the `+2b` term with the wrong stratum — i.e. substituting the
recorded `e_z` for `e(G[D(v)])` — yields
`2d² − 45d − 2(E45[d] − E45[m]) + 2a − 2b`, which *is* linearly independent of
the frozen basis (rank 8 → 9) and appears to cut all three routes hard. It is
not even an identity: on the path P₄ its general-n form sums to −4 instead of 0,
and it is nonzero on 2,016,606 of the 2,097,152 labelled graphs at n = 7.
Pinned by `tests/test_pair_lift_mixed.py::test_mis_keyed_row_is_not_even_an_identity`.

## h-interval tables

For a hypothetical Ramsey(5,5,n) graph, a degree-d vertex contributes
`h(v) = e(F_v^-) − e(F_v^+) − ½d(n−2d)` with
`h ∈ [C(m,2) − E45(m) − E45(d) − ½d(n−2d), C(m,2) − e45(m) − e45(d) − ½d(n−2d)]`,
m = n−1−d, because F_v^+ is a Ramsey(4,5,d) graph and F_v^- is the complement
of a Ramsey(4,5,m) graph. Theorem 2 forces Σ_v h(v) = 0. Computed exactly:

| n | degrees | h_min per degree |
|---|---|---|
| 43 | 18–24 | −4, −17/2, −13, −29/2, −13, −17/2, −4 |
| 44 | 19–24 | −5, −9, −11, −11, −9, −5 |
| 45 | 20–24 | −6, −15/2, −8, −15/2, −6 |
| 46 | 21–24 | −5, −5, −5, −5 |
| 47 | 22–24 | −3, −5/2, −3 |
| 48 | 23–24 | −1, −1 |
| 49 | 24 | 0 |
| 50 | (empty) | R(5,5) ≤ 50 immediately |

All h-intervals verified against the 656 real graphs at n = 42 (0 violations).

## Theorem 3 (deficiency budget)

Define the extremality deficiency of v as
`s(v) = [e(F_v^-) − (C(m,2) − E45(m))] + [E45(d) − e(F_v^+)] ≥ 0`
(how far v's neighborhood pair sits below the edge-extremal configuration).
Then `Σ_v s(v) = Σ_v (−h_min(d(v)))`, giving the *total deficiency budgets*:

| n | 49 | 48 | 47 | 46 | 45 |
|---|---|---|---|---|---|
| max Σ s(v) | **0** | 48 | 141 | 230 | 360 |

Interpretation: neighborhoods of hypothetical Ramsey(5,5,n) graphs are forced
to be nearly edge-extremal Ramsey(4,5) graphs, with an exact global budget.
This is the quantitative engine behind the published upper-bound proofs and
the precise measure of the hardness gradient from 49 down to 45.

## Theorem 4 (n = 49 forcing; classical result, re-derived and verified)

Any Ramsey(5,5,49) graph F is 24-regular with h(v) ≡ 0, hence:
- every neighborhood induces one of the **two** 132-edge Ramsey(4,5,24)
  graphs X1, X2 (both 11-regular, 176 triangles, Σdeg² = 2904);
- every dual induces the complement of X1 or X2 (12-regular, 144 edges);
- therefore every edge has exactly 11 common neighbors (λ = 11) and every
  non-edge exactly 12 (μ = 12):

**F must be a strongly regular graph SRG(49, 24, 11, 12)** — conference-graph
parameters — with clique and independence number ≤ 4.

Consistency checks performed (all pass, showing this configuration is
remarkably self-consistent): triangle/K4 divisibility, the second-moment
identity Σdeg² = 2904, eigenvalue interlacing (both X1, X2 have non-principal
spectrum exactly within [−4, 3], endpoints attained).

## Theorem 5 (type matching — new in this session's presentation)

In any Ramsey(5,5,49) graph, the neighborhood type and the dual type of each
vertex must MATCH: `N(v) ≅ X_i  ⇒  D(v) ≅ complement(X_i)`.

*Proof.* The bipartite adjacency M between N(v) and D(v) satisfies
`M·Mᵀ = R(X_a)` and `Mᵀ·M = C(X_b)`, both fully determined integer matrices
(diagonal 12; off-diagonal 10/11/12 minus co-degrees as forced by λ = 11,
μ = 12). Since spec(MMᵀ) = spec(MᵀM) always, the traces of powers must agree.
Exact integer computation gives
`tr(R(X1)²) = 21984 ≠ 22032 = tr(C(X2)²)` (and symmetrically), eliminating
both mixed combinations. ∎

The two diagonal combinations pass every counting, divisibility, spectral,
and Gram-consistency test attempted. Their elimination (known to follow from
McKay–Radziszowski 1997) requires the gluing computation; a strengthened SAT
encoding (with the additional exact row-constraints
`e(D[M-row of w]) = 22 + ℓ_X(w)` derived from applying the 132-edge forcing
to neighbors of v0) is preserved in `src/n49_gluing_sat.py` but was NOT
run to completion (compute budget; a weaker variant did not resolve within
50 minutes).

## Feasibility analysis for the open frontier n = 45

The published R(5,5) ≤ 46 proof needed catalogs of Ramsey(4,5,m) graphs with
deficiency ≤ ~5 and ~80 CPU-years. At n = 45 the budget rises to 8 per vertex
(and a fifth degree class, d = 20, opens). Observed catalog growth per unit
of edge-deficiency at order 23: ×43–66 per step
(122: 2 → 121: 119 → 120: 7,800 → 119: 332,778). Extrapolating conservatively
(×40/step), catalogs to deficiency 8 at order 23 reach ~10^10–10^13 graphs.

**Conclusion (evidence-based):** determining R(5,5) ≤ 45 by the current
method is compute-infeasible by 3–5 orders of magnitude; a genuinely new
structural idea — e.g., a second identity family tightening the deficiency
budget, or forcing statements analogous to Theorems 4–5 at deficiency > 0 —
is the actual bottleneck. This is where research effort should aim.

## Dependency ledger

- Unconditional (given standard published values R(4,5) = 25 and the
  completeness of the order-24 census and extremal classes): Theorems 1–5,
  all tables.
- Verified in-session: all graph data used; identity and intervals against
  656 real Ramsey graphs; the two 132-edge extremal graphs and their spectra.
- Not verified in-session: completeness of the published censuses (the hard
  published computations); the final elimination of SRG(49,24,11,12)
  candidates (published: McKay–Radziszowski 1997; superseded by ≤ 48, ≤ 46).

## Reproduction

- `src/structural_constraints.py` — end-to-end pipeline (verification, tables, budgets,
  n = 49 chain, trace-power proof).
- `src/n49_gluing_sat.py` — optional, unresolved gluing SAT (do not run
  casually; expect hours+).
- `data/structural_tables.json` — all computed numbers.
- Data: `data/` (validated in gate 1; SHA-256s in `data/VALIDATION.json`)
  (McKay, https://users.cecs.anu.edu.au/~bdm/data/ramsey.html).

## m=3 identity cone result (2026-08-17; gate recorded, machine-checked)

Campaign `higher_order_identity_positive_deficiency_m3` completed with the
frozen negative disposition **`M3_NO_CUT_IN_FROZEN_CONE`**. Artifact:
`data/higher_identity_m3.json` (SHA-256 `8295ab4d...`; 3,215-state cone over
degrees {20..24}, exact q-histogram windows for all 53 catalog classes and
sound outer windows elsewhere; every number exact rational arithmetic).

What is PROVED (given the recorded catalog trust roots and the m=3 identity,
itself exhaustively machine-verified on all graphs through order 6 and
replayed on all 656 known R(5,5,42) graphs and complements):

1. The frozen continuous cone over the (d, a, b, excess-balance, g-interval)
   state space admits NO exact affine certificate meeting any of the three
   quantified success routes. Each route carries both an exact rational
   certificate attaining its optimal bound and an exact rational primal
   obstruction witness attaining the same value:
   - `total_deficiency`: bound 360 (certificate `alpha=8, beta=1/2`), witness
     states {1324: 155/4, 1620: 25/4} with value 360; the 315 target is
     unreachable in this relaxation.
   - `degree20_count`: bound 14310/349 (certificate `alpha=318/349,
     beta=6/349, delta=1/1047`), witness {85: 38745/1396, 153: 18495/1396,
     2657: 1395/349} with the same value; a bound < 1 is unreachable.
   - `deficiency_ge8_count`: bound 45 (certificate `alpha=1, beta=0`),
     witness {1407: 45} with value 45; a bound < 1 is unreachable.
2. The route bounds are the exact LP optima of the frozen cone, certified
   twice by disjoint machinery (the producer's exact verifiers and an
   independent exact-simplex audit; see `notes/cone_lp_robustness_2026-08-17.md`).
3. The complete evidence chain was re-verified by the disjoint checker
   `src/check_m3_certificate.py` (stdlib + graph6 parsing only): 58 input
   hashes byte-identical, 656/656 replay residuals zero under its own kernel,
   all 8,500,211 catalog graphs reswept into the same 53 class histograms, all
   3,215 states rebuilt exactly, every certificate/witness/status re-derived:
   `M3 EVIDENCE VERIFIED: M3_NO_CUT_IN_FROZEN_CONE`.

What this does NOT say: no theorem about all possible uses of the m=3
identity is claimed; the negative is specific to the frozen cone, its frozen
objective registry, and the continuous relaxation. No `R(5,5)` bound is
claimed or changed. The `required_local_family` route stays explicitly
`UNAVAILABLE_NO_COVER_CERTIFICATE` (no hash-pinned complete gluing-cover
manifest exists in this repository).

Measured consequence for the next step (`notes/cone_lp_robustness_2026-08-17.md`):
the cone is LP-robust — envelope sharpening, exact g-signs, census-class
tightening, and parity/integrality refinements all leave routes 1 and 3
unchanged and move route 2 by less than 10 units; the only measured flip
channel inside reach of new mathematics is degree-class exclusion. The
separately reviewed m=4 plan (`notes/higher_order_identity_m4_plan_2026-08-17.md`)
adds the last nontrivial T-family row (per-state t, diamond, K3+K1, i4(dual)
intervals) under the same exact-certificate + disjoint-checker discipline.

## 2026-08-18 - THE M=4 ROW CLOSES THE T-FAMILY LINE AT n=45 (exact negative)

The m=4 campaign (`notes/higher_order_identity_m4_plan_2026-08-17.md`) is
complete: identity kernel exhaustively verified, catalog motif vector tables
built, envelope LP certified, cone extended, cuts searched, disposition
`M4_NO_CUT_IN_FROZEN_BASIS` recorded in `bench/spec.json
next_research_campaign.status`.

The exact statement (trust-rooted on the named catalog completeness claims,
with the independent checker replaying every step):

1. The frozen continuous cone over the nine-field state space
   (d, a, b, excess-balance, g-interval, h-interval) - the m=2/m=3 rows
   reused as-is plus the m=4 row
   h_v = 3(47-2d_v) t(N_v) + 4 diamond(N_v) - 6 s(K3+K1, N_v)
         - 12 i4(y-class) with plain stratum classes -
   admits NO exact affine certificate meeting any of the three quantified
   success routes. Each instantiable route carries an exact rational
   certificate and an exact rational obstruction witness at the same value:
   - `total_deficiency`: bound 360, witness states {1406: 175/4,
     1780: 5/4}; the 315 target unreachable in the frozen basis.
   - `degree20_count`: bound 14310/349, witness {85: 38745/1396,
     153: 18495/1396, 2657: 1395/349}; a bound < 1 unreachable.
   - `deficiency_ge8_count`: bound 45, witness {1251: 45}; a bound < 1
     unreachable.
   - `required_local_family` stays UNAVAILABLE_NO_COVER_CERTIFICATE (no
     hash-pinned complete gluing-cover manifest exists in this repository).
2. The three accepted certificates are EXACTLY the m=3 certificates: every
   accepted certificate carries epsilon = zeta = 0, so the m=4 row adds no
   cut strength in the frozen basis. This is a theorem of the recorded
   relaxation (verified twice), not a numerical trend.
3. The envelope LP (vertices of {s >= 0 : R1, R2 exact; R3-R5 two-sided}
   over the eleven induced 4-motifs) is enumerated exactly and provably
   completely (the basis cut is singular-row elimination, verified against
   the C(17,9) reference); all 130 strata classes were certified (38
   catalog + 92 envelope_lp) with ZERO fallbacks and ZERO audit violations
   across 8,500,211 streamed graphs.
4. The full chain was re-derived by the disjoint checker
   `src/check_m4_certificate.py` (stdlib + graph6 parsing only):
   61 input hashes byte-identical, 656/656 replay residuals zero, n=49
   anchors (row 1584 = 12*132, observed i4 multiset [138,144] above the
   forced mean 132), 8,500,211 graphs reswept to the same 53 class
   histograms, 3,215 states rebuilt exactly, all search evidence
   re-verified against the REBUILT states:
   `M4 EVIDENCE VERIFIED: M4_NO_CUT_IN_FROZEN_BASIS`.

What this does NOT say: no theorem about all possible uses of the m=4
identity is claimed; the negative is specific to the frozen basis, the
frozen objective registry, and the continuous relaxation. No R(5,5) bound
is claimed or changed. The checker trusts the PRODUCER's envelope-LP
windows for the 92 missing classes (caps + schema only there; exact
soundness for present classes) - scoped in its docstring; a future
ACCEPTED-tier disposition under any follow-on campaign requires the
plan-letter R1-R5 LP recompute before the evidence label holds.

Consequence: the T-family identity line (m=2, 3, 4; m>=5 void) is closed
at n=45. The Engström mixed row and its level-2 pair lift have since also
closed (`MIXED_NO_CUT_IN_FROZEN_BASIS`, `MIXED_PAIR_LIFT_CLOSED`). The
remaining exclusion lane is non-local certified gluing. The complementary
construction lane is addressed next.

## Theorem 6 (no order-45 Cayley Ramsey graph; 2026-08-22)

No undirected Cayley graph on 45 vertices has both clique number and
independence number at most four.

*Proof.* Let `G` be a group of order 45. Its Sylow-5 subgroup is normal because
its number divides 9 and is 1 modulo 5. Conjugation by a Sylow-3 subgroup has
image whose order divides both 9 and `|Aut(C5)|=4`, so it is trivial. Thus
`G = C5 × H`, where the order-9 group `H` is `C9` or `C3 × C3`. Hence the only
two groups are `C45` and `C15 × C3`, both abelian.

Neither has an involution, so the 44 nonidentity elements form 22 inverse
pairs. An undirected Cayley graph selects an arbitrary subset `S` of those
pairs. Translation makes any K5 contain zero; therefore `Cay(G,S)` has a K5
iff four elements of `S` form a K4 under pair differences. Its complement is
`Cay(G, G\\({0}∪S))`. Exhausting one representative of every complementary
pair of connection sets gives `2^21` cases per group. In every one, the
producer found and recorded a K5 on one side. A disjoint derivation converts
the same condition to a 22-variable CNF and checks all `2^22` assignments by
exact integer bitsets, leaving zero survivors for both groups. ∎

Exact evidence:

| group | normalized K5 patterns | CNF clauses | complement classes | survivors |
|---|---:|---:|---:|---:|
| `C45` | 12,749 | 25,498 | 2,097,152 | 0 |
| `C15 × C3` | 12,650 | 25,300 | 2,097,152 | 0 |

`src/cayley_r55_search.py` commits all 4,194,304 monochromatic witnesses to
two SHA-256 coverage digests. `src/check_cayley_r55.py --resweep` independently
repeats the full search with a different K4 kernel and reproduces both digests.
The compact CNF checker is the default proof path and takes about three seconds.

**Conference corollary.** No Cayley graph has strongly regular parameters
`(45,22,10,11)`. Both groups have a character of order three. Since `S=-S`,
each inverse pair contributes either 2 or `ω+ω²=-1`, making that character's
Cayley eigenvalue an integer. But every nonprincipal eigenvalue of an
`srg(45,22,10,11)` solves `x²+x−11=0`, whose roots are irrational.

**Scope.** The theorem eliminates every regular group action, including every
circulant construction, but not all vertex-transitive graphs of order 45 and
not arbitrary graphs. It changes no bound on `R(5,5)`. A targeted literature
search found no published exact statement of this order-45 Cayley closure;
that absence is not a novelty proof.

## Theorem 7 (order-3 automorphism parity in the conference lane)

Let `A` be the adjacency matrix of any `srg(45,22,10,11)`, and let an
automorphism of order three have `f` fixed vertices. Then

`f ≡ 3 (mod 6)`.

In particular, the graph has no fixed-point-free automorphism of order three.
This strictly strengthens the Cayley conference corollary: translation by an
element of order three would be fixed-point-free.

*Proof.* The SRG relation is

`A² + A − 11I = 11J`.

Over `K=Q(ω)`, decompose the permutation representation of the automorphism
into its `1,ω,ω²` eigenspaces. If there are `c=(45−f)/3` three-cycles, the
`ω`-eigenspace `W` has dimension `c`. It is `A`-invariant, and `J` vanishes on
it, so `A|W` satisfies `x²+x−11=0`. This polynomial is irreducible over
`Q(ω)=Q(√−3)` because its discriminant 45 is not a square in that field.
Therefore `W` is a vector space over the quadratic extension
`K[x]/(x²+x−11)` and `c` is even. Hence `f=45−3c ≡ 3 (mod 6)`. ∎

The statement is an elementary in-repo derivation, not a novelty claim.
Maksimović's 2018 orbit-matrix enumeration studies
`srg(45,22,10,11)` graphs with `S3` automorphisms and is the closest located
prior-art lane; the exact congruence above was not source-verified there in
this session.

## Theorem 8 (no order-three automorphism in the Ramsey conference lane)

Let `G` be a strongly regular graph with parameters `(45,22,10,11)` and with
no clique or independent set of order five. Then `G` has no automorphism of
order three. Equivalently, `3` does not divide `|Aut(G)|`.

*Proof for `f=3`.* Every fixed-graph degree is `1 mod 3`; on three vertices it
must equal one. A graph cannot have three odd degrees.

*Proof for `f=9`.* Let `H` be the fixed graph. Its degrees lie in `{1,4,7}`.
A degree-one vertex and any nonneighbor share at most one fixed neighbor, but
their total common-neighbor count is 11 and all nonfixed orbits contribute
multiples of three, requiring `2 mod 3` fixed common neighbors. Thus degree
one is impossible; complementation excludes degree seven. Hence `H` is
4-regular. Adjacent fixed pairs have exactly one fixed common neighbor and
nonadjacent pairs exactly two, so `H = srg(9,4,1,2) ≅ L_2(3)`.

There are twelve 3-cycles. The trace calculation in Theorem 7 forces six to
induce triangles and six to induce independent triples. A fixed vertex has
four fixed neighbors and therefore meets exactly six orbit columns. A triangle
orbit's fixed-neighbor set must be independent (else it completes a K5), while
an independent orbit's fixed-nonneighbor set must be a clique (else it
completes an independent five-set). Since `L_2(3)` has independence and clique
number three, the six triangle columns all have size three and the six
independent columns all have size six. Exact multiplicity equations force each
of the six independent triples and each complement of the six triangles of
`L_2(3)` to occur once. This gives a unique `9×12` incidence matrix `B`, up to
the fixed canonical ordering, with `BB^T=3(I+J)`.

For the equitable quotient

`Q = [[A, 3B], [B^T, T]]`,

the fixed-to-orbit SRG equations, diagonal `diag(T)=(2^6,0^6)`, and row sums
`(19^6,16^6)` have exact rank 77 in the 78 symmetric entries of `T`. Their
entire affine family has entries `T[0,7]=τ` and `T[0,6]=3−τ`; integer bounds
give `τ∈{0,1,2,3}`. The cycle-to-cycle equation

`3B^T B + T^2 + T − 11I = 33J`

leaves exactly `τ=1,2`.

In both matrices `T[0,3]=1`. Independently rotate orbit 3 relative to orbit 0
so that this single matching has phase zero. Fixed vertices `{2,3,6}` are
independent and miss both orbit columns 0 and 3; positions 0 and 1 in those
orbits are nonadjacent. These five vertices are independent, contradiction.

*Proof for `f≥15`.* Put `c=(45−f)/3`. Theorem 7 makes `c` even. On each
nontrivial order-three eigenspace, irreducibility of `x²+x−11` gives the two
nonprincipal adjacency eigenvalues multiplicity `c/2` each. Hence on the
orbit-constant subspace they have multiplicity `22−c` each, and the trace of
the adjacency quotient there is

`22 + (22−c)(r+s) = 22 − (22−c) = c`.

Each 3-orbit induces either a triangle, contributing `2` to the quotient
trace, or an independent triple, contributing zero. Therefore exactly `c/2`
of the 3-orbits are triangles.

Choose a triangle orbit `O`. Let `s` fixed vertices be adjacent to it. For the
other `c−1` three-orbits, let `t_i∈{0,1,2,3}` be the number of neighbors in
that orbit of one vertex of `O`. Equal orbit sizes make the reciprocal quotient
entry also `t_i`. The degree equation gives

`Σ_i t_i = 20−s`.

Two adjacent vertices in `O` have ten common neighbors. Their fixed common
neighbors contribute `s`, the third vertex of `O` contributes one, and a
degree-`t_i` cyclic bipartite graph between two 3-orbits contributes
`C(t_i,2)`. Thus

`Σ_i C(t_i,2) = 9−s`,

so `s≤9` and `Σ_i t_i² = 38−3s`. No Ramsey hypothesis entered these
identities.

For `f≥15`, `c≤10`. Cauchy–Schwarz would require

`(20−s)² ≤ (c−1)(38−3s) ≤ 9(38−3s)`.

But the difference at `c=10` is `s²−13s+58`, whose discriminant is `−63`;
it is positive for every real `s`. Decreasing `c` only increases the gap
because `38−3s=Σ_i t_i²≥0`. This contradiction excludes
`f=15,21,27,33,39` for every `srg(45,22,10,11)`, Ramsey-good or not.
Theorem 7 leaves only `f=3,9,15,21,27,33,39` for a nonidentity order-three
automorphism (`f=45` is the identity), so the three cases above exhaust all
possibilities in the Ramsey-good lane. ∎

**Corollary.** No Ramsey-good `srg(45,22,10,11)` is vertex-transitive. A
transitive automorphism group on 45 vertices has order divisible by three,
and Cauchy's group theorem would supply an automorphism of order three.

Exact artifacts: `src/order3_srg_frontier.py`,
`src/check_order3_srg_frontier.py`, `tests/test_order3_srg_frontier.py`, and
`data/order3_srg_frontier.json` (schema 2). The producer uses only
integer/Fraction arithmetic. The checker independently constructs `L_2(3)`
as a rook graph, repeats the complete `f=9` census and quotient calculation,
derives the spectral trace from `(45,22,10,11)`, and derives both triangle
moments from `k=22` and `λ=10` before checking every integral support.

**Closest source-verified prior art.** Maksimović, *Symmetry* 10 (2018), 212,
Table 4, found seven `Z3` orbit matrices and 288
`srg(45,22,10,11)` graphs from the sole viable `S3` orbit distribution
`(d1,d2,d3,d6)=(1,4,4,4)`. Restricting that action to `Z3` gives nine fixed
vertices and twelve 3-orbits, exactly the `f=9` lane above. The published
[adjacency data](https://www.math.uniri.hr/~mmaksimovic/srg45.txt) expose the
same cycle structure. The paper does not state the Ramsey-good exclusion or
the `f≥15` moment argument. This is a scope comparison, not a novelty
certification: https://doi.org/10.3390/sym10060212.

**Scope.** This excludes order-three symmetry—and therefore vertex
transitivity—only inside the Ramsey-good strongly-regular construction lane.
Non-strongly-regular hypothetical Ramsey(5,5,45) graphs remain open. No
Ramsey-number bound is changed.

## Theorem 9 (odd-prime automorphisms are impossible in the Ramsey conference lane; 2026-08-23)

Let `G` be a strongly regular graph with parameters `(45,22,10,11)` and with
no clique or independent set of order five. Then

`|Aut(G)|` is a power of two.

Equivalently, `G` has no automorphism of odd prime order. Theorem 8 supplies
the prime three case. The remaining primes are excluded as follows.

### Common cyclotomic reduction

Write the nonprincipal eigenvalues as

`r,s = (−1±3√5)/2`,

each with multiplicity 22. If an automorphism of odd prime order `p` has `c`
moving `p`-cycles, each nontrivial character space has dimension `c` and is
annihilated by `J`. For `p≠5`, `Q(ζ_p)` does not contain `Q(√5)`, so
`x²+x−11` remains irreducible and `c` is even. This leaves

- `p=7`: `(f,c)=(3,6),(17,4),(31,2)`;
- `p=11`: `(f,c)=(1,4),(23,2)`;
- `p=13,17,19`: only `c=2`, giving `f=19,11,7`;
- `p≥23`: no nonidentity case, since at most one `p`-cycle fits but `c` must
  be even.

The prime five field is exceptional because `√5∈Q(ζ_5)`, so it needs a
separate trace argument.

### Order five

On the four nontrivial character spaces, Galois symmetry gives multiplicities
`a_1=a_4`, `a_2=a_3`, and `a_1+a_2=c`; it does **not** force `c` even. On the
orbit-constant subspace, both nonprincipal adjacency eigenvalues have
multiplicity `22−2c`, so the quotient trace is `2c`.

Every moving 5-orbit induces an invariant graph of degree `0,2`, or `4`. If
their counts are `n_0,n_2,n_4`, the trace equation is

`2n_2+4n_4 = 2c = 2(n_0+n_2+n_4)`,

and hence `n_0=n_4`. Ramsey avoidance forbids both an independent orbit
(`d=0`) and a complete orbit (`d=4`), so in the Ramsey lane every moving orbit
is a `C5`.

There is a stronger Ramsey-free cut for `f≥10`. Choose any moving orbit, with
internal degree `d∈{0,2,4}`, fixed support `s`, and quotient degrees `t_i` into
the other `c−1` moving orbits. Degree and the diagonal quotient equation give

`Σt_i = 22−d−s`,  
`Σt_i² = 66−d−d²−5s`.

For `f≥10`, `c≤7`. At the weakest value `c=7`, the Cauchy gaps for
`d=0,2,4` are respectively

`s²−14s+88`, `s²−10s+40`, `s²−6s+48`,

with discriminants `−156,−60,−156`. All are positive for every real `s`.
Smaller `c` only increases the gap because the square sum is nonnegative.
Thus `f=10,15,20,25,30,35,40` are impossible for every
`srg(45,22,10,11)`, without a Ramsey hypothesis.

Only `f=0,5` remain.

*Fixed-point-free case.* There are nine `C5` orbits. Let `B` be their symmetric
equitable quotient, so `diag(B)=2`, `B1=22·1`, and

`B²+B = 11I+55J`.

Set `M=2B+I−5J`. Then

`M1=0`, `diag(M)=0`, `M²=45I−5J`,

and every off-diagonal entry lies in `{-5,−3,−1,1,3,5}`. Each row therefore
has sum zero and square sum 40. Up to permuting the other eight orbits, its
off-diagonal entries have exactly one of three forms:

`(−3,−3,−1,−1,1,1,3,3)`,  
`(−3,−1,−1,−1,−1,1,1,5)`,  
`(−5,−1,−1,1,1,1,1,3)`.

A complete row-by-row integer enumeration tests 56,755 prefixes. In the three
first-row cases, respectively, `72,64,8` second rows survive; `132,96,48`
three-row prefixes survive; and no fourth row can be added. Hence no such
quotient exists.

*Five fixed points.* Every fixed degree is `2 mod 5`, so the fixed graph is a
`C5`. Each fixed vertex is adjacent to four of the eight moving orbits. If
those fixed-neighbor supports are viewed as eight blocks on five points, every
point has replication four and every pair occurs in exactly two blocks.

For a fixed triple `T`, let `m(T)` be the number of blocks containing it. If
`m(T)=0`, choose an edge of the fixed `C5` in `T`; its two blocks both omit the
third point. Their ten moving vertices contain a triangle or an independent
four-set by `R(3,4)=9`, which extends to a `K5` or `I5`. If `m(T)≥2`, choose a
nonedge of the fixed `C5` in `T`; inclusion-exclusion shows that `m(T)` blocks
contain the remaining point and omit that nonedge. Two such blocks, together
with `R(4,3)=9`, again extend to a `K5` or `I5`.

Thus Ramsey avoidance would require `m(T)=1` for all ten triples. This is
impossible. If a block has size five, the other seven blocks have size at most
two and cannot supply the second occurrence of all ten pairs. Otherwise every
block has size at most four. For any pair, its two blocks must then have sizes
three and four, so every one of the ten pairs occurs in exactly one size-three
block; but size-three blocks account for pairs in multiples of three. This
excludes `f=5` and completes order five.

### Order seven

For any two fixed vertices, their moving-orbit supports intersect in at most
one orbit: two shared 7-orbits would already contribute 14 common neighbors.
The same holds for the complementary supports because the complement has the
same SRG parameters.

- `f=3`: every fixed degree is one, violating the handshake lemma.
- `f=17`: support sizes are `1,2,3`. At most one size-one support, at most one
  size-three support, and at most one copy of each of the six size-two
  supports can occur. The simultaneous intersection constraints have exact
  maximum six (attained by the six size-two supports), not 17.
- `f=31`: the two moving orbits have internal degrees two and four. Original
  and complementary intersection bounds leave, up to complementation, one
  fixed vertex adjacent to both orbits, none adjacent to neither, cross degree
  three, and support sizes 17 and 15. Let `B,C` be their disjoint residual
  support groups, of sizes 16 and 14. For the fixed adjacency matrix `H` and
  incidence matrix `X`,

  `H²+H = 11(I+J)−7XX^T`.

  Taking `v=1_B−1_C` gives `||v||²=30`, `Σv=2`, and
  `X^Tv=(16,−14)`, so the right-hand quadratic form is

  `11(30+2²)−7(16²+14²)=−2790`.

  But every eigenvalue of `H²+H` is at least `−1/4`, so the same form is at
  least `−30/4`. Contradiction.

### Order eleven

For `f=23`, fixed support sizes are `0,1,2`. Original and complementary
intersection bounds allow at most one vertex of support size two and at most
one of support size zero. At least 21 fixed vertices therefore have singleton
support. One of the two support classes has size at least 11, and any two
vertices in it share eleven moving common neighbors; they cannot be adjacent
because `λ=10`. This already gives an independent five-set.

For `f=1`, the fixed vertex meets exactly two of the four moving orbits.
Ramsey avoidance restricts an invariant 11-orbit to degree four or six, and
trace 20 forces two of each. The fixed-to-orbit and diagonal quotient
equations leave, up to swapping the last two orbits, the unique moving
quotient

```text
4 6 4 7
6 4 7 4
4 7 6 5
7 4 5 6
```

The two orbits adjacent to the fixed vertex each have internal degree four
and cross degree six. They have exactly

`C(5,2)² C(11,6) = 10²·462 = 46,200`

cyclic realizations. Exact bitset enumeration finds a `K4` in 44,000 and, in
the remaining 2,200, an `I5`. A `K4` extends with the fixed vertex to a `K5`.
The coverage digest is
`b888c63a7b77419c3e6df970acd3759391da073a3357f3318f11f0611d5538a6`.

### Larger primes and group conclusion

For `p=13,17,19`, the sole spectral candidates have fixed graphs respectively
9-regular on 19 vertices, 5-regular on 11 vertices, and 3-regular on seven
vertices. Their degree sums `171,55,21` are odd. For `p≥23`, parity of the
number of moving cycles excludes the only possible single cycle. No prime
larger than 43 divides the order of a subgroup of `S_45`.

Together with Theorem 8, every odd prime order is excluded. If an odd prime
divided `|Aut(G)|`, Cauchy's theorem would supply an element of that order.
Therefore `|Aut(G)|` is a power of two. ∎

Exact artifacts: `src/odd_prime_srg_frontier.py`,
`src/check_odd_prime_srg_frontier.py`,
`tests/test_odd_prime_srg_frontier.py`, and
`data/odd_prime_srg_frontier.json` (schema 1). The checker imports no producer
symbol: it reconstructs the order-five quotient search from row multisets,
derives the order-seven arithmetic independently, derives both order-eleven
quotients, and resweeps all 46,200 cyclic blocks with a second graph builder
and separate clique kernels. It composes the independently checked
order-three artifact rather than copying that proof.

As an external-data check, `src/validate_published_s3.py` parses Maksimović's
published GAP/GRAPE file without executing GAP. The source payload has
1,470,604 bytes and SHA-256
`a53366d919f5d29876d07031a65c1bd722d05e88b59b39a04e9988c9162ebec0`.
All 288 records verify as `srg(45,22,10,11)`, all admit the specified
order-three action with nine fixed points and twelve 3-cycles, and every
record contains both a `K5` and an `I5`. This corroborates Theorem 8's `f=9`
lane but is not used in any proof.

**Closest source-verified prior art.** Maksimović's 2018 paper classifies the
`S3` lane used by Theorem 8. Her 2023 survey/construction paper, Table 2, lists
the known `srg(45,22,10,11)` full automorphism groups; in particular eight
known graphs have group `Z10`, so order-five symmetry certainly exists outside
the Ramsey-good lane. The parameter set remains unclassified. Neither source
states the Ramsey-good odd-prime exclusion above. This is a scope comparison,
not a novelty certification:
https://doi.org/10.3390/sym10060212 and
https://doi.org/10.3390/sym15020408.

**Scope and next symmetry frontier.** The theorem applies only to a
Ramsey-good `srg(45,22,10,11)`. It leaves involutions—and therefore nontrivial
2-groups—open. It does not constrain a non-strongly-regular hypothetical
Ramsey(5,5,45) graph and changes no Ramsey-number bound.
