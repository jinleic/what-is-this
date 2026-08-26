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
All 288 records verify as `srg(45,22,10,11)` and admit both specified `S3`
generators: the order-three action has nine fixed points and the involution
has five. Every record contains both a `K5` and an `I5`. The order-three check
corroborates Theorem 8's `f=9` lane; the involution check is the positive
control expanded after Theorem 11. Neither external-data check is used in a
proof.

**Closest source-verified prior art.** Maksimović's 2018 paper classifies the
`S3` lane used by Theorem 8. Her 2023 paper, Table 2, lists
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

## Theorem 10 (involution fixed-count reduction for every conference SRG; 2026-08-24)

Let `G` be any strongly regular graph with parameters `(45,22,10,11)`, and let
`τ` be a nonidentity involution. Then the number `f` of fixed vertices belongs
to

`{1,5,9,13}`.

No clique or independence hypothesis is used.

*Proof.* Write the nonfixed vertices as `c=(45-f)/2` transposition pairs
`P_i={x_i,τ(x_i)}`. On the rational anti-invariant basis
`b_i=e_{x_i}-e_{τ(x_i)}`, adjacency has a symmetric integral matrix `C`.
If the invariant `2×2` block between two transposition pairs is

```text
a b
b a
```

then `C_ij=a-b∈{-1,0,1}`. On the diagonal,
`C_ii=-ε_i`, where `ε_i` records whether `P_i` is an edge. Since `J` vanishes
on the anti-invariant subspace, the SRG equation

`A²+A-11I=11J`

becomes

`C²+C=11I`.                                                     (1)

The polynomial `x²+x-11` is irreducible over `Q`. Its two roots therefore
have equal multiplicity on this rational representation, so `c` is even and
`tr(C)=-c/2`. Consequently exactly `c/2` pair orbits are edges. The diagonal
part of (1) gives

`Σ_{j≠i} C_ij²=11`

for every row. Thus the off-diagonal support graph `S` of `C` is 11-regular,
forcing `c≥12`. Hence `f≡1 (mod 4)`, `f≤21`, and initially

`f∈{1,5,9,13,17,21}`.                                         (2)

For `c=12`, `S=K_12`. Switch signs by `C↦DCD` so that the first row is
positive. Choose a second row from the same diagonal class. Its off-diagonal
equation with the first row requires a sum of ten signs to equal one,
impossible by parity. Hence `f=21` does not occur.

It remains to exclude `c=14`. Let `E` and `N` be the seven indices with
diagonal respectively `-1` and `0`, and let `Z` be the zero graph, the
off-diagonal complement of `S`. It is 2-regular. Reduce (1) modulo two and
write `ε` for the indicator of `E`, `ν=1+ε` for the indicator of `N`.
The row support degree is eleven, so `C1=ν`. Applying (1) to `1` gives
`Cν=ε`, and therefore `Cε=1`. Every vertex consequently has odd cross-part
degree in `S`, hence even cross-part degree in `Z`; because `Z` is 2-regular,
that degree is zero or two. Each zero cycle is therefore wholly inside one
part or alternates between the parts.

Up to permutations preserving `E,N`, the cycle lengths have exactly twelve
forms. If `x` vertices of each part lie on alternating cycles, then
`x∈{0,2,3,4,7}`. The forms are

```text
x=0:  E,N each 7 or 4+3                              (4 forms)
x=2:  cross 2; E,N each 5                            (1 form)
x=3:  cross 3; E,N each 4                            (1 form)
x=4:  cross 4 or 2+2; E,N each 3                     (2 forms)
x=7:  cross 7, 5+2, 4+3, or 3+2+2                   (4 forms).
```

The off-diagonal part of (1), again modulo two, imposes the
sign-independent rule

`|S(i)∩S(j)| ≡ 1`

exactly when `i,j` lie in the same part and `ij∈S`; in every other case the
intersection must be even. Each of the twelve zero-cycle forms violates this
rule. Their respective numbers of violating pairs are

`14,24,24,34,22,28,34,42,28,32,28,36`.

Thus `c=14`, and hence `f=17`, is impossible. Together with (2), only
`f∈{1,5,9,13}` remain. ∎

The finite artifact also supplies a redundant signing exhaustion. For each
zero form, switching fixes a lexicographic spanning tree of `S` positive,
which selects one representative of every switching class because `S` is
connected. Row-prefix search then checks `10,008,588` sign assignments across
the twelve `c=14` forms and finds no completion. A separately written checker
generates the forms by integer partitions and repeats the search with recursive
tail assignment.

### Fixed-support interface and the next exact gate

Let `H` be the fixed graph and let the columns of the binary matrix `X` be the
fixed-neighbor supports of the transposition pairs. The fixed-fixed block of
the SRG equation is

`H²+H+2XXᵀ=11(J+I)`.                                         (3)

Thus all fixed degrees are even, adjacent fixed pairs have an even number of
fixed common neighbors, and nonadjacent pairs have an odd number. Equivalently,

`H²+H+I=J (mod 2)`.

Equation (3) also gives the exact support-design parameters

```text
r(v)=(22-d_H(v))/2,
m(u,v)=(10-c_H(u,v))/2  if uv is an edge,
m(u,v)=(11-c_H(u,v))/2  otherwise.
```

Under the Ramsey-good hypothesis, an edge-pair support is triangle-free and a
nonedge-pair nonsupport has no independent triple. These implications do not
close the surviving cases. At `f=5`, exhaustive enumeration of all `2^10`
labelled fixed graphs leaves twelve labelled 5-cycles, one isomorphism type.
An explicit 20-column support design nevertheless survives: take all five
rotations of one subset of each size `1,2,3,4`, assigning even sizes to edge
pairs and odd sizes to nonedge pairs. Every point has replication ten and
every pair multiplicity five. The next genuine gate is therefore the signed
incidence completion of these support designs, not another fixed-count cut.

Exact artifacts: `src/involution_srg_frontier.py`,
`src/check_involution_srg_frontier.py`,
`tests/test_involution_srg_frontier.py`, and
`data/involution_srg_frontier.json` (schema 1). The producer and checker use
only the Python standard library and share no code.

**Closest source-verified prior art.** Behbahani--Lam's general fixed-point
bound, quoted as Theorem 1 in Maksimović (2018), gives only `f≤25` for these
parameters. Maksimović's classified `S3` and `Z6` actions realize involutions
with `f=5` and `f=13`, respectively, in known graphs, so those surviving
fixed-count models are not vacuous. No
involution-only orbit-matrix classification for `(45,22,10,11)` was located.
This is a scope comparison, not a novelty certification:
https://doi.org/10.1016/j.disc.2010.10.005,
https://doi.org/10.3390/sym10060212, and
https://doi.org/10.3390/sym15020408.

**Scope.** The fixed-count theorem holds for every `srg(45,22,10,11)`, but it
does not exclude involutions with `f=1,5,9,13`. It does not constrain
non-strongly-regular hypothetical Ramsey(5,5,45) graphs and changes no
Ramsey-number bound.

## Theorem 11 (complete fixed-five balanced-support relaxation; 2026-08-24)

Assume the involution fixes five vertices. By Theorem 10 and (3), the fixed
graph is `C5`. There are ten internally adjacent transposition pairs and ten
internally nonadjacent pairs. Their fixed supports have even and odd
cardinality, respectively. On the five fixed points every point has
replication ten and every point-pair has multiplicity five.

Let `𝓔` be the group of the sixteen even subsets of the fixed points under
symmetric difference. For `E∈𝓔`, let `x_E` count the even supports. Complement
each odd support and let `z_T` count its even complement `T∈𝓔`. Define the sign
character

`χ_E(U)=(-1)^|E∩U|`.

The constant Fourier coefficients of `x,z` are both ten. Replication ten
implies equality of the five singleton coefficients:

`\hat x({v})=\hat z({v})`.

Indeed, complementing an odd support negates every singleton sign. For a pair
`u,v`, the sum of `χ_S({u,v})` over all twenty supports is

`20-2r(u)-2r(v)+4m(u,v)=20-20-20+20=0`.

Complementation preserves pair signs, hence

`\hat x({u,v})=-\hat z({u,v})`.

The constant, five singleton, and ten pair character classes exhaust the
sixteen characters of `𝓔`. Fourier inversion is therefore exact. Its kernel
depends only on `|E△T|` and gives

```text
4z_T = sum_{|E△T|=2} x_E
       - sum_{|E△T|∈{0,4}} x_E.
```

Since `sum_E x_E=10`, put

`A_T=sum_{|E△T|∈{0,4}} x_E`.

Then

`z_T=(5-A_T)/2`.                                           (4)

The graph on `𝓔` joining subsets at symmetric-difference distance four is the
Clebsch graph, with spectrum `5,1^10,(-3)^5`. Hence its closed-neighborhood
matrix has spectrum `6,2^10,(-2)^5`; on the affine hyperplane
`sum_E x_E=10`, (4) is visibly an involution.

Consequently an even-support multiset has a nonnegative integral odd partner
if and only if every `A_T` lies in `{1,3,5}`. The Fourier sign change is an
involution, so applying (4) in the reverse direction shows that every valid
`x_E,z_T` lies in `{0,1,2}`.

Enumerating the resulting

`sum_{d=0}^5 C(16,d) C(16-d,10-2d)=996,216`

ternary vectors gives exactly **7,872** balanced, moment-feasible support
multisets on a labelled fixed `C5`. Quotienting by the explicit action of
`Aut(C5)=D5` gives exactly **844** orbits. Their orbit sizes are

```text
size 1:   2
size 5: 110
size 10: 732.
```

Burnside gives the same count: the identity fixes 7,872 designs, each of the
four nonidentity rotations fixes 2, and each of the five reflections fixes
112, so

`(7872+4·2+5·112)/10=844`.

Every support system realized by an SRG completion lies in this list, but the
converse is not asserted: 7,872 and 844 are upper bounds on the realizable
labelled supports and support orbits. ∎

The producer proves (4) directly and exhausts the 996,216 ternary vectors.
The disjoint checker does not import the producer and does not assume the
multiplicity bound. It builds the even and odd `16×16` incidence-moment
matrices, inverts the odd matrix exactly over the rationals, and recursively
covers all `C(25,15)=3,268,760` weak compositions of ten even blocks, pruning
only when no assignment of the remaining multiplicity can make an odd count
nonnegative. It independently recovers all 7,872 balanced designs and
generates the ten automorphisms of `C5` by filtering all 120 point
permutations.

### Signed completion coordinates

Equation (4) closes the balanced-support relaxation, not graph completion.
The remaining relations have a compact exact Seidel form. Let `N` be the
`5×20` fixed
incidence matrix, put `R_vi=1-2N_vi`, and let
`Q_F=J-I-2A(C5)`. For two transposition pairs let `p_ij,q_ij` be the parallel
and crossed adjacency bits. Define

```text
W_ii = 1-2a_i,             T_ii = 2a_i-1,
W_ij = 2-2(p_ij+q_ij),     T_ij = 2(q_ij-p_ij).
```

Here `a_i` records the internal pair edge. The conference Seidel equation
splits into

```text
Q_F R + R W = -J,
W² + 2RᵀR = 45I - 2J,                              (5)
T² = 45I.                                           (6)
```

Every row of `W` has eight nonzero off-diagonal entries and every row of `T`
has eleven. The relation blocks are recovered exactly by

```text
b_ij = 1-W_ij/2,  c_ij=-T_ij/2,
p_ij=(b_ij+c_ij)/2,  q_ij=(b_ij-c_ij)/2,
```

and their binary domain is equivalent to

`W_ij²+T_ij²=4` for `i≠j`.                           (7)

Thus (5)--(7), followed by the orbit-level `K5/I5` conditions, are the next
finite gate for the 844 balanced-support candidates.

As a source-verified positive control, all 288 published `S3` records have the
specified involution with five fixed vertices and satisfy (4)--(7). A
hash-bound dependency check verifies that their 26 support types occur among
the 844 balanced candidates. All 288 contain both a `K5` and an `I5`. This
checks the signed model on real completions but neither claims completeness of
the published family nor supplies a Ramsey-good completion.

Exact artifacts: `src/involution_f5_support_census.py`,
`src/check_involution_f5_support_census.py`,
`tests/test_involution_f5_support_census.py`, and
`data/involution_f5_support_census.json` (schema 1). Terminal disposition:
`F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS`.

**Scope.** Theorem 11 is an exact census of the necessary balanced-support
relaxation and therefore an upper bound on realizable support types. Signed
completion of its 844 candidates remains open, as do fixed counts `1,9,13`
and every non-strongly-regular graph. No Ramsey-number bound changes.

## Theorem 12 (forced fixed-window `R(3,3)` filter; 2026-08-24)

Retain the notation of Theorem 11. For each edge `e` of the fixed `C5`, let
`n(e)` be the unique fixed nonedge disjoint from `e`, and define

```text
M_e = {i : e is contained in S_i and S_i is disjoint from n(e)}.
```

The multiset cardinality `m_e=|M_e|` counts transposition orbits, not individual
vertices.

**Claim.** Every Ramsey-good completion satisfies

`m_e<=2` for all five fixed edges `e`.                         (8)

**Proof.** Each orbit in `M_e` contributes both of its vertices. Those vertices
are adjacent to both endpoints of `e` and nonadjacent to both endpoints of
`n(e)`. If `m_e>=3`, choose three orbits, giving six vertices. Since
`R(3,3)=6`, they contain a triangle or an independent triple. In the first
case the triangle together with the endpoints of `e` is a `K5`. In the second
case the independent triple together with the endpoints of `n(e)` is an
`I5`. Both contradict Ramsey-goodness. Therefore (8) holds. ∎

On the labelled fixed cycle, the five ordered window pairs are

```text
edge 01 versus nonedge 24,
edge 12 versus nonedge 03,
edge 23 versus nonedge 14,
edge 34 versus nonedge 02,
edge 04 versus nonedge 13.
```

Applying (8) to every representative in Theorem 11 rejects exactly **139** of
the 844 `D5` orbits, with orbit-size histogram `1:0, 5:29, 10:110`.
They account for **1,245** of the 7,872 labelled balanced supports. The exact
remaining necessary Ramsey-support frontier is therefore

```text
6,627 labelled supports,
705 D5 orbits,
orbit sizes 1:2, 5:81, 10:622.                              (9)
```

The rejected representatives have 176 violating windows in total; every
violating multiplicity is exactly three. The five-window multiplicity
histogram over all 844 representatives is

```text
m=0: 750,  m=1: 1816,  m=2: 1478,  m=3: 176.
```

The condition is invariant under the explicit ten-element `D5` action. The
producer checks all 844 canonical representatives and weights them by their
already certified orbit sizes. The disjoint checker hash-binds the complete
Theorem 11 artifact, exhausts all `2^15` labelled graphs on six vertices to
verify `R(3,3)=6`, uses `C5` as the five-vertex sharpness control, and
reconstructs every decision and aggregate in (9).

### Four-vertex relation consequence

Suppose `m_e=2`, with selected transposition orbits `i,j`. Their four vertices
must contain neither a triangle nor an independent triple. Let `a_i,a_j` be
the internal edge bits, and let `p_ij,q_ij` be the parallel and crossed bits
from Theorem 11. Exhausting their four possible values gives

```text
(a_i,a_j)   allowed (p_ij,q_ij)       allowed W_ij/2
(0,0)       (0,1),(1,0),(1,1)         -1,0
(0,1)       (0,1),(1,0)                0
(1,0)       (0,1),(1,0)                0
(1,1)       (0,0),(0,1),(1,0)          0,+1.               (10)
```

Across the 705 survivors, (10) supplies exactly **1,227** early relation
domains: 337 of type `{-1,0}`, 553 forced zeros, and 337 of type `{0,+1}`.
They are recorded with deterministic orbit indices obtained by expanding the
ten even supports in mask order and then the ten odd supports in mask order.

### Additional exact signed actions

The balanced moments imply

```text
R 1 = 0,       R R^T = 20 I_5.
```

Put `u=R^T 1_5`, so `u_i=5-2|S_i|`. Transposing the first equation in (5) and
using `Q_F 1_5=0` gives

`W u = -5 1_20`.                                             (11)

Multiplying the second equation in (5) by `1_20` gives
`W^2 1_20=5 1_20`. The same equation proves that `W` is invertible without
using (12): if `Wv=0`, applying `R` to
`2R^T Rv=45v-2Jv` and using `RR^T=20I` and `R1=0` gives `Rv=0`.
It follows that `45v=2(1^T v)1`; summing coordinates then forces `1^T v=0`
and hence `v=0`. Equations (11) and `W^2 1=5 1` now give
`W(W1+u)=0`, so

`W 1_20 = -u`.                                               (12)

Write `x_ij=W_ij/2` off the diagonal and let `epsilon_i` be the internal
edge bit, so `d_i=W_ii=1-2epsilon_i`. Since every row has eight nonzero `x_ij`,
(12) determines both sign degrees:

```text
deg_i(W=+2) = 3 + floor(|S_i|/2),
deg_i(W=-2) = 5 - floor(|S_i|/2).                            (13)
```

The total fixed-support incidence is 50 and exactly ten supports are odd.
Summing (13) therefore gives exactly 40 unordered `W=+2` relations, 40
unordered `W=-2` relations, and 110 complementary `T`-supported relations.

Finally, `U=span{1_20, im(R^T)}` has dimension six. The first equation in (5)
fixes the action of `W` on `U`, where its square is `5I` and its trace is zero.
On `U^perp`, the second equation gives `W^2=45I`; the total trace is zero
because the diagonal has ten entries of each sign. Hence every completion has

```text
char_W(t) = (t^2-5)^3 (t^2-45)^7,
char_T(t) = (t^2-45)^10.                                    (14)
```

Equations (10)--(14) are exact completion constraints and integrity checks.
They do not assert that any of the 705 supports has a signed completion.

Exact artifacts: `src/involution_f5_ramsey_filter.py`,
`src/check_involution_f5_ramsey_filter.py`,
`tests/test_involution_f5_ramsey_filter.py`, and
`data/involution_f5_ramsey_filter.json` (schema 1). Terminal disposition:
`F5_RAMSEY_R33_SUPPORT_FILTER_EXACT_6627_LABELLED_705_D5_ORBITS`.

**Scope.** Theorem 12 is a necessary Ramsey filter inside the fixed-five
strongly regular branch. Signed completion of the 705 survivors remains open,
as do fixed counts `1,9,13`, Ramsey-good strongly regular graphs with no
nontrivial involution, and all non-strongly-regular graphs. No Ramsey-number
bound changes.

## Computational status after Theorem 12 (2026-08-24)

The finite signed completion search has been exhausted, without yet promoting
its negative result to a theorem. For each of the 705 representatives, exact
meet-in-the-middle enumeration produces every row of `W` satisfying (5), the
sign degrees, and (10). There are 48,865,656 such rows in total. A symmetric
`W` is exactly a 20-clique in the resulting partite compatibility graph:
compatibility of rows `i,j` means agreement at `(i,j)` and their prescribed
inner product in (11). The diagonal entries of (11) are automatic from the
eight nonzero off-diagonal `W` entries.

The exact stages are:

```text
705 D5 supports / 6,627 labelled supports
    -- W square --> 155 / 1,401
    -- T square --> 108 /   986
    -- K5/I5  -->   0 /     0.
```

For a fixed `W`, the 11-regular support of `T` is connected. Conjugating `T`
by a diagonal sign matrix corresponds to swapping vertices inside
transposition orbits, so every switching class has a representative in which
the edges of a fixed spanning tree are positive. Exhausting that gauge yields
73,336 signed completions. Direct reconstruction finds both a `K5` and an
`I5` in every completion.

This establishes a complete **candidate-level** negative for the fixed-five
branch, conditional on the finite engines. It is not yet a theorem-grade
negative: the complete VeriPB/CakePB proof bundle is pending. The sound
certificate control for source representative 0 derives its DFS clauses
inside the proof (none is inserted as an unproved formula axiom), contains no
unchecked deletion, and is accepted by both pinned checkers. Coverage is
1 of 550 `W`-negative supports; 549 other `W` negatives and all 155 deeper
signed/Ramsey cases remain uncertified. Accordingly,
`data/involution_f5_ramsey_square_census.json` retains
`PENDING_VERIPB_CAKEPB`, and no theorem or Ramsey-number bound is claimed from
the exhaustive search alone.
