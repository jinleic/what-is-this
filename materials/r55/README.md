# R(5,5) ≤ 45 campaign

**Canonical target.** Prove that no Ramsey(5,5,45)-graph exists — a graph on 45
vertices with no clique of size 5 and no independent set of size 5 — i.e.

$$\mathcal{R}(5,5,45)=\emptyset \iff R(5,5)\le 45.$$

Secondary (constructive) direction: search for a Ramsey(5,5,43)-graph, which
would prove $R(5,5)\ge 44$ and refute the standing conjecture $R(5,5)=43$.
Either outcome at either order is a publishable exact result.

## Status quo (source-verified 2026-08-13, `notes/status_2026-08-13.md`)

- $43 \le R(5,5) \le 46$, confirmed current: Radziszowski's Dynamic Survey DS1
  **revision #18 (2026-04-24)** lists it verbatim as "the still open case".
  Lower: Exoo 1989. Upper: Angeltveit–McKay, arXiv:2409.15709 (v2 2025-09-01).
  The `math-open-problems-2026.md` figure "$\le 48$" is their **2018** result,
  superseded. Not scooped: a full arXiv ti:Ramsey sweep (2025-01→2026-08)
  found zero claimed improvements in either direction; AlphaEvolve
  (arXiv:2603.09172) improved nine other Ramsey lower bounds but not
  $R(5,5)\ge43$. (arXiv:2508.16699 *estimates* $R(5,5)=45$ heuristically —
  explicitly not a proof.)
- McKay–Radziszowski conjecture $R(5,5)=43$; 656 Ramsey(5,5,42)-graphs known
  (verified here: 328 in the file + 328 complement classes, disjoint under
  canonical labeling, none self-complementary) and no larger ones.
- Competition risk: Angeltveit is actively computing next door
  ($R(K_5,K_5{-}e)=30$, arXiv:2602.11459, 2026-02); a CMU group
  (Braginsky–Yellenki–Zhu w/ Mackey) is SAT-attacking srg(45,22,10,11), whose
  existence would prove $R(5,5)=46$ — 141/313 encodings UNSAT so far. Ruling
  that srg out is a natural sub-goal of any $\le45$ proof.
- AM state the ≤46 method cannot reach 45 without "new theoretical insights";
  Gasarch's [REPORTED] estimate for 45 by brute extension: ~3000 CPU-years.

## Why this frame is forced (PROVED, elementary)

In a Ramsey(s,t,n)-graph, the neighborhood of any vertex induces a
Ramsey(s−1,t)-graph and the non-neighborhood induces a Ramsey(s,t−1)-graph,
so $n - R(s,t-1) \le \deg(v) \le R(s-1,t)-1$. With $R(4,5)=25$:

- In a hypothetical Ramsey(5,5,45)-graph every degree lies in $[20,24]$;
  every neighborhood induces a Ramsey(4,5,d)-graph, $d\in[20,24]$, and every
  non-neighborhood induces the complement of a Ramsey(4,5,44-d)-graph.
- This is why the Ramsey(4,5) catalogs are the raw material of any attack,
  and why the Angeltveit–McKay proof glues pairs of (4,5)-graphs along a
  vertex ("pointed gluing" — exact terminology to be confirmed from the
  paper, see `notes/angeltveit_mckay_2409.md`).

## Structural constraints from the excess identity (PROVED + machine-checked, 2026-08-15)

Full writeup: `notes/structural_constraints_2026-08-15.md`. One-command
re-verification: `src/structural_constraints.py` (prints `ALL CHECKS PASSED`,
~2 min; exercises all 656 R(5,5,42) graphs and the boundary extremal classes;
for the 24-vertex census it verifies the exact count and edge distribution and
property-checks a deterministic sample plus all graphs with $e\ge130$ — full
census property validation is gate 1's `src/validate_catalogs.py`).

- **Excess identity** (PROVED, elementary double count): for every graph,
  $\sum_v [e(F_v^-) - e(F_v^+) - \tfrac{1}{2}d(v)(n-2d(v))] = 0$. Verified to
  hold exactly on all 656 known Ramsey(5,5,42) graphs.
- **Deficiency budgets** (PROVED, given the catalog trust root): the total
  distance-from-edge-extremality of all neighborhoods in a hypothetical
  Ramsey(5,5,n) graph is $\le 0/48/141/230/360$ at $n=49/48/47/46/45$ —
  the exact quantitative form of the hardness gradient, and the number any
  new inequality must beat to shrink the gate-4 search space.
- **n=49 forcing chain** (PROVED, given the trust root; classical result
  re-derived): any Ramsey(5,5,49) graph is an $\mathrm{srg}(49,24,11,12)$
  whose every neighborhood induces one of the exactly 2 edge-maximal
  Ramsey(4,5,24) graphs (both 11-regular, 176 triangles, non-principal
  spectrum exactly in $[-4,3]$).
- **Type-matching theorem** (PROVED, new presentation): $N(v)\cong X_i$
  forces $D(v)\cong \overline{X_i}$ — mixed combinations die by an exact
  integer Gram-trace mismatch ($\operatorname{tr}R^2 = 21984 \ne 22032$).
  The two diagonal combos survive every counting/spectral/divisibility test
  attempted; their elimination (published: MR 1997) is a gluing computation.
  A strengthened SAT encoding is preserved in `src/n49_gluing_sat.py`
  (NOT run to completion: a weaker variant was inconclusive after 50 min;
  do not run casually).
- **Catalog growth** (NUMERICAL extrapolation): observed $\times43$–$66$
  census growth per unit edge-deficiency at order 23; deficiency-8 catalogs
  (what a 46-style proof at $n{=}45$ would consume) extrapolate to
  $10^{10}$–$10^{13}$ graphs. Quantifies AM's "new theoretical insights
  needed" — the budget, not raw compute, is the object to attack.


## m=3 identity cone campaign (COMPLETE 2026-08-17: exact negative)

`M3_NO_CUT_IN_FROZEN_CONE`, triple-verified. Full record:
`notes/structural_constraints_2026-08-15.md` (final section), artifact
`data/higher_identity_m3.json`, audits `notes/cone_lp_robustness_2026-08-17.md`.
In one line: the exact 3,215-state cone over the m=3 neighborhood identity
row admits no accepted cut — each route has an exact optimal certificate
(360 / 14310/349 / 45) *and* an exact obstruction witness at the same value;
the disjoint checker `src/check_m3_certificate.py` re-derives the entire
chain from raw catalogs (`M3 EVIDENCE VERIFIED`, ~44 min). No R(5,5) bound is
claimed; the negative is specific to the frozen cone and routes.

Next: `notes/higher_order_identity_m4_plan_2026-08-17.md` (m=4 T-family row,
last nontrivial member; separately reviewed, all findings applied).

## m=4 identity cone campaign (COMPLETE 2026-08-18: exact negative)

`M4_NO_CUT_IN_FROZEN_BASIS`, double-verified. Full record:
`notes/structural_constraints_2026-08-15.md` (m=4 final section), artifact
`data/higher_identity_m4.json` (schema v2). In one line: the frozen
continuous basis `{balance, g, h}` — the m=2/m=3 rows plus the m=4 row
(per-state t, diamond, K3+K1, i4(dual) intervals over the R1-R5 motif
relation system, missing classes certified by exact envelope-LP
enumeration) — proves the SAME three route bounds as m=3 alone
(360 / 14310/349 / 45): every accepted certificate carries epsilon = zeta
= 0, so the m=4 row adds no cut strength in this basis. The T-family
identity line (m=2, 3, 4; m>=5 void) is closed at n=45; the Engstrom
mixed-identity tier stays open for a separately reviewed campaign.
The disjoint checker `src/check_m4_certificate.py` re-derives the entire
chain from raw catalogs with its own kernel (`M4 EVIDENCE VERIFIED`,
~40 min); the envelope windows it trusts for missing strata are scoped
exactly in its docstring.

## Go/no-go gates (in order; failing a gate stops the campaign)

1. **DONE 2026-08-13 — catalogs validated independently.**
   `src/check_ramsey.py` (dependency-free graph6 parser + clique checker,
   self-tested against brute force, C5, Paley(17)) confirmed, across
   3,785,907 graphs in 32 catalog files plus 8,241,382 graphs in the 77
   r45extreme files: every graph has its claimed Ramsey property, every
   published count matches exactly, all degrees respect the provable window.
   Records: `data/VALIDATION.json`, `data/VALIDATION_extreme.json`. Extra:
   - nauty `shortg`: 352,366 canonical classes in `r45_24.g6` — duplicate-free
     (same-author tool caveat: nauty is McKay's);
   - edge histogram of `r45_24.g6` reproduces the paper's Table exactly
     (132:2, 131:3, 130:32, 129:147, 128:843, 127:3401 ⇒ |A| = 4,428), and
     the r45extreme counts equal the paper's B/C/D strata sizes
     (119:332,778; 120:7,800; 121:119; 122:2; 113:30,976; 114:133; 107:31);
   - completeness spot-verified at small orders by from-scratch enumeration
     (`geng` + own filter): n=6:114, n=7:627, n=8:5,588, n=9:81,321 — all equal
     the catalog files;
   - R(5,5,42) complement-pair structure verified (656 = 328 + 328, disjoint).
   **Not verified (open trust root):** completeness of the R(4,5,24) catalog
   and of the edge-extremal censuses at n≥10 — that is the 2016/2024
   enumeration itself, which gate 2 starts to attack.
2. **DONE 2026-08-15 — one stratum of the ≤ 46 proof reproduced** from an
   independent implementation, held-out protocol (published file opened only
   at the final canonical diff). Warm-up: Lemma 5.1's census
   R(5,4,17,e≤59) = **7,147** graphs reproduced exactly, all with a vertex of
   degree ≥ 8. Main target: the **D3 census** $\mathcal{R}(4,5,21,e{=}107)$
   — reproduced as **31 iso classes, canonical set equal to the published
   file**, via a min-degree case split (δ≤9: one-vertex extensions of the
   complete order-20 classes e∈{98,99,100}; δ≥10: audited gluing with proved
   within-case pair windows; d=13 case empty by arithmetic). 209 raw
   solutions, every one independently re-verified; 77 min wall on 8 workers.
   Both frozen engines were measured infeasible at this stratum (CPU-years /
   ~8,000 CPU-h) — the case split is what made it tractable. A terminal
   min-degree leak in `glue_census2.c` was found, fixed, and regression-tested
   (identical output on the full (16,71) stratum at the true bound).
   Method + proofs: `notes/d3_case_split_2026-08-15.md`; driver:
   `src/case_split_21_107.py`; artifacts: `outD/case_split_21_107/`.
3. **DONE 2026-08-16 — held-out benchmark fixed and frozen (v2)**:
   [`bench/spec.json`](bench/spec.json) + [`bench/README.md`](bench/README.md).
   9 dev items + 2 held-out ((18,85)=74, (21,107)=31), all reproduced MATCH by
   fail-closed pipelines; baselines frozen with wall-clock **and** search-size
   metrics (C-DFS nodes per pair + Cadical accum conflicts/decisions for every
   SAT part — (21,107): 4,937 s, 71.18B nodes, 404.3M conflicts). All
   benchmark logic SHA-pinned (engine, sat_pair, harness, checker, D3 driver).
4. **STOPPED 2026-08-16 — no reproducible improvement transferred.**
   Protocol and dated erratum in `bench/README.md`; complete experiment ledger
   in `bench/spec.json`. Candidate 1 (static DFS→SAT handoff) was rejected by
   held-out A/B. Candidate 2 (adaptive handoff) was deprioritized after exact
   tail SAT measurement showed only an estimated ~6 % perfect-oracle ceiling.
   Candidate 3 implemented a proved intrinsic R(3,5,m)/R(4,4,m) local
   edge-window pair prefilter. Its fixed low-overhead variant preserved all
   nine dev MATCH results and reduced DFS nodes 1.62 %, but deferrals stayed
   191 and same-session full-dev ABBA improved wall only 43.45→43.15 s
   (−0.7 %). On untouched held-out (18,85), both arms reproduced 74 classes
   MATCH while wall regressed 252.0→266.9 s (+5.9 %); nodes fell only 0.014 %
   and deferrals stayed 3,148. Candidate 3 is rejected. Exact hashes, replay
   scripts, both candidate variants, and A/B metrics:
   `outD/candidate3/measurement.json`. Any restart requires a qualitatively
   stronger structural inequality or decomposition, not another scheduling
   retune or this local prefilter.

## Trust discipline (inherits `math/README.md`; additions)

- AI (this repo's agents included) stays **outside the trust base**: results
  count only as exact enumeration reproduced independently, or SAT/ILP with
  machine-checkable certificates (LRAT/VeriPB), verified by an independent
  checker. An optimizer or LLM saying "UNSAT/none found" is not a result.
- Symmetry breaking invalidates naive DRAT certification; certificate
  strategy (VeriPB / re-derivable symmetry lemmas) is a first-class design
  constraint, not an afterthought (see `notes/prior_art.md`).
- Every catalog file carries its SHA-256 in `data/VALIDATION.json`; published
  counts are quoted from the fetched source page, never from memory.

## Layout

| path | contents |
|---|---|
| `data/` | mirrored McKay catalogs + `download.sh` + `VALIDATION.json` |
| `src/check_ramsey.py` | trusted checker: graph6 + $K_s$/$I_t$ (selftest inside) |
| `src/validate_catalogs.py` | gate-1 runner against published counts + degree theorems |
| `src/structural_constraints.py` | excess identity, h-tables, budgets, n=49 srg chain (self-verifying) |
| `src/n49_gluing_sat.py` | optional unresolved n=49 diagonal gluing SAT (hours+; needs `python-sat`) |
| `notes/` | source-verified research notes (status, paper extraction, prior art, structural constraints) |
| `data/structural_tables.json` | machine-generated numbers behind the structural note |

## Running

```sh
cd math
./.venv/bin/python r55/src/check_ramsey.py --selftest
./.venv/bin/python r55/src/validate_catalogs.py
./.venv/bin/python r55/src/structural_constraints.py   # ALL CHECKS PASSED, ~2 min

# m=3 identity campaign (full reproduction; single-process, in order)
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*m3*.py' -v   # 111 tests
./.venv/bin/python r55/src/m3_deficiency_cone.py \
  --output r55/data/higher_identity_m3.json            # rebuild cone artifact
./.venv/bin/python r55/src/search_m3_cuts.py \
  --analysis r55/data/higher_identity_m3.json          # cut search -> NO_CUT_IN_FROZEN_CONE
./.venv/bin/python r55/src/check_m3_certificate.py \
  r55/data/higher_identity_m3.json                     # disjoint checker, ~44 min

# m=4 identity campaign (full reproduction; single-process, in order)
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*m4*.py' -v   # 111 tests
./.venv/bin/python r55/src/m4_deficiency_cone.py \
  --output r55/data/higher_identity_m4.json            # rebuild v2 cone artifact
./.venv/bin/python r55/src/search_m4_cuts.py \
  --analysis r55/data/higher_identity_m4.json          # cut search -> NO_CUT_IN_FROZEN_BASIS
./.venv/bin/python r55/src/check_m4_certificate.py \
  r55/data/higher_identity_m4.json                     # disjoint checker, ~40 min
```

The last m=4 command must print `M4 EVIDENCE VERIFIED:
M4_NO_CUT_IN_FROZEN_BASIS`; the last m=3 command must print `M3 EVIDENCE
VERIFIED: M3_NO_CUT_IN_FROZEN_CONE`.
Rebuilding an artifact from scratch reproduces it byte-identically only if
the frozen inputs are unchanged (they are hash-pinned inside the artifact;
the m=4 search records make identical rebuilds deterministic, not
byte-identical, since the discovery LP's proposal path is outside the trust
base and acceptance depends only on the exact re-verified certificates).

## Mixed Engström tier at n=45 (closed, 2026-08-21)

```bash
# mixed tier (full reproduction; single-process, in order)
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*mixed*.py' -v
./.venv/bin/python r55/src/mixed_deficiency_cone.py \
  --output r55/data/engstrom_identity.json             # rebuild v3 cone artifact, ~75 min
./.venv/bin/python r55/src/search_mixed_cuts.py \
  --analysis r55/data/engstrom_identity.json           # -> MIXED_NO_CUT_IN_FROZEN_BASIS
./.venv/bin/python r55/src/verify_mixed_independent.py \
  r55/data/engstrom_identity.json                      # independent verifier, sub-second
./.venv/bin/python r55/src/check_mixed_certificate.py \
  --sweep --jobs 32 r55/data/engstrom_identity.json    # disjoint checker + 8.5M resweep, ~6 min
./.venv/bin/python r55/src/pair_lift_mixed.py \
  r55/data/engstrom_identity.json                      # level-2 pair lift -> CLOSED
```

Terminal lines: `MIXED EVIDENCE VERIFIED: MIXED_NO_CUT_IN_FROZEN_BASIS` and
`LIFT VERDICT: MIXED_PAIR_LIFT_CLOSED`.

**Recorded negative.** Levels 1 and 2 of the first exact basis at n=45 are both
closed. The mixed Engström K=4 row is geometrically new (endpoint rank 6 → 8)
but has zero cut strength; the level-2 pair-coupled lift moves the three route
bounds from `360, 14310/349, 45` only to `360, 41, 45` against acceptance edges
`315, 1, 1`. Collapsing every catalog-missing motif window to a point still
leaves a factor-≈39 gap, so data quality is not the obstruction. The campaign's
`excess_balance` row *is* the degree-squared identity
`Σ_v d(v)² = n·e + Σ_v e(G[N(v)]) − Σ_v e(G[D(v)])` (Theorem 2, equivalent form),
which retires every handshake-style repair. Remaining routes, in order:
VeriPB-certified gluing search, then srg(45,22,10,11) complementary SAT.
**No bound on R(5,5) is claimed by any of this.**

## Cayley construction lane at n=45 (exactly closed, 2026-08-22)

Instead of adding another local inequality, this lane attacks the opposite side:
construct a Ramsey(5,5,45) graph with a regular group action. It exhausts the
class and proves the narrow theorem:

> **No undirected Cayley graph on 45 vertices is Ramsey(5,5).**

Every group of order 45 is `C45` or `C15 × C3`: its Sylow-5 subgroup is normal,
the order-9 conjugation action on `Aut(C5)=C4` is trivial, and a group of order
9 is `C9` or `C3 × C3`. Each group has 22 inverse pairs. Therefore all
undirected Cayley presentations comprise `2 × 2^22 = 8,388,608` connection
sets, or 4,194,304 complementary classes. No degree window, catalog, floating
point computation, or unproved R(4,5) input enters this result.

The producer attaches a canonical monochromatic K5 witness to every class. A
disjoint checker also derives the equivalent 22-variable CNF: 25,498 clauses
for `C45` and 25,300 for `C15 × C3`. Exact integer-bitset truth tables leave
zero of `2^22` assignments in either group. `--resweep` independently repeats
all 4,194,304 witness searches with a second clique kernel and reproduces both
coverage SHA-256 digests.

```bash
./.venv/bin/python -m unittest r55.tests.test_cayley_r55_search -v
./.venv/bin/python r55/src/cayley_r55_search.py \
  --output r55/data/cayley_r55_45.json                  # ~3.5 min
./.venv/bin/python r55/src/check_cayley_r55.py \
  r55/data/cayley_r55_45.json                          # compact exact theorem proof, ~5 s
./.venv/bin/python r55/src/check_cayley_r55.py --resweep \
  r55/data/cayley_r55_45.json                          # disjoint witness replay, ~3 min
```

Default terminal: `CAYLEY THEOREM VERIFIED: NO_CAYLEY_RAMSEY_5_5_45 ... witness_replay=NOT_RUN`. With `--resweep`, the
stronger artifact-level terminal is `CAYLEY EVIDENCE VERIFIED: ... witness_replay=VERIFIED`.

The conference-graph subcase has an even shorter obstruction. Either group has
an order-3 character; inverse-closed connection sets give that character an
integer Cayley eigenvalue, whereas `srg(45,22,10,11)` requires nonprincipal
eigenvalues solving `x²+x−11=0`. More generally, the order-three frontier
below now excludes every order-three automorphism of a Ramsey-good
`srg(45,22,10,11)`. Thus a successful conference/SAT construction, if one
exists, is necessarily non-Cayley, non-vertex-transitive, and has automorphism
group order coprime to three. This says nothing about non-strongly-regular
Ramsey(5,5,45) graphs and changes no bound on R(5,5).


## Order-three conference frontier (exactly closed, 2026-08-22)

The Cayley obstruction extends through the entire non-Cayley order-three
conference lane:

> **A Ramsey-good `srg(45,22,10,11)` has no automorphism of order three.**

The spectral argument first forces `f ≡ 3 (mod 6)`, leaving the nonidentity
fixed-point counts `3,9,15,21,27,33,39`.

- `f=3`: the fixed graph would be 1-regular on three vertices.
- `f=9`: the fixed graph is uniquely `L_2(3)`. Ramsey avoidance uniquely fixes
  its twelve fixed-to-orbit columns. Exact rank-77 quotient algebra leaves two
  matrices (`τ=1,2`), and both force the same independent five-set after one
  legal orbit rotation.
- `f≥15`: write `c=(45−f)/3≤10`. Spectral trace forces `c/2` triangle
  3-orbits. For one with `s` fixed neighbors and quotient degrees `t_i`,
  degree and adjacent common-neighbor counts give `Σt_i=20−s`,
  `ΣC(t_i,2)=9−s`, and `Σt_i²=38−3s`. Cauchy–Schwarz would require
  `(20−s)²≤(c−1)(38−3s)`, but at `c=10` the gap polynomial
  `s²−13s+58` has discriminant `−63`. This part excludes all five counts
  for every `srg(45,22,10,11)`, without a Ramsey hypothesis.

Consequently the automorphism group has order coprime to three. In particular,
no Ramsey-good conference graph on 45 vertices is vertex-transitive.

```bash
./.venv/bin/python -m unittest r55.tests.test_order3_srg_frontier -v
./.venv/bin/python r55/src/order3_srg_frontier.py \
  --output r55/data/order3_srg_frontier.json
./.venv/bin/python r55/src/check_order3_srg_frontier.py \
  r55/data/order3_srg_frontier.json
```

Terminal disposition:
`ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG`. The result is exact,
solver-free, and independently replayed. It excludes symmetry only inside the
strongly-regular construction lane; non-strongly-regular hypothetical
Ramsey(5,5,45) graphs remain open, and no Ramsey-number bound changes.

## Odd-prime conference frontier (exactly closed, 2026-08-23)

The symmetry argument now closes every odd prime, not only three:

> **The automorphism group of a Ramsey-good `srg(45,22,10,11)` is a
> 2-group.**

The new proof excludes orders five, seven, eleven, and every prime at least
thirteen; the independently checked order-three theorem supplies the remaining
odd prime.

- **Order five.** A quotient-diagonal moment excludes fixed counts
  `10,15,…,40` for every `srg(45,22,10,11)`, without a Ramsey hypothesis.
  Ramsey trace then leaves `f=0,5` with every moving orbit a `C5`. For `f=0`,
  the required `9×9` integral Seidel quotient has three possible first-row
  types; a complete 56,755-prefix enumeration finds none. For `f=5`, an
  eight-block fixed-incidence design and `R(3,4)=R(4,3)=9` give a
  solver-free `K5/I5` bridge.
- **Order seven.** Cyclotomic parity leaves `f=3,17,31`. Handshake, a
  support/complement packing bound, and an exact fixed-block Rayleigh
  contradiction eliminate the three cases.
- **Order eleven.** The `f=23` support classes force an independent set. For
  `f=1`, quotient algebra is unique up to one swap; all
  `10²·C(11,6)=46,200` cyclic two-orbit blocks contain a `K4` or an `I5`.
- **Larger primes.** Orders `13,17,19` fail by fixed-graph handshake, while
  `p≥23` is incompatible with the required even number of moving cycles.

By Cauchy's group theorem no odd prime divides the automorphism-group order.
Involutions and nontrivial 2-groups are the remaining symmetry frontier.

```bash
./.venv/bin/python -m unittest r55.tests.test_odd_prime_srg_frontier -v
./.venv/bin/python r55/src/odd_prime_srg_frontier.py \
  --output r55/data/odd_prime_srg_frontier.json
./.venv/bin/python r55/src/check_odd_prime_srg_frontier.py \
  r55/data/odd_prime_srg_frontier.json
```

Independent published-data corroboration also parses Maksimović's 288
`S3`-invariant records, verifies their stated `f=9` order-three action, and
finds both a `K5` and an `I5` in every graph:

```bash
./.venv/bin/python r55/src/validate_published_s3.py --verify-existing
```

Terminal dispositions:
`AUTOMORPHISM_GROUP_IS_A_2_GROUP_FOR_RAMSEY_SRG` and
`PUBLISHED S3 DATA VERIFIED`. The result remains confined to the
strongly-regular construction lane. It says nothing about a
non-strongly-regular Ramsey(5,5,45) graph and changes no Ramsey-number bound.