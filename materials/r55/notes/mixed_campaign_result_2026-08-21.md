# Engström mixed identity at $n=45$: terminal campaign record

**Date:** 2026-08-21
**Campaign:** `higher_order_identity_positive_deficiency_mixed`
**Disposition:** `MIXED_NO_CUT_IN_FROZEN_BASIS` (accepted routes: 0)
**Artifact:** [`data/engstrom_identity.json`](../data/engstrom_identity.json) (schema 3)

## 1. Result

Let $\mathcal{C}_{45}$ be the finite state relaxation of a hypothetical
Ramsey$(5,5,45)$ graph used throughout this repository: $3{,}215$ states
$(d,a,b)$ with $d\in\{20,\dots,24\}$, each carrying exact integer intervals for
the m=2 excess balance, the m=3 row $g$, the m=4 row $h$, and — newly added by
this campaign — the Engström mixed degree-4 row $F$.

**Terminal result.** Over the full eight-coefficient basis
$(\alpha,\beta,\gamma,\delta,\varepsilon,\zeta,\eta,\theta)$ spanned by
$\{\text{balance}, g, h, F\}$, every instantiable route of the frozen registry
is rejected by an **exact rational primal witness**:

| route | exact upper bound $45\alpha$ | exact witness value | status |
|---|---|---|---|
| `total_deficiency` | $360$ | $360$ | `REJECTED_BY_EXACT_WITNESS` |
| `degree20_count` | $14310/349$ | $14310/349$ | `REJECTED_BY_EXACT_WITNESS` |
| `deficiency_ge8_count` | $45$ | $45$ | `REJECTED_BY_EXACT_WITNESS` |
| `required_local_family` | — | — | `UNAVAILABLE_NO_COVER_CERTIFICATE` |

This fulfils **stop condition 2** of `bench/spec.json:next_research_campaign`:
the recorded `first_exact_basis` (MR97 $m=2,3,4$ plus the Engström mixed
identity) is now exhausted at $n=45$. No claim about $R(5,5)$ is made or implied.

## 2. The sharp quantitative finding (the actual contribution)

The negative disposition is not the interesting part; the *mechanism* is. Three
exact measurements, all reproducible from the artifact:

**(a) The mixed row is geometrically new.** Over the $3{,}215$ states, the
endpoint matrix rank is exactly

$$
\operatorname{rank}[\,\mathbf{1}, \text{bal}, g_{lo}, g_{hi}, h_{lo}, h_{hi}\,]=6,
\qquad
\operatorname{rank}[\,\cdots, f_{lo}\,]=7,
\qquad
\operatorname{rank}[\,\cdots, f_{lo}, f_{hi}\,]=8 .
$$

So $f_{lo},f_{hi}$ are **not** affine combinations of the frozen endpoints:
the naive "coordinate collapse under projection" hypothesis of the plan's
OPEN RISK section is **refuted**. The row adds two genuinely new dimensions.

**(b) Yet it adds exactly zero cut strength.** For all three routes the exact
dual optimum *and* the exact primal optimum are bit-identical with and without
the two $F$ columns/rows:

| route | dual (m4 basis) | dual (mixed basis) | primal (m4) | primal (mixed) |
|---|---|---|---|---|
| `total_deficiency` | $360$ | $360$ | $360$ | $360$ |
| `degree20_count` | $14310/349$ | $14310/349$ | $14310/349$ | $14310/349$ |
| `deficiency_ge8_count` | $45$ | $45$ | $45$ | $45$ |

Every optimal certificate carries $\eta=\theta=0$. This is the same phenomenon
the $m=4$ tier produced ($\varepsilon=\zeta=0$), now measured for a
**non-separable** row whose support is disjoint from the frozen basis — i.e. the
zero-strength outcome is not explained by shared support or by rank deficiency.

**(c) The reason: the interval abstraction, not the algebra.** The mechanism is
that the new aggregate constraints are *strictly inactive* at the optimal
vertices. At the exact optimal witnesses:

| route | $\sum_s w_s f_{lo}(s)$ | $\sum_s w_s f_{hi}(s)$ |
|---|---|---|
| `total_deficiency` | $-18522745/2$ | $68122015/2$ |
| `degree20_count` | $-2245549320/349$ | $9446922090/349$ |
| `deficiency_ge8_count` | $-11622240$ | $3360420$ |

Both slacks are enormous because the per-state $F$ windows are enormous:

| row | min width | median width | mean width | max width |
|---|---|---|---|---|
| $g$ | $0$ | $930$ | $986$ | $1{,}458$ |
| $h$ | $0$ | $41{,}128$ | $41{,}808$ | $72{,}078$ |
| $F$ | $864$ | $802{,}396$ | $764{,}953$ | $1{,}318{,}306$ |

The median $F$ window is **wider than the median $|F|$ magnitude**
($802{,}396$ vs $664{,}870$; ratio $1.207$). The mixed row's information is
destroyed by the *per-coordinate interval relaxation*, not by the identity and
not by the projection.

**(d) The discarded structure is measurable.** On the complete catalog class
$(22,114)$ ($133$ graphs), the $x$-side contributions of the $h$ row and the $F$
row have Pearson correlation

$$\rho(h_x, F_x) = 0.9464 ,$$

with $h_x\in[1686,2056]$ and $F_x\in[234048,239156]$. The box relaxation
replaces this near-degenerate 2-dimensional point cloud with its
$370\times5108$ bounding box, i.e. it discards a 95 %-correlated joint
constraint. **This is the precise, measured reason the frozen basis is
exhausted, and it identifies the next campaign** (Section 5).

## 3. What was built

| file | role |
|---|---|
[`src/mixed_subgraph_identities.py`](../src/mixed_subgraph_identities.py) | Task 1 kernel: eleven raw induced 4-vertex counters, corrected complement-pair map, `mixed_vertex_row`, `mixed_residual` |
[`src/mixed_deficiency_cone.py`](../src/mixed_deficiency_cone.py) | Task 2 producer: fifteen-aggregate catalog streaming, envelope-LP/degree-histogram outer windows, eleven-field `State`, v3 artifact writer + validator |
[`src/search_mixed_cuts.py`](../src/search_mixed_cuts.py) | Task 3 search: exact eight-coefficient dual certificates, exact primal witnesses, frozen route registry, atomic artifact update |
[`src/verify_mixed_independent.py`](../src/verify_mixed_independent.py) | Independent (stdlib-only, zero producer imports) re-verification of the search layer and of all $3{,}215$ $F$ intervals |
[`tests/test_mixed_subgraph_identities.py`](../tests/test_mixed_subgraph_identities.py) | Kernel contracts incl. exhaustive $n\le6$ (33,867 graphs) and the six corrected-math mutation pins |
[`tests/test_mixed_deficiency_cone.py`](../tests/test_mixed_deficiency_cone.py) | Cone contracts (11 tests) incl. frozen-projection equality and fixture-pair containment |
[`tests/test_search_mixed_cuts.py`](../tests/test_search_mixed_cuts.py) | Search contracts (13 tests) incl. both new sign gates, least-$\alpha$ canonicality, and the $F$-wiring pin |

## 4. Verification evidence (all observed, this session)

1. **Producer run.** `mixed_deficiency_cone.py` streamed $8{,}500{,}211$ graphs
   in $4348.2$ s (budget $7200$ s), reporting
   `MIXED REPLAY: 656/656 residual zero`,
   `MIXED N49 REPLAY: extremal=2 combos=4 values=[0, 144, 288, 432] zero_corners=1`,
   `MIXED ENVELOPES: catalog violations=0 outer-window audits=ALL`,
   `MIXED STATES: 3215 states (asserted equal to frozen m4 projection)`.
2. **Focused suites.** `test_mixed_deficiency_cone.py` 11/11 OK (192 s);
   `test_search_mixed_cuts.py` 13/13 OK.
3. **Search run.** `MIXED CUT SEARCH COMPLETE status=MIXED_NO_CUT_IN_FROZEN_BASIS accepted=0`.
4. **Independent verifier.** Re-hashed 64 inputs, re-derived $8{,}500{,}211$
   catalog graph counts across 53 classes, recomputed all $3{,}215$ $F$
   intervals from the artifact's own windows with independently written affine
   arithmetic, re-derived every certificate/witness/status →
   `MIXED EVIDENCE VERIFIED: MIXED_NO_CUT_IN_FROZEN_BASIS`.
5. **Corruption battery.** Twelve targeted mutations (state $f_{lo}$,
   $\alpha$ lowered, $\alpha$ padded, primal weight, trust-root URL,
   $g$ endpoint vs frozen m4, histogram count, $n{=}49$ replay field, route
   status, deleted state, schema version, swapped $f$ pair, input SHA-256) —
   **all rejected**; untouched control passes.

### Remaining verification gap (stated, not hidden)

The plan's Task 4 disjoint checker is **not** complete. Specifically not yet
independently redone: (i) the $8.5$M-graph catalog re-sweep with a second,
independently written motif kernel; (ii) independent recomputation of the $g$
and $h$ endpoints from the R1–R5 relation system (currently anchored by exact
equality against the frozen `higher_identity_m4.json` projection); (iii)
independent re-derivation of the envelope-LP windows for the 77 missing stratum
classes (the m4 campaign scoped this identically). The disposition above is
therefore **producer-verified + independently re-verified at the search and
$F$-interval layer + cross-anchored to the frozen m4 artifact**, which is
strictly stronger than the m4 tier's own recorded evidence but short of the
full Task 4 tier.

## 5. Reproduction

```bash
cd /Users/jinleic/jinleic-workspace/math
# Task 1 + 2 + 3 contracts
./.venv/bin/python -m unittest discover -s r55/tests -p 'test_*mixed*.py' -v
# Task 2 producer (about 75 min; rewrites the artifact from scratch)
./.venv/bin/python r55/src/mixed_deficiency_cone.py \
  --output r55/data/engstrom_identity.json
# Task 3 search (seconds; refuses to run twice)
./.venv/bin/python r55/src/search_mixed_cuts.py \
  --analysis r55/data/engstrom_identity.json
# Independent verification (sub-second)
./.venv/bin/python r55/src/verify_mixed_independent.py \
  r55/data/engstrom_identity.json
```

## 6. Where the obstruction really is: window quality is NOT the bottleneck

After the disposition landed, four further exact experiments were run to decide
which follow-on campaign is worth executing. Three plausible candidates —
including the joint-hull design this note originally proposed — are **ruled out
by measurement**. Recorded here so nobody spends a campaign on them.

### 6.1 Where the width lives

Decomposing every state's $h$ and $F$ window width by side and by certification
method (summed over all $3{,}215$ states):

| bucket | summed width | share |
|---|---|---|
| catalog-certified classes | $35{,}649{,}683$ | $1.4\%$ |
| envelope-LP classes (77 of 130) | $2{,}558{,}092{,}854$ | $98.6\%$ |
| trivial fallback | $0$ | $0\%$ |

So the catalog-side correlation structure of §2(d) — real, $\rho=0.9464$ —
governs only $1.4\%$ of the total width. **The joint-hull relaxation over
catalog clouds cannot matter.** That kills candidate 1.

### 6.2 How loose the envelope really is

Calibrating the envelope LP against catalog ground truth on the five present
extremal classes (same class, same key, truth vs LP):

| class | key | catalog truth | envelope LP | looser by |
|---|---|---|---|---|
| $(22,114)$ | `i4` | $[67,99]$ | $[0,1285]$ | $40\times$ |
| $(22,114)$ | `k3k1` | $[227,364]$ | $[0,5341]$ | $39\times$ |
| $(24,132)$ | `k3k1` | $[528,528]$ | $[0,6153]$ | $6153\times$ |
| $(20,100)$ | `p3_k1` | $[746,746]$ | $[0,3188]$ | $3188\times$ |

Aggregate over the 25 measured (class, key) pairs: true width $836$ versus
envelope width $92{,}888$ — the envelope LP is **$111\times$ looser than
reality**. That looked like the dominant lever.

### 6.3 The decisive experiment: it still would not matter

Recomputing every $h$ and $F$ endpoint from the artifact's own window table
while shrinking all envelope-class windows toward their midpoints by a factor
$k$ (the $g$ row held at its frozen m3 values; at $k=1$ the recomputation
reproduces all $3{,}215$ recorded $h$ and $F$ intervals **exactly**, which is
itself an independent re-derivation of both rows):

| $k$ | `total_deficiency` (edge $\le315$) | `degree20_count` (edge $<1$) | `deficiency_ge8_count` (edge $<1$) |
|---|---|---|---|
| $1$ (as certified) | $360$ | $41.003$ | $45.000$ |
| $0.1$ | $360$ | $39.352$ | $38.632$ |
| $0.01$ | $360$ | $39.260$ | $38.387$ |
| $0$ (windows collapsed to points) | $360$ | $39.251$ | $38.363$ |

Strong duality holds throughout (primal $=$ dual at $k=1$ and $k=0$; at $k=0$
the primal still places $39.251$ of the $45$ vertex-weight on degree-20 states).

**Conclusion.** Even a *perfect* motif-window certification — every one of the
77 catalog-missing stratum classes pinned to an exact point — leaves the three
routes at $360$, $39.25$, $38.36$ against edges of $315$, $1$, $1$. The gap is a
factor of $\approx39$, not a few percent. This rules out candidate 2 (structural
tightening of missing-class windows, e.g. Zykov/Turán caps such as
$i_4(z)\le900$ at order 22) and candidate 3 (exact joint projection of the
envelope polytope, separately probed at only $1.5\times$ area tightening in the
$(h_z,F_z)$ plane). No data-quality campaign can close a $39\times$ gap.

*Caveat, stated precisely:* the $k<1$ rows are a sensitivity analysis, not a
theorem — the shrunken windows are centred at current midpoints and are not
themselves sound certifications. The rigorous content is the $k=1$ row. But the
robustness is overwhelming: three orders of magnitude of width reduction buys
$1.75$ units on a gap of $38$.

### 6.4 What the obstruction actually is

The relaxation has exactly five aggregate conditions ($\sum w=45$,
$\sum w\,\mathrm{bal}=0$, and the three row-sum sign pairs) over $3{,}215$
nonnegative weights. An optimal vertex therefore has support $\le5$, and five
linear conditions simply cannot separate the degree-20 states from the rest:
the LP is free to place $\approx39$–$41$ of the $45$ vertex-weight on them. The
per-vertex interval relaxation discards **all** vertex-to-vertex coupling — it
never records that neighbourhoods of adjacent vertices intersect.

Two cheap candidate repairs were tried and also fail: the degree/edge handshake
coupling $\sum_s w_s d_s = 2e$ combined with the averaged per-vertex inequality
$e \ge d_v + e(N_v) + e(D_v)$ is slack by $\approx200$ units at every state
(because $e\approx495$ while the right-hand side is $\approx250$), so it prunes
nothing.

### 6.5 The one direction the evidence supports

What is needed is **more constraints, not better data**: a relaxation that
couples pairs of vertices. Concretely, a level-2 (Sherali--Adams / Lasserre)
lift carrying pair variables $w_{s,s'}^{\text{adj}}$ and
$w_{s,s'}^{\text{non}}$ with marginal consistency to the present $w_s$, plus
the exact per-pair facts this repository can already certify:
for adjacent $u,v$ the common neighbourhood $N(u)\cap N(v)$ is $K_3$-free
(else $G$ has a $K_5$ with $u,v$), codegrees satisfy
$\sum_{uv\in E}|N(u)\cap N(v)|=3T$, and each MR97/Engström row can be applied
to the pair-restricted counts. That layer sees exactly the structure the
current relaxation is blind to, and it is the cheapest rigorous step toward the
recorded fallback tiers (VeriPB-certified gluing, then
$\mathrm{srg}(45,22,10,11)$ complementary SAT).

Falsification criterion for that campaign: if the level-2 lift's exact optimum
on `degree20_count` stays above $1$, the pair-coupled relaxation of this basis
is also closed, and the gluing/SAT tiers become the only remaining routes.



## 7. The level-2 pair lift was executed. It is closed.

Section 6.5 named a falsification criterion. It has been run, and the criterion
fired: **the pair-coupled relaxation of this basis is also closed.**

Implementation: `r55/src/pair_lift_mixed.py`, tests
`r55/tests/test_pair_lift_mixed.py` (18 tests). Same trust discipline as the
rest of the tier: HiGHS only *proposes* dual multipliers; every accepted number
is rationalised and re-verified in exact `Fraction` arithmetic by checking dual
feasibility $A_{eq}^{\top}y_{eq}+A_{ub}^{\top}y_{ub}\ge c$ componentwise, so any
multiplier vector that passes certifies a sound bound whether or not it is
optimal.

### 7.1 Verdict

| route | frozen exact | level-2 exact | delta | edge |
|---|---|---|---|---|
| `total_deficiency` | $360$ | $360$ | $0$ | $\le 315$ |
| `degree20_count` | $14310/349$ | $41$ | $-1/349$ | $<1$ |
| `deficiency_ge8_count` | $45$ | $45$ | $0$ | $<1$ |

Maximised over every integer edge count for which the programme is feasible,
which is exactly $e\in[454,536]$ (a mild narrowing of the trivial $[450,540]$
implied by $d\in\{20,\dots,24\}$). Two routes do not move at all; the third
improves by $1/349$. No route reaches its acceptance edge. Disposition:
`MIXED_PAIR_LIFT_CLOSED`. **No bound on $R(5,5)$ is claimed.**

The lift carried what §6.5 asked for: ordered degree-pair variables
$M[d,d']$ (adjacent) and $N[d,d']$ (non-adjacent), symmetric and nonnegative,
with marginals $\sum_{d'}M[d,d']=d\,n_d$ and $\sum_{d'}N[d,d']=(44-d)n_d$; the
per-degree-class neighbour-degree identity (below); the per-state window
$20d_s\le e+\varphi(s)\le 24 d_s$; and three codegree facts —
$\lambda\le R(3,5)-1=13$ on edges, $\lambda\le\min(d,d',d+d'-30)$ on non-edges
from $|D(u)\cap D(v)|=43-d-d'+\lambda\le R(5,3)-1=13$, and
$\lambda\ge\max(0,d+d'-43)$ — tied to the vertex layer by
$\sum_v e_x(v)=3T$ and $\sum_{u\not\sim v}\lambda=\sum_v\binom{d_v}{2}-3T$.

### 7.2 Why: `excess_balance` *is* the neighbour-degree identity

The sharper finding, and the reason the lift could not work. For any finite
simple graph, pure double counting gives, for every vertex $v$,

$$\sum_{u\in N(v)}\deg(u)=e+e(G[N(v)])-e(G[D(v)]),$$

since $e=d_v+e_x+e_D+e_{xz}$ and the left side equals $d_v+2e_x+e_{xz}$.
Summing over $v$ and using $\sum_v\sum_{u\in N(v)}\deg(u)=\sum_u\deg(u)^2$:

$$\sum_v \deg(v)^2 = 45e+\sum_v e(G[N(v)])-\sum_v e(G[D(v)]).$$

Verified exhaustively on all $2{,}131{,}018$ labelled graphs up to $n=7$ plus
random graphs at $n=12,20,45$, and independently re-verified by a second agent
over all graphs on $n=4,5,6$ ($209{,}184$ per-vertex checks, zero violations).

Adjoining it as a new aggregate row is the obvious move. **It is already
there.** In this repository the $z$ stratum is the *complement* of the
anti-neighbourhood, so the recorded $e_z=E_{45}[m]-b$ is
$e(\bar G[D(v)])$ and the true interior count is $e_y=\binom{m}{2}-e_z$. With
that substitution the state functional collapses identically:

$$2d^2-45d-2e_x+2e_y \;=\; \texttt{excess\_balance}(s)\qquad\text{for all }3215\text{ states.}$$

Exact rank over the states is $8$ for the frozen basis and stays $8$ when the
row is adjoined. So `excess_balance` — the campaign's $m=2$ row — *is* the
degree-squared/handshake identity in disguise. That retires the entire family of
handshake-, assortativity- and degree-correlation repairs at the aggregate
level, including the one this note itself proposed in §6.4. Only the
*per-degree-class* form is strictly finer than the aggregate, and the table in
§7.1 shows what that finer form buys: $1/349$.

### 7.3 Recorded trap, and a retraction

During this work I first substituted the recorded $e_z$ directly for
$e(G[D(v)])$. That produces the row

$$r(s)=2d^2-45d-2\big(E_{45}[d]-E_{45}[m]\big)+2a-2b,$$

which **is** linearly independent of the frozen basis (rank $8\to9$) and which
appears to cut hard: $360\to4905/14$, $14310/349\to225/11$, $45\to4905/112$,
all with matching exact dual certificates and exact primal witnesses, and
tightening further to $349/20/43$ under integrality of the degree distribution.
**Every one of those numbers is retracted.** The row is not a valid constraint:
it misreads the stratum convention. Exact certificates and matching witnesses
prove optimality *for the program you wrote*; they say nothing about whether the
program is sound.

Two guards are now in the repository so this cannot recur silently.
`test_invalid_ez_row_is_the_trap` pins both halves — the trap row is
independent of the basis *and* differs from `excess_balance` — and
`test_nonadjacent_bound_is_an_upper_bound` pins the sign of the non-adjacent
codegree bound, because reading $\lambda\le d+d'-30$ as a lower bound renders
every $e$ infeasible and so manufactures a false refutation of
$R(5,5)=45$. Both failure modes were hit in this session and both were caught
by the same rule: a claim is not a result until it is re-derived from the
definitions in exact arithmetic.

### 7.4 What remains

With the pair-coupled layer closed, the recorded residual tiers are the only
remaining routes, in the order already established: VeriPB-certified gluing
search, then $\mathrm{srg}(45,22,10,11)$ complementary SAT. Nothing in the
vertex- or pair-local relaxation of this basis can open the three routes; the
evidence now covers levels $1$ and $2$, and the $\approx39\times$ margin of
§6.3 explains why data quality cannot substitute for a genuinely non-local
argument.