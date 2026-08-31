# NA-Compiler pre-statement — certified-optimal small-instance transport baselines

Written 2026-08-29 after first-hand reads (same day) of:
- arXiv:2604.25478v2 (RSQASM evaluation frame; read via arXiv API + full PDF text).
- arXiv:2505.22715v1 (Routing-aware placement, QMAP implementation; full PDF text).
- MQT QMAP source `main`, mqt.qmap 3.9.0 (pinned wheel, commit-less; exact files listed
  in `campaigns/.../provenance`), docs `na_zoned_compiler.md`, and the shipped
  evaluator `eval/na/zoned/eval_ids_relaxed_routing.py` (all read raw from GitHub main
  2026-08-29).

---

## 1. The exact constraint model (pinned from code, not from paper prose)

Architecture (one storage zone SLM + one entanglement zone, AOD shuttling):

- **Sites.** Storage SLM is an r×c grid of traps; entanglement zone is a pair of SLMs
  with row-aligned paired traps (site_separation 12 µm × 10 µm in MQT's shipped
  `square_architecture.json`; storage separation 4 µm × 4 µm). Distances are Euclidean
  between exact trap coordinates. We scale row/column counts down for ≤10-qubit work,
  keeping the separations (hence physical timings) identical to MQT's shipped arch.
- **Initial placement.** QMAP's `VertexMatchingPlacer::makeInitialPlacement` fills the
  storage SLM row-major from row 0, column 0 (`reverseInitialPlacement` optional). Both
  QMAP and our ILP start from the same row-major initial placement unless stated
  otherwise. This removes initial-mapping freedom as a confounder.
- **Two-qubit gate layers.** ASAP scheduling exactly as documented in arXiv:2505.22715
  §III-A and implemented by QMAP's `ASAPScheduler`: gates enter the earliest layer in
  which both operands are still free; layers alternate 1Q/2Q. This layering is a
  deterministic function of the circuit's gate order and is computed identically by us
  and by QMAP, so all comparisons are placement+routing comparisons, never
  scheduling comparisons. 1Q gates do not move atoms and are ignored in transport cost.
- **AOD rearrangement (one "move batch").** A batch moves a set of atoms from sites to
  sites. Strict legality (the QMAP `IndependentSetRouter::isCompatibleMovement`
  relation, read in source) is pairwise over moves i→(xs_i,ys_i), i→(xt_i,yt_i):
  1. row preservation: (xs_i == xs_j) iff (xt_i == xt_j);
  2. column preservation: (ys_i == ys_j) iff (yt_i == yt_j);
  3. x non-crossing/order: (xs_i < xs_j) iff (xt_i < xt_j);
  4. y non-crossing/order: (ys_i < ys_j) iff (yt_i < yt_j).
  A batch is legal iff all its move pairs are compatible. (QMAP additionally has a
  RELAXED method allowing same-row order changes at a merge cost; we compare against
  STRICT, the conservative default, and document RELAXED as out of scope for gate A.)
- **Ghost-spots.** The paper (arXiv:2505.22715 §II-A) states the ghost-spot constraint
  (pick-up/drop-off affects whole AOD rows/columns). QMAP's shipped strict router does
  NOT implement it (grep over `IndependentSetRouter` found no ghost-spot check — this is
  a documented observation, not an inference). Primary gate A therefore uses the
  QMAP-comparable constraint set (the four conditions above), with a `--ghost` flag
  adding the stricter pick-up/drop-off row/column disjointness as a secondary run.
- **Minimum separation.** No additional pairwise minimum-distance constraint exists in
  QMAP's strict router beyond legality of start and target sites (sites are discrete
  traps). Separation is enforced by zone geometry, not by the router. Pinned
  accordingly.
- **Move-batch parallelism.** All moves inside one batch happen in parallel; batches of
  one layer transition execute sequentially; batches across different layer transitions
  never interleave (QMAP code-generator emits per-transition groups sequentially).

## 2. The one duration model (primary), and the runner-up

Both QMAP's output and our optimum are scored under **one** primary model, QMAP's own
official evaluator arithmetic (`eval/na/zoned/eval_ids_relaxed_routing.py`, read today):

- trap transfer (each `load` or `store` event): **+15 µs** (from QMAP
  `Architecture.hpp`: `timeAtomTransfer = 15`).
- one move batch with max atom travel distance d (µm): jerk-limited profile
  - t(d) = 2·(4d/j)^(1/3) µs for d ≤ 110 µm, with 2·(4·110/j)^(1/3) = 200 µs, i.e.
    j = 32·110/200³ µm/µs³;
  - t(d) = 200 + (d − 110)/1.1 µs for d > 110 µm.
- batch duration = t(d_max over atoms in the batch); transition duration = sum over its
  batches; run "transport cost" = sum over all transitions of batch durations +
  15 µs × (number of load/store events). Gate and 1Q durations (cz 0.36 µs,
  single-qubit 52 µs in the shipped arch) are identical for every solution of a given
  instance, so the gap metric is transport-only; we still record total schedule time.
- Secondary descriptor (reported, not gated): RSQASM frame's linear model,
  t = 2·t_trans + D/υ with t_trans = 20 µs, υ = 0.55 µm/µs (arXiv:2604.25478 §III-E),
  and QMAP paper's t = (d/2750)^(1/2) (arXiv:2505.22715 §VI-B). Different models change
  gap *magnitudes*, not gap *signs*; this sensitivity is part of the result.

## 3. Gate A — ILP certified optima vs MQT QMAP (exact pass/fail)

**Instances (all generators committed in `src/na_compiler/circuits.py`; provenance is generated
code + pinned dependency versions, no external downloads):**
- QFT family: Qiskit `QFT(n)` for n ∈ {4, 6, 8, 10}, decomposed at
  optimization_level=0; instance = its ordered 2-qubit gate sequence (cz pairs with
  included swaps), qiskit 2.5.2 pinned.
- 3-regular family: `networkx.random_regular_graph(3, n, seed)` for n ∈ {4, 6, 8, 10},
  seeds {0, 1, 2}, edges in sorted order; ASAP layering = greedy edge coloring in edge
  order (identical semantics to QMAP's ASAP scheduler).
- Architecture per instance: scaled copy of MQT's `square_architecture.json`
  (separations kept; storage c ≥ 2n, r ≥ 2, entanglement ≥ widest layer).

**A1 (primary, always run): routing-stage certificate.** Take QMAP's compiled
placement trajectory (parsed from its `.naviz` output). For each transition, compute
with an exact solver the minimum number of legal move batches (this is exactly
conflict-graph coloring: min batches = chromatic number χ of the strict conflict graph)
and, separately, the minimum transport duration over all legal partitions into that
many batches with free internal batch ordering. QMAP's router is a greedy
maximal-independent-set heuristic (source-verified), so its batch count ≥ χ.
- Certified output: χ per transition via exhaustive solver proof — Z3 returns
  `unsat` for a (χ−1)-coloring and `sat` for χ (both logs archived); duration optimum
  via HiGHS MIP with archived dual bound (log `gap 0.0`, `status Optimal`).
- **Gate A pass:** for every instance in the table, a certificate file exists with
  both solver logs and gap = (QMAP transport cost − certified optimum)/certified
  optimum ≥ 0, reproducible by `campaign_runner.py --verify-certificates`. Gap = 0
  corroborates QMAP's router on that transition; gap > 0 is a certified refutation
  instance for the router heuristic.
- **Gate A fail (fiction):** any reported optimum without both solver logs; any gap
  computed across different duration models; any instance where the parser silently
  dropped events (parser asserts balanced load/store-to-move correspondence).

**A2 (secondary, small): joint placement+rauting optimum.** Additionally give the
solver both QMAP's initial placement AND per-layer placement freedom (storage sites ≥
entanglement-capacity only), minimize transport duration; run only for n ≤ 6 and ≤ 3
2Q-layers, Z3-answered, always labeled A2 and never pooled with A1 rows. Comparison
meaning: end-to-end optimality gap of QMAP's full pipeline on tiny inputs.

## 4. Gate B — exhaustive hardest permutation instances (exact pass/fail)

- State: n atoms (n ≤ 8) injected onto distinct sites of an S×S grid, S ≤ 4; target =
  specified permutation of occupied sites. Moves = legal batches per §1 (four pairwise
  conditions on integer coordinates; no zones — single-zone monolithic grid so hardness
  is model-pure). Objective: minimum number of batches start→target.
- Algorithm: IDA*/BFS over batch-set branches with memohash canonicalization
  (symmetry-reduced states); per-instance optimality by exhaustive exhaustion of depth
  d−1 (unsat = no solution with d−1 batches exists; the search tree IS the certificate,
  archived as node counts + visited states).
- Enumeration protocol: for n ≤ 6 on grids ≤ 3×3 the target space is exhaustively
  enumerated (all P(S²,n) reachable states from a fixed start) and the *diameter-hard*
  instances (argmax min-batches) are published. For n ∈ {7,8} on 4×4, exhaustive
  enumeration of P(16,n) > 4·10^7 states exceeds in-warehouse budget; protocol = BFS
  from identity to full depth where the frontier stabilizes (documented per-n frontier
  caps), then certified hardest among *reached* states only, explicitly labeled
  "hardest-found", never "hardest".
- **Gate B pass:** per-instance certificates (node/visit counts) for the exhaustive n ≤
  6 range, published hardest-found table for n ∈ {7,8} with frontier sizes. **Fail:**
  any "hardest" claim without the exhaustion log; extrapolating n ≤ 6 diameters to n ∈
  {7,8}.

## 5. "Certified" — precise meaning

A number is CERTIFIED iff (a) an exact solver (Z3 SAT/SMT complete decision, or HiGHS
MIP with `Status: Optimal` + nonzero best dual bound logged) proves optimality of that
quantity under the §1 constraint set and §2 duration model; (b) the solver's decision
log for the tight instance (model with bound k−1 unsat) is archived verbatim; (c) the
instance generator and architecture JSON are content-hashed into the run manifest.
Everything else is labeled NUMERICAL, [DERIVED], or [INFERENCE]. Labels never improve.
No sampled shortcut may be promoted to a universal claim; "hardest" vs "hardest-found"
are distinct labels with distinct proof obligations (§4).

## 6. Non-goals

No hardware experiments, no fidelity-model unification beyond §2, no AOD RELAXED-mode
certificates (out of scope, documented), no cross-tool comparison beyond MQT QMAP in
gate A (RSQASM tools are comparison-context only).

---

## Revision 1 — gate-B complete rows beyond the 30k cap (written 2026-08-30, BEFORE any run)

**Trigger.** The frozen gate-B campaigns (`20260830T002708Z_6907d87f_deb7e0f2bbe9`) hit
the 30k-state BFS cap on 4×4 at n ∈ {5,6}, so their `hardest` entries are
"hardest-found" (frontier-capped), never "hardest". Revision 1 closes rows for at
least two of the capped n ∈ {5,6} 4×4 instances by replacing the state-cap BFS with
IDA*/DFS-with-admissible-pruning over the SAME frozen model (§1 legality relation
unchanged, §4 semantics unchanged, §2/§ evaluator duration model unchanged — arXiv
2505.22715 §VI-B/QMAP evaluator arithmetic in `src/na_compiler/config.py`). **Every
frozen artifact is left intact; this revision adds one new module (`src/idastar.py`)
and new campaigns. The original BFS code, its certificates, and all pre-Revision-1
campaigns are untouched.**

**R1.1 Target quantity (unchanged from §4).** For each instance (n, 4×4) with
identity start `state = (0, 1, ..., n−1)` and the FIXED grid, the depth-minimal
batch count

  d* = min{ k : identity → permutation π reachable in exactly k legal move
            batches on the fixed grid }

where "legal batch" is unchanged from §4: any set of singleton moves
atom → currently-free target site with pairwise-distinct targets, all pairwise
compatible under QMAP's strict 4-condition relation (§1, filters 1–4, no ghost
spots, no minimum-separation constraint). The move-space and the reachability
relation are IDENTICAL to §4; only the search algorithm differs.

**R1.2 Heuristic (pinning before the run).** For a state s with occupied sites
O(s), define

  h(s) = ceil( max_{atoms a in s} d1(site(a), goal-site(a)) / D_max )

where d1 is the Manhattan distance in grid steps between an atom's current site
and its goal site on the 4×4 grid, and D_max = maximum number of grid steps any
SINGLE atom can traverse in ONE legal batch on the 4×4 grid. D_max is computed
exactly by table: for each grid step-length pair (|Δx|, |Δy|) ≤ (3,3),
D_max(|Δx|, |Δy|) = max over all legal batch placements of the reachable grid-step
distance in one batch for an atom moving that (Δx, Δy). Because a batch moves each
atom by at most D_max steps in the max-manhattan sense, and every batch reduces
each atom's remaining Manhattan distance-to-goal by at most that bound, the number
of batches needed is ≥ ⌈max_a d1(a) / D_max⌉. Monotonicity: every edge
(s, batch, s′) changes h by at most 1 downward (a single batch covers at most
D_max of any one atom's remaining distance); h(goal)=0. Hence h is both
admissible and consistent for the IDA*/A* use below.

**R1.3 Admissibility argument (explicit).** Every legal batch is a set of disjoint
permutation moves on grid sites; each atom travels by (Δx, Δy) with
|Δx| ≤ 3, |Δy| ≤ 3 and, critically, each atom's move strictly reduces d1(a →
goal(a)) by at most D_max in one batch — formal LEGALITY LEMMA: for the true
optimal solution σ = (b_1, ..., b_d*), the first batch b_1 maps s → s′ with
h(s′) ≥ h(s) − 1; induction over batches gives h(s) − h(s′) ≤ d*, i.e.
d* ≥ h(s) for all s. This is admissible because NO batch can reduce any atom's
remaining Manhattan distance by more than D_max in one batch (a batch is a
single simultaneous displacement per atom), and the max over atoms is a lower
bound on what any ONE batch must accomplish for the slowest-to-close atom.
Consistency follows from the same per-batch decrement bound (h(s) ≤ h(s′) + 1).

**R1.4 What a "complete row" certifies.** A row for instance (n, 4×4) is COMPLETE
iff the algorithm (i) exhibits a witness: an explicit sequence of d* ≥ h legal
batches, each batch legal under §4 (pairwise strict-legality + free distinct
targets at application time), verified programmatically against
`na_compiler.conflict.compatible` + the §4 free-target rule, and (ii) proves
d* − 1 is UNSAT by exhaustive recursive search over all legal-batch choices up to
depth d* − 1, with pruning applied ONLY by the R1.2 admissible-consistent
heuristic (prune f = g + h > bound) and, where used, a transposition table that
is never a heuristic shortcut (stores exact DFS exhausted-depth bounds only).
Any such closed row hence certifies d* = optimal #batches under the §4 batch
semantics AND, immediately, optimal transport cost d* → the cost model in §2:
the duration is monotone non-decreasing in the number of batches under the
pinned evaluator (each batch costs ≥ 0, and a witness batch set executed in a
d*-length sequence with per-batch max-distances is the duration-optimal
realization of the min-batch schedule because an extra batch can never reduce
total duration). Per §5, all via exhaustive search + explicit optimality, no
solver dependency, no sampling, no heuristic acceptance.

**R1.5 Honesty: what remains open after Revision 1.** Any instance not closed by
R1.4 remains labeled "hardest-found" (open). Frontier-rate and node-count data
are REPORTED, not CERTIFIED. No extrapolation from n=5 to n=6,7,8. The original
§4 "hardest" vs "hardest-found" distinction is unchanged: a closed d* row
certifies ONLY that this state's min-batch count is d*, never that the state is
the diameter-instancemost-hard. Complete-vs-partial labeling in the revision
campaigns uses the code `complete=True` only if R1.4(i)+(ii) both hold.

**R1.6 Runtime/instrumentation.** Runs are nice -n 10 single-core; wall clock is
self-capped at 2 h total across all instances in this revision; every run
reports nodes expanded, time per depth level, and witness batches for
reproducibility. If the 2 h cap fires, closed rows are reported as complete,
the rest as frontier-capped with the exact node counts, per pre_statement §4
honesty. Revision 1 does NOT modify §1, §2, §3, §4, §5, §6 in any respect that
affects model comparability; frozen campaigns keep their original labels.

### Revision 1 cost-search erratum (written 2026-08-30 AFTER batch-count BFS, BEFORE cost runs)

R1.4's sentence that an extra batch can never reduce duration was over-broad:
nonnegative edge costs alone do not compare different paths whose per-batch
distances differ. It is retracted for the duration column. Batch-count
certificates already run are unaffected: full BFS exhausted P(16,5) and
P(16,6), so their optimal-batch/diameter proofs do not depend on that sentence.

For the duration column, the cost of a legal batch B is pinned, before any cost
run, to exactly the §2 evaluator arithmetic applied to that batch:

  w(B) = t(max_{(a→u)∈B} EuclideanDistance(site(a), u) in µm)
         + 2·15 µs·|B|.

The 4×4 sites use the frozen 4 µm × 4 µm storage-site spacing. Thus the second
term is one `load` plus one `store` for every atom moved by the batch. A path's
cost is sum_B w(B). `opt_cost_us` is the GLOBAL minimum of that additive cost
over all legal paths identity→target (batch count unconstrained). It is proved
by A* over the exact same finite state graph; the min-batch scalar is proved
separately by the already-pinned exhaustive BFS.

Cost-search heuristic:

  h_cost(s) = 30 µs·|M(s)|
              + max_{a∈M(s)} t(EuclideanDistance(site_s(a), goal(a))),

where M(s) is the set of misplaced atoms. Admissibility: every misplaced atom
must move in at least one future batch, so remaining load/store cost is at least
30·|M|. For any atom a, the sum of its future segment lengths is at least its
straight-line remaining distance (triangle inequality); the pinned t(d) is
increasing and subadditive on this grid's cubic branch, and each batch duration
dominates t of a's segment in that batch, so remaining batch-duration sum is at
least t(straight-line distance). Taking the maximum atom avoids double-counting
shared batch duration, while adding the independent load/store lower bound is
valid. Consistency follows from the same argument on one edge: an edge moving k
atoms pays 30k plus t(d_max), which dominates the drop in both heuristic terms.
Therefore A* finalizing the target certifies globally optimal `opt_cost_us`;
no assumption about its batch count is used.
