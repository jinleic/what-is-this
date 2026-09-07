# Progress ledger

Newest first. One entry per session. Every claim links to its verification.

---

## 2026-09-06 — calibrated magic-cat search stopped without rank-five witness

After the AES XOR assay failed to show local headroom, the next admitted
candidate was `chi(cat_8) <= 5`, which would improve QPG's Clifford+T
simulation exponent. The
[registered public-SA assay](stabrank-witness/campaigns/20260906T025240Z_725b49e8_31de21aa6d69/)
used one chain/thread at a time and fixed settings/seeds. Exact calibration
passed for cat4/rank2, cat6/rank3 and cat8/rank6. Seven rank-five seeds
completed without a witness; seed 7 hit the 60-second cell limit.
**FROZEN-INCONCLUSIVE**, 73.105411 CPU s, 72.049315 wall s; no retries.

The verifier independently reconstructs affine-support/quadratic-phase
stabilizers over Q(i), checks the full QPG target and rejects corrupt controls.
The frozen cat8/rank6 witness replayed all 256 coordinates in a standalone
process after closure. Q(i) suffices because certificate vectors omit each
state's nonzero normalization; exact cat6 coefficients and full replay
adjudicate the normalization advisory.

Post-hoc exact normal-equation checks on the seven already returned rank-five
tuples match the numerical residuals. Best squared distance: **205/2016 to
one particular span**, not a global floor or a lower bound on stabilizer rank.
Rank A-for-assay becomes **C for the unchanged SA lane**. Residual variation
does not establish LLM headroom. No large search or external submission.

[`Cross-field evidence and dossiers`](../data/frontier-scout-20260906/).

## 2026-09-06 — Portfolio supervision: compiler provenance adjudicated

**CLAIM VERIFIED, no new compilation or batching search.** The frozen
headroom runner directly constructs QMAP with strict routing, A-star
placement, `log_level="error"` and `max_nodes=100000`, matching its
preregistration. The advisory about the generic adapter's 10,000,000-node
default concerns a different call path. Source/result hashes and all 517
frozen payload checksums match; the raw exact-Fraction clique replay again
returned PASS on 12 circuits, 242 transitions and 438 optimum batches
(0.026308 CPU s). The frozen snapshot was not modified.

The test remains descriptive: fixed classical rules and a null policy for
zero-denominator gap closure, no adapted AI rule or generalization claim.
The existing strict-placement batching search remains stopped.
[`Audit evidence`](../data/portfolio-supervisor-20260906-na.json).

## 2026-09-05 — proposed batching search stopped by a classical-baseline ceiling

The [preregistered headroom assay](na-compiler/campaigns/20260906T012803Z_3c409686_c44470e4c9b4/preregistration.json)
closed **FROZEN-NEGATIVE on its investment gate**: all 12 fixed circuits,
242 strict-placement transitions and 438 high-level batches were already optimal
under QMAP, DSATUR, largest-first and smallest-last. No LLM candidates were
generated; 0/12 circuits had the required 10% residual gap (6/12 needed).

The [results](na-compiler/campaigns/20260906T012803Z_3c409686_c44470e4c9b4/results/results.json)
include every circuit and transition, raw NAViz, physical/high-level accounting,
subset-DP witnesses and Z3 UNSAT@k−1 / SAT@k decisions. Independent stdlib-only
[replay](na-compiler/campaigns/20260906T012803Z_3c409686_c44470e4c9b4/verify_cliques.py)
reconstructed all raw traces using exact Fractions and supplied a matching clique
lower bound for each partition (242/242). Assay: 1.033274 CPU s; replay: 0.025658
CPU s. A tool-free Fable review accepted the narrow claim; its overbroad argument
about spurious transition splits was not relied on (see review-adjudication.json).

This rules out batch-count improvement **on these fixed transitions**, not on
different placements or relaxed routing. The two n=8 regular graphs are isomorphic.
No hardware/novelty/AI-superiority claim, no larger run, and no GB9 or autonomous
supervisor restart. Portfolio decisions and unresolved literature leads:
[`data/research-resume-20260905.json`](../data/research-resume-20260905.json).

## 2026-09-05 — FSS reporting branch stopped on reference-integrity failure

[Primary Table 1](https://arxiv.org/html/2603.19062v3#S3.T1) gives
N576 `0.4453 [0.4452,0.4455]` and N900 `0.4674 [0.4673,0.4676]`.
The FSS row in RESULTS and the 2026-08-30 completion entry below instead
use `0.4608` and `0.4702`; their N576 CI-miss and four-of-five-overlap
statements are unsupported. Main read the primary table after an
independent read-only provenance review.

The same review found stale RUNNING wording in the FSS landing state and
a resource-scope error: the preregistered estimate is about50 core-hours
for the five-size Gate-A 200k ladder; about221 belongs to the two
additional Gate-B sizes. Those are historical estimates, not new launch
budgets. [Integrity evidence](../data/fss-reference-integrity-20260905.json)
records the affected statements and source authority.

**FSS work stopped pending integrity adjudication.** No reporting repair,
frozen-payload change, modern status-file fabrication, new scientific
verdict, or FSS workload was performed. The existing ν label remains
INCONCLUSIVE-UNDER-BUDGET; the reporting errors must not support further
promotion. GB9 supervision and its resource guard are independent.

Quarantine applied the same day, without touching frozen payloads or the ν
verdict: the FSS RESULTS row and landing state now carry an explicit
integrity hold naming the correct Table 1 values, the ≈50 core-hour Gate-A
scope, and a refusal to approve the 200k ladder off the current published
comparisons. `latest_campaign_verdict` reads
`UNVERIFIED-LEGACY-NO-STATUS-FILE` rather than a fabricated terminal verdict:
the legacy layout writes no `status.json`, and the evidence that the run
finished is indirect. Correcting the published numbers, and any FSS work at
all, still waits on adjudication.

## 2026-09-05 — MSD stale rebuild handoff replaced by the official-source route

The [MSD README](msd/README.md) and [semantic state](msd/state.json)
now route official-author reproduction to the
[canonical QLOPS handoff](qlops/README.md#official-source-revision-7--frozen-readiness-smoke-only).
The historical source-frame floor diagnosis remains intact, but no new
surrogate rebuild is queued. R7's eight-gate smoke PASS is frozen
REHEARSAL launch readiness only; the author artifact is d3/d9, not a
reproduction of the paper's d7/c≈300 headline.

The [frozen R7 summary](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/results/summary.md)
and [terminal status](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/status.json)
ground the handoff. MSD's historical data, preregistrations, negative
controls and verdict strings are unchanged; the paper claim remains
UNTESTED by its surrogate. No new campaign or computation was run.

## 2026-09-05 — all pinned shadow families classified algebraically for every k

The [complete pinned-family comparison](shadows/README.md#analytic-closure--complete-pinned-family-comparison)
adds disjoint-pair factorization and exact Clifford layer identities to
the chain recurrence below. Within these families, only the pure chain
and lex CNOT mesh at k4 strictly improve on U_ct. Their full histograms
agree for every k. Uniform-letter meshes and CZ mesh tie U_ct; every
near-perfect matching loses for k≥3; U_ct then CZ loses for k≥2.
The linked table states all small-k ties.

- Independent symbolic reviews validated both remaining proof slices.
- Main's exact checks passed for 96 layer maps on complete GF(2) bases
  at k=1..12, all three pair transvections, and matching direct sums at
  k=1..8; combined execution was under one second.
- [Evidence](../data/shadows-pinned-family-analytic-20260905.json) records
  the reviewed source hashes and exercised checks. No token or matching-
  family enumeration, sampling, solver, or new campaign was run.
- This is a scoped all-k theorem, not an arbitrary-Clifford obstruction.
  Largest empirical k remains12, exhaustive matching search remains
  k≤8, and all frozen payloads and experimental verdicts are unchanged.
  The separate noisy-variance gate remains FROZEN-CERTIFIED.

## 2026-09-05 — pure cx_chain asymptotic question closed by an exact transfer proof

The [canonical derivation](shadows/README.md#analytic-closure--pure-cx_chain-for-every-k)
settles the pinned pure open CNOT chain without another finite-k campaign.
Writing `A_k=9^k*w_cx(k)`, the exact transfer identity gives
`A_1=3, A_2=17, A_(k+2)=3*A_(k+1)+8*A_k`.
An integer forcing recurrence proves the complete comparison with U_ct:
**k=4 is the sole strict win; k=1,2,3,5 tie; every k≥6 loses.**
The weight ratio tends to zero at exponential base `(3+sqrt(41))/10<1`;
the reciprocal shadow-norm disadvantage grows.

- Two independent hand reviews checked the map/transfer construction and
  the all-k induction. The argument is not a fitted recurrence.
- Main's exact four-state polynomial check matched all nine stored k=4..12
  weight rows and all four complete k=9..12 size histograms. No token
  enumeration, sampling, solver, or new campaign was run.
- [Verification evidence](../data/shadows-cx-chain-analytic-20260905.json)
  records exact differences, matrix identities and input byte hashes.
  The legacy extension's short manifest identifier hashes its prereg/source
  inputs, not the output JSON; no raw-output ledger is inferred from it.
- Scope is the pinned pure open chain with uniform full-support X/Y/Z
  tokens. Largest empirical k remains12; family(b), inserted identities,
  rings and noise receive no new claim. Frozen campaigns and previous
  verdicts, including noisy-variance FROZEN-CERTIFIED, remain unchanged.

## 2026-09-04 — QLOPS official source pinned; R7 smoke PASS closed REHEARSAL; 71m full run unstarted

The stale no-official-source premise is retired. The
[`FujitsuResearch/Zero-level_CCZ_Distillation`](https://github.com/FujitsuResearch/Zero-level_CCZ_Distillation)
repository is pinned at commit
`1b59e223590492e224bd8623a4e0bcba59029e01`, tree
`9a5d89401fd40554635fb0e9d0b4da818670c2bc`, and codeload tar sha256
`d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`.
Its executable artifact is d=3/d=9; the paper's headline grown result is
d=7. These are never equated: paper d=7 and the physical c≈300 claim
remain **NOT-REPRODUCED**.

- **Instrument history:** R5
  `qlops/campaigns/20260904T052209Z_2a0a2723_203c9b1d1900/` is
  frozen/closed **SUPERSEDED** because its ungrown-p=0 shipped oracle was
  unsupported; it has no accepted result. R6
  `qlops/campaigns/20260904T054022Z_5f040b2a_1f64afb2c983/` is
  frozen/closed **CRASHED** after a results-directory orchestration bug;
  the regression was fixed before R7, and R6 has no accepted result.
  Neither predecessor is evidence.
- **Canonical readiness smoke:** campaign
  `qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/`, gate
  `zero-level-author-repro-r7`, preregistration sha256
  `99b3985b261dc3de561917736ea25c1e552283a67b2353df50f137a954bb0668`,
  is frozen and closed **REHEARSAL**. The smoke verdict is **PASS**, all
  eight gates: builder-only ungrown p=0 recovery invariants, zero-error
  DEM, zero detector/observable bits, and fixed-schedule replay; plus real
  shipped grown p=0.001 rebuilt equality, fault-id invariants,
  accepted-only scalar/batch bitwise agreement, and replay.
- Canonical evidence:
  [summary](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/results/summary.md),
  [smoke JSON](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/results/smoke.json),
  [manifest](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/manifest.json),
  [status](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/status.json),
  and [sha256 ledger](qlops/campaigns/20260904T054537Z_818bca48_18f5ac3b849c/sha256s.txt).
  `sha256sum -c` passed all 31 listed payloads. This smoke has
  launch-readiness weight only, not scientific-result weight.
- The full author-count schedule remains pending and **unstarted**:
  71,000,000 shots across 12 cells. After GB9 releases its 26 threads,
  follow the canonical exact two-command handoff in
  [`qlops/README.md`](qlops/README.md#official-source-revision-7--frozen-readiness-smoke-only).
  It mints a fresh `zero-level-author-repro-r7` campaign, substitutes its
  emitted run ID, consumes the immutable hash-verified R7 source tar, sets
  the OMP/OpenBLAS/MKL/NumExpr/vecLib thread caps to one, and launches via
  `/usr/bin/nice -n 10`; it does not depend on the ephemeral `/tmp` tar.
  The R7 environment is Python 3.14.3, stim 1.16.0, pymatching 2.4.0,
  and numpy 2.5.2.
- Existing conclusions are unchanged: Gate A remains **REPRODUCED**
  (81 checks, zero hard failures), and Gate B's magic-state axis still
  falsifies QLOPS comparability. No part of this readiness smoke promotes
  the zero-level physical result.

## 2026-09-04 — shadows noisy-variance gate closed FROZEN-CERTIFIED: exact 132-cell no-flip regression (gate noisy-variance-r7)

The channel-resolved real-versus-unitary variance question named as the
live next step on 2026-09-03 is answered — by exact confirmation of known
formulas, not discovery.

- Campaign `shadows/campaigns/20260904T034040Z_bed4a462_be2fd1106eed/`
  (gate `noisy-variance-r7`, pre_statement.md Revision 7) closed
  **FROZEN-CERTIFIED** / `NO-FLIP-CONFIRMED-IN-SCOPE` at
  2026-09-04T03:41:22Z. Closed-form `fractions.Fraction` arithmetic over
  the frozen grid — no sampling, no RNG, single core; all six frozen
  artifacts re-hash clean against the dir's `sha256s.txt`.
- Grid: d ∈ {4, 8, 16} × two channels (depolarizing 60 cells, amplitude
  damping 72 cells) × three anchors (GHZ projector, maximally-mixed
  Pauli-Z, GHZ Pauli-X) = 132 unique cells. All 31 preregistered
  regression identities hold (Eq. (37) == direct second-moment ratio and
  Eq. (44) variance ratio ≥ second-moment ratio in every row);
  `flip_cells: []`; no ratio ≤ 1 anywhere.
- Extrema with provenance cells: second-moment ratio (Eq. 37) min
  **1955/1539** (depolarizing d=16, p=1, GHZ projector), max **17/9**
  (depolarizing d=16, p=1, maximally-mixed Pauli-Z); variance ratio
  (Eq. 44) min **6815/4671** (depolarizing d=16, p=3/4, GHZ projector),
  max **2** (depolarizing d=4, p=1, GHZ Pauli-X).
- Limitation, frozen before the run: the paper proves the no-flip result
  analytically (arXiv:2608.18935v1, Corollary 3.11 / Eq. (44)). This
  campaign is regression-grade exact confirmation of the formula chain on
  the frozen grid plus margin quantification — not a discovery search,
  not a proof beyond the grid, and no evidence on whether the k=4 bounded
  `cx_chain` edge persists. All earlier bounded-family facts stand
  unchanged.
- Provenance: prereg sha256
  `8f9d51ba6adad68fc415526b7f37c1848ef9e17d0d04ca25b11155d39ad0fbf2`;
  module snapshot sha256
  `156b09c0cdaf75fb46053df91427693520e83af640a430db8b23a3f659ada5a0`;
  results sha256
  `fcb733b705b8535acb829f105aa5cf89b83fa077b723687befd7be1a2ddc82d2`;
  summary sha256
  `36885050aa43400e92cde27991fd1cbb718a03ec06c8e43a3bfd3591f4675227`.
  Superseded predecessor `20260904T031748Z_cbbe7c0b_fbc8a23e2ed9` (r4
  instrument, closed SUPERSEDED before any result was used; Revisions
  5–7 hardened the instrument only). Post-close scoped re-check:
  `shadows_noisy_variance.py --check` reproduces identities=31 failed=0,
  FROZEN-CERTIFIED, 132 rows, min ratios 1955/1539 and 6815/4671.
- Ledger sync: RESULTS.md shadows row and evidence ladder now carry the
  R7 entry; README follow-up section records the closure with the frozen
  scope; state.json now carries `latest_campaign`
  `20260904T034040Z_bed4a462_be2fd1106eed` / `FROZEN-CERTIFIED` (exactly
  what `campaign.py state refresh` derives from the run dir). No claim
  beyond the grid is made anywhere.

## 2026-09-04 — while GB9 decodes: qlops common-unit correction; shadows claims bounded

Independent read-only review found and corrected provenance errors in two
CPU-light targets; the 26-thread GB9 decoder was not interrupted.

**`qlops/` Revision 4 — common-unit accepted-output accounting.** Revision 3's
4.1×–632× band divided incompatible time units: zero-level depth 24 is circuit
layers, while the Litinski denominator is in syndrome-extraction cycles. The
frozen 193.3×–29,946.5× row is now explicitly
`[DERIVED; LEGACY-MIXED-UNITS]`, not a published magnitude. The corrected
deterministic comparison uses:

- Litinski accepted CCZ equivalent =
  `7 × Table-6 unit qubits × reported cycles`; Table 6 states that its cycle
  count already includes postselection.
- Zero-level accepted CCZ =
  `(22 + 3(2d²−1)) × 3 syndrome rounds / acceptance`, for the paper-simulated
  d=3 and d=7 variants.

The worst sourced arm is d=7 with acceptance 0.30:
**33.76693322683706×–5086.233185367411×**. Every arm still clears the
preregistered 2× falsifier, so the Gate-B comparability verdict survives.
The acceptance endpoints are approximate reported values cross-applied only as
sensitivity arms. Figures 10/11 do contain graphical pointwise error bars; the
paper gives no definition/confidence level, fit-coefficient uncertainty,
residuals, or decoder. Artifact
`qlops/src/evidence/zero_level_sensitivity.json`, sha256
`1558bf3ebc493f1c3eff02bba13117c647537b78235a4c2d832ac340881e01a5`.
Seven behavioral tests pass, including exact combination-set and worst-arm
guards; regenerated Gate-A and Gate-B JSON is byte-identical to the frozen
campaign.

**`shadows/` — finite evidence bounded exactly.** The k=4 improvement belongs
to the pure `cx_chain` and `cx_mesh_lex` maps, not their composition with
`U_ct`; the literal composite is an excluded map with weight 1/81. The frozen
extension proves EDGE-VANISHES only at k=9..12. Across those rows, the absolute
weight gap narrows while w(`cx_chain`)/w(`U_ct`) decreases, so the relative
and shadow-norm disadvantage widens. `largest_k_completed=12`; incremental
row times 0.468/1.535/5.042/16.356 s grow by about 3.25×, not 2×. **Decision:
do not run a further finite-k campaign.** More rows cannot change the
already-closed k≤8 Gate-C verdict and cannot settle the unresolved asymptotic
law; no claim is made beyond k=12, and family (b) was not searched past k=8.
The authoritative results table now includes the missing shadows row. The live
next question is the channel-resolved real-Clifford variance table
(arXiv:2608.18935).

## 2026-09-02 — GB9 launched: the paper's actual 17x (beam64_32res, num_results=32), paired 1e8, RUNNING

- Decoder built and proven before any sample. The paper's Algorithm 3 step 3
  inserts the seed-BP solution and returns it only if num_results=1; the
  harness Python reference returned on seed convergence regardless — the
  one deviation, fixed (`qldpc-dec/src/qldpc_dec/beam_search.py`). Weight is
  now a sequential index-order sum (np.dot delegates to BLAS, unmirrorable).
  `beam8.cpp` gained `--num-results`; built to a NEW binary `beam_nr_cpp`
  (sha `297db977a57b…`); the frozen `beam8_cpp` (`a480050b2041…`) is untouched.
- Equivalence gate: `beam64_32res` **2000/2000 identical at p=1e-3, 300/300
  at p=3e-3** vs the Python reference (21 s/shot, sharded over processes).
  K=1 regression on the six frozen gate points: zero mismatches and the C++
  prediction files byte-identical to the frozen binary's (six of six).
- The authors' code (github ionq-publications/BeamSearchDecoder,
  `beam_search.hpp`) agrees with the pseudocode on seed insertion, duplicate
  counting, sequential weight, strict `<`, best-so-far after max_rounds, but
  carries two undisclosed heuristics (converged-branch skip; Tanner-degree
  <=2 exclusion) — pre-registered as first suspects if 17x fails.
- Cost: num_results=32 is 188 ms/shot/thread (42 BP runs, 465 iterations per
  shot; every shot enters the search). Sizing from frozen GB7 rates
  (`qldpc-dec/src/evidence/gb9_sizing.json`, sha `5777193c…`): N = 1e8
  pinned — P(refute 17x | K=1-level) = 0.88, P(exclude K=1 baseline | true
  17x) = 0.915; the +-30% band needs ~5.6e8 shots and is disclosed as
  unreachable. ~220 h of decoding.
- Revision GB9 recorded 2026-09-02T10:56Z; `campaign.py init` minted
  `qldpc-dec/campaigns/20260902T114747Z_69eba724_61e9fe30ba5d` (prereg sha
  `45d7a47f…`) at 11:47Z; runner `src/gateb_gb9_32res.py` launched under
  the process supervisor (shard-wise decoding, resumable). Hard voids:
  stream byte-identity with GB7; beam8 predictions byte-identity with GB7;
  beam8 mask identity. No terminal verdict yet.

## 2026-09-02 — RETRACTION: "beam64 17x refuted" was a prereg mis-attribution; beam64_640iters is CONSISTENT with its published 7.0x (Revision GB8, no sampling)

- Trigger: GB7 refuted the 17x point but its interval's upper edge (12.65x)
  overlapped the GB5a band [11, 22] — against 11x/17x/22x the measured 8.26x
  [5.97, 12.65] contains 11x and excludes 17x and 22x, so the band was not
  excluded. Before sizing a band-exclusion campaign,
  the paper (arXiv:2512.07057) was re-read first-hand, v1 (2025-12-08) and
  v2 (2025-12-17), plus an independent second extraction; both agree.
- Finding: Table 1 has FOUR rows. Section III, para. 2: "beam32_340iters and
  beam64_640iters reduce the logical error rate by factors of 5.6x and 7.0x,
  respectively. Notably, the most advanced configuration,
  beam64_32res_640iters, achieves a 17x improvement over BP-OSD." Our beam64
  rung (64, 40, 30, 20, num_results=1, pinned in every frozen manifest and
  in the C++ binary) is `beam64_640iters` → published **7.0x**. The 17x is
  `beam64_32res_640iters` (num_results=32: collect 32 valid solutions,
  return the minimum-weight one), never run. The original Gate B section of
  the pre-statement (frozen 2026-08-29) transcribed the abstract's "beam
  width of 64 achieves 17x" onto the wrong row; the GB5a band [11, 22] and
  GB7's 1/17 target inherited it. Not a version drift. (A first draft of
  this entry attributed the transcription to Revision GB1, which pins beam8
  only; corrected the same day by an independent read-only audit.)
- Re-target (`qldpc-dec/src/gateb_gb8_retarget.py`, reads only the frozen
  GB7 summary, hash-checked; replays GB7's rule under the targets as run and
  must reproduce both frozen outcomes or abort): with 7.0x, reciprocal
  0.1429 ∈ frozen CI95 [0.0791, 0.1674] → beam64 outcome
  `NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE`, identical to
  beam32's. Measured 8.26x [5.97, 12.65] contains 7.0x. Artifact
  `qldpc-dec/src/evidence/gb8_retarget.json`, sha256
  `d432ea399707112075e2e66b763b66fe41ad9fe05b1343e9c6003a3a949c8dd1`.
- Verdict: beam64_640iters **CONSISTENT with the published point value**.
  The paper publishes no shots, failure counts or error bars, so nothing
  tighter than the point value can be claimed. 17x: NOT TESTED. The
  2026-09-01 entry below stands as the frozen record; its heading carries an
  in-place retraction annotation (original text preserved), its body is
  unchanged.
- Recorded source conflict (unresolved): abstract "same LER" vs Section III
  "1.3x" for beam8 vs bp30+osd. The 2026-08-29 pre-statement and GB7 use the abstract reading (GB5a
  measured 1.024 [0.660, 1.593]). Under 1.3x the implied factors vs beam8 are
  4.31x (inside our beam32 interval) and 5.38x (below our beam64 lower bound
  5.97x — measured benefit exceeds implied). No verdict drawn from it.
- Unchanged: every frozen artifact, GB7's rule and FROZEN-NEGATIVE label,
  all counts, cells and intervals; GB5a's beam64 INCONCLUSIVE (its CI
  contains 1/7.0 as it contains everything; the [11, 22] band is void).
  Source constants `RUNGS['beam64']` / `PUBLISHED_RECIPROCALS['beam64']`
  keep their values (the verifier's frozen-decision replay pins them) and
  now carry a GB8 pointer comment.
- Acceptance re-run with six new `gb8.*` re-derivation checks and the
  Revision-GB8 presence check: `qldpc-dec/scratch/verify_gb8.log`,
  QLDPC_GATEB_ACCEPTANCE_PASS, **87 PASS / 0 FAIL / 0 SKIP**.
- Step 2 of the planned band-exclusion campaign (3e8-4e8 shots) is NOT
  launched: the question it would answer no longer exists. What would test
  the 17x is a new `beam64_32res` rung (binary `--num-results` path + fresh
  equivalence gate + its own pre-registration).

## 2026-09-01 — Gate B decided: beam32 5.6x reproduced, beam64 17x refuted (GB7 paired 1e8) — beam64 headline RETRACTED 2026-09-02 (Revision GB8, see above)

- GB7 closed **FROZEN-NEGATIVE** (verdict keyed on the beam32 gap hypothesis,
  which is refuted): `qldpc-dec/campaigns/20260901T145247Z_b7ea9ac4_ac3f6689e03b/`.
  100,000,000 shots of the frozen bp30+osd stream; the first 2e7 shots, all
  three beam prediction files, and all three failure masks byte-matched the
  frozen GB5a paired artifacts (prefix failures 40/8/5 recovered exactly).
  Failures per 1e8: beam8 223, beam32 41, beam64 27.
- beam32 (paired vs beam8, four-cell bootstrap 100,000 draws): R32 = 0.1839,
  CI [0.1325, 0.2400]; the published 1/5.6 = 0.1786 is INSIDE. Measured
  factor **5.44x [4.17, 7.55]** vs published 5.6x — the published beam32
  expansion is quantitatively reproduced. Width effect: b=193 fixed, c=11
  reversed, exact McNemar 4.0e-44. Outcome
  `NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE`.
- beam64: R64 = 0.1211, CI [0.0791, 0.1674]; the published 1/17 = 0.0588 is
  BELOW the interval. Measured factor **8.26x [5.97, 12.65]** vs published
  17x — **GAP_CONFIRMED**. b=202, c=6, McNemar 5.2e-52. beam64's failures are
  a strict subset of beam32's (14 fixed, 0 reversed, McNemar 1.2e-4): the
  32->64 step is real but ~1.5x, not the ~3x the published ladder implies.
- Every determinism assumption was backstopped: a deliberately unaligned
  rehearsal was voided by the byte check pre-launch, the runner now refuses
  non-chunk-aligned shot counts, and the frozen GB5a `gateb_ladder.py`
  snapshot was shown byte-equivalent to the streaming rewrite at multi-chunk
  scale before any GB7 sample was drawn.
- Independent acceptance over all three frozen campaigns (raw masks
  recomputed from predictions, 529 + 128 + 30 frozen hashes re-checked,
  decisions re-bootstrapped from re-derived cells):
  `qldpc-dec/scratch/verify_final_gateb.log` —
  QLDPC_GATEB_ACCEPTANCE_PASS, 80 PASS / 0 FAIL / 0 SKIP.
- Verifier fix made during closure: its "control plane frozen" check was an
  existence test that could only pass post-close; replaced by a
  lifecycle-aware check that visibly SKIPs pre-close and, post-close,
  re-hashes every file in `sha256s.txt` (strictly stronger). No prereg rule
  or frozen artifact was touched.

## 2026-09-01 — Gate B ladder frozen INCONCLUSIVE; mechanism companion certified; GB7 escalation running

- The GB5a 2e7-shot ladder closed **FROZEN-INCONCLUSIVE**: BP+OSD 41,
  beam8 42 (ratio 1.024, CI [0.660,1.593]), beam32 10 (0.244, [0.102,0.452]),
  beam64 5 (0.122, [0.025,0.263]). No published band was recovered or
  excluded — independent-arm Poisson noise binds both sides of every ratio,
  which cannot separate the measured 4.1x/8.2x reductions (vs BP+OSD) from the published
  5.6x/17x. Frozen campaign: `qldpc-dec/campaigns/20260901T040018Z_55e4ac01_f81265905380/`.
- The embedded paired instrument (same-stream decode) found a real width
  effect: beam32 resolves 35 of beam8's 40 shared failures (3 reverse,
  exact McNemar 6.7e-8); beam64 36 (1 reverse, 5.5e-10). On identical shots
  R32 = 8/40 = 0.200 vs published 0.1786; R64 = 5/40 = 0.125 vs 0.0588.
- Control-plane activation was post-hoc for GB5a (it ran inside a sibling
  target's claim); frozen manifest merged, 529 artifacts hashed,
  `FROZEN-INCONCLUSIVE` recorded; the accidental empty mint dir was closed
  REHEARSAL. Disclosed in `qldpc-dec/pre_statement.md` Amendment.
- GB6 dense paired companion at p=3e-3 (2e5 shared shots, argrescale
  circuit) closed **FROZEN-CERTIFIED**: failures 654 (BP+OSD) / 464 (beam8)
  / 145 (beam32) / 80 (beam64); M1-M3 PASS; verdict EXPANSION_EFFECT_PRESENT —
  the width effect is structural, not sparse-count noise:
  `qldpc-dec/campaigns/20260901T135902Z_71377291_ddeb5cde09fe/`.
- GB7 paired escalation (1e8 shots of the GB5a frozen stream prefix +
  8e7 more, revisions GB1-GB7, `src/gateb_gb7_paired.py`) is RUNNING:
  `qldpc-dec/campaigns/20260901T145247Z_b7ea9ac4_ac3f6689e03b/`. Decision
  rule: exact McNemar width effect AND paired-ratio interval excluding 1.0
  AND excluding the published reciprocal (1/5.6, 1/17) => GAP_CONFIRMED.
  Shot/prediction/mask prefix must byte-match the frozen GB5a artifacts;
  mismatch voids the run. Pre-launch verification: cross-writer byte
  equivalence (frozen campaign snapshot vs streaming rewrite) at
  multi-chunk scale, 6/6 equivalence reports re-hashed, independent read-only
  review with NO_HIGH_CONFIDENCE_ISSUES.

## 2026-08-31 — GB3 overturns the historical OSD0 timing attribution

- The 2026-08-30 timing forensics below used the now-invalid scalar-mean BP
  prior. Under that prior, OSD-CS10 exhibited an expensive fallback mode on
  many shots while OSD0 bypassed that cost, leading to the incorrect
  inference that the paper's Table-II mean represented OSD0 semantics.
- The authoritative corrected-prior OSD-CS10 campaign instead reproduces both
  means and both tails: mean/p99.9 ratios versus the paper are 1.082/1.115 at
  p=5e-4 and 1.038/1.065 at p=1e-3.
- A focused 3,000-decode diagnostic on the first shots of those same seeded
  streams kept the exact heterogeneous channel but changed only to OSD0. It
  produced mean/p99.9 = **1.208/5.986 ms** at p=5e-4 and
  **1.875/6.147 ms** at p=1e-3 — respectively 0.34/0.022 and 0.18/0.021
  times the published values. OSD0 therefore does not reproduce the
  published tail protocol under the corrected prior.
- Correction: retract the claim that the published timing row is BP+OSD0.
  The evidence supports the published OSD-CS10/max-iter-30 semantics once the
  full prior vector is supplied. Gate B independently retains exactly that
  OSD-CS10/order-10/max-iter-30 denominator for Fig. 2 comparability.
- The live harness now uses one `IONQ_BPOSD_DECODER_CONFIG` value for its
  timing and ordinary BP+OSD LER arms and records the resolved configuration
  in manifests/results. An end-to-end smoke verified the exact five fields;
  the frozen Gate-B manifest independently records the same values.
- Diagnostic artifact:
  `qldpc-dec/scratch/gatea_vector_prior_osd0_probe.json`; authoritative timing
  amendment:
  `qldpc-dec/campaigns/20260831T111233Z_755b3791_641ea5506f3b/`.


## 2026-08-31 — Corrected Gate A closed: timing reproduced, accuracy descriptive

- Repeated the paper-pinned BP+OSD timing protocol after isolating all other
  decoder campaigns: one compiled decoder and 10,000 individual
  `perf_counter` calls per point, with the full 8784-value heterogeneous
  channel asserted exactly.
- At p=5e-4, mean/p99.9 are **3.842/303.767 ms** versus paper
  3.55/272.5 ms (ratios **1.082/1.115**). At p=1e-3 they are
  **10.995/307.796 ms** versus 10.59/289.0 ms (ratios **1.038/1.065**).
  Both satisfy the frozen symmetric factor-2 criterion: **REPRODUCED**.
- The source summary omitted the p=5e-4 verdict because its reference key was
  `"5e-04"` while the measurement key was `"0.0005"`. Revision GA3 disclosed
  the bug after measurement, required exact key-set equality, and froze a
  hash-checked amendment without rerunning or selecting samples.
- Corrected accuracy is descriptive because no matching published target
  exists: all six X/Z rows at p=3e-4, 5e-4, and 1e-3 observed **0/100,000**
  failures; each per-round 95% Wilson upper bound is 3.2011e-6.
- The primary p=3e-4 X row was invalid because the old loader silently used
  the derived Z circuit. GA2 replaced only that row, under the same seed, with
  a separately derived and provenance-hashed X circuit/DEM.
- Frozen evidence: timing source
  `qldpc-dec/campaigns/20260831T110541Z_e3e01d3d_774ea68fcaa7/`, authoritative
  GA3 verdict
  `qldpc-dec/campaigns/20260831T111233Z_755b3791_641ea5506f3b/`, and
  authoritative GA2 accuracy
  `qldpc-dec/campaigns/20260831T103828Z_160bc8cc_05a5ebaec745/`.

## 2026-08-31 — Corrected Gate B beam8 band closed inconclusive

- The corrected p=1e-3 Z campaign completed 1,000,000 deterministic shots per
  arm. Its denominator kept the published OSD-CS10/order-10/max-iter-30
  configuration. BP+OSD had **4 failures** (3.3333e-7 per round); bit-exact
  C++ beam8 had **1 failure** (8.3333e-8 per round).
- The beam8/BP+OSD point ratio is **0.25**. The frozen 100,000-draw
  independent-binomial bootstrap interval is **[0, 2]**, which overlaps the
  published-equivalence band [0.87, 1.15]. The pre-stated beam8-band decision
  is **INCONCLUSIVE**, not confirmation or refutation.
- Beam32 and beam64 remain unrun because the frozen ladder gated them on
  beam8 passing. The 20 concurrent BP shard times are aggregate throughput,
  not Gate-A single-shot latency evidence.
- All 20 BP shard ranges, failure masks, hashes, counts, and the beam
  prediction/mask hashes were independently recomputed. The source samples
  are byte-identical to the invalidated scalar-prior campaign; only the BP
  decoder channel changed.
- Frozen evidence:
  `qldpc-dec/campaigns/20260831T091453Z_fd47e8d5_f3b51514a09e/`. The original
  scalar-prior run remains at
  `qldpc-dec/campaigns/20260831T075902Z_a34b4014_b5d11010cd55/` with outcome
  `INVALIDATED_GB3_scalar_prior`; its 66 BP failures support no ratio claim.

## 2026-08-31 — Campaign closure made fail-closed and immutable

- A NumPy scalar boolean made the original Gate-B JSON close fail after all
  results existed. GB4 disclosed the recovery; verdict booleans are now
  converted to Python `bool` before serialization.
- Campaign writes now reject non-finite JSON before opening a destination,
  use deterministic gzip metadata, refuse an existing campaign directory,
  require the manifest before results/close, and reject every mutation after
  close. Config is deep-copied at creation, re-hashed before manifest freeze,
  and canonical manifest fields cannot be overridden through extras.
  Regression checks covered strict JSON, external/direct config mutation,
  reserved-field overrides, create collisions, deterministic compressed
  bytes, and lifecycle transitions.

## 2026-08-31 — Gate C BP comparators restored under GB3

- Added `gatec_exact.run_campaign`, which freezes pre-run config/source
  snapshots, circuit/DEM hashes, exact BP channel assertions, per-shot CSVs,
  compressed results, and immutable-close inventory.
- Reran both Z and X memory at p=0.01/0.03, 2000 shots per point. Every
  exact-ML, beam8, BP+OSD, and NMS count matches the previous campaigns.
  All four corrected per-shot CSVs are byte-identical to their historical
  counterparts.
- Reason: these code-capacity DEMs have 36 mechanisms with one uniform prior
  q=2p/3, so the former scalar mean equals the full vector exactly. This
  restores the Gate-C BP rows; it does not rescue heterogeneous circuit-level
  Gate-A/B runs.
- Corrected rates remain:
  - Z p=0.03/0.01: exact ML 0.0145/0.0005, beam 0.0325/0.0040,
    BP+OSD 0.0330/0.0055, NMS 0.0260/0.0020;
  - X p=0.03/0.01: exact ML 0.0145/0.0005, beam 0.0370/0.0045,
    BP+OSD 0.0365/0.0025, NMS 0.0230/0.0010.
- Frozen evidence: Z
  `qldpc-dec/campaigns/20260831T093503Z_2f2ed071_884b7d216c6f/`; X
  `qldpc-dec/campaigns/20260831T093513Z_3039c8d0_4a350b3b0b06/`.

## 2026-08-31 — Corrected Gate D gross calibration closed

- Reran the [[144,12,12]] derived p=0.003 gross calibration with the corrected
  forward pre-transition AIS path product, 13 candidate classes, T=64, K=64,
  q0=0.02, B=2500, and the frozen base seed.
- Result: **200/200** beam seeds syndrome-consistent, **200/200** decisions
  certified, **199/200** certified decisions observable-correct. All 200 chose
  candidate class 0. Margins: minimum 46.52 nats, mean 71.63; median wall
  time 13.43 s/shot.
- Verification read all 200 JSONL records: shot indices are exactly 0…199,
  aggregates reproduce the summary, the compressed result has one matching
  record, and every copied source/input hash matches the pre-run manifest.
- Scope remains deliberately narrow: this is a candidate-set calibration,
  not a global-ML proof. Because no shot selected a nonzero class, it still
  does not exercise contested-class arbitration. The prior endpoint-weight
  artifacts remain retracted rather than being overwritten.
- Frozen evidence:
  `qldpc-dec/campaigns/20260831T083647Z_51cd3b47_26199282c466/`
  (manifest, copied source and inputs, 200 per-shot records, summary,
  compressed result, immutable-close inventory).

## 2026-08-31 — BP+OSD prior bug corrected; GB3 replacement launched

- First-hand source audit of the paper-pinned `stimbposd==0.1.0` tag (commit
  `7921f5eb1b358ff616f9822280c9961e83df06cb`) found that
  `stimbposd.BPOSD` passes the merged DEM prior vector as
  `error_channel=list(priors)`. The local adapter instead averaged that vector
  and passed one scalar `error_rate`.
- This is not an equivalent parameterization. The pinned p=1e-3 Z DEM has
  8784 merged columns, 9 distinct probabilities, and range
  [0.000533333333333148, 0.003721584619834077]. A deterministic two-column
  regression now observes configured priors [0.1, 0.2], not the former
  [0.15, 0.15], and the full configured channel equals the merged vector
  element-for-element.
- Impact at this ledger point: every historical local BP+OSD accuracy/runtime
  value used the wrong scalar channel. Gate A's BP+OSD claims and the BP
  comparator rows in the Gate-C exact campaigns were retracted pending the
  corrected reruns now closed above. Exact-ML, beam, and NMS outputs were
  unaffected. The original Gate-B million-shot campaign could not decide the
  gate.
- `qldpc-dec/pre_statement.md` Revision GB3 discloses that this correction was
  found after interim Gate-B output existed. No gate threshold or stochastic
  input changed. At this ledger point, corrected campaign
  `qldpc-dec/campaigns/20260831T091453Z_fd47e8d5_f3b51514a09e/` was running
  with source snapshots and an exact prior-channel assertion; it is closed
  as the inconclusive beam8-band result above.
- The corrected campaign's independently regenerated sample files are
  byte-identical to the original frozen streams: BP sha256
  `fd0e018d…f06788`, beam sha256 `ba9d9696…e4c5fe`. A 200-shot smoke campaign
  closed immutably with all 8784 configured priors verified.
- Separate backend validation on the actual million-shot beam stream checked
  199 uniformly spaced shots plus every C++ logical-failure shot: 0/200
  Python/C++ prediction mismatches; both backends fail the same selected shot
  (campaign index 280665).

## 2026-08-31 — Gate C exact reference completed in X memory

- Extended the [[36,4,4]] exact code-capacity runner from Z-only to an
  explicit `--basis X|Z` contract. Z memory uses X errors with `(Hz,Lz)`;
  X memory uses the Hadamard-dual Z errors with `(Hx,Lx)`.
- Regression proof: the refactored Z circuit string is byte-for-byte equal to
  the frozen 2026-08-30 implementation. X DEM dimensions are 18 detectors ×
  36 mechanisms with 4 observables and every prior equal to q=2p/3.
- Frozen X campaign, 2000 shots/point:
  - p=0.03: exact ML **0.0145** vs beam8 0.0370 / bp30+osd 0.0365 /
    nms-ens24 0.0230;
  - p=0.01: exact ML **0.0005** vs beam8 0.0045 / bp30+osd 0.0025 /
    nms-ens24 0.0010.
- Exact ML therefore beats every heuristic arm at both points in both CSS
  sectors. This remains same-family small-code code-capacity evidence, not a
  [[144,12,12]] circuit-level claim.
- Evidence:
  `qldpc-dec/campaigns/20260831T085112Z_df905ed6_9b30cfe313ed/`
  (pre-run manifest, source hashes/snapshots, two aggregate JSONs, two
  per-shot CSVs, compressed results, immutable-close inventory).

## 2026-08-31 — Gate B sparse-ratio bootstrap tail bug fixed pre-outcome (GB2)

- `qldpc_dec.bootstrap.ratio_ci()` dropped all non-finite bootstrap ratios.
  In sparse data this conditions away positive-numerator/zero-denominator
  draws and can replace an unbounded upper endpoint with a falsely finite one.
- Deterministic reproduction with the frozen 100000-draw design:
  beam/BP counts 5/1 in 1e6 returned **[0.667, 9.0]**; corrected handling
  returns **[1.0, +inf]**. Failure count, not total shots, controls this edge
  probability.
- Fix retains positive/zero as `+inf`, discards only undefined zero/zero
  draws, and uses an infinity-safe linear percentile. Dense 50/50 check stays
  finite at [0.672, 1.487]; finite-only percentiles match NumPy exactly.
- Revision GB2 was appended to `qldpc-dec/pre_statement.md` before any
  million-shot shard result existed. Samples, decodes, target band, 100000
  draws, bootstrap seed, and three-way verdict rule are unchanged.
- The active process already imported the old helper. Its raw counts remain
  evidence; after immutable close, a separate campaign amendment will apply
  GB2 and supersede only the original ratio CI/verdict.

## 2026-08-31 — Gate D gross evidence retracted; corrected AIS closes the small-code exact check

- Root cause: `AISEngine.run()` used
  `logW = log gamma_1(u_T) - log gamma_0(u_0)` after intermediate
  Metropolis moves. AIS weights are a path product: each
  `gamma_beta_t/gamma_beta_{t-1}` ratio must be evaluated at the same
  pre-transition state. Markov moves do not telescope pathwise.
- Deterministic reproduction on a one-bit model: the old estimator returned
  **1.8310097** for exact Z=**1.1111111** (relative bias **+64.79%**, K=200k).
  An independent read-only Fable audit reached the same diagnosis and
  recurrence.
- Fix: accumulate
  `-delta_beta * (E(u_{t-1}) + log p0(u_{t-1}))`, then apply the
  beta_t-invariant sweep. Energy and base-state Hamming weight are updated
  incrementally on accepted flips, avoiding a full mechanism scan per stage.
- Consequence: the prior [[144,12,12]] Gate-D pilot/calibration entries
  (49/49 and 200/200 certified) used biased weights and are **RETRACTED**.
  They support no current certificate claim; historical artifacts remain
  untouched.
- Corrected checks:
  - one-bit exact partition: relative error -0.061% at T=1 and -0.0028% at
    T=32 (K=400k); cached energies/Hamming weights match full recomputation;
  - d=5 mini exact check: **260/260 AIS decisions match exact ML** and
    260/260 satisfy the paired-bootstrap criterion;
  - strict BB [[36,4,4]] code-capacity check, p=0.03, all 16 logical classes
    represented on every shot: **191/200** AIS point decisions match exact
    Fraction ML; **185/200** satisfy the criterion; **185/185 certified
    decisions match exact ML**, 0 certified-wrong, and all 9 point-decision
    errors are uncertified.
- Frozen corrected evidence:
  `qldpc-dec/campaigns/20260831T083154Z_71fb5839_0b7c2262d96d/`
  (pre-run manifest, source hashes/snapshots, full per-shot exact masses,
  summary, immutable-close inventory). Scope is explicitly small-code
  same-family evidence, not a [[144,12,12]] headline verdict.
- Next Gate-D action: rerun the gross 200-shot calibration with corrected AIS;
  until then Gate D is closed only at the exact small-code validation scope.

## 2026-08-31 — qldpc-dec absolute LER normalization corrected (Gate B ratios unchanged)

- First-hand check of IonQ's official
  `ionq-publications/beamsearchdecoder/simulation_functions.py` confirms:
  `x_logical_error_rate = (x_num_errors/x_num_shots)/distance` and likewise
  for Z, with `distance_dictionary[(144,12)] = 12`.
- Harness bug: `qldpc_dec.run_gate.run_point()` stored raw
  `any-logical failure/shots` under `ler` without `/12`; absolute circuit-level
  LER labels from that path were therefore **12× high**. Frozen artifacts are
  not edited.
- Source now records both `raw_logical_failure_rate` and source-defined
  per-round `ler = raw/12`, with both Wilson intervals. Controlled behavioral
  check: 12 failures/1200 shots => raw 0.01, LER 0.0008333333333333334,
  and every CI endpoint scales exactly by 1/12.
- Frozen Gate-A accuracy companion at p=1e-3 is therefore: raw shot-failure
  rates **7e-5 X / 4e-5 Z** (7/4 failures in 1e5), per-round LER
  **5.833e-6 X / 3.333e-6 Z**, total **9.167e-6**. Zero-failure points stay
  zero. Gate-B beam/BP ratios are invariant because `/12` cancels; the
  p=3e-3 companion's 0.1995/0.0015 numbers are raw shot-failure rates.
- Scope: this corrects the 12-round IonQ circuit harness only. The [[36,4,4]]
  code-capacity exact reference is one-shot and already labels `fails/shots`
  appropriately.

## 2026-08-31 — qldpc-dec pinned Gate B million-shot campaign launched (pre-data)

- Exact target frozen before decoding in `qldpc-dec/pre_statement.md`
  Revision GB1: p=1e-3, Z basis, 1e6 shots **per arm**, beam8/BP+OSD target
  band [0.87,1.15].
- Original independent arm streams retained:
  `sampling_seed(20260829, p, "Z", arm)`. BP+OSD shards a single already
  sampled stream into 20 contiguous 50k slices; beam8 remains single-thread.
- Beam backend is `src/beam_cpp/beam8_cpp`, already 0/1000 bit-exact against
  the frozen Python arm at both p=1e-3 and p=3e-3. BP arm remains the pinned
  max_iter=30/min-sum/parallel/OSD-CS10 order-10 decoder.
- Decision frozen pre-data: 100000-draw independent-binomial bootstrap;
  CI inside band => REPRODUCED, outside => NOT-REPRODUCED, crossing =>
  INCONCLUSIVE. Raw shot-failure probability and per-round LER (=raw/12)
  are both recorded; their ratio is identical.
- New runner `qldpc-dec/scratch/gateb_pinned.py` smoke-closed end-to-end at
  200 shots (2 shards): bit-packed seeded streams matched ordinary Stim
  sampling exactly; BP+OSD 56.50 ms/shot, C++ beam8 3.18 ms/shot; expected
  zero-failure/UNDERDETERMINED smoke verdict. Full-run estimate ≈1 wall-hour
  with 16 BP workers + one single-thread beam process on this 28-core host.

## 2026-08-30 — msd floor EXPLAINED from first principles; output-side armoring provably cannot fix it

- Advisory pointed at `gate_a.py:237-239`: each observable is
  `output_x_records[k]` XOR the RAW `source_x_records[q]` over LX[k] —
  a bare parity of the distance-2 [[8,3,2]] source block. So every
  observable has ≥2 independent single-fault paths and output-side
  armoring closes only one. Verified verbatim in source, then tested by
  enumeration instead of argument.
- Enumerated every weight-1 undetectable logical mechanism in the
  unpatched Gate-A DEM (no decomposition, gauge detectors allowed):
  exactly **7 bare mechanisms**, whose observable patterns are precisely
  the 7 nonzero subsets of {0,1,2} — matching source wires 0..6 under
  LX = {0,1,2,3}/{0,1,4,5}/{0,2,4,6} (wire 7 is in no support and
  contributes 2 locations only). [NUMERICAL]
- **The floor is quantitatively the bare mass**, i.e. the measured
  "logical error rate" is the probability of these weight-1 paths:

  | p | bare mass (enumerated) | frozen measured p_L | gap |
  |---|---|---|---|
  | 1e-4 | 0.038781 | 0.038356 | 1.1% |
  | 3e-4 | 0.114423 | 0.110333 | 3.7% |
  | 1e-3 | 0.360266 | 0.319565 | 12.7% (union overcount at large mass) |

  This supersedes the earlier "slope 0.992 ⇒ distance-1 frame" reading
  with an exact mechanism-level account: not just the exponent, the
  coefficient too.
- Fault-location decomposition of those 7 mechanisms (via
  `explain_detector_error_model_errors`): **522 locations on source
  DATA wires, 775 on the three output reps, 167 elsewhere**. Therefore
  output-only armoring removes at most 53% of the locations — it can
  change the COEFFICIENT by ≈2× at best and **cannot change the
  exponent**. The measured patch outcome (slope 0.9902, and 2.7× worse
  because it added anchorless-frame channels) is exactly this prediction.
- Route status: the output-patch family is now closed **by construction**,
  not merely by one failed build. The only remaining unblock is
  fault-tolerant readout of the SOURCE frame (measuring the source X
  parities through a distance-carrying structure) — a full rebuild whose
  decoder and scheduling would be [INFERENCE] against a paper that names
  no decoder. Standing conclusion unchanged: discharge the qlops
  `ZERO_LEVEL_CCZ` REPORTED tag from the paper's own numbers, not from a
  surrogate.
- Scope caveat retained: even a future slope→2 result would bind only
  this reconstruction, since arXiv:2605.21867 specifies no decoder.

## 2026-08-30 — na-compiler cost rows upgraded to GLOBAL µs optima (and a heap bug I wrote, caught by arithmetic)

- Advisory (correct) flagged that the frozen cost rows only proved
  "no cheaper schedule at the minimum batch count". That scope is
  genuinely insufficient: each batch pays 2·15 µs per moved atom, so a
  6-batch schedule of short moves can undercut a 4-batch schedule
  containing one long move — µs-optimality is NOT implied by
  batch-count optimality.
- Settled it directly: `cdijkstra.c` (in the new campaign) runs Dijkstra
  over the FULL reachable 4×4 state space with the frozen edge weight,
  batch generation a line-for-line replica of the frozen `cbfs.c gen()`.
  Results: **n=5 → 539.043122292 µs** (claim 539.043122398901, Δ 1.1e-7)
  and **n=6 → 666.354656219 µs** (claim 666.354655402228, Δ 8.2e-7),
  both inside the 2⁻²⁰ µs fixed-point quantum; witnesses have 4 and 5
  batches (= BFS diameters) and replay through the frozen verifier to
  2.3e-13 / 4.5e-13 µs. **Both rows are global optima; the
  consistent-A* dependency is dropped.** n=6 sweep: 5,756,682 settled of
  5,765,760 reachable, 343 s single-core.
- **Bug I wrote, and how it surfaced.** The first cdijkstra build had a
  sift-down comparing children against the vacated root instead of the
  relocated last element — the heap invariant was destroyed, pops came
  out unordered. Symptoms were loud once I looked: 533M "settles"
  against 5.77M reachable states (92×), a heap pinned at 8.5M entries
  for 9 hours, and an n=3 answer of 1003.54 µs / 8 batches. After the
  fix the same n=3 target solves at 288.780458 µs / 3 batches, matching
  an independent Python Dijkstra (288.78045748) to 1e-6. I killed the
  bad 9-hour run rather than let it produce a number. Lesson recorded:
  a settle count exceeding the state count is an invariant violation,
  not slow convergence — check it as a first-class assertion.
- Frozen: `na-compiler/campaigns/20260830T021851Z_6e74b176_4683ec3d3949`
  (cdijkstra.c + both result JSONs + logs + manifest carrying the bug
  record and the exact-replay verification).

## 2026-08-30 — Gate D calibration complete: 200/200 certified, first certified-ML FAILURE observed and exercised

- 200-shot calibration run (same settings, median 10.9 s/decode, ≈41
  min wall): **200/200 shots certified** (0 beam-seed inconsistencies),
  decision class 0 in all 200 — the certificate pipeline is robust at
  gross generator count and zero stall/seed failures across 200
  consecutive runs. Margins: min 47.9 / median 73.9 / max 88.5 nats.
- **The interesting shot is 151**: certified AND its observable prediction
  was WRONG (199/200 certified-and-correct; certified_obs_fail_rate
  0.005). This is exactly the paper's predicted behavior and the codified
  caveat: certification asserts ML-within-candidates, NOT decoding
  success — ≈14% of certified decisions still fail in the paper's
  d=5 surface experiments at p=0.05 (irreducible ML floor). A certified
  failure at the noise floor of this DERIVED circuit is therefore the
  certificate doing its job, not a bug candidate. Still recorded as the
  first such event in this harness.
- Caveat unchanged and now quantified: all 200 decisions class 0 ⇒ the
  certificate has not yet arbitrated a contested class in 250 total
  decoded shots. Exercising cross-class competition needs either higher
  p or adversarial syndromes; that (plus a global-ML check via the Gate C
  exact reference at small scale) remains the open slice before the
  certificate can be attached to headline gate shots.
- Artifacts: `qldpc-dec/scratch/gateD_cert/gross_result.json` (200-shot
  aggregate overwrites the 50-shot one; the 50-shot run's per-shot JSONL
  rows are preserved in `gross_shots.jsonl` = 256 lines total),
  `gross200.log`, owner reruns `mini_owner_150.json` /
  `gross_owner_6.json`.

## 2026-08-30 — beam8 C++ port adjudicated BIT-EXACT (both p, 1000-shot formal diffs)

- Owner-run formal equivalence (fixed pipeline, 1000 shared shots, seed
  20260830): p=3e-3 (DERIVED rescale circuit, 34 obs):
  **0/1000 mismatches**, C++ 14.5 ms/shot vs Python 2698 twin-run =
  **186×**. p=1e-3 (canonical IonQ circuit sha 106f4547…, 12 obs):
  **0/1000 mismatches**, C++ 3.31 ms/shot vs Python 548 twin-run =
  **166×**. Both reports committed under `src/beam_cpp/scratch/`.
- The Port's earlier 9/10-mismatch report (`equiv_p0.003_n10`) predates
  its three bug fixes and is superseded (kept in scratch as the bug
  ledger). The port is now adjudicated bit-exact at both operating
  points, single-threaded.
- Consequence (updates the 2026-08-30 cost correction): pinned-p beam
  arm cost for a ≥1e6-shot resolution drops to ≈0.9 single-core hours
  (3.3 ms/shot); the pinned Gate B's real remaining cost is the BP+OSD
  arm (279 ms/shot ⇒ ≈78 single-core hours per 1e6 shots).

## 2026-08-30 — Gate D pilot: certification SURVIVES at gross generator count (100% certified, 75-nat margins)

- Subagent `GateDPilot` implemented the arXiv:2608.25545 spacetime
  sampling route (914 LOC, `qldpc-dec/src/gateD_cert/`): bit-packed GF(2)
  kernel basis + greedy sparsifier (mean generator weight 4.05 on the
  gross DEM — matches the paper's device-i target ~4), lightened
  representatives, CRN-AIS with concentrated Bernoulli(q0=0.02) bridge,
  paired-bootstrap certificate (B=2500, Bonferroni δ=0.05/12).
- Mini sanity (rotated surface d=5, code-capacity p=0.10, 2 classes,
  exact ML by 2^12 coset enumeration per class): AIS decisions
  **260/260 = 100% agreement with exact ML, 100% certified, 0 true
  violations** (acceptance ≥90%) [NUMERICAL]. Owner rerun at 150 fresh
  syndromes: 150/150 agreement, block success reproduced 0.6267 both
  engines.
- Gross pilot (BB_144 Z-memory DERIVED circuit, 936 dets × 8784 mechs,
  **7842 generators** = 4× beyond anything in the paper): **49/49
  decodable shots certified = 100%**, mean margin **75.2 nats** (min
  55.2), median wall 11.4 s/decode single-core nice-10 — ~4× faster than
  the 3.5–4 min linear-scaling estimate. 1/50 shots excluded: beam seed
  not syndrome-consistent (shot 21), an honest dependence of the
  pipeline on seed validity. Owner rerun on 6 fresh shots: 6/6
  certified, margins ≈75.3 nats, same wall.
- The riskiest unknown from the CertDec spec — AIS variance drowning
  class gaps at ~7.8k generators — did NOT materialize: margins are
  huge because lightened inside-class minima sit far below single-logical
  shifts. Honest caveat recorded: all 49 decisions were class 0 (no
  cross-class competition exercised in this sample), so the certificate
  currently asserts ML-within-candidates around the beam center; we
  have not yet seen it arbitrate a contested class. A 200-shot
  calibration (bg job, ~45 min) is now running to exercise cross-class
  decisions and calibrate certified-fraction vs LER.
- Scope discipline: certificates assert ML-within-candidates under the
  assumed DEM; the single-logical candidate set was never validated
  against a global-ML optimum (the paper validated only the 2-logical
  set on BB72 code-capacity). Gate C's exact reference ([[36,4,4]]) can
  later serve as the global check at small scale.
- Artifacts: `qldpc-dec/scratch/gateD_cert/` (mini_result.json,
  gross_result.json + per-shot jsonl, DESIGN.md, owner reruns, precomputed
  sparse basis npz — one-time setup ~2.5 min, reusable). Code sha noted
  in campaign freeze to follow with the 200-shot result.

## 2026-08-30 — msd unblock attempt: expanded d3/d7 output patches are NOT distance-carrying (negative result, frozen as evidence)

- Subagent ran the named unblock from the 2026-08-30 reclassification
  (expand output patches so observables carry code distance). Verdict:
  **the patch does not show quadratic suppression anywhere reachable, and
  is strictly worse than the base surrogate at every sub-saturation p.**
  Per the handoff's no-iteration rule it stopped instead of tuning.
- Primary campaign (d=3, five points 1e-6…1e-4, 2e7 shots/pt, base seed
  20260627, acceptance/sampling conventions exactly as the frozen
  harness): slope **0.9902** (target ≈2), fitted c **957.4** (target
  [150, 600]), rmse_log 0.014 (tight), p_L(1e-4) = 0.103 vs base 0.0384
  — 2.7× WORSE. d=7 variant is even worse (slope 0.47, p_L(1e-4)=0.283):
  extending distance amplifies the patch's own weight-1 frame-flip
  channels [NUMERICAL].
- Root cause [DERIVED then verified]: the patch removes the base's
  weight-1 readout channels (source/output single-MX flips now
  comparison-detected at weight ≥2) but ADDS new weight-1 channels: the
  pre-injection window carries an anchorless random X frame (undetected
  Z copy to the whole codeword), and each temporal-comparison chain has
  an unanchored first round on both output and source legs. A random X
  frame has no deterministic round-1 anchor — the hole is structural in
  the teleportation-based layout, not a tuning artifact.
- Mechanistic surprises recorded by the agent (all verified numerically):
  (1) spatial parity "mx(o)=mx(a)=mx(b)" is WRONG — X-side fan-out
  produces free Z-side terms, per-wire MX values are random in the ideal
  circuit, only 3-wire X parity is deterministic (v2 detectors fired 50%
  on ideal shots; discarded); (2) v1/v3 reproduce Gate-A's teleportation
  relation exactly at p=0 (zero detectors/observables over 4000 shots);
  (3) chunked sampler calls do NOT bit-match a single call — rerun
  single-call; (4) driver accounting bug (counting nonzero rows) fixed,
  superseded artifacts renamed, stem preserved in falsified/ dir.
  Frozen default-circuit hashes re-verified 4× against campaign values.
- Consequence for the target: the teleportation-only unblock route is
  now CLOSED with evidence — no patch of this family (X-frame repetition
  + temporal comparison) can show c·p². An expanded patch variant needs
  deterministic spatial distance (a genuine 2D surface-code patch with
  postselection on a nontrivial code), i.e. a different construction, not
  an iteration on this one. Remaining open route for discharging the
  qlops REPORTED tag: the paper's own full numbers, never a surrogate.
- Acceptance criteria: default (no-flag) behavior unchanged (circuit
  hash equality verified); owner spot-checked artifacts by direct JSON
  re-read (fit slope 0.9902, c 957.41, 5/5 points; d7 slope 0.4741,
  c 24.1) and venv import of expanded.py. Artifacts:
  `msd/scratch/expanded_d3_..._shots20000000_seed20260627.json` (+ d7 3×3e7),
  code `msd/src/expanded.py` (2cb1b4002953) + `msd/src/run_expanded.py`.

## 2026-08-30 — Gate C CLOSED: exact degenerate-ML beats beam8/BP+OSD/nms where checkable ([[36,4,4]])

- Subagent built `qldpc-dec/src/gatec_exact/`: exact coset enumeration over
  the FULL affine syndrome space of the smallest same-family BB instance
  ([[36,4,4]], A=x³+y+y², B=y+x+x³), exact rational class masses (Fractions),
  ML argmax, per-run brute-force equality assertion. 2,000 shots/point,
  p ∈ {0.01, 0.03}, runtime 135 s total single-core.
- Headline [NUMERICAL]: exact ML LER **0.0145** vs nms-ens24 0.0260 /
  beam8 0.0325 / bp30+osd 0.0330 at p=0.03; **0.0005** vs 0.0020/0.0040/0.0055
  at p=0.01. The ML-vs-heuristic headroom is real even where beam8 crushes
  BP+OSD: on this instance ML halves the LER of the best heuristic.
  Regret data: bp30+osd's worst wrong-class decision at p=0.01 carries
  e^18.9 ≈ 1.6e8 likelihood penalty vs ML.
- Owner verification (independent): module brute-force assertion re-run on
  owner seed 20260830 — PASS; dx=dz=4 recomputed; fresh-seed spot check
  (400 shots @ p=0.03, seed 987654321): ML 0.020 < nms24 0.0225 < beam8
  0.030 < bp+osd 0.0425 — ordering reproduced at small sample.
- Scope kept explicit: code-capacity Z-memory sector; circuit-level exact
  inference documented unreachable even at n=36 (2^1.6e4 affine space,
  high treewidth) — the harness's first explicitly bounded-scope exact
  reference; NO claims on the [[144,12,12]] headline gates.
- Frozen: `qldpc-dec/campaigns/20260830T184327Z_53a194ad_e663d5477de1`.

## 2026-08-30 — beam8 C++ port: 184–262× faster, bit-exact on shared shots (agent), formal 1k diff finalizing

- Subagent `BeamCppPort` built `qldpc-dec/src/beam_cpp/beam8.cpp` (~850
  lines: stim DEM parser, bit-packed shot IO, masked min-sum BP + beam
  search, deterministic — the Python decoder has NO RNG, seeds only feed
  the stim sampler). Three real bugs found and fixed during bit-exactness
  work: pair-sum gather order (fancy-index copy ⇒ pairwise), stim shot
  packing (u64-LE 120-byte rows, not flat d>>3), observable repacking.
  Evidence [NUMERICAL]: BP state bit-identical to numpy on BB_144 shot 0
  (nu and sum_llr arrays equal at bit level, same convergence iters=13,
  same min|sum_llr| at same arg); 10-shot shared set p=3e-3 10/10
  per-shot observable predictions identical; C++ **9.0 ms/shot** vs
  Python 1650 ms/shot measured twin-run (184×), vs frozen 2362 ms/shot
  reference (262×). p=1e-3: 3.254 ms/shot decode-only, 1000 shots all
  seed-converged.
- Consequence for the parked pinned Gate B: the ≈650 core-h Python beam
  arm collapses to ≈2.5 core-h at C++ speed — the pinned p=1e-3
  ≥1e6-shot resolution becomes plainly runnable. Awaiting the formal
  1000-shot per-shot diff (agent's process was still running; Main
  adopted it) before adjudicating the port bit-exact.
- **SUPERSEDED same day by the adjudication entry above**: the formal
  diffs landed 0/1000 at BOTH p. The correct figures are beam arm
  **3.31 ms/shot ⇒ ≈0.9 single-core-h per 1e6 shots** (the ≈2.5 here was
  extrapolated from the 9.0 ms/shot 10-shot p=3e-3 measurement), and
  BP+OSD **279 ms/shot ⇒ ≈78 single-core-h**, which is the gating cost.

## 2026-08-30 — historical p=3e-3 companion (BP side retracted by GB3)

- This 2,000-shot companion used the scalar BP prior later invalidated by
  GB3. Its raw output was BP 399/2000 versus beam8 3/2000, but the BP rate,
  0.00752 ratio, confidence interval, and inferred p-dependence are all
  **RETRACTED**. They compare beam8 to a semantically wrong denominator and
  do not pass or fail the p=1e-3 beam8 band.
- The unaffected beam arm remains historical pipeline/performance evidence:
  3/2000 failures and 2.36 s/shot for the Python implementation at p=3e-3.
  That cost motivated the bit-exact C++ port; it is not an accuracy verdict.
  Campaign:
  `qldpc-dec/campaigns/20260830T125840Z_9d7bc086_87ca0ed14eac`.
- Circuit provenance unchanged and explicit: derived-by-rescale from the
  committed IonQ p=1e-3 Z circuit, sha256 ea2de75c77969fcb…, labelled
  DERIVED, campaign labelled `B-companion-highp`.

## 2026-08-30 — campaign-directory inventory corrected (Main self-correction)

- Earlier in this session I described
  `qldpc-dec/campaigns/20260830T005616Z_3402651f_5b4e76ede0ee` as
  "residue from the aborted hub-broker launch (manifest only, no
  process)". **That was wrong.** It is the FROZEN Gate A accuracy
  campaign (manifest + results.json.gz + summary.json + inventory.json,
  all six points). The partial pre-freeze snapshot is the *other* dir,
  `...005856Z_48dbfa91...`. The later Gate-A closure entry states this
  correctly; the launch entry's characterization is superseded by this
  note.
- Full campaign-dir inventory taken (files/bytes per dir). Result:
  `shadows/campaigns/2026-08-29T232903Z_0a8203f2_ed782ad007ad` is not
  stray at all — it is my own owner rerun of shadows gates A+B+C (5
  artifacts, independent reproduction) and had simply never been listed;
  it is now in the RESULTS campaign ladder alongside the producer run.
- Genuinely empty (0 files, 0 bytes) and awaiting the owner's deletion
  decision: `shadows/campaigns/2026-08-29T232653Z_feb5c5e5` and
  `.../2026-08-29T232754Z_7bf6c2a3`, both created by my own cancelled
  re-check invocations. No campaign data exists in either; they are
  removable with zero evidence loss, but nothing is deleted without
  explicit approval.

## 2026-08-30 — Gate-B companion at elevated p: beam8 is not "equal", it is ~43× better there

- New companion runner `qldpc-dec/src/gateb_companion.py` (explicitly NOT
  the pinned gate: different p, campaign labelled
  `B-companion-highp`). Rationale: at the pinned p=1e-3 both arms gave 0
  failures in 1e4 shots, and resolving that needs ≈1e6 shots/arm; at
  elevated p failures are resolvable in ~1e2–1e3 shots, which tells us
  whether the equal-accuracy claim's direction is even plausible before
  committing ~71 core-h.
- Smoke (p=3e-3, 200 shots, campaign
  `20260830T124629Z_f989fcb1_33103ff47457`): bp30+osd LER **0.215**
  (43/200) vs beam8 **≈0.005** (1/200) → ratio 0.023, CI95 [0, 0.081],
  decisively BELOW the pinned band [0.87, 1.15]. So in this harness the
  beam8-vs-BP+OSD accuracy relation is **strongly p-dependent**: claimed
  equal at p=1e-3, but beam8 ≈43× better at p=3e-3. A 2000-shot run is
  now in flight for real error bars.
- Scope discipline recorded up front: the elevated-p circuit does not
  exist upstream and is DERIVED by uniform rescale of every noise
  argument of the committed IonQ p=1e-3 circuit (same documented
  technique as the in-repo p=3e-4 point), sha256
  ea2de75c77969fcb…; this companion cannot pass or fail the pinned
  p=1e-3 Gate B and is never to be quoted as if it did.

## 2026-08-30 — arXiv:2608.25545 owner-read: gives qldpc-dec an exact-ML reference standard

- Read first-hand (arXiv API abstract, 2026-08-26; Krishnamoorthy,
  Gerhardt, Knaute, … Piatkowski). Content that bears on this repo:
  degenerate ML decoding is cast as inference in a strictly positive MRF
  over check variables (generalizing the random-bond Ising mapping of the
  surface code to arbitrary CSS codes, including spacetime/circuit-level
  noise). Two decoders: annealed importance sampling that **attaches an
  optimality certificate to each decision** (paired bootstrap; with WISH,
  an exact optimality proof), and a region-based Bethe decoder that
  "reproduces exact ML decoding on every tested surface-code instance at
  millisecond cost" and, with elimination clusters, makes **exact
  degenerate ML decoding of [[72,12,6]] feasible**. Evaluated on exactly
  our objects: surface codes plus BB [[72,12,6]] and [[144,12,12]] under
  code-capacity, phenomenological, and circuit-level noise. `[REPORTED]`
  at abstract level — full text not yet read.
- Why it changes this target's ceiling: `qldpc-dec/` currently measures
  heuristics against each other (BP+OSD vs beam vs NMS). An exact-ML
  reference on [[72,12,6]] converts that into an **optimality-gap**
  measurement — how far from ML is each fast decoder, and on which
  syndromes — which is the fair-baseline question the target exists for.
  It also supplies a per-shot distrust flag, the natural instrument for
  the tail analysis this repo already ran on the OSD mean gap.
- **Proposed Gate D (not yet pinned, not yet run):** on [[72,12,6]] under
  code-capacity noise, measure the decision-level optimality gap of
  bp30+osd against a certified-ML reference, reporting (i) fraction of
  syndromes where BP+OSD's logical class differs from ML, (ii) the LER
  penalty that difference costs, (iii) whether the disagreement set
  concentrates on the certificate-flagged syndromes. Requires the full
  paper read + a decision on reimplementing the region-based decoder vs
  a brute-force coset-sum reference at small n; costed before any launch.

## 2026-08-30 — na-compiler optimal-COST rows verified (and a Main verifier bug)

- NaCompilerExt returned what my freeze was missing: optimal **cost**
  rows with explicit witnesses (my earlier campaign had batch counts /
  diameters only). Rows: n=5 target 0-2-1-4-3 → **539.043122398901 µs**
  in 4 batches; n=6 target 1-0-4-8-6-3 → **666.354655402228 µs** in 5
  batches. It also removed its own duplicate campaign dirs on request and
  left the owner-frozen campaign untouched.
- Owner verification (`na-compiler/src/verify_cost_rows.py`, independent
  of the agent's code path, using the FROZEN primitives
  `conflict.compatible` + `config.move_batch_duration_us` + 2×15 µs
  transfers): every batch legal, replay reproduces each intermediate
  witness state and lands exactly on target, and cost reproduces to
  **2.3e-13 / 4.5e-13 µs**. Batch counts equal the diameters I had
  already established by exhaustive BFS (4 and 5) ⇒ the witnesses are
  batch-optimal, and their cost is the pinned-model cost of an optimal
  schedule. Frozen:
  `na-compiler/campaigns/20260830T124308Z_409b4d89_2afb5e487a3b/`.
- **My error, recorded:** the first verifier run reported both rows as
  illegal/mismatched. The fault was mine — I assumed move pairs meant
  (source_site, target_site); they are (atom_index, target_site), which
  exact witness replay settles unambiguously. Under the correct semantics
  the agent's rows are exact. Had I trusted my own first output I would
  have wrongly accused a correct result; the replay-against-frozen-
  primitives discipline is what caught it.
- Scope note kept honest: cost-OPTIMALITY of the value (no cheaper
  4-/5-batch schedule) rests on the agent's consistent-A* argument, spot-
  validated by exact Dijkstra on three n=3 targets — labelled [DERIVED],
  not owner-proved. Owner-proved here: legality, replay, cost arithmetic,
  and batch-count optimality.

## 2026-08-30 — deadline/liveness audit of the long campaigns (advisory check)

- Concern raised: a bash launch without `timeout: 0` could be killed by a
  default deadline mid-run and then be mis-recorded as "running".
  Audited all four long launches — each carried `timeout: 0`, and the
  settled wall times prove no deadline truncated them: Gate A accuracy
  22,737 s (6.3 h), Gate A timing 1,056 s, Gate B pass 1 5,962 s (died on
  the `ratio_ci` bug, not a deadline — full traceback captured), Gate B
  relaunch 6,018 s. All four exceed any default cutoff by design.
- Liveness discipline was applied before every completion claim, matching
  the bar used on the FssBb handoff: pid + argv checked with `ps`,
  campaign dir growth checked on disk, and the ledger claim written only
  after reading the frozen `summary.json`/`inventory.json` (or, for
  fss-bb, `thresholds.json` + `fss_fit.json`). One false alarm during
  Gate A (no `results.json.gz` after 18 min) was diagnosed as expected
  behaviour — the runner appends per POINT, and a point takes ≈32 min —
  and recorded rather than acted on.
- Current state at audit time: zero campaign processes alive
  (`ps` count 0), zero "RUNNING/in flight" claims left in RESULTS.md;
  every row points at a frozen artifact. The only live peer is the
  `NaCompilerExt` agent, told to stand down to avoid duplicate campaign
  dirs after its work was owner-verified and frozen.

## 2026-08-30 — msd verdict RECLASSIFIED: surrogate frame is distance-1, not a refutation

- Owner diagnostic (single core, 1e6–2e6 shots/point, seed 20260830): the
  surrogate's three logical observables each fail at ≈1.6–1.7% at
  p=1e-4 (obs0 0.01668, obs1 0.01671, obs2 0.01613; any-obs pL 0.03834 —
  independently reproducing the frozen prime campaign's 0.038356).
- **Exponent measured directly across two decades:** per-observable
  failure vs p has log-log slope **0.992** with a stable coefficient
  (obs0/p = 166.9, 167.8, 166.0, 162.3 at p = 1e-5, 3e-5, 1e-4, 3e-4).
  A distance-1 (unencoded) frame responds linearly; a distillation
  circuit with the paper's claimed suppression responds quadratically.
  The committed reconstruction is therefore **structurally incapable** of
  exhibiting c·p² at any p — its unexpanded d=3/d=7 output patches (an
  [INFERENCE] simplification declared in the target README from the
  start) leave the output frame bare.
- **Verdict reclassified** accordingly: SURROGATE_REFUTATION_AT_LOW_P
  understates and mis-frames what was measured. The correct statement is
  **SCOPE-LIMITED-NON-TEST**: the reconstruction measures its own
  distance-1 floor (≈167·p), not the paper's protocol, so arXiv:2605.21867's
  c≈300 is neither supported nor refuted by this target as currently
  built. Failures are also positively correlated across observables
  (P(no fail) 0.96166 observed vs 0.95129 independent-model), consistent
  with common-cause faults in a shared bare frame.
- Concrete unblock, now the target's named next action: expand the
  d=3/d=7 surface-code patches on the three output representatives so
  the frame carries distance; only then is the c·p² window reachable.
  Everything else in the target (encoder from Fig. 4a, postselection
  estimator, campaign machinery) is reusable as-is.

## 2026-08-30 — na-compiler Gate B blocker CLOSED (owner-verified)

- The 30k-state completeness cap that blocked honest 4×4 rows is gone.
  New module `na-compiler/src/na_compiler/idastar.py` (admissible
  Manhattan/D_max heuristic + IDA* + bidirectional search + a compiled
  bitstate BFS enumerator) closes both previously capped instances by
  FULL state-space exhaustion, not sampling:
  - n=5 on 4×4: 524,160 states visited = P(16,5) exactly, optimal batch
    diameter **4**, frontiers [1, 341, 29847, 353858, 140113, 0].
  - n=6 on 4×4: 5,765,760 states = P(16,6) exactly, diameter **5**,
    frontiers [1, 499, 68670, 2042809, 3641150, 12631, 0], 268 s.
  Owner reran both from scratch and reproduced every number digit-for-
  digit before freezing; campaign
  `na-compiler/campaigns/20260830T123706Z_4ed93f70_9bc5ba162924/`
  (manifest + results + the exact C source + inventory).
- Revision 1 of `na-compiler/pre_statement.md` was appended pre-run by
  the agent and states the heuristic, its admissibility/consistency
  argument, and what a COMPLETE row certifies (witness batch sequence
  verified against the frozen legality code + exhaustive UNSAT proof at
  d*−1, pruning only by the admissible heuristic).
- Process note: the owning agent (`NaCompilerExt`) delivered a stale stub
  as its final artifact and died before freezing; its work was real and
  verifiable on disk, so Main verified independently and froze the
  campaign. Artifact-vs-claim divergence recorded rather than papered
  over.

## 2026-08-30 — source audit carried over from `cs/`: three unknown sources, one load-bearing

Forwarded by the owner while auditing `../cs/`'s frontier-scan coverage. Not a `physics/`
campaign entry and nothing here is claimed as verified — this is a **source note** so the
relevant targets do not rediscover it later.

- **ECCC TR26-160** (19 Aug 2026), Louis Golowich, Itzhak Tamo, Guanyu Zhu, *"Improved
  Transversal Non-Clifford Gates from Cup Products"*. Abstract, as read from the ECCC listing
  (**REPORTED** — the full report is not yet owner-read): constructs quantum codes with
  low-weight stabilizers supporting transversal, i.e. low-depth, implementations of the
  non-Clifford $C^{r-1}Z$ gate **for every constant $r\ge3$**, in particular obtaining
  length-$n$ quantum LDPC codes with that property.
- **Why it matters here.** Two live targets are directly adjacent: `msd/` (zero-level CCZ
  reproduction against arXiv:2605.21867, currently `SURROGATE_REFUTATION`) and `qldpc-dec/`
  (BB-code decoder reproduction harness). A low-overhead transversal non-Clifford route is the
  competing mechanism to magic-state distillation, so it bears on `msd/`'s framing, and
  "quantum LDPC codes with low-weight stabilizers" is `qldpc-dec/`'s object.
- **How it was found, and the transferable lesson.** `../cs/`'s scan claimed ECCC coverage as
  `TR26-001..154`. ECCC numbers are **not chronological** — 11 of 19 adjacent pairs in
  TR26-142..161 are date-inverted, and TR26-161 is dated 12 Aug while TR26-154 is dated 28 Aug.
  A numeric ceiling silently excluded seven reports dated 12–28 August. **Applies to this repo
  too: state source coverage by DATE against the dated listing, never by report number.**
- **No target opened, no gate changed.** Owner read required before this influences any
  `msd/` or `qldpc-dec/` claim.

**Same-session source audit, carried over from `../cs/`.** All 20 arXiv ids referenced anywhere
in `physics/` were version-checked against the arXiv API: 11 sit at v2+, but **none was revised
after 2026-08-29**, so no claim in this repo rests on a superseded version. `msd/`'s refutation
target `2605.21867` is still v1. A competitive-field sweep then found two papers with **zero
prior mentions in `physics/`**:

- **`arXiv:2608.25545` (2026-08-26), *"Certified decoding of quantum LDPC codes"* —
  LOAD-BEARING for `qldpc-dec/`, owner-read at abstract level (`REPORTED`).** It works on
  **exactly this target's objects**: the bivariate bicycle codes **[[72,12,6]] and
  [[144,12,12]]**, under code-capacity, phenomenological, and circuit-level noise. Four points
  that matter here. (i) It states flatly that *"the workhorse decoder BP+OSD sidesteps degeneracy
  heuristically and offers no guarantees"* — and BP+OSD replay is `qldpc-dec/`'s harness. (ii) It
  attaches *"to every decision a certificate of optimality: a paired bootstrap test, or, composed
  with constant-factor estimators such as WISH, an exact optimality proof."* (iii) Region-based
  Bethe free energy *"reproduces exact ML decoding on every tested surface-code instance at
  millisecond cost"*, and elimination clusters make *"exact degenerate ML decoding of the
  [[72,12,6]] bivariate bicycle code feasible"* — i.e. a **reference standard** against which a
  reproduction can be validated rather than merely replayed. (iv) *"the certificate flags exactly
  the syndromes on which any fast decoder should be distrusted"*, which is a per-shot instrument
  for precisely the kind of tail analysis this repo just did on the OSD mean gap. Directly
  relevant to the entry below: an OSD-postprocessing protocol mismatch is a symptom of BP+OSD
  having no degeneracy guarantee, which is this paper's thesis.
- **`arXiv:2606.24170` (2026-06-23), *"Low Spatial Cost CCZ Magic State Factory"* — adjacent to
  `msd/`, not colliding.** Reconstructs an eight-to-three CCZ distillation protocol as compact
  joint-measurement surface-code architecture, preserving single-fault detection and
  leading-order suppression *"with lower spatial cost than the design of Gidney and Fowler"*.
  `msd/`'s target is **zero-level** distillation (Itogawa et al.), a different mechanism, so this
  is design-space work next door rather than a competing claim.

**Method note, and why it is worth repeating here.** This sweep is the same technique that in
`../cs/` turned one missed paper into five unknown works, one of which used that repo's own
instrument on its own object. The two finds above were reached the same way: version-check every
cited source, then keyword-sweep the target's object and check each hit against the repo.
**Neither find changes any `physics/` verdict, and neither may, until owner-read in full.**

---

## 2026-08-30 — shadows gate-C k-extension: edge VANISHES at each enumerated k=9..12

- ShadowsExt completed the Revision-2 campaign
  (`shadows/campaigns/2026-08-30T114124Z_ba4f8430_65d2435684e3/`, 23.8 s,
  exact Fractions): regression gate 35/35 exact matches versus the frozen
  k≤8 artifact — PASS; k=9..12 are all **EDGE-VANISHES**, with w_cx<w_ct
  exactly. From k=9 to 12, the absolute weight gap narrows while w_cx/w_ct
  decreases, so the relative/shadow-norm disadvantage widens. Eq. (4)
  equals enumeration at every completed k. The pure `cx_chain` k=4
  improvement is absent at every enumerated k=6..12; this is a bounded
  finite-k result, not an asymptotic claim.
  Contract JSON copied at shadows/src/ext_k_contract.json.
- Process notes recorded honestly: the agent first yielded empty (caught,
  restarted); then caught a direction inversion in MY context contract
  before freezing (same bug class as the target's Revision 1) — pinned
  (a) PERSISTS iff w_M(k) > w_ct(k) per Revision-1 semantics; also
  frozen the carrier-identity scope note (pure cx_chain is the k=4
  improvement carrier; literal cx_chain∘U_ct composition is a different
  map, w=1/81 at k=4, excluded).

## 2026-08-30 — historical scalar-prior OSD diagnostic (attribution retracted)

- This 3,000-decode diagnostic used the scalar-mean BP prior later invalidated
  by Revision GB3. Under that wrong prior, OSD0 had mean 3.54 ms, pure BP had
  mean 3.36 ms, and OSD-CS10 had mean 65.8 ms with 28% of shots in the
  approximately 150–380 ms fallback slice.
- The measurements remain useful evidence that the scalar-prior decoder
  exercised the costly CS10 path far more often than OSD0. The original
  inference that the paper's Table-II mean was an OSD0 or OSD-skipping number
  is **RETRACTED**:
  corrected heterogeneous-prior OSD-CS10 subsequently reproduced both means
  and p99.9 tails, while corrected-prior OSD0 fell far below every published
  tail. See the 2026-08-31 correction at the top of this ledger.

## 2026-08-30 — historical scalar-prior bimodality (attribution retracted)

- Four thousand timed decodes under the scalar-mean prior later invalidated by
  GB3 were sharply bimodal: median 3.19 ms; 71.7% below 5.9 ms with mean
  2.6 ms; 28.3% in an approximately 150–380 ms fallback mode with mean
  225.6 ms. The p99.9 of 320.8 ms matched the invalidated timing campaign's
  326.7 ms tail.
- Those measurements isolate the behavior of the wrong-prior run: it
  exhibited frequent expensive OSD fallback. They do not establish that the
  paper skipped OSD. The former OSD-skipping attribution is
  **RETRACTED** because corrected-prior OSD-CS10 reproduces the published
  mean and tail while corrected-prior OSD0 does not. The live
  `timing_forensics.py` now uses the corrected vector prior; the numbers above
  are retained only as the historical scalar-prior record.

## 2026-08-30 — qldpc-dec Gate B first pass: UNDERDETERMINED at 1e4 shots

- Beam8-vs-BP+OSD at p=1e-3 Z, 1e4 shots/arm (frozen
  `campaigns/20260830T092524Z_c3272f02_91994ec12b0f/`): **0 failures on
  both arms** — both decoders below the resolvable LER at this budget, so
  the [0.87, 1.15] equal-accuracy band check is UNDERDETERMINED, not
  passed. Verdict string records the outcome; the target's own Gate-B
  falsification requires resolvable arms.
  - Enabled by two code fixes made en route (bg_2 crashed on the
    same-statistics point): `ratio_ci` degeneracy guard (0/0 → full
    support rather than IndexError) and an explicit
    UNDERDETERMINED branch in cmd_gate_b; smoke + compile checks ran
    before relaunch. The crash itself was informative: BP+OSD LER at
    p=1e-3 Z on the committed DEMs sits below 1e-4 with 1e5 shots —
    consistent with the accuracy campaign (4e-5 at 1e5).
- Next (budget-sized before launch): the paper's equal-accuracy claim
  lives at LER ≈ 2e-4 (Fig. 2); resolving the band at 95% needs enough
  failures that the ratio CI is narrow — with both arms ≈0 failures at
  1e4, the required shot count is ≥1e6/arm (~53 core-h beam8 numpy +
  ~18 core-h bp30+osd). Parked as a sized decision, not silently
  launched, because it exceeds the single-session core budget.

## 2026-08-30 — fss-bb 20k ladder FROZEN; paper-CI overlap on 4 of 5 sizes

- Detached run completed and froze as
  `fss-bb/campaigns/20260829T235018Z_6df86846_e8697f77b182/`
  (thresholds.json, fss_fit.json, results.csv 41 rows, manifest; dem/
  kept). Process 91617 exited cleanly.
- Threshold ladder vs paper (arXiv:2603.19062 Table 1): 12×6 p*=0.3707
  CI [0.36985, 0.37064] OVERLAPS paper 0.3701 [0.3697, 0.3703];
  18×9 0.4387 vs 0.4386 (tight CI, consistent); 24×12 0.4456 [0.4452,
  0.4463] vs paper 0.4608 — CI-MISS at 20k shots ⇒ INCONCLUSIVE-UNDER-
  BUDGET per Amendment A1 (not a refutation); 30×15 0.4675 vs 0.4702;
  36×18 0.4706 [0.4705, 0.4708] vs 0.4706 — matches to 3 decimals.
- Point-estimate anomaly resolved: the pre-freeze quote "0.3707 CI
  (0.3699, 0.3706)" was a rounded-underscore artifact; the frozen file's
  bracket [0.37032, 0.37070] shows the point is the Illinois upper edge
  and the bootstrap CI is the authoritative interval — self-consistency
  hold against `thresholds.json` passes.
- FSS fit (all windows w0.04/0.06/0.08 agree): p*_inf = 0.4980 CI
  [0.4973, 0.4988], ν = 1.412 CI [1.391, 1.434]. Paper: 0.488±0.001,
  ν=1.18±0.01. **ν mismatch lives at 20k budget:** 1.41 vs 1.18 is
  outside the paper's CI, but per A1 with 20k shots (3.2× wider CIs) and
  our CI half-width ≈0.02, this is INCONCLUSIVE-UNDER-BUDGET on the ν
  question; the 200k-shot full ladder (221 core-h, parked) is the only
  repo route to a real ν verdict. Falsifiable statement of record: under
  this budget, no claim either way on ν=1.18-vs-4/3.

## 2026-08-30 — Gate B beam-search campaign launched

- `run_gate gateB 10000` launched as bg_2 (arm beam8 vs bp30+osd baseline;
  1e4 shots/point per the measured beam8 526 ms/shot numpy cost ⇒ ≈1.5-h
  first point; full ladder sized at launch). Checks: beam8 LER ratio band
  [0.87, 1.15] vs bp30+osd at p=1e-3; ratio CIs from ratio-wilson. This
  campaign targets the beam8-equal-accuracy claim only; beam32/64
  ratio-band checks need ≈30+ core-h each and are queued behind an
  explicit compute window (numpy port runs ~8× slower than the paper's
  compiled decoder, so ratio-verdicts on LER remain fair while ms
  comparisons stay excluded).

## 2026-08-30 — qldpc-dec Gate A TIMING arm closed: split verdict

- Timing campaign `qldpc-dec/campaigns/20260830T071750Z_00327d68_a05388ffcf42/`
  (10,000 single-decode timings per point, Z basis): p=1e-3 mean 67.8 ms /
  p99.9 326.7 ms; p=5e-4 mean 37.7 ms / p99.9 315.5 ms. Against arXiv:
  2512.07057 Table II/III (ldpc-stack config, M3): p999 ratio 1.13× and in
  the factor-2 band ≈ REPRODUCED; mean ratio 6.4× — OUTSIDE factor 2 ⇒
  `within_factor2: false` at p=1e-3.
- Reading (honest, scope-limited): the tail distribution reproduces the
  published anchor; the mean does not. Descriptive, not gate-breaking by
  itself: pre_statement pinned factor-2 per-quantity, so mean-ms at 1e-3
  is a REPRODUCTION MISS on that one quantity. Candidate causes for the
  mean gap (ldpc-package version/threads, batched-vs-single decode
  protocol, on-node contention with the FSS ladder) are named hypotheses
  for the next diff session — none confirmed; do not promote either way
  without them.
- RESULTS row updated: Gate A timing = tail reproduced, mean missed;
  beam-search (Gate B) and GARI (Gate C) campaigns remain the next
  campaigns under the same harness.

## 2026-08-30 — qldpc-dec Gate A campaign CLOSED and verified

- Full 6-point accuracy companion (BP+OSD bp30+osd/osd_cs-10, 1e5 shots
  per basis-point, seeds logged per point) completed in ≈6.4 core-h:
  LER = 0 at 3e-4 and 5e-4 (both bases, CI95 upper 3.84e-5), 7e-5 X /
  4e-5 Z at 1e-3 (total-XZ ≈ 1.1e-4). Campaign frozen at
  `qldpc-dec/campaigns/20260830T005616Z_3402651f_5b4e76ede0ee/`
  (manifest/results/summary/inventory; process exit 0).
  Bookkeeping note: the first manifest written at 00:56 belonged to the
  aborted hub-broker launch; the OS reassigned that process's campaign
  object to the surviving run, so the FROZEN directory is `...005616Z...`
  and `...005856Z...` holds the partial pre-freeze snapshot — both
  consistent with the ledger trail above.
- Interpretation against the pre-statement: Gate A accuracy was declared
  a companion, NOT a reproduction criterion (no published LER curve
  exists); the numbers are consistent with the decoder operating in its
  nominal regime on the committed DEMs (single-obs decodable at 3e-4/5e-4
  within 1e5 shots; failures appearing at 1e-3, heavier in X). No
  universal threshold claim — these are point estimates with logged
  seeds. Timing arm (`gateA_timing`) remains the actual Gate A
  reproduction criterion and is queued next.

## 2026-08-29 — qldpc-dec Gate A campaign launched (owner instruction "continue")

- User directed continuation of all parked work; Gate A accuracy
  (BP+OSD, 1e5 shots × {3e-4, 5e-4, 1e-3} × {X, Z}, seeds base 20260829)
  is now running as bg_1 under hub name `qldpcec-gatea` (pid 97581,
  campaign `campaigns/20260830T005856Z_48dbfa91_5b4e76ede0ee`). An earlier
  hub-broker launch aborted mid-request, leaving residue dir
  `...005616Z_3402651f...` (manifest only, no process, empty results) —
  quarantine candidate, do not treat as a campaign.
- False-alarm analysis, recorded to prevent repeat: ~18 min in there is no
  `results.json.gz` yet — by design. `run_point` decodes each point's full
  1e5-shot batch in one harness call and appends only at point end; at the
  smoke-measured 64.1 ms/shot a point takes ≈107 min, so the 6-point
  ladder legitimately runs ≈10.7 single-core hours (matches QldpcDec's
  feasibility estimate). RSS stable, sampler phases cloud-free.
- Expected completion: ≈11:30 UTC next day. On completion: point-wise
  LERs vs pre_statement bands (gate A honesty check), timing rows vs
  IonQ Table II/III anchors, then RESULTS row promotion.

## 2026-08-29 — fss-bb handoff: detached 20k ladder verified running

- FssBb hit its request cap and handed off cleanly. Main verified the
  claim rather than trusting it: `hub ps` shows `fssbb-gateA-20k`
  pid 91617 detached (uptime 1h01m at check), `ps -p 91617` confirms the
  exact argv (`campaign_runner.py --sizes 12x6,18x9,24x12,30x15,36x18
  --shots 20000 --tag gateA-20k --seed-base 20260829 --n-boot-thresh 5000
  --n-boot-fss 500`), and `campaigns/20260829T235018Z_6df86846/dem/`
  exists on disk. Freeze protocol: runner exits and renames the dir to
  `<UTC>_<uuid>_<sha12>` with manifest/results/thresholds/fss_fit/dem.
- Amendment A1 (recorded in fss-bb/pre_statement.md BEFORE data use): 20k
  CIs are ~3.2× wider than the paper's 200k CIs, so a CI miss is
  INCONCLUSIVE-UNDER-BUDGET, never a refutation. Main endorses this rule.
- **Open check for freeze time (Main flag):** FssBb's pre-freeze message
  quotes 12×6 as p* = 0.3707 with CI (0.3699, 0.3706) — the point estimate
  lies outside its own interval, so at least one of the two numbers is
  mis-transcribed. Do not promote that row until `thresholds.json` is read
  from the frozen campaign and the estimate/CI pair is checked for
  self-consistency.
- Also delivered: FSS module validated against synthetic data with known
  exponents (T1 exact recovery, T2 out-of-class bias <1%), two pre-sweep
  bugs found and disclosed (sector-failure pairing; degenerate single-init
  fit → multi-start), feasibility bench (decode ~N², 0.69→441 ms/shot).

## 2026-08-29 — msd Gate-A-prime closed; Main steering-arithmetic error corrected

- Msd ran Gate-A-prime (1e7 shots × {1e-4, 3e-4, 1e-3}, artifact
  20260829T234229Z_ce1346fd_40d024611d98): p_L = 0.038356 / 0.110333 /
  0.319565 [NUMERICAL]. Owner spot rerun at p=1e-4, 2e5 shots: p_L
  0.038218, z = −0.31 vs the worker point — statistically identical.
- Main's steering text for the prime run mis-multipled the paper's
  implied p_L (wrote 3e-3 for c=300·p² at p=1e-4; correct is 3e-6). Msd
  caught the error before any promotion used the wrong band; measured gaps
  vs the paper curve are 1.28e4× / 4.09e3× / 1.07e3× at the three p
  values, with falling c_eff (3.8e6 → 1.2e6 → 3.2e5) ⇒ measured exponent
  ≈ p⁰·⁹², floor-like, not quadratic at any reachable p. Verdict string
  SURROGATE_REFUTATION_AT_LOW_P unchanged; target band corrected by the
  implementation agent in pre_statement rev2 (original frozen statement
  and rev1 untouched). The error is Main's; the catch is Msd's.

## 2026-08-29 — na-compiler gates delivered (owner-verified in part)

- `NaCompiler` completed: QMAP C++ sources + RSQASM read raw from GitHub;
  constraint model pinned from `IndependentSetRouter::isCompatibleMovement`
  + QMAP's own evaluator. Gate A campaign certified qft3/qft4/reg3_4
  optima (2888.28/12369.39/5301.13 µs) vs MQT QMAP routing-aware placement
  (gaps +0.9188/+0.4444/+0.6645, Z3 UNSAT@k−1 + SAT@k + HiGHS dual bound
  per instance). Owner rerun of the qft3 row (0.42 s) matched the frozen
  table. Gate B: 3×3 BFS complete through P(9,4)=3024; 4×4 frontier
  requires bitstate/IDA* before honest complete rows (documented
  blocker). During integration Main noticed the na-compiler row had been
  accidentally dropped from RESULTS.md by an earlier line-numbered cut
  (the qldpc-dec/qlops shadow-cast); row restored with the same-line edit
  and verified against the on-disk table — incident recorded here, no data
  loss (the frozen campaigns were untouched).

## 2026-08-29 — fss-bb 20k ladder incident (self-reported) + relaunch

- First 20k-shot ladder run (bg_48) died at the 24×12 stage: the FssBb
  agent archived stale campaign residue INTO campaigns-smoke/crashed/
  while the runner was live, breaking a dem/ path the runner assumed
  stable — an operational error, self-reported, not a physics result.
- Fixes adopted: runner recreates dem/ per size (robust to external dir
  moves); campaign dirs strictly hands-off until freeze; fresh relaunch
  bg_3 with unchanged seed lineage (identical results modulo wall time),
  ETA ~9–10 h from 22:59 UTC.
- Pre-freeze NUMERICAL signal so far: 12×6 threshold 0.3707 vs paper
  0.3701; 18×9 0.4387 vs 0.4386 — both agreeing. Unchanged budget ruling:
  20k variant now; 200k fullbudget ladder remains parked pending explicit
  user approval.

## 2026-08-29 — msd Gate-A v1 closed (surrogate grid = saturation regime)

- `Msd` completed the zero-level-CCZ Gate-A v1: figure-derived 22-wire
  Clifford surrogate ([8,3,2] encoder transcribed from Fig. 4a + identity
  T/T† surrogate + postselection estimator; ideal 256-shot static check
  passes exactly), frozen 4-point × 1e5-shot grid run and frozen
  (campaign ...d9fe3c0d...; owner rerun ...9e9627ec... reproduced the fit
  identically: slope 0.5531, c = 13.578, accepted-shot p_L 0.18→0.93).
- Main scope correction: the frozen grid sits in the saturation regime
  (slope ≠ 2), so the paper's asymptotic p_L ≈ 300 p² claim is UNTESTED by
  v1, not refuted. The pre-statement verdict string
  REFUTE_FOR_THIS_RECONSTRUCTION is retained as the frozen-protocol
  verdict, now annotated with the regime finding.
- Gate-A-prime requested via hub (pre_statement Revision 1): p ∈
  {1e-4, 3e-4, 1e-3} × 1e7 shots/point directly probes the quadratic
  window; verdict either branch (surrogate floor vs ≈3e-3 at p=1e-4) is a
  clean first RESULT-grade datum for this target.

## 2026-08-29 — shadows gates closed (owner-verified); FssBb ladder launched

- `Shadows` closed all three gates with exact arithmetic (194 s, no RNG):
  Gate A exhaustive census (11,520 labels) confirms the npj-QI lemma and
  sharpens it — every two-qubit Clifford contracts exactly 0 or 4 of the 9
  size-2 Paulis (bimodal histogram, nothing in between); Gate B re-derives
  Eq. (3)–(7) exactly (k≤10 exhaustive, k≤16 integer-exact) plus the
  1/(2^k+1) baseline; Gate C finds a scope-limited exact improvement: the
  pure `cx_chain` and `cx_mesh_lex` maps have weight 361/6561 versus
  353/6561 for `U_ct` at k=4 (shadow-norm ratio 353/361, ~2.2% better at
  that single k). The literal `cx_chain`∘`U_ct` composite is a different
  excluded map with weight 1/81 at k=4; no asymptotic claim is made. Owner
  rerun of `run_campaign.py` reproduced the artifacts (193.0 s).
  Pre-statement had been frozen before implementation;
  one handedness revision is documented (Revision 1), never silently
  edited.
- `FssBb` reproduced the N=144 paper-budget fragment (p* = 0.37032,
  CI [0.37004, 0.37101] vs paper 0.3701 [0.3697, 0.3703]; CIs overlap) and
  independently re-derived p*_inf = 0.4896 on their own ladder
  [DERIVED]. Full 5-size 200k-shot Gate A would cost ≈221 core-hours and
  stays parked pending explicit user approval; on Main's instruction the
  20k-shot fallback ladder (≈6–12 core-hours, all 5 sizes) is now running
  as FssBb background job bg_48, seeds pinned per-point, results due
  later tonight.
- RESULTS.md updated with verdict rows for qldpc-dec (ACTIVE, smoke),
  shadows (A+B+C closed), qlops (A+B closed), fss-bb (fragment + ladder
  running); campaign-ladder section now lists all frozen snapshots.

## 2026-08-29 — first worker verdicts integrated

- `Qlops` completed Gate A and Gate B. Owner rerun of
  `qlops/src/gate_a_reproduce.py` and `gate_b_sweep.py` passed (both exit 0):
  81 Gate-A checks, 0 hard failures; QLOPS Table-5 rows reproduce to ≤4e-8
  after the paper's Table-2a `t_r` one-slot correction. The printed table
  has a 33–43% `t_r` row-shift error for d≥17; the correction is recorded,
  not silently applied. Latency sensitivity max = 1.82× (<2×), while the
  7-T ↔ zero-level-CCZ magic-state swap gives 193×–29,944× spread and
  falsifies comparability on that axis. Artifacts:
  `qlops/campaigns/20260829T223540Z_dcf693ab/`.
- `QldpcDec` completed the receiver-side BB decoder harness and smoke
  validation. Gate A was corrected after first-hand reads: no published
  BP+OSD LER curve exists at the originally assumed p values; the published
  gate quantities are runtime (arXiv:2512.07057 Table II/III). The smoke
  timing anchor reproduced p99.9 = 272.2 ms vs 272.5 ms (0.1%); BP+OSD,
  beam-search, and NMS arms decode end-to-end. Three smoke campaigns and
  frozen DEM provenance are in `qldpc-dec/`; full campaigns remain pending
  because measured costs are ≈10.7 core-hours for Gate A and ≈15–30 hours
  for the numpy beam arm.
- These are **smoke/benchmark verdicts**, not universal QEC theorems.
  RESULTS.md now records both rows; the five remaining workers/long runs
  stay in progress.

## 2026-08-29 — six targets opened; six implementation agents launched

- All six shortlist targets opened with per-target READMEs (gates + exact
  pass/fail criteria + pre-campaign requirements): `qldpc-dec/`, `msd/`,
  `shadows/`, `na-compiler/`, `fss-bb/`, `qlops/`. RESULTS.md now carries
  one row per target (all `BENCHMARK (launching)`, no results yet).
- Six owning implementation agents launched in parallel (`QldpcDec`,
  `Msd`, `Shadows`, `NaCompiler`, `FssBb`, `Qlops`), each with a strict
  single-folder ownership contract: they may edit only their own target
  folder; this session (Main) remains the single writer of README.md,
  RESULTS.md, PROGRESS.md (SSOT consistency across targets).
  Shared-infrastructure rule: FssBb must coordinate the decoder-harness
  interface with QldpcDec via hub messaging; if unresolved, it builds a
  local minimal path inside its own folder.
- Toolchain ready before launch: `physics/.venv` (stim, ldpc, numpy,
  scipy, highspy, z3-solver), smoke-imports pass.
- Deliverable convention for this launch wave: per-target pre_statement.md
  + working src/ + one end-to-end smoke artifact each (Shadows' gate-A
  census and Qlops' gate-A reproduction are expected to produce verdict-
  grade evidence already).

## 2026-08-29 — date-precision correction (advisory, retraction)


- The Published/Event cleanup had introduced unverified day-precision dates.
  Verified against primaries today (arXiv API + nature.com): arXiv:2511.09551
  published **2025-11-12** (was shown as 2025-11-14 — wrong, retracted);
  Quantum Echoes Nature (10.1038/s41586-025-09526-6, "Observation of
  constructive interference at the edge of quantum ergodicity", Nature 646,
  825–830) published **2025-10-22** (the earlier 2025-06-25 in this ledger
  was wrong — retracted); arXiv:2508.20699 published **2025-08-28**
  (confirmed). `physics/docs/QC_FRONTIER_2026-08-29.md` is now
  month-precision throughout, with day precision retained only where
  verified against the source (2508.20699 and the two above).
- The 13,000× Quantum-Echoes figure was pinned down first-hand from the
  Nature page: tensor-network contraction ≈3.2 yr vs 2.1 h experimental
  collection on Frontier — a collection-time ratio, not a factor over all
  classical simulation methods.
- Note: the older drift-fix entry below still cites the pre-correction dates
  (2025-06-25 / 2025-11-14); those numbers there are superseded by this
  entry and were wrong at the time — the drift-fix entry's stated edge
  (pub vs event separation) was correct; its cited days were not.

## 2026-08-29 — five-scout parallel frontier scan

- Launched 5 read-only scout subagents (qEC, neutral-atom compilation,
  classical shadows, BP+OSD/qLDPC decoding, beyond-the-four). All
  completed; JSON outputs kept as session artifacts. New ranked shortlist
  with falsifiable gates written to
  [`docs/QC_SCOUT_SCAN_2026-08-29.md`](docs/QC_SCOUT_SCAN_2026-08-29.md);
  the README candidate table updated to match.
- Key convergence: QecScout and DecoderScout independently nominated the
  same target — a neutral cross-paper decoder reproduction harness for BB
  codes ([[144,12,12]]): Relay-BP (arXiv:2506.01779), beam search
  (arXiv:2512.07057), GARI (arXiv:2510.14060) all publish record claims on
  the same code with private DEMs, and arXiv:2603.19062 just showed
  independent reproduction methodology is itself publishable (erasure
  channel).
- New shortlist: 1 `qldpc-dec/`, 2 `msd/` (zero-level CCZ gate, hours on
  one core), 3 `shadows/` (11,520-Clifford exhaustive census),
  4 `na-compiler/` (certified-optimal ILP baselines), 5 `fss-bb/` (ν
  anomaly), 6 `qlops/` (resource-estimate audit). `sdec/` folded into
  `qldpc-dec/`; OTOC(2) census rejected as too crowded; `qmaqcma/`,
  `scars/`, `mbl2d/` parked from v1.
- Evidence caveat recorded: scout VERIFIED-PRIMARY tags are subagent
  verdicts, not owner reads; every gate statement requires a first-hand
  primary read before a campaign starts.
- Toolchain prepared: `physics/.venv` with stim, ldpc (BP+OSD), numpy,
  scipy, highspy, z3-solver; smoke-import verified.

## 2026-08-29 — naming/date drift fixed (advisory)

- Folder names in `README.md` (`qec/`, `decoder/`, `qma/`, `mbps/`) disagreed
  with `docs/QC_FRONTIER_2026-08-29.md` (`qldpc-dec/`, `sdec/`, `qmaqcma/`,
  `scars/`, `mbl2d/`) — the exact naming-drift class `../math/README.md`
  documents as a past incident. Both now use one canonical set; README examples
  updated. The scan doc's date column mixed publication dates with event
  dates; separate **Published** / **Event** columns were introduced.
  [Correction 2026-08-29, later the same day: the specific day values cited
  in this entry (Nature Echoes 2025-06-25; STOC arXiv 2025-11-14; Willow
  2024-12-09) were partly wrong — see the date-precision correction entry
  above; the Structural fix stands.]

## 2026-08-29 — physics/ folder created

- Created `physics/` per user request (PhD focus: quantum computing; goal:
  push the frontier on open physics/QC problems).
- Wrote tracking contract mirroring `../math/README.md` (SSOT: RESULTS.md +
  PROGRESS.md + per-target READMEs + immutable campaigns).
- Ran a frontier scan over 2025–2026 quantum-computing developments (qLDPC
  experiments, Willow/Quantum Echoes, QCMA oracle separation, magic-state
  distillation, many-body scars / 2D MBL); primaries for the headline claims
  read and verified: arXiv:2604.15427 (TNBP cannot simulate Quantum Echoes),
  arXiv:2511.09551 (STOC 2026 QMA/QCMA oracle separation, via Quanta writeup),
  arXiv:2602.09385 (second oracle separation). Scan output:
  [`docs/QC_FRONTIER_2026-08-29.md`](docs/QC_FRONTIER_2026-08-29.md).
- No target opened yet; candidate shortlist recorded with first falsifiable
  gates and an overlap-warning against `../math/qec/`.
