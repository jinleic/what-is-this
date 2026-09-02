# Pre-statement — Campaign C: minimal attribution-clean cap restoration,
# tile 3 at panel cap 2^19, then band resume tiles 4–49

Agent `KgCampC`, UTC 2026-09-01. **COMPUTE PENDING — this file is committed
before any numerical evaluation, probe, or band run of this campaign.** The
only prior files of this directory are the inherited instrument copies
 in commit `821edf5`, whose checksums are recorded in `checksums_addendum.sha256`
at freeze time (and match Campaign B's ledger byte-for-byte for the two
inherited files).

Predecessors: Campaign A
`campaigns/20260831T231500Z_KgEvenBoxRerun_campA_corr_env` (frozen commit
`8d2cb26`) and Campaign B
`campaigns/20260901T051620Z_KgCampB_subdiv_env` (frozen commit `f744dd2`).
Campaign B's quantified obstruction
(`obstruction_tile3.md`/`obstruction_tile3.json`, owner-arithmetic-checked) is
accepted as-is: tile 3 of band [4.083, 6.0] stays OPEN at A's panel-cap
131072; the pre-registered mq = t² subdivision recovers exactly 0 on the
frontier cell (crossing root s* = 45.703407790989350730 outside [0,8]); and
the BINDING constraint is the panel cap — A's envelope at cap 2^19 already
certifies the frontier cell's leaf at margin
+4.46309091764170239642742026145e-6 (single-cell certified run, B's
decomposition D3). The kgcloseB3 retraction STANDS; no new PASS has been
earned; nothing here bears on the paper's correctness in either direction.

## 0. Provenance of the single knob (not domain-shopping — a disclosed weakening restored)

Campaign A's ORIGINAL pre-registration (its Addendum 1) fixed the panel cap
at **524288 = 2^19**. A's Addendum 2 later weakened it to 131072 = 2^17 as a
pure cost measure ("Recomputed cap 131072 … the projected band drops to a few
hours"), with the disclosed tightness loss on the kernel slab. Campaign B's
decomposition then measured that on the actual failing ring cell the weakened
cap loses strictly (deficit −1.7629e-5 at 2^17) while the ORIGINAL cap 2^19
wins strictly (+4.4631e-6, certified). Restoring A's original pre-registered
value 2^19 is therefore reverting a DISCLOSED WEAKENING, not shopping the
domain for a pass: the c-tiling, root cover, box floor, depth, panel ladder
rungs (4096→16384→65536→262144→524288), S_MAX, tail, precision, verdict
semantics, and DFS order are all byte-shared with A; exactly ONE knob changes
so the outcome is attributable to it.

## 1. What is re-used, byte-identical (declared now)

From Campaign A's frozen `code/` (sha256 at that campaign's `checksums.sha256`
and in Campaign B's `checksums_addendum.sha256`, byte-verified on copy):
`gate_corr_env.py` (sha 94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30)
— envelope cores `e_panel_range`, `dist0_box`, `box_intersects_disk`,
`root_boxes`, `split_box`, `can_split`, tiling `c_tiles` ratio 1.008, ladder
4096×4, floor `BOX_FLOOR_Q = 1/1024`, `BOX_DEPTH_MAX = 40`, budget
semantics, `strict_margin_lower`, S_MAX = 8, precision 256 — and
`startup_controls.py`
(sha cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3), the
full battery: kernel-admission incl. cap-at-1 counterfactual, 2-D vs 3-D
odd-radius counterfactual, K_o/K_e identities, cited-tail anchor with the
byte-verified d3h sha 80e945589b…, margin-predicate trio, envelope-admissibility
on p*. **No byte of either is edited.** The ONLY delta of this campaign's
instrument (`gate_c2p19.py`, written fresh in `code/` but importing A's
modules as-is) is the module-level knob value:

- A: `PANEL_CAP = 131072` (its frozen file fixes this constant);
- C: effective cap **524288**, applied by patching `gce.PANEL_CAP = 524288`
  right after import and BEFORE any envelope evaluation. `panel_ladder`
  reads the module global, so the ladder becomes
  4096 → 16384 → 65536 → 262144 → 524288 (×4, unique, capped once) — the
  same ladder function, longer by one rung; no other code path differs.
  A frozen equi-run check in the Phase-1 driver asserts
  `panel_ladder(c) == [4096, 16384, 65536, 262144, 524288]` on this band
  before any cell evaluation.

Shared dependency `20260831T082425Z_kg_direct_d4_restart/code/core.py` is
imported, not copied (as in A and B). python-flint landmine discipline
(owner rule): every interval used as an interval is built in midpoint/radius
or explicit-endpoint/outward form; `arb(lo, hi)` midpoint/radius constructor
is never used as an interval constructor. All arithmetic inside
`Prec(256)`; 320 bits only inside the inherited startup battery, exactly as
frozen.

## 2. PHASE 1 — tile 3 re-run ENTIRELY at cap 2^19 (single knob, clean attribution)

Tile: the pair (c_L, c_R) = (4.1485893120, 4.181778026496) = exactly
`gce.c_tiles('4.083','6.0')[2]`, asserted adjacent to the A/B pair. Root
cover 132 (asserted), DFS order/split rule/bookkeeping byte-shared with A's
`certify_tile` (C drives it via the imported module with only the cap knob
patched — see §1). Budget: 260,000 distinct (box, panel-count) envelope
evaluations for Phase 1 (A's per-band budget, unchanged in meaning).
Stop rule: identical to A's (an open leaf at panel cap and floor depth ⇒
tile OPEN ⇒ Phase-1 FAIL). Optional but default: the run uses A's DFS from
the 132 roots (not a resume from A's stack), so the tile verdict is earned
end-to-end under the C cap.

## 3. Startup gate (before ANY tile compute) — full battery re-run

The inherited `startup_controls.py` runs **16 check families** (16/16 green in A and B); this campaign
re-runs the full battery in-place and instances it in `logs/startup_controls.json`:

1. `kernel_unit_norm`; 2. `kernel_e0_is_sqrt32` (e(0) = √(3/2) admitted);
3. `kernel_admitted_by_envelope`; 4. `kernel_envelope_tight` (slack 1.87e-96);
5. `cap1_counterfactual_rejected` (deficit 0.224744871391589… certified);
6. `dist0_exact_value`; 7. `three_d_inflation_is_looser` (3-D-inflated 1.6542
   vs corner-exact 1.3243);
8. `two_d_envelope_corner_exact_at_s0`; 9. `d3h_cert_sha256` (80e945589b…);
10. `Ko_Ke_identities`; 11. `B3_import_is_upper`; 12.
    `margin_straddle_rejected`; 13. `margin_tight_accepted`;
14. `margin_above_rejected`; 15. `pstar_closed_form_window`;
16. `envelope_admits_pstar`.

ANY failure ⇒ STOP, freeze nothing else, report (a failed control means the
inherited instrument is broken — attribution to the knob would be void).

## 4. Reconciliation gate (before any cap-2^19 number is reported)

Before Phase 1's verdict, C recomputes A's frozen tile-3 open margin under
the C instrument **with the 2^17 cap** on A's frozen frontier cell
a0 ∈ [−638/768, −637/768], a2 ∈ [429/768, 430/768]. Gate: the margin must
reproduce −1.76291787028287693244667476546e-5 as an Arb ball of radius
< 2^-40 around A's frozen value (B achieved bit-exact reproduction of A's
frozen display string; the reproduced margin must satisfy
`(margin − A_frozen).rad() < 2^-40`, which for an exact-ball difference of
zero means bit-identical). **If the gate fails, STOP: the campaigns disagree
and knob attribution is void** — freeze the disagreement exactly (inputs,
both margins, the failing comparison) as the campaign's obstruction, run
nothing further, and report.

Only after the gate passes: the same cell at cap 2^19 is computed (this
number is fresh compute, not a reconciliation) and reported as the first
Phase-1 leaf datum, alongside B's D3 certified ball
+4.46309091764170239642742026145e-6 ± 3.26e-36 for cross-campaign
consistency (comparison labelled COMPUTATIONAL-EVIDENCE cross-check; the
C number itself is MACHINE-VERIFIED from this campaign's run).

## 5. PHASE-1 VERDICT — threshold fixed now

Tile 3 PASSES iff **every disk-intersecting leaf of tile 3, including the
frozen frontier cell and A's certificate-corner leaf, certifies with a
margin ball whose LOWER end is strictly positive** (a straddling ball is
NOT a pass; the margin predicate is the inherited
`d(c_U).lower() − U(c_L).upper() > 0` and its ball lower end is what must
be > 0). Tile 3 FAILS iff any leaf is OPEN at cap 524288 (reason "panel cap
reached", or budget/stack exhaustion per the byte-shared rules).

Reported on PASS and FAIL alike: the per-leaf worst margin (Arb ball, full
rule-15 subinterval form), leaf count certified / total, achieved max depth
reached during evaluation (with floor leaves distinguished from non-floor
open leaves), max panels used across leaves, per-rung panel histogram, and
`time.process_time` CPU (rule 17e: in-process only).

## 6. PHASE 2 (only on Phase-1 PASS) — resume tiles 4–49 at cap 2^19

- **Tiles:** 4–49 of the same ratio-1.008 tiling, adjacent continuation from
  tile 3's right edge c_R = 4.181778026496 through 6.0, final tile clamped
  to exactly 6.0 (the same 49-tile tiling; adjacency asserted per tile).
  Tiles 1–2 PASSes stand as A's (not re-litigated — outside this campaign's
  budget).
- **Stop-at-first-failing-tile rule:** inherited verbatim. A tile that stops
  OPEN freezes the exact frontier (open cell, margin ball, depth, panels,
  remaining stack) and the band resumes NO further — no tile beyond the
  failing one is evaluated.
- **Global CPU budget (pre-registered): 40 CPU-hours** (time.process_time,
  sum over tiles), measured against a wall-clock ceiling of 72 h. Rationale
  (COMPUTATIONAL-EVIDENCE): A's tile-1/2 runs cost ~2,665/3,455 s wall at
  cap 2^17; the cap-2^19 ladder adds one rung ≈ 4× the 2^17 rung on
  non-skipped panels; mid-band tiles are cheaper than kernel-adjacent ones
  (A's tile-3 cost 863 s). 40 CPU-h is a hard budget: on exhaustion the
  band freezes at its last completed tile with the exact frontier.
- **Stop rule (pre-registered):** stop when (a) the first failing tile
  opens, or (b) tile 49 certifies (band complete), or (c) the global CPU
  budget binds mid-tile. In case (c) the interrupted tile is reported as
  OPEN with reason 'budget exhausted' — it consumes no verdict — and the
  campaign freezes.
- **Incremental freeze:** after each tile completes (pass or fail), its
  serialised row is written to `logs/inproc_tile_NNN.json` (rule-15 fields)
  and a running `logs/band_resume_progress.json` accumulates all completed
  rows + panel histogram + CPU seconds. An interrupted run loses at most
  the tile in flight; completed tiles are never recomputed.
- **Resumability at tile granularity:** `--start-tile N` (default 4 per the
  phase plan; any integer 1..49 accepted) evaluates tiles N..49 in order;
  a stop file records the next unstarted tile index. A resume skips
  completed tiles by reading any existing `inproc_tile_NNN.json` rows
  (their checksums in the addendum ledger identify them) and continues
  from the stored next-tile pointer. Per-tile c-pairs are regenerated from
  the byte-shared `c_tiles` and asserted equal to A-style adjacency —
  the tiling is deterministic, so a resumed run sees the same tiles.
- Per-tile artifacts follow A's format exactly (roots 132 asserted, evals,
  panel histogram, margin balls, open fields where applicable).

A completed band = all 49 tiles certified = band PASS (first under a sound
envelope since the retraction; the kgcloseB3 retraction stands regardless).
Any failing tile ⇒ band FAILURE TO CERTIFY at that tile with the exact
frontier, exactly as A filed it.

## 7. PHASE 3 — only if Phase 1 FAILS: deriving B's D4 soundly (or freezing the blocker)

Campaign B's D4 (disk-feasible |e|-sup) was correctly NOT adopted: it is a
diagnostic because coefficient-space feasibility of the even part B∩disk was
applied WITHOUT re-deriving the odd radius jointly — the paper's
`|o(s)| ≤ r_o√K_o` route requires, for a box straddling the disk boundary,
that the full box's |e|-sup be covered, and a disk-restricted even sup is
only sound jointly with an odd-slope bound re-derived under the same
restriction. If Phase 1 FAILS, Phase 3 attempts:

1. Quote the paper's 2-D even-box form verbatim (lines 1972–1981, frozen in
   Campaign A's pre-statement §1 and re-checked first-hand against
   `scratch/paper_full.txt` this session — the same text both A and B used).
2. Derive the joint bound: for a unit cubic with even part in
   B∩disk (sup re-stricted accordingly) AND |‖o‖₂| ≤ √(1−|‖e‖₂|²) under the
   same restriction, produce an outward-certified per-panel
   (E⁺\*, W⁺\*) pair. Clean-room derivation; no universal quantifier is
   sampled (repo rule 7).
3. Counterfactual battery on the joint form (pre-registered): it must ADMIT
   the paper's kernel point (e(0) = √(3/2) at the kernel box), and it must
   be shown LOOSER-than-or-equal-to every variant it replaces in the
   permissive direction (the joint bound on B∩disk must be ≤ the
   corresponding full-box bound for boxes entirely inside the disk, and
   their difference certified ≥ 0 on demo boxes).
4. Only if (1)–(3) all pass is any D4-style mechanism allowed into a
   verdict-bearing re-run of tile 3 — in a NEW pre-statement addendum, not
   retro-fitted here.
5. If the derivation does not close: freeze the blocker exactly (which
   step failed, what would additionally be needed), as the campaign's
   obstruction file, with the Phase-1 FAIL accounted as-is. NOT adopted
   then: B's D4 diagnostic stays a diagnostic; tile 3 remains OPEN; no
   Phase-2 resume runs.

## 8. Refutation gate (inherited verbatim)

A certified statement J.lower() > d.upper() on any leaf of any tile = a
candidate refutation of the paper's Lemma D.2 — STOP the whole campaign,
escalate to Main via hub with the exact cell and intervals, and write
nothing further. (Predecessor campaigns never fired this; neither may we
casually.)

## 9. Labels

- MACHINE-VERIFIED: any margin/sign/interval claim from outward Arb
  endpoints at the stated domains (256-bit outward), scripts frozen in
  `code/`.
- COMPUTATIONAL-EVIDENCE: timings (time.process_time /
  resource.getrusage only, rule 17e), cross-campaign comparisons, and any
  projected costs.
- REPORTED: paper Table-1 numbers, Campaign A tile-1/2 PASSes, Campaign B
  diagnostics.
- FAILURE TO CERTIFY is a result and is reported honestly; PASS requires
  strictly positive margin-ball lower ends on every leaf of every in-scope
  tile.

## 10. Scope (rule 7, fixed now)

Swept here (conditional on phases passing): Phase 1 = exactly tile 3 of
[4.083, 6.0] under A's byte-shared envelope with the single knob PANEL_CAP
restored to 2^19 (132-root cover, full even disk, all compatible odd parts,
[0,8] panels + closed tail, 256-bit outward Arb); Phase 2 = tiles 4–49 of
the same band under the identical knob state, with pre-registered global
CPU budget and stop rules. NOT swept: tiles 1–2 re-runs (A's PASSes stand;
B reconciled them) except as budget-leftover orientation; any cap value
other than 2^19 on verdict-bearing tiles (the 2^17 reconciliation cell run
is a gate, not a verdict); the mq = t² subdivision mechanism (B's Phase-1
instrument — NOT adopted here; its correct domain was shown empty for this
band's ring cells); any coefficient-space envelope change (Phase 3 is
derivation-gated, verdict-adopted only via a new addendum); bands outside
[4.083, 6.0]; the never-run bands ([1.0,1.3], [1.45,1.75], [1.75,3.5],
[6.0,12.0]); the paper's C1/splice/C3 certificates; degree > 3 families;
non-unit cubics; and Constructions outside the unit-cubic threshold family.
The kernel-admission/cap/3-D counterfactual battery sweeps only its 16
planted asserts.

## 11. Rule-15 discipline

Every certified subinterval (leaf cell, margin, tail, c-pair) is reported
with exact fmpq endpoints or full Arb ball display (magnitude + radius), so
the owner can re-derive each number from the frozen artifacts. The frontier
cell coordinates are exact fmpq (−638/768 etc. — inherited from A's frozen
strings at 30 display digits, re-derived to exact fractions this session
from the tiling arithmetic: 1/12-root grid ÷ 2^depth ÷ 1/1024 floor).

## 12. Campaign C result (frozen below, appended after runs)


## CAMPAIGN C RESULT (frozen; checksums in `checksums_addendum.sha256`; appended after runs)

**Phase 1: PASS — tile 3 of band [4.083, 6.0] certifies COMPLETELY at the single
restored knob PANEL_CAP = 2^19. Phase 2 was DEFERRED BY OWNER DECISION ON COST
GROUNDS (quantified below); it is NOT an instrument failure. No band verdict is
claimed for [4.083, 6.0] beyond tile 3.**

Gates (all before any 2^19 verdict number):

| gate | result | evidence |
|---|---|---|
| startup battery 16/16 | PASS, MACHINE-VERIFIED | `logs/startup_controls.json` (sha 211d8e9a…, byte-identical to A/B ledgers) |
| reconciliation (cap 2^17 vs A's frozen margin) | PASS — bit-identical, radius 0, diff −4.19e-35 < 2^-40 | `logs/phase1_gates.json` field `reconciliation_2p17` |
| frontier cell at cap 2^19 | PASS margin [4.46309091764170239642742026145e-6 ± 3.26e-36] | same file, `frontier_cell_2p19` |

Phase-1 tile verdict (MACHINE-VERIFIED, Arb 256-bit outward; artifact
`logs/inproc_tile_003.json`):

* c-pair (4.1485893120, 4.181778026496) — exactly A/B's tile 3 (adjacency
  asserted from the byte-shared `c_tiles('4.083','6.0')`).
* **484/484 disk-intersecting leaves certified, ZERO open**; DFS completed with
  the stack empty (1120 boxes evaluated, 4140 envelope evaluations).
* **Worst certified margin +5.06389732008514710540159866096e-7 ± 2.76e-37**
  (strictly positive lower end) at leaf
  a0 ∈ [−0.8203125, −0.817708333333333333333333],
  a2 ∈ [+0.5703125, +0.572916666666666666666667], panels 16384.
* Panel histogram (by ladder rung): 4096: 1120, 16384: 908, 65536: 784, 262144: 682, 524288: 646
  ; CPU 12,330 s (`time.process_time`, rule 17e).
* Single-knob attribution: every other instrument byte is A's
  (`gate_corr_env.py` sha 94c07adf…; `startup_controls.py` sha cf7d57fb…);
  the ladder became [4096, 16384, 65536, 262144, 524288] — the pre-registered
  assertion checked it before any compute.
* A's frozen frontier leaf (a0 ∈ [−638/768, −637/768], a2 ∈ [429/768, 430/768])
  certifies at cap 2^19 with margin +4.46309091764170239642742026145e-6,
  agreeing with Campaign B's D3 certified ball to 3.26e-36 (within B's own
  radius): the −1.7629e-5 deficit was purely panel-cap slack, as B's
  decomposition concluded.

Refutation gate: never fired — all 484 leaf margins are strictly positive; no
certified J > d interval was produced anywhere.

**Phase-2 deferral (owner decision, pre-compute — `logs/phase2_deferral.json`):**
tiles 4-49 at cap 2^19 are INFEASIBLE within the pre-registered 40 CPU-h
budget: optimistic 78.8 CPU-h (measured tile-3 cost halved per tile), like-tile-3
157.6 CPU-h, measured-multiplier 265 CPU-h — min 1.97× over budget. Zero Phase-2
CPU was spent. The instrument is NOT invalidated: this is a budget-scoping
decision on a measured cost model. A's tile-1/2 PASSes stand as A's; A's tile-3
OPEN stands as a statement about A's own cap-2^17 instrument (not retracted);
the band [4.083, 6.0] has NO overall PASS claim from this campaign — tile 3
certified under the 2^19 knob is the entire executed scope.

Refutation gate never fired; all leaf margins strictly positive; no J > d
interval anywhere.

Rule-7 scope: as pre-statement §10, executed subset — the tile-3 verdict under
the 2^19 knob; the cap-2^17 reconciliation cell; the cap-2^19 frontier datum.
NOT swept: tiles 1-2 re-runs, tiles 4-49, other caps as verdict knobs, other
bands, the mq = t² mechanism (B's, not adopted), joint coefficient-space
envelopes (Phase 3 never triggered — Phase 1 passed), degree > 3 families,
non-unit cubics, constructions outside the unit-cubic threshold family.

Named next action: Campaign D — fresh pre-statement for tiles 4-49 at cap 2^19
with budget ≥ 160 CPU-h and tile-granular checkpointing (mechanism proven in
`logs/c_tile3_state.json`), or a jointly-derived (even, odd) coefficient-space
tightening first. Queued behind the three never-run band campaigns with Main;
`cs/kg` is released.

Disclosed defects (rule 5): two detached `nohup` launches of the Phase-1
driver were throttled ~40× by the macOS QoS/Python.app mechanism Campaign A
documented; both were killed, produced zero verdicts, and their `sample`
captures are frozen (`logs/sample_*.txt`) as the monitoring record. The first
chunk engine's margin min/max update carried a comparison bug (an `if/else`
truthiness slip on an `arb`), repaired before chunk 2; the affected state file
was discarded and recomputed from the corrected path, so every frozen number is
machine-recomputed under the corrected code. The empty `phase1_stdout.log` is
the killed second nohup launch's output; the verdict-bearing run executed in
the session kernel with the same frozen code bytes, checkpointing per chunk to
`logs/c_tile3_state.json` (final: stack empty, 1120 boxes, 484 certified).
