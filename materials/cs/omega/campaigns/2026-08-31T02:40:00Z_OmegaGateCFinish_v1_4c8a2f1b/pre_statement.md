# PRE-STATEMENT — gate C V1 re-run under the corrected v11 interval core
# (agent OmegaPenaltyProbe, 2026-09-01T02:40Z)

Written and committed BEFORE the first compute of this campaign. Base:
continues frozen campaigns
`2026-08-30T22:58:26Z_OmegaGateCFinish_f5269412f0b1` (pre-statement
b41f384; manifest amendments dc4da1f, 5e0068a),
`2026-08-31T07:57:00Z_OmegaGateCDiag_4859327` (pre-statement 4859327),
`2026-08-31T09:18:39Z_OmegaGateCCandidateWitness_e9c6414`,
`2026-08-31T09:51:45Z_OmegaGateCEntropyRepair_1ec82a2`,
`2026-08-31T11:30:16Z_OmegaGateCEntropyRepair_v10_4cdd3c85`,
`2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5`, and
`2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959`. Live sources are
REUSED byte-identically, never rebuilt; no source file is edited in this
campaign unless the diagnostic branch below explicitly opens one.

## 0. Status and reconciliation of the inherited ticket (read first)

The assigned ticket carries an inherited state:
"V1 gate replay: Om_raw=3.272425321778391 vs frozen raw 2.3715538358350803,
residual +0.9008714859433109 ... FAIL". **This state is STALE and was
already retracted by the owner** (omega/README.md lines 590-599, agent
OmegaGateCDiag 2026-08-31, frozen in
`campaigns/2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/report.md`): the
advertised failing source (SHA-256 4504365a…) raises
`NameError: name 's_const' is not defined` at line 191 BEFORE its first
aggregate checkpoint, so `3.272425321778391` / `+0.9008714859433109` have no
reproducible producer. The ticket's diagnosis chain built on that number
(penalty-vs-eps split, 2.9e5 amplification question, stage bisection for a
0.9009 residual) inherits the stale premise. This campaign does NOT
re-adjudicate the dead branch; it runs the handoff's prescribed order on
the LIVE corrected state, whose actual open question is different and
narrower (see §1).

Facts this campaign adopts as frozen anchors (all machine-verified in their
own campaigns; re-asserted at startup here where consumed):

- A. LIVE REPLAY V1 PASS (pre-v11 core): in GateCDiag, the repaired
  gate_c_endgame.py at SHA-256 baaecb8b… (== live src copy, byte-identical
  as of write time) produced V1_frozen_machine_equality with all seven R
  branches, R=2.8170035674609757, M_low=2.0942543887102634, raw
  2.3715538358350807 == frozen 2.3715538358350807, residual vs registered
  anchor 2.3715538358350803 equal to 4.440892098500626e-16 < 5e-13: V1
  PASS. Also the narrower lemma1-off replay machine-equals the frozen
  include_lemma1=false record (raw 2.3715518061863814; shift
  2.029648699330977e-6), answering the inherited (a)/(b)/(c): (b), not
  epsilon amplification, not (c).
- B. CORE SWAP SINCE THAT PASS: the GateCDiag run executed under
  interval_core.py SHA-256 67ca5fbf… (pre-v11). The owner later found the
  outward-rounding defect in `from_two` (upper endpoint took `b.lower()`;
  README "Interval-enclosure suspension") and the live core is now v11
  SHA-256 ce70f959… with the owner-directed fixes (`from_two` uses
  `b.upper()`; endpoint-based hulls; 1/e Arb overlap; ent/nent guards;
  production v12 and the two-rung replay already ran on it). The live
  gate_c_slope_pass.py (c6b8cef2…, ==HEAD==asrun) embeds the v11-dependent
 exact-Arb leaf endpoints. **No V1 re-run of the slope replay has ever
 been executed under v11.** The frozen-control comparison in this
 campaign's decision tree provides the fresh cross-core anchor.
- C. FROZEN PROTOCOL CONTROL (Amendment 2): the frozen protocol
  reproduced standalone at Om_raw 2.3715538358351957, delta 1.155e-13 vs
  frozen raw 2.3715538358350803 — INSIDE the 5e-13 V1 gate. Fresh
  re-computation under v11 is part of this campaign's first pass (the
  two are computed in the same process, no code edit between).
- D. question status: the 21-dim kernel question is UNCHANGED and OPEN
  (GateCDiag full-box enclosure straddled zero, width 1.168202911124e-5 =
  73.86x the imported 1.5816497000997742e-7 signal, MACHINE-VERIFIED;
  CandidateWitness five-point ray all overlap zero, OPEN; v12 five rungs
  strictly worsen the repaired fixed-m objective, OPEN). Gate C stays
  PARTIAL regardless of anything below.

## 1. The ONE pass of this campaign (pre-registered decision tree)

Single script `run_v1_v11.py` (this directory), tree fixed NOW:

- Step 1 — FROZEN-PROTOCOL CONTROL on the v11 core: byte-run the frozen
  stage_b_rung2 protocol (`campaigns/2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/stage_b_rung2.py`
  with `import sys; sys.argv.append('without')` variant in-process second)
  against the live src (interval_core v11, vxxz24 31657c92…). Registers
  ctrl_with_raw, ctrl_without_raw, ctrl_M_low, ctrl_R_sum, ctrl value
  endpoint.
  * C-gate: |ctrl_with_raw − 2.3715538358350803| <= 5e-13 (the 5e-13 V1
    threshold, same instrument as V1). Expected hit: the rung-2 corrected
    replay recorded 2.3715538358350807 under v11 (two-rung campaign
    rung2_with_replay.json). If the fresh control drifts OUTSIDE 5e-13,
    STOP: report an environment/regression finding with both numbers; no
    further branch, no source edit, no diagnosis from a moving control.
- Step 2 — REPLAY V1 on the same process/core: run the live replay point
  pass (radius 0) + gate_c_endgame.aggregate exactly as the as-run
  GateCDiag script does (reuse of baaecb8b… logic; the checkpoint writer
  is the only new code and touches no computation). Registers
  replay_R_detail (7 branches), replay_R_sum, replay_M_low, replay_raw.
  * V1-gate (pre-registered, unchanged): |replay_raw − 2.3715538358350803|
    <= 5e-13 AND replay seven branches/R/M float-equal the GateCDiag
    pre-v11 records (V1_frozen_machine_equality). PASS => proceed to
    step 3. FAIL/branch-mismatch => branch B below.
- Step 3 (only on step-2 PASS) — LP/remainder stages, exactly as
  pre-registered in pre_statement_gateC2.md §2c/§4 and amended by
  pre_statement_repairC.md §3: box pass radius 1e-7, branch-safe full-box
  slope bundle, scaled midpoint LP (u-form, |V·u|_inf <= 1), ANY-y
  rigorous width, width/signal ratio, decision vocab
  {PASS-a, PASS-b, FAILURE-TO-CERTIFY} as written there. No new decision
  rule is introduced here; the existing frozen logic runs as-is.
  * Scope guardrail, carried verbatim: even total success in
    D = {delta : A·delta = 0 exactly, ||delta||_inf <= 1e-7} caps >= 12.8x
    short of the published 2.37155181; this route can NEVER produce a
    record. The LP/ANY-y outcomes certify facts about the relaxation over
    the named domains only; the 21-dimensional kernel question is
    UNCHANGED and OPEN. 'No feasible improving direction found among the
    directions tried' is NOT 'none exists' — six probes cannot certify a
    21-dimensional space. Gate C stays PARTIAL regardless.

### Pre-registered resolution tree (fixed before any run)

- C-gate PASS and V1 PASS => branch A: the defect named in the handoff is
  confirmed ABSENT from the live tree; the +0.9009 residual is finally
  instrumentally dead under the corrected core; record the fresh triple
  (control, replay, agreement) with labels; LP stage runs.
- C-gate PASS and V1 branch-mismatch/FAIL => branch B: a genuine
  v11-specific replay aggregation divergence exists. Then and only then,
  ONE diagnostic instrument (pre-registered here): re-run the replay's
  glob-level Lemma-1 aggregation comparing per-shape eps contributions
  against the frozen protocol in the SAME process, per-part first
  mismatch reported (file: gate_c_slope_pass.py glob eps block lines
  ~443-457, the named implicated component). No other edit, no second
  suspicion-driven pass; results go to the manifest as the state for the
  next campaign.
- C-gate FAIL => branch Z (expected NOT to hit): environment regression;
  stop as in Step 1.

## 2. Anchors (asserted at startup)

FROZEN_RAW = 2.3715538358350803; FROZEN_CERT = 2.3715538358544617;
defect arm = 1.9380875689560958e-11; FROZEN_R = 2.8170035674609757;
FROZEN_M = 2.0942543887102634; FROZEN R_detail 7 numbers as in
stage_b_rung2_with_lemma1.json; ctrl delta tolerance 5e-13; V1 tolerance
5e-13. All flo point floats machine-compared, never re-typed from prose.

## 3. Budget, environment, protocol

Single nice -n 10 process, OMP/BLAS threads = 1, MID = 300, __debug__ on
(asserts live). Two point passes + one box pass + LP: worst-case ~5 min
elapsed. No subdivision, no sampling, no probes outside the pre-registered
D/radius of gateC2. No formatter/linter/test-suite runs. Freeze protocol:
artifacts land here; checksums.sha256 written at freeze; single git commit
per phase.

## 4. Falsifiable outcomes

- A: replay_raw within 5e-13 of frozen under v11 core + all seven
  branches float-equal pre-v11 records + control inside gate + LP stage
  decision line (whatever it honestly states: PASS-a, PASS-b, or
  FAILURE-TO-CERTIFY with the width/signal ratio).
- B: fresh branch-mismatch table (per-block first divergence point),
  which becomes the named next campaign's starting state.
- Z: control drift outside 5e-13 => regression finding, stop.

## 5. Rule 7 sentence (mandatory; the scope of any certified sentence
## from the step-3 branch)

The swept set is and remains D = {delta in R^45 : A·delta = 0 exactly,
|delta_i|_inf <= r} around the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector (A = frozen 27x45 margin matrix, rank 24
exact over Q; kernel basis V exact rational from
kernel_basis_V_exact_rational.json, NOT a float SVD basis; r = 1e-7),
under the transcribed 3-region single-p_comp program at max_level 3,
q = 5. What was swept: one radius-0 point pass, one un-subdivided full
box with exact Arb endpoints IC.const(center_i) +/- IC.const(r), and the
LP/ANY-y instruments over D. What was NOT swept: region-1/2 glob blocks;
any part/split/lambda/region_prop or other parameter block (held at p*
as points); multi-block/cross-block direction classes; subdivision; Monte
Carlo; probes outside D; any other radius; the omega x single_mat_size
Schoenage-line 2-coordinate class; the published 2.37155181 and the
record 2.371177 (never searched); any other rung, regime, or construction
differing from the enclosed laser-method program (including 2026
asymptotic-rank/centroid-based improvements). Even total success in D
caps 12.8x short of the published 2.37155181; no record claim may appear
anywhere. Release-side commits this campaign may only touch cs/omega/.
