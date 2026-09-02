# Manifest — gate C V1 re-run under the corrected v11 interval core
# (agent OmegaPenaltyProbe, 2026-09-01; campaign 2026-08-31T02:40:00Z_…)

Headline: **BRANCH A. V1 PASSES under the v11 core — replay raw
2.3715538358350807, residual 4.440892098500626e-16 against the frozen raw
anchor 2.3715538358350803 (tolerance 5e-13) — and all seven frozen R
branches, R, and M are float-exact against both the frozen rung-2 record
and the pre-v11 GateCDiag replay records. The LP/remainder stage ran per
the pre-statement §1 step 3; rigorous verdict: FAILURE TO CERTIFY (ANY-y
enclosure straddles zero over D, width/signal 73.85977508485684). The
21-dimensional kernel question is UNCHANGED and OPEN; gate C stays
PARTIAL.**

## 0. Status and why this campaign exists (DIAGNOSTIC-ONLY completion of
## the inherited handoff)

The assigned ticket said: "V1 gate replay: Om_raw=3.272425321778391 …
residual +0.9008714859433109 → FAIL… the defect is decisively in the
slope replay path… zero the lam_sum penalty … band the residual." That
premise was already DEAD before this campaign started:

- omega/README.md lines 590-599 (owner retraction, 2026-08-31): the
  advertised failing source (SHA 4504365a…) raises
  `NameError: name 's_const' is not defined` at line 191 BEFORE its first
  aggregate checkpoint. `3.272425321778391` and `+0.9008714859433109`
  have no reproducible producer; they remain REPORTED history only.
- GateCDiag (2026-08-31, frozen
  `campaigns/2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/`) had already
  repaired the replay, PASSED V1 (delta 4.440892098500626e-16), answered
  the inherited penalty diagnostic with the exact narrower lemma1-off
  control (shift 2.029648699330977e-6, MACHINE-VERIFIED against the
  frozen include_lemma1=false record; no ~3e5 amplification; the
  (a)/(b)/(c) answer is (b)), and run the full-box LP/ANY-y instrument
  (74x signal, straddles zero, OPEN). The pre-statement of THIS
  campaign re-derives all of this and pre-registers the fresh question
  that remained genuinely open: **no slope-replay V1 had ever been
  executed under the corrected v11 interval core** (the owner's
  interval-enclosure suspension landed after GateCDiag; v12 and the
  two-rung replay ran under v11 but the slope replay did not).

## 1. What ran (one pass; single nice -n 10 process; threads=1; MID=300)

`run_v1_v11.py` (this directory), implementing the pre-registered tree:

1. **FROZEN-PROTOCOL CONTROL on the v11 core** (C-gate): in-process
   byte-run of the frozen
   `campaigns/2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/stage_b_rung2.py`
   with and without the Lemma-1 charge against live deps.
   - ctrl with-lemma raw **2.3715538358350807** → delta vs anchor
     **+4.440892098500626e-16** ≤ 5e-13: **C-GATE PASS**
     (COMPUTATIONAL-EVIDENCE label on the fresh run; the identical value
     is already MACHINE-VERIFIED as the frozen raw in
     stage_b_rung2_with_lemma1.json).
   - ctrl without-lemma raw 2.3715518061863814, R_sum
     2.817007818061672 — both float-identical to the frozen
     without-lemma record.
2. **REPLAY V1 on the same process/core** (MACHINE-VERIFIED):
   - replay raw **2.3715538358350807**, delta **+4.440892098500626e-16**
     ≤ 5e-13: **V1 PASS**.
   - All seven R branches float-equal the pre-v11 GateCDiag records AND
     the frozen with-lemma record:
     R_comp[0,3]=0.24530661807075219, R_comp[1,3]=0.24715974210930786,
     R_comp[2,3]=0.2456557897939231, R_comp[0,2]=0.5883309785276519,
     R_glob[0]=0.49684932826618977, R_glob[1]=0.4968521822710635,
     R_glob[2]=0.49684892842208755.
   - R_sum=2.8170035674609757 (== frozen), M_low=2.0942543887102634
     (== frozen), lemma1-off raw 2.3715518061863814 (== frozen without
     record). The replay aggregation is **byte-stable across both
     cores** for every consumed float.
3. **LP/remainder stage per gateC2 §2c/§4 as amended** (only on the V1
   pass, per the handoff; all ran):
   - Box pass r=1e-7 on the region-0 glob dist[0] block; all 45 leaf
     endpoints asserted exact (`IC.const(center) ± IC.const(radius)`)
     under the construction prec envelope (MACHINE-VERIFIED). The
     binary64 counterfactual repeats the GateCDiag plant exactly: 24
     lower and 25 upper endpoints would round INWARD under float
     subtract/add (identical coordinate sets) — the old path could not
     enclose the registered box.
   - Full-box slope bundle: width_max 4.640574461822503, L1 width sum
     113.66562530936552 (pre-v11: 4.640574461822501 / 113.66562530936517
     — float-display-level agreement; both enclosures straddle with the
     same sign structure).
   - Scaled midpoint LP (u-form, |V·u|≤1, exact rational JSON basis with
     the integer A·V=0 asserts re-run in-script): min
     -1.5773329015168567e-7 (COMPUTATIONAL-EVIDENCE only; pre-v11
     -1.577332901516852e-7).
   - **ANY-y rigorous instrument** (MACHINE-VERIFIED, Arb end-to-end):
     fitted arbitrary y ⇒ residual L1 upper 58.41014555619974; swing
     upper 5.841014555620028e-6; enclosure
     [-5.841014555620028e-6, +5.841014555620028e-6]; width
     1.1682029111240058e-5; zero-y counterfactual 1.37871966251969e-5
     is strictly wider (the instrument demonstrably responds).
   - Width / signal 1.5816497000997742e-7 = **73.85977508485684**
     (pre-v11 73.85977508485647; identical to display precision; the
     ratio label stays INFERENCE, the enclosure MACHINE-VERIFIED).

**Decision (frozen vocab, gateC2 §4): FAILURE TO CERTIFY — the ANY-y
slope enclosure straddles zero over D; the 21-dimensional question
remains OPEN.** Not PASS-a, not PASS-b, not "no improving direction":
six LP probes cannot certify a 21-dimensional space.

## 2. Regression anchor (pre-statement §0.C, all in one process, no edits)

The frozen-protocol control and the failing-till-now replay AGREE on the
same machine, same data, same anchor, same core, same process — inside
5e-13 on every consumed float (C-gate 4.44e-16 and V1 4.44e-16; seven
branches exact). This is the fresh cross-core evidence the pre-statement
asked for: the defect-class named in the handoff (a replay aggregation
defect producing a 0.9-scale residual) is **instrumentally dead** on the
live tree; the +0.9009 state it references was retired by the owner
retraction and is now double-confirmed non-reproducible.

## 3. Where the remaining gap lives (unchanged, restated)

The full-box slope bundle width (L1 113.67 → enclosure width 1.168e-5)
is 73.9x the 1.58e-7 signal; the branch-straddle uncertainty at the
base point (R interval span ~4.25e-6 across the seven components; Lemma-1
residual intervals straddling zero) is the binding constraint, exactly
as frozen in GateCDiag/CandidateWitness. Next named action (unchanged
from the frozen README): subdivision (branch-and-bound) of the box so
each subbox's branch structure is pinned, or a tighter certified
gradient enclosure; both must be a new pre-registered campaign. Even
total success in D caps ≥12.8x short of the published 2.37155181 —
this route can NEVER produce a record.

## 4. Contents

- pre_statement.md — pre-registration, committed ec16efe BEFORE first compute
- run_v1_v11.py — the one-pass script as run. Three runner-only defects
  were hit and fixed BEFORE the first completed pass (disclosed, not
  backdated): an `exec` namespace dict bug, a box-endpoint assert placed
  outside the construction prec envelope, and a missing timing variable.
  None touch aggregation or control math; aborted attempts produced no
  verdict (stderr tracebacks retained in the logs); every re-run
  restarted at step 1, so the completed pass re-derived control+V1
  together in one process.
- run_v1_v11_checkpoint.json — the machine checkpoint (labels inline)
- run_v1_v11_stdout.log / run_v1_v11_stderr.log — raw run transcripts
  (the completed pass has an EMPTY stderr; the aborted attempts'
  tracebacks were overwritten by each re-run's redirect and survive only
  in this manifest's description — disclosed as a freeze-hygiene gap)
- checksums.sha256 — sha256 of every file in this directory (freeze)
- (the decision block inside the checkpoint is the authoritative verdict)

## 5. Rule 7 sentence (mandatory)

This campaign swept exactly: (i) one radius-0 point pass and one
un-subdivided full box over D = {delta in R^45 : A·delta = 0 exactly,
|delta_i|_inf ≤ 1e-7} around the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector (A = frozen 27x45 margin matrix, rank 24
exact over Q; kernel basis V exact rational, integer-asserted in-run;
NOT a float SVD basis), with the scaled midpoint LP and the ANY-y
instrument evaluated over that un-subdivided D; and (ii) the same
frozen stage_b_rung2 protocol re-run at the anchor point as control —
all under the transcribed 3-region single-p_comp program at max_level 3,
q = 5. NOT swept: region-1/2 glob blocks; any part/split/lambda/
region_prop or other parameter block (held at p* as points); multi-block
or cross-block direction classes; any subdivision of the box; Monte
Carlo; any other radius; probes outside D; the omega × single_mat_size
Schoenage 2-coordinate class; the published 2.37155181 and the record
2.371177 (never searched; nothing here bears on them); any other rung,
max_level/q regime, or construction differing from the enclosed
laser-method program (including the 2026 asymptotic-rank/centroid-based
improvements). Even total success in D caps 12.8x short of the published
2.37155181; no record claim may appear anywhere.
