# Manifest — gate C stage-3 campaign (agent OmegaSlope, 2026-08-30)

Campaign: certified slope (forward-mode interval AD) through the transcribed
VXXZ24 tree + certified second-order remainder, deciding the 21-dim
nullspace question. Pre-statement committed BEFORE first compute:
`cs/omega/pre_statement_gateC2.md` (git 83e3f09, 2026-08-30T17:5xZ).

## Outcome: FAILURE TO CERTIFY (machinery incomplete, reported per rule 5 / pre-statement §4 FAIL mode)

The certified slope pass (gate_c_slope_pass.py + gate_c_slope_core.py)
BUILDS and RUNS end to end: 45-coordinate interval-AD through all 135
level-3 parts, 1377 level-2 parts, glob stage, p_comp chains and the
Lemma-1 penalty terms — two passes (point r=0 + box r=1e-7) in ~110 s
each, all arithmetic in Arb endpoint intervals inside IC.prec().

The VALUE re-aggregation plugged on top of it is INCOMPLETE and produces
WRONG numbers (V1 failed: point re-aggregation 41.9962 vs frozen certified
2.3715538358350803, delta ~+39.62). Root cause, identified and isolated
after the failure: the replay never synthesizes the level-2 PartLv2
`mat_size_contribution` (~2.094 of the ~2.094-point mat_size), so M_low
collapses from 2.09425 to 0.14985 and Omega_raw inflates by the missing
2.094/M ≈ 14x plus R-side gaps. The slope BUNDLES it produced are
consistent with this wrong value definition (gradient scale ~100s vs the
true ~2.5), i.e. the failure is in the aggregation layer, not the AD core;
but under the repo's discipline the two cannot be cleanly separated at
this budget, so NOTHING from those runs is a certified or even
floating-point claim about Omega. Zero numbers from past endgame_test.json
enter any README or index claim.

Per pre-statement §4 FAIL mode: reported with the exact step where the
width/VALue blows up and with no domain move, no re-centering, no budget
extension. Open question status: UNCHANGED AND OPEN (the predecessor's
scoped negative stands; the 21-dim nullspace question remains open).

## Contents
- pre_statement_gateC2.md copy — the pre-registration (original: cs/omega/pre_statement_gateC2.md)
- gate_c_slope_core.py — interval-AD core: Slope(val Ival, grad array-of-Ival),
  Leibniz/quotient/entropy/nent chain rules in Arb.
- gate_c_slope_pass.py — the certified replay of the transcribed tree
  (Workspace/GlobalStage/Part/PartZero/PartLv2 evaluate hooks in order).
- gate_c_endgame.py — endgame driver: value re-aggregation with branch
  pinning, V1–V4 validations (percent-pass/float-grad/kernel-grad/ZC),
  box pass, LP over D (scaled + unscaled), decision logic.
- endgame_partial_checkpoint.json — the raw checkpoint state of the final
  (failing) run; kept for provenance, NOT evidence.
- raw_stdout_last_run.log — stdout of the final run.
- interval_core.py, vxxz24_float.py — byte-copies of the frozen interval
  machinery AS USED.
- checksums.sha256 — sha256 of every script.

## What the next agent must fix (the named, exact gap)
1. PartLv2 post: build t.mat_size_contribution (Slope form of
   PartLv2.evaluate_post's `inner = entropy_vec([s0, s0, 1-2 s0])
   + 2 ln q (1 − 2 s0)`; contributions rotated/frac-multiplied exactly as
   the transcription does) — currently the replay sets only
   num_block_contribution; ALSO the msc loop must skip PartZero's zeros.
2. Re-aggregate and re-run V1; ONLY on a V1 delta < 5e-13 proceed to the
   LP/decision stages that are already coded.
