# Manifest — gate C campaign 2026-08-30T11:24:41Z (OmegaGateC)

Gate C, VXXZ24 rung centered (per Main steering; record 2.371177 has no
published parameters as of 2026-08-30 — artifact audit re-run by Main —
so no other centering is possible).

## Contents
- `pre_statement.md` — PRE-REGISTERED before first B&B call: box (r=1e-7,
  all 6759 coordinates), branching/screen rule, 63-box budget, 240-min
  wall cap, stage ordering, escalation triggers, rule-7 sentence.
- `precision_sweep_with_lemma1.json` / `precision_sweep_without_lemma1.json`
  — Part 1: 5-precision decomposition of the 2.0258545e-6 gap. BOTH arms.
- `part1_gap_decomposition.md` — Part 1 analysis + REVISED structural
  attribution (§Final supersedes intermediate reading).
- `stage_a_screen.json` + `stage_a_screen_record.md` — stage (a) float
  screen (mandated first; anchor clean; probes recorded including
  infeasible-at-exact-equality ones; NOT certified).
- `stage_b_boxes.json` + `stage_b_boxes_record.md` — stage (b) interval
  box screen: B0 anchor + 14 signed structured boxes at r=1e-7,
  per-box certified endpoints with interval widths as first-class
  numbers; zero-crossing guard clean everywhere.
- `vxxz24_float.py`, `interval_core.py`, `stage_b_rung2.py` — byte-copies
  of the src/ machinery AS RUN (includes this campaign's edits:
  G.__init__ object-array branch, ParamManager.set_value box coupling,
  IntervalTree.get box wiring with exact-dyadic endpoints, ent_vec
  pre-registered zero-crossing guard with ZC_FLAG).
- `gate_c_precision_sweep.py`, `gate_c_stage_a_screen.py`,
  `gate_c_stage_b_boxes.py` — this campaign's scripts as run.
- `checksums.sha256`, `tool_versions.txt` — provenance.

## Provenance / evidence labels used
- Anchors (B0 = rung-2 frozen endpoint): MACHINE-VERIFIED.
- Part 1 precision sweep endpoints and widths: MACHINE-VERIFIED.
- Stage (a) float screen: COMPUTATIONAL-EVIDENCE (float64, no enclosure).
- Stage (b) per-box endpoints: MACHINE-VERIFIED (outward-rounded interval
  aggregation; per-box params below).
- Published VXXZ24 2.37155181: CITED-DEPENDENCY (osf.io/7wgh2, file read
  first-hand by the 05:20:00Z campaign; sha of .mat recorded there).

## Inputs
- `scratch/rmmcode/data/K100_2.37155181.mat`
  sha256 df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0722c3e25871da93d
  (6759 float64 params, 3-region single-p_comp program, max_level 3, q 5).

## Determinism / reproduction
Single process, OMP/OPENBLAS threads = 1, `nice -n 10`, no RNG. Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 10 \
      ./.venv/bin/python gate_c_stage_b_boxes.py <outdir>/stage_b_boxes.json 1800
from the frozen directory (scripts are self-locating via absolute sys.path
into cs/omega/src; the byte-copies here are the versions that ran).
