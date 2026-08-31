# Rung 2 (VXXZ24 omega<=2.371552) — structural delta assessment (SCOPING)

Written 2026-08-30T04:10Z after scoping diff review; per Main timebox this
is the scoping deliverable, not the run itself. Run deferred if the diff
exceeds a small delta.

## Diffs assessed (VXXZ24 = scratch/rmmcode/, Alman25 = scratch/alman_code/)

- evaluation/GlobalStage.m: ~200 line diff. Key changes:
  3 regions -> 6 regions; `p_comp` (one function, Z-compat) ->
  `p_compY`/`p_compZ` (Y and Z compats, both used); `part_frac` -> `term_frac`;
  `parts` -> `terms`; extra symmetry constraint in omega-mode.
- evaluation/Workspace.m: ~196 line diff. Introduces explicit
  `num_retain_glob{r}` groups (1x6 instead of implicit); adds explicit
  Schonhage line constraint `value >= power*log(q+2)`.
- evaluation/PartInfo.m ↔ TermInfo.m: ~374 line diff (renaming + expanded
  asymmetric hashing in the 6-region version).
- evaluation/PartInfoLv2.m ↔ TermInfoLv2.m: same 3-shape + special cases.
- evaluation/PartInfoZero.m ↔ TermInfoZero.m: same zero-shape logic.
- src/verify/VerifyOmega.m + GetFeasibility.m: numerically identical drivers.
- src/autograd/: identical.
- Parameter file: K100_2.37155181.mat — 6759 float64 parameters, readable
  directly (verified by loadmat).

## Structure: SAME overall program, different lo竞技 dimension

Both rungs implement the same replacement-laser program structure (global
stage across 6 regions in Alman25 version, 3 regions in VXXZ24 version),
with VXXZ24 being the structurally simpler earlier version. The interval
machinery (interval_core.py, alman25_float.py dual-mode) works unchanged;
only the parameter-vector layout (6759 groups vs 24855) and the specific
zero-shape/level-2 term set differ.

## Estimated work

Transcription: 2-3 hours of careful work (mostly mapping 3-region p_comp).
Verification: same per-block containment + taint test procedures, applied
to the 6759-param vector.

## Decision (timebox respected)

Given remaining runway, NOT starting rung 2 tonight. Recommendation to Main:
freeze rung 1 at the tightened certified result and run rung 2 as a fresh
campaign when runway allows; the scoping notes here are the starting layout.
