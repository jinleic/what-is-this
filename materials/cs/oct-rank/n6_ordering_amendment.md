# N6 ordering amendment — plant control before TF path

Date: 2026-09-04. Run: `20260904T052006Z_e5f2adce_0e7798f95678`.
Supersedes nothing scientifically; fixes an execution-order defect found in the
launched driver bytes (`n6_qi_real13_descent.py` SHA-256
`0b4be99ebba9d7d5323a2fe7c202b4aecc66c72dec8aad3f86519babe6a2b609`).

## Defect

`n6_prereg.md` section 2 (lines 54-66) requires: "Before TF tracking, run a
deterministic positive plant through the same continuation code", and "Any
failed positive control or false substitution accept is INVALID INSTRUMENT;
TF is not run."

The launched driver computed `target_chart = exact_chart(start_q, tf, "TF")`
before entering the plant control block. Launch 1 therefore stopped with
`INCONCLUSIVE-CHART` (TF exact realified chart rank 322/384, 0.595799 CPU s)
before the prereg-required full-path plant control ever ran.

## Disposition of launch 1

Launch-1 result `n6_results.json` (SHA-256
`c7bbe657dff2fce474a04874e72dc4f941aa6f72d9bf0c328b78c864ff0aba5d`) is a
lifecycle/launch defect record, **nonadjudicative** for N6 science: no plant
control, no continuation, and no certification ran. It is preserved byte
for byte under `attempt1_launch_defect/` together with the launch-1 progress
checkpoint, the supervisor-captured stdout (`n6_launch1.stdout`), and a
disposition note (`n6_launch1.stderr.note`) recording that no separate stderr
stream capture exists for launch 1; no attempt-1 artifact is overwritten.
The directory additionally holds hash-identical duplicate copies of the two
records created during preservation.

## Amendment

In `run()`, the TF anchor chart computation and its progress checkpoint are
moved to after the full plant control block completes successfully
(`progress("plant-control-pass", ...)`), preserving every other gate,
constant, control, and stopping rule of the bound preregistration. The plant
control now executes exactly as registered: any plant failure remains
INVALID-INSTRUMENT, and the TF path — including its chart — cannot start
before the plant passes.

## Independent singularity evidence (pre-relaunch)

A standalone auditor importing nothing from the driver
(`n6_independent_chart_audit.py`) rebuilt TF, reconstructed the frozen T5
point, and confirmed by exact Q(i) elimination: substitution 192/192 (13-term,
conjugate, 12-term), singleton corruption rejects, and
`exact_complex_rank_full = 161` (realified 322/384, kernel dimension 86).
The TF anchor singularity is intrinsic to the point — no column selection
yields a chart — so the amended relaunch is expected to reproduce a
controlled `INCONCLUSIVE-CHART` after a passing plant control, unless the
plant itself exposes a defect.
