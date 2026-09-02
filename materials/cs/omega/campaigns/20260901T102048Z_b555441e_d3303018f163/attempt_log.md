# Attempt log (freeze-hygiene: one entry per attempt, distinct files)

- attempt 1 (scratch/branchpin_attempt1.py): ABORTED at control C2 (exit 1).
  Registered C2 expectation ("eps == 0 under drift") mis-modelled the
  Ival.pos() hull semantics; traceback in branchpin_attempt1_stderr.log.
  No verdict produced. Timing: 2026-09-01T10:28Z.
- AMENDMENT 1 to C2 gate applied (documented in manifest §attempt history
  BEFORE attempt 2 started); no engine or interval-core code changed.
- attempt 2 (scratch/branchpin_attempt2.py == branchpin.py here): EXIT=0,
  verdict FALSIFIED (G1/G2 straddle). 2026-09-01T10:30Z.
- authoritative replay in THIS dir: EXIT=0, decision fields identical
  (branchpin_checkpoint.json). 2026-09-01T10:31Z.
