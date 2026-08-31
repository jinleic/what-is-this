Campaign manifest — gate C4 (peeling-constant probe) + C1(a)/C2 lower halves

Frozen: 2026-08-30T12:15:25Z (UTC timestamp in directory name)
Agent: OctRankGateC
Machine: arm64 (Apple M3 Ultra), darwin

Tool versions (exact):
  python 3.14.3 (cs/.venv/bin/python)
  python-flint 0.9.0 (fmpq, fmpq_mat; exact rational arithmetic)
  sympy 1.14.0 (exact symbolic expansion; Sturm sequences)
  BLAS threads pinned: OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1
  No floating-point arithmetic anywhere in the three new scripts.

Pre-registration:
  pre_statement.md, section "Gate C addendum C4/C3 — n=8 peeling-constant
  probe [DECIDED AHEAD, 2026-08-30]" — committed BEFORE the probe ran.
  The probes S1 (pencil floor), S2 (stopping points j=0..6), S4
  (kernel-dimension tightness) and the NOT-SWEPT gap S3 are pre-registered
  there, including the pass/fail rule.

Contents (all FILES copied as run; outputs are the raw stdout+exit codes):
  gate_c_tau_lower_exact.py        subtask 1 lower half: exact symbolic
                                   re-derivation (coefficient-wise, no
                                   floats) of every polynomial identity the
                                   tau chain consumes; verdict PASS.
  gate_c_sharpness_lower_exact.py  subtask 2 lower half: self-contained
                                   exact machine proof R_R(T_C) >= 3
                                   (determinant multiplicativity +
                                   coefficient matching + Sturm zero-root
                                   certificate for t^2+1); verdict PASS.
  gate_c4_peel_n8.py               THE C4 PROBE: exact 12-term witness for
                                   (I_8, J_8), stopping-point arithmetic,
                                   key_bound tightness; verdict NEGATIVE
                                   (exit 0 as pre-registered).
  gate_c_sharpness_n2_n4.py        prior gate C artifact (upper witnesses,
                                   n=2/n=4), copied for self-containment.
  gate_c_tau_upper.py              prior gate C artifact (tau Krawczyk
                                   replay), copied for self-containment.
  *.out                            raw run outputs as frozen.
  code_hashes.sha256               SHA256 of all five scripts.
  manifest.md                      this file.

Headline results (evidence labels per repo contract):
  1. R_R(tau) = 7 — upper: MACHINE-VERIFIED (prior campaign, re-run PASS
     here); lower: chain identities re-derived exactly here
     (MACHINE-VERIFIED) + Lean replay [CITED-DEPENDENCY, gate A].
  2. (5/2)n - 2 sharp at n=2 and n=4: n=2 fully self-contained
     (MACHINE-VERIFIED both directions: exact 3-term witness + own exact
     lower-bound proof); n=4 upper MACHINE-VERIFIED (8-term witness), lower
     peel(2)+pencil(6) [CITED-DEPENDENCY: Lean replay].
  3. C4 probe: NEGATIVE, all three probes MACHINE-VERIFIED exact:
     (S1) pencilRank == 12 exactly on the consumed class;
     (S2) lower == 18 at every stopping point j = 0..6;
     (S4) every key_bound inequality exactly tight on the 12-term witness.
     Rule-7 scope in gate_c4_peel_n8.py "rule7_scope" and in VERDICT.md.
     No improvement of 18 exists inside the peeling+Thm-2 machinery.
  S3 (3-slice residual rank >= 14?) NOT swept — open, would improve 18.

Anchor-first discipline honored: the 3- and 8-term upper witnesses and the
tau Krawczyk replay were re-run and PASSED before any new computation was
trusted; the C4 probe reuses the campaign-verified 2x2 witness as its base
case (verified entrywise in fmpq inside gate_c4_peel_n8.py as well).

No README/RESULTS/PROGRESS edits by this campaign other than the appended
"Current state" section in cs/oct-rank/README.md (append-only contract).
