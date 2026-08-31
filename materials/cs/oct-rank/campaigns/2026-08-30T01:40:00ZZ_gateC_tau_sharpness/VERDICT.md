# Gate C frozen verdicts — tau R_R = 7 re-verified; sharpness C2 done

UTC freeze: 2026-08-30T01:40Z. Contents:
- `gate_c_tau_upper.py` — independent python-flint (fmpq + arb outward) replay
  of the upstream tau_r7 Krawczyk certificate.
- `gate_c_lower_rank4.py` — exact flattening/slice scaffolding (NOT a CP lower
  bound; scope stated in its docstring and in pre_statement.md C1(a)).
- `gate_c_sharpness_n2_n4.py` — exact CP upper witnesses for T_C (rank 3) and
  T_H (rank 8); the matching lower bounds are the paper's pencil argument,
  machine side = Lean replay (Gate A artifact) [CITED-DEPENDENCY].

## Verdicts (all MACHINE-VERIFIED in this artifact's scripts)

1. R_R(tau_H) <= 7 — Krawczyk certificate re-verified end-to-end in exact
   rationals, all three hypotheses PASS with outward-rounded arb confirmation:
   - ||Y J_S - I||_inf = 6.124e-15 (exact rational, < 1e-6 consistency gate)
   - K = 0.005351823996997164 (exact: 9777124417.../1826877046...) < 1
   - Banach margin = 9.946e-06 > 0 (= (1-K)*rho - off_max at rho = 1e-5)
   - worst per-equation margin = 9.946e-06 > 0
   - matches paper's announced ~9.9e-6 margin. [REPRODUCED]
   - Evidence label: MACHINE-VERIFIED re-verification; mathematics remains
     CITED-DEPENDENCY of arXiv:2608.16649 (+ its Lean artifacts).

2. Sharpness of (5/2)n - 2 at n = 2: exact 3-term CP witness for T_C
   (Gauss identity, integer entries) — upper bound 3 exact.
   Lower bound 3 = pencil Thm 2 at n=2 [CITED-DEPENDENCY; Lean replay].
   => R_R(T_C) = 3 confirmed (upper here, lower via Lean).
   [REPRODUCED] for the witness; label MACHINE-VERIFIED for the entrywise
   identity; the EQUALITY statement is conditional on the CITED lower bound.

3. Sharpness at n = 4: exact 8-term CP witness for T_H (Hadamard 4-term +
   sparse 4-term scheme; factors in {0, ±1/4, ±1, ±2}) — upper bound 8 exact.
   Lower bound 8 = peel(2) + pencil(6) [CITED-DEPENDENCY; Lean replay].
   => R_R(T_H) = 8 confirmed (same conditional structure).
   [REPRODUCED]; witness identity MACHINE-VERIFIED entrywise in fmpq.

4. Flattening scaffolding (explicitly NOT a bound): slices of tau_H have
   matrix rank exactly 4, pairwise non-proportional; all 128 homogeneous
   sign patterns have flattening rank exactly 4 (no sign-flip drops the
   necessary condition below 4). Scope honesty: these are necessary
   conditions only; no CP lower bound is claimed from them.

## How to run

    /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c_tau_upper.py
    /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c_sharpness_n2_n4.py
    /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c_lower_rank4.py

Certificates read from ../scratch/upstream_ref/certs/tower_cert/tau_r7
(symbolic link inside scratch/upstream_ref/verify/certs/). Environment:
python 3.14.3, flint 0.9.0, single-threaded; no randomness (deterministic).

## Also stored here

- `gate_b_rank24_results.json`, `gate_b_rank24_best.npz` — copy of the
  completed random-start rank-24 sweep summary (canonical copy lives in the
  gate B ladder campaign; see campaigns/2026-08-30T01:42:00ZZ_gateB_ladder).
