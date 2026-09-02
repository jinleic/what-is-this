# Protocol correction 9 — v9 freeze after artifact inspection (2026-08-31)

Status: ZERO real compute (py_compile + prepare-only + mock self-test only).

## v8 artifacts that contradicted the round-6 report — all confirmed and fixed

1. **Base plant STILL positional.** asrun_v8 lines 733-734 called
   `interval_dual(matrix, b_intervals, planted, residual_matrix,
   residual_b)` — residual_matrix was consumed as mass. The v9 source and
   FROZEN ARTIFACT now call
   `interval_dual(matrix, b_intervals, planted, target_mass,
   residual_matrix=residual_matrix, residual_b=residual_b)`. Verified by
   grep on the artifact bytes (see below).

2. **Feasibility conversion still bypassed the round-trip.**
   `mass_aware_feasible` had a private
   `Fraction.from_float(float(mid_symbol))`; pmax and the second p
   conversion path bypassed `exact_dyadic_fraction`. Now centralized:
   Fraction inputs pass through; interval points route through
   `exact_dyadic_fraction` (binary64 round-trip or hard abort); non-point
   pmax is RECORDED as non-point and REJECTED (None slot + status
   `nonpoint=True`), never silently rounded.

3. **The alleged y=0 identity plant never called interval_dual.** It built
   the target formula and a manual copy of that same formula — a tautology
   that could not catch a production-algebra or signature defect. The new
   plant calls the ACTUAL `interval_dual(B3_dual, b3_iv, zeros(rows), m3,
   residual_matrix=B3_dual, residual_b=b3_iv)` and compares the returned
   interval to an independently constructed
   `m*log(n) - m*log(m)` enclosure.

## Corrected wrong-formula diagnostic (audit round-7)

   bare(-log m) - correct(-m log m) = (m-1)*log(m)

At m = 3/4: (1/4)*log(3/4) ~ 0.07192 (not the previously reported
m|log m| ~ 0.2158, which is the -m log m magnitude). The plant asserts
|wrong_gap| > 2^-6 and separation |wrong_gap| - sep_upper > 2^-8.

## Artifact self-inspection (performed on asrun_v9 BYTES)

- `grep -n 'interval_dual\(' gate_c_entropy_repair.py.asrun_v9` →
  4 hits: 1 definition + 3 call sites, ALL with `target_mass` in
  argument position 4 and keyword `residual_matrix=`/`residual_b=`
  (lines 715/753 in-source; artifact-identical).
- Conversion policy grep: `exact_dyadic_fraction(` used in
  `mass_aware_feasible._to_fraction`, the pmax loop, and production p
  conversion; the private `Fraction.from_float(float(mid_symbol))` is
  GONE from `mass_aware_feasible`.
- `mock_self_test_v7.json` verbatim:
  - `dual_identity_uses_actual_interval_dual: true`
  - `dual_identity_abs_separation_upper: 9.029685343895233e-17`
  - `wrong_formula_gap_value_m_minus_1_log_m: 0.07192051811294523`
  - `wrong_formula_gap_separates_from_certified: true`
  - `roundtrip_converter_rejects_over53bit: true`
  - `over53bit_changed_by_old_float: [1, 1]` vs
    `over53bit_true_value: [1180591620717411303425, 1180591620717411303424]`
  - `honest_fallback_used: false` (still reported as-is)
