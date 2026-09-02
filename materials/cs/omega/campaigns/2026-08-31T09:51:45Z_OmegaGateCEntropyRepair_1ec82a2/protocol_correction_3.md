# Protocol correction 3 — re-audit response (2026-08-31, static only)

Status: ZERO compute. Only a prepare-only smoke and a tiny mock-object
self-test (`--mock-self-test`) were executed; no real point evaluation and no
entropy solve over the actual VXXZ24 workspace.

Response to each rejection defect in commit `65c0cc8`:

1. **Wrapper recursion.** `wrap_class_method` now captures the original class
   methods `_ORIGINAL_PART`/`_ORIGINAL_GLOB` BEFORE assignment, and
   `patched_part_evaluate_post`/`patched_glob_evaluate_post` call those
   captured originals. The mock test proves a single invocation returns.

2. **CERT_LOG never populated.** `installed_repair_certificate` appends
   metadata exactly once per certificate (line 593). The mock self-test
   routes through `installed_repair_certificate` and asserts
   `len(CERT_LOG) == 1` plus `stored_dist_max_feasible is False`
   (i.e., the fallback path is genuinely exercised, so the log is
   nonvacuous). `evaluate_repaired` additionally asserts non-zero
   glob (>=3) and part (>=1) counts.

3. **Callback-before-R ordering.** `evaluate_repaired` installs
   `W.gm_penalty = patched_gm_penalty` and
   `W.part_penalty = patched_part_penalty` BEFORE calling
   `W.evaluate_point`; the manual R loop and the Workspace c_viol
   hashing constraints therefore read the identical
   `hash_penalty_term` field installed by the wrapper. The mock test
   asserts the installed field and `patched_part_penalty` read-back are
   endpoint-identical.

4. **feasible_primal_lower tuple order.** The status string is the SECOND
   element; both call sites in `entropy_certificate` and the plants now
   unpack `(outcome, status)` in that order and test `status`.
   `_plant_outcome` is also validated so a list outcome can never be
   confused with a status string.

5. **exact_rational_primal rewritten to the audited construction.** The
   function now takes the FULL 0/1 matrix and its exact Fractions rhs,
   augments B with a normalization row (all ones, rhs 1), selects exact
   independent rows via RREF of the augmented matrix, picks an invertible
   PIVOT column set, seeds all free coordinates with the exact dyadic
   midpoint value 1/n, solves pivot coordinates with SymPy over Fractions,
   and validates normalization + EVERY full B row + q_i >= 0 exactly. On
   success the Arb lower bound is `entropy_lower_arb(q)`; otherwise the
   caller falls back honestly to H(p). Raw interval containment is never
   accepted as feasibility; the status string is recorded in metadata and
   H(p) fallback explicitly notes "raw containment never accepted".

6. **Duplicate feasible_primal_lower.** The live definition is the single
   executable one; the byte-identical second pair remains but is now
   wrapped in the string literal `_REFACTORED_DUPLICATE_feasible_primal_lower`
   with an UNREACHABLE HISTORICAL CODE header (owner-only deletion), so no
   executable shadow is left.

7. **_entropy_arb dead helper.** Preserved inside `_DEAD_entropy_arb` string
   literal labelled UNREACHABLE HISTORICAL CODE (owner-only deletion); never
   called. The live pair is `entropy_lower_arb` (exact Fractions -> Arb) and
   `_entropy_iv` (interval -> Arb).

## Audited sequence executed

1. `.venv/bin/python -m py_compile omega/src/gate_c_entropy_repair.py` — OK.
2. `.venv/bin/python omega/src/gate_c_entropy_repair.py /tmp/mock_selftest.json --mock-self-test`
   with thread pools pinned to 1 — ran to completion WITHOUT recursion and
   produced `mock_self_test.json` with:
   - `no_wrapper_recursion: true`,
   - `cert_log_nonvacuous: true`,
   - `callback_before_R_installed: true`,
   - `projection_plants_rejected: ["exact_negative", "interval_overlap"]`,
   - `honest_fallback_used: true`.
3. No objective evaluation, no entropy solve over the real workspace, no
   `results.json`, and no decision: the campaign status remains
   HELD for owner re-audit.
