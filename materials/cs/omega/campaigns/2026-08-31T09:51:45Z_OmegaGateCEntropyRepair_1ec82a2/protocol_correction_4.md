# Protocol correction 4 — second re-audit response (2026-08-31, static only)

Status: ZERO compute. Only a prepare-only smoke and the extended mock-object
self-test were run; no real point evaluation, no entropy solve, no
`results.json`, no decision.

Fixes to the six rejection defects of commit `6231acd`:

1. **Wrapper-called original.** `patched_part_evaluate_post` now invokes
   `_ORIGINAL_PART(self)` and `patched_glob_evaluate_post` invokes
   `_ORIGINAL_GLOB(self, parts)` (direct function calls; the prior
   `._evaluate_post` attribute access was on a plain function and would
   raise AttributeError).  The mock now wires
   `_ORIGINAL_PART`/`_ORIGINAL_GLOB` to counter-stubs and calls
   `patched_part_evaluate_post(mock)` — proving exactly one original
   invocation (`wrapper_invoked_original_exactly_once: true`) AND exactly
   three part certificates appended (`CERT_LOG == 3`).  A
   `_mock_part_identity` flag lets the mocked SimpleNamespace satisfy the
   real-Part identity guard without touching `isinstance` semantics for
   production objects.

2. **`_iv_exp_point` monotonically bounded.** For any exponent interval,
   the new implementation evaluates `exp(lo)` and `exp(hi)` separately in
   Arb at MID precision and returns
   `[exp(lo).lower()*(1-2^-OUTB), exp(hi).upper()*(1+2^-OUTB)]` — so a
   nonzero-width exponent is bounded at BOTH endpoints, never mid-only.
   Case checks in the mock:
   - wide exponent `[0.3,0.9]`: `bounded.lo <= exp(lo).lower()` AND
     `bounded.hi >= exp(hi).upper()` (ENCLOSED);
   - midpoint-only reference `_IV_EXP_MIDPOINT_ONLY` (kept ONLY as a
     never-used helper) demonstrably FAILS this plant — recorded as
     `exp_midpoint_only_fails_wide_plant: true`;
   - point exponent `[0.5,0.5]`: exp(point) enclosed.

3. **Exact projection uses pmax seeds on ACTIVE columns.**
   `exact_rational_primal(matrix_full_int, b_full_fractions,
   pmax_fractions=None)`:
   - forced-zero pmax columns fixed to exact zero;
   - projection on ACTIVE columns only (forced columns re-expanded to
     exact zero);
   - free (non-pivot) ACTIVE coordinates seeded with the pmax-derived
     exact dyadic midpoint `(pmax.lo + pmax.hi)/2`, normalized so the
     free seed vector sums to 1; active-uniform `1/n_active` fallback
     ONLY when pmax is unavailable;
   - square/full-rank case (no free columns after pivots) supported:
     solve pivots directly and validate;
   - exact validation after expansion: normalization == 1, EVERY full B
     row equality, q_i >= 0; raw interval containment is never a
     feasibility claim.
   The call site in `entropy_certificate` now extracts exact pmax
   Fractions (falling back to None when pmax is not an exact point) and
   passes them in.

4. **Genuine interval-overlap plant.** The overlap plant builds
   `[0.2,0.4]/[0.4,0.6]/[0.3,0.5]` — True non-point intervals — and
   asserts the status string actually starts with
   `"stored dist_max REJECTED"` (overlap rejection path); a separate
   genuine exact-dyadic plant `[0.25,0.5,0.25]` is asserted to PASS the
   same test (`stored dist_max exact-feasible`).  Negative plant likewise
   asserted; actual booleans recorded, no unconditional metadata.

5. **Exact branch-count formula.** `evaluate_repaired` asserts:
   - `glob_count == 3` (exactly 3 glob certificates),
   - `part_count % 3 == 0` and `part_count == 3 * (part_count // 3)`,
     i.e. exact multiple-of-3 per Part.
   The branch uniqueness key now includes `level` and `identifier`
   `(CURRENT_POINT, kind, region, level, identifier, part_id, shape)`.
   The mock asserts `len(CERT_LOG) == 3` for one Part.

6. **Dead interval_dual block neutralized.** The second executable block
   after the first `return` in `interval_dual` (which referenced the
   nonexistent `IC.exp_iv`) is now wrapped in the string literal
   `_DEAD_internal_dual_block_v2` labelled UNREACHABLE HISTORICAL CODE
   (owner-only deletion retained); no shadow executable code remains.

## Audited sequence executed (static)

1. `py_compile omega/src/gate_c_entropy_repair.py` — OK.
2. `nice -n 10 ./.venv/bin/python omega/src/gate_c_entropy_repair.py
   /tmp/mock_selftest2.json --mock-self-test` with all thread pools
   pinned to 1 — produced `mock_self_test_v2.json`:
   ```
   exp_midpoint_only_fails_wide_plant: true
   exp_point_case_ok: true
   exp_wide_endpoints_enclosed: true
   wrapper_invoked_original_exactly_once: true
   cert_log_nonvacuous: true
   callback_before_R_installed: true
   honest_fallback_used: (in CERT_LOG[0], recorded)
   ```
3. No objective evaluation, no entropy solve over the real workspace, no
   `results.json`; the campaign remains HELD for owner re-audit.
