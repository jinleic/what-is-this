# Protocol correction 5 — third re-audit response (2026-08-31, static only)

Status: ZERO compute. Only `py_compile`, `--prepare-only`, and the extended
mock self-test ran. No real point evaluation, no entropy solve over the
VXXZ24 workspace, no `results.json`, no decision.

## Retraction of an overstated claim (audit round-3 item 6)

My previous report claimed `honest_fallback_used: true`. **That was wrong.**
The frozen artifact `mock_self_test_v2.json` recorded
`honest_fallback_used: false`, and the new run records the same. The mock's
`dist_max` is exactly feasible, so the **fallback path was NOT exercised** by
the mock; the exact-feasible path was. The narrative is corrected here and
the flag is reported verbatim from the artifact. No claim is made that the
H(p) fallback has been exercised in a test.

## Fixes to the six defects of commit `84cdebd`

1. **Projection indexing bug.** `seed_vector` was indexed by ACTIVE column
   id inside `for j in free` although it had FREE-position length. The
   rewritten `exact_rational_primal` builds a per-FREE-POSITION `seeds`
   list and consumes it via `enumerate(free)` / `seeds[position]` in both
   the RHS reduction and the expansion, eliminating the index mismatch.

2. **pmax midpoints now always used; no free-seed renormalization.** The
   caller builds `pmax_seeds` for EVERY coordinate as the exact dyadic
   midpoint `(pmax.lo + pmax.hi)/2` of the finite pmax interval —
   regardless of width — and passes them in. The projection assigns each
   free coordinate its own pmax seed value and does NOT renormalize the
   free seeds; the pivot coordinates carry the normalization constraint
   (verified by the exact `sum(q) == 1` check after the solve). Seeds are
   explicitly documented as seeds, never a feasibility claim.

3. **Structural forced set passed explicitly.** `exact_rational_primal`
   now takes `forced_columns` and pins exactly those coordinates to 0.
   The caller passes the structurally proven `forced_columns` (induced by
   zero target marginals), not a pmax-zero heuristic.

4. **Production overlap plant replaced.** The production plant now builds
   genuine non-point intervals (`centre ± 2^-20` per coordinate, with an
   assertion that at least one is non-point), ASSERTS
   `status_ov.startswith("stored dist_max REJECTED")`, and records the
   real booleans plus both status strings in metadata. The old
   clip-of-`p` construction and hardcoded `True` are gone.

5. **Exact derived counts.** `prepare()` now computes
   `EXPECTED_PART_CERTS = 3 * sum(isinstance(t, F.Part) for level in
   ws.parts for t in level)` and `EXPECTED_GLOB_CERTS = 3` from the
   CONSTRUCTED Workspace, persists both in `protocol_guards`
   (`real_part_objects_in_workspace: 63`,
   `expected_part_certificates: 189`, `expected_glob_certificates: 3`),
   and `evaluate_repaired` asserts exact equality for both at every
   point. The divisibility tautology is removed. The mock now ALSO
   invokes the GlobalStage wrapper with a stub original, so the "+3 glob"
   half of the claim is tested rather than asserted.

6. **Narrative corrected** — see the retraction above.

## New direct projection plants (audit closing request)

All executed in the mock, all asserted:

| plant | construction | result |
|---|---|---|
| noncontiguous free indices, nonuniform pmax | `B=[[1,1,0,0],[0,0,1,1]]`, pivots `{0,2}` → free `{1,3}`; pmax `[1/16,7/16,3/8,1/8]` | succeeds; `q=[1/16,7/16,3/8,1/8]`; sum==1 and BOTH full rows exact |
| infeasible pmax repaired | pmax `[1/16,7/16,1/2,1/2]` (sums 3/2) | succeeds; exact normalization + all full rows |
| negative pivot | pmax `[0,9/16,0,1/2]` forces `q0=-1/16` | REJECTED (`ok=False`, `q=None`) |
| structural forced column | `forced_columns=(3,)` | succeeds with `q[3]==0` exactly; sum==1; all full rows exact |

## Mock self-test result (verbatim from `mock_self_test_v3.json`)

```
exp_wide_endpoints_enclosed: true
exp_midpoint_only_fails_wide_plant: true
exp_point_case_ok: true
wrapper_invoked_original_exactly_once: true
glob_wrapper_invoked_original_exactly_once: true
mock_part_certificates: 3
mock_glob_certificates: 3
cert_log_nonvacuous: true
callback_before_R_installed: true
honest_fallback_used: false      <-- fallback NOT exercised; reported as-is
interval_overlap_plant_rejected: true
interval_overlap_plant_status: "stored dist_max REJECTED (overlap=True negatives=False)"
exact_negative_plant_rejected: true
projection_noncontiguous_free_ok: true
projection_repairs_infeasible_pmax: true
projection_rejects_negative_pivot: true
projection_forced_columns_pinned_zero: true
```

No objective evaluation and no entropy solve were performed; the campaign
remains HELD for owner re-audit.
