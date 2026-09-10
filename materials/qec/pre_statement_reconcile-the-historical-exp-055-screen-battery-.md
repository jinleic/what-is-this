# FREEZE — historical EXP-055 reference reconciliation

Cycle: 20260909T094859Z_2e813b. Actor: GPT-6-Astra-Main; independent reviewer: institute-opus-reviewer (Opus expensive tier, read-only).

## Frozen question

Can the complete current EXP-055 screen battery be reconciled through EXP-063's proof-preserving monotone rebind under the already-admitted Liang row-a/row-b exact-reference battery, preserving every threshold/source and revalidating every carried physical witness? Then inventory and scope, but do not search, the 21 remaining n=270 k=8 classes for a separate allocation. No repeat Liang row-b searches.

The fixed inventory is all 24 current screen shards: 3x3, 5x3, 7x3, 9x3, 9x5, 15x3, 7x7, 9x7, 21x3, 25x3, 17x5, 9x9, 15x5, 27x3, 31x3, 11x9, 33x3, 15x7, 21x5, 35x3, 13x9, 39x3, 15x9, 27x5. These contain 10,574 records. The first 22 have fingerprint 7daefcaa10015fbe7a3312710ae13a45f9d3ad8d527e407e82f658e175655125; the last two have 4dd928285865f2798e19ddb9106f1bf1f2ed5384fb73b806d75899bff710ec6d. The current fingerprint must be freshly derived from validated references, not hardcoded. The existing aggregate's scope remains unchanged: reconciliation must not promote the unresolved n=270 frontier into a closed screen.

## Acceptance

- C1: Snapshot/hash all inputs and old archive bindings before mutation. Validate admitted references including row a and row b from their existing certificates/records without launching a search. Every current shard and the existing aggregate carry the freshly derived reference fingerprint and current validator version afterward.
- C2: Strict aggregate/record validation passes on all 24 shards. Record identities, orbit totals, thresholds, threshold-source strings, scientific verdicts and carried witness supports/bounds remain identical. In particular domination_threshold(270,8) stays (20, EXP-070 [[270,8,20]] Liang row a). Every carried witness_support and ceiling_witness_support, including open/no-reference rows, is physically rechecked for exact weight, commutation and non-stabilizer membership. No lower bound is inferred from a witness or timeout.
- C3: Rebinding is repeatable within the same validator version: immutable preimages remain hash-bound, and a second identical rebind is byte-idempotent. Historical transport source bindings resolve to the exact recorded source bytes; current transported shards bind current source bytes. No archive is overwritten and no strict validator check is weakened. Repair only demonstrated protocol obstructions necessary for this objective (e.g. archive version collisions or historical transport bindings), with a failing regression and passing replay.
- C4: Fresh audit repeats the complete verification after execution. Inventory all 21 k=8 residuals by lattice/class identity with their existing bounds/supports, threshold/source, prior effort and a bounded proposed next experiment with explicit stop/acceptance criteria. Do not execute that next experiment. Preserve 147 k=12 undecided and 210 k=20/24 unpromoted survivors.
- C5: Required read-only Opus review examines primary evidence and scoped source changes before the terminal verdict. Retain its response and resolve material findings. All monitor completion receipts, commands, stdout/stderr, hashes and deviations are retained in the run before freeze/close.

## Plan, resource and authority boundary

Use campaign.py init exactly once after this file is written. All evidence goes in the minted run; existing closed campaigns, certificates and preregs are read-only. Permitted changes are EXP-063, narrowly necessary EXP-055 historical-binding validation, focused regression tests, the current screen shards/aggregate and new append-only rebind archives, target owner fields and math/PROGRESS.md plus math/qec/README.md. No institute/lifecycle/budget edits, no state status edits, no math/RESULTS.md, no services, installs or unmanaged jobs.

Stable scientific-objective task id: exp055-reference-reconciliation. First monitored start is recorded by institute_monitor.py; every later attempt preserves that original start and the fixed deadline 2026-09-09T12:00:00Z. All mini execution uses local scripts/institute_monitor.py run, run-bound, host mini-pro.local, one lane, one BLAS/OpenMP thread, real progress log/JSON nominated as outputs, stall bound 1800 seconds. Heavy verification and focused tests execute on the mini, never via an unmonitored SSH science command. No existing native search is relaunched. Projection: under 60 minutes for baseline/rebind/full validation plus under 30 minutes for fresh audit, review and termination margin; total before the fixed deadline and far below 72 hours and the allocated 36,828 seconds. If the frozen gate cannot finish, stop honestly rather than shrinking it or resetting task names/deadlines. Full measured actor wall time is recorded separately from solver/job sums.

## Exact verification command

From the target directory, under the same bounded monitor and environment OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1:

```text
.venv/bin/python -B campaigns/<minted-run>/reconcile_verify.py audit
```

The audit must independently compare current artifacts with the pinned inputs, derive the current reference fingerprint, perform C1-C4 checks and emit audit.json with all_checks_pass=true only when every check succeeds. Run it once after production and once fresh for AUDIT; logs and summary files distinguish the two invocations. Focused protocol regressions use `.venv/bin/python -m pytest tests/test_exp063_reference_rebind.py -q` plus directly affected EXP-055 screen-proof tests. No project-wide suite or new solver work.

## Verdict mapping

- FROZEN-CERTIFIED: all C1-C5 hold, with exact receipt-bound evidence and no unresolved material review finding. This certifies reconciliation, not a new code, better threshold, residual distance or n=270 closure.
- FROZEN-NEGATIVE: a verified scientific counterexample refutes the frozen preservation/witness premise. Tool failures and archive-layout defects are not scientific negatives.
- FROZEN-INCONCLUSIVE: valid partial scientific evidence exists but some scientific acceptance obligation remains unproved within the deadline; no promotion of unchecked artifacts.
- CRASHED: execution/environment/tool failure prevents producing the required scientific evidence. Administrative restriction or inability to obtain the mandated actor/reviewer maps to REJECTED, never FROZEN-NEGATIVE.
