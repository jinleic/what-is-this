# Revision-9 shipped-generation follow-up — hmz FREEZE

## Revision 9 — direct sampling of the four F-Z4 shipped artifacts

Written before campaign init and before any new scientific execution for cycle
`20260910T013633Z_25f67e`. Gate:
`optional-follow-up-new-cycle-would-need-funding-`.
Stable scientific objective: `zero-level-shipped-fz4-four-cell-r9`.
This is the funded follow-up, not another Revision-8 rebuilt-oracle run.
The prior preregistrations and closed runs remain immutable.

### Question and reference semantics

Are the four shipped grown F-Z4 circuits statistically compatible with the
embedded author rows, when sampled DIRECTLY with the pinned author's
postselection and accepted-only decoding semantics? The historical rebuilt
`grown@0.0006` acceptance excess (z=+3.79) motivates, but does not select or
replace, the four-cell test. Report that cell's signed acceptance residual
and z against the author row; compare to the preserved R8 row descriptively.
Agreement is finite statistical compatibility, NOT proof that these exact
files generated the historical author rows. Repeating a seeded stream does
not eliminate Monte-Carlo fluctuation as an explanation of a discrepancy.

Primary repository: https://github.com/FujitsuResearch/Zero-level_CCZ_Distillation/tree/1b59e223590492e224bd8623a4e0bcba59029e01
Commit `1b59e223590492e224bd8623a4e0bcba59029e01`; tree
`9a5d89401fd40554635fb0e9d0b4da818670c2bc`; tar size 1,072,583 bytes and
sha256 `d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`.
Use the canonical R7 retained tar (no download substitution). Preserve the
full existing member-pin table, captured tar bytes, code, prereg and reference
rows inside this new run. Pinned members must pass safe extraction and blob
checks before sampling; the same validated bytes are the archived bytes.

| grown p | shipped git blob SHA1 | author shots | LER | acceptance | seed |
|---|---|---:|---:|---:|---:|
| 0.0008 | 27fa0d4283ee05c991098e0d4d3a9efc3854d274 | 5000000 | 0.0002219381364258208 | 0.3568562 | 2202 |
| 0.0006 | ea17f178c42b46c117d6826ddd87c14b861fdba7 | 5000000 | 0.00011960131162771752 | 0.4615334 | 2203 |
| 0.0004 | 3a28a475a9f8893e55638df97ccebba96c81125c | 5000000 | 5.422499759000011e-05 | 0.5975104 | 2204 |
| 0.0001 | cecca6083c5575da45071615fcef6a5001356804 | 10000000 | 3.866977783985162e-06 | 0.8792396 | 2206 |

For each row, open the pinned shipped `.stim` bytes, verify their blob and
sha256, instantiate that exact circuit, and use THAT object for both
`detector_error_model(decompose_errors=True)` and
`compile_detector_sampler(seed=seed)`. Never replace it with a rebuild.
Record raw-byte, canonical/flattened circuit, DEM and mask hashes. Rebuild
only to recover `Stim_builder.postselct_numbers()` from the pinned driver
prefix (parameterizing only its error-rate assignment), to record factual
flattened inequality, and to check detector/observable index definitions.
Require three observables/fault IDs and detector count matching the mask
length; compare the ordered DETECTOR and OBSERVABLE_INCLUDE definitions and
archive their signatures. The mask's reference is explicitly the pinned
reader semantics, not an inferred historical generator.

The upstream driver calls `output_circuit_text` BEFORE reading its file.
This amendment deliberately does NOT execute output methods: doing so would
overwrite the shipped generation with the rebuilt generation and defeat
the funded question. The matching historical generator is not recovered.
Full-circuit rebuild equality is neither a premise nor a success criterion;
the expected F-Z4 inequality is recorded truthfully. The rebuild is never
sampled for these four primary or replay cells.

Accepted iff no selected detector fires. Reject shots BEFORE decoding;
`decode_batch` sees accepted syndromes only. Error iff any of the three
predicted observable bits differs. Chunk size exactly 10,000; a single
continuous seeded sampler per cell; do not reseed chunks. One primary row
per cell, fixed seeds above, no optional extra shots or seed shopping. Save
complete ordered chunk counts and syndrome+observable digests. A completed
row is never resampled as primary; a mismatch in its identity is a refusal.

### Ordered verification and scientific gates

1. Validate environment pins: stim 1.16.0, pymatching 2.4.0, numpy 2.5.2;
   record Python, machine/platform and host. Preserve campaign manifest and
   prereg hash binding, source pins, code snapshots and live-run checks.
2. Execute the complete existing eight-gate smoke battery first: builder-only
   ungrown p=0 zero-error/invariant/replay tests and genuine shipped grown
   p=0.001 equality/fault-id/accepted-only scalar-batch/replay tests.
3. Before primary cells, additionally sample each of the FOUR direct shipped
   split cells for 20,000 smoke shots, chunk 5,000, seed 3302; require exact
   accepted-row scalar/batch agreement and same-seed ledger replay. These
   are readiness evidence, never added to primary author-count estimates.
4. Execute all four primary cells at the exact table budgets and fixed
   schedules. No substitution of an easier smoke or subset objective.
5. After primary evidence exists, obtain independent `institute-opus-reviewer`
   review of source binding, masks/decoding, counts, statistical contract,
   scope and planned closure. Preserve model, agent/session handle, reviewed
   paths and findings in REVIEW.md. Missing/failed/unresolved veto blocks
   scientific closure. No scientific verdict is chosen before that barrier.
6. Fresh audit re-extracts the pinned tar, revalidates every row and code/input
   identity, recomputes ALL eight author comparisons and the descriptive R8
   signature, re-executes the smoke battery, and fully replays ALL FOUR
   primary cells with the original seeds/chunks on the same mini. Require
   exact counts and EVERY chunk digest. Replay is verification, not new
   independent statistical evidence. Review-driven verification stays under
   the original objective/deadline and cannot change the primary schedule.

All eight comparisons use the unchanged combined analytic binomial sigma
`hypot(sigma_author, sigma_repro)`, with LER sigma
`sqrt(l*(1-l)/(shots*acceptance))` and acceptance sigma
`sqrt(a*(1-a)/shots)`. These are the pinned `std_calc.py` formulas. Retain
`|z| <= 3.53`, the upward-rounded dependence-valid Bonferroni 1% bound for
the original 24-comparison family, conservatively rather than relaxing it
for the eight new comparisons. No c-value, fit/slope or signature-only gate.
The z criterion is the inherited approximate binomial-normal instrument;
selection based on prior R8 evidence and reuse of author rows mean no new
claim of a single combined sequential familywise error guarantee.

### Verdict mapping (only after independent review and fresh audit)

With all four primary cells complete, clean source/sampling provenance,
all required readiness/replay checks passed, and independent review resolved:
- any defined |z| >= 5 OR at least three values beyond 3.53:
  `FROZEN-NEGATIVE`, refuting compatibility under THIS four-cell instrument;
- otherwise undefined statistics or one/two mild failures:
  `FROZEN-INCONCLUSIVE`;
- all eight defined comparisons |z| <= 3.53: `FROZEN-CERTIFIED`, certifying
  four-cell finite statistical agreement only, not historic generation.

A source/blob/mask/definition/DEM/decoder/smoke/integrity failure maps to
administrative `REJECTED`; unavailable execution, deadline/scheduling
obstruction or crash maps to `CRASHED` (or administrative `REJECTED` when
appropriate), never a scientific negative. Missing or failed required review
leaves evidence preserved and scientific closure blocked. Freeze and close
remain local `campaign.py` operations after closing artifacts exist.

### Scope, bounded placement and exact commands

Exclusive lane `jinleic@mini-pro.local`; NEVER mini-0. One low-priority process
at a time, thread environment caps one, `nice -n 10`. All scientific
execution (including arithmetic and replay) goes through local
`scripts/institute_monitor.py run`. No services, raw nohup or hand-written
monitor records. Scientific objective origin is fixed to the first measured
worker timestamp `2026-09-10T01:44:12Z`; all attempts use fixed mini deadline
`2026-09-10T04:15:00Z`. All preparation, science, review, audit, termination,
freeze/close and accounting must finish by `2026-09-10T04:45:00Z`, inside
current target allocation `seconds_max=14400`, `seconds_used=2522.423` and
one remaining funded run. The separate measured accounting origin is the
coordinator-observed batch dispatch `2026-09-10T01:43:58.643Z`. No retry,
phase, host or name resets any window. Before launch, confirm the full
remaining primary + replay + review + termination projection fits. Historical
R8 wall evidence supports reserving up to 3,600 seconds each for four-cell
primary and full replay, plus preparation/review/closure within these bounds;
if that cannot fit, preserve an obstruction, not a narrowed gate.

Let RUN be the single campaign.py-minted id and REMOTE be
`/Users/jinleic/institute/qlops-r9/physics/qlops`. Upload only this run's
manifest/plan/prereg, source tar, reference R8 evidence, and source snapshots.
The exact science command, launched from the workspace root, is:

```sh
python3 scripts/institute_monitor.py run --domain physics --target qlops --run RUN --host mini-pro.local --user jinleic --task zero-level-shipped-fz4-four-cell-r9 --task-started-utc 2026-09-10T01:44:12Z --remote-cwd /Users/jinleic/institute/qlops-r9/physics/qlops --output campaigns/RUN/primary-progress.jsonl --output campaigns/RUN/results/analysis.json --deadline-utc 2026-09-10T04:15:00Z --stall-seconds 600 --cap-seconds 10800 --json -- env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /usr/bin/nice -n 10 uv run --no-project --python 3.13 --with stim==1.16.0 --with pymatching==2.4.0 --with numpy==2.5.2 python src/zero_level_shipped_r9.py primary --run-dir campaigns/RUN --source-tar source_tar.gz
```

The exact fresh verification command uses that identical monitor binding,
origin, deadline, caps and environment, with outputs
`campaigns/RUN/audit-progress.jsonl` and `campaigns/RUN/results/audit.json`,
and final argv:

```sh
python src/zero_level_shipped_r9.py audit --run-dir campaigns/RUN --source-tar source_tar.gz
```

Fetch real completion receipts and progress and compare transferred hashes.
Heartbeat the local campaign during compute; freshly verify every monitor
attempt terminal before final lane release. Write PLAN, ITERATE, REVIEW,
AUDIT and RESULTS before freeze. Physics owner directly writes RESULTS.md
and PROGRESS.md with full-entry exact deduplication; no other domain edits.

Historical Gate A (81 checks reproduced) and Gate B (magic-state
comparability falsified) are unchanged, not rerun. This run does not certify
all twelve R8 cells or reinterpret any prior verdict. The executable grown
artifact is d=9; paper d=7 and physical c~300 remain NOT-REPRODUCED. No d=7
reconstruction, external generator search or second campaign is authorized.
