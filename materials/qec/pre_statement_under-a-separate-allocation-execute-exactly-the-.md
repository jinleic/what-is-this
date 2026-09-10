# FREEZE — exact 21-record n=270 k=8 residual witness search

Cycle: `20260910T001346Z_7c9821`. Actor: top-level GPT-6 Astra. Required independent reviewer: `institute-opus-reviewer` (Opus expensive tier), read-only, after primary evidence. No delegated science actor and no new institute cycle.

## Gate and input authority

Can the existing `experiments/exp069_n270_frontier.py::deepen_record(record, tries=20000)` supply independently verified weight-at-most-20 logical witnesses for **all 21** exact records in `campaigns/20260909T095538Z_4809c2f8_816c74c3bcbd/RESIDUAL_SCOPE.md` and `residuals.json` under this bounded allocation?

Input `residuals.json` SHA-256: `554a3704c125a003aa8eca83d7e21490d655299071d0ede5f685871a8f5e851e`, verified against its closed run pin. Its 17 records on 15x9 and four on 27x5 are the complete selection, in stored order. Preserve every identity, A/B support, orbit, k=8, n=270, threshold 20 and source `EXP-070 [[270,8,20]] Liang row a`. Prior weight-22 witnesses are upper bounds, not exact distances.

Pinned existing producer SHA-256: `02ecb6817eb10804c0665017bd4fe492e97d8b633b01208562c47466c6707d8e` (EXP-069), `e4a0e63f3faae5745ba553ac633be5b334fd764b327c3f2b7b7b9f0c9e06a60a` (EXP-055), `4a1eda59a2f162a87d77998bdb672e3f872e14535b651585e693cdb81b181e29` (EXP-056). Pin the copied source dependency tree before launch. Do not modify production algorithms.

## Resources and stop rules

One mini-pro.local science lane, one BLAS/OpenMP thread. Stable scientific-objective task id: `n270-k8-residual21-deepen20000`; reuse it for all validation/recovery attempts, never reset its window or relaunch completed search. Original allocation starts **2026-09-10T00:14:12Z** and ends **2026-09-10T01:14:12Z** (3,600 seconds including input preparation, verification, review and termination). The earlier initial brief-reading interval is included separately in full cycle usage if its dispatch timestamp is available. Monitor hard deadline is tightened to **2026-09-10T01:13:12Z**, leaving 60 seconds for monitor/watchdog termination. Every remote science or audit command must run through local `scripts/institute_monitor.py run`, with this task/run binding, original start/deadline, real progress outputs and explicit 180-second stall bound.

At most one `deepen_record(tries=20000)` call per selected record, with a 120-second whole child-process cap including import/build/verification; no retry of timed-out or completed records. The seed is the existing deterministic identity/threshold seed, hence up to the first 2,000 randomized orders repeat the previous search. This is **not 20,000 independent new tries**; at most 18,000 additional orders per record, with early target stopping allowed by the existing algorithm. Requested tries must not be reported as measured completed tries.

Historical linear projection is 1,569.484 seconds for 21 calls, not a benchmark. Before launch require remaining time to contain that estimate plus 420 seconds audit/review and 60 seconds termination. After the first three records, use 1.25 times the maximum observed whole-child duration (at most the 120-second cap) for remaining-record projection, plus 420 seconds audit/review and 60 seconds termination. If the full remainder cannot fit the original deadline, stop search and retain unattempted obligations as undecided; do not change the gate. Before each further record require its full 120-second slot plus the same reserves. A physical-check counterexample or structural input mismatch stops the batch for audit. Operational errors are not scientific refutations.

## Acceptance and exact verification command

Store each call's full input identity, deterministic seed, returned supports/bound, measured duration, outcome and physical verification; write real append-only progress on record transitions and the three-record projection. Validate all original weight-22 supports before searching. Independently reconstruct the checks and verify n/k, support uniqueness/range, actual weight, zero syndrome, non-stabilizer membership and weight-preserving X/Z transport using the existing NumPy and bitset GF(2) paths. Re-run verification fresh after search without redoing randomized search.

Exact science replay command (cwd the run directory, with the run-local pinned `source/` tree):

```
<target>/.venv/bin/python -B residual_campaign.py verify
```

`<target>` is `/Users/jinleic/jinleic-workspace/math/qec` on mini-pro. The command writes `verification.json`; the same argv with `--output verification-fresh.json` supplies the fresh audit. Both are launched via the local monitor using the frozen stable objective and deadline. Verification must reject duplicate/missing identities, threshold/source drift, duplicate/out-of-range or invalid supports, bound/weight mismatch and false domination. Exercise in-memory negative controls; no extra search. Hash-bind inputs, source, outputs and completion receipts. Script evidence belongs to this experiment, not a new public API or project-wide test claim.

A record is dominated only when a physically reverified upper bound is <=20. Larger bounds, timeouts, no witness or deadline-deferred records stay undecided. Whole-gate acceptance requires exactly all 21 dominated, no missing obligation, and independent Opus review without unresolved correctness blockers. Verified partial bounds remain useful but do not certify the gate.

## Verdict mapping

- `FROZEN-CERTIFIED`: all 21 exact obligations dominated by independently verified weight<=20 witnesses, all frozen gates verified and independent Opus review accepted.
- `FROZEN-INCONCLUSIVE`: sound bounded execution/partial physical evidence, but at least one of the 21 lacks an accepted weight<=20 witness (including timeout/deferred). No lower bound or scientific negative follows.
- `FROZEN-NEGATIVE`: only an independently established scientific refutation of the frozen physical premise/gate. Search failure or witness weight>20 never qualifies.
- `CRASHED`: operational unavailability, failed launch, tool/model restriction, execution failure or inability to obtain required verification/review. Preserve exact run and evidence, do not retry completed search.
- Administrative `REJECTED`: admission/policy refusal or invalid evidence that cannot be repaired inside this exact gate and deadline. Never substitute an easier gate.

## Scope and publication boundaries

No `resolve_screen`, no new initial screen, no Liang rerun, no SAT/MIP/native solver or exact lower-bound work, no k12 records, no k20/k24 survivor promotion, no threshold/source change, no manuscript/package edits, no publication. Existing 147 k12 undecided, 210 unpromoted survivors and closed n<=234 result are unchanged. Retain search results run-locally; do not rebind or mutate canonical screen/protocol artifacts in this allocation. Owner state/ledgers may report independently verified run-local progress with that explicit boundary. All frozen/closed prior runs remain byte-unchanged.

Only target-local prereg/run evidence, permitted target owner fields, math/PROGRESS.md, target README and the exact cycle usage record may be written. Single math target: coordinator writes both ledgers with full-entry exact deduplication. Do not modify institute lifecycle/budget/priority, monitors by hand, or orchestration plans. Close only after PLAN, ITERATE, fresh AUDIT, independent Opus review and RESULTS exist; freeze then close via campaign.py; never write into the sealed run.
