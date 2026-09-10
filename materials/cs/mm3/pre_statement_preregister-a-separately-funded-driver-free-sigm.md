# FREEZE — H-MM3-SIGMA2-PAIRS-1

Cycle 20260910T050037Z_700a86; actor Mm3Sigma2, GPT-6 Astra only.
Exactly one new run; exclusive lane jinleic@mini-pro.local. No other host.
Stable scientific-objective monitor task H-MM3-SIGMA2-PAIRS-1 across every
chunk, retry, audit and review-driven replay; never reset its first start.

## Objective and hard scope

Decide ALL 608 already registered sigma^2 positive-pair floor instances at
T=d for the locked sixth source S1_smirnov_repo_139. Ordered keys use side
U,V,W and increasing (a,b); the certified side counts are 128/256/224.
Use the sealed pure-Python encoder and in-sandbox python-sat, clause-evaluate
every SAT model, and independently replay UP/RUP for every UNSAT. This is the
funded census itself, not a substitute registration, probe or symmetry gate.
No historical sigma^2 decisions are available as cross-instrument controls:
all decisions here need fresh evidence. The completed sigma^1 ledger and the
original 520-line sigma^0 prefix are hash-checked read-only and preserved,
never reused by symmetry or merged into this ledger.

## Cited, pinned dependencies

Workspace root is /Users/jinleic/jinleic-workspace. The source-lock dependency
is cs/mm3/campaigns/20260904T035648Z_a4bf0c37_d790a569af43/VERDICT.md.
The certified encoder/registration dependency is
cs/mm3/campaigns/20260905T100124Z_4e1069d8_da4f5d1086e0/autonomy/finalize.json.
Byte-identical copies of the following inputs and source modules are taken
from cs/mm3/campaigns/20260910T014717Z_5f1b5fb1_136d5deadc31, whose sigma^1
RESULTS.md, AUDIT.md, REVIEW.md and certificate method are cited dependencies:

- inputs/triple.json SHA256 804b0d193bd6756ed070710ee30b09f652ef7bf802c3d462245d1bfec967b951.
- Canonical triple serialization json.dumps(dict(zip('UVW', triple)), sort_keys=True): cf9d004ce93cfecfcd79f3a3adee39adc2215bf1d2b84c6977231c227338332c.
- inputs/registration.json SHA256 61b9b82db3274897442ec08ed764851342ef8a80f90cabb04c9300721ff115c6.
- sigma^2 instances_sha256 9cddb483857a1328ac03301b40a9be73de510667d85dd2fd0d0032015c868330; ordered pair_keys_sha256 c00162793e5e6afae9e7b03e86251f84609fbe2465a2a8359d90472f4137310d.
- inputs/prefix520.jsonl SHA256 0adc781ea347e10308c787240c45ca1c89b23f09b14c2d417f7865ddf7f40e09.
- sealed/encoder.py SHA256 e536bf2ab307a24a356631e079ccf1cfe1364993cd4925816b1b9af8363e5973.
- sealed/independent_encoder.py SHA256 4f95c925fb7d50ba8f394c0e2f043328ddd999a4b932d23eead6357b281649c5.
- sealed/rup_checker.py SHA256 60046cd55aff7608cad7e3ee2ba1a55e70d55845ba6f259837a6fd45c92dfb16.

Before first launch, seal all run-specific source, copied inputs (including a
preservation-only sigma1-ledger.jsonl), solver environment specification and
resource wrappers by SHA256 in seal.json. Do not import any legacy main or
restart the stopped supervisor or unchanged orbit_min_run.py disk-growth driver.

## Falsifiable acceptance

1. Both primary and independent encoders freshly reconstruct the sigma^2
   pair universe exactly as the pinned registration: 608 ordered keys and
   exact key,d,vars,clauses,DIMACS-SHA256 equality for every accepted decision.
   Sigma^2 is selected explicitly; no transfer of sigma^0/1 SAT/UNSAT values.
2. Solve each floor CNF with python-sat==1.8.dev24, Glucose3, at T=d in a
   deny-default sandbox. Pin package version and native extension hash before
   each solve chunk. Dependency setup is outside the science sandbox but
   inside its monitored finite wrapper. No external solver executable.
3. Accept SAT only after a solver-independent evaluator satisfies every CNF
   clause. Accept UNSAT only after the separate pure-Python checker derives
   contradiction by original-CNF UP or replays every RUP addition, ending in
   the empty clause. Empty raw Glucose logs do not certify anything: direct
   UP must separately succeed and its one-step empty conclusion is retained.
   Deletions may be ignored soundly; RAT-only/unverified steps fail closed.
4. Retain exact resumable JSONL ledger fields key,d,cnf_sha256,sat,lb, with
   lb=d for SAT and lb=d+1 for verified UNSAT, plus ordinal, vars, clauses,
   certificate path/hash and SHA256 chain. Genesis is SHA256 of UTF-8
   H-MM3-SIGMA2-PAIRS-1 plus newline. Each next digest hashes prior digest
   bytes concatenated with sorted compact JSON of the record without chain.
   Resume only a validated exact prefix; reconstruct tallies from bound
   records, never trust summary counters. Incomplete final lines fail closed.
5. Preserve compact SAT models or UP/RUP refutations; regenerate CNFs rather
   than retaining large DIMACS files. Replay every accepted certificate with
   the independent encoder in a separate process; then repeat the exact
   verification command fresh for AUDIT. Positive and negative model/RUP
   controls and missing-middle/reordered/edited-ledger rejection must pass.
   Hash-check the preserved sigma^0 and sigma^1 inputs on every invocation.
6. After primary evidence exists, obtain mandatory institute-opus-reviewer
   review, wait for it, and preserve model, session/agent handle, evidence
   scope and findings in REVIEW.md. No unresolved blocking finding permits
   scientific closure. The reviewer may request one bounded diagnostic only
   within this objective, original deadline, resource limits and same lane.

## Bounded execution and stopping

The coordinator did not retain an exact worker dispatch timestamp. Use the
known cycle origin 2026-09-10T05:00:37Z as a conservative objective start, NOT
an invented dispatch time. Fixed overall cap is 2026-09-10T06:00:37Z (one
hour from cycle origin), including preparation, all compute, review, audit,
transfer and termination. The earliest worker clock sampled independently
is 2026-09-10T05:04:17.455232Z. Record timing basis honestly; charge a
conservative full-wall envelope from cycle origin rather than omit earlier
preparation. No deadline reset on retry. All monitor launches bind
--task-started-utc 2026-09-10T05:00:37Z and
--deadline-utc 2026-09-10T05:57:37Z, leaving three minutes for close.

Each invocation has at most 300 seconds real wall and 8 MiB new evidence:
285-second outer timer plus bounded TERM/KILL cleanup, 235-second science
soft alarm, at most 128 next-prefix instances, 7 MiB prospective payload
allowance reserving space for blocks/logs/summaries. Existing sealed
resource_guard.py additionally enforces nice 10, one numeric thread,
CPU pause/resume 40/38%, 100 GiB free reserve, 2 GiB group RSS, 270 charged
seconds and 8 MiB sampled output growth. Pauses and dependency preparation
count against the finite outer wall timer. No new primary launch after
05:48:37Z; audit/review-driven replay only if it and termination fit before
05:57:37Z. Maximum aggregate scientific invocation wall is 2000 seconds.
Historical sigma^1 timings support a projection well below this allocation;
the full 608-instance objective is never silently replaced by a subset.
Verify live host identity/occupancy before impactful work; never overlap an
unrelated heavy predecessor. Retain exclusive lane ownership until every
planned computation and replay is terminal in fresh monitor evidence.

## Exact verification command

From the remote run directory, inside the sealed sandbox and resource guard:

`python3 sigma2_pairs.py replay --output audit.json`

The local launch for the minted run is:

`python3 /Users/jinleic/jinleic-workspace/scripts/institute_monitor.py run --domain cs --target mm3 --run <minted-run> --host mini-pro.local --user jinleic --task H-MM3-SIGMA2-PAIRS-1 --remote-cwd /Users/jinleic/institute/mm3/<minted-run> --output evidence/audit.json --output control/audit.log --task-started-utc 2026-09-10T05:00:37Z --deadline-utc 2026-09-10T05:57:37Z --stall-seconds 300 --cap-seconds 3600 -- python3 sealed/outer.py replay audit.json`

Primary uses the same launch binding and wrapper with `solve chunk-NN.json`;
the sealed sandbox command adds --limit 128. The primary replay is named
replay-primary.json; the fresh audit is audit.json. Every actual command,
return, host identity and completion receipt is retained in the run.

## Terminal mapping and nonclaims

- FROZEN-CERTIFIED: all 608 decisions satisfy every acceptance criterion,
  fresh audit passes, independent Opus review has no unresolved blocker.
- FROZEN-INCONCLUSIVE: a permitted resource/time boundary leaves incomplete
  coverage or unreplayed proof evidence; retain exact coverage and the full
  unresolved objective. Certified partial statements still require fresh
  replay and nonblocking independent review.
- CRASHED: software/environment/execution failure prevents valid completion.
  Administrative unavailability or control/model restriction may be REJECTED.
- FROZEN-NEGATIVE is not an outcome of this finite-decision gate; no operational
  obstruction or bounded incomplete coverage is scientific refutation.
- Missing/failed/unresolved vetoing review or impossible safe closure means
  UNFINISHED for coordinator recovery, never a fabricated terminal verdict.

No new addition-count record, complete orientation rows, universal no-54,
global optimum or full sixth-source orientation census follows from pair
coverage. The previous exact ladder and every closed predecessor remain
unchanged. Write PLAN.md, ITERATE.md, AUDIT.md, REVIEW.md, RESULTS.md and
literal proposed entries under exact cs/RESULTS.md and cs/PROGRESS.md headings
in LEDGER_ENTRIES.md before campaign freeze. Main alone writes those ledgers.
Freeze then close locally; never mutate the sealed run afterward. Write the
cycle usage JSON outside the run; update only permitted owner state fields,
append a target README current-state note, and run campaign.py state refresh.
Skip formatters, linters, builds and project-wide tests, not scientific gates.
This explicit workflow supersedes older local Fable/CLI delegation, commit
and institute.py instructions: only the mandatory approved Opus task review
is permitted; all control-plane commands remain local.
