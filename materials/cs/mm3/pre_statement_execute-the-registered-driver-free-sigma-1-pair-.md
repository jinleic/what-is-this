# FREEZE — H-MM3-SIGMA1-PAIRS-1

Cycle: 20260910T013633Z_25f67e. Actor: MM3Sigma1, GPT-6 Astra only.
Preparation accounting starts 2026-09-10T01:44:05Z; absolute campaign deadline
2026-09-10T02:44:05Z (3600 s including preparation, every computation, review,
audit, transfer and termination). One run; one lane: jinleic@mini-0.local.
Stable monitor task: H-MM3-SIGMA1-PAIRS-1 for every attempt/replay.

## Frozen question and scope

Execute the registered driver-free sigma^1 pair census for the locked sixth
source S1_smirnov_repo_139: ALL 608 positive pair instances in the registered
ordered sigma^1 universe, at T=d, using the sealed pure-Python encoder and
in-sandbox python-sat. This is not a substitute resource/probe gate. All 52
retained sigma^1 decisions (attempt-3 lines 609–660) are mandatory independent
instrument controls. Preserve the separate 520-line sigma^0 attempt-4 prefix
byte-for-byte; it is not reused as new sigma^1 evidence or promoted.

Prior input directory (workspace absolute):
/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/20260904T035648Z_a4bf0c37_d790a569af43
Inputs: sixth_source_triple.json; branchA_pair_checkpoint.attempt3.jsonl;
branchA_pair_checkpoint.jsonl. File SHA256s respectively:
804b0d193bd6756ed070710ee30b09f652ef7bf802c3d462245d1bfec967b951;
1b1f54ec1baede85bf213cb926fbeb9188b6110ca0fa0ba34292d0698d3f699a;
0adc781ea347e10308c787240c45ca1c89b23f09b14c2d417f7865ddf7f40e09.
Canonical triple hash (json.dumps U/V/W, sort_keys=True):
cf9d004ce93cfecfcd79f3a3adee39adc2215bf1d2b84c6977231c227338332c.

Encoder dependency: byte-identical encoder_fidelity_gate.py from
/Users/jinleic/jinleic-workspace/cs/mm3/campaigns/20260905T100124Z_4e1069d8_da4f5d1086e0
SHA256 e536bf2ab307a24a356631e079ccf1cfe1364993cd4925816b1b9af8363e5973.
Use its pure functions, never its old main or the legacy orbit driver.
Fresh replay also uses the byte-identical independent
review/review_encoder_replay.py from that run, and the frozen registration
JSON as data. The latter independently re-encodes every accepted instance.

## Falsifiable acceptance

1. Regenerate the sigma^1 positive-pair universe in side U,V,W and increasing
   (a,b) order; exactly 608 entries (224/128/256), matching the pinned
   registration's key,d,vars,clauses,CNF hashes. No inferred sigma symmetry.
2. Solve at T=d using python-sat Glucose3 with proof logging in a deny-default
   sandbox; pin the installed python-sat version and native module hash before
   solving. Dependency preparation may use uv outside the science sandbox.
   No external SAT executables, legacy supervisor or unchanged disk-growth code.
3. Evaluate every SAT model against every regenerated CNF clause. UNSAT is
   accepted only if a separate pure-Python UP/RUP checker replays a refutation:
   direct unit-propagation to contradiction, or each added proof clause by
   negating it and unit-propagating against the active clause database,
   including an empty-clause conclusion. RAT-only steps are not accepted.
   Deletion steps may be ignored (retaining proven clauses is sound).
4. Append an exact resumable JSONL ledger with mandatory
   (key,d,cnf_sha256,sat,lb), where lb=d for SAT and d+1 for checked UNSAT;
   bind each line to ordinal, model/proof evidence and a SHA256 chain. Resume
   only the validated longest exact prefix; recompute all tallies from lines,
   not unbound checkpoint counters. Disagreement with any retained control
   blocks promotion and stops solving. Never count an unverified result.
5. All 52 controls must agree in key,d,CNF hash,sat,lb; all original prefix
   bytes remain intact. Fresh replay must reconstruct the ordered universe,
   re-encode accepted CNFs with the independent encoder, re-evaluate models and
   independently replay UNSAT evidence. Include positive and negative SAT/RUP
   checker controls plus missing/reordered/edited-ledger rejection controls.

## Bounded execution and stopping

Each supervised compute invocation has at most 300 s real wall and at most
8 MiB newly written evidence, with a soft stop before the boundary so complete
records can be retained. Use no more than 2000 s of scientific invocation wall
in total; reserve remaining wall for preparation, independent review, fresh
audit and close. Stop launching science after 02:32:05Z unless a bounded audit
or review-driven replay fits before 02:41:05Z with termination/close reserve.
The entire 608-instance objective remains registered if coverage is partial.
Chunks use deterministic next-prefix order, at most 128 pairs each, with
pre-write byte accounting and no retained full DIMACS files (regenerate by
hash). Proof/model objects are retained compactly; never delete historical
outputs to recover space. All scientific computations run on assigned mini
through local scripts/institute_monitor.py run, original task start/deadline
binding retained. Enforce nice 10, all numeric thread variables 1, the existing
CS resource_guard.py CPU pause/resume 40/38%, 100 GiB free reserve, 2 GiB group
RSS cap, 8 MiB guard growth limit. A finite outer 300 s timer includes guard
admission and pauses. No concurrent heavy job on mini-0; no other host.

## Exact verification command

From the remote run directory, under the sealed sandbox and resource guard:
`python3 sigma1_pairs.py replay --output audit.json`

The local launch is bound to the minted run, with task
H-MM3-SIGMA1-PAIRS-1, --task-started-utc 2026-09-10T01:44:05Z and
--deadline-utc 2026-09-10T02:41:05Z; the sealed invocation wrapper additionally
caps each invocation at 300 s. Primary chunks use
`python3 sigma1_pairs.py solve --output chunk-N.json --limit 128`.
Paths, source/solver hashes, actual argv, per-attempt receipt, and all command
outcomes are recorded in the run. The fresh replay is a new read-only pass over
primary evidence, not a re-solve or new objective.

## Mandatory review and verdict mapping

After primary evidence exists, request the independent institute-opus-reviewer
(expensive-tier approved Opus) and wait before selecting a scientific verdict.
Preserve its exact identity/handle, evidence scope and findings in REVIEW.md.
No missing, failed, or unresolved vetoing review permits scientific closure.
Resolve a repair only inside this fixed objective, original wall and lane.

- FROZEN-CERTIFIED: all 608 finite pair decisions satisfy every criterion,
  all 52 controls agree, fresh replay passes, independent review has no
  unresolved blocker. This certifies only this sigma^1 pair census.
- FROZEN-INCONCLUSIVE: permitted resource/time boundary leaves fewer than 608
  verified decisions, or bounded proof replay cannot establish every decision;
  exact retained coverage and unresolved objective are reported. Any retained
  verified subresult requires successful replay/review. No complete-census claim.
- CRASHED: execution/software/environment/model/tool failure prevents valid
  execution or verification. Administrative unavailability may be REJECTED.
- FROZEN-NEGATIVE is NOT used for operational failure or partial coverage;
  no negative scientific branch is expected for this finite-decision gate.

No new addition-count record, orientation-row completeness, global optimum,
universal no-54, or sigma^0/2 census follows, even if all 608 pairs complete.
Write PLAN, ITERATE, AUDIT, REVIEW, RESULTS and exact staged cs/RESULTS.md and
cs/PROGRESS.md append proposals before campaign freeze/close. Shared ledgers
remain untouched. Record measured full usage externally after close. Retain
mini-0 through all planned/review-driven replay; explicitly hand off to
QECResidual147 and Main only after every monitor attempt is freshly terminal.
