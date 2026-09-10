# Pre-statement — finish binding Liang row b: run prefix175 initial

Cycle `20260909T064241Z_f8f3a9` · target `math/qec` · actor GPT-6 Astra (Main).
Frozen before `campaign.py init`; this file is the hmz freeze artifact and the
campaign preregistration.

## Fixed context (inherited, not re-litigated)

Run `20260908T194313Z_2eb90dad_d3c7ccce6e21` banked the base pair
(`liang270b_base_initial` 1,212.8 s and `liang270b_base_replay` 1,326.0 s,
both exhaustive through 16) plus `liang270b_prefix55_initial` (4,173.2 s,
exhaustive through 18). Run `20260909T022016Z_81295935_492d78f41de6`
(FROZEN-INCONCLUSIVE) banked prefix55 replay (4,059.9 s), prefix56
initial+replay (4,188.2 s / 4,176.5 s) and prefix135 initial+replay
(4,147.1 s / 4,128.9 s), every record `validate_run_record` PASS fresh
locally under the pinned prefix solver (binary sha256 `12dd8b65…e61e769`).
**8 of the 12 two-root records are banked and valid.** The fresh prefix
partition is {55, 56, 135, 175, 269} (complete). Row a is already
CERTIFIED_EXACT [[270,8,20]] and admitted to `EXACT_REFERENCES`
(`exp070_odd_exact`), so the k=8 domination threshold at n=270 is 20
regardless of this run.

Outstanding for row b: exactly **4 solver runs** — prefix175 initial+replay
and prefix269 initial+replay (observed 4,060–4,188 s each, worst case
4,400 s, single-threaded, nice ≥ 10) — then `assemble` + `validate`, and on
CERTIFIED_EXACT the EXACT_REFERENCES admission of row b.

## Gate question

Can Liang Table III row b be bound completely and honestly this window: the
4 remaining prefix runs each passing `validate_run_record`, then
`assemble --target liang270b` writing the certificate and printing
`CERTIFIED_EXACT`, `validate --target liang270b` re-validating, and the row-b
admission to `EXACT_REFERENCES` passing `_validated_exact_reference` with
`domination_threshold(270,8)==20` unchanged — under the pinned solver, the
strict no-`allow_stale_thresholds` contract, and the pair-stop time-box?

## Bounded objective (from the cycle assignment; hard ceilings)

- **R1 — the 4 prefix runs.** Two lanes, one monitored solver run per mini
  at a time (mini-pro + mini-0, both verified free): `nice -n 15 uv run
  --no-sync python experiments/exp070_n270_connected_cluster.py run --target
  liang270b --mode prefix --prefix <P> --run-id <initial|replay> --solver
  results/partial_runs/exp070_n270_cluster/solver/dist-m4ri-src/src/dist_m4ri_prefix
  --solver-source-commit 538d119f6e98b3782415eb78f7ce322a076f8101 --timeout
  7200`. Lane A (mini-pro): prefix175 initial → prefix175 replay. Lane B
  (mini-0): prefix269 initial → prefix269 replay. Every launch goes through
  the local bounded supervisor `scripts/institute_monitor.py run` (task id
  `liang270b-distance`, run-bound, fixed cycle compute deadline
  `2026-09-09T11:15:00Z`, stall bound 7,200 s, real progress outputs =
  the run record and solver stdout); no raw nohup. **Pair-stop time-box:**
  never start a prefix initial unless initial AND replay both fit before
  that deadline at 4,400 s each; never start a replay unless it fits.
  All attempts retain the original objective start `2026-09-09T02:20:16Z`
  and its 72-hour ceiling; this allocation does not reset that window.
- **R2 — assemble + validate (only if all 12 records exist).** Locally:
  `assemble --target liang270b` writes
  `results/certificates/exp070_270_8_20_liang270b_distance.json` and prints
  `CERTIFIED_EXACT`; `validate --target liang270b` re-validates. Then admit
  row b to `EXACT_REFERENCES` in `experiments/exp055_odd_lattice_sweep.py`
  following the row-a pattern (`certificate_kind="exp070_odd_exact"`, local
  (15,9) terms A=[[0,0],[6,1],[6,2]], B=[[0,0],[4,4],[14,8]]) and re-check
  `_validated_exact_reference` + unchanged `domination_threshold(270,8)==20`.
  If the 12-record set is incomplete in-window: skip, and name the missing
  prefix pairs in AUDIT/RESULTS.

## Falsifiable acceptance criteria

- **C1 (runs):** every completed `liang270b` run record of this run passes
  `validate_run_record` (COMPLETE, returncode 0, hash-bound solver build,
  NO_LOGICAL_THROUGH_CAP exhaustive through 18, invocation self-consistent,
  parent nice ≥ 10, stdout/stderr hashes match).
- **C2 (full bind):** all 12 `liang270b` records present and valid;
  `assemble` prints `CERTIFIED_EXACT` (d = d_X = d_Z = 20) and the written
  certificate passes `validate`; the EXACT_REFERENCES admission edit passes
  `_validated_exact_reference` and leaves `domination_threshold(270,8)==20`.
- **C3 (premises):** before work, local↔remote sha256 parity of the four
  drivers (exp055/exp056/exp067/exp069/exp070 chain as used), the solver
  binary, and `solver_build.json` on every host used; the solver binary
  matches the pinned sha256 `12dd8b65…e61e769`; both minis free of other
  heavy user jobs (transient single-core macOS system daemons such as
  mediaanalysisd on a 10-core host are recorded, not treated as
  obstruction); records synced back into the local run dir area before
  audit.

## Verdict mapping (campaign.py terminal verdicts)

- C1 ∧ C2 complete (row b certified and admitted) → **FROZEN-CERTIFIED**.
- C1 holds for ≥1 completed run but C2 incomplete, with the incomplete
  prefixes named explicitly → **FROZEN-INCONCLUSIVE**.
- Scientifically valid counterevidence refuting the frozen exact-distance
  gate → **FROZEN-NEGATIVE**.
- Operational premise failure, unavailable tools/hosts, or execution failure
  preventing valid evidence → **CRASHED**, or administrative **REJECTED**
  before science. These are never scientific negatives.

## Scope boundary

In scope: the 4 liang270b prefix runs on the minis; record validation;
assemble/validate iff the 12-record set completes; the conditional
EXACT_REFERENCES admission of row b. Out of scope: the 21 k=8 residuals
(witness 22 vs threshold 20), 147 k=12 undecided, 210 k=20/24 unpromoted
survivors; any re-screen or rebind (the threshold does not move); any change
to solver, drivers, or validator semantics; any other target; any further
provisioning (both minis were provisioned and provenance-verified in prior
cycles and are re-verified by hash parity here).

## Resource envelope and compute placement

Heavy compute on the minis under `institute_monitor.py run` (detached
supervisor + caffeinate + mini-local watchdog), one single-threaded solver
run per host at a time, nice 15 (exp070 refuses < 10), per-run
`--timeout 7200`, fixed monitor deadline `2026-09-09T11:15:00Z` and stall
bound 7,200 s. Audit/review/close by `2026-09-09T12:15:00Z`; projected total
owner wall time under 20,000 seconds, below remaining 21,626 seconds.
The stable objective `liang270b-distance` began `2026-09-09T02:20:16Z`;
its original 72-hour end remains `2026-09-12T02:20:16Z`, with termination margin.
Control-plane commands stay local. Remote layout path-identical:
`/Users/jinleic/jinleic-workspace/math/qec`. Heartbeats via
`campaign.py heartbeat` during compute.

## Exact verification command

Run locally at the end (re-derives everything from the frozen evidence):

```
python3 campaigns/<run_id>/audit_verify.py
```

which checks: (a) premise hashes (drivers, solver, solver_build) at freeze
time against the values recorded in `sha256s.txt`; (b) `validate_run_record`
on every liang270b run record present (the 8 banked records plus every
record this run produces); (c) the completeness accounting against the
12-record two-root protocol with any missing pair listed explicitly; (d)
when the certificate exists, `validate_exact_certificate_payload` plus the
exp055 `_validated_exact_reference` admission check and
`domination_threshold(270,8)==20`. Plus `shasum -a 256 -c sha256s.txt`
inside the frozen run dir.

## Independent review and evidence retention

After primary evidence exists, request the sanctioned read-only Opus reviewer
before choosing the verdict. No actor self-review substitutes for that review.
Copy completed records, stdout/stderr, input/build/source provenance and any
certificate into this campaign before freeze; retain existing canonical files
without rewriting banked evidence. The fresh verifier must check full
completeness and certificate/admission, not merely process exit.
