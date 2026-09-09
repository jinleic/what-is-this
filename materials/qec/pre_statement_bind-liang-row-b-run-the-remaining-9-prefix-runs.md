# Pre-statement — bind Liang row b: run the remaining 9 prefix runs

Cycle `20260909T021553Z_fa9cef` · target `math/qec` · agent QecRowB.
Frozen before `campaign.py init`; this file is the hmz freeze artifact and the
campaign preregistration.

## Fixed context (inherited, not re-litigated)

Run `20260908T194313Z_2eb90dad_d3c7ccce6e21` (FROZEN-INCONCLUSIVE) banked:
`liang270b_base_initial.json` COMPLETE (2026-08-25, 1,212.8 s, exhaustive
through 16), `liang270b_base_replay.json` COMPLETE (1,326.0 s, exhaustive
through 16), and `liang270b_prefix55_initial.json` COMPLETE (4,173.2 s,
NO_LOGICAL_THROUGH_CAP exhaustive through weight 18) — all
`validate_run_record` PASS under the pinned prefix solver (binary sha256
`12dd8b65…e61e769`). The fresh root prefix partition is
{55, 56, 135, 175, 269} (complete). Row a is already CERTIFIED_EXACT
[[270,8,20]] and admitted to `EXACT_REFERENCES` (`exp070_odd_exact`), so the
k=8 domination threshold at n=270 is 20 regardless of this run.

Outstanding for row b under the two-run (initial+replay) stop rule: **9
solver runs** — prefix55 replay; prefixes 56, 135, 175, 269 initial+replay
(each ~3,700–4,300 s, single-threaded, nice ≥ 10) — then `assemble` +
`validate`.

## Gate question

Can Liang Table III row b be bound honestly within the ~4 h wall window:
the 9 remaining prefix runs each passing `validate_run_record`, in prefix
pairs completed in partition order (55 first), then
`assemble --target liang270b` writing the certificate and
`validate --target liang270b` re-printing `CERTIFIED_EXACT` — under the
pinned solver, the strict no-`allow_stale_thresholds` contract, and the
pair-stop time-box?

## Bounded objective (from the cycle assignment; hard ceilings)

- **R1 — the 9 prefix runs.** Two lanes, one detached solver run per mini at
  a time (mini-pro + mini-0, both free): `nice -n 15 uv run --no-sync python
  experiments/exp070_n270_connected_cluster.py run --target liang270b --mode
  prefix --prefix <P> --run-id <initial|replay> --solver
  results/partial_runs/exp070_n270_cluster/solver/dist-m4ri-src/src/dist_m4ri_prefix
  --solver-source-commit 538d119f6e98b3782415eb78f7ce322a076f8101 --timeout
  7200`. Pair order: 55 (replay only), then 56, 135, 175, 269
  (initial+replay). **Pair-stop time-box:** the compute window closes ≈ 4 h
  after init (~06:20Z); never start a prefix initial unless its initial AND
  replay both fit in the window at the 4,400 s worst case each; never start
  a replay unless the replay alone fits. Stop cleanly at the cap and bank
  what is complete.
- **R2 — assemble + validate (only if all 12 records exist).**
  `assemble --target liang270b` writes
  `results/certificates/exp070_270_8_20_liang270b_distance.json` and prints
  `CERTIFIED_EXACT`; `validate --target liang270b` re-validates. Then admit
  row b to `EXACT_REFERENCES` in `experiments/exp055_odd_lattice_sweep.py`
  following the row-a pattern (`certificate_kind="exp070_odd_exact"`, local
  (15,9) terms A=[[0,0],[6,1],[6,2]], B=[[0,0],[4,4],[14,8]]) and re-check
  `_validated_exact_reference` + unchanged `domination_threshold(270,8)==20`.

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
  heavy jobs; records synced back into the local run dir before audit.

## Verdict mapping (campaign.py terminal verdicts)

- C1 ∧ C2 complete (row b certified and admitted) → **FROZEN-CERTIFIED**.
- C1 holds for ≥1 completed run but C2 incomplete, with the incomplete
  prefixes named explicitly → **FROZEN-INCONCLUSIVE**.
- Any C3 premise fails, or host obstruction prevents banking any run →
  **FROZEN-NEGATIVE**.

## Scope boundary

In scope: the 9 liang270b prefix runs on the minis; record validation;
assemble/validate iff the 12-record set completes; the conditional
EXACT_REFERENCES admission of row b; environment provisioning of mini-0
(path-identical tree + the identical homebrew `libpng16.16.dylib` artifact,
sha256 `0af10df8…dc33` — the solver binary itself is never modified).
Out of scope: the 21 k=8 residuals (witness 22 vs threshold 20), 147 k=12
undecided, 210 k=20/24 unpromoted survivors; any re-screen or rebind (the
threshold does not move); any change to solver, drivers, or validator
semantics; any other target.

## Resource envelope and compute placement

Heavy compute detached on the minis (`nohup caffeinate -s … > run.log 2>&1
&`), one single-threaded solver run per host at a time, nice 15 (exp070
refuses < 10), per-run `--timeout 7200`. Control plane strictly local.
Remote layout path-identical: `/Users/jinleic/jinleic-workspace/math/qec`.
mini-pro is already provisioned (tree hash-parity + libpng repaired last
cycle); mini-0 is provisioned this run and hash-verified before use.
Heartbeats via `campaign.py heartbeat` during compute.

## Exact verification command

Run locally at the end (re-derives everything from the frozen evidence):

```
python3 campaigns/<run_id>/audit_verify.py
```

which checks: (a) premise hashes (drivers, solver, solver_build) at freeze
time against the values recorded in `sha256s.txt`; (b) `validate_run_record`
on every liang270b run record present (base pair, prefix55 initial, and
every record this run produces); (c) when the certificate exists,
`validate_exact_certificate_payload` plus the exp055
`_validated_exact_reference` admission check and
`domination_threshold(270,8)==20`; (d) the explicit list of prefix pairs NOT
completed. Plus `shasum -a 256 -c sha256s.txt` inside the frozen run dir.
