# Pre-statement — n=270 closure: deep target-20 re-screen of the 1,657 k=8 residuals + Liang row b binding

Cycle `20260908T193320Z_ba7b32` · target `math/qec` · agent QecRunner2.
Frozen before `campaign.py init`; this file is the hmz freeze artifact and the
campaign preregistration.

## Fixed context (inherited, not re-litigated)

Run `20260908T175247Z_e40ae2c8_8d516c12d880` (FROZEN-INCONCLUSIVE) certified
fresh: the (15,9)/(45,3) collapse, the validated n=270 screen (strict
validator-v11), and Liang Table III **row a** as CERTIFIED_EXACT
[[270,8,20]] admitted to `EXACT_REFERENCES` (battery fingerprint
`4dd928285865…`, k=8 domination threshold at n=270 = 20). Its named residuals
are this run's gate:

- **1,657 k=8 undecided** classes (1,474 in 15x9 + 183 in 27x5) in the rebound
  shards `results/partial_runs/exp055_screen/{15x9,27x5}.json` (sha256
  `0f2d0670…4713e` / `53abd152…f6cff`), carried witnesses only at weights
  22–36 vs. threshold 20. Record thresholds in the shards already read 20.
- **Liang row b unbound**: `liang270b_base_initial.json` COMPLETE (1,212.8 s,
  exhaustive through 16), `liang270b_base_replay.json` stale-RUNNING
  (2026-08-25 orphan), all five prefix pairs absent.

## Gate question

Can the two named residuals be advanced honestly within the bounded envelope —
(a) the deferred deep target-20 re-screen of the 1,657 k=8 classes, and
(b) the Liang row b two-root binding — without breaking the strict
validator-v11 publication contract?

## Bounded sub-objectives (from the cycle assignment; hard ceilings)

- **S1 (≤ 2 h) — deep target-20 re-screen.** For every k=8 undecided record in
  the rebound shards, run the EXP-069 deepen route
  (`resolve_candidate_record` k=8 branch → `deepen_record(tries=2000)`):
  deterministic information-set witness reduction whose seed incorporates the
  record identity **and the live threshold 20**, so this is a fresh witness
  search, not a replay of the threshold-18 trajectory. Mechanically:
  `experiments/exp069_n270_frontier.py resolve --lattice 15x9 --tries 2000`
  and `--lattice 27x5 --tries 2000` (the live battery forces a fresh
  protocol id → fresh initial screen → fresh resolve), which writes and
  strictly validates the published shards, then
  `experiments/exp063_reference_rebind.py run --lattices 15x9,27x5`
  (monotone; archives originals; strict v11 on the rebound shards).
- **S2 (≤ 2 h) — Liang row b binding.** Complete the stale base replay
  (`exp070 run --target liang270b --mode base --run-id replay --force`), then
  the five prefix pairs in partition order. **Stop after the first 2 full
  solver runs** (worst case ≈ 1,212 s + 4,303 s): if the bind is not
  complete by then, record the partial state and do not start a third run.

## Falsifiable acceptance criteria

- **C1 (S1 done):** both published shards re-validated: aggregate checks plus
  per-record checks **strict** (no `allow_stale_thresholds`); every record
  verdict transition is `undecided → dominated_by_deep_reduction` (or stays
  `undecided`); every newly dominated k=8 record carries a physically
  verified witness of weight ≤ 20 (HX-kernel, HZ-rank increment); the exact
  resolved/remaining undecided counts are reported.
- **C2 (S2 progress):** every completed `liang270b` run record passes
  `validate_run_record` (COMPLETE, returncode 0, hash-bound solver,
  NO_LOGICAL_THROUGH_CAP at its cap). Full bind additionally requires:
  base initial+replay + 5 prefix initial+replay, then
  `assemble --target liang270b` and `validate --target liang270b` printing
  `CERTIFIED_EXACT`.
- **C3 (premises):** the rebound shards validate strict under the live
  battery before any work; the prefix solver binary on the compute host
  matches the pinned sha256 `12dd8b65…e61e769`; local↔remote hashes of the
  four experiment drivers match.

## Verdict mapping (campaign.py terminal verdicts)

- C1 ∧ full bind of row b (C2 complete incl. CERTIFIED_EXACT) →
  **FROZEN-CERTIFIED**.
- C1 holds, or C2 partial (some run records validated), with named residuals
  remaining → **FROZEN-INCONCLUSIVE**.
- Any C3 premise fails → **FROZEN-NEGATIVE**.

## Scope boundary

In scope: the 1,657 k=8 re-screen; the liang270b base replay + prefix runs
within the 2-run stop rule; strict validation and publication of the shards.
Out of scope (next cycles): the 147 k=12 exact fallback; exact certification
of the 210 k=20/24 survivors; registering row b into `EXACT_REFERENCES`
unless C2 completes (if it does, the registration edit follows the row-a
pattern and is re-validated); any change to the solver, drivers, or
validator semantics. The k=12/k=20/24 records ride the same resolve pipeline
but are not acceptance targets.

## Resource envelope and compute placement

All compute on **mini-pro.local only** (shared with a Lean build: solver jobs
single-threaded, nice ≥ 10, and the screen pipeline capped at 6 threads;
checked via `pgrep -fl "lean|lake"`). Control plane strictly local. Remote
jobs detached: `nohup caffeinate -s <cmd> > run.log 2>&1 &`. Remote layout is
the path-identical `/Users/jinleic/jinleic-workspace/math/qec` on mini-pro
(hash-verified against local before work). Exp-070 runs refuse to launch
below nice 10 by construction.

## Exact verification command

Run locally at the end (re-derives everything from the frozen evidence):

```
python3 campaigns/<run_id>/audit_verify.py
```

which checks: (a) local↔remote hash parity of the four drivers and both
shards at freeze time (values pinned in `sha256s.txt`); (b) the published
shards validate strict (aggregate + records) under the live battery;
(c) every k=8 record newly dominated by the deep re-screen has a physically
verified weight-≤20 witness; (d) every completed liang270b run record passes
`validate_run_record`; (e) when present, the row-b certificate validates.
