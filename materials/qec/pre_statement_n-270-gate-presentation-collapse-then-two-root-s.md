# Preregistration — n=270 gate: presentation collapse then two-root screening

Cycle: `20260908T174625Z_a9a475` · Target: `math/qec` · Agent: QecRunner
Frozen: 2026-09-08T17:53Z (before `campaign.py init`; this file is the hmz
freeze artifact and the campaign prereg).

## Gate question

Can the n=270 odd-lattice BB frontier be (a) collapsed from the two isomorphic
noncyclic presentations (15,9) and (45,3) to a single screened presentation,
and (b) screened — 5,024 noncyclic (15,9) classes plus 688 cyclic (27,5)
classes — under the EXP-055 v3 screen protocol with the EXP-069 resolve layer,
and (c) can Liang et al.'s two published Table III [[270,8,20]] twisted-torus
rows be bound as local exact references by the pinned EXP-070 prefix-partitioned
connected-cluster protocol ("two-root" style rooted exhaustion)?

## Disclosed starting state (no pre-staged claim is assumed; all is re-audited)

The following artifacts exist from un-campaigned work of 2026-08-24/25
(untracked scripts `experiments/exp069_n270_frontier.py`,
`experiments/exp070_n270_connected_cluster.py`; checkpoints under
`results/partial_runs/`). This run treats them as *claims to re-verify*, not
as accepted evidence:

- EXP-069 initial screens: COMPLETE for (15,9) [5,024 classes, 282.4 s] and
  (27,5) [688 classes, 26.3 s], protocol hash `0ff81aa98d0d`.
- EXP-069 resolved_v2 screens: COMPLETE for (15,9) [2,031 residuals, 1,215.8 s]
  and (27,5) [234 residuals, 166.9 s]; final shards
  `results/partial_runs/exp055_screen/{15x9,27x5}.json` written.
- EXP-070 pinned dist-m4ri runs (single-threaded, nice ≥ 10, hash-bound prefix
  build): liang270a base (dmin=1..wmax=16) initial+replay COMPLETE; liang270a
  prefix branches 55/62/135/150/159 at weight 18 initial+replay COMPLETE;
  liang270b base initial COMPLETE; liang270b base replay stale-RUNNING
  (orphaned checkpoint from 2026-08-25); liang270b prefix runs ABSENT;
  liang270a monolithic `original_initial` run UNDECIDED_BUDGET (7,260 s;
  superseded by the prefix partition; not required by `assemble`).
- Final shard verdicts as pre-staged: 15x9 — 2,993 dominated_by_witness,
  187 survivor, 1,844 undecided; 27x5 — 454 dominated_by_witness,
  23 survivor, 211 undecided.

## New work this run

1. Fresh EXP-069 transport audit (the collapse): exhaustive bijection,
   inverse, and homomorphism checks of the CRT isomorphism
   Z_15×Z_9 → Z_45×Z_3 over all 135² pairs; candidate-class counts and
   represented-pair totals on both presentations; exact HX/HZ matrix transport
   for both Liang-bound rectangular presentations; published quotient-lattice
   transport for both Table III rows.
2. Fresh validator-v11 (aggregate + record level) validation of both final
   screen shards.
3. EXP-070 `assemble` + `validate` for liang270a: solver-free certificate
   assembly (parity law, BB duality, pure-second-block exhaustive enumeration,
   root prefix partition, weight-20 upper witness, 12 pinned run bindings) and
   independent re-validation. If any step raises, the obstruction is recorded
   as-is.
4. Reference promotion + monotone rebind: register the certified liang270a
   [[270,8,20]] row in EXP-055 (`EXACT_REFERENCES` + `FRONTIER_REFERENCE_KINDS`
   dispatch entry, following the exp067 pattern), then run the EXP-063
   monotone rebind on the two n=270 shards only. Thresholds only rise; every
   carried witness is rechecked; originals are archived by the tool.
5. Honest residual accounting (see verdict mapping).

## Falsifiable acceptance criteria

- **C1 (collapse certified):** the fresh transport audit reports
  `all_coordinates_bijective`, `inverse_exact`, `homomorphism_exhaustive`,
  candidate classes 5,024 on (15,9) and 5,024 on (45,3), equal represented
  pair totals, and `HX_exact`/`HZ_exact` matrix transport for both Liang rows;
  both published transports report quotient order 135, surjective quotient map,
  vanishing lattice relations, and k=8 in both rectangular presentations.
- **C2 (screen evidence valid):** both final shards pass
  `_validate_screen_shard_aggregates` and `_validate_screen_shard_records`
  (validator v11) re-run fresh this run; candidate totals 5,024/688; every
  record has a terminal verdict.
- **C3 (Liang binding):** EXP-070 certificate for liang270a assembles and
  validates `CERTIFIED_EXACT` with d=d_X=d_Z=20. (Both rows = full C3; one row
  = partial C3. liang270b is known-infeasible this cycle: its missing base
  replay + 5 prefix branches × initial/replay are ≈ 11 pinned solver runs of
  1,200–4,300 s each, ≫ the 7,200 s cycle budget.)
- **C4 (closure standard):** after the rebind, zero survivor/undecided rows
  across the two n=270 shards. Anticipated NOT met: residual k=8 classes with
  verified witnesses only at weight ≥ 22, 147 k=12 classes whose exact
  fallback exceeded the bounded screen budget, and 210 survivor classes
  (k=20/24, decided d > threshold) awaiting exact certification.

## Verdict mapping (outcome → campaign.py terminal verdict)

- C1 ∧ C2 ∧ C3(both rows) ∧ C4 → **FROZEN-CERTIFIED** (n=270 closure achieved;
  `current_gate` advances).
- C1 ∧ C2, with C3 partial/obstructed and/or C4 unmet, all residuals named with
  evidence → **FROZEN-INCONCLUSIVE** (gate advanced: collapse + screen
  certified; closure not achieved; no overclaim).
- C1 fails (the (15,9)/(45,3) collapse is refuted) OR C2 fails (pre-staged
  screen evidence does not survive fresh validation) OR a COMPLETE liang270a
  run record fails validation → **FROZEN-NEGATIVE** (a checked premise of the
  documented gate is false).
- Harness/tooling failure before any check completes → **CRASHED**.

## Scope boundary

No new solver runs (no dist-m4ri, no CP-SAT, no CDCL launches) beyond replay
of the already-COMPLETE pinned run records; no liang270b compute; no survivor
promotion; no k=12 exact fallback; no deep-reduction re-screen of the residual
k=8 classes; no edits outside `math/qec/`, this cycle's usage record, and the
domain ledgers (`math/PROGRESS.md`, `math/qec/README.md`); the 22 shards of the
n≤234 closure are not rebound. Solver-free re-validation only; every reused
pre-staged artifact is disclosed above and re-audited fresh.

## Exact verification command

On the synced copy on mini-pro.local (`~/qec-work/math/qec`), detached:

```
nohup caffeinate -s uv run --no-sync python campaigns/<run_id>/audit_verify.py \
  > campaigns/<run_id>/audit_verify.stdout.log 2>&1 &
```

`audit_verify.py` re-runs, fresh: (1) the EXP-069 transport audit and every C1
assertion; (2) validator-v11 aggregate+record validation of both rebound
shards; (3) EXP-070 `validate_exact_certificate_payload` on the liang270a
certificate (which itself re-validates all 12 pinned run records, the hash-bound
solver build, the parity law, the duality, the pure-second enumeration, the
prefix partition, and the weight-20 witness); (4) rebind-archive consistency.
It writes a JSON summary with per-artifact SHA-256; the summary and stdout log
are copied back into the sealed run dir.

## Resource envelope

Cycle budget: 1 run, 7,200 s wall (qec target cap). Heavy compute
(transport audit enumeration, certificate assembly/validation, shard
validation, rebind witness recheck) runs on mini-pro.local detached via
`nohup caffeinate -s`; control plane (campaign.py, state.json, ledgers) stays
local. No claim is made beyond what the frozen evidence supports; timeout or
obstruction is recorded, never converted into a bound.
