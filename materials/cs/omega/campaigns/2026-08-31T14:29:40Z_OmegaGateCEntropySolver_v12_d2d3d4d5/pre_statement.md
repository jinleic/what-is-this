# PRE-STATEMENT — OmegaGateCEntropySolver v12

## Scope and non-claim

This campaign certifies only the repaired generalized-entropy/absorbed-defect
objective on the frozen fixed-m affine ray inherited from v11. The exact guard
shows that the ray preserves the base mass, which is approximately
`0.9999999999997855`, not `1`. Therefore even a strict descent is not a feasible
construction, does not refute constrained local optimality, and does not improve
an omega bound.

No v12 execution has occurred at registration time.

## Parent failure and evidence used to design v12

The v11 production aborted before its base checkpoint because it asserted that
the repaired raw interval enclosed `FROZEN_RAW` from the pre-repair dependency
objective. A bounded diagnostic proved the intervals disjoint: repaired raw
`[2.3715518061851393, 2.3715522159185425]`; legacy anchor
`[2.37155383583508, 2.3715538358350807]`.

A second bounded v11 diagnostic found only 53 of 192 base entropy certificates
passed the exact `1/2^32` KKT gate under L-BFGS-B alone. A pre-stated two-stage
experiment then used each L-BFGS-B point to initialize SciPy `least_squares` on
`B (m softmax(B^T y)) - b` with analytic Jacobian
`B (diag(q) - q q^T/m) B^T`. All 192 base certificates passed the unchanged
outward interval checker; maximum certified KKT upper bound was
`1.4456822252392066e-10`. This was base-only diagnostic evidence, not rung or
decision evidence.

## Frozen v12 change

Release source SHA-256 is
`a5b1ba39bb86b0652e7c73a7618d32b1e6b4bee5e13b8e16bc4b8df0bc92dc91`.
It changes only the multiplier solver, repaired-raw invariant, campaign-local
v12 filenames, solver metadata, and loop completion:

1. Run the inherited L-BFGS-B solve.
2. Warm-start `least_squares(method="trf", x_scale="jac")` on the KKT residual,
   using the analytic generalized-softmax Jacobian, tolerances `1e-15`, and
   at most 4,000 function evaluations.
3. Keep solver status as `COMPUTATIONAL-EVIDENCE`; accept a multiplier only
   when the independently recomputed outward Arb KKT upper bound is at most
   exact `1/4294967296`.
4. Recompute `(target-R)/M` through the repaired fields and require exact
   endpoint identity with `base["raw"]`.
5. Record `FROZEN_RAW` only as legacy context for the old objective. It is not
   a v12 invariant and cannot accept or reject the repaired computation.
6. Evaluate every registered divisor `[16, 8, 4, 2, 1]`; retain the first
   accepted strict rung without terminating the loop early.

All four code dependencies and all three numeric inputs are copied into this
campaign and checked against hard-coded SHA-256 digests before loading.

## Decision rule

For each point, all 192 entropy certificates must pass the exact Arb KKT gate.
A rung is strict only when both base and rung pass every required certificate
and the outward upper endpoint of the common-denominator full-objective
difference is strictly below zero. The decision records all five attempted
divisors and the first strict one, if any. If no rung meets both conditions, the
campaign remains OPEN.

## Registered checks and outputs

Before release:

- compile the source in memory without writing bytecode;
- run `--prepare-only` to exercise frozen hashes and 270/270 exact coordinate
  guards without any entropy or objective point evaluation;
- run `--mock-self-test` to exercise wrapper, interval, generalized-mass,
  projection, cache, and counterfactual plants without a real point;
- verify all campaign checksums and the direct campaign-local no-overwrite
  output perimeter.

Reserved static outputs are `prepare_only_v12.json`,
`prepare_only_v12_stdout.log`, `mock_self_test_v12.json`, and
`mock_self_test_v12_stdout.log`. Reserved production outputs are
`results_v12_production.json` and `run_v12_production_stdout.log`. Existing
outputs and `.tmp` paths must be refused, never overwritten.

## Registered release command

```sh
cd cs && set -o noclobber && nice -n 10 env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 ./.venv/bin/python omega/campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/gate_c_entropy_solver.py.asrun_v12 omega/campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/results_v12_production.json > omega/campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/run_v12_production_stdout.log 2>&1
```

The source applies a 900-second production alarm. The shell runs at niceness
10 with all named numerical thread pools pinned to one and suppresses bytecode
writes.

## Static-check result

The frozen source compiled in memory. The initial checksum ledger passed
10/10. Prepare-only passed with 270/270 exact coordinate reconstructions,
189 expected Part certificates plus three expected GlobalStage certificates,
and no entropy or objective point evaluation. The mock self-test passed every
boolean plant with no real point. Pre-existing direct-campaign output and an
off-campaign `/tmp` output were both refused before prepare; the former was
retained and the latter was not created. Production remains held pending the
focused read-only source audit.

## Focused release audit

A fresh read-only Fable audit of the exact release source returned
`PASS — no release blocker found`. It raised two residual uncertainties, both
resolved by owner reads of the frozen dependencies before release:

1. `interval_core.Ival.pos()` is the absolute-value interval operation, not a
   one-sided positive-part clamp. A negative interval maps to `[-hi,-lo]`; an
   interval straddling zero maps to `[0,max(-lo,hi)]`. The KKT maximum therefore
   bounds both residual signs.
2. The candidate-witness dependency defines its historical `MAT` and `CAMP0`
   names at import but reads them only inside its uncalled `main()`. The v12
   `prepare()` path reads `W.MAT` and `W.CAMP0` only after v12 has redirected
   both to hard-hash-verified local payloads.

The source hash is unchanged. Owner release is granted for the exact registered
command and reserved outputs.

## Production result

The exact released command exited `0`. It evaluated the base and every divisor
`[16,8,4,2,1]` in 18.4318 internal seconds. All 1,152 entropy certificates
(192 at each of six points) passed the unchanged exact Arb KKT gate, and all
1,152 least-squares solves reported success. The largest outward KKT upper
bound was `1.4456822252392066e-10 < 1/4294967296`.

Every full-objective difference was strictly positive:

| divisor | outward full difference, moved minus base |
|---:|---|
| 16 | `[6.475289389514407e-9, 6.475289390521749e-9]` |
| 8 | `[1.2954036629144681e-8, 1.2954036630152025e-8]` |
| 4 | `[2.5911541940662355e-8, 2.5911541941669702e-8]` |
| 2 | `[5.182659614230212e-8, 5.182659614330947e-8]` |
| 1 | `[1.0365687869954472e-7, 1.0365687870055208e-7]` |

Thus no registered move descends on this named fixed-m ray. Formal status is
**OPEN**, not local optimality: the points violate the registered `sum=1`
constraint, and five positive ray samples do not exclude other feasible or
infeasible directions. Result SHA-256 is
`0403e652099a217d421c9bb1aa7fd56f07fb6daff72ddf14273d6a6dff063122`;
combined stdout/stderr SHA-256 is
`d365a5644fb06d8fb929c8e1368e9d9cb37a3b9d17c5d02b25147bec4550fb48`.
