# Reproduction commands

All commands below are intended to be run from the workspace root.  They name
real files present in the 2026-08-21 campaign snapshot.

```sh
cd ~/jinleic-workspace
PY="$PWD/scratch/kobon-audit/venv/bin/python"
"$PY" -c 'import pysat; print("pysat", pysat.__version__)'
```

Use `scratch/kobon-audit/venv/bin/python`, not whichever `python` happens to be
on `PATH`.  `python-sat`/`pysat` is required by the engine, CNF builders, the
self-test, and `max_packing.py`.  The two appendix obstruction checkers use only
the standard library.

## A. Maiorana `n=14`, 54-triangle certificate

### Pinned solution 1, fail closed

```sh
"$PY" math/kobon/verification/maiorana_verify.py \
  scratch/kobon/n14/maiorana/sol1_lines_rational.json
printf 'exit=%s\n' "$?"
```

Success is exit `0` and the final JSONL record

```json
{"stage": "CERTIFIED", "K_gen_14_lower_bound": 54, "encoding_convention": "campaign engine K_gen"}
```

The checked source SHA-256 must be
`83fc26666d73cd39791018a72aba71caeac4a3f7efb27109ed2955b4540b2523`.
The archived successful transcript is
`math/kobon/verification/maiorana_verify.log`.

### All 15 pinned arrangements

The sweep script intentionally resolves its `maiorana/` input directory next to
its executed source.  Run the original source-location copy:

```sh
"$PY" scratch/kobon/n14/sweep_all_verify.py
printf 'exit=%s\n' "$?"
```

Success is exit `0`, fifteen rows with `"pass": true`, and

```json
{"stage": "SUMMARY", "solutions_verified": 15, "all_pass": true}
```

The archived transcript is
`math/kobon/verification/sweep_all_verify.log`.  The copy of the script in
`math/kobon/verification/` is byte-identical but retains that source-relative
input contract.

## B. Exact finite checks for the `n=14` and `n=18` obstructions

These do not launch a SAT solver.

```sh
"$PY" math/kobon/verification/appendix_independent_check.py
printf 'n14_exit=%s\n' "$?"

"$PY" math/kobon/verification/appendix_check_n18.py
printf 'n18_exit=%s\n' "$?"
```

The expected final lines are respectively:

```text
CERTIFICATE OK: Theorem M holds for the full Q=3 branch at n=14.
CERTIFICATE OK: n=18 Q=3 legs (Aprime)+(Bprime) hold.
```

Any assertion failure exits nonzero.  The scripts certify the finite
cycle/placement legs; the accompanying hand proofs and hypotheses are in
`scratch/kobon/n14/escape_obstruction.md` and
`scratch/kobon/n18/obstruction_n18.md`.

## C. Rebuild an `n=14` concurrency cube or monolith

The builder reads `N14_TARGET` at import time and writes
`scratch/kobon/n14/<NAME>.cnf`.  Its exact registered names are:

```text
n14monolith
n14c-q0tp1  n14c-q1tp1d  n14c-q1tp1s
n14c-q0tp2d n14c-q0tp2s  n14c-q0k4
```

Inspect the target-dependent signature lattice without creating a CNF:

```sh
N14_TARGET=54 "$PY" scratch/kobon/n14/concurrency_cases.py enum
N14_TARGET=54 "$PY" scratch/kobon/n14/concurrency_cases.py orbits
```

Build one registered cube at target 54:

```sh
N14_TARGET=54 "$PY" scratch/kobon/n14/concurrency_cases.py build n14c-q0tp1
```

Build the assumption-free monolith at any target, for example 55 or 56:

```sh
N14_TARGET=55 "$PY" scratch/kobon/n14/concurrency_cases.py build n14monolith
N14_TARGET=56 "$PY" scratch/kobon/n14/concurrency_cases.py build n14monolith
```

`n14monolith` always writes `scratch/kobon/n14/n14monolith.cnf`; rename or hash
it before building a different target.  A cube command overwrites an existing
same-name CNF, so preserve an authoritative artifact before rebuilding it.  A
successful build prints the variable count, clause count, byte count, elapsed
time, and (for a cube) its fixed parallel classes, planted triples, multipoints,
and per-line capacities.

For comparison, the current target-specific monolith headers are recorded in
`experiments.md`; they are not inferred from filenames.

## D. Re-run exact max packing on any witness JSON

Run the source-location script, whose workspace-root calculation is based on
that location:

```sh
"$PY" scratch/kobon/max_packing.py \
  scratch/kobon/discoveries/n20_known116.json \
  --from 116 --time-budget 900
```

For a Maiorana witness:

```sh
"$PY" scratch/kobon/max_packing.py \
  scratch/kobon/n14/maiorana/sol1_lines_rational.json \
  --from 54 --time-budget 900
```

Generic form:

```sh
"$PY" scratch/kobon/max_packing.py FILE.json \
  --from KNOWN_COUNT --time-budget SECONDS
```

The program writes `FILE.maxpacking.json` next to the input and will overwrite
an existing result of that name.  Do not treat exit `0` as an optimality proof:
require both the terminal `RESULT` record to contain
`"proved_optimal": true` and the written file to contain
`"optimality_proved": true`.  On a time-budget stop it exits normally but
writes a non-proved result.

The archival copy can be invoked if the engine path is supplied explicitly:

```sh
PYTHONPATH="$PWD/math/kobon" "$PY" \
  math/kobon/verification/max_packing.py FILE.json \
  --from KNOWN_COUNT --time-budget SECONDS
```

## E. Engine self-test

```sh
"$PY" math/kobon/test_engine.py
printf 'exit=%s\n' "$?"
```

Success is exit `0` and the standard `unittest` summary ending in `OK`.  The
suite exercises exact geometry, bounded-face and per-line `TC` constraints,
multipoint charging, cardinalities, and SAT assumptions; it therefore requires
`pysat` and its bundled/available solver backend.

## F. Proof-producing UNSAT promotion (when a live run settles)

A Kissat log alone is discovery evidence.  For a proof-producing run, preserve
the exact CNF and complete DRAT, then run the campaign's independent checker:

```sh
scratch/kobon-audit/tools/drat-trim/drat-trim INSTANCE.cnf PROOF.drat -w
printf 'drat_trim_exit=%s\n' "$?"
```

Only `s VERIFIED` for that exact pair promotes the branch.  Never promote an
in-flight/truncated DRAT or a solver-only `s UNSATISFIABLE` line.

## G. Exact face square-penalty counterexample

This replay uses exact rational geometry and no SAT solver:

```sh
"$PY" scratch/kobon/square_penalty_counterexample.py \
  --out /tmp/square_penalty_counterexample.json
printf 'square_exit=%s\n' "$?"
```

Success is exit `0` and a JSON object containing:

```json
{
  "status": "VERIFIED_COUNTEREXAMPLE",
  "triangular_faces": 30,
  "claimed_lhs_3F": 90,
  "claimed_rhs": 89,
  "sector_defect_sum": 27,
  "shared_edges_with_two_multipoint_endpoints": 26,
  "sector_excess_identity_rhs": 1
}
```

The release-local copy is also runnable from
`math/kobon/release/kobon-2026-08/`:

```sh
"$PY" math/kobon/release/kobon-2026-08/scripts/square_penalty_counterexample.py
```

## H. Coefficient-2/3 endpoint certificate at `n=7`

Regenerate the essential all-degeneracy direct-gap CNF:

```sh
"$PY" scratch/kobon/square_penalty_sat.py 7 \
  --coefficient 2/3 \
  --out /tmp/c23_n7_violation.cnf
sha256sum /tmp/c23_n7_violation.cnf
```

The expected CNF SHA-256 is
`075ca17768ad70e4ff69673543c85d629a69737656b36211c949b8e0a32676d9`.
Check the bundled proof:

```sh
scratch/kobon-audit/tools/drat-trim/drat-trim \
  math/kobon/release/kobon-2026-08/certificates/c23_n7_violation.cnf \
  math/kobon/release/kobon-2026-08/certificates/c23_n7_violation.drat -w
```

Success is `s VERIFIED`. The hash-bound metadata and checker-core counts are
in `math/kobon/release/kobon-2026-08/certificates/c23_n7_certificate.json`.

## I. Hereditary, vertex-sector, shared-ray, and Pappus probes

Replay the exact-rational sector and shared-ray implementation audit:

```sh
"$PY" scratch/kobon/sector_bound_audit.py \
  --samples-per-mode 100 \
  --out /tmp/sector_bound_audit.json
```

Success is exit `0`, `"status": "PASS"`, 3,600 arrangements, 93,000 line-pair
checks, 83,115 face-sector assignments, 17,325 shared-pair face pairs, 1,591
same-ray cases, a four-sector sharpness control, and zero violations.

Replay the guarded dual-Pappus relaxation-gap probe:

```sh
"$PY" scratch/kobon/pappus_relaxation_probe.py \
  --out /tmp/pappus_relaxation_probe.json
```

Success is exit `0`, the base abstract model SAT with a false Pappus
conclusion, and the same assumptions UNSAT after the three projective
conclusion clauses are added.

Regenerate the live combined \(n=12,T=39\) discovery instance:

```sh
"$PY" scratch/kobon/gap_faces.py 12 39 \
  --sector-bounds \
  --sub-bound 11:32 \
  --sub-bound 10:25 \
  --out /tmp/n12_gap_sector_deletion_t39.cnf
```

The generator must report 374,381 variables / 830,030 clauses and reproduce
SHA-256 `5435cb9958710878b1a68bc112e353a46d9a031a723e1afaa31191a9a9546909`.
The cited smaller-\(n\) inputs and promotion rules are pinned in
`scratch/kobon/n12_gap_sector_deletion_t39.metadata.json`. This CNF is a sound
discovery lane, not a certificate. An UNSAT solver result requires a complete
DRAT and the promotion gate in section F.

Regenerate the opt-in shared-ray solver-diversity lane:

```sh
"$PY" scratch/kobon/gap_faces.py 12 39 \
  --sector-bounds \
  --shared-ray-bounds \
  --sub-bound 11:32 \
  --sub-bound 10:25 \
  --out /tmp/n12_gap_shared_ray_deletion_t39.cnf
```

The generator must report 374,381 variables / 841,910 clauses and reproduce
SHA-256 `0aa9e81e1806d5ff7413669d7357454c9e08a496d2172d8b4a512b27115b7e7d`.
`scratch/kobon/shared_ray_bound_experiment.json` records the proof, exact audit,
mixed A/B benchmarks and promotion gate. This lane is discovery-only.
 
## J. Independent faces-only certificate at `n=10`

Regenerate the target-26 CNF and compare it byte-for-byte:

```sh
"$PY" scratch/kobon/faces_only.py 10 26 \
  --out /tmp/ladder_n10_t26_faces.cnf
cmp scratch/kobon/ladder_n10_t26_faces.cnf \
  /tmp/ladder_n10_t26_faces.cnf
```

Then check the complete proof:

```sh
scratch/kobon-audit/tools/drat-trim/drat-trim \
  scratch/kobon/ladder_n10_t26_faces.cnf \
  scratch/kobon/ladder_n10_t26_faces.drat -w -t 200000
```

Success is `s VERIFIED`. The expected checker statistics and SHA-256 pins are
in `scratch/kobon/ladder_n10_t26_faces_certificate.json`; the saved transcript
is `scratch/kobon/ladder_n10_t26_faces.dratcheck.log`.

## K. Endpoint-closure and multipoint-incidence frontier

Replay the expanded exact audits:

```sh
"$PY" scratch/kobon/sector_bound_audit.py \
  --samples-per-mode 100 \
  --out /tmp/sector_bound_audit.json
"$PY" scratch/kobon/reified_face_audit.py \
  --out /tmp/reified_face_audit.json
"$PY" scratch/kobon/chirotope_gp_audit.py \
  --samples-per-mode 100 \
  --out /tmp/chirotope_gp_audit.json
```

All three print `"status": "PASS"`. The endpoint audit reports 1,607
colliding endpoints, 84 colliding same-ray cases, 42 single-line shared
segments, 1,591 shared-pair segments, and zero violations. The fixed-witness
audit accepts exactly 25 and 38 faces; the \(n=12\) witness has eight
multipoint-incident faces. The chirotope audit checks 195,000 determinant
signs and 1,386,000 Grassmann--Plücker relations with zero mismatch, while
reporting `NO_GAP_FOUND_THROUGH_N7_PROPAGATION_ONLY`.

Regenerate the primary endpoint lane:

```sh
"$PY" scratch/kobon/gap_faces.py 12 39 \
  --endpoint-closure \
  --k4-bound \
  --simple-bound 37 \
  --sub-bound 11:32 \
  --sub-bound 10:25 \
  --out /tmp/n12_gap_endpoint_mi_k4_deletion_t39.cnf
```

Expected header: 375,037 variables / 1,005,366 clauses. Expected SHA-256:
`6dca104aa120ee1f209aada8aca07a08bf3341f1335e6e99f4ae1f11bc2ba76b`.

Regenerate the high-risk \(K(9)=21\) hereditary variant by adding
`--sub-bound 9:21`. Expected header: 842,317 variables / 1,917,926 clauses.
Expected SHA-256:
`34406004c932ce77d33f76f8b7d19c7b0f5a2e3e241e94651846269bf750f2df`.

`scratch/kobon/frontier_endpoint_research.json` binds the proofs, source
survey, exact outputs, A/B controls, hashes, and live-process names. Both CNFs
are discovery inputs only. SAT requires exact rational realization; UNSAT
requires a proof-producing rerun and independent proof checking.
