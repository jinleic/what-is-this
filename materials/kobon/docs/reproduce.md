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
