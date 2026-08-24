# Verification scripts

This directory is a byte-for-byte copy of the small, inspectable verifiers used
by the campaign.  The authoritative inputs remain under `scratch/kobon/`; the
commands below are also collected in `../docs/reproduce.md`.  Solver logs and
DIMACS files are **not** proof certificates unless a complete DRAT has been
checked independently.

| File | What it checks | Success / failure contract | Expected terminal record |
|---|---|---|---|
| `maiorana_verify.py` | Fail-closed exact-rational certification of Maiorana solution 1: pinned SHA-256, 14 distinct lines, declared multipoints, 54-triangle census, cevian-aware empty-interior test, and `engine.verify_selection`. | Exit `0` only after all stages; explicit check failures exit `1` with a JSON `FAIL` record. I/O/import exceptions also exit nonzero. | Final JSON line: `{"stage": "CERTIFIED", "K_gen_14_lower_bound": 54, "encoding_convention": "campaign engine K_gen"}`. |
| `sweep_all_verify.py` | Recomputes the exact-rational 54-triangle census for all 15 pinned Maiorana arrangements and checks every file against `SHA256SUMS.txt`. | Exit `0` iff all 15 rows pass; exit `1` if any row fails. Unhandled input errors exit nonzero. | Fifteen solution rows followed by `{"stage": "SUMMARY", "solutions_verified": 15, "all_pass": true}`. |
| `appendix_independent_check.py` | Pure-stdlib exact-integer check of both `n=14`, `Q=3` direction-cycle obstruction legs, with feasible controls. | Normal completion exits `0`; any failed invariant raises `AssertionError` and exits nonzero. | Final line: `CERTIFICATE OK: Theorem M holds for the full Q=3 branch at n=14.` |
| `appendix_check_n18.py` | Pure-stdlib exact-integer check of the `n=18`, `Q=3` paired and triad direction-cycle legs, orbit counts, reflection counts, and feasible controls. | Normal completion exits `0`; any failed invariant raises `AssertionError` and exits nonzero. | Final line: `CERTIFICATE OK: n=18 Q=3 legs (Aprime)+(Bprime) hold.` |
| `max_packing.py` | Builds the exact overlap graph of all nondegenerate line triples and uses SAT to find a maximum independent set; SAT at `K` supplies a packing and UNSAT at `K+1` proves optimality for that fixed arrangement. | A normal run exits `0` even if its time budget expires. **Exit code alone is not a certificate:** require the final `RESULT` row to have `"proved_optimal": true` and the written JSON to have `"optimality_proved": true`. Parse/import/solver exceptions exit nonzero. | Proved runs end with a `RESULT` row and a `written` row; immediately before them, the decisive row is `{"stage":"UNSAT","K": optimum+1,"optimum": optimum}` (spacing may differ). |
| `../release/kobon-2026-08/scripts/square_penalty_counterexample.py` | Exact rational construction of a generic chart of simplicial `A(12,1)`; cross-checks all triangular faces with strict-straddle and direct-gap predicates, then verifies the complete vertex, edge, sector-run, and square-penalty census. | Exit `0` only after every identity and `engine.verify_selection` pass; any mismatch raises `AssertionError`. | JSON object with `"status": "VERIFIED_COUNTEREXAMPLE"`, `"claimed_lhs_3F": 90`, `"claimed_rhs": 89`, and `"sector_excess_identity_rhs": 1`. Persisted output: `../release/kobon-2026-08/verification/square_penalty_counterexample.json`. |
| `../release/kobon-2026-08/scripts/square_penalty_sat.py` + `../release/kobon-2026-08/certificates/c23_n7_*` | Encodes a strict violation of the coefficient-`2/3` face penalty at `n=7`, with zero selected faces allowed, essentiality enforced, and all parallel/concurrency patterns free; the bundled DRAT proves UNSAT. | Regeneration must match the manifest-bound CNF; `drat-trim CNF DRAT -w` must return `s VERIFIED`. | Included transcript: 74,507 variables, 152,419 clauses, 59,128 formula clauses / 211,032 lemmas in core, `s VERIFIED`. |
| `../release/kobon-2026-08/verification/sector_bound_audit.py` | Exact-rational audit of the four/two/zero selected-face sector bounds over deterministic generic, parallel, concurrent and mixed arrangements. | Exit `0` only with zero violations; any violation raises `AssertionError`. | `"status": "PASS"` with 3,600 arrangements, 93,000 line-pair checks and 6,204 multipoint-pair checks. |
| `../release/kobon-2026-08/verification/pappus_relaxation_probe.py` | Shows that the abstract order model admits one guarded non-Pappus pattern, then adds the three projective conclusion clauses and resolves the same assumptions. | Exit `0` iff the base assumptions are SAT and the cut version is UNSAT. | `"status": "PASS"`, both Boolean verdict fields `true`. |
| `maiorana_verify.log` | Captured successful transcript for the pinned solution-1 verifier. | Transcript only; it has no exit code. | Seven JSONL records ending in `CERTIFIED`. |
| `sweep_all_verify.log` | Captured successful transcript for the 15-solution sweep. | Transcript only; it has no exit code. | Sixteen JSONL records ending in `SUMMARY`, `all_pass: true`. |

## Relocation note

The copies intentionally preserve the executed source bytes.  Two scripts use
source-relative paths:

- `sweep_all_verify.py` expects a sibling directory named `maiorana/`; the
  executed source copy lives at `scratch/kobon/n14/sweep_all_verify.py`, beside
  `scratch/kobon/n14/maiorana/`.
- `max_packing.py` derives the workspace root from its source location; run the
  executed source copy at `scratch/kobon/max_packing.py`, or set
  `PYTHONPATH="$PWD/math/kobon"` when invoking this archival copy.

`maiorana_verify.py` takes its witness path explicitly. The two appendix
checkers are pure-stdlib and location independent. The square-penalty replay
and endpoint generator import the bundled release `scripts/engine.py` and
require the same PySAT-capable Python environment as the campaign engine.
