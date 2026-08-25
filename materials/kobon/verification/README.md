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
| `../release/kobon-2026-08/verification/sector_bound_audit.py` | Exact-rational audit of four/two/zero sector bounds, shared-pair rays, single-line endpoint closure, and shared-edge opposite-side parity over deterministic generic, parallel, concurrent and mixed arrangements. | Exit `0` only with zero violations and a valid four-sector sharpness control; any mismatch raises `AssertionError`. | `"status": "PASS"` with 3,600 arrangements, 1,607 endpoint collisions, 84 colliding same-ray cases, 42 single-line shared segments, 1,591 shared-pair segments, and zero violations. |
| `../release/kobon-2026-08/verification/reified_face_audit.py` | Fixes every `P/C/X/S` literal from the exact \(K(10)=25\) and \(K(12)\ge38\) rational witnesses, compares the complete direct-face set, and checks reified endpoint/perturbation/K(4)/hereditary CNF acceptance. | Exit `0` iff both face sets match exactly, both fixed CNFs are SAT, and the multipoint-incidence lower bounds hold. | `"status": "PASS"`; 25 and 38 exact faces, with 0 and 8 multipoint-incident faces. |
| `../release/kobon-2026-08/verification/chirotope_gp_audit.py` | Replays the dual determinant sign map with parallels/concurrencies, all zero-aware rank-three Grassmann--Plücker relations, and a base-relaxation one-sided-product probe. | Exit `0` iff exact signs/relations have zero mismatch and any found base violation is rejected. | `"status": "PASS"`; 195,000 triples / 1,386,000 relations; `NO_GAP_FOUND_THROUGH_N7_PROPAGATION_ONLY`. |
| `../release/kobon-2026-08/verification/pappus_relaxation_probe.py` | Shows that the abstract order model admits one guarded non-Pappus pattern, then adds the three projective conclusion clauses and resolves the same assumptions. | Exit `0` iff the base assumptions are SAT and the cut version is UNSAT. | `"status": "PASS"`, both Boolean verdict fields `true`. |
| `../release/kobon-2026-08/certificates/ladder_n10_t26_faces_{certificate.json,cnf,dratcheck.log}` | Binds the byte-reproducible faces-only target-26 CNF to the complete external 14.05 GB DRAT and its independent checker result. | Regenerated CNF must compare byte-identically; checking the hash-pinned DRAT with `drat-trim -w -t 200000` must return `s VERIFIED`. | 131,314 variables / 792,462 clauses; 52,808 formula clauses and 40,395,481 lemmas in core; 3,464,683,123 resolution steps; `s VERIFIED`. |
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
