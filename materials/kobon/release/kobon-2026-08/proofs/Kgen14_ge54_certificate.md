# Certificate: K_gen(14) >= 54 — VERIFIED 2026-08-20

## Source

- Repository: https://github.com/rufio72/kobon_triangles_k14
- Commit: `e47c7cfd9661e54d29befb83118bf83d5f804028` (2026-08-12)
- File: `sol1/lines_rational.json`
- File SHA-256: `83fc26666d73cd39791018a72aba71caeac4a3f7efb27109ed2955b4540b2523`
  (matches repository `SHA256SUMS.txt` manifest)
- Archived copy: `scratch/kobon/n14/sol1_lines_rational.json` (byte-identical)
- Engine-format lines + the certified 54-triangle selection:
  `scratch/kobon/n14/sol1_engine_check.json`

## Verifier independence matrix (all exact rational arithmetic)

| # | Implementation | Result |
|---|---|---|
| 1 | Repo authors' `verify_direct_exact.py` (re-implemented equivalently) | 54 = claimed 54 |
| 2 | Independent cevian-aware open-interior counting (custom, this workspace) | 54 |
| 3 | Campaign-audited `math/kobon/engine.py::verify_selection(n=14, sel, minimum=54)` | `True ok` |
| 4 | Distinct-line and concurrency census | 14 distinct lines; N3 = 2 triple points, no other degeneracy, Q = 0 |
| 5 | Full sweep of sol2..sol15 with SHA manifest checks | 54 each, hashes match, lines distinct |

Duplicates/candidate caveats: none observed; selection size 54 with pairwise
interior-disjoint open triangles in homogeneous exact rational geometry
(implementation #3's convention, same core as the campaign's Theorem 1/TM
lower-bound certificates).

## Scope statement (advisory-corrected)

- This certificate proves **K_gen(14) >= 54** under this campaign's broad
  convention (parallel lines, multi-points, crossed selections all permitted;
  the witness itself is non-simple with N3=2, Q=0).
- Equality `K(14) = 54` in the standard face/straight-line convention
  additionally uses the OEIS-published upper bound U=54 (BBL lineage); that
  upper-bound chain was not re-proved here and its convention scope is
  documented in `scratch/kobon/litrefs.md` (simple affine pseudoline faces).
- The broad upper frontier is therefore target **55**: a model for T=55 would
  violate the current literature upper bound; its CNF UNSAT (or the
  capacity-theoretic argument) remains open as the next upper-bound object.
- Consistency note (no overclaim): the witness is NON-SIMPLE with two
  shared-line triple points (0,4,10) and (0,11,12), Q=0. Concurrency is NOT
  proven required for 54 under this campaign: Theorem M closes only the
  no-concurrency Q=3 branch, Blanc caps simple faces at 53, and the broad
  no-concurrency Q0/Q1/Q2 crossed-selection cubes at target 54 remain
  capacity-admissible (solvers still running).

## Reproduction

One fail-closed command (exit 0 on success, exit 1 with JSON diagnostic on
any failed check):

```
scratch/kobon-audit/venv/bin/python \
    scratch/kobon/n14/maiorana_verify.py \
    scratch/kobon/n14/sol1_lines_rational.json
```

`scratch/kobon/n14/maiorana_verify.log` is the verbatim exit-0 transcript of
exactly this command (seven JSONL stage lines ending `CERTIFIED`).
