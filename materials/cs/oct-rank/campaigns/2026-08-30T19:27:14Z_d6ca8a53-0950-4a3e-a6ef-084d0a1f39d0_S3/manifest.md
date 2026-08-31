# manifest.md — S3 campaign (agent OctRankS3, 2026-08-30)

**Frozen 2026-08-30.** Do not edit after launch. All scripts byte-copied
below; SHA256 in `code_hashes.sha256`.

## Pre-registration

`pre_statement.md` — routes, budget (4 h wall-clock; 30 min Route-E cap),
adjudication rule, anchor suite, and the "failure to reach 14" semantics —
committed BEFORE the first S3 computation. The anchor scripts ran first
and PASSED before any Route was attempted (rule 14, harness validity).

## Files

| file | role |
|---|---|
| `pre_statement.md` | Pre-registration (read first) |
| `s3_anchor.py` / `.out` | Anchors 1-3 exact (n=2=>3, n=4=>8, tau=7 CONTROL); exact certificate parsing fix documented inline |
| `s3_routeA_strassen.py` | SUPERSEDED first-draft of Route A (rule-5 retraction: bookkeeping draft with a wrong derived pairing; kept as a byte-copy) |
| `s3_routeAB.py` / `.out` | Route A corrected: Strassen commutator rank, both formula twins, 173 triples |
| `s3_routeB_AFT.py` / `.out` | Route B: AFT/substitution pivots, double-peel, flattening caps |
| `s3_routeC_substitution.py` / `.out` | Route C: residual-pencil floor = 12-exact (C4) => 13 ceiling |
| `s3_routeD_flattening.py` / `.out` | Route D: flattening ranks 3/8/8 (supersedes the first version's mode-(c) slicing error; retraction in VERDICT) |
| `s3_routeE_groebner.py` / `.out` | Route E: pre-registered BUDGET-FAIL; linear-feasibility diagnostic only |
| `s3_routeF_subalgebra.py` / `.out` | Route F part 1: H' = span{e0..e3} quaternion subalgebra; block-diagonal / off-diagonal action split (exact) |
| `s3_routeF2_blocks.py` / `.out` | Route F part 2: complement-block identification (per-block structure) |
| `s3_routeF3_global.py` / `.out` | Route F part 3: global T-conjugation exact; upper 14 for the (1, i, j) triple; lower 13; the 13-or-14 question open |
| `VERDICT.md` | Outcome, tables, retractions, obstruction, rule-7 scope |
| `code_hashes.sha256` | SHA256 of all frozen files |

## Reproduce

    PY=/Users/jinleic/jinleic-workspace/cs/.venv/bin/python
    cd oct-rank/campaigns/2026-08-30T19:27:14Z_d6ca8a53-0950-4a3e-a6ef-084d0a1f39d0_S3
    for s in s3_anchor s3_routeAB s3_routeC_substitution s3_routeB_AFT \
             s3_routeD_flattening s3_routeE_groebner s3_routeF_subalgebra \
             s3_routeF2_blocks s3_routeF3_global; do
      $PY $s.py
    done

All 9 exit 0 (verified at freeze). Environment: CPython 3.14.3,
python-flint 0.9.0 (fmpq exact), sympy 1.14.0 (Route E diagnostic only);
single-threaded; no floats in any load-bearing check (the Route-E
diagnostic is labelled necessary-only in its own output).

## Escalation status

No claim here improves or contradicts a published result: the best
certified bound on the 3-slice family remains 13 (gap to 14 = 1); the
window 18 <= R_R(T_O) <= 25 is untouched. The Route-A twin (16 = 8 + 8)
is explicitly labelled NOT a certified universal bound (sweep-only + twin
contradicted by the tau control at n=4 + Strassen-form provenance
unverified), so per the campaign rules NOTHING required escalate-to-Main
before writing: this artifact records a NEGATIVE with a named, quantified
obstruction. Each labelled note is laid out for the owner's audit in
`VERDICT.md` and `../../scratch/s3_three_slice_notes_HUMAN-AUDIT-PENDING.md`.
