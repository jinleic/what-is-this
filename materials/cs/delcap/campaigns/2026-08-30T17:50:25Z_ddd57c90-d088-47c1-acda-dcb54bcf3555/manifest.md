# Manifest — Orbit-Certificate Campaign (Open Item 1: (3,10) + q=3 ladder)

Campaign directory: `2026-08-30T17:50:25Z_ddd57c90-d088-47c1-acda-dcb54bcf3555`
Owner: DelcapNextGate. Frozen: 2026-08-30T19:5xZ (this manifest written after
all 20 pre-registered rows completed; nothing in this directory was edited
after the run — the resume-patch to `gate_c_orbit.py` happened BETWEEN runs
and is itself part of the frozen script, hash below).

## Contents (byte-copies in place, sha256 in `checksums.sha256`)

| file | role |
|---|---|
| `pre_statement.md` | committed BEFORE the first computation: open item named, group/object line (Main's check), falsifiable outcome, adjudication rule, budget, box, anti-numerology clause. One post-commit addition: the GROUP-OBJECT LINE requested by Main mid-run — inserted verbatim with its confirmation evidence, and the falsifiable-outcome section was re-worded around it; no adjudication rule, budget, box or stop-rule was changed. |
| `gate_c_orbit.py` | the sparse orbit-certificate pipeline, sha256 `d17c2b43…` (see note below on the resume patch) |
| `orbit_rows.jsonl` | 20 rows: (3,10) x d in {1/2,1/5,1/10,1/20} + q=3 ladder n=6..9 x the same 4 d's. Each row: certified [primal, dual] per symbol, width, Tavakoli LB1/LB+/UB, verdicts. |
| `crosscheck.json` | dense-vs-orbit cross-validation at (3,6),(3,7) d=1/2: OVERLAP OK both (orbit lower endpoint matches dense to 14 digits) |
| `inf_branch_check.json` | the frozen invalid-certificate branch re-verified VERBATIM at q=2,n=10,d=1/20: raw per-word snap at 2^-40 zeros exactly 1 reachable word (mass 9.766e-14 < 9.095e-13), `cert_dual = +inf` (no claim) — the branch Bitu ucak preserves upper-bound trustworthiness is live in the orbit path |
| `log.txt` | full run log (selftest anchors A1-A6, inf check, crosschecks, 20 row lines) |
| `run_stdout.log` | raw stdout of all runs, including the two earlier stalled runs and their diagnostics |
| `log_prev_run.txt` | preserved log of the pre-fix stall (kept for the diff record, rule 17a spirit) |
| `tool_versions.txt` | python 3.14.3, python-flint 0.9.0, numpy 2.5.2, mpmath 1.3.0, Arb 400 bits |

## Two bugs found and fixed DURING the campaign (recorded, not hidden)

1. **`snap` remainder loop**: the frozen one-at-a-time largest-remainder loop
   spins ~1e11 iterations when the float vector's sum sits far from 1
   (partial-alphabet marginals). Replaced by the cycle-decomposed equivalent
   (add `rem//n` to all, then +1 to the top `rem mod n` fractionals — provably
   identical modulo cycling symmetry, verified equal on 200 random + 60
   adversarial trials against the frozen loop byte-for-byte).
2. **ylist-only normalization** of the dual D' candidate rehearsal produced a
   NEGATIVE width once (dual below primal — an invalid certificate shape).
   Root cause: the per-word marginal was snapped over only the words the
   representatives emit; the missing relabeled words carry real mass, so the
   normalized D' was sub-stochastic. Fix: extend the G-invariant float
   marginal to the FULL output alphabet through the orbit identity before
   snapping. Recorded here because it is exactly the failure shape the frozen
   SOLVER_AUDIT warns about (an invalid certificate that looks clean), caught
   by the width-sign check, not by printing.
3. **granularity of the +inf check (amendment)**: the raw-snap +inf branch
   lives at PER-WORD granularity; my first re-verification snapped
   orbit-AGGREGATED masses, which never zero out — the check as first coded
   returned a finite dual and would have FAILED its own assert. Fixed to the
   frozen per-word construction (verified: 1 zeroed word, cert_dual = +inf).

## Resume protocol

Runs are crash-resumable: rows already in `orbit_rows.jsonl` are skipped
(see `resume skip` lines in the log). The `(3,10,d=1/2)` row was produced by
the pre-resume script; all remaining rows by the resume-enabled script; the
resume patch changes ONLY the row-loop skip logic, not any certificate
computation.

## Rule-7 sweep sentence (mandatory)

**Swept set exactly:** `q = 3`; rows `(q,n,d)` = `(3,10)` at `d` in
`{1/2, 1/5, 1/10, 1/20}` and `(3,n)` at `n` in `{6,7,8,9}` at the same four
`d` — 20 rows, each by exact-rational orbit certificates (primal = exact
mutual information of an explicit rational input; dual = max over input
g-orbit representatives of KL(W(.|x)||D') for an exact rational G-invariant
D'), Arb 400 bits, snap 2^30. **Not swept:** `q >= 4` at `n > 5`; any `n > 10`;
Morozov-Duman `m > 23`; the Lambda-search gap at `m` in `{22,23}`; Pinto-
Ribeiro `n` in `{29,31}`; the `n -> infinity` capacity `C(BDC_d)` — every row
here is a finite-`n` theorem about `C_{3,n}(d)` and says nothing asymptotic
(rules 7/9/14).
