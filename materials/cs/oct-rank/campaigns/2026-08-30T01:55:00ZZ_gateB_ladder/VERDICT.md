# Gate B frozen record — rank-24 attack: FAILED (no rank-24 candidate)

UTC freeze: 2026-08-30T01:55Z. No claim about R_R(T_O) >= 24 follows from
anything in this artifact (README rule; pre_statement.md Gate B rule).
This is COMPUTATIONAL-EVIDENCE about the declared search only.

## Budget actually consumed (declared in pre_statement, 2026-08-29)

- Engine: CP-ALS (2000 sweeps) + LM/trf-lsmr (300 rounds), float64, scipy
  1.18.1/numpy 2.5.2, threads pinned 1, `nice -n 10`.
- Random sweep: 256 starts at r=24 (seed stream PCG64 [20260829, 24, 0..255]).
  Wallclock 26.8 min, well within the 8 h/rank cap; ranks 23..20 NOT run
  (redirected per owner approval 2026-08-30).
- Structured probes (owner-directed): 5 compress-merges, 4 alphabet snaps,
  1 continuation path (7 scaled-s rows + 1 true-24 endpoint).
  Total probes: 256 + 5 + 4 + 8 = 273; every probe holds rel-trajectory +
  factor-norm records (joint reporting).

## Results — the informative contrast (joint residual/norm rows)

| rank | probe | best rel | factor norms at best | verdict |
|------|-------|----------|----------------------|---------|
| 25 | control (cert x0, LM 41 rds) | 1.05e-15 | (0.67, 13.72, 8.70) | exact cert exists |
| 25 | random start (ALS 2000 + LM 287) | 1.85e-14 | (2.30, 9.02, 9.46) | converges, bounded |
| 24 | random, 256 starts | 2.63e-03 | (9.34, 30.42, 18.23) — start 139 | NO convergence; med. rel ~3e-2 |
| 24 | compress (merge cols 11,23 -> true 24 cols) | 1.76e-02 | (1.32, 36.73, 22.08) | NO convergence |
| 24 | alphabet snap ({0,±1/2,±1,±2}) | rel ~3.0 | — | fails outright |
| 24 | continuation term 0 (path, 25-col) | 1.2e-15..1.5e-15 | norms explode mid-path: 13.7 → 1135.7 (s=0.01) | re-convergence, NOT rank-24 |
| 24 | continuation s=0, TRUE 24 columns | 5.19e-01 (start 0.563) | (0.66, 13.72, 8.67) | LM runs away, bounded |

- ALL 256 random LM runs consumed the full 300 rounds (max_nfev), i.e. the
  optimizer never signalled xtol convergence at r=24.  The r=25 control
  converged in 287 rounds with status=3.
- Divergence-without-convergence signature present (start 204: resid
  1.36e-01, maxnorm 388.6; start 236: 1.15e-01 / 95.6).
- Continuation path (term 0 scaled by s, 25-col bookkeeping at every row):
  machine-precision residual at every s, but the solution manifold returns
  to a rank-25-type point (norms reset at s=0).  At the true-24 endpoint
  (column removed) LM DEGRADES (0.563 → 0.519) with bounded norms.  This is
  the "least informative clean negative" of the three owner-specified
  signatures — a 25-col configuration is never labelled rank 24 anywhere.

## Escalations

- One sub-1e-13 event (merge-mode s=0.9, rel 1.75e-15) escalated to Main
  BEFORE recording (owner replied; correct bookkeeping: 25th term still
  present at s=0.9 → re-convergence, not rank-24).  Owner-approved
  reframing to border-rank probe followed; result above.

## Appendix — continuation probe, additional terms (owner-approved)

**Appendix (post-freeze, owner-approved, ≤5 terms):** `gate_b_continuation_terms_appendix.jsonl`
— same protocol deflating terms 5, 11, 17, 23, 24. Rows labelled with the
same 25-col/24-col bookkeeping; adds basin-diversity evidence. FROZEN at
2026-08-30T02:05Z; hashes appended to `code_hashes.sha256` (appendix block).

Appendix result (16 rows, all LM status=3): terms 5/11/17/23/24 all show
the same clean-negative signature. 25-col rows re-converge when the term's
scale is shrunk (terms 0 and 5: machine precision; terms 11/17/23/24 land
in a degraded basin at rel 0.49-1.27, showing the certificate's columns
are NOT freely interchangeable — the certificate point sits in a narrow
basin). True-24 endpoints degrade (rel 0.49..1.57) with norms staying
O(14). No rank-24 candidate; no divergence-at-machine-precision (that
signature only ever appeared on 25-col rows, i.e. re-convergence).

## Conclusion

GATE B VERDICT: FAILED (no rank-24 candidate found) under 273 probes
(256 random + 17 structured).  The window 18 ≤ R_R(T_O) ≤ 25 remains open
at both ends.  No statement about R_R(T_O) ≥ 24 is made or implied.

## Files

- gate_b_search3.py — random sweep engine (256 seeds r=24; also ran the
  r=25 1-start control).
- gate_b_rank24_results.json / .npz — full per-start record (residual,
  factor norms, trajectories; 256 rows).
- gate_b_rank25_results.json / .npz — control record.
- gate_b_deflate.py — compress-merge + alphabet + scale-continuation tools.
- gate_b_continuation.py + gate_b_continuation_term0.jsonl — the
  border-rank-probe rows (joint residual/norm, phase-labelled).
- gate_b_ladder_report.py — regenerates the summary table from the jsons.
- code_hashes.sha256 — sha256 of everything above.
- gate_b_continuation_terms_appendix.jsonl — appendix rows (terms 5/11/
  17/23/24), FROZEN 2026-08-30T02:05Z.
