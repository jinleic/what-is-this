# Campaign notes — gate A closure + gate B stage (a) complete, stage (b) incomplete

Run: 2026-08-29..30, machine M3 Ultra, python 3.14.3, python-flint 0.9.0,
sympy 1.14, single-user session; BLAS threads pinned at import (scipy default
was used in stage (a); no BLAS-heavy ops in the program — entropies only).

## Gate A — verdict:record artifact not released (dated finding, 2026-08-29)

Primary reads complete (all four papers, full texts; quotes in
2608.16884_main.tex and 2404.16349_numerical.tex):
- arXiv:2608.16884v1 publishes the reformulated combination-loss program
  (Eq. 11), the H^max certificate scheme (Lemma 1), and states q=5, l*=4,
  ~7e6 parameters. NO parameter table. NO arXiv ancillary files (assets/ has
  only logos). Verification section verbatim: "We are preparing a repository in
  which we will release the verification code and our discovered solution."
- Search ledger (2026-08-29): author GitHub accounts (orbitingflea, EmilienDupont,
  marvineisenberger[404], bkozl[404]); google-deepmind/alphaevolve_results
  (May-2025 notebook only); GitHub/OSF text search for "2.371177": no repo.
=> Gate A verdict: as of 2026-08-29, the parameters attaining omega < 2.371177
   are not published at any precision; the bound is not independently
   checkable from the paper as published. NOT a claim of wrongness or
   impropriety; author repository explicitly promised.

## Corrigenda established by this campaign (framing correction)

All three predecessors DID release code + parameters on OSF:
- VXXZ24 osf.io/7wgh2 -> data/K100_2.37155181.mat (6759 float64 params)
- Alman25 osf.io/mw5ak -> data/W1.00_2.371339.mat (24855 float64 params)
- DWZ23 osf.io/dta6p -> fmm-code.zip with MATLAB MCOS object serialization
  (unreadable without MATLAB; rung NOT evaluated — stopped per directive).
The ladder is therefore reproducible in floating point; what was never done,
before this campaign, is a rigorous interval enclosure of any rung.

## Gate B — stage (a) COMPLETE for rung 1 (Alman25, omega <= 2.371339)

Transcription of the released verifier (~1.3k MATLAB LOC across 11 files) to
src/alman25_float.py. Float64 evaluation at the released point:
- registered parameter groups == file params == 24855 (exact match)
- omega in file 2.3713389005434182; published 2.371339; |delta| 9.95e-8
  (display rounding, consistent)
- max c violation 1.137e-10; max ceq violation 2.390e-11; Schönage line
  value - 4*ln(7) = -6.2e-15 (constraint requires >= 0 → violated by half-ulp
  noise)
- verdict: REPRODUCED within their refine tolerance 1.1e-9 (order-of-magnitude
  inside). Per pre_statement Addendum A2, stage (b) may now be trusted about
  the BOUND if it fails; per Addendum A3 any certified endpoint must carry the
  explicit absorbed slack eps.

## Gate B — stage (b) INCOMPLETE (infrastructure, not a bound result)

src/alman25_interval.py runs the same evaluation tree under python-flint
interval semantics, but the R_sum/M/eps aggregation produces -4.27 for
R_sum_lower, which is arithmetically impossible for this instance (the true
retained sum is ~+2.8): the interval plumbing mixes numpy-array and arb
objects in the hashing aggregation and p_compY/Z handlers. Per Addendum A2's
interpretation rule and Main's directive: this is TRANSCRIPTION-stage evidence
about MY PLUMBING, not about the published bound. No certified endpoint is
claimed. Tasks to finish (all mechanical, ~1-2 sessions):
1. route p_compY/Z hashing marginals through the interval Vec path (they are
   currently numpy objects in the aggregation);
2. fix the globals-side min-over-dims retained count to use interval lowers;
3. charge eps_total (Lemma-1 residual, already verified at 3.7e-6 order) into
   the P_true upper bounds and re-check the Schönage line with R at certified
   lows;
4. report omega_cert(PASS/FAIL @ publish rounding) + explicit eps per A3.

## Provenance / evidence labels

- All quotes: MACHINE-COPIED from primary TeX/HTML files in this campaign dir.
- float64 reproduction numbers: [REPRODUCED] of the released program at the
  released point (notанием [REPRODUCED] per repo rule 13 — same inputs, same
  stated pipeline).
- Stage (b) numbers in alman25_interval.py output: NOT EVIDENCE (plumbing bug);
  ignored entirely.
- Predecessor-release discovery: HUMAN-AUDITED (OSF pages + zip contents
  fetched and inspected; sha256 listed in stage_a_result.json).
