# Campaign notes — gate B rung 2 (VXXZ24 K100_2.37155181) certified enclosure

(Framing fixed per Main owner re-check 2026-08-30, 4-part statement. Owner
confirmed at 200 bits: cert − raw = 1.9381e-11; cert − published =
+2.0258545e-6; no-arm 2.3715518062057623 with margin 3.7942377e-9; rung
order gap 2.1375219e-4; double-count test clean.)

Run: 2026-08-30T05:2xZ. Machine: M3 Ultra, python 3.14.3, python-flint 0.9.0.
Single process, `nice -n 10`, OMP/OPENBLAS threads = 1. Reproduces from:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 10 \
      /Users/jinleic/jinleic-workspace/cs/.venv/bin/python stage_b_rung2.py
    # (append "without" for the no-Lemma-1 arm; inputs read from
    #  ../../scratch/rmmcode/data/K100_2.37155181.mat)

Both arms re-verified byte-identical on re-run from this frozen directory
(with: raw 2.3715538358350807; without: raw 2.3715518061863814).

## Stage (a) — float64 transcription, [REPRODUCED]

Transcription of the released VXXZ24 verifier (`rmmcode`, osf.io/7wgh2) to
`src/vxxz24_float.py` (dual float/interval mode). MATLAB sources read
first-hand this session. Deltas vs rung implemented exactly: 3 regions,
single `p_comp` (shared dim = region), `min(ptr.shape)==0` special case,
VXXZ24 registration order (interleaved num_retain per region), explicit
global-stage region_prop symmetry linear constraints ([0,1,-1] and [1,-1,0]).

At the released point (`params_sha256`
df75ae3acaa5b1388bd2f17cce266aec169d874c9fea1eb0722c3e25871da93d):

- parameter count: **6759 == 6759** (registered == file, exact)
- omega stored in file: 2.3715518061833021 (file name 2.37155181 = display
  rounding, consistent)
- max inequality violation c_max = 4.6349590904e-12
- max equality violation ceq_max = 3.5953524880e-11
- max linear-constraint violation = 0.0 (bit-exact)
- Schoenage line `value - 4 ln 7` = **0.0 exactly** (better than rung 1's
  -6.2e-15; the released vector satisfies the line at float64 exactly)
- value = 7.783640596221253 == target to all 15 digits

All residuals are inside the harness's own refine tolerance 1.1e-9 by 4+
orders. VERDICT: REPRODUCED; stage (b) licensed (pre_statement A2).

## Stage (b) — interval enclosure

Validations executed before any margin was computed (Addendum A3/A4):

- Per-block containment (elementwise): global num_block 9/9, penalty 3/3;
  documented float-noise tolerance 64 ulp because the float64 reference
  carries pairwise-summation noise (observed max 1.74 ulp at block r1d2;
  transcription errors would appear 9 orders above this).
- Lagrange constraints: interval path emits LOG-FORM pairs
  ln(dmv)-(tv-1); all 1512 checked form-to-form consistent with the float
  exp-form residuals (|exp-form| = dmv*|e^{-L}-1| <= |L|): 1512/1512.
- Precision sweep (proves interval path live, no float leak): num_block
  widths 1.084e-19 at prec 300/200/128/96 -> 1.084e-18 at 64; max width
  2.168e-19 -> 4.770e-18. Monotone growth; the 1e-19 widths are genuine
  enclosures, not zeros.
- Orchid test (zero-width check): value line width 0.0 is CORRECT — value
  sums pm point intervals only (same structural fact as rung 1).
- DOUBLE-COUNT TEST (CORRECTION_coincidence.md protocol), run first:
  - WITHOUT Lemma-1 charge: raw = 2.3715518061863814
  - WITH Lemma-1 charge:    raw = 2.3715538358350807
  - delta = 2.0297e-6 = the charge (lemma1_eps_max 3.0767e-6, weighted),
    so the Lemma-1 gap is present EXACTLY ONCE inside the raw endpoint
    (via the hashing-constraint P_true charge) and must not be added to
    eps again. Rung-2 semantics therefore mirror the corrected rung 1.

## Certified result (Addendum A3 option (a) semantics)

    WITH Lemma-1 (the pre-committed semantics — P_true must be charged):
      omega_cert_upper_raw = 2.3715538358350807
      eps_abs (defects only) = 1.9380875689560958e-11
      omega_cert_with_absorbed_defects = 2.3715538358544617

    WITHOUT Lemma-1 (diagnostic upper ecosystem line):
      omega_cert_upper_raw = 2.3715518061863814
      omega_with_absorbed_defects = 2.3715518062057623 < 2.37155181 (pub)

Interval components (with arm): R_sum_low = 2.8170035674609757 (widths
1e-19 each), M_low = 2.0942543887102634, target hi = 7.783640596221253.

## Verdict — the honest one (4-part statement per Main; lead with the
positive result)

**PART 4 FIRST — THE POSITIVE RESULT.** The two-rung ladder ORDER survives
rigorous arithmetic: rung-1 certifies to 2.3713400836689 and rung-2 to
2.3715538358544617 — gap 2.1375e-4 — so Alman25's rung genuinely improves
on VXXZ24's even when both are enclosed under identical semantics with
explicit slack. The first two-rung rigorous comparison in this line, and it
holds under both treatments tested (Lemma-1 charged and diagnostic-off).

**PART 1 — WHAT IS ESTABLISHED.** A rigorous interval enclosure at VXXZ24's
released parameter point, under the pre-committed Lemma-1 semantics, gives
omega <= 2.3715538358544617. That is a valid rigorous upper bound. It sits
2.0258545e-6 ABOVE their published 2.37155181, so it does not reproduce
their six-decimal rounding. (The fixed A3 threshold verdict here is an
HONEST FAIL of the rounded-reproduction test, reported weaker, not
massaged.)

**PART 2 — WHY, WITH THE ATTRIBUTION EXPLICIT.** The shortfall is entirely
the ln-space dual residual of the SHIPPED Lagrange multipliers (worst
3.0767e-6 in glob region 2; 2.1502e-6 at r0; 1.1433e-6 at r1, ~1500
per-part residuals below ~7e-8) measured against interval widths of
1.08e-19 (individual R entries: R_comp widths ~1e-10..1e-9 dominated by
Lemma-1 eps, num_block entropy widths 1e-19). The arithmetic is nowhere
near the bottleneck — five orders of magnitude tighter than the gap. The
no-arm diagnostic shows that with tight multipliers the published value IS
recoverable: without-arm certified 2.3715518062057623 (margin 3.7942377e-9
BELOW the published rounding). So this is a statement about the LOOSENESS
OF THE RELEASED CERTIFICATE DATA, NOT about the validity of the bound, and
NOT about any error in the paper.

**PART 3 — THE GENERAL FINDING.** Both rungs behave the same way: the cost
of rigour is 1.084e-6 on Alman25's rung and 2.026e-6 on VXXZ24's, dominated
in both cases by the Lemma-1 residual of the shipped multipliers rather
than by arithmetic. The published constants in the combination-loss line
are float-optimal points whose accompanying dual certificates are loose at
the ~1e-6 level; rigorously enclosing the released data costs about one to
two parts in 1e6. First quantitative measurement of the method's artifact
quality.

## Two-rung rigorous comparison (the payoff — see PART 4 above)

    rung 2 (VXXZ24)   omega <= 2.3715538358544617
    rung 1 (Alman25)  omega <= 2.3713400836689
    certified gap     = 2.1375218556e-4 (also < 2.371866, DWZ23, unenclosed)

Both certified endpoints are MACHINE-VERIFIED under identical semantics;
this answer, unlike either published number, carries its certificate.



## Provenance / evidence labels

- Stage (a) reproduction numbers: [REPRODUCED] (same inputs, same stated
  pipeline: released .mat params through transcribed verifier).
- Stage (b) certified endpoints: MACHINE-VERIFIED (interval arithmetic,
  outward-rounded, per-block containment asserted in-repo; deterministic
  re-run from the frozen scripts reproduced both arms byte-identically).
- Lemma-1 dual-residual characterization: COMPUTATIONAL-EVIDENCE at the
  dual-certificate standard of pre_statement A3 (the shipped multipliers
  are what they are; no full re-optimization within this run — see the
  timebox verdict below).
- VXXZ24/Alman25 published numbers: CITED-DEPENDENCY (papers read
  first-hand per gate A ledger).

## Optional multiplier re-solve — timebox verdict (Main's optional extra)

Attempted the cheap probe only: shifting each region's lam_sum by the mean
residual (a certified closed-form re-centring) reduces the worst global
residual from 2.150e-6/1.143e-6/3.077e-6 to only
2.120e-6/1.124e-6/3.065e-6 — the residual is genuinely SCATTERED across
shapes, not a constant offset. A full fix requires re-deriving the Lemma-1
dual (GetLambda max-entropy construction) per (region, part) group with its
own enclosure proof — a fresh campaign, not a cheap win. NOT attempted
within this run (pre-existing budget); recorded as the named follow-up
that would plausibly upgrade rung-2's verdict to "reproduced with
recomputed multipliers".
