# Part 1 — decomposition of the 2.0258545e-6 gap (run 2026-08-30, this campaign)
## REVISED ATTRIBUTION (same day, after first-order screen; see §Final at bottom)

Method: the SAME interval code path (`src/stage_b_rung2.py` aggregation, re-run
through `src/gate_c_precision_sweep.py` which calls `run_at_prec` with
`interval_core.prec()` re-entered at each working precision MID = 300, 200,
128, 96, 64 bits), at the released parameter point p* (sha df75ae3a…).
Both arms. Raw JSON in this directory:
precision_sweep_with_lemma1.json / precision_sweep_without_lemma1.json.

## Gap = certified endpoint − published 2.37155181, per precision

| MID bits | R_sum width (with arm) | gap_to_published (with arm) | gap (without arm) |
|---|---|---|---|
| 300 | 1.57e-89 | 2.0258544615181506e-06 | −3.794237812826395e-09 |
| 200 | 1.99e-59 | 2.0258544615181506e-06 | −3.794237812826395e-09 |
| 128 | 9.40e-38 | 2.0258544615181506e-06 | −3.794237812826395e-09 |
| 96 | 3.03e-28 | 2.0258544615181506e-06 | −3.794237812826395e-09 |
| 64 | 1.73e-18 | 2.0258544615181506e-06 | −3.794237812826395e-09 |

Ancillary widths same at every precision: M width 0.0 exactly, value width
0.0 exactly; eps_abs = 1.9380875689560958e-11 (float64 input, not interval).

**Primary conclusion (unchanged by revision):** the certified endpoint does
NOT move as precision drops 236 bits. The 2.0258545e-6 gap is NOT
outward-rounding accumulation. Rounding contribution to the gap measured
directly: ≤ 1.7e-18 R-width / 2.094 (M) ≈ 8.3e-19 of omega at 64 bits —
degenerate. So the gap is **structural**, not arithmetic.

**Revised attribution of the structural part.** First measurement I failed
to pin down on the first pass, now corrected: the gap is NOT the Lemma-1
dual residual of the SHIPPED multipliers (that was the 05:20:00Z campaign's
diagnosis and it is an over-statement when read without qualification).
Measured directly at p*: the *certified* enclosure at the shipped point
already sits inside the feasible polytope — the first-order screen
(stage a, below) reads omega = 0.24317919629188728 at p* with the
Schonhage line value − target = 0.0 exactly, and max_ineq = 4.63e-12,
max_eq = 3.60e-11. So at the SHIPPING point the residual-charge
interpretation collapses; the gap decomposes as:

  gap = 2.0258545e-6  =  0 (rounding, bounded by 1.7e-18 in R,
                            bounded above by 8.3e-19 after M division)
                      + 2.0258545e-6 (Lemma-1 charge, i.e. the +2·ε_IV
                        correction on the global-stage hashing penalty in
                        the with-arm; cf. frozen report glob_r0 2.1502e-6,
                        glob_r1 1.1433e-6, glob_r2 3.0767e-6 multiplied by
                        their region_props and weighted into the raw
                        endpoint)

**and the charge is an artifact of the CERTIFICATION SEMANTICS at the
shipping point, not of the released certificate data** — the without-arm
run proves the published rounding IS reachable (−3.7942e-9 margin at
every precision) under the transcribed code path. What the released
multipliers are: they are GOOD ENOUGH for the paper's own float64
verifier (which does not charge Lemma-1 residuals — it runs the exp-form
residuals at 1.1e-9 tolerance and passes at 3.6e-11), but they are NOT
tight enough for our stricter ln-space certified standard. That is a
statement about the DIFFERENCE IN STANDARDS, not a defect of the release.

Evidence labels: MACHINE-VERIFIED (all interval quantities recomputed
inside flint Arb at each stated precision from frozen scripts; the
precision sweep JSON in this directory re-runs byte-identical).
STRUCTURAL attribution of the 2.0258545e-6 gap: REPORTED-then-REPRODUCED
here at 5 precisions in this campaign; the interpretation quoted from
campaign 05:20:00Z is commented on above under exactly that label.

## §Final (authoritative, supersedes the intermediate reading above and
##  the 05:20:00Z Part 2 wording where they disagree)

gap(precision) = 2.0258544615e-6 CONSTANT across MID ∈ {300,200,128,96,64}.
Not rounding. Attribution by the double-arm experiment, measured here:
  • with Lemma-1 charge (pre-committed certified semantics): +2.0259e-6
  • without (diagnostic): −3.7942e-9, i.e. published rounding recovered.
The 2.0259e-6 is exactly the ln-space Lemma-1 residual of the released
multipliers (glob-region ε's 2.15e-6 / 1.14e-6 / 3.08e-6 before region_prop
weighting), charged because our certified standard is STRICTER than the
released float64 verifier's own standard (1.1e-9 exp-form tolerance). The
release is not defective relative to its own standard; our certified
standard absorbs a real slack the released verifier never had to.

## B0 delta note (requested by Main, 2026-08-30)

This campaign's B0 re-run of the rung-2 certified endpoint reads
2.3715538358544612 vs the frozen 05:20:00Z value 2.3715538358544617:
delta 5.0e-16, ours the SMALLER. Cause: the frozen value was carried
through one more outward-rounded float64 presentation step; both are the
same certified quantity computed under identical semantics, the re-run is
the tighter outward rounding, and the frozen number is not retracted —
the movement is in the safe direction (bound got tighter). MACHINE-VERIFIED
both runs; nothing else moved (gap to published identical to 2.0259e-6).
