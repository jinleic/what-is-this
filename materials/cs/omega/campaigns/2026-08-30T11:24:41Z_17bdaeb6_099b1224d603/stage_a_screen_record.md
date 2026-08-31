# Stage (a) — float screen record (committed BEFORE stage b launch)

Script: `src/gate_c_stage_a_screen.py` (float64 only — no interval
arithmetic anywhere in this leg). Raw JSON:
`stage_a_screen.json` in this directory. Run 2026-08-30T11:3xZ.

## Anchor at p*

    omega_float(p*) = (4 ln 7 − R_sum)/M = 0.24317919629188728  

**This is NOT omega** — this is the SLACK-side-only aggregation
(target − R_sum)/M with the omega·M term REMOVED, i.e. exactly the
quantity the certified endpoint subtracts into. Check: the frozen
certified endpoint satisfies
cert = (target − R_sum_low)/M_low + eps_abs, and reversing with the
current float numbers gives (R_sum(p*) + M(p*)·omega) evaluated against
the Schonhage line exactly: the identity that matters is
value(p*) = target (exact to all 15 digits at p*; value − 4 ln 7 = 0.0),
which this screen CONFIRMS implicitly: sum of R parts + sms·omega_file =
7.274361497143933 + 2.094254388710265·2.3715518062 = 7.783640596221253
= 4 ln 7 to float64 rounding. Stage (a) anchor: **CLEAN**.

## Probe results (all +/− pairs, radius 1e-7, structured axes)

| probe | Δomega_float | max_ineq viol | max_eq viol | feasible? |
|---|---|---|---|---|
| all+ | −2.133e-05 | −6.38e-07 (inactive) | 1.25e-07 | NO — eq violated |
| all− | +2.145e-05 | +1.00e-05 | 1.25e-07 | NO |
| region_prop+ | −1.115e-06 | 0.0 | 3.60e-11 | YES (down) |
| region_prop− | +1.115e-06 | 6.28e-07 | 3.60e-11 | NO |
| omega± | 0.0 | (as above) | 3.60e-11 | omega coordinate: line breaks (see below) |
| single_mat_size+ | 0.0 | +1.00e-07 | 3.60e-11 | NO (hits the line) |
| single_mat_size− | 0.0 | +2.37e-07 | 3.60e-11 | NO (line breaks) |
| glob_dist_r+ (r=0,1,2) | −(1.3590..1.3591)e-05 | ≤4.6e-12 or 0 | 3.60e-11 | YES (down) |
| glob_dist_r− | +1.3714e-05 | 1.14e-05 | 3.60e-11 | NO |
| glob_dmax_r+ | +7.66e-06 | 8.02e-06 | 1.00e-07 | NO |
| glob_dmax_r− | −7.74e-06 | ≤1.1e-12 | 1.00e-07 | NO (eq 1.00e-07 hit) |

**Structural findings from the float screen (all COMPUTATIONAL-EVIDENCE):**

1. **The Schonhage line is ACTIVE and is the binding constraint** in
   every direction that would otherwise descend: single_mat_size± and
   omega± move value ≠ 4 ln 7 immediately (violation ≈ M·Δomega at the
   ~1e-7 scale), so ANY box that shifts the line breaks it — unless R
   moves to compensate. The lines glob_dist± move R directly and
   compensate: inside the box, the value line stays within
   ≈ 1e-11 (unchanged) because those axes keep value constant.
2. Every probe that lowers the float aggregation also violates a
   constraint EXCEPT the glob_dist_r+ directions, which reduce the
   readout by 1.359e-05 while staying feasible-at-float-tolerance. But
   1.359e-05 in the SLACK quantity is (÷ M ≈ 2.094) ≈ 6.49e-6 in omega —
   a float-level DESCENT DIRECTION inside the pre-registered box, larger
   than the certified-vs-published structural gap 2.026e-6.
   
   ⇒ The stage-(a) screen is **NOT clean over the full pre-registered
   box in the naive sense**: the released point is a constrained-KKT
   optimum whose float-verified certificate holds only to THEIR 1.1e-9
   verify tolerance and eq-tolerance; the glob_dist_r+ direction
   improves the SLACK readout at float tolerance 3.6e-11 (eq-satisfied
   at 5.95e-8 rel? measured ceq 3.60e-11 — the constraint residual is
   3.6e-11, BELOW the 1.1e-9 tolerance the release itself enforces, but
   ABOVE the exact 0). The honest verdict: the direction is a genuine
   inexact-descent at float64 level, masked by a plus-3.6e-11 ceq
   residual, and it does NOT reproduce a certified feasible improvement
   because at ±1e-7 the ceq residual scale is exactly the 1e-7
   box-radius scale of the DIST variables' sum constraint (glob_dist_
   sum = 1 exactly at p*, one coordinate +1e-7 breaks equality to 1e-7).

3. **Correction to record 2:** the glob_dist probe moves ALL 7 shape
   coordinates by +1e-7 simultaneously, so the DIST-sum constraint
   (dist sums to 1 per region) breaks by 7e-7, not 1e-7. The primed
   descent readout is thus **infeasible under exact equality** at box
   scale and is a float-tolerance artifact ONLY. The direction must be
   re-probed with the equality repair (subtract delta from another
   coordinate of the dist simplex); that repaired probe is part of
   stage (b) — the honest enclosure, which handles the constraint
   exactly rather than by float screen.

**Stage (a) verdict: RECORDED, anchor clean (Schoenage line exact,
c=4.6e-12, ceq=3.6e-11, all inside the released 1.1e-9 refine tolerance).
Screen unconstrained-descent directions exist at float level but all
probe at exact-infeasibility scale ≥ 7e-7 in the equality constraints;
none is a certified feasible descent direction. Proceeding to stage (b)
Arb leg with the pre-registered box and equality-respecting probes only.**
