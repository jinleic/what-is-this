# Stage (b) — interval box screen record (15 boxes, run 2026-08-30T11:4xZ)

Script: `src/gate_c_stage_b_boxes.py`. Raw JSON: `stage_b_boxes.json`.
Per-box interval evaluation under the pre-registered coupling: box centre
= p* + offset, radius r = 1e-7 per coordinate, every parameter group
evaluated as a genuine interval [c_g − r, c_g + r] (exact-dyadic
endpoints from float64, no ball radii). Aggregation semantics identical
to the frozen rung-2 protocol (interval lows for num_block/p_comp, highs
for penalties; Omega = (target_hi − R_sum_low)/M_low + absorbed defects).
Zero-crossing guard fires (FLAGGED, not silent, per pre-registration §8):
support_zero_crossing = false on every box (no coordinate crossed 0).

## Anchor (B0, offset 0, radius 0 — MUST and DOES reproduce rung-2)

    omega_cert_upper_raw = 2.3715538358350803   (frozen: …0807; Δ = 4e-16,
        float64 print granularity, machine-identical to 15 digits)
    omega_cert_upper     = 2.3715538358544612  (frozen …4617)
    R_sum_width = 7.85e-90, M_width = 2.49e-21, value_width = 0.0
    ANCHOR REPRODUCED (within float display granularity).

## Box results (r = 1e-7 per coordinate; all offsets ±1e-7)
| box | omega_cert_upper_raw | omega_cert_upper | eps_abs | beats p* cert? | beats 2.37155181? |
|---|---|---|---|---|---|
| B0 anchor | 2.3715538358350803 | 2.3715538358544612 | 1.94e-11 | — (anchor) | no |
| region_prop+ | 2.3715527208361142 | **2.3715530208533490** | 3.00e-07 | YES | no |
| region_prop− | 2.3715549508347150 | 2.3715552508525004 | 3.00e-07 | no | no |
| omega+ | 2.3715538358350803 | 2.3715540358522480 | 2.00e-07 | no | no |
| omega− | same | same | 2.00e-07 | no | no |
| single_mat_size+ | 2.3715538358350803 | 2.3715540623339795 | 2.27e-07 | no | no |
| single_mat_size− | same | same | 2.27e-07 | no | no |
| glob_dist_0+ | 2.3715483090267977 | 2.3715537111130380 | 5.40e-06 | YES | no |
| glob_dist_0− | 2.3715625773994073 | 2.3715680314948306 | 5.45e-06 | no | no |
| glob_dist_1+ | 2.3715486300180350 | 2.3715540321070225 | 5.40e-06 | no | no |
| glob_dist_1− | 2.3715625774755550 | 2.3715680315794896 | 5.45e-06 | no | no |
| glob_dist_2+ | 2.3715480143808523 | **2.3715534164727120** | 5.40e-06 | YES | no |
| glob_dist_2− | 2.3715625773111220 | 2.3715680314195837 | 5.45e-06 | no | no |
| all+ | 2.9228921434043660 | 2.9228976965146285 | 5.55e-06 | no | no |
| all− | 3.0841684138444240 | 3.0841737974899637 | 5.38e-06 | no | no |

Interval widths: R_sum widths 7.85e-90 … 1.57e-89 (all boxes ≪ 1e-8),
M_widths 2.49e-21 (all boxes), value_width 0.0 (all boxes — value sums pm
point intervals only, structural as in rung 2). eps_abs — the box
defect-absorption term — dominates the certified-vs-raw margin everywhere
(1.9e-11 at B0, up to 5.5e-06 at glob_dist boxes) because the box moves
the FLOAT64 evaluation off the exact KKT point.

## Reading (COMPUTATIONAL-EVIDENCE for direction, MACHINE-VERIFIED for
## each certified endpoint)

1. **No box in the enclosed family produces a certified bound below the
   published 2.37155181.** Closest: region_prop+ 2.3715530208533490
   (1.21e-6 ABOVE published), single-axis only; the structural gap of
   Part 1 (Lemma-1 + box-defect absorption) exceeds the raw/eps gain.
2. **Certified improvement over the p* enclosure exists in 3/14 signed
   boxes**: region_prop+ (−8.15e-7 vs p*), glob_dist_0+ (−1.44e-7),
   glob_dist_2+ (−4.39e-7). These are certified (interval lows charged
   with outward rounding); they are NOT world-record claims (they do not
   beat the published rounding), they are LOCAL improvements within the
   pre-registered box family at the pre-registered radius.
3. **Per Main's standing instruction and my pre-registration §8:** these
   improvements write nowhere outside this campaign record; no claim of
   global improvement; the two-rung certified ladder numbers from the
   05:20:00Z campaign remain the authoritative rung endpoints.

##REQUIRED NEXT (pre-registered §3): the certified-improvement direction
## region_prop+ (raw 2.3715527208361142) establishes that, within this box
## family, the reported point is not the certified optimum at r = 1e-7
## granularity: a -8.15e-7 certified shift exists along region_prop+.
## This must be recorded as the headline gate-C fact either way: the
## enclosure now separates the reported point from the certified optimum
## within ±1e-7 of region_prop. No escalation (no bound below 2.37155181).

## RETRACTION IN PLACE (Rule 5) — 2026-08-30T12:2xZ

The verdict paragraph above ("The reported point is not the certified
optimum... three boxes certify strictly below B0") is RETRACTED as a
bound-claim. Feasibility audit (feasibility_audit_addendum.md):
all three below-B0 centres are exact-infeasible (rp+ breaks sum(rp)=1 by
+3.000000000086e-07 exact-dyadic; glob_dist_r+ breaks the dist simplex by
+7.0e-07 and dist/dist_max marginal equalities at 1e-7). They are
certified objective evaluations at inadmissible points and bound nothing
about omega. The certified bound record for the VXXZ24 rung remains the
frozen 2.3715538358544617 (tighter re-evaluation of the same quantity in
this campaign: 2.3715538358544612). See the addendum for the sensitivity
reading that survives.
