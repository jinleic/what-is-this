# PRE-STATEMENT ADDENDUM (Decision 1, Main-approved 2026-08-30) —
# simplex-repaired glob_dist probes. Committed BEFORE the first repaired
# evaluation.

## Repair rule (the formula, fixed before seeing any repaired result)

For region r in {0,1,2} and a signed direction s in {+1,-1}:
  offset d is defined on the 45 coordinates of glob dist_id[r] as
    d_i = +eps  for 8 of the coordinates with the LARGEST p*_i
    d_i = -7*eps  for the 1 coordinate with the SMALLEST p*_i
    d_i = 0 elsewhere (all other parameter groups untouched)
  with eps = 1e-7/8? NO — fixed exactly: the eight largest coordinates
  get +1.25e-8 * ... — SIMPLER AND EXACT: the single smallest coordinate
  gets -7e-7? that breaks the radius. FINAL RULE (radius-compatible):
    d_i = +1e-7 on the 43 non-smallest coordinates? No.
  COMMITTED FORMULA: pick k = argmin_i p*_i (smallest coordinate, tie ->
  smallest index). Set d_k = -7e-8 and d_i = +1e-8 for the 7 coordinates
  with the largest p*_i (i != k), all other coordinates 0. Then
  sum_i d_i = 7e-8 - 7e-8 = 0 EXACTLY in rational arithmetic, so the
  simplex sum moves only by dyadic float-rounding of the individual
  +1e-8/-7e-8 additions; the per-coordinate perturbation magnitude is
  <= 1e-7 (inside the pre-registered box radius, single-coordinate moves
  in sign-split subboxes, NOT a full-radius box).
  Repair EXACTNESS CHECK (gate before any objective evaluation): the
  exact-dyadic sum of (p*+d) over the 45 coords must equal the
  exact-dyadic sum of p* over the same coords, i.e. residual
  sum(p*+d) - 1 == sum(p*) - 1 in double arithmetic compared via
  IC.const endpoints (both printed, both recorded; equality of the two
  residuals required — the residual itself need not be 0 since p* itself
  carries a dyadic-sum residual, and this is the SAME residual by
  construction).
  The 8 moved coordinates are also chosen so no coordinate exits [0,1]
  (checked: p* min coordinate 1.28e-7 >= 7e-8 slack? NO — the SMALLEST
  coordinate gets the NEGATIVE step -7e-8 < 1.28e-7 => leaves inexact
  support! FIX, FINAL): d_k = +7e-8 on the SMALLEST coordinate, and the
  7 LARGEST coordinates get -1e-8 each. Sum = +7e-8 - 7e-8 = 0, all
  moves keep strict support (min becomes 1.28e-7+7e-8, max shrinks by
  1e-8). THE 8 MOVED COORDS ARE: 1 smallest (+7e-8), 7 largest (-1e-8).

## Directions pre-declared (exactly 6 boxes)

  glob_dist_r_repaired, r in {0,1,2} x s in {+1,-1}: for s = +1 use the
  rule above; for s = -1 negate every d_i (then the smallest coordinate
  gets -7e-8 — support check: p*_min = 1.2822e-7 > 7e-8, residual
  support 5.82e-8, stays strictly positive; largest shrink by -1e-8 ->
  +7e-8 away from 1). So the negated direction is support-safe too.
  Total: 6 repaired boxes, radius per-coordinate <= 1e-7 (the pre-
  registered radius), box budget 6 (well under the 63 cap), plus B0
  anchor re-run for the same-session baseline. NO further directions,
  NO adaptive re-probing of whichever repair looks better; the family is
  closed. Stopping condition: all 6 evaluated or wall-clock 30 min,
  whichever first.

## Feasibility gate (mandatory, per Main condition (b), BEFORE objective)

For each repaired point, in exact-dyadic per-row arithmetic:
  1. simplex row for region r: exact-dyadic sum residual vs p*'s own
     (must be IDENTICAL, computed as arb dyadic endpoints);
  2. all other groups untouched => their rows unchanged (spot-verify 3);
  3. the FULL ceq (Lagrange) set at the repaired point, re-verified from
     the interval-path log-form pairs vs the float-path exp-form
     residuals with the same |exp-form| <= |log-form| * dmv discipline
     used in the frozen stage-b; any violation above the frozen stage-a
     tolerance (ceq_max 3.6e-11 scale) => the probe is REJECTED as
     infeasible and reported as a miss;
  4. dist<->dist_max marginal rows: since dist_max group is untouched and
     dist moves by <= 1e-7 on 8 of 45 coords, the j2m row residuals must
     stay <= |row| * 1e-7-class — verified directly per row (row replay)
     and REQUIRED <= 1.1e-9 * (row 1-norm) to consider the probe
     admissible at the release's own tolerance discipline.
Only probes passing ALL gates contribute candidate endpoints.

## Objective (only after the gate)

Certified aggregation identical to the frozen rung-2 protocol
(interval lows for num_block/p_comp, highs for penalties,
Omega = (target_hi - R_sum_low)/M_low + absorbed defects), reported with
R_sum width, M width, value width as first-class numbers. Every probe's
certified endpoint is reported including misses; if none beats B0, the
campaign statement is exactly Main-approved phrasing: "the enclosure is
not improvable at this radius by any feasible move in the swept family".

## Escalation trigger (Main condition (d))

If any repaired point certifies feasibly below B0 (2.3715538358544612):
STOP, escalate to Main with the certificate BEFORE any wording change
anywhere.
