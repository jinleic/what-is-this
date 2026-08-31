# Decision-1 repaired probes — RESULT (run 2026-08-30T12:4xZ, after
# pre_statement_addendum_decision1.md was committed)

All 6 pre-registered repaired probes ran the full feasibility gate.
Raw JSON: `repaired_probes.json` (this directory). Script as run:
`src/gate_c_repaired_probes.py` (copy sealed below).

## Outcome: ALL SIX REJECTED at the feasibility gate — none reached
## objective evaluation.

Per-probe (marginal = worst |A d - A dmax| over the region's 3 margins;
ceq = max |Lagrange residual| at the repaired point):

| probe | simplex residual base vs repaired | worst marginal | ceq_max | verdict |
|---|---|---|---|---|
| glob_dist_r0_repaired+ | −2.1449508835758024e-13 (IDENTICAL) | 7.000036723187053e-08 | 3.5953524879506205e-11 | REJECT (marginal) |
| glob_dist_r0_repaired− | SAME | 7.000000000000000e-08 | 3.5953524879506205e-11 | REJECT (marginal) |
| glob_dist_r1_repaired+ | +3.597122599785507e-14 (IDENTICAL) | 7.000019317665584e-08 | 3.5953524879506205e-11 | REJECT (marginal) |
| glob_dist_r1_repaired− | SAME | 7.000002211210443e-08 | 3.5953524879506205e-11 | REJECT (marginal) |
| glob_dist_r2_repaired+ | +2.9620750296999176e-13 (IDENTICAL) | 6.999999999999999e-08 | 3.5953524879506205e-11 | REJECT (marginal) |
| glob_dist_r2_repaired− | SAME | 7.000097472856237e-08 | 3.5953524879506205e-11 | REJECT (marginal) |

The repair rule worked exactly as designed on the simplex row: the
repaired residual EQUALS p*'s own residual bitwise in every case (the
dyadic perturbations cancel by construction). But the dist<->dist_max
"share-marginals" rows move by exactly ±7e-8 — the size of the summed
perturbation — because dist_max is held fixed and the marginals are
LINEAR in the moved coordinates. There is no degree of freedom that
moves dist alone: every marginal-sharing row pins one linear functional
of dist to dist_max's, so ANY strictly-inside dist move at the 1e-7
scale is infeasible unless dist_max moves with it in the exact opposite
pattern, which the repair rule does not allow (dist_max was declared
untouched).

## Campaign statement (Main-approved phrasing, satisfied structurally)

"Within the pre-registered box, at this radius, by any feasible move in
the swept family, the enclosure is NOT improvable." Now stronger than
requested: the gate killed the family at feasibility, so no objective
comparison was even needed. Combined with the rp no-go (asymmetry rows
freeze rp), the certified local-optimality picture at ±1e-7 is:

- region_prop: frozen by 3 exact rows — no feasible ±1e-7 move exists;
- glob dist (single-region moves): every move breaks either the simplex
  row (unrepaired) or the 27 marginal-sharing rows (repaired); a
  coordinated dist+dist_max move would need the repair rule generalized
  to BOTH blocks (regenerating matching marginals), which was NOT
  pre-registered and is a new campaign;
- omega/single_mat_size/single coordinates: break the Schoenage line at
  the 1e-7 scale (stage-a screen), unless R moves to compensate — same
  coordinated-move class.

The honest gate-C conclusion is the NEGATIVE, in exactly the shape Main
approved: the released point is certified-locally-optimal at the ±1e-7
scale within every feasible single-block move family tried; the 2.026e-6
gap to the published rounding cannot be closed by any in-family
feasible move, and the infeasible-direction sensitivity (−8.15e-7 max)
shows where a coordinated dist+dist_max re-optimization campaign would
have to look. Machine labels: feasibility gates MACHINE-VERIFIED
(exact-dyadic per-row); no certified objectives were produced in this
run (all probes rejected pre-objective per the mandatory gate order).

## POST-RECORD SHARPENING (Main, 2026-08-30T12:5xZ — adopted into README)

Main's sharpening accepted and propagated to the README: the gate-C
verdict is LOCAL OPTIMALITY BY CONSTRAINT RIGIDITY, not by comparison.
No feasible alternative was evaluated-and-found-worse; there ARE none in
the swept single-block families (rp pinned by 3 rows; glob dist jointly
frozen by simplex-vs-marginals — marginals linear in dist, dist_max held
fixed; omega/sms break the Schoenhage line without compensating R). The
released point is the UNIQUE feasible point of each swept family inside
the box; the objective was never binding because the families died at
feasibility. This is a statement about the PARAMETERIZATION, not the
landscape. Explicitly NOT established: "no better nearby point exists" —
coordinated multi-block moves (dist+dist_max in opposite linear
patterns; omega/sms with compensating R) are entirely unexamined and are
the live class. Rung-vs-record distinction kept: VXXZ24 rung at
2.37155181; record 2.371177 unpublished as of 2026-08-30, never searched.
