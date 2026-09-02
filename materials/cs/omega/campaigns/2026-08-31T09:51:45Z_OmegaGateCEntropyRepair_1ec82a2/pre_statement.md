# PRE-STATEMENT — OmegaGateCEntropyRepair

Written 2026-08-31T09:51:45Z before any repaired entropy-certificate objective
evaluation. Parent witness campaign:
`2026-08-31T09:18:39Z_OmegaGateCCandidateWitness_e9c6414` (all five direct
point differences overlap; OPEN).

## Fixed target and point order

This is a certificate repair for the same direction, not a new search. Use the
exact integer basis coefficients, direction d, center, and exact rational
coordinates already frozen in the parent campaign. Assert exact `A*d=0`,
`sum(d)=0`, positivity, and radius guards at startup. Evaluate in this fixed
order: base, then `R/16`, `R/8`, `R/4`, `R/2`, `R`, with
`R=0x1.ad7f29abcaf48p-24`. Stop immediately after saving the first strict
witness; otherwise evaluate all six points.

## Branch equation and certified entropy dual

For each level-3 Part entropy branch and each glob entropy branch, let B be the
0/1 matrix obtained by vertically stacking the three checked-in `j2m.mats`
(the columns use the branch's frozen split/shape order), let p be the exact
point distribution, and let `b=B*p`. The constrained entropy maximum is

    Hmax(b) = max { H(x) : x>=0, sum(x)=1, B*x=b }.

For any real multiplier y on any exact independent row subset of B, define
`theta=B_ind^T*y` and

    U(y) = log(sum_i exp(theta_i)) - y^T b_ind.

By nonnegativity of KL divergence, `Hmax(b) <= U(y)` for arbitrary y; solver
success is never a premise. Evaluate U with Arb, a fixed max-shifted log-sum-exp,
and y interpreted as exact binary64 dyadics. The lower bound is H(p), because p
itself is primal feasible by construction. If a stored `dist_max` happens to
pass exact interval nonnegativity, normalization, and all B-marginal equality
tests, its entropy lower bound MAY replace H(p); otherwise it is not consumed.
Thus the repaired penalty interval is

    scale * [max(0,Hprimal_low-H(p)_high),
             max(0,U(y)_high-H(p)_low)],

where `scale=region_prop` for glob and
`scale=part_frac*region_prop` for Part. No `hm+2*eps` residual enclosure from
the parent is consumed.

The old max-residual update used `min(current_lower,residual_lower)` from a
zero seed, which conservatively forced the lower endpoint to zero and created
the ~2.03e-6 raw width. It remains frozen as history but is superseded here.

## Multiplier generation and solver-independent acceptance

Generate candidate multipliers only on an exact linearly independent row set
of B (row indices selected deterministically by exact SymPy RREF of `B.T`).
Use float L-BFGS solely to minimize the convex dual and return a candidate y.
The certificate consumes only the Arb value U(y).

Separately compute the Arb KKT residual `B*q(y)-b`, where
`q_i=exp(theta_i)/sum exp(theta)`. A branch multiplier is ACCEPTED only if the
maximum absolute residual upper endpoint is <= `2^-32`. Every entropy branch
consumed by a strict witness MUST be accepted; any rejection forces the
campaign decision to OPEN, even though its arbitrary-y upper bound remains
valid. Record solver status/iterations but never promote them above
COMPUTATIONAL-EVIDENCE.

Counterfactual plant, before any real decision: for the base glob region-0
branch add exact `1/8` to the first independent multiplier. The planted
candidate MUST fail the `2^-32` KKT acceptance threshold and have a strictly
larger certified dual upper bound than the unplanted candidate. Failure to
catch the plant aborts the campaign.

## Full objective and decision

Recompute all seven R components with interval
`max(0,min(candidate_0,candidate_1,candidate_2))`, the interval minimum M, and
Arb c/ceq/line defects. Compare each point to base with the same registered
shared-target cross product

    ((T-R_new)M_base-(T-R_base)M_new)/(M_new M_base)

plus repaired defect difference. Exact denominator positivity is mandatory.

- WITNESS only if every used entropy branch is accepted and
  `upper(full_new_minus_base)<0` strictly; freeze that first rung and stop.
- OPEN if intervals overlap, if every strict test fails, or if any required
  branch multiplier is rejected. Strictly positive lower at one rung rejects
  only that rung.
- No midpoint LP sign, optimizer success flag, or paper/record comparison is
  consumed.

## Rule 7 scope sentence

This campaign evaluates exactly the base and the five already-registered
points `p_star+eps*d` for `eps=R/{16,8,4,2,1}` on one exact integer direction,
changing only the 45 region-0 glob `dist[0]` coordinates of the frozen VXXZ24
`K100_2.37155181` vector under the transcribed 3-region single-p_comp program
at max_level=3, q=5; it changes only the entropy certificate, not the points or
objective. It does not search another direction, negative/intermediate steps,
a box, another region/parameter block, multi-block/cross-block classes,
another q/max-level regime, unpublished record parameters, or another
construction. A witness would refute local optimality only on this named ray;
no paper or record claim is in scope.

## Pending command and cap

Run only after the shared heavy slot is released:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 nice -n 10 \
    ./.venv/bin/python omega/src/gate_c_entropy_repair.py \
      omega/campaigns/2026-08-31T09:51:45Z_OmegaGateCEntropyRepair_1ec82a2/results.json

Cap: base plus five fixed points, at most 2000 L-BFGS iterations per entropy
branch, 900 wall seconds, no adaptive points or threshold changes.
