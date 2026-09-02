# PRE-STATEMENT — OmegaGateCCandidateWitness

Written 2026-08-31T09:18:39Z before any candidate-witness objective
evaluation. Parent campaign:
`2026-08-31T07:57:00Z_OmegaGateCDiag_4859327` (V1 PASS; full-box
FAILURE TO CERTIFY / OPEN). This campaign asks only whether its numerical LP
candidate can be converted into one exact feasible improving witness.

## Exact direction and startup anchors

Use the authoritative exact integer basis V from
`kernel_basis_V_exact_rational.json` (45x21; entries -1/0/1) and exact
coefficient vector

    u = [1,-1,0,-1,1,0,1,0,0,1,-1,-1,-1,1,0,-1,0,0,-1,0,1].

The direction is computed at startup by exact integer multiplication d=V*u
and MUST equal

    d = [0,1,1,-1,-1,-1,1,0,0,-1,-1,1,0,1,1,-1,0,0,1,-1,0,
         -1,1,0,-1,1,0,0,1,-1,1,-1,-1,1,0,1,-1,0,0,1,-1,0,-1,1,0].

Startup MUST assert with integer arithmetic that A*d=0 in all 27 rows,
sum(d)=0, rank(V)=21 (CITED-DEPENDENCY from the frozen exact-basis
certificate), and every d_i is in {-1,0,1}. The convention is columns of V
are kernel coordinates, d=V*u; A has margin constraints as rows and the 45
region-0 glob dist[0] probabilities in the frozen shape order as columns.

## Predeclared exact dyadic step ladder

Let R be the exact binary64 dyadic `0x1.ad7f29abcaf48p-24` (display 1e-7).
Evaluate, in this fixed order,

    eps in [R/16, R/8, R/4, R/2, R].

Every moved coordinate is the exact Arb number
`IC.const(p_i_star) + eps*d_i`; no binary64 addition may construct the moved
point. The ladder is fixed before results. The midpoint first-variation
candidate at R is -1.577332901516852e-7 (COMPUTATIONAL-EVIDENCE), so its
linear prediction at R/16 is -9.858330634480325e-9; this is the reason the
smallest rung is still large enough to separate from the inherited
approximately 1.94e-11 endpoint defect scale while reducing quadratic terms.
This reason is a hypothesis, not a premise.

Exact startup guards MUST establish for every rung: all 45 moved probabilities
are strictly positive; sum(new)-sum(base)=eps*sum(d)=0 exactly; A*(eps*d)=0
exactly; max_i |eps*d_i| <= R; every non-target parameter is byte-identical to
the frozen center. At R, the exact positivity ceiling among negative
direction coordinates is 1.6166317499672092e-5, so R is only
0.0061857006088138706 of that ceiling (both INFERENCE from exact dyadics,
recomputed and asserted at startup).

## Full direct objective and common-expression comparison

For base and each moved point, rebuild the full transcribed 3-region,
single-p_comp, max_level=3, q=5 workspace with exact point Ivals. Recompute all
level-3 Part, level-2, and three glob R components. For every R component use
the interval `max(0,min(candidate_0,candidate_1,candidate_2))`, with lower
endpoint `max(0,min lows)` and upper endpoint `max(0,min highs)`; similarly
use the interval minimum of all M branches. Thus no endpoint argmin is pinned.
Use the primary-source Part/glob penalty formula `prop*(Hmax-H+2*eps)` and
compute c, ceq, and line-slack defect terms from their Arb intervals, never
from float residuals.

Compute the raw difference with the shared-target cross product

    raw_new - raw_base =
      ((T-R_new)*M_base - (T-R_base)*M_new)/(M_new*M_base),
      T = 4*log(7),

then add the interval difference of the full absorbed-defect terms. This is
the consumed common-subexpression comparison; independently rounded endpoint
subtraction is report-only. M lower endpoints and all denominator intervals
MUST be strictly positive.

The frozen base anchor MUST enclose raw 2.3715538358350803 and certified
2.3715538358544617. Failure of any anchor/guard aborts without a witness.

## Decision rule fixed before compute

- **WITNESS / scoped counterexample to local optimality:** for any registered
  eps, `upper(full_new_minus_base) < 0` strictly. Freeze the first such rung
  and all attempted earlier rungs. This establishes only that the frozen
  center is not locally optimal within the named one-dimensional exact
  kernel ray.
- **OPEN:** no rung has strict negative upper endpoint. Overlap is OPEN, not
  no-improvement. A strict positive lower endpoint at a rung only rejects
  that one registered step.
- No float LP feasibility or midpoint sign is consumed by the decision.

## Rule 7 scope sentence

This campaign sweeps exactly five points `p_star + eps*d` for the predeclared
eps ladder above on one exact integer direction d in the 21-dimensional
kernel of the frozen 27x45 margin matrix, changing only the 45 region-0 glob
`dist[0]` coordinates of the VXXZ24 `K100_2.37155181` released vector under
the transcribed 3-region single-p_comp program at max_level=3, q=5. It does
not search any other kernel direction, negative steps, intermediate step,
box, region-1/2 block, part/split/lambda/region_prop or other parameter block,
multi-block/cross-block class, q/max-level regime, unpublished record
parameters, or construction outside this enclosed laser-method program. A
witness, if obtained, refutes local optimality only at this released center
and cannot establish or improve any published/record exponent; paper and
record claims are out of scope.

## Pending command and resource cap

Run only after the shared heavy slot is handed over:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 nice -n 10 \
    ./.venv/bin/python omega/src/gate_c_candidate_witness.py \
      omega/campaigns/2026-08-31T09:18:39Z_OmegaGateCCandidateWitness_e9c6414/results.json

Hard cap: six full direct point evaluations (one base plus five rungs), no
adaptive steps, no retries with changed arithmetic, and 900 wall seconds.
