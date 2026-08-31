# SUBSPACE CERTIFICATION — CHECKPOINT REPORT (T+~35 min of the 60-minute
# budget; outcome: FAILURE TO CERTIFY, reported as such per Main condition 3)

Domain attempted (Amendment A): D = { delta : A·delta = 0 (exact kernel,
27 margin rows, rank 24 proven over Q), |delta|_inf <= 1e-7 } on the
region-0 glob dist block (45 coordinates). Everything below is
computed, nothing is asserted.

## OUTCOME: the certified linear-plus-remainder route does NOT close with
## the current evaluator. Reported as FAILURE TO CERTIFY.

Reason, with the numbers:

1. NAIVE INTERVAL ENCLOSURE OVER THE ENCLOSING BOX (radius 1e-7 on the
   45 dist coords, all other coordinates points; the box CONTAINS D):
   two-sided certified endpoints
     Omega_LOW = 2.3715431556360005
     Omega_UP  = 2.3715624359361276
     WIDTH     = 1.9280300127100247e-05
   Point check at p* with the same code path: both endpoints
   2.3715518061863823, width EXACTLY 0.0 (so the box widening is real
   interval growth, not noise — the earlier point-evaluation defect is
   not present here; glob num_block widths at radius 1e-7 are
   6.656e-06 / 6.616e-06 / 6.616e-06 for region 0 and 3.614e-20 for the
   untouched regions, i.e. the box is LIVE and correctly localized).

2. THE SIGNAL WE MUST RESOLVE (first-order, float, NOT certified):
   coordinate gradient of the certified endpoint (45 central differences
   at h = 1e-9 through the same certified aggregation; 114.6 s):
     G range [-2.5865745101327775, -1.0278962125909175], |G|_1 = 79.070
   LP over D (well-scaled form, see the trap below):
     min_{delta in D} G·delta = -1.5816497000997742e-07
     max = +1.5816497000997742e-07 (symmetric, as D is symmetric)
   So the best conceivable first-order improvement inside D is
   1.58e-07 in omega.

3. THE GAP BETWEEN 1 AND 2 IS THE FAILURE: the interval enclosure is
   1.93e-05 wide where the signal is 1.58e-07 — a factor of 122 too
   loose. A certified cap needs the enclosure tighter than the signal;
   dependency loss in the entropy terms (each of the 45 coordinates
   entering H(dist) and the p_comp weights independently) dominates.
   Subdivision cannot fix it: widths scale ~linearly in the radius, so
   resolving 1e-9 would need radius ~1e-10, i.e. (1e3)^45 subboxes.
   Closing this requires a DERIVATIVE-AWARE enclosure (AD/slope
   arithmetic through the transcribed tree, mean-value or centred form),
   which is new machinery, not a bounding step — correctly outside the
   60-minute budget, so the run stops here per condition 3.

## THE INSTRUMENT THAT DOES WORK, recorded for the follow-up campaign

For ANY y in R^27 and every delta with A·delta = 0:
    G·delta = (G − A^T y)·delta  >=  −R·|| G − A^T y ||_1
so a rigorous lower bound on the linear term needs NO dual feasibility
check whatsoever (the equality constraint kills the y-term identically).
With the LP's own y: −R·||G − A^T y||_1 = −1.5816497000997758e-07,
matching the LP optimum to 1.6e-22 — and with y rounded to a 1e-6
rational grid the bound is −1.5816622065217566e-07, i.e. rounding costs
1.25e-12 and stays rigorous. Trivial y = 0 gives the useless
−7.907046373212268e-06. So: the linear term is cheaply certifiable once
a certified gradient ENCLOSURE over the box exists; the missing piece is
exactly that enclosure, plus the second-order remainder.

## AN INSTRUMENT TRAP CAUGHT IN PASSING (worth propagating)

scipy.optimize.linprog(method="highs") with bounds of magnitude 1e-7
returned a "solution" violating those bounds by 100% (two coordinates at
-2e-07 against a -1e-07 lower bound), because HiGHS's default primal
feasibility tolerance (1e-7) equals the bound magnitude. The objective
value it reported (-1.8050059225060977e-07) was correspondingly wrong —
it exceeded the true optimum -1.5816497000997742e-07 by 14%. Fix used:
scale to |u| <= 1 with delta = R·u; the scaled solve returns bound
violation exactly 0.0 and kernel residual exactly 0.0. Same shape as the
session's other instrument failures: the checking tolerance was the size
of the thing being checked.

## STATUS AFTER THIS CHECKPOINT

- Gate C's scoped negative (parent campaign) STANDS as the final result:
  no feasible improvement among the pre-registered single-block probe
  directions at +-1e-7; region_prop admits no feasible motion at all
  (exact 3-row pin, det = -3, feasible set a single point); whether a
  feasible improving direction exists in the exact 21-dimensional
  nullspace is OPEN and remains open.
- What this checkpoint ADDS, all MACHINE-VERIFIED except where marked:
  (i) the exact kernel (rank 24 / dim 21 over Q, A·V = 0 identically);
  (ii) the three structural dependencies named, with {dep0, dep1, dep2}
  proven to be a basis and the fourth relation shown to be dep1 + dep2;
  (iii) a two-sided certified enclosure over the enclosing box —
  Omega in [2.3715431556360005, 2.3715624359361276] for EVERY point of
  that box, hence for every point of D: this IS a rigorous statement,
  just far too weak to decide the 1.58e-07 question;
  (iv) the first-order bound (float, COMPUTATIONAL-EVIDENCE) that no
  point of D can improve the endpoint by more than 1.58e-07 — 12.8x
  short of the 2.0258544615181506e-06 gap to the published 2.37155181,
  which independently confirms the pre-statement's scope-honesty
  paragraph: this domain cannot reach the published value even on total
  success.


## WHAT THE POLYTOPE ROUTE BOUGHT BEYOND RIGOUR (Main, owner-verified;
## added at close of campaign)

Best conceivable gain over the true domain D versus the inscribed cubes
that the earlier parameterizations would have certified:

| domain | best conceivable gain | shortfall to published 2.37155181 |
|---|---|---|
| inscribed cube, SVD basis (r_c = 3.3631e-8) | ~1.34e-9 | ~810x |
| inscribed cube, rref basis (r_c = 4.7619e-9) | ~1.90e-10 | ~10,677x |
| **polytope D = kernel ∩ box (this route)** | **1.5816497000997742e-07** | **12.808x** |

So dropping the cube raised the POTENTIAL by roughly two orders of
magnitude as well as removing the basis artifact: the open question is
now "12.8x short" rather than "four orders short" — the difference
between a dead end and a live open question. Owner re-verified at 200
bits: width ratio 121.900x, shortfall 12.808x, HiGHS overshoot 14.12%,
subdivision cost (1e3)^45 = 1e135 boxes.

## THE LAGRANGIAN INSTRUMENT IS STRONGER THAN THE CAUTION IT ANSWERS
## (recorded at Main's instruction)

Main's condition required exact verification of DUAL FEASIBILITY, on the
grounds that an approximately-feasible y gives a non-rigorous bound.
The concern evaporates structurally: on the kernel A·delta = 0
identically, so yᵀA·delta = 0 for ANY y whatsoever; hence
G·delta = (G − Aᵀy)·delta and Hölder gives
    G·delta >= −R·|| G − Aᵀy ||_1   for ARBITRARY y,
with optimizing y merely tightening the bound. There is no dual
feasibility to check because the equality constraint annihilates the
y-term by construction. Rational-grid y (1e-6) costs 1.25e-12 and stays
rigorous. Whoever builds the gradient/slope enclosure gets a rigorous
linear bound for free — this is the reusable instrument.

## REPO-STANDARD REMEDY REGISTERED (Main broadcast to Mm3GateC and
## DelcapGateC)

Solver tolerance equalling the quantity being checked: HiGHS default
primal feasibility tolerance 1e-7 against box radius 1e-7 returned a
point violating its own bounds by 100% with a plausible objective 14.12%
off. REMEDY: non-dimensionalize (delta = R·u, |u| <= 1) rather than
tightening tolerances and hoping. Mirror-image risk named by Main: an
ILP "INFEASIBLE" that is really "infeasible within tolerance".

## FINAL STATE (Main's recording, 2026-08-30T15Z; AD build NOT approved)

Result: a scoped negative — no feasible improvement among the six
pre-registered single-block probes at ±1e-7; region_prop admits no
feasible motion at all (det A = −3, unique point (1/3,1/3,1/3),
owner-verified). Assets, not results: exact rational kernel basis; exact
rank-24 proof over Q; three named structural dependencies with the
redundancy identified (dep1 + dep2 = the second total relation);
the Lagrangian bound instrument; the rigorous-if-weak enclosure
omega ∈ [2.3715431556360005, 2.3715624359361276] valid over ALL of D;
and a measured statement of the missing machinery (certified
gradient/slope enclosure through the transcribed tree + second-order
remainder). Rationale for stopping: new multi-hour machinery whose total
success caps 12.8x short of the published value cannot produce a record
and would improve only our own enclosure — a clean complete negative
beats an extended incomplete positive.
