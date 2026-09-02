# PRE-STATEMENT — gate C stage 5: analytic branch-pinning of the R Glob[0]
# minimizer identity (Lemma-1 residual pin) over the dist[0] subdomain D.
# Agent OmegaNext, written 2026-09-01T~09:30Z, committed to git BEFORE the
# first compute of this campaign (Main strict-order rule: no parameter-
# dependent probe before this commit; only tool/help inspection and reads of
# frozen artifacts happened before this file).

## 0. Selecting the technology

Main's directive 2026-09-01 names two registered routes:
(a) certified-gradient feasibility (a tighter certified gradient enclosure),
(b) analytic branch-pinning (deriving and interval-certifying branch-order /
    minimizer inequalities over explicit subdomains, then reducing
    dimension only where those identities are proved).

Uniform branch-and-bound over D is FROZEN INFEASIBLE (parent campaign
2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c: k* = 36 of 45 at 90% L1
width, ~2^223 subboxes; the obstruction explicitly does not bind "an
analytic branch-pinning argument"). This campaign selects (b) — analytic
branch-pinning — and applies it to the single branch that carries the
dominant hull cost: `R_glob[0]`, whose candidate-interval span of
3.2123059392596964e-05 is ~203x the first-order signal
1.5816497000997742e-07 (stage1_decomposition.json,
`branch_straddle_1_3.R_glob[0].max_candidate_span`), whereas every other
branch's span is <= 2.06e-6. NOT selected: (a) certified-gradient
feasibility — the parent's Gaussian/Hessian route is the one whose
successor (slope-bundle + ANY-y width/signal 73.86) just measured SPREAD,
and no concrete certified-gradient technology with concentrates-width
structure was registered on this target; branch-pinning has a concrete,
checkable identity target (below). No re-run of B&B occurs.

## 1. The analytic identity (derived from recorded frozen evidence only)

Written with the transcribed program's own symbols (src/vxxz24_float.py,
GlobalStage): for region r the R aggregation candidates are
  cand_r,tt = nb_r,tt - (pc_r if tt == r else pen_r),
where nb_r,tt is the entropy of the tt-th margin of dist[.], pc_r is the
r-th p_comp term and pen_r = prop_r * (H(dm_r) - H(d_r) + 2*eps_r); the
Lemma-1 residual eps_r = max_t (ln(dm_r,t) - (lam_sum_r - 1 +
sum_d lam_margin_r,d[shp_t,d])).pos() is built ONLY from dm (dist_max
block), lam_sum, lam_margin — three parameter blocks that rule 16 holds at
p* (POINTS, never boxed). The frozen CandidateWitness base record exposes
this arithmetically: at r = 0 (radius 0) the R_glob[0] candidate interval
for cand_0,tt equals for tt = 1, 2 exactly:
  cand_0,1 - cand_0,0 = prop_0*(H(dm0) - H(d0) + 2*eps_0) + pc_0
                       - (nb-gap terms that vanish at the released point
                          because the three margins of the released dist0
                          are equal to display tolerance),
i.e. the straddle 1.433472463174823e-06 at r = 0 decomposes as
  2*prop_0*eps_0 + [pen entropy-mismatch + pc terms].
Frozen values: glob_r0 eps_hi = 2.1502086942121845e-06,
2*prop_0*eps_0 = 1.433472462808123e-06, matching the candidate span to
3.7e-16 (float render granularity). THE IDENTITY IS REGISTERED AS THE
CAMPAIGN'S TARGET; it is verified IN-RUN as control C1 (below), not
assumed.

Over D (which moves ONLY the 45 dist[0] coordinates), eps_0 is CONSTANT:
both of its argument families are p*-frozen parameter blocks. So the
boxed growth of the R_glob[0] candidate span from 1.433e-6 (r = 0) to
3.212e-5 (r = 1e-7) is entirely the movement of
  nb_0,tt (entropy of the t-th J2M margin of the moving dist0),
  pc_0   (the p_comp[0] functional of dist0),
over the box — while the epS term stays frozen at 2.1502086942121845e-06
(point-exact via IC.const). Branch-pinning therefore reduces to ONE
question: can the hulled minimizer identity of R_glob[0] over D be
pinned by an inequality on the FIXED side, so that the hull over three
candidates collapses to one?

## 2. AMENDMENT 1 (committed PRE-COMPUTE, before any runner execution —
BnBGate-precedent amend discipline): the minimizer identity direction

The first committed version of this file (754b66d) registered the pinned
candidate as cand_0,0 (the p_comp candidate). Deeper reads of the frozen
records (still NO compute done) contradict that: stage_b WITH_lemma1
records R_glob[0] = 0.49684932826618977 — economically the PENALTY
candidate's low (CandW cand1 = 0.49684932826618944, cand2 = ...5137) —
while stage_b WITHOUT_lemma1 records R_glob[0] = 0.49685076173807374,
exactly the p_comp candidate (CandW cand0). The Lemma-1 charge
(2*prop_0*eps_0 = 1.4334724628e-06) is what FLIPS the winner from
cand_0,0 to the penalty candidate at p*. Boxing dist[0] (radius 1e-7,
BnBGate stage-1.2) moves all three candidate lows DOWN: cand0 by
1.538e-5, cand1/cand2 by 1.142e-5 — the p_comp candidate falls FASTER
(its pc term inherits num-block motion), and at the box the argmin flips
back to tt=0 while all three interval candidates OVERLAP
(possible_minimizers {0,1,2}, span 3.2123e-5). THE CORRECTED PIN
QUESTION, registered now: certify over the whole box D the ORDER
INEQUALITIES between the INTERVAL candidates
  G1: cand_0,1.lo - cand_0,0.hi (lower-bound gap cand1 vs cand0),
  G2: cand_0,2.lo - cand_0,0.hi,
each candidate interval being the OUTWARD interval of nb_0,tt - rhs_tt
over the box (nb interval-heaved, rhs interval-heaved). PASS requires
BOTH gaps strictly positive at their outward-rounded endpoints (cand_0,0
uniquely the possible minimizer over D; the R_glob[0] hull collapses to
{0}; the 3.212e-5 hull width collapses to single-candidate width). C1 is
retained but re-aimed: at radius 0 the frozen data gives cand1.lo -
cand0.lo = -1.433471884e-6 = -2*prop_0*eps_0 (match 5.8e-13); the
identity to certify is
  cand_0,tt - cand_0,0  contains/is  -(2*prop_0*eps_0) + delta_nb,
with delta_nb the interval-certified nb-gap term; C1 asserts the radius-0
reproduction of the frozen candidates within 1e-12. Everything else
  (domain Sec 3, machinery Sec 4, controls Sec 5, backend/tolerances/cap
Sec 6, adjudication Sec 7 adapted to G1/G2, scope Sec 8, not-swept
Sec 9) stands as originally registered, with the PIN FAILURE margin
redefined as the numeric value of min(G1, G2) should it fail. No other
telemetry changes; NO compute has occurred; this amendment is committed
before the first runner execution (manifest attempt-log attests order).

## 3. Domain (unchanged, rule 16)

D = { delta in R^45 : (margin matrix A)·delta = 0 exactly,
      |delta_i|_inf <= 1e-7 } around the region-0 glob dist[0] block of
the VXXZ24 K100_2.37155181 released vector (MAT sha256 of the params
vector recorded in the frozen rung-2 campaign: df75ae3a-prefix record;
assert in-run). A = frozen 27x45 0/1 margin matrix (rank 24 over Q; the
exact kernel basis is not needed for THIS pin because the pin is proved
on the whole axis-aligned box D, strictly larger than the kernel
parametrization). The box endpoints are EXACT ARB:
IC.const(center_i) +/- IC.const(radius), center/radius binary64, radius
1e-7 (RADIUS constant identical to the parent). max_level = 3, q = 5.
Every other parameter block (region_prop, dist_max[0..2], lam blocks of
all parts and glob, part split_dist blocks, part region_props, omega,
single_mat_size, num_retain_*) is held at p* — point-evaluated, never
boxed, per rule 16. The certified subdomain theorem (if PASS) is that
the R_glob[0] hull collapses to the UNIQUE possible minimizer over ALL
of D (per Amendment 1's corrected order-inequality target G1/G2).

## 4. Machinery (pre-declared build; no existing instrument does this)
- The pin pass reuses the parent's certified interval tree (interval_core
  v11, `ce70f95922f3d5098e892e1b7047570ebb5267826b51f2776c91b56d9934adb5`,
  hash-checked at startup; vxxz24_float
  `31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890`) with
  SlopeG DISABLED — this is a pure interval run, no gradient machinery.
- P0 (nb-gap identity): the candidate_0,tt entries at tt = 1, 2 differ
  from candidate_0,0 by (nb_0,tt - nb_0,0) + (pc_0 - pen_0); that nb-gap
  vanishes only to float tolerance at p* (frozen stage-b detail shows
  6e-13-level identity), a fact that must become an INTERVAL statement:
  nb_0,tt - nb_0,0 over the box is certified directly (both terms are
  entropies of margin projections of the SAME moving dist0, and their
  box-enclosure difference is interval-adjacent). If the certified
  nb-gap width makes LB <= 0, the pin FAILS at exactly this term.
- P1 (frozen-argument invariance): assert in-run that dm0's and
  lam(isum/margin)'s IntervalTree leaves carry point intervals
  (width 0), i.e. that no dist-coupled term contaminates eps_0. The
  engine's evaluate path (Workspace.evaluate with dist0 leaves
  interval-boxed and everything else IC.const points) is the SAME code
  path the parent's stage (b) used for the frozen point run, so
  point-invariance is a leaf-width assertion, not a new instrument.
- P2 (candidate order at the centre): assert at the CENTRE of the box
  (radius-0 point pass, the frozen CandidateWitness base point) that
  cand_0,0 < cand_0,1 and cand_0,0 < cand_0,2 hold as FLOAT statements
  and the lb/ub margins match the frozen record ([0.4968507617380737,
  0.4968507617380738] cand_0 vs [0.49684932826618944,
  0.4968507617386526]; already frozen discipline).
- P3 (the pin inequality over the box): interval-certify, with outward
  rounding only (buffered `low()/up()` per the v11 core), the four
  interval quantities { nb_0,0, nb_0,1, nb_0,2, pc_0 } over the FULL
  dist0 1e-7 box, and a derivative-free bound on their BOX-EXCHANGE:
  nb-gap width and pc width both scale with the width of the J2M
  projected marginals, so a DIRECT interval pass over the box (no
  subdivision) produces:
    LB_box := EPS_lo - pen_0,hi_box + nb_0,0,hi_box - nb_0,1,lo_box
    (with pen_0,hi_box = prop*(H(dm0) - H(d0,lo_box) + 2*eps_0) and the
    analogous bt = 2 term).
  The v11 core's R3 containment checks (per-block containment against
  stage-(a) float) run at the box as they do at the point.
- P4 (the subdomain closure): if LB_box > 0, the certified subdomain
  theorem registers: "for all dist0 in D (1e-7, kernel A*delta == 0
  implied), the R_glob[0] branch contributes min(cand_0,tt) = cand_0,0
  uniquely; the hulled R_glob[0] slot of the parent's LP/J1 any-y signal
  collapses from 3 to 1 candidate and the hull-cost term
  3.2123059392596964e-05 vanishes in favour of single-candidate width."
  This is a LOCAL statement about the branch structure (NOT a pass on
  the 21-dim kernel question). The pre-declared end value added: the
  R_glob[0] branch's hull-overhead 3.212e-5 vanishes into a pinned
  candidate, and the natural successor's ANY-y width also loses the
  (prop-scaled) branch-straddle share that stage-1.3 attributed to
  R_glob[0].
- P5 (quarantine note, NOT consumption): the parent's 12.808x guardrail
  means no possible outcome of this pin makes the 21-dim kernel
  decision (the LP bound caps total attainable improvement at 12.8x
  BELOW the gap to any record). The pin is a STRUCTURAL increment: it
  either validates or kills one named hull-cost term of the currently
  frozen quantified obstruction. Nothing here touches the published
  2.37155181 or the record 2.371177, neither is reachable by any
  outcome of this campaign, and no record claim may appear anywhere.

## 5. Controls (registered obstacle-rejection gates)

- C1 (identity check, plant-accept): the P0 identity at the base point
  (radius-0) must reproduce, in v11 interval arithmetic, the frozen
  CandidateWitness base record within Machine tolerance — large margin
  (>= 2*prop_0*eps_0 interval's low must be <= the observed span's
  interval upper, and the observed interval must have width order
  1e-16 relative to itself). A C1 FAIL means the recorded identity
  mis-carries the mechanism: STOP before P3.
- C2 (structural perturbation rejection): bool-disturbed oracle. Shift
  EVERY lam_margin entry of glob region 0 by +delta_eps_user = +1e-3
  (a sign-flip-scale perturbation of the eps dual terms, broadcasting
  the SAME +1e-3 to all 9+1 lam entries) and re-certify: the pin
  boundary LB must CROSS SIGN (from possibly-positive to strictly
  negative), i.e. the must-FAIL side of the C2 gate. If the pin still
  reports LB > 0 under this drift, the pin is insensitive to the
  quantity it claims to control and is a false positive (REJECT
  verdict, FROZEN-NEGATIVE on the pin claim).
- C3 (sign-flip control, false-plant rejection): evaluate the pin at
  prop_0 located at region PROP INDEX 1 (a permutation-adjacent index
  0-free query) and confirm the machinery returns the corresPONDING
  region's dual lookup (glob_r1 eps = 1.143325024764798e-06 from the
  frozen with-lemma1 record) — not the region-0 value; i.e. the pin
  machinery must reject a wrong-index plant. A C3 FAIL means the
  machinery is reading the wrong lam blocks and every PASS is invalid.
- GATE Z (engine equivalence): re-derive the frozen point aggregation
  values (R_detail and Om_raw at radius 0, with and without lemma1)
  within the pre-registered replay-noise tolerance of the parent BnBGate
  campaign (max abs delta <= 1e-12 per R branch, <= 1e-6 on Omega). This
  is the same GATE-Z tolerance the parent's stage-1 measured (1.2e-14
  drift); it is re-registered here unchanged.
- PLANT-PASS (known-true object): before P3, the C1 pass at radius 0 IS
  the known-true certificate acceptance (frozen-matching interval
  bracket). The C2 drift and C3 wrong-index query are the known-wrong
  rejections. Both rejections are MANDATORY for a valid PASS outcome.

## 6. Arithmetic backend, tolerances, cap, kill logic

- Backend: python-flint arb, v11 endpoint-pair Ival semantics, MID = 300
  bits, OUTB = 64 (unchanged from the parent; every decision inequality
  uses the buffered `.low()/.up()` outward endpoints).
- Tolerances: GATE-Z: R-branch abs diff <= 1e-12, Omega abs diff <=
  1e-6. C1: span interval must satisfy |observed - predicted| <= 1e-12
  (i.e. the C1 identity is float-robust over the 1e-16 float-render
  granularity). C2: LB's cert interval must contain 0 or sit negative
  (LB_hi < 0) under the drift; the PASS boundary requires LB_lo > 0
  (strict positivity of the outward-rounded low). C3: the C3-extracted
  lambda-must trip the wrong-index check; if the machinery returns
  glob_r0's eps for a glob_r1 query, FAIL.
- Compute cap: 4 hours wall (checkpoint + write-first protocol per
  repo discipline; per-attempt stderr files under DISTINCT filenames per
  the 2026-09-01 freeze hygiene rule), single process, threads = 1, nice
  -n 10. Kill: any C1/C2/C3/GATE-Z failure stops the campaign with the
  failing control named verbatim; the outcome is exactly CERTIFIED (all
  controls + LB > 0 over D), FALSIFIED (any control fails or LB's
  interval brackets 0 with an identified term), or OPEN (budget).
- No subdivision, no Monte Carlo, no float-precision games (exact-Arb
  box endpoints), no float SVD basis (the pin needs no kernel basis at
  all).

## 7. Outcome adjudication (fixed now)

- PASS (CERTIFIED): C1-C3 + GATE-Z all pass AND
  LB_box = EPS_lo - pen_0,HI + nb_0,0,hi - nb_0,tt,lo > 0 outward for
  tt in {1, 2}. Deliverable: the subdomain theorem
  "R_glob[0]'s minimizer over the full 1e-7 dist0 box D is pinned to
  the tt = 0 (p_comp) candidate", plus the recorded numeric margin
  (LB_box vs the hull-cost 3.212e-5 it removes). This is certified
  LOCAL progress; the global 21-dimensional kernel question is
  UNCHANGED AND OPEN, and the 12.808x guardrail keeps every record
  unreachable via any outcome.
- FAIL (FALSIFIED): any control fails OR LB_box <= 0. Deliverable: the
  exact failing term and its numeric margin (the interval bracket of
  the specific inequality that cannot be pinned), registered as the
  successor of the BnB obstruction — the analytic route to THIS hull
  term is now also measured-false, frozen as quantified negative
  evidence for future route choice.
- OPEN: budget exhausted at any point; records the exact step reached
  and the next needed control. No verdict is inferred from a partial.

## 8. Rule-7 sentence (the scope of every statement here)

Registered swept set (if PASS): the R_glob[0] branch minimizer
inequality over D = {delta in R^45 : A·delta = 0 exactly,
|delta|_inf <= 1e-7} at the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector under the transcribed 3-region
single-p_comp program at max_level = 3, q = 5, with all other parameter
blocks (regions 1/2 glob dist and dist_max, region_prop, all lam blocks,
all part-level split_dist blocks, omega, single_mat_size, num_retain)
at p* POINT values and the interval engine v11 (ce70f959...). This is a
branch-structure theorem about ONE hull-cost term of the FROZEN BnB
obstruction — NOT a decision on the 21-dimensional kernel question
(which remains UNCHANGED and OPEN), not a statement about the published
2.37155181, not about the record 2.371177 (parameters unpublished,
never searched), not about any other rung, any other max_level/q
regime, multi-block or cross-block direction classes, or any
construction differing from the enclosed laser-method program
(including the 2026 asymptotic-rank/centroid-based improvements). Even
total success in D caps 12.808x short of the published 2.37155181, so
this route can NEVER produce a record and no record claim may appear
anywhere.

## 9. Explicit not-swept list

- All other branches of R (only R_glob[0]'s identity is pinned; the six
  other branch hulls are untouched).
- Regions 1 and 2 glob blocks, all part-level blocks, all lam/district
  blocks beyond what C3's instrument reads (C3 is an extra QUERY, not a
  sweep), the omega x single_mat_size Schoenage 2-coordinate class.
- Any subdivision of D, any radius other than 1e-7, any direction
  class outside the axis-aligned D box.
- The published 2.37155181 and record 2.371177: untouched; nothing here
  bears on the validity of the VXXZ24 paper.
- No reruns of the parent's B&B, no LP, no slope passes (the P0-P3
  instrument is a pure interval pass; SlopeG machinery is not
  exercised).
