# ADDENDUM — RETRACTION of the constraint-rigidity local-optimality
# conclusion (Main, 2026-08-30T13:0xZ); rule 5 inline; supersedes the
# "Gate-C verdict" bullet in README Part 2 and the sharpening paragraph of
# repaired_probes_record.md wherever they conflict.

## What is retracted, and where the error entered

The universal claim — "the released point is the UNIQUE feasible point of
each swept single-block family inside the box; local optimality by
constraint rigidity" — is RETRACTED. The error was introduced by Main's
sharpening instruction (the rigidity framing), which I implemented
faithfully; the record shows it entered there. My original wording
("within every feasible single-block family TRIED") was the correct
scope, and six rejected probe directions never licensed the universal
form.

## Why it is wrong — Main's dimension argument, made exact (item 3)

Main's inferred figures were 45 coords / 27+1 rows / nullspace >= 17.
Exact count (SVD rank of the recorded equality machinery for ONE glob
dist block):

- coordinates: 45 (dist_id[r] block);
- margin rows: 3 margins x 9 rows = 27 CANDIDATE rows, but rank of the
  stacked margin matrix is **24** (3 duplicate directions across the
  three margins — singular values show exactly three zeros; each margin's
  own 9 rows are full rank 9, and t1's span adds 0 independent rows over
  {t0, t2}: the three margins overlap by 3 rows total);
- simplex row (ones): rank of margin-span + ones = 24 — the simplex row
  is DEPENDENT on the margin span (it lies in it), adds nothing;
- therefore: 24 INDEPENDENT equality rows on 45 coordinates;
- **exact nullspace dimension = 45 - 24 = 21** (>= Main's 17; the true
  number is 21).

## The decisive structural facts (computed, not asserted)

1. Every nullspace vector is automatically orthogonal to the simplex
   direction: max |v·1| over an SVD basis = 5.55e-16 — the simplex row is
   preserved by ANY kernel move, i.e. the simplex-vs-marginals
   "joint freeze" I wrote is FALSE (the simplex row is a linear
   combination of the margin rows, which is why the repaired probes
   cleared it bitwise — that was a trace of the dependency, not of a
   freeze).
2. All 21 kernel dimensions are two-sided w.r.t. the stored point:
   each basis vector has both strictly positive and strictly negative
   entries, so for the released dist (support min 1.28e-7 > 0) there
   exist small eps with p* ± eps·v strictly inside the box AND strictly
   inside the open simplex AND preserving all 27 marginals — i.e. the
   feasible set inside the box around dist is at least 21-dimensional,
   and my "+7e-8/-7×(−1e-8)" repair pattern was neverin that kernel.
3. The certified objective is NOT constant on the kernel: the marginals
   (hence the global num_block terms and the penalty's H(margin)
   parts) are invariant, but the joint entropies H(dist), H(dist_max is
   fixed), the p_comp complete-split terms and the part-level blocks
   are not; measured float-level slack sensitivity along kernel basis
   vectors at eps=1e-7 (direction-finding only, NOT certified):
   -2.5e-9 .. +6.8e-9 both signs occur — there are kernel directions
   that increase the slack term (float level) and directions that
   decrease it. Whether any direction certifies feasibly below B0 is
   EXACTLY the open question the rigidity claim wrongly closed.

## What survives (unchanged, re-verified by Main)

- The rp no-go: 3 coords pinned by 3 independent rows — zero-dimensional
  feasible set, exactly-feasible |Δr| ~5.5e-17. Genuine rigidity proof.
- The six probes were infeasible as measured (worst marginal move
  7.000097472856237e-08 = 63.64x the 1.1e-9 gate, Main-verified).
- No box certified below published 2.37155181 — no record claim existed.
- Part 1 precision-invariance (bit-identical endpoints 300/200/128/96/64,
  rounding <= 8.3e-19, stricter-standard-not-defect): UNTOUCHED.
- Both accessor bugs and their fixes, closed by three-instrument
  agreement at 3.5258684860650646e-11: stand.

## The honest lesson (Main's sentence, adopted verbatim)

Main caught three shared-instrument audit failures this session — kg's
finite-difference check against its own wrong function, mm3's three
solvers over one modelling premise, our accessor comparing every row to
one b — and then made a dimension-counting error of exactly the same
family: "six directions rejected" was accepted as evidence about a space
because a mechanism story sounded structural. A mechanism that sounds
structural is not a proof that the mechanism is exhaustive.

## README correction (item 2)

README Part 2 wording corrected to the scoped claim; attribution of the
retracted universal claim to Main's instruction recorded there in
writing. Sealed records untouched (nothing unsealed).


## FOLLOW-UP — Main's owner-verification of the rp no-go and the RANK
## DICHOTOMY (adopted as the campaign's standard for any "frozen" claim)

Main re-derived the rp no-go independently: A = [[1,1,1],[0,1,-1],[1,-1,0]],
det A = -3 exactly (my float det −2.9999999999999996; det is exactly −3
over the integers), rank 3 on 3 coordinates → nullspace dimension EXACTLY
0, feasible set a single point, and A x = (1,0,0)ᵀ solves to x =
(1/3,1/3,1/3) — matching the stored triple's role as the unique feasible
point up to dyadic rounding. Owner-verified; recorded in the index.

THE RANK DICHOTOMY (the campaign's general test, supersedes the
"mechanism sounds structural" pattern that produced the retraction):
  - rank == coordinate count → the feasible set is a point, NO direction
    exists; proof complete, zero probes needed.
  - rank < coordinate count → a positive-dimensional family of feasible
    directions exists; NO finite number of rejected probes says anything
    about it; the subspace must be certified.
  - No middle case. Six rejected directions + a structural-sounding story
    was the exact error retracted above.

Each gate-C block placed in the dichotomy (Main's request; my counts):
  1. region_prop: FIRST BRANCH. det A = −9/3? exact −3 ⇒ rank 3 = coords 3
     ⇒ feasible set = {(1/3,1/3,1/3)}. NO-GO PROVEN, no probes needed.
  2. glob_dist (per region): SECOND BRANCH. 45 coords, 27 margin rows of
     rank 24 (3 dependent; spans overlap), simplex row dependent on the
     margin span (adds 0) ⇒ 24 independent rows ⇒ nullspace dimension
     EXACTLY 21 (Main's ≥17 lower bound made exact). All 21 kernel basis
     vectors are orthogonal to the simplex direction (max |v·1| = 5.55e-16)
     and two-sided at p*, so feasible kernel moves exist; the certified
     objective is NOT constant on the kernel (float slack sensitivity
     along kernel basis vectors at eps=1e-7 spans −2.5e-9..+6.8e-9, both
     signs). OPEN: subspace certification (interval directional
     derivatives + certified Hessian-Lipschitz remainder) is the real
     gate-C question; NOT attempted (needs fresh pre-statement).
  3. omega / single_mat_size: SECOND BRANCH (checked, not assumed). NO
     linear equality rows touch either coordinate (verified: 0 rows in the
     recorded ledger); K is pinned by BOUNDS (lb=ub=1), not a row and not
     part of this box. The acting constraint is the Schoenhage line
     value = sms·omega + R = 4 ln 7 (float-exact at p*), ONE nonlinear
     equation on TWO box coordinates ⇒ a codim-1 curve: sampling the box
     at 21 sms values puts 17 of them on admissible line points with
     omega' = M·w/sms' inside the omega side of the box (feasible family,
     dimension 1, though omega' varies only ~±9e-8 over the sms box
     fraction ~9e-8/(2e-7) of the box). NOT a no-go; whether any line
     point also keeps R-side rows satisfied is part of the OPEN
     multi-block question (the line can also be held by moving R, which
     is the compensating-R class).


## RANK CORRECTION + EXACT-KERNEL UPGRADE (Main ORDER, 2026-08-30T13:3xZ —
## certification run STOPPED before compute; basis replacement executed)

Main's blocker accepted in full: the SciPy SVD basis is float, so (i)
"kernel" held only to 5.55e-16 (feasible-in-program-tolerance but breaks
the certified chain) and (ii) — the fatal part — float rank 24 was never
proven exact; a lower true rank would mean MISSED kernel directions, i.e.
incomplete coverage presented as complete — the retracted error one level
up.

EXECUTED (exact, this day):
1. Coefficient extraction: every j2m entry is EXACTLY 0 or 1 (verified:
   distinct values {0.0, 1.0}, all integer — exact dyadics; no
   rationalization needed).
2. Exact rank/nullspace over Q via sympy rref on the 27x45 margin matrix
   (saved margin_matrix_M.npy): EXACT RANK 24, EXACT NULLSPACE DIMENSION
   21 — Main's number CONFIRMED EXACTLY, the "exactly 21" in the
   pre-statement stands (rule 5 correction not needed: float 24 == exact
   24).
3. A · V_exact == 0 verified IDENTICALLY over Q (sympy zero matrix, not a
   tolerance). Exact rational basis saved:
   scratch/kernel_basis_V_exact.npy (float rendering) +
   scratch/kernel_basis_V_exact_rational.json (numerator/denominator
   pairs, 45x21). Float residual of the exact basis: 0.0.
4. Simplex-row dependency proven EXACTLY: rank(margins + ones) = 24 adds
   0 — the simplex row lies in the margin span over Q (previously a float
   observation). Exact-basis float eval gives simplex orthogonality
   EXACTLY 0.0 (the rationals cancel in pairs).
5. Feasibility re-check on the exact basis: all 21 dims two-sided;
   inscribed cube side r_c = 4.7619e-9 (margin factor 21.0 — the exact
   basis is less isotropic than the SVD one; STILL an inscribed subset,
   see Main's second problem below).
6. Certification compute: STOPPED before the 60-minute remainder budget
   per the order. Direction-finding only (NOT certified, NOT an outcome):
   worst single-basis slack delta at an exact-basis box corner is
   −6.33e-9 (j=2, decreasing slack) — the feasibility gate question,
   whether kernel directions can IMPROVE the certified endpoint at
   1e-7-scale, remains genuinely open and now has an exact basis to
   answer it on.

## MAIN'S SECOND PROBLEM ADOPTED INTO THE RECORD — domain naming

The inscribed cube C is a CONSERVATIVE INNER approximation: a certified
cap over C does NOT cap the region between C and the box boundary.
Adopted wording discipline: any certified result must NAME its domain
exactly (e.g. "over the inscribed cube C of side r_c in kernel
coordinates, margin factor 21"), and the region between C and the box is
UNCERTIFIED unless a circumscribed-polytope remainder is later bounded
(NOT inside the 60-minute checkpoint).

## ON RE-CENTERING — rule adopted verbatim

Re-centering may NOT convert a failure-to-certify into a success. Any
re-centered result is a SEPARATE certified statement with its own named
domain, reported alongside the original attempt. Domain-shopping after
seeing the outcome is the box-shifting failure mode.

## FLAT-GRADIENT OBSERVATION (Main)

All |g_j| <= 0.04 (worst −0.0398) suggests a flat profile: a cap over C
will be driven by the SUM of many small contributions across 21
dimensions — precisely the regime where a loose Hessian bound swamps the
linear term. Recorded so the checkpoint report states this rather than
tightening anything under time pressure.

## THE GENERAL LESSON (Main, third instance of the shape)

An approximation standing in for the set it represents: six probes for a
21-dim space; an inscribed cube for its circumscribing region; a float
SVD basis for an exact kernel. The instrument must be able to establish
COVERAGE, not merely be consistent with it.


## THE THREE ROW DEPENDENCIES NAMED (Main's request, 2026-08-30T13:4xZ —
## exact, over Q, from rref of the 27x45 margin matrix; row (t,j) = margin
## t, value j; f_{t,j} = functional summing dist coords with margin-t = j)

- dep0 — GRAND TOTAL across margins t=0 and t=1:
  −Σ_j f_{0,j} + Σ_j f_{1,j} = 0 (weights t0: −1×9, t1: +1×9, t2: 0).
  Structural: each margin's rows partition all 45 coords, so every
  margin-total equals Σ dist = 1·(total). Verified: t0-total − t1-total
  = 0 row; t0-total − t2-total = 0 row.
- dep1 — FIRST MOMENT, weighted −(8+j)/(8−j)/(8−j):
  Σ_j −(8+j)·f_{0,j} + Σ_j (8−j)·f_{1,j} + Σ_j (8−j)·f_{2,j} = 0.
  Size check on the coordinate counts (col sums 9−j each): lhs
  Σ(8+j)(9−j) = 480 = rhs Σ(8−j)·2(9−j) = 480. EXACT.
- dep2 — the same first-moment family, variant with the t2 sign flipped
  and offset 7: Σ_j (7+j)·f_{0,j} = Σ_j (8−j)·f_{1,j} + Σ_j (7−j)·f_{2,j}
  (435 = 435 EXACT). Together dep1 and dep2 span the two first-moment
  relations implied by the triangular (i+j−like) geometry;
  no duplicate rows anywhere (all 27 rows pairwise distinct — checked).

VERDICT: all three dependencies are STRUCTURAL (combinatorial identities
of the joint-to-margin geometry — totals shared by all three margin
partitions, plus the two first-moment relations of the triangular
structure). NO duplicated row, NO mis-transcription. The third
dependency beyond Main's generic three-block model (rank 25/nullspace
20) is the second first-moment relation, which the generic model lacks.
This explains the exact discrepancy rank 24 vs 25 concretely.


## DEPENDENCY ACCOUNTING CLOSED (Main, 2026-08-30T14:0xZ) — the four vs
## three gap and its resolution

Gap found by Main and CONFIRMED EXACT by me: my report named dep0
(t0-total = t1-total), dep1, dep2, AND verified t0-total = t2-total —
four relations, but rank 24 permits exactly THREE independent
left-nullspace vectors. Resolution (Main, verified by me coefficient-
wise over Q): dep1 + dep2 = −Σ_j f_{0,j} + Σ_j f_{2,j}, i.e. the second
total relation is the combination dep1 + dep2. {dep0, dep1, dep2} is a
BASIS of the left-nullspace; the count closes with the proven rank 24
by an argument independent of the rref (second instrument).

GEOMETRIC FACT AS CORRECTED: the three margin relations are ONE shared
total relation (t0 = t1) plus a 2-dimensional first-moment family whose
sum reproduces the second total relation (t0 = t2). The "two totals +
two moments" picture was an over-count; the moment span CONTAINS the
second total relation. Main's size-check re-verification: 480 = 480 and
435 = 435 both hold.

PROCESS RULE (adopted per Main, recorded here): "OWNER-SUPPLIED FIGURES
ARE INPUTS, NOT EVIDENCE" — any claim resting on one must be
independently re-derived before it is written down. Applied to my own
dependency report: I announced four relations without checking their
span count against the proven rank; caught by Main's accounting, now
closed. LESSON LINE (Main, adopted): "A round number in a scaling
factor is a basis artifact until proven geometry" — generalizes to any
suspiciously round constant in a pipeline: a property of the
representation until proven otherwise.

DUAL-CERTIFICATE DISCIPLINE (Main condition, recorded into the
certification plan): the LP dual point must be verified for DUAL
FEASIBILITY in exact (rational) or certified arithmetic — not merely
plugged into the duality expression. Float duals get rounded to exact
rationals and every dual constraint verified; if rounding breaks dual
feasibility, perturb conservatively in the safe direction (a slightly
suboptimal exactly-dual-feasible point gives a valid, possibly loose,
rigorous bound — a valid loose bound beats an invalid tight one).

CLOCK STARTED: 60 minutes, Hessian remainder only. LP-dual step runs
first (minutes). Report at checkpoint or completion.
