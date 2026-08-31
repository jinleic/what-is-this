# Pre-statement — gate C stage 3: certified gradient (forward-mode AD/slope)
# + certified second-order remainder, deciding the 21-dim nullspace question.
#
# Agent OmegaSlope, written 2026-08-30T17:5xZ, BEFORE the first computation
# of this campaign (only instrument-familiarization reads and pure-numpy
# reproduction of FROZEN predecessor numbers happened before this file was
# committed; no new numeric evidence was produced).
#
# Campaign base: continues cs/omega (frozen predecessor campaign
# campaigns/2026-08-30T11:24:41Z_17bdaeb6_099b1224d603, outcome FAILURE TO
# CERTIFY; its assets are reused EXACTLY as frozen). New frozen campaign
# directory: campaigns/2026-08-30T17:5xZ_OmegaSlope_<hash>/ (producer-
# generated name fixed at freeze time; this file is committed there first).

## 1. The question (verbatim from omega/README.md:339-356 — unchanged)

Is there a feasible improving direction in the EXACT 21-dimensional dist
nullspace around the VXXZ24 K100_2.37155181 released parameter vector?
Domain fixed by the predecessor and NOT moved (rule 16):

    D = { delta : A·delta = 0 exactly, |delta|_inf <= 1e-7 }

acting on the region-0 glob dist block (45 coordinates), A = the recorded
27x45 margin matrix (exact rank 24 over Q; kernel basis V, 45x21, exact
rational, kernel_basis_V_exact_rational.json). delta = V·c,
c in R^21, D = {Vc : |Vc|_inf <= 1e-7}.

Program: the transcribed 3-region single-p_comp verifier
(src/vxxz24_float.py, max_level 3, q = 5) under the frozen certified
aggregation semantics (interval lows for num_block/p_comp, highs for
penalties, Omega = (4 ln 7)^up − R_sum_low, all divided by M_low, plus the
FROZEN point defects 1.9380875689560958e-11 / M_low absorbed — the frozen
campaign absorbs them at the box CENTRE; this campaign must handle the
defect movement over the box, see §4c).

## 2. Machinery to be built (named pre-compute)

(a) Forward-mode interval automatic differentiation (slope arithmetic)
through the SAME transcription: each leaf is an (Ival value, gradient
interval vector over the 45 dist coordinates) pair or a certified slope
enclosure over the box; all field ops (+,-,*,/, ent, nent, ln-like
residuals) get interval derivative rules in Arb with explicitly controlled
radii. The arb(lo,hi) midpoint/RADIUS landmine is respected: intervals are
built from arb midpoint + radius we control, never from float endpoints.

(b) Certified second-order remainder over D (= the whole box on the dist
block; valid and conservative — the kernel constraint can only shrink the
achieved set): |R(delta)| <= (1/2)·||delta||^2-sup over the box × a
certified bound on the Hessian/operator norm of the interval gradient
(enclosure of sup over the box of ||J_G||). If the Hessian enclosure is
too loose for the remainder to sit below the signal, that is the REPORTED
finding (failure to certify at step (b)), not a patched tolerance.

(c) Combination: certified slope enclosure S(delta) with
    Omega(p*+delta) in Omega0 + S(delta) + [-Rmax, +Rmax]
for EVERY delta in D. S encloses the first-order term restricted to the
kernel: use ANY-y instrument G·delta = (G − A^T y)·delta >= −1e-7·
||G − A^T y||_1 (A·delta = 0 identically kills the y-term; NO dual
feasibility needed). Decide:
  * certificate (a): if min over D of the certified enclosure >= Omega0's
    certified endpoint 2.3715538358544617 (with the defect treatment),
    then NO feasible improving direction exists in the exact 21-dim
    kernel at radius 1e-7 — certified local optimality, named domain, one
    sentence;
  * certificate (b): if a delta in D with certified Omega < endpoint −
    (width) exists, exhibit it exactly and its certified decrease;
  * failure mode: report the exact ratio (achieved enclosure width) /
    (1.5816497000997742e-07 signal) = how much tighter the machinery must
    be, against the predecessor's 122x baseline.

## 3. Budget

Wall clock: 4 hours for the AD build + 1 hour for the Hessian/remainder
step (hard checkpoint, write-first protocol, state dumped to the campaign
directory before stopping). No subbox sampling of D beyond what the
certified statements themselves cover. Precision: MID = 300 bits working
(the frozen default; drop to 200 only for speed with re-verifying B0
bit-identity at 300 first). No project-wide test/format runs. Threads = 1.

## 4. Falsifiable outcome / adjudication rule (fixed NOW)

- PASS-a (no improving direction): certified over D, width of the total
  enclosure (slope hull + remainder) STRICTLY BELOW the distance from the
  float-gradient LP optimum to Omega0 — i.e. the machinery resolves the
  1.58e-07 signal — AND the min of the enclosure over D is >= Omega0.
  Then the sentence is: "no feasible improving direction exists in the
  exact 21-dimensional nullspace at |delta|_inf <= 1e-7" — with the
  caveat that the certificate covers D, not any larger set.
- PASS-b (improving direction exists): an exact delta (rational
  coordinates given exactly) in D with certified Omega(p*+delta) +
  Rmax < Omega0. Reported with its certified decrease.
- FAIL: report the achieved (width / signal) ratio and the exact step
  where the width blows up; the open question stays open with a quantified
  distance. No domain move, no re-centering, no budget extension.

## 5. Scope honesty (binding, verbatim commitment)

Even total success inside D caps 12.8x short of the published
2.37155181 (float-gradient LP bound: max achievable improvement in D is
1.58e-07 while the gap Omega0-to-published is 2.0258544615181506e-06:
ratio 12.808). **This route can NEVER produce a record and no record
claim may appear anywhere.** This target does not claim to improve omega;
a certified enclosure is not a new bound. Nothing written here bears on
the validity of the VXXZ24 paper (escalation rule: anything that WOULD
bear on it goes to Main by hub first and enters no file until Main
answers).

## 6. What is NOT searched (pre-declared)

- The other two glob dist blocks (regions 1, 2): same machinery applies
  but the pre-registered domain of THIS campaign is the region-0 block —
  the same block the 21-dim kernel was proven on. Multi-block classes
  remain open.
- All other parameter coordinates (lam multipliers, splits, part-level
  dists, region_prop, omega, single_mat_size): held at p* values as
  POINTS, not boxed. The known 2-coordinate Schoenage-line class
  (omega × single_mat_size) is NOT swept.
- The published 2.37155181 and the record 2.371177: not touched, not
  searched; any statement here is about OUR certified enclosure value.
- No subdivision of D; no Monte Carlo; no probes outside D.

## 7. Rule-7 sentence (the scope of any certified sentence this campaign)

The certified statements of this campaign apply to the exact swept set
D = {delta in R^45 : A·delta = 0 exactly, |delta|_inf <= 1e-7} (A = the
frozen 27x45 margin matrix, rank 24 exact over Q, kernel basis V exact
rational) around the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector, under the transcribed 3-region
single-p_comp program at max_level 3, q = 5, with the fixed construction
and parameterization enclosed; they do NOT apply to the matrix
multiplication exponent record omega < 2.371177 (parameters unpublished as
of 2026-08-30, never searched), to any point outside D, to any other rung,
to any other max_level/q regime, to multi-block or cross-block direction
classes, or to any construction differing from the enclosed laser-method
program (including the 2026 asymptotic-rank / centroid-based improvements
outside every statement here).
