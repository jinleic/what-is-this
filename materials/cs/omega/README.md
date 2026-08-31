# `omega/` — rigorous enclosure of the matrix multiplication exponent record

**Status: BENCHMARK (no campaign run yet).** The current world record on the
matrix multiplication exponent $\omega$ is produced by a machine-learning-driven
numerical optimization. No rigorous interval enclosure of it is published.

## The claim

Owner-read (**V-owner**, 2026-08-29) — arXiv:**2608.16884**, Dupont,
Eisenberger, Kozlovskii, Mehrabian, Ruiz, See, Zhou, Alman, Vassilevska
Williams, Balog, *"Improving the matrix multiplication exponent with modern
optimization and AlphaEvolve"* (2026-08-17), abstract verbatim:

> The current best bounds on the matrix multiplication exponent $\omega$ are
> obtained through a refinement of the laser method called combination loss
> analysis (Duan et al., 2022; Williams et al., 2024; Alman et al., 2025). In
> this note, we address the optimization problem at the core of this approach
> and propose several improvements. First, we reformulate the optimization
> problem allowing us to solve it in a larger setting than was previously
> possible. Second, we leverage recent advances in machine learning to design a
> new optimization algorithm for this problem. Finally, we refine the resulting
> optimization algorithm with AlphaEvolve. Our combined approach yields an
> upper bound of $\omega$ < 2.371177, improving the previous best bound of
> 2.371339.

The most-cited constant in algorithm design now rests on the numerical optimum
of a reformulated combination-loss program, found by an ML optimizer refined
with AlphaEvolve.

## Why this is a certificate target, and why it is unowned

An upper bound on $\omega$ from combination loss analysis is mathematically a
*feasible point* plus an inequality: given the parameter vector, the bound
follows by evaluating an explicit expression built from entropies, multinomial
coefficients and a constrained maximization. A feasible point found in floating
point does **not** by itself establish the inequality — the evaluation must be
outward-rounded, and every inner maximization must be bounded above rigorously
rather than sampled. That is a bounded, well-posed `python-flint`/Arb job.

Nobody is positioned to do it: DeepMind's contribution is the optimizer, not a
kernel-grade enclosure, and the theory co-authors (Alman, Vassilevska Williams,
Zhou) are authors of the record itself, so an independent audit is by
definition external. The genre — "we certify an existing constant" — carries no
conference incentive, which is exactly why the niche is open.

**Framing discipline (repo rules 9 and 12):** this target does **not** claim to
improve $\omega$. A certified enclosure is not a new bound. The deliverable
answers one question: *is the published record rigorously established as
stated?*

## Gates

1. **Gate A (reconstruct and enclose the record, hours–days).**
   Reconstruct the combination-loss optimization exactly as reformulated in
   arXiv:2608.16884; take the paper's reported optimal parameters; re-evaluate
   the resulting $\omega$ bound in Arb with outward rounding, bounding every
   inner maximization rigorously from above.
   **Pass:** the certified interval's upper endpoint is $<2.371177$ — report the
   margin. **Refutation:** the enclosure does not close below $2.371177$;
   report the certified upper endpoint and the exact step where the
   float-to-interval gap opens. That would mean the record is not rigorously
   established as stated — a significant finding requiring owner re-audit
   before any claim leaves this repo.
2. **Gate B (the predecessor ladder — where the certificates are actually
   reachable, days).**
   **Correction, 2026-08-30, from agent `Omega`'s primary reads — this replaces
   an earlier wrong statement in this README that "the method's published
   constants are all float-only".** They are not. All three predecessors
   released code *and* parameters on OSF:

   | rung | bound | artifact |
   |---|---|---|
   | Duan–Wu–Zhou (2022) | $2.371866$ | `osf.io/dta6p` — `power8_dup_2.371866.mat`, MATLAB MCOS objects (may be unreadable without MATLAB) |
   | Vassilevska Williams–Xu–Xu–Zhou (2024) | $2.37155181$ | `osf.io/7wgh2` — `rmmcode.zip`, `data/K100_2.37155181.mat` = flat 6759-dim float64 vector + MATLAB verify script with a $10^{-9}$ refine mode |
   | Alman et al. (2025) | $2.371339$ | `osf.io/mw5ak` — `code_matrix_mult.zip`, `data/W1.00_2.371339.mat` = flat 24855-dim float64 vector + verify script |

   So the predecessor rungs are **reproducible in floating point but have never
   been enclosed in interval arithmetic**. That is the real gap, and it is a
   different and weaker claim than "float-only". Gate B: re-evaluate the
   released parameter vectors in Arb with outward rounding, bounding every
   inner max-entropy rigorously from above. A success is the **first rigorous
   interval enclosure of any rung of the combination-loss ladder** — the
   publishable core of this target regardless of how $2.371177$ resolves.
   Priority order: Alman25's $2.371339$ first (the immediate predecessor and
   the universally cited number), then VXXZ24's $2.37155181$.

   **Mandatory two-stage transcription control.** The released verifiers are
   MATLAB and there is no MATLAB on this machine, so the evaluation code must
   be transcribed to Python. A transcription bug can masquerade as a failed
   enclosure, and a false "the record does not certify" is the worst output
   this repo could emit. Therefore: (a) transcribe to plain float64 Python
   first and reproduce the published value from the released vector to the
   paper's own tolerance ($\sim10^{-9}$), recording reproduced vs published;
   (b) only then switch the *same* code path to Arb. An enclosure failure after
   a clean stage (a) is evidence about the bound; an enclosure failure without
   a clean stage (a) is evidence about the transcription and must be reported
   as such.
3. **Gate C (certified local optimality, days).**
   Interval-certified branch-and-bound over a box around the reported optimum.
   **Pass:** a certified improvement — a genuine new bound, reported only with
   its full certificate. **Negative result:** a certified cap for the box, i.e.
   the reported point is optimal within the searched region — publishable as
   the first rigorous local-optimality statement for the current record.

## Pre-campaign requirements

- Owner first-hand read of the **full** arXiv:2608.16884 plus the three
  predecessor papers it names. Extract: the exact reformulated program, the
  reported optimal parameters **and whether they are published to enough digits
  to be re-evaluated at all**, the constraint set, and every inner maximization.
- If the parameters are not published at sufficient precision and no code or
  parameter file is released, gate A's finding becomes *"the record is not
  independently checkable from the paper as published"*. Report that as the
  result; do not reverse-engineer a substitute point and present it as a check.
- One-page pre-statement: precision, rounding mode, the rigorous upper-bounding
  strategy per inner maximization, and the pass/fail threshold $2.371177$.

## Disjointness

No `../math/` or `../physics/` target touches matrix multiplication or the laser
method. In-repo neighbour `../mm3/` concerns the *additive* complexity of a fixed
rank-23 $3\times3$ algorithm — a finite circuit-optimality question, not the
asymptotic exponent. Share only exact-arithmetic tensor utilities.

## Layout

- `src/`, `campaigns/`, `scratch/` — created at first use.

## Current state (agent `OmegaRung2`, 2026-08-30; wording fixed per Main
owner re-check, 4-part framing)

**POSITIVE RESULT FIRST — the two-rung ladder ORDER survives rigorous
arithmetic.** Enclosing both released parameter points under identical
pre-committed semantics (interval arithmetic with outward rounding, Lemma-1
dual-certificate residual charged once through the hashing constraints,
float defects absorbed explicitly), rung 1 (Alman25) certifies to
**omega <= 2.3713400836689** and rung 2 (VXXZ24) to
**omega <= 2.3715538358544617**, a gap of **2.1375e-4**. Alman25's rung
genuinely improves on VXXZ24's even when both are enclosed rigorously with
explicit slack — the first two-rung rigorous comparison in the
combination-loss line. MACHINE-VERIFIED, deterministic (byte-identical
re-runs). Both certified endpoints also sit below DWZ23's 2.371866
(gaps >= 3.12e-4).

### Rung 2 (VXXZ24, K100_2.37155181) — what was done

Mandatory two-stage protocol (pre_statement A2/A4) applied to the released
OSF vector (osf.io/7wgh2, `data/K100_2.37155181.mat`, sha256 `df75ae3a…93d`):

- **Stage (a), [REPRODUCED], MACHINE-VERIFIED.** Transcription of the
  released VXXZ24 verifier to `src/vxxz24_float.py` (3 regions, single
  p_comp, VXXZ24 registration order). Parameter count 6759 == 6759 exact;
  c_max = 4.6349590904e-12, ceq_max = 3.5953524880e-11, linear violations
  0.0, Schoenage line `value − 4 ln 7` = 0.0 EXACTLY, value =
  7.783640596221253; all inside the harness's own refine tolerance 1.1e-9.
- **Stage (b), MACHINE-VERIFIED.** Interval enclosure with outward rounding
  (`src/stage_b_rung2.py` + `interval_core.py`): per-block containment 9/9
  num_block + 3/3 penalty + 1512/1512 Lagrange form-consistency; precision
  sweep live (num_block widths 1.084e-19 → 4.770e-18 as precision drops);
  double-count test run first — toggling the Lemma-1 charge moves the raw
  endpoint by the charge itself (2.0297e-6), so it is present exactly once
  and never added twice (CORRECTION_coincidence.md protocol).

### Rung 2 verdict (frozen campaign
`campaigns/2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/`; owner re-check
confirmed all numbers at 200 bits)

**Part 1 — WHAT IS ESTABLISHED.** A rigorous interval enclosure at VXXZ24's
released parameter point (6759-wide float64 match, stage-(a) reproduction
clean, see below) under the pre-committed Lemma-1 semantics gives

    omega <= 2.3715538358544617
           = 2.3715538358350807 (raw interval endpoint)
             + 1.9380875689560958e-11 (c + ceq + Schoenage defects / M_low)

That is a valid rigorous upper bound. It sits 2.0258545e-6 ABOVE VXXZ24's
published 2.37155181, so it does not reproduce their six-decimal rounding.

**Part 2 — WHY, WITH ATTRIBUTION.** The shortfall is entirely the ln-space
dual residual of the SHIPPED Lagrange multipliers — worst 3.0767e-6 in glob
region 2, 2.1502e-6 at region 0, 1.1433e-6 at region 1, all ~1500 smaller
per-part residuals below ~7e-8 — measured against interval widths of
1.08e-19. The arithmetic is nowhere near the bottleneck: it is five orders
of magnitude tighter than the gap. The diagnostic arm (Lemma-1 charge off)
certifies 2.3715518061863814 (+1.94e-11), i.e. with tight multipliers the
published value IS recoverable, with margin 3.794e-9. **This is a statement
about the looseness of the released certificate data, NOT about the
validity of the VXXZ24 bound, and NOT about any error in the paper** — the
enclosure says nothing negative about VXXZ24's mathematics, which derives
from exact max-entropic dist_max rather than from certificates the released
float vector happens to carry.

**Part 3 — THE GENERAL FINDING (first measurement of its kind).** Both
rungs behave identically: the cost of rigour is 1.084e-6 on Alman25's rung
and 2.026e-6 on VXXZ24's, in both cases dominated by the Lemma-1 residual
of the shipped multipliers rather than by arithmetic. The published
constants in the combination-loss line are float-optimal points whose
accompanying dual certificates are loose at the ~1e-6 level; rigorously
enclosing the released data costs about one to two parts in 1e6. This is a
quantitative, reusable statement about the artifact quality of the method's
public releases, measured here for the first time.

### Two-rung rigorous ladder (the payoff — achieved)

Identical rigorous semantics on both rungs (interval + Lemma-1 dual
certificates + absorbed defects):

| rung | certified endpoint | published | published met? |
|---|---|---|---|
| DWZ23 (2022) | (MCOS-opaque, not evaluated) | 2.371866 | — (below it though) |
| VXXZ24 (2024) **rung 2** | 2.3715538358544617 | 2.37155181 | NO (+2.03e-6) |
| Alman25 (2025) **rung 1** | 2.3713400836689 | 2.371339 | NO (+1.08e-6) |

- **Two-sided ladder statement:** rung 1's certified endpoint is below
  rung 2's certified endpoint by 2.1375218556e-4 — Alman25's improvement
  over VXXZ24 survives symmetric rigorous treatment. Also both rungs sit
  below DWZ23's 2.371866 (gaps ≥ 3.12e-4). First two-rung rigorous
  comparison in the combination-loss line. MACHINE-VERIFIED.
- Rung 1 numbers carried from `campaigns/2026-08-30T03:20:00Z_c9d4e2a8`
  (CORRECTION_coincidence.md); unchanged here.

### Honest gaps

- Neither rung's published six-decimal rounding is reproduced at the
  SHIPPED multipliers under the repo's dual-certificate standard (rung 2 by
  +2.026e-6, rung 1 by +1.08e-6). Per Part 2/3 this reflects released
  certificate looseness, not bound invalidity; multiplier re-optimization
  (diagnostically sufficient to recover rung 2's rounding, margin 3.79e-9)
  was NOT attempted — a follow-up campaign could upgrade "not reproduced
  from the released data" to "reproduced with recomputed multipliers".
- DWZ23 rung not attempted (MCOS .mat opaque; per gate-B priority note).
- `interval_core.py` gained a `float_module` parameter (rung 1 regression:
  24855-param tree builds and evaluates, both orders of evaluation). All
  arithmetic MUST run inside `IC.prec()` — probing outside it fails, as
  expected for ball arithmetic at default 53-bit precision.

### How to run

See `campaigns/2026-08-30T05:20:00Z_e97c35ae_70c28fc61780/REPRODUCTION_COMMAND.md`.
Deterministic: both arms re-ran byte-identical from the frozen directory.

## Current state (agent `OmegaGateC`, 2026-08-30) — gate C, VXXZ24 rung

Frozen campaign: `campaigns/2026-08-30T11:24:41Z_17bdaeb6_099b1224d603/`
(pre-statement, scripts as run, raw outputs, checksums, manifest;
pre_statement committed before the first B&B call). Centered on the
**VXXZ24 rung** per Main steering: the record omega < 2.371177
(arXiv:2608.16884) has NO published parameter vector (artifact audit
re-run 2026-08-30: google-deepmind/alphaevolve_results carries no
parameter table; GitHub searches return 0), so there is no reported
optimum of the record to enclose. The reported optimum enclosed here is
the VXXZ24 released vector `K100_2.37155181.mat` (osf.io/7wgh2, sha
df75ae3a…, 6759 float64 params), whose float transcription
`src/vxxz24_float.py` was already parameter-count-verified exact.

### Part 1 — the 2.0258545e-6 gap is NOT rounding (5-precision sweep)

Same interval code path re-run at MID = 300/200/128/96/64 bits, both
arms (Lemma-1 charged and diagnostic-off), at p*. The certified endpoint
is bit-identical at all five precisions: 2.3715538358544617 (with arm;
gap +2.0258544615181506e-06) and the without arm lands
−3.794237812826395e-09 below the published rounding at every precision.
R_sum widths fall 1.57e-89 → 1.73e-18 across the sweep while the gap
does not move: the rounding contribution to the gap is measured
<= 8.3e-19 of omega at 64-bit working precision. **The gap is
structural**, with the refined attribution: it is the ln-space Lemma-1
residual of the shipped multipliers charged by our strict certified
standard (glob-region residuals 2.1502e-6 / 1.1433e-6 / 3.0767e-6
before region_prop weighting), which the released float64 verifier
never had to charge (its own exp-form ceq residual at p* is 3.60e-11
against a 1.1e-9 tolerance). The release is not defective relative to
its own standard; the diagnostic without arm recovers the published
rounding with 3.79e-9 margin at every precision. This supplements (does
not contradict) the 05:20:00Z Part 2 wording.

### Part 2 — certified box screen (the gate-C deliverable)

Pre-registered box: B = {p* + offset : |offset_i| <= 1e-7 in EVERY one
of the 6759 coordinates}. Stage (a) float screen ran first (clean
anchor; all probes recorded including exactly-infeasible ones — see
stage_a_screen_record.md). Stage (b) evaluated B0 + 14 signed structured
boxes, each a genuine interval box (exact-dyadic endpoints,
outward-rounded aggregation identical to the frozen rung-2 protocol;
zero-crossing guard clean on every box):

- **B0 anchor reproduced** the frozen rung-2 endpoint: raw
  2.3715538358350803 / certified 2.3715538358544612 (float display
  granularity of …0807/…4617; R_sum width 7.85e-90, M width 2.49e-21,
  value width 0.0). MACHINE-VERIFIED.
- **No box certifies below the published 2.37155181** — no
  world-record improvement, nothing escalated, consistent with
  pre-statement §8.
- **RETRACTED (feasibility audit; campaign feasibility_audit_addendum.md):
  the earlier report that three boxes "certify strictly below B0" does NOT
  improve any bound.** Feasibility audit (exact-dyadic check of the linear
  equality rows) shows all three below-B0 centres — region_prop+
  (2.3715530208533490, -8.15e-7), glob_dist_2+ (2.3715534164727120,
  -4.39e-7), glob_dist_0+ (2.3715537111130380, -1.44e-7) — are
  EXACT-INFEASIBLE points of the program: region_prop+ breaks the
  sum(rp)=1 row by +3.000000000086e-07; each glob_dist_r+ breaks its
  region's dist-simplex sum by +7.0e-07 and the dist/dist_max marginal
  equalities at the 1e-7 scale. They are certified objective EVALUATIONS
  at inadmissible points and bound nothing about omega. The certified
  bound for the VXXZ24 rung remains 2.3715538358544617 (frozen; this
  campaign's tighter re-evaluation of the same quantity:
  2.3715538358544612). The probes stand as upper-sensitivity evidence
  only.
- **RETRACTED (2026-08-30T13Z, Main; record: campaign
  retraction_addendum_rigidity.md).** The "LOCAL OPTIMALITY BY CONSTRAINT
  RIGIDITY" universal conclusion announced below the fold of this entry
  was MY (OmegaGateC's) implementation of a Main sharpening instruction,
  and the instruction was WRONG — the error entered at that instruction
  and the record shows it. Exact dimension count (SVD rank, one glob
  dist block): 45 coordinates, 27 margin rows of RANK **24** (3 dependent
  directions; each margin's own 9 rows full-rank but the three margins
  overlap by 3), simplex row DEPENDENT on the margin span (adds 0
  rank) — therefore **24 independent equality rows and a 45-24 = 21
  dimensional EXACT nullspace** (Main's 17 was the lower bound). Every
  kernel basis vector is orthogonal to the simplex direction (max |v·1|
  = 5.55e-16) and two-sided, so p* ± eps·v stays feasible (box + open
  simplex + all marginals exact) for small eps: the feasible set inside
  the box around dist is >= 21-dimensional, and the certified objective
  is NOT constant on the kernel (float-level slack sensitivity along
  kernel basis vectors at eps=1e-7 spans −2.5e-9..+6.8e-9, both signs).
  "Six directions rejected" was accepted as evidence about a space; a
  mechanism that sounds structural is not a proof that the mechanism is
  exhaustive.
- **What SURVIVES of gate C, as the scoped result** (MACHINE-VERIFIED
  unless noted): (a) the rp no-go — region_prop (3 coords) pinned by 3
  INDEPENDENT rows (sum=1, Y−Z=0, X−Y=0), zero-dimensional feasible
  set, exactly-feasible |Δr| ~5.5e-17, 10 orders below the 1.94e-11
  defect absorption: this rigidity proof is exact and stands; (b) the
  six pre-registered repaired probes are infeasible as measured (worst
  marginal move 7.000097472856237e-08 = 63.64x the 1.1e-9 gate);
  (c) no box certified below the published 2.37155181 (never a record
  claim); (d) Part 1's precision-invariance of the certified endpoint at
  MID 300/200/128/96/64 (rounding contribution <= 8.3e-19; the 2.026e-6
  gap is the stricter-standard Lemma-1 cost, not a defect). The gate-C
  status is **PARTIAL**: whether a feasible improving direction exists
  in the exact 21-dimensional dist nullspace (or a multi-block class) is
  OPEN. Subspace-certification attempt (Main-approved, pre-registered
  Amendment A, run 2026-08-30, outcome FAILURE TO CERTIFY and reported
  as such): over the true domain D = {delta : A.delta = 0 exactly,
  |delta|_inf <= 1e-7} the naive interval enclosure is
  omega in [2.3715431556360005, 2.3715624359361276] (width 1.928e-5,
  rigorous but 122x looser than the 1.5816497000997742e-07 first-order
  signal), so the cap cannot be decided without a certified
  gradient/slope enclosure through the transcribed tree plus a
  second-order remainder — new machinery, not built. Assets produced:
  exact rational kernel basis with rank-24 proof over Q, three named
  structural margin dependencies (one total relation plus a
  2-dimensional first-moment family whose sum reproduces the second
  total relation), and a reusable rigorous linear-bound instrument
  (for ANY y, G.delta >= -R*||G - A^T y||_1 on the kernel, no dual
  feasibility needed). Even total success in D caps 12.8x short of the
  published 2.37155181, so this route can never produce a record.
- Run distinction kept as before: statements concern the VXXZ24 rung at
  2.37155181, NOT the record 2.371177 (parameters unpublished as of
  2026-08-30, never searched).

### (RETRACTED TEXT, kept for the record per rule 5 — the
### rigidity-framing bullets removed from the verdict above; see
### retraction_addendum_rigidity.md for the full text and the error's
### provenance)

### Rule 7 sentence (mandatory)

The certified statements above apply to the exact swept set
B = {p* + offset : |offset_i| <= 1e-7} in R^6759 around the VXXZ24
K100_2.37155181 released parameter vector, under the transcribed
3-region single-p_comp program at max_level 3, q = 5, with the fixed
construction and parameterization enclosed; NOT to the matrix
multiplication exponent record omega < 2.371177, which was not and could
not be searched because its parameter vector is unpublished as of
2026-08-30; NOT to any point outside the box, any other rung's released
vector, or any other max_level/q regime; and NOT to any construction
differing from the enclosed laser-method program — including the 2026
asymptotic-rank / centroid-based construction improvements
(arXiv:2605.21738, arXiv:2608.27434), which move the landscape rather
than move within it and are outside every statement here.

### What was NOT swept (explicit)

- All 2^6759 corners of B; only B0 + 14 structured signed axes of a
  63-box budget were evaluated; no depth-3 interval subtree.
- Equality-repaired dist-simplex probes: the raw glob_dist probes move
  all 7 shape coordinates at once and break the simplex equality at
  7e-7 — their certified endpoints are valid enclosures of their own
  boxes, but the step they represent is not a feasible point of the
  exact program. Feasible-repair probes remain unexplored.
- Multiplier re-optimization (fresh dual derivation — the named
  follow-up from 05:20:00Z, unchanged).
- The Alman25 rung-1 box; the DWZ23 rung (MCOS-opaque); the record
  2.371177 (no parameters); any deeper B&B tightening.

### How to run (this campaign)

See `campaigns/2026-08-30T11:24:41Z_17bdaeb6_099b1224d603/manifest.md`.
Deterministic, OMP/OPENBLAS threads = 1, `nice -n 10`; scripts as run
are byte-copies in the frozen directory with sha256 checksums.

## Current state (agent `OmegaSlope`, 2026-08-30) — gate C stage 3: the
## certified slope/second-order machinery

Frozen campaign: `campaigns/2026-08-30T21:59:00Z_OmegaSlope_7d90d14bc58f/`
(pre-statement `pre_statement_gateC2.md` committed to git 83e3f09 BEFORE
the first computation; scripts as run as byte-copies with sha256
checksums; raw checkpoint + stdout of the final failing run kept).

**Outcome: FAILURE TO CERTIFY — MACHINERY INCOMPLETE, reported as such per
the pre-statement's pre-registered FAIL mode. The 21-dimensional
nullspace question is UNCHANGED AND OPEN.** No decision (neither local
optimality nor an improving direction) was reached; nothing below is a
claim about omega.

### What was BUILT and what WORKS (all MACHINE-VERIFIED unless noted)

- `src/gate_c_slope_core.py`: a certified forward-mode interval-AD core
  (Slope = {value interval, gradient-slope array of 45 Ivals}) with
  Leibniz product, quotient, and entropy (`-t ln t`) chain rules, all in
  Arb endpoint intervals with buffered outward rounding
  (`arb(lo,hi)` midpoint/radius landmine respected). ~110 s per pass.
- `src/gate_c_slope_pass.py`: the certified REPLAY of the transcribed
  3-region single-p_comp program (max_level 3, q 5) carrying Slope
  payloads through every consumed quantity: glob dist/dist_max entropies
  and marginals, the Lemma-1 eps residual terms (both glob and part
  level), part_frac chains, 135 level-3 parts' num_block/penalty/p_comp,
  1377 level-2 contributions, and the value/Schoenage line slopes. Two
  passes run end to end (point r=0 and box r=1e-7).
- `src/gate_c_endgame.py`: endgame driver — value re-aggregation with
  branch pinning, validations V1–V4, scaled + unscaled LP over D, and
  the decision logic from the pre-statement.

### What FAILED, with the exact step and attribution

- Validation **V1 FAILED**: the point (r=0) value re-aggregation atop the
  slope pass returns Omega_raw = 41.9962 against the frozen certified
  2.3715538358350803 (delta +39.62). Root cause isolated: the replay
  never builds the level-2 PartLv2 `mat_size_contribution` — the joint
  ~2.094-point of mat_size — so M_low collapses 2.09425 → 0.14985 and the
  aggregate endpoint inflates ~14x, with R-side orderings correspondingly
  perturbed (branch pattern differs from the frozen aggregation's).
  Gradient scales derived from that wrong value definition (~100s per
  coordinate) are equally attached to the wrong field and are NOT the
  missing certified gradient enclosure. Nothing from those runs is
  evidence about Omega at any evidence level.
- Budget note: the value-aggregation gap is small, named and exactly
  localizable (see manifest "What the next agent must fix"); the budget
  ended before it could be repaired, and per the pre-statement §4 FAIL
  mode no repair-and-continue, no budget extension, and no partial
  numbers were allowed to propagate.

### State of the open question (unchanged)

- The predecessor's scoped negative STANDS: no feasible improvement among
  the six pre-registered single-block probes at ±1e-7; region_prop admits
  no feasible motion (exact 3-row pin); the certified box enclosing D is
  [2.3715431556360005, 2.3715624359361276] (width 1.928e-5), rigorous but
  122x looser than the 1.5816497000997742e-07 first-order signal.
- Whether a feasible improving direction exists in the exact
  21-dimensional dist nullspace at |delta|_inf <= 1e-7 remains OPEN.
  The instrument required — certified slope (mean-value form) plus a
  certified second-order remainder over D — is now ~90% built and
  exactly one value-aggregation repair away from its first real
  decision run.

### Rule 7 sentence (mandatory)

Nothing in this section is certified beyond the scripts' own internal
arithmetic: the swept set of THIS attempt was D = {delta : A·delta = 0
exactly, |delta|_inf <= 1e-7} around the region-0 glob dist block of the
VXXZ24 K100_2.37155181 released vector (A = the frozen 27x45 margin
matrix, rank 24 exact over Q; kernel basis exact rational), under the
transcribed 3-region single-p_comp program at max_level 3, q = 5, with
the fixed construction and parameterization enclosed — and even over D
no decision was reached (machinery incomplete, see above). The
statements here do NOT apply to the matrix multiplication exponent
record omega < 2.371177 (parameters unpublished as of 2026-08-30, never
searched), to any point outside D, to any other rung's released vector,
to any other max_level/q regime, to multi-block or cross-block direction
classes, or to any construction differing from the enclosed laser-method
program (including the 2026 asymptotic-rank / centroid-based
construction improvements, outside every statement here). Even total
success in D caps 12.8x short of the published 2.37155181, so this route
can NEVER produce a record and no record claim may appear anywhere.

### What was NOT swept (explicit)

- No certified decision over D was produced: neither PASS-a (no
  improving direction) nor PASS-b (an exhibited improving direction) was
  reached; the completion distance is the named value-aggregation repair
  (manifest) plus a re-run of the endgame.
- Regions 1 and 2 glob dist blocks; the multi-block/cross-block classes;
  the omega x single_mat_size Schoenage-line 2-coordinate class; all
  other parameter coordinates (lam multipliers, splits, part-level
  dists, region_prop) — held at p* points, not boxed.
- The published 2.37155181 and the record 2.371177: not touched, not
  searched; nothing written anywhere in this campaign bears on the
  validity of the VXXZ24 paper (no escalation was triggered and none was
  needed — no claim about the paper was written).
- No subdivision of D, no Monte Carlo, no probes outside D, no float SVD
  basis substitution (the exact rational kernel basis was used as-is).

## Current state (agent `OmegaGateCFinish`, 2026-08-30) — gate C stage 3:
## the pre-registered PartLv2 repair applied; V1 FAILED at +0.90; FAILURE
## TO CERTIFY

Frozen campaign:
`campaigns/2026-08-30T22:58:26Z_OmegaGateCFinish_f5269412f0b1/`
(pre-statement `pre_statement_repairC.md` committed to git b41f384
BEFORE the first repair computation; scripts as run as byte-copies with
sha256; full diagnostic chain in manifest.md).

**Outcome: FAILURE TO CERTIFY, reported per the pre-registered FAIL
mode. The 21-dimensional nullspace question is UNCHANGED AND OPEN.** The
named repair was applied and VERIFIED to do exactly what it claimed:
M_low restored 0.14985 → 2.0942543887102634 (frozen 2.0942543887101843,
agreement 8e-17) and level-2 aggregation parity vs the primary source
`vxxz24_float.PartLv2.evaluate_post` holds at worst 6.9e-18 over all 1377
PartLv2 parts × 6 slots, with an embedded control proving the parity
harness detects the named defect (4.4e-2 ≫ 5e-13). Two FURTHER defects in
the predecessor's replay (never named by the repair note) were found and
fixed on the way: the part-level penalty's lam_sum entered with a flipped
sign (−lam_sum instead of +lam_sum; fixed, residual dropped 2.3 → 2.05e-6)
and the endpoint-slope numerator dropped R from the M' chain term
(Leibniz error, value part was −T/M ≈ −3.72 instead of Ω ≈ 2.37; fixed).

**V1 result as a number (pre-registered gate: |Om_raw − 2.3715538358350803|
< 5e-13):** the repaired point re-aggregation returns
Om_raw = 3.272425321778391, residual **+0.9008714859433109 — V1 FAILED**
(predecessor's defective value: 41.9962, residual +39.6247 — the repair
closed the level-2/M_low defect and part of the penalty defect but NOT
the whole aggregation gap). Per the pre-registration: STOP, no LP, no
decision run, nothing propagated to any index.

**Named next blocker (exact, from the parity evidence):** the glob-region
Lemma-1 eps term. The slope replay's max-over-shapes eps matches the
frozen protocol's value numerically (e.g. 3.076685e-6 at shape (8,0,0),
r=2), but the frozen interval tree's evaluate-post dot-path evaluates the
same term to ~1.78e-15-scale, so the two papers of the aggregation use
different eps computation paths; the divergence site in
`gate_c_slope_pass.py`'s glob penalty block (lines ~437-456) vs
`vxxz24_float`'s interval evaluate path is NOT identified — that single
block is the completion distance to a V1 re-run.

### Rule 7 sentence (mandatory)

The attempted swept set of THIS campaign was D = {delta : A·delta = 0
exactly, |delta|_inf ≤ 1e-7} around the region-0 glob dist block of the
VXXZ24 K100_2.37155181 released vector (A = the frozen 27x45 margin
matrix, rank 24 exact over Q; kernel basis exact rational) under the
transcribed 3-region single-p_comp program at max_level 3, q = 5 — and
even over D no decision was reached (V1 failed on a residual aggregation
defect). The statements here do NOT apply to the matrix multiplication
exponent record omega < 2.371177 (parameters unpublished as of
2026-08-30, never searched), to any point outside D, to any other rung's
released vector, to any other max_level/q regime, to multi-block or
cross-block direction classes, or to any construction differing from the
enclosed laser-method program (including the 2026 asymptotic-rank /
centroid-based improvements, outside every statement here). Even total
success in D caps 12.8x short of the published 2.37155181, so this route
can NEVER produce a record and no record claim may appear anywhere.

### What was NOT swept (explicit)

- No certified decision over D: neither PASS-a (no improving direction)
  nor PASS-b (an exhibited improving direction); LP/decision stages not
  run (V1 gate failed pre-LP by pre-registration).
- Regions 1, 2 glob dist blocks; multi-block/cross-block classes; the
  omega x single_mat_size Schoenage-line 2-coordinate class; all other
  parameter coordinates (held at p* points).
- The published 2.37155181 and the record 2.371177: not touched, not
  searched; nothing written here bears on the validity of the VXXZ24
  paper (no escalation triggered, none needed).
- No subdivision of D, no Monte Carlo, no probes outside D, no float SVD
  basis substitution (exact rational kernel basis re-verified and used).
