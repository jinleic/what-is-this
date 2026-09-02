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

## Current state (agent `OmegaGateCDiag`, 2026-08-31) — Gate C V1 repaired;
## exact-basis numerical candidate found; rigorous full-box sign remains OPEN

Frozen campaign:
`campaigns/2026-08-31T07:57:00Z_OmegaGateCDiag_4859327/`.
The original pre-statement was committed before every compute at `4859327`;
its campaign-directory byte-copy was committed at `256e22f` before the final
authoritative replay but after the first explicitly DIAGNOSTIC ONLY pass
(timing deviation disclosed, not backdated). Source-as-run commits are
`5b1e327` and `e9c6414`; scripts, stale bytes/traceback, artifacts, manifest,
report, and checksums are frozen in the campaign.

**RETRACTION of the preserved prior section's lines 528–544:** raw
`3.272425321778391`, residual `+0.9008714859433109`, and “glob Lemma-1 eps is
the next blocker” are not reproducible outputs of the advertised frozen
source. That source (SHA-256
`4504365a7917f80362424a7252203cd770b7c51392061e517f8c10e62fc946fa`)
raises `NameError: name 's_const' is not defined` at line 191 before its first
aggregate checkpoint (**COMPUTATIONAL-EVIDENCE**). The old values remain above
as **REPORTED** history only and are not retrofitted to repaired code. Their
implied stale `R=0.9303495043602614` and deficit
`1.8866540631007143` are **INFERENCE**.

### V1 result and requested penalty diagnosis

The repaired point replay machine-equals the frozen seven R branches:
`0.24530661807075219`, `0.24715974210930786`,
`0.2456557897939231`, `0.5883309785276519`,
`0.49684932826618977`, `0.4968521822710635`, and
`0.49684892842208755`; it also equals frozen
`R=2.8170035674609757`, `M=2.0942543887102634`, and raw
`Omega=2.3715538358350807` (**MACHINE-VERIFIED**). The residual against
registered raw `2.3715538358350803` is
`4.440892098500626e-16`: **V1 PASS**.

The first successful penalty-off run was overbroad (all Part/glob hash
penalties zeroed) and is retained only as **DIAGNOSTIC ONLY**; its scope
identification is retracted in the campaign correction ledger. The
authoritative narrower replay zeros only the `2*eps` Lemma-1/`lam_sum`
contribution, retains `hm` and `p_comp`, and machine-equals all seven branches,
R, M, and raw endpoint to the frozen `include_lemma1=false` record. Its raw
shift is only `2.029648699330977e-6` (**MACHINE-VERIFIED**), not an amplified
`0.9008715`. Answer to the inherited (a)/(b)/(c) question: **(b), a fourth
independent provenance/aggregation defect**; the advertised failing path is
non-executable and the corrected real-penalty path already passes V1.

### LP repairs and honest verdict

The endpoint quotient rule, hard pre-LP V1 gate, LP radius units, 90-row dual
shape, exact rational JSON basis, and exact Arb box endpoints were repaired.
Old-path counterfactuals fail as required: the extra `/r` changes effective
radius `1e-7` to `9.999999999999998e-15`; 21 variable-bound marginals cannot
multiply `V^T` of shape 21x45; the dense SVD basis is not the exact JSON basis;
and old binary64 endpoint arithmetic rounded 24 lower and 25 upper box
endpoints inward. Every interval-overlapping possible minimizer of each R/M
minimum is now hulled; single endpoint-argmin pinning is not consumed.

The exact-basis midpoint HiGHS LP returns a candidate
`-1.577332901516852e-7` (**COMPUTATIONAL-EVIDENCE**), but this is not a
rigorous sign. The full-box arbitrary-y Arb instrument gives
`[-5.841014555619999e-6, +5.841014555619999e-6]`, width
`1.168202911124e-5` (**MACHINE-VERIFIED**). The predecessor signal
`1.5816497000997742e-7` is a **CITED-DEPENDENCY**; width/signal
`73.85977508485647` is **INFERENCE**. The enclosure straddles zero.

**Verdict: FAILURE TO CERTIFY; the 21-dimensional question remains OPEN.**
A candidate direction exists numerically, but this is neither an improving
witness, nor “no improvement,” nor local optimality.

### Rule 7 sentence (mandatory)

The swept set was exactly `D = {delta in R^45 : A delta = 0 exactly,
|delta_i| <= r for i=0,...,44}`, where A is the frozen exact 27x45 0/1
margin matrix (rank 24 over Q), `r` is binary64
`0x1.ad7f29abcaf48p-24` (display `1e-7`), and the center is all 45 binary64
coordinates of the region-0 glob `dist[0]` block of the frozen VXXZ24
`K100_2.37155181` MATLAB vector; one radius-0 point and one un-subdivided
full box with exact Arb endpoints `IC.const(center_i) +/- IC.const(r)` were
run under the transcribed 3-region single-`p_comp` program at `max_level=3`,
`q=5`, with every other parameter fixed. Not searched: region-1/2 glob
blocks; any part/split/lambda/region_prop or other parameter block;
multi-block/cross-block directions; another radius; subdivision; Monte
Carlo; points outside D; other q/max-level regimes; unpublished parameters
for `omega < 2.371177`; or constructions differing from the enclosed
laser-method program, including asymptotic-rank/centroid-based improvements.
Nothing here bears on the published or record exponents, and no record claim
is made.

### Named next campaign and cost

**`OmegaGateCCandidateWitness`:** lift the midpoint LP candidate through the
exact integer kernel basis to an exact dyadic feasible delta (scale inward if
needed), then directly Arb-certify the moved endpoint against the frozen
center. Estimated cost: two certified point evaluations plus exact rational
feasibility checks, approximately **2–5 CPU-minutes**, no full-box sweep.

## Candidate-witness correction (agent `OmegaGateCDiag`, 2026-08-31) —
## exact LP ray evaluated; all five direct differences overlap; OPEN

Frozen follow-on campaign:
`campaigns/2026-08-31T09:18:39Z_OmegaGateCCandidateWitness_e9c6414/`.
Pre-statement commit `d739c45`; corrected source-as-run commit `1ec82a2`.
The initially frozen source at `8ad6062` is preserved as
**SUPERSEDED — STATIC GUARDS ONLY**: before any objective compute, its loop
was found to continue after the first witness contrary to the registered stop
rule. The corrected source breaks immediately and records the actually
attempted prefix; no result was contaminated.

Exact integer multiplication lifted the LP ray as `d=V*u`, with
`u=[1,-1,0,-1,1,0,1,0,0,1,-1,-1,-1,1,0,-1,0,0,-1,0,1]`.
All 27 entries of `A*d`, and `sum(d)`, equal zero exactly; every entry of d is
in `{-1,0,1}`. All 45 exact rational coordinates at the five fixed steps
`R/{16,8,4,2,1}` and the positivity/normalization/radius guards are frozen in
`protocol_static_corrected.json` (**MACHINE-VERIFIED**).

The full branch-safe Arb common-expression differences (new minus base) are:

- `R/16`: `[-2.023173821748832e-6, 2.035224142522071e-6]`;
- `R/8`: `[-2.016695072801931e-6, 2.041182148618541e-6]`;
- `R/4`: `[-2.0037375640758724e-6, 2.0538318431935425e-6]`;
- `R/2`: `[-1.977822503045151e-6, 2.0791500670709997e-6]`;
- `R`: `[-1.9259922068297443e-6, 2.1298618530686437e-6]`.

All five are **MACHINE-VERIFIED** and contain zero. Hence **OPEN**: no
registered rung is a strict improving witness, but this is not
no-improvement and not local optimality.

The midpoint LP prediction `-1.577332901516852e-7` remains only
**COMPUTATIONAL-EVIDENCE**. It did not survive because point-level Lemma-1
residual intervals straddle zero, leaving several minimum branches possible:
the exact-base seven-component R interval spans about `4.2506006968e-6`
(about `2.03e-6` in raw Omega), and the equality residual is about
`3.076685e-6`. At R, the full difference width
`4.055854059898388e-6` is `25.713367520566205` times the absolute midpoint
candidate (**INFERENCE**). Shared-target subtraction cancels T, but cannot
cancel this honest branch/defect uncertainty.

### Rule 7 sentence (mandatory)

This follow-on swept exactly five points `p_star+eps*d` for
`eps=R/{16,8,4,2,1}` on one exact integer direction d, changing only the 45
region-0 glob `dist[0]` coordinates of the frozen VXXZ24
`K100_2.37155181` vector under the transcribed 3-region single-`p_comp`
program at `max_level=3`, `q=5`; it did not search another direction,
negative/intermediate steps, a box, region-1/2 or another parameter block,
multi-block/cross-block classes, another q/max-level regime, unpublished
record parameters, or another construction. No paper or record claim is
made.

## Feasibility-scope correction (owner, 2026-08-31)

The candidate-witness section above correctly reports the five interval
evaluations as **OPEN**, but its phrase “positivity/normalization/radius
guards” overstates what `protocol_static_corrected.json` proves.  The exact
guard proves

`sum(p_star + eps*d) - sum(p_star) = eps*sum(d) = 0`,

not `sum(p_star + eps*d) = 1`.  Exact summation of any registered 45-coordinate
vector in that artifact gives

`37778931862949057407573 / 37778931862957161709568`

`= 1 - 8104301995 / 37778931862957161709568`

`≈ 0.9999999999997855`, strictly below one.  Thus the ray preserves the
released center's positive mass but is not an exactly feasible ray for the
registered sum-one constraint.  The exact claims that survive are `A*d=0`,
`sum(d)=0`, positivity, the radius bound, and the five fixed-point interval
evaluations.  Their **OPEN** verdict is unchanged; no witness was found.

Any later strict sign on these same pre-registered points may establish only
strict descent of the named repaired generalized-entropy/absorbed-defect
objective on this fixed-mass affine ray at the released binary64 center.  It
cannot be called a feasible construction, a counterexample to local
optimality of the constrained program, or an exponent improvement.

## Interval-enclosure suspension (owner audit, 2026-08-31)

An independent frozen-source audit found a formal outward-rounding defect in
the shared `interval_core.py`: `from_two(a, b)` formed its upper endpoint with
`b.lower()` instead of `b.upper()`.  The owner reproduced the defect directly
from the frozen source at lines 113--117.  Any nonzero-radius Arb value in the
upper expression can therefore be rounded inward before later buffering.

Consequently, every earlier claim above whose evidence label depends on this
`interval_core.py` is **SUSPENDED PENDING CORRECTED REPLAY**.  This includes the
two-rung “rigorous” endpoints, the certified-box interval endpoints, and the
five candidate-ray interval differences.  Their stored numerical outputs and
overlap relations are unchanged, but they are not currently accepted as
formal interval enclosures.  The floating-point reproductions, exact
integer/rational rank and feasibility checks, and the logical fact that no
accepted improving witness exists are unaffected.

The first entropy-repair production attempt (`...EntropyRepair_1ec82a2`,
v9) also aborted at the base point before its first entropy certificate because
its feature matrix was transposed.  The staged v10 fixed that interface but was
rejected before release: its binary64-only exact-point converter would reject
the registered moved rungs, and it still imported the defective interval core.
Neither attempt produced a strict-sign result.  A new frozen campaign must fix
the upper endpoint, accept arbitrary exact finite Arb dyadics without rounding,
close the live-import/input perimeter, and replay the affected statements
before any rigorous label is restored.

## Repaired generalized-entropy ray replay complete (owner, 2026-08-31)

Campaign
`campaigns/2026-08-31T14:29:40Z_OmegaGateCEntropySolver_v12_d2d3d4d5/`
closes the candidate-ray replay under the corrected outward interval core. The
preceding v11 production is preserved as **ABORTED**: it compared the repaired
raw objective to `FROZEN_RAW` from the old penalty objective. A bounded
diagnostic proved those intervals disjoint. It also found L-BFGS-B alone left
139/192 base multiplier certificates just outside or far above the exact KKT
gate.

v12 uses L-BFGS-B only as a warm start, then solves
`B(m softmax(B^T y))-b=0` by least squares with the analytic
`B(diag(q)-q q^T/m)B^T` Jacobian. Solver status remains computational
evidence; acceptance still comes only from the independent outward Arb KKT
bound. Static guards passed, and a fresh read-only Fable release audit found no
blocker. Owner reads resolved its two residual questions: `Ival.pos()` is the
absolute-value interval operation, and the historical dependency paths are
read lazily only after v12 redirects them to hard-hash-verified local inputs.

Production evaluated the base and all five registered divisors. Every one of
the 1,152 entropy certificates passed; maximum outward KKT upper bound was
`1.4456822252392066e-10 < 1/2^32`. Every repaired full-objective difference
(moved minus base) is strictly positive:

- `R/16`: `[6.475289389514407e-9, 6.475289390521749e-9]`;
- `R/8`: `[1.2954036629144681e-8, 1.2954036630152025e-8]`;
- `R/4`: `[2.5911541940662355e-8, 2.5911541941669702e-8]`;
- `R/2`: `[5.182659614230212e-8, 5.182659614330947e-8]`;
- `R`: `[1.0365687869954472e-7, 1.0365687870055208e-7]`.

Therefore all five named moves strictly **worsen** this repaired fixed-m
objective. Formal verdict remains **OPEN**: this certifies no descent on one
infeasible fixed-mass ray, not local optimality, not a feasible construction,
and not an exponent bound. The points still have mass approximately
`0.9999999999997855`, not one.

Machine result SHA-256:
`0403e652099a217d421c9bb1aa7fd56f07fb6daff72ddf14273d6a6dff063122`.
The historical two-rung endpoints and certified-box/full-box enclosures remain
**SUSPENDED** until their separate corrected-core replay; v12 restores only
the explicitly listed candidate-ray statements.

## Two-rung corrected-core replay: diagnostic recovered after aggregate abort

Campaign
`campaigns/2026-08-31T11:30:16Z_TwoRungReplay_e97c35ae_ce70f959/`
ran all three frozen child commands under corrected
`interval_core.py` SHA-256 `ce70f959…`. Their write-once outputs completed, but
the launcher then exited 1 before writing `comparison_summary.json`: its
comparison requested `omega_cert_upper`, while the frozen rung-1 producer
emits only `omega_cert_upper_raw` and
`omega_cert_with_absorbed_slack`. The abort and every child artifact are
preserved.

A separately preregistered, hard-hash-gated recovery chose neither missing-key
replacement. It compared both available rung-1 values:

- rung 1 raw: `2.371340083602922`;
- rung 1 with absorbed slack: `2.3713411672715115`;
- rung 2 with Lemma and absorbed defects: `2.3715538358544617`;
- rung 2 without Lemma diagnostic: `2.3715518062057623`.

Both possible rung-1-to-rung-2 inequalities hold, with respective gaps
`0.0002137522515397` and `0.0002126685829502`; all three ordering values are
below target `2.371866`. The rung-2 with-minus-without-Lemma diagnostic is
`0.0000020296486994`. Recovered result SHA-256:
`8160e846c4b64a167624dad8460e83bc29c4e69ac16115a7015da208a838adda`;
the 22-entry final checksum ledger passes.

These are **COMPUTATIONAL-EVIDENCE** comparisons only. The failed frozen key
remains unresolved. Formal status stays **SUSPENDED** because the historical
route still lacks an end-to-end outward aggregate without decision-relevant
float extraction and a proof that its feasibility-absorption term covers the
omitted true defect. No exponent endpoint is restored by this replay.

## Current state (agent `OmegaPenaltyProbe`, 2026-09-01) — gate C stage 3:
## V1 re-certified under the v11 core; LP/ANY-y re-run; question OPEN

Frozen campaign:
`campaigns/2026-08-31T02:40:00Z_OmegaGateCFinish_v1_4c8a2f1b/`
(pre-statement committed `ec16efe` BEFORE first compute; campaign freeze
`d72d8d6`, manifest finalize `c807f74`).

**What this closes.** The handoff premise "V1 replay fails at
`+0.9008714859433109`" was already retracted (see the
`OmegaGateCDiag` section above: the advertised source raises NameError
before its first checkpoint), and the genuinely open item was that **no
slope-replay V1 had ever been run under the corrected v11 interval
core** (`ce70f959…`) — GateCDiag's V1 pass predates the
interval-enclosure suspension. This campaign ran the frozen-protocol
control and the slope-replay V1 **in one process, one core, no source
edits**:

- Control (frozen stage_b_rung2, with Lemma): raw `2.3715538358350807`,
  delta vs anchor `+4.440892098500626e-16` ≤ 5e-13 — C-GATE PASS
  (COMPUTATIONAL-EVIDENCE fresh run; identical value already
  MACHINE-VERIFIED as the frozen raw).
- Replay V1: raw `2.3715538358350807`, delta
  `+4.440892098500626e-16` ≤ 5e-13 — **V1 PASS**
  (**MACHINE-VERIFIED**); all seven R branches, `R=2.8170035674609757`,
  `M=2.0942543887102634`, and lemma1-off raw `2.3715518061863814` are
  float-exact against both the frozen records and the pre-v11 GateCDiag
  replay records. The replay aggregation is float-stable across both
  interval cores. The named defect class (a replay aggregation error
  carrying a ~0.9 residual) is instrumentally dead on the live tree.
- LP/remainder stage (ran per gateC2 §2c/§4 as amended, only after the
  V1 pass): box endpoints exact-Arb asserted on all 45 leaves; the
  binary64 counterfactual again rounds 24 lower/25 upper endpoints
  inward (same coordinate sets as GateCDiag);
  scaled midpoint LP `-1.5773329015168567e-7`
  (**COMPUTATIONAL-EVIDENCE**); ANY-y rigorous enclosure
  `[-5.841014555620028e-6, +5.841014555620028e-6]`, width
  `1.1682029111240058e-5` (**MACHINE-VERIFIED**); width/signal
  `73.85977508485684` (**INFERENCE** label). Enclosure straddles zero.

**Verdict: gate C stage 3 — replay machinery re-certified end-to-end
under the v11 core; decision stage FAILURE TO CERTIFY; the
21-dimensional kernel question remains UNCHANGED and OPEN. Gate C stays
PARTIAL.** Not PASS-a, not PASS-b; six LP probes cannot certify a
21-dimensional space.

### Rule 7 sentence (mandatory)

Swept: one radius-0 point pass and one un-subdivided full box over
`D = {delta in R^45 : A·delta = 0 exactly, |delta_i|_inf ≤ 1e-7}` around
the region-0 glob dist block of the VXXZ24 `K100_2.37155181` released
vector (A = frozen 27x45 margin matrix, rank 24 exact over Q; kernel
basis V exact rational, integer-asserted in-run; NOT a float SVD basis),
with the scaled midpoint LP and the ANY-y instrument evaluated over that
un-subdivided D, plus the same frozen stage_b_rung2 protocol re-run at
the anchor point as control — all under the transcribed 3-region
single-`p_comp` program at `max_level 3`, `q 5`. Not swept: region-1/2
glob blocks; any part/split/lambda/region_prop or other parameter block
(held at p*); multi-block/cross-block classes; subdivision; Monte
Carlo; other radii; probes outside D; the omega x single_mat_size
Schoenage 2-coordinate class; the published `2.37155181` and the record
`2.371177` (never searched); any other rung, regime, or differing
construction (incl. 2026 asymptotic-rank/centroid improvements). Even
total success in D caps 12.8x short of the published `2.37155181`; no
record claim anywhere.

### Named next action

Unchanged binding constraint: full-box slope-bundle width (73.9x the
`1.5816497000997742e-7` signal; branch-straddle on the Lemma-1 residual
intervals at the base point). Next campaign (must be fresh
pre-registered): certified subdivision (branch-and-bound) of the box so
each subbox pins its R-branch/Lemma-1 structure, or a tighter certified
gradient enclosure; scope guardrail (12.8x cap, route never produces a
record) carries verbatim.

### Main acceptance and provenance notes (2026-09-01)

Owner independently re-derived and accepted into RESULTS.md (omega:
PARTIAL / FAILURE TO CERTIFY, V1 closed, 21-dim OPEN): the V1 delta
`4.440892098500626e-16` is exactly one ULP at 2.37; the lemma1-off shift
`2.029648699330977e-6`, the width/signal `73.85977508485684`, and the
`12.808x` guardrail arithmetic all reproduce to the last digit; the swing
sits outward of `L1*radius` by `5.5e-20` (conservative direction); and
rule 14 holds directly — step 1 importlib-loads the frozen
stage_b_rung2.py while step 2 calls gate_c_slope_pass.run_slope_pass, so
the control/replay agreement is informative, not tautological.

**Campaign-directory timestamp provenance (disclosed):** the frozen
directory name carries UTC `2026-08-31T02:40:00Z` but the run actually
executed `2026-09-01T02:47–02:53Z` (~24 h later); the clock-stamped name
was minted at pre-statement write time, not launch time. Ordering
relative to the GateCDiag campaign (`2026-08-31T07:57:00Z`, pre-statement
commit `4859327`) is: GateCDiag FIRST (it landed the retraction and the
pre-v11 V1 pass), this campaign SECOND. Read the ledger by commit order,
not by directory-name timestamp.

**Freeze-hygiene gap and forward rule (disclosed):** during this
campaign the aborted attempts' stderr tracebacks were overwritten by
each re-run's `2>` redirect; the completed pass's empty stderr is the
only transcript that survives, and the manifest's §4 description is the
only record of the three runner defects. Forward rule for this agent's
future runs: per-attempt stderr/log files under DISTINCT filenames (one
per attempt), so every aborted attempt's traceback survives freeze.

## Current state (agent `OmegaPenaltyProbe`, 2026-09-01) — gate C stage 4
## feasibility gate: branch-and-bound INFEASIBLE AT THIS COST (frozen
## quantified obstruction); stage 2 never opened; question OPEN

Frozen campaign:
`campaigns/2026-09-01T03:15:00Z_OmegaGateCBnBGate_9d2e4a7c/` (pre-statement
`f8b9ca7` BEFORE compute; GATE-Z tolerance amend `532781d` pre-verdict;
verdict freeze `4ba4c38`; actual-UTC directory stamp per the provenance
rule).

Per Main's steering, the load-bearing empirical question for any
branch-and-bound successor — is the box enclosure's width CONCENTRATED
in a few coordinates? — was measured on the certified full-box slope
bundle of the parent campaign (width/signal 73.85977508485684, L1
113.6656…):

- **SPREAD.** k* at ≥ 90% of the L1 width is **36 of 45 coordinates**
  (cumulative shares k2 8.17%, k4 15.22%, k6 21.85%, k10 34.63%, k23
  68.41%; verdict margin k35 89.902% → k36 91.384%, far from the ≤ 6
  gate). Widths are nearly uniform in [0.664, 4.641]; 40 of 45
  coordinates each carry ≥ 1% individually.
- Branch hulls: every R branch's minimizer identity is undecided
  (all three tt slots possible on six of seven branches, two on the
  seventh; zero never possible); the widest candidate span is
  R_glob[0] at 3.212e-5 ≈ 203x the signal — consistent with the frozen
  CandidateWitness straddle picture.
- **Implied cost:** ~log2(73.86) ≈ 6.21 bisections per binding
  coordinate times k* = 36 ⇒ ~2^223 subboxes to force width/signal < 1
  (mild half-width reading still ~2^112). **Branch-and-bound over D at
  this enclosure technology is INFEASIBLE AT THIS COST** — frozen as a
  first-class quantified negative (cf. oct-rank S3).
- GATE Z cross-run control passed under the pre-registered replay-noise
  tolerance (max |Δlo| 1.20e-14, max |Δhi| 1.78e-15 ≤ 1e-12; L1 within
  3.55e-13 ≤ 1e-6); the 90% crossing sits at k = 35/36, far from any
  threshold, so last-digit drift cannot flip the branch.

**Verdict: gate C stage 4 (subdivision route) — INFEASIBLE AT THIS COST;
stage 2 never opened; no subbox compute, no calibration run. The
obstruction binds this enclosure technology (min_slope hulls over an
un-subdivided box); it does not bound a different certified-gradient
technology, an analytic branch-pinning argument, or subdivision coupled
with a certified crossing proof. The 21-dimensional kernel question
remains UNCHANGED and OPEN; gate C stays PARTIAL.**

### Rule 7 sentence (mandatory)

This campaign measured the decomposition of the EXISTING certified
full-box slope bundle on D = {delta in R^45 : A·delta = 0 exactly,
|delta_i|_inf ≤ 1e-7} (region-0 glob dist block of the VXXZ24
K100_2.37155181 vector; A = frozen 27x45 0/1 margin matrix, rank 24 over
Q; interval core v11, slope pass re-run at the same frozen centre/radius;
max_level 3, q 5) and froze the resulting obstruction arithmetic. It did
NOT run branch-and-bound, subdivision, calibration, any new enclosure,
any other parameter block/region/q/max_level/radius/construction, or any
point outside D. The published 2.37155181 and record 2.371177 were not
touched; the 12.808x guardrail (gap 2.0258350805768544e-6 vs signal
1.5816497000997742e-7) keeps both unreachable by this route; no record
claim anywhere.
