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
