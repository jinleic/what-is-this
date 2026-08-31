# Gate C pre-statement — certified local optimality, VXXZ24 rung — PRE-REGISTERED

Campaign: `2026-08-30T11:24:41Z_17bdaeb6_099b1224d603`. Written and committed
BEFORE the first branch-and-bound call (write-first protocol; file timestamp
predates all B&B outputs in this directory).

## 0. Centering (Main steering 2026-08-30, supersedes assignment default)

Gate C is centered on the **VXXZ24 rung** (`omega <= 2.37155181`), released
vector osf.io/7wgh2 `data/K100_2.37155181.mat`
(sha256 `df75ae3acaa5b1388bd2f17cc e266aec169d874c9fea1eb0722c3e25871da93d`,
6759 float64 params), because our transcription `src/vxxz24_float.py` is
already built and parameter-count-verified exact (6759 == 6759), giving a
known-good stage-(a) anchor. The Alman25 rung would be the repeat target if
runway allows.

**Scope honesty (mandatory).** The current record omega < 2.371177
(arXiv:2608.16884) has NO published parameter vector (artifact audit re-run
by Main 2026-08-30: google-deepmind/alphaevolve_results holds no parameter
table; GitHub searches return 0). There is no "reported optimum" of the
record to put a box around. Any gate-C result here is a statement about the
**rung-2 predecessor optimum**, is a repo-first rigorous local-optimality
statement, and must NEVER be written as a statement about the record. The
rule-7 sentence below says this explicitly.

## 1. The reported optimum (centre), sourced by file and line

- Centre vector p* := the 6759-dim float64 vector stored in
  `scratch/rmmcode/data/K100_2.37155181.mat`, key `params`, loaded via
  `scipy.io.loadmat` and flattened (`np.asarray(...).flatten()`), exactly as
  in `src/stage_b_rung2.py:68-70`. Source: osf.io/7wgh2 (VXXZ24 release).
- The released harness's own optimum readout:
  `omega_id` param value = `pm.cur_x[pm.start[ws.omega_id]]` — file stores
  2.3715518061833021 (campaign notes of run 05:20:00Z; re-verified below
  before any B&B call). The filename constant 2.37155181 is the display
  rounding.
- No re-derivation of the centre is attempted; p* is a fixed bit-exact
  float64 array, hashed before use.

## 2. The box (per-coordinate, exact, pre-committed)

- Rung-1/rung-2 coupling convention: perturbation is applied to ALL float64
  parameters of p*, i.e.
  `p(x) = p* + r * x` elementwise (numpy semantics, no exception list).
- Per-coordinate radius: **r = 1e-7** (float64 magnitude of the vector's
  smallest meaningful digits; 6759 coordinates; the box is
  `[p*_i - 1e-7, p*_i + 1e-7]` in every coordinate — a hypercube in R^6759).
- Chosen BEFORE search as the smallest radius that still plausibly
  dominates the float64 evaluation granularity:
  - SMALLEST use: any coordinate exactly 1 ulp from p*_i is inside the box
    (1 ulp of p*_i <= 2^-52 * |p*_i| <= 4.4e-16 for the range of
    coordinates involved, 60+ orders below r=1e-7). Stage-(a) sampled
    directions are evaluated at both boundary corners +-r * NaN-injection:
    `pm.set_value(p* + r * np.sign(direction) * 1.0)` per coordinate.
  - LARGEST use: r is well above the float64 evaluation granularity so the
    bound cannot be an artifact of interval under-width, yet well below the
    1.1e-9 refine tolerance? NO — r = 1e-7 is 100x the refine tolerance
    1.1e-9: this is a box around the ENGULAR neighbourhood, not the
    tolerance neighbourhood. Rationale: the Lemma-1 residual of the shipped
    multipliers is ~3.08e-6 (glob r2), so the certified-enclosure cost at
    p* itself is larger than 1.1e-9; the box must be large enough that the
    bound degrades observably if the transcribed optimum hides a descent
    direction larger than the residual scale. 1e-7 is roughly 3% of the
    residual scale and orders of magnitude above float64 granularity —
    large enough to be informative, small enough for the relaxed interval
    widths to stay manageable.
- **Per-coordinate, explicit:** radius is uniform r = 1e-7 in every one of
  the 6759 coordinates. No coordinate is shrunk; no coordinate is
  re-expanded. The swept set is exactly
  `B = {p* + r*x : x in [-1,1]^6759}, r = 1e-7` in float64 embedding.

## 3. Branching rule (pre-committed)

- BoxType `Box(offsets: 6759-dim float64 array of half-widths per coordinate,
  signature)`. Initial box: all-coordinate widths 2e-7 (r = 1e-7 each side).
- Valuation: certified Omega upper = `float(enclosure_centre(box))` from the
  SAME interval code path (stage_b_rung2-style evaluate -> cert aggregation)
  evaluated at the box's midpoint, with `eps.width = r` propagated as a
  widening of every Ival endpoint:
  - `stage(b-reuse)`: at a box with centre c, we form the interval tree by
    `pm.set_value(c)` then assert every Ival endpoint is widened by the
    L∞ perturbation weight pre-computed per group as `w_g = r * L1_sensitivity(g)`
    where the sensitivity weights are TRANSCRIBED from the interval tree's
    own evaluation gradients (NOT re-run: cache from the p* evaluation).
    This is the standard "interval enclosure at c + [dep] * r" relaxation.
  - Concretely for each registered parameter group g with coordinates
    c_g, the Ival fed to the interval path becomes
    `Ival(c_g - w_g, c_g + w_g)` (endpoints are arb, outward-rounded),
    where `w_g = r` (the SAME uniform radius — a deliberate over-relaxation
    for diagonal boxes: anything smaller than per-coordinate r would be a
    first-order claim, and we pre-commit to the honest over-approximation).
- Branching: pick the coordinate with LARGEST |offset_i| (tie: smallest
  index), split at the midpoint into two half-boxes; recurse depth-first.
- **Depth cap: 3 levels of branching** (initial box + at most 3 splits
  along the walk before a subbox must be decided, not split again).
- **Box budget: 63 subboxes total across the whole run** (2^6 - 1 upper
  bound on a 6-level complete walk would be 63; we cap the total number of
  pushed boxes at 63 regardless of tree shape). Counted and logged to the
  checkpoint file as each box is opened — no budget overrun possible.
- **Wall-clock cap: 240 minutes** from campaign start (pre-registered).
  On expiry the run stops and reports the last certified aggregate.

## 4. Stage ordering (README-mandated)

1. **Anchors FIRST.** Before any box other than B0 is touched, the B&B
   evaluator must reproduce, to 4+ significant digits of both endpoints:
   - omega <= 2.3715538358544617 (rung 2, WITH Lemma-1 charge, at p*),
   - the no-Lemma-1 diagnostic readout 2.3715518062057623 (not a pass/fail
     gate — recorded only, since it is the run without the charge);
   - rung 1 anchor 2.3713400836689 requires the Alman25 code path
     (`stage_b_rung1_certified.py` semantics); it is RE-EXECUTED here only
     as a smoke check that the shared machinery is intact before B&B, not
     re-derived. Anchor failure => STOP, report the failure, no B&B.
2. **Stage (a) float B&B** — same box, same branching, float64 aggregation
   via `vxxz24_float.max_violation` + the float R/M/Omega reconstruction in
   `stage_b_rung2.py` (float path). This is the cheap screen that must come
   out CLEAN (enclosure closes at or below the p* level, at worst slightly
   above due to float noise <= 1e-11).
3. **Stage (b) Arb B&B** — the SAME code path with Ivals (interval_core).
   Every aggregate reported is an interval; widths are first-class numbers
   in the result JSON.

An enclosure failure in stage (b) after a clean stage (a) is evidence about
the bound. An enclosure failure in stage (b) WITHOUT a clean stage (a) is
evidence about the transcription and must be reported as such.

## 5. Part 1 first — the 2.0258545e-6 gap decomposition (BEFORE B&B)

Re-run the identical int code path (`stage_b_rung2.run(include_lemma1=True)`
via `interval_core.prec()` with `MID` swept over multiple precisions; both
arms). Report gap = cert_upper - 2.37155181 at each precision together with
R_sum widths and M widths. Gap shrink with precision => outward-rounding
accumulation; gap flat => structural (transcription or VXXZ24's own
rounding), and Main is escalated before any conclusion. This decomposition
RUNS BEFORE the B&B and calibrates how much room a genuine improvement has
to beat.

## 6. Evidence labels (pre-committed)

- Anchors: MACHINE-VERIFIED (interval arithmetic, per-block containment
  asserted, deterministic re-run).
- Precision sweep: MACHINE-VERIFIED widths; gap-vs-precision CLAUIM is
  MACHINE-VERIFIED only if every precision point is re-run from the frozen
  scripts. COMPUTATIONAL-EVIDENCE if any number is copied from the earlier
  campaign.
- B&B float leg: COMPUTATIONAL-EVIDENCE at best; never labelled stronger.
- B&B Arb leg: MACHINE-VERIFIED if and only if per-block containment
  asserts fired for every evaluated box; any box failing containment is
  reported as a MISS, never dropped.
- Every box evaluated appears in the report, including boxes whose
  enclosure is too weak to close (misses). No survivor-only reporting.

## 7. Rule 7 (pre-committed, exact sentence)

A certified cap or improvement applies to the exact swept box
`B = p* + [-1e-7, +1e-7]^6759` around the VXXZ24 K100_2.37155181 released
parameter vector under the transcribed 3-region single-p_comp program at
max_level 3 and q = 5 — it is NOT a claim about the matrix multiplication
exponent record omega < 2.371177, which was not and could not be searched
because its parameter vector is unpublished as of 2026-08-30, and it is
NOT a claim about any point outside the stated box, any other released
rung vector, or any other max_level/q regime.

## 8. Escalation triggers (pre-committed)

- Any box certifies omega strictly below 2.37155181 (an improvement over
  the published rung-2 rounding): freeze, hand to Main, write nothing.
- Gap decomposition shows the 2.026e-6 gap NOT shrinking with precision:
  Message Main with the sweep table, stop that thread, do not conclude
  anything about the paper.
- Any per-block containment assert fires during B&B: stop, record as MISS,
  do not tighten semantics post hoc.
