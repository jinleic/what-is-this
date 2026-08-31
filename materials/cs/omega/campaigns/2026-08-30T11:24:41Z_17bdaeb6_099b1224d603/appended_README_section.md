
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
- **The reported point is NOT the certified optimum of the swept
  family:** three of the fourteen signed boxes certify strictly below
  B0 — region_prop+ (2.3715530208533490, −8.15e-7 vs p*), glob_dist_2+
  (2.3715534164727120, −4.39e-7), glob_dist_0+ (2.3715537111130380,
  −1.44e-7). Certified endpoint improvements over the p* enclosure
  within the pre-registered box family, NOT improvements of the
  published bound (all stay above 2.37155181); recorded here only as
  gate-C evidence per pre-statement §8.
- The winning direction region_prop+ is the symmetry coordinate: the
  released vector has region_prop X/Y/Z-symmetric, and the certified
  enclosure prefers a tilted region_prop at the 1e-7 scale.

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
