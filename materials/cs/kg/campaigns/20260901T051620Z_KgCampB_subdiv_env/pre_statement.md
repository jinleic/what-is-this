# Pre-statement — Campaign B: pre-registered mq = t² panel subdivision, validated
# on Campaign A's identified failing cell (tile 3 of band [4.083, 6.0])

Agent `KgCampB`, UTC 2026-09-01. **COMPUTE PENDING — this file is committed
before any numerical evaluation of this campaign.** The only prior files of
this directory are the skeleton commit `01bd720` and inherited instrument
copies whose checksums are recorded in `checksums.sha256` at freeze time.

Predecessor: Campaign A
`campaigns/20260831T231500Z_KgEvenBoxRerun_campA_corr_env` (frozen at commit
`8d2cb26`): band [4.083, 6.0] under the corrected 2-D even-box envelope —
tile 1 PASS, tile 2 PASS, tile 3 OPEN at a disk-ring leaf, band **FAILURE TO
CERTIFY** at 3/49 by its pre-registered stop-at-first-failing-tile rule;
tiles 4–49 NOT RUN. The kgcloseB3 retraction (rule 5) STANDS; this campaign
earns no PASS retroactively and bears on the paper's correctness in neither
direction.

## 0. What is re-used, byte-identical (declared now)

From Campaign A's frozen `code/` (sha256 at `checksums.sha256` of that
campaign): `gate_corr_env.py` (envelope cores `e_panel_range`, `dist0_box`,
`box_intersects_disk`, `root_boxes`, `split_box`, `can_split`, tiling
`c_tiles` ratio 1.008, ladder 4096×4 capped once at 131072, floor
`BOX_FLOOR_Q = 1/1024`, `BOX_DEPTH_MAX = 40`, budget semantics,
`strict_margin_lower`, S_MAX = 8, precision 256) and `startup_controls.py`
(the full six-part battery: kernel admission incl. the cap-at-1
counterfactual, 2-D vs 3-D odd-radius counterfactual, K_o/K_e identities,
cited-tail anchor with the byte-verified d3h sha 80e945589b…, margin-predicate
trio, envelope-admissibility on p*). Both are copied verbatim into this
campaign's `code/`; no byte of either is edited. Shared dependency
`20260831T082425Z_kg_direct_d4_restart/code/core.py` is imported, not copied.
python-flint landmine discipline (owner rule): every interval is built in
midpoint/radius or outward endpooint form explicitly; `arb(lo, hi)` is never
used as an interval.

## 1. The pre-registered mechanism being validated (Phase 1)

**Lineage (rule 16 provenance).** kgcloseB3 `pre_statement.md` Addendum 2
pre-registered, as its named next-campaign first change: *"the fix that would
finish it is panel-subdivision AT THE CROSSING"* — the crossing of the fold
case boundary M·q = t² (M = |e|+|o|, q = ||e|−|o||, t = c_L·s). It was
deliberately NOT adopted mid-flight in kgcloseB3 and not present in Campaign
A's instrument (whose per-panel bound applies the second hinge
unconditionally). Adopting it in this new campaign is the pre-registered
legitimate route.

**Definition of the subdivided envelope V(B, c_L; n).** Panels: uniform n
panels on [0, S_MAX=8], exact fmpq endpoints, skip rule unchanged (panel
skipped iff √K(b) ≤ c_L·a). With odd2 := (1 − l²)₊ outward (l = dist(0, B)
2-D, min'd at 1) and W_P := √odd2·√K_o(b) as in Campaign A:

1. **Certified second-hinge kill (per panel/subpanel).** Let [E_lo, E_hi] be
   the corner-exact even range (A's `e_panel_range`, unchanged), E⁺ its
   outward |·| upper, and E_min the outward distance of [E_lo, E_hi] from 0
   (0 if the interval contains 0). For any compatible unit cubic, pointwise
   on the panel: q = ||e|−|o|| = |e²−o²|/(|e|+|o|) ≤ max(E⁺², W²)/E_min
   whenever E_min > 0. Kill condition: **max(E⁺², W²) ≤ t_L·E_min**
   (t_L = c_L·a); when it certifies, the second hinge is 0 on that
   subpanel: G^V = (E⁺ + W − t_L)₊. When the even interval contains 0 or the
   inequality fails, G^V = (E⁺ + W − t_L)₊ + (|E⁺ − W| − t_L)₊ (A's bound —
   always valid by Lemma D.4 monotonicity).
2. **Panel subdivision at the fold-case crossing.** The mq = t² boundary is
   the set {s ≥ 0 : odd2·K_o(s) = (c_L·s)²}, i.e. s > 0 roots of
   odd2·(5/2 − s² + s⁴/6) = c_L². For s ∈ [√3, ∞) the left side is strictly
   increasing, so all roots are bracketed outward in Arb by monotone bisection
   from a doubled-bracket start (abort if the bracket fails to double past the
   root). Every root inside a panel splits that panel at its outward bracket
   endpoints (no gap, no overlap; endpoints exact fmpq where the bracket has
   converged exactly, else outward arb endpoints used only as partition
   bounds — the Gaussian weight is evaluated with outward `gauss_cdf` on each
   subpanel, so the outward split cannot lose mass). The case rule of item 1
   is applied per subpanel.
3. **Soundness insurance (assert, abort on failure).** V ≤ U pointwise
   capability: for every panel, the per-subpanel weighted sum cannot exceed
   A's whole-panel bound (E⁺_sub ≤ E⁺_panel, t_L_sub ≥ t_L_panel, kill only
   removes), verified numerically at fixed tolerance 2^-200 on planted cells;
   if any V evaluation exceeds the corresponding U evaluation, ABORT.

**Exact-reduction property (declared so the measurement is interpretable):**
on a cell whose crossing s* lies outside [0, 8], item 2 inserts no split and
V differs from A's envelope only by item 1's kill term. The mechanism's
recovery is then exactly the kill's certified mass, measured below.

## 2. Frozen caps, knobs, budget (unchanged from A; rule 16)

Tile: the pair (c_L, c_R) = (4.1485893120, 4.181778026496) — tile 3 exactly
as A left it, generated by the byte-shared `c_tiles('4.083', '6.0')` ratio
1.008 tiling (A's frozen values; the assignment's displayed 4.1817780265 is
this number rounded). Root cover 132, floor 1/1024, depth 40, panel ladder
4096→16384→65536→131072, cap 131072, precision 256, S_MAX 8, closed tail.
Budget: 260,000 envelope evaluations for Phase 1. Stop rule: identical to A's
(open leaf at panel cap and floor depth ⇒ tile OPEN ⇒ Phase-1 FAIL). DFS
order, split rule, and all bookkeeping are byte-shared with A's
`certify_tile`.

## 3. PHASE-1 VERDICT — threshold fixed in advance

Re-run tile 3 from the 132 roots under V. **PHASE-1 PASS ⇔ every
disk-intersecting leaf, including A's frozen frontier cell
a0 ∈ [−0.830729166666666666666667, −0.829427083333333333333333],
a2 ∈ [+0.558593750000000000000000, +0.559895833333333333333333],
certifies with margin ball d(c_R).lower − U^V.upper having strictly positive
lower end (> 0, not merely straddling).** PHASE-1 FAIL ⇔ any leaf OPEN.

- On PASS: proceed to Phase 2 (§5). On FAIL: STOP — Phase 2 does not run —
  and freeze the quantified obstruction per §4 and §6.
- Reported either way: the new open-leaf margin (Arb ball, rule 15 exact
  subintervals), the recovered fraction of A's −1.7629e-5 deficit
  (certified Arb difference of frontier-cell margins), certified-leaf count
  out of total, achieved depth, and per-rung panel histogram.

## 4. Slack-decomposition battery (frozen now; runs after Phase 1 either way)

All on A's frozen frontier cell only, all outward Arb 256-bit, all certified
differences are rule-15 named balls:

- **D1 (baseline recompute).** A's envelope U on the cell at cap 131072 —
  must reproduce A's frozen open-margin −1.7629e-5 to width, else STOP and
  escalate (reconciliation before anything else).
- **D2 (Phase-1 mechanism).** V at cap 131072. Decomposition entry
  "fold/second-hinge": D1.upper − D2.upper (≥ 0 by the soundness insurance).
- **D3 (panel-sup marginal).** U at caps 2^19 and 2^21 (single cell; ~10 s
  and ~40 s). Entry "panel-sup": D1.upper − U(2^21).upper, plus U(2^19) and
  U(2^21) values; a measured per-quadrupling improvement series (orientation
  for the "implied panel count", clearly labelled COMPUTATIONAL-EVIDENCE
  where extrapolated).
- **D4 (|e|-hinge coefficient-space diagnostic — NOT a Phase-1 verdict
  mechanism).** The even sup evaluated over the **feasible set
  B ∩ closed unit disk** instead of B (the sliver: sup of the affine
  e(s) = a0 + a2·ψ2(s) over box∩disk — attained at box corners inside the
  disk, box-edge/circle intersection points, or the disk-corner direction
  (±1, ±ψ2)/√(1+ψ2²) if feasible; all candidates evaluated outward). Entry
  "|e|-hinge (coefficient space)": U(2^19).upper − U^disk(2^19).upper.
  This term measures, certified, the excess A's box-corner sup carries over
  the box's disk-feasible part. **Firewall: a nonnegative D4 is a diagnostic
  input to the obstruction freeze only; it cannot flip the Phase-1 verdict,
  which is defined solely by V (§3).**
- **D5 (tail).** Tail_8 closed-form ball (~3.6e-13): identical in every
  variant, so the decomposition's tail differential is exactly 0; value
  reported.
- **D6 (t_L/within-panel skew).** Contained in D3's panel-sup series;
  reported as the 2^19→2^21 delta.

Every decomposition value is an outward Arb ball; the freeze lists them with
their sum and the residual (frontier margin at best measured knob state), so
the owner can re-add them from the JSON.

## 5. PHASE 2 (only on Phase-1 PASS) — band continuation protocol

Tiles 4–49 of the same ratio-1.008 tiling, continued from tile 3's right edge
under V, A's knobs and stop-at-first-failing-tile rule, budget 260,000
band-global envelope evaluations, per-tile artifacts `inproc_tile_NNN.json`,
final rule-15 table: every tile's c-pair and worst certified margin ball; any
failing tile reported with its frontier exactly as A did. A completed band
would be reported as a PASS only if all 49 tiles certify; anything else is
FAILURE TO CERTIFY at the first failing tile with the exact frontier. The
refutation gate is inherited verbatim: a certified J > d interval on any leaf
stops the campaign and escalates to Main.

## 6. Obstruction freeze (on Phase-1 FAIL)

`obstruction_tile3.md` + `obstruction_tile3.json`, frozen with checksums,
containing: residual deficit (Arb ball); the §4 decomposition table with
exact certified subinterval values; the certified statement of the floor
cell's refinement exhaustion (depth 40/floor 1/1024 reached — no box
refinement available within A's knobs); the measured/implied panel-count
analysis (D3 series; if the measured improvement per quadrupling is bounded
by a geometric tail below the deficit, the NONE-finite-cap conclusion is
filed as COMPUTATIONAL-EVIDENCE with the certified D4 mechanism recovery as
the machine-verified named route — the "implied depth/panel count to close"
field states explicitly that no panel/depth knob closes the cell and the
closing mechanism must be coefficient-space); the named next campaign.
Rule-14 semantics unchanged: a negative envelope margin is a failure of the
upper envelope at the caps, not a lower bound on J−d and not a paper
refutation.

## 7. Scope (rule 7, whole campaign)

Swept here: exactly the c-tile (4.1485893120, 4.181778026496) over the full
2-D even-coefficient disk (132-root cover, all compatible odd parts) with the
subdivided envelope V at A's frozen knobs, plus the §4 single-cell
decomposition battery on the one named frontier cell; and, only if Phase 1
passes, tiles 4–49 of [4.083, 6.0] under V. NOT swept: tiles 1–2 re-runs (A's
PASSes stand, not re-litigated), any c outside [4.083, 6.0], the never-run
bands [1.0,1.3], [1.45,1.75], [1.75,3.5] (Campaign C, separate pre-statement),
[6.0,12.0], the paper's C1/splice/C3 certificates, degree > 3 families, and
any coefficient-space envelope change (D4 is a diagnostic, never a verdict
mechanism in this campaign). Startup controls sweep only their six planted
assert families.

## 8. Evidence labels

MACHINE-VERIFIED: margin/sign/interval claims from outward Arb endpoints at
the stated domains, scripts frozen in `code/`. COMPUTATIONAL-EVIDENCE:
timings, wall/CPU duty (resource.getrusage/time.process_time only, rule 17e),
and any extrapolated series. REPORTED: paper Table-1 numbers. FAILURE TO
CERTIFY is a result and is reported honestly; no PASS is claimed anywhere
unless every leaf of every in-scope tile carries a strictly positive
certified margin ball.

## CAMPAIGN B RESULT (frozen; checksums in `checksums_addendum.sha256`)

**Phase 1: FAIL — tile 3 remains OPEN under the pre-registered mq = t² panel
subdivision; Phase 2 was NOT run.** The quantified obstruction is frozen in
`obstruction_tile3.md` / `obstruction_tile3.json`.

Startup controls: 16/16 PASS (`logs/startup_controls.json`, byte-frozen; the
full A battery re-ran green in this campaign, including the kernel-admission
assert with e(0) = √(3/2) admitted at slack 1.87e-96, the cap-at-1
counterfactual (certified deficit 0.224744871391589049…), the 3-D-inflation
counterfactual (1.65423992 vs corner-exact 1.32426407), the byte-verified
d3h certificate sha 80e945589b…, and the three-way margin-predicate trio).

Phase-1 instrument runs (MACHINE-VERIFIED, Arb 256-bit outward):

| run | artifact | result |
|---|---|---|
| 1A frontier cell under V | `logs/phase1_frontier_cell.json` | OPEN at panel cap 131072; margin −1.76291787028287693244667476546e-5 ± 4.20e-35; recovery vs A exactly 0 |
| 1B full tile 3 under V | `logs/phase1_tile3.json` (stdout `logs/phase1_tile3_stdout.log`) | OPEN; 819 evals, 95 leaves certified at stop; same frontier cell/depth 14/panel cap as A; tile certified_margin_min +6.85420926243409215516837143105e-7 |
| soundness insurance | inside both runs | V ≤ U asserted on 4 planted cells × 3 ladder rungs before any verdict |
| decomposition battery | `logs/decomposition_tile3.json` | see below |

Deficit accounting on the frontier cell (all rule-15 Arb balls):

* fold/second-hinge (D2, the Phase-1 mechanism): recovery **exactly 0** —
  the mq = t² crossing sits at s\* = 45.703407790989350730 (certified
  bracket), outside [0, 8]; V ≡ A's envelope there.
* panel-sup (D3): +2.20922696204705e-5 (2^17→2^19), +5.52295164708941e-6
  (2^19→2^21). The panel ladder cap, not the envelope case split, binds:
  at cap 2^19 the cell's margin turns +4.46309091764170239642742026145e-6.
* |e|-hinge coefficient-space excess (D4, diagnostic only, firewall held):
  the disk-feasible even sup is smaller by 4.94055685681106e-4 = 28.0× the
  deficit; margin would be +4.76426506978277e-4. This requires a JOINT
  (even, odd) disk-feasibility re-derivation of the envelope to be sound —
  named as Campaign C's first change; NOT adopted here.
* tail (D5): 3.6353172540086755822e-13, identical in every variant.

**Phase-1 verdict per §3: FAIL** (open leaf at panel cap; margin ball's
lower end < 0). Phase 2 does not run per the pre-registered threshold. Band
[4.083, 6.0] remains FAILURE TO CERTIFY at tile 3/49; A's two tile PASSes
stand untouched; the kgcloseB3 retraction stands; nothing here bears on the
paper's correctness.

Disclosed defects of this campaign (rule 5): an early V implementation
dropped the additive odd term from the first hinge, producing an invalid
"certified" cell PASS (+7.69e-5); caught by the frozen reconciliation gate
before any verdict was recorded, retracted, and fixed — the frozen V keeps
A's first hinge verbatim and zeroes only the certified second hinge. The
probe that exposed the mechanism's degeneracy (`s*` ≈ 45.7) was re-derived
machine-side before the Phase-1A run. No other defects known.

Rule-7 scope: this campaign swept exactly tile 3 of [4.083, 6.0] (132-root
cover, full even disk, all compatible odd parts) under V at A's caps, plus
the single-cell decomposition battery (caps 2^17/2^19/2^21 on that one cell)
; tiles 1-2, 4-49, and every other c-domain were NOT swept; the subdivision
mechanism was not deployed on any fresh band; no coefficient-space envelope
change was taken to verdict.

**Named next action (Campaign C, separate pre-statement):** joint (even,
odd) disk-feasibility envelope for the ring cells — re-derive r_o conditioned
on the even part being disk-feasible — plus a panel-cap raise to 2^19 (the
frontier cell alone at 2^19 costs ~68 s CPU, COMPUTATIONAL-EVIDENCE); then
re-run tile 3, and only on a PASS resume tiles 4-49 under the validated
mechanism with A's stop rule.
