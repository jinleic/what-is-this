# Pre-statement — Campaign A: corrected-envelope recertification of band [4.083, 6.0]

Agent `KgEvenBoxRerun`, 2026-08-31. **Status: COMPUTE PENDING — this file is the
first file of the campaign and is committed before any numerical evaluation,
probe, or band run of this campaign.** The predecessor campaign
`campaigns/20260830T084500Z_kgcloseB3/` (git d3ebb62, dad0eaf, f35886c,
f4e3b33) had its two band PASSes retracted to UNVERIFIED; the retraction
stands recorded (rule 5) and is not edited. History is not rewritten; a
clearing can only come from a new machine-verified campaign — this one.

## 0. The two predecessor defects being corrected (owner-verified; restated as facts)

1. **|e|-cap defect (unsound in the permissive direction).** The predecessor
   `j_integrand_bound` capped |e| ≤ 1, but at the paper's own D.4 reproducing
   kernel `p0(s) = (3−s²)/√6` (quoted below, line 2036), |e(0)| = √(3/2) ≈
   1.2247 > 1. A cap excluding the paper's own kernel point can only
   understate the integrand. Owner-verified in sympy at
   (a0, a2) = (√(2/3), −1/√3): |e(0)| = √(2/3) + (1/√3)/√2 = √(3/2).
2. **3-D circumradius inflation where the paper specifies a 2-D even-box
   envelope** (paper lines 1972–1981, quoted below).

## 1. Primary-source quotes (read first-hand this session,
`cs/kg/scratch/paper_full.txt`; verbatim, lines 1972–1981):

> The certificate applies Lemma D.4 with the Cauchy–Schwarz (Christoffel)
> envelopes
> |o(s)| ≤ r_o √K_o(s),
> K_o = ψ₁² + ψ₃² = (5/2)s² − s⁴ + (1/6)s⁶,
> |e(s)| ≤ |e_ā̄(s)| + ρ √K_e(s),
> K_e = ψ₀² + ψ₂² = 3/2 − s² + (1/2)s⁴,
> where ā = (ā₀, ā₂) is the center of an even-coefficient box of circumradius
> ρ, and r_o = √(1−l²) with l = max(0, ‖ā‖ − ρ) bounds the odd-coefficient
> norm of any unit cubic whose even part lies in the box. In this way a
> single evaluation certifies the box against all compatible odd parts, and
> the search space is the two-dimensional even disk rather than S³ (by
> J(c,−p) = J(c,p) one may take a₀ ≥ 0).

and lines 2030–2040 (kernel point, tail discussion):

> +2.539×10⁻⁵
> Table 1: Certificate (C2): the eight sub-bands of the dual-envelope
> branch-and-bound (envelope_band.py, precision 256 bits). Each row certifies
> J(c,p) ≤ d(c) for all unit cubics p and all c in its range; "sub-bands" is
> the geometric c-tiling within the row.
> D.4 Certificate (C3): the high-c tail c ≥ 12
> For c ≥ 12 substitute λ = 1/c ∈ (0, 1/12] and x = cs. Both sides of (56)
> then carry a factor λ: after dividing it out, the target becomes
> D(λ) ≥ I(λ,p) with D(λ) = 3ν/2 · 1/(1+√(1+λ²)),
> I(λ,p) = ∫₀^∞ [(|p(λx)|−x)₊ + (|p(−λx)|−x)₊] φ(λx) dx.
> As λ↓0 this degenerates, and it is worth seeing why. For large c the active
> condition |p(s)| ≥ c|s| selects a small neighbourhood of the origin, where
> J(c,p) = φ(0)p(0)²/c + O(c⁻³). Among unit cubics the largest possible value
> of p(0)² is the evaluation-kernel norm K≤3(0,0) = Σ_{j≤3} ψ_j(0)² = 3/2,
> attained by the normalized reproducing kernel at the origin,
> p0(s) = K≤3(s,0)/√K≤3(0,0) = (3−s²)/√6 = √(2/3) ψ₀ − (1/√3) ψ₂.
> Since φ(0) = ν/2, this gives sup_p J(c,p) ∼ 3ν/(4c), while also
> d(c) ∼ 3ν/(4c): the dual bound is asymptotically tight, with equality
> direction p0.

Supporting text already read first-hand in the same file: fold identity (57)
and f(x,y) = (x+y−t)₊ + (|x−y|−t)₊, t = cs (lines 1963–1966); Lemma D.4
monotonicity (lines 1967–1971); d(c) = (3ν/2)(√(1+c²) − c) (line 1923);
basis ψ₀=1, ψ₁=s, ψ₂=(s²−1)/√2, ψ₃=(s³−3s)/√6 (line 942).

## 2. The corrected 2-D even-box envelope, re-derived from the paper text
(re-used from NO predecessor code; every predecessor instrument was read, but
the formulas below are derived from the quotes in §1)

Notation: for a closed axis-aligned rectangle B ⊂ R²,
B = [α₀,α₀']×[α₂,α₂'] (an "even box"), s ∈ [0, S_MAX], panel P = [a, b],
c ∈ [c_L, c_U]:

- (E1) **Even bound (corner-exact 2-D).** e(s) = a₀ψ₀(s) + a₂ψ₂(s) is affine
  in (a₀, a₂, ψ₂(s)); an affine function on a box attains extrema at
  vertices. Let E_P = [E_lo, E_hi] be the interval enclosure of
  {a₀ψ₀(s) + a₂ψ₂(s) : (a₀,a₂) ∈ B, s ∈ P}, computed in outward-rounded Arb
  from the 8 vertices (a₀/α₀' × a₂/α₂' × ψ₂(a/b)). Then
  **|e(s)| ≤ E⁺_P := max(E_hi, −E_lo) ≤ E_P.upper-hull-outward**, and this
  holds for EVERY unit cubic whose even part lies in B. There is no ρ
  inflation term: the paper's |e(s)| ≤ |e_ā̄(s)| + ρ√K_e(s) is its boxed
  version — with ρ the circumradius of B, |e_ā̄(s)| + ρ√K_e(s) ≥
  max-corner-bound — and using the corner-exact 2-D form directly is strictly
  tighter and equivalent in soundness. NO cap at 1 is applied anywhere
  (counterfactual proof: the cap would clip the paper's own kernel at s = 0).
- (E2) **Odd bound (paper line 1981 verbatim form):** r_o = √(1−l²) with
  l = max(0, dist(0, B)₂) — the distance from the origin to the box, NOT a
  3-D circumradius. Then for every unit cubic with even part in B,
  ‖o‖₂² = 1 − ‖e‖₂² ≤ 1 − l² so |o(s)| ≤ √(1−l²) √K_o(b) on P (K_o monotone
  increasing on [0,∞): exact algebra K_o(s₂)−K_o(s₁) = (s₂−s₁)(s₁+s₂)·q/6
  with q = s₁⁴+s₁²s₂²+s₂⁴−6s₁²−6s₂²+15, and q = (x²+y²)² − x²y² − 6(x²+y²)+15
  ≥ (3/4)r² − 6r + 15 ≥ 3 > 0 for r = x²+y² ≥ 0 — machine-checked below).
  W_P := √(1−l²) √K_o(b) is the paper's own odd-radius formula applied to
  the 2-D box; no 3-D inflation is applied.
- (E3) **Integrand bound (Lemma D.4 direct).** With t = c_L s on panel
  P = [a, b] and threshold lower bound t_L = c_L a ≤ t:
  f(|e(s)|, |o(s)|; t) ≤ f(E⁺_P, W_P; t_L) =: G_P
  where f is the paper's folded integrand; splitting hinges:
  G_P = (E⁺_P + W_P − t_L)₊ + (|E⁺_P − W_P| − t_L)₊,
  each computed as an outward Arb ball from ±, nonnegative part clipped at
  zero (a hinge (x−t_L)₊ is monotone nondecreasing in x, so replacing x by
  its outward upper is sound by Lemma D.4 monotonicity).
- (E4) **Panel weights and total.** U(B, c_L) = Σ_P [Φ(b) − Φ(a)]·G_P +
  Tail_8, with Φ the Gaussian CDF evaluated in outward Arb, exact fmpq panel
  endpoints, and Tail_8 = closed-form tail majorant (paper lines 1982–1984)
  computed in outward Arb (value ≈ 3.635317254e-13 < 3.7e-13 as the paper
  REPORTED).
- (E5) **Band pairing** (paper line 1984): J(c,p) ≤ J(c_L,p) and
  d(c) ≥ d(c_U); certify d(c_U).lower() − U(B, c_L).upper() > 0 in Arb
  (strict endpoint separation). A full band passes iff every disk-
  intersecting leaf of every c-tile passes.

**Soundness claim (machine-checked by startup asserts below):** every
substitution (E1)–(E4) is nondecreasing under Lemma D.4 or an outward
interval operation, so U(B, c_L) ≥ J(c_L, p) for every unit cubic p with
even part in B, hence U(B, c_L) ≥ J(c, p) for every c ∈ [c_L, c_U] by the
band pairing.

## 3. Startup asserts (counterfactual discipline) — MUST pass before any band run

At Arb 320 bits (displayed anchors at 300-bit accuracy), abort on any failure:

1. **Kernel-admission assert (the retraction's core).** At the kernel point
   box `(a0,a2) = (√(2/3), −1/√3)` (a point box):
   - `a0² + a2² = 1` (ball radius < 1e-90);
   - `e(0) = a0 − a2/√2` has outward upper ≥ √(3/2) − 2⁻²⁶⁰ and lower ≤
     √(3/2) + 2⁻²⁶⁰, where √(3/2) = 1.2247448713915890490986420373529…;
   - the envelope's E_P at the point box/panel-union [0,0] satisfies
     E_P ≥ |e(0)| with slack < 2⁻²⁶⁰ (no cap, no 3-D inflation);
   - **counterfactual (must be REJECTED):** a cap `min(E_P, 1)` would give
     E_P^cap = 1 < √(3/2) — deficit 0.2247448713915890490986420373… is
     computed and asserted > 0 to show the cap is detectably unsound.
2. **Odd-radius assert.** On a demo box B = [0.8, 0.9]×[−0.6, −0.5]:
   l = dist(0,B) = √0.89 outward; r_o = √(1−0.89) computed; and the 3-D
   circumradius-counterfactual r_3D = circumradius of the box in R³ (corner
   to center+ext) is asserted STRICTLY GREATER than r_o, showing the 3-D form
   is permissively different (looser).
3. **K_o/K_e identities:** ψ₀²+ψ₂² = K_e and ψ₁²+ψ₃² = K_o exactly at
   s ∈ {0, 1, 1.7, 2, 3} (bit-equal Arb enclosures), matching paper lines
   1976/1980.
4. **Cited-tail anchor (independent of this campaign's envelope).** With N=251
   and ‖D³H‖ ≤ 14.44243664663976457 imported CITED from the authors' frozen
   `d3h_certificate.json` (sha256
   80e945589b796408f5de5add837272037160d163467050b354f15ab9d6822d2e, verified
   byte-identical this session), tail = B3/√(10·N⁵) must reproduce
   4.5756855097934878…e-6, agreeing with the frozen `4.57568550979349e-6`
   display string to its stated precision, and = 4.57569e-6 to 6 significant
   digits. This tail is consumed as a cited constant of Gate A's retraction
   record, is independent of the Gate-B envelope, and anchors provenance of
   the instrument chain (its malfunctioning is what would falsify reuse).
5. **Margin semantics counterfactuals.** On planted balls: (a) d.lower() −
   U.upper() > 0 is required for PASS; (b) a straddling ball must FAIL; (c)
   an upper bound strictly above d must FAIL. (Same battery as predecessor.)
6. **Envelope-admissibility assert:** the analytic even anchor
   `p* = (4/5)ψ₀ − (3/5)ψ₂` (unit norm) must satisfy U(B*, c_L) ≥ J(c_L, p*)
   at the point/box c=1.30 with a positive slack, and the closed-form
   J(1.30, p*) must be inside the independently frozen restart-campaign
   window [0.37483055458876357002453330452976048, …050] (provenance anchor).
   [This controls that the new envelope is never looser than the exact
  一点的 computation on an even-only cubic.]

## 4. Frozen band, envelope, and run parameters (Rule 7 + Rule 16: fixed BEFORE the run)

- **Band:** `c ∈ [4.083, 6.0]` (the retracted band; numbers `4.083` and `6.0`
  as exact decimals, tiled by adjacent ratio-1.02 pairs [c_L, min(1.02·c_L,
  6.0)], final tile clamped to exactly 6.0 — identical tiling to predecessor
  campaigns, so margins are comparable).
- **Even-disk cover:** exact 12×12 root partition of [−1,1]², 132 boxes
  intersecting the closed unit disk (pre-anchored count), binary longer-axis
  splits, maximum depth 40, minimum side 1/1024 (matching the frozen restart
  instrument). No below-floor refinement.
- **Panel ladder:** start 4096 for 4 < c_L ≤ 8, times 4, unique, capped once
  at 32768. (4096 for this band per the frozen instrument's band table).
- **Precision:** 256 Arb bits for the envelope, 320-bit analytic startup
  asserts. Python-flint (Arb) outward rounding; no floats in trusted paths.
- **Budget:** 260,000 distinct (box, panel-count) envelope evaluations per
  named band; budget exhaustion = FAILURE TO CERTIFY + exact frontier (same
  rule as all predecessor campaigns).
- **S_MAX = 8**, tail closed-form (frozen core semantics), no below-8 panels.
- **Verdict rule:** tile PASS iff every disk-intersecting leaf in the tile
  has d(c_U).lower() − U(B, c_L).upper() > 0 in Arb. Band PASS iff every
  tile passes. Any leaf with margin ≤ 0 at floor depth and panel cap means
  the band is OPEN at that leaf → band verdict is FAILURE TO CERTIFY naming
  the subinterval and the exact margin ball.
- **Refutation gate (unchanged):** a certified statement J.lower() > d.upper()
  on any leaf of any tile = candidate refutation of the paper's Lemma D.2 —
  STOP, escalate to Main via hub with the exact cell and intervals; nothing
  else written. (Predecessor campaigns never fired this; neither may we
  casually.)

## 5. Anchored tail scope statement (rule 7 sentence, restated)

This campaign re-certifies exactly the named band `c ∈ [4.083, 6.0]` over the
full even-coefficient disk with all compatible odd parts via the paper's 2-D
even-box/odd-radius envelope at the review's specified parameters, using
interval panels on [0,8] plus the closed tail; it does not re-derive the
paper's C1 (Lemma D.1), the splice (S), or the C3 high-c λ-chart; it does not
touch c < 4.083 or c > 6.0; the only consumed external constant is the
authors' frozen ‖D³H‖ ≤ 14.44243664663976457 (CITED, sha256-verified, used
only in the cited-tail provenance assert), which has no logical role in this
campaign's envelope verdicts.

## 6. Evidence labels

- MACHINE-VERIFIED: any margin/sign claim from outward Arb endpoints on the
  stated coefficient/c/panel/tail domain, scripts frozen in `code/`.
- COMPUTATIONAL-EVIDENCE: float/mpmath values used only as orientation.
- REPORTED: paper Table-1 numbers, including the +2.539e-5/−3.571e-5 margins
  quoted from Table 1 (not re-run).
- FAILURE TO CERTIFY is a result and is reported honestly; PASS requires the
  positive strict endpoint margin on EVERY leaf of EVERY tile of the band.

## ADDENDUM 1 (recorded BEFORE any full-band run; per rule 16 a new pre-registered
statement reported alongside the original §4 tiling, not a silent change)

**The §4 tiling ratio 1.02 is arithmetic-impossible for this band — changed to 1.008
before the band run.**

Orientation achieved with the CORRECTED envelope before any band sweep (thin point-
box at the paper's kernel direction, panels 4096→524288):

- U(B, 4.083) → 0.142747244 at n=524288 vs independent closed-form
  J(4.083, p0) = 0.14273947723612726… (mpmath 50-dps and exact-Arb closed form
  agree): the corrected envelope is asymptotically EXACT at the kernel point.
- The wall is therefore the tile pair, not the envelope: with c_U = 1.02·c_L the
  kernel slab needs J(4.083,p0) ≤ d(4.16466), i.e. slack −1.0645e-3 — uncertifiable
  at ANY panel count with ANY correct envelope. At ratio 1.008 the same kernel slab
  certifies with margin +5.45e-4 (n=131072, certified); a moderate panel ladder with
  cap ≥ 131072 suffices there.

Parameters frozen for the band run (supersede the single conflicting §4 sentence):

1. c-tiling: adjacent geometric ratio **1.008** tiles [c_L, min(1.008·c_L, 6.0)],
   final tile clamped exactly to 6.0 — still adjacent, still gapless, still
   deterministic; tile-adjacency assert unchanged.
2. Panel ladder: start 4096, ×4, unique, **capped once at 524288** (2^19).
3. Budget: 260,000 distinct envelope evaluations per named band (unchanged).
4. Everything else in §4 stands unchanged (root cover 132, floor 1/1024, depth 40,
   tail, S_MAX, precision, verdict rule, refutation gate).

Honest disclosure: ratio 1.02 was listed in §4 as "identical tiling to predecessor"
for margin comparability. The pre-run diagnosis above (made with the corrected
envelope only, before any band-level compute) showed that choice cannot certify
tile 1 of the named band, so per rule 16 the new tiling is pre-registered here
with its rationale; the §4 text stands as originally committed.

## ADDENDUM 2 (recorded BEFORE any verdict-bearing tile compute; the Addendum-1
run produced NO tile verdict in ~50 minutes of wall time and was stopped; its
discarded partial compute is disclosed here and cited nowhere)

**Panel cap recomputation: 524288 → 131072, with pre-registered cost/tightness
accounting.** The Addendum-1 band launch ([4.083,6.0], ratio 1.008, cap 2^19)
spent >50 min inside tile 1/49 without a verdict: a single full-ladder envelope
evaluates in ~15.5 s at cap 2^19 (measured out-of-band: 9.33 s at n=524288 alone),
and kernel-adjacent root boxes split to the 1/1024 floor, so tile 1 alone
projects to hours and the band to days. Recomputed cap 131072 (2^17):

1. Tightness loss, measured on the worst cell (kernel thin-box at c_L=4.083,
   c_U=1.008·c_L): U(524288)=0.142747244 vs U(131072)=0.142769295 — the margin
   falls from +5.44843e-4 to +5.22e-4 (still comfortably strictly positive).
2. Single-ladder cost at cap 131072 ≈ 0.07+0.29+1.16+4.68 ≈ 6.2 s (vs 15.5 s);
   the projected band drops to a few hours.
3. Budget, tiling, floor, depth, precision, verdict rule, refutation gate:
   unchanged from Addendum 1.

The Addendum-1 statement stands as originally committed; this addendum
supersedes exactly its panel-cap sentence.

## ADDENDUM 3 (environment-only, recorded before the third launch of the band;
no mathematical parameter changes)

The Addendum-2 launch showed ~2% CPU duty (6 s CPU in 4.5 min wall) despite an
otherwise idle machine — macOS App-Nap throttling of background `Python.app`
processes, cross-confirmed with a sibling agent (their fresh censuses identical
symptom). Fix: relaunch under `caffeinate -i` (prevent idle nap) WITHOUT nice,
same command line, same code bytes, same panel cap. The two earlier partial
launches produced ZERO tile verdicts; both are disclosed as discarded smoke and
cited nowhere. Addendum-2 parameters stand.

## ADDENDUM 4 (environment-only; before the fourth launch)

Caffeinate did not change the ~2% CPU duty (measured: 5 s CPU / 4 min wall on the
Addendum-3 launch). Fourth launch wrapper: `taskpolicy -c utility` (macOS QoS
class utility, no nice). Same code bytes, same parameters; Addendum-2/3 params
stand. No verdicts existed in the three earlier partial launches.

**Fourth-launch status (recorded at +52 min):** alive, PID 84917, but duty still
~2.4% CPU (1 s CPU / 52 s wall) — taskpolicy -c utility did not lift the QoS
throttle either. Expect multi-hour wall per tile at this duty; verdicts pending.
This state note is a Monitoring record, not a verdict of any kind.

## ADDENDUM 5 (environment-only monitoring record, before the fifth launch)

Foreground probe (same instrument bytes, no wrapper): a full-ladder envelope at
cap 131072 runs in **2.34 s** — the QoS throttle applies to *launched detached
processes*, not to the python interpreter per se. Fifth launch: plain command,
no nice/caffeinate/taskpolicy, thread envs pinned to 1 (as in every launch).
The instrument, parameters, and verdict rules are byte-identical to Addendum 2;
this addendum records the launcher only. Addendum 4's launched process was
stopped (still ~2.4% duty) before this relaunch; it produced no verdicts.

## CAMPAIGN A RESULT (frozen in logs/band_4p083_6_result.json, commit 0f494f3)

**Band [4.083, 6.0]: FAILURE TO CERTIFY** at tile 3/49 under the corrected
envelope, per the pre-registered stop rule (§4: an open leaf at panel cap and
floor depth earns FAILURE TO CERTIFY with the exact frontier). A FAIL here is a
result; the kgcloseB3 retraction stands.

Per-tile verdicts (all margins Arb 256-bit outward):

| tile | c-pair | verdict | worst certified margin | evals |
|---|---|---|---|---|
| 1 | (4.083, 4.1156640) | PASS | +1.04499016763434518123845756638e-6 ± 4.27e-36 | 3,354 |
| 2 | (4.1156640, 4.1485893) | PASS | +4.71904853878616890385690164962e-6 ± 3.07e-36 | 3,476 |
| 3 | (4.1485893, 4.1817780) | **OPEN** | open leaf margin −1.76291787028287693244667476546e-5 ± 4.20e-35 | 819 |
| 4–49 | — | NOT RUN (pre-registered stop at first failing tile) | — | — |

Tile-3 frontier (exact frozen coordinates): even box
a0 ∈ [−0.830729166666666666666667, −0.829427083333333333333333],
a2 ∈ [+0.558593750000000000000000, +0.559895833333333333333333],
depth 14, panels 131072 (cap), 95/251 leaves certified, 124 boxes stacked when
stopped. The frontier cell sits on the disk ring (‖a‖ ≈ 1), where
r_o = √(1−dist²) → 0: the |e|-hinge panel-sup is the binding slack, and the
cell survives as OPEN — a failure of THIS instrument at the named caps, not a
lower bound on J−d and not a paper refutation (rule-14 semantics unchanged).

Runner disclosure: all detached launches (six) were throttled ≈40× by an
illustrative ps-accounting artifact at ~2.4% —/Main's later instrument ruling:
ps deltas are untrustworthy for these Python.app children; the in-process
getrusage/sample instrument is authoritative. The verdict-bearing run was
re-executed in-session (unthrottled foreground kernel), tile by tile, with
frozen per-tile artifacts `inproc_tile_*.json` and identical instrument bytes
(frozen code/gate_corr_env.py, sha-verified below). Per-tile wall times:
2,665 s (tile 1), 3,455 s (tile 2), 863 s (tile 3).
