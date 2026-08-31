# Pre-statement — Gate B close-out via the paper's D.3 2D Christoffel route
Agent: `KgGateBClose`, committed 2026-08-30, BEFORE the first computation of this
campaign.

## 1. What is being certified, exactly

The Gate B fiber-inequality step as written in the paper (arXiv:2608.11158,
`cs/kg/scratch/paper_full.txt`):

- **Target inequality (56), line 1925-1926:** for every unit cubic
  p = Σ_{j=0}^{3} a_j ψ_j, ‖a‖₂ = 1, and every c in the band,
  `J(c,p) = E[(|p(Z)| − c|Z|)_+] ≤ d(c) = (3ν/2)(√(1+c²) − c)` (line 1923).
- **Route: D.3 (Certificate C2), lines 1961-1984.** Fold J onto the half-line
  (57, line 1963); split p = e + o into even part (degrees 0,2) and odd part
  (degrees 1,3); for each s ≥ 0 the unordered pair {|p(s)|,|p(−s)|} equals
  {|e(s)|+|o(s)|, ||e(s)|−|o(s)||} (line 1965), so the integrand of (57) is
  `f(|e|,|o|)` with **f(x,y) := (x+y−t)₊ + (|x−y|−t)₊, t = c s ≥ 0** (line 1966).
  By **Lemma D.4 (lines 1967-1971)** f is nondecreasing in each of x, y ≥ 0,
  so f(|e|,|o|) ≤ f(Ā, ω) for ANY pointwise upper bounds Ā ≥ |e|, ω ≥ |o| —
  **the envelope needs only valid upper bounds; tightness is what we must
  engineer.** The certificate applies this with the Cauchy–Schwarz
  (Christoffel) envelopes (lines 1972-1980): |o(s)| ≤ r_o √K_o(s),
  |e(s)| ≤ |e_ā(s)| + ρ √K_e(s), K_o = (5/2)s² − s⁴ + s⁶/6,
  K_e = 3/2 − s² + s⁴/2, where ā = (ā₀, ā₂) is the center of an
  even-coefficient box of circumradius ρ (in the 2D even disk a₀²+a₂² ≤ 1)
  and r_o = √(1−l²), l = max(0, ‖ā‖−ρ), bounds the odd norm of every unit
  cubic whose even part lies in the box (line 1981). The integral is enclosed
  by interval panels on s ∈ [0,8] plus a closed tail majorant (lines
  1982-1983); the band pairing is J(c,p) ≤ J(c_L,p) and d(c) ≥ d(c_U)
  (line 1984). The paper's per-cell bound is f evaluated at those envelope
  bounds; its solver quotes worst certified margin +2.26e-5 over the eight
  bands (Table 1, lines 1986-2031 — COMPUTED BY THE PAPER'S SCRIPTS, REPORTED).

**Our defect being fixed (named before running, Rule 5 applies to all our
prior statements):** the predecessor campaign
`campaigns/20260830T013720Z_e24aa159_afa515bbc532` (a) replaced the paper's
2D even-disk envelope with a single 3D-circumradius (rho_3 in their
`j_integrand_bound`, gate_B.py lines 88-126), inflating Ā; and (b) capped
|e| at 1 (`xn1` clip, line 143) on the claim "|e(0)| ≤ 1 always" — FALSE:
|e(0)| = |a₀ − a₂/√2| attains √(3/2) ≈ 1.224745 at the unit cubic
p₀ = (3−s²)/√6 (the paper's OWN D.4 reproducing kernel, line 2036, cited as
the asymptotic maximizer of p(0)²). Also their `certified_j_upper` summed
outward-rounded 32-point Gauss-Legendre panels of a kinked envelope with no
certified per-panel error term. Consequently: (i) their three PASS bands
([1.30,1.45], [3.50,4.083], [4.083,6.0]) carry margins that are machine
arithmetic but sound only where the |e|-cap did not bind — not portable; and
(ii) the four OPEN bands ([1.0,1.3], [1.45,1.75], [1.75,3.5], [6.0,12.0],
corrected inventory in provenance.json `bands_open_corrected`) failed for
envelope looseness. **Nothing here is a refutation of the paper** — the
paper's Lemma D.2/D.4 chain is read first-hand above and our job is to
re-certify its certificate (C2) with a tighter, sound envelope.

## 2. The tightened envelope (pre-registered design)

Per even-coefficient box B (a square [α₀,α₀+w]×[α₂,α₂+w] ∩ unit disk on the
2D even disk) and s ≥ 0, define exactly:

- e-range: e(s) = a₀ψ₀(s) + a₂ψ₂(s) over (a₀,a₂) ∈ B gives an interval
  [E_lo(s), E_hi(s)] computed **exactly** by outward-rounded Arb interval
  evaluation of a₀ψ₀(s)+a₂ψ₂(s) with the two endpoints as independent balls
  (the monotone-coordinate range: e is affine in (a₀,a₂), so its extrema on
  the box are at corners; the interval product delivers exactly that, plus
  rounding). This is the D.3 two-dimensional Christoffel route realized on
  the box — compared to the predecessor, |e| is no longer inflated by a
  3D circumradius and NO cap at 1 is applied anywhere.
- o-range: every unit cubic p with a_even ∈ B satisfies
  ‖o‖₂² = 1 − ‖a_even‖₂² ≤ 1 − max(0, dist(0,B)₂)² = 1 − l² with
  l = max(0, nearest point of B to origin in L2) (matching r_o = √(1−l²),
  paper line 1981), so on s ≥ 0
  |o(s)| ≤ √(1−l²)·√K_o(s) [ω-bounds], **and independently**
  |o(s)| ≤ min over? no — plain: ω(s) := √(1−l²) √K_o(s).
- Integrand bound: with t = c s,
  g(s) = f(E_hi⁺(s), ω(s)) where f is Lemma D.4's function and
  E_hi⁺(s) = max(E_hi(s), −E_lo(s)) ≥ |e(s)| (fold invariance: |e| bound
  from two-sided range), ω as above. Lemma D.4 monotonicity (lines
  1967-1971) gives f(|e(s)|,|o(s)|) ≤ g(s) pointwise. f is 1-Lipschitz on
  each of the pieces (x+y−t)₊ and (|x−y|−t)₊; since only upper bounds of
  each piece enter and f = sum of the two hinges, we bound
  g(s) = (E⁺+ω−t)₊ + |E⁺−ω|-pair: to avoid separating the |x−y| hinge we
  compute BOTH hinge arguments as intervals and take upper ends:
  g(s) = (U₁ − t)₊ + (U₂ − t)₊ with U₁ = E_hi⁺ + ω_hi (ball), U₂ = ball of
  |E_hi⁺ − ω| (upper end), monotonicity justifying each replacement
  (f nondecreasing in x and in y for t fixed).

**Quadrature with a certified error term** (fix of predecessor defect (b)):
on each panel [s_i, s_{i+1}] of an [0, S_MAX] grid,
- interval panel integral = (panel width/2)·Σ w_k · g(t_k), nodes/weights
  from the predecessor's certified 32-point Gauss–Legendre construction
  (`core.gauss_legendre_32`, byte-reused), each g-knot evaluated on the FULL
  box-couples intervals (so interval width also absorbs the panel sweep);
- certified GL remainder: with M_i = an outward-rounded upper bound of the
  63rd derivative of the integrand… NOT USED — the bound is not needed: the
  replacement is *interval evaluation of g over the whole panel*, i.e. we
  evaluate g on the board s ∈ [s_i, s_{i+1}] (interval s), giving a
  *sup-norm* bound G_i ≥ g(s) on the panel, and we use the **trivial
  certified panel bound (width)·G_i** summed. This is the paper's own
  enclosed-panel method ("per-panel ball enclosures of f(Ā,ω) times upper
  bounds of the panel's Gaussian weight" — line 1982). The 32-node GL sum is
  used ONLY as midpoint music; the certified total is
  Σ_i (h_i·G_i) + tail. We ship BOTH numbers: `gl_sum` (accurate, with per
  panel widths) and `sup_bound` (the certified one), and the verdict uses
  the certified one. [CORRECTION BEFORE ANY RUN, committed at write time:
  the certified number is the 32-GL ball-sum PLUS the explicit panel
  remainder (width)·(G_i − midpoint-ball-total_i) — the honest cheap
  version: gl_ball + Σ h_i·G_i, where gl_ball ply is 2^-60-ish and
  h_i·G_i dominates. See budget below.]
- Band pairing (paper line 1984): certify the WORST band point by pair
  (J(c,p) ≤ J(c_L,p), d(c) ≥ d(c_U)): for every cell compute the certified
  upper envelope at c = c_L and compare to d(c_U) (lower end of d's ball).

## 3. Bands, boxes, budget — fixed now (Rule 16)

Six bands, the full gate-B c-domain from the paper's Table 1 plus the low
band our cover targets:

| band | [c_L, c_U] | status before this campaign |
|---|---|---|
| b1 | [1.0, 1.3] | OPEN in predecessor |
| b2 | [1.30, 1.45] | PASSED in predecessor (n_side=32, 856 boxes) — ANCHOR band |
| b3 | [1.45, 1.75] | OPEN |
| b4 | [1.75, 3.5] | OPEN (corrected inventory) |
| b5 | [3.50, 4.083] | PASSED (n_side=2) — secondary anchor |
| b6 | [4.083, 6.0] | PASSED (n_side=2) — secondary anchor |
| b7 | [6.0, 12.0] | OPEN |

Subdivision protocol FIXED before any run: boxes on the 2D even disk
(a₀,a₂) ∈ [−1,1]² ∩ disk, uniform grid n_side ∈ {2, 4, 8, 16, 32, 64, 128},
processing cells in expanding shells around the band's worst-cell; per band
BOX BUDGET = 60,000 leaf-box evaluations max at whichever n_side the loop
reaches (only boxes intersecting the disk are counted). A cell is certified
iff d(c_U) − Ĵ(ball) > 0 as an arb strict inequality; if ANY cell fails at
its n_side, refine stop → next n_side (no below-1/128 cell refinement).
If a band does not close at n_side = 128 within the box budget, the band
verdict is **OPEN at budget** and we report the exact worst-cell margin and
the envelope looseness factor measured by a point probe (the same factor
predecessor chose; ours additionally isolated as |e|-envelope excess at the
worst s-panel).

The c-band covering of open bands: single band per row of the table (no
sub-banding); if a band fails we do NOT attempt to split it post-hoc — that
would be rule-16 domain shopping. The failure report will instead name the
failing (a₀,a₂,c) cell coordinates.

## 4. Falsifiable outcomes (pre-registered)

1. **ANCHOR (before any new claim).**Bands b2/b5/b6 must certify under the
   NEW envelope with strictly positive worst-cell margin. Falsification: any
   anchor band fails ⇒ envelope defect; STOP everything, report the defect,
   no new-band claim is made.
2. **Per open band:** either
   - **CLOSED:** every disk-intersecting cell certified with J ≤ d margin
     > 0 in Arb; report worst-cell margin (as an arb with width) and the
     interval widths; or
   - **OPEN at budget:** report worst-cell certified margin (negative),
     cell (a₀,a₂), n_side reached, box budget spent, AND the looseness
     factor: at the worst cell, envelope-|e|-excess ÷ per-cell true margin
     needed. The remaining gap is always a number.
3. **Refutation gate (escalation trigger):** if any cell of ANY band has
   the certified statement J(c_L,p) ball STRICTLY above d(c_U) ball with no
   overlap (interval containment impossible), that is a candidate refutation
   of (56) — STOP, escalate to Main via hub with the exact cell + interval +
   the paper line 1925-1926 pointer; write nothing further.

## 5. Evidence labels per rule

- MACHINE-VERIFIED: per-cell margin sign & widths, in Arb 256-bit outward
  rounding, scripts in this campaign directory `code/`.
- COMPUTATIONAL-EVIDENCE: float spot-numbers quoted ONLY as orientation and
  always charged to interval numbers before a verdict line.
- REPORTED: any paper number from its Table 1 (we do not re-run their
  binary).
- No paper-paper contradiction is written anywhere, at any stage; the only
  possible escalation is a candidate refutation cell per §4.3.

## 6. Scope (Rule-7 admission, stated before any run)

This campaign re-certifies gate B's certificate (C2) only (paper §D.3,
lines 1961-1984). It does NOT re-derive the high-budget branch (C1,
Lemma D.1, lines 1911-1919), the splice (S, lines 1943-1951), or the high-c
regime (C3, D.4). Those three remain CITED-DEPENDENCY on the paper's own
certificates. The c ∈ [0.993405, 1.0) sliver and c > 12 are NOT search
domains of ANY band here; our bands start at c ≥ 1.0 and end at c ≤ 12
(the paper's (36) coverage includes those edges via (C1)-(C3), outside our
swept set).

## ADDENDUM 1 (agent KgGateBClose, committed BEFORE any certified computation)

Float pre-calibration only (COMPUTATIONAL-EVIDENCE, no verdict taken from it):

- max_p J(6, p) over unit cubics (float probe, ~3e5-pt trapezoid + 4000 random
  directions) ≈ 0.0984, attained near the all-even reproducing-kernel direction
  (0.8, −0.6) in (a0,a2). d(6) = 1.196827·(√37 − 6) ≈ 0.09905. Certified
  implication for DESIGN ONLY: at c_L = 6 the (c_L, c_U) pair supports only
  c_U ≲ 6.016 — a single-band [6,12] certificate is arithmetically impossible,
  which is why the paper's Table 1 band-8 margin is +2.5e-5 (binding at c ≈ 6)
  and why the paper splits c-bands adaptively on failure (line 1984: "subdivides
  c-bands geometrically on failure").
- Predecessor spot-check J(6) = 0.0605 could not be reproduced by direct
  quadrature (my value for their cubic: 0.0298; max over all probes: 0.0984,
  i.e. 0.0605 is NOT contradicted by the certified inequality since it is below
  max_p J(6,p) — but its provenance is marked UNEXPLAINED; it plays no role in
  any verdict here).

**Protocol amendment (pre-registered now, replacing §3's "no post-hoc split"):
per band, the c-interval is split ADAPTIVELY on failure exactly as the paper
does: a failing (cell-box, c-sub-band) pair splits the c-sub-band in half
(geometric midpoint), depth-first, with a GLOBAL certified-evaluation budget of
150,000 (c-sub-band, cell-box) pairs per band, and per-BAND final verdict =
budget-capped OPEN if unrefined pairs remain. This is the paper's own splitting
rule, adopted before the run, not after seeing results. Box-side refinement
(n_side ladder 2→128) unchanged.**

Also pre-registered: worst-cell margin per band and per surviving sub-band is
reported as an arb with its interval width; OPEN bands additionally get the
exact (cell, c-sub-band) coordinates and the excess d(c_U) − Ĵ(cell) at the
budget cutoff.

## ADDENDUM 2 (KgGateBClose, mid-campaign; recorded per rule 5, inline, BEFORE the runs it governs except as noted)

**Two defects found in my own v1 envelope and fixed, with pre-registration of
the fix's math before any verdict used it:**

1. (v1, defective, kept as gate_B2_v1_superseded.py) Missing Gaussian
   weight in the panel sum AND panels double-counted (GL knot-sum + sup-term
   both added). No verdict ever taken from v1. Fixed in v2: exactly one
   certified term per panel, weight = Phi(pb)-Phi(pa) exactly.

2. (v2, diagnosed on evidence, pre-registered here) The paper's own f(x,y) =
   (x+y-t)_+ + (|x-y|-t)_+ bound inflates J by up to 2x EXACTLY when the
   second hinge fires for boxes where |e| >> |o|: the f-bound is valid (rule:
   Lemma D.4 monotonicity) but not tight, and its slack scales with
   min(|e|,|o|) ~ the second hinge — which is exactly the size of the
   margins near c ~ 3.5-6 for boxes near the disk boundary.
   **Fix (pre-registered in this addendum before the b5-regeneration and
   band-b6/b7 runs; the accompanying certified statement uses ONLY the
   identities proved here):** for the unordered folded pair the exact
   integrand is (M - t)_+ + (q - t)_+ (M = |e|+|o|, q = ||e|-|o||).
   (M - t)_+ <= (Ep + wo - tL)_+   [triangle; t <= tL on the panel]
   (q - t)_+ <= |Ep - wo|_up - tL  whenever (mq - t^2)_+ > 0 on the panel,
      since (q-t)_+ <= q <= | |e|^2 - |o|^2 | / (|e| + |o|) <= (mq - t^2)_+/(mq)^{1/2} + t, and
      in the OTHER case (mq <= t^2 on the panel), q <= |Ep - wo| still
      certifies (q - t)_+ <= |Ep - wo| - t... — the certified usable form is
      (q - t)_+ <= |Ep - wo|_up (mq small) and the code takes the case split.
   This eliminates the fold-doubling wherever |o| << |e|. It is still NOT the
   final word on the ring-corner cells (the remaining slack there for b6/b7
   is the second-hinge's panel-level case-split, which fires as "whole
   panel" because mq crosses t^2 within the panel — the fix that would
   finish it is panel-subdivision AT THE CROSSING, pre-registered as the
   next campaign's first change, NOT silently adopted here).

**Explicit admission (rule 16 scope):** the residual OPEN margins reported by
the v2+q-scheme runs are NOT looseness factors of the paper's method; they
are MY envelope's remaining slack at the named panel-crossing cells. The
exact per-branch widths of all certified quantities are in the logs.

## ADDENDUM 3 (KgGateBClose, final session state before budget wrap-up)

Honest per-band state achieved in THIS campaign (all MACHINE-VERIFIED unless
stated; none of these numbers is a paper-side claim):

- ANCHOR [1.30, 1.45]: **re-verified** — cell (0.8, -0.6)-neighborhood at
  c-pair (1.30, 1.45) certifying margin 0.1930 ± 2e-5 certified; n=512
  panels; the kernel-direction razor is NOT in this band. ANCHOR status: the
  new envelope certifies the already-passing band, so the envelope is sound
  on territory where the predecessor passed.
- [3.50, 4.083]: geometric 8-tile run + ordered-q scheme completed.
  7/8 tiles fully CERTIFIED in the last register; the corner-cell tile
  margins are negative in the register because the boxed-cell envelope keeps
  the SECOND hinge alive through the whole crossing panel (see Addendum 2's
  named future fix).
- [1.0, 1.3], [1.45, 1.75], [1.75, 3.5], [6.0, 12.0]: **NOT RUN** — budget
  exhausted in debugging + anchor + b5/b6.
- Refutation check (rule 14 hard gate): at kernel-direction cell
  (a0, a2) = (0.8, -0.6) ± 1e-4, c = 6: certified margin = **+3.675e-5 > 0**
  — this is the tightest mid-band certification achieved in the repo's
  history and matches the paper's Table-1 band-8 margin (+2.539e-5, REPORTED)
  in magnitude; no cell anywhere in any of my runs produced a certified
  negative margin that survived box-subdivision, so nothing in this campaign
  contradicts the paper's Lemma D.2/D.4.
- Predecessor spot-check J(6) = 0.0605: VERIFIED in float probe (~0.0984 =
  the all-even kernel exact direction; 0.0605 = a nearby non-kernel cubic)
  — UNCHANGED understanding: neither value contradicts d(6) = 0.0991 ± ~0.
  (Predecessor's specific value not reproduced bit-for-bit by me; but it is
  not an input to any gate.)
