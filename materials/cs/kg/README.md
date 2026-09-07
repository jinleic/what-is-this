# `kg/` — the Grothendieck constant $K_G$

**Status: BENCHMARK (no campaign run yet).** Independently re-verify, and then
attempt to improve, the 2026 bounds on the Grothendieck constant. Selected
because both directions of the current record reduce to a **finite collection of
one-dimensional Gaussian-integral inequalities discharged by rigorous interval
arithmetic** — the exact machinery this workstation already runs at scale
(`python-flint`/Arb; cf. `../../math/uc/` 488,465,854-node Arb certificate).

## The constant

$K_G$ is the infimum of $K$ in Grothendieck's inequality; equivalently the
integrality gap of the canonical SDP relaxation of $\max_{x,y\in\{\pm1\}^{n}}
x^{\top}Ay$. Its exact value is open.

Bound history, as read first-hand in arXiv:2608.11158 §1 (**V**, owner-read
2026-08-29):

| bound | value | source as cited there |
|---|---|---|
| lower | $\pi/2$ | Grothendieck |
| lower | $K_{\mathrm{DR}}=1.6769\ldots$ | Davie; Reeds (independent; best for four decades) |
| lower | $K_{\mathrm{DR}}+10^{-26}$ | Heilman |
| lower | $K_{\mathrm{DR}}+10^{-12}$ | Jones–Malavolta |
| **lower** | $\mathbf{6\pi/11=1.7135\ldots}$ | **arXiv:2608.11158** |
| upper | $\pi/(2\log(1+\sqrt2))=1.7822\ldots$ | Krivine |
| upper | strict $<$ Krivine, unquantified | Braverman–Makarychev–Makarychev–Naor 2011 |
| upper | $\approx10^{-5}$ explicit improvements | Heilman; Li et al. |
| **upper** | $\mathbf{\pi/(2\log(1+\sqrt2))-3.47\times10^{-4}=1.7818\ldots}$ | **arXiv:2608.11158** |

The abstract states the weaker $-10^{-4}$; §1 item 1 states $-3.47\times10^{-4}$.
Both are recorded; the gates target the stronger stated constant.

**Owner-verified arithmetic (MACHINE-VERIFIED, Arb, `ctx.prec=400`):**

```
6*pi/11                    = [1.713595992867159948252351 +/- 6.37e-26]
pi/(2 log(1+sqrt 2))       = [1.782213978191369111774413 +/- 4.53e-25]
Krivine - 1e-4             = [1.782113978191369111774413 +/- 4.53e-25]
interval width             = [0.06851798532420916352206252 +/- 3.37e-27]
both endpoints in (1.7,1.8): True   ->  tenths digit of K_G is 7
6*pi/11 > 1.67696 (beats K_DR):     True
```

So the paper's headline corollary — the tenths digit of $K_G$ is $7$ — follows
from its two bounds by rigorous arithmetic. **This certifies the corollary, not
the bounds.** The bounds themselves are `CITED-DEPENDENCY` until gates A and B
discharge them.

## Why this is a certificate target, not a survey target

Both directions terminate in interval-checkable finite inequalities
(arXiv:2608.11158 §2, owner-read):

* **Upper bound.** Krivine preprocessing gives $M(\gamma)\le1\Rightarrow
  K_G\le\pi/(2\gamma)$, where $M(s)=\sum_{n\ge1}|a_n|s^n$ is the majorant of
  $H^{-1}(z)=\sum a_nz^n$. The paper proves the *sufficient* condition
  $$\gamma+\Delta_H<b_1,\qquad H(t)=b_1t+\sum_{m\ge3}b_mt^m,\quad
    \Delta_H:=\sum_{m\ge3}|b_m|,\quad \gamma=\tfrac{\pi}{2K_G^{UB}},$$
  so the inverse series is never needed: certify a **lower bound on $b_1$** and
  an **upper bound on $\Delta_H$**. The $b_n$ come from Hermite coefficients of
  $f,g$ (one-dimensional Gaussian integrals); the $\Delta_H$ tail is closed by a
  single boundary-norm estimate on $D^3H$, $D=t\,d/dt$. The paper states
  outright: *"Rigorous interval arithmetic verifies the resulting finite
  collection of inequalities."*
* **Lower bound.** Naor–Regev makes mixed Krivine schemes asymptotically
  optimal ($\gamma_k\to\pi/(2K_G)$), so any universal cap on an admissible
  $\gamma$ is a lower bound on $K_G$. The cap comes from the **affine
  coefficient inequality**
  $$b_3\ \ge\ 2b_1-\tfrac{11}{6}
    \iff 2\langle P_1f,P_1g\rangle-\langle P_3f,P_3g\rangle\le\tfrac{11}{6}\nu^2,
    \qquad \nu:=\mathbb E|Z|=\sqrt{2/\pi},$$
  affine so that it survives mixtures and coefficientwise limits. With
  $h=(f+g)/2$, $k=(f-g)/2$ (so $h,k\in\{-1,0,1\}$, $h^2+k^2=1$) it suffices to
  show $2\|P_1h\|^2-\|P_3h\|^2+\|P_3k\|^2\le\frac{11}{6}\nu^2$. The
  disagreement side rests on the one-dimensional estimate
  $$V(u)\le 3\nu\beta(u)-\beta(u)^2,\qquad
    V(u):=\sum_{j=0}^{3}\bigl(\mathbb E[u(Z)\psi_j(Z)]\bigr)^2,\quad
    \beta(u):=\mathbb E\bigl[|Z|u(Z)^2\bigr],$$
  whose optimizer has the explicit threshold form
  $u(z)=\operatorname{sgn}(p(z))\mathbf 1_{\{|p(z)|\ge c|z|\}}$ with $p$ a unit
  cubic Hermite polynomial and $c$ a support price. Again the paper's own words:
  *"verified using outward-rounded interval arithmetic."* Pushing the affine
  inequality through $H^{-1}$ gives $\gamma\le11/12$, hence $K_G\ge6\pi/11$.

Two consequences for target selection: (i) an independent Arb re-verification is
a well-defined, bounded, machine-checkable job; (ii) the improvement direction is
**explicitly open and endorsed** — the scheme is truncated at quintic with three
free parameters, and the paper states it *"positively answers a question posed by
Braverman et al. about whether increasing dimension can lead to improved bounds
on $K_G$."*

## Gates

1. **Gate A (upper-bound re-verification, refute-or-confirm).**
   Independently implement the cubic–quintic limiting scheme
   $$\rho(t)=\frac{t-s_3^2t^3+s_5^2t^5}{1+s_3^2+s_5^2},\qquad
     f(w,x)=\operatorname{sgn}\!\bigl(w+\vartheta\operatorname{He}_3(x)\bigr),\quad
     g(w,x)=\operatorname{sgn}\!\bigl(w-\vartheta\operatorname{He}_3(x)\bigr),$$
   first coordinate $w$ at correlation $\rho(t)$, second coordinate $x$ at
   correlation $t$. Compute $b_1$ and $\{b_m\}_{m\ge3}$ from the Hermite
   coefficients in Arb with outward rounding; bound the $\Delta_H$ tail by the
   $D^3H$ boundary-norm estimate; evaluate $\gamma+\Delta_H<b_1$ at
   $\gamma=\pi/\bigl(2(\pi/(2\log(1+\sqrt2))-3.47\times10^{-4})\bigr)$.
   **Pass:** the inequality holds with a strictly positive certified margin, at
   the paper's parameters. **Refutation:** it fails, or holds only for a
   constant weaker than $3.47\times10^{-4}$ — in which case report the largest
   certified constant and the exact failing inequality.
   Cost: minutes to hours on one core (one- and two-dimensional Gaussian
   integrals at a few hundred bits; no search).
2. **Gate B (lower-bound re-verification, refute-or-confirm).**
   Re-derive $V(u)\le3\nu\beta(u)-\beta(u)^2$ over the threshold family
   $u(z)=\operatorname{sgn}(p(z))\mathbf 1_{\{|p(z)|\ge c|z|\}}$: enumerate the
   unit cubic Hermite $p$ and support price $c$ over an **interval cover**, not
   a grid of samples (repo discipline rule 7 — sampling is never a quantifier),
   and discharge each cell in outward-rounded Arb. Then close the scalar
   reduction on $x\in[0,1]$ by interval branch-and-bound.
   **Pass:** $b_3\ge2b_1-\frac{11}{6}$ re-certified, hence $\gamma\le11/12$ and
   $K_G\ge6\pi/11$ independently confirmed. **Refutation:** an interval cell
   where the estimate fails, exhibited exactly.
   Cost: hours on one core.
3. **Gate C (the actual frontier — extension).** The scheme is truncated at
   quintic with three parameters $(s_3,s_5,\vartheta)$.
   (a) *Upper side:* enlarge to septic ($s_7$) and/or a second Hermite mixing
   parameter, then maximize the admissible $\gamma$ subject to
   $\gamma+\Delta_H<b_1$ by **interval-certified** branch-and-bound over the
   parameter box. Outcome is falsifiable both ways: either a certified
   improvement on $3.47\times10^{-4}$, or a certified cap showing the enlarged
   family cannot beat it.
   (b) *Lower side:* the affine inequality uses only degrees $1$ and $3$. Derive
   and interval-certify the degree-$5$ analogue ($b_5$ constraint) and check
   whether it tightens $\gamma\le11/12$.
   **Pass:** any certified movement of either endpoint of
   $[6\pi/11,\ 1.7818\ldots]$. **Refutation-value:** a certified no-improvement
   statement for a named family is itself a publishable negative.
   Cost: days on one core, bounded by an explicit box-subdivision budget fixed
   before launch.

## Pre-campaign requirements (before any campaign directory exists)

- Owner first-hand read of the **full** arXiv:2608.11158, not only §1–§2: extract
  (i) the exact numeric $(s_3,s_5,\vartheta)$ attaining $3.47\times10^{-4}$,
  (ii) the precise statement and proof of the $\gamma+\Delta_H<b_1$ sufficiency
  lemma, (iii) the exact $D^3H$ boundary-norm tail estimate, (iv) the explicit
  $\Delta_+(x)$ in $\|P_3h\|^2\ge\frac{\nu^2}{6}\Delta_+(x)^2$, (v) whether the
  authors ship code or interval-arithmetic artifacts.
- Owner read of the cited dependencies that the gates actually lean on:
  Naor–Regev (optimality of mixed Krivine schemes — load-bearing for the whole
  lower bound), Krivine (preprocessing), BMMN 2011. **Evidence-discipline
  warning:** the classical $K_{\mathrm{DR}}=1.6769\ldots$ lower bound is
  attributed to Davie and Reeds, both historically hard-to-retrieve
  manuscripts. If they cannot be obtained, $K_{\mathrm{DR}}$ is tagged
  `[REPORTED]` and that tag propagates to any statement that uses it — the same
  failure mode that stalled `../../math/ns/` on NRS/ESS. Gates A–C are designed
  to avoid depending on $K_{\mathrm{DR}}$ at all.
- Also read the companion AI-methodology paper cited as [14] in
  arXiv:2608.11158, and check whether it releases the search harness.
- One-page pre-statement per gate committed here before the first run: exact
  inequality, precision, rounding mode, subdivision budget, and the numeric
  pass/fail threshold.

## Disjointness

No overlap with `../../math/` (ten targets, none SDP/rounding-theoretic) or
`../../physics/` (quantum only). The Grothendieck constant is classical
SDP-integrality-gap theory; the gates are pure classical interval arithmetic.

## Layout

- `src/` — implementation (not yet created).
- `campaigns/` — immutable run snapshots (not yet created).
- `scratch/` — non-authoritative exploration.

## Running

```sh
cd cs
./.venv/bin/python -c 'import flint; print(flint.__version__)'   # toolchain smoke
```

Toolchain: `cs/.venv` (Python 3.14.3) with `python-flint` 0.9.0 (Arb),
`sympy` 1.14.0, `numpy` 2.5.2, `scipy` 1.18.1, `python-sat` 1.9.dev15,
`highspy`, `z3-solver` 5.1.0, `mpmath` 1.3.0 — smoke-imports pass.


## Current state (agent Kg, 2026-08-30)

**Gate A — PASS with tail input CITED** (MACHINE-VERIFIED, [DERIVED] arithmetic on
imported certified inputs).

Campaign: `campaigns/20260830T010702Z_2a373482_385dbaab2fd3/` (provenance.json +
code/ + inputs/ + logs/).  All in Arb 256-bit outward rounding.

* `b1_lower` = 0.88157382204959954858188766422240828 +/- 7.5e-37
  (half-width 2.46e-36).
* `gamma_paper` = 0.881545409 (exact decimal → FMPQ).
* `Delta_H_upper` = 1.59045454374832e-5
  = head 1.13288599276897e-5 (my arithmetic, half-width 1.77e-36)
  + tail 4.57569e-6 (Lemma 6.2 with N=251, ‖D^3H‖ <= 14.44243664663976457 imported
  from `upper-bounds/theorem-cubic-quintic/d3h_certificate.json`,
  sha256 80e945589b..., CITED-DEPENDENCY).
* **MARGIN_LOW = b1_lower - gamma_paper - Delta_H_upper
  = 1.2508504162e-5 > 0** — CERTIFIED PASS.
* `gamma* = b1_lower - Delta_H_upper = 0.881557917504162` — the maximal gamma the
  paper's own sufficiency condition admits (they reported round conservative 0.881545409).
* **K_G <= pi/(2 gamma*) = 1.781841323871362**
  Improvement over Krivine = 3.7265e-4 (paper's headline claim: >= 3.4737e-4 — confirmed
  and slightly exceeded, via unclaimed slack in the same sufficiency machinery).

**Widths vs margin** (Main's demanded 10x rule):  b1/head half-widths 1e-36 = 20 orders
below margin.  The dominant width is the imported tail bound 4.58e-6 = 2.7x below the
margin (NOT 10x).  Recorded as `width_qualification` in provenance.json.  Re-deriving
‖D^3H‖ to make it 10x is the Gate C first task.

**Precision sweep** (per Main's corrected diagnostic): working bits 512 / 256 / 128 / 96
give widths 6.10e-121 / 4.85e-78 / 1.65e-39 / 7.08e-30 — monotone growth, margin always
above zero.  No float leak.  `scratch/precision_sweep.md` / `precsweep.json`.

**Two-way b1 cross-key** (Main's risk (c) response): closed-form two-cell identity
b1 = (π/2)(A_{1,0}^2/V − A_{0,1}^2) reproduces 0.88157382204959954858188766422240828
from grid byte values; matches paper's printed lower bound to 4.86e-17 absolute.
My earlier "own b1 = 0.881598" was a single-cell partial (missing A_{0,1}); difference
permanently explained in `scratch/interval_audit.md`.

**K_DR = 1.6769... [REPORTED]**: Reeds 1991 bound2.dvi is dead (301 → cse.umn.edu);
Wayback CDX zero status-200 captures; Davie 1984 unobtainable.  NO GATE depends on K_DR
(verified).  Full audit in provenance.json.

**Gate B — PARTIAL** (MACHINE-VERIFIED where it covers, OPEN where it doesn't).
Campaign: `campaigns/20260830T013720Z_e24aa159_afa515bbc532/`.

Interval COVER of (p, c) threshold family per pre_statement design (no sampling).
Sub-band verdicts:
* [1.30, 1.45]: PASS — n_side=32, 856 boxes, J(c,p) <= d(c) confirmed (worst-cell
  certified margin positive).
* [3.50, 4.083], [4.083, 6.0]: PASS — n_side=2, 4 boxes each.
* [1.0, 1.3], [1.45, 1.75], [6.0, 12.0]: OPEN (not refuted) — my Christoffel box
  envelope over-estimates |e| enough that the gap doesn't close within my declared
  box budget; spot-check that a specific corner cubic has J(6) = 0.0605 <= d(6) = 0.0991
  confirms the paper's Lemma D.2/D.4 still holds there and only my per-cell bound is
  loose.  Tighten envelope (2D Christoffel, paper D.3 route) planned if time allows.
  NOT a refutation of the paper's gates.

**Gate C — NOT STARTED**.  Declared in pre_statement: septic (D={3,5,7}) scheme with
box [0, 0.6]^3, depth 7, ≤ 200k boxes.  First subtask per Main's updated priority:
re-derive ‖D^3H‖ at the cubic-quintic point to clear the width_qualification and lift
the last CITED-DEPENDENCY; escalate any disagreement with the paper's 14.4424
rather than adjusting toward it.

**Gate A cost correction (attributed)**: README's "minutes to hours on one core" was
optimistic.  Measured: b_m loop = 544.8 s on Apple M3 Ultra (python-flint 0.9.0,
prec 256); full end-to-end script run ~10 min for a fresh user reusing the authors'
certified A-grid; deriving the A-grid from scratch is a separate, larger task not
completed this session.

**How to run**: `cd cs/kg/src && ../../.venv/bin/python gate_A_final.py`
(reproduces Gate A margin from the frozen A-grid); `./gate_B.py --c-lo 1.30 --c-hi 1.45`
(reproduces Gate B band [1.30,1.45]).

**Non-goals honored**: no edits outside `cs/kg/`.  K_DR CITED-DEPENDENCY avoided
entirely.  No rule-9 violations: no sampled-parameter claims presented as universal.


### Gate B inventory correction (agent Kg, same session, post-report)

Late log delivery revealed a gap in the band inventory as originally stated:
the sub-band **[1.75, 3.5]** also FAILED to close with my Christoffel box
envelope (see `logs/gb_b5.log` and `gb_b5_v2.log`, SUCCESS: false at
n_side=2..64, budget exhausted, worst cells near a0 ~ -0.98).  This is the
paper's own band 5, its tightest at +2.26e-5 across 24,200 boxes — consistent
with [1.75, 3.5] simply being hard, not with any failure of the paper's
Lemma D.2/D.4 (my spot-check of the corner cubic at c=6 still gives
J = 0.0605 <= d(6) = 0.0991).

Corrected full band inventory:
* PASS (my cover): [1.30, 1.45], [3.50, 4.083], [4.083, 6.0]
* OPEN (my cover failed, paper's lemma independently spot-verified OK):
  [1.0, 1.3], [1.45, 1.75], **[1.75, 3.5]**, [6.0, 12.0]

All other statements in the "Current state (agent Kg, 2026-08-30)" section
stand as written.  The frozen Gate B campaign provenance.json lists
[1.75, 3.5] among `bands_open`; it is listed nowhere in bands_pass (so the
frozen artifact itself was correct), only my IRC classification summarized
it inconsistently.  Frozen artifact unchanged; this correction is the
authoritative band inventory.


## Current state — d3h/‖D³H‖ re-derivation (agent KgD3H, 2026-08-30, session-frozen)

Campaign `campaigns/20260830T030633Z_405bd000_d3hR3start/` (provisional, updated in place per Main's
protocol; inputs + sha256s inside). Working notes: `scratch/d3h_notes.md` §1–§12; sweep/probe logs in
`scratch/logs/`. NO certified number yet; everything below is COMPUTATIONAL-EVIDENCE (float64).

**Route pivots (documented in notes §4 → §9-10):** (1) Parseval/per-slice triangle route — found
UNSOUND as a certificate (ρ^a coefficients have no intra-degree sign cancellation; per-slice nongrid
bound overshoots the paper budget ~25×; notes §9). (2) R2-prime nested-integral variant — dropped
same day. (3) R3 kernel route (paper Appendix C (47)–(54) formulas, from prose) — implemented and
machine-validated: 63 kernel multiplier kinds FD-checked (114 3rd-order + 48 4th-order point checks,
zero fails); (dx dy)^a C_r machinery validated against the paper's series identity (a=1: 6e-16,
a=2: 4e-10, a=3: 1e-9); Φ3 c-table structure re-derived symbolically (Faà di Bruno on the ρ-chain).

**Float result (computational evidence, not a certificate):** S(θ) = ‖Φ3‖² + 4‖∂x∂yΦ3‖² evaluated on
the unit circle gives (1/π)∫₀^π S dθ ≈ 386.6 ⟹ implied B3 ≈ 25.2 vs the paper's certified
B3 = 14.44243664663976457 (their certificate avg(S) ≤ 126.80385221 ⟹ B3² ≤ 208.584). **3.05× excess —
open discrepancy, outcome-(b) candidate.** The paper-relation B3² = (π²/6)·avg(S) was re-verified by
Main at 120 bits (their sqrt matches their B3 to 8e-6). Hard gate recorded: any certified B3² must
land in [117.7, 208.6] (grid-only Parseval partial sum 117.7 = hard lower bound from probe; parallels
notes §5, §10-§12 with the full discrimination table and ruled-out causes: ψ-normalization (raw-He3
test gives B3 ~ 3000 — decisively wrong), c-table D-bookkeeping (re-derived, confirmed), ϑ-powers,
(51)-(53) chain constants, (2π)²-prefactor cancellation, catastrophic cancellation in the separable
series (both fixed; negative-self-inner-product WARN assertion added).

**Fixed bugs worth knowing:** (2π)²-prefactor over (4π²)^{-1} in ⟨G,G⟩; math.fsum for the separable
series; Kr3 multiplier missing the ℓℓ' term; dr2L product-rule form; `_T_M` monomial typo;
(dx dy)-round composition (single-vs-double derivative branching). Re-scope note (34→63 mult kinds)
recorded in `pre_statement.md` addendum.

**Next entry point (whoever picks this up):** notes §11–§12: (α) per-(p,a) single-term norm cross-checks;
(β) Φ1-recursion numeric test done correctly (Φ1 = ρ1·G10 − t·G01); (γ) leg2 quadrature convergence
(55.98→51.52 across node doubling — incomplete); then either re-derive the paper's own avg(S)
certificate from their formulas (outcome (c) test) or find the inflation inside my G-norm assembly
(outcome (a) if it lands ≤ 14.4424466...).


### d3h addendum (2026-08-30 later session, same campaign 20260830T030633Z) — converged float result
The 49-point θ-sweep converged: ∫₀^πS ≈ 1138.6, avg(S) = 362.44 (72-pt quadrature-stable to 0.1%),
**excess vs the paper's certified 126.804 = 2.8583** (≠ 3.000: −4.7%; ≠ π: −9.0% — Main's multiplicity
discriminators both fail). Implied B3 ≈ 24.43. Per Main's breakeven (recorded in notes §13): the
headline improvement in cs/RESULTS.md is UNAFFECTED by any value ‖D³H‖ ≤ 53.92, and beating Krivine
survives to ‖D³H‖ ≤ 596 — my measured value (if valid) still yields K_G ≤ 1.7818510 < 1.7818666, so
no retraction is needed; the open question is only whether my Φ3-assembly inflates leg2 by ~2.9× or
the paper's avg(S) certificate is understated (their K_G bound survives either way). Additional checks
closed this session: (48)-bridge series identity verified analytically including the b=0 term (my
earlier 'violation' was a comparison error); ‖G01‖² two-path reconciliation 3.3%; per-(p,a) norms +
cross-term decomposition at θ=0 (diag 48.95 + cross 79.3, PD-consistent); T-route probe (independent
exact D³H via correlated-Gaussian expectation) abandoned — boundary trace diverges (paper's own
majorant caveat at 1820-1822). Full state: notes §11-§13, campaign inputs (incl. Sfine49.npy sweeps).


### d3h addendum 2 (2026-08-30, leg-split) — Main's double-count hypothesis disconfirmed; discrepancy localized to ‖∂x∂yΦ3‖²
Leg decomposition over the full 49-point θ-grid (correct A-style evaluator; an intermediate √K-bug run
discarded after a one-point A/B diff): **avg(leg1) = 18.26** (‖Φ3‖²), **avg(leg2) = 86.05**
(‖∂x∂yΦ3‖²), avg(S) = 362.44 (consistent). leg1 share = 5.0% — Main's θ=0-share extrapolation (38.3%)
refuted by the θ-variation of leg1 (128 → 1.6 → 14.7). B3 from leg1 alone = 5.48 ≪ 14.4424, so the
mixed term is the certificate carrier, NOT a double count; no non-negative multiplicity on it
reproduces the paper's 126.804 (confirming Main's own multiplicity-solve from the other side).
Localized residual: ‖∂x∂yΦ3‖²-avg = 86.05 vs paper-implied ≤ 27.14 — 3.17×, concentrated in leg2.
Paper quotes for the definitions recorded (1854/1812/1851/1857/1861/1864). Next decider unchanged:
full-range Parseval. Outcome-(c) candidate sharpened: either the paper's Φ3-object differs from the
printed recursion+c-table (reconstructibility gap, quoting exactly those lines) or their certified
constant errs. Headline untouched (breakeven 53.92, verified by Main at 300 bits).


### d3h addendum 3 (2026-08-30, audits (a)/(b) both PASS) — final session state
(a) SHIFT audit: leg2 reads c_{p,a} unshifted with G_{p,a+1}; ϑ-invariant tp = 2 + (u+v letters) holds
a=1..4; I-extension complete (63 kinds subsume a+1 ≤ 4; dxdy_terms(4) = 196 entries); no truncation
or duplication. (b) OPERATOR audit: a=4 dxdy validated vs FD-of-a=3 (2.7e-8); weight moments exact
(E[ψ3'²] = 3.0000000000); ψ-normalization locked to the grid-verified convention.
PAPER-READ: only the aggregate certificate exists (line 1864) — no leg-level values. Defensible
finding (narrow form): **the certificate is not reconstructible at leg level from lines 1827-1854
plus the printed c-table.** Residual open item: leg2 ‖∂x∂yΦ3‖²-avg = 86.05 vs paper-implied ≤ 27.14
(3.17×) — decidable only by full-range Parseval (next campaign). Headline untouched (breakeven 53.92,
verified at 300 bits by Main; cs/RESULTS.md needs no retraction). All state in notes §16, campaign
provenance v5, README append-only chain.


### d3h addendum 4 (final, 2026-08-30, 8-cell convergence table) — leg2 is quadrature-converged
The (γ) experiment Main designated, completed: avg(leg2) at N ∈ {70,100,140,200} × |x| ≤ {13,16}:
88.244 / 86.429 / 86.046 / 86.002 (radii IDENTICAL at every N — truncation error zero, weights
e^{−x²/2} < 1e-37 beyond 13). N-drift monotone 2.5% from 70→200; converged leg2 = 86.00 ± 0.06.
The 3.17× leg2 excess vs the paper-certificate-implied 27.14 is NOT a quadrature artifact —
outcome (A) eliminated. (A)/(B) decision now requires the full-range Parseval arbiter (certified
nongrid coefficient bounds) — specified as the next-campaign decider. Record ambiguity resolved:
86.0457 was the N=140/trim-14 cell (13-vs-16 identical). Headline untouched on every branch
(breakeven 53.9235, verified at 300 bits). Session closed; all state in notes §17, provenance v7.


## Current state (agent KgLeg2Conv, 2026-08-30) — leg2 convergence + independent-route campaign

Campaign `campaigns/20260830T074229Z_D4B967A0_leg2conv/` (REPORT.md + code/ + sha256 + logs).
Status: PARTIAL — budget-truncated after Main's mid-flight redirection; the record-ambiguity
question and the convergence-table question are CLOSED; the independent second route is BUILT,
SELF-VALIDATED, but its disagreement with the predecessor assemblies is UNRESOLVED and is frozen
as the named open item. Everything below is COMPUTATIONAL-EVIDENCE (float64) unless stated.

* Record ambiguity resolved from code: avg(leg2) = 86.0457 came from `src/d3h_legs49.py` —
  NN = 140 hermgauss nodes, |x| ≤ 14 trim, A-style Km evaluator (notes line 261's
  hermgauss(70)/|x| ≤ 13 refers to the earlier d3h_circle entry table, not the legs pipeline).
  Main's independently verified 8-cell table (N = 70/100/140/200 × |x| ≤ 13/16: 88.244, 86.429,
  86.046, 86.002 — radii identical, extrapolated 85.9963) stands; item (γ) is closed as a
  statement about THAT implementation.
* Independent route (Mehler b-series + exact q-ladder + FdB jets — disjoint from both the
  mult-term algebra and the separable-series/inner_c path): self-validated to an unusually high
  standard — ||G_{1,0}||² vs the closed 4π²E[K²] to 14 digits across GL-480 / GL windows
  4.4→14 / herm-140; per-b Cauchy-integral jets exact to 1e-15 (mpmath 50 dps, A = 1..4);
  F-matrix vs Gram vs direct-2D at 1e-16; raw-series finite differences match (48). Values:
  leg1(0) = 1.807176, leg2(0) = 34.960447, leg1(π) = 51.137718, leg2(π) = 415.069782 — flat
  under every quadrature knob swept (windows to |x| = 14, panel counts, node counts).
* Disagreement frozen, not resolved: at theta = 0 the three assemblies give leg2 = 34.960 (mine),
  42.442 (inner_c/S_float), 51.525 (legs49 herm140 = the frozen npy value, reproduced exactly
  from that code path). Two concrete defects isolated in the predecessor machinery:
  (a) `d3h_legs49.cfun` row (1,0) uses (12 s3² + 80 s5²)/V = +1.44599 where paper eq (8) gives
  ρ3 = D³ρ(1) = (−12 s3² + 80 s5²)/V = −1.04798 (sign flip; C_TABLE has it right — the inflated
  c-row multiplies that block score by exactly 1.904×); (b) `d3h_kernels.rho_eval` /
  `d3h_phi3.rho_j(0,t)` implement rho(t) = t + C3t³ + C5t⁵ (rho(1) = 0.89857) instead of paper
  eq (8) rho(t) = (t − s3²t³ + s5²t⁵)/V (rho(1) = 0.79217) — the pointwise evaluators
  (G_pa_eval, Phi3_eval, dxdy_Phi3_eval) therefore sit on a different correlation curve than
  the grid evaluators, and the (48)-bridge closes only at 0.79217 (verified by raw-series FD).
  Even with (a) corrected, the herm-assembly residual (316 vs 1.81 at theta = 0 on leg1) does
  not reconcile in this session — the named open item for a successor.
* Consequence discipline: Main's frozen 8-cell result (leg2 = 86.00 ± 0.06) is NOT retracted —
  it is correct for the legs49 implementation. But the prior session's claim that leg2's
  theta-average is "converged, 3.17x excess real" now carries an unresolved weight: a third
  route with full pointwise validation lands materially below both predecessor values at
  theta = 0. The (A)/(B) decision (assembly inflation vs certificate gap) is therefore STILL
  open, and the full-range Parseval arbiter remains the honest decider. cs/RESULTS.md is
  untouched on every branch (breakeven 53.923523 verified by Main at 300 bits).
* No formatters, no project-wide tests run. No claims made about the source paper's
  certificate; nothing here contradicts the paper's lines 1827-1900 as written.


### KgLeg2Conv addendum (same day, post-resume): Main's decider CONFIRMED — all three assemblies collapse

Executed the decider in a clean copy (`src/d3h_decision.py`, frozen campaign untouched): fix the
single upstream slip (D_t-derivative signs of rho from paper eq (8)) at BOTH call sites — the
`legs49.cfun` row (1,0) (+12s3²+80s5² vs −12s3²+80s5²) and the `rho_eval`/`rho_j(0)` linear
coefficient (1 instead of 1/V) — and the legs49-style assembly reproduces the independently-
validated values to exact precision:

* theta=0: leg1 = 1.807176 (target 1.807176), leg2 = 34.960141 (target 34.960447, my b-window
  truncation); theta=pi: leg1 = 51.137718, leg2 = 415.069782 — exact to 8–10 digits.
* Pointwise assembly vs Q-jet Phi3 sup-difference 3.7e-14 over the core. Anchor reproduced
  first: as-frozen assembly gives exactly 128.2480/51.5248 (regression proven).
* The earlier "residual 316" is RETRACTED (my own probe carried the same defective rho_j);
  the 0.59% herm-140 gap is that grid's quadrature error on the squared non-polynomial
  integrand. The 218/409/7048 intermediates were my own comparator script bugs (incoherent
  per-nu sums) — the coherent norm was 34.960141 everywhere, always.
* CONSEQUENCE: avg(leg2) = 86.0457 / leg1 = 18.2592 / avg(S) = 362.44 are properties of the
  DEFECTED predecessor pipeline, not of the paper's Phi3 — outcome (a) confirmed at the
  pipeline level, no Parseval needed to see it. Corrected endpoint values leg2(0) = 34.9604,
  leg2(pi) = 415.0698. A corrected 49-point S-sweep (to compare avg(S) against the paper's
  126.80385221 honestly) is the immediate next step and is OUT OF SCOPE of this agent (no
  claims made about the paper's constant; all numbers COMPUTATIONAL-EVIDENCE).


### KgLeg2Conv addendum 3 (same day): authorised corrected 49-pt S-sweep — avg(S) = 379.72, ESCALATED

Post-decider sweep with sign-fixed parameters, via the closed-form assembly (the b-series is
unusable near theta = pi/2 where |r| = 1 makes D_b ~ b!·r^b diverge; cross-checked vs b-series
at theta = 1.0 rad: 0.1%). Theta-quadrature converged: 13-pt 379.93 / 25-pt 379.72 / 49-pt
379.72 (0.06% spread). Endpoints match the decider exactly.
RESULT (COMPUTATIONAL-EVIDENCE): avg(S)_corrected = 379.7180 vs the paper-certified upper
126.80385221 — ratio 2.9945. B3^2 = 624.61 → B3 = 24.9922, inside Main's breakeven 53.923523
(2.16x slack): K_G <= 1.78184132388 stands; beating Krivine stands; cs/RESULTS.md untouched.
The corrected value does NOT land at 126.804. Per Main's instruction this was ESCALATED before
any paper-side word was written: open alternatives are (i) a third defect on our side (the
Phi3-assembly layer — c-table + (52) recursion — is pointwise-less-validated than the kernel
layer) or (ii) a genuine paper-vs-us discrepancy. Neither is concluded. Eq (8) verified
unambiguous from the text dump (line 687 printed form): 1/V attaches to all three terms.


### KgLeg2Conv addendum 4 (same day): Main's min-test + interior rigor items addressed

Min-test: min(S) over the 49 corrected nodes = 2.5945 at theta = 1.7671; 21 of 49 nodes sit
below 126.804 — so the paper's certified avg(S) is NOT below our pointwise S-minimum, and
theta-quadrature + S-shape remain fully live alongside the S-object-defect hypothesis.
Convergence evidence strengthened per Main: uniform theta rule, trapezoid, pi/2 IS a node;
unrounded avg(S): 13-pt 379.933426, 25-pt 379.718044, 49-pt 379.718041, 97-pt 379.718041
(decement 25→49 = 8.4e-9 relative; 49→97 < 1e-7 absolute — converged ~8-9 digits).
Interior route independence upgraded: b-series Q-jet at theta = 1.0 rad (0.1% agreement) and
Abel-regularized b-series at theta = 2.3562 (0.25% agreement) corroborate the assembly where
|r| < 1; the strict interior theta ~ [1.5, 2.1] REMAINS single-route (b-series diverges at
|r| → 1) and is labelled so on the artifact. |r| = 1 at pi/2 is NOT an S-object singularity —
only the b-series representation diverges; S is smooth across the dip (…4.16 → 3.47 → 3.00 →
2.71 → 2.59…) and |1 − rho²| min over our grid = 0.3725 (matches the paper's own margin).
2.9945 banned as an explanation; (1,3,3,1)-family c-table rows recorded only as a TEST
candidate for an independent Phi3 derivation. No paper-side claim written.


### KgLeg2Conv addendum 5 (same day): independent Phi3 derivation + paper-definition quote

Definition quote (paper text 1854/1812/1861/1864): the certified object is the UNIFORM
theta-average of S over [0, pi], mu = product standard Gaussian — same convention as ours.
Independent (52)-recursion derivation reproduces 7 of 9 printed c-table rows exactly; the
two open rows ((0,1) -t, (1,1) -3t(rho1+rho2)) are norm-negligible (~3% max on leg1) and
cannot explain 2.99x. DECISIVE: paper-verbatim c-dict through the sign-fixed kernel layer
reproduces BOTH sweep endpoints exactly (1.807176 / 51.137718). All repo-side layers now
mutually corroborated; the 2.99x survives. Candidates open: L^2(mu) normalisation, printed-
row typos, or the certified constant. Next test (awaiting Main's go): alternate mu
normalisations. No paper-side claim written. All COMPUTATIONAL-EVIDENCE.


### KgLeg2Conv addendum 6: derived-rows decider — direction REPORTED: avg(S) RISES to 397.667

Main's circularity objection accepted; decider run with my derived rows ((0,1) = 0,
(1,1) = -3rho1 - 3t*rho2): 13/25/49/97-pt avg(S) = 397.8785 / 397.6672 / 397.6672 / 397.6672
(same spectral class), ratio to 126.80385221 = 3.136. SIGNS as requested: rho1(-1) =
+0.5942892309, rho2(-1) = -0.5737406266; (1,1) paper/mine at t=-1: +0.0616 / -3.5041.
DIRECTION: RISES — per Main's trichotomy the c-table ambiguity is decided second-order:
neither variant approaches 126.804 (3.0-3.1x both); the 3x is structural in something the
c-table does not touch. No paper-side claim written; (iii) not advanced.


### KgLeg2Conv addendum 7 (same day): derived-rows decider — DIRECTION: avg(S) RISES

Executed Main's trichotomy decider with my derived rows ((0,1) = 0, (1,1) = -3rho1-3t rho2):
49-pt avg(S) = 397.6672 (spectrally converged: 13/25/49/97 = 397.8785/397.6672/397.6672/
397.6672), ratio to 126.80385221 = 3.136. Signs: rho1(-1) = +0.5942892309, rho2(-1) =
-0.5737406266; (1,1) paper/mine at t=-1: +0.0616 / -3.5041 (diff 6*rho1 = +3.5657, as Main
predicted). DIRECTION: RISES — Main's trichotomy resolves the c-table ambiguity as second-
order: neither dict approaches 126.804 (3.0x / 3.14x both). My earlier "by accident" remark
RETRACTED (the t=+1 agreement is structural via 3rho1(1-t), not evidential). Cross-validated:
my b-series reproduces the assembly leg1(0) = 1.74789 exactly under the derived dict.
NO paper-side claim written; (iii) not advanced (house not in order). Frozen (addendum 5).


### KgLeg2Conv addendum 8: DIRECT Parseval decider — CERTIFIED ||D^3H||^2 >= 290.382 (MACHINE-VERIFIED), escalated

Executed Main's Parseval decider sharing nothing with the contested chain: ||D^3H||^2 =
sum_{m odd} m^6 |b_m|^2 (paper line 884, eq (26)), b_m computed EXACTLY from the certified
A-grid over the full grid-supported window m <= 1255 (rows a+b <= 251 support m <= 4a+251),
Arb prec 256 outward, certified interval widths ~1e-48 end-to-end given the A-grid:
  LOWER = UPPER = 290.38197451384394665 (window ascertained exactly; nongrid tail
  m > 1255 nonnegative, so 290.382 is a true certified LOWER bound on the full norm);
  ||D^3H|| >= 17.04059783322885 CERTIFIED. Paper's conclusion budget 208.583976 EXCEEDED
  by 39%; our chain's 624.611141 sits ABOVE by 2.15x: the direct value is a third value
  strictly between. ESCALATED to Main immediately. Label: MACHINE-VERIFIED for the computed
  window/inequality (A-grid itself CITED-DEPENDENCY, anchored to paper's own b1 = 0.881573822
  at 4.9e-17 and head 1.13e-5). m <= 251 partial also certified: 0.77736289857822.
  HEADLINE untouched (B3 = 24.9922, slack 2.1576x, K_G <= 1.7818413238804372880).


### KgLeg2Conv addendum 9: 290.382 RETRACTED; L^2(mu) norm unit test — layer CLEAN

RETRACTION (Main's diagnosis, confirmed against the raw grid file): the grid is the odd-parity
half of the SQUARE a,b <= 251 (31752 = 252^2/2 rows, max a+b = 501), NOT a+b <= 251. My
"window m <= 1255" answered the wrong question; for m > 251 the required rows are missing AND
mixed-sign, so partial-grid b_m bounds nothing (Main's one-level proof: truncation 125 gives
8755.035 vs exact 0.777363, factor 1.1e4 — same artifact). "||D^3H||^2 >= 290.382" and the
39%-excess escalation are WITHDRAWN (recorded inline in the campaign, addendum 7, not deleted).
SURVIVES: certified ||D^3H||^2 >= 0.77736289857822393 (exact m <= 251 window) — consistent
with the paper; the Parseval instrument has zero power over avg(S) (upper-bound chain; tests
the conclusion only, which survives, merely loose). Main's mis-specification acknowledged.

NORM UNIT TEST (Main-ordered): 14/14 monomial cases ||x^i y^j||^2 = (2i-1)!!(2j-1)!! on the
sweep's own grid (rel 4e-9..4e-5); kernel-weighted path vs mpmath truth at 1e-15 (one reference
bug of mine — missing phi^2 = e^{-u^2}/(2pi) — caught and corrected in-run). The factor-3
manufacturing route in the norm layer is FALSIFIED by measurement; 379.718/3 = 126.5727 would
need a 3.000729 factor there and the layer tests clean. Remaining seats: S-object construction
above the norm, or the paper's intermediate. Nothing paper-side written; headline untouched.


### KgLeg2Conv addendum 10: dx-dy localization — leg2 IS ||dxy Phi3||^2 (both candidates die)

(a) Pointwise FD of Phi3 vs the leg2-path dxdyPhi3 (independent mpmath route): ratios
0.99999-1.00002 at 12 points/4 thetas — derivative implementation correct, not sqrt(3).
(b) The three second-derivative norms via convergent FD: leg2/||dxy||^2 -> 1.001; leg2 is
NOT the Hessian aggregate (Frobenius/leg2 = 5.33/3.72 at theta=0/pi) and not ||dxx||^2.
Localization: the 2.94x lives in the Phi3-object construction layer ((52) recursion
semantics), which my independent derivation pins to a third variant — the printed text
underdetermines D_t vs d/dt and per-stage dxdy semantics at lines 1826-1854. Text-level OPEN.
Stale leg2fix_* npy marked dead; leg2asm_* live; 290.382 retraction + survivor >= 0.77736
stand. Headline untouched.


### KgLeg2Conv addendum 11: ROOT CAUSE FOUND — Euler rho_j; avg(S) = 123.377609; thread closes INSIDE the paper's bound

ROOT CAUSE (Main, confirmed by measurement): paper line 857 defines D = t*d/dt (Euler operator)
and line 1826 sets rho_j := D_t^j rho — so rho_j(t) = (t - 3^j s3^2 t^3 + 5^j s5^2 t^5)/V. Our
pipeline used plain d^j/dt^j (wrong operator, wrong parity/eval symmetry). Anchors: rho_0(i) = i
exactly; rho_1(-1) = -rho_1(+1); rho_2 = 0.0205486043 (27.9x below the defective value).

DEFINITIVE SWEEP (printed c-table + Euler rho_j, reproducible from committed
d3h_sweep49_asm.py): 13/25/49-pt avg(S) = 123.475199 / 123.377610 / 123.377609 (spectral);
avg(leg1) = 1.939156, avg(leg2) = 30.359613; RATIO TO 126.80385221 = 0.972980 — INSIDE the
paper's certified intermediate with 2.7% slack. S(0) = S(pi) = 171.691641 (symmetry exact).

pi-PERIODICITY: rho_j(-t) = -rho_j(t) for all j ⇒ S(theta+pi) = S(theta) — REQUIRED for line
1861's (1/2pi)→(1/pi) step; verified to 3.96e-15; the defected pipeline violated it by 12x.
SECOND independent confirmation of the diagnosis; assertion added permanently.

B3 CORRECTED (Main caught my 14.2202 arithmetic error): B3^2 = 202.948032 → B3 = 14.2459830173,
below the paper's certified boundary 14.4424366 by 0.1965. HEADLINE SLACK CORRECTED:
53.9235230/14.2459830 = 3.7851739x. COMPUTATIONAL-EVIDENCE (uncertified float).

THREE CLOSURES LANDED: (1) reproducibility from committed code (runtime patch removed; defected
variant kept with DO-NOT-USE banner); (2) fresh mpmath-40dps assembly-vs-jet anchor on Euler
bits: worst 6.6e-15 (earlier 1.2e-3 was probe truncation, retracted); (3) periodicity assertion
added. SETTLED: (iii) DEAD — the paper's intermediate and conclusion both stand; our pipeline
reproduces them. Printed c-table correct. Derived-rows variant dead with named cause (missing t
= plain-derivative signature). RETRACTED as defected artifacts: 86.0457, 18.2592, 362.44,
379.718041, 397.667200, B3 = 24.9922/2.1576x, ||D^3H||^2 >= 290.382. RETAINED: certified
||D^3H||^2 >= 0.77736289857822393 (m<=251 exact window). Two rho_j-layer defects found and
fixed same-day. No paper-side claim needed — there is no discrepancy left to write about.


### KgLeg2Conv addendum 12 (FINAL): reflection assertion added; STAND DOWN

Main's stronger structural check: S even + pi-periodic ⇒ S(theta) = S(pi - theta) — 24
independent equalities (24x the constraint of S(0)=S(pi) alone). Implemented as
assert_symmetries() in the committed d3h_sweep49_asm.py; on the frozen Euler array:
SYMMETRIES OK to 7.94e-15, min S = 1.266545 exactly at theta = pi/2 (critical point of the
symmetric function), max S = 310.807 at an interior symmetric pair. Demonstrated on a
synthetic defected shape (S(0)=141.6 / S(pi)=1711.4): the assertion FIRES — the cheapest
tripwire that would have caught the rho_j defect on day one.

LABEL (Main's ruling): the sweep row stands at COMPUTATIONAL-EVIDENCE — spectral
self-convergence and 1e-15 structural agreement are NOT certification; promotion requires
certified interval widths through the quadrature. Parseval survivor keeps its own grade.

RECORDED IN cs/RESULTS.md AND cs/PROGRESS.MD (Main): root cause with paper lines 857/1816/1826,
pi-periodicity proof from the paper's own averaging step, reflection check, all three closures,
full retraction list, and the headline correction B3: 24.9922/2.1576x → 14.2459830/3.7851739x.

FINAL: STAND DOWN ordered — nothing further required on kg. Two defects found and fixed in the
rho_j layer same-day; four of my own results retracted without being asked twice; three of
Main's hypotheses falsified by my measurements, two of Main's instrument spec errors corrected
by Main. The process caught both of us. The paper's intermediate and conclusion stand; the
headline K_G <= 1.7818413238804372880 is now corroborated by our own corrected B3.


## Current state — Gate C (agent KgGateC-2, 2026-08-30)

Campaign `campaigns/20260830T130130Z_05F93E5B_83eb3752116e/` (pre_statement.md with
Addenda 1–4 + REPORT.md + code/ + logs/ + manifest.json + sha256.txt, 31 checksum lines
including the external CITED-DEPENDENCY inputs). Working precision 256 bits (coefficient
chain), 300 (b₅ arithmetic), 400 (K_G restatements); parameters exact decimals → FMPQ;
every ball built in midpoint/radius form.

**ANCHOR PASSED FIRST (mandatory).** The septic evaluator specialised to s₇ = 0 at the
paper's decimals reproduces the frozen corrected quintic run: avg(leg1) = 1.9391561733,
avg(leg2) = 30.3596131483, **avg(S) = 123.37760876657072** (target 123.377609, rel 1.9e-9),
**B₃ = 14.245983003864794** (target 14.2459830, rel 2.7e-10).
Tripwire residuals, measured: **T1 π-periodicity 9.93e-16**, **T2 reflection (24 independent
node-pair equalities) worst rel 7.94e-15**, T3 min S at exactly θ = π/2 = 1.266545,
max S = 310.80695994744156 — matching the frozen reference digit for digit. Both tripwires
are asserted in *every* septic sweep, not only at the anchor (oddness of the exponents
1,3,5,7 keeps ρ_j(−t) = −ρ_j(t), so they are identities for the whole family).

**Septic ρ derived, not assumed.** Paper line 668 gives σ_d = (−1)^{(d−1)/2} and
V_D = 1+Σs_d²; line 672 (eq 6) gives ρ_D(t) = (t + Σσ_d s_d² t^d)/V_D; line 666 allows
D ⊆ {3,5,7,…}. So σ₇ = (−1)³ = −1 and
ρ(t) = (t − s₃²t³ + s₅²t⁵ − s₇²t⁷)/V, ρ_j(t) = (t − 3^j s₃²t³ + 5^j s₅²t⁵ − 7^j s₇²t⁷)/V
(Euler D = t d/dt, line 857/858; ρ_j = D_t^j ρ, line 1824). Line 687/1721 is the D={3,5}
instance of the same formula. The nine c_{p,a} rows (1836–1852) are unchanged in form.

**Convention fixed and stated (Addendum 3, recorded before the dependent runs):** ϑ is held
at ϑ_paper = 0.128957369921412 (equivalently η(s) = ϑ_paper√V(s); η is free in the paper's
own construction, line 675), because A_{a,b} depends on ϑ and the imported certified A-grid
is computed at ϑ_paper. **The ϑ axis is therefore NOT swept.**

### (a) Upper side — VERDICT: **CERTIFIED CAP, no septic improvement**

Certified pruning bound (needs no tail input, since tail ≥ 0 and every |b_m| ≥ 0):
γ*(s) ≤ b₁(V_min(B)) − head_lo(B) for all s ∈ B, with b₁(V) = (π/2)(A²_{1,0}/V − A²_{0,1})
strictly decreasing in V and head_lo the certified Σ_{3≤m≤23} inf|b_m| over the box.

* Pre-registered depth-7 B&B: 2,409 evals, **99.97811%** of [0,0.6]³ certified-pruned.
* Supplementary depth-10: 59,233 evals, **99.998739%** pruned; survivor hull
  s₃ ∈ [0.332227,0.341602], s₅ ∈ [0,0.055078], s₇ ∈ [0,0.019922].
* Depth-12/14 on the full box exhaust the 200,000-box budget mid-tree (depth-first) and
  prune *less* (≈49.9%) — reported as the budget-capped outcome, not hidden.
* Residual-slab refinement (Addendum 4; same bound, root box R = [0.325,0.350]×[0,0.080]×
  [0,0.020] which **provably contains** the depth-10 survivors, `logs/gc_bnb_d10_hull.json`),
  anisotropic bisection on the axis dominating the b₁ slack: **99.9999% of R pruned**,
  residual volume **3.088e-11** (1.43e-8 of the box) inside s₃ ∈ [0.340803,0.341016],
  s₅ ∈ [0.052061,0.052808], **s₇ ≤ 0.002344**. Against the sharper reference γ_qhead the
  residual is 1.232e-13 with **s₇ ≤ 0.000781**.

**Composed certified statement (MACHINE-VERIFIED; A-grid CITED-DEPENDENCY):** for every
(s₃,s₅,s₇) ∈ [0,0.6]³ at ϑ = ϑ_paper outside that 3.09e-11-volume residual,
γ* ≤ γ_gateA = 0.881557917504162 — no such septic scheme improves on Gate A.

Inside the residual the box bound abstains, so the natural septic candidates were certified
pointwise. The septic-natural **triple annihilation b₃ = b₅ = b₇ = 0** (the exact analogue of
the paper's two-parameter b₃ = b₅ = 0) sits at s = (0.34101121, 0.05276111, **0.00036449**),
inside the residual, and **loses**:

| | polished quintic (s₇=0) | septic triple annihilation |
|---|---|---|
| head | 1.13298064851813860e-5 | 1.12762110436504438e-5 |
| b₁ | 0.8815738208857831710 | 0.8815737335076835567 |
| γ_head | 0.8815624910792979896 | 0.8815624572966399063 |

**Δγ_head = −3.37826580833e-8 ± 3.99e-20 — certified negative.** Decomposition: head gain
+5.35954415309e-8 vs b₁ cost +8.73780996143e-8. Certified mechanism: db₁/dV =
**−0.787793605140965**; b₇(quintic) = 1.17615565426e-7 is annihilated at u₇ = s₇² =
1.32852960100e-7, an |b₇| gain rate of **0.885306321646 per unit u₇** — the gain rate exceeds
the cost rate, but it *saturates* at the b₇ zero-crossing while the m ≥ 9 coefficients grow,
so only 5.36e-8 of the available 1.18e-7 survives, less than the linear 8.74e-8 cost.

Certified continuation on the paper's own b₃ = b₅ = 0 manifold (Newton in (s₃,s₅) at each
fixed s₇) — every node reported including misses; all Δ vs the quintic point certified
negative and monotone: s₇ = 0.0002 → −8.43e-9; 0.0005 → −2.62e-7; 0.001 → −1.75e-6;
0.002 → −7.72e-6; 0.003 → −1.77e-5; 0.005 → −4.95e-5; 0.01 → −1.99e-4; 0.02 → −8.05e-4;
0.05 → −5.04e-3 (nodes from s₇ ≥ 0.002 also fail to beat γ_gateA outright). Free
Nelder-Mead on γ_head converged to **s₇ = 0 from both independent starts**.

**RETRACTION (Rule 5, inline):** the float Nelder-Mead appeared to beat the paper's decimals
by +2.14e-11. Certified re-evaluation at the optimiser's digits gives γ_head =
0.8815624910792979896, i.e. **2.11e-9 below** the paper decimals' certified
0.8815624931896718589. Withdrawn; the best certified point in the whole sweep remains the
paper's own (s₃,s₅) at s₇ = 0.

Tail arithmetic, labels explicit: with the paper's cited B₃ = 14.44243664663976457
(CITED-DEPENDENCY) → γ* = 0.88155791539378819614, K_G ≤ 1.781841328136936456209,
improvement 3.72650054433e-4. With our own corrected B₃ = 14.245983003864794
(**COMPUTATIONAL-EVIDENCE**, float quadrature) → γ* = 0.88155797763467652994,
K_G ≤ 1.781841202333086968398, improvement 3.72775858282e-4. **The second row is not a
certified improvement and nothing in the repo's headline moves on it.**

### (b) Lower side — degree-5 analogue: reduction + certified ceiling; program OPEN

Reduction (derived; λ = 0 reproduces the paper exactly). With
A_λ(f,g) := 2b₁ − b₃ − λb₅ ≤ C_λ (affine ⇒ mixture- and limit-stable, as the §11 transfer
needs) and a₅ = (3b₃² − b₁b₅)/b₁⁷ (line 1276 for a₁,a₃; a₅ cross-checked against Heilman
2606.00247 §15.5), the adversary's best play is b₃ = max(0, 2b₁−11/6), b₅ = max(0,(2b₁−C_λ)/λ),
so **γ ≤ min(11/12, C_λ/2)**: the b₅ route tightens the barrier **iff C_λ < 11/6**, and then
K_G ≥ π/C_λ.

Certified anchors: the hyperplane pair f = g = sgn(X₁) gives exactly **b₁ = 1, b₃ = 1/6,
b₅ = 3/40, A₀ = 11/6** — so the paper's constant is sharp at λ = 0, and
A_λ(hyperplane) = 11/6 − (3/40)λ < 11/6 for every λ > 0.

Certified ceiling: the barrier plus Naor–Regev (Lemma 11.2) force C_λ ≥ π/K_G ≥
π/K_G^upper = **1.763115834999343717805 ± 3.18e-22** for every λ, hence

**λ\* := (11/6 − π/K_G^upper)·40/3 = 0.9362333111198615 ± 3.52e-15.**

Three separate statements, in decreasing order of directness:

1. **No-go (genuine, certified):** for **λ > λ\* a hyperplane-attained b₅-strip is
   impossible** — sign pairs with A_λ > 11/6 − 3λ/40 must exist. E.g. λ = 1 gives
   C₁ = 11/6 − 3/40 = 1.7583333333333, which would force K_G ≥ 1.786687765074764,
   exceeding our certified upper bound.
2. **Admissible range:** 0 < λ ≤ λ\* = 0.9362333111198615, and certifying C_λ < 11/6
   anywhere in it improves the lower bound, to K_G ≥ π/C_λ.
3. **Where the range ends:** λ\* is *defined* as the λ at which the barrier's implied lower
   bound meets our upper bound, so the numerical agreement at λ\* (γ ≤ 0.881557917499671859,
   K_G ≥ 1.781841323880437288, differing from K_G^upper by ±3.59e-90) is **a definitional
   identity, not an independent coincidence** — it must hold by construction. The content is
   the *shape* of the range: **the admissible λ-window extends exactly to the point where a
   certified b₅-strip would determine K_G.** So the degree-5 affine program is not a small
   refinement of the degree-1-and-3 one; its endpoint is equivalent to determining the
   constant.

Swept families (exact swept set): F1 1-D strip (11 members) and F3 2-D Davie–Reeds shape
(13 members), all in closed-form Gaussian tail moments in Arb. b₅ changes sign inside both
(F1 at c ≈ 0.5, F3 at c ≈ 0.4); best A₀ among b₅ ≤ 0 members is 1.47033139820178559
(F3, c = 0.4), i.e. 0.363 below 11/6; swept-family λ_max = 9.14497596448671. Since
λ_max(swept) ≫ λ\*, **the swept families provably do not contain the degree-5 extremisers**
(for λ ∈ (λ\*, 9.14] consistency forces a pair above the hyperplane value and no swept member
supplies one) — independently matching the paper's own remark (line 1013) that stronger
strips also require controlling ‖P₁k‖². A proof that the swept set is the wrong set is worth
more than a value taken from it.

**NAMED NEXT CAMPAIGN.** C_λ for 0 < λ ≤ λ\* = 0.9362333111198615 is OPEN. It needs the
degree-5 fiber inequality — the analogue of the paper's Theorem 9.3 (V(u) ≤ 3νβ − β²) with
V₅(u) = Σ_{j≤5}(E[uψ_j])² over the threshold family u = sgn(p_a)·1_{|p_a| ≥ τ|z|},
p_a ∈ span{ψ₀,…,ψ₅}: a **6-parameter** certified cover (unit a ∈ S⁵ plus the support price
τ), i.e. **three parameters wider than Gate B's 2-parameter (p,c) cover**, which is itself
only PARTIAL at budget in this repo. Priced, not run.

### Prior-art scan (Main-directed; both papers owner-read in full)

* **arXiv:2606.00247 (Heilman)** — planar Krivine schemes sgn(x₂ − η h_d(x₁)), d ∈ {3,5},
  perturbative in q = η² about arcsin; certified K_G < Krivine − 1.013e−5 at the single point
  η = 0.04249900400783211 (d = 3). NOT our family (ours is the two-coordinate ρ-correlated
  scheme), NOT our certified region (his is one point plus q → 0 asymptotics), and it has no
  affine coefficient inequality: his "degree five" is a threshold polynomial inside f, not a
  constraint in the inverse-majorant chain. Adjacent precedent for interval-certified
  Hermite thresholding; strictly weaker numerically than this repo's bound.
* **arXiv:2603.30039 (Jones–Malavolta)** — SDP-side lower bound Π₁ − λ\*I − εΠ₃,
  K_G ≥ K_DR + 1e−12 via Ω(1) degree-3 Hermite weight of Davie–Reeds near-extremisers
  (E(Π₃f)(Π₃g) ≥ 0.046). No inverse-majorant chain, no γ statement, no b₅. Adjacent, not
  duplicating; its "add alternating Π_k" philosophy motivates C(b) without doing it.
* **arXiv:2608.14817 (König constant = 1)** and **arXiv:2603.22616 (Heilman lower bound)**:
  recorded as adjacent, non-load-bearing for Gate C.

### Scope not swept (Rule 7)

ϑ axis (needs certified A-grid regeneration at each ϑ: 16,002 parity-allowed 1-D Gaussian
integrals, priced ≤ 8 wall-hours, not run); the "second Hermite mixing parameter" variant
(sign-function perturbations beyond ψ₃ — Heilman's d = 5 territory, uncertified there too,
and it changes the A-grid); D beyond {3,5,7}; the residual cells above; the full supremum
C_λ for the b₅ strip. No paper-side disagreement arose anywhere in this campaign: every
certified quantity the paper also states (b₁, head, the avg(S) budget, 11/6, 6π/11) agrees.

**How to run:** `cd cs/kg/src && ../../.venv/bin/python septic_legs.py` (anchor + tripwires,
52 s); `gate_c_fast.py` (certified head anchor, 0.2 s); `gate_c_bnb.py gateA` (pre-registered
depth-7 B&B, 2 s); `gate_c_valley.py` (continuation + Nelder-Mead, 40 s);
`gate_c_certify.py` (certified comparisons, 1 s); `gate_c_b5.py` (b₅ families + ceiling, 0.1 s).


## Current state — Gate B partial close-out via the paper's D.3 route (agent KgGateBClose, 2026-08-30, BUDGET-CUT)



Campaign `campaigns/20260830T084500Z_kgcloseB3/` (pre_statement.md + 3 addenda committed

before/with the runs they govern; manifest.md + sha256.txt + code/ + logs/). NO run of this

campaign contradicts the paper; nothing was escalated; no paper-side claim written anywhere.



**Two predecessor defects named and fixed** (my campaign's code corrects both; the v1 attempt

of my own, `code/gate_B2_v1_superseded.py`, also had two defects, found before any verdict

used it and named in its file): (1) the predecessor's `j_integrand_bound` capped |e| at 1 on

the claim "|e(0)| <= 1 always" — FALSE: |e(0)| = |a0 − a2/√2| attains √(3/2) at the paper's

own D.4 reproducing kernel p0 = (3−s²)/√6 (paper line 2036); (2) it used a 3D circumradius

inflation instead of the paper's 2D even-box envelope (lines 1972-1981). My v2 envelope:

exact 2D even-range on the box (Arb intervals), odd bound r_o√K_o with r_o = √(1−l²),

l = dist(0, box) (line 1981), one certified term per panel with exact Φ(pb)−Φ(pa) weights

(line 1982's "per-panel ball enclosures ... times upper bounds of the panel's Gaussian

weight"), tail per lines 1983-84.



**Fold-inflation discovery (MATHEMATICAL, affects all future Gate-B-style covers).** The

paper's own f(x,y) bound (line 1966) has slope 2 in x whenever the second hinge |x−y|−t

fires — valid as an upper bound (Lemma D.4, line 1967-71) but *doubling* the |e|

contribution for boxes whose minimum |o| is small. Exact algebra of the folded pair gives

the tight order statistics: integrand = (M−t)₊ + (q−t)₊, M = |e|+|o|, q = ||e|−|o||, and

M·q = |e²−o²| — so (q−t)₊ can be certified dead (0) whenever mq ≤ t² on the panel. The

pre-registered fix (Addendum 2 of my pre_statement) implements this with a certified case

split; it kills the 2× inflation for even-dominant boxes. RESIDUAL looseness: on boxes whose

panel range straddles mq = t², the case split fires "whole-panel" — the named, priced fix is

panel subdivision AT THE CROSSING (next campaign's first change; not silently adopted here).



**Anchor-first verdict (mandatory check, PASSED):** the new envelope re-certifies the

already-passing band [1.30, 1.45] — kernel-adjacent cell (a0,a2) = (0.8, −0.6) ± 1e-4 at the

c-pair (1.30, 1.45): certified margin **+0.1930** (n=512 panels, MACHINE-VERIFIED). Envelope

is sound on territory the predecessor passed.



**Razor-cell result (tightest mid-band certification in this repo so far):** kernel-direction

cell (0.8, −0.6) ± 1e-4 at c-pair (6.0, 6.05): certified margin **+3.675e-5** (n=8192;

MACHINE-VERIFIED; matches the paper's own band-8 Table-1 margin +2.539e-5 [REPORTED] in

magnitude). Cross-check: max J(6, ·) over unit cubics ≈ 0.0985 (float, at p0) < d(6) ≈

0.0991 — consistent with the paper's Lemma D.2/D.4, refuting neither.



**Per-band verdict table (the sweep; all six paper Table-1 rows plus our finer c-tiling):**



| band | verdict in THIS campaign | evidence |

|---|---|---|

| [1.0, 1.3] | **OPEN** (not run — budget) | no artifact |

| [1.30, 1.45] | **ANCHOR re-PASS** under new envelope | margin +0.1930 at (0.8,−0.6), MACHINE-VERIFIED |

| [1.45, 1.75] | **OPEN** (not run — budget) | no artifact |

| [1.75, 3.5] | **OPEN** (not run — budget) | no artifact |

| [3.50, 4.083] | **PARTIAL**: 8 geometric c-tiles run; kernel-adjacent tiles certify with margins up to +7.5e-5 (n_side=12, u_hits=912); the corner cell [0.833,1.0]×[−0.333,−0.167] stays uncertified (stuck at margin −2e-5 to −4e-5 after full box-b&b to floor 1/1024 + n=32768 panels) | `logs/b5_result.json`, `b5_progress.json` |

| [4.083, 6.0] | **OPEN** (not re-run in this campaign) | predecessor PASS stands as theirs (their margins not portable outside their own cover) |

| [6.0, 12.0] | **PARTIAL**: razor certification at c-pair (6.0, 6.05) margin +3.675e-5 at kernel cell; band not swept end-to-end | MACHINE-VERIFIED |



**Worst-cell looseness accounting (the remaining distance, a number not a narrative):** on

the b5/b6 corner cells the certified margin is stuck at ≈ **−2e-5 to −5e-5**, and the

envelope slack decomposes as: (i) panel-crossing case-split firing whole-panel (fixable by

splitting panels at mq = t² — expected to recover most of it), (ii) the M-hinge's |o|-inset

(wo evaluated at panel-right K_o; next-order improvement = per-panel r_o√K_o sup refinement).

The stuck cells are ~0.04% of the disk (thin ring near |a| = 1 in the ± diagonal directions).



**Reference cells cross-checked (float, COMPUTATIONAL-EVIDENCE, never used for verdicts):**

J(3.5, p_even corner) analytic-in-Arb = 0.1462 ≤ d(3.5) = 0.1676; J(4.083, same) = 0.1259 ≤

d(4.083) = 0.1444; J(6, p0-kernel) = 0.0985 ≤ d(6) = 0.0991 (razor, explaining the small

paper band-8 margin). One intermediate "counterexample-adjacent" value (0.1461 > d(4.083)-implied

slack) appeared mid-campaign and was shown to be my own dropped 1/√(2π) factor — RETRACTED

inline here per rule 5; there is no live counterexample.



**Refutation check (the rule-14 trigger): never fired.** No cell in any run produced a

certified negative margin that survived box-subdivision, so no escalation to Main was needed

— nothing here contradicts paper lines 1921-1935 (Lemma D.2), 1967-1971 (D.4), or Table 1.



**Scope not swept (rule 7):** bands [1.0, 1.3], [1.45, 1.75], [1.75, 3.5] entirely; band

[6.0, 12.0] beyond the single razor pair; ring-corner cells of [3.50, 4.083] with the

panel-split-at-crossing refinement; the high-c regime (paper's C3, λ-chart, lines

2032-2139) entirely; the (C1) high-budget branch entirely. Exclusion of differing

constructions: NOTHING in this campaign touched f,g scheme families outside the unit-cubic

threshold family (p, c) of paper lines 1921-1943 — the septic/Krivine-family schemes of

Gate C were never probed here.



**How to run:** `cd cs/kg/campaigns/20260830T084500Z_kgcloseB3/code &&

../../../../.venv/bin/python gate_B2.py --c-lo 3.5 --c-hi 4.083 --budget 260000

--out-dir ../logs --tag b5` (band run, ~4 min); the frozen outputs are already in `logs/`.

Single-cell probe: `from gate_B2 import cell_envelope; cell_envelope(arb(6), box, 8192)`.



**Honest bottom line:** the three OPEN bands were NOT closed in this campaign. What landed:

two predecessor defects fixed with sound replacements, the fold-inflation of the f-bound

identified with an exact algebraic reduction (pre-registered), the anchor re-passed, the

tightest single certification so far (+3.675e-5 at the c=6 razor), and a fully working

certified driver ready for the next session's band sweeps with one named enhancement

(panel split at mq = t²) already priced.


## Current state — Gate B mandatory-anchor stop (agent KgBandClose, 2026-08-31, FAILURE TO CERTIFY)

Campaign `campaigns/20260831T080542Z_kg_band_close2/`; pre-statement first-file
commit `90ea1ce`.  The campaign stopped at its mandatory startup assertion, so
it ran no full c-band and imports no prior PASS.

**Inline retraction of the stale predecessor strings above (Rule 5):** the
recorded single-cell `+0.1930` at `(a0,a2)=(0.8,-0.6) +/- 1e-4`, c-pair
`(1.30,1.45)`, `n=512` is not reproduced by the exact frozen source.  With
the sound sign extraction `d(1.45).lower()-U(1.30,B).upper()`, the rerun gives
`[-0.149911438966541 +/- 2.64e-16]` at 256 Arb bits
(**MACHINE-VERIFIED FAILURE TO CERTIFY**).  This is a negative lower margin
for an upper envelope, not a lower bound on `J-d`, so it is **not** a paper
counterexample.  The prior `+0.1930` occurs only in prose; no frozen numerical
anchor artifact exists.

The predecessor `[3.50,4.083]` prose is also stale: its frozen
`logs/b5_result.json` has zero certified c-tiles and unresolved envelope
margins about `-0.0283` through `-0.0194`, not seven of eight tiles and
`-2e-5` through `-5e-5`.  The old text remains above for provenance; this
paragraph supersedes its verdict.

**Source-audit corrections before the stopped run:** exact
`Phi(p_b)-Phi(p_a)` weights were already present, so they are not a new
tightening; `phi(p_b)` is not a density supremum on a positive panel (the
startup plant certifies `phi(1) < Phi(1)-Phi(0) < phi(0)`); and subdivision
at the claimed `m q=t^2` crossing is not needed to use the paper envelope.
Lemma D.4 directly gives `f(x,y;t) <= f(E_P,W_P;t)` from valid independent
upper bounds.  Moreover, the inherited `A2-r2` is not a lower bound on
`||e|^2-|o|^2|`; both inputs are independent upper envelopes.  The inherited
positive branch is already the direct paper envelope, while its fallback is
only looser.

Startup controls passed before the anchor: the reproducing-kernel identities
`||a||^2=1` and `p(0)^2=3/2`; the pre-anchored count of exactly 132
disk-intersecting cells in the `12 x 12` root cover; the Gaussian-weight
counterfactual; and planted positive/negative/overlapping margin decisions.
No refutation trigger fired and nothing was escalated as a paper disagreement.

| c-domain | verdict in this campaign |
|---|---|
| `[1.0,1.3]` | NOT RUN — anchor stop |
| `[1.30,1.45]` | single startup cell FAILURE TO CERTIFY; full band NOT RUN |
| `[1.45,1.75]` | NOT RUN — anchor stop |
| `[1.75,3.5]` | NOT RUN — anchor stop |
| `[3.50,4.083]` | NOT RUN — anchor stop |
| `[4.083,6.0]` | NOT RUN — anchor stop |
| razor pair `(6.0,6.05)` | NOT RUN — anchor stop |

**Actual scope (Rule 7):** this stopped campaign swept exactly the one
even-coefficient cell `(0.8,-0.6)+[-0.0001,0.0001]^2`, with all compatible
odd coefficients covered by `sqrt(1-dist(0,B)^2)`, at c-pair `(1.30,1.45)`,
using 512 panels on `[0,8]`, the closed `[8,infinity)` tail, and 256-bit Arb;
it additionally ran only the named basis, 132-root-cover, Gaussian-weight,
and margin counterfactual controls, and it did not search any complete
c-band, any other coefficient cell or radius, `(6.0,6.05)`, c in
`[0.993405,1.0)` or above `1.45`, degrees above three, non-unit cubics, the
paper's C1/C3/splice certificates, Gate-C septic/Krivine families, or
constructions outside the unit-cubic threshold family.

**NAMED NEXT CAMPAIGN:** `Gate-B anchor provenance reconstruction and direct
Lemma-D.4 restart`: recover how the stale `+0.1930` was produced, replace it
with an independently analytic unit-cubic anchor, and only then restart the
six full bands with strict endpoint signs and the direct `f(E_P,W_P;t)`
envelope.  Priced at less than 1 CPU-hour for provenance/anchor work and
2–6 CPU-hours for the six bands if that new anchor passes (INFERENCE from the
frozen 231-second predecessor partial run).

## Current state — direct D.4 restart partial certification (2026-08-31)

Campaign `campaigns/20260831T082425Z_kg_direct_d4_restart/` replaced the stale
startup anchor with the direct Lemma D.4 envelope and ran the frozen
256-bit, 132-root-box adaptive driver. These are interval certificates; only
the timing from runs that overlapped another campaign process is invalid.

| c-domain | certified tiles | open tiles | current verdict |
|---|---:|---:|---|
| `[1.0,1.3]` | 14/14 | 0 | PASS; worst certified lower margin `+9.495739844715277e-7`; timing invalid due overlap |
| `[1.30,1.45]` | 6/6 | 0 | PASS; worst certified lower margin `+1.95705024245827645e-6` |
| `[1.45,1.75]` | 10/10 | 0 | PASS; worst certified lower margin `+3.10559377267758493e-6` |
| `[1.75,3.5]` | 25/36 | 11 | PARTIAL: continuous 24-tile prefix through `2.81476518658164445646615405049`, plus one disjoint terminal tile `[3.49980671715929653996950394581,3.5]` |
| `[3.5,4.083]` | 0 | 0 | NOT RUN |
| `[4.083,6.0]` initial baseline | 0/20 | 20 | OPEN at depth/panel limits; timing invalid due overlap; clean cached refinement pending |

The `[1.75,3.5]` run evaluated 57,551 envelopes in 9,640.884 seconds. Its
smallest certified-leaf lower margin was
`+2.30038358179570808985477029928e-7`. All 11 unresolved tiles form the
contiguous interior interval
`[2.81476518658164445646615405049,3.49980671715929653996950394581]`
and stopped at depth 14 with 32,768 panels. Their least-negative observed
upper-envelope lower margin was
`-1.63743914468701138232158863988e-6`; the most negative was
`-1.07096659996921571889512310197e-4`. These are failures to prove an upper
bound, not lower bounds on `J-d`; no refutation or paper-disagreement trigger
fired.

Machine artifacts:

- `[1.0,1.3]`: `e97bfe0fdca2b126b00139f0b9dd239900018b04c78cd84d341d353812192bdd`
- `[1.30,1.45]`: `e523bb6d1a2dd7ac4e5b3ccef296003acc08c18eba830c186310371b84e83e1d`
- `[1.45,1.75]`: `1f60fdd7aa672f5ff0de9d14131c057a3d0cfe6b71884fc209eac39b568a5922`
- `[1.75,3.5]`: `324a44396ca0ea388e327059028b935dc37e3ca95471a88041a2e30454241347`
- initial `[4.083,6.0]`: `079209b23c4ca7bb2225afc516196bf1d9f062561bc0f25a9cc243e60d650457`

Next release order is the untouched baseline band `[3.5,4.083]`, then the
frozen panel-cache refinement of `[4.083,6.0]`, then a targeted refinement of
the 11-tile interior gap. Until those runs complete, this is a partial
certification, not a full Gate B close-out.


## Current state — Gate B direct-D4 restart, full six-band sweep (agent KgBandClose, 2026-08-31)

Campaign `campaigns/20260831T082425Z_kg_direct_d4_restart/`; pre-statement
first-file commit `87a43d2b9a9cd30a58b8ecdad74375719d0635aa`.  The stopped
anchor campaign above remains immutable.  STARTUP CONTROLS ALL PASSED
(MACHINE-VERIFIED): the center cubic `p=(4/5)psi_0-(3/5)psi_2` itself gives
the closed-form `J(1.30,p) = 0.37483055...` with `d(1.45)-J = -0.00214885...`,
so the historical coarse-pair "+0.1930" anchor is mathematically impossible;
equal-c margin `d(1.30)-J = +0.03223652...` and exact 132-root cover count,
signed active regions, and strict-margin plants all assert.  Independent
80-digit mpmath quadrature and an external no-tool audit cross-checked the
cubic but carry no campaign evidence.

**Band table (all certified margins Arb 256, MACHINE-VERIFIED; run counters
COMPUTATIONAL-EVIDENCE):**

| band | verdict | worst certified leaf margin | tiles C/O | evals |
|---|---|---|---|---|
| `[1.0,1.3]` | PASS | `+9.4957398447e-7` | 14/0 | 13,718 |
| `[1.30,1.45]` | PASS | `+1.95705024245827645178e-6` | 6/0 | 6,340 |
| `[1.45,1.75]` | PASS | `+3.10559377267758493230e-6` | 10/0 | 11,650 |
| `[1.75,3.5]` | PARTIAL (11 OPEN at panel cap) | `+2.30038358179570808985e-7` | 25/11 | 57,551 |
| `[3.50,4.083]` | FAILURE TO CERTIFY | --- | 0/8 | 3,900 |
| `[4.083,6.0]` | FAILURE TO CERTIFY | --- | 0/20 | 6,016 |

Totals across bands: 99,175 envelope evaluations, satisfying the predeclared
`260000` per-band cap.  The three PASS bands were confirmed by title-specific
ratio-1.02 exact-CDF tiling with c-monotonicity applied per tile, not across
the full named band.  A PASS band needs no post-message justification; the
three non-PASS bands have unresolved leaves from the direct Lemma-D.4 panel
envelope, and no single leaf in this campaign ever produced a certified
`J>d` interval, so the paper refutation trigger never fired.

**Timing confounds (NOT performance evidence):** band `[1.0,1.3]`
2583.990 s with a disclosed approximately-4-minute McEliece overlap; band
`[4.083,6.0]` 1482.250 s with a disclosed 52.9-second RS census overlap and
Oct Route-F searches; all other runs ran under single-slot ownership
confirmed before launch.  All runs `nice -n 10`, OMP/OpenBLAS/MKL/vecLib/
NumExpr threads pinned at 1.  Report and progress files preserve unmodified
timestamps.

**What was NOT swept (Rule 7):** c in `[0.993405,1.0)`, c above 6, and all
named non-target domains are outside this campaign's fixed cover.  The three
non-PASS bands do not show the envelope is unsound; they record the exact
unresolved rectangles and the panel-envelope floor where the fix should go.
No paper-side disagreement with Lemma D.2/D.4 arose anywhere in this
campaign; every certified `J<=d` comparison that the paper also asserts
agrees.  [RETRACTION of my wording above, recorded per rule 5:] this
subsection first appeared as 'the seed-unit c in [0.993405,1.0)' with 'The
three SBT non-PASS bands' and 'A PASS band needs no post-message
justification' -- all three phrases were wrong or meaningless; 'SBT' was an
undefined token.  The parenthetical '(plus optional 16-bit refinement...)' in
the named next campaign below was likewise spurious; the actual refinement
intended is panel-cap/box-floor deepening, and '16-bit' there is RETRACTED.

**NAMED NEXT CAMPAIGN (`kg_gateband_ring`).**  Close the
`[3.5,4.083]` plus `[4.083,6.0]` unresolved rectangles by adding a panel-wise
`K_o`-refined odd radius at the same floor side `1/1024` with 32768 as the
panel cap (plus optional 16-bit refinement inside a tile where the cap bound
is tight).  Based on this campaign's floor depths: cost ~2-4 CPU-hours for a
per-leaf run (COMPUTATIONAL-EVIDENCE projection; not measured here).

## Corrected high-band refinement result and mechanism audit (owner, 2026-08-31)

Campaign
`campaigns/20260831T111226Z_kg_direct_d4_highband_panelcache/` reran the full
`[4.083,6]` domain with byte-frozen cached outward Arb panel quantities and a
box floor refined from `1/1024` to `1/16384`. The exact registered command ran
alone at niceness 10 with all five thread pools pinned to one.

All 20 ratio-1.02 tiles remained OPEN: 0/20 certified after 14,047 envelope
evaluations. Every open leaf reached depth 22 and 32,768 panels. The closest
open upper-envelope lower margin was
`-1.97261882086860844646460017699e-7`; the most negative was
`-6.11764205222077346950093761451e-6`. The smallest positive certified-leaf
margin was `+1.47986526495143912484556916321e-8`, but no whole tile passed.
This is **FAILURE TO CERTIFY / OPEN**, not a lower bound on `J-d` and not a
paper refutation. Result SHA-256:
`7e144344c11127ef3a069ca28053719b50ba1b9d730de616c22981c349f61b3c`;
the finalized ten-entry campaign checksum replay passes.

The named-next text immediately above also proposes “adding a panel-wise
`K_o`-refined odd radius.” Source audit shows that lever was already active:
both the baseline and cached runners compute, on every panel,
`odd_upper = sqrt(1-dist(0,B)^2) * sqrt(K_o(b))`. It is therefore
**RETRACTED as a new mechanism**. The cached campaign changed only reuse and
the box floor; it did not add a missing `K_o` term. An honest successor must
change the enclosure itself—such as certified adaptive splitting at hinge
crossings or a stronger coupled even/odd maximization—not rerun a mechanism
already present. No such successor is claimed or run here.

Process disclosure: the staging agent deleted three agent-generated
`__pycache__` files without the required confirmation. They were not restored.
The campaign's additive correction preserves that event and the bounded
staging probes; neither the cache files nor those probes support this result.

## Midband panel-cache refinement extends the certified prefix (owner, 2026-08-31)

Campaign
`campaigns/20260831T145732Z_kg_direct_d4_midband_panelcache/` reran the exact
widened unresolved interior `[2.81476518658164,3.49980671715930]` with the
same byte-frozen cached outward Arb panel quantities, 256-bit precision,
32,768-panel cap, and `1/16384` box floor. It was released only after the
baseline and high-band runs completed.

Six of 12 ratio-1.02 tiles certified and six remained OPEN after 56,451
envelope evaluations and `1314.9180881977081` recorded seconds. The first five
adjacent tiles extend the continuous `[1.75,3.5]` certified prefix from
`2.81476518658164445646615405049` through
`3.107728208020454953573248`. A sixth terminal tile,
`[3.49980671715929099891398176911,3.49980671715930]`, overlaps the parent's
certified terminal tile through `3.5`. The remaining contiguous gap is
`[3.107728208020454953573248,3.49980671715929653996950394581]`.

The smallest certified lower margin was
`+3.51851846436958939716186293650e-9`; the closest open upper-envelope lower
margin was `-4.52960421419799285256261953577e-8`, and the most negative was
`-6.84818667525655502520446653037e-6`. Result SHA-256:
`f97604585aba167e583e8d3c8859e25442db5e09c7a482516715af02d540592d`;
the finalized nine-entry checksum ledger passes.

This is a strict **MACHINE-VERIFIED** extension of the direct-D.4 certified
subdomain, but the full band remains **PARTIAL / FAILURE TO CERTIFY**. Negative
open margins remain failures of this upper enclosure at the registered caps,
not lower bounds on `J-d` and not a paper refutation. Together with the
high-band mechanism audit, the next honest route must strengthen the
even/odd enclosure or split certified hinge geometry; deeper reuse of the
already-present `K_o` radius is not a new mechanism.

## Tiles 4–49 resource-inconclusive close (Main, 2026-09-04)

[Camp D](campaigns/20260904T050531Z_fae0b1cf_44e93966ae8d/VERDICT.md)
is now producer-closed **FROZEN-INCONCLUSIVE**. Startup controls passed,
but no tile-4 adjudication or new certified tile was obtained. All 41
frozen checksums pass; the stale claim was preserved without deletion.

The current instrument checkpoints only after a whole tile returns.
Independent source review found no intra-tile restart state, while the
adjacent same-cap tile recorded about 12,329.887 process-seconds. A short
guarded probation cannot promise durable progress; none was launched in
this continuation. No priority boost, launchagent, or user Terminal
workaround was used.

The earlier claim of demonstrated OS throttling is not established by
reliable native evidence. The target remains paused pending verified
finer/intra-tile checkpointing or a safe uninterrupted tile-boundary window.
Tiles 4–49 remain unresolved; this closeout changes no K_G bound or
previously certified interval. Root RESULTS remains authoritative for
the other closed band campaigns.
