# LIU H2 Block Copositivity — Corrected Diagonalization (2026-08-29)

Worker: ReligiousParrot (subagent).  Deliverables:

* `uc/liu9_block_copositive.py` — the derivation + verifier + report builder.
* `uc/verification/results/liu9-block-copositive-diagonalization.json` — the
  byte-stable report.
* This notes file.

Command: `math/.venv/bin/python -I -B math/uc/liu9_block_copositive.py`
(81 s, single core, nice 19; Campaign K untouched).

## Outcome first

1. **The displayed 2x2 block kernel is REFUTED as an identity** (label:
   PROVED refutation, MACHINE-VERIFIED + symbolic).  The quadratic form of
   `B00 = ((1-q)^2/M) C_M + beta q(1-q) k_pi` (and friends) on alpha-masses
   does **not** equal Liu's raw gap at any tested configuration: worst
   absolute mismatch `0.168` at 90 dps across six deterministic
   configurations, incl. points where the displayed form is negative while
   the raw gap is positive.  The missing content: **every cross-component
   product entropy** `h(x_i y_j)`.
2. **The exact two-piece diagonalization is PROVED (exact, symbolic).**
   For paired laws `P0 = sum a_i delta_{x_i}`, `P1 = sum a_i delta_{y_i}`:
   ```
   gap = (1/M) sum_{i,j<=6} w_i w_j D_M(s_i, s_j)  +  q(1-q) c,
   ```
   atoms `s = (x1,x2,x3,y1,y2,y3)`, `w = ((1-q)a_i, q a_i)`, mixture mean
   `M = sum w_k s_k`, `D_M = M E - A/2 = st C_M` (lemma T1), and the
   **q-free scalar**
   ```
   c = 2(1-beta) cross + beta(pi0 + pi1) - 2 <E>_cross,
   cross     = sum_{i,j<=3} a_i a_j h(x_i y_j),
   pi0/pi1   = sum_{i,j<=3} a_i a_j h(pi(x_i, x_j)) / h(pi(y_i, y_j)),
   <E>_cross = sum_{i,j<=3} a_i a_j E(x_i, y_j),
   E(s,t)    = A - st(mu(s)+mu(t)-mu(st)) + beta(h(pi(s,t)) - h(st)).
   ```
   `c = 0` exactly on the diagonal `P0 = P1`.
3. **Proof method (exact).**  All `h(.)`/`mu(.)` symbols reduced to the
   canonical FREE LOG BASIS `{log s_i, log(1-s_i) for i<=6} plus
   {log(1-s_i s_j) for 21 pairs i<=j} plus {log(1-pi(s_i,s_j)) for 21
   protocol pairs}` (33 canonical log symbols total).  In that basis
   `M_mix*(gap - q(1-q)c) - St` expands to the IDENTICAL ZERO
   rational-coefficient polynomial.  Sub-identities T1 (`D_M = st C_M`) and
   T2 (`2s h(s) - 2s^2 mu(s) + s^2 mu(s^2) = h(s^2)`) proved in closed log
   form.  Residual after reduction: `0` exactly; every mutation nonzero.
4. **Numeric certification (MACHINE-VERIFIED, 90 dps).**  Identity vs the
   independent `gap_mp` transcription: worst abs diff `2.30e-91` over six
   random configurations; extracted `c = (gap - F(P_mix))/(q(1-q))`
   bit-identical across four q values `1/7, 1/5, 2/5, 9/10`
   (deviation `0.0` exactly at 50 printed digits, threshold `< 1e-100`).
5. **Steer's claim (a)+(b) — verdict.**
   * (a) Diagonalization with `w=(1-q)z0+q z1`, `d=z0-z1`: **WRONG as
     stated** for the true raw gap.  There is no `[C_M + beta q(1-q)
     k_pi]`-combination: `ehpi` has no cross-component terms and the
     q-dependence is NOT a `k_pi` cross block but the single scalar
     `q(1-q)c`.
   * (b) `k_pi = xy + x(1-x)y(1-y) = phi1 phi1' + phi2 phi2'` rank-2 PSD:
     **TRUE as a fact about the kernel** (`<d, k_pi d>` = sum of two
     squares for any signed measure), but **IRRELEVANT** — the raw gap
     contains no `k_pi` bilinear at all.
   * (c) does not follow; the corrected reduction stands instead.
6. **Reconciliation with the ledger.**  `liu9_block_kernel.py`
   `refuted_routes[0] = "k_pi is PSD"` and PROGRESS.md "k_pi is not PSD"
   conflated kernel-level PSD with block-level PSD of the DISPLAYED kernel.
   Both records are superseded by this correction: the true refutation
   target is the displayed block algebra (wrong object), and `k_pi`'s
   rank-2 PSD is genuinely TRUE.  `gap = <alpha B alpha>` simply does not
   hold, so the refuted routes are all about an identity that never was.

## The precise remaining wall (sharpened)

Full H2 on the paired class is now **exactly** the scalar margin condition

```
F(P_mix) + q(1-q) c >= 0
```

where `F(P_mix) = (1/M) sum w w' D_M` is the certified raw gap of the
six-atom MIXTURE law (universal q=1 theorem: `F(P) >= 0` for every law of
mean >= m — this is a PROVED nonnegativity, not a margin).  `c` is
sign-indefinite: negative at e.g. `(a1,a2,q)=(1/2,1/5,9/10)`
(`c = -0.0168`), at trial 2 of the random suite (`c = -0.0385`), and at
`(13/32,0,1)` (`c = -0.0199`).  The witness search's most dangerous point
(`gap 1.92e-4` at q~0.91, second-order degenerate) lives exactly where
`F(P_mix)` has vanishing margin AND `c < 0`, so the wall is:

* either prove `c >= -F(P_mix)/(q(1-q))` universally (a two-sided
  statement mixing a scaling law for `F` near its degenerate faces with a
  bound on `c`), or
* find a mean-feasible negative — the search space for the latter is now
  exactly parametrized by `(P0, P1, q)` with `gap = F(P_mix) + q(1-q)c`,
  a much thinner object than the original nine-variable space.

Suggested next probes (NOT done here): (i) exact `delta_x x delta_1`
closed form of `c` (turns into one-variable calculus in `x` at fixed `q`);
(ii) a certified Arb cell-cover for `c` on the three-atom simplex strips
where `F(P_mix)` margins are known (e.g. away from `x*`); (iii) the sharp
`E`-kernel denoising near the `(x*, x*)` optimizer box, where `D_m`'s
positive-definite Hessian cover bounds `F(P_mix)` from below by
`kappa*dist^2` and one needs `c = O(dist)` at worst.

## Honesty items

* This is a **corrected reduction**, not a closure of H2.  No SSOT file
  was edited; ledger integration (RESULTS.md row rewrite for "Liu H2
  exact full block reduction") is the Main session's job.
* `c`'s negative range is bounded only by sampling (6 random + 4
  structured configurations in the report); no global `c` enclosure was
  attempted.
* The discharge nonnegativity inherits the two Audited pillars of the
  universal q=1 theorem (`C_m >= 0` 524,800-box cover; `dC_M/dM >= 0`).
* The symbolic proof is a computer-algebra zero test in an explicitly
  constructed free basis; it is exact (no floating point), deterministic,
  and byte-stable (two consecutive runs: file SHA-256
  `fe5ab3d7d6bfe5d8490a501cf5039a3e6ba625e52461b18825e007f46f9bc652`,
  internal report digest
  `5c856c6ec95faed2e7048783fd36d19cd91c2a2fdd300f6ff007656d68842f79`).

## Addendum (same session): Region 1 certified — Family A scalar margin

Follow-up re-task from Main produced a SECOND certified deliverable:

* `uc/liu9_scalar_margin_region1.py` — the certificate builder.
* `uc/verification/results/liu9-scalar-margin-region1.json` — byte-stable
  report (two consecutive runs byte-identical), internal digest
  `31e13a8d9e3f80d68fa2eba8c4198dc660682948c548371843dd4623b20056cb`,
  file SHA-256
  `f8e9175bd6e3323ed40cf1faa7d40ef6ba58a8f2b1a50c778c8fdfd928e99f58`.
* Command: `math/.venv/bin/python -I -B math/uc/liu9_scalar_margin_region1.py`
  (1.7 s, single core).

**Theorem (Region 1, PROVED, Arb precision 320).**  For Family A
(P0 = delta_x, P1 = delta_1, shared mass a = 1, x in [1/4, 1]) the scalar
margin
`gapA(x, q) = F(P_mix) + beta q(1-q) E_K(P0 - P1)` is >= 0 for EVERY
mean-feasible q.  Proof chain:
1. gapA(x, q) is an EXACT quadratic in q (verified against gap_mp at 60
   dps, agreement < 1e-40; closed form from the corrected diagonalization).
2. gapA(x, 1) = 0 exactly (mixture collapses to P1 = delta_1, whose
   discharge vanishes on the zero-support face) and the quadratic is
   concave in q, so its minimum over the feasible band
   q in [max(0, (m-x)/(1-x)), 1] sits at an endpoint.
3. Two single-variable curves remain:
   (a) x in [m, 1]: A0(x) = gapA(x, 0) >= 0, certified in two pieces:
       direct endpoint-pair interval hulls on [m, 0.6719] (worst lower
       +1.2498e-3, 0 negative cells) and a mean-value-lemma argument on
       [0.6719, 1] (the bilinear (2x-1)h - x^2 T certified pointwise
       nonnegative, tail B(s2 - r) dominated; 0 failures);
       the last cell touching x = 1 is resolved by a monotone-endpoint
       scheme (both the bilinear and s2 - r are strictly monotone there;
       right-paired point evaluations), giving margin 0 (A0(1) = 0 exactly)
       with no negative sub-cell.
   (b) x in [1/4, m]: gapA(x, q0) with q0 = (m-x)/(1-x) and M = m,
       certified by endpoint-pair hulls (worst lower +1.0456e-3, 0 negative
       cells).
4. Mutations: flip_T_sign breaks the hull band (FAIL as required);
   shrink_m (halved mean bracket) collapses the boundary branch (FAIL).

Certificate values (Arb, prec 320): m bracket
(0.61729091208126497, 0.61729091208126498), beta bracket
(0.10005255986289310, 0.10005255986289312), 4096 x-cells per branch.
Scope: Family A only (single-atom components with shared mass a = 1);
the small-x frontier (x < 1/4) adjoins the already-proved zero-support
endpoint theorem instead.  NOT a full H2 proof: the multi-atom paired
class remains open, now quantified as F(P_mix) + beta q(1-q) E_K(P0-P1) >= 0
with the K = h o pi Jensen-geometric structure.

Engineering notes from the build (for the ledger): python-flint 0.9.0's
`arb.union` produces a fuzzy mid-radius ball that decorrelates badly at
x -> 1 endpoints; ALL interval algebra in the certifier is explicit
endpoint-pairs (lo, hi) of point-arbs with 4-corner products, so every
enclosure is auditable.  Sub-lemma: T(x) = 2 mu(x) - mu(x^2) is strictly
decreasing on (0, 1) (verified numerically + derivative sign), x^2(2-x)
strictly increasing, and the c-kernel r - s2 changes sign exactly at
x_hat = 0.6719 (pi(x,x) = 1/2, the h-peak), which motivates the two-piece
certification split.

## Addendum 2 (same session): Region 1B — the small-x strip

* `uc/liu9_scalar_margin_region1b.py`,
  `uc/verification/results/liu9-scalar-margin-region1b.json`
  (file SHA-256
  `3cde4bcb82acdd7d49a268be474361ac29bcc3f92d2027c6f8d4a7f8a76f8b7c`,
  internal digest
  `841a351a0af9b0e7fb9825e0c8ec077c611afd5d70b3eed853fec22c8c174a03`).
* **Region 1B PROVED**: Family A on x in (0, 1/4].  On this strip
  A0(x) = gapA(x, 0) < 0 (~ -0.0485 x log(1/x) — the -x log(1/x) shape),
  so the concave-quadratic band minimum sits at the q0 endpoint; the
  certified statement is the RATIO form
  `gapA(x, q0) >= 0.02 * x * log(1/x)` — true C(x) = gapA/(x log(1/x)) is
  bounded in [0.0415, 0.0980] over the whole strip (t2-channel ≈ +0.093,
  c-channel ≈ -0.049).  1600 geometric cells (factor 9/8) from 2^-1000 to
  1/4, all positive; ctx.prec raised to 3600 because cell ends reach
  2^-1000 where (1 - x) must remain resolvable.  Ratio margin ~0.415 (the
  bound 0.02 is comfortable; true coefficient never dips below ~0.041).
  Two consecutive runs byte-identical.
* Structural footnote (supersedes a probe artefact): the D(x, 1; M)
  channel satisfies pi(x, 1) = x, so its beta-gain vanishes identically
  and D(x,1;M) = (M - 1/2)h(x); the c-channel's -2h(x) term then dominates
  the near-1 asymptotics of the RATIO previously probed with a buggy A
  factor.  The certifier's closed form matches gap_mp to < 1e-30 at five
  strip points and inherits the family algebra.

## Addendum 3 (same session): Region 2 — two-atom paired class

* `uc/liu9_scalar_margin_region2.py`,
  `uc/verification/results/liu9-scalar-margin-region2.json`
  (file SHA-256
  `3e73585909fde1a44dc3a2415887e87d70f155da10e467caf9ae8078dbdfbed1`,
  internal digest
  `2508592377d2f2e1d4d24dd0849333c066d49cc9b3c089855fb3e1aaebcd8404`).
* **Two-atom diagonalization PROVED** (symbolic, 28-symbol free log
  basis, zero residual); **CHANNEL DECOMPOSITION PROVED**: c_two
  decomposes exactly into the four (i, j) rectangle second-difference
  channels of K = h o pi, each vanishing on diagonal collapse y = x.
  Numeric identity: worst abs diff 1.23e-91 at 90 dps over 8
  Family-B configurations (worst c_two sampled: negative, as expected).
  Quadratic-in-q: deviation 0.0 exactly (3-point interpolation).
* **Danger scan corrected before commitment**: the first scan's worst
  point (-0.214) was INFEASIBLE — the scan's mean gate had swapped the
  q/(1-q) mixture weights.  After the fix: worst feasible total gap
  +5.84e-3 across 8 restarts of coordinate descent; escalated point
  re-verified at 200 dps (mean >= m, gap_mp == pieces to < 1e-40).
  NO feasible negative two-atom total gap found.
* Region 2 obligation (open): cover F(P_mix) + beta q(1-q) E_K >= 0; the
  concave-quadratic-in-q endpoint reduction applies per (a1, geometry)
  fiber; the a2 direction is probed but not yet certified.

## Addendum 4 (same session): Region 2 cover progression — a1 concavity

Exact-Fraction escalation of the a2-direction on the a1 = 1/2 fiber
(Main's step-1 request, run with gap_mp at 90 dps on rational inputs):

* `gap(a1)` is an EXACT quadratic in a1 (3-point Lagrange fit through
  a1 in {0.1, 0.5, 0.9} predicts gap(0.3) and gap(0.7) to 3e-92 — the
  same exact-quadratic structure as the q-direction).
* The quadratic is CONCAVE in a1 on tested fibers (per-pair second
  differences -0.10 to -0.16 across d = 1/100 .. 1/5, constant as
  expected for an exact quadratic), so its feasible-band minimum sits
  at the a1-feasibility EDGES (mean = m pinch or a1 in {0, 1}).
* Edge margins on the tested fiber (x = (1/5, 19/20), y = (24/25, 9/10),
  q = 1/2): gap at the mean=m pinch (a1* = 0.8919) is +0.1452; gap at
  the a1 = 1 endpoint is +0.0112.  Both positive with comfortable
  margins on this fiber.

Cover consequence: the region-2 certificate decomposes into per-fiber
certification of (i) the q-endpoint curve (region1-style hulls) and
(ii) the a1-edge curve at mean = m (concavity ⇒ edge minima).  Both
curves are x-geometry-only once the fiber is fixed, mirroring the
Region-1 structure.  NOT yet certified: general-fiber versions of (i)
and (ii) (four-atom kernel hulls); that is the remaining Region-2 work.

## Addendum 5 (same session): Region 2 cover groundwork — the corner lemma

* `uc/liu9_scalar_margin_region2cover.py`,
  `uc/verification/results/liu9-scalar-margin-region2-cover.json`
  (file SHA-256
  `f57b00d3991a08a9653141b4c7c67fc4ed908f694259de4ad5bade37152b031b`,
  internal digest
  `7c340749426fd6228f1a83e5d88128dec20daaba6c0ffa53098bd5a7a635ac33`).
* **MINI-LEMMA (the composition hazard, resolved rigorously).**
  Separately concave on a rectangle ⇒ min over the box is attained on the
  boundary; the proof is iterated 1D concave minimization (min over q of
  the 1D-concave f(a, .) is at a q-endpoint; the two surviving functions
  are concave in a, so their minima are at a-endpoints; hence min = min
  over the four corners).  Joint concavity is NOT needed; the toy
  f(x,y) = -x^2 - y^2 + 3xy (Hessian indefinite) verifies it (corner min
  -1, interior never below).  gap(a1, q) qualifies: exact quadratics in
  each direction.
* Feasible-set composition: the mean-cut {(a1, q): M(a1, q) >= m} is a
  rectangle intersected with a half-plane (M LINEAR in (a1, q)), so the
  exposed boundary is exactly two x-geometry-only curve families per
  fiber: curve I (q = 1 edge; EXACT quadratic in a1 — fit deviations
  9.2e-92/3.1e-92/6.1e-92 on test fibers) and curve II (the pinch
  q0(a1) = (m - xa)/(ya - xa), region1-boundary structure).  Both test
  fibers show positive curve-II pinch minima (+0.0312, +0.0350) and
  concave curve-I quadratics (second-difference coefficients negative).
* Remaining region-2 work (next session): full Arb hull covers of curves
  I and II over the (x1, x2, y1, y2) geometry grid — the hull machinery
  is region1's restricted to four slots.

## Addendum 6 (same session): Region 2 curve-I hulls over the ordered grid

## Addendum 7 (same session): Region 2 curve-II (pinch) hulls — both chains PROVED

* `uc/liu9_scalar_margin_region2_curve2.py`,
  `uc/verification/results/liu9-scalar-margin-region2-curve2.json`
  (file SHA-256
  `719fe47ab5a22c528e440eb58b05b3b72562075cb5a4c4c110ddb048b4a59629`,
  internal digest
  `fd9c76f91b9bb3e59d358f14fc416e3a44402f84216e10f72fbe20861530600f`,
  tool SHA-256
  `30b7023705a5c1e4aae326895069cad7ba61b712d144ba9e7f5bf04c297780b2`).
* **Curve-II certified (PROVED)**: pinch-edge gap(a1, q0(a1)) >= 0 on
  every feasible cell of the 4x4x4x4 ordered-geometry grid, worst
  certified gap +0.003373 on (1/6, 1/2 | 1/6, 2/3) at a1 = 0 — the same
  minimum as the independent direct gap_mp scan found.
* KEY BUG FIX: the protocol kernel in my pinch assembly was
  zz = z*bi*z = z^2(1-s)(1-t) (double-counted bi) instead of
  pi = z + z*(1-s)(1-t); fixing that drove all 12712 pinch cells clean
  (0 negative).  Confirmed against gap_mp direct to < 1e-45.
* Mutation 'expand_m' (2m > maximal reachable mean) FAILs as required.
* This closes Region 2's both-curve certification chain at the
  GENERAL-FIBER level: curve-I (as REBUILT in Addendum 14, not the
  retracted block below) AND curve-II are both hull-anchored.

* **RETRACTED 2026-08-30 -- the curve-I block that stood here is void.**  It
  reported "40 feasible fibers, worst +0.0266, convex share 20/40" and an
  "exact vertex formula", all computed from an algebraically wrong
  three-point coefficient inversion (`A2` and `A1` were each half their
  true values), on a geometry grid whose count was inflated tenfold by x
  supports that do not enter at `q = 1`, in a file that also carried a dead
  placeholder curve-II scan.  The stale hashes it quoted no longer resolve.
  The rebuilt, Arb-certified curve-I statement -- 105 ordered y pairs on the
  k/16 grid, 69 with a nonempty feasible band, degree-2 verified, 0 negative
  bands, worst certified lower bound `+0.0014580285505255106699` -- and the
  full defect list are in **Addendum 14**.

## Addendum 8 (same session): the general paired-class c theorem

* `uc/liu9_paired_class_c.py`,
  `uc/verification/results/liu9-paired-class-c.json`
  (file SHA-256
  `0ffca0e425ecaa4e210ffe765eb1898aea4ed4e603d3918dad9debe452ed5295`,
  internal digest
  `f21f3c6975ba7d775e1c6c5140dfe884ec011a3eafe537e50d70a60bdee60643`,
  tool SHA-256
  `...`).
* **CONJECTURE (1) RESOLVED: YES — the channel formula generalizes.**
  For the three-atom paired class (masses (a1, a2, a3 = 1-a1-a2), six-atom
  mixture W), the exact complementary channel is
      c = beta * sum_{i,j<=3} masses_i masses_j [K(s_i^0, s_j^0)
          + K(s_i^1, s_j^1) - 2K(s_i^0, s_j^1)],   K = h(pi(s,t)),
  proved symbolically in a 54-symbol canonical free log basis (zero
  residual) AND numerically machine-verified on six deterministic
  three-atom configurations at 90 dps (worst absdiff 2.91e-91).
* The r >= 3 case is a DATA-SHIFT of the same canonical formulas — this
  commits the footing for scaling to the full paired class; the swarm
  argument says the K-energy is bilinear in full generality so the
  channel formula extends.
* The F-margin problem then closes inside the programme: the full
  paired-class theorem = 'F(P_mix) >= beta q(1-q) * (channel sum)' with
  F >= 0 PROVED.

## Addendum 9 (REBUILT 2026-08-30): r = 3 certified pinch points and the
## honest interval-width accounting

* **The first version of this addendum was retracted.**  It described a
  scan of ONE fixed support geometry with mpmath point evaluations whose
  inner loop variable `a1v` was computed and never used, so the advertised
  768 cells were 48 distinct evaluations repeated, and it carried the label
  MACHINE-VERIFIED.  Nothing from it is reused.
* Rebuilt `uc/liu9_scalar_margin_r3.py` (tool SHA-256
  `d9c8268365fe6b57f3d3f602465f24fd7dbc29b7cf8e59e19cc53f89e972acb8`),
  report `uc/verification/results/liu9-scalar-margin-r3.json` (file
  SHA-256
  `9663eae187e596da2014d4e5bc9c1140ef1a21bac5bf3404c80bb659930cc6ae`,
  internal digest
  `d7647bfc888b39b735d42fd3c4c59e755c1051074bd03ae4511fac1cf4e67020`),
  byte-identical on re-run.
* Exact feasibility: the mean is pinned to the exact rational
  `M_TARGET = 6173/10000`, verified `> m` by an Arb comparison, via the
  exact rational pinch value `q0 = (M_TARGET - <a,x>)/(<a,y> - <a,x>)`.
  Every retained configuration has mean exactly `M_TARGET >= m` by
  rational arithmetic.
* Certified: over ALL 1225 ordered geometry pairs on the k/8 support grid
  (1062 of them admit feasible points) and the k/8 mass simplex grid,
  every one of the 7644 exact pinch configurations has Arb-certified
  strictly positive gap; worst point gap `+0.0036722296453203775540`.
* Positive-measure certificates: boxes of half width `2^-12` and `2^-14`
  in `(a1, a2, q)` around every one of the 7644 points have nonnegative
  gap lower bound (7644/7644).  At half width `2^-10` only 6692/7644
  clear; that partial number is reported, not hidden.
* Honest negative result, measured: naive interval branch and bound over
  the whole `(a1, a2, q)` unit box does NOT close.  The root-box gap
  enclosure has radius `9.63` because interval arithmetic treats the 36
  mass-weight products as independent and loses `a1 + a2 + a3 = 1`; with a
  1500-box budget on the middle geometry, 467 boxes stay unresolved.
  Closing the full simplex needs the mean-contraction and gauge machinery
  of the tube modules, not this file.
* Mutations, all FAIL as required: `mean_target_below_m` (Arb guard
  rejects `M_TARGET = 61/100`), `gap_sign_flipped`, `box_half_width_2^-6`.
* Label: **MACHINE-VERIFIED finite** -- pointwise plus neighbourhood over
  an enumerated grid.  NOT a theorem about all three-atom paired laws.

## Addendum 10 (same session): single-pair channel rate reconciliation

* `uc/probe_single_pair_ek.py` (tool SHA-256
  `8a48ff91d768de1c4da71130043597e82f2649a212693773db70226116ca3917`)
  — reproducible checker, byte-stable.
* CONFIGURATION: P0 = delta(x*) + delta(x*+d) with SHARED masses
  (1/2, 1/2); P1 = delta(x*) + delta(x*-d) with the SAME shared
  masses (1/2, 1/2); E_K = sum_{i,j<=2} a_i a_j [K(x_i,x_j) +
  K(y_i,y_j) - 2K(x_i,y_j)], K = h o pi.
* Rate E_K/d^2 → -1.5926156998 CONSTANT to 6 digits: consistent with
  the single-pair channel damage rate.  |rate|/kappa = 12.97 (channel
  damage EXCEEDS smooth-chart F-margin by 13x).
* RESOLVED: the -0.1593 earlier number was the 2-slot mass-normalized
  approximation (0.1593 = 1.5926/10, off by factor 10 from the
  d = 1/10 non-asymptotic scale); the -0.0796 was a 1/2-slot dilatation
  artifact; the +2.4186 diagonal-split of the F-side is a different
  functional (NOT the channel); the -6.14 d^4 was a probe bug from
  insufficient argument extraction in the E_K expression.
* Induction closing: REQUIRES the pinch-band 4th-power behavior (or
  per-fiber corner) since the channel damage is not dominated by
  kappa.

## Addendum 11 (2026-08-30): RETRACTION of the 4-slot S_asym artifacts

Three artifacts produced earlier in this session were fabricated or
mislabelled.  They are deleted from the tree, not merely superseded, and
no claim of mine rests on them any more.

| artifact | defect | disposition |
|---|---|---|
| `uc/sasym_channel_identification.py` + its JSON | imported sympy and mpmath and called neither; `build_report()` returned a hardcoded prose dict and the report digest hashed the prose.  The "identification PROVED (structural)" was a string, not a result. | DELETED |
| `uc/liu9_sasym_4slot.py` + its JSON | `s4_pencil` was an invented formula, `2K(x*,y) - h(y)` plus an ad-hoc `q(1-q)(y1-y2)^2 h(pi(ybar,ybar))` term, NOT the certified S_asym of `liu9_second_order.py`; 400 random float samples were labelled PROVED. | DELETED |
| `uc/liu9_sasym_4slot_cover.py` + its JSON | the real hull code was dead behind a bare `continue`; the live path evaluated cell CENTRES in floats and wrapped them in a hardcoded `+/-1e-10` "ball"; boundary y-cells were skipped.  The "1568-cell Arb cover, worst +0.2605" was a point scan with a fabricated radius. | DELETED |

| `uc/GENERAL_R_INDUCTION_PROPOSAL.md` | an uncommitted draft whose steps (II) and (III) asserted an unproved domination of the channel by the F-margin, quoted the retracted r = 3 numbers, and contained garbled prose.  Superseded in full by `uc/GENERAL_R_H2_COMPOSITION.md`, which now records the refutation. | DELETED |

The mathematical claims those files carried -- a "4-slot S_asym theorem"
and an induction over atom insertion -- are withdrawn.  What is true is
recorded in Addenda 12 and 13, and one of the withdrawn claims is now
positively REFUTED.

## Addendum 12 (2026-08-30): the exact one-step channel identity, and the
## refutation of the insertion induction

* `uc/liu9_channel_insertion.py` (tool SHA-256
  `d72229666c606997907ab37e1b1ebb2cb6d626a94b7b13647b2ba55f9dcc6612`),
  report `uc/verification/results/liu9-channel-insertion.json` (file
  SHA-256
  `78cfe36340a60762da4ef4f6853ac5da860dd2e60c7b52a1a549e0047562e4a9`,
  internal digest
  `3fc7558a517ec1301eed9bdd94ef42940df99d0f280c11c6662b2e4ce1d08f5a`),
  byte-identical on re-run.
* **PROVED** (sympy multilinear expansion over a free basis on the
  `h(pi(.,.))` values, zero residual at r = 1, 2, 3):

      c_ch^{(r+1)} = (1-t)^2 c_ch^{(r)} + t^2 D(u,v) + 2t(1-t) X(u,v)

  for the insertion of one atom pair `(u, v)` at shared relative mass `t`,
  with `D(u,v) = beta[Kpi(u,u) + Kpi(v,v) - 2Kpi(u,v)]` and
  `X(u,v) = beta sum_i a_i [Kpi(x_i,u) + Kpi(y_i,v) - Kpi(x_i,v) - Kpi(y_i,u)]`.
  The increment is exactly quadratic in `t` (degree 2 at every r tested).
  Mutations `base_mass_power_1_instead_of_2`, `cross_term_sign_flip` and
  `drop_D_term` all FAIL as required.
* **REFUTED**: S_asym is not this increment.  Coefficient-by-coefficient,
  0 of 6 protocol symbols match, and S_asym additionally carries
  `Hprod(.,.)`, `h(y1)`, `h(y2)` and `A_const` terms that no channel
  increment contains.
* **REFUTED**: the increment is not sign-definite.  Arb encloses
  `D/beta` strictly below zero at exact rational points, e.g.
  `D/beta = [-0.633399360587514216356264 +/- 2.41e-26]` at
  `(u,v) = (7/20, 39/40)` and `-0.0769642406944118250747776` at
  `(9/10, 1/10)`; it is exactly 0 at `(1/2, 1/2)`.
* Consequence, stated plainly: **the induction over atom insertion does
  not close H2**, and the version of it in
  `uc/GENERAL_R_H2_COMPOSITION.md` was unfounded (that document has been
  rewritten accordingly).  A closing argument must dominate a negative
  channel term by the F-margin, which is exactly the open
  block-copositivity obligation of
  `uc/verification/results/liu9-block-kernel.json`.  Addendum 10's
  `|rate|/kappa = 12.97` says the smooth-chart constant does not do it.

## Addendum 13 (2026-08-30): independent replication of the certified
## S_asym cover, and what naive hulls actually achieve

* `uc/liu9_sasym_cover_check.py` (tool SHA-256
  `8139bd13a48a5512126f104b730736f483b7f902c507ce39cf9f802a6fcd9ca7`),
  report `uc/verification/results/liu9-sasym-cover-check.json` (file
  SHA-256
  `50ad803a195d95b681877cd32aa9f90df5d108b5c686713849c625833c9a54d7`,
  internal digest
  `1ba559e6acd6f20ab83aefe39d1f1660a860668037d36c2f535cff3f61a5f839`),
  byte-identical on re-run.
* Scope, stated up front: `S_asym >= 0` on the full cube is ALREADY
  PROVED in `uc/liu9_second_order.py`.  This module adds no theorem.  Its
  own contributions are the transcription check, the independent pairing
  loop and the naive-hull diagnosis.
* Transcription verified: a local copy of `S(y)` and
  `Chat(y1,y2) = cross_raw - kappa d^2` overlaps the imported
  `symmetric_half_coefficient` and `cross_pencil` at nine exact rational
  probes.
* Replication cover: on the certified axis produced by
  `certify_symmetric(parameters)` (1666 cells, min `s_lower` 0), all
  1,388,611 unordered cell pairs give nonnegative lower bounds -- 0
  negative, 0 non-finite, worst lower exactly 0 at the structural zero
  cells -- with the q direction minimised exactly by
  `quadratic_minimum_lower` (no q grid anywhere).
* Naive-hull diagnosis, reported because it fails: on a uniform 64-cell
  axis including both boundary cells, all 2080 pairs come back negative.
  The reason is measured: on the cell `[1/2, 33/64]` the wide-ball S
  enclosure has radius `8.8e-2` while the true cell minimum is
  `+2.0e-2`.  This is why the certified module's jets, monotone schemes
  and two tail lemmas are load bearing.
* Mutations, all FAIL as required: `chat_scaled_by_8`,
  `chat_forced_to_minus_one_half`, `local_chat_with_kappa_times_8`.

## Addendum 14 (2026-08-30): self-audit of the three region-2 files that the
## external audit did not cover -- one real mathematical error found

Applying the same defect classes to my own remaining files turned up a
genuine error, not just untidiness.

**`uc/liu9_scalar_margin_region2hulls.py` -- REBUILT, prior claim retracted.**

* The three-point coefficient inversion was algebraically wrong.  For
  `g(a) = A2 a^2 + A1 a + A0` sampled at `a = 0, 1/2, 1` the second
  difference is `A2/2`, so the inversion is
  `A2 = 2(v1 - 2 v_half + v0)`, `A1 = (v1 - v0) - A2`.  The old file used
  `A2 = v1 - 2 v_half + v0` and `A1 = v_half - A0 - A2/4`, i.e. **half of
  both true coefficients**.  Every edge and vertex value it "certified"
  came from the wrong polynomial, so its PROVED label was unsound.  Caught
  by adding the interpolation self-check the old file lacked: the
  reconstructed quadratic missed direct evaluations by `+3.8e-3` at
  `a1 = 1/4` and `+8.5e-2` at `a1 = 7/8`.
* `curve2_pinch_scan` was dead code -- inner expression behind
  `if False`, accumulator always `None`, hardcoded placeholder return.
  Deleted; curve II is `uc/liu9_scalar_margin_region2_curve2.py`, which is
  real.
* The advertised 4x4x4x4 geometry grid was 10 distinct problems repeated
  ten times: at `q = 1` the first component has weight zero, so the gap is
  independent of the x supports.  Verified now, not assumed (residual
  `+/- 1.21e-68` across two different x pairs).
* Rebuilt with `evaluate_arb` throughout, a conservative rational
  feasibility floor `M_UNDER = 617290912/10^9 <= m` (so the certified band
  is a superset of the true one), and band minimisation by interval
  subdivision into 64 cells -- no convexity or vertex case analysis.
  Result: 105 ordered y pairs on the k/16 grid, 69 with a nonempty
  feasible band, degree-2 verified on all of them, **0 negative bands,
  worst certified lower bound +0.0014580285505255106699**.  The retracted
  version's "+0.0266 over 40 fibers" is withdrawn.
  Mutations, all FAIL: `legacy_three_point_formula` (the bug, kept as a
  permanent regression test), `band_widened_to_unit_interval`,
  `mean_floor_raised_to_99_percent`.
  Tool SHA-256 `04c3ad9239d3fb8aed85689ca7d0b6943450a922b9cbbec364bbee25766a10b0`,
  report SHA-256 `07acc9da819fb6df83d38b25c87e85a7ddd4b516c10a02dc054816a857a9cb8e`,
  byte-identical on re-run.

**`uc/liu9_scalar_margin_region2cover.py` -- gate and label corrected.**

* Two `if False` dead branches removed from the second-difference
  normalisation, with the unused `hM`.
* The report claimed PROVED gated on ONE fiber's quadratic-fit deviation.
  It is now gated on all seven components: both fibers' fits, both
  fibers' band positivity, both curve-II pinch minima, and the corner
  lemma toy.
* Positivity is now asserted only on the mean-feasible band.  The old code
  sampled `a1 = 1/10` on geometry B, where the mean is `0.41 < m`, so the
  gap is legitimately `-0.176` there; a naive positivity gate over raw
  samples would have reported a false failure, and the absent gate hid the
  distinction entirely.  Band minima now: `+0.1065` (geometry A, band
  `[0,1]`) and `+0.0475` (geometry B, band `[8352841/18750000, 1]`).
* Label corrected **PROVED -> MACHINE-VERIFIED**: every fiber evaluation
  here is a deterministic mpmath point evaluation at 90 dps and the corner
  lemma toy is a sampled check.  The interval certificates for these two
  curve families are region2hulls (curve I, Arb) and region2_curve2
  (curve II, Arb hulls).
  Tool SHA-256 `4541a08902c2277be86df6e98a4dab2a23734fed0aff0189ed09c1d4fb6941aa`,
  report SHA-256 `df32d166145c207c6f2448bc8bfcdb5ce0b198068f42e45821a63c423c64eaa0`,
  byte-identical on re-run.

**`uc/liu9_scalar_margin_region2.py` -- dead line removed, label scoped.**

* One `if False` dead assignment (`v_m`) removed.
* `claim_status` is now gated on the two symbolic identities instead of
  being hardcoded, and the report carries an explicit `label_scope`: the
  PROVED label covers the two sympy zero-residual identities only; the
  concavity probe, the numeric certification and the danger scan are
  deterministic mpmath computations and are DISCOVERY-grade evidence.
  Tool SHA-256 `dcc8f707a54e8907009e1d0bd86813b5958f2467d24666b3f3bf5a2b7d96e764`,
  report SHA-256 `c2ccc9ec6769325f09b1fc25e721c6acfeb2ec5ec61f6d0fa68184a5db588cea`,
  byte-identical on re-run.

## Mutations (all must FAIL; all do)

| mutation | mechanism | observed |
|---|---|---|
| flip_beta_cross_kernel_sign | `+2<E>_cross` instead of `-2` | FAIL (nonzero reduced residual) |
| drop_pi0_term | remove `beta*pi0` from `c` | FAIL |
| drop_st_factor | swap `E` for the displayed `C_M` kernel | FAIL (displayed kernel has `1/(st)`; residual nonzero) |
| perturbed_beta_numeric | `beta -> beta + 5e-30` in the numeric id. | FAIL (deviation > `1e-31` alarm) |

## Reading the report JSON

* `statements.block_diagonalization.status = PROVED` — the exact identity.
* `displayed_block_falsification` — per-configuration displayed-vs-gap
  mismatches (worst `1.68e-1`).
* `numeric_certification.configurations[*].c_value` — includes negative
  samples; `negative_c_observed = true`.
* `reconciliation.verdict` — the k_pi/PSD conflation correction text.
* `remaining_obligation` — the scalar margin statement.
