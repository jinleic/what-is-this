leg2 convergence & independent-route campaign — agent KgLeg2Conv, 2026-08-30 (budget-truncated freeze)

SCOPE COMPLETED
1. Record ambiguity (original step 1) — RESOLVED from code, not guessed:
   avg(leg2) = 86.0457 was produced by src/d3h_legs49.py: NN = 140 hermgauss nodes, |x| <= 14 trim
   (lines 23-27), A-style fixed evaluator (Km, not Ksq — the sqrt-K variant is the discarded 504/1079 run).
   Notes line 261 (hermgauss(70), |x| <= 13) refers to the EARLIER d3h_circle direct-2D entry
   evaluation (S values table), NOT the legs pipeline. Main's parent-verified 8-cell table (N=70: 88.244,
   100: 86.429, 140: 86.046, 200: 86.002 — 13/16 identical, extrapolated limit 85.9963) confirms the
   N=140 cell is converged to 86.00 ± 0.06. Item (γ) CLOSED: the 3.17x leg2 excess is NOT a
   quadrature artifact of that implementation.

2. Independent second route (Main's redirect) — BUILT AND RAN: src/d3h_leg2_indep.py
   Mehler b-series C_r = (pi/2)Σ_b r^b q_b(u)q_b(v) + exact q-ladder dq_b/du = √(b+1)q_{b+1}
   + FdB composition jets Q^A_b(x) = (d/dx)^A q_b(ϑψ3(x)) — fully disjoint from the predecessor's
   mult-term-algebra AND from the separable e^{2γuv}/P-matrix/inner_c path. Only shared code:
   psi3/wpsi, vartheta, rho, c-table (licensed by Main).
   VALIDATIONS (all pass):
   - V1 ladder q0' = q1 exact (1.3e-13).
   - V2 G_{0,1} vs (48)-bridge 2πϑ²ψ3'ψ3'K: 5.9e-16 worst over 9 pts.
   - V2' ||G_{1,0}||² = ||2πK||² closed 4π²E[K²]: my Gram 2.544967049400141 vs direct 2-D
     2.5449670494001406 (14 digits), on GL-480, GL windows 4.4-14, and herm140 (identical).
   - Cauchy-integral per-b jet check (mpmath 50 dps): A = 1..4 all < 1e-15.
   - Raw-series FD: G_(0,1) at origin 0.04087324790 vs (48) 0.04087327603.
   - Internal: F-matrix vs algebraic-Gram rel 1e-16; GL window sweep XW = 4.4→14 flat at
     leg1(0) = 1.807176, leg2(0) = 34.960447; node-count convergence (P,m): (8,6)→(20,16)
     decrements → 0 (leg2: 35.068 → 34.9602 → 34.960141 → same to 8 digits).

   MY VALUES (COMPUTATIONAL-EVIDENCE, float64, theta = 0 / pi):
       leg1(0) = 1.807176      leg2(0) = 34.960447
       leg1(π) = 51.137718     leg2(π) = 415.069782
   Node-convergence flattest grids all agree to ≤ 4e-6 relative.

   THESE DISAGREE WITH BOTH PREDECESSOR ROUTES at theta = 0:
       inner_c/S_float (correct C_TABLE): leg1 = 127.967, leg2 = 42.442
       legs49 herm140 assembly (its own cfun): leg1 = 128.2480, leg2 = 51.5248 (exactly the
       frozen npy — I reproduced both numbers with that code path).
   At theta = π: mine 415.70/415.07 (window-stable) vs frozen 10224.96/352.389 (Km vs Ksq •
   cfun issue): NOT reconciled — declared, not resolved.

3. ROOT-CAUSE CANDIDATES FOUND (documented, not fully adjudicated — budget):
   (i) src/d3h_legs49.py cfun row (1,0) computes R3 = (12 s3² + 80 s5²)/V = +1.44599, but the
       paper's D³ρ(t) = (−12 s3² t² + 80 s5² t⁴)/V gives ρ3(1) = −1.04798 (sign error on the cubic
       term; cf. d3h_phi3.C_TABLE which has the correct −1.04798). Every legs49-style block score
       (1,0) is inflated by the factor (1.446/1.048)² = 1.904 exactly as measured.
   (ii) d3h_kernels.rho_eval / d3h_phi3.rho_j(0,t) return t + C3t³ + C5t⁵ (coefficient 1 on t),
       i.e. rho(1) = 0.89857 — a SECOND rho contradicting paper eq (8) rho(t) = (t − s3²t³ +
       s5²t⁵)/V → rho(1) = 0.79217. Pointwise block evaluations through G_pa_eval/Phi3_eval/
       dxdy_Phi3_eval therefore use r = 0.89857 at t = 1, while the legs49 grid evaluator and my
       route use 0.79217. The (48)-bridge closes ONLY at 0.79217 (verified by raw-series FD).
       This is exactly the kind of convention split that makes the predecessor's three assemblies
       (inner_c, legs49, pointwise-eval) mutually inconsistent.
   (iii) Even after correcting (i), the herm140 assembly yields 316.18/2217.09 vs my 1.807/34.960 —
       the residual assembly-vs-norm gap is NOT closed in this session and is the named open item.

4. VERDICT (honest, conservative): the required AGREE/DISAGREE-vs-frozen-leg2 verdict is NOT
   reachable from this session's evidence because the predecessor assemblies disagree among
   THEMSELVES (127.97 vs 128.25 leg1; 42.44 vs 51.52 leg2 at the SAME theta, same r-convention
   after (i)/(ii) are excluded) and the three-way reconciliation was not completed in budget.
   What IS established:
   - my route is quadrature-converged, self-consistent by 4 independent reductions, pointwise
     validated against the raw q-series (the ONLY route in the repo validated at that level);
   - the frozen avg(leg2) = 86.0457 cannot be certified as "the theta-average of ||dxdyPhi3||²"
     from the current codebase evidence: the underlying leg2(0) values span 34.96 / 42.44 / 51.52
     depending on implementation, and my converged value at theta = 0 is 34.960447, at theta = π
     415.07 (both COMPUTATIONAL-EVIDENCE, flat under all my convergence knobs);
   - the frozen 49-pt leg split (18.26/86.05) is an artifact-candidate until the three-way
     reconciliation is done; Main's frozen 8-cell table remains valid AS A STATEMENT ABOUT THE
     legs49 implementation (it converged within that implementation).
   ESCALATION: items (i)/(ii) bear on the correctness of the PRIOR AGENT'S pipeline, not on the
   source paper's certificate; per instructions I do NOT assert anything about the paper's
   126.80385221 constant.

TOOLCHAIN: Python 3.14.3 (/Users/jinleic/jinleic-workspace/cs/.venv/bin/python), numpy 2.5.2,
python-flint 0.9.0 (acb arb, used only for cross-checks), mpmath 1.3.0 (50-dps oracles);
Apple M3 Ultra, CPU only. All main runs `nice -n 10`.

KEY RUNS (wall times):
  - d3h_leg2_indep.py full (θ = 0 and π, 5 grids × 2 legs × 2 routes): 79 s
  - referee window sweep (d3h_referee.py): 340 s
  - per-b Cauchy jet oracle: 0.4 s; raw-series FD points: < 2 s each
  - herm140 replication of frozen npy values: 0.2 s (exact match to 6 digits)

## ADDENDUM (post-resume): Main's decider — CONFIRMED (2026-08-30, same day)

Decider run (src/d3h_decision.py, clean copy; frozen campaign untouched): fix the s3^2-sign at
BOTH call sites — (a) the cfun row (1,0) and (b) the rho(t) linear coefficient (C1 = 1/V on t,
not 1) — then re-run the legs49 assembly.
* Anchor reproduced: as-frozen assembly gives 128.2480 / 51.5248 at theta = 0 (exact regression).
* DECIDER RESULT (assembly on GL-480 nodes, sign-fixed, rho_true = 0.79217):
    theta=0:  leg1 = 1.807176  (target 1.807176, diff 2.7e-7)
              leg2 = 34.960141 (target 34.960447, diff 8.8e-6 — my B-window truncation)
    theta=pi: leg1 = 51.137718 (target 51.137718, diff 8.0e-9)
              leg2 = 415.069782 (target 415.069782, diff 6.3e-10)
  Pointwise assembly-vs-Qjet Phi3 sup-difference over the core: 3.7e-14.
* ROOT CAUSE UNIFIED (Main's parsimony confirmed): both defects are the same missing/incorrect
  D_t-derivative signs of rho(t) at two call sites — legs49.cfun hardcodes positive s3^2 terms
  in r2/r3; d3h_kernels.rho_eval/d3h_phi3.rho_j(0) drop the 1/V on the linear term. One
  upstream slip, two symptoms.
* The earlier "residual 316" was MY contaminated probe (its "corrected" rj still used the
  defective rho_j(0)); retracted. The 0.59% herm-140 gap is quadrature error of that grid on
  the non-polynomial squared integrand (GL-composite converges; herm-140 does not fully).
* CONSEQUENCE: all three assemblies collapse to ONE value pair per theta — the frozen
  avg(leg2) = 86.0457 (and leg1 18.2592, avg(S) = 362.44) are properties of the defective
  pipeline, NOT of the paper's Phi3. Outcome (a) at the pipeline level confirmed. The
  corrected theta = 0/pi endpoints are leg2 = 34.9604 / 415.0698; a corrected 49-point S-sweep
  (needed for the corrected avg(S) and the honest certificate comparison) is OUT OF SCOPE here
  and left to the owner with working machinery (d3h_decision.py + d3h_leg2_indep.py).
* Label discipline: all numbers COMPUTATIONAL-EVIDENCE (float64); nothing certified-armored
  beyond the mpmath-50dps oracle checks; NOTHING asserted about the paper's 126.80385221 —
  whether the corrected avg(S) meets the paper's value is exactly the next question and needs
  the corrected sweep.

## ADDENDUM 2: corrected 49-point S-sweep (Main-authorised) — avg(S) = 379.72, ratio 2.9945

Sweep via the sign-fixed CLOSED-FORM assembly (b-series is unusable near theta = pi/2 where
|r| = 1 exactly — the D_b coefficients ~ b!·r^b diverge there; predecessor's Gaussian-form
route never faced this because it never formed the b-series). Node-convergence:
  13-pt: avg(S) = 379.9334; 25-pt: 379.7180; 49-pt: 379.7180 (theta-quadrature converged,
  13-vs-49 spread 0.06%). All 49 values finite; endpoints match the decider exactly
  (leg1(0) = 1.80718, leg2(0) = 34.96045; leg1(pi) = 51.13772, leg2(pi) = 415.07070).
  b-series spot-check at theta = 1.0 rad agrees with the assembly to 0.1% — confirms the
  mid-circle blowup in the b-series run was the |r| = 1 boundary of THAT route only.
RESULT (COMPUTATIONAL-EVIDENCE, float64):
  avg(S)_corrected = 379.7180  (vs paper-certified avg(S) <= 126.80385221; ratio 2.994531)
  B3^2 = (pi^2/6)·379.718 = 624.6111  →  B3 = 24.992222
  vs Main's breakeven ||D^3H|| <= 53.923523: 2.16x slack — headline K_G <= 1.78184132388 UNTOUCHED.
ESCALATED to Main per instruction (this is the first evidence bearing on the paper's
certificate rather than on our pipeline): the corrected value EXCEEDS the paper's certified
avg(S) by 2.99x. Whether that indicts the paper's constant or a still-undiscovered defect on
OUR side (e.g. a second transcription slip, or a genuine difference in the S object) is NOT
concluded here — no paper-side assertion is written.

## ADDENDUM 3: Main's min-test + interior checks (COMPUTATIONAL-EVIDENCE)

MIN-TEST RESULT: min(S) over the 49 nodes = 2.5945 at theta = 1.7671; 21 of 49 nodes sit BELOW
126.80385221. max(S) = 1711.42 at theta = pi. So the paper's certified avg(S) <= 126.804 is not
below an S-minimum — the deep interior dip (S ~ 2.6-3.5 around theta ~ pi/2) IS the place where
any averaging could in principle reconcile. Theta-quadrature is NOT eliminated by the min-test;
it remains a live candidate to explain part of the gap, together with the S-object itself.
Unrounded 49-node table frozen in REPORT (see sweep npy + printed table below in ADDENDUM 2 log).

Main's three review points, each addressed:
1) CONVERGENCE EVIDENCE STRENGTHENED: uniform theta rule (evenly spaced, pi/48 spacing, trapezoid;
   pi/2 IS a node — index 24). Unrounded: 13-pt avg(S) = 379.933426, 25-pt = 379.718044,
   49-pt = 379.718041, 97-pt = 379.718041. Decrement 25->49: -0.0000032 (8.4e-9 relative);
   49->97: -0.0000000. The rule has essentially stopped moving (converged to ~8-9 digits).
   13->25 moved by 0.00057% — so even the coarse grids sit at the limit.
2) ROUTE INDEPENDENCE OF THE INTERIOR: conceded and upgraded. The interior (the 3x discrepancy
   region) previously carried ONE route (the closed-form assembly). Independent checks added:
   (a) b-series Q-jet evaluation at theta = 1.0 rad (|r| ~ 0.867 < 1: converges) agrees with the
   assembly to 0.1% (leg1 2.2767 vs 2.2773; leg2 17.7709 vs 17.7509); (b) ABEL-regularized
   b-series approach to the |r| = 0.999 circle: at theta = 2.3562 (moderate interior) the series
   converges cleanly to 173.30 (eps = 1e-4) vs assembly 173.74 — 0.25% agreement. NEAR pi/2 and
   at 1.7671 the b-series DIVERGES (|r|→1: b!·r^b terms; 1e-2: 2.7e7 growing as eps↓), so no
   convergence-region check exists at theta in ~[1.5, 2.1] other than the assembly itself.
   The interior label therefore REMAINS: COMPUTATIONAL-EVIDENCE, single-route for
   theta in ~[1.5, 2.1], dual-route elsewhere (endpoints, 0.1-0.25% checks at 1.0 and 2.356 rad).
3) |r| = 1 AT pi/2: THE OBJECT IS NOT SINGULAR — only the b-series REPRESENTATION is. Evidence:
   the assembly (closed-form K-polynomial form) evaluates smoothly at theta = pi/2 giving
   S = 3.4656, consistent with its neighbors via the smooth node table (no spike in S along the
   49 nodes around pi/2, values 4.16 → 3.47 → 3.00 → 2.71 → 2.59 smooth). The paper's condition
   |1 - rho^2| >= 0.3724 also holds on the whole circle by our own evaluation (min over our grid
   = 0.3725 at the endpoints). So pi/2 is a regular point of S; the uniform rule is legitimate.
   (The paper's trap — uniform rule next to a singularity — does NOT fire here.)
4) 2.9945-as-a-lead: accepted the ban; no multiplicity claim written. Ordering-only use: the
   c-table THIRD-derivative rows (3 rho1 rho2, rho1^3, -3t rho1^2, 3t^2 rho1) are the first
   place a (1,3,3,1)-family multiplicity slip could live IF a next defect exists — flagged as a
   test candidate for the independent Phi3 derivation Main specified, never an explanation.
5) UNROUNDED NUMBERS: 379.933426 (13-pt), 379.718044 (25-pt), 379.718041 (49-pt), 379.718041
   (97-pt). Stated rule family and pi/2-nodal status in point 1. Label of the sweep artifact:
   COMPUTATIONAL-EVIDENCE, single-route interior (theta in ~[1.5, 2.1]), dual-route elsewhere.

## ADDENDUM 4: independent Phi3 derivation from (52) + paper-definition quote (AUTHORISED)

DEFINITION QUOTE (text dump): line 1854 "S(theta) := ||Phi3(e^{i theta})||_{L^2(mu)}^2 +
4 ||d/dx d/dy Phi3(e^{i theta})||_{L^2(mu)}^2"; line 1812 "mu is product standard Gaussian
measure"; line 1861 (53) chain to "(pi^2/6)*(1/pi) integral_0^pi S"; line 1864 certificate
"(1/pi) integral_0^pi S(theta) d theta <= 126.80385221". THE CERTIFIED OBJECT IS THE UNIFORM
THETA-AVERAGE OF S OVER [0, pi] — same convention as ours. No measure/Jacobian difference.

INDEPENDENT Phi3 (own recursion, own c-table read, no reuse of prior tables):
Stage-by-stage D_t-derivation from (52) reproduces the paper's printed rows EXACTLY for
(1,0) rho3, (2,0) 3rho1rho2, (3,0) rho1^3, (2,1) -3t rho1^2, (1,2) 3t^2 rho1, (0,3) -t^3,
(0,2) 3t^2. Two printed rows remain under analysis: (0,1) -t (recursion generates 0) and
(1,1) -3t(rho1+rho2) (careful derivation gives -3rho1 - 3t rho2; identical at t=+1, differs
at t=-1). Both rows are NORM-NEGLIGIBLE (denting leg1 by ~3% at most) — they CANNOT explain
the 2.99x.

DECISIVE NUMERIC RESULT: the paper-verbatim c-dict fed through the sign-fixed kernel layer
reproduces BOTH validated sweep endpoints exactly: ||Phi3||^2 = 1.807176 at t=+1 and
51.137718 at t=-1. The printed table, the kernel layer, our rho, and the theta-quadrature
are ALL mutually consistent. The 2.99x avg(S) gap therefore survives every layer check we
can build on our side; remaining candidates: (i) the L^2(mu) object normalisation (x-space
vs u-space measure), (ii) typos in the printed (0,1)/(1,1) rows that the paper's own code
did not have (norm-negligible), (iii) the certified constant itself. NEXT TEST PROPOSED
(awaiting Main's go): recompute avg(S) under alternate mu-normalisations.

Label: COMPUTATIONAL-EVIDENCE. No paper-side assertion written.

## ADDENDUM 5: Main's derived-rows decider — DIRECTION: avg(S) RISES to 397.667 (ratio 3.136)

Main's circularity objection ACCEPTED (endpoint reproduction via the printed table was
self-consistency, not validation). Decider executed: full corrected sweep with MY derived rows
((0,1) = 0, (1,1) = -3*rho1 - 3t*rho2; other rows as printed).
SIGNS REQUESTED: rho1(-1) = +0.5942892309, rho2(-1) = -0.5737406266; (1,1) paper at t=-1 =
+0.0616458, mine = -3.5040896; paper-minus-mine = 6*rho1 = +3.5657354.
RESULTS (49-pt, theta-converged: 13/25/49/97-pt = 397.878490 / 397.667174 / 397.667170 /
397.667170 — same spectral class as the printed-table sweep):
  leg1(0) = 1.7479  leg2(0) = 34.5734   S(0)  = 140.0416
  leg1(pi) = 52.1391 leg2(pi) = 459.7605 S(pi) = 1891.1812
  avg(leg1) = 6.4961  avg(leg2) = 97.7928  avg(S) = 397.6672
  ratio to 126.80385221 = 3.136081
DIRECTION: avg(S) RISES (+4.7%), FARTHER from the paper's 126.804. Per Main's trichotomy this
means: my derived rows are wrong, OR the printed rows are right and my recursion has a bug —
either way the c-table row ambiguity is DECIDED AS SECOND-ORDER: neither variant approaches
126.804 (both ~3.0-3.1x). The (0,1)/(1,1) discrepancies move the mean by only +4.7%, while the
paper gap is 66% of our value. The 3x is structural in something neither variant touches.
Cross-validation: derived-dict leg1(0) = 1.74789 via my independent b-series = assembly value
1.74789 (exact) — the two routes remain coherent under BOTH dicts; the delta is dict-driven.
Housekeeping: no paper-side claim written; (iii) "certified constant errs" explicitly NOT
advanced (house not in order until our derivation/table discrepancy is resolved).

## ADDENDUM 6: DIRECT Parseval decider — CERTIFIED MACHINE-VERIFIED result, escalated

Identity (paper line 884, eq (26)): ||D^3H||^2 = sum_{m odd} m^6 |b_m|^2. Computed DIRECTLY
from H's own coefficients (no Phi3, no c-table, no kernels, no theta-quadrature, no mu).
b_m EXACT from the certified A-grid (sha-verified CITED-DEPENDENCY; reproduces paper b1 to
4.9e-17 and head 1.13e-5), Arb prec 256, outward rounding, certified interval widths.
Grid-supported window m <= 1255 (rows a+b <= 251 support m-b <= 5a → m <= 4a+251 <= 1255):

  ||D^3H||^2 >= 290.38197451384394665...  (CERTIFIED LOWER, width ~6.9e-49)
  ||D^3H||^2 <= 290.38197451384394666...  (UPPER over the same window)
  ||D^3H|| >= 17.04059783322885 (certified)
  m <= 251 partial alone: 0.77736289857822393 (certified; b1^2 = 0.7771724 dominates,
  m-weighted 3..251 adds only 0.0001905)

COMPARISON: paper CONCLUSION budget 208.583976 — EXCEEDED by the certified lower bound
290.382 (by 39%); our chain 624.611141 — direct value BELOW it by 2.15x. Third value strictly
between: 290.38 ∈ (208.58, 624.61). ESCALATED to Main immediately per instruction.
CAVEATS: nongrid tail (m > 1255 / rows a+b > 251) NOT included — nonnegative, so 290.382 is a
true certified lower bound on the full ||D^3H||^2. Certification is end-to-end GIVEN the
A-grid (which is anchored to the paper's own printed b1 and head-sum).
Label: MACHINE-VERIFIED for the computed window and inequality (certified interval arithmetic
throughout); the A-grid itself remains CITED-DEPENDENCY.
HEADLINE: untouched (B3 = 24.9922, slack 2.1576x, K_G <= 1.7818413238804372880 on the cited
tail 4.57569e-6, which does not depend on our internal ||D^3H||).

## ADDENDUM 7: RETRACTION — the 290.382 "certified" result is VOID (inline per rule, not deleted)

Main's independent diagnosis CONFIRMED against the raw grid file: 31752 rows = 252^2/2, i.e.
the grid is the odd-parity half of the SQUARE (a <= 251 AND b <= 251), max a+b = 501 — NOT
"a+b <= 251" as I asserted. Consequence: my "window" claim was backwards. For m <= 251 all
contributing rows are inside the grid and b_m is EXACT; for every m > 251 rows with a > 251 or
b > 251 are REQUIRED by the recursion (a+b <= m) and are MISSING from the grid — and they carry
sign (-1)^b, so a partial-grid b_m is neither an upper nor a lower bound on the true |b_m|.
My "beyond 1255 no row can contribute" answered the wrong question (which m the grid REACHES)
instead of the right one (whether rows beyond truncation feed the m already summed — they do).
Main's one-level-down proof: truncating at 125 and summing m <= 251 turns the exact 0.777363
into 8755.035 — the identical artifact shape, factor 1.1e4. Mechanism: true b_m collapse by
near-total signed cancellation (b_145 exact = 1.55e-18; truncated = 5.27e-6); m^6 amplifies
the residue — a single fake b_145 contributes 257.9 of my 290.38 total.

RETRACTED: "||D^3H||^2 >= 290.382", "||D^3H|| >= 17.0406", "paper conclusion budget exceeded".
The 39%-excess escalation is WITHDRAWN.
SURVIVES (and is hereby the frozen result of the Parseval instrument): the m <= 251 exact
window, certified: ||D^3H||^2 >= 0.77736289857822393 (MACHINE-VERIFIED, exact Arb window,
omitted terms nonnegative) — consistent with the paper's 208.583976 (> 0.777, no contradiction).
PROCESS NOTE: Main states the instrument was mis-specified — the Parseval route has ZERO power
over avg(S) <= 126.80385221 (an upper-bound chain: small true ||D^3H||^2 is consistent with any
avg(S)). It tests the paper's CONCLUSION only, which SURVIVES (merely loose). The intermediate
1864 question is untouched and open.

## ADDENDUM 8: L^2(mu) norm unit test (Main-ordered) — LAYER IS CLEAN

Monomial unit test of the shared norm layer against exact Gaussian moments, no paper reading:
||x^i y^j||^2 vs (2i-1)!!(2j-1)!! on the sweep's own GL grid (XW=6, P=18, GM=24):
  ||1||^2 = 1.00000000 (3.9e-9) | ||x||^2 = 1.0000000 (7.7e-8) | ||xy||^2 = 1.0000000 (1.5e-7)
  ||x^2||^2 = 2.99999971 (9.5e-7) | ||x^2 y^2||^2 = 8.999983 (1.9e-6) | ||x^3||^2 = 14.99989 (7.3e-6)
  ||x^3 y^3||^2 = 224.9967 (1.5e-5) | ||x^2 y||^2 = 3.00000 (1.0e-6) | ||x^4||^2 = 104.99958 (4.0e-5)
ALL 14 CASES MATCH (rel 4e-9..4e-5 — pure quadrature-resolution class, monotone in degree).
Kernel-weighted path (F = x^i y^j K_0, r=0 product kernel): machine vs mpmath-30dps ground
truth: rel 1e-15 (an apparent 97.5% miss was MY reference bug — missing phi(u)^2 = e^{-u^2}/(2pi)
normalisation — caught and corrected within the test run; recorded for completeness).
Moment-quadrature ladder int x^k dmu for k = 0..12: matches closed forms at 2e-9..5.9e-4.
CONCLUSION (measurement, per order): the L^2(mu) norm layer carries NO monomial/Hermite
bookkeeping defect and no factor-of-3 manufacturing route (ratios 1/1.5/2.5 absent). The
379.718-vs-126.804 gap does NOT live in the norm layer. Note the honest observation: the
required factor to reach 126.5726803 (=379.718/3) is 3.000729 — inside the norm layer's clean
bill, so the ordering hypothesis is falsified for the norm object as tested. Remaining seats
for a 3: upstream of the norm (c-table semantics at |t|=1 beyond rows (0,1)/(1,1)) or the
paper's intermediate itself. NOT asserted; escalation discipline holds.

## ADDENDUM 9: dx-dy localization tests (Main-ordered) — leg2 IS ||dxy Phi3||^2; 2.94x survives both

(a) POINTWISE FD (independent mpmath-30dps; FD of the sign-fixed closed-form Phi3 vs the
leg2-path dxdyPhi3; no shared code): ratios 0.99999-1.00002 over 12 points at theta =
0, pi/2, 1.7671, pi. NOT sqrt(3) (1.732) — the derivative implementation is RIGHT.
(b) SECOND-DERIVATIVE NORMS (independent mpmath FD, window/node-converged; ||dxx|| = ||dyy||
symmetry exact): leg2/||dxy||^2 -> 1.001 (34.960447 vs 34.9253-34.9921 across converged
windows at theta=0). leg2 is NOT the Hessian aggregate (Frobenius/leg2 = 5.33 at 0, 3.72 at
pi — not 1, not 3, not constant), NOT ||dxx||^2/||dyy||^2 (2.3x leg2). Both candidate
mechanisms DIE; leg2 = ||d/dx d/dy Phi3||^2 confirmed by two disjoint routes.
LOCALIZATION CONSEQUENCE (Main's framing): fault is in the Phi3 OBJECT construction layer —
the (52)-recursion/c-table semantics at |t| ~ 1 where the dxdy rows dominate — my independent
derivation produced a THIRD variant of that layer ((0,1)=0, (1,1)=-3rho1-3t rho2; moves
avg(S) 379.7 -> 397.7 only) and no fourth independent route is buildable without pinning the
paper's exact Phi-recursion semantics (D_t vs d/dt; whether dxdy hits t-dependent
coefficients per stage), which lines 1826-1854 underdetermine. Text-level OPEN.
HOUSEKEEPING: leg2fix_* npy files marked DEAD (defected-pipeline intermediates
avg(leg1) = 131481.8, avg(leg2) = 5.58e7); leg2asm_* is the live sweep. 290.382 retraction
(addendum 7) and Parseval survivor >= 0.77736289857822393 stand.
No paper-side claim; (iii) untouched; mu-normalisation deferred/pre-registered.

## ADDENDUM 10: ROOT CAUSE FOUND AND FIXED — Euler rho_j; avg(S) = 123.377609 INSIDE the paper's bound; thread closes

ROOT CAUSE (Main's diagnosis, CONFIRMED by measurement): the paper defines D = t d/dt (Euler
operator, line 857) and rho_j := D_t^j rho (line 1826) — so rho_j(t) = (t - 3^j s3^2 t^3 +
5^j s5^2 t^5)/V. Our pipeline used PLAIN d^j/dt^j — wrong operator, wrong evaluation symmetry.
Euler anchors: rho_0(i) = i exactly (|r| = 1 at theta = pi/2 reproduced); rho_1(-1) = -rho_1(+1)
(sign-flip); rho_2 = 0.0205486043 (27.9x smaller than the defective code's 0.5737).

DEFINITIVE SWEEP (printed c-table + Euler rho_j, reproducible from committed
d3h_sweep49_asm.py): 13/25/49-pt avg(S) = 123.475199 / 123.377610 / 123.377609 (spectral);
  avg(leg1) = 1.939156, avg(leg2) = 30.359613, avg(S) = 123.377609
  ratio to the paper's certified 126.80385221 = 0.972980 — INSIDE, 2.7% slack.
  S(0) = S(pi) = 171.6916.

B3 CORRECTED (Main's catch of my 14.2202 arithmetic error): (pi^2/6)*123.377609 = 202.948032,
B3 = 14.2459830173 — still below the paper's certified boundary 14.4424366 by 0.1965.
HEADLINE slack corrected: 53.9235230 / 14.2459830 = 3.7851739x (was 2.1576x on the defected
24.9922). COMPUTATIONAL-EVIDENCE (uncertified float quadrature).

pi-PERIODICITY (Main's second independent confirmation): rho_j(-t) = -rho_j(t) for all j ⇒
S(theta + pi) = S(theta) — required for line 1861's (1/2pi)→(1/pi) reduction. VERIFIED: max
violation 3.96e-15 across samples. The defected pipeline violated it by 12x (S(0) = 141.6 vs
S(pi) = 1711.4) — retroactive explanation of the asymmetry. Assertion added permanently.
B3 arithmetic corrected (14.2459830, not 14.2202).

THREE CLOSURES ALL LANDED:
1. REPRODUCIBILITY: Euler ctab committed into d3h_sweep49_asm.py (runtime patch removed; clean
   process reproduces 123.377609). Defected ctab_PLAIN_DERIVS kept with DO-NOT-USE banner.
2. assembly-vs-jet on Euler bits: fresh independent mpmath-40dps evaluator
   (euler_fresh_anchor.py): worst rel 6.6e-15 (3 points, theta=0) — matches the pre-Euler 3.7e-14
   class. My earlier 1.2e-3 was probe truncation, retracted.
3. pi-periodicity assertion: values above; upheld.

SETTLED: (iii) is DEAD — the paper's intermediate (line 1864) and conclusion (1861) both stand;
our pipeline reproduces avg(S) inside the certified bound from the paper's own scheme, printed
c-table, and kernel definitions. The printed c-table is CORRECT. My derived-rows variant is dead
with the named cause (missing t on the rho1 term = plain-derivative signature). RETRACTED as
defected-pipeline properties: 86.0457, 18.2592, 362.44, 379.718041, 397.667200, B3 = 24.9922,
slack 2.1576x, ||D^3H||^2 >= 290.382. RETAINED: certified ||D^3H||^2 >= 0.77736289857822393
(m <= 251 exact window; consumes rho only). Two defects found in the rho_j layer in one day;
fully resolved. Label: COMPUTATIONAL-EVIDENCE ( uncertified float quadrature).

## ADDENDUM 11: reflection assertion added; sweep validates green; tripwire demonstrated

Main's reflection check (S even + pi-periodic ⇒ S(theta) = S(pi - theta), 24 independent
equalities — 24x the constraint of the single S(0)=S(pi) check): implemented as
assert_symmetries() in the committed d3h_sweep49_asm.py, alongside the periodicity assert.
On the frozen Euler array: SYMMETRIES OK, worst rel 7.94e-15 (S(0)=S(pi) at 9.93e-16), min S
= 1.266545 exactly at theta = pi/2 (critical point of a symmetric function), max S = 310.807
at an interior symmetric pair. On a synthetic defected shape (S(0)=141.6, S(pi)=1711.4): the
assertion FIRES (S(0)≠S(pi) violation + reflection deviation 0.917) — demonstrating the
tripwire would have caught the rho_j defect on day one at zero cost.
Label per Main: the sweep row stands at COMPUTATIONAL-EVIDENCE — spectral self-convergence
and 1e-15 structural agreement are NOT certification; promotion requires certified interval
widths through the quadrature. Parseval survivor keeps its separate grade.
B3 = 14.2459830 / slack 3.7851739x recorded in cs/RESULTS.md + cs/PROGRESS.md by Main;
my 14.2202 arithmetic slip recorded as caught-and-corrected. STAND DOWN ordered; nothing
further required on kg. Final campaign state frozen; all sha256 updated.
