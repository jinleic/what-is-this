# Pre-statement — oct-rank campaign (agent OctRank), committed 2026-08-29

**Committed BEFORE any computation.** Thresholds and budgets below may not be
adjusted after seeing results. Any re-scoping is recorded in an appended
"Re-scoping log", keeping the original statement visible.

## Primary sources (first-hand reads completed before this statement)

* arXiv:2608.16649v1, Hardik Jain, "Bounds on the real tensor rank of octonion
  multiplication" (submitted 2026-08-17). Read in full from the arXiv HTML.
* github.com/hxrdxkxvxd/octonion-rank, commit `816a01e16acb5849ba5ad9858e80cbb5a39b967d`
  (2026-08-19 "Added arXiv link"), SHA256 hashes of all cert/script files
  recorded in the campaign artifact.
* Fiduccia & Zalcstein 1977: NOT read first-hand → `[REPORTED]` (via paper § 1
  which quotes their Thm 6 and Example 4). **No gate depends on it.**
  Cariow & Cariowa: same, `[REPORTED]`.

## Facts fixed by the paper (target of re-verification)

* Lower bound chain: substitution peeling (LM-style) ×6 leaves pencil
  (L_u, L_v), u,v independent; C = L_u^{-1} L_v satisfies an irreducible
  quadratic M² − 2⟨u,v⟩M + N(u)N(v)I = 0 with negative discriminant; pencil
  bound (Thm 2): rank(I,C) ≥ n + n/2; at n=8: 6 + 12 = 18.
* `τ` = 3-slice quaternion family (I, L_i, L_j) on ℝ⁴: lower 1+6 = 7; upper
  by the same Krawczyk route at ρ = 10⁻⁵, paper announces worst margin
  ≈ 9.9 × 10⁻⁶ (> 0), K ≈ 0.005352 (this second number from repo README).
* Upper bound for T_O: rank-25 CP residual F: R⁶⁰⁰ → R⁵¹²; 88 coordinates
  frozen (S = 512 free indices); certificate = (A, B, C, Y, S) dyadic
  rationals; hypotheses over the box B_∞(x₀, ρ), ρ = 10⁻⁶:
  (i) K := max row-sum of |I − Y J_S| + ρ|Y|W₁ + ρ²|Y|W₂, K < 1;
  (ii) ‖Y F(x₀)‖_∞ ≤ (1−K)·ρ;
  (iii) per-equation inclusion radii: ρ − (off_e + ρ·rs0_e + ρ²·rs1_e +
  ρ³·rs2_e) > 0 with rs0 = row sums of |I − Y J_S|, rs1, rs2 = row sums of
  |Y|W₁, |Y|W₂ restricted to S.
  Paper's announce-pass numbers: K ≈ 0.3196, off_max ≈ 1.3×10⁻¹⁴,
  (1−K)ρ ≈ 6.8×10⁻⁷, worst per-equation margin 6.8×10⁻⁷.

## Gate A — re-verify the published rank-25 certificate [DECIDED AHEAD]

**Fixed quantities (same as paper; I am re-verifying, not inventing):**
ρ = 1/10⁶ exactly. All certificate data read as exact dyadic rationals
(float64 → Fraction, denominators must be powers of two; assert this).
Hypotheses to pass, all evaluated in EXACT rational arithmetic:
A1. Structure tensor T_O equals the standard Cayley–Dickson octonion product
    (independent construction in my code, entries in {−1,0,1}); L_{e0} = I.
A2. Y·J_S ≈ I consistency: ‖Y·J_S − I‖_∞ < 10⁻⁶ (exact rational check;
    this is an archive-consistency check, not a hypothesis of Prop 6).
A3. K < 1 with K as above.
A4. off_max ≤ (1−K)·ρ (exact; paper: factor ≈ 5×10⁷ to spare).
A5. Every per-equation margin strictly positive (512 margins).
**Pass verdict for Gate A** = A1–A5 all hold in my independent python-flint
implementation, AND my independent re-derivation of the Krawczyk radius
(see below) confirms that a ρ of 10⁻⁶ works, AND Lean accepts the lower
bound theorems (Lean replay separate). **Refutation verdict** = A3, A4 or A5
fails in outward-rounded interval re-evaluation, or the Lean statements are
materially weaker than the paper's prose. That triggers the escalation rule
(scored immediately to Main via hub, before any file records it).

**Independent Krawczyk-radius re-derivation (Gate A, second half).** Using
my OWN interval arithmetic (python-flint `arb`, outward rounding), with the
certificate's Y, x₀ and the exact Jacobian Hessian-free derivative enclosure
|J_S(x) − J_S(x₀)| ≤ ρ W₁ + ρ² W₂ (paper's degree-2 enclosure), I
recompute:
  K_dim := max row-sum of |I − Y J_S| + ρ|Y|W₁ + ρ²|Y|W₂ (must be < 1),
  margin_e := ρ − (|YF(x₀)|_e + ρ rs0_e + ρ² rs1_e + ρ³ rs2_e) (must be > 0),
and additionally estimate the largest ρ* = max{ρ: hypotheses hold} by bisection
with OUTWARD rounding at every product/sum (each bound computed as an `arb`
ball whose errors are rounded up into the comparison). Pass: ρ* established
> 10⁻⁶ with margin; the value near 6.8×10⁻⁷ per-equation worst is
[REPRODUCED] if I recover ≈ that figure, [DERIVED] for my bisection number.
Precision discipline: every comparison is ball > ball-zero with ERROR SIGN
EXPLICIT (outward); a tie (ball straddles zero) is recorded as
NOT-CERTIFIED-BY-THIS-SCALE and escalated upward precision (fixed 512-bit
effective precision steeply exceeds dyadic input granularity — not expected).

## Gate B — rank-24 attack [DECIDED AHEAD]

**HARD SOUNDNESS RULE (from README, accepted):** real tensor rank is a
semialgebraic condition over R. NO enumeration/sampling of integer or rational
witnesses can ever establish a LOWER bound. A failed rank-24 search yields NO
claim about rank > 24. Only a SUCCESS certified by Krawczyk/interval-Newton
yields an upper-bound claim, and even that must pass the same 3-hypothesis
scheme at fixed ρ before being claimed.

**Search budget (fixed before launch):**
* Engine: CP-ALS (alternating least squares on the 8×8×8 tensor) followed by
  Gauss–Newton / Levenberg–Marquardt refinement of the rank-2r least-squares
  residual, in float64, with jax-free numpy/scipy; threads pinned to 1.
* r shown in the search: 24 first, then 23, 22, 21, 20 (do not go below 20).
* Seed count per r: 256 distinct seeded starts (LCG/PCG64, fixed base seed
  20260829 per r), each start = random factor triple drawn from
  N(0,1)/√n scaled by (‖T‖F^{1/3});
* Iteration cap per start: 2000 ALS sweeps (counter-counted), then up to 300
  Gauss–Newton rounds;
* Stop-early: relative residual < 10⁻¹⁴ (assumed success candidate) or 8 h
  wall-clock per r, whichever first; global wall-clock cap 24 h then halt;
* Numerial-success criterion going into Krawczyk: min over starts of F-norm
  residual of ≤ 1e-13; the best start is then SNAPPED to dyadic rationals and
  a fresh Y preconditioner is computed from the float64 Jacobian restricted to
  good free columns (88 pivots chosen by a greedy criterion I will fix in
  code: maximum available 88 pivots, the exact rule is the standard
  "complete-to-invertible" rule = any set making J_S invertible with
  condition number as small as I can get by a greedy sweep);
* Certification at r=24 uses the SAME 3-hypothesis test at ρ = 10⁻⁶ first;
  if it fails at 10⁻⁶, attempt ρ = 10⁻⁵, 10⁻⁴ (bigger boxes sometimes help
  when x₀ is very close to x⋆) — all in exact arithmetic, recorded. A
  success at ANY of these radii at r=24 is a certified upper bound
  R_R(T_O) ≤ 24 (which improves the 2026 state of the art of 25);
  escalation: message Main BEFORE writing the claim (per assignment).
* If every start fails (min residual stays > 1e-13), Gate B terminates with
  verdict "FAILED (no rank-24 candidate found)" and NO claim about rank.

## Gate C — τ and sharpness [DECIDED AHEAD]

C1. Re-verify R_R(τ) = 7 in both directions:
    (a) lower: independent reproduction of the peel+pencil argument at n=4
      (human-audited algebra, plus my own exact-arithmetic check of the
      needed identities on the τ slices — exact rational arithmetic);
    (b) upper: replay of the τ Krawczyk certificate in my python-flint
      implementation with ρ = 10⁻⁵ fixed, same 3 hypotheses; paper: worst
      margin ≈ 9.9×10⁻⁶ (i.e. anchored at the boundary ρ ≈ (1−K)⁻¹·off_max
      much larger than needed), threshold FIXED at ρ = 10⁻⁵ and margin > 0.
    Pass = both directions hold ⇒ R_R(τ) = 7 re-verified.
C2. Sharpness of (5/2)n − 2 at n = 2 and n = 4 by exact computation:
    n = 2 (C): R_R(T_C) = 3. Exact check by hand: slices (I, J), J² = −I,
    pencil bound with n=2 gives ≥ 2 + 1 = 3 and the known rank-3 witness
    gives ≤ 3; both directions are 4×4-matrix algebra — I check the
    existence of an actual rank-3 decomposition over Q by exact linear
    algebra and the lower bound by the same pencil argument specialized.
    n = 4 (H): R_R(T_H) = 8. Same: lower by peel(2)+pencil(6) = 8; upper by
    exact integer witness decomposition of the 4-slice quaternion tensor
    into 8 rank-one terms — constructed and CHECKED in exact arithmetic.
    These two sharpness facts are `[REPRODUCED]` from classical ranks with
    my own certificate checks.
C3. n = 8 constant improvement probe (HUMAN-AUDIT REQUIRED BEFORE ANY
    CLAIM): analyze whether the peeling step could stop later (peel only 5
    slices, leaving a 3-slice family whose rank I can lower-bound by more
    than 12 + extra) or whether the pencil bound at n=8 can exceed 12
    (e.g. via #{linearly independent columns of D} arguments in
    Quaternionic structure). Any candidate improvement is written only
    under `scratch/`, tagged HUMAN-AUDIT-PENDING, and messaged to Main; no
    campaign artifact claims a better constant without independent human
    sign-off recorded there.

## Precision and environment

* python-flint 0.9.0 (`fmpq_mat`, `fmpz_mat` exact; `arb` balls, outward).
* venv: /Users/jinleic/jinleic-workspace/cs/.venv/bin/python (3.14.3).
* BLAS threads pinned to 1 (OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1).
* All certs loaded from a copy I take of the upstream repo at commit
  `816a01e…`, with SHA256 hashes recorded; upstream is never mutated.

## Ordering of work

1. pre_statement.md (this file) — committed first, before any run.
2. Gate A: independent T_O + inequality recheck; independent radius
   derivation; Lean 4 replay (install toolchain v4.29.1 if needed, `lake
   build`, `#print axioms` audit); freeze verdict(s) in campaigns/.
3. Gate B: search + certificate attempt at r = 24, budget as fixed above.
4. Gate C: τ both directions, C2 sharpness, C3 scratch note.

## What is NOT claimed by this campaign (stated ahead of time)

* No failure of any rank-24 (or lower) search implies anything about
  R_R(T_O) ≥ 24 or any lower claim. It only means this method produced no
  certified witness at that rank within the budget.
* No claim of an improvement over 18 on the lower side is made by me; any
  candidate goes to scratch + Main at HUMAN-AUDIT-PENDING.
* No claim that 25 (or 24 if reached) is optimal.
* The rank-25 upper bound will be carried as MACHINE-VERIFIED
  **re-verification** only; the underlying mathematics remains
  CITED-DEPENDENCY of arXiv:2608.16649 and its Lean artifacts.

## Gate C addendum C4/C3 — n=8 peeling-constant probe [DECIDED AHEAD, 2026-08-30]

(Committed BEFORE the C4 computations are run. C3 of this statement and
"gate C4" (README honest-gaps item) are the same work item; this addendum
fixes the plan. The same run batch also implements the lower halves already
committed under C1(a) and C2 — no re-scoping.)

Object: the constant 18 in 18 ≤ R_R(T_O) as produced by the
substitution(peeling) + Thm-2(pencil) machinery of arXiv:2608.16649
(Lean: Oct18.lean chain). Goal: decide whether any parametric choice inside
this machinery lifts the constant. Expected outcome is a negative, filed with
an exact rule-7 scope statement. Any candidate improvement is NOT claimed
here: freeze evidence, message Main via hub, stop.

Hypotheses of the machinery being probed (as consumed by the Lean chain):
H1. T_O slices are L_{e_p}, p = 0..7; all corrected slice families arising in
    the induction are Φ(w_a) with the w's linearly independent in ℝ⁸, so no
    corrected slice vanishes (Lmat8_eq_zero_iff).
H2. The substitution lemma (LM Prop 3.1) applies at each peel: pivot slice
    nonzero, over the field ℝ (Lean: substitution_pivot_last).
H3. Every 2-slice residual (L_u, L_v), u,v independent, satisfies the matrix
    quadratic M² − 2⟨u,v⟩M + N(u)N(v)I = 0 with M = L_ū L_v (composition
    identity + polarization; Lean Oct18.lean steps (4)–(6)).
H4. The discriminant is strictly negative: ⟨u,v⟩² < N(u)N(v) (strict
    Cauchy–Schwarz for independent u,v), so the normalized
    J := c⁻¹(C − a·1) satisfies J² = −1.
H5. The pencil floor consumes ONLY the normalized complex structure: by
    direct read of the Lean source (Pencil.lean), `key_bound` takes exactly
    `J * J = −1` plus the two decompositions — no other structure on C.

Pre-registered search space (all four probed; every value reported, hit or
miss):
S1. Pencil floor > 12? Can pencilRank(1, C) exceed n + n/2 = 12 for any C
    the chain consumes? Route: build an exact 12-term witness for
    (I_8, J_8), J_8 = ⊕₄ [[0,1],[−1,0]], by direct-summing the exact 3-term
    (I_2, J_2) witness already machine-verified in this campaign
    (gate_c_sharpness_n2_n4.py); verify entrywise in fmpq. Universal lower
    ≥ 12 is Lean key_bound; standard real canonical-form conjugacy
    [INFERENCE] places every J with J² = −1 in the GL-orbit of J_8, and
    pencil rank is invariant under the conjugation (u ↦ Pu, v ↦ vP⁻¹).
    Expected verdict: floor exactly 12 on the full consumed class ⇒ no gain.
S2. Stopping point: peel j slices (0 ≤ j ≤ 6), bound the residual (8−j)-slice
    L-family by the same induction engine at k = 6−j:
    r ≥ j + (12 + (6−j)) = 18 for every j. Machine-check the arithmetic for
    all j = 0..6 and report all. Expected verdict: constant 18 at every
    stopping point ⇒ no gain.
S3. Three-slice residuals: peel 5 to a 3-family (L_u, L_v, L_w); machinery
    yields 5 + 13 = 18. Is rank ≥ 14 for 3-slice L-families true? NOT swept
    in this campaign — no computation here addresses it; recorded as the
    open gap that would lift 18 to 19 if established.
S4. Kernel-dimension (#{independent columns of D}) arguments: verify on the
    S1 witness that every inequality inside key_bound is TIGHT
    (U D = 0, ker D = n/2 = 2 exactly, ker U = r − n = 4 exactly), so no
    D-side slack exists to exploit. Expected verdict: all tight ⇒ no gain.

Pass/fail (decided ahead): S1/S2/S4 all "no gain" AND S3 untouched ⇒ verdict
NEGATIVE, filed with the rule-7 sentence: the negative covers ONLY the
substitution/peeling framework terminating in a Thm-2 pencil at n = 8; NOT
swept: any non-peeling lower-bound route for T_O, and any rank bound on
3-slice L-families (L_u, L_v, L_w) beyond the peel-implied 13 — a bound ≥ 14
would improve 18 and is not excluded by this work. Any actual improvement
candidate (S1 > 12 on the consumed class, S2 < 18 at some j, S3 ≥ 14
evidence): freeze exact artifact + hub message Main + STOP, no claim in
README.
