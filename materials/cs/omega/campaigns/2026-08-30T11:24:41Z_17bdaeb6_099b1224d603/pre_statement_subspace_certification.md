# PRE-STATEMENT — SUBSPACE CERTIFICATION (gate C stage 2)
# Campaign base: 2026-08-30T11:24:41Z_17bdaeb6_099b1224d603 (Main GO,
# 2026-08-30T13:1xZ, four conditions). Committed BEFORE any certification
# computation.

## 0. Scope honesty, FIRST (Main condition 1)

The certified endpoint at p* is 2.3715538358544617; the published rung-2
value is 2.37155181 (gap 2.02585446e-6). Measured float sensitivities of
the slack model along the 21 exact kernel directions at eps = 1e-7 are
−2.5e-9..+6.8e-9 per-direction worst-case, and the certified LINEAR
model's maximal swing on the largest fully-feasible coefficient cube
(side r_c = 3.36306043896845e-08, inscribed in the 1e-7 box in every
coordinate simultaneously) is ±1.06e-8, i.e. a best possible certified
improvement of ~1.34e-9. **Reaching the published value at this
sensitivity would need eps ~ 8.1e-5 — three orders outside the
pre-registered 1e-7 box.** THIS CAMPAIGN CANNOT AND WILL NOT AFFECT THE
PUBLISHED 2.37155181 VALUE, the record 2.371177, or any published number.
Its value is the CERTIFICATION (first rigorous local statement over a
full feasible subspace in this line), not any number it produces.

## 1. Object certified

The certified-aggregation endpoint Omega(c) along the EXACT 21-dimensional
kernel K of the one-block equality machinery (basis V: 45x21, SVD of the
rank-24 margin system, simplex-orthogonal to machine precision; saved
scratch/kernel_basis_V.npy, sha-verified at load). Coefficient space: c in
R^21, point p(c) = p* + Σ_j c_j V[j,:] restricted to
|V c|_inf <= 1e-7 (per-coordinate box radius) AND strict simplex support.
Feasible coefficient domain: D = {c : max_i |(Vc)_i| <= 1e-7}
⊇ cube C = {c : |c_j| <= r_c = 1e-7 / 2.9734820950964815 =
3.36306043896845e-08} (row-sum bound, computed and recorded).

## 2. Certified decomposition

Omega(c) = Omega0 + L(c) + R(c), where
- Omega0 = B0 certified endpoint = 2.3715538358544612 (re-evaluated in
  this campaign; frozen 2.3715538358544617),
- L(c) = certified linear model: interval directional derivatives
  D_j = ∂Omega/∂c_j evaluated by INTERVAL central differences inside
  IC.prec() through the same dual-mode tree (MID=300), each bounded
  [D_j.lo, D_j.hi] outward; linear model interval
  L(c) = Σ D_j c_j (interval product, outward),
- R(c) = certified remainder over the WHOLE cube C: from interval
  Hessians. Entropy terms: d²(−t ln t)/dt² = −1/t (certified interval
  enclosure per coordinate over [t_lo, t_hi] on C — support min on C
  computed and recorded); p_comp terms: low-degree polynomials of the
  splits (exact polynomial Hessians, bounded over C by interval
  substitution). Lipschitz constant Jac(c) for the gradient over C:
  ||H(c)||_1 bound (interval), giving
  |R(c)| <= (1/2) * 21 * r_c² * max||H||_1 = recorded number.
- Certified endpoint range on C:
  Omega0 + [min_c L(c) − |R|max, max_c L(c) + |R|max] with c ranging over
  the FULL cube C (c_j = ±r_c independently, 21 dims — certified over the
  whole cube, not sampled).

## 3. Coverage (condition 2, hard requirement)

The certificate MUST cover C = the full 21-dimensional cube. No
sub-sampling, no "most promising directions". If the remainder bound
cannot be closed over C, the outcome is REPORTED AS FAILURE TO CERTIFY
(the scoped negative from the parent campaign stands as final) — not
patched by shrinking the domain.

## 4. Checkpoint (condition 3)

HARD STOP for the remainder-bounding step at T+60 minutes from first
Hessian-bounding run. If the Hessian-Lipschitz bound is not closed by
then: stop, report failure-to-certify + scoped negative as final, write
the checkpoint state to disk first (write-first protocol).

## 5. Outcome handling (condition 4)

Both outcomes are results and BOTH escalate to Main before any README or
index wording:
- certified cap over C (Omega range ≥ B0 or bounded below by B0 within
  tolerance): first rigorous local-optimality-over-subspace statement in
  this line;
- certified feasible improvement (Omega range strictly < B0 − 0): improves
  our own enclosure (by ≤ ~1.34e-9 per §0) — still not a published-value
  claim.

## 6. Instruments (independence discipline)

- Kernel basis: SVD (scipy) of the margin matrix — independent of the
  accessor and of any probe machinery;
- directional derivatives and Hessian bounds: interval evaluation through
  the dual-mode tree, per-direction, inside IC.prec();
- cross-check: the float gradient (saved scratch/kernel_grads_float.npy,
  computed 2026-08-30T13:2xZ, all |g_j| <= 0.04) must be INSIDE its
  interval counterpart for every j — a disagreement is a stop-and-report;
- no further probe-direction sampling anywhere in this campaign.


## AMENDMENT A (Main-approved 2026-08-30T13:5xZ, BEFORE compute restart) —
## parameterization changed from coefficient cube to the natural polytope

Rationale (recorded): the rref basis is unnormalized — row-sums of |V|
land at ≈ the nullspace dimension (margin factor EXACTLY 21), a basis
artifact, not geometry; the inscribed cube r_c = 1e-7/21 covers 7× less
than the SVD-basis cube and ~all of the feasible region is lost.

NEW PARAMETERIZATION: delta = V·c in exact-kernel coordinates; certified
domain D = { delta : A·delta = 0 exactly (kernel), |delta|_inf ≤ 1e-7 }
= kernel ∩ box — the natural set. Certified enclosure:
  Omega(p* + delta) ∈ Omega0 + [LPmin, LPmax] + [−Rmax, +Rmax] for ALL
  delta ∈ D, where
- [LPmin, LPmax] = range of the certified linear term g·delta over D,
  obtained by LP min/max with CERTIFIED duals: solve float (scipy), then
  verify weak duality in interval arithmetic on the dual certificate
  (y for the 24 equalities, mu± for box bounds) per dual-feasible
  direction (2 LPs, one per endpoint);
- Rmax = box-wide certified Hessian–Lipschitz remainder over
  |delta|_inf ≤ 1e-7 (valid and conservative: the kernel constraint can
  only shrink the achieved remainder);
- named domain in every result sentence: D in delta-coordinates (not a
  cube; no inscribed-loss disclosure needed — D IS the feasible region
 .kernel∩box).

ALL OTHER CONDITIONS UNCHANGED: whole-subspace coverage (condition 2) is
now over D directly; T+60 hard stop for the remainder step (condition 3)
starts at the first Hessian-bounding run; both outcomes escalate before
wording (condition 4); scope-honesty paragraph (§0) unchanged and
binding — best possible improvement ~1.34e-9 at OLD-basis numbers; at
the polytope domain the linear swing is bounded by
r_box * Σ_j |g_j| · (exact-basis row-sum normalization) — recomputed and
recorded at LP time; reaching published 2.37155181 remains impossible at
this sensitivity.

MOMENT-STRUCTURE FACT (input to any future campaign): the three margin
row-dependencies are — shared grand totals (pairwise) and the two
first-moment identities Σ_j (8−j)-type — all structural (no duplicated
or mis-transcribed row; 27 rows pairwise distinct). Recorded in
retraction_addendum_rigidity.md.
