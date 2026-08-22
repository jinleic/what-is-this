# UNRESOLVED: the full translation-invariant scalar Kac--Ward family at order eight

## Abstract

**[COMPUTATION]** The full translation-invariant scalar family has 30 directed-direction-pair weights, or 25 coordinates on a generic five-edge gauge chart. We formed the exact trace equations at orders `k=4,6,8` simultaneously on the different free boxes `3x3x2` and `2x2x3`. We also enumerated all 32 zero/nonzero patterns of the five weights fixed by the wave-4 gauge tree, scanned every resulting polynomial system for variables that occur only linearly, and recorded exact Jacobian ranks over three finite fields. **These calculations do not decide the family.** Both a ten-box and the minimal two-box modular Gröbner computation timed out after 3600 seconds; a timeout is not a no-go theorem. A numerical point satisfying the ten-box construction equations fails an independent exact `3x3x3` prediction at `k=8`, falsifying that point but not the variety. Thus the full family remains **UNRESOLVED**.

Artifacts: `experiments/e30_kac_ward_full.py`, `results/kac_ward/full_family.json`, and `tests/test_kac_ward_full.py`.

## 1. Normalisation and exact two-dimensional control

Let

`P_G(v) = sum_{E even} v^|E|`

be the exact even-subgraph polynomial. The scalar Kac--Ward target is

`det(I-v Lambda) = P_G(v)^2`.

**[EXTERNAL — Kac--Ward/Whitney]** This identity holds for finite planar graphs with the signed planar half-angle weights.

**[COMPUTATION]** Before studying 3D, the experiment reconstructs `det(I-v Lambda)` exactly in `Z[zeta_8]` and compares all coefficients with the square of the repository's independent exact even-subgraph enumeration. On free square lattices `3x3`, `3x4`, and `4x4`, the maximum coefficient discrepancy and every non-rational cyclotomic component are exactly zero. This is a setup control, not evidence for a 3D identity.

## 2. Full translation-invariant scalar equations

Write the six cubic directions as `d in {+x,-x,+y,-y,+z,-z}`. A nonbacktracking translation-invariant scalar transition assigns one complex number `U(d,d')` to each allowed ordered pair `d' != -d`, hence 30 weights. Direction-state similarity

`U(d,d') -> g(d)^(-1) U(d,d') g(d')`

has effective dimension five. On a chart where the following tree weights are nonzero, they can be normalised to one:

`U(+x,+y), U(+x,-y), U(+x,+z), U(+x,-z), U(+y,-x)`.

This leaves 25 variables. Newton's identity gives the necessary trace equations

`F_(B,k)(U) := Tr_B(U^k) + 2 k [v^k] log P_B(v) = 0`.

**[LEMMA]** Any scalar determinant identity on every free box must satisfy these equations on each particular free box and at every order.

**[COMPUTATION]** The minimal stored construction system is

`F_(3x3x2,k) = F_(2x2x3,k) = 0`, for `k in {4,6,8}`.

It consists of six explicit polynomials over `Q` in 25 named variables. Their term counts are `7,41,350,7,37,179`. The complete expressions and a SHA-256 of their canonical SymPy representation are in `data.minimal_unresolved_system` of the JSON artifact, so the system can be handed unchanged to a stronger elimination engine.

**[COMPUTATION]** The minimal construction variety is explicitly nonempty after reduction
modulo each of `3,5,7,11`: the JSON stores one complete 25-coordinate solution over each
finite field, and direct substitution gives six zero residues. This is stronger than a sampled
Jacobian rank and rules out a modular emptiness certificate at those primes. It does **not**
prove a characteristic-zero solution: an `F_p` point need not lift to `Qbar`. Independent
`3x3x3` substitution rejects the stored witnesses at one or more orders for `p=5,7,11`; the
`p=3` witness passes that one holdout, which remains only a finite-field fact.

## 3. Linear elimination and its branch guards

**[COMPUTATION]** Across a larger span-separating set of free boxes, row reduction leaves 20 independent equations through `k=8`. Greedy exact elimination found six pivots before expression growth became prohibitive. The eliminated variables were

`u_mx_my, u_my_mz, u_mx_mx, u_mx_mz, u_my_my, u_mz_mz`.

The corresponding divisions require nonzero pivot guards. On that guarded chart, the system reduces to nine equations in 19 variables, with as many as 311 terms in one remaining equation. Each failed pivot guard creates a separate zero branch. **[COMPUTATION]** Continuing to a seventh pivot exceeded a 3600-second resource limit. **[NOT PROVED]** This partial elimination neither establishes consistency nor emptiness on the guarded chart, and says nothing by itself about the pivot-zero branches.

## 4. Zero-pattern stratification actually covered

The earlier generic gauge chart omitted points at which a fixed tree weight vanishes. We stratified those five particular coordinates by a bit string, in tree order, with `0=vanishes` and `1=normalised to one`. This yields all `2^5=32` patterns.

**[COMPUTATION]** For every pattern, the experiment substitutes the five values exactly into the ungauged six-polynomial system, records the active-variable count, and counts variables that occur with degree at most one in every remaining equation. It then evaluates the exact symbolic Jacobian at deterministic points over `F_101`, `F_1009`, and `F_10007` and row-reduces exactly modulo each prime. All 32 branch records, including their zero pattern and method, appear in the JSON artifact.

**[LEMMA — scope of the rank diagnostic]** The rank at a sampled point is only local differential information about the polynomial map at that point. Because the sampled point need not lie on the variety, it proves neither existence nor emptiness. Agreement across three primes is a cross-check, not a rational equality theorem and not a branch decision.

Consequently **no branch was discarded by the Jacobian diagnostic**. Every one of these 32 records is honestly marked `unresolved`. Also, these are the 32 patterns of one selected five-edge tree; they are not all `2^30` zero patterns of the original weights. A complete projective/gauge stratification would require covering all direction-state support graphs and their valid spanning-tree charts.

## 5. Modular elimination resource wall

**[COMPUTATION]** Two exact modular Gröbner attempts were made over `F_101`, both with grevlex order:

1. 30 equations in 25 variables from ten free boxes through `k=8`: timeout after 3600 seconds.
2. The minimal six equations in 25 variables from `3x3x2` and `2x2x3`: timeout after 3600 seconds.

**[NOT PROVED]** Neither timeout is evidence of a nonempty variety, and neither is a no-go certificate. No rational Gröbner call was made after the modular prefilter failed to decide the system. No Nullstellensatz certificate `sum_i g_i F_i=1` was found.

## 6. Independent prediction from a numerical survivor

To ensure that construction-box fitting was not mistaken for a solution, we sought a point numerically and reserved a box not used in the fit.

**[NUMERICAL]** Complex float64 least squares found a generic-chart point with maximum scaled residual `6.37e-17` on ten construction boxes at `k=4,6,8`. This is not an exact algebraic assignment. On the independent free `3x3x3` box, the `k=4` and `k=6` trace-equation residuals remain at floating roundoff, but the `k=8` residual is

`-24.39243951057324 + 98.57593668106370 i`,

of absolute value about `101.549`.

The right-hand side used here comes from exact even-subgraph enumeration of the holdout box. Therefore this particular numerical branch candidate is falsified at `k=8`. **[NOT PROVED]** Failure of one fitted point does not imply that the algebraic construction variety is empty or that every point fails the holdout.

## 7. Exact scope and conclusion

**[COMPUTATION]** Established here:

- exact 2D control identities on `3x3`, `3x4`, and `4x4` free squares;
- the explicit six-polynomial, 25-variable two-box system over `Q`;
- explicit `F_3`, `F_5`, `F_7`, and `F_11` solutions of the six construction equations, with independent `3x3x3` residues;
- enumeration of all 32 zero/nonzero patterns of the previous five-edge gauge tree;
- exact per-pattern degree scans and modular Jacobian ranks at three primes;
- a numerical construction point and its independent exact-target `3x3x3` failure at `k=8`.

**[EXTERNAL/UNCHANGED]** The prior exact no-go for the strict cubic-covariant two-weight slice at `k=8` remains valid.

**[NOT PROVED]** This work does not prove existence or nonexistence in the full 30-weight translation-invariant scalar family, does not decide arbitrary support/zero patterns among all 30 weights, and does not address projective weights, position-dependent weights, spin-structure sums, matrix-valued weights, or any order above eight.

**Conclusion: the full translation-invariant scalar Kac--Ward family at `k>=8` remains UNRESOLVED.** The reusable advance is the explicit minimal system, branch inventory, guarded partial elimination record, and precisely measured resource wall.

## 8. Reproduction

```text
.venv/bin/python experiments/e30_kac_ward_full.py
.venv/bin/python tests/test_kac_ward_full.py
```

Both commands print a final `PASS` on success. The experiment overwrites `results/kac_ward/full_family.json` using exact arithmetic except for the section explicitly labelled numerical.
