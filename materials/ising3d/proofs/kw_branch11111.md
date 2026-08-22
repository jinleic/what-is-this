# Exact degree-two and thin-torus investigation of Kac--Ward branch `11111`

Artifacts: `experiments/e80_kw_branch11111.py`,
`results/kac_ward/branch11111.json`, and `tests/test_kw_branch11111.py`.

## 1. Finite exact system and scope

**[COMPUTATION]** The calculation works on the generic gauge-tree chart `11111`:
the five selected tree entries are set to one, leaving 25 active coordinates.
It uses the exact finite-box equations at orders `4,6,8` on the 14 shapes

```text
permutations of (2,2,1), permutations of (3,2,1),
(3,3,2), (2,2,3), (2,2,2), (3,2,2), (2,3,2).
```

This gives 42 raw integer polynomials. Exact `Fraction` elimination of their
coefficient vectors gives a 21-dimensional `Q` span; the selected raw indices
are stored in the JSON. This is the same finite construction catalogue used by
the prior degree-one search, not an assertion about all boxes or all orders.

## 2. Multiplier-degree-two rational Nullstellensatz

For every one of the 21 independent raw constraints, multiply by every
monomial in the 25 active variables of total degree at most two. There are

```text
1 + 25 + binomial(25 + 1, 2) = 351
```

multipliers and therefore `21 * 351 = 7371` integer Macaulay columns. Their
canonical sparse hash is

```text
f35b6b5f6c40f7d990d5cd3fd24d15f1c07fe7ab6cfe506d35a0f3322ccc66c5.
```

**[THEOREM]** For this 42-constraint finite system, no rational
Nullstellensatz certificate with multiplier degree at most two exists.

**Exact proof.** Sparse Gaussian elimination over `Fraction` constructs 7365
nonzero pivots over `Q`, with pivot digest

```text
99e55bbf5ee2d62929ce343ff143b999fa31404c06564725fb8654785fbc08e8.
```

The six exact zero reductions occur at columns

```text
1739, 2457, 3861, 4898, 4914, 5099.
```

The JSON records each complete rational relation (each has 4 or 6 nonzero
terms), so they give the characteristic-zero upper bound
`rank <= 7371 - 6 = 7365`. The constructed pivots give the matching lower
bound. Reducing the target constant with the same exact basis leaves the
nonzero residual `1` at the all-zero 25-exponent monomial. Thus `1` is not in
the column span, proving the stated bounded absence result.

**[UNRESOLVED]** This theorem excludes only the displayed degree bound and
finite 42-constraint catalogue. It neither proves that the branch has a
characteristic-zero point nor that its full ideal is empty over `Q` or `C`.

## 3. Explicit elimination through the thin prime torus

Let the nine thin-radical binomials be

```text
u_mx_my*u_my_px + 1
u_mx_py*u_my_mx*u_py_px + 1
u_mx_mx*u_px_px - 1
u_mx_mz*u_mz_px*u_pz_mx + 1
u_mx_pz*u_mz_mx*u_pz_px + 1
u_mz_mz*u_pz_pz - 1
u_my_mz*u_mz_py*u_py_pz*u_pz_my + 1
u_my_pz*u_mz_my*u_py_mz*u_pz_py + 1
u_my_my*u_py_py - 1.
```

**[COMPUTATION]** The experiment independently solves these nine relations as
Laurent substitutions for the nine displayed leftmost variables. For example,
`u_mx_my=-1/u_my_px`,
`u_mx_py=-1/(u_my_mx*u_py_px)`, and
`u_mx_mx=1/u_px_px`. The remaining 16 coordinates are an explicit torus
chart. Direct exact substitution verifies that all 27 primary thin-box
constraints vanish on this chart.

**[COMPUTATION]** The 15 genuinely three-dimensional construction equations
are then reduced to Laurent polynomials in those 16 coordinates. Five
order-four rows vanish identically; the ten remaining cleared numerators have
exact `Q`-linear span rank 6, with independent raw constraint indices
`28,29,32,35,38,41`. The JSON stores each Laurent denominator, numerator
support size, maximum total degree, and a canonical sparse hash.

**[UNRESOLVED]** This completed coordinate elimination did not produce a
univariate consequence or a zero-dimensional elimination ideal. The six-rank
support calculation is structural information, not an emptiness certificate
or a construction of a complex point.

## 4. Reversal-conjugate fourth-root complex slice

Impose additionally the conditional relation

```text
U(-d',-d) = conjugate(U(d,d'))
```

and restrict every coordinate to the exact finite subset
`mu_4={1,i,-1,-i}` of `Q(i)`. The condition is an extra chart assumption; it
is not part of the 30-weight source parameterization.

**[LEMMA]** On the thin torus all six oriented elementary plaquette holonomies
are `-1`. Hence their sum is `-6`. Under the conditional reversal-conjugation
relation, the two orientations in each plane are conjugates, so each plane
has real holonomy part `-1`.

**[COMPUTATION]** Reversal and the five normalized tree entries leave 11 free
`mu_4` orbits, hence `4^11=4,194,304` raw phase assignments. The nine thin
relations become exact congruences modulo four. Unit-pivot row reduction has
pivot columns `(1,3,6)` and free columns `(0,2,4,5,7,8,9,10)`, reducing the
scan deterministically to 65,536 assignments. All 65,536 satisfy the nine
radical relations under exact Gaussian-integer evaluation.

**[COMPUTATION]** Every one of those 65,536 assignments fails at least one of
the 42 construction equations under exact Gaussian-integer evaluation. There
are zero construction-passing candidates, so no point is reported as
surviving and no `3x3x3` holdout candidate test is applicable. The complete
ordered first-failure digest is

```text
18ca797acc9d05c91c1b071f5d323253a182aa262dca7f5b6f9b9005e83d607d.
```

The standalone test independently rebuilds the phase chart, congruence
reduction, every finite candidate, and that digest.

**[UNRESOLVED]** The finite `mu_4` exhaustion is not a theorem about arbitrary
reversal-conjugate complex weights, much less arbitrary complex branch-`11111`
weights. It does not establish `Q`- or `C`-emptiness of the full branch.

## 5. Reproduction

```text
timeout 600 .venv/bin/python experiments/e80_kw_branch11111.py
timeout 600 .venv/bin/python tests/test_kw_branch11111.py
```

**[COMPUTATION]** Both commands use only exact integer, rational, Laurent, or
Gaussian-integer arithmetic for all deciding statements and print a final
`PASS` on successful reproduction.
