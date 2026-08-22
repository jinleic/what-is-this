# Boundary stratification of the Kac--Ward branch `11111` anchor-section census

Artifacts: `experiments/e146_kw_boundary.py`,
`results/kac_ward/boundary_lifts.json`, and `tests/test_kw_boundary.py`.
This note continues `proofs/kw_components.md` (wave-14) and
`proofs/kw_cubic.md` / `proofs/kw_branch11111.md`; it analyses exactly the
part of the anchor thin-torus section census that wave-14 left unanalysed:
the **non-unit strata**.

## 1. Scope, chart, and exact computational convention

**[COMPUTATION]** The setup is that of `proofs/kw_components.md` §1 and is
used verbatim: the gauge-tree chart `11111` (five `DIRECTION_GAUGE_TREE`
weights set to 1), the 42-equation construction catalogue (orders `4,6,8` on
the fourteen shapes), the thin-torus substitutions, diagonals set to one, and
the six primitive quotient numerators at indices `28,29,32,35,38,41` (term
counts `33, 277, 155, 7, 155, 155`).  Two distinct uses of the word *chart*
occur below and are kept separate:

* **C13-coordinate charts** — the `C(13,6) = 1716` choices of six of the
  thirteen non-diagonal variables along which one may Hensel-lift;
* **the thin-torus weight chart** — the rational parametrization of the
  twenty-five weights by the sixteen free coordinates, defined only where the
  nine denominator monomials are units (this is the sense of
  "chart-degenerate" in `proofs/kw_components.md`).

All arithmetic is exact integer residue arithmetic over `Z/5Z` or `Z/625Z`.
The *anchor section* holds the seven trailing thin coordinates
(`u_pz_mx, u_pz_py, u_pz_my, u_mz_px, u_mz_mx, u_mz_py, u_mz_my`) at the
integers `(1,1,1,1,1,4,4)` and solves in the six remaining coordinates
(`u_py_px, u_py_pz, u_py_mz, u_my_px, u_my_mx, u_pz_px`).

## 2. The full-grid census

**[COMPUTATION]** (full `5^6` census).  Over the COMPLETE `5^6 = 15625`-point
anchor-section grid (solve coordinates in all of `F_5`, not just `F_5^x`), the
six primitive numerators have **exactly 1553** solutions; the unit stratum
(six solve coordinates in `F_5^x`) contributes exactly the three points
`p0, p1, p2` of wave-14.  Moreover the full 42-equation system was
thin-reduced equation by equation: **32 of the 42 numerators vanish
identically on the thin torus** (all primary-shape equations and all order-4
full-shape equations), and the ten nonzero numerators (the five full shapes
`(3,3,2),(2,2,3),(2,2,2),(3,2,2),(2,3,2)` at orders 6 and 8, indices
`28,29,31,32,34,35,37,38,40,41`) have **the same 1553-point zero set,
set-equal** to the primitive zero set.  The boundary decomposes into 43
distinct zero-coordinate patterns, the largest strata being `u_py_px = u_pz_px
= 0` (256 points), `u_py_px = u_my_px = 0` (112), and `u_py_mz = 0` alone
(exactly one point, `q = (1,3,0,3,2,3)`); there is exactly one point with all
six solve coordinates zero.  The census, the 42-numerator reduction, and the
set equality are recomputed independently in the standalone test.

## 3. Rank stratification of the 1553 solutions

**[COMPUTATION]** (rank histogram).  For each of the 1553 solutions the FULL
`42 x 13` Jacobian rank mod 5 was computed exactly (13 variables = 6 solve +
7 held; the 32 identically-zero rows contribute nothing, so the rank is that
of the ten nonzero rows).  At every solution this equals the rank of the
`6 x 13` primitive Jacobian.  The histogram (rank -> count) is

```text
0 -> 49     1 -> 349     2 -> 364     3 -> 630     4 -> 157     6 -> 4.
```

* **Maximal rank 6** (= the number of defining equations) is attained at
  exactly four points: the three unit points `p0, p1, p2` and the boundary
  point `q = (1,3,0,3,2,3)`.  At each of them the `6 x 13` primitive Jacobian
  has rank 6, so each is a **smooth complete-intersection point of the
  primitive variety `S`** (local dimension `13 - 6 = 7`); the solve-6 minors
  have determinants `2, 4, 1` (units) and `3` (boundary point), reproducing
  and completing the wave-14 record.
* The **other 1549 solutions are singular** (rank `<= 4 < 6`): rank `< 6`
  forces EVERY `6 x 6` minor of EVERY one of the `C(13,6) = 1716`
  C13-coordinate charts to be singular mod 5.  This is the exact
  chart-degeneracy certificate for the whole singular stratum; the
  rank-0 stratum (49 points, tangent space all of `F_5^13`) contains the
  all-zero solve point.

**[COMPUTATION]** (exhaustive chart census on the maximal-rank stratum).  For
each of the four maximal-rank points all 1716 C13-coordinate charts were
tested for a nonsingular `6 x 6` minor (primitive rows and full 42 rows; the
two counts agree).  Nonsingular-chart counts out of 1716:

```text
p2 = (3,2,3,1,1,1): 1153      p0 = (1,2,3,3,1,1): 1225
p1 = (2,4,1,3,3,1): 1023      q  = (1,3,0,3,2,3):  588.
```

So `q` is emphatically NOT chart-degenerate in the C13 sense — it is a smooth
point of `S` in 588 of the 1716 coordinate charts — which is exactly why its
exclusion needs the unit theorem of §4 rather than a Jacobian argument.

## 4. The unit theorem

**[THEOREM]** (all free coordinates of a thin-torus point are units).  Each of
the nine thin-torus substitutions is, by construction, the relation

```text
u_mx_my  * u_my_px                                 = -1
u_mx_py  * u_my_mx * u_py_px                       = -1
u_mx_mx  * u_px_px                                 =  1
u_mx_mz  * u_mz_px * u_pz_mx                       = -1
u_mx_pz  * u_mz_mx * u_pz_px                       = -1
u_mz_mz  * u_pz_pz                                 =  1
u_my_mz  * u_mz_py * u_py_pz * u_pz_my             = -1
u_my_pz  * u_mz_my * u_py_mz * u_pz_py             = -1
u_my_my  * u_py_py                                 =  1
```

between the twenty-five weights, with right-hand side the constant unit
`+-1`.  In any commutative ring `R` (in particular `Z/625Z`, `Z_5`, `Q_5`,
`Z`), if a product of ring elements equals a unit then every factor is a
unit.  The left-hand factors jointly involve **all sixteen free coordinates**
(each appears in at least one denominator monomial — a machine-checked
certificate recorded in the artifact).  Hence:

> Every point of `R^25` satisfying the nine thin-torus relations has all
> sixteen free coordinates invertible in `R`.

*Proof.*  Immediate from the nine displayed product identities and the
unit-factor property; the coverage of all sixteen free coordinates is
verified exactly in the producer and re-verified in the standalone test
(`union_covers_all_16_free_coordinates: true`).  `[THEOREM]`

**[COMPUTATION]** (empirical face).  Over the anchor section the theorem says:
every thin-torus chart point reduces mod 5 into the unit stratum.  Verified
pointwise: all **1550** boundary census points have at least one vanishing
thin denominator (monomial `= 0 mod 5`), and none of the three unit points
has any.

## 5. Hensel lifts of the maximal-rank stratum

**[THEOREM]** (lifts, conditional on nonsingularity).  Each of the four
maximal-rank points has nonsingular solve-6 minor, so by the multivariate
Hensel lemma it lifts uniquely in the solve-6 chart (seven held coordinates
pinned at the integers `(1,1,1,1,1,4,4)`) to a point of the six-primitive
system modulo `5^4 = 625`, computed by exact Newton iteration through mod 25,
125, 625.  `[COMPUTATION]` for the values:

```text
L0 = (226, 197,  98,  63, 161, 216)   (wave-14 anchor lift, reproduced)
L1 = (157, 459, 431, 228,  38, 496)   (reproduced)
L2 = (393,  72, 218, 596, 446, 496)   (reproduced)
Lq = (276, 278,   0, 288, 282,  83)   (the boundary lift; u_py_mz = 0 EXACTLY)
```

At `Lq` all **ten** nonzero thin-reduced numerators (not only the six
primitives) vanish mod 625, so `Lq` is a genuine mod-625 point of the full
42-numerator system; nevertheless `u_py_mz = 0` exactly, the denominator of
the relation for `u_my_pz` is a zero divisor mod 625, and **no 25-weight
thin-torus chart point lies over `Lq`** (the reconstruction of
`reconstruct_values25` is impossible — the assert `den % 5 != 0` fires by
theorem).  The same holds for any other C13-chart lift of `q` or of any of
the 1549 singular boundary points: a lifted coordinate that starts `= 0
mod 5` stays `= 0 mod 5` (Newton corrections are multiples of the running
modulus), and every solve coordinate occurs in some denominator monomial
(§4).

## 6. Holdout outcomes

**[COMPUTATION]** (controls, unit lifts).  The three unit lifts reconstruct
into the weight chart, satisfy all 42 construction residues `= 0 mod 625`,
and their `3x3x3` open-box holdout residues at orders `4,6,8,10,12` — computed
by BOTH the walk-monomial symbolic expansion and the `108 x 108`
matrix-power route, cross-validated equal — reproduce wave-14 exactly:

```text
L0 : (0, 0, 250,  55, 266)
L1 : (0, 0, 456, 350, 132)
L2 : (0, 0, 430, 480, 286).
```

First failing order 8 at every unit lift.

**[THEOREM]** (boundary holdout non-evaluability).  At `Lq` (and at every
conceivable lift of every boundary census point) the `3x3x3` holdout is
**not evaluable**: the holdout equation is a statement about the twenty-five
weights, the thin-torus chart map is undefined (denominator a zero divisor)
at any 13-tuple with a coordinate `= 0 mod 5`, and by the unit theorem no
other chart parametrization reaches these points.  The artifact records the
status `NOT_EVALUABLE_CHART_DEGENERATE` with the exact vanishing denominators
(`u_my_pz` for `Lq`) and the kept-zero coordinate (`u_py_mz`).

## 7. The boundary exclusion theorem (main statement)

**[THEOREM]** (anchor-section boundary exclusion).  The statement holds over
any commutative ring `R`, in particular over `R = Z/625Z`, `Z_5`, and `Q_5`.
Let `X` be a point of `R^25` satisfying the nine thin-torus relations whose
seven held thin coordinates are the integers `(1,1,1,1,1,4,4)` and whose ten
nonzero thin-reduced construction numerators vanish.  Then:

1. all sixteen free coordinates of `X` are units (§4);
2. reducing mod 5, the six solve coordinates of `X` are nonzero, hence the
   reduction is one of `p0, p1, p2` (§2 census);
3. if moreover the seven held coordinates equal the integers exactly, `X` is
   one of `L0, L1, L2` (Hensel uniqueness, §5), and each of these fails the
   `3x3x3` holdout at orders 8, 10, 12 mod 625 (§6).

Consequently **none of the 1550 non-unit solutions of the anchor-section
census — the 1549 singular ones and the smooth boundary point `q` — is the
mod-5 reduction of ANY thin-torus chart point over `Z/625Z` or `Q_5`**, and
the wave-14 scope note ("the 1549 remaining boundary solutions ... are
unanalysed here") is discharged: they are not shadow components, they are
outside the chart altogether, provably and uniformly in the ring.

## 8. Cross-artifact controls and provenance

The standalone test `tests/test_kw_boundary.py` (run:
`PYTHONPATH=src .venv/bin/python tests/test_kw_boundary.py`) imports no
producer module.  It independently rebuilds the 42-equation catalogue and all
42 thin-reduced numerators, recomputes the complete 15625-point census (1553,
unit 3, 42-numerator set-equal), the FULL rank histogram and per-point ranks
(the three unit points and `q` at rank 6, the all-zero solve point at rank
0), the boundary Hensel lift `Lq = (276,278,0,288,282,83)` with all ten
numerator residues zero and `u_py_mz` kept exactly zero, the anchor unit lift
`(226,197,98,63,161,216)` with holdout residues `(0,0,250,55,266)`
cross-validated walk == matrix, the unit-theorem certificate (constant unit
right-hand sides, monomial denominators, all-16 coverage, 1550/1550 boundary
vanishing denominators, 0/3 unit), and asserts agreement with
`results/kac_ward/boundary_lifts.json` and with the wave-14 records in
`results/kac_ward/components.json`.

Producer run time ~45 s, test run time ~42 s, both `time.process_time()`
gated, RSS well under 6 GB; no floating point decides any statement.

## 9. Limits and open questions `[UNRESOLVED]`

* The theorem is **section-local** (held coordinates pinned at the integers
  `(1,1,1,1,1,4,4)`) and **chart-local** (thin-torus locus, gauge branch
  `11111`).  The `4^7 - 1` other thin-torus sections remain open, as in
  wave-14.
* The boundary exclusion says nothing about points of the construction
  variety with **vanishing non-gauge weights** (zero-weight branches outside
  the thin-torus parametrization) or with vanishing gauge weights — the same
  caveat as the `full_weight_finite_system` docstring ("all zero-entry
  branches remain outside the computation and must not be claimed as
  excluded").  Whether such zero-weight points of the 25-weight variety exist
  over the section is a separate question.
* The characteristic-zero component structure of the full 42-ideal over the
  whole thin torus — and hence the final `Q`/`Q_5` decision of branch
  `11111` — remains `[UNRESOLVED]`; the present note removes the anchor
  section's boundary as a possible source of missing components.
