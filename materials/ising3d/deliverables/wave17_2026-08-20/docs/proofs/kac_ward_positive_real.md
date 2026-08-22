# Positive-real scalar Kac--Ward no-go from one cube

Artifacts: `experiments/e71_kw_positive_real.py`, `results/kac_ward/positive_real.json`, and `tests/test_kw_positive_real.py`.

## 1. Exact claim and domain

**[DEFINITION — existing parameterization]** Let the six direction states be

```text
+x, -x, +y, -y, +z, -z.
```

The full translation-invariant scalar ansatz assigns a number `U(d,d')` to each ordered pair with `d' != -d`. There are 30 such pairs: six straight transitions and 24 orthogonal turns. An immediate reversal is absent. The finite-box matrix has the form `Lambda(v)=v U`, where the one common Ising edge variable is `v=tanh(K)`. Thus finite ferromagnetic coupling `K>0` gives `0<v<1`. There are no separate anisotropic edge variables in the full-family parameterization used here.

**[DEFINITION — domain decided here]** The strict positive-real domain is

```text
0 < v < 1,
U(d,d') in R and U(d,d') > 0 for all 30 allowed ordered pairs.
```

The proof actually covers the closed enlargement `U(d,d')>=0`. In the strict domain every transition weight is nonzero and `U=conjugate(U)` pointwise.

**[SCOPE — no invented conjugation assumption]** The 30-weight source parameterization does **not** impose a relation between distinct coordinates such as

```text
U(-d',-d) = conjugate(U(d,d')).
```

That path-reversal relation is therefore not assumed below. If a physical convention additionally imposes it, its positive-real points form a subset of the domain excluded by the theorem. Conversely, positivity of the Ising edge variable `v` alone says nothing about the signs or phases of `U`; complex half-angle weights such as those in planar Kac--Ward are outside this positive-real corner-weight theorem.

**[THEOREM]** No assignment in the nonnegative-real 30-weight domain satisfies the coefficientwise scalar Kac--Ward identity on every finite free cubic box. In fact, the single `k=4` trace equation on the free `2x2x2` box has no such assignment.

## 2. Relation to the full-family chart and branch map

**[DEFINITION — existing gauge]** Direction-state similarity acts by

```text
U(d,d') -> g(d)^(-1) U(d,d') g(d').
```

The existing gauge tree, in its stored order, is

```text
+x->+y, +x->-y, +x->+z, +x->-z, +y->-x.
```

Its five underlying edges form a spanning tree on the six direction states. On the strictly positive domain a positive choice of the six `g(d)` normalizes all five entries to one. Hence every strict positive-real point lies on branch `11111` of the stored map. The other 31 labels contain at least one zero among these selected entries and are boundary strata, not points of the strict domain.

**[COMPUTATION]** The source branch artifact has 32 labels from `00000` through `11111`, with `0=vanishes` and `1=normalized to one`. It proves 17 branches empty over `Q` and leaves 15 unresolved over `Q/C`. The experiment in this wave reproduces the tree and all 32 labels exactly; it does not alter those characteristic-zero statuses.

**[THEOREM — stronger boundary statement]** After substituting the five tree entries according to any one of the 32 labels and requiring every remaining coordinate to be nonnegative real, the same rational Positivstellensatz argument below excludes the specialized equation. Thus all 32 selected nonnegative strata are excluded. This is a semialgebraic result and does not claim that the 15 unresolved polynomial ideals are empty over `Q` or `C`.

## 3. Exact `2x2x2`, order-four equation

For a finite free box `B`, coefficientwise equality

```text
det(I-v U_B) = P_B(v)^2
```

requires, by the formal logarithm/Newton trace identity,

```text
F_(B,k)(U) := Tr_B(U^k) + 2 k [v^k] log P_B(v) = 0.
```

**[COMPUTATION — exact integers]** Independent even-subgraph enumeration on the free cube gives

```text
P_2x2x2(v) = 1 + 6 v^4 + 16 v^6 + 9 v^8,
[v^4] log P_2x2x2(v) = 6.
```

Introduce the six oriented elementary-plaquette products

```text
M_xy+ = u_px_py u_mx_my u_py_mx u_my_px,
M_xy- = u_px_my u_mx_py u_py_px u_my_mx,
M_xz+ = u_px_pz u_mx_mz u_pz_mx u_mz_px,
M_xz- = u_px_mz u_mx_pz u_pz_px u_mz_mx,
M_yz+ = u_py_pz u_my_mz u_pz_my u_mz_py,
M_yz- = u_py_mz u_my_pz u_pz_py u_mz_my.
```

(The ordering of commuting factors is immaterial.) Exact closed-walk enumeration yields

```text
Tr_2x2x2(U^4) = 8 (M_xy+ + M_xy- + M_xz+ + M_xz- + M_yz+ + M_yz-).
```

Therefore

```text
F_(2x2x2,4)
 = 48 + 8 (M_xy+ + M_xy- + M_xz+ + M_xz- + M_yz+ + M_yz-).       (1)
```

Every displayed coefficient and the constant 48 are integers reconstructed directly by `full_weight_finite_system(..., gauge_fix=False)`. No floating-point value, critical-coupling benchmark, optimization, or interval truncation enters (1).

**[LEMMA — immediate inequality]** If all 30 coordinates are nonnegative, every `M` in (1) is nonnegative, so

```text
F_(2x2x2,4) >= 48 > 0.
```

This already contradicts the necessary equation `F_(2x2x2,4)=0`.

## 4. Rational Positivstellensatz certificate

The inequality proof has an explicit rational preordering identity, not merely a sign inspection.

Let `T` be the preordering generated by the 30 inequalities `u_i>=0`. Each plaquette monomial is a squarefree product of four generators and hence belongs to `T`. Moreover

```text
1/6 = (1/3)^2 + (1/6)^2 + (1/6)^2.                              (2)
```

Combining (1) and (2) gives the exact identity in `Q[u]`

```text
-1 = -(1/48) F_(2x2x2,4)
     + (1/6) (M_xy+ + M_xy- + M_xz+ + M_xz- + M_yz+ + M_yz-).  (3)
```

The first term on the right belongs to the ideal generated by the required equality `F=0`; by (2), every remaining term is an explicit rational sum-of-squares multiplier times a squarefree product of inequality generators. Thus the right side lies in `ideal(F)+T`. At a common point of `F=0` and `u_i>=0`, (3) would assert `-1>=0`, which is impossible.

**[THEOREM — rational certificate]** Identity (3) is a rational Positivstellensatz infeasibility certificate for the full nonnegative orthant. The JSON stores the ideal multiplier, each source monomial, its odd-generator product, its even square-monomial part, and all rational squares in (2). The standalone test reconstructs the two sides as SymPy polynomials over `QQ` and verifies exact equality to `-1`.

### Branchwise specialization

For branch `b`, substitute each selected tree weight by its bit `0` or `1`. A zero kills any plaquette monomial containing it; a one removes that factor. Write the surviving specialized monomials as `M_(b,j)`. The constant never changes, and all surviving coefficients remain 8, so

```text
F_b = 48 + 8 sum_j M_(b,j),
-1 = -(1/48) F_b + (1/6) sum_j M_(b,j).                          (4)
```

Every `M_(b,j)` remains a squarefree product of nonnegative active coordinates, so (2) turns every term in (4) into an explicit preordering term. The artifact stores a separate specialized polynomial and certificate for each of the 32 labels. The test freshly substitutes the bits and rechecks every one of the 32 identities; no branch result is inferred only by numerical sampling or symmetry.

## 5. Exact sign, parity, and conditional conjugation consequences

The cube equation also supplies necessary constraints outside the positive cone.

**[LEMMA — real signed parity]** For arbitrary real, possibly signed, nonzero `U`, equation (1) requires

```text
M_xy+ + M_xy- + M_xz+ + M_xz- + M_yz+ + M_yz- = -6.            (5)
```

Consequently at least one oriented plaquette product is negative. A product of four nonzero real factors is negative exactly when an odd number of its factors is negative. Thus every real signed solution must carry an odd sign parity on at least one oriented elementary plaquette. The all-positive sign sector has even (zero) negative parity on every plaquette and is excluded.

**[LEMMA — conditional reversal conjugation]** Under the additional, explicitly conditional relation

```text
U(-d',-d) = conjugate(U(d,d')),
```

the two orientations in each plane are complex conjugates. With representatives `H_xy=M_xy+`, `H_xz=M_xz+`, and `H_yz=M_yz+`, equation (5) becomes

```text
Re(H_xy) + Re(H_xz) + Re(H_yz) = -3.                            (6)
```

In particular, positive-real holonomies cannot satisfy (6).

**[LEMMA — conditional unit-modulus corollary]** If one further assumes every corner weight has modulus one, then each `|H_plane|=1` and `Re(H_plane)>=-1`. Equality (6) forces all three lower bounds to saturate, hence

```text
H_xy = H_xz = H_yz = -1.
```

Unit modulus is **not** an assumption of the positive-real no-go. This corollary records exactly what the order-four coefficient would demand from a phase-valued, reversal-conjugate family; it does not exclude such a complex family at higher orders.

## 6. Result and limitations

**[THEOREM]** The strict positive-real scalar family is empty: its only gauge-tree chart is `11111`, and the ungauged rational identity (3) already excludes it on the free cube at order four. More strongly, the closed nonnegative orthant and all 32 of its stored tree specializations are empty for this equation, with exact margin 48.

**[UNRESOLVED]** This result does not prove that the 15 source branches lacking characteristic-zero certificates are empty over `Q` or `C`. Signed-real weights can make plaquette products negative, and complex phases can meet the order-four sign requirement. The theorem also does not address position-dependent weights, matrix-valued weights, spin-structure sums, projective constructions, or equality at only one fixed numerical value of `v` rather than coefficientwise polynomial equality.

**[COMPUTATION]** No positive-real candidate or numerical optimization was run, because the exact rational certificate decides the whole stated domain. Therefore no high-precision `3x3x3` holdout is applicable.

## 7. Reproduction

```text
.venv/bin/python experiments/e71_kw_positive_real.py
.venv/bin/python tests/test_kw_positive_real.py
```

**[COMPUTATION]** Each command prints a final `PASS` if and only if all of its checks pass. The experiment writes the provenance/data/checks envelope to `results/kac_ward/positive_real.json`; the standalone test reloads it and independently rebuilds the cube equation and all 32 branch certificates.
