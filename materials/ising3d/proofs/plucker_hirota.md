# Exact boundary-tensor Pluecker/Hirota certificate

## 1. Result and scope

[THEOREM] For the open `2 x 2 x 2` repository lattice with terminal order

\[
(0,0,0),\ (0,0,1),\ (1,1,1),\ (1,1,0),
\]

[THEOREM] For this ordered terminal set, the standard four-terminal Pfaffian Pluecker residual of the normalized Walsh transform of the exact fixed-boundary Ising partition tensor is

\[
R_{\rm cube}(v)=4v^6(1-v^2)^6.
\]

[THEOREM] For the open `2 x 2 x 3` repository lattice with the four terminals in cyclic order around the rectangular `x=0` surface patch,

\[
(0,0,0),\ (0,0,2),\ (0,1,2),\ (0,1,0),
\]

[THEOREM] For this ordered terminal set, the same residual is

\[
R_{\rm slab}(v)=4v^8(1-v^2)^8(1+v^2)(3+6v^2-v^4).
\]

[THEOREM] Both residuals are strictly positive for `0 < v < 1`; hence neither is the zero polynomial, and each is an exact finite counterexample to this particular four-channel matchgate identity on the declared boundary tensor and terminal order.

[UNRESOLVED] These counterexamples are not a no-go for an arbitrary determinant or Pfaffian representation, a different local basis change, another terminal order, a higher Pluecker/Dodgson system, a full octahedron recurrence, or an all-size representation of the three-dimensional Ising partition function.

## 2. Exact boundary tensor

[LEMMA] Let `G=(V,E)` be one of the finite open lattices returned by `ising.lattices.hyperrect`, let `B=(b_1,b_2,b_3,b_4)` be its ordered terminal set, and let bit `1` denote spin `+1` as in the repository API.  The producer constructs the integer-polynomial tensor

\[
T_\tau(v)=\sum_{\sigma_{V\setminus B}}
\prod_{\{x,y\}\in E}(1+v\sigma_x\sigma_y),
\qquad \tau\in\{\pm1\}^4.
\]

[LEMMA] If `v=tanh K`, then the physical fixed-boundary partition function is

\[
Z_\tau(K)=(\cosh K)^{|E|}T_\tau(\tanh K),
\]

[LEMMA] Therefore the omitted factor is common to every tensor entry and cannot create or remove a homogeneous quadratic identity.

[COMPUTATION] `e209_plucker_tensor.py` builds every `T_tau` by exact spin enumeration with Python integers; the largest declared case has only `2^12=4096` full spin states, and no benchmark value or floating-point arithmetic enters.

## 3. Walsh-current lemma

[LEMMA] For `S` contained in the four terminals, define

\[
C_S(v)=2^{-|V|}\sum_{\tau\in\{\pm1\}^4}
\left(\prod_{i\in S}\tau_i\right)T_\tau(v).
\]

[LEMMA] The transform equals

\[
C_S(v)=\sum_{F\subseteq E:\ \partial F=S}v^{|F|},
\]

[LEMMA] Here every interior vertex has even degree in `F` and the odd-degree terminal set is exactly `S`.

[LEMMA] Expanding the edge product labels each selected edge set `F` by `v^{|F|}`.  Summing a spin at vertex `x` gives zero unless `deg_F(x)` has the parity prescribed by the Walsh character, and otherwise gives `2`.  The product of the `|V|` surviving spin sums cancels the normalization `2^{-|V|}`.

[COMPUTATION] Every artifact case compares this fixed-spin/Walsh construction coefficient-by-coefficient with a second exact edge-parity dynamic program.  All five comparisons pass, all odd-cardinality current entries vanish, global spin flip pairs tensor entries, and every tensor entry at `v=0` equals `2^{|V|-4}`.

## 4. Identity and sign discipline

[LEMMA] The Pfaffian of a `4 x 4` antisymmetric matrix has the sign pattern

\[
\operatorname{Pf}(A)=a_{12}a_{34}-a_{13}a_{24}+a_{14}a_{23}.
\]

[LEMMA] Accordingly, the tested four-terminal matchgate/Pluecker relation is

\[
C_\varnothing C_{1234}-C_{12}C_{34}+C_{13}C_{24}-C_{14}C_{23}=0.
\]

[COMPUTATION] The coefficient vector is ordered as

\[
(C_\varnothing C_{1234},C_{12}C_{34},C_{13}C_{24},C_{14}C_{23})
\]

[COMPUTATION] This ordered coefficient vector equals `(1,-1,1,-1)` before any three-dimensional tensor is evaluated.

[COMPUTATION] The bounded search class is exactly

\[
(1,c_{12|34},c_{13|24},c_{14|23}),\qquad
c_{12|34},c_{13|24},c_{14|23}\in\{-2,-1,0,1,2\}.
\]

[COMPUTATION] The declared class therefore contains `5^3=125` normalized integer vectors and nothing else.

[COMPUTATION] Independent edge variables on the planar `2 x 3` training rectangle select the unique vector `(1,-1,1,-1)` in this class.  The `3 x 3` planar rectangle, the elementary cube, and the slab consume no coefficient-selection role.

## 5. Controls, holdouts, and certificates

| role | lattice and ordered terminals | exact residual | status |
|---|---|---|---|
| [COMPUTATION] planar training | `2 x 3`, four corners cyclic | `0` with independent edge variables | selects the unique bounded coefficient vector |
| [COMPUTATION] fresh planar holdout | `3 x 3`, four corners cyclic | `0` both multivariately and after uniform specialization | passes without fitting |
| [COMPUTATION] scope counterexample | `2 x 2 x 2`, four terminals cyclic on one face | `0` | three-dimensional geometry alone does not force failure |
| [THEOREM] minimal skew-cube certificate | `2 x 2 x 2`, two opposite parallel terminal edges | `4 v^6(1-v^2)^6` | nonzero for `0<v<1` |
| [THEOREM] surface-slab holdout | `2 x 2 x 3`, corners of one rectangular surface patch | `4 v^8(1-v^2)^8(1+v^2)(3+6v^2-v^4)` | nonzero for `0<v<1` |

[THEOREM] At the exact audit point `v=1/2`, the two nonzero values are

\[
R_{\rm cube}(1/2)=\frac{729}{65536},\qquad
R_{\rm slab}(1/2)=\frac{2329155}{268435456}.
\]

[LEMMA] Cube positivity follows directly from `4v^6(1-v^2)^6>0` on `0<v<1`.

[LEMMA] Slab positivity follows because all displayed factors are positive and, with `z=v^2` in `(0,1)`, `3+6z-z^2 > 3 > 0` after using `6z-z^2=z(6-z)>0`.

## 6. Multivariate cube witness

[COMPUTATION] The cube uses one independent variable per repository edge, in the exact order returned by `hyperrect((2,2,2), periodic=False).bonds`.  A product of two current monomials is encoded in base three, with digits in `{0,1,2}` recording the edge exponents.

[COMPUTATION] The full multivariate cube residual has `128` nonzero monomials and canonical digest

```text
bcd4ef8a2f5d5628b5ea684a092687e81f9f527a3bbea59d53e11ce1b19410a4
```

[COMPUTATION] Its first minimum-degree witness has base-three code `79572`, coefficient `2`, total degree `6`, and edge-exponent vector

```text
[0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]
```

[COMPUTATION] This exponent vector is in repository bond order.

[COMPUTATION] Uniformly setting all edge variables to `v` reproduces coefficient-for-coefficient `4v^6(1-v^2)^6`; the degree-six coefficient `4` is the sum of two distinct degree-six multivariate witnesses of coefficient `2`.

## 7. Exact bounded-class conclusion

[THEOREM] Within the explicitly declared set of `125` normalized four-channel integer coefficient vectors, no vector passes both the independent-edge planar `2 x 3` training tensor and the skew-terminal elementary cube.

[THEOREM] This is an ansatz-bounded finite no-go: the planar training tensor leaves only `(1,-1,1,-1)`, and its cube residual has the displayed nonzero monomial witness.

[UNRESOLVED] Coefficients outside the declared box were not searched, although the standard Pfaffian sign vector was also fixed directly by the `4 x 4` Pfaffian identity.

[UNRESOLVED] No inference is made from finitely many planar controls to a theorem about every planar Ising boundary tensor.

## 8. Counterexamples to overstatement

[COMPUTATION] The cofacial cyclic terminal order on the same elementary cube has zero residual.  It is therefore a concrete counterexample to the statement that every boundary tensor drawn from a three-dimensional lattice must fail this relation.

[COMPUTATION] The nonzero result is in the fixed normalized Walsh/current basis.  It does not imply that the raw spin tensor, or a tensor after another invertible local transform, has no determinant/Pfaffian realization.

[COMPUTATION] The multivariate witness prevents the cube failure from being an artifact of cancellations caused by uniform edge specialization, but it says nothing about identities with additional channels or larger minors.

[UNRESOLVED] No thermodynamic limit, critical coupling, free energy, magnetization, or critical exponent is inferred.  The benchmark `K_c=0.221654626` was not used to select, fit, or evaluate anything.

## 9. Reproduction

[COMPUTATION] The exact producer commands are:

```bash
nice -n 10 .venv/bin/python experiments/e209_plucker_tensor.py
nice -n 10 .venv/bin/python experiments/e210_plucker_search.py
nice -n 10 .venv/bin/python experiments/e211_plucker_certificate.py
```

[COMPUTATION] The standalone clean-room verifier command for the lead is:

```bash
nice -n 10 .venv/bin/python tests/test_plucker_hirota.py
```

[COMPUTATION] The machine-readable artifact is `results/integrability/plucker_hirota.json`; its top-level keys are exactly `meta`, `data`, and `checks`, and every stored passing check is computed by `e211_plucker_certificate.py`.
