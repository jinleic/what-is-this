# Exact spectral non-Gaussianity of the open `3 x 3` layer

**[COMPUTATION] Artifacts.** The reproducible experiment is
`experiments/e54_spectral_3x3.py`, its exact certificate is
`results/spectral/spectral_3x3.json`, and the independent standalone check is
`tests/test_spectral_3x3.py`.

**[THEOREM] Finite `3 x 3` extension of Theorem S.** For each of the two exact
parameters

```text
t = tanh(K*/2) = 1/3, 1/2
```

the positive 512-dimensional symmetrized transfer operator of the open `3 x 3`
Ising layer is not a nine-mode positive-spectrum Gaussian subset-product
operator. This is an exact theorem about these two named finite rational
operators. In particular, the obstruction now holds for a layer graph having a
degree-four vertex.

## 1. Exact rational operator

**[LEMMA] Rational parametrization.** Put

```text
q = exp(2K) = (1+t^2)/(2t),
P_t[k,l] = t^hamming(k,l),
e_k = (b_k-eps)/2,
R = P_t diag(q^e_k) P_t.
```

Here `b_k` is the sum of the twelve in-layer bond products in computational
basis state `k`, and `eps` is their common parity. The symmetrized transfer
operator differs from `R` by one positive global scalar. Therefore all
multiplicative eigenvalue ratios used below are unchanged.

**[LEMMA] Single integer clearing.** If `t=a/b` is rational, multiplying `R` by

```text
D = b^(2n) q.denominator^max(max e_k,0)
             q.numerator^max(-min e_k,0)
```

produces an integer symmetric matrix `A=D R`. The implementation applies the
local two-by-two factors of `P_t` in `O(n 4^n)` integer operations rather than
forming a cubic dense product. For the two cases, the exact data are

| `t` | `q=exp(2K)` | exponent range | integer scale `D` |
|---:|---:|---:|---:|
| `1/3` | `5/3` | `[-6,6]` | `4412961507515625` |
| `1/2` | `5/4` | `[-6,6]` | `16777216000000` |

**[COMPUTATION] Positive definiteness.** Exact inertia at shift zero is zero in
all symmetry sectors at both couplings, certifying that every eigenvalue is
positive.

## 2. Exact `D4 x` spin-flip projectors

**[LEMMA] Central projectors.** Let `U_(g,f)` be the computational-basis
permutation for `g in D4` and global spin flip `f in C2`. For a real irreducible
character `chi_rho`, dimension `d_rho`, and spin-flip sign `s`, the operator

```text
P_(rho,s) = (d_rho/16) sum_(g,f) chi_rho(g) s^f U_(g,f)
```

is the orthogonal projector onto the corresponding `D4 x C2` isotypic
component. For the two-dimensional irrep `E`, an axial reflection `a` gives the
exact refinement

```text
P_(E,s,+/-) = P_(E,s) (I +/- U_a)/2.
```

The transfer operator commutes with both refined rows. This reduces each
120-dimensional `E` isotypic pencil to two independently counted
60-dimensional pencils without discarding multiplicity.

**[COMPUTATION] Projector certificate.** On all 512 computational-basis columns,
the twelve common-denominator-32 refined projectors satisfy, over the integers,

```text
M_alpha^2 = 32 M_alpha,
M_alpha M_beta = 0  (alpha != beta),
sum_alpha M_alpha = 32 I,
M_alpha^T = M_alpha.
```

All 66 unordered orthogonality pairs were checked. The matrix entries also
satisfy `A[g i,g j]=A[i,j]` for every one of the sixteen `D4 x C2` actions;
hence every projector commutes with `A` exactly. The certified sector
dimensions are

```text
A1+ 51, A1- 51, A2+ 19, A2- 19,
B1+ 33, B1- 33, B2+ 33, B2- 33,
E(axis+,+) 60, E(axis-,+) 60,
E(axis+,-) 60, E(axis-,-) 60.
```

Their sum is exactly 512.

**[LEMMA] Generalized sector inertia.** Select exact independent integer columns
of each projector numerator and place them in `C_alpha`. Define

```text
B_alpha = C_alpha^T A C_alpha,
G_alpha = C_alpha^T C_alpha.
```

The columns are an exact basis of the projector range, so `G_alpha` is positive
definite. The restriction of `A-sigma I` is congruent to
`B_alpha-sigma G_alpha`. Consequently the number of full-matrix eigenvalues
below a rational shift `sigma` is the sum of the negative inertias of these
twelve pencils. Fraction-free Bareiss/LDL computes every such inertia over the
integers; a zero pivot is rejected rather than interpreted numerically.

## 3. Forced-value absence certificates

**[LEMMA] Gaussian forced value.** Suppose a positive `2^9`-element Gaussian
spectrum has subset-product form

```text
lambda_0 { product_(j in S) u_j : S subset {1,...,9} },  u_j > 1.
```

If its three lowest distinct values are `lambda_0 < lambda_1 < lambda_2`, then
`lambda_1/lambda_0` and `lambda_2/lambda_0` are its two smallest one-particle
factors. Their two-element subset forces

```text
lambda_forced = lambda_1 lambda_2 / lambda_0
```

to occur in the spectrum. Positive rational enclosures `[l_i,h_i]` give the
outward enclosure `[l_1 l_2/h_0, h_1 h_2/l_0]`.

**[COMPUTATION] Exact endpoint data.** Float64 is used only to propose each
initial bracket. Both endpoints are accepted only after exact sector-inertia
counts establish the bracketing invariant. The following decimals are display
values for exact rational endpoints stored in the JSON.

| `t` | `lambda_0` | `lambda_1` | `lambda_2` | forced-value window | exact counts |
|---:|---:|---:|---:|---:|---:|
| `1/3` | `[4.175582455794e-4,4.183941980230e-4]` | `[7.595515344851e-4,7.610721581778e-4]` | `[1.398067167484e-3,1.400866100752e-3]` | `[2.538046816600e-3,2.553320878955e-3]` | `8,8` |
| `1/2` | `[3.527937268310e-6,3.535000205784e-6]` | `[2.284036968313e-5,2.288609614896e-5]` | `[2.772196507428e-5,2.777746450385e-5]` | `[1.791173674059e-4,1.801953025412e-4]` | `11,11` |

**[THEOREM] `t=1/3` certificate.** The three exact enclosures are disjoint. At
both ends of the outward forced-value window, the per-sector exact inertia
vectors are identical and their sums are `8=8`. Hence the window contains no
eigenvalue, contradicting the Gaussian forced-value lemma.

**[THEOREM] `t=1/2` certificate.** The three exact enclosures are disjoint. At
both ends of the outward forced-value window, the per-sector exact inertia
vectors are identical and their sums are `11=11`. Hence the window contains no
eigenvalue, contradicting the Gaussian forced-value lemma.

**[COMPUTATION] Spectrum bookkeeping.** At both couplings, `lambda_0` lies in
`A1` with negative spin-flip sign, `lambda_1` lies in `A1` with positive
spin-flip sign, and `lambda_2` has multiplicity two across the two refined
positive-spin-flip `E` rows. This explains why the third and fourth ordered
eigenvalues coincide without invalidating the strict inequalities
`lambda_0 < lambda_1 < lambda_2` used by the forced-value lemma.

## 4. Controls and independent dense bookkeeping

**[COMPUTATION] One-dimensional controls.** The independently certified open
chain with `n=9` gives forced-window counts `10,11` at both `t=1/3` and `t=1/2`.
Thus each positive-width predicted window is occupied, as required of the
integrable control. Differing inertia counts prove that some eigenvalue lies in
the window; they do not by themselves prove exact algebraic equality with
`lambda_1 lambda_2/lambda_0`.

**[LEMMA] Exact dyadic residual inertia witness.** Let `M=A-sI`. A floating
eigensolver proposes `Q` and a diagonal `D`, after which the script rounds `Q`
to the exact dyadic matrix `Q_Z/2^36` and `D` to integers. It then computes
`Q_Z^T Q_Z`, `Q_Z D Q_Z^T`, both Frobenius residual bounds, and the strict
spectral safety inequality entirely with Python integers. If

```text
||M-Q D Q^T||_F < (1-||Q^T Q-I||_F) min_i |D_ii|,
```

then `Q` is invertible, the model has the inertia of `D` by congruence, and the
homotopy to `M` cannot cross zero. Thus the resulting dense inertia count is
exact even though float64 supplied the discarded proposal.

**[COMPUTATION] Full 512-dimensional check.** At `t=1/2`, the coarse exact
integer-`A` window `[3000000000,3100000000]` contains the narrower predicted
window. The full dense dyadic witness gives counts `11,11`, equal to the sums of
all twelve Bareiss sector counts. Both stored strict safety-margin numerators
are positive. This independently certifies the spectrum bookkeeping without
using the sector decomposition to decide the dense counts.

**[UNRESOLVED] Direct dense Bareiss route.** A direct full `512 x 512`
fraction-free Bareiss LDL count at shift `3000000000` did not finish before the
1200-second wall; the recorded outcome is therefore only a lower bound of 1200
seconds for that route. The exact dyadic residual certificate replaces it for
the completed dense comparison; no claim of a completed full Bareiss count is
made.

**[COMPUTATION] Reproducibility.** The recorded sector runs took about 32.95
seconds at `t=1/3` and 33.99 seconds at `t=1/2`; the full dense residual witness
took about 34.18 seconds. On the same Apple M4 Pro run,
`.venv/bin/python tests/test_spectral_3x3.py` independently rebuilt the
`t=1/2` matrix and all projectors, reproduced exact counts `11,11`, validated
the artifact, and ended in `PASS` after 38.46 seconds. Timings are observations,
not mathematical inputs.

## 5. Scope

**[UNRESOLVED] Coupling and size scope.** The two finite certificates do not
prove non-Gaussianity on an interval of coupling, for every coupling, for every
layer size, or in the thermodynamic limit.

**[UNRESOLVED] Operator-class scope.** The theorem uses the same
positive-spectrum subset-product hypothesis as Theorem S. It does not exclude
embeddings in larger Gaussian systems, altered parity-projected constructions,
or every possible auxiliary-degree-of-freedom representation.

**[UNRESOLVED] Control equality scope.** Window occupancy in the chain controls
is a consistency certificate only. An exact equality certificate for the
predicted algebraic value was not constructed here.
