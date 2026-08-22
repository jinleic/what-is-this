# Spectral non-Gaussianity at additional rational couplings and at `2 x 4`

**Status tags.** **[LEMMA]** proved below; **[COMPUTATION]** exact finite computation with a
reproducible certificate; **[UNRESOLVED]** not decided here.

Artifacts:

- experiment: `experiments/e39_spectral_couplings.py`;
- exact data: `results/spectral/couplings.json`;
- standalone independent check: `tests/test_spectral_couplings.py`.

This note extends the exact finite-operator calculation behind Theorem S without changing
`proofs/spectral_gaussianity.md`. The multiplicative obstruction used here is exactly the one proved
there: if the three lowest eigenvalues of a positive-spectrum Gaussian operator are distinct, then
`lambda_1 lambda_2 / lambda_0` must be an eigenvalue.

## 1. Exact operator and integer clearing

For rational `t = tanh(K*/2)`, put

```text
q = exp(2K) = (1+t^2)/(2t),
P_t[k,l] = t^hamming(k,l),
R = P_t diag(q^((b_k-eps)/2)) P_t.
```

Here `b_k` is the in-layer Ising bond energy and `eps` is its common parity. As in Theorem S, the
symmetrized transfer operator is a positive scalar times `R`; hence the scalar cancels from all
ratios in the obstruction.

**[LEMMA] (single integer clearing).** If `t=a/b` and `q=c/d` are reduced positive rationals, then a
positive integer `D` can be chosen so that `A=D R` is an integer symmetric positive-definite matrix.
The inertia of `R-sigma I` equals that of `A-D sigma I`.

*Proof.* Write `P_t=P_num/b^n`, where `P_num` is the `n`-fold Kronecker power of
`[[b,a],[a,b]]`. The finitely many integer exponents `e_k=(b_k-eps)/2` have a minimum and maximum;
`d^max(max e_k,0) c^max(-min e_k,0)` clears every denominator in `q^e_k`. Multiplication by this
factor and by `b^(2n)` therefore makes `R` integral. Both clearing factors are positive, so
Sylvester inertia is unchanged after scaling the shifted matrix by them. []

The implementation forms this exact product in `O(n 4^n)` integer operations by applying the local
two-by-two Kronecker factors successively. No floating-point matrix entry participates in a
certificate.

## 2. Exact symmetry reduction for the 256-dimensional cases

The open rectangle has commuting row reflection, column reflection, and global spin complement.
The open chain has reflection and global spin complement. Each operation is an exact permutation
that commutes with `A`.

**[LEMMA] (character-sector inertia).** Let the columns of `C_chi` be nonzero, unnormalised integer
orbit vectors for one character of these commuting involutions. Put

```text
B_chi = C_chi^T A C_chi,       G_chi = C_chi^T C_chi.
```

The orbit supports are disjoint, so `G_chi` is positive diagonal. The number of eigenvalues of the
restriction of `A` below a rational shift `s` equals the negative inertia of
`B_chi-s G_chi`; summing over all characters gives the full count.

*Proof.* Normalised orbit vectors give the symmetric restriction
`G_chi^(-1/2) B_chi G_chi^(-1/2)`. Subtracting `s I` is congruent to
`B_chi-s G_chi`, so Sylvester's law gives equal inertia. The character sectors are mutually
orthogonal and their recorded dimensions sum to the full dimension. []

For `2 x 4`, the eight sector dimensions are
`44,32,28,32,28,32,28,32`, summing to 256. For the chain control they are
`72,64,56,64`, also summing to 256. The JSON records exact checks that every proposed permutation
commutes with the constructed matrix.

Each sector count uses fraction-free symmetric Bareiss/LDL elimination. A zero pivot is never
silently accepted: the evaluated shift is reported as degenerate, and bisection tries a different
interior rational. At a final endpoint, any movement is outward only.

Float64 diagonalisation is used solely to propose small initial brackets. It proves nothing:
every accepted lower endpoint and upper endpoint has its own exact rational inertia count. The
three enclosures then come from exact rational bisection to relative target `10^-7` for every layer
case. The chain consistency control uses target `10^-4`; its endpoints and counts are still exact.

## 3. Exact finite certificates

The table gives display decimals only. Full rational endpoints are in
`results/spectral/couplings.json`.

| finite operator | `t` | `exp(2K)` | dimension | certified `lambda_0` enclosure | certified `lambda_1` enclosure | certified `lambda_2` enclosure | rigorous relative-gap lower bounds | forced-value window | inertia counts |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| open `2 x 3` layer | `1/5` | `13/5` | 64 | `[9.064793345888e-3, 9.064794037477e-3]` | `[9.431479674461e-3, 9.431480394026e-3]` | `[3.709626369509e-2, 3.709626652531e-2]` | `0.0388789`, `0.745757` | `[3.859686779337e-2, 3.859687662750e-2]` | `4,4` |
| open `2 x 3` layer | `2/5` | `29/20` | 64 | `[1.583836333295e-3, 1.583836454132e-3]` | `[5.434439505190e-3, 5.434439919805e-3]` | `[7.651299307106e-3, 7.651299890854e-3]` | `0.708556`, `0.289736` | `[2.625304090715e-2, 2.625304691599e-2]` | `8,8` |
| open `2 x 3` layer | `1/2` | `5/4` | 64 | `[2.087711792756e-4, 2.087711952036e-4]` | `[1.427652427451e-3, 1.427652536372e-3]` | `[1.716134906883e-3, 1.716135037814e-3]` | `0.853766`, `0.168100` | `[1.173554696209e-2, 1.173554964815e-2]` | `8,8` |
| open `2 x 4` layer | `1/3` | `5/3` | 256 | `[1.025537205153e-3, 1.025537283395e-3]` | `[2.017393635357e-3, 2.017393789272e-3]` | `[3.079201883056e-3, 3.079202117980e-3]` | `0.491652`, `0.344832` | `[6.057275909359e-3, 6.057277295759e-3]` | `8,8` |
| open chain `n=8` control | `1/3` | `5/3` | 256 | `[9.121474479417e-4, 9.122187094610e-4]` | `[2.299613740292e-3, 2.299793397615e-3]` | `[2.619924520044e-3, 2.620129217639e-3]` | `0.603317`, `0.122191` | `[6.604572305233e-3, 6.606120413124e-3]` | `9,10` |

For all five rows, the displayed intervals are pairwise disjoint in increasing order. The relative
gaps in the table are exact lower bounds computed as `(l_1-h_0)/l_1` and `(l_2-h_1)/l_2`, not
floating estimates.

**[COMPUTATION] (additional `2 x 3` couplings).** At each of `t=1/5`, `t=2/5`, and `t=1/2`, the
three lowest eigenvalues of the exact rational 64-dimensional open-layer operator are distinct.
The exact inertia counts at both outward-rounded ends of the forced-value window agree (`4=4`,
`8=8`, and `8=8`, respectively). Thus each window contains no eigenvalue, whereas a Gaussian
operator on six modes would be forced to have one. Consequently each of these three finite transfer
operators is not fermionic Gaussian on six modes, in any choice of Majorana modes satisfying the
positive-spectrum hypothesis of Theorem S.

**[COMPUTATION] (`2 x 4`).** At `t=1/3`, the three lowest eigenvalues of the exact rational
256-dimensional open `2 x 4` layer operator are distinct. Both exact outward-window inertia counts
are 8. Therefore the forced window contains no eigenvalue, and this finite operator is not
fermionic Gaussian on eight modes under the same hypothesis. This case completed; it is not an
`UNRESOLVED-RESOURCE` result.

**[COMPUTATION] (chain control, consistency only).** For the open chain `n=8` at `t=1/3`, the exact
counts differ, `9` versus `10`. This proves only that **some** eigenvalue lies in the
positive-width forced-value window. It does not prove exact equality with
`lambda_1 lambda_2/lambda_0`, and it is not used as an exact proof of Gaussianity.

## 4. Resource record and independent verification

On the recorded Apple M4 Pro run, the `2 x 4` row used 66 exact inertia shifts and completed in
about 7.28 seconds; its slowest full-spectrum count took about 0.38 seconds. The `n=8` chain control
used 36 shifts and completed in about 38.8 seconds. These timings are observations, not mathematical
claims. The configured hard walls were 3000 seconds per inertia count and 14400 seconds total for
the `2 x 4` case, so neither wall was approached.

`tests/test_spectral_couplings.py` does not import the experiment. It rebuilds the `2 x 3`, `t=1/5`
matrix by the direct dense rational sum, clears denominators independently, repeats full-matrix
exact bisection without symmetry sectors, and obtains the decisive counts `4,4`. It then recomputes
those two counts by the slower Fraction-valued LDL algorithm used in `e38`, again obtaining `4,4`.
The same test parses every exact rational stored in the JSON and checks interval arithmetic,
distinctness, outward containment, count/verdict consistency, and sector-dimension sums for all
rows.

## 5. Scope and limitations

- Each absence result above is a **[COMPUTATION]-grade exact theorem about one named finite rational
  operator**. There is no numerical tolerance in the spectral claim.
- These computations do not prove non-Gaussianity for an interval of coupling, for every coupling,
  for all layer sizes, or in the thermodynamic limit. In particular, three additional rational
  points do not establish a coupling-independent theorem.
- The result inherits the precise operator class and positive-real-spectrum hypothesis of the
  spectral-Gaussianity lemmas in `proofs/spectral_gaussianity.md`. It does not exclude embeddings
  into larger Gaussian systems or every parity-projected construction discussed there.
- A differing inertia count in the control direction is consistency evidence only. Exact presence
  of the predicted algebraic value would require an equality certificate and is **[UNRESOLVED]**
  here.
