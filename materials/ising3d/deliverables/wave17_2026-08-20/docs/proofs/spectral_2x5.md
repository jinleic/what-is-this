# Exact spectral non-Gaussianity of the open 2x5 layer

**Status tags.** **[THEOREM]** proved here; **[LEMMA]** proved here;
**[COMPUTATION]** finite exact computation; **[UNRESOLVED]** explicitly outside the result.

Script: `experiments/e67_spectral_2x5.py`  
Artifact: `results/spectral/spectral_2x5.json`  
Independent test: `tests/test_spectral_2x5.py`

## 1. Statement and scope

**[THEOREM S25].** Let `R` be the symmetrized transfer operator of the ten-spin open
`2x5` Ising layer, with the conventions below. At each of

- `t = tanh(K*/2) = 1/3`, hence `exp(2K) = 5/3`, and
- `t = tanh(K*/2) = 1/2`, hence `exp(2K) = 5/4`,

`R` does **not** have the positive spectrum of a Gaussian operator on ten fermionic modes.
Equivalently, its 1024 eigenvalues are not one positive ten-generator subset-product multiset.

This is a finite-size, two-coupling theorem. It does not solve the three-dimensional Ising
model and does not imply a coupling interval, an all-size result, an embedding theorem, or a
thermodynamic-limit statement.

The two rational values were fixed as simple carry-over and stress-test points before the
calculation. No critical-coupling benchmark was used to tune, select, or fit either value.

## 2. Exact operator

Number sites in row-major order. The open `2x5` graph has 10 sites, 13 bonds, and degree-3
vertices. Put

```
P_t[k,l] = t^hamming(k,l),
q = exp(2K) = (1+t^2)/(2t),
b_k = sum_(ij in E) z_i z_j.
```

All `b_k` have the parity of the bond count, here odd. With `eps = 1`, define

```
R = P_t diag(q^((b_k-eps)/2)) P_t.
```

**[LEMMA] (rational symmetric representative).** Up to a common positive scalar, `R` is
the symmetrized layer transfer operator and has the same eigenvalue ratios. For rational `t`,
a positive integer multiple `A = scale R` is an exact symmetric integer matrix.

The proof is the same elementary factorization used by the earlier spectral certificates:
`exp(K*A_spin/2)` is a positive scalar times `P_t`, while the bond diagonal contributes
`q^((b_k-eps)/2)` after removing one common scalar. The script clears all denominators before
any spectral decision. It obtains

| `t` | `q` | integer scale | exponent range |
|---|---:|---:|---:|
| `1/3` | `5/3` | `198583267838203125` | `[-7,6]` |
| `1/2` | `5/4` | `335544320000000` | `[-7,6]` |

The exact inertia below zero is 0 at both couplings, independently confirming positive
definiteness in the same sector bookkeeping used below.

## 3. Complete exact symmetry split

The rectangle admits commuting row reflection `r`, column reflection `c`, and global spin flip
`f`. They generate `C2 x C2 x C2`. For signs `rho,kappa,phi` in `{+1,-1}`, the character
projector is

```
P_(rho,kappa,phi) = (1/8) sum_(a,b,d in {0,1})
                    rho^a kappa^b phi^d r^a c^b f^d.
```

**[COMPUTATION].** On all 1024 computational-basis columns, integer projector numerators
satisfy exactly:

1. self-adjointness;
2. `N_chi^2 = 8 N_chi`;
3. `N_chi N_psi = 0` for all 28 unequal pairs;
4. `sum_chi N_chi = 8 I`;
5. integral traces summing to 1024; and
6. commutation with `A` under all eight group actions.

The complete sector dimensions are

| row sign | column sign | flip sign | dimension |
|---:|---:|---:|---:|
| `+` | `+` | `+` | 152 |
| `+` | `+` | `-` | 136 |
| `+` | `-` | `+` | 120 |
| `+` | `-` | `-` | 120 |
| `-` | `+` | `+` | 120 |
| `-` | `+` | `-` | 136 |
| `-` | `-` | `+` | 120 |
| `-` | `-` | `-` | 120 |

Their sum is 1024. Thus every full-spectrum inertia below is the exact sum over a complete,
mutually orthogonal projector family; no unexamined block is omitted.

For each character, exact independent integer columns of the projector numerator form a basis
matrix `C_chi`. The restricted generalized pencil is

```
B_chi - sigma G_chi,
B_chi = C_chi^T A C_chi,
G_chi = C_chi^T C_chi.
```

Because distinct orbit vectors have disjoint support, every `G_chi` is positive diagonal. Hence
the number of generalized eigenvalues below `sigma` is the negative inertia of this pencil.

## 4. Exact integer-congruence inertia certificate

The 1024-dimensional dense check is optional and was not used. The sector checks are internally
complete and are certified as follows.

**[LEMMA] (relative congruence certificate).** Let `M = d B - n G` for the rational shift
`sigma=n/d`, `d>0`. Let an arbitrary integer square matrix `Q` be proposed and compute exactly

```
C = Q^T M Q = D + E,
```

where `D` is diagonal and `E` has zero diagonal. Suppose every `D_ii` is nonzero. Put
`S_ii=sqrt(abs(D_ii))`, `J_ii=sign(D_ii)`, and `F=S^(-1) E S^(-1)`. If `||F||_2 < 1`, then `Q`
is invertible and

```
inertia(M) = inertia(C) = inertia(D).
```

*Proof.* We have `C=S(J+F)S`. If `J+sF` were singular for some `0<=s<=1`, a nonzero `x`
would obey `x=-s J F x`, impossible because `J` is orthogonal and `s||F||_2<1`. Therefore
`J+sF` stays nonsingular from `s=0` to `s=1`, so its inertia is constant. Thus `C` is
nonsingular. If `Q` were singular then `Qx=0` for some nonzero `x`, forcing `Cx=0`, a
contradiction. Sylvester congruence now gives the result. []

The implementation proves the norm hypothesis without floating point or square roots. Set
`k_i=floor(log2(abs(D_ii)))`. Then

```
||F||_2^2 <= ||F||_F^2
            = sum_(i != j) E_ij^2/(abs(D_ii) abs(D_jj))
            <= sum_(i != j) E_ij^2/2^(k_i+k_j).
```

The final dyadic sum is accumulated as one integer numerator over one power-of-two denominator,
and the stored strict integer margin is positive for every sector at every shift. Float64 only
proposes `Q` (the numerator of a 36-bit dyadic approximation to generalized eigenvectors) and the
initial rational brackets. All matrices `B,G,Q,C`, all norm bounds, all margins, and every inertia
count used in the theorem are recomputed with exact Python integers. A poor floating proposal can
only make the strict inequality fail; it cannot make an incorrect count pass.

## 5. The Gaussian forced-value criterion

**[LEMMA] (multiplicative Gaussian test).** If a positive `2^m`-element Gaussian spectrum has
three distinct lowest values

```
lambda0 < lambda1 < lambda2,
```

then `lambda1*lambda2/lambda0` is also an eigenvalue.

Indeed, after dividing by the lowest value, such a spectrum is the multiset of products of
subsets of `m` factors `u_i>=1`. Distinctness forces the next two values to be the two smallest
single factors. Their product is another allowed subset product. This is a necessary condition
on the spectrum alone and is unaffected by a change of spin or Majorana basis.

If exact brackets are `lambda_i in [l_i,h_i]`, positivity gives the outward interval

```
I_forced = [l1*l2/h0, h1*h2/l0].
```

Equal exact inertia counts at both nonsingular endpoints certify that this entire closed interval
contains no eigenvalue.

## 6. Exact results

Every one of the three brackets has relative width `2/999`. Endpoint inertias are respectively
`(0,1)`, `(1,2)`, and `(2,3)`, so each bracket isolates exactly one of the first three
eigenvalues. The brackets are pairwise disjoint.

### `t=1/3`

| value | certified interval for `R` |
|---|---:|
| `lambda0` | `[1.389247757450e-4, 1.392029034242e-4]` |
| `lambda1` | `[2.600732781010e-4, 2.605939453244e-4]` |
| `lambda2` | `[3.672147350310e-4, 3.679498996657e-4]` |
| forced value | `[6.860685916549e-4, 6.901973785554e-4]` |

The exact sector-summed inertias at the forced endpoints are **8 and 8**. In fact, the
per-sector vectors agree entry by entry:

```
(2,1,0,1,0,2,0,2) at both endpoints.
```

Thus the forced interval is empty.

### `t=1/2`

| value | certified interval for `R` |
|---|---:|
| `lambda0` | `[7.838678270298e-7, 7.854371319888e-7]` |
| `lambda1` | `[5.137126632228e-6, 5.147411170030e-6]` |
| `lambda2` | `[5.691432990394e-6, 5.702827250635e-6]` |
| forced value | `[3.722463682925e-5, 3.744865611068e-5]` |

The exact sector-summed inertias at the forced endpoints are **12 and 12**, again with identical
per-sector vectors:

```
(2,3,0,2,0,3,0,2) at both endpoints.
```

Thus this forced interval is also empty. The multiplicative lemma proves Theorem S25 at both
couplings.

## 7. Chain controls

The same exact integer operator construction and the complete
`C2(chain reflection) x C2(spin flip)` split were applied to the open ten-site chain. Its four
sector dimensions are `272,256,240,256`, summing to 1024. The forced-window counts are

| control | lower count | upper count |
|---|---:|---:|
| chain `n=10`, `t=1/3` | 11 | 12 |
| chain `n=10`, `t=1/2` | 11 | 12 |

Both windows are occupied, as expected for the free-fermionic control. This is only a consistency
check: a nonzero-width occupied interval proves that some eigenvalue lies in the interval, not
that it equals the forced expression exactly.

## 8. Resources and independent verification

**[COMPUTATION].** On the recorded Darwin/Apple-arm run, the four cases took approximately
`31.50 s`, `32.77 s`, `80.92 s`, and `70.08 s`; total wall time was `215.58 s`. Peak process RSS
reported by `getrusage` was `149258240` bytes. The artifact records matrix construction,
projector verification, sector construction, congruence preparation, and exact-count time for
each case, plus every per-sector strict integer margin.

Run:

```bash
.venv/bin/python experiments/e67_spectral_2x5.py
.venv/bin/python tests/test_spectral_2x5.py
```

The standalone test does not import `e67` or `e54`. It independently rebuilds the `t=1/3`
1024-by-1024 integer matrix, all eight projectors, every exact projector identity and commutator,
the complete sector bases, and independent integer congruence witnesses at zero, all six
eigenvalue endpoints, and both forced-window endpoints. It reproduces forced counts `(8,8)` and
the stored per-sector counts.

## 9. Escape routes left open

**[UNRESOLVED].** This theorem excludes one full positive-spectrum ten-mode Gaussian
subset-product representation of these two finite transfer operators. It does not exclude:

- parity-projected or paired Gaussian descriptions that are not the full 1024-value subset-product
  spectrum;
- realization as a restriction of a larger Gaussian operator;
- isolated algebraic coincidences at other couplings;
- Gaussian notions with non-positive or complex spectrum, outside the lemma's hypotheses;
- an unrelated exact method for the infinite three-dimensional model.

No finite computation here is promoted to all layer sizes or to an interval of couplings.
