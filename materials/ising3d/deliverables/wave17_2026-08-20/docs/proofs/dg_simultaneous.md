# Simultaneous Dolan–Grady and tridiagonal deformation: an exact bounded no-go

**[THEOREM] Status.** In the precisely enumerated projective chart below, the simultaneous
Dolan–Grady variety and the more general simultaneous tridiagonal/Askey–Wilson variety have no
noncommuting points.  Their complete solution sets consist only of proportional generators
`B' = delta A'`, which commute and therefore do not form an integrable tridiagonal pair.

**[COMPUTATION] Provenance.** The exact calculation is implemented by
`experiments/e36_dg_simultaneous.py`, recorded in
`results/deformation/dg_simultaneous.json`, and locked by
`tests/test_dg_simultaneous.py`.  It uses Python integers, `Fraction`, and SymPy polynomials over
`QQ`; no floating-point arithmetic occurs.

## 1. Four mandatory wave-4 controls

**[COMPUTATION] Control normalization.** The frozen wave-4 script writes
`ad_B^3(A) = 16 lambda_frozen ad_B(A)`.  Thus its `lambda_frozen=1` is the direct constant
`lambda=16` used below.

| tag | exact control | result |
|---|---|---|
| [COMPUTATION] | cycle `C6` | 12 DG2 Pauli equations; `lambda_frozen=1`; zero residual terms; tridiagonal parameters `(beta,gamma,rho)=(2,0,16)` from 96 equations |
| [COMPUTATION] | path `P6` plus bond `(1,4)` | 14 DG2 Pauli equations; no floating constant; exactly 2 quartic obstruction terms; the 84-equation tridiagonal system is also inconsistent |
| [COMPUTATION] | torus `3x3` F3 point | `A_F3=sum_v X_v[1-(Z_(v+x)Z_(v-x)+Z_(v+y)Z_(v-y))/2]`; its complete DG2 residual at direct `lambda=16` has zero terms |
| [COMPUTATION] | DG1 witness at the F3 point | `Y0 Z1 X2 X3` has coefficient `3` in `ad_A_F3^3(B)` and coefficient `0` in `ad_A_F3(B)` |

**[COMPUTATION] Ordering guarantee.** `run_all()` executes and asserts all four controls before it
constructs the new ansatz.  A failed control terminates the experiment with a final `FAIL` line and
a nonzero exit status.

## 2. The smallest common symmetry-orbit basis

**[LEMMA] Basis enumeration.** Let

```text
O_X = sum_v X_v,
O_E = sum_<uv> Z_u Z_v,
O_C = sum_v X_v (Z_(v+x)Z_(v-x) + Z_(v+y)Z_(v-y)).
```

Each displayed operator is one orbit sum under translations semidirect `D4`.  On the `3x3` torus
the group order is `72`, and the three orbit sizes are respectively `9`, `18`, and `18`, with
representatives `X0`, `Z0 Z1`, and `X0 Z1 Z2`.  On the `4x4` torus the group order is `128`, and
the orbit sizes are `16`, `32`, and `32`.

**[LEMMA] Minimality proof.** The Pauli supports of `O_X`, `O_E`, and `O_C` are pairwise disjoint,
so the three orbit sums are linearly independent.  The successful F3 pair uses `O_X` and `O_C` in
its first generator and `O_E` in its second generator.  Therefore every common Pauli-orbit basis
containing that pair contains at least these three independent orbit sums; the displayed
three-dimensional basis is minimal.

**[COMPUTATION] Simultaneous chart.** Both generators independently range over this basis:

```text
A' = O_X + alpha O_E + chi O_C,
B' = delta O_X + O_E + eta O_C.                         (1)
```

The common basis gives `6` raw pair coefficients.  Independent generator scalings fix the
coefficient of `O_X` in `A'` and the coefficient of `O_E` in `B'` to one, leaving the four exact
variables `(alpha,chi,delta,eta)`.  The wave-4 F3 point is
`(alpha,chi,delta,eta)=(0,-1/2,0,0)`, so (1) strictly contains it and also deforms `B`.

**[COMPUTATION] Exact scope.** The search covers every rational or complex point in this
four-dimensional chart, with arbitrary direct Dolan–Grady constants or arbitrary tridiagonal
parameters.  It does not cover the projective hyperplanes where the fixed `O_X` or `O_E`
coefficient vanishes, non-invariant coefficients, Pauli orbits beyond `O_X,O_E,O_C`, or generator
pairs outside this common local basis.

## 3. Simultaneous Dolan–Grady relations

**[COMPUTATION] Equations.** The imposed relations are

```text
DG1: ad_A'^3(B') = mu     ad_A'(B'),
DG2: ad_B'^3(A') = lambda ad_B'(A').                    (2)
```

Here `mu` and `lambda` are direct constants, so the conventional Dolan–Grady point is `16`.
Define the three coefficient-vector cross products

```text
p = alpha delta - 1,
q = chi delta - eta,
r = alpha eta - chi.                                    (3)
```

The two coefficient vectors `(1,alpha,chi)` and `(delta,1,eta)` are proportional exactly when
`p=q=r=0`.

**[COMPUTATION] Orbit reduction.** Exact Pauli expansion and translation/`D4` compression give:

| relation | full Pauli equations | orbit equations | distinct polynomials up to rational scale |
|---|---:|---:|---:|
| [COMPUTATION] DG1 | 3384 | 57 | 30 |
| [COMPUTATION] DG2 | 3384 | 57 | 30 |
| [COMPUTATION] both | 6768 row slots | 114 row slots | 60 |

**[LEMMA] Low-support DG certificate.** Three coefficient rows decide the noncommuting locus before
any Gröbner calculation:

```text
DG2, Y0 Z1 Z2 Z3 X6:       24 q.
DG1, X0 Y1 X2 Z3, q=0:      8 chi p.
DG2, Y0 Z1 Z2 Z3, chi=eta=0: 48 p.                       (4)
```

**[LEMMA] Consequence of (4).** The first row forces `q=0`.  If `p=0`, then `q=0` also gives
`r=0`, hence the generators are proportional.  On a putative nonproportional point, therefore,
`p != 0`; the second row then forces `chi=0`, and `q=0` forces `eta=0`.  The third row finally
forces `p=0`, a contradiction.  Thus (2) has no nonproportional point in the chart.

**[COMPUTATION] Linear elimination before Gröbner.** Each relation is linear in its floating
constant.  The script eliminates `mu` and `lambda` first by all exact `2x2` minors of the orbit-row
pairs `(third commutator, first commutator)`.  This produces `80` distinct DG1 polynomials and `80`
distinct DG2 polynomials, `160` in their union.

**[COMPUTATION] Gröbner certificate.** Only after that linear elimination, the grevlex basis over
`QQ[alpha,chi,delta,eta]` is

```text
p^2, p q, q^2, p r, q r, r^2.                            (5)
```

Saturating (5) separately by `p`, `q`, and `r` gives the unit basis `[1]` on all three charts.
These three charts cover the nonproportional locus, so this is an exact empty-variety certificate
over every characteristic-zero field, including `Q` and `C`.

**[THEOREM] DG decision.** The complete Dolan–Grady solution variety in (1) is

```text
alpha delta = 1,    eta = chi delta,
B' = delta A'.                                           (6)
```

The constants `mu` and `lambda` are arbitrary on (6), because both sides of both relations vanish.
This is a commuting degeneracy, not a noncommuting Dolan–Grady pair.

## 4. The strictly more general tridiagonal relations

**[COMPUTATION] Equations.** The same deformed pair is tested against

```text
TD1: [A', A'^2 B' - beta A'B'A' + B'A'^2
          - gamma(A'B'+B'A') - rho B'] = 0,
TD2: [B', B'^2 A' - beta B'A'B' + A'B'^2
          - gamma*(B'A'+A'B') - rho* A'] = 0.             (7)
```

The shared `beta` and independent `(gamma,rho,gamma*,rho*)` make (7) a five-parameter family;
Dolan–Grady is the point `(2,0,16,0,16)`.

**[COMPUTATION] Orbit reduction.** Exact expansion gives:

| relation | full Pauli equations | orbit equations | distinct polynomials up to rational scale |
|---|---:|---:|---:|
| [COMPUTATION] TD1 | 16344 | 263 | 134 |
| [COMPUTATION] TD2 | 16344 | 263 | 134 |
| [COMPUTATION] both | 32688 row slots | 526 row slots | 266 |

**[LEMMA] First tridiagonal stage.** On any nonproportional point, three exact rows force
`beta=2`:

```text
Y0 Z1 X2 X3 X4:       4(beta-2) q,
Y0 Z1 Z2 Z3 X4 X5:    4(beta-2) r,
Y0 Z1 X2 X3, q=r=0:   4(beta-2) p.                       (8)
```

Indeed, if `q` is nonzero the first row decides `beta`; if `q=0` and `r` is nonzero the second does.
If `q=r=0` but the point is nonproportional, then `p != 0`; the identities
`eta=chi delta` and `r=chi p` give `chi=eta=0`, so the third row decides `beta`.

**[LEMMA] Second tridiagonal stage.** At `beta=2`, the rows

```text
TD2, Y0 X1 Z2 Z4 Z6:       24 q,
TD1, X0 Y1 X2 Z3, q=0:      8 chi p                         (9)
```

force `q=0`, and on a nonproportional point force `p != 0` and then `chi=eta=0`.

**[LEMMA] Terminal tridiagonal certificate.** With `beta=2` and `chi=eta=0`, two TD2 rows are

```text
Y1 Z3:          -8 gamma* p,
Y0 Z1 Z2 Z3:   -4(gamma* - 12) p.                         (10)
```

On the saturated chart `p != 0`, the support-two row gives `gamma*=0`, while the support-four row
gives `gamma*=12`.  Eliminating the linear variable `gamma*` first leaves the exact nonzero
constant `48`; the ensuing Gröbner basis is `[1]`.

**[THEOREM] Tridiagonal decision.** Equations (8)–(10) exclude every nonproportional point of (1).
Conversely, every proportional pair (6) makes all outer commutators in (7) vanish for arbitrary
five parameters.  Hence the complete tridiagonal variety in the chart is again precisely the
commuting component (6), and no Askey–Wilson deformation survives.

## 5. Size independence and interpretation

**[COMPUTATION] Two-size exact check.** The single fixed rational tuple

```text
(alpha,chi,delta,eta) = (1/2,-1/2,2,-1)
```

was evaluated without refitting on both `3x3` and `4x4` torus layers.  It obeys `B'=2A'`; its
commutator and both DG residuals (at arbitrarily selected direct constants `7` and `11`) have zero
terms on both sizes.  Its two tridiagonal residuals at the arbitrarily selected parameters
`(beta,gamma,rho)=(3,5,7)` and `(beta,gamma*,rho*)=(3,11,13)` also have zero terms on both sizes.

**[COMPUTATION] Nonintegrable cancellation check.** The fixed F3 coefficients
`(0,-1/2,0,0)` cancel DG2 at direct `lambda=16` on both `3x3` and `4x4`, but no DG1 constant exists
on either size.  Thus even a size-independent cancellation of one relation is not integrability.

**[THEOREM] Size conclusion.** Any size-independent noncommuting coefficient tuple that worked on
both sizes would in particular give a noncommuting point of the certified `3x3` variety, which is
empty.  The two-size solution check therefore confirms only the explicitly classified commuting
degeneracy and makes no integrability claim.

## 6. Exact scope of the no-go

**[THEOREM] Covered statement.** The no-go is complete for the four-dimensional chart (1) on the
minimal three-orbit common basis and includes arbitrary relation constants, every algebraic branch,
and both Dolan–Grady and general tridiagonal relations.

**[COMPUTATION] Uncovered directions.** No conclusion is drawn about the two omitted projective
hyperplanes, larger translation/point-group orbit bases, coefficients that break those symmetries,
non-Pauli generators, or alternative factorizations of the transfer operator.  The result is a
bounded exact obstruction, not a theorem that the three-dimensional Ising model is unsolvable.
