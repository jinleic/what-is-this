# Exact characteristic-zero pair-product ranks for the five-site tree

**[THEOREM] Status.** At each of the four named physical points

\[
 (t,w)\in\{(1/3,5/3),(1/4,17/8),(2/5,29/20),(1/5,13/5)\},
 \qquad w=\frac{1+t^2}{2t},
\]

the five-site tree core has

\[
 r(T)=417,\qquad r_{\rm all}(T)=445
\]

exactly in characteristic zero.  Thus the characteristic-zero hypothesis CZ in
`proofs/allsize_gaussian.md` is discharged by the named witness
`(t,w)=(1/3,5/3)`; the other three rows are independent exact witnesses.

**[THEOREM] Restored U2 scope.** The previously withdrawn composite counts are now
exact: `r(T \uplus P_2)=3893` at all four points above and
`r(T \uplus P_4)=35597` at `(1/3,5/3)`.  These are precisely the named-point rows
computed in `proofs/allsize_gaussian.md` §10; no value at an unlisted coupling is
included in this statement.

---

## 1. [LEMMA] Exact definitions

**[LEMMA] The graph and operator.** Number the vertices of `T` by `0,...,4` and
use

\[
 E(T)=\{01,02,03,34\}.
\]

Vertex `0` has degree three, and its branches have lengths `1,1,2`.  In the
computational spin basis, put

\[
 P(t)=\bigotimes_{v=0}^4\begin{pmatrix}1&t\\t&1\end{pmatrix},\qquad
 D_T(w)_{\sigma,\sigma}=w^{a_T(\sigma)},
\]

where `a_T(sigma)` is the number of aligned edges of `T`, and define

\[
 R_T(t,w)=P(t)D_T(w)P(t).
\]

This is the `R_G` operator of `proofs/allsize_gaussian.md` §2.  A nonzero scalar
multiple has the same pair-product collision counts.

**[LEMMA] Slot counts.** If
`Lambda_T=(lambda_1,...,lambda_32)` is the spectrum with algebraic-multiplicity
slots retained, then

\[
\begin{aligned}
 r(T)&=\#\{\lambda_i\lambda_j:i<j\},\\
 r_{\rm all}(T)&=\#\{\lambda_i\lambda_j:i\le j\}.
\end{aligned}
\]

The braces count distinct displayed algebraic values after the allowed unordered
slot pairs are formed.  In particular, two equal eigenvalues remain two slots,
so their square occurs in `r(T)`.

**[COMPUTATION] Withdrawn values being repaired.** The withdrawn U2 paragraph in
`proofs/allsize_gaussian.md` asserted `r(T)=417`, `r_all(T)=445`, and then—under
INJ—`r(T \uplus P_2)=3893` and `r(T \uplus P_4)=35597`.  The old finite-field
calculation proved only the lower bounds `r(T)>=417`, `r_all(T)>=445`; inserting
those lower bounds into Lemma 5' did not give an upper bound.  The exact
characteristic-zero computation below supplies the missing direction.

---

## 2. [LEMMA] A factored characteristic-zero counting identity

**[LEMMA] Integral scaling.** Write `t=a/b` and `w=c/d` in lowest terms.  Then

\[
 M=b^{10}d^4R_T(t,w)
\]

is the exact integer matrix obtained by replacing each local `P` factor by
`[[b,a],[a,b]]` and each aligned/non-aligned edge factor by `c/d`.  Since the
scale is nonzero at every named point, `M` and `R_T` have identical equality
relations among pair products.

**[LEMMA] Root-set components.** Factor exactly over `Q`

\[
 \chi_M(x)=\prod_a f_a(x)^{e_a}
\]

with the `f_a` distinct, monic, irreducible and separable.  Define

\[
\begin{aligned}
 S_a(z)&=\operatorname{Res}_x(f_a(x),z-x^2),\\
 X_{ab}(z)&=\operatorname{Res}_x\!\left(f_a(x),
 x^{\deg f_b}f_b(z/x)\right),\\
 W_a(z)^2&=X_{aa}(z)/S_a(z),
\end{aligned}
\]

where `W_a` is the unique monic exact polynomial square root.  The roots of
`S_a`, `X_ab`, and `W_a` are respectively `alpha^2`, `alpha beta`, and
`alpha_i alpha_j` for `i<j`, with `alpha_i` roots of `f_a` and `beta` roots of
`f_b`.

**[LEMMA] Exact two-sided count.** Let `L_r` be the squarefree lcm of every `W_a`,
every `X_ab` with `a<b`, and `S_a` exactly when `e_a>=2`.  Let `L_all` add every
remaining `S_a`.  Then

\[
 r(T)=\deg L_r,\qquad r_{\rm all}(T)=\deg L_{\rm all}.
\]

For the **upper-bound direction**, every allowed slot-pair product is a root of
the relevant lcm, hence `r(T)<=deg L_r` and
`r_all(T)<=deg L_all`.  For the **lower-bound direction**, every resultant root
is an actual product of spectral roots, and the squarefree lcm has exactly its
degree in distinct roots, hence `r(T)>=deg L_r` and
`r_all(T)>=deg L_all`.  This proves equality without reducing modulo a prime.

---

## 3. [COMPUTATION] Exact characteristic-zero certificate

**[COMPUTATION] Spectral factorization.** At every named point, exact SymPy
factorization over `Q` gives the degree/exponent pattern

\[
 (1,2),(1,2),(3,1),(3,1),(11,1),(11,1).
\]

Thus the 32 slots contain 30 distinct eigenvalue values.  The artifact stores all
factor coefficients, irreducibility flags, matrix hashes, component-polynomial
hashes, squarefree degrees, overlap degrees, and cumulative lcm degrees.

**[COMPUTATION] Exact sandwiches.** The resultant/lcm calculation gives the same
closed sandwiches at all four named points:

| `t` | `w` | exact factor pattern | characteristic-zero sandwich for `r` | characteristic-zero sandwich for `r_all` |
|---:|---:|---|---:|---:|
| `1/3` | `5/3` | `1^2,1^2,3,3,11,11` | `417 <= r <= 417` | `445 <= r_all <= 445` |
| `1/4` | `17/8` | `1^2,1^2,3,3,11,11` | `417 <= r <= 417` | `445 <= r_all <= 445` |
| `2/5` | `29/20` | `1^2,1^2,3,3,11,11` | `417 <= r <= 417` | `445 <= r_all <= 445` |
| `1/5` | `13/5` | `1^2,1^2,3,3,11,11` | `417 <= r <= 417` | `445 <= r_all <= 445` |

**[COMPUTATION] Equivalent gcd degrees.** Since the exterior and symmetric slot
polynomials have degrees `C(32,2)=496` and `C(33,2)=528`, the exact values are

\[
 \deg\gcd_{\mathbb Q}(C_2,C_2')=496-417=79,
 \qquad
 \deg\gcd_{\mathbb Q}(C_{\rm all},C_{\rm all}')=528-445=83.
\]

These are characteristic-zero equalities, not modular one-sided bounds.

---

## 4. [THEOREM] U2 restored by exact two-sided squeezes

**[THEOREM] The `m=2` rows.** At each of the four named points, the existing exact
`F_1000003` calculation gives the one-sided characteristic-zero lower bounds

\[
 r(T\uplus P_2)\ge3893,\qquad r_{\rm all}(T\uplus P_2)\ge4005.
\]

The now-exact core counts give the unconditional Lemma 5' upper bounds

\[
\begin{aligned}
 r(T\uplus P_2)
 &\le(3^2-2^2)445+2^2\,417=3893,\\
 r_{\rm all}(T\uplus P_2)&\le3^2\,445=4005.
\end{aligned}
\]

Therefore `r(T \uplus P_2)=3893` and `r_all(T \uplus P_2)=4005` exactly.  Equality
in Lemma 5' certifies INJ(2), rather than assuming it.

**[THEOREM] The `m=4` row.** At `(t,w)=(1/3,5/3)`, the existing exact
`F_1000003` pair polynomial gives `r(T \uplus P_4)>=35597`; its coefficient hash
is `29e9f1fd1f6994b2df69191c4a9ac0715d5282383ea0f325aa6c212da0085ca3`.
The exact-core Lemma 5' bound has the opposite direction:

\[
 r(T\uplus P_4)\le(3^4-2^4)445+2^4\,417=35597.
\]

Hence `r(T \uplus P_4)=35597` exactly, and equality certifies INJ(4).

**[THEOREM] Exact excess identity at the certified rows.** For `m=2,4`,

\[
\begin{aligned}
 r(T\uplus P_m)
  &=(3^m-2^m)445+2^m417,\\
 r(T\uplus P_m)-(3^{5+m}-2^{5+m})
  &=(445-3^5)3^m-((445-417)-2^5)2^m\\
  &=202\,3^m+4\,2^m.
\end{aligned}
\]

The exact excesses are `1834` at `m=2` and `16426` at `m=4`.

---

## 5. [THEOREM] Consequence for F2 and the one remaining hypothesis

**[THEOREM] Unconditional `m=2,4` F2 specializations.** Let `G` and its cut edge
class `X` be as in Theorem F2 of `proofs/allsize_gaussian.md`, with
`G\X=T\uplus P_m`.  For `m=2` and `m=4`, all finite-witness hypotheses are now
discharged.  Lemmas 7 and 8 there give a nonzero resultant hypersurface such that,
outside it,

\[
 r(G)\ge 3^{5+m}-2^{5+m}+202\,3^m+4\,2^m.
\]

Thus the generic lower bounds are `3893` for `n=7` and `35597` for `n=9`.  The
hypersurface degree bound remains
`4(2n+|E|)S(S-1)`, `S=C(2^n,2)`.  This is a lower bound transported from the
exact decoupled point; exact equality at a generic three-parameter point is not
claimed.

**[UNRESOLVED] Remaining general-F2 hypothesis.** CZ is no longer a hypothesis:
the exact `(1/3,5/3)` core is a sufficient named specialization for every chain
length.  INJ(`n-5`) is the sole remaining hypothesis for chain lengths other than
the now-certified `m=2,4`.  No all-`m` injectivity theorem is proved here.

**[UNRESOLVED] Fully isotropic specialization.** The generic three-parameter
resultant theorem does not locate `y=w`; it remains unresolved whether that
fully isotropic value lies on the exceptional hypersurface.

**[UNRESOLVED] Coupling scope.** Exact core counts are proved at the four displayed
physical points only.  Constancy for arbitrary positive `t`, arbitrary `(t,w)`,
or a symbolic function-field point is not claimed.

---

## 6. [COMPUTATION] Artifacts and reproduction

**[COMPUTATION] Certificate paths.** The machine-readable certificate is
`results/spectral/rank_char0.json`; it cites the prior modular artifact
`results/spectral/allsize_gaussian.json` by SHA-256.  The characteristic-zero
producer used exact integer/rational/SymPy arithmetic, recorded process-CPU
usage `208.659257 s`, peak RSS `100941824` bytes, and 21 passing checks.

**[COMPUTATION] Independent verification.** `tests/test_rank_char0.py` imports none
of the three producer modules.  It rebuilds `R_T` by an explicit Fraction
triple product, clears the actual denominator lcm, uses explicit companion,
exterior-square, square, and Kronecker matrices instead of resultants, and then
rebuilds the five decisive modular composite rows from raw graph/parameter
inputs.  The final recorded run passed all 21 checks and peaked at `152748032` bytes.

```bash
.venv/bin/python experiments/e166_rank_char0_spectrum.py
.venv/bin/python experiments/e167_rank_char0_products.py
.venv/bin/python experiments/e168_rank_char0.py
.venv/bin/python tests/test_rank_char0.py
```
