# Parametric spectral-Gaussianity obstruction for the `2x3` layer

**Status tags.** **[THEOREM]** proved for the four explicit rational subintervals below;
**[LEMMA]** proved here; **[COMPUTATION]** finite exact computation;
**[UNRESOLVED]** not proved; **[EXTERNAL]** cited from an existing manifest entry.

Script: `experiments/e48_parametric_gaussianity.py`  
Certificate: `results/spectral/parametric_gaussianity.json`  
Standalone check: `tests/test_parametric_gaussianity.py`

## 1. Outcome: the requested connected interval is unresolved

**[UNRESOLVED]** The requested statement for every rational or real
`t in [1/5,1/2]` has **not** been proved. The interval-propagation argument below yields four
nonzero but disjoint rational neighbourhoods. They do not chain across the gaps, so this note does
not silently promote sampled couplings to a continuum theorem.

Numerical reconnaissance, used only to understand the failed step and never as a certificate,
shows the forced value approaching and apparently crossing individual layer eigenvalue branches
inside `[1/5,1/2]`. Consequently, the particular necessary Gaussian relation used by Theorem S may
hold accidentally at isolated couplings even though the full spectrum remains non-Gaussian. A
proof of the full requested interval would then need a different necessary Gaussian relation at
such crossings. No exact crossing-existence statement is claimed here.

The exact positive results are:

**[THEOREM S-parametric, four components].** For the `2x3` open layer, its full 64-eigenvalue
transfer-operator spectrum is not the subset-product spectrum of one positive Gaussian operator on
six fermionic modes for every `t` in each of

\[
\begin{aligned}
 I_1&=[19999999/10^8,\;20000001/10^8],\\
 I_2&=[99999997/(3\cdot10^8),\;100000003/(3\cdot10^8)],\\
 I_3&=[39999999/10^8,\;40000001/10^8],\\
 I_4&=[499999999/10^9,\;500000001/10^9].
\end{aligned}
\]

These are genuine all-parameter statements, not grid-point computations. Their centres are
`1/5, 1/3, 2/5, 1/2`; their exact half-widths are respectively
`10^-8, 10^-8, 10^-8, 10^-9`.

## 2. Exact parametrisation

As in `experiments/e38_gaussianity_certificate.py`, put

\[
 t=\tanh(K^*/2),\qquad q(t)=e^{2K}=\frac{1+t^2}{2t}.
\]

After deleting a positive global scalar, the symmetric layer operator is

\[
 R(t)=P_t\,\operatorname{diag}(q(t)^{e_s})P_t,
 \qquad (P_t)_{uv}=t^{d_H(u,v)}.
\]

For the seven-bond `2x3` graph, `e_s` ranges from `-4` to `3`. The script instead uses the
positive scalar multiple

\[
 M(t)=(2t)^7q(t)^4R(t).
\]

**[LEMMA] (integer-polynomial representative).** Every entry of `M(t)` is an integer polynomial
of degree at most 26 with nonnegative coefficients.

*Proof.* A state with exponent `e_s` contributes

\[
 t^{d_H(u,s)+d_H(s,v)}(2t)^{7-(e_s+4)}(1+t^2)^{e_s+4}.
\]

Here `0 <= e_s+4 <= 7`. Expanding the final factor gives only nonnegative integer
coefficients. Summing over states preserves that property. The script constructs these coefficient
lists directly and records that every coefficient is nonnegative. Multiplying by a positive scalar
preserves all spectral ratios and therefore the Gaussian forced relation. `[]`

## 3. Forced value

**[LEMMA] (Gaussian forced value).** If a positive `2^m`-point Gaussian spectrum has sorted
values `lambda_0 < lambda_1 < lambda_2 <= ...`, then

\[
 f=\frac{\lambda_1\lambda_2}{\lambda_0}
\]

is also an eigenvalue.

*Proof.* Write the spectrum as `lambda_0` times all subset products of factors `u_j >= 1`.
Distinctness of the first three values excludes a zero mode and a repeated smallest factor. Thus
`lambda_1/lambda_0` and `lambda_2/lambda_0` are the two smallest one-factor products, while their
two-factor product must also occur. `[]`

This is the same necessary condition as in Theorem S. It is basis-independent because it uses only
the eigenvalue multiset.

## 4. Exact interval-propagation lemmas

Let `J=[c-h,c+h]` be one of the four rational intervals and write `mu_j(t)` for the sorted
eigenvalues of `M(t)`.

**[LEMMA] (Weyl-Frobenius motion; Weyl inequality is [EXTERNAL]).** If

\[
 L_J\ge \sup_{t\in J}\|M'(t)\|_F,
\]

then for every `t in J` and every `j`,

\[
 |\mu_j(t)-\mu_j(c)|\le L_J|t-c|\le \rho_J:=L_Jh.
\]

*Proof.* For real symmetric matrices, the min-max principle gives the Weyl estimate
`|lambda_j(A)-lambda_j(B)| <= ||A-B||_2`. The fundamental theorem of calculus and
`||X||_2 <= ||X||_F` give

\[
 \|M(t)-M(c)\|_2\le |t-c|\sup_{s\in J}\|M'(s)\|_F.
\]

No numerical eigensolver is used in this deduction. `[]`

**[LEMMA] (exact rational derivative majorant).** If `M_uv(t)=sum_k a_uvk t^k`, all
`a_uvk >= 0`, and `0<t<=b`, then

\[
 \sup_{0<t\le b}\|M'(t)\|_F^2
 \le \sum_{u,v}\left(\sum_{k\ge1}ka_{uvk}b^{k-1}\right)^2=:S_b.
\]

The right side is rational when `b` is rational. The script takes the least multiple of
`10^-15` whose square is at least `S_b`; this is an exact outward rational `L_J`. Integer square
comparison verifies the rounding direction. `[]`

Suppose exact Sylvester-inertia bisection at `c` encloses the first three eigenvalues by
`[l_j,u_j]`. The Weyl lemma gives, simultaneously on all of `J`,

\[
 \mu_j(t)\in[l_j-\rho_J,u_j+\rho_J].
\]

**[LEMMA] (standing hypotheses and forced-value modulus).** If

\[
 l_0-\rho_J>0,\quad l_1-u_0-2\rho_J>0,\quad l_2-u_1-2\rho_J>0,
\]

then `mu_0(t)<mu_1(t)<mu_2(t)` throughout `J`, and positivity plus monotonic interval arithmetic
gives

\[
 f(t)\in F_J:=
 \left[
 \frac{(l_1-\rho_J)(l_2-\rho_J)}{u_0+\rho_J},
 \frac{(u_1+\rho_J)(u_2+\rho_J)}{l_0-\rho_J}
 \right].
\]

Thus `F_J` is an exact Lipschitz-style modulus enclosure for the forced value. `[]`

**[LEMMA] (inertia propagation).** If the exact inertia counts of `M(c)` are equal at the two
ends of

\[
 [\inf F_J-\rho_J,\;\sup F_J+\rho_J],
\]

then no `mu_j(t)` equals `f(t)` for any `t in J`.

*Proof.* Equal endpoint counts, with nondegenerate rational shifts, prove that `M(c)` has no
eigenvalue in that closed tested window (the implementation moves a degenerate endpoint only
outward). If `mu_j(t)=f(t)` for some `t`, Weyl motion would put `mu_j(c)` within `rho_J` of
`f(t) in F_J`, hence inside the tested window, contradiction. `[]`

Combining the last lemma with the Gaussian forced-value lemma proves the four-component theorem.

## 5. Exact coupling table

The `K` values here are display-only evaluations at 50 decimal digits; all certified endpoints are
the exact rational `t` values in the first column. On `(0,1)`,

\[
 K(t)=\frac12\log\frac{1+t^2}{2t}
\]

is decreasing, so the endpoint order reverses.

| centre | exact certified `t` interval | corresponding `K` interval (50-dps display) | centre inertia counts |
|---|---|---|---|
| `1/5` | `[19999999/100000000, 20000001/100000000]` | `[0.477755699436795773, 0.477755745590641927]` | `4,4` |
| `1/3` | `[99999997/300000000, 100000003/300000000]` | `[0.255412799882995603, 0.255412823882995603]` | `8,8` |
| `2/5` | `[39999999/100000000, 40000001/100000000]` | `[0.185781769164517566, 0.185781787267965842]` | `8,8` |
| `1/2` | `[499999999/1000000000, 500000001/1000000000]` | `[0.111571775057104879, 0.111571776257104879]` | `8,8` |

The requested `t` interval `[1/5,1/2]` corresponds to
`K in [0.111571775657104878..., 0.477755722513718181...]`.

**Observation only.** The benchmark `K_c ~= 0.221654626` lies inside that *requested but
uncertified-as-a-whole* `K` interval. It is not inside any of the four narrow certified components;
its corresponding `t` is approximately `0.36332099224`. The benchmark was not used to choose,
fit, or tune any certificate.

## 6. Sylvester certificates and control

**[COMPUTATION]** At each centre, all spectral enclosures and window counts are exact. The matrix is
split into four reflection-character sectors of dimensions `24,16,12,12`. In each sector,
fraction-free symmetric Bareiss/LDL elimination counts negative pivots in the generalized pencil
`B-sigma G`; Sylvester's law of inertia is the underlying **[EXTERNAL]** standard theorem. Weyl's
eigenvalue perturbation inequality, used in section 4, is likewise **[EXTERNAL]**; both are cited as
standard matrix-analysis results rather than newly proved from first principles here. The full
rational endpoints, derivative squares, outward square-root bounds, Weyl radii, positivity bounds,
distinctness margins, forced-value enclosures, tested windows, and counts are retained in the JSON
artifact rather than rounded into this note.

The result artifact uses the required `provenance/data/checks` envelope. Float64 proposes initial
brackets only. Every proposed bracket endpoint is accepted only after an exact inertia count, so
float64 decides no claim.

**[COMPUTATION] (negative control).** The same construction at `t=1/3` for the open `n=4` chain
returns different exact counts, `5` and `6`, around the forced-value window. Therefore the
machinery does **not** certify absence for this free-fermion control. As always, differing counts
show only that *some* eigenvalue lies in a positive-width window; they are not an exact equality
certificate.

The standalone test independently rebuilds the polynomial matrix and performs full-matrix exact
Fraction LDL counts (no experiment import and no symmetry reduction) for two stored
subintervals. It also validates the artifact schema, all embedded checks, and the chain control.

## 7. Honest scope and remaining escape routes

The four-component theorem excludes, on those precise intervals, representing the full physical
64-point spectrum as one six-mode positive Gaussian subset-product spectrum, in any Majorana
basis. It does **not** exclude:

1. a restriction of a Gaussian operator on more than six modes, whose larger subset-product
   spectrum contains the physical spectrum only as a selected submultiset;
2. non-physical-parity projections or a pair of Gaussian sectors in a basis whose fermion parity
   is not the physical spin flip;
3. Gaussianity at parameters outside the four certified components;
4. non-Gaussianity throughout `[1/5,1/2]` by some stronger criterion—the failure here is the
   inability of this single forced-value certificate to cross its zero-margin points, not evidence
   that the full spectrum becomes Gaussian.

Accordingly the full connected-interval upgrade remains **[UNRESOLVED]**, while each stated
rational component is a proved continuum result.
