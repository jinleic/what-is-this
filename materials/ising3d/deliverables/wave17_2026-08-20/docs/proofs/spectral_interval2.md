# Forced-value crossings obstruct connected Theorem S for the open `2x3` layer

Artifacts: `experiments/e81_spectral_interval2.py`,
`results/spectral/interval2.json`, and
`tests/test_spectral_interval2.py`.

**Status tags.** **[FALSIFIED]** for universal forced-eigenvalue absence on
any connected interval containing one of the certified crossings;
**[THEOREM]** for the four exact c-jump crossing certificates and the three
safe closed intervals; **[LEMMA]** for the multiplicity-safe forced-value
condition; **[COMPUTATION]** for finite exact inertia calculations;
**[UNRESOLVED]** for uniqueness/minimal polynomials and exhaustive crossing
classification; **[EXTERNAL, abstract-only]** for standard finite-dimensional
continuity/real-closed-field facts specified below. No external file or
citation is claimed or used as a numerical authority.

## 1. Result

Let `mu_0(t) <= ... <= mu_63(t)` be the sorted eigenvalues, with
multiplicity, of the positive polynomial representative

\[
 M(t)=(2t)^7q(t)^4R(t),\qquad
 q(t)=\frac{1+t^2}{2t},\qquad
 R(t)=P_t\operatorname{diag}(q(t)^{e_s})P_t,
\]

for the open `2x3` layer.  Positive scalar multiplication does not affect the
condition below.  Define

\[
 f(t)=\frac{\mu_1(t)\mu_2(t)}{\mu_0(t)}.
\]

**[THEOREM — exact crossing obstruction].** Each row in the following table
contains at least one real algebraic parameter \(t_r^*\) at which
\(f(t_r^*)\) is an eigenvalue of \(M(t_r^*)\).  The displayed bounds are
exact rational localizing brackets for this existential statement; their width
is \(1/16{,}777{,}216{,}000<10^{-10}\).

| `r` | exact rational bracket for an attained forced value | exact c-count change |
|---:|---|---|
| 1 | \([2248027753/8388608000,\;4496055507/16777216000]\) | \(4\to5\) |
| 2 | \([4881706207/16777216000,\;152553319/524288000]\) | \(5\to6\) |
| 3 | \([2522992499/8388608000,\;5045984999/16777216000]\) | \(6\to7\) |
| 4 | \([553214907/1677721600,\;5532149071/16777216000]\) | \(7\to8\) |

Their decimal locations, for orientation only, are approximately
`0.26798579135`, `0.29097236440`, `0.30076414338`, and
`0.32974178019`.  Decimals decide no part of the theorem.

Consequently, a statement asserting

\[
 f(t)\notin\operatorname{Spec}M(t)
\]

for every `t` in a connected closed interval that contains any one of these
parameters is **false**.  This is a genuine obstruction, not merely a failure
of a Lipschitz certificate.

The positive corrected statements obtained here are:

\[
\begin{aligned}
J_{1/4}&=[524287/2097152,\;524289/2097152], &&|J_{1/4}|=1/1048576,\\
J_{3/10}&=[50331643/167772160,\;50331653/167772160], &&|J_{3/10}|=1/16777216,\\
J_{1/3}&=[8388605/25165824,\;8388611/25165824], &&|J_{1/3}|=1/4194304.
\end{aligned}
\]

**[THEOREM — safe intervals].** On every `t` in each displayed closed rational
interval, \(f(t)\notin\operatorname{Spec}M(t)\), hence the full 64-value
spectrum cannot be a positive six-mode Gaussian subset-product spectrum.
The widest new interval has width \(1/1048576\), strictly exceeding the
largest wave-8 component width \(2\cdot10^{-8}\).

The theorem does **not** say that the forced relation fails everywhere outside
these three intervals, nor that the four brackets exhaust all crossings.

## 2. Exact matrix representative and positivity

For a layer state `s`, set

\[
 e_s=\frac{\sum_{(i,j)\in E}\sigma_i\sigma_j-1}{2}\in\{-4,-3,\ldots,3\}.
\]

The script expands every entry of `M` directly.  It obtains an integer
polynomial of degree at most 26 with nonnegative coefficients.  It also
checks the exact character-sector dimensions under row reflection, column
reflection, and global spin flip:

\[
(14,6,6,6,10,10,6,6).
\]

For `0<t<1`,

\[
(P_t)_{uv}=t^{d_H(u,v)}
 =\bigotimes_{j=1}^{6}\begin{pmatrix}1&t\\t&1\end{pmatrix}_{uv}.
\]

Each one-site factor has determinant \(1-t^2\ne0\).  Since the diagonal
factor `diag(q(t)^e_s)` is strictly positive, `R(t)` and `M(t)` are positive
definite.  In particular \(\mu_0(t)>0\), so `f(t)` is well-defined and
continuous.

## 3. The forced-value lemma without a distinctness hypothesis

**[LEMMA — multiplicity-safe Gaussian forced value].** Let a positive
`2^m`-point Gaussian spectrum be

\[
\left\{c\prod_{j=1}^{m}r_j^{\epsilon_j}:
 \epsilon_j\in\{0,1\}\right\},\qquad c>0,\quad r_j\ge1.
\]

If \(\lambda_0\le\lambda_1\le\lambda_2\) are its first three values,
counted with multiplicity, then

\[
\frac{\lambda_1\lambda_2}{\lambda_0}
\]

is also a value in the spectrum.

*Proof.* The smallest value is `c`.  If at least two modes attain the smallest
factor `r`, then \(\lambda_1=\lambda_2=cr\) and the forced value is
\(cr^2\), realized by the product of two distinct minimal modes.  If the
smallest factor is unique and greater than one, the next smallest factor is
`r_2`, and the forced value is \(cr_1r_2\).  If a smallest factor is one,
then either the first two values are both `c`, or the forced value is the next
one-mode product; in either case it is already in the spectrum. `[]`

This repairs an unnecessary restriction in the previous interval argument:
strict inequalities among the three lowest eigenvalues are not needed for the
necessary Gaussian condition.

## 4. Exact point certificates and the c-jump theorem

For a rational `t`, the program isolates each of
\(\mu_0,\mu_1,\mu_2\) using exact Sylvester inertia counts of the generalized
character-sector pencils.  If

\[
\mu_j(t)\in[l_j,u_j],
\]

then exact positive interval arithmetic gives

\[
 f(t)\in\left[\frac{l_1l_2}{u_0},\frac{u_1u_2}{l_0}\right].
\]

The count

\[
 c(t)=\#\{j:\mu_j(t)<f(t)\}
\]

is certified when the exact inertia counts at outward rational endpoints of
this forced-value window agree.  At the two rational endpoints of each table
row, the windows are clear, and their exact counts are the displayed distinct
integers.

**[LEMMA — c-jump implies forced equality].** If two points `a<b` have clear
forced-value windows and `c(a) != c(b)`, then some \(t^*\in(a,b)\) satisfies
\(f(t^*)\in\operatorname{Spec}M(t^*)\).

*Proof.* Ordered eigenvalues of a continuous real-symmetric matrix family are
continuous; this is the standard Weyl-continuity fact **[EXTERNAL,
abstract-only]**.  So are \(f\) and every difference
\(\mu_j-f\).  If none of those differences vanished on `[a,b]`, all their
signs, hence the integer `c`, would be locally constant and therefore
constant on that connected interval.  The certified endpoint inequality is a
contradiction. `[]`

Every condition in the endpoint certificates is integer or `Fraction`
arithmetic.  Float64 only suggested initial brackets later accepted or
rejected by exact inertia.

Why may the crossing be called algebraic?  The eigenvalues and the forced
relation form a first-order semialgebraic condition over \(\mathbb Q\): it is
expressible from the characteristic polynomial of the rational-polynomial
matrix, root ordering, and \(ad=bc\).  The c-jump lemma proves that this
condition has a real solution in each rational bracket.  The real algebraic
numbers form a real closed field, so the standard transfer/quantifier-
elimination fact **[EXTERNAL, abstract-only]** gives an algebraic solution in
that same bracket.  This experiment deliberately does **not** compute its
minimal polynomial, a resultant, a Sturm chain, or uniqueness within a
bracket.

The exact grid scan on `[1/4,7/20]` records the runs

\[
4\quad\to\quad5\quad\to\quad6\quad\to\quad7\quad\to\quad8,
\]

with coarse transition brackets
`[267/1000,268/1000]`, `[290/1000,291/1000]`,
`[300/1000,301/1000]`, and `[329/1000,330/1000]`.  The fine brackets above
are 24 exact rational bisections of these coarse witnesses.

## 5. Closed-interval absence propagation

Let `J=[c-h,c+h]` and let

\[
L_J\ge\sup_{t\in J}\lVert M'(t)\rVert_2,
\qquad \rho=L_Jh.
\]

Nonnegative coefficient lists yield exact rational upper bounds for both the
Frobenius norm and the absolute row-sum norm of `M'`; their minimum bounds the
spectral norm.  Weyl gives

\[
|\mu_j(t)-\mu_j(c)|\le\rho.
\]

The center enclosures then imply, without any distinctness assumption,

\[
 f(t)\in
 \left[
 \frac{(l_1-\rho)(l_2-\rho)}{u_0+\rho},
 \frac{(u_1+\rho)(u_2+\rho)}{l_0-\rho}
 \right]=:F_J.
\]

If an eigenvalue at `t` equalled `f(t)`, Weyl would place an eigenvalue of
`M(c)` in `F_J+[-rho,rho]`.  Equal exact inertia counts at outward endpoints
of that center window rule this out on all of `J`.

For example, on `J_{1/4}` the exact derivative upper bound selected by the
program is

\[
L_J=\frac{181723715527428907}{2500000000000000},\qquad
\rho=\frac{181723715527428907}{5242880000000000000000},
\]

and the two exact center-window counts are both `4`.  The standalone verifier
rebuilds `M` independently from `R=P_t\operatorname{diag}(q^e)P_t`, checks the
polynomial expansion agrees at `t=1/4`, reconstructs these bounds, and uses a
full unsectorized 64-by-64 fraction-free inertia calculation to reproduce the
`4,4` result.

## 6. All-interval one-dimensional control

**[LEMMA — open-chain control].** The open `n=4` chain is Gaussian throughout
`[1/5,1/2]`, and therefore its forced value is present throughout that
interval.

Use Jordan--Wigner Majorana strings

\[
\gamma_{2j}=X_0\cdots X_{j-1}Z_j,
\qquad
\gamma_{2j+1}=X_0\cdots X_{j-1}X_jZ_j.
\]

The exact Pauli-support checks in the script verify

\[
\gamma_{2j}\gamma_{2j+1}\sim X_j,
\qquad
\gamma_{2j+1}\gamma_{2j+2}\sim Z_jZ_{j+1},
\]

and all distinct Majorana strings anticommute.  Now

\[
P_t\propto\prod_j e^{a(t)X_j},\quad\tanh a(t)=t,
\qquad
D_t\propto\prod_j e^{K(t)Z_jZ_{j+1}},\quad e^{2K(t)}=q(t).
\]

Thus `P_tD_tP_t` is a product of exponentials of Majorana bilinears.  Their
commutators close on Majorana bilinears, so these products form the
finite-dimensional Gaussian/Spin group; a positive member has a positive
spectrum.  Combining this with the lemma in section 3 proves forced-value
presence for every control parameter.  The standard canonical-form statement
for positive Gaussian operators is the only finite-dimensional
**[EXTERNAL, abstract-only]** fact used in this paragraph.

As a direct machinery control, exact forced windows were repeatedly tightened
at `t=1/5,1/3,1/2`; each remained occupied, so the absence procedure correctly
refused to reject the Gaussian chain.  An explicit two-mode spectrum
`{1,A(t),B(t),A(t)B(t)}`, with `A=1+t^2` and `B=A^3`, independently verifies
symbolically that the multiplicity-safe forced formula returns the existing
subset product `A(t)B(t)`.

## 7. Scope and unresolved remainder

1. **[FALSIFIED]** The particular forced-value absence criterion cannot prove
   non-Gaussianity on a connected interval containing any certified bracket:
   it is actually false at an attained parameter inside each bracket.
2. **[THEOREM]** The three stated safe intervals are honest continuum no-go
   statements and are wider than each wave-8 component.
3. **[UNRESOLVED]** No resultant/minimal polynomial/Sturm certificate is
   supplied for a unique named root in any bracket, and the scan does not
   prove that only four crossings occur in `[1/5,1/2]`.
4. **[UNRESOLVED]** The forced relation is necessary, not sufficient.  Its
   attainment at `t_r^*` does not make the `2x3` layer Gaussian and does not
   negate independent spectral or algebraic obstructions.
5. No critical-coupling benchmark was used to choose a point, a crossing
   bracket, or an interval.
