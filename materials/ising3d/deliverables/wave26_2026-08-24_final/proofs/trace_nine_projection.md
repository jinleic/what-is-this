# Finite trace-nine projection for the open `2x4` layer

Artifacts: `experiments/e251_trace_nine_projection.py`,
`results/spectral/trace_nine_projection.json`, and the clean-room verifier
`tests/test_trace_nine_projection.py`. The frozen contract is
`checkpoints/wave26_research_plan.md`.

## 1. Result and boundary

**[THEOREM — finite trace-nine projection].** Let `q=(1+t^2)/(2t)>1` and use the
open-`2x4` polynomial transfer representative and normalized traces fixed in e248. If its full
256-point spectrum is a positive full eight-mode subset-product spectrum, then `q` is a root of one
nonzero polynomial over `Q` of degree at most

\[
182400.
\]

Consequently the trace-four-through-nine incidence can hold at no more than 182,400 complex `q`
values at which the normalized targets are defined, and at only finitely many physical couplings.

**[EXACT COMPUTATION — nonvanishing].** At the exact physical point `q=2`, reduction modulo the
prime `p=2147483647` gives the established 96-dimensional `F4,...,F8` quotient and an invertible
`96x96` multiplication-by-`F9` matrix. Its nonzero determinant makes `F9` a unit in the finite
field quotient; the length-96 Bezout equality excludes infinity and proves the cleared
characteristic-zero resultant nonzero at `q=2`. A second physical point/prime certificate and a
synthetic eight-mode singular control are stored in the artifact.

**[UNRESOLVED].** The primitive `F9` norm is not reconstructed. Its physical roots are not isolated,
and no coefficient branch receives a complete shifted-Stieltjes disposition. This theorem proves
finiteness, not emptiness. It neither lists exceptional couplings nor extends beyond the finite open
`2x4` layer.

## 2. The 96-dimensional coefficient algebra

For positive modes write

\[
\tau_i=z_i+z_i^{-1},\qquad u_i=\tau_i^2\ge4,\qquad
Q(u)=\prod_{i=1}^8(u-u_i).
\]

As proved and independently verified in e248,

\[
Q(0)=r_1^2,\qquad Q(2)=r_2,\qquad Q(3)=r_3/r_1
\]

solve `a0,a2,a1` in the monic degree-eight polynomial `Q`. The remaining variables are
`a7,a6,a5,a4,a3`. The Lucas factors for traces four through nine have degrees

\[
2,2,2,3,4,4,
\]

and their resultant equations have total degrees in the five free coefficients

\[
2,2,2,3,4,3.                                      \tag{2.1}
\]

The last degree drops from four to three after the three fixed evaluations are substituted. The
first five equations have parameter-independent leading forms. Their 35-element grevlex basis has
pure leading powers `(2,2,3,5,9)` and 96 standard monomials. Counting those literal monomials by
degree gives

\[
(1,5,12,19,22,19,12,5,1),                         \tag{2.2}
\]

which also follows from the complete-intersection Hilbert series

\[
(1+t)^3(1+t+t^2)(1+t+t^2+t^3).
\]

The maximum standard degree is eight. Multiplication by the cubic `F9` therefore requires reductions
only through total degree eleven.

The no-projective-zero argument in e248 is load-bearing: homogenizing `F4,...,F8` adds no point at
infinity, so the affine quotient has scheme length 96 over `Q(q)` and at every physical `q>1`.

## 3. Exact modular nonvanishing

The producer reconstructs the abstract Lucas incidence without importing e248 code. It computes each
small resultant as the determinant of multiplication by `Q mod f` in `Q[u]/(f)`, where
`deg f<=4`. After substituting an exact rational `q0` modulo a safe prime, it performs these steps:

1. evaluate all nine landed integer trace polynomials by modular Horner arithmetic;
2. evaluate every normalized target and reject zero residues of `q`, `q^2-1`, `H1`, or `H2` before
   an inversion;
3. compute the specialized five-equation grevlex basis and require exactly the generic 35 leading
   monomials and 96 standard monomials;
4. reduce all 96 products `F9*b` by a degree-descending finite-field reducer;
5. form the `96x96` multiplication matrix and compute its determinant by exact modular Gaussian
   elimination.

The primes satisfy `p^2<2^63`, so every vectorized product fits the guarded integer representation.
A nonzero modular determinant makes `F9` a unit in the specialized finite algebra, also after
extension to the algebraic closure; hence the six rows have no affine common zero in that fiber.
The matched initial ideal has length `96=2*2*2*3*4`, the full Bezout product of `F4,...,F8`, so
that fiber has no residual point at infinity. Therefore the six cleared homogenized rows have no
projective common zero over the finite-field algebraic closure, and their universal resultant is
nonzero at the specialized `q0`.

All target and coefficient denominators are units modulo the chosen prime. The finite-field rows
are the literal reductions of the uncleared characteristic-zero equations. The corresponding
cleared rows differ by the nonzero scalars `d_k(q0)D(q0)^(m_k)`, which are also units modulo `p`;
multihomogeneity therefore changes the resultant only by a nonzero factor. If the cleared
resultant at the exact rational `q0` were zero in characteristic zero, its good-prime reduction
would also be zero. The observed finite-field nonvanishing proves that characteristic-zero
specialization nonzero, and hence that the resultant is not the zero polynomial in `q`. In the
equivalent finite-flat quotient picture, this is the nonvanishing multiplication norm of `F9`.

The positive control chooses the literal modes

\[
u_i=4,5,6,7,8,9,10,11.
\]

It builds `Q` and all nine targets from these integers. The known five free coefficients annihilate
all six rows, and multiplication by `F9` is singular. This catches a determinant implementation that
would return nonzero independently of the incidence.

## 4. Projection degree

Let the canceled numerator/denominator degree pairs of `Q(0),Q(2),Q(3)` be

\[
(22,12),\quad(32,22),\quad(40,30).
\]

Let `D` be the product of the three denominators. Then

\[
\deg D=12+22+30=64.
\]

After putting the three fixed coefficients over `D`, every coefficient of `D Q(u)` has `q`-degree at
most

\[
\max\{22+64-12,\ 32+64-22,\ 40+64-30\}=74.       \tag{4.1}
\]

For a Lucas factor of degree `m`, write the target as `n_k/d_k`. Clearing
`Res(f_k,Q)-n_k/d_k` by `D^m d_k` gives a coefficient-degree bound

\[
\delta_k=\max(74m+\deg d_k,\ 64m+\deg n_k).
\]

For `k=4,...,9`, using `m=(2,2,2,3,4,4)` gives

\[
\delta=(192,200,224,296,384,392).                 \tag{4.2}
\]

The safe fourth degree in the last clearing step is retained even though the last affine equation
has coefficient-variable degree three. This avoids assuming that the cancellation also lowers every
`q` denominator.

For six homogeneous equations of degrees `(2,2,2,3,4,3)` in six projective variables, the
resultant is multihomogeneous in each row's coefficients. Its row multidegrees are the products of
the other five equation degrees:

\[
(144,144,144,96,72,96).                            \tag{4.3}
\]

Equations (4.2) and (4.3) therefore give

\[
\deg_q R\le
192(144)+200(144)+224(144)+296(96)+384(72)+392(96)
=182400.                                            \tag{4.4}
\]

The cleared resultant may contain target-denominator factors, but those are strictly positive for
`q>1`: the raw trace polynomials have nonnegative coefficients and positive configuration sums, and
`q>0`, `q^2-1>0`. Extra clearing factors can enlarge the candidate set but cannot invalidate the
upper bound. The modular witness proves `R` is nonzero, completing the theorem.

## 5. Lucas-product gates do not shortcut the projection

Set `x=u-4>=0`. The six factors become

\[
\begin{aligned}
f_4&=x^2+4x+2,& f_5&=x^2+3x+1,& f_6&=x^2+4x+1,\\
f_7&=x^3+5x^2+6x+1,&
f_8&=x^4+8x^3+20x^2+16x+2,&
f_9&=x^4+7x^3+15x^2+10x+1.
\end{aligned}
\]

Coefficientwise nonnegative identities give the necessary product inequalities

\[
\begin{gathered}
f_4-f_6=1,
\quad f_7-f_5=x(x^2+4x+3),
\quad f_8-f_9=f_7,
\quad f_4^2-f_8=2,\\
f_4f_6-f_5^2=2x^3+8x^2+6x+1,\\
f_5f_7-f_6^2=x^5+7x^4+14x^3+6x^2+x,\\
f_6f_8-f_7^2=2x^5+16x^4+42x^3+40x^2+12x+1.
\end{gathered}                                      \tag{5.1}
\]

Multiplying (5.1) over eight modes yields seven cheap necessary inequalities among the normalized
trace targets. For the literal open-`2x4` targets, the producer clears their positive denominators
and shifts `q=1+y`. Every numerator has strictly positive coefficients. Thus all seven gates hold
throughout `q>1`; none can replace the `F9` projection. This is a closed shortcut class, not positive
evidence for a mode representation.

## 6. Verification and resource scope

The producer and verifier share no implementation code. The producer uses small quotient norms,
heap reduction, and vectorized finite-field elimination. The verifier uses direct resultants, a
separate highest-term reducer, and scalar Gaussian elimination. It rebuilds every witness/control
matrix, the Hilbert count, degree arithmetic, Lucas gates, complete theorem/scope record, data hash,
and source hashes.

Every heavy stage is single-process. The 1.24-GiB e248 trace recurrence is consumed by content hash
rather than rerun. The coarse degree bound would require 182,401 interpolation nodes before a tighter
valuation is proved. Direct interpolation with a fresh five-variable Groebner basis at each node is
therefore rejected by the frozen resource gate. A future full norm needs a fixed-leading-form lift or
fraction-free border template, followed by a proved denominator, degree, and coefficient-height
envelope. Agreement across a few primes is never a substitute for height-bounded CRT reconstruction.

Recent Singular, Sage, msolve, and FLINT source was surveyed for modular standard bases,
multiplication matrices, polynomial-matrix determinants, and exact real-root tools. None is installed
in this project environment, and e251 adds no dependency. Those codebases inform the next algorithmic
choice; no external numerical or theorem result is imported into the proof above.
