# Selected-row Callen compression: the Walsh rank and the pair-row boundary

Artifacts: `experiments/e246_callen_walsh_compression.py`,
`results/correlations/callen_walsh_compression.json`, and the clean-room
verifier `tests/test_callen_walsh_compression.py`.  The producer output is
created only after all exact, provenance, scope, CPU, and RSS gates pass.

## 1. Result and exact scope

**[THEOREM — every even local degree at positive coupling].** Let a selected pivot have
`d=2m` neighbours with integer `m>=1`, let `0<v=tanh(K)<1`, and let `c_j(v)`
be the exact degree-`d` Callen coefficient for any odd `j`-subset.  With

\[
 E^\circ=\{U\subseteq N: |U|\text{ is even},\ U\ne N\},
 \qquad
 O=\{W\subseteq N: |W|\text{ is odd}\},                 \tag{1.1}
\]

put

\[
 L_d(U,W)=c_{|U\mathbin\triangle W|}(v),
 \qquad U\in E^\circ,\quad W\in O.                       \tag{1.2}
\]

Then

\[
 \operatorname{rank}L_d
   =2^{d-1}-\frac12\binom d m,
 \qquad
 \dim\ker L_d^{\mathsf T}
   =\frac12\binom d m-1.                                 \tag{1.3}
\]

Every `lambda` in the left kernel satisfies

\[
 \lambda(\varnothing)=0.                                  \tag{1.4}
\]

Thus, in this fixed selected-pivot/far-mark block, no linear combination of
selected Callen rows that cancels every leakage column can have a nonzero
coefficient on the pair row `U=empty`.

**[FINITE EXACT COMPUTATION].** At the exact control value `v=1/3`, the
producer and an independent reconstruction verify all coefficients, Walsh
layers, and matrix ranks for `d=2,4,6`.  For `d=6`, the selected raw matrix is
`31 x 32` of rank `22`; after deleting the empty row it is `30 x 32` of rank
`21`.  Their left nullities are both `9`.  Cubic inversion gives `19` row
classes and `16` column classes.  The folded ranks are respectively `13` and
`12`, and their left nullities are both `6`.  Six explicit folded kernel
vectors are given in Section 6.

**[SCOPE].** The result is exactly: six exact inversion-invariant
higher-template relations and the all-even-degree obstruction to pair-row
leakage cancellation in the fixed selected-pivot/far-mark Callen block; no
susceptibility closure and no contact-inclusive or radius-expanded no-go.
In particular, (1.4) is a boundary for this matrix, not a proof that the pair
correlation cannot be determined by a larger system.

**[EXTERNAL LEMMA].** None is needed.  The Walsh diagonalization and every
rank statement used below are proved directly.

## 2. The selected Callen block

Write

\[
 \sigma_A=\prod_{a\in A}\sigma_a,
 \qquad C(A)=\langle\sigma_A\rangle .
\]

Fix a pivot `x`, its neighbour set `N`, and a marked site `z` outside
`{x} union N`.  The degree-`d` Callen expansion is

\[
 \tanh\!\left(K\sum_{i\in N}\sigma_i\right)
   =\sum_{\substack{A\subseteq N\\|A|\text{ odd}}}
      c_{|A|}(v)\,\sigma_A .                              \tag{2.1}
\]

For an even `U subset N`, use `sigma_z sigma_U` as the factor independent of
`sigma_x`.  The corresponding row is

\[
 R_U(z):=
 C(\{x,z\}\cup U)
 -\sum_{\substack{A\subseteq N\\|A|\text{ odd}}}
       c_{|A|}(v)C(\{z\}\cup(U\mathbin\triangle A))=0.    \tag{2.2}
\]

Relabel the odd leakage support as `W=U triangle A`.  Its coefficient is
`c_|U triangle W|`, which is precisely (1.2).  The row `U=empty` has
`C({x,z})` on the left and is therefore the pair row.  The selected block
retains all even rows except `U=N`.

The far-mark condition is load-bearing only for the stated contact-free
interpretation: `z` never cancels against the pivot or a neighbour spin.
Symmetric differences among neighbour factors in (2.2) remain, as they must,
but they are not marked-site contacts.

## 3. Exact Callen coefficients and their Walsh spectrum

For a fixed odd set `A` of size `j`, Walsh interpolation gives

\[
 c_j(v)=2^{-d}\sum_{s\in\{-1,+1\}^d}
  \tanh\!\left(K\sum_{i=1}^d s_i\right)\prod_{i\in A}s_i. \tag{3.1}
\]

Permutation symmetry makes this depend only on `j`.  Conversely, Walsh
inversion of (3.1) is exactly (2.1).

Define the odd kernel on the group `(Z/2Z)^d` by

\[
 k(A)=\begin{cases}
 c_{|A|}(v),&|A|\text{ odd},\\
 0,&|A|\text{ even}.
 \end{cases}                                               \tag{3.2}
\]

For `Q subset N`, write `chi_Q(A)=(-1)^|Q intersect A|`.  Evaluating the
Walsh expansion at the spin assignment which is `-1` on `Q` and `+1` on its
complement gives the full convolution eigenvalue

\[
 \widehat k(Q)
 =\sum_A k(A)\chi_Q(A)
 =\tanh\!\bigl(K(d-2|Q|)\bigr).                           \tag{3.3}
\]

Because `0<v=tanh(K)<1` implies `K>0`, the last expression is zero if and only
if

\[
 |Q|=m.                                                    \tag{3.4}
\]

There are exactly `binom(d,m)` zero Walsh eigenvalues.  No numerical
approximation or physical critical-coupling input occurs here.

## 4. Rank of the full and selected blocks

Convolution by the odd kernel `k` swaps the even and odd coordinate
subspaces.  In the even/odd ordering its two off-diagonal blocks are
transposes, since `k(A triangle B)=k(B triangle A)`.  The full convolution has
one eigenvalue (3.3) per Walsh character, so its rank is

\[
 2^d-\binom d m.
\]

The two off-diagonal blocks have equal rank.  Hence the full even-by-odd block
has rank

\[
 r_d=2^{d-1}-\frac12\binom d m.                           \tag{4.1}
\]

It remains to justify that deleting the `U=N` row does not lower this rank.
Let `lambda` be a function on the even subsets and extend it by zero on odd
subsets.  Its Fourier transform satisfies

\[
 \widehat\lambda(Q)=\widehat\lambda(Q^c),                 \tag{4.2}
\]

because `chi_Q(U)=chi_Qc(U)` whenever `|U|` is even.  A full-block left-kernel
vector obeys `lambda*k=0`; by (3.3), its Fourier transform is supported on the
middle layer `|Q|=m`.  Conversely, any middle-layer Fourier data satisfying
(4.2) gives a full-block left-kernel vector.  Complementation has no fixed
point among subsets, so this kernel has dimension

\[
 t_d=\frac12\binom d m.                                   \tag{4.3}
\]

Fourier inversion on the middle layer gives

\[
 \begin{aligned}
 \lambda(\varnothing)
   &=2^{-d}\sum_{|Q|=m}\widehat\lambda(Q),\\
 \lambda(N)
   &=2^{-d}\sum_{|Q|=m}(-1)^{|Q|}\widehat\lambda(Q)
     =(-1)^m\lambda(\varnothing).                         \tag{4.4}
 \end{aligned}
\]

A left dependency among the selected rows is exactly a full-block dependency
extended by `lambda(N)=0`.  The evaluation functional at `N` is nonzero on
the full kernel: the middle character `chi_Q` restricted to even subsets has
value `(-1)^m` there.  Therefore `lambda(N)=0` removes exactly one of the
`t_d` kernel dimensions.  The selected left nullity is `t_d-1`, and its rank
is

\[
 (2^{d-1}-1)-(t_d-1)=2^{d-1}-t_d=r_d,                    \tag{4.5}
\]

which proves (1.3).  Equation (4.4) simultaneously proves (1.4).

There is also a useful corrected deletion criterion.  Every selected left
kernel vector already has `lambda(empty)=0`, so deleting the empty row leaves
the left nullity unchanged and lowers the rank by exactly one.  This is why
the `d=6` ranks are `22` and `21`, not two copies of `22`.

## 5. Exact finite controls at `v=1/3`

From (3.1), with exact rational arithmetic,

| degree | exact nonzero Callen coefficients |
|---:|---|
| `2` | `c1=3/10` |
| `4` | `c1=177/680`, `c3=-27/680` |
| `6` | `c1=4143/17680`, `c3=-27/1040`, `c5=243/17680` |

The direct matrix controls are:

| `d` | full even-by-odd rank | selected `E^o x O` shape/rank/nullity | delete-empty shape/rank/nullity |
|---:|---:|---:|---:|
| `2` | `1` | `1 x 2 / 1 / 0` | `0 x 2 / 0 / 0` |
| `4` | `5` | `7 x 8 / 5 / 2` | `6 x 8 / 4 / 2` |
| `6` | `22` | `31 x 32 / 22 / 9` | `30 x 32 / 21 / 9` |

For each degree the computation also sums every Walsh layer and obtains
`tanh(K(d-2|Q|))` exactly, with precisely the middle layer zero.  These finite
checks control the implementation; the proof of (1.3) is the general argument
in Sections 3 and 4, not extrapolation from the table.

## 6. Cubic inversion and six explicit higher-template relations

For `d=6`, order the neighbour directions as

\[
 (x_+,x_-,y_+,y_-,z_+,z_-).
\]

Let `iota` be the single cubic inversion which
simultaneously swaps every plus direction with its minus direction.  Among the
`31` selected even rows, seven are fixed by `iota`; hence Burnside gives

\[
 (31+7)/2=19
\]

row classes.  An `iota`-fixed subset is a union of complete axis pairs and is
therefore even, so no odd column is fixed.  The `32` odd columns give `16`
classes.

The folded convention is the following.  A folded row is the sum of the raw
rows in one row orbit, while one representative is retained from each pair of
identical column equations.  Exact elimination at `v=1/3` gives

| folded row set | shape | rank | left nullity |
|---|---:|---:|---:|
| all selected proper-even classes | `19 x 16` | `13` | `6` |
| delete the empty class | `18 x 16` | `12` | `6` |

Here is an explicit basis for that six-dimensional kernel.  For a middle
three-set `Q`, let \([Q]=\{Q,N\setminus Q\}\) be its central complement pair.  The
seven inversion orbits of such pairs, ordered with fixed orbits first, are

\[
\begin{array}{c|l|c}
 i&\mathcal P_i&w_i\\ \hline
1&\{[\{x_+,y_+,z_+\}]\}&1\\
2&\{[\{x_-,y_+,z_+\}]\}&1\\
3&\{[\{x_+,y_-,z_+\}]\}&1\\
4&\{[\{x_-,y_-,z_+\}]\}&1\\
5&\{[\{x_+,x_-,y_+\}],[\{x_+,x_-,y_-\}]\}&2\\
6&\{[\{x_+,y_+,y_-\}],[\{x_-,y_+,y_-\}]\}&2\\
7&\{[\{x_+,x_-,z_+\}],[\{y_+,y_-,z_+\}]\}&2.
\end{array}                                                \tag{6.1}
\]

For even `U`, define

\[
 \Phi_i(U)=\sum_{[Q]\in\mathcal P_i}(-1)^{|Q\cap U|},
 \qquad
 \lambda_\alpha(U)=\sum_{i=1}^7\alpha_i\Phi_i(U).         \tag{6.2}
\]

The summand is independent of the representative of `[Q]` because `U` is
even.  The orbit sum makes `lambda_alpha` inversion invariant.  Every `Q` in
(6.2) is on the zero layer (3.4), so each `lambda_alpha` annihilates the full
leakage convolution.  Its endpoint values are

\[
 \lambda_\alpha(\varnothing)=
   \alpha_1+\alpha_2+\alpha_3+\alpha_4
       +2\alpha_5+2\alpha_6+2\alpha_7,
 \qquad
 \lambda_\alpha(N)=-\lambda_\alpha(\varnothing).         \tag{6.3}
\]

Thus the exact weighted constraint is

\[
 \alpha_1+\alpha_2+\alpha_3+\alpha_4
       +2\alpha_5+2\alpha_6+2\alpha_7=0.                 \tag{6.4}
\]

The following six integer vectors form a basis of (6.4):

\[
\begin{array}{c|rrrrrrr}
 r&\alpha_1&\alpha_2&\alpha_3&\alpha_4&\alpha_5&\alpha_6&\alpha_7\\ \hline
1&-1& 1&0&0&0&0&0\\
2&-1&0& 1&0&0&0&0\\
3&-1&0&0& 1&0&0&0\\
4&-2&0&0&0& 1&0&0\\
5&-2&0&0&0&0& 1&0\\
6&-2&0&0&0&0&0& 1.
\end{array}                                                \tag{6.5}
\]

Let `lambda_r=lambda_alpha` for row `r` of (6.5).  Combining the Callen rows
(2.2), the leakage coefficient of every odd `W` is

\[
 \sum_{U\in E^\circ}\lambda_r(U)c_{|U\triangle W|}=0.      \tag{6.6}
\]

Consequently the six promised higher-template relations are, explicitly,

\[
 \boxed{\quad
 H_r(z):=\sum_{U\in E^\circ}
       \lambda_r(U)C(\{x,z\}\cup U)=0,
 \qquad r=1,2,3,4,5,6.\quad}                              \tag{6.7}
\]

Equations (6.1), (6.2), and the six rows of (6.5) specify every coefficient in
(6.7) as an integer.  The artifact additionally records all `19` folded row
coordinates of each vector.  The verifier reconstructs those coordinates,
checks their rank is `6`, and multiplies them by the independently rebuilt
`19 x 16` matrix.  Every residual is exactly zero.

By (6.3), `lambda_r(empty)=lambda_r(N)=0`.  Therefore (6.7) contains only the
`|U|=2` and `|U|=4` left-hand templates: it gives six exact higher-template
relations, but it does not contain or solve the pair row.

## 7. The periodic `3 x 3 x 3` far-mark control

Take the pivot to be the origin of `(Z/3Z)^3`.  Its six sites at displacement
`plus or minus e_i` are distinct.  The closed neighbourhood therefore has

\[
 1+6=7
\]

sites, leaving

\[
 27-7=20                                                    \tag{7.1}
\]

far marks.  The producer enumerates all `27` coordinate triples and records
the `20`-element complement.  The verifier rebuilds it with a separate triple
loop.  Every recorded mark is unequal to the pivot and to all six neighbour
variables, so the marked-site pivot-contact and neighbour-contact counts are
both exactly zero.  Thus (2.2) has no contact terms for any of these marks.

## 8. Computation, provenance, and hard gates

**[COMPUTATION].** The producer imports `e240_callen_identity_system.py` only
for the exact rational Callen coefficients.  It evaluates those rational
functions with `Fraction`, replays all local Walsh states, constructs the raw
and folded matrices over `Fraction`, performs exact Gaussian elimination, and
writes JSON only after every check succeeds.  It uses `v=1/3` solely as a
rational finite control, with `benchmark_used=false`.

The verifier imports neither the producer nor `e240`.  It reconstructs each
coefficient independently from (3.1) as a formal rational function over
`Q(v)`, verifies the stored numerator and denominator by polynomial cross
multiplication, then rebuilds all matrices, orbit partitions, six vectors,
and the far-mark set.  Matrix and certificate SHA-256 values are recomputed.

The artifact's exact source-hash map is a hard gate for:

- `experiments/e246_callen_walsh_compression.py`;
- `tests/test_callen_walsh_compression.py`;
- `experiments/e240_callen_identity_system.py`; and
- `experiments/e243_callen_termwise_closure.py`.

Thus the current `e240` coefficient source and current `e243` scope
predecessor are both bound, not merely named.  Changing either source makes
the clean-room verifier fail until a new artifact is produced from the
reviewed sources.

The producer records process CPU and peak RSS with a strict RSS limit below
`2 GiB`.  On Darwin, both producer and verifier use
`mach_task_basic_info.resident_size_max` through `task_info` flavor `20`, not
the inherited `ru_maxrss`; on Linux they use `getrusage` with the documented
KiB-to-byte conversion.  The measurement method itself is recorded and
verified.  Timing and memory measurements are metadata only; every
mathematical claim uses integers, exact fractions, or formal rational
functions.

## 9. What is not proved

The theorem and the six relations do **not** establish any of the following:

1. pair-correlation closure or susceptibility closure;
2. a no-go for systems that include marked-site contacts;
3. a no-go after enlarging the radius, using more pivots, or adding other
   observable templates;
4. a no-go for nonlinear or auxiliary-variable compressions; or
5. an exact critical coupling or thermodynamic solution.

The positive result is the six relations (6.7).  The negative result is only
the exact pair-row coefficient obstruction (1.4) inside the fixed
selected-pivot/far-mark block (1.2).
