# Replica-column trace nine for the open `2x4` layer

Artifacts: `experiments/e248_replica_trace_nine.py`,
`results/spectral/replica_trace_nine.json`, and the clean-room verifier
`tests/test_replica_trace_nine.py`.

## 1. Result and honest boundary

**[THEOREM — replica-column identity].** Let `B_L(q)` be the polynomially shifted transfer
representative for the open `2xL` layer used in e238. For every `L,k>=1`,

\[
 H_{k,L}(q):=\operatorname{tr}(B_L(q)^k)
\]

is the equality enumerator of `k` cyclic layer replicas. It is computed exactly by a column
recurrence with `4^k` states, degree `k(5L-2)`, and coefficient sum `2^(2Lk)`. The recurrence never
forms a dense symbolic `2^(2L) x 2^(2L)` matrix power.

**[EXACT COMPUTATION].** For the open `2x4` layer, the producer constructs all nine integer
polynomials `H_1,...,H_9`. The largest step uses 262,144 states, degree 162, and two packed uint64
arrays. A fixed final-column state has at most `2^54` partial configurations; summing 128 rows at a
time is bounded by `2^61<2^64`. Only the outer sum of the `2^11` block totals uses Python integers.

\[
9,31,43,67,79,103,115,139,151.
\]

The resulting normalized traces define the first correctly dimensioned eight-mode Lucas incidence:
`Q(0),Q(2),Q(3)` leave five mode-polynomial coefficients, `F4,...,F8` form a zero-dimensional
coefficient ideal over `Q(q)`, and `F9` is the first safe equation that can cut the coupling.
The equation degrees are

\[
2,2,2,3,4,3.
\]

The parameter-independent leading forms have a 35-element grevlex Groebner basis, no projective
zero, and exactly 96 standard monomials. Their pure-power bounds are `(2,2,3,5,9)` and the maximum
standard degree is eight.

**[THEOREM — coefficient quotient].** Over `Q(q)`, the affine scheme cut out by
`F4,...,F8` has length exactly 96. The same conclusion holds at every physical `q>1`, where all
recorded target denominators are strictly positive. This is a scheme-length theorem; an explicit
lifted Groebner basis in the displayed 96 monomials is not claimed.

**[UNRESOLVED].** This front does not materialize the `96x96` multiplication norm of `F9`, isolate
its `q>1` roots, or apply shifted Stieltjes/Hermite disposition. Therefore it does not prove that the
open `2x4` physical exceptional set is empty. It advances the exact frontier from a fixed-coupling
pair-product obstruction to the complete trace-nine input and a finite coefficient quotient, while
leaving the univariate projection and semialgebraic branch stage open.

**[RESOURCE OBSERVATION].** An optional `--quotient-probe` attempted to substitute the physical
`Q(q)` targets and build the actual coefficient quotient. The outer command reached its 900-second
wall under background QoS before returning a quotient record. It wrote no artifact and supplies no
mathematical negative result. The verified trace/leading-ideal payload above is independent of that
aborted optional stage.

## 2. Derivation of the column recurrence

The open `2xL` layer has `n=2L` sites and `m=3L-2` in-layer edges. For a spin state `s`, the shifted
diagonal exponent in e238 is `m` minus the number of disagreeing in-layer edges. The off-diagonal
kernel contributes `n-d_H(s,t)`, the number of equal sites in consecutive replicas. Expanding a
cyclic trace therefore gives

\[
H_{k,L}(q)=
\sum_{s_0,\ldots,s_{k-1}}
q^{\sum_r \#\{\text{satisfied in-layer edges of }s_r\}
 +\sum_r \#\{\text{equal sites of }s_r,s_{r+1}\}},       \tag{2.1}
\]

with `s_k=s_0`. This is an integer polynomial with nonnegative coefficients and maximum degree
`k(m+n)=k(5L-2)`.

Cut the `2 x L x k` replica graph between spatial columns. A column state contains two spins in each
of `k` replicas, hence `4^k` possibilities. Let `ell(t)` count its `k` rung equalities and its `2k`
directed cyclic-replica equalities. If `v_j(t)` is the partial polynomial after `j` columns, then

\[
\begin{aligned}
v_1(t)&=q^{\ell(t)},\\
v_{j+1}(t)&=q^{\ell(t)}
 \sum_s q^{2k-d_H(s,t)}v_j(s),\\
H_{k,L}(q)&=\sum_t v_L(t).                                \tag{2.2}
\end{aligned}
\]

The horizontal kernel is the Kronecker power

\[
\begin{pmatrix}q&1\\1&q\end{pmatrix}^{\otimes 2k}.
\]

Applying one two-entry butterfly per bit gives (2.2) in `O(L k 4^k deg H)` packed additions rather
than a dense `4^k x4^k` multiply. The producer validates the identity against a literal dense
`ZZ[q]` construction on open `2x2` through trace four and reproduces the landed e238 open-`2x3`
seven-trace content digest exactly.

## 3. Why trace eight is underdetermined

Assume a positive full eight-mode subset-product spectrum. Write each mode ratio as
`z_i^2>0`, divide by the unique positive determinant centre, and set

\[
\tau_i=z_i+z_i^{-1},\qquad u_i=\tau_i^2\ge4,\qquad
Q(u)=\prod_{i=1}^8(u-u_i).
\]

If `L_0(tau)=2`, `L_1(tau)=tau`, and
`L_k=tau L_(k-1)-L_(k-2)`, then
`L_k(z+z^{-1})=z^k+z^{-k}` and the centred power traces factor as

\[
r_k=\prod_{i=1}^8 L_k(\tau_i).                           \tag{3.0}
\]

The recurrence gives

\[
\begin{array}{c|c}
k&\text{factor used in the norm row}\\ \hline
4&u^2-4u+2\\
5&L_5/L_1=u^2-5u+5\\
6&L_6/L_2=u^2-4u+1\\
7&L_7/L_1=u^3-7u^2+14u-7\\
8&u^4-8u^3+20u^2-16u+2\\
9&L_9/L_1=u^4-9u^3+27u^2-30u+9 .
\end{array}                                               \tag{3.0a}
\]

Because the mode count is even, products of these factors are the corresponding monic resultants
with `Q`; this gives `F4,...,F9`. The producer now derives (3.0a) from the recurrence, and the
clean-room verifier derives it again rather than copying the six factors. For physical `q>1`, the
raw trace polynomials have nonnegative coefficients and positive configuration sum, while
`q>0` and `q^2-1>0`; hence `r1`, `r2`, and every displayed target denominator are nonzero.

For eight modes write

\[
Q(u)=u^8+a_7u^7+a_6u^6+a_5u^5+a_4u^4+a_3u^3+a_2u^2+a_1u+a_0.
\]

As in e232/e238, positivity and determinant centering give

\[
Q(0)=r_1^2,\qquad Q(2)=r_2,\qquad Q(3)=r_3/r_1.
\]

These determine `a0,a2,a1` and leave `a7,...,a3`: five coefficient variables. The Lucas factors for
traces four through nine give six residual equations. Their exact degrees after substitution are
`2,2,2,3,4,3`.

At fixed `q`, `F4,...,F8` are five equations in five coefficient variables. With `q` also free they
describe a one-dimensional incidence; inequalities such as shifted-Hankel positivity do not lower
its complex dimension. Consequently a trace-eight finite-root claim would be invalid. `F9` is the
first additional basis-free scalar equation that can produce a finite `q` projection. This
dimension correction is load-bearing.

The leading forms of `F4,...,F8` are independent of all trace targets. Their exact Groebner basis has
35 elements and leading pure powers bounded by `(2,2,3,5,9)`. Enumerating monomials not divisible by
any leading monomial gives 96 standard monomials, agreeing with the complete-intersection Bezout
product

\[
2\cdot2\cdot2\cdot3\cdot4=96.                            \tag{3.1}
\]

Here is the filtered/projective step that promotes the leading calculation. Let
`K=Q(q)` and homogenize the five residuals in a new variable `x0`, preserving degrees
`2,2,2,3,4`. On `x0=0`, the equations are exactly the five recorded leading forms. Their
homogeneous ideal is zero-dimensional in affine five-space, so their only affine common zero is
the origin and they have no projective common zero. Therefore the projective closure of the full
five-equation affine scheme has no point at infinity. It cannot have a positive-dimensional
component: any such projective component would meet the hyperplane `x0=0`. Hence the five
hypersurfaces form a zero-dimensional complete intersection in `P^5`. Projective Bezout gives
scheme length equal to the product (3.1). For physical `q>1`, `q`, `q^2-1`, and every raw positive
trace in the recorded denominators are nonzero, so specialization preserves the same leading
system and no-infinity argument.

Thus the next exact object is a `96x96`, not `192x192`, multiplication matrix for `F9`.

## 4. Verification and resource scope

The producer uses exact integers, rational functions, and Groebner bases over `Q`; no floating-point
root or fitted critical value enters. Darwin peak RSS is measured by
`mach_task_basic_info.resident_size_max`. It fails before writing on every semantic, provenance,
CPU, or RSS gate.

The verifier imports neither the producer nor its helpers. It reverses the column-state bit layout,
rebuilds every integer trace coefficient, independently reconstructs all nine normalized rational
targets, and derives the abstract Lucas incidence and 96 standard monomials from fresh symbols. The
open-`2x3` digest is an inherited positive regression, not a rerun of its degree-971 elimination.
