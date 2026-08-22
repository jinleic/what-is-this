# Physical-parity projected Gaussianity of the `2 x 3` layer

**Status tags.** **[THEOREM]** means proved for every case in its stated scope;
**[LEMMA]** means proved here; **[COMPUTATION]** means an exact finite computation.

Script: `experiments/e40_parity_sector.py`  
Artifact: `results/spectral/parity_sector.json`  
Standalone independent test: `tests/test_parity_sector.py`

## 1. The parity-pair question

Let

\[
P=\prod_{v=1}^{6}X_v
\]

be physical spin flip.  Theorem S rules out the full `2 x 3` layer transfer
operator as one six-mode fermionic Gaussian, but leaves open a parity-projected
interpretation.  The precise escape route tested here is this spectral one:
the two physical spectra are, in one order or the other, the even- and
odd-cardinality halves of one product spectrum

\[
 \left\{\lambda_*\prod_{k\in S}u_k:S\subseteq\{1,\ldots,m\}\right\}.
\]

The physical sectors both have dimension 32, so this ansatz has `m=6`.  The
argument below is stated for all `m >= 3`.

## 2. Forced cross-sector relation

**[LEMMA] (parity-projected pair relation).**  Let `m >= 3`, let
`lambda_* > 0`, and let

\[
  1<u_{(1)}<u_{(2)}\leq\cdots\leq u_{(m)}.
\]

Thus there are no zero modes (`u_k=1`) and the smallest factor is simple.
Define the two multisets, with subset multiplicity retained,

\[
 E=\left\{\lambda_*\prod_{k\in S}u_k:|S|\ {
m even}\right\},\qquad
 O=\left\{\lambda_*\prod_{k\in S}u_k:|S|\ {
m odd}\right\}.
\]

If `e_0` is the minimum of `E` and `o_0 < o_1` are the first two entries of
`O` sorted with multiplicity, then

\[
 e_0=\lambda_*,\qquad o_0=\lambda_*u_{(1)},\qquad
 o_1=\lambda_*u_{(2)},
\]

and consequently

\[
 \boxed{\frac{o_0o_1}{e_0}=\lambda_*u_{(1)}u_{(2)}\in E.}
\]

*Proof.*  The empty subset is even and has product 1.  Every other even subset
has at least two elements.  For every pair `i != j`, sortedness gives
`u_i u_j >= u_(1) u_(2)`: among all distinct index pairs, the two smallest
factors minimize the product.  An even subset of size at least four has still
larger product because every omitted factor is greater than 1.  In particular,
no two-element even product can undercut `u_(1)u_(2)`, and no nonempty even
product can undercut 1.  Hence `min E=lambda_*`.

The singleton odd products are the factors themselves, so their two smallest
entries are `u_(1)` and `u_(2)`.  Every nonsingleton odd subset has at least
three elements.  Its product is at least
`u_(1)u_(2)u_(3)>u_(2)`, because all factors are greater than 1.  Thus no
three-or-more-element odd product can intervene below the second singleton.
The strict inequality `u_(1)<u_(2)` makes the first two odd entries distinct,
including when `u_(2)` has higher multiplicity.  Finally, the subset `{1,2}`
is even, so its product gives the boxed member of `E`.  []

### What the sector spectra themselves must certify

It is useful not to assume the factor hypotheses silently.  Start more weakly
with canonically oriented factors `u_k >= 1`.

* The empty subset gives `e_0=lambda_*`.  If any zero mode has `u_k=1`, its
  one-element subset gives `o_0=lambda_*` as well.  Therefore the strict,
  certified ordering `e_0<o_0` excludes every zero mode.
* Once zero modes are excluded, if the smallest factor occurs at least twice,
  two singleton states give the first two odd entries the same value.
  Therefore `o_0<o_1` certifies that the smallest factor is simple.
* Under these two strict inequalities, the lemma applies.  Thus the three
  relevant values `e_0,o_0,o_1` must have the ordered, pairwise-disjoint
  pattern `e_0<o_0<o_1`, and the boxed value is forced in the even sector.

This is the exact degeneracy boundary of the present criterion.  If `e_0` and
`o_0` are not certified distinct, a zero mode has not been excluded.  If the
first two odd entries are not certified distinct, simplicity of the smallest
factor has not been established.  A repeated positive smallest factor admits
additional multiplicity arguments, but it is outside the stated lemma; the
simple-factor inference cannot be made from overlapping first-two enclosures.
Likewise, overlapping rational enclosures are merely non-decisive and are
never interpreted as equality.  These limitations parallel the three-lowest-
eigenvalue caveat in Theorem S.

## 3. Why both physical assignments must be tested

**[LEMMA] (assignment symmetry).**  Identifying physical spin flip with a
fermionic `Z_2` grading has an unavoidable sign convention.  Either

\[
(P=+1,P=-1)=(E,O)
\]

or

\[
(P=+1,P=-1)=(O,E).
\]

There is no spectral reason to choose one convention in advance.  Therefore a
physical-parity obstruction must reject both assignments.

*Proof.*  If `F` is the fermionic parity involution, both `F` and `-F` commute
with the same Gaussian and exchange only the names of its two eigenspaces.
An identification of physical `P` with the grading can consequently be
`P=F` or `P=-F`.  These are exactly the two assignments displayed above.  []

For each assignment, the even-sector minimum must be the global product
minimum `lambda_*`.  Thus an assignment with `o_0<e_0` is already impossible.
When the required minimum ordering holds, absence of `o_0o_1/e_0` from the
assigned even sector supplies the cross-sector contradiction.

## 4. Exact rational physical-parity restrictions

Take `t=tanh(K*/2)=1/3`.  As in the construction behind Theorem S, put

\[
 q=\exp(2K)=\frac{1+t^2}{2t}=\frac53,
 \qquad (P_t)_{ab}=t^{\operatorname{Ham}(a,b)}.
\]

If `b_a` is the intralayer Ising bond sum of computational state `a` and
`epsilon=b_a mod 2`, then, up to one common positive scalar, the symmetrized
transfer operator is

\[
 R=P_t\,\operatorname{diag}\left(q^{(b_a-\epsilon)/2}\right)P_t.
\]

Every entry of `R` is rational and `R` is symmetric positive definite.
The symmetrizing similarity is generated by `A=sum_v X_v`; since physical
`P` commutes with `A`, this similarity leaves the computational-basis matrix
of `P` unchanged.


Let `c(a)` be bitwise complement.  In the computational basis, physical `P`
is the permutation `e_a -> e_c(a)`, hence a permutation matrix (in particular,
a signed permutation matrix with all nonzero signs `+1`).  It is fixed-point
free.  Hamming distance obeys
`Ham(c(a),c(b))=Ham(a,b)`, and global spin reversal leaves every bond product
unchanged, so `b_c(a)=b_a`.  It follows algebraically that `PRP=R`, equivalently
`PR=RP`.  **[COMPUTATION]** The script also verifies the permutation property,
involution property, absence of fixed points, and all 4096 rational matrix
equalities `R[c(a),c(b)]=R[a,b]`.

Choose one representative `a=0,...,31` from each orbit `{a,c(a)}` and set

\[
 b_a^{\pm}=e_a\pm e_{c(a)}.
\]

These basis vectors have entries only in `{0,+1,-1}`.  If `B_+` and `B_-`
are the matrices with these columns, then

\[
 B_s^T B_s=2I,\qquad B_+^TB_-=0.
\]

Because the two subspaces are invariant, the coordinate restriction is

\[
 R_s=(B_s^TB_s)^{-1}B_s^TRB_s=\frac12B_s^TRB_s,
\]

an exact rational `32 x 32` symmetric matrix.  No irrational `1/sqrt(2)`
normalization is needed.  More explicitly, for a rational shift `sigma`,

\[
 B_s^T(R-\sigma I)B_s=2(R_s-\sigma I).
\]

The left side is the restricted quadratic form in the unnormalized basis, and
multiplication by the positive scalar 2 changes no inertia.  Hence exact
inertia counts for `R_s-sigma I` are exactly the spectral counts in the
physical sector.  This is the promised congruence-invariance argument: an
unnormalized rational basis preserves the inertia of the restricted form, and
its simple Gram matrix supplies the rational coordinate operator.

## 5. Exact Sylvester-inertia certificate

**[LEMMA].**  Exact symmetric elimination of a rational symmetric matrix
`R_s-sigma I`, when all pivots are nonzero, factors it by invertible rational
congruences into a diagonal matrix.  Each elimination step subtracts a rational
multiple of one row and the same multiple of the corresponding column, so it
is an invertible congruence.  Therefore the number of negative final pivots is
exactly the number of eigenvalues of `R_s` below `sigma`.

The implementation never assigns a count to a shift with a zero pivot.
Rational bisection isolates eigenvalues; the final forced interval is rounded
outward if an endpoint has a zero pivot.  Equal counts at its two endpoints
prove that the whole interval contains no target-sector eigenvalue.

**[COMPUTATION] (exact finite certificate).**  The script isolates the three
lowest eigenvalues of both physical sectors (more than the lemma needs) with
relative rational width at most `10^-8`.  The following decimals are display
only; every endpoint stored in the JSON is an exact fraction.

| physical sector | lowest enclosure | second enclosure | third enclosure |
|---|---:|---:|---:|
| `P=+1` | `[0.004537261830673439, 0.004537261859178012]` | `[0.031113904864936665, 0.031113905092973244]` | `[0.04847791128889378, 0.04847791174496694]` |
| `P=-1` | `[0.009649610643311775, 0.009649610695319496]` | `[0.016327359583079192, 0.016327359687094634]` | `[0.02005460503812154, 0.020054605142136982]` |

All displayed enclosures within each sector are disjoint.  The two assignment
checks are:

| physical assignment | forced interval in assigned `E` | target-sector counts | exact conclusion |
|---|---:|---:|---|
| `P+ -> E`, `P- -> O` | `[0.03472417235327969, 0.03472417297979301]` | `2, 2` in `P+` | no eigenvalue in the outward interval |
| `P+ -> O`, `P- -> E` | `[0.014629806051694562, 0.014629806329676049]` | `1, 1` in `P-` | no eigenvalue in the outward interval |

For the first assignment, the exact enclosure order is
`e_0<o_0<o_1`; it excludes zero modes, certifies a simple smallest factor, and
activates the lemma.  Equal target counts then contradict its forced relation.

For the swapped assignment, the exact sector minima instead obey
`o_0<e_0`.  This already contradicts `min E=lambda_*<min O` under the
no-zero-mode product ansatz.  The script nevertheless runs the requested
swapped cross-product construction and finds its interval empty as recorded in
the second row.  That additional absence is an exact computation, but the
minimum-order contradiction is the logically prior rejection because the
lemma's oriented-even hypothesis has already failed.

## 6. One-dimensional Gaussian control

**[COMPUTATION].**  The identical exact construction was run on the open
six-site chain, whose transfer operator is the free-fermionic control.  For
`P+ -> E` the forced interval is

\[
[0.0403862550979781,0.04038625586843526]
\]

and target-sector inertia changes `1 -> 2`.  For the swapped assignment the
interval is

\[
[0.015592306756822173,0.015592307039393873]
\]

and inertia again changes `1 -> 2`.  Thus each positive-width window contains
at least one eigenvalue.  This is only a consistency check: differing counts
prove the presence of *some* eigenvalue in the window, not exact equality to
the forced algebraic value.

## 7. Theorem S′ and scope

**[THEOREM S′].**  At the exact coupling
`t=tanh(K*/2)=1/3` (`exp(2K)=5/3`), the `2 x 3` open-layer Ising transfer
operator cannot realize a fermionic Gaussian parity pair whose grading is the
physical spin flip `P=product_v X_v`, under the positive product-spectrum and
distinct-lowest-value hypotheses of the lemma.  Equivalently, its two physical
`P`-sector spectra cannot be the even- and odd-subset-product halves of one
six-mode Gaussian spectrum, in either assignment.

*Proof.*  Exact commutation gives the two rational physical restrictions.  The
assignment `P+ -> E`, `P- -> O` has the three required ordered and disjoint
lowest values, but exact equal inertia counts show that its forced even-sector
value is absent.  The swapped assignment violates the necessary ordering of
the even and odd minima (and its mechanically constructed cross-value window
is also exactly empty).  The assignment-symmetry lemma exhausts the two ways
physical `P` can label Gaussian parity.  []

**Scope, stated narrowly.**

* This is a theorem about one exact coupling and the finite `2 x 3` open layer;
  it is not an all-coupling or all-size theorem.
* It covers only the **physical** parity operator `P=product_v X_v`.  A
  hypothetical Gaussian structure whose fermion parity is a different
  commuting `Z_2` symmetry is not excluded.
* It does not exclude restrictions of a Gaussian to a non-parity invariant
  subspace, nor arbitrary submultisets of a larger Gaussian spectrum.
* It assumes the positive real product-spectrum form with no zero modes and a
  simple smallest factor, as certified by the relevant lowest-sector ordering
  in the viable assignment.  Degenerate cases not satisfying those spectral
  hypotheses are not silently claimed.
* The chain's differing counts are a control only; no exact-presence theorem is
  inferred from a finite-width interval.
