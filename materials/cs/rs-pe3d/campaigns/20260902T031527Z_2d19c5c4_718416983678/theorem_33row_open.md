# Tied three-row Vandermonde product circuits through the ambient cap
## gate H-33ROW-OPEN, run `20260902T031527Z_2d19c5c4_718416983678`

Agent `RsPe3dH2`, 2026-09-02 UTC. Preregistration:
`prereg/H_33ROW_OPEN_PREREG_2026-09-02.md`, path-scoped source commit
`c1ff0216e7e0bdddf6ab89d082db3fcef6e1aa86`, SHA-256
`adfccb116b4ebde7ed6847fbbb7afca7424ef444694a566bc224fe9f3f8e2dc8`.
The in-run `pre_statement.md` is byte-identical to that source (both `cmp` and
SHA-256 checked); the preregistration was not amended.

Let

\[
a_x=(1,x,x^2)^{\mathsf T},\qquad b_y=(1,y,y^2)^{\mathsf T},\qquad
h(x,y)=a_x\otimes b_y=(x^iy^j)_{0\le i,j<3}.
\]

Each axis has at least four distinct evaluation points.  Thus each factor
represents the rank-three uniform matroid, its spark is four, the product spark
is four, and the product ambient rank is nine.  The results below separate proved
statements from exhaustive finite censuses.  In particular, a rank/minor
predicate is exact but is not called a structural parameterization.

## 1. Two elementary restrictions

If a circuit of size greater than four had four cells in one matrix row, those
four cells would be a proper B-fiber circuit.  Hence every row and column has
degree at most three.  If a nonfiber circuit used at most three values on one
axis, the corresponding factor columns would be independent.  Dual isolation
of those factor columns would make each opposite-axis group a relation, hence
would expose a four-fiber proper subcircuit.  Therefore every nonfiber circuit
has both projections of size at least four.

For later use, let a full-support relation on a support `S` have coefficients
`gamma_xy`, and for each represented A-value put

\[
z_x=\sum_{y:(x,y)\in S}\gamma_{xy}b_y\in F^3.
\tag{1}
\]

A row has at most three cells, so factor spark four and full support imply
`z_x != 0`.  The tensor relation is

\[
\sum_x a_xz_x^{\mathsf T}=0.
\tag{2}
\]

If the A-projection has `r` points, every coordinate vector of the family
`(z_x)` belongs to the nullspace of the `3 x r` Vandermonde matrix.  For
`r >= 4`, with

\[
w_x=\left(\prod_{x'\ne x}(x-x')\right)^{-1},
\]

that nullspace is

\[
\{(w_xq(x))_x:\deg q\le r-4\}.
\tag{3}
\]

Indeed, for every polynomial `f` of degree at most `r-2`, Lagrange
interpolation says that `sum_x w_x f(x)` is the coefficient of degree `r-1`
in its interpolant, hence is zero.  Taking `f(x)=x^kq(x)` for `k=0,1,2`
shows that every vector displayed in (3) is in the Vandermonde nullspace; both
spaces have dimension `r-3`, proving equality.

## 2. Complete layers four and five

**Theorem 2.1 (size four).**  The only size-four circuits are factor-circuit
fibers.  On `|X|=|Y|=n` points their number is

\[
\boxed{N_4=2n\binom n4}.
\tag{4}
\]

Every four-subset in one fiber transports a factor circuit and is therefore a
circuit.  Conversely, let a four-set not lie in a fiber and fix one of its
cells `(x_0,y_0)`.  Among the other three cells, place every cell with
`y=y_0` into a set `R`; nonfiberhood gives `|R| <= 2`, and their distinct
`x`-values are roots of a polynomial `f(x)` of degree at most two with
`f(x_0) != 0`.  If `R` is empty, add to `f` the `x`-value of any other cell
whose `x` differs from `x_0` (one exists because the set is not an A-fiber).
At most two cells remain; their `y`-values, all different from `y_0`, are
roots of a polynomial `g(y)` of degree at most two with `g(y_0) != 0`.
Thus `f(x)g(y)` is a bidegree-`(2,2)` form that vanishes on the other three
cells but not on the fixed cell.  Repeating this for every cell gives four
dual isolators, so the four product columns are independent.  Hence only the
two fiber orientations are circuits; they are disjoint and yield (4).

**Theorem 2.2 (size five).**  There are no size-five circuits.

*Proof.*  A factor has no size-five circuit, so there is no fiber case.  In a
mixed case both projections are at least four.  Suppose first that the
A-projection has four points.  Its Vandermonde nullspace is one-dimensional
and has full support.  Equations (2)--(3) therefore say that all four nonzero
vectors `z_x` are projectively equal.  The only possible row-degree sequence
at size five is `(2,1,1,1)`.  The three singleton rows consequently use the
same B-value `y_0`, because distinct vectors `(1,y,y^2)` are not parallel.  The
double row would have to express a multiple of `b_{y_0}` as a combination of
two distinct B-columns with both coefficients nonzero.  This contradicts
independence of every three distinct B-columns (and still contradicts full
support if `y_0` is one of the two).  The case of B-projection four is
symmetric.  The only remaining profile is all-distinct `(5,5)`, but applying
the full-support relation in both tensor directions gives the standard
rank-sum bound `3+3 <= 5`, a contradiction.  Thus no case remains.  ∎

## 3. Complete size-six classification

### 3.1 No hidden repeated-index profile

**Lemma 3.1.**  Every size-six circuit has profile `(4,4)` or `(6,6)`.

*Proof.*  A size-six fiber contains a four-fiber proper subcircuit, so the
circuit is nonfiber.  Section 1 gives `4 <= r <= 6` for its A-projection.

If `r=4`, equations (2)--(3) again make all `z_x` projectively equal.  The two
possible row-degree sequences are `(3,1,1,1)` and `(2,2,1,1)`.  Every
singleton row must use one common B-value `y_0`.  A degree-two row cannot
express a multiple of `b_{y_0}` with two nonzero coefficients, by independence
of every three B-columns.  Hence `(2,2,1,1)` is impossible.  In the other
case, the degree-three row expresses `b_{y_0}` using three B-columns.  Full
support and three-column independence force those three values to be distinct
from `y_0`; the four involved B-columns carry their unique four-circuit
relation.  The support is exactly three cells in column `y_0` and three cells
in the remaining row, with their intersection omitted.  It has profile
`(4,4)` and is a crossing.

If `r=5`, (3) has `q` of degree at most one.  Thus the projective points
`[z_x]` lie in a projective line in `P^2`.  The row-degree sequence is
`(2,1,1,1,1)`.  The four singleton rows give points `[1:y:y^2]` where that
line meets the nonsingular conic.  A line meets this conic in at most two
distinct points, so the singleton rows use at most two B-values.  The double
row adds at most two more, hence the B-projection has size at most four.  It is
at least four by Section 1, so it is four.  Applying the already proved
`r=4` case after swapping the axes would force the A-projection also to have
size four, a contradiction.

If `r=6`, every row has degree one.  Projection sizes four and five have just
been excluded on the other axis, so the B-projection also has size six.  These
are all possibilities.  ∎

### 3.2 Crossing channel

For four-subsets `T subset X`, `Z subset Y` and centers `x_0 in T`,
`y_0 in Z`, define

\[
S(T,x_0,Z,y_0)=
((T\setminus\{x_0\})\times\{y_0\})\cup
(\{x_0\}\times(Z\setminus\{y_0\})).
\tag{5}
\]

The two factor circuit relations can be scaled so that their coefficients at
the omitted center cancel.  The resulting relation has exactly support (5).
Every five-subset is independent by Theorem 2.2, so (5) is a circuit.
Conversely, the `r=4` part of Lemma 3.1 recovered exactly this shape and its
two four-circuit projections.  The row and column of degree three recover the
centers and the two projection sets, so the parameterization is injective.
Consequently

\[
\boxed{N_{44}^{(6)}=16\binom{|X|}{4}\binom{|Y|}{4}}.
\tag{6}
\]

This channel and its count are field-free once the two factor circuit spectra
are fixed.

### 3.3 All-distinct Möbius channel

**Theorem 3.2 (M6).**  Six cells with distinct `x`- and `y`-coordinates form a
circuit if and only if they lie on the graph of a unique nonconstant
linear-fractional map

\[
x=M(y)=\frac{t(y)}{r(y)},\qquad M\in\mathrm{PGL}(2,F),
\tag{7}
\]

with no selected pole.

*Proof.*  Order the cells `(x_j,y_j)` and let `lambda_j` be their nonzero
relation coefficients.  With

\[
w_j=\left(\prod_{k\ne j}(y_j-y_k)\right)^{-1},
\]

(3), applied to the three A-coordinate slices, gives quadratics
`q_0,q_1,q_2` such that

\[
\lambda_j=w_jq_0(y_j),\quad
x_jq_0(y_j)=q_1(y_j),\quad
x_j^2q_0(y_j)=q_2(y_j).
\tag{8}
\]

No `q_0(y_j)` vanishes.  The degree-at-most-four polynomial
`q_1^2-q_0q_2` has all six `y_j` as roots, hence is zero.  Unique
factorization in `F[y]` then gives, up to nonzero units,

\[
q_0=sr^2,\qquad q_1=srt,\qquad q_2=st^2,
\tag{9}
\]

with `gcd(r,t)=1`.  Since all three `q` have degree at most two, both `r` and
`t` have degree at most one.  They are linearly independent because the
selected `x_j` are distinct.  Equations (8)--(9) give (7), and no selected
`r(y_j)` is zero.

Conversely, a nondegenerate `(1,1)` graph is isomorphic to `P^1`, and
`O(2,2)` restricts to `O(4)`.  Equivalently, its bilinear equation multiplies
`H^0(O(1,1))` into the four-dimensional kernel of the nine-dimensional
restriction map, whose image is the complete five-dimensional degree-four
space.  Evaluation at any five distinct parameter values is nonsingular
Vandermonde evaluation, while six values have rank five.  Thus six graph
points form a circuit.  Three paired points determine a unique projective
linear transformation, proving uniqueness.  ∎

Put

\[
k_M=|\{y\in Y:M(y)\in X\}|.
\]

Uniqueness makes the graph classes disjoint, so

\[
\boxed{N_{66}^{(6)}(X,Y)=\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{6}}.
\tag{10}
\]

Combining Lemma 3.1, (6), and (10) gives the complete size-six count

\[
\boxed{N_6=16\binom{|X|}{4}\binom{|Y|}{4}
      +\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{6}}.
\tag{11}
\]

The first term is combinatorial; every field dependence in the all-distinct
term is exactly the induced PGL incidence profile.

## 4. Certified higher-layer statements

### 4.1 The entire size-seven layer is empty

**Theorem 4.1.**  There are no size-seven circuits.

*Proof.*  A size-seven fiber contains a four-fiber proper subcircuit, so the
circuit is nonfiber.  Write `r` for its A-projection size and continue to use
(1)--(3).  The projection bounds give `4 <= r <= 7`.

If `r=4`, all projective vectors `[z_x]` are equal.  The possible row-degree
sequences are `(3,2,1,1)` and `(2,2,2,1)`.  A singleton row identifies the
common vector with one conic point `[b_{y_0}]`, but a degree-two row cannot
express that point using two distinct B-columns with both coefficients
nonzero.  Both sequences are impossible.

If `r=5`, write the degree-at-most-one vector polynomial in (3) as
`Q(x)=u+xv`.  The row-degree sequence is either `(3,1,1,1,1)` or
`(2,2,1,1,1)`.  If `u,v` are linearly independent, the projective map
`x mapsto [u+xv]` is injective, so the four or three singleton rows give that
many distinct points where a line meets the conic `[1:y:y^2]`; a line has at
most two such intersections.  If `u,v` are dependent, all singleton rows use
one `y_0`.  Four of them violate the degree-three ceiling, while in the
three-singleton case either degree-two row would have to express
`b_{y_0}` using two distinct B-columns, again impossible.  Thus `r=5` is
impossible.

If `r=6`, the row-degree sequence is `(2,1,1,1,1,1)`.  In (3) write
`Q=(q_0,q_1,q_2)` with each `q_i` quadratic.  At each of the five singleton
rows, `[Q(x)]=[1:y:y^2]`, so the degree-at-most-four polynomial
`q_1^2-q_0q_2` has five distinct roots and vanishes identically.  Hence the
value at the double row also lies on the same nonsingular conic.  But that
value is a linear combination, with both coefficients nonzero, of the two
distinct conic points belonging to the double row.  Their secant line meets a
nonsingular conic only at its two endpoints, so no third projective conic
point can be such a combination.  This is impossible.

Finally, if `r=7`, every row is a singleton.  The symmetric versions of the
three cases just proved exclude B-projection sizes four, five, and six, so the
support is all-distinct.  For seven distinct `y_j`, (3) gives
`q_0,q_1,q_2` of degree at most three.  The polynomial
`q_1^2-q_0q_2` has degree at most six and seven roots, hence vanishes.  The
factorization (9) and the degree-three ceiling force
`deg r,deg t <= 1`.  All seven points lie on one Möbius graph and contain a
size-six circuit, contrary to minimality.  No case remains.  ∎

### 4.2 A proved size-eight rational-degree-two family

Let `R(y)=P(y)/Q(y)` be a reduced rational map of exact degree two, and choose
eight distinct finite parameter values at which `Q` is nonzero.  Its graph is
a divisor of bidegree `(1,2)` (up to the axis convention).  It is isomorphic
to `P^1`, and

\[
O(2,2)|_{\operatorname{graph}R}\simeq O_{P^1}(6).
\]

The graph equation leaves a two-dimensional kernel in
`H^0(O(2,2))`, so restriction is onto the seven-dimensional complete
degree-six space.  Seven distinct parameter values evaluate independently;
eight have rank seven.  Hence every such eight-point graph support is a
circuit.  This statement allows repeated images `R(y)`; the dedicated matching
subfamily additionally requires the eight images to be distinct.

Two reduced degree-at-most-two maps agreeing at five finite paired points are
equal: after cross multiplication their difference has degree at most four
and five roots.  Therefore this family is injectively parameterized.  If

\[
D_R(X,Y)=\{y\in Y:Q(y)\ne0,\ R(y)\in X\},
\]

then its exact incidence count is

\[
\boxed{N_{R_2}^{(8)}(X,Y)=
 \sum_{\substack{R\text{ reduced rational}\\\deg R=2}}
 \binom{|D_R(X,Y)|}{8}}.
\tag{12}
\]

No exhaustion claim follows.  Indeed, the dedicated all-distinct sweeps below
found many size-eight circuits and no member of this subfamily for those
particular finite point sets.  Complete intersections/pencils of `(2,2)`
curves are a genuine residual mechanism.

### 4.3 Exact size-nine and size-ten predicates

For nine points, the following is exact:

\[
S\text{ is a circuit}\iff
\det E_S=0\quad\text{and}\quad
\operatorname{rank}E_{S\setminus\{s\}}=8\ \text{for every }s\in S.
\tag{13}
\]

The determinant condition says precisely that a nonzero bidegree-`(2,2)` form
vanishes on all nine points.  The deletion clauses are indispensable
minimality conditions; a common curve alone overaccepts Möbius graphs.

Every ten-set is dependent in the nine-dimensional ambient space.  Thus

\[
S\text{ is a size-ten circuit}\iff
\det E_{S\setminus\{s\}}\ne0\ \text{for every }s\in S.
\tag{14}
\]

Equations (13)--(14) are exact decision predicates, not structural
parameterizations or closed counting formulas.

## 5. Exact finite evidence

All decisions used Python integer arithmetic modulo `p`; no floating point was
used.  Every full-grid subset was evaluated by two independently written
finite-field rank routines, which agreed on every subset.  The circuit test
was rank `m-1` together with independence of every deletion (and, at size ten,
all nine-deletion determinants nonzero).  Thus the two rank implementations
induce the same overall and per-signature circuit populations.

### 5.1 Full grids

| grid | p | size 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `4 x 4` | 7 | 8 | 0 | 16 | 0 | 66 | 0 | 16 |
| `4 x 4` | 11 | 8 | 0 | 16 | 0 | 66 | 0 | 16 |
| `4 x 4` | 13 | 8 | 0 | 16 | 0 | 66 | 0 | 16 |
| `5 x 5` | 7 | 50 | 0 | 400 | 0 | 6030 | 36378 | 181960 |
| `5 x 5` | 11 | 50 | 0 | 400 | 0 | 5950 | -- | -- |
| `5 x 5` | 13 | 50 | 0 | 400 | 0 | 5814 | -- | -- |

The dashes are deliberately outside the preregistered mandatory matrix.
Every size-six grid circuit is a crossing: the direct profile-`(4,4)` set is
exactly the constructed set, with counts `16` and `400` as (6) predicts.
There is no room for the all-distinct channel at `n<6`.

At `n=4`, the size-eight population has three label-invariant signatures with
counts `18,24,24`; size ten has one signature with count `16`.  At `n=5`,
size eight has seven signatures.  Its profile counts are

| p | `(4,4)` | `(4,5)` | `(5,4)` | `(5,5)` |
|---:|---:|---:|---:|---:|
| 7 | 1650 | 1800 | 1800 | 780 |
| 11 | 1650 | 1800 | 1800 | 700 |
| 13 | 1650 | 1800 | 1800 | 564 |

The first three profile support sets are byte-for-byte equal across all three
primes.  The `(5,5)` sets are not: pairwise intersections are `236`, while the
left/right-only populations are nonzero.  This is direct evidence against a
prime-independent extrapolation.  At `n=5,p=7`, size nine has 11 signatures
and profile counts `(4,5):7200`, `(5,4):7200`, `(5,5):21978`; size ten has 26
signatures and profile counts `(4,4):400`, `(4,5):20400`, `(5,4):20400`,
`(5,5):140760`.  Full signature strings and counts are retained in
`controls_results.json`.

### 5.2 Dedicated all-distinct sweeps

For every size-six matching, the direct circuit set equaled the independently
constructed PGL graph set:

| n | p=7 | p=11 | p=13 |
|---:|---:|---:|---:|
| 6 | 12 | 4 | 2 |
| 7 | 588 | 72 | 42 |

At `n=7`, all `7!=5040` full matchings at each prime were tested and the
size-seven circuit count was zero, as Section 4.1 proves.

At `n=8`, all `8!=40320` matchings were tested at each prime.  There were
`560` size-eight circuits over GF(11) and `416` over GF(13).  The common-curve
dependency predicate accepted `2264` and `1380` matchings respectively, so
minimality matters.  The exact rational-degree-two fitter found zero members
on these particular point sets; hence all `560/416` circuits are residual to
(12), not evidence that (12) is empty in general.  A separate GF(101)
`x=y^2` plant is a genuine rational-degree-two size-eight circuit.

At `n=9`, all `9!=362880` matchings were tested at each prime.  Predicate
(13) yielded `13248` circuits over GF(11) and `14214` over GF(13), while the
common-curve condition alone accepted `71160` and `50210`.  No PGL or simple
curve count is inferred.

### 5.3 Controls, runtime, and audit correction

The final corrected run exercises both A- and B-fiber ACCEPTs; a five-set
containing a fiber and an independent generic five-set REJECT; crossing and
Möbius ACCEPTs; inserted-crossing and explicit non-Möbius REJECTs; a
rational-degree-two fit/circuit ACCEPT and a degree-three non-fit/independence
REJECT; dependent-but-nonminimal Möbius-curve REJECT; ambient-ten ACCEPT and
nonminimal REJECT; a corrupted actual Kronecker column REJECT; and the
nine-coordinate Kronecker-versus-six-coordinate-concatenation guard.

All `143` final assertions passed.  The final invocation used `618.124829`
CPU seconds and `706.9142427500337` wall seconds against hard caps of
`7200/8000` seconds, under `nice -n 15`, with all four thread caps equal to
one and both bytecode guards active.  Including the preserved first pass and
the stopped partial rerun, cumulative use was `1256.758397` CPU seconds and
`1434.3297473749844` wall seconds, still below the registered caps.

The first complete computational pass was numerically correct but its plant
battery used one generic `fiber_accept` instead of separately naming both
fiber orientations and omitted the preregistered generic-five REJECT.  A proof
audit also found that the rational-degree-two ACCEPT checked circuitness but
not the fitter and had no explicit non-degree-two REJECT, and that the
Möbius channel lacked a named non-Möbius six-set REJECT.  That first complete
output is preserved as `controls_results_preplant_audit.json`, with its final
checkpoint as `checkpoint_preplant_audit.json`.  A superseded rerun was stopped
as soon as the latter two omissions were found; its partial checkpoint is
preserved as `checkpoint_cancelled_plant_audit.json`.  The source was corrected
and the complete matrix, not merely the affected plants, was rerun.  No
parameter, expected census count, or analytic claim was changed, and the
preregistration was not amended.  `defect_log.json` records this provenance.

## 6. Exact boundary of certification and verdict

Certified analytically:

1. size four consists exactly of fibers, with (4);
2. size five is empty;
3. size six consists exactly of crossings and unique Möbius graphs, with the
   complete count (11);
4. the entire size-seven circuit layer is empty;
5. the rational-degree-two size-eight family and its injective incidence count
   (12);
6. the exact curve-plus-minimality size-nine predicate (13); and
7. the exact ambient deletion-minor size-ten predicate (14).

The **first unresolved structural layer is size eight outside the
rational-degree-two family**.  Already the all-distinct `n=8` sweeps contain
`560/416` residual circuits, and the full `5 x 5` grids show a
field-dependent `(5,5)` population.  Sizes nine and ten have exact predicates
but no structural parameterization or closed incidence count.

Therefore not every requested layer 4--10 is analytically classified.  Under
the preregistered rule the whole-gate verdict is **FROZEN-INCONCLUSIVE**, while
the numbered statements above remain certified subresults.  The strongest
exhaustive ranges are all subsets of `4 x 4` at sizes 4--10 over
GF(7), GF(11), GF(13); all subsets of `5 x 5` at sizes 4--8 over those primes;
all `5 x 5` size-nine and size-ten subsets over GF(7); and the dedicated
matching ranges stated in Section 5.2.
