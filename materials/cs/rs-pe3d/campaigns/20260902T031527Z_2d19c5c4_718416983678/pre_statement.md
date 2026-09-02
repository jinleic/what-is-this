# PREREG-RSPE3D-13-33ROW-OPEN — gate **H-33ROW-OPEN**

Written 2026-09-02 UTC by agent `RsPe3dH2`, after
H-3ROW-DETERMINANTAL was frozen, closed FROZEN-CERTIFIED, independently
verified by Main, and committed. This file is path-scoped committed before any
new parameter-dependent computation. Binding order: preregistration -> commit
-> `campaign.py init` -> byte-identical in-run copy with source commit and
SHA-256 -> controls/censuses -> analytic record -> `freeze` -> exactly one
`close --verdict`. Write only under `cs/rs-pe3d/`; do not edit any frozen run
or root ledger.

The candidate claims below were derived analytically before this registration
but are not authority. Counts outside proved constructions remain census data.
A rank/minor predicate is exact but is not called a structural closed form
unless its graph/geometry and injectivity are proved.

## 1. Setting

Let both factors be 3-row Vandermonde/GRS matrices

\[
a_x=(1,x,x^2)^{\mathsf T},\qquad
b_y=(1,y,y^2)^{\mathsf T}
\]

at distinct evaluation points. Product columns are the nine bidegree-$(2,2)$
monomials

\[
h(x,y)=(x^i y^j)_{0\le i,j<3}\in F^9.
\]

Thus $d_A=d_B=d=4$, the ambient dimension is nine, and circuit sizes lie in
$4,\ldots,10$. For a circuit larger than four, every row and column degree is
at most three. Every nonfiber circuit has both projections at least four.

For a support $S$, write $E_S=[h(x,y)]_{(x,y)\in S}$. The exact fallback
predicate is

\[
S\text{ circuit}\iff \operatorname{rank}E_S=|S|-1
\text{ and every one-column deletion has rank }|S|-1.
\tag{C}
\]

At size ten, dependence is automatic and (C) is exactly the nonvanishing of
all ten $9\times9$ deletion determinants.

## 2. P1 — layers 4 and 5

1. **Size 4:** reproduce T-DGE. The only circuits are A- or B-factor four-fibers,
   count $2n\binom n4$ when both index sets have size $n$.
2. **Size 5:** prove there are no circuits for pure 3-row Vandermonde factors.
   The $d+1$ theorem reduces this layer to fixed-coordinate lifts of factor
   size-five circuits and nonfiber two-factor channels. A rank-three MDS factor
   has only four-element circuits, so the lift term is empty. Thm-CRIT excludes
   the tied all-distinct channel for $d\ge4$, while the crossing channel has
   size $d_A+d_B-2=6$. Repeated mixed profiles are excluded by the certified
   projection/rank argument. Direct census must give the empty set.

## 3. P2 — size six

### 3.1 Crossing channel

The complete Theorem-X converse gives exactly the profile-$(4,4)$ crossing
family, with injective count

\[
\boxed{N_{44}^{(6)}=16\binom{|X|}{4}\binom{|Y|}{4}}.
\tag{X6}
\]

Every constructed support must be verified as a circuit and the measured
profile-$(4,4)$ set must equal it.

### 3.2 All-distinct Möbius channel

Preregister the following candidate theorem and prove it independently.
For six distinct B values $y_j$, set

\[
w_j=\left(\prod_{k\ne j}(y_j-y_k)\right)^{-1}.
\]

The kernel of the $3\times6$ B Vandermonde matrix is

\[
\{(w_jq(y_j))_j:\deg q\le2\}.
\]

A product relation yields quadratics $q_0,q_1,q_2$ with

\[
\lambda_j=w_jq_0(y_j),\qquad
x_jq_0(y_j)=q_1(y_j),\qquad
x_j^2q_0(y_j)=q_2(y_j).
\]

For a circuit $q_0$ is nonzero on all six nodes. The polynomial
$q_1^2-q_0q_2$ has degree at most four and six roots, hence vanishes
identically. Prove the UFD factorization forces a nonconstant linear-fractional
map $x=M(y)$; conversely restriction of bidegree-$(2,2)$ forms to a
nondegenerate $(1,1)$ graph is the complete degree-four univariate space, so
six graph points have rank five and all five-subsets are independent.

Therefore the candidate exact theorem/count is

\[
\boxed{S\text{ all-distinct size 6 circuit}\iff S\subset\operatorname{graph}M
\text{ for a unique }M\in\mathrm{PGL}(2,F)},
\tag{M6}
\]

\[
\boxed{N_{66}^{(6)}(X,Y)=\sum_{M\in\mathrm{PGL}(2,F)}\binom{k_M}{6}},
\qquad k_M=|\{y\in Y:M(y)\in X\}|.
\tag{MC6}
\]

Verify (M6) as support sets on dedicated all-distinct sweeps at $n=6,7$ over
at least three primes.

### 3.3 Other profiles

Exhaust every measured profile and component/degree signature. A six-set
satisfying the degree/projection bounds is a circuit iff its six columns have
rank five, because P1 leaves no dependent proper subset. Record the explicit
condition that every $6\times6$ maximal minor of the $9\times6$ evaluation
matrix vanishes. Promote a further field-free template only with an analytic
parameterization, injectivity, and exact count. Do not infer absence or a
formula from two/three-prime stability alone.

## 4. P3 — sizes 7 through 10 and the curve hierarchy

For every measured support record profile, label-invariant component/degree
signature, exact rank, and deletion-minor predicate. Direct and independently
recomputed minor-predicate sets must agree overall and per signature.

For an all-distinct $k$-matching, dependency means its evaluation functionals
on

\[
H^0(\mathbb P^1\times\mathbb P^1,\mathcal O(2,2))
=\operatorname{span}\{x^iy^j:0\le i,j<3\}
\]

are dependent. State the following only as far as proved:

- size 7: extend the Vandermonde-kernel argument. If dependence forces a
  Möbius graph, the support contains a size-six circuit and is not a circuit;
  certify emptiness of the all-distinct size-seven circuit channel if the proof
  closes;
- size 8: a reduced rational degree-two graph
  $x=P_2(y)/Q_2(y)$ pulls $\mathcal O(2,2)$ back to degree six, so eight
  distinct graph points are a candidate circuit when every seven-subset is
  independent. Prove this subfamily and uniqueness from five paired points if
  possible. Do not assert it exhausts size-eight matchings: a pencil/complete
  intersection of two $(2,2)$ curves may give additional rank-seven sets;
- size 9: an all-distinct matching is a circuit exactly when the $9\times9$
  bidegree evaluation determinant vanishes and every $8$-deletion has rank
  eight. This says the nine points lie on a nonzero $(2,2)$ curve with no
  dependent proper subset; do not turn it into a PGL count;
- size 10: dependence is ambient; circuitness is exactly ten nonzero
  nine-deletion determinants.

A “common curve” condition without the deletion/minimality clauses is the
preregistered wrong predicate and must overaccept a planted Möbius graph
containing a size-six circuit.

## 5. Exact computation matrix

All arithmetic is exact modulo $p$, stdlib only, with true nested Kronecker
columns.

### Mandatory full-grid sweeps

1. $n=4$, $p\in\{7,11,13\}$: exhaust every subset at sizes 4--10.
2. $n=5$, $p\in\{7,11,13\}$: exhaust sizes 4--8.
3. $n=5$, $p=7$: additionally exhaust sizes 9--10. This stage is mandatory
   unless the hard cap is reached; if it cannot complete, preserve the partial
   artifact and close inconclusive with the last completed size.

The $n=5$, $p=11,13$ sizes 9--10 are out of the mandatory matrix; do not infer
prime stability for those layers.

### Dedicated all-distinct sweeps

1. At $n=6,7$ and $p\in\{7,11,13\}$, exhaust all size-six matchings and require
   exact set equality with the independently constructed PGL family.
2. At $n=7$ and the same primes, exhaust all size-seven matchings and test the
   claimed empty circuit channel.
3. At $n=8$, $p\in\{11,13\}$, exhaust all $8!$ full matchings, record exact
   size-eight circuits, construct the reduced rational-degree-two graph
   subfamily independently if the uniqueness proof closes, and record any
   residual rather than fitting it.
4. At $n=9$, $p\in\{11,13\}$, exhaust all $9!$ full matchings using the exact
   determinant/deletion predicate and record the curve-circuit population.

For every full-grid configuration/size, record the universe, circuit count,
profiles, and signature counts. Across common primes/sizes, compare exact
support sets, not totals only.

## 6. Plants and falsification

1. Known size-four A- and B-fiber ACCEPTs.
2. A five-set containing a four-fiber REJECT and a generic five-set REJECT.
3. Known crossing and PGL/Möbius size-six ACCEPTs.
4. Insert one cell around a known crossing: the larger set contains a proper
   circuit and must REJECT.
5. Wrong curve predicate: seven/eight points on a Möbius graph are dependent
   but contain a size-six circuit; dependency alone must REJECT.
6. Ambient size-ten plant: every ten-set is dependent, but separately test all
   nine-deletions; include one nonminimal REJECT and an ACCEPT if the census
   population is nonempty.
7. Corrupt one actual Kronecker column while retaining a label-valid known
   support; actual rank/minors must reject.
8. True nine-coordinate Kronecker versus six-coordinate block-concatenation
   guard.

Every positive promoted family needs a known-true ACCEPT and a planted REJECT
in-run. A mismatch writes and preserves `broken_results.json`, is disclosed,
and reruns every affected family after a justified fix.

## 7. Runtime, promotion, and lifecycle

Set `sys.dont_write_bytecode=True` before other imports, execute with
`PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, and
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
Hard budget: CPU 7,200 seconds (two hours), wall 8,000 seconds.

The whole gate is **FROZEN-CERTIFIED** only if every requested layer 4--10 is
analytically classified by template/predicate with all claimed constructions
proved and the mandatory/dedicated set equalities passing. Otherwise close
**FROZEN-INCONCLUSIVE**, while explicitly marking each proved layer/subfamily
as certified and naming the first unresolved layer/template plus the strongest
exhaustive range. Freeze before exactly one close. After close update only
`SCOPE_NOTE_H_MIX_CROSSING.md` and `state.json`; Main owns root ledgers.
