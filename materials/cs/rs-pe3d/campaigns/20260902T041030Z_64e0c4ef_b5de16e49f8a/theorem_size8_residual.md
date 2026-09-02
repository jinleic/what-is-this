# The all-distinct size-eight residual is a reduced `(2,2)` complete intersection
## gate H-33ROW-SIZE8-RESIDUAL, run `20260902T041030Z_64e0c4ef_b5de16e49f8a`

Agent `RsPe3dH2`, 2026-09-02 UTC. Preregistration:
`prereg/H_33ROW_SIZE8_RESIDUAL_PREREG_2026-09-02.md`, path-scoped source
commit `96de32f1d7e90326f4dad02747981f37cca4416b`, SHA-256
`cf5d4f72778ffec1e9cf678918c3e7aa4bf5b8b30b36188d98de7c984d303a11`.
The in-run `pre_statement.md` is byte-identical to that source by `cmp` and
SHA-256. The preregistration was not amended.

Let

\[
V=H^0(\mathbb P^1\times\mathbb P^1,\mathcal O(2,2)),\qquad \dim V=9,
\]

with affine evaluation columns

\[
h(x,y)=(x^iy^j)_{0\le i,j<3}.
\]

For an all-distinct eight-set `S`, put

\[
W_S=I_S(2,2)=\ker(E_S^{\mathsf T})\subset V.
\]

The preceding frozen gate proved the complete smaller layers: an all-distinct
proper subset of `S` can be a circuit only when six points lie on a unique
nondegenerate `(1,1)`/Möbius graph; there are no size-five or size-seven
circuits. This fact supplies the minimality clause below.

## 1. Classification theorem

**Theorem CI8.** An all-distinct eight-set is a circuit if and only if exactly
one of the following holds.

1. **Degree-two graph.** The eight points lie on a reduced irreducible divisor
   of class `(1,2)` or `(2,1)`, respectively a graph
   `x=R(y)` or `y=R(x)` of a reduced rational map of exact degree two.
2. **Reduced complete intersection.** `W_S` is two-dimensional and has no
   fixed curve component. Any basis `F,G` of `W_S` has no common component and
   its scheme-theoretic intersection is exactly the eight distinct points of
   `S`.

The two cases are disjoint. In the second case `S` determines the pencil
`P(W_S)` uniquely, and that fixed-component-free pencil determines `S`
uniquely as its base locus.

Thus the previous gate's unexplained `560` and `416` circuits are not missing
rational graphs. They are precisely reduced complete intersections of two
bidegree-`(2,2)` curves supported by the chosen finite grids.

### 1.1 From a circuit to a pencil

Circuitness gives

\[
\operatorname{rank}E_S=7,
\qquad \dim W_S=9-7=2.
\tag{1}
\]

The projectivization of `W_S` is therefore a pencil of `(2,2)` curves. Rank
seven alone is not enough: minimality is equivalent here to the absence of a
six-point Möbius subcircuit.

Choose a basis `F,G` of the pencil. Let `H` be their maximal fixed divisor, of
class `(a,b)`, and write

\[
F=HF',\qquad G=HG'.
\]

The residual forms `F',G'` have class `(2-a,2-b)` and no common component.
Two such divisors have intersection length

\[
(2-a,2-b)^2=2(2-a)(2-b).
\tag{2}
\]

Every selected point outside `H` lies in that residual intersection.

### 1.2 Exhausting every possible fixed divisor

There are only the following proper nonzero classes.

- If `H` has class `(a,0)` or `(0,b)`, it contains at most `a` or `b`
  all-distinct matching points. Equation (2) bounds the total by five for
  `(1,0)/(0,1)` and by two for `(2,0)/(0,2)`, never eight.
- If `H` has class `(1,1)`, at most two points lie off `H`. A reducible
  `(1,1)` divisor is one fiber from each ruling and contains at most two
  matching points, impossible. A nondegenerate `(1,1)` divisor therefore
  contains at least six selected points, which are a forbidden Möbius
  subcircuit. This is exactly the rank-seven-but-nonminimal channel seen in
  the census.
- If `H` has class `(1,2)`, the residual system has class `(1,0)`. Its complete
  two-dimensional linear system is basepoint-free, so all eight points lie on
  `H`. If `H` were reducible, its components would be `(1,1)+(0,1)` or a union
  of ruling fibers. Avoiding six points on the `(1,1)` component bounds the
  first possibility by six matching points, and a union of fibers contains at
  most three. Hence `H` is reduced and irreducible. Its equation is linear in
  `x`, so it is the graph `x=R(y)` of an exact degree-two rational map.
- Class `(2,1)` is symmetric and gives `y=R(x)`.
- Class `(2,2)` leaves only the one-dimensional constant residual space and
  cannot support the two-dimensional pencil `W_S`.

This proves that a circuit with a fixed component is exactly one of the two
degree-two graph classes.

### 1.3 The fixed-component-free case

If `F,G` have no common component, Bezout on the quadric surface gives

\[
(2,2)\cdot(2,2)=2\cdot2+2\cdot2=8.
\tag{3}
\]

The eight distinct selected points already lie in this intersection, so they
exhaust it with multiplicity one. Thus `S` is the full reduced complete
intersection.

Conversely, for a reduced complete intersection of two `(2,2)` divisors, the
Koszul resolution twisted by `(2,2)` is

\[
0\longrightarrow\mathcal O(-2,-2)
\longrightarrow\mathcal O^{\oplus2}
\longrightarrow\mathcal I_S(2,2)\longrightarrow0.
\]

The relevant `H^0` and `H^1` of `O(-2,-2)` vanish, so
`dim I_S(2,2)=2` and `rank E_S=7`. If six selected points lay on a
nondegenerate `(1,1)` curve `L`, the restriction of each generator to
`L isomorphic to P^1` would be a degree-four section with six distinct zeros;
each restriction would vanish identically, forcing `L` to divide both
generators. This contradicts the no-common-component hypothesis. There is no
smaller all-distinct circuit, so all seven-deletions are independent and `S`
is a circuit.

For a reduced `(1,2)` or `(2,1)` graph, `O(2,2)` restricts to `O(6)` and the
restriction kernel is the graph equation times a two-dimensional complementary
linear system. Seven distinct parameters evaluate independently and eight
have rank seven, so those graph supports are also circuits. An irreducible
`(1,2)` and an irreducible `(2,1)` curve intersect in only
`(1,2)\cdot(2,1)=5` points, so an eight-set cannot belong to both cases.
This completes both directions and disjointness.

## 2. Exact count principle; no unjustified orbit quotient

For finite `X,Y`, let `D_R(X,Y)` be the finite parameters in `Y` whose
images under `R` land in `X`, and put

\[
\mathcal A_R(X,Y)=
\{U\subseteq D_R(X,Y): |U|=8,\ |R(U)|=8\}.
\]

The image-cardinality clause is essential because a degree-two map need not be
injective on all of `D_R`. Five paired points determine a
degree-at-most-two rational map uniquely, so the two all-distinct graph counts
are the injective incidence sums

\[
N_{12}=\sum_{\deg R=2}|\mathcal A_R(X,Y)|,\qquad
N_{21}=\sum_{\deg R=2}|\mathcal A_R(Y,X)|.
\tag{4}
\]

Let `P_CI(X,Y)` be the set of two-dimensional subspaces
`W subset V` whose base scheme is reduced, fixed-component-free, consists of
eight points in `X times Y`, and has distinct projections. The pencil/base-
locus uniqueness proved above is a bijection, not a many-to-one orbit map.
Hence

\[
\boxed{N_8(X,Y)=N_{12}+N_{21}+|P_{CI}(X,Y)|.}
\tag{5}
\]

No quotient by `PGL(2) x PGL(2)`, and no closed automorphism-orbit formula, is
claimed: that would require a separate stabilizer analysis. Formula (5) is the
correct exact structural count principle. Its complete-intersection term can
be field- and point-set-dependent, exactly as `560` versus `416` demonstrates.

## 3. Resultants certify the complete base locus

Write a form as a homogeneous quadratic in `x` with coefficients in `F[y]`.
For a basis `F,G`, the exact `4 x 4` Sylvester resultant in `x` is identically
zero exactly when the primitive forms share a component of positive x-degree.
The homogeneous gcd of all coefficient forms detects a y-only component;
the degree deficit detects the y-infinity factor. Swapping the variables gives
an independent symmetric route. The instrument required these two routes to
agree on every rank-seven pencil (`2192+1344=3536` exact comparisons), and
unit-tested common `(1,1)`, common `(1,2)`, x-only, y-only, and coprime pairs.

For a fixed-component-free basis, both resultants are nonzero of degree at most
eight. Every projected coordinate of `S` is a root. Since the projections have
eight distinct finite values, normalization gives

\[
\operatorname{Res}_x(F,G)\doteq
\prod_{y\in\pi_Y(S)}(t-y),\qquad
\operatorname{Res}_y(F,G)\doteq
\prod_{x\in\pi_X(S)}(t-x).
\tag{6}
\]

All `560+416=976` complete intersections passed both polynomial identities,
coefficient for coefficient. This also proves that the eight affine points
exhaust the base scheme. The plane-cubic slogan that eight points determine a
ninth does **not** transfer here: on `P^1 x P^1`, two `(2,2)` divisors have
intersection number eight, not nine.

## 4. Exact exhaustive census and fingerprints

All `8!=40320` matchings on `X=Y={1,...,8}` were exhausted independently at
each prime. Two separately written finite-field rank algorithms agreed on all
`80640` matchings.

| field | rank-seven pencils | circuits | rank-seven noncircuits | common component | graph `(1,2)` | graph `(2,1)` | complete intersections |
|---|---:|---:|---:|---:|---:|---:|---:|
| GF(11) | 2192 | 560 | 1632 | 1632 | 0 | 0 | 560 |
| GF(13) | 1344 | 416 | 928 | 928 | 0 | 0 | 416 |

The following were exact support-set identities, with sorted-set SHA-256
digests retained in `controls_results.json`:

\[
\text{circuits}=\text{graph}_{12}\ \dot\cup\
\text{graph}_{21}\ \dot\cup\ \text{CI}=\text{CI},
\]

\[
\text{rank-seven noncircuits}=\text{six-Möbius-containing sets}
=\text{common-component non-graph pencils}.
\]

Thus the bare pencil/rank-seven hypothesis overaccepts `1632` supports over
GF(11) and `928` over GF(13). Minimality is load-bearing, not cosmetic.

Every residual circuit at a given prime has one and the same complete
fingerprint:

- profile `(8,8)`;
- bipartite template: eight isolated `K_2` components, i.e. a perfect matching;
- both degree sequences `(1,1,1,1,1,1,1,1)`;
- evaluation rank seven, rank drop one, and left-kernel dimension two;
- all eight seven-deletion ranks equal seven;
- all nine `8 x 8` row minors vanish, zero mask `511`;
- no common component and neither rational-degree-two graph orientation;
- both elimination resultants have degree eight and equal the normalized node
  polynomials in (6).

That fingerprint occurs exactly `560` times over GF(11) and `416` over GF(13),
with no second residual signature.

## 5. Plants, defect disclosure, and runtime

The deterministic GF(101) complete-intersection ACCEPT is

\[
F=(xy-1)(xy-2),\qquad G=(y-x-1)(y-x-3),
\]

with eight finite base points

`(1,2),(16,19),(22,23),(27,30),(71,74),(78,79),(82,85),(99,100)`.

It has distinct projections, no common component, rank seven, all deletion
ranks seven, all nine maximal minors zero, and both exact node resultants.
The two oriented graph ACCEPTs `x=y^2` and `y=x^2` passed exact degree-two
fits, common-component detection, and circuitness. The wrong-degree
`x=y^3` plant had rank eight and rejected. Six identity-Möbius points plus
`(7,8),(8,7)` had rank seven and all nine maximal minors zero, but exposed its
six-subcircuit and common `(1,1)` component and rejected. Ambient-ten ACCEPT
and nonminimal REJECT, all common-component detector controls, a corrupted CI
column (rank rose to eight), and the nine-coordinate Kronecker versus
six-coordinate concatenation guard all passed.

The first invocation failed at Python parse time because one surgical edit
left escaped quote characters in a dictionary literal. No module code or
parameter-dependent computation executed. The defective source and a manual
broken record are preserved as `run_size8_residual_syntax_error.py` and
`broken_results_syntax_error.json`. After removing every such escape, the
source was parsed without execution and the entire preregistered family was
run from the beginning. `defect_log.json` records the defect; no prereg
amendment occurred.

All `72` final assertions passed. Runtime was `18.887597` CPU seconds and
`19.621288208989426` wall seconds against hard caps of `5400/6000`, under
`nice -n 15`, with all thread caps one and both bytecode guards active.

## 6. Verdict

H-CI8 is proved in both directions; the fixed-divisor cases are exhausted;
the pencil and reduced base locus determine one another; exact common-
component, minimality, graph-orientation, and resultant support sets agree at
both primes; every one of the `560/416` former residuals is explained by the
single complete-intersection fingerprint; and every plant passes after the
disclosed full rerun. There is no unexplained residual in the registered
all-distinct size-eight matching channel; repeated-index size-eight templates
were outside this gate and are not silently claimed closed.

Under the preregistered rule the gate verdict is **FROZEN-CERTIFIED**. This is
the final campaign on `cs/rs-pe3d` for this session; after closure the target
rests with larger sizes governed by the already frozen exact deletion
predicates rather than an asserted unproved orbit count.
