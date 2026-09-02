# PREREG-RSPE3D-14-33ROW-SIZE8-RESIDUAL — gate **H-33ROW-SIZE8-RESIDUAL**

Written 2026-09-02 UTC by agent `RsPe3dH2`, after H-33ROW-OPEN was frozen,
closed FROZEN-INCONCLUSIVE, independently verified by Main, and path-scoped
committed. This is the final sequential campaign on this target for this
session. Binding order: preregistration -> path-scoped commit ->
`campaign.py init` -> byte-identical in-run copy recording source commit and
SHA-256 -> controls/census -> analytic record -> `freeze` -> exactly one
`close --verdict`. Write only under `cs/rs-pe3d/`; never alter a frozen run or
root ledger.

No new parameter-dependent computation preceded this registration. The prior
frozen run supplied the two anchors to explain: among all `8!` matchings on
`X=Y={1,...,8}`, there are 560 size-eight circuits over GF(11) and 416 over
GF(13), while its one oriented reduced-rational-degree-two fitter found zero.
Those are inherited evidence, not authority; this run must reproduce them.

## 1. Setting and exact fallback

Use the tied three-row Vandermonde product

\[
h(x,y)=(x^iy^j)_{0\le i,j<3}\in F^9,
\qquad V=H^0(\mathbb P^1\times\mathbb P^1,\mathcal O(2,2)).
\]

For an all-distinct eight-set `S`, let `E_S` be its `9 x 8` evaluation matrix
and

\[
W_S=I_S(2,2)=\ker(E_S^{\mathsf T})\subset V.
\]

The prior gate proved that the only smaller all-distinct circuit is a
six-point Möbius graph and that the entire size-seven circuit layer is empty.
Consequently

\[
S\text{ circuit}\iff \operatorname{rank}E_S=7
\text{ and no six-subset of }S\text{ lies on a }(1,1)\text{ graph}.
\tag{C8}
\]

Rank seven means `dim W_S=2`: projectively, a pencil of `(2,2)` curves. Rank
seven alone is the preregistered wrong predicate because a six-point Möbius
circuit plus two generic points can also have rank seven.

## 2. Candidate structural theorem

Preregister the following candidate and prove it rather than treating it as
an assumption.

**H-CI8.** Every all-distinct size-eight circuit belongs to exactly one of:

1. **degree-two graph:** all eight points lie on a reduced irreducible divisor
   of class `(1,2)` or `(2,1)`, equivalently `x=R(y)` or `y=R(x)` for a
   reduced rational map of exact degree two; or
2. **complete intersection:** `W_S` has no fixed/common curve component, and
   any basis `F,G` of `W_S` cuts out `S` as the full reduced complete
   intersection of two `(2,2)` curves.

Conversely every such graph support and every reduced all-distinct complete
intersection is a circuit. The two classes are disjoint. Thus, for finite
point sets,

\[
N_8=N_{(1,2)\text{-graph}}+N_{(2,1)\text{-graph}}+N_{\mathrm{CI}}.
\tag{D8}
\]

For the registered `n=8` point sets, transpose symmetry plus the predecessor's
zero oriented fit predicts both graph terms are zero, so H-CI8 predicts that
all 560/416 circuits are complete intersections.

### Proof route to discharge

Let `H` be the maximal fixed divisor of the pencil, of class `(a,b)`, and
write `W_S=H W'`, where `W'` has no fixed component. Two members of class
`(c,d)` without a common component intersect in length `2cd` on
`P^1 x P^1`. Use all-distinctness and minimality to exhaust the possibilities:

- `(a,0)` or `(0,b)` contains at most `a` or `b` matching points, while the
  residual intersection has length `2(2-a)(2-b)`, too small for eight;
- `(a,b)=(1,1)` leaves residual class `(1,1)`, hence at most two points off
  `H`; six points on a nondegenerate `(1,1)` component form a forbidden
  Möbius subcircuit, while a reducible `(1,1)` divisor contains at most two
  matching points;
- `(1,2)` and `(2,1)` leave a basepoint-free two-dimensional complementary
  linear system, so all eight points lie on `H`; reducibility cannot hold
  without either a six-point `(1,1)` subcircuit or too few matching points,
  leaving exactly the two reduced degree-two graph classes;
- `(2,2)` cannot support a two-dimensional quotient pencil.

If there is no fixed component, Bezout gives

\[
(2,2)\cdot(2,2)=8,
\]

so the eight distinct points exhaust the complete intersection. Conversely,
the Koszul resolution gives `dim I_S(2,2)=2`; six points on a `(1,1)` curve
would force that curve to divide both generators, contradicting no common
component, so (C8) gives circuitness. Unlike the plane-cubic Cayley-Bacharach
mnemonic, there is no ninth base point here: the intersection number is eight.
The support determines its pencil uniquely as `W_S`, and a fixed-component-
free pencil determines its reduced base locus uniquely. Do not replace this
proved pencil bijection by a count modulo automorphisms unless stabilizers and
injectivity are also proved.

## 3. Exact common-component and resultant certificates

For a basis `F,G` represented as bidegree-at-most-two affine coefficient
arrays, detect a common projective component in two independent exact ways:

1. regard them as homogeneous binary quadratics in `x` over `F[y]`; the
   degree-two Sylvester resultant is identically zero exactly when they share
   a component of positive x-degree (including the x-infinity component);
2. compute the homogeneous gcd/content of all x-coefficient forms in y,
   including a separate y-infinity multiplicity check, to detect y-only
   components.

A pencil is fixed-component-free iff neither test fires. Unit-test the
detector on planted common `(1,1)`, common `(1,2)`, x-fiber, y-fiber, and
coprime pairs.

For every fixed-component-free circuit, normalize the two nonzero elimination
resultants and require

\[
\operatorname{Res}_x(F,G)\doteq\prod_{y\in\pi_Y(S)}(t-y),\qquad
\operatorname{Res}_y(F,G)\doteq\prod_{x\in\pi_X(S)}(t-x).
\tag{R8}
\]

This is a basis-independent certificate up to nonzero scalar that all eight
finite, distinct projected points exhaust the base scheme.

## 4. Mandatory exact census and fingerprints

Use stdlib-only exact arithmetic modulo `p` and true nine-coordinate row-major
Kronecker columns.

For each `p in {11,13}` and `X=Y={1,...,8}`:

1. exhaust all `8!=40320` matchings;
2. compute ranks with two independent implementations;
3. materialize exact sets for rank-seven pencils, direct circuits (rank seven
   plus all seven-deletions rank seven), six-Möbius-subset nonminimal sets,
   each oriented exact-degree-two rational fit, fixed-component-free complete
   intersections, and common-component pencils;
4. require exact set equalities

\[
\text{direct}=\text{graph}_{12}\ \dot\cup\ \text{graph}_{21}
\ \dot\cup\ \text{CI},
\]

\[
\text{CI}=\{S:\operatorname{rank}E_S=7,
W_S\text{ has no fixed component}\},
\]

and require every rank-seven noncircuit to have a six-point Möbius subcircuit
and a common component not belonging to either degree-two graph class;
5. pin direct counts 560/416, pin both oriented graph counts to zero, and hence
   pin CI counts 560/416;
6. verify (R8) on every CI support.

For every direct/residual circuit record: profile `(8,8)`; perfect-matching
bipartite template; both degree sequences; rank drop `8-rank(E_S)`; left
kernel dimension; all eight seven-deletion ranks; the zero/nonzero mask of the
nine `8 x 8` row minors; common-component classification; normalized
resultant degrees; and a label-invariant pencil fingerprint. Record full
histograms, not samples. Record rank-seven wrong-predicate overacceptance and
all common-component/nonminimal populations without fitting them.

## 5. Plants and falsification

Run all plants in-run, exactly:

1. **CI ACCEPT:** deterministically search GF(101) for the lexicographically
   first `c_1<c_2`, `s_1<s_2` such that
   `F=(xy-c_1)(xy-c_2)` and
   `G=(y-x-s_1)(y-x-s_2)` have exactly eight finite, reduced intersection
   points with distinct x- and y-coordinates. Assert no common component,
   rank seven, every deletion rank seven, and both resultants (R8).
2. **Both graph ACCEPTs:** use `x=y^2` and its transpose on parameters 1--8
   over GF(101); require exact-degree-two fit in the intended orientation,
   common-component pencil, and circuitness.
3. **Wrong-degree REJECT:** `x=y^3` on parameters 1--8 over GF(101) is an
   exact degree-three graph but has rank eight and must fail every size-eight
   circuit/degree-two predicate.
4. **Rank-seven-alone / inserted-lower-circuit REJECT:** six identity Möbius
   points plus two generic all-distinct points; require rank seven and all
   `8 x 8` row minors zero, but a six-subcircuit, common `(1,1)` component,
   and failure of circuit minimality.
5. **Ambient plant:** a ten-set is dependent; separately include a known
   circuit ACCEPT and a nonminimal REJECT using exact nine-deletion tests.
6. **Common-component detector controls:** planted common `(1,1)`, common
   `(1,2)`, x-only, y-only, and fixed-component-free pairs must classify
   correctly.
7. **Corruption and layout:** corrupt one actual CI column and require the
   circuit predicate to reject; assert true nine-coordinate Kronecker output
   differs from six-coordinate block concatenation.

Every promoted positive class has a known-true ACCEPT and a planted REJECT. A
mismatch writes and preserves `broken_results.json`, is disclosed, and reruns
every affected family after a justified fix.

## 6. Runtime, promotion, and lifecycle

Set `sys.dont_write_bytecode=True` before other imports, execute with
`PYTHONDONTWRITEBYTECODE=1`, `nice -n 15`, and
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`.
Hard budget: CPU 5,400 seconds (90 minutes), wall 6,000 seconds.

Close **FROZEN-CERTIFIED** only if H-CI8 is proved, the common-component and
resultant certificates pass every circuit at both primes, all exact set
identities and plants pass, and no unexplained residual remains. Otherwise
close **FROZEN-INCONCLUSIVE**, preserving the full fingerprint census and
naming the first exact gap. Freeze before exactly one close. After close,
update only `SCOPE_NOTE_H_MIX_CROSSING.md` and `state.json`; Main owns root
ledgers. The target then rests for this session.
