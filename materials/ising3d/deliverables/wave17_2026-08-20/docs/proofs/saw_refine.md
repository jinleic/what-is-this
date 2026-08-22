# First-backtrack SAW refinement and a chained cubic-lattice certificate

Artifacts: `experiments/e92_saw_refine.py`,
`tests/test_saw_refine.py`, and `results/bounds/saw_refine.json`.

## Scope and notation

Let \(c_r\) denote the number of rooted nearest-neighbour self-avoiding
walks (SAWs) of length \(r\) in \(\mathbb Z^3\), and let
\(\mu=\lim_{r\to\infty}c_r^{1/r}\).  Let \(a_r\) count rooted \(r\)-step
SAWs avoiding the fixed vertex \(e_1\).  Cubic-lattice symmetry makes the
choice of that one neighbour immaterial.

For \(r\geq0\), let \(b_r\) count rooted \(r\)-step SAWs avoiding the
corner pair

\[
\{e_1,e_1+e_2\}.\tag{1}
\]

Thus one forbidden vertex is a neighbour of the origin and the other is the
adjacent square corner.  This is the exact two-vertex geometry forced below;
it is not an assumed asymptotic approximation.

For \(m\geq2\), write \(T_m\) for the number of \(m\)-step SAWs whose last
two directed steps are perpendicular.  For \(m\geq3\), write \(D_m\) for
the number whose last three chronological directed steps are
\(s,t,-s\), where \(s\perp t\).  Finally, let
\(A_m^{\rm turn}\) count the \(m\)-step SAWs avoiding \(e_1\) that end in a
turn.

The first-backtrack injection from `proofs/kc_interval2.md`, §3, is used as
an input and is reproved as the first family below.

## A three-family invalid-pair injection

**[THEOREM] Multi-neighbour first-backtrack inequality.**  For every
\(m,n\geq2\),

\[
 c_{m+n}\leq c_m\bigl(c_n-a_{n-1}\bigr)
             -T_m b_{n-2}-D_m\bigl(a_{n-1}-b_{n-2}\bigr).\tag{2}
\]

*Proof.*  Pair an \(m\)-step SAW \(\omega\) with an \(n\)-step SAW
\(\eta\), translating \(\eta\) to start at \(\omega_m\) when
concatenating.  Valid pairs are in bijection with \((m+n)\)-step SAWs.
It is therefore enough to construct disjoint injected families of invalid
pairs.

Fix \(\omega\), translate its endpoint to \(0\), and put

\[
 \sigma=\omega_m-\omega_{m-1},\qquad
 \tau=\omega_{m-1}-\omega_{m-2}.
\]

The following families have fixed initial segments, so their tails recover
their preimages and each construction is injective.

1. **F1 (the original first backtrack).**  Force
   \(\eta_1=-\sigma=\omega_{m-1}\).  The remaining \((n-1)\)-step tail
   starts at \(\omega_{m-1}\), avoids \(0\), and has exactly
   \(a_{n-1}\) choices after a translation and a cubic symmetry.  The
   concatenation revisits \(\omega_{m-1}\).  This gives
   \(c_m a_{n-1}\) invalid pairs.

2. **F2 (the square-corner family).**  Suppose \(\sigma\perp\tau\), so
   \(\omega\) is turn-ending.  Force
   \[
     \eta_1=-\tau,\qquad \eta_2=-\tau-\sigma=\omega_{m-2}.
   \]
   The tail after \(\eta_2\) must avoid the already used vertices
   \(0\) and \(-\tau\).  Translating by \(\sigma+\tau\), these become
   \(\sigma+\tau\) and \(\sigma\), an oriented cubic-lattice copy of
   (1).  Hence it has \(b_{n-2}\) choices.  The concatenation revisits
   \(\omega_{m-2}\).  F1 and F2 are disjoint because their first steps are
   \(-\sigma\) and \(-\tau\), respectively.  Thus F2 contributes
   \(T_m b_{n-2}\) pairs disjoint from F1.

3. **F3 (the fold-return family).**  Suppose the final three steps of
   \(\omega\) are \(s,t,-s\), with \(s\perp t\).  In the above notation,
   \(\sigma=-s\), \(\tau=t\), and \(\omega_{m-3}=-\tau\).  Force
   \(\eta_1=\omega_{m-3}\), then choose an \((n-1)\)-step self-avoiding
   tail from that vertex avoiding \(0\).  There are \(a_{n-1}\) choices
   and every result revisits \(\omega_{m-3}\).  Its overlap with F2 is
   exactly the subfamily whose next vertex is
   \(\omega_{m-2}=-\sigma-\tau\): after that forced next step, the
   remaining tail avoids the two vertices \(0\) and \(-\tau\), so the
   overlap has \(b_{n-2}\) choices by the same corner geometry.  Remove
   this overlap from F3.  The residual F3 family has
   \(D_m(a_{n-1}-b_{n-2})\) members.  It is disjoint from F1 because the
   F1 first step is \(-\sigma\ne-\tau\), and it is disjoint from F2 by
   the overlap removal.

The union of these three families is contained in the invalid pairs and has
cardinality
\[
 c_m a_{n-1}+T_m b_{n-2}+D_m(a_{n-1}-b_{n-2}).
\]
Subtracting it from the \(c_m c_n\) total pairs proves (2). \(\square\)

**[LEMMA] Turn/corner identity.**  For every \(m\geq2\),

\[
 T_m=24b_{m-2}.\tag{3}
\]

*Proof.*  Reverse a turn-ending walk.  It becomes a walk whose first two
oriented steps are perpendicular.  There are \(6\cdot4=24\) ordered choices
for this prefix.  After translating the vertex after that prefix to the
origin and applying a signed coordinate permutation, the remaining tail
must avoid exactly the corner pair (1).  Conversely, a choice of oriented
prefix and a \(b_{m-2}\)-tail reconstructs one reversed turn-ending SAW.
\(\square\)

Consequently the F2 part of (2) may equivalently be written as the symmetric
correction \(24b_{m-2}b_{n-2}\).  The fold correction remains useful because
it adds a further disjoint invalid subfamily.

**[COMPUTATION]** Direct exact visited-set enumeration gives

\[
\begin{aligned}
(b_0,\ldots,b_{10})={}&(1,5,24,117,559,2690,12817,61263,291093,\\
&1385747,6570755),\\
(c_{10},a_{10})={}&(8809878,6989025).
\end{aligned}
\]

The producer checks (2) for every \(m+n\leq10\).  The standalone verifier
re-enumerates all series through depth \(10\), exhaustively enumerates every
pair of components for \(m+n\leq9\), and confirms both the valid-pair
bijection and coverage by F1--F3.  At \((m,n)=(2,2)\) and \((3,2)\), the
three injected families exhaust the invalid pairs exactly.  These finite
checks corroborate the injection; they are not the proof of its all-size
statement.

## Chaining the all-\(m\) inequality

**[THEOREM] Chained first-backtrack consequence.**  For every \(n\geq1\),

\[
 \mu^n\leq c_n-a_{n-1}.\tag{4}
\]

*Proof.*  Set \(X_n=c_n-a_{n-1}\).  The F1 portion of the preceding proof
(the first-backtrack theorem) gives
\(c_{m+n}\leq c_mX_n\) for every \(m\geq1\).  Induction with
\(m=(k-1)n\) yields

\[
 c_{kn}\leq c_nX_n^{k-1}\qquad(k\geq1).\tag{5}
\]

Ordinary SAW concatenation gives submultiplicativity of \(c_r\), so Fekete's
lemma gives \(c_r^{1/r}\to\mu\), including along the subsequence \(r=kn\).
Taking \(kn\)-th roots in (5) and letting \(k\to\infty\) proves (4).
\(\square\)

The corner and fold terms in (2) are a genuine strict finite refinement, but
they are not automatically chainable: a repeated use would require a
uniform lower bound on the turn/fold population of the growing prefix.  No
such unproved density assertion is used here.

## Two-sided triple refinement

**[THEOREM]**  For \(\ell\geq1\), \(m\geq3\), and \(n\geq2\),

\[
\begin{aligned}
 c_{\ell+m+n}\leq c_\ell\big[&
  (c_m-a_{m-1})(c_n-a_{n-1})\\
  &-(T_m-A_{m-1}^{\rm turn})b_{n-2}\big].\tag{6}
\end{aligned}
\]

*Proof.*  In a triple \((\omega,\eta,\zeta)\), take the F1 invalid family at
the first junction (A), the F1 family at the second junction (B), and the
F2 corner family at the second junction (C).  Their sizes are

\[
 c_\ell a_{m-1}c_n,\qquad c_\ell c_m a_{n-1},\qquad c_\ell T_m b_{n-2}.
\]

The families B and C are disjoint because their forced first \(\zeta\)-steps
are the distinct reversals of the final and penultimate \(\eta\)-steps.
The intersections A\(\cap\)B and A\(\cap\)C have sizes
\(c_\ell a_{m-1}a_{n-1}\) and
\(c_\ell A_{m-1}^{\rm turn}b_{n-2}\), respectively: under A, remove the
forced first step of \(\eta\), and its remaining \((m-1)\)-tail is an
\(a\)-walk, ending in a turn in the latter intersection.  Inclusion--exclusion
for A, B, C, followed by subtraction from all \(c_\ell c_m c_n\) triples,
gives (6). \(\square\)

**[COMPUTATION]** Equation (6) is checked for every
\(\ell+m+n\leq10\) in its stated range \(m\geq3,n\geq2\); the iterated F1
triple inequality is checked for every positive split in that total range.
For the enumerated values, \(T_m>A_{m-1}^{\rm turn}\) whenever the correction
is tested, so (6) is strictly smaller than the pairwise F1 right side in each
such finite case.  No all-\(m\) strictness claim is made.

## Explicit lower families for \(a_r\)

Let \(\mathcal M=\{-e_1,e_2,e_3\}\).  The monotone family has \(3^r\)
walks and avoids \(e_1\).  Define
\[
 s(0)=1,\qquad s(j)=2\,3^{j-1}\quad(j\geq1),
\]
and
\[
 G_r=\sum_{k=0}^{r-2}3^k s(r-2-k),\qquad
 D_r=\sum_{k=0}^{r-3}3^k s(r-3-k),\qquad
 L_r=3^r+4G_r+2D_r.\tag{7}
\]

**[THEOREM]** \(a_r\geq L_r\) for every \(r\geq0\).

*Proof.*  Besides \(\mathcal M\), use four one-defect families: after a
monotone prefix, take \(-e_2\) after either \(e_3\) or \(-e_1\), or take
\(-e_3\) after either \(e_2\) or \(-e_1\); use a monotone suffix whose first
step is restricted away from the single direction that would revisit the
pre-defect path.  Each position of the defect contributes
\(3^k s(r-2-k)\).  Two analogous three-step fold defects,
\((-e_1,e_2,e_1)\) and \((-e_1,e_3,e_1)\), contribute
\(3^k s(r-3-k)\).  Coordinate monotonicity before and after the defect
proves self-avoidance; all paths have nonpositive first coordinate and hence
avoid \(e_1\).  The families are disjoint because they have either no defect
or exactly one nonmonotone step with different type and prescribed preceding
local pattern.  Summing gives (7). \(\square\)

**[COMPUTATION]** The exact formula gives
\[
 L_{29}=741377533262625,\qquad L_{35}=644233352324156721,
\]
and exhaustive classification independently agrees with all seven family
counts through \(r=10\).

A stronger finite family is useful at \(r=35\).  Put
\(h(x,y,z)=-x+y+z\).  Let \(\mathcal B\) be the set of 11-step SAWs from
height \(0\) to height \(7\), with every proper non-origin vertex having
height strictly between \(0\) and \(7\).  Let \(\mathcal H\) be the set of
13-step SAWs from the origin for which every non-origin vertex has positive
height.

**[COMPUTATION]** Independent direct enumerations give
\[
 |\mathcal B|=729000,\qquad |\mathcal H|=142016661.\tag{8}
\]

**[THEOREM] Finite height-slab certificate.**

\[
 a_{35}\geq |\mathcal B|^2|\mathcal H|
          =75473476338501000000.\tag{9}
\]

*Proof.*  Concatenate a block in \(\mathcal B\), a translated second block
in \(\mathcal B\), and a translated tail in \(\mathcal H\).  The first
block has interior heights \(1,\ldots,6\) and endpoint height \(7\); the
second has interior heights \(8,\ldots,13\) and endpoint height \(14\);
the tail's post-junction vertices have height greater than \(14\).  Hence
pieces have disjoint interiors and meet only at their prescribed splice
endpoints.  The resulting path is self-avoiding and uniquely recovers its
three fixed-length components.  All non-origin vertices have positive
height, whereas \(h(e_1)=-1\), so it avoids \(e_1\).  This gives the
product injection into the walks counted by \(a_{35}\). \(\square\)

## External count and exact endpoint comparison

**[EXTERNAL]** Schram--Barkema--Bisseling, Table I, supplies
\[
 c_{36}=2941370856334701726560670.
\]
The cached source has SHA-256
`898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12`;
both producer and clean-room verifier hash it before using the count.

Combining (4) at \(n=36\) with (9) gives the following finite certificate.

**[THEOREM]** With
\[
 M=2941295382858363225560670,
\]
one has
\[
 \mu^{36}\leq M.\tag{10}
\]
Moreover, the exact integer comparison
\[
 M^{36}<c_{36}^{36}\tag{11}
\]
proves that (10) strictly improves the incumbent \(c_{36}^{1/36}\) bound.
The negative difference \(M^{36}-c_{36}^{36}\) has bit length \(2917\) and
SHA-256 magnitude digest
`3700e3c69abd918126e9828947ab28a6de1fe69483da77e3b3b178c163ed3e02`.

Using the standard high-temperature SAW domination recorded in the prior
SAW-bound proof, (10) implies the strictly stronger symbolic Ising bound
\[
 K_c\geq\operatorname{atanh}(M^{-1/36})
      >\operatorname{atanh}(c_{36}^{-1/36}).\tag{12}
\]
The strict comparison in (12) rests on (11), not on a decimal calculation.

**[COMPUTATION]** Directed interval arithmetic at 90 decimal digits formats
the lower endpoint in (12), rounded down to 40 places, as
`0.2122120570061123902353626109344172245745`.  The previous displayed floor
was `0.2122119011661678393310862783954278184914`.

## Bounded negative results and the missing ratio lemma

**[FALSIFIED]** For the explicitly enumerated direct-pair screening set
\(m\in\{2,\ldots,14\}\cup\{30,\ldots,36\}\), \(1\leq n\leq8\), no finite
right-hand side from F1 plus every available F2/F3 correction satisfies the
integer screening comparison against \(c_{36}^{1/36}\).  The smallest
shortfall is \(41\), at \((m,n)=(2,1)\), where the candidate right-hand side
is \(150\).  This is an honestly bounded finite screening result, not an
all-pairs theorem and not a substitute for the chained proof above.

**[COMPUTATION]** At the representative three-block scale
\((\ell,m,n)=(30,4,2)\), (6) yields the finite right-hand side
\(3962496266420283709722030\), whereas the chained \(n=36\) certificate
has right-hand side \(M\).  Thus the verified triple correction does not
beat the chained certificate at this scale.

**[UNRESOLVED]** A sufficient (deliberately strong) route from (4) toward
\(\mu^2\leq c_{32}/c_{30}\) would be
\[
 a_{29}\geq c_{30}-\left\lfloor\frac{c_{32}}{c_{30}}\right\rfloor
          =270569905525454674592.\tag{13}
\]
It would make \(\mu^{30}\leq\lfloor c_{32}/c_{30}\rfloor\), hence imply
the desired weaker square inequality.  The explicit defect family supplies
only \(L_{29}=741377533262625\).  The exact gap factor is
\[
\frac{270569905525454674592}{741377533262625},
\]
and no claim closes that gap.

## Reproduction

From the repository root:

```sh
timeout 600 .venv/bin/python experiments/e92_saw_refine.py
timeout 600 .venv/bin/python tests/test_saw_refine.py
```

The latest producer run took approximately 85 seconds and the clean-room
verifier approximately 91 seconds; the displayed timeout exceeds twice each
measured wall time.
