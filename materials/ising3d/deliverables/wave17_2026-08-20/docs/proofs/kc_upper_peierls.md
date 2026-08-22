# Quantitative three-dimensional Peierls bound with an exact area-24 head

## 1. Result and honest comparison

For the nearest-neighbour ferromagnetic Ising model on \(\mathbb Z^3\), write
\(K=\beta J\).  The argument below proves the following valid, but weak, upper
bound.

**[THEOREM]**
\[
  \boxed{K_c\le K_0:=\frac{848449}{500000}=1.696898.}
\]

This does **not** improve the certified incumbent
\[
  K_c\le 0.2527310098586630030260020266135701299926.
\]
The incumbent therefore remains unchanged.  The new mathematical progress is
the exhaustive rooted outer-contour census through area \(24\), including its
independent small-area cross-check, and an explicit negative certificate for
this head-plus-tail strategy.  The exact contour head is combined with a fully
proved exponential tail, and the resulting majorant fails at the preceding
millionth \(1.696897\).  More decisively, this tail converges only when
\(K>\tfrac12\log(11^{11}/10^{10})\), so no enlargement of a *finite* exact head
can make this particular tail proof compete with the incumbent.

**[COMPUTATION]** The exact rooted outer-contour counts obtained here are

\[
\begin{array}{c|rrrrrrrrrr}
 A&6&8&10&12&14&16&18&20&22&24\\ \hline
 N(A)&1&0&6&0&45&12&332&240&2538&3040.
\end{array}
\]

All entries are exact finite integers.  No all-area formula and no asymptotic
claim is inferred from them.

## 2. Finite-volume Peierls inequality

### 2.1 Contours, including cubical self-contacts

Let \(\Lambda=[-L,L]^3\cap\mathbb Z^3\) contain the origin and impose fixed
plus spins outside \(\Lambda\).  At zero field the unnormalised weight is

\[
 w_\Lambda(\sigma)=
 \exp\!\left(K\sum_{\langle x,y\rangle:\{x,y\}\cap\Lambda\ne\varnothing}
                  \sigma_x\sigma_y\right).
\]

Suppose \(\sigma_0=-1\), and let \(C=C_0(\sigma)\) be the nearest-neighbour
connected component of minus spins containing \(0\).  It is finite because the
exterior spins are plus.  Let \(C_\infty\) be the infinite component of
\(\mathbb Z^3\setminus C\), and fill all finite holes:

\[
 P=P(C):=\mathbb Z^3\setminus C_\infty.
\]

Thus \(P\) and its complement are both face-connected, \(P\) is finite, and
\(0\in P\).  Associate a closed unit cube centred at every site of \(P\).  For
each nearest-neighbour edge with one endpoint in \(P\) and one outside, select
the perpendicular dual unit plaquette.  Their set is the **outer contour**
\(\Gamma(C)=\partial^*P\), and \(|\Gamma|\) is its number of plaquettes.

At a dual edge, the four surrounding primal cells occur cyclically.  Membership
in \(P\) changes an even number of times around that cycle.  Therefore every
dual edge is incident to \(0,2\), or \(4\) selected plaquettes.  This is the
closedness condition used here.  Incidence four is allowed: the proof does not
silently discard cubical self-contacts.

**[LEMMA]** \(\partial^*P\) is connected when plaquettes sharing a dual edge are
adjacent.

**Proof.** Regard \(\partial^*P\) as the primal edge cut
\(\delta P\), using the one-to-one correspondence between a primal edge and
its perpendicular dual plaquette.  Because both \(P\) and \(P^c\) are
face-connected, \(\delta P\) is an inclusion-minimal nonempty cut (a bond):
after deleting any proper subset of \(\delta P\), an undeleted crossing edge,
together with paths inside \(P\) and \(P^c\), still connects the whole lattice.

Suppose that the dual plaquettes split into two or more shared-edge components,
and let \(E_1\) be the nonempty proper subset of primal edges dual to one
component.  Every elementary primal square meets \(E_1\) evenly.  Indeed, its
four crossing dual plaquettes meet at one dual edge; all selected ones at that
edge belong to the same shared-edge component, and the full cut \(\delta P\)
meets every cycle evenly.  Every finite lattice cycle is a mod-two sum of
elementary squares: commute consecutive steps in different coordinate
directions, adding the square between the two orders, until opposite steps
cancel.

It follows that \(E_1\) meets every finite cycle evenly.  Fix a base vertex and
label any vertex by the parity of \(E_1\)-edges on a path from the base to that
vertex.  The cycle parity makes this label path-independent, and an edge lies
in \(E_1\) exactly when its endpoint labels differ.  Hence \(E_1\) is itself a
nonempty cut.  This contradicts the inclusion-minimality of the bond
\(\delta P\), because \(E_1\subsetneq\delta P\).  Therefore the dual
plaquettes have one shared-edge component. \(\square\)

Finally, if \(n=|P|\) and \(b(P)\) is the number of nearest-neighbour pairs with
both endpoints in \(P\), face counting gives

\[
  |\partial^*P|=6n-2b(P).                                      \tag{2.1}
\]

In particular every contour area is even.

### 2.2 The flip map and the direction of every inequality

For a possible outer contour \(\Gamma\) enclosing the origin, let \(E_\Gamma\)
be the set of configurations for which \(\sigma_0=-1\) and the outer contour of
\(C_0(\sigma)\) is exactly \(\Gamma\).  Its bounded cell set \(P=P(\Gamma)\) is
determined by \(\Gamma\).  Flip **all** spins in \(P\):

\[
 (F_\Gamma\sigma)_x=
 \begin{cases}
 -\sigma_x,&x\in P,\\
 \sigma_x,&x\notin P.
 \end{cases}                                                    \tag{2.2}
\]

Every boundary edge of \(P\) has its inside endpoint in \(C\), hence spin
\(-1\), and its outside endpoint is adjacent to \(C\), hence spin \(+1\).
(A filled finite hole cannot have a face-neighbour in \(C_\infty\), or the two
complement components would be the same.)  Thus every boundary bond is
unsatisfied before the flip and satisfied afterwards.  Bonds wholly inside
\(P\) have both endpoints flipped and are unchanged; exterior bonds are also
unchanged.  Consequently, for every \(\sigma\in E_\Gamma\),

\[
  w_\Lambda(\sigma)
   =e^{-2K|\Gamma|}\,w_\Lambda(F_\Gamma\sigma).                 \tag{2.3}
\]

For fixed \(\Gamma\), \(F_\Gamma\) is an involution and hence injective on
\(E_\Gamma\).  Its image need not be all configurations.  Dropping that image
restriction is precisely the safe overcount:

\[
\begin{aligned}
 \mu^+_{\Lambda,K}(E_\Gamma)
 &=\frac{e^{-2K|\Gamma|}}{Z^+_{\Lambda,K}}
   \sum_{\sigma\in E_\Gamma}w_\Lambda(F_\Gamma\sigma)\\
 &\le e^{-2K|\Gamma|}.                                        \tag{2.4}
\end{aligned}
\]

The events \(E_\Gamma\) partition \(\{\sigma_0=-1\}\).  Summing (2.4), and then
enlarging the finite-box contour family to all lattice contours enclosing the
origin, only increases a nonnegative sum.  If \(N(A)\) denotes the number of
such rooted outer contours of area \(A\), then uniformly in \(\Lambda\),

\[
 \mu^+_{\Lambda,K}(\sigma_0=-1)
 \le \sum_{\Gamma\ni 0}e^{-2K|\Gamma|}
 =\sum_{A\ge6}N(A)e^{-2KA}.                                   \tag{2.5}
\]

If the last sum is \(<1/2\), then

\[
 \langle\sigma_0\rangle^+_{\Lambda,K}
 =1-2\mu^+_{\Lambda,K}(\sigma_0=-1)>0                         \tag{2.6}
\]

uniformly in the volume.  Taking the monotone plus-boundary thermodynamic limit
preserves the strict positive lower bound.  Therefore the plus state has
positive spontaneous magnetisation at \(K\), and by the definition and
ferromagnetic monotonicity of the ordered regime, \(K_c\le K\).

This proves the complete Peierls inequality chain.  Notice that images for
different contours may overlap; (2.4) is used separately and the subsequent
sum is an upper bound, so such overlap cannot invalidate the argument.

## 3. Exact contour census through area 24

### 3.1 Contours and filled polycubes

The bounded side of an outer contour is a finite face-connected union \(P\) of
unit cells whose complement is face-connected.  Conversely, every such \(P\)
produces its outer contour \(\partial^*P\).  The contour determines its bounded
side, so this correspondence is injective even in the presence of cubical
self-contacts.

The enumeration first quotients translations but not rotations or reflections.
For a translation class with \(n\) cells, exactly \(n\) translations place one
of its cells at the fixed origin.  No two do so in the same way for a finite
set.  Thus a canonical translation class of size \(n\) contributes \(n\), not
one, to \(N(A)\).

### 3.2 Discrete isoperimetry and a rigorous cell cutoff

**[LEMMA: discrete Loomis--Whitney isoperimetry]** Every finite
\(P\subset\mathbb Z^3\) with \(n=|P|\) satisfies

\[
  |\partial^*P|\ge6n^{2/3}.                                    \tag{3.1}
\]

**Proof.** Let \(p_{xy},p_{xz},p_{yz}\) be the cardinalities of the three
coordinate projections of \(P\).  Each nonempty column parallel to the missing
coordinate direction has a first and a last cell, and therefore at least two
boundary faces normal to that direction.  Hence

\[
 |\partial^*P|\ge2(p_{xy}+p_{xz}+p_{yz}).                       \tag{3.2}
\]

For completeness, here is the finite Loomis--Whitney step.  Let
\(a(x,y),b(x,z),c(y,z)\) be the indicator functions of the three projections.
Every point of \(P\) contributes one to \(abc\), so

\[
 n\le\sum_{x,y,z}a(x,y)b(x,z)c(y,z).
\]

For fixed \(x\), Cauchy--Schwarz in \((y,z)\), followed by Cauchy--Schwarz in
\(x\), gives

\[
\begin{aligned}
 \sum_{x,y,z}abc
 &\le \|c\|_2\sum_x
       \sqrt{\sum_y a(x,y)^2}\sqrt{\sum_z b(x,z)^2}\\
 &\le \|a\|_2\|b\|_2\|c\|_2
 =\sqrt{p_{xy}p_{xz}p_{yz}}.
\end{aligned}
\]

Thus \(n^2\le p_{xy}p_{xz}p_{yz}\).  AM--GM in (3.2) now yields (3.1).
\(\square\)

Cubing (3.1) gives \(A^3\ge216n^2\).  At \(A\le24\),

\[
  216n^2\le24^3=216\cdot8^2,
\]

so \(n\le8\).  Conversely \(216\cdot9^2>24^3\).  Therefore enumerating all
polycubes through eight cells is a **complete**, not heuristic, census of every
outer contour with area at most 24.

### 3.3 Canonical BFS and exact counts

Start from the singleton.  From every canonical \(n\)-cell set, adjoin each of
its vacant face-neighbours, translate the three coordinate minima to zero, sort
the coordinate triples, and deduplicate exact tuples.  This generates every
translation class: every finite connected graph with at least two vertices has
a vertex whose deletion leaves it connected (take a leaf of a spanning tree),
so induction supplies a predecessor.  Canonicalisation removes translations
only, leaving every lattice orientation distinct.

The exact numbers of fixed polycubes modulo translation are

\[
\begin{array}{c|rrrrrrrr}
 n&1&2&3&4&5&6&7&8\\ \hline
 \#&1&3&15&86&534&3481&23502&162913.
\end{array}
\]

For every candidate with area at most 24, an exact flood fill in a one-cell
expanded bounding box checks that the complement is connected.  No candidate
in this area range has a cavity.  Boundary faces are counted directly, and each
translation class is weighted by its number \(n\) of origin placements.  This
gives

\[
\begin{array}{c|rrrrrrrrrr}
 A&6&8&10&12&14&16&18&20&22&24\\ \hline
 \text{classes modulo translation}&1&0&3&0&15&3&83&48&504&505\\
 N(A)&1&0&6&0&45&12&332&240&2538&3040.
\end{array}                                                     \tag{3.3}
\]

In particular, \(N(6)=1\) is the unit cube, \(N(8)=0\), and \(N(10)=6\) are
the six rooted orientations of the two-cell domino.

### 3.4 Independent direct plaquette search through area 10

A second engine does not grow cells.  It represents plaquette and edge centres
in doubled integer coordinates.  Starting from every possible positive-\(x\)-ray
root, it maintains the set of dual edges having odd selected-plaquette
incidence.  If this set is nonempty, the lexicographically first odd edge is
chosen and each unselected incident plaquette is tried.  If it is empty before
the area cap, every edge-adjacent plaquette starts a possible connected
extension.  Exact frozensets remove duplicate growth histories.

This search is exhaustive for a connected closed target \(S\).  At a partial
set with an odd edge, even incidence of \(S\) guarantees an unselected target
plaquette at that edge.  At a proper closed partial set, connectedness of \(S\)
guarantees an edge-adjacent target plaquette.  In either case the required
branch is present.  Completed sets are filtered by an exterior flood fill; the
recovered bounded cells must contain the origin and have boundary exactly equal
to the selected plaquettes.

**[COMPUTATION]** The direct search visited 240439 distinct partial plaquette
sets and found 67 closed sets before the enclosure filter.  It independently
returns

\[
 (N(6),N(8),N(10))=(1,0,6),
\]

agreeing with the polycube census.

## 4. A proved all-area tail

### 4.1 A root on a fixed ray

Let \(\Gamma=\partial^*P\) have area \(A\) and contain the origin cell.  The
intersection \(P\cap\{(j,0,0):j\in\mathbb Z\}\) contains \(0\) and is finite.
Let \(m\ge0\) be the largest nonnegative integer for which every
\((j,0,0)\), \(0\le j\le m\), lies in \(P\).  Then the positive face of
\((m,0,0)\), centred at \((m+\tfrac12,0,0)\), belongs to \(\Gamma\).  Likewise,
if \(r\le0\) is the smallest integer for which every \((j,0,0)\),
\(r\le j\le0\), lies in \(P\), the negative face of \((r,0,0)\), centred at
\((r-\tfrac12,0,0)\), belongs to \(\Gamma\).

The contour is edge-connected by the lemma in Section 2.  A simple adjacency
path between these two plaquettes has at most \(A-1\) steps, and one
edge-sharing step changes a plaquette centre's \(x\)-coordinate by at most one.
Their coordinate separation is \(m-r+1\ge m+1\), so \(0\le m\le A-2\).
Thus the selected positive-ray root has fewer than \(A\) possible positions.
We use the convenient overcount \(A\).

### 4.2 Connected plaquette animals via the universal cover

The graph whose vertices are dual plaquettes and whose adjacency is shared-edge
adjacency has degree 12: a plaquette has four edges, and at each edge there are
three other incident plaquettes.  These twelve neighbours are distinct.

Fix an ordering of the twelve neighbour slots.  A connected \(A\)-plaquette set
containing a specified root has a deterministic ordered depth-first spanning
tree: run depth-first search with that slot order and record the rooted tree,
including at every vertex the occupied child slots.  The record reconstructs
the animal, since following a slot from a known plaquette determines its
neighbour.  It is therefore injective.  Ignoring relations in the plaquette
graph embeds these records among the \(A\)-vertex rooted subtrees of the
infinite 12-regular tree; let their number be \(B_A\).

For a planted branch let \(T(z)\) count its vertices, and let \(R(z)\) count a
rooted subtree.  The root has 12 available child slots and every non-root vertex
has 11 after reserving its parent slot, so as formal power series

\[
 T=z(1+T)^{11},\qquad R=z(1+T)^{12}.                             \tag{4.1}
\]

For reference, the Lagrange coefficient identity used here follows directly
from formal residues: substituting \(z=u/\phi(u)\) in
\(\operatorname{Res} f(T(z))z^{-m-1}\,dz\), and collecting the derivative term,
gives

\[
 [z^m]f(T(z))=\frac1m[u^{m-1}]f'(u)\phi(u)^m
 \quad\text{when }T=z\phi(T).
\]

Apply it to \(f(u)=(1+u)^{12}\), \(\phi(u)=(1+u)^{11}\), and \(m=A-1\).  Then
\(B_1=1\) and, for \(A\ge2\),

\[
 B_A=\frac{12}{A-1}\binom{11A}{A-2}.                            \tag{4.2}
\]

No asymptotic lattice-animal constant is imported here.

### 4.3 Explicit exponential majorant and geometric tail

Put

\[
 \lambda:=\frac{11^{11}}{10^{10}}
 =\frac{285311670611}{10000000000}=28.5311670611.
\]

The binomial theorem at \(t=1/10\) gives, term by term,

\[
 \binom{11A}{A-2}10^{-(A-2)}
 \le(1+1/10)^{11A},
\]

and hence the exact inequality

\[
 \binom{11A}{A-2}\le\frac{\lambda^A}{100}.                     \tag{4.3}
\]

Combining the at-most-\(A\) root choices with (4.2)--(4.3),

\[
 N(A)\le A B_A
 \le\frac{12A}{100(A-1)}\lambda^A.                             \tag{4.4}
\]

For the unenumerated range \(A\ge26\), \(A/(A-1)\le26/25\), so

\[
 N(A)\le\frac{78}{625}\lambda^A.                               \tag{4.5}
\]

All areas are even by (2.1).  Writing \(y=e^{-2K}\) and assuming
\(\lambda y<1\), the entire tail is therefore bounded by the exact geometric
sum

\[
 \sum_{\substack{A>24\\A\text{ even}}}N(A)y^A
 \le \frac{78}{625}
       \frac{(\lambda y)^{26}}{1-(\lambda y)^2}.                \tag{4.6}
\]

This improves substantially on the elementary depth-first-walk estimate
\(12^{2A}=144^A\), while remaining fully explicit and proved.

## 5. Exact rational optimisation and certificate

From (3.3) and (4.6), define

\[
\begin{aligned}
 H(y)={}&y^6+6y^{10}+45y^{14}+12y^{16}+332y^{18}
          +240y^{20}+2538y^{22}+3040y^{24},\\
 T(y)={}&\frac{78}{625}\frac{(\lambda y)^{26}}{1-(\lambda y)^2}.
                                                                    \tag{5.1}
\end{aligned}
\]

Both are increasing for \(0\le\lambda y<1\).  At
\(K_0=848449/500000\), exact alternating Taylor sums through orders 60 and 61,
followed by outward rounding, give

\[
 \frac{3358096131959342840851141423768071469788221606373}{10^{50}}
 \le e^{-2K_0}\le
 \frac{3358096131959342840851141423768071469788221606374}{10^{50}}.
                                                                    \tag{5.2}
\]

Why (5.2) is rigorous: \(2K_0=848449/250000\); after order 60 the terms
\((2K_0)^j/j!\) decrease, and the alternating-series theorem places
\(e^{-2K_0}\) between the odd partial sum \(S_{61}\) and the even partial sum
\(S_{60}\).  Every term and both decimal comparisons are evaluated as exact
fractions.

Substituting the **upper** rational endpoint of (5.2) in (5.1) gives the outward
intervals

\[
\begin{aligned}
 H(y)&<0.0000000014340426259949935895833213413351,\\
 T(y)&<0.4999604373722501023948481567957125622441,\\
 H(y)+T(y)&<0.4999604388062927283898417463790339035792<\tfrac12.
                                                                    \tag{5.3}
\end{aligned}
\]

The exact rational values and the positive exact rational margin are stored in
`results/bounds/peierls_upper.json`; (5.3) is only their outward decimal
rendering.  The margin obeys

\[
 \tfrac12-H(y)-T(y)
 >0.0000395611937072716101582536209660964208.                    \tag{5.4}
\]

Equations (2.5), (5.1), and (5.3) prove the theorem in Section 1.

**[COMPUTATION]** Applying the same exact lower enclosure to the preceding
millionth gives

\[
 H(e^{-2(1.696897)})+T(e^{-2(1.696897)})
 >0.5000088171994883731344704944018403797616.                   \tag{5.5}
\]

Thus \(1.696898\) is the least point on the \(10^{-6}\) grid certified by this
specific head-plus-tail majorant.  No claim of a least rational over all
possible denominators is made.

**[THEOREM] Certified limitation of this tail.** Formula (4.6) requires
\(\lambda e^{-2K}<1\), or \(K>\tfrac12\log\lambda\).  In particular this
threshold is greater than 1: \(\lambda>9\), while
\(e=\sum_{j\ge0}1/j!<1+1+\sum_{j\ge2}2^{-(j-1)}=3\), so
\(\tfrac12\log\lambda>\log3>1\).  Since the incumbent is below 1, replacing any
finite number of terms of this geometric tail by exact counts cannot make this
tail argument improve the incumbent.  A genuinely sharper all-area contour
bound, not merely a longer finite census, would be required.

## 6. External sources and relation to prior work

- **[EXTERNAL, abstract only]** R. Peierls, *On Ising's model of
  ferromagnetism*, Math. Proc. Cambridge Philos. Soc. **32** (1936),
  DOI `10.1017/S0305004100019174`.  Crossref supplied the abstract and metadata;
  no full text from this route is used in the proof.
- **[EXTERNAL, metadata only]** R. B. Griffiths, *Peierls Proof of Spontaneous
  Magnetization in a Two-Dimensional Ising Ferromagnet*, Phys. Rev. **136**
  (1964), A437--A439, DOI `10.1103/PhysRev.136.A437`.  This is historical
  context, not a hidden proof dependency.
- **[EXTERNAL, full source]** S. Friedli and Y. Velenik, *Statistical Mechanics
  of Lattice Systems: A Concrete Mathematical Introduction*, Chapter 3,
  Section 3.7.2, especially Lemma 3.37 and Exercise 3.20.  The freely supplied
  final-draft chapter is cached at
  `sources/fulltext/friedli_velenik_ising_ch3.pdf`, SHA-256
  `4809d134deb0a97bb37bcc306345d3ace659f7543bbcd8cfedde53c5dbe47045`.
  It gives the modern rigorous finite-volume flip argument and points explicitly
  to the higher-dimensional extension.  The constants used here are proved
  independently above.
- **[EXTERNAL, full source; comparison only]** C. Bonati, *The Peierls argument
  for higher dimensional Ising models*, Eur. J. Phys. **35** (2014) 035002,
  arXiv:1401.7894, cached SHA-256
  `2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36`.
  Its specialised three-dimensional count gives the previously recorded
  \(K_c\le0.7488385\ldots\), still weaker than the infrared bound but stronger
  than the deliberately elementary all-plaquette-animal tail (4.6).

The source metadata and access status are recorded in `sources/manifest.yaml`.

## 7. Scope, reproducibility, and claim classification

- **[THEOREM]** Equations (2.2)--(2.6) prove the plus-boundary Peierls
  inequality with the outer contour of the origin minus component.
- **[LEMMA]** The cubical closedness, edge-connectedness, parity of area,
  discrete isoperimetric inequality, universal-cover animal count, and explicit
  tail (4.6) are proved above.
- **[COMPUTATION]** The table (3.3) is an exhaustive exact computation for the
  finite scope \(A\le24\), made complete by (3.1).  It is independently checked
  through \(A=10\) by direct plaquette closure search.
- **[THEOREM plus exact rational computation]** The strict rational inequality
  (5.3) proves \(K_c\le848449/500000\).
- **[COMPUTATION / negative result]** This bound is weaker than the incumbent;
  no improvement is claimed.  The benchmark critical-coupling estimate was not
  used for selection, fitting, or validation.

Reproduce from the repository root:

```text
.venv/bin/python experiments/e41_peierls_upper.py
.venv/bin/python tests/test_peierls_upper.py
```
