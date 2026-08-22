# Exact modular lower certificate and symmetry container for the open 4x4 layer

## Scope

Let \(\Lambda=\{0,1,2,3\}^2\) be the finite **open** square grid, with
24 nearest-neighbour edges, and define

\[
 A=\sum_{v\in\Lambda}X_v,\qquad
 B=\sum_{\{u,v\}\in E(\Lambda)}Z_uZ_v,
 \qquad \mathfrak g=\operatorname{Lie}_{\mathbb Q}\langle A,B\rangle.
\]

Every statement below concerns this one finite layer only.  It makes no claim
about a different layer, all layer sizes, or the thermodynamic 3D Ising model.

## Literal-word modular lower certificate

**[LEMMA]** Let \(w_1,\ldots,w_r\) be nested commutator words in integral
Pauli matrices.  If their coordinate matrix has rank \(r\) modulo an odd
prime \(p\), then the same words are linearly independent over \(\mathbb Q\).

Indeed, a nonzero \(r\)-by-\(r\) minor modulo \(p\) is an integer minor that
is nonzero over \(\mathbb Q\).  Dividing each commutator by \(2\) changes no
span over either \(\mathbb Q\) or an odd \(\mathbb F_p\).

**[COMPUTATION]** At
\(p=2{,}147{,}483{,}647\), a sparse ordered-monomial echelon computation
recorded 3,609.162794 total seconds and stopped its closure loop at the
3,600-second guard during depth 16, after constructing 1,794 independent
**literal** nested raw Pauli words.  Therefore

\[
 \dim_{\mathbb Q}\mathfrak g\geq 1794.
\]

This is a lower bound only; it is not a characteristic-zero dimension equality.
The stored JSON checkpoint gives the parent/action/depth recipe for every
accepted word, so the minor can be reconstructed without importing the
producer experiment.

**[UNRESOLVED]** The closure did not saturate, so no second-prime closure was
run and no two-prime equality claim is made.  The single good-prime minor is
already a rigorous lower bound, but it is not an exact-\(\mathbb Q\) squeeze.

**[COMPUTATION]** The last fully completed breadth-first depth is 15.  Its
exact cumulative rank profile is

\[
(2,3,5,7,11,16,26,41,67,105,169,271,448,741,1240),
\]

and its four Pauli grading-block ranks are

\[
(393,221,391,235).
\]
At the stopped partial depth 16, 634 of the 998 parent/action expansions
from the 499-word depth-15 frontier had been processed; 554 were independent.
The current grading-block ranks are \((393,494,391,516)\).  The completed-depth
checkpoint consists of 1,240 words; the partial certificate has 1,794 words,
maximum literal depth 16, and tree checksum

\[
\sum_{i=1}^{1794}i(p_i+2)(a_i+2)(d_i+1)
 \equiv 898915470\pmod{1{,}000{,}000{,}007}.
\]

The producer records the literal recipes themselves, not merely this checksum;
the checksum is diagnostic rather than the certificate.

**[COMPUTATION]** The sparse run discovered 378,843 exact \(D_4\)-orbits of
Pauli support, used a conservative sparse-basis estimate of 4,312,089,696
bytes, and reached 10,666,885,120 bytes RSS.  An independent support BFS also
exceeded 250,000 \(D_4\)-orbits, reaching 250,001, in 2.386095 seconds.  This
excludes the dense all-support route at that cap; it does not estimate the
final reachable support size.

## Exact finite-group sector data

**[THEOREM]** The commuting subgroup

\[
 C_2(\text{row reflection})\times C_2(\text{column reflection})
 \times C_2(\text{global spin flip})
\]

has exact character-projector ranks, in character order
\(000,001,010,011,100,101,110,111\),

\[
(8384,8192,8128,8192,8128,8192,8128,8192).
\]

This follows directly from the finite character-trace formula; no numerical
linear algebra is used.

**[THEOREM]** Resolving the full \(D_4\times C_2(\text{global flip})\)
symmetry gives the following exact multiplicity-space dimensions:

\[
\begin{array}{c|ccccc}
 &A_1&A_2&B_1&B_2&E\\\hline
 +&4324&3940&4060&4188&8128\\
 -&4224&3968&3968&4224&8192.
\end{array}
\]

Here the two-dimensional \(E\) representation contributes twice its
multiplicity to the physical-module dimension.  Since both generators commute
with this finite symmetry, the natural representation gives the exact
same-basis rational containment

\[
 \mathfrak g\subseteq
 \bigoplus_{k\in K}\operatorname{End}_{\mathbb Q}(\mathbb Q^k),
\]

where

\[
 K=(4324,3940,4060,4188,8128,4224,3968,3968,4224,8192).
\]

Its dimension is

\[
 \sum_{k\in K}k^2=268{,}591{,}168.
\]

This is an exact rational symmetry container, but it is far too broad to
squeeze the 1,697-dimensional lower certificate.

## What the sector data do not prove

**[CONJECTURE]** If the dynamic commutant-centre projectors of the generated
pair \((A,B)\) were exactly these ten multiplicity spaces and if their
non-scalar actions were linked as full split special-linear factors, the
resulting one-scalar candidate would have dimension

\[
1+\sum_{k\in K}(k^2-1)=268{,}591{,}159.
\]

This is only a symmetry-guided candidate, not a factor classification.  The
3x3 theorem demonstrates that symmetry-isotypic spaces can be further refined
or linked by the **dynamic** commutant centre, so the 4x4 factor dimensions
cannot be inferred from these ranks.

**[UNRESOLVED]** No exact rational dynamic commutant-centre calculation,
linked-sector projector construction, or tight same-basis \(\mathbb Q\) upper
container was materialized for 4x4.  In particular, neither
\(\dim_{\mathbb Q}\mathfrak g\) nor a Levi type, centre, radical, or Killing
rank is claimed.

## Clean-room verification

```bash
timeout 9000 .venv/bin/python tests/test_algebra_4x4.py
```

The verifier imports neither `experiments/e97_algebra_4x4.py` nor its closure
helpers.  It independently rebuilds the finite-group traces, the capped
support BFS, all stored literal words and their modular pivots, and the
completed depth-15 closure using a reversed seed/frontier/action order.  The
observed final line was `WALL_SECONDS=3731.169` followed by `PASS`.
