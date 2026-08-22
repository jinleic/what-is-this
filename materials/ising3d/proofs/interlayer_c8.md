# Exact eighth interlayer coefficient through \(v^{12}\)

Artifact: `results/interlayer/c8_series.json`  
Producer: `experiments/e126_interlayer_c8.py`  
Clean-room verifier: `tests/test_interlayer_c8.py`

[COMPUTATION] All displayed finite-series coefficients and every intermediate
calculation in these three files use Python integers or `fractions.Fraction`;
no floating point calculation contributes to a mathematical claim.

## 1. Conventions and statement

Put
\[
v=\tanh K,\qquad q=K_z,\qquad w=\tanh q,
\]
and use the repository convention
\[
f_{3\mathrm D}(K,w)=f_{2\mathrm D}(K)+
 \sum_{m\geq1}c_{2m}^{\mathrm{tot}}(v)w^{2m}.
\]
[LEMMA] The separate vertical high-temperature prefactor is
\[
\log\cosh(\operatorname{atanh}w)=-\tfrac12\log(1-w^2),
\]
so \([w^8]\log\cosh(\operatorname{atanh}w)=1/8\).  Below,
\(c_{8}^{\mathrm{res}}=c_{8}^{\mathrm{tot}}-1/8\) denotes the coefficient
with that one-bond prefactor removed.

[COMPUTATION] Exact rectangular finite-lattice inversion gives
\[
\boxed{
 c_{8}^{\mathrm{res}}(v)=
 2v^2+243v^4+13594v^6+\frac{1046015}{2}v^8
 +14894834v^{10}+341377607v^{12}+O(v^{14}).
}
\]
Equivalently,
\[
\boxed{
 c_{8}^{\mathrm{tot}}(v)=\frac18+2v^2+243v^4+13594v^6
 +\frac{1046015}{2}v^8+14894834v^{10}+341377607v^{12}+O(v^{14}).
}
\]
This is a finite, coefficientwise formal-series result only.

For auditability, the direct-coupling calculation preceding the exact change
of variable is
\[
\begin{aligned}
 c_{8,q}(v)={}&-\frac{17}{2520}+\frac{2}{315}v^2
 +\frac{1367}{35}v^4+\frac{1798474}{315}v^6\\
 &+\frac{201906161}{630}v^8+\frac{3454667434}{315}v^{10}
 +\frac{29200871981}{105}v^{12}+O(v^{14}).
\end{aligned}
\]

## 2. Structural controls at eighth order

[THEOREM] `proofs/interlayer_allorders.md`, Theorems 1 and 3, applies to this
open stack: odd interlayer coefficients vanish; at order \(2m\), only
consecutive gap patterns with even gap multiplicities contribute; and no more
than \(m+1\) layers occur.  For \(m=4\), this permits exactly the eight
compositions
\[
(4),\ (3,1),\ (2,2),\ (2,1,1),\ (1,3),\ (1,2,1),\ (1,1,2),\ (1,1,1,1).
\]
Here a part \(a_j\) denotes multiplicity \(2a_j\) in gap \(j-1\).  The
producer derives, rather than inserts, its bulk placement factor: it explicitly
enumerates the labeled assignments of the repeated gap multiset and divides
their number by \(8!\).  This gives
\(1/\prod_j(2a_j)!\).  The patterns use respectively \(2,3,3,4,3,4,4,5\)
layers; hence a height-six \(c_8\) weight must vanish.

[COMPUTATION] The producer independently evaluated the complete Bell(8)
partition sum and found zero primary cumulants \(c_1,c_3,c_5,c_7\) through
\(v^{12}\), zero odd \(v\)-degrees in \(c_{8,q}\), and an identically zero
height-six \(c_8\) rectangle weight.  Separately, the clean-room verifier
literally evaluates the layer-labelled Bell/moment polynomial for every one of
the 3,432 normalized eight-gap multisets, without pre-filtering for contiguity
or even multiplicity.  Its exactly eight nonzero profiles are precisely the
eight displayed compositions; it also exhausts all 4,140 partitions of eight
labels.

[THEOREM] The all-order ladder-cycle identity from
`proofs/interlayer_allorders.md` gives
\([v^2w^{2m}]c_{2m}^{\mathrm{res}}=2\) for every \(m\ge1\).  It was not
inserted into the calculation: the independent \(w\)-column FLM reconstruction
returns
\[
[v^2w^8]c_{8}^{\mathrm{res}}=2.
\]

## 3. Generated-cumulant route and the \(q\)-to-\(w\) conversion

[LEMMA] For vertical bonds
\(Y_{g,x}=\sigma_{x,g}\sigma_{x,g+1}\), the finite-volume expansion is
\[
\frac1{|\Lambda|}\log\langle e^{qV}\rangle_0
 =\sum_{n\ge1}\frac{q^n}{n!|\Lambda|}
   \sum_{b_1,\ldots,b_n}\kappa_0(Y_{b_1},\ldots,Y_{b_n}).
\]
Every partition block factorizes into zero-field 2D moments after XOR
(parity) reduction of repeated in-plane positions.  This is the established
all-order formula used in `experiments/e75_interlayer_allorders.py`; the
present producer evaluates it at \(n=8\) with exact boundary-mask layer
moments and rectangular finite-lattice inversion.

[COMPUTATION] The eight pattern rows have respectively:

| composition | ordered gap assignments | placement | surviving Bell partitions | collected raw moment monomials | distinct connected monomials |
| --- | ---: | ---: | ---: | ---: | ---: |
| \((4)\) | 1 | \(1/40320\) | 379 | 379 | 56,932 |
| \((3,1)\) | 28 | \(1/1440\) | 107 | 107 | 8,718 |
| \((2,2)\) | 70 | \(1/576\) | 83 | 83 | 4,458 |
| \((2,1,1)\) | 420 | \(1/96\) | 35 | 22 | 270 |
| \((1,3)\) | 28 | \(1/1440\) | 107 | 107 | 8,718 |
| \((1,2,1)\) | 420 | \(1/96\) | 35 | 31 | 621 |
| \((1,1,2)\) | 420 | \(1/96\) | 35 | 22 | 270 |
| \((1,1,1,1)\) | 2520 | \(1/16\) | 15 | 8 | 27 |

Thus the exact connected inventory has 796 admissible partitions and 80,014
distinct slot-indexed monomials before the harmless relabeling-orbit collapse.
Only \(G\) (two-point), \(U\) (connected four-point), \(W\) (connected
six-point), and \(W8\) (connected eight-point) 2D correlation objects occur.
The producer verifies the collection two ways (substitute before versus after
raw collection), even slot incidence, anchor connectivity, reflection pairing,
and regression to the already certified sixth-order 13-term formula.

[COMPUTATION] In addition to its counts, `c8_series.json` stores a SHA-256
digest for each complete placement-scaled coefficient map and the collection
digest
`86587088de826e47e404feec36e115be6afa1a418eb4673a94ae5e877e7d28dc`.
The standalone verifier separately rebuilds all eight maps, their atom-kind
histograms and Möbius-by-block-class maps, and rejects the artifact unless all
nine digests agree.  This inventory is **not** asserted to be a manageable
scalar closed form.

[LEMMA] Formal integration of
\(d(\operatorname{atanh}w)/dw=(1-w^2)^{-1}\) gives, exactly,
\[
[w^8]q^2=\frac{44}{105},\qquad [w^8]q^4=\frac{22}{15},\qquad
[w^8]q^6=2,\qquad [w^8]q^8=1.
\]
Therefore the conversion checked by exact series arithmetic is
\[
c_{8}^{\mathrm{tot}}
 =c_{8,q}+2c_{6,q}+\frac{22}{15}c_{4,q}+\frac{44}{105}c_{2,q}.
\]
At \(v=0\), this sends
\(-17/2520\) to \(1/8\), as required by the vertical \(\log\cosh\)
prefactor; after removing it, the residual constant term is zero.

## 4. Independent open-slab even-subgraph route

[LEMMA] For an open \(a\times b\) layer, let \(P_S(v)\) be the integer
high-temperature polynomial whose selected in-plane edges have odd boundary
mask \(S\).  With \(S_0=S_h=\varnothing\), the direct even-subgraph column
of an \(a\times b\times h\) slab is
\[
 C_{h,d}(v)=
 \sum_{\substack{S_1,\ldots,S_{h-1}\\
                  \sum_{j=1}^{h-1}|S_j|=d}}
 P_{S_1}(v)\,P_{S_1\triangle S_2}(v)\cdots
 P_{S_{h-2}\triangle S_{h-1}}(v)\,P_{S_{h-1}}(v).
\]
This is a raw even-subgraph identity: the vertical mask on each cut is exactly
the odd boundary which its two adjacent layers must supply.  Odd \(d\) columns
vanish because every layer cut is crossed an even number of times.

[COMPUTATION] The second route constructs these columns through \(d=8\),
computes \([w^d]\log(C_h/C_{h,0})\) with the exact Newton recursion, then
performs three-dimensional box Möbius inversion over the 28 ordered open
rectangles satisfying \((a-1)+(b-1)\le6\) and heights \(1\) through \(6\).
It agrees coefficientwise with the converted generated-cumulant residual
through every displayed degree \(v^0,\ldots,v^{12}\), including every exact
height component.  Heights one and six have zero \(w^8\) weight.

[LEMMA] A connected even subgraph spanning an \(a\times b\) in-plane
rectangle needs at least \(2((a-1)+(b-1))\) in-plane edges.  Consequently the
span-six rectangle family is complete through \(v^{12}\).  The producer also
computes the otherwise unnecessary span-seven rectangles \(8\times1\) and
\(7\times2\); all their exact-height cumulant weights vanish through
\(v^{12}\).

[COMPUTATION] A third route transforms the exact joint spin density of states
of every required open 3D box to \(P(v,w)\), takes the same logarithm, and
performs the same 3D Möbius inversion through \(v^6\).  It uses 50 boxes and
returns
\[
0+0v+2v^2+0v^3+243v^4+0v^5+13594v^6,
\]
which agrees with both prior routes through the full overlap.  Its exact
height-two, -three, -four, and -five contributions are respectively
\[
0,\qquad -5v^4-212v^6,\qquad -40v^4-1696v^6,
\qquad 2v^2+288v^4+15502v^6.
\]

[COMPUTATION] For the separately scoped height-six spin-DOS control, the
verifier performs 3D Möbius inversion for each of
\(1\times1\times6\), \(2\times1\times6\), \(3\times1\times6\),
\(4\times1\times6\), and \(2\times2\times6\), and checks every
\(w^2,w^4,w^6,w^8\) weight individually.  All are zero through \(v^6\).
This is a five-shape control-family check, not a claim that every possible
open-box shape has been enumerated.

## 5. Reproduction and resources

[COMPUTATION] Reproduce the producer and its clean-room verifier with

```sh
PYTHONPATH=src .venv/bin/python experiments/e126_interlayer_c8.py
PYTHONPATH=src .venv/bin/python tests/test_interlayer_c8.py
```

The result artifact records a producer wall budget of 1800 seconds and an
8-GiB RSS preflight ceiling.  Its measured producer wall was 335.151 seconds
and its in-process Darwin `ru_maxrss` was 262,012,928 bytes; these are
observed measurements, while the 8-GiB ceiling is explicitly a preflight
resource wall.  The clean-room test measured 212.3 seconds and
`ru_maxrss=267,583,488` bytes under its declared 900-second wall budget.

## 6. Scope

[UNRESOLVED] The result supplies only the stated finite high-temperature
coefficients about \(K=K_z=0\).  It gives no convergence radius, no analytic
continuation, no critical-point evaluation, no sign theorem, and no solution
of the three-dimensional Ising model.  In particular, the \(G,U,W,W8\)
inventory is a coefficientwise formal reduction to 2D connected correlations,
not a closed evaluation of their infinite-lattice sums.
