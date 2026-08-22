# Independent sixth-order slab certificate and a closed correlation form for \(c_4\)

Artifact: `results/interlayer/closed_form.json`  
Producer: `experiments/e84_interlayer_closed.py`  
Clean-room verifier: `tests/test_interlayer_closed.py`

All formal-series coefficients below are computed using Python `int` and
`fractions.Fraction`; no floating-point arithmetic enters either certificate.
Write \(v=\tanh K\) for the in-plane variable, \(q=K_z\) for the direct
interlayer coupling, and \(w=\tanh q\).

## 1. Results and scope

[COMPUTATION] An independent open-slab even-subgraph finite-lattice calculation
reproduces the previously primary/generated-cumulant value of the residual
\([w^6]\) coefficient through every in-plane order \(v^0,\ldots,v^{12}\):

\[
\boxed{
 [w^6]\log P(v,w)
 =2v^2+138v^4+\frac{13682}{3}v^6+109718v^8
   +2061698v^{10}+32654478v^{12}+O(v^{14}).
}
\]

Here \(P(v,w)\) is the high-temperature even-subgraph polynomial, so this is
exactly `c6_w_residual`, i.e. the total coefficient after removal of the
separate vertical \(\log\cosh K_z\) prefactor.  The two routes agree
coefficientwise, including the three formerly single-route entries at
\(v^8,v^{10},v^{12}\).

[THEOREM] Let the zero-field two-dimensional layer correlations be

\[
 G_{ij}=\langle\sigma_i\sigma_j\rangle_{2\mathrm D},\qquad
 U_{ijkl}=\langle\sigma_i\sigma_j\sigma_k\sigma_l\rangle_{2\mathrm D}
 -G_{ij}G_{kl}-G_{ik}G_{jl}-G_{il}G_{jk}.
\]

Then, as a coefficientwise high-temperature identity,

\[
\boxed{
\begin{aligned}
 c_{4,q}={}&\frac1{24}\sum_{r,s,t}U_{0rst}^{\,2}
 +\frac12\sum_{r,s,t}G_{0r}G_{st}U_{0rst}\\
 &+\frac34\sum_{r,s,t}G_{0r}G_{st}G_{0s}G_{rt}.
\end{aligned}}
\tag{1}
\]

The total fourth coefficient in the repository variable \(w\) is consequently

\[
\boxed{
 c_{4,w}^{\rm tot}=c_{4,q}+\frac13\sum_rG_{0r}^{\,2}.}
\tag{2}
\]

The sums in (1) and (2) are first finite-volume identities.  In the infinite
lattice they mean formal high-temperature series: at every fixed coefficient
only finitely many connected graphs contribute.  No convergence claim is made
near the critical point.

[UNRESOLVED] Neither result gives a radius of convergence, evaluates an
infinite sum at a critical coupling, or solves the three-dimensional Ising
model.

## 2. Exact finite-support certificate for \(c_6\)

### 2.1 Planar boundary-polynomial transfer

For an open \(a\times b\) layer and boundary mask \(S\), define

\[
 P_S(v)=\sum_{E\subseteq E_{a\times b}:\,\partial E=S}v^{|E|},
 \qquad m_S(v)=\frac{P_S(v)}{P_\varnothing(v)}.
\]

The producer obtains every \(P_S\) by an edge-by-edge transfer indexed by the
odd-degree mask.  Since \(P_\varnothing(0)=1\), the displayed division is a
formal division in \(\mathbb Z[[v]]\), performed coefficientwise through
\(v^{12}\).  Thus the computation obtains exact finite-rectangle two-point
and four-point layer moments without a numerical transfer matrix.

For an open stack of \(h\) layers, let \(S_g\) be the selected vertical-edge
mask on the gap from layer \(g\) to \(g+1\), with \(S_0=S_h=\varnothing\).
The exact normalized slab polynomial is

\[
\frac{P_h(v,w)}{P_\varnothing(v)^h}
 =\sum_{S_1,\ldots,S_{h-1}}
   w^{\sum_g|S_g|}\prod_{\ell=0}^{h-1}m_{S_\ell\triangle S_{\ell+1}}(v).
\tag{3}
\]

A nonzero planar factor in (3) has even boundary cardinality.  Starting from
\(S_0=\varnothing\), it follows inductively that every \(S_g\) has even
cardinality.  Hence all odd \(w\)-columns vanish identically before the
calculation begins.

### 2.2 Fast but exact \(w^6\) extraction

Let

\[
 Q_d=\sum_{|S|=d}m_S^2,\qquad
 F_{24}=\sum_{\substack{|A|=2\\|B|=4}}m_A m_{A\triangle B}m_B,
\]

and

\[
 H_{222}=\sum_{|B|=2}\left(\sum_{|A|=2}m_A m_{A\triangle B}\right)^2.
\]

These are finite exact sums of integer truncated series.  Enumeration of the
possible weights of the gap masks in (3) gives the following \(w^6\) columns:

\[
 r^{(2)}_6=Q_6,
 \qquad
 r^{(3)}_6=2Q_6+2F_{24},
\]

and

\[
 r^{(4)}_6=3Q_6+4F_{24}+2Q_2Q_4+H_{222}.
\tag{4}
\]

The last term is the only three-gap case: all three masks have weight two, and
its factorization is the identity

\[
 \sum_{|A|=|B|=|C|=2}
 m_A m_{A\triangle B}m_{B\triangle C}m_C
 =\sum_{|B|=2}\left(\sum_{|A|=2}m_A m_{A\triangle B}\right)^2.
\]

Thus (4) is a rearrangement of the full vertical-mask sum, not a truncation or
sampling.  With
\(R=1+r_2w^2+r_4w^4+r_6w^6+O(w^8)\), the exact logarithm is

\[
 [w^6]\log R=r_6-r_2r_4+\frac13r_2^3.
\tag{5}
\]

### 2.3 Why the box family is complete

[LEMMA] Every connected even subgraph contributing to
\([v^d w^6]\log P\) with \(d\le12\) fits in an open box with at most four
layers and with in-plane bounding rectangle \(a\times b\) satisfying

\[
 (a-1)+(b-1)\le6.
\tag{6}
\]

*Proof.* Each nonempty vertical cut has an even positive number of vertical
edges.  Six vertical edges can therefore occupy at most three gaps, hence at
most four layers.  If a connected even subgraph spans an in-plane coordinate
cut, it crosses that cut a positive even number of times.  It crosses each of
the \(a-1\) horizontal and \(b-1\) vertical bounding-box cuts at least twice,
using at least \(2[(a-1)+(b-1)]\) in-plane edges.  For \(d\le12\), this proves
(6). \(\square\)

The producer applies ordinary three-dimensional rectangular Möbius inversion
to all 112 open boxes satisfying (6) and heights one through four.  The log
selects connected polymer weights, and the inversion removes every proper
subbox embedding.  The lemma proves that a box omitted by this family has
either vertical degree at least eight or in-plane degree at least fourteen.
Its \([v^0,\ldots,v^{12}]w^6\) weight therefore vanishes exactly.  This is an
extrapolation-free finite-support identification: enlarging a side beyond the
bounding box cannot create a missing graph at the certified coefficient.

The exact vertical-height components of the certified residual are

\[
\begin{array}{c|rrrrrr}
\text{height}&v^2&v^4&v^6&v^8&v^{10}&v^{12}\\ \hline
2&0&0&56/3&494&7528&96904\\
3&0&-20&-560&-5472&16160&1680036\\
4&2&158&5102&114696&2038010&30877538.
\end{array}
\]

Their columnwise sum is the boxed \(c_6\) series in Section 1.  The full
per-box certificate is retained in
`results/interlayer/closed_form.json` together with SHA-256
`1f4ce1702faab90e3d805013aefd5dc77e4d2f12cb1275c01e7ac2ea2e7e439a`.

## 3. Derivation of the \(c_4\) correlation identity

Work first in a finite translation-invariant layer volume; all sums below are
finite.  Pin the first in-plane position at \(0\), and abbreviate

\[
\begin{aligned}
 M&=\langle\sigma_0\sigma_r\sigma_s\sigma_t\rangle,\\
 A&=G_{0r}G_{st},\qquad
 B=G_{0s}G_{rt},\qquad
 C=G_{0t}G_{rs},\\
 U&=M-A-B-C.
\end{aligned}
\]

The direct-coupling cumulant expansion of independent layers has two allowed
fourth-order gap patterns.  The all-four-on-one-gap contribution and the
two-on-each-of-two-adjacent-gaps contribution are, respectively,

\[
 \kappa_{\rm same}=M^2-A^2-B^2-C^2,
 \qquad
 \kappa_{\rm adjacent}=AM-A^2.
\tag{7}
\]

Their placement weights are \(1/4!\) and \(1/(2!2!)=1/4\), so

\[
 c_{4,q}=\frac1{24}\sum_{r,s,t}\kappa_{\rm same}
          +\frac14\sum_{r,s,t}\kappa_{\rm adjacent}.
\tag{8}
\]

Equation (7) follows directly from the fourth joint-cumulant partition
formula.  For one gap, the three pair partitions give \(A^2,B^2,C^2\).  For
the adjacent-gap pattern, outer-layer parity kills the two crossed pair
partitions, leaving only \(AM-A^2\).  This is the exact layer-cumulant
origin of both terms; no coefficient has been fitted.

Substitute \(M=A+B+C+U\) into (7):

\[
\begin{aligned}
 \kappa_{\rm same}
 &=U^2+2U(A+B+C)+2(AB+AC+BC),\\
 \kappa_{\rm adjacent}&=A(U+B+C).
\end{aligned}
\tag{9}
\]

Relabelling the dummy indices \(r,s,t\) leaves the finite sum invariant and
gives

\[
 \sum AU=\sum BU=\sum CU,
 \qquad
 \sum AB=\sum AC=\sum BC.
\tag{10}
\]

Putting (9) and (10) in (8) yields

\[
 c_{4,q}=\frac1{24}\sum U^2+\frac12\sum AU+\frac34\sum AB,
\]

which is (1).  Thus the identity is finite-volume algebra.  Passing to the
infinite lattice coefficientwise is legitimate in the formal
high-temperature expansion: a fixed coefficient has finite graph support.

For the \(w\) statement, write

\[
 F(q)=c_{2,q}q^2+c_{4,q}q^4+O(q^6),
 \qquad q=\operatorname{atanh}w=w+\frac13w^3+O(w^5).
\]

Then

\[
 [w^4]F(\operatorname{atanh}w)=c_{4,q}+\frac23c_{2,q}.
\]

The second cumulant is
\(c_{2,q}=\frac12\sum_rG_{0r}^2\), proving (2).

## 4. Exact coefficient verification of the closed form

The same integer boundary-polynomial transfer computes every finite-rectangle
\(G\), \(M\), and \(U\) needed in (1) through \(v^{12}\).  Planar FLM over
all 28 rectangles with \((a-1)+(b-1)\le6\) gives

\[
\begin{aligned}
 c_{2,q}={}&\frac12+2v^2+18v^4+118v^6+778v^8+4978v^{10}+31398v^{12},\\
 c_{4,q}={}&-\frac1{12}+\frac23v^2+51v^4+\frac{2914}{3}v^6
 +\frac{41113}{3}v^8+\frac{481894}{3}v^{10}+1679961v^{12},\\
 c_{4,w}^{\rm tot}={}&\frac14+2v^2+63v^4+1050v^6+14223v^8
 +163950v^{10}+1700893v^{12}.
\end{aligned}
\]

The three terms in (1), in their displayed order, are

\[
\begin{array}{c|rrrrrrr}
&v^0&v^2&v^4&v^6&v^8&v^{10}&v^{12}\\ \hline
\frac1{24}\sum U^2
&1/6&8/3&38&1336/3&14374/3&143080/3&447778\\
\frac12\sum AU
&-1&-20&-356&-4652&-53628&-558084&-5400268\\
\frac34\sum AB
&3/4&18&369&5178&62541&671022&6632451.
\end{array}
\]

They sum coefficientwise to the stated \(c_{4,q}\).  As elementary transfer
witnesses, on the open \(3\times3\) rectangle the producer records

\[
 G_{0,1}=v+v^3+v^5-v^7-13v^9-7v^{11}+O(v^{13}),
\]

and

\[
 U_{0,1,3,4}=-8v^4-16v^6+104v^{10}+232v^{12}+O(v^{13}).
\]

These witnesses are finite-rectangle data only; the infinite-volume values
above are the separately certified FLM contractions.

## 5. Reproduction

```sh
timeout 600 .venv/bin/python experiments/e84_interlayer_closed.py
timeout 600 .venv/bin/python tests/test_interlayer_closed.py
```

Both commands print `PASS` only after exact checks succeed.  The test imports
neither the producer nor its helpers: it independently rebuilds the boundary
transfer, the slab Möbius inversion, and the \(G/U_4\) contractions; it also
matches the normalized \(2\times2\times4\) slab columns against the separate
integer spin-DOS transform.
