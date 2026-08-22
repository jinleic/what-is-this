# All-order structure of the interlayer expansion

Artifact: `results/interlayer/allorders.json`.
Generator: `experiments/e75_interlayer_allorders.py`.
Independent control: `tests/test_interlayer_allorders.py`.

This note proves the all-order statements behind the interlayer expansion of a stack of
two-dimensional Ising layers: every odd coefficient vanishes, the even coefficients are given by an
explicit set-partition cumulant formula over two-dimensional spin moments, exactly \(2^{m-1}\)
layer-gap patterns contribute at order \(2m\), at most \(m+1\) consecutive layers participate, and
the "adjacent-gap correction" first seen in \(c_4\) is the first member of that pattern family.
It also proves one genuine all-order coefficient identity, \([v^2w^{2m}]=2\).

**Scope.** Nothing here solves the three-dimensional Ising model. Theorems 1–6 are structural or
low-order-exact; the interlayer series is not summed, no sign or positivity property is proved, and
no radius of convergence is claimed. Section 9 lists the non-results explicitly.

## 1. Setup and conventions

Let \(\Lambda_2\) be a finite graph (the layer) with vertex set \(V\) and edge set \(E_\parallel\),
and let \(L\ge1\). The stack is
\[
 \Lambda=V\times\{0,\dots,L-1\},
\]
with in-plane bonds \(\{(x,l),(y,l)\}\) for \(\{x,y\}\in E_\parallel\) and every layer \(l\), and
interlayer bonds \(\{(x,l),(x,l+1)\}\) for every \(x\in V\) and \(0\le l\le L-2\). The zero-field
Hamiltonian in reduced couplings is
\[
 -\beta H
 =\sum_{l}\sum_{\{x,y\}\in E_\parallel}K_{l,xy}\,\sigma_{x,l}\sigma_{y,l}
 \;+\;q\sum_{l=0}^{L-2}\sum_{x\in V}\sigma_{x,l}\sigma_{x,l+1},
 \qquad \sigma\in\{\pm1\}^{\Lambda}.
\]
Write \(Z_\Lambda(K,q)=\sum_\sigma e^{-\beta H}\), \(f_\Lambda(K,q)=|\Lambda|^{-1}\log Z_\Lambda\),
and, for the simple cubic case, \(v=\tanh K\) in plane and \(w=\tanh q\) between layers. The
interlayer expansion is
\[
 f(K,q)=f_{2\mathrm D}(K)+\sum_{n\ge1}c_n(K)\,q^n .
\]
The reference measure \(\mu_0\) is the product over layers of independent two-dimensional
zero-field Gibbs measures (i.e. \(q=0\)). For a gap index \(g\) and an in-plane site \(x\) put
\[
 Y_{g,x}=\sigma_{x,g}\,\sigma_{x,g+1},
 \qquad
 V=\sum_{g,x}Y_{g,x},
\]
so that \(Z_\Lambda(K,q)=Z_\Lambda(K,0)\,\langle e^{qV}\rangle_{0}\). Two-dimensional moments are
written \(M_l(S)=\langle\prod_{x\in S}\sigma_{x,l}\rangle_{0}\); with a single layer state these do
not depend on \(l\), and \(G(x,y)=M(\{x,y\})\).

## 2. Theorem 1 — all odd coefficients vanish

**Theorem 1 (finite volume).** Fix any finite zero-field stack as in Section 1, with arbitrary
in-plane graph and arbitrary layer-dependent in-plane couplings, whose interlayer bonds join only
layers of opposite parity. Then
\[
 Z_\Lambda(K,q)=Z_\Lambda(K,-q)\quad\text{identically in }q,
\]
hence \(f_\Lambda(K,\cdot)\) is even and real-analytic near \(0\), and every odd Taylor coefficient
of \(f_\Lambda\) at \(q=0\) is exactly zero. The same holds in \(w=\tanh q\).

*Proof.* Define the involution \(T\) on spin configurations by
\(T\sigma_{x,l}=(-1)^l\sigma_{x,l}\). It is a bijection of \(\{\pm1\}^\Lambda\). For an in-plane
bond both endpoints share the layer index, so
\(T\sigma_{x,l}T\sigma_{y,l}=(-1)^{2l}\sigma_{x,l}\sigma_{y,l}=\sigma_{x,l}\sigma_{y,l}\): every
in-plane term is invariant, for every choice of the couplings \(K_{l,xy}\). For an interlayer bond
the two endpoints lie in layers of opposite parity, so
\(T\sigma_{x,l}T\sigma_{x,l'}=(-1)^{l+l'}\sigma_{x,l}\sigma_{x,l'}=-\sigma_{x,l}\sigma_{x,l'}\):
every interlayer term changes sign. Therefore \((-\beta H)(T\sigma;K,q)=(-\beta H)(\sigma;K,-q)\),
and summing over \(\sigma\) gives \(Z_\Lambda(K,q)=Z_\Lambda(K,-q)\). Each \(Z_\Lambda\) is a finite
sum of exponentials, hence entire in \(q\) and strictly positive, so
\(f_\Lambda=|\Lambda|^{-1}\log Z_\Lambda\) is analytic in a complex neighbourhood of \(q=0\) and
even; an even analytic function has vanishing odd Taylor coefficients. Finally \(q\mapsto\tanh q\) is
odd and analytic with \(\tanh'(0)=1\), so \(w\mapsto f_\Lambda(K,\operatorname{atanh}w)\) is even
too. \(\square\)

**Remarks on the exact scope.**

1. *No in-plane bipartiteness is used.* Only the layer-parity two-colouring of the interlayer bonds
   matters. Check `finite slab interlayer evenness` verifies this on a four-layer slab over an
   in-plane triangle (a frustrated, non-bipartite layer): all of
   \([q^1],[q^3],[q^5],[q^7],[q^9]\) of \(Z\) are exactly \(0\).
2. *Zero field is essential.* A term \(h\sum\sigma\) is not invariant under \(T\).
3. *The parity colouring is essential.* Gluing the stack periodically with an odd number of layers
   creates interlayer bonds inside one parity class and destroys the argument. Check
   `odd periodic layer ring breaks evenness` shows a three-layer periodic ring with
   \([q^3]Z=8\neq0\). So Theorem 1 is sharp in this hypothesis: free stacking of any height and
   periodic stacking of even height are covered, periodic odd height is not.

**Theorem 1′ (thermodynamic limit).** Let \(\Lambda_k\) be any sequence of stacks satisfying the
hypotheses of Theorem 1 and suppose \(f_{\Lambda_k}(K,q)\to f(K,q)\) pointwise for \(q\) in a
symmetric neighbourhood \((-\delta,\delta)\). Then \(f(K,q)=f(K,-q)\) on that neighbourhood. If in
addition \(f(K,\cdot)\) is \(n\) times differentiable at \(0\), then all odd derivatives up to order
\(n\) vanish; if \(f(K,\cdot)\) is analytic near \(0\), then \(c_{2m+1}=0\) for every \(m\ge0\).

*Proof.* Evenness passes to pointwise limits: \(f(K,q)=\lim f_{\Lambda_k}(K,q)=\lim
f_{\Lambda_k}(K,-q)=f(K,-q)\). An even function that is \(n\) times differentiable at \(0\) has
vanishing odd derivatives there, because \(f^{(j)}(0)=(-1)^jf^{(j)}(0)\). Analyticity gives the
statement for every order. \(\square\)

The regularity hypothesis in Theorem 1′ is stated, not proved. Evenness itself is unconditional;
only the passage to individual Taylor coefficients needs regularity of the limit. This is exactly
the honest scope: Theorem 1 is a theorem about every finite volume, Theorem 1′ is a theorem about
whatever limit one is entitled to take.

## 3. Theorem 2 — general even-coefficient formula

For random variables \(Y_1,\dots,Y_n\) the joint cumulant is
\[
 \kappa_0(Y_1,\dots,Y_n)
 =\sum_{\pi\in\Pi_n}(-1)^{|\pi|-1}(|\pi|-1)!\prod_{B\in\pi}\Big\langle\prod_{i\in B}Y_i\Big\rangle_0 .
\]

**Theorem 2.** With \(V=\sum_{g,x}Y_{g,x}\) as in Section 1,
\[
 f_\Lambda(K,q)-f_\Lambda(K,0)
 =\frac1{|\Lambda|}\log\langle e^{qV}\rangle_0
 =\sum_{n\ge1}\frac{q^n}{n!\,|\Lambda|}
 \sum_{b_1,\dots,b_n}\kappa_0\!\left(Y_{b_1},\dots,Y_{b_n}\right),
\]
the inner sum running over all ordered \(n\)-tuples of interlayer bonds. Moreover each joint
cumulant factorizes through two-dimensional moments:
\[
 \kappa_0\!\left(Y_{g_1,x_1},\dots,Y_{g_n,x_n}\right)
 =\sum_{\pi\in\Pi_n}(-1)^{|\pi|-1}(|\pi|-1)!
 \prod_{B\in\pi}\ \prod_{l}\ M_l\!\left(S_{l,B}\right),
\]
where \(S_{l,B}\) is the *parity-reduced* set of in-plane positions that the bonds of \(B\) place on
layer \(l\): position \(x\) belongs to \(S_{l,B}\) iff the number of \(i\in B\) with
\(x_i=x\) and \(l\in\{g_i,g_i+1\}\) is odd. If \(|S_{l,B}|\) is odd then \(M_l(S_{l,B})=0\).

*Proof.* The cumulant generating function of \(V\) under \(\mu_0\) is
\(\log\langle e^{qV}\rangle_0=\sum_{n\ge1}q^n\kappa_n(V)/n!\), convergent for all \(q\) because
\(V\) is a bounded random variable on a finite space. Joint cumulants are multilinear, so
\(\kappa_n(V)=\sum_{b_1,\dots,b_n}\kappa_0(Y_{b_1},\dots,Y_{b_n})\). For the factorization, fix a
partition \(\pi\) and a block \(B\). Then \(\prod_{i\in B}Y_{g_i,x_i}
=\prod_{i\in B}\sigma_{x_i,g_i}\sigma_{x_i,g_i+1}\) is a product of spins; grouping by layer and
using \(\sigma^2=1\) gives \(\prod_l\prod_{x\in S_{l,B}}\sigma_{x,l}\). Distinct layers are
independent under \(\mu_0\), so the block moment is \(\prod_lM_l(S_{l,B})\). Finally the zero-field
layer measure is invariant under a global spin flip within that layer, so any odd product of layer
spins has zero mean. \(\square\)

**Bulk form.** In the thermodynamic limit with in-plane translation invariance, fix the pattern
multiset of gap indices and use translation invariance to pin one position:
\[
 c_{2m}=\sum_{P}\ \frac1{\prod_g m_g(P)!}
 \sum_{x_2,\dots,x_{2m}}
 \kappa_0\!\left(Y_{g_1,0},Y_{g_2,x_2},\dots,Y_{g_{2m},x_{2m}}\right),
\]
where \(P\) runs over the admissible normalized gap patterns of Theorem 3, \(g_1\le\dots\le g_{2m}\)
is the non-decreasing arrangement of \(P\), and \(m_g(P)\) is the multiplicity of gap \(g\). The
combinatorial factor is the number \((2m)!/\prod_gm_g!\) of ordered gap assignments with multiset
\(P\), divided by \((2m)!\). Interchanging the (infinite) position sum with the coefficient
extraction requires absolute convergence of that sum; every coefficient computed in the artifact is
instead obtained from exact finite rectangles followed by finite-lattice Möbius inversion with a
bound-stability check, so no convergence hypothesis enters the computed numbers.

**Algorithm.** Theorem 2 with Theorem 3 is directly executable and is implemented in
`cumulant_polynomial` and `pattern_rectangle_sum`:

1. enumerate the \(2^{m-1}\) admissible gap patterns as compositions of \(m\);
2. for each pattern enumerate the set partitions of the \(2m\) bond indices;
3. discard a partition as soon as one block has odd bond multiplicity in some gap (Lemma 3.1);
4. factor each surviving block into one two-dimensional moment per layer, parity-reducing repeated
   in-plane positions;
5. weight by \((-1)^{|\pi|-1}(|\pi|-1)!\) and combine equal moment monomials over \(\mathbb Z\);
6. sum over in-plane position tuples and finite-lattice invert.

Verification: `generic generator versus literal definition m=1,2,3` compares the generated
polynomial against the literal partition sum on a non-Gaussian rational moment oracle for all seven
patterns with \(m\le3\); `generic generator reproduces c2 and c4 formulas` recovers
\(\kappa_2=G_{ij}^2\) and both fourth-order formulas that were established independently in
`proofs`/`results` for \(c_4\).

## 4. Theorem 3 — which layer-index patterns contribute

Call \(g_1,\dots,g_n\) the gap indices of a bond tuple, \(m_g\) their multiplicities, and the
*support* the set of occupied gaps.

**Lemma 3.1 (block parity).** Let \(\pi\) be a partition contributing a nonzero term to the
factorized cumulant of Theorem 2. Then every block \(B\in\pi\) contains an even number of bonds in
*every* gap.

*Proof.* Let \(g^\ast=\max_ig_i\) and let \(c_B(g)\) be the number of bonds of \(B\) in gap \(g\).
Layer \(g^\ast+1\) is met only by gap-\(g^\ast\) bonds, and parity reduction removes spins in pairs,
so \(|S_{g^\ast+1,B}|\equiv c_B(g^\ast)\pmod 2\); a nonzero term forces \(c_B(g^\ast)\) even.
Descending, layer \(l\) is met exactly by the bonds in gaps \(l-1\) and \(l\), so
\(|S_{l,B}|\equiv c_B(l-1)+c_B(l)\pmod2\). Given \(c_B(l)\) even by induction, a nonzero term forces
\(c_B(l-1)\) even. \(\square\)

**Lemma 3.2 (gap parity).** If some \(m_g\) is odd then the cumulant vanishes identically.

*Proof.* \(m_g=\sum_{B\in\pi}c_B(g)\). By Lemma 3.1 every surviving partition has all \(c_B(g)\)
even, so \(m_g\) is even. Contrapositively, if some \(m_g\) is odd then no partition survives and
every term of the cumulant is zero. \(\square\)

**Lemma 3.3 (interval support).** If the support is not an interval, the cumulant vanishes.

*Proof.* Let \(g^\ast\) be an unoccupied gap with occupied gaps on both sides. Bonds in gaps
\(<g^\ast\) involve only layers \(\le g^\ast\); bonds in gaps \(>g^\ast\) involve only layers
\(\ge g^\ast+1\). Under \(\mu_0\) distinct layers are independent, so these two nonempty families of
random variables are independent. A joint cumulant of a collection that splits into two independent
sub-collections, both nonempty, is zero. \(\square\)

The vanishing is *formal*: a block moment always factorizes into one moment per layer, so in the
polynomial ring generated by the symbols \(M_l(S)\) the mixed cumulant is identically zero, and the
generator returns the empty polynomial. Check `inadmissible patterns vanish formally` verifies this
for \((0,0,2,2)\), \((0,2,2,4)\), \((0,0,0,0,2,2)\), \((0,0,1,1,3,3)\) and for the parity failures
\((0,0,0,1)\), \((0,1,1,2)\), \((0,0,0,1,1,1)\); the standalone test repeats it independently.

**Theorem 3 (selection rule and pattern count).** Fix \(n=2m\). Modulo vertical translation, every
gap pattern with a not-identically-zero cumulant has the form: occupied gaps
\(0,1,\dots,k-1\) with multiplicities \(2a_1,\dots,2a_k\), where \((a_1,\dots,a_k)\) is a
composition of \(m\) into \(k\ge1\) positive parts. Consequently

* at most \(\sum_{k\ge1}\binom{m-1}{k-1}=2^{m-1}\) patterns contribute, and no other pattern is
  excluded by Lemmas 3.1–3.3;
* a contributing pattern occupies \(k+1\le m+1\) consecutive layers, with \(m+1\) attained only by
  \(a=(1,1,\dots,1)\);
* odd \(n\) admits no pattern at all, which re-proves \(c_{\mathrm{odd}}=0\) inside the cumulant
  framework, independently of Theorem 1.

*Proof.* By Lemma 3.2 every occupied gap multiplicity is even, so \(m_g=2a_g\) with \(a_g\ge1\); by
Lemma 3.3 the occupied gaps form an interval, normalized to \(\{0,\dots,k-1\}\); and
\(\sum_ga_g=m\). This is precisely a composition of \(m\). No such multiset is excluded by the two
lemmas, since both of their hypotheses fail for it. Counting
compositions of \(m\) into \(k\) parts gives \(\binom{m-1}{k-1}\), and \(\sum_k\binom{m-1}{k-1}
=2^{m-1}\). A pattern with \(k\) occupied gaps touches layers \(0,\dots,k\), i.e. \(k+1\) layers,
and \(k\le m\) because each part is \(\ge1\). For odd \(n\), \(\sum_g m_g=n\) with every \(m_g\)
even is impossible. \(\square\)

Verification. `layer-gap selection rule m=1..4` checks the composition construction against
exhaustive search over normalized gap multisets for \(m\le4\), the count \(2^{m-1}\) and the extent
\(m+1\) through \(m=7\); `block parity lemma m=1,2,3` checks the non-trivial equivalence of Lemma
3.1 ("even at every layer, blockwise") with "even in every gap, blockwise" over all \(203\)
partitions of every pattern with \(m\le3\); the exhaustive pattern search and the inadmissibility of
patterns such as \((0,0,0,1)\), \((0,1,1,2)\), \((0,0,2,2)\) are re-verified independently in
`tests/test_interlayer_allorders.py`.

The counts of *surviving partitions* per pattern at \(m\le3\) are recorded in the artifact
(`surviving_partition_counts_m1_m2_m3`): \(1\) for \((0,0)\); \(4\) and \(2\) at \(m=2\); and
\(31,11,11,5\) at \(m=3\). For a single-gap pattern the surviving partitions are exactly the
partitions of \([2m]\) into even blocks, so this count is \(1,4,31,379,\dots\)

**Theorem 3′ (reflection pairing).** The layer reflection \(l\mapsto L-1-l\) preserves \(\mu_0\) and
maps the pattern of the composition \((a_1,\dots,a_k)\) to the pattern of \((a_k,\dots,a_1)\).
Hence the two patterns contribute identically to \(c_{2m}\); a palindromic composition is
self-paired.

*Proof.* Reflection is a graph automorphism of the stack that fixes every layer's in-plane structure
and permutes the layers, so it maps \(\mu_0\) to itself and sends gap \(g\) to gap \(L-2-g\).
Relabelling the sum over positions accordingly identifies the two pattern contributions. \(\square\)

Verification: `layer-reflection pairing of gap patterns` confirms this on the computed series for
\(m\le3\); at \(m=3\) the contributions of \((0,0,0,0,1,1)\) and \((0,0,1,1,1,1)\) coincide
coefficient by coefficient, both equal to \(-\frac43v^2-\frac{166}3v^4-\frac{3092}3v^6\).

## 5. Theorem 4 — the adjacent-gap correction, formalized

At \(m=2\) the compositions of \(2\) are \((2)\) and \((1,1)\), giving the patterns \((0,0,0,0)\)
and \((0,0,1,1)\) with weights \(1/4!\) and \(1/(2!\,2!)=1/4\). Applying Theorem 2 to the two
patterns with distinct in-plane positions \(i,j,k,l\) and equal layer states:

**Theorem 4.** With \(G_{ij}=M(\{i,j\})\) and \(M_{ijkl}=M(\{i,j,k,l\})\),
\[
 \kappa_0\big(Y_{0,i},Y_{0,j},Y_{0,k},Y_{0,l}\big)
 =M_{ijkl}^2-G_{ij}^2G_{kl}^2-G_{ik}^2G_{jl}^2-G_{il}^2G_{jk}^2,
\]
\[
 \kappa_0\big(Y_{0,i},Y_{0,j},Y_{1,k},Y_{1,l}\big)
 =G_{ij}G_{kl}M_{ijkl}-G_{ij}^2G_{kl}^2,
\]
and
\[
 c_4=\frac1{4!}\sum_{j,k,l}\kappa_{\mathrm{same}}(0,j,k,l)
 +\frac1{4}\sum_{j,k,l}\kappa_{\mathrm{adj}}(0,j,k,l).
\]
More generally, at order \(2m\) the multi-gap patterns are the \(2^{m-1}-1\) compositions with
\(k\ge2\) parts, each entering with weight \(1/\prod_i(2a_i)!\), and the composition
\((1,\dots,1)\) needs exactly \(m+1\) layers. The adjacent-gap term of \(c_4\) is therefore the
first member of a family forced by Theorem 3, not a finite-slab artefact: a two-layer slab can never
see it, and no truncation to \(\le m\) layers can produce \(c_{2m}\) correctly.

*Proof.* Both displayed formulas are the specializations of Theorem 2. For the same-gap pattern,
each block must have even size (Lemma 3.1 with one gap), so \(\pi\) runs over the four even-block
partitions of \([4]\): the single block gives \(M_{ijkl}\) on both layers, i.e. \(M_{ijkl}^2\), with
Möbius coefficient \(+1\); each of the three pairings gives \(G\cdot G\) on both layers with Möbius
coefficient \(-1\). For the adjacent pattern the three layers are \(0,1,2\); the bonds \(i,j\) meet
layers \(0,1\) and the bonds \(k,l\) meet layers \(1,2\). Blocks must be even in each gap, so only
two partitions survive: the single block, giving \(G_{ij}\) on layer \(0\), \(M_{ijkl}\) on layer
\(1\) and \(G_{kl}\) on layer \(2\); and the split \(\{i,j\}|\{k,l\}\), giving
\(G_{ij}\cdot G_{ij}G_{kl}\cdot G_{kl}\) with Möbius coefficient \(-1\). Because layers carry equal
states this is \(G_{ij}G_{kl}M_{ijkl}-G_{ij}^2G_{kl}^2\). The weights come from Theorem 2's bulk
form. The general statement is Theorem 3 plus the weight computation. \(\square\)

Verification of the extent statement at \(m=3\): `vertical extent four isolates the maximal c6
pattern` shows the independent even-subgraph four-layer finite-lattice weight at \(w^6\) equals the
\((0,0,1,1,2,2)\) cumulant contribution exactly, and `vertical extent bound m+1 verified to m=6`
shows all weights of slabs taller than \(m+1\) layers vanish at \(w^{2m}\) for \(m\le6\).

## 6. Exact coefficients, including a fresh \(c_6\)

Everything below is an exact rational finite-lattice-inverted prefix for the isotropic simple cubic
stack, with a bound-stability check (adding every rectangle with one extra in-plane span unit changes
nothing). Direct coupling \(q\):
\[
 c_2=\tfrac12+2v^2+18v^4+118v^6+778v^8,
\]
\[
 c_4=-\tfrac1{12}+\tfrac23v^2+51v^4+\tfrac{2914}3v^6+\tfrac{41113}3v^8,
\]
\[
 c_6=\tfrac1{45}+\tfrac4{45}v^2+\tfrac{304}5v^4+\tfrac{144236}{45}v^6 .
\]
In the repository variable \(w=\tanh q\), using the exact change of variable
\(q=\operatorname{atanh}w\):
\[
 [w^2]=c_2,\qquad [w^4]=c_4+\tfrac23c_2,\qquad [w^6]=c_6+\tfrac43c_4+\tfrac{23}{45}c_2,
\]
\[
 [w^2]f=\tfrac12+2v^2+18v^4+118v^6+778v^8,\quad
 [w^4]f=\tfrac14+2v^2+63v^4+1050v^6+14223v^8,
\]
\[
 [w^6]f=\tfrac16+2v^2+138v^4+\tfrac{13682}3v^6 .
\]
The \(c_2\) and \(c_4\) rows reproduce the previously certified artifact
`results/interlayer/c4_series.json` exactly, including the same-gap component
\(-\frac1{12}-\frac43v^2-17v^4-\frac{464}3v^6-\frac{3527}3v^8\) and the adjacent-gap component
\(2v^2+68v^4+1126v^6+14880v^8\). The \(c_6\) row is new here and was produced twice independently
inside this task:

1. the generic cumulant generator of Theorem 2 over exact two-dimensional moments with rectangular
   finite-lattice inversion;
2. an independent even-subgraph route: exact \([w^d]P(v,w)\) columns of open \(w\times h\times L\)
   slabs built from layer boundary masks, an exact Newton \(\log\) recursion, and a three-dimensional
   finite-lattice inversion over \(L\le4\).

A third route, full spin enumeration of open boxes through the integer joint density of states
(`ising.interlayer.anisotropic_box_even_subgraph`), reproduces \([w^2],[w^4],[w^6]\) through \(v^4\)
in `tests/test_interlayer_allorders.py`.

A fourth, fully external confirmation was obtained without reading any other artifact: the
coefficients above were sent verbatim to the concurrent agent `InterlayerC6`, which had computed
\(c_6\) by its own route, and it reported exact agreement on every coefficient, on the \(w\)
conversion, and on the pattern split (its three-layer total equals the sum of the two mirror
patterns, consistent with Theorem 3′).

Two exact normalization controls hold: at \(v=0\) the coefficients are
\(\frac12,-\frac1{12},\frac1{45}\), which are \([q^{2m}]\log\cosh q\); and in the \(w\) variable they
are \(\frac12,\frac14,\frac16\), i.e. \([w^{2m}]\!\left(-\frac12\log(1-w^2)\right)\). The pattern
decomposition of \(c_6\) in \(q\) is
\[
 (0,0,0,0,0,0):\ \tfrac1{45}+\tfrac{34}{45}v^2+\tfrac{202}{15}v^4+\tfrac{7406}{45}v^6,
\]
\[
 (0,0,0,0,1,1)=(0,0,1,1,1,1):\ -\tfrac43v^2-\tfrac{166}3v^4-\tfrac{3092}3v^6,
\qquad
 (0,0,1,1,2,2):\ 2v^2+158v^4+5102v^6 .
\]
All four patterns are nonzero, so the selection rule of Theorem 3 is not merely an upper bound on
the contributing set at \(m=3\).

## 7. Theorem 5 — an all-order coefficient identity at second in-plane order

**Theorem 5.** For the simple cubic stack with in-plane variables \(v_x,v_y\) and interlayer
variable \(w\), and \(\Phi=f-f_{2\mathrm D}-\log\cosh K_z\) the residual pressure per site,
\[
 [v_x^2w^{2m}]\,\Phi=[v_y^2w^{2m}]\,\Phi=1\quad\text{for every }m\ge1,
 \qquad\text{hence}\qquad [v^2w^{2m}]\,\Phi=2 .
\]
Moreover the entire contribution comes from clusters of exactly \(m+1\) layers, saturating the extent
bound of Theorem 3.

*Proof.* Use the standard even-subgraph representation
\(Z=2^{|\Lambda|}\prod_e\cosh K_e\cdot P\), \(P=\sum_{S\ \mathrm{even}}\prod_{e\in S}t_e\) with
\(t_e\in\{v_x,v_y,w\}\), so that \(\Phi=\lim|\Lambda|^{-1}\log P-\lim|\Lambda|^{-1}\log P_{2\mathrm D}\)
and the \(w\)-dependent part of \(\Phi\) is \(\lim|\Lambda|^{-1}\log(P/P|_{w=0})\).

*Step 1: every nonempty even subgraph uses at least two in-plane edges.* A nonempty even subgraph
contains a cycle. Along a cycle the in-plane position changes only at in-plane edges; vertical edges
at a fixed in-plane position span a path graph (the vertical line), which is acyclic, so the cycle
uses at least one in-plane edge, and with exactly one it could not return to its starting in-plane
position. Hence at least two.

*Step 2: the even subgraphs with exactly two in-plane edges are precisely the ladder cycles.* Let
\(S\) be even with exactly two in-plane edges. By Step 1, \(S\) has exactly one nonempty connected
component. Every vertex of \(S\) has even degree, at most \(2\) from vertical edges (subgraph of a
path) and at most \(2\) from in-plane edges (only two exist). A vertex of degree \(4\) would carry
both in-plane edges, say \(\{(x,l),(y,l)\}\) and \(\{(x,l),(y',l)\}\); then \((y,l)\) has in-plane
degree \(1\), so its vertical degree is odd, i.e. \(1\); following that vertical run to its far end
gives a vertex with vertical degree \(1\) and no in-plane edge left, so odd total degree — a
contradiction. Hence all degrees are \(2\) and \(S\) is a single cycle. Such a cycle consists of an
in-plane edge \(\{(x,l),(y,l)\}\), a vertical run at \(y\) from \(l\) to \(l'\), the in-plane edge
\(\{(x,l'),(y,l')\}\) over the *same* in-plane pair, and the vertical run back at \(x\). Its weight
is \(t_{xy}^2w^{2|l-l'|}\).

*Step 3: count and logarithm.* Writing \(P=1+X\), the terms of \(X^j\) with \(j\ge2\) carry at least
\(4\) in-plane edges by Step 1, so
\([t_{xy}^2w^{2m}]\log P=[t_{xy}^2w^{2m}]X\). By Step 2 the cycles contributing \(v_x^2w^{2m}\) are
indexed by an in-plane \(x\)-direction edge and the lower layer, one per site of the stack; hence the
per-site count is \(1\), and likewise for \(v_y\). Each such cycle occupies layers \(l,\dots,l+m\),
i.e. exactly \(m+1\) layers. The two-dimensional subtraction removes no \(w\)-dependent term.
\(\square\)

Verification: `all-order [v^2 w^(2m)] ladder identity` gives exactly \(2\) for \(m=1,\dots,6\) from
the independent even-subgraph finite-lattice route, `all-order identity bound stability` confirms the
values are unchanged when the in-plane box family is enlarged, `vertical extent bound m+1 verified to
m=6` confirms the concentration on \(m+1\) layers, and the standalone test re-derives \(m\le3\) from
full spin enumeration. The identity is elementary, but it is a genuine all-\(m\) statement about the
three-dimensional series, and it is consistent with the certified rows above:
\([v^2]\) equals \(2\) in \([w^2]f\), \([w^4]f\) and \([w^6]f\).

## 8. What is verified where

| Statement | Status | Evidence |
| --- | --- | --- |
| \(Z_\Lambda(q)=Z_\Lambda(-q)\), all odd coefficients zero, finite volume | [THEOREM] | Section 2; checks `finite slab interlayer evenness`, `odd vertical columns vanish` |
| Thermodynamic-limit evenness; odd coefficients zero under stated regularity | [THEOREM] (conditional statement, hypothesis explicit) | Theorem 1′ |
| Layer-parity hypothesis is necessary | [THEOREM] + exact witness | check `odd periodic layer ring breaks evenness`, \([q^3]Z=8\) |
| Set-partition formula for \(c_n\), factorized through 2D moments | [THEOREM] | Section 3; checks `generic generator versus literal definition m=1,2,3` |
| Selection rule; \(2^{m-1}\) patterns; extent \(m+1\) | [THEOREM] | Section 4; checks `layer-gap selection rule m=1..4`, `block parity lemma m=1,2,3`, `inadmissible patterns vanish formally`, `vertical extent bound m+1 verified to m=6` |
| Reflection pairing of patterns | [THEOREM] | Theorem 3′; check `layer-reflection pairing of gap patterns` |
| \(c_4\) same-gap and adjacent-gap formulas; general multi-gap weights | [THEOREM] | Section 5; check `generic generator reproduces c2 and c4 formulas` |
| \(c_2,c_4\) exact prefixes through \(v^8\) | [COMPUTATION] | checks `stored c2 prefix reproduced`, `stored c4 prefix and gap split reproduced`, `exact q to w change of variable` |
| \(c_6\) exact prefix through \(v^6\), three routes | [COMPUTATION] | checks `two-route agreement through w^6`, `c6 decoupled normalization`, `finite-lattice bound stability`; third route in the standalone test |
| \([v^2w^{2m}]\Phi=2\) for all \(m\ge1\) | [THEOREM] | Section 7; checks `all-order [v^2 w^(2m)] ladder identity`, `all-order identity bound stability` |

## 9. Non-results

* **No sign or positivity theorem.** The Möbius weights alternate; the computed \(c_4\) same-gap
  component is negative and the adjacent-gap component positive, and \(c_6\) pattern components have
  both signs. Nothing here proves a definite sign for \(c_{2m}\), for its pattern components, or for
  the coefficients of the \(v\) series. [UNRESOLVED]
* **No radius of convergence.** No bound is proved on the interlayer series in \(q\), in \(w\), or on
  the in-plane series in \(v\); in particular nothing here is valid at or near the isotropic critical
  point. The Theorem 2 bulk form carries an explicit absolute-convergence hypothesis, and the
  computed coefficients are finite-order and finite-lattice certified rather than resummed.
  [UNRESOLVED]
* **No closed form.** Theorem 2 reduces \(c_{2m}\) to two-dimensional multipoint moments, which are
  not evaluated in closed form here. The pattern count grows like \(2^{m-1}\) and the number of
  surviving partitions grows superexponentially (\(1,4,31,379,\dots\) already for the single-gap
  family), so the result is structural.
* **The 3D Ising model is not solved.**

## 10. Reproduction

```
.venv/bin/python experiments/e75_interlayer_allorders.py
.venv/bin/python tests/test_interlayer_allorders.py
```

Both print `PASS` only if every exact check passes; the experiment refuses to write the artifact
otherwise. All arithmetic is `int`/`fractions.Fraction`; no floating point is used anywhere in this
front, so no precision or error bound has to be recorded.
