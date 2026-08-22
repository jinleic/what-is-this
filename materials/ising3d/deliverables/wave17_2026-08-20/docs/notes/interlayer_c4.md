# Fourth-order interlayer coefficient from exact layer cumulants

## Variables, normalization, and the necessary change of variable

[LEMMA] Let a square layer contain `A` spins, let the stack contain `M > 4` layers before the thermodynamic limit, and divide `log Z` by `AM`; the restriction `M > 4` excludes fixed-height order-four wrapping graphs while leaving every local order-four cluster unchanged.

[LEMMA] Write
\[
q=K_z,\qquad w=\tanh K_z,\qquad v=\tanh K,
\]
and define `q` as the direct Boltzmann coupling and `w` as the high-temperature interlayer variable used by the wave-5/6 artifacts.

[LEMMA] The fourth-cumulant formula applies first to `q`, because with
\[
X_{\ell i}=\sigma_{\ell i}\sigma_{\ell+1,i},\qquad
H_{\rm int}=\sum_{\ell,i}X_{\ell i},
\]
the decoupled-layer expectation gives
\[
\frac1{AM}\log\left\langle e^{qH_{\rm int}}\right\rangle_0
 =\frac{q^2}{2AM}\kappa_2(H_{\rm int})
  +\frac{q^4}{24AM}\kappa_4(H_{\rm int})+O(q^6).
\]

[LEMMA] Since `q = atanh(w) = w + w^3/3 + O(w^5)`, coefficients in the two conventions obey
\[
c_{2,w}^{\rm tot}=c_{2,q},\qquad
c_{4,w}^{\rm tot}=c_{4,q}+\frac23c_{2,q}.
\]

[LEMMA] Separating the vertical-bond prefactor uses
\[
\log\cosh(\operatorname{atanh}w)
=-\frac12\log(1-w^2)=\frac12w^2+\frac14w^4+O(w^6),
\]
so
\[
c_{2,w}^{\rm res}=c_{2,q}-\frac12,
\qquad
c_{4,w}^{\rm res}=c_{4,q}+\frac23c_{2,q}-\frac14.
\]

## Symbolic fourth-cumulant factorization

[LEMMA] For four zero-mean random variables the complete set-partition formula is
\[
\kappa_4(Y_1,Y_2,Y_3,Y_4)
=\langle Y_1Y_2Y_3Y_4\rangle
-\langle Y_1Y_2\rangle\langle Y_3Y_4\rangle
-\langle Y_1Y_3\rangle\langle Y_2Y_4\rangle
-\langle Y_1Y_4\rangle\langle Y_2Y_3\rangle.
\]

[LEMMA] Put
\[
G_{ij}=\langle\sigma_i\sigma_j\rangle_{2D},
\qquad
M_{ijkl}=\langle\sigma_i\sigma_j\sigma_k\sigma_l\rangle_{2D}.
\]
For four vertical bonds on one layer gap, independence of the two layers gives every term separately as
\[
\begin{aligned}
\langle X_iX_jX_kX_l\rangle_0&=M_{ijkl}^2,\\
\langle X_iX_j\rangle_0\langle X_kX_l\rangle_0&=G_{ij}^2G_{kl}^2,\\
\langle X_iX_k\rangle_0\langle X_jX_l\rangle_0&=G_{ik}^2G_{jl}^2,\\
\langle X_iX_l\rangle_0\langle X_jX_k\rangle_0&=G_{il}^2G_{jk}^2.
\end{aligned}
\]

[LEMMA] Therefore the requested two-independent-layer identity is
\[
\boxed{\kappa^{(0)}_{ijkl}
=M_{ijkl}^2-G_{ij}^2G_{kl}^2-G_{ik}^2G_{jl}^2-G_{il}^2G_{jk}^2.}
\]

[LEMMA] A full three-dimensional stack has one additional connected order-four geometry: two bonds `i,j` on gap `(ell,ell+1)` and two bonds `k,l` on the adjacent gap `(ell+1,ell+2)`.

[LEMMA] For that adjacent-gap geometry, layer factorization gives every term separately as
\[
\begin{aligned}
\langle X_{\ell i}X_{\ell j}X_{\ell+1,k}X_{\ell+1,l}\rangle_0
 &=G_{ij}M_{ijkl}G_{kl},\\
\langle X_{\ell i}X_{\ell j}\rangle_0
 \langle X_{\ell+1,k}X_{\ell+1,l}\rangle_0
 &=G_{ij}^2G_{kl}^2,\\
\langle X_{\ell i}X_{\ell+1,k}\rangle_0
 \langle X_{\ell j}X_{\ell+1,l}\rangle_0&=0,\\
\langle X_{\ell i}X_{\ell+1,l}\rangle_0
 \langle X_{\ell j}X_{\ell+1,k}\rangle_0&=0.
\end{aligned}
\]

[LEMMA] The two cross-gap pair moments vanish because each leaves an outer layer with one spin insertion, and hence
\[
\boxed{\kappa^{(1)}_{ijkl}
=G_{ij}G_{kl}M_{ijkl}-G_{ij}^2G_{kl}^2
=G_{ij}G_{kl}(M_{ijkl}-G_{ij}G_{kl}).}
\]

[LEMMA] The connected 2D four-spin function used by the linked-cluster interpretation is
\[
U_{ijkl}
=M_{ijkl}-G_{ij}G_{kl}-G_{ik}G_{jl}-G_{il}G_{jk}.
\]
Writing
\[
P=G_{ij}G_{kl},\qquad Q=G_{ik}G_{jl},\qquad R=G_{il}G_{jk},
\]
the two exact cumulant inputs become
\[
\boxed{\kappa^{(0)}=U^2+2U(P+Q+R)+2(PQ+PR+QR),\qquad
\kappa^{(1)}=P(U+Q+R).}
\]

[COMPUTATION] The experiment obtains `G`, `M`, and therefore `U=M-P-Q-R` as exact `Fraction` series from the same edge-boundary polynomials; it sums the algebraically equivalent `M` formulas to avoid constructing a second convention for the connected function.

[LEMMA] The gap-count patterns are exhaustive: four bonds on one gap give `kappa^(0)`; a `2+2` split on adjacent gaps gives `kappa^(1)`; a `2+2` split on disjoint gaps has zero joint cumulant by independence; and every `3+1`, `2+1+1`, or `1+1+1+1` nonwrapping pattern leaves an outer layer with odd spin parity.

[LEMMA] There are six placements of two ordered entries on each of two adjacent gaps, so after summing over the `M` gaps and the `M` adjacent gap pairs, dividing by `AM`, and using translation invariance to anchor `i=0`, the full-stack direct-coupling coefficient is
\[
\boxed{
c_{4,q}(K)=
\frac1{24}\sum_{j,k,l}\kappa^{(0)}_{0jkl}
+\frac14\sum_{j,k,l}\kappa^{(1)}_{0jkl}.}
\]

[LEMMA] A two-layer slab contains only the first boxed term, so an `L x L x 2` calculation alone is not a valid check of the full three-dimensional `c4`; a height-three contribution or a vertical finite-lattice inversion is required.

## Exact two-dimensional inputs and finite-lattice reduction

[LEMMA] On an open square rectangle, let `P_S(v)` count in-plane edge subsets whose odd-degree boundary is the vertex mask `S`; the exact high-temperature identities are
\[
P_S(v)=\sum_{E:\,\partial E=S}v^{|E|},\qquad
\left\langle\prod_{i\in S}\sigma_i\right\rangle=\frac{P_S(v)}{P_\varnothing(v)}.
\]

[COMPUTATION] `e51_interlayer_c4.py` enumerates every `P_S` with exact integers and performs every ratio, product, cumulant subtraction, and finite-lattice inversion with `Fraction` arithmetic.

[LEMMA] For one finite rectangle define
\[
S_2=\sum_{i,j}G_{ij}^2,\qquad
S_4=\sum_{i,j,k,l}M_{ijkl}^2,\qquad
T=\sum_{i,j,k,l}G_{ij}G_{kl}M_{ijkl}.
\]
The fully ordered sums reduce the three same-gap pair partitions to three identical copies of `S_2^2`, and therefore the finite-rectangle quantities inverted by the first route are
\[
c_{2,q}^{\rm box}=\frac12S_2,\qquad
c_{4,q}^{(0),\rm box}=\frac1{24}(S_4-3S_2^2),\qquad
c_{4,q}^{(1),\rm box}=\frac14(T-S_2^2).
\]

[LEMMA] Rectangular Möbius inversion removes every contribution embeddable in a proper subrectangle; after this linked-cluster subtraction, a contributing projected even multigraph spanning an `a x b` rectangle crosses each of its `(a-1)+(b-1)` coordinate cuts at least twice and thus uses at least `2((a-1)+(b-1))` in-plane edges.

[LEMMA] Four vertical edges cross each horizontal layer cut evenly and can span at most two adjacent layer gaps, so height at most three and in-plane span at most six are complete through `v^12`.

[COMPUTATION] The cumulant route used all 28 oriented open rectangles with `(a-1)+(b-1) <= 6`, and the independent anisotropic route used their 84 height-one, height-two, and height-three slabs.

## Independent anisotropic construction

[LEMMA] The independent route constructs the full slab even-subgraph polynomial without inserting the cumulant formula: for two layers its `w^d` column is
\[
P_d^{(2)}(v)=\sum_{|S|=d}P_S(v)^2,
\]
and for three layers it is
\[
P_d^{(3)}(v)=
\sum_{|S|+|T|=d}P_S(v)P_{S\triangle T}(v)P_T(v).
\]

[LEMMA] Because a graph boundary has even cardinality, every odd-cardinality `P_S` vanishes, so the constructed `w^1` and `w^3` columns and hence `c3` vanish identically.

[LEMMA] With odd columns zero, each finite slab supplies
\[
[w^2]\log P=\frac{P_2}{P_0},\qquad
[w^4]\log P=\frac{P_4}{P_0}-\frac12\left(\frac{P_2}{P_0}\right)^2,
\]
which is then inverted over all three box dimensions.

[COMPUTATION] The two-layer slab weights agree exactly with the same-gap cumulant after the `q`-to-`w` conversion, and the height-three spanning weights agree exactly with the adjacent-gap cumulant at every coefficient through `v^12`.

## Exact coefficients

[COMPUTATION] The reproduced and extended second-order coefficient is
\[
c_{2,q}(v)=c_{2,w}^{\rm tot}(v)
=\frac12+2v^2+18v^4+118v^6+778v^8+4978v^{10}+31398v^{12}+O(v^{14}).
\]

[COMPUTATION] The same-gap, bilayer part of the direct-coupling fourth cumulant is
\[
c_{4,q}^{(0)}(v)=
-\frac1{12}-\frac43v^2-17v^4-\frac{464}{3}v^6
-\frac{3527}{3}v^8-\frac{23024}{3}v^{10}-41539v^{12}+O(v^{14}).
\]

[COMPUTATION] The adjacent-gap, genuinely three-layer part is
\[
c_{4,q}^{(1)}(v)=
2v^2+68v^4+1126v^6+14880v^8+168306v^{10}+1721500v^{12}+O(v^{14}).
\]

[COMPUTATION] The full-stack direct-coupling coefficient is therefore
\[
\boxed{
c_{4,q}(v)=
-\frac1{12}+\frac23v^2+51v^4+\frac{2914}{3}v^6
+\frac{41113}{3}v^8+\frac{481894}{3}v^{10}+1679961v^{12}
+O(v^{14}).}
\]

[COMPUTATION] In the repository wave convention `w=tanh K_z`, the requested full-stack coefficient including the single-bond prefactor is
\[
\boxed{
c_{4,w}^{\rm tot}(v)=
\frac14+2v^2+63v^4+1050v^6+14223v^8+163950v^{10}+1700893v^{12}
+O(v^{14}).}
\]

[COMPUTATION] After separating `log cosh(K_z)`, the wave-convention residual is
\[
\boxed{
c_{4,w}^{\rm res}(v)=
2v^2+63v^4+1050v^6+14223v^8+163950v^{10}+1700893v^{12}
+O(v^{14}).}
\]

[COMPUTATION] For comparison, a bilayer-only high-temperature calculation gives
\[
c_{4,w}^{\rm bilayer,res}(v)
=-5v^4-76v^6-657v^8-4356v^{10}-20607v^{12}+O(v^{14}),
\]
which differs from the full stack precisely by the nonzero height-three contribution.

## Verification and scope

[COMPUTATION] Route 1, the exact 2D moment/cumulant finite-lattice calculation, and route 2, the direct anisotropic slab even-subgraph finite-lattice calculation, agree coefficient by coefficient through `v^12` for `c2`, both fourth-order components, and full `c4`.

[COMPUTATION] Adding every rectangle with one extra in-plane span unit changes neither route through `v^8`, and both routes reproduce the wave-6 value `[v^8]c2=778` exactly.

[COMPUTATION] `tests/test_interlayer_c4.py` independently recomputes the full-stack result through `v^6` from exact spin density-of-states transforms of open three-dimensional boxes and obtains `c4_w^res = 2v^2+63v^4+1050v^6` and `c4_q = -1/12+(2/3)v^2+51v^4+(2914/3)v^6` exactly.

[COMPUTATION] At `K=0`, the arrays reduce to `c2_q=1/2`, `c4_q=-1/12`, `c4_w^tot=1/4`, and `c4_w^res=0`, exactly matching `log cosh K_z`.

[COMPUTATION] The machine-readable exact arrays, all order-by-order route comparisons, box lists, span bounds, and passing checks are stored in `results/interlayer/c4_series.json`.

[UNRESOLVED] This finite high-temperature prefix does not determine the convergence radius in `K_z`, does not justify evaluation at the isotropic critical point, and does not constitute an exact solution of the three-dimensional Ising model.
