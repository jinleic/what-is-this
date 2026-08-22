# Interlayer-coupling expansion about exactly decoupled 2D layers

## Conventions and finite-volume identity

[THEOREM] Let one square layer have `A` sites, let there be `M >= 3` layers with periodic boundary conditions in the layer direction, and put `N = AM`; the restriction `M >= 3` avoids the repository's length-two parallel-bond degeneracy while the thermodynamic limit is taken.

[THEOREM] With normalized expectation `⟨·⟩_0` in the product of `M` independent zero-field 2D layers, `X_(ell,r) = sigma_(ell,r) sigma_(ell+1,r)`, and `w = tanh K_z`, the exact finite-volume identity is

\[
 Z_{A,M}(K,K_z)
 = (\cosh K_z)^{AM}\, Z_{2,A}(K)^M\,Q_{A,M}(w),
 \qquad
 Q_{A,M}(w)=\sum_{S\subseteq B_z}w^{|S|}
 \left\langle\prod_{b\in S}X_b\right\rangle_0.
\]

[THEOREM] Consequently,

\[
 \phi_{A,M}(K,K_z)=\log\cosh K_z+\phi_{2,A}(K)
       +\frac1{AM}\log Q_{A,M}(w).
\]

[THEOREM] This note uses the convention forced by that displayed decomposition:

\[
 c_n^{\rm res}(K)=[w^n]\{\phi_{3D}-\phi_{2D}-\log\cosh K_z\}.
\]

[THEOREM] A separately named coefficient includes the single-bond prefactor:

\[
 c_n^{\rm tot}(K)=[w^n]\{\phi_{3D}-\phi_{2D}\}.
\]

## The coefficients at first and second order

[THEOREM] Independent spin reversal in each zero-field layer makes an expectation vanish whenever some layer contains an odd number of spin insertions.

[THEOREM] A one-bond term contains one spin in each of two layers, so

\[
 c_1^{\rm res}(K)=\frac1{AM}\sum_b\langle X_b\rangle_0
 =\langle\sigma\rangle_{2D}^2=0.
\]

[THEOREM] Since `log cosh(atanh w)` is even in `w`, the total convention also has `c_1^tot(K)=0`.

[LEMMA] Writing `Q = 1 + a_1 w + a_2 w^2 + ...`, spin reversal gives `a_1=0`, and therefore

\[
 [w^2]\log Q=a_2
 =\sum_{b<c}\langle X_bX_c\rangle_0.
\]

[LEMMA] Two distinct interlayer bonds on different layer pairs contribute zero: disjoint pairs leave four layers with one insertion each, and adjacent pairs leave each of the two outer layers with one insertion.

[LEMMA] The only nonzero second-order terms have both bonds between the same pair of adjacent layers and at distinct in-plane sites `r != r'`; layer independence then gives

\[
 \left\langle X_{(\ell,r)}X_{(\ell,r')}\right\rangle_0
 =\left\langle\sigma_r\sigma_{r'}\right\rangle_{2D}^2
 =G_A(r-r')^2.
\]

[THEOREM] The combinatorial convention is an unordered pair of distinct interlayer bonds: for each of `M` layer gaps, `A` choices of the first site and every nonzero displacement count each unordered pair twice, so

\[
 \sum_{b<c}\langle X_bX_c\rangle_0
 =\frac{MA}{2}\sum_{r\ne0}G_A(r)^2.
\]

[THEOREM] Dividing by `AM` and taking the translation-invariant thermodynamic limit, whenever the sum converges, gives the residual coefficient

\[
 \boxed{c_2^{\rm res}(K)=\frac12\sum_{r\ne0}G(r)^2.}
\]

[THEOREM] The omitted `r=0` term is exactly the coefficient already isolated in the bond prefactor because

\[
 \log\cosh K_z=-\frac12\log(1-w^2)
 =\frac{w^2}{2}+\frac{w^4}{4}+\frac{w^6}{6}+\cdots,
 \qquad G(0)=1.
\]

[THEOREM] Thus the requested all-site formula is correct for the **total** coefficient, not for `c_2` after `log cosh K_z` has been written separately:

\[
 \boxed{c_2^{\rm tot}(K)=\frac12\sum_{r}G(r)^2
 =\frac12+c_2^{\rm res}(K).}
\]

[THEOREM] At `K=0`, `G(r)=delta_(r,0)`, hence `c_2^res(0)=0`, `c_2^tot(0)=1/2`, and the full result is exactly `phi_3D(0,K_z)=log(2 cosh K_z)`; this control would fail if the `r=0` term were counted in both conventions.

[EXTERNAL] At fourth and higher order, the logarithm selects connected bond cumulants, and `c_4` includes connected 2D four-spin functions rather than only products of `G`; fixed 2D multipoint correlations are in principle available from the exact fermionic/Pfaffian solution, but `c_4` was not computed here.

## Independent exact series cross-check

[COMPUTATION] The two-dimensional route in `exact_2d_overlap_series(6)` used exact integer polynomial row transfer on a width-seven periodic cylinder of height thirteen, with the source six rows from either open boundary; it did not use finite-lattice Möbius inversion.

[LEMMA] Through degree six, a horizontal wrap in one member of a squared-correlation graph costs at least the circumference seven in the combined degree, a vacuum wrap also costs seven, and an excursion to an omitted open-boundary edge costs more than six; therefore this finite cylinder returns the infinite-square-lattice truncation exactly through `v^6`.

[COMPUTATION] The resulting exact correlation series is

\[
 \sum_rG(r)^2=1+4v^2+36v^4+236v^6+O(v^8).
\]

[COMPUTATION] The two coefficient conventions are consequently

\[
 c_2^{\rm tot}(v)=\frac12+2v^2+18v^4+118v^6+O(v^8),
 \qquad
 c_2^{\rm res}(v)=2v^2+18v^4+118v^6+O(v^8).
\]

[COMPUTATION] The independent three-dimensional route transformed the exact anisotropic integer density of states of every required open box into its bivariate even-subgraph polynomial `P(v,w)`, extracted `[w^2] log P`, and performed oriented rectangular finite-lattice Möbius inversion.

[LEMMA] A connected even graph with exactly two vertical edges crosses each horizontal layer cut evenly, so it has vertical extent exactly one; if its in-plane bounding rectangle is `a x b`, it uses at least `2((a-1)+(b-1))` in-plane edges.

[COMPUTATION] The ten boxes satisfying `(a-1)+(b-1) <= 3` therefore suffice through `v^6`, and adding five boxes with one unit of extra span leaves every coefficient unchanged.

[COMPUTATION] The anisotropic finite-lattice result is exactly

\[
 [w^2]\{\phi_{3D}-\phi_{2D}-\log\cosh K_z\}
 =2v^2+18v^4+118v^6+O(v^8),
\]

[COMPUTATION] The 2D-correlation and anisotropic-3D finite-lattice coefficient arrays agree as exact `Fraction` values at all degrees through six, including the first three nonzero orders `v^2`, `v^4`, and `v^6`.

## Numerical evaluation and the critical divergence

[COMPUTATION] The numerical calculation did **not** use the Wu/McCoy--Tracy--Barouch axial or diagonal Toeplitz determinants, because axial and diagonal values alone do not determine the required sum over every displacement.

[COMPUTATION] Instead it evaluated all finite-torus correlations with the algebraically exact symmetric row-transfer formula

\[
 T_{s,t}=\exp\!\left[K\left(\frac{E_h(s)}2+s\!\cdot\!t+\frac{E_h(t)}2\right)\right],
 \qquad
 G_L(x,y)=\frac{\operatorname{Tr}(S_0T^yS_xT^{L-y})}{\operatorname{Tr}T^L}.
\]

[COMPUTATION] Binary64 arithmetic was used for periodic tori `L=3,...,10`; an independent 70-decimal mpmath evaluation at `L=5` and exact `K_c^(2D)` differed by `3.55e-15`, while the largest observed `G(0)-1` or inversion-symmetry residual was `3.33e-16`.

[COMPUTATION] For subcritical couplings, a three-size `L=6,8,10` geometric-tail extrapolation gave the following **total** coefficient; each quoted uncertainty is twice the extrapolation correction plus roundoff and is an empirical finite-size estimate, not a rigorous interval.

| `K` | `c_2^tot(K)` | empirical absolute error | `c_2^res=c_2^tot-1/2` |
|---:|---:|---:|---:|
| 0.10 | 0.521767390652106 | 1.07e-8 | 0.021767390652106 |
| 0.15 | 0.554682092226023 | 1.50e-6 | 0.054682092226023 |
| 0.20 | 0.614589998224902 | 7.02e-5 | 0.114589998224902 |
| 0.25 | 0.726419529490255 | 1.86e-3 | 0.226419529490255 |
| 0.30 | 0.955057398940567 | 3.97e-2 | 0.455057398940567 |

[COMPUTATION] The growing uncertainty is itself a finite-size warning, so no near-critical infinite-volume number was inferred from the `L <= 10` tori.

[EXTERNAL] At the exact square-lattice critical point

\[
 K_c^{2D}=\frac12\log(1+\sqrt2),
 \qquad G(r)\sim A(\widehat r)|r|^{-\eta},\quad\eta=\frac14.
\]

[THEOREM] Positivity of the critical asymptotic implies

\[
 \sum_{|r|\le R}G(r)^2\asymp
 \sum_{n\le R}n\,n^{-1/2}\asymp R^{3/2},
\]

so the infinite-lattice sum and `c_2` are not finite at `K=K_c^(2D)`.

[EXTERNAL] On the disordered side, the exact 2D correlation-length exponent is `nu=1`; the scaling prediction is therefore

\[
 \sum_rG(r)^2\sim \xi^{2-2\eta}
 \sim (K_c^{2D}-K)^{-\nu(2-2\eta)}
 =(K_c^{2D}-K)^{-3/2}.
\]

[COMPUTATION] At criticality the torus data grow from `sum_r G_L(r)^2 = 5.94065333949` at `L=3` to `37.8028232436` at `L=10`; a log-log fit on `L=7,...,10` gives exponent `1.51263`, and varying to a `1/L`-corrected fit gives `1.48792`, recorded as `1.51263 +/- 0.02471`, consistent with the exact prediction `3/2`.

[THEOREM] The reconciliation is that the infinite-system coefficient is finite for fixed `K<K_c^(2D)`, diverges on approaching the critical point from below, and is already infinite at the critical point; finite boxes replace that infinity by the `L^(3/2)` growth measured above.

[EXTERNAL] For `K>K_c^(2D)`, `G(r)` tends to the nonzero squared spontaneous magnetization, so the all-space sum also diverges and an ordinary thermodynamic Taylor expansion at `K_z=0` is not available in the same form.

## Scope and stopping point

[THEOREM] At finite volume the subset identity is exact at every `w`, and every formal coefficient is a finite sum of connected, exactly defined 2D correlation functions.

[COMPUTATION] In the thermodynamic disordered phase this work computes only the `w^2` coefficient and its exact high-temperature prefix through `v^6`.

[COMPUTATION] The expansion is centered at `K_z=0`; its radius of convergence for fixed `K<K_c^(2D)` was not determined.

[COMPUTATION] Nothing here proves that the isotropic point `K_z=K` lies inside that unknown disk, supplies an analytic continuation to it, or constitutes an exact solution of the isotropic three-dimensional model.

[COMPUTATION] Three nonzero `v` coefficients of `c_2` cannot support a D-finite conclusion; no recurrence or differential equation was fitted, and the infinite displacement sum need not inherit the fixed-separation holonomicity of individual 2D correlators.

[COMPUTATION] This route gives no rigorous statement about the anisotropic-to-isotropic crossover beyond the local expansion coefficients at the decoupled-layer point.

[COMPUTATION] Machine-readable exact coefficients, numerical precision, empirical finite-size errors, fit windows, and all passing checks are in `results/interlayer/interlayer_expansion.json`, generated by `experiments/e31_interlayer.py`.
