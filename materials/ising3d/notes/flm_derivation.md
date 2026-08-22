# Finite-lattice derivation for the simple-cubic Ising free energy

## Conventions and exact objects

Throughout,

\[
K=\beta J,\qquad v=\tanh K,\qquad x=e^{-2K},\qquad
\phi=\lim_{N\to\infty}N^{-1}\log Z.
\]

For a finite free box with side lengths \(A=(a_1,\ldots,a_d)\), let \(L(A)\) be the logarithm of the normalized finite-volume polynomial. The high-temperature choice is \(L(A)=\log P_A(v)\). The low-temperature choice is \(L(A)=\log\Xi_A(x)\), with the exterior fixed to \(+1\). All coefficients in the implementation are Python integers or `fractions.Fraction`; no numerical coefficient recognition is used.

## Bounding-box inversion and the finite-difference correction

Let \(W(R)\) collect connected clusters whose smallest axis-aligned bounding box has side lengths \(R=(r_1,\ldots,r_d)\). Such a cluster has

\[
\prod_{i=1}^d(a_i-r_i+1)
\]

translations inside \(A\). Therefore

\[
L(A)=\sum_{R\leq A}\prod_{i=1}^d(a_i-r_i+1)W(R). \tag{1}
\]

Define the backward difference \(\Delta_i f(A)=f(A)-f(A-e_i)\), extending \(f\) by zero when any side length is non-positive. In one coordinate,

\[
\Delta_a(a-r+1)^+={\bf 1}_{a\geq r},\qquad
\Delta_a^2(a-r+1)^+=\delta_{a,r}.
\]

Consequently,

\[
\Delta_1\cdots\Delta_d L(A)=\sum_{R\leq A}W(R), \tag{2}
\]

whereas

\[
W(A)=\Delta_1^2\cdots\Delta_d^2L(A). \tag{3}
\]

This corrects an important possible misreading of the first-difference formula: one must **not** sum the cumulative quantity in (2) over all boxes. With a rectangular cutoff, the partial bulk sum is the single terminal mixed first difference in (2). With the order-bound simplices used here, the implementation performs the equivalent Möbius inversion of (1),

\[
W(A)=L(A)-\sum_{R<A}\prod_i(a_i-r_i+1)W(R), \tag{4}
\]

and then sums each admissible \(W(A)\) once.

The code keeps ordered boxes in (4). Since the isotropic polynomial is invariant under permutations of the sides, transfer calculations are cached by the sorted shape, with the longest side chosen as the transfer direction.

## High-temperature polynomial and exact logarithm

For a free box, `box_broken_bond_poly` returns integers \(c_q\) such that

\[
Z_A=e^{K n_b}\sum_qc_qx^q,
\qquad \sum_qc_q=2^N.
\]

Since

\[
\frac{e^K}{\cosh K}=1+v,
\qquad x=\frac{1-v}{1+v},
\]

comparison with \(Z_A=2^N(\cosh K)^{n_b}P_A(v)\) gives

\[
P_A(v)=2^{-N}\sum_qc_q(1-v)^q(1+v)^{n_b-q}. \tag{5}
\]

Equation (5) is evaluated with integer polynomial arithmetic. For every coefficient actually requested, the numerator is asserted divisible by \(2^N\), and the quotient is asserted non-negative. This is also an implementation-level check of the broken-bond convention. The constant coefficient is asserted to be one.

For \(A(t)=1+\sum_{n\geq1}a_nt^n\) and \(\log A(t)=\sum_{n\geq1}\ell_nt^n\), comparing coefficients in \(A'=(\log A)'A\) gives

\[
\ell_n=a_n-\frac1n\sum_{k=1}^{n-1}k\ell_ka_{n-k}. \tag{6}
\]

The implementation uses (6) with `Fraction` and truncates after every operation.

### Proven HT order bound

A connected even subgraph is Eulerian. In coordinate \(i\), an Euler circuit that visits both extremes of a box of side \(a_i\) must travel at least \(a_i-1\) steps from the minimum to the maximum and at least \(a_i-1\) steps back. Hence it has at least

\[
2\sum_{i=1}^d(a_i-1) \tag{7}
\]

edges. Therefore a box can be omitted through order \(v^n\) when

\[
2\sum_i(a_i-1)>n.
\]

After the finite-lattice interaction series is summed, the complete HT free energy is

\[
\phi_{\rm HT}=\log2+d\log\cosh K+\sum_AW_A(v)
=\log2-\frac d2\log(1-v^2)+\sum_AW_A(v). \tag{8}
\]

## Low-temperature polynomial and surface bound

With a frozen \(+1\) exterior, `box_broken_bond_poly(..., plus_boundary=True)` gives

\[
\Xi_A(x)=\sum_qc_qx^q.
\]

The all-up state is the unique state with no broken bond, so \(c_0=1\); the code asserts this for every box. Fixed exterior bonds make a droplet activity independent of its position relative to a finite-box face. The bulk series is

\[
\phi_{\rm LT}=3K+\sum_AW_A(x). \tag{9}
\]

There is one transfer-engine edge case: for the one-site, one-layer box, the first and last transfer-boundary steps coincide and the engine returns the down-spin exponent as five. The physical box has six frozen-exterior bonds. The series wrapper calls the engine, asserts that exact degenerate output, and replaces only this case by \(\Xi_{1,1,1}=1+x^6\).

### Proven LT order bound

Let a connected set of flipped sites span \(a,b,c\). Its projection onto the \(xy\) plane is a connected grid set spanning \(a\) and \(b\), and therefore contains at least \(a+b-1\) cells. Each occupied \(xy\) column has at least two boundary faces normal to \(z\). Applying the same argument cyclically gives

\[
|\partial S|\geq2[(a+b-1)+(a+c-1)+(b+c-1)]
=4(a+b+c)-6. \tag{10}
\]

A connected Mayer cluster of droplets has connected union. The sum of its droplet surfaces is no smaller than the exterior surface of that union, so (10) remains a lower bound on its total power of \(x\). Thus boxes with

\[
4(a+b+c)-6>n
\]

cannot affect the series through \(x^n\).

## Independent two-dimensional Onsager gate

The square-lattice free energy has the exact double-integral form

\[
\phi_{\rm sq}=\log2+\frac12\left\langle
\log\left[\cosh^2(2K)-\sinh(2K)(\cos p+\cos q)\right]
\right\rangle_{p,q}, \tag{11}
\]

where the brackets denote the normalized integral over \([0,2\pi]^2\). In terms of \(v\),

\[
\phi_{\rm sq}=\log2-\log(1-v^2)+\frac12\langle\log Q\rangle,
\]

\[
Q=(1+v^2)^2-2v(1-v^2)(\cos p+\cos q). \tag{12}
\]

The code expands \(\log Q\) as a bivariate polynomial in \(v\) and \(S=\cos p+\cos q\). It integrates each coefficient exactly from

\[
\langle\cos^{2r}p\rangle=\binom{2r}{r}4^{-r},
\qquad \langle\cos^{2r+1}p\rangle=0. \tag{13}
\]

This calculation shares no finite-box partition polynomials with FLM. It agrees with the 2D FLM result at all twelve nonzero powers \(v^2,v^4,\ldots,v^{24}\):

\[
\begin{aligned}
\phi_{\rm sq}-\log2={}&v^2+\frac32v^4+\frac73v^6+\frac{19}4v^8
+\frac{61}5v^{10}+\frac{75}2v^{12}+\frac{911}7v^{14}\\
&+\frac{3923}8v^{16}+\frac{17629}9v^{18}+\frac{81743}{10}v^{20}
+\frac{388323}{11}v^{22}+\frac{1880881}{12}v^{24}+O(v^{26}).
\end{aligned}
\]

## Derived three-dimensional series

The simple-cubic HT interaction contribution is

\[
\sum_AW_A(v)=3v^4+22v^6+\frac{375}{2}v^8+1980v^{10}
+24044v^{12}+319170v^{14}+\frac{18059031}{4}v^{16}+O(v^{18}).
\]

Including \(3\log\cosh K\),

\[
\begin{aligned}
\phi_{\rm HT}-\log2={}&\frac32v^2+\frac{15}4v^4+\frac{45}2v^6
+\frac{1503}8v^8+\frac{19803}{10}v^{10}+\frac{96177}4v^{12}\\
&+\frac{4468383}{14}v^{14}+\frac{72236127}{16}v^{16}+O(v^{18}).
\end{aligned}
\]

At order \(v^4\), there are three elementary coordinate planes per site and one plaquette of each orientation per site. No lower nonempty even graph exists, so `log P` gives weight \(+1\) to each plaquette and the interaction coefficient is exactly three. The test also enumerates all 375 four-cycles of a \(5^3\) torus (three per site); length five prevents a length-four winding cycle.

The simple-cubic LT result is

\[
\phi_{\rm LT}=3K+x^6+3x^{10}-\frac72x^{12}+15x^{14}-33x^{16}
+\frac{313}{3}x^{18}-\frac{561}{2}x^{20}+O(x^{22}).
\]

The leading coefficient is one because there is one single-site flip per site and it breaks six bonds. The \(-7/2\), rather than \(-3/2\), at \(x^{12}\) is derived, not assumed: the logarithm subtracts an overlapping pair at the same site with weight \(-1/2\) and six nearest-neighbor incompatible monomer pairs per site with total weight \(-3\).

## Validation and achieved orders

- The 2D FLM/Onsager gate passes through \(v^{24}\), twelve nonzero orders.
- The reported 3D HT series reaches \(v^{16}\), using every ordered box with \(\sum_i(a_i-1)\leq8\) (165 ordered boxes). The strong cutoff check adds all 55 ordered boxes at span budget nine and leaves every coefficient through \(v^{16}\) unchanged. Of those, 54 are evaluated through their finite-box transfer polynomials and cancel exactly. The remaining \(4^3\) box exceeds the transfer engine's 62-site limit, but (7) proves its weight starts no earlier than \(v^{18}\), so its retained weight is exactly zero. Extending the requested order from 14 to 16 also preserves the full prefix.
- Reaching \(v^{18}\) requires the free \(4^3\) box. It has 64 sites, while the public exact transfer engine rejects more than 62 sites because it stores coefficients in `int64`; this is the first order blocked by the available engine.
- The 3D LT series reaches \(x^{20}\), using all 20 ordered boxes with \(a+b+c\leq6\). Adding all 15 ordered boxes with side sum seven leaves every coefficient through \(x^{20}\) unchanged, and each added weight is identically zero at the retained orders.
- An exact \(4^3\) torus comparison is not reachable for the same 64-site `int64` guard (and the torus transfer matrix has a 16-spin cross-section). No finite-torus numerical comparison is reported or implied.
- The JSON artifacts record exact rational strings and precision as exact integer/`Fraction` arithmetic. No published critical benchmark was used to choose or fit any coefficient.
