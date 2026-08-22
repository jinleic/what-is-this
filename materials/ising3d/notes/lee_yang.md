# Lee–Yang zeros from exact finite lattices

## Scope and conventions

This track studies the field variable

\[
z=e^{-2H_f},\qquad H_f=\beta h,
\]

at fixed ferromagnetic coupling \(K=\beta J\).  A spin configuration is graded by

- \(k=\#\{i:\sigma_i=-1\}\), the number of **down** spins;
- \(q=\#\{(ij):\sigma_i\ne\sigma_j\}\), the number of unsatisfied bonds;
- \(x=e^{-2K}\).

Since \(M=N-2k\) and \(\sum_{(ij)}\sigma_i\sigma_j=n_b-2q\),

\[
\begin{aligned}
Z(K,H_f)
 &=\sum_{\sigma}e^{K(n_b-2q)+H_f(N-2k)}\\
 &=e^{Kn_b+H_fN}\sum_{k=0}^{N}\sum_{q=0}^{n_b}c[k,q]x^qz^k.
\end{aligned}
\]

The exact combinatorial object is therefore

\[
C(x,z)=\sum_{k=0}^{N}A_k(x)z^k,
\qquad A_k(x)=\sum_q c[k,q]x^q\in\mathbb Z[x],
\]

with nonnegative integer \(c[k,q]\).  At a rational coupling fugacity \(x=p/r\), multiplication by the irrelevant common factor \(r^{n_b}\) gives the literal integer polynomial

\[
\widetilde C_{p/r}(z)=\sum_k\widetilde A_kz^k,
\qquad
\widetilde A_k=\sum_q c[k,q]p^q r^{n_b-q}\in\mathbb Z.
\]

This distinction matters: at a generic real \(K\), \(A_k(e^{-2K})\) is not an integer, while the bivariate coefficient table \(c[k,q]\) is exact.  All root calculations in this track start either from that exact table or from a CRT-reconstructed integer \(\widetilde C_{p/r}\).

Global spin reversal preserves \(q\) and sends \(k\mapsto N-k\).  Hence the exact palindrome is

\[
c[k,q]=c[N-k,q],\qquad A_k(x)=A_{N-k}(x).
\]

At \(z=1\), summing over \(k\) removes the field grading:

\[
\sum_k c[k,q]=\operatorname{dos\_bonds}[n_b-q].
\]

Both identities are checked exactly, coefficient by coefficient.

## Exact engines and lattice reach

Three independent paths were used.

1. `field_dos_from_enumeration` reduces `joint_dos(lat, with_field=True)` to the exact \(c[k,q]\) table.  It was used for open and periodic cubic \(2\times2\times2\), \(3\times3\times2\), and \(3\times3\times3\) lattices, and open and periodic square \(4\times4\) and \(5\times5\) controls.
2. `transfer_field_dos` extends the layer transfer matrix with a second polynomial axis for \(k\).  It produced exact open and periodic \(6\times6\) square tables.  The periodic case is an exact trace; neither the shared transfer-matrix module nor its conventions were changed.
3. `rational_field_polynomial_transfer` keeps the field-polynomial axis but evaluates \(x=p/r\) exactly modulo a sequence of 30-bit primes.  A coefficient bound determines when the CRT modulus product is sufficient.  This reached the open \(4\times4\times4\) cube: \(N=64\), \(n_b=144\), degree 64 in \(z\).  Six exact integer polynomials were computed, for \(x=7/10,2/3,9/14,5/8,3/5,1/2\).

A periodic \(4^3\) trace was not attempted: the direct trace representation carries both current and initial \(2^{16}\)-state indices and would scale as \(4^{16}\) before polynomial axes.  The largest reached lattice is therefore the **open \(4\times4\times4\) cube (64 sites)**, not a periodic cube.

The full integer coefficient tables and the six degree-64 polynomials are in `results/lee_yang/lee_yang_analysis.json`.

## The Lee–Yang theorem, stated with its hypotheses

Let \(G=(V,E)\) be a finite graph with pairwise Ising couplings \(K_{ij}=\beta J_{ij}\ge0\).  Give each site a possibly complex field \(H_i\), and put \(z_i=e^{-2H_i}\).  After factoring the all-up Boltzmann weight, the multivariate field polynomial is

\[
P_G(\{z_i\})=
\sum_{S\subseteq V}
\left(\prod_{(ij)\in\delta S}e^{-2K_{ij}}\right)
\left(\prod_{i\in S}z_i\right),
\]

where \(\delta S\) is the set of edges with one endpoint in \(S\).

**Multivariate Lee–Yang zero-free statement.**  If all \(K_{ij}\ge0\), then

\[
|z_i|<1\quad\text{for every }i
\quad\Longrightarrow\quad
P_G(\{z_i\})\ne0.
\]

For a connected positive-coupling graph, the stronger closed-polydisk form says that \(|z_i|\le1\) for all \(i\), with at least one strict inequality, is also zero-free.  If zero couplings disconnect the positive-coupling graph, that stronger statement is applied component by component.  The boundary on which every \(|z_i|=1\) is not claimed to be zero-free.

For a uniform field \(z_i=z\), spin reversal gives

\[
P_G(z)=z^N P_G(1/z).
\]

The zero-free open unit disk excludes \(|z|<1\); reciprocity excludes \(|z|>1\).  Therefore every zero of the finite-volume uniform-field polynomial satisfies

\[
|z_j|=1.
\]

The ferromagnetic sign is essential: the argument does not extend to antiferromagnetic or mixed-sign couplings.

### Corollary for a real field: no nonzero real critical field

For real \(h>0\), \(H_f=\beta h>0\) and \(z=e^{-2H_f}\in(0,1)\).  For real \(h<0\), \(z>1\), which is zero-free by spin reversal.  Thus every finite-volume partition function is nonzero for real \(h\ne0\), and a consistent branch of \(\log Z_\Lambda\) is analytic in each complex half-plane \(\Re H_f>0\) and \(\Re H_f<0\).

For a finite-range, or sufficiently summable, ferromagnetic Ising interaction along a van Hove sequence, the pressure limit exists.  Uniform finite-volume bounds make the analytic pressures a normal family; convergence on the real half-axis and the Vitali/normal-family argument extend the thermodynamic limit analytically throughout each half-plane.  Consequently the reduced free energy is real analytic for every real \(h\ne0\).

**Corollary.**  Under those ferromagnetic and thermodynamic-limit hypotheses, there is no phase transition at a nonzero real magnetic field.  A real-field singularity can occur only at \(h=0\).

This conclusion depends on the hypotheses.  The circle theorem alone says nothing for mixed-sign/antiferromagnetic couplings, nonsummable interactions for which the pressure construction fails, or a different field coupling.  It also does not say that the complex-field edge is analytic: the Yang–Lee edge is precisely a thermodynamic-limit singularity at nonzero *imaginary* field.

## Numerical circle verification

Roots were computed with 100 decimal digits.  Palindromy reduces an even-degree polynomial to a real Chebyshev polynomial in

\[
t=\frac{z+z^{-1}}2=\cos\theta.
\]

Float64 Chebyshev roots are used only as isolating guesses.  Each \(t\) is refined at 100 digits, mapped to a complex initial guess, deliberately perturbed in the radial direction, and then independently refined against the original complex polynomial.  The final radius is therefore measured by complex Newton iteration; it is not imposed by returning \(e^{i\theta}\).  For odd \(N\), the exact palindromic factor \(z+1\) is removed first and its root is restored.

Across all open and periodic required lattices, all seven coupling values for the \(N\le36\) tables, and all six exact 64-site rational polynomials:

- every polynomial returned exactly \(N\) roots;
- the largest observed \(\max_j\lvert |z_j|-1\rvert\) was \(9.1045\times10^{-92}\);
- the largest normalized polynomial residual was \(6.2279\times10^{-101}\);
- the acceptance threshold was \(10^{-40}\).

A circle violation was treated as an algorithmic error condition, not as physics.

## First zero, local density, and finite-size behavior

For roots \(z_j=e^{i\theta_j}\) on the positive arc, \(0<\theta_1<\theta_2<\cdots\le\pi\).  The finite-volume near-edge density estimator recorded in the JSON is

\[
\rho_{j+1/2}=\frac{1}{N(\theta_{j+1}-\theta_j)}
\]

at the interval midpoint.  The first four intervals are retained for every lattice and every coupling.

For the open cubic sequence \(L=2,3,4\), the exact rational-grid results are:

| \(x\) | \(K=-\tfrac12\log x\) | \(\theta_1(L=2)\) | \(\theta_1(L=3)\) | \(\theta_1(L=4)\) | first \(L=4\) density |
|---:|---:|---:|---:|---:|---:|
| 7/10 | 0.1783374720 | 1.1783358153 | 0.6382832388 | 0.4223254135 | 0.1009661 |
| 2/3 | 0.2027325541 | 1.0801525467 | 0.5350689078 | 0.3236162072 | 0.1009637 |
| 9/14 | 0.2209163761 | 1.0144851528 | 0.4690918944 | 0.2632257504 | 0.1016685 |
| 5/8 | 0.2350018146 | 0.9675592288 | 0.4237333096 | 0.2233997475 | 0.1026093 |
| 3/5 | 0.2554128119 | 0.9051093485 | 0.3660944200 | 0.1755016657 | 0.1045994 |
| 1/2 | 0.3465735903 | 0.6921383367 | 0.2032758352 | 0.0719716887 | 0.1233986 |

At the separately evaluated benchmark \(K=0.221654626\), the exact \(c[k,q]\) tables give \(\theta_1=1.0119432504\) for the open \(L=2\) cube and \(0.4665943500\) for open \(L=3\).  No value at this exact \(K\) was synthesized for \(L=4\); \(x=9/14\) is explicitly a nearby independent grid point, not a stand-in fitted to the benchmark.

A naive fit \(\theta_1\propto L^{-y_{\mathrm{eff}}}\) on the open \(L=2,3,4\) sequence gives \(y_{\mathrm{eff}}=1.48,1.74,1.94,2.11,2.36,3.25\) in the table order.  These are diagnostics, not exponent estimates:

- for \(K<K_c\) (high temperature), \(\theta_1\) should approach a nonzero edge, so a zero-offset power law is inappropriate asymptotically;
- at \(K_c\), one expects critical magnetic scaling;
- for \(K>K_c\), zeros pinch the origin and the first-order coexistence regime ultimately gives volume scaling \(\theta_1\sim L^{-d}\).

The trend toward a power near 3 at the largest sampled \(K\) is qualitatively consistent with the last statement, but \(L\le4\) is too small for a quantitative exponent claim.

## High-temperature edge and the exponent \(\sigma\)

For \(K<K_c\), the thermodynamic zero density is expected to start at a positive angle and behave as

\[
g(\theta)\sim(\theta-\theta_{\mathrm{edge}})^\sigma,
\qquad \theta\downarrow\theta_{\mathrm{edge}}.
\]

Integrating this form motivates the finite-size quantile model

\[
\theta_j(L)=\theta_{\mathrm{edge}}+
C\left(\frac{j-1/2}{N}\right)^{1/(\sigma+1)}.
\]

The fit used the first four positive zeros from open cubes \(L=2,3,4\).  Sensitivity ranges vary the window through \(j_{\max}=3,4,5\) and repeat without \(L=2\).  They are systematic envelopes, not confidence intervals.

| \(K\) | provisional \(\theta_{\mathrm{edge}}\) | edge sensitivity range | naive \(\sigma\) | \(\sigma\) sensitivity range |
|---:|---:|---:|---:|---:|
| 0.1783374720 | 0.2772 | 0.2020–0.3326 | 0.471 | 0.328–0.649 |
| 0.2027325541 | 0.1816 | 0.1126–0.2352 | 0.438 | 0.300–0.598 |
| 0.2209163761 | 0.1244 | 0.0604–0.1762 | 0.412 | 0.280–0.560 |

The decreasing provisional edge as \(K\) approaches the quoted critical coupling is the supported qualitative result.  The exponent is **not** supported quantitatively.  Modern calculations place the three-dimensional Yang–Lee edge exponent near \(0.08\); for example, An, Mesterházy, and Stephanov report \(\sigma=0.0742(56)\).  The small-cube estimates above are far away and strongly window dependent.

More importantly, applying the identical pipeline to periodic square controls \(L=4,5,6\) yields positive naive values \(\sigma=0.0525,0.0657,0.0748\) at the three high-temperature grid points, rather than the exact two-dimensional \(\sigma=-1/6\).  That failed control shows that these volumes have not reached the edge asymptotic regime and that formal least-squares errors severely understate the uncertainty.  The responsible conclusion is therefore:

> The data provisionally locate a positive high-temperature edge and its movement toward zero, but they do not provide a defensible three-dimensional estimate of \(\sigma\).  The uncertainty is dominated by uncontrolled finite-size and boundary corrections, not root-finding precision.

No published critical coupling or exponent was used as a fitted constant.

## Reproduction

```text
.venv/bin/python tests/test_lee_yang.py
.venv/bin/python experiments/e04_lee_yang.py
```

The test prints `PASS` after checking the circle to better than \(10^{-40}\), the exact palindrome, the \(z=1\) zero-field reduction, root count, and the 64-site reach.  The experiment writes `results/lee_yang/lee_yang_analysis.json` with precision, timestamp, exact coefficients, all positive zero angles, local densities, scaling fits, and check records.

## Sources

1. T. D. Lee and C. N. Yang, “Statistical Theory of Equations of State and Phase Transitions. II. Lattice Gas and Ising Model,” *Physical Review* **87**, 410–419 (1952), [doi:10.1103/PhysRev.87.410](https://doi.org/10.1103/PhysRev.87.410).
2. G. Harcos, “The Lee–Yang Circle Theorem,” a concise proof of the multivariate zero-free form, [PDF](https://users.renyi.hu/~gharcos/lee-yang.pdf).
3. M. E. Fisher, “Yang-Lee Edge Singularity and \(\phi^3\) Field Theory,” *Physical Review Letters* **40**, 1610–1613 (1978), [doi:10.1103/PhysRevLett.40.1610](https://doi.org/10.1103/PhysRevLett.40.1610).
4. X. An, D. Mesterházy, and M. A. Stephanov, “Functional renormalization group approach to the Yang-Lee edge singularity,” *JHEP* **07** (2016) 041, [arXiv:1605.06039](https://arxiv.org/abs/1605.06039), especially Tables II–III.
