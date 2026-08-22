# Onsager--Kaufman control derivation

This note fixes the conventions used by `ising.onsager` and explains why the two-dimensional
transfer matrix closes on a finite-dimensional free-fermion algebra. The purpose is an
**algebraic control experiment** for the three-dimensional project: here every step can be
checked against exact finite-lattice integers.

The coupling convention is the repository convention

\[
K_x=\beta J_x,\qquad K_y=\beta J_y,
\]

and the reduced free energy is

\[
\phi(K_x,K_y)=\lim_{mn\to\infty}\frac{1}{mn}\log Z_{m,n}.
\]

The formulas below are for ferromagnetic couplings and an \(m\times n\) torus. A row has
\(n\) spins, its bonds carry \(K_x\), and the transfer direction has \(m\) rows with bonds
carrying \(K_y\). High-precision partition functions and integrals used `mpmath` at 60 decimal
digits; the explicit `expm`/`logm` Clifford diagnostics used SciPy `complex128`. The
machine-readable record is `results/onsager/onsager_2d.json`.

## 1. The symmetric row transfer matrix

Write \(X_j,Y_j,Z_j\) for Pauli matrices on row site \(j\). For two adjacent rows, the
one-site vertical Boltzmann matrix in the \(Z\)-spin basis is

\[
A_y=
\begin{pmatrix}
e^{K_y}&e^{-K_y}\\
e^{-K_y}&e^{K_y}
\end{pmatrix}.
\]

Define the dual coupling \(K_y^*\) by

\[
\tanh K_y^*=e^{-2K_y}.
\]

Comparison of diagonal and off-diagonal entries, together with
\(\det A_y=2\sinh(2K_y)\), gives the exact identity

\[
A_y=\sqrt{2\sinh(2K_y)}\,e^{K_y^*X}.
\]

Consequently

\[
V_1=A_y^{\otimes n}
  =[2\sinh(2K_y)]^{n/2}
    \exp\!\left(K_y^*\sum_{j=1}^{n}X_j\right),
\]

while the interactions inside one periodic row are

\[
V_2=\exp\!\left(K_x\sum_{j=1}^{n}Z_jZ_{j+1}\right),
\qquad Z_{n+1}=Z_1.
\]

The symmetric positive transfer matrix is

\[
V=V_2^{1/2}V_1V_2^{1/2},
\qquad Z_{m,n}=\operatorname{Tr}V^m.
\]

The symmetric placement is a similarity-equivalent version of the usual nonsymmetric transfer
matrix, but it makes the principal matrix logarithm unambiguous and Hermitian.

## 2. Jordan--Wigner Majoranas and closure of the algebra

The code uses an \(X\)-string Jordan--Wigner convention. For \(j=1,\ldots,n\), set

\[
 c_{2j-1}=\left(\prod_{\ell<j}X_\ell\right)Z_j,
 \qquad
 c_{2j}=-\left(\prod_{\ell<j}X_\ell\right)Y_j.
\]

These matrices are Hermitian and obey the Clifford relations

\[
\{c_a,c_b\}=2\delta_{ab}I.
\]

The minus sign in the definition of \(c_{2j}\) is deliberate. Direct Pauli multiplication gives

\[
\boxed{X_j=-i c_{2j-1}c_{2j}},
\qquad
\boxed{Z_jZ_{j+1}=-i c_{2j}c_{2j+1}}
\quad (j<n).
\]

Thus each generator appearing in \(V_1\) or in a non-wrapping part of \(V_2\) is a Majorana
bilinear. Commutators of bilinears are again linear combinations of bilinears. Equivalently, the
operators \(c_ac_b\), \(a<b\), furnish the spin representation of the Lie algebra
\(\mathfrak{so}(2n)\). Products and logarithms of their exponentials therefore remain in the
corresponding Gaussian, or matchgate, group (apart from the scalar prefactor in \(V_1\)). This is
the precise algebraic reason the two-dimensional calculation closes.

### The periodic-boundary qualification

For a periodic row, the last bond must not be silently treated as an ordinary local
Jordan--Wigner bond. Let

\[
P=\prod_{j=1}^n X_j
\]

be the spin-flip operator, which is also fermion parity in this convention. The wrap bond obeys

\[
-i c_{2n}c_1=-P Z_nZ_1.
\]

In the parity sector \(P=p\), \(p=\pm1\), it follows that

\[
Z_nZ_1=p\,i c_{2n}c_1.
\]

Therefore each fixed-parity extension \(V_p\) is Gaussian, but the unreduced periodic spin
matrix on the full \(2^n\)-dimensional space is generally **not one Gaussian operator**. It is a
direct sum of two Gaussian parity blocks with opposite fermionic boundary conditions. This
qualification is essential: omitting it gives the right local identities but the wrong finite
torus signs.

The numerical diagnostic makes this distinction visible. At \(K_x=0.31\), \(K_y=0.47\), and
\(n=3,4,5,6\), both parity-resolved logarithms have maximum off-bilinear Pauli coefficient

\[
1.246120990682649\times10^{-15}.
\]

By contrast, the unresolved physical matrices have an off-bilinear coefficient at least
\(0.34209468266573206\), carried by the Jordan--Wigner wrap term. The latter is not a numerical
failure; it is the parity-boundary obstruction just derived.

A “Majorana bilinear” means weight at most two in the Majorana generators. Its Pauli-string image
can contain a long Jordan--Wigner \(X\) string, so it need not have ordinary Pauli locality two.
The code classifies support using the actual matrices \(i c_ac_b\), not the number of nonidentity
Pauli factors.

## 3. Extracting the quadratic generator and its energies

For each parity-resolved explicit matrix, the code computes

\[
M_p=\operatorname{logm}(V_p)
   =\alpha_p I+\frac{i}{4}c^{\mathsf T}A_pc,
\]

where \(A_p\) is a real antisymmetric \(2n\times2n\) matrix. Orthogonality of Pauli strings gives
an especially direct extraction formula:

\[
\alpha_p=2^{-n}\operatorname{Tr}M_p,
\qquad
(A_p)_{ab}=2\,2^{-n}\operatorname{Tr}\!\left[(i c_ac_b)M_p\right],
\quad a<b.
\]

The Hermitian matrix \(iA_p\) has eigenvalues

\[
-\gamma_1,\ldots,-\gamma_n,+\gamma_1,\ldots,+\gamma_n.
\]

The full Pauli-basis expansion of `logm(V_p)` was performed for every \(n=3,4,5,6\). The largest
coefficient outside the scalar and bilinear basis is the \(1.2461\times10^{-15}\) number above,
and the largest difference between extracted and analytic single-particle energies is
\(6.439\times10^{-15}\).

The traces can be reconstructed from the \(2n\)-dimensional generator rather than the
\(2^n\)-dimensional transfer matrix:

\[
\operatorname{Tr}V_p^m
 =e^{m\alpha_p}\prod_{r=1}^n2\cosh\frac{m\gamma_r}{2}.
\]

The parity-inserted trace also retains the orientation of the antisymmetric generator:

\[
\operatorname{Tr}(P V_p^m)
 =(-1)^n\operatorname{sgn}\operatorname{Pf}(A_p)
  e^{m\alpha_p}\prod_{r=1}^n2\sinh\frac{m\gamma_r}{2}.
\]

Projection with \((I+pP)/2\), followed by summing \(p=+1,-1\), reconstructs the spin-torus
partition function. For \(n=m=3,4,5,6\), the largest relative difference from the exact integer
transfer-matrix evaluation was \(5.66194641473\times10^{-14}\). This check uses double-precision
`scipy.linalg.logm`; the separate Kaufman-versus-integer check below uses 60-digit arithmetic.

## 4. Dispersion relation and the four Kaufman products

The two parity sectors impose the two fermionic momentum grids

\[
q_r^{(o)}=\frac{(2r+1)\pi}{n},
\qquad
q_r^{(e)}=\frac{2r\pi}{n},
\qquad r=0,\ldots,n-1.
\]

The one-particle dispersion is

\[
\cosh\gamma(q)
 =\cosh(2K_x)\cosh(2K_y^*)
  -\sinh(2K_x)\sinh(2K_y^*)\cos q.
\]

All nonzero-momentum energies use the nonnegative principal `acosh`. The periodic zero mode must
instead retain its sign:

\[
\boxed{\gamma_0=2(K_y^*-K_x)}.
\]

Define

\[
\begin{aligned}
C_o&=\prod_{r=0}^{n-1}2\cosh\!\left(\frac{m\gamma(q_r^{(o)})}{2}\right),&
S_o&=\prod_{r=0}^{n-1}2\sinh\!\left(\frac{m\gamma(q_r^{(o)})}{2}\right),\\
C_e&=\prod_{r=0}^{n-1}2\cosh\!\left(\frac{m\gamma(q_r^{(e)})}{2}\right),&
S_e&=\prod_{r=0}^{n-1}2\sinh\!\left(\frac{m\gamma(q_r^{(e)})}{2}\right).
\end{aligned}
\]

The empirically fixed finite-torus convention is

\[
\boxed{
Z_{m,n}=\frac12[2\sinh(2K_y)]^{mn/2}
\left(C_o+S_o+C_e-S_e\right).
}
\]

The plus-parity block supplies the half-integer/odd grid, and the minus-parity block supplies the
integer/even grid. Above the transition \(\gamma_0>0\). Below it \(\gamma_0<0\), so \(S_e\)
changes sign while \(C_e\) does not. Replacing \(\gamma_0\) by `abs(gamma_0)` loses this phase
information and gives the wrong low-temperature partition function.

### How the convention was fixed

The finite integers, not a published thermodynamic benchmark, were used as ground truth. The
experiment propagates exact bivariate coefficients

\[
c_{m,n}(q_x,q_y)\in\mathbb Z_{\ge0}
\]

and evaluates

\[
Z_{m,n}=e^{mn(K_x+K_y)}
\sum_{q_x,q_y}c_{m,n}(q_x,q_y)e^{-2K_xq_x-2K_yq_y}.
\]

For every shape, diagonal collapse

\[
c(q)=\sum_{q_x+q_y=q}c(q_x,q_y)
\]

agreed coefficient-by-coefficient with
`ising.transfer_matrix.torus_broken_bond_poly((n,m))`. This also explains why a separate
bivariate propagation is needed for \(K_x\ne K_y\): the repository function intentionally stores
only the total broken-bond count \(q\), which is sufficient only after setting the two couplings
equal.

Eight cases were checked at 60 digits: square tori of sizes 3, 4, 5, and 6, and rectangular tori
\(3\times4\), \(4\times5\), \(5\times6\), and \(6\times3\), including four anisotropic
couplings. The largest relative error was

\[
3.6386727631733\times10^{-60}.
\]

Nearby conventions were explicitly falsified. At generic anisotropic and low-temperature points,
using all plus signs, changing the odd-sinh sign, swapping which momentum grid receives the minus
sign, or forcing \(\gamma_0\ge0\) produced relative errors respectively of order
\(0.6293\), \(0.6606\), \(0.03131\), and \(0.4452\). Thus the surviving convention is exactly
`C_odd + S_odd + C_even - S_even` with signed \(\gamma_0\).

## 5. Thermodynamic limit and normalization of the double integral

For fixed \(n\) and large \(m\), each sector product has
\(\log[2\cosh(m\gamma/2)]\sim m\gamma/2\). Taking both dimensions to infinity gives

\[
\phi(K_x,K_y)
 =\frac12\log[2\sinh(2K_y)]
  +\frac{1}{4\pi}\int_0^{2\pi}\gamma(q)\,dq.
\]

Using

\[
\cosh\gamma(q)
 =\frac{\cosh(2K_x)\cosh(2K_y)-\sinh(2K_x)\cos q}{\sinh(2K_y)}
\]

and the elementary angular identity

\[
\frac{1}{2\pi}\int_0^{2\pi}\log(A-B\cos\theta)\,d\theta
 =\log\frac{A+\sqrt{A^2-B^2}}{2},
\qquad A\ge |B|,
\]

one obtains

\[
\boxed{
\phi=\log 2+\frac{1}{8\pi^2}
\int_0^{2\pi}\!\int_0^{2\pi}
\log\!\left[
\cosh(2K_x)\cosh(2K_y)
-\sinh(2K_x)\cos\theta_1
-\sinh(2K_y)\cos\theta_2
\right]d\theta_1d\theta_2.
}
\]

Thus the normalization stated in the assignment is the one consistent with the finite integer
data; no extra factor is required. Three independent checks fix the additive and multiplicative
normalization:

1. At \(K_x=K_y=0\), the logarithm vanishes and \(\phi=\log2\), the free-spin value.
2. At the isotropic critical point, numerical integration gives
   \[
   \phi_c=0.929695398341610214985384973665878121765882375151802167582343,
   \]
   agreeing to \(7.78\times10^{-62}\) with the independent evaluation
   \(\frac12\log2+2G/\pi\), where \(G\) is Catalan's constant.
3. At criticality the finite Kaufman values approach the integral as follows:

| \(L=m=n\) | \(\log Z_{L,L}/L^2-\phi_c\) |
|---:|---:|
| 64 | \(1.5623335267914998981\times10^{-4}\) |
| 128 | \(3.9057428579247241253\times10^{-5}\) |
| 256 | \(9.7643003482629668656\times10^{-6}\) |

The residual falls by a factor approaching four, as expected for the leading \(L^{-2}\) torus
correction at criticality.

The implementation analytically integrates one angle. To avoid subtracting nearly equal numbers,
it uses

\[
\cosh(2K_x)\cosh(2K_y)-\sinh(2K_x)-\sinh(2K_y)
=\frac{[\sinh(2K_x)\sinh(2K_y)-1]^2}
 {\cosh(2K_x)\cosh(2K_y)+\sinh(2K_x)+\sinh(2K_y)}.
\]

This form remains stable exactly on the critical line.

## 6. Location and nature of the singularity

The determinant inside the double integral is minimized at
\(\theta_1=\theta_2=0\). The identity just displayed proves that its minimum vanishes if and
only if

\[
\boxed{\sinh(2K_x)\sinh(2K_y)=1.}
\]

The same condition is \(K_x=K_y^*\), hence \(\gamma_0=0\). It is therefore simultaneously the
zero of the Kac--Ward determinant and the closing of the fermion gap. No critical constant was
fit: in the isotropic case the equation itself gives

\[
K_c=\frac12\operatorname{arsinh}(1)
   =\frac12\log(1+\sqrt2).
\]

Expanding the determinant around zero momentum gives a mass squared plus a positive quadratic
form in \(\theta_1,\theta_2\). Integrating its logarithm produces the familiar
\(t^2\log|t|\) singular part of \(\phi\), so the second derivative diverges logarithmically. As a
numerical check, centered second differences along the isotropic line at step sizes
\(0.002,0.001,0.0005\) were

\[
15.98131166298,\qquad17.74628135550,\qquad19.51133362917,
\]

and increase under each halving.

## 7. Reproduction

From the repository root, run

```text
.venv/bin/python experiments/e02_onsager_2d.py
.venv/bin/python tests/test_onsager.py
```

Both commands print a final `PASS`. The experiment overwrites
`results/onsager/onsager_2d.json` with the timestamp, 60-digit precision setting, all finite-size
values, sector diagnostics, and individual check outcomes.

## Primary historical references

- L. Onsager, “Crystal Statistics. I. A Two-Dimensional Model with an Order-Disorder
  Transition,” *Physical Review* **65** (1944), 117–149,
  <https://doi.org/10.1103/PhysRev.65.117>.
- B. Kaufman, “Crystal Statistics. II. Partition Function Evaluated by Spinor Analysis,”
  *Physical Review* **76** (1949), 1232–1243,
  <https://doi.org/10.1103/PhysRev.76.1232>.
