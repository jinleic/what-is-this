# Claimed exact solutions and decisive checks

This note treats every formula below as an **unverified claim**.  Throughout,

\[
\Phi \equiv \lim_{N\to\infty}\frac{1}{N}\log Z=-\beta f
\]

is the reduced free energy per site.  Equation numbers are those printed in the linked arXiv PDFs, not renumbered here.  The comparison value used below is \(K_c=0.221654626(5)\), quoted in Perk's comment on p. 3 of [arXiv:2202.03136](https://arxiv.org/pdf/2202.03136).

## 1. Z.-D. Zhang (2007): four-fold-integral conjecture

**Primary claim:** [arXiv:0705.1045](https://arxiv.org/pdf/0705.1045), *Conjectures on exact solution of three-dimensional (3D) simple orthorhombic Ising lattices*.

### Claimed free energy

Zhang defines the extra coupling in Eq. (3.6) as

\[
K'''=\frac{K'K''}{K}.
\]

The orthorhombic reduced free energy printed as Eq. (3.37) is

\[
\begin{aligned}
\Phi={}&\log 2+\frac{1}{2(2\pi)^4}
\int_{-\pi}^{\pi}d\omega'\int_{-\pi}^{\pi}d\omega_x
\int_{-\pi}^{\pi}d\omega_y\int_{-\pi}^{\pi}d\omega_z\\
&\quad\times\log\!\Bigl[
 \cosh(2K)\cosh\!\bigl(2(K'+K''+K''')\bigr)
 -\sinh(2K)\cos\omega'\\
&\hspace{42mm}
 -\sinh\!\bigl(2(K'+K''+K''')\bigr)
  \bigl(w_x\cos\omega_x+w_y\cos\omega_y+w_z\cos\omega_z\bigr)
\Bigr].
\end{aligned}
\tag{0705.1045, Eq. 3.37}
\]

For the isotropic cubic case \(K'=K''=K\), hence \(K'''=K\), the argument becomes

\[
\cosh(2K)\cosh(6K)-\sinh(2K)\cos\omega'
-\sinh(6K)(w_x\cos\omega_x+w_y\cos\omega_y+w_z\cos\omega_z).
\]

The weight prescription is not a closed function determined independently of already-known data.  Appendix A gives \(w_x=1\) and

\[
w_y=w_z=\pm\sqrt{\sum_{i=0}^{\infty}b_i\kappa^{2i}},
\qquad \kappa=\tanh K,
\tag{0705.1045, Eq. A.2}
\]

where the displayed \(b_i\) were fitted from known high-temperature coefficients.  Under “Ansatz 1,” the paper then declares \(w_y=w_z=0\) at every finite temperature because the truncated radicand becomes negative; see the text immediately after Eq. (A.2), especially pp. 137–138 of the arXiv PDF.  This is the finite-temperature choice used for the subsequent critical-point and thermodynamic claims.

There is a source-level normalization conflict, not merely an OCR ambiguity:

* Eq. (3.37) has \(2(2\pi)^4\) in the denominator.
* The isotropic restatement Eq. (3.62), despite still displaying four integrations, has \(2(2\pi)^3\).
* Appendix Eq. (A.1) again has \(2(2\pi)^4\).

These exponents were checked visually in the arXiv PDF.  Consequently this claim is marked `uncertain: true` in the JSON.  The machine expression uses the mutually agreeing Eqs. (3.37) and (A.1), specializes to \(w_x=1,w_y=w_z=0\), and analytically performs the now-trivial \(\omega_y,\omega_z\) integrations:

```text
mp.log(2) + mp.quad(lambda op: mp.quad(lambda ox: mp.log(mp.cosh(2*K)*mp.cosh(6*K) - mp.sinh(2*K)*mp.cos(op) - mp.sinh(6*K)*mp.cos(ox)), [-mp.pi, mp.pi]), [-mp.pi, mp.pi]) / (2*(2*mp.pi)**2)
```

For completeness, following the literal \((2\pi)^3\) normalization printed in Eq. (3.62) would change only the last denominator of that reduced two-dimensional integral to `2*(2*mp.pi)`.  No choice can make Eqs. (3.62) and (A.1) identical as printed.

Wu, McCoy, Fisher, and Chayes refer to the corresponding published-paper formulas as Eq. (49), reducing to Eq. (74) in the isotropic case; see pp. 2–3 of [arXiv:0811.3876](https://arxiv.org/pdf/0811.3876).  The arXiv manuscript's corresponding numbered formulas are Eqs. (3.37), (3.62), and (A.1).

### Claimed critical point

The source states the general condition

\[
\sinh(2K)\,\sinh\!\bigl(2(K'+K''+K''')\bigr)=1
\tag{0705.1045, Eq. 3.41}
\]

and, after using Eq. (3.6) and setting \(K'=K''=K\), **also states the cubic condition directly**:

\[
\boxed{\sinh(2K_c)\sinh(6K_c)=1.}
\tag{0705.1045, Eq. 3.60}
\]

Thus the cubic condition both follows from the more general relation and appears explicitly in the source.  The paper also gives \(K_c=0.24060591\ldots\) in Eq. (3.56) and

\[
e^{-2K_c}=\frac{\sqrt5-1}{2},
\qquad
K_c=-\frac12\log\!\left(\frac{\sqrt5-1}{2}\right).
\]

A fresh 80-decimal calculation of that exact expression gives

\[
K_c=0.240605912529801723748879456712184\ldots,
\]

or **\(0.240605912529802\)** rounded to 15 significant decimal digits.  It exceeds \(0.221654626\) by \(0.018951286529802\), about \(8.55\%\).

Directly evaluable critical-point expression:

```text
-mp.log((mp.sqrt(5)-1)/2)/2
```

### Published objections

1. **Unknown weights followed by an incompatible ad hoc choice.**  Wu et al., [arXiv:0811.3876](https://arxiv.org/pdf/0811.3876), pp. 2–3, observe that the four-fold integral contains undetermined \(w_x,w_y,w_z\).  Appendix Eq. (A.2) is fitted to the known high-temperature series, but the paper then sets \(w_y=w_z=0\) just before the published Eq. (50) and uses that different choice thereafter.  The critical relation therefore inherits the unsupported weight choice.
2. **The first high-temperature discrepancy is already \(\kappa^2\).**  Zhang's own Eq. (A.12) gives the known result
   \[
   \lambda=2\cosh^3K\,[1+3\kappa^4+22\kappa^6+\cdots],
   \]
   whereas the finite-temperature expression Eq. (A.13) gives
   \[
   \lambda=2\cosh^3K\,[1+\tfrac12\kappa^2+\tfrac{87}{8}\kappa^4+\cdots].
   \]
   Perk identifies this exact \(\kappa^2\) mismatch in Remark 2.5 of [arXiv:1209.0731](https://arxiv.org/pdf/1209.0731).  Section 3.1 further calls Eq. (A.2) reverse engineering of known coefficients rather than a prediction.
3. **Analyticity rules out two different high-temperature series.**  Perk's Theorems 2 and 3 in arXiv:1209.0731 prove uniform analyticity of correlations and reduced free energy around \(K=0\), using the Suzuki identity (9)–(12).  Equation (17) supplies explicit nonzero lower bounds \(|\tanh K|<0.131652497\ldots\) and \(|K|<0.130899693\ldots\).  Therefore a special “infinitesimal-temperature-distance” series and a different finite-temperature series cannot describe the same thermodynamic free energy.
4. **The magnetization fails at its first possible low-temperature term.**  Wu et al., p. 4, visually and unambiguously print Zhang's Eq. (103) expansion as
   \[
   1-6x^8-12x^{10}-18x^{12}-\cdots,
   \]
   versus the exact cubic expansion
   \[
   1-2x^6-12x^{10}+14x^{12}-\cdots,
   \qquad x=e^{-2K}.
   \]
   Thus the first disagreement is at \(x^6\), not merely at \(x^{12}\).  The signs \(-18\) and \(+14\) were checked directly on the PDF page and are **not uncertain**.  The same section notes that the claimed magnetization exponent \(3/8\) remains \(3/8\) after the formula is reduced to two dimensions, contradicting Yang's exact \(1/8\).

### Cheapest decisive test

`test_id = "zhang2007_ht_k2"`.  Let

\[
R(K)=\frac{\exp(\Phi(K))}{2\cosh^3 K},\qquad \kappa=\tanh K.
\]

Extract \([\kappa^2]R\) either symbolically from Eq. (A.13) or from a few high-precision values near \(K=0\).  The finite-temperature claim predicts \(1/2\); exact graph enumeration predicts \(0\) (the first correction is \(3\kappa^4\)).  One coefficient, with no critical extrapolation or large lattice, is decisive.

## 2. Degang Zhang (2021): operator-algebra/screw-boundary claim

**Primary claim:** [arXiv:2110.11233](https://arxiv.org/pdf/2110.11233), *Exact Solution for Three-Dimensional Ising Model*.

The paper uses

\[
H_1=\frac{J_1}{k_BT},\qquad H_2=\frac{J_2}{k_BT},\qquad H=\frac{J}{k_BT},
\qquad H^*=\frac12\log\coth H=\operatorname{atanh}(e^{-2H}),
\]

in and immediately below Eq. (2).

### Claimed anisotropic critical line

The claimed critical line is

\[
\boxed{\sinh(2H)\sinh(2H_1+2H_2)=1.}
\tag{2110.11233, Eq. 32}
\]

For the isotropic simple-cubic case \(H_1=H_2=H=K\), this is

\[
\sinh(2K_c)\sinh(4K_c)=1.
\]

The paper prints \(H_c=0.30468893\).  Solving Eq. (32) independently at high precision gives

\[
K_c=0.304688931718003115768401685584199\ldots,
\]

or **\(0.304688931718003\)** to 15 significant decimal digits.  It exceeds \(0.221654626\) by \(0.083034305718003\), about \(37.46\%\).

### Full claimed integral, exactly transcribed

The partition function per atom is \(\lambda_\infty\), and Eqs. (35)–(36) print

\[
\boxed{
\log\lambda_\infty
=\frac12\log[2\sinh(2H)]
 +\frac{1}{2\pi}\lim_{m\to\infty}\int_0^\pi\xi_m(\omega)\,d\omega
}
\tag{2110.11233, Eq. 35}
\]

with

\[
\boxed{
\begin{aligned}
\cosh\xi_m(\omega)=D(\omega)={}&
 \cosh(2H_1)\cosh(2H_2)\cosh(2H^*)\\
&-\sinh(2H_1)\cosh(2H_2)\sinh(2H^*)\cos\omega\\
&-\cosh(2H_1)\sinh(2H_2)\sinh(2H^*)\cos(m\omega)\\
&+\sinh(2H_1)\sinh(2H_2)\cosh(2H^*)\cos[(m-1)\omega].
\end{aligned}}
\tag{2110.11233, Eq. 36}
\]

The source actually prints **\(D(\omega)\)**, not \(D_m(\omega)\); `D_m` below is only a code alias making the explicit \(m\)-dependence visible.  Equation (37) then states

\[
F=-Nk_BT\log\lambda_\infty.
\tag{2110.11233, Eq. 37}
\]

A verbatim coding transcription is:

```python
def Hstar(H):
    return mp.atanh(mp.exp(-2*H))

def D_m(omega, m, H1, H2, H):  # source notation: D(omega)
    Hs = Hstar(H)
    return (
        mp.cosh(2*H1)*mp.cosh(2*H2)*mp.cosh(2*Hs)
        - mp.sinh(2*H1)*mp.cosh(2*H2)*mp.sinh(2*Hs)*mp.cos(omega)
        - mp.cosh(2*H1)*mp.sinh(2*H2)*mp.sinh(2*Hs)*mp.cos(m*omega)
        + mp.sinh(2*H1)*mp.sinh(2*H2)*mp.cosh(2*Hs)*mp.cos((m-1)*omega)
    )

def xi_m(omega, m, H1, H2, H):
    return mp.acosh(D_m(omega, m, H1, H2, H))

def log_lambda_m(m, H1, H2, H):
    return (mp.log(2*mp.sinh(2*H))/2
            + mp.quad(lambda omega: xi_m(omega, m, H1, H2, H),
                      [0, mp.pi])/(2*mp.pi))
# Eq. (35) requires the limit of log_lambda_m over integer m -> infinity.
```

The order printed in Eq. (35) is important: integrate for each \(m\), then take \(m\to\infty\).  The integrand contains the oscillatory \(\cos(m\omega)\) and \(\cos[(m-1)\omega]\) terms.

For \(J_1=J_2\), the paper itself supplies a directly evaluable post-limit reduction,

\[
\cosh\xi_\infty(\omega)=
\cosh(2H_1)\cosh(2H^*_{2D})
-\sinh(2H_1)\sinh(2H^*_{2D})\cos\omega,
\tag{2110.11233, Eq. 40}
\]

\[
H^*_{2D}=H^*-H_1.
\tag{2110.11233, Eq. 41}
\]

Thus the JSON's isotropic \(K\)-only expression for \(\log\lambda_\infty\) is directly evaluable as written:

```text
(lambda Hs: mp.log(2*mp.sinh(2*K))/2 + mp.quad(lambda w: mp.acosh(mp.cosh(2*K)*mp.cosh(2*(Hs-K)) - mp.sinh(2*K)*mp.sinh(2*(Hs-K))*mp.cos(w)), [0, mp.pi])/(2*mp.pi))(mp.atanh(mp.exp(-2*K)))
```

### Published objections

1. **The operator identification is false.**  Perk, [arXiv:2202.03136](https://arxiv.org/pdf/2202.03136), Eq. (2) on p. 2, reconstructs the identification used in Zhang's Eqs. (16)–(17):
   \[
   \mathcal A_{p,1}=\mathcal L^p_{1,2}=\sigma_p^z\sigma_{p+m}^z
   \equiv A_{p,1}=L^p_{1,2}=L_{p,p+m}
   =\sigma_p^z\sigma_{p+1}^x\cdots\sigma_{p+m-1}^x\sigma_{p+m}^z.
   \]
   The two operators commute and have the same set of eigenvalues, but those \(\pm1\) eigenvalues are assigned differently to their common eigenvectors.  Equality therefore does not follow.  This invalidly replaces
   \[
   H_y=\sum_{\tau=1}^{mn}\sigma_\tau^z\sigma_{\tau+m}^z
   \]
   by the nonlocal Jordan–Wigner string
   \[
   A_m=\sum_{\tau=1}^{mn}\sigma_\tau^z
   \left(\prod_{j=1}^{m-1}\sigma_{\tau+j}^x\right)\sigma_{\tau+m}^z,
   \]
   converting the interacting 3D problem into a free-fermion one.  Perk's Eqs. (6)–(7) give an explicit \(m=2\) Kronecker-product contradiction.
2. **Boundary conditions cannot repair the series.**  Perk's Eq. (1) applies the Peierls–Bogolyubov bound
   \[
   |F[\mathcal H_2]-F[\mathcal H_1]|\leq\|\mathcal H_2-\mathcal H_1\|.
   \]
   Moving the periodic bonds to screw bonds changes the free energy per site by at most \(2|J_1|/m\), which vanishes as \(l,m,n\to\infty\).  Therefore the different thermodynamic-limit series that Zhang labels screw and periodic cannot both be correct.
3. **The critical line violates rotational invariance.**  Perk rewrites Eq. (32) as his Eq. (3).  For permutations of \((J_1,J_2,J)=(0.5,1.0,1.5)\), the same physical anisotropic lattice gives \(\beta_c=0.3503982204\), \(0.3046889317\), or \(0.2937911957\), depending only on which axis is called \(J\).  A rotation cannot change the critical temperature.
4. **The claimed series fails first at \(k^4\).**  With \(k=\tanh H\), Zhang's Eq. (44) gives
   \[
   \lambda_\infty=2\cosh^3H\,[1-3k^4-62k^6-1036k^8-20838k^{10}-\cdots],
   \]
   while the same paper's Eq. (45) quotes the established cubic result
   \[
   \lambda_\infty^{p}=2\cosh^3H\,[1+3k^4+22k^6+192k^8+2046k^{10}+\cdots].
   \]
   The first nonzero coefficient is already wrong in sign.  Perk's boundary-condition bound removes the paper's proposed explanation for the discrepancy.

### Cheapest decisive test

`test_id = "dzhang2021_ht_k4"`.  Form

\[
R(K)=\frac{\lambda_\infty(K)}{2\cosh^3K},\qquad k=\tanh K,
\]

and extract \([k^4]R\).  The claim predicts \(-3\) (Eq. 44); exact plaquette enumeration predicts \(+3\) (Eq. 45).  This needs only the first nontrivial graph coefficient and avoids the oscillatory \(m\to\infty\) quadrature entirely.

## 3. Rosengren (1985): exact-critical-point-only conjecture

This is not a complete free-energy solution, so `free_energy_python` is `null`.  It is included because it is an explicit exact-\(K_c\) claim discussed in the surveyed source.

Zhang's arXiv manuscript quotes Rosengren's conjecture on p. 35 of [arXiv:0705.1045](https://arxiv.org/pdf/0705.1045) as

\[
v_c=\tanh K_c=(\sqrt5-2)\cos\frac{\pi}{8}=0.218098372\ldots,
\]

with \(K_c=0.22165863\ldots\).  The displayed equation is **unnumbered in the source**; p. 35 and reference [230] are the exact locator, and no equation number is invented here.  Direct evaluation gives

\[
K_c=\operatorname{atanh}\!\left[(\sqrt5-2)\cos\frac{\pi}{8}\right]
=0.221658637208698872081822726113508\ldots,
\]

or **\(0.221658637208699\)** to 15 significant decimal digits.

Directly evaluable expression:

```text
mp.atanh((mp.sqrt(5)-2)*mp.cos(mp.pi/8))
```

The same p. 35 passage reports Fisher's objection that the basis of the guess is obscure: a no-backstep weighted-walk argument suggested the factor \(\sqrt5-2\), while \(\cos(\pi/8)\) was selected to match then-current 1981–1984 numerical estimates.  Independently, the conjecture is \(0.000004011208699\) above the later \(0.221654626(5)\) value quoted by Perk in arXiv:2202.03136, a difference of about 18.1 ppm and vastly larger than the quoted uncertainty.

**Cheapest decisive test:** `test_id = "rosengren1985_kc"`.  Evaluate the closed form once at 50 decimal digits and compare it with the interval \([0.221654621,0.221654631]\) represented by \(0.221654626(5)\).  The intervals are disjoint; no free-energy integration is involved.

## Uncertainty ledger

* `zhang2007`: `uncertain: true` solely because the primary source prints incompatible \((2\pi)^3\) and \((2\pi)^4\) normalizations for the same four-fold isotropic integral and supplies mutually incompatible fitted and finite-temperature weights.  Both printed alternatives are recorded above.
* `dzhang2021`: `uncertain: false`.  Eqs. (32), (35), (36), (40), (41), (44), and (45) were read directly from the PDF.  The source notation is \(D(\omega)\), despite its explicit \(m\)-dependence.
* `rosengren1985_kc`: `uncertain: false` for the quoted formula and decimal; the source simply does not assign its displayed equation a number.
* The low-temperature signs in Wu et al. p. 4 are not guessed from OCR: the PDF visibly reads \(-18x^{12}\) for Zhang's expansion and \(+14x^{12}\) for the exact expansion.
