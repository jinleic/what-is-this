# The infrared upper bound on `K_c` is exactly saturated

## 0. Result

For the nearest-neighbour ferromagnetic Ising model on `Z^3` with `K=beta*J`,
this note re-derives, re-certifies and then **exactly delimits** the incumbent
upper endpoint

\[
\boxed{K_c\le\frac{I_3}{2}
=0.2527310098586630030260020266135701299926 },
\qquad
I_3=\frac1{(2\pi)^3}\int_{[-\pi,\pi]^3}\frac{d^3k}{3-\sum_i\cos k_i},
\]

the decimal rounded upward.

**There is no strict improvement.** The wave-8 deliverable is the negative
result together with its certificate:

> **[THEOREM D, master barrier]** For every `K <= I_3/2` the explicit spectral
> profile `Ghat_*(k)=1/(I_3*lambda(k))`, `lambda(k)=3-sum_i cos k_i`, satisfies
> *every* constraint that reflection positivity supplies in the thermodynamic
> limit — spin normalisation `G(0)=1`, positive definiteness, GKS-I positivity
> `G(z)>=0`, and the pointwise infrared ceiling `Ghat<=1/(2K*lambda)` — and it
> has **zero** long-range order. Hence no consequence of that constraint system
> alone can prove long-range order at or below `I_3/2`, and `I_3/2` is the exact
> optimum of the method class.

Theorem D subsumes optimized test functions, block-spin averages, and weakened
anisotropic sub-models. Two further exact barriers (Theorems C and E) close the
same routes constructively, and Section 7 quantifies precisely where the witness
survives: on the whole current uncertainty interval
`[0.2122119011661678393310862783954278184914, 0.2527310098586630030260020266135701299926]`.

The benchmark `0.221654626` selected no test function, threshold, or rounding.
It appears only as a comparison in Section 7.

Artifacts: `experiments/e72_upper_infrared.py`,
`results/bounds/upper_infrared.json`, `tests/test_upper_infrared.py`.

## 1. Normalisation, and the two factors of two

**[EXTERNAL THEOREM, inherited with its audit trail]** In the ordered-pair
convention

\[
 H=-\sum_{x,y}j(y-x)\sigma_x\sigma_y,\qquad
 E(k)=\sum_z j(z)\bigl(1-\cos(k\cdot z)\bigr),
\]

Fröhlich–Simon–Spencer Gaussian domination gives, on an even torus `T_L` and for
every nonzero torus momentum,

\[
 \widehat G_L(k)\le\frac1{2\beta E(k)} .                     \tag{1}
\]

The source, the Ising-specific statement and the access audit are recorded in
`proofs/kc_bounds.md`, Section B2; the audited manifest
(`sources/manifest.yaml`, SHA-256
`f94175fc4fc536a74485ca6a9e6cc5947d35335b680a35fff380109d5c3cc024`) contains no
further reflection-positivity full text, so nothing stronger is imported here.

Each undirected bond occurs twice in the ordered sum, so reproducing
`-J*sum_{<x,y>} sigma_x sigma_y` forces

\[
 j(e_i)=j(-e_i)=\frac J2
 \quad\Longrightarrow\quad
 E(k)=J\lambda(k),\qquad \lambda(k)=3-\gamma(k),\qquad
 \gamma(k)=\sum_{i=1}^3\cos k_i ,
\]

and (1) becomes exactly

\[
 \widehat G_L(k)\le\frac1{2K\lambda(k)},\qquad k\ne0 .        \tag{2}
\]

Both factors of two are therefore pinned: the ordered-pair bond factor `1/2` and
the `2` of Gaussian domination. Using `j(\pm e_i)=J` instead would advertise the
false candidate `I_3/4`, which is why the normalisation is re-derived rather than
quoted.

Parseval with `sigma_x^2=1` gives the exact finite-volume sum rule

\[
 1=\underbrace{\frac{\widehat G_L(0)}{|T_L|}}_{M_L^2}
   +\frac1{|T_L|}\sum_{k\ne0}\widehat G_L(k) .                \tag{3}
\]

With (2), and because the `d=3` singularity is integrable so that even-torus
Riemann sums converge to `I_3`,

\[
 \liminf_{L\to\infty}M_L^2\ge1-\frac{I_3}{2K},
 \qquad\text{hence}\qquad K_c\le\frac{I_3}{2}.               \tag{4}
\]

## 2. Certified constants

With `W_sc=3I_3` and the audited Glasser–Zucker identity,

\[
 I_3=\frac{\sqrt6}{96\pi^3}
 \Gamma\!\left(\tfrac1{24}\right)\Gamma\!\left(\tfrac5{24}\right)
 \Gamma\!\left(\tfrac7{24}\right)\Gamma\!\left(\tfrac{11}{24}\right).
\]

`mpmath.iv` at `iv.dps=100`, outward rounded, gives

\[
\begin{split}
0.505462019717326006052004053227140259985129014817420892188993487886028773451173816800537247069896037925625
\le I_3\le{}\\
0.505462019717326006052004053227140259985129014817420892188993487886028773451173816800537247069896038011349,
\end{split}
\]

\[
\begin{split}
0.252731009858663003026002026613570129992564507408710446094496743943014386725586908400268623534948018962813
\le \tfrac{I_3}2\le{}\\
0.252731009858663003026002026613570129992564507408710446094496743943014386725586908400268623534948019005675 .
\end{split}
\]

The 100-digit enclosure (width `< 1e-100`) lies inside the incumbent 80-digit
enclosure and rounds upward to the same 40-place endpoint. Every inequality
between constants below is decided on these directed endpoints.

## 3. Two exact identities that control every refinement

### 3.1 Test functions and the block variance

Let `p` be a real polynomial with `p(3)=1`, and let `q:Z^3->R` be the finitely
supported real-space kernel with `qhat(k)=p(gamma(k))`. Because
`qhat(0)=p(3)=1`, Parseval applied to the nonnegative weight `p(gamma)^2` splits
the zero mode off exactly:

\[
 M_L^2
 =\sum_z (q\!*\!q)(z)\,G_L(z)
  -\frac1{|T_L|}\sum_{k\ne0}p(\gamma(k))^2\widehat G_L(k) .   \tag{5}
\]

The infrared bound (2) upper-bounds the second term by
`(2K)^{-1}<p^2/lambda>`. It says **nothing** about the first term, the block
variance. Every "optimized test function" proposal lives or dies on that first
term.

**[LEMMA 3.1, admissibility]** If `(q*q)(z)>=0` for all `z != 0`, then GKS-I
(`G(z)>=0`) gives

\[
 \sum_z (q\!*\!q)(z)G_L(z)\;\ge\;(q\!*\!q)(0)=\lVert q\rVert_2^2 .   \tag{6}
\]

This is the exact hypothesis under which the method closes; nonnegative kernels
`q>=0` are the special case used below.

### 3.2 The energy identity

`lambda = 3-gamma` has real-space kernel `3*delta_0 - (1/2)*sum_{\pm e_i}
delta_{\pm e_i}`, so Parseval gives the exact torus identity

\[
 \frac1{|T_L|}\sum_{k\ne0}\lambda(k)\widehat G_L(k)=3\bigl(1-G_L(e)\bigr),
                                                              \tag{7}
\]

`G_L(e)` being the nearest-neighbour correlation (the `k=0` term drops because
`lambda(0)=0`). Combining (7) with (2) yields the strongest nearest-neighbour
input the method itself can produce:

\[
 G_L(e)\ \ge\ 1-\frac{1}{6K}\Bigl(1-\frac1{|T_L|}\Bigr)
 \ \xrightarrow[L\to\infty]{}\ 1-\frac1{6K}.                  \tag{8}
\]

**[COMPUTATION]** Identity (7) is verified in exact rational arithmetic on the
eight-site `L=2` Ising torus at `w=e^{-2K}=3/5` (equivalently `v=tanh K=1/4`):
both sides equal `33939/15844` with `G(e)=4531/15844`. The experiment computes
the left side from a Boltzmann enumeration and the standalone test recomputes the
right side from an independent even-subgraph high-temperature enumeration.

**[COMPUTATION, hypothesis audit]** On that same degenerate torus the ceiling
(2) is *false*: at `w=3/5` the exact mode value `280125/269348 = 1.0400...`
exceeds `1/(2K\lambda)=0.9788...` for `lambda=2`, and all seven nonzero modes
violate it. The even-torus hypothesis behind (1) is load bearing and cannot be
dropped when shrinking the volume.

## 4. Exact barriers for test functions (Theorem C)

Let `C` denote the lattice Green function of `lambda`,

\[
 C(z)=\frac1{(2\pi)^3}\int_{[-\pi,\pi]^3}\frac{e^{ik\cdot z}}{\lambda(k)}d^3k,
 \qquad C(0)=I_3 .
\]

**[LEMMA 4.1, positivity of `C`]** For `0<=r<1`,
`1/(3-r\gamma)=\frac13\sum_{n\ge0}(r\gamma/3)^n` converges uniformly, and the
Fourier coefficient of `gamma^n` at `z` is `2^{-n}N_n(z)`, `N_n(z)` counting
`n`-step `\pm e_i` walks `0 -> z`. Hence
`C_r(z)=\frac13\sum_n r^n P_n(0\to z)` with
`P_n=N_n/6^n`, all terms nonnegative. Letting `r\uparrow1` with dominated
convergence (`1/(3-r\gamma)\le1/\lambda\in L^1`) gives

\[
 C(z)=\frac13\sum_{n\ge0}P_n(0\to z)\ \ge\ 0,                 \tag{9}
\]

with strict positivity for every `z` because some `n` admits a walk. Truncating
(9) gives **exact rational lower bounds**; at 30 steps, for instance,

\[
 C(0)\ge\frac{238615648580500269379}{511745184538734624768},
 \qquad
 C(e_1)\ge\frac{68033920400922061123}{511745184538734624768}.
\]

Independently, `sum_n P_n(0\to0)=3C(0)=W_sc`; the exact rational partial sums
`27413/20736`, `139555725223369/101559956668416`,
`238615648580500269379/170581728179578208256` increase and stay below the
certified lower endpoint of `W_sc`, and the last one is reproduced in the
standalone test by the multinomial formula
`N_{2n}(0)=\sum_{i+j+k=n}(2n)!/(i!^2j!^2k!^2)`.

**[THEOREM C, general test-function barrier]** Let `q` be finitely supported,
real, `q\ne0`, with `(q*q)(z)>=0` for `z\ne0` — i.e. exactly the hypothesis of
Lemma 3.1. Then

\[
 \Bigl\langle \frac{p(\gamma)^2}{\lambda}\Bigr\rangle
 =\sum_z (q\!*\!q)(z)\,C(z)
 \;\ge\;(q\!*\!q)(0)\,C(0)
 =I_3\lVert q\rVert_2^2 ,                                     \tag{10}
\]

so by (5)–(6) the resulting threshold is

\[
 R(q)=\frac{\langle p(\gamma)^2/\lambda\rangle}
           {2\lVert q\rVert_2^2}\;\ge\;\frac{I_3}2 ,
\]

with equality iff `(q*q)(z)=0` for every `z\ne0`. For finitely supported real
`q` the latter means `|qhat|^2` is constant, i.e. `qhat` is a monomial times a
constant, i.e. `q` is supported on a single site: the incumbent constant test.
**The very positivity that licenses the block-variance lower bound forces the
threshold back up to `I_3/2`.**

*Quantitative form.* If two distinct sites carry mass, say `q(z)q(w)>0` with
`z\ne w`, then

\[
 R(q)-\frac{I_3}2\ \ge\ \frac{q(z)q(w)\,C(z-w)}{\lVert q\rVert_2^2}
 \ \ge\ \frac{q(z)q(w)}{\lVert q\rVert_2^2}\cdot
 \frac13\sum_{n\le N}P_n(0\to z-w)\ >\ 0 ,
\]

an exact rational lower bound. For `q=\delta_0+\delta_{e_1}` the excess equals
`C(e_1)\ge68033920400922061123/511745184538734624768>0`.

### 4.1 The one-parameter family, in closed form

For `p_t(\gamma)=1+t\gamma`, `t>=0` (kernel `1` at the origin and `t/2` at each
neighbour, so `q>=0`), the moments
`<gamma^0,...,gamma^4>=1,0,3/2,0,45/8` and the recursion
`<gamma^n/\lambda>=3<gamma^{n-1}/\lambda>-<gamma^{n-1}>` give
`||q||_2^2=1+\tfrac32t^2` and

\[
 R(t)=\frac{I_3+2t(3I_3-1)+t^2(9I_3-3)}{2\left(1+\tfrac32t^2\right)},
 \qquad
 R(t)-\frac{I_3}2
 =\frac{t\bigl[(12I_3-4)+(15I_3-6)t\bigr]}{2\,(2+3t^2)} .      \tag{11}
\]

The certified interval gives `12I_3-4>2.0655` and `15I_3-6>1.5819`, both
increasing in `I_3`, so (11) is strictly positive for `t>0`: the optimum is
`t=0` with value `I_3/2`.

### 4.2 The energy-augmented family: exact first-order cancellation

Refusing to discard the off-diagonal terms, one may keep the six nearest
neighbours of `q*q` and lower-bound them by (8) with `g=1-1/(6K)`; all remaining
off-diagonal coefficients stay nonnegative and are dropped. The order criterion
becomes `D(t)>0` with

\[
 D(t)=2K\Bigl(1+\tfrac32t^2+6tg\Bigr)-\bigl[I_3+2t(3I_3-1)+t^2(9I_3-3)\bigr],
\]

and the exact polynomial identity (verified symbolically and by exact rational
sampling) is

\[
 \boxed{\,D(t)=(2K-I_3)(1+6t)+3t^2\,(K+1-3I_3)\,}.             \tag{12}
\]

Since `K+1-3I_3<0` for every `K\le I_3/2` (certified: `I_3>2/5`), we get
`D(t)\le0` throughout `K\le I_3/2`, with `D\equiv0` only in the limiting case
`t=0, 2K=I_3`. The linear-in-`t` gain from the nearest-neighbour energy is
`12Kg-2(3I_3-1)=6(2K-I_3)`, which **vanishes identically at the threshold**.
The saturation value in (8) at `K=I_3/2` is

\[
 1-\frac1{6K}\Big|_{K=I_3/2}=1-\frac1{W_{\rm sc}}
 =0.3405373295509991428262731844329\ldots ,
\]

exactly the nearest-neighbour value of the master witness of Section 5 (and,
by (9), Pólya's three-dimensional return probability `1-1/\sum_nP_n(0\to0)`).
The cancellation in (12) is therefore not a numerical accident: it is the
witness saturating the constraint it is tested against.

### 4.3 What the discarded term costs

Pointwise minimisation of only the infrared integral suggests smaller numbers:
`p=\gamma/3` gives `<p^2/\lambda>=I_3-1/3` and `(3I_3-1)/6=0.08606434...`, and
`p=(2\gamma^2-3)/15` gives `I_3-2/5` and `(5I_3-2)/10=0.05273100...`. Both are
recorded in the artifact as **not bounds**: their kernels have mixed signs, so
(6) is unavailable and the block variance in (5) is simply missing. Their
absurdity relative to any known behaviour is a useful alarm, but the rigorous
reason for rejection is the missing term, not the size of the number.

## 5. Master barrier (Theorem D)

Define the constraint system `S(K)` as everything the route of Sections 1–3
provides in the thermodynamic limit about a translation-invariant correlation
function `G` with spectral measure `mu = M^2 delta_0 + h(k)\,d^3k/(2\pi)^3`:

1. `G(0)=1` (spin normalisation / Parseval sum rule);
2. `mu>=0`, i.e. `G` positive definite (the correlation matrix is a Gram matrix);
3. `G(z)>=0` for all `z` (GKS-I);
4. `0<=h(k)<=1/(2K\lambda(k))` for a.e. `k` (the infrared ceiling (2)).

**[THEOREM D]** For every `K` with `2K<=I_3`, the profile

\[
 h_*(k)=\frac1{I_3\lambda(k)},\qquad M_*^2=0,
 \qquad\text{i.e.}\qquad
 G_*(z)=\frac{C(z)}{I_3},
\]

belongs to `S(K)` and has zero long-range order.

*Proof.* (1) `G_*(0)=C(0)/I_3=1`. (2) `h_*>=0`, and there is no negative part,
so `G_*` is positive definite; since the total mass is `G_*(0)=1`, `mu_*` is a
probability measure and `|G_*|\le1`. (3) `G_*(z)=C(z)/I_3\ge0` by Lemma 4.1.
(4) `h_*(k)=1/(I_3\lambda)\le1/(2K\lambda)` exactly when `2K\le I_3`.
Finally `M_*^2=0`. `QED`

**[COROLLARY]** No consequence of `S(K)` alone can prove `M^2>0` for any
`K\le I_3/2`. In particular `I_3/2` is the exact optimum of the method class:
every re-weighted sum rule, optimized test function, block-spin average, or
momentum-space rearrangement applied to the isotropic model is a consequence of
items 1–4 and therefore cannot certify a smaller threshold. Theorem C is the
constructive shadow of this fact, and the sub-model route, which replaces `G` by
the correlation function of a different Hamiltonian, is closed separately in
Section 6. Relaxing item 4 to permit further singular components away from the
zero mode only enlarges the feasible set, so it cannot break the barrier either.
The one case not covered is a genuinely finite-volume argument, because at finite
`L` the witness is infeasible; that case is analysed and left open in Section 7.

Two properties explain why the witness is so hard to exclude. Its
nearest-neighbour value saturates (8) exactly, and its susceptibility is
infinite:

\[
 \sum_z G_*(z)=\frac1{I_3}\sum_zC(z)
 =\frac1{3I_3}\sum_{n\ge0}\sum_zP_n(0\to z)=\infty
\]

by Tonelli on (9). The witness is a *critical-looking* profile with no order:
Gaussian domination cannot separate it from an ordered state below the point
where the sum rule saturates.

**Scope.** Theorem D is a barrier for `S(K)`, not for the model. Adding genuinely
new exact input — random-current representations, Simon–Lieb separator
inequalities, higher-order correlation identities, or a quantified uniform
finite-volume magnetisation floor — changes the feasible set and is not covered.
Nothing here suggests the true `K_c` is close to `I_3/2`; the benchmark sits far
below.

## 6. Anisotropic sub-models (Theorem E)

The only way to certify order in the isotropic model by studying another model is
to *weaken* couplings, because order in a sub-model implies order in the full one
(GKS-II monotonicity, as recorded in `proofs/kc_bounds.md`, Section B1).
Strengthening couplings bounds a different model.

**[THEOREM E]** Let `0<\alpha_i\le1` and give direction `i` the coupling
`\alpha_iK`. The infrared criterion for that reflection-positive sub-model reads
`K>T(\alpha)/2` with

\[
 T(\alpha)=\frac1{(2\pi)^3}\int_{[-\pi,\pi]^3}
 \frac{d^3k}{\sum_i\alpha_i(1-\cos k_i)} .
\]

Each integrand is nonincreasing in every `\alpha_i` and strictly decreasing on a
positive-measure set, so `T(\alpha)\ge T(1,1,1)=I_3` with equality iff
`\alpha=(1,1,1)`. Hence the best threshold available from any weakened
anisotropic sub-model is exactly `I_3/2`. `QED`

**[COMPUTATION]** Writing `A=\sum_i\alpha_i` and expanding as in Lemma 4.1,
`T(\alpha)=A^{-1}\sum_nR_n(\alpha)` with `R_n` the exact return probability of
the walk whose `\pm e_i` step has probability `\alpha_i/(2A)`. For
`\alpha=(1,1,\tfrac12)`, `(1,\tfrac12,\tfrac12)` and `(1,1,\tfrac1{10})` the
experiment verifies **termwise** `R_n(\alpha)/A\ge R_n(1,1,1)/3` for all
`n\le20` in exact rational arithmetic, and the partial sums are
`0.5636...`, `0.7100...`, `0.7725...` against the isotropic `0.4580...`. The
anisotropic route is therefore closed by proof, not by absence of a source.

## 7. Finite volume: exact data, rejected numbers, and the open route

Let `C_L` be the zero-mode-removed torus Green function,
`C_L(z)=|T_L|^{-1}\sum_{k\ne0}e^{ikz}/\lambda(k)`. For
`L\in\{2,3,4,6\}` all momentum cosines are rational, so `C_L` is exactly
rational:

| `L` | `C_L(0)` | `min_z C_L(z)` | argmin | `C_L(0)/2` |
|---:|---|---|---|---|
| 2 | `29/96` | `-11/96` | `(1,1,1)` | `29/192` |
| 3 | `88/243` | `-11/243` | `(1,1,1)` | `44/243` |
| 4 | `1517/3840` | `-49/1280` | `(2,2,2)` | `1517/7680` |
| 6 | `1289503/2993760` | `-13861/598752` | `(3,3,3)` | `1289503/5987520` |

Each row satisfies `\sum_z C_L(z)=0` and the defect equation
`3C_L(z)-\tfrac12\sum_{\pm e}C_L(z\pm e)=\delta_{z,0}-|T_L|^{-1}` exactly; the
standalone test verifies the defect equation at **every** site of every listed
torus, which pins the normalisation independently of the Fourier code.

Two honest consequences.

* **[REJECTED]** The numbers `C_L(0)/2` are all below the incumbent endpoint.
  They are **not** `K_c` bounds. The magnetisation criterion needs a positive
  limit inferior of `M_L^2`, and `M_L^2>0` holds trivially at finite `L`; a
  single torus certifies nothing about `Z^3`.
* **[OPEN ROUTE]** Because `min_z C_L<0`, the saturating witness of Theorem D is
  *infeasible* at finite even `L`: the finite-volume GKS constraint is not
  vacuous, so Theorem D does not literally cover finite-volume arguments. But
  `C_L(0)` increases toward `I_3` and `|min_z C_L|` decreases
  (`0.1146, 0.0453, 0.0383, 0.0231`), so the extra information vanishes with
  `L`. An improvement along this line therefore requires a **quantified
  uniform-in-`L` magnetisation floor**, which is absent from the audited
  manifest. That is the identified open route, and it is stated as such rather
  than claimed.

**Barrier window.** Adding susceptibility finiteness excludes the witness only
where finiteness is actually proved. Elementarily `c_n\le6\cdot5^{n-1}`, so the
self-avoiding-path majorant `chi\le\sum_nc_nv^n` is finite for `v<1/5`, i.e. for

\[
 K<\operatorname{atanh}(1/5)=0.2027325540540821909890065577321745682859\ldots ,
\]

and the repository's certified lower endpoint pushes the same conclusion to
`0.2122119011661678393310862783954278184914`. Hence the witness survives on

\[
 [\,0.2122119011661678393310862783954278184914,\;
   0.2527310098586630030260020266135701299926\,],
\]

which is exactly the current certified uncertainty interval. Any improvement of
the upper endpoint must supply new exact input inside that window. The benchmark
`0.221654626` lies inside it; it is quoted for comparison only.

## 8. Classification and reproduction

| item | status |
|---|---|
| infrared inequality (1) and its magnetisation consequence (4) | rigorous theorem, inherited with the audit trail of `proofs/kc_bounds.md` |
| ordered-pair normalisation and both factors of two | rigorous, re-derived here |
| Watson gamma identity | rigorous published identity, inherited |
| 100-dps interval for `I_3`, `I_3/2`, `1-1/W_sc`, `atanh(1/5)` | computer-assisted rigorous, outward directed rounding |
| positivity `C(z)>=0` and its rational truncations | rigorous theorem plus exact rational certificates |
| Theorem C, its equality case, and identities (11), (12) | rigorous theorems; identities verified symbolically and by exact rational sampling |
| Theorem D (master barrier) | rigorous theorem for the stated constraint system `S(K)` |
| Theorem E (anisotropic sub-models) | rigorous theorem; termwise exact verification |
| exact torus table and defect equation | exact rational computation |
| `C_L(0)/2` values, pointwise-optimized numbers `0.08606...`, `0.05273...` | explicitly rejected, not bounds |
| uniform-in-`L` finite-volume improvement | unresolved, open route |
| `0.221654626` | comparison only; no logical role |

Run from the repository root:

```text
.venv/bin/python experiments/e72_upper_infrared.py
.venv/bin/python tests/test_upper_infrared.py
```

The artifact `results/bounds/upper_infrared.json` uses the
`provenance / data / checks` envelope, records the SHA-256 of every audited local
input, and both commands print a terminal `PASS` only when all checks pass.
