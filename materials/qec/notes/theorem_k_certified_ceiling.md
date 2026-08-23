# Theorem K — the reciprocal-pole logical isomorphism and certified BB distance bounds

**Status (2026-08-21):** proved in the repository's matrix convention and
machine-verified on 27 sourced odd-lattice BB instances. The underlying
principal-code logical-space isomorphism is **published** (Eberhardt–Steffan,
arXiv:2407.03973v1, Corollaries 2.11–2.12); do not claim it as new. A novelty
search did not find the explicit general-BB *minimum-pole-weight ceiling* or its
use as a solver-free exhaustive rejection oracle, but both are short corollaries
of the published isomorphism and should be presented as such, not oversold as a
deep new theorem.

## 0. Setting and the load-bearing bar

Let $\ell,m$ be odd, $G=\mathbb Z_\ell\times\mathbb Z_m$, and
$$R=\mathbb F_2[G]
  =\mathbb F_2[x,y]/(x^\ell-1,y^m-1).$$
For $a,b\in R$, the CSS bivariate-bicycle code is
$$H_X=[A\ B],\qquad H_Z=[B^{\mathsf T}\ A^{\mathsf T}],\qquad n=2\ell m,$$
where $A,B$ are the binary circulants of $a,b$. Put
$S_Z=\operatorname{rowspace}H_Z$ and let $\bar{\cdot}$ be the involution
$x\mapsto x^{-1},y\mapsto y^{-1}$.

The code stores polynomial coefficient vectors as **rows**. Consequently
```python
I = nullspace(HX.T)
```
is a *left*-annihilator coefficient ideal. A physical column vector in the
right kernel of $H_X$ is obtained only after applying the reciprocal map. This
is the convention that must never be elided:
$$I=\operatorname{Ann}_{\rm left}(a,b),\qquad J=\bar I
  =\operatorname{Ann}_{\rm right}(a,b).$$
The distinction is invisible when $I=\bar I$ (for example the published
$[[90,8,10]]$ code) and load-bearing otherwise. On the published
$(7,7)\ [[98,6,12]]$ code every basis vector of raw $I$ fails the equation
$H_X(u,0)^{\mathsf T}=0$, while every basis vector of $J=\bar I$ satisfies it.
`test_bar_convention_on_noninvariant_code` locks this regression.

## 1. Published rate and logical-space input

Since $|G|$ is odd, Maschke's theorem makes $R$ semisimple:
$$R\cong\prod_{\chi\in X}\mathbb F_\chi,$$
with $X$ the Frobenius orbits of characters. Multiplication is diagonal, so if
$Z=\{\chi:a(\chi)=b(\chi)=0\}$, then
$$I\cong\bigoplus_{\chi\in Z}\mathbb F_\chi,\qquad
  k=2\dim_{\mathbb F_2}I
   =2\sum_{\chi\in Z}[\mathbb F_\chi:\mathbb F_2].$$
This is published in equivalent forms by:

* Panteleev–Kalachev, arXiv:1904.02703, Proposition 1
  ($k=2\deg\gcd(a,b,x^\ell-1)$ in the cyclic case);
* Lin–Pryadko, arXiv:2306.16400, Eq. (47);
* Wang–Mueller, arXiv:2408.10001v4, Eq. (11), for coprime lattices via
  $\pi=xy$;
* Postema–Kokkelmans, arXiv:2502.17052v4, Theorem 2.6;
* Eberhardt–Steffan, arXiv:2407.03973v1, Corollaries 2.11–2.12
  (if $\ell,m$ are odd, all BB codes are principal).

EXP-055 reproduces $k$ independently by annihilator dimension, direct matrix
ranks, and the published gcd formula where $\gcd(\ell,m)=1$.

## 2. The exact reciprocal-pole isomorphism

**Theorem K (matrix-convention form of the published principal-code
isomorphism).** Define
$$\mathcal P=J\oplus J\subseteq\mathbb F_2^{2\ell m},
  \qquad J=\bar I.$$
Then
$$\mathcal P\subseteq\ker H_X,\qquad
  \mathcal P\cap S_Z=0,\qquad
  \dim\mathcal P=k.$$
Therefore the quotient map restricts to an isomorphism
$$\boxed{\mathcal P\;\cong\;\ker H_X/S_Z.}$$
Every logical class has a unique reciprocal-pole representative $(u,v)$ with
$u,v\in J$.

**Proof.** In the character decomposition, $I$ is supported on $Z$ and
$J=\bar I$ on $Z^{-1}$. Physical column multiplication by $H_X$ evaluates the
reciprocal character, so both blocks of every $(u,v)\in J\oplus J$ are killed:
$\mathcal P\subseteq\ker H_X$. Meanwhile
$$S_Z=\{(\lambda\bar b,\lambda\bar a):\lambda\in R\}.$$
At a character $\chi\in Z^{-1}$,
$(\bar b(\chi),\bar a(\chi))=(b(\chi^{-1}),a(\chi^{-1}))=(0,0)$; hence $S_Z$
has zero component on every block that supports $\mathcal P$, and
$\mathcal P\cap S_Z=0$. Finally
$\dim\mathcal P=2\dim J=2\dim I=k=\dim(\ker H_X/S_Z)$ by the published rate
law. The injective quotient map is therefore surjective. $\square$

This proof also states precisely why using raw $I$ was wrong: $I$ is supported
on $Z$, while the physical pole and the zero blocks of $S_Z$ are supported on
$Z^{-1}$.

## 3. Exact variational formula and solver-free bounds

**Corollary K1 (exact pole-coset formula).**
$$\boxed{
 d(P)=\min_{0\ne p\in\mathcal P}\ \min_{s\in S_Z}\operatorname{wt}(p+s).
}$$
This is just the definition of distance after replacing arbitrary logical
representatives by the exact pole transversal $\mathcal P$.

**Corollary K2 (pole ceiling).** Taking $s=0$ and one pole block zero gives
$$\boxed{
 d(P)\le d_{\rm pole}:=
 \min\{\operatorname{wt}(u):0\ne u\in J\}.
}$$
Bar preserves Hamming weight, so the number equals
$\min\operatorname{wt}(I\setminus\{0\})$, but the *physical witness* is in
$J=\bar I$. For $k\le32$, exhaustive enumeration costs
$O(2^{k/2}\ell m)$ bit operations — microseconds, no solver, decoder, or
sampling.

**Corollary K3 (certified coset reduction).** Pick any nonzero pole
$p\in\mathcal P$ and any stabilizer $s\in S_Z$. Since
$\mathcal P\cap S_Z=0$, $p+s$ is automatically a nontrivial logical. Row
reducing selected pole vectors modulo $S_Z$, under several information-set
orders, therefore returns explicit self-certifying logical witnesses. Every
returned weight is a rigorous upper bound; randomisation changes only tightness,
never soundness.

This second route is much tighter in practice. On the published
$[[90,8,10]]$, the raw pole ceiling is $20$ while the reduced witness has weight
exactly $10$. In the Pareto screen, it turns many 30–90 second CP-SAT calls into
$\approx0.01$ s witness checks.

## 4. Validation and limitations

`experiments/exp055_odd_lattice_sweep.py literature` audits 27 sourced
odd-lattice instances:

* the reciprocal-pole isomorphism holds **27/27**: $\mathcal P\subseteq\ker H_X$,
  $\mathcal P\cap S_Z=0$, and $\dim\mathcal P=k$ by independent rank tests;
* 25 rows reproduce their printed $k$; two transcribed Wang–Mueller App. C rows
  ($(5,9)$ and $(7,11)$) give internal $k=0$ by both our matrix and gcd routes
  against printed 4 and 6, so they are excluded rather than reconciled;
* all 25 reproduced rows satisfy `pole ceiling >= reported d`, but this is only
  a **sanity check**: Wang–Mueller uses BP-OSD `distance_upperbound`, and
  Postema–Kokkelmans labels its table distances Monte-Carlo estimates;
* five rows are independently exact-certified here (Bravyi $[[90,8,10]]$,
  three small Postema rows, and Wang–Mueller's $[[126,12,10]]$ exactified by
  EXP-055): zero pole-ceiling violations, with slack min/median/max $2/10/22$.

The pole ceiling is sound but loose because setting $s=0$ ignores cancellations
inside a logical coset. The reduced-witness route repairs much of that looseness
but remains an upper-bound search, not an exact distance algorithm. Exactness
still requires exhaustive coset search, CP-SAT/SAT, matching bounds, or a proof.
No BP-OSD estimate is admitted as a domination threshold in EXP-055.

The semisimple character proof requires both $\ell,m$ odd. On even lattices the
ring has nilpotents and this pole decomposition is not valid without a local
Artinian refinement.

## 5. Retraction record

An intermediate EXP-055 draft used raw $I$ as a physical kernel and proposed a
member fallback when $I\cap\bar I=0$. That fallback was **wrong**; the new
regression on $(7,7)$ caught it before a final PDF or ledger entry was shipped.
The earlier restricted statement using $I\cap\bar I$ happened to be sound
because that subspace is bar-invariant, but it was unnecessarily weak. The
surviving result is the exact bar-aware isomorphism with $J=\bar I$ above.

Machine sources:

* `results/processed/exp055_literature_validation.json`
* `results/processed/exp055_odd_lattice_sweep.json`
* `results/processed/exp055_odd_lattice_screen.json`
* `results/certificates/exp055_odd_lattice_survivors.json`
* `tests/test_exp055_odd_lattice.py`
