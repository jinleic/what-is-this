# Theorem K — a certified solver-free distance ceiling for BB codes on odd lattices

**Status (2026-08-21):** proved; validated against 25 reproduced published
instances with zero violations (EXP-055). Novelty check (read-only literature
scout, 2026-08-21): this ceiling was **NOT FOUND** in the literature in this
form; the nearest published items are Eberhardt–Steffan's explicit logical bases
for odd lattices (arXiv:2407.03973v1, Cor. 2.11–2.12) and Rabeti–Mahdavifar's
solver-free ceilings for the Frobenius-restricted UB subclass
(arXiv:2605.14173v1, Cor. 4/6). Everything in §1 is *published* and is reproduced
here only to fix notation.

## 0. Setting

$\ell, m$ odd, $G=\mathbb Z_\ell\times\mathbb Z_m$, $R=\mathbb F_2[G]
=\mathbb F_2[x,y]/(x^\ell-1,\,y^m-1)$. For $a,b\in R$ the CSS bivariate-bicycle
code is
$$H_X=[A\;B],\qquad H_Z=[B^{\mathsf T}\;A^{\mathsf T}],\qquad n=2\ell m,$$
with $A,B$ the circulant matrices of $a,b$. Write $S_Z=\operatorname{rowspace}
H_Z$, and let $\bar{\cdot}$ be the involution $x\mapsto x^{-1},y\mapsto y^{-1}$,
so that $A^{\mathsf T}$ is the circulant of $\bar a$.

Because $|G|$ is odd, $\gcd(2,|G|)=1$ and Maschke's theorem makes $R$
**semisimple**:
$$R\;\cong\;\prod_{\chi\in X}\mathbb F_\chi ,$$
where $X$ indexes the Frobenius orbits of characters $G\to\bar{\mathbb F}_2^{\,*}$
and $\dim_{\mathbb F_2}\mathbb F_\chi$ is the orbit length. Multiplication by $a$
is diagonal in this decomposition, acting on the $\chi$-block as $a(\chi)$. Put
$$Z(a)=\{\chi: a(\chi)=0\},\qquad Z=Z(a)\cap Z(b),\qquad
I=\operatorname{Ann}(a)\cap\operatorname{Ann}(b),$$
and $Z^{-1}=\{\chi^{-1}:\chi\in Z\}$. All sets are Frobenius-closed;
$|\cdot|$ counts characters, i.e. $\mathbb F_2$-dimensions.

## 1. Published input (reproduced, not claimed)

**Lemma 1.** $\operatorname{Ann}(a)=\bigoplus_{\chi\in Z(a)}\mathbb F_\chi$, hence
$I=\bigoplus_{\chi\in Z}\mathbb F_\chi$ and
$$k \;=\; 2\dim_{\mathbb F_2} I \;=\; 2|Z| .$$

*Proof.* Multiplication by $a$ is diagonal, so its kernel is the sum of the
blocks where $a(\chi)=0$. $\square$

This is the odd-lattice rate law. It is published in several equivalent forms:
Panteleev–Kalachev arXiv:1904.02703 Prop. 1 ($k=2\deg\gcd(a,b,x^\ell-1)$, cyclic
case), Lin–Pryadko arXiv:2306.16400 Eq. (47), Wang–Mueller arXiv:2408.10001v4
Eq. (11) (coprime case, via $\pi=xy$), Postema–Kokkelmans arXiv:2502.17052v4
Thm. 2.6, and Eberhardt–Steffan arXiv:2407.03973v1 Cor. 2.11–2.12 ("if $\ell$ and
$m$ are odd, all BB codes are principal"). EXP-055 reproduces it by three
independent routes on 25 published instances.

## 2. Where the trivial operators live

**Lemma 2.** $I\times I\subseteq\ker H_X$, and
$$(I\times I)\cap S_Z=\bigoplus_{\chi\in Z\setminus Z^{-1}}(S_Z)_\chi ,
\qquad \dim\big((I\times I)\cap S_Z\big)=|Z\setminus Z^{-1}| .$$

*Proof.* For $u,v\in I$, $H_X(u,v)^{\mathsf T}=au+bv=0$, giving the inclusion.
Next, $S_Z=\{(\bar b\lambda,\ \bar a\lambda):\lambda\in R\}$, so at a character
$$(S_Z)_\chi=\operatorname{span}_{\mathbb F_\chi}\big\{(\,b(\chi^{-1}),\,
a(\chi^{-1})\,)\big\},$$
using $\bar a(\chi)=a(\chi^{-1})$. Hence $(S_Z)_\chi=0$ exactly when
$\chi^{-1}\in Z$, and is a line otherwise. Since $(I\times I)_\chi=
\mathbb F_\chi^{\,2}$ for $\chi\in Z$ and $0$ otherwise, and a line inside
$\mathbb F_\chi^2$ is contained in it, the intersection is the sum of $(S_Z)_\chi$
over $\chi\in Z$ with $\chi^{-1}\notin Z$. $\square$

Two consequences worth stating separately. First, $\dim S_Z=\ell m-|Z^{-1}|
=\ell m-|Z|$ and $\dim\ker H_X=\ell m+|Z|$, which re-derives $k=2|Z|$
independently of Lemma 1. Second, $I\times I$ spans the logical space **iff**
$Z$ is inversion-closed; otherwise it covers a subspace of dimension
$2|Z|-|Z\setminus Z^{-1}|$.

## 3. The ceiling

**Theorem K.** Let $I_0=I\cap\bar I$. Then
$I_0=\bigoplus_{\chi\in Z\cap Z^{-1}}\mathbb F_\chi$,
$(I_0\times I_0)\cap S_Z=0$, and consequently **every** nonzero $u\in I_0$ makes
$(u,0)$ a nontrivial logical operator. Hence
$$\boxed{\;d(P)\;\le\;\min\{\operatorname{wt}(u)\;:\;0\neq u\in I_0\}\;}$$
whenever $I_0\neq 0$.

*Proof.* $\bar I=\bigoplus_{\chi\in Z^{-1}}\mathbb F_\chi$ because the involution
permutes blocks by $\chi\mapsto\chi^{-1}$; intersecting gives the stated $I_0$,
whose index set $Z\cap Z^{-1}$ is inversion-closed. Apply Lemma 2 to the pair
$(I_0,I_0)$: its trivial part is supported on
$(Z\cap Z^{-1})\setminus(Z\cap Z^{-1})^{-1}=\varnothing$. So no nonzero element
of $I_0\times I_0$ lies in $S_Z$; in particular $(u,0)\notin S_Z$ while
$(u,0)\in\ker H_X$, i.e. $(u,0)$ is a nontrivial logical of weight
$\operatorname{wt}(u)$. $\square$

**Cost.** $O(2^{\dim I_0}\cdot \ell m)$ bit operations — microseconds for the
$k\le 32$ regime, since $\dim I_0\le\dim I=k/2$. No solver, no decoder, no
sampling. This is the only distance bound in the programme that costs less than
a rank computation.

**Corollary K1 (member-level fallback).** If $Z\cap Z^{-1}=\varnothing$ the
theorem is vacuous, but for a *fixed* pair the smallest weight of a nontrivial
$(u,0)$, $u\in I$, is still an upper bound on $d$, decidable by one rank test per
candidate $u$. EXP-055 uses this whenever $I_0=0$; on the 27-instance literature
battery it closes all 11 gaps left by Theorem K.

**Scope, stated honestly.** The character argument needs semisimplicity, i.e.
$\ell,m$ both odd; on even lattices $R$ has nilpotents and Lemma 2's blockwise
computation fails. Corollary K1 is convention-free and applies to any BB code,
but it is a per-instance rank test rather than a closed-form certificate.

## 4. What it is good for, and what it is not

The ceiling is an **upper** bound, so it certifies *rejection*: if
$\min\operatorname{wt}(I_0\setminus 0)\le D$ then the code provably cannot beat
distance $D$. In EXP-055's Pareto screen this rejects candidates with **zero**
solver calls, which is what makes an exhaustive odd-lattice screen affordable.

It is **not** tight. Measured slack against the 25 reproduced published
instances: minimum $2$, median $12$, maximum $42$ — roughly a factor $2$, the
same looseness recorded for the earlier idempotent-ideal bound in
`failed_routes.md` (FR-024). The reason is structural: $\ker H_X$ is strictly
larger than $I\times I$ whenever $Z$ is not inversion-closed, and the true
minimum is often carried by syzygy pairs $(u,v)$ with $\bar au=\bar bv\neq0$
rather than by the pole ideal. A tight algebraic bound would need the coset
minimum in $\ker H_X/S_Z$, which is exactly the open item D of
`next_breakthroughs.md`.

## 5. Machine record

* `experiments/exp055_odd_lattice_sweep.py literature` —
  `results/processed/exp055_literature_validation.json`: 27 sourced instances,
  25 reproduced, $k$ agreeing by three independent routes (annihilator, matrix
  rank, published gcd formula where the lattice is coprime), **zero ceiling
  violations**, slack $2$/$12$/$42$ (min/median/max).
* Two rows (both arXiv:2408.10001v4 App. C Table 4, $(5,9)$ and $(7,11)$) give
  $k=0$ by **both** of our independent routes against the printed $k=4$ and
  $k=6$; recorded as `unreproduced` and excluded from every count. We did not
  read that appendix ourselves, so a transcription error on our side is the
  likeliest explanation and is stated as such.
* `experiments/exp055_odd_lattice_sweep.py run` —
  `results/processed/exp055_odd_lattice_sweep.json`: all 65 odd lattices with
  $\ell m\le180$, $4.23\times10^9$ weight-$\le3$ pairs, $k$ by two routes with
  zero mismatches, 273 idempotence tests with zero violations (Theorem J-G1's
  mechanism re-verified on lattices outside our catalogue).
