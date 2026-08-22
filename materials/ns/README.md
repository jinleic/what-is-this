# Canonical target — Clay Navier–Stokes, Fefferman alternative (C)

> **Domain scope.** The obstruction map in this directory is for **(C), on
> $\mathbb{R}^3$**. It does **not** cover **(D)** on the torus $\mathbb{R}^3/\mathbb{Z}^3$.
> Every argument here — Lemma 1's Liouville step, the self-similar and DSS reductions,
> and Lemma 4's Leray-projection covariance — runs on the Euclidean dilation
> $x\mapsto\lambda x$ and the $\mathbb{R}^3$ Leray multiplier. That dilation is **not a
> self-map of the fixed torus** (it carries a $1$-periodic field to a $1/\lambda$-periodic
> one), so for $\lambda\ne1$ continuous self-similarity and DSS have no torus analogue as
> written. Liouville also changes: a periodic harmonic field is *constant*, not zero.
> The cited rigidity theorems (NRS, Tsai, Chae, Chae–Wolf, ESS) are all stated on
> $\mathbb{R}^3$. Transferring them to a torus singularity would plausibly go by rescaling
> around the singular point to a whole-space limit — **not attempted or claimed here.**
> Since (C) alone suffices for the prize, this is a scope statement, not a defect.

This directory holds the **actual** target of the project. Everything here is tied to
the official CMI problem statement. Nothing in `../ccf/` is on this path; see
[Scope discipline](#scope-discipline).

## The statement being attacked

Fefferman, *Existence and smoothness of the Navier–Stokes equation*
(<https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf>):

$$\partial_t u_i + \sum_j u_j\partial_{x_j}u_i = \nu\Delta u_i - \partial_{x_i}p + f_i,
\qquad \nabla\cdot u = 0,\qquad u(x,0)=u^{\circ}(x),\qquad \nu>0,\; n=3.$$

**(C)** There exist a smooth divergence-free $u^{\circ}$ on $\mathbb{R}^3$ obeying the
Schwartz-decay condition (4) and a smooth $f$ obeying (5), for which **no** pair $(p,u)$
satisfying (1),(2),(3),(6),(7) exists on $\mathbb{R}^3\times[0,\infty)$.
**(D)** The periodic analogue on $\mathbb{R}^3/\mathbb{Z}^3$.

Per CMI rules §5(b), *"a resolution in either direction will be evaluated by the
standard evaluation procedure"* — so a breakdown construction is fully prize-eligible.
See [`../docs/CMI_RULES.md`](../docs/CMI_RULES.md).

## Active object — the blowup-ansatz obstruction map

Before constructing anything, establish exactly which *shapes* a finite-time breakdown
for constant-viscosity NS is permitted to have. Each entry is either proved here or
cited with its exact hypotheses; the value of the map is that its holes are the only
places a construction can live.

### Lemma 0 (scaling-invariant ansätze). Immediate.

NS has the symmetry $u\mapsto \lambda u(\lambda^2 t,\lambda x)$. Writing $\tau=T-t$ and
demanding that $u=\tau^{-\alpha}U(x/\tau^{\beta})$ be *invariant* under it gives
$\lambda u(\lambda^2t,\lambda x) = \lambda^{1-2\alpha}\tau^{-\alpha}U\!\big(\lambda^{1-2\beta}x/\tau^{\beta}\big)$,
which equals $u$ for all $\lambda>0$ only if $\alpha=\beta=1/2$. So *scaling-invariant*
self-similarity is Leray's ansatz and nothing else. Lemma 1 addresses the broader
**first-kind** family, where $\alpha$ is not pinned by the symmetry.

### Lemma 1 (exponent rigidity, first-kind ansatz). Proved here.

*Let $\nu>0$, $\alpha>0$, $\alpha+\beta=1$, and suppose*
$$u(x,t) = \tau^{-\alpha}\,U(y),\qquad p(x,t)=\tau^{-2\alpha}P(y),\qquad
y=\frac{x}{\tau^{\beta}},\quad \tau=T-t,$$
*solves the unforced 3D Navier–Stokes equations exactly on $\mathbb{R}^3\times(t_0,T)$,
with $U$ smooth, divergence-free, and $U(y)\to0$ as $|y|\to\infty$. Then $\beta=1/2$ or
$U\equiv0$.*

**Proof.** With $\nabla_x=\tau^{-\beta}\nabla_y$ and $\Delta_x=\tau^{-2\beta}\Delta_y$,
$$\partial_t u = \tau^{-\alpha-1}\big[\alpha U + \beta\,y\cdot\nabla U\big],\qquad
(u\cdot\nabla)u = \tau^{-2\alpha-\beta}(U\cdot\nabla)U,$$
$$\nabla_x p = \tau^{-2\alpha-\beta}\nabla P,\qquad
\nu\Delta_x u = \nu\,\tau^{-\alpha-2\beta}\Delta U.$$
The hypotheses $\alpha+\beta=1$ and $\gamma=2\alpha$ make the three inviscid exponents
coincide: $\alpha+1 = 2\alpha+\beta = \gamma+\beta$. Dividing by $\tau^{-\alpha-1}$ and
using $\alpha-\beta = 1-2\beta$,
$$\alpha U + \beta\,y\cdot\nabla U + (U\cdot\nabla)U + \nabla P
\;=\; \nu\,\tau^{\,1-2\beta}\,\Delta U. \tag{$\star$}$$
The left side is $t$-independent. If $\beta\ne1/2$ then $1-2\beta\ne0$, so
$\tau\mapsto\tau^{1-2\beta}$ is non-constant; evaluating $(\star)$ at $t_1\ne t_2$ in
$(t_0,T)$ and subtracting,
$$\nu\big[\tau_1^{1-2\beta}-\tau_2^{1-2\beta}\big]\,\Delta U = 0 .$$
The bracket is nonzero and $\nu>0$, so $\Delta U\equiv0$; a harmonic function on
$\mathbb{R}^3$ vanishing at infinity is $0$ by Liouville. Hence $U\equiv0$, and then
$\nabla P\equiv0$. $\blacksquare$

> **The hypotheses are load-bearing, and $\alpha+\beta=1$ is assumed, not derived.**
> An earlier draft claimed "matching $\partial_t u$ against $(u\cdot\nabla)u$ forces
> $\alpha+\beta=1$". That is false as stated: powers of $\tau$ are linearly independent,
> so the equation only requires that the profile coefficients sharing each *distinct*
> exponent sum to zero — and terms may cancel in groups. For instance if
> $(U\cdot\nabla)U+\nabla P\equiv0$ (with $\gamma=2\alpha$) that pair forms its own
> group and imposes no exponent relation on it; the degenerate case $\alpha=\beta=0$ is
> stationary rather than blowup. A full classification of every exponent grouping is a
> separate statement and is **not attempted here**.
>
> $\alpha+\beta=1$ is the standard *first-kind* condition — exactly what makes the
> inertial terms balance the time derivative — and by Lemma 0 the scaling-invariant
> case forces the stronger $\alpha=\beta=1/2$ anyway. So the hypothesis costs little in
> practice, but it is a hypothesis and is recorded as one.

### Lemma 2 (the critical case $\beta=1/2$). Cited, hypotheses verified.

$\beta=1/2$ is exactly Leray's 1934 backward self-similar ansatz
$u=\lambda(t)U(\lambda x)$, $p=\lambda^2 P$, $\lambda=1/\sqrt{2a(T-t)}$, with profile
equation
$$-\nu\Delta U + aU + a\,(y\cdot\nabla)U + (U\cdot\nabla)U + \nabla P = 0,\qquad \nabla\cdot U=0.$$

| Result | Hypothesis | Conclusion |
|---|---|---|
| Nečas–Růžička–Šverák, Acta Math. **176** (1996) 283–294, [doi 10.1007/BF02551584](https://doi.org/10.1007/BF02551584), Thm 1 (`sources/NRS-1996-acta.pdf`) | weak profile solution ($W^{1,2}_{loc}$, div-free, distributional), $U\in L^3(\mathbb{R}^3)$ | $U=0$ |
| Tsai, ARMA **143** (1998) 29–51, Thm 1 ([PDF](https://personal.math.ubc.ca/~ttsai/publications/leray.pdf)) | weak profile $U\in L^q$, $q\in(3,\infty]$ | $U$ constant; $U=0$ for $q<\infty$ |
| Tsai, Thm 2 | **distributional** NS solution of self-similar form with a *local* energy estimate $\operatorname*{ess\,sup}\int_B|u|^2+\nu\iint_B|\nabla u|^2<\infty$ on some ball | $u=0$ |

Tsai's Theorem 2 is the operative one: it assumes only distributional NS plus a local
energy bound, which any physically meaningful solution satisfies. **Exact backward
self-similar blowup for constant-$\nu$ 3D Navier–Stokes is closed.**

### Lemma 3 (asymptotically self-similar). Corrected.

Chae, [arXiv:math/0604234](https://arxiv.org/abs/math/0604234), Math. Ann.
[doi 10.1007/s00208-007-0082-6](https://doi.org/10.1007/s00208-007-0082-6). Primary read.

**Thm 1.4** — $p\in[3,\infty)$, $v\in C([0,T);L^p(\mathbb{R}^3))$ classical. If there is
$\bar V\in L^p$ with
$\lim_{t\uparrow T}(T-t)^{\frac{p-3}{2p}}\lVert v(\cdot,t)-(T-t)^{-1/2}\bar V(\cdot/\sqrt{T-t})\rVert_{L^p}=0$,
then $\bar V=0$ and $v$ continues in $C([0,T+\delta);L^p)$.

**Thm 1.5** (localized) — same hypotheses, with either (i) $q\in[3,\infty)$ and finite
$R$, or (ii) $q\in[2,3)$ for every $R>0$, and
$\lim_{t\uparrow T}(T-t)^{\frac{q-3}{2q}}\sup_{t<\tau<T}\lVert v(\cdot,\tau)-(T-\tau)^{-1/2}\bar V((\cdot-z)/\sqrt{T-\tau})\rVert_{L^q(B(z,R\sqrt{T-t}))}=0$.
Then $\bar V=0$ and $v$ is Hölder near $(z,T)$. Remark 1.7 allows $q\in[2,3)$.

> **Correction.** An earlier revision said "Chae explicitly leaves weaker pointwise
> convergence open" **for Navier–Stokes**. Wrong attribution: that wording is Chae's
> Remark 1.4, which concerns **Euler** Theorem 1.2. Chae's NS statements carry no such
> reservation. For NS the coverage is global convergence $q\in[3,\infty)$ and localized
> $q\in[2,\infty)$.

**Do not confuse forward with backward.** Nontrivial *forward* self-similar and DSS
solutions genuinely exist — Jia–Šverák; Tsai; Bradshaw–Tsai
[arXiv:1610.01386](https://arxiv.org/abs/1610.01386); Chae–Wolf
[arXiv:1801.08060](https://arxiv.org/abs/1801.08060) in the local Leray class, from
$L^2_{loc}$ data. Forward self-similarity is evolution *away* from singular data and is
**not** a blowup mechanism. Citing it as evidence for blowup is a category error.

### Verification status

| Source | Retrieved | Status |
|---|---|---|
| Fefferman (CMI) | **primary, read here** | verified verbatim |
| Tsai, ARMA 143 (1998) | **primary** ([author PDF](https://personal.math.ubc.ca/~ttsai/publications/leray.pdf)) | verified verbatim |
| Chae, arXiv:math/0604234 | **primary** | verified verbatim |
| Chae–Wolf, arXiv:1610.09464 | **primary** (ar5iv) | verified verbatim |
| KNŠŠ, arXiv:0709.3599 | **primary** (ar5iv) | verified |
| Nečas–Růžička–Šverák (1996) | **primary** (Acta Math. backfile on [Project Euclid](https://projecteuclid.org/journals/acta-mathematica/volume-176/issue-2/On-Lerays-self-similar-solutions-of-the-Navier-Stokes-equations/10.1007/BF02551584.full); local `sources/NRS-1996-acta.pdf`) | verified verbatim 2026-08-15 |
| Escauriaza–Seregin–Šverák (2003) | **primary** (authors' IMA preprint 1904, [UMN Conservancy](https://conservancy.umn.edu/items/1ded2377-3cb2-43af-b8ff-f8d52ce5e3fe); local `sources/ESS-2003-ima1904.pdf`; published RMS **58**:2 (2003) 211–250, [doi 10.1070/RM2003v058n02ABEH000609](https://doi.org/10.1070/RM2003v058n02ABEH000609)) | verified verbatim 2026-08-15 |

Tsai Thm 1, verbatim: *"If a weak solution $U$ of (1.3) belongs to $L^q(\mathbb{R}^3)$,
for some $q\in(3,\infty]$, then it must be constant (and hence identically zero if
$q<\infty$)."*
Thm 2, verbatim: *"Suppose $u$ is a weak solution of (1.1) satisfying the local energy
estimates (1.4) in the cylinder $Q_1(0,T)=B_1(0)\times(T-1,T)$. If $u$ is of the form
$(1.2)_1$, then $u$ is identically zero."* His §2 states Thm 2 needs **only**
distributional NS plus the local energy estimate — no Leray–Hopf, suitable-weak or
pressure hypothesis. The project's "weak enough that nothing reasonable escapes" reading
**is supported**, but strictly for *exact, unforced* Leray self-similarity.

NRS Thm 1 (p. 291), verbatim: *"Let $U$ be a weak solution of (1.3) belonging to
$L^3(\mathbf{R}^3)$. Then $U\equiv0$ in $\mathbf{R}^3$."* Their weak solution (p. 286) is
$U\in W^{1,2}_{loc}$, divergence-free, distributional — no decay, pressure, or global
energy hypothesis; any $a>0$, $\nu>0$. Proof: decay estimates (Lemma 3.2,
$|\nabla^kU|=O(|y|^{-3-k})$) plus a maximum principle for
$\Pi=\frac12|U|^2+P+ay\cdot U$ (Lemma 3.3). NRS themselves flag (p. 284) that local
energy bounds (their (1.5)) do **not** imply $U\in L^3$, so their theorem alone leaves
locally-energy-finite self-similar singularities open — precisely the gap Tsai Thm 2
closes. The NRS→Tsai division of labor recorded here is the published one.

**ESS, verified at source.** Norm definition (preprint p. 4):
$\lVert f\rVert_{s,l,Q_T}=\big(\int_0^T\lVert f(\cdot,t)\rVert_s^l\,dt\big)^{1/l}$,
ess-sup in $t$ for $l=\infty$ — so the title's $L_{3,\infty}(Q_T)$ is the mixed class
$L^\infty_tL^3_x$, **not** Lorentz weak-$L^3$; the caveat below is confirmed at source.
Thm 1.3, verbatim: *"Assume that $v$ is a weak Leray-Hopf solution to the Cauchy problem
(1.1) and (1.2) in $Q_T$ and satisfies the additional condition (1.13)
[$v\in L_{3,\infty}(Q_T)$]. Then, $v\in L_5(Q_T)$ and hence it is smooth and unique in
$Q_T$."* Blow-up form (p. 4): if $]0,T_\star[$ is the maximal interval of existence of
the smooth solution and $T_\star<+\infty$, then
$\limsup_{t\uparrow T_\star}\int_{\mathbb{R}^3}|v(x,t)|^3\,dx=+\infty$ — the paper's own
display is a **limsup**, matching the precision note kept in `../PROGRESS.md`. Thm 1.4 is
the local version: on $Q=B\times]{-1},0[$ with $v\in L_{2,\infty}(Q)\cap
L_2(-1,0;W^1_2(B))$, $p\in L_{3/2}(Q)$ and $\lVert v\rVert_{3,\infty,Q}<\infty$, $v$ is
Hölder continuous on $\overline{Q(1/2)}$.

A corollary worth stating alone: **any finite-time singularity must have
$\limsup_{t\uparrow T}\lVert u(t)\rVert_{L^3}=\infty$.** Every surviving route below is a
way of arranging that.

### Type I / Type II

Koch–Nadirashvili–Seregin–Šverák, [arXiv:0709.3599](https://arxiv.org/abs/0709.3599),
primary: a singularity at $T$ is **Type I** iff $\sup_x|u(x,t)|\le C(T-t)^{-1/2}$;
**Type II** is by definition anything else ("slow blow-up"). Leray's lower bound
$\sup|u|\ge\varepsilon_1(T-t)^{-1/2}$ always holds.

> **Axisymmetric Type I is *not* unconditionally excluded.** KNŠŠ Thm 6.2 needs the
> Type-I bound **plus** exterior cylindrical decay $|u(x,t)|\le C/\sqrt{x_1^2+x_2^2}$ for
> radius $\ge R_0$ (e.g. a mild solution with fast-decaying data); Thm 6.1 assumes a
> global $|u|\le C/\sqrt r$. Treating axisymmetric Type I as closed outright would
> overclaim.

Seregin, [arXiv:2304.04045](https://arxiv.org/abs/2304.04045), studies a Type II scenario
via Euler scaling. An analytical program, **not** a construction.

## Consequence for this project

The profile-based program is **structurally inapplicable** to the Clay statement, not
merely unbridged: Lemma 1 forces $\beta=1/2$, Tsai Thm 2 closes $\beta=1/2$ under
hypotheses weak enough that no reasonable construction escapes, and Lemma 4 shows
admissible forcing cannot rescue it. Anything descended from Chen–Hou or the DeepMind
unstable-profile work is off this path by structure, whatever its merits as inviscid
model mathematics.

The surviving routes are tabulated once, in [Verdict](#verdict) at the end — not here,
to keep a single authority.

## The DSS reduction — derived, machine-checked, and then closed

Verified symbolically in [`dss.py`](dss.py) (all assertions pass). The singularity is
placed at the origin by translation — for Fefferman (C) the data start at $t=0$ and the
singularity sits at $(x_0,T)$ with $T>0$, so one first shifts $\xi=x-x_0$, $s=t-T$; then
$t\in[0,T)$ is $s\in[-T,0)$, and below $t$ denotes that shifted time. With
$$y=\frac{x}{\sqrt{-t}},\quad \tau=-\log(-t),\quad
u=(-t)^{-1/2}v(y,\tau),\quad p=(-t)^{-1}P(y,\tau)$$
turns unforced NS into, exactly,
$$\partial_\tau v + \tfrac12 v + \tfrac12 (y\cdot\nabla)v + (v\cdot\nabla)v + \nabla P
= \nu\Delta v,\qquad \nabla\cdot v=0. \tag{$\dagger$}$$
The $\lambda$-DSS condition leaves $y$ invariant, shifts $\tau\mapsto\tau-2\log\lambda$,
and the prefactor cancels, so
$$\textbf{DSS}\iff v(y,\tau)=v(y,\tau-S),\qquad S=2\log\lambda.$$

| ansatz | reduced problem | status |
|---|---|---|
| exact self-similar | **steady** solution of $(\dagger)$ — Leray at $a=1/2$ | closed (NRS, Tsai) |
| $\lambda$-DSS | **$S$-periodic** solution of $(\dagger)$ | closed — see below |

Lemma 1 has no DSS analogue, and that is structural: it plays a $t$-independent left
side against a $t$-dependent coefficient, whereas $(\dagger)$ is autonomous in $\tau$
and periodicity is a constraint, not a contradiction. Chae–Wolf's rigidity being
perturbative in $\lambda\to1^+$ is likewise explained — small $\lambda$ means small
period $S$, i.e. nearly steady.

### But exact DSS is closed for Clay-admissible data

$L^3$ is the scale-invariant norm: $\lVert\lambda u(\lambda\,\cdot)\rVert_{L^q}^q
= \lambda^{q-3}\lVert u\rVert_{L^q}^q$, invariant exactly at $q=3$. Hence exact DSS makes
$\lVert u(\cdot,t)\rVert_{L^3}$ **periodic** in $\tau$. Finite at a single time slice —
which Schwartz data guarantees — makes it bounded on the whole interval, so
$u\in L^\infty_tL^3_x$, and Escauriaza–Seregin–Šverák then removes the singularity.

> **For every $\lambda$.** The large-$\lambda$ exact-DSS hole exists only *outside* $L^3$,
> and such solutions cannot arise from Fefferman's condition (4). This is Chae–Wolf
> Remark 1.2. **Exact backward DSS is not a Clay route.**

### Two retractions

> **Retracted (energy).** A draft argued that Chae–Wolf Thm 1.1's bound
> $|u|\le C_*/(\sqrt{-t}+|x|)$ forces $u\sim1/|x|$ and hence infinite energy. Invalid:
> that is an **upper** bound and a *conclusion* of their theorem, so faster-decaying
> $L^2$ fields satisfy it too; it cannot bound decay from below. (The draft also called
> $\int r^2r^{-2}dr$ a logarithmic divergence; it is linear.) In fact the DSS scaling
> gives $E(\lambda^2t)=\lambda E(t)$ and
> $E(t)=(-t)^{1/2}\lVert v(\cdot,\tau)\rVert_{L^2}^2$, so if $v\in L^2$ the energy is
> finite and $\to0$ at the singular time. Finite energy and DSS are compatible.
>
> **Retracted (open route).** A draft then called Schwartz-decaying periodic DSS an open
> Clay route. It is not — the $L^3$ + ESS argument above closes it.

### The genuine hole — localized / asymptotic DSS

Chae–Wolf **Remark 1.2**, verbatim: *"If $u\in C((-\infty,0);L^3(\mathbb{R}^3))$, and
discretely self-similar, then $u\in L^\infty(-\infty,0;L^3(\mathbb{R}^3))$. Thus, in case
$p=3$, by using the result in [5], we get the full regularity $u$ in $\bar Q$."* That is
exactly the $L^3$+ESS argument reconstructed above — derived independently here, then
confirmed against the source.

**Thm 1.3**: for every $C_*>0$ there is $\lambda_*>1$ depending on $C_*$ *alone* ($\nu$
normalised to 1) such that $\lambda\in(1,\lambda_*)$ plus $|u|\le C_*/(\sqrt{-t}+|x|)$
forces $u\equiv0$. **Remark 1.4**: if $C_*$ is small enough, *all* $\lambda$ are closed.

**Is $\lambda\approx1$ essential or technical? Essential to the method.** The proof is a
compactness argument along $\lambda_j\to1$: rescale, extract local $L^3$ convergence and
uniform vorticity bounds via Arzelà–Ascoli, pick $k_j$ with $\lambda_j^{k_j}\to\mu$ so the
DSS relation passes to *continuous* self-similarity in the limit, then apply Tsai. It
works by degenerating DSS into exact self-similarity, which is only available near
$\lambda=1$. The authors do not claim it is removable.

**Why localization reopens the route.** Thm 1.5's convergence lives on **shrinking balls**
$B_{R\sqrt{t_*-t}}$ only. It controls nothing outside them, hence gives **no global $L^3$
bound**, so ESS cannot be invoked. Remark 1.6 adds that for $p>3$ the weight
$(t_*-t)^{(p-3)/2p}\to0$, so the displayed convergence is not guaranteed as ordinary local
$L^p$ convergence. Exact DSS supplies global scale invariance and therefore the bound;
localization destroys precisely that.

> **Sharp separator:** $\sup_{t<t_*}\lVert u(t)\rVert_{L^3(\mathbb{R}^3)}<\infty$. With
> it, ESS closes the route. Without it, localized/asymptotic DSS is **not** closed by this
> mechanism — the far field can carry growing $L^3$ mass while the local asymptotic
> relation still holds.
>
> **The hole is real, not cosmetic.** Chae–Wolf themselves state that nonexistence for
> nontrivial DSS remains open in the non-$L^3$ profile case.

### Forcing — Lemma 4

Fefferman (C) lets the claimant choose $f$ subject only to (5); only (A)/(B) say
*"we take $f(x,t)$ to be identically zero."* So an unforced map covers a subclass. (5)
with $\alpha=m=K=0$ gives $|f|\le C$. Verified in [`forcing.py`](forcing.py):

**Recentring — this matters.** Fefferman (C) starts the data at $t=0$ with a putative
singularity at $(x_0,T)$, $T>0$. The NS scaling symmetry is centred at a spacetime point,
so the correct relation is centred at $(x_0,T)$, **not** at $t=0$:
$$u(x,t)=\lambda\,u\big(x_0+\lambda(x-x_0),\;T+\lambda^2(t-T)\big).$$
An earlier draft used $u(x,t)=\lambda u(\lambda x,\lambda^2t)$, which centres the symmetry
on the *initial slice* and formally describes a singularity at $t=0$, before the data
exists. Put
$$\xi=x-x_0,\quad s=t-T,\quad w(\xi,s)=u(x_0+\xi,T+s),\quad \tilde f(\xi,s)=f(x_0+\xi,T+s),$$
so $t\in[0,T)$ becomes $s\in[-T,0)$. Translation is an NS symmetry, so $w$ solves the same
system with force $\tilde f$, and every formula below is stated for $w,\tilde f$.

**Continuous self-similarity.** Every retained term carries $\tau^{-(\alpha+1)}\to\infty$,
so a bounded $f$ is subcritical. That alone does **not** give $f\equiv0$: with arbitrary
admissible forcing the *pressure* need not share the self-similar scaling.

> **Retracted.** A draft argued that exact self-similarity of $u$ makes the whole left
> side $\tau^{-3/2}L(y)$, hence $f=\tau^{-3/2}L(x/\sqrt\tau)$, hence $f\equiv0$. False.
> **Counterexample:** $u\equiv0$ is exactly self-similar; take any admissible $\phi$ and
> set $f=\nabla\phi$, $p=\phi$. The equation holds with $u\equiv0$ and $f\not\equiv0$.
> Only the divergence-free part of $f$ is pinned down.

**The covariance is *derived*, not assumed.** With
$N(u,p)=\partial_tu+(u\cdot\nabla)u+\nabla p-\nu\Delta u$,
$\tilde w=\lambda w(\lambda\xi,\lambda^2s)$ and $\tilde q=\lambda^2q(\lambda\xi,\lambda^2s)$,
every term carries exactly $\lambda^3$ (checked termwise), so
$N(\tilde w,\tilde q)(\xi,s)=\lambda^3N(w,q)(\lambda\xi,\lambda^2s)$. Exact covariance of
$w$ gives $\tilde w=w$, and subtracting the two forms of the equation leaves
$$\lambda^3\tilde f(\lambda\xi,\lambda^2s)-\tilde f(\xi,s)=\nabla(\tilde q-q),$$
a pure gradient. The Leray projector $\mathbb{P}=I-\nabla\Delta^{-1}\mathrm{div}$ has a
degree-0 homogeneous symbol, so it kills gradients **and commutes with dilations**:
$$(\mathbb{P}\tilde f)(\xi,s)=\lambda^3(\mathbb{P}\tilde f)(\lambda\xi,\lambda^2 s).$$

- **Continuous case:** this holds for every $\lambda>0$; letting $\lambda\to0^+$ gives
  $\lvert\mathbb{P}\tilde f\rvert\le\lambda^3\lVert\mathbb{P}\tilde f\rVert_\infty\to0$.
- **DSS case** (single $\lambda>1$): rewrite as
  $(\mathbb{P}\tilde f)(\xi/\lambda^n,s/\lambda^{2n})=\lambda^{3n}(\mathbb{P}\tilde f)(\xi,s)$.
  The sampled times $s/\lambda^{2n}\to0^-$ **approach the singularity and stay inside
  $[-T,0)$**; iterating the other way would send $s\to-\infty$, off the interval of
  existence. Unboundedness then contradicts Lemma 4.1.

Either way $\mathbb{P}\tilde f\equiv0$: $f$ is a pure gradient, absorbable into the
pressure, and the equation is **effectively unforced**. Tsai's Thm 2 carries no pressure
hypothesis, so it applies with the shifted pressure $q-\phi$.

> **Scope — precise.**
> **Closed:** *exactly* self-similar or *exactly* DSS $u$ with any admissible $f$. Such a
> $u$ forces $\mathbb{P}f\equiv0$, so the rigidity theorems apply unchanged.
> **Not closed:** *asymptotically* self-similar / localized DSS with forcing — there $u$
> is only approximately covariant, the derivation above yields nothing, and $f$ need not
> be DSS at all.
> **Not closed:** wholly non-self-similar forced blowup. No forced 3D NS blowup
> construction exists; known forced results are Euler or Boussinesq with forces far
> rougher than (5), and Albritton–Brué–Colombo is weak non-uniqueness, not breakdown.

### Lemma 4.1 — $\mathbb{P}f$ is uniformly bounded. Proved.

The closures above need $\sup_t\lVert\mathbb{P}f(\cdot,t)\rVert_\infty<\infty$. This was
previously *asserted*; here is the proof (also written out in
[`forcing.py`](forcing.py), whose sympy checks cover only the constants).

$\mathbb{P}$ is the multiplier $P(\xi)=I-\xi\otimes\xi/|\xi|^2$, the orthogonal projection
onto $\xi^\perp$, so $\lVert P(\xi)\rVert_{op}=1$ for $\xi\ne0$. By Fourier inversion,
$$\lVert \mathbb{P}f(\cdot,t)\rVert_\infty \le (2\pi)^{-3}\lVert P(\xi)\hat f(\xi,t)\rVert_{L^1_\xi} \le (2\pi)^{-3}\lVert \hat f(\cdot,t)\rVert_{L^1}.$$
Split at $|\xi|=1$.

*Low, $|\xi|\le1$.* $|\hat f(\xi,t)|\le\lVert f(\cdot,t)\rVert_{L^1}$. By (5) with
$\alpha=0,m=0,K=4$, $|f(x,t)|\le C_{004}(1+|x|+t)^{-4}\le C_{004}(1+|x|)^{-4}$, and
$\int_{\mathbb{R}^3}(1+|x|)^{-4}dx=\tfrac{4\pi}{3}$ since $4>3$ — uniformly in $t$. So
$\int_{|\xi|\le1}|\hat f|\,d\xi\le\tfrac{4\pi}{3}A_0$, $A_0:=\sup_t\lVert f(\cdot,t)\rVert_{L^1}$.

*High, $|\xi|\ge1$.* Some coordinate has $|\xi_j|\ge|\xi|/\sqrt3$. Taking $\alpha=4e_j$ and
using $(i\xi)^\alpha\hat f=\widehat{\partial^\alpha f}$,
$|\xi_j|^4|\hat f|\le\lVert\partial_j^4f(\cdot,t)\rVert_{L^1}$, so
$|\hat f(\xi,t)|\le 9A_4|\xi|^{-4}$ with
$A_4:=\sup_t\max_{|\alpha|=4}\lVert\partial^\alpha f(\cdot,t)\rVert_{L^1}$, finite by (5)
with $|\alpha|=4$, $K=4$. Since $\int_{|\xi|\ge1}|\xi|^{-4}d\xi=4\pi$,
$\int_{|\xi|\ge1}|\hat f|\,d\xi\le36\pi A_4$.

Hence $\sup_t\lVert\mathbb{P}f(\cdot,t)\rVert_\infty\le(2\pi)^{-3}\big[\tfrac{4\pi}{3}A_0+36\pi A_4\big]<\infty$. $\blacksquare$

> The lemma uses (5)'s **full** strength. $\mathbb{P}$ is *not* bounded on $L^\infty$, so
> $|f|\le C$ alone would not do — the decay *and* four derivatives are both load-bearing.

### Verdict

| Route | Status |
|---|---|
| Exact self-similar, unforced | closed — Lemma 1 + Tsai Thm 2 |
| *Exactly* self-similar / DSS **with** admissible $f$ | closed — Lemma 4: $\mathbb{P}f\equiv0$, reduces to unforced |
| Exact backward DSS in $C_tL^3$ | closed — Chae–Wolf Rem. 1.2, via ESS (verified 2026-08-15) |
| Exact backward DSS outside $L^3$ | open, but unreachable from Schwartz data |
| **Localized / asymptotic DSS with unbounded $L^3$, i.e. $\limsup_{t\uparrow T}\lVert u(t)\rVert_{L^3}=\infty$** | **open — the live route** |
| Type II | open; one scenario studied, no construction |
| Non-self-similar forced blowup | open; neither constructed nor excluded |
| **Asymptotic / localized profile routes with forcing** | **open — Lemma 4 does not reach them** |

**2026-08-15:** both primaries retrieved and read — the two `[UNVERIFIED]` tags are
discharged (see Verification status above); the DSS closure now rests entirely on
verified sources. **Next deliverable:** decide whether unboundedness of
$\lVert u(t)\rVert_{L^3}$ — $\limsup=\infty$, *not* necessarily a full limit — is
compatible with a localized asymptotically-DSS ansatz. That is the one question standing
between the live route and a construction attempt.

## Scope discipline

- The canonical target is **(C)** above, on $\mathbb{R}^3$. (D) on the torus is out of
  scope — see the Domain scope note at the top. Progress gates refer to (C) and nothing else.
- `../ccf/` is **technique validation only, and currently paused**: a 1D inviscid model
  with no viscosity bridge and no domain bridge. It produced reusable machinery (an
  exactly-invertible Hilbert-transform basis) and one hard-won lesson about validating
  operators on the decay class you actually need. It is not progress on this problem and
  is not reported as such.
- No experiment enters this directory unless it is tied to a named lemma or hole in the
  obstruction map above.
