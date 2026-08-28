# Union-closed sets conjecture — active target

> **MACHINE-VERIFIED FINITE CERTIFICATE (2026-08-21).** Campaign
> [`cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8`](campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/)
> replayed all eight traces and certified the explicit relaxed functional
> \(\Phi_{\rm rel}\ge0\) on the feasible five-parameter family at exact rational
> \(t=0.3820660112501052\ge\psi+10^{-4}\): 488,465,854 nodes, no residuals.
> Fresh direct isolated and hardened clean-stage replays on 2026-08-26 both
> passed all eight slices; their canonical report hashes are
> `6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`
> and
> `57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.
> A third secure source-independent replay then passed all eight traces from
> byte zero, with a separately implemented interval-AD arithmetic layer,
> 343-file/runtime-image pins, exact tallies, and zero residuals.  Its trusted
> report-lock raw SHA-256 is
> `af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
> composite canonical SHA-256 is
> `4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.
> This is not an end-to-end machine proof of the union-closed consequence.
> Support reduction, the Margin Lemma, and the set-family bridge are ordinary
> human-audited mathematics; see [`AUDIT.md`](AUDIT.md) and
> [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The overall result is a
> certificate-backed candidate theorem pending external review.

**Conjecture (Frankl, 1979).** For any union-closed family
$\mathcal{F}\subseteq2^{[n]}$ with $\mathcal{F}\neq\{\emptyset\}$ there exists
$i\in[n]$ contained in at least $|\mathcal{F}|/2$ of the sets.

Claims below are labeled by evidence type. Machine verification covers only the
named executable facts; human and cited steps remain explicit.

## Why this target

Selected over Navier–Stokes on one criterion: **every primary is on arXiv**, so
no step depends on a paywalled source. The NS attempt stalled with NRS and ESS
publisher-blocked. Here the entire post-Gilmer literature is retrievable, and
the state of the art is *a constant defined by an optimisation problem*, which
is machine-checkable.

## Current frontier — closure-defect scalar route

The authoritative current verdict is in [`../RESULTS.md`](../RESULTS.md);
sections below preserve the derivation and its recorded corrections.

**PROVED/CERTIFIED (2026-08-25): Gate B is resolved in its stated scope, and
the scalar closure-defect route fails.** Campaign IV first refuted the
closure-blind arbitrary-family target $A_+>0$ with a non-union-closed $n=6$
family. The complete declared $n=7$ block census then refuted
$$
A_+(\mathcal F)+\frac1{50}\varepsilon_\vee(\mathcal F)>0,\qquad
\varepsilon_\vee=\Pr[X\vee Y\notin\mathcal F].
$$
Its certified normalized 45-row extremizer $\mathcal B$ is the union of the
seven cells
$$
(0,1),(0,3),(1,0),(1,2),(1,5),(2,0),(2,1)
$$
for a $2+5$ coordinate partition. Every coordinate count is 18, the exact
defect is $64/81$, and a 256-bit Arb recurrence proves
$$
A_+(\mathcal B)
\le-0.0286491867944683\ldots<-\frac7{250}.
$$

Gate B asked whether
$$
c_{\rm cl}^\star=
\sup_{\substack{\mathcal F\ {\rm cap/Reimer}\\A_+(\mathcal F)<0}}
\frac{-A_+(\mathcal F)}{\varepsilon_\vee(\mathcal F)}
$$
is finite. Put $\mathcal F_k=\mathcal B^{\boxtimes k}$ on disjoint
seven-coordinate blocks. A fixed-order Bellman induction proves exact
tensorization
$$
A_+(\mathcal F_k)=kA_+(\mathcal B),
\qquad
\varepsilon_\vee(\mathcal F_k)=1-\left(\frac{17}{81}\right)^k.
$$
Every $\mathcal F_k$ is normalized, has size $45^k$, meets the cap with each
coordinate count $18\cdot45^{k-1}=(2/5)45^k$, and satisfies Reimer strictly
because $45^5<2^{28}$. Therefore
$$
\boxed{\frac{-A_+(\mathcal F_k)}
{\varepsilon_\vee(\mathcal F_k)}
>\frac{7k}{250}\longrightarrow+\infty},
\qquad
\boxed{c_{\rm cl}^\star=+\infty}.
$$
No fixed scalar coefficient can repair this $A_+$ relaxation over all
cap/Reimer families.

**Independence audit (2026-08-26).** The finite input to that conclusion no
longer rests on one arithmetic stack. A second standalone checker,
`uc/gate_b/verify_gate_b_dyadic.py`, imports nothing outside the Python
standard library, works in $2^{-160}$ outward-rounded dyadic intervals with
rational transcendental tail bounds, and makes no Bellman clamp decision at all
because it relaxes every nondegenerate action interval to $[0,1]$. It proves
$A_+(\mathcal B)<-1/40$, which already gives ratio $>k/40\to\infty$. The
pristine system interpreter reproduces its certificate byte-identically. A
non-factored audit of the 2,025-row square and a 292-case product
falsification search — 205 of them with unequal block dimensions — found no
counterexample to the tensorization lemmas, and the shared algebra
($\sum P_{ac}=1$, shift invariance of $D$, and
$\max_s[sD+h(s)]=\log_2(1+2^D)$) is now verified symbolically.

**Second-base audit (2026-08-26).**  The conclusion no longer depends on
\(\mathcal B\).  The separate 25-row \(3+3\) family \(\mathcal D\), with cells
\((0,0),(0,1),(1,0),(1,2),(2,1)\), has exact closure success \(181/625\).
A 160-bit dyadic checker evaluates all 720 orders without symmetry reduction
and proves \(A_+(\mathcal D)<-1/80\); a 256-bit Arb checker independently proves
\(A_+(\mathcal D)<-17/1250\).  Its powers are normalized, meet the cap exactly,
satisfy Reimer by \(25^5<2^{24}\), and have ratio \(>k/80\to\infty\).  A
non-factored 625-row square audit found zero Bellman product gap on five hostile
global orders.

**Exact-arithmetic audit (2026-08-27).**  The finite input is now an exact
rational fact about the *true* objective.  `uc/gate_b/verify_gate_b_rational.py`
uses only the Python standard library, keeps the exact feasible action interval
\([s^*,U]\) — no clamp classification and no relaxation — and returns two-sided
enclosures over **all** 720 and **all** 5040 coordinate orders:
\(A_+(\mathcal D)\in[-0.013672107732177773561859915661,
-0.013672107732177773561859912998]\) and
\(A_+(\mathcal B)\in[-0.028649186794468317816026981933,
-0.028649186794468317816026976455]\), widths below \(6\times10^{-27}\),
re-deriving the sharp \(-17/1250\) and \(-7/250\) bounds and containing both Arb
values.  Its certified primitives are a bit-by-bit \(\log_2\), ceiling
square-root towers for \(2^x\), and a concavity case split.  The all-order value
lists partition into 10 classes of size 72 and 21 of size 240, recovering the
automorphism data rather than assuming it.  The same arithmetic re-attacks the
product lemmas on 35 products over all 15,010 global orders with worst
fixed-order Bellman discrepancy exactly zero.  Finally the growth rate is
pinned: \(c_{\rm cl}^\star(n)>n/250-7/250\), while \(Q\ge0\) and \(C_+\ge0\)
give \(-A_+\le\log_2m\le n\), so the ratio is \(\Theta(n)\) whenever
\(\varepsilon_\vee\ge\varepsilon_0>0\).

This does not refute a union-closed statement: all three bases are explicitly
non-UC, and every power defect exceeds its base's. The separate local-stability
problem restricted to \(\varepsilon_\vee\to0\) remains open, but that
restriction was not present in the authoritative Gate B supremum.

**Local frontier and the size ceiling (2026-08-27).** The defect multiplies
through success, so no Cartesian power reaches small defect and the theorem's
mechanism says nothing about the local regime. Three new facts:

*Growth closed, with an improved constant.* Incidence is at most
\(n\lfloor2m/5\rfloor\) while Reimer demands at least \(m\log_2m/2\), so every
admissible family obeys \(\log_2m\le4n/5\) (exactly: \(m^5\le2^{4n}\)). With
\(Q,C_+\ge0\) this gives \(-A_+\le4n/5\) unconditionally, and a
cloned-coordinate base at \(n=8\) raises the certified lower slope from
\(1/250\) to \(3/640\), an improvement of exactly \(75/64\). So
\((3/640)n-3/80\le\Lambda(n)\le4n/5\): the numerator's growth is \(\Theta(n)\)
and only the denominator is open. The clone is legitimate because **no lemma in
the product argument uses separation**, and Reimer for every power reduces to
the single base-level inequality \(m^m\le2^{2I}\), i.e. to \(I\ge R_m\).

*Frontier moved to eight coordinates.* The complete \(S_2\times S_6\) class at
\(n=8\) (2,272 canonical families) yields, after exact rational re-evaluation,
\(\min\{\varepsilon_\vee:A_+<0\}\le1144/1875=0.6101\) among separating families
and \(\le139/245=0.5673\) among all admissible ones, down from
\(1336/2025=0.6598\); and the first repair ratios here to beat the published
\(0.0362591\) — \(0.0370512\) separating, \(0.0580391\) otherwise. Sizes 70 and
75 are *infeasible* at \(n=7\) (\(196<R_{70}=215\)): the eighth coordinate
supplies the missing incidence, and in the 70-row families it is an exact
duplicate of the first, so coordinate cloning is what raises the ceiling.

*Why low defect needs a dominant set.* The cap forces
\(\mathbb E|X\vee Y|\ge\frac85\bar s\) while a successful union has size at most
\(M\), so \(\varepsilon_\vee\ge(\frac85\bar s-M)/(n-M)\) and the ratio is
\(O(n)\) whenever \(M<\frac85\bar s\). Reaching the local regime therefore
requires a set of size \(\ge\frac45\log_2m\). Separately, 25 and 45 are the
largest admissible sizes at \(n=6,7\) and **no set can be added to any of the
three \(n\le7\) bases**, so adding the missing unions is empty there. Emptying
the regime by proof would still give a \(2/5\) Frankl bound, about \(0.018\)
beyond the published frontier.

The complete statement, proof, ten exactly certified bases, append-only
searches, and paper draft are in [`gate_b/`](gate_b/); the current verified
state alone is in [`gate_b/CURRENT_STATUS.md`](gate_b/CURRENT_STATUS.md).

Reproduce from `math/`:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/test_gate_b.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/verify_gate_b.py
OMP_NUM_THREADS=1 /Library/Developer/CommandLineTools/usr/bin/python3 -B uc/gate_b/verify_gate_b_dyadic.py
OMP_NUM_THREADS=1 /Library/Developer/CommandLineTools/usr/bin/python3 -B uc/gate_b/verify_gate_b_n6_dyadic.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/verify_gate_b_n6_arb.py
OMP_NUM_THREADS=1 /Library/Developer/CommandLineTools/usr/bin/python3 -B uc/gate_b/verify_gate_b_rational.py --bases n6,n7 --orders all
OMP_NUM_THREADS=1 /Library/Developer/CommandLineTools/usr/bin/python3 -B uc/gate_b/audit_tensorization_exact.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/audit_n6_square_direct.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/audit_tensorization.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/audit_square_direct.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/gate_b/search_n8_block_symmetric.py --k 1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/shapley_n7_block_symmetric.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/shapley_n7_falsifier_cert.py
```


## The ladder and current evidence status

| constant | source | status |
|---|---|---|
| $0.01$ | Gilmer, [2211.09055](https://arxiv.org/abs/2211.09055) | proved |
| $\psi=\frac{3-\sqrt5}{2}=0.3819660112501051$ | Chase–Lovett, Alweiss–Huang–Sellke, Pebody, Sawin | proved analytically |
| $0.3820660112501052$ | this project | **candidate theorem**: human-audited reduction and entropy bridge plus replayed Arb certificate |
| $c^*=0.3823455333667027$ | Yu [2212.00658](https://arxiv.org/abs/2212.00658), Cambie [2212.12500](https://arxiv.org/abs/2212.12500) | claimed; finite-dimensional verification gap remains |
| $>c^*$, non-explicit | Liu [2306.08824](https://arxiv.org/abs/2306.08824) Thm 6 | proved, no explicit value |
| $0.382709087918741$ | Liu, Thm 13 | conditional: Hypothesis 1 now has a human-audited candidate proof; Hypothesis 2 remains open |

Cambie states the gap in his own paper, verbatim:

> *"As such, the remaining problem is a minimisation problem for which the
> statement can be checked by a computer. This implies that an exact rigorous
> calculus proof is missing."*

Liu's Theorem 13, verbatim:

> *"Under the positive-semidefiniteness hypothesis in Section V-A and the
> hypothesis of the global minimizer structure in Section V-B, the constant in
> the union-closed sets conjecture can be improved to $c'$ in (93)."*

Liu checked Hypothesis 1 on a float64 grid (roundoff-scale extremal
eigenvalue) and Hypothesis 2 with roughly \(10^5\) MATLAB restarts.  The
human-audited residual-kernel proof in [`../LIU_H1/AUDIT.md`](../LIU_H1/AUDIT.md)
now supplies a candidate resolution of Hypothesis 1, with independent exact
algebra checks. Hypothesis 2 remains open, so secondary claims that
\(0.38271\) was already proved remain overstatements.

## Result 0 — the OR-entropy kernel has exactly one positive square

[`decomposition.py`](decomposition.py), [`reduction.py`](reduction.py). This is
the structural fact the rest rests on. Work in $x=1-p$; since
$p\oplus q=1-(1-p)(1-q)$ and $h(1-z)=h(z)$, the iid kernel is $H(xy)$ with
$H=\ln2\cdot h$. Put

$$T(u)=u+(1-u)\ln(1-u)=\sum_{n\ge2}\frac{u^n}{n(n-1)},$$

so $T(xy)$ is a **positive-definite kernel** (positive power-series coefficients).
Then, from $H(u)=-u\ln u+u-T(u)$ and $-xy\ln(xy)=-xy\ln x-xy\ln y$:

$$\boxed{\;H(xy)=g(x)g(y)-a(x)a(y)-T(xy),\qquad a(x)=-x\ln x,\;\; g(x)=x(1-\ln x)\;}$$

**Verified as an exact algebraic tautology in sympy**, and numerically on 200
random *signed* measures to $10^{-40}$ — so it is a form identity, not a
coincidence on probability measures. Gram spectra on grids of 40/120/300 points
show exactly **one** positive eigenvalue at every resolution, explaining the
single positive mean-zero eigenvalue [`kernel.py`](kernel.py) found numerically:
it is the functional $G(\mu)=\int g$.

### Result 0′ — the Gilmer term in closed form

Rearranged, with $B=\int(1-p)\,d\mu=1-\text{mean}$ and
$W(x,y)=xy-xT(y)-yT(x)+T(xy)=\sum_{n\ge2}\frac{(x-x^n)(y-y^n)}{n(n-1)}$ (PSD):

$$\boxed{\;Q(\mu)=2B\,L(\mu)-E(\mu),\qquad E=\tfrac1{\ln2}\iint W\,d\mu\,d\mu\ \ge0\;}$$

equivalently $\ln2\cdot E=\sum_{n\ge2}\frac{(M_1-M_n)^2}{n(n-1)}$ with
$M_n=\int(1-p)^n d\mu$ — every term a squared **moment defect** $M_1-M_n\ge0$.
So the full two-strategy functional is

$$F(\mu)=\big[2(1-\alpha)B-1\big]L+\alpha C-(1-\alpha)E.$$

**This is sharp for the pure iid problem — $\alpha=0$ only.** At $\alpha=0$ the
diagonal $W(x,x)/\ln2=2xh(x)-h(x^2)$ collapses $F(\delta_p)$ to $h(x^2)-h(x)$, so
$F\ge0\iff(1-p)^2\ge p\iff p\le\psi$, with equality exactly at the golden point
$x=1/\varphi$ where $x^2=1-x$. The barrier $\psi$ *is* the statement that the PSD
defect $E$ first overtakes $(2B-1)L$.

> **The $\alpha=0$ hypothesis is essential and an earlier version of this section
> dropped it.** For $\alpha>0$ the coupled term does **not** vanish at a point
> mass: $s^*(p,p)=1/2$ on $(1/4,1/2)$, so $h(s^*)=1$ near $\psi$ and
> $$F(\delta_p)=(1-\alpha)h\big((1-p)^2\big)+\alpha h\big(s^*(p,p)\big)-h(p),
> \qquad F(\delta_\psi)=\alpha\big[1-h(\psi)\big]=0.0405813\,\alpha>0.$$
> So point masses do not cap at $\psi$ once $\alpha>0$.
>
> **A second correction, to the fix itself.** I then claimed the point-mass
> threshold exceeds $c^*$ for *every* $\alpha>0$. False: at $p=c^*$ the pure-iid
> part is already negative, $h((1-c^*)^2)-h(c^*)=-5.8929\times10^{-4}<0$ because
> $c^*>\psi$, and $F(\delta_{c^*},\cdot)$ is **affine** in $\alpha$, so it stays
> negative for all small $\alpha$. The exact crossover is
> $$\alpha_{\rm crit}=\frac{h(c^*)-h\big((1-c^*)^2\big)}{1-h\big((1-c^*)^2\big)}
> =0.0144054585140081,$$
> and *point-mass threshold $>c^*$ $\iff$ $\alpha>\alpha_{\rm crit}$*. Cambie's
> $\alpha=0.0356069$ exceeds it, so at **that** $\alpha$ point masses are not
> binding and the 2-atom obstruction $\{1,b\}$ is; for $\alpha<\alpha_{\rm crit}$
> point masses bind strictly below $c^*$. My test had sampled
> $\alpha\in\{0,0.0356069,0.1,0.3\}$ and skipped the entire interval where the
> claim fails. Now bracketed on both sides in
> [`reduction.py`](reduction.py) §6 ($0.0144$ below, $0.0145$ above), with the
> crossover asserted rather than sampled.

Consequently the $\psi$ theorem is equivalent to the pure moment inequality

$$\sum_{n\ge2}\frac{(M_1-M_n)^2}{n(n-1)}\;\le\;(2M_1-1)\int H(x)\,d\mu
\qquad(M_1\ge1/\varphi),$$

with no entropy on the left. Checked against 4000 random laws: no violations.

**What it does not give.** Dropping the variance ($E\le\frac1{\ln2}\int W(x,x)d\mu$,
Jensen, tight only at point masses) yields $\kappa^*=\max_x[2x-h(x^2)/h(x)]=
0.3049467$ attained at $x=0.8835$, i.e. the constant $(1-\kappa^*)/2=0.3475266<\psi$.
The maximiser sits at $p=0.117$, far below the mean bound, so the pointwise route
wastes the constraint. **The variance term is load-bearing** — which is why a
support reduction, not an inequality, is required.

## Result 0″ — support reduction, corrected: at most 3 *pair-orbits*

[`decomposition.py`](decomposition.py), [`support3.py`](support3.py).

> **Retraction.** An earlier version claimed **Theorem B**: the infimum is
> attained on laws with $\le3$ *marginal* atoms. That required
> $\Psi=(1-\alpha)\Psi_Q+\alpha C-L$ to be concave. It is not. My own argument —
> if $M_i$ couples $\mu_i$ then a convex combination couples the combination, so
> $C(\theta\mu_1+(1-\theta)\mu_2)\le\theta C(\mu_1)+(1-\theta)C(\mu_2)$ — is the
> definition of **convexity**, and I read it as concavity. Monge–Kantorovich
> confirms it independently: $C=\sup_\phi2\int\phi\,d\mu$ is a supremum of linear
> functionals. Witness: $C(\delta_{0.3})=1$, $C(\delta_{0.6})=0.970951$, but the
> half-half mixture gives $0.970951$, **below** the chord $0.985475$ by
> $1.45\times10^{-2}$; worst over 1500 random segments, $-0.452$. Bauer does not
> apply, and the "$\le3$ marginal atoms", "same dimension as Yu's 5 parameters",
> and "the $k=2..6$ searches confirm a theorem" claims are all withdrawn. The
> searches remain valid as evidence; they were never a proof.

The repair keeps the symmetric coupling itself as the variable, which turns the
convex term linear. Let $\Delta=\{(p,q):0\le p\le q\le1\}$, $c(p,q)=h(s^*(p,q))$,
and for $\nu\in\mathcal P(\Delta)$ let $\mu_\nu=\int(\delta_p+\delta_q)/2\,d\nu$
be the induced marginal. Put
$\Phi(\nu)=(1-\alpha)Q(\mu_\nu)+\alpha\int c\,d\nu-L(\mu_\nu)$.

**Theorem B′.**
$\inf\{F(\mu):\int p\,d\mu\le t\}=\inf\{\Phi(\nu):\nu\in\mathcal P(\Delta),\ \mathbb E_{\mu_\nu}[p]\le t\}$,
*and the right side is attained at some $\nu$ with at most **3 atoms** — at most 3
unordered pairs, hence a marginal with at most **6** atoms.*

*Proof.* Symmetric probability measures on $[0,1]^2$ correspond to
$\nu\in\mathcal P(\Delta)$ with marginal $\mu_\nu$. For any $\nu$,
$\int c\,d\nu\ge C(\mu_\nu)$ (it is *some* symmetric coupling of its own
marginal), so $\Phi(\nu)\ge F(\mu_\nu)$; conversely an optimal coupling attains
equality. Hence the infima agree. Now $\int c\,d\nu$ is **linear** in $\nu$, and
so are $\mu_\nu$ and therefore $L$, the mean, and $G$. By Result 0,
$Q(\mu_\nu)=G^2/\ln2+\Psi_Q(\mu_\nu)$ with $\Psi_Q=-[A^2+\iint T]/\ln2$ concave —
minus a sum of squares of *linear* functionals — hence concave in $\nu$. So
$\Phi=(1-\alpha)G^2/\ln2+[\text{concave}]$, with $G$ the only convex direction.
Freeze $G=\gamma$: on $\mathcal K_\gamma$ (convex, weak-\* compact, $\Phi$
continuous) Bauer's minimum principle puts the minimum at an extreme point, and
$\mathcal K_\gamma$ is cut from the positive cone by three moment conditions —
mass, the mean inequality, $G$ — so its extreme points carry $\le3$ atoms.
$\qquad\blacksquare$

**Scope, plainly.** This is *weaker* than the retracted claim: 3 pairs are 6
coordinates $+$ 2 free weights $=$ **8 parameters**, against Yu's 5, so B′ does
**not** match the dimension of the literature's family. It is a valid finite
reduction where the literature has an invalid one, and nothing more. The
obstruction law happens to need only **2** pair-orbits —
$\nu^*=\{b,b\}@0.842245+\{b,1\}@0.157755$, mean exactly $c^*$ — so B′ covers the
binding instance comfortably.

Checked operationally: the reformulation is lossless ($\Phi$ at the LP-optimal
coupling equals $F$ to $10^{-9}$), and searches over $k=2,\dots,5$ pair-orbits
never beat $k\le3$, attaining exactly $1$ at $t=c^*$ for every
$\alpha\in\{0,0.0356,0.2,0.5,1\}$.

## Result 0‴ — the margin lemma: scale-free, and tight at the obstruction

[`margin_lemma.py`](margin_lemma.py). This is what makes the certification
well-posed. Recall $E=\frac1{\ln2}\iint W$ with $W$ PSD; write
$W=\langle\psi(x),\psi(y)\rangle\ln2$ for its feature map, so

$$\|\psi(x)\|^2=W(x,x)/\ln2=\rho(p)h(p),\qquad
\rho(p):=2(1-p)-\frac{h(2p-p^2)}{h(p)}.$$

Crucially $W(0,\cdot)=W(1,\cdot)=0$, so $\psi$ vanishes **exactly** at the two
sink points and nowhere else. Put $\sigma(\mu)=\int\sqrt{\rho(p)h(p)}\,d\mu$.

**Margin Lemma.** *For $L(\mu)>0$,*
$$F(\mu)\ \ge\ L(\mu)\,\Lambda(\mu),\qquad
\Lambda:=2(1-\alpha)B-1+\alpha\frac{C}{L}-(1-\alpha)\frac{\sigma^2}{L},$$
*with equality whenever $\mu$ has at most one atom outside $\{0,1\}$.*

*Proof.* $E=\|\int\psi\,d\mu\|^2\le(\int\|\psi\|\,d\mu)^2=\sigma^2$ by the
triangle inequality; substitute and divide by $L$. Equality needs the integrand
collinear, which is automatic if only one atom has $\psi\ne0$. $\;\blacksquare$

Two properties do all the work.

**Scale-free, so the degeneracy is gone.** $\sigma^2$ is *quadratic* in the
non-sink mass while $L$ is *linear*, so $\sigma^2/L\to0$ at the sink and $\Lambda$
tends to $2(1-\alpha)(1-t)-1>0$. On the family that defeats every bound on $F$ —
$\mu=(1-u-w)\delta_0+u\delta_1+w\delta_{0.1165}$ with the mean pinned at $t$ —
$F$ and $L$ collapse through ten orders of magnitude while $\Lambda$ sits at
$0.1918625$. And $\sigma^2\le\int\rho h\,d\mu\le(\max\rho)L$ gives
$\sigma^2/L\le\max_p\rho=0.3049467$ always: $\Lambda$ is bounded, never singular.

**Tight where it counts.** The obstruction $\mu^*=a\delta_1+(1-a)\delta_b$ has one
atom outside $\{0,1\}$ (namely $b$; $p=1$ *is* a sink point), so
$\Lambda(\mu^*)=F(\mu^*)/L(\mu^*)$ exactly — verified to $10^{-12}$.

> Using **Jensen** ($E\le\int\rho h\,d\mu$) instead loses $1.70\times10^{-2}$
> there, reporting $\Lambda=-1.94\times10^{-2}$ against the true $0$; and it is
> *negative on the whole* sink-plus-one-atom family, since $\int\rho h\,d\mu$ is
> only linear in $w$. Recorded because I tried it first.

**Theorem B″.** Put $\Phi:=L\Lambda=2(1-\alpha)BL-L+\alpha\int c\,d\nu
-(1-\alpha)\sigma^2$ on pair-orbit measures. Since $2BL=(B+L)^2-B^2-L^2$,
$$\Phi=(1-\alpha)(B+L)^2+[\text{concave in }\nu]+[\text{linear in }\nu],$$
the concave part being $-(1-\alpha)(B^2+L^2+\sigma^2)$. So $\Phi$ has exactly one
convex direction, $B+L$; freezing it, Bauer plus three moment conditions gives
$\inf\Phi$ on $\le3$ orbits. Hence

$$F\ge0\text{ on all feasible }\mu\ \Longleftarrow\ \Lambda\ge0
\text{ on feasible }\nu\text{ with}\le3\text{ orbits},$$

an 8-parameter problem with no inner LP, no degeneracy, and no slack at the
binding case — **superseded by Theorem B‴ below, which reaches 5 parameters**.
Search evidence (not a certificate):

| $t$ | $t-\psi$ | $\min\Lambda$ | $F/L$ at the minimiser |
|---|---|---|---|
| $\psi$ | $0$ | $+6.185390\times10^{-4}$ | same |
| $\psi+10^{-4}$ | $10^{-4}$ | $+4.555637\times10^{-4}$ | same |
| $\psi+2\times10^{-4}$ | $2\times10^{-4}$ | $+2.925862\times10^{-4}$ | same |
| $c^*$ | $3.80\times10^{-4}$ | $-0.0000000000$ | $0$ |

The search table is historical pre-certificate evidence. Campaign I has since
completed the branch-and-bound; current acceptance and limitations are in
`AUDIT.md` and `REPRODUCIBILITY.md`.

### Result 0⁗ — Theorem B‴: the reduction lands on **two** pair-orbits

B′ and B″ each had to freeze a functional chosen to absorb a *convex* direction
($G$, or $B+L$), leaving three moment conditions and three orbits. Freezing
$B=1-\text{mean}$ instead removes the bilinear term outright — and $B$ is the
*only* place the mean enters $Q$, which is exactly what Theorem A′ makes visible.

**Theorem B‴ (HUMAN-AUDITED).** *$\inf\{F(\mu):\text{mean}\le t\}$ is attained
at a pair-orbit measure with at most **2** atoms — a marginal with at most 4 atoms.*

*Proof.* On the slice $\{\nu:B(\nu)=\beta\}$, Theorem A′ gives $Q=2\beta L-E$, so
$$\Phi\big|_{B=\beta}=\big(2(1-\alpha)\beta-1\big)L(\nu)+\alpha\!\int\! c\,d\nu-(1-\alpha)E(\nu).$$
$L$ and $\int c\,d\nu$ are linear in $\nu$ and $E$ is a **PSD** quadratic form in
$\nu$, so $-E$ is concave and $\Phi|_{B=\beta}$ is concave. Bauer plus **two**
moment conditions — mass and $B$ — gives $\le2$ atoms. Union over $\beta$.
$\;\blacksquare$

**Verified with a control**, since I have twice mistaken a convexity direction:
worst $\Phi(\text{mid})-\text{chord}$ over 1955 random segments with $B$ held
fixed is $+4.10\times10^{-9}$ (concave); with $B$ free it is $-1.95\times10^{-1}$
(not concave). The control establishes the test can detect failure.

The binding obstruction uses exactly 2 orbits, so B‴ is **tight in orbit count**,
and the certification target is **5 parameters** $(p_1,q_1,p_2,q_2,w_1)$.

> A cheaper-looking route fails and is recorded: bounding $B\ge1-t$ *inside* the
> product to get the concave $J=\kappa_tL+\alpha\!\int\!c\,d\nu-(1-\alpha)\sigma^2$
> loses too much — $\min J=-0.0317$ at $t=\psi$, attained at $B=0.8065$, far from
> the bound $1-t=0.618$. The product must be handled by slicing, not bounding.

## Result 1 — the one-step barrier, exactly

[`onestep.py`](onestep.py) — verified symbolically and numerically.

$\psi$ is *exactly* the threshold of the iid argument, for an elementary reason:
$1-(2p-p^2)=(1-p)^2$ and $h(x)=h(1-x)$, so

$$h(2p-p^2)\ge h(p)\iff (1-p)^2\ge p\iff p\le\tfrac{3-\sqrt5}{2}.$$

Sawin's two branches meet at $1$ at $\psi$ by exact identities:
$2\psi-\psi^2=1-\psi$, and $\varphi(1-\psi)=1$ since $1-\psi=1/\varphi$.

Searching $k$-atom laws ($k\le4$) confirms the reach is exactly $\psi$, with
extremiser the **point mass** $\delta_t$: $\inf G=0$ at $t=\psi-10^{-4}$ and
$-1.55\times10^{-4}$ at $t=\psi+10^{-4}$.

## Result 2 — closed form for $c^*$

[`coupling.py`](coupling.py). Derived here, not transcribed:

$$\boxed{\;c^*=1-\frac{1-b}{2-h(b)},\qquad h(b)\big(2-h(b)\big)=h\big((1-b)^2\big)\;}$$

with $b=0.329454738503036972$ the larger root, giving
$c^*=0.382345533366702721$ — matching Cambie's reported value to $10^{-15}$.
The derivation: at the obstruction all three entropy terms coincide,
$(1-a)h((1-b)^2)=h(b)$ and $1-a=1/(2-h(b))$, and eliminating $a$ gives the
displayed equation. Liu's (7)–(10) contain the same relations.

**Mechanism.** The extremal law is $P_{SR}(1,b)=P_{SR}(b,1)=a$,
$P_{SR}(b,b)=1-2a$. The atom at $p=1$ is an **entropy sink**: the adversary pays
nothing there ($h(1)=0$) and *every* coupling gains nothing ($s^*(1,\cdot)=1$).
Negative correlation forces $(1,b)$ pairs, killing the coupled gain. Slack is
then exactly zero **for every $\alpha$** — and this one *is* a genuine universal,
not a sample: $Q=C=L$ at the obstruction, so
$F=(1-\alpha)L+\alpha L-L\equiv0$ identically in $\alpha$. Machine-checked at
$\alpha\in\{0,0.0356,0.2,0.5,1\}$ as corroboration of the identity, not as its
basis.

## Result 3 — Yu's reduction has an invalid step

[`concavity.py`](concavity.py), [`yu_gap.py`](yu_gap.py),
[`tests/test_yu_counterexample.py`](tests/test_yu_counterexample.py).

Yu reduces the problem to five parameters via Krein–Milman, resting on (verbatim):

> *"Note that $P_{pq}\mapsto g(P_{pq},\alpha)$ is concave, since by [4, Lemma 5],
> $P_p\mapsto \mathbb{E}_{(p,q)\sim P_p^{\otimes2}}h(p+q-pq)$ is concave, and
> $P_{pq}\mapsto P_p$ is linear."*

**Both the concavity claim and the inference it supports are false.**

*Concavity.* With $x=1-p$ the functional is $F(\mu)=\iint h(xy)\,d\mu\,d\mu$, and
concavity is equivalent to $h(xy)$ being **conditionally negative semidefinite**.
It is not, even on the feasible set $\{\mathbb{E}[p]\le t\}$ and even with
$\mathbb{E}[h(p)]\ge0.8$. Exact defect for
$\mu_1=(1-t)\delta_0+t\delta_1$, $\mu_2=(1-w)\delta_0+w\delta_{p_2}$:

$$\tfrac{w}{4}\Big[w\,h(2p_2-p_2^2)+2h(p_2)(t-w)\Big],$$

$=+0.03804$ at $w=0.531$, $p_2=0.2179$, $t=c^*$. A non-degenerate witness:
$\mu_1=\delta_{0.306787}$, $\mu_2=0.10457\delta_1+0.89543\delta_{0.310216}$,
defect $+0.002741$. Both hand-verified.

*The inference step.* For $P=\sum_i\gamma_iQ_i$ with $Q_i$ extreme, Yu asserts
$g(P,\alpha)\ge\sum_i\gamma_ig(Q_i,\alpha)$. Take

$$Q_1=\delta_{(3/20,\,3/20)},\qquad
Q_2=(1-\beta)\delta_{(1/10,\,1/10)}+\beta\,\delta_{(1,1)},\quad
\beta=\frac{t-1/10}{9/10},\qquad \gamma=\tfrac12 .$$

Both have exactly Yu's extreme-point shape ($\beta=0$ with $a=3/20\le t$;
$\beta>0$ with $a=1/10\le t<b=1$), and $P$ lies in $\mathcal{P}_B$. Then

$$g(P)=0.570142\ <\ 0.594662=\gamma g(Q_1)+(1-\gamma)g(Q_2),$$

a defect of $-0.02476$ at $t=c^*$ and $-0.02470$ at $t=\psi$. Verified at 40
digits and reproduced by hand.

**Scope — this refutes the proof, not the result.** On the same instance the
*ratio* conclusion Yu actually needs still holds
($1.2233\ge\min(1.3988,1.0447)$). And the conclusion appears true: an
unrestricted search over $k$-atom laws, $k\le7$, with the inner minimisation over
couplings solved **exactly by LP**, finds no violation at $t=c^*$, and the
5-parameter family attains the same infimum to $10^{-13}$
([`robust.py`](robust.py)). So $c^*$ is very likely the true threshold and
currently has no valid proof.

## Result 4 — no linearisation of the quadratic term

*The kernel is not PSD* ([`kernel.py`](kernel.py)). $h(xy)$ has
$\lambda_{\max}=+0.684$ and $15$ negative eigenvalues; explicit 2-point witness
$x\in\{0.3205,0.9993\}$. So the quadratic term admits **no linearisation**
$Q(\mu)\ge2B(\mu_0,\mu)-B(\mu_0,\mu_0)$. Exact series identity used:
$H(u)=u\ln\frac1u+u-\sum_{k\ge2}\frac{u^k}{k(k-1)}$, so
$H(xy)$ is rank-2 indefinite $+$ rank-1 PSD $-$ a positive combination of
$x^ky^k$.

## Result 5 — what dual certificates can and cannot do

[`dual_barrier.py`](dual_barrier.py).

> **Retraction.** An earlier version of this section claimed a theorem, *"linear
> dual certificates cannot exceed $\psi$"*, obtained by demanding copositivity of
> a kernel $Q$ and evaluating it at $\delta_0,\delta_1,\delta_s$. That argument
> **conflated two different requirements** and is wrong. Both horns fail, and
> the accompanying table was measuring $\lambda(t)\ge1$ at the single point
> $s=t$ — nothing to do with duality. That is the entire source of the spurious
> *"$\alpha=0$ recovers $\psi$"* coincidence, which I reported as a punchline.

Monge–Kantorovich duality on the transportation polytope gives, for any $\phi$ with

$$\textbf{(A)}\qquad \phi(s)+\phi(r)\le h\big(s^*(s,r)\big)\quad\forall s,r\in[0,1],$$

the linear bound $C(\mu)\ge2\int\phi\,d\mu$. What happens next depends on
**which cone copositivity is demanded on**, and the two readings behave completely
differently.

### (a) Global multiplier — the class is empty

Absorbing $\mathbb{E}_\mu[s]\le t$ into a multiplier $\gamma\ge0$ demands the
inequality for *all* probability measures, so $\delta_s$ is admissible for
**every** $s$. But $s^*(0,0)=0$ and $h(0)=0$ force $\phi(0)\le0$ via **(A)**,
and $Q(0,0)=2\alpha\phi(0)-\gamma t\ge0$ then forces $\gamma=0$. With $\phi\le0$,

$$Q(s,s)\le(1-\alpha)h(2s-s^2)-h(s),$$

which at $s=0.45$ equals $(1-\alpha)(0.88347)-0.99277<0$ for every $\alpha\ge0$.
**No certificate of this form exists at any $t>0$** — the route is vacuous, not
capped.

### (b) Decoupled certificates — genuinely capped at $\psi$

**Proposition.** *Suppose a certificate (i) bounds the iid term by Sawin's sharp
$I(\mu)\ge\lambda(t)\,\mathbb{E}_\mu h$, (ii) bounds the coupled term by a linear
functional $2\int\phi\,d\mu$ with $\phi$ obeying* **(A)**, *and (iii) discharges the
mean constraint with a multiplier. Then $t\le\psi$.*

*Proof.* After (i)–(ii) the bound is **linear** in $\mu$, so (iii) is lossless and
the requirement is $G(s)+\gamma(s-t)\ge0$ for all $s$, where
$G(s)=[(1-\alpha)\lambda(t)-1]h(s)+2\alpha\phi(s)$. At $s=0$: **(A)** gives
$\phi(0)\le0$ while the requirement gives $2\alpha\phi(0)\ge\gamma t$, so
$\gamma=0$. At $s=1$: $h(1)=0$ gives $\phi(1)\ge0$, and **(A)** at $(1,1)$ gives
$\phi(1)\le0$, so $\phi(1)=0$; since $s^*(s,1)=1$ for every $s$, **(A)** at $(s,1)$
gives $\phi\le0$ throughout. Then $G(t)\ge0$ forces
$\lambda(t)\ge\frac1{1-\alpha}\ge1$, and $\lambda$ is strictly decreasing with
$\lambda(\psi)=1$. $\qquad\blacksquare$

The multiplier is applied **only to a linear functional**, where Lagrangian
duality is tight — precisely what the retracted version got wrong. Since
$\lambda(c^*)=0.999385920316<1$, this route cannot reach $c^*$.

### (c) Mean-constrained copositivity — CLOSED (Theorem C, the sink pair)

The retraction left this open. It is not: it closes, by a mechanism far simpler
than either of the above.

**Theorem C (sink-pair barrier).** *Let $\alpha\in(0,1]$ and let $\phi$ satisfy*
**(A)**. *If*
$$(1-\alpha)Q(\mu)+2\alpha\!\int\!\phi\,d\mu-L(\mu)\ \ge\ 0\quad\text{for every }\mu
\text{ with }\textstyle\int p\,d\mu\le t,$$
*then $t\le\psi$ — even if $\phi$ is allowed to depend on $\int p\,d\mu$.*

*Proof.* Take the **sink mixture** $\nu_t=(1-t)\delta_0+t\delta_1$, feasible with
mean exactly $t$. Every entropy term vanishes on it: $h(0)=h(1)=0$, $0\oplus0=0$,
and $0\oplus1=1\oplus1=1$, so $Q(\nu_t)=L(\nu_t)=0$ and the hypothesis reduces to
$$2\alpha\big[(1-t)\phi(0)+t\phi(1)\big]\ge0.$$
But $s^*(0,0)=0$ and $s^*(1,1)=1$, both with $h=0$, so **(A)** forces
$\phi(0)\le0$ *and* $\phi(1)\le0$. A nonnegative combination of two nonpositive
numbers is $\ge0$ only if both vanish, so $\phi(0)=\phi(1)=0$. Since
$s^*(p,1)=1$ for **every** $p$, **(A)** at $(p,1)$ gives $\phi\le0$ throughout.
Then $\delta_t$ is feasible and forces $(1-\alpha)h(2t-t^2)\ge h(t)$, i.e.
$\lambda(t)\ge\frac1{1-\alpha}\ge1$, i.e. $t\le\psi$. $\qquad\blacksquare$

Both $\nu_t$ and $\delta_t$ have mean exactly $t$, so the argument never leaves
the slice of mean $t$ — which is why mean-dependent potentials do not help.

**The mechanism.** $\{0,1\}$ is a two-point set on which the *entire* entropy
functional vanishes, and $\nu_u$ is feasible for every $u\le t$. So the tight set
is a whole **curve**, and any linear surrogate for $C$ must be exact along all of
it — which pins it to zero at both endpoints, after which **(A)** drags it
nonpositive everywhere. Point masses then cap the iid term at $\psi$.

The same curve kills the cheaper idea of linearising the convex direction:
$G(\nu_u)=1-u$ sweeps $[1-t,1]$, so no single tangent to $G^2$ is exact on the
tight set. **The support reduction must freeze $G$, not linearise it** — which is
exactly what Theorem B′ does.

Together: every linear treatment of the coupled term caps at $\psi$, so a support
reduction is unavoidable, and Theorem B′ provides one — at the cost of carrying
the coupling as a variable, which is precisely what Theorem C says is mandatory.

## Claimed full proofs — all fail

Checked first-hand, none accepted, all preprint-only.

| paper | verdict |
|---|---|
| Demontis [2405.03731](https://arxiv.org/abs/2405.03731) | gap in Thm 3: *"As we suppose $Y$ has maximum cardinality we must conclude $Y\cup\{i\}=R$"* — maximality in $D-D_i$ does not constrain $R\in D_i$ |
| Agama [1711.02665](https://arxiv.org/abs/1711.02665) v6 | Thm 5.1 takes $l\to\infty$ for a fixed finite universe; Def 3.1's $n\to\infty$ is ill-defined |
| Scandone [2302.03484](https://arxiv.org/abs/2302.03484) | body not retrievable; `[UNRETRIEVED]` |
| Hachimori–Kashiwabara [2504.13454](https://arxiv.org/abs/2504.13454) | **legitimate but narrow**: formalises average rarity for *ideal* families only; the full statement is `theorem frankl_conjecture ... := sorry` |

## Files

| file | contents |
|---|---|
| `gate_b/` | Gate B exact definitions, Cartesian-power proof that $c_{\rm cl}^\star=\infty$, the $\Theta(n)$ growth corollary, standalone Arb / dyadic / exact-rational verifiers, machine-readable candidates and certificates, checkpointed $n=8$ search and exact tensorization audit, tests, experiment ledger, and paper draft |
| `entropy.py` | base-2 $h$, $\varphi$, $\psi$, Sawin's $\lambda$, one-step functional |
| `onestep.py` | $\psi$ exact (sympy) + extremiser search |
| `coupling.py` | closed form for $c^*$; obstruction; exact-LP coupling search |
| `robust.py` | unrestricted $k$-atom search; does the family attain the infimum? |
| `kernel.py` | $h(xy)$ indefinite; series identity |
| `concavity.py` | concavity fails on the feasible set |
| `concavity_nondeg.py` | it fails away from the degenerate corner too |
| `yu_gap.py` | counterexample on genuine extreme points |
| `dual_barrier.py` | dual-route triage: empty / capped / **Theorem C** (one claim retracted) |
| `decomposition.py` | Theorem A/B′ human derivations plus exact/sampled controls; not collector-executed |
| `reduction.py` | Theorem A′ human series proof plus exact/sampled controls |
| `support3.py` | pair-orbit landscape; exact obstruction family; margins |
| `margin_lemma.py` | **human-audited** Margin Lemma and historical sampled controls |
| `cert.py` | rigorous **Arb** enclosures + B&B; **no certificate obtained** |
| `relax_probe.py` | three candidate dimension reductions, **all refuted**; $q_2\to1$ asymptotics |
| `diag_v2.py` | corner bound `phi_corner` (point coefficient via mean-monotonicity); soundness-checked vs mpmath |
| `diag_resolution.py` | $\Phi_{\rm lo}$ vs FULL box width at the obstruction; flips positive at $\approx2\times10^{-4}$ |
| `diag_localize.py` | do corner-B&B residuals localise? budget-limited, unresolved |
| `diag_classify.py` | true $\Phi$ at residual centres by distance bucket (sampled, not certified) |
| `lemma_rh.py` | **conjecture (numerical)** $rh(x)\le2\min(x,1-x)$; candidate corner discharge rule |
| `tests/test_yu_counterexample.py` | 40-digit check of that counterexample |
| `lemma_rh_proof.py` | human endpoint/sign proof plus exact symbolic and Arb subchecks |
| `thmB3_proof.py` | **human-audited** support-reduction writeup with cited standard inputs; executable controls are sampled |
| `bridge_uc.py` | human universal bridge identities plus exact finite implementation schemas |
| `bound_kkt.py` | centered/MVT lower bounds, five-dimensional gradient enclosures |
| `diag_exhaust.py` | exact-rational mean contractor + corner bound consumed by `cert3.py` |
| `arbcore.py` | shared Arb helpers: $h$, $s^*$, $\sqrt{rh}$ enclosures, global $rh$ cover |
| `cert2.py` | w-eliminated 4-D certifier; certified $\rho$ cover; ratio rule; exact endpoint monotonicity proof |
| `cert3.py` | composite 5-parameter certifier: contractor → corner → ratio → centered → pinned face; **emits one-byte-per-node proof traces** |
| `cert3_replay.py` | independent DFS replay: re-derives every box, re-proves every claimed rule, rejects residuals/tampering |
| `cert3_par.py` | immutable campaign producer/supervisor/worker (v2): frozen snapshot, atomic no-replace commits, trace binding |
| `cert3_collect.py` | independent checker: pinned inventories, pristine-snapshot listing, hash gates, **mandatory replay of all 8 traces** |

## What would close the gap

**Exactly one step remains**, and it is a bounded computation with no structural
gap in front of it.

**Status (2026-08-12).** As posed, the computation did NOT go through in 5-D
interval arithmetic: near the obstruction the corner bound needs boxes of full
width $\approx2\times10^{-4}$ (a uniform grid at that cell size is
$\approx3\times10^{18}$ boxes), and wide boxes elsewhere lose to $L_{\rm lo}=0$
at atom-boxes touching $\{0,1\}$.  Adaptive B&B is NOT ruled out.  See
PROGRESS.md, session 4.

**The target.** By **Theorem B‴**, prove $\Lambda\ge0$ on pair-orbit measures with
$\le2$ orbits and $\sum w_i(p_i+q_i)/2\le t$: **five** parameters
$(p_1,q_1,p_2,q_2,w_1)$ with $p_i\le q_i$. All three earlier obstacles are gone:

| obstacle | status |
|---|---|
| inner coupling LP | gone — $\int c\,d\nu$ is **linear** in $\nu$ |
| sink degeneracy | gone — $\Lambda$ is **scale-free**, bounded by $\max\rho=0.3049467$ |
| slack at the binding case | gone — $\Lambda$ is **exact** when $\le1$ atom is outside $\{0,1\}$ |
| dimension | 8 → **5**, and tight in orbit count (the obstruction uses 2) |

**Interval arithmetic is straightforward here.** $p\oplus q$ and $s^*$ are monotone
in both arguments, so their ranges on a box are exact; $h$ is enclosed by
monotonicity about $1/2$; $\rho h=2(1-p)h(p)-h((1-p)^2)$ needs only $h$; and near
the sink corner, where $L$'s enclosure straddles $0$, use the proved bound
$\sigma^2/L\le\max_{\text{box}}\rho$ instead of dividing — which is small there
precisely because $\rho\to0$ at both sink points.

**Margins to certify** ([`margin_lemma.py`](margin_lemma.py)):

| $t$ | $t-\psi$ | $\min\Lambda$ |
|---|---|---|
| $0.3820660112501051$ | $10^{-4}$ | $4.555637\times10^{-4}$ |
| $0.3821660112501051$ | $2\times10^{-4}$ | $2.925862\times10^{-4}$ |
| $0.3822660112501052$ | $3\times10^{-4}$ | $1.296\times10^{-4}$ |

The margin is linear in $c^*-t$ with slope $1.6298$ to four digits, so $t$ trades
directly against branch-and-bound cost.

**Answered during this work:** B″'s 8 parameters were not optimal. Theorem B‴
reaches 5 and is tight in orbit count, since the binding obstruction uses exactly
2 orbits. Whether 5 can be cut further — e.g. a lemma forcing $q_2=1$, which holds
at the obstruction — remains open.


## The certification attempt — built, sound, and it does **not** converge

[`cert.py`](cert.py). Reported here because a failed attempt with numbers is worth
more than an untried plan.

> **Soundness retraction.** My first version of this file was **not** a rigorous
> certifier. It extracted Arb bounds with `.lower()`/`.upper()`, converted them to
> Python floats, and then did all the range algebra in binary floating point.
> Round-to-nearest can move a lower bound *up* or an upper bound *down*, so any
> "certificate" it produced would have been worthless — and the random-point
> validation I had written could not detect it, because the defect was in the
> arithmetic, not the formulas. Rewritten with every quantity an Arb ball end to
> end, and every pruning decision taken from a **certified** Arb sign.

What is verified:

- flint's comparisons are conservative: `arb(0,1) >= 0` is `False`, and
  `arb(0.5,0.5) >= 0` is `False` because $0$ lies in that ball. So clearing a box
  requires a certified sign, not a hopeful one.
- $\max_p rh(p)\le0.234318694$, proved by covering $[0,1]$ with 20000 Arb balls
  (the float search reports $0.234229480$, so the bound brackets it from above).
- Enclosures for $h$, $s^*$, $rh$ contain the mpmath truth on 1200 random
  intervals, 0 failures; and the $\Phi$ enclosure contains the dps-50 mpmath
  reference at test points to within $10^{-30}$ — at the obstruction the point-box
  enclosure reproduces $+3.83871\times10^{-4}$ exactly.
- $\Phi$ needs **no division**, which is what makes the sink corner reachable at
  all: $\sigma^2$ is quadratic in the non-sink mass while $L$ is linear.

**The result: 0 boxes cleared out of 60000 processed**, at both $t=\psi$ and
$t=\psi+10^{-4}$. No certificate. The bottleneck is enclosure quality in 5
dimensions against a $\approx4\times10^{-4}$ margin: naive ball arithmetic has
width $O(\mathrm{diam})$, so resolving the margin needs $\mathrm{diam}\sim10^{-4}$,
i.e. $\sim10^{8}$ boxes — out of reach for this implementation.

What would plausibly close it. **Three candidate reductions were checked and all
died** — see [`relax_probe.py`](relax_probe.py) and `PROGRESS.md` correction 7:
aggregate monotonicity does *not* give $M=t$ (the coordinates are coupled), the
Pareto set of the orbit surface is generically 2-D not 1-D, and $\Phi$ is not
concave in $w$ (the obstruction minimises at interior $w=0.8422$). No
single-variable derivative is one-signed on the feasible set. **Dimension remains
5 and the $\sim10^8$-box estimate stands.**

What survived is a boundary structure at $q_2\to1$, which is where the obstruction
sits: $rh(1-\varepsilon)=\varepsilon^2/\ln2+O(\varepsilon^3)$, so
$\sqrt{rh}$ is **linear** there, and the only singular quantity is
$\eta=h(q_2)\sim\varepsilon\log_2(1/\varepsilon)$, entering $\Phi$ **affinely with
positive coefficient**. So the ranked plan is:

1. **Split the domain at $q_2=1-\varepsilon_0$.** On the slab, use the
   $\varepsilon\log_2(1/\varepsilon)$-vs-$\varepsilon$ rate gap to pin $q_2=1$
   rigorously — 4 parameters, and $h(q_2)=\sqrt{rh(q_2)}=0$ kills terms.
   Verified on the feasible range $w\ge0.842245414678$: moving $q_2$ off $1$
   strictly raises $\Phi$, 0 violations. $\varepsilon_0$ is **not yet computed**
   and the $1/\ln2$ constant is numeric, not proved.
2. **Centred / mean-value forms on the bulk $q_2\le1-\varepsilon_0$**, where $h''$
   is bounded, for $O(\mathrm{diam}^2)$ instead of $O(\mathrm{diam})$. This is the
   only place that route is defensible — it cannot work *at* the obstruction,
   because $h'(q_2)\to-\infty$ there.
3. Compiled arithmetic rather than Python-level Arb calls.

The blocker on step 1 is that $\varepsilon_0$ must be explicit and the slab
argument must survive jointly in $(w,q_2)$, not just along the obstruction family.

Until this is carried out, $\psi=\frac{3-\sqrt5}{2}$ remains the largest constant
with a complete proof, and that fact is not reported as such anywhere in the
literature. What has changed across this work: the gap was a missing *idea* (a
support reduction, since Yu's is invalid), and it is now a bounded *computation*
with an explicit margin and a bound that is tight at the binding instance.
