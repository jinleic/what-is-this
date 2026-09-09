# The analytic input, regime by regime

**Owner of this file:** the audit of *external* theorems against this project's
parameters. Proofs and reductions live in [`THEOREMS.md`](THEOREMS.md); this file
never restates them, it only records which published statement covers which
regime and why. Verdicts are `APPLIES`, `PARTIAL`, or `FAILS BECAUSE …`; an input that a
cited theorem fully covers is marked `DISCHARGED` and its use is load-bearing.

Snapshot date 2026-08-28. Every entry was traced to arXiv or publisher text; two
entries are flagged as unverified at source and must not be relied on.

## 1. What is actually needed

Section 15 proves unconditionally that every class of modulus
$M\le\sqrt N/2$ — more precisely every progression of length at least
$P_0\approx2\sqrt{2\,\mathrm{top}}$ — contains a non-witness. Section 16.2 proves
the single-position route cannot pass $\sqrt{2N}$. So the open interval is

$$\sqrt{2N}\ \le\ M\ \le\ N,$$

above which a class meets the block at most once and prescribes rather than
Seven analytic inputs now appear: I--III concern the survivor/class route,
IV--VII the compensation-good run route.

| # | Input | Statement | Section |
|---|---|---|---|
| I | prime localisation | some prime $p\sim P$ has $ap^{-1}\bmod q$ in a prescribed arc | 15.6, hypothesis IP |
| II | survivor bound | $\Psi$-type upper bound $S(N,M,a)<N/M$ for an individual modulus | 16.1, (45) |
| III | window primes | no $L$ consecutive residues $\bmod M$ are prime-free in $[N,2N)$ | 16.3, (46) |
| IV | level congruence | primes in $[y,2y)$ hitting $\ge q^r2^{-r}$ prescribed residues $\bmod q^r$, $q^r\le y^{1-\delta}$ | 17.6, **DISCHARGED** by Vinogradov |
| V | consecutive smooth pair | positive density with both $P^+(n),P^+(n+1)<x^{1/2}$ | 19.1, **DISCHARGED** by Pascadi--Yang transfer |
| VI | selected-prime $q$-adic mask | the three centered errors (73) are $o(x)$ | 20.4, **OPEN** |
| VII | mixed exceptional-spectrum large sieve | control coefficients $hk\lambda-h'k'\lambda'$ generated after Cauchy | 20.5, **OPEN; explicitly absent in Pascadi** |

III is weaker than I, and II is independent of both. The census in
`data/survivor_census.json` shows II and III holding with room to spare at
$N=10^7$; neither is a theorem.

## 2. Input I — inverse-prime sums (Kloosterman with primes)

**2026-09-08 update:** the two missing reductions now have a written candidate
for $9\le m\le C\sqrt N$, $M\le N^{17/32-\varepsilon}$, with fixed
$C,\varepsilon>0$ and sufficiently large $N$. It is **review-pending**:
the required Opus verdict was not obtained, so no theorem is certified here.
The unrestricted all-class prime-spike deduction is not valid.

The relevant modulus is $q=M/d_i$, not automatically $M$. Moving intervals
are replaced by a fixed rectangle on $P<p\le3P/2$; nonprimitive frequencies
are handled by truncation and effective-modulus reduction, not uniform
all-frequency cancellation. The table retains earlier primitive-sum
comparisons; their exponent balances alone are not every-class theorems.

| Bound | Hypotheses | Saving at $P\approx q$ | Verdict |
|---|---|---|---|
| Bourgain, *GAFA* 15 (2005) 1–34, explicit $\delta=0.0005\varepsilon^4$ by Baker, *Acta Arith.* 156 (2012) 351–372 | **$q$ prime**, $(ab,q)=1$, $X\ge q^{1/2+\varepsilon}$ | $q^{-\delta}$, $\delta\sim10^{-3}\varepsilon^4$ | APPLIES for prime $M$ throughout $[\sqrt{2N},N]$; yields $\theta\le\frac12+O(\delta)$ — a genuine but numerically minute crossing |
| Fouvry–Shparlinski, *Acta Arith.* 150 (2011) 285–314, Thm 3.1 | $q\ge2$, $(A,q)=1$, $X^{3/4}\le q\le X^{4/3}$; prefix convention in §3.1 | $(X^{15/16}+X^{2/3}q^{1/4})q^\rho$, every fixed $\rho>0$ | **SOURCE CHECKED 2026-09-08**; supports the restricted $17/32$ candidate below. Independent proof review remains pending. |
| Korolev, eq. (7) of arXiv:1911.09981 | any composite $q$, $X\ge q^{7/10+\varepsilon}$ | $X^{32/37}q^{7/74+\varepsilon}$, i.e. $q^{-3/74}$ at $X\approx q$ | APPLIES for composite $M$; true fixed point $\theta\le14/27=0.5185$ (the recorded $37/71=0.521$ evaluates at $X=q$) — dominated by Fouvry–Shparlinski above |
| Korolev, *Res. Number Theory* 6:24 (2020) Thm 1 | any $q$, $(ab,q)=1$, $q^{3/4+\varepsilon}\le X\le(q/2)^{3/2}$, $\Lambda$-weighted, inhomogeneous | $(q^{3/4}/X)^{1/7}$ or $(q^{2/3}/X)^{3/35}$ | APPLIES with an exact printed saving; at the floor length dominated by Fouvry–Shparlinski above |
| Changa–Korolev, *Math. Notes* 108:1 (2020) 87–93, DOI 10.4213/mzm12693, Thm 1 | $(a,q)=1$, $0<\delta_0<1/100$, $0<\varepsilon\le1/6$, $q^{1/2+\varepsilon}\le X\le q^{3/2}$ | $T_q(X)\ll Xq^{\delta_0}\Delta$; piecewise saving recorded below | Primitive input source-checked 2026-08-29, not consumed by the current F–S-only candidate. Its wider range may help smaller windows but does not independently certify eviction. |
| Kowalski–Michel–Sawin, *Ann. of Math.* 186 (2017) 413–500 | **$q$ prime**, general bilinear forms, $M=N>q^{3/7}$ | power | Engine, not directly our shape; §1.5.2 states composite moduli open |
| Bourgain–Garaev, *Izv. Math.* 78 (2014) 656–707 | **$q$ prime**, multilinear, $\prod|I_j|>q^{1/2+\varepsilon}$ | power | Engine; also records that one published short-sum claim has a proof "in doubt" |

**The Changa–Korolev saving in full (read at source, 2026-08-29).** In
$T_q(X)\ll Xq^{\delta_0}\Delta$ the envelope $\Delta$ is piecewise in $X/q$,
with $\alpha=(29-5\sqrt{33})/8=0.034648\ldots$ and $\beta=0.649625\ldots$:

$$
\Delta=\begin{cases}
(Xq^{-1/2})^{-\alpha}, & q^{1/2+\varepsilon}\le X\le q^{\beta};\\
(Xq^{-5/8})^{-4/19}, & q^{\beta}\le X\le q^{21/26};\\
(Xq^{-1/2})^{-1/8}, & q^{21/26}\le X\le q^{9/10};\\
(Xq^{-3/4})^{-1/3}, & q^{9/10}\le X\le q^{12/13};\\
X^{-1/16}, & q^{12/13}\le X\le q^{3/2}.
\end{cases}
$$

Writing the net saving as $q^{-s(X)}$ with, in the five pieces,
$s=\alpha\log_qX-\tfrac\alpha2$, $s=\tfrac4{19}\log_qX-\tfrac5{38}$,
$s=\tfrac18\log_qX-\tfrac1{16}$, $s=\tfrac13\log_qX-\tfrac14$, and
$s=\tfrac1{16}\log_qX$: $s$ is continuous at all four interior breakpoints —
exact identity at $\beta$, $21/26$, $9/10$, $12/13$ (which pins
$\beta=(95\sqrt{33}-511)/(2(95\sqrt{33}-519))=0.6496250\ldots$) — validating
the reading against the PDF text layer's $q^{-5/8-4/19}$ bracket-flattening.
Every piece has slope $<1$, so $s$ increases in $\log_qX$ throughout
$[q^{1/2+\varepsilon},q^{3/2}]$: at $X=q$, $s=1/16$; the maximum sits at the
range top $X=q^{3/2}$ with $s=3/32$.

**Floor-band admissibility.** At $X\asymp_C\sqrt N$, F–S requires
$q_i\gtrsim_C N^{3/8}$ from $X^{3/4}\le q_i$; C–K requires
$q_i\gtrsim_C N^{1/3}$ from $X\le q_i^{3/2}$. Below their respective
caps these inputs do not cover cofactor-dominated classes, a previously
identified failure regime. The candidate's $q\gg_C N^{2/5}$ clears the
F–S cap by $2/5-3/8=1/40$, explaining the growth of its cutoff $H$.

**Current source-to-target reduction.** Fouvry--Shparlinski's
[primary PDF](https://www.impan.pl/shop/publication/transaction/download/product/82801)
states Theorem 3.1 on printed p.297; §3.1, printed p.296, explicitly permits
prefix sums with the same bound up to a constant. Subtracting prefixes at
$P$ and $3P/2$ is therefore allowed when both range premises hold.

The previously missing gaps are now explicit:

1. **Moving interval:** use the fixed real interval
   $[(N+i)/(dP),\,2(2N+i)/(3dP))$, of width $(N-i)/(3dP)$.
   It lies inside every required $J_p$ on the shorter prime band. At
   $i=o(N)$ its main mass is $N/(6M\log P)$, with no moving Fourier weight.
2. **Nonprimitive frequencies and admissible positions:** the candidate
   chooses $H=\lfloor q/(3P/2)^{3/4}\rfloor$. For $h\le H$,
   $r=q/(h,q)\ge(3P/2)^{3/4}$ and the numerator is primitive modulo $r$.
   The discarded frequencies are paid for by the discrepancy term
   $\Pi/(H+1)$, not assumed to cancel. Five consecutive positions give
   $\prod d_i\mid24M$, selecting $d\le(24M)^{1/5}$; such five positions
   exist when $m\ge9$.

For $P\asymp_C\sqrt N$, outside the elementary small-$M$ range, this gives
$q\gg_C N^{2/5}$ and $H\gg_C N^{1/40}$. The candidate's relative errors,
apart from logarithms and a sufficiently small source loss $q^\rho$, have
powers
$-3/160-\varepsilon/5$, $-\varepsilon$, and
$-1/384-5\varepsilon/4$. Both endpoints of the published range and the
excluded-prime correction are checked in the written derivation.

**Conclusion for I.** The candidate covers arbitrary moduli and all classes
only with the stated $9\le m\le C\sqrt N$ restriction. The Opus review
attempts timed out, so it remains unaccepted. The former assertion that
closing two technical gaps automatically certifies every $m$ is withdrawn:
even $q,h=q/2$ disproves uniform all-frequency saving, and $m=2$,
$d=M=2^s$, $q=1$ gives an asymptotically empty floor-prime route at
$1/2<\theta<17/32$. Neither obstruction disproves class eviction by
another mechanism or Erdős #389.

Evidence, post-closure integrity disclosure and the exact remaining review gate:
[cycle reconciliation](../../docs/knowledge-system/institute/cycles/20260908T155849Z_f7f167/reconciliation.md).


## 3. Input II — survivors (smooth numbers in one progression)

Target: $\Psi(x,y;q,a)$ with $x=2N$, $y=\sqrt{4N}$, so $u=\log x/\log y\to2$,
and $q\in[\sqrt{2N},x^{2/3}]$, i.e. $q\approx y^{1+2\delta}$. Needed:
$S<N/M$, i.e. an effective constant below $1$ where the Dickman heuristic gives
$\rho(2)=0.30685$.

| Bound | Hypotheses | Verdict |
|---|---|---|
| Balog–Pomerance, *Proc. AMS* 115 (1992) 33–43 (as rendered in Granville's survey, eq. 4.7): $\Psi(x,y;a,q)\ll\Psi_q(x,y)/\varphi(q)$ | $y\ge q^{3/4+\varepsilon}$, $x\ge y^2$, individual $q$, no shape condition | **APPLIES to our exact corner** for every $q\le x^{2/3-\varepsilon}$, right direction, no averaging. FAILS TO CLOSE only because $\ll$ hides an unspecified absolute constant |
| Shiu, *J. reine angew. Math.* 313 (1980) 161–170 (Brun–Titchmarsh for multiplicative $f$) | $x^\kappa\le y\le x$, $k<y^{1-\alpha}$, $(a,k)=1$, individual $k$ | **APPLIES** for every $M\le N^{1-\alpha}$ on the interval $[N,2N)$; FAILS TO CLOSE on the same absolute constant. $u=2$ is its worst corner |
| Wright, arXiv:2508.17217v2 (2025), Thm 3.1: parts 1–2 for bounded / poly-log-growing $f\in\mathcal M_Q$ ($Q$-smooth-supported), $C_1=\alpha\kappa/41\le1/164$, $C_2=\frac5{656}\bigl(\frac{\alpha\kappa(1-\alpha\kappa/20)}{20}\bigr)$; plus a commented-out indicator corollary with $\rho$-exponent $1+C_1(1/\eta-1)+o(1)$ | as Shiu; **the printed part-1 display is false as printed**: it omits the $\exp(\sum_{p\le Q}f(p)/p)$ its own proof carries (proof line (563) ends $\cdots\exp(\sum_{p\le Q}f(p)/p)$, then the display drops it). Counterexample at $u=2$: $f=\mathbf 1_{Q\text{-smooth}}$, $k=1$, $y\approx x$, $Q=\sqrt x$ satisfies every hypothesis, true window count $0.3069\,y$, bound $\rho(2)^{C_1}y\log\log Q/\log x=o(y)$. Same defect in the v1/v2 abstracts and PDF p. 6 | **CHECKED at source (2026-08-29): PARTIAL.** Proof-faithful part 1 at $k=M$, $y=\sqrt{4N}$, $u=2$, $Q=\sqrt{4N}$: margin $\frac NM\cdot\rho(2)^{C_1}\cdot\tfrac12 e^{B_1}(1+o(1))=\frac NM\cdot0.649\ldots$ against a true $\frac NM\cdot\rho(2)=0.307\,\frac NM$ — it *halves* Shiu's margin ($e^{B_1}=1.299$, i.e. $4.23\times$ truth, because the prime sum truncates at $Q=\sqrt x$) but $(45)$ is still unclosed: margin $2.11\times$ truth, and the constant remains the same ineffective $\ll$. Wright's independent additions: part 2 handles unbounded $f(p)$ (Shiu cannot), and the commented-out indicator corollary is the only published-source mechanism with the *right* shape for (45) — $\Psi$-per-class $\ll y\,\rho(u_y)^{1+o(1)}$. Supersedes Balog–Pomerance only as the shortest effectivization route; the concrete task is to un-comment, complete, and effectivize that corollary |
| Granville, *Acta Math.* 170 (1993) 255–273 | lower bound $\Psi\gg x/q$ for $y\ge q^{3/4+\varepsilon}$; asymptotic only for $y\ge q^{1+\varepsilon}$ **and** $\log x/\log q\to\infty$ | Lower bound is the wrong direction; asymptotic FAILS because our $\log x/\log q\approx2$ |
| Harman, *Acta Arith.* 91 (1999) 279–289 | $x^\varepsilon<q<x^{1/2-\varepsilon}$ | FAILS: needs $x>q^{2+2\varepsilon}$, we have $x<q^2$ for every $q>\sqrt{2N}$ |
| Soundararajan, arXiv:0707.0299; Harper, *JNT* 132 (2012) 182–199; Banks–Shparlinski, *BLMS* (2022) | $q\le y^{4\sqrt e-\varepsilon}$ and $\log x/\log q\to\infty$ | FAIL on the limit: our $\log x/\log q\approx2$. Dartyge's survey states explicitly that this condition "doesn't cover the case $q\approx x^\alpha$ even for very small $\alpha$" |
| Fouvry–Tenenbaum 1991/1996; Drappeau, *Compos. Math.* 151 (2015); Drappeau–Granville–Shao, *Mathematika* 63 (2017); Pascadi, *Compos. Math.* 161 (2025) and arXiv:2505.00653 | averages over $q\le x^{3/5}$, $x^{66/107}$, $x^{5/8}$; and $y\le x^{1/c}$ | FAIL twice: averaged over moduli, and their $y$-range excludes $y=x^{1/2}$. A single class can be empty without violating any of these means |

**Why no structural obstruction exists here.** The only known mechanism that
breaks per-class upper bounds for smooth numbers is concentration into a
subgroup, which requires $y$ below the least non-residue; the record there is
$q^{1/(4\sqrt e)}=q^{0.152}$. Our $y\approx q^{1+2\delta}$ sits far above it, so
the primes up to $y$ generate $(\mathbb Z/M)^*$ and the subgroup mechanism
cannot operate. Combined with the $3.26$-fold margin, this is why the honest
description of II is *one effective constant*: Shiu's margin is $e^{B_1}=1.299$
per $\rho(2)$ unit ($4.23\times$ truth at $u=2$), Wright's proof-faithful part 1
halves it to $2.11\times$ by truncating the prime sum at $Q=\sqrt x$, and the
only shape-correct mechanism on record — Wright's commented-out indicator
corollary, $\Psi\ll y\,\rho(u_y)^{1+o(1)}$ per class — is exactly one effective
constant away from (45). Un-comment it, complete it, make it effective.


### 3.1 The commented-out indicator corollary at $u=2$ — REFUTED-FOR-PURPOSE (CHECKED at source and certified, 2026-08-29)

**Source.** arXiv:2508.17217v2 comment block, source lines 612–689 of
`/tmp/wrt/shiu_for_journal.tex`: parent commented smooth-count theorem
(l. 630–660, bound $\frac{y\rho(u_{x/Z})^{C_1}\rho(u_{Z/Q})\log^2Q}{\phi(k)\log x}$
with $Z=y/k^{1+\nu}$, $\nu<\kappa/2$) and corollary (Hilcor, l. 663–689), both
commented out. Conventions verified at source: $u=\log x/\log Q$ (l. 104);
$C_1=\alpha\kappa/41$ with $0<\alpha<\frac12$, $0<\kappa<\frac12$ (l. 58–62,
137–140; $C_1\le1/164$); $u_x$, $u_{x/Z}$, $u_{Z/Q}$ all abbreviate
$\log(\cdot)/\log Q$ with $x=x_sQ^u$-type products scaled by $Q$
(l. 639, 653; cf. $u_x$ in the l. 150–155 commented duplicate).

**The corner (correcting the §10 item-3 wording).** Input II has
$x=2N$, $Q=y=\sqrt{4N}$, so $\eta=\log y/\log x\to\frac12$, $u_x\to2$ — not
$\eta\to1$ — and $u_y=\log y/\log Q=1$. The smooth-count theorem forces
$Z=y/k^{1+\nu}\le y=Q$; with $k=1$ the two split indices are
$u_{x/Z}=\frac{\log(x/Z)}{\log Q}=\frac{\log(x/y)}{\log Q}=u_x-1\to1$ and
$u_{Z/Q}=\frac{\log(Z)}{\log Q}=u_y=1$: **both Dickman factors of the parent
bound sit at $\rho(1)=1$** before any corollary algebra begins.

**Findings at the corner (all certified in exact rationals and 300-bit Arb
balls by `e389/wright_indicator_corollary_margin.py`; byte-stable stdout;
run `./.venv/bin/python -I -B e389/wright_indicator_corollary_margin.py`
from `math/`).**

1. **Universal vacuity of form 1.** The printed exponent
   $C_1-(1-C_1)\eta=-\frac{161}{328}+o(1)<0$ at $\eta=\frac12$, and the
   factor $\rho(u_x)^{C_1-(1-C_1)\eta}$ is an upper bound over all $\eta$:
   its exponent ranges in $[-1+2C_1,0]$ and its base $\rho(2(1-\eta))$ in
   $[\rho(2),1)$ (using $\rho(v)=1$ on $[0,1]$, $\rho(v)=1-\ln v$ on $[1,2]$
   — both classical and exact), so the factor's supremum over $\eta$ is
   $\rho(2)^{2C_1-1}=3.21227658592331260215\ldots>1$: **form 1 is vacuous
   everywhere in the corner's $\eta$-range, certified, not merely at
   $\eta=\frac12$** (there it reads $\rho(2)^{-161/328}=1.7858\ldots$).
2. **Zero margin of form 2.** At the corner $u_y=\log y/\log Q=1$ and
   $\rho(1)=1$, so $\rho(u_y)^{1+C_1(1/\eta-1)+o(1)}=1$: form 2 reduces
   exactly to the class bound $y$, margin exactly zero — it never sees the
   Dickman gain, which lives at $u_x=2$ (the $\rho(2)$-exponent slot),
   not at $u_y=1$.
3. **The parent bound is elementarily vacuous at $Q=y$.** The $c$-split's
   dyadic-sum step evaluates $\sum_{c\le Q,\,Q\text{-smooth}}1/c\ge\sum_{c\le Q}1/c
   =H_Q=\log Q+\gamma$ (every $c\le Q$ is $Q$-smooth at $Q=y$), so the
   parent bound per unit $y$ is $\frac{\log^2Q}{\log x}
   =\frac12\log Q\,(1+o(1))$ — growing, with no constant to drive below $1$.
   This is a *split-of-unity washout*: forcing $Z\le y=Q$ makes the two
   $\rho$-factors $\rho(u_{x/Z})^{C_1}\rho(u_{Z/Q})=\rho(1)^{1+C_1}=1$
   identically in $\nu$.
4. **Printed exponent algebra slip.** Carrying the honest split chain through
   the proof's own $u^{-u(1+o(1))}$ heuristic at $\eta=\frac12$ gives total
   exponent $C_1(1-\eta)+\eta=\frac{165}{328}$ on $u^{-u}$-type factors —
   the printed chain $C_1-(1-C_1)\eta$ differs by exactly
   $2\eta(1-C_1)=\frac{163}{164}$ (the $\eta$-power of the second
   $\rho$-factor converted to a $(1-\eta)$-power and combined with the first
   instead of being dropped), and at $u=2$ the honest chain carries the
   split-origin prefactor $\eta^{-\eta u}(1-\eta)^{-C_1(1-\eta)u}
   \ge2^{(1+C_1)u/2}$, which the printed proof drops: at $\eta=\frac12$ the
   honest bound at $u=2$ is $2^{1+C_1}\rho(2)^{(1+C_1)/2}=1.1085782816\ldots y>
   $\frac{y}{2\sqrt2}$ — still vacuous, but by a *bounded* factor, and the honest
   exponent is *positive* ($165/328$), unlike the printed negative one.
5. **The mechanism cannot reach our corner at any saturation.** Effectivity
   of the honest repaired shape $2^{(1+C_1)u/2}\rho(u)^{(1+C_1)/2}$ is
   $C_1$-independent and reads $\rho(u)=2^{-u}$, whose unique root
   $u^\ast\in(2.21,2.22)$ is certified by rigorous midpoint quadrature of
   $\rho(u)=\rho(2)-\int_2^u\frac{1-\ln(t-1)}{t}dt$ with cell-wise
   second-derivative majorants (break-even first at
   $u=u^\ast\approx2.2124705858847$; and by $u=3$ the ratio has fallen to
   $\rho(3)/2^{-3}\le0.048608388/0.125<\frac12$, certified) — **the corner
   $u=2$ sits strictly below the break-even point; the repaired shape is
   vacuous on all of $[2,3]$ and only non-vacuous at $u>2.2125$ there.**
6. **Effectivity gap, quantified both ways.** To kill the class bound
   $N/M=y/(2\sqrt2)$ at $M=\sqrt{2N}$ the total smooth count needs a
   $\rho(2)$-exponent $e_{\rm req}=
   \frac{\ln(2\sqrt2)}{\ln(1/(1-\ln2))}=0.880084778655420\ldots$; the honest
   mechanism ceiling at the corner is $(1+C_1)/2=\frac{165}{328}
   =0.503048780487805\ldots$; shortfall $0.377035998167615\ldots$ in
   $\rho(2)$-exponent. Anchored instead to the best Wright-faithful bound
   (corrected Thm 3.1 part 1 with its restored $\exp(\sum_{p\le Q}1/p)$,
   f = indicator of $\mathcal S(Q)$: $e^{B_1}\rho(2)^{C_1}\frac{\log Q}{\log x}
   =0.6447752017\ldots$ per $y$, $2.101\times$ truth — and
   $1.8237\times$ the class bound), the gap closes to
   $e_{\rm need}=0.514709211536745\ldots$ vs ceiling $\frac{165}{328}$:
   **$0.0117$ $\rho(2)$-exponent units** and a multiplier of
   $e_{\rm need}/\frac{165}{328}=1.0232\times$ the honest ceiling
   (vs $1.7495\times$ against $e_{\rm req}$): any effective completion of
   Wright's method needs a $\rho(2)$-exponent at the $u_x=2$ slot
   $2.3\%$ larger than its own $\eta=\frac12$ ceiling supplies — the $u_y$
   slot reads $\rho(1)=1$ identically and can contribute nothing.
7. **Cross-checks.** printed $C_1=\alpha\kappa/41$, $\alpha,\kappa<\frac12$
   gives $C_1\le1/164$ (exact); $e^{B_1}=1.29887332140903\ldots$ matches the
   §3 record $e^{B_1}=1.299$; the part-1 ratio $2.1013\times$ truth matches
   the §3 record $2.11\times$ up to the certified refinement
   $\rho(2)^{C_1}=0.9928$; Changa–Korolev constants re-verified at 300-bit
   Arb ($\alpha=(29-5\sqrt{33})/8$, $\beta=(95\sqrt{33}-511)/(2(95\sqrt{33}-519))$,
   identity $\alpha(\beta-\frac12)=\frac4{19}(\beta-\frac58)$ to $<10^{-60}$,
   $\beta$ re-derivation via $\beta=(\alpha/2-\frac5{38})/(\alpha-\frac4{19})$
   likewise).

**Repair path for §10 item 3.** What needs proving for (45) at the corner is
an effective $\rho(2)$-exponent $e>0.8801$ ON A $\rho(u_x)$ SLOT, or
equivalently rescuing the mechanism above the break-even $u^\ast\approx2.2125$
(which our $u=2$ corner cannot supply — $u_x\to2$ is the Input-II identity,
not a free parameter). The corollary's own instrument cannot be repaired
either way: the $u_y$-slot reads $\rho(1)=1$ identically, the $u_x$-slot
carries the negative printed exponent, and the parent $c$-split's
$H_Q$-cost is structural at $Q=y$. A repaired argument must change the
$c$-split at $\eta=\frac12$ (e.g. accumulate the $c$-sum with its own
$\rho$-gain rather than against $H_Q$, or re-split $n=cd$ with
$Z<Q$-smooth decomposition) — a genuinely new input, not a comment-block
completion. The §10 item-3 phrasing "complete the $\eta\to1$ endpoint"
misdescribes the corner: the corollary's printed proof line "We will assume
that $\eta$ is a fixed constant less than 1, since the case of $\eta=1+o(1)$
follows from Hildebrand's result" shows the corollary never reaches even its
own intended endpoint. $u=2$ needs an $\eta=\frac12$ theorem, which Wright's
comment block does not contain even implicitly.

**Verdict: REFUTED-for-purpose (checked, certified, and closed against the
only two readings of the printed claims).** The corollary gives at $u=2$
$\eta=\frac12$ either $3.21y$/$1.79y$ (form 1, vacuous, universally in
$\eta$) or exactly $y$ (form 2, zero margin, $\rho(1)=1$), against the
needed effective constant below $N/M=y/(2\sqrt2)=0.35355\ldots y$. (45)
remains **exactly one effective constant** away — one $3.26$-faceted margin
of room, with **no** shape-correct mechanism on record: the §3
"Un-comment it, complete it, make it effective" instruction is hereby
**withdrawn**; Wright's comment block cannot be completed to our corner
effectively because its $c$-split loses the Dickman gain before the
corollary's exponents are even formed.


### 3.2 Shiu's own per-class theorem at $u=2$ — the constant is now metered: **NO at every depth** (CHECKED at source pp. 161–169, certified, 2026-08-29)

**Source.** Full text of P. Shiu, *J. reine angew. Math.* **313** (1980)
161–170, read this session from GDZ page images (iiif
`PPN243919689_0313:00000165..173` = printed pp. 161–169).  Wright's tex
quotes Shiu's Lemma 2 as $\Phi(x,y,z;k,a)\ll\frac{y}{\phi(k)\log z}+z^2$
(tex l. 361); the printed Lemma 2 (p. 164) is an **honest inequality with
constant 1**: $\Phi(x,y,z;k,a)\le\frac{y}{\phi(k)\log z}+z^2$, $z\ge2$,
$k<y\le x$ — Selberg's upper sieve per a single residue class (Halberstam–
Richert p. 104).  Notation map verified on the scans: Shiu's
$(\alpha,\beta)$ = (modulus range $k<y^{1-\alpha}$, length range
$x^\beta<y\le x$) correspond to Wright's $(\alpha,\kappa)$.  For
$f=\mathbf 1_{\mathcal S(Q)}$: $A_1=A_2=1$ exactly ($f(p^l)=1$, $f(n)\le1$).

**The corner, and a coverage correction.**  $x=2N$, $y=Q=\sqrt{4N}$,
$k=M$.  Shiu's hypothesis (2.2) ($k<y^{1-\alpha}$, *any fixed*
$\alpha>0$) forces $M<y=2\sqrt N$: **the theorem is silent on the upper
half of the (45) band** $M\in[\sqrt{2N},(2N)^{2/3}]$, i.e. everything past
$M=2\sqrt N$ — including $\theta>1/2$ — has no Shiu theorem at all.  What
follows meters the covered half $M\in[\sqrt{2N},\,2\sqrt N)$.

**The per-class meter (all certified by
`e389/shiu_effective_margin.py`, fmpq + Arb ≥ 300-bit, byte-stable;
run `./.venv/bin/python -I -B e389/shiu_effective_margin.py` from `math/`).**
Shiu's proof drives $z=y^{\alpha/10}$ ((5.1)); with the maximal legal
$\alpha_{\max}=\log(y/M)/\log y$ from (2.2), the sieve half-level obeys
$$\log z'=\frac{1}{2}\log z=\frac{\log(y/M)}{20},$$
which collapses across the band: at the frontier $M=\sqrt{2N}$ it is
$\frac{\log 2}{40}=0.0173286795128\ldots$.  Amortizing the class-I window
bound $y\,B(z,M)/(\phi(M)\log z)+z'^{3/2}M/y$ over $[N,2N)$
(Hildebrand's $\psi$-subadditivity, junk additive-nonnegative; no $o(1)$
decay claim is made at finite $N$) gives the per-class meter, with the
$b$-sum $B(z,M)=\sum_{b\le z,\,(b,M)=1}f(b)/b$; in the band $z\le2^{1/20}<2$
so only $b=1$ contributes and $B=1$ **exactly**:
$$R_{\min}=\frac{M/\phi(M)}{\log z'}=\frac{(40/\log2)\,(M/\phi(M))}
{\log_2(4N/M^2)}\;\xlongequal{\,M=\sqrt{2N}\,}\;57.7078016\ldots\,(M/\phi(M)).$$
Since $H_z\ge\log z$, **$R\ge2(M/\phi(M))\ge2$ universally** — via the
Lemma-4-form class-IV sum $\sum_{r\ge2}(r+1)e^{-r\log r/10}=23.03\ldots$
(form (5.8)) or its direct $\frac12\log2$-per-$r$ cousin (both certified);
and via $(2.3)$/Lemma 3: $\Lambda_3\le1$ by exact telescoping
$\sum_{n\ge2}\frac1{n(n-1)}=1$, $A_3\le2(\zeta(3/2)-1)=3.22475\ldots$.

**Two further floors, both $\gg1$.**  (a) *Balanced sieve*: re-picking $z$
to balance Lemma 2's $z'^2$ junk against the main term gives at best
$\log z'\le\frac12\log(y/M)+O(\log\log)$, hence
$R\ge\frac{4}{\log2}\,(M/\phi(M))=5.7707802\ldots(M/\phi(M))$ at the
frontier — oracle-level Selberg depth still cannot reach $1$.  (b) *Idealized
depth*: at $z'\to\infty$ the $b$-sum's exact Möbius expansion
$B(z,M)=\sum_{d\mid M}\frac{\mu(d)}d H_{\lfloor z/d\rfloor}
=\frac{\phi(M)}M\bigl(\log z+\gamma+\sum_{p\mid M}\frac{\log p}{p-1}\bigr)+o(1)$
(certified against a 327 379-term exact-rational direct sum for
$M=10^{40}$, $z'=327379$: **exactly equal**, $B=5.7486\ldots$, and matching
the asymptotic to $10^{-6}$) shows the cascade's realized constant tends to
$\mathbf{2\,(M/\phi(M))\ge2}$: the depth machinery carries the
$(\phi(M)/M)$ density as a second factor and **never returns $R<1$ at any
sieve depth**.  The ledger's "Shiu margin $e^{B_1}=1.299$ per $\rho(2)$
unit ($4.23\times$ truth)" survives as the exp-form normalization only; the
proof-faithful constant at idealized depth is $2(M/\phi(M))$, and the
finite-$N$ realized margin is $57.7\,(M/\phi(M))$.

**Three-way constant split.**  (i) *Printed-literal*: Lemma 2's $1$ and
$z^2$; Lemma 1's $3$ and $2$ (incl. (4.1)'s $2y/\log^2y$); Lemma 4's
$-\frac1{10}r\log r$ with $2A_1r^{1/4}$ and $r\le\frac{\log z}{4\log\log z}$;
(5.1) $z=y^{\alpha/10}$; (5.3)–(5.8) factors $2$ (log-split), $z^{-1/4}$,
$z^{-1/8}$, $f(n)\le n^{\alpha\beta/80}$, $r_0=\bigl[\frac{\log z}
{\log(\log x\log\log x)}\bigr]$, $A_5=A_1^{20/(\alpha\beta)}$; p. 169
assembly $z^2<\frac{y}{\phi(k)\log z}$.  (ii) *Derivable, Arb-certified*:
$A_1=A_2=1$; $\Lambda_3\le1$ (exact); $A_3\le2(\zeta(3/2)-1)$; band
$b$-sum $=1$ (exact); class-IV sums both forms; $\alpha_{\max}$,
$\log z'_{\max}=\frac{\log(y/M)}{20}$; the R-meter with floors
$57.708\,(M/\phi)$ [Shiu depth], $5.7708\,(M/\phi)$ [balanced],
$2\,(M/\phi)$ [idealized].  (iii) *Genuinely ineffective*: **none** — every
$\ll$ and $O(1)$ in the $f=\mathbf 1_{\mathcal S(Q)}$ chain realizes
explicitly; the single external link, (4.1)'s PNT form
$\sum_{p\le y}\frac1{\log p}\le\frac{2y}{\log^2y}$, is effective through
published explicit bounds (Rosser–Schoenfeld / Dusart) and is labeled
HUMAN-AUDITED here.

**Korolev 2019 cross-check (NO / orthogonal).**  arXiv:1911.09981 estimates
Kloosterman sums over *primes* to *composite* modulus
$q$, non-trivial for $q^{3/4+\varepsilon}\le X\ll q^{3/2}$ (the §2 row
records eq. (7) with the $X\ge q^{7/10+\varepsilon}$ premise — flagged, not
re-audited).  His $(X,q)$ live on the Input-I side; nothing in his saving
structure speaks to a per-class sieve *level* for smooth-indicator
progressions, so there is no bridge by which it could dodge the $Q=y$
depth collapse; it serves only the prime-side (Weyl-sum, Input III) route.
*Premise note for Main:* the §2 row's tail "$X\ge q^{7/10+\varepsilon}$"
differs from the arXiv abstract's $q^{3/4+\varepsilon}$ range — likely the
eq.-specific refined form; re-check at his eq. (7)
before any $\theta$-updating consumes it.

**Verdict: NO — the effective Shiu margin never enters (45) with constant
below 1.**  At Shiu's own depth, $R\ge57.708\,(M/\phi(M))\ge57.7$ across the
covered half-band; at optimal Selberg balance $R\ge5.77\,(M/\phi(M))$; at
idealized depth $R\to2(M/\phi(M))\ge2$; above $M=2\sqrt N$ there is no
theorem.  The binding obstruction is the **per-class Selberg sieve level**
(Lemma 2's $1/\log z'$ per-class shape with $\log z'=\frac{\log(y/M)}{20}$
collapsing across the band) — an $e_{\rm req}$-analog in *sieve-level
units*: closing (45) through this cascade would need $\log z'>M/\phi(M)$,
against $\log(y/M)/20$ available (shortfall factor $\frac{40}{\log2}$ at
the frontier of the covered half).  Constant below 1 via Shiu's cascade is
not a missing-effective-constant state to be repaired but a **shape
limitation**: no choice of depth, no effectivization of any lemma, and no
Korolev-style input closes it.

## 4. Input III — window primes

(46) needs only that some class among $L$ consecutive ones holds a prime of
$[N,2N)$. The relevant sums are Weyl sums $\sum_{p\sim N}e(hp/M)$, $|h|\le M/L$.

| Bound | Hypotheses | Verdict |
|---|---|---|
| Vinogradov/Vaughan: $\sum_{p\le x}e(\alpha p)\ll(x^{4/5}+x^{1/2}q_\alpha^{1/2}+xq_\alpha^{-1/2})\log^cx$ | $|\alpha-a/q_\alpha|\le q_\alpha^{-2}$ | PARTIAL: the $x^{4/5}$ term caps the route at $M\ll LN^{1/5}$ — weaker than Section 15's own theorem |
| GRH via the explicit formula | — | PARTIAL: error $x^{1/2}\log^2x$ caps at $M\ll L\sqrt N/\log^3N$ — the same $\sqrt N$ wall, with the factor $L$ |
| Montgomery, *Michigan Math. J.* 17 (1970): $\psi(x;q,a)=x/\varphi(q)+O(x^{1/2+\varepsilon}q^{-1/2})$ | conjectural, individual $q$ | APPLIES conditionally to every $q\le x^{1-2\varepsilon}$ — the only standard hypothesis that reaches our whole corner. Friedlander–Granville, *Ann. of Math.* 129 (1989) and *Compos. Math.* 81 (1992), bar only its *uniform* version, not any individual modulus |
| Li, arXiv:2505.09629 (2025): prime minorant with level $10/19$ to smooth moduli | $M$ smooth, $M\le N^{10/19-\varepsilon}$ | PARTIAL and genuinely new: covers $\theta\in(\frac12,\frac{10}{19})$ for smooth moduli only |
| Linnik/Xylouris, *Acta Arith.* 150 (2011) 65–91, $L_0\le5.18$ | $x\ge q^{L_0}$ | FAILS at $\theta>\frac12$: we have $x\approx q^{1/\theta}\le q^2$ |
| Siegel–Walfisz | $q\le(\log x)^A$ | FAILS trivially |

Note that GRH is *weaker than Section 15's unconditional theorem* here: at
$q\approx\sqrt N$ the main term $\mathrm{li}(x)/\varphi(q)\approx\sqrt N/\log N$
is below the GRH error $\sqrt N\log^2N$, so the conditional route reaches only
$q\le\sqrt N/\log^{2+}N$ while (38) reaches $\sqrt N/2$.

## 5. Input IV — the level congruence of Section 17 — DISCHARGED

Unlike I–III this input is *not* missing: a classical theorem covers it exactly,
which is how Section 17 proves R(1). It is recorded here because it is the only
place in the project where an external analytic theorem is load-bearing for a
result labelled PROVED, and because the same theorem is listed as PARTIAL for
input III above — the difference is instructive.

**Needed.** For a band $r\in\{1,2\}$, a prime $q$ and a cofactor $t$ with
$q\nmid t$, an asymptotic lower bound for
$$\#\{\,p\ \text{prime},\ p\in[y,2y):\ pt\bmod q^r\in S_r(q)\,\},$$
where $S_r(q)$ is the lower-half level set (50), of size $\ge q^r2^{-r}(1-1/q)$.
Parameters: $y=(2x)^{\alpha}$ with $\alpha\ge\frac12+\eta$, and
$q^{r}\le y^{1/(1+\eta)}$ by construction (49).

| Bound | Hypotheses | Verdict |
|---|---|---|
| Vinogradov (1937), in Vaughan's form: $\sum_{n\le y}\Lambda(n)e(n\theta)\ll(y\,q'^{-1/2}+y^{4/5}+\sqrt{y q'})\log^4y$ | $\lvert\theta-a/q'\rvert\le q'^{-2}$, $(a,q')=1$ | **APPLIES.** Here $\theta=ht/q^r$ is exactly rational with reduced denominator $q'\in[q,q^r]$, so with $q\ge(2x)^{1/4}$ and $q^r\le y^{1/(1+\eta)}$ all three terms are $O(y^{1-\delta})$, $\delta=\delta(\eta)>0$ |
| Mertens' second theorem | — | APPLIES; supplies the $q$-sums and, alone, Lemma 17.2's $\log2$ |
| Prime number theorem | — | APPLIES; supplies the main term of the $p$-count |
| Bombieri–Vinogradov | average over $q\le x^{1/2-\varepsilon}$ | **NOT USED, and would not help.** The level set is $\approx q^r/2$ classes, and summing a per-class maximum over that many classes costs a factor $q^r/2$, which no level-of-distribution theorem returns |

**Why the same bound is PARTIAL for input III and APPLIES here.** In §4 the sum
has length $N$ and modulus $M\approx\sqrt N$, so the main term $N/M$ is
$\approx\sqrt N$ while the error $\sqrt{NM}$ is $\approx N^{3/4}$ — the
instrument is being asked for cancellation on a scale where none is available.
In Section 17 the modulus sits *below* the sum length by a fixed power,
$q^r\le y^{1-\delta}$, which is the regime Vinogradov's method was built for.
The binding quantity is the ratio $\log q'/\log y$, not the absolute size of
$q'$; input III forces that ratio to $\tfrac12$, Section 17 keeps it below $1$
by construction. This is the reusable lesson: an equidistribution demand made
against a *long* sum is cheap even at large modulus, and Section 17 obtains its
long sum by putting the large prime — not the smooth cofactor — in the free
variable.

**Unverified at source:** none. Vinogradov's bound is standard (Vaughan,
*The Hardy–Littlewood Method*, 2nd ed., Thm 3.1); Mertens and PNT likewise.

## 6. Input V — consecutive smooth pairs — DISCHARGED

**Needed.** A positive lower density of consecutive $w,w+1$ with both largest
prime factors below their $\sqrt{2\,\cdot}$ thresholds. This makes both
coordinates free of the short-cofactor class $A$ and proves the strict
Bonferroni inequality (56).

**Cited results.**

| Source | Statement | Verdict |
|---|---|---|
| Yang, arXiv:2607.16032, Theorem 1.5 | positive density with $P^+(n)<P^+(n+1)<x^{41/107+\varepsilon}$ | **APPLIES** directly |
| Pascadi, arXiv:2505.00653v2, Theorem 1.5 | absolute-value smooth Bombieri--Vinogradov through moduli $x^{5/8-\varepsilon}$ | **APPLIES** to Yang's two-prime construction |
| Section 19.1 | substitutes $5/16+5/16=5/8$, leaving exponent $3/8$ | **DISCHARGES V** with $3/8+\varepsilon$ |

The transfer is explicit. Choose primes in the disjoint exponent ranges
$(5/16-3\delta,5/16-2\delta]$ and
$(5/16-2\delta,5/16-\delta]$. Their product is below Pascadi's level, and the
remaining cofactor is below $x^{3/8+5\delta}$. Pascadi controls the summed
progression error; Mertens and $\Psi(x,x^{1/C})\sim\rho(C)x$ give a positive
main term. `data/smooth_pair_transfer.json` checks every exponent inequality.

**Correction.** The prime--prime parameterization $pa+1=p'b$ is still binary:
its modulus $pp'>2w$ exceeds the range, so Selberg, Rosser--Iwaniec and GRH do
not produce that lower bound. Chen switching supplies a $P_2=rs$ partner, but
class $A$ then requires the extra imbalance $r>2bs$ when $w+1=brs$; standard
Chen gives no such inequality. The error was treating failure of this
parameterization as failure of the event. Yang and Pascadi instead count the
smooth complement by Type-I/Type-II dispersion.

## 7. Input VI — structured $q$-adic half-lift dispersion — OPEN

For $q\Vert n+1$, first-level compensation requires

$$
n\bmod q^2\in\mathcal H_q
=\{-1+qr\bmod q^2:\ q/2<r<q\}.
$$

The dense pair theorem permits $q\le x^{3/8+\varepsilon}$, so $q^2$ reaches
$x^{3/4+2\varepsilon}$. Pascadi's modulus-size range is
$x^{5/8-\varepsilon}$; it reaches $q^2$ only through $q\le x^{5/16}$, leaving
a prime strip of width $1/16+\varepsilon$ and a top exponent gap $1/8$.

In the actual construction $n+1=p_1p_2\ell$. The selected primes pass their
first level exactly when

$$
D_{p_1,p_2}(\ell)=
\mathbf1_{\{2(p_2\ell\bmod p_1)>p_1\}}
\mathbf1_{\{2(p_1\ell\bmod p_2)>p_2\}}=1 .
$$

The centered expansion has four different costs:

| frequency class | quotient modulus | generic $n$-lift | verdict |
|---|---:|---:|---|
| zero | $1$ | base $p_1p_2$ | **APPLIES:** Pascadi main term |
| first/second one-sided | $p_1$ or $p_2$ | limiting $15/16$ | **PARTIAL:** Drappeau--Shparlinski saves $1/32$ individually; naive pair sum loses $19/32$ |
| both nonzero | primitive $p_1p_2$ | limiting $5/4$ | **FAILS as a generic modulus:** beyond the $n$-range |

The exact target is (74):
$\mathcal E_{10}+\mathcal E_{01}+\mathcal E_{11}=o(x)$. Pascadi's triple
convolution allows arbitrary scalar outer weights and inner coefficients
independent of the modulus. The mask is neither: it couples the factorization
$p_1p_2$ to $(n+1)/(p_1p_2)$. Parseval reaches order $x$ on the one-sided
terms and $x^{9/8}$ on the bilinear term, where the trivial bound is already
order $x$; it supplies no power saving.

The tail is separate. Repeated selected primes are $O(x^{11/16+o(1)})$, and
repeated factors $q>x^\zeta$ are $o(x)$ (Section 21). Prime factors of
$\ell$, ultra-small factors of $n$, and small repeated primes remain open.
The finite census has $427/2{,}938=0.145337$ full survival; it is evidence,
not proof.

## 8. Input VII — mixed exceptional-spectrum large sieve — OPEN

The exact mask-preserving Poisson identity is Theorem 23.1. In Pascadi's first
dispersion sum, complementary variables $k_1\equiv k_2\pmod r$ produce

$$
\widehat\Phi\!\left(
\frac{rh-ak_1+bk_2}{rH}\right)e_r(-h\overline{k_1}),
$$

with $a,b$ the two mask frequencies. The $h$-interval retains length
$H=x^{1/4+o(1)}$ but its center $(ak_1-bk_2)/r$ depends on the same
complementary variables used later as Kloosterman variables.

| dispersion quantity | exponent |
|---|---:|
| native Poisson length | $1/4$ |
| balanced mask period | $5/16$ |
| spectral level / unmasked support | $1/2$ |
| masked moving-center support | $7/8$ |
| post-Cauchy support inflation | $3/8$ |
| support gap after the native additive-saving benchmark | $1/8$ |

Lowering the outer level to $3/5$ makes the inflation $2/5$ and the benchmark
residual $1/5$, so the older Drappeau range does not evade the moving center.

Thus the earlier $1/16$ comparison is only pre-Poisson. The exact missing
theorem is (94): the paired moving-center extension of Pascadi's
`lem:expo-bound-convo`, plus its one-mask cross variant, with

$$
\mathcal M=\mathcal L=x^{3/8+o(1)},\quad
\mathcal N=H=x^{1/4+o(1)},\quad K=R=x^{5/8+o(1)}.
$$

Additive reciprocity transfers the primitive $r^2$ phase to an inverse-square
trace modulo the smooth complementary $k$; Drappeau--Shparlinski's nonlinear
theorem requires a prime modulus. Even an optimistic generic square-moduli
benchmark gives $x^{17/16+o(1)}$ one-sided and $x^{23/16+o(1)}$ bilinear
bounds, both dominated by the trivial order-$x$ estimate. An unbalanced Chen
$P_2$ either introduces a short-cofactor failure or adds a third $q$-adic mask.

Pascadi explicitly states that no corresponding mixed
additive/multiplicative exceptional-spectrum large sieve is known. A search
through August 2026 found no unconditional theorem covering (94).

## 9. Regime map

With $M=N^\theta$ and $L=\lceil m/2\rceil$:

| Regime | Status | Instrument |
|---|---|---|
| $M\le\sqrt N/2$ | **PROVED**, all $N$, all classes | Theorem 15.3 (Bertrand + one inversion) |
| $\sqrt N/2<M\le\sqrt{2N}$ | proved for non-smooth classes; sharp obstruction for $m$-smooth $M$ | Theorem 15.3 vs Proposition 16.2 |
| $\sqrt{2N}<M\lesssim L\sqrt N/2$ | **CERTIFIED FINITE** per instance | covering criterion, $O(L\log L)$ per prime |
| $\theta\le17/32-\varepsilon$, $9\le m\le C\sqrt N$, fixed $C,\varepsilon>0$ | **CANDIDATE; independent review pending**, all sufficiently large $N$, arbitrary $M$ and classes | F–S 2011 + fixed rectangle + five-position gcd bound + truncated discrepancy; §15.6 |
| $\theta\le\frac12+O(10^{-3}\varepsilon^4)$, prime reduced modulus | primitive-sum comparison only; class/position hypotheses still required | Bourgain–Baker |
| $\theta\in(\frac12,\frac{10}{19})$, smooth $M$ | conditional on a cited 2025 preprint; subsumed by the $17/32$ row on its overlap | Li's prime minorant |
| $\theta<1$, any $M$ | **CERTIFIED FINITE** at $N=10^7$ for $m\ge27$ | survivor census, multi-position |
| $\theta<1$, all $N$ | **OPEN** | (45) with an effective constant, or (46) |

## 10. What to do next, in order

1. Input VII: prove the paired moving-center estimate (94) and its one-mask
   cross variant. Zero frequency is closed. The exact post-Cauchy support is
   $x^{7/8}$ against spectral level $x^{1/2}$; coefficient $\ell^1$ and
   $\ell^2$ masses are already controlled by Theorem 23.1.
2. Once (74) closes, finish Input VI's separate tail: factors of $\ell$,
   ultra-small factors of $n$, and small repeated primes. Selected-prime
   repetition and all repeated factors above $x^\zeta$ are already proved
   negligible.
3. In parallel at the class frontier, effectivize Wright's commented-out
   indicator corollary (LOCALIZATION §3 row; arXiv:2508.17217v2, source
   comment block) at $u=2$: un-comment, complete the $\eta\to1$ endpoint, and
   drive the constant against $1$. It is the only mechanism on record with the
   shape $\Psi\ll y\,\rho(u_y)^{1+o(1)}$ per class that (45) needs.
   *Margin note (2026-08-29, §3.1 audit):* the corollary is REFUTED-for-purpose
   at our corner — the (45) target is $\eta=\frac12$ (with $y=Q=\sqrt{4N}$,
   $u_y=1$, $\rho(1)=1$), not $\eta\to1$; both printed forms are vacuous or
   zero-margin there, the honest mechanism breaks even only at
   $u^\ast\approx2.2125>2$ (certified), and the needed correction is a
   $\rho(2)$-exponent $2.3\%$ beyond the method's own ceiling (details and
   repair directions in §3.1; reproduction artifact
   `e389/wright_indicator_corollary_margin.py`). What survives of this item is
   the route it names: an effective per-class bound with the *shape*
   $\Psi\ll c\,y\,\rho(2)^{e}$, $e>0.8801$ certified-necessary at
   $M=\sqrt{2N}$ — from a $c$-split that preserves (rather than washes out)
   the Dickman gain at $\eta=\frac12$, or from a different instrument.
4. Failing that, fall back to the Weyl-sum route of input III; only one class
   in each window of $L$ must hold a prime.
5. Input I: obtain an independent Opus verdict on the existing restricted
   candidate in campaign `20260908T155952Z_9f8455ad_30f0ba9e47a9`, including
   the floor-count and counting-convention corrections. Do not rerun the
   finite identities or promote the candidate before that verdict.
   Only afterward consider extending the $m\ge9$ argument to smaller
   windows; the $m=2,q=1$ prime-spike obstruction must remain explicit.
