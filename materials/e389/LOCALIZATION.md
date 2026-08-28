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
predicts. Three inequivalent inputs would close it. In increasing weakness:

| # | Input | Statement | Section |
|---|---|---|---|
| I | prime localisation | some prime $p\sim P$ has $ap^{-1}\bmod q$ in a prescribed arc | 15.6, hypothesis IP |
| II | survivor bound | $\Psi$-type upper bound $S(N,M,a)<N/M$ for an individual modulus | 16.1, (45) |
| III | window primes | no $L$ consecutive residues $\bmod M$ are prime-free in $[N,2N)$ | 16.3, (46) |
| IV | level congruence | primes in $[y,2y)$ hitting $\ge q^r2^{-r}$ prescribed residues $\bmod q^r$, $q^r\le y^{1-\delta}$ | 17.6, **DISCHARGED** by Vinogradov |

III is weaker than I, and II is independent of both. The census in
`data/survivor_census.json` shows II and III holding with room to spare at
$N=10^7$; neither is a theorem.

## 2. Input I — inverse-prime sums (Kloosterman with primes)

Completion (15.6, eq. 41) turns IP into $\sum_{p\sim P}e_q(h\,a\,p^{-1})$. Our
parameters: $q=M\in[\sqrt{2N},N]$, $P$ free in $[\sqrt{4N+2m},N/d_i]$, minimal
$P\approx2\sqrt N$, so $P\approx q^{1/(2\theta)}$ with $M=N^\theta$ — at
$\theta\approx\tfrac12$ this is $P\approx q$, comfortably inside the published
ranges. A saving $q^{-\delta}$ uniform in $h\not\equiv0$ evicts
$M\le N^{1/(2(1-\delta))-\varepsilon}$ (15.6, eq. 43).

| Bound | Hypotheses | Saving at $P\approx q$ | Verdict |
|---|---|---|---|
| Bourgain, *GAFA* 15 (2005) 1–34, explicit $\delta=0.0005\varepsilon^4$ by Baker, *Acta Arith.* 156 (2012) 351–372 | **$q$ prime**, $(ab,q)=1$, $X\ge q^{1/2+\varepsilon}$ | $q^{-\delta}$, $\delta\sim10^{-3}\varepsilon^4$ | APPLIES for prime $M$ throughout $[\sqrt{2N},N]$; yields $\theta\le\frac12+O(\delta)$ — a genuine but numerically minute crossing |
| Fouvry–Shparlinski, *Acta Arith.* 150 (2011) 285–314 | any $q$, $q^{3/4+\varepsilon}\le X\le q^{4/3-\varepsilon}$ | $X^{15/16}+X^{2/3}q^{1/4}$ | PARTIAL: covers composite $M$ at $\theta\approx\frac12$; $\delta=1/16$ at $X\approx q$ |
| Korolev, eq. (7) of arXiv:1911.09981 | any composite $q$, $X\ge q^{7/10+\varepsilon}$ | $X^{32/37}q^{7/74+\varepsilon}$, i.e. $q^{-3/74}$ at $X\approx q$ | APPLIES for composite $M$; $\delta=3/74=0.0405$ gives $\theta\le0.521$ |
| Korolev, *Res. Number Theory* 6:24 (2020) Thm 1 | any $q$, $(ab,q)=1$, $q^{3/4+\varepsilon}\le X\le(q/2)^{3/2}$, $\Lambda$-weighted, inhomogeneous | $(q^{3/4}/X)^{1/7}$ or $(q^{2/3}/X)^{3/35}$ | APPLIES with an exact printed saving — the cleanest citable composite-modulus input |
| Changa–Korolev, *Math. Notes* 108 (2020) 87–93 | any $q$, $q^{1/2+\varepsilon}\le X\ll q^{3/2}$ | not stated in the abstract | UNVERIFIED — full text inaccessible; would cover the whole corner if the saving is a power of $q$ |
| Kowalski–Michel–Sawin, *Ann. of Math.* 186 (2017) 413–500 | **$q$ prime**, general bilinear forms, $M=N>q^{3/7}$ | power | Engine, not directly our shape; §1.5.2 states composite moduli open |
| Bourgain–Garaev, *Izv. Math.* 78 (2014) 656–707 | **$q$ prime**, multilinear, $\prod|I_j|>q^{1/2+\varepsilon}$ | power | Engine; also records that one published short-sum claim has a proof "in doubt" |

Two technical gaps stand between these bounds and a theorem here, and this
project has not closed either:

1. **The $p$-dependent range.** In eq. (41) the interval $Z_p$ moves with $p$,
   so the $p$-sum is not the bare $\sum_p e_q(h a p^{-1})$ but carries
   coefficients $\hat Z_p(h)$ of unbounded variation. A smooth-weight version of
   the cited bounds, or a decomposition of $(P,2P]$ into ranges of constant
   $\lfloor N/p\rfloor$, is required.
2. **Degenerate $h$.** For $(h,q)=d>1$ the sum reduces to modulus $q/d$ and the
   admissible-length hypothesis must be rechecked at the shorter effective
   length.

**Conclusion for I.** Prime moduli are covered today up to
$\theta=\frac12+O(10^{-3})$; composite moduli up to $\theta\le0.521$ via
Korolev. Reaching every class that meets the block ($\theta\to1$) needs
$\delta\to\frac12$, i.e. square-root cancellation in a sum of length $\approx q$
over primes — beyond current technology.

## 3. Input II — survivors (smooth numbers in one progression)

Target: $\Psi(x,y;q,a)$ with $x=2N$, $y=\sqrt{4N}$, so $u=\log x/\log y\to2$,
and $q\in[\sqrt{2N},x^{2/3}]$, i.e. $q\approx y^{1+2\delta}$. Needed:
$S<N/M$, i.e. an effective constant below $1$ where the Dickman heuristic gives
$\rho(2)=0.30685$.

| Bound | Hypotheses | Verdict |
|---|---|---|
| Balog–Pomerance, *Proc. AMS* 115 (1992) 33–43 (as rendered in Granville's survey, eq. 4.7): $\Psi(x,y;a,q)\ll\Psi_q(x,y)/\varphi(q)$ | $y\ge q^{3/4+\varepsilon}$, $x\ge y^2$, individual $q$, no shape condition | **APPLIES to our exact corner** for every $q\le x^{2/3-\varepsilon}$, right direction, no averaging. FAILS TO CLOSE only because $\ll$ hides an unspecified absolute constant |
| Shiu, *J. reine angew. Math.* 313 (1980) 161–170 (Brun–Titchmarsh for multiplicative $f$) | $x^\kappa\le y\le x$, $k<y^{1-\alpha}$, $(a,k)=1$, individual $k$ | **APPLIES** for every $M\le N^{1-\alpha}$ on the interval $[N,2N)$; FAILS TO CLOSE on the same absolute constant. $u=2$ is its worst corner |
| Wright, arXiv:2508.17217 (2025), Shiu with an explicit $\rho(u)^{C_1}$ factor, $C_1=41/(\alpha\kappa)$ | as Shiu, $f$ smooth-supported | TO CHECK (preprint): right shape, but the constant at $u=2$ was not verified below the threshold |
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
description of II is *one effective constant*, and why making Balog–Pomerance
or Shiu explicit at $u=2$ is the most promising single task in this file.

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

## 6. Regime map

With $M=N^\theta$ and $L=\lceil m/2\rceil$:

| Regime | Status | Instrument |
|---|---|---|
| $M\le\sqrt N/2$ | **PROVED**, all $N$, all classes | Theorem 15.3 (Bertrand + one inversion) |
| $\sqrt N/2<M\le\sqrt{2N}$ | proved for non-smooth classes; sharp obstruction for $m$-smooth $M$ | Theorem 15.3 vs Proposition 16.2 |
| $\sqrt{2N}<M\lesssim L\sqrt N/2$ | **CERTIFIED FINITE** per instance | covering criterion, $O(L\log L)$ per prime |
| $\theta\le0.521$, composite $M$ | conditional on a cited bound, modulo the two gaps of §2 | Korolev inverse-prime sums |
| $\theta\le\frac12+O(10^{-3})$, prime $M$ | same | Bourgain–Baker |
| $\theta\in(\frac12,\frac{10}{19})$, smooth $M$ | conditional on a cited 2025 preprint | Li's prime minorant |
| $\theta<1$, any $M$ | **CERTIFIED FINITE** at $N=10^7$ for $m\ge27$ | survivor census, multi-position |
| $\theta<1$, all $N$ | **OPEN** | (45) with an effective constant, or (46) |

## 7. What to do next, in order

1. Make Balog–Pomerance or Shiu **effective at $u=2$** and compare the constant
   against $1$ (needs a factor $3.26$ of slack, or $1.63$ against index-two
   concentration). This is elementary in method — Friedlander's four-factor
   decomposition, distilled as Lemma 1 of Harman (1999) — and closes input II
   for every individual modulus up to $N^{2/3}$ in one stroke.
2. Failing that, take input III and the Weyl-sum route, exploiting that only one
   class in each window of $L$ must hold a prime.
3. Only then return to input I, whose two technical gaps (§2) must be closed
   before any citation becomes a theorem here.
