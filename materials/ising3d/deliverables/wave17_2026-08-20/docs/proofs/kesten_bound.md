# Audit of the proposed finite Kesten ratio bound

## Result

**[UNRESOLVED]** The proposed finite inequality

\[
\mu^2\le \frac{c_{n+2}}{c_n}\qquad\text{for every }n
\tag{1}
\]

was **not** found in, or derivable from, the sources obtained for this audit.
The auditable Kesten theorem is the ratio **limit**

\[
\lim_{n\to\infty}\frac{c_{n+2}}{c_n}=\mu^2,
\tag{2}
\]

not the finite one-sided statement (1). Therefore no ratio value is promoted to
an upper bound on \(\mu\), and the certified simple-cubic Ising lower endpoint
remains

\[
\boxed{
K_c\ge
\operatorname{atanh}\!\left(
2941370856334701726560670^{-1/36}
\right)
\ge
0.2122119011661678393310862783954278184914 .
}
\tag{3}
\]

**[COMPUTATION, CONDITIONAL ONLY]** If (1) were supplied with an auditable
proof, the exact \(n=34\) data would instead imply

\[
K_c\ge
\operatorname{atanh}\sqrt{\frac{c_{34}}{c_{36}}}
=
\operatorname{atanh}\sqrt{
\frac{132853629626823234210582}
     {2941370856334701726560670}}
>
0.2158152282695194762796979344535594144773,
\tag{4}
\]

which is strictly stronger than both the full incumbent endpoint in (3) and the
task's short display `0.21221190116616784`. Equation (4) is recorded as a
negative-result certificate, not as an Ising theorem.

## 1. Source audit

### 1.1 Kesten's primary paper

**[EXTERNAL]** H. Kesten, *On the Number of Self-Avoiding Walks*, Journal of
Mathematical Physics **4** (1963), 960--969,
[doi:10.1063/1.1704022](https://doi.org/10.1063/1.1704022), has the publisher and
Crossref abstract:

> It is shown that \(\chi_{n+2}/\chi_n-\beta^2\) ... tend[s] to zero as
> \(n\to\infty\), where \(\beta=\lim_{n\to\infty}\chi_n^{1/n}\).

**[UNRESOLVED]** The AIP article and direct PDF routes were blocked by a
Cloudflare/paywall page. OpenAlex reports closed access and no repository full
text. General author/repository searches did not locate a freely available
primary copy. Consequently no page-level assertion beyond the abstract is made
from the primary article, and no unobserved proof is trusted.

### 1.2 Exact Madras--Slade statement and theorem number

**[EXTERNAL]** N. Madras and G. Slade, *The Self-Avoiding Walk*, Chapter 7,
*Pattern theorems*, [doi:10.1007/978-1-4612-4132-4_7](https://doi.org/10.1007/978-1-4612-4132-4_7):

- the official downloadable preview, pp. 229--230, introduces Kesten's result
  as the ratio limits (7.1.1)--(7.1.3);
- the legal Google Books limited preview, book id `xo32YMglDOcC`, p. 248,
  displays **Theorem 7.3.4(a)**:
  \(\lim_{N\to\infty}c_{N+2}/c_N=\mu^2\);
- the theorem's hypothesis is nearest-neighbour self-avoiding walks on
  \(\mathbb Z^d\), in the chapter's standing notation. Part (b) treats fixed
  nonzero endpoints with the parity restriction stated there; parts (c)--(d)
  treat polygons and bridges.

**[EXTERNAL]** The suggested reference “Theorem 1.2.5” is not the ratio theorem.
Google Books search shows `(1.2.5)` as an equation in the book's subadditivity
section on p. 9. The ratio result is Theorem 7.3.4(a).

**[UNRESOLVED]** The complete Springer chapter is subscription-only and was not
cached. The obtained official two-page preview and the page displayed by the
legal limited preview are enough to identify the result and its numbering, but
they do not supply a finite one-sided theorem because none is stated there.

### 1.3 Freely available corroborating full texts

**[EXTERNAL]** R. Bauerschmidt, H. Duminil-Copin, J. Goodman, and G. Slade,
*Lectures on self-avoiding walks*, Clay Mathematics Proceedings **15** (2012),
Section 1.3, printed p. 4, equation (1.17), state precisely

\[
\lim_{n\to\infty}\frac{c_{n+2}}{c_n}=\mu^2.
\]

The authors' complete PDF is cached and hashed below.

**[EXTERNAL]** G. R. Grimmett, *Harry Kesten's work in probability theory*,
[arXiv:2004.03861](https://arxiv.org/abs/2004.03861), Section 5, Theorem 5.1(a),
records the quantitative form

\[
\left|\frac{c_{n+2}}{c_n}-\mu^2\right|
\le A n^{-1/3}
\tag{5}
\]

for some constant \(A\). This is a two-sided error estimate. The theorem gives
neither the sign required by (1) nor a numerical \(A\) that could turn the
\(n=34\) ratio into a useful explicit bound.

### 1.4 Logical consequence of the audit

**[LEMMA]** A limit does not determine the side from which a finite term
approaches it. In particular, (2) alone is consistent with
\(c_{n+2}/c_n-\mu^2\) having either sign at \(n=34\). Likewise, (5) yields only

\[
\mu^2\le \frac{c_{n+2}}{c_n}+A n^{-1/3},
\]

with an unavailable numerical constant. Neither statement proves (1).

**[UNRESOLVED]** Upgrading (4) requires either (i) a page-level theorem and
reproducible proof of the finite one-sided sign, with hypotheses covering
nearest-neighbour SAWs on \(\mathbb Z^3\), or (ii) a complete independent proof.
A citation to the ratio limit is not sufficient.

### 1.5 Cached-source hashes

**[COMPUTATION]** SHA-256 hashes are recomputed by both the experiment and its
standalone test:

| manifest key | obtained scope | local file | SHA-256 |
|---|---|---|---|
| `schram2011_exact_saw` | full | `sources/fulltext/schram_barkema_bisseling2011.pdf` | `898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12` |
| `madras_slade1996_ch7_preview` | official two-page preview | `sources/fulltext/madras_slade1996_ch7_preview.pdf` | `4f3db221eb1a4fb8162b4f1fed325c3f41ad9d7912550deb51ee8a0200bc8135` |
| `bauerschmidt2012_saw_lectures` | full | `sources/fulltext/bauerschmidt_duminil_copin_goodman_slade2012_saw_lectures.pdf` | `f41468f516eea8d4ba5866d8c444636fa579ca57bcb27210204718515a2530cd` |
| `grimmett2020_kesten_work` | full | `sources/fulltext/grimmett2020_kesten_work.pdf` | `1ab2dc69831fcdf3ca26a05506cb6766dffd6d24cbcb8fc427e9555b0f4729b7` |
| `schram2017_bcc_fcc` | full | `sources/fulltext/schram_barkema_bisseling_clisby2017_bcc_fcc.pdf` | `002b53281c491a0589b43139918db1cd22d0c35edfa4f9f9060e744175ebfc2f` |

**[EXTERNAL]** The Kesten primary entry is `abstract_only`, with `local_path` and
`sha256` deliberately null. No unavailable primary text is represented as
cached.

## 2. Exact SAW data

**[EXTERNAL COMPUTATION]** R. D. Schram, G. T. Barkema, and R. H. Bisseling,
*Exact enumeration of self-avoiding walks*, J. Stat. Mech. (2011) P06019,
[doi:10.1088/1742-5468/2011/06/P06019](https://doi.org/10.1088/1742-5468/2011/06/P06019),
Table I, gives the rooted nearest-neighbour simple-cubic counts used here:

| \(n\) | exact \(c_n\) |
|---:|---:|
| 30 | 270569905525454674614 |
| 31 | 1274191064726416905966 |
| 32 | 5997359460809616886494 |
| 33 | 28233744272563685150118 |
| 34 | 132853629626823234210582 |
| 35 | 625248129452557974777990 |
| 36 | 2941370856334701726560670 |

The complete source already used by the incumbent bound is cached under the
first hash in Section 1.5. No fitted value of \(\mu\) is used.

**[COMPUTATION]** Exact cross multiplication gives

\[
\frac{c_{36}}{c_{34}}<\frac{c_{34}}{c_{32}},
\qquad
\frac{c_{35}}{c_{33}}<\frac{c_{33}}{c_{31}}.
\tag{6}
\]

These decreasing same-parity ratios are sanity observations only. No
all-\(n\) monotonicity is assumed.

## 3. Conditional arithmetic certificate

### 3.1 Algebraic chain

**[LEMMA, CONDITIONAL]** The high-temperature SAW domination proved in
`proofs/kc_bounds.md` gives \(v_c\ge1/\mu\). If (1) were available at \(n=34\),
then

\[
\mu\le\sqrt{\frac{c_{36}}{c_{34}}},
\qquad
v_c\ge\sqrt{\frac{c_{34}}{c_{36}}},
\]

and monotonicity of `atanh` would give (4). This implication is valid; only its
premise (1) is unresolved.

### 3.2 Exact adjacent rational brackets

**[COMPUTATION]** Set \(S=10^{90}\) and

\[
\begin{aligned}
M_-&=4705309340807700171373809606438049928499981287517794153022627972588591592036873355061175466,\\
M_+&=M_-+1,\\
V_-&=212525878230217019425992579085858173341341762357141355950356551076118513866301486370354404,\\
V_+&=V_-+1.
\end{aligned}
\]

The experiment verifies with exact integers

\[
\begin{aligned}
c_{34}M_-^2&\le c_{36}S^2<c_{34}M_+^2,\\
c_{36}V_-^2&\le c_{34}S^2<c_{36}V_+^2.
\end{aligned}
\tag{7}
\]

Thus

\[
\frac{M_-}{S}\le\sqrt{\frac{c_{36}}{c_{34}}}<\frac{M_+}{S},
\qquad
\frac{V_-}{S}\le\sqrt{\frac{c_{34}}{c_{36}}}<\frac{V_+}{S}.
\tag{8}
\]

**[COMPUTATION]** `mpmath.iv` directed rounding at `iv.dps=90`, applied
monotonically to the second bracket in (8), gives

\[
\begin{aligned}
0.21581522826951947627969793445355941447733069891440514825753509862465429330062968004479858423464
\le{}&\operatorname{atanh}\sqrt{\frac{c_{34}}{c_{36}}}\\
<{}&
0.21581522826951947627969793445355941447733069891440514825753509862465429330062968004479858552328.
\end{aligned}
\tag{9}
\]

Rounding (9) outward to 40 decimal places gives the conditional interval

\[
[0.2158152282695194762796979344535594144773,
  0.2158152282695194762796979344535594144774].
\tag{10}
\]

### 3.3 Both parities

**[COMPUTATION]** The exact conditional values for the last two ratios of each
parity are:

| \(n\) | parity | conditional \(\sqrt{c_{n+2}/c_n}\) | conditional outward \(K\) interval (40 places) |
|---:|:---:|---:|---:|
| 31 | odd | 4.7072467119533859654361621981160196218793... | [0.2157236226903125885980814718938995436536, 0.2157236226903125885980814718938995436537] |
| 32 | even | 4.7065932992108327990202003454274486526995... | [0.2157545094680394956202703487289228194543, 0.2157545094680394956202703487289228194544] |
| 33 | odd | 4.7058920678834378421923174829346842356235... | [0.2157876666201136904202289494755504846036, 0.2157876666201136904202289494755504846037] |
| 34 | even | 4.7053093408077001713738096064380499284999... | [0.2158152282695194762796979344535594144773, 0.2158152282695194762796979344535594144774] |

Exact comparison of the rational squared arguments, not floating-point
selection, identifies \(n=34\) as the strongest of these four conditional
candidates.

## 4. Certified endpoint after rejection of the premise

**[THEOREM from EXTERNAL EXACT COMPUTATION]** Without (1), ordinary
submultiplicativity still gives \(\mu\le c_{36}^{1/36}\), and the established SAW
correlation majorant gives the incumbent (3).

**[COMPUTATION]** With the same \(S=10^{90}\), exact integer bisection brackets
\(c_{36}^{-1/36}\) between the adjacent rationals with numerators

```text
209082681600312329099111432193182679138662355052911528634811788862371711498508872381089533
209082681600312329099111432193182679138662355052911528634811788862371711498508872381089534
```

and denominator \(S\). Directed evaluation gives

\[
\begin{aligned}
0.21221190116616783933108627839542781849140948350233611567987896959990030939401211091157036711736
\le{}&\operatorname{atanh}(c_{36}^{-1/36})\\
<{}&
0.21221190116616783933108627839542781849140948350233611567987896959990030939401211091157036843668.
\end{aligned}
\tag{11}
\]

The lower endpoint rounded downward to 40 places is exactly the incumbent in
(3). The conditional interval (9) is disjoint and lies above (11), but it is not
certified because its premise was not audited.

## 5. Search for longer exact simple-cubic counts

**[EXTERNAL]** Schram, Barkema, Bisseling, and Clisby,
[arXiv:1703.09340](https://arxiv.org/abs/1703.09340), Introduction, printed p. 2,
say that the exact simple-cubic series through \(N_{\max}=36\) was still the
record as of 2017. Their new exact tables concern BCC and FCC lattices; later
entries labelled as extrapolations are not exact simple-cubic counts.

**[UNRESOLVED]** No freely auditable post-2011 exact simple-cubic \(c_n\) with
\(n>36\) was located in this audit. This is a documented search outcome, not a
theorem that no such count exists.

## 6. Reproduction

**[COMPUTATION]** From the repository root:

```text
.venv/bin/python experiments/e55_kesten_bound.py
.venv/bin/python tests/test_kesten_bound.py
```

Both programs print `PASS` as their final line. The experiment writes
`results/bounds/kesten_bound.json`; the standalone test independently uses
`math.isqrt` for all square-ratio brackets, recomputes the inverse 36th-root
bracket, evaluates directed intervals at 100 dps, verifies every cached-source
SHA-256, and checks that the conditional result is not promoted.
