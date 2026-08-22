# Improved certified bounds for the simple-cubic Ising critical coupling

## Result

**[THEOREM from EXTERNAL EXACT COMPUTATION]** For the nearest-neighbour
ferromagnetic Ising model on \(\mathbb Z^3\), with \(K=\beta J\),

\[
\boxed{
0.2122119011661678393310862783954278184914
\;\le K_c\le\;
0.2527310098586630030260020266135701299926 .
}
\]

The lower endpoint improves the incumbent certified value
`0.2074277114992039908436804465100887760027`; the upper endpoint is unchanged.
**[COMPUTATION]** The outward-decimal interval width decreases from
`0.04530329835945901218232158010` to `0.04051910869249516369491574822`;
the lower endpoint moves upward by `0.004784189666963848487405831885`.
Every displayed endpoint is rounded outward. The benchmark `0.221654626` lies
inside this interval, but is used only as a comparison-only sanity check and has
no role in selecting, fitting, or rounding either endpoint.

## 1. Improved lower endpoint from the exact 36-step SAW count

**[EXTERNAL COMPUTATION]** R. D. Schram, G. T. Barkema, and R. H. Bisseling,
*Exact enumeration of self-avoiding walks*, J. Stat. Mech. (2011) P06019,
[doi:10.1088/1742-5468/2011/06/P06019](https://doi.org/10.1088/1742-5468/2011/06/P06019),
Section IV first paragraph and Table I, report the exact simple-cubic count

\[
c_{36}=2\,941\,370\,856\,334\,701\,726\,560\,670.
\]

Their Section II derives the length-doubling inclusion–exclusion identity used
for the enumeration. The complete five-page source is cached as
`sources/fulltext/schram_barkema_bisseling2011.pdf`, SHA-256
`898d8333e8d70acd279e29f226ff108550b132463aa662f86df590b0ad72cb12`.

**[THEOREM]** Write \(n=36q+r\), \(0\le r<36\). Repeated
submultiplicativity gives
\[
c_n=c_{36q+r}\le c_{36}^{q}c_r .
\]
All of the finitely many \(c_r\), \(0\le r<36\), are listed exactly in the same
table. Therefore the high-temperature self-avoiding-path majorant proved in
`proofs/kc_bounds.md` obeys, for \(v\ge0\),
\[
\chi(v)\le\sum_{n\ge0}c_n v^n
\le
\frac{\sum_{r=0}^{35}c_r v^r}{1-c_{36}v^{36}}
\qquad (c_{36}v^{36}<1).
\]
This proves susceptibility finiteness directly for every
\(v<c_{36}^{-1/36}\), with no asymptotic fit and no use of Fekete's lemma.
Consequently
\[
v_c\ge c_{36}^{-1/36},\qquad
K_c\ge\operatorname{atanh}\!\left(c_{36}^{-1/36}\right).
\]
Equivalently, the usual connective-constant formulation follows from Fekete:
\(\mu=\inf_n c_n^{1/n}\le c_{36}^{1/36}\), so
\(v_c\ge1/\mu\ge c_{36}^{-1/36}\). In either form, the input is an **upper**
bound on growth, never the fitted estimate of \(\mu\).

**[COMPUTATION]** The corresponding directed enclosure of the rigorous growth
upper bound is
\[
4.782796893296140857201842295667809952862374008807568163493888101353618066774602108321
\le c_{36}^{1/36}\le
4.782796893296140857201842295667809952862374008807568163493888101353618066774602125188,
\]
improving the incumbent finite-count value \(4.8899016751\ldots\).

**[COMPUTATION]** Exact-integer bisection first brackets
\(c_{36}^{-1/36}\) between adjacent rationals with denominator \(10^{90}\).
`mpmath.iv` at `iv.dps=80` then gives

\[
\begin{aligned}
0.2122119011661678393310862783954278184914094835023361156798789695999003093940121097166
\le{}&\operatorname{atanh}(c_{36}^{-1/36})\\
\le{}&
0.2122119011661678393310862783954278184914094835023361156798789695999003093940121120886.
\end{aligned}
\]

Rounding the lower enclosure downward to 40 decimal places yields the displayed
certified endpoint.

### 1.1 Audit of the stronger published \(4.7114\) claim

**[EXTERNAL, REJECTED AS A CERTIFICATE]** The publisher abstract for
D. MacDonald et al., *Self-avoiding walks on the simple cubic lattice*,
J. Phys. A **33** (2000), 5973–5983,
[doi:10.1088/0305-4470/33/34/303](https://doi.org/10.1088/0305-4470/33/34/303),
ends with the exact sentence:

> “A new, improved rigorous upper bound for the connective constant
> \(\mu<4.7114\) is obtained.”

Only the abstract was obtainable: the full text is paywalled, no repository
copy was found through Crossref, OpenAlex, or Unpaywall, and consequently no
page/section/equation, theorem statement, finite matrix, or numerical
certificate behind `4.7114` could be audited. The sentence is an abstract
numerical conclusion, not a theorem or corollary visible in the material
obtained. Therefore `4.7114` is **not used** in the certified interval. This is
deliberately stricter than trusting the word “rigorous” in an abstract.
If its proof were audited, `4.7114` would imply the stronger comparison value
`0.2155275067875380527127805211792797739268`. We decline that value solely
because its proof is unavailable for audit, not because we have evidence that
the published claim is false.

For contrast, Schram et al. explicitly label their nearby
\(\mu=4.6840401(50)\) as a preliminary direct-fit **estimate** (Section IV);
that estimate is also not used. The present endpoint uses only their exact
integer \(c_{36}\), the submultiplicative theorem, and directed rounding.

## 2. Why the implemented Simon–Lieb route saturates below it

For the free finite set \(B\), the modified Simon–Lieb criterion used in the
frozen experiment is exactly

\[
\kappa_B(v)=v\sum_{x\in B}q_B(x)
 \langle\sigma_o\sigma_x\rangle_{B,v},
\qquad q_B(x)=6-\deg_B(x),
\]

and \(\kappa_B(v)<1\) implies \(K<K_c\).

**[LEMMA] Symmetry grouping does not sharpen this criterion.** The stabilizer
of the rooted `4x4x32` box only permutes equal summands in the positive boundary
sum. Orbit grouping changes storage and arithmetic cost, but neither
\(\kappa_B\) nor its root.

**[LEMMA] The proposed tail replacement cannot sharpen the exact sum.** The
transfer computation already includes every end-cap and axial-tail term
exactly, and all terms are nonnegative. Replacing an exact positive tail by a
proved upper bound can only leave \(\kappa_B\) unchanged or increase it; it
cannot certify a larger root. Deleting or down-weighting a crossing term is not
licensed by the Simon–Lieb theorem.

**[COMPUTATION]** An independent binary64 spin-transfer diagnostic at the
certified rational \(v=0.20234669789378\) decomposes the strongest box as

| contribution | value | fraction of \(\kappa\) |
|---|---:|---:|
| four lateral faces | `0.9999967239960587` | `0.9999967239965325` |
| two end caps | `0.0000032760034679322187` | `0.0000032760034679` |
| total | `0.9999999999995262` | `1` |

Thus `99.99967239965325%` of \(\kappa\) is lateral and only
`0.00032760034679%` comes from the two caps.

The largest rooted symmetry orbit has representative `(0,1,15)`, multiplicity
2, \(q=1\), correlation `0.2533595100378708`, and contribution
`0.10253292047229834`; the ten largest orbits contribute
`0.45749264778602244` of \(\kappa\). Thus the dominant terms are lateral
boundary sites near the origin, not the distant caps. This diagnoses the
`4x4` transverse boundary as the bottleneck and explains why changing length
from 16 to 32 barely moves the root.

**[COMPUTATION, NOT A BOUND]** Fitting the three numerical roots for lengths
8, 16, and 32 to \(v_L=v_\infty+a\rho^L\) gives
\(K_\infty\approx0.20517892018755\). This fit is descriptive only, but it is a
useful failure certificate: the observed axial convergence is far below even
the old SAW endpoint. No endpoint uses this fit.

**[LEMMA]** For the tested free rectangular regions, the alternative
Griffiths/Simon finite-region formulation reduces to the same nonnegative
boundary mass. Calling the same crossing sum a different finite-region
inequality is not a genuinely sharper criterion. A non-rectangular set or an
additional correlation inequality could still be a new route, but no such
proved criterion was found here.

## 3. Upper endpoint and the genuinely different Peierls attempt

### 3.1 Incumbent theorem and its source of slack

**[THEOREM]** The incumbent upper endpoint is produced by the
Fröhlich–Simon–Spencer reflection-positivity/infrared bound, with the physical
bond normalization fixed in `proofs/kc_bounds.md`:

\[
\widehat G(k)\le
\frac{1}{2K(3-\cos k_1-\cos k_2-\cos k_3)}.
\]

The exact spin sum rule then proves long-range order when

\[
K>\frac{I_3}{2},\qquad
I_3=\frac1{(2\pi)^3}\int_{[-\pi,\pi]^3}
\frac{d^3k}{3-\sum_i\cos k_i},
\]

so \(K_c\le I_3/2=W_{\rm sc}/6\). The Watson integral is already evaluated by
an exact Glasser–Zucker gamma identity and directed rounding. Its numerical
quadrature or special-function evaluation is therefore not what controls
sharpness. The slack lies in the Gaussian-domination inequality applied
separately at every nonzero momentum before the exact spin sum rule.

The imported 80-digit enclosure gives the unchanged outward upper endpoint

\[
K_c\le0.2527310098586630030260020266135701299926.
\]

### 3.2 Peierls contours with exact finite surface counts

**[THEOREM]** Claudio Bonati, *The Peierls argument for higher dimensional
Ising models*, [arXiv:1401.7894](https://arxiv.org/abs/1401.7894), equations
(17)–(23), proves the three-dimensional estimate below. The full source is
cached as `sources/fulltext/bonati2014_peierls.pdf`, SHA-256
`2bd9daa39aad312c9a1141b404a9de19020c85e97d0b779644ddfe71567acd36`.
With \(x=9e^{-4K}\),

\[
\frac{\langle N_-\rangle}{N}
\le
F(x):=\frac{3x^3(9-11x+4x^2)}{16(1-x)^3}.
\]

If \(F(x)<1/2\), the plus-boundary thermodynamic state has positive
magnetization. Hence the equality \(F(x)=1/2\) yields a rigorous candidate
upper bound.

**[COMPUTATION]** Exact-integer bisection brackets the unique root by adjacent
rationals with denominator \(10^{60}\):

\[
\frac{450170115456100294887244830222126761823921842208355006402280}{10^{60}}
<x_*<
\frac{450170115456100294887244830222126761823921842208355006402281}{10^{60}}.
\]

Directed interval evaluation of \(-\tfrac14\log(x_*/9)\) gives

\[
K_c\le 0.7488385776610204938868995703971670322475,
\]

which is rigorous but far weaker than the infrared bound.

**[COMPUTATION]** The repository surface enumerator independently gives the
exact plus-boundary broken-bond prefix for the `3x3x4` box:

\[
1+36x^6+75x^{10}+555x^{12}+250x^{14}+2102x^{16}
+5697x^{18}+8412x^{20}+27806x^{22}+58192x^{24}
+122144x^{26}+276078x^{28}+\cdots .
\]

These are exact finite-volume configuration counts, not counts of
infinite-volume connected contours surrounding a fixed site. They therefore
cannot replace the Peierls theorem's infinite tail. A finite exact prefix can
improve an endpoint only together with a proved tail domination in the same
contour normalization; no such domination was established here. Promoting the
prefix alone would be an unproved improvement and is explicitly rejected.

## 4. Classification and reproduction

| claim | classification |
|---|---|
| \(c_{36}=2941370856334701726560670\) | **[EXTERNAL COMPUTATION]**, complete primary source cached with SHA-256 |
| \(\mu\le c_{36}^{1/36}\) and \(K_c\ge\operatorname{atanh}(1/\mu)\) | **[THEOREM]**, submultiplicativity and SAW domination |
| new lower decimal | **[COMPUTATION]**, exact rational bracket and directed rounding at 80 dps |
| abstract claim \(\mu<4.7114\) | **[EXTERNAL, REJECTED]**, proof/full text not obtained; not used |
| Simon–Lieb orbit/cap diagnosis | **[COMPUTATION]**, descriptive binary64 transfer |
| Simon–Lieb grouping/tail conclusions | **[LEMMA]**, positivity and exact-sum argument |
| infrared upper endpoint | **[THEOREM]** plus imported directed rounding |
| Peierls formula | **[THEOREM]**, Bonati equations (17)–(23) |
| Peierls root | **[COMPUTATION]**, exact-rational bracket plus directed rounding |
| finite surface prefix | **[COMPUTATION]**, exact integers; finite scope only |

Reproduce from the repository root:

```text
.venv/bin/python experiments/e33_bounds_improve.py
.venv/bin/python tests/test_bounds_improve.py
```

The experiment writes `results/bounds/improved_bounds.json`, including every
candidate, theorem, certificate, comparison with the incumbent, and explicit
rigor scope.
