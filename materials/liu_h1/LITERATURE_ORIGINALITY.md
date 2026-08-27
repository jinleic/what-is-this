# Literature and originality audit: Liu Hypothesis 1

**Full-audit cutoff:** 2026-08-25 (local date)
**Target:** Jingbo Liu, arXiv:2306.08824v1, Section V-A, and the stronger unprojected residual-kernel theorem used locally  
**Search record:** [`literature/search/query_log.json`](literature/search/query_log.json)  
**Latest refresh:** the arXiv `all:"union-closed"` feed was queried again at
2026-08-26T15:10:15Z.  It still reported 103 records, with no record posted
after the full-audit cutoff and no Liu-H1 collision; see
[`arxiv_refresh_2026-08-26.json`](literature/search/arxiv_refresh_2026-08-26.json).

**Source index:** [`literature/metadata/sources.json`](literature/metadata/sources.json)  
**Citation screening:** [`literature/search/citation_screening.json`](literature/search/citation_screening.json)

## Status vocabulary

- **MACHINE-VERIFIED:** a finite file, checksum, or symbolic/computational assertion was checked by a machine in the local project.
- **HUMAN-AUDITED:** a theorem statement, hypothesis, formula, or citation context was read and compared by a human-level audit.
- **COMPUTATIONAL-EVIDENCE:** numerical evidence only; not a proof of an infinite-dimensional sign claim.
- **CITED-DEPENDENCY:** a claim is taken from the cited primary source and was not independently reproved in this literature audit.
- **OPEN:** unresolved, not exhaustively knowable from the searched corpus, or an absence-of-found-prior-art statement.
- **FAILED:** a named source or check could not be accessed or completed.

## Executive conclusion

- **HUMAN-AUDITED — no exact collision found.** None of the accessible sources screened through the cutoff proves the exact endpoint statement
  \[
  (x,y)\longmapsto h_2\!\left(xy+x(1-x)y(1-y)\right)
  \]
  to be conditionally negative semidefinite on signed measures annihilating \(1,x,x(1-x)\), and none states the stronger ordinary negative-semidefiniteness of the residual kernel \(R\) below.

- **OPEN — this is not proof of novelty.** The search result supports cautious wording such as “we are not aware of an earlier proof”; it does **not** justify an unqualified “first proof” or establish legal or historical priority. Two indexed citing items were not fully accessible, and unpublished, private, or unindexed work cannot be excluded.

- **HUMAN-AUDITED — the closest theorem-level prior art is substantial but not subsuming.** Alweiss–Huang–Sellke prove the unperturbed product-entropy case \(h_2(xy)\); Liu proves a perturbative version for \(xy+f(x)f(y)\) when the scale in \(f\) is sufficiently small; Liu then leaves the scale-one case \(f(x)=x(1-x)\) to numerical checks. These are the correct immediate predecessors to cite.

- **HUMAN-AUDITED — several proof ingredients are classical.** Nonnegative power-series/Schur closure and the positive-definiteness of \(1/(1-\langle u,v\rangle)\) on a Hilbert ball predate the local argument. Any originality claim should therefore be confined to the problem-specific entropy regrouping, the explicit Lorentz factorization and diagonal-cone check for \(D_\theta\), and the exact scale-one endpoint theorem—not to those general kernel facts.

- **HUMAN-AUDITED — consequence remains scoped.** Proving Section V-A removes only Liu's positive-semidefiniteness hypothesis. Liu's Section V-B global-minimizer hypothesis remains separate, so the explicit value \(0.382709087918741\) remains conditional.

## 1. Exact target and normalization

### 1.1 Liu's variables

**CITED-DEPENDENCY.** Let \(h_2(p)=-p\log_2p-(1-p)\log_2(1-p)\). In Liu's Section V-A, after writing \(x=1-s\) and \(y=1-t\), the scale-one protocol has
\[
I(x,y)=xy+x(1-x)y(1-y)=xy\bigl(1+(1-x)(1-y)\bigr).
\]
The relevant quadratic form is
\[
Q_H(\nu)=\iint_{[0,1]^2}h_2(I(x,y))\,d\nu(x)\,d\nu(y).
\]
The primary sources are [arXiv:2306.08824v1, Section V-A](https://arxiv.org/html/2306.08824#S5.SS1) and Liu's archived code [`frankl3.m`](literature/code/frankl3.m).

**HUMAN-AUDITED.** The signed tangent space is
\[
\mathcal N_3=
\left\{\nu:\int d\nu=\int x\,d\nu=\int x(1-x)\,d\nu=0\right\}.
\]
Because \(x(1-x)=x-x^2\), this is equivalently the annihilator of \(\{1,x,x^2\}\). In the original \(s=1-x\) variable it is equivalently the annihilator of \(\{1,s,s(1-s)\}\) or \(\{1,s,s^2\}\).

**HUMAN-AUDITED.** The codimension terminology needs care:

- Liu optimizes over **probability measures** with two additional affine moment constraints, so the paper describes a codimension-two constrained probability set.
- Its ambient linear tangent space of signed measures also includes the mass-zero condition and is therefore codimension three.
- Liu's MATLAB projector has three columns—\(1\), the grid coordinate, and \(s(1-s)\)—and consequently removes three linear directions.

**CITED-DEPENDENCY.** Liu's required endpoint assertion is
\[
Q_H(\nu)\le 0\qquad(\nu\in\mathcal N_3),
\]
or equivalently that \(-h_2(I(x,y))\) is positive semidefinite on that tangent space.

### 1.2 Stronger local residual theorem being compared

**HUMAN-AUDITED.** Put
\[
u=1-s,\qquad v=1-t,\qquad A=uv,\qquad B=st,\qquad z=A(1+B),
\]
\[
G(q)=q+(1-q)\log(1-q),
\]
and
\[
R(s,t)=A\bigl[B-(1+B)\log(1+B)\bigr]-G(z).
\]
On \(\mathcal N_3\), the projectively vanishing rank-one terms give
\[
(\log 2)Q_H(\nu)=\iint R(s,t)\,d\nu(s)\,d\nu(t).
\]
The local theorem is stronger:
\[
\sum_{i,j}c_i c_jR(s_i,s_j)\le0
\]
for every finite set of points \(s_i\in[0,1]\) and arbitrary real coefficients \(c_i\), with no moment constraint. By continuity, the corresponding finite-signed-measure formulation follows in the local proof audit. The canonical local implementation is [`../uc/liu_R_nsd.py`](../uc/liu_R_nsd.py); this report does not re-audit that proof.

**HUMAN-AUDITED.** The local proof method is an exact Taylor regrouping
\[
-R=G(A)-AB\log(1-A)
+B^2\int_0^1(1-\theta)
\frac{A}{(1+\theta B)\,[1-A(1+\theta B)]}\,d\theta.
\]
The first two terms are handled by nonnegative power series. For the last term, the denominator kernel
\[
D_\theta(s,t)=(1+\theta st)\,[1-(1-s)(1-t)(1+\theta st)]
\]
is shown to have a degree-three Lorentz factorization
\[
D_\theta(s,t)=a_\theta(s)a_\theta(t)-\langle b_\theta(s),b_\theta(t)\rangle,
\qquad \|b_\theta(s)/a_\theta(s)\|<1,
\]
so \(1/D_\theta\) is a positive-definite Hilbert-ball kernel by its geometric tensor-power expansion.

## 2. Publication, version, proceedings, and erratum trail

| Status | Record | What was established | Limitation |
|---|---|---|---|
| **CITED-DEPENDENCY** | [arXiv:2306.08824](https://arxiv.org/abs/2306.08824) | The arXiv API returns exactly one entry, `2306.08824v1`, published and last revised 2023-06-15. The saved record is [`arxiv_2306.08824.atom.xml`](literature/metadata/arxiv_2306.08824.atom.xml). | An arXiv version count does not rule out unpublished revisions. |
| **CITED-DEPENDENCY** | [CISS DOI 10.1109/CISS59072.2024.10480167](https://doi.org/10.1109/CISS59072.2024.10480167) | IEEE and Crossref identify a six-page proceedings paper, pages 1–6, published 2024-03-13. Its abstract says that additional evaluation results are in arXiv:2306.08824. | The IEEE PDF was subscription-denied; no claim of verbatim identity with the 16-page arXiv PDF is made. See [`ieee_xplore_observation.json`](literature/metadata/ieee_xplore_observation.json). |
| **HUMAN-AUDITED** | [Illinois Experts publication record](https://experts.illinois.edu/en/publications/improving-the-lower-bound-for-the-union-closed-sets-conjecture-vi/) | The institutional record lists the CISS contribution and points evaluation details to the arXiv item; no separate journal article is listed for this title. | Institutional profiles can be incomplete. |
| **OPEN** | Crossref relation data, IEEE record, arXiv history, title/DOI web searches | No correction, erratum, corrigendum, retraction, later arXiv revision, or later journal version was found by the cutoff. Crossref's `relation` object is empty. | This is an absence-of-found-record statement, not proof that none exists. |
| **FAILED** | IEEE proceedings PDF | Direct IEEE PDF endpoints returned access denial in this audit session. | The open arXiv v1 PDF is archived at [`liu_2023_arxiv_2306.08824v1.pdf`](literature/pdfs/liu_2023_arxiv_2306.08824v1.pdf). |

## 3. What Liu actually established before the scale-one hypothesis

### 3.1 The unperturbed product kernel

**CITED-DEPENDENCY.** Alweiss, Huang, and Sellke, [“Improved Lower Bound for Frankl's Union-Closed Sets Conjecture”](https://doi.org/10.37236/12232), Lemma 5, prove concavity of
\[
\mu\longmapsto \iint H(xy)\,d\mu(x)d\mu(y)-\int H(x)\,d\mu(x)
\]
on probability measures with fixed mean. Its second variation says
\[
\iint H(xy)\,d\nu(x)d\nu(y)\le0
\]
for signed \(\nu\) annihilating \(1\) and \(x\). This is the \(f=0\) or scale-zero endpoint of Liu's family, not the scale-one kernel.

**HUMAN-AUDITED.** Their proof uses two integrations by parts and the nonnegative power-series expansions of \(1/(1-xy)\) and \(\log(1/(1-xy))\). It neither has Liu's additional \(x(1-x)y(1-y)\) term nor the third moment constraint.

### 3.2 Liu's small-scale perturbation theorem

**CITED-DEPENDENCY.** Liu's Lemma 11 considers
\[
I_f(x,y)=xy+f(x)f(y),\qquad f(x)=\ell x p(x),\qquad p(1)=0,
\]
and proves the required concavity for **sufficiently small** \(\ell>0\), on the tangent space annihilating \(1,x,f(x)\). The method perturbs the Alweiss–Huang–Sellke integration-by-parts argument and bounds adverse series terms by favorable ones.

**HUMAN-AUDITED.** For Hypothesis 1, \(p(x)=1-x\) and \(\ell=1\). Liu gives no threshold shown to reach \(1\); Section V-A instead switches to grid eigenvalues and finite coefficient-matrix eigenvalues. Thus Lemma 11 does not subsume the endpoint.

### 3.3 Liu's numerical endpoint evidence

**COMPUTATIONAL-EVIDENCE.** [`frankl3.m`](literature/code/frankl3.m) samples spacing \(0.0004\), builds the scale-one entropy Gram matrix, applies a rank-three projector, symmetrizes, and reports a residual eigenvalue at approximately \(10^{-14}\).

**COMPUTATIONAL-EVIDENCE.** [`frankl7.m`](literature/code/frankl7.m) forms a finite coefficient matrix and reports minimum eigenvalues of truncations, with the paper reporting checks through finite cutoff \(L\).

**HUMAN-AUDITED.** Neither computation controls all finite point sets or the infinite coefficient tail. They establish Liu's priority for the conjecture and its numerical evidence, not an exact proof.

## 4. Theorem-level comparison with potentially subsuming literature

| Status | Source/result | Domain and function space | Kernel/formula | Constraints/projector | Proof method | Exact relation to Liu H1 |
|---|---|---|---|---|---|---|
| **CITED-DEPENDENCY** | Alweiss–Huang–Sellke 2024, Lemma 5 | Probability measures on \([0,1]\); signed second variations | \(H(xy)\) | \(\int1\,d\nu=\int x\,d\nu=0\), linear codimension 2 | Integration by parts; nonnegative series | Closest exact sign theorem, but only the unperturbed \(f=0\) kernel. Does not subsume scale one. |
| **CITED-DEPENDENCY** | Liu 2023, Lemma 11 | Probability measures with fixed first and \(f\)-moments | \(H(xy+f(x)f(y))\), \(f=\ell xp(x)\) | Tangent annihilates \(1,x,f\), linear codimension 3 | Perturbative integration by parts and series domination | Same structural family, but only sufficiently small \(\ell\); no result that \(\ell=1\) is included. |
| **COMPUTATIONAL-EVIDENCE** | Liu 2023, Section V-A | Uniform grid and finite Taylor matrices | Exact H1 kernel at \(f=x(1-x)\) | Rank-three grid projector; finite coefficient truncation | Floating-point eigenvalues | Exact conjecture, but numerical only. Establishes priority for the statement, not its proof. |
| **CITED-DEPENDENCY** | Fuglede–Topsøe 2004, Theorems 1–2 | Nonnegative scalars and probability distributions | Jensen–Shannon difference kernel, generated pointwise by \(a\log\frac{2a}{a+b}+b\log\frac{2b}{a+b}\) | Ordinary CND zero-sum coefficients; no polynomial moment projector | Negative-type integral representation/Hilbert embedding | “Entropy” and “negative definite” are genuine, but the formula is a divergence of two inputs, not entropy of Liu's polynomial coupling probability. |
| **CITED-DEPENDENCY** | Cuturi–Fukumizu–Vert 2005, Proposition 1 | Positive measures/densities under addition | Differential entropy as a semigroup-negative-definite function; exponentiated Jensen divergence PD | Ordinary zero-sum condition on measures in the semigroup | Pointwise scalar ND plus Schoenberg exponentiation | Different objects and composition law; does not imply the Liu kernel sign. |
| **CITED-DEPENDENCY** | Hein–Bousquet 2005 | Probability measures with densities | Integrals of homogeneous Hilbertian pointwise metrics/divergences | Ordinary Hilbertian negative type | Pointwise embedding and integration | Related entropy/divergence kernels, but no product-polynomial \(I(x,y)\), no order-three projector. |
| **CITED-DEPENDENCY** | Schoenberg 1938; Micchelli 1986 | General metric kernels; radial/interpolation CPD theory | Abstract negative type and CPD-of-order-\(m\) kernels | Zero-sum or polynomial annihilators of prescribed order | Hilbert embedding; distance-matrix/interpolation theory | Supplies vocabulary and closure principles. No theorem located whose hypotheses recognize Liu's nonradial kernel. |
| **CITED-DEPENDENCY** | Schur 1911 / Pólya–Szegő 1925, as surveyed by Guillot et al. 2026 | PSD matrices or kernels of every finite size | Entrywise products and convergent power series with nonnegative coefficients | No projection | Schur products and closure of the PSD cone | Covers \(G(A)\), \(-AB\log(1-A)\), and geometric-series pieces after they are isolated; it does not make their original signed entropy combination sign-definite. |
| **CITED-DEPENDENCY** | Agler–McCarthy 2000, Theorem 3.1 and universal kernel; McCullough–Quiggin | Hilbert-ball points / complete Pick kernels | \(k(u,v)=1/(1-\langle u,v\rangle)\) up to rank-one scaling | No moment projection | Complete Pick representation; geometric tensor powers | Covers the reciprocal step once the problem-specific factorization \(D_\theta=aa-\langle b,b\rangle\) and \(\|b/a\|<1\) are established. It does not supply those facts or the entropy identity. |
| **CITED-DEPENDENCY** | Boppana 2023; Wakhare 2025 | One scalar variable \(x\in[0,1]\) | \(\alpha_kH(x^k)\ge x^{k-1}H(x)\) | None | Derivatives, Rolle/Descartes, real-root counts | Important binary-entropy inequalities in the same union-closed program, but not a two-variable Gram-kernel theorem. |
| **HUMAN-AUDITED** | Local residual theorem | Finite Gram sets, equivalently continuous finite-signed-measure form on \([0,1]\) after the proof's limit step | \(R(s,t)\) above | None; ordinary NSD | Exact Taylor remainder plus explicit Lorentz/Gram factorization | Stronger than H1 and not stated in any screened source. Originality remains **OPEN**, not proved by this table. |

## 5. Citation and survey trail through the cutoff

**CITED-DEPENDENCY.** Bibliographic indexes disagree in scope: IEEE, Crossref, and OpenAlex expose two DOI-linked citations, OpenCitations exposes three rows (one without a citing DOI), while Google Scholar's merged CISS cluster displayed ten citing records. This is an indexing difference, not a mathematical discrepancy. Raw snapshots are in [`literature/metadata/`](literature/metadata/).

**HUMAN-AUDITED.** The complete screened ledger is [`citation_screening.json`](literature/search/citation_screening.json). The accessible citing works fall into four groups:

1. **HUMAN-AUDITED — background/future-direction citations.** Das–Wu (arXiv:2412.03862v3) cite Liu in the history and suggest adapting improved couplings to kth frequencies; Hachimori–Kashiwabara's combinatorial papers cite the constant-bound literature in their introductions; none addresses H1.
2. **HUMAN-AUDITED — scalar entropy work.** Wakhare's 2025 journal article studies univariate binary-entropy inequalities and root counts; it neither states nor implies the bivariate kernel theorem.
3. **HUMAN-AUDITED — secondary overstatements.** Lu–Raz (arXiv:2405.10639v2) call \(0.38271\) “proven by Liu,” and Colbert's 2025/2026 paper says Liu improved the bound to that value without the primary source's conditional caveat. Neither provides a proof of either numerical hypothesis. Primary-source wording controls.
4. **HUMAN-AUDITED — current rigorous-context signal.** Das–Wu describe the contemporary rigorous numerical frontier as approximately \(0.3823455\) and treat Liu's further coupling improvement as a direction to deploy, which is consistent with the fact that the explicit \(0.382709\) value remains conditional.

**OPEN.** Two Google Scholar citing records could not be fully screened: Wouter Antvelink's 2025 TU Delft undergraduate thesis returned a campus-only/closed response, and Jonas Gebendorfer's 2026 ResearchGate paper returned access controls. Their indexed titles/abstract snippets do not announce an entropy-kernel proof, but their inaccessibility prevents a definitive content claim.

**HUMAN-AUDITED.** A recent 2026 article by van der Hout–Roos was also checked. It discusses general Frankl results and conjectures, mentions improvements only generically as slightly above 38%, does not cite Liu, and contains no negative-semidefinite entropy-kernel theorem.

## 6. Originality classification by component

| Status | Component | Conservative classification | Required attribution |
|---|---|---|---|
| **OPEN** | Exact scale-one conditional NSD of \(h_2(xy+x(1-x)y(1-y))\) | No earlier proof found; plausible theorem-level novelty, but not certified novelty. | State that this is Liu's Section V-A hypothesis and cite both arXiv v1 and CISS DOI. |
| **OPEN** | Stronger ordinary NSD theorem \(R\preceq0\) without projection | No exact antecedent found; potentially the strongest new theorem. | Explain explicitly why projectively vanishing terms make it imply H1. Avoid claiming exhaustive priority. |
| **HUMAN-AUDITED** | Statement, grid test, and finite coefficient tests at scale one | Not novel locally; Liu has priority. | Cite Section V-A, `frankl3.m`, and `frankl7.m`; label these numerical. |
| **HUMAN-AUDITED** | Unperturbed product-entropy kernel | Prior art. | Cite Alweiss–Huang–Sellke, Lemma 5. |
| **HUMAN-AUDITED** | Sufficiently-small perturbation family | Prior art. | Cite Liu, Lemma 11, and distinguish “sufficiently small” from \(\ell=1\). |
| **HUMAN-AUDITED** | Nonnegative power-series and Schur-product PSD steps | Classical, not a novelty claim. | Cite Schur/Pólya–Szegő or a standard source; the archived 2026 survey gives precise modern attribution. |
| **HUMAN-AUDITED** | \(1/(1-\langle c(s),c(t)\rangle)\) Gram-series kernel | Classical Hilbert-ball/complete-Pick construction. | Cite Agler–McCarthy and/or identify the direct tensor-power expansion. |
| **OPEN** | Explicit \(D_\theta\) Lorentz factorization and diagonal time-cone verification | No matching factorization found; plausible problem-specific novelty. | Do not describe the general complete-Pick kernel as new. |
| **OPEN** | Exact Taylor regrouping of \(-R\) into three PSD pieces | No matching identity found; plausible problem-specific novelty. | Present the identity exactly and separate it from standard closure facts. |

## 7. Priority-collision decision

**HUMAN-AUDITED.** A priority collision was defined narrowly: an earlier source had to state the same scale-one kernel with the same three annihilating directions and prove its conditional sign, or prove the unprojected residual theorem from which it follows immediately.

**OPEN — no collision located.** No screened source meets that test. The nearest results stop at either:

- the product kernel \(H(xy)\);
- sufficiently small perturbation scale;
- finite numerical grids/truncations;
- a different entropy divergence or semigroup;
- a scalar entropy inequality; or
- a general kernel-closure theorem that handles only one isolated substep.

**OPEN.** This decision must be revisited if the inaccessible thesis/paper becomes available, if Liu or another author posts a later revision, or if a subject expert identifies an unindexed theorem with the exact formula.

## 8. Manuscript-ready cautious wording

### Recommended originality paragraph

> **HUMAN-AUDITED wording.** Liu introduced the scale-one kernel and supported its conditional negative semidefiniteness numerically in Section V-A of arXiv:2306.08824v1 (published in CISS 2024). The closest rigorous predecessors are the unperturbed product-entropy kernel of Alweiss, Huang, and Sellke and Liu's sufficiently-small-perturbation lemma. We are not aware of a previous proof at the scale-one endpoint. Here we prove a stronger statement: after removing terms annihilated by Liu's three tangent constraints, the resulting residual kernel is negative semidefinite without any projection.

### Recommended method-attribution paragraph

> **HUMAN-AUDITED wording.** The proof uses classical closure of positive-definite kernels under nonnegative sums, limits, and Schur products. Its reciprocal-kernel step is an explicit instance of the Hilbert-ball kernel \((1-\langle u,v\rangle)^{-1}\), familiar from complete Nevanlinna–Pick theory. The problem-specific content is the exact entropy Taylor decomposition and the Lorentz factorization that places the normalized polynomial feature map strictly inside the unit ball.

### Recommended consequence paragraph

> **HUMAN-AUDITED wording.** This settles only Liu's positive-semidefiniteness hypothesis in Section V-A. The separate nine-parameter global-minimizer hypothesis in Section V-B is not addressed; consequently, Liu's explicit value \(0.382709087918741\) remains conditional.

### Wording to avoid

- **HUMAN-AUDITED:** Avoid “the first proof” unless a broader expert priority check is completed; use “we are not aware of an earlier proof.”
- **HUMAN-AUDITED:** Avoid “Liu proved \(0.38271\)”; Liu's primary source says that value is under numerically verified hypotheses.
- **HUMAN-AUDITED:** Avoid calling Schur closure, nonnegative Taylor coefficients, or the Hilbert-ball reciprocal kernel new.
- **HUMAN-AUDITED:** Avoid calling Liu's projector merely codimension two without explaining the probability-mass convention; the signed tangent space and numerical projector have codimension three.

## 9. Archived evidence and reproducibility

**MACHINE-VERIFIED.** Every file under [`literature/pdfs/`](literature/pdfs/) has PDF magic `%PDF-`; SHA-256 digests for all archived PDFs and primary MATLAB files are in [`literature/manifest.sha256`](literature/manifest.sha256).

**CITED-DEPENDENCY.** Core archived primary or theorem-level sources include:

- [`liu_2023_arxiv_2306.08824v1.pdf`](literature/pdfs/liu_2023_arxiv_2306.08824v1.pdf)
- [`alweiss_huang_sellke_2024_entropy_kernel.pdf`](literature/pdfs/alweiss_huang_sellke_2024_entropy_kernel.pdf)
- [`agler_mccarthy_2000_complete_pick.pdf`](literature/pdfs/agler_mccarthy_2000_complete_pick.pdf)
- [`schoenberg_1938_metric_spaces.pdf`](literature/pdfs/schoenberg_1938_metric_spaces.pdf)
- [`fuglede_topsoe_2004_jensen_shannon.pdf`](literature/pdfs/fuglede_topsoe_2004_jensen_shannon.pdf)
- [`cuturi_et_al_2005_semigroup_kernels.pdf`](literature/pdfs/cuturi_et_al_2005_semigroup_kernels.pdf)
- [`hein_bousquet_2005_hilbertian_metrics.pdf`](literature/pdfs/hein_bousquet_2005_hilbertian_metrics.pdf)
- [`guillot_et_al_2026_entrywise_preservers.pdf`](literature/pdfs/guillot_et_al_2026_entrywise_preservers.pdf)
- [`boppana_2023_binary_entropy.pdf`](literature/pdfs/boppana_2023_binary_entropy.pdf)
- [`wakhare_2025_entropy_derivatives.pdf`](literature/pdfs/wakhare_2025_entropy_derivatives.pdf)

**CITED-DEPENDENCY.** Archived citation/survey PDFs include Das–Wu, Cambie, Lu–Raz, Colbert, Hu–Shi–Zhou, two Hachimori–Kashiwabara papers, and van der Hout–Roos; exact filenames and URLs are indexed in [`sources.json`](literature/metadata/sources.json).

**FAILED.** The CISS proceedings PDF, Antvelink thesis, and Gebendorfer paper were not archived because the serving sites denied access. Those failures and the exact attempted source URLs are retained in the query and screening ledgers.
