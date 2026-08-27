# Literature and originality audit: the explicit-constant and support-reduction route

**Search date:** 2026-08-25  
**Latest refresh:** the arXiv `all:"union-closed"` feed was queried again at
2026-08-26T15:10:15Z.  It still reported 103 records, no record posted after
the full-audit cutoff, and no theorem-level collision; see
[`arxiv_refresh_2026-08-26.json`](literature/arxiv_refresh_2026-08-26.json).

**Scope:** theorem-level prior art for an explicit union-closed frequency constant above
\(\psi=(3-\sqrt5)/2\), dependent-coupling formulations, finite-support reductions,
and computer-assisted certificates.  
**Machine-readable record:** [`literature/sources.json`](literature/sources.json),
[`literature/search_protocol.json`](literature/search_protocol.json), and
[`literature/citation_graph.json`](literature/citation_graph.json). The saved PDFs,
source bundles, and author code are content-addressed there by SHA-256.

> **Integration update (2026-08-26).** The independent mathematical audit
> found no fatal proof gap, and a fresh eight-slice direct replay passed all
> 488,465,854 events with exact expected hashes and tallies.  Thus the
> report's conditional classification is now instantiated: within this bounded
> search, the result is apparently the first explicit *certified* improvement
> over \(\psi\), not the first claimed improvement.  External review and
> universal priority remain open.

## Status vocabulary

The labels in this report have deliberately narrow meanings.

- **MACHINE-VERIFIED**: a deterministic local operation was actually run in this audit
  (for example, hashing a saved source). It does not mean that this audit replayed the
  local 488-million-box certificate.
- **HUMAN-AUDITED**: the cited theorem statement and the relevant proof step were read
  in the version named here.
- **COMPUTATIONAL-EVIDENCE**: floating-point optimization, plotting, gridding, or a
  high-precision numerical check, without an exhaustive proof certificate.
- **CITED-DEPENDENCY**: the claim is imported from the cited primary source or local
  artifact and was not independently re-proved here.
- **OPEN**: an assumption, proof-completeness issue, access gap, or search-coverage gap
  remains.
- **FAILED**: the displayed inference is false or the primary source itself records a
  failure. This label is not applied merely because a paper is unpublished.

## Executive conclusion

1. **HUMAN-AUDITED — quantitative record.** The last explicit constant in the audited
   literature with several complete independent proofs is
   \(\psi=(3-\sqrt5)/2=0.381966011250105151\ldots\), proved in the November 2022 papers
   of [Sawin v3, Theorem 1](https://arxiv.org/pdf/2211.11504v3),
   [Chase--Lovett v1](https://arxiv.org/pdf/2211.11689v1),
   [Alweiss--Huang--Sellke v4](https://arxiv.org/pdf/2211.11731v4), and
   [Pebody v1](https://arxiv.org/pdf/2211.13139v1). AHS was subsequently published as
   [EJC 31(3), P3.35](https://doi.org/10.37236/12232).

2. **FAILED / OPEN — the published decimal beyond \(\psi\).** Yu's published
   \(0.38234\) conclusion depends on Proposition 1's finite-support reduction. Its proof
   invokes AHS concavity outside the fixed-mean slice on which AHS proves it; a valid
   pair of Yu extreme points gives a concavity defect
   \(-0.024762714427119718\ldots\). Thus the proof of the reduction is invalid, although
   this does not by itself disprove the final functional inequality. Cambie's sharper
   \(c^*\approx0.382345533366703\) has a valid two-atom *upper obstruction* for Sawin's
   two-strategy route, but its matching lower verification is expressly numerical and
   graphical, not an exhaustive interval certificate. See the detailed audit below.

3. **CITED-DEPENDENCY / OPEN — Liu.** [Liu v1, Theorem 6](https://arxiv.org/pdf/2306.08824v1)
   gives an analytic, non-explicit strict improvement using conditionally-IID coupling,
   but it imports Proposition 2 and the optimizer/equality description of the
   Yu--Cambie baseline. [Theorem 13](https://arxiv.org/pdf/2306.08824v1#page=16) prints
   \(0.382709087918741\) only under two expressly numerical hypotheses: a projected
   positive-semidefiniteness assertion and a nine-parameter global-minimizer structure.
   Neither decimal is an unconditional, certified explicit theorem in the audited chain.

4. **CITED-DEPENDENCY — local result.** The local primary artifact
   [`PROOF.md`](PROOF.md) targets the exact rational
   \[
   t_{\rm cert}=0.3820660112501052
   =\frac{955165028125263}{2500000000000000}
   \ge \psi+10^{-4}.
   \]
   It claims a fixed-slice reduction to two pair-orbits and an Arb branch-and-bound
   certificate with replay traces. This literature audit did **not** perform the expensive
   replay; validity is conditional on the independent proof/certificate audit owned by the
   main integration task.

5. **HUMAN-AUDITED — bounded originality classification.** Conditional on that local
   proof audit passing, the explicit \(t_{\rm cert}>\psi\) theorem and the exhaustive
   replayable certificate are **apparently new** in the searched literature. The
   two-pair-orbit reduction is best classified as a **strengthening/new application**, not
   as an invention of the extreme-point method: fixed-moment concavity and few-atom
   extremizers already occur in AHS and Liu, and the abstract moment-set theorem is
   classical. The new-looking part is the pair-orbit formulation of the *exact* symmetric-
   coupling functional, the choice of the fixed \(B=1-\mathbb E p\) slice that restores
   concavity, and the resulting two-orbit/five-parameter certificate domain.

No conclusion in this report is based on title search alone.

## Dated search protocol and coverage

**HUMAN-AUDITED.** The following protocol was executed on 2026-08-25; exact query
strings and returned counts are preserved in
[`search_protocol.json`](literature/search_protocol.json).

1. The arXiv API query `all:"union-closed"`, sorted by submission date, returned 103
   records. Both pages (100 + 3 records) were screened. Entropy/coupling papers,
   explicit-constant claims, support/cardinality reductions, formalizations,
   computer-assisted work, and competing full-proof claims were opened beyond the
   title/abstract level.
2. The latest arXiv records were then pinned for all core papers and near-neighbours.
   This caught material changes such as Cambie v2 (2025-02-16), Sawin v3
   (2023-06-19), AHS v4 (2024-07-08), and Wakhare v2 (2025-01-14).
3. Crossref was searched for versions of record and DOI metadata; OpenAlex was queried
   separately for every arXiv and published node's citing works; zbMATH Open was queried
   by exact title, arXiv identifier, or DOI. Public author pages and cited code repositories
   were inspected.
4. Exact-phrase searches covered `pair-orbit`, fixed-moment/extreme-point terminology,
   the decimals \(0.3823455333667027\) and \(0.382709087918741\), and
   computer-assisted/interval/formalized variants. Hits were followed to primary papers,
   source bundles, or author repositories.
5. Primary PDFs and source/code files that could be fetched without credentials were
   saved under [`literature/papers/`](literature/papers/),
   [`literature/sources/`](literature/sources/), and
   [`literature/code/`](literature/code/). Their bytes and SHA-256 hashes are
   **MACHINE-VERIFIED** in [`sources.json`](literature/sources.json).

**OPEN — coverage limitations.** MathSciNet was not available without credentials;
zbMATH Open, Crossref, and OpenAlex supply the public metadata layer. OpenAlex splits
some preprint/published records and undercounts citations. Semantic Scholar returned
HTTP 429. The IEEE CISS PDF, Winkler and Bauer PDFs, the SIAM/ICM 2026 article, and the
MDPI version-of-record PDF were blocked or returned wrappers; the corresponding DOI
metadata and accessible arXiv/PMC versions were retained. Private manuscripts,
unindexed repositories, and sources not surfaced by the stated searches remain
outside the claim. Non-arXiv work after 2026-08-25 is also outside coverage;
the separately recorded arXiv refresh through 2026-08-26 found no later record.

## Primary-source table

| Primary source | Version and date audited | Exact theorem/assumptions relevant here | Support/coupling/computation | Audit status and durable copy |
|---|---|---|---|---|
| [Gilmer, *A constant lower bound*](https://arxiv.org/abs/2211.09055v2) | arXiv:2211.09055v2, 2022-11-28 | Theorems 1--2: iid samples; coordinate marginals at most \(0.01\); union-closed frequency \(0.01\). | Entropy tensorization/conditioning; no dependent coupling. | **HUMAN-AUDITED**. [PDF](literature/papers/arxiv-2211.09055v2.pdf). |
| [Sawin, *An improved lower bound*](https://arxiv.org/abs/2211.11504v3) | arXiv:2211.11504v3, 2023-06-19 | Theorem 1 proves \(\psi\). Section 2 is explicitly titled “Proof Sketch” and sketches existence of some \(\delta>0\), without an explicit value. | For the iid functional, a minimizer is reduced to \(\{v,1\}\). The dependent term uses a greedy max-entropy coupling and a convex combination. | **HUMAN-AUDITED** for \(\psi\); **OPEN** as a complete beyond-\(\psi\) theorem. [PDF](literature/papers/arxiv-2211.11504v3.pdf). |
| [Chase--Lovett, *Approximate union closed conjecture*](https://arxiv.org/abs/2211.11689v1) | arXiv:2211.11689v1, 2022-11-21 | Proves the \(\psi\) theorem and an approximate-union-closed extension. | Short entropy proof; cites the one-variable inequality as rigorously proved by AHS. | **HUMAN-AUDITED**. [PDF](literature/papers/arxiv-2211.11689v1.pdf). |
| [Alweiss--Huang--Sellke](https://doi.org/10.37236/12232) | arXiv:2211.11731v4, 2024-07-08; EJC version of record, 2024-09-20 | Theorems 1--2 prove \(\psi\). AHS minimizes \(F(\mu)=\mathbb E h(XY)-\mathbb E h(X)\) on the **fixed-mean** set \(M_\phi\). | Lemmas 5--6: concavity on \(M_\phi\), support at most 2, then \(\{0,x\}\). Appendix uses two floating tables with stated margins; source contains `numerics.py`, not interval arithmetic. | **HUMAN-AUDITED** theorem and proof; appendix outputs are **COMPUTATIONAL-EVIDENCE**. [published PDF](literature/papers/ahs-EJC-v31i3p35-published.pdf), [v4 PDF](literature/papers/arxiv-2211.11731v4.pdf), [code](literature/code/ahs-numerics-v4.py). |
| [Pebody, *Extension of a Method of Gilmer*](https://arxiv.org/abs/2211.13139v1) | arXiv:2211.13139v1, 2022-11-23 | Theorems 1--2 prove \(\psi\). Theorem 10 preserves mean and entropy while not increasing the iid product-entropy objective. | Repeated mass merging gives at most one nonzero support point, i.e. \(\{0,v\}\). | **HUMAN-AUDITED**. [PDF](literature/papers/arxiv-2211.13139v1.pdf). |
| [Yu, *Dimension-Free Bounds*](https://doi.org/10.3390/e25050767) | arXiv:2212.00658v2, 2023-05-05; *Entropy* 25(5), 767, 2023-05-08 | Theorem 1 is a general coupling criterion. Proposition 1/Corollary 1 claim reduction of Sawin's two-strategy functional to two symmetric pair-orbits; the paper then reports \(t=0.38234\), \(\alpha=0.035\), and numerical ratio \(1.00000889\). | Krein--Milman/Carathéodory plus a purported global concavity inference; cited code URL currently resolves to an empty template and the arXiv source has no code. | **FAILED** finite-support proof step; numerical value is **COMPUTATIONAL-EVIDENCE**, not a certificate. [v2 PDF](literature/papers/arxiv-2212.00658v2.pdf), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10217025/). |
| [Cambie, *Better bounds*](https://arxiv.org/abs/2212.12500v2) | arXiv:2212.12500v2, 2025-02-16 | Question 2 is the iid/max-entropy two-strategy inequality. Theorem 3 claims \(c^*\sim0.3823455\); §2.3 gives the two-atom ceiling obstruction and defining equations. | Marginal support claimed at most 3; §3.4 narrows to two two-variable cases using Yu. Final verification uses Maple optimization and plots. | Ceiling obstruction: **HUMAN-AUDITED**. Matching lower bound/proof completeness: **OPEN** with **COMPUTATIONAL-EVIDENCE**. [PDF](literature/papers/arxiv-2212.12500v2.pdf), [pinned worksheets](literature/code/), [computation PDFs](literature/papers/cambie-FinalComputation24.pdf). |
| [Liu, *Conditionally IID Coupling*](https://doi.org/10.1109/CISS59072.2024.10480167) | arXiv:2306.08824v1, 2023-06-15; CISS version of record, 2024-03-13 | Theorem 6 claims a non-explicit strict improvement over \(c^*\). Theorem 9 reduces one conditionally-IID problem to binary \(W\). Theorem 12 gives three paired atoms for sufficiently small \(l\). Theorem 13 gives \(0.382709087918741\) under two hypotheses. | Fixed-mean concavity for Theorem 9; projected PSD and nine-dimensional optimization at \(l=1\) checked by grids, eigenvalues, and random-start MATLAB. | Theorems 9/12: **HUMAN-AUDITED** in their stated domains. Theorem 6 has a **CITED-DEPENDENCY** on the Yu--Cambie baseline. Theorem 13 is **OPEN** and **COMPUTATIONAL-EVIDENCE**. [arXiv PDF](literature/papers/arxiv-2306.08824v1.pdf), [author code](literature/code/liu-frankl5.m). |
| [Boppana, *A Useful Inequality*](https://arxiv.org/abs/2301.09664v1) | arXiv:2301.09664v1, 2023-01-23 | Gives a calculus proof of the one-variable binary-entropy inequality used at the \(\psi\) barrier. | Analytic, one dimensional; no dependent-coupling optimization. | **HUMAN-AUDITED** as a non-overlapping proof ingredient. [PDF](literature/papers/arxiv-2301.09664v1.pdf). |
| [Winkler, *Extreme Points of Moment Sets*](https://doi.org/10.1287/moor.13.4.581); [Pinelis](https://arxiv.org/abs/1204.0249v1) | Winkler 1988; Pinelis arXiv v1 (2012), journal 2016 | Generalized-moment extreme measures have finite atomic representations; one probability normalization plus one moment yields the familiar two-atom ceiling under the required topological hypotheses. | Abstract convex/moment-set machinery; it does not supply concavity of the UC functional. | **HUMAN-AUDITED** via metadata, Pinelis Theorem 12, and local direct-perturbation comparison. [Pinelis PDF](literature/papers/arxiv-1204.0249v1.pdf). |
| [Bauer, *Minimalstellen…*](https://doi.org/10.1007/BF01898615); [Bru--de Siqueira Pedra](https://arxiv.org/abs/1610.03411v1) | Bauer 1958; Bru--Pedra arXiv v1 (2016) | Bauer's principle places extrema of suitable semicontinuous convex/concave functions on extreme points. | General principle only. | **CITED-DEPENDENCY**; not a novelty claim for the local work. [Bru--Pedra PDF](literature/papers/arxiv-1610.03411v1.pdf). |
| [Ho, *A generalization of Boppana's entropy inequality*](https://arxiv.org/abs/2601.19327v1) | arXiv:2601.19327v1, 2026-01-27; Lean repository pinned 2026-03-13 | Formalizes the generalized one-variable inequality and derives an approximate \(k\)-union-closed result. It does not certify an original-family constant above \(\psi\). | Lean 4 formalization of the analytic inequality. | **CITED-DEPENDENCY** for formalization scope; no overlap with the five-dimensional dependent-coupling certificate. [PDF](literature/papers/arxiv-2601.19327v1.pdf), [Lean source](literature/code/ho-entropy-2601.19327.lean). |

## Theorem-by-theorem comparison

### 1. The iid barrier and its complete proofs

**HUMAN-AUDITED.** The four November 2022 routes agree on the exact constant
\(\psi\), but their reductions are not interchangeable:

- [AHS Theorem 1 and Lemmas 5--6](https://arxiv.org/pdf/2211.11731v4#page=2)
  use concavity only after fixing the mean and reduce a marginal law to two points,
  one of them zero.
- [Sawin Lemma 3](https://arxiv.org/pdf/2211.11504v3#page=4) directly varies an
  optimizer of the iid functional under \(\mathbb E p\le u\) and gets support
  \(\{v,1\}\) (the complement convention accounts for the apparent difference).
- [Pebody Theorem 10](https://arxiv.org/pdf/2211.13139v1#page=4) preserves both
  first moment and entropy while merging nonzero atoms until only \(\{0,v\}\)
  remains.
- [Chase--Lovett](https://arxiv.org/pdf/2211.11689v1) packages the sharp binary
  entropy inequality into an approximate-union-closed theorem.

**HUMAN-AUDITED.** AHS's appendix is computer-assisted, but the same one-variable
inequality also has analytic proofs in Sawin and
[Boppana](https://arxiv.org/pdf/2301.09664v1). Therefore the status of \(\psi\)
does not depend on trusting the saved NumPy tables.

### 2. Sawin's dependent-coupling sketch

**HUMAN-AUDITED / OPEN.** [Sawin §2](https://arxiv.org/pdf/2211.11504v3#page=8)
combines the iid entropy with a greedily chosen max-entropy coupling. It sketches a
compactness/stability argument yielding some \(\delta>0\), and ends: “It would be
interesting to modify this argument to obtain an explicit value of \(\delta\).” The
paper's numbered Theorem 1 remains the \(\psi\) theorem. This is conceptual prior art
for the local coupling formulation, not prior art for an explicit certified decimal.

### 3. Yu's criterion and the failed finite-support inference

**HUMAN-AUDITED.** [Yu Theorem 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10217025/)
defines a general max-correlation coupling optimization \(\Gamma(t)\) and proves that
\(\Gamma(t)>1\) implies the union-closed frequency bound. The simpler Sawin mixture is
then encoded by \(\widehat\Gamma(t)\).

**FAILED.** In [Yu §4, proof of Proposition 1](https://arxiv.org/pdf/2212.00658v2#page=8),
a symmetric pair law on a finite grid with \(\mathbb E p\le t\) is decomposed into
extreme laws \(Q_i\), and the proof asserts
\[
 g\!\left(\sum_i\gamma_iQ_i,\alpha\right)
 \ge \sum_i\gamma_i g(Q_i,\alpha).
\]
The cited [AHS Lemma 5](https://arxiv.org/pdf/2211.11731v4#page=3) proves concavity
only on \(M_\phi\), where every measure has the *same* mean. Yu's extreme points in
the displayed decomposition can have different means.

A particularly small witness uses Yu's own extreme-point shapes. Let
\[
 t=0.382345533366702721,\quad \alpha=0.0356069,\quad
 Q_1=\delta_{(0.15,0.15)},
\]
\[
 Q_2=(1-\beta)\delta_{(0.10,0.10)}+\beta\delta_{(1,1)},
 \qquad \beta=\frac{t-0.10}{0.90},\qquad \gamma=\tfrac12.
\]
Then \(Q_1\) is a \(\beta=0\) extreme law with mean \(0.15<t\), and \(Q_2\) is a
boundary extreme law with mean exactly \(t\). Direct 80-digit evaluation gives
\[
 g(\gamma Q_1+(1-\gamma)Q_2,\alpha)
 -\gamma g(Q_1,\alpha)-(1-\gamma)g(Q_2,\alpha)
 =-0.0247627144271197180262\ldots .
\]
This is **COMPUTATIONAL-EVIDENCE** with a margin far from roundoff; the structural
reason for failure is the missing fixed-mean premise. The local script
[`yu_gap.py`](yu_gap.py) independently reproduced \(-2.476271\times10^{-2}\) in this
audit. The witness invalidates the proof step, not the possibility that a different
proof of the reduced or unreduced inequality exists.

### 4. Cambie's exact ceiling versus lower verification

**HUMAN-AUDITED.** [Cambie §2.3](https://arxiv.org/pdf/2212.12500v2#page=4) defines
\(b\) as the larger root in \((0,1)\) of
\[
 h(b)(2-h(b))=h(2b-b^2),
 \qquad a=\frac{1-h(b)}{2-h(b)},
 \qquad c^*=a+(1-a)b.
\]
It obtains \(b\approx0.329454738503037\),
\(a\approx0.0788772927059232\), and
\(c^*\approx0.3823455333667034\). The corresponding two-atom law makes both
strategy terms equal to the original entropy, so it is a valid upper ceiling for this
specific two-strategy formulation.

**OPEN / COMPUTATIONAL-EVIDENCE.** The other direction is not an interval proof:
Cambie v2 says the four-variable check is “slightly less rigorous,” identifies finite
precision and local/global-minimum issues, and says an even more rigorous analysis is
unnecessary. Section 3.4 reduces to two two-variable cases but sends the final claim to
Maple worksheets and “graphical confirmation.” The pinned
[`FinalComputation24`](literature/papers/cambie-FinalComputation24.pdf) likewise says
it solves optimization problems and plots the two-dimensional graph. No interval
bounds, exhaustive subdivision, proof trace, or checker are present.

**OPEN.** Question 2 is stated for means *less than* \(c\), whereas the proof of
Theorem 3 begins by assuming every frequency is “at most” \(c\) and later says the
conditional means are “less than” \(c\). This endpoint mismatch is harmless only with
an additional argument excluding equality. It matters for the local rational target,
which explicitly certifies the closed domain \(\mathbb E p\le t_{\rm cert}\).

### 5. Liu's conditionally-IID route

**HUMAN-AUDITED.** [Liu Definition 2 and Theorem 6](https://arxiv.org/pdf/2306.08824v1)
introduce protocols in which the two sequences are iid conditional on shared randomness.
This is genuinely broader than the iid component and distinct from Sawin's max-entropy
coupling. Theorem 6 is qualitative: it does not print a usable \(\delta\).

**CITED-DEPENDENCY / OPEN.** Its proof starts from Proposition 2, including the exact
optimality/equality description of \(c^*\), which Liu attributes to the Yu--Cambie
finite optimization and “independent numerical optimization.” Consequently the
analytic perturbation is important prior art, but the full theorem chain inherits the
baseline proof-completeness issue identified above.

**HUMAN-AUDITED.** [Liu Theorem 9](https://arxiv.org/pdf/2306.08824v1#page=10) is the
closest sound methodological predecessor to the local support reduction: after fixing
the mean, concavity and Krein--Milman reduce a conditionally-IID law to a mixture of two
product laws (binary auxiliary variable \(W\)). The component laws themselves remain
arbitrary.

**HUMAN-AUDITED / OPEN.** [Theorem 12](https://arxiv.org/pdf/2306.08824v1#page=14)
reduces each of two component laws to three paired atoms only for sufficiently small
\(l>0\). Section V sets \(l=1\) and replaces the needed PSD theorem by:

- a grid of spacing \(0.0004\) with a reported minimum projected eigenvalue
  \(-2.3685\times10^{-14}\);
- coefficient-matrix eigenvalue checks through truncation 90; and
- a nine-variable MATLAB search with random initializations.

The saved [`frankl5.m`](literature/code/liu-frankl5.m) uses unseeded `rand`, local
`sqp`, and 4000 starts in the pinned file. These are **COMPUTATIONAL-EVIDENCE**, not a
global certificate. [Theorem 13](https://arxiv.org/pdf/2306.08824v1#page=16) itself
correctly labels both the PSD assertion and global-minimizer structure as hypotheses.

### 6. The local route being classified

**CITED-DEPENDENCY.** According to [`PROOF.md`](PROOF.md), the local argument uses
Cambie's entropy-to-union-closed implication but supplies a different proof of the
functional inequality:

1. recast every symmetric coupling as a probability measure \(\nu\) on
   \(\Delta=\{(p,q):0\le p\le q\le1\}\);
2. decompose the iid OR-entropy as
   \(Q(\mu)=2B(\mu)L(\mu)-E(\mu)\), where \(E\) is a positive-semidefinite quadratic
   moment-defect term;
3. fix \(B(\nu)=\beta^*\) at a global minimizer, making \(BL\) linear on that slice
   and the exact objective concave;
4. use Bauer plus the one-moment extreme-point bound (and a direct three-set
   perturbation) to obtain at most two pair-orbits, hence at most four marginal atoms;
5. certify a lower relaxation on the resulting five parameters
   \((w,p_1,q_1,p_2,q_2)\) by Arb interval branch-and-bound and replay traces.

The functional, slice choice, and certificate are local contributions; Bauer,
Krein--Milman, and the one-moment atomic bound are cited classical tools.

## Support-reduction comparison

| Work | Optimization domain and fixed constraints | Concavity/extreme-point mechanism | Claimed support | Status | Relation to local method |
|---|---|---|---|---|---|
| [AHS Lemmas 5--6](https://arxiv.org/pdf/2211.11731v4#page=3) | Marginal laws \(\mu\), fixed mean \(\phi\) | iid objective concave on the fixed-mean slice | \(\le2\), then \(\{0,x\}\) | **HUMAN-AUDITED** | Known fixed-slice template; no dependent coupling. |
| [Pebody Theorem 10](https://arxiv.org/pdf/2211.13139v1#page=4) | Marginal laws; preserve mean and entropy | Pairwise mass merging | At most one nonzero atom | **HUMAN-AUDITED** | Strong iid-only reduction with two preserved scalars. |
| [Sawin Lemma 3](https://arxiv.org/pdf/2211.11504v3#page=4) | Marginal laws, \(\mathbb E p\le u\) | First variation plus convex/concave shape | \(\{v,1\}\) | **HUMAN-AUDITED** | Iid-only optimizer geometry. |
| [Yu Proposition 1](https://arxiv.org/pdf/2212.00658v2#page=8) | Symmetric pair laws, \(\mathbb E p\le t\) | Applies fixed-mean concavity across different means | Two pair-orbits (five-dimensional) | **FAILED** proof step | Same desired parameter count, but the global concavity premise is false. |
| [Cambie §§3.2--3.4](https://arxiv.org/pdf/2212.12500v2#page=7) | Marginal/coupling law for Question 2 | Mass-moving and curvature cases; final use of Yu | Marginal support \(\le3\), then two special joint forms | **OPEN** | Valuable structural restrictions, but final global check remains numerical and inherits Yu. |
| [Liu Theorem 9](https://arxiv.org/pdf/2306.08824v1#page=10) | Convex hull of conditionally-iid product laws; fixed mean | Fixed-mean concavity plus Krein--Milman | Binary mixture of two product laws | **HUMAN-AUDITED** | Closest sound precedent for fixing the mean before taking extremes; component laws not finitely supported. |
| [Liu Theorem 12](https://arxiv.org/pdf/2306.08824v1#page=14) | Two component laws, two additional linear constraints | Conditional PSD/concavity | Three paired atoms in each component | **HUMAN-AUDITED** only for sufficiently small \(l\); **OPEN** at \(l=1\) | Different conditionally-IID objective and larger parameterization. |
| [Local Theorem B-triple-prime](PROOF.md#3-support-reduction-theorem-b) | All symmetric pair-orbit laws; feasible mean inequality; then fixed \(B=\beta^*\) | Exact PSD-defect decomposition restores slice concavity; Bauer/Winkler/direct perturbation | Two pair-orbits, marginal \(\le4\), five parameters | **CITED-DEPENDENCY** pending main proof audit | **Strengthening/new application**: repairs the missing fixed-slice premise while retaining the full symmetric-coupling cost. |

**HUMAN-AUDITED.** The abstract “one moment implies at most two atoms” conclusion is
known from [Winkler](https://doi.org/10.1287/moor.13.4.581) and
[Pinelis](https://arxiv.org/pdf/1204.0249v1). Therefore any originality claim must be
attached to the UC functional's pair-orbit representation and slice concavity, not to
Bauer/Krein--Milman/Winkler themselves.

## Certificate-method comparison

| Source | What the machine does | Global/exact guarantees in the primary artifact | Status versus the local certificate |
|---|---|---|---|
| [AHS Appendix A](https://arxiv.org/pdf/2211.11731v4#page=7) | NumPy evaluates two finite tables; analytic derivative bounds propagate table margins. | Paper assumes output accuracy \(10^{-3}\) and says interval arithmetic *can* make it rigorous; downloaded code uses ordinary doubles. Independent analytic proofs cover the same \(\psi\) inequality. | **COMPUTATIONAL-EVIDENCE**; much smaller one-dimensional task, no proof trace. |
| [Yu v2](https://arxiv.org/pdf/2212.00658v2#page=5) | Numerical optimization of the claimed five-dimensional reduction. | Reports ratio \(1.00000889\); cited code is presently unavailable and the reduction proof fails. | **FAILED** chain plus **COMPUTATIONAL-EVIDENCE**. |
| [Cambie v2 and worksheets](https://github.com/StijnCambie/UCconjecture/tree/3fd901fb43d1da6aab05f941e0c8d61269c9951d) | Maple local minimization, high working precision, and plots for two cases. | No interval enclosure, exhaustive box cover, or independent checker; plots are used as global confirmation. | **COMPUTATIONAL-EVIDENCE** only. |
| [Liu §V and code](https://arxiv.org/pdf/2306.08824v1#page=15) | Dense floating eigenproblems, finite coefficient truncations, and random-start local nonlinear optimization. | Theorem 13 retains both conclusions as hypotheses. | **OPEN** and **COMPUTATIONAL-EVIDENCE**. |
| [Ho v1 / Lean](https://github.com/boonsuan/entropy-inequality/tree/93db124699c1ccc06a7b1333596c039bb37bad65) | Lean formalizes a one-variable generalized Boppana inequality. | Formalization is source-available, but does not address dependent coupling or an explicit \(>\psi\) constant. | **CITED-DEPENDENCY**; non-overlapping formal work. |
| [Local campaign and replay](PROOF.md#6-the-interval-certificate-5-d-branch-and-bound) | Arb ball arithmetic, branch-and-bound over five parameters, committed per-node traces, frozen collector/replayer. | Local artifacts claim an exhaustive cover at the exact rational target and replay of all eight traces. This audit checked source metadata, not the expensive replay. | If independently replayed successfully: **apparently new** proof-certificate application for this theorem. Currently **CITED-DEPENDENCY** here. |

## Citation graph and later reporting

The public citation graph is preserved with exact OpenAlex work IDs and query URLs in
[`citation_graph.json`](literature/citation_graph.json). OpenAlex counts are lower
bounds because it splits arXiv and published nodes.

```text
Gilmer (0.01)
   └── AHS / Chase--Lovett / Pebody / Sawin (psi, complete)
                                      └── Sawin dependent-coupling sketch (>psi, non-explicit)
                                              ├── Yu (finite optimization; failed reduction + numerics)
                                              └── Cambie (c* ceiling exact; lower side numerical/graphical)
                                                        └── Liu (conditionally IID)
                                                              ├── Theorem 6: qualitative >c*
                                                              └── Theorem 13: 0.382709..., two hypotheses

AHS/Boppana one-variable inequality
   ├── Wakhare (analytic generalization; reports Liu imprecisely)
   └── Ho (Lean formalization/generalization; no >psi UC constant)

Local route
   ├── imports Cambie's entropy-to-UC bridge
   ├── uses classical Bauer/Winkler/Pinelis extreme-point tools
   └── supplies a different fixed-slice pair-orbit reduction and interval certificate
```

The theorem-status propagation in later sources is mixed:

- **HUMAN-AUDITED.** [Cambie's 2023 survey v1](https://arxiv.org/pdf/2306.12351v1)
  accurately describes Yu's value as numerical and Cambie's conclusion as using
  three-dimensional plots, numerical minimization, and a precise two-atom case; its note
  added records Liu's later method.
- **FAILED.** [Lu--Raz v2](https://arxiv.org/pdf/2405.10639v2#page=2) says the current
  \(0.38271\) was “proven by Liu.” That conflicts with Liu's own Theorem 13, whose
  two hypotheses are explicit.
- **FAILED.** [Wakhare v2](https://arxiv.org/pdf/2312.14743v2#page=2) calls
  \(0.38237\) the current best and cites Liu. Liu prints \(0.382709\ldots\), not
  \(0.38237\), and prints it conditionally. The number is best treated as a digit
  transposition plus an omitted qualification, not an independent theorem.
- **OPEN.** OpenAlex shows that [Samotij's ICM 2026 invited lecture](https://doi.org/10.1137/25M180874X)
  cites Yu and published AHS, but the full text returned HTTP 403. Its treatment of the
  decimal could not be audited.
- **HUMAN-AUDITED.** [Ho v1 (2026)](https://arxiv.org/pdf/2601.19327v1) formalizes
  the Boppana inequality and still describes the optimized original-family constant as
  \(\psi\); it contains no subsequent dependent-coupling certificate.

No inspected citing work supplies a missing interval proof of Yu/Cambie, proves Liu's
nine-parameter global-minimizer hypothesis, or presents an explicit certified constant
between \(\psi\) and \(c^*\).

## Screened competing “full proof” claims

These records matter to a novelty audit, but they are not part of the audited
entropy/coupling route and were not accepted from titles alone.

- **FAILED.** [Scandone, arXiv:2302.03484v2](https://export.arxiv.org/api/query?id_list=2302.03484)
  is an entropy-route full-proof claim whose current arXiv metadata contains the author's
  comment: “There is a mistake in the proof of Proposition 1.2, at the end of page 4.”
- **OPEN.** [Demontis, arXiv:2405.03731v1 / DOI 10.33140/CRSM.03.02.04](https://doi.org/10.33140/CRSM.03.02.04)
  and [Schrader, arXiv:2501.03302v1](https://arxiv.org/abs/2501.03302v1) claim the
  full conjecture by unrelated routes. Both had zero OpenAlex citations at the cutoff;
  neither is used by the 2026 sources inspected here as a settled theorem. This report
  does not adjudicate those entire proofs, so they remain **uncertain**, not evidence that
  the local result is known.

## Bounded originality classification

The requested classification is intentionally component-wise.

| Local component | Classification through 2026-08-25 | Basis and qualification |
|---|---|---|
| Exact explicit theorem at \(t_{\rm cert}=0.3820660112501052>\psi\) | **Apparently new**, conditional on local proof/certificate replay | No audited primary source gives a complete explicit theorem above \(\psi\). Yu/Cambie claim a larger value but lack a valid complete chain; Liu's printed value is explicitly conditional. **OPEN** outside the stated search coverage. |
| Pair-orbit folding of a symmetric coupling and the max-entropy cost | **Known / reformulation** | Symmetric couplings and the same max-entropy cost are explicit in Yu Proposition 1 and Cambie Question 2. |
| General use of fixed-moment concavity, Bauer/Krein--Milman, and \(\le2\)-atom extreme points | **Known** | AHS, Liu Theorem 9, Winkler, and Pinelis are direct precedents. |
| Exact decomposition \(Q=2BL-E\) with PSD moment-defect \(E\), used to choose the fixed \(B\)-slice | **Apparently new / uncertain** | No matching decomposition was found in the primary texts, citing works, or exact-formula searches. Absence from a bounded search is not proof of universal novelty. |
| Two-pair-orbit/five-parameter support theorem for the exact all-symmetric-coupling functional | **Strengthening and apparently new application** | It reaches Yu's desired two-orbit count with the missing fixed-slice premise restored, is stronger in finite support than Liu Theorem 9 for a different feasible class, and is not Cambie's three-marginal-atom argument. The abstract extreme-point theorem itself is known. |
| Arb exhaustive certificate with committed traces and independent replay for an explicit \(>\psi\) target | **Apparently new application**, not a new numerical method | Interval branch-and-bound and Arb are standard; no prior UC source audited here provides an exhaustive certificate/replay chain for a beyond-\(\psi\) constant. Status remains **CITED-DEPENDENCY** until the main replay finishes. |
| The sharp two-strategy ceiling \(c^*\) and its two-atom obstruction | **Known** | Cambie §2.3, and Yu/Liu's summaries, already give the equations and obstruction. The local target is deliberately below it. |
| An unconditional theorem at Liu's \(0.382709\ldots\) | **Not claimed; OPEN** | Liu Theorem 13 keeps two hypotheses. The local UC certificate does not reach or establish that value. |

## Submission guidance from the audit

1. **HUMAN-AUDITED.** State “first explicit *certified* constant above \(\psi\) found
   in the literature search,” not “first claimed explicit constant above \(\psi\).” Yu and
   Cambie indisputably claimed larger decimals earlier.
2. **HUMAN-AUDITED.** Describe the support theorem as a new pair-orbit application and
   repair/strengthening of the finite-dimensional route; cite AHS Lemmas 5--6, Liu
   Theorem 9, Winkler, and Pinelis prominently.
3. **HUMAN-AUDITED.** Separate Cambie's valid ceiling obstruction from the incomplete
   lower verification. Do not say Cambie's numerical value is false.
4. **HUMAN-AUDITED.** Quote Liu Theorem 13's two hypotheses whenever mentioning
   \(0.382709087918741\); do not repeat the Lu--Raz/Wakhare shorthand.
5. **OPEN.** Make the final novelty sentence conditional on successful independent replay
   and on the explicitly bounded search protocol. Avoid an unqualified “largest known”
   claim until that replay and the proof audit are both closed.
