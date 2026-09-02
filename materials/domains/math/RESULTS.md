# Union-closed project: authoritative results index

**Current through 2026-08-29.** This is the entry point for current results.
[`PROGRESS.md`](PROGRESS.md) remains the newest-first laboratory ledger, and
[`uc/campaigns/README.md`](uc/campaigns/README.md) remains the campaign
inventory. Evidence labels are literal: **MACHINE-VERIFIED** is a named
executable fact; **HUMAN-AUDITED** is ordinary checked mathematics;
**COMPUTATIONAL-EVIDENCE** is non-proof numerical evidence;
**CITED-DEPENDENCY** is imported; **OPEN** is unresolved; **FAILED** means the
wording or inference is not established.

## Headline result: certificate-backed UC candidate

**Campaign K is CERTIFIED (2026-08-31), target rational
$t=0.3822660112501052\ge\psi+3\times10^{-4}$.** All 16/16 slices committed
COMPLETE with `residual = stack = budget_time = 0`, and the frozen collector
replayed all sixteen traces — **805,542,368 boxes re-derived from their exact
dyadic roots with every claimed discharge re-proved in Arb** (117,615 s
single-core) — then printed `COMPOSITE CERTIFICATE` and exited 0. **The
certified constant therefore strengthens from $\psi+10^{-4}$ to
$\psi+3\times10^{-4}$.** Verdict transcript:
[`collect_output.txt`](uc/campaigns/cert3_20260822T015019Z_41f2e119129446d3b1a936d66bbe7683_b49a21ec80ba/collect_output.txt)
(sha256 `e4e799a8261b8ffde7d8993d570eb7399816e7873433dd26bf2ebb732e454602`);
the earlier pre-completion rejection is preserved beside it as
`collect_output_2026-08-24_rejection.txt`. See the campaign row in
[Open problems](#open-problems-and-adopted-next-targets) and
[`uc/campaigns/README.md`](uc/campaigns/README.md).

**MACHINE-VERIFIED finite claim.** Campaign I's frozen collector replayed all
eight one-byte-per-node traces in Arb arithmetic and printed:

> COMPOSITE CERTIFICATE: Phi >= 0 on the feasible 5-parameter family at
> certified rational t = 0.3820660112501052 >= exact psi + 0.0001 over w in
> [1/2,1]. The proved orbit-swap symmetry covers w in [0,1].

Here `Phi` is the explicit relaxed pair-orbit functional
\(\Phi_{\rm rel}\), not the original infinite-dimensional functional.  The
campaign contains 488,465,854 nodes, with
`residual = stack = budget_boxes = budget_time = 0` for every slice. Exact
external pins, hashes, environment, and replay commands are in
[`uc/REPRODUCIBILITY.md`](uc/REPRODUCIBILITY.md); the line-by-line mathematical
audit is [`uc/AUDIT.md`](uc/AUDIT.md).
The fresh 2026-08-26 direct replay independently re-executed every trace from
the raw pinned inputs: eight exit-zero processes, exact stdout hashes/tallies,
488,465,854 total nodes, zero residuals, and identical pre/post input maps.
Summary SHA-256:
`6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`.
A separate hardened replay from fresh read-only staging also passed all eight
slices and aggregated to 488,465,854 nodes with zero residuals; composite
SHA-256:
`57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.

The common Python-rule limitation is now separately closed.  A secure
source-independent verifier replayed all eight traces from byte zero, derived
centered bounds by interval AD, rediscovered every face proof, and reproduced
the exact 488,465,854-node tally with zero residuals.  It shares the trace and
Arb primitives, and its remaining OS/custody assumptions are explicit.  The
externally trusted report-lock raw SHA-256 is
`af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
the accepted composite canonical SHA-256 is
`4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.

**HUMAN-AUDITED analytic chain.** Exact decomposition, compact optimal-coupling
attainment, fixed-mean support reduction, the Margin Lemma, and the
\(L=0\) endpoint case prove
\[
\Phi_{\rm rel}\ge0\Longrightarrow\inf_{\mathbb E p\le t}F(\mu)\ge0.
\]
The revised manuscript also reconstructs the sequential Bernoulli-coupling and
entropy-chain argument from Cambie's Question 2/Section 4. A separate 256-bit
Arb check proves the strict first-coordinate endpoint margin. Cambie is now a
source citation rather than an opaque implication.

**Candidate theorem, not end-to-end machine verification.** Combining those
human arguments with the replayed finite certificate gives the candidate
frequency bound \(0.3820660112501052\ge\psi+10^{-4}\). A bounded literature
audit found this apparently new as an explicit certified improvement over
\(\psi\), but external mathematical review and publication priority remain
open. The arithmetic certificate does not machine-check the surrounding
functional analysis or set-family proof.

## Second result: candidate proof of Liu's Hypothesis 1

**HUMAN-AUDITED theorem with MACHINE-VERIFIED algebra.** The complete audit is
[`LIU_H1/AUDIT.md`](LIU_H1/AUDIT.md), and the submission manuscript is
[`LIU_H1/paper/main.pdf`](LIU_H1/paper/main.pdf).  For
\[
R(s,t)=(1-s)(1-t)\,[st-(1+st)\ln(1+st)]
-[z+(1-z)\ln(1-z)],\quad
z=(1-s)(1-t)(1+st),
\]
the ordinary Taylor/Lorentz/Gram proof establishes
\[
\iint R(s,t)\,d\nu(s)d\nu(t)\le0
\]
for every finite real signed Borel measure \(\nu\).

On Liu's constrained subspace, the entropy and \(R\) quadratic forms agree.
The unprojected theorem is stronger than, not equivalent to, that restricted
identity; it implies Liu's Section V-A Hypothesis 1. The proof's finite
polynomial identities are independently checked without CAS by
[`LIU_H1/verification/independent_exact_checker.py`](LIU_H1/verification/independent_exact_checker.py).
Calculus, convergence, endpoint, integration, and signed-measure extension
remain human-audited.

Two fresh read-only falsification passes found no blocking or major
mathematical defect.  Their four minor findings are repaired in the current
12-page manuscript: the Schur-product citation locator, an explicit
fixed-domain PSD pullback sequence at the axes, joint
\((r,s,t)\)-continuity/measurability before parameter integration, and a
complete inventory of the three infinite Gram series.  The exact reports are
preserved under [`LIU_H1/verification/logs/`](LIU_H1/verification/logs/).

The degree-8/12/16 Arb computations provide certified negative
\emph{one-sided upper bounds} for compressed \(R\) eigenvalues via a finite
Loewner majorant; they do not prove the continuum theorem. Float64 spectra and
the reported generalized ratio \(0.9999999906\) are computational evidence,
not an exact zero mode or a proof of exact tightness.

**OPEN scope.** Liu's independent Section V-B nine-parameter
global-minimizer hypothesis remains open. Therefore
\(0.382709087918741\) remains conditional. The literature audit found no
earlier proof of Hypothesis 1 in the searched sources, but absence from a
bounded search does not establish universal priority.


## Third result: Gate B scalar closure-defect repair is impossible

**PROVED/CERTIFIED, with two combinatorially distinct bases and three
arithmetically disjoint certificate stacks.**  For the normalized 45-row family
\(\mathcal B\) from the \(n=7\) census, the 256-bit Arb verifier proves
\(A_+(\mathcal B)<-7/250\) and counts join-success probability \(17/81\).
A second standalone checker that imports nothing outside the Python standard
library, uses \(2^{-160}\) outward-rounded dyadic intervals, and relaxes every
nondegenerate Bellman action interval to \([0,1]\) independently proves
\(A_+(\mathcal B)<-1/40\).  Either bound suffices.

A separate 25-row \(n=6\) family \(\mathcal D\), reconstructed from the five
\(3+3\) cells \((0,0),(0,1),(1,0),(1,2),(2,1)\), gives a second infinite
construction.  Its exact closure success is \(181/625\).  The dyadic checker
evaluates all \(6!=720\) orders without a symmetry quotient and proves
\(A_+(\mathcal D)<-1/80\).  A different 256-bit Arb evaluator gives the stronger
consistent finite bound \(A_+<-0.0136\), and a direct plain-mask check of the
625-row square found zero Bellman product gap on five hostile orders.

A third stack removes every remaining numerical concession.  An exact rational
evaluator, standard library only, maximizes over the **exact** feasible action
interval — no clamp classification, no relaxation to \([0,1]\), no interval
library — and returns two-sided enclosures over **all** \(720\) and **all**
\(5040\) coordinate orders:
$$
A_+(\mathcal D)\in[-0.01367210773217777356185991566,\,
 -0.01367210773217777356185991299],
$$
$$
A_+(\mathcal B)\in[-0.02864918679446831781602698193,\,
 -0.02864918679446831781602697645],
$$
each of width below \(6\times10^{-27}\), so it re-derives the sharp bounds
\(A_+(\mathcal D)<-17/1250\) and \(A_+(\mathcal B)<-7/250\) at once.  Both Arb
values lie strictly inside these enclosures.  Its three certified primitives
are a bit-by-bit binary logarithm, ceiling square-root towers for \(2^x\), and
a concavity case split for the constrained entropy maximum.  The same
arithmetic re-attacks the product lemmas on 35 products over all 15,010 global
orders, with zero intersection failures and a worst fixed-order Bellman
discrepancy of exactly zero, replacing the previous float64 evidence.

For the Cartesian powers \(\mathcal F_k=\mathcal B^{\boxtimes k}\) and
\(\mathcal G_k=\mathcal D^{\boxtimes k}\), a fixed-order Bellman induction
proves
$$
A_+(\mathcal F_k)=kA_+(\mathcal B),\quad
\varepsilon_\vee(\mathcal F_k)=1-\left(\frac{17}{81}\right)^k,
$$
$$
A_+(\mathcal G_k)=kA_+(\mathcal D),\quad
\varepsilon_\vee(\mathcal G_k)=1-\left(\frac{181}{625}\right)^k.
$$
Both sequences stay normalized and cap/Reimer-admissible; the \(k\)-free
Reimer witnesses are \(45^5<2^{28}\) and \(25^5<2^{24}\).  Therefore either
sequence proves
$$
\boxed{c_{\rm cl}^\star=+\infty}.
$$
This resolves Gate B as stated and rules out every fixed scalar coefficient.
It does not refute a union-closed statement.  Both defects tend to one, so the
separate local-stability question with \(\varepsilon_\vee\to0\) remains open.

The divergence is also quantified, and the numerator's growth is now closed on
both sides with an improved constant.  Write \(\Lambda(n)\) for the supremum of
\(-A_+\) over admissible families on at most \(n\) coordinates.  Incidence is at
most \(n\lfloor2m/5\rfloor\) while Reimer demands at least \(m\log_2m/2\), so
**every** admissible family obeys the size ceiling \(\log_2m\le4n/5\) — exactly,
as the integer test \(m^5\le2^{4n}\) — and with \(Q,C_+\ge0\),
$$\frac{3}{640}\,n-\frac{3}{80}\ \le\ \Lambda(n)\ \le\ \frac{4n}{5},$$
so \(\Lambda(n)=\Theta(n)\).  The lower constant improved this round by exactly
\(75/64=1.1719\ldots\), from \(1/250\) to \(3/640\), using a cloned-coordinate
base at \(n=8\); no lemma in the product argument uses separation, and one
base-level Reimer check \(m^m\le2^{2I}\) certifies every power.  Only the
denominator is still open: under a positive defect floor \(\varepsilon_0\) the
ratio grows exactly of order \(n\), and the only escape is
\(\varepsilon_\vee\to0\).

### The local regime is now bracketed (2026-08-27)

**PROVED/CERTIFIED where stated; the local question stays OPEN.**  The escape
above is not merely unaddressed by the theorem, it is unreachable by its
mechanism: the closure defect multiplies through success, so every Cartesian
power satisfies \(\varepsilon_\vee\ge\varepsilon_\vee(\text{base})\) and no
power of either base has defect below \(444/625\).  Writing
\(e^\star=\min\{\varepsilon_\vee(\mathcal F):\mathcal F\ \text{admissible},\
A_+(\mathcal F)<0\}\):

Ten bases are registered and exactly certified; one additional finite symbolic
clone family is certified by exact all-order reweighting. Current records:

| record | value | base | separating |
|---|---|---|---|
| **lowest defect with \(A_+<0\)** | \(\mathbf{62/225=0.275556}\) | `n10 m15` (e27 descent; Main-audited exact witness at order \((8,1,2,3,7,0,5,4,6,9)\), \(A_+\in[-0.0034369252227487344,-0.0034369252227487344]\); supersedes the \(14/45\) normalized-core clone) | no |
| lowest defect, separating | \(1144/1875=0.610133\) | `n8lo` | yes |
| highest repair ratio | \(147/2536=0.0579653\) | `n8clone_hi` | no |
| highest repair ratio, separating | \(531/11360=0.0467430\) | `n8best` | yes |
| **best asymptotic slope** | \(\mathbf{3/640=0.0046875}\) | `n8clone_hi` | no |
| **best asymptotic slope, separating** | \(\mathbf{177/40000=0.004425}\) | `n8best` | yes |

Normalization (active *and* separating) is a search and reporting condition in
[`uc/gate_b/DEFINITIONS.md`](uc/gate_b/DEFINITIONS.md), not a condition in the
displayed supremum, so both conventions are reported.  Published values before
this round were defect \(1336/2025\), ratio \(0.0362591\) and slope \(1/250\);
all three are improved, and the slope twice — once with a cloned base and once,
by \(10.6\%\), without leaving the separating class.  Block-symmetric coverage at
\(n=8\) is complete: all four classes \(k\in\{1,2,3,4\}\) enumerated,
52,494,332 raw masks, 27,431 canonical families, 150 screen-negative, **147
certified negative**.

* **CERTIFIED (2026-08-29):** \(e^\star\le62/225=0.27555\ldots<2/5\) in the
  displayed cap/Reimer class, from the e27-seeded \(n{=}10\), \(m{=}15\)
  family \((0,32,64,128,176,256,272,304,432,527,591,623,719,991,1023)\)
  (independent Main audit: admissible, incidence \(60\ge R_{15}=30\),
  \(15^{15}\le2^{120}\), all degrees at cap 6, \([10]\) present, exact
  negative fixed order). Lowest-defect **admissible** family known:
  \(56/225=0.248889\) at the same corner — its every fixed order is
  exactly positive on a 6,001-order scan, so defect alone does not give
  negativity. The separating record remains \(e^\star\le1144/1875\);
  the prior \(14/45\) normalized-core clone witness stands superseded in
  the displayed convention and retained as the normalized-slice record.
  Under active separation the new records are **not** certified
  (normalized \(=\) FALSE on both new families, reported per the
  both-conventions discipline).
* **CERTIFIED, independent earlier base:** the \(n=7\) witness is a 45-row
  family with all degrees at cap 18 and exact enclosure
  \(A_+\in[-0.00085518533372314171100789579,-0.00085518533372314171100789029]\).
  Its automorphism group is trivial, so all 5040 orders carry distinct
  enclosures. Its powers give a third infinite construction with
  \(\varepsilon_\vee=1-(689/2025)^k\).
* **CERTIFIED (2026-08-27):** the first family to beat the published repair
  ratio.  A separating 75-row family on eight coordinates has
  \(A_+\le-0.027744001\) at defect \(468/625\), ratio \(0.0370512\ldots>0.0362591\ldots\);
  dropping separation, which `DEFINITIONS.md` records as a reporting condition
  rather than a condition in the supremum, a 70-row family reaches
  \(0.0580391\ldots\).  Sizes 70 and 75 are *infeasible* at \(n=7\)
  (\(7\cdot28=196<R_{70}=215\)): the eighth coordinate supplies the missing
  incidence, and in the 70-row families it is an exact duplicate of the first.
* **PROVED (2026-08-27):** two defect floors.  The cap forces
  \(\mathbb E|X\vee Y|\ge\frac85\bar s\) while a successful union has size at
  most \(M=\max_A|A|\), so \(\varepsilon_\vee\ge(\frac85\bar s-M)/(n-M)\); and
  downset counting gives \(\varepsilon_\vee\ge1-m^{-2}\sum_AN(A)^2\), attained
  with equality.  Hence the ratio is \(O(n)\) on every family *without* a
  dominant set, and reaching the local regime **requires** a set of size at
  least \(\frac85\bar s\ge\frac45\log_2m\).
* **PROVED (2026-08-27):** 25 and 45 are the largest admissible sizes at
  \(n=6,7\); all three \(n\le7\) bases attain them with every degree exactly at
  the cap, and **no set can be added to any of them**, checked over all \(2^n\)
  candidates.  Adding the missing unions — the only route that lowers the defect
  and raises \(\log_2m\) at once — is therefore empty at the certified bases.
* **PROVED (elementary):** failures come in pairs, so a non-union-closed family
  has \(\varepsilon_\vee\ge2/m^2\) and \(-A_+/\varepsilon_\vee\le m^2\log_2m/2\);
  reaching defect \(\varepsilon\) forces \(m\ge\sqrt{2/\varepsilon}\).
* **NUMERICAL:** minimizing \(A_+\) under a hard defect cap at \(m=45\) gives a
  monotone curve crossing zero between caps \(0.64\) and \(0.66\); the least
  \(A_+\) found rises to \(+0.0795\) at cap \(0.35\).  All 38 minimum-defect
  witnesses at \(n\in\{6,7\}\), \(m\ge15\) have \(A_+>0\), from \(+0.0588\) to
  \(+0.5607\).  The searched minimum defect of an admissible family never fell
  below \(2/9\) at any \(n\le8\), against a \(2/m^2\) combinatorial floor three
  orders of magnitude smaller.  So as the defect falls \(-A_+\) falls faster and
  the ratio *collapses*, from \(0.0362591\) at \(\mathcal B\) to \(0.0012962\) at
  the new base (certified lower bound \(27/21376\)).
* **BLOCKED, precisely:** a family with \(\varepsilon_\vee=0\) is union-closed
  with all degrees \(<m/2\), i.e. a Frankl counterexample.  So proving a positive
  universal defect floor at cap \(2/5\) — what would empty the local regime —
  implies every union-closed family has an element in more than \(2m/5\) of its
  sets, about \(0.018\) beyond the published frontier \(\psi=0.381966\ldots\).
  Frankl is verified for \(n\le12\) **[REPORTED]**, so no admissible family in
  the searched range has zero defect.
* **METHODOLOGICAL, diagnosed:** the four complete block-symmetric
  \(n=8\) classes contain 52,494,332 raw masks, 27,431 canonical admissible
  families and 150 screen-reported negatives. Exact re-evaluation certifies 147,
  finds four screen negatives exactly non-negative, and finds one family the
  screen called non-negative to be exactly **negative**. This is *not* float64
  imprecision: on the same order set
  the float and exact evaluators agree to \(10^{-14}\), and the exact orbit route
  agrees with all \(8!\) orders exactly.  The census averages over one order per
  \(S_k\times S_{8-k}\) pattern, which represents the objective only when the
  automorphism group *is* that block group — true at \(k=1\), false at \(k=2\),
  where \(|\mathrm{Aut}|=1440\) matches \(|S_2\times S_6|\) but only 120 elements
  preserve the block partition, so 28 block orders cover 12 of 28 true orbits.
  Any subset average lies within \(\alpha(\max_\pi C_{+,\pi}-\min_\pi C_{+,\pi})\)
  of the truth, at most \(5.2\times10^{-3}\) here, so the re-ranking threshold is
  \(+0.01\).  Both frozen targets for the first new bases were written from
  census values and both were rejected by the verifier before correction.

### The local \(2/5\) barrier is crossed; the zero limit remains open (2026-08-28)

**PROVED/CERTIFIED where stated; search coverage remains heuristic.**

* **Sharp clone reachability (PROVED).** For multiplicities \(r_i\), the
  collapsed first-appearance order has exact Plackett--Luce probability
  \[
  \Pr_{\mathbf r}(\pi)
  =\prod_k\frac{r_{\pi_k}}{\sum_{j\ge k}r_{\pi_j}}.
  \]
  Every cloned fixed-order value equals its collapsed original value. Geometric
  multiplicities concentrate on any chosen order, so the clone infimum/supremum
  equal the fixed-order minimum/maximum. A finite clone multiset is negative
  iff one original fixed order is negative.
* **A row-changing core reaches the searched defect floor exactly.** The
  corrected `--score min-a` search found the normalized core
  \[
  \mathcal V=(0,1,2,4,5,8,10,43,64,190,192,193,245,254,255)
  \subseteq2^{[8]}.
  \]
  It has degrees \((6,6,6,6,4,5,6,6)\), incidence \(45\ge R_{15}=30\), full
  set, and defect \(14/45\). Its uniform objective is positive,
  \(A_+\in[0.230029313333,0.230029313334]\), but order
  \((6,1,2,0,3,4,5,7)\) has exact enclosure
  \([-0.003540262562,-0.003540262561]\).
* **Finite negative family below \(2/5\) (CERTIFIED EXACTLY).** Ratio \(R=41\)
  gives multiplicities
  \[
  (2825761,4750104241,115856201,68921,1681,41,194754273881,1).
  \]
  Exact reweighting of all 40320 fixed-order enclosures gives a finite symbolic
  family with dimension \(199623130728\), \(m=15\), incidence
  \(1197738780965\), unchanged defect \(14/45\), and
  \[
  A_+\in[-0.000084288238,-0.000084288237]<0.
  \]
  Certificate:
  [`gate_b_subtwofifths_clone_rational_v1.json`](uc/gate_b/certificates/gate_b_subtwofifths_clone_rational_v1.json),
  verdict `PROVED_SUB_TWO_FIFTHS_NEGATIVE_EXACT_RATIONAL`. Sampled orders select
  candidates only; the consumed order and final average are exact two-sided
  rational enclosures.
* **Both conventions are explicit.** The core \(\mathcal V\) is active and
  separating; the cloned family repeats columns. Thus \(14/45\) is the displayed
  cap/Reimer record. The separating negative frontier remains \(1144/1875\).
* **The previous \(14/45\) obstruction remains valid, but is family-specific.**
  The earlier normalized core
  \((0,2,4,6,8,9,16,22,31,32,96,105,112,125,127)\) has every fixed order at
  least \(0.145914214376\), and its one genuinely new defect-free extension has
  all 40320 orders positive. Every defect-free extension sequence of that core
  remains positive; changing rows is what exposed the new negative order.
* **Reimer is essential and discontinuous.** Chase--Lovett [REPORTED] refute
  cap + normalization + dominant-set floors and have normalized Reimer deficit
  tending to \(h(\psi)/2-\psi>0.0977433528\) while
  \(\varepsilon_\vee\to0\). Hence no pair-defect-continuous approximate Reimer
  inequality with vanishing normalized error exists. At exact zero defect,
  Reimer's cited average-set-size theorem and the dominant set are automatic,
  so a positive floor would still prove the \(2/5\) frequency theorem.
* **Open endpoint.** The new family has fixed defect \(14/45\), not a sequence
  tending to zero. Further progress needs a still lower-defect row-changing core
  with negative fixed order, a separating amplification, or a Reimer-essential
  floor using stronger support information.

The complete statement, proof, ten registered exact bases, the new symbolic
clone certificate, direct-square audits, the 292-case float and 35-case exact
tensorization falsification searches, exact \(n=8\) class re-ranking, finite
controls, and paper draft are indexed in
[`uc/gate_b/README.md`](uc/gate_b/README.md); the current verified state alone is
in [`uc/gate_b/CURRENT_STATUS.md`](uc/gate_b/CURRENT_STATUS.md).

## Theorem and lemma status

Commands are from `math/`. Executable controls prove only their named finite
or algebraic claims; they do not automatically verify the adjacent human
proof.

| Name | One-line statement | Status | Machine checker / primary source | Run from `math/` |
|---|---|---|---|---|
| Campaign I certificate | $\Phi_{\rm rel}\ge0$ on the feasible five-parameter family at exact rational $t=0.3820660112501052\ge\psi+10^{-4}$. | **MACHINE-VERIFIED finite** | Frozen traces/replayer plus external lock and independent structural audit; see [`uc/REPRODUCIBILITY.md`](uc/REPRODUCIBILITY.md) | Commands and observed status in the reproducibility guide |
| Entropy-to-UC bridge | The closed-domain functional inequality implies a present element of frequency at least $t$. | **HUMAN-AUDITED** universal statement; **MACHINE-VERIFIED** exhaustively for every family on $n\le4$ coordinates | Self-contained proof in [`uc/AUDIT.md`](uc/AUDIT.md) and the revised paper (Cambie v2 is the source); exhaustion report [`entropy-bridge-exhaustive.json`](uc/verification/results/entropy-bridge-exhaustive.json) | `./.venv/bin/python -I -B uc/verification/entropy_bridge_exhaustive.py --max-coordinates 4`; 256-bit endpoint subcheck: `./.venv/bin/python -I -B uc/verification/cambie_bridge_strictness.py` |
| Exact notation bridge | $s^\star$ formulas agree; symmetrization/folding preserves marginal and cost. | **HUMAN-AUDITED** with exact finite schemas | [`uc/bridge_uc.py`](uc/bridge_uc.py) | `./.venv/bin/python -I -B uc/bridge_uc.py` |
| Lemma 1 (PSD series) | $T(u)=u+(1-u)\ln(1-u)=\sum_{n\ge2}u^n/[n(n-1)]$, hence $T(xy)$ is PSD. | **HUMAN-AUDITED** with exact controls | [`uc/decomposition.py`](uc/decomposition.py) | Script checks finite/exact consequences |
| Theorem A | $H(xy)=g\otimes g-a\otimes a-T(xy)$. | **HUMAN-AUDITED** with exact SymPy identity | [`uc/decomposition.py`](uc/decomposition.py) | `./.venv/bin/python -I -B uc/decomposition.py` |
| Theorem A′ | $Q=2BL-E$ with $E$ an explicit PSD moment-defect form. | **HUMAN-AUDITED** with exact/sampled controls | [`uc/reduction.py`](uc/reduction.py) | `./.venv/bin/python -I -B uc/reduction.py` |
| Defect identity and iid sharpness | $F=[2(1-\alpha)B-1]L+\alpha C-(1-\alpha)E$; at $\alpha=0$, the point-mass threshold is $\psi$. | **HUMAN-AUDITED** | [`uc/reduction.py`](uc/reduction.py) | Same command |
| Theorem B′ (superseded) | Keeping the symmetric coupling as variable gives attainment on at most three pair-orbits. | **HUMAN-AUDITED**; superseded | [`uc/decomposition.py`](uc/decomposition.py) | Script runs controls only |
| Theorem B″ (superseded) | The relaxed functional needs at most three pair-orbits after freezing one convex direction. | **HUMAN-AUDITED**; superseded | [`uc/margin_lemma.py`](uc/margin_lemma.py) | Script runs sampled controls |
| Theorem B‴ | The exact pair-orbit functional attains its minimum on at most two pair-orbits. | **HUMAN-AUDITED** | [`uc/AUDIT.md`](uc/AUDIT.md), [`uc/thmB3_proof.py`](uc/thmB3_proof.py); standard cited compactness/Bauer inputs | Executable controls are sampled, not the proof |
| Margin Lemma | For $L>0$, $F\ge L\Lambda$; $L=0$ is handled separately. | **HUMAN-AUDITED** | [`uc/AUDIT.md`](uc/AUDIT.md), [`uc/margin_lemma.py`](uc/margin_lemma.py) | Script controls are supporting evidence |
| Liu H1 residual theorem | The continuous kernel \(R\) is NSD for every finite signed Borel measure and implies Liu Section V-A on the constrained subspace. | **HUMAN-AUDITED**; algebra **MACHINE-VERIFIED** | [`LIU_H1/AUDIT.md`](LIU_H1/AUDIT.md), [`LIU_H1/paper/main.pdf`](LIU_H1/paper/main.pdf) | `./.venv/bin/python -I -B LIU_H1/verification/independent_exact_checker.py` |
| Historical $\operatorname{rh}$ corner conjecture | The dense scan proposed $\operatorname{rh}(x)\le2\min(x,1-x)$; that script alone is not a proof. | **NUMERICAL** | [`uc/lemma_rh.py`](uc/lemma_rh.py) (1,999-point discovery scan) | `./.venv/bin/python uc/lemma_rh.py` |
| Rigorous $\operatorname{rh}$ bound | $\operatorname{rh}(x)\le2x$ on $[0,1/2]$ by exact identities plus an Arb cover; $h\le1$ gives $\operatorname{rh}(x)\le2(1-x)$ on the right half. | **PROVED** | [`uc/lemma_rh_proof.py`](uc/lemma_rh_proof.py) | `./.venv/bin/python uc/lemma_rh_proof.py` |
| One-step barrier | The pure iid one-step inequality holds exactly through $p\le\psi$ and fails immediately above it. | **PROVED** | [`uc/onestep.py`](uc/onestep.py) | `./.venv/bin/python uc/onestep.py` |
| Theorem C (sink-pair / dual cap) | Any mean-constrained certificate that treats the coupled term through a linear Kantorovich potential is capped at $t\le\psi$, even with mean-dependent potentials. | **PROVED** | [`uc/dual_barrier.py`](uc/dual_barrier.py) | `./.venv/bin/python uc/dual_barrier.py` |
| Yu Krein–Milman-step refutation | Explicit extreme points violate the concavity inequality used in Yu Section 4; this refutes that proof step, not the proposed constant. | **PROVED** | [`uc/yu_gap.py`](uc/yu_gap.py) and independent 40-digit regression [`uc/tests/test_yu_counterexample.py`](uc/tests/test_yu_counterexample.py) | `./.venv/bin/python uc/yu_gap.py && ./.venv/bin/python uc/tests/test_yu_counterexample.py` |
| Closed form for $c^*$ | For the two-orbit obstruction, $h(b)(2-h(b))=h((1-b)^2)$ and $c^*=1-(1-b)/(2-h(b))\approx0.3823455333667027$. | **PROVED** (closed-form relations); decimal is a checked approximation | [`uc/coupling.py`](uc/coupling.py) | `./.venv/bin/python uc/coupling.py` |
| $n=7$ closure-defect falsifier | A normalized 45-row cap/Reimer family has $\varepsilon_\vee=64/81$, $A_+\le-0.0286491867944683\ldots$, and $A_++\frac1{50}\varepsilon_\vee\le-0.0128467176586658\ldots<0$; it is explicitly non-UC. | **PROVED/CERTIFIED FINITE** | Complete declared block-class producer [`uc/shapley_n7_block_symmetric.py`](uc/shapley_n7_block_symmetric.py) plus independent 256-bit Arb checker [`uc/shapley_n7_falsifier_cert.py`](uc/shapley_n7_falsifier_cert.py) | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python -B uc/shapley_n7_falsifier_cert.py` |
| Gate B Cartesian-power obstruction | Two normalized cap/Reimer sequences independently satisfy \(A_+(\mathcal F_k)=kA_+(\mathcal B)\), \(\varepsilon_\vee(\mathcal F_k)=1-(17/81)^k\), and \(A_+(\mathcal G_k)=kA_+(\mathcal D)\), \(\varepsilon_\vee(\mathcal G_k)=1-(181/625)^k\); each has diverging \(-A_+/\varepsilon_\vee\), hence \(c_{\rm cl}^\star=\infty\). | **HUMAN-AUDITED; two bases MACHINE-CERTIFIED** | Mathematical proof [`uc/gate_b/PROOF.md`](uc/gate_b/PROOF.md), Arb checker [`uc/gate_b/verify_gate_b.py`](uc/gate_b/verify_gate_b.py), dependency-free \(n=7\) checker [`uc/gate_b/verify_gate_b_dyadic.py`](uc/gate_b/verify_gate_b_dyadic.py), and all-720-order \(n=6\) checker [`uc/gate_b/verify_gate_b_n6_dyadic.py`](uc/gate_b/verify_gate_b_n6_dyadic.py) | Run the three verifier commands in [`uc/gate_b/README.md`](uc/gate_b/README.md) |
| Gate B product lemmas | Fixed-order \(Q_\pi\) and \(C_{+,\pi}\) are additive under disjoint-block products and closure success multiplies. | **PROVED**; symbolic algebra, a 292-case float search, a 35-case/15,010-order **exact rational** search with worst gap exactly zero, and two non-factored-square searches **MACHINE-CHECKED** | [`uc/gate_b/PROOF.md`](uc/gate_b/PROOF.md), [`uc/gate_b/audit_tensorization_exact.py`](uc/gate_b/audit_tensorization_exact.py), [`uc/gate_b/audit_tensorization.py`](uc/gate_b/audit_tensorization.py), [`uc/gate_b/audit_square_direct.py`](uc/gate_b/audit_square_direct.py), [`uc/gate_b/audit_n6_square_direct.py`](uc/gate_b/audit_n6_square_direct.py) | Run the four audit commands in [`uc/gate_b/README.md`](uc/gate_b/README.md) |
| Gate B exact rational base certificate | Exact two-sided rational enclosures of the true \(A_+\), over the exact feasible action set and all \(720\)/\(5040\) coordinate orders, give \(A_+(\mathcal D)<-17/1250\) and \(A_+(\mathcal B)<-7/250\) with width \(<6\times10^{-27}\) and contain both Arb values. | **CERTIFIED EXACTLY** (standard library only; no interval class, no clamp classification, no relaxation) | [`uc/gate_b/verify_gate_b_rational.py`](uc/gate_b/verify_gate_b_rational.py), [`uc/gate_b/certificates/gate_b_unbounded_rational_v1.json`](uc/gate_b/certificates/gate_b_unbounded_rational_v1.json) | `nice -n 10 /Library/Developer/CommandLineTools/usr/bin/python3 -B uc/gate_b/verify_gate_b_rational.py --bases n6,n7 --orders all` |
| Gate B sub-\(2/5\) clone certificate | A normalized eight-coordinate \(m=15\) core at defect \(14/45\) has an exact negative fixed order. The sharp Plackett--Luce clone law and ratio-41 multiplicities give a finite cap/Reimer family with the same defect and \(A_+\in[-0.000084288238,-0.000084288237]\). | **PROVED / CERTIFIED EXACTLY**; displayed convention, non-separating clone | [`uc/gate_b/PROOF.md`](uc/gate_b/PROOF.md) Propositions 25--26; [`audit_clone_limit.py`](uc/gate_b/audit_clone_limit.py); [`gate_b_subtwofifths_clone_rational_v1.json`](uc/gate_b/certificates/gate_b_subtwofifths_clone_rational_v1.json) | `nice -n 19 ./.venv/bin/python -B uc/gate_b/audit_clone_limit.py` |
| Gate B growth rate | \(c_{\rm cl}^\star(n)>n/250-7/250\), and \(-A_+\le\log_2m\le n\), so on families with \(\varepsilon_\vee\ge\varepsilon_0\) the ratio is \(\Theta(n)\). | **PROVED**; instantiations exact-integer **MACHINE-CHECKED** | [`uc/gate_b/PROOF.md`](uc/gate_b/PROOF.md) corollary, `asymptotic` block of the rational certificate | `nice -n 10 ./.venv/bin/python -B uc/gate_b/test_gate_b.py` |
| Liu H2 zero-support boundary layer | With $y_0=1/32$, every mean-feasible nine-parameter point satisfies $\partial\,\mathrm{gap}/\partial b_j=w_jG_j$ and $\mathrm{gap}(V)\ge\mathrm{gap}(V\vert_{\text{layer}:=0})+m_0\sum_{j}w_jb_j$ with $m_0\ge0.1876317632$; the layer reduces to the face $b_j=0$, where $\lvert h'''\rvert\le1022.94$ for supports in $[y_0,1-y_0]$, and admits raw-gap $\kappa\le0.8730627$ at component-mean deviation $1/10$, above the raw-gap smooth ceiling $0.3871250878$ (the $0.7004675$ printed by `liu9_tube.py` is in $(\Phi-1)/\mathrm{dist}^2$ units; $\mathrm{gap}=EHX(\Phi-1)$ and $EHX(P^\ast)=p\,h(x)=0.5526667300$). | **HUMAN-AUDITED** derivation with **MACHINE-VERIFIED** Arb constants, 217,600 interval cells, and a 3,240-point falsification pass | [`uc/liu9_boundary_layer.py`](uc/liu9_boundary_layer.py), report [`liu9-boundary-layer.json`](uc/verification/results/liu9-boundary-layer.json), mutation tests [`uc/verification/test_liu9_boundary_layer.py`](uc/verification/test_liu9_boundary_layer.py) | `./.venv/bin/python -I -B uc/liu9_boundary_layer.py` |
| Liu H2 mirror layer | For $t=1-b_j\le t_0$ the singular coefficient is $1-2(1-\beta)W_1-2\beta A_1$; a tube forces $W_1\le\rho^2/g(1-t_0)$ with $g(b)=[b(b-x)]^2$, so $G_j\ge M\log(1/t)-K_1$ and raising the mirror atoms to $b_j=1$ gains $m_1\sum_jw_jt_j$ while raising the mean. | **PROVED** at $(\rho,t_0)=(1/10,1/64)$, raw-gap $\kappa\le0.393693$; **REFUTED** at $(1/10,1/32)$ | [`uc/liu9_mirror_layer.py`](uc/liu9_mirror_layer.py), report [`liu9-mirror-layer.json`](uc/verification/results/liu9-mirror-layer.json), mutation tests [`uc/verification/test_liu9_mirror_layer.py`](uc/verification/test_liu9_mirror_layer.py) | `./.venv/bin/python -I -B uc/liu9_mirror_layer.py` |
| Liu H2 smooth chart | On the active-mean chart both $\mathrm{gap}$ and $\mathrm{dist}^2$ vanish to second order, so the mean-value Hessian form gives $\mathrm{gap}\ge\kappa\,\mathrm{dist}^2$ with $\kappa$ the largest value making the pencil $\mathrm{Hess\,gap}-\kappa\,\mathrm{Hess\,dist}^2$ PSD. No third derivative is used. | **PROVED** on $\lvert s-x\rvert,\lvert d\rvert\le1/2048$ for every $q\in[1/4,3/4]$, with $\kappa\ge234627/1048576$ = 57.8% of the ceiling; **OPEN** at tube-relevant radius | [`uc/liu9_smooth_chart.py`](uc/liu9_smooth_chart.py), report [`liu9-smooth-chart.json`](uc/verification/results/liu9-smooth-chart.json), audit [`uc/LIU_H2_INGREDIENT_I.md`](uc/LIU_H2_INGREDIENT_I.md) | `./.venv/bin/python -I -B uc/liu9_smooth_chart.py` |
| Liu H2 residual frontier | All 79 residual complement boxes are covered by the zero-face hypothesis $2(1-\beta)\mu-1\ge0.11105875229>0$ (70) or are mirror-only and entirely outside the tube (9). | **MACHINE VERIFIED**, and it corrects an earlier count: 61 boxes are coordinate-identical across the two localized runs, not nine | [`uc/liu9_survivors.py`](uc/liu9_survivors.py), report [`liu9-survivor-classification.json`](uc/verification/results/liu9-survivor-classification.json) | `./.venv/bin/python -I -B uc/liu9_survivors.py` |
| Liu H2 endpoint stratum | Inner core $\lvert d\rvert\le1/32$ has $q$-uniform $\kappa=1/3$; the endpoint annulus $1/32\le\lvert d\rvert\le1/4$ with $\min(q,1-q)\le1/4096$ has $\kappa=1/20$, and the seam holds under $\rho^2\le p^2\varepsilon_{\rm sm}^2q_\ast(1-q_\ast)$. | **PROVED** on the pure-$d$ chart; **CONDITIONAL** on an interior chart at $\varepsilon_{\rm sm}=1/32$; $q\mapsto1-q$ alone **REFUTED** | [`uc/liu9_qdegenerate.py`](uc/liu9_qdegenerate.py), report [`liu9-qdegenerate.json`](uc/verification/results/liu9-qdegenerate.json), tests [`uc/verification/test_liu9_qdegenerate.py`](uc/verification/test_liu9_qdegenerate.py) | `./.venv/bin/python -I -B uc/liu9_qdegenerate.py` |
| Liu H2 centered chart form | A third-order interval jet bounds each Hessian entry by its exact centre value plus the third-derivative enclosure times the box radius, legitimate because the two boundary strata are separately closed and every chart support stays in $[0.659,0.722]$. | **PROVED**: inflation constant $811.29\to72.20$, certified radius $8\times$ larger at $1/256$ with $\kappa\ge4119063/33554432$; $1/128$ **REFUTED** | [`uc/liu9_chart_centered.py`](uc/liu9_chart_centered.py), report [`liu9-chart-centered.json`](uc/verification/results/liu9-chart-centered.json), tests [`uc/verification/test_liu9_chart_centered.py`](uc/verification/test_liu9_chart_centered.py) | `./.venv/bin/python -I -B uc/liu9_chart_centered.py` |
| Liu H2 chart cover to radius $1/32$ | The pencil $\mathrm{Hess\,gap}-\kappa\,\mathrm{Hess\,dist}^2$ is positive semidefinite on $\lvert s-x\rvert,\lvert d\rvert\le1/32$ for every $q\in[1/4096,4095/4096]$ at $\kappa\ge4119063/33554432$, over $393{,}216$ exhaustively abutting cells; weakest determinant margin $2.48\times10^{-5}$. A cover suffices because PSD is pointwise and the box is convex and contains the centre, so no radial quadrature is needed; Step A (the PSD statement) is pointwise; Step B, which turns it into $\mathrm{gap}\ge\kappa\,\mathrm{dist}^2$, inherits a deficit $\ge-7.361\times10^{-69}$ that is an **all-$q$ enclosure, not a sampled bound**: the centre value and $\partial_s$ come from a $q$-free jet, and $\partial_d$ is exactly zero for every $q$ because each term carries a $q$-free multiplier times $(1-q)(-q)+q(1-q)$, the zero polynomial in exact `Fraction`. The residual is Liu's numerically determined root of (87)-(90); the binding constraint was the $q$ range, not the radius, and it needs octave cells of relative width $1/64$ because $\mathrm{gap.hdd}$ and $\mathrm{dist.hdd}$ carry an exact factor $q(1-q)$. | **MACHINE VERIFIED** with exact-`Fraction` coverage | [`uc/liu9_chart_cover.py`](uc/liu9_chart_cover.py), 26 mutation tests [`test_liu9_chart_cover.py`](uc/verification/test_liu9_chart_cover.py), report [`liu9-chart-cover.json`](uc/verification/results/liu9-chart-cover.json) | `./.venv/bin/python -I -B uc/liu9_chart_cover.py --radius 1/32 --delta 1/512 --q-min 1/4096 --q-pieces 64` |
| Liu H2 nine-variable extension | The chart result does **not** extend to an unrestricted ambient nine-variable neighbourhood: $(p-\varepsilon,\varepsilon,q,x,y,0,x,y,0)$ at $(y,\varepsilon,q)=(1/32,1/1024,1/2)$ gives raw gap $-5.5338\times10^{-4}$, with exact `Fraction` bounds putting the linear coefficient below $-1/2$. That point is mean-INFEASIBLE, so Hypothesis 2 is untouched; the mean-preserving control was positive at all 65 scanned $y$ (0 negative) while the ambient split was negative at 46 -- that scan is now **superseded** by the universal enclosure in [`uc/liu9_transverse.py`](uc/liu9_transverse.py), which also explains the scan's tiny minimum near $x$ as a **double root**. | **REFUTED** (ambient). The feasible half-space is no longer conditional: **superseded for the one-new-support-per-component family by the two transverse rows below**, which prove it on the whole interval and, via the decoupling identity, on all of $(0,1)^2\times[0,1]$ | [`uc/liu9_ninevar.py`](uc/liu9_ninevar.py), 8 tests [`test_liu9_ninevar.py`](uc/verification/test_liu9_ninevar.py), report [`liu9-ninevar.json`](uc/verification/results/liu9-ninevar.json) | `./.venv/bin/python -I -B uc/liu9_ninevar.py` |
| Liu H2 transverse first variation | The mean-preserving transverse coefficient $L(y)=2pK(x,y)-h(y)-(y/x)A-\kappa y^2(y-x)^2$ satisfies $L(y)\ge-4.562\times10^{-136}$ for **every** $y\in(0,1)$, with a **double root** at $y=x$. $L(x)=0$ is structural (a rational-function identity in exact `Fraction`, entropy values as free symbols); $L'(x)=0$ needs Liu's eq. (87) *and* the $\beta$-system (89)-(90); $L''(x)=1.2289613682$. Five exactly abutting regions: a $y\log(1/y)$ layer with $\mu(u)\in[1-u,1]$ by termwise series comparison plus one tail cell containing $y=0$; $11{,}777$ bulk cells (weakest $7.32\times10^{-3}$); a $4{,}096$-cell window with $L''\ge0.4975$. Replaces the previous 65-point scan with a universal enclosure. | **PROVED** (deficit is the *square* of Liu's root residual) | [`uc/liu9_transverse.py`](uc/liu9_transverse.py), 35 mutation tests [`test_liu9_transverse.py`](uc/verification/test_liu9_transverse.py), independent derivation [`independent_transverse_check.py`](uc/verification/independent_transverse_check.py), report [`liu9-transverse.json`](uc/verification/results/liu9-transverse.json) `sha256 153fb77bbe2ae4eff027c9d3d80c9e1d905c4480db288be386649c21671534f8` | `./.venv/bin/python -I -B uc/liu9_transverse.py` |
| Liu H2 asymmetric decoupling | Because the objective's mean is the $q$-weighted mixture, unequal supports $y_1\ne y_2$ are mean-feasible with correction $\bar y=(1-q)y_1+qy_2$; the perturbation measure is then exactly the $q$-convex combination of the two symmetric ones, and every term is linear or bilinear in the measure, so $L_{\text{asym}}(y_1,y_2,q)=(1-q)L(y_1)+qL(y_2)$. The positivity bound extends to all of $(0,1)^2\times[0,1]$ with **no** additional certification. | **PROVED** (measure algebra); identity **MACHINE VERIFIED** to $9.68\times10^{-32}$ over 8 configurations | [`uc/liu9_transverse.py`](uc/liu9_transverse.py) | `./.venv/bin/python -I -B uc/liu9_transverse.py` |
| Liu H2 exact second order | For simultaneous mean-preserving insertion, the pencil is exactly $\varepsilon\bar L+\varepsilon^2P_2$: epsilon changes masses only, so no $\varepsilon^2\log(1/\varepsilon)$, higher term, or remainder exists. $S=2L+Q\ge0$ on $[0,1]$ over 1,666 cells, with $Q/y^2\ge0.091066$ at zero and $S''\ge1.406691$ at the double root. | **PROVED** | [`uc/liu9_second_order.py`](uc/liu9_second_order.py), 21 mutations, independent verifier, report [`liu9-second-order.json`](uc/verification/results/liu9-second-order.json) `sha256 0b1b5bfc11313af7c97224e06b14532fd3da647141bc6b1be7045cde536d2958` | `./.venv/bin/python -I -B uc/liu9_second_order.py --output uc/verification/results/liu9-second-order.json` |
| Liu H2 simultaneous insertion interaction | The second-order distance contains $\varepsilon^2q(1-q)(y_1-y_2)^2$ and the cross term is not sign-definite. A divided-difference Arb cover plus exact quadratic minimization in $q$ proves $S_{\rm asym}\ge0$ on $[0,1]^2\times[0,1]$: 1,388,611 cell pairs, 544,770 with negative interaction. Hence the exact pencil is nonnegative for every $0\le\varepsilon\le1/2$. | **PROVED**; dropping the distance cross factor and first-order decoupling are rejected mutations | Same module/report | Same command |
| Liu H2 strict mean-feasible normal | For mean excess $\delta\ge0$, exact bilinearity gives $P(\varepsilon,\delta)=P(\varepsilon,0)+\delta[h(x)+2\varepsilon\bar R]+\delta^2[K(x,x)-\kappa x^2]$. Arb proves the linear coefficient $\ge0.235013$ and quadratic coefficient $\ge0.633528$ on all 1,666 cells; component-mean cross terms cancel exactly. | **PROVED:** composition covers the full half-space $\mathrm{mean}\ge px$, not only the active boundary | `liu9_second_order.py`, 21 mutations | Same command |
| Liu H2 piecewise local tube | Relabel every two-way atom split so the smaller fragment has $\varepsilon\le1/2$; zero-mass, coincident-support, permutation and inactive-law coordinates are quotient gauges. The second-order theorem composes with the zero-support, mirror, chart and q-degenerate certificates. | **PROVED:** $\rho=1/1701$, near-maximal q-degenerate seam | SSOT [`uc/liu9_tube.py`](uc/liu9_tube.py), 7 integration tests, report [`liu9-piecewise-tube.json`](uc/verification/results/liu9-piecewise-tube.json) `sha256 4651ad293d5205a9ba416a3474fbc0b4ee6e553c762d4790f2ce7ccc83612742` | `./.venv/bin/python -I -B uc/liu9_tube.py --rhos 1/1701 --skip-complement` |
| Liu H2 whole insertion-segment pencil | At $y_1=y_2=1$, $\varepsilon=px$, numerator and $E H(X)$ both vanish on the $\{0,1\}$ law, so raw gap is structurally zero while the pencil is $\le-0.007245231$. | **REFUTED** (pencil only); **not** a raw-gap counterexample | `liu9_second_order.py` endpoint certificate | Same command |
| Liu H2 conservative complement frontier | The deterministic $\rho=1/1728$ cover remains a conservative decomposition of the complement of the proved $1/1701$ tube. Its old 49,767-box unfinished frontier is historical computational evidence; the former q=1 subproblem is now closed analytically rather than by subdivision. | **SUPERSEDED as a q=1 blocker; full block inequality remains OPEN** | [`liu9-complement-residual-rho1-1728-100k.json`](uc/verification/results/liu9-complement-residual-rho1-1728-100k.json) `sha256 494afc93c37eefb07cd0a8e0c982746db020e1387e0c93ab7b731fadd4524b69` | `./.venv/bin/python -I -B uc/liu9_residual.py --rho 1/1728 ...` |
| Liu H2 near-maximal q seam | $q_*=1/2254$, interior $q_{\min}=1/4096$, cutoff $1/32$ certify $\rho=1/1701$. The current seam proof ceiling is $5.8815684\times10^{-4}$, so this nearly saturates the method and refutes an orders-of-magnitude gain from $q_*$ alone. | **PROVED** | [`liu9-qdegenerate-maximal.json`](uc/verification/results/liu9-qdegenerate-maximal.json) `sha256 353131d9124516841fe32da2a8f09829d401f9304890a604d91ff544ca693107`, 5 mutations | `./.venv/bin/python -I -B uc/liu9_qdegenerate.py --q-star 1/2254 --tube-radius 1/1701 ...` |
| Liu H2 universal q=1 theorem | The exact endpoint blocker is closed with raw gap $\ge0.01232508028745$. More strongly, for every probability law $P$ of mean $M\ge m$, the size-biased identity $F(P)=M\mathbb E_{Q\otimes Q}C_M$ and a pointwise certificate $C_m\ge0$, $\partial_M C_M\ge0$ prove $F(P)\ge0$. The proof uses 524,800 bulk boxes, exact small-support and infinite-tail inequalities, and a positive-definite optimizer Hessian cover. | **PROVED**, arbitrary laws (stronger than the 230-box quotient) | [`uc/liu9_endpoint_support.py`](uc/liu9_endpoint_support.py), [`uc/liu9_size_biased.py`](uc/liu9_size_biased.py); reports `sha256 6510dc6c8b55b0132ee6993cee7f0141a11046013c1de5e0f06b51b1b71df843`, `sha256 eb645792526578d115da8356f7d9378812164273900a60607691a8f32661f19a`; independent endpoint `sha256 7b97e2d4e3936c248ef70f76431fbcab082186cf781a810341c9e59317221837` | `./.venv/bin/python -I -B uc/liu9_endpoint_support.py`; `./.venv/bin/python -I -B uc/liu9_size_biased.py --certify` |
| Liu H2 exact full block reduction | The inward-q dependence is the exact polynomial $G_0+(1-q)G_1+(1-q)^2G_2$. The $m$-frozen shortcut is false for a mean-feasible $(\delta_x,\delta_1)$ family, although its true raw gap remains positive. **Superseded display note (2026-08-29):** the $2\times2$ block-kernel display $B_{00}=((1-q)^2/M)C_M+\beta q(1-q)k_\pi$, $B_{11}=(q^2/M)C_M+\beta q(1-q)k_\pi$, $B_{01}=(q(1-q)/M)C_M-\beta q(1-q)k_\pi$ is **not** an identity of the raw gap — its quadratic form on $\alpha$-masses misses every cross-component product entropy $h(x_iy_j)$ (worst machine-verified mismatch $0.168$ at 90 dps; points exist where the displayed form is negative while the raw gap is positive). The earlier `refuted_routes[0] "k_pi is PSD"` conflated block-level PSD of that display with kernel-level PSD; the kernel fact $k_\pi=xy+x(1-x)y(1-y)=\phi_1\otimes\phi_1+\phi_2\otimes\phi_2$ (rank-2 PSD) is TRUE and recorded, but no $k_\pi$ bilinear enters the raw gap. The exact structure is the two-piece diagonalization $gap=(1/M)\sum_{i,j\le6}w_iw_jD_M(s_i,s_j)+q(1-q)c$ (exact symbolic zero residual in the 33-symbol free log basis; $D_M=st\,C_M$; $c=2(1-\beta)\,\mathrm{cross}+\beta(\pi_0+\pi_1)-2\langle E\rangle_{\mathrm{cross}}$ is $q$-free and $c=0$ on the diagonal $P_0=P_1$). Since the first piece is exactly the certified $q=1$ gap $F(P_{\rm mix})\ge0$ of the mixture law, **the entire remaining paired-class H2 obligation is the scalar margin $F(P_{\rm mix})+q(1-q)c\ge0$** with $c$ sign-indefinite (recorded negatives $-0.0385$, $-0.0168$). | **PROVED exact diagonalization + PROVED refutation of the displayed block identity**; remaining scalar margin OPEN | [`uc/liu9_block_copositive.py`](uc/liu9_block_copositive.py) (T1/T2/T3 symbolic; 4/4 mutations fail), report [`liu9-block-copositive-diagonalization.json`](uc/verification/results/liu9-block-copositive-diagonalization.json) inner sha256 `5c856c6ec95faed2e7048783fd36d19cd91c2a2fdd300f6ff007656d68842f79` (byte-stable across two runs, file `fe5ab3d7d6bfe5d8490a501cf5039a3e6ba625e52461b18825e007f46f9bc652`), notes [`uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md`](uc/LIU9_BLOCK_COPOSITIVE_2026-08-29.md); independent Main audit at 200 dps reproduced trial-0 `gap_mp` to $2.5\times10^{-31}$ (30-digit recorded inputs account for the shift; in-run identity residual $9.2\times10^{-92}$) | `./.venv/bin/python -I -B uc/liu9_block_copositive.py`, then `./.venv/bin/python -I -B uc/verification/test_liu9_qendpoint_lift.py` for the unchanged endpoint/q stack |
| Liu H2 mean-feasible witness search | Staged discovery-only search over the mean-feasible family (structured corners, 330,000 biased random draws across standard+deep budgets, constrained SLSQP polish of the best points): no negative raw gap. Minimum accepted value is the EHX-to-zero boundary `5.5162e-15` (EHX itself `5.5e-15`, a zero-support face point); best strictly interior raw gap `7.6127e-5`. Polish converges only to the proved zero-support faces. Gaps evaluated at the declared 50-decimal mpmath context (workdps), not the ambient default. Independent checker re-derives every recorded gap at 90 decimals (worst relative agreement `9.8e-17`), re-checks exact mean feasibility (20/20), and escalates any point below $-10^{-12}$ to 200 decimals; zero escalations, zero witnesses. Digest-tamper and sign mutations verified; two-run smoke determinism verified after replacing a PYTHONHASHSEED-dependent jitter seed. | **COMPUTATIONAL EVIDENCE**, not certificate | [`uc/liu9_block_witness.py`](uc/liu9_block_witness.py) + [`uc/verification/independent_block_witness_check.py`](uc/verification/independent_block_witness_check.py); search `sha256 adf08e4e36d882cbb22129fdc31b7016a7c60fe8c2488d2834a412450eea43db`, check `sha256 69e074b757dca3889ac9f9231349a06137a794bda37bf2d3f9148524c4abb989` | `./.venv/bin/python -I -B uc/liu9_block_witness.py --stage deep && ./.venv/bin/python -I -B uc/verification/independent_block_witness_check.py` |
| Liu H2 q=1 feasible search | Feasible starts in all 230 active boxes; 4,965 feasible raw-gap evaluations, 4,465 invalid points rejected before interpretation, minimum 0 at binary boundary, no feasible negative. Known ambient negative reproduced and rejected by mean deficit. | **COMPUTATIONAL EVIDENCE**, not certificate | [`liu9-qone-feasible-search.json`](uc/verification/results/liu9-qone-feasible-search.json) `sha256 690d0b48c7d7df8d773815e5a7198b37d7545f37223bdd381e05c022a4fb8690`; independent verifier | — |
| Liu H2 protocol identity, attribution | Liu's eq. (87) $x^{*2}+x^{*2}(1+\bar x^{*2})=1$ **is** $\mathrm{prot}(x,x)=1-x^2$, stated directly as a defining relation; expanding gives $x^4-2x^3+3x^2-1=0$, which **never appears in the paper**. Unremarked there and in all ten citing works: entropy symmetry forces $h(\mathrm{prot}(x,x))=h(x^2)$, so $K(x,x)=h(x^2)$ for **every** $\beta$ and $pK(x,x)=h(x)$. Root uncatalogued in OEIS. | **PROVED** (exact integer polynomials); novelty is **absence after stated search** | primary source `LIU_H1/literature/pdfs/liu_2023_arxiv_2306.08824v1.pdf`, verified in-repo | `pdftotext` + [`uc/liu9_transverse.py`](uc/liu9_transverse.py) |

The load-bearing UC chain is the machine-verified finite certificate embedded
in the human-audited A/A′, B‴, Margin, endpoint, and entropy-bridge arguments.
No single command verifies the entire chain.

## Campaign ladder

Campaign identities, budgets, completion counts, verdicts, and stop reasons
come from [`uc/campaigns/README.md`](uc/campaigns/README.md). Totals below mean
the exact sum of `tallies.processed` in the committed `result_slice*.json`
records linked by each campaign directory. These are **PROVED record facts**,
not a mathematical certificate unless the verdict says **CERTIFIED**.

| Campaign | Budget per slice | Committed slices | Total processed boxes | Verdict | Why it ended there |
|---|---:|---|---:|---|---|
| [E `…5bde6bea1eca`](uc/campaigns/cert3_20260815T213603Z_db655540567b41aa9fc545c59b11520d_5bde6bea1eca/) | 8 h | 6 COMPLETE / 8 records | 158,700,894 | **NO CERTIFICATE** | Slices 6 and 7 reached the wall with pending stacks 12 and 18; residual was zero. |
| [F `…eb00229a822b`](uc/campaigns/cert3_20260816T053855Z_9267788923da4b79bfa76b97c904a4e3_eb00229a822b/) | 12 h | 7 COMPLETE / 8 records | 196,113,302 | **NO CERTIFICATE** | Slice 7 exposed 1,447,470 `sink` residuals because `ratio_plausible` hid the Margin-Lemma ratio rule. |
| [G `…4c2b69644e5f`](uc/campaigns/cert3_20260816T175635Z_cabbeb698e6e4e7780ba8f3ea96346eb_4c2b69644e5f/) | 24 h | 7 COMPLETE / 8 records | 187,614,798 | **NO CERTIFICATE** | Removing that gate eliminated the residual barrier; slice 7 then stopped only at the wall with stack 16 and residual zero. |
| [H `…2f23a58ebdb8`](uc/campaigns/cert3_20260817T212001Z_4e6268eb7ce44aa3883ef15ea401c7c7_2f23a58ebdb8/) | 24 h | 5 COMPLETE / 8 records (0–3, 7) | 366,265,482 | **NO CERTIFICATE** | The influence-weighted splitter completed the first slice 7, but slices 4–6 hit the wall with stacks 25/27/30 and residual zero. |
| [I `…2f23a58ebdb8`](uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/) | 72 h | 8 COMPLETE / 8 records | 488,465,854 | **CERTIFIED 2026-08-21** | All stacks, residuals, and budget flags were zero; the frozen collector replayed all eight traces and printed `COMPOSITE CERTIFICATE`. |
| [J `…57a8908ba7b8`](uc/campaigns/cert3_20260820T181524Z_08ad738a1cf24a08bfc4189e65f9ef44_57a8908ba7b8/) | 72 h | 6 COMPLETE / 8 records (0–4, 7) | 584,645,062 | **NO CERTIFICATE** | Slices 5 and 6 hit the 72 h wall after 142,886,584 / 156,713,996 boxes with residual 0 but stack=`budget_time` 31 / 29. The frozen collector rejected at admission with exactly six problems and did not replay. |


To recompute any row, replace `<campaign-id>` by the linked directory name:

```sh
./.venv/bin/python -c 'import glob,json; r=[json.load(open(p)) for p in glob.glob("uc/campaigns/<campaign-id>/result_slice*.json")]; print(len(r), sum(x["verdict"]=="COMPLETE" for x in r), sum(x["tallies"]["processed"] for x in r))'
```

Campaign J targeted rational
$t=0.3821660112501052\ge\psi+2\times10^{-4}$; the exact relation is checked in
its [`launch.json`](uc/campaigns/cert3_20260820T181524Z_08ad738a1cf24a08bfc4189e65f9ef44_57a8908ba7b8/launch.json).
Its six COMPLETE counts are 21,891,091 / 21,637,061 / 24,563,611 /
37,877,105 / 85,793,971 / 93,281,643. Slices 5 and 6 were time-limited only;
all 584,645,062 processed boxes across the campaign had residual zero.

## Reproduction

Do **not** use a live source checkout or the root `requirements-freeze.txt` as
the certificate pin. Campaign I's own `launch.json` pins CPython 3.14.3,
python-flint 0.9.0 (Arb), mpmath 1.3.0, numpy 2.5.2, scipy 1.18.0, sympy
1.14.0, and macOS 26.5.2 arm64, together with every executable/review source
SHA-256. Use the existing interpreter
`/Users/jinleic/jinleic-workspace/math/.venv/bin/python`.

The exact one-line collector invocation is:

```sh
cd /Users/jinleic/jinleic-workspace/math && ./.venv/bin/python -B uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/snapshot/cert3_collect.py uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/launch.json
```

Allow about **25 h single-core wall time**; this expectation and the command are
recorded in
[`data/uc-certificate-psi-plus-1e-4-2026-08-21/BUNDLE.md`](../data/uc-certificate-psi-plus-1e-4-2026-08-21/BUNDLE.md).
Success means exit 0 after all eight `REPLAY slice=i PASS` lines and the exact
`COMPOSITE CERTIFICATE` statement above. Every other outcome is rejection.
The immutable evidence bundles are:

- `/Users/jinleic/jinleic-workspace/data/uc-certificate-psi-plus-1e-4-2026-08-21/`
- `/Users/jinleic/jinleic-workspace/data/uc-certificate-psi-plus-1e-4-2026-08-21.zip`
- `/Users/jinleic/jinleic-workspace/data/uc-certificate-psi-plus-1e-4-2026-08-21-v2.zip`

## Open problems and adopted next targets

These rows come from the “Adopted — next target after H” section of
[`PROGRESS.md`](PROGRESS.md) and the as-executed
[`uc/campaign_i_plan.md`](uc/campaign_i_plan.md). Projections remain explicitly
**NUMERICAL**, even when later campaign records provide actual counts.

| Target | Current status | Evidence label and boundary |
|---|---|---|
| Close Campaign J at rational $0.3821660112501052\ge\psi+2\times10^{-4}$. | **CLOSED — NO CERTIFICATE:** six slices COMPLETE; slices 5 and 6 hit the 72 h wall with stack 31/29, residual zero. The collector rejected exactly those two verdict/stack/time triplets. | **PROVED record facts**; verbatim negative-control output in J's `collect_output.txt`. |
| Certify rational $0.3822660112501052\ge\psi+3\times10^{-4}$ (campaign K). | **CERTIFIED 2026-08-31 — `COMPOSITE CERTIFICATE`, exit 0.** All 16/16 slices COMPLETE with residual=stack=`budget_time`=0 (committed `processed` = 805,542,368 boxes). Slices 10, 12, 13 — lost with the daemon broker on 2026-08-22 — were relaunched 2026-08-27 under the same launch manifest (which cannot overwrite a committed artifact) and finished 2026-08-28. The frozen collector then replayed **all sixteen traces, 16/16 PASS, 805,542,368 boxes re-derived from their exact dyadic roots with every claimed discharge re-proved in Arb**, 117,615 s single-core (slice 12 alone: 166,211,875 boxes in 25,328 s), and printed the certificate for $\Phi\ge0$ on the feasible 5-parameter family over $w\in[1/2,1]$, extended to $w\in[0,1]$ by the proved orbit-swap symmetry. Verdict transcript [`collect_output.txt`](uc/campaigns/cert3_20260822T015019Z_41f2e119129446d3b1a936d66bbe7683_b49a21ec80ba/collect_output.txt) (sha256 `e4e799a8261b8ffde7d8993d570eb7399816e7873433dd26bf2ebb732e454602`); the pre-completion rejection is preserved verbatim as `collect_output_2026-08-24_rejection.txt` and at [`collect_rejection_precompletions.txt`](/Users/jinleic/jinleic-workspace/data/uc-certificate-psi-plus-3e-4-2026-08-29/collect_rejection_precompletions.txt). Recorded process fact: the first collector run died after ~33 h without persisting its verdict (terminal-only stdout), so the checker was re-run under supervision with the transcript captured to disk. | **MACHINE-VERIFIED finite certificate.** Checker `6a2ae5ea06…`, campaign `20260822T015019Z_41f2e119…`, code `b49a21ec80ba…`, launch `2170d56f651c…`. The claim is the explicit-constant frequency statement at the stated rational $t$ — not Frankl's $1/2$, and independent of Liu's Hypothesis 2. |
| Go to or beyond the two-orbit obstruction $c^*$. | Earlier closure-blind routes and the proposed $1/50$ scalar repair are **REFUTED** by the certified $n=6/n=7 families. Gate B is now **RESOLVED**: Cartesian powers of the certified normalized $n=7$ base remain cap/Reimer and satisfy $A_+(\mathcal F_k)=kA_+(\mathcal B)$, $\varepsilon_\vee(\mathcal F_k)=1-(17/81)^k$, and ratio $>7k/250\to\infty$. Thus no fixed scalar closure-defect coefficient works. The larger UC-specific programme remains open because the base is non-UC; a separate small-defect question with $\varepsilon_\vee\to0$ is also open. | **PROVED/CERTIFIED:** exact proof and evidence in `uc/gate_b/`. **COMPLETE EXACT + NUMERICAL:** declared $n=6/n=7$ classes and the new $n=8$ $S_1\times S_7$ class. No union-closed counterexample or improvement beyond the certified $0.3820660112501052$ frequency constant is claimed. |
| Prove Liu's Hypothesis 2 (Section V-B global minimizer), on which the conditional constant $0.3827090879\ldots$ depends. | **OPEN only on an endpoint layer of width $2^{-14}$ (2026-09-01).** The continuum is eliminated analytically, and $R\ge0$ and $A\ge0$ are PROVED on all of $[0,1]^2$ while $B\ge0$ is PROVED on $99.9756\%$ of it plus the whole diagonal and both axes. H2 for the paired class reduces **exactly** from 9 variables $+\,q$ $+$ the mean constraint to **two two-variable inequalities on $[0,1]^2$**, valid for arbitrarily many atoms. With $h$ the natural-log binary entropy, $\pi(s,t)=st(1+(1-s)(1-t))$, $P_2=(1-\beta)h(st)-\frac{t h(s)+s h(t)}{2m}$, $Q_2=\beta h(\pi)\ge0$, $R=P_2+Q_2$, and $\varphi(s)=\sqrt{Q_2(s,s)}$, set $A=P_2+\varphi\otimes\varphi$ and $B=Q_2-\varphi\otimes\varphi$, so $R=A+B$. Then $T=\langle\mu,B\mu\rangle+\langle\nu,B\nu\rangle+2\langle\mu,A\nu\rangle+(\langle\mu,\varphi\rangle-\langle\nu,\varphi\rangle)^2$ and $\mathrm{gap}=\frac{(1-q)^2I_{00}+q^2I_{11}}{M}+q(1-q)T$ for **every** $q$, whence $\boxed{A\ge0\ \wedge\ B\ge0\Rightarrow\text{H2}}$. This is Diananda's $N+\mathrm{PSD}$ with a rank-one PSD part, made explicit; $\varphi$ is **forced** (scaling it up breaks $B$ on the diagonal, down breaks $A$ at $x^*$). **NO LONGER REQUIRED:** the activation lemma $\min\{M\ge m\}=\min\{M=m\}$ and its off-grid gap, the active-face $q^*$ elimination, the $\Delta=0$ case split, every box cover of $[0,1]^7$, $3\times3$ copositivity, 6-dimensional branch-and-bound. **Why covers were doomed:** $R$ has a *nondegenerate interior zero* at $(x^*,x^*)$, $x^*=0.6907875939249880141505$ — $R$ and $\nabla R$ both enclose $0$, Hessian eigenvalues $0.2568180929472731604836$ and $1.249573351$ — so no interval subdivision can ever certify it. The channel is second order there, $c\sim\beta K_{12}\Delta^2$ with $K_{12}=-1.59260964939$, and **cancels exactly** in the critical symmetric direction, leaving margin $H_{11}+H_{12}=0.2568180929$ independent of $\beta\lvert K_{12}\rvert=0.1593446723$. **$R\ge0$ on $[0,1]^2$ is PROVED** — 400-bit Arb, exhaustive three-stratum cover (near-zero strip $\min(s,t)\le10^{-8}$ with $R\ge st\cdot0.03708$; local $x^*$ square of half-width $9/500$ inside the certified $\rho=3\lambda/C_3=0.026969$, value and both gradients enclosing $0$, $\lambda_{\min}\ge0.2568180929472731604836$, $C_3\le28.5676185$; $(1,1)$ corner $R\ge u(0.189954\log(1/u)-\frac{\log2+1}{2m})$), plus $14$ rational rectangles cleared by branch-and-bound — $7{,}462$ processed $=3{,}738$ accepted $+3{,}724$ split, **no unresolved cell**, worst bound $+1.0135\times10^{-18}$; exhaustiveness machine-asserted. I re-ran it and recomputed its hash, and its constants are conservative in the right direction (my $\max\lvert\partial^3R\rvert=5.8785$ below their $10.1002$). Hence $I_{00},I_{11}\ge0$ unconditionally and **H2 narrows to the single sandwich $0\le B\le R$** (the right half is $A\ge0$). **NEW EXACT STRUCTURE:** $x^*$ is ALGEBRAIC of degree 4 — the root in $(0,1)$ of $x^4-2x^3+3x^2-1=0$, which is exactly $\pi(x,x)=1-x^2$, so $h(\pi)=h(x^2)$ by $h(u)=h(1-u)$ and $D(x^*,x^*)=m\,h(x^{*2})-x^*h(x^*)=0$ identically; and $m=x^*h(x^*)/h(x^{*2})$ in closed form (verified to $10^{-80}$). The tight point is therefore exact algebra, not a solved decimal. **REFUTED with exact witnesses:** entrywise nonnegativity ($\Psi(0,\tfrac12,\tfrac12,0)=-0.0553708106423719494349154$, since the cross-protocol terms cancel identically giving $\Psi=P_2(s,\tau)+P_2(t,\sigma)+Q_2(s,t)+Q_2(\sigma,\tau)$), the shortcut $L_{st}\le0$ for $B\ge0$ ($17{,}906/249{,}001$ points violate it, worst $+0.0784896$ at $(0.17,0.998)$). The older framings — the displayed $2\times2$ block kernel, `liu9-block-kernel.json`'s `remaining_obligation`, `\Psi\ge0$ on $[0,1]^4$, and $M_T$ copositive on $[0,1]^6$ — are all **SUPERSEDED** and must not be cited as the live obligation. | **MACHINE-VERIFIED reduction** (Arb 480 bits, six chain claims, residuals $\le3.2\times10^{-68}$, 4 mutations caught, two runs byte-identical, `23cb0c56\ldots`; independently re-derived symbolically with exact zero residuals in `liu9-psi-reduction.json`, 4/4 mutations). **PROVED (400-bit Arb, `liu9-h2-twovar.json`, byte-identical, 3 mutations/lemma, exhaustive machine-asserted partition, no unresolved cell):** $R\ge0$ and $A\ge0$ on all of $[0,1]^2$; the $x^*$ centre is a certified Arb isolating ball for the quartic root (220 bisections, $f(lo)=-5.27\times10^{-70}$, $f(hi)=+1.02\times10^{-69}$, $f'>0$), with $m=x^*h(x^*)/h(x^{*2})$ cross-checked against the binding solver ($\pm2.56\times10^{-69}$, contains 0). **PARTIAL (`liu9-h2-phi.json`):** $B\ge0$ on $[a,1-a]^2$, $a=2^{-14}$, area $67092481/67108864$, zero unresolved cells across 797,282 cells in three strata. **OPEN:** the residual layer, exact area $16383/67108864$, where the only available bound is the vacuous $\Phi\ge-(\log2)^2$; it needs an endpoint asymptotic remainder lemma, since $h'$ is singular at entropy arguments $0,1$. Matched asymptotics are derived (the both-small regime is a perfect square $(2st\log(s/t))^2$) but **not yet certified**. |
