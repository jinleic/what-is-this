# Union-closed Campaign I proof and certificate audit

> **Integration update (2026-08-26).** After this line-by-line audit, the
> manuscript now states the exact-to-relaxed functional bridge, handles
> optimal-coupling attainment and \(L=0\), proves the entropy-to-UC bridge
> directly, and separates human from machine evidence.  A fresh direct replay
> then processed all eight traces (488,465,854 nodes, zero residuals), with
> every exit code, stdout hash, and tally matching; independent structural
> reports before and after have identical input maps.  Evidence:
> [`verification/results/direct-isolated-replay/summary.json`](verification/results/direct-isolated-replay/summary.json),
> canonical SHA-256
> `6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`.
> The verifier was additionally hardened after independent review: external
> lock bytes are hash-parsed from one descriptor, dependency module hashes are
> pinned, unauthenticated resume is disabled, PASS schemas constrain execution
> fields, and an executable eight-slice aggregator validates composite reports.
> A second full replay from fresh read-only staging also passed all eight
> slices under those hardened gates; the composite report SHA-256 is
> `57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.
> Discrepancy line references below preserve the pre-correction snapshot that
> exposed each issue.

## 0. Scope, evidence convention, and verdict

This audit reconstructs the implication chain for the exact rationals

\[
\alpha_0=\frac{356069}{10^7}=0.0356069,
\qquad
t_{\rm cert}=\frac{955165028125263}{2500000000000000}
=0.3820660112501052.
\]

The immutable campaign examined here is

`uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8/`.

In citations below, `I/` means that exact directory and `S/` means its
`snapshot/` subdirectory. Line ranges always refer to the frozen Campaign I
copy, not a later development copy. The principal prose statements compared
against it are `uc/PROOF.md:1-455` and `uc/paper/main.tex:213-1110,1552-1646`.

Status words have the following deliberately narrow meanings.

- **MACHINE-VERIFIED**: a named exact-symbolic, interval-arithmetic, or replay
  gate passed. This status covers only what that gate actually computes; it
  does not authenticate a citation or prove the mathematical interpretation
  of hand-written verifier code.
- **HUMAN-AUDITED**: the displayed proof has been checked line by line here,
  including its domains, signs, and boundary cases, but is not a formal proof
  object.
- **COMPUTATIONAL-EVIDENCE**: a numerical sample, optimization seed, result
  record, timing, or archived transcript supports a claim but is not by itself
  a universal proof.
- **CITED-DEPENDENCY**: a named external theorem is used rather than proved in
  the frozen sources.
- **OPEN**: not established by this chain or still requiring a submission
  correction.
- **FAILED**: the stated wording or check does not establish what it says. A
  failed wording claim can coexist with a sound mathematical resolution.

### Audit verdict

**HUMAN-AUDITED:** after supplying three omitted but sound lemmas—existence of
an optimal self-coupling, the entropy-zero case, and the sequential
union-closed bridge—I found no circular inference, uncovered feasible region,
reversed load-bearing inequality, or unsound reuse of a one-sided interval
bound in the mathematical chain. The exact five-parameter functional actually
certified is a lower relaxation of the exact pair-orbit functional, in the
correct direction.

**MACHINE-VERIFIED, narrow scope:** the archived frozen collector transcript
reports successful replay of all eight Campaign I traces and
`488,465,854` processed nodes, and a targeted rerun during this audit passed
`S/lemma_rh_proof.py` (44-cell Arb cover) and `S/bridge_uc.py`. The full
488-million-node replay was not repeated as part of this document; the
archived replay occurrence is evidenced by `I/collect_output.txt:1-25` and is
independently reproducible from the frozen artifacts.

**Not an end-to-end machine proof:** the decomposition, support reduction,
Bochner argument, verifier-to-inequality correspondence, endpoint derivative
estimate, Shannon bridge, and imported functional-analysis results contain
human-audited or cited obligations. Moreover, replay calls the same discharge
functions as the worker, so it is independent reconstruction of the DFS proof
object, not an independent implementation of the interval inequalities
(`S/cert3_replay.py:22-29,187-331`). No overall union-closed claim in this
audit is labeled end-to-end machine-checked.

## 1. Exact statement and the two functionals that must not be conflated

Let \(h(u)=-u\log_2u-(1-u)\log_2(1-u)\), continuously extended by
\(h(0)=h(1)=0\). For a regular Borel probability measure
\(\mu\in P([0,1])\), define

\[
L(\mu)=\int h(p)\,d\mu(p),\qquad
Q(\mu)=\iint h(p+q-pq)\,d\mu(p)d\mu(q),
\]

and

\[
s^*(p,r)=\max\!\left(p,r,\min(p+r,\tfrac12)\right).
\]

For the set \(\Pi(\mu,\mu)\) of all self-couplings, put

\[
C(\mu)=\min_{\pi\in\Pi(\mu,\mu)}\int h(s^*(p,r))\,d\pi(p,r),
\qquad
F(\mu)=(1-\alpha)Q(\mu)+\alpha C(\mu)-L(\mu).
\]

These are the theorem quantities in `uc/paper/main.tex:240-258`. The minimum
exists; Section 4.2 below supplies the compactness argument omitted at the
point where the paper writes `min`.

For a pair-orbit law \(\nu\in P(\Delta)\),
\(\Delta=\{0\le p\le q\le1\}\), let

\[
\mu_\nu=\int\frac{\delta_p+\delta_q}{2}\,d\nu(p,q),
\quad
M_\nu=\int\frac{p+q}{2}\,d\nu,
\quad
C_\nu=\int h(s^*(p,q))\,d\nu.
\]

The exact pair-orbit functional is

\[
\Phi_{\rm ex}(\nu)
=(1-\alpha)Q(\mu_\nu)+\alpha C_\nu-L(\mu_\nu).
\]

The decomposition defines a positive-semidefinite energy \(E\) and

\[
\operatorname{rh}(p)=2(1-p)h(p)-h((1-p)^2),
\quad
S_\nu=\int\frac{\sqrt{\operatorname{rh}(p)}+
\sqrt{\operatorname{rh}(q)}}2\,d\nu.
\]

The **different** functional evaluated by `cert3.py` is

\[
\boxed{
\Phi_{\rm rel}(\nu)=
[2(1-\alpha)(1-M_\nu)-1]L(\mu_\nu)
+\alpha C_\nu-(1-\alpha)S_\nu^2.}
\]

The load-bearing comparison is

\[
\Phi_{\rm ex}(\nu)-\Phi_{\rm rel}(\nu)
=(1-\alpha)(S_\nu^2-E(\mu_\nu))\ge0. \tag{1}
\]

This is stated correctly in `uc/PROOF.md:160-169` and follows from the Margin
Lemma feature-map triangle inequality (`S/margin_lemma.py:9-37`). It is not
stated explicitly in the paper between the exact \(\Phi_\alpha\) of
`uc/paper/main.tex:615-768` and the certificate rules at
`uc/paper/main.tex:934-1045`. That omission is the most important
verifier-to-theorem presentation defect found in this audit: code variables
called `C` are \(C_\nu\), not the optimal-coupling value \(C(\mu_\nu)\).

The exact logical route is therefore

\[
\Phi_{\rm rel}\ge0\text{ on all feasible two-orbit laws}
\Rightarrow
\Phi_{\rm ex}\ge0\text{ at an exact global minimizer}
\Rightarrow
\inf_{\mathbb E_\mu p\le t}F(\mu)\ge0. \tag{2}
\]

It is not the generally false pointwise assertion
\(F(\mu_\nu)\ge\Phi_{\rm rel}(\nu)\): \(C_\nu\) may exceed
\(C(\mu_\nu)\). Equality of the two global optimization problems, plus the
support theorem, is essential.

**Status:** the comparison (1) is **HUMAN-AUDITED**; the interval evaluation
of \(\Phi_{\rm rel}\) at every accepted trace leaf is **MACHINE-VERIFIED**. Evidence:
`S/reduction.py:15-51`, `S/thmB3_proof.py:8-20,129-207,271-280`,
`S/margin_lemma.py:9-37`, `uc/PROOF.md:160-169`, and the code formula
`S/cert3.py:1065-1093`.

## 2. Dependency DAG

| ID | Component and exact conclusion | Depends on | Status | Principal evidence |
|---|---|---|---|---|
| D0 | Entropy, \(s^*\), regular Borel measures, and exact \(\alpha,t\) | — | **HUMAN-AUDITED** definitions; exact target relation **MACHINE-VERIFIED** | `uc/paper/main.tex:213-258`; `S/arbcore.py:8-16,29-74`; `S/cert3_par.py:198-224` |
| D1 | Cambie formula equals repository clipped median | D0 | **HUMAN-AUDITED** universal three-case proof; finite schema **MACHINE-VERIFIED** | `S/bridge_uc.py:7-20,68-109`; targeted `-B S/bridge_uc.py` run |
| D2 | Every self-coupling symmetrizes and folds to a pair-orbit law, with equal marginal and cost | D0, D1 | **HUMAN-AUDITED**; finite coefficient schemas **MACHINE-VERIFIED** | `S/bridge_uc.py:22-38,133-166,169-245`; `uc/paper/main.tex:597-613` |
| D3 | \(\inf F=\inf\Phi_{\rm ex}\), with both minima attained | D2, D6 | **HUMAN-AUDITED** after adding compact optimal-coupling lemma | Section 4.2; `S/thmB3_proof.py:26-58,138-146` |
| D4 | Kernel decomposition \(H(xy)=g\otimes g-a\otimes a-T(xy)\) | D0 | **HUMAN-AUDITED** exact algebra; an exact SymPy assertion exists but is not collector-run | `S/decomposition.py:17-42,220-248`; `uc/paper/main.tex:304-367` |
| D5 | \(Q=2BL-E\), \(E\ge0\), moment-square expansion | D4 | **HUMAN-AUDITED** with uniform-convergence justification | `S/reduction.py:15-45`; `uc/paper/main.tex:387-450` |
| D6 | Compactness, continuity, and global attainment on \(P(\Delta)\) | D0, D5 | **CITED-DEPENDENCY** for Riesz/Banach–Alaoglu/Stone–Weierstrass; hypotheses **HUMAN-AUDITED** | `S/thmB3_proof.py:26-146,283-303` |
| D7 | Fixed-\(B\) concavity | D5, D6, \(0\le\alpha\le1\) | **HUMAN-AUDITED** exact quadratic identity; sampled controls only **COMPUTATIONAL-EVIDENCE** | `S/thmB3_proof.py:149-193,430-531` |
| D8 | Bauer minimizer exists at an extreme point | D6, D7 | **CITED-DEPENDENCY**, hypotheses **HUMAN-AUDITED** | `S/thmB3_proof.py:195-207`; Bauer DOI `10.1007/BF01898615`; Bru–Pedra arXiv:1610.03411, Lemma 3.3 |
| D9 | Extreme one-moment slice has support at most two | D8 | **HUMAN-AUDITED** direct perturbation; Winkler/Pinelis is corroborative, not load-bearing | `S/thmB3_proof.py:209-265`; `uc/paper/main.tex:710-728` |
| D10 | Exact global minimizer is a two-orbit law | D3, D6–D9 | **HUMAN-AUDITED** | `S/thmB3_proof.py:267-280` |
| D11 | \(E\le S^2\), hence \(\Phi_{\rm ex}\ge\Phi_{\rm rel}\); handle \(L=0\) separately | D5 | **HUMAN-AUDITED** | `S/margin_lemma.py:15-37`; Sections 4.5–4.6 |
| D12 | \(\operatorname{rh}\) collars, \(\rho\) caps, and endpoint derivative cap | D5 | Collar theorem **MACHINE-VERIFIED** by targeted frozen run plus **HUMAN-AUDITED** signs; derivative cap **HUMAN-AUDITED** | `S/lemma_rh_proof.py:42-207,266-424`; `S/cert2.py:51-65,193-313`; `S/cert3.py:126-171` |
| D13 | Every code 0–8 discharge is a valid infeasibility proof or lower bound for \(\Phi_{\rm rel}\); splits preserve coverage | D0, D11, D12 | Rule derivations **HUMAN-AUDITED**; accepted leaf arithmetic reported **MACHINE-VERIFIED** | Section 6; `S/cert3.py:1406-1593` |
| D14 | Trace reconstruction covers all eight roots and rejects residual/pending work | D13 | **MACHINE-VERIFIED** by frozen, structural, and secure independent byte-zero replays | `S/cert3_replay.py:187-331`; `verification/results/independent-arithmetic/secure-full/composite.json` |
| D15 | Campaign I establishes \(\Phi_{\rm rel}\ge0\) on the feasible five-parameter domain | D13, D14 | **MACHINE-VERIFIED** independently at the Python formula/control layer, conditional on the stated trace/Arb/runtime/custody trust | secure composite and report lock; archived Campaign I records |
| D16 | \(F(\mu)\ge0\) for all \(\mathbb E p\le t_{\rm cert}\) | D2–D15 | **HUMAN-AUDITED** implication from machine arithmetic and analytic lemmas | Equations (1)–(2); `uc/PROOF.md:95-169` |
| D17 | Finite union-closed consequence | D16 | **HUMAN-AUDITED** self-contained reconstruction here; current paper labels it **CITED-DEPENDENCY** | Section 4.1; Cambie arXiv:2212.12500v2, Q2 and Section 4 |

There is no edge from the certificate back into D4–D12. In particular, the
support theorem uses the exact energy \(E\), not the relaxed certificate, and
the \(\rho\) endpoint theorem is proved independently of any branch-and-bound
clear. No circularity was found.

## 3. Quantifier, domain, continuity, and sign ledger

| Object or assertion | Exact quantifier/domain | Required assumption and boundary convention | Status and evidence |
|---|---|---|---|
| \(h\) | Every \(u\in[0,1]\) | Base 2; \(0\log0=0\); continuous and bounded | **HUMAN-AUDITED**; `S/arbcore.py:29-43`, `S/thmB3_proof.py:64-88` |
| Analytic \(\alpha,t\) theorem | Every fixed \(\alpha,t\in[0,1]\) | \(1-\alpha\ge0\) is needed for concavity and (1) | **HUMAN-AUDITED**; `S/thmB3_proof.py:8-20,175-193` |
| Certified parameters | One exact \(\alpha_0\), one exact \(t_{\rm cert}\) | Arb balls enclose these rationals; KKT decimals are merely nonnegative shift choices | **MACHINE-VERIFIED** arithmetic / **HUMAN-AUDITED** interpretation; `I/launch.json:1`, `S/arbcore.py:12`, `S/cert3.py:489-493` |
| \(\mu\) | Every regular Borel probability on \([0,1]\) with \(\int p\,d\mu\le t\) | Compact domain makes all continuous entropy kernels integrable | **HUMAN-AUDITED**; `uc/paper/main.tex:213-225,240-258` |
| Self-coupling \(\pi\) | Every regular Borel probability on \([0,1]^2\) with both marginals \(\mu\) | Nonempty via \(\mu\otimes\mu\); weak-* compact; symmetric cost | **HUMAN-AUDITED**; Section 4.2 |
| \(\nu\) | Every regular Borel probability on compact \(\Delta\) | Weak-* topology; marginal map is continuous linear | **CITED-DEPENDENCY** plus **HUMAN-AUDITED** mapping; `S/thmB3_proof.py:26-58,108-132` |
| Two-orbit parameters | \(0\le p_i\le q_i\le1\), \(0\le w\le1\) | Coincident points and weights 0 or 1 are allowed, so one-orbit laws are included | **HUMAN-AUDITED**; `uc/paper/main.tex:572-595` |
| Feature map | Every \(p\in[0,1]\), into \(\ell^2\) | Uniformly summable tail gives a bounded strongly measurable map; hence Bochner integrable | **HUMAN-AUDITED**; `S/reduction.py:19-31`; `S/margin_lemma.py:15-37` |
| Margin quotient | \(L(\mu)>0\) only | \(L=0\) must not be divided by; handled separately in Section 4.6 | **HUMAN-AUDITED**; `S/margin_lemma.py:25-37` |
| \(\rho=\operatorname{rh}/h\) | Only \(0<p<1\) | Endpoint collars avoid division by an interval containing zero entropy | **HUMAN-AUDITED** and targeted theorem **MACHINE-VERIFIED**; `S/cert2.py:193-282` |
| Shift \(\lambda\) | Each manifest value, all \(\lambda\ge0\) | Correct direction is \(\Phi_\lambda=\Phi+\lambda(M-t)\le\Phi\) on \(M\le t\) | **HUMAN-AUDITED**; `S/bound_kkt.py:139-177`; `S/cert3.py:1096-1186` |
| Corner coefficient | Only clear when \(\kappa_*=2(1-\alpha)(1-m_*)-1>0\) | Positivity licenses multiplying a lower entropy bound | **HUMAN-AUDITED**; `S/diag_exhaust.py:163-207` |
| Pin coefficient | Only pin when its strip enclosure has `coeff_lo > 0` | Negative entropy-derivative term then has the stated sign | **HUMAN-AUDITED**; `S/cert3.py:377-427` |
| Union-closed family | Every finite \(\mathcal F\ne\{\varnothing\}\) | Choose a present element first; absent ground-set elements may be discarded | **HUMAN-AUDITED**; Section 4.1; primary Cambie Section 4 |

The exact target comparison is particularly close and therefore must not be
replaced by a binary-float comparison. With

\[
r=t_{\rm cert}-10^{-4}
=\frac{954915028125263}{2500000000000000},
\]

exact rational arithmetic gives

\[
r^2-3r+1=
-\frac{673679581180831}{6250000000000000000000000000000}<0.
\]

Since \(0\le r\le3/2\), this proves \(r\ge(3-\sqrt5)/2\). This is exactly the
predicate in `S/cert3_par.py:198-210` and `S/cert3_collect.py:128-136`.
**Status: MACHINE-VERIFIED** exact arithmetic, with no float comparison in the
acceptance gate.

## 4. Detailed analytic audit

### 4.1 The union-closed bridge, including the primary-source wording defect

Cambie's v2 Question 2 says **expectation less than \(c\)**, not “at most
\(c\)”; see the fixed primary source
[arXiv:2212.12500v2, Question 2](https://arxiv.org/html/2212.12500v2#Thmtheorem2).
Cambie's Section 4 then begins its contradiction by saying every element occurs
in “at most” a \(c\)-fraction. Thus the primary source itself switches from
`< c` to `<= c`. The local claim that Q2's hypothesis itself says “at most” is
**FAILED wording** in `S/bridge_uc.py:40-46`, `uc/PROOF.md:25-32`, and the
corresponding paper discussion. The closed-domain certificate and the strict
first-coordinate argument below repair the mathematics.
Section 4 also says that \(Q_i\) is independent of the other two parameter
variables and then calls \(P_i,Q_i\) dependent. Those clauses are inconsistent;
the construction and Q2 require \(P_i,R_i\) to be the dependent pair and
\(Q_i\) to be independent. The reconstruction below uses the mathematically
consistent reading.


A self-contained bridge is as follows.

1. Let \(A\) be uniform on a finite union-closed family \(\mathcal F\), and
   order a present element first. For every positive-probability prefix \(a\),
   define
   \[
   p_i(a)=\Pr(i\in A\mid A_{<i}=a),\qquad P_i=p_i(A_{<i}).
   \]
   Zero-probability prefixes may be assigned arbitrarily. The tower property
   gives \(\mathbb E P_i=\Pr(i\in A)\), the frequency of element \(i\).

2. Let \(B\) be an independent uniform copy. Then
   \(Q_i=p_i(B_{<i})\) is iid with \(P_i\), and conditional on both prefixes
   the OR bit has parameter \(P_i+Q_i-P_iQ_i\).

3. Construct a dependent uniform copy \(C\) sequentially. Given Bernoulli
   transition parameters \(p,r\), the possible OR probabilities form
   \([\max(p,r),\min(p+r,1)]\). Projecting \(1/2\) onto this interval gives
   \[
   s^*=\min(\max(\tfrac12,p,r),\min(p+r,1)).
   \]
   Use the four conditional masses
   \[
   \Pr(1,1)=p+r-s^*,\quad
   \Pr(1,0)=s^*-r,\quad
   \Pr(0,1)=s^*-p,\quad
   \Pr(0,0)=1-s^*.
   \]
   They are nonnegative exactly because \(s^*\) lies in that feasible
   interval, and their marginals are Bernoulli \(p\) and Bernoulli \(r\).
   Induction over prefixes (equivalently the product rule) therefore keeps
   both \(A\) and \(C\) uniform on \(\mathcal F\). The induced
   \((P_i,R_i)\), where \(R_i=p_i(C_{<i})\), is a self-coupling of the law of
   \(P_i\). This reconstructs the max-entropy Bernoulli coupling in
   [Cambie, Section 4](https://arxiv.org/html/2212.12500v2#S4).

4. Chain rule and conditioning on the finer pair of input prefixes give
   \[
   H(A\cup B)\ge\sum_i\mathbb E h(P_i+Q_i-P_iQ_i),
   \]
   \[
   H(A\cup C)\ge\sum_i\mathbb E h(s^*(P_i,R_i)),
   \qquad
   H(A)=\sum_i\mathbb E h(P_i).
   \]
   The inequality direction is correct because conditioning on
   \((A_{<i},B_{<i})\) or \((A_{<i},C_{<i})\) gives at least as much
   information as conditioning only on the corresponding OR prefix.

5. If every element frequency is at most \(t_{\rm cert}\), D16 applies to
   every coordinate law and every induced parameter coupling. Summing yields
   \[
   (1-\alpha)H(A\cup B)+\alpha H(A\cup C)\ge H(A).
   \]
   Union closure makes both unions \(\mathcal F\)-valued, so each entropy is
   at most \(\log_2|\mathcal F|=H(A)\).

6. The inequality is strict at coordinate 1. Its law is \(\delta_u\), where
   the chosen element has \(u>0\). Put
   \[
   G(u)=(1-\alpha)h(2u-u^2)+\alpha h(s^*(u,u))-h(u).
   \]
   For \(0<u\le\psi\), the iid gap is nonnegative (split at
   \(p_0=1-1/\sqrt2\) and use entropy symmetry), while the dependent gap is
   strictly positive: \(s^*(u,u)=2u\) for \(u\le1/4\), and \(s^*=1/2\) for
   \(1/4\le u<1/2\). For \(\psi\le u\le t_{\rm cert}\), \(s^*=1/2\) and
   \[
   G'(u)=(1-\alpha)(2-2u)h'(2u-u^2)-h'(u)<0,
   \]
   because \(2u-u^2>1/2\), \(u<1/2\), and \(0<\alpha<1\). This is the sign
   proof in `S/reduction.py:347-377`. A targeted 256-bit Arb evaluation at the
   exact endpoint produced
   \[
   G(t_{\rm cert})=
   [0.0012927838426659849883486212008270239776312579875
   \;\mathbin{+/-}\;4.45\times10^{-76}]>0.
   \]
   Hence \(G(u)>0\) throughout \((0,t_{\rm cert}]\), including the boundary.
   The entropy sum is therefore strict, contradicting the upper bound.

**Resolution and status.** The historical paper still labels the final bridge
as a **CITED-DEPENDENCY** (`uc/paper/main.tex:260-278`), but the implication is
fully reconstructed above at **HUMAN-AUDITED** level. It now depends only on
standard finite Shannon chain rule/conditioning facts and the displayed
sequential coupling, not on Cambie's numerical optimization. The targeted
endpoint sign is **MACHINE-VERIFIED**; the sequential proof is not machine
checked. This also shows that the stronger assumption “all frequencies
\(\le t\)” leads to a contradiction, so the `<`/`<=` mismatch does not leave
a boundary gap.

**2026-08-27 update — the finite content of this bridge is now exhaustively
machine-checked.** `verification/entropy_bridge_exhaustive.py` rebuilds steps
1--3 above from the corollary statement alone, importing neither
`S/bridge_uc.py` nor any campaign or certificate module, and runs them on
**every** nonempty family of subsets of \([n]\) for \(n\le4\): 65,808 families
and 1,631,880 coupled prefix states.  In exact rational arithmetic it confirms
the four conditional masses, both Bernoulli marginals, the OR parameter
\(s^*\), prefix-by-prefix uniformity of \(A\) and \(C\),
\(\mathcal L(P_i)=\mathcal L(R_i)\), \(\mathbb E P_i=\) the frequency of
element \(i\), and the chain rule \(\mathsf H(A)=\sum_i\mathbb E h(P_i)\).  In
256-bit Arb it certifies the two conditioning inequalities, each enclosure
being provably nonnegative or a structural tie inside \(\pm2^{-200}\).  All
5,096 enumerated union-closed families satisfy the \(t_{\rm cert}\)
conclusion, the extreme case being exactly \(1/2\).  The sequential proof for
arbitrary finite families remains **HUMAN-AUDITED**: this is an exhaustive
finite control, and `verification/test_entropy_bridge_exhaustive.py` (4/4)
demonstrates that it rejects a broken \(s^*\) clip, a constant prefix
probability, and an inflated union bit.  Report
`verification/results/entropy-bridge-exhaustive.json`, canonical SHA-256
`6681f8faf13d9344f8e7bf3a5d7785ed289f1e82c5e9badd1e7967c03e3ebfef`.

A second independent read-only referee re-derived this bridge and the
Section 4.2/Theorem B\('''\) chain on 2026-08-27 and found **no proof defect
and no blocker**.  Eight minor exposition defects were repaired in the
manuscript; see
`verification/results/uc-bridge-referee-2026-08-27.md`.

### 4.2 Coupling symmetrization, folding, and existence of the minimum

The two clipped-median formulas agree by three exhaustive cases in
`S/bridge_uc.py:7-20,78-97`. The targeted frozen script checked 4,225 exact
rational pairs and all three cases, but the universal status comes from the
case proof, not the grid.

For arbitrary regular Borel couplings, let \(\tau(p,r)=(r,p)\). If
\(\pi\in\Pi(\mu,\mu)\), then
\(\bar\pi=(\pi+\tau_\#\pi)/2\) is symmetric, has both marginals \(\mu\), and
preserves every symmetric cost. Pushing \(\bar\pi\) through the continuous
fold map \((p,r)\mapsto(\min(p,r),\max(p,r))\) gives \(\nu\in P(\Delta)\).
Conversely, \((\nu+\tau_\#\nu)/2\) unfolds a pair-orbit law. Marginals and
cost agree term by term. The exact finite coefficient schema is
`S/bridge_uc.py:133-166`; the arbitrary-measure argument is
**HUMAN-AUDITED**.

The paper writes a minimum without spelling out its attainment
(`uc/paper/main.tex:249-250,597-605`). This is sound: \(P([0,1]^2)\) is weak-*
compact; the two fixed-marginal conditions are weak-* closed because marginal
test integrals use continuous functions; \(\Pi(\mu,\mu)\) is nonempty via
\(\mu\otimes\mu\); and the continuous bounded cost \(h\circ s^*\) attains its
minimum. Consequently, for every \(\nu\),
\(\Phi_{\rm ex}(\nu)\ge F(\mu_\nu)\), while an optimal coupling for every
\(\mu\) folds to equality. This proves \(\inf\Phi_{\rm ex}=\inf F\).

**Status:** **HUMAN-AUDITED, resolved omission**. It uses the same
**CITED-DEPENDENCY** compactness foundation as D6. No finite-support
assumption is introduced.

### 4.3 Exact decompositions and integrability

The \(T\)-series follows by integrating
\(-\log(1-u)=\sum_{m\ge1}u^m/m\) on \([0,1)\) and using the telescoping
endpoint sum. Its nonnegative coefficients give

\[
\iint T(xy)\,d\eta(x)d\eta(y)
=\sum_{n\ge2}\frac{(\int x^n\,d\eta)^2}{n(n-1)}\ge0
\]

for every finite signed regular measure \(\eta\). The kernel identity is a
direct pointwise algebraic identity, not a conclusion from the three-atom
SymPy sample (`S/decomposition.py:17-42,220-248`).

Similarly,

\[
W(x,y)=\sum_{n\ge2}\frac{(x-x^n)(y-y^n)}{n(n-1)}
\]

converges absolutely and uniformly because each summand has absolute value at
most \(1/[n(n-1)]\). Uniform convergence licenses integration against finite
signed measures and yields

\[
\ln2\,E(\mu)=\sum_{n\ge2}\frac{(M_1-M_n)^2}{n(n-1)}\ge0,
\qquad Q=2BL-E.
\]

This supplies the convergence/integrability sentence that is explicit in
`S/thmB3_proof.py:165-173` but only implicit in parts of the paper.

**Status:** **HUMAN-AUDITED** exact analysis; the frozen scripts' sampled
values and matrix spectra are only **COMPUTATIONAL-EVIDENCE**. Evidence:
`S/decomposition.py:17-42,204-268`, `S/reduction.py:15-45,152-257`, and
`uc/paper/main.tex:304-450`.

### 4.4 Compactness, continuity, fixed-\(B\) concavity, and Bauer

- \(\Delta\) is compact metric. Riesz identifies regular Borel measures with
  \(C(\Delta)^*\), and Banach–Alaoglu makes its dual unit ball weak-* compact.
  Positivity and mass-one conditions are weak-* closed, so \(P(\Delta)\) is
  compact (`S/thmB3_proof.py:26-58`).
- Entropy, \(s^*\), \(T\), \(W\), and the orbit cost are continuous on their
  closed domains. The marginal map is continuous linear. Product test
  functions are Stone–Weierstrass dense; hence
  \(\mu_j\otimes\mu_j\to\mu\otimes\mu\), so \(E\) and
  \(\Phi_{\rm ex}\) are continuous (`S/thmB3_proof.py:61-132`).
- The feasible set \(H_t\) is closed, nonempty, and compact, so a global
  minimizer \(\nu^*\) exists (`S/thmB3_proof.py:135-146`).
- On the slice \(B(\nu)=B(\nu^*)=\beta^*\),
  \[
  \Phi_{\rm ex}=[2(1-\alpha)\beta^*-1]L+\alpha C_\nu-(1-\alpha)E.
  \]
  The first terms are linear. The PSD quadratic \(E\) is convex by the exact
  polarization identity in `S/thmB3_proof.py:175-181`; therefore
  \(-(1-\alpha)E\) is concave exactly when \(\alpha\le1\). The slice stays
  feasible because fixing \(B\) fixes the mean (`S/thmB3_proof.py:149-193`).
- The weak-* dual is locally convex Hausdorff. Bauer applied to the continuous
  convex function \(-\Phi_{\rm ex}\) gives an extreme minimizer
  (`S/thmB3_proof.py:195-207`).

**Status:** hypotheses and signs **HUMAN-AUDITED**. Riesz,
Banach–Alaoglu, Stone–Weierstrass, and Bauer are **CITED-DEPENDENCY**. Exact
quotes and source metadata are in `S/thmB3_proof.py:42-47,112-115,195-201,
283-303`. The hard-coded quote strings checked at
`S/thmB3_proof.py:326-334,534-553` do not machine-authenticate those external
sources.

### 4.5 Extreme-point support reduction

If an extreme point of the one-moment slice had three support points, choose
three disjoint positive-mass neighborhoods. Their two-vectors
\((\nu(A_i),\int_{A_i}g\,d\nu)\) are linearly dependent. The resulting
nonzero signed restriction \(\eta\) has zero mass and zero \(g\)-moment, and
for sufficiently small \(\varepsilon\), both
\(\nu\pm\varepsilon\eta\) remain nonnegative slice members. This contradicts
extremality. The density bound is explicit:
\(\varepsilon<1/\max_i|a_i|\) suffices.

This direct proof covers atomic and non-atomic measures because any three
points of the topological support have disjoint neighborhoods of positive
measure. It does not assume that a measure already has finite support.

**Status:** **HUMAN-AUDITED** at `S/thmB3_proof.py:251-265` and
`uc/paper/main.tex:710-728`. Winkler's theorem via Pinelis
(`S/thmB3_proof.py:220-249`) is corroborative, not logically load-bearing;
the trust base can omit Winkler/Pinelis if the direct proof is retained.
Bauer remains load-bearing.

### 4.6 Margin Lemma, feature-map integrability, and \(L=0\)

The feature map may be chosen componentwise as

\[
\mathbf v_n(p)=
\frac{(1-p)-(1-p)^n}{\sqrt{\ln2\,n(n-1)}}.
\]

The uniformly bounded square-summable tail makes
\(p\mapsto\mathbf v(p)\) continuous into \(\ell^2\), hence strongly
measurable and Bochner integrable against every probability measure. Its norm
satisfies \(\|\mathbf v(p)\|^2=\operatorname{rh}(p)\), so

\[
E(\mu)=\left\|\int\mathbf v\,d\mu\right\|^2
\le\left(\int\|\mathbf v\|\,d\mu\right)^2=S^2.
\]

The coefficient \(1-\alpha\ge0\) then gives (1). Equality with one non-sink
support point is correct because all nonzero feature vectors in the integral
are the same vector direction. Evidence: `S/margin_lemma.py:9-37` and
`uc/paper/main.tex:778-825`.

The Margin Lemma quotient assumes \(L>0\). If \(L=0\), nonnegativity and
strict interior positivity of \(h\) force
\(\mu=a\delta_0+(1-a)\delta_1\). Then \(Q=C=E=S=F=0\); every pair-orbit
representing this marginal is supported on \((0,0),(0,1),(1,1)\), for which
\(h(s^*)=0\), so \(\Phi_{\rm ex}=\Phi_{\rm rel}=0\). No division or limiting
argument is needed.

**Status:** **HUMAN-AUDITED, resolved boundary omission**. Source support:
`S/reduction.py:48-54`, `S/margin_lemma.py:20-37`, and
`S/arbcore.py:29-43,77-104`.

### 4.7 The \(\operatorname{rh}\), \(\rho\), and derivative endpoint lemmas

For \(0\le x\le1/2\), the exact identity

\[
\ln2\,[2x-\operatorname{rh}(x)]
=x^2\ln(2/x)-x(2-x)\ln(1-x/2)
\]

has two nonnegative terms, strictly positive for \(x>0\). The frozen proof
also supplies an explicit near-zero remainder and an exact contiguous Arb
cover of \([1/16,1/2]\). A targeted `python -B S/lemma_rh_proof.py` run during
this audit printed 44 accepted subintervals, minimum lower bound
`0.002407365711405873...`, and `PROOF COMPLETE: machine-checked`. The right
half follows from \(0\le h\le1\). Thus
\(\operatorname{rh}(x)\le2\min(x,1-x)\), including equality-compatible
endpoints. Evidence: `S/lemma_rh_proof.py:42-207,266-424`.

For \(0<z<1\),

\[
h(z)-zh'(z)=-\log_2(1-z)>0,
\]

so \(z/h(z)\) increases and, by entropy symmetry,
\((1-z)/h(z)\) decreases. The exact SymPy identities are asserted at import
in `S/cert2.py:51-65`; the sign inference is human. Therefore the endpoint
collar caps used in `rho_upper` are valid. Interior cells use

\[
\rho(z)\le 2(1-u)-
\frac{\inf_{[u,v]}h(2z-z^2)}{\sup_{[u,v]}h(z)},
\]

which is dependency-safe because it combines one-sided nonnegative bounds.
The 20,000-cell global cap and the per-box minimum of valid caps are in
`S/cert2.py:202-282`.

For the \(q_2\to1\) derivative, write \(x=1-q_2\). Then

\[
\ln2\,\operatorname{rh}(1-x)
=(1-x)^2\ln(1-x)+(1-x^2)\ln(1+x)=F(x).
\]

On \(0\le x\le X=1/16\), elementary integral remainder bounds give
\(F(x)\ge c x^2\), where \(c=(2-3X+X^2)/2>0\). Also

\[
F'(x)=2[(1-x)(-\ln(1-x))-x\ln(1+x)],
\]

and
\(0\le F'(x)\le2x/(1-X)\): the lower sign follows from
\((1-x)(-\ln(1-x))\ge x(1-x)\) and
\(x\ln(1+x)\le x^2\); the upper bound drops the nonnegative second term and
uses \(-\ln(1-x)\le x/(1-X)\). Hence

\[
\left|\frac{d}{dx}\sqrt{\operatorname{rh}(1-x)}\right|
\le\frac1{(1-X)\sqrt{c\ln2}}=\texttt{DQ_ENDPOINT}.
\]

The estimate is sound. It is, however, **HUMAN-AUDITED**, not fully
machine-proved by the import-time assertions: `S/cert3.py:126-171` asserts
three symbolic derivative/value identities and a numerical Arb range for the
constant, but not every Taylor/order inequality. The message “Endpoint
series/remainder proof PASS” at `S/cert3.py:1634-1639` is therefore overstated.

## 5. Verifier-to-theorem parameter correspondence

### 5.1 Frozen parameters

| Mathematical/code parameter | Frozen Campaign I value | Exact use and evidence | Audit result |
|---|---:|---|---|
| \(\alpha\) | `0.0356069` | Hard-coded Arb enclosure `S/arbcore.py:12`; manifest `I/launch.json:1` | Encloses exact rational; **HUMAN-AUDITED** correspondence |
| \(t\) | `0.3820660112501052` | Manifest and `_arb(t)` at `S/cert3.py:1421-1424` | Exact decimal rational; target relation **MACHINE-VERIFIED** |
| Root | \([0,1]^4\times[1/2,1]\) | `S/cert3_par.py:154-157,553-555` | Full ordered domain after orbit swap; **HUMAN-AUDITED** |
| Slices | Eight width-\(1/16\) dyadic intervals | `I/launch.json:1`; collector exact-Fraction check `S/cert3_collect.py:197-234` | No gap/overlap; **MACHINE-VERIFIED** protocol check |
| Ordinary/collar floor | `0.001` / `0.000125` | Worker overrides `cert3.py`'s unused default: `S/cert3_par.py:61-68,567-581` | Correct manifest-to-worker mapping |
| Face floor/budget | `0.001` / `100000` | `S/cert3_par.py:65-67,573-575`; replay `S/cert3_replay.py:117-132,289-299` | Failure only prevents a clear; sound |
| Box budget | `2000000000` | `S/cert3_par.py:67,570` | Zero budget tally required; sound |
| Time budget | `259200` s = 72 h per slice | Every run in `I/launch.json:1` | `uc/PROOF.md:279-282` saying 24 h is stale |
| Native/work precision | 160 / 80 bits | `S/cert3_replay.py:32-34,135-161`; manifest work precision | Correct import/setup order |
| Full face \(\lambda\) family | `0,1,1.3408514938246723,1.3682158100251758,1.3955801262256793,1.626,2.2` | `I/launch.json:1`; `S/cert3.py:449-493` | KKT accuracy is heuristic only; every represented \(\lambda\ge0\) is valid |
| Center-rule \(\lambda\) family | `0` and `1.3682158100251758` | `S/cert3.py:1430-1435`; replay `S/cert3_replay.py:148-161` | Narrower than the paper's phrase “the lambda family,” but sound |
| `q2pin` | `15/16` | `S/cert3.py:38,377-427` | Matches derivative lemma domain |
| Center gates | width \(\le1/32\), corner \(\ge-0.02\) | `S/cert3.py:1032-1035,1495-1528` | Heuristic gates can only skip a rule |
| Code/checker/launch hashes | `2f23a58e...ce05`, `95d09322...81ec`, `6414affc...65cc` | `I/collect_output.txt:1,10-12`; `I/launch.json:1` | Must be externally bound together; see Section 7.3 |

A lightweight environment check during this audit returned CPython `3.14.3`,
flint `0.9.0`, mpmath `1.3.0`, NumPy `2.5.2`, SciPy `1.18.0`, SymPy `1.14.0`,
and `macOS-26.5.2-arm64-arm-64bit-Mach-O`, matching `I/launch.json:1`.
This is **COMPUTATIONAL-EVIDENCE** about the current interpreter; the collector
itself does not enforce the environment field.

### 5.2 Interval primitives

| Theorem quantity | Exact verifier implementation | Soundness correspondence | Status |
|---|---|---|---|
| \(h([u,v])\) | `h_pt`, `h_encl`, `S/arbcore.py:29-43` | Endpoint values exact; endpoint extrema plus 1 if crossing \(1/2\) | **HUMAN-AUDITED**; Arb arithmetic machine |
| min/max and \(s^*\) | `amin`, `amax`, `sstar_encl`, `S/arbcore.py:46-74` | Outward hull on incomparable balls; monotone two-corner range | **HUMAN-AUDITED** |
| \(\operatorname{rh}\), \(\sqrt{\operatorname{rh}}\) | `rh_encl`, `sqrt_rh_enc`, `S/arbcore.py:77-104` | Natural enclosure, clamp uses proved \(\operatorname{rh}\ge0\), global upper cap only tightens | **HUMAN-AUDITED** plus Arb |
| Global \(\operatorname{rh}\) cap | `_rh_global`, `S/arbcore.py:107-130` | 20,000 covering cells; process-local lazy cache | **MACHINE-VERIFIED** when recomputed; no disk-cache trust |
| Global/per-box \(\rho\) cap | `_rho_global`, `rho_upper`, `S/cert2.py:221-282` | Minimum only among independently valid upper caps | **HUMAN-AUDITED** plus Arb |
| Float endpoint semantics | `_arb_float_point`, `_decimal_fraction`, `S/diag_exhaust.py:140-160,210-242` | Every proof interval encloses the shortest-decimal rational; contractor computes that rational exactly | **HUMAN-AUDITED** |

The independent replay map exposed one conversion boundary:
\path{cert2.rho_upper} calls `arb(float)` for atom endpoints, whereas the
declared box convention elsewhere is `arb(repr(float))`.  For Campaign I this
does not create a different atom interval.  Atom coordinates start at \(0,1\),
only midpoint bisection changes them, and the ordinary/collar floors bound
their depths by \(10/13\).  An exhaustive exact-Fraction regression checks
every \(k/2^d\), \(0\le d\le13\), and finds
`Fraction(repr(float(k/2**d))) == Fraction(k,2**d)`.  The mean contractor
changes only \(w\), whose ratio-rule window is explicitly constructed through
`repr`.  The standalone verifier nevertheless uses the outward decimal
conversion consistently, removing this campaign-specific invariant from its
arithmetic implementation.

### 5.3 Every trace discharge and transition rule

In this table, “machine” means the archived trace replay recomputed the named
function in Arb/Fraction arithmetic. The assertion that the function is a
valid theorem bound is separately human-audited.

| Trace code | Rule and exact theorem implication | Preconditions and sign checks | Function/lines | Status |
|---:|---|---|---|---|
| 0 | From orbit lower means \(l_i\), \(M\ge l_2+w(l_1-l_2)\); exact solution of the necessary inequality contracts \(w\) or proves no feasible point | Uses an outward upper endpoint for \(t\); decimal cutoffs stepped outward | `mean_contract`, `S/diag_exhaust.py:210-286` | Bound **HUMAN-AUDITED**; event arithmetic/replay **MACHINE-VERIFIED** |
| 1 | Prepared ordered rectangle or corner mean enclosure is infeasible | `phi_corner` returns `None` only for empty ordered part or `mean.lower()>t` | `S/diag_exhaust.py:163-185`; replay `S/cert3_replay.py:238-244` | **HUMAN-AUDITED + MACHINE-VERIFIED** |
| 2 | \(\Phi_{\rm rel}\ge\kappa_*L_{lo}+\alpha C_{lo}-(1-\alpha)S_{hi}^2\ge0\) | \(m_*=\min(M_{hi},t)\); requires \(\kappa_*>0\); all three aggregates nonnegative | `phi_corner`, `S/diag_exhaust.py:163-207` | **HUMAN-AUDITED + MACHINE-VERIFIED** |
| 3 | Drop \(\alpha C_\nu\ge0\), use \(M\le t\) and \(S^2\le L\int\rho\,d\mu\), then prove \(\kappa_t-(1-\alpha)\rho_{hi}\ge0\) | Own contracted \(w\)-window; convex combination of one-sided caps maximized at its two endpoints | `ratio_rule`, `S/cert2.py:285-313` | Corrected inequality direction **HUMAN-AUDITED**; replay **MACHINE-VERIFIED** |
| 4 | Full five-variable mean-value lower bound for \(\Phi_\lambda\) | Prepared rectangle; all atom intervals strict interior; one certified affine \(s^*\) branch; \(\lambda\ge0\) | `interval_gradient`, `S/bound_kkt.py:337-400`; `centered5_best`, `S/cert3.py:1096-1186` | Gradient formulas **HUMAN-AUDITED**; Arb leaf **MACHINE-VERIFIED** |
| 5 | MVT in \((p_1,q_1,w)\), natural interval anchor in \((p_2,q_2)\) | Differentiated orbit strict interior; anchor encloses all of the second-orbit rectangle | `centered5_mixed`, `S/cert3.py:1189-1274` | **HUMAN-AUDITED + MACHINE-VERIFIED** |
| 6 | Apply code 5 after exact orbit swap \((A,B,w)\mapsto(B,A,1-w)\) | Arb subtraction gives outward decimal-semantic weight endpoints | `orbit_swap_box`, `centered5_mixed_swap`, `S/cert3.py:1277-1299` | Invariance **HUMAN-AUDITED**; endpoint enclosure/replay **MACHINE-VERIFIED** |
| 7 | One-dimensional MVT in \(w\), with a natural interval anchor over all atom coordinates | No atom derivative, so finite at entropy sinks; exact SymPy weight derivative identity asserted at import | `S/cert3.py:1302-1320,1323-1394` | Identity **MACHINE-VERIFIED**, bound **HUMAN-AUDITED**, leaf **MACHINE-VERIFIED** |
| 8 | If \(\partial_{q_2}\Phi_\lambda\le0\) on the entire strip, then \(\Phi_\lambda(q_2)\ge\Phi_\lambda(1)\); a complete full-face B&B proves the latter nonnegative | \(q_{2,lo}\ge15/16\), \(p_{2,hi}\le q_{2,lo}\), strip coefficient positive, \(\lambda\ge0\); **no feasibility restriction after pinning** | `q2_pin`, `S/cert3.py:377-427`; `face_best_bound`, `face_bb`, `S/cert3.py:331-355,990-1014` | Derivative estimate/bound **HUMAN-AUDITED**; nested face arithmetic/replay **MACHINE-VERIFIED** |
| 9 | Uncleared terminal residual | All coordinates at applicable floors | `S/cert3.py:1547-1562`; replay always rejects at `S/cert3_replay.py:300-303` | **FAILED** as proof event; Campaign I count is zero |
| 16+\(j\) | Split coordinate \(j\); children share a strict representable midpoint and cover the contracted parent exactly | Choice by influence is cost-only; replay need not verify heuristic ordering or floor | `S/cert3.py:982-987,1038-1060,1563-1566`; replay `S/cert3_replay.py:304-317` | Coverage **HUMAN-AUDITED + MACHINE-VERIFIED** |

The \(q_2\)-pin inequality direction was checked explicitly. Its exact
derivative is

\[
(\kappa v/2+\alpha v)h'(q_2)
-(1-\alpha)vL-(1-\alpha)Sv(\sqrt{\operatorname{rh}})'(q_2)
+\lambda v/2.
\]

`q2_pin` upper-bounds the first negative term using lower nonnegative
coefficients and the decreasing endpoint value of \(h'\), drops the
nonpositive \(-(1-\alpha)vL\), and bounds the remaining possible positives by
`DQ_ENDPOINT` (`S/cert3.py:409-427`). Thus derivative `<= 0` implies the face
value at 1 is a lower bound, not an upper bound. The face code correctly
forbids a feasibility-restricted face bound because pinning can increase the
mean (`S/cert3.py:331-355,1000-1003`).

### 5.4 Campaign rule-use tallies

Summing the eight committed one-line result records gives:

| tally | total |
|---|---:|
| processed | 488,465,854 |
| split | 244,232,923 |
| infeasible (codes 0 and 1 combined) | 43,494,129 |
| corner | 66,762,785 |
| ratio | 2,422,979 |
| centered5 | 83,634,710 |
| mixed | 14,335,937 |
| mixed-swap | 3,250,523 |
| weight-only | 97,084 |
| face | 30,234,784 |
| residual / final stack / box budget / time budget | 0 / 0 / 0 / 0 |

Every positive discharge category was therefore exercised by Campaign I. The
terminal total is `244,232,924 = split + 1`, and
`processed = terminals + split`, both exact for the aggregate and separately
for every slice. Evidence is the eight files named by `I/launch.json:1`, each
at line 1:

- `I/result_slice0_2fb154187bc54147855004d5dafaf214.json:1`
- `I/result_slice1_af8b2ece3eb947339f8a4e273c46e3bd.json:1`
- `I/result_slice2_26a971c6ac644309946431b3900b780a.json:1`
- `I/result_slice3_89cda161a3584c8aa9c1631258e37fde.json:1`
- `I/result_slice4_5136a0e6ac8c4078a5f8ceb73e34b777.json:1`
- `I/result_slice5_7ab1c417449946a7b48c2d66c672bf93.json:1`
- `I/result_slice6_30d759beee9b4026928e3e9181f6afaa.json:1`
- `I/result_slice7_e8b0b30179f441129e30d9afa266cc78.json:1`

These counts are **COMPUTATIONAL-EVIDENCE** about rule use; the archived
collector's successful recomputation of their corresponding trace events is
the narrow **MACHINE-VERIFIED** claim at `I/collect_output.txt:2-25`.

## 6. Boundary-case matrix

| Boundary or degeneracy | Required conclusion | Verifier/proof handling | Status |
|---|---|---|---|
| \(h(0)=h(1)=0\) | No `log(0)` evaluation; integrals continuous | Exact endpoint branches in `S/arbcore.py:29-43`; analytic continuity `S/thmB3_proof.py:82-88` | **HUMAN-AUDITED** |
| \(L=0\) | Functional is zero without division | Endpoint-support argument in Section 4.6 | **HUMAN-AUDITED, resolved** |
| \(\operatorname{rh}=0\) at 0 and 1 | Square root finite; no derivative assumed | Clamp and endpoint-safe `sqrt_rh_enc`, `S/arbcore.py:77-104`; centered atom derivatives return `None` and fall back | **HUMAN-AUDITED + machine arithmetic** |
| \(\rho=\operatorname{rh}/h\) at 0 or 1 | Never divide by zero-containing entropy interval | Collar bounds in `S/cert2.py:193-282` | **HUMAN-AUDITED** |
| \(p_i=q_i\) | Ordered triangle includes diagonal | Closed inequalities; pair contributes the same endpoint twice | **HUMAN-AUDITED** |
| Raw box overlaps \(p_i>q_i\) | Do not lose ordered part | Raise only `q_lo` to `max(q_lo,p_lo)`, yielding a rectangular superset; empty intersections are infeasible | `S/bound_kkt.py:123-136`; **HUMAN-AUDITED** |
| \(p+r=1/2\), \(p+r=1\), \(\max(p,r)=1/2\) | Cost continuous across branch ties | Min/max formula; non-smooth centered rules either certify one branch or use a Lipschitz cap/skip | `S/bridge_uc.py:15-20`; `S/bound_kkt.py:312-334`; **HUMAN-AUDITED** |
| \(M=t\) | Feasible boundary retained; shift is exactly zero | Contractor uses permissive `t_upper`; all comparisons are non-strict where required | `S/diag_exhaust.py:245-286`; **HUMAN-AUDITED** |
| \(\kappa_*=0\) | Corner lower-L multiplication cannot be used | Rule explicitly requires `coeff_star > 0`, otherwise returns negative/open | `S/diag_exhaust.py:185-188`; **MACHINE-VERIFIED branch** |
| \(w=0\) or 1 | One orbit may be vestigial; all laws still covered | Parameterization includes endpoints; mixed/weight/ratio rules remain valid; q2 pin with \(v=0\) is equality | `S/cert3.py:377-427,1189-1394`; **HUMAN-AUDITED** |
| \(w=1/2\) | Orbit-swap fundamental domains meet without a gap | Exact dyadic slice endpoint and involution | `S/cert3_par.py:154-157`; **HUMAN-AUDITED** |
| \(w<1/2\) | Must not be omitted | Swap orbit labels and send \(w\mapsto1-w\) | `S/cert3.py:1277-1299`; **HUMAN-AUDITED** |
| Coincident two orbit atoms | One-orbit minimizer included; continuity retained | Weights/points may coincide; no division by their separation | `S/thmB3_proof.py:8-20,460-512`; **HUMAN-AUDITED** |
| \(q_{2,lo}=15/16\) | Endpoint derivative cap applies inclusively | `dq_upper` accepts `>=`; `q2_pin` rechecks branch strip | `S/cert3.py:166-171,377-401`; **HUMAN-AUDITED** |
| \(q_2\to1\) | Entropy derivative diverges but \((\sqrt{\operatorname{rh}})'\) stays bounded | Separate analytic `DQ_ENDPOINT`; integration to endpoint uses continuity | Section 4.7; **HUMAN-AUDITED** |
| Pin may make \(M>t\) | Never use a feasible-only face candidate | `allow_feasible=False` in `face_bb` | `S/cert3.py:331-355,990-1003`; **MACHINE-VERIFIED path** |
| Decimal endpoint after cancellation | Binary `nextafter` alone is not trusted | Exact `Fraction(repr(float))` cutoff and outward decimal stepping | `S/diag_exhaust.py:210-286`; **HUMAN-AUDITED** |
| Split midpoint rounding | Children must cover and make progress | Shared Python-float midpoint; replay requires `lo < mid < hi` | `S/cert3_replay.py:304-317`; **MACHINE-VERIFIED** |
| All widths at floor | Cannot silently accept | Emit residual code 9; collector/replay rejects; Campaign I has zero | `S/cert3.py:1547-1562`; **MACHINE-VERIFIED** |
| \(\alpha=1\) in support theorem | Concave energy term vanishes | Allowed analytically; Campaign uses interior \(\alpha_0\) | `S/thmB3_proof.py:18-20,190-193`; **HUMAN-AUDITED** |
| \(t=0,1\) in support theorem | Feasible set nonempty for analytic quantifier | \(\delta_{(0,0)}\) works; certificate fixes interior target | `S/thmB3_proof.py:138-146`; **HUMAN-AUDITED** |
| First union-closed coordinate \(u=0\) | Strict contradiction would fail | Choose a present element first; \(\mathcal F\ne\{\varnothing\}\) guarantees \(u>0\) | Section 4.1; **HUMAN-AUDITED** |
| First coordinate \(u=t_{\rm cert}\) | Resolve Q2 `<` versus local `<=` | Exact derivative sign plus positive 256-bit Arb endpoint value | Section 4.1; **MACHINE-VERIFIED endpoint + HUMAN-AUDITED interval** |

## 7. Trace replay, campaign acceptance, and trust boundaries

### 7.1 What replay proves

`replay_trace` starts from the exact dyadic root, reapplies the exact mean
contractor, and consumes one byte per raw DFS node. It re-evaluates only the
claimed terminal rule; optional rules need not have failed before a split.
It rejects unknown events, a false infeasibility claim, residual code 9, a
nonstrict split, trailing bytes, a nonempty final stack, and any mismatch in
the committed tallies (`S/cert3_replay.py:187-331`). This is sufficient: each
terminal is independently valid and each split exactly covers its parent.
Rule order is a performance choice, not a theorem premise.

Replay imports `cert2`, `cert3`, `arbcore`, `bound_kkt`, `diag_exhaust`, and
`entropy` from the frozen directory and verifies their module locations
(`S/cert3_replay.py:19-29,50-63`). It recomputes the global covers at 160 bits,
switches to 80-bit work precision, and recomputes/checks the lambda family
(`S/cert3_replay.py:104-161`).

**Scope limitation:** replay calls the same `mean_contract`, `phi_corner`,
`ratio_rule`, centered, pin, and `face_bb` functions as the producer
(`S/cert3_replay.py:224-299`). It is independent of the producer's trace
claims and search-control path, but not independent source code or an
independent mathematical implementation. Common-mode rule bugs are excluded
by the human source audit, not by replay alone.

### 7.2 Standalone arithmetic implementation and adversarial review

`verification/independent_arithmetic_replay_secure.py` reimplements every
terminal predicate without importing or executing the frozen Python
arithmetic.  Its interval automatic-differentiation layer derives the
five-variable, mixed, weight-only, and face gradients from the displayed
functional; the \(q_2\) endpoint lemma and nested face B\&B are separately
implemented.  It shares only the authenticated outer trace partition and the
python-flint/Arb primitive backend.

The first standalone version failed closed at trace event 303,061, opcode 5.
Wrapping adjacent exact interval endpoints in a second Arb hull obscured the
certified \(q_1^{lo}=p_1^{hi}\) branch tie and widened one gradient enough to
lose a \(2.2\times10^{-6}\) margin.  Preserving the declared endpoint balls
recovered a positive independent bound
\(2.2277095\times10^{-6}\), close to but not copied from the frozen
\(2.2278896\times10^{-6}\) result.  This was an independent-verifier strength
failure, not a false campaign terminal; it is now a permanent regression.

An adversarial code review then found a real acceptance flaw: the first
version's checkpoint HMAC key was public, so a forged EOF checkpoint could
inject public expected tallies and obtain `PASS`.  That version's ongoing
byte-zero run is retained only as diagnostic evidence.  The secure version:

- forbids resumed checkpoints from producing full `PASS`;
- discloses and requires `initial_trace_offset=0`, `resume_used=false`;
- binds and rechecks the exact source, inputs, target, precision, runtime seal,
  eight roots, tallies, topology, and live no-resume commands;
- requires an out-of-band trusted raw digest for a report lock fixing all
  reports, the launch record, and the live process/image attestation.

The 343-file prelaunch runtime seal is supplemented by a mid-run `ps`/`lsof`
inventory: all eight exact commands were live, 62 non-system images per worker
matched the seal, two additional Python framework images were explicitly
content-pinned, and `/usr/lib/dyld` was the sole unsealed system image.  A
second adversarial recheck found no remaining material false-`PASS` path under
the explicitly retained trusted assumptions: custody of the out-of-band
report-lock digest, a nonmalicious OS/kernel/process-inspection/hash stack, and
absence of a hostile same-inode file-rewrite race.

The secure eight-slice byte-zero replay is now **MACHINE-VERIFIED** under that
trust boundary.  All eight processes exited zero after independently
recomputing \(488{,}465{,}854\) events, \(244{,}232{,}923\) splits, and
\(244{,}232{,}931\) leaves with zero residuals and exact locked tallies.
The externally trusted report-lock raw SHA-256 is
`af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
the accepted composite canonical SHA-256 is
`4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.
Exact per-slice hashes, runtimes, commands, runtime pins, adversarial tests,
and the remaining OS/custody assumptions are in `uc/REPRODUCIBILITY.md`.

### 7.3 Campaign I acceptance record

| slice | exact \(w\)-interval | processed | split | result status | archived collector replay |
|---:|---|---:|---:|---|---|
| 0 | \([1/2,9/16]\) | 21,135,611 | 10,567,805 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:2` |
| 1 | \([9/16,5/8]\) | 20,827,065 | 10,413,532 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:3` |
| 2 | \([5/8,11/16]\) | 23,486,265 | 11,743,132 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:4` |
| 3 | \([11/16,3/4]\) | 34,186,855 | 17,093,427 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:5` |
| 4 | \([3/4,13/16]\) | 71,405,987 | 35,702,993 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:6` |
| 5 | \([13/16,7/8]\) | 117,709,435 | 58,854,717 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:7` |
| 6 | \([7/8,15/16]\) | 121,456,775 | 60,728,387 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:8` |
| 7 | \([15/16,1]\) | 78,257,861 | 39,128,930 | COMPLETE; all four zero tallies | PASS, `I/collect_output.txt:9` |

The collector checks exact run IDs and roots, a pristine 15-name snapshot,
source/result/trace hashes, trace size equal to processed count, normal worker
exit 0, COMPLETE, positive work, and
`stack=residual=budget_boxes=budget_time=0`; it then replays all traces before
printing `COMPOSITE CERTIFICATE` (`S/cert3_collect.py:139-285,306-358,
361-484,487-516`). The archived output names code hash
`2f23a58ebdb8b14275284e2ffefca7862373f7baf5d354a06e607cff8d78ce05`
and checker hash
`95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec`
(`I/collect_output.txt:1,10-25`).

### 7.4 Acceptance limitations that must remain in the claim envelope

1. **Review proofs are hash-pinned, not run by the collector.** The collector's
   replay imports only the executable rule modules. It does not execute
   `bridge_uc.py`, `decomposition.py`, `reduction.py`, `thmB3_proof.py`,
   `margin_lemma.py`, or `lemma_rh_proof.py`; it merely checks their filenames
   and hashes (`S/cert3_collect.py:42-64,427-484`; imported-module set
   `S/cert3_replay.py:22-29,50-63`). Therefore `COMPOSITE CERTIFICATE` is not a
   machine check of the analytic implication chain.

2. **The generic collector does not hard-code the accepted code digest.** It
   hard-codes the 15 filenames, then checks actual hashes against hashes
   supplied by `launch.json` and checks the supplied combined hash
   (`S/cert3_collect.py:237-280`). It never requires the externally audited
   `2f23...ce05` digest. A self-consistent different snapshot and rewritten
   manifest could be accepted by this generic checker. The specific Campaign I
   claim is bound only when the submission/verifier independently requires the
   three printed code/checker/launch hashes. SHA-256 supplies integrity after
   that binding, not authorship or authenticity.

3. **The collector does not compare the current runtime to
   `launch.environment`.** Workers did compare `environment_record()` to the
   launch at start and after execution (`S/cert3_par.py:399-439,532-588,
   659-704`), and the current environment matched in this audit. But
   `S/cert3_collect.py:139-285` contains no environment comparison, so a fresh
   replay uses the current external libraries. This is a reproducibility/trust
   limitation, not an inequality error.

4. **The historical campaign did not content-hash external binaries.** A new
   observed-runtime seal now hashes the resolved CPython executable, every
   imported module file, every regular file in the python-flint distribution,
   and records `otool -L` dependencies: 257 files, report
   `verification/results/runtime-environment-seal.json`, canonical SHA-256
   `53cb75c28ffcc806dc547602f2f82611f50147c11dc8ca4ca821ccdc454a95eb`.
   This narrows byte drift but is not a sandbox or authenticated execution
   attestation; dyld-shared-cache system libraries, the kernel, firmware, and
   hardware remain identified only by platform metadata.

5. **No persistent numerical cache is trusted.** The \(\operatorname{rh}\) and
   \(\rho\) global caps are process-local lazy values recomputed from Arb
   covers (`S/arbcore.py:107-130`, `S/cert2.py:221-253`). `-B`, exact snapshot
   listings, and module-location checks exclude snapshot bytecode caches. The
   remaining startup/import behavior is explicitly in the CPython trust base.

6. **The archived text is evidence, not a signature.** `collect_output.txt` and
   `campaigns/CERTIFICATE.txt` are byte-identical transcripts
   (`I/collect_output.txt:1-25`, `uc/campaigns/CERTIFICATE.txt:1-25`). The
   closeout script structurally checked that text but deliberately did not
   replay (`uc/closeout.sh:1-35,68-85`). A third party obtains machine evidence
   by running the frozen collector under the required hashes/environment.

## 8. Active fault search

| Requested failure mode | Audit test | Finding and resolution |
|---|---|---|
| Circularity | Followed D0–D17 and checked exact versus relaxed uses | No cycle. Support reduction is for \(\Phi_{\rm ex}\); only afterward is \(\Phi_{\rm rel}\) introduced. |
| Missing feasible region | Checked ordering, orbit swap, all \(w\) endpoints, contractor, and trace stack | No gap. Each orbit may be internally ordered; orbit-label swap covers \(w<1/2\); exact dyadic slices cover \([1/2,1]\). |
| Wrong inequality direction | Re-derived (1), lambda shift, corner coefficient, ratio rule, and pin derivative | All load-bearing directions are correct. Full-face rather than feasible-face pinning is essential and present. |
| Float/rational conversion hazard | Compared target, roots, mean cutoff, split, and swap semantics | Target relation and contractor are exact rational/outward. Arb encloses shortest decimals. Shared split endpoints preserve decimal-semantic coverage. |
| Unsound interval dependency | Inspected `rh_encl`, direct \(\rho\), aggregate \(w\), and natural interval formulas | Natural dependency loss only widens bounds. The formerly unsound one-sided slope expression is not used; endpoint maximum is used at `S/cert2.py:301-309`. |
| One-sided cap misuse | Traced every `rho_upper` combination | Caps are combined only with nonnegative weights and maximized at the two \(w\) endpoints. No subtraction of an upper cap occurs. |
| Cache trust | Checked lazy caps, `-B`, pristine listing, and module locations | No disk numerical cache is accepted. CPython startup and installed extension modules remain trusted. |
| Replay independence | Compared worker and replay call graphs | “Independent” is limited to trace/control-flow reconstruction; mathematical functions are shared. Claim envelope corrected in Section 7. |
| Stale contradictory wording | Compared proof, paper, manifests, and prior result records | Multiple non-load-bearing contradictions listed in Section 9. |
| Campaign acceptance | Reconciled all result records, exact tree invariants, output hashes, and collector gates | Static records are internally consistent. Generic collector still needs external code-hash anchoring. |

## 9. Discrepancies, severity, and resolution

| ID | Severity | Status | Discrepancy | Risk | Resolution / exact envelope |
|---|---|---|---|---|---|
| U-01 | **HIGH** | **RESOLVED in manuscript** | The earlier paper did not distinguish the exact pair-orbit functional from the certificate's relaxation. | A reader could infer an invalid pointwise bound for \(F\). | The certificate section now defines both functionals and proves \(\Phi_{\rm ex}-\Phi_{\rm rel}=(1-\alpha)(S^2-E)\ge0\), including \(L=0\). |
| U-02 | **HIGH** | **RESOLVED in manuscript/status documents** | “Machine-checked functional inequality” and broad downstream wording overstated the collector's scope. | Submission could claim stronger verification than exists. | The theorem is now “certificate-backed”; the proof-status box and every root status call the result hybrid and not end-to-end machine checked. |
| U-03 | **MEDIUM** | **RESOLVED** | Primary Q2 says \(\mathbb E p<c\); earlier local summaries said “at most.” | Boundary/strictness gap if left unexplained. | The closed-domain theorem, monotonicity argument, and positive independent Arb endpoint check cover \(\le t\). |
| U-04 | **MEDIUM** | **RESOLVED for the accepted replay; residual trust stated** | The historical collector hard-codes names rather than source digests. | `COMPOSITE CERTIFICATE` alone does not identify the implementation. | The external campaign lock, structural auditor, direct summary, and clean-room verifier pin and recheck the frozen code/checker/launch bytes. The collector's standalone limitation remains explicit. |
| U-05 | **MEDIUM** | **RESOLVED at the declared trust boundary** | The historical collector does not enforce `launch.environment`. | A replay could silently use another interval/library implementation. | The clean-room verifier enforces the exact interpreter/dependency versions and top-level module hashes. System-site/native-library trust and lack of sealed attestation remain explicit limitations. |
| U-06 | **MEDIUM** | **RESOLVED human lemma** | Optimal-coupling existence was unstated where `min` was introduced. | Equality of the \(F\) and pair-orbit optimizations otherwise had a missing attainment step. | The revised paper proves compactness of \(\Pi(\mu,\mu)\) and continuity of the cost before using `min`. |
| U-07 | **MEDIUM** | **RESOLVED human lemma** | The Margin Lemma assumes \(L>0\), but the chain did not isolate \(L=0\). | Division by zero at sink laws. | The revised proof proves all exact and relaxed terms vanish on endpoint-supported laws. |
| U-08 | **MEDIUM** | **RESOLVED claim scope** | Frozen `S/cert3.py` prints an endpoint remainder “proof PASS,” although its import assertions do not execute every order inequality. | The automation behind the \(q_2\)-pin could be overstated. | The final audit and paper classify the derivative/order argument as **HUMAN-AUDITED**; Arb evaluates only the resulting constant. The historical source label is not used as theorem evidence. |
| U-09 | **LOW** | **RESOLVED in manuscript** | Earlier wording said orbit-swap invariance was “asserted exactly,” while the helper includes sampled value controls. | Verification-mode misdescription. | The final paper proves invariance as a termwise human identity and describes replay's exact swapped-bound checks separately. |
| U-10 | **LOW** | **RESOLVED in manuscript** | The earlier Campaign G caption gave stack 24 for slice 7. | Historical record inconsistency only. | The final caption records stack 16, matching the committed result and `uc/PROOF.md`. |
| U-11 | **LOW** | **RESOLVED in proof ledger** | A generic campaign sentence formerly said 24 h per Campaign I slice. | Campaign I wall-time misstatement. | `uc/PROOF.md` now records 72 h/slice, matching `I/launch.json`. |
| U-12 | **LOW** | **RESOLVED by provenance labeling** | Frozen `S/margin_lemma.py` ends “Certified branch-and-bound ... NOT done.” | It can contradict Campaign I if quoted without chronology. | The audit labels it a frozen pre-campaign review-script summary; only Campaign I results and the two fresh replays determine current acceptance. |
| U-13 | **LOW** | **RESOLVED in manuscript** | Earlier paper wording blurred the lambda sets used by centered and face rules. | Reproduction ambiguity, no soundness effect. | The certificate section now states that centered rules use only \(0,\lambda_{\rm seed}\), while the face rule uses all seven values. |
| U-14 | **LOW** | **RESOLVED in manuscript** | The earlier one-line \(q_2\)-pin description omitted ordering, coefficient-sign, and full-face preconditions. | Essential applicability conditions were hidden. | The final rule inventory states \(q_{2,\rm lo}\ge15/16\), \(p_{2,\rm hi}\le q_{2,\rm lo}\), the positive strip coefficient, derivative sign, and unrestricted full-face verification. |
| U-15 | **LOW** | **RESOLVED to bounded originality wording** | “Largest/first explicit constant with a complete proof” exceeded certificate evidence. | Priority cannot be established from arithmetic. | The dedicated literature audit supports only “apparently first explicit certified improvement” within a dated bounded search; universal priority remains open. |
| U-16 | **LOW** | **RESOLVED interpretation** | Cambie Section 4 says \(Q_i\) is independent of the other parameter variables and then calls \(P_i,Q_i\) dependent. | The two clauses cannot both hold and obscure which self-coupling Q2 uses. | The proof reads the second pair as \(P_i,R_i\), as required by Q2 and the sequential construction in Section 4.1. |
| U-17 | **MEDIUM** | **RESOLVED in proof ledger** | `PROOF.md` wrote \(\|\psi\|^2=\operatorname{rh}\,h\) while defining \(\operatorname{rh}\) as the full diagonal defect. | The notation multiplied by entropy twice and obscured the Margin Lemma. | It now states \(\|\psi\|^2=W(x,x)/\ln2=\operatorname{rh}=\rho h\). |
| U-18 | **MEDIUM** | **RESOLVED in manuscript** | The manuscript folded symmetric couplings to pair-orbit laws but did not state the reverse unfolding before claiming equality of the optimizations. | Only one inequality direction was visible in the submitted proof. | The paper now explicitly splits off-diagonal orbit mass equally between \((p,q)\) and \((q,p)\), keeps diagonal mass diagonal, and verifies both marginals, cost, and inverse folding. |
| U-19 | **LOW** | **RESOLVED in manuscript/script prose** | The scale-free paragraph claimed convergence to \(\kappa_t\) for arbitrary sink collapse. | In general \(C/L\ge0\) need not tend to zero and the limiting mean may be \(m<t\). | With non-sink mass \(\varepsilon\), the proof now gives \(\sigma^2\le\varepsilon\rho_{\max}L\) and \(\liminf\Lambda\ge2(1-\alpha)(1-m)-1\ge\kappa_t\); equality is reserved for the displayed mean-\(t\), zero-cost path. |
| U-20 | **LOW** | **RESOLVED in manuscript/proof script** | The direct extreme-point perturbation proof said “support points” without specifying topological support. | Positive mass of the chosen neighborhoods was implicit. | Both sources now name topological support and invoke its definition before choosing disjoint positive-mass Borel neighborhoods. |

The fresh read-only independent analytic referee found no blocker in the human
chain and produced U-17--U-20.  Its scope, checked transitions, exact findings,
post-repair commands, and residual assumptions are preserved in
[`verification/results/uc-analytic-referee-2026-08-26.md`](verification/results/uc-analytic-referee-2026-08-26.md).

No rejected discrepancy was found to reverse the certified inequality or expose
an uncovered feasible box.  All manuscript/status discrepancies above are
resolved.  The remaining open items are declared trust boundaries and external
mathematical review, not hidden certificate failures.

## 10. Residual trust base

### 10.1 Analytic trust

- Standard finite Shannon entropy facts: chain rule, conditioning reduces
  entropy, and entropy of a variable supported on \(N\) values is at most
  \(\log_2N\). They are used in the human bridge of Section 4.1.
- Riesz representation, Banach–Alaoglu, and Stone–Weierstrass as cited in
  `S/thmB3_proof.py:42-47,112-127,283-303`.
- Bauer's maximum principle as cited in `S/thmB3_proof.py:195-207`.
- Winkler/Pinelis are not required if the direct perturbation proof remains.
- Cambie v2 is corroborative for the bridge construction. Within this audit,
  the bridge is no longer an opaque black-box implication, but it remains a
  human proof rather than a machine proof.

### 10.2 Machine-arithmetic and execution trust

- Correctness of python-flint/Arb outward ball arithmetic and elementary
  functions, including its native library.
- CPython 3.14.3 language, import, `-B`, startup-hook, file-I/O, and exception
  semantics.
- The frozen Python source bytes under externally fixed code/checker/launch
  hashes.
- mpmath/SciPy behavior used to generate deterministic nonnegative lambda
  choices; KKT correctness is not assumed. SymPy correctness for asserted
  identities; NumPy and diagnostics are not universal proof premises.
- SHA-256 collision resistance for integrity binding.
- At campaign commit time, Darwin `renamex_np(RENAME_EXCL)`, `fsync`, and file
  system semantics (`S/cert3_par.py:258-310,565-719`). A fresh mathematical
  replay needs the trace bytes but does not need to trust the original timing
  or search heuristics.
- No claim of protection against a malicious interpreter, malicious same-version
  native library, forged external hash publication, or concurrent hostile file
  system outside these assumptions.

## 11. Explicit final claim envelope

### Established within this audit

1. **HUMAN-AUDITED:** for exact
   \(\alpha_0=356069/10^7\), if the Campaign I relaxed functional is
   nonnegative on its feasible five-parameter two-orbit domain, then
   \(F(\mu)\ge0\) for every regular Borel \(\mu\) on \([0,1]\) with
   \(\mathbb E p\le t_{\rm cert}\). This uses exact decomposition,
   compactness/continuity, fixed-\(B\) concavity, Bauer, the direct support
   perturbation, and \(E\le S^2\).

2. **MACHINE-VERIFIED by an independent arithmetic/control implementation,
   under the named trust base:** all eight exact Campaign I traces were
   replayed from byte zero with empty stacks, no residual event, exact matching
   tallies, and independently rediscovered face proofs.  This establishes
   \(\Phi_{\rm rel}\ge0\) over
   \([0,1]^4\times[1/2,1]\) subject to ordering and \(M\le t_{\rm cert}\);
   orbit-label symmetry covers \(w\in[0,1]\).  The secure verifier imports no
   frozen arithmetic and derives gradients by interval AD, but shares the
   authenticated traces and Arb primitives.  It is bound to verifier hash
   `e67058...91a9a`, runtime-seal raw hash `58c493...2a0c`, trusted report-lock
   raw hash `af8648...e730`, and composite canonical hash `4acbd3...3243`.

3. **HUMAN-AUDITED:** the sequential entropy bridge proves that the functional
   inequality on the closed mean domain implies a present element of every
   finite union-closed \(\mathcal F\ne\{\varnothing\}\) has frequency
   strictly greater than \(t_{\rm cert}\), hence at least
   \(t_{\rm cert}\). The strictness proof includes the
   **MACHINE-VERIFIED** positive Arb endpoint value for the first coordinate.

4. **MACHINE-VERIFIED exact rational comparison:**
   \(t_{\rm cert}\ge(3-\sqrt5)/2+10^{-4}\).

### Not established

- **OPEN:** an end-to-end formal or independently implemented machine proof of
  the analytic chain.
- **OPEN:** any valid certificate at a larger target, including the later
  incomplete campaigns.
- **OPEN:** optimality of \(t_{\rm cert}\), Frankl's \(1/2\), Cambie's reported
  \(0.3823455\ldots\), or Liu's conditional \(0.3827090\ldots\).
- **OPEN in this document:** novelty/priority claims; those require the
  separate literature record.
- **FAILED if stated without qualifications:** “the overall UC implication
  chain is end-to-end machine-checked.” The correct description is a replayed
  Arb certificate for the reduced functional embedded in a human-audited and
  partly cited analytic proof.
