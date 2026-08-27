# A certificate-backed candidate constant 0.3820660112501052 for union-closed families

**Finite-certificate status — read this first.** Campaign I
(`cert3_20260818T212601Z_425f109c…_2f23a58ebdb8`) committed all eight
slices: 488,465,854 boxes, zero residual, zero stack. Its frozen collector
replayed every trace and printed **`COMPOSITE CERTIFICATE`** (checker SHA-256
`95d09322a7821f41…`; verbatim output at
[`campaigns/CERTIFICATE.txt`](campaigns/CERTIFICATE.txt)). This machine result
is the five-dimensional interval certificate only. The reduction to it and
the entropy-to-union-closed bridge are ordinary human proofs audited in
[`AUDIT.md`](AUDIT.md); the overall result is a candidate theorem pending
external mathematical review, not an end-to-end machine-checked proof.

**Fresh replay (2026-08-26).** Eight isolated direct replay processes all
exited zero and reproduced the exact expected stdout hashes/tallies; pre/post
input maps matched. Machine-readable summary:
[`verification/results/direct-isolated-replay/summary.json`](verification/results/direct-isolated-replay/summary.json),
canonical SHA-256
`6ef126ced337b035e5c5f22a130b1be3697666e60f16586b75540ef6fc734b16`.
The hardened fresh-copy replay independently passed all eight slices as well;
its composite SHA-256 is
`57ca3f52c054bd9bccfefc349f440674c74fe7ef671d06c103ef3f52dc33d53d`.

## Abstract

We present a certificate-backed, machine-assisted candidate proof that every
nontrivial finite union-closed family has an element belonging to at least
$0.3820660112501052\cdot|\mathcal F|$ of its sets. Exact identities decompose
the OR-entropy kernel ($Q=2BL-E$). A human-audited support-reduction theorem
places the minimum on two pair-orbits, and a human-audited Margin Lemma lowers
the exact functional to a five-parameter surrogate. Arb interval
branch-and-bound certifies that surrogate at the exact rational target
$t=0.3820660112501052\ge\psi+10^{-4}$, where
$\psi=(3-\sqrt5)/2$. The computation writes one byte per DFS node; the replay
re-derives every box and re-proves the named interval rule. The implication
from the functional inequality to union-closed families is reconstructed
directly by a sequential Bernoulli-coupling and entropy-chain argument in the
revised manuscript, following Cambie's Question 2 and Section 4.

The machine-verified layer does not include compactness, support reduction,
the Bochner-triangle Margin Lemma, or the set-family bridge. Sampled controls
and numerical KKT seeds are evidence or heuristics only. No claim of formal
verification, peer review, publication priority, or Frankl's constant $1/2$
is made.

## The claim, precisely

Let $h(u)=-u\log_2u-(1-u)\log_2(1-u)$ (binary entropy, $h(0)=h(1)=0$),
$\alpha=0.0356069$ (Cambie's $\alpha$), and
$t_{\mathrm{cert}}=0.3820660112501052\ge\psi+10^{-4}$.

* **Machine-verified finite claim** (`cert3.py` + traces + replay): the explicit
  lower functional $\widehat\Phi$ is nonnegative on every feasible
  five-parameter two-orbit law at $t_{\mathrm{cert}}$.
* **Human-audited reduction:** exact algebra, compactness, fixed-mean
  concavity, Bauer/extreme-point support reduction, and the Margin Lemma give
  $F(\mu)\ge\widehat\Phi(\nu)$ at a minimizing two-orbit representation. The
  case $L(\mu)=0$ is handled separately and gives equality zero.
* **Human-audited set-family bridge:** the sequential max-entropy Bernoulli
  coupling, entropy chain rule, and a rigorous positive point-mass margin
  convert $F\ge0$ into the union-closed frequency conclusion. Cambie
  arXiv:2212.12500v2 is cited as the source of this construction, not used as
  an unaudited black box.

Together these form a candidate proof of frequency
$\ge0.3820660112501052$. The comparison with $\psi$ is rational-exact: for
$r\in[0,3/2]$, $r\ge\psi\iff r^2-3r+1\le0$
(`cert3_par.target_relation_holds`).

## Plain-English explainer

A family of sets is *union-closed* if the union of any two of its members is
again a member. Frankl's 1979 conjecture says that in any such family (other
than the trivial $\{\emptyset\}$) some element appears in at least half the
sets. It is open.

Gilmer's 2022 idea: choose a set uniformly at random from the family and
study the entropies of the events "element $i$ is present". His
one-strategy argument gives $0.01$; refinements by Chase–Lovett,
Alweiss–Huang–Sellke, Pebody, and Sawin reach
$\psi=(3-\sqrt5)/2\approx0.3819660112501051$, the benchmark with a complete
analytic proof independent of later finite-dimensional numerical claims.
All of these entropy arguments assert that a certain functional of a
probability distribution is nonnegative. Going further requires handling
distributions with many atoms. Yu reduced the problem to five parameters, but
that proof has an invalid concavity step —
explicit counterexamples are in [`README.md`](README.md) (Result 3). Cambie
reformulated the remaining task as an entropy functional inequality and
stated the gap himself, verbatim: *"an exact rigorous calculus proof is
missing"* — the remaining problem is a minimisation checkable by computer.

This work supplies a candidate repair with a reproducible finite certificate.
The OR-entropy kernel decomposes exactly as one positive square minus a
positive-definite remainder (Theorem A), giving $Q=2BL-E$ with $E$ a sum of
squared moment defects (Theorem A′). The human Theorem B‴ places a minimizer
on two pair-orbits. The human Margin Lemma bounds $E$ by the Bochner triangle
inequality in a form that stays finite at the entropy sinks. The interval
certificate then checks the resulting five-dimensional lower functional.

Anything sampled is labeled NUMERICAL and is never a proof step. The KKT seed
$\lambda^\star$ only selects cheap nonnegative multipliers. The certificate
does not machine-check the surrounding human arguments.

## Evidence layers

| step | evidence level | artifact |
|---|---|---|
| $s^\star$ identity; symmetrize-and-fold | human identity + exact finite schemas | `bridge_uc.py` |
| Theorem A kernel decomposition | human series argument + exact SymPy identity | `decomposition.py` |
| Theorem A′: $Q=2BL-E$ | human series argument + exact controls | `reduction.py` |
| Theorem B‴: attainment on ≤2 pair-orbits | **human-audited**, standard cited functional analysis | `thmB3_proof.py`, `AUDIT.md` |
| Margin Lemma: $F\ge L\Lambda$ | **human-audited** Bochner triangle inequality | `margin_lemma.py`, `AUDIT.md` |
| $\operatorname{rh}$ and endpoint bounds | human inequalities + exact symbolic/Arb subchecks | `lemma_rh_proof.py`, `cert2.py`, `cert3.py` |
| 5-D interval certificate and node traces | **machine-verified** Arb arithmetic | `cert3.py`, campaign traces |
| arithmetic replay of every node | **machine-verified**, shares frozen rule implementation | `cert3_replay.py` |
| independent trace topology/tally audit | **machine-verified structural check** | `verification/structural_trace_audit.py` |
| campaign identity and immutable-input pins | **machine-verified external lock** | `verification/campaign-lock.json` |
| entropy-to-union-closed bridge | **human-audited**, with independent Arb strictness check | paper §1, `verification/cambie_bridge_strictness.py` |

Each branch-and-bound rule (mean contractor, corner bound, ratio rule, KKT-shifted center rules, $q_2$-pin face rule, split) is an independently sound lower bound or infeasibility proof in Arb ball arithmetic; skipping a rule can only send a box to the split branch. Campaign parameters (v2 manifests): min box width $10^{-3}$, collar floor $1.25\times10^{-4}$, box budget $2\times10^9$, 24 h/slice, 160-bit covers, 80-bit work.

Tamper tests: the adversarial suite of eight trace attacks — all rejected by the replayer. Determinism: campaigns E and F ran slice 0 five hours apart on different snapshots — 16,231,621 nodes each, byte-identical traces (SHA-256 `9eed0cc1382f5485…`).

## Bridge to union-closed families

Cambie's Question 2 is stated for expectations strictly below $c$, not “at
most $c$.” The interval certificate proves the stronger closed domain
$\mathbb Ep\le t_{\mathrm{cert}}$. The revised paper now gives the bridge
directly: it constructs the dependent uniform sample coordinate by coordinate,
applies the entropy chain rule, and proves strictness at the first present
coordinate. The strict point-mass endpoint is independently enclosed at
$+0.0012927838426659849\ldots$ by 256-bit Arb. This layer remains an ordinary
human proof, not part of the branch-and-bound certificate.

## Relation to the literature

| constant | source | audited status |
|---|---|---|
| $0.01$ | Gilmer, arXiv:2211.09055 | proved |
| $\psi=\frac{3-\sqrt5}{2}$ | Chase–Lovett, AHS, Pebody, Sawin | complete analytic proofs |
| $c^*=0.3823455333667027$ | Yu; Cambie | claimed; Yu's support step is invalid and Cambie's matching lower check remains numerical/graphical |
| $>c^*$, non-explicit | Liu, Theorem 6 | analytic perturbation but imports the unresolved Yu--Cambie optimizer description |
| $0.382709087918741$ | Liu, Theorem 13 | explicitly conditional; local H1 candidate does not resolve H2 |
| $0.3820660112501052$ | this project | candidate theorem: human-audited chain plus replayed finite Arb certificate |

The full bounded search through 2026-08-25 and an arXiv-feed refresh on
2026-08-26 found this apparently the first explicit *certified* improvement
over \(\psi\), not the first claimed improvement. See
[`LITERATURE_ORIGINALITY.md`](LITERATURE_ORIGINALITY.md); universal priority
and external review remain open.

## How to verify

Use [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), not a campaign-local replay
marker. It pins the exact Campaign I launch, fifteen sources, eight results,
and eight traces through [`verification/campaign-lock.json`](verification/campaign-lock.json);
checks the interpreter and dependency versions; supplies isolated `-I -B`
per-slice replay commands; and defines accepted machine-readable reports.

The independent structural command is:

```sh
./.venv/bin/python -I -B uc/verification/structural_trace_audit.py \
  --campaign uc/campaigns/cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8 \
  --output uc/verification/results/structural-traces.json \
  --nice 10 --cpu-limit-seconds 600
```

The structural command checks all 488,465,854 trace bytes, hashes, opcodes, DFS
topology, termination, and tallies without importing certificate arithmetic;
by itself it is not an interval proof.  A separate secure arithmetic
implementation has now replayed every trace from byte zero without importing
the frozen rule modules.  It uses interval AD, rediscovers every face proof,
and is bound to a 343-file runtime seal and live process/image inventory.  Its
trusted report-lock raw SHA-256 is
`af86480901c2c497739916bb410f1e117e4eaf3eefb82575ce494b7d8c04e730`;
the accepted composite canonical SHA-256 is
`4acbd3b935bda7e51ed387e42e0598debf22f4a85b7a24ad976c3eaca4f23243`.

The historical collector printed `COMPOSITE CERTIFICATE` on 2026-08-21, but
does not itself enforce those external provenance/independence gates; the lock,
secure verifier, runtime evidence, and trusted report-lock digest do.

## Campaign I record

Launched 2026-08-18 21:26:01 UTC, 72 h/slice, influence-weighted split
(width × influence; on campaign F's hardest slab this cut 36,511 boxes
to 131 — 279× — and the obstruction root from 231 to 109 boxes). Code
SHA-256 `2f23a58ebdb8…78ce05` (the same code hash as its 24-hour
predecessor H, which committed 5/8 slices before its wall);
launch SHA-256 `6414affcaa54bd5e…df65cc` (full values in its
`launch.json`). The earlier campaigns chart the search: E (5 COMPLETE
slices before the 8 h wall) found nothing mathematical in the way,
F located the true residual barrier (1,447,470 sink residuals on slice 7
— the near-total-sink slab the ratio rule's dead zone left untouched),
G removed it (7/8, slice 7 out of wall clock with residual 0), and I
closed the target with time to spare.

| slice | processed | verdict |
|---|---|---|
| 0 | 21,135,611 | COMPLETE |
| 1 | 20,827,065 | COMPLETE |
| 2 | 23,486,265 | COMPLETE |
| 3 | 34,186,855 | COMPLETE |
| 4 | 71,405,987 | COMPLETE |
| 5 | 117,709,435 | COMPLETE |
| 6 | 121,456,775 | COMPLETE |
| 7 | 78,257,861 | COMPLETE |
| total | **488,465,854** | collector output: **COMPOSITE CERTIFICATE** |

## Where each number comes from

* $\psi$, $t_{\mathrm{cert}}$, $\alpha$, the proof chain §§1–6 and its
  PROVED/CITED/NUMERICAL labels: `PROOF.md`; protocol, collector gates,
  replay, adversarial suite, determinism, trust base: `PROOF.md` §§7–10.
* Frankl's conjecture, the literature ladder, Yu's invalid step, Cambie's
  verbatim caveat: `README.md`.
* Collector one-liner: `campaigns/README.md`; campaign I identifiers and
  pinned environment: its `launch.json`.
