# A machine-checked constant 0.3820660112501052 for union-closed families

**Status — read this first.** **CERTIFIED 2026-08-21.** Campaign I (`cert3_20260818T212601Z_425f109c…_2f23a58ebdb8`) committed all eight slices — 488,465,854 boxes, zero residual, zero stack — and its frozen collector replayed every trace inside its own acceptance gate and printed **`COMPOSITE CERTIFICATE`** (checker SHA-256 `95d09322a7821f41…`; verbatim output at [`campaigns/CERTIFICATE.txt`](campaigns/CERTIFICATE.txt)). The result of [`PROOF.md`](PROOF.md) stands; the campaign-I record below is the authoritative account (see *How to verify*).

## Abstract

Every union-closed family $\mathcal F\neq\{\emptyset\}$ contains an element
belonging to at least $0.3820660112501052\cdot|\mathcal F|$ of its sets,
conditional on a single cited implication (Cambie, arXiv:2212.12500v2,
Question 2 and Section 4), which turns a two-strategy entropy functional
inequality into the union-closed bound. We prove that functional inequality
and machine-check the proof end to end: exact identities decompose the
OR-entropy kernel ($Q=2BL-E$, with $E$ a sum of squared moment defects); a
support-reduction theorem shows the minimum is attained on two pair-orbits
(five parameters); a margin lemma controls the defect by Cauchy–Schwarz;
and an interval branch-and-bound in Arb ball arithmetic certifies
$\Phi\ge0$ at $t=0.3820660112501052\ge\psi+10^{-4}$ (rational-exact: the
certificate's $t-\psi$ strictly exceeds $10^{-4}$, never equals it, since
$t$ is rational and $\psi$ irrational), exceeding
$\psi=\frac{3-\sqrt5}{2}\approx0.3819660112501051$, the largest previously
completely proved constant. Every pruning decision is a certified Arb sign.
All sampled heuristics are labeled NUMERICAL and are never proof steps.
The computation writes one-byte-per-node proof traces; an independent
replayer re-proves every node; a frozen collector re-proves the campaign
before acceptance. The pipeline is deterministic.

## The claim, precisely

Let $h(u)=-u\log_2u-(1-u)\log_2(1-u)$ (binary entropy, $h(0)=h(1)=0$),
$\alpha=0.0356069$ (Cambie's $\alpha$), and
$t_{\mathrm{cert}}=0.3820660112501052\ge\psi+10^{-4}$.

* **Machine-checked here** (`cert3.py` + trace + replay): for
  every probability measure $\mu$ on $[0,1]$ with $\mathbb Ep\le
  t_{\mathrm{cert}}$,
  $$F(\mu)=(1-\alpha)\,Q(\mu)+\alpha\,C(\mu)-L(\mu)\ \ge\ 0,$$
  where $L=\int h\,d\mu$, $Q=\iint h(p+q-pq)\,d\mu d\mu$, and $C$ is the
  symmetric-coupling term with cost $h(s^\star)$, $s^\star(p,r)=
  \max(p,r,\min(p+r,\tfrac12))$.
* **Cited, not reproved** (Cambie, arXiv:2212.12500v2, Q2 + §4): if that
  inequality holds for all $\mu$ with $\mathbb Ep\le c$, then every
  union-closed family $\neq\{\emptyset\}$ has an element with frequency
  $\ge c$.

Together: frequency $\ge 0.3820660112501052$ for every union-closed family
$\neq\{\emptyset\}$. The comparison $t_{\mathrm{cert}}\ge\psi$ is itself
rational-exact: for $r\in[0,3/2]$, $r\ge\psi\iff r^2-3r+1\le 0$
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
$\psi=(3-\sqrt5)/2\approx0.3819660112501051$ — until now the largest
explicit constant with a complete proof. All are statements of the form "a
certain functional of a probability distribution is nonnegative". Going
further requires handling distributions with many atoms. Yu reduced the
problem to five parameters, but that proof has an invalid concavity step —
explicit counterexamples are in [`README.md`](README.md) (Result 3). Cambie
reformulated the remaining task as an entropy functional inequality and
stated the gap himself, verbatim: *"an exact rigorous calculus proof is
missing"* — the remaining problem is a minimisation checkable by computer.

This work closes that gap as a checked computation. The OR-entropy kernel
decomposes exactly as one positive square minus a positive-definite
remainder (Theorem A), giving $Q=2BL-E$ with $E$ a sum of squared moment
defects (Theorem A′). A support-reduction theorem (B‴, tight in orbit
count) shows the minimum of $F$ is attained on two "pair-orbits" — five
parameters $(p_1,q_1,p_2,q_2,w)$. A margin lemma bounds $E$ by
Cauchy–Schwarz in a form that stays finite where the functional degenerates
(at the entropy sinks $0$ and $1$). The rest is a five-dimensional
branch-and-bound in Arb interval arithmetic with guaranteed signs.

Two honesty notes. Anything sampled is labeled NUMERICAL and is never a proof step — the KKT seed $\lambda^\star$ that picks cheap multipliers is numerical, and correctness never depends on it. And the step from the functional inequality to union-closed families is Cambie's cited result, not ours; the certificate covers the functional inequality end to end.

## What is machine-checked

| step | artifact |
|---|---|
| $s^\star$ identity; symmetrize-and-fold (exact bridge to Cambie's form) | `bridge_uc.py` |
| Theorem A: kernel has exactly one positive square (SymPy tautology) | `decomposition.py` |
| Theorem A′: $Q=2BL-E$, $E$ a squared-moment-defect series | `reduction.py` |
| Theorem B‴: attainment on ≤2 pair-orbits (5 parameters) | `thmB3_proof.py` |
| Margin Lemma: $F\ge L\Lambda$, scale-free at the sinks | `margin_lemma.py` |
| $\operatorname{rh}$ bounds; endpoint/derivative lemmas | `lemma_rh_proof.py`, `cert2/3.py` |
| 5-D interval certificate, one-byte-per-node traces | `cert3.py` |
| independent replay of every node | `cert3_replay.py` |
| immutable campaign, atomic commits, hashed snapshot | `cert3_par.py` |
| acceptance only by replaying all 8 traces | `cert3_collect.py` |

Each branch-and-bound rule (mean contractor, corner bound, ratio rule, KKT-shifted center rules, $q_2$-pin face rule, split) is an independently sound lower bound or infeasibility proof in Arb ball arithmetic; skipping a rule can only send a box to the split branch. Campaign parameters (v2 manifests): min box width $10^{-3}$, collar floor $1.25\times10^{-4}$, box budget $2\times10^9$, 24 h/slice, 160-bit covers, 80-bit work.

Tamper tests: the adversarial suite of eight trace attacks — all rejected by the replayer. Determinism: campaigns E and F ran slice 0 five hours apart on different snapshots — 16,231,621 nodes each, byte-identical traces (SHA-256 `9eed0cc1382f5485…`).

## The conditional step

The implication "functional inequality at $c$ ⟹ union-closed element with frequency $\ge c$" is Cambie's (arXiv:2212.12500v2, Question 2 and Section 4) and is used as a citation. His hypothesis is stated for $\mathbb Ep\le c$ ("at most"); the certificate proves the closed domain $\mathbb Ep\le t_{\mathrm{cert}}$, which covers the boundary. Everything downstream of that citation is machine-checked; it is the single non-machine step. The residual trust base (PROOF.md §10): CPython semantics with pinned-source execution (`python -B`), flint/Arb, the pinned library versions, Darwin `renamex_np` (RENAME_EXCL) atomicity, SHA-256, and the published citations (Bauer 1958; Winkler 1988 via Pinelis 2016; van Neerven; Bru–de Siqueira Pedra; Cambie).

## Relation to the literature

| constant | source | status per README.md |
|---|---|---|
| $0.01$ | Gilmer, arXiv:2211.09055 | proved |
| $\psi=\frac{3-\sqrt5}{2}=0.3819660112501051$ | Chase–Lovett, Alweiss–Huang–Sellke, Pebody, Sawin (arXiv:2211.11504v3 Thm 1) | proved; previous largest with a complete proof |
| $c^*=0.3823455333667027$ | Yu (arXiv:2212.00658), Cambie (arXiv:2212.12500) | value likely correct; published proofs invalid (Yu's concavity step refuted by explicit counterexample) |
| $>c^*$, non-explicit | Liu, arXiv:2306.08824 Thm 6 | proved, no explicit value |
| $0.382709087918741$ | Liu, Thm 13 | conditional on two unproved hypotheses |
| $0.3820660112501052$ | this work | machine-checked modulo the cited Cambie implication |

This note claims the largest **explicit** constant with a complete
proof; Liu's non-explicit theorem already exceeds $c^*$.

## How to verify

One command, from the repository's `math/` directory, using the campaign's
own frozen collector (as checked in at [`campaigns/README.md`](campaigns/README.md)):

```sh
cd math
./.venv/bin/python -B uc/campaigns/<ID>/snapshot/cert3_collect.py \
    uc/campaigns/<ID>/launch.json
```

with `<ID> = cert3_20260818T212601Z_425f109c15b64a6198785c6cebbbdaab_2f23a58ebdb8`.

The collector trusts nothing it did not re-derive: hardcoded file
inventories (the manifest's own inventory is not trusted), a pristine
snapshot listing, exact dyadic roots, unique runs, exit code 0, COMPLETE
verdicts, positive work, `stack`=`residual`=`budget_boxes`=`budget_time`=0,
consistent trace digests/sizes — then it **replays all eight traces**
(`cert3_replay.py` re-runs the mean contractor, re-proves each claimed rule
in Arb, follows the recorded splits, requires an empty stack, no residuals,
matching tallies) and only then prints `COMPOSITE CERTIFICATE`. On
2026-08-21 it printed exactly that (verbatim at
[`campaigns/CERTIFICATE.txt`](campaigns/CERTIFICATE.txt)); checker SHA-256
`95d09322a7821f41850ef1ce77345e5e93eb4fa570cc0e7f58f817716cf681ec`.
Environment pinned by `launch.json`: CPython
3.14.3, flint 0.9.0, mpmath 1.3.0, numpy 2.5.2, scipy 1.18.0, sympy
1.14.0, macOS arm64. The eight workers cover the exact dyadic partition of
$w\in[\tfrac12,1]$ (the exact orbit-swap invariance makes $w\ge\tfrac12$
cover the full cube).

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
