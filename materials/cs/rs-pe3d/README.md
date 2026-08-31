# `rs-pe3d/` — high-dimensional product expansion for Reed–Solomon tensor codes

**Status: BENCHMARK (no campaign run yet).** Lowest duplication risk in the
whole scan: the conjecture is ten days old, and the only existing attempt at it
is an AI-generated proof that **the authors themselves state they have not
verified**.

## The claim

Owner-read (**V-owner**, 2026-08-29) — ECCC **TR26-150**, Mitali Bafna & Nikhil
Vyas, *"Private PCPs from Product Expansion"* (2026-08-19). Abstract, verbatim
on the two load-bearing sentences:

> A proof obtained by ChatGPT 5.6 Sol establishes constant 2-dimensional
> product expansion whenever both rates are bounded away from one, crossing the
> tight sum-of-rates-below-one barrier in the product expansion theorem of
> Polishchuk and Spielman. We conjecture the analogous statement in higher
> dimensions.

Context from the same abstract: their private-PCP construction gets
$\sqrt n$-query private PCPs of size $O(n\log n)$ unconditionally; **conditional
on the high-dimensional product expansion conjecture for Reed–Solomon codes**
the PCPs have length $n^{1+o(1)}$, $n^{o(1)}$ queries, and fractional privacy
and soundness gap $n^{-o(1)}$. The codes are tensor products of Reed–Solomon
codes whose *"key innovation is to choose the evaluation domains as
multiplicative subgroups of pairwise coprime orders of $\mathbb F_q^\star$."*

Scout-extracted from the full text (**V-scout**, owner read pending): the
formal statement is **Conjecture 4.2**, carrying an explicit **balance
parameter $\beta$** and expansion $\rho(r,\eta,\beta)$; the proved
2-dimensional **Theorem 2.1/4.11** gives $\rho=\Omega(\varepsilon^6/K^3)$ but
under a $|S_1|\approx|S_2|$ balance hypothesis **which the informal
higher-dimensional conjecture drops**; the AI candidate proof is of the form
$\rho=\varepsilon^{\exp(k)}$ and is unverified.

## Why this is the right shape of target

The conjecture is about a **finite, explicitly constructible object**: for fixed
$r$, prime power $q$, multiplicative subgroups $S_1,\dots,S_r\le\mathbb F_q^\star$
of pairwise coprime orders, and Reed–Solomon degree profile, the expansion
constant $\rho$ is the optimum of a finite program over an explicitly presented
$\mathbb F_q$-linear space. Exact $\mathbb F_q$ construction plus exact-rational
LP/ILP gives $\rho$ *exactly*, with a rational dual certificate. Nobody has
computed a single exact $\rho(3,\eta,\beta)$ value. Two possible outcomes, both
valuable: the first exact 3-dimensional data points, or a small-parameter
violation that refutes the conjectured form — and since the near-linear PCP
length depends on it, a refutation is load-bearing.

The balance-hypothesis gap is the sharpest wedge: the proved 2D theorem needs
$|S_1|\approx|S_2|$; the conjecture as stated for higher $r$ drops it.
**Deliberately unbalanced small instances are therefore the highest-value
probe** — they test exactly the step the AI proof glosses.

## Gates

1. **Gate A (first exact 3D expansion constants, days).**
   Fix $r=3$. **Correction, 2026-08-30, from agent `RsPe3d`:** this README
   originally said "choose $q$ with $12\mid q-1$" and listed
   $q\in\{13,25,37,\dots\}$. That condition was an imprecise proxy and is
   replaced by the actual requirement — **$q-1$ must admit three pairwise
   coprime divisors $>1$**, since the construction needs multiplicative
   subgroups $S_1,S_2,S_3\le\mathbb F_q^\star$ of pairwise coprime orders.
   ($q=31$ satisfies this, with $q-1=30=2\cdot3\cdot5$ giving orders
   $(2,3,5)$, even though $12\nmid30$.) Also noted: $q=25$ is a prime power,
   so $\mathbb F_{25}$ is a perfectly valid field, but $|\mathbb F_{25}^\star|=24=2^3\cdot3$
   admits no three pairwise coprime divisors $>1$, so it is correctly dropped.
   Adopted census: $q\in\{31,61,109,181,241,421,601\}$, all coprime subgroup
   triples among the divisors of $q-1$, rate slack
   $\eta\in\{1/128,1/64,1/32\}$, and **both balanced and deliberately
   unbalanced** $(|S_1|,|S_2|,|S_3|)$ profiles. For each instance, construct the
   tensor-code spaces in exact $\mathbb F_q$ arithmetic and compute the
   expansion ratio exactly by LP/ILP (HiGHS with exact-rational verification of
   the returned dual).
   **Deliverable:** the first table of exact $\rho(3,\eta,\beta)$ values.
   **Refutation:** an instance whose exact $\rho$ is below any $\rho$ the
   conjectured form permits at that $(\eta,\beta)$ — exhibited with its dual
   certificate. **Pass:** the conjecture survives every instance in the sweep,
   with the measured $\rho$ trend reported against $\beta$.
2. **Gate B (isolate the balance hypothesis, days).**
   Hold $q,\eta$ fixed and vary $\beta$ monotonically from balanced to maximally
   unbalanced, measuring exact $\rho(\beta)$. **Deliverable:** whether $\rho$
   degrades gracefully (conjecture plausible without balance) or collapses at a
   threshold $\beta^\*$ (conjecture needs the hypothesis the 2D theorem has and
   the conjecture omits — a precise, publishable correction to the conjecture's
   statement rather than to its truth).
3. **Gate C (audit the AI proof, days).**
   Transcribe the unverified AI candidate proof's claimed $\rho=\varepsilon^{\exp(k)}$
   dependence, instantiate it at the gate-A parameters, and compare against the
   exact measured values. **Refutation:** measured $\rho$ inconsistent with the
   candidate's form — the first evidence on an explicitly unverified AI proof.
   **Pass:** consistent, which is evidence for, never proof of, the conjecture.

## Pre-campaign requirements

- Owner first-hand read of the full TR26-150: the formal Conjecture 4.2 with
  its quantifiers, Definition 4.1, Theorem 2.1/4.11 and the exact role of the
  balance hypothesis, and the precise statement of the AI candidate proof's
  claim. Quote all of it verbatim into the pre-statement — the gate is
  meaningless if the conjecture is paraphrased.
- Confirm the expansion definition used (which norm/weight, which quotient) so
  a "refutation" cannot be dismissed as a definitional mismatch. This is the
  single largest risk on this target and must be closed before any run.
- One-page pre-statement per gate: instance enumeration, the exact LP, the dual
  verification procedure, and the numeric refutation threshold.

## Disjointness

TR26-150 is motivated by quantum PCPs, but every object here is **classical**:
Reed–Solomon tensor codes over $\mathbb F_q$, classical expansion, classical LP.
No quantum code, decoder, or circuit appears, so there is no overlap with
`../physics/` or `../math/qec/`. If the work drifts toward quantum codes'
distance or decoding, it must stop and be re-scoped.

## Layout

- `src/`, `campaigns/`, `scratch/` — created at first use.

## Current state (agent RsPe3d, 2026-08-29)

- **Primary read:** ECCC TR26-150, Bafna--Vyas, *Private PCPs from Product
  Expansion*, including Definition 4.1, Theorem 2.1/4.11, and Conjecture 4.2.
  The formal expansion quantity and the AI claim transcription are fixed in
  `pre_statement.md`.
- **Gate A — COMPUTATIONAL-EVIDENCE, finite window only:** at
  $(q,s,\eta,t)=(31,(2,3,5),1/32,(1,1,1))$, every point support of weight
  at most $3$ was exhaustively tested by exact $\mathbb F_{31}$ row reduction.
  The normalized reduced-basis candidate class has 55 vectors; all 55
  HiGHS MILPs are `kOptimal`, and each returned decomposition passed exact
  sum/weight checks plus exhaustive exact line-support checks below its
  returned cost. The minimum ratio in this candidate class is
  $\mathrm{wt}/\delta=3/5$ (witness weight $3$, $\delta=5$). This is an
  upper bound for the global $\rho_{\rm inst}$, not a global exact value:
  heavier supports and non-basis values remain **OPEN**. The same complete
  weight-$\le3$ support census, without optimization, is recorded for all
  eight feasible $t$ profiles at this pilot.
- **Gate B — OPEN:** no fixed-$q,\eta$ balance sweep over alternative subgroup
  triples has yet been closed; imbalance dependence remains unanswered.
- **Gate C — COMPUTATIONAL-EVIDENCE:** under the pre-registered
  $\varepsilon=\eta,\ k=3,\ \exp(k)=e^3$ mapping, the proved universal bound
  $\rho_{\rm inst}\ge 1/(3N)=1/90$ and the exact inequality
  $\eta^{e^3}<\eta^{20}<1/90$ for $\eta\in\{1/32,1/64,1/128\}$ give the
  verdict “consistent-with-form on proved interval.” The global exact
  $\rho_{\rm inst}$ is **OPEN**, so this is not an exact-value comparison.
- **Unit findings:** the exhaustive tiny unit
  $(q,s,t)=(5,(2,2,2),(1,1,1))$ has exact $\rho=1/3$ and 5^7−1 nonzero
  vectors; see `campaigns-smoke/unit-anchors.json`. A line-sum equality I
  derived as a unit premise fails on the feasible substitution
  $(13,(2,3,4),(2,3,4))$: an exact split has cost 23 versus the
  one-direction cost 24. This is a scratch/unit premise failure, not a
  TR26-150 claim; see `scratch/anchor1_counterexample.json`.
- **Frozen artifacts:** Gate-A pilot
  `campaigns/2026-08-30T00-46-32Z_d7452695/`; Gate-C audit
  `campaigns/2026-08-30T00-49-26Z_98cae3e5_gateC/`.
- **Reproduction:** from this directory, run
  `PYTHONPATH=$PWD ../.venv/bin/python -m src.campaign` (the repository
  interpreter path is `/Users/jinleic/jinleic-workspace/cs/.venv/bin/python`);
  verify a frozen witness file with
  `python -m src.verify campaigns/<run>/witnesses.json`.

## Current state (agent RsPe3d, 2026-08-30) — exact finite-window table

Exact per-witness values now frozen at
`campaigns/2026-08-30T01-04-49Z_exacttable/` (includes
`witnesses.json`, `payload.json`, and `optimality.json`).
Window, identical for every row: **weight $\le3$ support census over all
point sets, candidates = normalized reduced basis of each nonzero
$V\cap F^S$; per-witness $\delta$ proven exact by exhaustive line-tuple
sweep (all lower costs absent, $\ge1$ at cost).**

| q | s | t | ratio | weight | delta | label |
|---|---|---|-------|--------|-------|-------|
| 13 | (2,2,4) | (1,1,1) | 1/2 | 2 | 4 | MACHINE-VERIFIED (per-witness finite) |
| 13 | (2,2,4) | (1,1,2) | 3/8 | 3 | 8 | MACHINE-VERIFIED (per-witness finite) |
| 13 | (2,2,4) | (1,1,4) | 1/4 | 1 | 4 | MACHINE-VERIFIED (per-witness finite) |
| 31 | (2,3,5) | (1,1,1) | 3/5 | 3 | 5 | MACHINE-VERIFIED (per-witness finite) |
| 31 | (2,3,10) | (1,1,1) | 3/5 | 3 | 5 | MACHINE-VERIFIED (per-witness finite) |
| 31 | (2,6,5) | (1,1,1) | 1 | 2 | 2 | MACHINE-VERIFIED (per-witness finite) |
| 31 | (3,10,2) | (1,1,1) | 3/5 | 3 | 5 | MACHINE-VERIFIED (per-witness finite) |
| 61 | (2,3,5) | (1,1,1) | 3/5 | 3 | 5 | MACHINE-VERIFIED (per-witness finite) |
| 61 | (3,4,5) | (1,1,1) | 1 | 3 | 3 | MACHINE-VERIFIED (per-witness finite) |

Scope: ratios are **best within the stated window**, not global minima;
global rho over heavier supports stays OPEN (same addendum policy).
Gate-B signal at q=61: (2,3,5) gives 3/5 while more balanced
(3,4,5) collapses to 1 — imbalance does NOT collapse rho here; it
improves (3/5 vs 1). No threshold beta* found in this small sweep.

## Current state (agent RsPe3dGateB, 2026-08-30) — Gate B CLOSED, verdict "neither"

**Outcome first.** Gate B's dichotomy does not apply. On the swept set the
window ratio $\rho^{\rm window}$ is **not a function of $\beta$ at all**, so
there is neither graceful degradation nor a collapse threshold $\beta^\ast$ —
the pre-registered **third outcome**. 13 of 13 instances landed as exact
rationals, every one MACHINE-VERIFIED in-window. Frozen at
`campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/`
(pre-statement, `payload.json`, `witnesses.json`, `optimality.json`,
`manifest.json`, `checksums.sha256`, scripts as run).

### The sweep (all rows: $\eta=1/32$, $\Lambda=\mathrm{Id}$, exact rationals)

| id | q | s | t | $\beta=\frac{\max s_i}{\min s_i}$ | $\rho^{\rm window}$ | wt | $\delta$ | coprime? |
|----|---|---|---|------|------|----|----|----|
| B13 | 241 | (3,4,5) | (1,1,1) | 5/3 | 1 | 3 | 3 | yes |
| B02 | 61 | (3,4,5) | (1,1,1) | 5/3 | 1 | 3 | 3 | yes |
| B07 | 13 | (2,2,4) | (1,1,1) | 2 | 1/2 | 2 | 4 | no |
| B08 | 13 | (2,2,4) | (1,1,2) | 2 | 3/8 | 3 | 8 | no |
| B09 | 13 | (2,2,4) | (1,1,4) | 2 | 1/4 | 1 | 4 | no |
| B01 | 61 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | 3 | 5 | yes |
| B03 | 31 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | 3 | 5 | yes |
| B12 | 241 | (2,3,5) | (1,1,1) | 5/2 | 3/5 | 3 | 5 | yes |
| B10 | 241 | (3,5,8) | (1,1,1) | 8/3 | 1 | 3 | 3 | yes |
| B05 | 31 | (2,6,5) | (1,1,1) | 3 | 1 | 2 | 2 | no |
| B04 | 31 | (2,3,10) | (1,1,1) | 5 | 3/5 | 3 | 5 | no |
| B06 | 31 | (3,10,2) | (1,1,1) | 5 | 3/5 | 3 | 5 | no |
| B11 | 241 | (3,5,16) | (1,1,1) | 16/3 | 1 | 3 | 3 | yes |

Every row: complete weight-$\le3$ point-support census, candidates = normalized
reduced basis of each nonzero $V\cap\mathbb F^S$, and $\delta$ closed **exactly**
per witness (all cheaper line-tuples proved absent, $\ge1$ present at the cost).
B11's census enumerated all **2,304,200** supports (3,449 s). Twelve rows were
additionally re-run in a fresh process: **zero ratio disagreements**; B11 is
single-pass with its witness independently verified and closure-proved. Rows
marked "coprime? no" violate Conjecture 4.2's pairwise-coprime hypothesis and
are shape probes inherited from gate A, not tests of the conjecture.

### Exclusion E — a CONCLUSION, independent of any hypothesis

$\rho^{\rm window}$ is not monotone — indeed not single-valued — in $\beta$. The
exact rows interleave: at $q=31$, $\beta=5/2\to3/5$, $\beta=3\to1$,
$\beta=5\to3/5$ (a peak); at $q=241$, $\beta=5/3\to1$, $\beta=5/2\to3/5$,
$\beta=8/3\to1$ (a valley). B11 seals it at the extreme: $\beta$ grows $3.2\times$
from $5/3$ to $16/3$ with $\rho$ going $1\to1$, while the *intermediate*
$\beta=5/2$ sits at $3/5$. So no monotone $\rho(\beta)$ fits, no $\beta^\ast$
exists to locate, and the grid-resolution uncertainty statement applies
vacuously (swept $\beta$: $5/3, 2, 5/2, 8/3, 3, 5, 16/3$).

### Pattern P — EVIDENCE, not a theorem

What does track the data: $\{2,3\}\subseteq\{s_1,s_2,s_3\} \iff
\rho^{\rm window}=3/5$, exact on all 13 rows (5 contain both and all are $3/5$;
8 do not and none is). **12 rows are in-sample** (P was found by inspecting
them); the in-sample separation figure $1/\binom{13}{5}=7.8\times10^{-4}$ is an
**IN-SAMPLE** number and is not a prospective p-value. **One row, B11, is
out-of-sample**: P predicted $1$, a $\beta$-degradation story predicted
$<1/4$ (the grid floor); observed $1$ — the grid *ceiling*. The prediction was
committed mid-census with a machine-checkable order certificate
(`prediction_order_proof.json`: pre-statement sha256, census 1,350,000/2,304,200,
`b11_result_file_exists=false`). **P survived its first prospective test; P is
not established.** It dies to any triple containing $\{2,3\}$ with
$\rho\ne3/5$, or any triple lacking the pair with $\rho=3/5$. Known limits: P
speaks only about *orders* at $t=(1,1,1)$; $\rho$ demonstrably varies with $t$ at
fixed orders (the three $q=13$ rows give $1/2, 3/8, 1/4$), so P cannot be the
whole story.

### Finding 2 — where the balance hypothesis actually lives (provenance)

First-hand read of TR26-150 (PDF sha256 `fb51d7dc8486…`) and of the linked AI
candidate proof (`e2d9b4b61f9e…`, recovered from the PDF's `/URI` annotation).
Balance is **present** in Conjecture 4.2 itself ("balance parameter $\beta\ge1$
… $\beta^{-1}\le s_i/s_j\le\beta$"), in Theorem 4.11 ($K^{-1}\le n_1/n_2\le K$,
$\rho=\Omega(\varepsilon^6/K^3)$), in Lemma 5.4 item 2 ($1/4\le s_i/s_j\le4$,
"4-balanced", used as $\rho(k,1/128,4)$), in §6's own summary of its hypotheses,
**and in the unverified AI proof** (its "$K$-comparable" Proposition 5;
$\rho=\frac{1}{2dK^d}\big(\frac{\varepsilon}{(C_0dK)^{10d}}\big)^{d/\gamma_d}$).
It is absent **only** from the informal Conjecture 1.4 ("$\rho$ being only a
function of $\varepsilon$ and $k$"). So the correction is a wording correction to
one informal restatement; **Conjecture 4.2 needs no correction on this axis**,
and the gate's founding premise — that the 3-D conjecture omits the hypothesis
the 2-D theorem has — is false as applied to Conjecture 4.2. Full quotes with
line numbers: `campaigns/…_gateB/finding2_balance_hypothesis_provenance.md`.
That file also carries a **rule-5 inline retraction**: I first reported the AI
form's balance cost as $K^{-4863}$ by guessing the ambiguous ASCII exponent as
$3^d6^{d-2}$; the proof's own recurrence $\gamma_1=1$,
$\gamma_e=1/(3\cdot6^{e-2})$ pins the exponent to $d/\gamma_d = 3d\cdot6^{d-2}$
($6$ at $d{=}2$, $54$ at $d{=}3$ — the $3^d$ reading gives $9$ and $162$ and
fails at both), so the correct figure is $\rho\sim K^{-1623}$ at $d=3$
[DERIVED, exact `Fraction` arithmetic]. Still an enormous negative power of the
balance parameter, which is finding 2's point.

### Rule 7 — the swept set exactly, and what was not swept

Swept: $\eta=1/32$; $\Lambda=\mathrm{Id}$; $t=(1,1,1)$ except the two extra
$q=13$ rows at $t=(1,1,2),(1,1,4)$; the 13 $(q,s)$ pairs tabulated above
($q\in\{13,31,61,241\}$); witnesses restricted to weight-$\le3$ supports and the
normalized reduced-basis candidate class.

**Not swept, plainly:** supports of weight $\ge4$; non-basis field values inside
$\dim>1$ intersections; other $t$; other $\eta$; **non-identity $\Lambda$** (the
conjecture's $\forall\Lambda_i$ quantifier is untouched); primes outside
$\{13,31,61,241\}$; the $q=421$ nineteen-triple census — **priced at 36 h**
census-only at the measured 650 supports/s, which is why it stays OPEN;
continuous $\beta$ (only the discrete grid); any asymptotic claim in $\beta$, $q$
or $N$. No global $\rho_{\rm inst}$ is claimed: each row certifies only the
exact interval $\rho_{\rm inst}\in[1/(3N),\ \rho^{\rm window}]$, e.g. B11
$[1/720,\ 1]$ of exact width $719/720$, B01 $[1/90,\ 3/5]$ of exact width
$53/90$ — both endpoints exact rationals, zero-width arithmetic throughout (no
float anywhere in the chain).

**Anchors before extension:** unit $(5,(2,2,2),(1,1,1))\to\rho=1/3$ exact
(78,124 nonzero $V$ vectors) reproduced by the sweep code **and** by a
structurally independent second route (V-basis coefficient enumeration +
ascending-cost closure, 536 s, agreeing at wt 2 / $\delta$ 6); all nine frozen
gate-A ratios reproduced exactly (3/5, 1, 3/5, 3/5, 1, 3/5, 1/2, 3/8, 1/4) with
zero mismatches.

**Reproduce:** `cd cs/rs-pe3d && PYTHONPATH=$PWD ../.venv/bin/python -c "from
src.gateb import run_instance; print(run_instance({'id':'B01','q':61,'s':(2,3,5),'t':(1,1,1),'expect':'3/5'})[0])"`;
B11 via the frozen `run_b11_fast.py.asrun` (3,449 s); verify every witness with
`python -m src.verify campaigns/2026-08-30T11-43-59Z_50595CC9_fffe84b0_gateB/witnesses.json`
(passes 13/13 with exhaustive lower-support closure); re-freeze with
`python -m src.gateb_freeze`.

## Current state (agent RsPe3dPatternP, 2026-08-30) — P2 battery: P survives
a second pre-registered battery; the owner coprimality hypothesis FALSIFIED;
window-emptiness gate found and prospectively confirmed

**Outcome first.** The second out-of-sample battery for pattern P landed at
`campaigns/2026-08-30T18-04-32Z_99AFA09B_patternP2/` (pre-statement committed
as the directory's FIRST file, git `bf541f0`, with per-instance P-predictions
fixed before computation; order proof `order_proof.json`). **No P miss
occurred.** Broken out by arm, per Main's interpretive correction, because
the arms are not equally informative:

- **PP arm (contains $\{2,3\}$; P predicts exactly $3/5$) — the
  discriminating arm: 9/9 hits.** PP1–PP9 at new primes
  $q\in\{13,43,97,109,67,157,601,193\}$ ($s=(2,3,k)$, $k\in$
  $\{4,7,8,9,11,13,16,25,32\}$, $\beta$ up to 16) all landed
  $\rho^{\rm window}=3/5$ exactly (wt=3, $\delta=5$, per-witness closure
  proven). Illustrative only, labeled as such: a uniform one-of-five null on
  this arm is $(1/5)^9=1/1953125\approx5.12\times10^{-7}$; this is not a
  rigorously defined null and is not a p-value.
- **PN arm (lacks $\{2,3\}$; P predicts only "not $3/5$") — the weak arm:
  7/7 landed hits, 3 untestable.** PN1 $(41,(2,4,5))$,
  PN2 $(71,(2,5,7))$, PN3 $(211,(3,5,7))$, PN6 $(43,(2,6,7))$,
  PN7 $(337,(3,4,7))$ gave $1/1$ and PN8 $(11,(2,2,5))$,
  PN9 $(17,(2,2,8))$ gave $1/2$ — all exact, wt/delta =
  $2/2,2/2,3/3,2/2,3/3,2/4,2/4$ with closure proven (PN8/PN9 completed in
  freeze addendum 2 after the first runner's crash at PN4; the battery's
  only genuinely new sub-1 value is PN8/PN9's $1/2$). PN4, PN5, PN10 have
  EMPTY weight-$\le3$ windows (below), so P makes no testable prediction
  there.
- **Status: P survived its second battery and is still not established.**
  Out-of-sample record: B11 + 16 landing rows here; the thin direction (PN)
  remains under-tested and is where a biconditional usually dies.

### Owner checks on the frozen table — independently confirmed (MACHINE-VERIFIED)

- **OC1:** $\rho=\mathrm{wt}/\delta$ in every row — all 13 frozen witnesses
  re-verified in exact $\mathbb F_q$ (`src.verify.verify_record`), and every
  recorded ratio equals `Fraction(weight, delta)` recomputed from its own
  witness. Zero mismatches.
- **OC2:** P is not a shadow of `pairwise_coprime` — confirmed within both
  classes: pairwise-coprime rows (B01,B02,B03,B10,B11,B12,B13): the three
  $(2,3,5)$ give $3/5$, the four others give $1$; non-coprime rows
  (B04,B05,B06,B07,B08,B09): $(2,3,10),(3,10,2)$ give $3/5$, the rest do
  not. (Note the class assignment is PAIRWISE: B04–B06 are mutually but not
  pairwise coprime — gcd chain 1 while a pair shares a factor.)
- **OC4 (description):** the owner's observation is exact as a description of
  the frozen 13: coprime rows take $\{1,3/5\}$ only; non-coprime rows take
  five values $\{1/4,3/8,1/2,3/5,1\}$.

### H-COP — the owner observation pre-registered as a hypothesis and FALSIFIED

Pre-registered in P2-4 with its own falsifiable content, adjudicated
separately from P: H-COP.a says coprime rows land in $\{1,3/5\}$ with no-pair
coprime rows exactly $1$; H-COP.b says non-coprime rows show $\ge3$ distinct
sub-1 values. **Verdict: H-COP.a holds on landed rows; H-COP.b FAILS — the landed
non-coprime rows take only $3/5$ (PP1, PP3, PP4, PP7, PP9), $1$ (PN1, PN2,
PN6) and $1/2$ (PN8, PN9; added in freeze addendum 2): TWO distinct sub-1
values, not the predicted $\ge3$.** The five-value richness did not reproduce on new rows.
Main's ledger: this falsifies the hypothesis Main reframed from Main's own
OC4 observation (owner structural record now 0-for-13). Caveat kept: the
larger non-coprime instances are window-blind, so the falsification is
restricted to the landed set.

### H-GATE — the window-emptiness gate (Main's min(s) hypothesis), prospectively confirmed

Finding promoted per Main: **at every instance ever run on this target, the
weight-$\le3$ window is non-empty iff $\min_i s_i\le3$.** All 13 frozen rows
($\min s\in\{2,3\}$) have non-empty windows; PN4 $(421,(4,5,7))$, PN5
$(61,(4,5,6))$ and PN10 $(281,(5,7,8))$ ($\min s=4,4,5$) have windows with
**zero** nonzero $V\cap\mathbb F^S$ supports, by complete censuses of
$457{,}450=140+9{,}730+447{,}580$, $288{,}100$ and $3{,}658{,}900$ supports
respectively (vectorized kernel over the $V$-basis; instrument validated by
exact agreement with the frozen census code on non-empty PN6 3486=3486 and
PP2 917=917; random sampling of $V$ at PN4 puts its minimum weight around
$134$ of $140$). Main's prediction was committed mid-census
(`addendum_HGATE_prediction.md` + `hgate_order_proof.json`, PN10 at 49%,
zero hits so far) and the prospective read came back EMPTY as predicted.
The confound is broken: product magnitude does not separate the classes
(frozen $(3,5,16)$ has product 240 and a non-empty window; PN4 has product
140 and an empty one) — $\min(s)$ does. Status: EVIDENCE on 17 consistent
battery+frozen rows with one prospective confirmation; the mechanism
(dual product-code weight structure) is OPEN; corner cases $s_i=1$ and
equal-order triples are untested.
**Significance (rule-15 named domain):** if H-GATE holds, the gate-B
programme's effective domain is exactly $\{$instances with an order
$\le3\}$, and every ratio this target has reported — the frozen 13 included
— is conditional on that predicate. Under H-GATE a min($s$)$\ge4$ instance
cannot be tested in this window at all; testing P at larger orders requires
a heavier window (support weight $\ge4$), which is OPEN everywhere.

### Rule 7 — the swept set exactly, and what was not swept

Swept: $t=(1,1,1)$; $\Lambda=\mathrm{Id}$; the 19 pre-registered $(q,s)$
pairs at $q\in\{11,13,17,41,43,61,67,71,97,109,157,193,211,281,337,421,601\}$
with $(s_1,s_2,s_3)$ exactly as frozen in the pre-statement's P2-2 table;
complete weight-$\le3$ point-support censuses per instance with the
normalized reduced-basis candidate class; per-witness $\delta$ exactly
closed at each landing best witness (all cheaper line-tuples absent,
$\ge1$ present); empty-window rows certified only as emptiness facts.

**Not swept, plainly:** supports of weight $\ge4$ (and the H-GATE-hardened
suspicion that they are where min($s$)$\ge4$ instances live); non-basis
field values inside $\dim>1$ intersections; other $t$; other $\eta$;
non-identity $\Lambda$ (the conjecture's $\forall\Lambda_i$ quantifier is
untouched); every $(q,s)$ pair outside the 19-row list; the frozen 13 rows
themselves (NOT re-run — non-goal); continuous $\beta$; any asymptotic claim
in $q$, $N$, $\beta$; and the global $\rho_{\rm inst}$ everywhere — every
number in this section is $\rho^{\rm window}$ (window best), never
$\rho_{\rm inst}$; each landing row certifies only the exact interval
$\rho_{\rm inst}\in[1/(3N),\rho^{\rm window}]$ (e.g. PP1 $[1/72,3/5]$,
PN3 $[1/315,1]$). Nothing here touches TR26-150's statements; no
escalation arose.

Evidence labels: all 16 landing rows MACHINE-VERIFIED (finite window; exact
$\mathbb F_q$ witness verification, regenerable via
`python -m src.verify campaigns/2026-08-30T18-04-32Z_99AFA09B_patternP2/witnesses.json`,
passes 14/14 core witnesses.json + 2 in witnesses_addendum2.json); the 3 empty-window rows MACHINE-VERIFIED as emptiness
censuses; no float anywhere in the chain. Inheritance disclosure: the
battery reuses the frozen gate-B instrument (census→MILP→closure code,
sha256 unchanged $a9ff7e0a\ldots$) and re-derives OC1/OC2/OC4 from the
frozen witnesses rather than payload strings; PN6/PN7 used a vectorized
census kernel validated against the frozen census before use.

**Reproduce:** `cd cs/rs-pe3d && PYTHONPATH=$PWD ../.venv/bin/python
campaigns/2026-08-30T18-04-32Z_99AFA09B_patternP2/p2_run.py.asrun` (scripts
are byte-copies; restore to /tmp or run in place with PYTHONPATH set);
single-row re-check via the same `src.gateb.run_instance` recipe with the
P2 instance; freeze artifacts and sha256 in
`campaigns/2026-08-30T18-04-32Z_99AFA09B_patternP2/{manifest.json,checksums.sha256}`.
