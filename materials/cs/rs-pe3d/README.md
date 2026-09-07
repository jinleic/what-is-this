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

## Current state (agent RsPe3dMechanism, 2026-08-31) — H-GATE mechanism RESOLVED at t=(1,1,1)

**Outcome first.** The H-GATE hypothesis (weight-≤3 window non-empty iff
min(s) ≤ 3) is now a **THEOREM on the t=(1,1,1) slice**, in the stronger
general form: **for each support S, V ∩ F^S = {0} iff |S| < min_i d(C_i)**,
where d(C_i) = s_i − t_i + 1 under hypotheses **H1** (1 ≤ t_i ≤ s_i; if the
paper ever allows t_i = 0 that is the zero code with d = ∞, handled
separately, never as s_i + 1), **H2** (distinct evaluation points — subgroup
elements qualify; *pairwise coprimality of the s_i is irrelevant to this
distance theorem*), **H3** (all Λ-entries nonzero; the diagonal is
weight-preserving). At t = (1,1,1): d(C_i) = s_i, so the min(s) form follows.
FIRST nonempty support weight = min_i d(C_i) exactly. The proof is the
owner's quotient skeleton, **independently re-derived and audited by this
agent** (five links: tensor right-exactness V = ker π_0⊗π_1⊗π_2 = Σ_i
A⊗C_i⊗A; MDS coset independence; a corrected point-isolation assignment step
(T_i taken as a SET — |T_i| ≤ w survives coordinate duplicates); a separator
functional ℓ_i ∈ Q_i^* correctly pulled back as ℓ_i∘π_i (raw coordinate
functionals do NOT descend); and the minimum-line-word sharpness giving the
biconditional both directions). **No novelty claim** — the mechanism is the
owner's; this campaign audited it, closed its two gaps (coordinate repeats,
functional descent), fixed its hypotheses (t_i ≤ s_i), and machine-anchored
it. Frozen at `campaigns/2026-08-31T08-26-54Z_hgateMechanism/`
(`pre_statement.md` committed first, pre-compute, sha256
`8a84ab12…`; `theorem_quotient_proof.md` sha256 `c965bfc0…`;
`controls_part1.py.asrun`; `checksums.sha256`; `manifest.json`).

### MACHINE-VERIFIED controls (exact F_q end to end, semantics audited)

dim V = N − Π(s_i−1) at (13,(2,2,4))→13, (31,(2,3,5))→22, (17,(4,4,4))→37;
kernel identity rank(Σ_i A⊗C_i⊗A) = dim V = N − Π t_i = 37 at (17,(4,4,4));
exhaustive below-d emptiness 16/30/1,830 supports (d = 2,2,3); every axis
line (w = d) carries a codeword; **full w=3 census at (41,(4,4,5)):
82,160 supports, ZERO nonzero** (pre-registered probe P2, confirmatory);
w ≤ 2 at (61,(4,5,5)): 5,050 empty; **s_i = 1 corners: (5,(1,2,2)) 14/14 and
(7,(1,2,3)) 41/41 supports non-empty** — the corner is INSIDE the ≤ arm,
machine-verified; A0 instrument revalidation PN6 3486 = 3486, PP2 917 = 917;
Λ-uniformity: dim V = 13 under three random nonzero Λ draws and Λ = Id at
(13,(2,2,4)); V equals the triple-sum space {x0(a1,a2)+x1(a0,a2)+x2(a0,a1)}
as a SET at (13,(2,2,4)); delta-exact closures: (17,(4,4,4)) line word
δ = 4, (5,(1,2,2)) line word δ = 2 (exhaustive ascending-cost sweeps).

### Rule-5 retractions, both preserved inline in the frozen files

(1) This campaign's own pre-proof derivation "V = ker D0D1D2" is RETRACTED
as an equality of spaces: the machine check at (13,(2,2,4)) found 28
difference-operator rows pairing nontrivially with V-basis vectors, and
dim ker(D0D1D2) = 1 ≠ dim V = 13 there (root cause: each summand x_i
depends on its TWO non-i axes, so D_i annihilates only the one summand
missing axis i; the three kernels intersect in the constants). The
dimension identity and the two correct characterizations (quotient kernel
identity, triple-sum equality) survive. (2) The pre-statement's
pre-registered adjudication "non-empty probe window ⟹ H-GATE min-form
REFUTED" was falsified PRE-COMPUTE by the audited proof and is retracted
inline in `pre_statement.md` §7 (a non-empty window can now only indicate
an instrument bug — halt for bug, not refutation). (3) Process: the E3
census (52.9 s) breached the one-heavy-slot resource contract during
KgBandClose's occupancy — recorded in the theorem file and `manifest.json`
per Main's instruction; a consolidated control replay was cancelled
mid-run on Main's order; all further CPU-bound work by this agent STOPPED;
the pre-registered P1 (43,744 supports) and P3 (166,750) full censuses
remain HELD and are unnecessary to the theorem.

### Rule 7 — the swept set exactly, and what was not swept

Swept: t = (1,1,1) with Λ = Id unless stated (L1 adds three random nonzero-Λ
draws at (13,(2,2,4))); exact-F_q checks on the ten instances named in the
manifest; completeness limited to all supports below d on the E1 instances,
the full w = 3 census on (41,(4,4,5)), w ≤ 2 on (61,(4,5,5)), all supports
w ≤ 3 on the two s_i = 1 corners and the two A0 anchors; delta closed
exactly (exhaustive sweeps) only at the D1/D2 words.
**Not swept, plainly:** any support of weight ≥ 4 anywhere (ρ^window stays
UNDEFINED — not zero — at min(s) ≥ 4 instances by the theorem; the theorem
makes no ratio claim); the pre-registered P1/P3 full censuses (HELD);
every t ≠ (1,1,1) — the frozen q=13 t-profile rows are only
COROLLARY-CONSISTENT with the general threshold min_i(s_i−t_i+1), an
INFERENCE from quoted README rows, not re-derived; t_i = 0 (out of scope,
d = ∞ convention); nonprime fields; Conjecture 4.2's ρ(r,η,β) content, all
η, and every global ρ_inst — untouched, as is the weight-≥4 window
everywhere. The theorem is HUMAN-AUDITED derivation MACHINE-ANCHORED by the
controls; it is not a machine enumeration beyond the completed slices
listed. Evidence labels: controls MACHINE-VERIFIED (exact F_q); theorem
HUMAN-AUDITED (independently re-derived; owner skeleton credited); all
frozen prior rows REPORTED (never re-run, per non-goals).

**Reproduce:** from this directory,
`PYTHONPATH=$PWD ../.venv/bin/python campaigns/2026-08-31T08-26-54Z_hgateMechanism/controls_part1.py.asrun`
(≈5 min, CPU-bound — schedule in a free heavy slot; the session transcript
holds the per-block as-run outputs (E3 52.9 s, corners <1 s, A0 70.6 s);
`sha256sum -c campaigns/2026-08-31T08-26-54Z_hgateMechanism/checksums.sha256`.
Next campaign (priced in `manifest.json`): **H-WINDOW4** — at (17,(4,4,4))
the theorem puts the FIRST nonempty supports exactly at weight 4; the open
question is whether the weight-≥4 window's optimal ratio drops below 1
(line words give wt/δ = 1). First slice w=4 only: C(64,4) = 635,376
supports ≈ 409 s at the measured 1,553 supports/s; w ≤ 6 ≈ 15.4 h.
Cheap add-on: upgrade the q=13 t-profile corollary INFERENCE to
machine-check (minutes).

## Owner correction — 2026-08-31 (RsPe3dMechanism campaign, per Main)

The owner review REJECTS five textual claims in the frozen section above and
in the frozen campaign files (`README` section
"H-GATE mechanism RESOLVED", `theorem_quotient_proof.md`, `manifest.json`).
**The theorem survives; the wording did not.** The frozen campaign is NOT
modified and NOT re-checksummed — cite it only THROUGH this correction.
The five corrections, originals preserved verbatim in the frozen files:

1. **Individual-support biconditional is FALSE.** Original frozen text said
   "V ∩ F^S = {0} iff |S| < min_i d(C_i)" for each individual S. Correct
   theorem: for EVERY S with |S| < d, V ∩ F^S = {0}; and there EXISTS a line
   support S of size d with V ∩ F^S ≠ {0}. Equivalently: the UNION window
   over supports |S| ≤ w is empty iff w < d. **No claim is made for each
   individual S with |S| ≥ d.** (This README's earlier sentence "FIRST
   nonempty support weight = min_i d(C_i) exactly" must be read in the
   exists-threshold sense for the union window, not as a per-S claim.)
2. **Proof link 2 (coset independence), wording fix:** an arbitrary relation
   Σ c_a π_i(e_a) = 0 yields c = Σ c_a e_a ∈ C_i with supp(c) ⊆ T — the
   frozen text's "⟺ 1_T ∈ C_i" is false in general (the relation vector c
   need not be the all-ones indicator; only its support is controlled).
   The independence conclusion stands (wt(c) ≤ |T| < d(C_i)).
3. **Minimum GRS word, wording fix:** it is the degree-(t_i − 1) polynomial
   vanishing at t_i − 1 evaluation points (not at s_i − 1). At t_i = 1 both
   readings coincide with the constant/one-point form, which is why the
   t = (1,1,1) slice's line-indicator statement survives verbatim. The
   minimum-distance FACT d(C_i) = s_i − t_i + 1 (MDS) is unaffected.
4. **Kernel dimension, arithmetic fix:** dim ker(π_0⊗π_1⊗π_2) =
   N − Π_i (s_i − t_i), not "N − Π_i t_i" as one frozen line says; at
   (17,(4,4,4)): 37 = 64 − 27 = 64 − (s_i − t_i)³. The machine check K1's
   measured value 37 was correct; the formula printed next to it was not.
5. **Triple-sum representation is NOT unique.** Only equality of spaces/sets
   holds (V = {x0(a1,a2) + x1(a0,a2) + x2(a0,a1)} as sets, machine-verified
   T1); the frozen sentence "every v ∈ V has a UNIQUE representation" (in
   `pre_statement.md` §1) is false as stated.

These errors affect WORDING, not the corrected minimum-distance proof: the
right-exactness kernel identity, MDS distance, union-window emptiness ⟺
w < d, the existence threshold at d, and all MACHINE-VERIFIED control
numbers stand. Campaign inventory pointer:
`campaigns/2026-08-31T08-26-54Z_hgateMechanism/` (frozen at git `dab17e7`,
read its `checksums.sha256`/`manifest.json` only through this correction;
biconditional statements anywhere in the frozen trio are superseded by
correction 1's two-part form). No compute was run for this correction.

**Consequence for the next-campaign note above (H-WINDOW4):** its first
sentence must also be read through correction 1 — at (17,(4,4,4)) the UNION
window is empty for w ≤ 3 and NON-empty at w = 4 (line words exist); whether
the w = 4 window's optimal ratio drops below 1 remains the open question,
unchanged in substance.

**Correction 6 (preregistration provenance, no compute).** The frozen
section's claims that `pre_statement.md` was "committed as the FIRST file,
pre-compute" are FALSE as git provenance: `git log --reverse` on the
campaign directory shows ONLY the freeze commit `dab17e7`, which adds
`pre_statement.md` together with the theorem, script, manifest, and
checksums in one commit. The file may have been written before compute per
session chronology, but no separate commit proves preregistration.
Evidence grade for every "pre-registered/pre-committed" claim in the frozen
section and campaign files: **session chronology / REPORTED** — not
commit-anchored. The frozen campaign remains preserved and unmodified;
`dab17e7` is cited only as the freeze commit. The two inline retractions
inside `pre_statement.md` remain part of its recorded content; their
TIMING relative to compute is now REPORTED-grade, not PROVEN by commits.

**Correction 7 (tuple typo in frozen pre_statement, no compute).** Frozen
`pre_statement.md` line 161 names the A0 anchor instance "PP2 (43,(2,3,4))";
the actual frozen PP2 is **(43,(2,3,7))** (as in the prior P2
pre-statement/run and as used by `controls_part1.py.asrun` line 72 and the
theorem file's control table). Tuple typo only — the 917 nonzero-support
replay was run at the CORRECT (43,(2,3,7)), so the frozen count 917 = 917
stands unchanged. Frozen file untouched; frozen campaign preserved.

---

## Campaign W (2026-09-01T20-30-00Z_win4_pn4) — CLOSED, verdict FROZEN-CERTIFIED (agent RsPe3dW2, 2026-09-01)

Weight-EXACTLY-4 window at PN4 (421,(4,5,7),t=(1,1,1),Λ=Id), complete parity-PRIMARY
census of all C(140,4) = 15,329,615 supports (exact F_421; the non-candidate side is
settled by rank = |S| directly, so no missed-candidate mode exists): **exactly 35
candidates, SET-EQUAL to the 35 direction-0 minimum lines** {o,o+35,o+70,o+105},
o ∈ [0,34] — 0 non-line supports, 0 missing, all dim(V∩F^S) = 1, per-hit lex indices
re-verified. Per-candidate rho: delta = 4 exactly for all 35 (single direction-0 line
closure each), so **ρ^window_w4 = 1** exactly (Fraction 4/4). Evidence labels:
census + set-equality + rho MACHINE-VERIFIED (exact enumeration, zero float);
H-GATE status at this instance: its ≤3-window content is the M theorem's corollary
(PN4 w≤3 empty re-anchored, 457,450/457,450), the w=4 nonemptiness matches Theorem T
part 2 (existence at |S| = d = 4 — NOT part 1, which is emptiness strictly below d),
and the observed 35-peak characterization is consistent with the owner's separate
H-MINLINE conjecture (recorded as INPUT-consistent; not folded into W). Pattern P
makes NO committed prediction at this window (its scope is ≤3); ρ^window_w4 = 1 is
the expected-floor retrospective observation; H-GATE is NOT falsified by anything
here and its single prospective read (PN10, P2-frozen) stands untouched. V-basis
ladder (w_ladder): B01-B10 exact at freeze; B11-B13 still running (frozen bytes
untouched; post-freeze rows land in census_ledger.jsonl per the addendum's named
next action). Two instrument defects in the NEW w=4 parity code (ledger set-content
crash; idx-counter line lost in a fixing edit) are disclosed with preserved evidence
ledgers (DEFECT*/DEFECT2*); neither touched any frozen byte, the M campaign, or the
mathematical formulation — the final run is defect-free end-to-end. Grandfathered
legacy campaign per the 2026-09-01 workspace control plane (no migration; frozen
pre-statement 8114a8f untouched; all additions additive + checksums_addendum_w2.sha256).
Full detail, option-(b) record, and rule-7 scope sentence: campaigns/2026-09-01T20-30-00Z_win4_pn4/freeze_addendum_w2.md.

---

## H-SPARK2-CONVERSE — CLOSED, FROZEN-CERTIFIED (2026-09-04)

Run `20260904T040116Z_00df70c5_9bf664de2def` closes the proportional-column
boundary left out of the frozen $d_A,d_B\ge3$ crossing converse. For
$d_A=2$, the size-$d_B$, profile-$(2,d_B)$ circuits are exactly the
nontrivial two-row partitions of a $B$-circuit across a proportional
$A$-column pair, with exact count

$$C_A(2)C_B(d_B)(2^{d_B}-2).$$

The circuit coefficients have the same sign on both rows; if
$a_{u'}=\lambda a_u$, the second-row coefficient is
$\tau\lambda^{-1}c_J(v)$ and recovering $\tau$ multiplies by $\lambda$.
Crossing-only holds at $d_B=2$ (two diagonals, with the 2-to-1 crossing
parameter coincidence) and $d_B=3$; the noncrossing remainder
$2^{d_B}-2-2d_B$ is positive exactly for $d_B\ge4$. The result is symmetric
when the spark-2 factor is $B$.

The final complete rerun passed 120 assertions: all 15 row/column pairs,
9,300 registered partition candidates, 536 relation-pattern checks, exact
rank-before-deletion promotion, full-rank noncircuit rejects, and 4,500
below-threshold independence checks. Crucially, the separately coded
$(4,2)$ loop used a scale-2 B duplicate at GF(13) and GF(11); all 70 circuits
per prime recovered one common $\tau$ by multiplying with $\lambda=2$, while
all 70 deliberately inverted recoveries were rejected. Full theorem,
provenance, failed-launch disclosure, audit, and checksums are in the run
directory. No frozen predecessor or `cs/RESULTS.md` was edited.

Immediate successor `20260904T051205Z_d03471d8_0604b708af38`
(H-C3-34-RESIDUAL-POINT) is RUNNING. Its preregistration and init are
commit-bound. Reviewed launch 2 passed 26 pre-geometry controls (including a
one-million-distinct-pair stream control) and completed the registered
one-million-pair GF(17) search with outcome E-empty: no qualifying reduced
all-distinct 12-point base set occurred in that fixed factorized family.
That finite emptiness is not a universal C3 conclusion; the analytic
Bezout residual-length-one lemma is a separate obligation, and a nonempty
construction family requires a bound amendment before any further search.

Existence language remains conditional: for the all-distinct
$K=2+\sum_i(r_i-1)$ law, “none below $K$” is unconditional, while the
simultaneous-PGL classification at $K$ supplies a circuit only when the
relevant incidence population is nonzero. This README does not claim spark
exactly $K$ in empty instances.

## H-C3-34 residual closeout correction (Main, 2026-09-04)

The POINT successor above is no longer RUNNING: it closed
**FROZEN-INCONCLUSIVE** with an E-empty one-million-pair sample. Its
residual-length-one lemma is separately analytically verified under the
stated proper-intersection and local-length-one hypotheses; the rational
residual point may lie outside the affine grid, and existence/minimality
remain separate obligations.

The registered
[FULLPOOL run](campaigns/20260904T152816Z_bb611cd1_455c4fdc03e3/VERDICT.md)
is now producer-closed **FROZEN-CERTIFIED: finite-pool E-empty**. All
**41,587,200** compatible pairs were exhausted, with 49/49 assertions passing.
All **185** size-12 intersections repeat a projection coordinate, so there
are **zero qualifying base sets**. A fresh repaired driver-free replay
matched the complete pool, histogram, candidate list and decisions; all
nine comparisons and four GCD controls passed. All 45 frozen checksums
verified. Failed launches, the stopped memory-defective replay, and the
misnamed attempt-5 artifact copy are preserved and disclosed.

This closes only the exact registered finite family, not universal C3,
characteristic-zero, or extension-field existence. The next gate needs a
nonempty construction family or an analytic existence/minimality argument.
The successful replay obeyed the low-impact policy: sampled host CPU at
most 39.25%, 20.2 MiB group RSS, 28 KiB new output, at least 128.388 GiB free.
[Root results](../RESULTS.md) hold the authoritative outcome and scope.

## Explicit characteristic-two residual circuits — CLOSED (Main, 2026-09-04)

The next nonempty construction is now certified:
[H-C3-34-CHAR2-EXPLICIT](campaigns/20260904T224223Z_9c13dba4_cbfbca8c94c3/VERDICT.md),
producer verdict **FROZEN-CERTIFIED**. For the twelve roots of
$T^{12}+T^3+T^2+1$ over **GF(2048)**, the points $(t^3,t^2)$ form the
reduced proper $(2,3)$ complete intersection
$x^2+y^3=x^2y^3+x+y+1=0$ with injective projections.

Full evaluation rank is **10**; all **66** ten-point subsets are independent
and all **12** eleven-point deletions are circuits. For omitted root $a$,
the coefficients are $(t+a)^2/t^4$ on the remaining points. The residual
point is the omitted point, possibly outside the retained coordinate grid.
Independent analytic review, 86 exact primary checks, 106 independent
bit-arithmetic checks and all 20 frozen checksums pass.

The final primary/replay computations together used 0.099366 CPU seconds
and 28 KiB new output; neither enumerated the field or full product grid.
The warning-bearing first attempt is preserved beside the source repair and
clean full rerun. This proves characteristic-two existence/minimality for
these supports, not global spark eleven, classification, or other
characteristics. A prospective extension to other characteristics requires
new registration and evidence; the earlier GF(17) finite-pool miss is intact.


<!-- cs-autonomy:20260905T045017-d18d5052195f:1 -->
## Autonomous rs-pe3d closeout — 2026-09-05 05:15 UTC

**FROZEN-INCONCLUSIVE**. REVIEW of H-C3-35-ALLCHAR-EXPLICIT: INCONCLUSIVE (no compute receipt; the worker wall budget was exhausted before the sealed verifier could run, so no verified/rejected verdict is issued). What was completed: (1) Full analytic proof audit by reading the registered theorem_allchar_explicit.md (sha 6dab404f...): every step checked by hand and found sound within scope — disc(P mod p) = D mod p for monic P (valid also at p | 12 since lc = 1); x-injectivity from g + x^2 f = x^4 - x - y - 1; y-injectivity and a+b != 0 from P(T) - P(-T) = -2T^3 (char 2 handled); B = {f=g=0} with no points at infinity (F, G checked) and Jacobian det = 12t^13 - 3t^4 - 2t^3 = t^2 P'(t) recomputed by hand, so all 12 points reduced and exhaust the length-12 intersection; C has unit rows e_0,e_2..e_10 and zero columns 1, 11 so rank 10 in every characteristic with left kernel span(f,g); ten-subset independence from [T]Q_ab = -(a+b)/(a^2 b^2) != 0; relation from [T^11]L_t = 1/P'(t), [T]L_t = 1/(t^2 P'(t)); converse for p in E. Only a cosmetic note: 'degree dividing 12!' can be sharpened to 11! since P(-1) = 0. Primary receipts compute-01/02 (rc 0, sealed_intact) and the 50/50 PASS stdout log were read directly. (2) An independent pure-Python verifier review_verify.py (sha ff477886..., no sympy/flint, primary driver not imported) was written, preregistered and sealed in review/evidence together with review_statement.md: R1 Sylvester/Bareiss discriminant and Res(P(T),P(-T)); R2 own mod-p gcd + distinct-degree factorisation for all 47 primes compared to the primary table; R3 direct Z[a,b] expansion of both ideal-membership witnesses; R4 C-matrix rank over Q and mod all probe primes, left kernel, tau; R5 fresh exact replay at a different non-exception prime p = 29 in GF(29^10) = GF(29)[u]/(h) (h = degree-10 factor of P, roots as Frobenius conjugates): 66 ten-subsets, 12 eleven-subsets, explicit lambda_a with full support, corrupted relation; R6 positive failure certificate at p = 13 and 10417770367 (exact double root, hence <= 11 points); R7 data comparison with core_allchar.json (sha 5a13c42c...). It has not been executed; the claim remains a candidate pending that run. Nothing is promoted.

Evidence: [`rs-pe3d/campaigns/20260905T050431Z_28f76694_6cffbbf64083`](campaigns/20260905T050431Z_28f76694_6cffbbf64083/autonomy/finalize.json). The label applies only to the registered claim; no broader frontier improvement is implied.

Next registered-work proposal: Re-run this review with a fresh wall budget: read only job_info, copy the already-written pure-Python verifier from scratch/autonomy/20260905T050609-677fcc187da7/1/review/evidence/review_verify.py and its statement verbatim, preregister within ~60 s, and execute it once (argv: review_out.json, path to the primary core_allchar.json; expected < 60 s, < 8 MiB). Finish verified iff verdict REVIEW-PASS on all R1-R7 checks including the GF(29^10) replay and the p = 13 failure certificate.


<!-- cs-autonomy:20260905T045017-73185087234c:1 -->
## Autonomous rs-pe3d closeout — 2026-09-05 05:28 UTC

**FROZEN-INCONCLUSIVE**. REVIEW of H-C3-36-FAMILY-PRIMEFIELD: INCONCLUSIVE by tool rule only (the sealed verifier exited rc 1 because of a reviewer-side artefact; no substantive discrepancy with the primary was found, and the tool requires an rc-0 receipt for 'verified'). Compute-01 (sealed_intact, 23.3 s, 12 KiB growth) ran review_verify.py (sha ce23987d..., pure Python: own mod-p elimination, own coset DP, own brute force; primary driver not imported; primary JSON read as data). Substantive results, all PASS: R2 brute force over all 372,736 admissible 12-subsets of GF(29)^* yields exactly 448 solutions in 16 scaling orbits, containing the three primary-reported sets and the explicit instance R={2,3,4,7,11,12,19,20,21,23,24,28}; R3 ALL 448 GF(29) solutions pass every exact check (c_1=c_11=0, c_0!=0, no +- pair, injective squares/cubes, 12x12 M over bidegree-(2,3) monomials rank 10, f,g in left kernel and independent, Jacobian 2x g_y+3y^2 g_x = t^2 P'(t) != 0, g(t^3,t^2)=P(t) for all t in GF(29), 66 ten-subsets rank 10, 12 eleven-subsets rank 10 with lambda_a=(t^2-a^2)/(t^2P'(t)) annihilating, full support, corrupted lambda rejected, Q_ab/Q_t coefficient formulas, residual off-grid); recomputed P=[16,0,13,7,13,8,19,21,22,13,23,0,1] matches. R4 the three reported solutions at each of p=37,41,43,47 pass all checks with matching P. R1 own DP reproduces N_29=448, N_37=108, N_41=305,840, N_43=3,360, N_47=2,507,506 and admissible counts 372,736/46,656/515,973,120/6,205,248/5,538,111,488, and the exact set of primes with max admissible <12 = {2,3,5,7,11,13,17,19,23,31} with maxima 1,1,2,2,5,4,8,6,11,10 equal to the primary table. R5 negative control: 8 random 12-subsets of GF(29)^* with sum r = sum 1/r = 0 containing a +- pair {a,-a}: the ten-subset omitting {a,-a} has rank exactly 9 and all other ten-subsets rank 10, confirming 'independent iff no +- pair'. The only failed flag, R1_dp_table_matches_primary, is because my DP caps |S| at 12 (line 'if k + dk > 12: continue'), so max_admissible saturates at 12 for p in {29,41,43,47} where the primary reports 14,20,14,23; the uncapped values equal n_cosets times per-coset max (1 for p=2 mod 3, 2 for p=1 mod 3) and were checked by hand only. Analytic audit by hand, sound within scope: [1]Q_ab=c_0/(ab), [T]Q_ab=(a+b)c_0/(a^2b^2) from c_1=0; Jacobian by chain rule; C has unit rows e_0,e_2..e_10 and zero columns 1,11 iff c_1=c_11=0 (rank 10, kernel span(f,g)); relation via sum_t t^k/P'(t)=[T^11](T^k mod P)=0 for k in {-2,0,2,..,12} with T^-2 = -(T^10+c_10T^8+..+c_2)/c_0 mod P; no points at infinity; mu_6-coset lemma (exponent differences +-1 in Z/6 => <=2 per coset when p=1 mod 3) giving p>=29 / p>=37; c_1=-c_0 sum 1/r, c_11=-sum r. Minor wording: 'e=11 not representable' holds within the bidegree-(2,3) box only (x y^4 is outside it), which is what is meant. Reviewer's verdict on substance: the registered claim and scope are supported (p=29 smallest prime field for THIS family; no claim on extension fields q=25,27, other supports, C3 or spark eleven); promotion should wait for the rc-0 rerun proposed below.

Evidence: [`rs-pe3d/campaigns/20260905T052113Z_c8edc589_cabb1348014b`](campaigns/20260905T052113Z_c8edc589_cabb1348014b/autonomy/finalize.json). The label applies only to the registered claim; no broader frontier improvement is implied.

Next registered-work proposal: Re-run the review verifier scratch/autonomy/20260905T052236-17e9b035bbae/1/review/evidence/review_verify.py with one edit: in coset_dp remove the '|S| <= 12' cap (or cap at 24) so max_admissible is uncapped, keep everything else verbatim, preregister within 60 s and execute once (argv: review_out.json, path to family_primefield.json; expected ~25 s, <16 MiB). Finish verified iff verdict REVIEW-PASS with rc 0.
