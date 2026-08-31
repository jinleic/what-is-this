# PRE-STATEMENT — P2 battery — rs-pe3d — second out-of-sample test battery for pattern P

Committed: 2026-08-30T18:04Z (UTC), by agent `RsPe3dPatternP`, BEFORE any
battery computation ran. This file was written into the fresh campaign
directory `campaigns/2026-08-30T18-04-32Z_99AFA09B_patternP2/` as its FIRST
file; when this directory was created it contained ONLY this file and nothing
was written to `scratch/`. Code hash at commitment time
(`src.campaign._code_hash()`):
`a9ff7e0ae1a716640c249f63baf4cd8ba774d6e7763822a7ce32d8fca6fde810`
(unchanged from the frozen gate-B code — the battery reuses the exact frozen
census/optimization/closure machinery).

Any later edit to this file must be a dated addendum; the original text stays
visible. Every prediction below is fixed now. Any instance added after
computation starts is a separate pre-registered statement reported alongside
(rule 16: no domain-shopping; the searched set is this list, reported with
every outcome including misses).

---

## P2-0. Where pattern P stands (inherited state, verified against frozen files)

Pattern P (gate-B GB9.1): at $t=(1,1,1)$, $\Lambda=\mathrm{Id}$, in the
weight-$\le3$ support census over all point sets with the normalized
reduced-basis candidate class of $V\cap\mathbb F^S$,

$$ \{2,3\}\subseteq\{s_1,s_2,s_3\}\;\Longleftrightarrow\;\rho^{\rm window}=3/5. $$

It held on all 13 frozen rows; 12 are IN-SAMPLE (P was found by inspecting
them), and exactly one row, B11 = $(241,(3,5,16),(1,1,1))$, was out-of-sample
(pre-registered GB9 with `prediction_order_proof.json`) and PASSED: orders
contain 3 but not 2, observed $\rho^{\rm window}=1$, not $3/5$.
The in-sample separation statistic $1/\binom{13}{5}=1/1287=7.77\times10^{-4}$
is IN-SAMPLE ONLY, quoted with that label, never as a prospective p-value.
P rests on exactly one out-of-sample test.

Rho conventions: **all numbers in this battery are $\rho^{\rm window}$** —
the min of $\mathrm{wt}(M)/\delta(M)$ over the stated weight-$\le3$
reduced-basis candidate class — never global $\rho_{\rm inst}$, which stays
OPEN (heavier supports $\ge4$, non-basis values, other $t$, other $\eta$,
non-identity $\Lambda$ all excluded). Every row additionally certifies only
the exact interval $\rho_{\rm inst}\in[1/(3N),\,\rho^{\rm window}]$.

## P2-1. Owner checks to confirm independently (before any new instance)

Three checks on the frozen 13-row table, to be re-derived from frozen data,
not assumed (rule 17: owner-supplied figures are inputs, not evidence):

1. **(OC1)** $\rho = \mathrm{wt}/\delta$ in every row — i.e. the recorded
   ratio string equals the exact `Fraction(weight, delta)` of its own
   witness, with the witness verified in exact $\mathbb F_q$ by
   `src.verify.verify_record`.
2. **(OC2)** P is not a shadow of `pairwise_coprime`: it holds separately
   within each coprimality class. Among the 7 pairwise-coprime rows, exactly
   the three $(2,3,5)$ rows give $3/5$ and $(3,4,5)\times2,(3,5,8),(3,5,16)$
   give not-$3/5$; among the 6 non-coprime rows, $(2,3,10),(3,10,2)$ give
   $3/5$ and $(2,2,4)\times3,(2,6,5)$ do not.
3. **(OC3)** Verification is confirmation of arithmetic only; P's in-sample
   status is unchanged by OC1/OC2 (recorded here so it cannot be upgraded).

**(OC4, owner observation, explicitly NOT a finding — pre-registered here as
a FALSIFIABLE hypothesis.)** On the 13 frozen points the 7 coprime rows take
only two distinct $\rho^{\rm window}$ values $\{1,3/5\}$ while the 6
non-coprime rows take five, $\{1/4,3/8,1/2,3/5,1\}$. Pronounced NOT a
finding; treated here as a candidate hypothesis to kill (owner structural
hypotheses are 0-for-9). Its second, independent prediction is fixed in P2-4
below. First confirm OC4's description exactly from the frozen rows.

## P2-2. The battery design (fixed before any computation)

19 NEW instances at NEW primes or NEW triples — none among the frozen 13
$(q,s,t)$ rows. All at $t=(1,1,1)$, $\Lambda=\mathrm{Id}$, $\eta$- slack
trivially satisfied ($t_i=1\le s_i$). Every triple verified NOW, pre-run, to
satisfy $s_i\mid q-1$, $q$ prime, $s_i\ge2$. Two families deliberately
weighted toward BREAKING P:

- **PP arm ("contains {2,3}"; P predicts $3/5$):** 9 instances. If ANY lands
  at $\rho^{\rm window}\ne3/5$, P is falsified in the forward direction.
- **PN arm ("lacks {2,3}"; P predicts $\rho^{\rm window}\ne3/5$):** 10
  instances, chosen at orders near the $3/5$-row family (small primes,
  including $(2,4,5),(2,5,7),(3,4,7),(4,5,7)$ — orders close to but unequal
  the winning $(2,3,5)$). If ANY lands at exactly $3/5$, P is falsified in
  the reverse direction.

Additional pre-registered design properties, committed now: every row is a
MINIMAL discriminator by construction — the arms differ from the frozen rows
in exactly one structural knob (new third order with $2,3$ kept; dropped
$3$ with small orders kept; non-coprime triple WITH $\{2,3\}$; larger $q$).

The instance list (P-predictions are FIXED below, BEFORE any computation):

| id | q | s | t | beta(s) | coprime? | P-prediction $\rho^{\rm window}$ |
|----|---|---|---|---------|----------|-----------------------------------|
| PP1 | 13 | (2,3,4) | (1,1,1) | 2 | no  | **3/5** |
| PP2 | 43 | (2,3,7) | (1,1,1) | 7/2 | yes | **3/5** |
| PP3 | 97 | (2,3,8) | (1,1,1) | 4 | no  | **3/5** |
| PP4 | 109 | (2,3,9) | (1,1,1) | 9/2 | no  | **3/5** |
| PP5 | 67 | (2,3,11) | (1,1,1) | 11/2 | yes | **3/5** |
| PP6 | 157 | (2,3,13) | (1,1,1) | 13/2 | yes | **3/5** |
| PP7 | 97 | (2,3,16) | (1,1,1) | 8 | no  | **3/5** |
| PP8 | 601 | (2,3,25) | (1,1,1) | 25/2 | yes | **3/5** |
| PP9 | 193 | (2,3,32) | (1,1,1) | 16 | no  | **3/5** |
| PN1 | 41 | (2,4,5) | (1,1,1) | 5/2 | no  | **not 3/5** |
| PN2 | 71 | (2,5,7) | (1,1,1) | 7/2 | no  | **not 3/5** |
| PN3 | 211 | (3,5,7) | (1,1,1) | 7/3 | yes | **not 3/5** |
| PN4 | 421 | (4,5,7) | (1,1,1) | 7/4 | yes | **not 3/5** |
| PN5 | 61 | (4,5,6) | (1,1,1) | 3/2 | no  | **not 3/5** |
| PN6 | 43 | (2,6,7) | (1,1,1) | 7/2 | no  | **not 3/5** |
| PN7 | 337 | (3,4,7) | (1,1,1) | 7/3 | yes | **not 3/5** |
| PN8 | 11 | (2,2,5) | (1,1,1) | 5/2 | no  | **not 3/5** |
| PN9 | 17 | (2,2,8) | (1,1,1) | 4 | no  | **not 3/5** |
| PN10 | 281 | (5,7,8) | (1,1,1) | 8/5 | yes | **not 3/5** |

Notes committed now:
- PN8/PN9 use REPEATED orders $(2,2,\cdot)$ — same shape as frozen rows
  B07/B08/B09 $(13,(2,2,4))$; deliberately different primes. Repeated-order
  triples cannot be pairwise coprime; they are shape probes, never tests of
  Conjecture 4.2 itself.
- PP8 ($q=601$) is the largest PRIME in the parent pre-statement's in-scope
  census list $\{31,61,109,181,241,421,601\}$ not yet swept at $t=(1,1,1)$;
  $q=601$ satisfies the corrected instance condition (three pairwise coprime
  divisors of $600$: $2,3,25$ — note $25\mid600$ and $\gcd(3,25)=1$,
  $\gcd(2,25)=1$).
- PP9 ($N=192$, census $1{,}179{,}808$ supports) and PP6/PP8 are the
  deliberate "larger q" stress rows; PN10 ($N=280$, census $3{,}658{,}900
  \approx$ B11's $2{,}304{,}200$) pushes the same scale outside the coprime
  $\{2,3\}$-free hypothesis space.
- battery census totals at 650 supports/s measured throughput: arm PP
  ~2.1 M supports ($\sim$54 min census), arm PN ~4.8 M ($\sim$124 min), plus
  per-candidate MILPs and per-witness closure.

## P2-3. Adjudication for P (fixed now)

- **P-falsification event (headline):** ANY PP row with
  $\rho^{\rm window}\ne3/5$, or ANY PN row with $\rho^{\rm window}=3/5$.
  Reported FIRST and prominently, with its exact witness certificate. P then
  has exactly the falsifying row(s) attached to its name in every summary.
- **P-survival event:** all 19 rows match their P2-2 predictions. Then P has
  TWO out-of-sample tests passed (this battery + B11) and 0 fails; status
  upgrades to "P survived two independent preregistrations; still not a
  theorem" — never stronger language than that.
- **Mixed:** some PP fail or PN succeed but not both patterns — report the
  exact rows; the hypothesis is falsified (one counterexample ends it, per
  gate-B's own what-would-kill-P list), regardless of how many rows agree.

The window is FIXED: weight-$\le3$ point supports, complete census;
candidates = normalized reduced basis of each nonzero $V\cap\mathbb F^S$;
best candidate tuned by exact MILP (HiGHS, 120 s cap) with every accepted
decomposition re-verified in exact $\mathbb F_q$; per-witness $\delta$ closed
exactly by the frozen exhaustive cheaper-line-tuple sweep (every cheaper cost
proved absent, $\ge1$ tuple present at the claimed cost) — the identical
standard as the frozen 13. All arithmetic exact rationals; no float anywhere
in the certified chain (the python-flint interval landmines and the HiGHS
1e-7 tolerance trap are inherited as standing cautions; none of this
pipeline touches $q-1$ scale ratios near tolerance, and no interval
arithmetic is constructed in this campaign at all — every quantity is an
exact integer or Fraction).

Sanity: any computed ratio $>1$ or $=0$ is a bug — STOP, not recorded.

## P2-4. Second independent prediction: the coprimality ratio-set separation (OC4)

The owner's observation is reframed here as the following FALSIFIABLE
statement — its first prospective test, with the instance set FIXED (the
same 19 rows):

$$\textbf{H-COP: } \text{At }t=(1,1,1),\ \Lambda=\mathrm{Id}\text{, in the same window:}
\ \ |\{\rho^{\rm window}\}|\ \text{over the pairwise-coprime instances}\ \le\ 2
\ \text{ values in } \{1,3/5\},\ \text{while non-coprime instances take}\
\ge3\ \text{ distinct values below } 1.$$

More precisely, the committed predictions per row-class on the battery:

- **H-COP.a (coprime arm:** PP2, PP5, PP6, PP8, PN3, PN4, PN7, PN10 — 8
  coprime rows): each row's $\rho^{\rm window}\in\{1,3/5\}$. Prediction:
  every one of the 8 lands in $\{1,3/5\}$, colliding with values P already
  assigns: coprime rows containing $\{2,3\}$ must be $3/5$ (overlaps
  P2-2); **coprime rows NOT containing $\{2,3\}$ must be exactly $1$** (the
  frozen coprime no-pair rows are all exactly $1$: B02, B10, B11, B13).
- **H-COP.b (non-coprime arm:** PP1, PP3, PP4, PP7, PP9, PN1, PN2, PN5,
  PN6, PN8, PN9 — 11 rows): the row set as a whole takes **at least 3
  distinct** $\rho^{\rm window}$ values strictly below $1$ (the frozen
  non-coprime rows took five: $1/4,3/8,1/2,3/5,1$; the observation is that
  coprimality, not the $\{2,3\}$ pair, drives multi-valuedness below 1).
  Note H-COP.b is a SET-level prediction: it is FALSIFIED iff all 11
  non-coprime rows land on fewer than 3 distinct sub-1 values.
- H-COP.a and P make OPPOSITE-ish demands only jointly with the frozen
  data: PN3/PN7/PN4/PN10 (coprime, no pair) must be $1$ under H-COP.a but
  merely "not $3/5$" under P; any of them landing at exactly $3/5$ kills
  BOTH P and H-COP simultaneously; any landing at e.g. $2/5$ or $4/5$
  kills H-COP.a while P survives. This discriminates the two hypotheses —
  the point of pre-registering both.
- Adjudication (own section, separate outcome): H-COP confirmed iff BOTH
  H-COP.a and H-COP.b hold on the 19 rows as fixed above; falsified
  otherwise; the two sub-branches reported separately either way. H-COP's
  outcome NEVER upgrades P (and vice versa): two labels, two verdicts.

## P2-5. Rule-7 statement (the swept set, fixed now; restated verbatim in the report)

Swept: $t=(1,1,1)$; $\Lambda=\mathrm{Id}$; $\eta$ unconstrained-and-moot
($t_i=1$); the 19 listed $(q,s)$ pairs with $q\in\{11,13,17,41,43,61,67,71,$
$97,109,157,193,211,281,337,421,601\}$; witnesses restricted to weight-$\le3$
supports and the normalized reduced-basis candidate class of
$V\cap\mathbb F^S$; per-witness $\delta$ exhaustively closed at the
best-only witness of each instance (all cheaper line-tuples absent, one
present at cost).

NOT swept, plainly: supports of weight $\ge4$; non-basis field values inside
$\dim>1$ intersections (reduced-basis candidates only); other $t$; other
$\eta$ whose constraint binds ($t_i>1$ anywhere); non-identity $\Lambda$
(the conjecture's $\forall\Lambda_i$ quantifier untouched); OTHER $(q,s)$
pairs not listed above — the frozen 13 rows themselves are NOT re-run (they
are complete, exact, MACHINE-VERIFIED and one of this campaign's
non-goals); continuous $\beta$; every asymptotic claim in $q$, $N$, $\beta$;
global $\rho_{\rm inst}$ (each row certifies only
$[1/(3N),\rho^{\rm window}]$); TR26-150 itself — nothing here touches the
conjecture's truth or statement (falsifiability stance, parent §1.9).

## P2-6. Budget, phase order, escalation (fixed now)

- Budget: single workstation; census priced above; per-candidate MILP cap
  120 s; instance-level outer deadline 12 h wall; checkpoints to
  `scratch/p2_checkpoint.json` after EVERY instance (write-first; two
  predecessors lost work to mid-flight failures).
- Phase order (COMMITTED):
  1. OC1/OC2 re-derivation from the frozen table (cheap; file reading +
     `verify_record` on frozen witnesses; zero new optimization).
  2. Smoke: unit anchor $(5,(2,2,2),(1,1,1))$ exhaustive
     $\rho^{\rm window}=1/3$ must reproduce with this campaign's code (it
     is MACHINE-VERIFIED frozen); then PP1 (smallest N) must reproduce its
     P2-2 prediction or the battery HALTS and escalates.
  3. Arm PP in id order; 4. Arm PN in id order; largest first if the
  per-instance outer deadline forces cutbacks — cutbacks recorded as OPEN,
  never failures.
- Escalation: any numeric disagreement with a frozen row (impossible by
  construction — no overlapping instances), any PP/PN row needing a
  statement about TR26-150 itself, or any anchor failure → hub Main with
  the exact discrepancy BEFORE any claim lands in any README. A P- or
  H-COP-falsifying row is NOT an escalation: it is a battery result (the
  repo's purpose), reported prominently.

## P2-7. Pre-commitment on reporting

Every instance is reported, including and especially misses. The number
comes first, the interpretation second (inherited from GB9.2's reporting
rule). No value anywhere in the battery may be re-derived with a different
window or candidate class if it fails to match a prediction (rule 16).

— RsPe3dPatternP, 2026-08-30T18:04Z. This file precedes every battery number.
