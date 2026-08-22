# A rate–distance ceiling for perturbed bivariate-bicycle codes, and its circuit-level cost

*Draft, 2026-08-13 (rev. 3: basis-independence, certified reversals, schedule-class separation).*

## Abstract

Perturbed bivariate-bicycle (PBB) codes augment the $X$-type stabiliser block of
a bivariate-bicycle (BB) code with $Z$-type support from two extra ring elements
$C,D$, producing genuinely non-CSS stabilisers. They were introduced as a way to
reach rate–distance regimes that the CSS structure forbids. We prove a
structural ceiling that rules this out for the majority of the published
catalogue, and measure what the construction costs at the circuit level.

Writing $H_X=[A\,B]$, $H_Z=[B^{\mathsf T}A^{\mathsf T}]$, $P=[C\,D]$, we show
that the pure-$Z$ centraliser of a PBB code is $\ker H_X$ **independently of $C$
and $D$**, that its pure-$Z$ stabiliser subgroup is
$\langle H_Z\rangle+\{uP: uH_X=0\}$, and hence that
$k_{\mathrm{PBB}}=k_{\mathrm{BB}}-\delta$ for an integer $\delta\ge0$. When
$\delta=0$ the nontrivial pure-$Z$ logical operators of the PBB code are
*identical* to those of its parent BB code, so
$d_{\mathrm{PBB}}\le d_Z(\mathrm{BB}(A,B))$; combined with the observation that
every BB code satisfies $d_X=d_Z$ (via an explicit qubit involution), the parent
CSS code weakly dominates the PBB code in $[[n,k,d]]$. This applies to 213 of
the 368 published PBB codes, including every member of the headline
$[[144,12,12]]$ family. For arbitrary $\delta$ we exhibit a *CSS shadow* on the
same qubits with the same $k$ whose $Z$-distance upper-bounds the PBB distance;
that is a ceiling rather than domination. A strict-provenance audit of all 155
$\delta>0$ catalogue codes then shows domination is **false in general**: in
two certified cases a $\delta=4$ perturbation of a $[[72,8,4]]$ parent yields a
genuinely non-CSS $[[72,4,6]]$ — the perturbation strictly *increases* the
distance (both distances exact, independently re-verified) — alongside 66
proved dominations and 87 instances beyond that audit's certification budget.
A CDCL closure campaign (EXP-036) then decided all 29 of the undecided
$n=108/144$ instances and 5 more at $n=180/360$: 29 further dominations and
**five further certified reversals** at $\delta\in\{2,4,6\}$ — running totals 95
dominations, **7 reversals**, 41 undecided (all at $n\in\{180,360\}$, in flight).  Every
reversal code remains dominated by a CSS code at equal $n$, so the CSS
envelope survives (machine-checked: **all seven** certified reversals — EXP-036's five plus EXP-027's two — have an explicit dominating certified-exact CSS code at equal $n$, `results/processed/exp036_envelope_check.json`); the mechanism "perturbations never buy distance" fails at
three lattice sizes.

At the circuit level we derive the correctness criterion for one-ancilla
syndrome extraction of arbitrary mixed stabilisers — the number of
anticommuting shared qubits at which one check precedes another must be even for
every pair — and show it *forces* depth 7 for weight-6 BB codes within the
translation-invariant schedule class, reproducing the depth-7 cycle of Bravyi
et al. in that class (the unrestricted depth-6 exclusion is certified
externally by ASC). For the $[[144,12,12]]$ family
we prove the depth separation is **basis-independent**: no stabiliser element
of symplectic weight $\le 7$ has a nonzero $X$-part (solver-certified with
completeness certificates for all 14 members), so *every* generating set
contains a weight-$\ge8$ generator — in the optimal-code-weight terminology of
Wei et al., $W\ge8$ for all 14 members versus $W=6$ exactly for the Gross code
— and every one-ancilla circuit needs $\ge 8$ layers regardless of which
generators it measures. Enumerating single-ancilla hooks in the $n{=}72$ pair,
58.3 % of them in the PBB `phase2_47` circuit corrupt both Pauli sectors at
once, against 0 % for the CSS $[[72,12,6]]$ circuit. In a
schedule-controlled matched simulation — minimum-depth translation-invariant
schedules exhaustively
enumerated, five sampled uniformly per code, slot maps pinned — all five PBB
point estimates exceed all five CSS ones (Mann–Whitney $p=0.0040$;
schedule-averaged ratio 3.55×, Welch $p=0.014$), although the predeclared
familywise pointwise test does not separate them. Decoding is 6–17× slower on
an idle machine (p50 10.3×, p95 16.6×, p99 6.3×).

We also report an incidental finding that constrains all such comparisons:
schedule choice alone moves the logical error rate by 2.8× *within* a single
code, so any qLDPC circuit-level claim that does not pin its schedule is
under-specified.

## 1. Introduction

*(Motivation: qLDPC overhead, the Gross code, the appeal of escaping the CSS
constraint. State the three contributions: the rate/distance ceiling theorem,
the scheduling criterion and depth separation, the matched end-to-end
simulation. Emphasise that we reproduce rather than merely cite.)*

## 2. Preliminaries

Binary symplectic representation; $R=\mathbb F_2[x,y]/(x^\ell-1,y^m-1)$; BB and
PBB definitions; the distinction between code distance, circuit distance and
decoder estimates.

## 3. Structure of PBB codes

Lemma 0 (validity $\iff$ $AC^{\mathsf T}+BD^{\mathsf T}$ symmetric — agreeing
with, and independently re-deriving, the original definition). Propositions 1–3
($\delta$ and the rate identity). Theorem 1 and Corollary 1 (distance ceiling
and parent domination). Proposition 4 ($d_X=d_Z$ for BB). Theorem 2 (CSS shadow).
Limitation: $\delta>0$ codes need not be dominated by *their own* shadow;
`phase2_58`, `phase2_60` are explicit counterexamples, and we say so.

## 4. Syndrome extraction for mixed stabilisers

Lemma C1 and its proof. The CP-SAT scheduler. Verification protocol: zero
detection events over 4000 noiseless shots, on $[[5,1,3]]$, $[[7,1,3]]$,
$[[9,1,3]]$ and all BB/PBB instances. Proposition C2 and Theorems C3, C4.
Proposition C5 on translation-invariant schedules, with the explicit note that
the restriction, not the code, is what blocks 83 catalogue entries.

## 5. Hook-error mechanism

The two-sector hook count. Table comparing CSS and PBB at $n=72$.

## 6. Matched circuit-level benchmark

Protocol, noise model, decoder, statistics. Resource table. LER table with
Clopper–Pearson intervals. Explicit quarantine of timing measured under shared
machine load, and the isolated re-measurement protocol.

## 7. Discussion

What generalises: any construction that leaves one Pauli sector's centraliser
untouched inherits a distance ceiling from its CSS parent. What does not: this
says nothing about non-CSS codes obtained by genuinely mixing both sectors, nor
about flag or cat-state extraction, which change the ancilla budget.

## 8. Reproducibility

Pinned commits, package versions, seeds, commands, and all twenty-one
failed-route entries recording claims or workflows that were falsified,
downgraded, or repaired during the work.

---

### Claim ledger for this draft

| # | Claim | Status |
|---|---|---|
| 1 | $k_{\rm PBB}=k_{\rm BB}-\delta$ | proved; 368/368 verified |
| 2 | $\delta=0\Rightarrow$ pure-$Z$ logicals identical to parent | proved; exact check on `12_6_0187` |
| 3 | every BB code has $d_X=d_Z$ | proved; 375/375 verified |
| 4 | parent CSS code dominates for $\delta=0$ | proved; 213/368 codes |
| 5 | CSS shadow, same $(n,k)$, $d_{\rm PBB}\le d_Z(Q')$, all $\delta$ | proved; $d_Z$ identity 6/6 exact. A **ceiling, not domination** — it does not bound $d_X(Q')$ |
| 6 | some $\delta>0$ codes escape their own shadow | exact; 2 explicit examples |
| 6b | **$\delta>0$ domination is false in general**: 2 catalogue codes have certified $d_{\rm PBB}=6>4=d_{\rm parent}$ ($[[72,8,4]]\to[[72,4,6]]$, $\delta=4$) | exact both sides (OPTIMAL, witnesses verified, double-checked); still dominated by CSS $[[72,12,6]]$ at equal $n$; audit: 66 dominated / 2 reversals / 87 undecided of 155 |
| 6c | **five further reversals** (EXP-036): `phase2_71`/`72` ($\delta=6$, $n=108$), `phase2_88`, `9_6_0183` ($\delta=2$, $n=108$, $[[108,4,\le8]]\to[[108,2,>8]]$ witness 10), `12_6_0217` ($\delta=4$, $n=144$, $[[144,8,\le8]]\to[[144,4,>8]]$ witness 10) | certified: verified witness upper bounds + replayed CDCL UNSAT lower bounds, CNF-hash-bound stamps, tie-safe protocol; running audit totals 90 / 7 / 58 of 155 |
| 7 | one-ancilla correctness criterion | proved; 0 detector firings on 3 toy codes |
| 8 | criterion forces depth 7 for weight-6 BB **within the TI schedule class** | exact; TI-model $T{=}6$ `INFEASIBLE` on 4 instances; achievement of 7 is class-free; unrestricted 6-exclusion is external (ASC) |
| 9 | PBB $[[144,12,12]]$ needs depth $\ge8$ **for every generating set** (Theorem C6) | proved + exact: no weight-$\le7$ element has nonzero $X$-part (14/14 INFEASIBLE, completeness-certified); independent rebuild checks ranks, witnesses, bounded combinations and coset controls for 14/14; a second monolithic encoding agrees on benchmark member `12_6_0193`; $\mathrm{rank}(V_7)=66<132$; Gross $W=6$ exact |
| 10 | 58.3 % of PBB hooks are two-sector ($n{=}72$ pair: `phase2_47` vs CSS $[[72,12,6]]$; no $n{=}144$ enumeration exists) | exact enumeration |
| 11 | all 5 sampled PBB schedules worse than all 5 sampled CSS schedules; schedule-mean ratio 3.55× | statistical over the exhaustively enumerated **translation-invariant** minimum-depth class; Mann-Whitney $p=0.004$, Welch $p=0.014$; predeclared Bonferroni pointwise test **FAILED**; no best-schedule claim |
| 12 | PBB decoding 6–17× slower (p50 10.3×, p95 16.6×, p99 6.3×) | statistical, isolated machine |
| 12b | cause is *not* differential BP convergence | measured and **refuted**: 0 % BP convergence for both codes |
| 13 | schedule choice alone moves LER 2.8× within each code | statistical; forced a redesign of our own benchmark |
| 14 | two-sector hooks cause the LER gap | **not established** — structural difference only; no fault-injection study run |
| 15 | one-ancilla min depth $=w_{\max}+[\text{exact parity-CSP obstruction}]$, two-value law | predictor exact on all 208 decided EXP-031 cases (150 unrestricted random sets + 47 BB and 11 PBB translation-invariant; proxy confusion matrix spans all 211 attempted); **FALSIFIED in the translation-invariant class** by 3 weight-9 PBB members (no TI schedule through depth 13), but EXP-033 finds unrestricted depth **9 OPTIMAL** for all three. Thus $T_{\rm unrestricted}=9$ while $T_{\rm TI}>13$ for the adversarial trio; unrestricted law supported, not universally proved |
| 16 | undecomposed-DEM mechanism distance of both $n{=}144$ circuits | certified $d_{\rm DEM\,mech}\ge4$ (weights $\le3$ exhaustive); validated weight-12 witnesses give $d_{\rm DEM\,mech}\le12$ both; **not** an exact physical-location circuit-distance claim; heuristics uninformative |
| 17 | **mapped** surface $d{=}12$ baseline at matched $k{=}12$, numeric $p{=}0.002$ | measured: $Z$-memory $2.5\times10^{-5}$ (5/200k, 95% CI $[0.8,5.8]\times10^{-5}$), $X$-memory $2.0\times10^{-5}$ (4/200k, CI $[0.5,5.1]\times10^{-5}$) LER/round on 3444 qubits vs Gross $7.4\times10^{-4}$ (CI $[3.5,11.4]\times10^{-4}$) on 288 (schedule-sample mean over 5 TI schedules, exact block$\to$round conversion; 12× fewer qubits; LER-ratio point estimate ~30×, interval-compatible with ≈6–140× at 4–5 events); **not identical noise locations or decoder**; all convention deltas enumerated |
| 18 | catalogue has no duplicates under proven symmetry group | exact; 368→368 classes, 14-family mutually distinct |
| 19 | benchmark target `12_6_0193` is genuinely non-CSS (Theorem E) | proved + exact (EXP-034): not CSS after row ops ($20{+}66{=}86<132$), any qubit permutation, any local-H assignment, **or any independent local Clifford** — nineteen parity masks plus one one-hot parity row XOR to the GF(2) contradiction $0{=}1$ (persisted, machine-rechecked; CP-SAT independently INFEASIBLE); row space indecomposable under row ops, permutations, and local Cliffords; $[[5,1,3]]$ stays linear-feasible, proving the certificate is not vacuous |
| 20 | independent target distance for `12_6_0193` | **proved + exact — Theorem F** (EXP-035): $d=12$ two-sided — exact MITM excludes weight $\le5$; complete weight-6 classification (72/72 weight-6 zero-syndrome vectors in $S$, equal to the stored pure-$Z$ checks — no weight-6 logical); weights 7–11 excluded by four UNSAT proofs from CaDiCaL on equisatisfiable CNFs (502/674/422/576 s single-threaded) over the machine-checked translation-orbit reduction (transported dual rank 24 on two GF(2) paths); independently re-verified weight-12 witness; the catalogue's `deep_milp` value is confirmed, not trusted; CP-SAT (7 workers, multi-hour budgets) proved none of the four — the frozen dual-backend protocol (INFEASIBLE/UNSAT proof pairs, fail-closed collector) is itself a finding |

### Retractions made during this work

Twenty-one failed-route entries are retained in the append-only ledger.  The
sixteen publication-critical corrections are summarized here:

* **FR-001** a distance refutation that used the wrong nontriviality detector;
* **FR-002** a scheduler that was conflict-free but measured the wrong operators;
* **FR-005** an unschedulability claim that was an artifact of the
  translation-invariant schedule class;
* **FR-006** an erratum we nearly filed against a paper that was correct;
* **FR-007** a benchmark that did not hold the circuit fixed across noise points;
* **FR-008** the BP-convergence explanation for the decoder slowdown, measured
  and refuted (0 % convergence for both codes).
* **FR-012** canonical artifacts could be clobbered by partial runs (three
  successive guard defects; final gate coverage-aware and exercised);
* **FR-013** a verdict field that reported counterexamples as POSITIVE;
* **FR-014** certified translation-invariant refutations misfiled as
  "undecided" — correcting this exposed a genuine schedule-class separation:
  the three codes have $T_{\rm TI}>13$ but EXP-033 proves
  $T_{\rm unrestricted}=9$ exactly.
* **FR-015** an aggregate coverage sentence that could be read as 741 tests
  across three lattices; the actual evidence is 247 perturbations across three
  $(6,6)$ parents and zero coverage on $(6,4)$ or $(4,6)$.
* **FR-016** Theorem C3 and EXP-016 silently dropped their
  translation-invariant schedule-class scope; every statement is now scoped,
  and the audit table records the failure rather than denying it;
* **FR-017** legacy EXP-012 called five parents "weaker" using catalogue
  distances in the wrong bound direction; only two survive EXP-027's strict
  certification and the other three are undecided.

* **FR-018** EXP-026's 1.16-million-decode Cartesian grid failed its
  feasibility gate before any production row completed; bounded, resumable,
  isolation-gated protocol v2 replaced it before results;
* **FR-019** EXP-025 consumed its stage-1 budget on 70 Gross schedules and
  completed zero PBB rows; the asymmetric data are development-only and a
  symmetric held-out v2 was frozen before any PBB result.
* **FR-020** EXP-024's frozen-v2 execution gate was dead code on the normal
  CLI, which could have benchmarked through the retired legacy loop; the
  production entry now fails closed before any structural work and the
  canonical benchmark remains unexecuted;
* **FR-021** EXP-033's noiseless evidence was not content-bound to its
  nondeterministically regenerated slot maps; the artifact now persists and
  tests the canonical-input SHA-256 and per-row slot hashes, and the checks
  were regenerated under that binding.
Each is recorded in `notes/failed_routes.md` with the corrected statement.
