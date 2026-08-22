# Structural ceilings, certified reversals, and circuit costs in perturbed bivariate-bicycle codes

**Date** 2026-08-13 · **Repository** `qec-codesign`

Every claim below carries an epistemic label:
**[proved]** · **[exact computation]** (GF(2) arithmetic, or a solver that
returned proven `OPTIMAL`/`INFEASIBLE`) · **[certified bound]** ·
**[statistical]** (Monte Carlo with an interval) · **[external]** (someone
else's result, reproduced or not) · **[conjecture]**.

---

## 1. What was asked and what was found

The assignment was to determine whether a code, syndrome circuit, decoder and
hardware mapping can be co-designed to beat bivariate-bicycle (BB) and surface
baselines under circuit-level noise, with the non-CSS **perturbed
bivariate-bicycle** (PBB) construction as the concrete target, and specifically
whether the non-CSS PBB $[[144,12,12]]$ has any end-to-end advantage over the
CSS Gross code $[[144,12,12]]$.

**Answer: it does not, and the reason is structural rather than accidental.**

1. **[proved]** A PBB code can never encode more logical qubits than its parent
   BB code: $k_{\rm PBB}=k_{\rm BB}-\delta$ with $\delta\ge0$.
2. **[proved]** When $\delta=0$ the nontrivial pure-$Z$ logical operators of the
   PBB code are *literally the same set* as the parent's, so
   $d_{\rm PBB}\le d_Z({\rm parent})=d({\rm parent})$. The parent CSS code
   therefore weakly dominates in $[[n,k,d]]$. This covers **213 of the 368**
   published PBB codes and **all 14** members of the headline $[[144,12,12]]$
   family.
3. **[proved, weaker]** For arbitrary $\delta$ there is an explicit CSS code on
   the same qubits with the same $k$ whose *$Z$-distance* upper-bounds
   $d_{\rm PBB}$.  This is a ceiling, **not** domination.  The strict-provenance
   audit over all 155 $\delta>0$ codes (EXP-027) settles the shape:
   **66 provably dominated by the parent; 2 certified reversals**
   (`phase2_58`, `phase2_60`: a $\delta=4$ perturbation of a $[[72,8,4]]$
   parent yields non-CSS $[[72,4,6]]$ — the perturbation **strictly increased
   the distance**, both distances exact/OPTIMAL, double-verified); **87
   undecided** at certification budget (immutable artifact).  **EXP-036**
   then turned the Theorem-F CDCL method on those 87: **24 further
   dominations and 5 further replay-stamped reversals** (`phase2_71`/`72`,
   `phase2_88`, `9_6_0183` at $n=108$; `12_6_0217` at $n=144$) — running
   totals **107 dominations / 7 reversals / 41 undecided**, the remainder all
   at $n\in\{180,360\}$ and in flight.  Every reversal is still dominated by
   a CSS code at equal $n$, so the CSS envelope stands (machine-checked: **all seven** certified reversals — EXP-036's five plus EXP-027's two — have an explicit dominating certified-exact CSS code at equal $n$, `results/processed/exp036_envelope_check.json`); what falls is
   "perturbations never buy distance" — now at three lattice sizes.
   $\delta=0$ domination (item 2) is untouched.
4. **[proved + exact computation, basis-independent]** For every catalogue PBB
   $[[144,12,12]]$ member, **every generating set** of the stabiliser group
   contains a generator of symplectic weight $\ge 8$ (EXP-023 / Theorem C6:
   no weight-$\le7$ element has a nonzero $X$-part — CP-SAT INFEASIBLE with
   completeness certificates; an independently rebuilt artifact rechecks all
   ranks, stored witnesses, bounded row combinations, and coset controls;
   pure-$Z$ subgroup rank exactly $66<132$).  A second monolithic encoding was
   run only for benchmark member `12_6_0193`, where it agreed.  Hence **any
   one-ancilla syndrome circuit needs $\ge 8$ two-qubit layers** regardless of
   which generators it measures; the Gross code's optimal max generator weight
   is exactly 6 ($\mathrm{rank}(V_5)=0$) and it achieves 7 layers — minimal
   for its published generators within the translation-invariant schedule
   class (Theorem C3 scope; the unrestricted depth-6 exclusion is external,
   ASC arXiv:2603.21499).  The $\ge8$-vs-7 separation is a property of the
   codes, not of their presentations.
5. **[exact computation, $n=72$ pair]** 58.3 % of single-ancilla hooks in the
   PBB `phase2_47` $[[72,12,6]]$ circuit corrupt both Pauli sectors at once;
   0 % do in the CSS $[[72,12,6]]$ circuit.  (No hook enumeration exists for
   the $n=144$ pair.)
6. **[statistical]** With the schedule controlled — minimum-depth
   **translation-invariant** schedules *exhaustively enumerated*, five drawn
   *uniformly at random* from each
   complete set, slot maps pinned and persisted — **all five** sampled PBB
   point estimates exceeded **all five** sampled CSS point estimates
   (Mann-Whitney one-sided $U=25/25$, $p=0.0040$), and the schedule-averaged
   rate was 3.55× higher (Welch $p=0.014$).  The **predeclared stricter test
   failed**: with Bonferroni-corrected simultaneous intervals at familywise
   95 % the smallest PBB lower bound ($1.30\times10^{-2}$) overlaps the largest
   CSS upper bound ($1.91\times10^{-2}$), so pointwise separation of every pair
   is **not** established.  We did not optimise, and did not identify, the best
   schedule for either code.
7. **[statistical, isolated]** At matched decoder settings on an idle machine
   the PBB circuit decodes 10.3× (p50) / 16.6× (p95) / 6.3× (p99) slower.  Our
   proposed cause — worse BP convergence — was **measured and refuted**: BP
   converges on 0.0 % of shots for *both* codes.

**The defensible summary.**  The rate and distance results (1-4) are theorems
with exact verification and do not depend on any simulation.  The circuit-level
results (5-7) are a *consistent body of evidence* pointing the same way, not a
single decisive measurement: the PBB code needs more depth and more gates by a
counting argument, produces a fault mode CSS cannot, was worse in every one of
five sampled schedules, and decoded an order of magnitude slower.  What we do
**not** have is pointwise familywise separation, an optimised schedule for
either code, or an isolated causal mechanism for the decoder slowdown.

Scoped precisely: for the **213 of 368** catalogue codes with $\delta=0$ —
including every member of the headline $[[144,12,12]]$ family — the non-CSS
degree of freedom provably buys nothing in $[[n,k,d]]$.  For all $\delta$ it
provably cannot raise $k$, and it is capped by a same-$(n,k)$ CSS code's
$Z$-distance.  For $\delta>0$, domination is **false in general**: EXP-027
certifies two catalogue reversals ($[[72,8,4]]\to[[72,4,6]]$) and the EXP-036
SAT closure certifies five more at $n=108/144$, alongside 90 proved
dominations and 58 still undecided ($n\in\{180,360\}$, in flight).  No
circuit-level result favoured the $[[144,12,12]]$ PBB
candidate: depth, gate count, hook structure, sampled LER, and latency moved
against it, while the certified DEM-mechanism distance intervals tied.  The
one-ancilla depth separation is now basis-independent (Theorem C6).

---

## 1b. The 2026-08-12 extension campaign (eleven parallel directions)

Eleven independent investigations were run concurrently by subagents against
pre-registered protocols, each with its own artifact and report under
`notes/agent_reports/`.  Summary, with verification status:

| # | direction | verdict | one-line result |
|---|---|---|---|
| EXP-023 | min-over-bases generator weight | **BREAKTHROUGH_CANDIDATE** | Theorem C6: all 14 PBB members need a weight-$\ge8$ generator in *every* basis; Gross $w^*{=}6$ exact |
| EXP-024 | cat/flag extraction (open Q4) | **NEGATIVE_STRUCTURAL**; benchmark **BLOCKED** | 2-ancilla 4+4 split: total-2q depth 7 **INFEASIBLE** (certified); TI depth-9 built, semantics verified; canonical Monte Carlo benchmark never executed — normal CLI fails closed until the frozen v2 driver exists (FR-020) |
| EXP-025 | schedule optimisation | **STOPPED before PBB data** | 70/70 Gross stage-1 rows, 0/69 PBB; asymmetric v1 is development-only; symmetric held-out v2 frozen unrun (FR-019) |
| EXP-026 | decoder co-design + shadow decoder | **STOPPED at feasibility gate** | 1.16M-decode grid failed its feasibility gate before any production row (FR-018); frozen v2 (finite BP+OSD/BP+LSD set, resumable, deadline-bounded) not yet run |
| EXP-027 | strict $\delta>0$ domination audit | **BREAKTHROUGH_CANDIDATE** | 66 dominated / **2 certified distance-increasing reversals** / 87 undecided |
| EXP-028 | out-of-catalogue reversal hunt | NEGATIVE, tightly scoped | no reversal among **247 perturbations across 3 parents**; every parent/PBB aggregate distance call returned `OPTIMAL` (some sector subproblems are `INFEASIBLE`, as expected); only $(6,6)$ reached before the reserve-triggered stop under the 90-min wall budget; no evidence for $(6,4)$ or $(4,6)$; coverage: 3 of 300 sampled parents reached, 247 of 844 enumerated out-of-catalogue perturbations tested — non-exhaustive |
| EXP-029 | undecomposed-DEM mechanism distance at $n{=}144$ | POSITIVE | certified $d_{\rm DEM\,mech}\ge4$ both circuits (weights ≤3 exhaustive/MITM); validated weight-12 witnesses give $d_{\rm DEM\,mech}\le12$ both; no exact physical-location circuit-distance claim; heuristic searches timed out |
| EXP-030 | mapped real surface-code baseline | POSITIVE | rotated $d{=}12$, matched-$k$: $Z$ $2.5\times10^{-5}$ (5/200k, 95% CI $[0.8,5.8]\times10^{-5}$), $X$ $2.0\times10^{-5}$ (4/200k, CI $[0.5,5.1]\times10^{-5}$) LER/round at $p{=}0.002$ on 3444 qubits vs Gross $7.4\times10^{-4}$ (CI $[3.5,11.4]\times10^{-4}$, 5-schedule t-interval) on 288; ratio point estimate ~30× but interval-compatible with ≈6–140× at 4–5 events; **not identical noise/decoder conventions** |
| EXP-031 | one-ancilla depth criterion | **NEGATIVE** (criterion falsified in the TI class) | exact predictor on all 208 decided cases (150 unrestricted random sets, 47 BB + 11 PBB translation-invariant); 3 weight-9 PBB members certified `INFEASIBLE` at every TI depth 9–13 |
| EXP-033 | unrestricted falsifier follow-up | POSITIVE | all three adversaries have unrestricted depth **9 OPTIMAL**: $T_{\rm TI}>13$ vs $T_{\rm unrestricted}=9$ — the failure is a cost of translation invariance, not of the codes |
| EXP-032 | catalogue dedup under proven symmetries | POSITIVE | 368 → 368 classes; the 14-family is mutually distinct |
| EXP-034 | target equivalence (`12_6_0193`) | **POSITIVE, fully decided** | genuinely non-CSS (Theorem E): not CSS after row ops (86<132), any qubit permutation, any local-H ($0{=}1$ certificate), or **any independent local Clifford** — 19 parity masks + 1 one-hot parity row XOR to $0{=}1$ over GF(2), a solver-free linear certificate; CP-SAT independently INFEASIBLE; row space indecomposable under the same group; $[[5,1,3]]$ linear-feasible control shows non-vacuousness |
| EXP-035 | independent target distance (`12_6_0193`) | **POSITIVE, exact — Theorem F** | $d=12$ certified two-sided: exact MITM excludes weight $\le5$; complete weight-6 classification (72/72 zero-syndrome vectors in $S$, equal to the stored pure-$Z$ checks — no weight-6 logical); weights 7–11 excluded by four UNSAT sector proofs (CaDiCaL/PySAT equisatisfiable CNF, 502/674/422/576 s) over the machine-checked translation-orbit reduction (union dual rank 24, two GF(2) paths); verified weight-12 witness; CP-SAT racing the same sectors produced no proof — CDCL, not CP, closed the problem |
| EXP-036 | SAT closure of the $\delta>0$ gap | **POSITIVE, in flight** | Theorem-F CDCL method on EXP-027's 87 undecided rows: all 29 at $n=108/144$ decided plus 17 at $n=180/360$ — 41 dominations + **5 new replay-stamped reversals** (`phase2_71`/`72`, `phase2_88`, `9_6_0183`, `12_6_0217`); running totals 107/7/41; a sound one-clause translation symmetry break gives 6.2x on decisive UNSAT proofs (premise machine-checked) while sector decomposition measured 8x slower and was abandoned (FR-022); every verdict = re-verified witness + replayed decisive UNSAT with CNF-hash-bound stamps; $n\in\{180,360\}$ waves running |
| EXP-037 | catalogue-wide CSS-envelope classification | **POSITIVE, in flight** | asks the programme's central question directly: does *any* of the 368 published PBB codes escape the CSS BB envelope at its own length? A row $Q=[[n,k,\cdot]]$ is envelope-dominated when some same-length CSS BB code has $k_C\ge k_Q$ and certified exact $d_C\ge U(Q)$, where $U(Q)$ is a *verified witness* upper bound. The asymmetry makes the sweep affordable: the PBB side needs only a cheap SAT witness, and the expensive UNSAT work is spent once per CSS code and amortised over every row at that length. **250 of 368 rows certified envelope-dominated so far, zero genuine escapees** — the 22 flagged rows all have no certified same-length CSS code with $k_C\ge k_Q$ yet (pool coverage, not a finding), 96 rows pending. Both domination dimensions are re-derived from rebuilt algebra at classification time ($k$ from the fresh stabilizer, $U$ by re-verifying the stored witness); the CSS dimension uses $k=n-\mathrm{rank}(H_X)-\mathrm{rank}(H_Z)$ cross-checked against the independent theory module. Bounded witness descent additionally yields exact PBB distances for free: reversals `phase2_71`, `phase2_72`, `phase2_88` now carry persisted UNSAT-at-5 records alongside their weight-6 witnesses, so each is **exactly $d=6$** against a parent $\le4$.  **First escape leads:** four $n=180$, $k=6$ rows (`15_6_0218/0219/0222/0229`) have verified witnesses at weight 14-16 while the best certified same-length CSS code at $k_C\ge6$ reaches only $d=12$. These are *leads, not claims*: an escape needs a PBB **lower** bound (UNSAT at 12, running as dedicated probes) plus completion of the $n=180$ CSS frontier (2 of 46 distinct parents certified so far). |
| EXP-038 | structural survival criterion (**Theorem G**) | **POSITIVE, proved + machine-verified** | explains *why* reversals are rare, and turns the explanation into a zero-search domination test.  The PBB's pure-$Z$ centralizer equals the parent's ($\ker[A\,B]$); its pure-$Z$ stabilizers are the parent's plus a dressing space $\Delta=\{\lambda[C\,D]:\lambda[A\,B]=0\}$; and $k_Q=k_P-\dim\Delta$ is **proved** (via $\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\Delta$) and verified with zero mismatches on all 368 rows.  Hence any surviving minimum-weight parent $Z$-logical certifies $d_Q\le d_Z(P)$ with **no search on the PBB** - the survivor is the certificate.  Counting form: $d_Q>d_Z(P)\Rightarrow k_Q\le k_P-t$ for $t$ the rank of the minimum-weight translation orbit.  Falsification gate passed: the criterion is **silent on all certified reversals**, and on each the margin $t-\dim\Delta$ is **exactly $0$** - reversals sit on the knife edge.  Explains the observed confinement of all 7 reversals to the $k$-halving class (58 rows), none among the 97 rows at $k_Q/k_P\in\{3/4,5/6\}$.  Runtime collapses from 785 s to 0.3-4 s per row once the parent's CSS distance is pooled |
| prior art | novelty refresh (2026-08-12) | partial scoop | ASC (2603.21499) certifies no-depth-6 for IBM BB codes; our analytic criterion, C6, reversals, and end-to-end negative remain unscooped |

Process defects found and fixed during the campaign, each with a failed-route
entry and regression test: FR-012 (canonical-artifact clobber by partial runs;
write routing now `canonical_route`-gated and exercised), FR-013 (verdict
polarity), FR-014 (certified TI refutations misfiled as "undecided"; correction
exposed the real schedule-class separation, resolved by EXP-033).
FR-015 (an aggregate coverage sentence that multiplied counts and implied
untested lattices were covered) was caught and fixed after the campaign wrote
its first reports; the 2026-08-13 adversarial audit then found one more
defect of the FR-005/FR-014 class — the TI schedule-class scope silently
dropped from Theorem C3 and EXP-016 prose — which is corrected throughout
this revision.  Post-campaign hardening added FR-018/019 (the non-executable
decoder and schedule grids), FR-020 (EXP-024's frozen-v2 gate was dead code on
the normal CLI; production now fails closed before any structural work), and
FR-021 (EXP-033's noiseless evidence was not content-bound to its slot maps;
the artifact now stores and tests the canonical-input SHA-256 and per-row
slot-map hashes, and the checks were regenerated under that binding).


---

## 2. Setting

$R=\mathbb F_2[x,y]/(x^\ell-1,y^m-1)$, $\dim=\ell m$, $n=2\dim$. With
$H_X=[A\;B]$, $H_Z=[B^{\mathsf T}A^{\mathsf T}]$, $P=[C\;D]$,

$$H=\begin{pmatrix}H_X & P\\ 0 & H_Z\end{pmatrix},\qquad v=(x|z),\qquad
\langle(x|z),(x'|z')\rangle_s=x\!\cdot\!z'+z\!\cdot\!x'.$$

$C=D=0$ recovers the CSS BB code, the **parent**.

Full proofs: `proofs/pbb_structure.md`. Implementation:
`src/qec_research/codes/{bicycle,pbb_theory}.py`.

---

## 3. Baseline reproduction (Gates 2, 3, 7)

**[exact computation]** All seven Bravyi BB instances were rebuilt from
polynomial shift operators with no code shared with `qldpc` or the reference
repositories. Commutation, rank, $k$, check weights, qubit degrees and
connectivity were each computed by **two independent GF(2) paths** (Python
`int` bitsets and packed `numpy`), agreeing everywhere.

| code | $n$ | $k$ | max check wt | max degree | components | exact $d$ | solver status |
|---|---|---|---|---|---|---|---|
| $[[72,12,6]]$ | 72 | 12 | 6 | 6 | 1 | $d_X=d_Z=6$ | `OPTIMAL` |
| $[[90,8,10]]$ | 90 | 8 | 6 | 6 | 1 | $d_X=d_Z=10$ | `OPTIMAL` |
| $[[108,8,10]]$ | 108 | 8 | 6 | 6 | 1 | $d_X=d_Z=10$ | `OPTIMAL` |
| $[[144,12,12]]$ | 144 | 12 | 6 | 6 | 1 | $d_X=d_Z=12$ | `OPTIMAL`, 563 s |
| $[[288,12,18]]$ | 288 | 12 | 6 | 6 | 1 | not attempted | — |

Distances come from a CP-SAT integer program per logical sector, using the
identity $d=\min_j\min\{\,\mathrm{wt}_s(v): v\in S^{\perp_s},\ \langle v,L_j\rangle_s=1\,\}$.
Only sectors that return `OPTIMAL` or `INFEASIBLE` count as exact; anything
else is reported as a bound.

**Convention control.** `qcode-discovery` builds monomial matrices as
`lift().T`, the transpose of ours — the ring involution $x^ay^b\mapsto
x^{-a}y^{-b}$. **[proved + exact computation]** the two conventions are related
by the qubit permutation $(i,j)\mapsto(-i,-j)$ and span identical row spaces
(`tests/test_conventions.py`), so all $n$, $k$, weight and distance statements
transfer.

---

## 4. The commutation condition, independently re-derived

**[proved]** The PBB validity condition is that
$M=AC^{\mathsf T}+BD^{\mathsf T}$ be **symmetric**, not zero:
$H_xH_z^{\mathsf T}=\begin{pmatrix}M&AB+BA\\0&0\end{pmatrix}$, and $AB=BA$
because $R$ is commutative, so validity is $M+M^{\mathsf T}=0$.

**[external, verified]** This agrees with the primary source. arXiv:2606.02418v1
§III.2 states: *"All rows commute if and only if $(AC^{\top}+BD^{\top})\bmod 2$
is symmetric over $\mathbb F_2$ (this is the only nontrivial commutativity
condition; commutativity between block 1 and block 2 is automatic from BB ring
commutativity)."* The reference implementation
(`qcode-discovery/evaluation/pbb_code.py:71-73`) tests the same predicate.

An intermediate summary in our own pipeline had paraphrased the condition as
"$AC^{\mathsf T}+BD^{\mathsf T}=0$", which is strictly stronger and would
exclude valid codes; we briefly recorded that as an erratum candidate against
the paper. **Reading the primary source retracted our claim, not theirs.**
Logged as `notes/failed_routes.md` FR-006. **[exact computation]** all 368
catalogue codes satisfy symmetry; our independent rebuild had zero commutation
failures.

## 5. Structural results

### 5.1 Rate

**[proved]** (Props 1–3). The pure-$Z$ centraliser of a PBB code is $\ker H_X$,
independent of $C,D$; the pure-$Z$ stabiliser subgroup is
$S_Z=\langle H_Z\rangle+\{uP: uH_X=0\}$; and with
$\delta=\dim S_Z-\operatorname{rank}H_Z$,

$$k_{\rm PBB}=k_{\rm BB}-\delta,\qquad 0\le\delta\le \dim-\operatorname{rank}H_X.$$

**[exact computation]** Catalogue $k$, this formula, and a direct symplectic
rank agree on **368/368** codes.
$\delta\in\{0{:}213,\;2{:}115,\;4{:}34,\;6{:}2,\;8{:}2,\;10{:}2\}$.

### 5.2 Distance

**[proved]** *Theorem 1.* If $\delta=0$ the nontrivial pure-$Z$ logicals of the
PBB code are exactly the parent's, hence $d_{\rm PBB}\le d_Z({\rm parent})$.

**[proved]** *Proposition 4.* The permutation $(\mathrm L,g)\leftrightarrow
(\mathrm R,-g)$ exchanges $H_X$ and $H_Z$, so **every** BB code has $d_X=d_Z$.
**[exact computation]** verified on 375/375 codes, and independently by exact
CP-SAT on four instances.

**[proved]** *Corollary 1.* $\delta=0\Rightarrow$ the parent CSS BB code has the
same $n$, the same $k$, and $d\ge d_{\rm PBB}$, with max check weight 6.

**[exact computation]** Sharpest single confirmation: for PBB `12_6_0187`
($[[144,12,12]]$, catalogue `trust_level: EXACT`), the minimum pure-$Z$ logical
weight is **12** (`OPTIMAL`, 1086 s) and the parent's is also **12**
(`OPTIMAL`, 272 s) — Theorem 1 predicts equality, and equality is what was
found. Controls in the same run reproduced $6$ and $10$ for
$[[72,12,6]]$ and $[[108,8,10]]$.

**[proved]** *Theorem 2.* For any $\delta$, the CSS shadow $Q'$ ($X$-checks
$H_X$, $Z$-checks $S_Z$) has $n(Q')=n(Q)$, $k(Q')=k(Q)$ and $d(Q)\le d_Z(Q')$.

**Honest limitation.** **[exact computation]** Theorem 2's $d_Z$ identity held
6/6 on the small $\delta>0$ codes tested, but 2 of those 6 — `phase2_58` and
`phase2_60`, both $[[72,4,6]]$ with $\delta=4$ — have $d_X(Q')=4<6$, so their
own shadow does **not** dominate them. They remain dominated by the certified
Bravyi $[[72,12,6]]$ (same $n$, same $d$, $k=12$ against 4), but that is a
different CSS code, and the general-$\delta$ statement is therefore weaker than
the $\delta=0$ one. **[legacy diagnostic, superseded]** An earlier
uncertified sweep (EXP-012, $n\le108$, 77 codes) tallied 49 parent-dominates /
5 parent-weaker / 23 undecided using catalogue labels without certified
bounds.  The strict-provenance audit (EXP-027) supersedes it: over all 155
$\delta>0$ rows, 66 dominations proved, 2 certified reversals, 87 undecided.
EXP-036 then closes 46 of the 87 (41 dominations, 5 reversals; running totals
107/7/41).  Two of EXP-012's five "parent weaker" cases — `phase2_58` and
`phase2_60` — survived certification as the first two reversals; whether the
remaining three coincide with EXP-036's new reversals has not been checked
against the legacy sweep and is not claimed.

---

## 6. Syndrome-extraction circuits

### 6.1 The correctness criterion

**[proved]** *Lemma C1.* In a one-ancilla-per-check circuit,
$\mathrm{C}\text{-}P_{a,j}\,\mathrm{C}\text{-}Q_{b,j}
=\mathrm{C}\text{-}Q_{b,j}\,\mathrm{C}\text{-}P_{a,j}\cdot\mathrm{CZ}_{a,b}$
when $P,Q$ anticommute. Hence the circuit measures the intended generators
**iff** for every pair of checks the number of anticommuting shared qubits at
which the first acts before the second is even.

**[exact computation]** A König edge-colouring of the Steane code violates this
on 5 pairs and fires **1402** detectors in 200 noiseless shots — our own first
scheduler did exactly that (`notes/failed_routes.md` FR-002). A CP-SAT
scheduler that enforces it fires **0** detectors in 4000 noiseless shots on
Steane $[[7,1,3]]$, Shor $[[9,1,3]]$ and the non-CSS perfect code $[[5,1,3]]$,
and on every BB and PBB circuit built here.

### 6.2 Depth

**[proved]** *Proposition C2.* $T\ge\max(\max_a|\mathrm{supp}\,g_a|,\ \max_j\deg j)$.

**[exact computation, TI schedule class]** *Theorem C3.* For BB $[[72,12,6]]$,
$[[108,8,10]]$, $[[144,12,12]]$, $[[288,12,18]]$: $T=6$ is proven `INFEASIBLE`
under Lemma C1 **within the translation-invariant (orbit) schedule class**
(`exp004` solves the orbit model; every `results/raw/exp004_schedules.json`
entry records `"mode": "orbit"`), and $T=7$ is achieved (achievement is
class-free).  Within that class this **derives** the depth-7 cycle of
arXiv:2308.07915 from the correctness criterion rather than assuming it.  Our
artifacts do not exclude an unrestricted depth-6 schedule; ASC
(arXiv:2603.21499) independently certifies that exclusion for the IBM BB
instances.  The audit that caught this silent restriction is §9.

**Independent cross-check (Gate 7).** A second, deliberately independent path
was run: rather than *deriving* a schedule, it *recovered* the published one
straight from the original artifact
(`third_party/BivariateBicycleCodes/decoder_setup.py:45-51`; layer bodies at
`:208-225`, `:227-240`, `:242-256`; measurement/reset slice at `:258-264`):

```
sX = ['idle', 1, 4, 3, 5, 0, 2]
sZ = [3, 5, 0, 1, 2, 4, 'idle']
```

— 8 effective time slices but exactly **7 CNOT layers**. The circuit built from
it uses **288 qubits = 144 data + 72 X-ancillas + 72 Z-ancillas** and exactly
**10 368 CX over 12 rounds = 864 per round**, with **0 detector events and 0
logical flips over 10 000 noiseless shots** (seed 230807915).

Every one of those numbers matches what our parity-criterion scheduler produced
from scratch: depth 7, 288 qubits, 864 two-qubit gates per round, 0 noiseless
detector firings. The two paths share no scheduling code — one solves a
constraint program, the other transcribes a published table — so the agreement
validates the circuit layer independently rather than being two front-ends onto
one library.

**[proved + exact computation]** *Theorem C4.* All 14 catalogue PBB
$[[144,12,12]]$ codes have max check weight $\ge 8$, so every one-ancilla
schedule needs $T\ge8$. The best member reaches exactly 8, which is therefore
optimal for it.

**Retracted claim.** An earlier draft asserted that 83 catalogue PBB codes were
unschedulable. **[exact computation]** the three smallest of them *are*
schedulable once the translation-invariance restriction is dropped
(`phase2_26` at depth 7, `phase2_21`/`phase2_44` at depth 9, each with 0
noiseless detector firings). The infeasibility proofs were sound but belonged
to a schedule class, not to the codes. See `notes/failed_routes.md` FR-005.
The $[[144,12,12]]$ depth separation was re-derived from Proposition C2 so it
does not depend on any schedule class.

---

## 7. Matched circuit-level comparison

Matched on: noise model and strength; rounds; decoder and all its
hyperparameters; observable set; shot budget and stopping rule; worker count.
Differing only in the code and its circuit.

**Resource accounting** (per syndrome round, $n=144$, both with 144 data +
144 ancilla qubits, both with 0 noiseless detector firings in 4000 shots):

| | CSS Gross | non-CSS PBB `12_6_0193` |
|---|---|---|
| max check weight | 6 | 8 |
| mixed checks | 0 | 72 |
| certified depth | **7** | **8** |
| two-qubit gates / round | **864** | **1008** |
| gate types | CX, CZ | CX, CY, CZ |
| $Z$-basis observables | 12 | 12 |

> **Methodological correction (load-bearing).**  Our first sweep (EXP-007)
> re-solved the scheduling CP-SAT model on every call, unseeded and
> multi-threaded, so the circuit was not held fixed across noise points.
> Schedule choice moves the logical error rate by about a factor of two for the
> *same* code: Gross at $p=0.002$, 12 rounds, identical decoder, gave
> $81/8000$ under one schedule and $203/40000$ under another, with disjoint
> intervals.  EXP-007 is therefore **downgraded to arbitrary-schedule
> exploratory evidence**, and the claim that minimal depth singles out a "fair"
> PBB circuit is **withdrawn**.  See `notes/failed_routes.md` FR-007.
>
> EXP-016 replaces it.  The minimum-depth translation-invariant schedules form
> a finite set that CP-SAT enumerates **exhaustively** (status `OPTIMAL`, no
> truncation): **8496** for the Gross code at depth 7 and **9968** for the PBB
> code at depth 8.  We draw a *uniform random sample* from the complete set,
> persist every sampled slot map in full, reuse one map for both the $p=0$
> validity check and the noisy run, and report **Bonferroni-corrected
> simultaneous** Clopper-Pearson intervals at familywise 95 % together with a
> t-interval on the schedule-averaged rate.  Sampled extrema are reported as
> extrema *of the sample*, never as global best or worst.

### 7.1 Schedule-controlled result (primary)

**[statistical]** EXP-016.  The minimum-depth translation-invariant schedules
were **exhaustively enumerated** (CP-SAT status `OPTIMAL`, no truncation):
**8496** for the Gross code at depth 7 and **9968** for the PBB code at depth 8.
Five were drawn *uniformly at random* from each complete set (RNG seed
20260811); every sampled slot map is persisted in full in the result file.  Each
circuit passed the $p=0$ check (0 detector firings, 0 observable flips in 4000
shots) using the *same* slot map as the noisy run.  $p=0.002$, 12 rounds, 8000
shots per schedule, identical decoder.

| code | depth | valid schedules | schedule-averaged LER | sampled range | within-code spread |
|---|---|---|---|---|---|
| CSS Gross | 7 | 8496 | $8.90\times10^{-3}$ | $[5.25,\,14.88]\times10^{-3}$ | 2.83× |
| non-CSS PBB | 8 | 9968 | $3.16\times10^{-2}$ | $[16.75,\,46.88]\times10^{-3}$ | 2.80× |

* **Schedule-averaged ratio 3.55×** (PBB worse).
* Welch two-sample $t$ on the schedule-level rates: $t=3.83$, $p=0.014$.
* Mann-Whitney one-sided: $U=25$ out of a maximum 25, $p=0.0040$ — **every**
  sampled PBB schedule is worse than **every** sampled CSS schedule.
* **The stricter test fails and we report it as failing.**  With
  Bonferroni-corrected simultaneous intervals (familywise 95 %, per-interval
  $\alpha=0.005$), the smallest PBB lower bound is $1.30\times10^{-2}$ while
  the largest CSS upper bound is $1.91\times10^{-2}$; these overlap.  At 8000
  shots per schedule we therefore do **not** establish pointwise separation of
  every pair at familywise 95 %, only the rank separation and the mean
  difference above.

**A finding in its own right.**  Schedule choice alone moves the logical error
rate by **2.8×** *within each code* — comparable to the 3.55× between the code
families.  Two valid, minimum-depth, formally verified schedules for the same
code are not interchangeable.  Any qLDPC co-design claim that does not pin and
publish its schedule is under-specified, and schedule optimisation is an
optimisation axis at least as valuable as the CSS/non-CSS choice.  We did *not*
optimise the schedule for either code; both were sampled uniformly, which is
fair but pessimistic for both.

### 7.2 Exploratory sweep (schedule not controlled)

**[statistical, exploratory]** EXP-007, 12 rounds, per-point Clopper-Pearson
95 % intervals.  Retained for completeness; **not** load-bearing.

| $p$ | code | shots | failures | LER / shot | 95 % CI | LER / round |
|---|---|---|---|---|---|---|
| 0.0015 | CSS Gross | 40000 | 23 | $5.75\times10^{-4}$ | $[3.65,\,8.63]\times10^{-4}$ | $4.79\times10^{-5}$ |
| 0.0015 | non-CSS PBB | 40000 | 143 | $3.58\times10^{-3}$ | $[3.01,\,4.21]\times10^{-3}$ | $2.98\times10^{-4}$ |
| 0.002 | CSS Gross | 40000 | 203 | $5.08\times10^{-3}$ | $[4.40,\,5.82]\times10^{-3}$ | $4.24\times10^{-4}$ |
| 0.002 | non-CSS PBB | 10400 | 249 | $2.39\times10^{-2}$ | $[2.11,\,2.71]\times10^{-2}$ | $2.02\times10^{-3}$ |
| 0.003 | CSS Gross | 10400 | 1171 | $1.13\times10^{-1}$ | $[1.07,\,1.19]\times10^{-1}$ | $9.91\times10^{-3}$ |
| 0.003 | non-CSS PBB | 10400 | 2472 | $2.38\times10^{-1}$ | $[2.30,\,2.46]\times10^{-1}$ | $2.24\times10^{-2}$ |

The direction agrees with EXP-016 at every point, but because the schedule was
not pinned, the magnitudes carry no weight.

**Schedule-independent fact.**  The detector error model is denser for the PBB
circuit under **every** schedule measured: **82 512-82 944** error mechanisms
against **67 032-67 752** for the Gross code at the same $p$ and rounds
(+22 to +24 %).  The DEM size does vary a little with the schedule within each
code (67 032-67 752 for Gross; 82 512-82 944 for PBB), but the two ranges are
disjoint, so the +23 % gap holds under every schedule measured.

### 7.3 Decoder cost, measured in isolation

**[statistical, isolated]** EXP-013.  The EXP-007 latency percentiles were
collected under shared machine load with the codes measured sequentially and
are quarantined (FR-004).  This replacement runs on an effectively idle machine
(1-minute load average 1.83 at start, 2.62 at end, on 28 logical CPUs; the
script refuses to start above a threshold), with `OMP_NUM_THREADS=1`, a single
worker, a discarded warm-up batch, and the two codes **interleaved in
alternating blocks** so any residual drift hits both equally.  Identical decoder
settings, $p=0.002$, 12 rounds.

| | CSS Gross | non-CSS PBB | ratio |
|---|---|---|---|
| DEM error mechanisms | 67 032 | 82 800 | 1.24× |
| decode latency p50 | 118.2 ms | 1220.2 ms | **10.3×** |
| decode latency p95 | 579.2 ms | 9614.2 ms | **16.6×** |
| decode latency p99 | 2159.0 ms | 13629.1 ms | **6.3×** |
| throughput | 4.89 shots/s/core | 0.43 shots/s/core | **11.4× slower** |

The latency penalty is an order of magnitude — far more than the 24 % increase
in detector-error-model size would suggest on its own.

**We tested the obvious explanation and it is wrong** (`notes/failed_routes.md`
FR-008).  EXP-017 instrumented `BpOsdDecoder.converge` directly: at 30 min-sum iterations belief propagation
converges on **0.0 % of shots for *both* codes**, so every shot falls through to
the ordered-statistics stage in both cases.  Differential BP convergence is
therefore **not** the cause.  What we can state is the measured scaling: the PBB
circuit has 1.24× the DEM error mechanisms and **1.37×** the detection events
per shot (163.4 versus 119.2), and its OSD stage costs 10.8× more.  The growth
is superlinear in both quantities but we have **not** isolated which term
dominates, and we do not claim to have.

**Decoder caveat.**  Because BP never converges at these settings, this
comparison is effectively OSD-0 applied to the two detector error models.  It
is matched — identical settings for both codes — but it is *not* the decoder
configuration of the original artifact (OSD-CS order 7 with far more BP
iterations).  A stronger BP stage could change the absolute numbers, and
possibly the ratio.  For a real-time decoder the p99 is what sizes the buffer,
so the gap matters, but it should be re-measured with a production decoder
before being used for hardware planning.

---

---

## 7b. Mechanism: mixed checks corrupt both Pauli sectors at once

What fault mechanism differs between circuits with identical $n$, $k$, and
catalogue $d$?  A single ancilla fault can behave qualitatively differently.

For a pure-$X$ check with the ancilla as control, a mid-sequence ancilla fault
propagates $X$ onto the remaining data targets only: the resulting data error
lives in one sector and is visible to the $Z$-checks.  For a **mixed** check the
same fault propagates $X$ onto the remaining $X$-targets *and* $Z$ onto the
remaining $Z$-targets, producing a weight-$w$ data error with support in both
sectors from a **single** circuit fault.

**[exact computation]** Enumerating every (check, cut point) pair of the
realised schedules — purely combinatorial, no simulation and no solver:

| | CSS-BB $[[72,12,6]]$ | non-CSS PBB $[[72,12,6]]$ `phase2_47` |
|---|---|---|
| syndrome depth | 7 | 8 |
| two-qubit gates / round | 432 | 504 |
| mixed checks | 0 | 36 |
| single-ancilla hooks | 360 | 432 |
| hooks corrupting **both** sectors | **0 (0.0 %)** | **252 (58.3 %)** |
| max hook weight | 5 | 7 |
| mean hook weight | 3.00 | 3.58 |

The PBB circuit has 20 % more hook locations, heavier hooks, and a majority of
them produce two-sector errors — a failure mode *absent by construction* from a
CSS circuit, and a property of mixed stabilizers as such rather than of this
particular code.

**Scope.**  This is a structural difference in the fault set, established by
exact enumeration.  It is a *plausible contributor* to the measured logical-error
gap, not a demonstrated cause: our decoder consumes the undecomposed hypergraph
DEM, which already represents two-sector faults as single error mechanisms, so
we have not shown that two-sector hooks are what the decoder actually fails on.
Establishing causation would need a fault-injection study that conditions
failures on hook type; we did not run one.

**[certified mechanism-model bounds, EXP-029]** On the full $n=144$, 12-round
circuits, exhaustive weight-1/2 search plus a completed exact weight-3
meet-in-the-middle search certifies
$d_{\rm DEM\,mech}\ge4$ for both.  Explicit XOR-validated weight-12 data-logical
witnesses give $d_{\rm DEM\,mech}\le12$ for both.  The model counts one
positive-probability instruction in the flattened
`decompose_errors=False` detector error model as one mechanism; it is not
silently promoted to exact physical-location circuit distance.  The bounds are
identical and therefore *establish* no circuit-distance separation; they do
not *exclude* one (the true mechanism distances could be anywhere in
$[4,12]$ independently).  They do rule
out any one-, two-, or three-mechanism undetectable logical hook in either
circuit.  The PBB DEM is denser (82,800 mechanisms, maximum detector degree 14)
than Gross (67,032, degree 9), but that structural difference is not a causal
explanation for the LER or latency gaps.

---

## 8. What this does *not* show

* It does not show that **no** non-CSS qLDPC code can beat CSS.  For PBB it
  proves parent domination only when $\delta=0$ (213/368 catalogue rows,
  including all 14 headline $[[144,12,12]]$ members).  Domination is **false in
  general** for $\delta>0$: two certified perturbations increase distance from
  4 to 6 while reducing $k$ from 8 to 4.  Those two are still dominated by the
  CSS $[[72,12,6]]$, but that is an envelope observation, not the parent
  theorem.
* It does not establish pointwise statistical separation of the logical error
  rates.  The predeclared familywise test failed; what we have is complete
  rank separation among five sampled schedules per code plus a mean
  difference.  EXP-025's schedule-optimised comparison is reported separately
  when complete.
* The one-ancilla depth separation is basis-independent (Theorem C6), but it
  does not rule out flag-qubit or cat-state gadgets trading ancillas and
  preparation/verification depth against data-coupling depth.  EXP-024 tests
  one explicit two-ancilla 4+4 construction; it is not a theorem over all
  multi-ancilla gadgets.
* The isolated decoder comparison uses matched BP+OSD settings, not the
  original artifact's full decoder configuration.  BP never converges at
  those settings, so the observed 6–17× gap (p50 10.3×, p95 16.6×, p99 6.3×)
  has no resolved cause.
  EXP-026 tests a predeclared co-design grid and two-stage shadow decoder.
* Circuit distance is bounded, not resolved exactly:
  $4\le d_{\rm DEM\,mech}\le12$ for both $n=144$ circuits in the explicitly
  stated undecomposed-DEM mechanism model.  No hook-induced distance
  degradation was demonstrated.
* Two-sector hooks are a structural difference in the fault set, not a
  demonstrated cause of the logical-error gap.
* A real rotated-surface-code baseline was simulated, but it is a **mapped**
  baseline: twelve independent $d=12$ patches, Stim/PyMatching conventions,
  and the same numeric $p$, not identical noisy locations or decoder.  It
  closes the matched-$k$ resource/reliability table but not a same-noise
  dominance claim.
* Track B was executed on $(12,12)$ and $(15,12)$ and found parent domination
  throughout its declared weight-$\le2$ held-out sample.  That does not exhaust
  higher-weight perturbations or controlled generalisations beyond PBB.
  Rare-event estimation and hardware-specific mapping were not executed;
  ordinary Monte Carlo sufficed at the measured $10^{-3}$–$10^{-1}$ rates.

## 9. Adversarial audit

| audit question | answer |
|---|---|
| Decoder estimate used as exact distance? | No. Code distances called exact have complete `OPTIMAL` sector searches and checked witnesses. EXP-029 reports only $4\le d_{\rm DEM\,mech}\le12$ in its named mechanism model. |
| Bound direction reversed? | **Yes, once — caught and corrected.** Legacy EXP-012 used catalogue `d_reported` as a PBB lower bound; its "weaker in 5" tally is relabelled a superseded diagnostic (§5.2). EXP-027's strict protocol treats catalogue $d$ only as an upper bound; its two reversals satisfy both certified directions and were rerun. |
| Code decomposable / duplicate? | BB references are single-component. EXP-032 finds 368/368 distinct only under its explicitly finite proven symmetry group. For the benchmark target, EXP-034 proves exactly (Theorem E): genuinely non-CSS under row operations, qubit permutations, local-H, and arbitrary independent local Cliffords — via a persisted 20-row GF(2) XOR contradiction with CP-SAT corroboration — and row-space **indecomposable** under that group. Catalogue-wide classification stays open and is labelled so. |
| Different noise models compared? | The Gross/PBB EXP-016 comparison uses one model identically. EXP-030's surface row deliberately uses Stim/PyMatching conventions and is labelled **mapped, non-identical-noise**, not a dominance test. |
| Ancillas or idles omitted? | One-ancilla Gross/PBB comparison counts 144 ancillas each and inserts idle noise at every schedule layer. Multi-ancilla EXP-024 has separate resource accounting. |
| Pseudo-threshold called a threshold? | No threshold is claimed. |
| Tuned on the test set? | EXP-016 settings and samples were fixed before decoding; EXP-025/026 have separate pre-registered development/held-out partitions. |
| $X$/$Z$/$Y$ correlations ignored? | Gross/PBB uses the undecomposed hypergraph DEM. The mapped surface baseline uses graphlike PyMatching and says so. |
| Zero failures reported as zero rate? | No. EXP-030's two zero-failure points are explicit 95% Clopper–Pearson upper bounds. |
| Latency compared under different load? | Caught and fixed: first timings quarantined (FR-004), re-measured on an idle machine with interleaved blocks and warm-up. No new latency is accepted under the current shared load. |
| Circuit held fixed across noise points? | Caught and fixed: EXP-007 did not (FR-007); it is downgraded. EXP-016 pins and persists every slot map. |
| Sample extrema reported as global optima? | No. EXP-016 extrema are sampled; the Bonferroni pointwise test failed. Schedule optimisation is a separate EXP-025 protocol. |
| Schedule class silently restricted? | **Yes, three times — all now corrected.** FR-005 (unschedulability was TI-only), FR-014/EXP-033 ($T_{\rm TI}>13$ vs $T_{\rm unrestricted}=9$), and the 2026-08-13 audit found Theorem C3's $T{=}6$ `INFEASIBLE` and EXP-016's "exhaustive" enumeration are also TI-class statements; both are now labelled as such. |
| Mechanism asserted without measurement? | BP convergence was measured and refuted as the latency cause (0% both); two-sector hooks and DEM density remain non-causal structural observations. |
| Partial run overwrote canonical evidence? | Caught repeatedly (FR-012). EXP-021/023/024/027/028/029 and both EXP-033 writers now route incomplete runs to `results/partial_runs/` and failures to `results/quarantine/`; checksum/path regressions cover the guards. EXP-032 alone wrote ungated, but only after a complete run and was hash-audited after the fact. |
| Own claims falsified or downgraded? | Twenty-one failed-route entries are retained with corrected statements; FR-013/014 test semantics, FR-015 coverage wording, FR-016/017 schedule-class/bound-direction defects, FR-018/019 the non-executable decoder/schedule grids, FR-020 the dead production gate on EXP-024's CLI, and FR-021 the unbound EXP-033 noiseless evidence. |

---

## 10. Artifacts

* Proofs — `proofs/pbb_structure.md`
* Hypothesis ledger — `notes/hypotheses.csv` (14 entries: 10 confirmed,
  3 falsified, 1 provisionally supported; regenerated from
  `experiments/write_hypotheses.py`)
* Failed routes — `notes/failed_routes.md` (21 entries)
* Experiments — `experiments/exp001`…`exp035` (numbered campaign scripts;
  canonical, partial, and quarantined outputs separated)
* Certified distances and witnesses — `results/certificates/`
* Raw and processed data — `results/raw/`, `results/processed/`
* Source and artifact provenance — `sources/manifest.yaml`,
  `third_party/manifest.yaml`, `artifacts/manifest.yaml`
