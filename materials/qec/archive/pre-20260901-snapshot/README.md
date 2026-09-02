# qec-codesign

Co-design study of quantum LDPC codes, fault-tolerant syndrome-extraction
circuits and decoders under circuit-level noise, centred on the non-CSS
**perturbed bivariate-bicycle** (PBB) construction of arXiv:2606.02418 and the
CSS bivariate-bicycle baselines of arXiv:2308.07915.

## Headline result

For the **213 of 368** published PBB codes with $\delta=0$ — including every
member of the headline $[[144,12,12]]$ family — perturbing a bivariate-bicycle
code provably buys nothing in $[[n,k,d]]$.  For all $\delta$ it provably cannot
raise $k$.  Parent domination for the remaining 155 codes is **false in
general**: **seven certified reversals** are now known (two at $n=72$, four at
$n=108$, one at $n=144$), and the EXP-036 SAT campaign is closing the last
undecided rows (41 of 87 still open as of 2026-08-16, all at $n\in\{180,360\}$).
No circuit-level result favoured the PBB candidate: depth, gates, hook
structure, sampled LER, and latency moved against it; the certified
mechanism-distance intervals tied.

* **[proved]** $k_{\rm PBB} = k_{\rm BB} - \delta$ with $\delta \ge 0$: the
  perturbation can only *remove* logical qubits.
* **[proved]** For $\delta = 0$ the nontrivial pure-$Z$ logical operators of the
  PBB code are *identical* to the parent CSS code's, so
  $d_{\rm PBB} \le d({\rm parent})$ and the parent weakly dominates in
  $[[n,k,d]]$. Covers **213 / 368** published PBB codes and **all 14** members
  of the headline $[[144,12,12]]$ family.
* **[proved, weaker]** For *arbitrary* $\delta$ there is a same-$(n,k)$ CSS
  "shadow" whose **$Z$-distance** caps $d_{\rm PBB}$ — a ceiling, not domination.
* **[exact, strict-provenance audit — EXP-027]** Over all 155 $\delta>0$
  catalogue codes: the parent provably dominates **66**; **2 are certified
  reversals** — `phase2_58`/`phase2_60`, where the $\delta=4$ perturbation of a
  $[[72,8,4]]$ parent yields a genuinely non-CSS $[[72,4,6]]$, i.e. **the
  perturbation strictly increased distance** (both distances exact/OPTIMAL,
  double-verified); **87 undecided** (parent distances uncertifiable at budget;
  literature values were *not* accepted as certificates). The two reversal
  codes are still dominated by the CSS $[[72,12,6]]$ at the same $n$ (more $k$,
  same $d$), so the CSS envelope stands; but "perturbations never buy
  distance" is now **refuted** for $\delta>0$. The $\delta=0$ domination
  theorem (Corollary 1) is untouched.
* **[exact, SAT closure — EXP-036]** The CDCL method behind Theorem F turned
  on EXP-027's 87 undecided rows: **46 decided so far** — 41 parent
  dominations and **5 new certified reversals** beyond `phase2_58`/`60`:
  `phase2_71`, `phase2_72` ($\delta=6$: $[[108,12,\le4]]$ parents,
  $[[108,6,>4]]$ PBBs), `phase2_88`, `9_6_0183` ($\delta=2$:
  $[[108,4,\le8]]$ parent, $[[108,2,>8]]$ PBB with witness 10), and
  `12_6_0217` ($\delta=4$: $[[144,8,\le8]]$ parent, $[[144,4,>8]]$ PBB with
  witness 10).  Every verdict is a verified-witness upper bound plus a
  replayed CDCL UNSAT lower bound under the frozen tie-safe comparison
  protocol, with CNF-hash-bound replay stamps
  (`experiments/exp036_delta_closure.py`; per-row identity-gated records;
  `verify` re-proves each decisive UNSAT from rebuilt matrices).  All 29
  $n=108/144$ rows are closed plus 17 at $n=180/360$; the remaining 41 run
  under a budget-laddered sharded schedule. A sound one-clause **translation symmetry break** (BB translations are code automorphisms acting transitively on each block, so any solution can be translated to anchor support at block index 0) cuts decisive UNSAT proofs by **6.2x** (77.8 s -> 12.5 s on the $[[144,12,12]]$ cap-11 instance); premise machine-checked in `tests/test_translation_symmetry_break.py`.  Sector decomposition of the nontriviality disjunction was exact but measured **8x slower** and was abandoned (FR-022).
* **[proved + machine-verified — EXP-038, Theorem G]** *Why* perturbation almost never
  buys distance, and a domination test that needs no search on the PBB.  The PBB's pure-$Z$
  centralizer is exactly the parent's $\ker[A\,B]$, and its pure-$Z$ stabilizers are the
  parent's plus the dressing space $\Delta=\{\lambda[C\,D]:\lambda[A\,B]=0\}$, with
  $k_Q=k_P-\dim\Delta$ **proved** (from $\operatorname{rank}H_Q=\operatorname{rank}H_P+\dim\Delta$)
  and verified with zero mismatches on all 368 rows.  So a surviving minimum-weight parent
  $Z$-logical certifies $d_Q\le d_Z(\mathrm{parent})$ outright, the survivor itself being the
  certificate; contrapositively $d_Q>d_Z(P)$ forces $k_Q\le k_P-t$ for $t$ the rank of the
  minimum-weight translation orbit.  This explains the confinement of **all seven** certified
  reversals to the $k$-halving class (58 of 155 rows; none among the 97 at $k_Q/k_P\in\{3/4,5/6\}$),
  each at margin **exactly zero** — they buy the smallest distance step parity allows ($+2$) at a
* **[proved + exact computation — EXP-039, Theorem H]** The quantifier move from one perturbation
  to an entire family.  $\Delta$ is an $R$-*submodule* of the $Z$-sector (the coefficient set
  $\{\lambda:\lambda[A\,B]=0\}$ is an ideal), so absorbing one minimum-weight parent logical absorbs
  its whole translation orbit.  Hence, with $T(P)=\dim\,(M(P)+S_Z)/S_Z$ the orbit-span dimension of
  *all* minimum-weight $Z$-logicals, computed *exactly* by a SAT enumeration terminating in certified
  UNSAT: $d_Q>d_Z(P)\Rightarrow k_Q\le k_P-T(P)$ for **every** perturbation $[C\,D]$ of $P$.
  $T=k_P$ closes the family outright.  Machine-certified: exact $T$ for **134/202** catalogue parents
  (exhaustively all $117$ at $n\le144$), **61 family-closed** including the Gross code itself
  ($T=12=k_P$) and **all eleven** catalogue $[[144,12,12]]$ parents; $249/368$ rows capped a priori
  with no per-row search.  All **7/7** certified reversals sit exactly on the bound
  ($k_Q=k_P-T$, slack zero, hypotheses replay-certified), and every reversal is CSS-dominated at equal
  $n$ against $162$ certified-exact CSS candidates (`exp036_envelope_check.json`, schema v2).
  **Theorem I (forced saturation, EXP-039/040).**  $\Delta$ is the image of the
  left kernel $L=\{\lambda:\lambda[A\,B]=0\}$, and $k_P=2\dim L$ on every BB parent
  (elementary rank-nullity; verified on all 202).  So $\dim\bar\Delta\le k_P/2$ always,
  and whenever $T\ge k_P/2$ — **133 of 134 certified parents** — Theorems H+I sandwich:
  $$d_Q>d_Z(P)\ \Rightarrow\ k_Q=k_P-T(P)\ \text{exactly}.$$  The trade is not just
  bounded, it is *priced to the qubit*: every distance increase pays precisely $T$
  logical qubits.  All 7 catalogue reversals and all 496 small-lattice strict
  increases (26,898 perturbations, 64 parents, zero upper-law violations, containment
  $\bar M=\bar\Delta$ verified) sit on the law.  Proofs: `proofs/pbb_structure.md`
  (Theorem H + I sections); tests: `tests/test_pbb_theorems.py`.
  Repository: `src/qec_research/codes/pbb_nogo.py`, `experiments/exp039_nogo_module.py`,
  `tests/test_pbb_nogo.py`; paper: `reports/paper_pbb_nogo.md` (claims locked by
  `tests/test_paper_claims.py`).  $n=180/360$ sweep continues; figures are monotone lower bounds.

* **[exact, in flight — EXP-037]** The programme's central question, asked directly and catalogue-wide: does *any* of the 368 published PBB codes escape the CSS BB envelope at its own length?  A row $[[n,k,\cdot]]$ is **envelope-dominated** when some same-length CSS BB code has $k_C\ge k$ and certified exact $d_C\ge U$, where $U$ is a verified-witness upper bound on the PBB distance.  The asymmetry makes the sweep affordable — cheap SAT witness on the PBB side, expensive UNSAT spent once per CSS code and amortised over every row at that length.  **250 of 368 rows certified envelope-dominated, zero genuine escapees so far**; the 22 flagged rows simply have no certified same-length CSS code with $k_C\ge k$ yet.  Both domination dimensions are re-derived from rebuilt algebra at classification time.
* **[proved + exact]** For every catalogue PBB $[[144,12,12]]$ member, **every
  generating set** of the stabiliser group contains a generator of symplectic
  weight $\ge 8$ (EXP-023: no element of weight $\le 7$ has a nonzero $X$-part,
  CP-SAT INFEASIBLE certificates; the pure-$Z$ subgroup has rank exactly 66 <
  132). Hence $\ge 8$ two-qubit layers for **any** one-ancilla schedule —
  the depth separation is basis-independent, a property of the code.
  The Gross code's optimal max generator weight is exactly 6 (rank $V_5=0$),
  and it achieves **7** layers; 7 is minimal for its published generators
  within the translation-invariant schedule class (unrestricted 6-vs-7 is
  closed only externally, by ASC arXiv:2603.21499).
* **[exact, $n=72$ pair]** **58.3 %** of single-ancilla hooks in the PBB
  `phase2_47` $[[72,12,6]]$ circuit corrupt
  *both* Pauli sectors; **0 %** do in the CSS circuit — a fault mode absent by
  construction from CSS. A structural difference; **not** a demonstrated cause
  of the logical-error gap (no fault-injection study was run).
* **[statistical]** Schedule-controlled matched simulation (minimum-depth
  **translation-invariant** schedules **exhaustively enumerated**, five
  sampled uniformly per code, slot
  maps pinned): all five PBB point estimates exceeded all five CSS ones
  (Mann-Whitney $p=0.004$), schedule-averaged ratio **3.55×** (Welch
  $p=0.014$).  The predeclared Bonferroni pointwise test **failed** — we report
  that, and make no best-schedule claim.
* **[statistical, isolated]** Decoder latency on an idle machine:
  **10.3× (p50)**, **16.6× (p95)**, **6.3× (p99)** worse for the PBB code
  (6–17× across the three percentiles).  Our
  proposed explanation (worse BP convergence) was **measured and refuted** — BP
  converges 0 % of the time for *both* codes — so the cause is unresolved.
* **[proved + exact, target-only — EXP-034]** The benchmark member
  `12_6_0193` is **genuinely non-CSS**: not CSS after arbitrary stabilizer
  row operations (dim $S_X$ + dim $S_Z$ = 86 < 132 = dim $S$), after any
  qubit permutation, after any local-Hadamard assignment (affine GF(2)
  system infeasible, certificate $0=1$), **and after any independent
  per-qubit local Clifford**: nineteen parity masks plus one one-hot parity
  row XOR to the GF(2) contradiction $0=1$, a solver-free linear
  certificate persisted and machine-rechecked (CP-SAT independently agrees,
  INFEASIBLE).  Because CSS-ness is permutation-invariant and permutation
  conjugates of local Cliffords are local Cliffords, the refutation covers
  the full row-operation × qubit-permutation × local-Clifford group.  Its
  row space is also **indecomposable** under that same group.  The
  $[[5,1,3]]$ control stays linear-feasible, so the certificate is not
  vacuous — the target's linear infeasibility is a strictly special
  structural property.
* **[proved + exact — EXP-035]** The benchmark target `12_6_0193`
  ($[[144,12,d]]$, non-CSS) has **independently certified exact distance
  $d = 12$**.  Lower bound: exact meet-in-the-middle exclusion of all
  weight-$\le5$ centralizer elements; a **complete weight-6 classification**
  (the triple-triple syndrome MITM enumerates exactly 72 weight-6
  zero-syndrome vectors, re-verifies each, and proves the set equals the 72
  stored pure-$Z$ checks, so no weight-6 logical exists); and weight
  $7$–$11$ excluded by **four UNSAT proofs** — one per orbit-representative
  sector of the machine-checked translation-orbit reduction of the
  $2^{24}-1$ logical classes, whose transported orbit union spans the full
  24-dimensional dual quotient (rank 24 on two independent GF(2) paths).
  Upper bound: an independently re-verified weight-12 witness.  The four
  sector proofs came from the **PySAT/CaDiCaL** backend (equisatisfiable
  CNF; totalizer cardinality + chained XOR) in **502/674/422/576 s**
  single-threaded — after OR-Tools CP-SAT, racing the same frozen sectors
  with 7 workers and multi-hour budgets, produced no proof (sector 0
  UNKNOWN across days).  The collector is fail-closed: frozen
  backend/status proof pairs, mandatory re-verification of any claimed
  solution through both GF(2) paths, and cross-backend contradiction
  fail-stops.  `results/certificates/pbb_12_6_0193_distance.json`.

### Findings we did not expect

1. **Schedule choice alone moves the logical error rate by 2.8× within each
   code** — comparable to the 3.55× between the two code families. Two valid,
   minimum-depth, formally verified schedules for the same code are not
   interchangeable. Any qLDPC co-design claim that does not pin and publish its
   schedule is under-specified. Our first sweep did not hold the schedule fixed
   and was downgraded (`notes/failed_routes.md` FR-007).
2. **Translation invariance can cost more than four syndrome layers.** Three
   weight-9 PBB $[[144,12,12]]$ members have no translation-invariant schedule
   through depth 13 (certified), yet EXP-033 finds an unrestricted depth-9
   `OPTIMAL` schedule for every one: $T_{\rm unrestricted}=9$ versus
   $T_{\rm TI}>13$. A symmetry ansatz that is convenient computationally can
   hide the actual circuit frontier.
3. **CDCL demolished CP-SAT on symplectic weight-bounded feasibility.** The
   four weight-$\le11$ orbit sectors of the `12_6_0193` distance problem are
   XOR-heavy Boolean systems (132 stabilizer parities + 1 detector parity +
   a cardinality cap over 144 OR-variables).  OR-Tools CP-SAT with native
   XOR constraints and 7 workers made no proof in multi-hour budgets
   (sector 0: days of UNKNOWN); CaDiCaL via PySAT on an equisatisfiable CNF
   (totalizer cardinality + chained XOR auxiliaries) proved all four UNSAT
   in 7–12 minutes each, single-threaded.  Exact qLDPC distance
   certification at this scale is a SAT problem, not a CP problem.

Full argument: [`reports/technical_report.md`](reports/technical_report.md).
Proofs: [`proofs/pbb_structure.md`](proofs/pbb_structure.md).

## A by-product worth its own line

The correctness criterion for one-ancilla syndrome extraction — *for every pair
of checks, the number of anticommuting shared qubits at which the first acts
before the second must be even* — **forces** depth 7 for weight-6 BB codes
within the **translation-invariant schedule class** (TI depth 6 is proven
`INFEASIBLE`; our artifacts do not exclude an unrestricted depth-6 schedule —
ASC, arXiv:2603.21499, independently certifies that exclusion). Within that
class this derives the depth-7 syndrome cycle of arXiv:2308.07915 from first
principles instead of assuming it.

## Layout

```
src/qec_research/
  gf2/            exact GF(2) linear algebra, two independent implementations
  symplectic/     stabilizer codes, centralizers, logicals, direct-sum detection
  codes/          BB and PBB construction; PBB structure theory (delta, shadows)
  distance/       certified exact distance (CP-SAT) and sector-restricted variants
  circuits/       mixed-stabilizer circuit compiler, parity-aware scheduler
  decoders/       BP+OSD on the undecomposed hypergraph DEM; parallel harness
experiments/      exp001 .. exp033, config-driven; canonical/partial outputs separated
tests/            validation gates and regression guards for falsified routes
proofs/           theorem statements with proofs and epistemic labels
notes/            open status, novelty matrix, hypothesis ledger, failed routes
results/          raw/, processed/, certificates/, partial_runs/, quarantine/
third_party/      pinned upstream artifacts (see third_party/manifest.yaml)
artifacts/        hashes and provenance for canonical research outputs
```

## Reproducing

```bash
uv sync --locked --extra dev
export PYTHONPATH=src

# validation gates
.venv/bin/python -m pytest tests/ -q
.venv/bin/python experiments/exp021_verify_equations.py
.venv/bin/python experiments/write_verified_results.py

# certified exact distance of the Gross code (about 10 min on 12 cores)
.venv/bin/python experiments/exp001_exact_distance_bb.py '[[144,12,12]]' 12 5400 12

# the delta=0 domination theorem over the whole catalogue (seconds)
.venv/bin/python experiments/exp006_parent_domination.py

# certified syndrome-extraction depth, BB vs PBB (translation-invariant model)
.venv/bin/python experiments/exp004_circuit_cost_survey.py 90

# schedule-controlled matched circuit-level benchmark (primary; supersedes exp007)
.venv/bin/python experiments/exp016_schedule_controlled.py 8000 0.002 12 26 5 20260811

# exploratory only, superseded by exp016 (FR-007):
# .venv/bin/python experiments/exp007_matched_benchmark.py 40000 0.0015,0.002,0.003 12 26 200

# decoder latency -- REQUIRES AN IDLE MACHINE, self-refuses otherwise
.venv/bin/python experiments/exp013_isolated_latency.py 300 12 0.002 4 4.0
```

The EXP-024 two-ancilla cat benchmark has **no canonical run**: its normal CLI
fails closed before any schedule or circuit work until the frozen supervised
v2 driver exists (FR-020).  Only `--structural-only` and `--validate-only` are
executable; `results/processed/exp024_cat_extraction.json` is intentionally a
zero-byte placeholder.

```bash
# EXP-034 target equivalence certificate (seconds; algebra only)
.venv/bin/python experiments/exp034_target_equivalence.py

# EXP-035 independent target-distance bounds; sectors resume and are load-gated
.venv/bin/python experiments/exp035_pbb_exact_distance.py assemble --seed 2026081934
```

Upstream artifacts are cloned under `third_party/` at the commits recorded in
`third_party/manifest.yaml` and are never modified.

## Epistemic conventions

Every claim in this repository is labelled **[proved]**, **[exact computation]**
(GF(2) arithmetic, or a solver returning proven `OPTIMAL`/`INFEASIBLE` — with
its schedule class or search scope named), **[certified bound]**,
**[statistical]** (with an interval), **[external]** or **[conjecture]**.  A
decoder-derived distance is never called exact.  A run with zero observed
failures is never reported as zero error rate.

Twenty-one failed-route entries from our own workflow are recorded, with their
corrections and regression guards, in
[`notes/failed_routes.md`](notes/failed_routes.md).
