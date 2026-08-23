# qec-codesign

Co-design study of quantum LDPC codes, fault-tolerant syndrome-extraction
circuits and decoders under circuit-level noise, centred on the non-CSS
**perturbed bivariate-bicycle** (PBB) construction of arXiv:2606.02418 and the
CSS bivariate-bicycle baselines of arXiv:2308.07915.

## Headline result

Two proved layers, one on each Pauli sector, plus a consistent circuit-level
negative.

**Z-sector (rate–distance).** The perturbation degree of freedom is priced
exactly: with $T(P)$ the orbit-span dimension of all minimum-weight parent
$Z$-logicals, $d_Q>d_Z(P)\Rightarrow k_Q\le k_P-T(P)$ for *every* perturbation,
and $=$ whenever $T\ge k_P/2$ (Theorems G/H/I). $63$ of $202$ catalogue parents
are family-closed, including the Gross code and all $11$ catalogue
$[[144,12,12]]$ parents. Every row of the published catalogue through $n\le180$
($318/318$) is CSS-envelope-dominated; only $n=360$ rows remain open.

**X-sector (collapse mechanism) — fully classified.** The single-row syzygy
channel that silently destroys $X$-distance is governed by one ideal
$I=\operatorname{Ann}_R(A,B)$: $\dim S=2\dim I^\infty$, so a parent is
demote-full iff $I$ is nilpotent, demote-immune iff $I^2=I$, mixed iff neither
(**Theorem J-G**, EXP-053: reproduces the whole $202$-parent classification in
$8.9$ s, $I^2=0$ on all $192$ demoting parents, $I^2=I$ on all $10$ immune).
Consequences: **odd$\times$odd lattices cannot demote at all** (semisimple $R$;
in our seven-instance baseline table only $[[90,8,10]]$ on $(15,3)$ has this
protection; the wider literature contains many more odd-lattice BB codes), and
**no pair with $\operatorname{wt}(A),\operatorname{wt}(B)\le3$ is ever mixed**
(**Theorem J-I**), while weight $4$ attains it (**Theorem J-J**: $(2,3)$,
$A=(1+x)(y+y^2)$, $B=Ay$, $240/255$ classes demote). Census:
$653{,}022{,}021$ weight-$\le3$ pairs over $18$ lattices, zero mixed.

**Odd-lattice extension (EXP-055).** The prior claim that $[[90,8,10]]$ was the
only odd$\times$odd BB code in print was false and is withdrawn: $27$ sourced
instances now form the validation battery. The published principal-code
structure becomes, in our bar convention,
$J^2\cong\ker H_X/S_Z$ with $J=\bar{\operatorname{Ann}_{\rm left}(A,B)}$;
verified **27/27**. The pole ceiling has zero violations on the **five**
independently exact-certified literature rows; all 25 reproduced reported
values also pass as a non-certifying sanity check. Exhaustive sweep: all $65$
odd lattices with $\ell m\le180$,
$4{,}229{,}823{,}962$ weight-$\le3$ pairs, zero $k$ mismatches. Exact screen
waves produced connected/indecomposable $[[30,8,4]]$, $[[54,8,6]]$ and
$[[126,12,10]]$ BB references. The first is a new constructor within the
checked BB tables but **globally dominated** by Grassl's $[[30,8,7]]$; the
latter two independently exactify Wang--Mueller BP-OSD parameters. None is an
end-to-end Pareto result.
Fixed-point screen through $n=126$: **715** symmetry classes / **10,247**
normalised pairs; all **598/598** classes with a locally exact reference are
dominated (588 explicit witnesses + 10 CP-SAT), 117 high-$k$ classes have no
reference, zero survivors/undecided. The $k$ census remains complete through
$n=360$; no distance-closure claim beyond $n=126$.

**Circuit level.** No circuit-level result favoured the PBB candidate: depth
(basis-independent $\ge8$ vs Gross $7$), gates, hook structure, sampled LER and
decoder latency all moved against it; certified mechanism-distance intervals
tied. A fault-tolerant one-ancilla mixed-stabilizer circuit for the flagship
remains unconstructed *and* unrefuted — the open end-to-end item.

Full ordered directions: [`notes/next_breakthroughs.md`](notes/next_breakthroughs.md).

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

* **[proved + machine-exact — EXP-046…054, Theorems J-C/J-E/J-E″/J-G/J-H/J-I/J-J]**
  The $X$-sector story, end to end. Candidate monotonicity $d_X(Q)\ge d_X(P)$ is
  **false** (J.5 refuted: four machine-verified witnesses at full $k$, the
  flagship's own sibling is exactly $[[144,12,6]]$ — EXP-048 certifies $133/133$
  siblings terminal-UNSAT, $99$ strict drops, $98$ dimension-preserving). The
  mechanism is *decidable in polynomial time* (Theorem J-C, one GF(2) rank test)
  and *exactly characterized*: no single-row demotion $\iff
  1\in L_{\rm pre}+\operatorname{Ann}_R(M)$ (Theorem J-E, Nakayama), with the
  demote fixed set equal to the stable image $I^\infty M$ (Theorem J-E″, proved
  and independently audited). **Theorem J-G reduces all of it to the ideal
  $I=\operatorname{Ann}_R(A,B)$** — $\dim S=2\dim I^\infty$; nilpotent ⟺
  demote-full, idempotent ⟺ immune, otherwise mixed — verified on $202/202$
  parents against the independent module route, and sharpened to $I^2=0$
  (all $192$ demoting) / $I^2=I$ (all $10$ immune). **J-G1**: both $\ell,m$ odd
  ⟹ $R$ semisimple ⟹ *no* demotion possible ($3{,}600$ parents, zero
  violations; the published $[[90,8,10]]$ on $(15,3)$ is the only protected
  baseline). **J-H/J-I**: writing $G=G_2\times G_{\rm odd}$, exact vanishing at
  a local factor $F_\chi[G_2]$ costs two support points per occupied $2$-coset,
  so **no weight-$\le3$ pair is ever mixed on any lattice**; the single-coset
  predicate separates the catalogue exactly ($10$ vs $192$). **J-J**: weight $4$
  realizes mixed — $(2,3)$, $A=(1+x)(y+y^2)$, $B=Ay$: $k_P=8$, $\dim S=4$,
  $240/255$ classes demote, verified by ideal chain, module chain and brute
  force. Census: $653{,}022{,}021$ weight-$\le3$ pairs over $18$ lattices
  (all catalogue lattices plus the previously underexplored $(12,12)$,
  $(15,12)$), **zero mixed**. X-side partition of the catalogue: $192$
  demotion-realized / $8$ **certified** X-monotone (exhaustive one-shot CP-SAT
  channel exclusion, EXP-051) / $2$ X-undecided with verified weight-$8$
  multi-row channels ($=202$). Notes:
  `notes/theorem_je_exact_decision.md`, `notes/theorem_je2_demote_trichotomy.md`,
  `notes/theorem_jg_ideal_invariant.md`; tests `tests/test_exp052_*.py`,
  `tests/test_exp053_ideal_invariant.py`.
* **[published structure reproduced + new exhaustive computation — EXP-055]**
  Odd-lattice BB logicals have the exact reciprocal-pole transversal
  $J^2\cong\ker H_X/S_Z$ (Eberhardt–Steffan; $J=\bar I$ is load-bearing).
  Machine audit: $27/27$ isomorphisms; pole ceiling zero violations on four
  locally exact references; 65 lattices / $4.23\times10^9$ pairs / zero $k$
  mismatches. The Pareto screen excludes all BP-OSD/Monte-Carlo estimates,
  uses explicit reduced logical witnesses first, then capped CP-SAT.
  Artifacts: `results/processed/exp055_*.json`; proof/convention note:
  `notes/theorem_k_certified_ceiling.md`; tests:
  `tests/test_exp055_odd_lattice.py`.

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
