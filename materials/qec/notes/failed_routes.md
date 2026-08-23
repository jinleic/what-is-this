# Failed routes and falsified claims

Append-only.  Every entry records a route that was tried and abandoned, why it
failed, and the artifact it left behind.  Falsified *own* claims are recorded
here with the same prominence as external ones.

---

## FR-001 — Wrong nontriviality detector for pure-Z logical operators
**Date** 2026-08-11 · **Track** C (distance) · **Status** FALSIFIED BY OUR OWN GATE

**What we claimed (briefly).**  An early run of `exp002_pbb_ceiling_audit.py`
reported that six catalogue PBB `[[144,12,12]]` codes — including entries
marked `trust_level: EXACT` — possessed a weight-6 nontrivial logical operator,
apparently refuting their distance claim.

**Why it was wrong.**  To decide whether a vector `v` in the centraliser is a
*nontrivial* logical we require `<v, L_j>_s = 1` for some logical basis element
`L_j`, and

    <(0|v), (L_x|L_z)>_s  =  v . L_x .

The detector must therefore be the **x-part** of a logical basis element.  The
audit instead paired the pure-Z candidate against *Z-type* representatives with
a plain dot product.  That pairing is degenerate: a genuine element of the
stabilizer group can satisfy `rep . v = 1`, so the constraint admits stabilizer
elements as false witnesses.

**How it was caught.**  Gate 2 / Gate 4: the explicit witness was re-tested with
`StabilizerCode.is_logical`, which reported `lies in stabilizer group: True`,
i.e. `is_logical = False`.  The claim was retracted before it propagated.

**Correct result after the fix.**  The minimum pure-Z logical weight of PBB
`12_6_0187` is **12** (CP-SAT `OPTIMAL`), matching the catalogue.  The parent
BB code's minimum pure-Z logical weight is also **12** (CP-SAT `OPTIMAL`) —
which *confirmed* Theorem 1 rather than refuting the catalogue.

**Artifacts left behind.**
* `src/qec_research/distance/sectors.py` — a correct sector-restricted
  minimum-weight solver with the pitfall documented at the top of the file.
* Controls in the same run: `[[72,12,6]] -> 6`, `[[108,8,10]] -> 10`, both
  matching known exact distances.

**Lesson encoded.**  Any coset-nontriviality constraint must use the symplectic
pairing against the *opposite* Pauli sector, never a same-sector dot product.

---

## FR-002 — Arbitrary edge colouring is not a valid syndrome-extraction schedule
**Date** 2026-08-11 · **Track** A (circuits) · **Status** FALSIFIED BY STIM

**What we tried.**  Schedule the two-qubit gates of a one-ancilla-per-check
syndrome circuit by König edge-colouring the check/qubit incidence graph,
reasoning that conflict-freedom (no qubit or ancilla used twice per layer) was
the only requirement.  This attains the combinatorial optimum
`Delta = max(max check weight, max qubit degree)` layers.

**Why it was wrong.**  Controlled-Paulis from *different* ancillas onto the
*same* data qubit do not commute when the Paulis anticommute:

    C-X_{a,j} C-Z_{b,j}  =  C-Z_{b,j} C-X_{a,j} . CZ_{a,b}

Reordering to a canonical order therefore emits one `CZ` between the two
ancillas for each anticommuting shared qubit at which the order is inverted.
Because `CZ^2 = I`, the circuit measures the intended operators **iff**

    #{ j in J(a,b) : a acts before b }  is even    for every pair (a,b),

where `J(a,b)` is the set of shared qubits with anticommuting Paulis.  A plain
edge colouring does not enforce this.

**How it was caught.**  Gate 5 (circuit semantics).  A König-coloured Steane
circuit has **5** parity-violating pairs and fires **1402** detectors in 200
noiseless shots.

**Artifacts left behind.**
* `src/qec_research/circuits/scheduling.py` — the parity criterion, an
  independent verifier (`verify_schedule`), and a CP-SAT scheduler that
  enforces it.  Validated to **0** noiseless detector firings on Steane
  `[[7,1,3]]`, Shor `[[9,1,3]]` and the non-CSS perfect code `[[5,1,3]]`.
* A by-product: the parity constraint *forces* depth 7 for every weight-6 BB
  code (depth 6 proven `INFEASIBLE`), independently reproducing the depth-7
  cycle of arXiv:2308.07915 from first principles.

**Lesson encoded.**  Conflict-freedom is necessary but not sufficient; the
anticommuting-overlap parity rule is the actual correctness condition, and it
strictly raises the achievable depth.

---

## FR-003 — Over-broad CSS-shadow refutation campaign
**Date** 2026-08-11 · **Track** C · **Status** ABANDONED (superseded, not wrong)

**What we tried.**  `exp003_css_shadow_domination.py` attempted, for all 318
catalogue codes with `n <= 216`, to decide *both* `d_Z(Q') >= D` and
`d_X(Q') >= D` by CP-SAT refutation at cap `D-1`.

**Why it was abandoned.**  At `D = 12` the capped search is still large; every
one of the first 22 codes returned `undecided` after the 300 s budget, at
~27 CPU-minutes each.  Extrapolated cost was far beyond value.

**What replaced it.**  The `delta = 0` special case needs *no solver at all*:
when `delta = 0` the CSS shadow **is** the parent BB code `BB(A,B)`, and
Proposition 4 (verified on 368/368 catalogue entries and all 7 Bravyi
instances) gives `d_X = d_Z` for every BB code.  Domination then follows
algebraically.  `exp006_parent_domination.py` settles 213/368 codes in 7
seconds.

**Lesson encoded.**  Look for the algebraic special case before spending
solver time; a theorem that removes a search is worth more than a faster search.

---

## FR-004 — Decoder timing collected under shared machine load
**Date** 2026-08-11 · **Track** D · **Status** METHODOLOGICAL DEFECT, QUARANTINED

**What happened.**  `exp007` reports decoder latency percentiles and measures
the CSS baseline first and the PBB code second.  A concurrent 20-process
subagent job and a 12-worker CP-SAT job were running during part of the
collection, so the two codes were not exposed to identical machine load.

**Consequence.**  Logical-failure *counts* are unaffected (they depend only on
the sampled faults and the decoder's deterministic output).  All *timing*
numbers from that run are biased and must not be used for a decoder comparison.

**Remedy.**  Latency and throughput are re-measured in a dedicated isolated
run with no other job on the machine, fixed worker count, and a warm-up batch.
Only those numbers are reported as timing evidence.

**Lesson encoded.**  Never interleave a timed A/B comparison with unrelated
compute; record machine load alongside every latency figure.

---

## FR-005 — "Unschedulable" claim was an artifact of the schedule class
**Date** 2026-08-11 · **Track** A · **Status** OUR OWN CLAIM FALSIFIED BY OUR OWN EXPERIMENT

**What we claimed.**  Theorem C4, first draft: 83 of 318 catalogue PBB codes
admit *no* one-ancilla syndrome-extraction schedule, on the strength of CP-SAT
returning `INFEASIBLE` at every depth from the combinatorial lower bound up to
(number of monomial directions) + 2.

**Why it was wrong.**  Those searches were all run over the *translation-
invariant* model, in which one time slot is shared by every translate of a
monomial direction.  That is a strict sub-class of one-ancilla schedules.  The
infeasibility proofs were sound; the generalisation was not.

**How it was caught.**  `exp009_general_schedule_nogo.py` re-ran the three
smallest such codes on the unrestricted model (one time variable per edge).
All three are schedulable, with CP-SAT status `OPTIMAL`:

| code | n | translation-invariant | general model | noiseless detectors |
|---|---|---|---|---|
| `phase2_26` | 36 | INFEASIBLE at every depth | depth 7 | 0 / 4000 shots |
| `phase2_21` | 36 | INFEASIBLE at every depth | depth 9 | 0 / 4000 shots |
| `phase2_44` | 36 | INFEASIBLE at every depth | depth 9 | 0 / 4000 shots |

**What survives.**  The headline depth separation does not depend on the
schedule class at all.  Every catalogue PBB `[[144,12,12]]` has maximum check
weight >= 8, so Proposition C2 (an ancilla takes part in at most one two-qubit
gate per layer) forces depth >= 8 for *any* one-ancilla schedule, against the
Gross code's 7.  Theorem C4 was restated on that basis and the schedule-class
observation demoted to Proposition C5.

**Lesson encoded.**  An `INFEASIBLE` certificate is only as general as the
model it was proved in.  State the model in the theorem, and prefer a bound
that follows from a counting argument over one that follows from a solver.

---

## FR-006 — We nearly filed an erratum against a paper that was right
**Date** 2026-08-11 · **Track** C · **Status** OUR OWN CLAIM RETRACTED

**What we claimed.**  That arXiv:2606.02418 states the PBB commutation
condition as `A C^T + B D^T = 0`, which is strictly stronger than the correct
symmetry condition, and that this was a paper-versus-code discrepancy worth
recording as an erratum candidate.

**Why it was wrong.**  The `= 0` wording came from a *summary of the paper
produced inside our own pipeline*, not from the paper.  Reading the primary
source settles it; arXiv:2606.02418v1 Sec. III.2 says, verbatim:

> "All rows commute if and only if (A C^T + B D^T) mod 2 is symmetric over F_2
> (this is the only nontrivial commutativity condition; commutativity between
> block 1 and block 2 is automatic from BB ring commutativity)."

That is exactly the condition we derived and exactly what the reference
implementation tests.  There is no discrepancy.

**How it was caught.**  Fetching and reading the HTML of the primary source
before letting the claim reach the technical report's conclusions.

**What survives.**  Our independent derivation of the condition (Lemma 0) is
still worth stating, because every downstream structural result (Props 1-3,
Theorems 1-2) is built on the same symplectic computation.  It is now labelled
as agreeing with the source rather than correcting it.

**Lesson encoded.**  A summary of a source is not the source.  Never attribute
a claim to a paper without quoting the paper.  This applies with full force to
summaries produced by our own tooling.

---

## FR-007 — EXP-007 did not hold the circuit fixed
**Date** 2026-08-11 · **Track** A/D · **Status** RESULT DOWNGRADED, EXPERIMENT REPLACED

**What was wrong.**  `exp007_matched_benchmark.py` called `build()` separately
for the p=0 sanity check and again inside the loop over every noise strength.
Each call re-solved the scheduling CP-SAT model, which at the time had no fixed
random seed and ran multi-threaded (`num_workers=12`).  CP-SAT returns whichever
feasible schedule its workers reach first, so the circuit under test could
differ between noise points, and neither the slot map nor a schedule hash was
archived.

**Why it matters.**  Schedule choice is not a nuisance parameter.  Rebuilding
the *same* Gross code at p=0.002, 12 rounds, with the same decoder settings but
a schedule obtained from a different CP-SAT worker count gave

    81 / 8000  = 1.0125e-2      versus      203 / 40000 = 5.075e-3

with disjoint Clopper-Pearson intervals — a factor of two from schedule
ordering alone.  Gate ordering changes how ancilla faults propagate, so this is
expected in hindsight; it was not controlled.

**How it was caught.**  A cross-check run (`exp015`) reproduced the Gross point
under nominally identical settings and disagreed with the headline sweep.  The
same run also revealed that the ancilla-idle-noise strip was a no-op, so the
noise convention was *not* the explanation.

**Consequences.**
* EXP-007's logical error rates are downgraded to **arbitrary-schedule
  exploratory evidence** and are no longer the basis of any claim.
* The claim "this is the fair/best PBB circuit because its depth is minimal" is
  withdrawn: minimal depth does not single out a circuit.
* `cpsat_schedule` now takes `random_seed` and pins the solver to a
  deterministic single-worker search, so a seed reproduces a schedule exactly.
* `exp016_schedule_controlled.py` replaces it: one solve per (code, seed), the
  slot map hashed and persisted, the *same* map reused for the p=0 check and
  the noisy run, and the comparison decided by the **distribution of LER over
  several distinct minimum-depth schedules** per code rather than by one
  circuit.

**What is unaffected.**  Everything algebraic (Props 1-4, Theorems 1-2,
Corollary 1), the depth bounds (Prop C2, Theorems C3/C4 — these are about
existence and minimality, not about which schedule is returned), the
hook-spectrum enumeration, and the independent artifact-derived cross-check of
the depth-7 Gross circuit.

**Lesson encoded.**  If a solver picks part of the experimental apparatus, its
seed is part of the protocol.  Persist the artifact it returns, not just the
command that produced it.

---

## FR-008 — "BP fails to converge" did not explain the decoder slowdown
**Date** 2026-08-11 · **Track** D · **Status** OUR OWN EXPLANATION FALSIFIED

**What we claimed.**  EXP-013 measured the non-CSS PBB circuit decoding 10.3x
(p50) / 16.6x (p95) / 6.3x (p99) slower than the CSS Gross circuit at matched
settings, while its detector error model is only 1.24x larger.  We wrote that
the gap came from decoder *behaviour*: belief propagation converging less often
on the mixed-check Tanner graph, pushing more shots into the ordered-statistics
stage.  That was an assertion with no measurement behind it.

**How it was caught.**  `exp017_bp_convergence.py` reads
`BpOsdDecoder.converge` after every call.  At 30 min-sum iterations, BP
converges on **0.0 % of shots for BOTH codes**.  Every shot already falls
through to OSD in both cases, so differential BP convergence cannot explain
anything.

| | CSS Gross | non-CSS PBB | ratio |
|---|---|---|---|
| BP converged | 0.0 % | 0.0 % | 1.00 |
| DEM error mechanisms | 67 032 | 82 800 | 1.24 |
| detection events / shot | 119.2 | 163.4 | 1.37 |
| mean decode latency | 183.2 ms | 1982.9 ms | 10.8 |

**What survives.**  The latency measurement itself (isolated machine,
interleaved blocks, warm-up discarded) stands.  What is withdrawn is the causal
story.  All we can say is that the cost grows superlinearly in the 1.24x
column count and 1.37x detection-event count, and that we have not isolated
which term dominates.

**A caveat this exposed.**  Because BP never converges at our settings, the
comparison is effectively OSD-0 applied to the two detector error models.  It
is matched -- identical settings for both codes -- but it is *not* the decoder
configuration of the original artifact (OSD-CS order 7, far more BP
iterations).  A stronger BP stage could change the absolute numbers and
possibly the ratio.  This is now stated in the report.

**Lesson encoded.**  "The obvious mechanism" is a hypothesis, not a finding.
If the instrumentation to test it is one attribute read away, read it before
writing the sentence.

---

## FR-009 — "Exact" distances that only proved a lower bound
**Date** 2026-08-11 · **Track** C · **Status** OUR OWN CERTIFICATES INVALIDATED AND REGENERATED

**What was wrong.**  `exact_distance_css` and `_min_pure_sector` both seeded
`best = upper_bound` from the caller-supplied search cap and then solved each
logical sector with `upper_bound = best - 1`.  When every sector returned
`INFEASIBLE`, `best` was still sitting at the seed, `best_sol` was still
`None`, and the routines returned that seed with `exact = True`.

That establishes only

    d >= D        (nothing of weight <= D-1 exists)

and never

    d <= D        (no operator of weight D was ever exhibited).

EXP-001 was invoked with the catalogue value as the cap, so all four persisted
BB certificates claimed `exact: true` while containing **no witness at all** —
verified by grepping: zero of four files had a witness field.

**Why it slipped through.**  The numbers were right.  Every value agreed with
the published literature and with earlier uncapped runs that *had* produced
genuine witnesses, so nothing looked wrong on the surface.  The defect was in
what the artifact *proved*, not in what it said.

**Fix.**
* `upper_bound` is now documented and treated as a **search cap only**; `best`
  starts at `None` in both routines and per-sector caps use `best` rather than
  `best - 1`, so a witness at the optimal weight can still be found.
* `exact` requires all sectors decided **and** a witness whose weight equals
  the reported value.
* A found `best` is no longer reported as a lower bound while any sector is
  undecided (a lighter operator could live in the undecided sector);
  `d_*_is_upper_bound_only` flags that case.
* EXP-001 now re-verifies every witness through the **independent pure-Python
  bitset path** — weight, membership in the opposite-type kernel,
  non-membership in the stabilizer row space, and `is_logical` on the full
  symplectic vector — and persists the witness support in the certificate.
  The certificate field is `CERTIFIED_EXACT`, true only if all of that passes.
* Two regression tests: hand each solver a cap *below* the true distance and
  require it to report no value and `exact = False`.

**Regenerated results** (all with independently verified witnesses):
`[[72,12,6]] d=6`, `[[90,8,10]] d=10`, `[[108,8,10]] d=10`, `[[144,12,12]] d=12`.
Values unchanged; what changed is that they are now actually certified.

**Lesson encoded.**  A solver flag named `exact` must be earned by the two-sided
argument, not by the absence of a counterexample below a cap the caller chose.
Seeding an optimiser with the answer you expect is how you get it back.


---

## FR-010 — Wrong centralizer constraint in the held-out domination test
**Date** 2026-08-12 · **Track** B · **Status** UNSOUND RESULTS DELETED BEFORE PUBLICATION

**What was wrong.**  `exists_logical_within` built its feasibility model with
`constraints = code.H`.  Membership in the symplectic centralizer is
$(\Lambda H)v = 0$, not $Hv = 0$; `exact_distance_symplectic` has always used
`lambda_swap(H)`.  On CSS blocks the two can coincide sector by sector, which
is why a smoke test passed — but on **mixed** PBB stabilisers the wrong
condition admits vectors that are not logical operators, so a `DOMINATED`
verdict was not sound.

**Fix.**  `constraints = lambda_swap(code.H)`; every witness is now verified
with `code.is_logical(v)` and a symplectic-weight recheck before the verdict is
accepted, raising rather than reporting if either fails; and a regression test
cross-checks the query against `exact_distance_symplectic` on a genuinely mixed
non-CSS code.  The 12 already-written verdicts were **deleted, not patched**.

*Removal path, recorded because the ledger must not assert an action it did not
take.*  The first attempt bundled `rm` with this very ledger entry in one shell
command and was refused by tool policy, so neither happened.  The deletion was
then performed separately and its effect verified by listing
`results/**/exp01[89]*`, which returned no matches.  The regenerated artifacts
below are written by a fresh run, not by editing the old files.

**Lesson encoded.**  CSS structure hides non-CSS bugs.  Worse, my first attempt
at the regression test asserted `(x_j and z_j)` on a single qubit — that tests
for $Y$ support, not for a mixed stabiliser, and would have let a CSS-only code
pose as the non-CSS test bed.  The guard needs the right predicate: the row has
$X$ support **and** the row has $Z$ support.

---

## FR-011 — "Exhaustive" enumeration over a family containing non-trinomials
**Date** 2026-08-12 · **Track** B · **Status** STAGE-1 ARTIFACTS INVALIDATED AND REGENERATED

**What was wrong.**  `canonical_A` generated `[(a1,0), (0,a2), (0,a3)]` over all
`a1` and all `a2 < a3`.  When `a1 = a2 = 0` the monomial $(0,0)$ appears twice,
XORs to zero, and the "trinomial" is actually a **monomial** — a different
family, with correspondingly degenerate codes.  `canonical_B` had the mirror
defect at `b1 = b2 = 0`.

Affected: $(12,12)$ 11/792 of $A$ and 11/792 of $B$, i.e. **17 303 of 627 264
pairs (2.8 %)**; $(15,12)$ 11/990 and 14/1260, i.e. **27 566 of 1 247 400
(2.2 %)**.  Every stage-1 pair count and $k$-histogram was therefore computed
over a family that was not the one named in the write-up.

**Fix.**  Filter `len(set(terms)) == 3` in both generators; record the true
family sizes in the artifact; regenerate both stage-1 files from scratch.

**Measured impact, after regenerating.**  On $(12,12)$ the family drops from
$792\times792$ to $781\times781$ and the enumerated pairs from 627 264 to
609 961.  The 17 303 removed pairs fell **entirely in the $k=0$ bucket**
(557 632 $\to$ 540 329); every $k>0$ class is bit-identical
($k=4$: 40 960, $k=8$: 19 584, $k=12$: 576, $k=16$: 7 608, $k=20$: 16,
$k=24$: 144, $k=32$: 702, $k=48$: 4, $k=64$: 37, $k=128$: 1).  So the error was
real and the word "exhaustive" was wrong, but no QEC-relevant number moved.
That is the honest statement: the defect is reported because it existed, not
because it changed a conclusion.

**Second defect, same fix.**  Stage 1 stored parent *indices* into the
canonical family and EXP-019 resolved them with `canonical_A(...)[ia]`.  Once
the family changes size, every stored index silently points at a different
code.  Representatives are now stored as **terms**, so an artifact cannot be
reinterpreted against a regenerated family.

**Lesson encoded.**  "Exhaustive over family F" is a claim about F, and it has
to be tested — a generator that emits members outside its own stated form makes
the count, the histogram, and the word *exhaustive* all wrong at once.  Never
key a saved artifact to a position in a list that a later edit can renumber.

---

## FR-012 - Canonical artifact writes unprotected against partial/failing runs (three successive defects)
**Date** 2026-08-12 · **Track** infrastructure (EXP-021) · **Status** FIXED; gate exercised in all three directions

**What was wrong.** (1) exp021 wrote its artifact unconditionally: a failing 60-row
smoke run overwrote the full-catalogue 368-row result. (2) The first fix gated on
pass/fail only; after E15/E16 were widened to scan the whole catalogue, a *passing*
60-row prefix would still have clobbered the canonical artifact. (3) The second fix
routed failing runs to `results/partial_runs/` while the report claimed they went to
`results/quarantine/` - the quarantine branch was unreachable.

**Fix.** Routing rule extracted to `src/qec_research/artifacts.py::canonical_route`
(canonical iff clean AND full coverage; clean-partial -> `results/partial_runs/`;
any failing run -> `results/quarantine/`), unit-tested over all four combinations
(`tests/test_conventions.py::test_artifact_routing`), and exercised end-to-end:
injected-failure run -> quarantine with canonical sha256 unchanged and exit 1;
passing 60-row run -> partial_runs with canonical sha256 unchanged; full 368-row
run -> canonical rewritten.

**Lesson.** A guard is verified by attempting the exact clobber it must prevent -
reading the guard's code is not verification. All three defects survived reviews
that looked at the diff but did not exercise the failure path.

---

## FR-013 - Verdict-polarity bug: counterexample path emitted POSITIVE
**Date** 2026-08-12 · **Track** infrastructure (EXP-031) · **Status** FIXED + regression-tested

**What was wrong.** exp031's mismatch branch set `verdict = "POSITIVE"` with a
reason string describing a counterexample. Root cause: the parent task spec
itself said 'any mismatch => POSITIVE-as-counterexample', conflating 'we learned
something' with the machine-readable field downstream integration keys on. Both
outcomes mapping to POSITIVE made the one decisive bit unreadable. The shipped
zero-mismatch artifact never exercised the branch, so no artifact is affected.

**Fix.** Verdict decision extracted to `decide_verdict()` (pure, importable):
mismatch -> NEGATIVE ('criterion FALSIFIED...'), zero-mismatch on >=150 decided
-> POSITIVE, else INCONCLUSIVE; counterexamples dominate small samples.
Regression tests in `tests/test_exp031_verdict.py` exercise all four paths,
including the exact historical bug (mismatch != POSITIVE).

**Lesson.** Verdict taxonomies are contracts: 'valuable outcome' and 'conjecture
status' are different axes and must be different fields. The spec author owns
spec-induced bugs.

---

## FR-014 - Certified refutations misfiled as 'undecided': the two-value law was actually falsified
**Date** 2026-08-12 · **Track** A/infrastructure (EXP-031) · **Status** FIXED; verdict corrected to NEGATIVE; regression-tested

**What was wrong.** exp031 built mismatch records only from `decided` instances
(those with a certified minimum). Three PBB members (pbb144-01/-03/-11, all
w_max = 9) carried certified INFEASIBLE at every depth 9..13 in the
translation-invariant class - the criterion's prediction (depth 10) and the
two-value law were certifiably REFUTED for them - yet they were filed as
'undecided' and the run reported 208/208 matches, verdict POSITIVE.
A certified infeasibility at the predicted depth refutes the prediction
regardless of where the true minimum lies.

**Fix.** `classify_prediction()` now three-values both laws: certified
INFEASIBLE at the predicted depth => matches=False (mismatch) even when
undecided; INFEASIBLE at both w_max and w_max+1 => two-value law False.
Timeouts/UNKNOWN never refute. Mismatch records are self-describing
(trace statuses, class, refutation flags). Verdict: NEGATIVE - 'criterion
FALSIFIED: 3 instances'. Regression tests cover refuted-undecided, UNKNOWN-
stays-undecided, and decided paths (tests/test_exp031_verdict.py, 8 pass).

**Scientific consequence (the good kind).** The corrected statement is
class-scoped: the criterion + two-value law hold on all 208 decided instances
(47 BB TI, 11/14 PBB TI, 150 random unrestricted), and are REFUTED in the
TI class for 3 weight-9 PBB members, which admit NO translation-invariant
schedule at any depth up to 13. Whether the law survives unrestricted for
those members is decided by EXP-033.

**Follow-up (EXP-033, 2026-08-13).**  The law survives this adversarial trio
once the symmetry restriction is removed: all three have verified unrestricted
depth-9 schedules, exactly attaining $w_{\max}=9$, while $T_{\rm TI}>13$.
Thus the corrected scientific result is a large schedule-class separation, not
an unrestricted counterexample.

**Lesson.** 'Undecided' must mean 'no certificate either way'. Any bucket that
can swallow a certified refutation converts falsification into silence.


---

## FR-015 - Coverage sentence multiplied aggregate counts and implied untested lattices
**Date** 2026-08-13 · **Track** B/infrastructure (EXP-028) · **Status** FIXED; metadata corrected; regression-tested

**What was wrong.**  EXP-028 correctly stored 247 total solver-tested
perturbations distributed across three reached parents, all on $(6,6)$.  Its
machine-readable `coverage_statement` rendered this as “3 parents x 247
perturbations on (6,6),(6,4),(4,6)”.  That wording could be read as 741 tests
and as positive coverage of the two lattices whose counters were exactly zero.
The per-lattice records and every solver result were correct; the aggregate
sentence was not.

**Fix.**  Coverage formatting is now a pure `format_coverage_statement()`
function derived from `search_accounting.per_lattice`.  The canonical JSON
metadata and agent report state: **247 perturbations across 3 parents; $(6,6)$
only; zero on $(6,4)$ and $(4,6)**.  Scientific counts, candidates, witnesses,
statuses, and the negative scoped verdict are unchanged.  Regression test
`test_coverage_statement_does_not_multiply_or_imply_empty_lattices` rejects the
old Cartesian-product wording and requires explicit zero-coverage strata.

**Lesson.**  A correct table can still be contradicted by its summary.  Coverage
sentences must be generated from the same per-stratum counters and must use
“across”, not multiplication notation, for aggregate totals.

---

## FR-016 - Theorem C3 and EXP-016 published without their translation-invariant scope
**Date** 2026-08-13 · **Track** A · **Status** CORRECTED throughout; third instance of the silent-class defect

**What we claimed.**  Theorem C3: "CP-SAT proves INFEASIBLE at $T=6$" for the
four weight-6 BB codes, presented as if it excluded *any* depth-6 one-ancilla
schedule, and "this derives the depth-7 cycle of arXiv:2308.07915".  EXP-016:
"minimum-depth schedules exhaustively enumerated", presented without naming the
enumerated class.  `checkpoints/verified_results.json` carried `thmC3` as
unqualified "exact computation".

**Why it was wrong.**  Both generating runs are orbit-restricted:
`exp004_circuit_cost_survey.py` passes `symmetry_orbits=orb` (its own field is
named `depth_is_certified_min_for_orbit_schedules` before being renamed on
output), and every `results/raw/exp004_schedules.json` entry records
`"mode": "orbit"`; `exp016`'s `env.protocol` says "translation-invariant"
explicitly.  Nothing in our artifacts excludes an unrestricted depth-6 Gross
schedule.  EXP-033 proved on three PBB members that the TI class can be more
than four layers weaker than unrestricted - exactly the gap this scoping hides.
This is the third instance of the silent-class defect (FR-005, FR-014), found
by the 2026-08-13 adversarial audit *after* the audit table had answered
"Schedule class silently restricted? No".

**What survives.**  Everything, once scoped: TI $T=6$ INFEASIBLE certificates;
depth-7 achievement (class-free); the 8-vs-7 separation (PBB $\ge8$ is
class-free via Prop C2 + Theorem C6); all EXP-016 statistics as statements
about the TI minimum-depth sample.  ASC (arXiv:2603.21499) independently
certifies the unrestricted no-depth-6 exclusion; that result is theirs, not
ours.

**Fix.**  TI qualifiers inserted at every statement of Theorem C3 and EXP-016
(README, both reports, proofs, equations, checkpoint, hypothesis generator,
verified-results generator); audit-table row corrected to "Yes, three times".

**Lesson.**  A solver's search space is part of every claim it certifies.
"INFEASIBLE" without the model's class is not a theorem statement.  Any future
scheduler artifact must carry `schedule_class` at top level, and prose must
quote it.

---

## FR-017 - Legacy EXP-012 "parent weaker in 5" carried an [exact computation] label
**Date** 2026-08-13 · **Track** B · **Status** RELABELLED legacy diagnostic; superseded by EXP-027

**What we claimed.**  "[exact computation] On $\delta>0$ codes with $n\le108$
(77 codes) the parent dominates in 49 cases, is weaker in 5, and 23 were
undecided within budget."

**Why it was wrong.**  Declaring the parent "weaker" needs a certified *lower*
bound on the PBB distance; EXP-012 used catalogue `d_reported` values
(upper-bound provenance at best) as if they were lower bounds - the
bound-direction error the strict audit protocol was later built to refuse.
Under EXP-027's strict provenance, `phase2_71` (one of the five) is
`CATALOGUE_UPPER_BOUND` with `pbb_d_lower_bound = 1` and verdict `UNDECIDED`.

**What survives.**  Two of the five - `phase2_58` and `phase2_60` - survived
certification and are the two double-verified distance-increasing reversals.
The other three are undecided, not refutations.  The corrected global tally is
EXP-027's: 66 dominations proved / 2 certified reversals / 87 undecided over
all 155 $\delta>0$ rows.

**Fix.**  Both passages relabelled "[legacy diagnostic, superseded]" with the
survivor split stated; audit-table row "Bound direction reversed?" now answers
"Yes, once - caught and corrected".

**Lesson.**  A catalogue value's provenance travels with it.  Any comparison
that flips which side of a bound a number certifies must be recomputed, not
relabelled.

---

## FR-018 - Full decoder Cartesian grid was computationally non-executable
**Date** 2026-08-13 · **Track** D (EXP-026) · **Status** FAILED FEASIBILITY GATE before production data; protocol v2 frozen

**What we pre-registered.**  EXP-026 v1 proposed 29 decoder configurations
(baseline, a $3\times3\times3$ BP/min-sum/OSD grid, and the shadow-assisted
decoder), on two codes, with 20,000 accuracy shots and 2,000 timed shots per
configuration/code pair: 1.16 million accuracy decodes before timing.

**Why it failed.**  Zero of the 58 production accuracy rows completed.  The
worker submitted a whole 20,000-shot configuration as one atomic task, so
interruption discarded every decoded shot.  Even a deliberately optimistic
baseline-equivalent estimate from the independently measured EXP-013 medians
is over 216 core-hours; higher OSD orders and iteration limits were not a
defensible lower-cost assumption, and attempted waves completed no rows.  In
the similarly configured EXP-024 PBB arm, all six workers accumulated more
than three CPU-hours without returning a 200-shot chunk; attaching a debugger
to one active worker placed it in `OsdDecoder::decode` sparse row reduction.
Shared-machine timing would also have repeated FR-004.  Thus v1 was
not a bounded experiment, and no canonical result or scientific decoder
verdict exists.

**What survives.**  The fixed circuits, noise model, detector semantics,
baseline configuration, shadow construction, two-shot API smoke, and
canonical-write guard all passed their non-production checks.  They remain
inputs to v2, not evidence for a performance claim.

**Fix.**  Before any production row existed, v2 was frozen around a finite
five-decoder set (baseline BP+OSD, shadow BP+OSD, and three BP+LSD
localisation widths), paired development and held-out faults, checkpointing
at most every 25 shots, and a supervised 30-second per-shot deadline whose
timeouts count as operational failures.  The whole latency pass is
isolation-gated.  V1 remains discoverable at
`results/partial_runs/exp026_decoder_codesign-smoke.json`.

**Lesson.**  Pre-registration does not rescue an experiment whose atomic work
unit exceeds its practical budget.  Benchmark protocols need a measured
feasibility gate, bounded per-item work, resumable checkpoints, and an
isolation gate *before* a large grid is called executable.

---

## FR-019 - Schedule-optimisation protocol spent its budget before the PBB arm
**Date** 2026-08-13 · **Track** A (EXP-025) · **Status** V1 STOPPED with asymmetric development data; symmetric v2 frozen before any PBB result

**What we pre-registered.**  EXP-025 v1 sampled 64 schedules per code plus
persisted/default extras, used 2,500 stage-1 shots per schedule, then planned
40,000 fresh shots for six finalists per code (plus the Gross default):
867,500 BP+OSD decodes in total.

**Why it failed.**  The run completed all 70 deduplicated Gross stage-1 rows
but zero of 69 PBB rows.  Gross alone consumed 19,814 wall-seconds.  The first
PBB row used the same unbounded atomic BP+OSD task that stalled in EXP-024;
stopping it returned no checkpointed shot.  Continuing would create an
asymmetric, multi-day search and could still lose a whole schedule to one
pathological decode.  The v1 data therefore cannot support a Gross/PBB
optimisation claim.

**What survives.**  The complete translation-invariant schedule catalogues,
uniform sampler, pinned schedule hashes, noiseless validation, Gross
development rows, and selection code remain useful.  The 70 Gross rows are
explicitly development/partial evidence and will not be reused in v2's
claim-bearing selection or evaluation.

**Fix.**  V2 was frozen before observing a PBB schedule result: fresh uniform
$K=16$ TI schedules per code plus the five EXP-016 schedules and Gross
default; 1,000-shot symmetric screening; three finalists per code; 10,000
fresh held-out shots; independent exact tests because equal Stim seeds across
different circuits are not paired physical faults; and the same resumable,
30-second-deadline operational decoder harness as EXP-026 v2.  Gross depth 7
and the 8,496/9,968 populations are labelled TI-class; PBB depth 8 remains
class-free exact by weight lower bound plus witness.

**Lesson.**  A symmetric comparison needs a symmetric feasibility check before
one arm consumes the budget.  Circuit seeds do not create paired trials when
the circuits have different fault-location orderings, and an optimisation
search must separate development selection from fresh held-out evaluation.

---

## FR-020 - EXP-024's frozen-v2 gate was dead code on the production CLI
**Date** 2026-08-13 · **Track** A (EXP-024) · **Status** CORRECTED before any production result

**What failed.**  `require_benchmark_v2_execution_gate()` existed, but the
normal CLI never called it.  `run()` performed the schedule work and then
entered the retired `_benchmark`/`collect` loop, which could write the
canonical artifact without using the frozen supervised, deadline-bounded,
resumable v2 protocol.  A reproduction command could therefore bypass the
protocol that the report claimed was mandatory.

**Impact.**  No benchmark conclusion survives and no production data were
generated: `results/processed/exp024_cat_extraction.json` remains a zero-byte
placeholder.  The separately named structural artifact remains valid within
its structural-only scope.

**Fix.**  The normal production entry now fails closed *before* loading the
target, solving a schedule, or building a circuit.  It first requires the
attested shared-harness gate and, even if that gate becomes complete, raises
until a real v2 driver replaces the legacy loop.  `--structural-only` and
`--validate-only` remain explicit non-production paths.  The regression test
plants spies on target loading and the legacy benchmark and proves neither is
reached; the direct CLI smoke fails at the missing gate and leaves the
canonical placeholder untouched.

**Lesson.**  Defining a gate is not enforcing it.  A claim-bearing entry point
must exercise the gate before expensive work, and a frozen protocol without an
implemented driver must hard-block rather than fall through to legacy code.

---

## FR-021 - EXP-033 noiseless evidence was not bound to its slot maps
**Date** 2026-08-13 · **Track** A (EXP-033) · **Status** CORRECTED and evidence regenerated

**What failed.**  The auxiliary Stim artifact named its canonical input only
by path and matched rows only by instance identifier.  EXP-033 uses
multi-worker CP-SAT and can persist a different valid depth-9 slot map on a
rerun.  The old test would still pass after such a rerun, so zero detector
firings and observable flips for an earlier map could be misattributed to the
currently persisted witness.

**Impact.**  The algebraic depth certificates and structural schedule checks
were unaffected, but the statement that *the persisted slot maps* had passed
the Stim check temporarily lacked provenance.  We therefore regenerated the
three noiseless checks rather than carrying the old zero counts forward.

**Fix.**  The auxiliary artifact now stores and tests the SHA-256 of the whole
canonical EXP-033 input and a canonical-JSON SHA-256 of every slot map.  Its
canonical route additionally requires the exact predeclared three-instance
set and `depth == certified_depth == w_max == 9` for every passing row.  The
current maps each produced zero detector firings and zero observable flips in
4,000 noiseless 12-round shots, and the focused regression recomputes all four
bindings from the live canonical artifact.

**Lesson.**  Matching labels and scalar summaries does not bind derived
evidence to a nondeterministically generated witness.  Persist and verify the
content identity of every claim-bearing input.

---

## FR-022 — Sector decomposition of the nontriviality disjunction is 8× slower
**Date** 2026-08-16 · **Track** C (distance) · **Status** CORRECT BUT ABANDONED ON MEASUREMENT

**What we expected.**  EXP-036 encodes "v is a *nontrivial* logical" as a
disjunction over the $2k$ logical functionals: at least one $\langle v,L_j
\rangle$ must be odd.  The hypothesis was that this disjunction is what makes
the hard $n=180/360$ rows intractable, since CDCL must reason about all
logical classes at once.  The disjunction decomposes exactly:

$$\exists v:\ \text{parity} \wedge \bigvee_j \langle v,L_j\rangle{=}1 \wedge \mathrm{wt}\le c
\iff \exists j\ \exists v:\ \text{parity} \wedge \langle v,L_j\rangle{=}1 \wedge \mathrm{wt}\le c,$$

so the full instance is UNSAT iff every sector is, and each sector is a
strictly more constrained (and embarrassingly parallel) subproblem.

**Why it failed.**  Measured on the $[[144,12,12]]$ CSS $X$-side at cap 11
(a genuine UNSAT): **monolithic 83.8 s versus sectors 670.7 s over 12
sectors** — an 8× regression.  CaDiCaL amortizes the shared parity and
cardinality structure across all logical classes within a single solve; the
decomposition throws that away and re-derives the common core 12 times.
Parallelizing the sectors does not recover it either: average sector cost
(~56 s) is not below the whole monolithic solve, so 12 cores would buy
nothing over 1.

**What survives.**  The decomposition is *correct* — statuses agreed and
bundle digests reproduce — and is retained as a tested, hash-bound
alternative (`decide_by_sectors`, `--decomposition sectors`), including the
soundness guard that exhausting a strict *subset* of sectors returns
`UNSAT_SUBSET` and can never pass a proof gate.  The default stays
monolithic.

**Lesson.**  "Split the disjunction" is not automatically a win for CDCL:
clause sharing across the disjuncts can be worth more than the extra
propagation each disjunct gains.  Measure before adopting a decomposition,
and keep the measurement in the ledger — the same afternoon, a *one-clause*
symmetry break on the identical instance bought 6.2× (77.8 s → 12.5 s) where
this restructuring cost 8×.

---

## FR-023 — Envelope check silently covered 5 of 7 reversals
**Date** 2026-08-17 — **Track** C (equivalence/envelope) — **Status** DEFECT FOUND AND FIXED

**What was wrong.** `exp036_envelope_check.py` iterated `E36.undecided_rows()` — EXP-027's
87 *undecided* rows — to collect reversals. EXP-027's own two certified reversals
(rows 331 `phase2_58`, 339 `phase2_60`) are not undecided, so they were **structurally
invisible** to the pass. The artifact's `all_reversals_css_dominated_at_equal_n: true`
therefore certified 5/7 while reading as 7/7, and its CSS candidate pool held no $n=72$
entry at all, so the $n=72$ pair could not have been checked even if collected.

**How it was caught.** Writing the $kd^2/n$ comparison table for
`reports/paper_pbb_nogo.md` §6 required a dominator for every reversal. The $n=72$ row had
no backing in the artifact — the only $n=72$ "6" in the repo was `d_z_parent` in EXP-039's
parent certificates, which is a $Z$-distance, not a certified code distance. Citing it
would have passed a $Z$-distance off as a `d_css_exact`.

**Fix.** Two reversal sources, both re-derived rather than trusted: EXP-036 tie-branch rows
(replay-stamped) and EXP-027's artifact rows, whose PBB upper bound is taken from EXP-037's
independently verified witness and required to agree with the immutable artifact. CSS
candidates now also come from EXP-037's pool (`exact` = verified witness at $U$ plus UNSAT
at $U-1$ on both sides — the same standard as EXP-036's tie branch), with each stored
witness re-verified through two GF(2) paths at load. Non-vacuity gate widened to the
*union* of both certified sets. Dominator selection canonicalised to the strongest by
$k d^2$, since `dominators[0]` made the reported comparison depend on iteration order.
Schema bumped `exp036-envelope-check-v1` -> `v2`; $7$ reversals, $162$ candidates, all
dominated.

**Lesson.** A fail-closed gate is only as wide as its input set. "All X are dominated" must
enumerate X from the *union of every source that can produce an X*, and a subset-scoped
iterator makes a universal claim quietly non-universal. Locked by
`test_envelope_check_covers_both_reversal_sources`.

---

## FR-024 — Idempotent-ideal minimum weight is not the immune parents' distance
**Date** 2026-08-21 · **Track** J (X-sector) · **Status** ABANDONED (bound true but loose)

**What we hoped.** Theorem J-G shows that for a demote-immune parent every
logical class lives in the idempotent ideal $e_IR=I^\infty$, a 2D cyclic code of
GF(2)-dimension $k_P/2$. Since $(u,0)$ with $u\in\bar e_IR$ lies in $\ker H_X$,
this gives $d_X(P)\le\min\operatorname{wt}(e_IR\setminus\{0\})$, and the ideal
has only $2^{k_P/2}\le256$ words — so we hoped for an *exact* algebraic distance
formula for the immune class, i.e. distance without a solver.

**Why it fails.** Measured on all 10 immune catalogue parents (enumerating the
whole ideal): minimum ideal weight is $18$ on the seven $(9,6)$ parents and $30$
on `phase2_109`, `15_6_0256` (15,6) and `30_6_0289` (30,6), against certified
parent distances $d_P\in\{6,8,10\}$. The bound is loose by $2$–$3\times$: the
true minimum-weight logicals are *syzygy pairs* $(u,v)$ with
$\bar au=\bar bv\neq0$, which live in $\ker H_X$ but not in $e_IR^2$, so
restricting to the idempotent block throws away exactly the vectors that realize
the distance.

**Lesson / what survives.** The inequality is valid and free, but too weak to be
useful; the live version is per-*factor* (bound the coset minimum of
$\ker H_X/S_Z$ inside each local ring $F_\chi[G_2]$ using classical cyclic-code
bounds), recorded as direction D in `notes/next_breakthroughs.md`. Numbers
reproduced by the probe in that note's D section; no artifact promoted.

---

## FR-025 — The cube-root-only trinomial criterion on odd lattices
**Date** 2026-08-22 · **Track** K (odd-lattice search) · **Status** FALSIFIED BEFORE IMPLEMENTATION

**What we almost claimed.** In characteristic two, if odd-order roots of unity
satisfy $1+u+v=0$, then $\{u,v\}$ must be the two primitive cube roots. The
proposed corollary was $3\nmid\ell m\Rightarrow k=0$ for every odd-lattice
weight-3 BB pair, which would have pruned lattices such as $(7,7)$ and $(25,7)$.

**Why it is false.** Every nonzero element of a finite field has odd
multiplicative order. In $\mathbb F_8$, if
$\alpha^3+\alpha+1=0$, then $1+\alpha+\alpha^3=0$ with
$\operatorname{ord}(\alpha)=7$. Direct rank computation gives
$\dim\operatorname{Ann}(1+x+y)=6$ on $(7,7)$ despite $3\nmid49$.
Wang–Mueller publish $[[98,6,12]]$ on $(7,7)$, and
Postema–Kokkelmans' trinomial classification includes Mersenne-prime channels
$3,7,31,\ldots$ rather than only 3.

**Correct result.** The general odd-lattice statement is the published
common-root formula $k=2|Z(a)\cap Z(b)|$. EXP-055 therefore uses field-free
annihilator ranks and scans **all** odd lattices; no divisibility prefilter is
used. On $3$-power lattices the cube-root picture remains a special case only.

**Artifacts.** `experiments/exp055_odd_lattice_sweep.py`;
`notes/theorem_k_certified_ceiling.md`; the 65-lattice census.

---

## FR-026 — Raw left-annihilator vectors used as physical logicals
**Date** 2026-08-22 · **Track** K (distance ceiling) · **Status** UNSOUND FALLBACK RETRACTED; BAR REGRESSION ADDED

**What was wrong.** `ann_basis` computes
`nullspace(poly_matrix(a).T)`: a **left**-annihilator coefficient space $I$.
An intermediate EXP-055 fallback treated $(u,0)$ with $u\in I$ as a physical
right-kernel vector. That is false unless the ideal is reciprocal-invariant.

**How it was caught.** The new witness test failed on the published
$(7,7)\,[[98,6,12]]$ code: every basis vector of raw $I$ violates
$H_X(u,0)^T=0$. Applying the involution fixes all of them:
$J=\bar I$ is the physical right-kernel pole. The bar-invariant
$[[90,8,10]]$ control had hidden the defect.

**Correct result.** On odd lattices
$$J\oplus J\cong\ker H_X/S_Z,\qquad J=\bar I,$$
the Eberhardt–Steffan principal-code isomorphism in the repository convention.
All pole ceilings and reduced witnesses now use $J$; every persisted witness is
rechecked in $\ker H_X\setminus S_Z$. The isomorphism passes 27/27 literature
audits. The earlier restricted $I\cap\bar I$ statement was sound but weak.

**Regression.** `test_bar_convention_on_noninvariant_code` and
`test_reduced_witness_is_a_real_logical` in
`tests/test_exp055_odd_lattice.py`.
