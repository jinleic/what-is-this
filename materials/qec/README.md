# qec-codesign

Co-design study of quantum LDPC codes, fault-tolerant syndrome-extraction
circuits and decoders under circuit-level noise, centred on the non-CSS
**perturbed bivariate-bicycle** (PBB) construction of arXiv:2606.02418 and the
CSS bivariate-bicycle baselines of arXiv:2308.07915.

## Current research card — exact gross-code circuit distance eleven

**Question and scope (2026-09-06).** Certify or refute \(d_{\rm circ}=11\)
for the **specific three-colour, depth-eight schedule in Figure 12 / Appendix
H.2 of Strikis–Browne–Beverland, arXiv:2603.05481v1**, not an arbitrary
gross-code schedule. The underlying code is [[144,12,12]], with 72 checks of
each type and one ancilla per check. This is the sole new research mainline;
MoE drafting remains a backup, not a concurrent experiment. Union-closed
manuscript materials remain an externally unreviewed review lane; no outreach
is authorized. Earlier PBB results below are unchanged.

**Sources and admission.** [Primary paper](https://arxiv.org/pdf/2603.05481v1),
[versioned source](https://arxiv.org/src/2603.05481v1),
[authors' implementation](https://github.com/PurePhys/LR-circuits), and
[archived circuits, version 18853601](https://doi.org/10.5281/zenodo.18853601).
The archived Stim bundle contains no gross-code circuit. The complete arXiv
source contains the Figure 10–12 vector PDFs; absence of a generator script
does not preclude reconstruction. The source analysis below recovers their
printed annotations and checks the resulting interpretations mathematically.
Source bytes and figure hashes are pinned in
`../../physics/qldpc-dec/circuits/sbb2603_05481/`. An independently checked,
unambiguous transcription of the graph, colours, CNOT times and staggered
round boundaries is an admission requirement, not an assumed equivalence.
The existing IonQ/beam-search circuit and decomposed DEM are **not substitutes**.

**Fault and observable contract.** Use Section II.C's circuit-noise support:
one fault per physical operation location; idle X/Y/Z, any of the 15
nonidentity post-CNOT Paulis, preparation eigenstate flip, measurement outcome
flip. Alternatives at one location are mutually exclusive and cost one,
including a two-qubit Pauli; DEM decomposition fragments do not each become
independent physical faults. Test both X- and Z-memory experiments and all
12 logical observables. Detectors compare repeated check outcomes, with
initial known-basis checks and terminal data-parity checks. **Correction:**
Section II.C specifies a perfect terminal syndrome round, but the
Section IV.B proof explicitly includes final data-measurement errors.
“Perfect” therefore does not by itself mean noiseless destructive readout.
Direct construction confirms the **generic** convention: the pinned
`stim_circuit.py` applies data `init_error` and terminal `meas_error`;
the documented caller sets both to \(p\), and its archived surface-code
output includes the terminal error before `M`. This does **not** identify
the slightly modified Figure 12 circuit's finite head/tail, staggered round
cut or boundary idle exposure. Those target-specific details remain unresolved.
The conjecture test used one QEC cycle. The generic logical-error protocol
sets \(m=d\), hence twelve rounds for Gross; the source does not report a
separate twelve-cycle Figure 12 experiment. Neither convention fixes
the modified circuit's finite boundaries.

**Baseline, missing ingredient and novelty.** The standard Bravyi et al.
gross SEC has a reported weight-ten circuit logical. The paper proves
\(d_{\rm circ}\le11\) for non-interleaved SECs and conjectures Figure 12
saturates it after 10,000 heuristic outer loops without a sub-eleven witness.
Its Theorem 1 identifies distance with the complete residual extended-code
distance for the stated non-interleaved model; applying that theorem requires
checking its hypotheses and lifting witnesses to the actual timed circuit.
The missing ingredient is a trustworthy physical-model derivation plus a
checkable exclusion of every weight-at-most-ten nontrivial undetected fault
configuration. [Webster–Jacob–Higgott v2](https://arxiv.org/abs/2603.22532v2)
provides exact SAT/MIP/enumeration baselines (Gurobi is its recommended exact
BB-circuit baseline), not a certificate for this particular circuit.
Morphing BB circuits reported at distance twelve change the circuit class;
this is not a global hardware-performance record attempt.
No later resolution was found in the targeted search; novelty remains
**unconfirmed**, not established by search absence.

**First gate and resources.** Admit the exact source first. Then compare an
independent binary-Pauli propagation model with direct Stim fault injection;
exercise known-answer circuits, witness replay and corrupt/missing-fault,
wrong-boundary and wrong-logical controls. Only after these gates pass, run
one serial, capped bound-search assay using existing QEC machinery: at most
60 seconds per SAT call, four calls, one thread, a sampled 2 GiB RSS guard,
and 100 MiB per-output-file cap. Solver wall time is capped at 240 seconds;
source and validation stages have separate preregistered caps. Darwin did
not support the attempted hard address-space limits. Freeze effective
parameters and inputs before execution through `scripts/campaign.py`.
Success is an admitted model and reproducible feasibility evidence; an exact
distance additionally requires matching physical upper/lower evidence.
Source ambiguity or a failed independent check stops search. Timeout is
inconclusive, never a lower bound. Redirect only to resolving the concrete
source/model blocker; no automatic compute escalation or topic expansion.

### n=270 deep re-screen: 1,636/1,657 k=8 residuals decided (2026-09-08, cycle campaign)

[Campaign `20260908T194313Z_2eb90dad_d3c7ccce6e21`](campaigns/20260908T194313Z_2eb90dad_d3c7ccce6e21/)
(cycle `20260908T193320Z_ba7b32`) executed the deferred deep target-20
re-screen and the bounded Liang row b continuation. **FROZEN-INCONCLUSIVE**
per the frozen prereg mapping — the re-screen sub-objective is certified
fresh; row b advanced but is not bound.

* The live threshold-20 battery forced a fresh screen protocol: both
  lattices re-enumerated, fresh initial screens, and 2,000-try deepening of
  every \(k=8\) undecided class with threshold-20-derived seeds.
  **1,636 of 1,657 decided** (15x9: 1,457; 27x5: 179), all with physically
  re-verified witnesses at weights 12–20; shards validate strict under
  validator-v11. Remaining: 21 \(k=8\) undecided (best witness 22), 147
  \(k=12\) undecided, 210 \(k=20/24\) survivors unpromoted.
* Liang row b: base replay COMPLETE (1,326 s, exhaustive through 16,
  `validate_run_record` PASS) ⇒ the base pair is closed; prefix partition
  {55, 56, 135, 175, 269} computed; prefix55 initial run under the two-run
  stop rule. Nine prefix runs + assemble/validate remain; row b is NOT in
  `EXACT_REFERENCES`.
* One environment obstruction was repaired: mini-pro lacked the homebrew
  `libpng16.16.dylib` the pinned solver links; the identical artifact was
  provisioned (solver binary hash unchanged).
* Full findings and claims-that-do-not-hold:
  [RESULTS](campaigns/20260908T194313Z_2eb90dad_d3c7ccce6e21/RESULTS.md) ·
  [AUDIT](campaigns/20260908T194313Z_2eb90dad_d3c7ccce6e21/AUDIT.md).

### n=270 screen certified; closure open (2026-09-08, cycle campaign)

[Campaign `20260908T175247Z_e40ae2c8_8d516c12d880`](campaigns/20260908T175247Z_e40ae2c8_8d516c12d880/)
(cycle `20260908T174625Z_a9a475`) executed the documented n=270 gate:
presentation collapse then two-root screening. **FROZEN-INCONCLUSIVE** per
the frozen prereg mapping — the gate advanced, closure did not close.

* The CRT isomorphism \(\mathbb Z_{15}\times\mathbb Z_9\to\mathbb Z_{45}\times\mathbb Z_3\)
  was re-audited exhaustively (bijection, inverse, 18,225-pair homomorphism);
  both presentations enumerate 5,024 classes / 192,834 pairs, so (45,3) is a
  duplicate of the same physical family and (15,9) is the sole screened
  noncyclic presentation ([fresh audit](campaigns/20260908T175247Z_e40ae2c8_8d516c12d880/transport_audit.json)).
* Liang Table III **row a** is now bound by the two-root protocol: the
  EXP-070 certificate was re-assembled from the 12 pinned run records and
  re-validated — base exhaustion through weight 16, five prefix branches at
  weight 18 (initial+replay pairs), even-parity exclusion of 17/19,
  exhaustive pure-second-block enumeration, weight-20 witness ⇒
  **CERTIFIED_EXACT \([[270,8,20]]\)**; admitted to `EXACT_REFERENCES`
  (`exp070_odd_exact`), raising the \(k=8\) threshold at \(n=270\) to 20.
  Row b remains a source lead (≈11 pinned runs exceeded the cycle budget).
* Both screen shards passed strict validator-v11 aggregate+record
  re-validation; the EXP-063 monotone rebind (repaired this run — v11 had
  made it unusable after any battery raise; disclosed opt-in tolerant mode
  for historical snapshots only) moved 223+28 \(k=8\) classes with verified
  weight-20 witnesses to dominated. Final: 15x9 = 3,216 dominated / 187
  survivor / 1,621 undecided; 27x5 = 482 / 23 / 183.
* **Closure residuals (no overclaim):** 1,657 \(k=8\) undecided (witnesses
  only at 22–36 vs. threshold 20; deep re-screen deferred), 147 \(k=12\)
  undecided (exact fallback over budget), 210 \(k=20/24\) survivors decided
  \(d\ge5>4\) but unpromoted. Pre-staged (2026-08-24/25) screen checkpoints,
  CP-SAT payloads and dist-m4ri runs were hash-bound and re-validated, not
  re-executed; the native solver remains trusted code.
* Full findings, claims-that-do-not-hold and evidence map:
  [RESULTS](campaigns/20260908T175247Z_e40ae2c8_8d516c12d880/RESULTS.md) ·
  [AUDIT](campaigns/20260908T175247Z_e40ae2c8_8d516c12d880/AUDIT.md).

### Authorized-assay successor prepared; awaiting budget (2026-09-08)

[Successor `20260908T124248Z_78b7187d_b9950942769a`](campaigns/20260908T124248Z_78b7187d_b9950942769a/)
is registered and hash-bound to the admitted model. **Zero solver calls are
authorized; none ran.** No distance bound is claimed. The preregistration is
not a grant: solving requires a fresh explicit user authorization bound to
this run, its canonical preregistration and its validated execution seal.

Verified without a solver: all **166** frozen parent files unchanged;
**175** representative primitive replays per basis plus the wrong-boundary,
wrong-logical, corrupt-reference and same-location composition controls, now
consuming the frozen `.stim` bytes and pinned builder metadata instead of a
rebuilt circuit; **21** durable-accounting scenarios on synthetic fixtures
with Python-only probes; the real CLI refusing `status`, `solve-next`,
`reconcile` and `_child` before any model, encoding or reservation work; and
the production acceptance path verifying fixed logical controls while
rejecting **8** malformed or over-cap witnesses and a modified adjudication.

Two accounting defects were reproduced and repaired before sealing
([pre-fix evidence](campaigns/20260908T124248Z_78b7187d_b9950942769a/validation/pre_fix_guard_regressions.json)):
a measured overrun under 0.5 seconds counted as compliant, and completed
receipts were reused without rechecking retained output hashes. A solver-free
regression now covers each. One full 60-second slot is charged durably before
launch and never refunded; a detached supervisor owns the watchdog and
survives launcher death; measured time is recorded unclipped.

Sealed disposition: **READY_PENDING_EXPLICIT_USER_AUTHORIZATION**
([readiness](campaigns/20260908T124248Z_78b7187d_b9950942769a/readiness.json),
[seal](campaigns/20260908T124248Z_78b7187d_b9950942769a/execution_seal.json)).
The run stays live and holds the `math/qec` claim; it is deliberately not
frozen, because no assay conclusion exists yet.

### Source-derived reconstruction admitted (2026-09-08)

[Methodology-control successor](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/)
resolves the reconstruction. Running the authors' **unmodified** constructor
produced six noiseless generic-LRC memory controls (X and Z at 1, 2, 3
rounds) with deterministic detectors and observables. Under true
preparation-to-measurement chronology the earlier periodic criterion
rejects **zero** check pairs; sorting that same valid circuit by cyclic
clock residue falsely rejects **216**. The earlier obstruction was a
cyclic-snapshot ordering error, not a source ambiguity.

The residual law is identical under \(X\leftrightarrow Z\): hook suffixes,
with the paper's **90** reduced classes reproduced. All three upper Figure 12
panels match Figure 11 at cyclic cut 0 and all three lower panels at cut 3,
with **no endpoint permutation**. The previously explored monomial cycles
mixed axial and diagonal spokes and are not used.

From the pinned block structure the physical ticks are \(\ell-1\) for
\(L_X, R_Z, R_X\) and \(\ell+1\) for \(L_Z\), whose half belongs to the
preceding calendar period. With one unknown row origin \(\beta\), the source
round convention \(Z_r<X_r<Z_{r+1}\) forces
\(\max(D_0)-8<\beta<\min(D_0)\); the observed interval \((-2,0)\) gives
**\(\beta=-1\) uniquely**. All 864 CNOTs and ancilla boundaries are
qubit-exclusive across neighbouring rounds and all 648 overlapping pairs
satisfy the ordering. The old two-check witness becomes differences
\(2,4\) — no obstruction.

Both one-cycle memories were built in an effective and a source-style
serialized representation with identical named outcome/detector/observable
definitions; **360 endpoint moves provably cross no operation or noise** on
the moved qubit. A separate proof gives exact boundary-move equivalence and
minimum-weight invariance to other idle/boundary choices, with its
assumptions and a necessity counterexample recorded.

[Admission and metric scope](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/source_admission.json):
the target is admitted as a **source-derived canonical representative**, not
an author-identical file. Any distance value is **one-cycle** and tied to the
stated detector sets; the earlier frozen obstruction is retracted as a
blocker, with its files untouched
([correction](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/methodology_correction.json)).

**Physical-model validation passed.** Each memory basis has 1,728 physical
locations and 14,400 Pauli/record-flip options. Independent binary propagation
and Stim agree on every detector/observable column. Repaired replay checks
pass 175 representative primitive cases per basis, fixed logical and
corrupted-readout controls, the wrong-boundary control, and same-location
composition/cancellation. All 14,400 columns per basis also agree across the
two boundary serializations after semantic relabelling
([replay evidence](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/replay_validation.json)).

The frozen [comparison report](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/fault_model.json)
localizes the surface-code control's raw-record differences to 112 of 1,259
primitive-fault columns, only on records not individually certified
deterministic. There are zero differences on the 12 certified records out
of 33, and every detector/observable column agrees. These counts establish
where the raw-frame mismatch occurs, not why the representations differ.
Validation adopts this scoped agreement on the tested circuits, not general
raw-record frame equivalence.

**The distance assay is invalid; no bound is accepted.** Main incorrectly
restarted the ladder after an outer-wrapper timeout and again after a replay
bug. Per-process counters reset instead of preserving the campaign-wide
four-call/240-second budget, and report exclusivity was checked only after
solver work. Individual call answers were not durably retained. The original
failed report is preserved; elapsed time is not evidence of UNSAT.
The campaign is locked against further solver calls. Both entrypoints now
fail before dispatch; replay-only repairs used no additional solver calls.
Further bound search requires explicit new user authorization and a separately
accounted budget—not another reset of this run.
Any future UNSAT answer must be labelled solver-asserted unless its proof is
independently checked; agreement from another solver is corroboration, not a
checked proof. A published known-answer control must match the circuit and
boundaries being modelled. Failure of a capped search to find a published
weight-10 witness is inconclusive, not by itself evidence of a model error.
[Qualified result](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/result_summary.json),
[incident audit](campaigns/20260908T015946Z_8066cf36_0e6e4247e2ba/assay_incident.json).

Terminal disposition: **FROZEN-INCONCLUSIVE**. All **166 frozen evidence-file
SHA-256 checks passed**; source/model admission does not make the failed
assay a distance result.

### Superseded interpretation obstruction (2026-09-07)

The section below is **retracted as an admission blocker** by the
2026-09-08 control above; the conditional computations themselves stand and
their run is frozen.

### Coordinate-certified reconstruction — explicit interpretation obstruction

[Source-only successor](campaigns/20260907T120059Z_2a5cca9c_a7b4eb60c485/)
recovers all **36 Figure 12 annotations** with unique coordinate-based
endpoint associations. Independent PDFKit decoding agrees with the Poppler
text/vector extraction; all values match the old transcription. Direct
Figure 10 vector checks also confirm all **144 check colours** and the
four displayed long-range endpoint displacements. No closed graph, rank,
static-witness or metadata sweep was repeated.

The interpretation checks use the standard CSS extraction gate pattern:
directed check CNOTs, ancilla preparation/measurement and idles, not
undisclosed additional data gates.

The new result is a **minimal conditional obstruction**, not a distance
bound. Under the literal Figure 10 monomial association, canonical
\(X(0,0)\) and \(Z(1,1)\) share exactly \(L(0,1)\) and \(R(1,0)\).
Figure 12 gives chronological ranks \((0,4)\) versus \((4,0)\) if its
integers are read in increasing time order. In a stationary depth-eight
cycle with separate unit preparation and measurement, every reused
weight-six ancilla is saturated. The rank differences \(-4,+4\) force a
collision or an odd ancilla-to-ancilla transfer for every relative origin.
[Independent polynomial/CNOT verification](campaigns/20260907T120059Z_2a5cca9c_a7b4eb60c485/two_check_validity_certificate.json)
and a [reset-aware adjoint check](campaigns/20260907T120059Z_2a5cca9c_a7b4eb60c485/reset_observable_certificate.json)
confirm this necessary validity failure. It excludes that interpretation,
**not the authors' intended circuit or the distance-eleven conjecture**.

A common residual-compatible endpoint map is unique within the declared
block-preserving model, but neither marker-based nor caption-based family
assignment has a non-interleaved clock completion. This was checked first
with family/colour phases, then with **144 independent per-check phases**
and replayable domain-removal proofs. These are explicit restricted
interpretation results, not an exhaustive impossibility claim.

The explicit global CSS duality \(P H_X\Pi=H_Z,\ P H_Z\Pi=H_X\) is verified;
arbitrary global X/Z naming alone is not an external blocker. Conditional
on the saturated chronological model, first-CNOT time \(s\) forces
preparation at \(s-1\) and measurement at \(s+6\). The amortized depth
definition alone does not identify the target's finite cut.

**Exact target admission remains unmet.** The precise unresolved input is a
joint endpoint/order/clock interpretation that addresses the saved
four-binding witness, not access to the LaTeX or raw annotation values.
[Full derivation, scope limits and executable evidence](campaigns/20260907T120059Z_2a5cca9c_a7b4eb60c485/source_resolution.txt).
No target primitive-fault model, physical replay/control gate or distance
assay was run; GB9, the decoder contract and frozen runs remain untouched.

### Authoritative-source follow-up — external blocker (2026-09-07)

[New source-discovery campaign](campaigns/20260907T012847Z_9687c876_4f597eeab895/)
found recoverable author-owned BB144 outputs in
[commit `4eca62a`](https://github.com/PurePhys/LR-circuits/commit/4eca62a251c0ee1f56f8a5d0d52da6b8ffb31671),
dated 2025-10-13. They include X/Z partner-order JSON, colouring arrays,
residual data and two PCMs. The immediate child removed the generated outputs
on 2025-10-17. The adjacent committed drivers select surface codes, not a
BB144 finite-memory experiment; no author statement identifies these
historical files as Figure 12.

The [successor erratum](campaigns/20260907T061036Z_058f56ba_0590281d38d2/erratum.json)
corrects the frozen report's nonexistent `input/PCM/BB_144.npz` entry:
the actual paths are `circuit-generator/example-pcms/BB_144_Hx.npz` and
`circuit-generator/example-pcms/BB_144_Hz.npz`. Closed reports, statuses and
hash manifests are unchanged. The successor's
[complete source provenance](campaigns/20260907T061036Z_058f56ba_0590281d38d2/source_provenance.json)
pins three repository revisions, including the full historical constructor
and order/residual implementation, with URLs, retrieval times and hashes.

The [complete-source admission decision](campaigns/20260907T061036Z_058f56ba_0590281d38d2/source_admission.json)
distinguishes known generic clocks, repeat assembly and noisy destructive
closure from the missing **target-specific** construction. The paper calls
the target a “slightly modified LRC” but supplies no modification recipe.
Captured generic versions differ in final-loop `TICK` placement, endpoint
idle suppression and observable indexing; none is identified as Figure 12.
Combining a stock constructor with the PDF is therefore not exact admission.

The historical code **does** define its array values as data-column IDs
ordered by edge colour. A
[new-candidate correspondence check](campaigns/20260907T012847Z_9687c876_4f597eeab895/historical_candidate_check.json)
used preserved canonical results without rebuilding them: both 72-row
support sets match. Each fixed Figure 10 colour class contains **four**
normalized local orders in each basis. This count is **diagnostic-only**:
no paper-level one-order invariant, contradiction, exclusion, residual
inequivalence, Figure 12 identity or distance bound follows from it.
Figure 12's own endpoint numerals remain semantically untyped, particularly
the lower row headed **“X-check CNOT residuals.”**

Earlier discovery covered all 39 reachable main-branch commits, public refs,
issues/discussions, forks, author-linked assets and publication/deposit
relationships. The official QCTiP 2026 abstract (book page 58) and QEC 2026
poster listing (#167, board #16) confirm presentation provenance but supply
no construction details. The QCTiP extended submission requires
authentication; no public deck/poster was located. These are unavailable
**leads**, not sources known to contain the missing construction. No
authenticated private retrieval or outreach was attempted.

[Archived source-only disposition and coverage](campaigns/20260907T120059Z_2a5cca9c_a7b4eb60c485/source_admission.json):
at that stage the exact target was unadmitted. The then-current obstruction
specifies which joint endpoint/monomial and numeric timing interpretation
must be resolved before instantiating the finite circuit. An executable
listing is not mandatory if an unambiguous mathematical derivation suffices.
No closed numerical check or target physical/distance stage was rerun.
Reopen for a concrete new source, clarification or source-derived proof
addressing those constraints, not another unchanged sweep. GB9, the decoder
contract and prior frozen evidence
remain untouched.

The 07:12 UTC repository delta after 06:51:11 UTC returned no new
default-branch commits or updated issues/PRs. The exact arXiv record still
exposed v1 at 07:20 UTC. These metadata-only observations did not reopen
the closed construction or numerical checks.
The 07:38 UTC Zenodo concept/latest-version delta also resolved record
18853601: publication date 2026-03-03, modification timestamp and
publisher-reported file descriptors unchanged. No archive was downloaded
or re-inventoried.

### Canonical reconciliation — source blocker, no distance assay

[Resumed campaign](campaigns/20260906T201730Z_6c85d571_e3816a3962de/)
reuses `physics/fss-bb/src/bb_codes.py:bb_code(12,6)`, not a replacement
Gross graph. Both parity-check matrices are \(72\times144\), commute, have
rank 66 in the builder and two independent GF(2) implementations, and have
row/column weights 6/3. Thus \(k=144-66-66=12\); the twelve logical pairs
have full-rank pairing. All entries agree with the Figure 10 graph under
the explicit map, for index \(6u+v\) and coordinates modulo \((24,12)\):

| Canonical object | Figure coordinates |
| --- | --- |
| X check | \((1-2u,1-2v)\) |
| Z check | \((4-2u,4-2v)\) |
| Left data | \((1-2u,4-2v)\) |
| Right data | \((4-2u,1-2v)\) |

[Canonical checks](campaigns/20260906T201730Z_6c85d571_e3816a3962de/canonical_reconciliation.json)
and [static witness replay](campaigns/20260906T201730Z_6c85d571_e3816a3962de/static_witness_replay.json)
are distinct evidence: three existing weight-twelve static logical witnesses
replay with zero syndrome and nonzero logical pairing. The published static
distance twelve and historical cap-eleven UNSAT records are inherited;
no new lower-bound solve or portable UNSAT-proof check was performed.

The new source obstruction is **not** the old 0/16 timing result.
[Exact rational vector extraction](campaigns/20260906T201730Z_6c85d571_e3816a3962de/residual_vector_correspondence.json)
identifies all 18 Figure 11 shaded residual regions and 36 data-node
centres. **Conditionally** sorting the Figure 12 endpoint integers gives
matching residuals for all three upper/filled panels, but none of the
lower/hollow panels, even modulo the full-check stabilizer. For example,
lower red predicts pairs \(\{S,N\}\), \(\{W,E\}\); Figure 11 shows
\(\{N,NE\}\), \(\{SW,W\}\). Both triples agree modulo the check.
Independent colour cycling or reversing the order cannot remove this mismatch.

A within-block monomial-slot cycle could repair those lower patterns, but
no pinned source authorizes that non-geometric relabeling. Figure 11's
colour-cycle permission preserves an **extended-code upper bound**, not
the identity of every Figure 12 physical gate. The Figure 10 versus
Figure 11/12 X/Z marker-name conflict also remains explicit. Neither an
invented relabeling, a colour-permuted variant nor the stock generator's
different stagger is admitted as the exact target.

Remaining external input: a corrected endpoint/monomial and clock binding
for Figure 12, with its modified-LRC finite boundary and idle exposure,
or the authors' corresponding machine-readable circuit. Logical bases can
be derived locally; matching an author's arbitrary basis ordering is not
a blocker. No primitive-fault model, physical known-answer/replay gate or
distance assay was run. **No new circuit-distance bound is claimed.**
The decoder contract, GB9 and prior frozen campaign remain untouched.

### Previous source-admission campaign — FROZEN-INCONCLUSIVE

[Registered campaign](campaigns/20260906T192206Z_955565c2_c9339b06c72c/):
**STOP distance search; redirect to exact source reconstruction.** The
versioned Zenodo archive contains 124 Stim files, all from four other code
families. The authors' complete, non-truncated repository tree at
`eded83cf054ae1c79607fbd2bae98748d95f43bf` has 37 files and only a surface-code
Stim example; its BB144 matrices do not specify the proposed SEC.

Two tentative figure readings were exercised, matching graphical markers
or X/Z caption labels. The old probe assumed the endpoint annotations were
absolute CNOT clock values. **That assumption was not established from the
source.** Its six basic graph checks and numerical collision counts describe
those assumed encodings only; the 0/16 result is **not an admitted physical
timing test and must not gate the resumed reconciliation**.

| Figure reading | Simultaneous data-qubit CNOT collision slots |
| --- | --- |
| Marker-matched | 48, 72, 168, 240, 288, 240, 168, 72 |
| Caption-matched | 240, 120, 240, 168, 144, 120, 96, 168 |

The frozen numerical output is retained, but the prior timing-based stop
inference is withdrawn. Endpoint/qubit labels, local CNOT ordering and
absolute staggered clock slots must be distinguished from primary
construction evidence. A canonical graph mismatch is likewise not
established by these counts; the resumed gate checks full coordinate
equivalence, ranks, source order semantics and direct boundary construction.

[Commands and logs](campaigns/20260906T192206Z_955565c2_c9339b06c72c/execution_portable.json),
[first probe](campaigns/20260906T192206Z_955565c2_c9339b06c72c/source_admission.json),
and [caption probe](campaigns/20260906T192206Z_955565c2_c9339b06c72c/caption_interpretation/source_admission.json):
0.035366 / 0.036094 seconds probe wall time; 29,540,352 / 29,736,960 bytes
peak RSS. Native `RLIMIT_AS` and `RLIMIT_DATA` setup failed before compute;
prospective amendments used a 30-second subprocess wall limit and sampled
2 GiB RSS supervision, **not a hard address-space limit**. Failed attempts
are retained. No primitive-fault model, known-answer physical-circuit
validation, logical witness replay, SAT call or distance bound was produced.

Next admission input: an exact authors' timed circuit, or a source-justified
mapping resolving the graph/colour/time conventions, with collision-free
unit-duration preparation/CNOT/measurement/idle placement and explicit
one-cycle/12-round boundaries, detectors and 12 logical observables per
basis. Then perform the independent physical-fault and witness gates before
any solver escalation. No external request has been sent. The decoder GB9
claim, earlier PBB results and legacy `state.json` routing are untouched;
this campaign uses the CLI-returned absolute run path for freeze/close
because legacy `math/qec` inventories use `experiments/`.

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

**Odd-lattice extension (EXP-055--068).** The prior claim that $[[90,8,10]]$
was the only odd$\times$odd BB code in print was false and is withdrawn:
$29$ sourced instances now form the validation battery. The published
principal-code structure becomes, in our bar convention,
$J^2\cong\ker H_X/S_Z$ with $J=\bar{\operatorname{Ann}_{\rm left}(A,B)}$;
verified **29/29**. The pole ceiling has zero violations on the **eight**
independently exact-certified literature rows; all 27 reproduced reported
values also pass as a non-certifying sanity check. Exhaustive sweep: all $65$
odd lattices with $\ell m\le180$,
$4{,}229{,}823{,}962$ weight-$\le3$ pairs, zero $k$ mismatches. Exact screen
waves produced connected/indecomposable $[[30,8,4]]$, $[[54,8,6]]$ and
$[[126,12,10]]$ BB references. The first is a new constructor within the
checked BB tables but **globally dominated** by Grassl's $[[30,8,7]]$; the
latter two exactify Wang--Mueller BP-OSD parameters. EXP-056 additionally
exactifies Wang--Mueller's $(3,27)$ $[[162,8,14]]$: 20/20 logical-class orbits
UNSAT and replayed, plus an explicit weight-14 witness and exact $X/Z$ duality.
EXP-057--067 push the fixed point further: exact $[[170,16,10]]$, seven
$[[186,10,14]]$ classes, thirteen $[[210,18,8]]$ classes,
$[[210,24,4]]$, $[[210,14,12]]$, two $[[210,10,16]]$ classes, and two
$[[234,8,18]]$ classes. The last two are Liang et al. Table III prior art:
EXP-067 uniquely maps their twisted-torus presentations to the two open
EXP-066 bundles, then independently exactifies both by rooted connected-cluster
enumeration. A separate $n=210$, $k=8$ residual is dominated by independent
weight-16 witnesses; its exact distance is not claimed. An explicit CRT
permutation collapses the three coprime $n=210$ presentations to one search;
parity/duality one-sector ratcheting replaces the stalled CP-SAT route. None is
an end-to-end Pareto result.

Exact fixed-point closure now reaches $n=234$: **22** nonempty-frontier
lattices, **4,862** symmetry classes / **150,581** normalised pairs, with all
**4,658/4,658** referenced classes dominated, 204 high-$k$ classes without a
reference, and zero survivors or undecided. EXP-066's full
$\operatorname{Aut}(\mathbb Z_{39}\times\mathbb Z_3)$ action has 576 verified
coordinate permutations and compresses 182 initial hard classes to 30 exact
bundles; candidate-local witnesses close 28 bundles (158 classes). EXP-067
closes the final two 12-class bundles: two original/block-swapped rooted
searches, each replayed, exhaust cap 16; even kernel weight and verified
weight-18 $X/Z$ witnesses give exact $d=18$. The race-prone multithreaded
solver coordinator was rejected before certification in favour of a
hard-digest-pinned single-thread entrypoint. EXP-068 then binds physical
witnesses to all 57 legacy fallback dominations and validator v11 rebuilds
identity/transport, threshold and witness proofs on every one of the 22
shards. The $k$ census remains complete through $n=360$; no distance-closure
claim is made beyond $n=234$. At $n=270$ the screen itself is now certified
(collapse of $(15,9)/(45,3)$; $5{,}024+688$ classes validated; Liang row a
bound exact $[[270,8,20]]$ and admitted as a reference; $3{,}698$ classes
dominated) but closure is OPEN: $1{,}657$ $k=8$ + $147$ $k=12$ undecided,
$210$ $k=20/24$ survivors unpromoted, Liang row b unbound — see the
[2026-09-08 cycle note](#n270-screen-certified-closure-open-2026-09-08-cycle-campaign)
above.

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
* **[published structure reproduced + new exhaustive computation — EXP-055--068]**
  Odd-lattice BB logicals have the exact reciprocal-pole transversal
  $J^2\cong\ker H_X/S_Z$ (Eberhardt–Steffan; $J=\bar I$ is load-bearing).
  Machine audit: $29/29$ isomorphisms; pole ceiling zero violations on eight
  locally exact references; 65 lattices / $4.23\times10^9$ pairs / zero $k$
  mismatches. EXP-056 exact-certifies $[[162,8,14]]$ by 20 replayed class-orbit
  UNSATs; EXP-067 independently exactifies Liang et al.'s two
  $[[234,8,18]]$ rows by eight rooted connected-cluster runs. EXP-068 binds
  physical witnesses to all 57 fallback dominations, and validator v11 checks
  every record across all 22 shards. The Pareto screen excludes all unpromoted
  BP-OSD/Monte-Carlo estimates and uses physical witnesses against locally
  validated exact references. Artifacts: `results/processed/exp055_*.json`,
  `results/certificates/exp056_wm_162_8_14_distance.json`,
  `results/certificates/exp067_234_8_18_*_distance.json`,
  `results/partial_runs/exp068_screen_witnesses/`; note:
  `notes/theorem_k_certified_ceiling.md`; tests:
  `tests/test_exp055_odd_lattice.py`, `tests/test_exp056_odd_distance.py`,
  `tests/test_exp067_n234_connected_cluster.py`,
  `tests/test_exp068_screen_proof_repair.py`.

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
