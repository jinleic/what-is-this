# pre_statement.md — `qlops/` gates A and B

Written before the first campaign run. Source paper read first-hand (arXiv HTML v2,
timestamp 22 Apr 2026 v2 line on the PDF page; DOI 10.1145/3797968, ACM TQC 7(2):14).
Authors: Kong, Zhang, Chen (Zhongguancun Lab / Tsinghua).

## Gate A — reproduce

**Quantity set** (each reproduced by pure arithmetic from the paper's own stated
inputs; full list with references in README "Transcription"):

1. Eq. (1): `p0 = 1-(1-pL)^(1/(k*d))` for the six GB codes of Table 4 -> Table 5
   `p0` column (6 numbers).
2. Eq. (2): `Q = k/((ceil(tr/tSEC)+d)*tSEC)` for the six GB rows of Table 5
   (neutral-atom QLOPS, 6 numbers) and the QLOPS-density column (12 numbers).
3. Table 5 superconducting rows: QLOPS and QLOPS density at matched p0, computed
   with Table 1 SEC times (0.86 us current / 0.40 us future, sums verified from
   components) and Table 2 reaction times (12 + 6 = 18 numbers). Requires the
   surface-code distance column as an input (not reproducible from prose; see
   "Known non-reproducible components").
4. Sec. 3.5 RSA-2048 numbers: Q_SC = 4.0314e7 (1411 logical, tSEC 1 us, tr 10 us,
   d=25), data-qubit count 1280*430 + 131*(2*25^2-1) = 714,019, ratio 56.4611;
   Q_atom = 6.8089e6 (6128 logical, tSEC 900 us), 0.7626 density; Eq. (3) ratio
   5.244; Eq. (4) ratio 2.4204; underestimation ratios ~110 / ~270.
5. Table 6 structural arithmetic: unit-qubit formula
   `2*(dX+4*dZ)*3*dX+4*dm` (13 distinct protocols, 18 rows), `Total = unit * n_units`
   integer factorization (18 rows), and inversion `cycles = 6*dm/(1-p_fail)` ->
   implied `p_fail` per row.

**Tolerance.** The paper prints no tolerance. Definition used: a number is
*PROVED-reproduced* if |computed/published - 1| <= 5% (relative). Justification:
every quantity above is arithmetic on constants the paper prints to 4-6 significant
digits; Monte-Carlo sampling error (pL measured from ~1e5+ shots) does not enter
because pL itself is taken as the published input. Expected deviations are
rounding-boundary artifacts of the `ceil` in Eq. (2) (observed in hand-checks:
<= 0.35%). A deviation > 5% therefore indicates a transcription/typo discrepancy in
the paper or in my reading, and must be recorded as a finding, per-number:

- reproduced within 5%: **PROVED-reproduced** (with exact delta in the artifact)
- not reproducible from stated inputs (missing constant / unstated convention):
  **NOT-REPRODUCED (blocked-on-missing-input)** -- recorded precisely, not guessed
- reproduced only after a correction the paper's own data implies: **PROVED-reproduced [with correction]**, correction documented.

**Gate A passes** if every reproducible number is PROVED-reproduced (the
non-reproducible set is reported separately); **fails** if any reproducible number
deviates > 5% and no documented correction resolves it.

**Known non-reproducible-from-prose components** (identified before coding):
- Fig. 1 exponential-fit constants: not printed in prose; only recoverable as
  feasibility ranges from Table 5's published distances.
- Table 5 "Distance" column (surface-code distance matching each p0): requires the
  paper's surface-code simulations; treated as input, consistency-checked.
- Table 6 "Distillation error" and "Number of cycles": produced by the paper's
  5-qubit density-matrix simulation under platform noise inputs that are not
  printed; cycles are checked structurally (6*dm/(1-p_fail) inversion) and p_out
  column is tagged [REPORTED].
- Table 6 "Total number of qubits" unit-count convention ("syndrome extraction
  cycle difference and post-selection rate"): not fully specified; checked as
  integer factorization.

## Gate B — sensitivity

**Axis 1 (decoder latency).** Multiply every reaction time `tr` (Tables 2a, 2b, 4)
by m in {0.5, 0.75, 1.0, 1.5, 2.0}; recompute every Eq. (2) QLOPS and density with
all other inputs fixed. Falsifier: any output shifts by more than 2x vs baseline
=> cross-paper comparability falsified on this axis.
Analytic note (to be confirmed numerically): for m<=2,
Q(m*tr)/Q(tr) = (ceil(m*r)+d)/(ceil(r)+d) < (2r+1+d)/(r+d) <= 2 (+epsilon), so a
> 2x shift is structurally impossible for Eq. (2) at fixed d; expected verdict:
robust, with ceil-boundary steps quantified.

**Axis 2 (magic-state protocol swap).** Baseline: Litinski 15-to-1 one-level
constants as used by the paper (Table 6: unit qubits 810..25,098, cycles 18.6..90.6,
per-protocol p_out). Replacement: zero-level CCZ constants of arXiv:2605.21867
(abstract-level, read 2026-08-29): c ~ 300 (p_L = 300 p^2), 22 physical qubits,
circuit depth 24, 3 logical qubits, space-time 22*24 = 528 qubit*cycles per CCZ.
These are tagged **[REPORTED]** until `../msd/` gate A reproduces c ~ 300 within 2x.
Test: under the paper's own rule ("infidelity of the output magic states smaller or
equal to the logical error rate" p_out <= p0), does the swapped protocol satisfy
p_out <= p0 for each of the six GB targets at the platform's physical error rates?
Resource axis: per-magic-state space-time factor 15-to-1 vs zero-level; falsifier:
output (resource or feasibility) shifts > 2x => comparability falsified on this
axis. Note the direction is expected to be opposite on the two axes (resources
shrink ~30x, accuracy limit worsens), which itself is the interesting result.

**Labels.** PROVED / CONDITIONAL / NUMERICAL / [REPRODUCED] / [DERIVED] /
[REPORTED] / [INFERENCE] as per repo contract. All runs low-priority, single core,
bounded; computation here is seconds of arithmetic.

## Revision 3 — zero-level CCZ: first-hand provenance + accounting sensitivity (2026-09-03, no sampling)
> **Revision 3's resource magnitude and “no error bars” statement are
> superseded by the Revision 4 correction below.**

Recorded after the gate-A/B campaign froze. Documentation and deterministic
re-derivation only: no new measurement, no campaign, no frozen artifact
touched. `gate_a_results.json` and `gate_b_results.json` regenerate
byte-identical to `campaigns/20260829T223540Z_dcf693ab/artifacts/` after this
revision (checked 2026-09-03), and the six behavioral tests still pass.

- **Transcription half of the `[REPORTED]` tag: DISCHARGED.** The
  arXiv:2605.21867 constants were taken from the abstract on 2026-08-29. The
  full text (v1, the only version) was read first-hand 2026-09-03 and every
  value is located in the body: c = 300 is a least-squares fit stated in
  Sec. IV and the Fig. 10 caption (the paper prints no numbered tables), over
  six p values in [1e-4, 1e-3] with 1e7–1e8 trials each, under a uniform
  circuit-level depolarizing model; no fit residuals, fit-coefficient
  uncertainty/CI, or decoder are stated. Revision 3 incorrectly said that no
  error bars existed; Figures 10 and 11 contain pointwise graphical error
  bars of unstated definition/confidence level. The 22 / 24 / 3 values appear
  in the abstract, Sec. I, Sec. III and Sec. V. Both simulated output variants
  (non-expanded and expanded to d = 7) are reported with the same 300 p².
- **Reproduction half: stays open, and the stated unblock route is closed.**
  `../msd/` cannot discharge it — its reconstruction has a distance-1 output
  frame (per-observable failure linear in p, log-log slope 0.992) and is a
  SCOPE-LIMITED-NON-TEST. The paper ships no Stim file, no repository and no
  ancillary files (Appendix A is figure-only), so reproduction needs a
  distance-carrying rebuild. Until then the physical claim remains
  `[REPORTED]`; only its provenance is now first-hand.
- **Superseded accounting record.** Revision 3 correctly recorded that 528
  is not a paper value, that 22 omits the three output patches, and that
  300p² is acceptance-conditional. It then incorrectly treated
  `22 physical qubits × 24 circuit layers` as if it shared the Litinski
  physical-qubit × syndrome-cycle basis. The then-used repo-derived legacy
  ratio also prices one zero-level attempt, not one accepted CCZ; QLOPS
  Table-6 cycles already include Litinski postselection.
- **Superseded sensitivity result.** Revision 3 incorrectly described the
  repo-derived 193×–29,947× legacy row as “published” and, because of the
  mixed time units, reported 4.1×–632× after patch/acceptance adjustments.
  Those historical numbers are retained here only as an audit trail. They
  have no common-unit magnitude interpretation and are replaced by Revision
  4 below. The qualitative observation that the resource result clears the
  preregistered 2× falsifier is rechecked on the corrected basis.
- **Falsifier unchanged** (>2× output shift ⇒ comparability falsified). No
  tolerance, axis or rule in Revisions 1–2 is modified.

## Revision 4 — common-unit accepted-output correction (2026-09-04, no sampling)

An independent read-only review found two factual defects in Revision 3.
This correction is deterministic, changes no preregistered gate or
falsifier, starts no campaign, and touches no frozen artifact.

- **The 4.1×–632× band mixed time units.** The frozen Gate-B numerator
  used 22 physical qubits × 24 *circuit layers*, while the Litinski
  denominator used physical qubits × *syndrome-extraction cycles*. The
  frozen 193.3×–29,946.5× row is therefore retained only as a
  `[DERIVED; LEGACY-MIXED-UNITS]` reproduction row, not as a published or
  physically comparable magnitude.
- **One accepted-output basis.** The zero-level paper states that depth 24
  corresponds to three surface-code syndrome rounds. QLOPS Table 6 states
  that its reported cycle counts already include Litinski postselection.
  Revision 4 therefore compares
  `7 × unit_qubits × reported_cycles` against
  `(22 + 3(2d²−1)) × 3 / acceptance` for the paper-simulated d∈{3,7}
  variants. All products, the seven-T conversion, static patch pricing,
  ratios and verdict are `[DERIVED]`.
- **Result.** The lowest sourced sensitivity arm is d=7 with acceptance
  0.30: **33.76693322683706×–5086.233185367411×**. Every enumerated
  combination remains above the unchanged 2× falsifier, so the Gate-B
  comparability verdict survives. The acceptance endpoints are approximate
  paper values and are applied across both variants only as sensitivity
  arms; the paper does not assign 0.30 or 0.40 to one named variant.
- **Uncertainty wording corrected.** Figures 10 and 11 contain graphical
  pointwise error bars, but their definition and confidence level are not
  stated. The paper reports no uncertainty/CI for fitted c=300 and no fit
  residuals; no decoder is named.
- **Evidence.** `src/evidence/zero_level_sensitivity.json`, sha256
  `1558bf3ebc493f1c3eff02bba13117c647537b78235a4c2d832ac340881e01a5`.
  Seven behavioral tests pass. Gate-A and Gate-B regenerated JSON files are
  byte-identical to the frozen campaign artifacts.

## Revision 5 — zero-level CCZ author-repository reproduction
 (control-plane registration, 2026-09-03; no sampling performed in this
 revision; no gate in Revisions 1–4 is modified)

**Goal.** Reproduce the c ≈ 300 zero-level CCZ Monte-Carlo study (Itogawa et
al., arXiv:2605.21867v1) from the authors' own repository, under the
`campaign.py` control plane, as gate **`zero-level-author-repro-r5`**. The
executable artifact of this gate is the adapter
`src/zero_level_repo_repro.py` (+ `src/test_zero_level_repo_repro.py`); the
campaign run dir is created by `campaign.py init` and is never created,
frozen, or closed by the adapter. It must be exactly under
`physics/qlops/campaigns`, its basename/timestamp/hash binder must match the
manifest minted by `campaign.py`, and `status.json` / `sha256s.txt` must remain
absent; captured manifest bytes and both markers are rechecked at every point
acceptance/write and final-artifact write. The paper's d = 7 grown claim remains
**NOT-REPRODUCED** regardless of this gate's outcome; no d = 7 circuit is
reconstructed anywhere.

**Exact source pins (primary).**
Repo `FujitsuResearch/Zero-level_CCZ_Distillation`, commit
`1b59e223590492e224bd8623a4e0bcba59029e01`, tree
`9a5d89401fd40554635fb0e9d0b4da818670c2bc`; codeload tar.gz advisory sha256
`d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`, size
1,072,583 bytes. Members are pinned by repo-relative path AND git blob
sha1: drivers `stim/ungrown/4_surface3_d3.py` =
`6205e13302cbe0285f3e704b9a66da1df3ee6173` and
`stim/grown/3_surface3_d3_expand.py` =
`592d1a8422ae2aad47b37d66e694c4eeb2b0f019`; builders
`stim/{ungrown,grown}/stim_builder.py` = `f6c25a55…` / `7d06f315…`;
surfaces `stim/{ungrown,grown}/surface_func.py` = `1dc06f18…` /
`3af84832…`; `std_calc.py` = `d2df580b…` / `0e0aaf3e…`;
`plot_LER_suc/CCZ_plot.py` = `35d58268…`; `LICENSE.txt` = `b9d9364d…`;
plus all 13 shipped circuit texts `stim/{ungrown,grown}/832_text_*.stim`
(full 40-hex table lives in the adapter and its tests, which also assert
the ticket's truncated prefixes). Tar size/sha256 mismatch, a missing or
blob-mismatched member, a duplicate member, or any traversal-shaped member
name is a refusal before anything is written. The source tar path is read
exactly once into an immutable byte capture; size and sha256 are validated on
that capture, pinned members are extracted through `io.BytesIO` over those same
bytes, and exactly those bytes are retained in the run dir. The source path is
never reopened after capture. A source tar inside the run dir is refused (no
provenance self-overwrite).

**Scope findings, preserved verbatim.**
- **F-Z1** — the PAPER reports the grown variant at d = 7 with 1e7–1e8
  trials per point.
- **F-Z2** — the REPOSITORY executes `distance_expand = 9` and embeds
  published counts taken at 1e6–1e7 shots per point. Repo d9 is never
  silently equated with paper d7; this gate reproduces the executable
  d3/d9 artifact only.
- **F-Z3** — the ungrown p = 1e-4 embedded row conflicts with its own
  comment: the comment says 5e6 shots, but the printed LER
  (3.1638657726672983e-06) times the printed acceptance
  (0.9166002 × 5e6 = 4,583,001) gives 14.5 "errors", impossible; at 1e7
  shots it is exactly 29 errors over 9,166,002 accepted. The row is
  therefore treated as 1e7-shot, and the driver's in-file
  `shots = 10 ** 7` leftover corroborates it. The authors' plotted std
  for that point (8.308705429272706e-07, CCZ_plot.py) is the 5e6-budget
  value, i.e. the published std list is itself inconsistent with the
  published LER there. The conflict is named, not silently reconciled.

**Executable schedule and budget.** Ungrown shots
[5e6, 5e6, 5e6, 5e6, 5e6, 1e7] at p = [1e-3, 8e-4, 6e-4, 4e-4, 2e-4,
1e-4]; grown [1e6, 5e6, 5e6, 5e6, 1e7, 1e7]; **71,000,000 author shots
total** (35m ungrown — five 5e6 rows plus the F-Z3-resolved 1e7 row —
plus 36m grown). One primary seeded run per point; there is no unseeded+seeded
duplication and no mid-point resume (a per-point JSON exists ⇒ it is
consumed, never resampled; identity mismatch ⇒ refusal in favor of a
fresh run dir).

**Sampling contract (mirrors the pinned drivers exactly).** The circuit
sampled is the shipped `832_text_<p>.stim` — the same artifact the
authors' drivers load. Masks are recovered ONLY by executing the pinned
driver's construction prefix (everything up to and including its single
`circuit = circuit_builder.build()` line, with the driver's own per-run
`error_rate = [...]` line parameterized — precisely how the authors
operated the repo) against the pinned builder, and reading
`Stim_builder.postselct_numbers()`; masks are never inferred from `.stim`
text. Always `detector_error_model(decompose_errors=True)`;
`matching.num_fault_ids == circuit.num_observables == 3` asserted for
every noisy arm (the noiseless p = 0 arm instead asserts a zero-error
DEM and needs no matcher). Accepted ⇔ no builder-postselected detector
fires; error ⇔ any of the 3 predicted observables differs. Rejected
shots are skipped BEFORE decode; only accepted syndromes reach the
matcher. `circuit.compile_detector_sampler(seed=seed)`; `sample()` is
never given a seed; the fixed chunk schedule (all chunks equal except
the final remainder; smoke chunk 5,000, full chunk 10,000) is part of
run identity. Full mode executes the ENTIRE smoke battery FIRST — no
certified result may bypass it.

**Smoke gates (both modes; any failure ⇒ REJECTED refusal before
scientific artifacts, never a FAIL report with exit 0).** Ungrown
noiseless: rebuilt-vs-shipped `flattened() == flattened()`; DEM carries
zero error mechanisms; 10,000 seeded shots yield acceptance 1.0, zero
actual-observable flip rows, and zero fired detectors. Grown p = 1e-3:
`num_fault_ids == 3`; 20,000 seeded shots; `decode_batch` bitwise
cross-checked against scalar `decode` on exactly the accepted rows
(mismatches ⇒ refusal; a `decode_batch ValueError` is an API semantic
failure ⇒ refusal). Both arms: same-seed fixed-schedule replay must
reproduce chunk digests bitwise.

**Full mode.** Smoke battery first; then the 12 author-count points one
at a time (ungrown then grown, p descending 1e-3 → 1e-4; fixed per-point
seeds 1101–1106 / 2201–2206; recorded schedule), each point writing an
atomic, resumable JSON. Resume acceptance requires the exact pinned recovery
record and a complete 10,000-shot chunk ledger (ordered `n`, 64-hex digest,
per-chunk count invariants, and sums matching the top-level counts), with
`reason = null` and `crosscheck = null`; a scalar-only row is refused before
new sampling. Final analysis requires exactly the 12 unique cells.
Shipped-vs-rebuilt equality is
checked at every point; inequality ⇒ refusal. Environment pins
stim == 1.16.0, pymatching == 2.4.0, numpy == 2.5.2 (mismatch ⇒ refusal:
seeded-stream identity is version-bound).

**Statistical contract.** All 12 LER and all 12 acceptance estimates are
compared to the embedded author rows with combined analytic binomial
sigma `hypot(σ_author, σ_repro)`, where σ_LER =
`sqrt(l·(1−l)/(shots·acceptance))` and σ_acc = `sqrt(a·(1−a)/shots)` —
the pinned `std_calc.py` formulas (author side at the authoritative row
budget; ungrown 1e-4 uses the 1e7 reading). Familywise 1% over all 24
comparisons: the dependence-valid **Bonferroni** bound is the
preregistered derivation — exact z = Φ⁻¹(1 − 0.01/48) =
3.529296088834735, rounded upward to **|z| ≤ 3.53** (exact Šidák
3.528022653482286 rounds up to the same gate; Šidák's independence
assumption does not hold for these 24 dependent statistics). There is NO
c ≤ 350 (or any other ad-hoc magnitude) gate. Acceptance z is always
defined; LER z is undefined only at zero accepted shots. Verdict
precedence over the defined statistics: (1) any |z| ≥ 5, or ≥ 3 values
beyond 3.53 ⇒ **FROZEN-NEGATIVE** (undefined statistics never mask a
decisive mismatch); (2) otherwise undefined statistics, or 1–2 mild
failures ⇒ **FROZEN-INCONCLUSIVE**; (3) all defined |z| ≤ 3.53 ⇒
**FROZEN-CERTIFIED**. Clean provenance with all 24 passing certifies
the **executable d3/d9 artifact only**. Any source, blob, mask, DEM,
decoder-API, circuit-inequality, or smoke-gate failure ⇒
**REJECTED/refusal before sampling**. Fits: the authors'
unweighted-through-origin c fits (least squares of LER = c·p², i.e.
Σp²·LER/Σp⁴ — their `curve_fit` on `a·x²` with no weights) are recorded
as 282.1598546872 (ungrown) and 346.5564346243 (grown) and recomputed
from the embedded rows; reproduced fits are reported with the linearly
propagated standard uncertainty sqrt(Σ[(p_i²/Σp_j⁴)² σ_i²]);
fit/slope is descriptive, never a gate. Zero-error-count clean provenance
⇒ freeze (`campaign.py freeze`)
remains the harness's act; the adapter never freezes or closes.

**Deterministic limitations (recorded in every artifact).** Seeded
detector streams are reproducible only for this stim build, this machine,
and this exact chunk schedule; NO stream identity is claimed, implied, or
testable across stim versions, machines, or schedules. Smoke repeats
re-use the same seed precisely to verify fixed-schedule bitwise count
identity within one environment.

**Tests.** `src/test_zero_level_repo_repro.py` is deterministic and
sampling-free: git blob hashing vectors and the full pin table;
author-row integer reconciliation, the 71 m schedule, and the F-Z3
conflict; the authors' published std lists reproduced exactly by the
formulas plus the corrected 1e7 sigmas (5.875141951920214e-07 LER,
8.743241581928296e-05 acceptance); through-origin fits (Σp⁴ denominator)
against both author constants; propagated σ² scaling; exact Bonferroni/
Šidák criticals; verdict boundaries and decisive-over-undefined
precedence; vectorized postselection/error counting with fake arrays and
a fake matcher (decoder sees only accepted rows; ValueError ⇒ refusal);
noiseless detector and actual-observable checks; bytecode-free pinned local
imports with global-flag restoration; chunk-schedule identity; tar traversal,
duplicate, mismatch, idempotent no-overwrite, and single-capture
extract/archive identity; campaign parent/basename/mint/stamp/live-marker
validation; strict recovery/chunk-ledger resume refusal and complete-row
resume; full-mode smoke-first wiring; CLI return codes.

## Revision 6 — builder-only zero-level smoke oracle replacement
 (2026-09-04; supersedes only the Revision-5 campaign; no result from that
 campaign is used; no gate in Revisions 1–4 is modified)

**R5 disposition and observed blocker.** The first real R5 smoke run,
`campaigns/20260904T052209Z_2a0a2723_203c9b1d1900`, did exactly what the
preregistered refusal was supposed to do: it stopped before sampling when its
ungrown p = 0 shipped-equality premise failed. The pinned current ungrown
driver is checked in with `error_rate = [1e-4]` and loads
`832_text_0.0001.stim`. R5 parameterized that current driver/builder path to
p = 0 and rebuilt a circuit with **106 detectors**, while the separate
historical `832_text_0.stim` has **58 detectors**. Thus those two artifacts
cannot be flattened-equality counterparts. Moreover, p = 0 is not one of the
12 author-sampled points: the six author points for each variant remain
p = 1e-3, 8e-4, 6e-4, 4e-4, 2e-4, and 1e-4. Treating
`832_text_0.stim` as a shipped oracle for the newly parameterized current
builder was therefore unsupported. The R5 run was frozen solely to preserve
this refusal evidence and then closed **SUPERSEDED**; it produced no sampling
result and no result from it is accepted or carried forward. Its directory
and contents remain immutable.

**Replacement identity and exact scope.** The replacement gate is
**`zero-level-author-repro-r6`**, adapter version `r6`, bound to this latest
Revision 6. A newly initialized R6 campaign is required; an R5 manifest cannot
resume under R6. Revision 6 changes only the ungrown p = 0 smoke oracle
described below. It does not reinterpret the observed R5 inequality as
success, does not special-case an equality failure, and does not weaken any
equality assertion for a circuit at an author-sampled point.

**Ungrown p = 0 is builder-only.** This smoke arm executes the pinned current
ungrown driver prefix and pinned builder with only the driver's
`error_rate = [...]` assignment parameterized to p = 0. The returned rebuilt
circuit is itself the sampled circuit, and
`Stim_builder.postselct_numbers()` from that same execution is the sole mask
source. `832_text_0.stim` is neither opened nor used as an oracle; no shipped
path or shipped blob is attributed to this recovery, and
`shipped_rebuilt_flattened_equal` is recorded as null, never fabricated as
true. Provenance explicitly records `recovery_mode = builder_only_p0`,
`sampled_circuit_source = dem_source =
rebuilt_from_pinned_driver_builder`, and `shipped_oracle_used = false`.
Builder-only recovery is executable only for the exact `(ungrown, p=0)`
smoke arm; requesting it for any author-sampled arm is itself a refusal.

The p = 0 smoke gates now require: a zero-error
`detector_error_model(decompose_errors=True)` from the rebuilt circuit;
10,000 fixed-seed shots from that rebuilt circuit with zero fired detector
bits, zero actual-observable flip rows, full acceptance, and zero errors; and
same-seed, same-chunk-schedule bitwise replay. This is a builder/DEM/sampler
API and determinism check only. It makes no shipped-circuit equality claim.

**Author-sampled equality remains mandatory.** The grown p = 0.001 smoke arm
still loads the authors' shipped `832_text_0.001.stim`, rebuilds the same point
through the pinned current driver/builder, and requires exact flattened
equality. Default/production recovery now raises an immediate REJECTED refusal
on inequality, before DEM construction or sampling. Its 20,000-shot
`num_fault_ids == 3`, accepted-row scalar-vs-batch decoder, and fixed-seed
replay gates are unchanged. Full mode still executes the complete smoke
battery first and then requires shipped-vs-rebuilt flattened equality at
**every one of all 12 author-sampled cells** before sampling that cell.

**Everything else is unchanged.** The 71,000,000-shot full schedule, all
12 author rows and seeds, 10,000-shot full chunks, per-point atomic resume
schema and complete chunk ledgers, all 24 LER/acceptance comparisons, combined
analytic binomial uncertainties, dependence-valid Bonferroni
`|z| <= 3.53` contract, decisive/inconclusive/certified verdict precedence,
fit reporting, environment/source/blob pins, single-capture tar extraction,
traversal/duplicate/no-overwrite defenses, campaign parent/mint/live-marker
checks, manifest and preregistration TOCTOU checks, and immutable
freeze/close division of responsibility remain exactly as in Revision 5.
Deterministic tests additionally pin that builder-only p = 0 works with no
shipped file present and records no equality claim, that builder-only mode
cannot bypass an author-sampled arm, that shipped inequality refuses inside
default recovery before DEM creation, and that the R6 smoke wiring and
provenance use the two distinct oracle semantics above.

## Revision 7 — smoke artifact-directory initialization
 (2026-09-04; supersedes only the Revision-6 campaign; no R6 result is used;
 no scientific gate or earlier result is modified)

The first R6 instrument run,
`campaigns/20260904T054022Z_5f040b2a_1f64afb2c983`, completed the in-memory
smoke-gate calculations but then crashed before writing `results/smoke.json`:
the smoke path invoked the exclusive atomic writer before creating its
`results/` parent directory. The resulting `FileNotFoundError` is an
instrument-orchestration bug, not a scientific verdict. That run was frozen
only to preserve its bound source/code/preregistration inputs and closed
**CRASHED**; it has no accepted smoke result and is never resumed.

The replacement gate is **`zero-level-author-repro-r7`**, adapter version
`r7`, bound to this latest Revision 7 and a newly initialized campaign.
After all smoke gates and the live-manifest/preregistration rechecks succeed,
the smoke orchestration now creates only the targeted `results/` directory
before atomically writing `smoke.json`, `summary.md`, and `inventory.json`.
A deterministic orchestration regression requires a fresh campaign-shaped
directory with no pre-existing `results/` directory and proves that all
three artifacts are created.

Nothing else changes: the Revision-6 builder-only p=0 semantics, mandatory
grown-p=0.001 and all-12-point shipped equality checks, 71,000,000-shot full
schedule, exact source/environment/control-plane pins, fixed seeds/chunks,
complete resume ledgers, 24-comparison Bonferroni gate, verdict precedence,
fit reporting, and executable-d3/d9 versus paper-d7 limitation remain frozen.

## Revision 8 — rebuilt-oracle amendment for the F-Z4 split cells
 (2026-09-08; resolves finding F-Z4 of run
 `20260908T174922Z_dc3b41ff_d04f0ae58ca0`; no earlier revision, verdict,
 or frozen artifact is modified; the REJECTED run stays frozen and is
 never reinterpreted)

**Finding F-Z4 (frozen, cited):** inside the pinned repository (commit
`1b59e22…2901`), the shipped grown circuits
`832_text_{0.0008,0.0006,0.0004,0.0001}.stim` differ from the pinned-driver
rebuild by an identical flattened delta (52 shipped-only / 32 rebuilt-only
lines; first divergence at flattened line 2086), while
`832_text_{0.001,0.0002}.stim` are flattened-equal to the rebuild.

**Mechanism diagnosis (this cycle, builder-only, from structure —
`scratch/diag-report.json` / `scratch/diag2-report.json`):** the split is
**(b) driver-version drift**, not (a) a serialization/flattening artifact
and not (c) a genuinely different construction.  Evidence: identical
968-qubit lattice with identical coordinate sets; after relabeling shipped
qubit indices onto rebuilt indices by `QUBIT_COORDS`, all 841 `DETECTOR`
definitions and all 3 `OBSERVABLE_INCLUDE`s are byte-identical; the delta is
confined to 5 of 36 tick-rounds (rounds 10–14), where the two generations
schedule and wire the reset, hookup, and measurement of three ancilla
qubits differently (net +13 CX, −3 R, −2 X_ERROR on the shipped side with
matching `DEPOLARIZE2` placement); noise strengths are correctly
parameterized per-p on both sides; the detector-error models share 14,307
of 14,466 lines (Jaccard 0.989); the interaction graphs are
WL-refinement-equal; and the noise-normalized cross-p comparison shows all
six rebuilds are one generation while the four split shipped files are
another.  Because the instruction multisets, measurement streams, and DEMs
genuinely differ, the split shipped artifacts can never serve as oracles
for the pinned driver — and because the lattices, detectors, and
observables are identical, the pinned-driver rebuild is a well-defined,
honest-broker reference for those cells.

**Amendment — rebuilt-oracle semantics for exactly four cells.**
`REBUILT_ORACLE_CELLS = {("grown","0.0008"), ("grown","0.0006"),
("grown","0.0004"), ("grown","0.0001")}`.  For these cells ONLY:

* The reference and sampled circuit is the circuit REBUILT from the pinned
  driver/builder at the cell's error rate (recovery mode
  `rebuilt_oracle_r8`); `flattened_equal=False` is expected, recorded
  factually, and is NOT a defect and NOT a refusal.
* The shipped artifact is demoted to provenance: it is opened once, its
  pinned blob sha1 and its factual flattened inequality are recorded, and
  it is never DEM'd, decoded, or sampled.
* The author row (shots, LER, acceptance), seed, 10,000-shot chunk
  schedule, and both z-comparisons for the cell are UNCHANGED; the row
  carries `oracle: "rebuilt"` and is comparable on the strength of its own
  clean recovery and sampling.

Every other cell keeps shipped-oracle semantics exactly as in Revisions
5–7 (shipped artifact sampled; flattened inequality ⇒ REJECTED refusal
before sampling).  The ungrown p=0 builder-only smoke arm is unchanged.

**Everything else is unchanged.**  The 71,000,000-shot 12-cell schedule
(run fresh, end-to-end, under this revision), all 12 author rows and
seeds, the full smoke battery executed first in full mode, per-point
atomic resume schema and complete chunk ledgers, all 24 LER/acceptance
comparisons, combined analytic binomial uncertainties, the
dependence-valid Bonferroni `|z| <= 3.53` familywise-1% contract,
decisive/inconclusive/certified precedence, fit reporting,
environment/source/blob pins, tar capture and extraction defenses,
campaign parent/mint/live-marker checks, manifest and preregistration
TOCTOU checks, and the executable-d3/d9 versus paper-d7 limitation remain
exactly as in Revisions 5–7.  Package environment pins are unchanged
(stim 1.16.0, pymatching 2.4.0, numpy 2.5.2; interpreter version is
recorded in run identity, not pinned); the Revision-8 execution host is
`mini-0.local` per the cycle host assignment (thread caps <= 6,
`nice -n 10`).

**Gate and adapter**: replacement gate `zero-level-author-repro-r8`,
adapter version `r8`, bound to this latest Revision 8 and a newly
initialized campaign.  A full-mode candidate FROZEN-CERTIFIED under this
revision certifies that the PINNED-DRIVER construction reproduces the
author schedule statistics for the executable d3/d9 artifact, with the
shipped split-generation artifacts documented as F-Z4 provenance; it does
not certify the four shipped split artifacts, and paper d=7 / physical
c≈300 remain NOT-REPRODUCED.
