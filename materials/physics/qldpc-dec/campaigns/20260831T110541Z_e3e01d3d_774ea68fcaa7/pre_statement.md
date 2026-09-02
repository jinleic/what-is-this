# Pre-statement: cross-paper decoder reproduction gates for the BB [[144,12,12]] code

Frozen: 2026-08-29, before any campaign run. Owner: qldpc-dec.

Basis. This statement is written AFTER a first-hand read of all four source
papers (2026-08-29) and their official code/artifact repositories. The task
description's gate A ("BP+OSD curve for p in {3e-4, 5e-4, 1e-3}") required
correction: NO paper read today publishes a BP+OSD curve at p=3e-4, and the
beam-search paper's accuracy claims are per-basis, not combined. Freezing
gates to unpublishable points would produce unfalsifiable comparisons, so the
gates are pinned to what the sources actually state. All 1e-3 points for
gates A/B/C coincide, by construction, with the beam-search paper's runtime
arms (Tables II/III), so the same campaign feeds all three gates.

## Published-side conventions (pinned, first-hand)

Circuit/DEM conventions apply to every gate:

- Circuit: `circuits/BB_144_144_12_memory_Z_p{p}_sr12_ionq.stim`, sha256
  given in README. Source: official IonQ repo
  (github.com/ionq-publications/beamsearchdecoder), directory StimCircuit/,
  files named `BB[[144,12,12]],memory_Z,error_rate={p},syndrome_rounds=12.stim`.
  Read off the circuit text: uniform DEPOLARIZE1/2(p) on every qubit after
  R/RX/M/MX/CX (measurement instruments carry p), 12 noisy SD rounds, no
  perfect final round, 936 detectors, 24 data + 72 X-type + 72 Z-type
  check-pair ancillas = 288 qubits, 12 Z-logical observables.
- DEM: `stim.Circuit.detector_error_model(decompose_errors=True,
  ignore_decomposition_failures=True)`. Direct measurement (2026-08-29):
  936 detectors, 8784 unsplit errors, decomposed component sizes {2: 864,
  3: 5328, 4: 864, 5: 864, 6: 864}. 864 Y-type errors remain as single
  3/4/6-detector graphlike hyperedges (stim did NOT split them). Gate A is
  therefore tagged [DERIVED]: an independent BP+OSD stack will not see the
  same decoder-side hyperedge convention.
- Sampling: deterministic seeds, stim `compile_detector_sampler
  (seed=...)`; seeds recorded in each campaign's seed ledger. sinter is NOT
  used as an orchestrator (it does not allow explicit seed pinning).
- LER definition (verbatim from the IonQ authors' simulation_functions.py):
  total LER = X-basis LER + Z-basis LER, where each basis LER =
  (errors/shots)/d with d=12; "errors" = any-shot logical failure of the
  12 observables in that basis's memory experiment.
- BP+OSD reference config (beam-search paper Sec. III, matching IonQ's
  bp30+osd decoder_dictionary): `ldpc.BpOsdDecoder(dem, error_rate=...,
  max_iter=30, bp_method="ms", schedule="parallel", osd_method="osd_cs",
  osd_order=10)`. Relay-BP paper's own baseline (Sec. "Flexible decoding")
  is the same decoder family at max_iter=10000, CS order 10 -- recorded as
  a second reproduction arm.

## Gate A -- reproduce the BP+OSD baseline curve for the [[144,12,12]] BB code

PUBLISHED values read from the papers directly: NONE. The beam-search
paper (arXiv:2512.07057, read 2026-08-29) reports the BP-OSD decoder only
as runtimes in Table II/III for [[144,12,12]] at p in {5e-4, 1e-3}, never
as a logical error rate number or curve read-out. The Relay-BP paper
(arXiv:2506.01779) shows BP+OSD+CS-10 as a curve in Fig. 2 but gives no
read-out table. The GARI paper (arXiv:2510.14060) explicitly notes:
"reported logical error rates for BPOSD vary substantially across papers,
possibly due to differences in circuit construction, noise model, or
decoder implementation" (verified first-hand, read 2026-08-29).

Therefore the gate quantity and tolerance must be constructed from
calibration points that ARE stated in the source papers, i.e. the decoding
time performance of the same reference pipeline:

- **Gate A quantity** (reproduce): the BP-OSD arm of beam-search paper
  Table II + Table III, basis Z, per-shot decoding wall-time, pipeline
  `ldpc>=2.4` BpOsdDecoder with max_iter=30, bp_method=min_sum, osd_order=10,
  osd_method=osd_cs, on the SAME DEMs generated from the IonQ circuits:
  - avg. time (Table II): 3.55 ms @ p=5e-4; 10.59 ms @ p=1e-3 (Apple M3)
  - p99.9 (Table III): 272.5 ms @ p=5e-4; 289.0 ms @ p=1e-3 (Apple M3)
- **Gate A accuracy companion** (no published reference; recorded as an
  internal, harness-consistency quantity, not a reproduction success
  criterion): total LER = X+Z LER at p in {5e-4, 1e-3} with 95% bootstrap
  CI, >=1e5 shots per point per basis, plus an extra point at p=3e-4.
- **Falsification criterion**: avg/p99.9 wall time agree within a factor
  2 of the Table II/III values on this M3-class machine (hardware gate,
  CONDITIONAL on same-process measurement protocol of the IonQ authors:
  10k shots of the identical compiled BpOsdDecoder object, one perf_counter
  call each, no batch/sinter dispatch). A factor-2 agreement ON THIS
  HARDWARE is the reproduction bar; hardware-dependent runtime identities
  are cross-checkable by throughput ratios (p99.9/avg) rather than
  absolute ms.
- **Self-consistency checks (required, not falsification criteria)**:
  (i) X and Z basis LERs agree within 3 sigma on >=1e5 shots (paper does
  not publish per-basis asymmetry beyond d/2 normalization);
  (ii) shots per point recorded with 95% bootstrap CI on the LER.

## Gate B -- beam search decoder reproduction (arXiv:2512.07057)

Source of truth: Fig. 2 (p-sweep of per-basis LER), Table I (decoder
params), Table II (avg. time), Table III (p99.9 time).

Pinned published quantities, all at p=1e-3, [[144,12,12]], 12 rounds (all
read directly from paper text/tables on 2026-08-29):

| decoder | per-basis LER ratio vs bp30+osd (Fig. 2) | avg time (ms, Table II) | p99.9 (ms, Table III) |
|---|---|---|---|
| bp30+osd   | 1.0 (baseline) | 10.59 | 289.0 |
| beam8_230iters  | ~=1 (same LER) | 2.318 | 11.01 |
| beam32_340iters | 5.6x lower  | 3.837 | 14.18 |
| beam64_640iters | 17x lower  | 6.698 | 17.70 |

- **Gate B quantity** (reproduce): same quantity as published, namely
  per-basis logical error rate ratios AND per-shot runtime for our own
  from-scratch Python beam-search decoder (Algorithm 3 of the paper,
  re-implemented in numpy from Appendix A/B equations), on the SAME DEMs
  (IonQ circuits, decomp DEM), with the published parameter table for
  beam8_230iters (beam_width=8, initial_iters=30, iters_per_round=20,
  max_rounds=10, num_results=1), beam32_340iters (32, 40, 30, 10, 1),
  beam64_640iters (64, 40, 30, 20, 1) from Table I exactly.
- **Ler ratio reference**: at p=1e-3, beam8_230iters's total LER must be
  within the 95% bootstrap CI of the harness-baseline bp30+osd LER
  (target: ratio in [0.87, 1.15] around 1.0); beam32 (5.6x) and beam64
  (17x) determinations are quantitative ratios: OUR beam32/bp30+osd ratio
  must match the published 5.6x ratio within [4.5, 7.0] to reproduce; the
  beam64 ratio must match within [11, 22] (published 17x). These
  tolerance bands are OUR choice (the paper does not publish CIs on the
  ratios, only the point estimates); they are sized as +-25-30% around
  the published figure, an honest reproduction band for a from-scratch
  reimplementation, and each batch ratio carries its own bootstrap
  uncertainty on top.
- **Falsification criterion**: beam8 LER ratio outside [0.87, 1.15], or
  beam32 ratio outside [4.5, 7.0], or beam64 ratio outside [11, 22] at
  p=1e-3 with >=1e5-shots-per-arm evidence, = the corresponding paper
  claim NOT reproduced under this harness. p99.9 speedups (26.2x for
  beam8, 20.4x for beam32) are SECONDARY reproduction targets under the
  same factor-2 bar as gate A (hardware-coupled, CONDITIONAL).
- **Config set**: at minimum arm beam8_230iters works E2E (first
  falsification target). beam32/beam64 run only after beam8 passes.

## Gate C -- GARI-NMS ensemble ensemble claim (arXiv:2510.14060)

Source of truth: Fig. 4(c) (D_Z curve, [[144,12,12]], uniform
depolarizing), Sec. III A "Logical error rates", and the abstract
headline. All read first-hand 2026-08-29.

- Published quantity: GARI-NMS-x24 on gross code, uniform depolarizing
  circuit-level noise, 12 rounds + final Z-basis data measurement,
  per-round LER = 6.70e-9 at p=1e-3. Uncertainty bound: "99% confidence
  ... at least 100 decoding failures" => relative statistical error of
  1/sqrt(100) = 10%, now frozen at +-1.93e-9 (the abstract's number).
- **Crucial comparisons constraints read from paper**: error bars 99% CI
  from >=100 failures; scipy data on Zenodo (DOI 10.5281/zenodo.17991975);
  code at github.com/astra-decoders/gari-nms. NOTE: the paper does NOT
  pin the shot count for the 1e-3 gross-code point; CIs from the
  >=100-failure rule implies at least ~1.5e10 shots (because 100 /
  6.7e-9). This is not feasible in this harness at the pinned shot
  budgets. [INFERENCE, labeled as such]
- **Gate C quantity** (reproduce as LOWER BOUND, not point): OUR ensemble
  NMS decoder (24 randomized serial-schedule decoders on the plain DEM,
  published NMS normalization 0.96875, max 400 iters, single Z-D matrix
  on the Z-basis memory) must reach per-round LER <= 1.93e-9 at p=1e-3
  with >=100 decoding failures within the campaign budget before the
  harness is credited with reproducing the GARI claim. Below that bar,
  the harness CLAIM is restricted to: "ensemble NMS reaches LER below
  the 1.93e-9 floor" (sufficient condition, not the published point).
  FALSIFICATION: our ensemble fails to beat 1.93e-9 per-round in >=1e9
  effective shots per basis.
- **Explicit non-goals** (honesty about the unavoidable gap): we do NOT
  reproduce (i) the GARI graph augmentation/rewiring (4-cycle-free
  hyperedge-rewritten decoder graph -- a different decoding problem
  representation), (ii) the >=1.5e10-scale shots needed for a pointed
  6.70e-9 determination to the paper's CI. Gate C as runnable HERE is
  the ensemble-NMS WITHOUT GARI on the plain decomp DEM. The
  falsification semantics are correspondingly weaker: failing to reach
  the floor does NOT refute the paper (their graph transform is the
  claimed innovation); succeeding DOES confirm that a 24-ensemble of
  random-schedule NMS on the plain DEM is already sufficient. Either
  way this is a meaningful, publishable-consistent result.

## Shot/cost budget summary

| gate | point(s) | shots/basis (min) | hardware cost estimate |
|---|---|---|---|
| A runtime (Table II/III) | p=1e-3 Z | 10,000 timed decode() calls | ~2 min at 10.59ms avg |
| A accuracy | 3e-4,5e-4,1e-3 X+Z | 1e5 per (basis, p) | ~E5 x 6 decode calls at avg ~5-10ms = ~1.5-3 h |
| B beam8 ratio | p=1e-3 X+Z | 1e5 per arm, 2 decoders | ~E5 x 2 x (1-10ms) = ~0.5-2 h |
| B beam32/64 ratio | p=1e-3 X+Z | 1e5 per arm | +0.5-2h per arm (slower decoders) |
| C ensemble-NMS | p=1e-3 X+Z | 1e9 scale needed for 1.93e-9 floor | NOT FEASIBLE in single-campaign; multi-session continuation, bounded per session |

Gate C at the published point is explicitly declared out of budget reach
for this harness; the deliverable is the lower-bound form plus the honest
statement of what is out of reach. If 1e9+ shots become available via
continuation across sessions in bounded chunks, the floor criterion
above applies.

## Sources read first-hand (2026-08-29)

- arXiv:2512.07057 (beam search; v2 list, IonQ public repo simulation
  code + stim circuits + notebook read first-hand)
- arXiv:2506.01779 (Relay-BP; v2 list, IBM relay repo README + code +
  testdata circuits read first-hand)
- arXiv:2510.14060 (GARI; v2 list, code/Zenodo pointers read first-hand)
- arXiv:2603.19062 (fair baselines; v3 list, read first-hand)

## Signature

These gates and tolerances are frozen as of 2026-08-29, before any
campaign run in this harness. Changes to a gate require a documented
revision entry in this file (append below, do not silently edit).

## Revision GB1 — pinned million-shot beam8 campaign (2026-08-31, pre-data)

No published target, tolerance, circuit, DEM, decoder parameter, sampling
stream, or LER-ratio definition changes.

- Increase the first beam8 determination from the minimum 1e5 to **1e6
  shots per arm**, basis Z, p=1e-3. Each arm retains its original independent
  `sampling_seed(20260829, p, basis, arm)` stream. BP+OSD shards one already
  sampled stream; sharding changes execution only.
- Execute beam8 with `src/beam_cpp/beam8_cpp`, single-threaded. This is a
  backend substitution only: formal per-shot comparisons against the frozen
  Python decoder found **0/1000 mismatches** at both p=1e-3 and p=3e-3
  (seed 20260830; reports frozen under `src/beam_cpp/scratch/`).
- Freeze the decision rule before decoding: independent-binomial percentile
  bootstrap, 100000 draws with a recorded derived seed. Ratio CI wholly
  inside [0.87,1.15] => REPRODUCED; wholly outside => NOT-REPRODUCED;
  overlap/crossing => INCONCLUSIVE. This is deliberately stricter than
  promoting a point estimate alone.
- Record both raw shot-failure probability and per-round LER (raw/12). Their
  beam8/BP+OSD ratio is identical, so this clarification cannot move Gate B.

## Revision GB2 — sparse-ratio bootstrap tail correction (2026-08-31, pre-data)

Recorded while the million-shot decode was still running, before any shard
result or arm failure count existed. No circuit, sample, decoder, shot count,
seed, point estimate, target band, or three-way decision rule changes.

- Bug: the frozen helper dropped every non-finite bootstrap ratio. That
  incorrectly removes positive-numerator/zero-denominator draws and truncates
  the upper tail when observed failures are sparse. Reproduction:
  5/1 failures in 1e6 shots with the frozen 100000-draw design returned the
  finite CI [0.667, 9.0], even though denominator-zero draws occur often
  enough that the 97.5th percentile is unbounded.
- Correction: retain positive/zero draws as `+inf`; discard only zero/zero
  draws, whose ratio is undefined. Encode an unbounded JSON endpoint as
  `null`. Keep the same independent-binomial draws, draw count, derived seed,
  and percentile rule.
- The already-running campaign imported the old helper at process start.
  Its samples, decoder outputs, failure masks, raw counts, and rates remain
  authoritative, but its ratio CI/verdict will be marked superseded. A
  separate immutable amendment will recompute only the CI/verdict from those
  raw counts under this pre-data correction; the original campaign directory
  will not be edited.

## Revision GB3 — BP+OSD heterogeneous-prior correction (2026-08-31, post-interim)

Recorded after the beam8 output and 14 of 20 BP+OSD shards from the
million-shot campaign existed. This is therefore a disclosed post-interim
protocol correction, not a pre-data amendment. No circuit, sample stream,
decoder iteration/OSD parameter, shot count, target band, bootstrap seed, or
decision rule changes.

- Source audit: the paper pins `stimbposd==0.1.0`. In its tagged source
  (commit `7921f5eb1b358ff616f9822280c9961e83df06cb`),
  `stimbposd.BPOSD` constructs `ldpc.BpOsdDecoder` with
  `error_channel=list(priors)`. The local adapter instead averaged the merged
  DEM priors and passed that scalar as `error_rate`.
- The distinction is material and structural: the pinned p=1e-3 Z DEM has
  8784 merged columns and 9 distinct prior probabilities, ranging from
  0.000533333333333148 to 0.003721584619834077. A scalar mean is not the
  published baseline. The adapter now passes the full vector exactly; an
  explicit scalar override remains experimental only.
- The original million-shot campaign is retained as invalidated evidence; its
  scalar-prior BP+OSD counts cannot decide Gate B. The corrected replacement
  independently regenerates the same two seeded sample streams and reruns both
  arms under the corrected sources. Acceptance requires both regenerated
  sample-file hashes to equal the original pre-decode hashes.
- Revision GB2 still governs sparse ratio tails. The replacement campaign
  snapshots the corrected bootstrap helper and emits the final Gate B verdict
  directly, so no GB2-only amendment is needed for the valid replacement.

## Revision GA1 — apply the heterogeneous-prior correction to Gate A (2026-08-31)

Recorded after the GB3 source correction, while the corrected descriptive
accuracy companion was running and before any corrected Gate-A timing sample
existed. At recording time, only its first row was available: X memory at
p=3e-4 had 0 failures in 100000 shots. Gate-A accuracy has no published
acceptance target, so that interim count cannot move a reproduction verdict.

- Retract the historical scalar-prior Gate-A accuracy and timing rows. Reuse
  the frozen circuits, DEMs, shot counts, seeds, decoder parameters, per-round
  normalization, published timing references, and factor-2 timing bar.
- The corrected accuracy companion uses the full merged-DEM heterogeneous
  channel and snapshots every relevant local source. It remains descriptive.
- Run the corrected 10000-shot single-decode timing arm only after other
  decoder campaigns leave the host, so campaign contention does not bias the
  hardware-coupled measurement. The timing artifact must record both a source
  snapshot and an exact configured-channel audit before it can decide Gate A.

## Revision GA2 — basis-correct the derived p=3e-4 companion (2026-08-31)

Recorded after both corrected 100000-shot p=3e-4 accuracy rows had completed,
and before the remaining four companion rows or any corrected timing sample
were available. Both observed rows had zero failures; accuracy remains
descriptive and has no reproduction acceptance target.

- Bug: for p=3e-4, `load_circuit(p, basis)` appended the derived **Z** circuit
  regardless of `basis`. The row labeled X therefore sampled and decoded a Z
  memory circuit; it is invalid. The Z row is basis-correct and remains valid.
- Correction: derive an X-memory p=3e-4 artifact by applying the same exact
  textual noise-argument rescale `(0.001)` to `(0.0003)` to the committed
  IonQ X-memory p=1e-3 circuit. The circuit sha256 is
  `300f8625ac81c9bf35f16b85b3bebf6aa087947650209387ab0f688168f372a2`;
  its 936-by-8784, 12-observable DEM sha256 is
  `d1d6fdb838d407155ec1a5940b8e1d32b7f2c32d5661dea84dad3feb8f4d4579`.
- Freeze a focused 100000-shot replacement of only the invalid X row, using
  the original X sampling seed and corrected vector prior. Do not mutate the
  in-flight six-row campaign. Any combined Gate-A accuracy table must take
  p=3e-4 X from the focused replacement and the other five rows from the
  original corrected-prior campaign.

## Revision GB4 — recover the NumPy-bool close failure (2026-08-31)

Recorded after all scalar-prior BP shards and all corrected-prior BP shards
had completed, while the corrected beam subprocess was still running. This is
an artifact-serialization correction only; it cannot alter samples, decoder
outputs, failure counts, rates, bootstrap draws, target band, or verdict rule.

- Bug: the ratio helper could return NumPy scalar booleans. Both arm results
  were appended successfully, but `json.dump(summary)` rejected
  `ci_inside_band` during campaign close, leaving a truncated summary and no
  immutable inventory.
- Correction: force percentile endpoints and verdict flags to built-in Python
  scalar types. For a process that imported the old code before launch, use
  `gateb_pinned.py --recover <campaign>` only after its parent exits. Recovery
  must validate both sample hashes, all shard ranges and aggregate fields,
  beam predictions and failure mask, snapshot the recovery sources, then
  write the summary and inventory. It refuses an already-closed or live
  campaign.
- The scalar-prior campaign is closed with outcome
  `INVALIDATED_GB3_scalar_prior`; no ratio is evaluated. A corrected-prior
  campaign uses the unchanged GB2/GB3 ratio computation.
