# EXP-010: Bravyi bivariate-bicycle circuit-level baseline

## Scope and result

This experiment independently reconstructed the CSS bivariate-bicycle (BB) syndrome circuit from the artifact accompanying Bravyi *et al.*, [arXiv:2308.07915](https://arxiv.org/abs/2308.07915), implemented it as a Stim Z-basis memory experiment, and decoded sampled detector data with `ldpc.BpOsdDecoder`.

All four required circuit-validation gates passed. The true artifact structure is **eight effective time slices containing seven CNOT layers**, not eight CNOT layers. For the Gross code, the implementation uses 144 data qubits and 144 ancillas, executes 864 CNOTs per syndrome round, and produced no detection events or logical flips in 10,000 noiseless shots. Stim's hypergraph logical-error search found a six-fault undetectable logical for `[[72,12,6]]`; this is a heuristic **upper bound**, not an exact circuit-distance certificate.

The complete machine-readable record is `results/raw/exp010_bb_circuit_baseline.json`.

## Primary-artifact audit and exact schedule

The four requested artifact files were read before implementation:

- `README.md:5-8` says `decoder_setup.py` constructs the check matrices, syndrome circuit, and decoding matrices. `README.md:10-13` identifies BP-OSD, and `README.md:15-21` defines the original result columns and failure criterion.
- `circuit_distance.py:1` describes that script as an upper-bound search. Its BP-OSD parameters are at `circuit_distance.py:47-52`, and the random logical searches are at `circuit_distance.py:89-146`. It does not define the syndrome schedule.
- `decoder_setup.py:45-51` gives the schedule arrays verbatim:
  - `sX = ['idle', 1, 4, 3, 5, 0, 2]`
  - `sZ = [3, 5, 0, 1, 2, 4, 'idle']`
- `decoder_setup.py:181-205` defines directions 0, 1, 2 from the three monomials of the left block and directions 3, 4, 5 from the right block, including the transposed mapping for Z checks.
- `decoder_setup.py:208-225` implements the first CNOT layer together with X-ancilla preparation; `decoder_setup.py:227-240` implements the five doubly interleaved middle layers; `decoder_setup.py:242-256` implements the last X-check CNOT layer together with Z measurement; and `decoder_setup.py:258-264` implements the final all-data-idle/X-measurement/Z-reset slice.
- `decoder_run.py:84-190` is the online noise sampler. It applies measurement flips, uniformly random single-qubit Pauli idle faults, basis-opposite reset faults, and one of 15 nonidentity two-qubit Paulis after a faulty CNOT.

`DEPTH7_SCHEDULE` stores the following pairs `(X direction, Z direction)`:

| CNOT layer | X-check direction | Z-check direction | Concurrent single-qubit action |
|---:|---:|---:|---|
| 0 | idle | 3 | prepare X ancillas |
| 1 | 1 | 5 | none |
| 2 | 4 | 0 | none |
| 3 | 3 | 1 | none |
| 4 | 5 | 2 | none |
| 5 | 0 | 4 | none |
| 6 | 2 | idle | measure Z ancillas |
| final non-CNOT slice | — | — | all data idle; measure X ancillas; reset Z ancillas |

Thus each check family participates in six CNOT layers, the middle five layers interleave X- and Z-check gates, and there are seven two-qubit layers total. The builder reconstructs the six monomial permutation matchings from `HX=[A B]`, reconstructs the transposed Z matchings, and asserts for every layer that no qubit occurs in more than one CNOT. For `[[144,12,12]]`, direct instruction counting gave 10,368 CNOTs over 12 rounds, or 864 per round (`6n`).

## Noise model and memory boundaries

`decoder_setup.py:35-40` sets initialization, idle, CNOT, and measurement error strengths all equal to `p`. The Stim implementation uses the equivalent Pauli channels:

- `DEPOLARIZE2(p)` after every CNOT;
- `DEPOLARIZE1(p)` at every scheduled data-idle location;
- `R; X_ERROR(p)` for Z reset and `RX; Z_ERROR(p)` for X reset;
- `X_ERROR(p); MZ` and `Z_ERROR(p); MX` for noisy measurement.

Data are indexed `0..n-1`, X ancillas `n..n+n/2-1`, and Z ancillas `n+n/2..2n-1`. A Z-memory experiment initializes and destructively measures data in Z. Round-zero Z-check results are detectors against the known product-state boundary; later Z detectors compare consecutive rounds. X detectors begin only at round 1 and compare consecutive X-check measurements. Terminal Z-check detectors compare the final data parity with the last Z-check measurement. The 12 observables are independent representatives of `ker(HX) / row(HZ)`, hence pure-Z logicals.

## Validation gates

Validation used `p=0.001` for detector-error-model construction and seed `230807915` for noiseless sampling.

| Gate | `[[72,12,6]]` | `[[144,12,12]]` | Outcome |
|---|---:|---:|---|
| (a) `detector_error_model(decompose_errors=False)` | 432 detectors, 12 observables, 16,164 DEM errors | 1,728 detectors, 12 observables, 67,752 DEM errors | PASS |
| (b) noiseless, 10,000 shots | 0 detection events; 0 logical flips | 0 detection events; 0 logical flips | PASS |
| (d) exact qubit count | 144 = 72 data + 36 X ancillas + 36 Z ancillas | 288 = 144 data + 72 X ancillas + 72 Z ancillas | PASS |

For gate (c), `circuit.shortest_graphlike_error()` returned `Failed to find any graphlike logical errors.` Stim's hypergraph-aware `search_for_undetectable_logical_errors` was then run with maximum intermediate detector-set size 6, maximum edge degree 6, symptom-degree-increasing moves disabled, and circuit-error canonicalization enabled. It returned an undetectable logical containing **6 physical fault mechanisms** for `[[72,12,6]]`.

**Circuit-distance conclusion:** `d_circuit <= 6`. This is a heuristic upper bound witnessed by a concrete fault set. It is not a proof that no fault set of size below six exists, so it is not an exact distance result.

## BP-OSD benchmark

Each row is an independently executed circuit simulation. A logical failure means that at least one of the 12 decoded pure-Z observables differs from Stim's sampled logical frame. The decoder parity-check matrix is the undecomposed Stim DEM fault matrix. The reported confidence interval is the two-sided, equal-tailed 95% Clopper-Pearson interval. Per-round rates and interval endpoints use

`p_round = 1 - (1 - p_shot)^(1 / rounds)`.

| Code | p | Rounds | Shots | Failures | LER/shot | 95% CP interval/shot | LER/round | Transformed 95% interval/round |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `[[72,12,6]]` | 0.001 | 6 | 20,000 | 13 | 6.500000e-4 | [3.461416e-4, 1.111263e-3] | 1.083627e-4 | [5.769858e-5, 1.852964e-4] |
| `[[72,12,6]]` | 0.002 | 6 | 20,000 | 281 | 1.405000e-2 | [1.246477e-2, 1.577870e-2] | 2.355494e-3 | [2.088334e-3, 2.647242e-3] |
| `[[72,12,6]]` | 0.003 | 6 | 20,000 | 1,493 | 7.465000e-2 | [7.104451e-2, 7.837971e-2] | 1.284730e-2 | [1.220729e-2, 1.351155e-2] |
| `[[72,12,6]]` | 0.005 | 6 | 20,000 | 9,112 | 4.556000e-1 | [4.486797e-1, 4.625332e-1] | 9.637892e-2 | [9.447455e-2, 9.830718e-2] |
| `[[144,12,12]]` | 0.001 | 12 | 2,000 | 0 | **<=1.842740e-3 (95% upper bound)** | [0, 1.842740e-3] | **<=1.536915e-4 (95% upper bound)** | [0, 1.536915e-4] |
| `[[144,12,12]]` | 0.002 | 12 | 2,000 | 3 | 1.500000e-3 | [3.094429e-4, 4.377320e-3] | 1.250860e-4 | [2.579057e-5, 3.655106e-4] |
| `[[144,12,12]]` | 0.003 | 12 | 2,000 | 57 | 2.850000e-2 | [2.165552e-2, 3.676887e-2] | 2.406600e-3 | [1.822790e-3, 3.116956e-3] |
| `[[144,12,12]]` | 0.005 | 12 | 2,000 | 1,271 | 6.355000e-1 | [6.139705e-1, 6.566315e-1] | 8.066288e-2 | [7.625584e-2, 8.522692e-2] |

The Gross-code run was reduced to 2,000 shots per point after a 20,000-shot attempt was cancelled while heavily CPU-contended and before it emitted point-level results. No result from the cancelled run is included. In particular, the zero-failure `p=0.001` row is reported only as a Clopper-Pearson upper bound, never as a zero logical-error rate.

### Decoder configuration

The executed decoder was `ldpc.BpOsdDecoder` with minimum-sum BP, variable min-sum scaling (`ms_scaling_factor=0`), parallel schedule, 200 maximum BP iterations, and OSD-0 fallback. Each of 20 forked worker processes used one decoder OpenMP thread.

This is deliberately recorded as a deviation from the original artifact. `decoder_run.py:67-72` uses minimum-sum BP, 10,000 maximum iterations, variable scaling, and combination-sweep OSD of order 7. The present OSD-0/200-iteration run satisfies the requested BP-OSD decoder contract but is not a decoder-parameter reproduction; its degradation is especially visible at higher `p`.

## Environment and reproducibility

- Host: Apple M3 Ultra, arm64, macOS 26.5.2.
- Python: 3.13.9.
- Logical CPU count reported by `os.cpu_count()`: 28.
- Per point: 20 worker processes, one OMP decoder thread each.
- `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `VECLIB_MAXIMUM_THREADS` were all set to 1.
- The host was shared during part of collection with a separate 26-worker experiment. Consequently, recorded wall times are **not isolated performance measurements**; shot counts, seeded samples, failure counts, and confidence intervals are unaffected.
- Stim seed base: `10082023008079015`. Each code was collected in a separate invocation. For p = 0.001, 0.002, 0.003, and 0.005, the 20 worker seeds were respectively `10082023008079015..10082023008079034`, `10082023008089015..10082023008089034`, `10082023008099015..10082023008099034`, and `10082023008109015..10082023008109034`. Exact per-worker seeds and shot splits are in the JSON.
- Noiseless validation seed: `230807915`.

Exact installed distributions:

| Package | Version |
|---|---:|
| stim | 1.16.0 |
| sinter | 1.16.0 |
| pymatching | 2.4.0 |
| ldpc | 2.4.1 |
| qldpc | 0.3.3 |
| numpy | 2.4.6 |
| scipy | 1.18.0 |
| networkx | 3.6.1 |
| python-sat | 1.9.dev13 |
| ortools | 9.15.6755 |

## What did and did not reproduce

### Reproduced

1. The primary artifact's direction order and asynchronous eight-slice/seven-CNOT-layer schedule were reconstructed from `decoder_setup.py:45-51,181-264` without using another scheduler.
2. Every generated two-qubit layer is disjoint. The Gross code independently gives 288 total qubits and 864 CNOTs per syndrome round.
3. The circuit-level Pauli channels match the artifact's uniform circuit-based noise primitives.
4. Both codes pass Stim DEM construction and 10,000-shot noiseless checks.
5. The six-fault `[[72,12,6]]` witness reproduces the published/artifact-style statement `d_circuit <= 6` as an upper bound.
6. For `[[72,12,6]]` at p = 0.001, the measured per-round rate is `1.083627e-4` with interval `[5.769858e-5, 1.852964e-4]`; the paper's Table 1 value `7e-5` lies inside this interval.

### Not reproduced or not tested

1. The paper's abstract reports a family threshold of 0.8%. Four finite-size points for two codes, all at p <= 0.005, do not estimate that asymptotic threshold. No threshold reproduction is claimed.
2. The paper reports about `2e-7` logical error per cycle for `[[144,12,12]]` at p = 0.001, motivating the nearly-one-million-cycle statement. The present 0/2,000 result only establishes the much weaker 95% upper bound `1.536915e-4` per round. It has insufficient statistical power to confirm or refute `2e-7`.
3. The decoder parameters differ materially: OSD-0/200 iterations here versus OSD-CS order 7/10,000 iterations in `decoder_run.py:67-72`.
4. The observable and boundary experiments differ. This run is a product-state Z-memory test with 12 pure-Z observables, noisy initial data reset, destructive final data measurement, and terminal detectors. The artifact evaluates both Pauli sectors, appends noiseless final syndrome cycles (`README.md:21`; `decoder_setup.py:513-520,568-578`), and declares failure for any nonidentity logical Pauli. Therefore the numerical tables are not a like-for-like reproduction of the published decoder experiment.
5. The original artifact uses its explicitly generated linearized single-fault matrices; this run uses Stim's undecomposed detector error model as the BP-OSD fault matrix. They implement closely related linearized decoding models but are not byte-identical constructions.

Accordingly, this work validates the schedule, qubit accounting, detector construction, noiseless behavior, and a six-fault circuit-distance upper bound. It does **not** validate the published 0.8% threshold or million-cycle lifetime claim.
