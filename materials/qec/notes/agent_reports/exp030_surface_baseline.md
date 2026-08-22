POSITIVE — Both real Stim/PyMatching d=12 memory circuits were sampled at all three pre-registered noise points, closing the p=0.002 matched-k=12 Gross/PBB/surface table with the noise-convention deltas explicit.

# EXP-030 real circuit-level rotated-surface-code baseline

## Verdict and headline

At matched k=12, twelve parallel d=12 rotated patches use **3444 active physical qubits** versus **288** for Gross: **3156 more (11.958x)**.
At p=0.002 their mapped k=12 LER/round is Z-memory 2.5e-05 [8.11748e-06, 5.83407e-05] and X-memory 2e-05 [5.44935e-06, 5.12071e-05].
This is a **mapped baseline, not an identical-noise comparison**; the complete placement differences are listed below.

## Complete three-way Pareto table at p=0.002, 12 rounds, k=12

Surface resources and two-qubit gate counts are multiplied by 12 independent patches; patch depth is unchanged because the patches run in parallel. Gross/PBB statistics are imported without rerunning them.

| code / memory | data | ancilla | total | depth/round | 2q gates/round | shots / fails | 12-round block LER (95% CI) | LER/round (95% CI) | CI type |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Gross [[144,12,12]] / Z | 144 | 144 | 288 | 7* | 864 | 40000 / 356 | 0.0089 [0.00416247, 0.0136375] | 0.000744709 [0.000347536, 0.00114363] | schedule-mean t CI |
| PBB 12_6_0193 [[144,12,12]] / Z | 144 | 144 | 288 | 8* | 1008 | 40000 / 1263 | 0.031575 [0.0158458, 0.0473042] | 0.00267012 [0.00133017, 0.00403016] | schedule-mean t CI |
| 12 x rotated surface d=12 / Z | 1728 | 1716 | 3444 | 7 exact fixed-circuit TICK count* | 6336 | 200000 / 5 | 0.000299959 [9.74054e-05, 0.000699864] | 2.5e-05 [8.11748e-06, 5.83407e-05] | CP95 on one patch, monotone-mapped x12 |
| 12 x rotated surface d=12 / X | 1728 | 1716 | 3444 | 7 exact fixed-circuit TICK count* | 6336 | 200000 / 4 | 0.000239974 [6.53903e-05, 0.000614313] | 2e-05 [5.44935e-06, 5.12071e-05] | CP95 on one patch, monotone-mapped x12 |

* Gross depth 7 is exact within the translation-invariant class; ASC arXiv:2603.21499 supplies the unrestricted depth-6 exclusion. PBB depth 8 is class-free exact here (weight-8 lower bound plus a valid depth-8 witness), while the five-schedule statistical population remains translation-invariant. Surface depth is the literal fixed-circuit count (84 TICKs / 12), not a global optimum. qLDPC depth counts only two-qubit layers; surface TICK depth includes H boundaries.

## Surface sweep: every simulated point

| p | memory | seed | shots | fails | one-patch 12-round block | mapped k=12 12-round block LER (CP95) | mapped k=12 LER/round (CP95) | wall (s) |
|---:|:---:|---:|---:|---:|---|---|---|---:|
| 0.001 | Z | 20260812 | 200000 | 0 | <=1.84442e-05 (CP95 upper; [0, 1.84442e-05]) | <=0.000221308 (CP95 upper; [0, 0.000221308]) | <=1.84442e-05 (CP95 upper; [0, 1.84442e-05]) | 1.39 |
| 0.001 | X | 20260813 | 200000 | 0 | <=1.84442e-05 (CP95 upper; [0, 1.84442e-05]) | <=0.000221308 (CP95 upper; [0, 0.000221308]) | <=1.84442e-05 (CP95 upper; [0, 1.84442e-05]) | 1.39 |
| 0.002 | Z | 20260814 | 200000 | 5 | 2.5e-05 [8.11748e-06, 5.83407e-05] | 0.000299959 [9.74054e-05, 0.000699864] | 2.5e-05 [8.11748e-06, 5.83407e-05] | 2.34 |
| 0.002 | X | 20260815 | 200000 | 4 | 2e-05 [5.44935e-06, 5.12071e-05] | 0.000239974 [6.53903e-05, 0.000614313] | 2e-05 [5.44935e-06, 5.12071e-05] | 2.34 |
| 0.003 | Z | 20260816 | 200000 | 33 | 0.000165 [0.000113581, 0.000231714] | 0.0019782 [0.00136212, 0.00277702] | 0.000165 [0.000113581, 0.000231714] | 3.58 |
| 0.003 | X | 20260817 | 200000 | 40 | 0.0002 [0.000142887, 0.000272333] | 0.00239736 [0.00171329, 0.00326311] | 0.0002 [0.000142887, 0.000272333] | 3.63 |

Zero-failure points are never reported as zero: the table substitutes the 95% Clopper–Pearson upper confidence bound and labels it as an upper bound. Block-to-round conversion exactly copies EXP-010: `-expm1(log1p(-block_ler) / 12)`. The k=12 surface block first applies `1-(1-single_patch_block_ler)^12`; interval endpoints are transformed monotonically in both steps.

Gross and PBB have no p=0.001 or p=0.003 measurements in the permitted source artifacts, so those cells are intentionally unavailable rather than fabricated. Their complete imported comparison point is p=0.002; the surface sweep covers all requested p values in both bases.

## Circuit-derived surface resources

The generated circuit contains 144 final-data-measurement targets and 143 other active QUBIT_COORDS targets, hence 287 active physical qubits per patch. Stim's `num_qubits=323` is only the sparse index span (maximum index plus one), not the active physical-qubit count.
Exact fixed-circuit counts per round are 7 TICKs, 528 CX pairs, 144 data DEPOLARIZE1 target locations, 144 ancilla DEPOLARIZE1 target locations, and 144 H target locations. The equality of ancilla DEPOLARIZE1 and H target counts reflects post-H noise, not ancilla-idle noise.

## Convention deltas — mapped baseline, not identical noise

The same numeric p therefore does **not** denote the same set or number of noisy circuit locations.

1. **data idle noise.** Stim generated: before_round_data_depolarization=p applies one DEPOLARIZE1(p) to every data qubit at the start of each syndrome round, independent of whether that qubit is idle in a particular CX layer. Custom qLDPC: DEPOLARIZE1(p) is applied separately at every two-qubit schedule layer to each data qubit idle in that layer; a data qubit participating in the layer instead receives the two-qubit gate's DEPOLARIZE2(p).
2. **ancilla idle noise.** Stim generated: there is no general ancilla-idle channel; ancillas receive noise after reset, before measurement, after an H when present, and after a CX when active. Custom qLDPC: every ancilla idle in each two-qubit schedule layer receives DEPOLARIZE1(p).
3. **single-qubit Clifford noise.** Stim generated: after_clifford_depolarization=p adds DEPOLARIZE1(p) after both H layers on the X-check ancillas. Custom qLDPC: ancillas are prepared and measured directly with RX/MX and no separate post-H noise location exists.
4. **two-qubit operations.** Stim generated: a four-layer nearest-neighbour CX schedule is used and DEPOLARIZE2(p) follows every CX. Custom qLDPC: seven (Gross) or eight (PBB) controlled-Pauli layers are used; DEPOLARIZE2(p) follows every CX/CY/CZ operation.
5. **ancilla reset, basis rotation, and measurement ordering.** Stim generated: ancillas are reset with R; X-check ancillas are rotated by H; all ancillas are measured and reset together with MR.  The measurement flip is X_ERROR(p) immediately before MR, and the next reset flip is X_ERROR(p) immediately after MR. Custom qLDPC: all check ancillas are freshly prepared with RX then Z_ERROR(p), interact as controls, and are measured separately with Z_ERROR(p) then MX; the next round performs a separate RX instead of reusing an MR reset.
6. **round-boundary and terminal ordering.** Stim generated: the data-round DEPOLARIZE1 follows the leading TICK and the reset error; MR also resets ancillas after the final measured round and therefore emits a terminal after-reset error on ancillas that are never used again. Custom qLDPC: there is no once-per-round data channel or combined MR boundary; idle channels occur inside interaction layers, and no unused reset follows the last ancilla measurement.
7. **final data measurement.** Stim generated: before_measure_flip_probability=p inserts the basis-opposite Pauli just before the final M/MX data readout. Custom qLDPC: the same type of basis-opposite probability-p flip is inserted before the final data readout; this location is aligned even though the preceding round noise is not.
8. **decoder and DEM mapping.** Stim generated: the DEM is graphlike-decomposed and decoded by PyMatching. Custom qLDPC: EXP-016 used BP+OSD on an undecomposed hypergraph DEM.
9. **logical-block mapping.** Stim generated: each circuit contains one logical qubit; twelve statistically independent patches are mapped to k=12 using 1-(1-r)^12 rather than simulated jointly. Custom qLDPC: each sampled shot contains the actual twelve-logical-qubit code block and fails when any decoded logical observable is wrong.
10. **depth metric.** Stim generated: depth/round is the literal TICK count, including the two H boundaries. Custom qLDPC: the imported depth is the number of two-qubit controlled-Pauli schedule layers and excludes reset/measurement basis-operation slices.

## Imported evidence and provenance

- `results/raw/exp013_isolated_latency.json`: Gross/PBB depth and two-qubit gates per round. No latency number is used in this baseline.
- `results/raw/exp016_schedule_controlled.json`: persisted slot maps (used to count 144 data indices and 144 check/ancilla indices), p=0.002 shots/failures, schedule-mean block LER, and schedule-mean t intervals. The five schedules per code were a fixed-seed uniform sample from completely enumerated minimum-depth translation-invariant schedules.
- Neither imported artifact was rerun or modified.
- EXP-030 performed no code-distance search. Gross distance 12 has an independent exact certificate. PBB distance 12 is catalogue metadata pending an independent full-symplectic certificate. Surface `d=12` is Stim's nominal generated-code parameter; no independent physical-location circuit-distance certification was performed.

## Reproduction

```bash
cd /Users/jinleic/jinleic-workspace/qec-codesign
PYTHONPATH=src .venv/bin/python experiments/exp030_surface_baseline.py
```

Pre-registered point order is p outermost and Z then X within p; seeds are 20260812+i. Total new shots: 1200000. Experiment wall time: 14.67 s.
Environment: Python 3.13.14, Stim 1.16.0, PyMatching 2.4.0, NumPy 2.5.2, SciPy 1.18.0; macOS-26.5.2-arm64-arm-64bit-Mach-O.
Shared-machine caveat: one process ran batched PyMatching decoding while other agents could be active; load average moved from [2.4638671875, 2.451171875, 2.67138671875] to [2.28662109375, 2.41162109375, 2.65283203125]. Wall times are operational metadata, not latency benchmarks.

## Caveats

- The surface k=12 rate and resource row assumes twelve independent, identical patches operated in parallel. It is an analytic monotone mapping of one-patch samples, not a jointly sampled 12-patch circuit.
- Z-memory is the direct basis match to EXP-016's qLDPC measurement; X-memory is reported separately and is not averaged into it.
- Gross/PBB confidence intervals are schedule-mean t intervals that include measured schedule-to-schedule variation. Surface intervals are binomial Clopper–Pearson intervals conditional on Stim's fixed generated schedule. They answer related but not identical uncertainty questions.
- PyMatching and BP+OSD are different decoders, and the two DEM treatments differ. Decoder quality is part of this mapped end-to-end baseline.
- No claim of matched target logical error, identical noise-location count, identical geometry/connectivity, or globally optimized wall-clock cycle duration is made.
