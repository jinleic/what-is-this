# `qldpc-dec/` — cross-paper decoder reproduction harness for BB codes

**Status: ACTIVE — Gates A/C closed at their stated scopes; pinned Gate B
million-shot beam8/BP+OSD determination running; prior gross Gate-D
certificate evidence retracted after an AIS-weight audit, with corrected
small-code exact validation closed (2026-08-31).**
One deterministic harness for bivariate-bicycle (BB) codes under circuit-level
noise, used to reproduce or refute published decoder-performance claims.
Receiver-side only — no code construction; algebraic qLDPC stays in `../math/qec/`.

## Contract (fixed before any campaign)

- Code: [[144,12,12]] gross code. **Primary circuits: published stim circuits
  from the official IonQ beam-search repo**
  (github.com/ionq-publications/beamsearchdecoder, `StimCircuit/`), committed
  under `circuits/` with sha256 (see `dem/PROVENANCE.json`). These are
  12-round memory experiments, uniform DEPOLARIZE1/2(p) (m.m. instruments at p),
  936 detectors, 12 logical observables, 288 qubits. The in-repo
  `src/qldpc_dec/code.py` generator (Bravyi et al. 2024 polynomial pair
  A=x³+y+y², B=y³+x+x², Nature 627, 778, arXiv:2308.07915) is a documented
  fallback, not the campaign circuit.
- Noise: one fixed circuit-level depolarizing DEM per (p, basis), derived as
  `stim.Circuit.detector_error_model(decompose_errors=True,
  ignore_decomposition_failures=True)` — same convention as the IBM relay
  repo loader (arXiv:2506.01779). DEMs frozen under `dem/` BEFORE any campaign;
  936 dets × 8784 error mechanisms each; 864 Y-type mechanisms stay as
  graphlike hyperedges (documented; makes the reproduction [DERIVED] per the
  parent README's terminology).
- Decoders (all three smoke-validated end-to-end):
  - BP+OSD baseline: `ldpc` 2.4.1 `BpOsdDecoder`, max_iter=30, ms, osd_cs
    order 10 (matches the IonQ paper's bp30+osd arm, arXiv:2512.07057 Sec. III).
  - Beam search (gate B): from-scratch numpy reimplementation of Algorithm 3 +
    Appendix A/B of arXiv:2512.07057 in `src/qldpc_dec/beam_search.py`.
  - Ensemble NMS (gate C): 24 randomized-serial-schedule NMS decoders,
    alpha=0.96875, max 400 iters, minimum-latency ensemble stopping
    (arXiv:2510.14060 configuration) in `src/qldpc_dec/nms.py`.
- Shots: ≥1e5 per point per basis (gate A); seeds derived from base seed
  20260829 via `src/qldpc_dec/seeds.py` (stable SHA-256 derivation), recorded
  per run in the campaign inventory; 95% bootstrap CIs via
  `src/qldpc_dec/bootstrap.py` (Wilson + binomial-percentile).
- Circuit-level rates follow the official IonQ script: first compute the raw
  any-logical shot-failure probability, then divide by 12 syndrome rounds for
  the reported per-round LER. `run_point()` records both quantities. A
  2026-08-31 audit found and corrected an earlier missing `/12`; frozen
  artifacts remain raw and are relabeled rather than rewritten.

## Gates

Pinned in `pre_statement.md` (2026-08-29, frozen before any campaign run).
Summary:

1. **Gate A (reproduce):** the *published* BP+OSD reference quantities for
   [[144,12,12]] are the decoding-time numbers of arXiv:2512.07057 Table II/III
   (mean 3.55/10.59 ms, p99.9 272.5/289.0 ms at p ∈ {5e-4, 1e-3}, single-core
   Apple M3). Falsification bar: within factor 2 on this machine under the
   same single-shot perf_counter protocol. Published BP+OSD *LER curves* at
   the task's p ∈ {3e-4, 5e-4, 1e-3} do NOT exist in any of the four sources
   (documented in pre_statement.md); LERs are recorded as an
   accuracy-companion quantity, not a reproduction criterion.
2. **Gate B (refute-or-confirm):** per-basis LER ratios of beam8/beam32/beam64
   vs bp30+osd at p=1e-3 (arXiv:2512.07057 Fig. 2): beam8 ≈ equal accuracy
   (band [0.87, 1.15]), beam32 = 5.6× (band [4.5, 7.0]), beam64 = 17× (band
   [11, 22]) lower LER. ≥2σ / band miss = refuted published claim.
   Absolute-ms runtime claims are hardware-conditioned (our beam search is
   numpy; only ratios are gate quantities).
3. **Gate C (refute-or-confirm):** GARI ensemble claim, LER
   (6.70±1.93)e-9 per round at p=1e-3 (arXiv:2510.14060 abstract; 99% CI from
   ≥100 failures ⇒ ≥1.5e10 shots for the point determination — outside this
   harness budget; pre-statement froze the reachable form: our 24-ensemble
   NMS on the plain DEM must reach ≤1.93e-9 per-round to confirm the weaker
   claim; the GARI graph transform itself is an explicit non-goal).

## Pre-campaign requirements (discipline) — DONE 2026-08-29

- [x] Owner first-hand reads: arXiv:2512.07057 (v2, 2025-12-17 version),
      2506.01779 (v2, 2025-08-22), 2510.14060 (v2, 2026-03-06), 2603.19062
      (v3, 2026-04-26) — full texts + official repos read.
- [x] Per paper: exact DEM construction, shots, runtime percentile
      definition, figure/table pinned (see `pre_statement.md`).
- [x] One-page pre-statement per gate with falsification criteria,
      committed BEFORE first campaign run.

## Smoke test (2026-08-29, both committed under campaigns-smoke/)

- `20260829T225407Z_53700f78_87fb945aabf3` — 1e3-shot end-to-end:
  BP+OSD raw shot-failure rate 0/1000 (raw CI95 ⊂ [0, 3.8e-3]), beam8 0/1000,
  NMS-ensemble 0/25 (2-decoder rehearsal), timing arm x50. Pipeline proof
  only; no gate verdicts.
- `20260829T230642Z_755286ce_1f387aff4429` — 5-shot validation that every gate
  point (p ∈ {3,5e-4,1e-3} Z + 1e-3 X, both decoder arms) decodes cleanly,
  including the derived p=3e-4 circuit.

Measured rehearsal costs (this machine, M3 Ultra, single core, nice 10):

| arm | ms/shot (batch) | source |
|---|---|---|
| bp30+osd (ldpc C++) | 64.1 | smoke run (1000 shots) |
| beam8 (numpy, ours) | 526 | smoke run (1000 shots) |
| BP+OSD single-shot timing | mean 68.1, p50 3.1, p99.9 275.2 ms | smoke timing arm x50 |

Gate A full cost at the pre-statement budget (≥1e5 shots × 6 points, BP+OSD
arm): ≈ 10.7 single-core hours. [NUMERICAL]

## How to run

```bash
cd physics
# 1e3-shot end-to-end rehearsal (writes campaigns-smoke/<UTC>_<uuid>_<hash>/)
.venv/bin/python -m qldpc_dec.run_gate smoke                    # from qldpc-dec/src
.venv/bin/python -m qldpc_dec.run_gate gateA_accuracy 100000    # gate A LER companion
.venv/bin/python -m qldpc_dec.run_gate gateA_timing 10000       # gate A timing (Table II/III)
.venv/bin/python -m qldpc_dec.run_gate gateB 100000             # gate B beam8 vs bp30+osd
```

Each run writes an immutable campaign snapshot
`campaigns/<UTC>_<uuid8>_<confighash12>/` with `manifest.json` (config, seeds,
package versions, circuit+DEM sha256), `results.json.gz`, `summary.json`,
`inventory.json`; never edited after close. Rehearsals go to `campaigns-smoke/`.

## Layout

- `pre_statement.md` — frozen gate statement (2026-08-29), read first.
- `circuits/` — committed stim circuits (IonQ-basis; provenance in `dem/`).
- `dem/` — frozen canonical DEMs + `PROVENANCE.json` (sha256 of circuits+DEMs).
- `src/qldpc_dec/` — harness package: `code.py` (BB constructor, fallback),
  `circuits.py` (circuit/DEM loader), `dem_matrices.py` ((H,A,p) extraction,
  raw + merged conventions), `harness.py` (BBCodeHarness; shape contract
  (num_shots,num_dets)→(num_shots,num_obs) shared with fss-bb),
  `bp_osd.py`, `beam_search.py`, `nms.py`, `bootstrap.py`, `seeds.py`,
  `runner.py` (campaign snapshots), `run_gate.py` (entry points).
- `campaigns/` — immutable campaign snapshots (18 as of 2026-08-30; see
  `../RESULTS.md` ladder).
- `src/beam_cpp/` — C++ port of beam8 (`beam8.cpp`, `equiv.py`,
  `benchmark.sh`, `README.md`): bit-exact to the numpy decoder on shared
  shots (0/1000 mismatches at p=3e-3, seed 20260830 — report
  `scratch/equiv_p0.003_n1000_s20260830.json`); 9.0–14.5 ms/shot at
  p=3e-3 vs 1650–2698 numpy twin-run = 186–262× (single thread).
  Deterministic (no RNG in the decoder; seeds only feed the stim sampler).
  NOTE: pre-fix equivalence artifacts (`equiv_p0.003_n10_*`, 9/10
  mismatches) are supereded — kept under scratch for the bug ledger.
- `src/gatec_exact/` — Gate C exact degenerate-ML reference
  ([[36,4,4]] full affine coset enumeration, exact Fractions;
  brute-force assertion per run; `benchmark.sh`).
- `src/gateD_cert/` — Gate D certified-decoding experiment (arXiv:2608.25545
  spacetime sampling route: sparse trivial-kernel basis + lightening +
  corrected forward CRN-AIS path weights + paired bootstrap;
  `run_pilot.py --mode mini|gross`, exact 16-class cross-check in
  `exact_crosscheck.py`).
- `scratch/` — non-authoritative exploration.

## Interface shared with fss-bb (agreed 2026-08-29)

`BBCodeHarness.decode(dem, shots[bool,(N,M)]) -> predictions[int,(N,K)]`;
campaign dir convention `<UTC>_<uuid8>_<hash12>`; fss-bb imports from
`physics/qldpc-dec/src` (single implementation; erasure DEMs and per-shot
priors handled on their side).

## Pinned reference numbers from the four papers (first-hand reads, 2026-08-29)

- arXiv:2512.07057 (beam search, IonQ): Table I decoder params (exact);
  Table II mean times (3.55/10.59 ms bp30+osd at p=5e-4/1e-3; beam8 1.627/
  2.318 ms); Table III p99.9 (272.5/289.0 ms bp30+osd; beam8 8.704/11.01 ms);
  Fig. 2 per-basis LER curves; claims "17× LER reduction" (beam64),
  "26.2× p99.9 reduction" (beam8), "5.6× + sub-ms p99.9" (beam32). [REPORTED
  → to be tagged REPRODUCED/REFUTED by gate A/B campaigns]
- arXiv:2506.01779 (Relay-BP, IBM): BP+OSD baseline = ldpc-package
  BP+OSD+CS-10 at 10,000 BP iterations (Sec. "Flexible decoding"); gross-code
  Relay-BP params: gamma0=0.125 first leg (80 iters), then legs of 60 iters
  with γ ∈ [−0.24, 0.66]; XZ-decoding matrices ~1k×9k. Relay-BP itself is
  NOT in this harness's gate set (recorded as candidate future arm).
- arXiv:2510.14060 (GARI): headline (6.70±1.93)e-9 per-round LER at p=1e-3
  under uniform depolarizing circuit noise, [[144,12,12]], 24-decoder
  ensemble (99% CI from ≥100 failures); GARI-NMS normalization 0.96875 (11/16)
  or 0.9921875 (127/128), max 400 iters; eq. (7) per-round conversion used in
  `nms.py`. Stim circuits from Gong et al. arXiv:2403.18901 (Fig. 4c) — NOTE:
  our DEM source is the IonQ circuit family instead; this is a declared
  deviation (gate C compares LER values under a same-dimensional 12-round
  memory DEM, not a byte-identical DEM; flagged in pre_statement.md).
- arXiv:2603.19062 (fair baselines, erasure channel): 200k shots/point,
  seed 12345, BP-OSD order 10 erasure-prior decoding; pseudo-thresholds
  0.3701→0.4706 for N=144→1296, FSS p*∞=0.488±0.001, ν=1.18±0.01. Erasure
  channel itself belongs to the fss-bb target; methodology (bootstrap seeds
  recorded, pseudo-threshold ≠ asymptotic) adopted here.

## Next actions (updated 2026-08-30 after gates A–D work)

1. ~~C++ port decision~~ DONE: `src/beam_cpp/` is bit-exact and 166–262×
   faster — formal per-shot diffs are **0/1000 mismatches at p=3e-3**
   (14.5 ms/shot vs 2698 numpy) and **0/1000 at p=1e-3** (3.31 ms/shot
   vs 548 numpy), seed 20260830, single thread. Beam-arm cost for 1e6
   shots collapses from ≈650 single-core-h (numpy) to **≈0.9
   single-core-h** (3.31 ms × 1e6 = 3310 s).
2. Pinned Gate B is now resolvable: `src/beam_cpp/beam8` arm at ≥1e6
   shots/arm p=1e-3 (pinned circuit `BB_144_144_12_memory_Z_p0.001_
   sr12_ionq.stim`) costs ≈0.9 core-h; the Python BP+OSD arm at 279
   ms/shot costs **≈78 single-core-h** per 1e6 shots (279e3 s) and is
   now the sole gating cost. Consider the per-shot protocol fix
   documented in Gate A before spending it.
3. Gate D gross calibration evidence is **RETRACTED**: `AISEngine.run()`
   incorrectly used an endpoint-only target/base ratio after intermediate
   Metropolis moves. Correct forward AIS evaluates every density-ratio
   increment at the pre-transition state. Verification: a one-bit partition
   estimate is within 0.003% of exact (T=32, K=400k); the corrected d=5 mini
   check is 260/260 exact; and a strict all-16-class [[36,4,4]] run is
   191/200 point-estimate agreement, 185/200 certified, with **185/185
   certified decisions exact and 0 certified-wrong**. The prior
   [[144,12,12]] 200/200 gross result used biased weights and supports no
   certificate claim; rerunning that calibration with corrected AIS is the
   next Gate-D slice.
4. Gate C extension: the [[36,4,4]] exact reference supports code-capacity
   scans (slope/c_eff vs p already computed); X-basis transpose is the
   natural next bounded slice.
5. Original gate-C floor run (≤1.93e-9 per-round GARI-class claim)
   remains parked: honest expectation stays "floor not reached in
   available budget".

## Related

- `../fss-bb/` (same harness interface, erasure channel + finite-size scaling).
- Sibling `../../math/qec/` — algebraic side, do not touch.
