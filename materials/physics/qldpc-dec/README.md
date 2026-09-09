# `qldpc-dec/` — cross-paper decoder reproduction harness for BB codes

**Status (2026-09-09): PARKED.** GB9's mini-0 relaunch
[`20260908T194719Z_e9e0a889_e5962c4be23b`](campaigns/20260908T194719Z_e9e0a889_e5962c4be23b/RESULTS.md)
was stopped by owner direction: 16 remaining shards at 30–45 h each exceed
the 72-hour per-task per-mini occupation limit. Administrative close
`CRASHED` (cancellation); no scientific verdict on the partial num_results=32
experiment. GB9 will not resume under this allocation. Gate B remains DECIDED
on the completed paired 1e8-shot instrument
(GB7, `campaigns/20260901T145247Z_b7ea9ac4_ac3f6689e03b`, terminal verdict
FROZEN-NEGATIVE on the beam32 gap hypothesis) and RE-TARGETED by Revision
GB8 (2026-09-02, no sampling). On one shared stream at p=1e-3 Z: beam8 223,
beam32 41, beam64 27 failures per 1e8. beam32 measured expansion **5.44x
[4.17, 7.55]** vs published 5.6x — REPRODUCED (McNemar 4e-44). beam64
measured **8.26x [5.97, 12.65]**; the paper's factor for the configuration
we ran (`beam64_640iters`, num_results=1) is **7.0x**, inside the interval —
**CONSISTENT with the published point value**. The "17x" is the paper's
`beam64_32res_640iters` (num_results=32), never completed: neither reproduced nor
refuted. The 2026-09-01 headline "beam64 17x refuted" is **RETRACTED**
(transcription error in the original 2026-08-29 Gate B pre-statement, corrected by Revision GB8); GB7's frozen `GAP_CONFIRMED` against
1/17 stands as arithmetic but tests a target the paper does not make. GB5a
2e7 independent-arm ladder FROZEN-INCONCLUSIVE (BP+OSD 41 / beam8 42 / beam32
10 / beam64 5). GB6 dense paired mechanism FROZEN-CERTIFIED (p=3e-3:
654/464/145/80 per 2e5, M1-M3 PASS, EXPANSION_EFFECT_PRESENT). Original
Gate-C floor campaign remains PAUSED (≥1.5e10 shots). Historical scalar-prior
BP evidence remains retracted; corrected Gate-C exact comparisons and Gate-D
calibrations are closed (2026-08-31).**
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
  repo loader (arXiv:2506.01779). Published p=5e-4 and p=1e-3 circuits are
  frozen under `circuits/`; the p=3e-4 X and Z circuits are derived separately
  by exact probability rescaling from the corresponding p=1e-3 basis circuit.
  Circuit/DEM hashes and derivations live in `dem/PROVENANCE.json`. Each DEM
  has 936 detectors and 8784 merged error mechanisms; 864 Y-type mechanisms
  remain graphlike hyperedges (documented; this makes the reproduction
  [DERIVED] under the parent README terminology).
- Decoders (all three smoke-validated end-to-end):
  - BP+OSD baseline: `ldpc` 2.4.1 `BpOsdDecoder`, full heterogeneous merged-DEM
    `error_channel`, max_iter=30, ms, osd_cs order 10 (matching
    `stimbposd==0.1.0`, the IonQ paper's pinned bp30+osd implementation).
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
   **Revision GB8 (2026-09-02): the 17× / [11, 22] target was mis-attributed.**
   Our beam64 rung is the paper's `beam64_640iters` (num_results=1), published
   at 7.0× (Section III); 17× is `beam64_32res_640iters` (num_results=32),
   never run. The band is void; no corrected band is registered.
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

- `20260829T225407Z_53700f78_87fb945aabf3` — historical 1e3-shot
  end-to-end rehearsal. Its BP+OSD arm used the now-invalid scalar mean prior;
  beam8 and NMS remain pipeline proofs only, with no gate verdict.
- `20260829T230642Z_755286ce_1f387aff4429` — 5-shot validation that every gate
  point (p ∈ {3,5e-4,1e-3} Z + 1e-3 X, both decoder arms) decodes cleanly,
  including the derived p=3e-4 circuit.

Historical rehearsal costs below used the invalid scalar-prior BP adapter;
they are retained only as bug evidence, not published-baseline measurements:

| arm | ms/shot (batch) | status |
|---|---|---|
| bp30+osd (ldpc C++) | 64.1 | retracted scalar-prior run |
| beam8 (numpy, ours) | 526 | unaffected rehearsal |
| BP+OSD single-shot timing | mean 68.1, p50 3.1, p99.9 275.2 ms | retracted scalar-prior run |

Corrected vector-prior costs and Gate-A results are recorded below.

## How to run

```bash
cd physics/qldpc-dec
# 1e3-shot end-to-end rehearsal (writes campaigns-smoke/<UTC>_<uuid>_<hash>/)
PYTHONPATH=src ../.venv/bin/python -m qldpc_dec.run_gate smoke
PYTHONPATH=src ../.venv/bin/python -m qldpc_dec.run_gate gateA_accuracy 100000
PYTHONPATH=src ../.venv/bin/python -m qldpc_dec.run_gate gateA_timing 10000
PYTHONPATH=src ../.venv/bin/python -m qldpc_dec.run_gate gateB 100000
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
- `campaigns/` — immutable campaign snapshots; see `../RESULTS.md` for the
  authoritative campaign ladder.
- `src/beam_cpp/` — C++ port of beam8 (`beam8.cpp`, `equiv.py`,
  `benchmark.sh`, `README.md`): bit-exact to the numpy decoder on shared
  shots (0/1000 mismatches at p=3e-3, seed 20260830 — report
  `scratch/equiv_p0.003_n1000_s20260830.json`); 9.0–14.5 ms/shot at
  p=3e-3 vs 1650–2698 numpy twin-run = 186–262× (single thread).
  Deterministic (no RNG in the decoder; seeds only feed the stim sampler).
  NOTE: pre-fix equivalence artifacts (`equiv_p0.003_n10_*`, 9/10
  mismatches) are superseded — kept under scratch for the bug ledger.
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

## Pinned reference numbers from the four papers (first-hand reads, 2026-08-29; beam-search LER factors re-read 2026-09-02, Revision GB8)

- arXiv:2512.07057 (beam search, IonQ): Table I decoder params (exact);
  Table II mean times (3.55/10.59 ms bp30+osd at p=5e-4/1e-3; beam8 1.627/
  2.318 ms); Table III p99.9 (272.5/289.0 ms bp30+osd; beam8 8.704/11.01 ms);
  Fig. 2 per-basis LER curves; Section III factors vs bp30+osd at p=1e-3:
  beam8 1.3× (abstract: "same"), beam32 5.6×, beam64_640iters 7.0×,
  beam64_32res_640iters (num_results=32) 17×; "26.2× p99.9 reduction"
  (beam8) and "5.6× + sub-ms p99.9" (beam32). No shots, failure counts or
  error bars are published. Gate-A timing is REPRODUCED; Gate-B below.
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

## Corrected Gate-A/B evidence (closed 2026-08-31)

All corrected BP rows pass an element-for-element assertion that the
`BpOsdDecoder` received the 8784-value heterogeneous merged-DEM channel. The
old scalar-mean campaigns remain immutable and retracted.

### Gate A timing — REPRODUCED

One compiled decoder and 10,000 individual `perf_counter` calls per point,
run after all other decoder campaigns:

| p | mean ours (ms) | paper (ms) | ratio | p99.9 ours (ms) | paper (ms) | ratio | factor-2 gate |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 5e-4 | 3.842 | 3.55 | 1.082 | 303.767 | 272.5 | 1.115 | pass |
| 1e-3 | 10.995 | 10.59 | 1.038 | 307.796 | 289.0 | 1.065 | pass |

The timed decoder is the published OSD-CS10/order-10/max-iter-30
configuration, not OSD0. A corrected-prior 3,000-decode OSD0 diagnostic
produced mean/p99.9 = 1.208/5.986 ms at p=5e-4 and 1.875/6.147 ms at p=1e-3,
far below the published tails. This retracts the historical inference, made
under the invalid scalar prior, that the paper's timing row used OSD0
semantics. Diagnostic: `scratch/gatea_vector_prior_osd0_probe.json`.

The source timing campaign is
`campaigns/20260831T110541Z_e3e01d3d_774ea68fcaa7/`. Its summary omitted the
5e-4 verdict because one reference key was formatted as `"5e-04"` while the
measurement key was `"0.0005"`; no timing row was affected. Revision GA3
prospectively froze a no-rerun, hash-checked verdict amendment at
`campaigns/20260831T111233Z_755b3791_641ea5506f3b/`.

### Gate A accuracy companion — DESCRIPTIVE

The authoritative six rows (X/Z at p=3e-4, 5e-4, and 1e-3) each observed
0 failures in 100,000 shots. Each arm's 95% Wilson upper bound is
3.2011e-6 per round. No source publishes a matching accuracy target, so this
is not a reproduction verdict. The first corrected campaign had silently
loaded the derived Z circuit for its p=3e-4 row labeled X; Revision GA2
replaced only that row using the separately derived X circuit and the same
frozen seed. Authoritative amendment:
`campaigns/20260831T103828Z_160bc8cc_05a5ebaec745/`.

### Gate B beam8 band — INCONCLUSIVE

At p=1e-3 Z memory, corrected OSD-CS10/order-10/max-iter-30 BP+OSD had
4 failures in 1,000,000 shots (3.3333e-7 per round), and bit-exact C++ beam8
had 1 failure in 1,000,000 shots (8.3333e-8 per round). The point ratio
beam8/BP+OSD is 0.25, but the frozen
100,000-draw independent-binomial bootstrap interval is [0, 2]. It overlaps
the published-equivalence band [0.87, 1.15], so the pre-stated three-way
decision is **INCONCLUSIVE**, not a confirmation or refutation. This closes
only the beam8 band, not all of Gate B. Beam32 and beam64 remain unrun because
the frozen protocol gated them on beam8 passing.
The 20-shard BP aggregate timings were concurrent throughput measurements,
not Gate-A latency evidence. Authoritative campaign:
`campaigns/20260831T091453Z_fd47e8d5_f3b51514a09e/`.

## Remaining actions

1. The original Gate-C floor campaign for the ≤1.93e-9 per-round GARI-class
   claim remains parked: the available budget cannot reach that floor.
2. The GB5a ladder (below) is the prospective revision that authorized the
   beam32/beam64 rungs; the earlier frozen ladder correctly stopped at beam8.
3. A future Gate-D extension should force nonzero candidate classes before
   making a contested-class or global-ML claim.


## Gate B ladder (GB5a) — FROZEN INCONCLUSIVE (2026-09-01)

`campaigns/20260901T040018Z_55e4ac01_f81265905380/` — 20,000,000 shots per
arm at p=1e-3 Z, pre-statement revisions GB1–GB5a, closed
**FROZEN-INCONCLUSIVE** (529 artifacts hashed into `sha256s.txt`; control
manifest merged post-hoc per the Amendment in `pre_statement.md` because the
ladder had run inside a live claim on a sibling target).

Independent-arm results (each vs the shared BP+OSD denominator, GB2-semantics
bootstrap, 100,000 draws):

| rung | failures/2e7 | raw rate (per round) | ratio vs BP+OSD | ratio CI95 | published band | verdict |
|---|---|---|---|---|---|---|
| BP+OSD | 41 | 1.708e-7 | — | — | — | — |
| beam8  | 42 | 1.750e-7 | 1.024 | [0.660, 1.593] | [0.87, 1.15] | INCONCLUSIVE |
| beam32 | 10 | 4.167e-8 | 0.244 | [0.102, 0.452] | [1/7, 1/4.5] | INCONCLUSIVE |
| beam64 | 5  | 2.083e-8 | 0.122 | [0.025, 0.263] | [1/22, 1/11] (void, GB8) | INCONCLUSIVE |

Ladder outcome: `INCONCLUSIVE_no_band_excluded_or_recovered_for_every_rung`.
No band is recovered and none is excluded; both arms of every ratio carry
independent Poisson noise, so the measured 4.1x (beam32) and 8.2x (beam64)
reductions cannot be separated from the published 5.6x / 7.0x (the beam64
band above was built on the mis-attributed 17x; its verdict is unchanged
because the CI contains 1/7.0 as it contains every candidate — Revision GB8).

Paired instrument on the SAME shared stream (bp30+osd seed), which is the
ship-embedded mechanism data for GB7's prefix continuity:

| pair | b (first fails, second ok) | c (second fails, first ok) | both fail | McNemar exact p |
|---|---|---|---|---|
| bp30+osd vs beam8  | 31 | 30 | 10 | 1.0 |
| bp30+osd vs beam32 | 37 | 4  | 4  | 1.0e-7 |
| bp30+osd vs beam64 | 37 | 1  | 4  | 2.8e-10 |
| beam8 vs beam32    | 35 | 3  | 5  | 6.7e-8 |
| beam8 vs beam64    | 36 | 1  | 4  | 5.5e-10 |
| beam32 vs beam64   | 3  | 0  | 5  | 0.25 |

Paired failure counts: 41 / 40 / 8 / 5 (BP+OSD / beam8 / beam32 / beam64).
On identical shots the R32 = 8/40 = 0.200 and R64 = 5/40 = 0.125; at 2e7 the
sample is too small to exclude the published reciprocals (0.1786, 0.0588),
which is precisely the question GB7 escalates.

## Gate B mechanism companion (GB6) — FROZEN CERTIFIED (2026-09-01)

`campaigns/20260901T135902Z_71377291_ddeb5cde09fe/` — 200,000 shots at
p=3e-3 Z on ONE shared seeded stream (argrescale circuit), four decoders,
revisions GB6 + GB5a. Failures: BP+OSD 654, beam8 464, beam32 145,
beam64 80. All three preregistered clauses PASS: M1 expansion monotonicity
(464 >= 145 >= 80 with the beam8/beam64 ratio CI excluding 1), M2 sign
stability (beam8/BP+OSD ratio < 1 with CI excluding 1, matching the frozen
p=1e-3 direction), M3 ladder sanity (80 <= 145). Verdict:
**EXPANSION_EFFECT_PRESENT** — the width effect is structural, not a
sparse-count artifact. Closed **FROZEN-CERTIFIED** (128 artifacts hashed).

## Gate B paired escalation (GB7) — FROZEN NEGATIVE on the beam32 gap; beam64 GAP_CONFIRMED against a mis-attributed 17x, re-targeted to 7.0x by GB8 (2026-09-01 / 2026-09-02)

`campaigns/20260901T145247Z_b7ea9ac4_ac3f6689e03b/` — pre-statement Revision
GB7 (recorded pre-decode). 100,000,000 shots of the frozen bp30+osd stream
(sampling_seed(20260829, 1e-3, "Z", "bp30+osd")), decoded by beam8/32/64 at
Table-I parameters. Continuity: the first 2e7 shots, all three beam
prediction files, and all three failure masks byte-match the frozen GB5a
paired artifacts (prefix failure counts 40/8/5 recovered exactly). Closed
**FROZEN-NEGATIVE** (30 artifacts hashed): the verdict label keys on the
beam32 gap hypothesis, which is refuted; the beam64 gap is confirmed and is
reported verbatim below.

Paired four-cell tables against beam8 (reference) on identical shots:

| rung | f(rung)/1e8 | both fail | beam8 fails, rung ok (b) | beam8 ok, rung fails (c) | McNemar exact p |
|---|---|---|---|---|---|
| beam32 | 41 | 30 | 193 | 11 | 4.0e-44 |
| beam64 | 27 | 21 | 202 | 6  | 5.2e-52 |
| beam64 vs beam32 | — | 27 | 14 | 0 | 1.2e-4 |

beam8 = 223/1e8. beam64's failures are a strict subset of beam32's (c=0).

Pre-frozen decision rule applied (four-cell multinomial percentile bootstrap,
100,000 draws, seed derive_seed(20260829, "gb7-paired-ratio", rung, 1e8)):

| rung | R = f(rung)/f(beam8) | R CI95 | published 1/factor | width effect | CI excludes 1 | published excluded | measured factor [CI95] | published | outcome |
|---|---|---|---|---|---|---|---|---|---|
| beam32 | 0.1839 | [0.1325, 0.2400] | 0.1786 | YES | YES | **NO** | **5.44x [4.17, 7.55]** | 5.6x | NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE |
| beam64 | 0.1211 | [0.0791, 0.1674] | 0.0588 | YES | YES | **YES** | **8.26x [5.97, 12.65]** | 17x | **GAP_CONFIRMED** |

Reading (as frozen, 2026-09-01): the beam32 rung reproduces the paper's 5.6x
quantitatively — the published value sits inside a paired interval only 0.11
wide. The beam64 rung excludes 1/17 at 95%. The in-harness width expansion
from 32 to 64 is real (beam64 fixes 14 of beam32's 41 failures with zero
reverses) but roughly 1.5x.

### Revision GB8 re-target (2026-09-02, no sampling) — RETRACTION of "17x refuted"

First-hand re-read of arXiv:2512.07057 (v1 and v2 identical): Table 1 has
four rows. Section III assigns 5.6x to `beam32_340iters`, **7.0x to
`beam64_640iters` (num_results=1)** and 17x to `beam64_32res_640iters`
(num_results=32, a different termination rule: collect 32 valid solutions,
return the minimum-weight one). Our beam64 rung is (64, 40, 30, 20,
num_results=1) = `beam64_640iters`. The original Gate B section of the
pre-statement (2026-08-29) transcribed the abstract's "beam width of 64
achieves 17x" onto the wrong row; GB5a's band [11, 22] and GB7's 1/17 target
inherited it.

Re-applying GB7's own three-condition rule to the frozen beam64 interval with
the corrected factor (`src/gateb_gb8_retarget.py` → `src/evidence/gb8_retarget.json`,
sha256 `d432ea39…c8dd1`; re-derived by the verifier):

| rung | paper row (num_results) | published factor | 1/factor | frozen R CI95 | published excluded | measured factor [CI95] | outcome |
|---|---|---|---|---|---|---|---|
| beam32 | beam32_340iters (1) | 5.6x | 0.1786 | [0.1325, 0.2400] | NO | 5.44x [4.17, 7.55] | NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE |
| beam64 | beam64_640iters (1) | **7.0x** | 0.1429 | [0.0791, 0.1674] | **NO** | 8.26x [5.97, 12.65] | NO_GAP_WIDTH_EFFECT_CONFIRMED_AND_INTERVAL_EXCLUDES_ONE |
| — | beam64_32res_640iters (32) | 17x | — | — | — | NOT TESTED | — |

Verdict: **beam64_640iters is CONSISTENT with the published point value**;
the paper publishes no shots, failure counts or error bars, so no tighter
statement is available. The 17x claim is untested (the C++ binary pins
num_results=1; a `beam64_32res` rung needs a binary change, a new equivalence
gate and its own pre-registration). Recorded source conflict, not resolved:
the abstract says beam8 has "the same" LER as BP-OSD, Section III says 1.3x;
the 2026-08-29 pre-statement and GB7 use the abstract reading (our GB5a: 1.024 [0.660, 1.593]). Under the
1.3x reading the paper-implied factors vs beam8 are 4.31x (inside our beam32
interval) and 5.38x (below our beam64 lower bound 5.97x — the measured
benefit exceeds the implied one); two curve read-outs with no stated
uncertainty support no verdict either way.

Scope: this is the paired instrument against beam8 on a single Z-basis
stream at p=1e-3 with the pinned IonQ circuit and DEM; the beam8-vs-BP+OSD
anchor remains the frozen GB5a 2e7 result (paired 40 vs 41, McNemar 1.0).
No frozen artifact, rule or verdict label changes. Independent acceptance
(raw artifacts re-derived, frozen hashes re-checked, decisions
re-bootstrapped, GB8 re-target re-derived from the frozen summary):
`scratch/verify_gb8.log`, QLDPC_GATEB_ACCEPTANCE_PASS, 87 PASS / 0 FAIL /
0 SKIP.

## Gate B — the paper's actual 17x: beam64_32res (GB9) — stopped, unresolved

`campaigns/20260902T114747Z_69eba724_61e9fe30ba5d/` — pre-statement Revision
GB9 (recorded 2026-09-02T10:56Z, before init and before any sample). Decoder
under test: arXiv:2512.07057 Table 1 row 4, `beam64_32res_640iters` = (64,
40, 30, 20, **num_results=32**), the configuration the paper credits with
17x (Section III) and which GB8 showed had never been run.

What was built and proven before launch:

- Python reference `src/qldpc_dec/beam_search.py` made faithful to Algorithm
  3 step 3 for num_results>1 (seed solution inserted, search continues); the
  weight is a sequential index-order float64 sum; num_results=1 unchanged.
- `src/beam_cpp/beam8.cpp` gained `--num-results=K`, built to a NEW binary
  `src/beam_cpp/beam_nr_cpp` (sha `297db977a57b…`); the frozen `beam8_cpp`
  (`a480050b2041…`) is untouched and still pinned by GB5a/GB6/GB7.
- Equivalence gate (`src/gb5_equiv_width.py --rung beam64_32res --binary
  beam_nr_cpp`): **2000/2000 identical at p=1e-3 and 300/300 at p=3e-3** vs
  the Python reference (21 s/shot; sharded over 12-14 workers). K=1
  regression on the six frozen gate points with the new binary: zero
  mismatches and C++ prediction files **byte-identical** to the frozen
  binary's on every point (`src/beam_cpp/evidence/gb5eq_*_beam_nr_cpp.json`).
- Cost measured: num_results=32 is 188 ms/shot/thread (7.9 ms/shot on 24
  threads) — 42 BP runs and 465 BP iterations per shot, since every shot
  enters the search once the seed solution only counts as one result.
- Sizing from frozen GB7 rates (`src/evidence/gb9_sizing.json`, sha
  `5777193c…`): N = 1e8 pinned — P(refute 17x | no better than
  beam64_640iters) = 0.88, P(exclude the K=1 baseline | true 17x) = 0.915;
  the +-30% band [13.1x, 22.1x] is unreachable at any feasible N (needs
  ~5.6e8) and is reported, not decisive.
- Design: the ENTIRE frozen GB7 stream (resampled; byte-identity is a hard
  void), beam8 re-decoded by the new binary (must equal GB7's frozen
  predictions byte-for-byte: a 1e8-shot K=1 regression), beam64_32res
  decoded shard-wise (20 x 5e6; assembled file byte-identical to a single
  run, demonstrated pre-launch), the frozen GB7 beam64_640iters mask read for
  the same-shot comparison.
- Known differences from the authors' own code (github
  ionq-publications/BeamSearchDecoder): a converged-branch skip and a
  Tanner-degree-<=2 exclusion, both absent from the published Algorithm 3 —
  pre-registered as the first suspects if 17x is not reproduced.

Decision rule, verdict mapping and reporting are frozen in Revision GB9.
The initial campaign closed CRASHED on 2026-09-08. The mini-0 relaunch
[`20260908T194719Z_e9e0a889_e5962c4be23b`](campaigns/20260908T194719Z_e9e0a889_e5962c4be23b/RESULTS.md)
was canceled during shard 04 on 2026-09-09 for projected runtime. Both are
failed execution records, not scientific negatives; the partial-sample rule
still forbids a scientific verdict. No new factor or interval is claimed.


## Related

- `../fss-bb/` (same harness interface, erasure channel + finite-size scaling).
  Sibling `../../math/qec/` — algebraic side, do not touch.
