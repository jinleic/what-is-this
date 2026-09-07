# `qlops/` — FTQC resource-estimate arithmetic audit

**Status: BENCHMARK (gate A complete: reproduced; gate B complete: latency
robust, magic-state-swap falsifies comparability; official zero-level source
pinned; R7 smoke PASS is frozen REHEARSAL only; the 71m-shot full run is
unstarted and paper d=7/c≈300 remain NOT-REPRODUCED).**
Reuse published fault-tolerant QC resource/benchmark estimates from their stated
formulas, reproduce their numbers, then run sensitivity sweeps. Pure arithmetic;
minutes per run.

## Reference

- **Primary (read first-hand 2026-08-29, arXiv HTML v2):** Kong, Zhang &
  Chen, *Benchmarking fault-tolerant quantum computing hardware via QLOPS*,
  arXiv:2507.12024v2 [quant-ph, 22 Apr 2026]; ACM TQC 7(2):14, DOI
  10.1145/3797968. Every formula/constant transcribed in
  `src/paper_data.py` (anchors = arXiv HTML LaTeX ids, e.g. `S2.E2` =
  Eq. (2), `S3.T5` = Table 5).
- CS-challenges white paper (arXiv:2601.20247): names evaluation as a
  first-class bottleneck for early FTQC. (Not yet read first-hand.)
- Litinski 2019b (arXiv:1905.06903v3, read 2026-08-29): 15-to-1 protocol
  behind the paper's Table 6 (space formula its Fig. 11; time cost
  6·d_m/(1−p_fail)).
- Itogawa et al., zero-level CCZ (arXiv:2605.21867v1, the only version):
  gate-B replacement constants. **Full text read first-hand 2026-09-03**;
  every constant is located in the body. The
  [official repository](https://github.com/FujitsuResearch/Zero-level_CCZ_Distillation)
  is pinned at commit
  `1b59e223590492e224bd8623a4e0bcba59029e01`, tree
  `9a5d89401fd40554635fb0e9d0b4da818670c2bc`, and codeload tar sha256
  `d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`.
  Its executable artifact is d=3/d=9, while the paper's headline grown
  result is d=7. Paper d=7 and the physical c≈300 claim remain
  `[REPORTED]` / **NOT-REPRODUCED**.

## Transcription (gate A input registry)

All constants live in `src/paper_data.py`, one entry per printed value:
Table 1 (SC hardware, current/future), t_SEC = 0.86/0.40 µs (component sums
re-verified exactly), Table 2a/2b (PyMatching t_r for d = 5..27), Table 3
(NA hardware), move model √(6Δx/a_p), Table 4 (GB codes: p_L, t_r, t_SEC),
Table 5 (QLOPS + density, 6 codes × 3 rows), §3.5 RSA-2048 constants,
Table 6 (18 Litinski 15-to-1 rows), plus the zero-level CCZ constants.
Float policy: exact ints/Fractions where the paper's inputs allow; Decimal
ceilings (never float ceil — see `src/test_reproduction.py::
test_rsa_exact_vs_float_ceil`); printed 4–6-sig-fig values reproduced to
their own printing precision.

## Gate A — verdict: REPRODUCED (81 checks, 0 hard failures)

Run: `python3 src/gate_a_reproduce.py <out.json>` (≈1 s, single core).
Artifacts: `campaigns/20260829T223540Z_dcf693ab/artifacts/`
(`gate_a_results.json`, `gate_a_stdout.txt`, frozen `src/` copies).
Tolerance: ≤5 % relative (the paper states none); every reproduced number
is in fact ≤0.02 % except two `~N` prose roundings.

| quantity | verdict |
| --- | --- |
| Eq. (1) p₀, 6 codes (S3.T5 col 2) | PROVED-reproduced (≤0.02 %) |
| Eq. (2) Q, NA rows, 6 codes | PROVED-reproduced (≤0.02 %) |
| QLOPS density, NA + N=2n convention | PROVED-reproduced (≤0.02 %) |
| SC matched rows, 12 sub-rows, qubit counts N=k(2d²−1) | PROVED-reproduced (≤0.02 %) |
| SC matched rows, current-hw d∈{17,19} (3 Q values) | PROVED-reproduced **[with correction]** — see finding F1 |
| RSA-2048: Q_SC=4.0314e7, 714 019 data qubits, 56.4611, Q_NA=6.8089e6, 0.7626, Eq. (3) 5.244, Eq. (4) 2.4204, ~110/~270 | PROVED-reproduced (≤4.3 %; `~270` is prose) |
| Table 6: 18 unit-qubit counts, 18 total = unit × integer | PROVED-reproduced (exact) |
| Table 6: distillation error, cycles | NOT-REPRODUCED-from-inputs [REPORTED] — produced by the paper's 5-qubit density-matrix sim under unstated noise inputs; cycles inverted to implied p_fail ∈ [0.00320, 0.07024] ([DERIVED]) |
| Fig. 1 fit constants | NOT-REPRODUCED-from-inputs — not printed; only visual read-offs (image in campaign artifacts) |

**Paper-internal findings (evidence in campaign artifacts):**

- **F1 — Table 2(a) t_r appears shifted by one slot for d ≥ 17.** The
  printed d=17 entry (9.7008e-5) is out-of-trend, and Table 5's current-SC
  QLOPS for d=17/19 reproduce the published values exactly (dev ≤ 4×10⁻⁸)
  only when t_r is taken from the *next* d row (d=19→4.8795e-5,
  d=21→6.9903e-5). Using entries as printed: deviations +43 %/+33 % (far
  outside 5 %). [DERIVED]
- **F2 — ceil boundary trap.** [[288,12,18]] current-SC row implies
  ⌈t_r/t_SEC⌉ = 259 while t_r/t_SEC = 2.2188e-4 / 0.86e-6 = 258.0 exactly;
  paper value sits one unit high (0.35 % effect, within tolerance). The
  same one-unit float-ceil error appears in `ceil(1e-5/1e-6) = 11`; our
  harness takes ceilings on exact rationals. [DERIVED]

## Gate B — sensitivity (sweep completed; verdicts below)

Run: `python3 src/gate_b_sweep.py <out.json>`; artifacts
`campaigns/20260829T223540Z_dcf693ab/artifacts/gate_b_results.json`.

- **Latency axis (0.5×–2× all published t_r, 100 evals): ROBUST.** Max
  shift 1.82× (< 2× falsifier), analytically bounded by
  (⌈m·r⌉+d)/(⌈r⌉+d) ≤ (2r+1+d)/(r+d) < 2 for m ≤ 2. Largest jumps at the
  ⌈·⌉ boundaries (e.g. SC [[288,12,18]] at m=0.5: 1.82×). No latency
  factor within ±2× falsifies the benchmark's outputs. [PROVED for the
  implemented formula; sampled over the paper's own 13 (k,d,t_r)
  configurations]
- **Magic-state swap (7-T Litinski 15-to-1 → zero-level CCZ, constants
  `[REPORTED]`, provenance first-hand): COMPARABILITY FALSIFIED on this
  axis.** Accuracy-wise, the zero-level rule p_L = 300p² violates the
  paper's own p_out ≤ p₀ matching rule for 16/18 Table-6 rows at p=1e-3
  and 11/18 at p=1e-4. QLOPS-style headline numbers are therefore not
  invariant under this protocol substitution. [REPORTED inputs; DERIVED
  verdict]

  The frozen campaign's **193.3×–29,946.5×** resource row is a historical
  repo derivation, not a published value. It divides 22 physical qubits ×
  24 circuit layers by Litinski physical qubits × syndrome-extraction
  cycles, so it has no common-unit magnitude interpretation. Revision 4
  preserves that row only to reproduce the frozen artifact.

  **Common-unit accepted-output sensitivity (Revision 4,
  `src/evidence/zero_level_sensitivity.json`, sha256
  `1558bf3ebc493f1c3eff02bba13117c647537b78235a4c2d832ac340881e01a5`).**
  The comparison basis is physical qubits × surface-code
  syndrome-extraction cycles per accepted CCZ:

  - Litinski: 7 × Table-6 unit qubits × the reported cycle count, which
    already includes postselection.
  - Zero level: (22 + three full rotated output-patch footprints) × three
    reported syndrome rounds / acceptance.

  | zero-level basis | acceptance | advantage vs accepted 7-T Litinski |
  |---|---:|---:|
  | + 3 output patches d=3 | 0.40 | 193.0× – 29,077.5× |
  | + 3 output patches d=3 | 0.30 | 144.8× – 21,808.1× |
  | + 3 output patches d=7 | 0.40 | 45.0× – 6,781.6× |
  | + 3 output patches d=7 | 0.30 | **33.8× – 5,086.2×** |
  | + 3 output patches d=7 | 0.90 | 101.3× – 15,258.7× |

  The paper simulates d=3 and d=7. Its 0.30–0.40 values are approximate
  range endpoints not assigned to a specific variant, so cross-pairings are
  explicitly sensitivity arms. The static full-patch footprint is
  `[DERIVED]`, not the paper's own logical-qubit × syndrome-round estimate.
  The worst lower endpoint, 33.77×, remains above the preregistered 2×
  falsifier; the verdict survives. Frozen Gate-A/B JSON still regenerates
  byte-identically.

## Official-source Revision 7 — frozen readiness smoke only

The official
[`FujitsuResearch/Zero-level_CCZ_Distillation`](https://github.com/FujitsuResearch/Zero-level_CCZ_Distillation)
source is pinned exactly at commit
`1b59e223590492e224bd8623a4e0bcba59029e01`, tree
`9a5d89401fd40554635fb0e9d0b4da818670c2bc`, and tar sha256
`d67dbe7482b391984da5e64aeff7668bdaee45c262e6fb2dfd41dfed3d7f3338`.
The repository executes the ungrown d=3 and grown d=9 artifacts. The
paper reports the grown headline at d=7; the two are never equated, no d=7
reconstruction is attempted, and paper d=7/c≈300 remain
**NOT-REPRODUCED**.

**Instrument history (neither predecessor is evidence):**

- R5 campaign
  [`20260904T052209Z_2a0a2723_203c9b1d1900`](campaigns/20260904T052209Z_2a0a2723_203c9b1d1900/status.json)
  is frozen/closed **SUPERSEDED**. It correctly refused before sampling
  because its ungrown p=0 shipped-oracle equality premise was unsupported;
  no R5 result is accepted.
- R6 campaign
  [`20260904T054022Z_5f040b2a_1f64afb2c983`](campaigns/20260904T054022Z_5f040b2a_1f64afb2c983/status.json)
  is frozen/closed **CRASHED**. Its in-memory smoke calculations completed,
  but a results-directory orchestration bug prevented artifact writing.
  The regression was fixed before R7; no R6 result is accepted.
- Canonical R7 campaign
  `20260904T054537Z_818bca48_18f5ac3b849c`, gate
  `zero-level-author-repro-r7`, preregistration sha256
  `99b3985b261dc3de561917736ea25c1e552283a67b2353df50f137a954bb0668`,
  is frozen/closed **REHEARSAL**. Its machine smoke verdict is **PASS**
  across all eight gates: four builder-only ungrown p=0
  invariant/zero-error/replay gates, and four real shipped grown p=0.001
  equality/fault-id/accepted-only scalar-vs-batch/replay gates. See the
  canonical [summary](campaigns/20260904T054537Z_818bca48_18f5ac3b849c/results/summary.md),
  [smoke result](campaigns/20260904T054537Z_818bca48_18f5ac3b849c/results/smoke.json),
  [manifest](campaigns/20260904T054537Z_818bca48_18f5ac3b849c/manifest.json),
  [status](campaigns/20260904T054537Z_818bca48_18f5ac3b849c/status.json),
  and [sha256 ledger](campaigns/20260904T054537Z_818bca48_18f5ac3b849c/sha256s.txt).
  `sha256sum -c` passed all 31 listed payloads. This is launch-readiness
  evidence only, with no scientific-result weight.

The frozen full schedule is 71,000,000 author-count shots across 12 cells.
It remains pending and **has not been launched**. After GB9 releases its
26 threads, run this sequence from `physics/qlops`. First mint a fresh R7
campaign:

```sh
python3 ../../scripts/campaign.py init \
  --gate zero-level-author-repro-r7 \
  --prereg pre_statement.md
```

Substitute the emitted `run` value for `<printed-run-id>`, then launch
against the immutable, hash-verified R7 tar with an explicit one-thread
environment and low priority:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
/usr/bin/nice -n 10 ../.venv/bin/python src/zero_level_repo_repro.py \
  --source-tar campaigns/20260904T054537Z_818bca48_18f5ac3b849c/source/source_tar.gz \
  --run-dir campaigns/<printed-run-id> \
  --mode full
```

The R7 smoke environment was Python 3.14.3, stim 1.16.0, pymatching
2.4.0, and numpy 2.5.2.

## How to run

```
python3 src/gate_a_reproduce.py out.json      # gate A, ~1 s
python3 src/gate_b_sweep.py out.json          # gate B axes
python3 src/zero_level_provenance.py          # Revision 4 sensitivity, <1 s
python3 src/test_reproduction.py              # 7 behavioral tests
```

## Layout

- `campaigns/20260829T223540Z_dcf693ab/` — frozen gate A+B artifacts
  (JSON results, stdout, code copies, Fig. 1 image).
- `campaigns/20260904T052209Z_2a0a2723_203c9b1d1900/` — frozen/closed
  SUPERSEDED R5 refusal; unsupported p=0 oracle, no accepted result.
- `campaigns/20260904T054022Z_5f040b2a_1f64afb2c983/` — frozen/closed
  CRASHED R6 instrument run; results-directory regression, no accepted result.
- `campaigns/20260904T054537Z_818bca48_18f5ac3b849c/` — canonical
  frozen/closed R7 REHEARSAL; eight-gate smoke PASS, readiness only.
- `scratch/` — non-authoritative downloads.
- `pre_statement.md` — pre-registered gates, tolerances, falsifiers
  (written before the run; only tolerance language unchanged).

## Open items

- The official source is now pinned, but the 71,000,000-shot, 12-cell R7
  full run remains pending and unstarted; the completed eight-gate smoke
  is REHEARSAL readiness evidence only. A successful full run can evaluate
  only the repository's executable d=3/d=9 artifact. The paper's d=7
  headline and physical c≈300 claim remain **NOT-REPRODUCED**, and no d=7
  reconstruction is attempted.
- Fig. 1 fit constants and Table 6 sim-level numbers remain
  NOT-REPRODUCED-from-stated-inputs (missing quantitative inputs, recorded
  per pre_statement.md policy — not guessed).

## Related

- `../msd/` is a closed SCOPE-LIMITED-NON-TEST for the zero-level c≈300
  claim, not an input source. Decoder-latency constants are transcribed from
  the QLOPS paper; `../qldpc-dec/` runtime measurements are not integrated.
