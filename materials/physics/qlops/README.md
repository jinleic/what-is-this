# `qlops/` — FTQC resource-estimate arithmetic audit

**Status: BENCHMARK (gate A complete: reproduced; gate B complete: latency robust,
magic-state-swap falsifies comparability; msd-dependent parts await `../msd/`).**
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
- Itogawa et al., zero-level CCZ (arXiv:2605.21867, abstract read
  2026-08-29): gate-B replacement constants, [REPORTED] until `../msd/`
  gate A verifies c≈300 within 2×.

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
| Table 6: distillation error, cycles | NOT-REPRODUCED-from-inputs [REPORTED] — produced by the paper's 5-qubit density-matrix sim under unstated noise inputs; cycles inverted to implied p_fail ∈ [4e-4, 3.7e-3] ([DERIVED]) |
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
  [REPORTED] pending `../msd/`): COMPARABILITY FALSIFIED on this axis.**
  Space-time per CCZ drops 193×–29 944× (ratio range 3.3e-5–5.2e-3,
  T→CCZ scaling ×7 stated); accuracy-wise the zero-level rule p_L=300p²
  violates the paper's own p_out ≤ p₀ matching rule at p=1e-4 for 11/18
  Table-6 targets (at p=1e-3, 16/16). QLOPS-style headline numbers are
  therefore not invariant under this protocol substitution. [REPORTED
  inputs; DERIVED verdict]

## How to run

```
python3 src/gate_a_reproduce.py out.json      # gate A, ~1 s
python3 src/gate_b_sweep.py out.json          # gate B axes
python3 src/test_reproduction.py              # 6 behavioral tests
```

## Layout

- `campaigns/20260829T223540Z_dcf693ab/` — frozen gate A+B artifacts
  (JSON results, stdout, code copies, Fig. 1 image).
- `scratch/` — non-authoritative downloads.
- `pre_statement.md` — pre-registered gates, tolerances, falsifiers
  (written before the run; only tolerance language unchanged).

## Open items

- `../msd/` gate A must verify c≈300 (2×) before the zero-level column of
  gate B is promoted from [REPORTED] to [REPRODUCED]; verdicts above are
  stated conditionally on those constants.
- Fig. 1 fit constants and Table 6 sim-level numbers remain
  NOT-REPRODUCED-from-stated-inputs (missing quantitative inputs, recorded
  per pre_statement.md policy — not guessed).

## Related

- Consumes magic-state constants from `../msd/`; decoder latency data from
  `../qldpc-dec/` runtime measurements (not yet integrated).
