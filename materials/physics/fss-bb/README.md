# `fss-bb/` — BB-code finite-size scaling and the ν anomaly

**Status: ACTIVE — Gate-A 20k-shot fallback campaign RUNNING, detached
(launched 2026-08-29 23:50 UTC, hub process `fssbb-gateA-20k`, pid 91617,
campaign dir `campaigns/20260829T235018Z_6df86846` — hands-off until the
runners exits and renames it to `<UTC>_<uuid>_<sha12>`; ETA ~8–10 h).
Survives this session; monitor: `hub logs fssbb-gateA-20k` / wait-for-exit.**
Approved scope (Main, 2026-08-29): 20,000 shots/eval — **caveat: the 20k vs 200k choice changes bootstrap CI width
(~3.2× wider), NOT the estimator chain**, which is the paper's verbatim
(see pre_statement.md). The 200k-shot variant (~50 h) stays PARKED pending
explicit user approval — that is the named next action after the 20k
result. Gate B remains parked with measured feasibility.

Independent finite-size-scaling study of bivariate-bicycle codes on the
erasure channel, attacking the open question left by arXiv:2603.19062: is the
BB-under-BP+OSD critical exponent ν ≈ 1.18 (vs random-bond Ising 4/3) a
finite-size artifact or a distinct effective universality class?

## Reference numbers (from arXiv:2603.19062v3 — read first-hand 2026-08-29,
HTML v3, all cited sections read directly)

Single-author paper (Tushar Pandey, Texas A&M; received 19 Mar 2026, v3 26
Apr 2026):

- Erasure channel, N = 144–1296: (L,M) = (12,6), (18,9), (24,12), (30,15),
  (36,18); A = x³+y+y², B = y³+x+x², Hx = [A|B], Hz = [Bᵀ|Aᵀ]; 200k
  shots/evaluation, seed 12345, bootstrap 95% CIs.
- Pseudo-thresholds (WER = 0.10, Sec 2.4, adaptive bracket-from-0.38 + Illinois
  regula falsi to bracket < 5e-4, ⩽10 evals; 5000-parametric-bootstrap CIs):
  **0.3701 → 0.4386 → 0.4453 → 0.4674 → 0.4706** (Table 1).
- FSS (Sec 2.5): WER ≈ f[(p−p*_inf)·N^{1/ν}], f = degree-3 polynomial,
  window |p−p*| < 0.06, 500 bootstrap datasets; **p*_inf = 0.488 ± 0.001**,
  **ν = 1.18 ± 0.01**; window sensitivity ⩽0.004 (0.04/0.08),
  sys ≲ 0.005; linearized companion (p*(N) ≈ p*_inf + c·N^{−1/ν})
  gives 0.490 ± 0.003.
- Their Sec 4.4 verdict on ν: gap to 4/3 "could be finite-size or
  corrections-to-scaling … or a different effective universality class … we
  leave it open" — hence Gate B.
- [REPORTED] (paper's own claims, not verified by us): BB N=1296 beats toric
  on 12× normalized overhead; uninformed MWPM = random guessing.

## Gates (exact criteria also in pre_statement.md)

1. **Gate A (reproduce):** re-run the erasure-channel FSS with an independent
   harness (own scripts, same sizes, 200k-shot floor, THEIR estimator
   definitions verbatim — see pre_statement.md for the verbatim citations).
   CONFIRM = ≥4/5 per-size p*(N) CI-overlap AND |p*_inf − 0.488| ≤ 0.006
   (their stat ± sys envelope) AND |Δν| ≤ 0.15; else REFUTE-claim flag.
2. **Gate B (decide):** add (42,21) N=1764 (K=8) and (48,24) N=2304 (K=16;
   exact GF(2) rank, corrected from an earlier Ki=12 guess) if one-core
   feasible within a 48 h nice-10 budget (feasibility measured, see
   campaigns-smoke/feasibility_bench.json); test whether ν drifts toward 4/3.
   |Δν| ≥ 0.1 with disjoint CIs = drift (artifact-leaning); stable = evidence
   for decoder-dependent universality. Hard stop conditions pre-committed.

## Status detail (2026-08-29 session)

- [x] arXiv:2603.19062v3 read first-hand; estimator definitions extracted
      verbatim into `pre_statement.md` (channel, decoder, pseudo-threshold
      search, FSS ansatz, windowing, CI method, sizes, shots, seed).
- [x] `src/` pipeline built: BB construction (K verified vs Table 1 for all
      five sizes), erasure sampler, per-shot-prior BP+OSD, logical-failure
      tests, FSS fitter, bootstrap CIs, Illinois pseudo-threshold search,
      DEM export, campaign runner.
- [x] **FSS module validated on synthetic data with known exponents**
      [REPRODUCED, artifact `campaigns-smoke/fss_module_validation.json`]:
      collapse estimator recovers (p0,ν0) = (0.55, 1.35) → 0.5499
      [0.5499, 0.5500], 1.349 [1.347, 1.351] (95% bootstrap CI); out-of-class
      (logistic, (0.61, 1.05)) → (0.6100, 1.050) — < 1% model bias. The
      estimator needed a deterministic multi-start init; the pre-fix
      single-init version FAILED this validation and was fixed before any
      real sweep (history in pre_statement.md / git).
- [x] Decode smoke, N=144, 1e3 shots: WER(0.20)=0.000, WER(0.30)=0.006,
      WER(0.35)=0.051, WER(0.37)=0.116, WER(0.40)=0.259, WER(0.45)=0.652;
      0.06-0.67 ms/shot both sectors — behavior concordant with
      Table-1 p\*(144)=0.3701.
- [x] Campaign-runner smoke (`campaigns/20260829T223700Z_c7e131a5_7f337723833e/`,
      N=144, 2000-shot floor, tag smoke): p\* = 0.3720, CI95 (0.3715, 0.3719)
      [NUMERICAL] — consistent with the paper value at the wider shot
      noise. Full artifacts (manifest, results.csv, thresholds.json,
      fss_fit.json, dem export).
- [x] **Gate-A fragment at the paper floor [REPRODUCED]**
      (`campaigns/20260829T224546Z_987904a9_40c474752126/`, 2026-08-29):
      N=144, 200k shots/eval, exact paper estimator —
      **p\* = 0.37032, 95% CI (0.37004, 0.37101)** vs their Table-1
      0.3701 [0.3697, 0.3703]: CIs overlap. Bracket closed to 2.0e-4 < 5e-4
      (Illinois), 6 evals ≤ 10, wall 326 s. Per-point seed lineage + versions
      pinned in manifest.json (python 3.14.3, numpy 2.5.2, scipy 1.18.1,
      stim 1.16.0, ldpc 2.4.1).
- [x] Gate-B feasibility measured and decided
      (`campaigns-smoke/feasibility_bench.json`): per-shot cost ≈ N²;
      Gate-A floor run ≈ 50 h core; Gate-B both new sizes ≈ 221 h ≫ 48 h
      ceiling ⇒ Gate B PARKED at paper floor, re-scope options with numbers
      recorded in pre_statement.md (recommended: 20k-shot variant).
- [x] Interface agreement with qldpc-dec recorded (below).
- [ ] **Full Gate-A campaign (5 sizes × 200k shots ≈ 50 h core) — PENDING
      budget decision**: a scheduled multi-session job; the 20k-shot
      alternative (≈ 5 h, ~3.2× wider CIs) is the fallback — see
      pre_statement.md feasibility table.
- [ ] Gate B PARKED pending Gate-A outcome + pre-statement amendment for
      the reduced-shot re-scope (options with numbers in pre_statement.md).
- [x] Cross-check against the paper's own data [DERIVED]: our
      `fit_linearized_fss` applied to their Table-1 ladder
      (0.3701/0.4386/0.4453/0.4674/0.4706, ν=1.18) returns
      p\*_inf = 0.4896 — their published linearized companion is
      0.490 ± 0.003 (Fig. 2). Our estimator reproduces their intermediate
      result from their numbers; the collapse estimator remains the
      independent arm of Gate A.

## Shared harness interface (agreed with qldpc-dec, 2026-08-29, via hub)

1. Single implementation lives in `../qldpc-dec/src/`; `fss-bb` imports it
   (never edits it). Until that harness lands, `fss-bb/src/decode.py` is the
   erasure-path implementation (documented above); planned merge = thin
   adapter, already shape-compatible.
2. Decode signature (their `BBCodeHarness`): `(dem: stim.DetectorErrorModel,
   shots: bool[num_shots, num_dets]) -> int[num_shots, num_obs]`.
   `fss-bb/src/decode.py:decode_dem_shots` implements exactly this contract.
3. Primary erasure path is NOT DEM-based (per-shot erasure priors 0.5/1e-10
   cannot live in a static DEM — paper Sec 2.2 updates posteriors per shot):
   direct sampling + `ldpc.BpOsdDecoder.update_channel_probs`. DEMs are
   exported for interop per convention: `campaigns/<UTC>_<uuid>_<hash>/dem/
   N{N}_p{p}.dem` (marginal-flattening semantics documented in
   `src/dem_export.py`; X/Z-erasure correlations are not DEM-representable —
   stated, not hidden).
4. Campaign naming: `<UTC-timestamp>_<uuid8>_<shorthash12>/` under each
   target's own `campaigns/` (matches physics/README.md contract; shorthash =
   sha256(results.csv+manifest.json)[:12]).

## Layout

- `pre_statement.md` — gate statement: verbatim estimator definitions,
  pass/fail criteria, feasibility budget, stop conditions.
- `src/bb_codes.py` — Hx/Hz construction, GF(2) tools, logical-operator
  quotient bases; asserts K against Table 1.
- `src/erasure_channel.py` — channel sampler + sector failure tests.
- `src/decode.py` — BpOsdConfig (paper-verbatim), per-shot-prior batch
  decode (timing-instrumented), DEM-shots adapter.
- `src/fss_fit.py` — FSS ansatz fit (multi-start), 500-dataset bootstrap,
  linearized companion, Illinois pseudo-threshold search + bootstrap.
- `src/dem_export.py` — marginal erasure DEM export (interop only).
- `scripts/validate_fss_synthetic.py` — synthetic-exponent validation
  (T1 exact / T2 model-bias); artifact in campaigns-smoke/.
- `scripts/bench_feasibility.py` — per-shot cost vs N incl. Gate-B sizes;
  artifact campaigns-smoke/feasibility_bench.json.
- `scripts/campaign_runner.py` — threshold/grid modes, artifacts per layout.
- `campaigns/` — frozen run snapshots (one smoke run present).
- `campaigns-smoke/` — validation artifacts + crashed-runner debris (named).

## How to run

```bash
cd physics/fss-bb
V=../.venv/bin/python
# validation (must pass before any real sweep)
OMP_NUM_THREADS=1 $V scripts/validate_fss_synthetic.py
# feasibility table
OMP_NUM_THREADS=1 nice -n 10 $V scripts/bench_feasibility.py
# smoke campaign (small)
OMP_NUM_THREADS=1 $V scripts/campaign_runner.py --sizes 12x6 --shots 2000 \
    --tag smoke --n-boot-thresh 500 --n-boot-fss 100
# Gate A campaign (paper floor; ~hours — see pre_statement budget)
OMP_NUM_THREADS=1 nice -n 10 $V scripts/campaign_runner.py \
    --shots 200000 --tag gateA --seed-base 20260829
```

## Related

- Shares the decoder harness with `../qldpc-dec/` (interface above; merge
  plan: swap `src/decode.py` internals for BBCodeHarness behind the same
  shape contract).
