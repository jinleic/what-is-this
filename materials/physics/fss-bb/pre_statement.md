# fss-bb — pre-statement (v0.1, 2026-08-29, before any real campaign run)

Owner: fss-bb agent. Source paper: arXiv:2603.19062v3, "Fair Decoder Baselines
and Rigorous Finite-Size Scaling for Bivariate Bicycle Codes on the Quantum
Erasure Channel", Tushar Pandey, received 19 Mar 2026, last revised 26 Apr
2026 (v3). Read first-hand 2026-08-29 (HTML v3, full text, all sections
cited below read directly — not abstracts or summaries).

---

**AMENDMENT A1 (2026-08-29, BEFORE any 20k data was used in a fit):**
Main approved the 20,000-shot/eval fallback for the full 5-size Gate-A run
(budget decision; the 200k-shot variant stays parked pending user
approval). EFFECT: bootstrap CI widths widen ~3.2× (1/sqrt(10) shot-noise
scaling); the estimator chain, windowing, bootstrap methods, seeds, and
pass/fail criteria above are UNCHANGED. Per-size CI-overlap checks against
the paper's Table-1 CIs may therefore fail on CI width alone — in that
case the gate verdict is recorded as INCONCLUSIVE-UNDER-BUDGET, never as a
refutation, and the 200k variant remains the decisive arbiter.

---

## Gate A — independent reproduction of the erasure-channel FSS

**Exactly reproduced quantities (their estimator, verbatim restatement):**

1. **Channel (Sec 2.2):** each qubit erased independently w.p. `p`; erased
   qubits receive independent Pauli X and Z errors with Pr = 0.5; decoding is
   BP-OSD with erasure-informed priors (erased qubits: channel probability
   0.5; non-erased: 1e-10); logical error if either X or Z sector fails;
   WER = fraction of trials with a logical error in either sector.
2. **Decoder (Sec 2.3):** "min-sum BP, OSD-CS post-processing, OSD order 10,
   maximum 50 BP iterations" (bposd package = ldpc BpOsdDecoder core; we call
   the core directly since bposd is not installed — the ldpc version pinned
   in the campaign manifest).
3. **Pseudo-threshold (Sec 2.4):** "the erasure rate at which WER = 0.10 for
   a given finite code size", located by adaptive search: step upward from
   p = 0.38 in increments of 0.04 until WER > 0.10 (bracketing); refine the
   bracket using regula falsi (false-position) with the Illinois
   modification; stop when bracket width < 5e-4 or after at most 10 total
   evaluations; report the final linear interpolant. Every evaluation uses
   200,000 shots. Uncertainty: 5,000-iteration parametric bootstrap
   resampling from the binomial at the bracketing points (95% CIs).
4. **FSS ansatz (Sec 2.5):** WER(p, N) ≈ f[(p − p*_inf) · N^(1/ν)], f
   approximated by a **degree-3 polynomial** in the scaling variable x =
   (p − p*_inf)·N^(1/ν); minimize residual sum of squares over
   (p*_inf, ν) using per-size WER data within window **|p − p*| < 0.06**
   (p* = that size's pseudo-threshold).
5. **Uncertainty (Sec 2.5):** bootstrap CIs from **500 resampled datasets**;
   fit-window sensitivity at |p − p*| < 0.04 and < 0.08 (they report p*_inf
   shifts ≲ 0.004, systematic ≲ 0.005); companion linearized fit
   p*(N) ≈ p*_inf + c·N^(−1/ν).

**Sizes/loads (Sec 2.1, Table 1):** family A = x³ + y + y², B = y³ + x + x²,
Hx = [A|B], Hz = [Bᵀ|Aᵀ], N = 2LM; (L,M) = (12,6) N=144 K=12,
(18,9) N=324 K=8, (24,12) N=576 K=16, (30,15) N=900 K=8, (36,18) N=1296
K=12. 200,000 shots per evaluation; seed 12345 (their runs); our runs use
recorded per-point seed lineage (sha256(seed_base|tag|LxM|p)).

**Their reported constants (to compare against):** p*(N) ladder 0.3701 →
0.4386 → 0.4453 → 0.4674 → 0.4706 (Table 1); FSS p*_inf = 0.488 ± 0.001
(95% boot CI), ν = 1.18 ± 0.01 (95% boot CI); window-shift systematic ≲
0.005 on p*_inf; linearized companion p*_inf = 0.490 ± 0.003.

**Gate A pass/fail criterion (falsification):** CONFIRM per-size if our 95%
bootstrap CI on p*(N) overlaps theirs for ≥ 4 of 5 sizes AND our p*_inf is
within their quoted stat+sys envelope (|p*_inf_ours − 0.488| ≤ 0.006) with ν
consistent within the (large) systematic uncertainty of both pipelines
(|Δν| ≤ 0.15). REFUTE-claim flag if either fails after seed-variation and
window checks. Either outcome is a deliverable result.

**Independence statement:** our channel sampler, decoder driver, logical-op
construction, threshold search, and FSS fitter are written from scratch in
`fss-bb/src/` from the paper's Sec-2 text; no code is imported from the
paper's repository (github.com/pandey-tushar/QLDPC, not consulted). The only
shared dependency is the decoder library (ldpc/bposd core) — same as theirs.

## Gate B — ν drift test at larger N (decision experiment)

**Question:** does ν drift toward 4/3 (finite-size artifact) or stay ≈ 1.18
(decoder-dependent universality) when sizes beyond N = 1296 enter the fit?
**Design:** add (42,21) N=1764 (K=8, rate 0.45%) and (48,24) N=2304 (K=16,
rate 0.69%) — K from exact GF(2) rank at build (a "Ki=12" guess in an
earlier draft was wrong; see campaigns-smoke/feasibility_bench.json) —
to the five original sizes; re-fit with the same
ansatz; report the ν shift with bootstrap CIs. Falsification criterion:
|Δν| ≥ 0.1 with non-overlapping CIs = drift detected (artifact-leaning);
stable ν within CIs = decoder-dependent-universality-leaning; the p*_inf
shift is reported either way.
**Feasibility analysis (one core, bounded stop conditions):** measured
per-shot decode wall-time from `scripts/bench_feasibility.py`
(campaigns-smoke/feasibility_bench.json), single core (OMP/BLAS pinned),
nice -n 10, at the near-threshold worst case p = 0.45:

| (L,M) | N | ms/shot | 200k eval | search ~6 evals |
|---|---|---|---|---|
| 12x6 | 144 | 0.69 | 0.04 h | 0.23 h |
| 18x9 | 324 | 3.06 | 0.17 h | 1.02 h |
| 24x12 | 576 | 11.3 | 0.63 h | 3.8 h |
| 30x15 | 900 | 38.7 | 2.15 h | 12.9 h |
| 36x18 | 1296 | 97.6 | 5.4 h | 32.5 h |
| 42x21 | 1764 | 221 | 12.3 h | 73.8 h |
| 48x24 | 2304 | 441 | 24.5 h | 147.0 h |

Consequences (decision, 2026-08-29):

- **Gate A** at the paper floor (200k shots, 5 sizes, ~6 evals each):
  ≈ 50 h single-core nice-10. Feasible only as a scheduled multi-session
  run; NOT launched today. A time-scaled alternative remains open: 20k
  shots/eval costs ≈ 5 h total and brackets p*(N) CIs ≈ 3.2× wider
  (1/sqrt(10) shot-noise scaling) — still able to test overlap with the
  paper's Table-1 CIs only marginally; the decision to run either variant
  is recorded as PENDING in README.md.
- **Gate B** with both new sizes at the paper floor: 73.8 + 147.0 = 221 h ≫
  the pre-committed 48 h ceiling. Gate B is therefore PARKED at the paper
  floor. Re-scoping options (require a pre-statement amendment BEFORE any
  fit uses the data): (a) drop N=2304, run N=1764 only (≈ 74 h — still over
  ceiling), (b) reduced-shot Gate B (20k/eval: N=1764 ≈ 7.4 h, N=2304 ≈
  14.7 h — inside ceiling together, CIs ~3.2× wider), (c) mixed-channel
  points where WER curves are steeper (fewer evals to bracket). Option (b)
  is the recommended default if Gate B is approved.
- Per-shot decode cost grows ≈ N^2 (single-core BP+OSD), so any campaign
  re-scope should prefer adding windows/shots at existing sizes over new
  sizes when the target is p*_inf precision vs ν resolution.

**Stop conditions (both gates):** hard data cap per campaign; seeds fixed
before launch; no adaptive shot top-ups; window set fixed in advance
(0.04/0.06/0.08) — window selection from data is forbidden.

## Module validation before real sweeps (DONE, artifact)

`campaigns-smoke/fss_module_validation.json` [REPRODUCED from synthetic
truth, not from the paper]: the collapse estimator recovers known synthetic
exponents (p0, ν0) = (0.55, 1.35) → (0.5499 [0.5499, 0.5500], 1.349
[1.347, 1.351]) inside 95% bootstrap CIs; out-of-class truth (logistic f,
(0.61, 1.05)) recovers (0.6100, 1.050) — sub-1% model bias. A previous
version of the estimator failed these tests (single-init local-minimum
degeneracy); the multi-start grid init was added to fix it — that history is
visible in this file's git record.

## Status

- [x] pre-statement committed before first real run
- [x] src/ pipeline built and smoke-tested (N=144, 2000-shot campaign
      smoke: p* = 0.3720 CI95 (0.3715, 0.3719) vs paper 0.3701 ± 0.00034 at
      200k shots — consistent within shot-noise)
- [x] **Gate-A fragment at the paper floor [REPRODUCED]**
      (campaigns/20260829T224546Z_987904a9_40c474752126, 2026-08-29):
      N=144, 200k shots, exact paper estimator chain — p* = 0.37032,
      95% CI (0.37004, 0.37101) vs their 0.3701 [0.3697, 0.3703]:
      CIs overlap. This is the first, smallest rung of Gate A only — NOT
      the full gate (p*_inf/ν need all five sizes).
- [x] FSS module validated on synthetic data with known exponents
- [x] Gate-B feasibility benchmark (campaigns-smoke/feasibility_bench.json)
- [ ] Gate A full campaign (5 sizes × 200k shots ≈ 50 h core) — PENDING
      budget decision (20k-shot fallback ≈ 5 h, ~3.2x wider CIs)
- [ ] Gate B (PARKED; blocked on Gate A outcome + feasibility re-scope
      amendment)
