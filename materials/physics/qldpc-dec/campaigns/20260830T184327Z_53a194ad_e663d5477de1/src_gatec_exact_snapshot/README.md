# `gatec_exact/` — exact degenerate-ML reference at small scale (Gate C companion)

**Status: COMPLETE (2026-08-30).** Exact per-shot degenerate-ML reference for the
**[[36,4,4]]** bivariate-bicycle code — the smallest instance of the Bravyi et al.
monomial family used by this harness — cross-checked against the three harness
decoders (beam8, bp30+osd, NMS-ensemble-24). Scratch outputs:
`qldpc-dec/scratch/gatec/`.

## Scope (read before quoting any number here)

- This is the **small-scale complement** the task charter calls for. It is **NOT**
  the [[144,12,12]] headline configuration: no number in this module says anything
  about the [[144,12,12]] gates, and none can. [[144,12,12]] circuit-level exact ML
  (936 detectors, ~8.8k mechanisms) is out of reach of this machinery by a factor
  of ~2^(detectors) in the affine-space dimension.
- arXiv:2608.25545 gives *certified* exact ML for [[72,12,6]] via region-based Bethe
  with elimination clusters; our scope is a cheaper, fully-controlled independent
  reference on a smaller instance of the same code family (same 3-term monomial
  pair, smaller torus).

## Code instance

- Grid (L, M) = (6, 3), monomials A = x^3 + y + y^2, B = y + x + x^3
  (same family as the [[72,12,6]] / [[144,12,12]] codes; "much smaller a, b,
  low distance" per the task). Parameters verified by exhaustive span
  enumeration: **[[36, 4, 4]]** — n = 36 data qubits (≤ 40 target), k = 4
  logicals (≤ 12 target), d_X = d_Z = 4.
- Z-memory sector only; the X-memory sector is the transpose (Lx vs Lz swap)
  and was not run separately.

## Noise model: code-capacity depolarizing, with exact marginalization

Committed-circuit noise is circuit-level depolarizing; exact ML over a
circuit-level DEM of even this small code is not tractable by the coset method
(DEM mechanisms ≈ 1.7e4 → affine dimension enormous; see "Cost wall" below).
The reference therefore conditions on the **code-capacity** limit, which is the
sector where the coset machinery is EXACT and complete:

- per data qubit: depolarizing channel with parameter p;
- Z-memory syndrome/observables depend only on the X-error component;
- marginalizing each qubit's channel over the Z-component leaves an
  independent X-error per qubit with probability **q = 2p/3**;
- stim DEM = 36 single-qubit X mechanisms, prior q each (verified numerically:
  `dem_to_matrices_raw` priors == q).

Points run: p in {1e-2, 3e-2} (circuit-equivalent), i.e. q in {6.67e-3, 2e-2}
— chosen so that all decoders make errors "sometimes" (task charter) rather
than the too-clean 1e-3 regime.

## Exact method (complete, no truncation)

For a sampled syndrome s (18 bits), the consistent error set
{e : H e = s} is an affine space of dimension 36 − 16 = **20**
(2^20 = 1,048,576 members) — enumerated COMPLETELY per shot as packed 36-bit
words:

1. nullspace span table E0 (build once): 2^20 uint64 words;
2. per shot: particular solution via joint [H | s] GF(2) elimination
   (the correct rhs transform; a raw-s back-substitution is a known-wrong
   shortcut rejected during validation), Es = E0 XOR e_p;
3. integer weight histograms n[c, w] over all 2^20 solutions split by the
   4-bit logical class (16 classes), with a hard assertion that the histogram
   totals 2^20 every shot (partition of the affine space);
4. exact rational class masses Z_c = sum_w n[c,w] r^w, r = q/(1-q) = 2p/(3-2p)
   via Python `Fraction` — the ML argmax uses **no floating point**;
5. ML class c* = argmax_c Z_c; a decoder that picks class c gets
   regret = ln(Z_{c*}/Z_c) >= 0 in nats (0 iff it matches ML).

Degeneracy is handled by construction: every stabilizer-equivalent solution
contributes its prior weight, so this is genuine degenerate ML (coset-sum
partition function), not minimum-weight decoding.

### Cost wall (documented honestly)

The same machinery on the **circuit-level** DEM (12-round memory, built and
verified during development: 432 detectors, 16,732 mechanisms, decomposed
form) would need an affine space of dimension ~1.6e4 — exhausting 2^(n-rank)
is impossible; the interaction graph of that DEM has substantial treewidth
(data errors persist across all 12 rounds, coupling every round's detectors),
so transfer-matrix elimination does not factor it either. The code-capacity
limit above is the largest instance where a COMPLETE exact class sum is
affordable. Per-shot cost here: ~7-10 ms single core (M3 Ultra, numpy 2.5);
2000 shots ≈ 20 s per noise point for the exact reference.

## Validation chain (all green, run in-source)

1. **Tiny-code exhaustive check**: 3-qubit repetition code, ALL 2^3 errors:
   class weights of the histogram machinery equal direct enumeration exactly.
2. **Brute-force low-weight check** (`validate_exact_vs_bruteforce`): for
   random syndromes of the real [[36,4,4]] instance (syndromes taken from
   random full-weight errors, not planted low-weight ones), the exact
   histogram rows for w <= 5 equal brute-force enumeration of all weight-<=5
   errors. Runs automatically at the start of every `run_gatec` point.
3. **Partition assertion**: every shot's histogram must sum to exactly 2^20
   (built into `class_histogram`).
4. **Code-parameter assertion**: an independent distance computation
   (`bb_distances` via exhaustive ker/rowspace coset enumeration) must return
   [[36,4,4]] at build time, or the runner aborts.
5. Syndrome consistency: the true error's class always has nonzero mass.

## Results (2026-08-30, 2000 shots/point, base seed 20260830)

See `scratch/gatec/gatec_p{0.01,0.03}_2000shots.json` — the summary tables
printed by `run_gatec.py` carry the headline numbers (per-decoder: logical
failures, ML-mismatch rate, mean/max regret in nats, ms/shot; plus
exact-ML's own logical failure count — ML minimizes P(wrong class); it is
not perfect).

## Files

- `exact_ml.py` — code construction (BB parity checks, exact distances,
  logicals), code-capacity circuit, DEM extraction, `ExactMLReference`
  (coset machinery + validators).
- `run_gatec.py` — CLI runner: builds instance, validates, samples, runs 3
  decoder arms, exact ML, regrets; writes JSON report + prints table.
- `benchmark.sh` — the full reproduction run (both p points, 2000 shots).
- scratch outputs: `../../scratch/gatec/` (JSON reports).

## Interface notes

- Decoders are the harness implementations, imported unmodified from
  `qldpc_dec.beam_search` / `qldpc_dec.bp_osd` / `qldpc_dec.nms`, bound to the
  same stim DEM (beam/NMS raw convention, BP+OSD merge convention — each
  decoder's native convention, same underlying DEM).
- Sampling seeds derived via `qldpc_dec.seeds.derive_seed` (repo seed policy).
- No campaign snapshots were created (freeze privileges are Main's); JSON
  reports land in `qldpc-dec/scratch/gatec/` per the charter.

## Honest limitations

- One sector (Z-memory) run; X-sector is the transpose and omitted.
- Code-capacity, not circuit-level: the optimality gaps measured here bound
  what an exact circuit-level reference would show only weakly; they are a
  clean small instance where "how far from optimal is each decoder" is
  answerable exactly, not a reproduction of any published circuit-level number.
- The [[72,12,6]] instance of the same family was scanned too; n=72 exceeds
  the <=40-data-qubit target and its affine dimension is 2^40 — out of reach
  by a factor of ~1e6 vs the [[36,4,4]] space.
