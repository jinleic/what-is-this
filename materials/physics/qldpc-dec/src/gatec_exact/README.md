# `gatec_exact/` — exact degenerate-ML reference at small scale (Gate C companion)

**Status: COMPLETE in both memory sectors with GB3-corrected BP+OSD
comparators (2026-08-31).** Exact per-shot degenerate-ML reference for the
**[[36,4,4]]** bivariate-bicycle code — the smallest instance of the Bravyi
et al. monomial family used by this harness — cross-checked against beam8,
BP+OSD, and NMS-ensemble-24. Authoritative corrected evidence: Z memory
`campaigns/20260831T093503Z_2f2ed071_884b7d216c6f/`; X memory
`campaigns/20260831T093513Z_3039c8d0_4a350b3b0b06/`.

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
- Both CSS sectors are now run independently: Z memory uses X errors with
  `(Hz, Lz)`; X memory uses the Hadamard-dual Z errors with `(Hx, Lx)`.

## Noise model: code-capacity depolarizing, with exact marginalization

Committed-circuit noise is circuit-level depolarizing; exact ML over a
circuit-level DEM of even this small code is not tractable by the coset method
(DEM mechanisms ≈ 1.7e4 → affine dimension enormous; see "Cost wall" below).
The reference therefore conditions on the **code-capacity** limit, which is the
sector where the coset machinery is EXACT and complete:

- per data qubit: depolarizing channel with parameter p;
- Z-memory syndrome/observables depend on the X-error component; X memory
  depends on the Z-error component;
- either marginalization leaves an independent opposite-Pauli error per qubit
  with probability **q = 2p/3**;
- each stim DEM has 36 single-qubit mechanisms at prior q (verified
  numerically by `dem_matrices_raw`).

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

## Results (2000 shots/point/sector, base seed 20260830)

All values are code-capacity block-failure rates, not 12-round circuit LERs:

| p | basis | exact ML | beam8 | bp30+osd | nms-ens24 |
|---:|:---:|---:|---:|---:|---:|
| 0.03 | Z | **0.0145** | 0.0325 | 0.0330 | 0.0260 |
| 0.03 | X | **0.0145** | 0.0370 | 0.0365 | 0.0230 |
| 0.01 | Z | **0.0005** | 0.0040 | 0.0055 | 0.0020 |
| 0.01 | X | **0.0005** | 0.0045 | 0.0025 | 0.0010 |

Exact ML beats every heuristic arm at both points in both sectors. GB3
reruns assert that the configured BP channel equals the complete merged-DEM
prior vector. All 36 code-capacity mechanisms have the same prior q, so the
former scalar mean happened to be exactly equivalent on this small-code
experiment. All four corrected per-shot CSVs are byte-identical to their
historical counterparts; this equivalence does not extend to the
heterogeneous 8784-column circuit-level DEM.

## Files

- `exact_ml.py` — BB construction, exact distances/logicals, basis-explicit
  code-capacity circuits, DEM extraction, and `ExactMLReference`.
- `run_gatec.py` — `--basis X|Z` point runner: validates, samples, runs the
  three heuristic arms, asserts the BP channel, and computes exact ML/regrets.
- `run_campaign.py` — immutable two-point campaign wrapper with source
  snapshots, circuit/DEM hashes, compressed results, and close inventory.
- `benchmark.sh` — both p points; set `BASIS=X` for the X-memory sector.
- Corrected Z campaign: `../../campaigns/20260831T093503Z_2f2ed071_884b7d216c6f/`.
- Corrected X campaign: `../../campaigns/20260831T093513Z_3039c8d0_4a350b3b0b06/`.

## Interface notes

- Decoders are the harness implementations, imported from
  `qldpc_dec.beam_search` / `qldpc_dec.bp_osd` / `qldpc_dec.nms`, bound to the
  same stim DEM. Corrected campaigns use and assert the merged DEM's complete
  BP error channel; exact/beam/NMS outputs also regression-match the prior
  campaigns byte-for-byte.
- Sampling seeds derived via `qldpc_dec.seeds.derive_seed` (repo seed policy).
- Campaign snapshots contain pre-run manifests, source hashes and snapshots,
  both aggregate reports, per-shot CSVs, and immutable-close inventories.

## Honest limitations

- Both CSS memory sectors are covered, but only under code-capacity noise.
  The measured optimality gaps constrain a clean small instance; they only
  weakly bound a circuit-level reference and reproduce no published
  circuit-level number.
- The [[72,12,6]] instance of the same family was scanned too; n=72 exceeds
  the <=40-data-qubit target and its affine dimension is 2^40 — out of reach
  by a factor of ~1e6 vs the [[36,4,4]] space.
