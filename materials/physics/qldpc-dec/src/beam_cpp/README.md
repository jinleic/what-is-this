# beam8_cpp — C++ port of the Python beam8 beam-search decoder

C++ port of `qldpc_dec.beam_search.BeamSearchDecoder` (the "beam8" arm of
Gate B; arXiv:2512.07057 Algorithm 3) at the pinned companion parameters
(beam_width=8, initial_iters=30, iters_per_round=20, max_rounds=10,
num_results=1). Decodes the same raw-convention DEM for a batch of shots and
writes per-shot predicted-observable bits; `equiv.py` diffs C++ vs Python
per-shot outcomes.

## What runs where

- `beam8.cpp` / binary `beam8_cpp` — the port. Inputs: stim textual DEM
  (straight-line `error(p) D.. L..` lines only — the canonical form of every
  committed BB_144 DEM in this repo), 1-D float64 `.npy` of `lam` values, and
  a bit-packed shots file. Output: little-endian binary of per-shot
  observable words. No RNG anywhere (the Python decoder is deterministic;
  the only seeded component, Stim's sampler, runs on the Python side).
- `equiv.py` — equivalence + timing harness. Generates shared seeded shots,
  runs both arms, diffs per-shot predictions, writes a JSON report.

## Equivalence contract (what had to match to the bit)

1. `lam_j = log((1-p)/clip(p,1e-12))` — np.log can differ from std::log by
   1 ulp, so the harness exports the exact float64 array the Python decoder
   uses (`--lam`); `--check-lam` reports any bit difference when computing
   in-binary.
2. numpy float64 pairwise-summation order (PW_BLOCKSIZE=128, unrolled-8
   accumulator tree) for every `np.sum`/`np.mean` in the pipeline, including
   the fancy-indexed `mu[ks].sum()` (gather-then-sum, never a raw slice —
   an error's edge pairs are scattered because pairs are grouped by
   detector).
3. Tie-break orders: `np.argsort(kind="stable")` min-rank check in A1;
   `np.argmin` first-minimum scans; stable `next_set.sort(key=score,
   reverse=True)` (equal scores keep creation order); `(path, val)` expansion
   order; first converged+valid child returns immediately (num_results=1).
4. Exact early-return / fallback structure of `decode()`, including
   converged-but-invalid decomposition results continuing into the beam.

## Build

    clang++ -O3 -std=c++17 -Wall -Wextra beam8.cpp -o beam8_cpp

(no external dependencies; any C++17 compiler works, g++ fine too)

## Run

    ./beam8_cpp --dem <dem.dem> --lam <lam.npy> --shots <shots.bin> \
                --out <preds.bin> [--limit N] [--repeats K] [--threads T] \
                [--beam-width W] [--initial-iters I] [--iters-per-round P] \
                [--max-rounds R] [--check-lam] [--verbose]

- shots file: 24-byte little-endian header `<ndet:u64> <nobs:u64> <n:u64>`,
  then `n` rows of detector bits, then `n` rows of observable bits; each row
  is ceil(count/64) little-endian u64 words, bit b of word w = element 64w+b
  (Stim's own packing; produced by `equiv.py`).
- output file: 12-byte header `<magic:u32=0x42385031 ("B8P1")>
  <nobs:u32> <n:u32>` then per-shot u64-LE observable words.
- stderr reports wall-ms and the aggregate convergence counters.

## Equivalence harness

Exact commands (from this dir; /Users/jinleic/jinleic-workspace/physics/.venv/bin/python is the venv python 3.14 with numpy/scipy/stim):

    PYTHONPATH=../../src /Users/jinleic/jinleic-workspace/physics/.venv/bin/python \
        equiv.py --p 3e-3 --shots 1000 --seed 20260830 --cpp-repeats 3

- `--p` picks the committed circuit: 1e-3 uses the committed IonQ circuit,
  3e-3 uses the committed DERIVED rescale artifact (the same files the frozen
  Gate-B companion ran on).
- one-off arms: re-run C++ only with `--cpp-only` (reuses cached Python predictions).

## Results

See `scratch/equiv_p3e-3_n1000_s20260830.json` (and the p=1e-3 twin) for the
full per-run record: per-shot mismatches (count + first ids), python/C++
ms/shot, and speedup. C++ timing excludes DEM parse; Python timing is the
batch wall-clock (`decode_batch`), matching how the frozen 2362 ms/shot
companion number was produced.
