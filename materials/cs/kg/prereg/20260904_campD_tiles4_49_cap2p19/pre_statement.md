# Pre-statement — Campaign D: close [4.083, 6.0] tiles 4–49 at cap 2^19

Agent `KgBand175`, UTC 2026-09-04. **COMPUTE PENDING.** This file is committed before fresh `campaign.py init`, controls, or any D-tile arithmetic. A byte-identical copy plus source commit and sha256 is recorded in the minted run.

## 0. Frozen lineage prefix and exact scope

This campaign certifies exactly tiles 4–49 of Campaign A/C's unchanged 49-tile table `c_tiles('4.083','6.0')` (`TILE_RATIO='1.008'`). Its left edge is Campaign C tile 3's exact right edge `c=4.181778026496`; tile 49 ends exactly at 6.0. Before D compute: assert 49 total tiles, exact cover, adjacency of every pair, tile-4 left edge equals tile-3 right edge, and final endpoint equals 6.0.

The prefix is frozen input, not recomputed and not charged to D's budget:

- A tile 1 PASS at cap 2^17: `20260831T231500Z_KgEvenBoxRerun_campA_corr_env/logs/inproc_tile_001.json`, sha256 `08f5e3bfaafefcc95006568467b2bc1ad6e7876a67f00e27a07fa793cd2326a5`;
- A tile 2 PASS at cap 2^17: sibling `inproc_tile_002.json`, sha256 `623598b887db1cfc0cea1bcd64aa0375b7d72417f09eae431c7f7c9e5b9441e6`;
- C tile 3 PASS at cap 2^19: `20260901T055808Z_KgCampC_pan2p19_campC/logs/inproc_tile_003.json`, sha256 `8b02e4eb24650edb715ad696c90e29a5f104f176a7d33dabeec06f1452af1c99` (484/484 leaves, worst margin `[5.06389732008514710540159866096e-7 +/- 2.76e-37]`).

A D PASS combines those standing prefix certificates with fresh tiles 4–49 and closes the whole [4.083,6.0] band. A D OPEN/budget frontier claims only the exact contiguous prefix through the last completed D tile and retracts none of A/C.

## 1. Frozen instrument; one knob

Copy Campaign A's frozen `gate_corr_env.py` byte-identically, sha256 `94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30`, and Campaign A/C's byte-identical `startup_controls.py`, sha256 `cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3`. Shared dependency remains `20260831T082425Z_kg_direct_d4_restart/code/core.py`. Assert hashes before import. Set `PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode=True` before frozen imports; no bytecode in frozen/live runs.

The only computational knob changed from A is `gce.PANEL_CAP: 131072 -> 524288 = 2^19`, immediately after import and before any envelope arithmetic. Then assert the high-band ladder equals `[4096,16384,65536,262144,524288]`. C used the same single knob. No envelope, tiling, root cover, box floor (`1/1024`), depth (`40`), `S_MAX=8`, tail, panel weights, margin predicate, DFS/stack/child order, disk filter, or precision is changed. All verdict arithmetic is 256-bit outward Arb.

## 2. Controls before tile 4

Controls run in order and freeze after each:

1. Startup battery: all 16 families PASS, `aborted:null`, and output byte-identical to A/B/C's `startup_controls.json`, expected sha256 `211d8e9afd6491af62d7bb89282e586f300f5a32535320caa31aa342fc63beac`.
2. Ladder/cover gate: cap exactly 524288; ladder exactly `[4096,16384,65536,262144,524288]`; 49-tile exact cover and prefix shas/adjacency all green.
3. Cap-2^17 reconciliation (non-verdict gate): on C's exact tile-3 frontier cell `a0 in [-638/768,-637/768]`, `a2 in [429/768,430/768]`, c-pair `(4.148589312,4.181778026496)`, direct n=131072 margin reproduces `-1.76291787028287693244667476546e-5` within ball radius `<2^-40`.
4. Known-true accept: same cell/c-pair, direct n=524288 margin lower endpoint >0 and agrees with C's frozen `+4.46309091764170239642742026145e-6` within `2^-40`. C gates source `phase1_gates.json`, sha256 `ae7e8c530e8447e3b030d8445f5b239f10db11a7b7d6a01a2832cf6dc03022ff`.
5. Planted reject: same frontier cell and `cL=4.148589312`, deliberately widen `cR_probe=1.05*cL`, direct n=524288; required margin upper endpoint <0. This checks that even restored cap 2^19 rejects an over-wide c-tile.

Any startup, identity, ladder/cover, known-accept, reconciliation, or planted-reject control failure stops before tile 4 and classifies the instrument as invalid; it closes **FROZEN-INCONCLUSIVE**, never FROZEN-NEGATIVE, with the exact failed ball. No parameter adapts to control outcomes.

## 3. D-tile DFS and checkpoints

Run tiles 4,5,...,49 in order. For each tile, call A's `certify_tile` unchanged after the single cap patch. Require all disk-intersecting leaves strictly positive under `d(c_U).lower()-U(c_L).upper()>0`; straddling is OPEN. Stop at first failing tile. Do not run a later tile after OPEN.

After every completed tile atomically freeze `logs/d_tile_NNN.json` with exact c-pair, root/box/leaf counts, depth/open cell, panel count, margin ball, remaining stack count, panel histogram, envelope evaluations, and `time.process_time` CPU. Atomically rewrite `logs/d_tiles4_49_progress.json` with all completed D rows, aggregate histogram/CPU, exact certified prefix, and `next_tile`. Resume only at the exact next tile; completed tiles are never recomputed. External interruption mid-tile produces no verdict artifact and resumes that tile.

Budget is checked at each tile boundary; an already-started exact tile completes atomically. Registered D budget: **200 CPU-hours = 720000 s** for fresh tiles 4–49 only; prefix/control CPU reported separately. Registered wall ceiling: **240 h = 864000 s** from manifest creation. If the boundary budget binds, freeze the exact boundary frontier (next tile/c-pair, completed prefix, remaining named tiles, histograms, CPU/wall used) and close FROZEN-INCONCLUSIVE. This exceeds Campaign C's registered-next-action minimum `>=160 CPU-h` and its like-tile-3 projection 157.6 CPU-h.

Runtime: sole writer `KgBand175`; foreground eval kernel; niceness 0; `RLIMIT_CPU` soft unlimited; OMP/OpenBLAS/MKL/vecLib/NumExpr threads pinned to 1. No competing writer.

## 4. Verdict and freeze

- `FROZEN-CERTIFIED`: all D tiles 4–49 PASS; combined standing A tiles 1–2 + C tile 3 + D tiles 4–49 certify the full [4.083,6.0] band.
- `FROZEN-NEGATIVE`: first terminal OPEN leaf under cap 2^19 after every instrument/control gate has passed; freeze exact obstruction and certified prefix only. This verdict is reserved for a genuine sound-instrument OPEN leaf.
- `FROZEN-INCONCLUSIVE`: any instrument/control failure, a registered CPU/wall boundary binding first, or an external interruption preventing a terminal result; freeze the failed control or exact prefix/frontier.

Before close: run-root `result.json`, per-tile table/aggregate histogram/CPU, prefix provenance, refutation-gate status, exact unswept scope, and a result appendix in this run copy with both pre/post hashes in `provenance.json`. Then `campaign.py freeze` and exactly one `campaign.py close --verdict ...`.
