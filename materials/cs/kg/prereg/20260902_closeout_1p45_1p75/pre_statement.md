# Pre-statement — corrected-envelope lineage close-out of [1.45, 1.75],
# primary cap 2^17, registered 2^19 fallback; separate successor campaign

Agent `KgBand175`, UTC 2026-09-02. **COMPUTE PENDING.** This file is committed
before the fresh `campaign.py init` and before any compute of this campaign.
`campaign.py init` will mint the run directory; a byte-identical copy of this
file plus its source commit + sha256 is recorded there before controls.

## 0. Why this campaign exists, and what it supersedes

This is a **lineage close-out, not a frontier gap**. The sound direct-D.4
instrument already holds a frozen certificate for this band: 10/10 old
ratio-1.02 tiles, worst margin +3.10559377267758493230e-6. That PASS stands
independently and is not superseded by anything here. This campaign
re-earns the band under the corrected envelope lineage (Campaign A's frozen
`gate_corr_env.py`, pointwise looser-or-equal upper bounds), per owner
instruction of 2026-09-02, under its own prereg/init/controls.

Precondition (input, not authority): the primary band campaign
`20260901T135742Z_1a48b071_282f6f2490bc` closed **FROZEN-CERTIFIED**
(87/87 tiles PASS on [1.75,3.5]). Two frozen artifacts of that campaign are
cited here: (a) its driver-written `logs/band_1p75_3p5_result.json`, sha256
`f42d0463f11f4d7014febe9c6f799b8754c73d436a7ba227173265ae270e8bf1` as
committed at `ff14541`; (b) its run-root summary `result.json`, sha256
`93588cc2a7ff510c4800cb66fe911b854db30db5acbd2452dd0a21f1ea3b5b22`. Item
(a) is copied byte-identical into this run's `logs/` before band compute
solely to satisfy the driver's registered close-out precondition guard. Its
content is never used as arithmetic input. This band shares the endpoint
band's left edge; the tiling endpoints are exact on both sides and no
double-counting occurs.

## 1. Frozen lineage and byte identity

Identical to the primary campaign's §1: Campaign A's frozen
`gate_corr_env.py` sha256
`94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30`;
Campaign A/C byte-identical `startup_controls.py` sha256
`cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3`;
Campaign B frozen `decompose_cell.py` sha256
`291eddc1d76b5a60338855012629077a7d5799db9c7434e07b18ae19df6204ff` and its
dependency `gate_subdiv.py` sha256
`6f1b8ee0f4bd11e3a3dec1c0a584dfe72975c9b91ab53ff2b26366cb2e0b5df5` (D1–D6
escalation template only); shared dependency
`campaigns/20260831T082425Z_kg_direct_d4_restart/code/core.py`. Run copies
are byte-identical, shas asserted in-run before any import or arithmetic.
`PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode = True` before any
frozen import; no bytecode inside any frozen directory or this run.

The only knobs set are the two inherited from A: panel cap 131072 = 2^17 for
the primary ladder, explicit fallback rung 524288 = 2^19 on an otherwise
terminal leaf, and the named c-range [1.45, 1.75]. No envelope is re-derived;
root cover, box floor 1/1024, depth max 40, S_MAX 8, closed tail, 256-bit
outward Arb (320 only inside the startup battery), panel weights, margin
predicate `d(c_U).lower() - U(B,c_L).upper() > 0`, DFS order, child reversal,
split rule, checkpointing, `--start-tile` resume, and stop-at-first-failing
-tile are exactly as in the primary campaign's frozen driver
(`20260901T135742Z_1a48b071_282f6f2490bc/code/gate_band175.py`, computational
core byte-identical; only the header docstring names this campaign).

## 2. Frozen band parameters

- Band: exact decimals `c ∈ [1.45, 1.75]`.
- Tiling: A's `c_tiles('1.45','1.75')`, `TILE_RATIO='1.008'`, final endpoint
  clamped exactly to 1.75. Assert `len(tiles)==24` and
  `assert_tile_cover('1.45','1.75',tiles)` before tile compute.
- Ladder: A's inherited `panel_ladder(c_left)`, asserted non-empty per tile;
  rungs are whatever A's frozen module returns for this band — no adaptation.
- Everything else (root cover 132, box floor, depth, S_MAX, precision,
  margin predicate, 260,000-evaluation budget, 60 CPU-hour campaign budget
  from `time.process_time` deltas, 72-hour wall ceiling from manifest,
  per-tile atomic freeze with exact frontier serialization, refutation gate
  with escalation to Main) is the primary campaign's §2 verbatim.
- Runtime disclosure: compute runs inside the worker agent's foreground
  eval kernel at niceness 0, `RLIMIT_CPU` soft raised to unlimited, all
  thread pools pinned to 1. Priority is scheduling-only.

Measured-rate note (orientation only, from the primary campaign): midband
tiles cost 1000-2100 s CPU; this band is expected to be cheaper, but no
verdict depends on that expectation.

## 3. Required controls after init, before band compute

The primary campaign's §3 battery, verbatim and in order: startup battery
(16/16, output byte-identical, sha256
`211d8e9afd6491af62d7bb89282e586f300f5a32535320caa31aa342fc63beac`);
instrument-identity reconciliation at (4.083, 4.115664) against frozen
+1.04499016763434518123845756638e-6 with 2^-40 ball tolerance; near-miss
discriminator on the OLD ratio-1.02 midband cell, padded 1e-20, 32768
negative then 131072 positive, D1–D6 escalation on failure; ratio-1.05
planted-failure control strictly negative at 131072. Each JSON frozen before
the next starts.

## 4. DFS and verdict rule

Tiles 1..24 in order; per-cell primary ladder; 2^19 fallback exactly on
unsplittable terminal leaves; per-tile freeze and progress rewrite as in the
primary. Verdict rule: PASS only if all 24 tiles finish strictly positive;
FROZEN-NEGATIVE at the first terminal OPEN leaf after the registered 2^19
rung or a control failure; FROZEN-INCONCLUSIVE if the budget binds first.
A PASS here supersedes nothing; an OPEN/INCONCLUSIVE here retracts nothing
from the standing direct-D.4 certificate.

## 5. Freeze and report

Before close: `result.json` at run root with the certified prefix, per-tile
table, aggregate histogram, CPU, frontier if any, refutation-gate status, and
exact scope not swept; result section appended to the run copy of this
pre-statement with both hashes recorded in `provenance.json`. Then
`campaign.py freeze` and exactly one `campaign.py close --verdict ...`.
