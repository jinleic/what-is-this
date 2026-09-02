# Pre-statement — corrected-envelope certification of [1.75, 3.5]
# primary cap 2^17, registered 2^19 fallback; conditional lineage close-outs

Agent `KgBand175`, UTC 2026-09-01. **COMPUTE PENDING.** This file is
committed before the fresh `campaign.py init` and before any compute of the
fresh campaign. `campaign.py init` will mint its run directory; a
byte-identical copy of this file and its source commit + sha256 will then be
recorded there before controls.

## 0. Protocol provenance and excluded scouting

A first no-compute run, `20260901T134705Z_9da59c74_a3bb7abdb738`, was frozen
and closed **REJECTED** before controls because the worker had run two
read-only scouting checks before that run's prereg commit. Its retained
`protocol_rejection.json`, `sha256s.txt`, and terminal `status.json` carry no
certificate weight.

Those excluded pre-commit checks were:

1. the byte-shared `c_tiles('1.75','3.5')` count/cover, which returned 87
   tiles, first pair (1.75, 1.764), last pair
   (3.4725124389279784651..., 3.5);
2. the owner's named near-miss SUPERSET cell at the frozen old-tiling c-pair
   (3.16988277218086405264471296000,
   3.23328042762448133369760721920), which returned
   `[-4.14781740358947041431994319933e-8 +/- 3.81e-38]` at 32768 panels and
   `[8.72566255993993200750351436424e-5 +/- 4.25e-35]` at 131072 panels.

They are disclosed as **EXCLUDED SCOUTING ONLY** and are not evidence of
this fresh campaign. No threshold or knob could adapt to them: before either
check, the assignment had already fixed the c-range, ratio-1.008 tiling,
primary cap 2^17, fallback cap 2^19, exact near-miss cell, expected
negative→positive flip, ratio-1.05 planted-failure direction, root cover,
box floor, depth, S_MAX, precision, margin predicate, DFS order, budget, and
verdict rule. The fresh campaign reruns every required control after init and
lets its own DFS earn every tile verdict.

Owner/Main scouting supplied independently — including one PASSing real
`certify_tile` at tile 75/87 (1781.8 s CPU, cap 2^17) — is likewise cost
orientation only, never a certificate or authority.

## 1. Frozen lineage and byte identity

This campaign consumes:

- Campaign A:
  `campaigns/20260831T231500Z_KgEvenBoxRerun_campA_corr_env/code/gate_corr_env.py`,
  sha256
  `94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30`;
- Campaign A/C's byte-identical `startup_controls.py`, sha256
  `cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3`;
- Campaign B's frozen `decompose_cell.py`, sha256
  `291eddc1d76b5a60338855012629077a7d5799db9c7434e07b18ae19df6204ff`,
  as the D1–D6 template for the Gate-N escalation branch;
- shared dependency
  `campaigns/20260831T082425Z_kg_direct_d4_restart/code/core.py`, imported
  exactly as A/B/C imported it.

The run copies are byte-identical and their shas are asserted in-run before
any module import or arithmetic. `PYTHONDONTWRITEBYTECODE=1` and
`sys.dont_write_bytecode = True` are set before importing any frozen campaign
module. No bytecode may be written into any frozen directory or the live run.

The only knobs set are:

1. panel cap 131072 = 2^17 (A's inherited default) for the primary ladder;
2. explicit panel count 524288 = 2^19 as the registered fallback rung on an
   otherwise terminal leaf that remains open at 2^17;
3. the named c-ranges [1.75,3.5], [1.0,1.3], [1.45,1.75].

No envelope is re-derived. `root_boxes`, `e_panel_range`, `cell_envelope`,
`dist0_box`, `box_intersects_disk`, `split_box`, `can_split`,
`strict_margin_lower`, panel weights, tail, and every Arb operation are called
from A's byte-identical module. The driver mirrors A's DFS loop only to add
checkpoints, budget checks, full frontier serialization, and the one
registered 2^19 fallback at the exact point where A would otherwise return a
terminal OPEN leaf. Root order, stack order, child reversal, split rule, and
stop-at-first-failing-tile are unchanged.

## 2. Frozen primary-band parameters

- Band: exact decimals `c ∈ [1.75,3.5]`.
- Tiling: A's `c_tiles('1.75','3.5')`, `TILE_RATIO='1.008'`, final endpoint
  clamped exactly to 3.5. Before tile compute assert `len(tiles)==87` and
  `assert_tile_cover('1.75','3.5',tiles)`.
- Root cover: `root_boxes(12)`, exactly 132 disk-intersecting roots asserted.
- Box floor: `BOX_FLOOR_Q = 1/1024`, unchanged.
- Maximum depth: `BOX_DEPTH_MAX = 40`, unchanged.
- `S_MAX = 8`, closed tail unchanged.
- Precision: 256-bit outward Arb; 320 bits only inside the inherited startup
  battery.
- Primary ladder on this band: `[1024,4096,16384,65536,131072]`, asserted
  from A's `panel_ladder` before each tile.
- Margin predicate: exactly
  `d(c_U).lower() - U(B,c_L).upper() > 0`; strictly positive lower endpoint
  required. A straddling ball is not a pass.
- Inherited envelope-evaluation budget: 260,000 distinct (box,panel-count)
  evaluations for each named band.
- Campaign hard budget: 60 CPU-hours total across controls, primary band, and
  conditional close-outs, using `time.process_time` deltas; 72-hour wall
  ceiling from the fresh manifest's `created_utc`.
- Runtime: niceness 15; OMP/OpenBLAS/MKL/vecLib/NumExpr threads pinned to one.
- Refutation gate: any certified `J.lower() > d.upper()` interval stops the
  campaign and escalates to Main; a negative envelope margin alone is never
  interpreted as a lower bound on `J-d`.

Budget is checked before every envelope evaluation. An already-started Arb
call completes atomically; if it crosses the limit, the measured overshoot is
reported and the current cell/result is frozen. The exact frontier includes:
tile and c-pair, current cell, margin ball, depth, panels, full remaining
stack (every box and depth), stack count, panel histogram, envelope-eval
count, CPU used/remaining, and wall used/remaining.

## 3. Required controls after init, before band compute

### 3.1 Startup battery

Run A's full startup battery. Required result: 16/16 green, `aborted:null`,
and output byte-identical to A/B/C's frozen
`logs/startup_controls.json` (expected sha256
`211d8e9afd6491af62d7bb89282e586f300f5a32535320caa31aa342fc63beac`).
Any mismatch stops the campaign before further compute.

### 3.2 Instrument-identity reconciliation

Re-evaluate A's frozen tile-1 worst cell at c-pair
(4.083,4.115664):

- a0 ∈ [−13/16, −77/96];
- a2 ∈ [55/96, 7/12];
- panels 131072.

Required margin: A's frozen
`+1.04499016763434518123845756638e-6`. Pass iff the recomputed margin is
strictly positive and
`abs(recomputed - frozen).upper() < 2^-40`. Inputs, upper ball, margin ball,
difference ball, panel count, and CPU are frozen. Failure stops and escalates.

### 3.3 Near-miss discriminator

Use the OLD ratio-1.02 frozen midband c-pair solely for this control:

- cL = 3.16988277218086405264471296000;
- cR = 3.23328042762448133369760721920;
- a0 ∈ [−0.849527994791666666666667,
         −0.849446614583333333333333];
- a2 ∈ [+0.527669270833333333333333,
         +0.527750651041666666666667].

As the owner did, pad each printed endpoint outward by 1e-20. This SUPERSET
box gives `U_superset >= U_true`, so a positive margin certifies the true
cell; the negative result is a discriminator observation, not a `J>d` claim.
Compute exactly at 32768 then 131072 panels. Required directions:

- 32768 margin upper endpoint < 0;
- 131072 margin lower endpoint > 0.

Both full Arb balls and CPU are frozen. If it does not flip, the panel-cap
classification is falsified: STOP, run and freeze the adapted Campaign-B
D1–D6 battery on this cell (D1 32768 baseline; D2 131072 target; D3
2^19/2^21 panel-sup series; D4 disk-feasible even-sup diagnostic with the
firewall held; D5 tail; D6 2^19→2^21 skew), escalate the exact obstruction to
Main, and write no band conclusion.

### 3.4 Planted-failure control

On the same padded cell and same cL, deliberately widen the **c-tile** to
ratio 1.05: `cR_probe = 1.05*cL`. Evaluate at 131072 panels. Required result:
strictly NEGATIVE, i.e. margin upper endpoint < 0. This tests that the strict
margin predicate rejects an over-wide tile even at the primary cap. A
nonnegative/straddling result fails the preregistered control and stops before
the band.

Controls run in order: startup → reconciliation → near-miss flip →
ratio-1.05 planted failure. Each JSON is frozen before the next starts.

## 4. Primary DFS and 2^19 fallback

Run tiles 1..87 in order. For each cell, call the inherited primary ladder in
order and stop at its first strictly positive margin exactly as A does. A
failed-but-splittable box is split with A's `split_box`, children pushed in
A's reversed order, disk filter unchanged. Only when a box is unsplittable
and its 2^17 result is `panel cap reached` does the driver evaluate the same
box once at the explicit fallback rung n=524288, using A's
`cell_envelope(cL,box,524288)` and `strict_margin_lower(cR,upper)` unchanged.

- fallback margin > 0: leaf certifies at 524288; DFS continues;
- fallback margin ≤ 0 or straddles: leaf OPEN; tile OPEN; freeze exact
  frontier; stop the band immediately;
- no intermediate 2^18 claim is made; the fallback is exactly the registered
  2^19 rung.

Per tile write `logs/inproc_tile_NNN.json` with c-pair, verdict, root/box/leaf
counts, maximum depth, worst margin ball and cell, max panels, per-rung panel
histogram including 524288, fallback count, envelope evaluations, CPU, and
frontier if any. After every tile atomically rewrite
`logs/band_progress.json` with all completed rows, aggregate panel histogram,
CPU, exact certified prefix, and `next_tile`.

The run is resumable by `--start-tile N`. Existing completed tile JSONs are
read and checksum-checked; completed tiles are never recomputed. Resume
regenerates the c-table from A's module and asserts every stored c-pair.
Stop at the first failing tile. If the budget binds mid-tile, freeze the exact
frontier described in §2 and claim only the last completed contiguous prefix.

Primary-band verdict rule:

- PASS only if all 87 tiles finish with every leaf strictly positive under
  the primary or registered fallback rung;
- FAILURE TO CERTIFY at the first terminal OPEN leaf, with exact c-prefix and
  frontier;
- budget-bound means partial prefix only, never a band claim.

## 5. Conditional corrected-lineage close-outs

Only if [1.75,3.5] completes PASS with CPU and wall budget left, run under the
identical knob state:

1. `c_tiles('1.0','1.3')`, assert 33 tiles and exact cover;
2. then `c_tiles('1.45','1.75')`, assert 24 tiles and exact cover.

Use the same primary/fallback, DFS, budget, checkpoint, resume, and stop rules.
These are lineage close-outs, not frontier gaps. The sound direct-D.4 frozen
PASSes already stand independently:

- [1.0,1.3]: 14/14 old ratio-1.02 tiles, worst margin
  +9.495739844715277e-7;
- [1.45,1.75]: 10/10 old ratio-1.02 tiles, worst margin
  +3.10559377267758493230e-6.

A corrected-lineage PASS supersedes nothing; an incomplete/OPEN close-out
retracts nothing from the direct-D.4 instrument. Report exactly what the new
DFS earned.

## 6. Freeze, report, and exactly one terminal verdict

Before close, write `result.json` and a result section appended to the run
copy of this pre-statement. Report:

- exact certified c-prefix of [1.75,3.5];
- every processed tile's c-pair, verdict, worst margin ball, panel histogram,
  fallback usage, envelope evals, and measured CPU;
- aggregate panel histogram and CPU;
- exact frontier if incomplete: cell, margin ball, depth, panels, full
  remaining stack;
- close-out tables only if run;
- refutation gate status;
- exact scope not swept.

Then run `campaign.py freeze`, followed by one
`campaign.py close --verdict ...`:

- `FROZEN-CERTIFIED`: primary [1.75,3.5] PASS (conditional close-outs may be
  PASS, incomplete, or unrun for lack of residual budget; their state is
  disclosed separately);
- `FROZEN-NEGATIVE`: first terminal OPEN leaf after the registered 2^19 rung,
  or a required control failure, with exact obstruction;
- `FROZEN-INCONCLUSIVE`: CPU/wall budget binds before tile 87 or an external
  interruption prevents a terminal leaf/tile conclusion, with exact prefix
  and frontier.

No claim about [4.083,6.0] is made. Campaign C's tile-3 PASS at cap 2^19 and
its tiles-4–49 cost deferral remain exactly as frozen.
