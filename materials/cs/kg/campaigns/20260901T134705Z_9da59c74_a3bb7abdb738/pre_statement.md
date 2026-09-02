# Pre-statement — Campaign D (KgBand175): corrected-envelope certification of
# the never-run band [1.75, 3.5] at panel cap 2^17, with 2^19 fallback rung;
# then lineage close-outs of [1.0, 1.3] and [1.45, 1.75] on budget left

Agent `KgBand175`, UTC 2026-09-01. **Status: COMPUTE PENDING — this file is
committed BEFORE any numerical evaluation, tile run, or band run of this
campaign.** No run dir of this campaign exists yet (`campaign.py init` runs
AFTER this commit and mints the run dir, which then receives a byte-identical
copy of this pre-statement plus its sha256 and the source commit hash).
The only prior files of the minted directory will be the inherited instrument
copies committed alongside this pre-statement (byte-frozen, checksums
recorded in `checksums_addendum.sha256` at freeze time and verified in-run
against the frozen predecessor ledgers).

Owner scouting evidence on record (supplied; reproduce it, do not cite it as
authority): all six frozen open cells of the OLD ratio-1.02 midband run
close at cap 2^17 under the corrected instrument (margins −4.5e-6..−4.1e-8
at 32768 panels → +8.0e-5..+8.7e-5 at 131072 panels, measured on outward
SUPERSET boxes so they lower-bound true margins); pointwise
corrected ≥ direct-D.4 confirmed by measurement; and Main's tile-probe of
c=[3.18110465298457505, 3.20655349020845165] PASSed the full `certify_tile`
at 2^17 in 1781.8 s CPU (scouting only). Per the owner's caveat, those
measurements visit DIFFERENT cells than this campaign's ratio-1.008 tiling,
so they size the problem — they do NOT prove the band closes. The verdict is
earned only by this campaign's own DFS.

## 0. Predecessors and what is inherited, byte-identical (declared now)

Predecessors (all in `cs/kg/campaigns/`):

- **Campaign A** `20260831T231500Z_KgEvenBoxRerun_campA_corr_env` (frozen
  commit `8d2cb26`): corrected 2-D even-box envelope for [4.083, 6.0];
  tiles 1–2 PASS, tile 3 OPEN at cap 2^17; band FAILURE TO CERTIFY at 3/49.
- **Campaign B** `20260901T051620Z_KgCampB_subdiv_env` (frozen commit
  `f744dd2`): mq = t² subdivision recovers exactly 0 (crossing s* = 45.7034
  outside [0, 8]); quantified obstruction frozen (`obstruction_tile3.*`);
  D1–D6 decomposition battery frozen in `code/decompose_cell.py`.
- **Campaign C** `20260901T055808Z_KgCampC_pan2p19_campC` (frozen commit
  `95922c8`): tile 3 of [4.083, 6.0] re-run ENTIRELY at the single restored
  knob PANEL_CAP = 2^19 — 484/484 leaves certified, worst margin
  +5.06389732008514710540159866096e-7; Phase 2 (tiles 4–49) DEFERRED by
  owner decision on quantified cost grounds (min 1.97× over C's 40 CPU-h
  budget; zero Phase-2 CPU spent). The kgcloseB3 retraction STANDS; C
  claims no band PASS beyond tile 3.

The band [1.75, 3.5] is a NEVER-RUN band (listed OPEN in the README's
corrected band inventory since the 2026-08-30 correction; out of scope for
B §7 and C §10). Its binding obstruction is NOT established to be the panel
cap — [1.75, 3.5] has no measured obstruction like B's tile-3 finding — so
the primary cap here is A's frozen 2^17, not C's 2^19: raising the cap
across the board would spend budget C explicitly declined, and the two
kernels' cap behavior differs (this band's ladder starts at 1024 vs the
[4.083,6.0] band's 4096 — `panel_ladder`'s start bracket differs by c_left).

## 1. What is re-used, byte-identical (declared before compute)

From Campaign A's frozen `code/` (sha256 recorded in A's
`checksums.sha256`, reproduced verbatim in B's and C's addendum ledgers;
byte-verified on copy and re-verified IN-RUN before any cell math):

- `gate_corr_env.py`
  `94c07adf14fbd860c54c5a59f5430dafe9125a86107c1f0b0f85da441d9eee30`
  — envelope cores `e_panel_range`, `dist0_box`, `box_intersects_disk`,
  `root_boxes`, `split_box`, `can_split`, tiling `c_tiles` with
  `TILE_RATIO = '1.008'`, ladder `panel_ladder` (start by c_left bracket,
  ×4, unique, capped once at PANEL_CAP), floor `BOX_FLOOR_Q = 1/1024`,
  `BOX_DEPTH_MAX = 40`, budget semantics, `strict_margin_lower`,
  `evaluate_cell`, `_update_leaf_extrema`, `certify_tile` (DFS from the
  132-root cover), `serialise_tile`, `S_MAX = 8`, `PREC = 256` outward Arb,
  tail closed-form, N_SIDE = 12.

From Campaign C's frozen `code/` (byte-identical to A's, whose
`cf7d57fbb3052cb09fc719cafda6eb6de6830ce826b5f4ad1d8211cb7e965ac3` is the
recorded sha for BOTH; sha-verified at freeze here):

- `startup_controls.py` — the full 16-family battery (see §3).

From Campaign B's frozen `code/`: `decompose_cell.py`
  `291eddc1d76b5a60338855012629077a7d5799db9c7434e07b18ae19df6204ff`
  — the D1–D6 decomposition battery. Used ONLY in the escalation branch of
  Gate N (§4: on a failed flip, the obstruction is re-frozen with this
  battery and Main is escalated, exactly as B did, before any band verdict
  direction is written).

Shared dependency `20260831T082425Z_kg_direct_d4_restart/code/core.py` is
IMPORTED (not copied) from the frozen restart dir under name `d4core`,
exactly as A, B and C did. python-flint landmine discipline (owner rule):
every interval used as an interval is built in midpoint/radius or explicit
outward-endpoint form; `arb(lo, hi)` is never used as an interval
constructor. All arithmetic inside `Prec(256)`; 320 bits inside the
inherited startup battery, exactly as frozen.

**The only knobs of this campaign** (nothing else touched anywhere):

1. Panel cap PRIMARY = `gce.PANEL_CAP` left at its inherited value 131072
   (2^17) — the module constant, no patch at import.
2. Panel cap FALLBACK = 524288 (2^19), a REGISTERED rung applied ONLY when
   a leaf has stayed OPEN through the whole 2^17 ladder ('panel cap
   reached'): patched per-leaf, run reaching either PASS or final OPEN at
   2^19, then immediately restored to 131072. Hierarchically this extends
   A's ladder 1024→4096→16384→65536→131072 with one more rung
   524288 for that leaf only — the same ladder function semantics C used
   (one patch point, asserted before use, restored after).
3. c-range: `c_tiles('1.75', '3.5')` for the primary band;
   `c_tiles('1.0', '1.3')` and `c_tiles('1.45', '1.75')` for Phase-C
   close-outs (pre-registered domains).

No new envelope derivation, no retiling, no box-floor change, no depth
change, no S_MAX change, no precision change, no margin-predicate change,
no DFS order change, no split rule change, no budget-semantics change, no
root-cover change.

## 2. Band definition and frozen knobs (rule 7 + rule 16: fixed BEFORE the run)

- **Primary band:** `c ∈ [1.75, 3.5]` as exact decimals, tiled by the
  byte-shared `c_tiles('1.75','3.5')` at ratio 1.008, final tile clamped
  exactly to 3.5. Pre-registered assertions (checked in-run BEFORE tile 1
  computes): `len(tiles) == 87`; `assert_tile_cover('1.75','3.5', tiles)`
  passes (first tile [1.75, 1.764], adjacency gapless through the final
  tile [3.47251..., 3.5]).
- **Root cover:** inherited `root_boxes(12)` — exactly 132 boxes asserted
  intersecting the closed unit disk; binary longer-axis splits; maximum
  depth 40; minimum side 1/1024; no below-floor refinement. UNCHANGED.
- **Panel ladder:** `panel_ladder(c_left)` per A's frozen band-table
  function with PANEL_CAP = 131072: start 1024 for this band (1.75 ≤ c_left
  ≤ 4 ⇒ start = 1024), so `panel_ladder(c_left) ==
  [1024, 4096, 16384, 65536, 131072]`, asserted per tile before it runs.
  The registered 2^19 fallback EXTENDS this to [1024, 4096, 16384, 65536,
  131072, 524288] for the specific open leaf only.
- **Precision:** 256 Arb bits (outward), inherited; no floats in trusted
  paths; 320-bit startup battery exactly as frozen.
- **S_MAX = 8**, tail closed-form (inherited `d4core.tail_account`).
- **Verdict rule (byte-shared):** tile PASS iff every disk-intersecting
  leaf in the tile has `d(c_U).lower() − U(B, c_L).upper() > 0` in Arb
  (strict endpoint separation; margin ball's lower end is what must be
  > 0 — a straddling ball is NOT a pass). Band PASS iff every tile passes.
  An open leaf at floor depth with panel cap exhausted (2^17 then
  registered 2^19 rung) ⇒ tile OPEN ⇒ band FAILURE TO CERTIFY at that
  leaf, freezing the exact frontier.
- **In-band envelope-eval budget (byte-shared):** 260,000 distinct
  (box, panel-count) evaluations per named band — the inherited A budget
  parameter, per band: [1.75, 3.5] gets its own 260,000; each Phase-C
  close-out band runs under its own fresh 260,000. (The direct-D.4 full
  six-band sweep averaged ≤ 57,551 per band; 260,000 is the frozen
  predecessor config, unchanged.)
- **CPU budget (pre-registered):** 60 CPU-hours (`time.process_time`, sum
  across this campaign's verdict-bearing runs) against a 72-hour wall
  ceiling. On either binding, the in-flight tile freezes as its next
  checkpoint (OPEN, reason 'budget exhausted', exact frontier with margin
  ball/depth/panels/remaining stack) and the campaign closes honestly with
  the certified prefix. Wall ceiling hit ⇒ same freeze semantics.
- **Stop-at-first-failing-tile:** inherited verbatim (A §4 bullet: tile
  OPEN ⇒ band stop; no tile past a failing one is evaluated). The 2^19
  rung is NOT a new tile — it is part of exhausting the caps on the open
  leaf before the stop rule fires.
- **Refutation gate:** inherited verbatim (A §4 / B §8 / C §8): a certified
  `J.lower() > d.upper()` on any leaf of any tile = candidate refutation of
  the paper's Lemma D.2 — STOP the whole campaign, escalate to Main via hub
  with the exact cell and intervals, freeze nothing further. Predecessors
  never fired it; neither may we casually.

## 3. Startup gate (before ANY tile compute) — full battery re-run

Campaign C's frozen `startup_controls.py` (byte-identical to A's) re-runs
the full 16-check-family battery in-place, instancing
`logs/startup_controls.json` (sha logged at freeze):

1. `kernel_unit_norm`; 2. `kernel_e0_is_sqrt32` (e(0) = √(3/2) admitted);
3. `kernel_admitted_by_envelope`; 4. `kernel_envelope_tight`; 5.
   `cap1_counterfactual_rejected` (deficit 0.224744871391589…); 6.
   `dist0_exact_value`; 7. `three_d_inflation_is_looser` (3-D-inflated
   1.6542 vs corner-exact 1.3243); 8. `two_d_envelope_corner_exact_at_s0`;
9. `d3h_cert_sha256` (80e945589b…); 10. `Ko_Ke_identities`; 11.
   `B3_import_is_upper`; 12. `margin_straddle_rejected`; 13.
   `margin_tight_accepted`; 14. `margin_above_rejected`; 15.
   `pstar_closed_form_window`; 16. `envelope_admits_pstar`.

ANY failure ⇒ STOP, freeze nothing else, report — a failed control means
the inherited instrument is broken and every downstream number would be
void. (Families 11/15/16 anchor the instrument against restart constants
independent of this campaign's band; they run unchanged as provenance
anchors.)

## 4. Gates — run in order; each freezes its own JSON before the next runs

### Gate R — instrument-identity reconciliation (MUST pass before any band compute)

Reproduce Campaign A's frozen tile-1 margin under this campaign's
instrument. A's frozen tile-1 row: c-pair (4.083, 4.1156640000000000000),
worst certified margin **+1.04499016763434518123845756638e-6 ± 4.27e-36**
at worst certified cell a0 ∈ [−0.8125, −0.802083333333333333333333],
a2 ∈ [+0.572916666666666666666667, +0.583333333333333333333333]
(131072 panels, exact coordinates from A's frozen
`band_4p083_6_result.json`). This campaign rebuilds the cell from the
frozen 30-digit printed endpoints pad padded outward by 1e-20 (the printed
endpoints carry radii ≤ 3.34e-25, so the padded box is a certified outward
SUPERSET — its margin is a LOWER bound of the boxed-over true margin), runs
the full 2^17 ladder per A's frozen `evaluate_cell`/`panel_ladder` for the
[4.083, 4.115664] c-pair, and takes the final margin.

Gate passes iff:
- the recomputed margin's lower end is strictly positive, AND
- `(recomputed_margin − arb('1.04499016763434518123845756638e-6')).rad()
  < 2^-40` — i.e. the recomputed ball contains the frozen midpoint within
  2^-40 (the frozen ball's radius is 4.27e-36 ≪ 2^-40, so this pins the
  recomputation to A's number to far better than the superset pad's
  resolution).

Note on direction: the superset pad can only DEPRESS the margin (larger box
⇒ larger U ⇒ smaller d−U). The frozen certified margin was earned on the
exact same cell at the same 131072 panels; A's own DFS reported it at the
EXACT box (fmpq endpoints). To make the gate decisive in BOTH directions,
the gate evaluates the cell TWICE: once with the exact fmpq endpoints
of A's frozen cell (±exact halves from the 30-digit display, which show
da to 24 decimals and pin the fmpq exactly: −0.8125 = −13/16 etc.), and
once with the 1e-20 padded superset. BOTH must reproduce the frozen
midpoint within 2^-40 (the exact-box one must match bit-near-exactly; the
superset one may differ only by the measured pad effect, also < 2^-40 —
the pad's measurable effect scales with the box width, so it is checked
empirically rather than via a universal claim). If either check fails ⇒
STOP, freeze the discrepancy exactly (inputs, both margins, the failing
comparison) as the campaign's obstruction, escalate to Main, run nothing
further — the campaigns disagree and knob attribution is void.

### Gate G — gap orientation (INFORMATIVE, no verdict impact, freeze-only)

Dump into `logs/gate_orientation.json`: (a) the frozen midband run's
`tiles_open` boundary rows; (b) this campaign's `c_tiles('1.75','3.5')`
tile table with per-tile `panel_ladder`; (c) the owner's six measured
superset cells and this campaign's re-measurement of the cited flip cell
(consumed in Gate N). No pass/fail; the band verdict in no way conditions
on this gate.

### Gate F — planted-probe control (MUST be REJECTED): a deliberately widened
ratio-1.05 probe tile on Gate N's exact cell must certify NEGATIVE

Pre-registered counterfactual: take the Gate-N frozen cell and evaluate a
deliberately WIDENED probe box — a0 and a2 both padded outward by a factor
1.05 of the box half-widths (a strictly larger box whose sup over e is
looser and whose dist(0, ·) may cross differently). We assert the planted,
deliberately-loose probe fails to certify (margin lower end < 0) at the
same panel count where the unfrozen true cell certifies (Pad effect: the
uncorrected 1.05-widened box over-states the odd radius bound W_P and the
even sup E⁺; the envelope's own soundness is not in question — the PLANTED
PROBE's failure is the control). If the probe PASSES where the true cell
passes (margin lower end > 0 on the widened probe at 131072 panels),
report the discrepancy, freeze both margins, and PROCEED ONLY IF the
planted probe's cell remains OPEN (the probe's widened-odd-radius
artefact). Threshold is pre-fixed: the probe must be NEGATIVE at 131072
panels. A positive probe ⇒ STOP, freeze, escalate to Main, run nothing
further (this would mean the margin predicate is insensitive to obviously
worse boxes — an unsound envelope, not a panel-sup matter).

(The literal ratio-1.05 mechanism: probe box := center ± 1.05 · half-width;
asserted strictly-larger in both axes.)

### Gate N — near-miss discriminative flip (decisive; the owner's
classification test)

On the frozen midband open cell from
`20260831T145732Z_kg_direct_d4_midband_panelcache/logs/band_2p814765_3p499807_panelcache_result.json`
tiles_open[1] — whose `open_cell` is
a0 ∈ [−0.849527994791666666666667, −0.849446614583333333333333],
a2 ∈ [+0.527669270833333333333333, +0.527750651041666666666667] with
c-pair (3.16988277218086405264471296000, 3.23328042762448133369760721920),
frozen open margin −4.52960421419799285256261953577e-8 at 32768 panels —
evaluate at 32768 then 131072 panels with this campaign's instrument on
the 1e-20-outward-padded superset box:

- the 32768 margin MUST be negative (sign match; magnitude within 2^-40
  of the frozen −4.52960421419799285256261953577e-8 modulo the pad
  effect, and within the pad's measured signature);
- the 131072 margin MUST be positive (owner measured +8.69e-5; only the
  SIGN and strictly-positive-lower-end are gated; magnitude reported as
  COMPUTATIONAL-EVIDENCE).

If the sign does NOT flip: the owner's panel-cap classification is
falsified. STOP. Freeze the obstruction with Campaign B's frozen
`decompose_cell.py` D1–D6 battery (re-run against this cell: D1 baseline at
cap 2^17, D2 subdivision kill recovery, D3 panel-sup series
2^17→2^19→2^21, D4 disk-feasible even-sup diagnostic with its firewall
held, D5 tail, D6 within-panel skew), escalate to Main via hub with the
exact cell + intervals + decomposition, and write NO band conclusion.
The campaign freezes at that point with no Phase-B run.

## 5. PHASE-B — the band run (only after ALL gates hold)

Per-tile, in tile order 1..87, DFS per A's frozen `certify_tile` from the
132 roots; ladder asserted per tile before the tile runs;
per-tile artifacts `logs/inproc_tile_NNN.json` (A's `serialise_tile` row,
rule-15 fields, plus `process_seconds_tile`,
per-tile envelope evals, per-tile panel histogram); a running
`logs/band_progress.json` accumulating every completed row + aggregate
panel histogram + envelope evals + CPU + next-tile pointer. An interrupted
run loses at most the tile in flight; completed tiles are never recomputed
(C's phase-2 contract, mechanism proven in C's `logs/c_tile3_state.json`).

Resumability at tile granularity: `--start-tile N` evaluates tiles N..87 in
order; a stop file records the next unstarted tile index; a resume skips
completed tiles by reading existing `inproc_tile_NNN.json` rows (their
shas in the addendum ledger identify them) and continues from the stored
pointer. Per-tile c-pairs are regenerated from the byte-shared `c_tiles`
and asserted equal to adjacency — the tiling is deterministic, so a
resumed run sees the same tiles.

2^19 fallback rung on any OPEN ('panel cap reached') leaf: after the 2^17
ladder exhausts on that leaf, patch `gce.PANEL_CAP = 524288`, re-assert
the extended ladder `panel_ladder(c_left) == [1024, 4096, 16384, 65536,
131072, 524288]`, re-run the leaf's panel-loop from the ladder start
(the panel-envelope is deterministic per n; the only difference is the
additional 524288 rung), record the leaf verdict under the extended
ladder in the tile JSON as `fallback_2p19` with its margin ball,
restore `gce.PANEL_CAP = 131072`, and continue. The tile's verdict uses
the network of leaf verdicts under their respective knob states.

CPU accounting: `time.process_time` deltas per tile, accumulated in
`band_progress.json`; at the 60 CPU-h budget the in-flight tile is frozen
as 'budget exhausted' (OPEN, exact frontier: cell, margin ball, depth,
panels, remaining stack), and the campaign closes with the exact certified
prefix. Stop-at-first-failing-tile is inherited verbatim: any OPEN leaf
after the 2^19 rung ⇒ tile OPEN ⇒ band FAILURE TO CERTIFY, exact frontier
frozen, no tile past it evaluated.

Refutation gate check: inherited verbatim, per leaf, on every margin.

## 6. PHASE-C — close-out bands, ONLY on a [1.75, 3.5] PASS with budget left

Both bands under the IDENTICAL knob state (module-inline 2^17 ladder at
each tile's own band table, one registered 2^19 fallback per open leaf,
byte-shared instrument, per-band fresh 260,000 envelope-eval budget):

- `[1.0, 1.3]` — `c_tiles('1.0', '1.3')` (expected 33 tiles at ratio
  1.008 — asserted in-run; the direct-D.4 baseline ran 14 ratio-1.02
  tiles, a different tiling, so the count is NOT expected to match — the
  assert is `assert_tile_cover` + documented tile count of THIS tiling);
- `[1.45, 1.75]` — `c_tiles('1.45', '1.75')` (expected 24 tiles; the
  direct-D.4 baseline ran 10 ratio-1.02 tiles — same caveat).

Both bands ALREADY hold valid PASS certificates under the sound direct-D.4
instrument (A-heredity `20260831T082425Z_kg_direct_d4_restart` records
[1.0,1.3] 14/14 PASS worst margin +9.495739844715277e-7 and [1.45,1.75]
10/10 PASS worst margin +3.10559377267758493230e-6). A PASS under a
pointwise-LOOSER upper-bound envelope (this campaign's corrected
envelope) is a VALID proof, so these close-outs are LINEAGE close-outs,
NOT frontier gaps — they supersede nothing; the frozen direct-D.4 PASSes
stand on their own instrument. They are run to bring the two bands under
the SAME corrected envelope as [1.75,3.5] so the three never-run
corrected-envelope bands form one consistent instrument, and to spend
budget C's owner decision declared available here (after the hard band,
which is the priority).

Phase-C aborts cleanly (no close-out band claim) if the primary band
fails, or if either close-out band opens at 2^17 + registered 2^19 rung —
a close-out is reported under its own band, never as a [1.75,3.5]
result.

## 7. Labels

- MACHINE-VERIFIED: any margin/sign/interval claim from outward Arb
  endpoints at the stated domains (256-bit outward), scripts frozen in
  `code/`.
- COMPUTATIONAL-EVIDENCE: timings (time.process_time / resource.getrusage
  only, rule 17e), wall-clock costs, cross-campaign comparisons (owner
  scouting probes, Main's tile probe), any projected cost.
- REPORTED: A's tile-1/2 margins; the frozen midband cells and margins; A's
  direct-D.4 PASS bands; C's tile-3 result.
- FAILURE TO CERTIFY is a result and is reported honestly; PASS requires a
  strictly positive margin-ball lower end on EVERY leaf of EVERY tile of
  the band under the registered knob state.

## 8. Scope (rule 7, fixed now)

Swept here (conditional on gates passing): band [1.75, 3.5] over the full
2-D even-coefficient disk (132-root cover, all compatible odd parts, [0,8]
panels + closed tail, 256-bit outward Arb) at primary cap 2^17 with the
one-leaf-registered 2^19 fallback rung per open leaf; bands [1.0, 1.3]
and [1.45, 1.75] under the identical knob state, only on a primary
[1.75,3.5] PASS and only with CPU budget left. NOT swept: any re-run of
[4.083, 6.0] tiles (C's tile-3 verdict stands under its own 2^19
instrument; tiles 4–49 are C's named next campaign at budget ≥ 160 CPU-h
— out of scope here); the c-domain outside [1.0, 1.3] ∪ [1.45, 1.75] ∪
[1.75, 3.5]; any cap value applied other than 2^17 (primary) and 2^19
(registered fallback, open-leaf only); any coefficient-space envelope
change; the mq = t² subdivision mechanism (B's Phase-1 instrument —
measured recovery exactly 0; NOT adopted here); degree > 3 families;
non-unit cubics; Constructions outside the unit-cubic threshold family;
the paper's C1/splice/C3 certificates; Gate-C septic/Krivine families.
The startup battery sweeps only its 16 planted asserts; Gate F sweeps
only its one planted widened-box cell; Gate N sweeps only its one
named midband cell at two panel counts.

## 9. In-run assert order (fixed now; enforced by the driver before tile 1)

1. Byte-verify sha256 of the run-dir copies of `gate_corr_env.py`,
   `startup_controls.py`, `decompose_cell.py` against A/B/C's frozen
   leadgers (ABORT on any mismatch). Assert `gce.PANEL_CAP == 131072`
   pristine at import.
2. `len(c_tiles('1.75','3.5')) == 87`, `assert_tile_cover` passes;
   for tile 1's c_left: `panel_ladder(c_left) == [1024, 4096, 16384,
   65536, 131072]`.
3. Startup battery: 16/16 green (`logs/startup_controls.json`).
4. Gate R (reconciliation): both sub-checks pass
   (`logs/gate_reconciliation.json`).
5. Gate F (planted probe): widened 1.05 box certifies NEGATIVE at 131072
   panels (`logs/gate_planted_probe.json`).
6. Gate N (near-miss flip): negative at 32768, positive at 131072
   (`logs/gate_near_miss_flip.json`).
7. Only after 1–6 all hold: Phase-B runs.
8. Phase-B completed at PASS: Phase-C runs with its own band-level gates.

## 10. Freezing & close

`PYTHONDONTWRITEBYTECODE=1` and `sys.dont_write_bytecode = True` before
importing any module from a frozen campaign dir; no `__pycache__` is ever
written inside the run dir (the run driver enforces this, and
`campaign.py freeze` re-checks pristine state). At campaign end:
`campaign.py freeze` (pins sha256s.txt over the run dir), then
`campaign.py close --verdict <exactly one>`. The verdict is chosen as:

- **FROZEN-CERTIFIED** if tile 87 certifies (band [1.75,3.5] PASS) AND
  Phase-C ran per §6 and closed both close-out bands — with all gates
  green and no refutation-gate fire; OR if tile 87 certifies but budget
  or gates leave Phase-C un-run: still FROZEN-CERTIFIED for the primary
  band, with the close-out deferral disclosed in the campaign result and
  in `state.json` (a partial close-out is NOT a secondary band claim).
- **FROZEN-NEGATIVE** if the primary band stopped OPEN anywhere (with the
  exact frontier), or the refutation gate fired (with escalation), or a
  registered gate failed (with the gate artifact staged). Any remaining
  Open or gate-failure supersedes a partial band claim.
- **FROZEN-INCONCLUSIVE** if the CPU (60 CPU-h) or wall (72 h) budget
  bound before tile 87 — freeze the exact frontier (certified prefix,
  in-flight tile state with its open-leaf margin ball/depth/panels/
  remaining stack) and close, with the certified prefix stated
  explicitly in the result, never left implicit.

Verdict semantics identical to A/B/C: partial results reported exactly;
no band PASS claimed beyond what the DFS earned; a partial band is a
partial band.
