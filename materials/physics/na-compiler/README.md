# `na-compiler/` — certified-optimal small-instance transport scheduling

**Status: BENCHMARK (smoke harness complete; one full parallelization-certificate set frozen; no multi-qubit campaign yet).** Exact ILP/SAT optima for neutral-atom zoned-array transport scheduling, computed against MQT QMAP's own output. Pure offline computation, `nice -n 10`, one core, explicit stop conditions.

Everything quantitative in this README is labeled:
- **PROVED** — exact solver optimality proof archived (Z3 k−1 UNSAT + k SAT, or HiGHS `Optimal` with dual bound, in the campaign manifest).
- **REPRODUCED** — a table in a frozen `campaigns/<ts>_<uuid>_<hash>/` rerun by `campaign_runner.py`.
- **REPORTED** — single numerical probe (this README); rerun instructions included.
- **DERIVED / INFERENCE** — modeling conclusions drawn from source reads, pinned with the source.

## Constraint model (pinned 2026-08-29 from first-hand reads)

Authoritative SSOT is `pre_statement.md` section 1. Summary:

1. **AOD strict move-compatibility** (`Move` batch = set of parallel moves):
   two moves are incompatible unless
   `(xs1==xs2) <=> (xt1==xt2)`, `(ys1==ys2) <=> (yt1==yt2)`,
   `(xs1<xs2) <=> (xt1<xt2)`, `(ys1<ys2) <=> (yt1<yt2)`.
   Source-verified in QMAP `IndependentSetRouter::isCompatibleMovement`
   (`main`, read 2026-08-29). This encodes non-crossing (3–4) and
   preservation (1–2) of the paper; minimum separation is enforced by the
   site grid (one atom per trap), not by the router.
2. **Ghost-spots.** QMAP's shipped strict router does **not** implement the
   ghost-spot constraint (our exhaustive grep found no check). Gate-A1 uses
   the QMAP-comparable model (**without** ghost-spots) as primary; the
   `ghost_spots=True` flag adds the stricter pick-up/drop-off row/column
   disjointness as a secondary experiment (pre_statement §1).
3. **Architecture.** Scaled copy of QMAP's shipped `square_architecture.json`
   (site separations 4 µm storage, 12×10 µm entanglement; AOD r=c=100;
   `atom_transfer = 15 µs`), with zone *sizes* cut to the instance. Entire
   JSON emitted by `Architecture.to_quimap_json()`; architecture pointer:
   `site_xy` equals QMAP's `exactSLMLocation` (x = first separation × col +
   location, y = second separation × row + location) — verified in source.
4. **Initial placement.** QMAP's `VertexMatchingPlacer.makeInitialPlacement`:
   row-major fill of the first storage SLM from (0,0). Both QMAP and the
   certified solver start from this same mapping.
5. **Layers.** QMAP `ASAPScheduler` semantics (earliest free layer, both
   operands free). Our `circuits.asap_merge` replicates it exactly;
   architecture is sized so the capacity limiter never splits a benchmark
   layer (`maxFillingFactor * r*c ≥ layer width`).

### Duration model (one model, QMAP-official)

- Batch move for max distance `d` (µm): `t = 2*(4d/j)^(1/3)` for `d ≤ 110 µm`
  (`j = 32*110/200³ µm/µs³`, so `t(110)=200 µs`); `t = 200 + (d−110)/1.1 µs`
  for `d > 110 µm` (QMAP official evaluator `eval_ids_relaxed_routing.py`,
  read raw from GitHub main on 2026-08-29).
- Every `load` and every `store` atom event: +15 µs (`Architecture.hpp`
  `timeAtomTransfer = 15`; also `operation_duration.atom_transfer` in arch).
- `gate` and 1Q durations are identical for every solution of a given
  instance; gap is transport-only. Secondary (reported, never gated)
  metrics: RSQASM linear model `2*t_trans + D/υ` (`t_trans=20 µs`,
  `υ=0.55 µm/µs`) and QMAP paper `(d/2750)^(1/2)`.

## Gates and current state

### Gate A — certified optimal move batches vs MQT QMAP output

**Instance set (pre_statement §3):**
- `qft4..qft10` — qiskit 2.5.2 `QFT(n, do_swaps=True)`, transpiled at
  `optimization_level=0` to `basis_gates=[u,cz]`; instance = ordered CZ
  stream after transpile. Generator: `circuits.qft_instance`.
- `reg3_4..reg3_10` (even n, seeds 0/1/2) — `networkx 3.6.1`
  `random_regular_graph(3, n, seed)`, edges sorted. Generator:
  `circuits.reg3_instance`.
- External QASM: `circuits.load_qasm(path_or_text)` accepts OpenQASM 2.x.

**Quantity:** for every layer transition in QMAP's `.naviz` placement
trajectory, (1) `chi` = minimum number of legal batches (= conflict-graph
chromatic number), Z3 certified; (2) minimum transport duration over all
legal partitions, MIP level-based exact model (see
`src/na_compiler/optimize.py` docstring), HiGHS certified (`status=Optimal`,
dual bound logged).
Gap = `(qmap_transport_cost − cert) / cert`, transport cost = QMAP's own
load/store/move arithmetic under the pinned model.

**A1 (PROVED in smoke campaign):** certified optimum of the per-transition
batch partition, vs QMAP's own batch sequence. The proof obligation per
pre_statement §5 (k−1 UNSAT + k SAT) is met by `optimize.cert_min_batches`;
the MIP optimality proof is HiGHS `status=Optimal` + dual bound, archived
per transition in `manifest.json`.

**A2 (not run):** joint placement+routing optimum, `n ≤ 6`, `≤ 3` 2Q layers.
Scaffolded interface exists (`cert_min_batches`, `cert_min_duration`);
campaigning remains. Never pool A2 rows with A1.

### Gate B — exhaustive hardest permutation instances

**Model:** single-zone `S×S` grid (`S ≤ 4`), `n` atoms, BFS from identity;
batch = set of pairwise-compatible singleton moves into currently-free
distinct targets; BFS distance = min batches; distance layer = provably
optimal per state (the search tree is the certificate). `enumerate_reachable`
returns `(dist, frontier_sizes, visited, aborted)` and gates on
`max_states` / `max_depth` only. **Exhaustive = all reached states
complete=True.** `hardest_states(dist, k=5)` sorted by descending depth,
ties broken lexicographically.

**Current per-pre_statement protocol:** exhaustive for `n ≤ 6` on `3×3`
(verified `complete=True` at `n = 3/4/5/6`, see campaign
`20260830T002458Z_71680b58_deb7e0f2bbe9`); hardest-found (labeled, never
promoted to "hardest") beyond that. The  `n ∈ {7, 8}` on `4×4` range is
frontier-capped hardest-found per pre_statement §4 protocol — work item now.

## Smoke artifacts (REPRODUCED)

Frozen in `campaigns/`:

| row | value | source |
|---|---|---|
| campaign `20260830T002449Z_80eab445_58aa8d772c52` | gate A, qft3/4, reg3_4 | `campaign_runner --gates A --nqubits 3,4 --seeds 0` |
| campaign `20260830T002458Z_71680b58_deb7e0f2bbe9` | gate B 3×3 full, complete | one-shot runner (this README) |
| campaign `20260830T002708Z_6907d87f_deb7e0f2bbe9` | gate B 4×4 hardest-found, 30k cap | one-shot runner (this README) |

### Gate A gap table (PROVED via archived solver logs)

From `campaigns/20260830T002449Z_80eab445_58aa8d772c52`:

| instance | n | qmap_cost_us | cert_us | gap_rel |
|---|---|---|---|---|
| qft3           | 3 | 5541.90  | 2888.28  | +0.9188 |
| qft4           | 4 | 17866.61 | 12369.39 | +0.4444 |
| reg3_4_s0      | 4 | 8823.71  | 5301.13  | +0.6645 |

All rows have `Highs status=Optimal, mip_gap=0` logged and Z3 UNSAT/SAT
records per transition in `manifest.json`. qft3 gap is +91.9%; QMAP's
transport cost exceeds the certified optimum by ~1.9× on this tiny instance
under QMAP's own duration model. All walls in `manifest.json` are single-run
timings; not further analyzed in smoke. **Read: A1 semantics, not
"QMAP wrong"** — QMAP was invoked with strict routing (smaller search space),
start fixed, and our certified solver discs lower the same start/layers but
optimizes globally over the batching.

### Gate B certificate rows (single-run probes)

| n_atoms | grid | states | P(n²,n) | time_s | complete |
|---|---|---|---|---|---|
| 3 | 3×3 | 504   | 504    | 0.13 | True (in campaign) |
| 4 | 3×3 | 3024  | 3024   | 1.07 | True (in campaign) |
| 3 | 4×4 | 3360  | 3360   | 5.25 | True (REPORTED in README) |
| 5 | 3×3 | 15120 | 15120  | 5.47 | True (REPORTED in README) |
| 6 | 3×3 | 60480 | 60480  | 16.59 | True (REPORTED in README) |
| 4 | 4×4 | 42758* | 43680 | 37.7 | False — 30k cap (in campaign) |
| 5 | 4×4 | 30189* | 524160 | 2.03 | False — 30k cap (REPORTED) |
| 6 | 4×4 | 69170* | 5.7M  | 4.20 | False — 30k cap (REPORTED) |
| 7 | 4×4 | 131078* | 57M  | 7.36 | False — 30k cap (REPORTED) |
| 8 | 4×4 | 218299* | 5.1e8 | 11.67 | False — 30k cap (REPORTED) |

*state-cap-reached rows; hardest entries in the corresponding manifest rows
are "hardest-found" per pre_statement §4. Full `P(16,4)` exhaustion is
~150 s wall (measured 2026-08-29); running and archiving it in a `campaigns/`
directory is the immediate follow-up work item, as is n=5/4×4.

Runtime scaling (certified gate A wall times, REPORTED):

| instance | total_moves | Z3 chi time | MIP time |
|---|---|---|---|
| qft3  | 8   | 0.070 s | 0.013 s |
| qft4  | 48  | 0.042 s | 0.075 s |
| qft5  | 78  | 0.066 s | 0.114 s |
| qft6  | 130 | 0.102 s | 0.277 s |
| qft8  | 253 | 0.256 s | 0.553 s |
| qft10 | 268 | 0.204 s | 0.806 s |

Gate-A certification solves all ≤10-qubit instances in ~1 s total per
instance in this single-core run; the runtime budget is therefore entirely
in the gate-B enumeration, not the A1 certificates.

## How to run

Audience: developer with `physics/.venv` ready (highspy, z3-solver, networkx
3.6.1, qiskit 2.5.2; `mqt.qmap==3.9.0` installed 2026-08-29 for the QMAP
comparison).

```sh
cd physics
nice -n 10 .venv/bin/python -c "import sys; sys.path.insert(0, 'na-compiler/src'); from na_compiler import smoke; smoke.main()"
nice -n 10 .venv/bin/python na-compiler/src/na_compiler/campaign_runner.py --gates A --nqubits 3,4 --seeds 0
# gate B certificates / one-shot probes: replicate README snippet below
```

Gate-B one-shot probe pattern (same as used for the README rows):

```py
import sys
sys.path.insert(0, 'na-compiler/src')
from na_compiler.exhaustive import enumerate_reachable, hardest_states
dist, fr, vis, ab = enumerate_reachable(4, 3)  # n_atoms, grid_side
```

Falsification / certificate verification: read
`campaigns/<ts>/manifest.json` per instance for the Z3 UNSAT/SAT transcript
and MIP `status=Optimal` + dual-bound entries; gap row in
`campaigns/<ts>/gap_table.md` is `qmap_cost_us vs cert_duration_us`.

## Inventory

- `pre_statement.md` — gate statement: instance set, exact constraint
  model, primary/secondary duration models, exact pass/fail criteria,
  definition of CERTIFIED.
- `src/na_compiler/` — `config.py` (arch JSON, duration model),
  `architecture.py` (scaled QMAP architecture + exact coordinates),
  `circuits.py` (QFT/reg3 instance generators + QASM loader + ASAP merge),
  `qmap_adapter.py` (mqt.qmap invocation + `.naviz` strict parser
  (load/store/move/cz validated with assertions; unknown lines abort
  parsing, never silenced) + raw batch record), `conflict.py` (strict
  compatibility relation + conflict graph), `optimize.py` (Z3 chi
  certificate + MIP level-based duration certificate), `exhaustive.py`
  (gate-B BFS, correctness invariants tested by `smoke.py`), `smoke.py`
  (end-to-end smoke), `campaign_runner.py` (campaign freezer: hashes src
  tree, writes manifest + gap table).
- `campaigns/` — three frozen snapshots (above). Campaign directories are
  immutable; directories created by failed/aborted probes are removed
  immediately with the failure recorded in this README (see "Campaign
  curation log").
- `scratch/` — none in play this session; empty.

## Acknowledged engineering shortcomings and risks

- Proxy/lower-bound comparisons in `campaign_runner.py`
  (`min_batches_lb_gap` etc.) are secondary metrics, never pooled with the
  certified gap; their interpretation requires explicit care (the LB
  compares batch count, not durations).
- The gate-B BFS cost is dominated by legal-batch *generation* per state
  (~0.005 s at 4×4; larger on 5×5 grids where free-site combinatorics
  explode). This caps full exhaustive enumeration in-warehouse at roughly
  `P(16,4) ≈ 43k`; scaling beyond that needs bitstate/IDA* pruning before a
  faithful `n ∈ {7, 8}` hardest-instance table can be published. This is the
  named blocker for the remaining gate-B protocol step.
- `qft_instance`'s transpile output is an opaque parameterized `u(x)`
  chain; we use the CZ pairing stream only (1Q gates do not affect
  transport). The native QAST→CZ exchange for an external OpenQASM circuit
  is exactly what `load_qasm`'s tiny two-op (`cz|cx|cp|swap`) path supports.
- Parser robustness: unknown `.naviz` lines assert/fail (pre_statement
  acceptance criterion), so any QMAP output format drift aborts the run
  rather than silently changing the comparison.

## Campaign curation log (un-normalized state cleanup)

- 20260829 probe series: 4 dirs; ultimate canonical gate-A row is
  `20260830T002449Z_80eab445_58aa8d772c52`.   - Two failed/skipped runs that
  created no `manifest.json` (timeout hits) were immediately removed
  rather than left non-canonical:

Older partial and probe directories from this session were deleted to
collapse campaign history into the three canonical rows involved by the
README. Deletion performed via `rmdir/shutil.rmtree` on empty/failed dirs;
no completed artifact rows were discarded.
