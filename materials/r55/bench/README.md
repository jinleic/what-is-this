# Gate-3 frozen benchmark (R(5,5) ≤ 45 campaign)

Frozen 2026-08-16, **v2 (node-count instrumented)**. Single source of truth:
[`spec.json`](spec.json) — machine, toolchain SHA-256s, per-item baselines
(wall, DFS nodes, SAT conflicts/decisions), experiment log, and the current
`status` field. Harness: [`../src/bench_stratum.py`](../src/bench_stratum.py)
(fail-closed: any canonical mismatch, or an empty result, exits nonzero).

## What one item does

`bench_stratum.py n e` reproduces the full edge-extremal stratum R(4,5,n,e)
end-to-end from validated catalogs (max-degree gluing at the proved
min-degree bound; C DFS engine + exact AllSAT tail for deferred pairs; every
graph re-verified by the trusted checker), then canonically diffs against the
published census member, which is unsealed only at the end.

## Items (v2 instrumented baselines, 8 workers)

| kind | (n,e) | classes | wall | DFS nodes | SAT conflicts |
|---|---|---|---|---|---|
| dev | (12,48) | 1 | 0.6 s | 4 | 0 |
| dev | (13,52) | 10 | 1.1 s | 11,963 | 0 |
| dev | (13,53) | 2 | 0.9 s | 8 | 0 |
| dev | (14,60) | 1 | 1.2 s | 148 | 0 |
| dev | (15,66) | 1 | 2.0 s | 84,652 | 0 |
| dev | (16,71) | 138 | 7.4 s | 49,309,494 | 27,703 |
| dev | (16,72) | 5 | 3.7 s | 9,147,704 | 1,131 |
| dev | (17,78) | 86 | 26.9 s | 475,709,989 | 714,376 |
| dev | (17,79) | 1 | 7.3 s | 98,625,550 | 43,044 |
| held-out | (18,85) | 74 | 314.5 s | 6,780,102,449 | 36,195,685 |
| held-out | (21,107) | 31 | 4,937 s | 71,184,533,784 | 404,295,016 |

Node-count semantics: `dfs_nodes` = C-engine DFS nodes summed over every
pair; `sat_conflicts` = Cadical195 accum conflicts over the exact SAT parts.
The (21,107) case-split driver additionally reports case-1 extension-SAT
conflicts/decisions (its case 1 is a pure SAT pipeline over 305,671 base
graphs — the SAT metric *is* its search-size metric).

All 11 strata have been reproduced with exact canonical-set equality; the
(21,107) instrumented rerun (2026-08-16, rc 0) reproduced the gate-2 census
bit-for-bit (209 raw, identical degree-class breakdown). Gate 3 is
**COMPLETE** — see `spec.json:status`. The (21,107) row includes case-1
conflicts (10,165,906 over 305,671 extension SATs); per-case detail in
`spec.json:held_out_items[0].baseline_detail`; its decomposition driver
`../src/case_split_21_107.py` is SHA-pinned in `spec.json:toolchain`.

## Gate-4 protocol

An improvement claim MUST: (a) keep every dev item MATCH; (b) improve
wall/deferral metrics on dev; (c) transfer in a same-session certified A/B on
a held-out item (both arms MATCH); (d) log the experiment in `spec.json`
(negative results included); (e) state any toolchain delta (SHA-256 bump +
dev-suite re-MATCH). Search-size metrics are required diagnostics, not a
substitute for the wall/deferral promotion gate.

**Protocol erratum, 2026-08-16:** criterion (b) previously said “wall or
search-size”, conflicting with the authoritative `spec.json:purpose`
wall/deferral rule. The wording above restores the frozen SSOT; held-out
transfer remains mandatory before any improvement claim.

## Experiment log

1. **cand1_static_node_limit_100k — REJECTED by held-out.** Dev suite −17 %
   (dominated by (17,78): 23.8→16.1 s), but (18,85) A/B: 302.5→323.7 s (+7 %).
   Both arms MATCH everywhere; correctness never at risk. Insight: the
   DFS→SAT handoff point is pair-dependent → candidate 2 = per-pair adaptive
   handoff policy.
2. **cand2_dfs_tail_concentration_and_measured_sat_cost — candidate 2
   deprioritized on measured evidence.** Two corrections then a measurement:
   (a) the v2 sparse CSVs are a *censored sample* — completed zero-solution
   pairs ((17,78) d=10: 13,064 completed rows vs 129 deferred) emit no row;
   fixed with opt-in `GLUE_CSV_ALL` (v2.1 build, pinned separately in
   `spec.json:toolchain.candidate2_instrumentation`; default mode
   byte-identical to v2, dev suite re-MATCHed with identical node counts).
   (b) The resulting `../outD/cand2_data/` (19 dev slices, 35,576 rows;
   held-out excluded) is a **DFS-tail concentration dataset**, not a handoff
   oracle: it holds DFS nodes only, deferred rows stay right-censored at
   budget+1, and SAT cost is not bounded by the DFS node cap.
   **Measured SAT cost on the candidate tail** ((17,78) d=10, top-2 %
   completed pairs = 261 pairs holding 51 % of completed-pair nodes):
   13.392442 s summed pair wall (mean 51.312 ms, median 8.126 ms,
   max 1.881356 s; 994,130 conflicts, 1,347,157 decisions; all 261 SAT
   solution counts equal DFS ground truth). Exact pair IDs and per-pair
   timings/stats: `../outD/cand2_data/17_78_d10_top2pct_sat.csv`
   (SHA-256 `521ca1f9…d21b`). Tail DFS 16.2 s and the perfect-predictor
   2.8 s / ~6 % saving are **estimates** from the measured full-slice
   average rate (376,004,258 nodes / 48.1 s), not an end-to-end forced-defer
   wall measurement; a feature-based predictor keeps less.
   **Scheduling lever largely exhausted →
   candidate 3: structural pair prefilter (R(3,5,m) neighborhood edge
   windows / second identity family) attacking node mass directly.**
   Baseline provenance unchanged: every frozen baseline remains attributed
   to the archived v2 engine (`spec.json:toolchain.engine`).
3. **cand3_local_edge_window_prefilter — REJECTED by held-out.** A proved
   pre-DFS relaxation intersects exact intrinsic R(3,5,m) neighborhood and
   R(4,4,m) anti-neighborhood edge windows. Small exhaustive replay:
   43,200 valid cones / 562 pair-edge cases, zero false rejections. The first,
   fixed-point implementation rejected 6,111/35,576 admitted dev rows and
   11,072,343 completed DFS nodes, but its overhead lost wall; it was
   superseded by a one-pass outer relaxation. The fixed fast variant rejected
   5,160 rows / 10,272,379 completed nodes (no positive or deferred dataset
   row), and all nine dev items MATCH: 632,889,512→622,617,100 DFS nodes
   (−1.62 %), 191→191 deferred. Same-session full-dev ABBA measured
   43.45→43.15 s (−0.7 %). On untouched held-out (18,85), both arms reproduced
   74 classes MATCH, but wall regressed **252.0→266.9 s (+5.9 %)**; nodes
   changed only 6,780,102,449→6,779,182,394 (−0.014 %) and deferrals stayed
   3,148. No transfer; rejected. Exact sources, binaries, hashes, replay
   scripts, A/B metrics, and corrections:
   `../outD/candidate3/measurement.json` (SHA-256 `7f9217dd…a24`).
   Live `../src/glue_census2.c` was byte-restored to v2.1 after archiving the
   rejected variants.
