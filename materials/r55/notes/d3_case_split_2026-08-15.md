# D3 census by min-degree case split (gate 2, 2026-08-15)

Goal: reproduce the AM D3 stratum census R(4,5,21,e=107) — published: 31 iso
classes — independently, after both frozen engines hit a wall at this stratum.

## Why the frozen engines stalled (measured, this session)

| engine | config | measured | projection for (21,107) |
|---|---|---|---|
| `srcB/glue2` (DFS) | as frozen | > 150 s/pair on a d=13 probe (timeout) | CPU-years |
| `src/glue_census2_c` | GLUE_DMIN=7, limit 1e6 | 190–284 ms/pair, ~100% deferred | 87 CPU-h wasted DFS, then SAT |
| `src/sat_pair.py` | dmin=7 | 18.4 s/pair (d=11), 3.4 s/pair (d=13) | ~8,000 CPU-h |

Root cause: at (21,107) the cone is dense (47–72 of ≤ 99 cross-slots), the
worst regime for both DFS and AllSAT over the raw spec constraints.
Also found: `run_stratum_c.py` defaults `GLUE_BIN` to v1 `glue_census_c`,
which never emits deferred-pair markers — the SAT fallback silently never
runs unless `GLUE_BIN` points to v2. Recorded here as an operational trap.

## The case split (proofs)

Let G be in the census, delta = min degree, Delta = max degree.

1. delta >= 7: deleting a min-degree vertex leaves an R(4,5,20) graph with
   107 - delta edges <= E(4,5,20) = 100 (order-20 extremal value, verified
   from the gate-1-validated archive).
2. delta <= 10: min <= mean = 214/21 < 10.2.
3. **Case 1 (delta <= 9):** G - u lies in the published complete classes
   r4520.{100,99,98}.g6 (1 + 822 + 304,848 graphs). G is a one-vertex
   extension: N(u) triangle-free-hitting constraints (no K3 of G-u inside
   N(u); N(u) hits every I4 of G-u), |N(u)| = 107 - e'. Exact iff. 20-var
   AllSAT per base graph, ~0.5 ms each.
4. **Case 2 (delta >= 10):** audited max-degree gluing spec with the
   within-case prune GLUE_DMIN=10. Summing per-vertex cone lower bounds
   (deg_G(a) >= 10 over A; deg_G(b) >= 10 over B) yields proved pair windows:
   - d=11: e(B)+3 <= e(A) <= e(B)+6  -> 244,323 of 1,543,605 pairs survive (15.8%)
   - d=12: e(B)+13 <= e(A) <= e(B)+15 -> 508 of 24,948 pairs
   - d=13: forces e(B) <= 3 < 5 = min over R(4,4,7): **case empty, proved**
5. Cases are disjoint (delta <= 9 vs >= 10) and exhaustive; the census is the
   deduplicated union. Published reference held out until the final labelg diff.

## Effect of the case prune (measured)

- glue_census2_c, GLUE_DMIN=10, limit 1e6: d=12/13 ~2 ms/pair with zero
  deferrals; d=11 ~124 ms/pair with ~21% deferrals.
- sat_pair with dmin=10: 99 ms/pair vs 18,400 ms/pair at dmin=7 — the
  within-case min-degree bound speeds the SAT tail by ~186x.
- Whole-stratum plan cost: case 1 ~1 min (8 workers), case 2 d=12 ~33 s,
  d=11 ~1–2 wall-hours at 8 workers. Versus CPU-years for the frozen configs.

## Run results (this session)

- Warm-up (Lemma 5.1): R(5,4,17, e<=59) = 7,147 graphs reproduced exactly
  from archive classes r4517.{77,78,79}; 0 property violations; every graph
  has a vertex of degree >= 8 as claimed.
- Case 1: 45 raw solutions, all at delta = 9 (none at 7 or 8), including
  both raw copies of the unique Delta=12 census graph; every output
  re-verified (K4-free, I5-free, e=107, min-degree invariant).
- Case 2, d=12: 0 solutions (the Delta=12 graph has min degree 9 and is
  case-1 territory — confirmed).
- Case 2, d=11: DFS pass 640 bucket-chunks in 43 min (8 workers): 98 raw
  solutions, 44,474 deferred pairs; parallel SAT tail resolved them in
  ~33 min, adding 66 more: 164 raw at delta >= 10, zero contract violations
  after the engine fix.
- **Final: 209 raw verified solutions -> 31 canonical iso classes; canonical
  set EQUAL to the published census (`D3 RESULT: MATCH`); 77 min wall.**
- Engine fix regressions: leaking bucket chunk emits 0 with the fixed binary
  (1 with the old, min degree 9); full (16,71) stratum at true bound
  GLUE_DMIN=5: old and fixed byte-identical (21,316 raw graphs). Pre-fix
  binary: `outD/bench/glue_census2_c.pre_fix`.

## Artifacts

- `src/case_split_21_107.py` — orchestrator (proofs in docstring; every
  emitted graph independently re-verified before entering the raw file).
- `outD/case_split_21_107/` — isolated outputs: raw g6, canon, pub.canon.
- `outD/bench/` — the engine benchmarks above.
- outA shard-file incident: my earlier default-outdir launch truncated
  outA/tmp_21_107/s{0..7}.{g6,csv}; snapshot inspection proved all 20
  original shard files were 0 bytes (mtime 08-13 07:07:25), and the 16
  truncated files were restored byte-identical with original mtimes from
  the 15:03:06 APFS snapshot. Zero evidentiary loss.
