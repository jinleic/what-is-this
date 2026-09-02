# q=4 total-output-orbit-mass certification report

## Headline and verdict

**MACHINE-VERIFIED:** all 24 preregistered finite-blocklength rows for `q=4`, `n=5..10`, and `d in {1/2,1/5,1/10,1/20}` are `CERTIFIED`. Every conservative interval lies strictly inside both published Tavakoli-Nguyen-Bose endpoints: the certified lower endpoint beats `LB+`, and the certified upper endpoint beats `UB`. Every selected dual is the exact total-output-orbit-mass candidate `ba_total_orbit_mass`; uniform was evaluated over the same full input alphabet but never selected.

This campaign closes the deferred q=4 restart. It does not make an asymptotic-capacity claim and does not adjudicate `q=3,n=11`.

Campaign: `campaigns/2026-08-31T09:59:17Z_b6cd7315-caf3-4225-bb78-e4a81bf9a0f7_q4-total-output-orbit-mass/`.

Frozen production source SHA-256: `c9c3a84213d46da9a38fdf57a3761adddff950817091d8068e781630376f08b1`.

## Certified table

The interval endpoints are conservative outward decimals backed by exact binary rationals and 400-bit Arb enclosures. Widths are outward upper bounds, in bits per symbol. Every row has verdict `CERT_LOWER_BEATS_LBplus+CERT_UPPER_BEATS_UB`.

| n | d | certified interval | width upper | selected dual |
|---:|:---:|:---|---:|:---|
| 5 | 1/2 | `[0.6664806108007938, 0.6664806978587714]` | `8.705797740570955e-08` | `ba_total_orbit_mass` |
| 5 | 1/5 | `[1.3277006864199778, 1.3277006950433863]` | `8.623408079458642e-09` | `ba_total_orbit_mass` |
| 5 | 1/10 | `[1.6380465841890821, 1.638046590354154]` | `6.165071504436472e-09` | `ba_total_orbit_mass` |
| 5 | 1/20 | `[1.812336819003943, 1.8123368449964559]` | `2.5992512288507346e-08` | `ba_total_orbit_mass` |
| 6 | 1/2 | `[0.6358515734372481, 0.6358521206535938]` | `5.47216345443326e-07` | `ba_total_orbit_mass` |
| 6 | 1/5 | `[1.2979227045312098, 1.2979227333476206]` | `2.8816410510323828e-08` | `ba_total_orbit_mass` |
| 6 | 1/10 | `[1.6189070156070138, 1.6189070393689304]` | `2.3761916190006316e-08` | `ba_total_orbit_mass` |
| 6 | 1/20 | `[1.801553474543007, 1.8015535054230871]` | `3.088007948575566e-08` | `ba_total_orbit_mass` |
| 7 | 1/2 | `[0.611799221875365, 0.6118040734281766]` | `4.851552811291828e-06` | `ba_total_orbit_mass` |
| 7 | 1/5 | `[1.2736606165054785, 1.2736607523151153]` | `1.3580963632724146e-07` | `ba_total_orbit_mass` |
| 7 | 1/10 | `[1.602961519916888, 1.6029616108940126]` | `9.097712407753996e-08` | `ba_total_orbit_mass` |
| 7 | 1/20 | `[1.7924499270995096, 1.792450151213695]` | `2.241141850599111e-07` | `ba_total_orbit_mass` |
| 8 | 1/2 | `[0.5923193309337754, 0.5923424560693681]` | `2.3125135592425783e-05` | `ba_total_orbit_mass` |
| 8 | 1/5 | `[1.253452581576217, 1.2534531084585532]` | `5.268823357441099e-07` | `ba_total_orbit_mass` |
| 8 | 1/10 | `[1.589402062820447, 1.5894023890079079]` | `3.261874603539729e-07` | `ba_total_orbit_mass` |
| 8 | 1/20 | `[1.7846066171994794, 1.7846067975948448]` | `1.803953648481094e-07` | `ba_total_orbit_mass` |
| 9 | 1/2 | `[0.5761595996416342, 0.5761607360557939]` | `1.136414159464924e-06` | `ba_total_orbit_mass` |
| 9 | 1/5 | `[1.2363290100118065, 1.2363309197159202]` | `1.9097041134344336e-06` | `ba_total_orbit_mass` |
| 9 | 1/10 | `[1.5776876705095761, 1.5776899006080805]` | `2.2300985039430338e-06` | `ba_total_orbit_mass` |
| 9 | 1/20 | `[1.7777426950856503, 1.7777433069921658]` | `6.119065150642508e-07` | `ba_total_orbit_mass` |
| 10 | 1/2 | `[0.562495856308736, 0.5625006203492582]` | `4.764040521997813e-06` | `ba_total_orbit_mass` |
| 10 | 1/5 | `[1.2216151596360123, 1.221624099427772]` | `8.939791759484224e-06` | `ba_total_orbit_mass` |
| 10 | 1/10 | `[1.567438365675021, 1.5674413275691748]` | `2.9618941532341816e-06` | `ba_total_orbit_mass` |
| 10 | 1/20 | `[1.7716609170039264, 1.771663270194063]` | `2.3531901360821368e-06` | `ba_total_orbit_mass` |

The width range is `6.165071504436472e-09` at `(n,d)=(5,1/10)` through `2.3125135592425783e-05` at `(8,1/2)`. The smallest serialized positive Arb lower-margin lower endpoint is approximately `0.0007402311927961205436222872684` bits/symbol, at `(5,1/20)`; the smallest upper-margin lower endpoint is approximately `0.08766315500354434284029482356` bits/symbol, also at `(5,1/20)`. The strict inequalities themselves are certified by the archived Arb balls, not by these decimal summaries.

## Exact checks and resource evidence

**MACHINE-VERIFIED:** all 24 census anchors passed; all input and output orbit-size sums equal the full alphabet sizes; exact input and output distributions sum to one and pass every expanded group-generator equality; every primal orbit calculation overlaps a direct full-input calculation; and both dual candidates agree between representative and direct full-input maxima. Aggregate checks recorded in the rows are:

- `11,182,080` direct full-input dual evaluations across the two fixed candidates;
- `3,585,643,216` positive-`W`/positive-`D` conditional terms checked for those duals;
- `1,792,821,608` direct full-input primal conditional entries evaluated;
- `16,773,120` expanded input-distribution generator checks;
- `44,728,272` expanded output-distribution generator checks;
- zero primal or dual full-word/representative overlap failures.

The runner's immutable-resume validator returned `PASS` with exactly 24 rows after completion. Independent owner checks verified all 24 newline-inclusive JSONL line hashes and established that the CSV is an exact projection of the 24 canonical JSONL rows.

**COMPUTATIONAL-EVIDENCE:** the production run used one low-priority process with all five thread pools pinned to one. Campaign runtime was `4651.466 s` wall / `4649.094 s` CPU; the maximum row-recorded RSS was `1,441,169,408` bytes. No row hit a resource stop.

## Evidence and provenance

Canonical row evidence is `orbit_rows.jsonl`; `orbit_rows.line_checksums.json` binds each canonical newline-inclusive row. `TABLE_q4_conservative_intervals.csv` is a derived presentation table. `summary.json` is the terminal runner summary. `manifest.json` and `checksums.sha256` are preserved pre-production/interim records; additive `manifest_final.json` and `checksums_final.sha256` close the campaign without rewriting that history.

The post-acceptance hold note was appended after an instruction to stop edits; its bytes remain preserved. Subsequent owner perimeter corrections, source refreezes, static guards, and the released anchor are likewise append-only in `pre_statement.md`. The production result uses the final frozen source hash above. Import-generated `__pycache__/` files remain preserved but are excluded from the evidence checksum ledgers because they are neither frozen source nor result artifacts.

## Rule-7 scope sentence

The completed sweep is exactly `q=4` with `(n,d)` in `{5,6,7,8,9,10} x {1/2,1/5,1/10,1/20}`, using `G=S_4 x C_2` on symbol values and word reversal, exact total input-orbit mass snap `2^30`, exact total output-orbit mass snap `2^30`, fixed dual candidates `(ba_total_orbit_mass, uniform)`, direct full-alphabet KL comparison, full finite channel support, and Arb 400-bit outward rounding. **NOT swept:** `q=3,n=11` or any q=3 row; `q>=5`; `q=4,n>=11`; other deletion probabilities; `ba_word`; non-G-invariant input or output candidates; the falsified `S_n` position/type reduction; the cause of Tavakoli-Nguyen-Bose rounding; Morozov-Duman, Pinto-Ribeiro, and every `n -> infinity` or asymptotic-capacity claim. Every certified row is only a finite-n theorem.

## Next campaign

Named next campaign: **q=3,n=11 total-output-orbit-mass extension**. It must be a separate prospective campaign because `q=3,n=11` was explicitly excluded here. Its source, exact row box, resource cap, and stop rules must be frozen before any result computation.
