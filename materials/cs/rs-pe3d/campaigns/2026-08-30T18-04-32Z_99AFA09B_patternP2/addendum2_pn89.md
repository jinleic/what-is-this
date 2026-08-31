# ADDENDUM 2 — PN8/PN9 completed (missed by the crashed first runner), 2026-08-30T~00:10Z

The first runner crashed at PN4; PN8 (11,(2,2,5)) and PN9 (17,(2,2,8)) had
NOT run when the freeze was cut. Both are now complete (scratch/
p2_pn89_results.json; witnesses appended in witnesses_addendum2.json):

- PN8 (11,(2,2,5)): rho^window = 1/2 (wt=2, delta=4, closure proved),
  30 candidates. P predicted "not 3/5" -> HIT. First genuinely NEW sub-1
  value produced by the battery.
- PN9 (17,(2,2,8)): rho^window = 1/2 (wt=2, delta=4, closure proved),
  48 candidates. P predicted "not 3/5" -> HIT.

Battery final: 19 rows = 9 PP (9 hits) + 7 PN landed (7 hits; PN8, PN9 =
1/2, the rest 1/1) + 3 empty-window (PN4, PN5, PN10). Zero P misses.

H-COP.b RE-ADJUDICATION with the two new rows: distinct sub-1 values among
landed non-coprime rows are now {3/5, 1/2} = TWO (was one). Needed >= 3.
The falsification of H-COP.b STANDS (2 < 3), but the falsifying gap
narrowed from 1 vs 3 to 2 vs 3 - the frozen table's five-value richness
still did not reproduce. H-COP.a unchanged (PN8/PN9 add no coprime rows).

H-GATE unchanged (both new rows have min(s)=2 <= 3, non-empty windows).

pn8_pn9 witnesses verified in exact F_q via src.verify (inside this
addendum generation); checksums for the two new files below.
