# PRE-STATEMENT ADDENDUM — P2 battery — Main min(s) prediction, committed BEFORE the remaining windows land

Committed 2026-08-30T22:59Z (UTC) by agent RsPe3dPatternP. Census state at
commit: PN10 (281,(5,7,8)) count-census is RUNNING at 1,800,000 of 3,658,900 supports (checkpoint
scratch/p2_pn10_census.json: done=1,800,000, nonzero=0) and has NOT finished; PN5, PN6, PN7 outcomes ARE already known
(PN5 empty, PN6 nonempty-42, PN7 nonempty-28). PN4 empty is known.

Main hypothesis H-GATE (Main, 2026-08-30, 0-for-12 record, input to be
tested): on this battery + frozen table, the emptiness gate of the
weight<=3 window is exactly min(s) <= 3 vs min(s) >= 4:
  - every instance with min(s) >= 4 has an EMPTY weight<=3 window;
  - every instance with some order <= 3 has a NON-EMPTY window.

Prospective part (not yet observed when this file was written): PN10 =
(281,(5,7,8)), min(s)=5. H-GATE predicts PN10 empty.
Retro-fitted part (already observed, labeled): PN5 (61,(4,5,6)) min=4 ->
empty (observed); PN6 (43,(2,6,7)) min=2 -> non-empty (observed 42 cands,
all wt=2); PN7 (337,(3,4,7)) min=3 -> non-empty (observed 28 cands, all
wt=3). PN4 (421,(4,5,7)) min=4 -> empty (observed).

Adjudication fixed now: H-GATE is CONFIRMED-in-window iff PN10 lands empty;
it is FALSIFIED iff PN10 lands non-empty (a single counterexample kills it,
same standard as P). The 14-row description behind H-GATE is labeled
IN-SAMPLE; only the PN10 read is prospective. Main caveats inherited: the
only-if direction rests on one prior empty point plus whatever this census
adds; no mechanism is derived from the dual product code here (OPEN).

This addendum file is written to scratch/ BEFORE the PN10 census outcome is
read (census at 49%, nonzero=0 so far). the outcome may be zero-non-
empty == simple; report either way.
