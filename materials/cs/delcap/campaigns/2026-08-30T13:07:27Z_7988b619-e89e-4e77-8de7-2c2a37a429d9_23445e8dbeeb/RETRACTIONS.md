# Retractions — Gate C (rule 5: recorded inline, never removed)

## R1 — "the LO-CVB rows disagree with the paper" — RETRACTED (cleared)

ORIGINAL CLAIM (mine, 2026-08-30, escalated to Main): 19 of 36 LO-CVB rows
"disagree" with Morozov-Duman's printed Table III, my certified value lying
below theirs by 5.0e-6 to 1.04e-5, and the cause is "last-digit rounding in a
printed table".

FALSIFIED IN TWO STAGES.

(a) The stated CAUSE was falsified by Main with exact decimal arithmetic. For
m=22,n=32 my certified 0.66389960 rounds to 0.66390 under BOTH round-half-up
and ceiling while the printed value is 0.66391; for m=5,n=1024 my certified
0.80277963 rounds to 0.80278 under both, printed 0.80279. The differences
1.040e-5 and 1.037e-5 exceed the 5dp ulp of 1.0e-5, and a printed-table
rounding artifact of my own UNSHIFTED value cannot exceed one ulp. Withdrawn.

(b) The 19-vs-17 PARTITION was then falsified by me. Sorting all 36 rows by
gap = printed - certified, my partition falls exactly at the sort boundary:
every row with gap < 5.04e-6 was labelled "reproduces", every row with
gap >= 5.04e-6 "strictly tighter". The split was a cut at the +-5e-6
comparison window I chose myself, and carries ZERO structural information. The
two classes are interleaved on every structural property: m in {5,22,23} on
both sides, full-support Lambda 11/19 vs 6/17, n from 1 to infinity on both
sides, |Lambda| from 5 to 24 on both sides. Correlations of gap are weak and
consistent with noise at n=36: n +0.40, m*n +0.28, is_full_set +0.20, m -0.08,
|Lambda| -0.03.

(c) A third hypothesis, Main's, was also falsified: that the gaps show a
SYSTEMATIC ADDITIVE finite-n term, inferred from two quoted gaps sitting
3.0e-8 apart. The gaps are a CONTINUUM — 35 distinct values among 36 rows at
1e-8 resolution — and those two are simply the top two of it (third:
1.0238e-5), a top-of-distribution proximity artifact.

(d) CORRECTION I VOLUNTEERED AGAINST MY OWN REPORT: the direction split is
36-0, not 19-0. In every one of the 36 rows the printed value exceeds the
certified value; minimum gap +7.1e-7.

WHAT REPLACED IT: a reproduction with a named printing convention. See
MD_TABLE_III_FINDING.md in this campaign. Delta band under ceiling printing
(3.990e-7, 7.10008e-7], determined by the min and max gaps rather than fitted,
non-empty precisely because the spread 9.689365e-6 is less than one 5dp ulp.
No published claim is contradicted.

## R2 — no other Gate C claim has been retracted.

The d=1/2 sandwich improvement, the paper-grid recomputation rows, the
Pinto-Ribeiro extension rows and the delta in {1/20,1/2,4/5} extension rows
stand as first reported.
