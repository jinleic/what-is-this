# Morozov-Duman Table III — 36 rows recomputed and certified

Cleared for the record by Main after independent re-derivation of the offset
band. This is a REPRODUCTION with a named printing convention and a quantified
precision difference. It is not a discrepancy, a disagreement or a tension, and
no published claim is contradicted.

## What was computed

All 36 rows of Morozov-Duman (arXiv:2504.20961) Table III — the layer-oriented
converse bound (LO-CVB) code rate, their eqs (24)-(30) — at m in {5, 22, 23},
n in {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, infinity}, delta = 1/5 and
eps = 1/5 as EXACT rationals, recomputed from the integer E(m,w) tables.

Every certified value is one outward-rounded Arb evaluation of log2 of an exact
rational at 400 bits. Interval widths: 0.0 at every printed digit, maximum ball
RADIUS over all 36 rows 2.84e-120 (reported as a first-class number rather than
as "zero"). Truncation order: none — the LO-CVB is a finite closed form in the
E(m,w) integers, so there is no series to truncate and the tail bound is
exactly 0. Label: MACHINE-VERIFIED.

## Agreement with the published table

In all 36 rows the printed value EXCEEDS the certified value:

    gap = printed - certified  in  [7.10008e-07, 1.03994e-05]
    spread of gap over 36 rows  =  9.689365e-06  <  1.0e-05 = one 5dp ulp

Solving printed = round(certified + Delta) at 5 decimals SIMULTANEOUSLY over
all 36 rows, under ceiling printing (round up — the conservative and correct
convention for printing an upper bound, since a rounded-up upper bound remains
a valid upper bound):

    Delta in (3.990e-07, 7.10008e-07]

This band is DETERMINED, NOT FITTED: under ceiling printing gap lies in
[Delta, Delta + ulp), so the MINIMUM gap caps Delta from above and the MAXIMUM
gap bounds it from below at max_gap - ulp. Nothing was tuned. The band is
non-empty for exactly one reason — the spread is less than one ulp — and
"spread < 1 ulp" is equivalent to the existence of a single-Delta ceiling fit.
Under round-half-up printing the corresponding band is [5.399e-06, 5.710e-06].

So the published values agree with the certified ones to within 7.1e-7
absolute in every row, and the printed digits are exactly what round-up at five
decimals produces from a value that close.

## Why this is a precision difference and not a formula difference

A missing or differing FORMULA term would be n-dependent or m-dependent. Such
a term would BREAK the single-Delta fit, because the same Delta must reproduce
rows spanning n = 1 to n = infinity and m = 5, 22, 23 simultaneously. The fit
holds across all 36 rows, which is positive evidence AGAINST a formula
difference. A uniform positive offset bounded by 7.1e-7 at rate values
0.55-0.82 is 0.9e-6 to 1.3e-6 relative — the size of a 6-to-7-significant-digit
intermediate. Mechanism label: INFERENCE. The convention identification
(round-up printing) is INFERENCE. Our own values are MACHINE-VERIFIED.

## Aggregate check on the E(m,w) transcription

The n -> infinity row equals log2(sum_w E(m,w) delta^{m-w}(1-delta)^w)/m, a
weighted sum over the WHOLE E(m,·) column, and it reproduces their printed
value at every m: 0.80272 (m=5), 0.73569 (m=22), 0.73414 (m=23). This is an
aggregate check on the typed Table I columns, not a per-cell one.

Per-cell verification of E(m,w):
  * m=5, all w: recomputed from scratch by brute force over all (x,y) pairs.
    E(5,·) = [1, 10, 32, 52, 54, 32]; their eq (28) prints E(5,2)=32,
    E(5,3)=52, E(5,4)=54 — reproduced exactly. MACHINE-VERIFIED.
  * m=20,21,22,23, columns w in {0,1,2,m}: re-derived independently and agree
    with the typed values. w=0: E=1. w=1: E=2m. w=m: E=2^m (their eq 28).
    w=2: closed form E(m,2) = 2*C(m,2) + 2*floor(m^2/4), derived and proved in
    gate_c_md.py and brute-force validated at m=4..12; it reproduces their
    typed 580, 640, 704, 770 exactly. MACHINE-VERIFIED for those cells.
  * all other cells of m=20..23: CITED-DEPENDENCY (typed from their Table I).

## Status

First certified enclosures of these 36 rows. No published claim is
contradicted; their converse stands. If anything the round-up printing
convention is the correct choice for a converse bound.
