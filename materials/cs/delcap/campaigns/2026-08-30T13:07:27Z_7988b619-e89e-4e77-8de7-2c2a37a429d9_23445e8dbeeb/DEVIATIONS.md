# Pre-registration deviations — Gate C

`pre_statement.md` was committed BEFORE any Gate C computation (its own header
records this, and a byte-identical copy is frozen in this campaign). It is NOT
edited retroactively. Every difference between the committed grid and what was
actually run is listed here, with cause. Rule 7: the swept set below is the
exact claim; nothing outside it is claimed.

## (A) Tavakoli-Nguyen-Bose sandwich rows

Committed: q in {2,3}, n in {3,5,10}, d in {1/20, 1/10, 1/5} (18 points), plus
extension to q=4, to q=3 n in {7,8}, and to d=1/2.

RAN, exactly:
  d = 1/2 (the ticket's target): (q,n) = (2,2) (2,3) (2,4) (2,5) (2,6) (2,7)
      (2,8) (3,2) (3,3) (3,4) (3,5) (4,3) (4,4)                = 13 rows
  paper grid: (q,n) = (2,3) (2,5) (2,10) (3,3) (3,5), each at
      d in {1/20, 1/10, 1/5}                                    = 15 rows
  28 distinct rows. Two rows — (3,3,1/2) and (3,5,1/2) — were additionally
  re-run in a SEPARATE process as a reproducibility check and agreed to every
  printed digit (0.594406224780 / 0.594406260821 and 0.512042 both times).

NOT RUN, with cause:
  * (q,n) = (3,10) at d in {1/20,1/10,1/5} — 3 committed points dropped.
    Cause: the exact channel is 3^10 = 59049 inputs by sum_k 3^k = 88573
    outputs, a dense exact-integer matrix of 5.2e9 entries — beyond this
    workstation. The pre-statement already flagged this point as at-risk.
  * extension (q,n) = (3,7), (3,8) at d=1/2 — not run. Cause: 3^7 x 3280 and
    3^8 x 9841 exact-rational certificate loops exceed the wall-clock budget.
    The q=3 ladder stops at n=5; the q=2 ladder at n=8 (n=10 on the paper grid).

METHOD DEVIATION (material, stated plainly):
  The pre-statement fixes the locating step as "mpmath 150-dps BA". These rows
  were produced with a float64 numpy BA locating step instead.
  Cause: at q=3, n=5 the mpmath-150dps locator cost ~290 s per row and, under
  heavy deletion (d=1/20), its iterates underflowed and collapsed onto a
  near-degenerate p, giving a valid but very loose dual.
  WHY THIS CANNOT AFFECT VALIDITY: the locator only proposes candidates. The
  certified lower bound is the exact mutual information I(p*) of an exact
  rational p* = m/M — a true MI, hence <= C for ANY p*. The certified upper
  bound is max_x KL(W(.|x) || D') for an exact rational full-support D' —
  hence >= C for ANY such D'. Both are one-shot outward-rounded Arb
  evaluations of logs of exact rationals; there is no interval iteration. A
  worse locator costs interval WIDTH and nothing else. Achieved widths are
  2.9e-10 to 6.7e-6 bits/symbol, tighter than the mpmath path achieved.

## (B) Pinto-Ribeiro C_{n,k} rows

Committed: (n,k) with k <= n <= 12, listed explicitly up to n=8.
RAN: n = 6,7,8,9 for all 1 <= k <= n — 30 rows.
NOT RUN: n in {10,11,12}. Cause: wall-clock.
Pinto-Ribeiro publish C_{n,k} only at n in {29,31} (GPU), out of certified
reach here — stated in the pre-statement as a structural limit, not discovered
late. Our n=6..9 rows therefore have NO published analogue: they are new
certified values, not a reproduction.

BECAUSE there is no published value to check them against, each row is
produced TWICE by different machinery and the two are compared:
  route 1 (gate_c_pr.py): float64 BA locator -> exact rational snap (2^30 in,
     2^40 out) -> gate_c_dhalf.cert_primal / cert_dual, Arb at 400 bits
  route 2 (delcap_cert, frozen): mpmath 150-dps BA locator -> exact rational
     snap (2^24) -> certify_ab, Arb at 300 bits
The routes share only the exact integer channel matrix subseq_matrix, which is
itself anchored by the frozen gate-B table. Both intervals must contain the
same true C_{n,k}, so the test is that they OVERLAP; a disjoint pair would mean
at least one certificate is wrong. RESULT: 22 of 30 rows have both routes
(route 2 was stopped at n=9,k=1 for cost — an mpmath row at n=9 costs
~10-80 min against ~1 s for route 1); all 22 OVERLAP, 0 disjoint. The
remaining 8 rows (n=9, k=2..9) are route-1 only and are labelled as such in
TABLE_pr_cnk_rows.csv (verdict FASTPATH_ONLY_no_mpmath_row).

## (C) Morozov-Duman LO-CVB rows

Committed: m in {5, 10, 15, 23}, E(m,w) recomputed from scratch for m <= 10,
anchored on their eq (28) values E(5,2)=32, E(5,3)=52, E(5,4)=54.
RAN: m in {5, 22, 23}.
DEVIATION: m in {10,15} dropped, m=22 added. Cause: their Table III prints
columns for m = 5, 22 and 23 only, so m=10 and m=15 have no published value to
compare against while m=22 does. Choosing the comparable columns is a
deviation from the committed set and is recorded as one.
NOT SWEPT (rule 7): the Lambda search is EXHAUSTIVE over all 2^6 = 64 subsets
at m=5. At m=22 and m=23 all 2^23 and 2^24 subsets are enumerated and ranked
in float64, and only the top 24 per row are certified in Arb. Since EVERY
Lambda yields a valid converse bound, incomplete certification can only make
the reported bound larger (weaker), never invalid. The reported value is the
min over CERTIFIED candidates — NOT a proven global optimum over Lambda.
