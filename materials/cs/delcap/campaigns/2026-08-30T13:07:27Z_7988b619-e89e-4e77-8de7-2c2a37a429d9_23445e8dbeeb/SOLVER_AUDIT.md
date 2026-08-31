# Solver / tolerance audit — Gate C

Prompted by a repo-wide warning from Main: a scipy HiGHS call elsewhere in the
repo returned a point violating its own bounds by 100% because the solver's
default primal feasibility tolerance (1e-7) EQUALLED the box radius being
enforced — an instrument that cannot resolve the quantity it is asked about.

## Finding: no solver anywhere in the Gate C certified path

grep for scipy / linprog / highs / optimize across gate_c_dhalf.py,
gate_c_md.py, gate_c_run.py and delcap_cert.py returns 0 hits. There is no
LP/QP/MILP in the certified path. Float appears in exactly two places, both
non-load-bearing: (a) a float64 numpy BA locating step, (b) float64 ranking of
Lambda subsets in the LO-CVB search. Both feed an exact-rational snap, after
which every certified number is an outward-rounded Arb evaluation of logs of
exact rationals at 400 bits.

Tolerance versus magnitude: certified per-symbol widths run 2.9e-10 to 6.7e-6.
No float tolerance is ever compared against anything at that scale. The
sandwich comparisons are Arb-endpoint inequalities (certified lower endpoint vs
the exact-Arb bound's upper endpoint) with margins of 0.02-0.13 bits/symbol —
four to six orders of magnitude above any interval width.

## A live invalid-certificate branch, found before it bit

The dual certificate is max_x KL(W(.|x) || D') for an exact rational D'. If
D'(y) = 0 while W(y|x) > 0 for some x, that KL is +infinity. Silently skipping
such a term would yield a too-SMALL "upper bound" — an invalid certificate that
looks perfectly clean.

This is not hypothetical. VERIFIED at q=2, n=10, d=1/20:
  * output alphabet 2047; BA output marginal minimum mass 9.765625e-14
  * snap resolution 2^-40 = 9.094947e-13, i.e. 9.3x COARSER than the mass it
    must represent
  * exactly 1 of 2047 entries snaps to zero, and that output IS reachable
    (W(y|x) > 0 for some x)
  * cert_dual on the raw snap returns +inf (no claim) — NOT a finite too-small
    number. The unsafe branch does not exist in the code.
  * with the +1 support bump: 8.529977690111 (valid, tight)
  * with uniform D': 9.672330910401 (valid, looser)
  * the row takes the min over candidates, so the tight valid one is used.

Same failure shape as the HiGHS trap — a resolution coarser than the quantity
being resolved — caught pre-emptively rather than after a retraction.
