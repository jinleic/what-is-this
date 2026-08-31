# ADDENDUM — feasibility audit of the three improving boxes (Main question
# (a)-(c)), run 2026-08-30T12:1xZ; supersedes the "box results" verdict
# wording wherever they conflict.

## Answer

(a) Feasibility status of the three certified-below-B0 boxes:

- region_prop+ (centre rp = p*_rp + 1e-7 each of the 3 coords):
  **INFEASIBLE** as an exact program point. The region_prop block is
  pinned by THREE linear equality rows (sum(rp)=1, [0,1,-1]·rp=0,
  [1,-1,0]·rp=0 — rows 12/85/86 of the recorded ledger; the asymmetry
  rows exist and are correctly recorded; my earlier "no -1 rows" scan
  was a bad filter). Uniform shift +1e-7 keeps Y−Z = 0 and X−Y = 0
  exactly, but the exact-dyadic sum of the shifted triple is
  1.0000003 (residual +3.000000000086e-07 against b=1). Measured in
  exact dyadic arithmetic (IC.const endpoints, no ball radii),
  not float eyeballing.

- glob_dist_0+ / glob_dist_1+ / glob_dist_2+ (all 7 shape coords of one
  region +1e-7): **INFEASIBLE**. Each moves the exact simplex-sum
  residual of that region's dist block by +7.0e-07 (7 coordinates ×
  1e-7), and it also breaks the dist↔dist_max marginal-equality rows
  (interval path equality residual at the 1e-7 scale) while dist_max is
  held fixed. Same exact-dyadic method.

So Main's suspicion was correct, and Main's instinct to ask was the
load-bearing move: (a) => all three boxes are objective evaluations at
inadmissible points. They are NOT feasible points; they bound nothing
about omega.

(b) n/a (no box is feasible).

(c) Downgrade performed: the README Part 2 bullet and the campaign box
record now say "certified objective evaluations at exact-infeasible
centres bound nothing about omega" instead of "certify strictly below
B0". Rule-5 retraction is inline (this file, plus README edit marked
RETRACTED in place, not deleted).

## Feasibility repair — cheap enough to settle in THIS campaign?

Partially, and the repair analysis changes the picture materially:

1. EXACT/rs-float level: for rp, ALL THREE rows (sum=1, Y−Z=0, X−Y=0)
   pin rp to a uniform triple r,r,r with 3r = 1 exactly over the
   rationals. Exactly-feasible float triples exist (r = 1/3 and
   r = nextafter(1/3,↑) both give r+r+r == 1.0 exactly and asym 0
   exactly) but the released p*_rp is NOT one of them: its stored triple
   (0.3333333333333307, 0.33333333333333465, 0.33333333333333465)
   violates X−Y by −3.941291737419306e-15 (sum exactly 1.0, Y−Z exactly
   0) — a rational inexactness the frozen eps_abs (1.94e-11) already
   absorbs. A MIXED ±1e-7 repair (x+2e-7, y−1e-7, z−1e-7) keeps the
   exact sum at 1.0 but violates X−Y by 3.0e-7 — so within this box
   there is NO non-uniform rp move that is exactly feasible; the rp
   coordinate is effectively Frozen by the three rows at ±(dyadic tail)
   scale. Feasible rp neighbourhood: |Δr| <= 5.5e-17-class only.
   Those near-micro moves were NOT probed and cannot matter: they are
   10 orders below the 1.94e-11 defect absorption.

2. For glob_dist: a simplex-repaired probe (e.g. 6 coords +1.16667e-7,
   1 coord −1e-7 to first order, exact-dyadic-verified sum) IS cheap
   (~4 s per box on this machine, same script path with a repaired
   offset vector). It was not part of the pre-registered probe family,
   so per pre-statement §8 it must go through a pre-statement addendum
   BEFORE running: this addendum pre-registers exactly that repaired
   family (which coordinates move, by how much, exact-dyadic sum
   residual acceptance threshold 0 in double float64, asym rows exact 0,
   dist-vs-dist_max equality residual <= 1.1e-9 — wait, these equality
   residuals are NONLINEAR in the certified standard; repaired probes
   must also be re-verified against the ceq set at the repaired point
   before any bound claim).

3. THE DEEPER POINT (why the three probes still told us something, and
   what we can honestly claim): the probes measured the SENSITIVITY of
   the certified endpoint to moving OUT of the feasible set. A
   certified local-optimality statement over the box requires, per
   direction, BOTH (i) an enclosure at feasible points in that
   direction and (ii) no descent. The current campaign establishes
   (ii)-style information only at infeasible points, so the honest
   gate-C statement is the WEAKER one:
     "Within the pre-registered box, at the only evaluated feasible
      point (p* itself, with its 1.94e-11-absorbed rational
      inexactness), the certified cap stands at 2.3715538358544612;
      at fourteen structured exact-infeasible centres the certified
      objective evaluation drops by up to 8.15e-7 — an upper-
      sensitivity statement, NOT a bound on omega."
   The first clause is a rigorous local statement about p*; the second
   is sensitivity information, not optimality.

## Not done here (would need Main approval + fresh pre-registration)
- simplex-repaired probe family NOT run: it extends the pre-registered
  probe family, which needs Main approval + fresh pre-registration
  first. (Feasibility-repair feasibility: ~1 min compute; deferred.)
- Any claim that the certified enclosure of the VXXZ24 rung improved to
  2.3715530208533490 — RETRACTED: the three centres are
  exact-infeasible, so no bound improved. The certified bound record
  for the VXXZ24 rung remains 2.3715538358544617 (frozen) with this
  campaign's tighter re-evaluation at 2.3715538358544612 (same
  quantity, tighter outward rounding).

## DECISION 2 EXECUTION — ParamManager ledger patch (2026-08-30T12:3xZ)

Bug 1 (b-replacement): `add_lincon_eq` REPLACED `_lin_beq` on every call
while `linear_violations()` zipped all 3174 rows against that single
trailing b. Fixed in src/ ONLY (`src/vxxz24_float.py`): `_lin_beq_list`
extended one-b-per-row, kept `_lin_beq` write for frozen-script reads.
The frozen campaign scripts are untouched (immutable snapshots); their
"linear violations 0.0" claims are RETAINED but their provenance now
rests on the independent exact-dyadic per-row replay at p* (which gives
0.0 for the old corrupt-as-recorded rows? no — see bug 2), not on the
accessor.

Bug 2 (row expansion, found DURING re-validation): `_lc_rows` read
`coeff[r, c]` with r = output-row index, silently corrupting EVERY
matrix-coefficient row set (all j2m "share-marginals" equalities, part
and global). Symptom used to catch it: with bug 1 fixed, the accessor
suddenly showed 1.478849e-04 at p* — the corrupt recorded rows are not
satisfied by the released vector, while the DIRECT marginal check over
all 63 Parts is <= 3.5258684860650646e-11 clean. Fix: MATLAB semantics
verbatim (ParamManager.m:142-152 — `lin_Aeq(rows, group) = coeff'`),
i.e. coeff is (n_group x nrow), transposed into per-output-row blocks.

Decision 2(b) re-validation, both rungs:
- rung 2 (this fix applies): corrected accessor = independent per-row
  replay = 3.525868e-11, agreement exact; vs the direct structural
  recomputation max |A d - A dmax| over all 63 Parts x 3 margins =
  3.5258684860650646e-11 — three instruments, one number. The stage-(a)
  headline numbers are UNCHANGED (c 4.6349590904e-12, ceq 3.5953524880e-11,
  value/target exact); the stage_b_rung2 anchor re-ran at
  2.3715538358544617 (unchanged; M-verified in this session's earlier
  runs) and the corrected `max linear violation` printout at the entry
  point now reads 3.525868e-11 inside the 1.1e-9 refine tolerance.
- rung 1 (alman25_float.py): the transcription records NO linear ledger
  (no AddLinearConstraintEq calls — its ParamManager has no such API), so
  the accessor n/a; rung-1's own frozen stage-(a) numbers came from the
  c/ceq violation path, not from this accessor, and are unaffected.

Transferable lesson (Main's request, one sentence): third instance this
session of a validation validating a component it shares with the thing
checked — kg's finite-difference check of a wrong function against its
own derivative, mm3's three solvers on one modelling premise, and here an
accessor comparing 3174 rows against one b; the fix in every case is an
instrument built from different parts than the thing under test.

Post-fix rp+ / glob_dist+ feasibility conclusions UNCHANGED: the row
data corruption was in the recorded coefficients, not in the published
parameter data; the feasibility audit used exact-dyadic per-row replay
of the CORRECTED constraint semantics (simplex sums and asymmetry rows
are group-uniform rows that were never corrupted) and stands.
