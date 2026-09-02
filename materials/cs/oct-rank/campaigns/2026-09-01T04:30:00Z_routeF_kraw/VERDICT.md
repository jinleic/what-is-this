# VERDICT.md — Route F / af-triple13-existence campaign (OctRankRouteF)

Campaign: `campaigns/2026-09-01T04:30:00Z_routeF_kraw/`. Pre-statement
committed BEFORE any compute (commit `739e63f`); refinement log R1-R7 also
committed before the certification runs ( Regulations-5 corrections: exact
trilinear remainder W_e, P-pos plant disambiguation/1e-8 scale, sound
quarantine criteria, tau equation-count correction 48-not-64). Instrument
`rf_krawczyk.py`, replay `rf_replay.py`, outputs frozen with SHA256 in
`checksums.sha256`. All timing numbers: in-process `time.process_time()`
(rule 17e); no ps-derived numbers anywhere.

## Headline

**FAILURE TO CERTIFY — the Route-F rank-13 Krawczyk attempt does not close.
rank((L_1, L_i, L_j)) stays OPEN in {13, 14}. No rank-13 witness was
produced; no impossibility argument was found. By the pre-registered
vocabulary (and Main's condition 4): a numerical rank-13 candidate is not a
rank-13 result; the absence of a 13-witness is NOT evidence for 14, and the
absence of an impossibility argument is NOT evidence for 13. The published
window 18 <= R_R(T_O) <= 25 is untouched; no published work is refuted.**

## What was decided BEFORE the run (precompute artifacts)

1. `pre_statement.md` (commit 739e63f): the handoff's "20-variable
   symmetric-real" instrument was diagnosed as existence-incapable (a
   gradient-system Krawczyk certifies critical points, not zeros; interval
   arithmetic cannot turn a 1e-6 residual into an exact zero) — accepted by
   Main. Substituted instrument: SQUARE-SLICE Krawczyk on the rank-13 CP
   system — 192 equations (3*8*8), 247 unknowns, 55 coordinates frozen at
   exact dyadic values, 247 - 55 = 192 = balanced. Margin algebra = the
   gate-A margin algebra (exact fmpq throughout, outward-arb final
   comparisons, explicit midpoint/radius control; arb(lo,hi) is
   midpoint/RADIUS — never an interval, and no radius was formed by float
   endpoint subtraction anywhere).
2. `REFINEMENT_LOG.md` (pre-compute): R1 exact remainder W_e (the Taylor
   remainder of the trilinear system is quadratic+cubic with an exact
   absolute-value bound — unit-tested); R2/R3 P-pos plant construction
   (S-only perturbation, scale 1e-8 registered before compute); R4 sound
   quarantine criteria (false certificate = containment on the globally
   infeasible tau-6 system); R5 pass-window honesty; R6 tau count 48 not
   64 (found at dry-run, before any ladder compute); R7 tau-6 budget note.
3. Fixed seeds: the two frozen candidates
   (rank13_candidate.npz rel 1.3252060344282049e-3 seed 2; 
   rank13_extension_candidate.npz rel 9.99874493968644e-5 seed 2), both
   SHA256-verified against the Route-AF frozen ledger BEFORE use.
4. Fixed 15-rung ladder 1e-2 .. 1e-12 (Main's 1e-2..1e-8 verbatim, extended
   downward only); fixed polish budget (5 rounds x 5000 nfev TRF-exact per
   pipeline); fixed slice rule (QRCP with declared fallback); fixed
   selection and tie rules; NO other seeds, NO box widening, NO re-slicing.

## Startup anchors (all re-verified in every run; run5 = freeze of record)

- [ANCHOR-TBL] 3-way octonion table byte-exactness: cd == gateA ==
  upstream (integer lists), e0 unit, norm multiplicativity at 20 fixed
  points. PASS [MACHINE-VERIFIED].
- [ANCHOR-BLK] the corrected certification-slice equivalence chain, all
  three slices entrywise exact: L-matrix of (L_1,L_i,L_j) IS
  blockdiag(Lq_p, S2 Lq_p S2) with S2 = diag(1,1,1,-1) [this is the exact
  content of the frozen s3_routeF2/F3 facts]; and Sbig = diag(...,-1 on e7)
  conjugates blockdiag(Lq_p, Lq_p) to blockdiag(Lq_p, S2 Lq_p S2). Both are
  diagonal(similitude) identities, Machine-verified — therefore
  rank(blockdiag(tau,tau)) = rank((L_1,L_i,L_j)): a certificate for the
  target certifies the Route-F triple. [MACHINE-VERIFIED entrywise]
  (Rule-5 note: the first anchor attempt asserted entrywise
  T L T == blockdiag(tau,tau) which is FALSE — the complement block carries
  the S2 conjugation; caught by the instrument's own startup anchor at run2,
  fixed as above, documented in RUN_NOTES.md.)
- Seed-reproduction anchors on the CORRECT (blockdiag) target:
  CP rel reproduces 1.3252060344282049e-3 bit-exactly (diff 0.0), EXT rel
  9.9987449458480189e-5 (diff 6.16e-14 < 1e-12 trigger). PASS
  [MACHINE-VERIFIED / float-reproduction diagnostic].

## Certification outcome (run5 of record; run4 bit-identical on all
## mathematical lines — determinism check: 95/95 lines equal)

MAIN system (T_F-as-blockdiag target, rank 13, seeds P0..P5 polished):
- Polish: CP pipeline 1.325e-3 -> 8.816619e-06 rel (5x5000 nfev); EXT
  pipeline -> 5.694968e-05; selection rule picks CP. max|factor| = 32.4
  (within the 1e4 bound).
- Slice: QRCP, rank(J) = 192/247, min float pivot 7.9e-7 (>= 1e-9:
  Slice-1 fires; float diagnostic only — the exact facts below carry the
  load).
- Exact facts: Y = J_S^-1 exact inverse verified (Y*J_S = I in fmpq),
  30.2 s CPU; exact rank(J_F) = 55/55 (FULL F-block transversality — no
  slice-handicap).
- Ladder (all 15 rungs reported): containment 0/192 at EVERY rung; certified
  exclusion (K disjoint from the box) grows 44/192 coords at rho=1e-4,
  190/192 at 1e-5, 192/192 (complete, certified no-root-in-box over the
  named box) for all rho <= 1e-6 through 1e-12.
- The two straddling coordinates at rho=1e-5 (exact re-computation, frozen
  factors): S-coord 38 = original column 240 = c-block (7,6); S-coord
  50 = original column 99 = b-block (4,8) — i.e. the c7-row-6 and b4-col-8
  factor coordinates are the last two whose bounce/remainder still straddle.
- Verdict: **FAILURE TO CERTIFY**. The polished residual 8.8e-6 relative
  (7.9 max-dyn) is ~4 orders too large for the certified window: the
  containment needs ||Y g0|| < rho while the quadratic margin needs
  rho <~ sigma_min(J_S)/(132 * factor-scale) ~ 1.9e-5-ish (float est) — the
  window is EMPTY at this residual. No rung passes; exclusion is complete
  and certified at 1e-6..1e-12, i.e. the instrument CERTIFIES there is no
  root in the (shrinking) boxes around the polished candidate — the
  candidate is not near any on-slice exact rank-13 decomposition, at the
  level of certified statements. This is NOT a rank-13-impossibility
  (off-slice roots are not excluded; nothing global is claimed).

## Controls (all through the identical instrument path; behavior contract)

- C-TAU-7 (known-certifiable case): CONTAINMENT at rho=1e-3 (48/48 coords,
  arb-out re-verified strict on all coords). The tau_r7 frozen certificate
  is independently re-certified as an EXISTENCE certificate by this
  instrument in 0.27 s of exact arithmetic. [MACHINE-VERIFIED]
- C-TAU-6 (known-impossible case, the quarantine test): NO containment at
  any rung; certified exclusion fires: 40/48 coords at rho=1e-2, complete
  48/48 at rho <= 1e-3 — a MACHINE-VERIFIED no-root-in-named-box statement
  for the tau-6 system (whose exact rank-7 status is the frozen
  machine-verified fact). The instrument correctly REJECTS the globally
  infeasible system. [MACHINE-VERIFIED]
- P-POS (positive synthetic rank-13 control): CONTAINMENT at rho=1e-4
  (192/192, arb strict) — the plant root (exact dyadic factors, S-only
  1e-8 perturbation per R2/R3) is found and certified end-to-end.
  [MACHINE-VERIFIED]
- P-WRNG (counterfactual, near-miss corrupted target): the T_F-certified
  structure does NOT contain on the corrupted table (entry [0,0,0]=2):
  containment 0/192 everywhere; certified exclusion from rho=1e-3 (191/192)
  and complete at <= 1e-4. The instrument does not certify proximity into
  existence. [MACHINE-VERIFIED]
- Quarantine criterion (R4): NOT triggered (no containment ever fired on
  the globally infeasible tau-6 system). Behavior contract: PASS
  (exit 0; `behavior_ok: true` in the frozen JSON).

## The instrument's calibrated discrimination (rule 14 summary)

  tau-6 (no root anywhere):  rejected, certified exclusion complete.
  MAIN polished (rel 8.8e-6): rejected, certified exclusion complete at
                              1e-6..1e-12; no root within certified boxes.
  tau-7 (true root nearby):   ACCEPTED at rho=1e-3.
  P-POS (true root nearby):   ACCEPTED at rho=1e-4.
  P-WRNG (near-miss wrong):   rejected with complete exclusion.
The instrument separates every case in the expected direction; the
lower-detection edge (exact root present but window does not open) is
exactly the failure branch the two rejections share, which is why the MAIN
outcome can only be reported as failure-to-certify, never as impossibility.

## Upward side (assignment item 4) — honest absence

No exact argument that rank 13 is IMPOSSIBLE for this slice class was found
(U1: the frozen theorem audit's negatives stand — no covering
additivity/multiplicativity theorem, Strassen route closed at 12 < 13;
U2: the pre-committed necessary-condition probe produced no obstruction and
no proof either way; U3: Koiran's exact characterization noted as the
natural full-problem system, out of budget). **Absence of an impossibility
argument is not evidence for 13**; equally, the certified local exclusions
are not evidence for 14. Both absences are recorded with exactly this mean-
ing: the question stays OPEN in {13, 14} with the same certified interval
[13, 14] as before this campaign.

## Defects disclosed (rule 5)

1. Target-convention error in the FIRST instrument run (targeted raw slices
   instead of the conjugated blockdiag target — seed anchors failed
   8.2e-1 vs 1.3e-3): caught by the pre-registered seed-anchor check, no
   result was produced from the wrong target, all target-independent
   diagnostics of run 1 are documented in RUN_NOTES.md. Fixed and
   re-anchored before any certification run.
2. First ANCHOR-BLK formulation was entrywise-false (T L T !=
   blockdiag(tau,tau); the complement block carries S2): caught by the
   instrument's own startup output at run2 HALT (exit 3, no data
   produced), corrected to the two-anchor chain (L already IS
   blockdiag(Lq, S2LqS2); Sbig gives the target).
3. tau-7 "natural" slice choice was singular (precondition failure logged
   at run3); re-run with the registered qrcp rule — the only change is the
   slice that had no result anyway.
4. Two missing imports (`csv`, `json`) crashed runs 2/3 after their MAIN
   sections completed identically; no mathematical output affected (run 4
   complete, run 5 = freeze of record with factor CSVs).
5. `rf_replay.py` initially asserted tau7 factors are an EXACT root (wrong
   semantic — they are certified-approximate decimals by frozen
   discipline); corrected to the box semantics (residual < 1e-3 AND
   containment). The retraction is inline in the script comments.
6. Pre-statement section 2.3 said P-pos "+1e-6" — superseded by R3 (1e-8)
   BEFORE compute; the pre-statement is byte-frozen so R3 governs. Also the
   C-tau "64 equations" arithmetic slip — corrected by R6 to 48 BEFORE any
   ladder compute. Both were pre-compute corrections (rule 16 complaint:
   no threshold/domain moved after any result was seen).

## Cost (in-process accounting, rule 17e; run5)

time.process_time() totals: instrument run5 = 415.5 s CPU single-thread
(2x192x192 exact inverses ~61 s of it; 50,000 TRF nfev ~330 s; margins
~15 s). Replay = 42.8 s CPU. NO ps-derived numbers used anywhere. nice -n
10 throughout; Main confirmed no slot serialization needed (no breach).

## Rule-7 scope sentence

Covered: the fixed conjugated Route-F tensor (slices
blockdiag(tau_p, tau_p), proven CP-rank-equivalent to (L_1, L_i, L_j) by
two machine-verified diagonal identities) at rank 13, from the two frozen
seeds via the declared 5-round TRF polish, through the declared QRCP square
slice (192 eqs / 55 frozen / full exact transversality), the fixed 15-rung
Krawczyk ladder, and the four declared controls (tau-7 certify, tau-6
reject-with-certified-exclusion, P-pos certify, P-wrong reject) — exact
fmpq arithmetic on every load-bearing number, outward-arb re-verification
of every passing comparison, in-process CPU accounting only. NOT covered,
hence not excluded by anything here: any root of the rank-13 CP system off
the frozen slice or outside the certified boxes (the 55 frozen coordinates
were NOT searched); any rank-13 witness anywhere; any rank >= 14 statement
about the triple; the full 8-slice T_O window 18 <= R_R(T_O) <= 25
(untouched, nothing here bears on it); border rank; complex decompositions;
other octonion triples; the absence claim does NOT extend beyond the two
declared seeds and their certified neighborhoods.

## Named next action

The declared instruments on the declared seeds are exhausted: no
re-slicing, no seed additions, no box widening is permitted or useful at
this residual scale. The deciding steps that remain inside the frozen
plan-with-a-new-domain discipline are (owner decision required, rule 16
re-scoping):
  N1 (downward, decisive if it works): exact 192x192 certified Krawczyk at
  a candidate whose EXACT residual is at arithmetic zero — requires a
  dramatically better polish (e.g. exact RREF-based CP completion from a
  partial witness, or a structured 4+9-block ansatz exploiting the
  blockdiag(tau, tau) shape with cross-block zeros); the window analysis in
  RUN_NOTES.md gives the exact inequality any such candidate must satisfy
  (||Y g0||_max < rho and 132 |factors| rho^2 < rho - ||Y g0||_max).
  N2 (upward, decisive): exact infeasibility of the rank-13 system
  (Gröbner/elimination over the 247-variable x - 55 = 192 polynomial
  system) — the Route-E budget-fail scale problem, un-budgeted here.
  N3 (model side): attempt a rank-13 witness for the FULL tau ⊥ tau
  (different tensor, 64 equations... owner-decide; the current campaign
  did not touch it).
Escalation: none needed — no published number is touched and no claim
about R_R(T_O) is made.
