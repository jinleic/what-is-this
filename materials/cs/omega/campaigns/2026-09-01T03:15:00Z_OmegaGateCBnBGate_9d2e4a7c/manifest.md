# Manifest — gate C stage 4 feasibility gate (BnB) — SPREAD branch:
# branch-and-bound over D is INFEASIBLE AT THIS COST (agent
# OmegaPenaltyProbe, 2026-09-01; campaign 2026-09-01T03:15:00Z_…,
# actual UTC start stamp)

Headline: **STAGE 1 verdict = SPREAD → QUANTIFIED OBSTRUCTION. No set of
≤ 6 coordinates carries ≥ 90% of the L1 width; the measured k* is 36 of
45 coordinates, each of 40 coordinates carries ≥ 1% of the total
individually, and the width is nearly uniform across [0.664, 4.641].
Uniform (or coordinate-restricted) branch-and-bound to force
width/signal < 1 from width/signal 73.86 needs ~2^(6.2·36) ≈ 2^223
subboxes (even the mildest honest half-width reading ~2^112). STAGE 2
does not open. Frozen as a first-class negative. The 21-dimensional
kernel question is UNCHANGED and OPEN; gate C stays PARTIAL.**

## Protocol conformance

- Pre-statement committed `f8b9ca7` BEFORE any compute; GATE-Z tolerance
  amendment committed `532781d` after attempt 1 measured a 1.2e-14
  last-digit replay drift but BEFORE any verdict compute (logged in
  pre-statement §2; the amended tolerance — max|live−frozen| ≤ 1e-12
  per endpoint, L1 within 1e-6 — is what both gates below test).
- Runner: `stage1_decompose.py`; attempts 1–3 under DISTINCT log files
  (freeze-hygiene rule): `stage1_attempt{1,2,3}_std{out,err}.log`.
  Attempts 1–2 aborted on runner defects (bit-equality gate written
  before the amendment; a scoping slip), tracebacks preserved in their
  own stderr files; no verdict was produced by any aborted attempt.
  Attempt 3 completed EXIT=0.

## GATE Z (cross-run control, per amended tolerance)

- live vs frozen L1: 113.66562530936517 vs 113.66562530936552 —
  agreement 3.55e-13 ≤ 1e-6.
- live vs frozen widths: max |Δlo| = 1.20e-14, max |Δhi| = 1.78e-15 —
  both ≤ 1e-12 (float-extraction replay noise across process images,
  as pre-registered; the DECISION quantities are ordering and
  concentration, both far from any threshold boundary).
- VERDICT margin: the 90% mark is crossed between k=35 (cum 89.902%)
  and k=36 (cum 91.384%); both are an order of magnitude beyond the
  ≤ 6 gate, so float last-digit drift cannot flip the branch.
- Conclusion: the frozen parent enclosure reproduces in-process;
  verdict computed from the FRESH live arrays bit-exactly.

## STAGE 1.1 — per-coordinate width distribution (MACHINE-VERIFIED,
## fresh live arrays; the decomposition of width L1
## 113.66562530936517)

- k* at ≥ 90%: **36 of 45** (cumulative shares: k2 8.17%, k4 15.22%,
  k6 21.85%, k10 34.63%, k23 68.41%, k36 91.38%).
- Coordinate widths span [0.6639861348561744, 4.640574461822501] with
  median 2.3324: 42 of 45 coordinates carry width ≥ 1 unit; 40 of 45
  carry ≥ 1% of the total individually; 23 of 45 carry ≥ 2%.
- Top-6 coords (0, 8, 9, 16, 17, 23) carry only 21.85%. There is no
  dominant set; the width is essentially UNIFORM across the 45
  coordinates.

## STAGE 1.2 — per-term / per-block attribution (from the same live
## aggregate; hull membership is the measurable pre-hull proxy)

- Every R branch requires hulling over ALL THREE candidate tt slots
  (n_possible_minimizers = 3 for R_comp[0/1/2,3], R_comp[0,2],
  R_glob[0], R_glob[1]; 2 for R_glob[2]) with zero_possible false —
  the minimum's branch identity is undecided everywhere, so every
  branch contributes hull width from multiple candidates.
- M: possible_minimizers multiple (recorded in checkpoint); the
  quotient-rule hull adds the M-branch spread on top.

## STAGE 1.3 — branch-straddle component (MACHINE-VERIFIED intervals)

- Max candidate-interval span per branch: R_comp[0,3] 3.826e-7,
  R_comp[1,3] 3.362e-7, R_comp[2,3] 3.604e-7, R_comp[0,2] 6.171e-7,
  **R_glob[0] 3.212e-5**, R_glob[1] 7.622e-7, R_glob[2] 2.051e-6.
  R_glob[0]'s span alone (3.2e-5) is ~203x the 1.5816e-7 signal —
  consistent with the frozen CandidateWitness finding (base-point
  R-branch span ~4.2506e-6 across the seven components; Lemma-1
  residual intervals straddle zero at the base point).
- Mechanism: with ALL slots overlapping at radius 1e-7, min_slope must
  hull 3 candidates per branch; subdivision of one coordinate reduces
  its own entropy widths but cannot decide the branch identity while
  the OTHER 44 coordinates keep their full width — precisely the
  uniformity the 1.1 distribution exhibits.

## The quantified obstruction (frozen)

- Width/signal now: 73.85977508485684 (= 1.1682029111240058e-5 /
  1.5816497000997742e-7, both frozen).
- Bisections per binding coordinate to shrink width ~74x (width ∝
  diameter): ~log2(73.86) ≈ 6.21.
- Binding-coordinate count at the 90% threshold: **k* = 36**.
- Implied subbox count: 2^(6.2×36) ≈ **1.55e67** (full reading);
  2^(6.2×18) ≈ **3.94e33** (mild half-width reading). Compare: this is
  beyond any feasible sweep by ~40+ orders of magnitude even against
  the repo's largest censuses.
- Verdict sentence, verbatim from the pre-statement: **branch-and-bound
  over D at this enclosure technology is INFEASIBLE AT THIS COST — a
  first-class frozen negative, not a failure of the loop.**
- Stage 2 never opened; no subbox compute; no calibration; the exact
  rational kernel basis was not needed and was not exercised.

## What this does NOT say (scope honesty)

- The obstruction binds THIS technology (min_slope hulls over an
  un-subdivided box at radius 1e-7, interval AD slopes). It does NOT
  bound: a different certified-gradient technology with width
  concentrated in few coordinates, an analytic argument that pins the
  branch structure globally (e.g. a certified sign of the Lemma-1
  residual displacement), or subdivision COUPLED with a proof that
  bisection of few coordinates crosses a branch-deciding threshold
  (ruled out HERE only empirically at radius 1e-7).
- Any success in D says nothing about the published 2.37155181 or the
  record 2.371177 (12.808x guardrail: gap 2.0258350805768544e-6 vs
  signal 1.5816497000997742e-7); neither is reachable by this route
  and no record claim may appear anywhere. Gate C stays PARTIAL; the
  21-dimensional kernel question remains UNCHANGED and OPEN.

## Rule 7 sentence (mandatory)

Measured set: the existing certified full-box slope bundle on
D = {delta in R^45 : A·delta = 0 exactly, |delta_i|_inf ≤ 1e-7} around
the region-0 glob dist block of the VXXZ24 K100_2.37155181 released
vector (A = frozen 27x45 0/1 margin matrix, rank 24 exact over Q;
interval core v11 ce70f959…, slope pass c6b8cef2…, point+box pass under
the transcribed 3-region single-p_comp program at max_level 3, q = 5,
radius 1e-7, exact-Arb endpoints). This campaign computed a
decomposition of THAT enclosure's L1 width and branch-straddle
structure and froze the resulting infeasibility arithmetic; it did NOT
run branch-and-bound, subdivision, any new enclosure, any other
parameter block, region, q/max_level, radius, construction, or any
point outside D; the published 2.37155181 and record 2.371177 were not
touched and are unreachable via this route (12.808x guardrail).

## Contents

- pre_statement.md — pre-registration (commit f8b9ca7) + GATE-Z amend
  (532781d, pre-verdict)
- stage1_decompose.py — the runner (as run in attempt 3)
- stage1_decomposition.json — the machine checkpoint (all tables +
  verdict)
- stage1_attempt{1,2,3}_std{out,err}.log — per-attempt logs (distinct
  filenames per the 2026-09-01 freeze-hygiene rule; attempts 1–2
  aborted with tracebacks preserved, no verdict produced)
- checksums.sha256 — written at freeze
