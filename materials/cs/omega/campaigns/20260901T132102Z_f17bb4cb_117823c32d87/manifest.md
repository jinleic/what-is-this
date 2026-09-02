# Manifest — OmegaPcVersusPen (gate C stage 6): certified subdomain pinning
# of the p_comp-vs-penalty ordering on the exact-kernel ray family
# — **BOX-CONTRADICTION** (agent OmegaPcPen, 2026-09-01; campaign
# 20260901T132102Z_f17bb4cb_117823c32d87)

Headline: **The registered successor pin — the p_comp-vs-penalty ordering
G1/G2 restricted to the exact 21-dimensional kernel slice of the dist[0]
box (a single certified ray family p* + rho·s, s = the frozen
CandidateWitness direction, sum(s)=0 and A·s=0 over Z, radius exactly
1e-7 = the box-wall-bound feasible radius) — FAILS with a quantified
margin: at the ray endpoint rho = 1e-7 both order inequalities take
outward intervals G1 in [-1.8897456060889468e-05, +1.6396874550250076e-05]
(width 3.5294330611139544e-05 = 223.15x the 1.5816497000997742e-07 signal)
and G2 in [-1.9088262419827114e-05, +1.6587298910286668e-05] (width
3.5675561330113786e-05 = 225.56x signal). Both straddle zero. The
eps_0 exonorated by the parent campaign is RE-CONFIRMED live
(bit-constant 2.1502086942121845e-06 at the ray endpoint too); the
straddle carrier is the pc_0-vs-pen_0 box motion itself. The 21-dim
kernel question remains UNCHANGED and OPEN; gate C stays PARTIAL.**

## Verdict numbers (all outward-rounded v11 interval endpoints)

| quantity | value (outward) |
|---|---|
| G1 @ e0 (rho=0, the p* point) | [-1.4334718839595232e-06, +5.788488742896399e-13] — the frozen 2*prop_0*eps_0 charge floor, negative |
| G1 @ e1 (rho=1e-7 on the ray) | [-1.8897456060889468e-05, +1.6396874550250076e-05], width 3.5294330611139544e-05 = 223.14884647917356 x signal |
| G2 @ e1 | [-1.9088262419827114e-05, +1.6587298910286668e-05], width 3.5675561330113786e-05 = 225.5591824653922 x signal |
| incumbent whole-box widths (falsified parent) | G1/G2 width 5.63e-05 = 355.9x signal |
| kernel-slice / whole-box width ratio | 3.5294e-05 / 5.6292e-05 = 0.6270 (62.7% — collapsing to the kernel ray shrinks the straddle only 1.59x) |
| R_glob[0] candidate hull span at e1 | 2.1108752192733893e-05 (parent whole-box: 3.2123059392596964e-05; 1.521x reduction, still 133.46x signal = 65.7% of the parent span) |
| eps_0 at the ray endpoint | 2.1502086942121845e-06 — BIT-IDENTICAL to the radius-0 value and to the parent's boxed value |
| Om_raw at the ray endpoint (uncertified link, reference only) | 2.3715605661603827 |
| M_low at the ray endpoint | 2.094252326247832 |

## PFIRST (the certified feasible ray radius) — exact-Fraction walls

- max |s_i| over the 45 ray coordinates = 1 exactly (s in {-1,0,+1}^45;
  the frozen CandidateWitness direction), so the BOX wall binds:
  R_box = 1e-7 / 1 = Fraction(944473296573929, 9444732965739290427392)
  EXACTLY (binary64 1e-7 = 0x1.ad7f29abcaf48p-24 is a dyadic).
- positivity wall (from the 45 exact dyadic dist0 coordinates):
  4771454744492531 / 295147905179352825856 = 1.6166317499672092e-05 —
  EXACTLY the frozen CandidateWitness positivity_ceiling (independent
  re-derivation, not copied), and NOT binding.
- r0 = 1e-7 exactly; the prereg's expected 1.0607145909782878e-06 was a
  transcription error (that number is R/16 of the CandidateWitness LP-ray
  steps, a different quantity) — corrected in Amendment 1 (commit 2bc7628)
  BEFORE the verdict compute, with the runner OPEN-at-E-FIRST gate
  enforcing the correction mechanically.

## PV (certified kernel membership of the swept endpoint)

delta = r0·s satisfies, by EXACT integer/rational arithmetic in-run:
A·delta = 0 on all 27 rows; sum(delta) = 0; every |delta_i| <= 1e-7
(= the certified box); every center_i + delta_i > 0 (positivity);
max |delta_i| = 1e-7 = the radius. The swept set is inside the exact
kernel slice D of pre_statement §1. A·s = 0 column-wise on the sha-pinned
exact rational kernel basis V (21 columns). A's rank re-proved exactly 24
over Q by Fraction elimination in-run.

## Protocol conformance

- Pre-statement committed a7c80de BEFORE any campaign existed (prereg sha
  443fc302…, byte-identical copy in this dir, PREREG_PROVENANCE.json
  binding source commit a7c80de + sha256).
- Amendment 1 committed 2bc7628 BEFORE the first verdict-bearing compute
  (attempt 4). It corrected ONLY the registered E-FIRST expectation
  (r0 = 1e-7 box-wall-bound; the prior number was a transcription error)
  and documented the e1-pass instrument detail; no inequality, control,
  or adjudication semantics changed.
- Amendment 2 committed 7305623 (display-doc: hash-text correction in this
  copy; the E0 machine gate always asserted the registered 31657c92 hash).
- No B&B run, no LP solver call, no slope pass, no subdivision, no Monte
  Carlo, no float-precision shortcuts — a pure interval campaign (v11
  core, hash-verified at startup, plus vxxz24_float, kernel-basis JSON,
  margin-matrix npy, and the frozen branchpin.py dependency-of-record,
  all hash-asserted in-run per E0).

## Attempt history (per-attempt logs under DISTINCT filenames — the
## 2026-09-01 freeze-hygiene rule; every attempt's stderr preserved)

- attempt 1 (`attempt1_stdout.log`, `attempt1_stderr.log`): E0 PASS;
  crashed with NameError T0 (runner defect, engine untouched). No verdict.
- attempt 2 (`attempt2_stdout.log`, `attempt2_stderr.log`): crashed
  NameError SCRATCH (same runner-defect class, fixed). No verdict.
- attempt 3 (`attempt3_stdout.log`, `attempt3_stderr.log`): OPEN at
  E-FIRST — exposed the prereg's wrong r0 expectation AND two runner
  comparison/indexing bugs in the fresh exact-wall block (inv_smax
  min/max logic, positivity index shadow). Exact recomputation produced
  R_exact = the frozen positivity ceiling to the bit. No verdict.
  (Runner repairs + E-FIRST correction committed as Amendment 1,
  2bc7628, BEFORE any verdict compute.)
- attempt 4 (`attempt4_stdout.log`, `attempt4_stderr.log`): crashed at
  frac_dyadic's q>0 assert (7 of the 45 ray offsets are exactly 0 — a
  zero-width dyadic radius is the correct object; zero case added). No
  verdict (crash predates the P2 aggregation read). Controls GATE-Z,
  C1, C2, C3, P1, PV all PASSED in this attempt already.
- attempt 5 (display-key fix run: `report["P"]` → `report["P2"]`
  post-aggregation typo, runner display only): EXIT=0, VERDICT
  BOX-CONTRADICTION. Logs for this physical run are
  `authoritative_stdout.log`/`authoritative_stderr.log` (the second
  in-dir run, post-Amendment-2) — see the authoritative entry below;
  the two runs' decision fields are identical (deterministic).
- authoritative re-run in THIS directory (`authoritative_stdout.log`,
  empty `authoritative_stderr.log`): EXIT=0, decision fields identical
  to attempt 5's (deterministic; `pcpen_checkpoint.json`).

## GATE-Z (engine equivalence, both frozen records) — PASS

- with-lemma1: all 7 R branches within 3.33e-16 of the frozen
  stage_b_rung2_with_lemma1.json values (gate 1e-12); Om_raw
  2.3715538358350816 vs frozen 2.3715538358350807 (8.9e-16 <= 1e-6);
  eps0_live bit-equals the frozen 2.1502086942121845e-06;
  M_low 2.094254388710263.
- without-lemma1: all 7 branches within 3.3e-16; Om_raw
  2.3715518061863823 vs frozen ...814 (8.9e-16).
- params sha256 df75ae3a…da93d verified in-run (E0 + PFIRST-params).

## C1 (identity-accept, radius 0) — PASS

- lhs = cand1.lo - cand0.lo (frozen CandW base record) = -1.4334718842490268e-06
  vs registered -2*prop_0*eps_0 = -1.4334724628081230e-06: identity gap
  5.785590961461637e-13 <= 1e-12.
- live candidates reproduce the frozen triple (all six endpoints within
  9e-13 of the frozen CandW values).

## C2 (wrong-index reject) — PASS

- Region-1 dual lookup returns eps_1 = 1.143325024764798e-06, bit-for-bit
  the frozen with-lemma1 glob_r1 record, distinct from eps_0 (1.0e-6
  apart). The machinery reads the correct lam blocks.

## C3 (dual-drift reject, +1e-3 broadcast on every glob-0 lam entry) — PASS

- Drifted eps boundary 0.004002150208693772 = eps_0 + 4*drift
  (coupling gap 4.41e-16 <= 1e-12): the eps term is FULLY dual-coupled
  (re-confirms the parent's C2 exactly).
- Branch response: R_glob[0] moves -2.666666666666373e-03 vs predicted
  -2*prop*4*drift = -2.666666666666667e-03 (gap 2.9e-16): the pc-side
  branch identity responds linearly to the dual drift.

## C1a (leaf-width-zero invariance at the ray endpoint) — PASS

- 0 nonzero-width leaves among ALL non-dist0 blocks (region_prop, dist_1,
  dist_2, dist_max_0..2, all lam_margin/lam_sum) in the ray pass:
  the boxed quantity is exactly the 45 dist0 coordinates.
- eps_0 at the ray endpoint BIT-EQUALS the radius-0 value — the parent
  campaign's exonoration REPRODUCES on the kernel slice itself.

## P (the ordering pin over the kernel ray family) — **BOX-CONTRADICTION**

Adjudication per pre_statement §8 as amended: the PASS branch required
BOTH G1@e1.lo > 0 AND G2@e1.lo > 0 outward. Both lows are NEGATIVE
(-1.890e-05, -1.909e-05) and both intervals straddle zero with widths
~3.53e-05 / 3.57e-05. The failing terms (registered FAIL deliverable):
- G1: interval [-1.8897456060889468e-05, +1.6396874550250076e-05],
  width 3.5294330611139544e-05, width/signal 223.14884647917356.
- G2: interval [-1.9088262419827114e-05, +1.6587298910286668e-05],
  width 3.5675561330113786e-05, width/signal 225.5591824653922.
- interpretation: the pc-vs-pen ORDER is not fixed even on the certified
  1-dimensional kernel subdomain at the registered radius — the negative
  lobe of G1 (|−1.89e-05|) is 119.5x the first-order signal. The order
  cannot be pinned on ANY subdomain that contains this ray family
  (which every 21-dim kernel neighbourhood does, since s is a kernel
  basis combination: the negative-lobe behaviour is in the interior of
  the kernel slice's feasible set, not an artifact of the box superset).

## What this freezes

1. **Certified local statement (fail-fast negative, registered FROZEN-
   NEGATIVE branch):** on the exact-kernel ray family
   {p* + rho·s : rho in [0, 1e-7]} (sum(s)=0, A·s=0 over Z, positivity
   and box certified), the R_glob[0] p_comp-vs-penalty candidate order is
   NOT pinned: both order inequalities straddle zero at the ray endpoint
   with widths 223.15x / 225.56x the registered first-order signal. The
   named mover (pc_0's box motion vs pen_0's) is REMOVED from the
   "pinnable by domain restriction" list: the parent campaign's whole-box
   straddle was NOT a box artifact — 62.7% of it survives on a single
   certified kernel ray, so no ordering pin via subdomain restriction of
   THIS kind can collapse the R_glob[0] hull. Any successor must either
   (a) certify ordering on a DIFFERENT branch (the other six: six have
   whole-box spans <= 2.06e-6 — two orders smaller than R_glob[0]'), or
   (b) attack the ordering with genuinely different machinery (e.g. a
   certified derivative sign of the difference functionals along s).
2. **Structural re-confirmation:** eps_0 is bit-constant over the kernel
   ray too (third independent confirmation: parent box, this ray) — the
   Lemma-1 residual is FULLY exonerated for every future analysis.
3. **Instrument debt retired:** the exact-Fraction feasible-radius
   instrument reproduces the frozen CandidateWitness positivity ceiling
   to the bit from first principles (independent re-derivation), and the
   PV membership proof (A·delta==0, sum==0, |delta|<=radius, positivity)
   is a reusable certified kernel-membership gate for any future
   subdomain campaign on this target.
4. The 21-dimensional kernel question: **UNCHANGED and OPEN** (this
   campaign measured ONE certified ray family inside the slice; the
   remaining directions are untouched and unmeasured).

## Rule 7 sentence (mandatory)

The measured set of THIS campaign: the R_glob[0] branch candidate-order
inequalities G1/G2 at two points of the certified kernel ray family
{p* + rho·s : rho in [0, 1e-7]} — rho = 0 (the released p*) and
rho = 1e-7 (exact dyadic ray endpoint, PV-certified kernel member) —
on the region-0 glob dist[0] block of the VXXZ24 K100_2.37155181 released
vector (params sha df75ae3a…), under the transcribed 3-region
single-p_comp program at max_level 3, q = 5, interval core v11
ce70f959…, vxxz24_float 31657c92…, kernel basis
7d69c1dd…, margin matrix 423dbafd…, every other parameter block at p*
point values (C1a-verified width-0), radius exactly binary64 1e-7. No
point outside the ray family and its rho=0 endpoint, no other branch of
R, no other region, no other parameter block, no other radius, no
subdivision, no Monte Carlo, no LP, no slope/gradient pass. The
published 2.37155181 and the record 2.371177 were not touched and are
unreachable via any outcome of this campaign (12.808x guardrail:
certified-gap 2.0258350805768544e-6 = 12.808x signal 1.5816497000997742e-7
— every outcome of this campaign lives at the ~1e-5 hull scale, orders
OUTSIDE the guardrail in the wrong direction, so no record claim may
appear anywhere). Gate C stays PARTIAL; the 21-dimensional kernel
question remains UNCHANGED and OPEN.

## Contents

- pre_statement.md — byte-identical copy of the committed prereg
  (commits a7c80de + AMC1 2bc7628 + AMC2 7305623)
- PREREG_PROVENANCE.json — source-commit binding for the prereg copy
- pcpen.py — the runner as run (authoritative in-dir bytes)
- pcpen_checkpoint.json — the authoritative in-dir run's checkpoint
  (decision fields identical to the scratch-equivalent fixed run)
- attempt1..4 stdout/stderr logs + authoritative stdout/stderr logs —
  per-attempt preservation per the freeze-hygiene rule (attempt history
  above; attempt 5's physical logs are the authoritative_* files)
- precompute_note.md — append-only runner-defect and correction ledger
- manifest.md / manifest.json / sha256s.txt / status.json — freeze
