# PRE-STATEMENT — gate C stage 3 DIAGNOSTIC CAMPAIGN (agent OmegaGateCDiag)

Written 2026-08-31T07:57Z, BEFORE the first diagnostic computation of this
agent. Base: continues frozen campaigns
`2026-08-30T21:59:00Z_OmegaSlope_7d90d14bc58f` and
`2026-08-30T22:58:26Z_OmegaGateCFinish_f5269412f0b1` (both FAILURE TO
CERTIFY — V1 residuals +39.62 and +0.9009 respectively). Committed to git
BEFORE any run that depends on it.

## 0. Status of this campaign: DIAGNOSTIC ONLY

Per the Main amendment (manifest.md Amendment 2, one-line test) and the
assigned steering, the FIRST compute of this campaign is a single
diagnostic pass, marked DIAGNOSTIC ONLY in the manifest:

  Zero the lam_sum penalty contribution in the slope replay and re-run V1
  against frozen RAW 2.3715538358350803.

Interpretation, fixed NOW:
- Residual collapses toward ~3e-6 ⇒ the penalty aggregation is the
  blocker; open that component, identify the defect, fix, re-run V1.
- Residual barely moves ⇒ the penalty is exonerated; run stage-by-stage
  frozen-value substitution bisection instead.

## 1. Anchors (fixed before any run; asserted at startup per repo rule #3)

All from the frozen authoritative records:
- FROZEN_RAW = 2.3715538358350803 (stage_b_boxes.json B0_anchor
  omega_cert_upper_raw; V1 pre-registered anchor, tolerance 5e-13).
- FROZEN R_sum_low = 2.817003567460976; FROZEN M_low =
  2.0942543887102634 (B0_anchor record).
- Standalone frozen-protocol control PASSES the same gate:
  2.3715538358351957 (delta 1.155e-13; Amendment 2). The defect is
  isolated to the slope replay path.
- Defective replay state: Om_raw = 3.272425321778391, residual
  +0.9008714859433109 (OmegaGateCFinish).
- Level-2 R branch is parity-exact (0.5883309785276518 vs frozen
  0.5883309785276386, diff 1.32e-14); glob/level-3 num_block lows MATCH.
- Since 0.9008714859433109 = 3.272425321778391 − 2.3715538358350803 and
  R_low(2.3715538358)·M_low(2.0942543887) with the SAME M branch
  (M at tt index 0, the M branch is ~2.0942 vs ~19.90 at branch A1) pins
  the total_R deficit: replay total_R must be frozen_total_R − 0.9009·M
  = 2.8170035675 − 1.8866342744 ≈ 0.93037, i.e. ~1.88662 SHORT against
  the frozen 2.81700. The diagnostic measures whether zeroing the
  lam_sum penalties explains this 1.8866.

## 2. The diagnostic (exact, fixed now)

Single script gate_c_diag_penalty.py:
- Runs the repaired slope replay point pass (radius 0) — code path
  UNCHANGED except a diag-zero toggle on the Part-level
  hash_penalty_term slopes (the only place lam_sum enters the
  aggregation; the glob blocks carry no lam_sum). Toggle read from
  env GATEC_DIAG_ZERO_PEN ∈ {0,1}.
- Recomputes aggregation candidates (level-3 per r, level-2, glob per
  r), reports per-block cand values both with the real penalties and
  with penalties zeroed, and the residual against FROZEN_RAW under
  penalty-zeroing.
- Anchors asserted at startup: FROZEN_RAW, FROZEN_R_sum, FROZEN_M_low.
- Output marked "DIAGNOSTIC ONLY — NOT A CERTIFIED PATH".

## 3. Interpretation table (fixed now)

| residual after zeroing | verdict |
|---|---|
| collapses to O(3e-6) | penalty agg = blocker (answer (b) with the penalty defect) |
| stays O(0.9) | penalty exonerated ⇒ stage bisection next |
| partial movement | (c) both; decompose exact deficit |

## 4. Rule 7 sentence (the scope of any statement in this campaign)

The attempted swept set of THIS campaign remains D = {delta : A·delta = 0
exactly, |delta|_inf ≤ 1e-7} around the region-0 glob dist block of the
VXXZ24 K100_2.37155181 released vector (A = frozen 27x45 margin matrix,
rank 24 exact over Q; kernel basis exact rational) under the transcribed
3-region single-p_comp program at max_level 3, q = 5 — and this campaign
attempts NO certified decision over D; every run here is a flagged
DIAGNOSTIC on the replay's own aggregation arithmetic (point pass,
r = 0), not a certified enclosure. Statements do NOT apply to the
matrix multiplication exponent record omega < 2.371177 (parameters
unpublished, never searched), to any point outside D, to any other
rung's released vector, to any other max_level/q regime, to multi-block
or cross-block direction classes, or to any construction differing from
the enclosed laser-method program. Even total success in D caps 12.8x
short of the published 2.37155181, so this route can NEVER produce a
record and no record claim may appear anywhere.

## 5. What is NOT swept (explicit)

- No certified decision over D: neither PASS-a nor PASS-b; LP/decision
  stages not run (V1 gate pre-registers them behind a V1 pass).
- Regions 1/2 glob dist blocks; multi-block classes; omega x
  single_mat_size class; all other parameter coordinates (held at p*).
- The published 2.37155181 and the record 2.371177: not touched, not
  searched; nothing here bears on the validity of the VXXZ24 paper.
- No float SVD basis substitution; no domain move; no subdivision of D;
  no Monte Carlo; no probes outside D.
