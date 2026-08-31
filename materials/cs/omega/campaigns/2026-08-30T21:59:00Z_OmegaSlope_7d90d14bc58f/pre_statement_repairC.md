# PRE-STATEMENT — gate C stage 3 REPAIR + decision run (agent OmegaGateCFinish)

Written 2026-08-30 BEFORE the first repair computation of this agent.
Base: continues frozen campaign cs/omega/campaigns/2026-08-30T21:59:00Z_OmegaSlope_7d90d14bc58f
(OmegaSlope, FAILURE TO CERTIFY — V1 failed on the named single defect).
The repair is executed in-place on the working copies in cs/omega/src/ and
re-frozen after; the frozen 21:59:00Z copies stay untouched as the defect
record. This file is committed to git and the campaign directory BEFORE the
first computation that depends on the repair.

## 1. The pre-registered repair (exact, fixed now — no alternatives explored)

The level-2 replay in src/gate_c_slope_pass.py never builds
`mat_size_contribution` for PartLv2 shapes. Derived FIRST-HAND from
vxxz24_float.PartLv2.evaluate_post (lines 576-650, file sha
31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890):

For every PartLv2 part t, in Slope form:
- shape_type "022": num_block_contribution = (0,0,0);
  mat_size_contribution = (0, 0, inner) with
  inner = entropy_vec([s0, s0, 1-2*s0]) + 2*ln(q)*(1-2*s0), s0 = split_0
  (POINT param, zero gradient — split_0 is NOT a dist coordinate).
- "112": num_block_contribution = (ln 2, ln 2, entropy_vec([s0,s0,1-2s0]));
  mat_size_contribution = ((1-2s0)ln q, 2*s0*ln q, (1-2s0)ln q).
- "013"/"031": nbc = (0,0,0); msc = (0, 0, ln 2 + ln q).
- "004": nbc = (0,0,0); msc = (0,0,0).
Then BOTH contributions undergo Rot3(., rotate_num) and BOTH are
multiplied by part_frac (a genuine Slope carrying dist gradients).
Known deviation from OmegaSlope's replay found while deriving: the
predecessor's "022" block put the entropy inner term into
num_block_contribution slot 2 (nbc = (0, 0, inner)) — the primary source
clearly puts inner into MAT_SIZE_CONTRIBUTION slot 2 and keeps
num_block_contribution = (0,0,0). The repair fixes BOTH defects
(missing msc, wrong nbc022) in one cutover; V1 below adjudicates.

Also per manifest: the ms loop must skip parts whose msc is absent OR
whose num_block_contribution is [None]*3 (PartZero zeros are explicit
(0,0,0) Slope constants — harmless either way; the loop must add
PartLv2 and PartZero contributions over the right slots, both of which
ARE Slopes after the repair).

Organizational rule (primary-source discipline, cs/README.md rule 1):
this block is derived from vxxz24_float.py first-hand today, not from the
note paraphrase; the parity check in §3 is the control.

## 2. Validations (fixed now)

- V1 (ABSOLUTE gate, pre-registered): the point (r=0) re-aggregation must
  produce Omega_raw within 5e-13 of the frozen RAW endpoint
  2.3715538358350803 (frozen rung-2 campaign; the certified endpoint is
  raw + 1.9380875689560958e-11). If V1 fails: STOP, report FAILURE TO
  CERTIFY with the residual; NO LP run, no partial numbers.
- V2 (consistency, not gating the LP): point-pass coordinate gradient
  bundle midpoint vs coord_grad_full.npy — REPORT ONLY (the reference
  gradient was derived from the WRONG value definition, so agreement is
  not expected; the correct gradient scale is O(2.5), the corrupt one
  O(100s)).
- V3 (control for the repair, rule 14): float parity of the repaired
  replay vs vxxz24_float.PartLv2.evaluate_post: for EVERY level-2 part,
  evaluate the MATLAB-transcribed float evaluate_post at p* and the
  repaired Slope replay's VALUE parts at r=0; require
  max |delta| <= 5e-13 per component over all 3 slots of both
  contributions (spread over all 1377 level-2 parts... the actual part
  count comes from the tree). This is a control the previous run FAILED
  silently: it must be able to DETECT the 022 defect (verify: with the
  022 defect left in — as a diagnostic toggle ONLY, not a run mode — the
  parity residual exceeds the tolerance).
- V4 box zero-crossing guard (IC.ZC_FLAG) clean on the box pass.
- M_low report before/after: 0.14985 (defective) target ~2.09425.
- V5 point defects budget check: reported |Omega_raw + defects/M_low -
  certified frozen endpoint| informational only (defects are the FROZEN
  c + ceq + Schoenage quantities by the frozen protocol; this campaign
  does not re-derive them and no certified endpoint is claimed here —
  the certified endpoint remains the frozen 2.3715538358544617).

## 3. Post-repair decision stages (run only if V1 passes)

Exactly as pre-registered in pre_statement_gateC2.md §2c/§4 (domain D, ANY-y
LP instrument, HiGHS u-form handling preserved, adversarial duals
(LP_duals), R2slack, 12.8x cap statement), with amendments fixed now:
- The float reference gradients V2/V3 of OmegaSlope CANNOT serve as
  controls after the repair (they encoded the wrong value definition).
  The repaired gradient control is the SUCCESS criterion of the LP stage
  itself: the certified slope bundle over the box must have
  ratio (enclosure width over D) / (1.5816497000997742e-07 signal)
  FINITE and REPORTED — that ratio is the pre-registered adjudicator.
- Decision logic unchanged: PASS-a iff min over D of the certified
  enclosure >= Omega0; PASS-b iff an exact delta in D with certified
  Omega < Omega0 is exhibited; otherwise FAIL with the width/signal
  ratio and the exact step where the width blows up. No domain move,
  no re-centering, no budget extension (rule 16).
- FROZEN RAW check: the certified endpoint at the BOX CENTRE from the
  frozen rung-2 re-evaluation is 2.3715538358544612; the certified
  endpoint FROZEN for V1 anchoring is 2.3715538358544617 (+5e-16 display
  granularity) and the RAW component frozen is 2.3715538358350803. The
  5e-13 tolerance is against the RAW value.

## 4. Falsifiable outcomes (fixed now)

- V1 pass + LP decision run -> the pre-registered three-way outcome
  (PASS-a certified local optimality over D at radius 1e-7 NAMED;
  PASS-b exact improving direction with certified decrease; FAIL with
  ratio + step).
- V1 fail -> FAILURE TO CERTIFY, residual stated, no LP, no numbers
  propagated to any index.
- Scope honesty (binding, from pre_statement_gateC2.md §5): even total
  success in D caps 12.8x short of the published 2.37155181; this route
  can NEVER produce a record; no record claim may appear anywhere; a
  certified enclosure is not a new bound; nothing here bears on the
  validity of the VXXZ24 paper.

## 5. Rule-7 sentence (the scope of any certified sentence this campaign)

Identical to pre_statement_gateC2.md §7: the swept set is
D = {delta in R^45 : A*delta = 0 exactly, |delta|_inf <= 1e-7} around the
region-0 glob dist block of the VXXZ24 K100_2.37155181 released vector,
A = the frozen 27x45 margin matrix (rank 24 exact over Q; kernel basis V
exact rational, NOT a float SVD basis), under the transcribed 3-region
single-p_comp program at max_level 3, q = 5; NOT the record
omega < 2.371177 (parameters unpublished, never searched), NOT any point
outside D, NOT other rungs, other max_level/q, multi-block classes, or
other constructions.
