# Manifest — OmegaBranchSix (gate C stage 7): certified subdomain pins
# for the SIX remaining R branches on the kernel-ray family
# — **FROZEN-CERTIFIED** (agent OmegaPcPen, 2026-09-01; campaign
# 20260901T133130Z_f4779759_b81dc930e87c; Main steering "one more attempt,
# then omega rests")

Headline: **ONE of the six remaining branches is CERTIFIED-EXCLUDED on the
exact-kernel ray family {p* + rho·s : rho in [0, 1e-7]}: R_glob[2]'s
candidate 2 is impossible as the R_glob[2] minimizer at BOTH ray
endpoints — the exclusion gap cand2.lo - argmin{cand0,cand1}.hi is
[+4.5483061219230014e-12, +2.0511278540081266e-06] outward at e0 AND e1
(positive low, certified). The other five branches
(R_comp[0,3]/[1,3]/[2,3]/[0,2], R_glob[1]) are Tier-3 (no pin): their
closest miss is R_glob[1]'s CandW-vs-argmin gap straddling at
[-2.539357667283998e-12, +7.622141438568529e-07] — a negative lobe only
16.1x the signal (2.539e-12 vs 1.5816e-7), two orders below the parent
campaign's R_glob[0] failures. The 21-dimensional kernel question
remains UNCHANGED and OPEN; gate C stays PARTIAL; the 12.808x guardrail
is restated — no record claim anywhere.**

## Verdict numbers (all outward-rounded v11 interval endpoints)

| quantity | value (outward) |
|---|---|
| R_glob[2] Tier-2 exclusion gap @ e0 | [+4.5483061219230014e-12, +2.0511278540081266e-06] — POSITIVE LOW, CERTIFIED |
| R_glob[2] Tier-2 exclusion gap @ e1 | [+4.5483061219230014e-12, +2.0511278540081266e-06] — POSITIVE LOW, CERTIFIED (same render; exact intervals agree to render precision) |
| implication | cand_2 (the pen candidate with the frozen 0.496850979549942 low, 2.05e-6 ABOVE the cand0/cand1 pair) cannot be the R_glob[2] minimizer anywhere on the ray family; the possible-minimizer set shrinks {0,1,2} → {0,1} |
| R_glob[1] Tier-2 e1 gap | [-2.539357667283998e-12, +7.622141438568529e-07] — straddles; closest miss, negative lobe 16.1x signal |
| R_comp[0,3] e1 gap lows | [-2.089e-07, +2.104e-07] — Tier-3 |
| R_comp[1,3] e1 gap lows | [-2.962e-07, +2.962e-07] — Tier-3 |
| R_comp[2,3] e1 gap lows | [-2.769e-07, +2.759e-07] — Tier-3 |
| R_comp[0,2] e1 gap lows | [-4.136e-07, +4.136e-07] — Tier-3 |

## Protocol conformance

- Pre-statement committed `24b22b7` BEFORE campaign init; run-dir
  byte-identical copy + PREREG_PROVENANCE.json carrying BOTH init and
  final hashes (per the successor rule provenance_correction_pcpen.md);
  zero amendments were needed — the runner matched the prereg on the
  first full run (after two runner-only crashes, logged per-attempt).
- E0 hash-gate: v11 core ce70f959…, vxxz24_float 31657c92…, kernel basis
  JSON 7d69c1dd…, margin matrix 423dbafd…, parent runner f07b70ae…
  (f07b70ae1127a986ab2a169819448df2c2d0f2e1b62ac41455dcff62fc5c8d43 =
  the parent pcpen.py as frozen) — all asserted in-run before any
  import; params prefix df75ae3a verified.
- No B&B, no LP, no slope pass, no subdivision, no Monte Carlo; pure v11
  interval passes: ONE radius-0 pass (GATE-Z/C1/C2/P1) + ONE drift pass
  (C3) + ONE ray-endpoint pass (P2) + exact-Fraction PV.

## Controls (all PASS, both directions)

- GATE-Z: both frozen records reproduced (7 R-branches within 3.3e-16;
  Om_raw 8.9e-16 <= 1e-6) — same frozen anchors as the parent.
- C1 (known-true accept plant): identity gap 5.79e-13 <= 1e-12; live
  R_glob[0] candidates match the frozen CandW triple within 9e-13.
- C2 (wrong-index reject, region 2 — fresh index): eps_2 live
  3.076684958606817e-06, bit-distinct from eps_0 (2.1502e-6) and eps_1
  (1.1433e-6); matches the in-run glob_r2 Lemma residual.
- C3 (dual-drift reject, region 1 — fresh plant): +1e-3 broadcast on
  glob-1 lam entries moves eps_1 by exactly 4*drift (coupling gap
  4.41e-16) and R_glob[1] by -2*prop*4*drift (predicted -2.667e-3,
  observed -0.002666666666666373, gap 2.9e-16).
- C1a: 0 nonzero-width leaves among all non-dist0 blocks at the ray
  endpoint; eps_0 bit-constant there (fourth confirmation).
- PV: A·delta = 0 (27 rows, exact), sum(delta) = 0, |delta_i| <= 1e-7,
  positivity — the swept endpoint is a certified kernel member.

## P (per-branch tiers, per pre_statement §2)

Tier arithmetic at BOTH endpoints; final tier = min(e0, e1):

| branch | t* | tier e0 | tier e1 | final | evidence (e1) |
|---|---|---|---|---|---|
| R_comp[0,3] | 2 | 3 | 3 | 3 | gaps ±2.09e-7 straddle |
| R_comp[1,3] | 2 | 3 | 3 | 3 | ±2.96e-7 |
| R_comp[2,3] | 1 | 3 | 3 | 3 | ±2.77e-7 |
| R_comp[0,2] | 1 | 3 | 3 | 3 | ±4.14e-7 |
| R_glob[1] | 2 | 3 | 3 | 3 | exclusion straddles at -2.54e-12 lo |
| R_glob[2] | 1 | 2 | 2 | **2** | **exclusion certified: +4.55e-12 lo** |

R_glob[2] detail: at both endpoints the three candidate intervals are
(0.49684892842286427..0.49685097954616997), (0.49684892842208755
..0.49685097954539326), and the singleton cand_2 = 0.49685097954994156 —
cand_2's LOW sits 4.55e-12 above the hull-top of {cand_0, cand_1} with
the opposite margin 2.05e-6, so the exclusion inequality
cand2.lo - max(cand0.hi-argmin).hi > 0 holds OUTWARD on the whole ray
family. This is a genuine, certified, subdomain-restricted structural
statement — the second branch ever pinned on this target, and the FIRST
Tier-2 exclusion.

## What this freezes

1. **Certified subdomain statement:** on the exact-kernel ray family
   (rho in [0, 1e-7]), R_glob[2]'s possible-minimizer set is {0, 1}
   (cand_2 excluded, certified width margin 4.55e-12..2.05e-6). Any
   successor certified-gradient/LP instrument over the kernel slice can
   drop cand_2 from the R_glob[2] hull — one candidate fewer in the
   hull-cost budget, on the branch family whose whole-box spans are
   already 13x-signal class.
2. **Closest-miss register:** R_glob[1]'s exclusion gap misses by a
   negative lobe of 2.54e-12 = 16.1x signal — quantified; a successor
   with a tighter nb/R_glob[1] hull could plausibly certify it.
3. **Tier-3 negatives:** the four comp branches' argmin identities
   straddle at ~2-4e-7 — recorded with intervals; the pc/pen prefix
   race carries them (eps exonerated per the parent campaign).
4. **Instrument:** the six-branch candidate-triple aggregation extends
   the parent's single-branch recorder; controls all re-passed with the
   plants re-aimed (region-2 lookup, region-1 drift) — reusable.
5. The 21-dimensional kernel question: **UNCHANGED and OPEN** — this
   campaign trimmed ONE candidate from ONE branch's hull on ONE ray
   family; it decides nothing about improving directions.

## Rule 7 sentence (mandatory)

The measured set of THIS campaign: the candidate-order inequalities of
the six branches R_comp[0,3], R_comp[1,3], R_comp[2,3], R_comp[0,2],
R_glob[1], R_glob[2] at two points (rho = 0 and the PV-certified exact
dyadic ray endpoint rho = 1e-7) of the certified kernel ray family
{p* + rho·s : rho in [0, 1e-7]} on the region-0 glob dist[0] block of
the VXXZ24 K100_2.37155181 released vector (params sha df75ae3a…) under
the transcribed 3-region single-p_comp program at max_level 3, q = 5,
interval core v11 ce70f959…, vxxz24_float 31657c92…, kernel basis
7d69c1dd…, margin matrix 423dbafd…, every other parameter block at p*
point values (width-0 in-run), radius exactly binary64 1e-7. No other
branch (R_glob[0] frozen falsified), no other radius, subdomain,
region, parameter block, subdivision, Monte Carlo, LP, or gradient
machinery. The published 2.37155181 and the record 2.371177 were not
touched and are unreachable via any outcome of this campaign — the
12.808x guardrail (certified gap 2.0258350805768544e-6 = 12.808x signal
1.5816497000997742e-7) survives every branch's verdict; the certified
R_glob[2] exclusion lives at the 2.05e-6 hull scale and produces NO
exponent improvement. Gate C stays PARTIAL; the 21-dimensional kernel
question remains UNCHANGED and OPEN; no record claim anywhere.

## Attempt history

- Two runner-only crashes inside attempt 1's build (NameError `cand`
  after the comp-triple rewrite clobbered the frozen min subtree;
  Ival-from-float in tier_calc) — fixed, no verdicts, engine untouched.
- attempt 1 (attempt1_stdout/stderr.log): EXIT=0, VERDICT
  FROZEN-CERTIFIED, pinned R_glob[2].
- authoritative in-dir rerun (authoritative_stdout/stderr.log): EXIT=0,
  decision fields byte-identical (branchsix_checkpoint.json vs
  checkpoint.json compared key-by-key sans timestamps — identical).

## Contents

- pre_statement.md (byte-identical prereg, both hashes ea864a71…),
  PREREG_PROVENANCE.json (dual-hash rule applied)
- branchsix.py (runner as run), branchsix_checkpoint.json (authoritative
  in-dir checkpoint), checkpoint.json (authoritative rerun's copy)
- attempt1 + authoritative stdout/stderr logs, attempt_log.md
- manifest.md / manifest.json / sha256s.txt / status.json
