# Manifest — gate C stage-3 repair campaign (agent OmegaGateCFinish, 2026-08-30)

Campaign: apply the OmegaSlope-named PartLv2 `mat_size_contribution` repair,
pass V1, then run the coded decision stages. Pre-statement committed to git
**b41f384** (2026-08-30T21:42Z) BEFORE the first repair computation:
`cs/omega/campaigns/2026-08-30T21:59:00Z_OmegaSlope_7d90d14bc58f/pre_statement_repairC.md`.

## Outcome: FAILURE TO CERTIFY (V1 FAILED; reported per the pre-registered
## FAIL mode — no LP run, no decision on the 21-dim nullspace question,
## which therefore remains UNCHANGED AND OPEN)

### V1 result (pre-registered: point re-aggregation must hit the frozen RAW
### endpoint 2.3715538358350803 within 5e-13)

    Om_raw (repaired replay, radius 0) = 3.272425321778391
    FROZEN_RAW                         = 2.3715538358350803
    V1 residual                        = +0.9008714859433109   → FAIL

    For reference, the DEFECTIVE predecessor value was 41.9962 (residual
    +39.6247); the repair moved the residual from +39.62 to +0.90 but did
    NOT close it.

### M_low before/after (pre-registered report item)

    M_low before repair (OmegaSlope): 0.14985
    M_low after repair:               2.0942543887102634
    frozen protocol M_low:            2.0942543887101843   (matches, 8e-17)
The named repair (PartLv2 msc + related) WORKED as specified: M_low is
restored exactly, and level-2 aggregation parity vs the frozen protocol
holds to 6.9e-18 over all 1377 PartLv2 parts x 6 slots (harness
parity_lv2_check.py in this directory; defect-detection control inside it
confirms the harness detects the named defect at 4.4e-2 >> 5e-13).

### V1 failure localization (next agent's starting point — evidence below is
### float-level DIAGNOSTIC, not certified)

Per-block R-candidate comparison vs the frozen protocol reproduced
standalone (Om_raw repro 2.3715538358351957, delta vs frozen 1.15e-13):

  - level-2 R branch (nb2): 0.5883309785276518 vs frozen 0.5883309785276386
    → MATCHES (the repaired part is correct).
  - glob num_block values: MATCH (e.g. r=0 nb 0.49586563875 exactly).
  - level-3 num_block lows: MATCH (0.234012078.. both).
  - MISMATCH 1 (dominant): level-3 hash_penalty. Worst per-part:
    slope 1.4728e-01 vs frozen 7.5e-09 for part shape (3,3,2)/(2,3,3).
    Root cause found and FIXED as a first candidate (lam_sum sign: the
    replay had g_val = −lam_sum − 1 + margins, the frozen protocol uses
    +lam_sum − 1 + margins — confirmed by numeric reconstruction:
    frozen diff −4.4891447 − (−4.489145) ≈ −3e-16 vs slope −3.525): after
    the fix the L3 block residual dropped from 2.3 down to 2.05e-6 total.
  - MISMATCH 2 (remaining, unresolved at budget): glob hash_penalty /
    Lemma-1 eps block. Both trees' eps MAX over shapes agrees numerically
    (3.0767e-6 at shape (8,0,0) r=2 etc.), but the frozen tree's
    contribution evaluates to ~1.78e-15-scale while the slope replay keeps
    3.08e-6 → the frozen evaluate-post dot-path computes eps differently
    from the standalone replication (suspect: in the interval-driven tree
    the `hi <= 0` guard and/or pos() of point-interval diffs behaves
    differently than the slope replay's re-derived max-over-shapes; the
    exact divergence site is NOT identified — this is the named next
    blocker). Full-tree parity harness: parity_full_tree.py (worst after
    the lam_sum fix: glob pen r2, 2.05e-6).

### What ran (all inside IC.prec(), MID=300, threads=1, OMP=1)

- Point pass (r=0) + box pass (r=1e-7) both ran end-to-end post-repair
  (~111 s each) with all 1377 PartLv2 parts carrying complete_split,
  num_block_contribution, mat_size_contribution.
- V1 point re-aggregation: FAILED as above. Per pre-statement §2/§4 and
  the parent pre-statement §4: STOP — no LP, no PASS-a/PASS-b, no partial
  numbers promoted. The 21-dim nullspace question is UNCHANGED AND OPEN.

### Repairs applied to src/ (working copies; frozen here as byte-copies)

Derived FIRST-HAND from vxxz24_float.PartLv2.evaluate_post (sha
31657c92a1f841805e7355d5b5390f76e0d7ca2d04b1a3f2ce9fc204d9d5d890,
byte-identical to the frozen campaign copy):
1. PartLv2 replay block rebuilt (previously: no msc at all; 022 inner
   wrongly in nbc slot 2 with coefficient ln q instead of 2 ln q; 022
   complete_split marginals summed into one vector; no part_frac
   multiply): now mirrors the primary exactly (022 nbc=(0,0,0),
   msc=(0,0,ent3+2lnq(1−2s0)); 112 msc = ((1−2s0)lnq, 2s0·lnq, (1−2s0)lnq);
   013/031 msc=(0,0,ln2+lnq); 004 zeros; distinct cs1/cs2/cs3 with the
   primary's EncodeCSD entries; Rot3 on both contributions then frac
   multiply on both; complete_split assigned once AFTER Rot3c).
2. gate_c_endgame.aggregate: level-2 R branch now pins its argument and
   appends the (nb2, zero-rhs) slope pair; endpoint-slope numerator fixed
   from −R'M − T·M' to −R'M − Ω·M' (Leibniz error; its value part was
   −T/M ≈ −3.72 instead of Ω ≈ 2.37).
3. Part-level penalty lam_sum sign fixed (−lam_sum → +lam_sum), matching
   the frozen protocol as reproduced standalone.
4. Dead `if False else` debris in the glob num_block path removed
   (clean cutover; behavior of the surviving branch unchanged).

### Contents
- pre_statement_repairC.md — the pre-registration (git b41f384)
- gate_c_slope_pass.py, gate_c_endgame.py — repaired scripts as run
- gate_c_slope_core.py, vxxz24_float.py, interval_core.py — unchanged
  machinery byte-copies (sha of core/vxxz equal the frozen OmegaSlope's)
- parity_lv2_check.py — V3 harness (per-part parity + defect control): PASS
- parity_full_tree.py — V6 full-tree parity harness: FAIL at 2.05e-6
  (glob pen eps block) — the named next blocker
- checksums.sha256 — sha256 of every file here

### Rule 7 sentence (mandatory)

Nothing in this campaign is certified beyond the scripts' internal
arithmetic: the attempted swept set was D = {delta : A·delta = 0 exactly,
|delta|_inf ≤ 1e-7} around the region-0 glob dist block of the VXXZ24
K100_2.37155181 released vector (A = the frozen 27x45 margin matrix,
rank 24 exact over Q — re-verified over Q this session; kernel basis
exact rational, max|A·V| over Q = 0, dim nullspace = 21, V spans it,
max|v·1| over Q = 0), under the transcribed 3-region single-p_comp
program at max_level 3, q = 5 — and even over D no decision was reached
(V1 failed on a residual aggregation defect: the repaired level-2 block
is parity-exact, but the Lemma-1 eps/penalty aggregation in the slope
replay still disagrees with the frozen protocol at the 2e-6 level). The
statements here do NOT apply to the matrix multiplication exponent record
omega < 2.371177 (parameters unpublished as of 2026-08-30, never
searched), to any point outside D, to any other rung's released vector,
to any other max_level/q regime, to multi-block or cross-block direction
classes, or to any construction differing from the enclosed laser-method
program (including the 2026 asymptotic-rank / centroid-based construction
improvements, outside every statement here). Even total success in D caps
12.8x short of the published 2.37155181, so this route can NEVER produce
a record and no record claim may appear anywhere.

### What was NOT swept (explicit)

- No certified decision over D was produced: neither PASS-a nor PASS-b;
  the LP/decision stages of gate_c_endgame.py were NOT run (V1 gate
  failed).
- Regions 1 and 2 glob dist blocks; multi-block/cross-block classes; the
  omega x single_mat_size Schoenage-line 2-coordinate class; all other
  parameter coordinates — held at p* points.
- The published 2.37155181 and the record 2.371177: not touched, not
  searched; nothing here bears on the validity of the VXXZ24 paper (no
  escalation was triggered — no claim about the paper was written).
- No subdivision of D, no Monte Carlo, no probes outside D, no float SVD
  basis substitution (the exact rational kernel basis was used as-is and
  re-verified over Q).

## AMENDMENT (2026-08-30, post-yield; Main-approved handoff question, top of note)

**MAIN'S ARITHMETIC OBJECTION — the (a)/(b)/(c) question MUST be answered
before any successor writes repair code:**

- The named blocker (glob Lemma-1 eps, frozen ~1.8e-15 vs replay ~3.08e-6,
  a 1.71e9x path discrepancy) has ABSOLUTE size 3.08e-6, but the V1
  residual to close is 0.9008715 — eps alone would need ~2.925e5
  downstream amplification to account for it. Hence:
  (a) eps is genuinely amplified ~3e5 through a division/logarithm
      (possible — but must be DEMONSTRATED, not assumed), OR
  (b) a FOURTH defect, still unidentified, carries most of the 0.9008715,
      OR
  (c) both.
  Fixing eps blindly risks a third FAILURE TO CERTIFY at the same residual.
- Cheap diagnostic prescribed by Main (bisect on STAGES, not code): take
  the frozen tree's own intermediate values and substitute them one stage
  at a time into the replay's aggregation; the first stage where a frozen
  substitution collapses the residual localizes the defect. ~a dozen
  stages, not a search over 1377 parts. Frozen intermediates are exposed
  by the float tree and endgame_partial_checkpoint.json.
- Magnitude context: the residual 0.9009 is ~38% of the raw value 2.3716 —
  a structural missing/double-counted term, NOT a precision-scale effect,
  consistent with the three defects already found (41.9962 → 3.2724 =
  M_low term; 3.2724 → 2.3716 = one more term of similar magnitude, e.g.
  a wrong R-branch/candidate structure, not eps at 1e-6).

**STATUS AS RECORDED BY MAIN:** gate C stays PARTIAL; the 21-dimensional
question stays OPEN and unchanged; the scoped negative over the six
pre-registered probes stands; the index gains no certified number from
this campaign. Three defects fixed and the 43.98x residual reduction are
recorded as genuine movement.

**Handoff order for a successor:** (1) answer (a)/(b)/(c) via stage-by-
stage frozen-value substitution; (2) only then repair; (3) re-run V1.

## AMENDMENT 2 (2026-08-30 post-yield; Main refinement of the handoff question)

**PROMOTED CONTROL (stronger than an evidence-list item):** the frozen
protocol was reproduced standalone in this session at
Om_raw = 2.3715538358351957 vs frozen raw 2.3715538358350803 —
**delta 1.155e-13, INSIDE the 5e-13 V1 threshold** — on the same machine,
same data, same anchor as the replay that fails. The frozen protocol
PASSES the very gate the replay fails: the defect is decisively isolated
to the slope replay path; environment, frozen data, and the anchor are
all ruled out.

**BLOCKER REFINEMENT — penalty vs eps split (Main):** the "Lemma-1
eps/penalty path" has two components and they must not be conflated:
- eps is O(3.08e-6); it would need ~2.9e5 amplification to carry the
  0.9009 residual — implausible without demonstration.
- The lam_sum PENALTY is naturally O(1), and one sign flip was already
  found and fixed in exactly that term. A penalty mis-aggregation of
  0.9009 requires NO amplification. Likely answer to the (a)/(b)/(c)
  question: (c) both, with the PENALTY carrying essentially all of the
  residual and eps a separate, much smaller, genuine bug.

**ONE-LINE TEST (one pass, settles which component to open):** zero the
penalty contribution in the replay and re-run V1.
- Residual collapses 0.9009 -> ~3e-6 ⇒ the penalty is the blocker,
  eps is the remainder.
- Barely moves ⇒ penalty exonerated; run the stage bisection instead.
(Run it as a DIAGNOSTIC TEST, not a fix — it must not enter any certified
path; per Main this is its thirteenth structural suggestion, the first
twelve were wrong, though this one derives from the parity data.)

**Handoff order (supersedes Amendment 1's order):**
1. Run the one-line penalty-zeroing test (diagnostic only).
2. Open the implicated component (penalty aggregation or stage bisection).
3. Re-run V1; only on a pass proceed to the LP/decision stages.

NOTHING ELSE CHANGES: gate C PARTIAL; the 21-dimensional question OPEN
and unchanged; no certified number entered any index; the scoped negative
over the six pre-registered probes stands; no record claim anywhere.
