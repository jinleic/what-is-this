# Manifest — OmegaGateCBranchPin (gate C stage 5): analytic branch-pinning
# of the R_glob[0] hull over D — **FALSIFIED** (agent OmegaNext, 2026-09-01;
# campaign 20260901T102048Z_b555441e_d3303018f163)

Headline: **The registered analytic branch-pin of R_glob[0] FAILS with a
quantified margin: both order inequalities G1/G2 straddle zero over the
full 1e-7 dist[0] box, interval [-2.9602e-05, +2.6690e-05] — 355.9x the
1.5816e-07 signal in width, 187.2x in the negative lobe. `eps_0` is
CERTIFIED CONSTANT over D (bit-identical to the frozen point value
2.1502086942121845e-06) and is therefore EXONERATED as the box-straddle
driver; the straddle is carried by the pc_0 (p_comp) hull moving ~1.54e-5
against the pen hull moving ~1.14e-5 over the box. The 21-dimensional
kernel question remains UNCHANGED and OPEN; gate C stays PARTIAL.**

## Protocol conformance

- Pre-statement committed `754b66d` BEFORE any compute (Main strict order);
  **Amendment 1** committed `68a65d9` BEFORE the first runner execution
  (corrects the pin direction: the Lemma-1 charge is what flips the
  minimizer TO the penalty candidate at p*; the registered pin target
  became order inequalities G1/G2 — both gaps strictly positive, outward
  — for the p_comp candidate to be the unique possible minimizer over D).
- Prereg sha256 (as committed and as byte-copied into this dir):
  - `754b66d`: f15596b94d1d2cb94f133743b57fca0cb5f7d4a9975d93848bfab017ca009aba
  - `68a65d9`: f681bc61a738470fe64f89f32521d6608896d937df57b18364b08a7212129209
- No derived-object file was touched before its commit. No B&B run, no LP,
  no slope pass, no subdivision — a pure interval pass (v11 core,
  hash-verified at startup: interval_core `ce70f959…`, vxxz24_float
  `31657c92…` — both asserted equal to the frozen parent values in-run).

## Attempt history (per-attempt logs under DISTINCT filenames, 2026-09-01
## freeze-hygiene rule; both attempt files are byte-copied here)

- attempt 1 (`scratch/branchpin_attempt1.py`; stdout/stderr kept at
  `scratch/branchpin_attempt1_stdout.log`, `…_stderr.log`): **ABORTED at
  control C2 with exit 1** — the registered C2 expectation
  ("eps0 == 0 exactly under drift") was WRONG about the hull semantics:
  `Ival.pos()` of a straddling interval [-a, b] is [0, max(a,b)], and the
  drifted residual set is strictly negative, so pos() returns the
  |residual| hull, not zero. Traceback preserved in
  `branchpin_attempt1_stderr.log` (copied here). No verdict was produced
  by attempt 1 (abort predates P).
- attempt 2 (`branchpin.py`, identical bytes in this dir; source copy
  `scratch/branchpin_attempt2.py`): completed EXIT=0, produced the FALSIFIED
  verdict below. The C2 gate was re-derived from the v11 semantics BEFORE
  attempt 2 ran (documented here; no engine code was touched between
  attempts — only the runner's registered expectation).
- Authoritative re-run of `branchpin.py` inside THIS directory: EXIT=0,
  all decision fields byte-identical to the scratch run (deterministic;
  `branchpin_authoritative_stderr.log` empty as expected).

## GATE-Z (engine equivalence) — PASS

- Max |R-branch deviation| vs frozen `stage_b_rung2_with_lemma1.json`:
  **3.33e-16** (gate 1e-12); vs `…without_lemma1.json`: **3.33e-16**.
- Omega_raw deviation: with **8.88e-16** (2.3715538358350816 vs frozen
  2.3715538358350807), without **8.88e-16** (2.3715518061863823 vs
  2.3715518061863814) — both ≤ 1e-6.
- params sha256 `df75ae3a…da93d` verified in-run.

## C1 (identity-accept, radius 0) — PASS

- `cand1.lo - cand0.lo` (frozen CandW base record) = **-1.4334718842490268e-06**
  vs registered `-2*prop_0*eps_0` = **-1.4334724628081230e-06**:
  identity gap **5.79e-13** ≤ 1e-12.
- Live-in-run candidates reproduce the frozen CandidateWitness base triple
  within 6e-15 (candidate_0 low matches the frozen 0.4968507617380737 to
  3e-16).
- Mechanism reading CONFIRMED: the Lemma-1 charge IS the pd sentence that
  flips the R_glob[0] minimizer to the penalty candidate at p* (the
  frozen with/without records differ by exactly 2*prop_0*eps_0 at
  R_glob[0]).

## C3 (wrong-index reject) — PASS

- Region-1 dual lookup returns eps_1 = **1.143325024764798e-06**, equal to
  the frozen with-lemma1 `glob_r1` record bit-for-bit and distinct from
  eps_0 (>1e-9 apart). The machinery reads the correct lam blocks.

## C2 (dual-drift reject, +1e-3 broadcast on every glob-0 lam entry) — PASS

- Drifted eps boundary: **0.004002150208693772** = eps_0 + 4*drift +
  4.4e-16 (coupling gap 4.41e-16 ≤ 1e-12): the boundary moves by EXACTLY
  4*drift (1 lam_sum + 3 lam_margin per shape), i.e. the eps term is
  FULLY dual-coupled — it responds to a planted dual perturbation at the
  analytic rate.
- Behavioural response: R_glob[0] moves by **-2.666666666666373e-03** =
  predicted -2*prop*4*drift = -2.666666666666667e-03 (gap 2.9e-16): the
  pc-side branch identity responds linearly to the dual drift, confirming
  the eps-coupling is LIVE in the aggregation (not a dead codepath).

## P (the pin over the full 1e-7 box) — **FALSIFIED**, margin quantified

One interval pass, dist[0] leaves = exact-Arb `IC.const(center_i) ±
IC.const(radius)`, everything else point constants:

| quantity | value (outward) |
|---|---|
| G1 = cand_1.lo - cand_0.hi | **[-2.9602186479446094e-05, +2.6689872684820378e-05]** |
| G2 = cand_2.lo - cand_0.hi | **[-2.9602189567791676e-05, +2.668987624020422e-05]** |
| hull span cand_0 (pc) | 3.2123059392635975e-05 (parent's stage-1.2: 3.2123059392596964e-05; dev 3.9e-17) |
| hull span cand_1/2 (pen) | 2.4168999771630385e-05 (parent: 2.4168999771589927e-05; dev 4.0e-17) |
| eps_0 over box | 2.1502086942121845e-06 — **BIT-IDENTICAL to the radius-0 value** |
| R_glob[0] box low | 0.4968353851809626 — **EQUAL to the parent's candidate_lows[0] exactly** |
| Om_raw over box | 2.371563781108452 (raw, uncertified link — reference only) |
| M_low over box | 2.094251842455663 |

Verdict per Amendment 1: G1 and G2 BOTH fail the strictly-positive test
with an interval straddling zero. **FALSIFIED** — the analytic pin of the
R_glob[0] branch identity via the eps-constancy mechanism CANNOT hold at
radius 1e-7: the p_comp candidate's interval top intermittently exceeds
the penalty candidates' interval lows by up to 2.96e-05 across the box.
The failing terms (named per the pre-registration FAIL boundary):
- `G1` and `G2` — the ORDER of the pc vs pen candidates is NOT fixed over
  D; the pc candidate is NOT the unique possible minimizer.
- The pin instrument's `eps_0` term is exonerated (bit-constant over D,
  and FULLY dual-coupled per C2) — the failure is carried by the relative
  BOX MOTION of `pc_0` (~15.4e-6 down) vs `pen_0` (~11.4e-6 down), i.e.
  the p_comp functional of the moving dist[0] distribution, not by the
  Lemma-1 residual.

## Cross-run reconstruction checks (free; frozen-number arithmetic)

- Box candidate spans reproduce the parent BnBGate stage-1.2 floats to
  3.9e-17/4.0e-17 and the R_glob[0] box low EXACTLY (0.4968353851809626).
- eps_0 over the box BIT-EQUALS the radius-0 eps_0 (both
  2.1502086942121845e-06) — the constancy claim is machine-verified, not
  inferred.

## What this freezes

1. **Certified local statement (fail-fast negative):** over the FULL
   un-subdivided D+ (the radius-1e-7 axis-aligned box on the region-0 glob
   dist[0] block, a SUPERSET of the kernel slice D with A·delta = 0 — the
   pin covered the superset box), the R_glob[0] branch minimizer identity
   is NOT pinnable
   by the eps-constancy analytic route: both order inequalities carry a
   −2.96e-05..+2.67e-05 straddle. This is the analytic-route successor of
   the BnBGate quantified obstruction: the branch-identity hull cost on
   R_glob[0] cannot be analytically collapsed by freezing its eps term.
2. **Control-passing instrument:** all four controls (GATE-Z, C1, C2, C3)
   reproduce frozen anchors to ≤8.9e-16 and exhibit the EXACTLY-predicted
   linear responses (4*drift boundary shift; -2*prop*4*drift branch
   response; correct-region dual lookup). The instrument accepts the
   known-true plant and rejects both planted perturbations.
3. **Structural finding:** `eps_0` is RADIUS-INVARIANT over D (bit-exact)
   — any future box analysis can treat the glob Lemma-1 residual as a
   frozen constant, removing it from the movers list; the movers that
   remain are `pc_0`, `nb_*` and (d_box) entropy terms. This narrows the
   parent obstruction's attribution: the R_glob[0] hull straddle is a
   p_comp-vs-penalty race, not an eps artifact.
4. The 21-dimensional kernel question: **UNCHANGED and OPEN** (this
   campaign made NO decision over D's kernel; it measured one branch's
   hull inside D's superset box).

## Rule 7 sentence (mandatory)

The measured set of THIS campaign: the R_glob[0] branch candidate-order
inequalities over D⁺ = the axis-aligned 1e-7 box on the region-0 glob
dist[0] block of the VXXZ24 K100_2.37155181 released vector (D⁺ ⊇ D =
{delta : A·delta = 0 exactly, |delta|_inf ≤ 1e-7}; A = frozen 27x45 0/1
margin matrix, rank 24 over Q; interval core v11 `ce70f959…`,
vxxz24_float `31657c92…`; transcribed 3-region single-p_comp program at
max_level 3, q = 5; every other parameter block at p* point values,
radius exactly binary64 1e-7). No point outside D⁺, no other branch of
R, no other region, no other parameter block, no other radius, no
subdivision, no Monte Carlo, no LP, no slope/gradient pass. The
published 2.37155181 and the record 2.371177 were not touched and are
unreachable via any outcome of this campaign (12.808x guardrail); no
record claim anywhere. Gate C stays PARTIAL.

## Contents

- pre_statement.md — byte-identical copy of the committed prereg
  (commits 754b66d + amend 68a65d9; sha256 f681bc61…9209)
- branchpin.py — the attempt-2 runner (as run, byte-identical to
  scratch/branchpin_attempt2.py; sha256 d7aca589c5dafbf7…)
- branchpin_checkpoint.json — the authoritative in-dir run's checkpoint
  (decision fields byte-identical to the scratch run; sha256 d0b0c260…)
- branchpin_attempt1_stderr.log — attempt-1's preserved traceback (C2
  expectation corrected before attempt 2; no verdict produced)
- branchpin_authoritative_stderr.log — empty (clean final run)
- checksums_dependencies.txt — dependency/workspace sha256s recorded at run
