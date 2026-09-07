# N3 pre-statement — gate `n3-technology-switch`: exact complex-relaxation probing
# (Q(i) CP-witness construction + Koszul/mode flattening battery + ordered-field
#  demotion of the elimination route)

Campaign: to be minted by `python3 scripts/campaign.py init --gate
n3-technology-switch --prereg n3_prereg.md` from cs/oct-rank AFTER this file
is committed (path-scoped commit, `git add <this file>` + `git commit --only
-- <this file>`). Owner agent: `OctRankSwitch`. Date of prereg: 2026-09-03.

## 0. Excluded technology validation (disclosure, per Main's steering of 2026-09-03)

Two /tmp probe runs of one file were executed BEFORE this commit and are
**excluded technology validation** only:

- `/tmp/n3probe/probe.py`, sha256
  `f2cc548dd2e5baa567e12dda36ff293e4a4f0a7b78fc8c6cb38ca84f98b4ce60`.
  Run 1 crashed on an internal assertion (probe's own `cj` copy had dropped
  the frozen len-1 branch — a probe code defect, fixed to a verbatim copy of
  the frozen `cd_ok`/`cj` recursion). Run 2 completed.
- The probes were used ONLY to find defects in the probe's own code
  (conjugation recursion, block-extraction form, eigenvector-basis
  fragility) and to confirm the frozen table conventions. **No scientific
  knob of this campaign — technology claim, witness, bound, instrument
  list, threshold, or stopping rule — was chosen, tuned, or validated by a
  probe result.** Every such knob below is fixed by the mathematical
  construction or by frozen prior evidence, and every decisive
  substitution/bound is re-run from scratch inside the campaign in exact
  arithmetic. A drafted second probe file was never written to disk nor
  executed. No campaign file existed before init; nothing was written into
  `cs/oct-rank/` pre-commit except this prereg.

## 1. Technology choice and justification (from frozen cost evidence)

**Technology: exact complex-relaxation probing.** Three exact, purely
linear-algebraic instruments over Q and Q(i), all deterministic, none
search-based:

- **C1 (witness construction):** an explicit exact Q(i)-rational CP
  decomposition of the certification target TF with as few terms as the
  construction yields (the M2(C)-module regrouping below gives exactly 12
  terms — a theorem of the construction, not a tuned parameter), verified
  by direct substitution into the exact integer tensor entries.
- **C2 (lower-bound battery):** exact fmpq ranks of a FIXED list of nine
  mode-flattening / Koszul-flattening linear maps of TF, giving a certified
  lower bound for rank_C(TF) via the Koszul lemma (rank(K_T) <= 2R, proved
  in section 4).
- **C3 (ordered-field demotion):** exact Q(i)-points of the N2 raw/E/E1
  systems, certifying those ideals are PROPER over Q, plus a pure-integer
  sizing table proving the naive real-Nullstellensatz/SOS route is out of
  scale at every useful degree.

**Why this technology and not the frozen ones** (frozen cost evidence):

1. N2 (FROZEN-NEGATIVE, `20260901T115801Z_8d6683fb_39c752e7eadd`): F4
   elimination died inside degree 9 on a 66-variable calibration system
   (8,080,476 x 74,499,764 Macaulay matrix, 22.65 GB, degree 8 alone
   3802 s wall / 3593 s CPU), 4x smaller than the 273-variable main
   systems. Elimination is out of reach by orders of magnitude — AND it is
   structurally mispointed: msolve Q-mode certifies {1} over Q only if the
   elimination ideal is the unit ideal, i.e. only if the system has NO
   complex point. This campaign certifies the converse fact directly
   (C3): the rank-13 CP system of TF HAS complex points, so the {1}-over-Q
   verdict was impossible ab initio and the only possible elimination
   outputs were mod-p artifacts, which N2's own prereg (section 3, clause
   (ii)) already admits cannot be airtight. The elimination technology has
   **no path to a rank-14 verdict at any budget**; that is a theorem-shaped
   route closure, not a budget observation. This campaign never invokes
   msolve, F4, or any Groebner engine.
2. N1 (FROZEN-INCONCLUSIVE, `20260901T103338Z_06de66a3_d0b7bc2dc6a5`):
   40-seed TRF + exact-Newton admission: 0/40 admitted, exact residuals
   1.4e-5..7.8e-4, four orders above the exact gate. The real-feasibility
   side is search-hostile. This campaign never invokes TRF, seeds, or
   numerical admission. **No TRF and no F4 code path is rerun anywhere.**
3. The complex side is the unattacked flank: rank_C(TF) <= rank_R(TF) <= 14
   (the frozen real 14-term upper bound is in particular a complex
   decomposition), so complex 13-witnesses exist in principle; the
   construction below produces one EXACTLY, by exact linear algebra, with
   zero search. Exact witnesses are then self-certifying by substitution
   (rule 14: the instrument can establish the claim, not merely be
   consistent with it).

This is candidate technology family (b) of the assignment ("an exact
infeasibility/feasibility certificate that avoids Groebner bases entirely"),
instantiated as: feasibility certificates over C (witnesses) + flattening
lower bounds + an ordered-field demotion proving which future certificates
must be real-specific.

## 2. Object, conventions, and frozen anchors

- **Object:** TF = blockdiag(tau, tau) as a 3x8x8 integer tensor, where
  tau = (I_4, L_i, L_j) is the quaternion 3-slice tensor on the
  Cayley-Dickson basis (1, i, j, k), product
  (a,b)(c,d) = (ac - conj(d)b, da + b conj(c)), table entry
  T[p][b][c] = (e_p * e_b)_c. The instrument carries a **verbatim semantic
  copy** of the frozen `cd_ok`/`cj` recursion from
  `campaigns/20260901T115801Z_8d6683fb_39c752e7eadd/n2_systems.py`
  (sha256 `52f6e65111d1fde64a96cae45c06704cfe0aad9ac30fc4f49576df6ca1f43c69`),
  including the len-1 conjugation branch.
- **Frozen rank-equivalence (consumed, not re-derived):** rank(TF) =
  rank((L_1, L_i, L_j)) via the two frozen entrywise diagonal-similitude
  identities (Route-F campaign `2026-09-01T04:30:00Z_routeF_kraw`).
- **Frozen certified facts consumed:** floor rank_R(TF) >= 13 (1 peel +
  pencil 12; S3/C4/RouteAF chain; Route A family cap 12 < 13 intact);
  upper rank_R(TF) <= 14 (blockwise 7+7); rank_R(tau) = 7; the published
  window 18 <= R_R(T_O) <= 25 (untouched by everything here).
- **Byte-level object anchors (asserted in-run, Stage A):**
  (a) tau has exactly 12 nonzero entries, TF exactly 24;
  (b) parsing the 48 polynomials of the frozen `s1a_tau_r7.txt` (sha256
  `16fe2bc45fe5813567722934e3acc11e7615af1efed29ff6778a907d79af2530`)
  reproduces -tau entrywise as constant terms;
  (c) parsing the first 192 polynomials of the frozen
  `s2e_tf_r13_ext.txt` (sha256
  `117fe71e08f335c268a49358ea14161c7193a9623d6c1c1d255eca97131d33f6`)
  reproduces -TF entrywise;
  (d) concat-vs-Kronecker guard: blockdiag(M, M) == M (x) I_2 verified
  entrywise on a fixed exact plant M (the 4x4 matrix with M[j][k] =
  (j*k) mod 5 - 2, j,k in 0..3);
  (e) dimension plants: every verifier refuses tensors of mismatched
  shape (3x8x8 vs 3x4x4 vs 8x8) with a named exception.

## 3. Arithmetic backend

Python 3.14.3 from `/Users/jinleic/jinleic-workspace/cs/.venv/bin/python`
(Clang 17.0.0 build). The instrument uses **only the standard library**
(`fractions`, `itertools`, `json`, `hashlib`, `resource`, `time`, `sys`).
Exact arithmetic: rational pairs (a, b) representing a + b*i with
`fractions.Fraction`; exact Gaussian elimination for Q(i) linear systems;
exact RREF for Q ranks. **No floats, no numpy, no BLAS in any
claim-relevant path.** Floats appear nowhere.

## 4. The fixed construction (theorems T1-T6, no free parameters)

All objects over Q(i) unless stated. E_conj = C^2_doub (x) C^2_std (x)
C^2_mult with index (d, s, m) -> 4d + 2s + m.

- **T1 (regrouping, exact linear solve).** Let
  X_i = [[0,1],[-1,0]], X_j = [[0,-i],[-i,0]] in M_2(Q(i)), and
  K_i = X_i (x) I_2, K_j = X_j (x) I_2 in M_4(Q(i)). Solve the linear
  system in the 16 entries of S0 in M_4(Q(i)):
  S0 L_i = K_i S0 and S0 L_j = K_j S0  (32 linear equations).
  By the module theorem (C (x)_R H_L ≅ M_2(C) acting on C^4 ≅ C^2 (x) C^2
  as standard (x) trivial, and the presentation
  <x, y | x^2 = y^2 = -1, xy + yx = 0> ≅ M_2(C)), the solution space has
  dimension exactly 4 and consists of the isomorphisms Hom_{M_2(C)}-space.
  Fixed selection rule (no search): take the nullspace basis in the order
  produced by exact Gauss elimination with first-nonzero-pivot rule;
  S0 := [n1 n2 n3 n4]. If det(S0) = 0, the single deterministic fallback
  is [n4 n3 n2 n1]; if that is also singular: verdict INCONCLUSIVE
  (construction anomaly), no retries. In-run asserts: nullity = 4 exactly
  (a different nullity is a recorded anomaly -> INCONCLUSIVE), det != 0,
  and BOTH conjugation identities entrywise.
- **T2 (rank-1 trinomials).** In M_2(Q(i)): p1 = I2 + i*X_i,
  p2 = I2 - i*X_i, p3 = I2 + i*X_j. Each satisfies
  det = x^2 + y^2 + z^2 evaluated at (x,y,z) in {(1,1,0),(1,-1,0),(1,0,1)}
  = 0 (the general identity det(xI + yX_i + zX_j) = x^2 + y^2 + z^2 holds
  because X_i^2 = X_j^2 = -I2 and X_iX_j = -X_jX_i), hence each p_r is
  rank <= 1; in-run factors p_r = u_r v_r^T exactly (first nonzero entry
  rule) and asserts the factorization entrywise.
- **T3 (alpha solve).** Solve the 12x3 exact linear system
  [vec(p1) vec(p2) vec(p3)] alpha_col = vec(slice) for the three W-slices
  I2, X_i, X_j; in-run asserts the residual is identically zero and that
  the three slices of W = (I2, X_i, X_j) lie in span{p1, p2, p3}
  (span check: the 3x3 coefficient matrix is exactly invertible — asserted
  by the unique-solution residual zero).
- **T4 (12-term witness).** In the conjugated frame with slices
  I_8, K_i (x)-lifted as I2_doub (x) X_i (x) I2_mult etc., the identity
  I_D (x) (u_r v_r^T) (x) I_M = sum_{d=0,1} sum_{m=0,1}
  (d-vec (x) u_r (x) m-vec)(d-vec (x) v_r (x) m-vec)^T gives exactly
  3*2*2 = 12 rank-1 terms. Pull back each factor by S0^{-1} to the
  original frame (v -> S0^{-1} v applied to both matrix factors). In-run:
  **substitution-verify all 3*64 = 192 entries of TF against the exact
  integer table** (conjugated-frame check first, then original frame).
  A single mismatch -> verdict INCONCLUSIVE (witness failure), full
  residual table frozen.
- **T5 (13-term split, raw-system point).** Split term r=0 (fixed):
  (alpha_r, P, Q) -> (alpha_r, P, Q + W) + (-alpha_r, P, W) with the fixed
  W = e_3 (0-based coordinate 3) of the conjugated 8-space; deterministic
  fallback order r = 1 then r = 2 if term 0's gauge coordinates fail
  T6's nonsingularity; no other fallback. Verifies 192/192 again and that
  all 13 term triples are nonzero. This is an exact Q(i) point of the RAW
  192-equation/247-unknown rank-13 CP system of TF.
- **T6 (gauge coordinates + E/E1 points).** For each of the 13 terms of
  T5's point, (a[0,s], c[0,s]) != (0,0) is asserted exactly; the 2x2
  gauge system a[0,s] d_s - c[0,s] g_s = 1, c[0,s] d_s + a[0,s] g_s = 0
  is solved exactly per term. The E-point (all 218 equations of the frozen
  `s2e_tf_r13_ext.txt` evaluated == 0) and the E1-point (all 220 equations
  of the frozen `s2e1_tf_r13_ext_pin.txt` evaluated == 0, using the
  pinned zero-term point: term index 0 = (0,0,0), terms 1..12 = the T4
  12-term witness) are both evaluated by an exact parser of the frozen
  files' restricted grammar (unit-coefficient monomials and integer
  constants), with the parser self-checked by stage A's constant
  extraction.
- **Consequence (the demotion theorem, proved in VERDICT.md from T5/T6):**
  a Q(i)-point of a Q-polynomial system implies its ideal is proper over Q
  (a unit ideal has an empty complex variety, by the Nullstellensatz).
  Hence the raw rank-13 CP system of TF, and the N2 gauge systems E and
  E1, all have proper elimination ideals over Q: **N2's clause-(i)
  ({1}-over-Q) verdict was impossible ab initio**; any future mod-p {1}
  outcome is a bad-reduction artifact carrying no char-0 content (this
  sharpens, and does not contradict, N2 prereg section 3); and any future
  rank >= 14 certificate MUST be real/order-specific (Krawczyk-type
  exclusion, SOS with real structure, or new mathematics) — elimination
  over Q is closed as a route, structurally, not just by budget.

## 5. The fixed lower-bound battery (nine instruments, no shopping)

LB battery, all exact fmpq ranks, all on the integer TF:

1. rank of mode-1 flattening (3 x 64); LB = that rank.
2. rank of mode-2 flattening (8 x 24); LB = that rank.
3. rank of mode-3 flattening (8 x 24); LB = that rank.
4. Koszul p=1, A = mode-1 (3), B = mode-2 (8), C = mode-3 (8): the map
   A (x) B* -> Lambda^2 A (x) C*, alpha (x) beta |->
   sum_i (alpha ^ e_i) (x) <beta, T(e_i, -, -)>; LB = floor(rank/2) from
   the **Koszul lemma**: if T = sum_{r=1}^R a_r (x) b_r (x) c_r then
   rank(K_T) <= sum_r rank(kappa_{a_r}) * rank(ell_r) <= 2R, since
   kappa_a: alpha |-> alpha ^ a has rank 2 (kernel span{a}) and
   ell_r: beta |-> beta(b_r) c_r has rank 1. (Lemma proved here; used for
   every Koszul instrument.)
5. Koszul p=1, A = mode-2 (8), B = mode-3 (8), C = mode-1 (3):
   kappa_b has rank 7 (dim Lambda^1 - dim Lambda^0 of span{b}), so
   LB = floor(rank/7).
6. Koszul p=1, A = mode-3 (8), B = mode-1 (3), C = mode-2 (8): LB =
   floor(rank/7).
7. Koszul p=2, A = mode-1 (3): map Lambda^2 A (x) B* -> Lambda^3 A (x) C*;
   summand rank <= 1 (kappa: Lambda^2 -> Lambda^3 in 3 dims, alpha |->
   alpha ^ a, rank 1); LB = rank.
8. Koszul p=2, A = mode-2 (8): summand rank <= C(7,1) = 7; LB =
   floor(rank/7).
9. Koszul p=2, A = mode-3 (8): LB = floor(rank/7).

LB := max over the nine. rank_C(TF) is reported as the exact interval
[LB, 12]; the statement "rank_C(TF) = 12 EXACTLY" is made **iff** LB = 12.
No other rank statement is made. If LB < 12 the interval is the honest
result and the verdict says so.

**Battery plants (both directions, run before the battery):**
- P-KOSZ (accept): T = sum_{r=1}^{3} e_r (x) f_r (x) g_r on R^3 (x) R^4 (x)
  R^4 (basis-indexed diagonal). Closed form by the lemma's equality
  analysis: the three kappa-images are the pairwise-disjoint 2-dim spaces
  span{e_r ^ e_s : s != r}, so instrument #4 must return exactly
  rank = 6. Any other value -> INVALID INSTRUMENT.
- P-KOSZ-R (reject): same plant with slice 3 identically zero: closed
  form rank = 4 (two disjoint images). Must return exactly 4 and differ
  from P-KOSZ.

## 6. Controls (all through the SAME code path, run BEFORE the construction)

- **A1 (accept known-valid):** the frozen exact T_C rank-3 witness and
  **A2** the frozen exact T_H rank-8 witness (both entrywise fmpq, from
  the gate-C campaign `2026-08-30T01:40:00ZZ_gateC_tau_sharpness`, re-built
  inside the instrument from the frozen construction code semantics with
  the frozen file hashes recorded) must PASS substitution verification on
  their own tensors T_C = (I_2, [[0,-1],[1,0]]) and T_H (4x4x4).
- **R1 (reject one-coordinate perturbation of the target):** T_cor = TF
  with entry [0][0][0] += 1. BOTH the T4 12-term witness and the T5
  13-term point must FAIL on T_cor (>= 1 mismatched entry, mismatch count
  recorded exactly; a PASS here = FALSE CERTIFICATE = quarantine +
  INVALID INSTRUMENT).
- **R2 (reject corrupted accept-control):** T_C with entry [0][0][0] += 1
  must FAIL the A1 witness.
- **G-DIM / G-KRON:** section 2(d),(e) dimension plants and Kronecker
  guard.

## 7. C3 sizing table (route negative, pure integer arithmetic)

For the fixed degree list d in {3, 4, 5, 6}, for a degree-d real
Nullstellensatz certificate on the raw 192-cubic system: report
M(d) = C(247 + d, d) (monomials of degree <= d in 247 variables),
N(d-3) = C(247 + (d-3), d-3) (per-multiplier coefficient block), the
dense coefficient-matrix RAM M(d) x (192*N(d-3)) x 8 bytes, and the
SOS-Gram RAM C(247 + floor(d/2), floor(d/2))^2 x 8 bytes. Compute exactly;
compare against the workstation RAM (reported via `sysconf`); freeze the
table. The claim frozen is ONLY: every degree >= 4 naive Ncert/SOS
certificate is out of scale by the tabulated factors on this workstation;
the degree-3 constant-multiplier family has no existence guarantee and is
recorded as impractical-without-new-structure. This is a sizing negative
about a route family, never a rank claim.

## 8. Resource cap, accounting, and stopping rules

- Total campaign CPU cap: 2 h (`time.process_time` printed at every stage
  boundary; hard self-check aborts with verdict INCONCLUSIVE-BUDGET if
  exceeded). Expected: minutes (all stages are polynomial-size exact
  linear algebra; the largest object is a 24x24 exact RREF and ~600 small
  substitutions).
- `PYTHONDONTWRITEBYTECODE=1` in the environment; `sys.dont_write_bytecode
  = True` before any import; `resource.setrlimit(RLIMIT_CPU, soft -> hard)`
  raised at the top of the instrument; `nice` increment 10 at launch;
  single process, no threads, no daemons, no blocking waits (the
  instrument is a single foreground run; wall-clock is bounded by the CPU
  cap).
- Stopping rules: every stage runs once; a failure records its full
  history and moves the verdict per section 9; **no parameter is ever
  adjusted after seeing a result** (the only fallbacks are the two
  deterministic branches fixed in T1 and T5).
- Reproduction: run command
  `PYTHONDONTWRITEBYTECODE=1 nice -n 10 ../../.venv/bin/python
  n3_instrument.py` from the run dir (venv path absolute:
  `/Users/jinleic/jinleic-workspace/cs/.venv/bin/python`); the instrument
  reads the three frozen N2 files listed in section 2 by RELATIVE path
  from the run dir
  (`../20260901T115801Z_8d6683fb_39c752e7eadd/...`) and records their
  sha256 in-run; environment: python 3.14.3, stdlib only.

## 9. Adjudication vocabulary (fixed now)

- **FROZEN-CERTIFIED**: stages A (anchors), C (controls A1,A2,R1,R2,
  plants, guards), T1-T6, battery+plants, and sizing all complete with
  their expected outcomes. Claims allowed, exactly:
  (i) rank_C(TF) in [LB, 12] with both ends exact and machine-verified
  (equality claimed iff LB = 12);
  (ii) rank_R(TF) >= 13 > 12 >= rank_C(TF): a certified strict
  real-over-complex rank gap for this tensor (floor consumed from the
  frozen chain; upper end from T4);
  (iii) the demotion theorem of T6 (proper ideals; N2 clause-(i) void;
  mod-p {1} = bad-reduction artifact only; upward route requires
  real-specific certificates);
  (iv) the C3 sizing negative.
  **The {13,14} frontier does NOT move** under any outcome: no rank_R
  claim of any kind is made in either direction; the verdict must state
  verbatim: "Absence of a 13-witness is NOT evidence for 14; absence of an
  impossibility argument is NOT evidence for 13." The published window
  18 <= R_R(T_O) <= 25 is untouched; no published work is refuted.
- **FROZEN-INCONCLUSIVE (INVALID INSTRUMENT)**: any control break
  (A1/A2 fail; R1/R2 pass; plant value off; guard failure; T1 nullity
  anomaly). Nothing reported as evidence.
- **FROZEN-INCONCLUSIVE (CONSTRUCTION ANOMALY)**: T1-T6 residual
  nonzero / gauge singularity with both fixed fallbacks exhausted.
- **FROZEN-INCONCLUSIVE (BUDGET/CRASHED)**: cap abort or script failure,
  with partial history frozen.
- status.json mapping: FROZEN-CERTIFIED -> FROZEN-CERTIFIED;
  INVALID INSTRUMENT / CONSTRUCTION ANOMALY / BUDGET / CRASHED ->
  FROZEN-INCONCLUSIVE with the label in the verdict.

## 10. Rule-7 prospective scope

Covered: TF = blockdiag(tau, tau) (with the frozen rank-equivalence to
(L_1, L_i, L_j)) at rank exactly 13's complex feasibility side; the nine
fixed flattening instruments and two plants; the fixed construction
T1-T6; the three frozen N2 system files; the C3 sizing table at degrees
3..6. NOT covered, hence decided by nothing here: rank_R(TF) in either
direction ({13,14} stays OPEN with floor 13 and upper 14 intact); any
real-specific certificate question; any rank other than 13's complex
analogue; T_O or any other slice family; border rank; tau's complex rank
(named candidate for a FOLLOW-UP gate, not committed here); the frozen
N1/N2 campaigns (consumed read-only); R_R(T_O).

## 11. Lifecycle and freeze plan

Commit this file (path-scoped) -> `campaign.py init --gate
n3-technology-switch --prereg n3_prereg.md` from cs/oct-rank ->
byte-identical copy of this file into the minted run dir +
PROVENANCE.md (source commit hash + prereg sha256, copied AFTER init,
BEFORE any in-run compute) -> copy the probe disclosure file
(/tmp/n3probe/probe.py + its sha256 + the two-run transcript summary of
section 0) into the run dir as `probe_disclosure/` -> write and run
`n3_instrument.py` (sections 2-8, single pass) -> freeze (`campaign.py
freeze`) with instrument, results JSON, witness JSON (exact Q(i)
coordinates), verdict notes -> `campaign.py close` with the section-9
mapping -> append the target README current-state section. Per Main's
steering: if this gate lands FROZEN-CERTIFIED, the next implied frontier
gate (candidates: exact rank_C(tau) via the same technology — the
tower-gap pattern with T_C (3 real / 2 complex) and tau as rungs; or a
witness-structured exact real-13 fusion attempt) will be preregistered
and run with fresh controls immediately after, in this same session.
