# N3 pre-statement v2 — gate `n3-technology-switch`: exact complex-relaxation
# probing (Q(i) CP-witness construction + Koszul/mode flattening battery +
# ordered-field demotion of the elimination route)

Campaign: to be minted by `python3 scripts/campaign.py init --gate
n3-technology-switch --prereg n3_prereg_v2.md` from cs/oct-rank AFTER this
file is committed (path-scoped). Owner agent: `OctRankSwitch`. Date:
2026-09-04.

**This v2 supersedes the invalid prereg of run
20260904T021928Z_66f9d4ee_138fcf06a676 (FROZEN-INCONCLUSIVE, INVALID-PREREG,
closed before any control/compute; sha256
29e08a181d5409777380c9f1c846f4cdc44d71efea83ef2ac34acba365194846, commit
b4239c9).** Corrections carried by v2, preserved in the closed run's
CLOSURE_NOTE.md: (1) Koszul p=2 divisors for dim A = 8 are C(7,2) = 21, not
C(7,1) = 7 — with 7, instruments 8/9 would report invalid lower bounds;
generally rank(kappa_a: Lambda^p A -> Lambda^(p+1) A, alpha |-> alpha ^ a)
= C(dim(A)-1, p); (2) T4's pullback is asymmetric: from S0 L S0^-1 = K the
term-level pullback of P Q^T is (S0^-1 P)(S0^T Q)^T, i.e. the second factor
pulls back by S0^T; (3) the invalid implication "rank_C <= rank_R <= 14 so
complex 13-witnesses exist in principle" is removed — complex 13-feasibility
is established ONLY by the explicit 12-term construction plus its split,
which is exactly what this campaign certifies; (4) T2's coordinates stated
correctly: the rank-1 test elements are xI + yX_i + zX_j at
(x,y,z) in {(1,i,0),(1,-i,0),(1,0,i)} where x^2+y^2+z^2 = 0 over Q(i).

## 0. Excluded technology validation (disclosure, per Main's steering)

One /tmp probe file ran twice pre-commit of the superseded prereg:
`/tmp/n3probe/probe.py` (sha256
`f2cc548dd2e5baa567e12dda36ff293e4a4f0a7b78fc8c6cb38ca84f98b4ce60`; run 1
crashed on the probe's own `cj` defect, run 2 completed), copied to the
closed run's `probe_disclosure/`. It is **excluded technology validation**:
it located code defects in the probe itself (conjugation recursion,
block-extraction form, eigenvector-basis fragility) and confirmed the
frozen table conventions. **No scientific knob of this campaign —
technology, witness, bound, instrument list, divisor, threshold, or
stopping rule — was chosen, tuned, or validated by any probe result**, and
no probe informed the v1 instruments that turned out defective (the two v1
defects above were found by Main's review of the prereg text itself). Every
knob below is fixed by mathematics or frozen prior evidence; every decisive
substitution/bound is rerun from scratch inside this campaign in exact
arithmetic.

## 1. Technology choice and justification (from frozen cost evidence)

**Technology: exact complex-relaxation probing.** Three exact, purely
linear-algebraic instruments over Q and Q(i), all deterministic, none
search-based:

- **C1 (witness construction):** an explicit exact Q(i)-rational CP
  decomposition of TF with the term count the construction yields (the
  M2(C)-module regrouping gives exactly 12 terms — a theorem of the
  construction, not a tuned parameter), verified by direct substitution
  into the exact integer tensor entries, then padded by an exact split to
  an exact 13-term point of the raw rank-13 CP system. **Complex
  13-feasibility of TF is an OUTPUT of this campaign, established only by
  these verified constructions** (the frozen real 14-term bound gives a
  complex 14-term decomposition and implies nothing at 13).
- **C2 (lower-bound battery):** exact fmpq ranks of a FIXED list of nine
  mode-flattening / Koszul-flattening linear maps of TF, giving a
  certified lower bound for rank_C(TF) via the Koszul lemma with the
  CORRECT divisors (section 5).
- **C3 (ordered-field demotion):** exact Q(i)-points of the raw N2 rank-13
  CP system and (conditionally, see T6) of the gauge system E, certifying
  those ideals are PROPER over Q; plus detection of any self-contradiction
  in the serialized E1 (pre-registered in T6); plus a pure-integer sizing
  table proving the naive real-Nullstellensatz/SOS route is out of scale
  at every useful degree.

**Why this technology and not the frozen ones:**

1. N2 (FROZEN-NEGATIVE, `20260901T115801Z_8d6683fb_39c752e7eadd`): F4
   elimination died inside degree 9 on a 66-variable calibration system
   (8,080,476 x 74,499,764 Macaulay matrix, 22.65 GB, degree 8 alone
   3802 s wall / 3593 s CPU), 4x smaller than the 273-variable main
   systems. Elimination is out of reach by orders of magnitude — AND it is
   structurally mispointed: msolve Q-mode certifies {1} over Q only if the
   elimination ideal is the unit ideal, i.e. only if the system has NO
   complex point. This campaign certifies the converse fact directly
   (C3): the raw rank-13 CP system of TF HAS an exact complex point, so a
   {1}-over-Q verdict on it is impossible ab initio and the only possible
   elimination outputs are mod-p artifacts, which N2's own prereg
   (section 3, clause (ii)) already admits cannot be airtight. The
   elimination technology has **no path to a rank-14 verdict at any
   budget**; that is a theorem-shaped route closure, not a budget
   observation. This campaign never invokes msolve, F4, or any Groebner
   engine.
2. N1 (FROZEN-INCONCLUSIVE, `20260901T103338Z_06de66a3_d0b7bc2dc6a5`):
   40-seed TRF + exact-Newton admission: 0/40 admitted. The real-feasibility
   side is search-hostile. This campaign never invokes TRF, seeds, or
   numerical admission. **No TRF and no F4 code path is rerun anywhere.**
3. The complex side is the unattacked flank: the only frozen complex-side
   fact is the trivial rank_C(TF) <= rank_R(TF) <= 14. Everything at 13
   and below over C is open input to this campaign and is settled, if at
   all, by exact construction (C1) and exact bounds (C2) — both
   self-certifying by substitution/rank computation (rule 14: the
   instrument can establish the claim, not merely be consistent with it).

This is candidate technology family (b) of the assignment
("an exact certificate route that avoids Groebner bases entirely"),
instantiated as: feasibility witnesses over C + flattening lower bounds +
an ordered-field demotion proving which future certificates must be
real-specific.

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
  entrywise on the fixed plant M[j][k] = (j*k) mod 5 - 2 (j,k in 0..3);
  (e) dimension plants: every verifier refuses tensors of mismatched
  shape (3x8x8 vs 3x4x4 vs 8x8) with a named exception.

## 3. Arithmetic backend

Python 3.14.3 from `/Users/jinleic/jinleic-workspace/cs/.venv/bin/python`.
The instrument uses **only the standard library** (`fractions`, `re`,
`itertools`, `json`, `hashlib`, `math`, `resource`, `time`, `sys`, `os`).
Exact arithmetic: rational pairs (a, b) = a + b*i over
`fractions.Fraction`; exact Gauss elimination/RREF over Q(i); exact RREF
over Q. **No floats, no numpy, no BLAS in any claim-relevant path.**

## 4. The fixed construction (theorems T1-T6, no free parameters)

All objects over Q(i) unless stated. Conjugated frame:
E_conj = C^2_doub (x) C^2_std (x) C^2_mult, index (d, s, m) -> 4d + 2s + m.

- **T1 (regrouping, exact linear solve).** Let
  X_i = [[0,1],[-1,0]], X_j = [[0,-i],[-i,0]] in M_2(Q(i)); both square to
  -I2 and anticommute, so the C-algebra they generate carries the
  presentation <x, y | x^2 = y^2 = -1, xy + yx = 0> ≅ C (x)_R H ≅ M_2(C).
  Set K_i = X_i (x) I_2, K_j = X_j (x) I_2 in M_4(Q(i)) and solve the
  linear system in the 16 entries of S0 in M_4(Q(i)):
  S0 L_i = K_i S0 and S0 L_j = K_j S0  (32 linear equations over Q(i)).
  By the module theorem (the left-multiplication representation of
  C (x)_R H_L on C^4 is the standard M_2(C)-module with multiplicity 2,
  and both generator pairs realize the same presentation), the solution
  space is the intertwiner space of two isomorphic M_2(C)-modules and has
  dimension exactly 4. Fixed selection rule (no search): take the
  nullspace basis in the order produced by exact RREF with
  first-nonzero-pivot rule; S0 := [n1 n2 n3 n4]. If det(S0) = 0, the
  single deterministic fallback is [n4 n3 n2 n1]; if that is also
  singular: verdict INCONCLUSIVE (construction anomaly), no retries.
  In-run asserts: nullity = 4 exactly (any other value -> CONSTRUCTION
  ANOMALY), det != 0, and BOTH conjugation identities entrywise
  (S0 L_i S0^-1 = K_i, S0 L_j S0^-1 = K_j).
- **T2 (rank-1 trinomials).** In M_2(Q(i)): p1 = I2 + i*X_i,
  p2 = I2 - i*X_i, p3 = I2 + i*X_j; these are xI + yX_i + zX_j at
  (x,y,z) = (1,i,0), (1,-i,0), (1,0,i). The identity
  det(xI + yX_i + zX_j) = x^2 + y^2 + z^2 holds because
  X_i^2 = X_j^2 = -I2 and X_iX_j + X_jX_i = 0 (so
  (yX_i + zX_j)^2 = -(y^2 + z^2) I2 and the eigenvalue argument gives the
  determinant), and x^2 + y^2 + z^2 = 0 at all three points, so each p_r
  has rank <= 1; in-run factors p_r = u_r v_r^T exactly (first nonzero
  entry rule) and asserts the factorization entrywise.
- **T3 (alpha solve).** Solve the 12x3 exact linear system
  [vec(p1) vec(p2) vec(p3)] alpha_col = vec(slice) for the three W-slices
  I2, X_i, X_j; in-run asserts the solution is unique (pivot columns
  {0,1,2}) and the residual identically zero, which certifies
  span{p1, p2, p3} contains all three W-slices.
- **T4 (12-term witness, ASYMMETRIC pullback).** In the conjugated frame,
  I_D (x) (u_r v_r^T) (x) I_M
  = sum_{d=0,1} sum_{m=0,1} (d-vec (x) u_r (x) m-vec)
      (d-vec (x) v_r (x) m-vec)^T
  gives exactly 3*2*2 = 12 rank-1 terms summing (with the alpha
  coefficients) to the conjugated-frame tensor with slices
  I2_d (x) I4 (x) I2_m, I2_d (x) K_i (x) I2_m, I2_d (x) K_j (x) I2_m.
  In-run verifies all 3*64 = 192 entries in the conjugated frame. Then
  pull back to the original frame: since L = S0^-1 K S0, a rank-1
  conjugated-frame matrix P Q^T pulls back to
  S0^-1 P Q^T S0 = (S0^-1 P)(S0^T Q)^T — i.e. **b_orig = S0^-1 b_conj and
  c_orig = S0^T c_conj (ASYMMETRIC pullback; the v1 wording "S0^-1 applied
  to both factors" was an algebra slip, corrected at remint)**. Coefficients
  alpha are slice coefficients and are unchanged. In-run verifies all 192
  integer TF entries in the original frame. Any single mismatch -> verdict
  INCONCLUSIVE (witness failure), full residual table frozen.
- **T5 (13-term split, raw-system point).** Split term r=0 (fixed):
  (alpha_r, P, Q) -> (alpha_r, P, Q + W) + (-alpha_r, P, W) with the fixed
  W = e_3 (0-based coordinate 3) of the conjugated 8-space; deterministic
  fallback order r = 1 then r = 2 if term 0 fails T6's gauge conditions;
  no other fallback. Verifies 192/192 substitution again and that all 13
  term triples are nonzero. This is an exact Q(i) point of the RAW
  192-equation/247-unknown rank-13 CP system of TF — the fact the demotion
  theorem consumes.
- **T6 (gauge coordinates + E/E1 handling, pre-registered).** For each of
  the 13 terms, (a[0,s], c[0,s]) != (0,0) is asserted exactly, and the
  2x2 gauge system a[0,s] d_s - c[0,s] g_s = 1, c[0,s] d_s + a[0,s] g_s
  = 0 is solved exactly per term; over Q(i) the determinant is
  a[0,s]^2 + c[0,s]^2, which can vanish for nonzero isotropic pairs
  (unlike over F_p, p = 3 mod 4). If any term is gauge-isotropic, the
  E-point is not certified and this is recorded exactly (the demotion
  theorem still stands on the raw system). The E-point evaluation parses
  ALL 218 equations of the frozen `s2e_tf_r13_ext.txt` and requires every
  one to evaluate to 0. For E1, the pinned point is term 0 = (0,0,0)
  (satisfying the pin equations a0_0 = 0, c0_0 = 0) with terms 1..12 =
  the T4 witness; **pre-registered prediction with exact handling:** the
  serialized E1 (220 equations) may pair the pin with a gauge pair on the
  same term (the v1 generator gauged all 13 terms), in which case the
  ideal of E1-as-serialized is the UNIT IDEAL by the explicit 3-polynomial
  combination (a0_0 = 0, c0_0 = 0, a0_0 d_0 - c0_0 g_0 - 1 = 0 give
  -1 = 0), E1-as-serialized admits no point over ANY field, and it
  carries no rank content (a {1} certificate on it would be vacuous about
  rank). The instrument detects and records this exactly (which equations
  fail at the pinned point, and the explicit combination); it is
  disclosed as a defect of the FROZEN N2 serialization (rule 5), not an
  anomaly of this campaign, and does not affect the FROZEN-CERTIFIED
  path — the demotion theorem for the raw system and E stands on its own.
- **Consequence (the demotion theorem, proved in VERDICT.md from T5):**
  a Q(i)-point of a Q-polynomial system implies its ideal is proper over
  Q (a unit ideal has an empty complex variety, by the Nullstellensatz).
  Hence the raw rank-13 CP system of TF has a proper elimination ideal
  over Q (and E likewise if its T6 point certifies): **N2's clause-(i)
  ({1}-over-Q) verdict was impossible ab initio**; any future mod-p {1}
  outcome is a bad-reduction artifact carrying no char-0 content (this
  sharpens, and does not contradict, N2 prereg section 3); and any future
  rank >= 14 certificate MUST be real/order-specific (Krawczyk-type
  exclusion, SOS with real structure, or new mathematics) — elimination
  over Q is closed as a route, structurally, not just by budget.

## 5. The fixed lower-bound battery (nine instruments, correct divisors)

LB battery, all exact fmpq ranks, all on the integer TF. The Koszul lemma:
if T = sum_{r=1}^R a_r (x) b_r (x) c_r then the Koszul flattening
K_T(alpha (x) beta) := sum_a (alpha ^ e_a) (x) T_a(beta), with T_a(beta)
the C-vector k |-> T[a][j][k] contracted against beta, satisfies
K_T = sum_r (kappa_{a_r} (x) ell_r) where
kappa_{a_r}: alpha |-> alpha ^ a_r has rank C(dim(A)-1, p) on
Lambda^p A (p = 1: dim A - 1; p = 2: C(dim(A)-1, 2)), and
ell_r: beta |-> beta(b_r) c_r has rank 1; hence
rank(K_T) <= R * C(dim(A)-1, p) and R >= floor(rank(K_T) /
C(dim(A)-1, p)). Divisors: p=1 dim A=3: 2; p=1 dim A=8: 7; p=2 dim A=3:
C(2,2)=1; **p=2 dim A=8: C(7,2)=21** (the v1 prereg's 7 was the fatal
proof error, corrected here).

1. rank of mode-1 flattening (3 x 64); LB = that rank.
2. rank of mode-2 flattening (8 x 24); LB = that rank.
3. rank of mode-3 flattening (8 x 24); LB = that rank.
4. Koszul p=1, (A,B,C) = (mode-1, mode-2, mode-3); LB = floor(rank/2).
5. Koszul p=1, (A,B,C) = (mode-2, mode-3, mode-1); LB = floor(rank/7).
6. Koszul p=1, (A,B,C) = (mode-3, mode-1, mode-2); LB = floor(rank/7).
7. Koszul p=2, (A,B,C) = (mode-1, mode-2, mode-3); LB = floor(rank/1).
8. Koszul p=2, (A,B,C) = (mode-2, mode-3, mode-1); LB = floor(rank/21).
9. Koszul p=2, (A,B,C) = (mode-3, mode-1, mode-2); LB = floor(rank/21).

LB := max over the nine. rank_C(TF) is reported as the exact interval
[LB, 12]; the statement "rank_C(TF) = 12 EXACTLY" is made **iff** LB = 12.
No other rank statement is made. If LB < 12 the interval is the honest
result and the verdict says so.

**Battery plants (both directions, run before the battery):**
- P-KOSZ (accept): T = sum_{r=1}^{3} e_r (x) f_r (x) g_r on R^3 (x) R^4 (x)
  R^4 (diagonal). Hand-verified closed form: K(e_p (x) f_j*) =
  (e_p ^ e_j) (x) g_j for j <= 2 and 0 for j = 3 (in the p=1, A=3
  instrument), so the image is the direct sum of three 2-dim pieces with
  distinct g-components: rank exactly 6. Any other value -> INVALID
  INSTRUMENT.
- P-KOSZ-R (reject): same plant with slice 3 identically zero: rank
  exactly 4. Must return exactly 4.

## 6. Controls (through the SAME code path)

- **A1 (accept known-valid):** the frozen exact T_C rank-3 witness and
  **A2** the frozen exact T_H rank-8 witness (both entrywise exact
  rational, from the gate-C campaign
  `2026-08-30T01:40:00ZZ_gateC_tau_sharpness`, construction semantics
  copied from
  `gate_c_sharpness_n2_n4.py` sha256
  `bcd6e32b08df5246694ecfa14467dea7a4cc925efd62e4755e0b9dae15ba934c`)
  must PASS substitution verification on their own tensors
  T_C = (I_2, [[0,1],[-1,0]]) and the 4x4x4 quaternion multiplication
  table T_H.
- **R1 (reject one-coordinate perturbation of the target):** T_cor = TF
  with entry [0][0][0] += 1. BOTH the T4 12-term witness and the T5
  13-term point must FAIL on T_cor (>= 1 mismatched entry, mismatch count
  recorded exactly; a PASS here = FALSE CERTIFICATE = quarantine +
  INVALID INSTRUMENT). Ordering note (fixed now): A1/A2/R2 and the plants
  run BEFORE the construction; R1 necessarily runs immediately after T4
  and T5, BEFORE the battery and adjudication — the quarantine semantics
  are identical either way.
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
SOS-Gram RAM C(247 + floor(d/2), floor(d/2))^2 x 8 bytes; compare against
the workstation RAM (`sysconf`); freeze the table. The claim frozen is
ONLY: every degree >= 4 naive Ncert/SOS certificate is out of scale by
the tabulated factors on this workstation; the degree-3
constant-multiplier family has no existence guarantee and is recorded as
impractical-without-new-structure. A sizing negative about a route
family, never a rank claim.

## 8. Resource cap, accounting, and stopping rules

- Total campaign CPU cap: 2 h (`time.process_time` printed at every stage
  boundary; hard self-check aborts with verdict INCONCLUSIVE-BUDGET if
  exceeded). Expected: minutes to tens of minutes; the largest single
  object is the 168 x 224 exact RREF of instruments 8/9.
- `PYTHONDONTWRITEBYTECODE=1` in the environment; `sys.dont_write_bytecode
  = True` before anything else; `resource.setrlimit(RLIMIT_CPU, soft ->
  hard)` raised at the top of the instrument; `nice -n 10` at launch
  (shell), no in-process renice; single process, no threads, no daemons,
  no blocking waits.
- Stopping rules: every stage runs once; a failure records its full
  history and moves the verdict per section 9; **no parameter is ever
  adjusted after seeing a result** (the only fallbacks are the two
  deterministic branches fixed in T1 and T5).
- Reproduction: run command
  `PYTHONDONTWRITEBYTECODE=1 nice -n 10
  /Users/jinleic/jinleic-workspace/cs/.venv/bin/python n3_instrument.py`
  from the run dir; the instrument reads the three frozen N2 files by
  RELATIVE path
  (`../20260901T115801Z_8d6683fb_39c752e7eadd/...`) and records their
  sha256 in-run against the section-2 pins; environment: python 3.14.3,
  stdlib only.

## 9. Adjudication vocabulary (fixed now)

- **FROZEN-CERTIFIED**: stages A (anchors), C (A1, A2, R2, plants,
  guards), T1-T4, R1, T5/T6 (with E and E1 handled exactly per their
  pre-registered conditionals), battery+plants, and sizing all complete
  with their expected outcomes. Claims allowed, exactly:
  (i) rank_C(TF) in [LB, 12] with both ends exact and machine-verified
  (equality claimed iff LB = 12);
  (ii) rank_R(TF) >= 13 > 12 >= rank_C(TF): a certified strict
  real-over-complex rank gap for this tensor (floor consumed from the
  frozen chain; upper end from T4) — and the certified complex
  13-feasibility (12-term witness + split) that the demotion theorem
  consumes;
  (iii) the demotion theorem of T5/T6 (raw-system ideal proper over Q;
  E conditional on its T6 point; E1 handled per its pre-registered
  structural case; N2 clause-(i) void; mod-p {1} = bad-reduction artifact
  only; upward route requires real-specific certificates);
  (iv) the C3 sizing negative.
  **The {13,14} frontier does NOT move** under any outcome: no rank_R
  claim of any kind is made in either direction; the verdict must state
  verbatim: "Absence of a 13-witness is NOT evidence for 14; absence of
  an impossibility argument is NOT evidence for 13." The published window
  18 <= R_R(T_O) <= 25 is untouched; no published work is refuted.
- **FROZEN-INCONCLUSIVE (INVALID INSTRUMENT)**: any control break
  (A1/A2 fail; R1/R2 pass; plant value off; guard failure; T1 nullity
  anomaly). Nothing reported as evidence.
- **FROZEN-INCONCLUSIVE (CONSTRUCTION ANOMALY)**: T1-T4 residual
  nonzero / S0 singular under both fixed selections / T5 failure with
  both fixed fallbacks exhausted.
- **FROZEN-INCONCLUSIVE (BUDGET/CRASHED)**: cap abort or script failure,
  with partial history frozen.
- status.json mapping: FROZEN-CERTIFIED -> FROZEN-CERTIFIED;
  INVALID INSTRUMENT / CONSTRUCTION ANOMALY / BUDGET / CRASHED ->
  FROZEN-INCONCLUSIVE with the label in the verdict.

## 10. Rule-7 prospective scope

Covered: TF = blockdiag(tau, tau) (with the frozen rank-equivalence to
(L_1, L_i, L_j)) at rank exactly 13's complex feasibility side; the nine
fixed flattening instruments with the corrected divisors and two plants;
the fixed construction T1-T6 with the asymmetric pullback; the three
frozen N2 system files; the C3 sizing table at degrees 3..6. NOT covered,
hence decided by nothing here: rank_R(TF) in either direction ({13,14}
stays OPEN with floor 13 and upper 14 intact); any real-specific
certificate question; any rank other than 13's complex analogue; T_O or
any other slice family; border rank; tau's complex rank (named candidate
for a FOLLOW-UP gate, not committed here); the frozen N1/N2 campaigns
(consumed read-only); R_R(T_O).

## 11. Lifecycle and freeze plan

Commit this file (path-scoped) -> `campaign.py init --gate
n3-technology-switch --prereg n3_prereg_v2.md` from cs/oct-rank ->
byte-identical copy of this file into the minted run dir + PROVENANCE.md
(source commit hash + this file's sha256, copied AFTER init, BEFORE any
in-run compute; provenance also binds the superseded run
20260904T021928Z_66f9d4ee_138fcf06a676 and both prereg hashes) -> probe
disclosure carried by reference to the closed run's `probe_disclosure/` ->
write and run `n3_instrument.py` (sections 2-8, single pass, from
scratch) -> freeze (`campaign.py freeze`) with instrument, results JSON,
witness JSON (exact Q(i) coordinates), verdict notes -> `campaign.py
close` with the section-9 mapping -> append the target README
current-state section. Per Main's steering: if this gate lands
FROZEN-CERTIFIED, the next implied frontier gate (candidates: exact
rank_C(tau) via the same technology — the tower-gap pattern with T_C
(3 real / 2 complex) and tau as rungs; or a witness-structured exact
real-13 fusion attempt) will be preregistered and run with fresh controls
immediately after, in this same session.
