# N3 pre-statement v3 — gate `n3-technology-switch`: exact complex-relaxation
# probing (Q(i) CP-witness construction + Koszul/mode flattening battery +
# ordered-field demotion of the elimination route)

Campaign: to be minted by `python3 scripts/campaign.py init --gate
n3-technology-switch --prereg n3_prereg_v3.md` from cs/oct-rank AFTER this
file is committed (path-scoped). Owner agent: `OctRankSwitch`. Date:
2026-09-04.

**v3 supersedes v2 (run 20260904T023254Z_a4dd077a_7130b640666a, closed
FROZEN-INCONCLUSIVE INVALID-PREREG before any authoritative launch; prereg
sha256 e85be086e96086e0dd2ec8fc6c0c1ae7ed855ae4ac8d70cce3fde40da1b2c5d1,
commit 3dda194), which superseded v1 (run 20260904T021928Z_66f9d4ee_
138fcf06a676, sha256 29e08a181d5409777380c9f1c846f4cdc44d71efea83ef2ac34
acba365194846, commit b4239c9). All three closure notes are in the frozen
run dirs.** v3 changes vs v2: (a) section 7 / section 9(iv) sizing claim
CORRECTED — the v2 blanket "every degree>=4 Ncert/SOS certificate is out
of scale by RAM" was false for the SOS branch at d=4,5; (b) the p=2
Koszul plant pair P-KOSZ2 / P-KOSZ2-R is registered in the controls;
(c) T5 registers the explicit true-TF substitution gate (all 192 entries
zero-mismatch asserted BEFORE any raw-system/properness record); (d) T6
registers direct gauge recomputation for the E1 witness terms (no index
reuse from the split system); (e) T1 registers the repaired PROVABLE
cyclic-vector S0 construction (the v1/v2 nullspace-basis rule produced
only singular candidates); (f) the eleven pre-run instrument defects of
the v2 run are enumerated in its PROVENANCE and bound here by reference.

## 0. Excluded technology validation (disclosure, per Main's steering)

One /tmp probe file ran twice pre-commit of the superseded v1 prereg:
`/tmp/n3probe/probe.py` (sha256
`f2cc548dd2e5baa567e12dda36ff293e4a4f0a7b78fc8c6cb38ca84f98b4ce60`;
preserved in the v1 run's probe_disclosure/). EXCLUDED technology
validation only: it located defects in the probe's own code and confirmed
frozen table conventions. **No scientific knob — technology, witness,
bound, instrument list, divisor, threshold, or stopping rule — was chosen,
tuned, or validated by any probe or by any aborted v2 start.** All v2
instrument defects were found by Main's static review or by crash-time
diagnostics BEFORE any adjudication, and are enumerated/repaired in the v2
run's PROVENANCE. Every decisive substitution/bound runs from scratch
inside THIS campaign.

## 1. Technology choice and justification (from frozen cost evidence)

**Technology: exact complex-relaxation probing.** Three exact, purely
linear-algebraic instruments over Q and Q(i), deterministic, no search:

- **C1 (witness construction):** an explicit exact Q(i)-rational CP
  decomposition of TF with the term count the construction yields
  (M2(C)-module regrouping: exactly 12 terms — a theorem of the
  construction, not a tuned parameter), verified by direct substitution
  into the exact integer tensor entries, then padded by an exact split to
  an exact 13-term point of the raw rank-13 CP system. **Complex
  13-feasibility of TF is an OUTPUT of this campaign** (the frozen real
  14-term bound gives a complex 14-term decomposition and implies nothing
  at 13).
- **C2 (lower-bound battery):** exact fmpq ranks of a FIXED list of nine
  mode-flattening / Koszul-flattening linear maps of TF, lower bounds via
  the Koszul lemma with divisors C(dim(A)-1, p) (section 5).
- **C3 (ordered-field demotion):** exact Q(i)-points of the raw N2
  rank-13 CP system and (conditionally, T6) of the gauge system E,
  certifying those ideals PROPER over Q; detection of the pre-registered
  possible self-contradiction in serialized E1; plus a sizing table whose
  claims are EXACTLY the corrected ones of section 7.

**Why this technology and not the frozen ones:**

1. N2 (FROZEN-NEGATIVE, `20260901T115801Z_8d6683fb_39c752e7eadd`): F4
   died inside degree 9 on a 66-variable calibration system (8,080,476 x
   74,499,764 Macaulay matrix, 22.65 GB, degree 8 alone 3802 s wall /
   3593 s CPU), 4x smaller than the 273-variable main systems. And it is
   structurally mispointed: {1}-over-Q requires an empty complex variety;
   this campaign certifies the raw system's ideal is PROPER by exhibiting
   an exact Q(i)-point, so N2's clause (i) was impossible ab initio and
   only mod-p artifacts were ever obtainable (which N2's own prereg
   section 3 already admits cannot be airtight). The elimination
   technology has **no path to a rank-14 verdict at any budget**.
   No msolve/F4/Groebner code is invoked here.
2. N1 (FROZEN-INCONCLUSIVE, `20260901T103338Z_06de66a3_d0b7bc2dc6a5`):
   40-seed TRF + exact-Newton: 0/40 admitted. This campaign never invokes
   TRF, seeds, or numerical admission. **No TRF and no F4 code path is
   rerun anywhere.**
3. The complex side is the unattacked flank: the only frozen complex-side
   fact is the trivial rank_C(TF) <= rank_R(TF) <= 14. Everything at 13
   and below over C is settled here, if at all, by exact construction
   (C1) and exact bounds (C2) — self-certifying by substitution and rank
   computation (rule 14).

This is assignment candidate family (b): an exact certificate route that
avoids Groebner bases entirely — feasibility witnesses over C + flattening
lower bounds + an ordered-field demotion proving which future 14-
certificates must be real-specific.

## 2. Object, conventions, and frozen anchors

As in v2, carried verbatim: TF = blockdiag(tau, tau) 3x8x8 integer tensor
with tau = (I_4, L_i, L_j) on the Cayley-Dickson basis (1, i, j, k),
(a,b)(c,d) = (ac - conj(d)b, da + b conj(c)), T[p][b][c] = (e_p * e_b)_c;
verbatim semantic copy of the frozen `cd_ok`/`cj` recursion
(`n2_systems.py` sha256
`52f6e65111d1fde64a96cae45c06704cfe0aad9ac30fc4f49576df6ca1f43c69`,
including the len-1 branch). Frozen rank-equivalence rank(TF) =
rank((L_1, L_i, L_j)) (Route-F); frozen floor rank_R(TF) >= 13 (1 peel +
pencil 12; Route A cap 12 < 13 intact); frozen upper rank_R(TF) <= 14;
frozen rank_R(tau) = 7; window 18 <= R_R(T_O) <= 25 untouched.
Byte-level anchors asserted in-run: (a) tau/TF nonzeros 12/24; (b) the 48
constants of `s1a_tau_r7.txt` (sha256 `16fe2bc4...af2530`) reproduce
-tau; (c) the first 192 constants of `s2e_tf_r13_ext.txt` (sha256
`117fe71e...31d33f6`) reproduce -TF; (d) concat-vs-Kronecker guard:
blockdiag(M,M) == I2 (x) M under the 4-block flat layout (identity SLOW),
fixed plant M[j][k] = (j*k) mod 5 - 2; (e) dimension plants refused.

## 3. Arithmetic backend

Python 3.14.3 venv interpreter; **stdlib only** (`fractions`, `re`,
`itertools`, `json`, `hashlib`, `math`, `resource`, `time`, `sys`, `os`).
Exact rational pairs (a, b) = a + b*i; exact RREF/Gauss over Q and Q(i).
No floats, no numpy, no BLAS in any claim-relevant path.

## 4. The fixed construction (theorems T1-T6, no free parameters)

Conjugated frame: E_conj = C^2_doub (x) C^2_std (x) C^2_mult, index
(d, s, m) -> 4d + 2s + m.

- **T1 (regrouping; PROVABLE cyclic-vector S0 construction — repaired at
  v3, disclosed).** X_i = [[0,1],[-1,0]], X_j = [[0,-i],[-i,0]] in
  M_2(Q(i)); both square to -I2 and anticommute; K_a = X_a (x) I_2. Solve
  NOTHING heuristically: construct S0 directly. Let
  y = sum_s e_s (x) e_s in the K-frame (cyclic for the M_2(C) (x) I_2
  rep by the vectorization bijection A |-> (A (x) I_2) y, so
  tgt = (y, K_i y, K_j y, K_i K_j y) is a basis of C^4). For the first x
  in the FIXED candidate list (e1, e2, e3, e4, e1+e2, e1+e3, e1+e4,
  e1+e2+e3, e1+e2+e3+e4) with B = (x, L_i x, L_j x, L_i L_j x) a basis of
  C^4 (det test), set S0 := tgt * B^{-1}. Then S0 L_a = K_a S0 holds on B
  (S0 x = y, S0 L_a x = K_a y) and hence on the generated algebra; S0 is
  invertible (both B and tgt are bases). In-run asserts: the nullspace of
  the 32x16 intertwiner system has dimension exactly 4 (cross-check; any
  other value -> CONSTRUCTION ANOMALY); det(S0) != 0; BOTH identities
  S0 L_i S0^-1 = K_i and S0 L_j S0^-1 = K_j entrywise. If no candidate x
  yields a basis: CONSTRUCTION ANOMALY, no retries (mathematically
  impossible for this rep — recorded should it ever occur).
- **T2 (rank-1 trinomials).** p1 = I2 + i*X_i, p2 = I2 - i*X_i,
  p3 = I2 + i*X_j (elements xI + yX_i + zX_j at (x,y,z) = (1,i,0),
  (1,-i,0), (1,0,i)); det(xI + yX_i + zX_j) = x^2 + y^2 + z^2 (X_i^2 =
  X_j^2 = -I2, X_iX_j + X_jX_i = 0), so all three p_r are rank <= 1;
  in-run factors p_r = u_r v_r^T exactly (first nonzero entry rule) and
  asserts the product entrywise.
- **T3 (alpha solve).** The 12x3 exact system [vec p1 vec p2 vec p3]
  alpha_col = vec(slice) for slices I2, X_i, X_j; in-run asserts unique
  solution (pivot columns {0,1,2}) and zero residual — certifying the
  three W-slices lie in span{p1, p2, p3}.
- **T4 (12-term witness; ASYMMETRIC pullback).** I_D (x) (u_r v_r^T)
  (x) I_M = sum_{d,m} (d-vec (x) u_r (x) m-vec)(d-vec (x) v_r (x)
  m-vec)^T gives exactly 3*2*2 = 12 rank-1 terms; with the alpha
  coefficients their sum is the conjugated-frame tensor with slices
  I2_d (x) [I4, K_i, K_j] (x) I2_m — verified on all 192 entries. Pull
  back: L = S0^-1 K S0, so P Q^T |-> (S0^-1 P)(S0^T Q)^T:
  **b_orig = S0^-1 b_conj, c_orig = S0^T c_conj** (alpha unchanged).
  In-run verifies all 192 integer TF entries in the original frame. Any
  mismatch -> INCONCLUSIVE (witness failure), residuals frozen.
- **T5 (13-term split; EXPLICIT true-TF gate).** Split term r=0 (fixed):
  (alpha_r, P, Q) -> (alpha_r, P, Q + W) + (-alpha_r, P, W), W = e_3 of
  the conjugated 8-space (nonzero, fixed); deterministic fallback r = 1
  then r = 2 if gauge conditions fail; no other fallback. THE INSTRUMENT
  THEN (i) verifies substitution of all 13 terms against the TRUE TF:
  all 3*64 = 192 entries must match with ZERO mismatches, the count is
  recorded, and the assertion fires BEFORE any raw-system or properness
  record; (ii) asserts all 13 term triples are nonzero. This point is an
  exact Q(i) solution of the raw 192-equation/247-unknown rank-13 CP
  system of TF.
- **T6 (gauge coordinates + E/E1 handling).** For the T5 13-term point,
  each term's (a[0,s], c[0,s]) != (0,0) is asserted exactly and the 2x2
  gauge system solved exactly per term (determinant a0^2 + c0^2 over
  Q(i); isotropic pairs recorded, not silently skipped: if any term is
  isotropic the E-point is not certified and the fact is recorded — the
  demotion theorem stands on the raw system alone). E-point: ALL 218
  equations of the frozen `s2e_tf_r13_ext.txt` evaluate to 0.
  E1-point: term 0 = (0,0,0) (pins a0_0 = c0_0 = 0), terms 1..12 = the
  T4 12-term witness, with gauge solutions recomputed DIRECTLY for these
  12 terms (no index reuse from the split system). **Pre-registered
  structural case:** if the serialized E1 pairs the pin with a gauge pair
  on the same term (the v1 generator gauged all 13 terms), that gauge
  polynomial evaluates to the constant -1 on the pin locus, so
  E1-as-serialized has the UNIT ideal by an explicit 3-polynomial
  combination and carries no rank content; the instrument detects and
  records this (failing indices + the combination) as a defect of the
  FROZEN N2 serialization (rule 5), not a campaign anomaly.
- **Consequence (demotion theorem, proved in VERDICT.md from T5):** a
  Q(i)-point of a Q-polynomial system implies its ideal is proper over Q
  (unit ideal => empty complex variety, Nullstellensatz). Hence the raw
  rank-13 CP system of TF has a proper elimination ideal over Q (E
  likewise if its T6 point certifies): **N2's clause-(i) {1}-over-Q
  verdict was impossible ab initio**; any future mod-p {1} outcome is a
  bad-reduction artifact with no char-0 content; any future rank >= 14
  certificate MUST be real/order-specific. Elimination over Q is closed
  structurally, not just by budget.

## 5. The fixed lower-bound battery (nine instruments + four plants)

Koszul lemma: for T = sum_r a_r (x) b_r (x) c_r, the flattening
K_T(alpha (x) beta) := sum_a (alpha ^ e_a) (x) T_a(beta),
T_a(beta) = (T[a][j][k])_{k} contracted against beta, satisfies
K_T = sum_r (kappa_{a_r} (x) ell_r), kappa_a: alpha |-> alpha ^ a of rank
C(dim(A)-1, p) on Lambda^p A, ell_r of rank 1; hence
R >= floor(rank(K_T) / C(dim(A)-1, p)). Divisors: p=1 dim A=3: 2;
p=1 dim A=8: 7; p=2 dim A=3: 1; **p=2 dim A=8: C(7,2)=21**.

Instruments (all exact fmpq ranks on the integer TF):
1. mode-1 flattening rank (3 x 64); LB = rank.
2. mode-2 flattening rank (8 x 24); LB = rank.
3. mode-3 flattening rank (8 x 24); LB = rank.
4. Koszul p=1, (A,B,C) = (1,2,3); LB = floor(rank/2).
5. Koszul p=1, (A,B,C) = (2,3,1); LB = floor(rank/7).
6. Koszul p=1, (A,B,C) = (3,1,2); LB = floor(rank/7).
7. Koszul p=2, (A,B,C) = (1,2,3); LB = floor(rank/1).
8. Koszul p=2, (A,B,C) = (2,3,1); LB = floor(rank/21).
9. Koszul p=2, (A,B,C) = (3,1,2); LB = floor(rank/21).

Plants (before the battery; hand-closed forms):
- P-KOSZ (p=1, A=3 accept): diagonal 3x4x4, rank exactly 6.
- P-KOSZ-R (p=1, A=3 reject): slice 3 zeroed, rank exactly 4.
- P-KOSZ2 (p=2, A=3 accept): same diagonal plant, rank exactly 3.
- P-KOSZ2-R (p=2, A=3 reject): slice 3 zeroed, rank exactly 2.
Any deviation -> INVALID INSTRUMENT.

LB := max over the nine. Reported result: rank_C(TF) in [LB, 12], both
ends exact; "rank_C(TF) = 12 EXACTLY" iff LB = 12. If LB < 12 the
interval is the honest result and the verdict says so.

## 6. Controls (through the SAME code path)

- **A1/A2 (accept):** the frozen exact T_C rank-3 witness and T_H rank-8
  witness (arrays verbatim from `gate_c_sharpness_n2_n4.py` sha256
  `bcd6e32b08df5246694ecfa14467dea7a4cc925efd62e4755e0b9dae15ba934c`)
  PASS substitution on their own tensors.
- **R1 (reject, target):** T_cor = TF with entry [0][0][0] += 1 must FAIL
  under BOTH the T4 12-term witness and the T5 13-term point (>= 1
  mismatch, count recorded; a PASS = FALSE CERTIFICATE = quarantine +
  INVALID INSTRUMENT). Ordering: A1/A2/R2 and plants run before the
  construction; R1 runs immediately after T4 and again after T5, before
  the battery and adjudication.
- **R2 (reject, control):** T_C with [0][0][0] += 1 must FAIL the A1
  witness.
- **G-DIM / G-KRON:** section 2(d),(e).

## 7. C3 sizing table (route negative, EXACT corrected claims)

For d in {3, 4, 5, 6} report exactly: M(d) = C(247+d, d),
N(d-3) = C(247+(d-3), d-3), dense Ncert coefficient RAM
M(d)*192*N(d-3)*8 bytes, SOS Gram dimension C(247+floor(d/2), floor(d/2))
and dense Gram RAM (dimension^2)*8 bytes; compare with workstation RAM
(sysconf). **Registered claims, and only these:** (i) the Ncert dense
coefficient route at d >= 4 exceeds this workstation's RAM by the
tabulated factor (~61.5 TB at d=4); (ii) the SOS route at d >= 6 exceeds
it (~52.96 TB at d=6); (iii) **SOS at d in {4,5} (~7.63e9 bytes dense)
and Ncert at d = 3 (~3.94e9 bytes dense) are NOT excluded by this
sizing** — they are recorded as unexcluded-by-RAM, with no runtime or
conditioning certificate claimed. This is a sizing statement about route
families, never a rank claim.

## 8. Resource cap, accounting, and stopping rules

- Total campaign CPU cap: 2 h (`time.process_time` at every stage
  boundary; abort -> INCONCLUSIVE-BUDGET). Expected minutes; largest
  object the 168 x 224 exact RREF (instruments 8/9).
- `PYTHONDONTWRITEBYTECODE=1` at launch; `sys.dont_write_bytecode = True`
  first line; `resource.setrlimit(RLIMIT_CPU, soft -> hard)` at top;
  `nice -n 10` at launch; single process, no threads, no daemons; no
  `py_compile` or any bytecode-writing tool is ever invoked inside the
  run dir (the v2 run's disclosed pyc debris lives outside it, in
  `scratch/n3_excluded_pyc/`).
- Stopping rules: every stage runs once; a failure records its history
  and moves the verdict per section 9; no parameter is adjusted after
  seeing a result (the only fallbacks are the fixed lists in T1 and T5).
- Reproduction: `PYTHONDONTWRITEBYTECODE=1 nice -n 10
  /Users/jinleic/jinleic-workspace/cs/.venv/bin/python n3_instrument.py`
  from the run dir; frozen N2 files read by relative path
  `../20260901T115801Z_8d6683fb_39c752e7eadd/...` with in-run sha256
  checks against section 2 pins; python 3.14.3, stdlib only.

## 9. Adjudication vocabulary (fixed now)

- **FROZEN-CERTIFIED**: stages A, C (A1, A2, R2, all four plants,
  guards), T1-T4, R1, T5 (including its explicit true-TF zero-mismatch
  gate), T6 (with the registered E/E1 conditionals), battery+plants, and
  sizing all complete with expected outcomes. Claims, exactly:
  (i) rank_C(TF) in [LB, 12], both ends exact and machine-verified
  (equality iff LB = 12);
  (ii) rank_R(TF) >= 13 > 12 >= rank_C(TF): a certified strict
  real-over-complex rank gap for this tensor (floor consumed from the
  frozen chain; upper end from T4), plus the certified complex
  13-feasibility (12-term witness + split);
  (iii) the demotion theorem of T5/T6 (raw ideal proper over Q; E
  conditional; E1 per its structural case; N2 clause-(i) void; mod-p {1}
  = bad-reduction artifact; upward route requires real-specific
  certificates);
  (iv) the section-7 sizing claims (i)-(iii) exactly as corrected there
  — no blanket out-of-scale claim.
  **The {13,14} frontier does NOT move** under any outcome: no rank_R
  claim in either direction; the verdict states verbatim: "Absence of a
  13-witness is NOT evidence for 14; absence of an impossibility argument
  is NOT evidence for 13." The published window 18 <= R_R(T_O) <= 25 is
  untouched; no published work is refuted.
- **FROZEN-INCONCLUSIVE (INVALID INSTRUMENT)**: control break (A1/A2
  fail; R1/R2 pass; any plant value off; guard failure; T1 nullity
  anomaly). Nothing reported as evidence.
- **FROZEN-INCONCLUSIVE (CONSTRUCTION ANOMALY)**: T1-T6 residual nonzero,
  S0 singular, no cyclic candidate, or T5 true-TF mismatch.
- **FROZEN-INCONCLUSIVE (BUDGET/CRASHED)**: cap abort or script failure
  with partial history frozen.
- status.json mapping: FROZEN-CERTIFIED -> FROZEN-CERTIFIED; all
  INCONCLUSIVE labels -> FROZEN-INCONCLUSIVE with the label in the
  verdict.

## 10. Rule-7 prospective scope

Covered: TF at rank exactly 13's complex feasibility side; the nine fixed
flattening instruments with correct divisors; the four registered plants;
the fixed construction T1-T6 with the asymmetric pullback and explicit
true-TF gates; the three frozen N2 system files; the section-7 sizing
table. NOT covered: rank_R(TF) in either direction ({13,14} stays OPEN,
floor 13 and upper 14 intact); any real-specific certificate question;
any other rank; T_O or other slice families; border rank; tau's complex
rank (named FOLLOW-UP candidate, not committed); the frozen N1/N2
campaigns (read-only); R_R(T_O).

## 11. Lifecycle and freeze plan

Commit this file (path-scoped) -> init --gate n3-technology-switch --prereg
n3_prereg_v3.md from cs/oct-rank -> byte-identical copy into the run dir +
PROVENANCE.md (source commit + sha256; BOTH superseded prereg hashes
bound) -> write n3_instrument.py implementing sections 2-8 -> ONE
authoritative single-pass execution -> freeze -> close per section 9 ->
append the target README current-state section. Per Main's steering: on
FROZEN-CERTIFIED, the next implied frontier gate (candidates: exact
rank_C(tau) via the same technology — the tower-gap pattern T_C 3/2, tau
7/6?; or a witness-structured exact real-13 fusion attempt) is
preregistered and run with fresh controls immediately after, same session.
