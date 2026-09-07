# N4 pre-statement — gate `n4-tower-complex-rank`: exact complex ranks on the
# quaternion tower (T_C, tau) with fresh controls; tower-gap pattern as
# COMPUTATIONAL-EVIDENCE

Campaign: to be minted by `python3 scripts/campaign.py init --gate
n4-tower-complex-rank --prereg n4_prereg.md` from cs/oct-rank AFTER this
file is committed (path-scoped). Owner agent: `OctRankSwitch`. Date:
2026-09-04.

## 0. Relation to N3 (FROZEN-CERTIFIED,
run 20260904T034332Z_f5a61843_0333458ed434)

N3 certified rank_C(TF) = 12 EXACTLY (Koszul LB 12; exact Q(i) 12-term
witness) and the demotion theorem. N4 applies the SAME technology — exact
Q(i) witness construction + exact flattening battery — to the SMALLER
tower tensors, with fresh in-run controls. N3's artifacts are consumed
READ-ONLY (construction pattern and discipline; all arithmetic rerun).
Technology justification is the frozen N3 VERDICT (this target's route
map); nothing here touches TRF or F4. The published window
18 <= R_R(T_O) <= 25 is untouched; no rank_R claim is made anywhere.

## 1. Registered questions and claims

Objects (frozen conventions, same cd_ok recursion):
- tau = (I_4, L_i, L_j), 3x4x4 integers (the quaternion 3-slice tensor;
  frozen rank_R(tau) = 7 machine-verified).
- T_C = (I_2, [[0,1],[-1,0]]), 2x2x2 (frozen rank_R(T_C) = 3 exact).
- T_H = the full 4x4x4 quaternion multiplication table (frozen
  rank_R(T_H) = 8 exact).

Claims to certify (any outcome reported honestly):
- (i) rank_C(tau) in [LB_tau, 6]: UB from a NEW exact 6-term Q(i)
  witness (the M2(C) regrouping at n=4: tau_C slices are
  I_std (x) I_mult, X_i (x) I_mult, X_j (x) I_mult with
  {X_i, X_j} = [[0,1],[-1,0]], [[0,-i],[-i,0]]; the W-tensor
  (I2, X_i, X_j) has its three slices in the span of the three exact
  rank-1 elements I2 +- iX_i and I2 + iX_j; doubling over the multiplicity
  factor gives exactly 3*2 = 6 rank-1 terms; asymmetric pullback
  b_orig = S0^-1 b_conj, c_orig = S0^T c_conj with S0 built by the
  registered provable cyclic-vector construction). LB from the fixed
  battery below. "rank_C(tau) = 6 EXACTLY" iff LB_tau = 6.
- (ii) rank_C(T_C) in [LB, 2]: exact 2-term Q(i) witness (I2 and X_i both
  lie in span{I2 + iX_i, I2 - iX_i}, each of the two elements rank-1);
  LB = mode-1 flattening rank. Expected exact 2.
- (iii) rank_C(T_H) in [LB, 8]: exact 8-term Q(i) witness (the algebra
  basis {I2, X_i, X_j, X_iX_j} of the conjugated slice space is spanned
  by the four rank-1 2x2 matrices e_a e_b^T, a,b in {1,2}, each doubled
  over the multiplicity factor); LB from the battery. The interval may
  remain open at the lower end; that is an honest outcome.
- (iv) IF (i) yields LB_tau = 6 AND (ii) is exact 2: the tower table
  T_C 3 > 2, tau 7 > 6, TF >= 13 > 12 is recorded as
  **COMPUTATIONAL-EVIDENCE for a real-over-complex gap-1 pattern** on
  this tower — explicitly labeled conjecture-grade; it is NOT a rank_R
  claim about any object and in particular does not decide
  rank((L_1, L_i, L_j)) in {13, 14}.

## 2. Fixed instruments and controls (fresh, in-run)

Battery on tau (exact fmpq ranks, divisors C(dim(A)-1, p)):
1. mode-1 flattening (3 x 16); LB = rank.
2. mode-2 flattening (4 x 12); LB = rank.
3. mode-3 flattening (4 x 12); LB = rank.
4. Koszul p=1, (A,B,C) = (1,2,3); LB = floor(rank/2).
5. Koszul p=1, (A,B,C) = (2,3,1); LB = floor(rank/3).
6. Koszul p=1, (A,B,C) = (3,1,2); LB = floor(rank/3).
7. Koszul p=2, (A,B,C) = (1,2,3); LB = floor(rank/1).
8. Koszul p=2, (A,B,C) = (2,3,1); LB = floor(rank/3).
9. Koszul p=2, (A,B,C) = (3,1,2); LB = floor(rank/3).
Battery_T_C: mode-1 flattening rank (2 x 4) only (LB 2 = registered UB).
Battery_T_H: instruments 1-4 and the A=mode-2/mode-3 p=1,p=2 variants,
divisors as registered (dim A = 4: C(3,1) = 3, C(3,2) = 3).

Plants (same closed forms as N3, fresh in-run): P-KOSZ 6, P-KOSZ-R 4,
P-KOSZ2 3, P-KOSZ2-R 2. Guards: Kronecker (I2 slow, 4-block layout,
plant M[j][k] = (j*k) mod 5 - 2 on the relevant size) and dimension
plants (genuinely mismatched shapes refused).

Accept controls: A1 frozen exact T_C rank-3 real witness and A2 frozen
exact T_H rank-8 real witness (arrays verbatim from
`gate_c_sharpness_n2_n4.py` sha256
`bcd6e32b08df5246694ecfa14467dea7a4cc925efd62e4755e0b9dae15ba934c`) PASS
substitution on their own tensors. Reject controls: R2 corrupted T_C
([0][0][0] += 1) FAILS the A1 witness; R1 corrupted tau ([0][0][0] += 1)
FAILS both the 6-term tau witness and the padded raw-system point.
Object anchors: tau constants vs frozen `s1a_tau_r7.txt` (sha256
`16fe2bc45fe5813567722934e3acc11e7615af1efed29ff6778a907d79af2530`)
entrywise; T_C/T_H tables recomputed from the frozen gate-C semantics.

## 3. Arithmetic backend, cap, discipline

Identical to N3 section 3/8: stdlib-only exact Fraction pairs over Q(i);
no floats in claim-relevant paths; `PYTHONDONTWRITEBYTECODE=1` +
`sys.dont_write_bytecode = True`; RLIMIT_CPU soft->hard with 7200 s cap
enforced in-process at every stage boundary (persisted budget-abort
marker; verdict FROZEN-INCONCLUSIVE BUDGET on breach); `nice -n 10`;
single deterministic pass; the only fallbacks are registered fixed lists.
Reproduction command printed in-run; frozen N2/gate-C files read by
relative path with in-run sha256 checks.

## 4. Adjudication vocabulary (fixed now)

- **FROZEN-CERTIFIED**: all controls clean, witnesses verified by
  substitution (T_C 8/8, tau 48/48, T_H 64/64 — zero mismatches
  asserted and recorded), batteries computed, registered conditionals
  handled. Claims: exactly the section-1 intervals with the exactness
  rule, plus (iv) IF its condition holds, labeled
  COMPUTATIONAL-EVIDENCE.
- **FROZEN-INCONCLUSIVE (INVALID INSTRUMENT)**: any control/plant/guard
  break.
- **FROZEN-INCONCLUSIVE (CONSTRUCTION ANOMALY)**: witness substitution
  mismatch, singular S0, no cyclic candidate in the fixed list, or
  alpha-solve non-uniqueness.
- **FROZEN-INCONCLUSIVE (BUDGET/CRASHED)**: cap abort / script failure.
- Frontier sentence (mandatory in every outcome): {13, 14} stays OPEN;
  "Absence of a 13-witness is NOT evidence for 14; absence of an
  impossibility argument is NOT evidence for 13."; 18 <= R_R(T_O) <= 25
  untouched; no published work refuted.

## 5. Rule-7 scope

Covered: rank_C of tau, T_C, T_H at the registered upper bounds 6, 2, 8
via the fixed constructions and batteries; the registered plants and
controls. NOT covered: rank_R of anything (all frozen real statements
consumed, none re-derived); T_O; border rank; any other triple or
family; TF's real rank; the N3 campaign (consumed read-only).

## 6. Lifecycle

Commit (path-scoped) -> init --gate n4-tower-complex-rank --prereg
n4_prereg.md from cs/oct-rank -> byte-identical copy + PROVENANCE.md
(commit + sha256; N3 run referenced read-only) -> n4_instrument.py ->
one authoritative pass -> freeze -> close per section 4 -> README
current-state append.
