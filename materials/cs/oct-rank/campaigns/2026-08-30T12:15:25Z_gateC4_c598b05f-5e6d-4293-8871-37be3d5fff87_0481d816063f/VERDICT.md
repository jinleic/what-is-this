# Gate C4 frozen verdict — peeling-constant probe NEGATIVE; C1(a)/C2 lower halves done

UTC freeze: 2026-08-30T12:15Z (campaign dir timestamp). Agent: OctRankGateC.

Pre-registration: pre_statement.md, section "Gate C addendum C4/C3 — n=8
peeling-constant probe [DECIDED AHEAD, 2026-08-30]" was appended BEFORE any
probe computation was run; probes, pass/fail rule, and escalation path were
fixed there. This artifact follows it.

## Subtask 3 — THE C4 PROBE: verdict NEGATIVE (as pre-registered)

Question: can the constant 18 in 18 <= R_R(T_O) be lifted by any parametric
choice INSIDE the substitution-peeling + Thm-2-pencil machinery of
arXiv:2608.16649 (Lean Oct18.lean: 6 peels + pencil floor 12)?

Probes (all exact fmpq/sympy, no floats — MACHINE-VERIFIED):

S1. Pencil floor. An exact 12-term entrywise-verified rank-1 decomposition
    of (I_8, J_8) is constructed by direct-summing four shifted copies of
    the campaign-verified 3-term (I_2, J_2) Gauss witness
    (gate_c_sharpness_n2_n4.py: A2=B2=[[1,1,0],[1,0,1]],
    C2=[[0,1,-1],[1,-1,-1]]; reverified in fmpq here). J_8^2 = -I_8 also
    verified entrywise. With the universal floor 12 = n + n/2 from Lean
    key_bound [CITED-DEPENDENCY, gate A replay] and the standard facts that
    every C with C^2 - 2aC + bI = 0, a^2 < b normalizes to a J with
    J^2 = -I, that every such J is R-similar to J_8 (real canonical form;
    irreducible quadratic t^2+1), and that pencil rank is invariant under
    the conjugation transporting (I, J) to (I, P J P^-1):
    pencilRank(1, C) = 12 EXACTLY on the whole class the chain consumes.
    => the pencil half of the 18 is saturated; no gain there.
    Note [INFERENCE-label]: the conjugacy/normalization steps are textbook
    real canonical-form facts, stated not re-proved here; the exact
    arithmetic (witness, J8^2 = -I) is MACHINE-VERIFIED.

S2. Stopping points. Peeling j slices (0 <= j <= 6) of the 8-slice
    L_{e_p}-family and finishing by the pencil bound gives
    j + (12 + (6 - j)) = 18 for EVERY j; verified numerically for all
    seven values and symbolically (jj + 12 + (6-jj) - 18 == 0).
    => the substitution framework has no better stopping point.

S4. Kernel-dimension (#{independent columns of D}) slack. On the 12-term
    witness, with Lean key_bound's matrices U (8x12), A, B, D = A + B J:
      U D = 0, U A = I_8, U B = J_8 (all entrywise exact),
      rank D = 4 = n/2 exactly, rank U = 8 = n exactly,
      dim ker U = 12 - 8 = 4 = r - n exactly.
    => every inequality in the core count is TIGHT at r = 12; there is no
    D-side slack for a sharper count of this shape.

Rule 7 (scope of the negative), verbatim from the script's rule7_scope:
  "This negative covers ONLY the substitution/peeling framework terminating
   in a Thm-2 n=8 pencil (the exact machinery of Oct18.lean: 6 peels +
   pencil floor 12). NOT swept: any non-peeling lower-bound route for T_O
   (higher-order flattening, Koszul/Young-type bounds, intersection-
   theoretic methods, symmetry-adapted bounds), any rank bound on 3-slice
   residual L-families (L_u, L_v, L_w) beyond the peel-implied 13 (S3
   above), and any improvement via non-similar (non-J^2=-I) normalizations
   — none of these were computed on."

S3 (open, NOT swept): does rank(L_u, L_v, L_w) >= 14 hold for independent
u, v, w in R^8? The machinery gives only 5 + 13 = 18 (peel 5 + n + n/2 + 1).
A >= 14 bound would improve 18 to 19 and is NOT excluded by this work.

Escalation clause: no improvement candidate appeared; nothing escalated to
Main for refutation risk. (An "IMPROVEMENT-CANDIDATE" exit during
development was a witness-construction bug — wrong factor orientation in
the 8x8 reduction — caught by the exact entrywise checks and fixed; it was
never a mathematical signal and was not escalated.)

## Subtask 1 — R_R(tau) = 7, lower half [gate C1(a)]

gate_c_tau_lower_exact.py, verdict PASS. Every polynomial identity the tau
chain consumes re-derived from scratch by exact expansion (sympy, Z[...]):

  Lmat_mul (L(p)L(q) = L(qmul p q)), Lmat_conj_self (L(qbar)L(q) = N(q) I),
  Lmat_quad (L(m)^2 - 2 m0 L(m) + N(m) I = 0),
  det(Lmat q) = N(q)^2,
  pencil_imaginary (|Im(conj(u) v)|^2 = l1^2 + l2^2 + 1 for
  u = (-l1,1,0,0), v = (-l2,0,1,0)),
  Lmat_sub_real (L_i - l1 I = L(-l1,1,0,0); L_j - l2 I = L(-l2,0,1,0)),
  scaled quadratic (C^2 - 2aC + bI = c^2(L(m)^2 - 2m0 L(m) + N(m) I)
  with a = c m0, b = c^2 N(m)),
  discriminant split (b - a^2 = c^2(l1^2 + l2^2 + 1); instantiated at
  m = conj(u) v),
plus a second independent route: 8 deterministic fmpq probe points
(seed 20260830) for Lmat_mul/Lmat_conj_self and an fmpq grid for
pencil_imaginary. All exact, all PASS, no floats anywhere.
Ordered-field inputs (squares >= 0; strict Cauchy-Schwarz for independent
vectors) are stated, not simulated. The logical skeleton is the gate-A
Lean replay (tau_rank_ge_seven, axioms [propext, Classical.choice,
Quot.sound]) [CITED-DEPENDENCY]. Together with the upper replay (re-run
PASS: K = 97771244.../18268770... < 1, margins 9.946e-06 > 0):
R_R(tau) = 7 re-verified in both directions.

## Subtask 2 — sharpness of (5/2)n - 2 at n = 2 and n = 4 [gate C2]

Paper's sharpness wording, verbatim (arXiv:2608.16649v1, section 2, after
the R_R(T_C) >= 3, R_R(T_H) >= 8, R_R(T_O) >= 18 display):
  "The first two are the exact ranks [7, 11], so on the
   classical cases the bound is sharp."
Abstract: "sharp for C and H". Values: (5/2)*2 - 2 = 3; (5/2)*4 - 2 = 8;
(5/2)*8 - 2 = 18.

n = 2: R_R(T_C) = 3 now FULLY SELF-CONTAINED machine-matter.
  Upper: exact 3-term witness (prior campaign, reused; entrywise PASS).
  Lower (new, own algebra — NOT the paper's Thm 2):
    rank <= 2 hypothesis forces x^2 + y^2 = W (x f10 + y f11)(x f20 + y f21)
    for W = det(U)det(V); coefficient matching gives W f10 f20 = 1,
    W f11 f21 = 1, W(f10 f21 + f11 f20) = 0; ring arithmetic (W != 0) gives
    (W f10 f21)^2 + 1 = 0 over R, impossible since t^2 + 1 has zero real
    roots (exact Sturm certificate; sanity aux: t^2 - 1 has two).
    Machine checks A-E all exact and PASS.
  => R_R(T_C) >= 3 MACHINE-VERIFIED by own computation; equality exact.

n = 4: R_R(T_H) = 8. Upper: exact 8-term witness (prior campaign, reused).
  Lower: peel(2) + pencil(6) [CITED-DEPENDENCY: Lean replay gate A]; the
  quaternion chain identities consumed by the same machinery are re-derived
  exactly in gate_c_tau_lower_exact.py. Equality conditional on the CITED
  lower bound, as in the prior campaign verdict.

## How to run (from cs/oct-rank/src/)

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c_tau_lower_exact.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c_sharpness_lower_exact.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
      /Users/jinleic/jinleic-workspace/cs/.venv/bin/python gate_c4_peel_n8.py

Environment: python 3.14.3, flint 0.9.0, sympy 1.14.0, single-thread,
deterministic (no randomness except the fixed-seed fmpq probes, seed
20260830). Raw outputs frozen in this directory (*.out with exit codes).

## Retraction note (Rule 5)

During development the S1 8x8 witness was first built with the wrong factor
orientation (A-rows/B-rows instead of B-columns/C-columns), which made the
entrywise checks FAIL and briefly produced an "IMPROVEMENT-CANDIDATE"
verdict. That construction was wrong; the corrected witness (documented
above) verifies exactly. The failure is retained here as a recorded
retraction of the intermediate false reading; it was a bug in my own
witness bookkeeping, not evidence about rank 18.
