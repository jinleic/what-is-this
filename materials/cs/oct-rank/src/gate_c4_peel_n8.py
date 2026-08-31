"""gate_c4_peel_n8.py — Gate C4 (pre_statement C3): n = 8 peeling-constant
probe. Decide whether the constant 18 in 18 <= R_R(T_O), as produced by
substitution-peeling + Thm-2 pencil, can be lifted by any parametric choice
INSIDE that machinery. Expected verdict: NEGATIVE, with rule-7 scope.

Everything below is EXACT rational / symbolic arithmetic (fmpq, sympy).
No floats anywhere.

Pre-registered probes (pre_statement.md, Gate C addendum, committed before
this run):

S1. Pencil floor > 12? Build an EXACT 12-term rank-1 decomposition of
    (I_8, J_8), J_8 = direct sum of four [[0,1],[-1,0]] blocks, by
    direct-summing four copies of the exact 3-term (I_2, J_2) witness of
    this campaign (gate_c_sharpness_n2_n4.py). Verify entrywise in fmpq.
    Universal lower half 12 <= pencilRank(1, C): Lean key_bound (gate A).
    Every C consumed by the chain satisfies C^2 - 2aC + bI = 0, a^2 < b,
    hence normalizes to a J with J^2 = -I [ordered-field fact] and every
    such J is R-similar to J_8 [real canonical form; stated, standard].
    Pencil rank is conjugation-invariant (u_s -> P u_s, v_s^T -> v_s^T P^-1
    transports a witness). So pencilRank == 12 EXACTLY on the whole
    consumed class <=> the 12 witness here + Lean 12 floor.

S2. Stopping-point independence: peeling j slices of the 8-slice L-family
    (0 <= j <= 6) leaves an L-family of (8 - j) independent vectors; the
    same induction engine at k = 6 - j gives
        r >= j + (12 + (6 - j)) = 18   for EVERY j.
    Verified here by exact integer arithmetic for all j and the identity
    j + (6 - j) == 6 (symbolic in j - trivially exact).

S4. Kernel-dimension (#{independent columns of D}) slack: on the S1
    witness, verify EXACTLY that every inequality in the Lean key_bound
    core is TIGHT at r = 12:
      (i)  U D == 0                      (exact matrix identity)
      (ii) rank D == 4 == n/2            (exact fmpq rank)
      (iii) rank U == n == 8, ker U dim == r - n == 4 (exact)
    => no D-side slack for a sharper count.

S3 (NOT swept, recorded): whether a 3-slice residual family
    (L_u, L_v, L_w) has rank >= 14. The machinery gives only 5 + 13 = 18.
    No computation here addresses S3; a >= 14 bound would improve 18.

Pass/fail (decided ahead): S1 floor==12 AND S2 all-18 AND S4 all-tight
=> NEGATIVE verdict for this probe, filed with the rule-7 sentence.
Any improvement candidate => freeze + hub Main + STOP (not expected).
"""
import json
import sys

from flint import fmpq, fmpq_mat

Q0 = fmpq(0)
Q1 = fmpq(1)
Qm1 = fmpq(-1)


def fmpq_mat_from(rows):
    nr, nc = len(rows), len(rows[0])
    return fmpq_mat([[fmpq(rows[i][j]) for j in range(nc)]
                     for i in range(nr)])


def rank_of(M):
    return M.rref()[1]



# ---------------------------------------------------------------- S1 -------
# Exact 3-term witness for (I_2, J_2), from the campaign-verified Gauss
# identity in gate_c_sharpness_n2_n4.py:
#   T[i,j,k] = sum_m A2[i][m] B2[j][m] C2[k][m]; slice_i = matrix over
#   (j, k). Hence with u_m = column m of B2 (over j), v_m = column m of
#   C2 (over k), f_I(m) = A2[0][m], f_J(m) = A2[1][m]:
#     I_2  = sum_m f_I(m) u_m v_m^T ,   J_2 = sum_m f_J(m) u_m v_m^T .
A2 = [[1, 1, 0], [1, 0, 1]]
B2 = [[1, 1, 0], [1, 0, 1]]
C2 = [[0, 1, -1], [1, -1, -1]]

I2sum = [[fmpq(0), fmpq(0)], [fmpq(0), fmpq(0)]]
J2sum = [[fmpq(0), fmpq(0)], [fmpq(0), fmpq(0)]]
for m in range(3):
    fI_m = fmpq(A2[0][m])
    fJ_m = fmpq(A2[1][m])
    for j in range(2):
        for k in range(2):
            I2sum[j][k] += fI_m * fmpq(B2[j][m]) * fmpq(C2[k][m])
            J2sum[j][k] += fJ_m * fmpq(B2[j][m]) * fmpq(C2[k][m])
s1_witness_2x2_exact = (
    I2sum == [[fmpq(1), fmpq(0)], [fmpq(0), fmpq(1)]]
    and J2sum == [[fmpq(0), fmpq(1)], [fmpq(-1), fmpq(0)]])

# ---- 8x8 witness: direct sum of four shifted copies of the 2x2 witness
# on coordinate blocks (2p, 2p+1). One I-slice and one J-slice.
terms8 = []
for p in range(4):
    offs = (2 * p, 2 * p + 1)
    for m in range(3):
        u8 = [fmpq(0)] * 8
        v8 = [fmpq(0)] * 8
        for t, off in enumerate(offs):
            u8[off] = fmpq(B2[t][m])     # u over j = B2[:,m]
            v8[off] = fmpq(C2[t][m])     # v over k = C2[:,m]
        terms8.append((fmpq(A2[0][m]), fmpq(A2[1][m]), u8, v8,
                       f"blk{p}_m{m}"))

I8 = [[fmpq(0)] * 8 for _ in range(8)]
J8 = [[fmpq(0)] * 8 for _ in range(8)]
for p in range(4):
    J8[2*p][2*p + 1] = fmpq(1)
    J8[2*p + 1][2*p] = fmpq(-1)
    I8[2*p][2*p] = fmpq(1)
    I8[2*p + 1][2*p + 1] = fmpq(1)

I8sum = [[fmpq(0)] * 8 for _ in range(8)]
J8sum = [[fmpq(0)] * 8 for _ in range(8)]
for wI, wJ, u8, v8, tag in terms8:
    for i in range(8):
        for j in range(8):
            prod = u8[i] * v8[j]
            if prod != 0:
                I8sum[i][j] += wI * prod
                J8sum[i][j] += wJ * prod
s1_12term_exact = (I8sum == I8 and J8sum == J8)
terms8_count = len(terms8)

# J8^2 == -I8 exact
J8m = fmpq_mat_from(J8)
I8m = fmpq_mat_from(I8)
J8sq = J8m * J8m
s1_J8sq = all((J8sq[i, j] == (Qm1 if i == j else Q0)
               for i in range(8) for j in range(8)))

# ---------------------------------------------------------------- S2 -------
s2_loop_values = []
for j in range(0, 7):
    # peel j slices: contributes j; residual L-family of (8-j) vectors gets
    # 12 + (number of further peels to reach a pencil) = 12 + (6 - j);
    # total = j + 12 + (6 - j) = 18 for every j.
    total = j + 12 + (6 - j)
    s2_loop_values.append({"j": j,
                           "residual_family_slices": 8 - j,
                           "further_peels_to_pencil": 6 - j,
                           "lower": total,
                           "equals_18": total == 18})
s2_all_18 = all(d["equals_18"] for d in s2_loop_values)
# pure arithmetic identity, symbolic in j:
import sympy as sp
jj = sp.symbols('j', integer=True, nonnegative=True)
s2_symbolic = sp.expand(jj + 12 + (6 - jj) - 18) == 0

# ---------------------------------------------------------------- S4 -------
# Extract the key_bound matrices from the 12-term witness and verify all
# inequalities tight. Slice framework: I = sum alpha_s u_s v_s^T,
# J = sum beta_s u_s v_s^T with alpha_s = wI, beta_s = wJ.
# U (8 x 12): columns u_s. A (12 x 8): rows alpha_s v_s^T.
# B (12 x 8): rows beta_s v_s^T. D = A + B J.
alphas = [fmpq(int(t[0])) for t in terms8]
betas = [fmpq(int(t[1])) for t in terms8]
us = [t[2] for t in terms8]
vs = [t[3] for t in terms8]

U = [[fmpq(0)] * terms8_count for _ in range(8)]
for s_i, u8 in enumerate(us):
    for i in range(8):
        U[i][s_i] = u8[i]
A = [[fmpq(0)] * 8 for _ in range(terms8_count)]
B = [[fmpq(0)] * 8 for _ in range(terms8_count)]
for s_i in range(terms8_count):
    for j in range(8):
        A[s_i][j] = alphas[s_i] * (vs[s_i])[j]
        B[s_i][j] = betas[s_i] * (vs[s_i])[j]


def fmpq_mat_rows(rows):
    return fmpq_mat(rows)


Um = fmpq_mat_rows(U)
Am = fmpq_mat_rows(A)
Bm = fmpq_mat_rows(B)
Dm = Am + Bm * J8m

UD = Um * Dm
s4_UD_zero = all((UD[i, j] == Q0 for i in range(8) for j in range(8)))
rankD = rank_of(Dm)
rankU = rank_of(Um)
UA = Um * Am
UB = Um * Bm
s4_UA_I = all((UA[i, j] == (Q1 if i == j else Q0)
               for i in range(8) for j in range(8)))
s4_UB_J = all((UB[i, j] == J8m[i, j] for i in range(8) for j in range(8)))

s4 = {
    "UD_is_zero_exact": bool(s4_UD_zero),
    "UA_is_I_exact": bool(s4_UA_I),
    "UB_is_J_exact": bool(s4_UB_J),
    "rank_D": int(rankD),
    "rank_D_equals_n_half": bool(rankD == 4),
    "rank_U": int(rankU),
    "rank_U_equals_n": bool(rankU == 8),
    "ker_U_dim": int(terms8_count - rankU),
    "ker_U_dim_equals_r_minus_n": bool(terms8_count - rankU == 4),
    "all_inequalities_tight": bool(rankD == 4 and rankU == 8
                                   and terms8_count - rankU == 4),
}

# ---------------------------------------------------------------- verdict --
negative = bool(s1_witness_2x2_exact and s1_12term_exact and s1_J8sq
                and s2_all_18 and s2_symbolic and s4["all_inequalities_tight"])
improvement_candidate = not negative

verdict = {
    "gate": "C4/C3.n8_peel_constant_probe",
    "constant_under_probe": 18,
    "S1_pencil_floor": {
        "claim": "pencilRank(1, C) == 12 exactly on the class the chain "
                 "consumes (C^2 - 2aC + bI, a^2 < b): lower 12 = Lean "
                 "key_bound [gate A]; upper 12 = exact witness here, "
                 "conjugation-invariance of pencil rank + real canonical "
                 "form J^2 = -I ~ J_8 [stated, standard]",
        "witness_2x2_base_exact": bool(s1_witness_2x2_exact),
        "witness_12term_for_I8_J8_entrywise_exact": bool(s1_12term_exact),
        "terms": terms8_count,
        "J8_squared_is_minus_I_exact": bool(s1_J8sq),
        "floor_is_12_exactly": bool(s1_witness_2x2_exact
                                    and s1_12term_exact and s1_J8sq),
    },
    "S2_stopping_points": {
        "claim": "peel j then finish by pencil: lower 18 for every "
                 "j in {0..6}; the substitution framework has no better "
                 "stopping point",
        "per_j": s2_loop_values,
        "symbolic_identity_j_plus_12_plus_6_minus_j_minus_18_is_zero":
            bool(s2_symbolic),
        "all_18": bool(s2_all_18),
    },
    "S4_kernel_dimension_tightness": s4,
    "S3_three_slice_residuals_NOT_SWEPT": {
        "question": "does rank(L_u, L_v, L_w) >= 14 hold for independent "
                    "u, v, w in R^8?",
        "machinery_gives": "5 + 13 = 18 (peel 5 + n + n/2 + 1 = 13)",
        "status": "NOT swept in this campaign; a >= 14 bound would improve "
                  "18 and is NOT excluded by this work",
    },
    "rule7_scope": "This negative covers ONLY the substitution/peeling "
                   "framework terminating in a Thm-2 n=8 pencil (the exact "
                   "machinery of Oct18.lean: 6 peels + pencil floor 12). "
                   "NOT swept: any non-peeling lower-bound route for T_O "
                   "(higher-order flattening, Koszul/Young-type bounds, "
                   "intersection-theoretic methods, symmetry-adapted "
                   "bounds), any rank bound on 3-slice residual "
                   "L-families (L_u, L_v, L_w) beyond the peel-implied 13 "
                   "(S3 above), and any improvement via non-similar "
                   "(non-J^2=-I) normalizations — none of these were "
                   "computed on.",
    "verdict": "NEGATIVE" if negative else "IMPROVEMENT-CANDIDATE",
    "all_exact_checks": {
        "s1_witness_2x2": bool(s1_witness_2x2_exact),
        "s1_witness_12": bool(s1_12term_exact),
        "s1_J8sq": bool(s1_J8sq),
        "s2_all_18": bool(s2_all_18),
        "s2_symbolic": bool(s2_symbolic),
        "s4_tight": bool(s4["all_inequalities_tight"]),
    },
}
print("GATE_C4_JSON " + json.dumps(verdict))
if improvement_candidate:
    # per pre-statement: freeze evidence, message Main, STOP — not expected.
    print("GATE_C4_VERDICT: IMPROVEMENT-CANDIDATE (escalation path)",
          file=sys.stderr)
    sys.exit(2)
print("GATE_C4_VERDICT: NEGATIVE (as predicted)")
sys.exit(0)
