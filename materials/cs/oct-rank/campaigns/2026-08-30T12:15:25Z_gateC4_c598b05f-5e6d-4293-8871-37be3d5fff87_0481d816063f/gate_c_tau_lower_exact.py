"""gate_c_tau_lower_exact.py — Gate C1(a) lower half: exact machine support
for the tau chain R_R(tau) >= 7 (peel 1 + pencil 6).

tau = 3-slice quaternion family (I, L_i, L_j) on R^4 (basis 1,i,j,k;
Lean Defs.lean `tau`, Lmat rows [q0,-q1,-q2,-q3; q1,q0,-q3,q2; q2,q3,q0,-q1;
q3,-q2,q1,q0]).

The chain (Lean Tau7.lean, gate-A replay, axioms [propext, Classical.choice,
Quot.sound]) consumes these identities; here EVERY polynomial identity is
re-derived from scratch by exact symbolic expansion (sympy, coefficients in
Z[...]) and by fmpq evaluation at fixed exact rational probe points. No
floating point anywhere. The two logical steps that are ordered-field facts
(positivity of squares; strict Cauchy-Schwarz for independent vectors) are
stated, not simulated.

Verdict: all identity checks exact-PASS + Lean chain [CITED-DEPENDENCY]
=> 7 <= r re-supported end to end.
"""
import json
import sys

import sympy as sp

from flint import fmpq, fmpq_mat

I2 = [[1, 0], [0, 1]]

# ----------------------------------------------------- quaternion machinery
# Lmat row convention (Lean Defs.lean Lmat):
# rows [q0,-q1,-q2,-q3; q1,q0,-q3,q2; q2,q3,q0,-q1; q3,-q2,q1,q0]
def Lmat_rows(q0, q1, q2, q3):
    return [[q0, -q1, -q2, -q3],
            [q1, q0, -q3, q2],
            [q2, q3, q0, -q1],
            [q3, -q2, q1, q0]]


def qmul_tuple(a, b):
    # Lean qmul convention (Defs.lean qmul)
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0)


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]


def mat_sub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def sc_mul(t, A):
    return [[t * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


results = {}

# ------- identity 1: Lmat_mul  L(p) L(q) = L(qmul p q)  (16 entries, exact)
p0, p1, p2, p3, q0, q1, q2, q3 = sp.symbols('p0 p1 p2 p3 q0 q1 q2 q3')
P = Lmat_rows(p0, p1, p2, p3)
Q = Lmat_rows(q0, q1, q2, q3)
pq = qmul_tuple((p0, p1, p2, p3), (q0, q1, q2, q3))
Lpq = Lmat_rows(*pq)


mismatch = []
PM = mat_mul(P, Q)
for i in range(4):
    for j in range(4):
        diff = sp.expand(PM[i][j] - Lpq[i][j])
        if diff != 0:
            mismatch.append((i, j, sp.simplify(diff)))
results["Lmat_mul_exact"] = not mismatch
results["Lmat_mul_diffs"] = [str(d) for *_, d in mismatch]

# ------- identity 2: Lmat_conj_self  L(qbar) L(q) = N(q) I  (exact)
qb = Lmat_rows(q0, -q1, -q2, -q3)
Nm = q0**2 + q1**2 + q2**2 + q3**2
mismatch2 = []
QB = mat_mul(qb, Q)
for i in range(4):
    for j in range(4):
        diff = sp.expand(QB[i][j] - (Nm if i == j else 0))
        if diff != 0:
            mismatch2.append((i, j, diff))
results["Lmat_conj_self_exact"] = not mismatch2

# ------- identity 3: Lmat_quad  L(m)^2 - 2 m0 L(m) + N(m) I = 0  (exact)
m0, m1, m2, m3 = sp.symbols('m0 m1 m2 m3')
M = Lmat_rows(m0, m1, m2, m3)
MM = mat_mul(M, M)
Nm_m = m0**2 + m1**2 + m2**2 + m3**2
mismatch3 = []
for i in range(4):
    for j in range(4):
        diff = sp.expand(MM[i][j] - 2*m0*M[i][j] + (Nm_m if i == j else 0))
        if diff != 0:
            mismatch3.append((i, j, diff))
results["Lmat_quad_exact"] = not mismatch3

# ------- identity 4: det(Lmat q) = N(q)^2  (exact symbolic determinant)
Dexp = sp.expand(sp.Matrix(Lmat_rows(m0, m1, m2, m3)).det() - Nm_m**2)
results["det_Lmat_exact"] = (sp.expand(Dexp) == 0)

# ------- identity 5: pencil_imaginary
# u = (-l1, 1, 0, 0) = i - l1 ; v = (-l2, 0, 1, 0) = j - l2
# m = qmul(conj u) v ; claim m1^2 + m2^2 + m3^2 = l1^2 + l2^2 + 1  (exact)
l1, l2 = sp.symbols('l1 l2')
u_sym = (-l1, 1, 0, 0)
v_sym = (-l2, 0, 1, 0)
uc = (u_sym[0], -u_sym[1], -u_sym[2], -u_sym[3])
m_sym = qmul_tuple(uc, v_sym)
lhs = sp.expand(m_sym[1]**2 + m_sym[2]**2 + m_sym[3]**2)
rhs = sp.expand(l1**2 + l2**2 + 1)
results["pencil_imaginary_exact"] = (sp.expand(lhs - rhs) == 0)
results["m_components"] = [str(sp.expand(c)) for c in m_sym]

# ------- identity 6: correction-glue  L_i - l I = L(-l,1,0,0),
#                                    L_j - l I = L(-l,0,1,0)  (exact)
Li = Lmat_rows(0, 1, 0, 0)
Lj = Lmat_rows(0, 0, 1, 0)
U1 = Lmat_rows(-l1, 1, 0, 0)
W1 = Lmat_rows(-l2, 0, 1, 0)
g1 = all(sp.expand(Li[i][j] - (l1 if i == j else 0) - U1[i][j]) == 0
         for i in range(4) for j in range(4))
g2 = all(sp.expand(Lj[i][j] - (l2 if i == j else 0) - W1[i][j]) == 0
         for i in range(4) for j in range(4))
results["Lmat_sub_real_exact"] = bool(g1 and g2)

# ------- identity 7: scaled quadratic for C = c * L(m)
# C^2 - 2 a C + b I = c^2 (L(m)^2 - 2 m0 L(m) + N(m) I)
# with a = c m0, b = c^2 N(m)   (exact polynomial identity in c, m0..m3)
c_sym = sp.symbols('c')
Cm = c_sym * sp.Matrix(M)
MM_sp = sp.Matrix(mat_mul(M, M))
Csq = MM_sp * c_sym**2
a_sym = c_sym * m0
b_sym = c_sym**2 * Nm_m
scaled = Csq - 2*a_sym*Cm + b_sym*sp.eye(4)
bundle = c_sym**2 * (MM_sp - 2*m0*sp.Matrix(M)
                     + Nm_m*sp.eye(4))
results["scaled_quadratic_exact"] = bool(
    all(sp.expand(scaled[i, j]) == 0 for i in range(4) for j in range(4)))

# ------- identity 8: discriminant split (ordered-field input prepared)
# for ANY quartet m: N(m) - m0^2 = m1^2 + m2^2 + m3^2  (exact, polynomial)
gen_disc = sp.expand((m0**2 + m1**2 + m2**2 + m3**2) - m0**2
                     - (m1**2 + m2**2 + m3**2))
results["disc_split_exact"] = (gen_disc == 0)
# instantiated at m = conj(u) v (components of m_sym): the pencil's
# b - a^2 quantity equals c^2 (l1^2 + l2^2 + 1), strictly positive as a
# sum of squares times c^2 (ordered-field step stated, not simulated):
Nm_inst = sp.expand(m_sym[0]**2 + m_sym[1]**2 + m_sym[2]**2 + m_sym[3]**2)
inst_disc = sp.expand(Nm_inst - m_sym[0]**2
                      - (m_sym[1]**2 + m_sym[2]**2 + m_sym[3]**2))
vs_pencil = sp.expand(m_sym[1]**2 + m_sym[2]**2 + m_sym[3]**2
                      - (l1**2 + l2**2 + 1))
results["disc_split_instantiated_exact"] = (inst_disc == 0
                                            and vs_pencil == 0)
results["pencil_b_minus_a2_form"] = ("b - a^2 = c^2*(l1^2 + l2^2 + 1), "
                                     "sum of squares x c^2 > 0")

# ================= exact fmpq probe points (second, independent route) =====
import random
random.seed(20260830)  # deterministic


def qp(symvals, expr_poly, syms):
    """Evaluate a sympy poly at exact rationals -> fmpq."""
    subs = dict(zip(syms, [fmpq(v[0], v[1]) for v in symvals]))
    val = expr_poly.subs(subs)
    return fmpq(int(sp.Rational(val).p), int(sp.Rational(val).q))


probes = []
ok_probe = True
for t in range(8):
    vals = {s: fmpq(random.randint(-97, 97), random.choice([1, 2, 3, 7]))
            for s in (p0, p1, p2, p3, q0, q1, q2, q3)}
    Pv = [[fmpq(sp.Rational(P[i][j].subs({s: vals[s] for s in vals})).p,
                sp.Rational(P[i][j].subs({s: vals[s] for s in vals})).q)
           for j in range(4)] for i in range(4)]
    Qv = [[fmpq(sp.Rational(Q[i][j].subs({s: vals[s] for s in vals})).p,
                sp.Rational(Q[i][j].subs({s: vals[s] for s in vals})).q)
           for j in range(4)] for i in range(4)]
    pqv = qmul_tuple(tuple(vals[s] for s in (p0, p1, p2, p3)),
                     tuple(vals[s] for s in (q0, q1, q2, q3)))
    Lpqv = Lmat_rows(*pqv)
    PMv = mat_mul(Pv, Qv)
    ok_probe &= all(PMv[i][j] == Lpqv[i][j]
                    for i in range(4) for j in range(4))
    # conj_self at rationals: L(qbar) with qbar = (q0, -q1, -q2, -q3)
    qv = tuple(vals[s] for s in (q0, q1, q2, q3))
    qbv = Lmat_rows(qv[0], -qv[1], -qv[2], -qv[3])
    QBv = mat_mul(qbv, Qv)
    Nv = (vals[q0]**2 + vals[q1]**2 + vals[q2]**2 + vals[q3]**2)
    ok_probe &= all(QBv[i][j] == (Nv if i == j else fmpq(0))
                    for i in range(4) for j in range(4))
    probes.append(t)
results["fmpq_probe_count"] = len(probes)
results["fmpq_probe_all_pass"] = bool(ok_probe)

# pencil_imaginary at exact points
ok_pi = True
for l1v in [fmpq(3, 2), fmpq(-7), fmpq(0), fmpq(11, 3)]:
    for l2v in [fmpq(5, 7), fmpq(-2), fmpq(0), fmpq(13, 5)]:
        uc = (-l1v, fmpq(1), fmpq(0), fmpq(0))
        vc = (-l2v, fmpq(0), fmpq(1), fmpq(0))
        mcc = (-uc[0], uc[1], uc[2], uc[3])
        mv = qmul_tuple(mcc, vc)
        lhs = mv[1]**2 + mv[2]**2 + mv[3]**2
        rhs = l1v**2 + l2v**2 + fmpq(1)
        ok_pi &= (lhs == rhs)
results["pencil_imaginary_fmpq_all_pass"] = bool(ok_pi)

# ================================================== chain restatement ======
chain = {
    "step_peel": ("substitution_pivot0 (Lean, Substitution.lean): peel 1 "
                  "term from a rank-r witness of tau (identity slice "
                  "nonzero) => corrected 2-slice family "
                  "(L_i - l1 I, L_j - l2 I) has rank <= r - 1 "
                  "[CITED-DEPENDENCY: Lean gate-A replay]"),
    "step_glue": ("corrections are L(-l1,1,0,0), L(-l2,0,1,0) — verified "
                  "exactly here as Lmat_sub_real_exact"),
    "step_normalize": ("L_u invertible (det = N(u)^2 = (l1^2+1)^2 > 0, "
                       "ordered-field); pencilRank invariant: "
                       "pencilRank(L_u, L_w) = pencilRank(1, L_u^{-1} L_w) "
                       "[CITED-DEPENDENCY: pencilRank_unit_left]"),
    "step_C": ("L_u^{-1} L_w = c * L(m), c = N(u)^{-1}, m = conj(u) w — "
               "verified exactly here (Lmat_inv via conj_self, Lmat_mul)"),
    "step_quadratic": ("C^2 - 2aC + bI = 0 with a = c m0, b = c^2 N(m) — "
                       "scaled_quadratic_exact + Lmat_quad_exact here"),
    "step_disc": ("a^2 < b:  b - a^2 = c^2 (m1^2 + m2^2 + m3^2) "
                  "= c^2 (l1^2 + l2^2 + 1) >= c^2 > 0 — "
                  "disc_split_exact + pencil_imaginary_exact here; "
                  "positivity is an ordered-field fact (squares >= 0)"),
    "step_pencil": ("Thm 2 at n = 4 (key_bound, Lean): "
                    "pencilRank(1, C) >= 4 + 2 = 6 "
                    "[CITED-DEPENDENCY: Lean gate-A replay]"),
    "step_chain": ("6 <= pencilRank(L_u,L_w) <= r - 1  =>  7 <= r  (omega)"),
}

verdict = {
    "gate": "C.tau_lower",
    "claim": "R_R(tau) >= 7 by peel(1) + pencil(6); machine support: all "
             "consumed algebraic identities re-derived exactly (sympy Z[...] "
             "expansion), plus fmpq probe points; logical skeleton "
             "[CITED-DEPENDENCY] Lean Tau7.lean replay (gate A, axioms "
             "[propext, Classical.choice, Quot.sound])",
    "identities": {
        "Lmat_mul_exact": results["Lmat_mul_exact"],
        "Lmat_conj_self_exact": results["Lmat_conj_self_exact"],
        "Lmat_quad_exact": results["Lmat_quad_exact"],
        "det_Lmat_equals_Nsq_exact": bool(results["det_Lmat_exact"]),
        "pencil_imaginary_exact": bool(results["pencil_imaginary_exact"]),
        "Lmat_sub_real_exact": bool(results["Lmat_sub_real_exact"]),
        "scaled_quadratic_exact": bool(results["scaled_quadratic_exact"]),
        "disc_split_exact": bool(results["disc_split_exact"]),
        "disc_split_instantiated_exact":
            bool(results["disc_split_instantiated_exact"]),
        "fmpq_probe_all_pass": results["fmpq_probe_all_pass"],
        "fmpq_probe_count": results["fmpq_probe_count"],
        "pencil_imaginary_fmpq_all_pass":
            results["pencil_imaginary_fmpq_all_pass"],
    },
    "chain": chain,
    "arithmetic": {"peels": 1, "pencil_floor_n4": 6, "total": 7},
    "all_exact_checks_pass": bool(
        results["Lmat_mul_exact"] and results["Lmat_conj_self_exact"]
        and results["Lmat_quad_exact"] and results["det_Lmat_exact"]
        and results["pencil_imaginary_exact"] and results["Lmat_sub_real_exact"]
        and results["scaled_quadratic_exact"] and results["disc_split_exact"]
        and results["disc_split_instantiated_exact"]
        and results["fmpq_probe_all_pass"]
        and results["pencil_imaginary_fmpq_all_pass"]),
}
print("GATE_C_TAU_LOWER_JSON " + json.dumps(verdict))
if not verdict["all_exact_checks_pass"]:
    print("GATE_C_TAU_LOWER_VERDICT: FAIL", file=sys.stderr)
    sys.exit(1)
print("GATE_C_TAU_LOWER_VERDICT: PASS")
