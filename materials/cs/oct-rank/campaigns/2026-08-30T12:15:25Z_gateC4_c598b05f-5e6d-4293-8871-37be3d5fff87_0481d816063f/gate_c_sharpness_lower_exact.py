"""gate_c_sharpness_lower_exact.py — Gate C2 lower halves, exact.

Paper sharpness wording (arXiv:2608.16649v1, section 2, after Thm 4's
R_R(T_C) >= 3, R_R(T_H) >= 8, R_R(T_O) >= 18 display):

    "The first two are the exact ranks [7, 11], so on the
     classical cases the bound is sharp."

(abstract: "sharp for C and H"). At n = 2: (5/2)n - 2 = 3; at n = 4: 8;
at n = 8: 18.

This file makes the n = 2 case FULLY SELF-CONTAINED machine-matter:

  LOWER bound R_R(T_C) >= 3 by own algebra (not the paper's pencil Thm 2).
  Hypothesis to refute: the 2-slice family (I, J), J = [[0,1],[-1,0]],
  has slice rank <= 2, i.e. there are u_s, v_s in R^2, scalars
  f_s0, f_s1 with
      I = sum_s f_s0 u_s v_s^T,   J = sum_s f_s1 u_s v_s^T   (s = 1, 2).
  Then for EVERY (x, y):
      xI + yJ = U diag(g1, g2) Vt,  g_s = x f_s0 + y f_s1,
  with U = [u1 u2], Vt rows v1^T, v2^T. Taking determinants (2x2 exact
  multiplicativity) and using det(xI + yJ) = x^2 + y^2:
      x^2 + y^2 = W (x f10 + y f11)(x f20 + y f21),  W = det(U) det(Vt).
  Matching coefficients gives the EXACT system
      E1: W f10 f20 = 1,  E2: W f11 f21 = 1,  E3: W (f10 f21 + f11 f20) = 0.
  Derivation (ring arithmetic over the coefficient identities):
      (E1)(E2): W^2 f10 f20 f11 f21 = 1;  with P = f10 f21, Q = f11 f20:
      f10 f20 f11 f21 = P Q (exact reassociation); E3 (W != 0 from E1):
      P + Q = 0.  Hence W^2 P (-P) = 1, i.e. (W P)^2 = -1.
  But t^2 + 1 has NO real root (exact Sturm count: 0). Contradiction.
  => R_R(T_C) >= 3.

  Machine-verified parts (all exact, sympy over Z[...] / Q[...]):
    [A] det(xI + yJ) - (x^2 + y^2) == 0 (polynomial identity).
    [B] det(U diag(g1,g2) Vt) - det(U) g1 g2 det(Vt) == 0 (generic; and at
        an exact rational probe point).
    [C] coefficient matching: each of E1, E2, E3 holds as a polynomial
        identity in the 12 symbols (u,v,a,b) GIVEN [A]+[B] — verified by
        expanding (x^2+y^2) - W g1 g2 and reading coefficients.
    [D] the derivation chain to (W P)^2 + 1 == 0 as exact reassociations.
    [E] Sturm sequence of t^2 + 1: exactly 0 real roots (exact algorithm
        over Q; this is the only analytic input, classical).
  Ordered-field inputs used and STATED: squares are >= 0; nonzero reals
  invertible; equality-multiplication in an integral domain. No floats.

  => R_R(T_C) >= 3 MACHINE-VERIFIED by own computation; with the exact
  3-term upper witness (gate_c_sharpness_n2_n4.py, this campaign):
  R_R(T_C) = 3, no external dependency.

For n = 4: R_R(T_H) = 8 — upper 8-term witness exact (existing script);
lower peel(2)+pencil(6) [CITED-DEPENDENCY: Lean replay, gate A]. Identity
support for the quaternion chain verified in gate_c_tau_lower_exact.py.

n = 8 (context): (5/2)*8 - 2 = 18 is exactly the constant probed by the
C4/C3 machinery — see gate_c4_peel_n8.py.
"""
import json
import sys

import sympy as sp

x, y = sp.symbols('x y')
I2 = sp.eye(2)
J2 = sp.Matrix([[0, 1], [-1, 0]])

# [A] det(x I + y J) = x^2 + y^2, exact
check_a = sp.expand((x*I2 + y*J2).det() - (x**2 + y**2)) == 0

# [B] determinant multiplicativity on the 2-term family, generic symbols
u1i, u1j, u2i, u2j, v1i, v1j, v2i, v2j, f10, f11, f20, f21 = sp.symbols(
    'u1i u1j u2i u2j v1i v1j v2i v2j f10 f11 f20 f21')
U = sp.Matrix([[u1i, u2i], [u1j, u2j]])           # columns u1, u2
Vt = sp.Matrix([[v1i, v1j], [v2i, v2j]])          # rows v1^T, v2^T
g1 = x*f10 + y*f11
g2 = x*f20 + y*f21
W = U.det() * Vt.det()
lhs_B = (U * sp.diag(g1, g2) * Vt).det()
check_b = sp.expand(lhs_B - W*g1*g2) == 0

# rational probe point for [B]
subs_pt = {u1i: 3, u1j: -1, u2i: 2, u2j: 5,
           v1i: 4, v1j: 7, v2i: -2, v2j: 1,
           f10: 1, f11: 2, f20: -3, f21: 1, x: sp.Rational(5, 3),
           y: sp.Rational(-2, 7)}
check_b_pt = sp.expand(sp.Rational(lhs_B.subs(subs_pt))
                       - sp.Rational((W*g1*g2).subs(subs_pt))) == 0
# [C] coefficient matching: (x^2 + y^2) - W g1 g2 == 0 forces E1, E2, E3
resid = sp.Poly(sp.expand((x**2 + y**2) - W*g1*g2), x, y)
E1 = sp.expand(W*f10*f20 - 1)          # x^2 coeff == 0  <=>  E1 == 0
E2 = sp.expand(W*f11*f21 - 1)          # y^2 coeff == 0  <=>  E2 == 0
E3 = sp.expand(W*(f10*f21 + f11*f20))  # xy  coeff == 0  <=>  E3 == 0
c_x2 = sp.expand(resid.coeff_monomial(x**2))
c_y2 = sp.expand(resid.coeff_monomial(y**2))
c_xy = sp.expand(resid.coeff_monomial(x*y))
check_c = (c_x2 + E1 == 0 and c_y2 + E2 == 0 and c_xy + E3 == 0)

# [D] derivation chain (pure ring identities; the FACTS are the
# coefficient equalities E1, E2, E3 == 0 verified in [C], with W != 0
# implied by E1):
P = f10*f21
Q = f11*f20
# s1: W^2 P Q == (W f10 f20)(W f11 f21)   [= 1 by E1, E2]
s1 = sp.expand(W**2 * P * Q - (W*f10*f20)*(W*f11*f21)) == 0
# s2: W^2 P (P + Q) == W^2 P^2 + W^2 P Q   [P + Q == 0 by E3, W != 0
#     => W^2 P^2 = - W^2 P Q]
s2 = sp.expand(W**2 * P * (P + Q) - (W**2 * P**2 + W**2 * P * Q)) == 0
# s3: (W P)^2 == W^2 P^2   [=> (WP)^2 = -W^2 P Q = -1]
s3 = sp.expand((W*P)**2 - W**2 * P**2) == 0
d_chain = bool(s1 and s2 and s3)

# [E] Sturm: t^2 + 1 has exactly 0 real roots (exact over Q)
t = sp.symbols('t')
chain_t = sp.sturm(t**2 + 1, t)


def sign_at_infinity(poly, at_neg):
    lc = sp.LC(poly, t)
    deg = sp.degree(poly, t)
    if at_neg:
        return sp.sign(lc) * (-1)**deg
    return sp.sign(lc)


def variations_at_infinity(chain, at_neg):
    signs = [sign_at_infinity(c, at_neg) for c in chain]
    return sum(1 for a, b in zip(signs, signs[1:])
               if a != 0 and b != 0 and a != b)


V_neg = variations_at_infinity(chain_t, True)
V_pos = variations_at_infinity(chain_t, False)
real_roots_t2p1 = V_neg - V_pos
check_e = (real_roots_t2p1 == 0)

# sanity: Sturm machinery is not degenerate — t^2 - 1 has exactly 2 roots
chain_s = sp.sturm(t**2 - 1, t)
check_e_sanity = (variations_at_infinity(chain_s, True)
                  - variations_at_infinity(chain_s, False) == 2)

n2 = {
    "claim": "R_R(T_C) >= 3 (self-contained); with exact 3-term upper "
             "witness: R_R(T_C) = 3",
    "A_det_xI_yJ_is_sum_sq_exact": bool(check_a),
    "B_det_multiplicativity_generic_exact": bool(check_b),
    "B_det_multiplicativity_rational_point": bool(check_b_pt),
    "C_coefficient_system_exact": bool(check_c),
    "D_derivation_chain_exact": bool(d_chain),
    "E_sturm_t2_plus_1_zero_real_roots": bool(check_e),
    "E_sturm_sanity_t2_minus_1_two_roots": bool(check_e_sanity),
    "sturm_variations": {"V_minus_inf": V_neg, "V_plus_inf": V_pos},
    "ordered_field_inputs_stated": [
        "squares are nonnegative; t^2 + 1 > 0 (Sturm certifies 0 roots)",
        "W != 0 follows from E1 (W f10 f20 = 1)",
        "multiplying equalities: integral-domain laws",
    ],
    "conclusion": "rank <= 2 hypothesis forces (W P)^2 = -1 over R; "
                  "impossible => R_R(T_C) >= 3. MACHINE-VERIFIED by own "
                  "computation (no pencil theory, no external lower dep).",
}


n4 = {
    "claim": "R_R(T_H) = 8; upper 8-term exact witness in "
             "gate_c_sharpness_n2_n4.py (this campaign), lower "
             "peel(2)+pencil(6) [CITED-DEPENDENCY: Lean replay gate A]; "
             "quaternion chain identities re-derived exactly in "
             "gate_c_tau_lower_exact.py (this batch)",
}

verdict = {
    "gate": "C.sharpness_lower",
    "paper_quote_verbatim": "The first two are the exact ranks [7, 11], so "
                            "on the classical cases the bound is sharp.",
    "paper_abstract_quote_verbatim": "sharp for $\\mathbb{C}$ and "
                                     "$\\mathbb{H}$",
    "bound_values": {"n2": "(5/2)*2 - 2 = 3", "n4": "(5/2)*4 - 2 = 8",
                     "n8": "(5/2)*8 - 2 = 18"},
    "n2_complex": n2,
    "n4_quaternion": n4,
    "all_exact_checks_pass": bool(
        check_a and check_b and check_b_pt and check_c and d_chain
        and check_e and check_e_sanity),
}
print("GATE_C_SHARPNESS_LOWER_JSON " + json.dumps(verdict))
if not verdict["all_exact_checks_pass"]:
    print("GATE_C_SHARPNESS_LOWER_VERDICT: FAIL", file=sys.stderr)
    sys.exit(1)
print("GATE_C_SHARPNESS_LOWER_VERDICT: PASS")
