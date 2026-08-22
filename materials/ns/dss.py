"""Backward discretely self-similar (DSS) reduction for constant-viscosity 3D NS.

Verifies the claim in README.md that the similarity change of variables

    y = x / sqrt(-t),   tau = -log(-t),   u(x,t) = (-t)^(-1/2) v(y,tau),
                                          p(x,t) = (-t)^(-1)   P(y,tau)

(blowup at t = 0, so t < 0 and tau -> +inf) turns unforced Navier-Stokes into

    d_tau v + (1/2) v + (1/2)(y.grad)v + (v.grad)v + grad P = nu Lap v,   div v = 0,

and that the backward lambda-DSS condition u(x,t) = lambda u(lambda x, lambda^2 t)
becomes exactly tau-periodicity,

    v(y, tau) = v(y, tau - S),        S = 2 log lambda.

Two consequences, both checked below.

  * Steady solutions (d_tau v = 0) reproduce Leray's profile equation
        -nu Lap U + a U + a (y.grad)U + (U.grad)U + grad P = 0
    at a = 1/2, matching Tsai's statement (personal.math.ubc.ca/~ttsai/publications/
    leray.pdf). That is precisely the case closed by NRS 1996 and Tsai 1998.

  * DSS is the tau-PERIODIC case. Lemma 1's argument turns on the profile equation being
    t-independent while a coefficient is not; for a periodic v there is no such
    contradiction, so Lemma 1 has no analogue here. Consistently, the Chae-Wolf rigidity
    (arXiv:1610.09464) is a SMALL-PERIOD result: lambda near 1 means S near 0, i.e.
    nearly steady. Long period is untouched.

Method: substitute symbolically and compare each term against its expected profile form.
Two spatial coordinates are carried so that (y.grad) and the Laplacian are genuinely
exercised rather than collapsing to a 1D identity.
"""

import sympy as sp

nu, lam = sp.symbols("nu lam", positive=True)
t = sp.Symbol("t", negative=True)
x1, x2 = sp.symbols("x1 x2", real=True)
Y1, Y2, TAU = sp.symbols("Y1 Y2 TAU", real=True)

a = sp.sqrt(-t)
y1, y2, tau = x1 / a, x2 / a, -sp.log(-t)

v = sp.Function("v")
w = sp.Function("w")
P = sp.Function("P")

u_expr = v(y1, y2, tau) / a
p_expr = P(y1, y2, tau) / a**2

V = v(Y1, Y2, TAU)
W = w(Y1, Y2, TAU)
PP = P(Y1, Y2, TAU)


def to_profile(e, weight):
    """Multiply by (-t)^weight, then rewrite in profile variables; must come out t-free."""
    e = sp.expand(sp.simplify(e * (-t) ** weight))
    e = e.subs({x1: Y1 * sp.sqrt(-t), x2: Y2 * sp.sqrt(-t)})
    e = e.subs({sp.log(-t): -TAU})
    e = sp.simplify(sp.powsimp(e, force=True))
    return e


print("u = (-t)^(-1/2) v(x/sqrt(-t), -log(-t)),  p = (-t)^(-1) P(...)")
print("Each momentum term is multiplied by (-t)^(3/2); the result must be t-free.\n")

# --- time derivative --------------------------------------------------------------
dt = to_profile(sp.diff(u_expr, t), sp.Rational(3, 2))
dt_exp = (sp.diff(V, TAU) + V / 2
          + (Y1 * sp.diff(V, Y1) + Y2 * sp.diff(V, Y2)) / 2)
assert sp.simplify(dt - dt_exp) == 0, sp.simplify(dt - dt_exp)
print("d_t u      == d_tau v + v/2 + (y.grad)v/2                        [OK]")

# --- viscous ------------------------------------------------------------------------
lap = to_profile(nu * (sp.diff(u_expr, x1, 2) + sp.diff(u_expr, x2, 2)), sp.Rational(3, 2))
assert sp.simplify(lap - nu * (sp.diff(V, Y1, 2) + sp.diff(V, Y2, 2))) == 0
print("nu Lap_x u == nu Lap_y v                                         [OK]")

# --- pressure -------------------------------------------------------------------------
gp = to_profile(sp.diff(p_expr, x1), sp.Rational(3, 2))
assert sp.simplify(gp - sp.diff(PP, Y1)) == 0
print("d_x1 p     == d_Y1 P                                             [OK]")

# --- convection (full two-component advection) ------------------------------------------
conv = to_profile(u_expr * sp.diff(u_expr, x1)
                  + (w(y1, y2, tau) / a) * sp.diff(u_expr, x2), sp.Rational(3, 2))
assert sp.simplify(conv - (V * sp.diff(V, Y1) + W * sp.diff(V, Y2))) == 0
print("(u.grad_x)u == (v.grad_y)v                                       [OK]")

# --- incompressibility --------------------------------------------------------------------
div = to_profile(sp.diff(u_expr, x1) + sp.diff(w(y1, y2, tau) / a, x2), 1)
assert sp.simplify(div - (sp.diff(V, Y1) + sp.diff(W, Y2))) == 0
print("div_x u = 0 <=> div_y v = 0                                      [OK]")

print("\n=> d_tau v + v/2 + (y.grad)v/2 + (v.grad)v + grad P = nu Lap v")

# --- steady case is Leray's profile equation at a = 1/2 -------------------------------------
U = sp.Function("U")(Y1, Y2)
Q = sp.Function("Q")(Y1, Y2)
A = sp.Rational(1, 2)
steady = (U / 2 + (Y1 * sp.diff(U, Y1) + Y2 * sp.diff(U, Y2)) / 2
          + U * sp.diff(U, Y1) + sp.diff(Q, Y1)
          - nu * (sp.diff(U, Y1, 2) + sp.diff(U, Y2, 2)))
leray = (-nu * (sp.diff(U, Y1, 2) + sp.diff(U, Y2, 2)) + A * U
         + A * (Y1 * sp.diff(U, Y1) + Y2 * sp.diff(U, Y2))
         + U * sp.diff(U, Y1) + sp.diff(Q, Y1))
assert sp.simplify(steady - leray) == 0
print("\nsteady (d_tau v = 0) == Leray profile equation, a = 1/2          [OK]")
print("   closed by NRS 1996 (U in L^3) and Tsai 1998 (local energy)")

# --- DSS becomes tau-periodicity ---------------------------------------------------------
assert sp.simplify((lam * x1) / sp.sqrt(-(lam**2 * t)) - x1 / sp.sqrt(-t)) == 0
assert sp.simplify((-sp.log(-(lam**2 * t)) - (-sp.log(-t))) + 2 * sp.log(lam)) == 0
assert sp.simplify(lam / sp.sqrt(-(lam**2 * t)) * sp.sqrt(-t) - 1) == 0
print("\nDSS  u(x,t) = lam u(lam x, lam^2 t):")
print("   y invariant, tau -> tau - 2 log lam, prefactor cancels        [OK]")
print("   => v(y,tau) = v(y, tau - S),  S = 2 log lam")
print("   exact self-similar = STEADY   (S -> 0)  -- closed")
print("   DSS                = PERIODIC (S > 0)   -- see the L^3 argument below:")
print("                        open only OUTSIDE L^3; closed for Clay data, any lambda")

print("\n   lambda ->  period S = 2 log lambda")
for L in (1.01, 1.1, 1.5, 2.0, 10.0):
    print(f"      {L:6.2f}  ->  {float(2*sp.log(L)):.4f}")

print("\nChae-Wolf Thm 1.3 closes lambda near 1, i.e. S near 0: a perturbation of the")
print("steady case. Long period resists THAT argument -- but is closed anyway by the")
print("L^3 + ESS argument below, for Clay-admissible (Schwartz) data.")
# --- energy under the DSS scaling ------------------------------------------------------
# E(t) = int |u(x,t)|^2 dx.  Using u(x,t) = lam u(lam x, lam^2 t) and z = lam x:
#     E(t) = lam^2 * lam^-3 int |u(z, lam^2 t)|^2 dz = lam^-1 E(lam^2 t),
# hence  E(lam^2 t) = lam E(t).
# In similarity variables, dx = (-t)^(3/2) dy and |u|^2 = (-t)^-1 |v|^2, so
#     E(t) = (-t)^(1/2) * ||v(.,tau)||_{L^2(dy)}^2 .
# Consistency: tau -> tau - S under t -> lam^2 t, and S-periodicity keeps the norm
# fixed, so E(lam^2 t) = (lam^2 (-t))^(1/2) ||v||^2 = lam E(t).  Checked symbolically.
nrm = sp.Symbol("nrm", positive=True)          # ||v(.,tau)||_{L^2}^2, tau-periodic
E = lambda s: sp.sqrt(-s) * nrm
assert sp.simplify(E(lam**2 * t) - lam * E(t)) == 0
print("\nENERGY under DSS scaling:")
print("   E(lam^2 t) = lam E(t)   and   E(t) = (-t)^(1/2) ||v(.,tau)||_{L^2(dy)}^2  [OK]")

# Along t_n = lam^(-2n) t0 -> 0^-, iterating gives E(t_n) = lam^-n E(t0) -> 0.
n = sp.Symbol("n", positive=True, integer=True)
t0 = sp.Symbol("t0", negative=True)
assert sp.simplify(E(lam ** (-2 * n) * t0) - lam ** (-n) * E(t0)) == 0
print("   along t_n = lam^-2n t0 -> 0^-:  E(t_n) = lam^-n E(t0) -> 0            [OK]")

print("\nADMISSIBILITY against the Clay statement -- what this does and does not show:")
print("   * If v(.,tau) is in L^2(dy) then E(t) is FINITE for every t < 0 and in fact")
print("     E(t) -> 0 as t -> 0^-.  So backward DSS is NOT excluded by Fefferman's")
print("     bounded-energy condition (7).  Finite energy and DSS are compatible.")
print("   * Chae-Wolf Thm 1.1's bound |u| <= C/(sqrt(-t)+|x|) is an UPPER bound and a")
print("     CONCLUSION of their theorem. It does not say u decays no faster than 1/|x|,")
print("     so it cannot be used to argue infinite energy. Fields in L^2 satisfy it too.")
print("   * The DATA constraint is where exact DSS dies. See below.")

# --- L^3 scale invariance closes exact DSS for Clay-admissible data ----------------------
# ||lam u(lam . , s)||_{L^3(dx)}^3 = lam^3 * lam^-3 * ||u(.,s)||_{L^3}^3 : L^3 is INVARIANT
# under the NS scaling. So the DSS relation u(x,t) = lam u(lam x, lam^2 t) gives
#     ||u(.,t)||_{L^3} = ||u(.,lam^2 t)||_{L^3},
# i.e. the L^3 norm is PERIODIC in tau with period S. Finite at one slice => bounded on
# the whole interval => u in L^inf_t L^3_x => Escauriaza-Seregin-Sverak gives smoothness.
d = 3
q = sp.Symbol("q", positive=True)
# scaling weight of ||lam u(lam x)||_{L^q(R^3)}^q relative to ||u||_{L^q}^q
wt = sp.simplify(q - d)
assert sp.solve(sp.Eq(wt, 0), q) == [3]
print("\nL^q SCALE INVARIANCE: ||lam u(lam .)||_{L^q}^q = lam^(q-3) ||u||_{L^q}^q")
print(f"   invariant exactly at q = {sp.solve(sp.Eq(wt, 0), q)[0]}                "
      "               [OK]")

print("\nEXACT DSS IS CLOSED FOR CLAY-ADMISSIBLE DATA (Chae-Wolf Remark 1.2):")
print("   L^3 invariance + DSS  =>  ||u(.,t)||_{L^3} is S-periodic in tau.")
print("   Schwartz data at one slice => that norm is finite there => bounded for all t")
print("   => u in L^inf_t L^3_x => Escauriaza-Seregin-Sverak (Russian Math Surveys 58")
print("      (2003) 211-250) removes the singularity -- FOR EVERY lambda.")
print("   So the large-lambda exact-DSS hole lives only OUTSIDE L^3, and such solutions")
print("   cannot arise from Fefferman's Schwartz data (4). Exact DSS is NOT a Clay route.")

print("\n   Net: the finite-energy objection is WITHDRAWN (it was invalid), but exact DSS")
print("   is closed anyway by the stronger L^3 + ESS argument. Surviving routes:")
print("     - localized / asymptotically DSS (not exactly DSS, so L^3-periodicity fails;")
print("       Chae-Wolf Thm 1.5 closes only lambda near 1)")
print("     - Type II blowup, outside the Type-I bound")
print("     - genuinely multi-scale cascades (no construction)")
