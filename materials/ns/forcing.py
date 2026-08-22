"""Lemma 4: an admissible force cannot alter an EXACT profile-type singularity.

DOMAIN SCOPE -- READ FIRST.  Everything here is for Fefferman alternative **(C)**, on
R^3.  It does NOT cover (D) on the torus R^3/Z^3.  Every argument runs on the Euclidean
dilation x -> lam x and the R^3 Leray multiplier, and that dilation is not a self-map of
R^3/Z^3 (it sends a 1-periodic field to a 1/lam-periodic one), so for lam != 1 it is not
a map of the fixed torus to itself.  Continuous self-similarity and DSS therefore have no
whole-space meaning there.  (D) needs its own periodic treatment, most plausibly by
rescaling around a putative singular point and passing to a whole-space limit -- not
attempted, not claimed.  (C) alone suffices for the prize.

WHY THIS FILE EXISTS.  Fefferman (C) lets the claimant choose the force, subject only to

    (5)  |d_x^alpha d_t^m f(x,t)| <= C_{alpha m K} (1 + |x| + t)^{-K}
         on R^3 x [0,inf), for any alpha, m, K,

whereas (A) and (B) say explicitly "we take f(x,t) to be identically zero".  So an
obstruction map built for the unforced equation does not automatically cover (C).  This
file establishes exactly how much of it does.

RESULT (Lemma 4).  Let f satisfy (5) and let u solve (C)'s system on R^3, with data at
t = 0 and a putative singularity at (x_0, T), T > 0.  Work in the recentred variables
xi = x - x_0, s = t - T, so t in [0,T) is s in [-T,0), and set w(xi,s) = u(x_0+xi, T+s),
ftil(xi,s) = f(x_0+xi, T+s).

  (a) If w is EXACTLY self-similar, w(xi,s) = lam w(lam xi, lam^2 s) for every lam > 0,
      then P ftil == 0.
  (b) If w is EXACTLY lam-DSS for a single lam > 1, then P ftil == 0.

  In both cases f = grad phi is a pure gradient, absorbable into the pressure, so the
  equation is EFFECTIVELY UNFORCED and the unforced rigidity theorems apply verbatim
  (Tsai's Theorem 2 imposes no pressure condition, so the shift p -> p - phi is free).

METHOD, and two retractions it replaces.

  Each term of N(u,p) = d_t u + (u.grad)u + grad p - nu Lap u carries exactly lam^3 under
  utilde = lam u(lam x, lam^2 t), ptilde = lam^2 p(lam x, lam^2 t).  Exact covariance of u
  gives utilde = u, so subtracting the two forms of the momentum equation leaves a PURE
  GRADIENT.  The Leray projector P = I - grad Lap^{-1} div has a degree-0 homogeneous
  symbol, so it annihilates gradients AND commutes with dilations:

      (P f)(x,t) = lam^3 (P f)(lam x, lam^2 t).

  Case (a): this holds for every lam > 0, so lam -> 0+ with P f bounded gives P f == 0
  in one step.  Case (b): iterate the single lam toward the singular time; the lam^{3n}
  growth contradicts boundedness.  Boundedness of P f is Lemma 4.1, proved below.

  RETRACTED #1.  An earlier draft concluded f == 0 in case (a), reasoning that exact
  self-similarity of u makes the whole left side tau^{-3/2} L(y).  FALSE: with arbitrary
  admissible forcing the PRESSURE need not share the self-similar scaling.  Counterexample:
  u == 0 is exactly self-similar, and f = grad phi with p = phi satisfies the equation
  with f not identically zero.  Only the divergence-free part of f is pinned down.

  RETRACTED #2.  An earlier draft ASSUMED an exactly-DSS force in case (b).  Fefferman
  grants no such thing; f is arbitrary subject to (5).  The covariance of P f is DERIVED
  above, not assumed.

SCOPE OF THE CONCLUSION.  Lemma 4 covers EXACT self-similarity and EXACT DSS only.  It
says nothing about ASYMPTOTICALLY self-similar or localized DSS solutions with forcing:
there u is only approximately covariant, the subtraction leaves more than a gradient, and
f need not be covariant at all.  Nor does it address wholly non-self-similar forced
blowup.  Both remain OPEN -- see README.md.
"""

import sympy as sp

nu, lam = sp.symbols("nu lam", positive=True)
al, be = sp.symbols("alpha beta", positive=True)
s = sp.Symbol("tau_s", positive=True)          # tau = T - t, taken to 0+
C = sp.Symbol("C", positive=True)              # the uniform bound from (5) with K=0

print("Fefferman (C) on R^3 permits a force f subject only to (5).")
print("(D) on the torus is OUT OF SCOPE: x -> lam x is not a self-map of R^3/Z^3.")
print("Taking alpha = m = 0, K = 0 in (5):  |f| <= C  uniformly.\n")

# --- (a) continuous self-similarity: also a Leray-projection argument -------------------
# WRONG EARLIER CLAIM. A previous draft argued: if u is exactly self-similar the whole
# left side is tau^-3/2 L(y) for a fixed L, so f = tau^-3/2 L(x/sqrt(tau)), and
# boundedness gives L == 0 hence f == 0. That is FALSE, because with an arbitrary
# admissible force the PRESSURE need not share the self-similar scaling.
#
# COUNTEREXAMPLE. u == 0 is exactly self-similar. Take any admissible scalar phi and set
# f = grad phi, p = phi. Then d_t u + (u.grad)u + grad p - nu Lap u = grad phi = f is
# satisfied with u == 0, yet f is not identically zero. So "f == 0" cannot be concluded.
# Only the DIVERGENCE-FREE PART of f is pinned down.
#
# CORRECT ARGUMENT, parallel to (b) below. Continuous self-similarity is invariance under
# u -> lam u(lam x, lam^2 t) for EVERY lam > 0 (Lemma 0). Each term of N(u,p) carries
# lam^3, so subtracting the two forms of the equation leaves a pure gradient, and applying
# the Leray projector P (degree-0 homogeneous symbol, so it kills gradients and commutes
# with dilations) gives
#       (P f)(x,t) = lam^3 (P f)(lam x, lam^2 t)      for every lam > 0.
# Send lam -> 0+ with P f bounded:  |(P f)(x,t)| <= lam^3 * ||P f||_inf -> 0.
# Hence P f == 0: f is a pure gradient, absorbable into the pressure, and the equation is
# EFFECTIVELY UNFORCED. Tsai's Theorem 2 has no pressure hypothesis, so it applies to
# u with the shifted pressure p - phi.
print("(a) CONTINUOUS SELF-SIMILARITY")
print("    Retained-term order tau^-(alpha+1):")
order = -(al + 1)
lim = sp.limit(s**order, s, 0, "+")
assert lim == sp.oo
print(f"      lim_(tau->0+) tau^({order}) = {lim}  for every alpha > 0       [OK]")
rel = sp.simplify(C / s ** (-(al + 1)))
assert sp.limit(rel, s, 0, "+") == 0
print(f"      bounded f is subcritical: C*tau^(alpha+1) -> "
      f"{sp.limit(rel, s, 0, '+')}                     [OK]")
print("    But subcriticality alone does NOT give f == 0 -- the pressure need not be")
print("    self-similar. Counterexample: u == 0 with f = grad phi, p = phi.")
bound = sp.Symbol("Pf_bound", positive=True)
decay = sp.limit(lam**3 * bound, lam, 0, "+")
assert decay == 0
print("    Leray route: (P f)(x,t) = lam^3 (P f)(lam x, lam^2 t) for EVERY lam > 0,")
print(f"    so |P f| <= lam^3 * ||P f||_inf -> {decay} as lam -> 0+           [OK]")
print("    => P f == 0: effectively unforced.  NOT f == 0.\n")

# --- (b) DSS: the force INHERITS covariance, it is not assumed ---------------------------
# The earlier draft of this file simply ASSUMED an exactly-DSS force. Fefferman grants no
# such thing: f is arbitrary subject to (5). The covariance must be derived.
#
# RECENTERING (this matters). Fefferman (C) starts the data at t = 0 and puts the
# putative singularity at some T > 0. The NS scaling symmetry is centred at a spacetime
# point, so the correct DSS relation is centred at (x_0, T), NOT at t = 0:
#
#       u(x,t) = lam u( x_0 + lam(x - x_0),  T + lam^2 (t - T) ).
#
# An earlier draft used u(x,t) = lam u(lam x, lam^2 t), which centres the symmetry on the
# INITIAL slice and formally describes a singularity at t = 0 -- before the data exists.
# Work in translated variables and the whole computation below is unchanged:
#
#       xi := x - x_0,   s := t - T,   w(xi,s) := u(x_0 + xi, T + s),
#       ftil(xi,s) := f(x_0 + xi, T + s),
#
# so t in [0,T) corresponds to s in [-T, 0), and DSS reads  w(xi,s) = lam w(lam xi, lam^2 s).
# Translation is a symmetry of NS, so w solves the same system with force ftil.
#
# Write N(u,p) = d_t u + (u.grad)u + grad p - nu Lap u, so the momentum equation for the
# translated fields is N(w,q) = ftil. Put wtil(xi,s) = lam w(lam xi, lam^2 s) and
# qtil(xi,s) = lam^2 q(lam xi, lam^2 s). Each term carries exactly lam^3 (verified below):
#       N(wtil, qtil)(xi,s) = lam^3 N(w,q)(lam xi, lam^2 s).
# If w is exactly lam-DSS then wtil = w, so
#       N(w, qtil) = lam^3 ftil(lam xi, lam^2 s),   while   N(w, q) = ftil.
# Subtracting, the two differ only by a gradient:
#       lam^3 ftil(lam xi, lam^2 s) - ftil(xi,s) = grad(qtil - q).
# Apply the Leray projection P = I - grad Lap^{-1} div, whose symbol I - xi(x)xi/|xi|^2 is
# homogeneous of degree 0 and therefore commutes with every dilation. Gradients die, so
#       (P ftil)(xi,s) = lam^3 (P ftil)(lam xi, lam^2 s).
# THAT is the covariance -- derived, for arbitrary admissible f, not assumed.
#
# ITERATION DIRECTION (also matters). Rewriting the relation as
#       (P ftil)(xi/lam^n, s/lam^{2n}) = lam^{3n} (P ftil)(xi,s),
# the sampled times s/lam^{2n} -> 0^- APPROACH the singularity and stay inside [-T,0);
# iterating the other way would send s -> -inf and leave the interval of existence.
# So if (P ftil)(xi,s) != 0 the left side is unbounded, contradicting Lemma 4.1.
x1, tt = sp.symbols("x1 t", real=True)
X, TT = sp.symbols("X T_", real=True)          # placeholders for the scaled arguments
u = sp.Function("u")
p = sp.Function("p")


def nse_terms(uu, pp, xv, tv):
    """The four momentum terms of N(u,p), in one space dimension."""
    return {
        "d_t u":     sp.diff(uu, tv),
        "(u.grad)u": uu * sp.diff(uu, xv),
        "grad p":    sp.diff(pp, xv),
        "nu Lap u":  nu * sp.diff(uu, xv, 2),
    }


# N(u,p) evaluated at the SCALED point, as a function of the placeholders.
base = nse_terms(u(X, TT), p(X, TT), X, TT)
# N(utilde, ptilde) at (x,t), then rewritten in terms of the same placeholders.
ut = lam * u(lam * x1, lam**2 * tt)
pt = lam**2 * p(lam * x1, lam**2 * tt)
tilde = nse_terms(ut, pt, x1, tt)

print("(b) DSS -- term weights under w -> lam w(lam xi, lam^2 s), recentred at (x_0,T);")
print("    the weight computation is translation-independent, so it is done locally:")
for name in base:
    lhs = tilde[name].doit()
    rhs = (lam**3 * base[name]).subs({X: lam * x1, TT: lam**2 * tt}).doit()
    assert sp.simplify(lhs - rhs) == 0, (name, sp.simplify(lhs - rhs))
    print(f"    {name:11s} carries exactly lam^3                            [OK]")
print("    => N(utilde,ptilde)(x,t) = lam^3 N(u,p)(lam x, lam^2 t)          [OK]")
print("    => for exactly DSS u:  P f  is exactly DSS with weight lam^3")
print("       (Leray symbol is degree-0 homogeneous, so it commutes with dilations)")

n = sp.Symbol("n", positive=True, integer=True)
e = sp.Symbol("e", positive=True)
assert sp.limit((1 + e) ** (3 * n), n, sp.oo) == sp.oo
print("\n    orbit toward the singularity, s_n = s_0/lam^(2n) -> 0^- inside [-T,0):")
print("      |P ftil(xi/lam^n, s_n)| = lam^(3n) |P ftil(xi,s_0)|")
for L in (1.01, 2.0, 10.0):
    print(f"      lam={L:5.2f}: lam^(3n) at n=10 is {float(L)**30:.3e}")
print("    -> oo for every lam > 1, so P ftil(xi,s_0) must vanish            [OK]")
print("    (iterating the other way would send s -> -inf, outside [-T,0))")
# --- Lemma 4.1: sup_t ||P f(.,t)||_inf < inf for f satisfying (5) -------------------------
# The DSS/self-similar closures need P f BOUNDED, uniformly in t. That is an analytic fact,
# not a symbolic identity, so the proof is written out here; sympy only checks the three
# constants that appear in it.
#
# LEMMA 4.1. Let f satisfy Fefferman's (5) on R^3 x [0,inf). Then
#     sup_{t>=0} ||P f(.,t)||_{L^inf(R^3)} < inf.
#
# PROOF. The Leray projector is the Fourier multiplier P(xi) = I - xi (x) xi / |xi|^2,
# the orthogonal projection onto xi^perp, so its operator norm is exactly 1 for xi != 0.
# With the convention fhat(xi) = int f e^{-i x.xi} dx and Fourier inversion carrying
# (2 pi)^{-3}, the inversion bound ||g||_inf <= (2 pi)^{-3} ||ghat||_{L^1} gives
#
#     ||P f(.,t)||_inf <= (2 pi)^{-3} || P(xi) fhat(xi,t) ||_{L^1_xi}
#                      <= (2 pi)^{-3} || fhat(.,t) ||_{L^1}.
#
# Split the frequency integral at |xi| = 1.
#
# LOW, |xi| <= 1.  Trivially |fhat(xi,t)| <= ||f(.,t)||_{L^1}. Apply (5) with
# alpha = 0, m = 0, K = 4:  |f(x,t)| <= C_{004}(1+|x|+t)^{-4} <= C_{004}(1+|x|)^{-4},
# and (1+|x|)^{-4} is integrable on R^3 because 4 > 3. The bound is uniform in t, so
#     A_0 := sup_t ||f(.,t)||_{L^1} <= C_{004} * int_{R^3} (1+|x|)^{-4} dx = C_{004} * 4pi/3.
# Hence int_{|xi|<=1} |fhat| dxi <= |B_1| * A_0 = (4pi/3) A_0.
#
# HIGH, |xi| >= 1.  For any xi there is a coordinate j with |xi_j| >= |xi|/sqrt(3).
# Taking alpha = 4 e_j and using  (i xi)^alpha fhat = FT(d^alpha f):
#     |xi_j|^4 |fhat(xi,t)| = |FT(d_j^4 f)(xi,t)| <= ||d_j^4 f(.,t)||_{L^1},
# so with A_4 := sup_t max_{|alpha|=4} ||d^alpha f(.,t)||_{L^1}, finite by (5) with
# |alpha| = 4 and K = 4 (same integrability, uniform in t),
#     |fhat(xi,t)| <= (sqrt(3))^4 A_4 |xi|^{-4} = 9 A_4 |xi|^{-4}.
# Then int_{|xi|>=1} |xi|^{-4} dxi = 4 pi int_1^inf r^{-2} dr = 4 pi, so
#     int_{|xi|>=1} |fhat| dxi <= 36 pi A_4.
#
# Therefore sup_t ||fhat(.,t)||_{L^1} <= (4pi/3) A_0 + 36 pi A_4 < inf, and
# sup_t ||P f(.,t)||_inf <= (2 pi)^{-3} [ (4pi/3) A_0 + 36 pi A_4 ] < inf.   QED
#
#
# TRANSLATION INVARIANCE. Lemma 4 is applied to ftil(xi,s) = f(x_0+xi, T+s). The bound
# above survives verbatim: A_0 and A_4 are SPATIAL L^1 norms, invariant under the shift
# x -> x_0 + xi, and (5) is uniform over t >= 0, hence uniform over s in [-T,0). So
#     sup_{s in [-T,0)} ||P ftil(.,s)||_inf <= (2 pi)^{-3} [ (4pi/3) A_0 + 36 pi A_4 ],
# the same constant. No new hypothesis is needed for the recentred statement.
# Note the lemma is genuinely about (5)'s FULL strength: P is not bounded on L^inf, so
# |f| <= C alone would not suffice. Decay plus four derivatives is what does it.
r = sp.Symbol("r", positive=True)
ball = sp.simplify(4 * sp.pi * sp.integrate(r**2, (r, 0, 1)))
assert ball == sp.Rational(4, 3) * sp.pi, ball
tail = sp.simplify(4 * sp.pi * sp.integrate(r**2 * r**-4, (r, 1, sp.oo)))
assert tail == 4 * sp.pi, tail
decay_l1 = sp.simplify(4 * sp.pi * sp.integrate(r**2 * (1 + r) ** -4, (r, 0, sp.oo)))
assert decay_l1 == sp.Rational(4, 3) * sp.pi, decay_l1
assert sp.sqrt(3) ** 4 == 9
print("    Lemma 4.1 (proved in the source comment): sup_t ||P f(.,t)||_inf < inf.")
print(f"      |B_1| = {ball},  int_(|xi|>=1)|xi|^-4 dxi = {tail},")
print(f"      int_(R^3)(1+|x|)^-4 dx = {decay_l1},  (sqrt 3)^4 = 9        [OK]")
print("      => ||P f||_inf <= (2pi)^-3 [ (4pi/3) A_0 + 36 pi A_4 ],  A_0,A_4 from (5)")
print("      P is NOT bounded on L^inf, so |f| <= C alone would not suffice;")
print("      the proof uses (5)'s decay AND its fourth derivatives.")

# --- scope statement --------------------------------------------------------------------
print("SCOPE -- Fefferman (C) on R^3 ONLY, stated precisely:")
print("  CLOSED    - EXACTLY self-similar or EXACTLY DSS u with any admissible f.")
print("              Such a u forces P f == 0, so the rigidity theorems apply unchanged.")
print("  NOT CLOSED- ASYMPTOTICALLY self-similar / localized DSS with forcing. There u is")
print("              only approximately covariant, so the argument above yields nothing")
print("              about f, and f need not be DSS at all.")
print("  NOT CLOSED- wholly non-self-similar forced blowup. Nothing constructed, nothing")
print("              excluded. No forced 3D NS blowup construction exists in the")
print("              literature; known forced results are Euler/Boussinesq with forces")
print("              far rougher than (5), and Albritton-Brue-Colombo is weak")
print("              non-uniqueness rather than breakdown.")
print("  OUT OF SCOPE - Fefferman (D) on the torus R^3/Z^3. The Euclidean dilation")
print("              x -> lam x is not a self-map of R^3/Z^3, so the covariance argument")
print("              has no torus analogue as written. (D) needs its own periodic")
print("              treatment; none is claimed here.")
print("\n  NRS and ESS primaries: retrieved and verified verbatim 2026-08-15")
print("  (ns/sources/NRS-1996-acta.pdf, ns/sources/ESS-2003-ima1904.pdf; see ns/README.md")
print("  'Verification status'). No [UNVERIFIED] tag remains on this map. All of")
print("  NRS/Tsai/Chae/Chae-Wolf/ESS are stated for f == 0; the bridge above is what")
print("  carries them to the forced EXACT case.")
