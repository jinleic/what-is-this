#!/usr/bin/env python3
"""THE MIXTURE THEOREM: Liu's two-protocol entropy inequality holds for EVERY
conditionally i.i.d. coupling, with the uniform constant M/m and no mean
constraint.  This is the exact input Liu's Proposition 3 needs, and it removes
both of his hypotheses (Section V-A positive semidefiniteness, Section V-B
minimiser structure) from the chain to the union-closed constant c' = 1 - m.

SETTING (natural log; x-variables, x = 1 - s as in Liu Section V-A).
    h(u) = -u log u - (1-u) log(1-u);  pi(x,y) = x y (1 + (1-x)(1-y));  K = h o pi.
A conditionally i.i.d. coupling of a law mu on [0,1] with itself is
    Gamma = int nu_u (x) nu_u dP_U(u),        mu = int nu_u dP_U(u),
for a probability space (Omega, P_U) and a Markov kernel u -> nu_u (Liu,
Definition 2 and the paragraph after Example 5: the pair (S_i, T_i) of
conditional means produced by a conditionally i.i.d. protocol is conditionally
i.i.d. given U^{i-1}).  With Example 5, f(x) = x(1-x), the protocol has
Pi_{s,t}(0,0) = pi(x,y).  Put M = <mu, x>, ehx = <mu, h> and

    numerator(Gamma) = (1-beta) <mu (x) mu, h(xy)> + beta int <nu_u (x) nu_u, K> dP_U(u).

    THEOREM (mixture).  For every conditionally i.i.d. coupling Gamma of every
    Borel probability law mu on [0,1]:
            numerator(Gamma)  >=  (M/m) * ehx.
    Equivalently gap := numerator - ehx >= ((M-m)/m) ehx.  In particular gap >= 0
    whenever M >= m (Liu H2 for every mixture), and for M >= 1-c > m the ratio
    numerator/ehx is >= (1-c)/m = (1-c)/(1-c') > 1 uniformly: this is the
    constant C > 1 that Proposition 3 requires for every c < c'.

PROOF.  Write w-sums for a finite mixture (int dP_U in general; every step
below is bilinearity, symmetry of the kernels, and Fubini with bounded
continuous integrands, so it is literally the same computation).  Let
E_kl = <nu_k (x) nu_l, h(xy)>, K_kl = <nu_k (x) nu_l, K>, H_k = <nu_k, h>,
M_k = <nu_k, x>, R = D_m/m = P2 + Q2 with P2 = (1-beta) h(xy) - (y h(x) +
x h(y))/(2m), Q2 = beta K, phi(x) = sqrt(Q2(x,x)), A = P2 + phi (x) phi,
B = Q2 - phi (x) phi, F_k = <nu_k, phi>, and P_kl = <nu_k (x) nu_l, P2>
= (1-beta) E_kl - (M_l H_k + M_k H_l)/(2m).

  (I)   sum_k w_k K_kk - sum_kl w_k w_l K_kl = (1/2) sum_kl w_k w_l (K_kk + K_ll - 2 K_kl)
        [uses sum_l w_l = 1], so
        gap = <mu(x)mu, (1-beta) h(xy) + beta K> - <mu,h> + (beta/2) sum_kl w_k w_l (K_kk+K_ll-2K_kl).
  (II)  <mu(x)mu, (1-beta) h(xy) + beta K> - <mu,h> = <mu(x)mu, R> + ((M-m)/m) <mu,h>
        [because <mu(x)mu, (y h(x) + x h(y))/(2m)> = M <mu,h>/m].
  (III) <mu(x)mu,R> + (beta/2) sum_kl w_k w_l (K_kk+K_ll-2K_kl)
        = sum_kl w_k w_l P_kl + beta sum_k w_k K_kk
        = sum_kl w_k w_l <nu_k(x)nu_l, A> + sum_k w_k <nu_k(x)nu_k, B> + [ sum_k w_k F_k^2 - (sum_k w_k F_k)^2 ].
  The last bracket is a variance, hence >= 0; A >= 0 and B >= 0 pointwise on
  [0,1]^2 (liu9-h2-twovar.json; liu9-h2-boundary.json with beta > 0), so
  (III) >= 0 and gap >= ((M-m)/m) <mu,h>.  QED.  R >= 0 is not even needed.
  The two-component case P_U = (1-q) delta_0 + q delta_1 is exactly the chain
  of liu9_h2_reduction.py / liu9_h2_general_lift.py: the variance is
  q(1-q)(F_0 - F_1)^2 and the channel is beta q(1-q) <(nu_0-nu_1)^(x)2, K>.

WHAT IS MACHINE-CHECKED HERE.
  X1  (I)+(II)+(III) as exact rational-function identities (m != 0) in the
      free pairing algebra for K = 2, 3, 4 components with symbolic weights (sum w = 1
      substituted) -- sympy residuals identically 0.  The general-K proof is the
      three displayed lines; K = 2,3,4 are the instances that exercise every term.
  X2  The same identities in Arb at random mixtures with K up to 6 components
      and up to 4 atoms per component, unequal atom counts, atoms at 0 and 1,
      to 1e-100; and numerator - (M/m) ehx = the (III) form to 1e-100.
  X3  Independent implementation: at K = 2 with three shared-mass atoms the
      numerator/ehx of this module agree with liu9_objective.evaluate_arb
      (Liu's audited transcription) to 1e-100; at K = 1 the numerator is
      (1-beta)<mu(x)mu,h(xy)> + beta<mu(x)mu,K> (Gilmer + protocol, no channel).
  X4  Witnesses.  Laws supported on {0,1}: numerator = ehx = 0 exactly (Liu's
      trivial equality case; the theorem is vacuous there, so the strictness of
      Proposition 3 must come from coordinates with ehx > 0, which is where
      C = M/m > 1 bites).  Liu's minimiser (mean exactly m): ratio = 1 = M/m.
      A five-component mixture with M > m: certified ratio > M/m > 1.
  X5  Mutations, each must be caught: drop the variance term; mixture-form
      protocol term K_kl in place of K_kk; drop the 1/2 in the channel sum;
      reversed bound numerator <= (M/m) ehx (refuted at a certified point);
      wrong protocol s+t-st against the audited evaluator; the 100-digit
      truncations declared to be the constants (refuted by W_truncated).
      Six semantic fault injections, five distinct residual modes: the
      mixture-form protocol term and the dropped 1/2 both add or remove exactly
      one channel, so they share the same error polynomial.

THE CONSTANTS.  Exactly as in liu9_h2_general_lift.py (whose `Constants` is
imported and whose module is pinned): x* is the unique root in (0,1) of
z^4 - 2z^3 + 3z^2 - 1, p* = h(x*)/h(x*^2), m* = p* x*, c' = 1 - m*, beta* by
the closed-form stationarity of Liu (89)-(90).  A >= 0 and B >= 0 are certified
at these EXACT constants only (the interior zero of R is exact and
nondegenerate), so every theorem-validity SIGN claim here (the minimiser
equality, the five-component ratio, the reversed-bound mutation) is evaluated
at Arb balls containing them; the 100-digit binding rationals are TEST points
for parameter-free identities (X2, X3, the wrong-protocol mutation) only, and
witness W_truncated deliberately uses them to record that at the truncations
the bound is certified FALSE at Liu's minimiser (by ~6e-102).  The certified
enclosure of c' is recorded; Liu's printed (93) is wrong in its last two digits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
from fractions import Fraction
from pathlib import Path
from typing import Callable, Sequence

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "FLINT_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import sympy as sp  # noqa: E402
from flint import arb, ctx  # noqa: E402

import liu9_objective  # noqa: E402
from liu9_h2_general_lift import Constants  # noqa: E402

PRECISION_BITS = 400
ctx.prec = PRECISION_BITS
DEFAULT_OUTPUT = HERE / "verification" / "results" / "liu9-h2-mixture-theorem.json"
IDENTITY_TOL = Fraction(1, 10 ** 100)
SEED = 20260903

DEPENDENCIES = {
    "liu9-h2-twovar": {
        "path": "verification/results/liu9-h2-twovar.json",
        "sha256": "a028ff17ae8492211a2ef9de5721de7f3a0afd8bd45bf7d8d9b77b7c90d2e6eb",
        "claim_status": "CERTIFIED_PROVED",
        "provides": "A >= 0 (and R >= 0) on all of [0,1]^2 with beta, m as Arb balls around the exact root-defined constants",
    },
    "liu9-h2-boundary": {
        "path": "verification/results/liu9-h2-boundary.json",
        "sha256": "e4d6d96df3ae435800cef5739fb1ae45b98febf539ab418df116ceb03bd481c9",
        "claim_status": "PROVED",
        "provides": "Phi >= 0 on all of [0,1]^2, hence B >= 0 (beta > 0)",
    },
    "liu9-h2-general-lift": {
        "path": "verification/results/liu9-h2-general-lift.json",
        "sha256": "eb4874d39904593bddc3f2ae638d1bad6756a19d4ebf6e06a60427b6bb7ac1ce",
        "claim_status": "PROVED",
        "provides": "the two-component general-measure theorem, L0 (kernel regularity), L4 (B >= 0 from Phi >= 0), and the exact-constant discipline; this module is its K-component extension",
    },
    "liu9_h2_general_lift.py": {
        "path": "liu9_h2_general_lift.py",
        "sha256": "b265e726d110098a09c00de331c4433454973ba741966b3018b77afcc225b5c3",
        "claim_status": None,
        "provides": "the shared `Constants` (exact balls from the 70-digit bracket of x*; TEST rationals) and its twovar cross-checks",
    },
    "liu9_binding.py": {
        "path": "liu9_binding.py",
        "sha256": "feb3a1aaded3cff0805d852c4349ff7142c55e24be4b2da3c0a27f8e7f58cac6",
        "claim_status": None,
        "provides": "the 70-digit rational bracket of x* and the 100-digit test solution of (87)-(90)",
    },
    "liu9_objective.py": {
        "path": "liu9_objective.py",
        "sha256": "f043b8c99cb26c62804c6f4f29190065a7a876e47d5e7c27c9e6f555c4098c6c",
        "claim_status": None,
        "provides": "h_arb and the audited nine-variable transcription evaluate_arb",
    },
}


# --------------------------------------------------------------------------- #
# Arb primitives (independent of liu9_h2_general_lift.py)
# --------------------------------------------------------------------------- #

def frac_arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def fmt(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def h(u: arb) -> arb:
    return liu9_objective.h_arb(u)


def pi_(x: arb, y: arb) -> arb:
    return x * y * (1 + (1 - x) * (1 - y))


def pi_wrong(x: arb, y: arb) -> arb:
    return x + y - x * y


Law = Sequence[tuple[Fraction, Fraction]]                 # (mass, atom)
Mixture = Sequence[tuple[Fraction, Law]]                  # (weight, component law)


def law_arb(nu: Law) -> list[tuple[arb, arb]]:
    return [(frac_arb(w), frac_arb(x)) for w, x in nu]


def pair1(nu: list[tuple[arb, arb]], f: Callable[[arb], arb]) -> arb:
    total = arb(0)
    for w, x in nu:
        total += w * f(x)
    return total


def pair2(a: list[tuple[arb, arb]], b: list[tuple[arb, arb]], k: Callable[[arb, arb], arb]) -> arb:
    total = arb(0)
    for w, x in a:
        for v, y in b:
            total += w * v * k(x, y)
    return total


def barycenter(G: Mixture) -> list[tuple[Fraction, Fraction]]:
    return [(w * m, x) for w, nu in G for m, x in nu]


def mean_exact(G: Mixture) -> Fraction:
    return sum(w * x for w, x in barycenter(G))


def numerator_ehx_arb(G: list[tuple[arb, list[tuple[arb, arb]]]], beta: arb,
                      protocol: Callable[[arb, arb], arb] = pi_) -> tuple[arb, arb]:
    """numerator(Gamma) and ehx for a finite mixture Gamma = sum_k w_k nu_k (x) nu_k (Arb data)."""
    mu = [(w * m, x) for w, nu in G for m, x in nu]
    kern = lambda x, y: h(protocol(x, y))  # noqa: E731
    ehxy = pair2(mu, mu, lambda x, y: h(x * y))
    ehpi = arb(0)
    for w, nu in G:
        ehpi += w * pair2(nu, nu, kern)
    ehx = pair1(mu, h)
    return (1 - beta) * ehxy + beta * ehpi, ehx


def numerator_ehx(G: Mixture, beta: arb, protocol: Callable[[arb, arb], arb] = pi_) -> tuple[arb, arb]:
    """numerator(Gamma) and ehx for a finite mixture with rational data."""
    return numerator_ehx_arb([(frac_arb(w), law_arb(nu)) for w, nu in G], beta, protocol)


def nonnegative_form(G: Mixture, beta: arb, m: arb) -> tuple[arb, arb, arb, arb]:
    """The (III) right-hand side: (sum w w <A>, sum w <B>, variance of F, total)."""
    comps = [(frac_arb(w), law_arb(nu)) for w, nu in G]

    def P2(x: arb, y: arb) -> arb:
        return (1 - beta) * h(x * y) - (y * h(x) + x * h(y)) / (2 * m)

    def phi(x: arb) -> arb:
        return (beta * h(pi_(x, x))).sqrt()

    def A(x: arb, y: arb) -> arb:
        return P2(x, y) + phi(x) * phi(y)

    def B(x: arb, y: arb) -> arb:
        return beta * h(pi_(x, y)) - phi(x) * phi(y)

    F = [pair1(nu, phi) for _, nu in comps]
    sumA = arb(0)
    for wk, nk in comps:
        for wl, nl in comps:
            sumA += wk * wl * pair2(nk, nl, A)
    sumB = arb(0)
    for wk, nk in comps:
        sumB += wk * pair2(nk, nk, B)
    meanF = arb(0)
    meanF2 = arb(0)
    for (wk, _), Fk in zip(comps, F):
        meanF += wk * Fk
        meanF2 += wk * Fk * Fk
    var = meanF2 - meanF * meanF
    return sumA, sumB, var, sumA + sumB + var


def abs_upper(value: arb) -> arb:
    return abs(value).upper()


# --------------------------------------------------------------------------- #
# X1: exact identities in the free pairing algebra, K components
# --------------------------------------------------------------------------- #

def mixture_identities(Kc: int, mutate: str | None = None) -> dict[str, object]:
    beta, m = sp.symbols("beta m")
    w = list(sp.symbols(f"w1:{Kc + 1}"))
    H = sp.symbols(f"H1:{Kc + 1}")
    Mk = sp.symbols(f"M1:{Kc + 1}")
    F = sp.symbols(f"F1:{Kc + 1}")

    def sym(name: str, k: int, l: int) -> sp.Symbol:
        a, b = sorted((k, l))
        return sp.Symbol(f"{name}{a + 1}{b + 1}")

    E = {(k, l): sym("E", k, l) for k in range(Kc) for l in range(Kc)}
    KK = {(k, l): sym("K", k, l) for k in range(Kc) for l in range(Kc)}
    subs = {w[-1]: 1 - sum(w[:-1])}

    Ssum = lambda f: sum(f(k) for k in range(Kc))                             # noqa: E731
    Dsum = lambda f: sum(f(k, l) for k in range(Kc) for l in range(Kc))       # noqa: E731

    Mtot = Ssum(lambda k: w[k] * Mk[k])
    ehx = Ssum(lambda k: w[k] * H[k])
    ehxy = Dsum(lambda k, l: w[k] * w[l] * E[(k, l)])
    if mutate == "mixture_protocol_term":
        ehpi = Dsum(lambda k, l: w[k] * w[l] * KK[(k, l)])
    else:
        ehpi = Ssum(lambda k: w[k] * KK[(k, k)])
    numerator = (1 - beta) * ehxy + beta * ehpi
    gap = numerator - ehx

    half = sp.Rational(1, 2) if mutate != "drop_half_in_channel" else 1
    channel = half * beta * Dsum(lambda k, l: w[k] * w[l] * (KK[(k, k)] + KK[(l, l)] - 2 * KK[(k, l)]))
    mu_mu_R = Dsum(lambda k, l: w[k] * w[l] * ((1 - beta) * E[(k, l)] + beta * KK[(k, l)]
                                               - (Mk[l] * H[k] + Mk[k] * H[l]) / (2 * m)))
    # (I)+(II): gap = <mu(x)mu,R> + channel + ((M-m)/m) ehx
    res_I_II = sp.expand((gap - (mu_mu_R + channel + (Mtot - m) / m * ehx)).subs(subs))

    P = {(k, l): (1 - beta) * E[(k, l)] - (Mk[l] * H[k] + Mk[k] * H[l]) / (2 * m) for k in range(Kc) for l in range(Kc)}
    A = {(k, l): P[(k, l)] + F[k] * F[l] for k in range(Kc) for l in range(Kc)}
    Bd = {k: beta * KK[(k, k)] - F[k] ** 2 for k in range(Kc)}
    var = Ssum(lambda k: w[k] * F[k] ** 2) - Ssum(lambda k: w[k] * F[k]) ** 2
    rhs_III = Dsum(lambda k, l: w[k] * w[l] * A[(k, l)]) + Ssum(lambda k: w[k] * Bd[k])
    if mutate != "drop_variance_term":
        rhs_III += var
    # (III): <mu(x)mu,R> + channel = sum w w A + sum w B + Var
    res_III = sp.expand((mu_mu_R + channel - rhs_III).subs(subs))
    # quantitative: numerator - (M/m) ehx = rhs_III
    res_Q = sp.expand((numerator - Mtot / m * ehx - rhs_III).subs(subs))
    ok = res_I_II == 0 and res_III == 0 and res_Q == 0
    short = lambda e: str(e) if len(str(e)) < 160 else str(e)[:160] + "..."  # noqa: E731
    return {"components": Kc, "ok": ok, "residual_I_II": short(res_I_II),
            "residual_III": short(res_III), "residual_quantitative": short(res_Q)}


def claim_X1(mutate: str | None = None) -> dict[str, object]:
    checks = [mixture_identities(Kc, mutate) for Kc in (2, 3, 4)]
    ok = all(c["ok"] for c in checks)
    return {
        "name": "X1_mixture_identities_free_algebra" if mutate is None else f"X1_mutated_{mutate}",
        "label": "MACHINE-VERIFIED" if ok else "FAILED",
        "statement": "for K = 2,3,4: gap = <mu(x)mu,R> + (beta/2) sum w_k w_l <(nu_k-nu_l)^(x)2,K> + ((M-m)/m) ehx, and <mu(x)mu,R> + channel = sum w_k w_l <nu_k(x)nu_l,A> + sum w_k <nu_k(x)nu_k,B> + Var_w(F); equivalently numerator - (M/m) ehx = that nonnegative form",
        "general_K": "the three displayed lines of the docstring proof; K = 2,3,4 exercise every term; the continuum case is the same computation with int dP_U (Fubini, bounded continuous integrands)",
        "checks": checks,
    }


# --------------------------------------------------------------------------- #
# X2/X3: Arb checks
# --------------------------------------------------------------------------- #

def random_law(rng: random.Random, n: int) -> list[tuple[Fraction, Fraction]]:
    weights = [rng.randint(1, 20) for _ in range(n)]
    total = sum(weights)
    out = []
    for wt in weights:
        roll = rng.random()
        x = Fraction(0) if roll < 0.05 else Fraction(1) if roll < 0.10 else Fraction(rng.randint(1, 999), 1000)
        out.append((Fraction(wt, total), x))
    return out


def random_mixture(rng: random.Random, Kc: int) -> list[tuple[Fraction, list[tuple[Fraction, Fraction]]]]:
    weights = [rng.randint(1, 20) for _ in range(Kc)]
    total = sum(weights)
    return [(Fraction(wt, total), random_law(rng, rng.randint(1, 4))) for wt in weights]


def claim_X2(rng: random.Random, beta: arb, m: arb, m_exact: Fraction, draws: int = 60) -> dict[str, object]:
    """Identity check at TEST parameters (the identities are parameter-free); the form's sign is a sanity sample."""
    worst = arb(0)
    min_form_lower = None
    below_threshold_seen = 0
    for _ in range(draws):
        G = random_mixture(rng, rng.randint(1, 6))
        num, ehx = numerator_ehx(G, beta)
        Mx = frac_arb(mean_exact(G))
        _, _, _, form = nonnegative_form(G, beta, m)
        residual = num - Mx / m * ehx - form
        worst = max(worst, abs_upper(residual))
        lower = form.lower()
        if min_form_lower is None or lower < min_form_lower:
            min_form_lower = lower
        if mean_exact(G) < m_exact:
            below_threshold_seen += 1
    tol = frac_arb(IDENTITY_TOL)
    ok = bool(worst <= tol)
    return {
        "name": "X2_identities_in_arb_random_mixtures",
        "label": "MACHINE-VERIFIED" if ok else "FAILED",
        "parameters": "TEST rationals (identity is parameter-free)",
        "draws": draws,
        "components_up_to": 6,
        "atoms_per_component_up_to": 4,
        "worst_abs_residual_numerator_minus_Mm_ehx_minus_form": fmt(worst, 12),
        "tolerance": str(IDENTITY_TOL),
        "min_lower_bound_of_nonnegative_form_over_sample": fmt(min_form_lower, 15) if min_form_lower is not None else None,
        "draws_with_mean_below_m": below_threshold_seen,
        "note": "the identity is certified at every draw; the form's lower bounds are a sanity sample (its nonnegativity is the theorem, from A >= 0, B >= 0 pointwise and Var >= 0)",
    }


def evaluate_nine(a: Sequence[Fraction], X: Sequence[Fraction], Y: Sequence[Fraction],
                  q: Fraction, beta: arb) -> tuple[arb, arb]:
    values = [frac_arb(a[0]), frac_arb(a[1]), frac_arb(q), frac_arb(X[0]), frac_arb(X[1]), frac_arb(X[2]),
              frac_arb(Y[0]), frac_arb(Y[1]), frac_arb(Y[2])]
    terms = liu9_objective.evaluate_arb(values, beta)
    return terms.numerator, terms.ehx


def claim_X3(rng: random.Random, beta: arb, protocol: Callable[[arb, arb], arb] = pi_, draws: int = 30) -> dict[str, object]:
    worst_num = arb(0)
    worst_ehx = arb(0)
    worst_k1 = arb(0)
    for _ in range(draws):
        a = [wt for wt, _ in random_law(rng, 3)]
        X = [x for _, x in random_law(rng, 3)]
        Y = [y for _, y in random_law(rng, 3)]
        q = Fraction(rng.randint(0, 100), 100)
        G = [(1 - q, list(zip(a, X))), (q, list(zip(a, Y)))]
        num, ehx = numerator_ehx(G, beta, protocol)
        num9, ehx9 = evaluate_nine(a, X, Y, q, beta)
        worst_num = max(worst_num, abs_upper(num - num9))
        worst_ehx = max(worst_ehx, abs_upper(ehx - ehx9))
        # K = 1: numerator is (1-beta)<mu(x)mu,h(xy)> + beta<mu(x)mu,K>
        nu = random_law(rng, rng.randint(1, 4))
        num1, _ = numerator_ehx([(Fraction(1), nu)], beta, protocol)
        mu = law_arb(nu)
        direct = (1 - beta) * pair2(mu, mu, lambda x, y: h(x * y)) + beta * pair2(mu, mu, lambda x, y: h(protocol(x, y)))
        worst_k1 = max(worst_k1, abs_upper(num1 - direct))
    tol = frac_arb(IDENTITY_TOL)
    ok = bool(worst_num <= tol) and bool(worst_ehx <= tol) and bool(worst_k1 <= tol)
    return {
        "name": "X3_agreement_with_audited_evaluator",
        "label": "MACHINE-VERIFIED" if ok else "FAILED",
        "parameters": "TEST rationals (agreement is parameter-free)",
        "draws": draws,
        "worst_abs_diff_numerator_vs_evaluate_arb": fmt(worst_num, 12),
        "worst_abs_diff_ehx_vs_evaluate_arb": fmt(worst_ehx, 12),
        "worst_abs_diff_K1_vs_direct": fmt(worst_k1, 12),
        "statement": "this module's numerator and ehx coincide with liu9_objective.evaluate_arb at three-atom shared-mass two-component mixtures (Liu's nine-parameter objective), and the one-component case is Gilmer + protocol with no channel",
    }


# --------------------------------------------------------------------------- #
# X4: witnesses
# --------------------------------------------------------------------------- #

def witnesses(C: Constants) -> list[dict[str, object]]:
    beta, m = C.beta, C.m
    out = []
    # trivial laws on {0,1} (exact zeros for every beta)
    for p in (Fraction(0), Fraction(1, 3), Fraction(1)):
        nu = [(1 - p, Fraction(0)), (p, Fraction(1))]
        G = [(Fraction(1, 2), nu), (Fraction(1, 2), [(Fraction(1), Fraction(1))])]
        num, ehx = numerator_ehx(G, beta)
        exact_zero = bool(num == 0) and bool(ehx == 0)
        out.append({"name": f"W_trivial_law_on_01_p_{p}", "label": "MACHINE-VERIFIED" if exact_zero else "FAILED",
                    "parameters": "EXACT balls", "numerator": fmt(num, 10), "ehx": fmt(ehx, 10),
                    "statement": "laws supported on {0,1}: numerator = ehx = 0 exactly (Liu's trivial equality case; h(0) = h(1) = 0, xy and pi(x,y) in {0,1})"})
    # Liu's minimiser at the EXACT balls: numerator - (M/m*) ehx = 0 is a theorem; the enclosure must contain 0
    p = C.m / C.x
    G_arb = [(arb(1), [(1 - p, arb(0)), (p, C.x)])]
    num, ehx = numerator_ehx_arb(G_arb, beta)
    resid = num - ehx                     # M = m* exactly, so (M/m*) ehx = ehx
    contains = bool(resid.contains(arb(0)))
    out.append({"name": "W_liu_minimiser_equality_exact", "label": "PROVED" if contains else "FAILED",
                "parameters": "EXACT balls", "numerator_minus_ehx_enclosure": fmt(resid, 12), "enclosure_contains_zero": contains,
                "proof_of_equality": "mean = m*; numerator - ehx = gap = I00/m* with D_{m*}(0,.) = 0 and D_{m*}(x*,x*) = m* h(x*^2) - x* h(x*) = 0 by Liu (87)-(88)",
                "statement": "at Liu's minimiser the bound numerator >= (M/m*) ehx is an equality (M = m*): the theorem is sharp; equality is proved by algebra, never numerically"})
    # the same configuration with the TEST rationals declared exact: certified NEGATIVE
    p_t = C.m_t / C.x_t
    G_t = [(Fraction(1), [(1 - p_t, Fraction(0)), (p_t, C.x_t)])]
    num_t, ehx_t = numerator_ehx(G_t, frac_arb(C.beta_t))
    resid_t = num_t - frac_arb(mean_exact(G_t)) / frac_arb(C.m_t) * ehx_t
    neg = bool(resid_t.upper() < 0)
    out.append({"name": "W_truncated_constants_certified_negative", "label": "MACHINE-VERIFIED" if neg else "FAILED",
                "parameters": "TEST rationals declared exact (beta_t, m_t, x_t)", "numerator_minus_Mm_ehx": fmt(resid_t, 12), "certified_negative": neg,
                "statement": "with the 100-digit truncations in place of the constants the bound is FALSE at Liu's minimiser (~ -6e-102): perturbing the constants cannot be assumed harmless (this joint truncation breaks the bound; raising m alone would not), which is why all sign claims use the exact balls"})
    # five-component mixture with M > m*: certified ratio above M/m* > 1
    G = [(Fraction(1, 5), [(Fraction(1, 2), Fraction(7, 10)), (Fraction(1, 2), Fraction(9, 10))]),
         (Fraction(1, 5), [(Fraction(1), Fraction(4, 5))]),
         (Fraction(1, 5), [(Fraction(1, 4), Fraction(1, 2)), (Fraction(3, 4), Fraction(1))]),
         (Fraction(1, 5), [(Fraction(1, 3), Fraction(3, 5)), (Fraction(1, 3), Fraction(7, 10)), (Fraction(1, 3), Fraction(19, 20))]),
         (Fraction(1, 5), [(Fraction(1), Fraction(3, 4))])]
    num, ehx = numerator_ehx(G, beta)
    Mx = mean_exact(G)
    ratio = num / ehx
    bound = frac_arb(Mx) / m
    ok = bool(ratio > bound) and bool(bound > 1)
    out.append({"name": "W_five_component_mixture_ratio_above_M_over_m", "label": "MACHINE-VERIFIED" if ok else "FAILED",
                "parameters": "EXACT balls", "mean": str(Mx), "ratio": fmt(ratio, 20), "M_over_m": fmt(bound, 20),
                "statement": "certified ratio > M/m > 1 at a five-component mixture with unequal atom counts and an atom at 1"})
    return out


# --------------------------------------------------------------------------- #
# X5: mutations
# --------------------------------------------------------------------------- #

def run_mutations(rng: random.Random, C: Constants, wits: list[dict[str, object]]) -> list[dict[str, object]]:
    beta, m = C.beta, C.m
    muts = []
    for mut in ("drop_variance_term", "mixture_protocol_term", "drop_half_in_channel"):
        r = claim_X1(mutate=mut)
        first_bad = next((c for c in r["checks"] if not c["ok"]), None)
        muts.append({"mutation": f"X1_{mut}", "caught": r["label"] == "FAILED",
                     "observed": "first nonzero residual at K=%s: %s" % (first_bad["components"], [first_bad[k] for k in ("residual_I_II", "residual_III", "residual_quantitative") if first_bad[k] != "0"][0][:120]) if first_bad else "no residual"})
    # reversed bound: claim numerator <= (M/m) ehx, refuted where the nonnegative form is certified positive
    G = [(Fraction(1, 2), [(Fraction(1), Fraction(4, 5))]), (Fraction(1, 2), [(Fraction(1), Fraction(9, 10))])]
    num, ehx = numerator_ehx(G, beta)
    excess = num - frac_arb(mean_exact(G)) / m * ehx
    muts.append({"mutation": "reversed_bound_numerator_le_Mm_ehx", "caught": bool(excess.lower() > 0),
                 "observed": f"numerator - (M/m*) ehx = {fmt(excess, 15)} certified > 0 at the two-point mixture (4/5, 9/10), exact balls"})
    r = claim_X3(rng, frac_arb(C.beta_t), protocol=pi_wrong, draws=5)
    muts.append({"mutation": "wrong_protocol_s+t-st_vs_evaluate_arb", "caught": r["label"] == "FAILED",
                 "observed": "worst |numerator diff| " + r["worst_abs_diff_numerator_vs_evaluate_arb"]})
    wt = next(w for w in wits if w["name"] == "W_truncated_constants_certified_negative")
    muts.append({"mutation": "declare_the_100_digit_truncations_to_be_the_constants", "caught": bool(wt["certified_negative"]),
                 "observed": "at Liu's minimiser with (beta_t, m_t, x_t): numerator - (M/m_t) ehx = " + str(wt["numerator_minus_Mm_ehx"]) + " certified < 0"})
    return muts


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #

def check_dependencies() -> tuple[bool, dict[str, object], dict[str, object]]:
    out, ok = {}, True
    twovar_constants: dict[str, object] = {}
    for name, spec in DEPENDENCIES.items():
        path = HERE / spec["path"]
        if not path.is_file():
            out[name] = {"path": spec["path"], "status": "MISSING"}
            ok = False
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        status = None
        if path.suffix == ".json":
            data = json.loads(path.read_text(encoding="utf-8"))
            status = data.get("claim_status")
            if name == "liu9-h2-twovar":
                twovar_constants = data.get("constants", {})
        match = digest == spec["sha256"] and status == spec["claim_status"]
        ok = ok and match
        out[name] = {"path": spec["path"], "sha256_expected": spec["sha256"], "sha256_on_disk": digest,
                     "claim_status_expected": spec["claim_status"], "claim_status_on_disk": status, "match": match,
                     "provides": spec["provides"]}
    return ok, out, twovar_constants


def constant_provenance(C: Constants, twovar_constants: dict[str, object]) -> dict[str, object]:
    cprime = C.cprime
    out = C.provenance(twovar_constants)
    out["sign_claims_use"] = ("EXACT balls for every theorem-validity sign claim (W_liu_minimiser_equality_exact, "
                              "W_five_component_mixture_ratio_above_M_over_m, the reversed-bound mutation, beta* > 0); "
                              "W_truncated_constants_certified_negative deliberately declares the TEST rationals exact; "
                              "X2, X3 and the wrong-protocol mutation are parameter-free checks at TEST rationals")
    out["c_prime"] = {
        "definition": "c' = 1 - m* = 1 - p* x*",
        "certified_enclosure_45": fmt(cprime, 45),
        "lower_bound_decimal_40": cprime.lower().str(40, radius=False),
        "exceeds_0.38270908791873": bool(cprime > arb("0.38270908791873")),
        "exceeds_sawin_yu_cambie_0.3823455": bool(cprime > arb("0.3823455")),
        "liu_printed_93_decimal": "0.382709087918741",
        "liu_printed_93_is_above_certified_value": bool(cprime < arb("0.382709087918741")),
        "note": "Liu's (93) prints 0.382709087918741; the root-defined constant is 0.38270908791873502993... so the printed value is wrong in its 14th-15th digits and must not be cited as the bound",
    }
    return out


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def build_report() -> dict[str, object]:
    deps_ok, deps, twovar_constants = check_dependencies()
    C = Constants()
    provenance = constant_provenance(C, twovar_constants)
    constants_ok = (C.bracket_ok and C.derivative_positive_on_bracket and C.binding_overlap
                    and all(provenance["twovar_uses_the_same_definitions"].values()) and bool(C.beta > 0))
    beta_t, m_t = frac_arb(C.beta_t), frac_arb(C.m_t)
    rng = random.Random(SEED)
    claims = [claim_X1(), claim_X2(rng, beta_t, m_t, C.m_t), claim_X3(rng, beta_t)]
    wits = witnesses(C)
    mutations = run_mutations(rng, C, wits)
    failed = [c["name"] for c in claims if c["label"] == "FAILED"] + [w["name"] for w in wits if w["label"] == "FAILED"]
    all_caught = all(mt["caught"] for mt in mutations)
    if not deps_ok:
        status = "DEPENDENCY-MISMATCH"
    elif not constants_ok:
        status = "CONSTANTS-MISMATCH"
    elif failed or not all_caught:
        status = "FAILED"
    else:
        status = "PROVED"
    theorem = {
        "statement": "at the exact constants (beta*, m*): for every Borel probability law mu on [0,1] and every conditionally i.i.d. coupling Gamma = int nu_u (x) nu_u dP_U (mu = int nu_u dP_U): (1-beta*)<mu(x)mu,h(xy)> + beta* int <nu_u(x)nu_u, h o pi> dP_U >= (M/m*) <mu,h>, M = <mu,x>",
        "label": status,
        "proof": [
            "(I) sum_k w_k K_kk - sum_kl w_k w_l K_kl = (1/2) sum_kl w_k w_l (K_kk + K_ll - 2K_kl)  [sum w = 1]",
            "(II) <mu(x)mu,(1-beta)h(xy)+beta K> - <mu,h> = <mu(x)mu,R> + ((M-m)/m)<mu,h>  [<mu(x)mu,(y h(x)+x h(y))/(2m)> = M<mu,h>/m]",
            "(III) <mu(x)mu,R> + channel = sum_kl w_k w_l <nu_k(x)nu_l,A> + sum_k w_k <nu_k(x)nu_k,B> + Var_w(<nu_k,phi>)",
            "A >= 0, B >= 0 pointwise on [0,1]^2 at the exact (beta*, m*) (pinned artifacts; beta* > 0 certified), Var >= 0; every pairing finite (kernels continuous on the compact square, general-lift L0); Fubini for a general P_U",
            "hence numerator - (M/m*) ehx >= 0; R >= 0 is not used",
            "perturbing the constants cannot be assumed harmless (W_truncated_constants_certified_negative refutes the 100-digit joint truncation); they are certified as balls, never as decimals",
        ],
        "consequences": [
            "Liu H2 for every mixture: gap >= 0 whenever M >= m",
            "uniform strictness for Liu's Proposition 3: for laws with mean(S) <= c < c' = 1 - m, i.e. M >= 1 - c > m, the ratio numerator/ehx is >= (1-c)/m > 1",
            "the trivial laws on {0,1} have numerator = ehx = 0; they do not affect Proposition 3 because H(X^n) > 0 forces some coordinate with ehx > 0",
        ],
    }
    report: dict[str, object] = {
        "artifact": "liu9-h2-mixture-theorem",
        "claim": theorem["statement"],
        "claim_status": status,
        "theorem": theorem,
        "definitions": {
            "h": "-u log u - (1-u) log(1-u), natural log", "pi": "x y (1 + (1-x)(1-y)) = Liu Example 5 with f(x) = x(1-x): Pi_{s,t}(0,0) with x = 1-s",
            "numerator": "(1-beta) ehxy + beta ehpi with ehxy the barycenter pairing of h(xy) and ehpi the P_U-average of the component pairings of h(pi)",
            "R": "D_m/m = P2 + Q2; P2 = (1-beta) h(xy) - (y h(x) + x h(y))/(2m); Q2 = beta h(pi)",
            "phi, A, B": "phi(x) = sqrt(Q2(x,x)); A = P2 + phi(x)phi; B = Q2 - phi(x)phi",
        },
        "claims": claims,
        "witnesses": wits,
        "mutations": mutations,
        "all_mutations_caught": all_caught,
        "failed_claims": failed,
        "dependencies": deps,
        "constants": provenance,
        "labels": {"PROVED": "complete argument", "MACHINE-VERIFIED": "exact symbolic residual 0 or Arb-certified comparison at the listed points",
                   "COMPUTATIONAL-EVIDENCE": "sample only", "OPEN": "not established"},
        "precision_bits": PRECISION_BITS,
        "natural_log": True,
        "seed": SEED,
        "single_core_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "FLINT_NUM_THREADS")},
        "sympy_version": sp.__version__,
        "tool": "liu9_h2_mixture_theorem.py (sympy exact identities + python-flint Arb)",
        "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "determinism": "no timestamps; fixed seed; sympy residuals as strings; Arb strings at fixed digit counts; byte-identical under the recorded single-thread environment (thread variables are set with setdefault, so a pre-set value is inherited and recorded)",
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON plus newline with report_sha256 omitted",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def write_report(report: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=1) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Liu two-protocol inequality for every conditionally i.i.d. coupling.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = build_report()
    write_report(report, args.output)
    print(f"MIXTURE_THEOREM {report['claim_status']}")
    for c in report["claims"]:
        print(f"  {c['name']:<46} {c['label']}")
    for w in report["witnesses"]:
        print(f"  {w['name']:<46} {w['label']}")
    for mt in report["mutations"]:
        print(f"  mutation {mt['mutation']:<44} caught={mt['caught']}")
    cp = report["constants"]["c_prime"]
    print(f"  c' certified enclosure {cp['certified_enclosure_45']}")
    print(f"  constants: twovar same definitions {report['constants']['twovar_uses_the_same_definitions']}")
    print(f"REPORT {args.output}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    return 0 if report["claim_status"] == "PROVED" else 3


if __name__ == "__main__":
    raise SystemExit(main())
