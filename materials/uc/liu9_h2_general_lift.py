#!/usr/bin/env python3
"""LIU H2 FOR GENERAL MEASURES: the paired-class theorem lifts to every pair of
Borel probability measures on [0,1], with no density argument.

STATEMENT.  Natural log throughout.  h(u) = -u log u - (1-u) log(1-u) with
h(0) = h(1) = 0; pi(s,t) = s t (1 + (1-s)(1-t)) (Liu's protocol); K = h o pi.
For Borel probability measures mu, nu on [0,1] and q in [0,1] put

    P   = (1-q) mu + q nu                                  (the mixture law)
    M   = <P, x>                                           (its mean)
    gap(mu,nu,q) = (1-beta) <P (x) P, h(xy)>
                   + beta [ (1-q) <mu (x) mu, K> + q <nu (x) nu, K> ]
                   - <P, h>.

At finitely supported mu, nu this is exactly `liu9_objective._formula`
(numerator minus ehx): ehxy is the mixture pairing, ehpi is the COMPONENT-wise
protocol pairing, ehx the mixture entropy.  Liu's mean constraint (the sole
nonlinear constraint of frankl5.m, `SOURCE_MEAN_LINE`) is M >= 1 - c = m.

THE CONSTANTS -- exact definitions, Liu (87)-(90).  Everything is stated at
    x*    the unique root in (0,1) of f(z) = z^4 - 2 z^3 + 3 z^2 - 1      (87)
    p*    = h(x*) / h(x*^2)                                                (88)
    m*    = p* x*,   c' = 1 - m*
    beta* = ((h'(x*) + h(x*)/x*)/p* - a') / (b' - a'),
            a' = h'(x*^2) 2x*,  b' = h'(pi(x*,x*)) d/dx pi(x,x)|_{x*}      (89)-(90)
These are the constants at which the pinned inputs certify R >= 0 and A >= 0
(`liu9-h2-twovar.json`, `constants` block: the same quartic, the same closed
form for m, the same beta, with the interior zero (x*,x*) handled by EXACT
algebra: R(x*,x*) = 0 because h(pi(x*,x*)) = h(x*^2), and grad R(x*,x*) = 0 by
the definition of beta*).  R >= 0 is a TIGHT statement: the zero at (x*,x*) is
nondegenerate, so perturbing the constants cannot be assumed harmless (lowering
m, or moving beta at fixed m, makes R negative near (x*,x*) by first order;
raising m keeps R >= 0 by L2).  Consequently this module keeps
two parameter sets and never confuses them:
  * EXACT balls: Arb enclosures of x*, p*, m*, beta*, c' propagated from a
    70-digit rational bracket of x* (re-certified here: f(lo) < 0 < f(hi) and
    f' > 0 on (0,1]) through the exact formulas above.  Every theorem-validity
    SIGN claim in this artifact (W1a, W2, W3, the sample, the sign mutations)
    is evaluated with these balls, so it holds at the exact constants.
  * TEST rationals: the 100-digit binding solution converted exactly from
    mpf.  They are within 1e-100 of the constants and are used ONLY as test
    points for parameter-free identities (which hold for every beta, m).  At
    a tight point they are NOT the constants: declaring them exact makes the
    theorem false by ~6e-102 at Liu's minimiser (witness W1b, a certified
    negative, found by the independent review and recorded here on purpose).

    THEOREM (general measures).  For all Borel probability measures mu, nu on
    [0,1] and all q in [0,1] with M >= m*:   gap(mu,nu,q) >= 0, with beta = beta*.

Liu's Hypothesis 2 (Section V-B, nine parameters: three shared masses, six
supports, q) is the special case mu = sum a_i delta_{x_i}, nu = sum a_i
delta_{y_i}; the paired class with arbitrarily many atoms (PROGRESS.md entry
LIU-H2-PAIRED-CLASS-PROVED-AND-AUDITED) is the case of shared masses.

WHY NO DENSITY ARGUMENT IS NEEDED.  The reduction chain C1/C5/C6 of
`liu9_h2_reduction.py` is an identity in the free algebra of pairings: it uses
only (i) bilinearity of (mu,nu) -> <mu (x) nu, k>, (ii) symmetry k(s,t) =
k(t,s) of every kernel, and (iii) the product rule <mu (x) nu, f(s) g(t)> =
<mu,f><nu,g>.  All three hold for arbitrary Borel probability measures as soon
as the kernels are bounded Borel functions on [0,1]^2 -- and every kernel here
is CONTINUOUS on the compact square (L0).  Integrating the certified pointwise
inequalities R >= 0, A >= 0, B >= 0 against product probability measures then
gives the theorem outright.  Weak density of finitely supported laws and the
(true) weak continuity of gap are recorded as a corollary but carry no weight,
so no modulus of continuity of h -- which is not Lipschitz at 0 or 1 -- ever
enters.

THE PROOF, CLAIM BY CLAIM (labels: PROVED / MACHINE-VERIFIED /
COMPUTATIONAL-EVIDENCE / OPEN).

  L0  Regularity.  h is continuous on [0,1] with 0 <= h <= log 2 (h'' =
      -1/(u(1-u)) < 0, h(u) = h(1-u), h(1/2) = log 2, u log u -> 0); pi maps
      [0,1]^2 into [0,1] because 1 - pi(1-a,1-b) = a(1-b) + b(1-a) +
      ab(a+b-ab) is a sum of nonnegative terms.  Hence h(xy), K, D_M, R, P2,
      Q2, phi = sqrt(Q2(s,s)), A = P2 + phi (x) phi, B = Q2 - phi (x) phi are
      continuous on [0,1]^2 and every pairing below is finite.  PROVED
      (elementary); the identities are checked by sympy.
  L1  Free-algebra form of C1.  With E_ab = <.,h(xy)>, K_ab = <.,K>, H_a =
      <.,h>, M_a = <.,x> for a,b in {0,1} (0 = mu, 1 = nu), I_ab = <.,D_M>,
      M = (1-q) M_0 + q M_1 > 0, c = beta <(mu-nu)(x)(mu-nu), K>:
          gap = ((1-q)^2 I_00 + q^2 I_11)/M + q(1-q) T,   T = 2 I_01/M + c.
      MACHINE-VERIFIED: sympy residual identically 0 (beta, m free symbols).
  L2  Monotonicity in the mean.  D_M/M - D_m/m = (t h(s) + s h(t))(M-m)/(2Mm)
      >= 0 pointwise when M >= m > 0.  PROVED (sympy identity + signs).
  L3  Rank-one splitting (C6) in the free algebra.  T_m := c + 2 <mu (x) nu, R>
      = <mu(x)mu, B> + <nu(x)nu, B> + 2 <mu(x)nu, A> + (<mu,phi> - <nu,phi>)^2.
      MACHINE-VERIFIED: sympy residual identically 0.
  L4  B >= 0 from Phi >= 0.  Phi = K_st^2 - K_ss K_tt = (K_st - sqrt(K_ss
      K_tt))(K_st + sqrt(K_ss K_tt)) and B = beta (K_st - sqrt(K_ss K_tt));
      with K >= 0 and beta > 0, Phi >= 0 forces B >= 0 (if the second factor
      vanishes then K_st = 0 = K_ss K_tt and B = 0).  PROVED (elementary).
  L5  Inputs, pinned by sha256 and claim_status, and checked to use THE SAME
      constant definitions: R >= 0 and A >= 0 on [0,1]^2 at (beta*, m*)
      (`liu9-h2-twovar.json`), Phi >= 0 on [0,1]^2, parameter-free
      (`liu9-h2-boundary.json`), the chain (`liu9-h2-reduction.json`,
      re-derived here symbolically).  `liu9_binding.py` and
      `liu9_objective.py` are pinned as well.
  (a) Product coupling.  For mu = sum a_i delta_{x_i} (n atoms) and nu = sum
      b_j delta_{y_j} (m atoms), the nm-atom PAIRED configuration with shared
      masses g_ij = a_i b_j, x-supports x_i, y-supports y_j has exactly the
      same gap and the same mean, because the paired formula depends on the
      configuration only through its marginals and sum_j b_j = sum_i a_i = 1.
      MACHINE-VERIFIED by sympy for (n,m) in {(1,1),(1,2),(2,1),(2,2),(2,3),
      (3,2),(3,3)} with the function values as free symbols (a polynomial
      identity in the masses), and by Arb against `liu9_objective.evaluate_arb`
      at the sizes it accepts.  The general-(n,m) argument is the same one-line
      computation.  PROVED.
  (b) General measures.  Compose L1-L5 at (beta*, m*):
          gap = (1-q)^2 <mu(x)mu, D_M/M> + q^2 <nu(x)nu, D_M/M>
                + q(1-q) [ 2 <mu(x)nu, D_M/M> + c ]
              >= (1-q)^2 <mu(x)mu, R> + q^2 <nu(x)nu, R> + q(1-q) T_m     (L2)
              >= q(1-q) [ <mu(x)mu,B> + <nu(x)nu,B> + 2<mu(x)nu,A> + (<mu,phi>-<nu,phi>)^2 ]   (R>=0; L3)
              >= 0                                                        (A >= 0, B >= 0 pointwise).
      PROVED, conditional only on the pinned artifacts.

SHARPNESS AND THE ROLE OF THE CONSTANTS (witnesses).
  W1a  gap = 0 at Liu's minimiser mu = nu = (1 - m*/x*) delta_0 + (m*/x*)
       delta_{x*} (mean exactly m*): PROVED by algebra (D_{m*}(0,.) = 0 and
       D_{m*}(x*,x*) = m* h(x*^2) - x* h(x*) = 0 by (87)-(88)); the Arb
       enclosure at the EXACT balls contains 0.  Equality is never decided
       numerically.
  W1b  the same configuration with the TEST rationals declared exact has a
       CERTIFIED NEGATIVE gap (~ -6e-102): the theorem is false for these
       100-digit truncations.  (The independent review's counterexample;
       the tight point is exact.)
  W2/W3  gap < 0 certified at laws with mean below m* (paired and non-paired):
       the mean constraint is load-bearing.

MUTATIONS (all must be caught): mixture-form protocol term in L1; dropped
square in L3; wrong coupling masses in (a); the unconstrained claim (refuted
by a certified negative gap); the wrong protocol s+t-st against the audited
evaluator; reversed L2 monotonicity (refuted at a certified point); the TEST
rationals declared to be the constants (refuted by W1b).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import re
import sys
from decimal import Decimal
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
from liu9_binding import certify_equation_parameters, solve_equation_parameters  # noqa: E402

PRECISION_BITS = 400
ctx.prec = PRECISION_BITS
DEFAULT_OUTPUT = HERE / "verification" / "results" / "liu9-h2-general-lift.json"
IDENTITY_TOL = Fraction(1, 10 ** 100)   # Arb agreement between evaluation routes at TEST points
SEED = 20260902

# Pinned inputs.  A mismatch fails the artifact instead of silently composing
# with a changed input.
DEPENDENCIES = {
    "liu9-h2-twovar": {
        "path": "verification/results/liu9-h2-twovar.json",
        "sha256": "a028ff17ae8492211a2ef9de5721de7f3a0afd8bd45bf7d8d9b77b7c90d2e6eb",
        "claim_status": "CERTIFIED_PROVED",
        "provides": "R >= 0 and A >= 0 on all of [0,1]^2 at the exact (beta*, m*) (400-bit Arb, exhaustive strata; exact algebra at (x*,x*))",
    },
    "liu9-h2-boundary": {
        "path": "verification/results/liu9-h2-boundary.json",
        "sha256": "e4d6d96df3ae435800cef5739fb1ae45b98febf539ab418df116ceb03bd481c9",
        "claim_status": "PROVED",
        "provides": "Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s)) h(pi(t,t)) >= 0 on all of [0,1]^2 (parameter-free; kappa >= 9/4)",
    },
    "liu9-h2-reduction": {
        "path": "verification/results/liu9-h2-reduction.json",
        "sha256": "c63ae94e4fbd2397e19ea9e2a4e2c155f1e02265852d0aacad5c97d70f078e3a",
        "claim_status": "REDUCTION-CERTIFIED",
        "provides": "the chain C1-C5 at the nine-variable level (re-derived here in the free pairing algebra, so cited, not load-bearing)",
    },
    "liu9_binding.py": {
        "path": "liu9_binding.py",
        "sha256": "feb3a1aaded3cff0805d852c4349ff7142c55e24be4b2da3c0a27f8e7f58cac6",
        "claim_status": None,
        "provides": "the 70-digit rational bracket of x* and the 100-digit test solution of (87)-(90); its formulas are re-implemented here and cross-checked",
    },
    "liu9_objective.py": {
        "path": "liu9_objective.py",
        "sha256": "f043b8c99cb26c62804c6f4f29190065a7a876e47d5e7c27c9e6f555c4098c6c",
        "claim_status": None,
        "provides": "h_arb and the audited nine-variable transcription evaluate_arb",
    },
}
TWOVAR_EXPECTED_CONSTANT_DEFINITIONS = {
    "xstar_defining_polynomial": "x^4-2*x^3+3*x^2-1",
    "m_closed_form": "xstar*h(xstar)/h(xstar^2)",
}


# --------------------------------------------------------------------------- #
# Arb primitives
# --------------------------------------------------------------------------- #

def exact_mpf_fraction(value: object) -> Fraction:
    sign, mantissa, exponent, _bits = value._mpf_  # type: ignore[attr-defined]
    signed = -mantissa if sign else mantissa
    if exponent >= 0:
        return Fraction(signed * (1 << exponent))
    return Fraction(signed, 1 << (-exponent))


def frac_arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def fmt(value: arb, digits: int = 30) -> str:
    return value.str(digits, radius=True)


def h(u: arb) -> arb:
    return liu9_objective.h_arb(u)


def hp(u: arb) -> arb:
    """h'(u) = log((1-u)/u) on (0,1)."""
    return ((1 - u) / u).log()


def pi_(s: arb, t: arb) -> arb:
    return s * t * (1 + (1 - s) * (1 - t))


def pi_wrong(s: arb, t: arb) -> arb:
    """Mutation only: the protocol s + t - st that Liu does NOT use."""
    return s + t - s * t


def hxy(s: arb, t: arb) -> arb:
    return h(s * t)


def K(s: arb, t: arb) -> arb:
    return h(pi_(s, t))


def D_M(s: arb, t: arb, M: arb, beta: arb) -> arb:
    return M * ((1 - beta) * hxy(s, t) + beta * K(s, t)) - (t * h(s) + s * h(t)) / 2


def R_kernel(s: arb, t: arb, beta: arb, m: arb) -> arb:
    return D_M(s, t, m, beta) / m


def abs_upper(value: arb) -> arb:
    return abs(value).upper()


def certified_below(value: arb, bound: arb) -> bool:
    return bool(value.upper() < bound)


def certified_above(value: arb, bound: arb) -> bool:
    return bool(value.lower() > bound)


# --------------------------------------------------------------------------- #
# The constants
# --------------------------------------------------------------------------- #

class Constants:
    """EXACT balls (x, p, m, beta, cprime) and TEST rationals (x_t, m_t, beta_t)."""

    def __init__(self) -> None:
        mp = solve_equation_parameters(100)
        certified = certify_equation_parameters(mp)
        self.root_lo, self.root_hi = certified.root_lo, certified.root_hi
        lo, hi = frac_arb(self.root_lo), frac_arb(self.root_hi)
        f = lambda z: z ** 4 - 2 * z ** 3 + 3 * z ** 2 - 1  # noqa: E731
        self.bracket_ok = bool(f(lo) < 0) and bool(f(hi) > 0)
        x = lo.union(hi)
        self.derivative_positive_on_bracket = bool((4 * x ** 3 - 6 * x ** 2 + 6 * x) > 0)
        self.x = x
        self.p = h(x) / h(x * x)
        self.m = self.p * x
        self.cprime = 1 - self.m
        a_prime = hp(x * x) * (2 * x)
        t = pi_(x, x)
        dt = 2 * x + 2 * x * (1 - x) ** 2 - 2 * x * x * (1 - x)
        b_prime = hp(t) * dt
        self.beta = ((hp(x) + h(x) / x) / self.p - a_prime) / (b_prime - a_prime)
        # cross-checks against the binding's own propagation
        self.binding_overlap = (bool(certified.x.overlaps(self.x)) and bool(certified.p.overlaps(self.p))
                                and bool(certified.mean.overlaps(self.m)) and bool(certified.beta.overlaps(self.beta)))
        # TEST rationals: the stored 100-digit mpf, converted exactly
        self.x_t = exact_mpf_fraction(mp.x)
        self.m_t = exact_mpf_fraction(mp.mean)
        self.beta_t = exact_mpf_fraction(mp.beta)
        self.c_t = exact_mpf_fraction(mp.c)

    def provenance(self, twovar_constants: dict[str, object]) -> dict[str, object]:
        x_res_t = self.x_t ** 4 - 2 * self.x_t ** 3 + 3 * self.x_t ** 2 - 1
        return {
            "definitions": {
                "x*": "unique root in (0,1) of z^4 - 2 z^3 + 3 z^2 - 1 (Liu (87)); f(0) = -1 < 0 < f(1) = 1 and f'(z) = 2z(2(z-3/4)^2 + 15/8) > 0 on (0,1]",
                "p*": "h(x*)/h(x*^2) (Liu (88))",
                "m*": "p* x*;  c' = 1 - m*",
                "beta*": "((h'(x*) + h(x*)/x*)/p* - a')/(b' - a'), a' = h'(x*^2) 2x*, b' = h'(pi(x*,x*)) (2x + 2x(1-x)^2 - 2x^2(1-x))|_{x*} (Liu (89)-(90), the binding's closed form)",
            },
            "exact_bracket_of_xstar": {"lo": str(self.root_lo), "hi": str(self.root_hi),
                                       "width": str(self.root_hi - self.root_lo),
                                       "f_lo_negative_and_f_hi_positive_certified": self.bracket_ok,
                                       "f_prime_positive_on_bracket_certified": self.derivative_positive_on_bracket},
            "exact_balls_45_digits": {"x*": fmt(self.x, 45), "p*": fmt(self.p, 45), "m*": fmt(self.m, 45),
                                      "beta*": fmt(self.beta, 45), "c'": fmt(self.cprime, 45)},
            "binding_propagation_overlaps_own_propagation": self.binding_overlap,
            "twovar_uses_the_same_definitions": {
                "polynomial": twovar_constants.get("xstar_defining_polynomial") == TWOVAR_EXPECTED_CONSTANT_DEFINITIONS["xstar_defining_polynomial"],
                "m_closed_form": twovar_constants.get("m_closed_form") == TWOVAR_EXPECTED_CONSTANT_DEFINITIONS["m_closed_form"],
                "beta_enclosure_overlaps": overlaps_enclosure_string(self.beta, twovar_constants.get("beta_enclosure")),
                "m_enclosure_overlaps": overlaps_enclosure_string(self.m, twovar_constants.get("m_enclosure")),
                "xstar_enclosure_overlaps": overlaps_enclosure_string(self.x, twovar_constants.get("xstar_enclosure")),
                "xstar_bracket_overlaps_twovar_bracket": brackets_overlap(self.root_lo, self.root_hi,
                                                                          twovar_constants.get("xstar_root_lo"), twovar_constants.get("xstar_root_hi")),
            },
            "test_rationals": {
                "source": "liu9_binding.solve_equation_parameters(100), stored mpf converted exactly (sign, mantissa, exponent)",
                "role": "test points for parameter-free identities ONLY; within 1e-100 of the constants; NOT the constants (see W1b)",
                "x_t": str(self.x_t), "m_t": str(self.m_t), "beta_t": str(self.beta_t),
                "one_minus_c_t_minus_m_t_is_zero": (1 - self.c_t - self.m_t) == 0,
                "quartic_residual_at_x_t": fmt(frac_arb(x_res_t), 12),
                "x_t_inside_exact_ball": bool(self.x.contains(frac_arb(self.x_t))),
                "beta_t_inside_exact_ball": bool(self.beta.contains(frac_arb(self.beta_t))),
                "m_t_inside_exact_ball": bool(self.m.contains(frac_arb(self.m_t))),
            },
            "sign_claims_use": "EXACT balls for every theorem-validity sign claim (W1a, W2, W3, the sample, the drop-mean-constraint and reversed-monotonicity mutations); W1b is deliberately the TEST rationals declared exact; parameter-free identity checks (arb_crosschecks, the wrong-protocol mutation) use TEST rationals",
        }


def parse_enclosure(text: object) -> tuple[Fraction, Fraction] | None:
    if not isinstance(text, str):
        return None
    match = re.fullmatch(r"\[(-?[0-9.]+) \+/- ([0-9.]+(?:e[+-]?[0-9]+)?)\]", text.strip())
    if not match:
        return None
    return Fraction(Decimal(match.group(1))), Fraction(Decimal(match.group(2)))


def overlaps_enclosure_string(ball: arb, text: object) -> bool:
    parsed = parse_enclosure(text)
    if parsed is None:
        return False
    mid, rad = parsed
    theirs = frac_arb(mid - rad).union(frac_arb(mid + rad))
    return bool(ball.overlaps(theirs))


def brackets_overlap(lo: Fraction, hi: Fraction, other_lo: object, other_hi: object) -> bool:
    try:
        olo, ohi = Fraction(str(other_lo)), Fraction(str(other_hi))
    except (ValueError, TypeError):
        return False
    return max(lo, olo) <= min(hi, ohi)


# --------------------------------------------------------------------------- #
# Measures and the functional
# --------------------------------------------------------------------------- #

Measure = Sequence[tuple[Fraction, Fraction]]   # (mass, support), masses sum to 1
ArbMeasure = list[tuple[arb, arb]]


def measure_arb(mu: Measure) -> ArbMeasure:
    return [(frac_arb(w), frac_arb(x)) for w, x in mu]


def pair1(mu: ArbMeasure, f: Callable[[arb], arb]) -> arb:
    total = arb(0)
    for w, x in mu:
        total += w * f(x)
    return total


def pair2(mu: ArbMeasure, nu: ArbMeasure, k: Callable[[arb, arb], arb]) -> arb:
    total = arb(0)
    for w, x in mu:
        for v, y in nu:
            total += w * v * k(x, y)
    return total


def mixture_arb(A: ArbMeasure, B: ArbMeasure, qa: arb) -> ArbMeasure:
    return [((1 - qa) * w, x) for w, x in A] + [(qa * v, y) for v, y in B]


def gap_general_arb(A: ArbMeasure, B: ArbMeasure, qa: arb, beta: arb,
                    protocol: Callable[[arb, arb], arb] = pi_) -> arb:
    P = mixture_arb(A, B, qa)
    kern = lambda s, t: h(protocol(s, t))  # noqa: E731
    ehxy = pair2(P, P, hxy)
    ehpi = (1 - qa) * pair2(A, A, kern) + qa * pair2(B, B, kern)
    ehx = pair1(P, h)
    return (1 - beta) * ehxy + beta * ehpi - ehx


def gap_general(mu: Measure, nu: Measure, q: Fraction, beta: arb,
                protocol: Callable[[arb, arb], arb] = pi_) -> arb:
    """The general-measure functional at finitely supported mu, nu (rational data)."""
    return gap_general_arb(measure_arb(mu), measure_arb(nu), frac_arb(q), beta, protocol)


def mean_exact(mu: Measure, nu: Measure, q: Fraction) -> Fraction:
    return sum((1 - q) * w * x for w, x in mu) + sum(q * v * y for v, y in nu)


def gap_chain_rhs(mu: Measure, nu: Measure, q: Fraction, beta: arb) -> arb:
    """Right-hand side of L1 evaluated directly from the D_M and K pairings."""
    A, B = measure_arb(mu), measure_arb(nu)
    qa = frac_arb(q)
    M = frac_arb(mean_exact(mu, nu, q))
    d = lambda s, t: D_M(s, t, M, beta)  # noqa: E731
    I00, I11, I01 = pair2(A, A, d), pair2(B, B, d), pair2(A, B, d)
    c = beta * (pair2(A, A, K) + pair2(B, B, K) - 2 * pair2(A, B, K))
    T = 2 * I01 / M + c
    return ((1 - qa) ** 2 * I00 + qa ** 2 * I11) / M + qa * (1 - qa) * T


def gap_paired_atoms(g: Sequence[Fraction], X: Sequence[Fraction], Y: Sequence[Fraction],
                     q: Fraction, beta: arb) -> arb:
    """The N-atom paired formula: the literal generalisation of liu9_objective._formula."""
    return gap_general(list(zip(g, X)), list(zip(g, Y)), q, beta)


def product_coupling(mu: Measure, nu: Measure) -> tuple[list[Fraction], list[Fraction], list[Fraction]]:
    g, X, Y = [], [], []
    for a, x in mu:
        for b, y in nu:
            g.append(a * b)
            X.append(x)
            Y.append(y)
    return g, X, Y


# --------------------------------------------------------------------------- #
# Symbolic claims (sympy, exact).  beta and m are free symbols here.
# --------------------------------------------------------------------------- #

def claim_L0_regularity() -> dict[str, object]:
    u = sp.Symbol("u", positive=True)
    a, b = sp.symbols("a b")
    hh = -u * sp.log(u) - (1 - u) * sp.log(1 - u)
    second = sp.simplify(sp.diff(hh, u, 2) + 1 / (u * (1 - u)))
    symmetric = sp.simplify(hh - hh.subs(u, 1 - u))
    at_half = sp.simplify(hh.subs(u, sp.Rational(1, 2)) - sp.log(2))
    lim0 = sp.limit(hh, u, 0, "+")
    lim1 = sp.limit(hh, u, 1, "-")
    pi_ab = (1 - a) * (1 - b) * (1 + a * b)          # pi(1-a, 1-b)
    complement = a * (1 - b) + b * (1 - a) + a * b * (a + b - a * b)
    pi_range = sp.expand(1 - pi_ab - complement)
    z = sp.Symbol("z")
    quartic_derivative = sp.expand(sp.diff(z ** 4 - 2 * z ** 3 + 3 * z ** 2 - 1, z) - 2 * z * (2 * (z - sp.Rational(3, 4)) ** 2 + sp.Rational(15, 8)))
    ok = all(expr == 0 for expr in (second, symmetric, at_half, lim0, lim1, pi_range, quartic_derivative))
    return {
        "name": "L0_regularity",
        "label": "PROVED" if ok else "FAILED",
        "statement": "h continuous on [0,1], 0 <= h <= log 2; pi([0,1]^2) subset [0,1]; hence h(xy), K, D_M, R, P2, Q2, phi, A, B continuous on [0,1]^2 and every pairing finite; f'(z) = 2z(2(z-3/4)^2+15/8) so the quartic has exactly one root in (0,1)",
        "sympy": {
            "h_second_derivative_plus_1_over_u_1mu": str(second),
            "h_minus_h_of_1mu": str(symmetric),
            "h_half_minus_log2": str(at_half),
            "limit_h_at_0_plus": str(lim0),
            "limit_h_at_1_minus": str(lim1),
            "one_minus_pi_1ma_1mb_minus_nonnegative_sum": str(pi_range),
            "quartic_derivative_minus_2z_2_z_minus_3_4_sq_plus_15_8": str(quartic_derivative),
        },
        "why_bounded": "concave (h'' < 0), symmetric, maximum h(1/2) = log 2, limits 0 at both ends; each term of a(1-b) + b(1-a) + ab(a+b-ab) is >= 0 on [0,1]^2 so 0 <= pi <= 1",
    }


def free_symbols() -> dict[str, sp.Symbol]:
    names = "q beta M0 M1 H0 H1 E00 E01 E11 K00 K01 K11 P01 F0 F1"
    return dict(zip(names.split(), sp.symbols(names)))


def claim_L1_identity(mutate: str | None = None) -> dict[str, object]:
    S = free_symbols()
    q, beta = S["q"], S["beta"]
    M = (1 - q) * S["M0"] + q * S["M1"]
    Pxy = (1 - q) ** 2 * S["E00"] + 2 * q * (1 - q) * S["E01"] + q ** 2 * S["E11"]
    Ph = (1 - q) * S["H0"] + q * S["H1"]
    if mutate == "mixture_protocol_term":
        ehpi = (1 - q) ** 2 * S["K00"] + 2 * q * (1 - q) * S["K01"] + q ** 2 * S["K11"]
    else:
        ehpi = (1 - q) * S["K00"] + q * S["K11"]
    gap = (1 - beta) * Pxy + beta * ehpi - Ph
    I00 = M * ((1 - beta) * S["E00"] + beta * S["K00"]) - S["M0"] * S["H0"]
    I11 = M * ((1 - beta) * S["E11"] + beta * S["K11"]) - S["M1"] * S["H1"]
    I01 = M * ((1 - beta) * S["E01"] + beta * S["K01"]) - (S["M1"] * S["H0"] + S["M0"] * S["H1"]) / 2
    c = beta * (S["K00"] + S["K11"] - 2 * S["K01"])
    T = 2 * I01 / M + c
    rhs = ((1 - q) ** 2 * I00 + q ** 2 * I11) / M + q * (1 - q) * T
    residual = sp.simplify(sp.cancel(gap - rhs))
    return {
        "name": "L1_free_algebra_chain_identity" if mutate is None else f"L1_mutated_{mutate}",
        "label": "MACHINE-VERIFIED" if residual == 0 else "FAILED",
        "statement": "gap = ((1-q)^2 I00 + q^2 I11)/M + q(1-q)(2 I01/M + c) with M = (1-q)M0 + qM1, I_ab = <.,D_M>, c = beta<(mu-nu)^(x)2, K>, as an identity in the pairings E_ab, K_ab, H_a, M_a",
        "residual": str(residual),
        "pairing_facts_used": [
            "<P(x)P, k> = (1-q)^2 <mu(x)mu,k> + 2q(1-q) <mu(x)nu,k> + q^2 <nu(x)nu,k> for symmetric k (bilinearity + symmetry)",
            "<mu(x)nu, t h(s)> = <mu,h><nu,x> and <mu(x)nu, s h(t)> = <mu,x><nu,h> (product rule / Fubini, bounded continuous factors)",
            "<mu(x)mu, (t h(s) + s h(t))/2> = <mu,x><mu,h>",
        ],
    }


def claim_L2_monotone() -> dict[str, object]:
    s, t, hs, ht, hst, Kst, M, m, beta = sp.symbols("s t h_s h_t h_st K_st M m beta")
    D = lambda mean: mean * ((1 - beta) * hst + beta * Kst) - (t * hs + s * ht) / 2  # noqa: E731
    residual = sp.simplify(D(M) / M - D(m) / m - (t * hs + s * ht) * (M - m) / (2 * M * m))
    return {
        "name": "L2_mean_monotonicity",
        "label": "PROVED" if residual == 0 else "FAILED",
        "statement": "D_M/M - D_m/m = (t h(s) + s h(t))(M - m)/(2 M m) >= 0 pointwise on [0,1]^2 whenever M >= m > 0 (h >= 0, s,t >= 0)",
        "residual": str(residual),
    }


def claim_L3_splitting(mutate: str | None = None) -> dict[str, object]:
    S = free_symbols()
    beta = S["beta"]
    Q00, Q01, Q11 = beta * S["K00"], beta * S["K01"], beta * S["K11"]
    c = Q00 + Q11 - 2 * Q01
    T_m = c + 2 * (S["P01"] + Q01)            # c + 2 <mu(x)nu, R>, R = P2 + Q2
    F0, F1 = S["F0"], S["F1"]
    rhs = (Q00 - F0 ** 2) + (Q11 - F1 ** 2) + 2 * (S["P01"] + F0 * F1)
    if mutate != "drop_square_term":
        rhs += (F0 - F1) ** 2
    residual = sp.expand(T_m - rhs)
    return {
        "name": "L3_rank_one_splitting" if mutate is None else f"L3_mutated_{mutate}",
        "label": "MACHINE-VERIFIED" if residual == 0 else "FAILED",
        "statement": "c + 2<mu(x)nu,R> = <mu(x)mu,B> + <nu(x)nu,B> + 2<mu(x)nu,A> + (<mu,phi> - <nu,phi>)^2 with A = P2 + phi(x)phi, B = Q2 - phi(x)phi, <mu(x)nu, phi(x)phi> = <mu,phi><nu,phi>",
        "residual": str(residual),
    }


def claim_L4_B_from_Phi(beta_ball: arb) -> dict[str, object]:
    Kst, Kss, Ktt = sp.symbols("K_st K_ss K_tt", nonnegative=True)
    root = sp.sqrt(Kss * Ktt)
    residual = sp.expand(Kst ** 2 - Kss * Ktt - (Kst - root) * (Kst + root))
    beta_positive = bool(beta_ball > 0)
    ok = residual == 0 and beta_positive
    return {
        "name": "L4_B_nonnegative_from_Phi",
        "label": "PROVED" if ok else "FAILED",
        "statement": "Phi = (K_st - sqrt(K_ss K_tt))(K_st + sqrt(K_ss K_tt)); B = beta (K_st - sqrt(K_ss K_tt)); K >= 0 and beta* > 0, so Phi >= 0 on [0,1]^2 implies B >= 0 on [0,1]^2",
        "factorisation_residual": str(residual),
        "beta_star_positive_certified": beta_positive,
        "degenerate_case": "if K_st + sqrt(K_ss K_tt) = 0 then K_st = 0 and K_ss K_tt = 0, so B = 0",
    }


def coupling_identity(n: int, m: int, mutate: str | None = None) -> dict[str, object]:
    """The nm-atom paired formula with masses a_i b_j equals the (n,m)-atom general functional."""
    q, beta = sp.symbols("q beta")
    a = sp.symbols(f"a1:{n + 1}")
    b = sp.symbols(f"b1:{m + 1}")

    def sym2(name: str, u: tuple[str, int], v: tuple[str, int]) -> sp.Symbol:
        lo, hi = sorted((u, v))
        return sp.Symbol(f"{name}_{lo[0]}{lo[1]}_{hi[0]}{hi[1]}")

    def sym1(name: str, u: tuple[str, int]) -> sp.Symbol:
        return sp.Symbol(f"{name}_{u[0]}{u[1]}")

    def pair2_sym(A: list, B: list, name: str) -> sp.Expr:
        return sum(wa * wb * sym2(name, ua, ub) for wa, ua in A for wb, ub in B)

    def pair1_sym(A: list, name: str) -> sp.Expr:
        return sum(w * sym1(name, u) for w, u in A)

    mu = [(a[i], ("x", i)) for i in range(n)]
    nu = [(b[j], ("y", j)) for j in range(m)]
    P = [((1 - q) * w, u) for w, u in mu] + [(q * w, u) for w, u in nu]
    gap_gen = (1 - beta) * pair2_sym(P, P, "hxy") \
        + beta * ((1 - q) * pair2_sym(mu, mu, "K") + q * pair2_sym(nu, nu, "K")) - pair1_sym(P, "h")
    mean_gen = pair1_sym(P, "x")

    if mutate == "wrong_coupling":
        g = {(i, j): (a[i] + b[j]) / (n + m) for i in range(n) for j in range(m)}   # still sums to 1
    else:
        g = {(i, j): a[i] * b[j] for i in range(n) for j in range(m)}
    atoms = [(g[(i, j)], ("x", i), ("y", j)) for i in range(n) for j in range(m)]
    mix = [((1 - q) * w, ux) for w, ux, _ in atoms] + [(q * w, uy) for w, _, uy in atoms]
    ehxy = pair2_sym(mix, mix, "hxy")
    comp_x = sum(w1 * w2 * sym2("K", u1, u2) for w1, u1, _ in atoms for w2, u2, _ in atoms)
    comp_y = sum(w1 * w2 * sym2("K", v1, v2) for w1, _, v1 in atoms for w2, _, v2 in atoms)
    ehpi = (1 - q) * comp_x + q * comp_y
    ehx = pair1_sym(mix, "h")
    gap_paired = (1 - beta) * ehxy + beta * ehpi - ehx
    mean_paired = pair1_sym(mix, "x")

    subs = {a[n - 1]: 1 - sum(a[:n - 1]), b[m - 1]: 1 - sum(b[:m - 1])}
    res_gap = sp.expand((gap_paired - gap_gen).subs(subs))
    res_mean = sp.expand((mean_paired - mean_gen).subs(subs))
    ok = res_gap == 0 and res_mean == 0
    return {
        "sizes": [n, m],
        "atoms_paired": n * m,
        "ok": ok,
        "gap_residual": str(res_gap) if len(str(res_gap)) < 200 else str(res_gap)[:200] + "...",
        "mean_residual": str(res_mean) if len(str(res_mean)) < 200 else str(res_mean)[:200] + "...",
    }


def claim_a_product_coupling(sizes: Sequence[tuple[int, int]], mutate: str | None = None) -> dict[str, object]:
    results = [coupling_identity(n, m, mutate) for n, m in sizes]
    ok = all(r["ok"] for r in results)
    return {
        "name": "a_product_coupling" if mutate is None else f"a_mutated_{mutate}",
        "label": "PROVED" if ok else "FAILED",
        "statement": "for mu with n atoms (masses a_i) and nu with m atoms (masses b_j), the nm-atom paired configuration g_ij = a_i b_j, x_i, y_j has the same gap and mean; polynomial identity in the masses modulo sum a = sum b = 1, with all function values free symbols",
        "general_argument": "sum_j a_i b_j f(x_i) = a_i f(x_i) sum_j b_j = a_i f(x_i), and likewise for every pairing; the paired formula depends only on the marginals",
        "symbolic_checks": results,
    }


# --------------------------------------------------------------------------- #
# Arb cross-checks at TEST points (parameter-free identities)
# --------------------------------------------------------------------------- #

def random_measure(rng: random.Random, n: int) -> list[tuple[Fraction, Fraction]]:
    weights = [rng.randint(1, 20) for _ in range(n)]
    total = sum(weights)
    atoms = []
    for w in weights:
        roll = rng.random()
        if roll < 0.05:
            x = Fraction(0)
        elif roll < 0.10:
            x = Fraction(1)
        else:
            x = Fraction(rng.randint(1, 999), 1000)
        atoms.append((Fraction(w, total), x))
    return atoms


def evaluate_nine(a: Sequence[Fraction], X: Sequence[Fraction], Y: Sequence[Fraction],
                  q: Fraction, beta: arb) -> arb:
    """gap from the audited evaluator: numerator - ehx at [a1,a2,q,b0,b2,b4,b1,b3,b5]."""
    values = [frac_arb(a[0]), frac_arb(a[1]), frac_arb(q),
              frac_arb(X[0]), frac_arb(X[1]), frac_arb(X[2]),
              frac_arb(Y[0]), frac_arb(Y[1]), frac_arb(Y[2])]
    terms = liu9_objective.evaluate_arb(values, beta)
    return terms.numerator - terms.ehx


def arb_crosschecks(rng: random.Random, beta_test: arb, protocol: Callable[[arb, arb], arb] = pi_,
                    draws: int = 40) -> dict[str, object]:
    """Three agreements, each certified to IDENTITY_TOL at every draw (TEST beta: identities are parameter-free)."""
    worst_paired3 = arb(0)
    worst_couple = arb(0)
    worst_chain = arb(0)
    for _ in range(draws):
        a = [w for w, _ in random_measure(rng, 3)]
        X = [x for _, x in random_measure(rng, 3)]
        Y = [y for _, y in random_measure(rng, 3)]
        q = Fraction(rng.randint(0, 100), 100)
        diff = gap_paired_atoms(a, X, Y, q, beta_test) - evaluate_nine(a, X, Y, q, beta_test)
        worst_paired3 = max(worst_paired3, abs_upper(diff))
        mu3 = random_measure(rng, 3)
        y = random_measure(rng, 1)
        g, Xc, Yc = product_coupling(mu3, y)
        diff = gap_general(mu3, y, q, beta_test, protocol) - evaluate_nine(g, Xc, Yc, q, beta_test)
        worst_couple = max(worst_couple, abs_upper(diff))
        x1 = random_measure(rng, 1)
        nu3 = random_measure(rng, 3)
        g, Xc, Yc = product_coupling(x1, nu3)
        diff = gap_general(x1, nu3, q, beta_test, protocol) - evaluate_nine(g, Xc, Yc, q, beta_test)
        worst_couple = max(worst_couple, abs_upper(diff))
        mu = random_measure(rng, rng.randint(1, 4))
        nu = random_measure(rng, rng.randint(1, 4))
        if mean_exact(mu, nu, q) > 0:
            diff = gap_general(mu, nu, q, beta_test, protocol) - gap_chain_rhs(mu, nu, q, beta_test)
            worst_chain = max(worst_chain, abs_upper(diff))
    tol = frac_arb(IDENTITY_TOL)
    ok = bool(worst_paired3 <= tol) and bool(worst_couple <= tol) and bool(worst_chain <= tol)
    return {
        "name": "arb_crosschecks_against_audited_evaluator",
        "label": "MACHINE-VERIFIED" if ok else "FAILED",
        "draws": draws,
        "parameters": "TEST rationals (the identities checked are parameter-free)",
        "tolerance": str(IDENTITY_TOL),
        "worst_abs_diff_paired3_vs_evaluate_arb": fmt(worst_paired3, 12),
        "worst_abs_diff_product_coupling_31_13_vs_evaluate_arb": fmt(worst_couple, 12),
        "worst_abs_diff_general_gap_vs_chain_rhs": fmt(worst_chain, 12),
        "statement": "the general functional coincides with liu9_objective.evaluate_arb (numerator - ehx) on three-atom paired inputs and on product-coupled (3,1)/(1,3) inputs, and L1 holds in Arb at unequal atom counts",
    }


# --------------------------------------------------------------------------- #
# Sign claims at the EXACT balls
# --------------------------------------------------------------------------- #

def sample_nonnegativity(rng: random.Random, C: Constants, draws: int = 300) -> dict[str, object]:
    feasible = certified_nonneg = undecided = negatives = 0
    min_lower = None
    min_at = None
    m_hi = C.m.upper()
    for _ in range(draws):
        mu = random_measure(rng, rng.randint(1, 4))
        nu = random_measure(rng, rng.randint(1, 4))
        q = Fraction(rng.randint(0, 100), 100)
        if not bool(frac_arb(mean_exact(mu, nu, q)) > m_hi):
            continue                         # only certified mean-feasible draws
        feasible += 1
        value = gap_general(mu, nu, q, C.beta)
        lower = value.lower()
        if min_lower is None or lower < min_lower:
            min_lower = lower
            min_at = {"mu": [[str(w), str(x)] for w, x in mu], "nu": [[str(w), str(y)] for w, y in nu], "q": str(q)}
        if certified_above(value, arb(0)):
            certified_nonneg += 1
        elif certified_below(value, arb(0)):
            negatives += 1
        else:
            undecided += 1
    return {
        "name": "random_sample_mean_feasible_nonnegativity",
        "label": "COMPUTATIONAL-EVIDENCE",
        "parameters": "EXACT balls",
        "draws": draws,
        "certified_mean_feasible": feasible,
        "certified_positive": certified_nonneg,
        "sign_undecided": undecided,
        "certified_negative": negatives,
        "min_certified_lower_bound": fmt(min_lower, 20) if min_lower is not None else None,
        "min_at": min_at,
        "statement": "sanity sample only; the theorem is proved by composition, not by this sample",
    }


def witnesses(C: Constants) -> list[dict[str, object]]:
    out = []
    half = frac_arb(Fraction(1, 2))
    # W1a: Liu's minimiser at the EXACT balls -- enclosure must contain 0; equality is PROVED by algebra
    p = C.m / C.x
    minimiser = [(1 - p, arb(0)), (p, C.x)]
    value = gap_general_arb(minimiser, minimiser, half, C.beta)
    contains_zero = bool(value.contains(arb(0)))
    out.append({
        "name": "W1a_liu_minimiser_exact_zero",
        "label": "PROVED" if contains_zero else "FAILED",
        "parameters": "EXACT balls",
        "mu_eq_nu": [["1 - m*/x*", "0"], ["m*/x*", "x*"]],
        "gap_enclosure": fmt(value, 20),
        "enclosure_contains_zero": contains_zero,
        "proof_of_equality": "mean = m* exactly; gap = I00/m* with I00 = p^2 D_{m*}(x*,x*) + 2p(1-p) D_{m*}(0,x*) + (1-p)^2 D_{m*}(0,0); D_{m*}(0,.) = 0 since h(0) = 0 and pi(0,.) = 0; D_{m*}(x*,x*) = m* h(x*^2) - x* h(x*) = 0 by h(pi(x*,x*)) = h(1 - x*^2) = h(x*^2) (Liu (87)) and m* = x* h(x*)/h(x*^2) (Liu (88))",
        "statement": "equality is a theorem, not a numerical fact; the enclosure is the consistency check",
    })
    # W1b: the same configuration with the TEST rationals declared exact -- certified NEGATIVE
    p_t = C.m_t / C.x_t
    min_t = [(1 - p_t, Fraction(0)), (p_t, C.x_t)]
    value_t = gap_general(min_t, min_t, Fraction(1, 2), frac_arb(C.beta_t))
    neg = certified_below(value_t, arb(0))
    out.append({
        "name": "W1b_truncated_constants_certified_negative",
        "label": "MACHINE-VERIFIED" if neg else "FAILED",
        "parameters": "TEST rationals declared exact (beta_t, m_t, x_t)",
        "mean_equals_m_t_exactly": mean_exact(min_t, min_t, Fraction(1, 2)) == C.m_t,
        "gap_enclosure": fmt(value_t, 20),
        "certified_negative": neg,
        "statement": "with the 100-digit truncations in place of the constants the theorem is FALSE (by ~6e-102): the interior zero of R at (x*,x*) is nondegenerate and exact, so perturbing (beta*, m*, x*) cannot be assumed harmless (this joint truncation lowers R near (x*,x*); raising m alone would keep R >= 0 by L2). Recorded from the independent review's counterexample. This is why every sign claim here uses the exact balls.",
    })
    # W2: unconstrained statement is false (paired, mean 1/2 < m*)
    d_half = [(Fraction(1), Fraction(1, 2))]
    value = gap_general(d_half, d_half, Fraction(1, 2), C.beta)
    out.append({
        "name": "W2_mean_constraint_load_bearing_paired",
        "label": "MACHINE-VERIFIED" if certified_below(value, arb(0)) else "FAILED",
        "parameters": "EXACT balls",
        "mu_eq_nu": [["1", "1/2"]],
        "mean": "1/2",
        "mean_below_m_certified": bool(half < C.m),
        "gap": fmt(value, 20),
        "statement": "gap < 0 certified: the theorem is false without M >= m*",
    })
    # W3: non-paired pair with mean below m*
    mu = [(Fraction(1, 3), Fraction(1, 5)), (Fraction(2, 3), Fraction(3, 5))]
    nu = [(Fraction(1), Fraction(1, 2))]
    value = gap_general(mu, nu, Fraction(1, 3), C.beta)
    mean3 = mean_exact(mu, nu, Fraction(1, 3))
    out.append({
        "name": "W3_mean_constraint_load_bearing_general",
        "label": "MACHINE-VERIFIED" if certified_below(value, arb(0)) else "FAILED",
        "parameters": "EXACT balls",
        "mu": [["1/3", "1/5"], ["2/3", "3/5"]],
        "nu": [["1", "1/2"]],
        "q": "1/3",
        "mean": str(mean3),
        "mean_below_m_certified": bool(frac_arb(mean3) < C.m),
        "gap": fmt(value, 20),
        "statement": "gap < 0 certified at a two-atom / one-atom pair with mean below m*",
    })
    return out


# --------------------------------------------------------------------------- #
# Mutations
# --------------------------------------------------------------------------- #

def run_mutations(rng: random.Random, C: Constants, wits: list[dict[str, object]]) -> list[dict[str, object]]:
    muts: list[dict[str, object]] = []

    r = claim_L1_identity(mutate="mixture_protocol_term")
    muts.append({"mutation": "L1_mixture_form_protocol_term", "caught": r["label"] == "FAILED",
                 "observed": "residual " + r["residual"]})

    r = claim_L3_splitting(mutate="drop_square_term")
    muts.append({"mutation": "L3_drop_square_term", "caught": r["label"] == "FAILED",
                 "observed": "residual " + r["residual"]})

    r = claim_a_product_coupling([(2, 3)], mutate="wrong_coupling")
    muts.append({"mutation": "a_wrong_coupling_masses_(a_i+b_j)/(n+m)", "caught": r["label"] == "FAILED",
                 "observed": "gap residual " + str(r["symbolic_checks"][0]["gap_residual"])[:120]})

    d_half = [(Fraction(1), Fraction(1, 2))]
    value = gap_general(d_half, d_half, Fraction(1, 2), C.beta)
    muts.append({"mutation": "drop_mean_constraint_(claim gap>=0 for all mu,nu,q)", "caught": certified_below(value, arb(0)),
                 "observed": f"gap(delta_1/2, delta_1/2) = {fmt(value, 15)} certified < 0 at the exact balls, mean 1/2 < m*"})

    caught = False
    worst = arb(0)
    beta_t = frac_arb(C.beta_t)
    for _ in range(5):
        mu3 = random_measure(rng, 3)
        y = random_measure(rng, 1)
        q = Fraction(rng.randint(0, 100), 100)
        g, Xc, Yc = product_coupling(mu3, y)
        diff = gap_general(mu3, y, q, beta_t, pi_wrong) - evaluate_nine(g, Xc, Yc, q, beta_t)
        if abs(diff).lower() > frac_arb(Fraction(1, 10 ** 6)):
            caught = True
        worst = max(worst, abs_upper(diff))
    muts.append({"mutation": "wrong_protocol_s+t-st_vs_evaluate_arb", "caught": caught,
                 "observed": f"largest |diff| {fmt(worst, 10)} (certified > 1e-6 at some draw)"})

    s = t = frac_arb(Fraction(1, 2))
    gapM = D_M(s, t, arb(1), C.beta) / arb(1) - R_kernel(s, t, C.beta, C.m)
    muts.append({"mutation": "L2_reversed_monotonicity_(claim D_M/M <= R for M > m)", "caught": certified_above(gapM, arb(0)),
                 "observed": f"D_1/1 - R at (1/2,1/2) = {fmt(gapM, 15)} certified > 0 at the exact balls"})

    w1b = next(w for w in wits if w["name"] == "W1b_truncated_constants_certified_negative")
    muts.append({"mutation": "declare_the_100_digit_truncations_to_be_the_constants", "caught": bool(w1b["certified_negative"]),
                 "observed": "W1b: gap at Liu's minimiser with (beta_t, m_t, x_t) = " + str(w1b["gap_enclosure"]) + " certified < 0"})
    return muts


# --------------------------------------------------------------------------- #
# Dependencies and report
# --------------------------------------------------------------------------- #

def check_dependencies() -> tuple[bool, dict[str, object], dict[str, object]]:
    out: dict[str, object] = {}
    ok = True
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
                     "claim_status_expected": spec["claim_status"], "claim_status_on_disk": status,
                     "match": match, "provides": spec["provides"]}
    return ok, out, twovar_constants


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def build_report() -> dict[str, object]:
    deps_ok, deps, twovar_constants = check_dependencies()
    C = Constants()
    provenance = C.provenance(twovar_constants)
    constants_ok = (C.bracket_ok and C.derivative_positive_on_bracket and C.binding_overlap
                    and all(provenance["twovar_uses_the_same_definitions"].values()))
    beta_test = frac_arb(C.beta_t)
    rng = random.Random(SEED)

    sizes = [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 2), (3, 3)]
    claims = [
        claim_L0_regularity(),
        claim_L1_identity(),
        claim_L2_monotone(),
        claim_L3_splitting(),
        claim_L4_B_from_Phi(C.beta),
        claim_a_product_coupling(sizes),
        arb_crosschecks(rng, beta_test),
    ]
    sample = sample_nonnegativity(rng, C)
    wits = witnesses(C)
    mutations = run_mutations(rng, C, wits)

    failed = [c["name"] for c in claims if c["label"] == "FAILED"] + [w["name"] for w in wits if w["label"] == "FAILED"]
    all_caught = all(mt["caught"] for mt in mutations)
    if not deps_ok:
        status = "DEPENDENCY-MISMATCH"
    elif not constants_ok:
        status = "CONSTANTS-MISMATCH"
    elif failed or not all_caught or sample["certified_negative"] > 0:
        status = "FAILED"
    else:
        status = "PROVED"

    composition = {
        "name": "b_general_measures",
        "label": status if status in ("PROVED", "FAILED") else "OPEN",
        "statement": "at the exact constants (beta*, m*): for all Borel probability measures mu, nu on [0,1] and q in [0,1] with (1-q)<mu,x> + q<nu,x> >= m*: gap(mu,nu,q) >= 0",
        "proof": [
            "L1: gap = (1-q)^2 <mu(x)mu, D_M/M> + q^2 <nu(x)nu, D_M/M> + q(1-q)[2 <mu(x)nu, D_M/M> + c]",
            "L2 with M >= m*: each D_M/M pairing >= the same pairing of R (integrand difference >= 0, product probability measures); (1-q)^2, q^2, q(1-q) >= 0",
            "R >= 0 pointwise at (beta*, m*) (liu9-h2-twovar.json): the first two terms are >= 0",
            "L3: 2 <mu(x)nu, R> + c = <mu(x)mu,B> + <nu(x)nu,B> + 2<mu(x)nu,A> + (<mu,phi>-<nu,phi>)^2",
            "A >= 0 pointwise at (beta*, m*) (liu9-h2-twovar.json); B >= 0 pointwise by L4 from Phi >= 0 (liu9-h2-boundary.json) and beta* > 0; a square is >= 0",
            "L0: every kernel is continuous on [0,1]^2, so every pairing exists, is finite, and the pairing of a nonnegative kernel against a product probability measure is nonnegative",
        ],
        "constants": "the theorem is about the root-defined (beta*, m*) only; W1b shows it fails for the 100-digit truncations, so perturbing the constants cannot be assumed harmless; they are certified as balls, never as decimals",
        "finite_support_special_case": "mu = sum a_i delta_x_i, nu = sum b_j delta_y_j with any atom counts and masses (claim (a) shows this is also the paired class via the product coupling)",
        "liu_nine_parameter_case": "shared masses a_i, three atoms each: liu9_objective._formula; the audited evaluator agrees with the general functional there (arb_crosschecks)",
        "weak_density_remark": {
            "label": "PROVED (elementary), not load-bearing",
            "statement": "finitely supported laws are weakly dense in the probability measures on [0,1], and gap(mu,nu,q) is jointly weakly continuous in (mu,nu) because every kernel is continuous on the compact square; hence the general theorem is also the weak closure of the finitely supported one. The direct proof above never takes a limit, so no modulus of continuity of h (which is not Lipschitz at 0, 1) is needed.",
        },
    }

    report: dict[str, object] = {
        "artifact": "liu9-h2-general-lift",
        "claim": composition["statement"],
        "claim_status": status,
        "definitions": {
            "h": "-u log u - (1-u) log(1-u), natural log, h(0) = h(1) = 0",
            "pi": "s t (1 + (1-s)(1-t))",
            "K": "h(pi(s,t))",
            "D_M": "M[(1-beta) h(st) + beta h(pi)] - (t h(s) + s h(t))/2",
            "R": "D_m/m = P2 + Q2, P2 = (1-beta) h(st) - (t h(s) + s h(t))/(2m), Q2 = beta h(pi)",
            "phi": "sqrt(Q2(s,s))",
            "A": "P2 + phi (x) phi",
            "B": "Q2 - phi (x) phi",
            "Phi": "h(pi(s,t))^2 - h(pi(s,s)) h(pi(t,t))",
            "gap": "(1-beta)<P(x)P,h(xy)> + beta[(1-q)<mu(x)mu,K> + q<nu(x)nu,K>] - <P,h>, P = (1-q)mu + q nu",
            "mean_constraint": "M = <P,x> >= m*",
        },
        "constants": provenance,
        "claims": claims,
        "composition": composition,
        "random_sample": sample,
        "witnesses": wits,
        "mutations": mutations,
        "all_mutations_caught": all_caught,
        "failed_claims": failed,
        "dependencies": deps,
        "labels": {
            "PROVED": "complete argument, elementary steps machine-checked where stated",
            "MACHINE-VERIFIED": "exact symbolic residual 0, or Arb-certified comparison at the listed points",
            "COMPUTATIONAL-EVIDENCE": "random sample, not a proof",
            "OPEN": "not established",
        },
        "precision_bits": PRECISION_BITS,
        "natural_log": True,
        "seed": SEED,
        "single_core_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "FLINT_NUM_THREADS")},
        "sympy_version": sp.__version__,
        "tool": "liu9_h2_general_lift.py (sympy exact identities + python-flint Arb)",
        "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "determinism": "no timestamps; fixed seed; sympy residuals as strings; Arb strings at fixed digit counts; byte-identical under the recorded single-thread environment (thread variables are set with setdefault, so a pre-set value is inherited and recorded)",
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON plus newline with report_sha256 omitted",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def write_report(report: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=1) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Liu H2: lift from the paired class to general measures.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report()
    write_report(report, args.output)
    print(f"GENERAL_LIFT {report['claim_status']}")
    for claim in report["claims"]:
        print(f"  {claim['name']:<44} {claim['label']}")
    print(f"  {'b_general_measures':<44} {report['composition']['label']}")
    for w in report["witnesses"]:
        print(f"  {w['name']:<44} {w['label']}  {w.get('gap_enclosure', w.get('gap'))}")
    rs = report["random_sample"]
    print(f"  sample: {rs['certified_mean_feasible']} certified mean-feasible of {rs['draws']}, certified positive {rs['certified_positive']}, "
          f"undecided {rs['sign_undecided']}, negative {rs['certified_negative']}, min lower bound {rs['min_certified_lower_bound']}")
    for mt in report["mutations"]:
        print(f"  mutation {mt['mutation']:<60} caught={mt['caught']}")
    cs = report["constants"]
    print(f"  constants: bracket ok {cs['exact_bracket_of_xstar']['f_lo_negative_and_f_hi_positive_certified']}, twovar same definitions {cs['twovar_uses_the_same_definitions']}")
    print(f"REPORT {args.output}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    return 0 if report["claim_status"] == "PROVED" else 3


if __name__ == "__main__":
    raise SystemExit(main())
