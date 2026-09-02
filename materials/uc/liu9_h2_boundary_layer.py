#!/usr/bin/env python3
"""Rigorous Arb certificate: Phi(s,t) >= 0 on the WHOLE unit square.

    Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s)) * h(pi(t,t))

with h(u) = -u log u - (1-u) log(1-u) (natural log), pi(s,t) = s t (1+(1-s)(1-t))
(Liu's protocol).  Phi >= 0 is equivalent to B >= 0 in the H2 reduction chain
(B = beta (h(pi_st) - sqrt(h(pi_ss) h(pi_tt))), beta > 0), and the residual
endpoint layer left open by liu9-h2-phi.json is a subset of the square, so this
artifact discharges that layer by inclusion.  Phi contains no parameter: beta, m,
xstar never enter (beta is only the positive prefactor of B).

The proof does NOT subdivide near the singular endpoints.  It rests on an exact
identity that turns the two-variable inequality into a statement about ONE
one-variable function, gamma(x) = log( h(pi(x,x)) / pi(x,x) ), whose sufficient
condition kappa >= 2 has a uniform margin (numerical minimum 2.394 at s=t=0.54)
and diverges to +infinity at every boundary.  The endpoint layer is therefore the
EASY part; the only two analytic estimates needed are lower bounds on |gamma'|
near 0 and near 1, both elementary.

Notation.  Lam(u) = h(u)/u = log(1/u) + mu(u) with mu(u) = -(1-u) log(1-u)/u.
pi_d(x) = pi(x,x) = x^2 (2 - 2x + x^2).  g = Lam o pi_d,  gamma = log g.
Subscripts: pi_st = pi(s,t), Lam_ss = Lam(pi(s,s)), etc.  a = 1-s, b = 1-t.

L0 (elementary; proofs in the docstrings of the checking functions)
 (a) u h'(u) - h(u) = log(1-u), hence Lam'(u) = log(1-u)/u^2 < 0 on (0,1).
 (b) pi_d'(x) = 2x(2-3x+2x^2) > 0 on (0,1] (discriminant -7), and
     pi_d''(x) = 4(1-3x+3x^2) > 0 (discriminant -3): pi_d and pi_d' increase.
 (c) g and gamma are strictly decreasing on (0,1), and
     gamma'(x) = log(1-pi_d) pi_d' / (pi_d h(pi_d)).
 (d) mu(u) = 1 - sum_{k>=2} u^{k-1}/(k(k-1)) on [0,1): decreasing, concave,
     0 <= mu <= 1, continuous on [0,1] with mu(1)=0.
 (e) pi_ss pi_tt - pi_st^2 = s^2 t^2 (s-t)^2   (polynomial identity, checked
     EXACTLY on a 6x6 rational grid; degree <= 4 in each variable).
 (f) M_mu := mu(pi_ss) + mu(pi_tt) - 2 mu(pi_st) <= 0.  By (e),
     pi_st <= sqrt(pi_ss pi_tt) <= (pi_ss+pi_tt)/2; mu decreasing then concave.
 (g) Lam >= 0 on (0,1].

L1 (exact identity).  With delta = (s-t)^2/(1+ab)^2 (so pi_ss pi_tt =
    pi_st^2 (1+delta)) and Lbar = (Lam_ss+Lam_tt)/2,

    Phi = pi_st^2 [ (Lam_ss-Lam_tt)^2/4 + (log(1+delta) - M_mu)/2 (Lam_st+Lbar) ]
          - s^2 t^2 (s-t)^2 Lam_ss Lam_tt.

    Proof.  Phi = pi_st^2 (Lam_st^2 - Lam_ss Lam_tt) - (pi_ss pi_tt - pi_st^2)
    Lam_ss Lam_tt;  Lam_st^2 - Lam_ss Lam_tt = (Lam_ss-Lam_tt)^2/4 +
    (Lam_st-Lbar)(Lam_st+Lbar);  2(Lam_st-Lbar) = log(pi_ss pi_tt/pi_st^2) - M_mu
    because Lam = log(1/.) + mu;  then (e).
    The middle term is >= 0 by (f), (g) and log(1+delta) >= 0, so with
    pi_st = s t (1+ab):

    (dagger)  Phi >= s^2 t^2 [ (1+ab)^2 (g(s)-g(t))^2/4 - (s-t)^2 g(s) g(t) ].

L2 (master condition).  For 0 < s < t < 1 put
      kappa(s,t) = (1 + (1-s)(1-t)) (gamma(s)-gamma(t)) / (t-s)  (> 0 by (c)).
    Then  Phi(s,t) >= s^2 t^2 (s-t)^2 g(s) g(t) (kappa^2/4 - 1);  in particular
    kappa >= 2 implies Phi >= 0.
    Proof.  With x = (gamma(s)-gamma(t))/2 > 0,  g(s)-g(t) = 2 sqrt(g(s)g(t))
    sinh x, so the bracket in (dagger) is (s-t)^2 g(s)g(t) [((1+ab) sinh x/(t-s))^2
    - 1], and sinh x >= x gives (1+ab) sinh x/(t-s) >= kappa/2.

L3 (boundary and symmetry).  Phi(s,s) = 0, Phi(0,t) = 0, Phi(s,1) = h(s)^2 >= 0
    (pi(s,1) = s, pi(1,1) = 1), Phi(s,t) = Phi(t,s).  Hence it suffices to
    certify kappa >= 2 on the open triangle 0 < s < t < 1.

L4 (endpoint derivative bounds).
 N0: for 0 < x <= x0 < exp(-(1+log 2)/2) = 0.4289...:
       |gamma'(x)| >= nu0(x0) := (2-3x0) / ( x0 (2 log(1/x0) + 1 - log 2) ).
     Proof.  -log(1-pi) >= pi;  pi_d'(x) >= 2x(2-3x0);  h(pi) <= pi(log(1/pi)+1)
     (mu <= 1);  pi_d(x) <= 2x^2 and u(log(1/u)+1) is increasing on (0,1), so
     h(pi) <= 2x^2 (2 log(1/x) - log 2 + 1);  finally x (2 log(1/x) + 1 - log 2)
     is increasing for x < exp(-(1+log 2)/2), so it is at most its value at x0.
 N1: for 3/4 <= x1 <= x < 1, with ell1 = log(1/(1-pi_d(x1))):
       |gamma'(x)| >= nu1(x1) := pi_d'(x1)/(1-pi_d(x1)) * ell1/(ell1+1).
     Proof.  pi <= 1;  h(pi) = h(1-pi) <= (1-pi)(log(1/(1-pi)) + 1);  on
     [3/4,1] the factor 2-3x+2x^2 increases (vertex 3/4) so pi_d' >= pi_d'(x1);
     ell/(ell+1) increases in ell and 1/(1-pi) increases in x.

L5 (interior derivative bound).  N(x) := -log(1-pi_d(x)) pi_d'(x) is positive and
    strictly increasing on (0,1) (product of positive increasing functions), so
    for [l,r] in (0,1):  inf_{[l,r]} |gamma'| >= N(l) / ( pi_d(r) max_{[pi_d(l),
    pi_d(r)]} h ),  where max h over an interval is log 2 if it contains 1/2 and
    the larger endpoint value otherwise (h increases on [0,1/2], decreases after).

Cell rules.  For a rational cell C = [s1,s2] x [t1,t2] with s1 < t2, every
(s,t) in C with 0 < s < t < 1 satisfies, by the mean value theorem and (c),
  (D)  kappa >= (1+(1-s2)(1-t2)) * inf_{xi in [s1,t2] cap (0,1)} |gamma'(xi)|,
       the infimum being the minimum of N0(min(t2,x0)) on (0,x0], N1(max(s1,x1))
       on [x1,1) and L5 on the middle piece, whichever apply;
  (Q)  if s2 < t1:  kappa >= (1+(1-s2)(1-t2)) (gamma(s2)-gamma(t1)) / (t2-s1).
A cell is accepted when the lower endpoint of either Arb bound exceeds THETA;
otherwise it is bisected along its longer side.  Cells with s1 >= t2 contain no
point with s < t and are discarded.  THETA = 9/4 > 2 is certified, which by L2
gives the quantitative bound Phi >= (17/64) s^2 t^2 (s-t)^2 g(s) g(t).

Invocation:  nice -n 19 ./.venv/bin/python -I -B uc/liu9_h2_boundary_layer.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "FLINT_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_thread_variable, "1")

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from flint import arb, ctx  # noqa: E402

PRECISION_BITS = 400
ctx.prec = PRECISION_BITS

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
LOG2 = TWO.log()

HALF = Fraction(1, 2)
X0 = Fraction(1, 8)          # N0 is applied on (0, min(t2, X0)]
X1 = Fraction(7, 8)          # N1 is applied on [max(s1, X1), 1)
THETA = Fraction(9, 4)       # certified lower bound for kappa (> 2 needed)
MAX_DEPTH = 60
MAX_CELLS = 400_000
RESIDUAL_TOLERANCE = arb("1e-100")
LAYER_A = Fraction(1, 16384)  # residual layer half-width of liu9-h2-phi.json


@dataclass(frozen=True)
class Cell:
    s_lo: Fraction
    s_hi: Fraction
    t_lo: Fraction
    t_hi: Fraction
    depth: int = 0


@dataclass
class CoverResult:
    ok: bool
    processed: int = 0
    accepted_d: int = 0
    accepted_q: int = 0
    discarded: int = 0
    max_depth: int = 0
    min_certified: arb | None = None
    min_cell: Cell | None = None
    layer_cells: int = 0
    failure: str | None = None
    failing_cell: Cell | None = None


# --------------------------------------------------------------------------- #
# Basic arithmetic helpers
# --------------------------------------------------------------------------- #

def frac_arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def format_arb(value: arb, digits: int = 40) -> str:
    return value.str(digits)


def format_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_min(values: list[arb]) -> arb:
    """Minimum of exact (zero-radius) Arb numbers."""
    result = values[0]
    for value in values[1:]:
        if not (value.rad() == 0 and result.rad() == 0):
            raise AssertionError("exact_min requires exact endpoints")
        if value < result:
            result = value
    return result


def exact_max(values: list[arb]) -> arb:
    result = values[0]
    for value in values[1:]:
        if not (value.rad() == 0 and result.rad() == 0):
            raise AssertionError("exact_max requires exact endpoints")
        if value > result:
            result = value
    return result


def h_arb(u: arb) -> arb:
    """Natural-log entropy for a ball strictly inside (0,1)."""
    return -(u * u.log() + (ONE - u) * (ONE - u).log())


def h_frac(u: Fraction) -> arb:
    if u == 0 or u == 1:
        return ZERO
    return h_arb(frac_arb(u))


def mu_arb(u: arb) -> arb:
    return -(ONE - u) * (ONE - u).log() / u


def lam_arb(u: arb) -> arb:
    return h_arb(u) / u


def pi_frac(s: Fraction, t: Fraction) -> Fraction:
    return s * t * (1 + (1 - s) * (1 - t))


def pid_frac(x: Fraction) -> Fraction:
    return x * x * (2 - 2 * x + x * x)


def pidp_frac(x: Fraction) -> Fraction:
    return 2 * x * (2 - 3 * x + 2 * x * x)


def gamma_point(x: Fraction) -> arb:
    """gamma(x) = log(h(pi_d(x))/pi_d(x)) at a rational 0 < x < 1."""
    if not (0 < x < 1):
        raise AssertionError("gamma is evaluated only at interior points")
    p = frac_arb(pid_frac(x))
    return lam_arb(p).log()


def dgamma_point(x: Fraction) -> arb:
    """Closed form of gamma'(x) (L0(c)); used only for sanity samples."""
    p = frac_arb(pid_frac(x))
    return (ONE - p).log() * frac_arb(pidp_frac(x)) / (p * h_arb(p))


def phi_direct(s: Fraction, t: Fraction) -> arb:
    return h_frac(pi_frac(s, t)) ** 2 - h_frac(pi_frac(s, s)) * h_frac(pi_frac(t, t))


# --------------------------------------------------------------------------- #
# L0 / L1 machine checks
# --------------------------------------------------------------------------- #

IDENTITY_POINTS: list[tuple[Fraction, Fraction]] = [
    (Fraction(1, 3), Fraction(2, 3)),
    (Fraction(1, 7), Fraction(5, 7)),
    (Fraction(1, 2), Fraction(3, 5)),
    (Fraction(27, 50), Fraction(11, 20)),
    (Fraction(1, 10), Fraction(9, 10)),
    (Fraction(1, 2 ** 20), Fraction(1, 3)),
    (Fraction(1, 2 ** 20), Fraction(1, 2 ** 19)),
    (Fraction(1, 2 ** 40), Fraction(1, 2 ** 15)),
    (Fraction(1, 3), Fraction(1, 1) - Fraction(1, 2 ** 20)),
    (Fraction(1, 1) - Fraction(1, 2 ** 20), Fraction(1, 1) - Fraction(1, 2 ** 21)),
    (Fraction(1, 2 ** 20), Fraction(1, 1) - Fraction(1, 2 ** 20)),
    (Fraction(1, 16384) - Fraction(1, 2 ** 30), Fraction(17, 32)),
    (Fraction(3, 16384), Fraction(16383, 16384) + Fraction(1, 2 ** 30)),
    (Fraction(9, 10), Fraction(1, 5)),
]


def identity_residual(s: Fraction, t: Fraction, mutate: str | None = None) -> arb:
    """Phi minus the right-hand side of L1 (must enclose zero)."""
    pst, pss, ptt = pi_frac(s, t), pi_frac(s, s), pi_frac(t, t)
    if mutate == "wrong_protocol":
        pst = s + t - s * t
    ast, ass, att = frac_arb(pst), frac_arb(pss), frac_arb(ptt)
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    delta = frac_arb((s - t) ** 2 / (1 + (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    quarter = 4 if mutate != "wrong_quarter" else 2
    bracket = (lss - ltt) ** 2 / quarter + (delta.log1p() - m_mu) / 2 * (lst + (lss + ltt) / 2)
    rhs = ast ** 2 * bracket - frac_arb(s * s * t * t * (s - t) ** 2) * lss * ltt
    return phi_direct(s, t) - rhs


def middle_term(s: Fraction, t: Fraction) -> tuple[arb, arb]:
    """The dropped middle term of L1 and M_mu, at a rational point."""
    pst, pss, ptt = pi_frac(s, t), pi_frac(s, s), pi_frac(t, t)
    ast, ass, att = frac_arb(pst), frac_arb(pss), frac_arb(ptt)
    lst, lss, ltt = lam_arb(ast), lam_arb(ass), lam_arb(att)
    delta = frac_arb((s - t) ** 2 / (1 + (1 - s) * (1 - t)) ** 2)
    m_mu = mu_arb(ass) + mu_arb(att) - 2 * mu_arb(ast)
    return (delta.log1p() - m_mu) / 2 * (lst + (lss + ltt) / 2), m_mu


def polynomial_identity_exact() -> dict[str, object]:
    """L0(e) checked exactly: a polynomial of degree <= 4 in s and <= 4 in t that
    vanishes on a 6x6 grid of distinct rationals is identically zero."""
    grid = [Fraction(k, 7) for k in range(1, 7)]
    worst = Fraction(0)
    for s in grid:
        for t in grid:
            lhs = pi_frac(s, s) * pi_frac(t, t) - pi_frac(s, t) ** 2
            rhs = s * s * t * t * (s - t) ** 2
            worst = max(worst, abs(lhs - rhs))
    return {
        "statement": "pi(s,s)*pi(t,t) - pi(s,t)^2 == s^2*t^2*(s-t)^2",
        "grid": "s,t in {1/7,...,6/7}",
        "degree_bound": "<= 4 in each variable, so a 5x5 grid already suffices",
        "max_abs_residual_exact": format_fraction(worst),
        "status": "PROVED" if worst == 0 else "FAILED",
    }


def derivative_identity_check() -> dict[str, object]:
    """L0(a): u h'(u) - h(u) == log(1-u), checked in Arb at rational points."""
    worst = ZERO
    for u in (Fraction(1, 5), Fraction(1, 2), Fraction(7, 9), Fraction(1, 2 ** 30), 1 - Fraction(1, 2 ** 30)):
        au = frac_arb(u)
        hp = ((ONE - au) / au).log()
        residual = au * hp - h_arb(au) - (ONE - au).log()
        if not residual.contains(ZERO):
            raise AssertionError("u h'(u) - h(u) = log(1-u) failed")
        worst = exact_max([worst, abs(residual).upper()])
    return {
        "statement": "u*h'(u) - h(u) == log(1-u); hence Lam'(u) = log(1-u)/u^2 < 0",
        "algebra": "u*h' = u log(1-u) - u log u; -h = u log u + (1-u) log(1-u); sum = log(1-u)",
        "max_abs_residual_upper": format_arb(worst, 12),
        "status": "PROVED",
    }


def identity_checks(mutate: str | None = None) -> dict[str, object]:
    worst = ZERO
    min_middle: arb | None = None
    max_m_mu: arb | None = None
    for s, t in IDENTITY_POINTS:
        residual = identity_residual(s, t, mutate)
        if not residual.contains(ZERO):
            raise AssertionError(f"L1 identity residual excludes zero at {s},{t}: {residual}")
        if not abs(residual).upper() < RESIDUAL_TOLERANCE:
            raise AssertionError(f"L1 identity residual too wide at {s},{t}: {residual}")
        worst = exact_max([worst, abs(residual).upper()])
        mid, m_mu = middle_term(s, t)
        if not mid.lower() >= 0:
            raise AssertionError(f"middle term not certified nonnegative at {s},{t}")
        if not m_mu.upper() <= 0:
            raise AssertionError(f"M_mu not certified nonpositive at {s},{t}")
        min_middle = mid.lower() if min_middle is None else exact_min([min_middle, mid.lower()])
        max_m_mu = m_mu.upper() if max_m_mu is None else exact_max([max_m_mu, m_mu.upper()])
    return {
        "points": len(IDENTITY_POINTS),
        "max_abs_residual_upper": format_arb(worst, 12),
        "residual_tolerance": "1e-100",
        "middle_term_min_lower_at_points": format_arb(min_middle, 20),
        "M_mu_max_upper_at_points": format_arb(max_m_mu, 20),
        "status": "MACHINE-VERIFIED at the listed rational points; the identity itself is the algebra in the module docstring (L1)",
    }


# --------------------------------------------------------------------------- #
# L4 / L5 derivative lower bounds
# --------------------------------------------------------------------------- #

def n0_validity(x0: Fraction) -> bool:
    """x0 < exp(-(1+log 2)/2), i.e. log(x0) + (1+log 2)/2 < 0, certified in Arb."""
    return bool((frac_arb(x0).log() + (ONE + LOG2) / 2).upper() < 0)


def nu0(x0: Fraction) -> arb:
    """Exact lower bound of |gamma'| on (0, x0] (L4, N0)."""
    if not n0_validity(x0):
        raise AssertionError("N0 hypothesis x0 < exp(-(1+log2)/2) is not certified")
    x = frac_arb(x0)
    return ((2 - 3 * x) / (x * (2 * (ONE / x).log() + ONE - LOG2))).lower()


def nu1(x1: Fraction) -> arb:
    """Exact lower bound of |gamma'| on [x1, 1) (L4, N1)."""
    if not x1 >= Fraction(3, 4):
        raise AssertionError("N1 hypothesis x1 >= 3/4 violated")
    one_minus_p = frac_arb(1 - pid_frac(x1))
    ell = (ONE / one_minus_p).log()
    return (frac_arb(pidp_frac(x1)) / one_minus_p * ell / (ell + ONE)).lower()


def h_max_upper(lo: Fraction, hi: Fraction) -> arb:
    """Exact upper bound of h on [lo,hi] in (0,1): h increases on [0,1/2] and decreases on [1/2,1]."""
    if lo <= HALF <= hi:
        return LOG2.upper()
    return exact_max([h_frac(lo).upper(), h_frac(hi).upper()])


def dgamma_lower_middle(l: Fraction, r: Fraction) -> arb:
    """L5: inf of |gamma'| on [l,r] in (0,1) via monotone numerator."""
    if not (0 < l <= r < 1):
        raise AssertionError("middle piece must be strictly interior")
    pl, pr = pid_frac(l), pid_frac(r)
    numerator = (-(ONE - frac_arb(pl)).log()) * frac_arb(pidp_frac(l))
    denominator = frac_arb(pr) * h_max_upper(pl, pr)
    return (numerator / denominator).lower()


def dgamma_lower(l: Fraction, r: Fraction) -> arb:
    """Exact lower bound of |gamma'(xi)| for xi in [l,r] cap (0,1), 0<=l<r<=1."""
    if not (0 <= l < r <= 1):
        raise AssertionError("bad derivative interval")
    pieces: list[arb] = []
    if l < X0:
        pieces.append(nu0(min(r, X0)))
    if r > X1:
        pieces.append(nu1(max(l, X1)))
    ml, mr = max(l, X0), min(r, X1)
    if ml < mr:
        pieces.append(dgamma_lower_middle(ml, mr))
    if not pieces:
        raise AssertionError("derivative interval not covered")
    return exact_min(pieces)


# --------------------------------------------------------------------------- #
# Cell cover of the triangle 0 < s < t < 1
# --------------------------------------------------------------------------- #

def cell_factor(cell: Cell, mutate: str | None) -> arb:
    if mutate == "drop_factor":
        return ONE
    return ONE + frac_arb((1 - cell.s_hi) * (1 - cell.t_hi))


def certify_cell(cell: Cell, theta: arb, mutate: str | None) -> tuple[str, arb] | None:
    factor = cell_factor(cell, mutate)
    d_bound = (factor * dgamma_lower(cell.s_lo, cell.t_hi)).lower()
    if d_bound > theta:
        return "D", d_bound
    if cell.s_hi < cell.t_lo:
        quotient = (gamma_point(cell.s_hi) - gamma_point(cell.t_lo)) / frac_arb(cell.t_hi - cell.s_lo)
        q_bound = (factor * quotient).lower()
        if q_bound > theta:
            return "Q", q_bound
    return None


def split_cell(cell: Cell) -> tuple[Cell, Cell]:
    depth = cell.depth + 1
    if cell.s_hi - cell.s_lo >= cell.t_hi - cell.t_lo:
        mid = (cell.s_lo + cell.s_hi) / 2
        return (
            Cell(cell.s_lo, mid, cell.t_lo, cell.t_hi, depth),
            Cell(mid, cell.s_hi, cell.t_lo, cell.t_hi, depth),
        )
    mid = (cell.t_lo + cell.t_hi) / 2
    return (
        Cell(cell.s_lo, cell.s_hi, cell.t_lo, mid, depth),
        Cell(cell.s_lo, cell.s_hi, mid, cell.t_hi, depth),
    )


def in_layer(cell: Cell) -> bool:
    """Does the cell meet the residual layer of liu9-h2-phi.json?"""
    a = LAYER_A
    return cell.s_lo < a or cell.t_lo < a or cell.s_hi > 1 - a or cell.t_hi > 1 - a


def cover_triangle(theta: Fraction, mutate: str | None = None, max_cells: int = MAX_CELLS) -> CoverResult:
    theta_arb = frac_arb(theta)
    result = CoverResult(ok=True)
    stack = [Cell(Fraction(0), Fraction(1), Fraction(0), Fraction(1))]
    while stack:
        cell = stack.pop()
        if cell.s_lo >= cell.t_hi:
            result.discarded += 1
            continue
        result.processed += 1
        result.max_depth = max(result.max_depth, cell.depth)
        if result.processed > max_cells:
            result.ok = False
            result.failure = f"cell budget {max_cells} exhausted"
            result.failing_cell = cell
            return result
        verdict = certify_cell(cell, theta_arb, mutate)
        if verdict is None:
            if cell.depth >= MAX_DEPTH:
                result.ok = False
                result.failure = f"max depth {MAX_DEPTH} reached"
                result.failing_cell = cell
                return result
            left, right = split_cell(cell)
            stack.append(right)
            stack.append(left)
            continue
        kind, bound = verdict
        if kind == "D":
            result.accepted_d += 1
        else:
            result.accepted_q += 1
        if in_layer(cell):
            result.layer_cells += 1
        if result.min_certified is None or bound < result.min_certified:
            result.min_certified = bound
            result.min_cell = cell
    return result


# --------------------------------------------------------------------------- #
# Sanity samples (DISCOVERY only; never load-bearing)
# --------------------------------------------------------------------------- #

def endpoint_samples() -> dict[str, object]:
    samples0 = [Fraction(1, 2 ** k) for k in (4, 8, 16, 32, 64)]
    samples1 = [1 - Fraction(1, 2 ** k) for k in (3, 8, 16, 32, 64)]
    bound0, bound1 = nu0(X0), nu1(X1)
    for x in samples0:
        if not abs(dgamma_point(x)).lower() >= bound0:
            raise AssertionError("N0 sample below its bound")
    for x in samples1:
        if not abs(dgamma_point(x)).lower() >= bound1:
            raise AssertionError("N1 sample below its bound")
    return {
        "nu0_at_X0": format_arb(bound0, 20),
        "nu1_at_X1": format_arb(bound1, 20),
        "samples_below_X0": {format_fraction(x): format_arb(abs(dgamma_point(x)), 12) for x in samples0},
        "samples_above_X1": {format_fraction(x): format_arb(abs(dgamma_point(x)), 12) for x in samples1},
        "note": "samples confirm the closed-form gamma' exceeds the analytic bounds; the bounds themselves are proved in L4",
    }


def kappa_point(s: Fraction, t: Fraction) -> arb:
    return (ONE + frac_arb((1 - s) * (1 - t))) * (gamma_point(s) - gamma_point(t)) / frac_arb(t - s)


def phi_quantitative_samples(theta: Fraction) -> dict[str, object]:
    """Check Phi >= (theta^2/4-1) s^2 t^2 (s-t)^2 g(s) g(t) at sample points (sanity)."""
    coefficient = frac_arb(theta * theta / 4 - 1)
    worst_ratio: arb | None = None
    for s, t in IDENTITY_POINTS:
        gs = lam_arb(frac_arb(pid_frac(s)))
        gt = lam_arb(frac_arb(pid_frac(t)))
        lower = coefficient * frac_arb(s * s * t * t * (s - t) ** 2) * gs * gt
        phi = phi_direct(s, t)
        if not (phi - lower).lower() >= 0:
            raise AssertionError(f"quantitative bound violated at {s},{t}")
        ratio = (phi / lower).lower()
        worst_ratio = ratio if worst_ratio is None else exact_min([worst_ratio, ratio])
        if not kappa_point(s, t).lower() > frac_arb(theta):
            raise AssertionError(f"kappa sample below theta at {s},{t}")
    return {
        "coefficient": format_fraction(theta * theta / 4 - 1),
        "min_ratio_Phi_over_bound_at_points": format_arb(worst_ratio, 15),
        "note": "sanity only; the bound is a theorem once the cover certifies kappa >= theta",
    }


# --------------------------------------------------------------------------- #
# beta > 0 from the binding solver (exact conversion, never via mpf())
# --------------------------------------------------------------------------- #

def exact_mpf_fraction(value: object) -> Fraction:
    sign, mantissa, exponent, _bits = value._mpf_  # type: ignore[attr-defined]
    signed = -mantissa if sign else mantissa
    if exponent >= 0:
        return Fraction(signed * (1 << exponent), 1)
    return Fraction(signed, 1 << (-exponent))


def binding_provenance() -> dict[str, object]:
    """Record beta, m, xstar from liu9_binding.solve_equation_parameters(100).

    All three are converted from the solver's stored binary mpf to EXACT
    rationals (never re-rounded through mpf()).  Only beta > 0 is used, as the
    prefactor of B; m and xstar are recorded for audit and do not enter Phi.
    """
    from liu9_binding import solve_equation_parameters  # noqa: E402

    parameters = solve_equation_parameters(100)
    beta = exact_mpf_fraction(parameters.beta)
    mean = exact_mpf_fraction(parameters.mean)
    xstar = exact_mpf_fraction(parameters.x)
    if not beta > 0:
        raise AssertionError("beta not positive")
    return {
        "source": "liu9_binding.solve_equation_parameters(100), fields beta, mean, x",
        "conversion": "mpf -> exact Fraction via the stored (sign, mantissa, exponent); no mpf() re-rounding",
        "beta_exact_binary_fraction_bits": beta.denominator.bit_length() - 1,
        "beta_decimal_40": format_arb(frac_arb(beta), 40),
        "m_decimal_40": format_arb(frac_arb(mean), 40),
        "xstar_decimal_40": format_arb(frac_arb(xstar), 40),
        "role": "beta is only the positive prefactor in B = beta*(h(pi_st) - sqrt(h(pi_ss) h(pi_tt))); Phi itself is parameter-free; m and xstar are unused",
        "status": "PROVED (beta > 0 as an exact rational)",
    }


# --------------------------------------------------------------------------- #
# Mutations
# --------------------------------------------------------------------------- #

def run_mutations(legit_cells: int) -> list[dict[str, object]]:
    mutations: list[dict[str, object]] = []
    budget = max(4 * legit_cells, 20_000)

    # M1: threshold above the numerical minimum 2.394 must NOT be certifiable.
    # The mutant claim is moreover REFUTED by a certified point with kappa < 12/5.
    m1 = cover_triangle(Fraction(12, 5), None, max_cells=budget)
    w1s, w1t = Fraction(27, 50), Fraction(27, 50) + Fraction(1, 10 ** 6)
    w1 = kappa_point(w1s, w1t)
    refuted1 = bool(w1.upper() < frac_arb(Fraction(12, 5)))
    mutations.append({
        "mutation": "theta_2.4_above_true_minimum",
        "theta": "12/5",
        "caught": (not m1.ok) and refuted1,
        "observed": m1.failure or "unexpectedly certified",
        "processed_cells": m1.processed,
        "failing_cell": cell_record(m1.failing_cell),
        "refutation_witness": {
            "point": [format_fraction(w1s), format_fraction(w1t)],
            "kappa_enclosure": format_arb(w1, 25),
            "certified_below_12/5": refuted1,
        },
    })

    # M2: dropping the (1+ab) factor must fail even at theta = 2.
    # Refuted by a certified point with |gamma'| < 2.
    m2 = cover_triangle(Fraction(2), "drop_factor", max_cells=budget)
    w2 = abs(dgamma_point(Fraction(9, 20)))
    refuted2 = bool(w2.upper() < TWO)
    mutations.append({
        "mutation": "drop_(1+(1-s)(1-t))_factor",
        "theta": "2",
        "caught": (not m2.ok) and refuted2,
        "observed": m2.failure or "unexpectedly certified",
        "processed_cells": m2.processed,
        "failing_cell": cell_record(m2.failing_cell),
        "refutation_witness": {
            "point": format_fraction(Fraction(9, 20)),
            "abs_dgamma_enclosure": format_arb(w2, 25),
            "certified_below_2": refuted2,
        },
        "reason": "the diagonal limit of the factor-free condition is |gamma'(x)|>=2, which fails at x=9/20",
    })

    # M3: the identity with the wrong protocol s+t-st must fail.
    caught3 = False
    try:
        identity_checks("wrong_protocol")
    except AssertionError as error:
        caught3 = True
        observed3 = str(error)[:160]
    mutations.append({
        "mutation": "identity_with_protocol_s+t-st",
        "caught": caught3,
        "observed": observed3 if caught3 else "identity unexpectedly held",
    })

    # M4: the identity with (Lam_ss-Lam_tt)^2/2 instead of /4 must fail.
    caught4 = False
    try:
        identity_checks("wrong_quarter")
    except AssertionError as error:
        caught4 = True
        observed4 = str(error)[:160]
    mutations.append({
        "mutation": "identity_quarter_replaced_by_half",
        "caught": caught4,
        "observed": observed4 if caught4 else "identity unexpectedly held",
    })

    # M5: N0 hypothesis outside its validity range must be rejected.
    caught5 = False
    try:
        nu0(Fraction(1, 2))
    except AssertionError as error:
        caught5 = True
        observed5 = str(error)
    mutations.append({
        "mutation": "N0_applied_at_x0=1/2_beyond_validity",
        "caught": caught5,
        "observed": observed5 if caught5 else "hypothesis check unexpectedly passed",
    })
    return mutations


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def cell_record(cell: Cell | None) -> dict[str, object] | None:
    if cell is None:
        return None
    return {
        "s": [format_fraction(cell.s_lo), format_fraction(cell.s_hi)],
        "t": [format_fraction(cell.t_lo), format_fraction(cell.t_hi)],
        "depth": cell.depth,
    }


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def build_report() -> dict[str, object]:
    polynomial = polynomial_identity_exact()
    derivative = derivative_identity_check()
    identities = identity_checks()
    endpoints = endpoint_samples()
    if not n0_validity(X0):
        raise AssertionError("X0 outside N0 validity")
    cover = cover_triangle(THETA)
    if not cover.ok:
        raise AssertionError(f"cover failed: {cover.failure} at {cell_record(cover.failing_cell)}")
    if not cover.min_certified > frac_arb(THETA):
        raise AssertionError("certified minimum does not exceed theta")
    quantitative = phi_quantitative_samples(THETA)
    binding = binding_provenance()
    mutations = run_mutations(cover.processed)
    all_caught = all(bool(m["caught"]) for m in mutations)
    if not all_caught:
        raise AssertionError("a mutation was not caught")
    if polynomial["status"] != "PROVED":
        raise AssertionError("polynomial identity failed")

    layer_area = Fraction(16383, 67108864)
    report: dict[str, object] = {
        "tool": "liu9_h2_boundary_layer.py (python-flint Arb)",
        "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "precision_bits": PRECISION_BITS,
        "natural_log": True,
        "single_core_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "FLINT_NUM_THREADS")},
        "claim": "Phi(s,t) = h(pi(s,t))^2 - h(pi(s,s))*h(pi(t,t)) >= 0 for all (s,t) in [0,1]^2, with pi(s,t)=s*t*(1+(1-s)*(1-t))",
        "claim_status": "PROVED",
        "parameter_independence": "Phi contains no parameter; beta, m, xstar do not enter the certificate.  beta>0 is recorded only to conclude B>=0.",
        "quantitative_form": {
            "statement": "for 0<s<t<1: Phi >= (theta^2/4 - 1) * s^2 t^2 (s-t)^2 g(s) g(t) with theta = 9/4, i.e. coefficient 17/64, where g(x) = h(pi(x,x))/pi(x,x) > 0",
            "sanity": quantitative,
        },
        "definitions": {
            "h": "-u*log(u) - (1-u)*log(1-u), h(0)=h(1)=0",
            "pi": "s*t*(1+(1-s)*(1-t))",
            "Lam": "h(u)/u = log(1/u) + mu(u)",
            "mu": "-(1-u)*log(1-u)/u",
            "pi_d": "pi(x,x) = x^2*(2-2x+x^2)",
            "g": "Lam(pi_d(x))",
            "gamma": "log g(x)",
            "kappa": "(1+(1-s)(1-t)) * (gamma(s)-gamma(t))/(t-s) for 0<s<t<1",
        },
        "lemmas": {
            "L0a_derivative_identity": derivative,
            "L0b_pi_d_monotone": {
                "statement": "pi_d'(x)=2x(2-3x+2x^2)>0 on (0,1]; pi_d''(x)=4(1-3x+3x^2)>0",
                "proof": "discriminants 9-16=-7<0 and 9-12=-3<0 with positive leading coefficients",
                "status": "PROVED",
            },
            "L0c_gamma_decreasing": {
                "statement": "g = Lam o pi_d and gamma = log g are strictly decreasing on (0,1); gamma' = log(1-pi_d) pi_d' /(pi_d h(pi_d))",
                "proof": "Lam' < 0 by L0a, pi_d' > 0 by L0b, chain rule",
                "status": "PROVED",
            },
            "L0d_mu_shape": {
                "statement": "mu(u) = 1 - sum_{k>=2} u^{k-1}/(k(k-1)); mu decreasing, concave, 0<=mu<=1 on [0,1]",
                "proof": "-(1-u)log(1-u) = u - sum_{k>=2} u^k/(k(k-1)) by multiplying the log series by (1-u); all derivative coefficients are negative on (0,1); mu(1)=0 by continuity; mu>=0 since -(1-u)log(1-u)>=0",
                "status": "PROVED",
            },
            "L0e_polynomial_identity": polynomial,
            "L0f_M_mu_nonpositive": {
                "statement": "mu(pi_ss)+mu(pi_tt)-2mu(pi_st) <= 0",
                "proof": "L0e gives pi_st <= sqrt(pi_ss pi_tt) <= (pi_ss+pi_tt)/2; mu decreasing then mu concave (L0d)",
                "status": "PROVED",
            },
            "L1_identity": identities,
            "L2_master": {
                "statement": "for 0<s<t<1, Phi >= s^2 t^2 (s-t)^2 g(s) g(t) (kappa^2/4 - 1); kappa>=2 implies Phi>=0",
                "proof": "drop the nonnegative middle term of L1; g(s)-g(t) = 2 sqrt(g(s)g(t)) sinh((gamma(s)-gamma(t))/2); sinh x >= x for x >= 0",
                "status": "PROVED",
            },
            "L3_boundary": {
                "statement": "Phi(s,s)=0; Phi(0,t)=0; Phi(s,1)=h(s)^2>=0 since pi(s,1)=s and pi(1,1)=1; Phi symmetric",
                "status": "PROVED",
            },
            "L4_endpoint_bounds": {
                "N0": {
                    "statement": "0<x<=x0<exp(-(1+log2)/2): |gamma'(x)| >= (2-3x0)/(x0(2log(1/x0)+1-log2))",
                    "x0_used": format_fraction(X0),
                    "validity_certified": n0_validity(X0),
                    "value": endpoints["nu0_at_X0"],
                },
                "N1": {
                    "statement": "3/4<=x1<=x<1: |gamma'(x)| >= pi_d'(x1)/(1-pi_d(x1)) * ell1/(ell1+1), ell1=log(1/(1-pi_d(x1)))",
                    "x1_used": format_fraction(X1),
                    "value": endpoints["nu1_at_X1"],
                },
                "samples": endpoints,
                "status": "PROVED (proofs in module docstring)",
            },
            "L5_interior_bound": {
                "statement": "N(x)=-log(1-pi_d(x))pi_d'(x) is positive increasing; inf_[l,r]|gamma'| >= N(l)/(pi_d(r) max_[pi_d(l),pi_d(r)] h)",
                "status": "PROVED",
            },
        },
        "cell_cover": {
            "region": "closed triangle 0<=s<=t<=1; only points with 0<s<t<1 need kappa (L3 handles the rest); cells with s_lo>=t_hi are discarded",
            "theta": format_fraction(THETA),
            "acceptance": "lower endpoint of the Arb bound (D) or (Q) strictly exceeds theta",
            "processed_cells": cover.processed,
            "accepted_by_D": cover.accepted_d,
            "accepted_by_Q": cover.accepted_q,
            "discarded_cells": cover.discarded,
            "max_depth": cover.max_depth,
            "min_certified_bound": format_arb(cover.min_certified, 30),
            "min_certified_cell": cell_record(cover.min_cell),
            "accepted_cells_meeting_residual_layer": cover.layer_cells,
            "numerical_minimum_of_kappa_DISCOVERY": "2.394 at s=t=0.54 (float scan, not load-bearing)",
        },
        "residual_layer_of_liu9_h2_phi": {
            "set": "([0,a) x [0,1] union (1-a,1] x [0,1] union [a,1-a] x ([0,a) union (1-a,1])) minus ({s=0} union {t=0} union {s=t})",
            "a": format_fraction(LAYER_A),
            "exact_area": format_fraction(layer_area),
            "status": "PROVED by inclusion in [0,1]^2; no separate endpoint remainder lemma is needed because kappa -> +infinity at every boundary (L4)",
        },
        "consequences": {
            "B_nonnegative": "B = beta*(h(pi_st) - sqrt(h(pi_ss) h(pi_tt))) >= 0 on [0,1]^2 (beta>0 below)",
            "binding_parameters": binding,
            "combined_with": [
                "uc/verification/results/liu9-h2-reduction.json (chain C1-C6)",
                "uc/verification/results/liu9-h2-twovar.json (A>=0 and R>=0 on [0,1]^2)",
            ],
            "yields": "T >= 0 hence gap >= 0 for every q, every mass vector, every admissible M >= m, arbitrarily many atoms, with NO support restriction",
        },
        "mutations": mutations,
        "all_mutations_caught": all_caught,
        "bound_directions": {
            "identity_middle_term": "dropped term is (log(1+delta)-M_mu)/2*(Lam_st+Lbar) with log(1+delta)>=0, -M_mu>=0 (L0f), Lam>=0 (L0g): dropping it lowers Phi",
            "D_rule": "MVT: (gamma(s)-gamma(t))/(t-s) = |gamma'(xi)| >= inf over [s_lo,t_hi]; (1+ab) >= 1+(1-s_hi)(1-t_hi); product of lower bounds of nonnegative quantities",
            "Q_rule": "gamma(s)>=gamma(s_hi), gamma(t)<=gamma(t_lo) (gamma decreasing), 0<t-s<=t_hi-s_lo; requires s_hi<t_lo so the numerator lower bound is positive",
            "N0": "numerator uses -log(1-pi)>=pi and pi_d'>=2x(2-3x0); denominator uses h(pi)<=pi(log(1/pi)+1), pi_d<=2x^2, u(log(1/u)+1) increasing, x(2log(1/x)+c) increasing below exp(-(1+log2)/2)",
            "N1": "numerator uses pi_d'>=pi_d'(x1) on [x1,1] (x1>=3/4) and ell/(ell+1) increasing; denominator uses pi<=1 and h(pi)<=(1-pi)(log(1/(1-pi))+1)",
            "L5": "numerator N is increasing so N(l) is its minimum; denominator pi_d h(pi_d) <= pi_d(r) * max h on the image interval",
            "sinh": "sinh x >= x for x >= 0 gives kappa_sinh >= kappa/2; used in L2 only",
        },
        "determinism": "no timestamps; deterministic depth-first traversal; Arb strings at fixed digit counts",
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON with report_sha256 omitted",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def write_report(report: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=1) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    default_output = _HERE / "verification" / "results" / "liu9-h2-boundary.json"
    parser = argparse.ArgumentParser(description="Certify Phi >= 0 on [0,1]^2 (Liu H2, B >= 0).")
    parser.add_argument("--output", type=Path, default=default_output)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report()
    write_report(report, args.output)
    cover = report["cell_cover"]
    print(f"PHI_NONNEGATIVE_ON_UNIT_SQUARE PROVED theta={cover['theta']} cells={cover['processed_cells']} depth={cover['max_depth']}")
    print(f"RESIDUAL_LAYER PROVED (by inclusion)")
    print(f"MUTATIONS all_caught={report['all_mutations_caught']}")
    print(f"REPORT {args.output}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
