#!/usr/bin/env python3
"""Rigorous Arb certificates for three two-variable Liu H2 lemmas.

The certified claims are, on the closed unit square,

    R(s,t) = P2(s,t) + beta*h(pi(s,t)) >= 0,
    Lam(s,t) = 2*P2(s,t) + beta*h(pi(s,s)) + beta*h(pi(t,t)) >= 0,
    A2(s,t) = P2(s,t) + phi(s)*phi(t) >= 0,

where phi(s)=sqrt(beta*h(pi(s,s))).

The proof is deliberately stratified.  Endpoint faces and neighborhoods are
handled analytically, the common nondegenerate zero (x*,x*) by a Taylor bound,
and the compact remainder by deterministic Arb branch-and-bound.  A box that
touches a zero is never passed to the global positivity cover.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable, Iterable

# Keep the certificate single-core even when imported by a launcher that has not
# already constrained numerical libraries.  This file itself imports no NumPy.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_variable, "1")

# `python -I` deliberately omits the script's own directory from sys.path, so the
# deferred `liu9_binding` import in binding_parameter_crosscheck() would fail
# under the invocation used everywhere else in this repository
# (`./.venv/bin/python -I -B uc/<module>.py`).  Added by Main as a purely
# mechanical reproducibility fix; no mathematics is affected.
import sys  # noqa: E402

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from flint import arb, ctx  # noqa: E402


PRECISION_BITS = 400
ctx.prec = PRECISION_BITS

ZERO = arb(0)
ONE = arb(1)
TWO = arb(2)
HALF = ONE / TWO
LOG2 = TWO.log()
SQRT2 = TWO.sqrt()

EPSILON = Fraction(1, 100_000_000)       # near s=0 or t=0
CORNER_DELTA = Fraction(1, 10_000)       # near (1,1), R and Lam
CORNER_DELTA_A2 = Fraction(1, 100_000_000)
LOCAL_HALF_A = Fraction(18, 1000)
LOCAL_HALF_B = Fraction(16, 1000)
LOCAL_HALF_A2 = Fraction(16, 1000)
TAYLOR_BOX_RADIUS = arb("0.03")
TAYLOR_BOX_RADIUS_A2 = arb("0.025")
MAX_BB_DEPTH = 48
MAX_BB_DEPTH_A2 = 64
MAX_BB_CELLS = 1_000_000


@dataclass(frozen=True)
class Box:
    s_lo: Fraction
    s_hi: Fraction
    t_lo: Fraction
    t_hi: Fraction
    depth: int = 0


@dataclass
class CoverStats:
    processed: int = 0
    accepted: int = 0
    split: int = 0
    max_depth: int = 0
    worst_lower: arb | None = None
    worst_box: Box | None = None


@dataclass(frozen=True)
class Parameters:
    x: arb
    p: arb
    m: arb
    beta: arb
    root_lo: Fraction
    root_hi: Fraction
    root_polynomial_lo: arb
    root_polynomial_hi: arb


@dataclass(frozen=True)
class TaylorCertificate:
    value: arb
    grad_s: arb
    grad_t: arb
    h_ss: arb
    h_st: arb
    h_tt: arb
    lambda_lower: arb
    third_partial_enclosures: tuple[arb, arb, arb, arb]
    max_third_upper: arb
    c3_upper: arb
    rho_formula: arb
    derivative_box_radius: arb
    chosen_half_width: Fraction
    chosen_distance_upper: arb
    coefficient_lower: arb
    point_c3_required_lower: arb


def arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator

def exact_mpf_fraction(value: object) -> Fraction:
    """Convert an mpmath mpf's stored binary value to an exact rational."""
    sign, mantissa, exponent, _bit_count = value._mpf_  # type: ignore[attr-defined]
    signed_mantissa = -mantissa if sign else mantissa
    if exponent >= 0:
        return Fraction(signed_mantissa << exponent)
    return Fraction(signed_mantissa, 1 << (-exponent))


def hull(*values: arb) -> arb:
    if not values:
        raise ValueError("hull requires at least one value")
    result = arb(values[0])
    for value in values[1:]:
        result = result.union(arb(value))
    return result


def min_enclosure(a: arb, b: arb) -> arb:
    """An enclosure no larger than either argument, used for minima."""
    if a <= b:
        return a
    if b <= a:
        return b
    return a.union(b).lower()


def format_arb(value: arb, digits: int = 50) -> str:
    return value.str(digits)


def format_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def polynomial(z: arb) -> arb:
    return z**4 - 2 * z**3 + 3 * z**2 - 1


def h_interior(u: arb) -> arb:
    """Natural-log binary entropy for an interval strictly in (0,1)."""
    return -(u * u.log() + (ONE - u) * (ONE - u).log())


def h_point_fraction(u: Fraction) -> arb:
    if u == 0 or u == 1:
        return ZERO
    if not 0 < u < 1:
        raise ValueError(f"entropy point outside [0,1]: {u}")
    return h_interior(arb_fraction(u))


def h_range_fraction(lo: Fraction, hi: Fraction) -> arb:
    """Exact range enclosure of natural-log binary entropy on [lo,hi]."""
    if not 0 <= lo <= hi <= 1:
        raise ValueError(f"entropy interval outside [0,1]: [{lo},{hi}]")
    result = hull(h_point_fraction(lo), h_point_fraction(hi))
    if lo <= Fraction(1, 2) <= hi:
        result = result.union(LOG2)
    return result


def hp(u: arb) -> arb:
    return ((ONE - u) / u).log()


def hpp(u: arb) -> arb:
    return -ONE / (u * (ONE - u))


def hppp(u: arb) -> arb:
    return (ONE - 2 * u) / (u**2 * (ONE - u) ** 2)


def pi_arb(s: arb, t: arb) -> arb:
    return s * t * (ONE + (ONE - s) * (ONE - t))


def pi_fraction(s: Fraction, t: Fraction) -> Fraction:
    return s * t * (1 + (1 - s) * (1 - t))


def qdiag_fraction(s: Fraction) -> Fraction:
    return pi_fraction(s, s)


def interval_fraction(lo: Fraction, hi: Fraction) -> arb:
    return hull(arb_fraction(lo), arb_fraction(hi))


def isolate_parameters() -> Parameters:
    """Isolate the unique root and propagate it through Liu's equations."""
    root_lo = Fraction(69, 100)
    root_hi = Fraction(691, 1000)
    for _ in range(220):
        midpoint = (root_lo + root_hi) / 2
        sign = polynomial(arb_fraction(midpoint))
        if sign < 0:
            root_lo = midpoint
        elif sign > 0:
            root_hi = midpoint
        else:
            raise RuntimeError("Arb could not determine a bisection sign")

    lo = arb_fraction(root_lo)
    hi = arb_fraction(root_hi)
    f_lo = polynomial(lo)
    f_hi = polynomial(hi)
    if not (f_lo < 0 < f_hi):
        raise AssertionError("root endpoints do not have opposite certified signs")

    # f'(z)=2z(2(z-3/4)^2+15/8)>0 for 0<z<=1, so the root is unique.
    x = lo.union(hi)
    hx = h_interior(x)
    hxx = h_interior(x * x)
    p = hx / hxx
    mean = p * x
    diagonal_protocol = pi_arb(x, x)
    diagonal_protocol_prime = (
        2 * x + 2 * x * (ONE - x) ** 2
        - 2 * x * x * (ONE - x)
    )
    a_prime = hp(x * x) * 2 * x
    b_prime = hp(diagonal_protocol) * diagonal_protocol_prime
    beta = ((hp(x) + hx / x) / p - a_prime) / (b_prime - a_prime)

    if not (x > 0 and x < 1 and p > 0 and mean > 0 and beta > 0 and beta < 1):
        raise AssertionError("derived parameter sign/range check failed")
    return Parameters(x, p, mean, beta, root_lo, root_hi, f_lo, f_hi)


PARAMETERS = isolate_parameters()
XSTAR = PARAMETERS.x
# Put the rational partition center just to the left of the certified root
# bracket.  This keeps mutant bisections from cutting through the bracket while
# making the center displacement itself rigorously controlled by the bracket
# width, rather than by a hand-copied decimal.
CENTER = 2 * PARAMETERS.root_lo - PARAMETERS.root_hi
P = PARAMETERS.p
M = PARAMETERS.m
BETA = PARAMETERS.beta
C0 = ONE - BETA - ONE / (2 * M)
C1 = ONE - ONE / (2 * M)

if not (C0 > 0 and C1 > 0):
    raise AssertionError("face coefficients are not certified positive")

def binding_parameter_crosscheck() -> dict[str, object]:
    """Compare the rigorous closed forms with the binding solver's stored mpfs."""
    from liu9_binding import solve_equation_parameters

    binding = solve_equation_parameters(100)
    binding_x_fraction = exact_mpf_fraction(binding.x)
    binding_mean_fraction = exact_mpf_fraction(binding.mean)
    binding_x = arb_fraction(binding_x_fraction)
    binding_mean = arb_fraction(binding_mean_fraction)
    x_difference = XSTAR - binding_x
    mean_difference = M - binding_mean
    if not (x_difference.contains(0) and mean_difference.contains(0)):
        raise AssertionError(
            "liu9_binding constants disagree with the rigorously isolated "
            f"quartic parameters: x difference {x_difference}; "
            f"mean difference {mean_difference}"
        )
    return {
        "solver": "liu9_binding.solve_equation_parameters(100)",
        "binding_x_exact_stored_mpf_fraction": format_fraction(binding_x_fraction),
        "binding_mean_exact_stored_mpf_fraction": format_fraction(
            binding_mean_fraction
        ),
        "binding_x_enclosure": format_arb(binding_x),
        "binding_mean_enclosure": format_arb(binding_mean),
        "xstar_minus_binding_x_enclosure": format_arb(x_difference),
        "closed_form_mean_minus_binding_mean_enclosure": format_arb(
            mean_difference
        ),
        "status": "AGREES: both certified Arb difference enclosures contain zero",
    }


def p2_point(s: Fraction, t: Fraction) -> arb:
    st = s * t
    return (
        (ONE - BETA) * h_point_fraction(st)
        - (
            arb_fraction(t) * h_point_fraction(s)
            + arb_fraction(s) * h_point_fraction(t)
        ) / (2 * M)
    )


def k_point(s: Fraction, t: Fraction) -> arb:
    return h_point_fraction(pi_fraction(s, t))


def r_point(s: Fraction, t: Fraction) -> arb:
    return p2_point(s, t) + BETA * k_point(s, t)


def lam_point(s: Fraction, t: Fraction) -> arb:
    return (
        2 * p2_point(s, t)
        + BETA * k_point(s, s)
        + BETA * k_point(t, t)
    )


def phi_point(s: Fraction) -> arb:
    if s == 0 or s == 1:
        return ZERO
    return (BETA * k_point(s, s)).sqrt()


def a2_point(s: Fraction, t: Fraction) -> arb:
    return p2_point(s, t) + phi_point(s) * phi_point(t)


def mu_interior(u: arb) -> arb:
    return -(ONE - u) * (ONE - u).log() / u


def rcheck_point(s: Fraction, t: Fraction) -> arb:
    if not (0 < s <= 1 and 0 < t <= 1):
        raise ValueError("Rcheck requires positive coordinates")
    sa, ta = arb_fraction(s), arb_fraction(t)
    rho = ONE + (ONE - sa) * (ONE - ta)
    st = sa * ta
    return (
        st.log() * (HALF - M - M * BETA * (ONE - sa) * (ONE - ta))
        + M * (ONE - BETA) * mu_interior(st)
        + M * BETA * rho * mu_interior(st * rho)
        - M * BETA * rho * rho.log()
        - (mu_interior(sa) + mu_interior(ta)) / 2
    ) / M


def r_range(box: Box) -> arb:
    s = interval_fraction(box.s_lo, box.s_hi)
    t = interval_fraction(box.t_lo, box.t_hi)
    h_st = h_range_fraction(box.s_lo * box.t_lo, box.s_hi * box.t_hi)
    h_s = h_range_fraction(box.s_lo, box.s_hi)
    h_t = h_range_fraction(box.t_lo, box.t_hi)
    p_lo = pi_fraction(box.s_lo, box.t_lo)
    p_hi = pi_fraction(box.s_hi, box.t_hi)
    h_pi = h_range_fraction(p_lo, p_hi)
    return (
        (ONE - BETA) * h_st
        - (t * h_s + s * h_t) / (2 * M)
        + BETA * h_pi
    )


def lam_range(box: Box) -> arb:
    s = interval_fraction(box.s_lo, box.s_hi)
    t = interval_fraction(box.t_lo, box.t_hi)
    h_st = h_range_fraction(box.s_lo * box.t_lo, box.s_hi * box.t_hi)
    h_s = h_range_fraction(box.s_lo, box.s_hi)
    h_t = h_range_fraction(box.t_lo, box.t_hi)
    q_s = h_range_fraction(qdiag_fraction(box.s_lo), qdiag_fraction(box.s_hi))
    q_t = h_range_fraction(qdiag_fraction(box.t_lo), qdiag_fraction(box.t_hi))
    return (
        2 * (ONE - BETA) * h_st
        - (t * h_s + s * h_t) / M
        + BETA * q_s
        + BETA * q_t
    )


def sqrt_nonnegative_range(value: arb) -> arb:
    """Monotone square-root range for a quantity known a priori nonnegative."""
    upper = value.upper()
    if upper < 0:
        raise AssertionError("nonnegative radicand has negative upper endpoint")
    lower = value.lower()
    sqrt_lower = lower.sqrt() if lower > 0 else ZERO
    sqrt_upper = upper.sqrt() if upper > 0 else ZERO
    return hull(sqrt_lower, sqrt_upper)


def phi_range(lo: Fraction, hi: Fraction) -> arb:
    diagonal_entropy = h_range_fraction(
        qdiag_fraction(lo), qdiag_fraction(hi)
    )
    return sqrt_nonnegative_range(BETA * diagonal_entropy)


def a2_range(box: Box) -> arb:
    s = interval_fraction(box.s_lo, box.s_hi)
    t = interval_fraction(box.t_lo, box.t_hi)
    h_st = h_range_fraction(box.s_lo * box.t_lo, box.s_hi * box.t_hi)
    h_s = h_range_fraction(box.s_lo, box.s_hi)
    h_t = h_range_fraction(box.t_lo, box.t_hi)
    p2 = (
        (ONE - BETA) * h_st
        - (t * h_s + s * h_t) / (2 * M)
    )
    return p2 + phi_range(box.s_lo, box.s_hi) * phi_range(
        box.t_lo, box.t_hi
    )


def logarithmic_product_lower(box: Box, multiplier: arb) -> arb | None:
    """Return the analytic entropy lower bound when its coefficient is positive.

    For positive ``s,t``, both lemmas use

        c0*log(1/(s*t)) - 1/m.

    The logarithm is smallest at the upper product endpoint.  Once that
    coefficient is positive, multiplying it by the lower product endpoint is
    the correct lower-bound direction.
    """
    product_lo = box.s_lo * box.t_lo
    product_hi = box.s_hi * box.t_hi
    coefficient = (
        C0 * (ONE / arb_fraction(product_hi)).log() - ONE / M
    ).lower()
    if not coefficient > 0:
        return None
    return (multiplier * arb_fraction(product_lo) * coefficient).lower()


def stronger_lower(a: arb, b: arb) -> arb:
    """Return the stronger of two certified scalar lower bounds."""
    if a >= b:
        return a
    if b >= a:
        return b
    return a.union(b).lower()


def phi_prime_range(lo: Fraction, hi: Fraction) -> arb:
    coordinate = interval_fraction(lo, hi)
    diagonal_protocol = interval_fraction(
        qdiag_fraction(lo), qdiag_fraction(hi)
    )
    qprime = (
        4 * coordinate - 6 * coordinate**2 + 4 * coordinate**3
    )
    return (
        BETA * hp(diagonal_protocol) * qprime
        / (2 * phi_range(lo, hi))
    )


def gradient_lower(box: Box, lemma: str) -> arb | None:
    """Centered first-order lower bound using Arb gradient ranges."""
    if box.s_hi == 1 or box.t_hi == 1:
        return None
    s = interval_fraction(box.s_lo, box.s_hi)
    t = interval_fraction(box.t_lo, box.t_hi)
    st = interval_fraction(
        box.s_lo * box.t_lo, box.s_hi * box.t_hi
    )
    protocol = interval_fraction(
        pi_fraction(box.s_lo, box.t_lo),
        pi_fraction(box.s_hi, box.t_hi),
    )
    h_s = h_range_fraction(box.s_lo, box.s_hi)
    h_t = h_range_fraction(box.t_lo, box.t_hi)
    p2_s = (
        (ONE - BETA) * hp(st) * t
        - (t * hp(s) + h_t) / (2 * M)
    )
    p2_t = (
        (ONE - BETA) * hp(st) * s
        - (s * hp(t) + h_s) / (2 * M)
    )
    pi_s = t * (2 - 2 * s - t + 2 * s * t)
    pi_t = s * (2 - 2 * t - s + 2 * s * t)
    if lemma == "R":
        grad_s = p2_s + BETA * hp(protocol) * pi_s
        grad_t = p2_t + BETA * hp(protocol) * pi_t
        midpoint_value = r_point(
            (box.s_lo + box.s_hi) / 2,
            (box.t_lo + box.t_hi) / 2,
        )
    elif lemma == "Lam":
        q_s = interval_fraction(
            qdiag_fraction(box.s_lo), qdiag_fraction(box.s_hi)
        )
        q_t = interval_fraction(
            qdiag_fraction(box.t_lo), qdiag_fraction(box.t_hi)
        )
        qprime_s = 4 * s - 6 * s**2 + 4 * s**3
        qprime_t = 4 * t - 6 * t**2 + 4 * t**3
        grad_s = 2 * p2_s + BETA * hp(q_s) * qprime_s
        grad_t = 2 * p2_t + BETA * hp(q_t) * qprime_t
        midpoint_value = lam_point(
            (box.s_lo + box.s_hi) / 2,
            (box.t_lo + box.t_hi) / 2,
        )
    elif lemma == "A2":
        grad_s = (
            p2_s
            + phi_prime_range(box.s_lo, box.s_hi)
            * phi_range(box.t_lo, box.t_hi)
        )
        grad_t = (
            p2_t
            + phi_range(box.s_lo, box.s_hi)
            * phi_prime_range(box.t_lo, box.t_hi)
        )
        midpoint_value = a2_point(
            (box.s_lo + box.s_hi) / 2,
            (box.t_lo + box.t_hi) / 2,
        )
    else:
        raise ValueError(lemma)
    if not (grad_s.is_finite() and grad_t.is_finite()):
        return None
    half_s = arb_fraction((box.s_hi - box.s_lo) / 2)
    half_t = arb_fraction((box.t_hi - box.t_lo) / 2)
    return (
        midpoint_value.lower()
        - grad_s.abs_upper() * half_s
        - grad_t.abs_upper() * half_t
    ).lower()


def r_global_lower(box: Box) -> arb:
    analytic = logarithmic_product_lower(box, ONE)
    if analytic is not None:
        return analytic
    raw = r_range(box).lower()
    centered = gradient_lower(box, "R")
    return raw if centered is None else stronger_lower(raw, centered)


def lam_global_lower(box: Box) -> arb:
    analytic = logarithmic_product_lower(box, TWO)
    if analytic is not None:
        return analytic
    raw = lam_range(box).lower()
    centered = gradient_lower(box, "Lam")
    return raw if centered is None else stronger_lower(raw, centered)


def a2_global_lower(box: Box) -> arb:
    analytic = logarithmic_product_lower(box, ONE)
    if analytic is not None:
        return analytic
    raw = a2_range(box).lower()
    centered = gradient_lower(box, "A2")
    return raw if centered is None else stronger_lower(raw, centered)


def pi_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    value = pi_arb(s, t)
    ps = t * (2 - 2 * s - t + 2 * s * t)
    pt = s * (2 - 2 * t - s + 2 * s * t)
    pss = -2 * t * (ONE - t)
    ptt = -2 * s * (ONE - s)
    pst = 2 - 2 * s - 2 * t + 4 * s * t
    psst = -2 + 4 * t
    pstt = -2 + 4 * s
    return value, ps, pt, pss, pst, ptt, psst, pstt


def p2_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    u = s * t
    hu = h_interior(u)
    h1u, h2u, h3u = hp(u), hpp(u), hppp(u)
    hs, ht = h_interior(s), h_interior(t)
    value = (ONE - BETA) * hu - (t * hs + s * ht) / (2 * M)
    ds = (ONE - BETA) * h1u * t - (t * hp(s) + ht) / (2 * M)
    dt = (ONE - BETA) * h1u * s - (s * hp(t) + hs) / (2 * M)
    dss = (ONE - BETA) * h2u * t**2 - t * hpp(s) / (2 * M)
    dtt = (ONE - BETA) * h2u * s**2 - s * hpp(t) / (2 * M)
    dst = (
        (ONE - BETA) * (h2u * s * t + h1u)
        - (hp(s) + hp(t)) / (2 * M)
    )
    dsss = (ONE - BETA) * h3u * t**3 - t * hppp(s) / (2 * M)
    dttt = (ONE - BETA) * h3u * s**3 - s * hppp(t) / (2 * M)
    dsst = (
        (ONE - BETA) * (h3u * s * t**2 + 2 * h2u * t)
        - hpp(s) / (2 * M)
    )
    dstt = (
        (ONE - BETA) * (h3u * t * s**2 + 2 * h2u * s)
        - hpp(t) / (2 * M)
    )
    return value, ds, dt, dss, dst, dtt, dsss, dsst, dstt, dttt


def k_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    value, ps, pt, pss, pst, ptt, psst, pstt = pi_derivatives(s, t)
    hv = h_interior(value)
    h1v, h2v, h3v = hp(value), hpp(value), hppp(value)
    ds = h1v * ps
    dt = h1v * pt
    dss = h2v * ps**2 + h1v * pss
    dtt = h2v * pt**2 + h1v * ptt
    dst = h2v * ps * pt + h1v * pst
    dsss = h3v * ps**3 + 3 * h2v * ps * pss
    dttt = h3v * pt**3 + 3 * h2v * pt * ptt
    dsst = (
        h3v * ps**2 * pt
        + h2v * (2 * ps * pst + pt * pss)
        + h1v * psst
    )
    dstt = (
        h3v * pt**2 * ps
        + h2v * (2 * pt * pst + ps * ptt)
        + h1v * pstt
    )
    return hv, ds, dt, dss, dst, dtt, dsss, dsst, dstt, dttt


def r_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    p2 = p2_derivatives(s, t)
    kernel = k_derivatives(s, t)
    return tuple(p2[index] + BETA * kernel[index] for index in range(10))


def qdiag_derivatives(s: arb) -> tuple[arb, arb, arb, arb]:
    q = s**2 * (ONE + (ONE - s) ** 2)
    q1 = 4 * s - 6 * s**2 + 4 * s**3
    q2 = 4 - 12 * s + 12 * s**2
    q3 = -12 + 24 * s
    return (
        h_interior(q),
        hp(q) * q1,
        hpp(q) * q1**2 + hp(q) * q2,
        hppp(q) * q1**3 + 3 * hpp(q) * q1 * q2 + hp(q) * q3,
    )


def lam_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    p2 = p2_derivatives(s, t)
    gs = qdiag_derivatives(s)
    gt = qdiag_derivatives(t)
    return (
        2 * p2[0] + BETA * (gs[0] + gt[0]),
        2 * p2[1] + BETA * gs[1],
        2 * p2[2] + BETA * gt[1],
        2 * p2[3] + BETA * gs[2],
        2 * p2[4],
        2 * p2[5] + BETA * gt[2],
        2 * p2[6] + BETA * gs[3],
        2 * p2[7],
        2 * p2[8],
        2 * p2[9] + BETA * gt[3],
    )


def phi_derivatives(s: arb) -> tuple[arb, arb, arb, arb]:
    entropy = qdiag_derivatives(s)
    w = BETA * entropy[0]
    w1 = BETA * entropy[1]
    w2 = BETA * entropy[2]
    w3 = BETA * entropy[3]
    root = w.sqrt()
    first = w1 / (2 * root)
    second = w2 / (2 * root) - w1**2 / (4 * root**3)
    third = (
        w3 / (2 * root)
        - 3 * w1 * w2 / (4 * root**3)
        + 3 * w1**3 / (8 * root**5)
    )
    return root, first, second, third


def a2_derivatives(s: arb, t: arb) -> tuple[arb, ...]:
    p2 = p2_derivatives(s, t)
    ps = phi_derivatives(s)
    pt = phi_derivatives(t)
    return (
        p2[0] + ps[0] * pt[0],
        p2[1] + ps[1] * pt[0],
        p2[2] + ps[0] * pt[1],
        p2[3] + ps[2] * pt[0],
        p2[4] + ps[1] * pt[1],
        p2[5] + ps[0] * pt[2],
        p2[6] + ps[3] * pt[0],
        p2[7] + ps[2] * pt[1],
        p2[8] + ps[1] * pt[2],
        p2[9] + ps[0] * pt[3],
    )


def maximum_abs_upper(values: Iterable[arb]) -> arb:
    result = ZERO
    for value in values:
        upper = value.abs_upper()
        if upper > result:
            result = upper
        elif not result >= upper:
            result = result.union(upper).upper()
    return result.upper()


def maximum_abs_lower(values: Iterable[arb]) -> arb:
    result = ZERO
    for value in values:
        lower = value.abs_lower()
        if lower > result:
            result = lower
    return result.lower()


def taylor_certificate(
    derivatives: Callable[[arb, arb], tuple[arb, ...]],
    half_width: Fraction,
    derivative_box_radius: arb = TAYLOR_BOX_RADIUS,
) -> TaylorCertificate:
    at_zero = derivatives(XSTAR, XSTAR)
    if not (at_zero[0].contains(0) and at_zero[1].contains(0)
            and at_zero[2].contains(0)):
        raise AssertionError("interior zero/stationarity enclosure misses zero")

    lambda_s = at_zero[3].lower() - at_zero[4].abs_upper()
    lambda_t = at_zero[5].lower() - at_zero[4].abs_upper()
    lambda_lower = min_enclosure(lambda_s, lambda_t).lower()
    if not lambda_lower > 0:
        raise AssertionError("Hessian lower eigenvalue is not positive")

    derivative_box = hull(
        XSTAR.lower() - derivative_box_radius,
        XSTAR.upper() + derivative_box_radius,
    )
    on_box = derivatives(derivative_box, derivative_box)
    third = tuple(on_box[6:10])
    if not all(value.is_finite() for value in third):
        raise AssertionError("third derivative enclosure is not finite")
    max_third = maximum_abs_upper(third)
    c3 = (2 * SQRT2 * max_third).upper()
    rho = (3 * lambda_lower / c3)
    if not (rho > 0 and rho.upper() < derivative_box_radius):
        raise AssertionError("Taylor radius does not lie in derivative box")

    center_ball = arb_fraction(CENTER)
    center_error = (XSTAR - center_ball).abs_upper()
    distance = (SQRT2 * (arb_fraction(half_width) + center_error)).upper()
    coefficient = (
        lambda_lower / 2 - c3 * distance / 6
    ).lower()
    if not (distance < rho.lower() and coefficient > 0):
        raise AssertionError("chosen local square does not fit positive Taylor ball")

    # Any valid scalar tensor-norm remainder constant must in particular
    # dominate each pure coordinate third derivative (take d=e_s or e_t).
    point_pure_third = (at_zero[6], at_zero[9])
    point_c3_required = maximum_abs_lower(point_pure_third).lower()
    return TaylorCertificate(
        value=at_zero[0],
        grad_s=at_zero[1],
        grad_t=at_zero[2],
        h_ss=at_zero[3],
        h_st=at_zero[4],
        h_tt=at_zero[5],
        lambda_lower=lambda_lower,
        third_partial_enclosures=third,
        max_third_upper=max_third,
        c3_upper=c3,
        rho_formula=rho,
        derivative_box_radius=derivative_box_radius,
        chosen_half_width=half_width,
        chosen_distance_upper=distance,
        coefficient_lower=coefficient,
        point_c3_required_lower=point_c3_required,
    )


def split_box_dimension(box: Box, dimension: str) -> tuple[Box, Box]:
    depth = box.depth + 1
    if dimension == "s":
        mid = (box.s_lo + box.s_hi) / 2
        return (
            Box(box.s_lo, mid, box.t_lo, box.t_hi, depth),
            Box(mid, box.s_hi, box.t_lo, box.t_hi, depth),
        )
    if dimension == "t":
        mid = (box.t_lo + box.t_hi) / 2
        return (
            Box(box.s_lo, box.s_hi, box.t_lo, mid, depth),
            Box(box.s_lo, box.s_hi, mid, box.t_hi, depth),
        )
    raise ValueError(dimension)


def split_box(box: Box) -> tuple[Box, Box]:
    dimension = (
        "s" if box.s_hi - box.s_lo >= box.t_hi - box.t_lo else "t"
    )
    return split_box_dimension(box, dimension)


def split_cover_box(box: Box) -> tuple[Box, Box]:
    """Isolate the logarithmic near-face region without squaring thin strips."""
    product_lo = box.s_lo * box.t_lo
    product_hi = box.s_hi * box.t_hi
    lower_coefficient = (
        C0 * (ONE / arb_fraction(product_lo)).log() - ONE / M
    )
    upper_coefficient = (
        C0 * (ONE / arb_fraction(product_hi)).log() - ONE / M
    )
    if lower_coefficient > 0 and not upper_coefficient > 0:
        # The coordinate with the smaller lo/hi ratio is responsible for the
        # near-face part.  Repeatedly splitting it peels off a certified thin
        # strip while leaving the other coordinate wide.
        if box.s_lo * box.t_hi <= box.t_lo * box.s_hi:
            return split_box_dimension(box, "s")
        return split_box_dimension(box, "t")
    return split_box(box)


def cover_boxes(
    roots: Iterable[Box],
    evaluator: Callable[[Box], arb],
    max_depth: int = MAX_BB_DEPTH,
    max_cells: int = MAX_BB_CELLS,
) -> CoverStats:
    stack = list(reversed(list(roots)))
    stats = CoverStats()
    while stack:
        box = stack.pop()
        stats.processed += 1
        stats.max_depth = max(stats.max_depth, box.depth)
        if stats.processed > max_cells:
            raise RuntimeError(f"branch-and-bound exceeded {max_cells} cells")
        enclosure = evaluator(box)
        lower = enclosure.lower()
        if lower > 0:
            stats.accepted += 1
            if stats.worst_lower is None or lower < stats.worst_lower:
                stats.worst_lower = lower
                stats.worst_box = box
            elif not stats.worst_lower <= lower:
                stats.worst_lower = stats.worst_lower.union(lower).lower()
            continue
        if box.depth >= max_depth:
            raise RuntimeError(
                "branch-and-bound reached depth cap with nonpositive lower "
                f"bound {format_arb(lower)} on {box}"
            )
        left, right = split_cover_box(box)
        stats.split += 1
        stack.append(right)
        stack.append(left)
    if stats.accepted == 0 or stats.worst_lower is None:
        raise AssertionError("empty global cover")
    if stats.processed != stats.accepted + stats.split:
        raise AssertionError("branch-and-bound accounting mismatch")
    return stats


def global_root_boxes(
    local_half: Fraction,
    corner_delta: Fraction = CORNER_DELTA,
) -> list[Box]:
    local_lo = CENTER - local_half
    local_hi = CENTER + local_half
    corner_start = Fraction(1) - corner_delta
    cuts = [EPSILON, local_lo, local_hi, corner_start, Fraction(1)]
    if cuts != sorted(cuts) or len(set(cuts)) != len(cuts):
        raise AssertionError("invalid stratum partition")
    boxes: list[Box] = []
    for s_lo, s_hi in zip(cuts[:-1], cuts[1:]):
        for t_lo, t_hi in zip(cuts[:-1], cuts[1:]):
            inside_local = (
                s_lo >= local_lo and s_hi <= local_hi
                and t_lo >= local_lo and t_hi <= local_hi
            )
            inside_corner = (
                s_lo >= corner_start and t_lo >= corner_start
            )
            if inside_local or inside_corner:
                continue
            boxes.append(Box(s_lo, s_hi, t_lo, t_hi))
    if len(boxes) != 14:
        raise AssertionError(f"unexpected global partition size: {len(boxes)}")
    return boxes


def nested_zero_cell(
    evaluator: Callable[[Box], arb],
    initial: Box,
    zero_kind: str,
    depth: int,
) -> tuple[arb, Box]:
    box = initial
    enclosure = evaluator(box)
    if not enclosure.contains(0):
        raise AssertionError("mutant start box does not contain the known zero")
    for _ in range(depth):
        left, right = split_box(box)
        if zero_kind == "xstar":
            if left.s_lo != right.s_lo:  # split in s
                split = left.s_hi
                if XSTAR.upper() < arb_fraction(split):
                    box = left
                elif XSTAR.lower() > arb_fraction(split):
                    box = right
                else:
                    raise AssertionError("root bracket straddles mutant split")
            else:  # split in t
                split = left.t_hi
                if XSTAR.upper() < arb_fraction(split):
                    box = left
                elif XSTAR.lower() > arb_fraction(split):
                    box = right
                else:
                    raise AssertionError("root bracket straddles mutant split")
        elif zero_kind == "corner11":
            box = right
        else:
            raise ValueError(zero_kind)
        enclosure = evaluator(box)
        if not enclosure.contains(0):
            raise AssertionError("nested zero cell lost the exact zero")
    if enclosure.lower() > 0:
        raise AssertionError("zero-containing mutant cell was incorrectly discharged")
    return enclosure, box


def identity_checks() -> dict[str, object]:
    maximum_factorization_residual = ZERO
    maximum_lam_residual = ZERO
    maximum_a2_residual = ZERO
    maximum_face_residual = ZERO
    count = 0
    grid = [Fraction(index, 10) for index in range(1, 10)]
    for s in grid:
        for t in grid:
            factorization = r_point(s, t) - arb_fraction(s * t) * rcheck_point(s, t)
            lam_identity = (
                lam_point(s, t)
                - 2 * r_point(s, t)
                - BETA * (k_point(s, s) + k_point(t, t) - 2 * k_point(s, t))
            )
            a2_identity = (
                lam_point(s, t)
                - 2 * a2_point(s, t)
                - (phi_point(s) - phi_point(t)) ** 2
            )
            a2_diagonal = a2_point(s, s) - r_point(s, s)
            if not (
                factorization.contains(0)
                and lam_identity.contains(0)
                and a2_identity.contains(0)
                and a2_diagonal.contains(0)
            ):
                raise AssertionError("kernel identity enclosure misses zero")
            maximum_factorization_residual = max(
                maximum_factorization_residual,
                factorization.abs_upper(),
            )
            maximum_lam_residual = max(
                maximum_lam_residual,
                lam_identity.abs_upper(),
            )
            maximum_a2_residual = max(
                maximum_a2_residual,
                a2_identity.abs_upper(),
                a2_diagonal.abs_upper(),
            )
            count += 1

    for u in [Fraction(index, 16) for index in range(17)]:
        face_residuals = (
            r_point(Fraction(1), u) - C1 * h_point_fraction(u),
            r_point(Fraction(0), u),
            lam_point(Fraction(0), u) - BETA * k_point(u, u),
            lam_point(Fraction(1), u)
            - 2 * C0 * h_point_fraction(u) - BETA * k_point(u, u),
            a2_point(Fraction(0), u),
            a2_point(Fraction(1), u) - C0 * h_point_fraction(u),
        )
        for residual in face_residuals:
            if not residual.contains(0):
                raise AssertionError("face identity enclosure misses zero")
            maximum_face_residual = max(
                maximum_face_residual,
                residual.abs_upper(),
            )

    quartic_protocol_residual = pi_arb(XSTAR, XSTAR) + XSTAR**2 - ONE
    r_at_zero = r_derivatives(XSTAR, XSTAR)
    lam_at_zero = lam_derivatives(XSTAR, XSTAR)
    a2_at_zero = a2_derivatives(XSTAR, XSTAR)
    if not (
        quartic_protocol_residual.contains(0)
        and r_at_zero[0].contains(0)
        and r_at_zero[1].contains(0)
        and r_at_zero[2].contains(0)
        and lam_at_zero[0].contains(0)
        and lam_at_zero[1].contains(0)
        and lam_at_zero[2].contains(0)
        and a2_at_zero[0].contains(0)
        and a2_at_zero[1].contains(0)
        and a2_at_zero[2].contains(0)
    ):
        raise AssertionError("exact xstar zero/stationarity enclosure misses zero")
    return {
        "grid_point_count": count,
        "r_equals_st_rcheck_max_residual_enclosure": format_arb(
            maximum_factorization_residual
        ),
        "lam_decomposition_max_residual_enclosure": format_arb(
            maximum_lam_residual
        ),
        "a2_rank_one_identity_max_residual_enclosure": format_arb(
            maximum_a2_residual
        ),
        "face_identity_max_residual_enclosure": format_arb(
            maximum_face_residual
        ),
        "xstar_quartic_protocol_residual_enclosure": format_arb(
            quartic_protocol_residual
        ),
        "xstar_exact_zero_derivation_lemma_a": (
            "The quartic is pi(x,x)+x^2-1=0, hence h(pi(x,x))=h(x^2). "
            "With p=h(x)/h(x^2) and m=p*x this gives R(x,x)=0. "
            "The definition of beta sets the diagonal derivative to zero; "
            "symmetry then gives R_s=R_t=0."
        ),
        "xstar_exact_zero_derivation_lemma_b": (
            "Lam(x,x)=2*R(x,x)=0 and, on the diagonal, "
            "Lam_s=Lam_t=2*R_s=0."
        ),
        "xstar_exact_zero_derivation_lemma_a2": (
            "A2(s,s)=P2(s,s)+phi(s)^2=R(s,s) on the whole diagonal. "
            "Also A2_s=R_s there because phi*phi' is half the diagonal "
            "derivative of Q2.  Thus A2(x,x)=A2_s=A2_t=0."
        ),
    }


def box_record(box: Box | None) -> dict[str, str] | None:
    if box is None:
        return None
    return {
        "s_lo": format_fraction(box.s_lo),
        "s_hi": format_fraction(box.s_hi),
        "t_lo": format_fraction(box.t_lo),
        "t_hi": format_fraction(box.t_hi),
        "depth": str(box.depth),
    }


def cover_record(stats: CoverStats) -> dict[str, object]:
    return {
        "processed_cells": stats.processed,
        "accepted_cells": stats.accepted,
        "split_cells": stats.split,
        "max_depth": stats.max_depth,
        "worst_certified_lower_enclosure": format_arb(stats.worst_lower),
        "worst_box": box_record(stats.worst_box),
    }


def taylor_record(certificate: TaylorCertificate) -> dict[str, object]:
    return {
        "value_enclosure": format_arb(certificate.value),
        "gradient_s_enclosure": format_arb(certificate.grad_s),
        "gradient_t_enclosure": format_arb(certificate.grad_t),
        "hessian_ss_enclosure": format_arb(certificate.h_ss),
        "hessian_st_enclosure": format_arb(certificate.h_st),
        "hessian_tt_enclosure": format_arb(certificate.h_tt),
        "lambda_min_lower_enclosure": format_arb(certificate.lambda_lower),
        "third_partial_enclosures": [
            format_arb(value) for value in certificate.third_partial_enclosures
        ],
        "max_third_partial_upper_enclosure": format_arb(
            certificate.max_third_upper
        ),
        "c3_upper_enclosure": format_arb(certificate.c3_upper),
        "rho_equals_3lambda_over_c3_enclosure": format_arb(
            certificate.rho_formula
        ),
        "derivative_box_radius_enclosure": format_arb(
            certificate.derivative_box_radius
        ),
        "proof_center": (
            "the certified Arb isolating ball for the unique quartic root xstar"
        ),
        "partition_square_center": format_fraction(CENTER),
        "partition_center_to_xstar_distance_upper_enclosure": format_arb(
            (XSTAR - arb_fraction(CENTER)).abs_upper()
        ),
        "distance_bound_direction": (
            "sqrt(2)*(square half-width + sup|xstar-partition center|)"
        ),
        "chosen_square_half_width": format_fraction(
            certificate.chosen_half_width
        ),
        "chosen_euclidean_distance_upper_enclosure": format_arb(
            certificate.chosen_distance_upper
        ),
        "quadratic_minus_remainder_coefficient_lower_enclosure": format_arb(
            certificate.coefficient_lower
        ),
    }


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def build_report() -> dict[str, object]:
    identities = identity_checks()
    binding_crosscheck = binding_parameter_crosscheck()

    log_epsilon = (ONE / arb_fraction(EPSILON)).log()
    strip_margin = (C0 * log_epsilon - ONE / M).lower()
    if not strip_margin > 0:
        raise AssertionError("near-zero strip lower bound is not positive")

    delta = arb_fraction(CORNER_DELTA)
    log_corner = (ONE / (2 * delta)).log()
    corner_a_coefficient = C1 - delta * (ONE + BETA) / 2
    corner_a_margin = (
        corner_a_coefficient * log_corner - (LOG2 + ONE) / (2 * M)
    ).lower()
    corner_b_coefficient = 2 - ONE / M - (ONE + BETA) * delta
    corner_b_margin = (
        corner_b_coefficient * log_corner - (LOG2 + ONE) / M
    ).lower()
    corner_b_entropy_coefficient = BETA * (2 - 2 * delta) - ONE / M
    if not (
        corner_a_coefficient > 0
        and corner_a_margin > 0
        and corner_b_coefficient > 0
        and corner_b_margin > 0
        and corner_b_entropy_coefficient < 0
    ):
        raise AssertionError("(1,1) corner lower bound failed")

    delta_a2 = arb_fraction(CORNER_DELTA_A2)
    log_corner_a2 = (ONE / (2 * delta_a2)).log()
    corner_a2_coefficient = (
        C0 - (ONE - BETA) * delta_a2 / 2
    )
    corner_a2_margin = (
        corner_a2_coefficient * log_corner_a2
        - (LOG2 + ONE) / (2 * M)
    ).lower()
    if not (corner_a2_coefficient > 0 and corner_a2_margin > 0):
        raise AssertionError("A2 (1,1) corner lower bound failed")

    taylor_a = taylor_certificate(r_derivatives, LOCAL_HALF_A)
    taylor_b = taylor_certificate(lam_derivatives, LOCAL_HALF_B)
    taylor_a2 = taylor_certificate(
        a2_derivatives,
        LOCAL_HALF_A2,
        TAYLOR_BOX_RADIUS_A2,
    )

    global_a = cover_boxes(global_root_boxes(LOCAL_HALF_A), r_global_lower)
    global_b = cover_boxes(global_root_boxes(LOCAL_HALF_B), lam_global_lower)
    global_a2 = cover_boxes(
        global_root_boxes(LOCAL_HALF_A2, CORNER_DELTA_A2),
        a2_global_lower,
        max_depth=MAX_BB_DEPTH_A2,
    )

    local_a_box = Box(
        CENTER - LOCAL_HALF_A, CENTER + LOCAL_HALF_A,
        CENTER - LOCAL_HALF_A, CENTER + LOCAL_HALF_A,
    )
    local_a_mutant_range, local_a_mutant_box = nested_zero_cell(
        r_range, local_a_box, "xstar", 14
    )
    corner_b_box = Box(
        Fraction(1) - CORNER_DELTA, Fraction(1),
        Fraction(1) - CORNER_DELTA, Fraction(1),
    )
    corner_b_mutant_range, corner_b_mutant_box = nested_zero_cell(
        lam_range, corner_b_box, "corner11", 14
    )
    local_a2_box = Box(
        CENTER - LOCAL_HALF_A2,
        CENTER + LOCAL_HALF_A2,
        CENTER - LOCAL_HALF_A2,
        CENTER + LOCAL_HALF_A2,
    )
    local_a2_mutant_range, local_a2_mutant_box = nested_zero_cell(
        a2_range, local_a2_box, "xstar", 14
    )
    corner_a2_box = Box(
        Fraction(1) - CORNER_DELTA_A2,
        Fraction(1),
        Fraction(1) - CORNER_DELTA_A2,
        Fraction(1),
    )
    corner_a2_mutant_range, corner_a2_mutant_box = nested_zero_cell(
        a2_range, corner_a2_box, "corner11", 14
    )

    mutant_strip_margin = (-C0 * log_epsilon - ONE / M).upper()
    if not mutant_strip_margin < 0:
        raise AssertionError("flipped strip-sign mutant unexpectedly passed")

    mutant_c3_a = (taylor_a.c3_upper / 16).upper()
    mutant_c3_b = (taylor_b.c3_upper / 16).upper()
    mutant_c3_a2 = (taylor_a2.c3_upper / 16).upper()
    if not (
        mutant_c3_a < taylor_a.point_c3_required_lower
        and mutant_c3_b < taylor_b.point_c3_required_lower
        and mutant_c3_a2 < taylor_a2.point_c3_required_lower
    ):
        raise AssertionError("understated C3 mutant was not refuted at xstar")

    tool_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report: dict[str, object] = {
        "artifact": "liu9-h2-twovar",
        "claim_status": "CERTIFIED_PROVED",
        "precision_bits": PRECISION_BITS,
        "tool": "liu9_h2_twovar_lemmas.py (python-flint Arb)",
        "tool_sha256": tool_hash,
        "constants": {
            "xstar_defining_polynomial": "x^4-2*x^3+3*x^2-1",
            "xstar_enclosure": format_arb(XSTAR),
            "xstar_root_lo": format_fraction(PARAMETERS.root_lo),
            "xstar_root_hi": format_fraction(PARAMETERS.root_hi),
            "xstar_isolating_interval_width_enclosure": format_arb(
                arb_fraction(PARAMETERS.root_hi - PARAMETERS.root_lo)
            ),
            "root_polynomial_lo_enclosure": format_arb(
                PARAMETERS.root_polynomial_lo
            ),
            "root_polynomial_hi_enclosure": format_arb(
                PARAMETERS.root_polynomial_hi
            ),
            "root_uniqueness": (
                "f'(z)=2z(2(z-3/4)^2+15/8)>0 on (0,1]; "
                "the certified endpoint signs therefore isolate the unique root"
            ),
            "p_enclosure": format_arb(P),
            "m_closed_form": "xstar*h(xstar)/h(xstar^2)",
            "m_closed_form_derivation": (
                "p=h(xstar)/h(xstar^2) and m=p*xstar"
            ),
            "m_enclosure": format_arb(M),
            "beta_enclosure": format_arb(BETA),
            "c0_equals_1_minus_beta_minus_1_over_2m_enclosure": format_arb(C0),
            "c1_equals_1_minus_1_over_2m_enclosure": format_arb(C1),
        },
        "binding_parameter_crosscheck": binding_crosscheck,
        "identities": identities,
        "lemma_a": {
            "claim": "R(s,t)>=0 on [0,1]^2",
            "status": "PROVED",
            "faces": {
                "s_or_t_zero": "R(0,t)=R(s,0)=0 exactly",
                "s_or_t_one": "R(1,t)=R(t,1)=c1*h(t)>=0 with certified c1>0",
            },
            "near_zero_strip": {
                "epsilon": format_fraction(EPSILON),
                "log_1_over_epsilon_enclosure": format_arb(log_epsilon),
                "lower_coefficient_enclosure": format_arb(strip_margin),
                "bound": (
                    "If min(s,t)<=epsilon, R>=s*t*(c0*log(1/(s*t))-1/m)"
                    ">=s*t*lower_coefficient>=0; endpoint faces are exact zeros"
                ),
            },
            "corner_1_1": {
                "delta": format_fraction(CORNER_DELTA),
                "log_1_over_2delta_enclosure": format_arb(log_corner),
                "log_coefficient_enclosure": format_arb(corner_a_coefficient),
                "lower_coefficient_enclosure": format_arb(corner_a_margin),
                "derivation": (
                    "Let W=1-s*t=u-a*b and V=1-pi(s,t)=W-a*b*s*t. "
                    "Then W>=u*(1-delta/2), V>=u*(1-delta), and both "
                    "log(1/W),log(1/V)>=log(1/u).  Also, with "
                    "A=a*log(1/a)+b*log(1/b), A<=u*(log(1/u)+log(2)) "
                    "and (1-b)h(a)+(1-a)h(b)<=A+u.  Applying h(z)>="
                    "z*log(1/z) to the positive terms and h(z)<="
                    "z*(log(1/z)+1) to the subtracted terms gives the bound."
                ),
                "bound": (
                    "For a=1-s,b=1-t,u=a+b in (0,2delta], "
                    "R>=u*(log_coefficient*log(1/u)-(log(2)+1)/(2m))"
                    ">=u*lower_coefficient>=0; u=0 is the exact corner zero"
                ),
            },
            "local_xstar": taylor_record(taylor_a),
            "global_branch_and_bound": cover_record(global_a),
            "mutations": [
                {
                    "name": "drop_xstar_local_stratum",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_local_coefficient_lower_enclosure": format_arb(
                        taylor_a.coefficient_lower
                    ),
                    "mutant_depth": local_a_mutant_box.depth,
                    "mutant_zero_cell_enclosure": format_arb(local_a_mutant_range),
                    "mutant_zero_cell": box_record(local_a_mutant_box),
                    "reason": "the nested cell still contains R(x*,x*)=0 and cannot have a strictly positive lower bound",
                },
                {
                    "name": "flip_near_zero_log_term_sign",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_margin_enclosure": format_arb(strip_margin),
                    "mutant_margin_enclosure": format_arb(mutant_strip_margin),
                    "reason": "the flipped-sign lower coefficient is certified negative",
                },
                {
                    "name": "understate_c3_by_factor_16",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_c3_upper_enclosure": format_arb(taylor_a.c3_upper),
                    "mutant_c3_enclosure": format_arb(mutant_c3_a),
                    "pointwise_required_c3_lower_enclosure": format_arb(
                        taylor_a.point_c3_required_lower
                    ),
                    "reason": "the mutant is below a certified pointwise third-derivative requirement at xstar",
                },
            ],
        },
        "lemma_b": {
            "claim": "Lam(s,t)>=0 on [0,1]^2",
            "status": "PROVED",
            "faces": {
                "s_or_t_zero": "Lam(0,t)=Lam(t,0)=beta*h(pi(t,t))>=0",
                "s_or_t_one": "Lam(1,t)=Lam(t,1)=2*c0*h(t)+beta*h(pi(t,t))>=0",
                "corner_zeros": "the four corners are exact zeros of the displayed face formulas",
            },
            "near_zero_strip": {
                "epsilon": format_fraction(EPSILON),
                "log_1_over_epsilon_enclosure": format_arb(log_epsilon),
                "lower_coefficient_enclosure": format_arb(2 * strip_margin),
                "bound": (
                    "If min(s,t)<=epsilon, Lam>=2*s*t*"
                    "(c0*log(1/(s*t))-1/m)>=s*t*lower_coefficient>=0"
                ),
            },
            "corner_1_1": {
                "delta": format_fraction(CORNER_DELTA),
                "log_1_over_2delta_enclosure": format_arb(log_corner),
                "log_coefficient_enclosure": format_arb(corner_b_coefficient),
                "lower_coefficient_enclosure": format_arb(corner_b_margin),
                "auxiliary_entropy_coefficient_enclosure": format_arb(
                    corner_b_entropy_coefficient
                ),
                "derivation": (
                    "Let W=1-s*t=u-a*b and A=a*log(1/a)+b*log(1/b). "
                    "Then W>=u*(1-delta/2), A<=u*(log(1/u)+log(2)), "
                    "and the subtracted weighted entropies are at most A+u. "
                    "For U_a=1-pi(1-a,1-a), "
                    "(2-2delta)*a<=U_a<=2a, so h(U_a)>="
                    "(2-2delta)*a*(log(1/a)-log(2)), and similarly for b. "
                    "The certified auxiliary entropy coefficient is negative, "
                    "so substituting the upper bound on A is the correct lower-"
                    "bound direction and simplifies to the displayed bound."
                ),
                "bound": (
                    "For a=1-s,b=1-t,u=a+b in (0,2delta], "
                    "Lam>=u*(log_coefficient*log(1/u)-(log(2)+1)/m)"
                    ">=u*lower_coefficient>=0; u=0 is the exact corner zero"
                ),
            },
            "local_xstar": taylor_record(taylor_b),
            "global_branch_and_bound": cover_record(global_b),
            "mutations": [
                {
                    "name": "drop_corner_1_1_stratum",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_corner_coefficient_lower_enclosure": format_arb(
                        corner_b_margin
                    ),
                    "mutant_depth": corner_b_mutant_box.depth,
                    "mutant_zero_cell_enclosure": format_arb(corner_b_mutant_range),
                    "mutant_zero_cell": box_record(corner_b_mutant_box),
                    "reason": "the nested cell still contains Lam(1,1)=0 and cannot have a strictly positive lower bound",
                },
                {
                    "name": "understate_c3_by_factor_16",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_c3_upper_enclosure": format_arb(taylor_b.c3_upper),
                    "mutant_c3_enclosure": format_arb(mutant_c3_b),
                    "pointwise_required_c3_lower_enclosure": format_arb(
                        taylor_b.point_c3_required_lower
                    ),
                    "reason": "the mutant is below a certified pointwise third-derivative requirement at xstar",
                },
                {
                    "name": "flip_near_zero_log_term_sign",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_margin_enclosure": format_arb(2 * strip_margin),
                    "mutant_margin_enclosure": format_arb(2 * mutant_strip_margin),
                    "reason": "the flipped-sign lower coefficient is certified negative",
                },
            ],
        },
        "lemma_a2": {
            "claim": "A2(s,t)=P2(s,t)+phi(s)*phi(t)>=0 on [0,1]^2",
            "status": "PROVED",
            "definition": "phi(s)=sqrt(beta*h(pi(s,s)))",
            "faces": {
                "s_or_t_zero": "A2(0,t)=A2(t,0)=0 exactly",
                "s_or_t_one": "A2(1,t)=A2(t,1)=c0*h(t)>=0 with certified c0>0",
            },
            "near_zero_strip": {
                "epsilon": format_fraction(EPSILON),
                "log_1_over_epsilon_enclosure": format_arb(log_epsilon),
                "lower_coefficient_enclosure": format_arb(strip_margin),
                "bound": (
                    "Dropping phi(s)*phi(t)>=0 gives A2>=P2>="
                    "s*t*(c0*log(1/(s*t))-1/m).  If min(s,t)<=epsilon "
                    "this is >=s*t*lower_coefficient>=0; endpoint faces "
                    "are exact zeros."
                ),
            },
            "corner_1_1": {
                "delta": format_fraction(CORNER_DELTA_A2),
                "log_1_over_2delta_enclosure": format_arb(log_corner_a2),
                "log_coefficient_enclosure": format_arb(
                    corner_a2_coefficient
                ),
                "lower_coefficient_enclosure": format_arb(corner_a2_margin),
                "derivation": (
                    "For a=1-s,b=1-t,u=a+b and W=1-s*t=u-a*b, "
                    "h(W)>=u*(1-delta/2)*log(1/u).  With "
                    "A=a*log(1/a)+b*log(1/b)<=u*(log(1/u)+log(2)), "
                    "the subtracted weighted entropies are at most A+u. "
                    "Dropping phi(s)*phi(t)>=0 therefore gives the bound."
                ),
                "bound": (
                    "For u in (0,2delta], A2>=P2>=u*"
                    "(log_coefficient*log(1/u)-(log(2)+1)/(2m))"
                    ">=u*lower_coefficient>=0; u=0 is the exact corner zero."
                ),
            },
            "local_xstar": taylor_record(taylor_a2),
            "global_branch_and_bound": cover_record(global_a2),
            "mutations": [
                {
                    "name": "drop_xstar_local_stratum",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_local_coefficient_lower_enclosure": format_arb(
                        taylor_a2.coefficient_lower
                    ),
                    "mutant_depth": local_a2_mutant_box.depth,
                    "mutant_zero_cell_enclosure": format_arb(
                        local_a2_mutant_range
                    ),
                    "mutant_zero_cell": box_record(local_a2_mutant_box),
                    "reason": "the nested cell contains A2(x*,x*)=0 and cannot have a strictly positive lower bound",
                },
                {
                    "name": "drop_corner_1_1_stratum",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_corner_coefficient_lower_enclosure": format_arb(
                        corner_a2_margin
                    ),
                    "mutant_depth": corner_a2_mutant_box.depth,
                    "mutant_zero_cell_enclosure": format_arb(
                        corner_a2_mutant_range
                    ),
                    "mutant_zero_cell": box_record(corner_a2_mutant_box),
                    "reason": "the nested cell contains A2(1,1)=0 and cannot have a strictly positive lower bound",
                },
                {
                    "name": "understate_c3_by_factor_16",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_c3_upper_enclosure": format_arb(taylor_a2.c3_upper),
                    "mutant_c3_enclosure": format_arb(mutant_c3_a2),
                    "pointwise_required_c3_lower_enclosure": format_arb(
                        taylor_a2.point_c3_required_lower
                    ),
                    "reason": "the mutant is below a certified pure-coordinate third derivative at xstar",
                },
                {
                    "name": "flip_near_zero_log_term_sign",
                    "status": "FAILED_AS_REQUIRED",
                    "sound_margin_enclosure": format_arb(strip_margin),
                    "mutant_margin_enclosure": format_arb(mutant_strip_margin),
                    "reason": "the flipped-sign lower coefficient is certified negative",
                },
            ],
        },
        "bound_directions": {
            "entropy_lower": "h(u)>=-u*log(u); the omitted term -(1-u)log(1-u) is nonnegative",
            "entropy_upper": "h(u)<=-u*log(u)+u; this is equivalent to -(1-u)log(1-u)<=u",
            "strip_lower": "all positive entropy terms not displayed are dropped, and h(s),h(t) in the subtracted term use upper bounds",
            "corner_lower": "positive h terms use lower bounds; the subtracted weighted entropies use h(a)+h(b)<=a*log(1/a)+b*log(1/b)+a+b and the coefficient multiplying their upper bound is negative",
            "hessian_lower": "lambda_min is bounded from below by min(Hss-|Hst|,Htt-|Hst|), with diagonal lower endpoints and the off-diagonal absolute upper endpoint",
            "third_derivative_upper": "C3=2*sqrt(2)*max|partial_ijk| uses Arb absolute upper endpoints on each reported whole derivative box (radius 0.03 for R/Lam, 0.025 for A2)",
            "taylor_lower": "Remainder magnitude is at most C3*||d||_2^3/6, hence f>=||d||_2^2*(lambda_min/2-C3*||d||_2/6)",
            "local_root_center": (
                "All local values, gradients, Hessians, and derivative boxes are "
                "evaluated at or around the Arb ball bounded by the certified "
                "quartic isolating endpoints, never a decimal approximation.  "
                "The rational partition center is derived just left of that "
                "bracket; its reported displacement upper bound is added to the "
                "square half-width before the Taylor distance is bounded above."
            ),
            "global_lower": (
                "Each accepted box uses a certified lower bound: either the "
                "analytic logarithmic-product bound, the lower endpoint of the "
                "raw Arb range, or f(mid)-sup|f_s|*halfwidth_s-"
                "sup|f_t|*halfwidth_t.  The strongest certified bound is used; "
                "entropy ranges use endpoint values and log(2) exactly when the "
                "argument interval crosses 1/2."
            ),
            "protocol_range": "pi is coordinatewise nondecreasing because pi_s>=t^2 and pi_t>=s^2; q(s)=pi(s,s) is increasing since q'=2s(2-3s+2s^2)>0",
            "square_root_range": (
                "beta*h(pi(s,s)) is nonnegative.  Since sqrt is increasing, "
                "phi interval lower=sqrt(max(0,radicand lower)) and upper="
                "sqrt(radicand upper); products use nonnegative interval factors. "
                "A2 corner and near-zero lower bounds only drop this positive term."
            ),
            "a2_gradient_lower": (
                "On strict interior boxes phi>0 and phi'=beta*q'(s)*"
                "h'(q(s))/(2*phi).  Arb encloses this derivative on the whole "
                "box before applying the centered first-order lower bound."
            ),
            "parameter_crosscheck": (
                "The stored binary mpf values returned by "
                "liu9_binding.solve_equation_parameters(100) are converted to "
                "exact rationals and then to Arb.  Agreement means the Arb "
                "enclosures of isolated-root minus binding x and closed-form "
                "mean minus binding mean both contain zero."
            ),
        },
        "coverage": (
            "Near-zero strips cover min(s,t)<=epsilon.  On [epsilon,1]^2, "
            "each lemma's local xstar square and (1,1) corner square are "
            "certified analytically (corner delta 1/10000 for R/Lam and "
            "1/100000000 for A2).  For each lemma, the sorted distinct cut "
            "list gives 16 rectangles; removing exactly the local and corner "
            "rectangles leaves 14 exact rational rectangles exhausted by "
            "deterministic Arb branch-and-bound."
        ),
        "report_sha256_scope": "SHA-256 of canonical sorted-key compact JSON with report_sha256 omitted",
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    return report


def write_report(report: dict[str, object], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(report, sort_keys=True, indent=1) + "\n"
    output.write_text(serialized, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    default_output = (
        Path(__file__).resolve().parent
        / "verification" / "results" / "liu9-h2-twovar.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=default_output)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report()
    write_report(report, args.output)
    print("LEMMA_A PROVED")
    print("LEMMA_B PROVED")
    print("LEMMA_A2 PROVED")
    print(f"REPORT {args.output}")
    print(f"REPORT_SHA256 {report['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
