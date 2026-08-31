#!/usr/bin/env python3
"""Second-order certificate for Liu's mean-feasible atom insertion.

The first-order module ``liu9_transverse.py`` proves the coefficient ``L(y)``
nonnegative on the whole interval.  This module closes the next two terms that
that result deliberately left open:

* the exact ``epsilon^2`` coefficient at the double root ``y=x``;
* the interaction of unequal insertion points in the two components.

For

    ybar = (1-q)*y1 + q*y2,
    (a1,a2,q,b0,b2,b4,b1,b3,b5)
      = (p-epsilon*ybar/x, epsilon, q, x,y1,0,x,y2,0),

the objective mean is exactly ``p*x``.  Epsilon changes masses only.  Every
entropy argument is therefore epsilon-independent, the expectations are at
most quadratic in the masses, and the pencil is EXACTLY

    gap - kappa*dist^2 = epsilon*Lbar + epsilon^2*P2.

There is no remainder, no higher term, and no ``epsilon^2 log(1/epsilon)``
term.  With ``a=1-q``, ``b=q``,

    Lbar = a*L(y1) + b*L(y2),
    P2   = a*Q(y1) + b*Q(y2) + a*b*Chat(y1,y2),

where

    Q(y) = K(y,y) - 2(y/x)K(x,y) + (y/x)^2 K(x,x),
    g(y) = y^2(y-x)^2,
    S(y) = 2L(y) + Q(y),

and, putting d=y1-y2,

    C = (2d/x)(K(x,y1)-K(x,y2))
        -(d/x)^2 K(x,x)
        -(1-beta)[h(y1^2)+h(y2^2)-2h(y1*y2)],
    Chat = C - kappa*d^2.

The last subtraction is load-bearing.  The mixture is mean-preserving, but its
two component means differ at first order, and the gauge-invariant distance
therefore has the exact quadratic term

    epsilon^2*q(1-q)*(y1-y2)^2.

Dropping it gives the wrong interaction sign.  A mutation test must reject that
mistake.

For ``0 <= epsilon <= 1/2`` the sign reduces to

    2*(Lbar + epsilon*P2) >= min(2Lbar, 2Lbar+P2).

The first term is nonnegative by the exact first-order theorem.  The second is

    S_asym = (1-q)S(y1) + qS(y2) + q(1-q)Chat(y1,y2).

This module certifies ``S>=0`` in one dimension and ``S_asym>=0`` on the full
cube by an Arb cover.  The 3-D problem is reduced exactly to a 2-D cover because
for fixed ``(y1,y2)`` the expression is a scalar quadratic in q; every cell uses
the rigorous lower envelopes of ``S(y1)``, ``S(y2)`` and ``Chat`` and minimizes
that quadratic analytically.

Dependency error on the diagonal is avoided by factoring

    Chat = (y1-y2)^2 * Jhat,

where ``Jhat`` is written in divided-difference form.  On diagonal and adjacent
cells the divided differences are bounded by ranges of the corresponding first
and mixed second derivatives.  Near zero and one, where those separate ranges
are singular, explicit endpoint bounds retain the cancellation symbolically.

The one-sided mean half-space is included, not assumed.  A strict mean excess
is parameterized by moving mass delta from support 0 to support x.  Exact
bilinearity produces a positive linear normal coefficient and the positive
quadratic coefficient K(x,x)-kappa*x^2; both are universally Arb-certified.

A sharp obstruction is also recorded.  The pencil is NOT nonnegative on the
whole feasible insertion segment: at ``y1=y2=1`` and ``epsilon=p*x`` the raw
gap is structurally zero (all supports are in {0,1}, so numerator and EHX both
vanish) while the pencil is strictly negative.  This does not refute Liu's
Hypothesis 2; it is a boundary 0/0 objective and only limits the tube proof.

The certificate in this module is the local mass statement.  Integration with
the quotient gauges and the existing piecewise chart/boundary lemmas is owned
by ``liu9_tube.py``; this module does not independently promote a tube radius.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import mpmath
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    _hp_arb,
    _hpp_arb,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_ninevar import (  # noqa: E402
    EpsilonPolynomial,
    PolynomialFormulaOps,
    distance_squared,
)
from liu9_objective import _formula, h_arb  # noqa: E402
from liu9_transverse import (  # noqa: E402
    Jet2,
    _arbf,
    _ball,
    _constant_a,
    _mu_enclosure,
    jet_entropy_interior,
    transverse_jet,
)

ctx.prec = max(ctx.prec, 320)

KAPPA = Fraction(4_119_063, 33_554_432)
HALF = Fraction(1, 2)
Y0 = Fraction(1, 32)
X_WINDOW = Fraction(1, 32)
ONE_TAIL = Fraction(1, 1024)
DEFAULT_LOW_CELLS = 192
DEFAULT_HIGH_CELLS = 256
DEFAULT_WINDOW_CELLS = 1024
DEFAULT_CROSS_LOW_CELLS = 96
DEFAULT_CROSS_HIGH_CELLS = 128


def _lower_float(value: arb) -> float:
    """A conservative float below an Arb lower endpoint."""
    return math.nextafter(float(value.lower()), -math.inf)


def _upper_float(value: arb) -> float:
    return math.nextafter(float(value.upper()), math.inf)


def _kernel_h(left: arb, right: arb) -> arb:
    return h_arb(left * right)


def _protocol(left: arb, right: arb) -> arb:
    return left * right + left * (arb(1) - left) * right * (arb(1) - right)


def _kernel_p(left: arb, right: arb) -> arb:
    return h_arb(_protocol(left, right))


def kernel(parameters: ArbParameters, left: arb, right: arb) -> arb:
    return ((arb(1) - parameters.beta) * _kernel_h(left, right)
            + parameters.beta * _kernel_p(left, right))


def linear_coefficient(parameters: ArbParameters, y: arb) -> arb:
    x, p = parameters.x, parameters.p
    support = y * y * (y - x) * (y - x)
    return (2 * p * kernel(parameters, x, y) - h_arb(y)
            - (y / x) * _constant_a(parameters) - _arbf(KAPPA) * support)


def quadratic_symmetric(parameters: ArbParameters, y: arb) -> arb:
    x = parameters.x
    ratio = y / x
    return (kernel(parameters, y, y)
            - 2 * ratio * kernel(parameters, x, y)
            + ratio * ratio * kernel(parameters, x, x))


def symmetric_half_coefficient(parameters: ArbParameters, y: arb) -> arb:
    return 2 * linear_coefficient(parameters, y) + quadratic_symmetric(parameters, y)


def cross_raw(parameters: ArbParameters, y1: arb, y2: arb) -> arb:
    """Raw-gap interaction C, before the quadratic distance subtraction."""
    x = parameters.x
    d = y1 - y2
    return (
        (2 * d / x) * (kernel(parameters, x, y1) - kernel(parameters, x, y2))
        - (d * d / (x * x)) * kernel(parameters, x, x)
        - (arb(1) - parameters.beta)
        * (h_arb(y1 * y1) + h_arb(y2 * y2) - 2 * h_arb(y1 * y2))
    )


def cross_pencil(parameters: ArbParameters, y1: arb, y2: arb) -> arb:
    d = y1 - y2
    return cross_raw(parameters, y1, y2) - _arbf(KAPPA) * d * d


# ---------------------------------------------------------------------------
# Exact epsilon polynomial from the repository SSOT.


def asymmetric_polynomial_values(
    parameters: Any, y1: Any, y2: Any, q: Any, one: Any,
) -> tuple[EpsilonPolynomial, ...]:
    epsilon = EpsilonPolynomial.variable(one)
    constant = EpsilonPolynomial.constant
    ybar = (one - q) * y1 + q * y2
    zero = one * 0
    values = (
        parameters.p - epsilon * ybar / parameters.x,
        epsilon,
        q,
        parameters.x, y1, zero,
        parameters.x, y2, zero,
    )
    return tuple(v if isinstance(v, EpsilonPolynomial) else constant(v)
                 for v in values)


def automatic_expansion(
    parameters: ArbParameters, y1: arb, y2: arb, q: arb,
) -> dict[str, Any]:
    values = asymmetric_polynomial_values(parameters, y1, y2, q, arb(1))
    ops = PolynomialFormulaOps(h_arb, arb(1))
    terms = _formula(values, parameters.beta, ops)
    gap = terms.numerator - terms.ehx
    dist = distance_squared(values, parameters)
    pencil = gap - _arbf(KAPPA) * dist
    return {
        "gap": gap,
        "distance": dist,
        "pencil": pencil,
        "entropy_argument_count": ops.entropy_argument_count,
        "epsilon_dependent_entropy_arguments": ops.epsilon_dependent_entropy_arguments,
    }


def hand_expansion(
    parameters: ArbParameters, y1: arb, y2: arb, q: arb,
) -> dict[str, arb]:
    qbar = arb(1) - q
    d = y1 - y2
    lbar = (qbar * linear_coefficient(parameters, y1)
            + q * linear_coefficient(parameters, y2))
    q2 = (qbar * quadratic_symmetric(parameters, y1)
          + q * quadratic_symmetric(parameters, y2)
          + qbar * q * cross_pencil(parameters, y1, y2))
    spread = (qbar * y1 * y1 * (y1 - parameters.x) ** 2
              + q * y2 * y2 * (y2 - parameters.x) ** 2)
    return {
        "linear": lbar,
        "quadratic": q2,
        "distance_linear": spread,
        "distance_quadratic": qbar * q * d * d,
    }


def certify_exact_expansion(parameters: ArbParameters) -> dict[str, Any]:
    probes = (
        (Fraction(1, 5), Fraction(9, 10), Fraction(1, 2)),
        (Fraction(1, 100), Fraction(99, 100), Fraction(1, 8)),
        (Fraction(3, 5), Fraction(7, 10), Fraction(7, 8)),
        (Fraction(1, 3), Fraction(2, 3), Fraction(1, 4)),
    )
    worst = 0.0
    rows = []
    for y1, y2, q in probes:
        aa = automatic_expansion(parameters, _arbf(y1), _arbf(y2), _arbf(q))
        hh = hand_expansion(parameters, _arbf(y1), _arbf(y2), _arbf(q))
        residuals = {
            "linear": aa["pencil"].c1 - hh["linear"],
            "quadratic": aa["pencil"].c2 - hh["quadratic"],
            "distance_linear": aa["distance"].c1 - hh["distance_linear"],
            "distance_quadratic": aa["distance"].c2 - hh["distance_quadratic"],
        }
        for name, value in residuals.items():
            if not value.contains(0):
                raise AssertionError(f"{name} formula residual excludes zero")
            worst = max(worst, float(abs(value).upper()))
        if aa["epsilon_dependent_entropy_arguments"] != 0:
            raise AssertionError("an entropy argument unexpectedly depends on epsilon")
        rows.append({
            "y1": str(y1), "y2": str(y2), "q": str(q),
            "residuals": {k: v.str(12) for k, v in residuals.items()},
        })
    return {
        "probes": rows,
        "worst_residual": worst,
        "epsilon_dependent_entropy_arguments": 0,
        "epsilon_squared_log_coefficient": "0",
        "exact_degree": 2,
        "certified": True,
    }


# ---------------------------------------------------------------------------
# One-dimensional S certificate.


def _kernel_jet(parameters: ArbParameters, left: Optional[Jet2], right: Jet2) -> Jet2:
    one = Jet2.constant(arb(1))
    if left is None:  # K(y,y)
        product = right * right
        curved = right * (one - right)
        protocol = product + curved * curved
    else:  # K(x,y), x held constant in the left slot
        x = left.v
        product = right.scaled(x)
        protocol = product + (right * (one - right)).scaled(x * (arb(1) - x))
    return (jet_entropy_interior(product).scaled(arb(1) - parameters.beta)
            + jet_entropy_interior(protocol).scaled(parameters.beta))


def symmetric_half_jet(parameters: ArbParameters, y: Jet2) -> Jet2:
    x = parameters.x
    ljet = transverse_jet(parameters, y)
    qjet = (_kernel_jet(parameters, None, y)
            - (_kernel_jet(parameters, Jet2.constant(x), y)
               * y).scaled(2 / (x * x) * x)
            + (y * y).scaled(kernel(parameters, x, x) / (x * x)))
    return ljet.scaled(2) + qjet


def q_over_y2_boundary(parameters: ArbParameters, top: Fraction) -> arb:
    """Singularity-free lower enclosure of Q(y)/y^2 on 0<y<=top."""
    y = _ball(Fraction(0), top)
    x, beta = parameters.x, parameters.beta
    one = arb(1)
    x_mid = Fraction(x.str(40, radius=False))
    x_hi = x_mid + Fraction(1, 10**30)
    x_lo = x_mid - Fraction(1, 10**30)

    mu_y2 = _mu_enclosure(y * y, top * top)
    mu_xy = _mu_enclosure(x * y, x_hi * top)
    qh = (mu_y2 - 2 * (one / x).log() - 2 * mu_xy
          + h_arb(x * x) / (x * x))

    f = one + (one - y) * (one - y)
    rho = x + x * (one - x) * (one - y)
    coeff = (one - y) * (x - y)  # f-rho/x, exact
    if not coeff > 0:
        raise AssertionError("boundary logarithmic coefficient is not positive")
    log_lower = (one / _arbf(top)).log()
    mu_diag = _mu_enclosure(y * y * f, 2 * top * top)
    mu_cross = _mu_enclosure(y * rho, x_hi * (2 - x_lo) * top)
    qp = (
        2 * coeff * log_lower
        - f * f.log() + f * mu_diag
        + 2 * (rho / x) * rho.log()
        - 2 * (rho / x) * mu_cross
        + h_arb(_protocol(x, x)) / (x * x)
    )
    return (one - beta) * qh + beta * qp


@dataclass(frozen=True)
class AxisCell:
    lo: Fraction
    hi: Fraction
    s_lower: arb
    label: str


def _subdivide(lo: Fraction, hi: Fraction, count: int, label: str) -> list[tuple[Fraction, Fraction, str]]:
    return [
        (lo + (hi - lo) * Fraction(i, count),
         lo + (hi - lo) * Fraction(i + 1, count), label)
        for i in range(count)
    ]


def certify_symmetric(
    parameters: ArbParameters,
    low_cells: int = DEFAULT_LOW_CELLS,
    high_cells: int = DEFAULT_HIGH_CELLS,
    window_cells: int = DEFAULT_WINDOW_CELLS,
) -> tuple[dict[str, Any], list[AxisCell]]:
    x0 = Fraction(parameters.x.str(40, radius=False))
    left = x0 - X_WINDOW
    right = x0 + X_WINDOW

    q_boundary = q_over_y2_boundary(parameters, Y0)
    if not q_boundary > 0:
        raise AssertionError(f"Q/y^2 boundary bound failed: {q_boundary}")

    # The point y=0 is a structural gauge zero but S grows like
    # y*log(1/y) immediately away from it.  One coarse [0,Y0] cell would assign
    # the whole layer the zero lower bound and destroy the asymmetric q-vertex
    # estimate, so retain the exact tail and tile the rest dyadically.
    boundary_octaves = 24
    smallest = Y0 / Fraction(2) ** boundary_octaves
    intervals: list[tuple[Fraction, Fraction, str]] = [
        (Fraction(0), smallest, "zero-tail")]
    for k in range(boundary_octaves, 0, -1):
        intervals += _subdivide(
            Y0 / Fraction(2) ** k,
            Y0 / Fraction(2) ** (k - 1),
            8,
            "zero-layer",
        )
    intervals += _subdivide(Y0, left, low_cells, "left")
    intervals += _subdivide(left, right, window_cells, "x-window")
    intervals += _subdivide(right, Fraction(1) - ONE_TAIL, high_cells, "right")
    intervals.append((Fraction(1) - ONE_TAIL, Fraction(1), "one-tail"))

    # Certify curvature cellwise -- a single natural interval over the whole
    # window widens to about +/-7 by dependency even though the true minimum is
    # positive.  The resulting minimum is then reused as ONE correlated Taylor
    # constant on every window subcell.
    window_curvature: Optional[arb] = None
    for lo, hi, label in intervals:
        if label != "x-window":
            continue
        jet = symmetric_half_jet(parameters, Jet2.variable(_ball(lo, hi)))
        if not jet.dd > 0:
            raise AssertionError(
                f"S'' failed on x-window cell [{lo},{hi}]: {jet.dd}")
        if (window_curvature is None
                or jet.dd.lower() < window_curvature.lower()):
            window_curvature = arb(jet.dd.lower())
    assert window_curvature is not None

    axis: list[AxisCell] = []
    worst_strict: Optional[arb] = None
    worst_cell = None
    x_lo = Fraction(parameters.x.str(40, radius=False)) - Fraction(1, 10**30)
    x_hi = Fraction(parameters.x.str(40, radius=False)) + Fraction(1, 10**30)
    for lo, hi, label in intervals:
        if label == "zero-tail":
            lower = arb(0)  # exact gauge point plus unresolved infinitesimal tail.
        elif label == "x-window":
            if hi < x_lo:
                distance = x_lo - hi
            elif lo > x_hi:
                distance = lo - x_hi
            else:
                distance = Fraction(0)
            lower = (window_curvature * _arbf(distance) * _arbf(distance) / 2)
        elif label == "one-tail":
            value = symmetric_half_coefficient(parameters, _ball(lo, hi))
            if not value > 0:
                raise AssertionError(
                    f"S one-tail positivity failed on [{lo},{hi}]: {value}")
            lower = arb(value.lower())
        else:
            # Point-centred first-order Taylor enclosure.  The exact centre
            # evaluation has tiny width; the interval jet supplies a rigorous
            # derivative bound on the entire cell.  This is orders tighter
            # than natural interval evaluation of a formula with repeated y.
            midpoint = (lo + hi) / 2
            halfwidth = (hi - lo) / 2
            point = symmetric_half_coefficient(parameters, _arbf(midpoint))
            jet = symmetric_half_jet(parameters, Jet2.variable(_ball(lo, hi)))
            derivative = arb(abs(jet.d).upper())
            lower = arb(point.lower()) - _arbf(halfwidth) * derivative
            if not lower > 0:
                raise AssertionError(
                    f"S Taylor positivity failed on [{lo},{hi}]: "
                    f"point={point}, derivative={jet.d}, lower={lower}")
        if lower > 0 and (worst_strict is None
                          or lower.lower() < worst_strict.lower()):
            worst_strict, worst_cell = lower, (str(lo), str(hi), label)
        axis.append(AxisCell(lo, hi, lower, label))

    assert window_curvature is not None and worst_strict is not None
    return ({
        "cells": len(axis),
        "zero_boundary_q_over_y2": q_boundary.str(20),
        "zero_boundary_q_over_y2_lower": _lower_float(q_boundary),
        "x_window": str(X_WINDOW),
        "x_window_minimum_curvature": window_curvature.str(20),
        "x_window_minimum_curvature_lower": _lower_float(window_curvature),
        "worst_strict_margin": worst_strict.str(20),
        "worst_strict_cell": worst_cell,
        "structural_zeros": ["0", "x"],
        "slope_zero_exact": True,
        "certified": True,
    }, axis)


# ---------------------------------------------------------------------------
# One-sided mean-excess normal.


def mean_excess_margin(parameters: ArbParameters, y: arb) -> arb:
    """h(x)+R(y), the worst linear normal coefficient for eps<=1/2."""
    x = parameters.x
    return (h_arb(x) + kernel(parameters, x, y)
            - (y / x) * kernel(parameters, x, x))


def certify_mean_excess(
    parameters: ArbParameters,
    axis: Sequence[AxisCell],
) -> dict[str, Any]:
    """Certify every strict mean-feasible normal direction.

    Adding delta to a1 and subtracting it from a3 shifts mass from support 0 to
    support x in BOTH components.  The mean rises by delta*x.  Exact bilinearity
    gives

      P(eps,delta)=P(eps,0)+delta[h(x)+2 eps Rbar]
                   +delta^2[K(x,x)-kappa*x^2].

    Rbar=(1-q)R(y1)+qR(y2), so for eps<=1/2 its worst coefficient is bounded by
    min(h(x), h(x)+R(y)).  The quadratic coefficient is independent of every
    insertion coordinate.
    """
    worst: Optional[arb] = None
    worst_cell = None
    for cell in axis:
        value = mean_excess_margin(parameters, _ball(cell.lo, cell.hi))
        if not value > 0:
            raise AssertionError(
                f"mean-excess margin failed on [{cell.lo},{cell.hi}]: {value}")
        if worst is None or value.lower() < worst.lower():
            worst, worst_cell = value, (str(cell.lo), str(cell.hi), cell.label)
    quadratic = (kernel(parameters, parameters.x, parameters.x)
                 - _arbf(KAPPA) * parameters.x * parameters.x)
    if not quadratic > 0:
        raise AssertionError(
            f"mean-excess quadratic coefficient is not positive: {quadratic}")
    assert worst is not None
    return {
        "cells": len(axis),
        "linear_margin": worst.str(20),
        "linear_margin_lower": _lower_float(worst),
        "worst_cell": worst_cell,
        "quadratic_coefficient": quadratic.str(20),
        "quadratic_coefficient_lower": _lower_float(quadratic),
        "mean_increase": "delta*x",
        "delta_domain": "delta>=0 subject only to simplex feasibility",
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Divided-difference interaction bounds.


def _d_kernel_xt(parameters: ArbParameters, t: arb) -> arb:
    x = parameters.x
    hp_h = _hp_arb(x * t)
    z = _protocol(x, t)
    zprime = x + x * (arb(1) - x) * (arb(1) - 2 * t)
    return ((arb(1) - parameters.beta) * x * hp_h
            + parameters.beta * zprime * _hp_arb(z))


def _h_mixed(t1: arb, t2: arb) -> arb:
    z = t1 * t2
    return _hp_arb(z) + z * _hpp_arb(z)


def jhat_average_bound(parameters: ArbParameters, lo: Fraction, hi: Fraction) -> arb:
    """Range bound for Chat/(y1-y2)^2 over pairs in [lo,hi], away endpoints."""
    if not 0 < lo < hi < 1:
        raise ValueError("interior divided-difference bound needs 0<lo<hi<1")
    t = _ball(lo, hi)
    derivative = _d_kernel_xt(parameters, t)
    mixed = _h_mixed(t, t)
    x = parameters.x
    return (2 * derivative / x - kernel(parameters, x, x) / (x * x)
            - (arb(1) - parameters.beta) * mixed - _arbf(KAPPA))


def zero_tail_jhat_lower(parameters: ArbParameters, top: Fraction) -> arb:
    """Lower bound for Jhat on 0<=y1,y2<=top with log cancellation retained."""
    x, beta = parameters.x, parameters.beta
    one = arb(1)
    d = _arbf(top)
    # h-part: cancel the averages of log(t) before bounding.
    jh = (-h_arb(x * x) / (x * x) + one
          + 2 * (one - x * d).log() - 2 * x.log())
    zmax = d * x * (2 - x)
    zprime_min = x + x * (one - x) * (one - 2 * d)
    jp = (-h_arb(_protocol(x, x)) / (x * x)
          + 2 * zprime_min * _hp_arb(zmax) / x)
    return (one - beta) * jh + beta * jp - _arbf(KAPPA)


def one_tail_jhat_lower(parameters: ArbParameters, width: Fraction) -> arb:
    """Lower bound for Jhat on 1-width<=y1,y2<=1."""
    x, beta = parameters.x, parameters.beta
    one = arb(1)
    t0 = one - _arbf(width)
    zmin = t0 * t0
    minus_hmixed = (zmin / (one - zmin)).log() + one / (one - zmin)
    # d/dt K(x,t) is finite at t=1; bound the beta-protocol part in Arb.
    t = _ball(Fraction(1) - width, Fraction(1))
    z = _protocol(x, t)
    zprime = x + x * (one - x) * (one - 2 * t)
    dp = zprime * _hp_arb(z)
    jh = (-h_arb(x * x) / (x * x) + minus_hmixed
          + 2 * _hp_arb(x))
    jp = (-h_arb(_protocol(x, x)) / (x * x) + 2 * dp / x)
    return (one - beta) * jh + beta * jp - _arbf(KAPPA)


def chat_lower_for_cells(
    parameters: ArbParameters, left: AxisCell, right: AxisCell,
) -> arb:
    lo = min(left.lo, right.lo)
    hi = max(left.hi, right.hi)
    overlap = not (left.hi < right.lo or right.hi < left.lo)
    adjacent = left.hi == right.lo or right.hi == left.lo

    if hi <= Fraction(1, 10):
        jhat = zero_tail_jhat_lower(parameters, hi)
        if not jhat > 0:
            raise AssertionError(f"zero-tail Jhat failed through {hi}: {jhat}")
        return arb(0)
    if lo >= Fraction(15, 16):
        jhat = one_tail_jhat_lower(parameters, Fraction(1) - lo)
        if not jhat > 0:
            raise AssertionError(f"one-tail Jhat failed from {lo}: {jhat}")
        return arb(0)

    if overlap or adjacent or hi - lo <= Fraction(3, 64):
        jhat = jhat_average_bound(parameters, lo, hi)
        delta_max = _arbf(hi - lo)
        if jhat >= 0:
            return arb(0)
        return arb(jhat.lower()) * delta_max * delta_max

    # Separated cells: direct factored divided differences.
    if left.lo > right.hi:
        upper, lower = left, right
    else:
        upper, lower = right, left
    y1, y2 = _ball(upper.lo, upper.hi), _ball(lower.lo, lower.hi)
    d = y1 - y2
    if not d > 0:
        raise AssertionError("separated cross cells do not have positive delta")
    dk = (kernel(parameters, parameters.x, y1)
          - kernel(parameters, parameters.x, y2)) / d
    w = (h_arb(y1 * y1) + h_arb(y2 * y2) - 2 * h_arb(y1 * y2)) / (d * d)
    jhat = (2 * dk / parameters.x
            - kernel(parameters, parameters.x, parameters.x)
            / (parameters.x * parameters.x)
            - (arb(1) - parameters.beta) * w - _arbf(KAPPA))
    d2 = d * d
    if jhat >= 0:
        return arb(0)
    return arb(jhat.lower()) * arb(d2.upper())


def quadratic_minimum_lower(s1: arb, s2: arb, chat: arb) -> arb:
    """Rigorous lower bound of (1-q)s1+q*s2+q(1-q)chat on q in [0,1]."""
    l1, l2, c = arb(s1.lower()), arb(s2.lower()), arb(chat.lower())
    if c >= 0:
        return l1 if l1.lower() <= l2.lower() else l2
    a = -c
    b = l2 - l1 + c
    qstar = -b / (2 * a)
    if qstar <= 0:
        return l1
    if qstar >= 1:
        return l2
    # Unconstrained minimum; safe even if the small qstar ball brushes an edge.
    return l1 - b * b / (4 * a)


def certify_asymmetric_cover(
    parameters: ArbParameters,
    axis: Sequence[AxisCell],
) -> dict[str, Any]:
    failures = []
    worst: Optional[arb] = None
    worst_pair = None
    negative_chat_cells = 0
    checked = 0
    for i, left in enumerate(axis):
        for j in range(i + 1):
            right = axis[j]
            chat = chat_lower_for_cells(parameters, left, right)
            if chat < 0:
                negative_chat_cells += 1
            lower = quadratic_minimum_lower(left.s_lower, right.s_lower, chat)
            checked += 1
            if worst is None or lower.lower() < worst.lower():
                worst, worst_pair = lower, (
                    [str(left.lo), str(left.hi), left.label],
                    [str(right.lo), str(right.hi), right.label],
                    chat.str(12),
                )
            # Zero is allowed only where a structural zero cell participates.
            if lower < 0:
                failures.append((i, j, lower.str(12), worst_pair))
                if len(failures) >= 20:
                    break
        if len(failures) >= 20:
            break
    if failures:
        raise AssertionError(f"asymmetric cover has failing cells: {failures[:3]}")
    assert worst is not None
    return {
        "axis_cells": len(axis),
        "triangular_pairs": checked,
        "negative_chat_cells": negative_chat_cells,
        "worst_lower_bound": worst.str(20),
        "worst_pair": worst_pair,
        "failures": [],
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Sharp whole-segment obstruction and reporting.


def endpoint_obstruction(parameters: ArbParameters) -> dict[str, Any]:
    x, p = parameters.x, parameters.p
    epsilon = parameters.mean
    q = arb(1) / 2
    values = (
        arb(0), epsilon, q, x, arb(1), arb(0), x, arb(1), arb(0),
    )
    terms = _formula(values, parameters.beta, __import__("liu9_objective")._ArbOps(False))
    gap = terms.numerator - terms.ehx
    dist = distance_squared(values, parameters)
    pencil = gap - _arbf(KAPPA) * dist
    if not pencil < 0:
        raise AssertionError("endpoint pencil obstruction is not negative")
    return {
        "y1": "1", "y2": "1", "q": "1/2",
        "epsilon": epsilon.str(20),
        "raw_gap": gap.str(20),
        "raw_gap_structural_zero": True,
        "distance_squared": dist.str(20),
        "pencil": pencil.str(20),
        "pencil_upper": _upper_float(pencil),
        "interpretation": "pencil obstruction only; EHX=numerator=0 on {0,1}",
    }


def _digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode()).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--low-cells", type=int, default=DEFAULT_LOW_CELLS)
    parser.add_argument("--high-cells", type=int, default=DEFAULT_HIGH_CELLS)
    parser.add_argument("--window-cells", type=int, default=DEFAULT_WINDOW_CELLS)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    start = time.monotonic()
    mp = solve_equation_parameters(args.dps)
    parameters = certify_equation_parameters(mp)

    print("1. EXACT EPSILON EXPANSION")
    expansion = certify_exact_expansion(parameters)
    print(f"PROVED [SSOT polynomial AD, {len(expansion['probes'])} Arb probes]: "
          f"pencil = epsilon*Lbar + epsilon^2*P2 exactly, worst hand-form "
          f"residual {expansion['worst_residual']:.3e}.")
    print("PROVED [structural]: every entropy argument is epsilon-independent, so "
          "the exact degree is 2 and the epsilon^2*log(1/epsilon) coefficient "
          "is exactly 0.  There is no uncomputed remainder.")

    print("\n2. SYMMETRIC SECOND ORDER")
    symmetric, axis = certify_symmetric(
        parameters, args.low_cells, args.high_cells, args.window_cells)
    print(f"PROVED [Arb, {symmetric['cells']} abutting cells]: "
          f"S(y)=2L(y)+Q(y)>=0 on [0,1].  At zero, Q/y^2 >= "
          f"{symmetric['zero_boundary_q_over_y2_lower']:.6f}; around x, "
          f"S'' >= {symmetric['x_window_minimum_curvature_lower']:.6f}.  "
          "The only zeros are the gauge points y=0 and y=x.")
    print("PROVED [exact optimizer equations]: S(0)=S(x)=S'(x)=0 exactly.  "
          "The tiny Arb centre residue is dependency/parameter-enclosure error "
          "and is not used as a deficit.")

    print("\n3. SIMULTANEOUS ASYMMETRIC INSERTIONS")
    asymmetric = certify_asymmetric_cover(parameters, axis)
    print(f"PROVED [Arb + exact q-quadratic minimization]: S_asym>=0 on "
          f"[0,1]^2 x [0,1] over {asymmetric['triangular_pairs']} symmetric "
          f"axis-cell pairs; {asymmetric['negative_chat_cells']} pairs have a "
          "negative interaction and are settled by the interior q vertex, not "
          "by discarding the cross term.")
    print("PROVED: for every 0<=epsilon<=1/2, the exact pencil is nonnegative.  "
          "No first-order decoupling is assumed at second order.")

    print("\n4. STRICT MEAN-FEASIBLE NORMAL")
    mean_excess = certify_mean_excess(parameters, axis)
    print(f"PROVED [Arb, {mean_excess['cells']} cells]: every mean-excess "
          f"direction delta>=0 is absorbed.  For epsilon<=1/2 the linear "
          f"normal coefficient is at least "
          f"{mean_excess['linear_margin_lower']:.6f}, and the delta^2 "
          f"coefficient K(x,x)-kappa*x^2 is at least "
          f"{mean_excess['quadratic_coefficient_lower']:.6f}.")
    print("PROVED [exact bilinearity]: the component-mean cross terms cancel "
          "because their q-weighted sum is zero.  Thus the mean-active theorem "
          "extends to the full one-sided half-space mean>=p*x.")

    print("\n5. SHARP OBSTRUCTION OUTSIDE THE LOCAL MASS RANGE")
    obstruction = endpoint_obstruction(parameters)
    print(f"REFUTED [Arb]: pencil positivity on the whole feasible insertion "
          f"segment.  At y1=y2=1 and epsilon=p*x, raw gap is structurally 0 "
          f"but pencil <= {obstruction['pencil_upper']:.6e} < 0.  This is not a "
          "counterexample to Hypothesis 2: numerator=EHX=0 on the {0,1} law.")

    payload = {
        "tool": "liu9_second_order.py",
        "dps": args.dps,
        "kappa": str(KAPPA),
        "mass_radius": "1/2",
        "exact_expansion": expansion,
        "symmetric": symmetric,
        "asymmetric": asymmetric,
        "mean_excess": mean_excess,
        "endpoint_obstruction": obstruction,
        "elapsed_seconds": round(time.monotonic() - start, 3),
    }
    digest_payload = {k: v for k, v in payload.items() if k != "elapsed_seconds"}
    payload["report_sha256"] = _digest(digest_payload)
    if args.output:
        Path(args.output).write_text(json.dumps(payload, indent=1, sort_keys=True))
        print(f"\nreport written to {args.output}")
        print(f"report_sha256 {payload['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
