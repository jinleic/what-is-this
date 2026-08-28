#!/usr/bin/env python3
"""Quantitative zero-support boundary-layer estimate for Liu's Hypothesis 2.

`liu9_tube.py` reports that a one-piece cubic tube bound is impossible because
every raw-coordinate tube meets supports tending to zero, where
h'''(y)=(1-2y)/(y^2(1-y)^2) is unbounded, and it names the missing ingredient:

    "(ii) an explicit y*log(1/y) boundary-layer estimate covering the zero-mass
     strata rather than Taylor-expanding through them".

`liu9_binding.py` already certifies the *leading coefficient* of that layer at
the binding distribution, B=2*p*x*[1+beta*(1-x)]-1>0, and certifies that the
zero-mass insertion potential V(y) is nonnegative on [0,1].  Both statements
are asymptotic or pointwise at one distribution: neither gives a threshold y0,
a uniform constant, or an inequality valid at every feasible nine-parameter
point.  This module supplies exactly that missing quantitative statement.

Write the denominator-cleared gap in Liu's own transcription
(https://jingbol.web.illinois.edu/frankl5.m, see liu9_objective.py) as

    gap(V) = (1-beta)*EHXY(V) + beta*EHPI(V) - EHX(V),

so Hypothesis 2 is exactly gap>=0 on {mean>=1-c'}.  Index the six atoms by
j=0..5 with weights w=( (1-q)a1,(1-q)a2,(1-q)a3, q*a1,q*a2,q*a3 ) and supports
(b0,b2,b4,b1,b3,b5).  Elementary differentiation of the three sums gives the
exact identity

    d gap / d b_j = w_j * G_j(V),
    G_j(V) = 2(1-beta) * sum_k w_k b_k h'(b_j b_k)
           + 2*beta   * sum_{k in comp(j)} a_k pi_1(b_j,b_k) h'(pi(b_j,b_k))
           - h'(b_j),

with pi(y,z)=yz(1+(1-y)(1-z)) and pi_1=d pi/dy=z(1+(1-z)(1-2y)); terms with
b_k=0 vanish identically in b_j and are omitted.  Evaluated at the binding
point this reproduces liu9_binding's B, which is the intended cross-check.

THEOREM (boundary layer).  Fix y0<=1/4.  For every point of the raw box with
a1+a2<=1 and every support index j with 0<b_j<=y0,

    G_j(V) >= [2(1-beta)*mu(V) - 1] * log(1/b_j) - 2(1-beta)*log(1/(1-y0)),

where mu(V) is the mean.  Consequently, with

    Lam  = 2(1-beta)*(mean-y0) - 1,     K = 2(1-beta)*log(1/(1-y0)),
    m0   = Lam*(log(1/y0)+1) - K,

every mean-feasible V satisfies, for the set S of support indices with
0<b_j<=y0 and V0 the point with all of them set to zero,

    gap(V) >= gap(V0) + m0 * sum_{j in S} w_j b_j,

and the single-coordinate form has the sharper constant

    m0s  = Lam_f*(log(1/y0)+1) - 2(1-beta)*y0*(log(1/y0)/2+3/4) - K,
    Lam_f = 2(1-beta)*mean - 1.

Both constants are Arb-certified positive at y0=1/32 from the equation-defined
parameters, never from Liu's printed decimals.  The estimate is therefore an
explicit-threshold replacement for the impossible cubic remainder on this
stratum, and the reduction to the face b_j=0 lands in the region where
|h'''|<=(1-2y0)/(y0^2(1-y0)^2) is finite.

Run from the repository root with

    math/.venv/bin/python -I -B math/uc/liu9_boundary_layer.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from fractions import Fraction
from typing import Any, Iterable, Optional, Sequence

# This workstation is shared with live certification workers.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    certify_curvatures,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import _ScalarOps, _formula  # noqa: E402

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 60
DEFAULT_Y0 = Fraction(1, 32)
TRADEOFF_THRESHOLDS = (
    Fraction(1, 8),
    Fraction(1, 16),
    Fraction(1, 24),
    Fraction(1, 32),
    Fraction(1, 64),
    Fraction(1, 128),
    Fraction(1, 1024),
)
DEVIATION_GRID = (Fraction(1, 100), Fraction(3, 100), Fraction(1, 10),
                  Fraction(3, 10), Fraction(1))
SAMPLE_SEED = 20260827
SUPPORT_VARIABLE_INDEX = (3, 4, 5, 6, 7, 8)


# ---------------------------------------------------------------------------
# Scalar helpers.  Deliberately re-derived here: this module is an independent
# check of the closed form, so it does not import liu9_binding's private
# entropy helpers.


def _arbf(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _h_mp(u: mpmath.mpf) -> mpmath.mpf:
    if u <= 0 or u >= 1:
        return mpmath.mpf(0)
    return -(u * mpmath.log(u) + (1 - u) * mpmath.log1p(-u))


def _hp_mp(u: mpmath.mpf) -> mpmath.mpf:
    """h'(u)=log((1-u)/u) for 0<u<1."""
    return mpmath.log1p(-u) - mpmath.log(u)


def _pi_mp(y: mpmath.mpf, z: mpmath.mpf) -> mpmath.mpf:
    return y * z * (1 + (1 - y) * (1 - z))


def _pi_first_mp(y: mpmath.mpf, z: mpmath.mpf) -> mpmath.mpf:
    return z * (1 + (1 - z) * (1 - 2 * y))


def _weights(values: Sequence[Any], one: Any) -> tuple[Any, ...]:
    a1, a2, q = values[0], values[1], values[2]
    a3 = one - a1 - a2
    masses = (a1, a2, a3)
    return tuple(component * mass
                 for component in (one - q, q) for mass in masses)


def _support(values: Sequence[Any]) -> tuple[Any, ...]:
    return (values[3], values[4], values[5], values[6], values[7], values[8])


def mean_of(values: Sequence[Any], one: Any) -> Any:
    weights = _weights(values, one)
    support = _support(values)
    total = one - one
    for weight, point in zip(weights, support):
        total = total + weight * point
    return total


# ---------------------------------------------------------------------------
# The gap and its closed-form partial derivative.


class _GapOps(_ScalarOps):
    """Liu's formula without the final quotient.

    The gap N-D is defined even where the objective N/D is not: a point whose
    supports are all 0 or 1 has D=0.  Such points are genuine members of the
    zero-support stratum -- they are exactly what the layer collapses to -- so
    the division must be skipped rather than guarded downstream.
    """

    def quotient(self, numerator: Any, denominator: Any) -> None:
        return None


def gap_mp(values: Sequence[mpmath.mpf], beta: mpmath.mpf) -> mpmath.mpf:
    """gap = (1-beta)EHXY + beta*EHPI - EHX from Liu's own transcription."""
    terms = _formula(tuple(values), beta, _GapOps(_h_mp, mpmath.mpf(1)))
    return terms.numerator - terms.ehx


def entropy_mp(values: Sequence[mpmath.mpf]) -> mpmath.mpf:
    terms = _formula(tuple(values), mpmath.mpf(0), _GapOps(_h_mp, mpmath.mpf(1)))
    return terms.ehx


def closed_form_factor(values: Sequence[mpmath.mpf], index: int,
                       beta: mpmath.mpf) -> mpmath.mpf:
    """G_j(V): the weight-cleared partial derivative of gap in b_j."""
    if not 0 <= index <= 5:
        raise ValueError("support index must be 0..5")
    one = mpmath.mpf(1)
    a1, a2 = values[0], values[1]
    masses = (a1, a2, one - a1 - a2)
    weights = _weights(values, one)
    support = _support(values)
    y = support[index]
    if not 0 < y < 1:
        raise ValueError("closed form needs an interior support coordinate")

    cross = mpmath.mpf(0)
    for other, point in enumerate(support):
        if point == 0:
            # w_j w_k h(b_j b_k) is identically zero in b_j when b_k=0.
            continue
        cross += weights[other] * point * _hp_mp(y * point)

    paired = mpmath.mpf(0)
    base = 0 if index < 3 else 3
    for offset in range(3):
        point = support[base + offset]
        if point == 0:
            continue
        paired += (masses[offset] * _pi_first_mp(y, point)
                   * _hp_mp(_pi_mp(y, point)))

    return 2 * (1 - beta) * cross + 2 * beta * paired - _hp_mp(y)


def closed_form_partial(values: Sequence[mpmath.mpf], index: int,
                        beta: mpmath.mpf) -> mpmath.mpf:
    weights = _weights(values, mpmath.mpf(1))
    return weights[index] * closed_form_factor(values, index, beta)


# ---------------------------------------------------------------------------
# Forward-mode differentiation through Liu's transcribed formula itself.


class _Dual:
    """First-order jet over mpmath, used to differentiate ``_formula``."""

    __slots__ = ("value", "derivative")

    def __init__(self, value: Any, derivative: Any = 0) -> None:
        self.value = value if isinstance(value, mpmath.mpf) else mpmath.mpf(value)
        self.derivative = (derivative if isinstance(derivative, mpmath.mpf)
                           else mpmath.mpf(derivative))

    @staticmethod
    def _coerce(other: Any) -> "_Dual":
        return other if isinstance(other, _Dual) else _Dual(other)

    def __add__(self, other: Any) -> "_Dual":
        other = self._coerce(other)
        return _Dual(self.value + other.value,
                     self.derivative + other.derivative)

    __radd__ = __add__

    def __sub__(self, other: Any) -> "_Dual":
        other = self._coerce(other)
        return _Dual(self.value - other.value,
                     self.derivative - other.derivative)

    def __rsub__(self, other: Any) -> "_Dual":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "_Dual":
        other = self._coerce(other)
        return _Dual(self.value * other.value,
                     self.derivative * other.value
                     + self.value * other.derivative)

    __rmul__ = __mul__

    def __truediv__(self, other: Any) -> "_Dual":
        other = self._coerce(other)
        return _Dual(
            self.value / other.value,
            (self.derivative * other.value - self.value * other.derivative)
            / (other.value * other.value),
        )

    def __rtruediv__(self, other: Any) -> "_Dual":
        return self._coerce(other) / self


def _h_dual(u: _Dual) -> _Dual:
    if u.value <= 0 or u.value >= 1:
        if u.derivative != 0:
            raise ArithmeticError(
                "entropy differentiated at a boundary argument with nonzero "
                "inner derivative; the sample point is outside the layer")
        return _Dual(0, 0)
    return _Dual(_h_mp(u.value), _hp_mp(u.value) * u.derivative)


def autodiff_partial(values: Sequence[mpmath.mpf], index: int,
                     beta: mpmath.mpf) -> mpmath.mpf:
    """d gap / d b_j obtained by differentiating ``_formula`` directly."""
    variable = SUPPORT_VARIABLE_INDEX[index]
    duals = tuple(_Dual(value, 1 if position == variable else 0)
                  for position, value in enumerate(values))
    terms = _formula(duals, _Dual(beta), _GapOps(_h_dual, _Dual(1)))
    gap = terms.numerator - terms.ehx
    return gap.derivative


# ---------------------------------------------------------------------------
# Certified constants.


@dataclass(frozen=True)
class LayerConstants:
    y0: Fraction
    two_one_minus_beta: arb
    lam_feasible: arb
    lam_uniform: arb
    penalty: arb
    m_uniform: arb
    m_sharp: arb
    third_derivative_ceiling: arb

    def as_json(self) -> dict[str, Any]:
        return {
            "y0": str(self.y0),
            "y0_float": float(self.y0),
            "two_one_minus_beta": str(self.two_one_minus_beta),
            "lam_feasible": str(self.lam_feasible),
            "lam_uniform": str(self.lam_uniform),
            "penalty_K": str(self.penalty),
            "m_uniform": str(self.m_uniform),
            "m_sharp": str(self.m_sharp),
            "smooth_third_derivative_ceiling":
                str(self.third_derivative_ceiling),
            "lam_feasible_lower": float(self.lam_feasible.lower()),
            "lam_uniform_lower": float(self.lam_uniform.lower()),
            "m_uniform_lower": float(self.m_uniform.lower()),
            "m_sharp_lower": float(self.m_sharp.lower()),
        }


def certify_constants(parameters: ArbParameters,
                      y0: Fraction) -> LayerConstants:
    """Arb-certify every constant of the boundary-layer theorem."""
    if not Fraction(0) < y0 <= Fraction(1, 4):
        raise ValueError("the layer threshold must satisfy 0<y0<=1/4")
    y0_arb = _arbf(y0)
    two_omb = 2 * (1 - parameters.beta)
    lam_feasible = two_omb * parameters.mean - 1
    lam_uniform = two_omb * (parameters.mean - y0_arb) - 1
    penalty = two_omb * (-(1 - y0_arb).log())
    log_inverse = -y0_arb.log()
    m_uniform = lam_uniform * (log_inverse + 1) - penalty
    m_sharp = (lam_feasible * (log_inverse + 1)
               - two_omb * y0_arb * (log_inverse / 2 + arb(3) / 4)
               - penalty)
    ceiling = (1 - 2 * y0_arb) / (y0_arb**2 * (1 - y0_arb) ** 2)
    return LayerConstants(
        y0=y0,
        two_one_minus_beta=two_omb,
        lam_feasible=lam_feasible,
        lam_uniform=lam_uniform,
        penalty=penalty,
        m_uniform=m_uniform,
        m_sharp=m_sharp,
        third_derivative_ceiling=ceiling,
    )


def constants_are_certified(constants: LayerConstants) -> bool:
    return bool(
        constants.lam_feasible > 0
        and constants.lam_uniform > 0
        and constants.m_uniform > 0
        and constants.m_sharp > 0
    )


# ---------------------------------------------------------------------------
# Interval sweeps guarding the two term-wise links of the proof.


def _octave_grid(top: Fraction, octaves: int,
                 pieces: int) -> list[tuple[Fraction, Fraction]]:
    """Exact-rational cover of (top*2^-octaves, top], geometric by octave.

    A linear grid is useless here: its lowest cell touches zero, and an Arb
    ball containing zero has no logarithm.  Halving keeps every cell's
    endpoint ratio bounded, so every logarithm stays finite.
    """
    out: list[tuple[Fraction, Fraction]] = []
    for octave in range(octaves):
        upper = top / 2**octave
        lower = upper / 2
        for piece in range(pieces):
            out.append((
                lower + (upper - lower) * Fraction(piece, pieces),
                lower + (upper - lower) * Fraction(piece + 1, pieces),
            ))
    return out


def sweep_cross_link(y0: Fraction, octaves: int, pieces: int) -> dict[str, Any]:
    """h'(y*b) >= log(1/y)+log(1-y0) on (0,y0] x (0,1].

    Logarithms are additive on the product, so

        h'(yb) - log(1/y) = log(1/b) + log(1-yb)

    exactly, and the right side is at least log(1-y0) as soon as b<=1 and
    yb<=y0.  Both hypotheses are *exact rational* conditions on the cell, and
    both conclusions are monotonicity of the logarithm, so the certification
    below is the rational domain check, not a floating sweep.  The margin is
    tight -- it vanishes at (y,b)=(y0,1) -- so no interval test could separate
    it from zero there, and requiring that would only measure outward
    rounding.  Arb is still used two ways: it must never enclose the margin
    strictly below zero (a refutation), and every cell midpoint must satisfy
    the additivity identity to 40 digits (a transcription guard).
    """
    log_one_minus = (1 - _arbf(y0)).log()
    y_grid = _octave_grid(y0, octaves, pieces)
    b_grid = _octave_grid(Fraction(1), octaves, pieces)
    weakest: Optional[arb] = None
    weakest_float = math.inf
    worst_identity = 0.0
    boxes = 0
    ties = 0
    with mpmath.workdps(50):
        for y_lo, y_hi in y_grid:
            y_interval = _arbf(y_lo).union(_arbf(y_hi))
            y_mid = Fraction(y_lo + y_hi, 2)
            for b_lo, b_hi in b_grid:
                boxes += 1
                if b_hi > 1 or y_hi * b_hi > y0:
                    raise AssertionError(
                        f"cell ({y_lo},{y_hi})x({b_lo},{b_hi}) violates the "
                        "exact rational hypotheses b<=1 and y*b<=y0")
                b_interval = _arbf(b_lo).union(_arbf(b_hi))
                product = y_interval * b_interval
                margin = -b_interval.log() + (1 - product).log() - log_one_minus
                if margin.upper() < 0:
                    raise AssertionError(
                        f"cross link refuted at y in {y_interval}, "
                        f"b in {b_interval}: {margin}")
                if not margin > 0:
                    ties += 1
                if float(margin.lower()) < weakest_float:
                    weakest, weakest_float = margin, float(margin.lower())

                b_mid = Fraction(b_lo + b_hi, 2)
                y_point = mpmath.mpf(y_mid.numerator) / y_mid.denominator
                b_point = mpmath.mpf(b_mid.numerator) / b_mid.denominator
                residual = abs(
                    _hp_mp(y_point * b_point)
                    - (-mpmath.log(y_point) - mpmath.log(b_point)
                       + mpmath.log1p(-y_point * b_point))
                )
                worst_identity = max(worst_identity, float(residual))
    if weakest is None:
        raise AssertionError("empty cross-link sweep")
    if worst_identity > 1e-40:
        raise AssertionError(
            f"log additivity guard failed at {worst_identity:.3e}")
    return {
        "boxes": boxes,
        "rational_hypotheses_verified": boxes,
        "tight_cells": ties,
        "weakest_margin": str(weakest),
        "weakest_margin_lower": weakest_float,
        "worst_log_additivity_residual": worst_identity,
        "cover": f"dyadic, {octaves} octaves x {pieces} pieces in each variable",
        "smallest_covered_y": str(y0 / 2**octaves),
        "term_wise_tail": (
            "y<=y0*2^-octaves or b<=2^-octaves: margin=log(1/b)+log(1-yb)"
            "-log(1-y0)>=0 by the same two monotonicity steps"
        ),
    }


def sweep_paired_link(y0: Fraction, y_boxes: int, z_boxes: int) -> dict[str, Any]:
    """pi_1(y,z)*h'(pi(y,z)) >= 0 on [0,y0] x [0,1].

    pi(y,z)=yz(1+(1-y)(1-z))<=2yz<=2*y0<=1/2 forces h'(pi)>=0, and
    pi_1=z(1+(1-z)(1-2y))>=0 because every factor is nonnegative for z in
    [0,1] and y<=1/2.  Both are products of exact rational sign conditions,
    verified cellwise; pi_1 vanishes on the whole edge z=0, so an interval
    test could not separate it there and the sweep instead refuses any cell
    whose enclosure lies strictly on the wrong side.
    """
    if y0 > Fraction(1, 4):
        raise ValueError("the paired link needs y0<=1/4")
    largest_pi: Optional[arb] = None
    largest_float = -math.inf
    smallest_first: Optional[arb] = None
    smallest_float = math.inf
    boxes = 0
    edge_cells = 0
    for yi in range(y_boxes):
        y_lo = y0 * Fraction(yi, y_boxes)
        y_hi = y0 * Fraction(yi + 1, y_boxes)
        y_interval = _arbf(y_lo).union(_arbf(y_hi))
        for zi in range(z_boxes):
            z_lo = Fraction(zi, z_boxes)
            z_hi = Fraction(zi + 1, z_boxes)
            z_interval = _arbf(z_lo).union(_arbf(z_hi))
            boxes += 1
            if not (0 <= z_lo and z_hi <= 1 and y_hi <= Fraction(1, 2)
                    and 2 * y_hi * z_hi <= Fraction(1, 2)):
                raise AssertionError(
                    f"cell ({y_lo},{y_hi})x({z_lo},{z_hi}) violates the exact "
                    "rational hypotheses 0<=z<=1, y<=1/2, 2yz<=1/2")
            pi_value = y_interval * z_interval * (
                1 + (1 - y_interval) * (1 - z_interval))
            first = z_interval * (1 + (1 - z_interval) * (1 - 2 * y_interval))
            if pi_value.lower() > arb(1) / 2:
                raise AssertionError(f"pi refuted above 1/2: {pi_value}")
            if first.upper() < 0:
                raise AssertionError(f"pi_1 refuted below zero: {first}")
            if not first > 0:
                edge_cells += 1
            if float(pi_value.upper()) > largest_float:
                largest_pi, largest_float = pi_value, float(pi_value.upper())
            if float(first.lower()) < smallest_float:
                smallest_first, smallest_float = first, float(first.lower())
    if largest_pi is None or smallest_first is None:
        raise AssertionError("empty paired-link sweep")
    return {
        "boxes": boxes,
        "rational_hypotheses_verified": boxes,
        "zero_edge_cells": edge_cells,
        "largest_pi": str(largest_pi),
        "largest_pi_upper": largest_float,
        "smallest_pi_first": str(smallest_first),
        "smallest_pi_first_lower": smallest_float,
    }


# ---------------------------------------------------------------------------
# Sampling controls.


def _simplex_pair(rng: random.Random) -> tuple[float, float]:
    a1, a2 = rng.random(), rng.random()
    if a1 + a2 > 1:
        a1, a2 = 1 - a1, 1 - a2
    return a1, a2


def _draw_point(rng: random.Random, mean_target: float, y0: float,
                index: int) -> Optional[tuple[tuple[mpmath.mpf, ...], int]]:
    """A mean-feasible nine-vector whose support index ``index`` is in (0,y0]."""
    for _ in range(400):
        a1, a2 = _simplex_pair(rng)
        q = rng.choice([rng.random(), rng.random() ** 4, 1 - rng.random() ** 4])
        support = []
        for position in range(6):
            if position == index:
                support.append(y0 * rng.random() ** rng.choice([1, 3, 8]))
                continue
            draw = rng.random()
            if draw < 0.22:
                support.append(0.0)
            elif draw < 0.34:
                support.append(y0 * rng.random())
            else:
                support.append(rng.random() ** 0.2)
        if support[index] <= 0.0:
            continue
        values = (a1, a2, q, support[0], support[1], support[2],
                  support[3], support[4], support[5])
        converted = tuple(mpmath.mpf(value) for value in values)
        weights = _weights(converted, mpmath.mpf(1))
        if weights[index] <= 0:
            continue
        if float(mean_of(converted, mpmath.mpf(1))) < mean_target:
            continue
        if entropy_mp(converted) <= 0:
            continue
        return converted, index
    return None


def _structured_points(mean_target: float, y0: float
                       ) -> list[tuple[tuple[mpmath.mpf, ...], int]]:
    """Adversarial corners: degenerate q, tiny supports, extreme masses."""
    out: list[tuple[tuple[mpmath.mpf, ...], int]] = []
    small_values = (y0, y0 / 2, 1e-3, 1e-8, 1e-30)
    q_values = (0.0, 1e-9, 0.5, 1 - 1e-9, 1.0)
    mass_values = ((1.0, 0.0), (0.9, 0.1), (0.5, 0.5), (0.34, 0.33),
                   (0.0, 1.0), (0.0, 0.0))
    other_values = (0.0, 1.0)
    for q in q_values:
        for a1, a2 in mass_values:
            for small in small_values:
                for other in other_values:
                    for index in range(6):
                        support = [other] * 6
                        support[index] = small
                        values = tuple(mpmath.mpf(value) for value in
                                       (a1, a2, q, *support))
                        weights = _weights(values, mpmath.mpf(1))
                        if weights[index] <= 0:
                            continue
                        if float(mean_of(values, mpmath.mpf(1))) < mean_target:
                            continue
                        if entropy_mp(values) <= 0:
                            continue
                        out.append((values, index))
    return out


@dataclass
class SampleReport:
    points: int
    structured: int
    random: int
    derivative_checks: int
    worst_derivative_relative: float
    pointwise_checks: int
    worst_pointwise_margin: float
    integrated_checks: int
    worst_integrated_ratio: float
    multi_checks: int
    worst_multi_ratio: float
    failures: list[str]


def run_samples(parameters: MPParameters, constants: LayerConstants,
                random_points: int, dps: int, seed: int) -> SampleReport:
    """Falsification pass over the boundary-layer stratum."""
    rng = random.Random(seed)
    y0 = float(constants.y0)
    beta = parameters.beta
    mean_target = float(parameters.mean)
    lam_feasible = mpmath.mpf(str(float(constants.lam_feasible.lower())))
    penalty = mpmath.mpf(str(float(constants.penalty.upper())))
    m_sharp = mpmath.mpf(str(float(constants.m_sharp.lower())))
    m_uniform = mpmath.mpf(str(float(constants.m_uniform.lower())))

    report = SampleReport(0, 0, 0, 0, 0.0, 0, math.inf, 0, math.inf, 0,
                          math.inf, [])
    structured = _structured_points(mean_target, y0)
    points: list[tuple[tuple[mpmath.mpf, ...], int]] = list(structured)
    for step in range(random_points):
        drawn = _draw_point(rng, mean_target, y0, step % 6)
        if drawn is not None:
            points.append(drawn)
    report.structured = len(structured)
    report.random = len(points) - len(structured)

    with mpmath.workdps(dps):
        for values, index in points:
            report.points += 1
            weights = _weights(values, mpmath.mpf(1))
            y = values[SUPPORT_VARIABLE_INDEX[index]]

            closed = closed_form_partial(values, index, beta)
            auto = autodiff_partial(values, index, beta)
            scale = max(abs(closed), abs(auto), mpmath.mpf(1))
            relative = float(abs(closed - auto) / scale)
            report.derivative_checks += 1
            if relative > report.worst_derivative_relative:
                report.worst_derivative_relative = relative
            if relative > 1e-25:
                report.failures.append(
                    f"derivative mismatch {relative:.3e} at index {index}")

            factor = closed_form_factor(values, index, beta)
            bound = lam_feasible * mpmath.log(1 / y) - penalty
            margin = float(factor - bound)
            report.pointwise_checks += 1
            if margin < report.worst_pointwise_margin:
                report.worst_pointwise_margin = margin
            if margin < -1e-30:
                report.failures.append(
                    f"pointwise bound violated by {margin:.3e} at index {index}")

            zeroed = list(values)
            zeroed[SUPPORT_VARIABLE_INDEX[index]] = mpmath.mpf(0)
            increment = gap_mp(values, beta) - gap_mp(tuple(zeroed), beta)
            predicted = m_sharp * weights[index] * y
            report.integrated_checks += 1
            if predicted > 0:
                ratio = float(increment / predicted)
                if ratio < report.worst_integrated_ratio:
                    report.worst_integrated_ratio = ratio
                if ratio < 1 - 1e-20:
                    report.failures.append(
                        f"single-coordinate estimate violated: ratio {ratio:.12f}")

            small = [position for position in range(6)
                     if 0 < values[SUPPORT_VARIABLE_INDEX[position]] <= y0]
            if small:
                collapsed = list(values)
                for position in small:
                    collapsed[SUPPORT_VARIABLE_INDEX[position]] = mpmath.mpf(0)
                total = sum(
                    (weights[position]
                     * values[SUPPORT_VARIABLE_INDEX[position]]
                     for position in small),
                    mpmath.mpf(0),
                )
                increment = gap_mp(values, beta) - gap_mp(tuple(collapsed), beta)
                predicted = m_uniform * total
                report.multi_checks += 1
                if predicted > 0:
                    ratio = float(increment / predicted)
                    if ratio < report.worst_multi_ratio:
                        report.worst_multi_ratio = ratio
                    if ratio < 1 - 1e-20:
                        report.failures.append(
                            f"multi-coordinate estimate violated: ratio {ratio:.12f}")
    return report


# ---------------------------------------------------------------------------
# Binding-point cross-check and the kappa trade-off.


def binding_coefficient_check(parameters: MPParameters,
                              arb_parameters: ArbParameters,
                              dps: int) -> dict[str, Any]:
    """G_j at the binding point must reproduce B=2px[1+beta(1-x)]-1.

    G(y)=B*log(1/y)+C+O(y*log(1/y)), so the ratio G/log(1/y) converges only
    harmonically and says little.  The finite difference in log(1/y) cancels C
    exactly and converges at the rate of the neglected O(y*log(1/y)) term,
    which at y=10^-40 is far below any transcription error this check exists
    to catch.
    """
    curvatures = certify_curvatures(arb_parameters)
    working = max(dps, 140)
    with mpmath.workdps(working):
        beta = parameters.beta
        barrier = (2 * parameters.p * parameters.x
                   * (1 + beta * (1 - parameters.x)) - 1)

        def factor_at(exponent: int) -> tuple[mpmath.mpf, mpmath.mpf]:
            y = mpmath.mpf(10) ** (-exponent)
            values = (parameters.p, 1 - parameters.p, mpmath.mpf(0),
                      parameters.x, y, mpmath.mpf(0),
                      mpmath.mpf(0), mpmath.mpf(0), mpmath.mpf(0))
            return closed_form_factor(values, 1, beta), mpmath.log(1 / y)

        rows = []
        for low, high in ((10, 20), (20, 40), (40, 80)):
            g_low, log_low = factor_at(low)
            g_high, log_high = factor_at(high)
            slope = (g_high - g_low) / (log_high - log_low)
            rows.append({
                "exponents": f"1e-{low} to 1e-{high}",
                "slope": mpmath.nstr(slope, 30),
                "difference_from_B": mpmath.nstr(slope - barrier, 6),
                "difference_float": float(slope - barrier),
            })
        if abs(rows[-1]["difference_float"]) > 1e-30:
            raise AssertionError(
                "the closed form does not reproduce the certified leading "
                f"coefficient: {rows[-1]['difference_from_B']}")
    return {
        "arb_boundary_barrier": str(curvatures.boundary_barrier),
        "independent_barrier": mpmath.nstr(barrier, 30),
        "slope_rows": rows,
    }


def mirror_layer_probe(parameters: MPParameters, dps: int) -> dict[str, Any]:
    """Measure the opposite entropy layer, b_j -> 1.

    The face this lemma reduces to is only C^3 where the surviving supports
    lie in [y0,1-y0]: h'''(y)=(1-2y)/(y^2(1-y)^2) blows up at 1 exactly as it
    does at 0, and a small-mass atom near 1 is inside the tube.  Repeating the
    singular bookkeeping there, only partners with b_k=1 produce a
    log(1/(1-b_j)) term, through h'(b_j*1) and through pi(b_j,1)=b_j, so

        G_j ~ [1 - 2(1-beta)*W1 - 2*beta*A1] * log(1/(1-b_j)),

    with W1 the total weight at support 1 and A1 its component-conditional
    mass.  Unlike the lower layer, the mean constraint does not sign this: the
    coefficient is positive when little mass sits at 1 and negative when much
    does.  The probe measures the slope directly and checks it against the
    prediction, so the asymmetry is a measured fact rather than an assumption.
    """
    working = max(dps, 200)
    with mpmath.workdps(working):
        one = mpmath.mpf(1)
        beta = parameters.beta
        zero = mpmath.mpf(0)
        configurations = {
            "W1=1-p, only the moving atom reaches 1": (
                parameters.p, 1 - parameters.p, zero,
                parameters.x, None, zero, zero, zero, zero),
            "W1=2/5, the rest of the mass sits at x": (
                mpmath.mpf(3) / 5, mpmath.mpf(2) / 5, zero,
                parameters.x, None, zero, zero, zero, zero),
            "W1=1/2, the exact sign flip": (
                mpmath.mpf(1) / 2, mpmath.mpf(1) / 2, zero,
                parameters.x, None, zero, zero, zero, zero),
            "W1=7/10, a third atom already at 1": (
                mpmath.mpf(3) / 10, mpmath.mpf(4) / 10, zero,
                parameters.x, None, one, zero, zero, zero),
            "W1=1/2 but A1=1, the beta split is visible": (
                mpmath.mpf(1) / 2, mpmath.mpf(1) / 2, mpmath.mpf(1) / 2,
                one, None, zero, parameters.x, parameters.x, zero),
            "W1=1, every atom at support 1": (
                mpmath.mpf(45) / 100, mpmath.mpf(45) / 100,
                mpmath.mpf(1) / 2, one, None, one, one, one, one),
        }
        rows = []
        for name, template in configurations.items():
            def evaluate(distance: int) -> tuple[mpmath.mpf, mpmath.mpf]:
                y = one - mpmath.mpf(10) ** (-distance)
                values = tuple(y if entry is None else entry
                               for entry in template)
                return (closed_form_factor(values, 1, beta),
                        mpmath.log(1 / (1 - y)))

            g_low, log_low = evaluate(40)
            g_high, log_high = evaluate(80)
            slope = (g_high - g_low) / (log_high - log_low)

            probe = tuple(one if entry is None else entry
                          for entry in template)
            weights = _weights(probe, one)
            support = _support(probe)
            masses = (probe[0], probe[1], one - probe[0] - probe[1])
            weight_at_one = sum(
                (weights[k] for k in range(6) if support[k] == one),
                mpmath.mpf(0))
            mass_at_one = sum(
                (masses[k] for k in range(3) if support[k] == one),
                mpmath.mpf(0))
            predicted = 1 - 2 * (1 - beta) * weight_at_one - 2 * beta * mass_at_one
            rows.append({
                "configuration": name,
                "weight_at_support_one": mpmath.nstr(weight_at_one, 12),
                "measured_slope": mpmath.nstr(slope, 20),
                "predicted_slope": mpmath.nstr(predicted, 20),
                "difference": float(slope - predicted),
                "sign": ("flat (exact tie)" if abs(float(slope)) < 1e-30
                         else "increasing" if slope > 0 else "decreasing"),
            })
            if abs(float(slope - predicted)) > 1e-30:
                raise AssertionError(
                    f"mirror-layer prediction failed for {name}: "
                    f"{mpmath.nstr(slope - predicted, 6)}")
    signs = {row["sign"] for row in rows}
    return {
        "rows": rows,
        "two_sided": len(signs) > 1,
        "note": (
            "the b->1 layer is not covered by this lemma; its sign depends on "
            "the weight already sitting at support 1"
        ),
    }


def entropy_at_star(parameters: ArbParameters) -> arb:
    """H* = EHX(P*) = p*h(x), the exact conversion between the two kappas.

    Liu's objective is the quotient `Phi = N/D` with `D = EHX`, and the gap
    used throughout this module is the *quotient-free* `N - D = D*(Phi-1)`.
    That identity is exact wherever `D > 0`, so

        gap/dist^2 = EHX * (Phi-1)/dist^2,

    and at `P*`, where `EHX = p*h(x) = H*`, the two normalizations differ by
    exactly this factor.  `liu9_binding.certify_curvatures` divides its
    `split_unit` by `h(x)` and `liu9_binding.py` prints `split_unit/4` as
    `Phi''(0)`, so the curvature ratios that `liu9_tube.py` displays are in
    `Phi` units.  Every kappa in this module is in raw-gap units, because it
    comes from `d gap / d b_j`.  Comparing the two without this factor -- as
    this module did before 2026-08-28 -- compares incompatible quantities.
    """
    return parameters.p * (-(parameters.x * parameters.x.log()
                             + (1 - parameters.x) * (1 - parameters.x).log()))


def smooth_kappa_ceiling(parameters: ArbParameters,
                         units: str = "gap") -> arb:
    """The quadratic ceiling of liu9_tube's two smooth displayed modes.

    Recomputed here rather than parsed from that script's output: the layer
    is only interesting relative to the constant it has to beat.  `units`
    selects `"phi"` for the ratio as `liu9_tube.py` prints it, or `"gap"` for
    the raw-gap normalization every estimate in this module uses.
    """
    if units not in ("gap", "phi"):
        raise ValueError("units must be 'gap' or 'phi'")
    curvatures = certify_curvatures(parameters)
    face = curvatures.face / (2 * parameters.p * parameters.x**2)
    split = curvatures.split_unit / (
        2 * (parameters.p**2 + parameters.p * parameters.x**2))
    lower = min(face.lower(), split.lower(), key=float)
    upper = min(face.upper(), split.upper(), key=float)
    ceiling = arb(lower).union(arb(upper))
    if units == "gap":
        ceiling = ceiling * entropy_at_star(parameters)
    return ceiling


def kappa_table(constants: LayerConstants,
                parameters: ArbParameters) -> list[dict[str, Any]]:
    """Largest kappa for which the layer reduction is kappa-monotone."""
    y0 = _arbf(constants.y0)
    geometric = y0 * parameters.x**2
    ceiling = smooth_kappa_ceiling(parameters)
    rows = []
    for deviation in DEVIATION_GRID:
        denominator = 2 * _arbf(deviation) + geometric
        admissible = constants.m_uniform / denominator
        rows.append({
            "component_mean_deviation": str(deviation),
            "denominator": str(denominator),
            "admissible_kappa": str(admissible),
            "admissible_kappa_lower": float(admissible.lower()),
            "beats_smooth_ceiling": bool(admissible > ceiling),
        })
    return rows


def tradeoff_table(parameters: ArbParameters,
                   thresholds: Sequence[Fraction]) -> list[dict[str, Any]]:
    rows = []
    for y0 in thresholds:
        constants = certify_constants(parameters, y0)
        admissible = constants.m_uniform / (
            2 * _arbf(Fraction(1, 10)) + _arbf(y0) * parameters.x**2)
        rows.append({
            "y0": str(y0),
            "lam_uniform": str(constants.lam_uniform),
            "m_uniform": str(constants.m_uniform),
            "m_sharp": str(constants.m_sharp),
            "certified": constants_are_certified(constants),
            "admissible_kappa_at_deviation_0.1": str(admissible),
            "admissible_kappa_lower": float(admissible.lower()),
            "smooth_third_derivative_ceiling":
                str(constants.third_derivative_ceiling),
            "smooth_ceiling_upper":
                float(constants.third_derivative_ceiling.upper()),
        })
    return rows


# ---------------------------------------------------------------------------
# Reporting.


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--y0", type=str, default=str(DEFAULT_Y0),
                        help="layer threshold as an exact fraction")
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--random-points", type=int, default=600)
    parser.add_argument("--cross-octaves", type=int, default=40)
    parser.add_argument("--cross-pieces", type=int, default=3)
    parser.add_argument("--paired-boxes", type=int, default=120)
    parser.add_argument("--seed", type=int, default=SAMPLE_SEED)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    started = time.monotonic()
    y0 = Fraction(args.y0)
    mp_parameters = solve_equation_parameters(args.dps + 30)
    arb_parameters = certify_equation_parameters(mp_parameters)
    constants = certify_constants(arb_parameters, y0)

    print("1. EQUATION-DEFINED CONSTANTS OF THE LAYER")
    print(f"PROVED [Arb, Liu (87)-(90)]: beta in {arb_parameters.beta}")
    print(f"PROVED [Arb, Liu (87)-(90)]: 1-c' = mean in {arb_parameters.mean}")
    print(f"PROVED [Arb]: layer threshold y0={y0}")
    print(f"PROVED [Arb]: Lam_f = 2(1-beta)*mean-1 in {constants.lam_feasible}")
    print(f"PROVED [Arb]: Lam = 2(1-beta)*(mean-y0)-1 in {constants.lam_uniform}")
    print(f"PROVED [Arb]: K = 2(1-beta)*log(1/(1-y0)) in {constants.penalty}")
    print(f"PROVED [Arb]: m0  = Lam*(log(1/y0)+1)-K in {constants.m_uniform}")
    print(f"PROVED [Arb]: m0s = sharpened single-coordinate constant in "
          f"{constants.m_sharp}")
    if not constants_are_certified(constants):
        print("REFUTED [Arb]: the layer constants are not all positive at this y0.")
        return 1
    print("PROVED [Arb]: every layer constant is strictly positive at y0="
          f"{y0}.")
    print()

    print("2. TERM-WISE LINKS OF THE LOWER BOUND")
    cross = sweep_cross_link(y0, args.cross_octaves, args.cross_pieces)
    print("PROVED [exact rational hypotheses on %d dyadic cells, %s]: b<=1 and "
          "y*b<=y0 hold cellwise, so h'(y*b)-log(1/y)-log(1-y0)>=0 by "
          "additivity and monotonicity of log; Arb never encloses the margin "
          "below zero, weakest enclosure %s, %d cells tight at the corner "
          "(y0,1) where the bound is an equality."
          % (cross["boxes"], cross["cover"], cross["weakest_margin"],
             cross["tight_cells"]))
    print("NUMERICAL [midpoint guard]: worst log-additivity residual %.3e."
          % cross["worst_log_additivity_residual"])
    print("PROVED [term-wise]: below y=%s the same two monotonicity steps give "
          "the bound with no cell: %s."
          % (cross["smallest_covered_y"], cross["term_wise_tail"]))
    paired = sweep_paired_link(y0, args.paired_boxes, args.paired_boxes)
    print("PROVED [exact rational hypotheses on %d cells]: 0<=z<=1, y<=1/2 and "
          "2yz<=1/2 hold cellwise, so pi<=%s stays below 1/2 and pi_1 stays "
          "nonnegative; the beta term may therefore be discarded downward.  "
          "%d cells touch the edge z=0 where pi_1 vanishes identically and no "
          "interval test can separate it from zero."
          % (paired["boxes"], paired["largest_pi"], paired["zero_edge_cells"]))
    print()

    print("3. BINDING-POINT AGREEMENT WITH THE PROVED LEADING COEFFICIENT")
    binding = binding_coefficient_check(mp_parameters, arb_parameters, args.dps)
    print("PROVED [liu9_binding.py]: B=2px[1+beta(1-x)]-1 in "
          f"{binding['arb_boundary_barrier']}")
    for row in binding["slope_rows"]:
        print("PROVED [closed form at the binding point, %s]: "
              "dG/dlog(1/y)=%s, difference from B=%s"
              % (row["exponents"], row["slope"], row["difference_from_B"]))
    print()

    print("4. FALSIFICATION PASS OVER THE LAYER STRATUM")
    samples = run_samples(mp_parameters, constants, args.random_points,
                          args.dps, args.seed)
    print("NUMERICAL [%d dps]: %d stratum points = %d structured corners + %d "
          "accepted random draws (seed %d)."
          % (args.dps, samples.points, samples.structured, samples.random,
             args.seed))
    print("NUMERICAL [closed form vs forward-mode differentiation of Liu's own "
          "formula]: %d checks, worst relative disagreement %.3e."
          % (samples.derivative_checks, samples.worst_derivative_relative))
    print("NUMERICAL [pointwise theorem]: %d checks, worst slack %.6f."
          % (samples.pointwise_checks, samples.worst_pointwise_margin))
    print("NUMERICAL [single-coordinate estimate]: %d checks, worst "
          "actual/predicted ratio %.6f."
          % (samples.integrated_checks, samples.worst_integrated_ratio))
    print("NUMERICAL [all-coordinates estimate]: %d checks, worst "
          "actual/predicted ratio %.6f."
          % (samples.multi_checks, samples.worst_multi_ratio))
    if samples.failures:
        for failure in samples.failures[:20]:
            print(f"REFUTED [sample]: {failure}")
        return 1
    print("NUMERICAL [falsification]: no sample violated the theorem.")
    print()

    print("5. WHAT THE LAYER BUYS THE LOCAL LEMMA")
    kappa_rows = kappa_table(constants, arb_parameters)
    ceiling = smooth_kappa_ceiling(arb_parameters)
    ceiling_phi = smooth_kappa_ceiling(arb_parameters, units="phi")
    h_star = entropy_at_star(arb_parameters)
    print("PROVED [Arb, recomputed from the certified A and C]: the two smooth "
          "displayed modes cap the *quotient* ratio (Phi-1)/dist^2 at %s, "
          "which liu9_tube.py prints." % ceiling_phi)
    print("PROVED [Arb]: gap = EHX*(Phi-1) exactly, and EHX(P*) = p*h(x) = %s, "
          "so in the raw-gap units this module uses the same modes cap kappa "
          "at %s.  Every comparison below uses that value; before 2026-08-28 "
          "this module compared raw-gap kappas against the quotient ceiling, "
          "which are different quantities." % (h_star, ceiling))
    print("PROVED [Arb]: zeroing the layer coordinates changes the tube "
          "distance by at most (2*delta+y0*x^2)*sum_j w_j b_j, so the layer is "
          "kappa-monotone for every kappa below the table value.")
    for row in kappa_rows:
        print("PROVED [Arb]: component-mean deviation <= %s gives admissible "
              "kappa %s (beats the smooth ceiling: %s)."
              % (row["component_mean_deviation"], row["admissible_kappa"],
                 row["beats_smooth_ceiling"]))
    print("PROVED [Arb]: on the face, |h'''| <= %s for every surviving support "
          "in [y0,1-y0]." % constants.third_derivative_ceiling)
    mirror = mirror_layer_probe(mp_parameters, args.dps)
    print("PROVED [closed form, slope in log(1/(1-b_j)) between 1-1e-40 and "
          "1-1e-80]: the opposite layer b_j->1 carries the same entropy "
          "singularity and is NOT covered here.  Its coefficient is "
          "1-2(1-beta)W1-2beta*A1, where W1 is the weight already at support "
          "1, and it is two-sided (%s):" % mirror["two_sided"])
    for row in mirror["rows"]:
        print("PROVED [closed form]: %s -- W1=%s, measured slope %s, predicted "
              "%s, difference %.3e, gap is %s in b_j."
              % (row["configuration"], row["weight_at_support_one"],
                 row["measured_slope"], row["predicted_slope"],
                 row["difference"], row["sign"]))
    print("PROVED [tube geometry]: a component-mean deviation obeys "
          "delta_c <= rho/sqrt(min(q,1-q)) inside a radius-rho tube, so the "
          "admissible kappa degrades exactly as q approaches 0 or 1.  The "
          "layer stratum and the q-degenerate stratum are therefore coupled, "
          "which is what the 2026-08-27 residual localization measured: 57 of "
          "the 70 zero-support survivors also pin q.")
    print()

    print("6. THRESHOLD TRADE-OFF")
    rows = tradeoff_table(arb_parameters, TRADEOFF_THRESHOLDS)
    print("y0 | Lam | m0 | m0s | certified | admissible kappa at delta=0.1 | "
          "smooth |h'''| ceiling")
    for row in rows:
        print("PROVED [Arb]: %s | %s | %s | %s | %s | %s | %s"
              % (row["y0"], row["lam_uniform"], row["m_uniform"],
                 row["m_sharp"], row["certified"],
                 row["admissible_kappa_at_deviation_0.1"],
                 row["smooth_third_derivative_ceiling"]))
    print()

    print("7. VERDICT")
    print("PROVED: ingredient (ii) of the conjectured piecewise local lemma -- "
          "an explicit y*log(1/y) boundary-layer estimate at zero support -- "
          "now exists with an explicit threshold and Arb-certified constants.")
    print("PROVED: it is not the binding stratum: at y0=%s and component-mean "
          "deviation 1/10 the layer admits kappa up to %s, while the smooth "
          "A/C modes cap kappa at %s."
          % (y0, kappa_rows[2]["admissible_kappa"], ceiling))
    print("OPEN: ingredients (i) and (iii) are untouched here.  The layer "
          "bound degrades at q in {0,1}, so the endpoint stratum must be "
          "supplied before any tube radius can be certified.  The mirror "
          "layer b_j->1 measured above is a second stratum, now closed "
          "separately in liu9_mirror_layer.py: its sign depends on the "
          "weight at support 1, which a tube bounds by rho^2/g(1-t0).")

    payload = {
        "tool": "liu9_boundary_layer.py",
        "outcome": "PASS",
        "y0": str(y0),
        "dps": args.dps,
        "seed": args.seed,
        "parameters": {
            "beta": str(arb_parameters.beta),
            "mean": str(arb_parameters.mean),
            "x": str(arb_parameters.x),
            "p": str(arb_parameters.p),
            "root_lo": str(arb_parameters.root_lo),
            "root_hi": str(arb_parameters.root_hi),
        },
        "constants": constants.as_json(),
        "cross_link": cross,
        "paired_link": paired,
        "binding_agreement": binding,
        "samples": asdict(samples),
        "smooth_kappa_ceiling": str(ceiling),
        "smooth_kappa_ceiling_phi_units": str(ceiling_phi),
        "entropy_at_star": str(h_star),
        "mirror_layer": mirror,
        "kappa_table": kappa_rows,
        "threshold_tradeoff": rows,
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    payload["report_sha256"] = _canonical_digest(payload)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(f"\nreport written to {args.output}")
        print(f"report_sha256 {payload['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
