#!/usr/bin/env python3
"""Centered-Hessian certificate for Liu's active-mean smooth chart.

The direct interval Hessian in :mod:`liu9_smooth_chart` loses about 811 times
its ``(s,d)`` box radius to dependency.  This verifier differentiates the same
chart once more and encloses each Hessian entry by its value at ``(x,0)`` plus
a mean-value bound for its variation.  The third-order jet is used only on this
smooth chart; every nonconstant entropy argument is first certified to stay
strictly inside ``(0,1)``.  The zero and mirror strata remain the responsibility
of ``liu9_boundary_layer.py`` and ``liu9_mirror_layer.py``.

For Shannon entropy ``h(u)=-u log(u)-(1-u)log(1-u)``, the derivative used here
is

    h'''(u) = (1-2u)/(u^2 (1-u)^2).

An earlier task draft stated the opposite sign.  The independent finite-
difference gate below caught that sign error; the formula is not trusted merely
because it appears in this source.

Run from the repository root with

    ./.venv/bin/python -I -B uc/liu9_chart_centered.py \
        --output uc/verification/results/liu9-chart-centered.json
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
from datetime import datetime, timezone
from fractions import Fraction
from typing import Any, Optional, Sequence

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
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_boundary_layer import _arbf, entropy_at_star  # noqa: E402
from liu9_objective import evaluate_mpmath  # noqa: E402
from liu9_smooth_chart import (  # noqa: E402
    Jet2 as NaiveJet2,
    distance_squared_jet as naive_distance_squared_jet,
    gap_jet as naive_gap_jet,
    smooth_chart_ceiling,
)

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 80
DEFAULT_Q_CELLS = 16
DEFAULT_BISECTION_STEPS = 24
DEFAULT_THIRD_SUBDIVISIONS = 2
Q_FLOOR = Fraction(1, 4)
RADIUS_GRID = tuple(Fraction(1, 2**power) for power in range(3, 12))
REFERENCE_NAIVE_INFLATION = 811.0
TUBE_CHART_REACH = Fraction(3, 20)
FINITE_DIFFERENCE_DIGITS = 20


# ---------------------------------------------------------------------------
# A third-order jet in the two genuine chart variables (s,d).


class Jet3:
    """Value and all symmetric partial derivatives through order three."""

    __slots__ = (
        "v", "gs", "gd", "hss", "hsd", "hdd",
        "tsss", "tssd", "tsdd", "tddd",
    )

    def __init__(
        self,
        v: arb,
        gs: Optional[arb] = None,
        gd: Optional[arb] = None,
        hss: Optional[arb] = None,
        hsd: Optional[arb] = None,
        hdd: Optional[arb] = None,
        tsss: Optional[arb] = None,
        tssd: Optional[arb] = None,
        tsdd: Optional[arb] = None,
        tddd: Optional[arb] = None,
    ):
        zero = arb(0)
        self.v = v
        self.gs = zero if gs is None else gs
        self.gd = zero if gd is None else gd
        self.hss = zero if hss is None else hss
        self.hsd = zero if hsd is None else hsd
        self.hdd = zero if hdd is None else hdd
        self.tsss = zero if tsss is None else tsss
        self.tssd = zero if tssd is None else tssd
        self.tsdd = zero if tsdd is None else tsdd
        self.tddd = zero if tddd is None else tddd

    @staticmethod
    def constant(value: object) -> "Jet3":
        return Jet3(value if isinstance(value, arb) else arb(value))

    @staticmethod
    def variable(value: arb, which: str) -> "Jet3":
        if which == "s":
            return Jet3(value, gs=arb(1))
        if which == "d":
            return Jet3(value, gd=arb(1))
        raise ValueError("variable must be 's' or 'd'")

    def _coerce(self, other: object) -> "Jet3":
        return other if isinstance(other, Jet3) else Jet3.constant(other)

    def __add__(self, other: object) -> "Jet3":
        rhs = self._coerce(other)
        return Jet3(
            self.v + rhs.v,
            self.gs + rhs.gs,
            self.gd + rhs.gd,
            self.hss + rhs.hss,
            self.hsd + rhs.hsd,
            self.hdd + rhs.hdd,
            self.tsss + rhs.tsss,
            self.tssd + rhs.tssd,
            self.tsdd + rhs.tsdd,
            self.tddd + rhs.tddd,
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet3":
        return Jet3(
            -self.v,
            -self.gs,
            -self.gd,
            -self.hss,
            -self.hsd,
            -self.hdd,
            -self.tsss,
            -self.tssd,
            -self.tsdd,
            -self.tddd,
        )

    def __sub__(self, other: object) -> "Jet3":
        return self + (-self._coerce(other))

    def __rsub__(self, other: object) -> "Jet3":
        return self._coerce(other) + (-self)

    def __mul__(self, other: object) -> "Jet3":
        rhs = self._coerce(other)
        return Jet3(
            self.v * rhs.v,
            self.gs * rhs.v + self.v * rhs.gs,
            self.gd * rhs.v + self.v * rhs.gd,
            self.hss * rhs.v + 2 * self.gs * rhs.gs + self.v * rhs.hss,
            self.hsd * rhs.v + self.gs * rhs.gd + self.gd * rhs.gs
            + self.v * rhs.hsd,
            self.hdd * rhs.v + 2 * self.gd * rhs.gd + self.v * rhs.hdd,
            self.tsss * rhs.v + 3 * self.hss * rhs.gs
            + 3 * self.gs * rhs.hss + self.v * rhs.tsss,
            self.tssd * rhs.v + self.hss * rhs.gd
            + 2 * self.hsd * rhs.gs + 2 * self.gs * rhs.hsd
            + self.gd * rhs.hss + self.v * rhs.tssd,
            self.tsdd * rhs.v + 2 * self.hsd * rhs.gd
            + self.hdd * rhs.gs + self.gs * rhs.hdd
            + 2 * self.gd * rhs.hsd + self.v * rhs.tsdd,
            self.tddd * rhs.v + 3 * self.hdd * rhs.gd
            + 3 * self.gd * rhs.hdd + self.v * rhs.tddd,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "Jet3":
        if self.v.contains(0):
            raise ArithmeticError("jet reciprocal across zero")
        inverse = 1 / self.v
        return compose_unary(
            self,
            inverse,
            -(inverse**2),
            2 * inverse**3,
            -6 * inverse**4,
        )

    def __truediv__(self, other: object) -> "Jet3":
        return self * self._coerce(other).reciprocal()

    def is_constant_endpoint(self) -> bool:
        derivatives = (
            self.gs, self.gd, self.hss, self.hsd, self.hdd,
            self.tsss, self.tssd, self.tsdd, self.tddd,
        )
        return bool(
            all(value == 0 for value in derivatives)
            and (self.v == 0 or self.v == 1)
        )


def compose_unary(
    argument: Jet3,
    value: arb,
    first: arb,
    second: arb,
    third: arb,
) -> Jet3:
    """Compose a scalar C3 function with a two-variable third-order jet.

    The mixed Faa di Bruno terms are written explicitly.  In particular,
    ``ssd`` contains all three second/first pairings
    ``u_ss*u_d + 2*u_sd*u_s`` and ``sdd`` contains their transposes.
    """
    return Jet3(
        value,
        first * argument.gs,
        first * argument.gd,
        second * argument.gs * argument.gs + first * argument.hss,
        second * argument.gs * argument.gd + first * argument.hsd,
        second * argument.gd * argument.gd + first * argument.hdd,
        third * argument.gs**3
        + 3 * second * argument.gs * argument.hss
        + first * argument.tsss,
        third * argument.gs**2 * argument.gd
        + second * (
            argument.hss * argument.gd
            + 2 * argument.hsd * argument.gs
        )
        + first * argument.tssd,
        third * argument.gs * argument.gd**2
        + second * (
            argument.hdd * argument.gs
            + 2 * argument.hsd * argument.gd
        )
        + first * argument.tsdd,
        third * argument.gd**3
        + 3 * second * argument.gd * argument.hdd
        + first * argument.tddd,
    )


def entropy_third_formula(value: arb) -> arb:
    """Point formula for Shannon entropy's third derivative."""
    return (1 - 2 * value) / (value**2 * (1 - value) ** 2)


def entropy_third_enclosure(value: arb) -> arb:
    """A tight monotonicity enclosure of h''' on an Arb interval.

    ``h'''(u)=u^-2-(1-u)^-2`` is strictly decreasing on ``(0,1)``.
    Endpoint evaluation therefore avoids the dependency inflation of directly
    substituting a wide interval into the rational expression.
    """
    if not (value > 0 and value < 1):
        raise ArithmeticError(f"entropy third derivative left (0,1): {value}")
    at_lower = entropy_third_formula(value.lower())
    at_upper = entropy_third_formula(value.upper())
    return at_upper.union(at_lower)


def jet_entropy(argument: Jet3) -> Jet3:
    """Shannon entropy, with endpoints only for identically constant 0 or 1.

    Value and derivative ranges are assembled from endpoint values and the
    sole interior extremum at ``1/2``.  Direct ball substitution into
    ``u*(1-u)`` can cross zero on a wide positive interval even though the
    exact product cannot; the monotonicity form prevents that avoidable
    dependency failure.
    """
    if argument.is_constant_endpoint():
        return Jet3.constant(0)
    if not (argument.v > 0 and argument.v < 1):
        raise ArithmeticError(
            f"entropy jet argument left the open unit interval: {argument.v}"
        )

    def value_at(point: arb) -> arb:
        return -(point * point.log() + (1 - point) * (1 - point).log())

    def first_at(point: arb) -> arb:
        return (1 - point).log() - point.log()

    def second_at(point: arb) -> arb:
        return -1 / point - 1 / (1 - point)

    lower = argument.v.lower()
    upper = argument.v.upper()
    value = value_at(lower).union(value_at(upper))
    first = first_at(upper).union(first_at(lower))
    second = second_at(lower).union(second_at(upper))
    half = _arbf(Fraction(1, 2))
    if argument.v.contains(half):
        value = value.union(value_at(half))
        second = second.union(arb(-4))
    third = entropy_third_enclosure(argument.v)
    return compose_unary(argument, value, first, second, third)

def protocol_diagonal(argument: Jet3) -> Jet3:
    """The within-component protocol at the chart's sole nonzero support.

    On this chart every other support in the same component is identically
    zero, so EHPI's only nonzero entropy term has
    ``pi(u,u)=u^2*(1+(1-u)^2)``.  Treating this as a monotone unary map avoids
    the false interval value above one produced by multiplying four dependent
    copies of a wide ``u`` box.
    """
    if not (argument.v > 0 and argument.v < 1):
        raise ArithmeticError(f"protocol support left (0,1): {argument.v}")

    def function(value: arb) -> arb:
        return value**2 * (1 + (1 - value) ** 2)

    def first_derivative(value: arb) -> arb:
        return 4 * value - 6 * value**2 + 4 * value**3

    def second_derivative(value: arb) -> arb:
        return 4 - 12 * value + 12 * value**2

    lower = argument.v.lower()
    upper = argument.v.upper()
    value = function(lower).union(function(upper))
    first = first_derivative(lower).union(first_derivative(upper))
    second = second_derivative(lower).union(second_derivative(upper))
    if argument.v.contains(_arbf(Fraction(1, 2))):
        second = second.union(arb(1))
    third = -12 + 24 * argument.v
    return compose_unary(argument, value, first, second, third)



# ---------------------------------------------------------------------------
# The exact active-mean chart and the two functionals.


def chart(
    parameters: ArbParameters,
    q: arb,
    s: Jet3,
    d: Jet3,
    r: arb,
) -> tuple[list[Jet3], list[Jet3], list[Jet3]]:
    """Return global weights, supports, and within-component masses."""
    mass_active = Jet3.constant(parameters.mean) / s
    mass_gauge = Jet3.constant(r)
    mass_zero = Jet3.constant(1) - mass_active - mass_gauge
    masses = [mass_active, mass_gauge, mass_zero]
    support_zero = s - Jet3.constant(q) * d
    support_one = s + Jet3.constant(1 - q) * d
    supports = [
        support_zero,
        Jet3.constant(0),
        Jet3.constant(0),
        support_one,
        Jet3.constant(0),
        Jet3.constant(0),
    ]
    weights = [Jet3.constant(1 - q) * masses[index] for index in range(3)]
    weights += [Jet3.constant(q) * masses[index] for index in range(3)]
    return weights, supports, masses


def gap_jet(
    parameters: ArbParameters,
    q: arb,
    s: Jet3,
    d: Jet3,
    r: arb,
) -> Jet3:
    """``(1-beta) EHXY + beta EHPI - EHX`` on the chart."""
    weights, supports, masses = chart(parameters, q, s, d, r)
    ehxy = Jet3.constant(0)
    for left in range(6):
        for right in range(6):
            ehxy += (
                weights[left]
                * weights[right]
                * jet_entropy(supports[left] * supports[right])
            )

    # The component factors are load-bearing and match liu9_objective._formula.
    # Every term involving offsets 1 or 2 has an identically zero protocol
    # support, hence zero entropy.  Keeping only offset 0 is exact, and exposes
    # the repeated-support protocol to the tight unary map above.
    ehpi = Jet3.constant(0)
    for component, component_weight in enumerate((1 - q, q)):
        active_support = supports[3 * component]
        inner = (
            masses[0]
            * masses[0]
            * jet_entropy(protocol_diagonal(active_support))
        )
        ehpi += Jet3.constant(component_weight) * inner

    ehx = Jet3.constant(0)
    for index in range(6):
        ehx += weights[index] * jet_entropy(supports[index])

    beta = Jet3.constant(parameters.beta)
    return (Jet3.constant(1) - beta) * ehxy + beta * ehpi - ehx


def distance_squared_jet(
    parameters: ArbParameters,
    q: arb,
    s: Jet3,
    d: Jet3,
    r: arb,
) -> Jet3:
    """The same tube distance squared as ``liu9_smooth_chart.py``."""
    _, supports, masses = chart(parameters, q, s, d, r)
    total = Jet3.constant(0)
    for component, component_weight in enumerate((1 - q, q)):
        base = 3 * component
        component_mean = Jet3.constant(0)
        spread = Jet3.constant(0)
        for offset in range(3):
            point = supports[base + offset]
            component_mean += masses[offset] * point
            deviation = point * (point - Jet3.constant(parameters.x))
            spread += masses[offset] * deviation * deviation
        centered_mean = component_mean - Jet3.constant(parameters.mean)
        total += Jet3.constant(component_weight) * (
            centered_mean * centered_mean + spread
        )
    return total


# ---------------------------------------------------------------------------
# Domain and centered-form construction.


@dataclass(frozen=True)
class DomainCertificate:
    radius: Fraction
    q_floor: Fraction
    checked: bool
    safe: bool
    support_range: arb
    product_range: arb
    protocol_range: arb
    support_h3_bound: arb
    entropy_argument_h3_bound: arb
    failure: Optional[str] = None

    def as_json(self) -> dict[str, Any]:
        return {
            "radius": str(self.radius),
            "q_floor": str(self.q_floor),
            "checked": self.checked,
            "safe": self.safe,
            "support_range": str(self.support_range),
            "support_lower": float(self.support_range.lower()),
            "support_upper": float(self.support_range.upper()),
            "product_range": str(self.product_range),
            "protocol_range": str(self.protocol_range),
            "support_h3_bound": str(self.support_h3_bound),
            "support_h3_bound_upper": float(self.support_h3_bound.upper()),
            "entropy_argument_h3_bound": str(self.entropy_argument_h3_bound),
            "entropy_argument_h3_bound_upper": float(
                self.entropy_argument_h3_bound.upper()
            ),
            "failure": self.failure,
        }


def _abs_upper(value: arb) -> arb:
    return abs(value).upper()


def _maximum_arb(values: Sequence[arb]) -> arb:
    if not values:
        return arb(0)
    return max(values, key=lambda item: float(item.upper()))


def _h3_abs_bound(value: arb) -> arb:
    return _abs_upper(entropy_third_enclosure(value))


def certify_support_domain(
    parameters: ArbParameters,
    radius: Fraction,
    q_floor: Fraction = Q_FLOOR,
) -> DomainCertificate:
    """Certify all nonconstant entropy arguments before any jet evaluation."""
    if radius <= 0:
        raise ValueError("radius must be positive")
    if not Fraction(0) < q_floor < Fraction(1, 2):
        raise ValueError("q_floor must lie strictly between 0 and 1/2")

    radius_arb = _arbf(radius)
    largest_split_coefficient = _arbf(1 - q_floor)
    support_delta = (1 + largest_split_coefficient) * radius_arb
    support_range = (
        parameters.x - support_delta
    ).union(parameters.x + support_delta)

    empty = arb(0)
    if not (support_range > 0 and support_range < 1):
        return DomainCertificate(
            radius,
            q_floor,
            True,
            False,
            support_range,
            empty,
            empty,
            empty,
            empty,
            "a nonzero chart support is not certified inside (0,1)",
        )

    support_lower = support_range.lower()
    support_upper = support_range.upper()
    product_range = (support_lower * support_lower).union(
        support_upper * support_upper
    )

    def protocol(value: arb) -> arb:
        return value * value * (1 + (1 - value) * (1 - value))

    protocol_range = protocol(support_lower).union(protocol(support_upper))
    argument_ranges = (support_range, product_range, protocol_range)
    if not all(value > 0 and value < 1 for value in argument_ranges):
        return DomainCertificate(
            radius,
            q_floor,
            True,
            False,
            support_range,
            product_range,
            protocol_range,
            empty,
            empty,
            "a nonconstant entropy argument is not certified inside (0,1)",
        )

    support_h3 = _h3_abs_bound(support_range)
    all_h3 = _maximum_arb([_h3_abs_bound(value) for value in argument_ranges])
    return DomainCertificate(
        radius,
        q_floor,
        True,
        True,
        support_range,
        product_range,
        protocol_range,
        support_h3,
        all_h3,
    )


def _require_domain(
    domain: Optional[DomainCertificate],
    radius: Fraction,
    q_floor: Fraction,
) -> DomainCertificate:
    if domain is None or not domain.checked:
        raise ValueError("centered enclosure requires a checked support domain")
    if domain.radius != radius or domain.q_floor != q_floor:
        raise ValueError("support-domain certificate does not match this box")
    if not domain.safe:
        raise ArithmeticError(domain.failure or "unsafe support domain")
    return domain


@dataclass(frozen=True)
class HessianBox:
    ss: arb
    sd: arb
    dd: arb

    def entries(self) -> tuple[arb, arb, arb]:
        return self.ss, self.sd, self.dd


def centered_entry(
    center: arb,
    derivative_s: arb,
    derivative_d: arb,
    radius: Fraction,
) -> arb:
    """Mean-value enclosure on an ``l-infinity`` box of the given radius."""
    gradient_l1_bound = _abs_upper(derivative_s) + _abs_upper(derivative_d)
    variation = _arbf(radius) * gradient_l1_bound
    return center + (-variation).union(variation)


def _centered_hessian(jet_at_center: Jet3, jet_on_box: Jet3,
                      radius: Fraction) -> HessianBox:
    return HessianBox(
        centered_entry(
            jet_at_center.hss, jet_on_box.tsss, jet_on_box.tssd, radius
        ),
        centered_entry(
            jet_at_center.hsd, jet_on_box.tssd, jet_on_box.tsdd, radius
        ),
        centered_entry(
            jet_at_center.hdd, jet_on_box.tsdd, jet_on_box.tddd, radius
        ),
    )


def _jet3_hessian(jet: Jet3) -> HessianBox:
    return HessianBox(jet.hss, jet.hsd, jet.hdd)


def _jet2_hessian(jet: NaiveJet2) -> HessianBox:
    return HessianBox(jet.hss, jet.hsd, jet.hdd)


def _interval_width(value: arb) -> float:
    return float(value.upper()) - float(value.lower())

def _update_supremum(current: arb, candidate: arb) -> arb:
    candidate_upper = _abs_upper(candidate)
    if candidate_upper > current:
        return candidate_upper
    if candidate_upper <= current:
        return current
    return current.union(candidate_upper).upper()


def _third_derivative_suprema(
    parameters: ArbParameters,
    q: arb,
    radius: Fraction,
    subdivisions: int,
) -> tuple[Jet3, Jet3]:
    """Enclose third partials on a finite cover of the full (s,d) box."""
    if subdivisions <= 0:
        raise ValueError("third-derivative subdivisions must be positive")
    gap_bounds = [arb(0) for _ in range(4)]
    distance_bounds = [arb(0) for _ in range(4)]
    for s_index in range(subdivisions):
        s_lower_offset = -radius + 2 * radius * s_index / subdivisions
        s_upper_offset = -radius + 2 * radius * (s_index + 1) / subdivisions
        s_box = (
            parameters.x + _arbf(s_lower_offset)
        ).union(parameters.x + _arbf(s_upper_offset))
        for d_index in range(subdivisions):
            d_lower = -radius + 2 * radius * d_index / subdivisions
            d_upper = -radius + 2 * radius * (d_index + 1) / subdivisions
            d_box = _arbf(d_lower).union(_arbf(d_upper))
            gap = gap_jet(
                parameters,
                q,
                Jet3.variable(s_box, "s"),
                Jet3.variable(d_box, "d"),
                arb(0),
            )
            distance = distance_squared_jet(
                parameters,
                q,
                Jet3.variable(s_box, "s"),
                Jet3.variable(d_box, "d"),
                arb(0),
            )
            for index, value in enumerate(
                (gap.tsss, gap.tssd, gap.tsdd, gap.tddd)
            ):
                gap_bounds[index] = _update_supremum(
                    gap_bounds[index], value
                )
            for index, value in enumerate(
                (
                    distance.tsss,
                    distance.tssd,
                    distance.tsdd,
                    distance.tddd,
                )
            ):
                distance_bounds[index] = _update_supremum(
                    distance_bounds[index], value
                )
    return (
        Jet3(
            arb(0),
            tsss=gap_bounds[0],
            tssd=gap_bounds[1],
            tsdd=gap_bounds[2],
            tddd=gap_bounds[3],
        ),
        Jet3(
            arb(0),
            tsss=distance_bounds[0],
            tssd=distance_bounds[1],
            tsdd=distance_bounds[2],
            tddd=distance_bounds[3],
        ),
    )


@dataclass(frozen=True)
class CellEnclosure:
    q_lower: Fraction
    q_upper: Fraction
    gap: HessianBox
    distance: HessianBox
    naive_gap: Optional[HessianBox]
    naive_distance: HessianBox
    naive_gap_failure: Optional[str]
    center_containment_checks: int
    full_overlap_checks: int


def build_cell_enclosure(
    parameters: ArbParameters,
    radius: Fraction,
    q_lower: Fraction,
    q_upper: Fraction,
    domain: Optional[DomainCertificate],
    q_floor: Fraction = Q_FLOOR,
    third_subdivisions: int = DEFAULT_THIRD_SUBDIVISIONS,
) -> CellEnclosure:
    """Build one q-cell after validating its matching support certificate."""
    _require_domain(domain, radius, q_floor)
    if not q_floor <= q_lower < q_upper <= 1 - q_floor:
        raise ValueError("q cell lies outside the certified q range")

    radius_arb = _arbf(radius)
    q = _arbf(q_lower).union(_arbf(q_upper))
    s_box = (parameters.x - radius_arb).union(parameters.x + radius_arb)
    d_box = (-radius_arb).union(radius_arb)

    gap_box, distance_box = _third_derivative_suprema(
        parameters, q, radius, third_subdivisions
    )
    gap_center = gap_jet(
        parameters,
        q,
        Jet3.variable(parameters.x, "s"),
        Jet3.variable(arb(0), "d"),
        arb(0),
    )
    distance_center = distance_squared_jet(
        parameters,
        q,
        Jet3.variable(parameters.x, "s"),
        Jet3.variable(arb(0), "d"),
        arb(0),
    )
    centered_gap = _centered_hessian(gap_center, gap_box, radius)
    centered_distance = _centered_hessian(
        distance_center, distance_box, radius
    )

    naive_gap_failure = None
    try:
        naive_gap_box = naive_gap_jet(
            parameters,
            q,
            NaiveJet2.variable(s_box, "s"),
            NaiveJet2.variable(d_box, "d"),
            arb(0),
        )
        naive_gap: Optional[HessianBox] = _jet2_hessian(naive_gap_box)
        if not all(
            math.isfinite(float(value.lower()))
            and math.isfinite(float(value.upper()))
            for value in naive_gap.entries()
        ):
            naive_gap = None
            naive_gap_failure = (
                "liu9_smooth_chart.py returned a non-finite Hessian enclosure"
            )
    except ArithmeticError as error:
        # This is the direct liu9_smooth_chart.py comparison, not part of the
        # centered proof.  At wide radii its dependent protocol interval may
        # leave (0,1) even though the semantic support-domain check is safe.
        naive_gap = None
        naive_gap_failure = str(error)
    naive_distance_box = naive_distance_squared_jet(
        parameters,
        q,
        NaiveJet2.variable(s_box, "s"),
        NaiveJet2.variable(d_box, "d"),
        arb(0),
    )
    naive_gap_center = naive_gap_jet(
        parameters,
        q,
        NaiveJet2.variable(parameters.x, "s"),
        NaiveJet2.variable(arb(0), "d"),
        arb(0),
    )
    naive_distance_center = naive_distance_squared_jet(
        parameters,
        q,
        NaiveJet2.variable(parameters.x, "s"),
        NaiveJet2.variable(arb(0), "d"),
        arb(0),
    )
    naive_distance = _jet2_hessian(naive_distance_box)

    containment_checks = 0
    overlap_checks = 0
    comparisons = (
        (centered_gap, _jet2_hessian(naive_gap_center), naive_gap),
        (
            centered_distance,
            _jet2_hessian(naive_distance_center),
            naive_distance,
        ),
    )
    for centered, naive_center, naive_full in comparisons:
        for index, (centered_value, center_value) in enumerate(zip(
            centered.entries(), naive_center.entries()
        )):
            containment_checks += 1
            if not centered_value.contains(center_value):
                raise AssertionError(
                    "centered Hessian lost the independent naive center value"
                )
            if naive_full is not None:
                overlap_checks += 1
                if not centered_value.overlaps(naive_full.entries()[index]):
                    raise AssertionError(
                        "centered and naive whole-box Hessians do not intersect"
                    )

    return CellEnclosure(
        q_lower,
        q_upper,
        centered_gap,
        centered_distance,
        naive_gap,
        naive_distance,
        naive_gap_failure,
        containment_checks,
        overlap_checks,
    )


# ---------------------------------------------------------------------------
# Conservative PSD pencil and rational-kappa bisection.


def pencil_is_psd(
    gap: HessianBox,
    distance: HessianBox,
    kappa: arb,
) -> bool:
    """Certify a 2x2 interval pencil by M11 and determinant lower bounds."""
    m11 = gap.ss - kappa * distance.ss
    m22 = gap.dd - kappa * distance.dd
    m12 = gap.sd - kappa * distance.sd
    if not m11 >= 0:
        return False
    off_diagonal = _abs_upper(m12)
    determinant_lower = (m11 * m22).lower() - off_diagonal * off_diagonal
    return bool(determinant_lower >= 0)


def candidate_holds(
    cells: Sequence[CellEnclosure],
    candidate: Fraction,
    ceiling: arb,
) -> tuple[bool, Optional[str], Optional[tuple[Fraction, Fraction]]]:
    """Check the ceiling and every q cell for one exact rational kappa."""
    value = _arbf(candidate)
    if not value <= ceiling.lower():
        return False, "above_ceiling", None
    for cell in cells:
        if not pencil_is_psd(cell.gap, cell.distance, value):
            return False, "pencil", (cell.q_lower, cell.q_upper)
    return True, None, None


@dataclass(frozen=True)
class RadiusCertificate:
    radius: Fraction
    q_floor: Fraction
    q_cells: int
    domain: DomainCertificate
    kappa: Fraction
    kappa_resolution: Fraction
    ceiling: arb
    centered_hss_width: float
    naive_hss_width: Optional[float]
    centered_inflation: float
    naive_inflation: Optional[float]
    naive_failure: Optional[str]
    third_subdivisions: int
    center_containment_checks: int
    full_overlap_checks: int
    next_failure: Optional[str]
    next_failing_q_cell: Optional[tuple[Fraction, Fraction]]

    def certified(self) -> bool:
        return bool(
            self.domain.checked
            and self.domain.safe
            and self.kappa > 0
            and _arbf(self.kappa) <= self.ceiling.lower()
        )

    def as_json(self) -> dict[str, Any]:
        return {
            "radius": str(self.radius),
            "radius_float": float(self.radius),
            "q_floor": str(self.q_floor),
            "q_cells": self.q_cells,
            "q_cell_width": float((1 - 2 * self.q_floor) / self.q_cells),
            "third_derivative_subdivisions_per_axis": self.third_subdivisions,
            "support_domain": self.domain.as_json(),
            "kappa": str(self.kappa),
            "kappa_float": float(self.kappa),
            "kappa_resolution": str(self.kappa_resolution),
            "kappa_ceiling": str(self.ceiling),
            "kappa_ceiling_upper": float(self.ceiling.upper()),
            "fraction_of_ceiling": (
                float(self.kappa) / float(self.ceiling.lower())
                if self.kappa > 0
                else 0.0
            ),
            "certified": self.certified(),
            "centered_hss_width": self.centered_hss_width,
            "naive_hss_width": self.naive_hss_width,
            "centered_inflation_constant": self.centered_inflation,
            "naive_inflation_constant": self.naive_inflation,
            "naive_failure": self.naive_failure,
            "centered_beats_naive": (
                self.centered_hss_width < self.naive_hss_width
                if self.naive_hss_width is not None
                else None
            ),
            "center_containment_checks": self.center_containment_checks,
            "full_overlap_checks": self.full_overlap_checks,
            "next_failure": self.next_failure,
            "next_failing_q_cell": (
                [str(value) for value in self.next_failing_q_cell]
                if self.next_failing_q_cell
                else None
            ),
        }


def _q_cell_count(radius: Fraction, q_floor: Fraction,
                  minimum: int) -> int:
    span = 1 - 2 * q_floor
    sized = (span.numerator * radius.denominator +
             span.denominator * radius.numerator - 1) // (
                 span.denominator * radius.numerator
             )
    return max(minimum, sized)


def certify_radius(
    parameters: ArbParameters,
    radius: Fraction,
    q_floor: Fraction = Q_FLOOR,
    minimum_q_cells: int = DEFAULT_Q_CELLS,
    bisection_steps: int = DEFAULT_BISECTION_STEPS,
    third_subdivisions: int = DEFAULT_THIRD_SUBDIVISIONS,
) -> RadiusCertificate:
    """Certify one radius, including its support-domain prerequisite."""
    if minimum_q_cells <= 0:
        raise ValueError("minimum_q_cells must be positive")
    if bisection_steps <= 0:
        raise ValueError("bisection_steps must be positive")
    if third_subdivisions <= 0:
        raise ValueError("third_subdivisions must be positive")
    domain = certify_support_domain(parameters, radius, q_floor)
    ceiling = smooth_chart_ceiling(parameters)
    if not domain.safe:
        return RadiusCertificate(
            radius,
            q_floor,
            0,
            domain,
            Fraction(0),
            Fraction(1, 2**bisection_steps),
            ceiling,
            math.inf,
            None,
            math.inf,
            None,
            "out_of_domain",
            third_subdivisions,
            0,
            0,
            "out_of_domain",
            None,
        )

    q_cells = _q_cell_count(radius, q_floor, minimum_q_cells)
    span = 1 - 2 * q_floor
    cells = []
    for index in range(q_cells):
        q_lower = q_floor + span * index / q_cells
        q_upper = q_floor + span * (index + 1) / q_cells
        cells.append(
            build_cell_enclosure(
                parameters,
                radius,
                q_lower,
                q_upper,
                domain,
                q_floor,
                third_subdivisions,
            )
        )

    centered_width = max(_interval_width(cell.gap.ss) for cell in cells)
    naive_boxes = [cell.naive_gap for cell in cells]
    if all(value is not None for value in naive_boxes):
        naive_width: Optional[float] = max(
            _interval_width(value.ss)
            for value in naive_boxes
            if value is not None
        )
        naive_failure = None
    else:
        naive_width = None
        naive_failure = next(
            cell.naive_gap_failure
            for cell in cells
            if cell.naive_gap_failure is not None
        )
    radius_float = float(radius)
    naive_inflation = (
        naive_width / radius_float if naive_width is not None else None
    )
    containment_checks = sum(
        cell.center_containment_checks for cell in cells
    )
    overlap_checks = sum(cell.full_overlap_checks for cell in cells)

    zero_holds, zero_reason, zero_cell = candidate_holds(
        cells, Fraction(0), ceiling
    )
    if not zero_holds:
        return RadiusCertificate(
            radius,
            q_floor,
            q_cells,
            domain,
            Fraction(0),
            Fraction(1, 2**bisection_steps),
            ceiling,
            centered_width,
            naive_width,
            centered_width / radius_float,
            naive_inflation,
            naive_failure,
            third_subdivisions,
            containment_checks,
            overlap_checks,
            zero_reason,
            zero_cell,
        )

    low = Fraction(0)
    high = Fraction(1, 2)
    failure_reason: Optional[str] = None
    failure_cell: Optional[tuple[Fraction, Fraction]] = None
    for _ in range(bisection_steps):
        middle = (low + high) / 2
        holds, reason, failing_cell = candidate_holds(cells, middle, ceiling)
        if holds:
            low = middle
        else:
            high = middle
            failure_reason = reason
            failure_cell = failing_cell

    holds, _, _ = candidate_holds(cells, low, ceiling)
    if not holds:
        raise AssertionError("reported rational kappa failed its own certificate")
    if low > 0 and not _arbf(low) <= ceiling.lower():
        raise AssertionError("reported kappa exceeds the smooth-chart ceiling")
    return RadiusCertificate(
        radius,
        q_floor,
        q_cells,
        domain,
        low,
        high - low,
        ceiling,
        centered_width,
        naive_width,
        centered_width / radius_float,
        naive_inflation,
        naive_failure,
        third_subdivisions,
        containment_checks,
        overlap_checks,
        failure_reason,
        failure_cell,
    )


# ---------------------------------------------------------------------------
# Independent high-precision finite-difference gate.


def _mp_gap(mp_parameters: Any, q: mpmath.mpf, s: mpmath.mpf,
            d: mpmath.mpf, dps: int) -> mpmath.mpf:
    active_mass = mp_parameters.mean / s
    values = (
        active_mass,
        mpmath.mpf(0),
        q,
        s - q * d,
        mpmath.mpf(0),
        mpmath.mpf(0),
        s + (1 - q) * d,
        mpmath.mpf(0),
        mpmath.mpf(0),
    )
    terms = evaluate_mpmath(values, mp_parameters.beta, dps=dps)
    return terms.numerator - terms.ehx


def _third_difference(
    function: Any,
    s: mpmath.mpf,
    d: mpmath.mpf,
    step: mpmath.mpf,
    derivative: str,
) -> mpmath.mpf:
    """Second-order centered finite-difference stencils of total order three."""
    h = step
    if derivative == "sss":
        return (
            function(s + 2 * h, d)
            - 2 * function(s + h, d)
            + 2 * function(s - h, d)
            - function(s - 2 * h, d)
        ) / (2 * h**3)
    if derivative == "ddd":
        return (
            function(s, d + 2 * h)
            - 2 * function(s, d + h)
            + 2 * function(s, d - h)
            - function(s, d - 2 * h)
        ) / (2 * h**3)
    if derivative == "ssd":
        def second_s(at_d: mpmath.mpf) -> mpmath.mpf:
            return (
                function(s + h, at_d)
                - 2 * function(s, at_d)
                + function(s - h, at_d)
            )
        return (second_s(d + h) - second_s(d - h)) / (2 * h**3)
    if derivative == "sdd":
        def second_d(at_s: mpmath.mpf) -> mpmath.mpf:
            return (
                function(at_s, d + h)
                - 2 * function(at_s, d)
                + function(at_s, d - h)
            )
        return (second_d(s + h) - second_d(s - h)) / (2 * h**3)
    raise ValueError(f"unknown third derivative {derivative}")


def _richardson_third_difference(
    function: Any,
    s: mpmath.mpf,
    d: mpmath.mpf,
    step: mpmath.mpf,
    derivative: str,
) -> mpmath.mpf:
    coarse = _third_difference(function, s, d, step, derivative)
    fine = _third_difference(function, s, d, step / 2, derivative)
    return (4 * fine - coarse) / 3


def _arb_midpoint_mp(value: arb) -> mpmath.mpf:
    text = value.mid().str(100)
    token = text.split()[0].lstrip("[")
    return mpmath.mpf(token)


def finite_difference_check(
    mp_parameters: Any,
    parameters: ArbParameters,
    dps: int = 120,
) -> dict[str, Any]:
    """Check all four third partials at three interior chart points."""
    rows = []
    worst_relative = mpmath.mpf(0)
    with mpmath.workdps(dps):
        points = (
            (Fraction(1, 3), mp_parameters.x + mpmath.mpf(1) / 100,
             mpmath.mpf(1) / 200),
            (Fraction(2, 5), mp_parameters.x - mpmath.mpf(1) / 120,
             -mpmath.mpf(1) / 180),
            (Fraction(7, 10), mp_parameters.x + mpmath.mpf(1) / 256,
             -mpmath.mpf(1) / 300),
        )
        step = mpmath.mpf("1e-7")
        for q_fraction, s_mp, d_mp in points:
            q_mp = (
                mpmath.mpf(q_fraction.numerator) / q_fraction.denominator
            )
            cache: dict[tuple[str, str], mpmath.mpf] = {}

            def evaluate(s_value: mpmath.mpf,
                         d_value: mpmath.mpf) -> mpmath.mpf:
                key = (mpmath.nstr(s_value, dps), mpmath.nstr(d_value, dps))
                if key not in cache:
                    cache[key] = _mp_gap(
                        mp_parameters, q_mp, s_value, d_value, dps
                    )
                return cache[key]

            q_arb = _arbf(q_fraction)
            s_arb = arb(mpmath.nstr(s_mp, dps))
            d_arb = arb(mpmath.nstr(d_mp, dps))
            jet = gap_jet(
                parameters,
                q_arb,
                Jet3.variable(s_arb, "s"),
                Jet3.variable(d_arb, "d"),
                arb(0),
            )
            actuals = {
                "sss": jet.tsss,
                "ssd": jet.tssd,
                "sdd": jet.tsdd,
                "ddd": jet.tddd,
            }
            derivative_rows = {}
            for name, actual_ball in actuals.items():
                estimate = _richardson_third_difference(
                    evaluate, s_mp, d_mp, step, name
                )
                actual = _arb_midpoint_mp(actual_ball)
                relative = abs(actual - estimate) / max(
                    mpmath.mpf(1), abs(actual), abs(estimate)
                )
                worst_relative = max(worst_relative, relative)
                derivative_rows[name] = {
                    "jet": mpmath.nstr(actual, 50),
                    "finite_difference": mpmath.nstr(estimate, 50),
                    "relative_error": float(relative),
                }
            rows.append({
                "q": str(q_fraction),
                "s": mpmath.nstr(s_mp, 50),
                "d": mpmath.nstr(d_mp, 50),
                "derivatives": derivative_rows,
            })

    tolerance = mpmath.mpf(10) ** (-FINITE_DIFFERENCE_DIGITS)
    if not worst_relative < tolerance:
        raise AssertionError(
            "third-order jet failed the 20-digit finite-difference gate: "
            f"worst relative error {worst_relative}"
        )
    return {
        "points": rows,
        "derivative_checks": len(rows) * 4,
        "required_digits": FINITE_DIFFERENCE_DIGITS,
        "worst_relative_error": float(worst_relative),
        "passed": True,
        "sign_selected": "(1-2u)/(u^2(1-u)^2)",
        "rejected_sign": "(2u-1)/(u^2(1-u)^2)",
    }


# ---------------------------------------------------------------------------
# Reporting.


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def tube_cell_estimate(radius: Fraction,
                       reach: Fraction = TUBE_CHART_REACH) -> float:
    """Continuous covering count used by the original 9.7e7 estimate."""
    return (float(reach / radius) ** 2) * float(Fraction(1, 2) / radius)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--q-cells", type=int, default=DEFAULT_Q_CELLS,
                        help="minimum cells; cells are refined to width <= radius")
    parser.add_argument("--bisection-steps", type=int,
                        default=DEFAULT_BISECTION_STEPS)
    parser.add_argument("--third-subdivisions", type=int,
                        default=DEFAULT_THIRD_SUBDIVISIONS,
                        help="subdivisions per (s,d) axis for third partials")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    started = time.monotonic()
    mp_parameters = solve_equation_parameters(args.dps + 50)
    parameters = certify_equation_parameters(mp_parameters)
    ceiling = smooth_chart_ceiling(parameters)

    print("1. THIRD-ORDER JET SANITY GATE")
    finite_difference = finite_difference_check(
        mp_parameters, parameters, max(120, args.dps + 40)
    )
    print(
        "MACHINE VERIFIED [Richardson finite differences]: all %d third "
        "partials at three interior points agree to at least %d digits; worst "
        "relative error %.3e."
        % (
            finite_difference["derivative_checks"],
            finite_difference["required_digits"],
            finite_difference["worst_relative_error"],
        )
    )
    print(
        "MACHINE VERIFIED [sign mutation]: the gate selects h'''(u)="
        "(1-2u)/(u^2(1-u)^2) and rejects the opposite sign from the earlier "
        "task draft."
    )
    print()

    print("2. SUPPORT DOMAIN AND ENTROPY THIRD-DERIVATIVE BOUNDS")
    certificates = []
    for radius in RADIUS_GRID:
        certificate = certify_radius(
            parameters,
            radius,
            Q_FLOOR,
            args.q_cells,
            args.bisection_steps,
            args.third_subdivisions,
        )
        certificates.append(certificate)
        domain = certificate.domain
        if domain.safe:
            print(
                "PROVED [Arb domain]: radius %s, q in [1/4,3/4] -> supports "
                "in [%s,%s], support |h'''| <= %s, and |h'''| <= %s over "
                "every nonconstant entropy argument."
                % (
                    radius,
                    domain.support_range.lower(),
                    domain.support_range.upper(),
                    domain.support_h3_bound,
                    domain.entropy_argument_h3_bound,
                )
            )
        else:
            print(
                "REFUTED [domain]: radius %s was not evaluated: %s"
                % (radius, domain.failure)
            )
    print()

    print("3. CENTERED FORM VERSUS THE NAIVE INTERVAL HESSIAN")
    for certificate in certificates:
        if not certificate.domain.safe:
            continue
        if certificate.naive_hss_width is None:
            print(
                "MACHINE VERIFIED [Arb, %d q cells]: radius %s -> worst "
                "H_ss width centered %.9g (constant %.6g); the direct "
                "liu9_smooth_chart.py enclosure is OUT OF DOMAIN: %s."
                % (
                    certificate.q_cells,
                    certificate.radius,
                    certificate.centered_hss_width,
                    certificate.centered_inflation,
                    certificate.naive_failure,
                )
            )
        else:
            print(
                "MACHINE VERIFIED [Arb, same %d q cells]: radius %s -> worst "
                "H_ss width centered %.9g (constant %.6g), naive %.9g "
                "(constant %.6g), centered beats naive: %s."
                % (
                    certificate.q_cells,
                    certificate.radius,
                    certificate.centered_hss_width,
                    certificate.centered_inflation,
                    certificate.naive_hss_width,
                    certificate.naive_inflation,
                    certificate.centered_hss_width
                    < certificate.naive_hss_width,
                )
            )
    total_containment = sum(
        certificate.center_containment_checks for certificate in certificates
    )
    total_overlap = sum(
        certificate.full_overlap_checks for certificate in certificates
    )
    print(
        "MACHINE VERIFIED [independent Jet2 intersection gate]: %d centered "
        "entries contain the independently evaluated naive centre entry, and "
        "%d whole-box centered/naive pairs overlap."
        % (total_containment, total_overlap)
    )
    smallest = certificates[-1]
    print(
        "MACHINE VERIFIED [measured inflation]: at radius %s the centered "
        "constant is %.6g versus the same-run naive %.6g and the previously "
        "reported approximately %.0f."
        % (
            smallest.radius,
            smallest.centered_inflation,
            smallest.naive_inflation,
            REFERENCE_NAIVE_INFLATION,
        )
    )
    print()

    print("4. CENTERED PENCIL CERTIFICATES")
    print(
        "PROVED [Arb centre pencil]: every reported kappa is constrained below "
        "the q-free raw-gap ceiling %s."
        % ceiling
    )
    for certificate in certificates:
        if not certificate.domain.safe:
            continue
        status = "PROVED" if certificate.certified() else "REFUTED"
        print(
            "%s [Arb centered pencil]: radius %s, q cells %d -> kappa %s "
            "(%.2f%% of ceiling), next failure %s%s."
            % (
                status,
                certificate.radius,
                certificate.q_cells,
                certificate.kappa,
                100 * (
                    float(certificate.kappa) / float(ceiling.lower())
                    if certificate.kappa > 0
                    else 0
                ),
                certificate.next_failure,
                (
                    " at q cell "
                    + str(certificate.next_failing_q_cell)
                    if certificate.next_failing_q_cell
                    else ""
                ),
            )
        )
    print()

    positive = [item for item in certificates if item.certified()]
    largest = max(positive, key=lambda item: item.radius) if positive else None
    if largest is None:
        cells_needed = math.inf
        print("5. VERDICT")
        print(
            "REFUTED [this centered enclosure]: no radius on the requested "
            "grid certified a positive kappa."
        )
    else:
        cells_needed = tube_cell_estimate(largest.radius)
        print("5. VERDICT")
        print(
            "PROVED [Arb centered pencil]: the largest requested grid radius "
            "with positive kappa is %s, with kappa >= %s for every q in "
            "[1/4,3/4]."
            % (largest.radius, largest.kappa)
        )
        print(
            "COMPUTATIONAL EVIDENCE [same covering model as the old 9.7e7 "
            "count]: reaching |s-x| about 0.15 with cells of radius %s needs "
            "about %.6g cells in (s,d,q)."
            % (largest.radius, cells_needed)
        )
    print(
        "OPEN: this certificate concerns only the active-mean smooth chart; it "
        "does not replace the separately certified zero-support, mirror, or "
        "q-degenerate strata and does not extend the chart to all nine "
        "variables."
    )

    radius_rows = [item.as_json() for item in certificates]
    payload: dict[str, Any] = {
        "tool": "liu9_chart_centered.py",
        "outcome": "PASS",
        "dps": args.dps,
        "q_floor": str(Q_FLOOR),
        "minimum_q_cells": args.q_cells,
        "bisection_steps": args.bisection_steps,
        "third_derivative_subdivisions_per_axis": args.third_subdivisions,
        "parameters": {
            "x": str(parameters.x),
            "p": str(parameters.p),
            "mean": str(parameters.mean),
            "beta": str(parameters.beta),
            "entropy_at_star": str(entropy_at_star(parameters)),
        },
        "finite_difference_gate": finite_difference,
        "smooth_chart_ceiling": str(ceiling),
        "smooth_chart_ceiling_upper": float(ceiling.upper()),
        "radius_table": radius_rows,
        "inflation_comparison": {
            "reference_old_constant": REFERENCE_NAIVE_INFLATION,
            "same_run_naive_constant": smallest.naive_inflation,
            "new_centered_constant": smallest.centered_inflation,
            "improvement_factor": (
                smallest.naive_inflation / smallest.centered_inflation
                if (
                    smallest.naive_inflation is not None
                    and smallest.centered_inflation > 0
                )
                else None
            ),
            "beats_naive": (
                smallest.centered_hss_width < smallest.naive_hss_width
                if smallest.naive_hss_width is not None
                else None
            ),
        },
        "sanity_gates": {
            "center_containment_checks": total_containment,
            "full_overlap_checks": total_overlap,
            "finite_difference_digits": FINITE_DIFFERENCE_DIGITS,
            "all_kappas_below_ceiling": all(
                item.kappa == 0 or _arbf(item.kappa) <= ceiling.lower()
                for item in certificates
            ),
            "all_certified_radii_checked_support": all(
                not item.certified()
                or (item.domain.checked and item.domain.safe)
                for item in certificates
            ),
        },
        "largest_positive_certificate": (
            largest.as_json() if largest else None
        ),
        "tube_reach": float(TUBE_CHART_REACH),
        "cells_for_tube_reach": cells_needed,
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
