#!/usr/bin/env python3
"""q-weighted endpoint estimate for Liu's Hypothesis 2.

PROVED (definitions).  Let ``x,p,beta`` be the equation-defined constants,
``m=p*x``, and ``P*=p delta_x+(1-p)delta_0``.  For a displacement ``d`` put
``Q_d=p delta_(x+d)+(1-p)delta_0``.  The exact mean-active path

    x0=x-q*d,  x1=x+(1-q)*d,
    V(q,d)=(p,0,q,x0,0,0,x1,0,0)

has mean ``m`` and endpoints ``(P*,Q_d)`` at ``q=0`` and ``(Q_-d,P*)`` at
``q=1``.  Write ``G=gap=(1-beta)EHXY+beta*EHPI-EHX`` and

    delta(Q)^2=(M(Q)-m)^2+integral [y(y-x)]^2 dQ(y),
    dist^2=(1-q)delta(P0)^2+q delta(P1)^2.

MACHINE VERIFIED THEOREM (tube-restricted endpoint chart).  Fix
``R=1/4``, a smooth cutoff ``0<eps_sm<=1/32``, an endpoint threshold ``q_*``,
and a tube radius ``rho`` satisfying the two certified seam inequalities
printed by this program.  For every ``d`` with
``eps_sm<=|d|<=R`` and every ``q in [0,q_*]`` (or, by the exact component
swap, ``q in [1-q_*,1]``),

    G(V(q,d)) >= kappa_end * dist(V(q,d))^2,

with the explicit raw-gap ``kappa_end`` in the report.  At the exact endpoint
the inactive law has zero weight in both expressions, so this is a limiting
q-weighted statement, not evaluation at an interior surrogate.  On the inner
core ``|d|<=eps_sm``, a separately factored cubic estimate is uniform for all
``q in [0,1]``.  The default headline is ``eps_sm=1/32``; all dependence on
``eps_sm`` is retained in the certificate.

PROVED (endpoint first variation).  If ``D`` is exactly the normalized
functional certified nonnegative by ``liu9_endpoint.py``, then

    d/dq G(V(q,d)) at q=0 = H* D(Q_d),   H*=p h(x).

On the restricted annulus this module certifies
``H*D(Q_d)>=k_D*delta(Q_d)^2`` and hence an explicit positive ``D>=d0``.
It does not infer a positive unrestricted infimum from the sign theorem.

PROVED (powers of q).  In reduced coordinates the quadratic objective term
and squared distance both contain ``q(1-q)``.  The pure-d cubic at ``d=0`` has
weight

    (1-q)(-q)^3+q(1-q)^3=q(1-q)(1-2q),

and on the whole inner chart ``G'''_d=q(1-q)B(q,d)``.  Arb bounds ``B`` before
any division, on dyadic cells whose entropy arguments stay in
``[1/32,31/32]``; there ``|h'''|<1022.94``.

PROVED (symmetry and gauges).  Liu's transcription is invariant under the
simultaneous map ``(q,P0,P1)->(1-q,P1,P0)`` with the shared masses unchanged.
It is not invariant under ``q->1-q`` with component supports fixed.  At
``q=0`` all P1 terms have zero objective and distance weight, and at ``q=1``
all P0 terms do; the endpoint ``d`` direction is therefore a joint gauge.

CONDITIONAL (integrated seam).  If an interior smooth-chart theorem is valid
for ``q in [q_int,1-q_int]`` and ``|d|<=eps_sm``, the regimes meet without a
q-gap exactly when ``q_int<=q_*``.  Inside ``dist<=rho``, the inequality
``rho^2<=p^2 eps_sm^2 q_*(1-q_*)`` forces every interior-q pure-d point into
that smooth cutoff.  The program proves this implication and reports whether
the supplied constants satisfy it; it does not claim an unavailable smooth
chart radius.

OPEN.  Arbitrary inactive three-atom laws outside the displayed ``Q_d`` chart,
and pure-d laws with ``|d|>R``, remain assigned to the complement/boundary
strata.  In particular, the unrestricted infimum of ``D(Q)/delta(Q)^2`` is
not certified here.

Run from the ``math`` directory with one-core variables set:

    ./.venv/bin/python -I -B uc/liu9_qdegenerate.py \
        --output uc/verification/results/liu9-qdegenerate.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from dataclasses import asdict, dataclass
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
    MPParameters,
    certify_curvatures,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import (  # noqa: E402
    _ScalarOps,
    _formula,
    evaluate_mpmath,
    h_arb,
)

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 120
Y0 = Fraction(1, 32)
ENDPOINT_RADIUS = Fraction(1, 4)
DEFAULT_SMOOTH_CUTOFF = Fraction(1, 32)
DEFAULT_Q_STAR = Fraction(1, 4096)
DEFAULT_INTERIOR_Q_MIN = Fraction(1, 4096)
DEFAULT_TUBE_RADIUS = Fraction(1, 4096)
ENDPOINT_RATIO = Fraction(1, 4)
ENDPOINT_KAPPA = Fraction(1, 20)
CORE_KAPPA = Fraction(1, 3)
Q_REMAINDER_CAP = Fraction(1, 32)
DEFAULT_ENDPOINT_PIECES = 64
DEFAULT_Q_CELLS = 32
DEFAULT_D_CELLS = 128
EXPERIMENT_D = "1e-20"
EXPERIMENT_Q = (
    "1e-1", "1e-3", "1e-6", "1e-9",
    "0.9", "0.999", "0.999999", "0.999999999",
)


def _arbf(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _ball(lower: Fraction, upper: Fraction) -> arb:
    if lower > upper:
        raise ValueError("reversed interval")
    return _arbf(lower).union(_arbf(upper))


def _abs_upper(value: arb) -> arb:
    upper = value.upper()
    negative_lower = -value.lower()
    return negative_lower if upper < negative_lower else upper


def _entropy(value: Any) -> Any:
    if isinstance(value, arb):
        if not (value > 0 and value < 1):
            raise ArithmeticError(
                f"entropy interval touches a singular endpoint: {value}")
        return h_arb(value)
    if value == 0 or value == 1:
        return value - value
    if not 0 < value < 1:
        raise ArithmeticError("entropy argument outside [0,1]")
    return -(value * mpmath.log(value)
             + (1 - value) * mpmath.log(1 - value))


def _hp(value: Any) -> Any:
    if isinstance(value, arb):
        return (1 - value).log() - value.log()
    return mpmath.log(1 - value) - mpmath.log(value)


def _hpp(value: Any) -> Any:
    return -1 / (value * (1 - value))


def _hthird(value: Any) -> Any:
    return (1 - 2 * value) / (value**2 * (1 - value)**2)


def _pi(left: Any, right: Any) -> Any:
    return left * right * (1 + (1 - left) * (1 - right))

def _constant_like(value: Any, constant: int) -> Any:
    if isinstance(value, arb):
        return arb(constant)
    if isinstance(value, mpmath.mpf):
        return mpmath.mpf(constant)
    return constant


def _zero_like(value: Any) -> Any:
    return _constant_like(value, 0)


def _one_like(value: Any) -> Any:
    return _constant_like(value, 1)


class Jet3:
    """Value and first three ordinary derivatives for forward mode."""

    __slots__ = ("value", "first", "second", "third")

    def __init__(self, value: Any, first: Any = 0, second: Any = 0,
                 third: Any = 0) -> None:
        self.value = value
        zero = _zero_like(value)
        self.first = zero + first
        self.second = zero + second
        self.third = zero + third

    @staticmethod
    def _coerce(other: Any, like: Any) -> "Jet3":
        return (other if isinstance(other, Jet3)
                else Jet3(_zero_like(like) + other))

    def __add__(self, other: Any) -> "Jet3":
        other = self._coerce(other, self.value)
        return Jet3(
            self.value + other.value,
            self.first + other.first,
            self.second + other.second,
            self.third + other.third,
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet3":
        return Jet3(-self.value, -self.first, -self.second, -self.third)

    def __sub__(self, other: Any) -> "Jet3":
        return self + (-self._coerce(other, self.value))

    def __rsub__(self, other: Any) -> "Jet3":
        return self._coerce(other, self.value) - self

    def __mul__(self, other: Any) -> "Jet3":
        other = self._coerce(other, self.value)
        return Jet3(
            self.value * other.value,
            self.first * other.value + self.value * other.first,
            self.second * other.value + 2 * self.first * other.first
            + self.value * other.second,
            self.third * other.value + 3 * self.second * other.first
            + 3 * self.first * other.second + self.value * other.third,
        )

    __rmul__ = __mul__


class _GapJetOps(_ScalarOps):
    def quotient(self, numerator: Any, denominator: Any) -> None:
        return None


def entropy_jet(argument: Jet3) -> Jet3:
    value = argument.value
    if bool(value == 0 or value == 1):
        if not (argument.first == 0 and argument.second == 0
                and argument.third == 0):
            raise ArithmeticError("nonconstant entropy jet at an endpoint")
        return Jet3(value - value)
    entropy = _entropy(value)
    hp = _hp(value)
    hpp = _hpp(value)
    hthird = _hthird(value)
    return Jet3(
        entropy,
        hp * argument.first,
        hpp * argument.first**2 + hp * argument.second,
        hthird * argument.first**3
        + 3 * hpp * argument.first * argument.second
        + hp * argument.third,
    )


def _gap_from_jets(values: Sequence[Jet3], beta: Any, one: Any) -> Jet3:
    terms = _formula(
        values,
        Jet3(beta),
        _GapJetOps(entropy_jet, Jet3(one)),
    )
    return terms.numerator - terms.ehx


def formula_d_jet(q: Any, d: Any, parameters: Any) -> Jet3:
    one = _one_like(q)
    zero = _zero_like(q)
    qbar = one - q
    left = parameters.x - q * d
    right = parameters.x + qbar * d
    return _gap_from_jets((
        Jet3(parameters.p), Jet3(zero), Jet3(q),
        Jet3(left, -q), Jet3(zero), Jet3(zero),
        Jet3(right, qbar), Jet3(zero), Jet3(zero),
    ), parameters.beta, one)


def formula_q_jet(q: Any, d: Any, parameters: Any) -> Jet3:
    one = _one_like(q)
    zero = _zero_like(q)
    q_jet = Jet3(q, 1)
    left = Jet3(parameters.x) - q_jet * d
    right = Jet3(parameters.x) + (1 - q_jet) * d
    return _gap_from_jets((
        Jet3(parameters.p), Jet3(zero), q_jet,
        left, Jet3(zero), Jet3(zero),
        right, Jet3(zero), Jet3(zero),
    ), parameters.beta, one)


def _diag_j_third(value: Any) -> Any:
    jet = Jet3(value, 1)
    return entropy_jet(jet * jet).third


def _diag_k_third(value: Any) -> Any:
    jet = Jet3(value, 1)
    return entropy_jet(_pi(jet, jet)).third


def _ell_third(value: Any, parameters: Any) -> Any:
    return (parameters.beta * parameters.p**2 * _diag_k_third(value)
            - parameters.p * _hthird(value))


def factored_third_quotient(q: Any, d: Any, parameters: Any) -> Any:
    """B(q,d) in G'''_d=q(1-q)B, without dividing by q."""
    one = _one_like(q)
    qbar = one - q
    left = parameters.x - q * d
    right = parameters.x + qbar * d
    cross = entropy_jet(Jet3(left, -q) * Jet3(right, qbar)).third
    j_part = (
        -qbar * q**2 * _diag_j_third(left)
        + 2 * cross
        + q * qbar**2 * _diag_j_third(right)
    )
    ell_part = (
        -q**2 * _ell_third(left, parameters)
        + qbar**2 * _ell_third(right, parameters)
    )
    return (1 - parameters.beta) * parameters.p**2 * j_part + ell_part


def cubic_form_at_star(parameters: Any) -> Any:
    x_variable = Jet3(parameters.x, 1)
    j111 = entropy_jet(x_variable * Jet3(parameters.x)).third
    return (_ell_third(parameters.x, parameters)
            + 2 * (1 - parameters.beta) * parameters.p**2 * j111)


def pure_values_mp(q: mpmath.mpf, d: mpmath.mpf,
                   parameters: MPParameters) -> tuple[mpmath.mpf, ...]:
    one = mpmath.mpf(1)
    zero = mpmath.mpf(0)
    return (
        parameters.p, zero, q,
        parameters.x - q * d, zero, zero,
        parameters.x + (one - q) * d, zero, zero,
    )


def pure_gap_mp(q: mpmath.mpf, d: mpmath.mpf,
                parameters: MPParameters, dps: int = DEFAULT_DPS) -> mpmath.mpf:
    terms = evaluate_mpmath(
        pure_values_mp(q, d, parameters), parameters.beta, dps=dps)
    return terms.numerator - terms.ehx


def pure_gap_reduced_mp(q: mpmath.mpf, d: mpmath.mpf,
                        parameters: MPParameters) -> mpmath.mpf:
    qbar = 1 - q
    left = parameters.x - q * d
    right = parameters.x + qbar * d
    p = parameters.p
    ehxy = p**2 * (
        qbar**2 * _entropy(left**2)
        + 2 * qbar * q * _entropy(left * right)
        + q**2 * _entropy(right**2)
    )
    ehpi = p**2 * (
        qbar * _entropy(_pi(left, left))
        + q * _entropy(_pi(right, right))
    )
    ehx = p * (qbar * _entropy(left) + q * _entropy(right))
    return (1 - parameters.beta) * ehxy + parameters.beta * ehpi - ehx


def pure_distance2_mp(q: mpmath.mpf, d: mpmath.mpf,
                      parameters: MPParameters) -> mpmath.mpf:
    qbar = 1 - q
    left = parameters.x - q * d
    right = parameters.x + qbar * d
    mean = parameters.p * parameters.x
    delta0 = ((parameters.p * left - mean)**2
              + parameters.p * (left * (left - parameters.x))**2)
    delta1 = ((parameters.p * right - mean)**2
              + parameters.p * (right * (right - parameters.x))**2)
    return qbar * delta0 + q * delta1


def distance_factor(q: Any, d: Any, parameters: Any) -> Any:
    """Exact B with dist^2=q(1-q)d^2 B."""
    qbar = 1 - q
    return parameters.p**2 + parameters.p * (
        parameters.x**2
        + 2 * parameters.x * (1 - 2 * q) * d
        + (q**3 + qbar**3) * d**2
    )

def _core(parameters: Any) -> Any:
    if hasattr(parameters, "core"):
        return parameters.core
    x = parameters.x
    return (
        (1 - parameters.beta) * _entropy(x * x)
        + parameters.beta * _entropy(_pi(x, x))
    ) / (x * _entropy(x))


def endpoint_raw_first_variation(d: Any, parameters: Any) -> Any:
    """H*D(Q_d), independently reduced from liu9_endpoint.py's definition."""
    x = parameters.x
    p = parameters.p
    y = x + d
    h_star = p * _entropy(x)
    return (
        2 * (1 - parameters.beta) * p**2
        * (_entropy(x * y) - _entropy(x * x))
        + parameters.beta * p**2
        * (_entropy(_pi(y, y)) - _entropy(_pi(x, x)))
        - p * (_entropy(y) - _entropy(x))
        - h_star * _core(parameters) * p * d
    )


def endpoint_raw_first_derivative(d: Any, parameters: Any) -> Any:
    x = parameters.x
    p = parameters.p
    y = x + d
    y_jet = Jet3(y, 1)
    pi_prime = _pi(y_jet, y_jet).first
    h_star = p * _entropy(x)
    return (
        2 * (1 - parameters.beta) * p**2 * x * _hp(x * y)
        + parameters.beta * p**2 * pi_prime * _hp(_pi(y, y))
        - p * _hp(y)
        - h_star * _core(parameters) * p
    )


def endpoint_delta2(d: Any, parameters: Any) -> Any:
    y = parameters.x + d
    return d**2 * (parameters.p**2 + parameters.p * y**2)


def endpoint_delta2_derivative(d: Any, parameters: Any) -> Any:
    y = parameters.x + d
    return (
        2 * d * (parameters.p**2 + parameters.p * y**2)
        + 2 * parameters.p * y * d**2
    )


def cubic_weight(q: Fraction, power: int = 3) -> Fraction:
    q = Fraction(q)
    return (1 - q) * (-q)**power + q * (1 - q)**power


def regime_cover_is_complete(endpoint_max: Fraction,
                             interior_min: Fraction) -> bool:
    endpoint_max = Fraction(endpoint_max)
    interior_min = Fraction(interior_min)
    return bool(
        0 <= endpoint_max <= Fraction(1, 2)
        and 0 <= interior_min <= Fraction(1, 2)
        and interior_min <= endpoint_max
    )

def seam_geometry_holds(parameters: ArbParameters, smooth_cutoff: Fraction,
                        q_star: Fraction, tube_radius: Fraction) -> bool:
    smooth_cutoff = Fraction(smooth_cutoff)
    q_star = Fraction(q_star)
    tube_radius = Fraction(tube_radius)
    if not (smooth_cutoff > 0 and 0 < q_star < Fraction(1, 2)
            and tube_radius > 0):
        return False
    right_lower = (
        parameters.p**2 * _arbf(smooth_cutoff)**2
        * _arbf(q_star) * (1 - _arbf(q_star))).lower()
    return bool(_arbf(tube_radius)**2 <= right_lower)


def _is_power_of_two(value: int) -> bool:
    return value > 0 and not value & (value - 1)


def validate_domain(parameters: ArbParameters, smooth_cutoff: Fraction,
                    endpoint_radius: Fraction = ENDPOINT_RADIUS) -> dict[str, str]:
    smooth_cutoff = Fraction(smooth_cutoff)
    endpoint_radius = Fraction(endpoint_radius)
    if not 0 < smooth_cutoff <= Fraction(1, 32):
        raise ValueError("smooth cutoff must satisfy 0<eps_sm<=1/32")
    if not smooth_cutoff < endpoint_radius <= Fraction(1, 4):
        raise ValueError("endpoint radius must satisfy eps_sm<R<=1/4")
    actual_lo = parameters.root_lo - endpoint_radius
    actual_hi = parameters.root_hi + endpoint_radius
    declared_lo = Fraction(7, 16)
    declared_hi = Fraction(31, 32)
    if actual_lo < declared_lo or actual_hi > declared_hi:
        raise ValueError("endpoint radius leaves the certified support domain")
    entropy_lo = declared_lo**2
    pi_hi = declared_hi**2 * (1 + (1 - declared_hi)**2)
    entropy_hi = max(declared_hi, pi_hi)
    if entropy_lo < Y0 or entropy_hi > 1 - Y0:
        raise AssertionError("entropy arguments leave [y0,1-y0]")
    return {
        "smooth_cutoff": str(smooth_cutoff),
        "endpoint_radius": str(endpoint_radius),
        "actual_support_lower": str(actual_lo),
        "actual_support_upper": str(actual_hi),
        "declared_support_lower": str(declared_lo),
        "declared_support_upper": str(declared_hi),
        "entropy_argument_lower": str(entropy_lo),
        "entropy_argument_upper": str(entropy_hi),
        "y0": str(Y0),
    }


def _radial_cells(inner: Fraction, outer: Fraction,
                  pieces: int) -> list[tuple[Fraction, Fraction]]:
    if not _is_power_of_two(pieces):
        raise ValueError("radial pieces must be a power of two")
    cells: list[tuple[Fraction, Fraction]] = []
    start = inner
    while start < outer:
        stop = min(2 * start, outer)
        width = (stop - start) / pieces
        for index in range(pieces):
            lower = start + index * width
            upper = start + (index + 1) * width
            cells.append((lower, upper))
            cells.append((-upper, -lower))
        start = stop
    return cells


@dataclass(frozen=True)
class EndpointSweep:
    boxes: int
    radial_pieces: int
    cover: str
    boundary_points: int
    weakest_margin: arb

    def as_json(self) -> dict[str, Any]:
        return {
            "boxes": self.boxes,
            "radial_pieces": self.radial_pieces,
            "cover": self.cover,
            "boundary_points": self.boundary_points,
            "weakest_margin": str(self.weakest_margin),
            "weakest_margin_lower": float(self.weakest_margin.lower()),
        }


def certify_endpoint_ratio(parameters: ArbParameters, smooth_cutoff: Fraction,
                           endpoint_radius: Fraction = ENDPOINT_RADIUS,
                           pieces: int = DEFAULT_ENDPOINT_PIECES,
                           ratio: Fraction = ENDPOINT_RATIO) -> EndpointSweep:
    """Certify H*D(Q_d)>=ratio*delta(Q_d)^2 on the endpoint annulus."""
    smooth_cutoff = Fraction(smooth_cutoff)
    endpoint_radius = Fraction(endpoint_radius)
    ratio = Fraction(ratio)
    validate_domain(parameters, smooth_cutoff, endpoint_radius)
    cells = _radial_cells(smooth_cutoff, endpoint_radius, pieces)
    weakest: Optional[arb] = None
    for lower, upper in cells:
        midpoint = (lower + upper) / 2
        cell = _ball(lower, upper)
        delta = _ball(lower - midpoint, upper - midpoint)
        point = (
            endpoint_raw_first_variation(_arbf(midpoint), parameters)
            - _arbf(ratio) * endpoint_delta2(_arbf(midpoint), parameters)
        )
        derivative = (
            endpoint_raw_first_derivative(cell, parameters)
            - _arbf(ratio) * endpoint_delta2_derivative(cell, parameters)
        )
        enclosure = point + derivative * delta
        if not enclosure.is_finite() or not enclosure > 0:
            raise AssertionError(
                "endpoint ratio failed on exact cell [%s,%s]: %s"
                % (lower, upper, enclosure))
        lower_bound = enclosure.lower()
        if weakest is None or lower_bound < weakest:
            weakest = lower_bound
    assert weakest is not None
    return EndpointSweep(
        boxes=len(cells),
        radial_pieces=pieces,
        cover=(f"[-{endpoint_radius},-{smooth_cutoff}] union "
               f"[{smooth_cutoff},{endpoint_radius}]"),
        boundary_points=4,
        weakest_margin=weakest,
    )


@dataclass(frozen=True)
class CubicSweep:
    q_cells: int
    d_cells: int
    boxes: int
    cover: str
    boundary_cells: int
    corner_points: int
    bound: arb

    def as_json(self) -> dict[str, Any]:
        return {
            "q_cells": self.q_cells,
            "d_cells": self.d_cells,
            "boxes": self.boxes,
            "cover": self.cover,
            "boundary_cells": self.boundary_cells,
            "corner_points": self.corner_points,
            "bound": str(self.bound),
            "bound_upper": float(self.bound.upper()),
        }


def certify_cubic_sweep(parameters: ArbParameters, smooth_cutoff: Fraction,
                        q_cells: int = DEFAULT_Q_CELLS,
                        d_cells: int = DEFAULT_Q_CELLS) -> CubicSweep:
    smooth_cutoff = Fraction(smooth_cutoff)
    validate_domain(parameters, smooth_cutoff)
    if not _is_power_of_two(q_cells) or not _is_power_of_two(d_cells):
        raise ValueError("cubic cell counts must be powers of two")
    maximum = arb(0)
    boundary = 0
    for q_index in range(q_cells):
        q_lo = Fraction(q_index, q_cells)
        q_hi = Fraction(q_index + 1, q_cells)
        q = _ball(q_lo, q_hi)
        for d_index in range(d_cells):
            d_lo = -smooth_cutoff + 2 * smooth_cutoff * Fraction(
                d_index, d_cells)
            d_hi = -smooth_cutoff + 2 * smooth_cutoff * Fraction(
                d_index + 1, d_cells)
            enclosure = factored_third_quotient(
                q, _ball(d_lo, d_hi), parameters)
            if not enclosure.is_finite():
                raise ArithmeticError("non-finite cubic enclosure")
            bound = _abs_upper(enclosure)
            if maximum < bound:
                maximum = bound
            if (q_index in (0, q_cells - 1)
                    or d_index in (0, d_cells - 1)):
                boundary += 1
    return CubicSweep(
        q_cells=q_cells,
        d_cells=d_cells,
        boxes=q_cells * d_cells,
        cover=f"[0,1] x [-{smooth_cutoff},{smooth_cutoff}]",
        boundary_cells=boundary,
        corner_points=4,
        bound=arb(0).union(maximum.upper()),
    )


@dataclass(frozen=True)
class RemainderSweep:
    q_cells: int
    d_cells: int
    boxes: int
    cover: str
    boundary_cells: int
    corner_points: int
    bound: arb

    def as_json(self) -> dict[str, Any]:
        return {
            "q_cells": self.q_cells,
            "d_cells": self.d_cells,
            "boxes": self.boxes,
            "cover": self.cover,
            "boundary_cells": self.boundary_cells,
            "corner_points": self.corner_points,
            "bound": str(self.bound),
            "bound_upper": float(self.bound.upper()),
        }


def certify_q_remainder(parameters: ArbParameters,
                        endpoint_radius: Fraction = ENDPOINT_RADIUS,
                        q_cap: Fraction = Q_REMAINDER_CAP,
                        q_cells: int = DEFAULT_Q_CELLS,
                        d_cells: int = DEFAULT_D_CELLS) -> RemainderSweep:
    if not _is_power_of_two(q_cells) or not _is_power_of_two(d_cells):
        raise ValueError("remainder cell counts must be powers of two")
    maximum = arb(0)
    boundary = 0
    for q_index in range(q_cells):
        q_lo = q_cap * Fraction(q_index, q_cells)
        q_hi = q_cap * Fraction(q_index + 1, q_cells)
        q = _ball(q_lo, q_hi)
        for d_index in range(d_cells):
            d_lo = -endpoint_radius + 2 * endpoint_radius * Fraction(
                d_index, d_cells)
            d_hi = -endpoint_radius + 2 * endpoint_radius * Fraction(
                d_index + 1, d_cells)
            enclosure = formula_q_jet(
                q, _ball(d_lo, d_hi), parameters).second
            if not enclosure.is_finite():
                raise ArithmeticError(
                    f"non-finite q-remainder enclosure at q={q_lo,q_hi}, "
                    f"d={d_lo,d_hi}: {enclosure}")
            bound = _abs_upper(enclosure)
            if maximum < bound:
                maximum = bound
            if (q_index in (0, q_cells - 1)
                    or d_index in (0, d_cells - 1)):
                boundary += 1
    return RemainderSweep(
        q_cells=q_cells,
        d_cells=d_cells,
        boxes=q_cells * d_cells,
        cover=f"[0,{q_cap}] x [-{endpoint_radius},{endpoint_radius}]",
        boundary_cells=boundary,
        corner_points=4,
        bound=arb(0).union(maximum.upper()),
    )


@dataclass(frozen=True)
class QDegenerateCertificate:
    smooth_cutoff: Fraction
    endpoint_radius: Fraction
    q_star: Fraction
    interior_q_min: Fraction
    tube_radius: Fraction
    hthird_ceiling: arb
    entropy_at_star: arb
    raw_second_coefficient: arb
    cubic_form: arb
    core_gap_factor_lower: arb
    core_distance_factor_upper: arb
    core_ratio_lower: arb
    endpoint_raw_d0: arb
    endpoint_normalized_d0: arb
    q_star_ceiling: arb
    seam_radius_squared_ceiling: arb
    endpoint_kappa_lower: arb
    endpoint_sweep: EndpointSweep
    cubic_sweep: CubicSweep
    remainder_sweep: RemainderSweep

    def as_json(self) -> dict[str, Any]:
        return {
            "smooth_cutoff": str(self.smooth_cutoff),
            "endpoint_radius": str(self.endpoint_radius),
            "q_star": str(self.q_star),
            "interior_q_min": str(self.interior_q_min),
            "tube_radius": str(self.tube_radius),
            "hthird_ceiling": str(self.hthird_ceiling),
            "entropy_at_star": str(self.entropy_at_star),
            "raw_second_coefficient": str(self.raw_second_coefficient),
            "cubic_form_at_star": str(self.cubic_form),
            "core_gap_factor_lower": str(self.core_gap_factor_lower),
            "core_distance_factor_upper": str(self.core_distance_factor_upper),
            "core_ratio_lower": str(self.core_ratio_lower),
            "endpoint_raw_d0": str(self.endpoint_raw_d0),
            "endpoint_normalized_d0": str(self.endpoint_normalized_d0),
            "q_star_ceiling": str(self.q_star_ceiling),
            "seam_radius_squared_ceiling": str(
                self.seam_radius_squared_ceiling),
            "endpoint_kappa_lower": str(self.endpoint_kappa_lower),
            "endpoint_sweep": self.endpoint_sweep.as_json(),
            "cubic_sweep": self.cubic_sweep.as_json(),
            "remainder_sweep": self.remainder_sweep.as_json(),
            "certified": self.certified(),
        }

    def certified(self) -> bool:
        return bool(
            self.hthird_ceiling < arb("1022.94")
            and self.raw_second_coefficient > 0
            and self.core_gap_factor_lower > 0
            and self.core_ratio_lower > _arbf(CORE_KAPPA)
            and self.endpoint_raw_d0 > 0
            and self.endpoint_normalized_d0 > 0
            and _arbf(self.q_star) <= self.q_star_ceiling
            and _arbf(self.tube_radius)**2
            <= self.seam_radius_squared_ceiling
            and self.endpoint_kappa_lower > _arbf(ENDPOINT_KAPPA)
            and regime_cover_is_complete(
                self.q_star, self.interior_q_min)
        )


def certify(parameters: ArbParameters,
            smooth_cutoff: Fraction = DEFAULT_SMOOTH_CUTOFF,
            q_star: Fraction = DEFAULT_Q_STAR,
            interior_q_min: Fraction = DEFAULT_INTERIOR_Q_MIN,
            tube_radius: Fraction = DEFAULT_TUBE_RADIUS,
            endpoint_pieces: int = DEFAULT_ENDPOINT_PIECES,
            q_cells: int = DEFAULT_Q_CELLS,
            d_cells: int = DEFAULT_D_CELLS) -> QDegenerateCertificate:
    smooth_cutoff = Fraction(smooth_cutoff)
    q_star = Fraction(q_star)
    interior_q_min = Fraction(interior_q_min)
    tube_radius = Fraction(tube_radius)
    validate_domain(parameters, smooth_cutoff)
    if not 0 < q_star <= Q_REMAINDER_CAP:
        raise ValueError("q_* must lie in (0,1/32]")
    if tube_radius <= 0:
        raise ValueError("tube radius must be positive")

    endpoint = certify_endpoint_ratio(
        parameters, smooth_cutoff, pieces=endpoint_pieces)
    cubic = certify_cubic_sweep(
        parameters, smooth_cutoff, q_cells=q_cells, d_cells=q_cells)
    remainder = certify_q_remainder(
        parameters, q_cells=q_cells, d_cells=d_cells)

    y0 = _arbf(Y0)
    hthird_ceiling = (
        (1 - 2 * y0) / (y0**2 * (1 - y0)**2)).upper()
    h_star = parameters.p * h_arb(parameters.x)
    raw_second = h_star * certify_curvatures(parameters).split_unit
    core_gap = (
        raw_second.lower() / 2
        - cubic.bound.upper() * _arbf(smooth_cutoff) / 6).lower()
    core_distance = (
        parameters.p**2
        + parameters.p * (parameters.x + _arbf(smooth_cutoff))**2).upper()
    core_ratio = (core_gap / core_distance).lower()

    raw_d0 = (
        _arbf(ENDPOINT_RATIO) * parameters.p**2
        * _arbf(smooth_cutoff)**2).lower()
    normalized_d0 = (raw_d0 / h_star.upper()).lower()
    q_ceiling = (raw_d0 / remainder.bound.upper()).lower()
    seam_squared = (
        parameters.p**2 * _arbf(smooth_cutoff)**2
        * _arbf(q_star) * (1 - _arbf(q_star))).lower()
    endpoint_distance_upper = (
        parameters.p**2
        + parameters.p * (parameters.x + _arbf(ENDPOINT_RADIUS))**2).upper()
    endpoint_kappa = (
        _arbf(ENDPOINT_RATIO) * parameters.p**2
        / (2 * endpoint_distance_upper)).lower()

    certificate = QDegenerateCertificate(
        smooth_cutoff=smooth_cutoff,
        endpoint_radius=ENDPOINT_RADIUS,
        q_star=q_star,
        interior_q_min=interior_q_min,
        tube_radius=tube_radius,
        hthird_ceiling=hthird_ceiling,
        entropy_at_star=h_star,
        raw_second_coefficient=raw_second,
        cubic_form=cubic_form_at_star(parameters),
        core_gap_factor_lower=core_gap,
        core_distance_factor_upper=core_distance,
        core_ratio_lower=core_ratio,
        endpoint_raw_d0=raw_d0,
        endpoint_normalized_d0=normalized_d0,
        q_star_ceiling=q_ceiling,
        seam_radius_squared_ceiling=seam_squared,
        endpoint_kappa_lower=endpoint_kappa,
        endpoint_sweep=endpoint,
        cubic_sweep=cubic,
        remainder_sweep=remainder,
    )
    if not certificate.certified():
        raise AssertionError("the supplied endpoint/seam constants are not certified")
    return certificate


@dataclass(frozen=True)
class ExperimentRow:
    q: str
    one_minus_q: str
    d: str
    ratio: str


def run_experiment(parameters: MPParameters, dps: int) -> tuple[ExperimentRow, ...]:
    rows = []
    with mpmath.workdps(dps):
        d = mpmath.mpf(EXPERIMENT_D)
        for q_text in EXPERIMENT_Q:
            q = mpmath.mpf(q_text)
            ratio = pure_gap_mp(q, d, parameters, dps) / pure_distance2_mp(
                q, d, parameters)
            rows.append(ExperimentRow(
                q=q_text,
                one_minus_q=mpmath.nstr(1 - q, 20),
                d=EXPERIMENT_D,
                ratio=mpmath.nstr(ratio, 50),
            ))
    return tuple(rows)


def autodiff_crosscheck(parameters: MPParameters, dps: int) -> dict[str, Any]:
    q_values = tuple(Fraction(value) for value in (
        0, Fraction(1, 4096), Fraction(1, 4), Fraction(1, 2),
        Fraction(3, 4), Fraction(4095, 4096), 1,
    ))
    d_values = (-Fraction(1, 4), -Fraction(1, 32), Fraction(0),
                Fraction(1, 32), Fraction(1, 4))
    worst_value = mpmath.mpf(0)
    worst_d_third = mpmath.mpf(0)
    worst_endpoint = mpmath.mpf(0)
    worst_distance = mpmath.mpf(0)
    checks = 0
    with mpmath.workdps(dps):
        for q_fraction in q_values:
            q = mpmath.mpf(q_fraction.numerator) / q_fraction.denominator
            for d_fraction in d_values:
                d = mpmath.mpf(d_fraction.numerator) / d_fraction.denominator
                d_jet = formula_d_jet(q, d, parameters)
                reduced = pure_gap_reduced_mp(q, d, parameters)
                factored = q * (1 - q) * factored_third_quotient(
                    q, d, parameters)
                worst_value = max(
                    worst_value,
                    abs(d_jet.value - reduced)
                    / max(abs(d_jet.value), abs(reduced), mpmath.mpf(1)),
                )
                worst_d_third = max(
                    worst_d_third,
                    abs(d_jet.third - factored)
                    / max(abs(d_jet.third), abs(factored), mpmath.mpf(1)),
                )
                if q == 0:
                    q_jet = formula_q_jet(q, d, parameters)
                    endpoint = endpoint_raw_first_variation(d, parameters)
                    worst_endpoint = max(
                        worst_endpoint,
                        abs(q_jet.first - endpoint)
                        / max(abs(q_jet.first), abs(endpoint), mpmath.mpf(1)),
                    )
                if 0 < q < 1 and d != 0:
                    direct_distance = pure_distance2_mp(q, d, parameters)
                    reduced_distance = (
                        q * (1 - q) * d**2
                        * distance_factor(q, d, parameters))
                    worst_distance = max(
                        worst_distance,
                        abs(direct_distance - reduced_distance)
                        / max(abs(direct_distance), mpmath.mpf(1)),
                    )
                checks += 1
    tolerance = mpmath.mpf(10) ** (-(dps - 25))
    if max(worst_value, worst_d_third, worst_endpoint,
           worst_distance) > tolerance:
        raise AssertionError("forward-mode cross-check failed")
    return {
        "checks": checks,
        "dps": dps,
        "worst_value_relative": float(worst_value),
        "worst_d_third_relative": float(worst_d_third),
        "worst_endpoint_first_relative": float(worst_endpoint),
        "worst_distance_relative": float(worst_distance),
    }


def _swap_components(values: Sequence[mpmath.mpf]
                     ) -> tuple[mpmath.mpf, ...]:
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    return (a1, a2, 1 - q, b1, b3, b5, b0, b2, b4)


def symmetry_crosscheck(parameters: MPParameters, dps: int) -> dict[str, Any]:
    with mpmath.workdps(dps):
        mp = mpmath.mpf
        values = (
            mp(3) / 10, mp(1) / 5, mp(1) / 5,
            mp(3) / 5, mp(7) / 10, mp(4) / 5,
            mp(13) / 20, mp(3) / 4, mp(17) / 20,
        )
        original = evaluate_mpmath(values, parameters.beta, dps=dps)
        swapped = evaluate_mpmath(
            _swap_components(values), parameters.beta, dps=dps)
        q_only_values = list(values)
        q_only_values[2] = 1 - q_only_values[2]
        q_only = evaluate_mpmath(
            tuple(q_only_values), parameters.beta, dps=dps)
        original_gap = original.numerator - original.ehx
        full_residual = abs(
            original_gap - (swapped.numerator - swapped.ehx))
        q_only_difference = abs(
            original_gap - (q_only.numerator - q_only.ehx))
        pure_residual = mpmath.mpf(0)
        for q_text in ("0.001", "0.2", "0.7", "0.999"):
            q = mp(q_text)
            d = mp(1) / 32
            pure_residual = max(
                pure_residual,
                abs(pure_gap_mp(q, d, parameters, dps)
                    - pure_gap_mp(1 - q, -d, parameters, dps)),
            )
        tolerance = mp(10) ** (-(dps - 20))
        if full_residual > tolerance or pure_residual > tolerance:
            raise AssertionError("component-swap symmetry failed")
        if q_only_difference < mp("1e-8"):
            raise AssertionError("q-only mutation was not detected")
        return {
            "dps": dps,
            "full_component_swap_residual": float(full_residual),
            "pure_q_d_swap_residual": float(pure_residual),
            "q_only_mutation_difference": float(q_only_difference),
            "shared_masses_unchanged": True,
        }


def endpoint_gauge_crosscheck(parameters: MPParameters, dps: int) -> dict[str, Any]:
    with mpmath.workdps(dps):
        mp = mpmath.mpf
        p, x = parameters.p, parameters.x

        def gap(values: Sequence[mpmath.mpf]) -> mpmath.mpf:
            terms = evaluate_mpmath(values, parameters.beta, dps=dps)
            return terms.numerator - terms.ehx

        q0_a = (p, 0, mp(0), x, 0, 0, mp("0.2"), mp("0.6"), mp("0.9"))
        q0_b = (p, 0, mp(0), x, 0, 0, mp("0.8"), mp("0.4"), mp("0.1"))
        q1_a = (p, 0, mp(1), mp("0.2"), mp("0.6"), mp("0.9"), x, 0, 0)
        q1_b = (p, 0, mp(1), mp("0.8"), mp("0.4"), mp("0.1"), x, 0, 0)
        q0_residual = abs(gap(q0_a) - gap(q0_b))
        q1_residual = abs(gap(q1_a) - gap(q1_b))
        tolerance = mp(10) ** (-(dps - 20))
        if max(q0_residual, q1_residual) > tolerance:
            raise AssertionError("inactive endpoint component changed the gap")
        return {
            "dps": dps,
            "q0_inactive_residual": float(q0_residual),
            "q1_inactive_residual": float(q1_residual),
            "distance_inactive_weight_exactly_zero": True,
        }


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--smooth-cutoff", type=str,
                        default=str(DEFAULT_SMOOTH_CUTOFF))
    parser.add_argument("--q-star", type=str, default=str(DEFAULT_Q_STAR))
    parser.add_argument("--interior-q-min", type=str,
                        default=str(DEFAULT_INTERIOR_Q_MIN))
    parser.add_argument("--tube-radius", type=str,
                        default=str(DEFAULT_TUBE_RADIUS))
    parser.add_argument("--endpoint-pieces", type=int,
                        default=DEFAULT_ENDPOINT_PIECES)
    parser.add_argument("--q-cells", type=int, default=DEFAULT_Q_CELLS)
    parser.add_argument("--d-cells", type=int, default=DEFAULT_D_CELLS)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)
    if args.dps < 80:
        parser.error("--dps must be at least 80")

    smooth_cutoff = Fraction(args.smooth_cutoff)
    q_star = Fraction(args.q_star)
    interior_q_min = Fraction(args.interior_q_min)
    tube_radius = Fraction(args.tube_radius)
    mp_parameters = solve_equation_parameters(args.dps + 30)
    parameters = certify_equation_parameters(mp_parameters)
    domain = validate_domain(parameters, smooth_cutoff)
    certificate = certify(
        parameters,
        smooth_cutoff=smooth_cutoff,
        q_star=q_star,
        interior_q_min=interior_q_min,
        tube_radius=tube_radius,
        endpoint_pieces=args.endpoint_pieces,
        q_cells=args.q_cells,
        d_cells=args.d_cells,
    )
    experiment = run_experiment(mp_parameters, args.dps)
    autodiff = autodiff_crosscheck(mp_parameters, args.dps)
    symmetry = symmetry_crosscheck(mp_parameters, args.dps)
    gauges = endpoint_gauge_crosscheck(mp_parameters, args.dps)
    quadratic_limit = (
        certificate.raw_second_coefficient
        / (2 * (parameters.p**2 + parameters.p * parameters.x**2)))

    print("PROVED [exact reduced algebra]: both the pure-d quadratic gap and "
          "distance carry q(1-q); the variance factor cancels in their ratio.")
    print("COMPUTATIONAL EVIDENCE [%d dps, direct Liu _formula, d=%s]:"
          % (args.dps, EXPERIMENT_D))
    for row in experiment:
        print("COMPUTATIONAL EVIDENCE: q=%s, 1-q=%s, gap/dist^2=%s"
              % (row.q, row.one_minus_q, row.ratio))
    print("MACHINE VERIFIED [Arb point limit]: the common d->0 raw-gap ratio "
          "is in %s." % quadratic_limit)
    print("PROVED [exact cubic weight]: (1-q)(-q)^3+q(1-q)^3="
          "q(1-q)(1-2q), and the transcribed EHXY cross sum has the same "
          "factor after exchanging its two arguments.")
    print("MACHINE VERIFIED [Arb, %d exact dyadic cells]: on |d|<=%s, "
          "|G'''_d/[q(1-q)]|<=%s; all entropy arguments lie in [%s,%s], "
          "where |h'''|<=%s.  No sign test divided by a ball touching zero."
          % (certificate.cubic_sweep.boxes, smooth_cutoff,
             certificate.cubic_sweep.bound, domain["entropy_argument_lower"],
             domain["entropy_argument_upper"], certificate.hthird_ceiling))
    print("MACHINE VERIFIED [directed Arb Taylor bound]: on the inner core, "
          "G/[q(1-q)d^2]>=%s, dist^2/[q(1-q)d^2]<=%s, hence "
          "G>=%s*dist^2 uniformly through both exact endpoints."
          % (certificate.core_gap_factor_lower,
             certificate.core_distance_factor_upper, CORE_KAPPA))
    print("PROVED [endpoint identity]: dG(V(q,d))/dq at q=0 equals "
          "H*D(Q_d), with D exactly as defined and certified in "
          "liu9_endpoint.py.")
    print("MACHINE VERIFIED [Arb centered forms on %d dyadic annulus cells]: "
          "for %s<=|d|<=%s, H*D(Q_d)>=%s*delta(Q_d)^2 and D(Q_d)>=%s."
          % (certificate.endpoint_sweep.boxes, smooth_cutoff,
             ENDPOINT_RADIUS, ENDPOINT_RATIO,
             certificate.endpoint_normalized_d0))
    print("MACHINE VERIFIED [Arb q remainder on %d dyadic cells]: "
          "|d^2G/dq^2|<=%s for 0<=q<=%s, so every q_*<=%s gives the "
          "endpoint estimate; the selected q_*=%s yields raw kappa>%s."
          % (certificate.remainder_sweep.boxes,
             certificate.remainder_sweep.bound, Q_REMAINDER_CAP,
             certificate.q_star_ceiling, q_star,
             certificate.endpoint_kappa_lower))
    print("PROVED [exact endpoints]: at q=0 P1, and at q=1 P0, has zero "
          "weight in EHXY, EHPI, EHX, and dist; d is a joint gauge and the "
          "endpoint assertion follows by the certified one-sided limit.")
    print("PROVED [transcribed symmetry]: only simultaneous "
          "(q,P0,P1)->(1-q,P1,P0) is guaranteed; the atom masses are shared "
          "and unchanged, while q->1-q alone is not a symmetry.")
    print("MACHINE VERIFIED [symmetry mutation guard]: full-swap residual %.3e; "
          "q-only mutation changes the gap by %.3e."
          % (symmetry["full_component_swap_residual"],
             symmetry["q_only_mutation_difference"]))
    print("MACHINE VERIFIED [forward mode through Liu _formula]: %d checks; "
          "worst value, d-third, endpoint-first, and distance residuals are "
          "%.3e, %.3e, %.3e, %.3e."
          % (autodiff["checks"], autodiff["worst_value_relative"],
             autodiff["worst_d_third_relative"],
             autodiff["worst_endpoint_first_relative"],
             autodiff["worst_distance_relative"]))
    print("PROVED [exact seam inequality]: endpoint ranges reach q_*=%s, the "
          "interior begins at %s, and rho^2<=p^2 eps_sm^2 q_*(1-q_*) "
          "holds for rho=%s; therefore no pure-d tube point is lost between "
          "the q regimes."
          % (q_star, interior_q_min, tube_radius))
    print("CONDITIONAL: integrating this seam with a smooth-chart theorem "
          "requires that theorem to be valid at the supplied eps_sm=%s and "
          "interior threshold %s; this module does not inflate its radius."
          % (smooth_cutoff, interior_q_min))
    print("OPEN: inf D(Q)/delta(Q)^2 for arbitrary inactive three-atom Q, "
          "and the finite-q remainder outside the stated pure-d annulus, "
          "remain assigned to the complement/boundary pieces.")

    payload: dict[str, Any] = {
        "tool": "liu9_qdegenerate.py",
        "outcome": "PASS",
        "claim_status": "MACHINE VERIFIED",
        "parameters": {
            "x": str(parameters.x),
            "p": str(parameters.p),
            "mean": str(parameters.mean),
            "beta": str(parameters.beta),
            "root_lo": str(parameters.root_lo),
            "root_hi": str(parameters.root_hi),
        },
        "domain": domain,
        "theorem": {
            "endpoint_domain": (
                f"q in [0,{q_star}] union [1-{q_star},1], "
                f"{smooth_cutoff}<=|d|<={ENDPOINT_RADIUS}"
            ),
            "inner_core": f"q in [0,1], |d|<={smooth_cutoff}",
            "tube_radius": str(tube_radius),
            "endpoint_kappa": str(ENDPOINT_KAPPA),
            "core_kappa": str(CORE_KAPPA),
            "equality": "d=0, or q in {0,1} modulo inactive-component gauges",
        },
        "certificate": certificate.as_json(),
        "quadratic_limit": str(quadratic_limit),
        "experiment": [asdict(row) for row in experiment],
        "autodiff_crosscheck": autodiff,
        "symmetry_crosscheck": symmetry,
        "endpoint_gauge_crosscheck": gauges,
        "structural_findings": {
            "q_factors_cancel": True,
            "quadratic_factor": "q*(1-q)",
            "cubic_factor_at_d_zero": "q*(1-q)*(1-2*q)",
            "valid_symmetry": "(q,P0,P1)->(1-q,P1,P0)",
            "invalid_symmetry": "q->1-q with P0,P1 fixed",
        },
        "seam": {
            "status": "PROVED for the supplied pure-d tube hypotheses; CONDITIONAL on an interior smooth-chart certificate",
            "endpoint_q_max": str(q_star),
            "interior_q_min": str(interior_q_min),
            "regime_cover_complete": regime_cover_is_complete(
                q_star, interior_q_min),
            "required_inequality": "rho^2 <= p^2*eps_sm^2*q_star*(1-q_star)",
            "rhs_lower": str(certificate.seam_radius_squared_ceiling),
        },
        "remaining_open": (
            "unrestricted inactive three-atom D(Q)/delta(Q)^2 and finite-q "
            "remainder outside the pure-d endpoint chart"
        ),
    }
    payload["report_sha256"] = _canonical_digest(payload)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print("MACHINE VERIFIED: report written to %s with canonical SHA-256 %s."
              % (args.output, payload["report_sha256"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
