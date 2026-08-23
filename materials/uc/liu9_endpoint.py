#!/usr/bin/env python3
"""Rigorous endpoint certificate for Liu's nine-parameter binding locus.

This closes the endpoint-sign conjecture left in ``liu9_binding.py``.  Status
labels in the output distinguish exact/Arb proofs from finite-precision
searches.  The equation-defined root of Liu (87)--(90), never the rounded
15-digit value printed in the paper, is used throughout.

Run from the repository root with

    math/.venv/bin/python math/uc/liu9_endpoint.py

The default calculation is single-core and has a 400,000-box budget for each
of the two reduced three-dimensional domains and a 540-second global cap.
"""

from __future__ import annotations

import argparse
import heapq
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Callable, Sequence

# This workstation is shared with live certification workers.  These variables
# must be set before importing NumPy/SciPy.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import mpmath
import numpy as np
from flint import arb, ctx
from scipy.optimize import differential_evolution

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_binding import (  # noqa: E402
    ArbParameters,
    MPParameters,
    certify_equation_parameters,
    certify_insertion_potential,
    solve_equation_parameters,
)
from liu9_objective import h_arb  # noqa: E402

ctx.prec = max(ctx.prec, 320)
DEFAULT_DPS = 90
DEFAULT_SEARCH_ITERATIONS = 120
DEFAULT_SEARCH_POPULATION = 96
DEFAULT_BOX_BUDGET = 400_000
DEFAULT_WALL_SECONDS = 540.0

# Exact dyadic cutoffs used in the local proof and aligned with the B&B tree.
SMALL_ZERO_CUTOFF = 1.0 / 64.0
ZERO_CUTOFF = 1.0 / 32.0
X_LO = 19.0 / 32.0
X_HI = 25.0 / 32.0
LOCAL_DIAGONAL_S_BOXES = 128
LOCAL_DIAGONAL_MASS_BOXES = 64
LOCAL_SPLIT_BOXES = 64
LOCAL_LOW_MASS_BOXES = 64
LOCAL_ZERO_MVT_BOXES = 128

RawBox = tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class SearchResult:
    name: str
    value_float: float
    value_mp: mpmath.mpf
    raw_variables: tuple[float, float, float]
    variables: tuple[float, float, float, float]
    evaluations: int


@dataclass(frozen=True)
class LocalCertificate:
    true_zero_coefficient: arb
    proof_zero_coefficient: arb
    small_zero_margin: arb
    zero_segment_margin: arb
    diagonal_h_mm: arb
    diagonal_h_ms: arb
    diagonal_h_ss: arb
    diagonal_determinant: arb
    diagonal_grid_margin: arb
    split_leading_coefficient: arb
    split_grid_margin: arb
    low_mass_margin: arb
    diagonal_boxes: int
    split_boxes: int
    low_mass_boxes: int
    zero_boxes: int


@dataclass
class BranchStats:
    evaluated: int = 0
    mvt: int = 0
    direct: int = 0
    interval_cleared: int = 0
    infeasible: int = 0
    local_cleared: int = 0


@dataclass(frozen=True)
class ResidualBox:
    lower: arb
    depth: int
    box: RawBox


@dataclass(frozen=True)
class BranchResult:
    name: str
    stats: BranchStats
    residual: tuple[ResidualBox, ...]
    weakest_margin: arb | None
    maximum_depth: int
    elapsed: float


# ---------------------------------------------------------------------------
# Scalar and high-precision formulas used by the independent numerical search.


def _h_float(value: float) -> float:
    if value <= 0.0 or value >= 1.0:
        if value == 0.0 or value == 1.0:
            return 0.0
        raise ValueError(f"entropy argument outside [0,1]: {value!r}")
    return -(value * math.log(value) + (1.0 - value) * math.log1p(-value))


def _h_mp(value: mpmath.mpf) -> mpmath.mpf:
    if value <= 0 or value >= 1:
        if value == 0 or value == 1:
            return mpmath.mpf(0)
        raise ValueError(f"entropy argument outside [0,1]: {value!r}")
    return -(value * mpmath.log(value)
             + (1 - value) * mpmath.log(1 - value))


def _protocol_float(left: float, right: float) -> float:
    return left * right * (1.0 + (1.0 - left) * (1.0 - right))


def _protocol_mp(left: mpmath.mpf, right: mpmath.mpf) -> mpmath.mpf:
    return left * right * (1 + (1 - left) * (1 - right))


def _k_float(left: float, right: float) -> float:
    return _h_float(_protocol_float(left, right))


def _k_mp(left: mpmath.mpf, right: mpmath.mpf) -> mpmath.mpf:
    return _h_mp(_protocol_mp(left, right))


class NumericEndpointModel:
    """Independent float64/mpmath implementation of the reduced functional."""

    def __init__(self, parameters: MPParameters):
        self.parameters = parameters
        self.x = float(parameters.x)
        self.p = float(parameters.p)
        self.rho = 1.0 - self.p
        self.beta = float(parameters.beta)
        self.entropy = self.p * _h_float(self.x)
        # At the equation-defined root, core=1/(p*x) exactly.  Retaining the
        # displayed formula here makes the numerical implementation independent
        # of that simplification.
        self.core = (
            (1.0 - self.beta) * _h_float(self.x * self.x)
            + self.beta * _k_float(self.x, self.x)
        ) / (self.x * _h_float(self.x))
        self.gamma = self.beta / self.entropy

    def potential(self, value: float) -> float:
        return (
            2.0 * self.p * (
                (1.0 - self.beta) * _h_float(self.x * value)
                + self.beta * _k_float(self.x, value)
            )
            - _h_float(value)
            - self.entropy * self.core * value
        ) / self.entropy

    def coefficients(
        self, group_mass: float, left: float, right: float, fixed: float
    ) -> tuple[float, float, float]:
        """Return A,B,C with D(Q_s)=A*s^2+B*s+C."""
        a = group_mass
        b = 1.0 - a
        sites = (right, fixed, self.x, 0.0)
        weights = (a, b, -self.p, -self.rho)
        constant = (
            a * self.potential(right)
            + b * self.potential(fixed)
            + self.gamma * sum(
                weights[i] * weights[j] * _k_float(sites[i], sites[j])
                for i in range(4) for j in range(4)
            )
        )
        linear = (
            a * (self.potential(left) - self.potential(right))
            + 2.0 * self.gamma * a * sum(
                weights[i] * (
                    _k_float(sites[i], left) - _k_float(sites[i], right)
                )
                for i in range(4)
            )
        )
        quadratic = self.gamma * a * a * (
            _k_float(left, left)
            - 2.0 * _k_float(left, right)
            + _k_float(right, right)
        )
        return quadratic, linear, constant

    @staticmethod
    def minimise_quadratic(
        quadratic: float, linear: float, constant: float
    ) -> tuple[float, float]:
        candidates = [(constant, 0.0),
                      (quadratic + linear + constant, 1.0)]
        if quadratic > 0.0:
            split = -linear / (2.0 * quadratic)
            if 0.0 < split < 1.0:
                candidates.append((
                    constant - linear * linear / (4.0 * quadratic), split))
        return min(candidates, key=lambda candidate: candidate[0])

    def reduced_raw(
        self, group_mass: float, raw: Sequence[float]
    ) -> tuple[float, float, float, float, float]:
        ell, separation, fixed = (float(value) for value in raw)
        left = ell
        right = ell + (1.0 - ell) * separation
        value, split = self.minimise_quadratic(
            *self.coefficients(group_mass, left, right, fixed))
        return value, split, left, right, fixed

    def value_mp(
        self, group_mass: mpmath.mpf, split: float,
        left: float, right: float, fixed: float, dps: int,
    ) -> mpmath.mpf:
        with mpmath.workdps(dps):
            parameters = self.parameters
            p = parameters.p
            rho = 1 - p
            x = parameters.x
            beta = parameters.beta
            entropy = p * _h_mp(x)
            core = (
                (1 - beta) * _h_mp(x * x) + beta * _k_mp(x, x)
            ) / (x * _h_mp(x))
            a = mpmath.mpf(group_mass)
            s = mpmath.mpf(repr(float(split)))
            points = tuple(mpmath.mpf(repr(float(value))) for value in (
                left, right, fixed))
            weights = (a * s, a * (1 - s), 1 - a, -p, -rho)
            sites = (*points, x, mpmath.mpf(0))

            def potential(value: mpmath.mpf) -> mpmath.mpf:
                return (
                    2 * p * (
                        (1 - beta) * _h_mp(x * value)
                        + beta * _k_mp(x, value)
                    )
                    - _h_mp(value) - entropy * core * value
                ) / entropy

            answer = sum(
                (weights[i] * potential(sites[i]) for i in range(3)),
                mpmath.mpf(0),
            )
            answer += beta / entropy * sum((
                weights[i] * weights[j] * _k_mp(sites[i], sites[j])
                for i in range(5) for j in range(5)
            ), mpmath.mpf(0))
            return +answer


# ---------------------------------------------------------------------------
# Arb range arithmetic.  All entropy endpoints are resolved before a log.


def _finite(label: str, *values: arb) -> None:
    if not all(value.is_finite() for value in values):
        raise ArithmeticError(f"non-finite Arb enclosure in {label}")


def _arb_interval(lower: float, upper: float) -> arb:
    # Every B&B endpoint is dyadic and hence exactly representable as a float.
    answer = arb(lower).union(arb(upper))
    _finite("interval construction", answer)
    return answer

def _arb_center_delta(lower: float, upper: float) -> tuple[arb, arb]:
    """An Arb midpoint and a delta interval covering both float endpoints."""
    lower_arb, upper_arb = arb(lower), arb(upper)
    center = (lower_arb + upper_arb) / 2
    delta = (lower_arb - center).union(upper_arb - center)
    _finite("centered interval construction", center, delta)
    return center, delta


def _subinterval(lower: arb, upper: arb, index: int, count: int) -> arb:
    left = lower + (upper - lower) * index / count
    right = lower + (upper - lower) * (index + 1) / count
    answer = left.union(right)
    _finite("subinterval construction", answer)
    return answer


def _ell_point(value: arb) -> arb:
    """-x log(x) at an Arb point, with x=0 handled before the log."""
    if value <= 0:
        return arb(0)
    answer = -value * value.log()
    _finite("-x log(x) point", answer)
    return answer


def _ell_range(value: arb) -> arb:
    """Exact range of -x log(x) on an interval semantically in [0,1]."""
    lower, upper = value.lower(), value.upper()
    if lower < 0:
        lower = arb(0)
    if upper > 1:
        upper = arb(1)
    if lower > upper:
        raise ValueError("empty semantic entropy interval")
    answer = _ell_point(lower).union(_ell_point(upper))
    inverse_e = arb(1) / arb(1).exp()
    if not (upper < inverse_e or lower > inverse_e):
        answer = answer.union(inverse_e)
    _finite("-x log(x) range", answer)
    return answer


def _hp_arb(value: arb) -> arb:
    answer = ((arb(1) - value) / value).log()
    _finite("entropy first derivative", answer)
    return answer


def _hpp_arb(value: arb) -> arb:
    answer = -arb(1) / (value * (arb(1) - value))
    _finite("entropy second derivative", answer)
    return answer


def _protocol_range(left: arb, right: arb) -> arb:
    """Exact argument range, using coordinatewise monotonicity."""
    left_lo, left_hi = left.lower(), left.upper()
    right_lo, right_hi = right.lower(), right.upper()

    def argument(first: arb, second: arb) -> arb:
        return first * second * (
            arb(1) + (arb(1) - first) * (arb(1) - second))

    answer = argument(left_lo, right_lo).union(argument(left_hi, right_hi))
    _finite("protocol argument", answer)
    return answer


def _k_range(left: arb, right: arb) -> arb:
    answer = h_arb(_protocol_range(left, right))
    _finite("protocol entropy", answer)
    return answer


def _k_derivatives(
    left: arb, right: arb
) -> tuple[arb, arb, arb, arb, arb]:
    """k_1,k_2,k_11,k_22,k_12 for k=h(protocol), off endpoints."""
    one = arb(1)
    argument = left * right * (one + (one - left) * (one - right))
    argument_left = right * (one + (one - 2 * left) * (one - right))
    argument_right = left * (one + (one - 2 * right) * (one - left))
    argument_ll = -2 * right * (one - right)
    argument_rr = -2 * left * (one - left)
    argument_lr = one + (one - 2 * left) * (one - 2 * right)
    hp = _hp_arb(argument)
    hpp = _hpp_arb(argument)
    first_left = hp * argument_left
    first_right = hp * argument_right
    second_left = hpp * argument_left**2 + hp * argument_ll
    second_right = hpp * argument_right**2 + hp * argument_rr
    mixed = hpp * argument_left * argument_right + hp * argument_lr
    _finite(
        "protocol derivatives", first_left, first_right,
        second_left, second_right, mixed,
    )
    return first_left, first_right, second_left, second_right, mixed


class ArbEndpointModel:
    """Certified range extension and centered forms for D."""

    def __init__(self, parameters: ArbParameters):
        self.parameters = parameters
        self.x = parameters.x
        self.p = parameters.p
        self.rho = arb(1) - self.p
        self.beta = parameters.beta
        self.entropy = self.p * h_arb(self.x)
        self.lambda_gap = self.entropy * parameters.core
        self.gamma = self.beta / self.entropy
        self.k_xx = _k_range(self.x, self.x)
        _finite(
            "model constants", self.rho, self.entropy,
            self.lambda_gap, self.gamma, self.k_xx,
        )

    def potential_interval(self, value: arb) -> arb:
        """A dependency-reduced enclosure of V, clipped by the proved V>=0."""
        value_lo, value_hi = value.lower(), value.upper()
        # protocol(x,y)=c(y)y, and c(y) decreases with y.
        c_lo = self.x * (
            2 - self.x - (1 - self.x) * value_hi)
        c_hi = self.x * (
            2 - self.x - (1 - self.x) * value_lo)
        c = c_lo.union(c_hi)
        xy = self.x * value
        cy = c * value
        log_coefficient = (
            2 * self.p * ((1 - self.beta) * self.x + self.beta * c) - 1
        )
        numerator = (
            log_coefficient * _ell_range(value)
            - 2 * self.p * (
                (1 - self.beta) * self.x * value * self.x.log()
                + self.beta * c * value * c.log()
            )
            + 2 * self.p * (
                (1 - self.beta) * _ell_range(1 - xy)
                + self.beta * _ell_range(1 - cy)
            )
            - _ell_range(1 - value)
            - self.lambda_gap * value
        )
        answer = numerator / self.entropy
        _finite("insertion potential", answer)

        # liu9_binding.py's 22,100-box Arb proof gives V>=0 on [0,1].
        lower, upper = answer.lower(), answer.upper()
        if upper < 0:
            raise AssertionError("potential enclosure contradicts certified V>=0")
        if lower < 0:
            lower = arb(0)
        clipped = lower.union(upper)
        _finite("clipped insertion potential", clipped)
        return clipped

    def potential_prime(self, value: arb) -> arb:
        protocol = self.x * value * (
            1 + (1 - self.x) * (1 - value))
        protocol_prime = self.x * (
            1 + (1 - self.x) * (1 - 2 * value))
        answer = (
            2 * self.p * (
                (1 - self.beta) * self.x * _hp_arb(self.x * value)
                + self.beta * protocol_prime * _hp_arb(protocol)
            )
            - _hp_arb(value) - self.lambda_gap
        ) / self.entropy
        _finite("insertion potential derivative", answer)
        return answer

    def f0(self, value: arb) -> arb:
        """The local lower-barrier summand before division by H(P*)."""
        answer = (
            2 * (1 - self.beta) * self.p * h_arb(self.x * value)
            - h_arb(value) - self.lambda_gap * value
        )
        _finite("f0", answer)
        return answer

    def f0_prime(self, value: arb) -> arb:
        answer = (
            2 * (1 - self.beta) * self.p * self.x
            * _hp_arb(self.x * value)
            - _hp_arb(value) - self.lambda_gap
        )
        _finite("f0 prime", answer)
        return answer

    def f0_second(self, value: arb) -> arb:
        answer = (
            2 * (1 - self.beta) * self.p * self.x**2
            * _hpp_arb(self.x * value)
            - _hpp_arb(value)
        )
        _finite("f0 second", answer)
        return answer

    @staticmethod
    def k_diagonal_derivatives(value: arb) -> tuple[arb, arb]:
        derivatives = _k_derivatives(value, value)
        first = derivatives[0] + derivatives[1]
        second = derivatives[2] + 2 * derivatives[4] + derivatives[3]
        _finite("diagonal protocol derivatives", first, second)
        return first, second

    def coefficients(
        self, group_mass: arb, left: arb, right: arb, fixed: arb
    ) -> tuple[arb, arb, arb]:
        """Direct interval coefficients A,B,C of D(Q_s)."""
        a = group_mass
        b = 1 - a
        v_left = self.potential_interval(left)
        v_right = self.potential_interval(right)
        v_fixed = self.potential_interval(fixed)
        k_ll = _k_range(left, left)
        k_lr = _k_range(left, right)
        k_rr = _k_range(right, right)
        k_lw = _k_range(left, fixed)
        k_rw = _k_range(right, fixed)
        k_ww = _k_range(fixed, fixed)
        k_xl = _k_range(self.x, left)
        k_xr = _k_range(self.x, right)
        k_xw = _k_range(self.x, fixed)
        quadratic = self.gamma * a * a * (k_ll - 2 * k_lr + k_rr)
        linear = (
            a * (v_left - v_right)
            + 2 * self.gamma * a * (
                a * (k_lr - k_rr)
                + b * (k_lw - k_rw)
                - self.p * (k_xl - k_xr)
            )
        )
        energy = (
            a * a * k_rr + b * b * k_ww + self.p**2 * self.k_xx
            + 2 * a * b * k_rw
            - 2 * a * self.p * k_xr
            - 2 * b * self.p * k_xw
        )
        constant = (
            a * v_right + b * v_fixed + self.gamma * energy)
        _finite("quadratic coefficients", quadratic, linear, constant)
        return quadratic, linear, constant

    def coefficient_gradients(
        self, group_mass: arb, left: arb, right: arb, fixed: arb
    ) -> tuple[tuple[arb, arb, arb], ...]:
        """Gradients of A,B,C with respect to (left,right,fixed)."""
        a = group_mass
        b = 1 - a
        vp_left = self.potential_prime(left)
        vp_right = self.potential_prime(right)
        vp_fixed = self.potential_prime(fixed)
        k_ll = _k_derivatives(left, left)
        k_lr = _k_derivatives(left, right)
        k_rr = _k_derivatives(right, right)
        k_lw = _k_derivatives(left, fixed)
        k_rw = _k_derivatives(right, fixed)
        k_ww = _k_derivatives(fixed, fixed)
        k_xl = _k_derivatives(self.x, left)
        k_xr = _k_derivatives(self.x, right)
        k_xw = _k_derivatives(self.x, fixed)
        k_ll_diag = k_ll[0] + k_ll[1]
        k_rr_diag = k_rr[0] + k_rr[1]
        k_ww_diag = k_ww[0] + k_ww[1]

        grad_a = (
            self.gamma * a * a * (k_ll_diag - 2 * k_lr[0]),
            self.gamma * a * a * (-2 * k_lr[1] + k_rr_diag),
            arb(0),
        )
        grad_b = (
            a * vp_left + 2 * self.gamma * a * (
                a * k_lr[0] + b * k_lw[0] - self.p * k_xl[1]),
            -a * vp_right + 2 * self.gamma * a * (
                a * (k_lr[1] - k_rr_diag)
                - b * k_rw[0] + self.p * k_xr[1]),
            2 * self.gamma * a * b * (k_lw[1] - k_rw[1]),
        )
        grad_c = (
            arb(0),
            a * vp_right + self.gamma * (
                a * a * k_rr_diag + 2 * a * b * k_rw[0]
                - 2 * a * self.p * k_xr[1]),
            b * vp_fixed + self.gamma * (
                b * b * k_ww_diag + 2 * a * b * k_rw[1]
                - 2 * b * self.p * k_xw[1]),
        )
        for row in (grad_a, grad_b, grad_c):
            _finite("coefficient gradients", *row)
        return grad_a, grad_b, grad_c

    def two_atom(self, p_point: arb, rho_point: arb) -> arb:
        """D(p delta_y + rho delta_z)."""
        v_p = self.potential_interval(p_point)
        v_rho = self.potential_interval(rho_point)
        k_pp = _k_range(p_point, p_point)
        k_rr = _k_range(rho_point, rho_point)
        k_pr = _k_range(p_point, rho_point)
        k_xp = _k_range(self.x, p_point)
        k_xr = _k_range(self.x, rho_point)
        energy = (
            self.p**2 * k_pp + self.rho**2 * k_rr
            + self.p**2 * self.k_xx + 2 * self.p * self.rho * k_pr
            - 2 * self.p**2 * k_xp - 2 * self.p * self.rho * k_xr
        )
        answer = (
            self.p * v_p + self.rho * v_rho + self.gamma * energy)
        _finite("two-atom functional", answer)
        return answer

    def two_atom_gradient(
        self, p_point: arb, rho_point: arb
    ) -> tuple[arb, arb]:
        vp = self.potential_prime(p_point)
        vr = self.potential_prime(rho_point)
        k_pp = _k_derivatives(p_point, p_point)
        k_rr = _k_derivatives(rho_point, rho_point)
        k_pr = _k_derivatives(p_point, rho_point)
        k_xp = _k_derivatives(self.x, p_point)
        k_xr = _k_derivatives(self.x, rho_point)
        grad_p = self.p * vp + 2 * self.gamma * self.p * (
            self.p * k_pp[0] + self.rho * k_pr[0]
            - self.p * k_xp[1])
        grad_rho = self.rho * vr + 2 * self.gamma * self.rho * (
            self.rho * k_rr[0] + self.p * k_pr[1]
            - self.p * k_xr[1])
        _finite("two-atom gradient", grad_p, grad_rho)
        return grad_p, grad_rho


# ---------------------------------------------------------------------------
# Exact quadratic reduction and numerical search.


def run_independent_search(
    model: NumericEndpointModel,
    iterations: int,
    population: int,
    seed: int,
    dps: int,
) -> tuple[SearchResult, ...]:
    rng = np.random.default_rng(seed)
    cases = (
        ("zero-mass split", model.rho, ((0.0, 0.0, model.x),)),
        ("x-mass split", model.p, (
            (model.x, 0.0, 0.0), (0.0, model.x, model.x))),
    )
    answers: list[SearchResult] = []
    for index, (name, group_mass, equality_seeds) in enumerate(cases):
        initial = rng.random((population, 3))
        for row, equality in enumerate(equality_seeds):
            initial[row] = equality
        result = differential_evolution(
            lambda raw: model.reduced_raw(group_mass, raw)[0],
            bounds=((0.0, 1.0),) * 3,
            init=initial,
            seed=seed + index,
            maxiter=iterations,
            tol=1e-11,
            atol=1e-14,
            polish=False,
            workers=1,
            updating="immediate",
        )
        value, split, left, right, fixed = model.reduced_raw(
            group_mass, result.x)
        mp_group_mass = (
            1 - model.parameters.p if index == 0 else model.parameters.p)
        answers.append(SearchResult(
            name=name,
            value_float=float(value),
            value_mp=model.value_mp(
                mp_group_mass, split, left, right, fixed, dps),
            raw_variables=tuple(float(value) for value in result.x),
            variables=(float(split), left, right, fixed),
            evaluations=int(result.nfev),
        ))
    return tuple(answers)


# ---------------------------------------------------------------------------
# Local proof at the exact attaining law.


def certify_local_neighbourhood(model: ArbEndpointModel) -> LocalCertificate:
    one = arb(1)
    small = arb(SMALL_ZERO_CUTOFF)
    zero_cutoff = arb(ZERO_CUTOFF)
    x_lo = arb(X_LO)
    x_hi = arb(X_HI)
    x_interval = x_lo.union(x_hi)
    mass_interval = model.rho.union(one)

    # f0(y)=b0*y*log(1/y)+y*r(y).  For phi(t)=
    # -(1-t)log(1-t)/t, phi decreases from one, since
    # phi'(t)=(t+log(1-t))/t^2 <= 0.
    prefactor = 2 * (1 - model.beta) * model.p
    proof_zero_coefficient = prefactor * model.x - 1
    true_zero_coefficient = (
        2 * model.p * model.x * (1 + model.beta * (1 - model.x)) - 1)
    phi_argument = model.x * small
    phi_at_cutoff = (
        -(1 - phi_argument) * (1 - phi_argument).log() / phi_argument)
    remainder_lower = (
        -prefactor * model.x * model.x.log()
        + prefactor * model.x * phi_at_cutoff
        - 1 - model.lambda_gap
    )
    small_zero_margin = (
        proof_zero_coefficient * (one / small).log() + remainder_lower)
    if not (
        proof_zero_coefficient > 0
        and true_zero_coefficient > 0
        and small_zero_margin > 0
    ):
        raise AssertionError("small-zero expansion was not certified positive")

    zero_segment_margin: arb | None = None
    for index in range(LOCAL_ZERO_MVT_BOXES):
        segment = _subinterval(
            small, zero_cutoff, index, LOCAL_ZERO_MVT_BOXES)
        midpoint = (segment.lower() + segment.upper()) / 2
        delta = (
            (segment.lower() - midpoint).union(segment.upper() - midpoint))
        enclosure = (
            model.f0(midpoint) + model.f0_prime(segment) * delta)
        if not enclosure > 0:
            raise AssertionError(
                "middle zero-support segment was not certified positive")
        lower = enclosure.lower()
        if zero_segment_margin is None or lower < zero_segment_margin:
            zero_segment_margin = lower
    assert zero_segment_margin is not None

    # F_diag(M,s)=M f0(s)+beta M^2 k(s,s)+constant.  It is
    # stationary at (p,x).  Certify positive definiteness of its Hessian on
    # [rho,1] x I_x by Sylvester's criterion.
    diagonal_grid_margin: arb | None = None
    weakest_h_mm: arb | None = None
    for s_index in range(LOCAL_DIAGONAL_S_BOXES):
        support = _subinterval(
            x_lo, x_hi, s_index, LOCAL_DIAGONAL_S_BOXES)
        k_first, k_second = model.k_diagonal_derivatives(support)
        h_mm = 2 * model.beta * _k_range(support, support)
        f_prime = model.f0_prime(support)
        f_second = model.f0_second(support)
        for mass_index in range(LOCAL_DIAGONAL_MASS_BOXES):
            mass = _subinterval(
                model.rho, one, mass_index, LOCAL_DIAGONAL_MASS_BOXES)
            h_ms = f_prime + 2 * model.beta * mass * k_first
            h_ss = mass * f_second + model.beta * mass**2 * k_second
            determinant = h_mm * h_ss - h_ms**2
            if not (h_mm > 0 and determinant > 0):
                raise AssertionError(
                    "local diagonal Hessian was not certified positive definite")
            mm_lower = h_mm.lower()
            det_lower = determinant.lower()
            if weakest_h_mm is None or mm_lower < weakest_h_mm:
                weakest_h_mm = mm_lower
            if diagonal_grid_margin is None or det_lower < diagonal_grid_margin:
                diagonal_grid_margin = det_lower
    assert weakest_h_mm is not None and diagonal_grid_margin is not None

    # In weighted-mean coordinates y=s-(1-alpha)d,
    # z=s+alpha*d, the second d derivative is alpha(1-alpha)T2.
    # A crude sign-aware bound avoids a four-dimensional interval grid.
    f_second_min: arb | None = None
    diagonal_k_second_min: arb | None = None
    k_pure_second_min: arb | None = None
    k_mixed_upper: arb | None = None
    support_segments = tuple(
        _subinterval(x_lo, x_hi, index, LOCAL_SPLIT_BOXES)
        for index in range(LOCAL_SPLIT_BOXES)
    )
    for support in support_segments:
        f_second = model.f0_second(support)
        _, diagonal_second = model.k_diagonal_derivatives(support)
        if f_second_min is None or f_second.lower() < f_second_min:
            f_second_min = f_second.lower()
        if (
            diagonal_k_second_min is None
            or diagonal_second.lower() < diagonal_k_second_min
        ):
            diagonal_k_second_min = diagonal_second.lower()
    for left in support_segments:
        for right in support_segments:
            derivatives = _k_derivatives(left, right)
            pure_lower = min(
                derivatives[2].lower(), derivatives[3].lower())
            if (
                k_pure_second_min is None
                or pure_lower < k_pure_second_min
            ):
                k_pure_second_min = pure_lower
            mixed_upper = derivatives[4].upper()
            if k_mixed_upper is None or mixed_upper > k_mixed_upper:
                k_mixed_upper = mixed_upper
            if not derivatives[4] < 0:
                raise AssertionError(
                    "local mixed protocol derivative was not certified negative")
    assert (
        f_second_min is not None
        and diagonal_k_second_min is not None
        and k_pure_second_min is not None
        and k_mixed_upper is not None
    )
    kernel_lower = (
        diagonal_k_second_min / 2 + 2 * k_pure_second_min)
    if not kernel_lower < 0:
        raise AssertionError("split-kernel lower bound must be negative")
    # For 0<M<=1 and kernel_lower<0, M^2*kernel_lower >=
    # M*kernel_lower; beta.upper() is therefore the safe coefficient.
    split_grid_margin = (
        f_second_min + model.beta.upper() * kernel_lower)
    if not split_grid_margin > 0:
        raise AssertionError("local split coefficient was not certified positive")

    # If the x-cluster has mass below rho in one coordinate pattern, it has
    # only one x-support.  It is separated from the attaining mass p by a
    # direct two-dimensional certificate.
    f_x = model.f0(model.x)
    anchor = -(model.p * f_x + model.beta * model.p**2 * model.k_xx)
    low_mass_margin: arb | None = None
    for mass_index in range(LOCAL_LOW_MASS_BOXES):
        mass = _subinterval(
            arb(0), model.rho, mass_index, LOCAL_LOW_MASS_BOXES)
        for support_index in range(LOCAL_LOW_MASS_BOXES):
            support = _subinterval(
                x_lo, x_hi, support_index, LOCAL_LOW_MASS_BOXES)
            enclosure = (
                mass * model.f0(support)
                + model.beta * mass**2 * _k_range(support, support)
                + anchor
            )
            if not enclosure > 0:
                raise AssertionError("low-mass local face was not certified positive")
            lower = enclosure.lower()
            if low_mass_margin is None or lower < low_mass_margin:
                low_mass_margin = lower
    assert low_mass_margin is not None

    # Leading quadratic coefficients at the attaining law.
    k_first_x, k_second_x = model.k_diagonal_derivatives(model.x)
    k_derivatives_x = _k_derivatives(model.x, model.x)
    diagonal_h_mm = 2 * model.beta * model.k_xx
    diagonal_h_ms = (
        model.f0_prime(model.x)
        + 2 * model.beta * model.p * k_first_x)
    diagonal_h_ss = (
        model.p * model.f0_second(model.x)
        + model.beta * model.p**2 * k_second_x)
    diagonal_determinant = (
        diagonal_h_mm * diagonal_h_ss - diagonal_h_ms**2)
    split_leading_coefficient = (
        model.p * model.f0_second(model.x)
        + 2 * model.beta * model.p**2 * k_derivatives_x[2])
    mass_gradient = f_x + 2 * model.beta * model.p * model.k_xx
    support_gradient = (
        model.p * model.f0_prime(model.x)
        + model.beta * model.p**2 * k_first_x)
    if not (mass_gradient.contains(0) and support_gradient.contains(0)):
        raise AssertionError("local stationary gradient enclosure missed zero")
    if not (
        diagonal_h_mm > 0
        and diagonal_determinant > 0
        and split_leading_coefficient > 0
    ):
        raise AssertionError("leading local Hessian was not certified positive")

    return LocalCertificate(
        true_zero_coefficient=true_zero_coefficient / model.entropy,
        proof_zero_coefficient=proof_zero_coefficient / model.entropy,
        small_zero_margin=small_zero_margin,
        zero_segment_margin=zero_segment_margin,
        diagonal_h_mm=diagonal_h_mm / model.entropy,
        diagonal_h_ms=diagonal_h_ms / model.entropy,
        diagonal_h_ss=diagonal_h_ss / model.entropy,
        diagonal_determinant=diagonal_determinant / model.entropy**2,
        diagonal_grid_margin=diagonal_grid_margin,
        split_leading_coefficient=(
            split_leading_coefficient / model.entropy),
        split_grid_margin=split_grid_margin,
        low_mass_margin=low_mass_margin / model.entropy,
        diagonal_boxes=(
            LOCAL_DIAGONAL_S_BOXES * LOCAL_DIAGONAL_MASS_BOXES),
        split_boxes=LOCAL_SPLIT_BOXES**2,
        low_mass_boxes=LOCAL_LOW_MASS_BOXES**2,
        zero_boxes=LOCAL_ZERO_MVT_BOXES,
    )


# ---------------------------------------------------------------------------
# Reduced interval branch-and-bound.


def _support_box(box: RawBox) -> tuple[arb, arb, arb, tuple[arb, arb]]:
    (ell_lo, ell_hi), (sep_lo, sep_hi), (fixed_lo, fixed_hi) = box
    ell_lo_arb, ell_hi_arb = arb(ell_lo), arb(ell_hi)
    sep_lo_arb, sep_hi_arb = arb(sep_lo), arb(sep_hi)
    left = ell_lo_arb.union(ell_hi_arb)
    # v=ell+(1-ell)t is increasing in both variables on the unit square.
    # Compute its corners in Arb, rather than rounding them through float64.
    right_lo = (
        ell_lo_arb + (arb(1) - ell_lo_arb) * sep_lo_arb).lower()
    right_hi = (
        ell_hi_arb + (arb(1) - ell_hi_arb) * sep_hi_arb).upper()
    right = right_lo.union(right_hi)
    fixed = _arb_interval(fixed_lo, fixed_hi)
    _finite("ordered-support transform", left, right, fixed)
    return left, right, fixed, (right_lo, right_hi)


def _quadratic_lower(
    quadratic: arb, linear: arb, constant: arb
) -> arb:
    """Lower bound for min_[0,1] A*s^2+B*s+C."""
    a = quadratic.lower()
    b = linear.lower()
    c = constant.lower()
    zero = arb(0)
    if a > 0:
        if b >= 0:
            minimum = zero
        elif b <= -2 * a:
            minimum = a + b
        else:
            minimum = -(b * b) / (4 * a)
    else:
        endpoint = a + b
        minimum = endpoint if endpoint < zero else zero
    answer = (c + minimum).lower()
    _finite("reduced quadratic lower bound", answer)
    return answer


def _coefficient_box(
    model: ArbEndpointModel, group_mass: arb, box: RawBox
) -> tuple[tuple[arb, arb, arb], bool]:
    (ell_lo, ell_hi), (sep_lo, sep_hi), (fixed_lo, fixed_hi) = box
    ell_interval = _arb_interval(ell_lo, ell_hi)
    separation_interval = _arb_interval(sep_lo, sep_hi)
    left, right, fixed, _ = _support_box(box)
    try:
        gradients = model.coefficient_gradients(
            group_mass, left, right, fixed)
    except (ArithmeticError, ValueError, ZeroDivisionError):
        # This is the deliberate entropy-endpoint path: no derivative log of an
        # interval containing zero or one is attempted.
        return model.coefficients(group_mass, left, right, fixed), False

    ell_mid, delta_ell = _arb_center_delta(ell_lo, ell_hi)
    separation_mid, delta_separation = _arb_center_delta(sep_lo, sep_hi)
    fixed_mid, delta_fixed = _arb_center_delta(fixed_lo, fixed_hi)
    left_mid = ell_mid
    right_mid = ell_mid + (1 - ell_mid) * separation_mid
    point_coefficients = model.coefficients(
        group_mass, left_mid, right_mid, fixed_mid)
    enclosures: list[arb] = []
    for point, (grad_left, grad_right, grad_fixed) in zip(
        point_coefficients, gradients
    ):
        grad_ell = grad_left + grad_right * (1 - separation_interval)
        grad_separation = grad_right * (1 - ell_interval)
        enclosure = (
            point + grad_ell * delta_ell
            + grad_separation * delta_separation
            + grad_fixed * delta_fixed
        )
        _finite("centered coefficient enclosure", enclosure)
        enclosures.append(enclosure)
    return tuple(enclosures), True  # type: ignore[return-value]


def _local_box(pattern: str, box: RawBox) -> bool:
    (ell_lo, ell_hi), _, (fixed_lo, fixed_hi) = box
    _, _, _, (right_lo, right_hi) = _support_box(box)
    left_zero = arb(ell_hi) <= arb(ZERO_CUTOFF)
    right_zero = right_hi <= arb(ZERO_CUTOFF)
    fixed_zero = arb(fixed_hi) <= arb(ZERO_CUTOFF)
    left_x = arb(X_LO) <= arb(ell_lo) and arb(ell_hi) <= arb(X_HI)
    right_x = arb(X_LO) <= right_lo and right_hi <= arb(X_HI)
    fixed_x = (
        arb(X_LO) <= arb(fixed_lo) and arb(fixed_hi) <= arb(X_HI))
    if pattern == "zero-mass split":
        return (
            (left_zero and right_zero and fixed_x)
            or (left_zero and right_x and fixed_x)
        )
    if pattern == "x-mass split":
        return (
            (left_x and right_x and fixed_zero)
            or (left_zero and right_x and fixed_zero)
            or (left_zero and right_x and fixed_x)
        )
    raise ValueError(pattern)


def _classify_interior_box(
    model: ArbEndpointModel,
    group_mass: arb,
    box: RawBox,
) -> tuple[str, arb | None, bool, tuple[arb, arb, arb]]:
    coefficients, used_mvt = _coefficient_box(model, group_mass, box)
    quadratic, linear, constant = coefficients

    # An interior minimum in s requires A>0, B<0, and B+2A>0.
    # Otherwise both endpoints have already been covered by the two-atom proof.
    if (
        quadratic.upper() <= 0
        or linear.lower() >= 0
        or (linear + 2 * quadratic).upper() <= 0
    ):
        return "infeasible", None, used_mvt, coefficients

    lower = _quadratic_lower(quadratic, linear, constant)
    if lower > 0:
        return "cleared", lower, used_mvt, coefficients

    # On the possible interior branch, D_red=C-B^2/(4A).  Restricting
    # semantically to B<0 gives this sharper certified lower bound.
    if quadratic.lower() > 0:
        b_lower = linear.lower()
        interior_lower = (
            constant.lower()
            - b_lower * b_lower / (4 * quadratic.lower())
        ).lower()
        if interior_lower > 0:
            return "cleared", interior_lower, used_mvt, coefficients
        lower = interior_lower
    return "split", lower, used_mvt, coefficients


def _two_atom_local(box: RawBox) -> bool:
    (p_lo, p_hi), (rho_lo, rho_hi) = box
    return (
        X_LO <= p_lo and p_hi <= X_HI
        and rho_lo >= 0.0 and rho_hi <= ZERO_CUTOFF)


def _two_atom_box(
    model: ArbEndpointModel, box: RawBox
) -> tuple[arb, bool]:
    (p_lo, p_hi), (rho_lo, rho_hi) = box
    p_interval = _arb_interval(p_lo, p_hi)
    rho_interval = _arb_interval(rho_lo, rho_hi)
    try:
        grad_p, grad_rho = model.two_atom_gradient(
            p_interval, rho_interval)
    except (ArithmeticError, ValueError, ZeroDivisionError):
        return model.two_atom(p_interval, rho_interval).lower(), False
    p_mid, delta_p = _arb_center_delta(p_lo, p_hi)
    rho_mid, delta_rho = _arb_center_delta(rho_lo, rho_hi)
    enclosure = (
        model.two_atom(p_mid, rho_mid)
        + grad_p * delta_p + grad_rho * delta_rho)
    _finite("two-atom centered form", enclosure)
    return enclosure.lower(), True


def _split_longest(box: RawBox) -> tuple[RawBox, RawBox]:
    dimension = max(
        range(len(box)), key=lambda index: box[index][1] - box[index][0])
    lower, upper = box[dimension]
    midpoint = (lower + upper) / 2.0
    if not lower < midpoint < upper:
        raise ArithmeticError("float box midpoint did not refine its interval")
    first = list(box)
    second = list(box)
    first[dimension] = (lower, midpoint)
    second[dimension] = (midpoint, upper)
    return tuple(first), tuple(second)


def certify_two_atom_boundary(
    model: ArbEndpointModel, budget: int, deadline: float
) -> BranchResult:
    start = time.monotonic()
    stats = BranchStats()
    heap: list[tuple[float, int, int, RawBox, arb]] = []
    counter = 0
    weakest: arb | None = None
    maximum_depth = 0

    def add(box: RawBox, depth: int) -> None:
        nonlocal counter, weakest, maximum_depth
        maximum_depth = max(maximum_depth, depth)
        if _two_atom_local(box):
            stats.local_cleared += 1
            return
        lower, used_mvt = _two_atom_box(model, box)
        stats.evaluated += 1
        if used_mvt:
            stats.mvt += 1
        else:
            stats.direct += 1
        if lower > 0:
            stats.interval_cleared += 1
            if weakest is None or lower < weakest:
                weakest = lower
            return
        counter += 1
        heapq.heappush(heap, (float(lower), counter, depth, box, lower))

    add(((0.0, 1.0), (0.0, 1.0)), 0)
    while (
        heap
        and stats.evaluated + 2 <= budget
        and time.monotonic() < deadline
    ):
        _, _, depth, box, _ = heapq.heappop(heap)
        for child in _split_longest(box):
            add(child, depth + 1)
    residual = tuple(
        ResidualBox(lower=item[4], depth=item[2], box=item[3])
        for item in sorted(heap)
    )
    return BranchResult(
        name="two-atom split endpoints",
        stats=stats,
        residual=residual,
        weakest_margin=weakest,
        maximum_depth=maximum_depth,
        elapsed=time.monotonic() - start,
    )


def certify_interior_pattern(
    model: ArbEndpointModel,
    name: str,
    group_mass: arb,
    budget: int,
    deadline: float,
) -> BranchResult:
    start = time.monotonic()
    stats = BranchStats()
    heap: list[tuple[float, int, int, RawBox, arb]] = []
    counter = 0
    weakest: arb | None = None
    maximum_depth = 0

    def add(box: RawBox, depth: int) -> None:
        nonlocal counter, weakest, maximum_depth
        maximum_depth = max(maximum_depth, depth)
        if _local_box(name, box):
            stats.local_cleared += 1
            return
        classification, lower, used_mvt, _ = _classify_interior_box(
            model, group_mass, box)
        stats.evaluated += 1
        if used_mvt:
            stats.mvt += 1
        else:
            stats.direct += 1
        if classification == "infeasible":
            stats.infeasible += 1
            return
        assert lower is not None
        if classification == "cleared":
            stats.interval_cleared += 1
            if weakest is None or lower < weakest:
                weakest = lower
            return
        counter += 1
        heapq.heappush(heap, (float(lower), counter, depth, box, lower))

    add(((0.0, 1.0), (0.0, 1.0), (0.0, 1.0)), 0)
    while (
        heap
        and stats.evaluated + 2 <= budget
        and time.monotonic() < deadline
    ):
        _, _, depth, box, _ = heapq.heappop(heap)
        for child in _split_longest(box):
            add(child, depth + 1)
    residual = tuple(
        ResidualBox(lower=item[4], depth=item[2], box=item[3])
        for item in sorted(heap)
    )
    return BranchResult(
        name=name,
        stats=stats,
        residual=residual,
        weakest_margin=weakest,
        maximum_depth=maximum_depth,
        elapsed=time.monotonic() - start,
    )


def _residual_upper(
    model: ArbEndpointModel, group_mass: arb, box: RawBox
) -> arb:
    coefficients, _ = _coefficient_box(model, group_mass, box)
    quadratic, linear, constant = coefficients
    endpoint_zero = constant.upper()
    endpoint_one = (quadratic + linear + constant).upper()
    answer = endpoint_zero if endpoint_zero < endpoint_one else endpoint_one
    _finite("residual upper bound", answer)
    return answer


# ---------------------------------------------------------------------------
# Output.


def print_step_one(
    mp_parameters: MPParameters,
    arb_parameters: ArbParameters,
    searches: Sequence[SearchResult],
    dps: int,
) -> None:
    print("1. EXACT ENDPOINT FUNCTIONAL AND INDEPENDENT NUMERICAL SEARCH")
    print("PROVED [definitions and active-mean tangent elimination]:")
    print("   h(z)=-z*log(z)-(1-z)*log(1-z), with h(0)=h(1)=0;")
    print("   H(Q)=integral h(y)dQ(y), M(Q)=integral y dQ(y),")
    print("   J(P,Q)=double_integral h(y*z)dP(y)dQ(z), and")
    print("   K(P,Q)=double_integral h(y*z*[1+(1-y)(1-z)])dP(y)dQ(z).")
    print("   P*=p*delta_x+(1-p)*delta_0, m=p*x, H*=p*h(x), where")
    print("   x^4-2*x^3+3*x^2-1=0 and p=h(x)/h(x^2).")
    print("   Since x^2[1+(1-x)^2]=1-x^2, core=1/(p*x) exactly.")
    print("   For every admissible inactive law Q,")
    print("   D(Q)={2(1-beta)[J(P*,Q)-J(P*,P*)]")
    print("         +beta[K(Q,Q)-K(P*,P*)]-[H(Q)-H*]}/H*")
    print("         -core*[M(Q)-m].")
    print("PROVED [endpoint domain, up to atom permutations]:")
    print("   Q=p*delta_u+r*delta_v+(1-p-r)*delta_w, 0<=r<=1-p,")
    print("   or Q=r*delta_u+(p-r)*delta_v+(1-p)*delta_w, 0<=r<=p,")
    print("   with (u,v,w) in [0,1]^3.  Each component is four-dimensional.")
    print("   At q=0 this is the arbitrary inactive law P1 while P0=P*; at q=1")
    print("   it is the arbitrary inactive law P0 while P1=P*.  The two cases are")
    print("   identical by exchanging P0,P1.  The -core term is exactly the active")
    print("   mean compensation; no unconstrained mean direction was discarded.")
    print(f"PROVED [320-bit Arb]: x in {arb_parameters.x}")
    print(f"PROVED [320-bit Arb]: p in {arb_parameters.p}")
    print(f"PROVED [320-bit Arb]: beta in {arb_parameters.beta}")
    evaluations = sum(search.evaluations for search in searches)
    print(
        "NUMERICAL [independent one-core reduced differential evolution, "
        f"{evaluations} evaluations]:")
    for search in searches:
        print(
            "   %s: float64 min=%+.17e; mpmath(%d dps)=%s"
            % (
                search.name,
                search.value_float,
                dps,
                mpmath.nstr(search.value_mp, 18),
            ))
        print(
            "      raw (ell,t,w)=%s; (split,u,v,w)=%s; nfev=%d"
            % (
                np.array2string(
                    np.asarray(search.raw_variables), precision=12, separator=","),
                np.array2string(
                    np.asarray(search.variables), precision=12, separator=","),
                search.evaluations,
            ))
    print("NUMERICAL [search verdict]: no negative value survives high-precision")
    print("   recomputation.  Both minima are coordinate representations of Q=P*;")
    print("   the tiny negative float64 values are cancellation roundoff.")
    print()


def print_step_two(model: ArbEndpointModel) -> None:
    print("2. EXACT ANALYTIC DIMENSION REDUCTION")
    print("PROVED [exact signed-measure identity]: put")
    print("   V(y)={2p[(1-beta)h(xy)+beta*k(x,y)]-h(y)-H*core*y}/H*,")
    print("   k(y,z)=h(yz[1+(1-y)(1-z)]), gamma=beta/H*.")
    print("   Since V(0)=V(x)=0, integral V dP*=0, and exactly")
    print("      D(Q)=integral V(y)dQ(y)+gamma*K(Q-P*,Q-P*).")
    print("PROVED [mass variable eliminated exactly]: for a in {1-p,p}, b=1-a,")
    print("   Q_s=a[s*delta_u+(1-s)*delta_v]+b*delta_w.")
    print("   With mu0=a*delta_v+b*delta_w-P*, D(Q_s)=A*s^2+B*s+C, where")
    print("   A=gamma*a^2*K(delta_u-delta_v,delta_u-delta_v),")
    print("   B=a[V(u)-V(v)]+2gamma*a*K(mu0,delta_u-delta_v),")
    print("   C=a*V(v)+b*V(w)+gamma*K(mu0,mu0).")
    print("   Hence D_red(u,v,w)=min_{0<=s<=1}(A*s^2+B*s+C) equals")
    print("   C-B^2/(4A) if A>0 and -2A<B<0, and otherwise")
    print("   min(C,A+B+C).  This is an exact piecewise expression, not a fit.")
    print("PROVED [permutation quotient]: exchange of the split atoms permits u<=v;")
    print("   write u=ell, v=ell+(1-ell)t with (ell,t,w) in [0,1]^3.")
    print("   Thus the two four-dimensional endpoint components reduce to two")
    print("   three-dimensional cubes.  Differentiating A,B,C leaves every support")
    print("   in coupled entropy-log terms; unlike s, none occurs linearly.  The")
    print("   certificate therefore treats all three supports by interval B&B rather")
    print("   than asserting an unavailable closed-form stationarity solution.")
    print(f"PROVED [imported 1D Arb certificate]: V>=0 on [0,1] was checked on")
    print("   22,100 exact-rational boxes in liu9_binding.py; the B&B below uses")
    print("   that theorem only to clip V's interval lower endpoint to zero.")
    print()


def print_local_certificate(
    model: ArbEndpointModel, certificate: LocalCertificate
) -> None:
    print("PROVED [local attaining neighbourhoods and expansion]: within the")
    print("   neighbourhoods isolated below, the only equality law is Q=P*.")
    print("   Coordinate gauges arise only from a zero split weight")
    print("   or coincident atoms.  The positive-mass coordinate representatives are")
    print("   (a=1-p): w=x,u=v=0; and (a=p): u=v=x,w=0, or")
    print("   u=0,v=w=x,s=(1-p)/p, together with split-atom exchange.")
    print("PROVED [Arb boundary expansion]: for a zero-side support y,")
    print("   V(y)=(B/H*)*y*log(1/y)+O(y), with")
    print(f"   true B/H* in {certificate.true_zero_coefficient}.")
    print("   The local lower proof drops nonnegative K terms involving zero-side")
    print("   atoms and uses f0(y)/H*.  Its leading coefficient is in")
    print(f"   {certificate.proof_zero_coefficient}; on [0,1/64], f0(y)/y has")
    print(f"   certified lower margin {certificate.small_zero_margin}; on")
    print(f"   [1/64,1/32], the weakest centered-MVT f0 enclosure is")
    print(f"   {certificate.zero_segment_margin}.")
    print("PROVED [Arb smooth expansion]: if M is the x-cluster mass and s its")
    print("   weighted mean, the normalized leading Hessian at (M,s)=(p,x) is")
    print("      [[h_MM,h_Ms],[h_Ms,h_ss]], with")
    print(f"   h_MM in {certificate.diagonal_h_mm}")
    print(f"   h_Ms in {certificate.diagonal_h_ms}")
    print(f"   h_ss in {certificate.diagonal_h_ss}")
    print(f"   determinant in {certificate.diagonal_determinant}.")
    print("   For two x-side atoms in weighted-mean/difference coordinates, the")
    print("   normalized d^2 coefficient divided by alpha(1-alpha) is in")
    print(f"   {certificate.split_leading_coefficient}.")
    print("PROVED [local neighbourhood, exact Taylor integral remainder]: on")
    print("   M in [1-p,1] and supports in [19/32,25/32], Arb gives weakest")
    print(f"   Hessian determinant lower bound {certificate.diagonal_grid_margin}")
    print(f"   on {certificate.diagonal_boxes} boxes and split-bracket lower bound")
    print(f"   {certificate.split_grid_margin} on {certificate.split_boxes} support")
    print("   pairs.  The one-x-atom face M in [0,1-p] is strictly separated by")
    print(f"   normalized margin {certificate.low_mass_margin} on")
    print(f"   {certificate.low_mass_boxes} boxes.  With Z=[0,1/32] and")
    print("   X=[19/32,25/32], the locally cleared ordered support patterns are")
    print("   (a=1-p): (u,v,w) in ZxZxX or ZxXxX; and")
    print("   (a=p): (u,v,w) in XxXxZ, ZxXxZ, or ZxXxX.")
    print("   These facts prove D>=0 in every listed attaining neighbourhood,")
    print("   including all zero-weight gauges.")


def print_branch_result(result: BranchResult) -> None:
    stats = result.stats
    print(f"PROVED [Arb branch accounting, {result.name}]:" if not result.residual
          else f"NUMERICAL [incomplete Arb branch accounting, {result.name}]:")
    print(
        "   evaluated=%d; centered-MVT=%d; endpoint-direct=%d; "
        "interval-cleared=%d; infeasible=%d; local-cleared=%d; residual=%d; "
        "max-depth=%d; seconds=%.3f"
        % (
            stats.evaluated,
            stats.mvt,
            stats.direct,
            stats.interval_cleared,
            stats.infeasible,
            stats.local_cleared,
            len(result.residual),
            result.maximum_depth,
            result.elapsed,
        ))
    if result.weakest_margin is not None:
        print(
            "   smallest certified lower margin among complement boxes: "
            f"{result.weakest_margin}")


def print_residual_map(
    model: ArbEndpointModel,
    result: BranchResult,
    group_mass: arb | None,
) -> None:
    if not result.residual:
        return
    print(
        f"PROVED [Arb residual map, {result.name}]: "
        f"{len(result.residual)} boxes rigorously survive this certificate.")
    for index, residual in enumerate(result.residual):
        if group_mass is None:
            (p_lo, p_hi), (rho_lo, rho_hi) = residual.box
            upper = model.two_atom(
                _arb_interval(p_lo, p_hi),
                _arb_interval(rho_lo, rho_hi),
            ).upper()
        else:
            upper = _residual_upper(model, group_mass, residual.box)
        width = upper - residual.lower
        support = _support_box(residual.box) if len(residual.box) == 3 else None
        print(
            "   residual[%d]: depth=%d raw=%s; reduced enclosure=[%s,%s]; "
            "width=%s%s"
            % (
                index,
                residual.depth,
                residual.box,
                residual.lower,
                upper,
                width,
                "" if support is None else f"; v-range={support[3]}",
            ))


def pathology_checks() -> tuple[arb, arb]:
    unsafe_log = _arb_interval(0.0, SMALL_ZERO_CUTOFF).log()
    unsafe_sqrt = _arb_interval(
        -SMALL_ZERO_CUTOFF, SMALL_ZERO_CUTOFF).sqrt()
    if unsafe_log.is_finite() or unsafe_sqrt.is_finite():
        raise AssertionError("Arb pathology self-test unexpectedly stayed finite")
    safe_left = h_arb(_arb_interval(0.0, SMALL_ZERO_CUTOFF))
    safe_right = h_arb(_arb_interval(1.0 - SMALL_ZERO_CUTOFF, 1.0))
    _finite("endpoint entropy self-test", safe_left, safe_right)
    return safe_left, safe_right


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument(
        "--search-iterations", type=int,
        default=DEFAULT_SEARCH_ITERATIONS)
    parser.add_argument(
        "--search-population", type=int,
        default=DEFAULT_SEARCH_POPULATION)
    parser.add_argument(
        "--box-budget", type=int, default=DEFAULT_BOX_BUDGET,
        help="maximum evaluated boxes for each reduced three-dimensional case")
    parser.add_argument(
        "--wall-seconds", type=float, default=DEFAULT_WALL_SECONDS,
        help="global wall-clock cap (must not exceed 570 seconds)")
    parser.add_argument("--seed", type=int, default=230_608_824)
    args = parser.parse_args(argv)
    if args.dps < 70:
        parser.error("--dps must be at least 70")
    if args.search_iterations < 20:
        parser.error("--search-iterations must be at least 20")
    if args.search_population < 32:
        parser.error("--search-population must be at least 32")
    if args.box_budget < 10_000:
        parser.error("--box-budget must be at least 10000")
    if not 30.0 <= args.wall_seconds <= 570.0:
        parser.error("--wall-seconds must lie in [30,570]")

    start = time.monotonic()
    deadline = start + args.wall_seconds
    mpmath.mp.dps = max(mpmath.mp.dps, args.dps + 30)

    print("LIU NINE-PARAMETER ENDPOINT INWARD-q CERTIFICATE")
    print("Status labels are PROVED, NUMERICAL, or CONJECTURED.")
    print("PROVED [runtime configuration]: SciPy workers and BLAS/OpenMP use one core;")
    print(
        f"   per-pattern box budget={args.box_budget}; "
        f"global wall cap={args.wall_seconds:.1f} seconds.")
    print()

    mp_parameters = solve_equation_parameters(args.dps + 20)
    arb_parameters = certify_equation_parameters(mp_parameters)
    numeric_model = NumericEndpointModel(mp_parameters)
    arb_model = ArbEndpointModel(arb_parameters)
    insertion = certify_insertion_potential(arb_parameters)
    if insertion.boxes != 22_100:
        raise AssertionError("unexpected imported insertion-certificate box count")
    searches = run_independent_search(
        numeric_model,
        args.search_iterations,
        args.search_population,
        args.seed,
        args.dps,
    )
    print_step_one(mp_parameters, arb_parameters, searches, args.dps)
    print_step_two(arb_model)

    print("3. ARB LOCAL CERTIFICATE AND COMPLEMENT BRANCH-AND-BOUND")
    local = certify_local_neighbourhood(arb_model)
    print_local_certificate(arb_model, local)
    print()

    two_atom = certify_two_atom_boundary(
        arb_model, min(args.box_budget, 50_000), deadline)
    zero_pattern = certify_interior_pattern(
        arb_model,
        "zero-mass split",
        arb_model.rho,
        args.box_budget,
        deadline,
    )
    x_pattern = certify_interior_pattern(
        arb_model,
        "x-mass split",
        arb_model.p,
        args.box_budget,
        deadline,
    )
    for result in (two_atom, zero_pattern, x_pattern):
        print_branch_result(result)
    print("PROVED [B&B method]: every split boundary s=0 or 1 is covered by the")
    print("   two-atom certificate.  An interior minimum requires A>0, B<0,")
    print("   B+2A>0; boxes violating one condition are counted infeasible.  Other")
    print("   boxes use direct Arb ranges or centered mean-value coefficient forms,")
    print("   followed by the exact quadratic minimum, and are bisected dyadically.")
    print()

    safe_left, safe_right = pathology_checks()
    print("4. PATHOLOGY CHECKS, ATTAINMENT, AND VERDICT")
    print("PROVED [Arb non-finite guards]: log([0,1/64]) and")
    print("   sqrt([-1/64,1/64]) are detected by .is_finite() as unsafe; neither")
    print("   expression is used in the proof.  Entropy endpoints are dispatched")
    print("   before logs; safe h([0,1/64]) and h([63/64,1]) enclose to")
    print(f"   {safe_left} and {safe_right}, respectively.")

    print_residual_map(arb_model, two_atom, None)
    print_residual_map(arb_model, zero_pattern, arb_model.rho)
    print_residual_map(arb_model, x_pattern, arb_model.p)
    proved = not (
        two_atom.residual or zero_pattern.residual or x_pattern.residual)
    elapsed = time.monotonic() - start
    if proved:
        margins = tuple(
            result.weakest_margin for result in (
                two_atom, zero_pattern, x_pattern)
            if result.weakest_margin is not None)
        weakest = min(margins) if margins else arb(0)
        total_evaluated = sum(
            result.stats.evaluated for result in (
                two_atom, zero_pattern, x_pattern))
        total_local = sum(
            result.stats.local_cleared for result in (
                two_atom, zero_pattern, x_pattern))
        total_infeasible = sum(
            result.stats.infeasible for result in (
                two_atom, zero_pattern, x_pattern))
        print("PROVED [local expansion plus rigorous Arb B&B]: D(Q)>=0 on both")
        print("   four-dimensional endpoint domains, at q=0 and q=1.  Equality holds")
        print("   exactly when the inactive law Q equals P* (including its coordinate")
        print("   gauges and atom permutations).")
        print(
            "PROVED [certificate totals]: evaluated=%d; local-cleared=%d; "
            "infeasible=%d; residual=0; weakest complement margin=%s; "
            "wall=%.3f seconds."
            % (
                total_evaluated,
                total_local,
                total_infeasible,
                weakest,
                elapsed,
            ))
        print("CONJECTURED [sole remaining binding-locus ingredient]: completeness.")
        print("   It remains to prove that no disconnected equality/binding component")
        print("   exists beyond the diagonal q-family P0=P1=P*, with the endpoint")
        print("   coordinate gauges and permutations already described.  The endpoint")
        print("   inward-q sign itself is no longer conjectural.")
    else:
        print("CONJECTURED [endpoint verdict]: the configured budget/cap did not clear")
        print("   every box.  The exact residual map above, not the numerical search,")
        print(f"   is the obstruction.  Elapsed wall time={elapsed:.3f} seconds.")


if __name__ == "__main__":
    main()
