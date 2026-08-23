#!/usr/bin/env python3
"""Certified local analysis of Liu's nine-parameter binding locus.

The mathematical value called c' in Liu's paper is defined by equations
(87)--(90); the printed decimals in (91)--(94) are approximations.  This
script keeps that distinction explicit.  It

* parameterises every algebraically forced equality stratum coming from
  Liu's two-atom optimizer and the q-invariant diagonal P0=P1;
* classifies every feasible local direction at its interior diagonal stratum;
* reduces every inward direction at the q=0,1 endpoint strata to an explicit
  four-variable first-variation functional and searches that functional on
  one core; and
* recomputes the paired-atom LOSS/GAIN ratio at the binding distribution.

Only statements labelled PROVED are algebraic or Arb-certified.  NUMERICAL
statements are finite-precision computations.  CONJECTURED statements are
interpretations that would require a global proof of Liu's Hypothesis 2.

Run from the repository root with

    math/.venv/bin/python math/uc/liu9_binding.py
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from fractions import Fraction
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

from liu9_objective import (  # noqa: E402
    LIU_BETA_DECIMAL,
    LIU_C_DECIMAL,
    evaluate_arb,
    evaluate_mpmath,
    h_arb,
    mean_value,
)

ctx.prec = max(ctx.prec, 320)
DEFAULT_DPS = 90
DEFAULT_SEARCH_ITERATIONS = 140
DEFAULT_SEARCH_POPULATION = 128
ROOT_DECIMAL_PLACES = 70


@dataclass(frozen=True)
class MPParameters:
    x: mpmath.mpf
    p: mpmath.mpf
    mean: mpmath.mpf
    c: mpmath.mpf
    beta: mpmath.mpf


@dataclass(frozen=True)
class ArbParameters:
    x: arb
    p: arb
    mean: arb
    c: arb
    beta: arb
    core: arb
    core_prime: arb
    core_second: arb
    root_lo: Fraction
    root_hi: Fraction


@dataclass(frozen=True)
class CurvatureCertificate:
    face: arb
    split_unit: arb
    boundary_barrier: arb
    gain: arb
    loss: arb
    loss_over_gain: arb
    gain_minus_loss: arb


@dataclass(frozen=True)
class InsertionCertificate:
    small_y_margin: arb
    left_mvt_lower: arb
    root_convexity: arb
    right_mvt_lower: arb
    endpoint_interval: arb
    boxes: int

@dataclass(frozen=True)
class EndpointSearch:
    name: str
    value_float: float
    value_mp: mpmath.mpf
    variables: tuple[float, ...]
    evaluations: int


def _h_mp(x: mpmath.mpf) -> mpmath.mpf:
    if x == 0 or x == 1:
        return mpmath.mpf(0)
    return -x * mpmath.log(x) - (1 - x) * mpmath.log(1 - x)


def _hp_mp(x: mpmath.mpf) -> mpmath.mpf:
    return mpmath.log((1 - x) / x)


def _hpp_mp(x: mpmath.mpf) -> mpmath.mpf:
    return -1 / (x * (1 - x))


def _h_arb(x: arb) -> arb:
    """Binary entropy for an Arb interval strictly inside (0,1)."""
    one = arb(1)
    return -(x * x.log() + (one - x) * (one - x).log())


def _h_arb_endpoint_safe(x: arb) -> arb:
    if x == 0 or x == 1:
        return arb(0)
    return _h_arb(x)


def _hp_arb(x: arb) -> arb:
    return ((arb(1) - x) / x).log()


def _hpp_arb(x: arb) -> arb:
    return -arb(1) / (x * (arb(1) - x))


def _arb_fraction(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _mp_fraction(value: Fraction) -> mpmath.mpf:
    return mpmath.mpf(value.numerator) / value.denominator


def _protocol_argument(x, y):
    return x * y * (1 + (1 - x) * (1 - y))


def _core_mp(x: mpmath.mpf, beta: mpmath.mpf) -> mpmath.mpf:
    numerator = ((1 - beta) * _h_mp(x * x)
                 + beta * _h_mp(_protocol_argument(x, x)))
    return numerator / (x * _h_mp(x))


def solve_equation_parameters(dps: int) -> MPParameters:
    """Solve Liu (87)--(90), not the rounded decimals (91)--(94)."""
    with mpmath.workdps(dps):
        one = mpmath.mpf(1)
        x = mpmath.findroot(
            lambda z: z**4 - 2 * z**3 + 3 * z**2 - 1,
            mpmath.mpf("0.690787593924988"),
        )
        p = _h_mp(x) / _h_mp(x * x)
        mean = p * x
        c = one - mean
        u = x * x
        t = _protocol_argument(x, x)
        du = 2 * x
        dt = 2 * x + 2 * x * (one - x) ** 2 - 2 * x * x * (one - x)
        a_prime = _hp_mp(u) * du
        b_prime = _hp_mp(t) * dt
        beta = ((_hp_mp(x) + _h_mp(x) / x) / p - a_prime) / (
            b_prime - a_prime)
        return MPParameters(*( +value for value in (x, p, mean, c, beta)))


class _SecondJet:
    """A scalar value with first and second derivatives, over Arb."""

    def __init__(self, value, first=0, second=0):
        self.value = arb(value)
        self.first = arb(first)
        self.second = arb(second)

    @staticmethod
    def coerce(value) -> "_SecondJet":
        return value if isinstance(value, _SecondJet) else _SecondJet(value)

    def __add__(self, other):
        other = self.coerce(other)
        return _SecondJet(
            self.value + other.value,
            self.first + other.first,
            self.second + other.second,
        )

    __radd__ = __add__

    def __neg__(self):
        return _SecondJet(-self.value, -self.first, -self.second)

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) + (-self)

    def __mul__(self, other):
        other = self.coerce(other)
        return _SecondJet(
            self.value * other.value,
            self.first * other.value + self.value * other.first,
            self.second * other.value
            + 2 * self.first * other.first
            + self.value * other.second,
        )

    __rmul__ = __mul__

    def reciprocal(self):
        return _SecondJet(
            1 / self.value,
            -self.first / self.value**2,
            (2 * self.first**2 - self.value * self.second) / self.value**3,
        )

    def __truediv__(self, other):
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other):
        return self.coerce(other) * self.reciprocal()


def _h_jet(x: _SecondJet) -> _SecondJet:
    value = _h_arb(x.value)
    first_entropy = _hp_arb(x.value)
    second_entropy = _hpp_arb(x.value)
    return _SecondJet(
        value,
        first_entropy * x.first,
        second_entropy * x.first**2 + first_entropy * x.second,
    )


def _core_jet(x: _SecondJet, beta: arb) -> _SecondJet:
    one = _SecondJet(1)
    protocol = x * x * (one + (one - x) * (one - x))
    return ((one - beta) * _h_jet(x * x) + beta * _h_jet(protocol)) / (
        x * _h_jet(x))


def certify_equation_parameters(mp_parameters: MPParameters) -> ArbParameters:
    """Bracket the algebraic x root and propagate it through Liu's equations."""
    scale = 10**ROOT_DECIMAL_PLACES
    with mpmath.workdps(ROOT_DECIMAL_PLACES + 30):
        scaled = mpmath.floor(mp_parameters.x * scale)
    root_lo = Fraction(int(scaled), scale)
    root_hi = Fraction(int(scaled) + 1, scale)
    lo = _arb_fraction(root_lo)
    hi = _arb_fraction(root_hi)
    x = lo.union(hi)

    def polynomial(z: arb) -> arb:
        return z**4 - 2 * z**3 + 3 * z**2 - 1

    derivative = 4 * x**3 - 6 * x**2 + 6 * x
    if not (polynomial(lo) < 0 < polynomial(hi)):
        raise AssertionError("the rational endpoints do not bracket Liu's x root")
    if not derivative > 0:
        raise AssertionError("the polynomial root was not proved unique in its bracket")

    hx = _h_arb(x)
    hxx = _h_arb(x * x)
    p = hx / hxx
    mean = p * x
    c = 1 - mean
    u = x * x
    protocol = _protocol_argument(x, x)
    du = 2 * x
    dt = 2 * x + 2 * x * (1 - x) ** 2 - 2 * x * x * (1 - x)
    a_prime = _hp_arb(u) * du
    b_prime = _hp_arb(protocol) * dt
    beta = ((_hp_arb(x) + hx / x) / p - a_prime) / (b_prime - a_prime)
    core_jet = _core_jet(_SecondJet(x, 1, 0), beta)
    if not core_jet.first.contains(0):
        raise AssertionError("Arb stationarity enclosure misses zero")
    return ArbParameters(
        x=x,
        p=p,
        mean=mean,
        c=c,
        beta=beta,
        core=core_jet.value,
        core_prime=core_jet.first,
        core_second=core_jet.second,
        root_lo=root_lo,
        root_hi=root_hi,
    )


def _f_xx_arb(x: arb, y: arb, beta: arb) -> arb:
    product = x * y
    j_xx = _hpp_arb(product) * y**2
    protocol = _protocol_argument(x, y)
    protocol_x = y * (1 + (1 - 2 * x) * (1 - y))
    protocol_xx = -2 * y * (1 - y)
    k_xx = (_hpp_arb(protocol) * protocol_x**2
            + _hp_arb(protocol) * protocol_xx)
    return (1 - beta) * j_xx + beta * k_xx


def _k_xy_arb(x: arb, y: arb) -> arb:
    protocol = _protocol_argument(x, y)
    protocol_x = y * (1 + (1 - 2 * x) * (1 - y))
    protocol_y = x * (1 + (1 - 2 * y) * (1 - x))
    protocol_xy = 1 + (1 - 2 * x) * (1 - 2 * y)
    return (_hpp_arb(protocol) * protocol_x * protocol_y
            + _hp_arb(protocol) * protocol_xy)


def certify_curvatures(parameters: ArbParameters) -> CurvatureCertificate:
    x, p, beta = parameters.x, parameters.p, parameters.beta
    face = parameters.mean * parameters.core_second

    x_jet = _SecondJet(x, 1, 0)
    diagonal_protocol = x_jet * x_jet * (
        1 + (1 - x_jet) * (1 - x_jet))
    protocol_diagonal_second = _h_jet(diagonal_protocol).second
    f_xx = _hpp_arb(x * x) * x**2
    split_unit = (
        2 * (1 - beta) * p * f_xx
        + beta * p * protocol_diagonal_second
        - _hpp_arb(x)
    ) / _h_arb(x)

    boundary_barrier = 2 * p * x * (1 + beta * (1 - x)) - 1
    gain = p * (2 * p * _f_xx_arb(x, x, beta) - _hpp_arb(x))
    loss = -2 * beta * p**2 * _k_xy_arb(x, x)
    if not (face > 0 and split_unit > 0 and boundary_barrier > 0
            and gain > 0 and gain - loss > 0):
        raise AssertionError("Arb did not certify the required positive curvatures")
    return CurvatureCertificate(
        face=face,
        split_unit=split_unit,
        boundary_barrier=boundary_barrier,
        gain=gain,
        loss=loss,
        loss_over_gain=loss / gain,
        gain_minus_loss=gain - loss,
    )

def _insertion_potential_arb(parameters: ArbParameters, y: arb) -> arb:
    """Reduced first variation for inserting one atom at y."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    entropy = p * _h_arb(x)
    lambda_gap = entropy * parameters.core
    protocol = _protocol_argument(x, y)
    return (
        2 * p * ((1 - beta) * h_arb(x * y) + beta * h_arb(protocol))
        - h_arb(y)
        - lambda_gap * y
    ) / entropy


def _insertion_potential_prime_arb(parameters: ArbParameters, y: arb) -> arb:
    x, p, beta = parameters.x, parameters.p, parameters.beta
    entropy = p * _h_arb(x)
    lambda_gap = entropy * parameters.core
    xy = x * y
    protocol = _protocol_argument(x, y)
    protocol_y = x * (1 + (1 - x) * (1 - 2 * y))
    return (
        2 * p * (
            (1 - beta) * _hp_arb(xy) * x
            + beta * _hp_arb(protocol) * protocol_y
        )
        - _hp_arb(y)
        - lambda_gap
    ) / entropy


def _insertion_potential_second_arb(parameters: ArbParameters, y: arb) -> arb:
    x, p, beta = parameters.x, parameters.p, parameters.beta
    entropy = p * _h_arb(x)
    xy = x * y
    protocol = _protocol_argument(x, y)
    protocol_y = x * (1 + (1 - x) * (1 - 2 * y))
    protocol_yy = -2 * x * (1 - x)
    return (
        2 * p * (
            (1 - beta) * _hpp_arb(xy) * x**2
            + beta * (
                _hpp_arb(protocol) * protocol_y**2
                + _hp_arb(protocol) * protocol_yy
            )
        )
        - _hpp_arb(y)
    ) / entropy


def _certify_mvt_segment(parameters: ArbParameters, start: Fraction,
                         end: Fraction, boxes: int) -> arb:
    """Prove V>0 on a closed segment by centered mean-value forms."""
    step = (end - start) / boxes
    weakest: arb | None = None
    weakest_lower = math.inf
    for index in range(boxes):
        lo = start + index * step
        hi = lo + step
        center = (lo + hi) / 2
        radius = step / 2
        interval = _arb_fraction(lo).union(_arb_fraction(hi))
        center_value = _insertion_potential_arb(
            parameters, _arb_fraction(center))
        derivative = _insertion_potential_prime_arb(parameters, interval)
        lower_form = center_value - abs(derivative) * _arb_fraction(radius)
        if not lower_form > 0:
            raise AssertionError(
                f"insertion-potential MVT failed on [{lo},{hi}]: {lower_form}")
        lower = float(lower_form.lower())
        if lower < weakest_lower:
            weakest, weakest_lower = lower_form, lower
    if weakest is None:
        raise AssertionError("empty insertion-potential MVT segment")
    return weakest


def certify_insertion_potential(
        parameters: ArbParameters) -> InsertionCertificate:
    """Prove the zero-mass atom insertion potential is nonnegative on [0,1]."""
    delta = Fraction(1, 10_000)
    delta_arb = _arb_fraction(delta)
    x, p, beta = parameters.x, parameters.p, parameters.beta
    entropy = p * _h_arb(x)
    lambda_gap = entropy * parameters.core
    # For 0<y<=delta, h(cy)>=cy*log(1/y) when 0<c<=1, while
    # h(y)<=y*log(1/y)+y.  The protocol argument is y*g(y), with g decreasing.
    g_delta = x * (1 + (1 - x) * (1 - delta_arb))
    logarithmic_coefficient = (
        2 * p * ((1 - beta) * x + beta * g_delta) - 1)
    small_y_margin = (
        logarithmic_coefficient * _arb_fraction(1 / delta).log()
        - 1
        - lambda_gap
    )
    if not small_y_margin > 0:
        raise AssertionError("small-y insertion bound was not positive")

    left_boxes = 15_000
    right_boxes = 7_000
    convexity_boxes = 100
    left_mvt = _certify_mvt_segment(
        parameters, delta, Fraction(17, 25), left_boxes)
    right_mvt = _certify_mvt_segment(
        parameters, Fraction(7, 10), Fraction(999, 1000), right_boxes)

    convexity_step = (Fraction(7, 10) - Fraction(17, 25)) / convexity_boxes
    weakest_convexity: arb | None = None
    weakest_convexity_lower = math.inf
    for index in range(convexity_boxes):
        lo = Fraction(17, 25) + index * convexity_step
        hi = lo + convexity_step
        interval = _arb_fraction(lo).union(_arb_fraction(hi))
        second = _insertion_potential_second_arb(parameters, interval)
        if not second > 0:
            raise AssertionError(
                f"insertion-potential convexity failed on [{lo},{hi}]: {second}")
        lower = float(second.lower())
        if lower < weakest_convexity_lower:
            weakest_convexity, weakest_convexity_lower = second, lower
    if weakest_convexity is None:
        raise AssertionError("empty insertion-potential convexity segment")

    endpoint_interval = _insertion_potential_arb(
        parameters,
        _arb_fraction(Fraction(999, 1000)).union(arb(1)),
    )
    if not endpoint_interval > 0:
        raise AssertionError("near-one insertion-potential interval was not positive")
    return InsertionCertificate(
        small_y_margin=small_y_margin,
        left_mvt_lower=left_mvt,
        root_convexity=weakest_convexity,
        right_mvt_lower=right_mvt,
        endpoint_interval=endpoint_interval,
        boxes=left_boxes + convexity_boxes + right_boxes,
    )


def _canonical_values(parameters: MPParameters, q: mpmath.mpf,
                      inactive: Sequence[mpmath.mpf] | None = None):
    r = (1 - parameters.p) / 2
    p0 = (parameters.x, mpmath.mpf(0), mpmath.mpf(0))
    p1 = p0 if inactive is None else tuple(inactive)
    return (parameters.p, r, q, *p0, *p1)


def print_parameterisation(mp_parameters: MPParameters,
                           arb_parameters: ArbParameters, dps: int) -> None:
    print("1. EQUATION-DEFINED BINDING VALUE AND LOCUS")
    print("PROVED [Liu (87)]: x is the unique root in the printed rational bracket")
    print(f"   {arb_parameters.root_lo} < x < {arb_parameters.root_hi},")
    print("   of x^4-2*x^3+3*x^2-1=0.")
    print("PROVED [Liu (87)--(90), direct algebra]:")
    print("   t=x^2[1+(1-x)^2]=1-x^2, p=h(x)/h(x^2), m=p*x, c=1-m;")
    print("   hence h(t)=h(x^2) and Phi=1 for every beta on the two-atom law")
    print("   P*=p*delta_x+(1-p)*delta_0.  The displayed beta formula imposes")
    print("   constrained stationarity rather than treating beta as a tenth variable.")
    print(f"PROVED [320-bit Arb]: x in {arb_parameters.x}")
    print(f"PROVED [320-bit Arb]: p in {arb_parameters.p}")
    print(f"PROVED [320-bit Arb]: c_equation in {arb_parameters.c}")
    print(f"PROVED [320-bit Arb]: beta_equation in {arb_parameters.beta}")
    print("NUMERICAL [paper/source decimal comparison]:")
    print("   c_printed-c_equation=" + mpmath.nstr(
        mpmath.mpf(LIU_C_DECIMAL) - mp_parameters.c, 18))
    print("   beta_source-beta_equation=" + mpmath.nstr(
        mpmath.mpf(LIU_BETA_DECIMAL) - mp_parameters.beta, 18))
    print("   The remaining sections use the equation-defined values; no false exact")
    print("   equality is asserted for either printed decimal.")
    print()

    print("PROVED [explicit coordinate strata, up to permutations]:")
    print("   Z(r): a=(p,r,1-p-r), 0<=r<=1-p, P*=(x,0,0);")
    print("   X(r): a=(r,p-r,1-p),   0<=r<=p,   P*=(x,x,0).")
    print("   For 0<q<1, set P0=P1=P* in either representation.")
    print("   With all masses positive this is a 2-dimensional (q,r) coordinate")
    print("   stratum.  At r=0 or an analogous zero-mass representation, its two")
    print("   unused support coordinates make a 3-dimensional gauge stratum.  Only q")
    print("   remains after quotienting repeated/zero-mass atom representations.")
    print("PROVED [endpoint enlargement]: at q=0, P1=(u1,u2,u3) is arbitrary;")
    print("   at q=1, P0 is arbitrary.  Each endpoint has maximal coordinate dimension")
    print("   four: (r,u1,u2,u3), or one unused active support plus the three inactive")
    print("   supports at a zero mass.  Permutations add discrete copies.")
    print("CONJECTURED [completeness]: these are all equality strata.  Proving that")
    print("   no disconnected equality component exists is essentially the global part")
    print("   of Hypothesis 2 and is not inferred from a local Hessian calculation.")
    print()

    q_values = tuple(mpmath.mpf(text) for text in ("0", "0.25", "0.5", "0.75", "1"))
    for q in q_values:
        values = _canonical_values(mp_parameters, q)
        terms = evaluate_mpmath(values, mp_parameters.beta, dps=dps)
        mean = mean_value(values)
        print(
            "NUMERICAL [mpmath %d dps, diagonal q=%s]: Phi-1=%s; mean-(1-c)=%s"
            % (
                dps,
                mpmath.nstr(q, 3),
                mpmath.nstr(terms.objective - 1, 8),
                mpmath.nstr(mean - mp_parameters.mean, 8),
            )
        )
    for inactive in (
        (mpmath.mpf("0.2"), mpmath.mpf("0.6"), mpmath.mpf("0.9")),
        (mpmath.mpf("0.8"), mpmath.mpf("0.4"), mpmath.mpf("0.1")),
    ):
        values = _canonical_values(mp_parameters, mpmath.mpf(0), inactive)
        terms = evaluate_mpmath(values, mp_parameters.beta, dps=dps)
        mean = mean_value(values)
        print(
            "NUMERICAL [mpmath %d dps, q=0 inactive P1=%s]: Phi-1=%s; "
            "mean-(1-c)=%s"
            % (
                dps,
                tuple(mpmath.nstr(value, 3) for value in inactive),
                mpmath.nstr(terms.objective - 1, 8),
                mpmath.nstr(mean - mp_parameters.mean, 8),
            )
        )
    print()


def print_reduced_hessians(mp_parameters: MPParameters,
                           arb_parameters: ArbParameters,
                           curvatures: CurvatureCertificate,
                           insertion: InsertionCertificate, dps: int) -> None:
    with mpmath.workdps(dps):
        core_prime = mpmath.diff(
            lambda z: _core_mp(z, mp_parameters.beta), mp_parameters.x)
        reduced_gradient = abs(mp_parameters.mean * core_prime)

    print("2. ALL LOCAL DIRECTIONS ON THE INTERIOR DIAGONAL STRATUM")
    print("PROVED [exact active-mean coordinates]: write")
    print("   x0=s-q*d, x1=s+(1-q)*d, p=m/s, and let r split the zero mass.")
    print("   Then the mean is identically m.  At (s,d)=(x,0), in coordinates")
    print("   (q,s,d,r), the reduced gradient is zero and the Hessian is")
    print("   diag(0,A,C*q*(1-q),0).")
    print("   The q null direction is exact P0=P1 invariance; the r null direction")
    print("   merely repartitions coincident zero atoms.")
    print("NUMERICAL [mpmath %d dps derivative check]: |reduced gradient|=%s"
          % (dps, mpmath.nstr(reduced_gradient, 8)))
    print(f"PROVED [Arb eigenvalue coefficient]: A in {curvatures.face}")
    print(f"PROVED [Arb eigenvalue coefficient]: C in {curvatures.split_unit}")
    print(f"PROVED [Arb stationarity enclosure]: core'(x) in {arb_parameters.core_prime}")
    for q_fraction in (Fraction(0), Fraction(1, 4), Fraction(1, 2),
                       Fraction(3, 4), Fraction(1)):
        split = curvatures.split_unit * _arb_fraction(
            q_fraction * (1 - q_fraction))
        if q_fraction in (0, 1):
            spectrum = f"[0, 0, 0, {curvatures.face}]"
        else:
            spectrum = f"[0, 0, {split}, {curvatures.face}]"
        print("PROVED [Arb reduced spectrum, q=%s]: %s"
              % (q_fraction, spectrum))
    print()

    print("PROVED [one-sided positive-mass box directions]: every zero-support atom")
    print("   activated into the domain contributes, to the denominator-cleared gap,")
    print("   a leading positive multiple of y*log(1/y).  Its common coefficient is")
    print(f"   B=2*p*x*[1+beta*(1-x)]-1 in {curvatures.boundary_barrier}.")
    print("   For 0<q<1 the component weights are positive, so this dominates every")
    print("   smooth O(y) or O(y^2) term.")
    print("PROVED [all zero-mass atom directions]: activating a zero-mass pair at")
    print("   (y0,y1) has reduced first variation (1-q)V(y0)+q*V(y1), where")
    print("   V(y)={2p[(1-beta)h(xy)+beta*h(xy(1+(1-x)(1-y)))]")
    print("         -h(y)-H(P*)*core*y}/H(P*).")
    print("   Direct algebra gives V(0)=V(x)=V'(x)=0.  Arb proves V(y)>=0 on [0,1],")
    print("   with equality only at 0 and x, as follows:")
    print(f"   (0,1e-4]: transformed small-y lower margin {insertion.small_y_margin};")
    print(f"   [1e-4,0.68]: weakest centered-MVT form {insertion.left_mvt_lower};")
    print(f"   [0.68,0.70]: weakest V'' enclosure {insertion.root_convexity} >0;")
    print(f"   [0.70,0.999]: weakest centered-MVT form {insertion.right_mvt_lower};")
    print(f"   [0.999,1]: direct range {insertion.endpoint_interval}.")
    print(f"   The interval proof used {insertion.boxes} exact-rational subboxes.")
    print("   Thus mass activation at unused support coordinates supplies no missing")
    print("   descent.  At q=0 or 1 the arbitrary inactive law is handled next.")
    print()


def _h_float_array(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    answer = np.zeros_like(values)
    interior = (values > 0) & (values < 1)
    selected = values[interior]
    answer[interior] = -(
        selected * np.log(selected) + (1 - selected) * np.log1p(-selected))
    return answer


class _EndpointFunctional:
    """The complete reduced inward-q derivative at a q=0 binding point."""

    def __init__(self, parameters: MPParameters):
        self.p = float(parameters.p)
        self.x = float(parameters.x)
        self.beta = float(parameters.beta)
        self.weights_p = np.array((self.p, 1 - self.p))
        self.points_p = np.array((self.x, 0.0))
        self.mean = self.p * self.x
        self.entropy_p = self._entropy(self.weights_p, self.points_p)
        self.j_pp = self._j_cross(
            self.weights_p, self.points_p, self.weights_p, self.points_p)
        self.k_pp = self._k(self.weights_p, self.points_p)
        self.lambda_phi = float(_core_mp(parameters.x, parameters.beta))
        self.lambda_gap = self.entropy_p * self.lambda_phi

    @staticmethod
    def _entropy(weights: np.ndarray, points: np.ndarray) -> float:
        return float(weights @ _h_float_array(points))

    @staticmethod
    def _j_cross(left_weights: np.ndarray, left_points: np.ndarray,
                 right_weights: np.ndarray, right_points: np.ndarray) -> float:
        return float(left_weights @ _h_float_array(
            np.outer(left_points, right_points)) @ right_weights)

    @staticmethod
    def _k(weights: np.ndarray, points: np.ndarray) -> float:
        product = np.outer(points, points)
        argument = product * (1 + np.outer(1 - points, 1 - points))
        return float(weights @ _h_float_array(argument) @ weights)

    def value(self, weights: np.ndarray, points: np.ndarray) -> float:
        derivative_gap = (
            2 * (1 - self.beta) * (
                self._j_cross(self.weights_p, self.points_p, weights, points)
                - self.j_pp)
            + self.beta * (self._k(weights, points) - self.k_pp)
            - (self._entropy(weights, points) - self.entropy_p)
        )
        derivative_mean = float(weights @ points) - self.mean
        return ((derivative_gap - self.lambda_gap * derivative_mean)
                / self.entropy_p)

    def zero_split(self, values: Sequence[float]) -> float:
        split, *points = values
        residual = 1 - self.p
        weights = np.array((self.p, residual * split,
                            residual * (1 - split)))
        return self.value(weights, np.asarray(points))

    def x_split(self, values: Sequence[float]) -> float:
        split, *points = values
        weights = np.array((self.p * split, self.p * (1 - split),
                            1 - self.p))
        return self.value(weights, np.asarray(points))


def _endpoint_value_mp(parameters: MPParameters, pattern: str,
                       values: Sequence[float], dps: int) -> mpmath.mpf:
    with mpmath.workdps(dps):
        split = mpmath.mpf(repr(float(values[0])))
        points = tuple(mpmath.mpf(repr(float(value))) for value in values[1:])
        p, x, beta = parameters.p, parameters.x, parameters.beta
        if pattern == "zero-mass split":
            weights = (p, (1 - p) * split, (1 - p) * (1 - split))
            p_points = (x, mpmath.mpf(0), mpmath.mpf(0))
        elif pattern == "x-mass split":
            weights = (p * split, p * (1 - split), 1 - p)
            p_points = (x, x, mpmath.mpf(0))
        else:
            raise ValueError(pattern)

        def entropy(ws, xs):
            return sum((w * _h_mp(z) for w, z in zip(ws, xs)), mpmath.mpf(0))

        def j_cross(wl, xl, wr, xr):
            return sum((
                wl[i] * wr[j] * _h_mp(xl[i] * xr[j])
                for i in range(3) for j in range(3)
            ), mpmath.mpf(0))

        def kval(ws, xs):
            return sum((
                ws[i] * ws[j] * _h_mp(_protocol_argument(xs[i], xs[j]))
                for i in range(3) for j in range(3)
            ), mpmath.mpf(0))

        entropy_p = entropy(weights, p_points)
        derivative_gap = (
            2 * (1 - beta) * (
                j_cross(weights, p_points, weights, points)
                - j_cross(weights, p_points, weights, p_points))
            + beta * (kval(weights, points) - kval(weights, p_points))
            - (entropy(weights, points) - entropy_p)
        )
        derivative_mean = sum((weights[i] * points[i] for i in range(3)),
                              mpmath.mpf(0)) - parameters.mean
        lambda_phi = _core_mp(x, beta)
        return +(derivative_gap / entropy_p - lambda_phi * derivative_mean)


def search_endpoint_directions(parameters: MPParameters, iterations: int,
                               population_size: int, seed: int,
                               dps: int) -> tuple[EndpointSearch, ...]:
    model = _EndpointFunctional(parameters)
    rng = np.random.default_rng(seed)
    searches: list[EndpointSearch] = []
    patterns: tuple[tuple[str, Callable[[Sequence[float]], float], np.ndarray], ...] = (
        (
            "zero-mass split",
            model.zero_split,
            np.array((0.5, model.x, 0.0, 0.0)),
        ),
        (
            "x-mass split",
            model.x_split,
            np.array((0.5, model.x, model.x, 0.0)),
        ),
    )
    for index, (name, function, equality) in enumerate(patterns):
        initial = rng.random((population_size, 4))
        initial[0] = equality
        initial[1] = equality
        initial[1, 0] = 0.125 if index == 0 else 0.875
        result = differential_evolution(
            function,
            bounds=((0.0, 1.0),) * 4,
            init=initial,
            seed=seed + index,
            maxiter=iterations,
            tol=1e-10,
            atol=1e-13,
            polish=False,
            workers=1,
            updating="immediate",
        )
        point = tuple(float(value) for value in result.x)
        searches.append(EndpointSearch(
            name=name,
            value_float=float(result.fun),
            value_mp=_endpoint_value_mp(parameters, name, point, dps),
            variables=point,
            evaluations=int(result.nfev),
        ))
    return tuple(searches)


def _endpoint_derivative_arb(parameters: ArbParameters, split: Fraction,
                             points: Sequence[Fraction]) -> arb:
    p, x, beta = parameters.p, parameters.x, parameters.beta
    split_arb = _arb_fraction(split)
    weights = (p, (1 - p) * split_arb, (1 - p) * (1 - split_arb))
    p_points = (x, arb(0), arb(0))
    q_points = tuple(_arb_fraction(point) for point in points)

    def entropy(ws, xs):
        return sum((w * _h_arb_endpoint_safe(z) for w, z in zip(ws, xs)), arb(0))

    def j_cross(wl, xl, wr, xr):
        return sum((
            wl[i] * wr[j] * _h_arb_endpoint_safe(xl[i] * xr[j])
            for i in range(3) for j in range(3)
        ), arb(0))

    def kval(ws, xs):
        return sum((
            ws[i] * ws[j] * _h_arb_endpoint_safe(
                _protocol_argument(xs[i], xs[j]))
            for i in range(3) for j in range(3)
        ), arb(0))

    entropy_p = entropy(weights, p_points)
    derivative_gap = (
        2 * (1 - beta) * (
            j_cross(weights, p_points, weights, q_points)
            - j_cross(weights, p_points, weights, p_points))
        + beta * (kval(weights, q_points) - kval(weights, p_points))
        - (entropy(weights, q_points) - entropy_p)
    )
    derivative_mean = sum((weights[i] * q_points[i] for i in range(3)), arb(0)) \
        - parameters.mean
    return derivative_gap / entropy_p - parameters.core * derivative_mean


def print_endpoint_analysis(parameters_mp: MPParameters,
                            parameters_arb: ArbParameters,
                            searches: Sequence[EndpointSearch]) -> None:
    print("3. q=0 AND q=1 ENDPOINT TANGENT CONES")
    print("PROVED [first variation, q=0; q=1 is symmetric]: for the inactive law Q,")
    print("   every inward direction, after using the active mean tangent, has")
    print("   D(Q)={2(1-beta)[J(P*,Q)-J(P*,P*)]+beta[K(Q,Q)-K(P*,P*)]")
    print("         -[H(Q)-H(P*)]}/H(P*) - core*[mean(Q)-m].")
    print("   Thus D(Q), not a sampled coordinate Hessian, covers every inward-q")
    print("   direction at every point of the four-dimensional endpoint stratum.")
    for split, points in (
        (Fraction(1, 2), (Fraction(1, 5), Fraction(3, 5), Fraction(9, 10))),
        (Fraction(1, 3), (Fraction(4, 5), Fraction(2, 5), Fraction(1, 10))),
    ):
        enclosure = _endpoint_derivative_arb(parameters_arb, split, points)
        if not enclosure > 0:
            raise AssertionError("sample endpoint derivative was not certified positive")
        print("PROVED [Arb sample endpoint gradient]: split=%s, Q=%s, D(Q) in %s"
              % (split, points, enclosure))
    evaluations = sum(search.evaluations for search in searches)
    print("NUMERICAL [one-core global differential-evolution search, %d evaluations]:"
          % evaluations)
    for search in searches:
        print("   %s: min float64 D=%+.6e; mpmath recomputation=%s; variables=%s; nfev=%d"
              % (
                  search.name,
                  search.value_float,
                  mpmath.nstr(search.value_mp, 12),
                  np.array2string(np.asarray(search.variables), precision=9,
                                  separator=","),
                  search.evaluations,
              ))
    print("NUMERICAL [search verdict]: no D(Q)<0 was found; every sharp candidate")
    print("   is a coordinate representation of Q=P*, where D(P*)=0 exactly.")
    print("CONJECTURED [global endpoint sign]: D(Q)>=0 on both four-dimensional")
    print("   endpoint strata, with equality only when Q represents P*.  The searches")
    print("   are strong local/global evidence, not an interval proof over all Q.")
    print("PROVED [endpoint face Hessian]: for an interior inactive Q, directions that")
    print("   keep q=0 have spectrum [0,0,0,0,A]: three inactive-support gauges, one")
    print("   repeated-zero-mass gauge, and the positive face eigenvalue A.  If D(Q)>0,")
    print("   the inward-q ray is not in the critical cone; at Q=P* it joins the exact")
    print("   diagonal null manifold and is governed by the split curvature above.")
    print()


def print_paired_retest(curvatures: CurvatureCertificate) -> None:
    print("4. q DIRECTION AND PAIRED-ATOM RETEST AT BINDING")
    print("PROVED [degenerate fibre algebra]: at P0=P1=P*, Phi is exactly independent")
    print("   of q on all of [0,1]; every pure-q derivative is zero.")
    print("PROVED [Arb, the only moving interior atom, unit displacement]:")
    print(f"   GAIN in {curvatures.gain}")
    print(f"   LOSS in {curvatures.loss}")
    print(f"   LOSS/GAIN in {curvatures.loss_over_gain}")
    print(f"   GAIN-LOSS in {curvatures.gain_minus_loss}")
    print("PROVED [comparison with the earlier non-binding probe]: LOSS/GAIN is <0.212")
    print("   here, not >4.2848.  The paired curvature is")
    print("   Phi''(0)=q*(1-q)*C >0 for every 0<q<1.")
    print("PROVED [Arb, q=1/2]: Phi''(0) in "
          + str(curvatures.split_unit / 4))
    print("   Hence the earlier paired-atom descent mechanism does NOT survive on the")
    print("   equation-defined binding manifold; it reverses sign and is increasing.")
    print()


def print_literal_decimal_check() -> None:
    """Separate a harmless but rigorous printed-decimal rounding issue."""
    c_printed = Fraction(LIU_C_DECIMAL)
    beta_source = Fraction(LIU_BETA_DECIMAL)
    mean_printed = 1 - c_printed
    # Thirty decimal places are ample because the source-beta crossing is
    # stationary; the x rounding error enters quadratically.
    x = Fraction("0.690787593925025407002335315473")
    p = mean_printed / x
    r = (1 - p) / 2
    values = (
        _arb_fraction(p),
        _arb_fraction(r),
        _arb_fraction(Fraction(1, 2)),
        _arb_fraction(x), arb(0), arb(0),
        _arb_fraction(x), arb(0), arb(0),
    )
    target = _arb_fraction(mean_printed)
    objective = evaluate_arb(values, _arb_fraction(beta_source)).objective
    mean_residual = mean_value(values) - target
    if objective is None or not objective < 1 or not mean_residual.contains(0):
        raise AssertionError("literal-decimal rational witness was not certified")

    print("5. PRINTED-DECIMAL CAVEAT (NOT THE EQUATION-DEFINED HYPOTHESIS)")
    print("PROVED [exact rationals + Arb]: if the 15-digit printed c and beta are")
    print("   interpreted literally as exact rationals, rounding alone gives a sub-1")
    print("   diagonal point with the mean exactly active:")
    print(f"   c={c_printed}; beta={beta_source}; q=1/2;")
    print(f"   x={x}; p={p}; a=(p,(1-p)/2,(1-p)/2);")
    print(f"   objective in {objective}; mean-(1-c) in {mean_residual}.")
    print("PROVED [scope]: this only says the rounded-up decimal 0.382709087918741 is")
    print("   about 5.97e-15 above the equation-defined value.  It is not a descent")
    print("   from the mathematical binding locus and is not a refutation of Liu's")
    print("   equation-defined Hypothesis 2.")
    print()


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--search-iterations", type=int,
                        default=DEFAULT_SEARCH_ITERATIONS)
    parser.add_argument("--search-population", type=int,
                        default=DEFAULT_SEARCH_POPULATION)
    parser.add_argument("--seed", type=int, default=230_608_824)
    args = parser.parse_args(argv)
    if args.dps < 70:
        parser.error("--dps must be at least 70")
    if args.search_iterations < 20:
        parser.error("--search-iterations must be at least 20")
    if args.search_population < 32:
        parser.error("--search-population must be at least 32")
    mpmath.mp.dps = max(mpmath.mp.dps, args.dps + 30)


    print("LIU HYPOTHESIS-2 BINDING-LOCUS ANALYSIS")
    print("Status labels are PROVED, NUMERICAL, or CONJECTURED.")
    print("PROVED [runtime configuration]: BLAS/OpenMP limits and SciPy workers are one.")
    print()

    mp_parameters = solve_equation_parameters(args.dps + 20)
    arb_parameters = certify_equation_parameters(mp_parameters)
    curvatures = certify_curvatures(arb_parameters)
    insertion = certify_insertion_potential(arb_parameters)
    searches = search_endpoint_directions(
        mp_parameters,
        args.search_iterations,
        args.search_population,
        args.seed,
        args.dps,
    )

    print_parameterisation(mp_parameters, arb_parameters, args.dps)
    print_reduced_hessians(
        mp_parameters, arb_parameters, curvatures, insertion, args.dps)
    print_endpoint_analysis(mp_parameters, arb_parameters, searches)
    print_paired_retest(curvatures)
    print_literal_decimal_check()

    print("FINAL VERDICT")
    print("PROVED [local equation-defined binding manifold]: every feasible tangent")
    print("   direction is an exact q/representation null, has positive Hessian")
    print("   curvature, meets the positive y*log(1/y) box-face barrier, or has the")
    print("   globally nonnegative Arb-certified atom-insertion potential V.")
    print("PROVED [paired-atom crux]: LOSS/GAIN<0.212 and Phi''(0)>0 at q=1/2;")
    print("   the earlier LOSS/GAIN>4.2848 descent does not survive at binding.")
    print("NUMERICAL [all endpoint inward rays]: the explicit four-variable D(Q)")
    print("   searches found no negative direction, with equality only at Q=P*.")
    print("CONJECTURED [critical manifold]: modulo coordinate gauges and permutations,")
    print("   the complete binding set is the one-dimensional diagonal q-family P0=P1=P*,")
    print("   enlarged in nine-coordinate space to 4D at q=0 and q=1 by the inactive law.")
    print("CONJECTURED [Hypothesis-2 verdict]: NO DESCENT FOUND and no certified")
    print("   refutation.  The local tangent cone is certified non-decreasing, but the")
    print("   global completeness of the binding locus and global endpoint sign D(Q)>=0")
    print("   remain the exact statements a proof must control.")


if __name__ == "__main__":
    main()
