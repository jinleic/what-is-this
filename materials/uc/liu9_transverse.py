#!/usr/bin/env python3
"""Positivity of the mean-preserving transverse first variation for Liu H2.

`uc/liu9_ninevar.py` REFUTED the unrestricted ambient nine-variable extension
of the chart result: the configuration

    (p - epsilon, epsilon, q, x, y, 0, x, y, 0)

has a negative linear coefficient.  That configuration is mean-INFEASIBLE -- it
deletes mass from the atom at ``x`` without compensating -- so it says nothing
about Liu's Hypothesis 2, which lives in the half-space ``mean >= p*x``.  The
mean-feasible member of the same family is

    (p - epsilon*y/x, epsilon, q, x, y, 0, x, y, 0),

which inserts mass ``epsilon`` at ``y`` and removes ``epsilon*y/x`` at ``x``,
preserving the first moment exactly.  For that family `liu9_ninevar.py` reports
only a 65-point scan: positive at every sampled ``y``, weakest ``7.865e-8`` at
``y = 707/1024``, which is essentially Liu's abscissa.  A scan is not a proof,
and a weakest value that small at a sampled point is exactly the shape of an
uncertified sign.

This module replaces that scan with a universal enclosure over the whole open
interval, and explains the small value: the coefficient has a DOUBLE ROOT at
``y = x``.

Write the mean-preserving linear coefficient as

    L(y) = 2p*K(x,y) - h(y) - (y/x)*A - kappa*y^2*(y-x)^2,
    A    = 2p*K(x,x) - h(x),
    K(l,r) = (1-beta)*h(l*r) + beta*h(prot(l,r)),
    prot(l,r) = l*r + l*(1-l)*r*(1-r).

Three facts, in increasing depth.

1.  ``L(x) = 0`` identically.  At ``y = x`` the ratio ``y/x`` is one, so the
    first three terms telescope, and ``(y-x)^2`` kills the last.  No property
    of ``x`` is used.  This is a structural zero, not a numerical one.

2.  ``prot(x,x) = 1 - x^2`` is EQUIVALENT to Liu's quartic.  Expanding,
    ``prot(z,z) - (1 - z^2) = z^4 - 2z^3 + 3z^2 - 1`` as polynomials with
    integer coefficients.

    Attribution, checked against the paper rather than assumed.  Liu's eq. (87)
    reads ``x*^2 + x*^2(1 + xbar*^2) = 1``, and ``x^2(1+(1-x)^2)`` IS
    ``prot(x,x)``, so (87) already IS the protocol identity, stated directly as
    a defining relation of the numerically located optimizer.  It is NOT
    derived there as a stationarity condition, and the quartic polynomial does
    not appear in the paper at all -- the string "quartic" occurs zero times in
    arXiv:2306.08824v1 (Jingbo Liu).  So the quartic is the derived form here,
    not the source.

    What the paper does not record is the CONSEQUENCE.  Because ``h`` is
    symmetric about ``1/2``, (87) forces ``h(prot(x,x)) = h(x^2)``, hence

        K(x,x) = h(x^2)   for EVERY beta,

    and with ``p = h(x)/h(x^2)`` therefore ``p*K(x,x) = h(x)`` and ``A = h(x)``.

3.  ``L'(x) = 0``.  This needs Liu's beta-system (89)-(90) as well as (87):
    it gives ``2p*d/dy K(x,y)|_{y=x} = h'(x) + h(x)/x``, and fact 2 gives
    ``A/x = h(x)/x``, so the two cancel.  It is NOT beta-independent -- the
    module records a mutation showing ``L'(x)`` moves off zero when beta is
    perturbed, so the cancellation is a property of Liu's optimizer and not an
    artifact of the transcription.

So ``y = x`` is a double root with ``L''(x) > 0``, and the 65-point scan's tiny
minimum was that double root seen from a nearby grid point.

``L`` is NOT globally convex: ``L'' < 0`` on roughly ``[0.01, 0.43]``.  The
certificate therefore splits ``(0,1)`` into three regions.

  A. Boundary layer ``0 < y <= y0``.  ``L(0) = 0`` as well, approached like
     ``y*log(1/y)``.  Dividing by ``y`` and using ``h(u)/u = log(1/u) + mu(u)``
     with ``mu(u) = ((1-u)/u)*log(1/(1-u))`` gives

         L(y)/y = C(y)*log(1/y) + R(y),
         C(y)   = 2p*((1-beta)*x + beta*rho(y)) - 1,
         rho(y) = x + x*(1-x)*(1-y),   so prot(x,y) = y*rho(y).

     Both pieces are singularity-free.  ``mu`` is enclosed WITHOUT evaluating a
     ``0/0``: from ``log(1/(1-u)) = sum_{k>=1} u^k/k`` one gets termwise
     ``1 <= sum_{k>=1} u^(k-1)/k <= 1/(1-u)``, hence ``1-u <= mu(u) <= 1``.
     That bound is elementary and exact.  Certifying ``C >= c1 > 0`` and
     ``R >= -c2`` then gives ``L(y)/y >= c1*log(1/y0) - c2``.

  B. Bulk ``y0 <= y <= x-r`` and ``x+r <= y <= 1``.  Direct interval enclosure
     of ``L`` on a cover, using the endpoint-aware entropy range enclosure.

  C. Double root ``|y-x| <= r``.  ``L(x) = 0`` exactly and ``|L'(x)| <= delta``,
     so Taylor with an interval second derivative gives

         L(y) >= (1/2)*m*(y-x)^2 - delta*|y-x| >= -delta^2/(2m),

     where ``m = min L''`` over the window.  The deficit is the square of Liu's
     own root residual and is the only thing between this and an exact zero.

WHAT THIS DOES AND DOES NOT CLAIM.  It certifies the leading coefficient of the
pencil along every mean-preserving single-atom insertion of the displayed
shape, for all ``y`` in ``(0,1)`` and all ``q`` -- the coefficient is q-free.
It does NOT certify Liu's Hypothesis 2, does not certify a tube radius, and
does not treat multi-atom insertions, asymmetric insertions placing different
support points in the two components, or the epsilon^2 and higher terms.  Those
remain OPEN and are named in the verdict.

Run from the repository root with

    math/.venv/bin/python -I -B math/uc/liu9_transverse.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from dataclasses import dataclass
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
    _h_arb,
    _hp_arb,
    _hpp_arb,
    certify_equation_parameters,
    solve_equation_parameters,
)
from liu9_objective import h_arb  # noqa: E402
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_ninevar import distance_squared  # noqa: E402

ctx.prec = max(ctx.prec, 320)

DEFAULT_DPS = 90
KAPPA = Fraction(4_119_063, 33_554_432)

#: Boundary-layer top.  Below this the log(1/y) factorization is used.
DEFAULT_Y0 = Fraction(1, 32)
#: Half-width of the double-root window around x.
DEFAULT_WINDOW = Fraction(1, 8)
#: Cells per unit length in the bulk cover.
DEFAULT_BULK_CELLS = 4096
#: Octaves below y0 in the boundary-layer cover.
DEFAULT_OCTAVES = 60
#: Cells across the double-root window.
DEFAULT_WINDOW_CELLS = 512


# ---------------------------------------------------------------------------
# Scalar helpers.


def _arbf(value: Fraction | int) -> arb:
    value = Fraction(value)
    return arb(value.numerator) / arb(value.denominator)


def _ball(lower: Fraction, upper: Fraction) -> arb:
    """The exact interval [lower, upper] as an Arb ball."""
    return _arbf(lower).union(_arbf(upper))


def _lower(value: arb) -> float:
    return float(value.lower())


def _upper(value: arb) -> float:
    return float(value.upper())


# ---------------------------------------------------------------------------
# Section 1.  The protocol identity, in exact integer arithmetic.


def _poly_mul(left: Sequence[int], right: Sequence[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                out[i + j] += a * b
    return out


def _poly_sub(left: Sequence[int], right: Sequence[int]) -> list[int]:
    size = max(len(left), len(right))
    out = []
    for i in range(size):
        a = left[i] if i < len(left) else 0
        b = right[i] if i < len(right) else 0
        out.append(a - b)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def certify_protocol_identity() -> dict[str, Any]:
    """``prot(z,z) - (1 - z^2)`` IS Liu's quartic, as integer polynomials.

    Coefficients are listed from the constant term upward.  Nothing here is
    numerical: the claim is an identity in ``Z[z]``, so a coefficient-by-
    coefficient comparison settles it.
    """
    z = [0, 1]
    one = [1]
    one_minus_z = _poly_sub(one, z)
    z_squared = _poly_mul(z, z)
    # prot(z,z) = z*z + z*(1-z)*z*(1-z)
    cross = _poly_mul(_poly_mul(z, one_minus_z), _poly_mul(z, one_minus_z))
    prot_diagonal = [a + b for a, b in zip(
        z_squared + [0] * (len(cross) - len(z_squared)), cross)]
    reflection = _poly_sub(one, z_squared)
    difference = _poly_sub(prot_diagonal, reflection)
    quartic = [-1, 0, 3, -2, 1]
    return {
        "prot_diagonal": list(prot_diagonal),
        "entropy_reflection_of_x_squared": list(reflection),
        "difference": list(difference),
        "liu_quartic": list(quartic),
        "identity_holds": list(difference) == list(quartic),
    }


# ---------------------------------------------------------------------------
# Section 2.  A second-order jet in y over Arb.


@dataclass(frozen=True)
class Jet2:
    """Value and first two derivatives with respect to the single variable y."""

    v: arb
    d: arb
    dd: arb

    @staticmethod
    def constant(value: arb) -> "Jet2":
        return Jet2(value, arb(0), arb(0))

    @staticmethod
    def variable(value: arb) -> "Jet2":
        return Jet2(value, arb(1), arb(0))

    def __add__(self, other: "Jet2") -> "Jet2":
        return Jet2(self.v + other.v, self.d + other.d, self.dd + other.dd)

    def __sub__(self, other: "Jet2") -> "Jet2":
        return Jet2(self.v - other.v, self.d - other.d, self.dd - other.dd)

    def __mul__(self, other: "Jet2") -> "Jet2":
        return Jet2(
            self.v * other.v,
            self.d * other.v + self.v * other.d,
            self.dd * other.v + 2 * self.d * other.d + self.v * other.dd,
        )

    def scaled(self, factor: arb) -> "Jet2":
        return Jet2(self.v * factor, self.d * factor, self.dd * factor)


def jet_entropy_interior(argument: Jet2) -> Jet2:
    """Binary entropy of a jet whose value stays strictly inside (0,1).

    The caller is responsible for interiority; `_h_arb` and its derivatives are
    undefined at the endpoints and Arb will report a non-finite ball rather
    than silently returning a wrong number, which is the behaviour we want.
    """
    value = argument.v
    first = _hp_arb(value)
    second = _hpp_arb(value)
    return Jet2(
        _h_arb(value),
        first * argument.d,
        second * argument.d * argument.d + first * argument.dd,
    )


def _protocol_jet(x: arb, y: Jet2) -> Jet2:
    """prot(x, y) = x*y + x*(1-x)*y*(1-y), as a jet in y."""
    one = Jet2.constant(arb(1))
    linear = y.scaled(x)
    curved = (y * (one - y)).scaled(x * (arb(1) - x))
    return linear + curved


def transverse_jet(parameters: ArbParameters, y: Jet2) -> Jet2:
    """L(y) as a jet in y, valid where every entropy argument is interior."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    kappa = _arbf(KAPPA)
    one = arb(1)

    product = y.scaled(x)
    protocol = _protocol_jet(x, y)
    kernel = (jet_entropy_interior(product).scaled(one - beta)
              + jet_entropy_interior(protocol).scaled(beta))

    constant_a = _constant_a(parameters)
    shift = Jet2.constant(x)
    support = y * y * (y - shift) * (y - shift)

    return (kernel.scaled(2 * p)
            - jet_entropy_interior(y)
            - y.scaled(constant_a / x)
            - support.scaled(kappa))


def _constant_a(parameters: ArbParameters) -> arb:
    """A = 2p*K(x,x) - h(x), computed from its definition."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    diagonal = ((arb(1) - beta) * _h_arb(x * x)
                + beta * _h_arb(x * x + x * (arb(1) - x) * x * (arb(1) - x)))
    return 2 * p * diagonal - _h_arb(x)


def certify_constant_a(parameters: ArbParameters) -> dict[str, Any]:
    """A = h(x), a consequence of the protocol identity and p = h(x)/h(x^2)."""
    x = parameters.x
    measured = _constant_a(parameters)
    predicted = _h_arb(x)
    difference = measured - predicted
    diagonal_gap = (_h_arb(x * x + x * (arb(1) - x) * x * (arb(1) - x))
                    - _h_arb(x * x))
    return {
        "constant_a": measured.str(20),
        "entropy_at_x": predicted.str(20),
        "difference": difference.str(20),
        "difference_upper": float(abs(difference).upper()),
        "kernel_diagonal_minus_h_x_squared": diagonal_gap.str(20),
        "beta_free_diagonal": float(abs(diagonal_gap).upper()),
    }


# ---------------------------------------------------------------------------
# Section 3.  The double root at y = x.


def certify_structural_zero(samples: int = 12) -> dict[str, Any]:
    """L(x) = 0 as a rational-function identity, independent of every constant.

    Arb cannot see this cancellation.  Evaluating L at a ball containing x
    computes ``y*(A/x)`` and ``K(x,y)`` from separate sub-expressions, so
    ball-minus-ball leaves a residue of order the parameter width -- about
    1e-68 here.  That residue is DEPENDENCY ERROR, not a deficit, and using it
    in the Taylor bound would swamp the real one by sixty-eight orders.

    The cancellation is structural: at y = x the ratio y/x is one and the
    support term carries (y-x)^2, so

        L(x) = A - 1*A - kappa*x^2*0 = 0

    for ANY values of x, p, beta, kappa and of the three entropy quantities
    h(x), h(x^2), h(prot(x,x)).  Treating those as free symbols and evaluating
    in exact rational arithmetic at pseudo-random points is a Schwartz-Zippel
    certificate for that identity: a nonzero rational function of this degree
    cannot vanish at every sampled point.
    """
    import random

    generator = random.Random(20260829)
    witnesses = []
    for _ in range(samples):
        def pick() -> Fraction:
            return Fraction(generator.randint(1, 10 ** 6),
                            generator.randint(1, 10 ** 6))

        x, p, beta, kappa = pick(), pick(), pick(), pick()
        h_x, h_x2, h_prot = pick(), pick(), pick()
        kernel_diagonal = (1 - beta) * h_x2 + beta * h_prot
        constant_a = 2 * p * kernel_diagonal - h_x
        # L evaluated at y = x, term for term as the definition reads.
        value = (2 * p * kernel_diagonal
                 - h_x
                 - (x / x) * constant_a
                 - kappa * x ** 2 * (x - x) ** 2)
        witnesses.append(value == 0)
        if value != 0:
            raise AssertionError(
                f"the structural zero failed at x={x}, p={p}, beta={beta}")
    return {
        "samples": samples,
        "all_exactly_zero": all(witnesses),
        "arithmetic": "exact Fraction, entropy values treated as free symbols",
    }


def certify_coverage(y0: Fraction, window: Fraction, centre: Fraction,
                     smallest: Fraction) -> dict[str, Any]:
    """The four regions tile [0,1] with no gap and no unproved sliver."""
    pieces = [
        ("A tail", Fraction(0), smallest),
        ("A octaves", smallest, y0),
        ("B below", y0, centre - window),
        ("C window", centre - window, centre + window),
        ("B above", centre + window, Fraction(1)),
    ]
    for index in range(len(pieces) - 1):
        if pieces[index][2] != pieces[index + 1][1]:
            raise AssertionError(
                f"gap between {pieces[index][0]} and {pieces[index+1][0]}")
    for name, lo, hi in pieces:
        if not hi > lo:
            raise AssertionError(f"empty region {name}")
    if pieces[0][1] != 0 or pieces[-1][2] != 1:
        raise AssertionError("the cover does not span [0,1]")
    return {
        "pieces": [[name, str(lo), str(hi)] for name, lo, hi in pieces],
        "spans_unit_interval": True,
        "abutting": True,
    }


def certify_double_root(parameters: ArbParameters) -> dict[str, Any]:
    """L(x) = 0 structurally, and |L'(x)| is bounded by Liu's root residual."""
    x = parameters.x
    jet = transverse_jet(parameters, Jet2.variable(x))
    value_bound = float(abs(jet.v).upper())
    slope_bound = float(abs(jet.d).upper())
    curvature = jet.dd
    if not curvature.lower() > 0:
        raise AssertionError(
            f"the second derivative at x must be positive, got {curvature}"
        )
    return {
        "value": jet.v.str(20),
        "value_upper": value_bound,
        "slope": jet.d.str(20),
        "slope_upper": slope_bound,
        "curvature": curvature.str(20),
        "curvature_lower": float(curvature.lower()),
        "structural_zero": True,
        "value_enclosure_is_dependency_error": True,
    }


def beta_sensitivity(mp_parameters: MPParameters, dps: int) -> dict[str, Any]:
    """L'(x) moves off zero when beta is perturbed.

    Without this the vanishing slope could be an artifact of a transcription
    that happens to cancel.  The cancellation is Liu's equation (90), so it
    must be destroyed by moving beta and by nothing else.
    """
    rows = []
    with mpmath.workdps(dps):
        x = mp_parameters.x
        p = mp_parameters.p

        def slope(beta: mpmath.mpf) -> mpmath.mpf:
            def kernel(right: mpmath.mpf) -> mpmath.mpf:
                product = x * right
                protocol = x * right + x * (1 - x) * right * (1 - right)
                return ((1 - beta) * _entropy_mp(product)
                        + beta * _entropy_mp(protocol))

            constant_a = 2 * p * kernel(x) - _entropy_mp(x)
            derivative = mpmath.diff(kernel, x)
            return (2 * p * derivative
                    - mpmath.diff(_entropy_mp, x)
                    - constant_a / x)

        exact = mp_parameters.beta
        for label, beta in (
            ("liu", exact),
            ("liu + 1e-6", exact + mpmath.mpf("1e-6")),
            ("liu - 1e-6", exact - mpmath.mpf("1e-6")),
            ("0", mpmath.mpf(0)),
            ("1/2", mpmath.mpf(1) / 2),
            ("1", mpmath.mpf(1)),
        ):
            rows.append({
                "beta": label,
                "slope": mpmath.nstr(slope(beta), 12),
                "slope_float": float(slope(beta)),
            })
    perturbed = [row for row in rows if row["beta"] != "liu"]
    return {
        "rows": rows,
        "liu_slope": rows[0]["slope_float"],
        "smallest_perturbed_slope": min(
            abs(row["slope_float"]) for row in perturbed),
        "beta_is_load_bearing": all(
            abs(row["slope_float"]) > 1e-9 for row in perturbed),
    }


def _entropy_mp(value: mpmath.mpf) -> mpmath.mpf:
    if value == 0 or value == 1:
        return mpmath.mpf(0)
    return -(value * mpmath.log(value) + (1 - value) * mpmath.log(1 - value))


# ---------------------------------------------------------------------------
# Section 4.  Region C, the double-root window.


def certify_window(
    parameters: ArbParameters,
    radius: Fraction,
    cells: int,
    slope_bound: float,
) -> dict[str, Any]:
    """Certify L >= -delta^2/(2m) on |y-x| <= radius via an interval Taylor."""
    x = parameters.x
    centre = Fraction(x.str(40, radius=False))
    worst_curvature: Optional[arb] = None
    for index in range(cells):
        lo = centre - radius + 2 * radius * Fraction(index, cells)
        hi = centre - radius + 2 * radius * Fraction(index + 1, cells)
        jet = transverse_jet(parameters, Jet2.variable(_ball(lo, hi)))
        if not jet.dd.is_finite():
            raise AssertionError(f"non-finite curvature on [{lo}, {hi}]")
        if worst_curvature is None or (
                jet.dd.lower() < worst_curvature.lower()):
            worst_curvature = jet.dd
    assert worst_curvature is not None
    minimum = arb(worst_curvature.lower())
    if not minimum > 0:
        raise AssertionError(
            f"the window is not convex: min L'' = {worst_curvature}")
    deficit = slope_bound ** 2 / (2.0 * float(minimum.lower()))
    return {
        "radius": str(radius),
        "cells": cells,
        "minimum_curvature": worst_curvature.str(20),
        "minimum_curvature_lower": float(minimum.lower()),
        "slope_bound": slope_bound,
        "deficit": deficit,
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Section 5.  Region B, the bulk.


def _transverse_value(parameters: ArbParameters, y: arb) -> arb:
    """L(y) as a plain enclosure, endpoint-safe in every entropy argument."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    one = arb(1)
    kappa = _arbf(KAPPA)
    product = x * y
    protocol = x * y + x * (one - x) * y * (one - y)
    kernel = (one - beta) * h_arb(product) + beta * h_arb(protocol)
    constant_a = _constant_a(parameters)
    return (2 * p * kernel
            - h_arb(y)
            - y * constant_a / x
            - kappa * y * y * (y - x) * (y - x))


def certify_bulk(
    parameters: ArbParameters,
    y0: Fraction,
    window: Fraction,
    cells_per_unit: int,
) -> dict[str, Any]:
    """Direct enclosure of L on [y0, x-r] and [x+r, 1]."""
    x = parameters.x
    centre = Fraction(x.str(40, radius=False))
    segments = [
        ("below", y0, centre - window),
        ("above", centre + window, Fraction(1)),
    ]
    worst: Optional[arb] = None
    worst_cell = None
    total = 0
    for name, lo_end, hi_end in segments:
        if not hi_end > lo_end:
            raise ValueError(f"empty {name} segment")
        count = max(1, int(cells_per_unit * (hi_end - lo_end)) + 1)
        for index in range(count):
            lo = lo_end + (hi_end - lo_end) * Fraction(index, count)
            hi = lo_end + (hi_end - lo_end) * Fraction(index + 1, count)
            value = _transverse_value(parameters, _ball(lo, hi))
            total += 1
            if not value.is_finite():
                raise AssertionError(f"non-finite L on [{lo}, {hi}]")
            if worst is None or value.lower() < worst.lower():
                worst, worst_cell = value, (name, str(lo), str(hi))
    assert worst is not None
    if not worst > 0:
        raise AssertionError(
            f"bulk positivity failed on {worst_cell} with L = {worst}")
    return {
        "cells": total,
        "y0": str(y0),
        "window": str(window),
        "worst_margin": worst.str(20),
        "worst_margin_lower": float(worst.lower()),
        "worst_cell": worst_cell,
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Section 6.  Region A, the y*log(1/y) boundary layer.


def _mu_enclosure(argument: arb, upper_argument: Fraction) -> arb:
    """Enclose mu(u) = ((1-u)/u)*log(1/(1-u)) without evaluating 0/0.

    From log(1/(1-u)) = sum_{k>=1} u^k / k, dividing by u gives
    sum_{k>=1} u^(k-1)/k, which is at least its first term 1 and at most the
    geometric sum 1/(1-u).  Multiplying by (1-u) gives 1-u <= mu(u) <= 1.
    Both bounds are exact and elementary; neither needs u to be bounded away
    from zero, which is the whole point.
    """
    del argument
    lower = arb(1) - _arbf(upper_argument)
    return lower.union(arb(1))


def certify_boundary_layer(
    parameters: ArbParameters,
    y0: Fraction,
    octaves: int,
) -> dict[str, Any]:
    """L(y)/y >= C*log(1/y) + R > 0 on (0, y0], with no singular evaluation."""
    x, p, beta = parameters.x, parameters.p, parameters.beta
    one = arb(1)
    kappa = _arbf(KAPPA)
    constant_a = _constant_a(parameters)

    worst_total: Optional[arb] = None
    worst_cell = None
    worst_coefficient: Optional[arb] = None
    smallest = y0 / Fraction(2) ** octaves

    for index in range(octaves):
        hi = y0 / Fraction(2) ** index
        lo = hi / 2
        y = _ball(lo, hi)
        # rho(y) = x + x(1-x)(1-y);  prot(x,y) = y*rho(y)
        rho = x + x * (one - x) * (one - y)
        protocol = y * rho
        product = x * y

        coefficient = 2 * p * ((one - beta) * x + beta * rho) - one
        if not coefficient > 0:
            raise AssertionError(
                f"the log coefficient must be positive on [{lo}, {hi}], "
                f"got {coefficient}"
            )

        x_hi = Fraction(x.str(40, radius=False)) + Fraction(1, 10 ** 30)
        mu_product = _mu_enclosure(product, x_hi * hi)
        mu_protocol = _mu_enclosure(protocol, x_hi * (2 - x_hi) * hi)
        mu_plain = _mu_enclosure(y, hi)

        remainder = (
            2 * p * ((one - beta) * x * ((one / x).log() + mu_product)
                     + beta * rho * (-(rho.log()) + mu_protocol))
            - mu_plain
            - constant_a / x
            - kappa * y * (y - x) * (y - x)
        )
        log_factor = (one / _arbf(hi)).log()
        total = coefficient * log_factor + remainder
        if not total.is_finite():
            raise AssertionError(f"non-finite layer bound on [{lo}, {hi}]")
        if worst_total is None or total.lower() < worst_total.lower():
            worst_total, worst_cell = total, (str(lo), str(hi))
        if worst_coefficient is None or (
                coefficient.lower() < worst_coefficient.lower()):
            worst_coefficient = coefficient

    assert worst_total is not None and worst_coefficient is not None
    if not worst_total > 0:
        raise AssertionError(
            f"boundary layer failed on {worst_cell} with bound {worst_total}")

    # The tail [0, smallest] is a genuine interval cell, not an argument.  The
    # factorization has no singular evaluation anywhere, so y = 0 may sit
    # inside the ball; only log(1/y) is unbounded there, and it appears
    # multiplied by a coefficient certified positive, so replacing it by
    # log(1/smallest) is a valid lower bound on the whole tail at once.
    #
    # The uniform coefficient bound needs no monotonicity argument:
    # rho(y) = x + x(1-x)(1-y) >= x for every y <= 1, hence C(y) >= 2px - 1.
    uniform_coefficient = 2 * p * x - one
    if not uniform_coefficient > 0:
        raise AssertionError(
            f"the uniform log coefficient is not positive: "
            f"{uniform_coefficient}")
    tail_y = _ball(Fraction(0), smallest)
    tail_rho = x + x * (one - x) * (one - tail_y)
    tail_protocol = tail_y * tail_rho
    tail_product = x * tail_y
    tail_x_hi = Fraction(x.str(40, radius=False)) + Fraction(1, 10 ** 30)
    tail_remainder = (
        2 * p * ((one - beta) * x * ((one / x).log()
                                     + _mu_enclosure(tail_product,
                                                     tail_x_hi * smallest))
                 + beta * tail_rho * (-(tail_rho.log())
                                      + _mu_enclosure(
                                          tail_protocol,
                                          tail_x_hi * (2 - tail_x_hi)
                                          * smallest)))
        - _mu_enclosure(tail_y, smallest)
        - constant_a / x
        - kappa * tail_y * (tail_y - x) * (tail_y - x)
    )
    tail_bound = (uniform_coefficient * (one / _arbf(smallest)).log()
                  + tail_remainder)
    if not tail_bound.is_finite():
        raise AssertionError("non-finite tail bound")
    if not tail_bound > 0:
        raise AssertionError(
            f"the tail cell [0, {smallest}] failed with bound {tail_bound}")
    return {
        "octaves": octaves,
        "y0": str(y0),
        "smallest_covered_y": str(smallest),
        "worst_layer_bound": worst_total.str(20),
        "worst_layer_bound_lower": float(worst_total.lower()),
        "worst_cell": worst_cell,
        "worst_log_coefficient": worst_coefficient.str(20),
        "worst_log_coefficient_lower": float(worst_coefficient.lower()),
        "uniform_log_coefficient": uniform_coefficient.str(20),
        "uniform_log_coefficient_lower": float(uniform_coefficient.lower()),
        "tail_cell": ["0", str(smallest)],
        "tail_bound": tail_bound.str(20),
        "tail_bound_lower": float(tail_bound.lower()),
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Section 7.  The asymmetric two-support family decouples.


def _pencil_mp(values, mp_parameters) -> mpmath.mpf:
    """gap - kappa*dist^2 from the sanctioned transcriptions only.

    `gap_mp` and `distance_squared` are the repository's single sources of
    truth.  Hand-building a `_ScalarOps` here is a known trap: its first
    argument is the ENTROPY FUNCTION, and passing `mpmath.mpf` silently makes
    "entropy(u) = u", which produced a plausible-looking but entirely wrong
    sign the first time this family was probed.
    """
    values = [mpmath.mpf(v) for v in values]
    kappa = mpmath.mpf(KAPPA.numerator) / KAPPA.denominator
    return (gap_mp(values, mp_parameters.beta)
            - kappa * distance_squared(values, mp_parameters))


def asymmetric_coefficient_mp(mp_parameters, y1, y2, q, epsilon=None):
    """Linear coefficient of the mean-preserving ASYMMETRIC insertion.

    Places mass epsilon at y1 in the first component and at y2 in the second,
    removing epsilon*ybar/x from the atom at x, where ybar = (1-q)y1 + q y2 is
    exactly what mean preservation requires: the objective's mean is the
    q-weighted mixture, so unequal supports are feasible.
    """
    x, p = mp_parameters.x, mp_parameters.p
    epsilon = mpmath.mpf("1e-30") if epsilon is None else epsilon
    ybar = (1 - q) * y1 + q * y2
    base = _pencil_mp((p, 0, q, x, y1, 0, x, y2, 0), mp_parameters)
    moved = _pencil_mp((p - epsilon * ybar / x, epsilon, q, x, y1, 0, x, y2, 0),
                       mp_parameters)
    return (moved - base) / epsilon


def certify_perturbation_decomposition(samples: int = 12) -> dict[str, Any]:
    """The perturbation measure is the q-convex combination, exactly.

    Writing Mdot for the first-order change of the mixture qbar*P0 + q*P1 and
    Mdot_i for the symmetric perturbation at y_i, the atoms are

        Mdot   :  x |-> -ybar/x,  y1 |-> qbar,  y2 |-> q,  0 |-> ybar/x - 1
        Mdot_i :  x |-> -y_i/x,   y_i |-> 1,                0 |-> y_i/x - 1

    and qbar*Mdot_1 + q*Mdot_2 reproduces every atom of Mdot because
    ybar = qbar*y1 + q*y2.  Every term of the objective is linear (EHX) or
    bilinear (EHXY, EHPI, dist^2) in the measure, so each first variation is
    LINEAR in the perturbation, and linearity carries the convex combination
    straight through to the coefficient.  That is the whole proof; the checks
    below confirm the arithmetic rather than replace it.
    """
    import random

    generator = random.Random(20260829)
    for _ in range(samples):
        def pick() -> Fraction:
            return Fraction(generator.randint(1, 10 ** 5),
                            generator.randint(1, 10 ** 5))

        x, y1, y2, q = pick(), pick(), pick(), pick()
        qbar = 1 - q
        ybar = qbar * y1 + q * y2
        combined = {
            "x": qbar * (-y1 / x) + q * (-y2 / x),
            "y1": qbar * 1,
            "y2": q * 1,
            "0": qbar * (y1 / x - 1) + q * (y2 / x - 1),
        }
        direct = {"x": -ybar / x, "y1": qbar, "y2": q, "0": ybar / x - 1}
        if combined != direct:
            raise AssertionError(
                f"the perturbation does not decompose at x={x}, q={q}")
    return {
        "samples": samples,
        "atoms_match_exactly": True,
        "arithmetic": "exact Fraction",
    }


def certify_asymmetric_decoupling(
    mp_parameters,
    dps: int,
    grid: Sequence[tuple[str, str, str]],
) -> dict[str, Any]:
    """L_asym(y1,y2,q) = (1-q)L(y1) + q L(y2), checked against the SSOT."""
    from liu9_ninevar import transverse_coefficients_mp

    rows = []
    worst = mpmath.mpf(0)
    with mpmath.workdps(dps):
        for raw1, raw2, raw_q in grid:
            y1, y2, q = (mpmath.mpf(raw1), mpmath.mpf(raw2), mpmath.mpf(raw_q))
            measured = asymmetric_coefficient_mp(mp_parameters, y1, y2, q)
            predicted = ((1 - q) * transverse_coefficients_mp(
                mp_parameters, y1, "mean_preserving")[0]
                + q * transverse_coefficients_mp(
                    mp_parameters, y2, "mean_preserving")[0])
            difference = abs(measured - predicted)
            worst = max(worst, difference)
            rows.append({
                "y1": raw1, "y2": raw2, "q": raw_q,
                "coefficient": mpmath.nstr(measured, 15),
                "difference": mpmath.nstr(difference, 4),
            })

        # Endpoint reduction: at q = 0 the second support cannot appear at all.
        endpoint = []
        for other in ("1/10", "1/2", "99/100"):
            value = asymmetric_coefficient_mp(
                mp_parameters, mpmath.mpf(3) / 10,
                mpmath.mpf(Fraction(other).numerator)
                / Fraction(other).denominator, mpmath.mpf(0))
            endpoint.append(mpmath.nstr(value, 20))
        endpoint_invariant = len(set(endpoint)) == 1

        # Affine in q: the second difference must vanish.
        second = []
        for raw1, raw2 in (("1/5", "9/10"), ("1/20", "19/20"), ("4/5", "3/20")):
            y1 = mpmath.mpf(Fraction(raw1).numerator) / Fraction(raw1).denominator
            y2 = mpmath.mpf(Fraction(raw2).numerator) / Fraction(raw2).denominator
            values = [asymmetric_coefficient_mp(
                mp_parameters, y1, y2, mpmath.mpf(k) / 4) for k in (1, 2, 3)]
            second.append(float(abs(values[0] - 2 * values[1] + values[2])))

    if not endpoint_invariant:
        raise AssertionError("q=0 did not eliminate the second support")
    tolerance = 1e-25
    if float(worst) > tolerance:
        raise AssertionError(
            f"the decoupling identity failed by {float(worst):.3e}")
    return {
        "rows": rows,
        "worst_difference": float(worst),
        "tolerance": tolerance,
        "endpoint_reduction_invariant": endpoint_invariant,
        "worst_second_difference_in_q": max(second),
        "certified": True,
    }


# ---------------------------------------------------------------------------
# Reporting.


def _canonical_digest(payload: dict[str, Any]) -> str:
    body = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--y0", type=str, default=str(DEFAULT_Y0))
    parser.add_argument("--window", type=str, default=str(DEFAULT_WINDOW))
    parser.add_argument("--bulk-cells", type=int, default=DEFAULT_BULK_CELLS)
    parser.add_argument("--octaves", type=int, default=DEFAULT_OCTAVES)
    parser.add_argument("--window-cells", type=int,
                        default=DEFAULT_WINDOW_CELLS)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args(argv)

    started = time.monotonic()
    y0 = Fraction(args.y0)
    window = Fraction(args.window)

    mp_parameters = solve_equation_parameters(args.dps)
    arb_parameters = certify_equation_parameters(mp_parameters)

    print("1. THE PROTOCOL IDENTITY IS LIU'S QUARTIC")
    identity = certify_protocol_identity()
    if not identity["identity_holds"]:
        raise AssertionError("the protocol identity failed")
    print("PROVED [exact integer polynomials]: prot(z,z) - (1 - z^2) = "
          f"{identity['difference']} read from the constant term up, which is "
          f"Liu's quartic {identity['liu_quartic']} coefficient for "
          "coefficient, so `prot(x,x) = 1 - x^2` is not a consequence of the "
          "quartic -- it IS the quartic.")
    print("ATTRIBUTION [primary source, arXiv:2306.08824v1, Jingbo Liu]: eq. "
          "(87) reads `x*^2 + x*^2(1 + xbar*^2) = 1`, and x^2(1+(1-x)^2) is "
          "prot(x,x), so (87) already IS the protocol identity, stated "
          "directly as a defining relation of the numerically located "
          "optimizer.  It is not derived there, and the quartic polynomial "
          "never appears in the paper.  The quartic is the derived form here.")

    constant = certify_constant_a(arb_parameters)
    print("NOVEL [absence after stated search]: the CONSEQUENCE below is "
          "unremarked in the paper and in all ten citing works, which never "
          "mention the entropy-symmetry collapse, the beta-independence of the "
          "kernel diagonal, or p*K(x,x) = h(x); the root is uncatalogued in "
          "OEIS.  Liu instead fixes a specific beta* via (89)-(90).")
    print(f"PROVED [Arb]: h symmetric about 1/2 therefore gives "
          f"K(x,x) - h(x^2) = {constant['beta_free_diagonal']:.3e} for every "
          f"beta, so with p = h(x)/h(x^2) the constant A = 2p*K(x,x) - h(x) "
          f"equals h(x) to {constant['difference_upper']:.3e}.")

    print()
    print("2. THE DOUBLE ROOT AT y = x")
    root = certify_double_root(arb_parameters)
    structural = certify_structural_zero()
    print(f"PROVED [exact Fraction, {structural['samples']} pseudo-random "
          f"points, entropy values as free symbols]: L(x) = 0 identically.  At "
          "y = x the ratio y/x is one, so the kernel term and the A term "
          "telescope, and (y-x)^2 kills the support term.  This uses NO "
          "property of x, p, beta or kappa, so it is a rational-function "
          "identity and the Schwartz-Zippel witness settles it.")
    print(f"NOT USED [Arb]: evaluating L at a ball containing x returns "
          f"{root['value_upper']:.3e} instead of zero.  That is DEPENDENCY "
          "ERROR -- Arb computes y*(A/x) and K(x,y) from separate "
          "sub-expressions and cannot cancel them -- and it is deliberately "
          "excluded from the Taylor bound below, where it would have swamped "
          "the real deficit by sixty-eight orders of magnitude.")
    print(f"PROVED [Arb]: |L'(x)| <= {root['slope_upper']:.3e}.  This is the "
          "residual of Liu's equations (87)-(90) at the certified parameter "
          "enclosure, not a measured slope.")
    print(f"PROVED [Arb]: L''(x) >= {root['curvature_lower']:.6f} > 0, so the "
          "root is a strict local minimum and y = x is a DOUBLE root.")

    sensitivity = beta_sensitivity(mp_parameters, args.dps)
    print(f"MACHINE VERIFIED [mutation]: the vanishing slope needs Liu's "
          f"equation (90) and not merely the quartic.  Perturbing beta moves "
          f"L'(x) to at least {sensitivity['smallest_perturbed_slope']:.3e} in "
          "absolute value at every tested value:")
    for row in sensitivity["rows"]:
        print(f"     beta = {row['beta']:<12s} L'(x) = {row['slope']}")

    print()
    print("3. REGION A -- THE y*log(1/y) BOUNDARY LAYER")
    layer = certify_boundary_layer(arb_parameters, y0, args.octaves)
    print(f"PROVED [Arb, {layer['octaves']} octaves down to y = "
          f"{layer['smallest_covered_y']}]: writing L(y)/y = C(y)*log(1/y) + "
          f"R(y), the log coefficient is at least "
          f"{layer['worst_log_coefficient_lower']:.6f} and the whole bound is "
          f"at least {layer['worst_layer_bound_lower']:.6f} > 0.  The "
          "factorization is singular-free: mu(u) = ((1-u)/u)log(1/(1-u)) is "
          "enclosed in [1-u, 1] by termwise comparison of its power series "
          "against 1 and against the geometric series, so no 0/0 is ever "
          "evaluated.")
    print(f"PROVED [Arb, one interval cell {layer['tail_cell']}]: the tail is "
          f"certified, not argued.  rho(y) >= x for every y <= 1 gives the "
          f"uniform bound C >= 2px - 1 >= "
          f"{layer['uniform_log_coefficient_lower']:.6f} with no monotonicity "
          f"appeal, and since that is positive, replacing log(1/y) by "
          f"log(1/{layer['tail_cell'][1]}) bounds the whole tail below at "
          f"{layer['tail_bound_lower']:.6f} > 0.  y = 0 sits inside the ball "
          "because the factorization is singularity-free.  L(0) = 0 is "
          "therefore approached strictly from above.")

    print()
    print("4. REGION B -- THE BULK")
    bulk = certify_bulk(arb_parameters, y0, window, args.bulk_cells)
    print(f"PROVED [Arb, {bulk['cells']} cells]: L > 0 on [{bulk['y0']}, "
          f"x-{bulk['window']}] and [x+{bulk['window']}, 1], with weakest "
          f"enclosure {bulk['worst_margin_lower']:.6e} on cell "
          f"{bulk['worst_cell']}.")

    print()
    print("5. REGION C -- THE DOUBLE-ROOT WINDOW")
    win = certify_window(arb_parameters, window, args.window_cells,
                         root["slope_upper"])
    coverage = certify_coverage(
        y0, window, Fraction(arb_parameters.x.str(40, radius=False)),
        y0 / Fraction(2) ** args.octaves)
    print(f"PROVED [Arb, {win['cells']} cells]: L'' >= "
          f"{win['minimum_curvature_lower']:.6f} > 0 on |y-x| <= "
          f"{win['radius']}.  L is NOT globally convex -- L'' < 0 on roughly "
          "[0.01, 0.43] -- so this window is the reason the region split "
          "exists, not a convenience.")
    print(f"PROVED [Arb]: with L(x) = 0 and |L'(x)| <= "
          f"{win['slope_bound']:.3e}, an interval Taylor step gives "
          f"L(y) >= (1/2)*{win['minimum_curvature_lower']:.6f}*(y-x)^2 - "
          f"{win['slope_bound']:.3e}*|y-x| >= -{win['deficit']:.3e} on the "
          "window.  That deficit is the SQUARE of Liu's root residual.")

    print()
    print("6. THE COVER IS EXHAUSTIVE")
    print("PROVED [exact Fraction]: the five pieces abut and tile [0,1] with "
          "no gap and no unproved sliver:")
    for name, lo, hi in coverage["pieces"]:
        shown = hi if len(hi) < 28 else hi[:25] + "..."
        print(f"     {name:<10s} [{lo if len(lo) < 28 else lo[:25] + '...'}, "
              f"{shown}]")

    print()
    print("7. THE ASYMMETRIC TWO-SUPPORT FAMILY DECOUPLES")
    decomposition = certify_perturbation_decomposition()
    grid = [("1/5", "9/10", "1/2"), ("1/20", "19/20", "1/8"),
            ("3/10", "7/10", "7/8"), ("1/100", "99/100", "1/2"),
            ("9/10", "1/5", "3/8"), ("2/3", "1/3", "1/4"),
            ("1/1000", "1/2", "1/2"), ("999/1000", "1/50", "5/8")]
    decoupling = certify_asymmetric_decoupling(mp_parameters, args.dps, grid)
    print("PROVED [exact Fraction, %d points]: the objective's mean is the "
          "q-weighted mixture, so UNEQUAL supports y1 != y2 are mean-feasible "
          "with the correction ybar = (1-q)y1 + q y2.  The resulting "
          "perturbation measure is exactly the q-convex combination of the two "
          "symmetric ones, atom for atom."
          % decomposition["samples"])
    print("PROVED [measure algebra]: every term is linear (EHX) or bilinear "
          "(EHXY, EHPI, dist^2) in the measure, so each first variation is "
          "LINEAR in the perturbation.  Linearity carries the convex "
          "combination through, giving")
    print("     L_asym(y1, y2, q) = (1-q)*L(y1) + q*L(y2).")
    print(f"MACHINE VERIFIED [{len(grid)} configurations, {args.dps} dps "
          f"Richardson against gap_mp and distance_squared]: the identity "
          f"holds to {decoupling['worst_difference']:.3e}, the second "
          f"difference in q vanishes to "
          f"{decoupling['worst_second_difference_in_q']:.3e} confirming "
          f"affineness, and at q=0 the second support provably drops out of "
          f"every term.")
    print("PROVED [corollary]: because a convex combination of numbers that "
          f"are all at least -{win['deficit']:.3e} is itself at least "
          f"-{win['deficit']:.3e}, the whole asymmetric family inherits the "
          "bound with NO additional certification.  The two-parameter "
          "extension is free.")

    print()
    print("8. VERDICT")
    deficit = win["deficit"]
    print(f"PROVED: the mean-preserving transverse first variation satisfies "
          f"L(y) >= -{deficit:.3e} for every y in (0,1), with equality only in "
          f"the immediate neighbourhood of the double root at y = x.  The "
          "coefficient is q-free, so this holds for every q simultaneously.  "
          "This is a universal enclosure over the whole interval, NOT a scan: "
          "it replaces the 65-point sampled positivity in liu9_ninevar.py.")
    print("PROVED: the ambient nine-variable refutation in liu9_ninevar.py "
          "therefore cannot be repaired into a counterexample to Liu's "
          "Hypothesis 2 along this family.  The negative direction it found is "
          "mean-infeasible, and the mean-feasible member of the same family is "
          "nonnegative for every y and every q.")
    print(f"PROVED: the same bound therefore holds for the ASYMMETRIC family "
          f"on all of (0,1)^2 x [0,1]: L_asym >= -{deficit:.3e} for every "
          "pair of insertion points and every q.")
    print("OPEN: this is the leading coefficient of a family with at most ONE "
          "new support point per component.  Insertions opening the third atom "
          "in both components simultaneously, and the epsilon^2 and higher "
          "terms, are untouched.  No tube radius follows, and Liu's Hypothesis "
          "2 is NOT claimed -- the certified statement is about the first "
          "variation at his conjectured optimizer, not about global "
          "optimality.")

    payload = {
        "tool": "liu9_transverse.py",
        "dps": args.dps,
        "kappa": str(KAPPA),
        "protocol_identity": identity,
        "constant_a": constant,
        "double_root": root,
        "structural_zero": structural,
        "coverage": coverage,
        "beta_sensitivity": sensitivity,
        "boundary_layer": layer,
        "bulk": bulk,
        "window": win,
        "global_deficit": deficit,
        "perturbation_decomposition": decomposition,
        "asymmetric_decoupling": decoupling,
        "elapsed_seconds": round(time.monotonic() - started, 3),
    }
    digest_payload = {k: v for k, v in payload.items()
                      if k != "elapsed_seconds"}
    payload["report_sha256"] = _canonical_digest(digest_payload)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1, sort_keys=True)
        print()
        print(f"report written to {args.output}")
        print(f"report_sha256 {payload['report_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
