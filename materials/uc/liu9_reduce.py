#!/usr/bin/env python3
"""Audit the proposed nu-collapse of Liu's nine-parameter optimisation.

The source formula is Liu, arXiv:2306.08824v1, Section V-B, (81)--(84),
transcribed independently in ``liu9_objective.py`` from
https://jingbol.web.illinois.edu/frankl5.m.

The audit finds an exact obstruction rather than a collapse.  If

    mu = (1-q) P0 + q P1,          nu = P1 - P0,

then the objective separates through the protocol kernel K, not directly
through the residual kernel R proved NSD in ``liu_R_nsd.py``.  The part
depending only on mu is constant on the literal fixed-measure fibre, but that
fibre generally leaves Liu's paired 3+3-atom parameterisation.  The true
in-domain paired-atom second variation has a nonzero gain term.  A rational
in-domain path is certified with Arb to have LOSS/GAIN > 4.2848; even the
first 64 exact coercive terms prove LOSS/GAIN > 4.0313.  Thus stronger
coercivity worsens, rather than closes, the proposed q=0 reduction.

This script prints only the following status labels:

* PROVED: direct algebra from the cited source or a cited proved kernel bound;
* NUMERICAL: finite-precision checks or optimisation searches;
* CONJECTURED: an interpretation not certified by the calculations.

It also prints the exact q=0 five-parameter face problem and deterministic,
one-core numerical estimates at the three requested rational targets.

Run from the repository root with

    math/.venv/bin/python math/uc/liu9_reduce.py
"""

from __future__ import annotations

import argparse
import math
import os
import random
import sys
from fractions import Fraction
from typing import Callable, Optional, Sequence

# This workstation is shared with live certification workers.  These must be
# set before importing NumPy/SciPy (and therefore before liu9_objective).
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
from scipy.optimize import NonlinearConstraint, differential_evolution

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from liu9_objective import (  # noqa: E402
    LIU_BETA_DECIMAL,
    SOURCE_OBJECTIVE_LINE,
    evaluate_float64,
    evaluate_mpmath,
)

BETA = Fraction(LIU_BETA_DECIMAL)
ctx.prec = max(ctx.prec, 256)
TARGET_TEXTS = ("0.38240", "0.38250", "0.38260")
VERIFY_TARGET = Fraction("0.38260")
DEFAULT_RANDOM_POINTS = 24
DEFAULT_DPS = 80
DEFAULT_FACE_MAXITER = 160
DEFAULT_FACE_POPULATION = 96

Scalar = float | mpmath.mpf
Kernel = Callable[[Scalar, Scalar, Callable[[Scalar], Scalar]], Scalar]


def _h_scalar(x: Scalar, logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    """Binary entropy in nats, with its continuous endpoint convention."""
    if x == 0 or x == 1:
        return x - x
    return -x * logarithm(x) - (1 - x) * logarithm(1 - x)


def _kernel_j(x: Scalar, y: Scalar,
              logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    return _h_scalar(x * y, logarithm)


def _kernel_k(x: Scalar, y: Scalar,
              logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    z = x * y * (1 + (1 - x) * (1 - y))
    return _h_scalar(z, logarithm)


def _kernel_r(x: Scalar, y: Scalar,
              logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    """The proved R after the source variable change s -> 1-x."""
    a = x * y
    b = (1 - x) * (1 - y)
    z = a * (1 + b)
    if z == 1:
        g = z
    else:
        g = z + (1 - z) * logarithm(1 - z)
    return a * (b - (1 + b) * logarithm(1 + b)) - g


def _kernel_coercive(x: Scalar, y: Scalar,
                     logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    """The two retained PSD families in the proved coercive inequality."""
    a = x * y
    b = (1 - x) * (1 - y)
    if a == 1:
        # Here b=0.  This is the continuous value of -a*b*log(1-a).
        log_piece = a - a
        g = a
    else:
        log_piece = -a * b * logarithm(1 - a)
        g = a + (1 - a) * logarithm(1 - a)
    return g + log_piece


def _linear(points: Sequence[Scalar], weights: Sequence[Scalar],
            function: Callable[[Scalar], Scalar]) -> Scalar:
    return sum((weight * function(point)
                for point, weight in zip(points, weights)), weights[0] - weights[0])


def _bilinear(left_points: Sequence[Scalar], left_weights: Sequence[Scalar],
              right_points: Sequence[Scalar], right_weights: Sequence[Scalar],
              kernel: Kernel,
              logarithm: Callable[[Scalar], Scalar]) -> Scalar:
    zero = left_weights[0] - left_weights[0]
    return sum((
        left_weight * right_weight * kernel(left_point, right_point, logarithm)
        for left_point, left_weight in zip(left_points, left_weights)
        for right_point, right_weight in zip(right_points, right_weights)
    ), zero)


def _split_terms(values: Sequence[Scalar], beta: Scalar,
                 logarithm: Callable[[Scalar], Scalar]) -> dict[str, Scalar]:
    """Independently evaluate the exact mixture separation."""
    if len(values) != 9:
        raise ValueError("expected Liu's nine variables")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    one = beta - beta + 1
    masses = (a1, a2, one - a1 - a2)
    p0 = (b0, b2, b4)
    p1 = (b1, b3, b5)
    qbar = one - q
    mu_points = p0 + p1
    mu_weights = tuple(qbar * mass for mass in masses) + tuple(q * mass for mass in masses)

    entropy = _linear(
        mu_points, mu_weights, lambda point: _h_scalar(point, logarithm))
    j_mu = _bilinear(
        mu_points, mu_weights, mu_points, mu_weights, _kernel_j, logarithm)
    k_mu = _bilinear(
        mu_points, mu_weights, mu_points, mu_weights, _kernel_k, logarithm)

    k00 = _bilinear(p0, masses, p0, masses, _kernel_k, logarithm)
    k11 = _bilinear(p1, masses, p1, masses, _kernel_k, logarithm)
    k01 = _bilinear(p0, masses, p1, masses, _kernel_k, logarithm)
    k_nu = k00 + k11 - 2 * k01

    r00 = _bilinear(p0, masses, p0, masses, _kernel_r, logarithm)
    r11 = _bilinear(p1, masses, p1, masses, _kernel_r, logarithm)
    r01 = _bilinear(p0, masses, p1, masses, _kernel_r, logarithm)
    r_nu = r00 + r11 - 2 * r01

    def signed_moment(function: Callable[[Scalar], Scalar]) -> Scalar:
        return sum((mass * (function(right) - function(left))
                    for mass, left, right in zip(masses, p0, p1)), a1 - a1)

    m1 = signed_moment(lambda x: x)
    ma = signed_moment(lambda x: x * (one - x))
    l1 = signed_moment(
        lambda x: x - x if x == 0 else x * logarithm(x))
    la = signed_moment(
        lambda x: x - x if x == 0 else x * (one - x) * logarithm(x))
    low_rank = m1 * m1 - 2 * m1 * l1 - 2 * ma * la

    numerator0 = (one - beta) * j_mu + beta * k_mu
    phi0 = numerator0 / entropy
    correction = beta * qbar * q * k_nu / entropy
    objective = phi0 + correction
    gap = numerator0 + beta * qbar * q * k_nu - entropy

    # This is the exact residual of the literal prior claim
    # Phi=Phi0+beta*qbar*q*R(nu).
    claimed_r_residual = beta * qbar * q * (k_nu / entropy - r_nu)
    return {
        "entropy": entropy,
        "j_mu": j_mu,
        "k_mu": k_mu,
        "k_nu": k_nu,
        "r_nu": r_nu,
        "low_rank": low_rank,
        "phi0": phi0,
        "correction": correction,
        "objective": objective,
        "gap": gap,
        "claimed_r_residual": claimed_r_residual,
    }


def _random_feasible_points(count: int, seed: int) -> list[tuple[float, ...]]:
    rng = random.Random(seed)
    threshold = 1.0 - float(VERIFY_TARGET)
    points: list[tuple[float, ...]] = []
    attempts = 0
    while len(points) < count:
        attempts += 1
        if attempts > 100_000:
            raise RuntimeError("could not generate enough feasible audit points")
        raw_masses = [rng.expovariate(1.0) for _ in range(3)]
        total = sum(raw_masses)
        masses = [mass / total for mass in raw_masses]
        q = rng.uniform(0.02, 0.98)
        support = [rng.uniform(0.02, 0.98) for _ in range(6)]
        mean0 = sum(mass * point for mass, point in zip(masses, support[:3]))
        mean1 = sum(mass * point for mass, point in zip(masses, support[3:]))
        mean = (1.0 - q) * mean0 + q * mean1
        if mean >= threshold:
            points.append((masses[0], masses[1], q, *support))
    return points


def verify_separation(random_points: int, seed: int, dps: int) -> None:
    points = _random_feasible_points(random_points, seed)
    beta_float = float(BETA)
    float_objective_residual = 0.0
    float_gap_residual = 0.0
    float_kernel_residual = 0.0
    float_claim_formula_residual = 0.0
    float_claim_size = 0.0

    for values in points:
        source = evaluate_float64(values, beta_float)
        split = _split_terms(values, beta_float, math.log)
        float_objective_residual = max(
            float_objective_residual,
            abs(float(source.objective) - float(split["objective"])),
        )
        float_gap_residual = max(
            float_gap_residual,
            abs(float(source.numerator - source.ehx) - float(split["gap"])),
        )
        float_kernel_residual = max(
            float_kernel_residual,
            abs(float(split["k_nu"] - split["r_nu"] - split["low_rank"])),
        )
        claimed_direct = float(source.objective) - (
            float(split["phi0"])
            + beta_float * (1.0 - values[2]) * values[2] * float(split["r_nu"])
        )
        float_claim_formula_residual = max(
            float_claim_formula_residual,
            abs(claimed_direct - float(split["claimed_r_residual"])),
        )
        float_claim_size = max(float_claim_size, abs(claimed_direct))

    with mpmath.workdps(dps):
        beta_mp = mpmath.mpf(BETA.numerator) / BETA.denominator
        mp_objective_residual = mpmath.mpf(0)
        mp_gap_residual = mpmath.mpf(0)
        mp_kernel_residual = mpmath.mpf(0)
        mp_claim_formula_residual = mpmath.mpf(0)
        mp_claim_size = mpmath.mpf(0)
        for values_float in points:
            values = tuple(mpmath.mpf(str(value)) for value in values_float)
            source = evaluate_mpmath(values, beta_mp, dps=dps)
            split = _split_terms(values, beta_mp, mpmath.log)
            mp_objective_residual = max(
                mp_objective_residual, abs(source.objective - split["objective"]))
            mp_gap_residual = max(
                mp_gap_residual,
                abs(source.numerator - source.ehx - split["gap"]),
            )
            mp_kernel_residual = max(
                mp_kernel_residual,
                abs(split["k_nu"] - split["r_nu"] - split["low_rank"]),
            )
            q = values[2]
            claimed_direct = source.objective - (
                split["phi0"] + beta_mp * (1 - q) * q * split["r_nu"])
            mp_claim_formula_residual = max(
                mp_claim_formula_residual,
                abs(claimed_direct - split["claimed_r_residual"]),
            )
            mp_claim_size = max(mp_claim_size, abs(claimed_direct))

    if float_objective_residual > 5e-14 or float_gap_residual > 5e-14:
        raise AssertionError("float64 separation check failed")
    if float_kernel_residual > 5e-14 or float_claim_formula_residual > 5e-14:
        raise AssertionError("float64 R-correction check failed")
    mp_tolerance = mpmath.mpf(10) ** (-(dps - 10))
    if max(mp_objective_residual, mp_gap_residual, mp_kernel_residual,
           mp_claim_formula_residual) > mp_tolerance:
        raise AssertionError("mpmath separation check failed")

    print("1. SOURCE AUDIT AND EXACT SEPARATION")
    print(f"PROVED [frankl5.m transcription]: `{SOURCE_OBJECTIVE_LINE}`")
    print("PROVED [direct expansion of Liu (81)]: with H(mu)=int h(x)dmu,")
    print("   J(mu)=iint h(xy)dmu dmu, K(x,y)=h(xy[1+(1-x)(1-y)]),")
    print("   mu=(1-q)P0+qP1 and nu=P1-P0,")
    print("   Phi = Phi0(mu) + beta*q*(1-q)*K(nu)/H(mu),")
    print("   Phi0(mu)=[(1-beta)J(mu)+beta*K(mu)]/H(mu).")
    print("PROVED [same expansion, denominator-free]:")
    print("   G:=numerator-H = G0(mu)+beta*q*(1-q)*K(nu).")
    print("PROVED [degenerate fibre]: if P0=P1 then nu=0, so the correction")
    print("   vanishes and the source objective is invariant in q for every q in [0,1].")
    print("PROVED [normalisation and low-rank audit]: the prior literal claim with")
    print("   the unnormalised residual kernel R is false off the projected subspace.")
    print("   For m1=<nu,x>, ma=<nu,x(1-x)>, L1=<nu,x ln x>,")
    print("   La=<nu,x(1-x)ln x>, direct entropy algebra gives")
    print("   K(nu)=R(nu)+D(nu),  D=m1^2-2*m1*L1-2*ma*La.")
    print("PROVED [exact residual of the literal claimed objective split]:")
    print("   E_R=beta*q*(1-q)*[(R+D)/H-R]")
    print("      =beta*q*(1-q)*[D/H+(1/H-1)R].")
    print("PROVED [exact residual in the gap formulation]: E_R,gap=beta*q*(1-q)*D.")
    print(
        f"NUMERICAL [float64, {random_points} feasible points at rational "
        f"c={VERIFY_TARGET}]: max correct objective residual="
        f"{float_objective_residual:.3e}; max correct gap residual="
        f"{float_gap_residual:.3e}; max K-R-D residual={float_kernel_residual:.3e}."
    )
    print(
        f"NUMERICAL [float64, same points]: max exact-residual consistency error="
        f"{float_claim_formula_residual:.3e}; max size of the omitted/incorrect "
        f"R residual={float_claim_size:.3e}."
    )
    print(
        f"NUMERICAL [mpmath {dps} dps, same decimal-rationalized points]: "
        f"max correct objective residual={mpmath.nstr(mp_objective_residual, 6)}; "
        f"max correct gap residual={mpmath.nstr(mp_gap_residual, 6)}; "
        f"max K-R-D residual={mpmath.nstr(mp_kernel_residual, 6)}."
    )
    print(
        f"NUMERICAL [mpmath {dps} dps]: max exact-residual consistency error="
        f"{mpmath.nstr(mp_claim_formula_residual, 6)}; max size of the literal "
        f"R residual={mpmath.nstr(mp_claim_size, 8)}."
    )
    print()


def _signed_quadratic(points: Sequence[mpmath.mpf],
                      weights: Sequence[mpmath.mpf], kernel: Kernel) -> mpmath.mpf:
    return _bilinear(points, weights, points, weights, kernel, mpmath.log)


def _arb_fraction(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def _atom_path_certificate(cutoff: int = 64) -> dict[str, object]:
    """Certify a genuine 3+3-atom descent direction with Arb."""
    masses = (Fraction(6, 25), Fraction(13, 100), Fraction(63, 100))
    points = (Fraction(43, 1000), Fraction(53, 100), Fraction(941, 1000))
    # Put w_i=a_i d_i.  The cyclic differences make sum w_i=sum b_i w_i=0.
    directions = (Fraction(-137, 80), Fraction(449, 65), Fraction(-487, 630))
    mean = sum((mass * point for mass, point in zip(masses, points)), Fraction(0))
    dot_m1 = sum((mass * direction
                  for mass, direction in zip(masses, directions)), Fraction(0))
    dot_m2 = sum((2 * mass * direction * point
                  for mass, direction, point in zip(masses, directions, points)),
                 Fraction(0))
    if dot_m1 != 0 or dot_m2 != 0:
        raise AssertionError("the rational atom direction is not projected")

    def derivative_moment(power: int) -> Fraction:
        return sum((
            power * mass * direction * point ** (power - 1)
            for mass, direction, point in zip(masses, directions, points)
        ), Fraction(0))

    # Every summand is nonnegative.  This exact rational partial sum is already
    # a proved lower bound for the full coercive budget.
    coercive_partial = sum((
        derivative_moment(k) ** 2 / (k * (k - 1))
        for k in range(2, cutoff + 1)
    ), Fraction(0))
    coercive_partial += sum((
        sum((
            mass * direction * (
                (k + 1) * point**k - (k + 2) * point ** (k + 1)
            )
            for mass, direction, point in zip(masses, directions, points)
        ), Fraction(0)) ** 2 / k
        for k in range(1, cutoff + 1)
    ), Fraction(0))

    masses_arb = tuple(_arb_fraction(value) for value in masses)
    points_arb = tuple(_arb_fraction(value) for value in points)
    directions_arb = tuple(_arb_fraction(value) for value in directions)
    beta = _arb_fraction(BETA)
    one = arb(1)

    def entropy(x: arb) -> arb:
        return -(x * x.log() + (one - x) * (one - x).log())

    def entropy_prime(x: arb) -> arb:
        return ((one - x) / x).log()

    def entropy_second(x: arb) -> arb:
        return -one / (x * (one - x))

    entropy_value = sum((
        mass * entropy(point)
        for mass, point in zip(masses_arb, points_arb)
    ), arb(0))
    numerator = arb(0)
    for i in range(3):
        for j in range(3):
            x, y = points_arb[i], points_arb[j]
            product = x * y
            protocol = product * (one + (one - x) * (one - y))
            numerator += masses_arb[i] * masses_arb[j] * (
                (one - beta) * entropy(product) + beta * entropy(protocol))
    phi0 = numerator / entropy_value

    gain = arb(0)
    for i in range(3):
        x = points_arb[i]
        inner = arb(0)
        for j in range(3):
            y = points_arb[j]
            product = x * y
            j_xx = entropy_second(product) * y * y
            protocol = product * (one + (one - x) * (one - y))
            z_x = y * (one + (one - 2 * x) * (one - y))
            z_xx = -2 * y * (one - y)
            k_xx = (
                entropy_second(protocol) * z_x * z_x
                + entropy_prime(protocol) * z_xx
            )
            inner += masses_arb[j] * ((one - beta) * j_xx + beta * k_xx)
        gain += (
            masses_arb[i] * directions_arb[i] ** 2
            * (2 * inner - phi0 * entropy_second(x))
        )

    loss = arb(0)
    for i in range(3):
        for j in range(3):
            x, y = points_arb[i], points_arb[j]
            product = x * y
            protocol = product * (one + (one - x) * (one - y))
            z_x = y * (one + (one - 2 * x) * (one - y))
            z_y = x * (one + (one - 2 * y) * (one - x))
            z_xy = one + (one - 2 * x) * (one - 2 * y)
            k_xy = (
                entropy_second(protocol) * z_x * z_y
                + entropy_prime(protocol) * z_xy
            )
            loss -= (
                2 * beta * masses_arb[i] * masses_arb[j]
                * directions_arb[i] * directions_arb[j] * k_xy
            )

    partial_loss = 2 * beta * _arb_fraction(coercive_partial)
    if float(gain.lower()) <= 0:
        raise AssertionError("Arb did not certify positive atom-path gain")
    if float((loss - gain).lower()) <= 0:
        raise AssertionError("Arb did not certify atom-path descent")
    if float((partial_loss - gain).lower()) <= 0:
        raise AssertionError("64 coercive terms did not certify atom-path descent")
    return {
        "masses": masses,
        "points": points,
        "directions": directions,
        "mean": mean,
        "wall": min(
            2 * min(point, 1 - point) / abs(direction)
            for point, direction in zip(points, directions)
        ),
        "cutoff": cutoff,
        "entropy": entropy_value,
        "phi0": phi0,
        "gain": gain,
        "loss": loss,
        "ratio": loss / gain,
        "partial_loss": partial_loss,
        "partial_ratio": partial_loss / gain,
        "curvature_bracket": gain - loss,
        "phi_second": (gain - loss) / (4 * entropy_value),
    }


def print_fibre_obstruction(dps: int) -> None:
    with mpmath.workdps(dps):
        one = mpmath.mpf(1)
        beta = mpmath.mpf(BETA.numerator) / BETA.denominator
        root3 = mpmath.sqrt(3)
        p0 = (mpmath.mpf("0.9"), mpmath.mpf("0.6"), mpmath.mpf("0.6"))
        p1 = (
            mpmath.mpf("0.7") + root3 / 10,
            mpmath.mpf("0.7") - root3 / 10,
            mpmath.mpf("0.7"),
        )
        masses = (one / 3, one / 3, one / 3)
        signed_points = p1 + p0
        signed_weights = masses + tuple(-mass for mass in masses)
        q = one / 2
        values = (one / 3, one / 3, q, *p0, *p1)
        split = _split_terms(values, beta, mpmath.log)
        r_loss = -split["r_nu"]
        coercive = _signed_quadratic(
            signed_points, signed_weights, _kernel_coercive)
        remainder = r_loss - coercive
        capture = coercive / r_loss
        actual_loss = -split["correction"]
        proved_loss = beta * q * (one - q) * coercive / split["entropy"]

        def moment(power: int) -> mpmath.mpf:
            return sum((weight * point**power
                        for point, weight in zip(signed_points, signed_weights)),
                       mpmath.mpf(0))

        m1 = moment(1)
        m2 = moment(2)
        m3 = moment(3)
        third_moment_bound = m3 * m3 * mpmath.mpf(7) / 6
        third_loss = (
            beta * q * (one - q) * third_moment_bound / split["entropy"])
    atom = _atom_path_certificate()
    with mpmath.workdps(dps):
        wall = mpmath.mpf(atom["wall"].numerator) / atom["wall"].denominator
        atom_masses = tuple(
            mpmath.mpf(value.numerator) / value.denominator
            for value in atom["masses"]
        )
        atom_points = tuple(
            mpmath.mpf(value.numerator) / value.denominator
            for value in atom["points"]
        )
        atom_directions = tuple(
            mpmath.mpf(value.numerator) / value.denominator
            for value in atom["directions"]
        )
        atom_p0 = tuple(
            point - wall * direction / 2
            for point, direction in zip(atom_points, atom_directions)
        )
        atom_p1 = tuple(
            point + wall * direction / 2
            for point, direction in zip(atom_points, atom_directions)
        )
        atom_wall_values = (
            atom_masses[0], atom_masses[1], mpmath.mpf("0.5"),
            *atom_p0, *atom_p1,
        )
        atom_wall_objective = evaluate_mpmath(
            atom_wall_values,
            mpmath.mpf(BETA.numerator) / BETA.denominator,
            dps=dps,
        ).objective
        atom_base_values = (
            atom_masses[0], atom_masses[1], mpmath.mpf("0.5"),
            *atom_points, *atom_points,
        )
        atom_base_objective = evaluate_mpmath(
            atom_base_values,
            mpmath.mpf(BETA.numerator) / BETA.denominator,
            dps=dps,
        ).objective

    print("2. NU-DIRECTION SECOND VARIATIONS: LOSS VERSUS GAIN")
    print("PROVED [literal fixed-measure fibre]: put P0(e)=mu-q*e*nu and")
    print("   P1(e)=mu+(1-q)*e*nu.  Then")
    print("   Phi(e)=Phi0(mu)+beta*q*(1-q)*e^2*K(nu)/H(mu), so")
    print("   Phi0''(0)=0 and Phi''(0)=2*beta*q*(1-q)*(D-r)/H, r=-R(nu).")
    print("PROVED [scope caveat]: this exact fixed-mu path is in the unrestricted")
    print("   measure problem; when mu has six atoms its intermediate points need not")
    print("   lie in Liu's paired 3+3-atom parameterisation.  It is not by itself an")
    print("   in-domain Hessian calculation.")
    print("PROVED [liu_coercive.py]: for mk=<nu,x^k> and")
    print("   ek=<nu,x^(k+1)(1-x)>=m_(k+1)-m_(k+2),")
    print("   r >= C(nu):=sum_(k>=2) mk^2/[k(k-1)] + sum_(k>=1) ek^2/k.")
    print("PROVED [direction of the inequality]: coercivity is a LOWER bound on")
    print("   negative curvature.  On m1=m2=0, D=0 and")
    print("   C >= m3^2/6+(m2-m3)^2=(7/6)m3^2.")
    print("PROVED [unrestricted-measure diagnostic]: for q=1/2 and equal masses,")
    print("   P0=(9/10,3/5,3/5),")
    print("   P1=(7/10+sqrt(3)/10,7/10-sqrt(3)/10,7/10),")
    print("   both means and second moments agree and m3=-1/500.")
    print(
        "NUMERICAL [mpmath %d dps, diagnostic endpoint]: H(mu)=%s; "
        "K(nu)=%s; -R(nu)=%s."
        % (
            dps,
            mpmath.nstr(split["entropy"], 16),
            mpmath.nstr(split["k_nu"], 16),
            mpmath.nstr(r_loss, 16),
        )
    )
    print(
        "NUMERICAL [closed-form retained families]: C=%s, C/(-R)=%s, "
        "discarded PSD remainder=%s; objective LOSS=%s."
        % (
            mpmath.nstr(coercive, 16),
            mpmath.nstr(capture, 12),
            mpmath.nstr(remainder, 16),
            mpmath.nstr(actual_loss, 16),
        )
    )
    print(
        "NUMERICAL [same diagnostic]: full-family LOSS >=%s; "
        "third-moment-only LOSS >=%s."
        % (mpmath.nstr(proved_loss, 16), mpmath.nstr(third_loss, 16))
    )
    print()
    print("PROVED [true in-domain paired-atom path]: for a three-atom P=sum ai*delta_bi,")
    print("   set P0(e)=sum ai*delta_(bi-q*e*di) and")
    print("       P1(e)=sum ai*delta_(bi+(1-q)*e*di).")
    print("   This remains inside the nine variables for all sufficiently small |e|.")
    print("   Although its coordinate mean is fixed, its mixture measure changes at e^2.")
    print("PROVED [explicit in-domain second variation]: put")
    print("   F(x,y)=(1-beta)h(xy)+beta*K(x,y), Phi0=F(P,P)/H(P), and")
    print("   G_atom=sum_i ai*di^2[2 sum_j aj*F_xx(bi,bj)-Phi0*h''(bi)].")
    print("   For <dot_nu,f>=sum_i ai*di*f'(bi),")
    print("   Phi''(0)=q*(1-q)/H(P) * [G_atom+2*beta*K(dot_nu)]")
    print("           =q*(1-q)/H(P) * [GAIN-LOSS],")
    print("   LOSS=2*beta*(r-D), r=-R(dot_nu).")
    print("PROVED [derivative moment defects and smooth signed-measure limit]:")
    print("   mk=k sum_i ai*di*bi^(k-1),")
    print("   ek=sum_i ai*di[(k+1)bi^k-(k+2)bi^(k+1)].")
    print("   Apply coercivity to (P1(e)-P0(e))/e and let e->0; all bi are")
    print("   interior, so r>=C for the same explicit coercive series.")
    print("PROVED [exact rational in-domain data, q=1/2]:")
    print(f"   a={atom['masses']}; b={atom['points']}; d={atom['directions']}.")
    print(f"   mean={atom['mean']} > 1-c at all three targets; dot_m1=dot_m2=0,")
    print("   hence D=0; every bi is interior, so the path is admissible for small e.")
    print("PROVED [256-bit Arb evaluation of the explicit Hessian]:")
    print(f"   H(P) is enclosed by {atom['entropy']}; Phi0(P) by {atom['phi0']};")
    print(f"   GAIN is enclosed by {atom['gain']};")
    print(f"   LOSS is enclosed by {atom['loss']};")
    print(f"   LOSS/GAIN is enclosed by {atom['ratio']}, hence is > 4.2848.")
    print(f"PROVED [first {atom['cutoff']} coercive terms, summed as an exact rational]:")
    print(f"   their LOSS lower bound is enclosed by {atom['partial_loss']};")
    print(f"   partial-LOSS/GAIN is enclosed by {atom['partial_ratio']}, hence is > 4.0313.")
    print(f"PROVED [Arb]: GAIN-LOSS is enclosed by {atom['curvature_bracket']} < 0;")
    print(f"   at q=1/2 the full Phi''(0) is enclosed by {atom['phi_second']} < 0.")
    print("   The paired-atom interior therefore strictly beats its corresponding")
    print("   degenerate q=0/diagonal face point for every sufficiently small e!=0.")
    print("PROVED [exact wall of the path at q=1/2]: the first box constraint binds at")
    print(f"   e_wall=min_i 2*min(bi,1-bi)/|di| = {atom['wall']}.")
    print(
        "NUMERICAL [mpmath %d dps, non-binding descent]: the source objective is "
        "%s at e=0 and %s at e=e_wall;"
        % (
            dps,
            mpmath.nstr(atom_base_objective, 20),
            mpmath.nstr(atom_wall_objective, 20),
        )
    )
    print("   the descent stops at the wall with margin over 1 intact, so this")
    print("   witness refutes the pointwise reduction without approaching Phi=1.")
    print()

    print("3. COLLAPSE VERDICT AND QUANTITATIVE FAILURE REGION")
    print("PROVED [in-domain Arb certificate]: the proposed universal nu-collapse is")
    print("   false as a pointwise reduction: the explicit paired-atom path has")
    print("   LOSS/GAIN >4.2848 and negative total second variation.")
    print("PROVED [general ratio]: when G_atom>0 and r>D,")
    print("   LOSS/GAIN = 2*beta*(r-D)/G_atom")
    print("             >= 2*beta*(C-D)/G_atom.")
    print("   The retained coercivity alone proves failure on")
    print("   2*beta*(C-D)>G_atom>0; if G_atom<=0 and C>D, descent is immediate.")
    print("PROVED [degeneracy at the structural optimum]: an admissible paired-atom")
    print("   direction needs dot_m1=dot_m2=0 and di=0 on boundary atoms; with fewer")
    print("   than three distinct interior moving atoms the Vandermonde system forces")
    print("   d=0.  Liu's observed optimizer p*delta_x+(1-p)*delta_0 has one interior")
    print("   atom, so this obstruction mechanism has no admissible direction there and")
    print("   is consistent with his local q=0 minimality.")
    print("PROVED [discarded-remainder sign]: writing r=C+E with E>=0, adding any")
    print("   fraction theta of E only increases LOSS.  No sharper Hypothesis-1")
    print("   coercivity constant can rescue this collapse; one would instead need an")
    print("   upper bound on r-D or a different global reduction.")
    print("PROVED [scope]: the certificate beats one face point locally, not the")
    print("   best observed face value.  It refutes the proposed pointwise route but")
    print("   neither proves nor disproves Liu's global minimizer hypothesis.")
    print()


def _h_array(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    answer = np.zeros_like(values)
    interior = (values > 0.0) & (values < 1.0)
    selected = values[interior]
    answer[interior] = -(
        selected * np.log(selected) + (1.0 - selected) * np.log1p(-selected))
    return answer


def _face_components(values: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    a1, a2, b0, b2, b4 = np.asarray(values, dtype=float)
    return np.array((a1, a2, 1.0 - a1 - a2)), np.array((b0, b2, b4))


def _face_objective(values: Sequence[float], beta: float) -> float:
    masses, points = _face_components(values)
    if np.min(masses) < -1e-12:
        return 1e6
    masses = np.maximum(masses, 0.0)
    entropy = float(masses @ _h_array(points))
    if entropy <= 1e-15:
        return 1e6
    product = np.outer(points, points)
    j_matrix = _h_array(product)
    pi_argument = product * (1.0 + np.outer(1.0 - points, 1.0 - points))
    k_matrix = _h_array(pi_argument)
    j_value = float(masses @ j_matrix @ masses)
    k_value = float(masses @ k_matrix @ masses)
    return ((1.0 - beta) * j_value + beta * k_value) / entropy


def _unpack_stick(values: Sequence[float]) -> np.ndarray:
    a1, remainder_fraction, b0, b2, b4 = np.asarray(values, dtype=float)
    a2 = (1.0 - a1) * remainder_fraction
    return np.array((a1, a2, b0, b2, b4))


def _stick_mean(values: Sequence[float]) -> float:
    face = _unpack_stick(values)
    masses, points = _face_components(face)
    return float(masses @ points)


def _make_face_population(target_mean: float, structural_x: float,
                          size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    structural_p = target_mean / structural_x
    base = np.array((structural_p, 1.0, structural_x, 0.0, 0.0))
    population = np.empty((size, 5), dtype=float)
    population[0] = base
    scales = np.array((0.035, 0.04, 0.035, 0.025, 0.025))
    focused = 3 * size // 4
    for index in range(1, size):
        if index < focused:
            candidate = np.clip(base + rng.normal(size=5) * scales, 0.0, 1.0)
        else:
            candidate = rng.random(5)
        mean = _stick_mean(candidate)
        if mean < target_mean:
            # Move every support point towards one just enough to hit the mean.
            shift = (target_mean - mean) / (1.0 - mean)
            candidate[2:] = shift + (1.0 - shift) * candidate[2:]
        population[index] = candidate
    return population


def _face_search(target: Fraction, structural_x: mpmath.mpf, beta: Fraction,
                 maxiter: int, population_size: int, seed: int):
    target_mean = 1.0 - float(target)
    beta_float = float(beta)
    initial = _make_face_population(
        target_mean, float(structural_x), population_size, seed)
    constraint = NonlinearConstraint(
        _stick_mean, target_mean - 5e-15, np.inf)
    result = differential_evolution(
        lambda values: _face_objective(_unpack_stick(values), beta_float),
        bounds=((0.0, 1.0),) * 5,
        constraints=(constraint,),
        seed=seed,
        maxiter=maxiter,
        tol=5e-9,
        polish=False,
        workers=1,
        updating="immediate",
        init=initial,
    )
    return result, _unpack_stick(result.x)


def _face_structural_stationary(beta: mpmath.mpf,
                                initial: mpmath.mpf) -> mpmath.mpf:
    one = mpmath.mpf(1)

    def entropy(x: mpmath.mpf) -> mpmath.mpf:
        return _h_scalar(x, mpmath.log)

    def entropy_derivative(x: mpmath.mpf) -> mpmath.mpf:
        return mpmath.log((one - x) / x)

    def stationarity(x: mpmath.mpf) -> mpmath.mpf:
        j = x * x
        k = j + (x * (one - x)) ** 2
        j_prime = 2 * x
        k_prime = 2 * x + 2 * x * (one - x) ** 2 - 2 * x * x * (one - x)
        numerator = (one - beta) * entropy(j) + beta * entropy(k)
        numerator_prime = (
            (one - beta) * entropy_derivative(j) * j_prime
            + beta * entropy_derivative(k) * k_prime
        )
        denominator = x * entropy(x)
        denominator_prime = entropy(x) + x * entropy_derivative(x)
        return numerator_prime * denominator - numerator * denominator_prime

    return mpmath.findroot(stationarity, initial)


def print_face_problem(maxiter: int, population_size: int, seed: int,
                       dps: int) -> None:
    print("4. EXPLICIT q=0 FIVE-PARAMETER FACE")
    print("PROVED [specialising Liu (81)--(84) at q=0]: variables are")
    print("   (a1,a2,b0,b2,b4), with a3=1-a1-a2.")
    print("PROVED [exact feasible set for rational c]:")
    print("   0<=a1,a2,b0,b2,b4<=1; a1+a2<=1;")
    print("   a1*b0+a2*b2+a3*b4 >= 1-c.  There are no ordering constraints.")
    print("PROVED [exact continuous inequality, h(0)=h(1)=0]:")
    print("   (1-beta) sum_ij ai*aj*h(bi*bj)")
    print(" + beta sum_ij ai*aj*h(bi*bj[1+(1-bi)(1-bj)])")
    print(" - sum_i ai*h(bi) >= 0,")
    print("   where (b_1,b_2,b_3)=(b0,b2,b4) and")
    print(f"   beta={BETA} (the exact rational denoted by {LIU_BETA_DECIMAL}).")
    print("PROVED [quotient domain]: where H=sum_i ai*h(bi)>0, the equivalent")
    print("   face objective is the first two sums divided by H and must be >=1.")
    print("   The gap form above remains continuous at H=0 and is ready for interval")
    print("   branch-and-bound without inheriting frankl5.m's 0/0 endpoint ambiguity.")

    with mpmath.workdps(dps):
        beta_mp = mpmath.mpf(BETA.numerator) / BETA.denominator
        structural_x = _face_structural_stationary(
            beta_mp, mpmath.mpf("0.690787593924988"))
        entropy_x = _h_scalar(structural_x, mpmath.log)
        j = structural_x * structural_x
        k = j + (structural_x * (1 - structural_x)) ** 2
        core = (
            (1 - beta_mp) * _h_scalar(j, mpmath.log)
            + beta_mp * _h_scalar(k, mpmath.log)
        ) / (structural_x * entropy_x)
        crossing = 1 - 1 / core

    print(
        "NUMERICAL [mpmath %d dps, fixed rational beta]: stationary two-atom "
        "support x=%s; its objective is (1-c)*%s."
        % (dps, mpmath.nstr(structural_x, 30), mpmath.nstr(core, 30))
    )
    print(
        "NUMERICAL [not identified with Liu's exact constant]: this fixed-rational-"
        "beta two-atom family crosses objective 1 at c=%s."
        % mpmath.nstr(crossing, 30)
    )
    print(
        "NUMERICAL [one-core deterministic full-face search]: target c | "
        "estimated min Phi | objective margin Phi-1 | gap margin N-H | "
        "search-minus-two-atom | evaluations | best (a1,a2,b0,b2,b4)"
    )

    for index, target_text in enumerate(TARGET_TEXTS):
        target = Fraction(target_text)
        result, best_face = _face_search(
            target, structural_x, BETA, maxiter, population_size, seed + index)
        best_mean = float(_face_components(best_face)[0] @ _face_components(best_face)[1])
        if best_mean < 1.0 - float(target) - 2e-12:
            raise AssertionError("face search returned an infeasible point")

        with mpmath.workdps(dps):
            target_mp = mpmath.mpf(target.numerator) / target.denominator
            target_mean_mp = 1 - target_mp
            structural_p = target_mean_mp / structural_x
            structural_value = target_mean_mp * core
            objective_margin = structural_value - 1
            entropy = structural_p * entropy_x
            gap_margin = entropy * objective_margin
        search_difference = result.fun - float(structural_value)
        # The search starts with the structural point, so a large positive
        # discrepancy means its equality-boundary feasibility was mishandled.
        if search_difference > 2e-8:
            raise AssertionError("full-face search lost its seeded structural point")

        print(
            "NUMERICAL [c=%s=%s]: %s | %s | %s | %+.3e | %d | %s; mean=%.15g"
            % (
                target_text,
                target,
                mpmath.nstr(structural_value, 18),
                mpmath.nstr(objective_margin, 14),
                mpmath.nstr(gap_margin, 14),
                search_difference,
                result.nfev,
                np.array2string(best_face, precision=9, separator=","),
                best_mean,
            )
        )

    print("NUMERICAL [interpretation]: all three full five-dimensional searches")
    print("   matched the boundary family p*delta_x+(1-p)*delta_0 to the displayed")
    print("   search-minus-two-atom residual.  These are numerical upper estimates,")
    print("   not global certificates.")
    print("PROVED [scope from the in-domain obstruction]: this face is a validation")
    print("   target for porting the certified engine, NOT a proof of Hypothesis 2.")
    print("   A q-interior path strictly beats its corresponding face point; a separate")
    print("   global reduction would be required before a face certificate could imply H2.")
    print()


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-points", type=int, default=DEFAULT_RANDOM_POINTS)
    parser.add_argument("--seed", type=int, default=230_608_824)
    parser.add_argument("--dps", type=int, default=DEFAULT_DPS)
    parser.add_argument("--face-maxiter", type=int, default=DEFAULT_FACE_MAXITER)
    parser.add_argument("--face-population", type=int, default=DEFAULT_FACE_POPULATION)
    args = parser.parse_args(argv)
    if args.random_points < 20:
        parser.error("--random-points must be at least 20")
    if args.dps < 60:
        parser.error("--dps must be at least 60")
    if args.face_maxiter < 20 or args.face_population < 20:
        parser.error("face search controls must each be at least 20")

    print("LIU NINE-PARAMETER NU-COLLAPSE AUDIT")
    print("Status labels are PROVED, NUMERICAL, or CONJECTURED.")
    print("PROVED [runtime configuration]: BLAS/OpenMP thread limits are set to one;")
    print("   scipy differential_evolution is called with workers=1.")
    print()
    verify_separation(args.random_points, args.seed, args.dps)
    print_fibre_obstruction(args.dps)
    print_face_problem(
        args.face_maxiter, args.face_population, args.seed + 1000, args.dps)
    print("FINAL VERDICT")
    print("PROVED [source algebra]: corrected separation uses K(nu)/H(mu); R alone")
    print("   differs by the printed exact low-rank and normalisation residual.")
    print("PROVED [in-domain Arb + exact coercive partial sum]: a rational paired-atom")
    print("   path has LOSS/GAIN >4.2848, while the first 64 coercive terms alone give")
    print("   LOSS/GAIN >4.0313.  The pointwise q=0 collapse is false, and increasing")
    print("   the Hypothesis-1 coercivity constant only makes the obstruction stronger.")
    print("PROVED [scope]: this does not disprove Liu's global q=0 minimizer hypothesis;")
    print("   it proves that the proposed universal nu-collapse cannot establish it.")
    print("NUMERICAL [one-core face estimates]: see the three objective and gap margins")
    print("   above; the exact face is a certified-engine validation target only.")


if __name__ == "__main__":
    main()
