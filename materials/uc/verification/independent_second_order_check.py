#!/usr/bin/env python3
"""Independent exact expansion check for the simultaneous insertion family.

This file uses the expectation definitions in ``liu9_objective._formula`` but
never imports ``liu9_second_order``.  Write

    h(t) = -t log(t) - (1-t) log(1-t),
    pi(u,v) = u*v + u*(1-u)*v*(1-v),
    K(u,v) = (1-beta)*h(u*v) + beta*h(pi(u,v)),
    w = 1-q,  ybar = w*y1 + q*y2,  r = ybar/x.

For

    (a1,a2,q,b0,b2,b4,b1,b3,b5)
      = (p-eps*ybar/x, eps, q, x,y1,0, x,y2,0),

the aggregate mean is identically ``p*x``.  Directly collecting the measure
expectations gives

    gap(eps) = eps*Gbar + eps**2*Rasym,
    dist2(eps) = eps*Sbar + eps**2*Vasym,
    pencil(eps) = gap(eps) - kappa*dist2(eps)
                = eps*Lbar + eps**2*Qasym,

where the zero constant uses the defining Liu identity
``p**2*K(x,x)-p*h(x)=0``, and

    A = 2*p*K(x,x) - h(x),
    g(y) = 2*p*K(x,y) - h(y) - (y/x)*A,
    s(y) = y**2*(y-x)**2,
    Gbar = w*g(y1) + q*g(y2),
    Sbar = w*s(y1) + q*s(y2),
    Lbar = Gbar - kappa*Sbar,
    Vasym = w*(y1-ybar)**2 + q*(y2-ybar)**2
           = w*q*(y1-y2)**2,
    Bq(y1,y2) = (1-beta)*(
        w**2*h(y1**2) + 2*w*q*h(y1*y2) + q**2*h(y2**2)
    ) + beta*(
        w*h(pi(y1,y1)) + q*h(pi(y2,y2))
    ),
    Rasym = r**2*K(x,x)
            - 2*r*(w*K(x,y1) + q*K(x,y2))
            + Bq(y1,y2),
    Qasym(y1,y2,q) = Rasym - kappa*w*q*(y1-y2)**2.

This is an exact polynomial, not a Taylor truncation.  Epsilon occurs only in
the three masses, which are affine; the two pair expectations are quadratic
in those masses, and the one-point expectation is linear.  Every argument of
``h`` is one of the fixed support expressions displayed above.  Consequently
there is no higher-order term and, in particular, no
``eps**2*log(1/eps)`` term: epsilon never opens an entropy argument at 0 or 1.

Run from the math directory with

    ./.venv/bin/python -I -B uc/verification/independent_second_order_check.py
"""

from __future__ import annotations

import heapq
import math
import os
import random
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Sequence

# This workstation is shared with live compute workers.  These must be set
# before importing anything that can transitively load a numerical library.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

HERE = Path(__file__).resolve().parent
UC = HERE.parent
if str(UC) not in sys.path:
    sys.path.insert(0, str(UC))

import mpmath  # noqa: E402

from liu9_binding import MPParameters, solve_equation_parameters  # noqa: E402
from liu9_boundary_layer import gap_mp  # noqa: E402
from liu9_ninevar import distance_squared  # noqa: E402


DPS = 110
KAPPA = Fraction(4_119_063, 33_554_432)
VERIFY_CONFIGURATIONS = 36
EPSILON_SCALES = ("1e-1", "1e-4", "1e-8", "1e-13", "1e-20")
SEARCH_BASE_SAMPLES = 210_000
VERIFY_SEED = 0x230608824
SEARCH_SEED = 0x20260828
SHORTLIST_SIZE = 48
MP_ZERO_TOLERANCE = mpmath.mpf("1e-90")
FLOAT_ZERO_TOLERANCE = 5e-14


@dataclass(frozen=True, slots=True)
class Expansion:
    raw_linear: Any
    raw_quadratic: Any
    distance_linear: Any
    distance_quadratic: Any
    pencil_linear: Any
    pencil_quadratic: Any

    def raw_gap(self, epsilon: Any) -> Any:
        return epsilon * self.raw_linear + epsilon * epsilon * self.raw_quadratic

    def dist2(self, epsilon: Any) -> Any:
        return (epsilon * self.distance_linear
                + epsilon * epsilon * self.distance_quadratic)

    def pencil(self, epsilon: Any) -> Any:
        return (epsilon * self.pencil_linear
                + epsilon * epsilon * self.pencil_quadratic)


@dataclass(frozen=True, slots=True)
class Candidate:
    label: str
    y1: float
    y2: float
    q: float
    epsilon: float
    epsilon_max: float
    raw_gap: float
    pencil: float
    reduced_pencil: float
    local: bool


@dataclass(frozen=True, slots=True)
class MPCandidate:
    source: Candidate
    raw_gap: mpmath.mpf
    distance_squared: mpmath.mpf
    pencil: mpmath.mpf
    reduced_pencil: mpmath.mpf
    mean_residual: mpmath.mpf
    reconstruction_residual: mpmath.mpf


def _mp_fraction(value: Fraction) -> mpmath.mpf:
    return mpmath.mpf(value.numerator) / value.denominator


def entropy_mp(value: mpmath.mpf) -> mpmath.mpf:
    """Binary entropy with the endpoint convention used by the SSOT."""
    if value == 0 or value == 1:
        return mpmath.mpf(0)
    if not 0 < value < 1:
        raise ValueError(f"entropy argument outside [0,1]: {value!r}")
    return -(value * mpmath.log(value)
             + (1 - value) * mpmath.log(1 - value))


def entropy_float(value: float) -> float:
    if value == 0.0 or value == 1.0:
        return 0.0
    if not 0.0 < value < 1.0:
        raise ValueError(f"entropy argument outside [0,1]: {value!r}")
    return -(value * math.log(value) + (1.0 - value) * math.log1p(-value))


def protocol(left: Any, right: Any) -> Any:
    return left * right + left * (1 - left) * right * (1 - right)


def kernel(left: Any, right: Any, beta: Any,
           entropy: Callable[[Any], Any]) -> Any:
    return ((1 - beta) * entropy(left * right)
            + beta * entropy(protocol(left, right)))


def derive_expansion(
    x: Any,
    p: Any,
    beta: Any,
    kappa: Any,
    y1: Any,
    y2: Any,
    q: Any,
    entropy: Callable[[Any], Any],
) -> Expansion:
    """Collect the two measure expectations by powers of epsilon."""
    w = 1 - q
    ybar = w * y1 + q * y2
    ratio = ybar / x
    k_xx = kernel(x, x, beta, entropy)
    k_xy1 = kernel(x, y1, beta, entropy)
    k_xy2 = kernel(x, y2, beta, entropy)
    kbar = w * k_xy1 + q * k_xy2
    hbar = w * entropy(y1) + q * entropy(y2)
    active_derivative = 2 * p * k_xx - entropy(x)

    raw_linear = 2 * p * kbar - hbar - ratio * active_derivative
    pair_inserted = (1 - beta) * (
        w * w * entropy(y1 * y1)
        + 2 * w * q * entropy(y1 * y2)
        + q * q * entropy(y2 * y2)
    ) + beta * (
        w * entropy(protocol(y1, y1))
        + q * entropy(protocol(y2, y2))
    )
    raw_quadratic = (
        ratio * ratio * k_xx - 2 * ratio * kbar + pair_inserted
    )

    distance_linear = (
        w * y1 * y1 * (y1 - x) * (y1 - x)
        + q * y2 * y2 * (y2 - x) * (y2 - x)
    )
    distance_quadratic = w * q * (y1 - y2) * (y1 - y2)
    return Expansion(
        raw_linear=raw_linear,
        raw_quadratic=raw_quadratic,
        distance_linear=distance_linear,
        distance_quadratic=distance_quadratic,
        pencil_linear=raw_linear - kappa * distance_linear,
        pencil_quadratic=raw_quadratic - kappa * distance_quadratic,
    )


def family_values(
    parameters: MPParameters,
    y1: mpmath.mpf,
    y2: mpmath.mpf,
    q: mpmath.mpf,
    epsilon: mpmath.mpf,
) -> tuple[mpmath.mpf, ...]:
    """Build the simultaneous insertion without any second-order helper."""
    zero = parameters.x * 0
    ybar = (1 - q) * y1 + q * y2
    return (
        parameters.p - epsilon * ybar / parameters.x,
        epsilon,
        q,
        parameters.x,
        y1,
        zero,
        parameters.x,
        y2,
        zero,
    )


def aggregate_mean(values: Sequence[Any]) -> Any:
    """Direct expectation of the six-atom support measure."""
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    a3 = 1 - a1 - a2
    return ((1 - q) * (a1 * b0 + a2 * b2 + a3 * b4)
            + q * (a1 * b1 + a2 * b3 + a3 * b5))


def epsilon_max(ybar: Any, p: Any, x: Any) -> Any:
    """Largest epsilon for which all three affine masses are nonnegative."""
    one = x * 0 + 1
    caps = [one]
    if ybar > 0:
        caps.append(p * x / ybar)  # a1 >= 0
    if ybar < x:
        caps.append(x * (1 - p) / (x - ybar))  # a3 >= 0
    return min(caps)


def _fmt(value: Any, digits: int = 26) -> str:
    if isinstance(value, float):
        return format(value, ".17g")
    return mpmath.nstr(value, digits)


def calibrate_zeros(parameters: MPParameters, kappa: mpmath.mpf) -> None:
    """Calibrate both independent transcriptions before any random check."""
    x = parameters.x
    q = mpmath.mpf("0.371")
    zero = mpmath.mpf(0)
    expansion = derive_expansion(
        x, parameters.p, parameters.beta, kappa, x, x, q, entropy_mp
    )
    base_values = family_values(parameters, x, x, q, zero)
    base_gap = gap_mp(base_values, parameters.beta)
    base_distance = distance_squared(base_values, parameters)
    base_pencil = base_gap - kappa * base_distance

    maximum_gap = abs(base_gap)
    maximum_distance = abs(base_distance)
    maximum_pencil = abs(base_pencil)
    maximum_gap_change = mpmath.mpf(0)
    for fraction in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1)):
        epsilon = parameters.p * _mp_fraction(fraction)
        values = family_values(parameters, x, x, q, epsilon)
        raw_gap = gap_mp(values, parameters.beta)
        dist2 = distance_squared(values, parameters)
        pencil = raw_gap - kappa * dist2
        maximum_gap = max(maximum_gap, abs(raw_gap))
        maximum_distance = max(maximum_distance, abs(dist2))
        maximum_pencil = max(maximum_pencil, abs(pencil))
        maximum_gap_change = max(maximum_gap_change, abs(raw_gap - base_gap))

    coefficient_zero = max(
        abs(expansion.raw_linear),
        abs(expansion.raw_quadratic),
        abs(expansion.distance_linear),
        abs(expansion.distance_quadratic),
        abs(expansion.pencil_linear),
        abs(expansion.pencil_quadratic),
    )
    x_float = float(parameters.x)
    p_float = float(parameters.p)
    beta_float = float(parameters.beta)
    kappa_float = float(kappa)
    float_expansion = derive_expansion(
        x_float,
        p_float,
        beta_float,
        kappa_float,
        x_float,
        x_float,
        float(q),
        entropy_float,
    )
    float_optimizer_gap = (
        p_float * p_float * kernel(
            x_float, x_float, beta_float, entropy_float
        )
        - p_float * entropy_float(x_float)
    )
    float_coefficient_zero = max(
        abs(float_expansion.raw_linear),
        abs(float_expansion.raw_quadratic),
        abs(float_expansion.distance_linear),
        abs(float_expansion.distance_quadratic),
        abs(float_expansion.pencil_linear),
        abs(float_expansion.pencil_quadratic),
    )
    print("\n[1] ZERO CALIBRATION")
    print("optimizer raw gap              =", _fmt(base_gap))
    print("optimizer distance squared      =", _fmt(base_distance))
    print("optimizer pencil                =", _fmt(base_pencil))
    print("max |raw gap| for y1=y2=x       =", _fmt(maximum_gap))
    print("max |gap(eps)-gap(0)| at y=x    =", _fmt(maximum_gap_change))
    print("max distance squared at y=x     =", _fmt(maximum_distance))
    print("max |pencil| at y=x             =", _fmt(maximum_pencil))
    print("max independent coefficient at x=", _fmt(coefficient_zero))
    print("float optimizer-gap calibration =", _fmt(float_optimizer_gap))
    print("float coefficient calibration   =", _fmt(float_coefficient_zero))
    if max(maximum_gap, maximum_distance, maximum_pencil,
           coefficient_zero) > MP_ZERO_TOLERANCE:
        raise AssertionError("zero calibration failed")
    if max(abs(float_optimizer_gap),
           float_coefficient_zero) > FLOAT_ZERO_TOLERANCE:
        raise AssertionError("float search transcription failed calibration")
    print("calibration verdict             = PASS")


def _seeded_verification_points(
    count: int,
    seed: int,
) -> list[tuple[mpmath.mpf, mpmath.mpf, mpmath.mpf]]:
    rng = random.Random(seed)
    denominator = 1 << 53
    points = []
    for index in range(count):
        y1 = mpmath.mpf(rng.getrandbits(53)) / denominator
        y2 = mpmath.mpf(rng.getrandbits(53)) / denominator
        q = mpmath.mpf(rng.getrandbits(53)) / denominator
        # Deterministically add endpoint clustering while retaining a seeded
        # interior point in every coordinate.
        if index % 4 == 1:
            y1 = y1**4
        elif index % 4 == 2:
            y2 = 1 - (1 - y2) ** 4
        elif index % 4 == 3:
            q = q**3
        points.append((y1, y2, q))
    return points


def verify_exact_reconstruction(
    parameters: MPParameters,
    kappa: mpmath.mpf,
) -> None:
    """Compare the independently collected polynomial with both SSOT calls."""
    points = _seeded_verification_points(VERIFY_CONFIGURATIONS, VERIFY_SEED)
    scales = tuple(mpmath.mpf(text) for text in EPSILON_SCALES)
    maximum_gap_residual = mpmath.mpf(0)
    maximum_distance_residual = mpmath.mpf(0)
    maximum_pencil_residual = mpmath.mpf(0)
    maximum_mean_residual = mpmath.mpf(0)
    maximum_zero_gap = mpmath.mpf(0)
    expansions: list[tuple[mpmath.mpf, mpmath.mpf, mpmath.mpf, Expansion]] = []

    for y1, y2, q in points:
        expansion = derive_expansion(
            parameters.x,
            parameters.p,
            parameters.beta,
            kappa,
            y1,
            y2,
            q,
            entropy_mp,
        )
        expansions.append((y1, y2, q, expansion))
        ybar = (1 - q) * y1 + q * y2
        feasible_maximum = epsilon_max(ybar, parameters.p, parameters.x)
        zero_values = family_values(parameters, y1, y2, q, mpmath.mpf(0))
        zero_gap = gap_mp(zero_values, parameters.beta)
        zero_distance = distance_squared(zero_values, parameters)
        maximum_zero_gap = max(maximum_zero_gap, abs(zero_gap))
        if abs(zero_distance) > MP_ZERO_TOLERANCE:
            raise AssertionError("distance did not calibrate to zero at epsilon=0")

        for scale in scales:
            epsilon = feasible_maximum * scale
            values = family_values(parameters, y1, y2, q, epsilon)
            raw_gap = gap_mp(values, parameters.beta)
            dist2 = distance_squared(values, parameters)
            centered_gap = raw_gap - zero_gap
            centered_distance = dist2 - zero_distance
            pencil = centered_gap - kappa * centered_distance
            expected_gap = expansion.raw_gap(epsilon)
            expected_distance = expansion.dist2(epsilon)
            expected_pencil = expansion.pencil(epsilon)
            maximum_gap_residual = max(
                maximum_gap_residual, abs(centered_gap - expected_gap)
            )
            maximum_distance_residual = max(
                maximum_distance_residual,
                abs(centered_distance - expected_distance),
            )
            maximum_pencil_residual = max(
                maximum_pencil_residual, abs(pencil - expected_pencil)
            )
            maximum_mean_residual = max(
                maximum_mean_residual,
                abs(aggregate_mean(values) - parameters.mean),
            )

    # Import the reference only after all independent coefficients and SSOT
    # reconstructions have been computed.
    from liu9_ninevar import transverse_coefficients_mp  # noqa: PLC0415

    maximum_first_order_disagreement = mpmath.mpf(0)
    for y1, y2, q, expansion in expansions:
        l1 = transverse_coefficients_mp(
            parameters, y1, "mean_preserving"
        )[0]
        l2 = transverse_coefficients_mp(
            parameters, y2, "mean_preserving"
        )[0]
        reference = (1 - q) * l1 + q * l2
        maximum_first_order_disagreement = max(
            maximum_first_order_disagreement,
            abs(expansion.pencil_linear - reference),
        )

    maximum_formula_residual = max(
        maximum_gap_residual,
        maximum_distance_residual,
        maximum_pencil_residual,
    )
    print("\n[2] EXACT EPSILON RECONSTRUCTION")
    print("working precision (decimal dps) =", DPS)
    print("seeded configurations            =", len(points))
    print("epsilon scales                   =", ", ".join(EPSILON_SCALES))
    print("SSOT evaluations                 =", len(points) * len(scales))
    print("max epsilon=0 raw-gap residual   =", _fmt(maximum_zero_gap))
    print("max aggregate-mean residual      =", _fmt(maximum_mean_residual))
    print("max raw-gap reconstruction resid =", _fmt(maximum_gap_residual))
    print("max distance reconstruction resid=", _fmt(maximum_distance_residual))
    print("max pencil reconstruction resid  =", _fmt(maximum_pencil_residual))
    print("maximum formula residual         =", _fmt(maximum_formula_residual))
    print("max first-order disagreement     =", _fmt(maximum_first_order_disagreement))
    if maximum_formula_residual > MP_ZERO_TOLERANCE:
        raise AssertionError("the exact epsilon reconstruction failed")
    if maximum_mean_residual > MP_ZERO_TOLERANCE:
        raise AssertionError("the simultaneous family did not preserve the mean")
    if maximum_first_order_disagreement > MP_ZERO_TOLERANCE:
        raise AssertionError("first-order coefficient disagrees with Lbar")
    print("reconstruction verdict           = PASS (exact quadratic; no log term)")
    print("first-order Lbar comparison      = PASS")


def _sample_geometry(
    rng: random.Random,
    index: int,
    x: float,
) -> tuple[float, float, float]:
    """Seeded mixture of interior, endpoint, diagonal, and q-boundary draws."""
    u, v, z = rng.random(), rng.random(), rng.random()
    mode = index % 12
    if mode == 0:
        return u, v, z
    if mode == 1:
        return u**4, v, z
    if mode == 2:
        return 1.0 - u**4, 1.0 - v**4, z
    if mode == 3:
        return u, 1.0 - v**8, z
    if mode == 4:
        return 1.0, v, z
    if mode == 5:
        return u, 1.0, z
    if mode == 6:
        return u, v, z**5
    if mode == 7:
        return u, v, 1.0 - z**5
    if mode == 8:
        delta = (v - 0.5) * 2e-5
        return min(1.0, max(0.0, u + delta)), u, z
    if mode == 9:
        delta = (u - 0.5) * 2e-5
        return min(1.0, max(0.0, x + delta)), x, z
    if mode == 10:
        return 0.0, v, z
    return u, 0.0, z


def _sample_epsilon_fraction(rng: random.Random, index: int) -> float:
    u = rng.random()
    mode = index % 6
    if mode == 0:
        return u
    if mode == 1:
        return u**4
    if mode == 2:
        return 1.0 - u**4
    if mode == 3:
        return math.sqrt(u)
    if mode == 4:
        return 0.5 + 0.5 * u
    return 1.0 - 10.0 ** (-2.0 - 10.0 * u)


def _push_shortlist(
    heap: list[tuple[float, int, Candidate]],
    score: float,
    serial: int,
    candidate: Candidate,
) -> None:
    """Keep the candidates with the smallest score in a bounded max-heap."""
    item = (-score, serial, candidate)
    if len(heap) < SHORTLIST_SIZE:
        heapq.heappush(heap, item)
    elif score < -heap[0][0]:
        heapq.heapreplace(heap, item)


def _float_candidate(
    label: str,
    y1: float,
    y2: float,
    q: float,
    epsilon: float,
    feasible_maximum: float,
    expansion: Expansion,
    local: bool,
) -> Candidate:
    raw_gap = expansion.raw_gap(epsilon)
    pencil = expansion.pencil(epsilon)
    reduced = expansion.pencil_linear + epsilon * expansion.pencil_quadratic
    return Candidate(
        label=label,
        y1=y1,
        y2=y2,
        q=q,
        epsilon=epsilon,
        epsilon_max=feasible_maximum,
        raw_gap=raw_gap,
        pencil=pencil,
        reduced_pencil=reduced,
        local=local,
    )


def _to_mp(value: float) -> mpmath.mpf:
    return mpmath.mpf(repr(value))


def reevaluate_candidate(
    candidate: Candidate,
    parameters: MPParameters,
    kappa: mpmath.mpf,
) -> MPCandidate:
    y1 = _to_mp(candidate.y1)
    y2 = _to_mp(candidate.y2)
    q = _to_mp(candidate.q)
    epsilon = _to_mp(candidate.epsilon)
    values = family_values(parameters, y1, y2, q, epsilon)
    raw_gap = gap_mp(values, parameters.beta)
    dist2 = distance_squared(values, parameters)
    pencil = raw_gap - kappa * dist2
    expansion = derive_expansion(
        parameters.x,
        parameters.p,
        parameters.beta,
        kappa,
        y1,
        y2,
        q,
        entropy_mp,
    )
    expected = expansion.pencil(epsilon)
    reduced = expansion.pencil_linear + epsilon * expansion.pencil_quadratic
    return MPCandidate(
        source=candidate,
        raw_gap=raw_gap,
        distance_squared=dist2,
        pencil=pencil,
        reduced_pencil=reduced,
        mean_residual=abs(aggregate_mean(values) - parameters.mean),
        reconstruction_residual=abs(pencil - expected),
    )


def _print_mp_candidate(title: str, candidate: MPCandidate) -> None:
    source = candidate.source
    print(title)
    print("  source       =", source.label)
    print("  y1           =", _fmt(_to_mp(source.y1), 18))
    print("  y2           =", _fmt(_to_mp(source.y2), 18))
    print("  q            =", _fmt(_to_mp(source.q), 18))
    print("  epsilon      =", _fmt(_to_mp(source.epsilon), 18))
    print("  epsilon_max  =", _fmt(_to_mp(source.epsilon_max), 18))
    print("  local sample =", source.local)
    print("  raw gap      =", _fmt(candidate.raw_gap))
    print("  dist2        =", _fmt(candidate.distance_squared))
    print("  pencil       =", _fmt(candidate.pencil))
    print("  reduced P/eps=", _fmt(candidate.reduced_pencil))
    print("  mean residual=", _fmt(candidate.mean_residual))


def search_seeded_family(
    parameters: MPParameters,
    kappa_mp: mpmath.mpf,
) -> tuple[MPCandidate, MPCandidate, MPCandidate, int, int]:
    """Fast float search followed by authoritative high-precision replay."""
    rng = random.Random(SEARCH_SEED)
    x = float(parameters.x)
    p = float(parameters.p)
    beta = float(parameters.beta)
    kappa = float(kappa_mp)
    raw_heap: list[tuple[float, int, Candidate]] = []
    pencil_heap: list[tuple[float, int, Candidate]] = []
    local_heap: list[tuple[float, int, Candidate]] = []
    raw_negative_screen = 0
    pencil_negative_screen = 0
    serial = 0
    local_domain_checked = 0

    def record(candidate: Candidate) -> None:
        nonlocal raw_negative_screen, pencil_negative_screen
        nonlocal serial, local_domain_checked
        serial += 1
        _push_shortlist(raw_heap, candidate.raw_gap, serial, candidate)
        _push_shortlist(pencil_heap, candidate.pencil, serial, candidate)
        if candidate.local:
            local_domain_checked += 1
            _push_shortlist(
                local_heap, candidate.reduced_pencil, serial, candidate
            )
        if candidate.raw_gap < -1e-12:
            raw_negative_screen += 1
        if candidate.pencil < -1e-12:
            pencil_negative_screen += 1

    for index in range(SEARCH_BASE_SAMPLES):
        y1, y2, q = _sample_geometry(rng, index, x)
        ybar = (1.0 - q) * y1 + q * y2
        feasible_maximum = float(epsilon_max(ybar, p, x))
        expansion = derive_expansion(
            x, p, beta, kappa, y1, y2, q, entropy_float
        )
        full_epsilon = feasible_maximum * _sample_epsilon_fraction(rng, index)
        local_cap = min(0.5, feasible_maximum)
        local_epsilon = local_cap * _sample_epsilon_fraction(rng, index + 1)
        record(_float_candidate(
            "seeded-full",
            y1,
            y2,
            q,
            full_epsilon,
            feasible_maximum,
            expansion,
            full_epsilon <= local_cap,
        ))
        record(_float_candidate(
            "seeded-local",
            y1,
            y2,
            q,
            local_epsilon,
            feasible_maximum,
            expansion,
            True,
        ))

    # Exact strata are included in addition to the seeded random population.
    anchors = (
        ("anchor-x-full", x, x, 0.5, p, False),
        ("anchor-x-local", x, x, 0.5, min(0.5, p), True),
        ("anchor-endpoint", 1.0, 1.0, 0.5, p * x, False),
        ("anchor-zero", 0.0, 0.0, 0.5, 1.0 - p, True),
    )
    for label, y1, y2, q, epsilon, _local in anchors:
        ybar = (1.0 - q) * y1 + q * y2
        feasible_maximum = float(epsilon_max(ybar, p, x))
        expansion = derive_expansion(
            x, p, beta, kappa, y1, y2, q, entropy_float
        )
        record(_float_candidate(
            label,
            y1,
            y2,
            q,
            epsilon,
            feasible_maximum,
            expansion,
            epsilon <= min(0.5, feasible_maximum),
        ))

    raw_replays = [
        reevaluate_candidate(item[2], parameters, kappa_mp)
        for item in raw_heap
    ]
    pencil_replays = [
        reevaluate_candidate(item[2], parameters, kappa_mp)
        for item in pencil_heap
    ]
    local_replays = [
        reevaluate_candidate(item[2], parameters, kappa_mp)
        for item in local_heap
    ]
    raw_minimum = min(raw_replays, key=lambda item: item.raw_gap)
    pencil_minimum = min(pencil_replays, key=lambda item: item.pencil)
    local_minimum = min(local_replays, key=lambda item: item.reduced_pencil)
    maximum_replay_residual = max(
        item.reconstruction_residual
        for item in raw_replays + pencil_replays + local_replays
    )
    maximum_mean_residual = max(
        item.mean_residual
        for item in raw_replays + pencil_replays + local_replays
    )

    print("\n[3] SEEDED MEAN-FEASIBLE SEARCH")
    print("seed                           =", SEARCH_SEED)
    print("seeded geometries              =", SEARCH_BASE_SAMPLES)
    print("full-range epsilon samples     =", SEARCH_BASE_SAMPLES)
    print("local-range epsilon samples    =", SEARCH_BASE_SAMPLES)
    print("explicit stratum anchors       =", len(anchors))
    print("total searched configurations  =", serial)
    print("total local-domain checked     =", local_domain_checked)
    print("float-screen raw gap < -1e-12  =", raw_negative_screen)
    print("float-screen pencil < -1e-12   =", pencil_negative_screen)
    print("max shortlist replay residual  =", _fmt(maximum_replay_residual))
    print("max shortlist mean residual    =", _fmt(maximum_mean_residual))
    _print_mp_candidate("\nraw-gap minimum (SSOT replay):", raw_minimum)
    _print_mp_candidate("\npencil minimum (SSOT replay):", pencil_minimum)
    _print_mp_candidate("\nlocal reduced-pencil minimum:", local_minimum)
    return (
        raw_minimum,
        pencil_minimum,
        local_minimum,
        raw_negative_screen,
        pencil_negative_screen,
    )


def verify_endpoint_obstruction(
    parameters: MPParameters,
    kappa: mpmath.mpf,
) -> None:
    """Check the sharp all-{0,1} endpoint, simplifying its zero mass exactly."""
    zero = mpmath.mpf(0)
    one = mpmath.mpf(1)
    q = mpmath.mpf("0.5")
    epsilon = parameters.p * parameters.x
    # At ybar=1 and eps=p*x, a1=p-eps/x=0 algebraically.  Writing that
    # simplified zero avoids mistaking a final-roundoff remnant for mass.
    values = (
        zero,
        epsilon,
        q,
        parameters.x,
        one,
        zero,
        parameters.x,
        one,
        zero,
    )
    raw_gap = gap_mp(values, parameters.beta)
    dist2 = distance_squared(values, parameters)
    pencil = raw_gap - kappa * dist2
    mean_residual = abs(aggregate_mean(values) - parameters.mean)
    expansion = derive_expansion(
        parameters.x,
        parameters.p,
        parameters.beta,
        kappa,
        one,
        one,
        q,
        entropy_mp,
    )
    expansion_residual = abs(expansion.pencil(epsilon) - pencil)
    unsimplified_a1 = family_values(
        parameters, one, one, q, epsilon
    )[0]

    print("\n[4] SHARP ENDPOINT PENCIL OBSTRUCTION")
    print("configuration: y1=y2=1, q=1/2 (arbitrary), epsilon=p*x")
    print("epsilon = p*x                    =", _fmt(epsilon))
    print("masses (a1,a2,a3)               = (0, p*x, 1-p*x)")
    print("unsimplified numerical a1 resid =", _fmt(unsimplified_a1))
    print("aggregate-mean residual          =", _fmt(mean_residual))
    print("raw gap                           =", _fmt(raw_gap))
    print("distance squared                  =", _fmt(dist2))
    print("pencil                            =", _fmt(pencil))
    print("expansion reconstruction residual =", _fmt(expansion_residual))
    print("raw-gap-zero verdict              = PASS")
    print("negative-pencil verdict           = PASS")
    if abs(raw_gap) > MP_ZERO_TOLERANCE:
        raise AssertionError("endpoint raw gap was not zero")
    if not pencil < 0:
        raise AssertionError("endpoint pencil was not negative")
    if mean_residual > MP_ZERO_TOLERANCE:
        raise AssertionError("endpoint did not preserve the mean")
    if expansion_residual > MP_ZERO_TOLERANCE:
        raise AssertionError("endpoint did not satisfy the exact expansion")


def main() -> int:
    started = time.monotonic()
    mpmath.mp.dps = DPS
    parameters = solve_equation_parameters(DPS)
    kappa = _mp_fraction(KAPPA)

    print("INDEPENDENT SIMULTANEOUS-INSERTION SECOND-ORDER CHECK")
    print("=" * 61)
    print("pencil formula: eps*Lbar + eps^2*Qasym (exact)")
    print("Qasym = Rasym - kappa*(1-q)*q*(y1-y2)^2")
    print("entropy arguments are epsilon-independent; no eps^2*log(1/eps)")
    print("x     =", _fmt(parameters.x, 52))
    print("p     =", _fmt(parameters.p, 52))
    print("beta  =", _fmt(parameters.beta, 52))
    print("kappa =", _fmt(kappa, 32))

    calibrate_zeros(parameters, kappa)
    verify_exact_reconstruction(parameters, kappa)
    _, _, local_minimum, _, _ = search_seeded_family(parameters, kappa)
    verify_endpoint_obstruction(parameters, kappa)

    print("\n[5] SAMPLED LOCAL STATEMENT")
    print("domain: epsilon <= min(1/2, epsilon_max(ybar))")
    print("minimum replayed reduced pencil =", _fmt(local_minimum.reduced_pencil))
    print("minimum replayed local pencil   =", _fmt(local_minimum.pencil))
    if local_minimum.reduced_pencil < 0:
        _print_mp_candidate("NEGATIVE LOCAL CONFIGURATION:", local_minimum)
        print("local sampled statement verdict = FAIL")
        return 1
    print("local sampled statement verdict = PASS (no negative sample)")
    print("\nFINAL VERDICT: all exact-expansion checks pass; the raw gap and")
    print("pencil were tracked separately, including the sharp endpoint obstruction.")
    print("elapsed seconds =", format(time.monotonic() - started, ".3f"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
