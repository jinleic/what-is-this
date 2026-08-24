#!/usr/bin/env python3
"""Certified log-cumulant obstructions for positive subset-product spectra.

For a positive full n-mode subset-product spectrum, determinant centering turns
its empirical log spectrum into a sum of independent Rademacher variables.  If
``kappa_r`` denotes its centered empirical cumulants, then

    p1 = kappa_2,
    p2 = -kappa_4/2,
    p3 = kappa_6/16,
    p4 = -kappa_8/272

are the power sums of the nonnegative squared mode energies.  Consequently the
Stieltjes/Hankel minor ``p1*p3-p2**2`` must be nonnegative.

This producer certifies that the minor is strictly negative for the open 2x2
and 2x3 bipartite layer graphs at the benchmark-independent physical point
``t=tanh(K*/2)=1/3``.  It also certifies the sign on a small explicit rational
interval about that point.  NumPy proposes an eigenbasis only: every spectral
claim is discharged by exact integer congruence bounds and exact outward-rounded
dyadic logarithm/cumulant arithmetic.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable

# Keep the exploratory eigensolver single-threaded on shared hosts.  Its output
# is only a proposal; exact arithmetic below is the certificate.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "spectral" / "bipartite_log_moments.json"
CPU_BUDGET_SECONDS = 300.0
RSS_CAP_BYTES = 2_000_000_000
POINT_T = Fraction(1, 3)
INTERVAL_RADIUS = Fraction(1, 1_000_000_000)
VECTOR_SCALE = 10**14
DYADIC_BITS = 160
DYADIC_SCALE = 1 << DYADIC_BITS
LOG_SERIES_TERMS = 48

Edge = tuple[int, int]
DyadicInterval = tuple[int, int]  # endpoints divided by DYADIC_SCALE


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    cpu = time.process_time() - started
    rss = max_rss_bytes()
    if cpu > CPU_BUDGET_SECONDS:
        raise TimeoutError(f"PREFLIGHT/RESOURCE: {stage}: process CPU {cpu:.1f}s")
    if rss > RSS_CAP_BYTES:
        raise MemoryError(f"PREFLIGHT/RESOURCE: {stage}: RSS {rss} bytes")


def _check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def chain_edges(sites: int) -> tuple[Edge, ...]:
    return tuple((site, site + 1) for site in range(sites - 1))


def grid_edges(rows: int, columns: int) -> tuple[Edge, ...]:
    edges: list[Edge] = []
    for row in range(rows):
        for column in range(columns):
            site = row * columns + column
            if row + 1 < rows:
                edges.append((site, (row + 1) * columns + column))
            if column + 1 < columns:
                edges.append((site, site + 1))
    return tuple(edges)


def bipartite_coloring(sites: int, edges: Iterable[Edge]) -> tuple[int, ...]:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        if left == right or not (0 <= left < sites and 0 <= right < sites):
            raise ValueError("edges must be simple and in range")
        adjacency[left].append(right)
        adjacency[right].append(left)
    colors = [-1] * sites
    for root in range(sites):
        if colors[root] >= 0:
            continue
        colors[root] = 0
        queue = [root]
        for vertex in queue:
            for other in adjacency[vertex]:
                if colors[other] < 0:
                    colors[other] = 1 - colors[vertex]
                    queue.append(other)
                elif colors[other] == colors[vertex]:
                    raise ValueError("graph is not bipartite")
    return tuple(colors)


def spin_bond_data(sites: int, edges: tuple[Edge, ...]) -> tuple[list[int], int, list[int]]:
    values: list[int] = []
    for state in range(1 << sites):
        spins = [1 - 2 * ((state >> (sites - 1 - site)) & 1) for site in range(sites)]
        values.append(sum(spins[left] * spins[right] for left, right in edges))
    parity = values[0] % 2
    if any((value - parity) % 2 for value in values):
        raise AssertionError("bond-energy parity failed")
    exponents = [(value - parity) // 2 for value in values]
    return values, parity, exponents


def build_rational_layer(
    sites: int, edges: tuple[Edge, ...], t: Fraction
) -> tuple[list[list[Fraction]], int, list[list[int]], Fraction, int]:
    """Build R=P_t diag(q^m) P_t and its exact minimal common scaling.

    The physical symmetric transfer operator differs from R only by a positive
    scalar.  Here q=(1+t^2)/(2t), P_t[i,j]=t^hamming(i,j), and all operations are
    over ``Fraction``.
    """

    if not 0 < t < 1:
        raise ValueError("physical t must lie in (0,1)")
    dimension = 1 << sites
    q = (1 + t * t) / (2 * t)
    _, parity, exponents = spin_bond_data(sites, edges)
    diagonal = [q**exponent for exponent in exponents]
    t_powers = [t**power for power in range(sites + 1)]
    p_matrix = [
        [t_powers[(row ^ column).bit_count()] for column in range(dimension)]
        for row in range(dimension)
    ]
    matrix = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    for row in range(dimension):
        for column in range(row, dimension):
            value = sum(
                (
                    p_matrix[row][state]
                    * diagonal[state]
                    * p_matrix[column][state]
                    for state in range(dimension)
                ),
                start=Fraction(0),
            )
            matrix[row][column] = value
            matrix[column][row] = value

    denominator = 1
    for row in matrix:
        for value in row:
            denominator = math.lcm(denominator, value.denominator)
    integral = [[int(value * denominator) for value in row] for row in matrix]
    if any(
        Fraction(integral[row][column], denominator) != matrix[row][column]
        for row in range(dimension)
        for column in range(dimension)
    ):
        raise AssertionError("integral scaling failed")
    return matrix, denominator, integral, q, parity


def matrix_digest(matrix: list[list[int]]) -> str:
    payload = ";".join(",".join(str(value) for value in row) for row in matrix)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _round_scaled(value: float, scale: int) -> int:
    return int(np.rint(value * scale))


def exact_congruence_certificate(
    integral: list[list[int]], denominator: int
) -> tuple[dict[str, object], list[tuple[Fraction, Fraction]]]:
    """Certify all ordered eigenvalues by an exact generalized min-max bound.

    Float64 ``eigh`` proposes V and sorted diagonal d.  With rational
    V=V_int/S, set G=V^T V and B=V^T M V for M=denominator*R.  Exact row sums give

        ||G-I||_2 <= f < 1,       ||B-diag(d)||_2 <= e.

    Hence G is positive definite, V is invertible, and the generalized min-max
    principle plus Weyl gives, for every ordered eigenvalue of R,

        (d_i-e)/((1+f) denominator) <= lambda_i
          <= (d_i+e)/((1-f) denominator).

    Degeneracies require no separation assumption.
    """

    dimension = len(integral)
    float_matrix = np.asarray(integral, dtype=np.float64)
    eigenvalues, eigenvectors = np.linalg.eigh(float_matrix)
    scale = VECTOR_SCALE
    vectors = [
        [_round_scaled(eigenvectors[row, column], scale) for column in range(dimension)]
        for row in range(dimension)
    ]
    diagonal_numerators = [_round_scaled(value, scale) for value in eigenvalues]
    if diagonal_numerators != sorted(diagonal_numerators):
        raise AssertionError("rounded proposed diagonal is not ordered")

    matrix_times_vectors = [
        [
            sum(integral[row][inner] * vectors[inner][column] for inner in range(dimension))
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    gram = [
        [
            sum(vectors[inner][row] * vectors[inner][column] for inner in range(dimension))
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    congruence = [
        [
            sum(
                vectors[inner][row] * matrix_times_vectors[inner][column]
                for inner in range(dimension)
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]

    scale_squared = scale * scale
    gram_error_numerator = max(
        sum(
            abs(gram[row][column] - (scale_squared if row == column else 0))
            for column in range(dimension)
        )
        for row in range(dimension)
    )
    congruence_error_numerator = max(
        sum(
            abs(
                scale * congruence[row][column]
                - (
                    diagonal_numerators[row] * scale_squared
                    if row == column
                    else 0
                )
            )
            for column in range(dimension)
        )
        for row in range(dimension)
    )
    gram_error = Fraction(gram_error_numerator, scale_squared)
    congruence_error = Fraction(congruence_error_numerator, scale_squared * scale)
    if not gram_error < 1:
        raise AssertionError("proposed eigenbasis is not certified invertible")

    intervals: list[tuple[Fraction, Fraction]] = []
    for numerator in diagonal_numerators:
        diagonal = Fraction(numerator, scale)
        lower = (diagonal - congruence_error) / (1 + gram_error) / denominator
        upper = (diagonal + congruence_error) / (1 - gram_error) / denominator
        if lower <= 0 or lower > upper:
            raise AssertionError("positive ordered eigenvalue enclosure failed")
        intervals.append((lower, upper))

    certificate = {
        "tag": "[COMPUTATION]",
        "method": (
            "float64 eigh proposes V and d only; exact integer products certify "
            "||V^T V-I||_2<=f and ||V^T M V-diag(d)||_2<=e by symmetric infinity norms, "
            "then generalized min-max and Weyl enclose every ordered eigenvalue"
        ),
        "degeneracy_handling": (
            "no eigenvalue separation or simplicity is assumed; f<1 proves V invertible and "
            "the ordered generalized-eigenvalue bound includes multiplicities"
        ),
        "vector_scale": scale,
        "eigenbasis_rows": vectors,
        "diagonal_scale": scale,
        "diagonal_numerators": diagonal_numerators,
        "gram_error_numerator": str(gram_error_numerator),
        "gram_error_denominator": str(scale_squared),
        "congruence_error_numerator": str(congruence_error_numerator),
        "congruence_error_denominator": str(scale_squared * scale),
        "gram_error_exact": str(gram_error),
        "congruence_error_exact": str(congruence_error),
        "ordered_eigenvalue_intervals": [
            {
                "lower": str(lower),
                "upper": str(upper),
                "lower_approx": float(lower),
                "upper_approx": float(upper),
            }
            for lower, upper in intervals
        ],
    }
    return certificate, intervals


# ---------------------------------------------------------------------------
# Exact outward-rounded fixed-denominator interval arithmetic.
# ---------------------------------------------------------------------------
def _floor_ratio(numerator: int, denominator: int) -> int:
    return numerator // denominator


def _ceil_ratio(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def dyadic_from_fraction(value: Fraction) -> DyadicInterval:
    numerator = value.numerator * DYADIC_SCALE
    return (
        _floor_ratio(numerator, value.denominator),
        _ceil_ratio(numerator, value.denominator),
    )


def dyadic_add(left: DyadicInterval, right: DyadicInterval) -> DyadicInterval:
    return left[0] + right[0], left[1] + right[1]


def dyadic_neg(value: DyadicInterval) -> DyadicInterval:
    return -value[1], -value[0]


def dyadic_sub(left: DyadicInterval, right: DyadicInterval) -> DyadicInterval:
    return dyadic_add(left, dyadic_neg(right))


def dyadic_scale(value: DyadicInterval, scalar: Fraction | int) -> DyadicInterval:
    scalar = Fraction(scalar)
    if scalar < 0:
        return dyadic_neg(dyadic_scale(value, -scalar))
    return (
        _floor_ratio(value[0] * scalar.numerator, scalar.denominator),
        _ceil_ratio(value[1] * scalar.numerator, scalar.denominator),
    )


def dyadic_mul(left: DyadicInterval, right: DyadicInterval) -> DyadicInterval:
    products = (
        left[0] * right[0],
        left[0] * right[1],
        left[1] * right[0],
        left[1] * right[1],
    )
    return (
        _floor_ratio(min(products), DYADIC_SCALE),
        _ceil_ratio(max(products), DYADIC_SCALE),
    )


def dyadic_power(value: DyadicInterval, exponent: int) -> DyadicInterval:
    if exponent < 0:
        raise ValueError("dyadic_power only supports nonnegative exponents")
    if exponent == 0:
        return DYADIC_SCALE, DYADIC_SCALE
    denominator = DYADIC_SCALE ** (exponent - 1)
    if exponent % 2:
        return (
            _floor_ratio(value[0] ** exponent, denominator),
            _ceil_ratio(value[1] ** exponent, denominator),
        )
    if value[0] <= 0 <= value[1]:
        lower_numerator = 0
    else:
        lower_numerator = min(value[0] ** exponent, value[1] ** exponent)
    upper_numerator = max(value[0] ** exponent, value[1] ** exponent)
    return (
        _floor_ratio(lower_numerator, denominator),
        _ceil_ratio(upper_numerator, denominator),
    )


def atanh_log_interval(z: Fraction) -> DyadicInterval:
    """Outward enclosure of 2*atanh(z) for 0<=z<=1/3."""

    if not 0 <= z <= Fraction(1, 3):
        raise ValueError("atanh series range reduction failed")
    z_interval = dyadic_from_fraction(z)
    z_squared = dyadic_mul(z_interval, z_interval)
    power = z_interval
    total: DyadicInterval = (0, 0)
    for index in range(LOG_SERIES_TERMS):
        total = dyadic_add(total, dyadic_scale(power, Fraction(2, 2 * index + 1)))
        power = dyadic_mul(power, z_squared)
    # The omitted positive tail is at most
    # 2*z^(2m+1)/((2m+1)(1-z^2)) <= 9*z^(2m+1)/(4(2m+1)).
    tail = dyadic_scale(
        (0, power[1]), Fraction(9, 4 * (2 * LOG_SERIES_TERMS + 1))
    )
    return total[0], total[1] + tail[1]


LOG_TWO_INTERVAL = atanh_log_interval(Fraction(1, 3))


def log_fraction_interval(value: Fraction) -> DyadicInterval:
    """Exact dyadic enclosure of log(value), using power-of-two reduction."""

    if value <= 0:
        raise ValueError("logarithm requires a positive rational")
    reduced = value
    exponent = 0
    while reduced >= 2:
        reduced /= 2
        exponent += 1
    while reduced < 1:
        reduced *= 2
        exponent -= 1
    z = (reduced - 1) / (reduced + 1)
    return dyadic_add(atanh_log_interval(z), dyadic_scale(LOG_TWO_INTERVAL, exponent))


def centered_log_cumulants(
    eigenvalue_intervals: list[tuple[Fraction, Fraction]]
) -> dict[str, object]:
    log_intervals = [
        (
            log_fraction_interval(lower)[0],
            log_fraction_interval(upper)[1],
        )
        for lower, upper in eigenvalue_intervals
    ]
    dimension = len(log_intervals)
    mean = dyadic_scale(
        (
            sum(value[0] for value in log_intervals),
            sum(value[1] for value in log_intervals),
        ),
        Fraction(1, dimension),
    )
    centered = [
        (value[0] - mean[1], value[1] - mean[0]) for value in log_intervals
    ]
    moments: dict[int, DyadicInterval] = {}
    for order in (2, 4, 6, 8):
        powers = [dyadic_power(value, order) for value in centered]
        moments[order] = dyadic_scale(
            (sum(value[0] for value in powers), sum(value[1] for value in powers)),
            Fraction(1, dimension),
        )

    kappa2 = moments[2]
    kappa4 = dyadic_sub(moments[4], dyadic_scale(dyadic_power(moments[2], 2), 3))
    kappa6 = dyadic_add(
        dyadic_sub(
            moments[6],
            dyadic_scale(dyadic_mul(moments[4], moments[2]), 15),
        ),
        dyadic_scale(dyadic_power(moments[2], 3), 30),
    )
    kappa8 = dyadic_sub(
        dyadic_add(
            dyadic_sub(
                dyadic_sub(
                    moments[8],
                    dyadic_scale(dyadic_mul(moments[6], moments[2]), 28),
                ),
                dyadic_scale(dyadic_power(moments[4], 2), 35),
            ),
            dyadic_scale(
                dyadic_mul(moments[4], dyadic_power(moments[2], 2)), 420
            ),
        ),
        dyadic_scale(dyadic_power(moments[2], 4), 630),
    )
    cumulants = {2: kappa2, 4: kappa4, 6: kappa6, 8: kappa8}
    powers = {
        1: kappa2,
        2: dyadic_scale(kappa4, Fraction(-1, 2)),
        3: dyadic_scale(kappa6, Fraction(1, 16)),
        4: dyadic_scale(kappa8, Fraction(-1, 272)),
    }

    p1, p2, p3, p4 = (powers[index] for index in range(1, 5))
    inequalities = {
        "hankel_n_p2_minus_p1_squared": dyadic_sub(
            dyadic_scale(p2, dimension), dyadic_power(p1, 2)
        ),
        "hankel_p1_p3_minus_p2_squared": dyadic_sub(
            dyadic_mul(p1, p3), dyadic_power(p2, 2)
        ),
        "hankel_p2_p4_minus_p3_squared": dyadic_sub(
            dyadic_mul(p2, p4), dyadic_power(p3, 2)
        ),
        "newton_2e2": dyadic_sub(dyadic_power(p1, 2), p2),
        "newton_6e3": dyadic_add(
            dyadic_sub(
                dyadic_power(p1, 3),
                dyadic_scale(dyadic_mul(p1, p2), 3),
            ),
            dyadic_scale(p3, 2),
        ),
        "newton_24e4": dyadic_sub(
            dyadic_add(
                dyadic_add(
                    dyadic_sub(
                        dyadic_power(p1, 4),
                        dyadic_scale(dyadic_mul(dyadic_power(p1, 2), p2), 6),
                    ),
                    dyadic_scale(dyadic_power(p2, 2), 3),
                ),
                dyadic_scale(dyadic_mul(p1, p3), 8),
            ),
            dyadic_scale(p4, 6),
        ),
    }
    return {
        "dyadic_denominator_power": DYADIC_BITS,
        "log_series_terms": LOG_SERIES_TERMS,
        "mean_log_interval": interval_record(mean),
        "centered_moments": {
            str(order): interval_record(value) for order, value in moments.items()
        },
        "cumulants": {
            str(order): interval_record(value) for order, value in cumulants.items()
        },
        "normalized_power_sums": {
            str(order): interval_record(value) for order, value in powers.items()
        },
        "inequalities": {
            name: interval_record(value) for name, value in inequalities.items()
        },
    }


def interval_record(value: DyadicInterval) -> dict[str, object]:
    if value[0] > 0:
        sign = "strictly_positive"
    elif value[1] < 0:
        sign = "strictly_negative"
    elif value == (0, 0):
        sign = "zero"
    else:
        sign = "contains_zero"
    return {
        "lower_numerator": str(value[0]),
        "upper_numerator": str(value[1]),
        "denominator": f"2^{DYADIC_BITS}",
        "sign": sign,
        "lower_approx": value[0] / DYADIC_SCALE,
        "upper_approx": value[1] / DYADIC_SCALE,
    }


def record_interval_sign(record: dict[str, object]) -> str:
    return str(record["sign"])


def q_of_t(t: Fraction) -> Fraction:
    return (1 + t * t) / (2 * t)


def power_interval_positive(
    lower: Fraction, upper: Fraction, exponent: int
) -> tuple[Fraction, Fraction]:
    if not 0 < lower <= upper:
        raise ValueError("positive power interval required")
    if exponent >= 0:
        return lower**exponent, upper**exponent
    return upper**exponent, lower**exponent


def uniform_matrix_perturbation(
    point_matrix: list[list[Fraction]],
    sites: int,
    edges: tuple[Edge, ...],
    lower_t: Fraction,
    upper_t: Fraction,
) -> Fraction:
    """Exact infinity-norm bound on R(t)-R(1/3) over a rational t interval."""

    if not 0 < lower_t <= POINT_T <= upper_t < 1:
        raise ValueError("uniform interval must contain the physical point inside (0,1)")
    dimension = 1 << sites
    _, _, exponents = spin_bond_data(sites, edges)
    # q(t)=(t+t^{-1})/2 is decreasing on (0,1).
    q_lower = q_of_t(upper_t)
    q_upper = q_of_t(lower_t)
    t_lower_powers = [lower_t**power for power in range(2 * sites + 1)]
    t_upper_powers = [upper_t**power for power in range(2 * sites + 1)]
    q_intervals = {
        exponent: power_interval_positive(q_lower, q_upper, exponent)
        for exponent in set(exponents)
    }
    hamming = [
        [(row ^ state).bit_count() for state in range(dimension)]
        for row in range(dimension)
    ]
    row_sums = [Fraction(0) for _ in range(dimension)]
    for row in range(dimension):
        for column in range(row, dimension):
            lower = Fraction(0)
            upper = Fraction(0)
            for state in range(dimension):
                total_hamming = hamming[row][state] + hamming[column][state]
                q_low, q_high = q_intervals[exponents[state]]
                lower += t_lower_powers[total_hamming] * q_low
                upper += t_upper_powers[total_hamming] * q_high
            point = point_matrix[row][column]
            if not lower <= point <= upper:
                raise AssertionError("natural entry interval missed the point matrix")
            deviation = max(point - lower, upper - point)
            row_sums[row] += deviation
            if column != row:
                row_sums[column] += deviation
    return max(row_sums)


def uniform_eigenvalue_intervals(
    point_intervals: list[tuple[Fraction, Fraction]], perturbation: Fraction
) -> list[tuple[Fraction, Fraction]]:
    intervals = [
        (lower - perturbation, upper + perturbation)
        for lower, upper in point_intervals
    ]
    if any(lower <= 0 for lower, _ in intervals):
        raise AssertionError("uniform Weyl enclosure no longer proves positivity")
    return intervals


def fraction_record(value: Fraction) -> dict[str, object]:
    return {"exact": str(value), "approx": float(value)}


def exact_mode_control(label: str, squared_energies: tuple[int, ...]) -> dict[str, object]:
    powers = {
        order: sum((Fraction(value) ** order for value in squared_energies), Fraction(0))
        for order in range(1, 5)
    }
    p1, p2, p3, p4 = (powers[index] for index in range(1, 5))
    inequalities = {
        "hankel_n_p2_minus_p1_squared": len(squared_energies) * p2 - p1 * p1,
        "hankel_p1_p3_minus_p2_squared": p1 * p3 - p2 * p2,
        "hankel_p2_p4_minus_p3_squared": p2 * p4 - p3 * p3,
        "newton_2e2": p1 * p1 - p2,
        "newton_6e3": p1**3 - 3 * p1 * p2 + 2 * p3,
        "newton_24e4": p1**4 - 6 * p1 * p1 * p2 + 3 * p2 * p2 + 8 * p1 * p3 - 6 * p4,
    }
    return {
        "tag": "[COMPUTATION]",
        "label": label,
        "squared_mode_energies": list(squared_energies),
        "normalized_power_sums": {str(key): str(value) for key, value in powers.items()},
        "cumulants": {
            "2": str(p1),
            "4": str(-2 * p2),
            "6": str(16 * p3),
            "8": str(-272 * p4),
        },
        "inequalities": {key: str(value) for key, value in inequalities.items()},
    }


def physical_case(
    label: str,
    sites: int,
    edges: tuple[Edge, ...],
    role: str,
    started: float,
) -> tuple[dict[str, object], list[tuple[Fraction, Fraction]], list[list[Fraction]]]:
    colors = bipartite_coloring(sites, edges)
    matrix, denominator, integral, q, parity = build_rational_layer(sites, edges, POINT_T)
    budget_tick(started, f"{label} exact layer construction")
    spectral, eigenvalue_intervals = exact_congruence_certificate(integral, denominator)
    budget_tick(started, f"{label} exact congruence certificate")
    moments = centered_log_cumulants(eigenvalue_intervals)
    budget_tick(started, f"{label} exact log intervals")
    row = {
        "tag": "[COMPUTATION]",
        "label": label,
        "role": role,
        "sites": sites,
        "dimension": 1 << sites,
        "edges": [list(edge) for edge in edges],
        "bipartite_coloring": list(colors),
        "point": {
            "t_tanh_half_Kstar": str(POINT_T),
            "q_exp_2K": str(q),
            "bond_energy_parity": parity,
        },
        "exact_matrix": {
            "definition": "R=P_t diag(q^((b-epsilon)/2)) P_t",
            "common_denominator": str(denominator),
            "integral_matrix_sha256": matrix_digest(integral),
            "dimension": 1 << sites,
            "positive_definiteness_proof": (
                "P_t is the Kronecker power of [[1,t],[t,1]] and is invertible for 0<t<1; "
                "the middle diagonal is strictly positive, so R is symmetric positive definite"
            ),
        },
        "spectral_certificate": spectral,
        "log_moments": moments,
    }
    return row, eigenvalue_intervals, matrix


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    mode_controls = [
        exact_mode_control("unequal_nonnegative_modes", (0, 1, 4, 9)),
        exact_mode_control("equal_nonzero_modes_hankel_boundary", (4, 4, 4)),
    ]
    _check(
        checks,
        "exact Rademacher normalization controls",
        mode_controls[0]["cumulants"] == {
            "2": "14",
            "4": "-196",
            "6": "12704",
            "8": "-1854496",
        }
        and mode_controls[1]["inequalities"]["hankel_p1_p3_minus_p2_squared"] == "0",
        "unequal control fixes kappa_2,kappa_4,kappa_6,kappa_8; equal modes saturate the localizing Hankel minor",
    )

    specs = (
        ("open_chain_P4", 4, chain_edges(4), "free_fermion_control"),
        ("open_chain_P6", 6, chain_edges(6), "size_matched_free_fermion_control"),
        ("open_grid_2x2", 4, grid_edges(2, 2), "bipartite_obstruction"),
        ("open_grid_2x3", 6, grid_edges(2, 3), "bipartite_obstruction"),
    )
    cases: list[dict[str, object]] = []
    point_intervals: dict[str, list[tuple[Fraction, Fraction]]] = {}
    point_matrices: dict[str, list[list[Fraction]]] = {}
    for label, sites, edges, role in specs:
        row, intervals, matrix = physical_case(label, sites, edges, role, started)
        cases.append(row)
        point_intervals[label] = intervals
        point_matrices[label] = matrix
        certificate = row["spectral_certificate"]
        _check(
            checks,
            f"{label} exact full-spectrum enclosure",
            Fraction(str(certificate["gram_error_exact"])) < 1
            and len(certificate["ordered_eigenvalue_intervals"]) == 1 << sites
            and all(
                Fraction(item["lower"]) > 0
                for item in certificate["ordered_eigenvalue_intervals"]
            ),
            f"dimension={1 << sites}, f={certificate['gram_error_exact']}, degeneracies allowed",
        )
        _check(
            checks,
            f"{label} normalized cumulant signs and Newton necessities",
            all(
                record_interval_sign(row["log_moments"]["normalized_power_sums"][str(order)])
                == "strictly_positive"
                for order in range(1, 5)
            )
            and all(
                record_interval_sign(row["log_moments"]["inequalities"][name])
                == "strictly_positive"
                for name in ("newton_2e2", "newton_6e3", "newton_24e4")
            ),
            "p1..p4 and the first three Newton numerators have certified positive lower endpoints",
        )

    by_label = {str(row["label"]): row for row in cases}
    _check(
        checks,
        "physical case coverage is exact",
        set(by_label) == {label for label, _, _, _ in specs},
        "P4, P6, open 2x2, and open 2x3; no uncomputed 2x4 claim",
    )
    _check(
        checks,
        "all physical cases are bipartite",
        all(
            all(
                row["bipartite_coloring"][left] != row["bipartite_coloring"][right]
                for left, right in row["edges"]
            )
            for row in cases
        ),
        "stored colorings cross every stored edge",
    )
    for label in ("open_chain_P4", "open_chain_P6"):
        inequalities = by_label[label]["log_moments"]["inequalities"]
        _check(
            checks,
            f"{label} low-order route non-violation",
            record_interval_sign(inequalities["hankel_p1_p3_minus_p2_squared"])
            == "strictly_positive"
            and record_interval_sign(inequalities["hankel_p2_p4_minus_p3_squared"])
            == "strictly_positive",
            "both tested adjacent 2x2 Hankel minors have positive certified lower endpoints",
        )
    for label in ("open_grid_2x2", "open_grid_2x3"):
        inequalities = by_label[label]["log_moments"]["inequalities"]
        _check(
            checks,
            f"{label} strict point obstruction",
            record_interval_sign(inequalities["hankel_p1_p3_minus_p2_squared"])
            == "strictly_negative",
            "certified upper endpoint of p1*p3-p2^2 is below zero at t=1/3",
        )

    lower_t = POINT_T - INTERVAL_RADIUS
    upper_t = POINT_T + INTERVAL_RADIUS
    uniform_rows: list[dict[str, object]] = []
    for label, sites, edges, _ in specs[2:]:
        perturbation = uniform_matrix_perturbation(
            point_matrices[label], sites, edges, lower_t, upper_t
        )
        intervals = uniform_eigenvalue_intervals(point_intervals[label], perturbation)
        moments = centered_log_cumulants(intervals)
        h1 = moments["inequalities"]["hankel_p1_p3_minus_p2_squared"]
        row = {
            "tag": "[COMPUTATION]",
            "label": label,
            "t_interval": {"lower": str(lower_t), "upper": str(upper_t)},
            "entrywise_to_operator_norm_bound": fraction_record(perturbation),
            "method": (
                "positive rational natural intervals enclose every entry of R(t); symmetry gives "
                "||R(t)-R(1/3)||_2<=||.||_infinity, Weyl enlarges every ordered point interval, "
                "then exact dyadic log/cumulant interval arithmetic is repeated"
            ),
            "log_moments": moments,
            "strict_obstruction": record_interval_sign(h1) == "strictly_negative",
        }
        uniform_rows.append(row)
        _check(
            checks,
            f"{label} strict uniform-interval obstruction",
            bool(row["strict_obstruction"]),
            f"p1*p3-p2^2 remains below zero for every t in [{lower_t},{upper_t}]",
        )
        budget_tick(started, f"{label} uniform interval")

    if len({str(item["name"]) for item in checks}) != len(checks):
        raise AssertionError("check names must be unique")
    if not all(bool(item["passed"]) for item in checks):
        failed = [str(item["name"]) for item in checks if not bool(item["passed"])]
        raise AssertionError(f"producer checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e228_bipartite_log_moments.py",
            "interpreter": sys.executable,
            "arithmetic": (
                "exact Fraction layer matrices; float64 eigenvectors/eigenvalues are proposals only; "
                "exact integer congruence row-sum certificates plus generalized min-max enclose the "
                "entire spectrum; logarithms and cumulants use exact outward-rounded 2^-160 intervals"
            ),
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
            "cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "rss_cap_bytes": RSS_CAP_BYTES,
        },
        "data": {
            "necessary_conditions": {
                "tag": "[LEMMA]",
                "hypotheses": (
                    "N=2^n positive slots lambda_S=a product_{i in S}u_i with a,u_i>0; "
                    "the empirical log spectrum is determinant-centered before cumulants are taken"
                ),
                "rademacher_form": (
                    "log(lambda_S)-N^-1 sum_T log(lambda_T)="
                    "sum_i epsilon_i E_i, E_i=(log u_i)/2, epsilon_i in {-1,+1} uniformly"
                ),
                "normalization": {
                    "p1": "kappa_2=sum E_i^2",
                    "p2": "-kappa_4/2=sum E_i^4",
                    "p3": "kappa_6/16=sum E_i^6",
                    "p4": "-kappa_8/272=sum E_i^8",
                },
                "newton": (
                    "with y_i=E_i^2>=0, k e_k=sum_{j=1}^k(-1)^(j-1)e_(k-j)p_j; "
                    "thus every e_k>=0 and e_k=0 for k>n"
                ),
                "hankel": (
                    "[p_(i+j)] and [p_(i+j+1)] are positive semidefinite Gram matrices; "
                    "in particular n*p2-p1^2>=0, p1*p3-p2^2>=0, and p2*p4-p3^2>=0"
                ),
                "decisive_cumulant_form": (
                    "J1=16*(p1*p3-p2^2)=kappa_2*kappa_6-4*kappa_4^2>=0; "
                    "the adjacent form J2=544*(p2*p4-p3^2)=kappa_4*kappa_8-68*kappa_6^2>=0 "
                    "is redundant and not selected as the theorem test"
                ),
            },
            "exact_mode_controls": mode_controls,
            "physical_representative": {
                "tag": "[LEMMA]",
                "point": "t=1/3, q=(1+t^2)/(2t)=5/3",
                "benchmark_used_for_selection": False,
                "K_c_used_for_selection_or_tuning": False,
                "reciprocity": (
                    "for a bipartite graph, J=(product_all Z)(product_one_color X) sends both "
                    "A=sum X and B=sum ZZ to their negatives, so the determinant-one symmetric "
                    "physical transfer operator is conjugate to its inverse; R differs only by a "
                    "positive scalar, removed by determinant centering"
                ),
            },
            "physical_cases": cases,
            "uniform_interval_certificates": uniform_rows,
            "finite_result": {
                "tag": "[THEOREM]",
                "statement": (
                    "For every t in the displayed rational interval about 1/3, the positive "
                    "reciprocal spectra of the open 2x2 and 2x3 bipartite Ising layer operators "
                    "are not full n-mode subset-product spectra, because their certified "
                    "p1*p3-p2^2 upper endpoints are strictly negative."
                ),
            },
            "route_scope": {
                "tag": "[UNRESOLVED]",
                "statement": (
                    "This is a finite-layer, finite-coupling-interval obstruction.  The P4 and P6 "
                    "open-chain controls do not violate the tested conditions.  No all-coupling, "
                    "all-size, 2x4, trace-reciprocity, localization, varying-conjugator, or "
                    "pair-polynomial-sequence claim is made."
                ),
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for item in payload["checks"]:
        print(f"[PASS] {item['name']}: {item['detail']}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
