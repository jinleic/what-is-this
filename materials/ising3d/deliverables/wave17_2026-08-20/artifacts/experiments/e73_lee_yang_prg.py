"""Calibration-frozen Lee--Yang PRG and an exact three-size resolution barrier.

The experiment consumes the exact zero data frozen by e34.  Estimator choices are
made and hashed using only the two-dimensional control.  The unchanged protocol
is then offered the open-cube data.  Because only L=2,3,4 exist, the exponent is
not reported; instead an exact interpolation certificate proves that two
monotone analytic-correction models with different exponents fit every available
cube angle.

A fixed-field transfer recurrence also extends 2x2 bars without materializing a
field polynomial.  Those anisotropic data are explicitly diagnostic and are not
substituted for missing isotropic cubes.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Callable

import mpmath as mp
import sympy as sp

from ising.lee_yang import (
    lee_yang_roots,
    positive_zero_angles,
    rational_field_polynomial_transfer,
)
from ising.transfer_matrix import layer_bonds


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e73_lee_yang_prg.py"
INPUT = ROOT / "results" / "lee_yang" / "scaling.json"
OUTPUT = ROOT / "results" / "lee_yang" / "prg.json"
DPS = 120
RECHECK_DPS = 160
SERIAL_DIGITS = 90
T_RADIUS = Fraction(1, 10**70)
THETA_RADIUS = Fraction(1, 10**60)
COSINE_EVEN_ORDER = 80
CALIBRATION_CENTERS = tuple(range(5, 15))
MIN_SUFFIX_CENTERS = 5
BST_OMEGA = 1
BST_DEPTH = 2
TARGET_SIGMA_2D = Fraction(-1, 6)
BARRIER_EDGE = Fraction(1, 50)
BARRIER_POWERS = (2, 4)
HIGH_TEMPERATURE_IDS = ("x_7_over_10", "x_2_over_3")
TUBE_LENGTHS = (2, 4, 8, 16, 32, 64)
SCOUT_SEGMENTS = 1024

mp.mp.dps = DPS


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def mp_text(value: object, digits: int = SERIAL_DIGITS) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def fraction_decimal(value: Fraction, digits: int = 50) -> str:
    with mp.workdps(digits + 20):
        return mp_text(mp.mpf(value.numerator) / value.denominator, digits)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def divide_by_z_plus_one(coefficients: list[int]) -> list[int]:
    degree = len(coefficients) - 1
    if degree % 2 != 1 or coefficients != coefficients[::-1]:
        raise ValueError("division requires an odd-degree palindrome")
    quotient = [0] * degree
    quotient[0] = coefficients[0]
    for index in range(1, degree):
        quotient[index] = coefficients[index] - quotient[index - 1]
    reconstructed = (
        [quotient[0]]
        + [quotient[index - 1] + quotient[index] for index in range(1, degree)]
        + [quotient[-1]]
    )
    if reconstructed != coefficients or quotient != quotient[::-1]:
        raise ArithmeticError("exact z+1 division failed")
    return quotient


def reduced_power_coefficients(coefficients: list[int]) -> list[int]:
    """Return Q(t), ascending in t, for the palindromic field polynomial."""

    values = (
        divide_by_z_plus_one(coefficients)
        if (len(coefficients) - 1) % 2
        else list(coefficients)
    )
    half = (len(values) - 1) // 2
    power = [values[half]]
    if half == 0:
        return power
    previous = [1]
    current = [0, 1]
    for order in range(1, half + 1):
        if order == 1:
            chebyshev = current
        else:
            chebyshev = [0] * (len(current) + 1)
            for index, value in enumerate(current):
                chebyshev[index + 1] += 2 * value
            for index, value in enumerate(previous):
                chebyshev[index] -= value
            previous, current = current, chebyshev
        if len(power) < len(chebyshev):
            power.extend([0] * (len(chebyshev) - len(power)))
        multiplier = 2 * values[half - order]
        for index, value in enumerate(chebyshev):
            power[index] += multiplier * value
    return power


def evaluate_fraction_polynomial(coefficients: list[int], value: Fraction) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(coefficients):
        total = total * value + coefficient
    return total


def cosine_bounds(value: Fraction) -> tuple[Fraction, Fraction]:
    """Exact alternating-series bounds for cos(value), with 0 <= value < 2."""

    if not 0 <= value < 2:
        raise ValueError("cosine bound is specialized to 0 <= value < 2")
    term = Fraction(1)
    partial = term
    even_partial = None
    for order in range(1, COSINE_EVEN_ORDER + 2):
        term *= -value * value
        term /= (2 * order - 1) * (2 * order)
        partial += term
        if order == COSINE_EVEN_ORDER:
            even_partial = partial
    if even_partial is None or COSINE_EVEN_ORDER % 2:
        raise AssertionError("the configured upper partial sum must end at even order")
    odd_partial = partial
    if not odd_partial <= even_partial:
        raise ArithmeticError("alternating cosine bounds were reversed")
    return odd_partial, even_partial


def certify_theta_interval(size_record: dict) -> dict:
    """Enclose theta_1 exactly using rational polynomial and cosine signs."""

    coefficients = [int(value) for value in size_record["A_k"]]
    power = reduced_power_coefficients(coefficients)
    root_record = size_record["positive_zero_angles"][0]
    t_center = Fraction(root_record["t_refined"])
    t_lower = t_center - T_RADIUS
    t_upper = t_center + T_RADIUS
    value_lower = evaluate_fraction_polynomial(power, t_lower)
    value_upper = evaluate_fraction_polynomial(power, t_upper)
    sign_change = value_lower * value_upper < 0

    polynomial = sp.Poly.from_list(list(reversed(power)), gens=sp.Symbol("t"), domain=sp.ZZ)
    sym_lower = sp.Rational(t_lower.numerator, t_lower.denominator)
    sym_upper = sp.Rational(t_upper.numerator, t_upper.denominator)
    root_count = int(polynomial.count_roots(sym_lower, sym_upper))

    theta_center = Fraction(root_record["theta"])
    theta_lower = theta_center - THETA_RADIUS
    theta_upper = theta_center + THETA_RADIUS
    cosine_upper_lower, _ = cosine_bounds(theta_lower)
    _, cosine_lower_upper = cosine_bounds(theta_upper)
    # cos is decreasing: cos(theta_upper) < t_root < cos(theta_lower).
    cosine_enclosure = cosine_lower_upper < t_lower and t_upper < cosine_upper_lower
    passed = sign_change and root_count == 1 and cosine_enclosure
    return {
        "theta_center_decimal_rational": fraction_text(theta_center),
        "theta_interval": [fraction_text(theta_lower), fraction_text(theta_upper)],
        "theta_radius": fraction_text(THETA_RADIUS),
        "t_interval": [fraction_text(t_lower), fraction_text(t_upper)],
        "t_radius": fraction_text(T_RADIUS),
        "exact_reduced_polynomial_sign_change": sign_change,
        "exact_root_count_in_t_interval": root_count,
        "exact_cosine_series_enclosure": cosine_enclosure,
        "cosine_series_orders": [COSINE_EVEN_ORDER + 1, COSINE_EVEN_ORDER],
        "passed": passed,
    }


def ratio_function(power: mp.mpf, sizes: tuple[int, int, int]) -> mp.mpf:
    left, middle, right = (mp.mpf(value) for value in sizes)
    return (left ** (-power) - middle ** (-power)) / (
        middle ** (-power) - right ** (-power)
    )


def solve_prg_power(
    sizes: tuple[int, int, int], values: tuple[mp.mpf, mp.mpf, mp.mpf], dps: int
) -> dict:
    """Solve the three-size shift-ratio equation by guarded bisection."""

    with mp.workdps(dps):
        observed = (values[0] - values[1]) / (values[1] - values[2])
        lower = mp.mpf("1e-12")
        upper = mp.mpf(16)

        def residual(power: mp.mpf) -> mp.mpf:
            return ratio_function(power, sizes) - observed

        f_lower = residual(lower)
        f_upper = residual(upper)
        if not f_lower * f_upper < 0:
            raise ArithmeticError(f"PRG power not bracketed for sizes {sizes}")
        for _ in range(600):
            midpoint = (lower + upper) / 2
            f_midpoint = residual(midpoint)
            if f_lower * f_midpoint <= 0:
                upper = midpoint
                f_upper = f_midpoint
            else:
                lower = midpoint
                f_lower = f_midpoint
            if upper - lower < mp.power(10, -dps + 20):
                break
        power = (lower + upper) / 2
        amplitude = (values[0] - values[1]) / (
            mp.mpf(sizes[0]) ** (-power) - mp.mpf(sizes[1]) ** (-power)
        )
        edge = values[0] - amplitude * mp.mpf(sizes[0]) ** (-power)
        return {
            "sizes": list(sizes),
            "observed_drop_ratio": mp_text(observed),
            "power_p": mp_text(power),
            "power_bracket": [mp_text(lower), mp_text(upper)],
            "edge_from_three_point_model": mp_text(edge),
            "ratio_residual": mp_text(abs(residual(power))),
        }


def solve_linear_system(matrix: list[list[mp.mpf]], vector: list[mp.mpf]) -> list[mp.mpf]:
    result = mp.lu_solve(mp.matrix(matrix), mp.matrix(vector))
    return [result[index] for index in range(len(vector))]


def correction_fit(
    centers: tuple[int, ...], powers: dict[int, mp.mpf], dimension: int
) -> dict:
    """Overdetermined p(h)=p_inf+a*h+b*h^2 with h=1/L."""

    design = [[mp.mpf(1), mp.mpf(1) / center, mp.mpf(1) / center**2] for center in centers]
    normal = [
        [sum(row[i] * row[j] for row in design) for j in range(3)]
        for i in range(3)
    ]
    rhs = [sum(row[i] * powers[center] for row, center in zip(design, centers)) for i in range(3)]
    coefficients = solve_linear_system(normal, rhs)
    residuals = [
        sum(coefficient * row[index] for index, coefficient in enumerate(coefficients))
        - powers[center]
        for row, center in zip(design, centers)
    ]
    p_infinity = coefficients[0]
    sigma = mp.mpf(dimension) / p_infinity - 1
    return {
        "centers": list(centers),
        "p_infinity": mp_text(p_infinity),
        "sigma": mp_text(sigma),
        "coefficients_p_inf_a_b": [mp_text(value) for value in coefficients],
        "max_abs_p_residual": mp_text(max(abs(value) for value in residuals)),
    }


def bulirsch_stoer_table(centers: tuple[int, ...], values: list[mp.mpf]) -> list[list[mp.mpf]]:
    """Neville/Bulirsch--Stoer extrapolation at h=0 with h=(1/L)^omega."""

    rows = [list(values)]
    for depth in range(1, BST_DEPTH + 1):
        previous = rows[-1]
        current = []
        for index in range(len(previous) - 1):
            ratio = (mp.mpf(centers[index + depth]) / centers[index]) ** BST_OMEGA
            current.append(previous[index + 1] + (previous[index + 1] - previous[index]) / (ratio - 1))
        rows.append(current)
    return rows


def build_calibration(source: dict, dps: int) -> dict:
    zero_data = source["data"]["two_dimensional_calibration"]["zero_data"]
    with mp.workdps(dps):
        theta = {int(size): mp.mpf(record["theta_1"]) for size, record in zero_data.items()}
        prg = {}
        powers: dict[int, mp.mpf] = {}
        for center in CALIBRATION_CENTERS:
            result = solve_prg_power(
                (center - 1, center, center + 1),
                (theta[center - 1], theta[center], theta[center + 1]),
                dps,
            )
            power = mp.mpf(result["power_p"])
            powers[center] = power
            result["sigma_effective_d2"] = mp_text(mp.mpf(2) / power - 1)
            prg[str(center)] = result

        central = correction_fit(CALIBRATION_CENTERS, powers, dimension=2)
        leave_one_out = []
        for omitted in CALIBRATION_CENTERS:
            centers = tuple(center for center in CALIBRATION_CENTERS if center != omitted)
            fit = correction_fit(centers, powers, dimension=2)
            fit["omitted_center"] = omitted
            leave_one_out.append(fit)
        suffix_windows = []
        last_start = CALIBRATION_CENTERS[-1] - MIN_SUFFIX_CENTERS + 1
        for start in range(CALIBRATION_CENTERS[0], last_start + 1):
            centers = tuple(range(start, CALIBRATION_CENTERS[-1] + 1))
            suffix_windows.append(correction_fit(centers, powers, dimension=2))

        all_fits = [central, *leave_one_out, *suffix_windows]
        sigmas = [mp.mpf(fit["sigma"]) for fit in all_fits]
        uncertainty = [min(sigmas), max(sigmas)]
        target = mp.mpf(TARGET_SIGMA_2D.numerator) / TARGET_SIGMA_2D.denominator

        table = bulirsch_stoer_table(
            CALIBRATION_CENTERS, [powers[center] for center in CALIBRATION_CENTERS]
        )
        rolling_bst = []
        for index, value in enumerate(table[BST_DEPTH]):
            p_value = value
            rolling_bst.append(
                {
                    "centers": list(CALIBRATION_CENTERS[index : index + BST_DEPTH + 1]),
                    "p_infinity": mp_text(p_value),
                    "sigma_d2": mp_text(mp.mpf(2) / p_value - 1),
                }
            )

        passed = uncertainty[0] <= target <= uncertainty[1]
        return {
            "claim_tag": "EMPIRICAL_CALIBRATION",
            "input_sizes": list(range(4, 16)),
            "theta_1_sequence": [
                {"L": size, "theta_1": mp_text(theta[size])} for size in range(4, 16)
            ],
            "local_prg_definition": (
                "solve (theta_(L-1)-theta_L)/(theta_L-theta_(L+1))="
                "((L-1)^(-p)-L^(-p))/(L^(-p)-(L+1)^(-p)); sigma=d/p-1"
            ),
            "crossing_form": (
                "A_left(p)=(theta_(L-1)-theta_L)/((L-1)^(-p)-L^(-p)); "
                "A_right(p)=(theta_L-theta_(L+1))/(L^(-p)-(L+1)^(-p)); "
                "the PRG root is the amplitude crossing A_left(p)=A_right(p)"
            ),
            "local_prg": prg,
            "correction_structure": {
                "h": "1/L",
                "omega": BST_OMEGA,
                "depth": BST_DEPTH,
                "model": "p_L=p_infinity+a/L+b/L^2",
                "reason": "open boundaries admit analytic 1/L boundary corrections",
            },
            "bulirsch_stoer_depth_2_rolling": rolling_bst,
            "central_overdetermined_fit": central,
            "leave_one_center_out": leave_one_out,
            "largest_size_suffix_windows_minimum_five_centers": suffix_windows,
            "uncertainty_protocol": (
                "hull of the full fit, every single-center deletion, and every suffix ending "
                "at L=14 with at least five effective-exponent centers"
            ),
            "empirical_sigma_interval": [mp_text(uncertainty[0]), mp_text(uncertainty[1])],
            "target_sigma_exact": fraction_text(TARGET_SIGMA_2D),
            "target_inside_empirical_interval": passed,
            "passed": passed,
            "scope": (
                "The interval is a robust finite-window envelope, not a confidence interval "
                "or a rigorous truncation bound."
            ),
        }


def protocol_definition() -> dict:
    return {
        "selected_using": "two-dimensional control only",
        "three_dimensional_benchmark_used": False,
        "raw_theta_sizes_required_count": 12,
        "raw_theta_sizes": list(range(4, 16)),
        "effective_centers_required": 10,
        "effective_centers": list(CALIBRATION_CENTERS),
        "ratio_estimator": "three-consecutive-size power-shift PRG",
        "dimension_inserted_only_after_p_extrapolation": True,
        "correction_h": "1/L",
        "bulirsch_stoer_omega": BST_OMEGA,
        "bulirsch_stoer_depth": BST_DEPTH,
        "central_fit": "unweighted least-squares p_L=p_inf+a/L+b/L^2 over all centers",
        "uncertainty": (
            "hull(full, every leave-one-center-out, every >=5-center suffix ending at largest center)"
        ),
        "calibration_pass_condition": "exact -1/6 lies in empirical uncertainty hull",
        "on_insufficient_sizes": "return not_identifiable; do not shorten or refit protocol",
    }


def invert_fraction_matrix(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    work = [row[:] + [Fraction(int(i == j)) for j in range(size)] for i, row in enumerate(matrix)]
    for pivot_index in range(size):
        pivot_row = next(
            row for row in range(pivot_index, size) if work[row][pivot_index] != 0
        )
        work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
        pivot = work[pivot_index][pivot_index]
        work[pivot_index] = [value / pivot for value in work[pivot_index]]
        for row in range(size):
            if row == pivot_index:
                continue
            factor = work[row][pivot_index]
            if factor:
                work[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(work[row], work[pivot_index])
                ]
    return [row[size:] for row in work]


def adversarial_model(
    theta_centers: list[Fraction], power: int, dimension: int = 3
) -> dict:
    """Exact interval certificate for edge+L^-p(c0+c1/L+c2/L^2)."""

    sizes = (2, 3, 4)
    matrix = [
        [Fraction(1), Fraction(1, size), Fraction(1, size * size)] for size in sizes
    ]
    inverse = invert_fraction_matrix(matrix)
    rhs = [Fraction(size**power) * (theta - BARRIER_EDGE) for size, theta in zip(sizes, theta_centers)]
    coefficients = [
        sum((inverse[row][column] * rhs[column] for column in range(3)), Fraction(0))
        for row in range(3)
    ]
    error_radii = [
        THETA_RADIUS
        * sum(
            (abs(inverse[row][column]) * sizes[column] ** power for column in range(3)),
            Fraction(0),
        )
        for row in range(3)
    ]
    intervals = [
        (coefficient - radius, coefficient + radius)
        for coefficient, radius in zip(coefficients, error_radii)
    ]
    c0_lower, c0_upper = intervals[0]
    c1_lower, c1_upper = intervals[1]
    c2_lower, c2_upper = intervals[2]
    c1_square_upper = max(abs(c1_lower), abs(c1_upper)) ** 2
    correction_discriminant_upper = c1_square_upper - 4 * c0_lower * c2_lower
    derivative_discriminant_upper = (
        (power + 1) ** 2 * c1_square_upper
        - 4 * power * (power + 2) * c0_lower * c2_lower
    )
    residuals = []
    for size, theta in zip(sizes, theta_centers):
        correction = sum(
            coefficients[index] * Fraction(1, size**index) for index in range(3)
        )
        predicted = BARRIER_EDGE + Fraction(1, size**power) * correction
        residuals.append(predicted - theta)
    interval_certificate = (
        c0_lower > 0
        and c2_lower > 0
        and correction_discriminant_upper < 0
        and derivative_discriminant_upper < 0
    )
    sigma = Fraction(dimension, power) - 1
    return {
        "model": "theta(L)=1/50+L^(-p)*(c0+c1/L+c2/L^2)",
        "p": power,
        "sigma_exact_d3": fraction_text(sigma),
        "sigma_decimal_d3": fraction_decimal(sigma),
        "nominal_coefficients_exact": {
            name: fraction_text(value)
            for name, value in zip(("c0", "c1", "c2"), coefficients)
        },
        "nominal_coefficients_decimal": {
            name: fraction_decimal(value)
            for name, value in zip(("c0", "c1", "c2"), coefficients)
        },
        "coefficient_intervals_covering_exact_angles": {
            name: [fraction_text(bounds[0]), fraction_text(bounds[1])]
            for name, bounds in zip(("c0", "c1", "c2"), intervals)
        },
        "nominal_exact_fit_residuals": [fraction_text(value) for value in residuals],
        "correction_quadratic_discriminant_upper_exact": fraction_text(
            correction_discriminant_upper
        ),
        "derivative_quadratic_discriminant_upper_exact": fraction_text(
            derivative_discriminant_upper
        ),
        "all_coefficient_box_models_above_edge_and_strictly_decreasing_for_L_positive": interval_certificate,
        "proof": (
            "C(u)>0 because c2>0 and disc(C)<0. Also Q(L)=p*c0*L^2+"
            "(p+1)*c1*L+(p+2)*c2>0 because its leading coefficient is positive "
            "and disc(Q)<0; theta'(L)=-L^(-p-3)Q(L)<0."
        ),
    }


def build_barrier(source: dict, angle_certificates: dict[str, dict[str, dict]]) -> dict:
    cube_data = source["data"]["three_dimensional"]["exact_open_cube_zero_data"]
    couplings = {}
    for grid_id in HIGH_TEMPERATURE_IDS:
        record = cube_data[grid_id]
        theta_centers = [
            Fraction(item["theta_1"]) for item in record["theta_1_sequence"]
        ]
        models = [adversarial_model(theta_centers, power) for power in BARRIER_POWERS]
        couplings[grid_id] = {
            "x": record["x"],
            "cube_sizes": [2, 3, 4],
            "theta_angle_certificates": angle_certificates[grid_id],
            "common_edge_exact": fraction_text(BARRIER_EDGE),
            "models": models,
            "sigma_separation_exact": fraction_text(
                abs(Fraction(3, BARRIER_POWERS[0]) - Fraction(3, BARRIER_POWERS[1]))
            ),
            "both_models_certified": all(
                model[
                    "all_coefficient_box_models_above_edge_and_strictly_decreasing_for_L_positive"
                ]
                and model["nominal_exact_fit_residuals"] == ["0", "0", "0"]
                for model in models
            ),
        }
    return {
        "claim_tag": "THEOREM",
        "statement": (
            "For each certified-high-temperature L=2,3,4 sequence, both p=2 "
            "(sigma=1/2) and p=4 (sigma=-1/4) analytic-correction models share "
            "edge 1/50, fit all three exact angles, stay above the edge, and decrease "
            "strictly for every real L>0."
        ),
        "general_interpolation_lemma": (
            "For any p>0 and edge e, the three values uniquely determine C(u) of "
            "degree at most two through C(1/L_i)=L_i^p*(theta_i-e)."
        ),
        "angle_input_radius_exact": fraction_text(THETA_RADIUS),
        "couplings": couplings,
        "conclusion": (
            "The current cube sizes do not identify sigma without an independently "
            "proved bound on correction coefficients or additional isotropic sizes."
        ),
    }


def layer_geometry(periodic_cross: tuple[bool, bool]) -> tuple[list[int], list[int]]:
    bonds = layer_bonds((2, 2), periodic_cross)
    broken = []
    down = []
    for state in range(16):
        broken.append(sum(((state >> left) ^ (state >> right)) & 1 for left, right in bonds))
        down.append(4 - state.bit_count())
    return broken, down


def fixed_field_value_float(
    theta: float, length: int, periodic_cross: tuple[bool, bool]
) -> float:
    broken, down = layer_geometry(periodic_cross)
    x = 2.0 / 3.0
    z = cmath.exp(1j * theta)
    weights = [(x ** broken[state]) * (z ** down[state]) for state in range(16)]
    vector = weights[:]
    for _ in range(1, length):
        transformed = vector
        for bit in range(4):
            transformed = [
                transformed[state] + x * transformed[state ^ (1 << bit)]
                for state in range(16)
            ]
        vector = [weights[state] * transformed[state] for state in range(16)]
        scale = max(abs(value) for value in vector)
        vector = [value / scale for value in vector]
    phased = cmath.exp(-0.5j * 4 * length * theta) * sum(vector)
    return phased.real


def fixed_field_value_mp(
    theta: mp.mpf, length: int, periodic_cross: tuple[bool, bool], dps: int
) -> tuple[mp.mpf, mp.mpf]:
    with mp.workdps(dps):
        broken, down = layer_geometry(periodic_cross)
        x = mp.mpf(2) / 3
        z = mp.exp(1j * theta)
        weights = [x ** broken[state] * z ** down[state] for state in range(16)]
        vector = weights[:]
        for _ in range(1, length):
            transformed = vector
            for bit in range(4):
                transformed = [
                    transformed[state] + x * transformed[state ^ (1 << bit)]
                    for state in range(16)
                ]
            vector = [weights[state] * transformed[state] for state in range(16)]
            scale = max(abs(value) for value in vector)
            vector = [value / scale for value in vector]
        phased = mp.exp(-mp.j * 2 * length * theta) * sum(vector)
        return mp.re(phased), abs(mp.im(phased))


def scout_first_sign_change(
    length: int, periodic_cross: tuple[bool, bool], upper: float
) -> tuple[float, float]:
    lower = 0.0
    previous = fixed_field_value_float(lower, length, periodic_cross)
    for step in range(1, SCOUT_SEGMENTS + 1):
        current_theta = upper * step / SCOUT_SEGMENTS
        current = fixed_field_value_float(current_theta, length, periodic_cross)
        if previous * current < 0:
            return upper * (step - 1) / SCOUT_SEGMENTS, current_theta
        previous = current
    if upper < math.pi:
        return scout_first_sign_change(length, periodic_cross, math.pi)
    raise ArithmeticError(f"no fixed-field sign change found for length {length}")


def bisect_fixed_field_root(
    length: int,
    periodic_cross: tuple[bool, bool],
    bracket: tuple[float, float],
    dps: int,
    target_width: str,
) -> dict:
    with mp.workdps(dps):
        lower = mp.mpf(repr(bracket[0]))
        upper = mp.mpf(repr(bracket[1]))
        f_lower, leak_lower = fixed_field_value_mp(lower, length, periodic_cross, dps)
        f_upper, leak_upper = fixed_field_value_mp(upper, length, periodic_cross, dps)
        if not f_lower * f_upper < 0:
            raise ArithmeticError("multiprecision transfer signs do not bracket a root")
        threshold = mp.mpf(target_width)
        for _ in range(700):
            midpoint = (lower + upper) / 2
            f_midpoint, _ = fixed_field_value_mp(midpoint, length, periodic_cross, dps)
            if f_lower * f_midpoint <= 0:
                upper = midpoint
                f_upper = f_midpoint
            else:
                lower = midpoint
                f_lower = f_midpoint
            if upper - lower <= threshold:
                break
        if upper - lower > threshold:
            raise ArithmeticError("fixed-field bisection did not reach target width")
        midpoint = (lower + upper) / 2
        return {
            "theta_midpoint": mp_text(midpoint),
            "theta_sign_bracket": [mp_text(lower), mp_text(upper)],
            "bracket_width": mp_text(upper - lower),
            "endpoint_real_values": [mp_text(f_lower), mp_text(f_upper)],
            "endpoint_imaginary_leakage_before_bisection": mp_text(max(leak_lower, leak_upper)),
        }


def build_tube_data() -> dict:
    boundary_specs = {
        "open_2x2_bar": (False, False),
        "one_transverse_periodic_cylinder": (True, False),
    }
    sequences = {}
    validations = []
    for name, periodic_cross in boundary_specs.items():
        records = []
        previous_upper = math.pi
        for length in TUBE_LENGTHS:
            scout = scout_first_sign_change(length, periodic_cross, previous_upper)
            primary = bisect_fixed_field_root(
                length, periodic_cross, scout, DPS, "1e-65"
            )
            recheck = bisect_fixed_field_root(
                length, periodic_cross, scout, RECHECK_DPS, "1e-75"
            )
            with mp.workdps(RECHECK_DPS):
                difference = abs(
                    mp.mpf(primary["theta_midpoint"])
                    - mp.mpf(recheck["theta_midpoint"])
                )
            record = {
                "length": length,
                "shape": [2, 2, length],
                "theta_1": recheck["theta_midpoint"],
                "theta_sign_bracket_160_dps": recheck["theta_sign_bracket"],
                "bracket_width_160_dps": recheck["bracket_width"],
                "120_vs_160_dps_midpoint_difference": mp_text(difference),
                "scout_bracket_float64": [repr(scout[0]), repr(scout[1])],
                "empirical_sign_bracket": True,
            }
            records.append(record)
            previous_upper = float(mp.mpf(recheck["theta_midpoint"]))

        exact_polynomial = rational_field_polynomial_transfer(
            (2, 2, 2), 2, 3, periodic=(*periodic_cross, False)
        )
        exact_theta = positive_zero_angles(
            lee_yang_roots(exact_polynomial, dps=DPS)
        )[0]
        direct_theta = mp.mpf(records[0]["theta_1"])
        validations.append(
            {
                "boundary": name,
                "length": 2,
                "exact_small_polynomial_degree": len(exact_polynomial) - 1,
                "absolute_difference": mp_text(abs(exact_theta - direct_theta)),
                "passed_to_1e-60": abs(exact_theta - direct_theta) < mp.mpf("1e-60"),
            }
        )
        sequences[name] = records

    return {
        "claim_tag": "EMPIRICAL_DIAGNOSTIC",
        "coupling": "x=2/3",
        "method": (
            "evaluate z^(-N/2) Z(z) at fixed z=exp(i theta) by a 16-state layer "
            "recurrence; normalize each layer by a positive scalar; bracket the first "
            "real sign change; no field-polynomial coefficient array for lengths 4..64"
        ),
        "precision": {
            "primary_mpmath_dps": DPS,
            "independent_recheck_mpmath_dps": RECHECK_DPS,
            "target_160_dps_bracket_width": "1e-75",
            "float64_scout_subintervals": SCOUT_SEGMENTS,
            "classification": "empirical multiprecision brackets, not directed-rounding intervals",
        },
        "small_polynomial_validations": validations,
        "sequences": sequences,
        "geometric_warning": (
            "Increasing only the transfer length at fixed 2x2 cross-section is a "
            "quasi-one-dimensional limit. These records cannot be inserted into the "
            "frozen isotropic-cube protocol and do not add L=5 or larger cubes."
        ),
    }


def build_three_dimensional_application(
    source: dict, protocol_hash: str, calibration_passed: bool
) -> dict:
    cube_data = source["data"]["three_dimensional"]["exact_open_cube_zero_data"]
    diagnostics = {}
    for grid_id, record in cube_data.items():
        values = tuple(mp.mpf(item["theta_1"]) for item in record["theta_1_sequence"])
        result = solve_prg_power((2, 3, 4), values, DPS)
        with mp.workdps(DPS):
            p_value = mp.mpf(result["power_p"])
            sigma_diagnostic = mp.mpf(3) / p_value - 1
        diagnostics[grid_id] = {
            "x": record["x"],
            "p_three_point": result["power_p"],
            "sigma_if_naively_d3": mp_text(sigma_diagnostic),
            "edge_three_point": result["edge_from_three_point_model"],
            "tag": "UNCONTROLLED_DIAGNOSTIC_NOT_AN_ESTIMATE",
        }
    available_sizes = [2, 3, 4]
    required_sizes = list(range(4, 16))
    sufficient = all(size in available_sizes for size in required_sizes)
    return {
        "protocol_sha256": protocol_hash,
        "calibration_passed_before_application": calibration_passed,
        "available_isotropic_cube_sizes": available_sizes,
        "required_isotropic_cube_sizes": required_sizes,
        "required_consecutive_isotropic_size_count": len(required_sizes),
        "protocol_instantiable": sufficient,
        "edge_exponent_estimate": None,
        "status": (
            "NOT_IDENTIFIABLE: the frozen calibrated protocol needs isotropic cubes "
            "L=4,...,15 (10 local PRG centers), but only L=2,3,4 exist."
        ),
        "three_size_diagnostics": diagnostics,
        "external_3d_exponent_used_for_selection_or_comparison": None,
    }


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    source = json.loads(INPUT.read_text(encoding="utf-8"))
    source_valid = (
        source.get("provenance", {}).get("script") == "experiments/e34_lee_yang_scaling.py"
        and source.get("checks")
        and all(item.get("passed") for item in source["checks"])
    )
    check("frozen_e34_input_valid", bool(source_valid), f"sha256={sha256_file(INPUT)}")

    calibration = build_calibration(source, DPS)
    calibration_recheck = build_calibration(source, RECHECK_DPS)
    with mp.workdps(RECHECK_DPS):
        calibration_difference = abs(
            mp.mpf(calibration["central_overdetermined_fit"]["sigma"])
            - mp.mpf(calibration_recheck["central_overdetermined_fit"]["sigma"])
        )
    check(
        "two_dimensional_calibration_passes",
        calibration["passed"] and calibration["target_inside_empirical_interval"],
        f"interval={calibration['empirical_sigma_interval']}; target=-1/6",
    )
    check(
        "calibration_reproduces_at_160_dps",
        calibration_difference < mp.mpf("1e-80"),
        f"central sigma difference={mp_text(calibration_difference)}",
    )

    protocol = protocol_definition()
    protocol_hash = canonical_sha256(protocol)
    check(
        "protocol_selected_without_3d_benchmark",
        protocol["selected_using"] == "two-dimensional control only"
        and protocol["three_dimensional_benchmark_used"] is False,
        f"protocol_sha256={protocol_hash}",
    )

    cube_data = source["data"]["three_dimensional"]["exact_open_cube_zero_data"]
    angle_certificates: dict[str, dict[str, dict]] = {}
    for grid_id in HIGH_TEMPERATURE_IDS:
        per_size = {}
        for size in (2, 3, 4):
            per_size[str(size)] = certify_theta_interval(
                cube_data[grid_id]["sizes"][str(size)]
            )
        angle_certificates[grid_id] = per_size
    all_angles_certified = all(
        certificate["passed"]
        for per_size in angle_certificates.values()
        for certificate in per_size.values()
    )
    check(
        "exact_theta_intervals_certified",
        all_angles_certified,
        f"six angle radii={fraction_text(THETA_RADIUS)}",
    )

    application = build_three_dimensional_application(
        source, protocol_hash, calibration["passed"]
    )
    check(
        "frozen_protocol_not_shortened_for_3d",
        not application["protocol_instantiable"]
        and application["edge_exponent_estimate"] is None,
        application["status"],
    )

    barrier = build_barrier(source, angle_certificates)
    barrier_passed = all(
        record["both_models_certified"]
        for record in barrier["couplings"].values()
    )
    check(
        "adversarial_correction_models_exactly_certified",
        barrier_passed,
        "p=2 and p=4 models fit both certified-high-temperature sequences",
    )
    check(
        "adversarial_exponents_diverge",
        all(
            record["sigma_separation_exact"] == "3/4"
            for record in barrier["couplings"].values()
        ),
        "sigma=1/2 versus sigma=-1/4",
    )

    tube_data = build_tube_data()
    tube_validation = all(
        record["passed_to_1e-60"]
        for record in tube_data["small_polynomial_validations"]
    )
    tube_precision = all(
        mp.mpf(record["bracket_width_160_dps"]) <= mp.mpf("1e-75")
        and mp.mpf(record["120_vs_160_dps_midpoint_difference"]) < mp.mpf("1e-60")
        for sequence in tube_data["sequences"].values()
        for record in sequence
    )
    check(
        "fixed_field_transfer_matches_small_exact_polynomials",
        tube_validation,
        "both boundary choices agree at length 2 to 1e-60",
    )
    check(
        "fixed_field_transfer_precision_recheck",
        tube_precision,
        "all length-64 sequences have <=1e-75 brackets and 120/160-dps agreement",
    )

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "frozen_input": str(INPUT.relative_to(ROOT)),
            "frozen_input_sha256": sha256_file(INPUT),
            "exact_arithmetic": (
                "integer field polynomials; Fraction interpolation, root signs, cosine "
                "series bounds, coefficient boxes, and discriminant inequalities"
            ),
            "multiprecision": {
                "primary_mpmath_dps": DPS,
                "recheck_mpmath_dps": RECHECK_DPS,
                "role": "PRG solves and fixed-field transfer diagnostics only",
            },
            "three_dimensional_benchmark_exponent_used": False,
        },
        "data": {
            "definitions": {
                "theta_1": "smallest positive Lee--Yang zero angle",
                "shift_model": "theta_1(L)=theta_edge+A*L^(-p)*(1+corrections)",
                "sigma_relation": "sigma=d/p-1 for isotropic d-dimensional volume L^d",
            },
            "frozen_protocol": protocol,
            "frozen_protocol_sha256": protocol_hash,
            "two_dimensional_calibration": calibration,
            "three_dimensional_application": application,
            "exact_identifiability_barrier": barrier,
            "fixed_field_strip_cylinder_diagnostic": tube_data,
            "assessment": (
                "THEOREM-LIKE NEGATIVE RESULT: the 2D-selected PRG/BST protocol calibrates "
                "within its robust empirical window, but cannot be instantiated on three "
                "cube sizes; exact adversarial correction certificates prove sigma is not "
                "identified. No thermodynamic 3D exponent is reported."
            ),
        },
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    if not checks or not all(item["passed"] for item in checks):
        failed = [item["name"] for item in checks if not item["passed"]]
        raise SystemExit(f"embedded checks failed: {failed}")
    print("PASS")


if __name__ == "__main__":
    main()
