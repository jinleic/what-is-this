#!/usr/bin/env python3
"""Independent finite-difference check of Liu's mean-preserving transverse family.

This module takes ``uc/liu9_objective.py::_formula`` as the authoritative
objective.  In its notation, put

    a3 = 1 - a1 - a2,
    rho_0 = 1 - q,  rho_1 = q,
    z_0 = (b0, b2, b4),  z_1 = (b1, b3, b5),
    w_{c,i} = rho_c * a_i,
    prot(u,v) = u*v + u*(1-u)*v*(1-v).

With ``h`` the binary entropy in nats, the three expectations transcribed by
``_formula`` are exactly

    ehxy = sum_{c,i,d,j} w_{c,i} w_{d,j} h(z_{c,i} z_{d,j}),
    ehpi = sum_c rho_c sum_{i,j} a_i a_j h(prot(z_{c,i}, z_{c,j})),
    ehx  = sum_{c,i} w_{c,i} h(z_{c,i}).

Thus the objective is the quotient

    objective = numerator / ehx,
    numerator = (1-beta)*ehxy + beta*ehpi.

The gap convention used here (and by the nine-variable audit) is

    gap = numerator - ehx = ehx*(objective - 1),

not its negative.  The quantity differentiated below is
``gap - (4119063/33554432)*dist^2``.  Its epsilon derivative is obtained only
by evaluations of the source objective followed by a fourth-order forward
finite difference and Richardson extrapolation; no closed form for the
coefficient is used.  The existing hand-form coefficient is imported only
after all independent grid values have been computed.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

# This workstation is shared with live compute workers.  Set these before any
# import that can transitively load NumPy/SciPy or a threaded numerical library.
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
from liu9_objective import evaluate_mpmath  # noqa: E402


DPS = 120
KAPPA = Fraction(4_119_063, 33_554_432)
Q = Fraction(1, 2)
EPSILON_STEP_TEXT = "1e-14"
Y_DERIVATIVE_STEP_TEXT = "1e-4"
ZERO_REFERENCE_THRESHOLD_TEXT = "1e-80"


@dataclass(frozen=True)
class CoefficientRow:
    label: str
    y: mpmath.mpf
    coefficient: mpmath.mpf
    richardson_correction: mpmath.mpf
    numerical_uncertainty: mpmath.mpf
    sign: str


def _mp_fraction(value: Fraction) -> mpmath.mpf:
    return mpmath.mpf(value.numerator) / value.denominator


def family_values(
    parameters: MPParameters,
    y: mpmath.mpf,
    epsilon: mpmath.mpf,
    q: mpmath.mpf,
) -> tuple[mpmath.mpf, ...]:
    """The requested nine-vector, without using liu9_ninevar helpers."""
    zero = mpmath.mpf(0)
    return (
        parameters.p - epsilon * y / parameters.x,
        epsilon,
        q,
        parameters.x,
        y,
        zero,
        parameters.x,
        y,
        zero,
    )


def distance_squared(values: Sequence[Any], parameters: MPParameters) -> Any:
    """Direct transcription of the nine-variable distance definition."""
    if len(values) != 9:
        raise ValueError("expected Liu's nine coordinates")
    a1, a2, q, b0, b2, b4, b1, b3, b5 = values
    one = q * 0 + 1
    zero = one - one
    masses = (a1, a2, one - a1 - a2)

    def component(points: Sequence[Any]) -> Any:
        mean = sum(
            (mass * point for mass, point in zip(masses, points)),
            zero,
        )
        spread = sum(
            (
                mass
                * (point * (point - parameters.x))
                * (point * (point - parameters.x))
                for mass, point in zip(masses, points)
            ),
            zero,
        )
        return (mean - parameters.mean) * (mean - parameters.mean) + spread

    return ((one - q) * component((b0, b2, b4))
            + q * component((b1, b3, b5)))


def pencil_value(
    parameters: MPParameters,
    y: mpmath.mpf,
    epsilon: mpmath.mpf,
    q: mpmath.mpf,
) -> mpmath.mpf:
    """Evaluate gap-kappa*dist^2 through the authoritative objective."""
    values = family_values(parameters, y, epsilon, q)
    terms = evaluate_mpmath(values, parameters.beta, dps=DPS)
    gap = terms.numerator - terms.ehx
    return +(gap - _mp_fraction(KAPPA) * distance_squared(values, parameters))


def fourth_order_forward_difference(
    parameters: MPParameters,
    y: mpmath.mpf,
    step: mpmath.mpf,
    q: mpmath.mpf,
) -> mpmath.mpf:
    """Five-point, one-sided derivative using only epsilon >= 0."""
    values = [
        pencil_value(parameters, y, index * step, q)
        for index in range(5)
    ]
    return (
        -25 * values[0]
        + 48 * values[1]
        - 36 * values[2]
        + 16 * values[3]
        - 3 * values[4]
    ) / (12 * step)


def finite_difference_coefficient(
    parameters: MPParameters,
    y: mpmath.mpf,
    q: mpmath.mpf,
) -> tuple[mpmath.mpf, mpmath.mpf, mpmath.mpf]:
    """Return the Richardson value, correction, and roundoff-aware uncertainty."""
    step = mpmath.mpf(EPSILON_STEP_TEXT)
    coarse = fourth_order_forward_difference(parameters, y, step, q)
    fine_step = step / 2
    fine = fourth_order_forward_difference(parameters, y, fine_step, q)
    extrapolated = fine + (fine - coarse) / 15
    correction = abs(extrapolated - fine)

    # A forward stencil subtracts O(1) values and divides by fine_step.  This
    # guard is deliberately conservative relative to mpmath's working epsilon;
    # the Richardson correction alone does not measure shared cancellation.
    roundoff_allowance = 512 * mpmath.eps / fine_step
    uncertainty = correction + roundoff_allowance
    return +extrapolated, +correction, +uncertainty


def exact_mean_preservation_check() -> tuple[Fraction, Fraction, Fraction]:
    """Check the epsilon coefficient and a complete witness using Fractions."""
    p = Fraction(17, 19)
    x = Fraction(7, 11)
    y = Fraction(5, 13)
    epsilon = Fraction(1, 1000)

    epsilon_coefficient = (-y / x) * x + y
    original_mean = p * x
    perturbed_mean = (p - epsilon * y / x) * x + epsilon * y
    if epsilon_coefficient != 0 or perturbed_mean != original_mean:
        raise AssertionError("the exact mean-preservation identity failed")
    return epsilon_coefficient, original_mean, perturbed_mean


def grid_points(parameters: MPParameters) -> list[tuple[str, mpmath.mpf]]:
    """A 32-point grid with endpoint and transverse clustering."""
    decimal_points = (
        "1e-8",
        "1e-6",
        "1e-4",
        "0.001",
        "0.005",
        "0.01",
        "0.03",
        "0.05",
        "0.1",
        "0.2",
        "0.3",
        "0.4",
        "0.5",
        "0.6",
    )
    points = [(text, mpmath.mpf(text)) for text in decimal_points]
    points.extend((
        ("x-1e-4", parameters.x - mpmath.mpf("1e-4")),
        ("x-1e-8", parameters.x - mpmath.mpf("1e-8")),
        ("x-1e-10", parameters.x - mpmath.mpf("1e-10")),
        ("x", parameters.x),
        ("x+1e-10", parameters.x + mpmath.mpf("1e-10")),
        ("x+1e-8", parameters.x + mpmath.mpf("1e-8")),
        ("x+1e-4", parameters.x + mpmath.mpf("1e-4")),
    ))
    points.extend(
        (text, mpmath.mpf(text))
        for text in (
            "0.72",
            "0.75",
            "0.8",
            "0.85",
            "0.9",
            "0.95",
            "0.99",
            "0.999",
            "0.9999",
            "0.999999",
            "0.99999999",
        )
    )
    if len(points) < 25:
        raise AssertionError("the requested grid needs at least 25 points")
    if any(not (0 < y < 1) for _, y in points):
        raise AssertionError("every grid point must be strictly between 0 and 1")
    if any(points[index][1] >= points[index + 1][1]
           for index in range(len(points) - 1)):
        raise AssertionError("the reporting grid must be strictly increasing")
    return points


def classify_sign(
    coefficient: mpmath.mpf,
    uncertainty: mpmath.mpf,
) -> str:
    if coefficient > uncertainty:
        return "POSITIVE"
    if coefficient < -uncertainty:
        return "NEGATIVE"
    return "ZERO_WITHIN_NUMERICAL_ERROR"


def compute_independent_grid(
    parameters: MPParameters,
    q: mpmath.mpf,
) -> list[CoefficientRow]:
    rows = []
    for label, y in grid_points(parameters):
        coefficient, correction, uncertainty = finite_difference_coefficient(
            parameters, y, q
        )
        rows.append(CoefficientRow(
            label=label,
            y=y,
            coefficient=coefficient,
            richardson_correction=correction,
            numerical_uncertainty=uncertainty,
            sign=classify_sign(coefficient, uncertainty),
        ))
    return rows


def five_point_y_derivatives(
    parameters: MPParameters,
    q: mpmath.mpf,
    step: mpmath.mpf,
) -> tuple[mpmath.mpf, mpmath.mpf]:
    """Fourth-order centered estimates of C'(x) and C''(x)."""
    values = [
        finite_difference_coefficient(parameters, parameters.x + offset * step, q)[0]
        for offset in (-2, -1, 0, 1, 2)
    ]
    first = (values[0] - 8 * values[1] + 8 * values[3] - values[4]) / (
        12 * step
    )
    second = (
        -values[0]
        + 16 * values[1]
        - 30 * values[2]
        + 16 * values[3]
        - values[4]
    ) / (12 * step * step)
    return first, second


def coefficient_y_derivatives(
    parameters: MPParameters,
    q: mpmath.mpf,
) -> tuple[mpmath.mpf, mpmath.mpf, mpmath.mpf, mpmath.mpf]:
    """Two Richardson levels on the centered y-derivative stencils."""
    step = mpmath.mpf(Y_DERIVATIVE_STEP_TEXT)
    coarse = five_point_y_derivatives(parameters, q, step)
    medium = five_point_y_derivatives(parameters, q, step / 2)
    fine = five_point_y_derivatives(parameters, q, step / 4)

    richardson_coarse = tuple(
        medium[index] + (medium[index] - coarse[index]) / 15
        for index in range(2)
    )
    richardson_fine = tuple(
        fine[index] + (fine[index] - medium[index]) / 15
        for index in range(2)
    )
    extrapolated = tuple(
        richardson_fine[index]
        + (richardson_fine[index] - richardson_coarse[index]) / 63
        for index in range(2)
    )
    corrections = tuple(
        abs(extrapolated[index] - richardson_fine[index])
        for index in range(2)
    )
    return (
        +extrapolated[0],
        +extrapolated[1],
        +corrections[0],
        +corrections[1],
    )


def compare_with_existing(
    parameters: MPParameters,
    rows: Sequence[CoefficientRow],
) -> dict[str, Any]:
    """Import and compare the existing hand form only after rows already exist."""
    from liu9_ninevar import transverse_coefficients_mp  # noqa: PLC0415

    zero_threshold = mpmath.mpf(ZERO_REFERENCE_THRESHOLD_TEXT)
    comparisons = []
    zero_comparisons = []
    for row in rows:
        reference = transverse_coefficients_mp(
            parameters, row.y, "mean_preserving"
        )[0]
        absolute = abs(row.coefficient - reference)
        if abs(reference) <= zero_threshold:
            zero_comparisons.append((absolute, row, reference))
            continue
        relative = absolute / abs(reference)
        comparisons.append((relative, absolute, row, reference))

    if not comparisons:
        raise AssertionError("there were no nonzero coefficients to compare")
    worst = max(comparisons, key=lambda item: item[0])
    maximum_zero_absolute = max(
        zero_comparisons,
        key=lambda item: item[0],
        default=(mpmath.mpf(0), None, mpmath.mpf(0)),
    )
    return {
        "worst_relative": worst[0],
        "worst_absolute": worst[1],
        "worst_row": worst[2],
        "worst_reference": worst[3],
        "maximum_zero_absolute": maximum_zero_absolute[0],
        "maximum_zero_row": maximum_zero_absolute[1],
        "maximum_zero_reference": maximum_zero_absolute[2],
    }


def _fmt(value: Any, digits: int = 24) -> str:
    return mpmath.nstr(value, digits)


def main() -> int:
    mpmath.mp.dps = DPS
    parameters = solve_equation_parameters(DPS)
    q = _mp_fraction(Q)

    print("INDEPENDENT LIU TRANSVERSE FINITE-DIFFERENCE CHECK")
    print("=" * 59)
    print("\n[1] OBJECTIVE AND GAP CONVENTION")
    print("objective = ((1-beta)*ehxy + beta*ehpi) / ehx")
    print("gap       = ((1-beta)*ehxy + beta*ehpi) - ehx")
    print("           = ehx*(objective-1), so the sign is not reversed")

    coefficient, original_mean, perturbed_mean = exact_mean_preservation_check()
    print("\n[2] EXACT MEAN-PRESERVATION CHECK")
    print("symbolic: (p-epsilon*y/x)*x + epsilon*y = p*x")
    print("epsilon coefficient (-y/x)*x+y =", coefficient)
    print("exact Fraction witness: original mean =", original_mean)
    print("exact Fraction witness: perturbed mean =", perturbed_mean)
    print("mean preservation: PASS")

    print("\n[3] NUMERICAL SETUP")
    print("mpmath decimal precision:", DPS)
    print("x    =", _fmt(parameters.x, 52))
    print("p    =", _fmt(parameters.p, 52))
    print("beta =", _fmt(parameters.beta, 52))
    print("kappa = 4119063/33554432 =", _fmt(_mp_fraction(KAPPA), 30))
    print("q = 1/2 (both component support triples coincide, so q is a gauge)")
    print("epsilon steps:", EPSILON_STEP_TEXT, "and", _fmt(mpmath.mpf(EPSILON_STEP_TEXT) / 2))
    print("method: five-point fourth-order forward difference + Richardson")

    # This complete table is deliberately computed before liu9_ninevar is
    # imported by compare_with_existing.
    rows = compute_independent_grid(parameters, q)
    print("\n[4] INDEPENDENT COEFFICIENT GRID")
    print("positive? is YES only when coefficient > numerical uncertainty")
    print(f"{'label':>12}  {'y':>25}  {'epsilon coefficient':>29}  {'positive?':>17}")
    for row in rows:
        positive = "YES" if row.sign == "POSITIVE" else "NO (numerical zero)"
        if row.sign == "NEGATIVE":
            positive = "NO (NEGATIVE)"
        print(
            f"{row.label:>12}  {_fmt(row.y, 19):>25}  "
            f"{_fmt(row.coefficient, 22):>29}  {positive:>17}"
        )

    minimum = min(rows, key=lambda row: row.coefficient)
    positive_rows = [row for row in rows if row.sign == "POSITIVE"]
    zero_rows = [row for row in rows if row.sign == "ZERO_WITHIN_NUMERICAL_ERROR"]
    negative_rows = [row for row in rows if row.sign == "NEGATIVE"]
    print("\ngrid size:", len(rows))
    print(
        "sign counts: positive =", len(positive_rows),
        ", numerical zero =", len(zero_rows),
        ", negative =", len(negative_rows),
    )
    print("grid minimum label:", minimum.label)
    print("grid minimum y:", _fmt(minimum.y, 52))
    print("grid minimum raw coefficient:", _fmt(minimum.coefficient, 30))
    print("grid minimum uncertainty:", _fmt(minimum.numerical_uncertainty, 12))
    print(
        "maximum epsilon-Richardson correction:",
        _fmt(max(row.richardson_correction for row in rows), 12),
    )
    if negative_rows:
        raise AssertionError("a negative independent grid coefficient was found")
    if len(zero_rows) != 1 or zero_rows[0].label != "x":
        raise AssertionError("only y=x should be numerically zero on this grid")

    first, second, first_correction, second_correction = (
        coefficient_y_derivatives(parameters, q)
    )
    print("\n[5] y-DERIVATIVES OF THE COEFFICIENT AT y=x")
    print("centered five-point stencils with two Richardson levels")
    print("base y step:", Y_DERIVATIVE_STEP_TEXT)
    print("C'(x)  =", _fmt(first, 32))
    print("C'(x) Richardson correction =", _fmt(first_correction, 12))
    print("C''(x) =", _fmt(second, 32))
    print("C''(x) Richardson correction =", _fmt(second_correction, 12))
    if abs(first) <= first_correction:
        print("first derivative verdict: numerically zero")
    else:
        print("first derivative verdict: resolved nonzero")

    # Import and compare with liu9_ninevar only after every independent number,
    # including the y-derivative estimates, has been computed.
    comparison = compare_with_existing(parameters, rows)
    worst_row = comparison["worst_row"]
    print("\n[6] FINAL COMPARISON WITH liu9_ninevar")
    print("relative disagreement = |finite_difference-reference|/|reference|")
    print("the reference-zero row y=x is excluded from that relative ratio")
    print("maximum relative disagreement:", _fmt(comparison["worst_relative"], 18))
    print("worst y label:", worst_row.label)
    print("worst y:", _fmt(worst_row.y, 52))
    print("worst absolute disagreement:", _fmt(comparison["worst_absolute"], 18))
    print("finite-difference uncertainty there:", _fmt(worst_row.numerical_uncertainty, 18))
    print("maximum absolute disagreement on reference-zero rows:",
          _fmt(comparison["maximum_zero_absolute"], 18))
    if comparison["worst_absolute"] > worst_row.numerical_uncertainty:
        print("*** DISAGREEMENT EXCEEDS FINITE-DIFFERENCE UNCERTAINTY ***")
        print("*** WORST Y:", _fmt(worst_row.y, 52), "***")
        raise AssertionError("existing coefficient disagrees with the independent check")
    print("comparison verdict: AGREES within finite-difference uncertainty")

    print("\nFINAL VERDICT: all nontrivial grid coefficients are positive;")
    print("y=x is zero within numerical error; the existing formula agrees.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
