"""Exact and numerical analysis of Liu's residual kernel R.

PROVED in this file: for 0 <= s,t <= 1,

    R(s,t) = (1-s)(1-t) F(st) - G((1-s)(1-t)(1+st))

is a negative-semidefinite kernel, where

    F(x) = x - (1+x) log(1+x),
    G(x) = x + (1-x) log(1-x).

PROVED consequence (using the exact reduction recorded in ``liu_kernel.py``):
this settles Liu's Hypothesis 1 in arXiv:2306.08824v1, Section V-A.
It does NOT settle Liu's Hypothesis 2 in Section V-B (the nine-parameter
optimisation), so Liu's constant 0.382709087918741 is not yet unconditional.

The proof below is an explicit Gram decomposition obtained from a
Lorentz-signature factorisation of a degree-three denominator.  NUMERICAL
Mercer diagnostics and PROVED finite Arb compressions are also reported.  A
status label on an output heading applies to every indented line below it until
the next status label.

Run from the repository root with

    math/.venv/bin/python math/uc/liu_R_nsd.py

The quadrature order defaults to 300 and may be changed within the requested
200--400 range with LIU_R_QUAD.  BLAS is constrained to one thread before NumPy
is imported.
"""

from __future__ import annotations

import os

# This diagnostic shares a workstation with live certification workers.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

from functools import lru_cache
from fractions import Fraction
from math import comb

import numpy as np
import sympy as sp
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_jacobi, eval_legendre

QUADRATURE_ORDER = int(os.environ.get("LIU_R_QUAD", "300"))
if not 200 <= QUADRATURE_ORDER <= 400:
    raise ValueError("LIU_R_QUAD must lie in the requested range 200..400")

FINGERPRINT_MODES = 5
CANDIDATE_DEGREE = 16
FIT_DEGREE = 12
THETA_QUADRATURE_ORDER = 64
ARB_PRECISION_BITS = 256
ARB_G_CUTOFF = 80


def f_np(x: np.ndarray) -> np.ndarray:
    """Return F(x) on the quadrature interior."""
    return x - (1.0 + x) * np.log1p(x)


def g_np(x: np.ndarray) -> np.ndarray:
    """Return G(x), including its continuous value G(1)=1."""
    x = np.asarray(x)
    out = np.array(x, copy=True)
    interior = x < 1.0
    out[interior] = (
        x[interior] + (1.0 - x[interior]) * np.log1p(-x[interior])
    )
    return out


def kernel_arrays(s: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return A, B, z, R and -R on a tensor product of the nodes."""
    left = s[:, None]
    right = s[None, :]
    a = (1.0 - left) * (1.0 - right)
    b = left * right
    z = a * (1.0 + b)
    r_kernel = a * f_np(b) - g_np(z)
    return a, b, z, r_kernel, -r_kernel


def weighted_correlation(
    first: np.ndarray, second: np.ndarray, weights: np.ndarray
) -> float:
    numerator = abs(float(np.dot(weights * first, second)))
    denominator = np.sqrt(
        float(np.dot(weights * first, first))
        * float(np.dot(weights * second, second))
    )
    return numerator / denominator


def weighted_span_residual(
    function: np.ndarray, design: np.ndarray, weights: np.ndarray
) -> float:
    root_weight = np.sqrt(weights)
    coefficients, *_ = np.linalg.lstsq(
        root_weight[:, None] * design, root_weight * function, rcond=None
    )
    residual = function - design @ coefficients
    return float(np.sqrt(np.dot(weights * residual, residual)))


def best_atom(
    function: np.ndarray,
    candidates: list[tuple[str, np.ndarray]],
    weights: np.ndarray,
) -> tuple[str, float]:
    correlation, name = max(
        (weighted_correlation(function, atom, weights), name)
        for name, atom in candidates
    )
    return name, correlation


def regression_r_squared(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    slope, intercept = np.polyfit(x, y, 1)
    prediction = slope * x + intercept
    residual = float(np.dot(y - prediction, y - prediction))
    centered = y - y.mean()
    total = float(np.dot(centered, centered))
    return float(slope), 1.0 - residual / total


def mercer_fingerprint() -> dict[str, object]:
    """Compute the requested Nystrom spectrum and candidate fits."""
    raw_nodes, raw_weights = leggauss(QUADRATURE_ORDER)
    s = (raw_nodes + 1.0) / 2.0
    weights = raw_weights / 2.0
    root_weight = np.sqrt(weights)
    a, b, z, r_kernel, minus_r = kernel_arrays(s)
    nystrom = root_weight[:, None] * minus_r * root_weight[None, :]
    nystrom = (nystrom + nystrom.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(nystrom)
    descending = eigenvalues[::-1]

    modes = np.arange(3, 21, dtype=float)
    logarithms = np.log(descending[2:20])
    root_slope, root_r2 = regression_r_squared(np.sqrt(modes), logarithms)
    geometric_slope, geometric_r2 = regression_r_squared(modes, logarithms)
    power_slope, power_r2 = regression_r_squared(np.log(modes), logarithms)

    u = 1.0 - s
    monomials = [(f"s^{degree}", s**degree) for degree in range(CANDIDATE_DEGREE + 1)]
    endpoint_polynomials = [
        (f"(1-s)s^{degree}", u * s**degree)
        for degree in range(CANDIDATE_DEGREE + 1)
    ]
    legendre_atoms = [
        (f"P_{degree}(2s-1)", eval_legendre(degree, raw_nodes))
        for degree in range(CANDIDATE_DEGREE + 1)
    ]
    jacobi_atoms: list[tuple[str, np.ndarray]] = []
    for alpha, beta in ((0, 1), (1, 0), (1, 1), (0, 2), (2, 0)):
        for degree in range(CANDIDATE_DEGREE + 1):
            jacobi_atoms.append(
                (
                    f"Jacobi({alpha},{beta})_{degree}",
                    eval_jacobi(degree, alpha, beta, raw_nodes),
                )
            )
    log_atoms = [
        ("log(s)", np.log(s)),
        ("log(1-s)", np.log(u)),
        ("s log(s)", s * np.log(s)),
        ("(1-s) log(1-s)", u * np.log(u)),
        ("log(1+s)", np.log1p(s)),
        ("(1-s) log(1+s)", u * np.log1p(s)),
    ]

    polynomial_span = np.column_stack([s**degree for degree in range(FIT_DEGREE + 1)])
    endpoint_span = np.column_stack(
        [u * s**degree for degree in range(FIT_DEGREE + 1)]
    )
    log_span = np.column_stack(
        [
            np.ones_like(s),
            u,
            np.log(s),
            np.log(u),
            s * np.log(s),
            u * np.log(u),
            np.log1p(s),
            u * np.log1p(s),
        ]
    )

    fits: list[dict[str, object]] = []
    for mode in range(FINGERPRINT_MODES):
        function = eigenvectors[:, -1 - mode] / root_weight
        function /= np.sqrt(np.dot(weights * function, function))
        largest = int(np.argmax(np.abs(function)))
        if function[largest] < 0.0:
            function = -function
        fits.append(
            {
                "mode": mode + 1,
                "monomial": best_atom(function, monomials, weights),
                "endpoint": best_atom(function, endpoint_polynomials, weights),
                "legendre": best_atom(function, legendre_atoms, weights),
                "jacobi": best_atom(function, jacobi_atoms, weights),
                "log": best_atom(function, log_atoms, weights),
                "poly_residual": weighted_span_residual(
                    function, polynomial_span, weights
                ),
                "endpoint_residual": weighted_span_residual(
                    function, endpoint_span, weights
                ),
                "log_residual": weighted_span_residual(function, log_span, weights),
            }
        )

    return {
        "s": s,
        "weights": weights,
        "a": a,
        "b": b,
        "z": z,
        "r": r_kernel,
        "minus_r": minus_r,
        "eigenvalues": eigenvalues,
        "descending": descending,
        "fits": fits,
        "root_slope": root_slope,
        "root_r2": root_r2,
        "geometric_slope": geometric_slope,
        "geometric_r2": geometric_r2,
        "power_slope": power_slope,
        "power_r2": power_r2,
    }


def series_coefficient(p: int, q: int) -> Fraction:
    """Exact coefficient c_(p,q) of A^p B^q in R."""
    if p == 1 and q >= 2:
        return Fraction(-1 if q % 2 == 0 else 1, q * (q - 1))
    if p >= 2 and 0 <= q <= p:
        return -Fraction(comb(p, q), p * (p - 1))
    return Fraction(0)


def exact_series_checks() -> None:
    """Check the all-index coefficient formula on exact symbolic truncations."""
    x, aa, bb = sp.symbols("x aa bb")
    f_series = sp.series(x - (1 + x) * sp.log(1 + x), x, 0, 14).removeO()
    for q in range(2, 13):
        fraction = series_coefficient(1, q)
        expected = sp.Rational(fraction.numerator, fraction.denominator)
        assert sp.expand(f_series).coeff(x, q) == expected

    g_truncation = sum(
        aa**p * (1 + bb) ** p / sp.Integer(p * (p - 1))
        for p in range(2, 9)
    )
    r_truncation = aa * f_series.subs(x, bb) - g_truncation
    for p in range(1, 9):
        for q in range(0, 13):
            fraction = series_coefficient(p, q)
            expected = sp.Rational(fraction.numerator, fraction.denominator)
            assert sp.expand(r_truncation).coeff(aa, p).coeff(bb, q) == expected


def reduced_feature_vector(p: int, q: int, degree: int) -> sp.Matrix:
    """Coefficients of (1-s)^(p-1) s^q after the common (1-s) factor."""
    return sp.Matrix(
        [
            sp.Integer((-1) ** (index - q) * comb(p - 1, index - q))
            if q <= index <= q + p - 1
            else sp.Integer(0)
            for index in range(degree + 1)
        ]
    )


def reduced_sos_block(max_f_degree: int, max_g_power: int) -> sp.Matrix:
    """Exact coefficient matrix of a finite block of -R after factoring uv."""
    degree = max(max_f_degree, 2 * max_g_power - 1)
    matrix = sp.zeros(degree + 1)
    for q in range(2, max_f_degree + 1):
        vector = sp.zeros(degree + 1, 1)
        vector[q] = 1
        matrix += sp.Rational((-1) ** q, q * (q - 1)) * vector * vector.T
    for p in range(2, max_g_power + 1):
        for q in range(p + 1):
            vector = reduced_feature_vector(p, q, degree)
            weight = sp.Rational(comb(p, q), p * (p - 1))
            matrix += weight * vector * vector.T
    return matrix


def exact_sos_hunt() -> dict[str, object]:
    """Produce one exact non-diagonal SOS block and its exact obstruction."""
    successful = reduced_sos_block(5, 4)
    lower, diagonal = successful.LDLdecomposition(hermitian=True)
    assert lower * diagonal * lower.T == successful
    pivots = [sp.factor(diagonal[index, index]) for index in range(diagonal.rows)]
    assert all(pivot > 0 for pivot in pivots)

    failed = reduced_sos_block(7, 6)
    failed_determinant = sp.factor(failed.det())
    assert failed_determinant == sp.Rational(-713423754787, 22394880000000)

    # If x=(1-s)s^2 and y=(1-s)^2s^2, then (1-s)s^3=x-y.
    small_matrix = sp.diag(sp.Rational(1, 2), sp.Rational(1, 2))
    difference = sp.Matrix([1, -1])
    small_matrix -= sp.Rational(1, 6) * difference * difference.T
    expected_small = sp.Matrix(
        [
            [sp.Rational(1, 3), sp.Rational(1, 6)],
            [sp.Rational(1, 6), sp.Rational(1, 3)],
        ]
    )
    assert small_matrix == expected_small

    return {"pivots": pivots, "failed_determinant": failed_determinant}


def exact_factorisation_checks() -> dict[str, object]:
    """Check every symbolic identity in the Lorentz/Gram proof."""
    aa, bb = sp.symbols("A B")
    g = lambda argument: argument + (1 - argument) * sp.log(1 - argument)
    phi = g(aa * (1 + bb)) + aa * g(-bb)
    second_derivative = aa / ((1 + bb) * (1 - aa * (1 + bb)))
    assert sp.simplify(phi.subs(bb, 0) - g(aa)) == 0
    assert sp.simplify(sp.diff(phi, bb).subs(bb, 0) + aa * sp.log(1 - aa)) == 0
    assert sp.simplify(sp.diff(phi, bb, 2) - second_derivative) == 0

    s, t, theta, lam = sp.symbols("s t theta lambda")
    denominator = sp.expand(
        (1 + theta * s * t)
        * (1 - (1 - s) * (1 - t) * (1 + theta * s * t))
    )
    polynomial = sp.Poly(denominator, s, t)
    coefficient_matrix = sp.zeros(4)
    for (row, column), coefficient in polynomial.terms():
        coefficient_matrix[row, column] = coefficient
    expected_matrix = sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, -theta - 1, 2 * theta, 0],
            [0, 2 * theta, -theta * (theta + 2), theta**2],
            [0, 0, theta**2, -theta**2],
        ]
    )
    assert (coefficient_matrix - expected_matrix).applyfunc(sp.expand) == sp.zeros(4)
    assert sp.factor(coefficient_matrix.det()) == -2 * theta**3

    characteristic_coefficients = [sp.factor(value) for value in coefficient_matrix.charpoly(lam).all_coeffs()]
    expected_characteristic = [
        sp.Integer(1),
        (theta + 1) * (2 * theta + 1),
        4 * theta**3 + 2 * theta - 1,
        -2 * theta * (theta**3 - theta**2 + theta + 1),
        -2 * theta**3,
    ]
    assert [sp.expand(value) for value in characteristic_coefficients] == [
        sp.expand(value) for value in expected_characteristic
    ]

    permutation = [1, 0, 2, 3]
    permuted = coefficient_matrix.extract(permutation, permutation)
    lower = sp.Matrix(
        [
            [1, 0, 0, 0],
            [-1 / (theta + 1), 1, 0, 0],
            [-2 * theta / (theta + 1), 2 * theta, 1, 0],
            [0, 0, -theta / (theta + 2), 1],
        ]
    )
    diagonal = sp.diag(
        -theta - 1,
        1 / (theta + 1),
        -theta * (theta + 2),
        -2 * theta**2 / (theta + 2),
    )
    assert (lower * diagonal * lower.T - permuted).applyfunc(sp.factor) == sp.zeros(4)

    permuted_monomials = sp.Matrix([s, 1, s**2, s**3])
    lorentz_coordinates = (lower.T * permuted_monomials).applyfunc(sp.factor)
    expected_coordinates = sp.Matrix(
        [
            -(2 * theta * s**2 - (theta + 1) * s + 1) / (theta + 1),
            1 + 2 * theta * s**2,
            s**2 * (theta + 2 - theta * s) / (theta + 2),
            s**3,
        ]
    )
    assert (lorentz_coordinates - expected_coordinates).applyfunc(sp.factor) == sp.zeros(4, 1)

    cubic = 2 - 2 * s + 2 * s**2 - s**3
    assert sp.discriminant(sp.diff(cubic, s), s) == -8
    assert cubic.subs(s, 0) == 2 and cubic.subs(s, 1) == 1
    diagonal_gap = sp.factor(1 - (1 - s) ** 2 * (1 + s**2))
    assert sp.expand(diagonal_gap - s * cubic) == 0

    return {
        "coefficient_matrix": coefficient_matrix,
        "characteristic_coefficients": expected_characteristic,
    }


def weighted_minimum_eigenvalue(kernel: np.ndarray, weights: np.ndarray) -> float:
    root_weight = np.sqrt(weights)
    matrix = root_weight[:, None] * kernel * root_weight[None, :]
    return float(np.linalg.eigvalsh((matrix + matrix.T) / 2.0)[0])


def numerical_factorisation_check(fingerprint: dict[str, object]) -> dict[str, float]:
    """Reconstruct -R from the exact integral identity in float64."""
    a = np.asarray(fingerprint["a"])
    b = np.asarray(fingerprint["b"])
    minus_r = np.asarray(fingerprint["minus_r"])

    base_g = g_np(a)
    base_log = a * b * (-np.log1p(-a))
    base = base_g + base_log
    theta_nodes, theta_weights = leggauss(THETA_QUADRATURE_ORDER)
    theta_nodes = (theta_nodes + 1.0) / 2.0
    theta_weights = theta_weights / 2.0
    integral = np.zeros_like(a)
    for theta, weight in zip(theta_nodes, theta_weights):
        denominator = (1.0 + theta * b) * (1.0 - a * (1.0 + theta * b))
        integral += weight * (1.0 - theta) * a / denominator
    reconstructed = base + b**2 * integral
    reconstruction_error = float(np.max(np.abs(reconstructed - minus_r)))
    assert reconstruction_error < 1.0e-12

    small_nodes, small_weights = leggauss(80)
    small_s = (small_nodes + 1.0) / 2.0
    small_weights = small_weights / 2.0
    small_a, small_b, *_ = kernel_arrays(small_s)
    base_g_minimum = weighted_minimum_eigenvalue(g_np(small_a), small_weights)
    base_log_minimum = weighted_minimum_eigenvalue(
        small_a * small_b * (-np.log1p(-small_a)), small_weights
    )
    worst_integrand_eigenvalue = 0.0
    for theta in (0.0, 0.25, 0.5, 0.75, 1.0):
        denominator = (1.0 + theta * small_b) * (
            1.0 - small_a * (1.0 + theta * small_b)
        )
        integrand = small_a * small_b**2 / denominator
        minimum = weighted_minimum_eigenvalue(integrand, small_weights)
        worst_integrand_eigenvalue = min(worst_integrand_eigenvalue, minimum)
    assert min(base_g_minimum, base_log_minimum, worst_integrand_eigenvalue) > -1.0e-12

    return {
        "reconstruction_error": reconstruction_error,
        "base_g_minimum": base_g_minimum,
        "base_log_minimum": base_log_minimum,
        "integrand_minimum": worst_integrand_eigenvalue,
    }


# Exact moment machinery for the Arb compression.  Q is the series index.
Q = sp.symbols("Q")


@lru_cache(maxsize=None)
def shifted_legendre_j(n: int) -> sp.Expr:
    """Rational J_n(Q)=integral_0^1 s^Q P_n(2s-1) ds."""
    numerator = sp.prod(Q - offset for offset in range(n))
    denominator = sp.prod(Q + offset for offset in range(1, n + 2))
    return sp.cancel(numerator / denominator)


@lru_cache(maxsize=None)
def shifted_legendre_u_moment(n: int) -> sp.Expr:
    """Rational integral of (1-s)s^Q P_n(2s-1)."""
    value = shifted_legendre_j(n)
    return sp.cancel(value - value.subs(Q, Q + 1))


def sympy_rational_to_fraction(value: sp.Expr) -> Fraction:
    value = sp.Rational(value)
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def alternating_tail_parts(pole_shift: int, power: int) -> tuple[Fraction, ...]:
    """Return rational, log(2), pi^2 coefficients of an alternating tail."""
    first_index = 2 + pole_shift
    assert first_index >= 1 and power in (1, 2)
    finite = sum(
        (
            Fraction(-1 if index % 2 else 1, index**power)
            for index in range(1, first_index)
        ),
        Fraction(0),
    )
    sign = -1 if pole_shift % 2 else 1
    rational = -sign * finite
    log_coefficient = Fraction(-sign) if power == 1 else Fraction(0)
    pi_squared_coefficient = Fraction(-sign, 12) if power == 2 else Fraction(0)
    return rational, log_coefficient, pi_squared_coefficient


def sum_alternating_rational_function(expression: sp.Expr) -> tuple[Fraction, ...]:
    """Exactly sum sum_(Q>=2) (-1)^Q expression(Q)."""
    partial_fractions = sp.apart(expression, Q)
    rational = Fraction(0)
    log_coefficient = Fraction(0)
    pi_squared_coefficient = Fraction(0)
    for term in sp.Add.make_args(partial_fractions):
        coefficient, remainder = term.as_coeff_Mul()
        exact_coefficient = sympy_rational_to_fraction(coefficient)
        if isinstance(remainder, sp.Pow):
            base = remainder.base
            power = -int(remainder.exp)
        else:
            base = remainder
            power = 1
        polynomial = sp.Poly(base, Q)
        assert polynomial.degree() == 1 and polynomial.LC() == 1
        pole_shift = int(polynomial.TC())
        tail = alternating_tail_parts(pole_shift, power)
        rational += exact_coefficient * tail[0]
        log_coefficient += exact_coefficient * tail[1]
        pi_squared_coefficient += exact_coefficient * tail[2]
    return rational, log_coefficient, pi_squared_coefficient


@lru_cache(maxsize=None)
def exact_f_moment_parts(n: int, m: int) -> tuple[Fraction, ...]:
    """Exact unnormalised Legendre matrix entry for A F(B)."""
    first = shifted_legendre_u_moment(n)
    second = shifted_legendre_u_moment(m)
    rational_function = -first * second / (Q * (Q - 1))
    return sum_alternating_rational_function(rational_function)


@lru_cache(maxsize=None)
def shifted_legendre_coefficients(n: int) -> tuple[int, ...]:
    """Integer coefficients of P_n(2s-1) in ascending powers of s."""
    return tuple(
        (-1) ** (n - k) * comb(n, k) * comb(n + k, k)
        for k in range(n + 1)
    )


def exact_feature_moments(p: int, q: int, degree: int) -> list[Fraction]:
    """Integrals of u^p s^q against shifted Legendre polynomials."""
    beta_moments = [Fraction(1, (p + q + 1) * comb(p + q, p))]
    for k in range(degree):
        beta_moments.append(
            beta_moments[-1] * Fraction(q + k + 1, p + q + k + 2)
        )
    return [
        sum(
            (
                Fraction(coefficient) * beta_moments[k]
                for k, coefficient in enumerate(shifted_legendre_coefficients(n))
            ),
            Fraction(0),
        )
        for n in range(degree + 1)
    ]


def to_fmpq(value: Fraction) -> fmpq:
    return fmpq(value.numerator, value.denominator)


def exact_truncated_g_matrix(degree: int, cutoff: int) -> fmpq_mat:
    """Exact matrix of sum_(p=2..cutoff) z^p/[p(p-1)]."""
    columns: list[list[Fraction]] = []
    weights: list[Fraction] = []
    for p in range(2, cutoff + 1):
        for q in range(p + 1):
            columns.append(exact_feature_moments(p, q, degree))
            weights.append(Fraction(comb(p, q), p * (p - 1)))

    rows = []
    weighted_rows = []
    for n in range(degree + 1):
        row = [to_fmpq(column[n]) for column in columns]
        weighted_row = [
            to_fmpq(column[n] * weight)
            for column, weight in zip(columns, weights)
        ]
        rows.append(row)
        weighted_rows.append(weighted_row)
    feature_matrix = fmpq_mat(rows)
    weighted_feature_matrix = fmpq_mat(weighted_rows)
    return feature_matrix * weighted_feature_matrix.transpose()


def arb_from_fraction(value: Fraction) -> arb:
    return arb(to_fmpq(value))


def arb_compression_bounds() -> dict[int, tuple[object, object]]:
    """Certify upper bounds using an exact finite Loewner majorant of R."""
    degree = 16
    ctx.prec = ARB_PRECISION_BITS
    log_two = arb(2).log()
    pi_squared = arb.pi() ** 2

    f_entries: list[arb] = []
    for n in range(degree + 1):
        for m in range(degree + 1):
            rational, log_coefficient, pi_coefficient = exact_f_moment_parts(n, m)
            entry = (
                arb_from_fraction(rational)
                + arb_from_fraction(log_coefficient) * log_two
                + arb_from_fraction(pi_coefficient) * pi_squared
            )
            f_entries.append(entry)
    f_matrix = arb_mat(degree + 1, degree + 1, f_entries)
    g_matrix = arb_mat(exact_truncated_g_matrix(degree, ARB_G_CUTOFF))
    comparison = f_matrix - g_matrix

    # sqrt(2n+1) P_n(2s-1) is orthonormal on [0,1].
    for n in range(degree + 1):
        for m in range(degree + 1):
            comparison[n, m] *= (arb(2 * n + 1) * arb(2 * m + 1)).sqrt()

    results: dict[int, tuple[object, object]] = {}
    for requested_degree in (8, 12, 16):
        dimension = requested_degree + 1
        entries = [
            comparison[row, column]
            for row in range(dimension)
            for column in range(dimension)
        ]
        compressed = arb_mat(dimension, dimension, entries)
        eigenvalues = compressed.eig(algorithm="rump")
        assert all(value.imag.contains(0) for value in eigenvalues)
        top = max(eigenvalues, key=lambda value: float(value.real.upper()))
        upper = top.real.upper()
        assert upper < 0
        results[requested_degree] = (top.real, upper)
    return results


def print_fingerprint(result: dict[str, object]) -> None:
    descending = np.asarray(result["descending"])
    eigenvalues = np.asarray(result["eigenvalues"])
    print("1. MERCER / EIGENFUNCTION FINGERPRINT")
    print(
        "NUMERICAL [float64 %d-node Gauss--Legendre Nystrom matrix]:"
        % QUADRATURE_ORDER
    )
    print("   largest eigenvalues of -R:")
    print("   " + "  ".join(f"{value:.12e}" for value in descending[:12]))
    print(
        "   smallest computed eigenvalue = %+.3e (roundoff-scale diagnostic)"
        % eigenvalues[0]
    )
    print(
        "   modes 3..20: log(lambda_k)=a%+.6f sqrt(k), R^2=%.6f"
        % (result["root_slope"], result["root_r2"])
    )
    print(
        "   comparison fits: geometric R^2=%.6f; power-law R^2=%.6f"
        % (result["geometric_r2"], result["power_r2"])
    )
    print(
        "   observed decay is root-exponential over those resolved modes; "
        "this is not an asymptotic theorem."
    )
    print("NUMERICAL [weighted L2 correlations; residuals are relative because ||phi||=1]:")
    print("   s^k(1-s) and (1-s)s^k are the same requested candidate family.")
    for fit in result["fits"]:
        monomial_name, monomial_corr = fit["monomial"]
        endpoint_name, endpoint_corr = fit["endpoint"]
        legendre_name, legendre_corr = fit["legendre"]
        jacobi_name, jacobi_corr = fit["jacobi"]
        log_name, log_corr = fit["log"]
        print(f"   mode {fit['mode']}: best single atoms")
        print(
            f"      polynomial {monomial_name} ({monomial_corr:.6f}); "
            f"endpoint {endpoint_name} ({endpoint_corr:.6f})"
        )
        print(
            f"      Legendre {legendre_name} ({legendre_corr:.6f}); "
            f"Jacobi {jacobi_name} ({jacobi_corr:.6f}); "
            f"log-type {log_name} ({log_corr:.6f})"
        )
        print(
            "      degree-12 span residuals: polynomial %.3e, "
            "(1-s) polynomial %.3e, log span %.3e"
            % (fit["poly_residual"], fit["endpoint_residual"], fit["log_residual"])
        )
    print(
        "CONJECTURED [only from the preceding fits]: the leading Mercer "
        "functions are not members of the tested classical closed-form families."
    )
    print()


def print_series_and_sos(sos: dict[str, object]) -> None:
    print("2. EXACT DOUBLE SERIES AND FINITE SOS HUNT")
    print("PROVED [Taylor series plus the binomial theorem, checked by SymPy]:")
    print("   R = sum_(p,q) c_(p,q) A^p B^q, with")
    print("     c_(1,q) = (-1)^(q+1)/[q(q-1)]                  (q >= 2),")
    print("     c_(p,q) = -binom(p,q)/[p(p-1)]       (p >= 2, 0 <= q <= p),")
    print("     c_(p,q) = 0 otherwise.")
    print("   Thus the exact positive set is {(1,q): q odd and q >= 3}.")
    print("   Every p>=2 coefficient is negative; at p=1 the even-q set is negative.")
    print("   Exact sample rows (Fraction arithmetic):")
    print(
        "     p=1, q=2..11: "
        + ", ".join(str(series_coefficient(1, q)) for q in range(2, 12))
    )
    for p in range(2, 6):
        print(
            f"     p={p}, q=0..{p}: "
            + ", ".join(str(series_coefficient(p, q)) for q in range(p + 1))
        )

    print("PROVED [explicit non-rank-one 2-by-2 regrouping]:")
    print("   Put x=(1-s)s^2 and y=(1-s)^2s^2, so (1-s)s^3=x-y.")
    print("   The q=2, q=3 and (p,q)=(2,2) terms of -R are")
    print("     (1/2)x@x + (1/2)y@y - (1/6)(x-y)@(x-y)")
    print("       = (1/4)(x+y)@(x+y) + (1/12)(x-y)@(x-y),")
    print("   an exact SOS that genuinely mixes two different (p,q) features.")
    print("PROVED [exact rational LDL^T, all pivots positive]:")
    print("   The joint prefix q<=5, p<=4 of -R is an SOS after factoring uv.")
    print("   Rational LDL pivots:")
    print("     " + ", ".join(str(pivot) for pivot in sos["pivots"]))
    print("PROVED [exact determinant obstruction to the natural next prefix]:")
    print("   For q<=7, p<=6 the reduced coefficient determinant is")
    print(f"     {sos['failed_determinant']} < 0,")
    print("   so that particular finite prefix is not an SOS.  This refutes only")
    print("   the natural prefix regrouping, not arbitrary infinite regroupings.")
    print()


def print_factorisation_proof(
    symbolic: dict[str, object], numerical: dict[str, float]
) -> None:
    print("3. SUBSTITUTION / EXPLICIT GRAM FACTORISATION")
    print(
        "NUMERICAL [%d-node theta Gauss--Legendre reconstruction]:"
        % THETA_QUADRATURE_ORDER
    )
    print(
        "   max absolute error in the claimed integral identity = %.3e"
        % numerical["reconstruction_error"]
    )
    print(
        "   weighted min eig G(A)=%+.3e; weighted min eig -AB ln(1-A)=%+.3e"
        % (numerical["base_g_minimum"], numerical["base_log_minimum"])
    )
    print(
        "   worst sampled integrand eigenvalue (theta=0,.25,.5,.75,1) = %+.3e"
        % numerical["integrand_minimum"]
    )
    print("   All diagnostics pass the requested 1e-12 threshold.")

    print("PROVED [two exact derivatives and Taylor's theorem with integral remainder]:")
    print("   Since G(-B)=-F(B), put phi(B)=G(A(1+B))+A G(-B)=-R.  SymPy gives")
    print("     phi(0)=G(A),        phi'(0)=-A ln(1-A),")
    print("     phi''(B)=A/[(1+B)(1-A(1+B))].")
    print("   Consequently")
    print("     -R = G(A)-AB ln(1-A)")
    print("          + integral_0^1 (1-theta) A B^2 / D_theta(s,t) dtheta,")
    print("     D_theta=(1+theta st)[1-uv(1+theta st)].")
    print("PROVED [nonnegative power-series Gram maps]:")
    print("   G(A)=sum_(n>=2) A^n/[n(n-1)] is PSD, and")
    print("   -AB ln(1-A)=sum_(n>=1) A^(n+1)B/n is PSD.")

    print("PROVED [exact SymPy coefficient extraction]:")
    print("   In the basis (1,s,s^2,s^3), D_theta has coefficient matrix")
    print(f"     {symbolic['coefficient_matrix']}")
    print("   det M(theta)=-2 theta^3, and its characteristic coefficients are")
    print("     1, (theta+1)(2theta+1), 4theta^3+2theta-1,")
    print("     -2theta(theta^3-theta^2+theta+1), -2theta^3.")
    print("PROVED [Descartes' rule, using symmetry of M(theta)]:")
    print("   For 0<theta<=1 the signs are (+,+,+/-,-,-), hence exactly one")
    print("   positive eigenvalue; det!=0 then gives inertia (1 positive, 3 negative).")
    print("   At theta=0 the inertia is (1 positive, 1 negative, 2 zero).")

    print("PROVED [exact rational LDL^T, giving an explicit Lorentz factor]:")
    print("   In the permuted basis (s,1,s^2,s^3), the pivots are")
    print("     -(theta+1), 1/(theta+1), -theta(theta+2), -2theta^2/(theta+2).")
    print("   Define")
    print("     z0=-(2theta s^2-(theta+1)s+1)/(theta+1),")
    print("     z1=1+2theta s^2,")
    print("     z2=s^2(theta+2-theta s)/(theta+2),   z3=s^3,")
    print("     a_theta=z1/sqrt(theta+1),")
    print("     b_theta=(sqrt(theta+1)z0,")
    print("              sqrt(theta(theta+2))z2,")
    print("              theta sqrt(2/(theta+2))z3).")
    print("   Then D_theta(s,t)=a_theta(s)a_theta(t)-b_theta(s).b_theta(t).")

    print("PROVED [exact cubic monotonicity check]:")
    print("   For 0<s<=1 and 0<=theta<=1,")
    print("     D_theta(s,s) >= s p(s),  p(s)=2-2s+2s^2-s^3.")
    print("   p'(s)=-3s^2+4s-2 has discriminant -8 and is negative; p decreases")
    print("   from 2 to 1.  Thus D_theta(s,s)>0, a_theta is nonzero with constant")
    print("   sign, and ||c_theta(s)||<1 for c_theta=b_theta/a_theta.")
    print("PROVED [geometric tensor-power Gram series]:")
    print("     1/D_theta(s,t) = 1/[a_theta(s)a_theta(t)]")
    print("       * sum_(n>=0) <c_theta(s)^tensor n, c_theta(t)^tensor n>.")
    print("   Cauchy--Schwarz gives absolute convergence.  Multiplication by")
    print("   A B^2=(1-s)s^2(1-t)t^2 and integration against 1-theta preserve PSD.")
    print("   At (0,0), B^2/D_theta extends continuously by zero; PSD follows by limits.")
    print("   Also AB[-ln(1-A)] extends by zero there: with m=max(s,t),")
    print("   st[-ln(s+t-st)] <= m^2[-ln m] -> 0.  Hence the full identity extends.")
    print("PROVED [the preceding explicit Gram decomposition]: R is NSD on [0,1].")
    print()


def print_arb_bounds(bounds: dict[int, tuple[object, object]]) -> None:
    print("4. RIGOROUS FINITE ARB COMPRESSIONS")
    print(
        "PROVED [exact Legendre moments, exact Fraction/fmpq arithmetic, "
        f"{ARB_PRECISION_BITS}-bit Arb]:"
    )
    print("   The comparison kernel is B_K=A F(B)-sum_(p=2..K) z^p/[p(p-1)]")
    print(f"   with K={ARB_G_CUTOFF}.  Since the omitted G tail is PSD, R <= B_K.")
    print("   Basis: sqrt(2n+1) P_n(2s-1), orthonormal on [0,1].")
    for degree in (8, 12, 16):
        enclosure, upper = bounds[degree]
        print(f"   degree <= {degree:2d}: lambda_max(B_K) in {enclosure.str(24)}")
        print(f"                 hence lambda_max(R) <= {upper.str(24)} < 0")
    print("PROVED [logic of finite sections]: each displayed inequality applies only")
    print("   to its stated polynomial subspace.  Finite compressions alone cannot")
    print("   prove an infinite-dimensional kernel claim; Section 3 supplies that proof.")
    print()


def main() -> None:
    print("LIU RESIDUAL KERNEL: EXACT NSD PROOF AND CERTIFIED COMPRESSIONS")
    print("Status labels are PROVED, NUMERICAL, or CONJECTURED.")
    print()

    fingerprint = mercer_fingerprint()
    print_fingerprint(fingerprint)

    exact_series_checks()
    sos = exact_sos_hunt()
    print_series_and_sos(sos)

    symbolic_factorisation = exact_factorisation_checks()
    numerical_factorisation = numerical_factorisation_check(fingerprint)
    print_factorisation_proof(symbolic_factorisation, numerical_factorisation)

    bounds = arb_compression_bounds()
    print_arb_bounds(bounds)

    print("FINAL STATUS")
    print("PROVED [exact integral Gram factorisation]: R <= 0 without projection.")
    print("PROVED [the exact reduction in liu_kernel.py]: Liu's Hypothesis 1 is settled.")
    print("PROVED [scope statement]: Liu's Hypothesis 2 (the nine-parameter Section V-B")
    print("   optimisation) remains open, so 0.382709087918741 is still conditional.")
    print("PROVED [Arb]: the degree-8, degree-12, and degree-16 compressions have")
    print("   strictly negative certified upper bounds as belt-and-braces checks.")


if __name__ == "__main__":
    main()
