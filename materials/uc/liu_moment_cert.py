"""Exact finite moment sections for Liu's Hypothesis 1.

On signed measures satisfying <mu,1> = <mu,s> = <mu,s^2> = 0, the
reduction studied here is

    ln(2) Q(mu)
      = sum_{j>=2} (-1)^(j+1) <mu, (1-s)s^j>^2 / (j(j-1))
        - sum_{k>=2} sum_{i=0}^k C(k,i) <mu, (1-s)^k s^i>^2 / (k(k-1)).

For a support-degree cutoff R, this script assembles the exact rational matrix
in the plain moments p_3,...,p_R, rigorously encloses its top eigenvalue with
Arb, and runs interval Cholesky checks.  It also checks the stronger conjecture
that the fully resummed reduced kernel is negative semidefinite without any
projection.

The support-degree sections are not proofs of the infinite inequality: they
omit rank-one terms of both signs.  The output records that obstruction rather
than hiding it behind a floating-point near-zero.

Run from the workspace root with

    math/.venv/bin/python math/uc/liu_moment_cert.py
"""

from __future__ import annotations

import os

# Keep every numerical backend used by this standalone script to one core.
for _thread_variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import random
from fractions import Fraction
from math import ceil, comb, gcd, lcm
from typing import NamedTuple

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat


DEGREES = (8, 12, 16, 24, 32)
ARB_PRECISION_BITS = 512
BOUND_DECIMAL_PLACES = 15
SELF_CHECK_DIGITS = 100

ctx.prec = ARB_PRECISION_BITS


class RankOneTerm(NamedTuple):
    """One retained weighted square of a moment g_{k,i}."""

    family: str
    k: int
    i: int
    weight: Fraction
    vector: tuple[int, ...]


def projected_coefficients(k: int, i: int, degree: int) -> tuple[int, ...]:
    """Coefficients of (1-s)^k s^i in the p_3,...,p_degree basis."""

    assert k + i <= degree
    coefficients = [0] * (degree - 2)
    for offset in range(k + 1):
        power = i + offset
        if power >= 3:
            coefficients[power - 3] += (-1) ** offset * comb(k, offset)
    return tuple(coefficients)


def retained_terms(degree: int) -> list[RankOneTerm]:
    """All terms in the reduction whose polynomial support degree is <= R."""

    terms: list[RankOneTerm] = []

    # The k=1 logarithmic family f_j=(1-s)s^j has support degree j+1.
    for j in range(2, degree):
        terms.append(
            RankOneTerm(
                "f",
                1,
                j,
                Fraction((-1) ** (j + 1), j * (j - 1)),
                projected_coefficients(1, j, degree),
            )
        )

    # The power-series family has 0 <= i <= k and support degree i+k.
    for k in range(2, degree + 1):
        for i in range(k + 1):
            if k + i <= degree:
                terms.append(
                    RankOneTerm(
                        "g",
                        k,
                        i,
                        Fraction(-comb(k, i), k * (k - 1)),
                        projected_coefficients(k, i, degree),
                    )
                )
    return terms


def assemble_exact(degree: int) -> tuple[list[list[Fraction]], list[RankOneTerm]]:
    """Assemble G_R exactly as a sum of rational rank-one matrices."""

    dimension = degree - 2
    matrix = [[Fraction(0) for _ in range(dimension)] for _ in range(dimension)]
    terms = retained_terms(degree)
    for term in terms:
        nonzero = [(index, value) for index, value in enumerate(term.vector) if value]
        for row, row_value in nonzero:
            for column, column_value in nonzero:
                matrix[row][column] += term.weight * row_value * column_value

    assert all(
        matrix[row][column] == matrix[column][row]
        for row in range(dimension)
        for column in range(dimension)
    )
    assert all(isinstance(value, Fraction) for row in matrix for value in row)
    return matrix, terms


def matrix_quadratic(matrix: list[list[Fraction]], vector: list[Fraction]) -> Fraction:
    """Evaluate x^T A x without leaving Fraction arithmetic."""

    return sum(
        matrix[row][column] * vector[row] * vector[column]
        for row in range(len(vector))
        for column in range(len(vector))
    )


def random_projected_measure(
    node_count: int = 16, seed: int = 7
) -> tuple[list[Fraction], list[Fraction]]:
    """A seeded random rational signed measure killing degrees zero through two."""

    rng = random.Random(seed)
    raw_weights = [0] * node_count
    # Each translated (1,-3,3,-1) stencil annihilates every quadratic on a
    # uniform grid.  Random integer combinations therefore satisfy all three
    # projection constraints exactly, without a floating-point projection.
    for start in range(node_count - 3):
        coefficient = rng.randint(-9, 9)
        for offset, stencil_value in enumerate((1, -3, 3, -1)):
            raw_weights[start + offset] += coefficient * stencil_value

    scale = sum(abs(value) for value in raw_weights)
    assert scale > 0
    nodes = [Fraction(index + 1, node_count + 1) for index in range(node_count)]
    weights = [Fraction(value, scale) for value in raw_weights]
    assert all(
        sum(weight * node**power for node, weight in zip(nodes, weights)) == 0
        for power in range(3)
    )
    return nodes, weights


def plain_moments(
    nodes: list[Fraction], weights: list[Fraction], degree: int
) -> list[Fraction]:
    return [
        sum(weight * node**power for node, weight in zip(nodes, weights))
        for power in range(3, degree + 1)
    ]


def retained_quadratic_direct(
    terms: list[RankOneTerm], nodes: list[Fraction], weights: list[Fraction]
) -> Fraction:
    """Evaluate retained moment squares directly, before basis expansion."""

    total = Fraction(0)
    for term in terms:
        moment = sum(
            weight * (1 - node) ** term.k * node**term.i
            for node, weight in zip(nodes, weights)
        )
        total += term.weight * moment * moment
    return total


def mp_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def direct_kernel_quadratic(
    nodes: list[Fraction], weights: list[Fraction]
) -> mp.mpf:
    """Evaluate ln(2) times Liu's entropy quadratic at high precision."""

    total = mp.mpf(0)
    mp_nodes = [mp_fraction(node) for node in nodes]
    mp_weights = [mp_fraction(weight) for weight in weights]
    for s, weight_s in zip(mp_nodes, mp_weights):
        for t, weight_t in zip(mp_nodes, mp_weights):
            z = (1 - s) * (1 - t) * (1 + s * t)
            entropy_natural = -z * mp.log(z) - (1 - z) * mp.log1p(-z)
            total += weight_s * weight_t * entropy_natural
    return total


def fmpq_matrix(matrix: list[list[Fraction]]) -> fmpq_mat:
    return fmpq_mat(
        [[fmpq(value.numerator, value.denominator) for value in row] for row in matrix]
    )


def primitive_integer_vector(values: list[fmpq]) -> tuple[int, ...]:
    """Clear denominators and normalize an exact rational vector."""

    common_denominator = 1
    for value in values:
        common_denominator = lcm(common_denominator, int(value.q))
    integers = [int(value.p) * (common_denominator // int(value.q)) for value in values]
    common_factor = 0
    for value in integers:
        common_factor = gcd(common_factor, abs(value))
    assert common_factor > 0
    integers = [value // common_factor for value in integers]
    first_nonzero = next(value for value in integers if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def exact_nullspace(
    matrix: list[list[Fraction]],
) -> tuple[int, list[tuple[int, ...]]]:
    """Return exact rank and a primitive integer basis of the nullspace."""

    exact = fmpq_matrix(matrix)
    dimension = exact.ncols()
    rank = exact.rank()
    if rank == dimension:
        return rank, []

    reduced, reduced_rank = exact.rref()
    assert reduced_rank == rank
    pivot_columns: list[int] = []
    for row in range(rank):
        pivot = next(column for column in range(dimension) if reduced[row, column] != 0)
        assert reduced[row, pivot] == 1
        pivot_columns.append(pivot)
    free_columns = [column for column in range(dimension) if column not in pivot_columns]

    basis: list[tuple[int, ...]] = []
    for free in free_columns:
        vector = [fmpq(0) for _ in range(dimension)]
        vector[free] = fmpq(1)
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -reduced[row, free]
        integer_vector = primitive_integer_vector(vector)
        assert all(
            sum(
                matrix[row][column] * integer_vector[column]
                for column in range(dimension)
            )
            == 0
            for row in range(dimension)
        )
        basis.append(integer_vector)
    return rank, basis


def arb_fraction(value: Fraction) -> arb:
    """An outward-rounded Arb enclosure of an exact rational."""

    return arb(value.numerator) / value.denominator


def arb_exact_negative_shift(
    matrix: list[list[Fraction]], shift: Fraction = Fraction(0)
) -> arb_mat:
    """Enclose -G + shift*I entry by entry from exact Fractions."""

    dimension = len(matrix)
    entries: list[list[arb]] = []
    for row in range(dimension):
        output_row: list[arb] = []
        for column in range(dimension):
            value = -matrix[row][column]
            if row == column:
                value += shift
            output_row.append(arb_fraction(value))
        entries.append(output_row)
    return arb_mat(entries)


def arb_positive_matrix(matrix: list[list[Fraction]]) -> arb_mat:
    return arb_mat([[arb_fraction(value) for value in row] for row in matrix])


def negate_arb_matrix(matrix: arb_mat) -> arb_mat:
    return arb_mat(
        [[-matrix[row, column] for column in range(matrix.ncols())]
         for row in range(matrix.nrows())]
    )


def interval_cholesky(matrix: arb_mat) -> tuple[bool, int | None, arb | None]:
    """Prove positive definiteness by interval Cholesky pivots."""

    dimension = matrix.nrows()
    assert dimension == matrix.ncols()
    lower = [[arb(0) for _ in range(dimension)] for _ in range(dimension)]
    for row in range(dimension):
        pivot = matrix[row, row]
        for prior in range(row):
            pivot -= lower[row][prior] * lower[row][prior]
        if not (pivot > 0):
            return False, row, pivot
        lower[row][row] = pivot.sqrt()
        for following in range(row + 1, dimension):
            value = matrix[following, row]
            for prior in range(row):
                value -= lower[following][prior] * lower[row][prior]
            lower[following][row] = value / lower[row][row]
    return True, None, None


def rigorous_top_eigenvalue(matrix: arb_mat) -> arb:
    """Isolate and identify the largest eigenvalue using Arb's Rump method."""

    eigenvalues = matrix.eig(algorithm="rump")
    assert len(eigenvalues) == matrix.nrows()
    assert all(value.imag.contains(0) for value in eigenvalues)
    ordered = sorted(eigenvalues, key=lambda value: float(value.real.mid()))
    top = ordered[-1].real
    assert all(top.lower() > value.real.upper() for value in ordered[:-1])
    return top


def decimal_upper_candidate(value: arb, places: int) -> Fraction:
    """Find a displayed decimal rigorously above an Arb interval."""

    scale = 10**places
    numerator = ceil(float(value.mid()) * scale) + 2
    while not (arb_fraction(Fraction(numerator, scale)) > value):
        numerator += 1
    while arb_fraction(Fraction(numerator - 1, scale)) > value:
        numerator -= 1
    return Fraction(numerator, scale)


def decimal_string(value: Fraction, places: int) -> str:
    scale = 10**places
    scaled = value * scale
    assert scaled.denominator == 1 and scaled.numerator >= 0
    numerator = scaled.numerator
    return f"{numerator // scale}.{numerator % scale:0{places}d}"


def reduced_kernel_float(node_count: int) -> tuple[float, float]:
    """Float64 spectrum of the unprojected, fully resummed reduced kernel."""

    s = np.arange(1, node_count + 1, dtype=float) / (node_count + 1.0)
    u = 1.0 - s
    y = np.outer(s, s)
    z = np.outer(u, u) * (1.0 + y)
    first = np.outer(u, u) * (y - (1.0 + y) * np.log1p(y))
    second = z + (1.0 - z) * np.log1p(-z)
    matrix = first - second
    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues = np.linalg.eigvalsh(matrix)
    return float(eigenvalues[-1]), float(eigenvalues[0])


def reduced_kernel_arb(node_count: int) -> arb_mat:
    """Arb enclosure of the unprojected reduced kernel on rational nodes."""

    nodes = [arb_fraction(Fraction(index + 1, node_count + 1))
             for index in range(node_count)]
    entries = [[arb(0) for _ in range(node_count)] for _ in range(node_count)]
    for row, s in enumerate(nodes):
        u = 1 - s
        for column in range(row + 1):
            t = nodes[column]
            v = 1 - t
            y = s * t
            z = u * v * (1 + y)
            first = u * v * (y - (1 + y) * (1 + y).log())
            second = z + (1 - z) * (1 - z).log()
            value = first - second
            entries[row][column] = value
            entries[column][row] = value
    return arb_mat(entries)


def cauchy_kernel_counterexample() -> tuple[Fraction, tuple[tuple[Fraction, ...], ...]]:
    """Exact two-point witness against PSD of 1/(1+theta*s*t)."""

    theta = Fraction(1, 2)
    points = (Fraction(1, 4), Fraction(3, 4))
    matrix = tuple(
        tuple(Fraction(1, 1) / (1 + theta * s * t) for t in points)
        for s in points
    )
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    assert determinant < 0
    return determinant, matrix


def claim(label: str, text: str) -> None:
    """Print every substantive output line with its epistemic label."""

    assert label in {"PROVED", "NUMERICAL", "CONJECTURED"}
    for line in text.splitlines():
        print(f"{label}: {line}")


def main() -> None:
    mp.mp.dps = SELF_CHECK_DIGITS
    nodes, weights = random_projected_measure()
    direct = direct_kernel_quadratic(nodes, weights)

    claim(
        "PROVED",
        "The support-degree convention is max supp(f_j)=j+1 and "
        "max supp(g_{k,i})=k+i; coordinates are p_3,...,p_R.",
    )
    claim(
        "PROVED",
        "Every G_R entry below is assembled with Python Fraction arithmetic from "
        "the stated rational rank-one weights.",
    )
    print()

    results = []
    for degree in DEGREES:
        matrix, terms = assemble_exact(degree)
        moments = plain_moments(nodes, weights, degree)
        matrix_value = matrix_quadratic(matrix, moments)
        direct_retained = retained_quadratic_direct(terms, nodes, weights)
        assert matrix_value == direct_retained

        dimension = degree - 2
        rank, null_vectors = exact_nullspace(matrix)
        float_matrix = np.array([[float(value) for value in row] for row in matrix])
        float_top = float(np.linalg.eigvalsh(float_matrix)[-1])

        arb_matrix = arb_positive_matrix(matrix)
        top = rigorous_top_eigenvalue(arb_matrix)
        assert top > 0

        unshifted_passed, failed_index, failed_pivot = interval_cholesky(
            arb_exact_negative_shift(matrix)
        )
        assert not unshifted_passed
        assert failed_index is not None and failed_pivot is not None and failed_pivot < 0

        epsilon = decimal_upper_candidate(top, BOUND_DECIMAL_PLACES)
        while True:
            shifted_passed, _, _ = interval_cholesky(
                arb_exact_negative_shift(matrix, shift=epsilon)
            )
            if shifted_passed:
                break
            epsilon += Fraction(1, 10**BOUND_DECIMAL_PLACES)
        assert arb_fraction(epsilon) > top

        positive_count = sum(term.weight > 0 for term in terms)
        negative_count = sum(term.weight < 0 for term in terms)
        claim(
            "PROVED",
            f"R={degree}: exact symmetric G_R is {dimension}x{dimension}, with "
            f"{len(terms)} retained formal terms ({positive_count} positive and "
            f"{negative_count} negative weights); exact rank-one self-check passed.",
        )
        if null_vectors:
            claim(
                "PROVED",
                f"R={degree}: exact FLINT rank is {rank}; primitive exact null "
                f"vector basis in (p_3,...,p_R) coordinates is {null_vectors}.",
            )
        else:
            claim(
                "PROVED",
                f"R={degree}: exact FLINT rank is {rank}={dimension}; there is no "
                "exact null vector to quotient out.",
            )
        claim(
            "NUMERICAL",
            f"R={degree}: float64 largest eigenvalue = {float_top:+.12e} "
            "(orientation only; the monomial basis becomes severely ill-conditioned).",
        )
        claim(
            "PROVED",
            f"R={degree}: Arb {ARB_PRECISION_BITS}-bit isolated top-eigenvalue "
            f"enclosure is {top.str(25)}; it is strictly positive, so G_R <= 0 is false.",
        )
        claim(
            "PROVED",
            f"R={degree}: interval Cholesky of -G_R FAILED at coordinate "
            f"p_{failed_index + 3}, with strictly negative pivot "
            f"{failed_pivot.str(12)}.",
        )
        claim(
            "PROVED",
            f"R={degree}: interval Cholesky of -G_R + eps I PASSED for exact "
            f"eps={decimal_string(epsilon, BOUND_DECIMAL_PLACES)}; hence "
            f"lambda_max(G_R) < {decimal_string(epsilon, BOUND_DECIMAL_PLACES)}.",
        )
        print()
        results.append((degree, matrix_value, top, epsilon))

    claim(
        "PROVED",
        "The seeded rational signed self-check measure has p_0=p_1=p_2=0 exactly "
        "(third-difference construction).",
    )
    claim(
        "NUMERICAL",
        "Its direct high-precision ln(2)*entropy-kernel quadratic is "
        f"{mp.nstr(direct, 25)} (mpmath, {SELF_CHECK_DIGITS} decimal digits).",
    )
    for degree, matrix_value, _, _ in results:
        matrix_mp = mp_fraction(matrix_value)
        claim(
            "PROVED",
            f"R={degree}: exact p^T G_R p equals the independently evaluated "
            "retained rank-one moment sum in Fraction arithmetic.",
        )
        claim(
            "NUMERICAL",
            f"R={degree}: p^T G_R p={mp.nstr(matrix_mp, 20)}, "
            f"|direct-kernel minus truncation|={mp.nstr(abs(direct - matrix_mp), 8)}.",
        )
    print()

    claim(
        "PROVED",
        "For these even cutoffs, retained positive terms are exactly "
        "+<f_j>^2/[j(j-1)] for odd 3 <= j <= R-1.",
    )
    for degree in DEGREES:
        claim(
            "PROVED",
            f"R={degree}: omitted positive terms are exactly "
            f"{{+<f_j>^2/[j(j-1)] : j odd, j >= {degree + 1}}} "
            f"(j={degree + 1},{degree + 3},{degree + 5},...).",
        )
    claim(
        "PROVED",
        "For an even R, omitted negative terms are exactly the f_j terms with even "
        "j>=R and the g_{k,i} terms with k>=2, 0<=i<=k, and k+i>R.",
    )
    claim(
        "PROVED",
        "Deleting a negative weighted square raises the quadratic form, but deleting "
        "a positive weighted square lowers it.  Because support truncation deletes "
        "both signs, G_R is not a one-sided bound for the infinite form.",
    )
    claim(
        "PROVED",
        "The k>=2 power-series negative family is confined to the band 0<=i<=k, "
        "while every positive term is g_{1,j}=f_j with odd j>=3, outside that "
        "band.  Thus that negative family admits no index-local term-by-term pairing "
        "with the positive terms.",
    )
    claim(
        "PROVED",
        "Although g_{2,j}=f_j-f_{j+1} is a polynomial identity, Liu's negative "
        "k=2 sum contains only i=0,1,2.  It contains no g_{2,j} weight for j>=3, "
        "so the proposed direct telescoping Cauchy-Schwarz budget is absent.",
    )
    claim(
        "PROVED",
        "The M_{k,i} are not free coordinates: all are linear images of the same "
        "plain moments p_r.  Any valid domination must therefore use these global "
        "linear relations rather than compare moment indices term by term.",
    )
    claim(
        "PROVED",
        "The certified positive eigenvalues show that none of the five requested "
        "support-degree sections proves G_R<=0; omitted higher-support negative "
        "terms are structurally relevant to any joint tail argument.",
    )
    claim(
        "CONJECTURED",
        "The infinite moment form is <=0 and the even-cutoff positive residue tends "
        "to zero.  The finite calculations here do not prove either statement.",
    )
    print()

    claim(
        "PROVED",
        "Let y=st, u=1-s, v=1-t, z=uv(1+y), and "
        "R(s,t)=uv[y-(1+y)ln(1+y)]-[z+(1-z)ln(1-z)].  On projected measures, "
        "the original entropy form equals the quadratic form of this reduced kernel.",
    )
    claim(
        "PROVED",
        "G(z)=z+(1-z)ln(1-z)=sum_{k>=2} z^k/[k(k-1)] is a PSD kernel: "
        "z=uv+a(s)a(t) is PSD, each entrywise power is PSD by Schur products, "
        "and all series coefficients are nonnegative.  Hence -G(z) is NSD.",
    )
    claim(
        "PROVED",
        "For F(y)=y-(1+y)ln(1+y), F''(y)=-1/(1+y) and F(0)=F'(0)=0, so "
        "F(y)=-y^2*integral_0^1 (1-theta)/(1+theta*y) dtheta exactly.",
    )
    determinant, cauchy_matrix = cauchy_kernel_counterexample()
    claim(
        "PROVED",
        "The hoped-for Schur route is closed: at theta=1/2 and points 1/4,3/4, "
        f"the kernel 1/(1+theta*s*t) has Gram matrix {cauchy_matrix} with exact "
        f"determinant {determinant}<0, so it is not PSD.",
    )
    print()

    claim(
        "CONJECTURED",
        "Stronger target: the fully resummed reduced kernel R is NSD on all signed "
        "measures, without Liu's codimension-three projection.",
    )
    for node_count in (100, 200):
        float_top, float_bottom = reduced_kernel_float(node_count)
        claim(
            "NUMERICAL",
            f"Unprojected R on {node_count} rational interior nodes: float64 "
            f"lambda_max={float_top:+.12e}, lambda_min={float_bottom:+.12e}; "
            "the top is at roundoff scale.",
        )
    for node_count in (16, 32):
        reduced = reduced_kernel_arb(node_count)
        top = rigorous_top_eigenvalue(reduced)
        passed, _, _ = interval_cholesky(negate_arb_matrix(reduced))
        assert top < 0 and passed
        claim(
            "PROVED",
            f"Unprojected R on {node_count} exact rational nodes: Arb top-eigenvalue "
            f"enclosure {top.str(25)} is strictly negative, and interval Cholesky "
            "of -R PASSED.",
        )
    claim(
        "PROVED",
        "Each Arb grid result covers only vectors supported on that finite node set; "
        "it cannot imply the continuum conjecture, whose spectrum may accumulate at zero.",
    )
    claim(
        "CONJECTURED",
        "The finite-grid evidence suggests the unprojected reduced-kernel statement "
        "is the clean proof target.  No continuum NSD proof is claimed here.",
    )


if __name__ == "__main__":
    main()
