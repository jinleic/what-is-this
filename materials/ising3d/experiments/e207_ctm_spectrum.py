"""Exact characteristic and pair-product algebra for the finite CTM front.

No eigenvalue approximation and no logarithm is used.  A four-slot spectrum is
a two-mode subset-product spectrum exactly when, after relabelling nonzero simple
slots, one complementary pair product is repeated.  The characteristic
polynomial of the exterior square records those six pair products without
choosing an eigenvalue order.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import permutations
from math import isqrt
from typing import Sequence

IntPolynomial = list[int]  # ascending powers of q
IntMatrix = list[list[int]]
PolynomialMatrix = list[list[IntPolynomial]]


def _trim_ascending(polynomial: Sequence[int]) -> IntPolynomial:
    out = [int(value) for value in polynomial]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [0]


def polynomial_add(left: Sequence[int], right: Sequence[int]) -> IntPolynomial:
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += int(value)
    for index, value in enumerate(right):
        out[index] += int(value)
    return _trim_ascending(out)


def polynomial_negate(polynomial: Sequence[int]) -> IntPolynomial:
    return [-int(value) for value in polynomial]


def polynomial_subtract(left: Sequence[int], right: Sequence[int]) -> IntPolynomial:
    return polynomial_add(left, polynomial_negate(right))


def polynomial_multiply(left: Sequence[int], right: Sequence[int]) -> IntPolynomial:
    out = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            out[left_degree + right_degree] += int(left_value) * int(right_value)
    return _trim_ascending(out)


def polynomial_scale(polynomial: Sequence[int], scalar: int) -> IntPolynomial:
    return _trim_ascending([int(scalar) * int(value) for value in polynomial])


def polynomial_divide_integer_exact(
    polynomial: Sequence[int], divisor: int
) -> IntPolynomial:
    if divisor == 0 or any(int(value) % divisor for value in polynomial):
        raise ArithmeticError("non-exact coefficient division")
    return _trim_ascending([int(value) // divisor for value in polynomial])


def matrix_multiply(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> IntMatrix:
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    if any(len(row) != inner for row in left) or any(len(row) != columns for row in right):
        raise ValueError("incompatible or ragged matrices")
    out = [[0] * columns for _ in range(rows)]
    for row in range(rows):
        for pivot in range(inner):
            value = int(left[row][pivot])
            if value == 0:
                continue
            for column in range(columns):
                out[row][column] += value * int(right[pivot][column])
    return out


def matrix_power_traces(matrix: Sequence[Sequence[int]], count: int) -> list[int]:
    dimension = len(matrix)
    if dimension == 0 or any(len(row) != dimension for row in matrix):
        raise ValueError("a nonempty square matrix is required")
    power = [[int(value) for value in row] for row in matrix]
    traces: list[int] = []
    for exponent in range(1, count + 1):
        traces.append(sum(power[index][index] for index in range(dimension)))
        if exponent < count:
            power = matrix_multiply(power, matrix)
    return traces


def charpoly_from_power_traces(traces: Sequence[int]) -> list[int]:
    """Monic characteristic coefficients [1,c1,...,cn], descending in z."""
    coefficients = [1]
    for order in range(1, len(traces) + 1):
        numerator = sum(
            coefficients[order - power] * int(traces[power - 1])
            for power in range(1, order + 1)
        )
        if numerator % order:
            raise ArithmeticError("Newton identity did not divide over Z")
        coefficients.append(-numerator // order)
    return coefficients


def characteristic_polynomial(matrix: Sequence[Sequence[int]]) -> list[int]:
    return charpoly_from_power_traces(matrix_power_traces(matrix, len(matrix)))


def pair_product_characteristic_polynomial(
    matrix: Sequence[Sequence[int]],
) -> list[int]:
    """Return charpoly(wedge^2 matrix), whose roots are slot products lambda_i lambda_j."""
    dimension = len(matrix)
    pair_count = dimension * (dimension - 1) // 2
    traces = matrix_power_traces(matrix, 2 * pair_count)
    pair_traces: list[int] = []
    for exponent in range(1, pair_count + 1):
        numerator = traces[exponent - 1] ** 2 - traces[2 * exponent - 1]
        if numerator % 2:
            raise ArithmeticError("exterior-square trace is not integral")
        pair_traces.append(numerator // 2)
    return charpoly_from_power_traces(pair_traces)


def _trim_fraction_descending(polynomial: Sequence[Fraction]) -> list[Fraction]:
    out = [Fraction(value) for value in polynomial]
    first = 0
    while first < len(out) - 1 and out[first] == 0:
        first += 1
    return out[first:] or [Fraction(0)]


def _fraction_remainder_descending(
    dividend: Sequence[Fraction], divisor: Sequence[Fraction]
) -> list[Fraction]:
    work = _trim_fraction_descending(dividend)
    divisor = _trim_fraction_descending(divisor)
    if divisor == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(work) < len(divisor):
        return work
    for offset in range(len(work) - len(divisor) + 1):
        factor = work[offset] / divisor[0]
        if factor:
            for index, value in enumerate(divisor):
                work[offset + index] -= factor * value
    remainder_length = len(divisor) - 1
    return _trim_fraction_descending(
        work[-remainder_length:] if remainder_length else [Fraction(0)]
    )


def polynomial_gcd_over_q_descending(
    left: Sequence[int], right: Sequence[int]
) -> list[Fraction]:
    first = _trim_fraction_descending([Fraction(value) for value in left])
    second = _trim_fraction_descending([Fraction(value) for value in right])
    while second != [0]:
        first, second = second, _fraction_remainder_descending(first, second)
    if first == [0]:
        return first
    leading = first[0]
    return [value / leading for value in first]


def derivative_descending(polynomial: Sequence[int]) -> list[int]:
    degree = len(polynomial) - 1
    return [int(value) * (degree - index) for index, value in enumerate(polynomial[:-1])]


def squarefree_gcd_descending(polynomial: Sequence[int]) -> list[Fraction]:
    return polynomial_gcd_over_q_descending(polynomial, derivative_descending(polynomial))


def integral_fraction_polynomial(polynomial: Sequence[Fraction]) -> list[int]:
    if any(value.denominator != 1 for value in polynomial):
        raise ArithmeticError("the monic gcd is not integral")
    return [int(value) for value in polynomial]


def determinant_bareiss(matrix: Sequence[Sequence[int]]) -> int:
    dimension = len(matrix)
    if dimension == 0 or any(len(row) != dimension for row in matrix):
        raise ValueError("a nonempty square matrix is required")
    if dimension == 1:
        return int(matrix[0][0])
    work = [[int(value) for value in row] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(dimension - 1):
        pivot_row = next(
            (row for row in range(pivot_index, dimension) if work[row][pivot_index]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, dimension):
            for column in range(pivot_index + 1, dimension):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // previous
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def leading_principal_minors(matrix: Sequence[Sequence[int]]) -> list[int]:
    return [
        determinant_bareiss([list(row[:order]) for row in matrix[:order]])
        for order in range(1, len(matrix) + 1)
    ]


def exact_rank(matrix: Sequence[Sequence[int]]) -> int:
    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(row_count):
            if row == rank or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(column_count)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def scale_invariant_second_moment(matrix: Sequence[Sequence[int]]) -> Fraction:
    traces = matrix_power_traces(matrix, 2)
    dimension = len(matrix)
    if traces[0] == 0:
        raise ZeroDivisionError("trace-normalized invariant is undefined")
    return Fraction(dimension * traces[1], traces[0] ** 2)


def polynomial_matrix_multiply(
    left: Sequence[Sequence[Sequence[int]]],
    right: Sequence[Sequence[Sequence[int]]],
) -> PolynomialMatrix:
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    out: PolynomialMatrix = [[[0] for _ in range(columns)] for _ in range(rows)]
    for row in range(rows):
        for pivot in range(inner):
            for column in range(columns):
                out[row][column] = polynomial_add(
                    out[row][column],
                    polynomial_multiply(left[row][pivot], right[pivot][column]),
                )
    return out


def polynomial_matrix_power_traces(
    matrix: Sequence[Sequence[Sequence[int]]], count: int
) -> list[IntPolynomial]:
    dimension = len(matrix)
    power: PolynomialMatrix = [
        [[int(value) for value in entry] for entry in row] for row in matrix
    ]
    traces: list[IntPolynomial] = []
    for exponent in range(1, count + 1):
        trace = [0]
        for index in range(dimension):
            trace = polynomial_add(trace, power[index][index])
        traces.append(trace)
        if exponent < count:
            power = polynomial_matrix_multiply(power, matrix)
    return traces


def polynomial_charpoly_from_traces(
    traces: Sequence[Sequence[int]],
) -> list[IntPolynomial]:
    coefficients: list[IntPolynomial] = [[1]]
    for order in range(1, len(traces) + 1):
        numerator = [0]
        for power in range(1, order + 1):
            numerator = polynomial_add(
                numerator,
                polynomial_multiply(coefficients[order - power], traces[power - 1]),
            )
        coefficients.append(
            polynomial_negate(polynomial_divide_integer_exact(numerator, order))
        )
    return coefficients


def polynomial_pair_product_charpoly(
    matrix: Sequence[Sequence[Sequence[int]]],
) -> list[IntPolynomial]:
    dimension = len(matrix)
    pair_count = dimension * (dimension - 1) // 2
    traces = polynomial_matrix_power_traces(matrix, 2 * pair_count)
    pair_traces: list[IntPolynomial] = []
    for exponent in range(1, pair_count + 1):
        numerator = polynomial_subtract(
            polynomial_multiply(traces[exponent - 1], traces[exponent - 1]),
            traces[2 * exponent - 1],
        )
        pair_traces.append(polynomial_divide_integer_exact(numerator, 2))
    return polynomial_charpoly_from_traces(pair_traces)


def _permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def polynomial_determinant_leibniz(
    matrix: Sequence[Sequence[Sequence[int]]],
) -> IntPolynomial:
    dimension = len(matrix)
    determinant = [0]
    for permutation in permutations(range(dimension)):
        term = [_permutation_sign(permutation)]
        for row, column in enumerate(permutation):
            term = polynomial_multiply(term, matrix[row][column])
        determinant = polynomial_add(determinant, term)
    return determinant


def polynomial_square_root_exact(polynomial: Sequence[int]) -> IntPolynomial:
    target = _trim_ascending(polynomial)
    valuation = next((index for index, value in enumerate(target) if value), None)
    if valuation is None or valuation % 2:
        raise ArithmeticError("polynomial has no integral square root")
    base = target[valuation:]
    if (len(base) - 1) % 2:
        raise ArithmeticError("polynomial degree is incompatible with a square")
    constant = isqrt(base[0])
    if constant * constant != base[0]:
        raise ArithmeticError("constant coefficient is not a square")
    root = [constant]
    root_degree = (len(base) - 1) // 2
    for degree in range(1, root_degree + 1):
        target_coefficient = base[degree] if degree < len(base) else 0
        known = sum(root[index] * root[degree - index] for index in range(1, degree))
        numerator = target_coefficient - known
        denominator = 2 * constant
        if denominator == 0 or numerator % denominator:
            raise ArithmeticError("square-root recurrence did not divide over Z")
        root.append(numerator // denominator)
    candidate = [0] * (valuation // 2) + root
    candidate = _trim_ascending(candidate)
    if polynomial_multiply(candidate, candidate) != target:
        raise ArithmeticError("polynomial square-root audit failed")
    return candidate


def divide_z_polynomial_by_linear_q_factor(
    coefficients: Sequence[Sequence[int]], root: Sequence[int]
) -> tuple[list[IntPolynomial], IntPolynomial]:
    """Synthetic division of a monic z-polynomial by z-root(q)."""
    if len(coefficients) < 2:
        raise ValueError("constant z-polynomial")
    quotient: list[IntPolynomial] = [[int(value) for value in coefficients[0]]]
    for coefficient in coefficients[1:]:
        quotient.append(
            polynomial_add(coefficient, polynomial_multiply(root, quotient[-1]))
        )
    remainder = quotient.pop()
    return quotient, _trim_ascending(remainder)


def polynomial_payload_sha256(payload: object) -> str:
    encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def two_dimensional_symbolic_certificate(
    matrix: Sequence[Sequence[Sequence[int]]],
) -> dict[str, object]:
    """Derive, rather than assume, the repeated complementary-product factor."""
    determinant = polynomial_determinant_leibniz(matrix)
    determinant_root = polynomial_square_root_exact(determinant)
    pair_charpoly = polynomial_pair_product_charpoly(matrix)
    first_quotient, first_remainder = divide_z_polynomial_by_linear_q_factor(
        pair_charpoly, determinant_root
    )
    second_quotient, second_remainder = divide_z_polynomial_by_linear_q_factor(
        first_quotient, determinant_root
    )
    _, third_remainder = divide_z_polynomial_by_linear_q_factor(
        second_quotient, determinant_root
    )
    return {
        "determinant_coefficients_ascending": determinant,
        "determinant_root_coefficients_ascending": determinant_root,
        "pair_charpoly_coefficient_degrees": [len(entry) - 1 for entry in pair_charpoly],
        "pair_charpoly_sha256": polynomial_payload_sha256(pair_charpoly),
        "first_remainder_coefficients_ascending": first_remainder,
        "second_remainder_coefficients_ascending": second_remainder,
        "third_remainder_coefficients_ascending": third_remainder,
        "quotient_after_two_factors_sha256": polynomial_payload_sha256(second_quotient),
    }
