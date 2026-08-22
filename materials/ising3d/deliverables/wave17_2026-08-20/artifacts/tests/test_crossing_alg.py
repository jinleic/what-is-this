"""Clean-room verifier for exact open-2x3 forced-crossing isolation.

The producer ``e96_crossing_alg.py`` is deliberately never imported.  This
program independently rebuilds the integer polynomial transfer representative,
its character-sector polynomials, the stored rational Krawczyk certificates,
and the exact inertia/branch-transport checks.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from typing import Iterable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "spectral" / "crossing_alg.json"
PRIOR = ROOT / "results" / "spectral" / "interval2.json"

T, X, Z = sp.symbols("t x z")
ROWS, COLS, SITES, DIMENSION, MAX_DEGREE = 2, 3, 6, 64, 26
T_RADIUS = Fraction(5, 10**11)
EIGEN_RADIUS = Fraction(2, 10**10)
INITIAL_RADIUS = Fraction(1, 10**10)
TUBE_RADIUS = Fraction(5, 10**9)
ORDER_MARGIN = Fraction(1, 10**7)


@dataclass(frozen=True)
class Basis:
    character: int
    vectors: tuple[dict[int, int], ...]
    gram: tuple[int, ...]


def q(value: str) -> Fraction:
    return Fraction(value)


def digest(poly: sp.Poly) -> str:
    content = json.dumps(
        [(list(exponents), str(coefficient)) for exponents, coefficient in sorted(poly.as_dict().items())],
        separators=(",", ":"),
    )
    return sha256(content.encode()).hexdigest()


def layer_coefficients() -> list[list[list[int]]]:
    """Independent raw-state construction of the degree-26 integer M(t)."""

    def at(row: int, col: int) -> int:
        return row * COLS + col

    bonds = [(at(row, col), at(row, col + 1)) for row in range(ROWS) for col in range(COLS - 1)]
    bonds += [(at(row, col), at(row + 1, col)) for row in range(ROWS - 1) for col in range(COLS)]
    energies = []
    for state in range(DIMENSION):
        spins = [1 - 2 * ((state >> (SITES - 1 - index)) & 1) for index in range(SITES)]
        energies.append(sum(spins[left] * spins[right] for left, right in bonds))
    powers = [(energy - (energies[0] % 2)) // 2 for energy in energies]
    assert min(powers) == -4 and max(powers) == 3
    shift, scalar_degree = 4, 7
    coefficients = [[[0] * (MAX_DEGREE + 1) for _ in range(DIMENSION)] for _ in range(DIMENSION)]
    for row in range(DIMENSION):
        for column in range(row, DIMENSION):
            polynomial = [0] * (MAX_DEGREE + 1)
            for middle, power in enumerate(powers):
                q_power = power + shift
                base = scalar_degree - q_power + (row ^ middle).bit_count() + (middle ^ column).bit_count()
                for binomial_index in range(q_power + 1):
                    polynomial[base + 2 * binomial_index] += 2 ** (scalar_degree - q_power) * math.comb(q_power, binomial_index)
            coefficients[row][column] = polynomial
            coefficients[column][row] = polynomial
    assert all(value >= 0 for row in coefficients for entry in row for value in entry)
    assert all(entry[index] == entry[MAX_DEGREE - index] for row in coefficients for entry in row for index in range(MAX_DEGREE + 1))
    return coefficients


def permutation(mapping: Sequence[int]) -> list[int]:
    answer = []
    for state in range(DIMENSION):
        new_state = 0
        for old, new in enumerate(mapping):
            new_state |= ((state >> (SITES - 1 - old)) & 1) << (SITES - 1 - new)
        answer.append(new_state)
    return answer


def clean_bases() -> tuple[Basis, ...]:
    def at(row: int, col: int) -> int:
        return row * COLS + col

    row_flip = permutation([at(ROWS - 1 - row, col) for row in range(ROWS) for col in range(COLS)])
    col_flip = permutation([at(row, COLS - 1 - col) for row in range(ROWS) for col in range(COLS)])
    spin_flip = [state ^ 63 for state in range(DIMENSION)]
    generators = (row_flip, col_flip, spin_flip)
    group = []
    for mask in range(8):
        element = list(range(DIMENSION))
        for bit, generator in enumerate(generators):
            if mask & (1 << bit):
                element = [generator[state] for state in element]
        group.append(element)
    representatives, visited = [], set()
    for state in range(DIMENSION):
        if state not in visited:
            orbit = {element[state] for element in group}
            representatives.append(min(orbit))
            visited |= orbit
    bases = []
    for character in range(8):
        vectors, gram = [], []
        for representative in representatives:
            vector: dict[int, int] = {}
            for group_index, element in enumerate(group):
                sign = -1 if (character & group_index).bit_count() % 2 else 1
                image = element[representative]
                vector[image] = vector.get(image, 0) + sign
            vector = {state: value for state, value in vector.items() if value}
            if vector:
                vectors.append(vector)
                gram.append(sum(value * value for value in vector.values()))
        bases.append(Basis(character, tuple(vectors), tuple(gram)))
    assert [len(basis.gram) for basis in bases] == [14, 6, 6, 6, 10, 10, 6, 6]
    return tuple(bases)


def reciprocal_table() -> tuple[sp.Expr, ...]:
    values: list[sp.Expr] = [sp.Integer(2), Z]
    for _ in range(2, 14):
        values.append(sp.expand(Z * values[-1] - values[-2]))
    return tuple(values)


def to_z(entry: Sequence[int], table: Sequence[sp.Expr]) -> sp.Expr:
    answer: sp.Expr = sp.Integer(entry[13])
    for offset in range(1, 14):
        answer += int(entry[13 + offset]) * table[offset]
    return sp.expand(answer)


def charpoly_z(coefficients: list[list[list[int]]], basis: Basis, table: Sequence[sp.Expr]) -> sp.Poly:
    block = []
    for left in basis.vectors:
        row = []
        for right in basis.vectors:
            value: sp.Expr = sp.Integer(0)
            for i, left_value in left.items():
                for j, right_value in right.items():
                    value += left_value * right_value * to_z(coefficients[i][j], table)
            row.append(value)
        block.append(row)
    determinant = sp.Poly((X * sp.diag(*basis.gram) - sp.Matrix(block)).det(method="domain-ge"), X, Z, domain=sp.ZZ)
    return determinant.primitive()[1]


def low_factor(charpoly: sp.Poly) -> tuple[sp.Poly, sp.Poly | None]:
    _content, factors = sp.factor_list(charpoly.as_expr(), X, Z)
    nonlinear = []
    linear = None
    for factor, multiplicity in factors:
        primitive = sp.Poly(factor, X, Z, domain=sp.ZZ).primitive()[1]
        assert multiplicity == 1
        if primitive.degree(X) == 1:
            linear = primitive
        else:
            nonlinear.append(primitive)
    assert len(nonlinear) == 1
    return nonlinear[0], linear


def restore(poly_z: sp.Poly, eigen_degree: int) -> sp.Poly:
    expression: sp.Expr = sp.Integer(0)
    for (x_power, z_power), coefficient in poly_z.as_dict().items():
        t_power = 13 * eigen_degree - 13 * x_power - z_power
        assert t_power >= 0
        expression += int(coefficient) * X**x_power * T**t_power * (1 + T * T) ** z_power
    return sp.Poly(sp.expand(expression), X, T, domain=sp.ZZ).primitive()[1]


Term = tuple[tuple[int, ...], int]


def embedded_terms(poly: sp.Poly, position: int) -> list[Term]:
    answer: list[Term] = []
    for (x_power, t_power), coefficient in poly.as_dict().items():
        exponents = [0] * 5
        exponents[0], exponents[position] = t_power, x_power
        answer.append((tuple(exponents), int(coefficient)))
    return answer


def partial(terms: Iterable[Term], coordinate: int) -> list[Term]:
    answer: list[Term] = []
    for exponents, coefficient in terms:
        if exponents[coordinate]:
            reduced = list(exponents)
            reduced[coordinate] -= 1
            answer.append((tuple(reduced), coefficient * exponents[coordinate]))
    return answer


def ev(terms: Iterable[Term], point: Sequence[Fraction], absolute: bool = False) -> Fraction:
    answer = Fraction(0)
    for exponents, coefficient in terms:
        term = abs(coefficient) if absolute else coefficient
        for value, power in zip(point, exponents, strict=True):
            term *= value**power
        answer += term
    return answer


def mm(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [[sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0)) for j in range(len(right[0]))] for i in range(len(left))]


def mv(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    return [sum((row[index] * vector[index] for index in range(len(vector))), Fraction(0)) for row in matrix]


def verify_krawczyk(record: dict[str, object], system: Sequence[list[Term]]) -> None:
    certificate = record["krawczyk"]
    assert isinstance(certificate, dict)
    center = [q(value) for value in certificate["center"]]
    radii = [q(value) for value in certificate["radii"]]
    preconditioner = [[q(value) for value in row] for row in certificate["preconditioner"]]
    upper = [value + radius for value, radius in zip(center, radii, strict=True)]
    jac = [[partial(system[row], column) for column in range(5)] for row in range(5)]
    hessian = [[[partial(jac[row][column], direction) for direction in range(5)] for column in range(5)] for row in range(5)]
    jacobian = [[ev(jac[row][column], center) for column in range(5)] for row in range(5)]
    f0 = [ev(polynomial, center) for polynomial in system]
    variation = [[sum((ev(hessian[row][column][direction], upper, True) * radii[direction] for direction in range(5)), Fraction(0)) for column in range(5)] for row in range(5)]
    correction = mm(preconditioner, jacobian)
    identity_error = [[Fraction(int(row == column)) - correction[row][column] for column in range(5)] for row in range(5)]
    propagated = mm([[abs(value) for value in row] for row in preconditioner], variation)
    q_matrix = [[abs(identity_error[row][column]) + propagated[row][column] for column in range(5)] for row in range(5)]
    image = [center[index] - value for index, value in enumerate(mv(preconditioner, f0))]
    containment = [
        (abs(image[row] - center[row]) + sum((q_matrix[row][column] * radii[column] for column in range(5)), Fraction(0))) / radii[row]
        for row in range(5)
    ]
    norm = max(sum((q_matrix[row][column] * radii[column] / radii[row] for column in range(5)), Fraction(0)) for row in range(5))
    assert all(value < q(certificate["coordinate_containment_upper_bound"]) for value in containment)
    assert norm < q(certificate["contraction_upper_bound"]) < 1


def two_terms(poly: sp.Poly) -> list[tuple[tuple[int, int], int]]:
    return [(tuple(map(int, powers)), int(coefficient)) for powers, coefficient in poly.as_dict().items()]


def partial_two(terms: Iterable[tuple[tuple[int, int], int]], coordinate: int) -> list[tuple[tuple[int, int], int]]:
    answer = []
    for exponents, coefficient in terms:
        if exponents[coordinate]:
            new = list(exponents)
            new[coordinate] -= 1
            answer.append((tuple(new), coefficient * exponents[coordinate]))
    return answer


def ev_two(terms: Iterable[tuple[tuple[int, int], int]], x_value: Fraction, t_value: Fraction, absolute: bool = False) -> Fraction:
    return sum(((abs(coefficient) if absolute else coefficient) * x_value**powers[0] * t_value**powers[1] for powers, coefficient in terms), Fraction(0))


def verify_implicit_derivative(poly: sp.Poly, center_x: Fraction, center_t: Fraction, expected: dict[str, object]) -> None:
    terms = two_terms(poly)
    px, pt = partial_two(terms, 0), partial_two(terms, 1)
    pxx, pxt = partial_two(px, 0), partial_two(px, 1)
    ptx, ptt = partial_two(pt, 0), partial_two(pt, 1)
    upper_x, upper_t = center_x + TUBE_RADIUS, center_t + T_RADIUS
    delta_x = ev_two(pxx, upper_x, upper_t, True) * TUBE_RADIUS + ev_two(pxt, upper_x, upper_t, True) * T_RADIUS
    delta_t = ev_two(ptx, upper_x, upper_t, True) * TUBE_RADIUS + ev_two(ptt, upper_x, upper_t, True) * T_RADIUS
    denominator = abs(ev_two(px, center_x, center_t)) - delta_x
    numerator = abs(ev_two(pt, center_x, center_t)) + delta_t
    assert denominator > 0
    actual = numerator / denominator
    derivative_upper = q(expected["derivative_upper_bound"])
    assert actual < derivative_upper
    motion = INITIAL_RADIUS + derivative_upper * T_RADIUS
    assert motion == q(expected["refined_motion_upper"])
    assert motion < EIGEN_RADIUS

def linear_factor_outside_window(
    linear_factor_z: sp.Poly | None, center_t: Fraction, lower: Fraction, upper: Fraction
) -> bool:
    """Check that an omitted linear factor has no root in a low-root window."""

    if linear_factor_z is None:
        return True
    linear_t = restore(linear_factor_z, 1)
    fixed_t = sp.Rational(center_t.numerator, center_t.denominator)
    polynomial_x = sp.Poly(linear_t.as_expr().subs(T, fixed_t), X, domain=sp.QQ)
    root = -polynomial_x.nth(0) / polynomial_x.nth(1)
    return not (
        sp.Rational(lower.numerator, lower.denominator)
        <= root
        <= sp.Rational(upper.numerator, upper.denominator)
    )


def inertia(lower: list[list[int]]) -> int | None:
    size = len(lower)
    work = [row[:] for row in lower]
    previous, negatives = 1, 0
    for pivot_index in range(size):
        pivot = work[pivot_index][pivot_index]
        if pivot == 0:
            return None
        if pivot * previous < 0:
            negatives += 1
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, row + 1):
                numerator = pivot * work[row][column] - work[row][pivot_index] * work[column][pivot_index]
                value, remainder = divmod(numerator, previous)
                assert remainder == 0
                work[row][column] = value
        previous = pivot
    return negatives


def matrix_at(coefficients: list[list[list[int]]], value: Fraction) -> tuple[list[list[int]], int]:
    numerator, denominator = value.numerator, value.denominator
    powers_numerator = [numerator**degree for degree in range(MAX_DEGREE + 1)]
    powers_denominator = [denominator ** (MAX_DEGREE - degree) for degree in range(MAX_DEGREE + 1)]
    return (
        [
            [sum(coefficient * powers_numerator[degree] * powers_denominator[degree] for degree, coefficient in enumerate(entry)) for entry in row]
            for row in coefficients
        ],
        denominator**MAX_DEGREE,
    )


def blocks_at(matrix: list[list[int]], bases: Sequence[Basis]) -> list[tuple[list[list[int]], tuple[int, ...]]]:
    answer = []
    for basis in bases:
        size = len(basis.gram)
        block = [[0] * size for _ in range(size)]
        for row in range(size):
            for column in range(row + 1):
                block[row][column] = block[column][row] = sum(
                    left_coefficient * matrix[left_state][right_state] * right_coefficient
                    for left_state, left_coefficient in basis.vectors[row].items()
                    for right_state, right_coefficient in basis.vectors[column].items()
                )
        answer.append((block, basis.gram))
    return answer


def counts(blocks: Sequence[tuple[list[list[int]], tuple[int, ...]]], scale: int, level: Fraction) -> list[int]:
    answer = []
    for block, gram in blocks:
        shifted = [[0] * len(gram) for _ in gram]
        for row in range(len(gram)):
            for column in range(row + 1):
                entry = level.denominator * block[row][column]
                if row == column:
                    entry -= level.numerator * scale * gram[row]
                shifted[row][column] = entry
        result = inertia(shifted)
        assert result is not None
        answer.append(result)
    return answer


def sqrt_upper(value: Fraction, denominator: int = 10**18) -> Fraction:
    integer = (value.numerator * denominator * denominator + value.denominator - 1) // value.denominator
    numerator = math.isqrt(integer)
    if numerator * numerator < integer:
        numerator += 1
    return Fraction(numerator, denominator)


def full_derivative(coefficients: list[list[list[int]]], upper: Fraction) -> Fraction:
    frobenius = Fraction(0)
    max_row = Fraction(0)
    for row in coefficients:
        row_sum = Fraction(0)
        for entry in row:
            derivative = sum((power * coefficient * upper ** (power - 1) for power, coefficient in enumerate(entry) if power), Fraction(0))
            frobenius += derivative * derivative
            row_sum += derivative
        max_row = max(max_row, row_sum)
    return min(sqrt_upper(frobenius), max_row)


def sector_coefficients(coefficients: list[list[list[int]]], bases: Sequence[Basis]) -> list[list[list[list[int]]]]:
    output = []
    for basis in bases:
        block = [[[0] * (MAX_DEGREE + 1) for _ in basis.gram] for _ in basis.gram]
        for row, left in enumerate(basis.vectors):
            for column, right in enumerate(basis.vectors):
                for i, left_coefficient in left.items():
                    for j, right_coefficient in right.items():
                        for power, coefficient in enumerate(coefficients[i][j]):
                            block[row][column][power] += left_coefficient * right_coefficient * coefficient
        output.append(block)
    return output


def sector_derivative(block: list[list[list[int]]], gram: Sequence[int], upper: Fraction) -> Fraction:
    frobenius = Fraction(0)
    for row in range(len(gram)):
        for column in range(len(gram)):
            derivative = sum((abs(power * coefficient) * upper ** (power - 1) for power, coefficient in enumerate(block[row][column]) if power), Fraction(0))
            frobenius += derivative * derivative / (gram[row] * gram[column])
    return sqrt_upper(frobenius)


def verify_inherited_brackets(data: dict[str, object]) -> None:
    prior = json.loads(PRIOR.read_text())
    old = prior["data"]["crossings"]
    new = data["crossings"]
    assert len(old) == len(new) == 4
    for old_record, new_record in zip(old, new, strict=True):
        assert old_record["fine_bracket_t"]["lower"] == new_record["fine_bracket"]["lower"]
        assert old_record["fine_bracket_t"]["upper"] == new_record["fine_bracket"]["upper"]
        left, right = old_record["c_counts_left_to_right"]
        assert right == left + 1
        inherited = new_record["inherited_c_jump_existence"]
        assert inherited["source"] == "proofs/spectral_interval2.md and results/spectral/interval2.json"
        assert inherited["c_counts_left_to_right"] == [left, right]


def main() -> None:
    artifact = json.loads(RESULT.read_text())
    data = artifact["data"]
    assert data["status"] == "PARTIAL_THEOREM"
    assert data["minimal_polynomial_status"] == "unresolved_after_resultant_wall"
    assert data["certified_bracket_crossing_count"] == 4
    assert data["full_scan_crossing_count"] is None
    assert data["full_scan_exhaustiveness"]["claim_tag"] == "[UNRESOLVED]"
    assert data["resultant_attempt"]["status"] == "wall"
    verify_inherited_brackets(data)

    coefficients = layer_coefficients()
    bases = clean_bases()
    table = reciprocal_table()
    degrees = {0: 14, 4: 9, 5: 9, 6: 6, 7: 6}
    polynomials: dict[int, sp.Poly] = {}
    linear: dict[int, sp.Poly | None] = {}
    for sector, degree in degrees.items():
        characteristic = charpoly_z(coefficients, bases[sector], table)
        low, linear_factor = low_factor(characteristic)
        restored = restore(low, degree)
        metadata = data["sector_polynomials"][str(sector)]
        assert characteristic.degree(X) == metadata["characteristic_x_degree"]
        assert characteristic.degree(Z) == metadata["characteristic_z_degree"]
        assert len(characteristic.terms()) == metadata["characteristic_terms"]
        assert low.degree(X) == metadata["low_factor_x_degree"]
        assert low.degree(Z) == metadata["low_factor_z_degree"]
        assert digest(restored) == metadata["restored_t_sha256"]
        polynomials[sector], linear[sector] = restored, linear_factor

    coefficient_blocks = sector_coefficients(coefficients, bases)
    for record in data["crossings"]:
        number = int(record["name"].rsplit("_", 1)[1])
        center = [q(value) for value in record["krawczyk"]["center"]]
        radii = [q(value) for value in record["krawczyk"]["radii"]]
        assert radii == [T_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS, EIGEN_RADIUS]
        lower, upper = q(record["fine_bracket"]["lower"]), q(record["fine_bracket"]["upper"])
        assert lower <= center[0] <= upper
        assert center[0] - T_RADIUS <= lower <= upper <= center[0] + T_RADIUS
        d_sector = int(record["ordered_sectors"][f"lambda{record['ordering']['local_eigenvalue_windows'][-1]['global_index']}"])
        system = [
            embedded_terms(polynomials[0], 1),
            embedded_terms(polynomials[5], 2),
            embedded_terms(polynomials[7], 3),
            embedded_terms(polynomials[d_sector], 4),
            [((0, 0, 1, 1, 0), 1), ((0, 1, 0, 0, 1), -1)],
        ]
        verify_krawczyk(record, system)

        matrix, scale = matrix_at(coefficients, center[0])
        sector_blocks_at_center = blocks_at(matrix, bases)
        ordering = record["ordering"]
        thresholds = [q(value) for value in ordering["thresholds"]]
        assert [sum(counts(sector_blocks_at_center, scale, value)) for value in thresholds] == ordering["global_counts"]
        for level, expected in zip(thresholds, ordering["sector_counts"], strict=True):
            assert counts(sector_blocks_at_center, scale, level) == expected
            assert counts(sector_blocks_at_center, scale, level - ORDER_MARGIN) == expected
            assert counts(sector_blocks_at_center, scale, level + ORDER_MARGIN) == expected

        windows = ordering["local_eigenvalue_windows"]
        for window in windows:
            level_low, level_high = q(window["lower"]), q(window["upper"])
            below, above = counts(sector_blocks_at_center, scale, level_low), counts(sector_blocks_at_center, scale, level_high)
            assert sum(below) == window["global_index"]
            assert sum(above) == window["global_index"] + 1
            assert above[window["sector"]] == below[window["sector"]] + 1
            for sector in range(8):
                if sector != window["sector"]:
                    assert above[sector] == below[sector]
            assert linear_factor_outside_window(linear[window["sector"]], center[0], level_low, level_high)

        assert full_derivative(coefficients, upper) * T_RADIUS < ORDER_MARGIN
        for key, certificate in ordering["implicit_branch_certificates"].items():
            window = next(item for item in windows if item["variable"] == key)
            sector = window["sector"]
            assert sector_derivative(coefficient_blocks[sector], bases[sector].gram, upper) * T_RADIUS < ORDER_MARGIN
            index = {"a": 1, "b": 2, "c": 3, "d": 4}[key]
            verify_implicit_derivative(polynomials[sector], center[index], center[0], certificate)
        motion = INITIAL_RADIUS + Fraction(1, 20) * T_RADIUS
        forced = ordering["forced_value_interval_on_fine_bracket"]
        forced_lower = (center[2] - motion) * (center[3] - motion) / (center[1] + motion)
        forced_upper = (center[2] + motion) * (center[3] + motion) / (center[1] - motion)
        assert all(center[index] > motion for index in (1, 2, 3))
        assert q(forced["selected_branch_motion_upper"]) == motion
        assert q(forced["lower"]) == forced_lower
        assert q(forced["upper"]) == forced_upper
        assert q(forced["target_band_lower"]) == thresholds[-2]
        assert q(forced["target_band_upper"]) == thresholds[-1]
        assert thresholds[-2] < forced_lower <= forced_upper < thresholds[-1]
        assert forced["strictly_inside_target_band"]
        assert record["exact_algebraic_representation"]["minimal_polynomial"] is None
        assert record["exact_algebraic_representation"]["transverse_sector_product_intersection"]

    print("OK: four exact Krawczyk systems, sector factors, inertia orderings, and branch transport certificates reverified")


if __name__ == "__main__":
    main()
