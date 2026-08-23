"""Standalone clean-room verifier for the finite CTM signature certificate.

This verifier imports none of experiments/e206_ctm_definition.py,
e207_ctm_spectrum.py, or e208_ctm_certificate.py.  It rebuilds the boundary
objects by a different enumeration order and obtains pair-product
characteristic polynomials from an explicit exterior-square matrix rather than
from the producer's power-trace identities.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "ctm_signature.json"
ALLOWED_TAGS = {
    "[THEOREM]",
    "[LEMMA]",
    "[COMPUTATION]",
    "[CONJECTURE]",
    "[EXTERNAL]",
    "[UNRESOLVED]",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def qtrim(polynomial: Sequence[int]) -> list[int]:
    out = [int(value) for value in polynomial]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out or [0]


def qadd(left: Sequence[int], right: Sequence[int]) -> list[int]:
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += int(value)
    for index, value in enumerate(right):
        out[index] += int(value)
    return qtrim(out)


def qneg(polynomial: Sequence[int]) -> list[int]:
    return [-int(value) for value in polynomial]


def qsub(left: Sequence[int], right: Sequence[int]) -> list[int]:
    return qadd(left, qneg(right))


def qmul(left: Sequence[int], right: Sequence[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            out[left_degree + right_degree] += int(left_value) * int(right_value)
    return qtrim(out)


def qpayload_sha256(payload: object) -> str:
    encoded = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def grid_vertices(dimension: int, radius: int) -> list[tuple[int, ...]]:
    return [tuple(vertex) for vertex in product(range(radius + 1), repeat=dimension)]


def grid_edges(
    dimension: int, radius: int
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    out: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for vertex in grid_vertices(dimension, radius):
        for axis in range(dimension):
            if vertex[axis] == radius:
                continue
            neighbour = list(vertex)
            neighbour[axis] += 1
            out.append((vertex, tuple(neighbour)))
    return out


def face_sites(dimension: int, radius: int, axis: int) -> list[tuple[int, ...]]:
    remaining = [coordinate for coordinate in range(dimension) if coordinate != axis]
    out: list[tuple[int, ...]] = []
    for local in product(range(radius + 1), repeat=dimension - 1):
        vertex = [0] * dimension
        for coordinate, value in zip(remaining, local):
            vertex[coordinate] = value
        out.append(tuple(vertex))
    return out


def boundary_polynomial(
    dimension: int, radius: int, face_states: Sequence[int]
) -> list[int]:
    assignments: dict[tuple[int, ...], int] = {}
    edge_list = grid_edges(dimension, radius)
    for axis, state in enumerate(face_states):
        for bit, vertex in enumerate(face_sites(dimension, radius, axis)):
            spin = (state >> bit) & 1
            if vertex in assignments and assignments[vertex] != spin:
                return [0] * (len(edge_list) + 1)
            assignments[vertex] = spin
    interior = [vertex for vertex in grid_vertices(dimension, radius) if vertex not in assignments]
    coefficients = [0] * (len(edge_list) + 1)
    for state in range(1 << len(interior)):
        spins = dict(assignments)
        for bit, vertex in enumerate(interior):
            spins[vertex] = (state >> bit) & 1
        broken = sum(spins[left] != spins[right] for left, right in edge_list)
        coefficients[broken] += 1
    return coefficients


def two_dimensional_polynomial_matrix() -> list[list[list[int]]]:
    matrix: list[list[list[int]]] = []
    for row in range(4):
        row_face = 1 | (row << 1)  # fixed + origin at bit zero
        entries: list[list[int]] = []
        for column in range(4):
            column_face = 1 | (column << 1)
            entries.append(boundary_polynomial(2, 2, (row_face, column_face)))
        matrix.append(entries)
    return matrix


def evaluate_scaled(polynomial: Sequence[int], denominator_power: int) -> int:
    return sum(
        int(coefficient) * 2 ** (denominator_power - degree)
        for degree, coefficient in enumerate(polynomial)
    )


def evaluate_matrix(
    matrix: Sequence[Sequence[Sequence[int]]], denominator_power: int
) -> list[list[int]]:
    return [
        [evaluate_scaled(entry, denominator_power) for entry in row] for row in matrix
    ]


def permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def qdeterminant(matrix: Sequence[Sequence[Sequence[int]]]) -> list[int]:
    determinant = [0]
    for permutation in permutations(range(len(matrix))):
        term = [permutation_sign(permutation)]
        for row, column in enumerate(permutation):
            term = qmul(term, matrix[row][column])
        determinant = qadd(determinant, term)
    return determinant


def zmul(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> list[list[int]]:
    """Multiply z-polynomials stored ascending, with coefficients in Z[q]."""
    out: list[list[int]] = [[0] for _ in range(len(left) + len(right) - 1)]
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            out[left_degree + right_degree] = qadd(
                out[left_degree + right_degree], qmul(left_value, right_value)
            )
    while len(out) > 1 and out[-1] == [0]:
        out.pop()
    return out


def zadd(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> list[list[int]]:
    out: list[list[int]] = [[0] for _ in range(max(len(left), len(right)))]
    for degree, coefficient in enumerate(left):
        out[degree] = qadd(out[degree], coefficient)
    for degree, coefficient in enumerate(right):
        out[degree] = qadd(out[degree], coefficient)
    while len(out) > 1 and out[-1] == [0]:
        out.pop()
    return out


def polynomial_matrix_charpoly_leibniz(
    matrix: Sequence[Sequence[Sequence[int]]],
) -> list[list[int]]:
    """Return det(zI-M), descending in z, by a bivariate Leibniz expansion."""
    dimension = len(matrix)
    determinant: list[list[int]] = [[0]]
    for permutation in permutations(range(dimension)):
        term: list[list[int]] = [[permutation_sign(permutation)]]
        for row, column in enumerate(permutation):
            constant = qneg(matrix[row][column])
            factor = [constant, [1]] if row == column else [constant]
            term = zmul(term, factor)
        determinant = zadd(determinant, term)
    return list(reversed(determinant))


def exterior_square_polynomial_matrix(
    matrix: Sequence[Sequence[Sequence[int]]],
) -> list[list[list[int]]]:
    pairs = list(product(range(len(matrix)), repeat=2))
    pairs = [(left, right) for left, right in pairs if left < right]
    out: list[list[list[int]]] = []
    for left, right in pairs:
        row: list[list[int]] = []
        for first, second in pairs:
            row.append(
                qsub(
                    qmul(matrix[left][first], matrix[right][second]),
                    qmul(matrix[left][second], matrix[right][first]),
                )
            )
        out.append(row)
    return out


def exterior_square_integer_matrix(
    matrix: Sequence[Sequence[int]],
) -> list[list[int]]:
    pairs = [(left, right) for left in range(len(matrix)) for right in range(left + 1, len(matrix))]
    return [
        [
            matrix[left][first] * matrix[right][second]
            - matrix[left][second] * matrix[right][first]
            for first, second in pairs
        ]
        for left, right in pairs
    ]


def integer_charpoly_leibniz(matrix: Sequence[Sequence[int]]) -> list[int]:
    polynomial_matrix = [[[int(value)] for value in row] for row in matrix]
    result = polynomial_matrix_charpoly_leibniz(polynomial_matrix)
    require(all(len(coefficient) == 1 for coefficient in result), "integer charpoly leaked q")
    return [coefficient[0] for coefficient in result]


def zdivide_linear_q(
    coefficients_descending: Sequence[Sequence[int]], root: Sequence[int]
) -> tuple[list[list[int]], list[int]]:
    quotient = [[int(value) for value in coefficients_descending[0]]]
    for coefficient in coefficients_descending[1:]:
        quotient.append(qadd(coefficient, qmul(root, quotient[-1])))
    remainder = quotient.pop()
    return quotient, qtrim(remainder)


def fraction_trim(polynomial: Sequence[Fraction]) -> list[Fraction]:
    out = [Fraction(value) for value in polynomial]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def fraction_divmod(
    dividend: Sequence[Fraction], divisor: Sequence[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    numerator = fraction_trim(dividend)
    denominator = fraction_trim(divisor)
    require(denominator != [0], "division by zero polynomial")
    if len(numerator) < len(denominator):
        return [Fraction(0)], numerator
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and numerator != [0]:
        shift = len(numerator) - len(denominator)
        factor = numerator[-1] / denominator[-1]
        quotient[shift] = factor
        for index, value in enumerate(denominator):
            numerator[index + shift] -= factor * value
        numerator = fraction_trim(numerator)
    return fraction_trim(quotient), numerator


def gcd_descending(polynomial: Sequence[int]) -> list[int]:
    degree = len(polynomial) - 1
    derivative = [int(polynomial[index]) * (degree - index) for index in range(degree)]
    first = [Fraction(value) for value in reversed(polynomial)]
    second = [Fraction(value) for value in reversed(derivative)]
    while fraction_trim(second) != [0]:
        _, remainder = fraction_divmod(first, second)
        first, second = second, remainder
    first = fraction_trim(first)
    first = [value / first[-1] for value in first]
    require(all(value.denominator == 1 for value in first), "nonintegral monic gcd")
    return [int(value) for value in reversed(first)]


def determinant_integer_leibniz(matrix: Sequence[Sequence[int]]) -> int:
    total = 0
    for permutation in permutations(range(len(matrix))):
        term = permutation_sign(permutation)
        for row, column in enumerate(permutation):
            term *= int(matrix[row][column])
        total += term
    return total


def leading_minors(matrix: Sequence[Sequence[int]]) -> list[int]:
    return [
        determinant_integer_leibniz([list(row[:order]) for row in matrix[:order]])
        for order in range(1, len(matrix) + 1)
    ]


def rational_rank(matrix: Sequence[Sequence[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rank + 1, len(work)):
            if work[row][column] == 0:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column, len(work[row])):
                work[row][index] -= factor * work[rank][index]
        rank += 1
    return rank


def second_moment(matrix: Sequence[Sequence[int]]) -> Fraction:
    trace = sum(matrix[index][index] for index in range(len(matrix)))
    trace_square = sum(
        matrix[row][column] * matrix[column][row]
        for row in range(len(matrix))
        for column in range(len(matrix))
    )
    return Fraction(len(matrix) * trace_square, trace * trace)


def three_dimensional_tensor() -> list[list[list[int]]]:
    """Independent route: enumerate 2^7 boundary assignments, not 16^3 leg triples."""
    faces = [face_sites(3, 1, axis) for axis in range(3)]
    boundary = sorted(set().union(*map(set, faces)))
    edge_list = grid_edges(3, 1)
    bulk = (1, 1, 1)
    tensor = [[[0] * 16 for _ in range(16)] for _ in range(16)]
    for boundary_state in range(1 << len(boundary)):
        fixed = {
            vertex: (boundary_state >> bit) & 1 for bit, vertex in enumerate(boundary)
        }
        leg_states: list[int] = []
        for face in faces:
            state = 0
            for bit, vertex in enumerate(face):
                state |= fixed[vertex] << bit
            leg_states.append(state)
        value = 0
        for bulk_spin in (0, 1):
            spins = dict(fixed)
            spins[bulk] = bulk_spin
            broken = sum(spins[left] != spins[right] for left, right in edge_list)
            value += 2 ** (len(edge_list) - broken)
        tensor[leg_states[0]][leg_states[1]][leg_states[2]] = value
    return tensor


def contract(tensor: Sequence[Sequence[Sequence[int]]], covector: Sequence[int]) -> list[list[int]]:
    return [
        [sum(tensor[row][column][state] * covector[state] for state in range(16)) for column in range(16)]
        for row in range(16)
    ]


def plus_edge_block(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    indices = [state for state in range(16) if state & 3 == 3]
    return [[matrix[row][column] for column in indices] for row in indices]


def transport_face_state(
    old_axis: int, state: int, permutation: Sequence[int]
) -> tuple[int, int]:
    new_axis = list(permutation).index(old_axis)
    old_face = face_sites(3, 1, old_axis)
    new_face = face_sites(3, 1, new_axis)
    positions = {vertex: bit for bit, vertex in enumerate(new_face)}
    result = 0
    for bit, vertex in enumerate(old_face):
        image = tuple(vertex[permutation[index]] for index in range(3))
        result |= ((state >> bit) & 1) << positions[image]
    return new_axis, result


def tensor_covariant(tensor: Sequence[Sequence[Sequence[int]]]) -> bool:
    for permutation in permutations(range(3)):
        maps = {
            (axis, state): transport_face_state(axis, state, permutation)
            for axis in range(3)
            for state in range(16)
        }
        for first in range(16):
            for second in range(16):
                for third in range(16):
                    image = [0, 0, 0]
                    for axis, state in enumerate((first, second, third)):
                        new_axis, new_state = maps[(axis, state)]
                        image[new_axis] = new_state
                    if tensor[first][second][third] != tensor[image[0]][image[1]][image[2]]:
                        return False
    return True

def permute_four_bits(state: int, permutation: Sequence[int]) -> int:
    out = 0
    for old_bit, new_bit in enumerate(permutation):
        out |= ((state >> old_bit) & 1) << new_bit
    return out



def fraction_payload(value: Fraction) -> dict[str, object]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    require(set(artifact) == {"meta", "data", "checks"}, "artifact top-level shape")
    require(artifact["meta"]["provenance"] == "experiments/e208_ctm_certificate.py", "provenance")
    require(artifact["meta"]["modular_arithmetic_used"] is False, "unexpected modular shortcut")
    require(artifact["meta"]["benchmarks_used_to_select_or_fit"] is False, "benchmark fitting")
    require(artifact["checks"] and all(check["passed"] for check in artifact["checks"]), "stored failed check")
    require(all(check["tag"] in ALLOWED_TAGS for check in artifact["checks"]), "untagged check")
    require(len({check["name"] for check in artifact["checks"]}) == len(artifact["checks"]), "duplicate check")
    for relative, expected_digest in artifact["meta"]["source_sha256"].items():
        require(sha256_file(ROOT / relative) == expected_digest, f"source digest {relative}")

    data = artifact["data"]
    require(data["finite_kernel_definition"]["tag"] == "[LEMMA]", "definition tag")
    require(data["conclusion"]["tag"] == "[UNRESOLVED]", "scope tag")
    require(data["conclusion"]["standard_3d_ising_solved"] is False, "solution overclaim")
    require(data["conclusion"]["infinite_volume_ctm_signature_decided"] is False, "limit overclaim")
    cross_checks = data["repository_enumerator_cross_checks"]

    matrix_q = two_dimensional_polynomial_matrix()
    matrix_2d = evaluate_matrix(matrix_q, 12)
    wedge_q = exterior_square_polynomial_matrix(matrix_q)
    pair_q_descending = polynomial_matrix_charpoly_leibniz(wedge_q)
    determinant_q = qdeterminant(matrix_q)
    root_q = [0] * 19
    for power, coefficient in enumerate((1, -6, 15, -20, 15, -6, 1)):
        root_q[6 + 2 * power] = coefficient
    first_quotient, first_remainder = zdivide_linear_q(pair_q_descending, root_q)
    second_quotient, second_remainder = zdivide_linear_q(first_quotient, root_q)
    _, third_remainder = zdivide_linear_q(second_quotient, root_q)

    control = data["two_dimensional_integrable_control"]
    symbolic = control["symbolic"]
    require(len(grid_edges(2, 2)) == control["edge_count"] == 12, "2D edge count")
    require(matrix_2d == control["scaled_matrix"], "2D matrix rebuild")
    require(
        2 * sum(sum(row) for row in matrix_2d)
        == cross_checks["two_dimensional_sum_of_fixed_corner_blocks"]
        == cross_checks["two_dimensional_open_3x3_scaled_partition"],
        "2D repository-enumerator cross-check",
    )
    require(qmul(root_q, root_q) == determinant_q, "2D symbolic determinant square")
    require(root_q == symbolic["determinant_root_coefficients_ascending"], "2D root coefficients")
    require(determinant_q == symbolic["determinant_coefficients_ascending"], "2D determinant coefficients")
    require(first_remainder == second_remainder == [0], "2D repeated pair factor")
    require(third_remainder != [0], "2D third factor unexpectedly vanishes")
    require(first_remainder == symbolic["first_remainder_coefficients_ascending"], "first remainder field")
    require(second_remainder == symbolic["second_remainder_coefficients_ascending"], "second remainder field")
    require(third_remainder == symbolic["third_remainder_coefficients_ascending"], "third remainder field")
    require([len(coefficient) - 1 for coefficient in pair_q_descending] == symbolic["pair_charpoly_coefficient_degrees"], "symbolic degree vector")
    require(qpayload_sha256(pair_q_descending) == symbolic["pair_charpoly_sha256"], "symbolic pair digest")
    require(qpayload_sha256(second_quotient) == symbolic["quotient_after_two_factors_sha256"], "symbolic quotient digest")

    charpoly_2d = integer_charpoly_leibniz(matrix_2d)
    pairpoly_2d = integer_charpoly_leibniz(exterior_square_integer_matrix(matrix_2d))
    minors_2d = leading_minors(matrix_2d)
    require(charpoly_2d == control["characteristic_polynomial_descending"], "2D charpoly")
    require(pairpoly_2d == control["pair_product_polynomial_descending"], "2D pair polynomial")
    require(gcd_descending(charpoly_2d) == control["characteristic_gcd_derivative_descending"] == [1], "2D simple spectrum")
    require(gcd_descending(pairpoly_2d) == control["pair_product_gcd_derivative_descending"] == [1, -46656], "2D pair gcd")
    require(minors_2d == control["leading_principal_minors"] and all(value > 0 for value in minors_2d), "2D positivity")
    require(control["distinct_pair_products"] == 5 and control["pair_product_slots"] == 6, "2D nonvacuity")

    tensor = three_dimensional_tensor()
    nonzero = sum(value != 0 for plane in tensor for row in plane for value in row)
    analogue = data["three_dimensional_honest_analogue"]
    require(len(grid_edges(3, 1)) == analogue["edge_count"] == 12, "3D edge count")
    require(nonzero == analogue["nonzero_compatible_entries"] == 128, "3D compatible tensor count")
    require(
        sum(value for plane in tensor for row in plane for value in row)
        == cross_checks["three_dimensional_boundary_tensor_sum"]
        == cross_checks["three_dimensional_open_2x2x2_scaled_partition"],
        "3D repository-enumerator cross-check",
    )
    require(tensor_covariant(tensor), "3D tensor coordinate covariance")
    require(analogue["intrinsic_characteristic_polynomial"] is None, "invented tensor spectrum")

    free = [1] * 16
    aligned = [int(state in (0, 15)) for state in range(16)]
    biased = [free[state] + aligned[state] for state in range(16)]
    require(all(value > 0 for value in free + biased), "positive covectors")
    require(
        all(
            covector[state] == covector[state ^ 15]
            for covector in (free, biased)
            for state in range(16)
        ),
        "spin-flip-invariant covectors",
    )
    require(
        all(
            covector[state] == covector[permute_four_bits(state, permutation)]
            for covector in (free, biased)
            for permutation in permutations(range(4))
            for state in range(16)
        ),
        "face-symmetric covectors",
    )
    require(
        any(
            free[left] * biased[right] != free[right] * biased[left]
            for left in range(16)
            for right in range(left + 1, 16)
        ),
        "nonproportional invariant covectors",
    )
    free_matrix = contract(tensor, free)
    aligned_matrix = contract(tensor, aligned)
    biased_matrix = contract(tensor, biased)
    free_block = plus_edge_block(free_matrix)
    biased_block = plus_edge_block(biased_matrix)

    spectralization = data["declared_free_face_spectralization"]
    charpoly_3d = integer_charpoly_leibniz(free_block)
    pairpoly_3d = integer_charpoly_leibniz(exterior_square_integer_matrix(free_block))
    minors_3d = leading_minors(free_block)
    require(free_block == spectralization["scaled_matrix"], "3D free-face block")
    require(charpoly_3d == spectralization["characteristic_polynomial_descending"], "3D charpoly")
    require(pairpoly_3d == spectralization["pair_product_polynomial_descending"], "3D pair polynomial")
    require(gcd_descending(charpoly_3d) == spectralization["characteristic_gcd_derivative_descending"] == [1], "3D simple spectrum")
    require(gcd_descending(pairpoly_3d) == spectralization["pair_product_gcd_derivative_descending"] == [1], "3D squarefree pair spectrum")
    require(minors_3d == spectralization["leading_principal_minors"] and all(value > 0 for value in minors_3d), "3D positivity")
    require(spectralization["distinct_pair_products"] == spectralization["pair_product_slots"] == 6, "3D pair verdict")

    obstruction = data["definition_level_obstruction"]
    require(obstruction["covectors"]["free"] == free, "free covector")
    require(obstruction["covectors"]["biased"] == biased, "biased covector")
    require(obstruction["covectors"]["aligned_counterexample"] == aligned, "aligned covector")
    ranks = {"free": rational_rank(free_matrix), "biased": rational_rank(biased_matrix), "aligned_counterexample": rational_rank(aligned_matrix)}
    require(ranks == obstruction["full_contraction_ranks"] == {"free": 16, "biased": 16, "aligned_counterexample": 8}, "contraction ranks")
    full_moments = {"free": fraction_payload(second_moment(free_matrix)), "biased": fraction_payload(second_moment(biased_matrix))}
    sector_moments = {"free": fraction_payload(second_moment(free_block)), "biased": fraction_payload(second_moment(biased_block))}
    require(full_moments == obstruction["full_scale_invariant_second_moments"], "full spectral ambiguity")
    require(sector_moments == obstruction["plus_edge_sector_scale_invariant_second_moments"], "sector spectral ambiguity")
    require(full_moments["free"] != full_moments["biased"], "contractions accidentally scale-equivalent")
    require(integer_charpoly_leibniz(biased_block) == obstruction["biased_plus_edge_characteristic_polynomial_descending"], "biased charpoly")

    print(
        "OK: clean-room exterior-square rebuild confirms the symbolic 2D product "
        "factor, the squarefree 3D free-face block, and the contraction ambiguity."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
