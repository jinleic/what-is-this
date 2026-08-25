#!/usr/bin/env python3
"""Clean-room verifier for the open-2x3x3 boundary matchgate defect.

This verifier imports neither e247 nor e241.  It reconstructs the open-box
edges by pairwise Manhattan distance, obtains the target signatures by a direct
retained-spin Walsh sum, obtains the planar control by another spin sum, and
enumerates boundary orders by factoradic unranking rather than the producer's
permutation iterator.
"""

from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "boundary_matchgate_defect.json"
PRODUCER = ROOT / "experiments" / "e247_boundary_matchgate_defect.py"
PROOF = ROOT / "proofs" / "boundary_matchgate_defect.md"
E241_PRODUCER = ROOT / "experiments" / "e241_box_genus_bound.py"
SHAPE = (2, 3, 3)
AUDIT_POINTS = (Fraction(1, 3), Fraction(1, 2))
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
X = sp.symbols("v")
FAILURES: list[str] = []

Coordinate = tuple[int, int, int]
Edge = tuple[Coordinate, Coordinate]
Poly = tuple[int, ...]
Signature = tuple[Poly, ...]

EXPECTED_FACTOR_DATA: dict[str, dict[str, object]] = {
    "R_v4_Q36": {
        "scalar": 4,
        "v_power": 4,
        "one_minus_v_power": 6,
        "one_plus_v_power": 6,
        "one_plus_v_squared_power": 2,
        "positive_remainder_coefficients_in_v_squared": [
            1,
            9,
            86,
            508,
            2389,
            8704,
            25555,
            60683,
            115407,
            178677,
            224201,
            214498,
            137351,
            59430,
            17373,
            3287,
            388,
            28,
            1,
        ],
        "first_ordering_count": 64,
        "all_four_subset_occurrences": 320,
    },
    "R_v4_Q38": {
        "scalar": 2,
        "v_power": 4,
        "one_minus_v_power": 6,
        "one_plus_v_power": 6,
        "one_plus_v_squared_power": 1,
        "positive_remainder_coefficients_in_v_squared": [
            1,
            17,
            190,
            1484,
            7961,
            34263,
            115785,
            315551,
            692059,
            1211057,
            1691429,
            1802207,
            1394971,
            751317,
            283139,
            72885,
            12784,
            1426,
            81,
            1,
        ],
        "first_ordering_count": 16,
        "all_four_subset_occurrences": 80,
    },
    "R_v6_Q28": {
        "scalar": 4,
        "v_power": 6,
        "one_minus_v_power": 8,
        "one_plus_v_power": 8,
        "one_plus_v_squared_power": 2,
        "positive_remainder_coefficients_in_v_squared": [
            1,
            30,
            199,
            990,
            3395,
            8616,
            17510,
            27180,
            30691,
            24110,
            13079,
            4406,
            785,
            76,
            4,
        ],
        "first_ordering_count": 32,
        "all_four_subset_occurrences": 160,
    },
    "R_v10_Q20": {
        "scalar": 8,
        "v_power": 10,
        "one_minus_v_power": 11,
        "one_plus_v_power": 11,
        "one_plus_v_squared_power": 1,
        "positive_remainder_coefficients_in_v_squared": [
            25,
            203,
            874,
            2266,
            3718,
            4018,
            2734,
            1178,
            293,
            47,
            4,
        ],
        "first_ordering_count": 8,
        "all_four_subset_occurrences": 40,
    },
}


def check(name: str, passed: bool, detail: str = "") -> None:
    suffix = f": {detail}" if detail else ""
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}{suffix}")
    if not passed:
        FAILURES.append(name)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        class TimeValue(ctypes.Structure):
            _fields_ = [
                ("seconds", ctypes.c_int32),
                ("microseconds", ctypes.c_int32),
            ]

        class MachTaskBasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_time", TimeValue),
                ("system_time", TimeValue),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        system = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        system.mach_task_self.restype = ctypes.c_uint32
        info = MachTaskBasicInfo()
        count = ctypes.c_uint32(
            ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_int32)
        )
        result = system.task_info(
            system.mach_task_self(),
            20,
            ctypes.byref(info),
            ctypes.byref(count),
        )
        if result != 0:
            raise OSError(f"mach task_info failed with status {result}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def peak_rss_measurement() -> str:
    if platform.system() == "Darwin":
        return "mach_task_basic_info.resident_size_max (task_info flavor 20)"
    return "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"


def trim(poly: Iterable[int]) -> Poly:
    coefficients = list(poly)
    if not coefficients:
        return (0,)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def add(left: Poly, right: Poly) -> Poly:
    return trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    )


def negate(poly: Poly) -> Poly:
    return tuple(-coefficient for coefficient in poly)


def subtract(left: Poly, right: Poly) -> Poly:
    return add(left, negate(right))


def multiply(left: Poly, right: Poly) -> Poly:
    coefficients = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            coefficients[left_index + right_index] += left_value * right_value
    return trim(coefficients)


def scale(poly: Poly, scalar: int) -> Poly:
    return trim(scalar * coefficient for coefficient in poly)


def evaluate(poly: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def exact_average(poly: Poly, denominator: int) -> Poly:
    if any(coefficient % denominator for coefficient in poly):
        raise ArithmeticError("Walsh average was not integral")
    return trim(coefficient // denominator for coefficient in poly)


def fraction_record(value: int | Fraction) -> dict[str, int]:
    rational = Fraction(value)
    return {
        "numerator": rational.numerator,
        "denominator": rational.denominator,
    }


def parity(vertex: Coordinate) -> int:
    return (vertex[0] ^ vertex[1] ^ vertex[2]) & 1


def clean_graph(shape: Sequence[int]) -> tuple[tuple[Coordinate, ...], tuple[Edge, ...]]:
    vertices = tuple(
        (x, y, z)
        for z in range(shape[2])
        for y in range(shape[1])
        for x in range(shape[0])
    )
    edges = tuple(
        sorted(
            (left, right) if left < right else (right, left)
            for left, right in itertools.combinations(vertices, 2)
            if sum(abs(left[axis] - right[axis]) for axis in range(3)) == 1
        )
    )
    return vertices, edges


def neighbour_map(
    vertices: Sequence[Coordinate], edges: Sequence[Edge]
) -> dict[Coordinate, tuple[Coordinate, ...]]:
    result = {vertex: [] for vertex in vertices}
    for left, right in edges:
        result[left].append(right)
        result[right].append(left)
    return {
        vertex: tuple(sorted(neighbours))
        for vertex, neighbours in result.items()
    }


def clean_boundary(
    vertices: Sequence[Coordinate], eliminated_parity: int
) -> tuple[tuple[Coordinate, ...], tuple[Coordinate, ...], dict[int, int], int]:
    retained = tuple(
        sorted(vertex for vertex in vertices if parity(vertex) != eliminated_parity)
    )
    by_face = {
        face: tuple(
            sorted(
                (vertex for vertex in retained if vertex[0] == face),
                key=lambda vertex: 3 * vertex[1] + vertex[2],
            )
        )
        for face in (0, 1)
    }
    selected = max(by_face, key=lambda face: len(by_face[face]))
    if sum(len(sites) == len(by_face[selected]) for sites in by_face.values()) != 1:
        raise AssertionError("target face maximum was not unique")
    boundary = by_face[selected]
    internal = tuple(vertex for vertex in retained if vertex not in boundary)
    return boundary, internal, {face: len(sites) for face, sites in by_face.items()}, selected


def elementary_even_star(spins: Sequence[int]) -> Poly:
    coefficients = []
    for degree in range(len(spins) + 1):
        if degree % 2:
            coefficients.append(0)
            continue
        coefficient = sum(
            math.prod(spins[index] for index in subset)
            for subset in itertools.combinations(range(len(spins)), degree)
        )
        coefficients.append(coefficient)
    return trim(coefficients)


def direct_signature(
    vertices: Sequence[Coordinate],
    edges: Sequence[Edge],
    eliminated_parity: int,
    boundary: Sequence[Coordinate],
) -> Signature:
    neighbours = neighbour_map(vertices, edges)
    eliminated = tuple(
        sorted(vertex for vertex in vertices if parity(vertex) == eliminated_parity)
    )
    retained = tuple(
        sorted(vertex for vertex in vertices if parity(vertex) != eliminated_parity)
    )
    retained_index = {vertex: index for index, vertex in enumerate(retained)}
    boundary_positions = tuple(retained_index[vertex] for vertex in boundary)
    spin_weights: list[Poly] = []

    for state in range(1 << len(retained)):
        signs = tuple(
            -1 if (state >> index) & 1 else 1
            for index in range(len(retained))
        )
        weight = (1,)
        for center in eliminated:
            factor = elementary_even_star(
                tuple(signs[retained_index[vertex]] for vertex in neighbours[center])
            )
            weight = multiply(weight, factor)
        spin_weights.append(weight)

    signature: list[Poly] = []
    for support in range(1 << len(boundary)):
        accumulator = (0,)
        for state, weight in enumerate(spin_weights):
            boundary_state = sum(
                ((state >> retained_position) & 1) << boundary_position
                for boundary_position, retained_position in enumerate(
                    boundary_positions
                )
            )
            sign = -1 if (boundary_state & support).bit_count() % 2 else 1
            accumulator = add(accumulator, scale(weight, sign))
        signature.append(exact_average(accumulator, 1 << len(retained)))
    return tuple(signature)


def direct_cycle_signature() -> Signature:
    spin_weights: list[Poly] = []
    for state in range(32):
        spins = tuple(-1 if (state >> index) & 1 else 1 for index in range(5))
        weight = (1,)
        for index in range(5):
            weight = multiply(
                weight,
                elementary_even_star((spins[index], spins[(index + 1) % 5])),
            )
        spin_weights.append(weight)
    signature = []
    for support in range(32):
        accumulator = (0,)
        for state, weight in enumerate(spin_weights):
            sign = -1 if (state & support).bit_count() % 2 else 1
            accumulator = add(accumulator, scale(weight, sign))
        signature.append(exact_average(accumulator, 32))
    return tuple(signature)


def brute_edge_signature(
    vertices: Sequence[Coordinate],
    edges: Sequence[Edge],
    eliminated_parity: int,
    boundary: Sequence[Coordinate],
) -> Signature:
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    boundary_index = {vertex: index for index, vertex in enumerate(boundary)}
    internal = {
        vertex
        for vertex in vertices
        if parity(vertex) != eliminated_parity and vertex not in boundary
    }
    constrained = {
        vertex
        for vertex in vertices
        if parity(vertex) == eliminated_parity or vertex in internal
    }
    accumulators: list[Poly] = [(0,)] * (1 << len(boundary))
    for subset in range(1 << len(edges)):
        incidence = 0
        for index, (left, right) in enumerate(edges):
            if (subset >> index) & 1:
                incidence ^= (1 << vertex_index[left]) | (1 << vertex_index[right])
        if any((incidence >> vertex_index[vertex]) & 1 for vertex in constrained):
            continue
        syndrome = sum(
            1 << boundary_index[vertex]
            for vertex in boundary
            if (incidence >> vertex_index[vertex]) & 1
        )
        monomial = (0,) * subset.bit_count() + (1,)
        accumulators[syndrome] = add(accumulators[syndrome], monomial)
    return tuple(accumulators)


def read_signature(record: dict[str, object]) -> Signature:
    entries = record["raw_entries"]
    dimension = int(record["boundary_legs"])
    expected_keys = {format(mask, f"0{dimension}b") for mask in range(1 << dimension)}
    if set(entries) != expected_keys:
        raise AssertionError("stored signature masks are incomplete")
    result: list[Poly] = []
    for mask in range(1 << dimension):
        row = entries[format(mask, f"0{dimension}b")]
        poly = trim(int(coefficient) for coefficient in row["coefficients_ascending"])
        terms = {
            str(degree): coefficient
            for degree, coefficient in enumerate(poly)
            if coefficient
        }
        degree = len(poly) - 1 if poly != (0,) else -1
        if row["nonzero_terms"] != terms or row["degree"] != degree:
            raise AssertionError("stored polynomial metadata disagrees with coefficients")
        result.append(poly)
    return tuple(result)


def unrank_permutation(rank: int, order: int = 5) -> tuple[int, ...]:
    if rank < 0 or rank >= math.factorial(order):
        raise ValueError("factoradic rank out of range")
    available = list(range(order))
    result = []
    remainder = rank
    for places in range(order, 0, -1):
        block = math.factorial(places - 1)
        index, remainder = divmod(remainder, block)
        result.append(available.pop(index))
    return tuple(result)


def physical_mask(order: Sequence[int], positions: Sequence[int]) -> int:
    mask = 0
    for position in positions:
        mask |= 1 << order[position]
    return mask


def four_residual(
    signature: Signature, order: Sequence[int], positions: Sequence[int]
) -> Poly:
    a, b, c, d = positions

    def pair(left: int, right: int) -> Poly:
        return signature[(1 << order[left]) | (1 << order[right])]

    pfaffian = add(
        subtract(
            multiply(pair(a, b), pair(c, d)),
            multiply(pair(a, c), pair(b, d)),
        ),
        multiply(pair(a, d), pair(b, c)),
    )
    target = multiply(
        signature[physical_mask(order, positions)], signature[0]
    )
    return subtract(target, pfaffian)


def divide_polynomials(dividend: Poly, divisor: Poly) -> tuple[Poly, Poly]:
    dividend_work = list(trim(dividend))
    divisor = trim(divisor)
    if divisor == (0,):
        raise ZeroDivisionError
    if len(dividend_work) < len(divisor):
        return (0,), trim(dividend_work)
    quotient = [0] * (len(dividend_work) - len(divisor) + 1)
    while len(dividend_work) >= len(divisor) and any(dividend_work):
        shift = len(dividend_work) - len(divisor)
        if dividend_work[-1] % divisor[-1]:
            break
        coefficient = dividend_work[-1] // divisor[-1]
        quotient[shift] = coefficient
        for index, divisor_coefficient in enumerate(divisor):
            dividend_work[shift + index] -= coefficient * divisor_coefficient
        while len(dividend_work) > 1 and dividend_work[-1] == 0:
            dividend_work.pop()
    return trim(quotient), trim(dividend_work)


def remove_factor(poly: Poly, factor: Poly) -> tuple[Poly, int]:
    exponent = 0
    quotient = poly
    while True:
        candidate, remainder = divide_polynomials(quotient, factor)
        if remainder != (0,):
            return quotient, exponent
        quotient = candidate
        exponent += 1


def clean_factor_projection(poly: Poly) -> tuple[str, dict[str, object], Poly]:
    v_power = next(index for index, coefficient in enumerate(poly) if coefficient)
    quotient = trim(poly[v_power:])
    quotient, one_minus = remove_factor(quotient, (1, -1))
    quotient, one_plus = remove_factor(quotient, (1, 1))
    quotient, one_plus_squared = remove_factor(quotient, (1, 0, 1))
    content = 0
    for coefficient in quotient:
        content = math.gcd(content, abs(coefficient))
    primitive = trim(coefficient // content for coefficient in quotient)
    if primitive[0] < 0:
        primitive = negate(primitive)
        content = -content
    if any(primitive[index] for index in range(1, len(primitive), 2)):
        raise AssertionError("residual remainder is not in ZZ[v^2]")
    coefficients_in_v_squared = list(primitive[::2])
    class_id = f"R_v{v_power}_Q{len(primitive) - 1}"
    projection = {
        "scalar": content,
        "v_power": v_power,
        "one_minus_v_power": one_minus,
        "one_plus_v_power": one_plus,
        "one_plus_v_squared_power": one_plus_squared,
        "positive_remainder_coefficients_in_v_squared": coefficients_in_v_squared,
    }
    return class_id, projection, primitive


def sympy_factorization(poly: Poly) -> dict[str, object]:
    expression = sum(
        sp.Integer(coefficient) * X**degree
        for degree, coefficient in enumerate(poly)
    )
    content, factors = sp.factor_list(sp.Poly(expression, X, domain=sp.ZZ))
    return {
        "content": int(content),
        "irreducible_factors": [
            {
                "coefficients_ascending": [
                    int(coefficient)
                    for coefficient in reversed(factor.all_coeffs())
                ],
                "exponent": int(exponent),
            }
            for factor, exponent in factors
        ],
        "reconstructs_residual": True,
    }


def target_symbolic_rebuild(
    signature: Signature,
) -> tuple[dict[str, object], dict[str, dict[str, object]], dict[tuple[int, ...], dict[str, object]]]:
    four_subsets = tuple(
        tuple(position for position in range(5) if position != omitted)
        for omitted in range(4, -1, -1)
    )
    orders = tuple(unrank_permutation(rank) for rank in range(120))
    residual_rows: dict[tuple[int, ...], tuple[Poly, ...]] = {
        order: tuple(four_residual(signature, order, subset) for subset in four_subsets)
        for order in orders
    }
    unique = set(
        residual
        for residuals in residual_rows.values()
        for residual in residuals
    )
    factors: dict[str, dict[str, object]] = {}
    residual_class: dict[Poly, str] = {}
    for residual in unique:
        class_id, projection, primitive = clean_factor_projection(residual)
        q_expression = sum(
            sp.Integer(coefficient) * X**degree
            for degree, coefficient in enumerate(primitive)
        )
        q_content, q_factors = sp.factor_list(
            sp.Poly(q_expression, X, domain=sp.ZZ)
        )
        projection["remainder_irreducible_over_ZZ"] = (
            int(q_content) == 1
            and len(q_factors) == 1
            and int(q_factors[0][1]) == 1
        )
        factors[class_id] = {
            "projection": projection,
            "residual": residual,
            "factorization_over_ZZ": sympy_factorization(residual),
        }
        residual_class[residual] = class_id

    first_histogram = Counter(
        residual_class[residual_rows[order][0]] for order in orders
    )
    all_histogram = Counter(
        residual_class[residual]
        for residuals in residual_rows.values()
        for residual in residuals
    )
    ordering_records = {
        order: {
            "first_nonzero_subset_positions": list(four_subsets[0]),
            "first_residual_class": residual_class[residual_rows[order][0]],
            "four_subset_class_sequence": [
                residual_class[residual] for residual in residual_rows[order]
            ],
        }
        for order in orders
    }
    summary = {
        "ordering_count": len(orders),
        "principal_pfaffians_checked": len(orders) * 16,
        "identically_zero_empty_and_pair_residuals": len(orders) * 11,
        "four_leg_residuals_checked": len(orders) * 5,
        "nonzero_four_leg_residuals": sum(
            residual != (0,)
            for residuals in residual_rows.values()
            for residual in residuals
        ),
        "passing_orderings": [
            list(order)
            for order, residuals in residual_rows.items()
            if all(residual == (0,) for residual in residuals)
        ],
        "first_residual_class_histogram": dict(sorted(first_histogram.items())),
        "all_four_subset_class_histogram": dict(sorted(all_histogram.items())),
    }
    return summary, factors, ordering_records


def pfaffian_from_matrix(
    matrix: Sequence[Sequence[Fraction]], positions: Sequence[int]
) -> Fraction:
    if not positions:
        return Fraction(1)
    if len(positions) == 2:
        return matrix[positions[0]][positions[1]]
    if len(positions) == 4:
        a, b, c, d = positions
        return (
            matrix[a][b] * matrix[c][d]
            - matrix[a][c] * matrix[b][d]
            + matrix[a][d] * matrix[b][c]
        )
    raise ValueError("five-leg verifier only expects sizes 0, 2, and 4")


def rational_rebuild(signature: Signature, value: Fraction) -> dict[str, object]:
    raw = tuple(evaluate(poly, value) for poly in signature)
    normalized = tuple(entry / raw[0] for entry in raw)
    even_position_masks = tuple(mask for mask in range(32) if mask.bit_count() % 2 == 0)
    matching = 0
    failing = 0
    passing: list[list[int]] = []
    failure_histogram: Counter[int] = Counter()

    for rank in range(120):
        order = unrank_permutation(rank)
        matrix = [[Fraction(0) for _ in range(5)] for _ in range(5)]
        for row in range(5):
            for column in range(row + 1, 5):
                pair = (1 << order[row]) | (1 << order[column])
                matrix[row][column] = normalized[pair]
                matrix[column][row] = -normalized[pair]
        order_failures = 0
        for position_mask in even_position_masks:
            positions = tuple(
                index for index in range(5) if (position_mask >> index) & 1
            )
            pfaffian = pfaffian_from_matrix(matrix, positions)
            expected = normalized[physical_mask(order, positions)]
            if pfaffian == expected:
                matching += 1
            else:
                failing += 1
                order_failures += 1
        failure_histogram[order_failures] += 1
        if order_failures == 0:
            passing.append(list(order))

    return {
        "v": fraction_record(value),
        "raw_empty_Walsh_entry": fraction_record(raw[0]),
        "normalized_signature": {
            format(mask, "05b"): fraction_record(entry)
            for mask, entry in enumerate(normalized)
        },
        "boundary_permutations_checked": 120,
        "even_principal_pfaffians_per_permutation": 16,
        "principal_pfaffians_checked": 1920,
        "matching_principal_pfaffians": matching,
        "failing_principal_pfaffians": failing,
        "passing_permutations": passing,
        "permutation_failure_count_histogram": {
            str(count): multiplicity
            for count, multiplicity in sorted(failure_histogram.items())
        },
    }


def dihedral_orders() -> set[tuple[int, ...]]:
    return {
        tuple((start + direction * offset) % 5 for offset in range(5))
        for start in range(5)
        for direction in (-1, 1)
    }


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text())
    meta = artifact["meta"]
    data = artifact["data"]

    expected_source_hashes = {
        "experiments/e247_boundary_matchgate_defect.py": file_sha256(PRODUCER),
        "tests/test_boundary_matchgate_defect.py": file_sha256(Path(__file__).resolve()),
        "proofs/boundary_matchgate_defect.md": file_sha256(PROOF),
        "experiments/e241_box_genus_bound.py": file_sha256(E241_PRODUCER),
    }
    check(
        "current source hash gate",
        meta["source_sha256"] == expected_source_hashes,
        "producer, verifier, proof, and inherited e241 source must match the artifact",
    )
    check(
        "certificate content hash",
        meta["certificate_sha256"] == canonical_sha256(data),
    )

    vertices, edges = clean_graph(SHAPE)
    geometry = data["geometry"]
    geometry_holds = (
        len(vertices) == 18
        and len(edges) == 33
        and len(set(edges)) == 33
        and {tuple(vertex) for vertex in geometry["vertices"]} == set(vertices)
        and {
            tuple(tuple(endpoint) for endpoint in edge)
            for edge in geometry["edges"]
        }
        == set(edges)
        and geometry["bipartition"]["sizes_even_odd"] == [9, 9]
        and all(parity(left) != parity(right) for left, right in edges)
    )
    check(
        "independent coordinate graph and 9+9 bipartition",
        geometry_holds,
        "pairwise Manhattan construction gives the same 18 vertices and 33 edges",
    )

    neighbours = neighbour_map(vertices, edges)
    rebuilt_signatures: list[Signature] = []
    stored_factor_classes = {
        record["class_id"]: record for record in data["residual_factor_classes"]
    }
    all_case_checks = True
    all_symbolic_checks = True
    all_rational_checks = True

    for eliminated_parity, case in zip((0, 1), data["checkerboard_cases"]):
        boundary, internal, face_counts, selected_face = clean_boundary(
            vertices, eliminated_parity
        )
        rebuilt = direct_signature(
            vertices, edges, eliminated_parity, boundary
        )
        stored = read_signature(case["boundary_signature"])
        rebuilt_signatures.append(rebuilt)

        expected_stars = {
            vertex: neighbours[vertex]
            for vertex in vertices
            if parity(vertex) == eliminated_parity
        }
        stored_stars = {
            tuple(row["center"]): tuple(tuple(vertex) for vertex in row["neighbour_order"])
            for row in case["local_stars"]["star_rows"]
        }
        all_case_checks &= (
            case["eliminated_color"]["parity_x_plus_y_plus_z_mod_2"]
            == eliminated_parity
            and case["retained_color"]["parity_x_plus_y_plus_z_mod_2"]
            == 1 - eliminated_parity
            and case["face_counts"]
            == {f"x={face}": count for face, count in face_counts.items()}
            and case["selected_face"]["coordinate"] == selected_face
            and tuple(tuple(vertex) for vertex in case["boundary_leg_order"])
            == boundary
            and tuple(tuple(vertex) for vertex in case["averaged_retained_coordinates"])
            == internal
            and stored_stars == expected_stars
            and Counter(map(len, stored_stars.values())) == {3: 4, 4: 4, 5: 1}
            and stored == rebuilt
            and rebuilt[0][0] == 1
            and all(
                rebuilt[mask] == (0,)
                for mask in range(32)
                if mask.bit_count() % 2
            )
        )

        symbolic, factors, ordering_records = target_symbolic_rebuild(rebuilt)
        stored_symbolic = case["symbolic_pfaffian_coverage"]
        all_symbolic_checks &= all(
            symbolic[key] == stored_symbolic[key]
            for key in (
                "ordering_count",
                "principal_pfaffians_checked",
                "identically_zero_empty_and_pair_residuals",
                "four_leg_residuals_checked",
                "nonzero_four_leg_residuals",
                "passing_orderings",
                "first_residual_class_histogram",
                "all_four_subset_class_histogram",
            )
        )
        stored_ordering_records = {
            tuple(row["permutation"]): row
            for row in stored_symbolic["ordering_records"]
        }
        all_symbolic_checks &= set(stored_ordering_records) == set(ordering_records)
        for order, rebuilt_row in ordering_records.items():
            stored_row = stored_ordering_records[order]
            all_symbolic_checks &= (
                stored_row["first_nonzero_subset_positions"]
                == rebuilt_row["first_nonzero_subset_positions"]
                and stored_row["first_residual_class"]
                == rebuilt_row["first_residual_class"]
                and stored_row["four_subset_class_sequence"]
                == rebuilt_row["four_subset_class_sequence"]
                and stored_row["coordinate_order"]
                == [list(boundary[index]) for index in order]
            )

        all_symbolic_checks &= set(factors) == set(EXPECTED_FACTOR_DATA)
        for class_id, rebuilt_factor in factors.items():
            expected = EXPECTED_FACTOR_DATA[class_id]
            projection = dict(rebuilt_factor["projection"])
            irreducible = projection.pop("remainder_irreducible_over_ZZ")
            expected_projection = {
                key: expected[key]
                for key in (
                    "scalar",
                    "v_power",
                    "one_minus_v_power",
                    "one_plus_v_power",
                    "one_plus_v_squared_power",
                    "positive_remainder_coefficients_in_v_squared",
                )
            }
            stored_factor = stored_factor_classes[class_id]
            all_symbolic_checks &= (
                projection == expected_projection
                and irreducible is True
                and stored_factor["positive_interval_factorization"]
                == expected_projection
                and stored_factor["first_ordering_count"]
                == expected["first_ordering_count"]
                and stored_factor["all_four_subset_occurrences"]
                == expected["all_four_subset_occurrences"]
                and tuple(stored_factor["residual"]["coefficients_ascending"])
                == rebuilt_factor["residual"]
                and stored_factor["factorization_over_ZZ"]
                == rebuilt_factor["factorization_over_ZZ"]
                and stored_factor["coefficient_positivity_certificate"]
                == {
                    "positive_scalar": True,
                    "all_remainder_coefficients_positive": True,
                    "conclusion": "strictly positive for 0<v<1",
                }
            )
            for point in AUDIT_POINTS:
                residual_value = evaluate(rebuilt_factor["residual"], point)
                expected_audit = {
                    "raw_residual_numerator": fraction_record(residual_value),
                    "normalized_residual": fraction_record(
                        residual_value / evaluate(rebuilt[0], point) ** 2
                    ),
                }
                all_symbolic_checks &= (
                    stored_factor["audit_points"][str(point)]
                    == expected_audit
                )

        for point in AUDIT_POINTS:
            rebuilt_audit = rational_rebuild(rebuilt, point)
            all_rational_checks &= (
                case["rational_controls"][str(point)] == rebuilt_audit
                and rebuilt_audit["matching_principal_pfaffians"] == 1320
                and rebuilt_audit["failing_principal_pfaffians"] == 600
                and rebuilt_audit["passing_permutations"] == []
            )

    check(
        "independent direct-spin target signatures",
        all_case_checks and rebuilt_signatures[0] == rebuilt_signatures[1],
        "both reflected color choices are rebuilt from 512 retained-spin assignments without e247",
    )
    check(
        "independent all-order symbolic residuals and ZZ factors",
        all_symbolic_checks,
        "factoradic 120-order enumeration reproduces all 600 positive residuals and four irreducible factor classes",
    )
    check(
        "independent exact rational Pfaffian controls",
        all_rational_checks,
        "v=1/3 and v=1/2 each give five exact failures in every target ordering",
    )

    control_vertices, control_edges = clean_graph((2, 2, 2))
    control_boundary = tuple(
        sorted(
            (
                vertex
                for vertex in control_vertices
                if parity(vertex) == 1 and vertex[0] == 1
            ),
            key=lambda vertex: 2 * vertex[1] + vertex[2],
        )
    )
    control_direct = direct_signature(
        control_vertices, control_edges, 0, control_boundary
    )
    control_edgesum = brute_edge_signature(
        control_vertices, control_edges, 0, control_boundary
    )
    small = data["small_control"]
    check(
        "independent small spin-sum and edge-subset control",
        len(control_edges) == 12
        and control_direct == control_edgesum
        and read_signature(small["direct_spin_signature"]) == control_direct
        and read_signature(small["edge_boundary_signature"]) == control_edgesum
        and small["edge_subsets"] == 4096
        and small["retained_spin_assignments"] == 16,
    )

    cycle = direct_cycle_signature()
    cycle_stored = data["planar_alternating_cycle_control"]
    cycle_rows = {}
    four_subsets = tuple(
        tuple(position for position in range(5) if position != omitted)
        for omitted in range(4, -1, -1)
    )
    cycle_zero = 120 * 11
    cycle_nonzero = 0
    for rank in range(120):
        order = unrank_permutation(rank)
        residuals = tuple(four_residual(cycle, order, subset) for subset in four_subsets)
        cycle_rows[order] = residuals
        cycle_zero += sum(residual == (0,) for residual in residuals)
        cycle_nonzero += sum(residual != (0,) for residual in residuals)
    cycle_passing = {
        order for order, residuals in cycle_rows.items() if all(residual == (0,) for residual in residuals)
    }
    planar_checks = (
        read_signature(cycle_stored["signature"]) == cycle
        and cycle_zero == 1520
        and cycle_nonzero == 400
        and cycle_passing == dihedral_orders()
        and {
            tuple(order)
            for order in cycle_stored["symbolic_coverage"]["passing_permutations"]
        }
        == cycle_passing
        and cycle_stored["symbolic_coverage"]["principal_pfaffians_checked"]
        == 1920
    )
    for point in AUDIT_POINTS:
        planar_checks &= (
            cycle_stored["rational_controls"][str(point)]
            == rational_rebuild(cycle, point)
        )
    check(
        "independent planar alternating-cycle positive control",
        planar_checks,
        "direct spin sum passes exactly the ten cyclic/reversed boundary orders under the same Pfaffian checker",
    )

    topology = data["topology_provenance"]
    check(
        "e241 topology provenance and terminology boundary",
        topology["shape"] == [2, 3, 3]
        and topology["vertices"] == 18
        and topology["edges"] == 33
        and topology["e241_genus_bound_numerator"] == 1
        and topology["e241_orientable_genus_lower_bound"] == 1
        and "no exact genus is asserted" in topology["conclusion"]
        and topology["inherited_source"]
        == {
            "path": "experiments/e241_box_genus_bound.py",
            "sha256": expected_source_hashes["experiments/e241_box_genus_bound.py"],
            "binding": "current source hash; no e241 result artifact is imported",
        }
        and "not called spin structures" in topology["terminology_boundary"],
    )

    scope = data["scope"]
    check(
        "scope remains fixed-boundary and non-global",
        data["parameter"]["domain"] == "0<v<1"
        and "fixed five-leg boundary contraction" in scope["precise_boundary"]
        and any("arbitrary hidden auxiliaries" in statement for statement in scope["not_proved"])
        and any("fermionic-SWAP" in statement for statement in scope["not_proved"])
        and any("number of Pfaffians" in statement for statement in scope["not_proved"])
        and any("all lengths L" in statement for statement in scope["not_proved"])
        and any("three-dimensional Ising solution" in statement for statement in scope["not_proved"]),
    )

    expected_checks = {
        "source_hashes_bound",
        "explicit_open_box_geometry_and_bipartition",
        "larger_opposite_color_face_selected",
        "exact_all_degree_star_contraction",
        "small_direct_spin_sum_equals_edge_boundary_enumeration",
        "complete_symbolic_ordering_and_pfaffian_coverage",
        "all_v_positive_obstruction",
        "exact_one_third_and_one_half_controls",
        "planar_alternating_cycle_positive_control",
        "e241_nonplanarity_provenance",
        "scope_boundary",
        "benchmark_not_used",
        "declared_resource_limits",
    }
    check(
        "producer checks are hard gates",
        {row["name"] for row in artifact["checks"]} == expected_checks
        and len(artifact["checks"]) == len(expected_checks)
        and all(row["passed"] is True for row in artifact["checks"]),
    )
    check(
        "producer metadata and resource gates",
        meta["producer"] == "experiments/e247_boundary_matchgate_defect.py"
        and meta["verifier"] == "tests/test_boundary_matchgate_defect.py"
        and meta["proof"] == "proofs/boundary_matchgate_defect.md"
        and meta["artifact"]
        == "results/integrability/boundary_matchgate_defect.json"
        and meta["arithmetic"]
        == "exact integers, fractions, and polynomials over ZZ[v]"
        and meta["benchmark_used"] is False
        and meta["process_cpu_seconds"] < meta["process_cpu_budget_seconds"]
        and meta["process_cpu_budget_seconds"] == CPU_BUDGET_SECONDS
        and meta["peak_rss_bytes"] < meta["rss_limit_bytes"]
        and meta["rss_limit_bytes"] == RSS_LIMIT_BYTES
        and meta["peak_rss_measurement"] == peak_rss_measurement(),
    )

    elapsed = time.process_time() - started
    peak_rss = peak_rss_bytes()
    check(
        "verifier resource walls",
        elapsed < CPU_BUDGET_SECONDS and peak_rss < RSS_LIMIT_BYTES,
        (
            f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; "
            f"peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes via {peak_rss_measurement()}"
        ),
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: open-2x3x3 fixed-boundary matchgate defect")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
