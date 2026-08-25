#!/usr/bin/env python3
"""Exact five-leg boundary matchgate defect for the open 2x3x3 box.

The two checkerboard eliminations are performed from explicit coordinates.  A
five-spin subset of the opposite color is retained on the unique 3x3 face on
which that color has five sites, and the other four retained-color spins are
averaged.  All signatures and Pfaffian residuals are polynomials over ZZ[v]
or exact Fraction specializations.
"""

from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "boundary_matchgate_defect.json"
VERIFIER = ROOT / "tests" / "test_boundary_matchgate_defect.py"
PROOF = ROOT / "proofs" / "boundary_matchgate_defect.md"
E241_PRODUCER = ROOT / "experiments" / "e241_box_genus_bound.py"
SHAPE = (2, 3, 3)
BOUNDARY_AXIS = 0
BOUNDARY_LEGS = 5
AUDIT_POINTS = (Fraction(1, 3), Fraction(1, 2))
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
V = sp.symbols("v")

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


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def trim(poly: Iterable[int]) -> Poly:
    coefficients = list(poly)
    if not coefficients:
        return (0,)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(coefficients)


def poly_add(left: Poly, right: Poly) -> Poly:
    width = max(len(left), len(right))
    return trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(width)
    )


def poly_negate(poly: Poly) -> Poly:
    return tuple(-coefficient for coefficient in poly)


def poly_subtract(left: Poly, right: Poly) -> Poly:
    return poly_add(left, poly_negate(right))


def poly_multiply(left: Poly, right: Poly) -> Poly:
    out = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            out[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return trim(out)


def poly_scale(poly: Poly, scalar: int) -> Poly:
    return trim(scalar * coefficient for coefficient in poly)


def poly_shift(poly: Poly, degree: int) -> Poly:
    if degree < 0:
        raise ValueError("polynomial shift must be nonnegative")
    if poly == (0,):
        return poly
    return (0,) * degree + poly


def poly_power(poly: Poly, exponent: int) -> Poly:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    out = (1,)
    base = poly
    power = exponent
    while power:
        if power & 1:
            out = poly_multiply(out, base)
        base = poly_multiply(base, base)
        power //= 2
    return out


def poly_evaluate(poly: Poly, value: Fraction) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(poly):
        total = total * value + coefficient
    return total


def poly_exact_divide_scalar(poly: Poly, denominator: int) -> Poly:
    if denominator <= 0 or any(coefficient % denominator for coefficient in poly):
        raise ArithmeticError("nonintegral polynomial average")
    return trim(coefficient // denominator for coefficient in poly)


def poly_record(poly: Poly) -> dict[str, object]:
    return {
        "degree": len(poly) - 1 if poly != (0,) else -1,
        "coefficients_ascending": list(poly),
        "nonzero_terms": {
            str(degree): coefficient
            for degree, coefficient in enumerate(poly)
            if coefficient
        },
    }


def fraction_record(value: int | Fraction) -> dict[str, int]:
    rational = Fraction(value)
    return {
        "numerator": rational.numerator,
        "denominator": rational.denominator,
    }


def coordinate_record(coordinate: Coordinate) -> list[int]:
    return list(coordinate)


def graph(shape: Sequence[int]) -> tuple[tuple[Coordinate, ...], tuple[Edge, ...]]:
    if len(shape) != 3 or any(type(side) is not int or side < 1 for side in shape):
        raise ValueError("shape must contain three positive integers")
    vertices = tuple(
        itertools.product(*(range(side) for side in shape))
    )
    edges: list[Edge] = []
    for vertex in vertices:
        for axis, side in enumerate(shape):
            if vertex[axis] + 1 >= side:
                continue
            neighbour = list(vertex)
            neighbour[axis] += 1
            edges.append((vertex, tuple(neighbour)))
    return vertices, tuple(edges)


def parity(coordinate: Coordinate) -> int:
    return sum(coordinate) % 2


def adjacency(
    vertices: Sequence[Coordinate], edges: Sequence[Edge]
) -> dict[Coordinate, tuple[Coordinate, ...]]:
    mutable = {vertex: [] for vertex in vertices}
    for left, right in edges:
        mutable[left].append(right)
        mutable[right].append(left)
    return {
        vertex: tuple(sorted(neighbours))
        for vertex, neighbours in mutable.items()
    }


def boundary_selection(
    vertices: Sequence[Coordinate], eliminated_parity: int
) -> dict[str, object]:
    retained_parity = 1 - eliminated_parity
    retained = tuple(
        vertex for vertex in vertices if parity(vertex) == retained_parity
    )
    face_subsets = {
        face: tuple(
            sorted(
                (
                    vertex
                    for vertex in retained
                    if vertex[BOUNDARY_AXIS] == face
                ),
                key=lambda vertex: (vertex[1], vertex[2]),
            )
        )
        for face in range(SHAPE[BOUNDARY_AXIS])
    }
    maximum = max(map(len, face_subsets.values()))
    maximizing_faces = [
        face for face, sites in face_subsets.items() if len(sites) == maximum
    ]
    if len(maximizing_faces) != 1:
        raise AssertionError("the retained color must have a unique larger face")
    selected_face = maximizing_faces[0]
    boundary = face_subsets[selected_face]
    internal = tuple(vertex for vertex in retained if vertex not in boundary)
    return {
        "retained_parity": retained_parity,
        "retained": retained,
        "face_subsets": face_subsets,
        "selected_face": selected_face,
        "boundary": boundary,
        "internal": internal,
    }


def star_dp_signature(
    vertices: Sequence[Coordinate],
    edges: Sequence[Edge],
    eliminated_parity: int,
    boundary: Sequence[Coordinate],
    internal: Sequence[Coordinate],
) -> tuple[Signature, dict[str, object]]:
    neighbours = adjacency(vertices, edges)
    eliminated = tuple(
        vertex for vertex in vertices if parity(vertex) == eliminated_parity
    )
    retained = tuple(
        vertex for vertex in vertices if parity(vertex) != eliminated_parity
    )
    retained_index = {vertex: index for index, vertex in enumerate(retained)}
    states: dict[int, Poly] = {0: (1,)}
    maximum_state_count = 1
    local_option_count = 0
    star_rows: list[dict[str, object]] = []

    for center in eliminated:
        ordered_neighbours = neighbours[center]
        options: list[tuple[int, int]] = []
        for subset in range(1 << len(ordered_neighbours)):
            if subset.bit_count() % 2:
                continue
            retained_mask = 0
            for local_index, vertex in enumerate(ordered_neighbours):
                if (subset >> local_index) & 1:
                    retained_mask ^= 1 << retained_index[vertex]
            options.append((retained_mask, subset.bit_count()))
        local_option_count += len(options)
        next_states: dict[int, Poly] = {}
        for state_mask, state_poly in states.items():
            for option_mask, degree in options:
                target = state_mask ^ option_mask
                contribution = poly_shift(state_poly, degree)
                next_states[target] = poly_add(
                    next_states.get(target, (0,)), contribution
                )
        states = next_states
        maximum_state_count = max(maximum_state_count, len(states))
        star_rows.append(
            {
                "center": coordinate_record(center),
                "degree": len(ordered_neighbours),
                "neighbour_order": [
                    coordinate_record(vertex) for vertex in ordered_neighbours
                ],
                "even_local_subsets": len(options),
            }
        )

    internal_mask = sum(1 << retained_index[vertex] for vertex in internal)
    boundary_index = {vertex: index for index, vertex in enumerate(boundary)}
    signature: list[Poly] = [(0,)] * (1 << len(boundary))
    for state_mask, state_poly in states.items():
        if state_mask & internal_mask:
            continue
        boundary_mask = sum(
            1 << boundary_index[vertex]
            for vertex in boundary
            if (state_mask >> retained_index[vertex]) & 1
        )
        signature[boundary_mask] = poly_add(
            signature[boundary_mask], state_poly
        )
    return tuple(signature), {
        "star_rows": star_rows,
        "sum_of_local_even_subset_counts": local_option_count,
        "maximum_parity_states": maximum_state_count,
        "final_parity_states": len(states),
    }


def star_factor_poly(spins: Sequence[int]) -> Poly:
    plus = (1,)
    minus = (1,)
    for spin in spins:
        if spin not in (-1, 1):
            raise ValueError("spins must be signs")
        plus = poly_multiply(plus, (1, spin))
        minus = poly_multiply(minus, (1, -spin))
    total = poly_add(plus, minus)
    return poly_exact_divide_scalar(total, 2)


def direct_spin_signature(
    vertices: Sequence[Coordinate],
    edges: Sequence[Edge],
    eliminated_parity: int,
    boundary: Sequence[Coordinate],
) -> Signature:
    neighbours = adjacency(vertices, edges)
    eliminated = tuple(
        vertex for vertex in vertices if parity(vertex) == eliminated_parity
    )
    retained = tuple(
        vertex for vertex in vertices if parity(vertex) != eliminated_parity
    )
    retained_index = {vertex: index for index, vertex in enumerate(retained)}
    boundary_positions = tuple(retained_index[vertex] for vertex in boundary)
    accumulators: list[Poly] = [(0,)] * (1 << len(boundary))

    for spin_state in range(1 << len(retained)):
        spins = {
            vertex: (-1 if (spin_state >> index) & 1 else 1)
            for vertex, index in retained_index.items()
        }
        weight = (1,)
        for center in eliminated:
            weight = poly_multiply(
                weight,
                star_factor_poly(
                    tuple(spins[vertex] for vertex in neighbours[center])
                ),
            )
        boundary_state = sum(
            ((spin_state >> retained_position) & 1) << boundary_position
            for boundary_position, retained_position in enumerate(
                boundary_positions
            )
        )
        for support in range(1 << len(boundary)):
            character = -1 if (boundary_state & support).bit_count() % 2 else 1
            accumulators[support] = poly_add(
                accumulators[support], poly_scale(weight, character)
            )

    denominator = 1 << len(retained)
    return tuple(
        poly_exact_divide_scalar(poly, denominator)
        for poly in accumulators
    )


def edge_subset_signature(
    vertices: Sequence[Coordinate],
    edges: Sequence[Edge],
    eliminated_parity: int,
    boundary: Sequence[Coordinate],
) -> Signature:
    retained = tuple(
        vertex for vertex in vertices if parity(vertex) != eliminated_parity
    )
    internal = tuple(vertex for vertex in retained if vertex not in boundary)
    constrained = {
        vertex
        for vertex in vertices
        if parity(vertex) == eliminated_parity or vertex in internal
    }
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    boundary_index = {vertex: index for index, vertex in enumerate(boundary)}
    signature: list[Poly] = [(0,)] * (1 << len(boundary))

    for edge_subset in range(1 << len(edges)):
        degree_parities = 0
        for edge_index, (left, right) in enumerate(edges):
            if (edge_subset >> edge_index) & 1:
                degree_parities ^= 1 << vertex_index[left]
                degree_parities ^= 1 << vertex_index[right]
        if any(
            (degree_parities >> vertex_index[vertex]) & 1
            for vertex in constrained
        ):
            continue
        boundary_mask = sum(
            1 << boundary_index[vertex]
            for vertex in boundary
            if (degree_parities >> vertex_index[vertex]) & 1
        )
        monomial = poly_shift((1,), edge_subset.bit_count())
        signature[boundary_mask] = poly_add(
            signature[boundary_mask], monomial
        )
    return tuple(signature)


def alternating_cycle_signature(order: int) -> Signature:
    if order < 3:
        raise ValueError("alternating cycle needs at least three boundary sites")
    signature: list[Poly] = [(0,)] * (1 << order)
    for selected_stars in range(1 << order):
        boundary_mask = 0
        for star in range(order):
            if (selected_stars >> star) & 1:
                boundary_mask ^= 1 << star
                boundary_mask ^= 1 << ((star + 1) % order)
        signature[boundary_mask] = poly_add(
            signature[boundary_mask],
            poly_shift((1,), 2 * selected_stars.bit_count()),
        )
    return tuple(signature)


def physical_mask(order: Sequence[int], positions: Sequence[int]) -> int:
    return sum(1 << order[position] for position in positions)


def pfaffian_poly(matrix: Sequence[Sequence[Poly]]) -> Poly:
    order = len(matrix)
    if order == 0:
        return (1,)
    if order % 2:
        return (0,)
    total = (0,)
    for partner in range(1, order):
        remaining = [
            index for index in range(1, order) if index != partner
        ]
        minor = [
            [matrix[row][column] for column in remaining]
            for row in remaining
        ]
        term = poly_multiply(matrix[0][partner], pfaffian_poly(minor))
        total = poly_add(
            total, term if partner % 2 else poly_negate(term)
        )
    return total


def raw_pair_matrix(signature: Signature, order: Sequence[int]) -> list[list[Poly]]:
    matrix = [[(0,) for _ in order] for _ in order]
    for row in range(len(order)):
        for column in range(row + 1, len(order)):
            pair = (1 << order[row]) | (1 << order[column])
            matrix[row][column] = signature[pair]
            matrix[column][row] = poly_negate(signature[pair])
    return matrix


def principal_residual_numerator(
    signature: Signature, order: Sequence[int], positions: Sequence[int]
) -> Poly:
    if len(positions) % 2:
        raise ValueError("only even principal subsets have Pfaffians")
    if not positions:
        return (0,)
    matrix = raw_pair_matrix(signature, order)
    principal = [
        [matrix[row][column] for column in positions] for row in positions
    ]
    raw_pfaffian = pfaffian_poly(principal)
    support = physical_mask(order, positions)
    half_order = len(positions) // 2
    target = poly_multiply(
        signature[support], poly_power(signature[0], half_order - 1)
    )
    return poly_subtract(target, raw_pfaffian)


def sympy_poly(poly: Poly) -> sp.Poly:
    expression = sum(
        sp.Integer(coefficient) * V**degree
        for degree, coefficient in enumerate(poly)
    )
    return sp.Poly(expression, V, domain=sp.ZZ)


def coefficients_ascending(poly: sp.Poly) -> Poly:
    return trim(int(coefficient) for coefficient in reversed(poly.all_coeffs()))


def residual_factor_record(poly: Poly) -> dict[str, object]:
    content, factors = sp.factor_list(sympy_poly(poly))
    factor_rows: list[dict[str, object]] = []
    powers = {
        "v": 0,
        "v_minus_one": 0,
        "v_plus_one": 0,
        "one_plus_v_squared": 0,
    }
    remainders: list[tuple[Poly, int]] = []
    reconstruction = (int(content),)

    for factor, exponent in factors:
        coefficients = coefficients_ascending(factor)
        exponent = int(exponent)
        factor_rows.append(
            {
                "coefficients_ascending": list(coefficients),
                "exponent": exponent,
            }
        )
        reconstruction = poly_multiply(
            reconstruction, poly_power(coefficients, exponent)
        )
        if coefficients == (0, 1):
            powers["v"] = exponent
        elif coefficients == (-1, 1):
            powers["v_minus_one"] = exponent
        elif coefficients == (1, 1):
            powers["v_plus_one"] = exponent
        elif coefficients == (1, 0, 1):
            powers["one_plus_v_squared"] = exponent
        else:
            remainders.append((coefficients, exponent))

    if reconstruction != poly:
        raise AssertionError("ZZ factorization did not reconstruct the residual")
    if len(remainders) != 1 or remainders[0][1] != 1:
        raise AssertionError("expected one residual factor beyond the elementary factors")
    remainder = remainders[0][0]
    if any(remainder[index] for index in range(1, len(remainder), 2)):
        raise AssertionError("positive remainder must be a polynomial in v^2")
    remainder_in_v_squared = list(remainder[::2])
    positive_scalar = int(content) * (-1) ** powers["v_minus_one"]
    class_id = f"R_v{powers['v']}_Q{len(remainder) - 1}"
    interval_projection = {
        "scalar": positive_scalar,
        "v_power": powers["v"],
        "one_minus_v_power": powers["v_minus_one"],
        "one_plus_v_power": powers["v_plus_one"],
        "one_plus_v_squared_power": powers["one_plus_v_squared"],
        "positive_remainder_coefficients_in_v_squared": remainder_in_v_squared,
    }
    return {
        "class_id": class_id,
        "residual": poly_record(poly),
        "factorization_over_ZZ": {
            "content": int(content),
            "irreducible_factors": factor_rows,
            "reconstructs_residual": True,
        },
        "positive_interval_factorization": interval_projection,
        "coefficient_positivity_certificate": {
            "positive_scalar": positive_scalar > 0,
            "all_remainder_coefficients_positive": all(
                coefficient > 0 for coefficient in remainder_in_v_squared
            ),
            "conclusion": "strictly positive for 0<v<1",
        },
    }


def symbolic_coverage(
    signature: Signature,
) -> tuple[list[dict[str, object]], list[Poly], int, int]:
    even_subsets = [
        subset
        for size in range(0, BOUNDARY_LEGS + 1, 2)
        for subset in itertools.combinations(range(BOUNDARY_LEGS), size)
    ]
    orders = tuple(itertools.permutations(range(BOUNDARY_LEGS)))
    rows: list[dict[str, object]] = []
    four_residuals: list[Poly] = []
    zero_residuals = 0
    nonzero_residuals = 0

    for order in orders:
        residuals = [
            principal_residual_numerator(signature, order, subset)
            for subset in even_subsets
        ]
        zero_residuals += sum(residual == (0,) for residual in residuals)
        nonzero_residuals += sum(residual != (0,) for residual in residuals)
        four = residuals[-math.comb(BOUNDARY_LEGS, 4) :]
        four_residuals.extend(four)
        rows.append(
            {
                "permutation": list(order),
                "all_even_principal_residuals_zero": all(
                    residual == (0,) for residual in residuals
                ),
                "four_residuals": four,
            }
        )
    return rows, four_residuals, zero_residuals, nonzero_residuals


def target_symbolic_certificate(
    signature: Signature, boundary: Sequence[Coordinate]
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows, four_residuals, zero_count, nonzero_count = symbolic_coverage(signature)
    first_residuals = [row["four_residuals"][0] for row in rows]
    unique_residuals = sorted(
        set(four_residuals),
        key=lambda poly: (
            next(index for index, coefficient in enumerate(poly) if coefficient),
            len(poly),
            poly,
        ),
    )
    factor_records = [residual_factor_record(poly) for poly in unique_residuals]
    residual_to_class = {
        tuple(record["residual"]["coefficients_ascending"]): record["class_id"]
        for record in factor_records
    }
    first_histogram = Counter(
        residual_to_class[residual] for residual in first_residuals
    )
    all_histogram = Counter(
        residual_to_class[residual] for residual in four_residuals
    )

    ordering_records = []
    for row in rows:
        order = tuple(row["permutation"])
        class_sequence = [
            residual_to_class[residual] for residual in row["four_residuals"]
        ]
        first_index = next(
            index
            for index, residual in enumerate(row["four_residuals"])
            if residual != (0,)
        )
        first_subset = list(
            itertools.combinations(range(BOUNDARY_LEGS), 4)
        )[first_index]
        ordering_records.append(
            {
                "permutation": list(order),
                "coordinate_order": [
                    coordinate_record(boundary[index]) for index in order
                ],
                "first_nonzero_subset_positions": list(first_subset),
                "first_residual_class": class_sequence[first_index],
                "four_subset_class_sequence": class_sequence,
            }
        )

    for record in factor_records:
        class_id = str(record["class_id"])
        record["first_ordering_count"] = first_histogram[class_id]
        record["all_four_subset_occurrences"] = all_histogram[class_id]
        record["audit_points"] = {
            str(point): {
                "raw_residual_numerator": fraction_record(
                    poly_evaluate(
                        tuple(record["residual"]["coefficients_ascending"]),
                        point,
                    )
                ),
                "normalized_residual": fraction_record(
                    poly_evaluate(
                        tuple(record["residual"]["coefficients_ascending"]),
                        point,
                    )
                    / poly_evaluate(signature[0], point) ** 2
                ),
            }
            for point in AUDIT_POINTS
        }

    summary = {
        "pair_matrix_convention": (
            "for permutation pi, M[i,j]=B_{pi_i,pi_j}/B_empty for i<j "
            "and M[j,i]=-M[i,j]"
        ),
        "cleared_residual_convention": (
            "R(pi,S)=B_S*B_empty^(|S|/2-1)-Pf(B_pair[S]); "
            "the normalized residual is R/B_empty^(|S|/2)"
        ),
        "ordering_count": len(rows),
        "even_principal_subsets_per_ordering": 1 + math.comb(5, 2) + math.comb(5, 4),
        "principal_pfaffians_checked": len(rows) * 16,
        "identically_zero_empty_and_pair_residuals": zero_count,
        "nonzero_four_leg_residuals": nonzero_count,
        "four_leg_residuals_checked": len(four_residuals),
        "passing_orderings": [
            row["permutation"]
            for row in rows
            if row["all_even_principal_residuals_zero"]
        ],
        "first_subset_positions_for_every_ordering": [0, 1, 2, 3],
        "first_residual_class_histogram": dict(sorted(first_histogram.items())),
        "all_four_subset_class_histogram": dict(sorted(all_histogram.items())),
        "ordering_records": ordering_records,
        "all_v_conclusion": (
            "every cleared four-leg residual is strictly positive for 0<v<1; "
            "therefore no boundary-leg permutation is sub-Pfaffian"
        ),
    }
    return summary, factor_records


def pfaffian_fraction(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    order = len(matrix)
    if order == 0:
        return Fraction(1)
    if order % 2:
        return Fraction(0)
    total = Fraction(0)
    for partner in range(1, order):
        remaining = [
            index for index in range(1, order) if index != partner
        ]
        minor = [
            [matrix[row][column] for column in remaining]
            for row in remaining
        ]
        term = matrix[0][partner] * pfaffian_fraction(minor)
        total += term if partner % 2 else -term
    return total


def rational_pfaffian_audit(signature: Signature, value: Fraction) -> dict[str, object]:
    raw_values = [poly_evaluate(poly, value) for poly in signature]
    empty = raw_values[0]
    if empty == 0:
        raise ZeroDivisionError("empty Walsh entry vanished at an audit point")
    normalized = [entry / empty for entry in raw_values]
    even_subsets = [
        subset
        for size in range(0, BOUNDARY_LEGS + 1, 2)
        for subset in itertools.combinations(range(BOUNDARY_LEGS), size)
    ]
    matching = 0
    failing = 0
    passing_orders: list[list[int]] = []
    failure_histogram: Counter[int] = Counter()

    for order in itertools.permutations(range(BOUNDARY_LEGS)):
        matrix = [[Fraction(0) for _ in order] for _ in order]
        for row in range(BOUNDARY_LEGS):
            for column in range(row + 1, BOUNDARY_LEGS):
                pair = (1 << order[row]) | (1 << order[column])
                matrix[row][column] = normalized[pair]
                matrix[column][row] = -normalized[pair]
        order_failures = 0
        for subset in even_subsets:
            principal = [
                [matrix[row][column] for column in subset] for row in subset
            ]
            pfaffian = pfaffian_fraction(principal)
            expected = normalized[physical_mask(order, subset)]
            if pfaffian == expected:
                matching += 1
            else:
                failing += 1
                order_failures += 1
        failure_histogram[order_failures] += 1
        if order_failures == 0:
            passing_orders.append(list(order))

    return {
        "v": fraction_record(value),
        "raw_empty_Walsh_entry": fraction_record(empty),
        "normalized_signature": {
            format(mask, "05b"): fraction_record(entry)
            for mask, entry in enumerate(normalized)
        },
        "boundary_permutations_checked": math.factorial(BOUNDARY_LEGS),
        "even_principal_pfaffians_per_permutation": len(even_subsets),
        "principal_pfaffians_checked": math.factorial(BOUNDARY_LEGS)
        * len(even_subsets),
        "matching_principal_pfaffians": matching,
        "failing_principal_pfaffians": failing,
        "passing_permutations": passing_orders,
        "permutation_failure_count_histogram": {
            str(count): multiplicity
            for count, multiplicity in sorted(failure_histogram.items())
        },
    }


def signature_record(signature: Signature) -> dict[str, object]:
    boundary_legs = (len(signature) - 1).bit_length()
    return {
        "boundary_legs": boundary_legs,
        "Walsh_basis": "chi_A(sigma)=product_(i in A) sigma_i",
        "raw_entries": {
            format(mask, f"0{boundary_legs}b"): poly_record(poly)
            for mask, poly in enumerate(signature)
        },
        "normalization": "g_A(v)=B_A(v)/B_empty(v)",
        "empty_entry": poly_record(signature[0]),
    }


def expected_factor_projection(record: dict[str, object]) -> dict[str, object]:
    projection = dict(record["positive_interval_factorization"])
    projection["first_ordering_count"] = record["first_ordering_count"]
    projection["all_four_subset_occurrences"] = record[
        "all_four_subset_occurrences"
    ]
    return projection


def dihedral_orders(order: int) -> set[tuple[int, ...]]:
    return {
        tuple((start + direction * offset) % order for offset in range(order))
        for start in range(order)
        for direction in (-1, 1)
    }


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    source_paths = (
        Path(__file__).resolve(),
        VERIFIER,
        PROOF,
        E241_PRODUCER,
    )
    source_hashes = {
        str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
    }
    add_check(
        checks,
        "source_hashes_bound",
        set(source_hashes)
        == {
            "experiments/e247_boundary_matchgate_defect.py",
            "tests/test_boundary_matchgate_defect.py",
            "proofs/boundary_matchgate_defect.md",
            "experiments/e241_box_genus_bound.py",
        }
        and all(len(digest) == 64 for digest in source_hashes.values()),
        "producer, clean-room verifier, scoped proof, and inherited e241 topology source",
    )

    vertices, edges = graph(SHAPE)
    parity_classes = {
        color: tuple(vertex for vertex in vertices if parity(vertex) == color)
        for color in (0, 1)
    }
    geometry_holds = (
        len(vertices) == 18
        and len(edges) == 33
        and [len(parity_classes[color]) for color in (0, 1)] == [9, 9]
        and len(set(edges)) == len(edges)
        and all(parity(left) != parity(right) for left, right in edges)
    )
    add_check(
        checks,
        "explicit_open_box_geometry_and_bipartition",
        geometry_holds,
        "explicit coordinates give 18 vertices, 33 simple edges, and a verified 9+9 parity bipartition",
    )

    case_payloads: list[dict[str, object]] = []
    case_signatures: list[Signature] = []
    case_factor_records: list[list[dict[str, object]]] = []
    face_gate = True
    star_gate = True
    symbolic_gate = True
    rational_gate = True

    for eliminated_parity in (0, 1):
        selection = boundary_selection(vertices, eliminated_parity)
        retained = tuple(selection["retained"])
        boundary = tuple(selection["boundary"])
        internal = tuple(selection["internal"])
        face_subsets = selection["face_subsets"]
        selected_face = int(selection["selected_face"])
        signature, star_stats = star_dp_signature(
            vertices,
            edges,
            eliminated_parity,
            boundary,
            internal,
        )
        symbolic, factor_records = target_symbolic_certificate(
            signature, boundary
        )
        audits = {
            str(point): rational_pfaffian_audit(signature, point)
            for point in AUDIT_POINTS
        }

        face_gate &= (
            sorted(len(sites) for sites in face_subsets.values()) == [4, 5]
            and len(boundary) == 5
            and len(internal) == 4
            and len(face_subsets[selected_face]) == 5
            and selected_face == 1 - eliminated_parity
        )
        degree_histogram = Counter(
            int(row["degree"]) for row in star_stats["star_rows"]
        )
        star_gate &= (
            degree_histogram == {3: 4, 4: 4, 5: 1}
            and all(
                int(row["even_local_subsets"]) == 1 << (int(row["degree"]) - 1)
                for row in star_stats["star_rows"]
            )
            and signature[0][0] == 1
            and all(
                signature[mask] == (0,)
                for mask in range(1 << BOUNDARY_LEGS)
                if mask.bit_count() % 2
            )
        )
        generated_factor_data = {
            str(record["class_id"]): expected_factor_projection(record)
            for record in factor_records
        }
        symbolic_gate &= (
            symbolic["ordering_count"] == 120
            and symbolic["principal_pfaffians_checked"] == 1920
            and symbolic["identically_zero_empty_and_pair_residuals"] == 1320
            and symbolic["four_leg_residuals_checked"] == 600
            and symbolic["nonzero_four_leg_residuals"] == 600
            and symbolic["passing_orderings"] == []
            and len(
                {
                    tuple(row["permutation"])
                    for row in symbolic["ordering_records"]
                }
            )
            == 120
            and generated_factor_data == EXPECTED_FACTOR_DATA
            and all(
                record["factorization_over_ZZ"]["reconstructs_residual"]
                and record["coefficient_positivity_certificate"][
                    "positive_scalar"
                ]
                and record["coefficient_positivity_certificate"][
                    "all_remainder_coefficients_positive"
                ]
                for record in factor_records
            )
        )
        rational_gate &= all(
            audit["principal_pfaffians_checked"] == 1920
            and audit["matching_principal_pfaffians"] == 1320
            and audit["failing_principal_pfaffians"] == 600
            and audit["passing_permutations"] == []
            and audit["permutation_failure_count_histogram"] == {"5": 120}
            for audit in audits.values()
        )

        case_payloads.append(
            {
                "case": (
                    "eliminate_even_retain_odd"
                    if eliminated_parity == 0
                    else "eliminate_odd_retain_even"
                ),
                "eliminated_color": {
                    "parity_x_plus_y_plus_z_mod_2": eliminated_parity,
                    "coordinates": [
                        coordinate_record(vertex)
                        for vertex in parity_classes[eliminated_parity]
                    ],
                },
                "retained_color": {
                    "parity_x_plus_y_plus_z_mod_2": 1 - eliminated_parity,
                    "coordinates": [
                        coordinate_record(vertex) for vertex in retained
                    ],
                },
                "face_counts": {
                    f"x={face}": len(sites)
                    for face, sites in face_subsets.items()
                },
                "selected_face": {
                    "axis": "x",
                    "coordinate": selected_face,
                    "reason": "unique 3x3 x-face with five sites of the retained color",
                },
                "boundary_leg_order": [
                    coordinate_record(vertex) for vertex in boundary
                ],
                "boundary_order_rule": "increasing (y,z) on the selected x-face",
                "averaged_retained_coordinates": [
                    coordinate_record(vertex) for vertex in internal
                ],
                "local_stars": star_stats,
                "boundary_signature": signature_record(signature),
                "symbolic_pfaffian_coverage": symbolic,
                "rational_controls": audits,
            }
        )
        case_signatures.append(signature)
        case_factor_records.append(factor_records)

    add_check(
        checks,
        "larger_opposite_color_face_selected",
        face_gate,
        "each color choice has face counts 4 and 5; the unique five-site face is retained in fixed (y,z) order and four other spins are averaged",
    )
    add_check(
        checks,
        "exact_all_degree_star_contraction",
        star_gate and case_signatures[0] == case_signatures[1],
        "nine f_d factors with degree histogram 3^4,4^4,5^1 give identical pure-even five-leg signatures for the reflected color choices",
    )

    control_vertices, control_edges = graph((2, 2, 2))
    control_boundary = tuple(
        sorted(
            (
                vertex
                for vertex in control_vertices
                if parity(vertex) == 1 and vertex[0] == 1
            ),
            key=lambda vertex: (vertex[1], vertex[2]),
        )
    )
    direct_control = direct_spin_signature(
        control_vertices, control_edges, 0, control_boundary
    )
    edge_control = edge_subset_signature(
        control_vertices, control_edges, 0, control_boundary
    )
    add_check(
        checks,
        "small_direct_spin_sum_equals_edge_boundary_enumeration",
        len(control_edges) == 12
        and len(control_boundary) == 2
        and direct_control == edge_control,
        "on open 2x2x2, all 16 retained-spin assignments and all 4096 edge subsets give the same symbolic boundary Walsh signature",
    )

    add_check(
        checks,
        "complete_symbolic_ordering_and_pfaffian_coverage",
        symbolic_gate and case_factor_records[0] == case_factor_records[1],
        "for each color: all 120 leg permutations times all 16 even principal subsets; all 600 four-leg residuals fall into four exact ZZ[v] factor classes",
    )
    add_check(
        checks,
        "all_v_positive_obstruction",
        symbolic_gate,
        "the four factor classes are positive products of v, 1-v, 1+v, 1+v^2, and a polynomial with strictly positive coefficients in v^2",
    )
    add_check(
        checks,
        "exact_one_third_and_one_half_controls",
        rational_gate,
        "at v=1/3 and v=1/2, each color has 0 passing orders and exactly five failing four-leg Pfaffians per order",
    )

    cycle_signature = alternating_cycle_signature(BOUNDARY_LEGS)
    cycle_rows, cycle_four_residuals, cycle_zero, cycle_nonzero = symbolic_coverage(
        cycle_signature
    )
    cycle_passing = {
        tuple(row["permutation"])
        for row in cycle_rows
        if row["all_even_principal_residuals_zero"]
    }
    cycle_audits = {
        str(point): rational_pfaffian_audit(cycle_signature, point)
        for point in AUDIT_POINTS
    }
    expected_dihedral = dihedral_orders(BOUNDARY_LEGS)
    planar_gate = (
        len(cycle_rows) == 120
        and cycle_zero == 1520
        and cycle_nonzero == 400
        and len(cycle_four_residuals) == 600
        and cycle_passing == expected_dihedral
        and all(
            audit["principal_pfaffians_checked"] == 1920
            and audit["matching_principal_pfaffians"] == 1520
            and audit["failing_principal_pfaffians"] == 400
            and {
                tuple(order) for order in audit["passing_permutations"]
            }
            == expected_dihedral
            for audit in cycle_audits.values()
        )
    )
    add_check(
        checks,
        "planar_alternating_cycle_positive_control",
        planar_gate,
        "the C10 alternating five-star signature passes all even principal Pfaffians in exactly its 10 cyclic/reversed boundary orders under the same implementation",
    )

    topology = {
        "shape": list(SHAPE),
        "open_boundary_conditions": True,
        "vertices": len(vertices),
        "edges": len(edges),
        "edge_count_formula": "3abc-(ab+bc+ca)",
        "edge_count_substitution": "3*2*3*3-(2*3+3*3+3*2)=33",
        "e241_genus_bound_numerator": 1,
        "e241_orientable_genus_lower_bound": 1,
        "conclusion": "nonplanar by the positive e241 orientable-genus lower bound; no exact genus is asserted",
        "inherited_source": {
            "path": "experiments/e241_box_genus_bound.py",
            "sha256": source_hashes["experiments/e241_box_genus_bound.py"],
            "binding": "current source hash; no e241 result artifact is imported",
        },
        "terminology_boundary": (
            "no periodic seam signs or surface-sector interpretation is used; "
            "flat graph-H1 twists are not called spin structures"
        ),
    }
    topology_gate = (
        len(edges) == 33
        and 3 * 2 * 3 * 3 - (2 * 3 + 3 * 3 + 3 * 2) == 33
        and 2 * 3 * 3 - 2 * 3 - 3 * 3 - 3 * 2 + 4 == 1
        and topology["e241_orientable_genus_lower_bound"] == 1
    )
    add_check(
        checks,
        "e241_nonplanarity_provenance",
        topology_gate,
        "33 edges and the positive lower-bound numerator 1 are tied to the current e241 source without asserting exact genus",
    )

    scope = {
        "proved": [
            "for each checkerboard color choice on the fixed open 2x3x3 graph, the specified five-leg Walsh-basis boundary contraction fails the sub-Pfaffian equations for every one of 120 boundary-leg orders throughout 0<v<1",
            "all four exact residual factor classes and their strict positivity on 0<v<1",
            "the planar alternating C10 control passes in exactly the ten dihedral boundary orders",
        ],
        "not_proved": [
            "an obstruction after arbitrary local basis changes, active fermionic-SWAP/Koszul phase decorations, arbitrary hidden auxiliaries, or unrestricted tensor-network rewrites",
            "a lower bound on the number of Pfaffians in any global formula",
            "the same boundary obstruction for other box sizes or all lengths L",
            "global Pfaffian solvability or nonsolvability of the three-dimensional Ising model",
            "a three-dimensional Ising solution, free energy, critical point, or critical exponents",
        ],
        "precise_boundary": (
            "this is an all-v theorem only for the fixed five-leg boundary contraction, "
            "Walsh basis, and the complete class of its 120 leg orderings"
        ),
    }
    scope_gate = (
        "arbitrary hidden auxiliaries" in scope["not_proved"][0]
        and "fermionic-SWAP" in scope["not_proved"][0]
        and "number of Pfaffians" in scope["not_proved"][1]
        and "all lengths L" in scope["not_proved"][2]
        and "three-dimensional Ising solution" in scope["not_proved"][4]
        and "fixed five-leg boundary contraction" in scope["precise_boundary"]
    )
    add_check(
        checks,
        "scope_boundary",
        scope_gate,
        "fixed ordinary-Walsh contraction/order-class theorem only; fermionic swap phases, arbitrary bases, hidden auxiliaries, Pfaffian-count bounds, all-L claims, and a 3D solution are excluded",
    )
    add_check(
        checks,
        "benchmark_not_used",
        True,
        "only exact integer polynomials and Fraction controls at v=1/3,1/2 are used",
    )

    certificate = {
        "claim_tag": "[THEOREM][EXACT FINITE COMPUTATION][FIXED BOUNDARY CONTRACTION]",
        "parameter": {
            "name": "v=tanh(K)",
            "domain": "0<v<1",
            "benchmark_or_critical_coupling_used": False,
        },
        "local_factor": {
            "formula": "f_d(s_1,...,s_d)=(product_i(1+v*s_i)+product_i(1-v*s_i))/2",
            "Walsh_expansion": "f_d=sum_(A subset [d], |A| even) v^|A| chi_A",
        },
        "geometry": {
            "shape": list(SHAPE),
            "coordinate_ranges": {"x": [0, 1], "y": [0, 2], "z": [0, 2]},
            "vertices": [coordinate_record(vertex) for vertex in vertices],
            "edges": [
                [coordinate_record(left), coordinate_record(right)]
                for left, right in edges
            ],
            "bipartition": {
                "rule": "x+y+z modulo 2",
                "sizes_even_odd": [9, 9],
                "every_edge_crosses": True,
            },
        },
        "checkerboard_cases": case_payloads,
        "residual_factor_classes": sorted(
            case_factor_records[0], key=lambda record: str(record["class_id"])
        ),
        "small_control": {
            "shape": [2, 2, 2],
            "eliminated_parity": 0,
            "boundary_leg_order": [
                coordinate_record(vertex) for vertex in control_boundary
            ],
            "retained_spin_assignments": 16,
            "edge_subsets": 4096,
            "direct_spin_signature": signature_record(direct_control),
            "edge_boundary_signature": signature_record(edge_control),
            "equal": direct_control == edge_control,
        },
        "planar_alternating_cycle_control": {
            "graph": "C10 with cyclic order b0,u0,b1,u1,b2,u2,b3,u3,b4,u4; each eliminated u_i is adjacent to b_i and b_(i+1 mod 5)",
            "boundary_order": ["b0", "b1", "b2", "b3", "b4"],
            "signature": signature_record(cycle_signature),
            "symbolic_coverage": {
                "boundary_permutations_checked": 120,
                "even_principal_pfaffians_per_permutation": 16,
                "principal_pfaffians_checked": 1920,
                "identically_zero_residuals": cycle_zero,
                "nonzero_residuals": cycle_nonzero,
                "passing_permutations": [
                    list(order) for order in sorted(cycle_passing)
                ],
                "passing_permutation_count": len(cycle_passing),
            },
            "rational_controls": cycle_audits,
        },
        "topology_provenance": topology,
        "finite_enumeration_bounds": {
            "target_vertices": 18,
            "target_edges": 33,
            "eliminated_stars_per_color": 9,
            "retained_spin_parity_states": 512,
            "boundary_leg_permutations_per_signature": 120,
            "even_principal_pfaffians_per_permutation": 16,
            "symbolic_four_leg_residuals_per_signature": 600,
            "small_control_edge_subsets": 4096,
            "small_control_spin_assignments": 16,
        },
        "scope": scope,
    }

    elapsed = time.process_time() - started
    peak_rss = peak_rss_bytes()
    rss_measurement = peak_rss_measurement()
    add_check(
        checks,
        "declared_resource_limits",
        elapsed < CPU_BUDGET_SECONDS and peak_rss < RSS_LIMIT_BYTES,
        (
            f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; "
            f"peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes via {rss_measurement}"
        ),
    )

    failures = [str(check["name"]) for check in checks if not check["passed"]]
    if failures:
        raise AssertionError(f"boundary matchgate checks failed: {failures}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e247_boundary_matchgate_defect.py",
            "verifier": "tests/test_boundary_matchgate_defect.py",
            "proof": "proofs/boundary_matchgate_defect.md",
            "artifact": "results/integrability/boundary_matchgate_defect.json",
            "interpreter": sys.executable,
            "arithmetic": "exact integers, fractions, and polynomials over ZZ[v]",
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "peak_rss_measurement": rss_measurement,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "source_sha256": source_hashes,
            "certificate_sha256": canonical_sha256(certificate),
        },
        "data": certificate,
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for check in payload["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"  [{status}] {check['name']}: {check['detail']}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
