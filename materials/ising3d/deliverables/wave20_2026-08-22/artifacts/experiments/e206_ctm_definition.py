"""Exact finite corner-boundary kernels for the Ising model.

The graph-theoretic definition is made before any spectralization.  For the open
orthant Q_d(r)={0,...,r}^d, one tensor leg records the spin assignment on each
coordinate face F_i={x_i=0}.  Compatible face assignments are held fixed, all
other spins are summed, and every open-grid bond is counted exactly once.

For d=2 this object is a matrix.  For d=3 it is honestly a rank-three tensor;
turning it into a matrix requires an additional boundary covector.
"""

from __future__ import annotations

from itertools import product
from math import gcd
from typing import Callable, Iterable, Sequence

PLUS = 1
MINUS = 0
MAX_BULK_SPINS = 8

IntPolynomial = list[int]  # ascending powers of q=e^{-2K}
Vertex = tuple[int, ...]


def orthant_vertices(dimension: int, radius: int) -> list[Vertex]:
    if dimension < 2:
        raise ValueError("a corner needs dimension at least two")
    if radius < 1:
        raise ValueError("radius must be positive")
    return [tuple(point) for point in product(range(radius + 1), repeat=dimension)]


def orthant_edges(dimension: int, radius: int) -> list[tuple[Vertex, Vertex]]:
    """Open nearest-neighbour bonds, with no periodic or duplicated bonds."""
    edges: list[tuple[Vertex, Vertex]] = []
    for vertex in orthant_vertices(dimension, radius):
        for axis in range(dimension):
            if vertex[axis] < radius:
                neighbour = list(vertex)
                neighbour[axis] += 1
                edges.append((vertex, tuple(neighbour)))
    return edges


def face_vertices(dimension: int, radius: int, axis: int) -> list[Vertex]:
    """Lexicographic bit order on F_axis={x_axis=0}; bit 1 means spin +1."""
    if not 0 <= axis < dimension:
        raise ValueError("face axis out of range")
    remaining = [coordinate for coordinate in range(dimension) if coordinate != axis]
    vertices: list[Vertex] = []
    for local in product(range(radius + 1), repeat=dimension - 1):
        vertex = [0] * dimension
        for coordinate, value in zip(remaining, local):
            vertex[coordinate] = value
        vertices.append(tuple(vertex))
    return vertices


def decode_face_state(
    dimension: int, radius: int, axis: int, state: int
) -> dict[Vertex, int]:
    vertices = face_vertices(dimension, radius, axis)
    if not 0 <= state < (1 << len(vertices)):
        raise ValueError("face state out of range")
    return {vertex: (state >> bit) & 1 for bit, vertex in enumerate(vertices)}


def merge_face_states(
    dimension: int, radius: int, states: Sequence[int]
) -> dict[Vertex, int] | None:
    """Merge compatible face assignments; return None on an overlap conflict."""
    if len(states) != dimension:
        raise ValueError("one state is required for each coordinate face")
    merged: dict[Vertex, int] = {}
    for axis, state in enumerate(states):
        for vertex, spin in decode_face_state(dimension, radius, axis, state).items():
            previous = merged.get(vertex)
            if previous is not None and previous != spin:
                return None
            merged[vertex] = spin
    return merged


def corner_entry_polynomial(
    dimension: int, radius: int, states: Sequence[int]
) -> IntPolynomial:
    """Return sum q^(number of broken bonds) for one boundary-tensor entry."""
    edges = orthant_edges(dimension, radius)
    fixed = merge_face_states(dimension, radius, states)
    if fixed is None:
        return [0] * (len(edges) + 1)

    bulk = [vertex for vertex in orthant_vertices(dimension, radius) if vertex not in fixed]
    if len(bulk) > MAX_BULK_SPINS:
        raise ValueError(
            f"refusing a {len(bulk)}-spin corner census; exact front limit is "
            f"{MAX_BULK_SPINS}"
        )

    coefficients = [0] * (len(edges) + 1)
    for bulk_state in range(1 << len(bulk)):
        spins = dict(fixed)
        for bit, vertex in enumerate(bulk):
            spins[vertex] = (bulk_state >> bit) & 1
        broken = sum(spins[left] != spins[right] for left, right in edges)
        coefficients[broken] += 1
    return coefficients


def fixed_corner_face_state(radius: int, arm_state: int, spin: int = PLUS) -> int:
    """Pack a d=2 face: origin at bit 0, arm sites at bits 1,...,radius."""
    if spin not in (MINUS, PLUS):
        raise ValueError("spin must be encoded by 0 or 1")
    if not 0 <= arm_state < (1 << radius):
        raise ValueError("arm state out of range")
    state = spin
    for offset in range(radius):
        state |= ((arm_state >> offset) & 1) << (offset + 1)
    return state


def corner_matrix_2d_polynomial(
    radius: int = 2, fixed_corner_spin: int = PLUS
) -> list[list[IntPolynomial]]:
    """Fixed-corner block of the canonical d=2 boundary kernel."""
    dimension = 1 << radius
    matrix: list[list[IntPolynomial]] = []
    for row in range(dimension):
        row_state = fixed_corner_face_state(radius, row, fixed_corner_spin)
        entries: list[IntPolynomial] = []
        for column in range(dimension):
            column_state = fixed_corner_face_state(radius, column, fixed_corner_spin)
            entries.append(corner_entry_polynomial(2, radius, (row_state, column_state)))
        matrix.append(entries)
    return matrix


def evaluate_polynomial_with_common_scale(
    coefficients: Sequence[int], numerator: int, denominator: int, degree_bound: int
) -> int:
    """Return denominator^degree_bound * P(numerator/denominator) over Z."""
    if denominator <= 0 or gcd(abs(numerator), denominator) != 1:
        raise ValueError("the rational parameter must be reduced with positive denominator")
    if len(coefficients) > degree_bound + 1:
        raise ValueError("polynomial exceeds the declared common degree bound")
    return sum(
        int(coefficient)
        * numerator**degree
        * denominator ** (degree_bound - degree)
        for degree, coefficient in enumerate(coefficients)
    )


def evaluate_polynomial_matrix_with_common_scale(
    matrix: Sequence[Sequence[Sequence[int]]],
    numerator: int,
    denominator: int,
    degree_bound: int,
) -> list[list[int]]:
    return [
        [
            evaluate_polynomial_with_common_scale(
                entry, numerator, denominator, degree_bound
            )
            for entry in row
        ]
        for row in matrix
    ]


def corner_tensor_3d_r1_scaled(
    numerator: int = 1, denominator: int = 2
) -> tuple[int, list[list[list[int]]]]:
    """The honest 16 x 16 x 16 radius-one octant tensor at rational q."""
    edge_count = len(orthant_edges(3, 1))
    scale = denominator**edge_count
    tensor = [[[0] * 16 for _ in range(16)] for _ in range(16)]
    for first in range(16):
        for second in range(16):
            for third in range(16):
                polynomial = corner_entry_polynomial(3, 1, (first, second, third))
                tensor[first][second][third] = evaluate_polynomial_with_common_scale(
                    polynomial, numerator, denominator, edge_count
                )
    return scale, tensor


def contract_third_leg(
    tensor: Sequence[Sequence[Sequence[int]]], covector: Sequence[int]
) -> list[list[int]]:
    """Contract a declared third-face boundary covector; this is extra data."""
    dimension = len(tensor)
    if len(covector) != dimension:
        raise ValueError("covector dimension does not match the tensor leg")
    return [
        [
            sum(
                int(tensor[row][column][state]) * int(covector[state])
                for state in range(dimension)
            )
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]


def shared_edge_sector_r1(
    matrix: Sequence[Sequence[int]], fixed_spin: int = PLUS
) -> list[list[int]]:
    """Fix F_0 intersect F_1 (the two-site z edge) before comparing spectra."""
    if len(matrix) != 16 or any(len(row) != 16 for row in matrix):
        raise ValueError("the radius-one face matrix must be 16 by 16")
    if fixed_spin not in (MINUS, PLUS):
        raise ValueError("spin must be encoded by 0 or 1")
    shared_vertices = set(face_vertices(3, 1, 0)) & set(face_vertices(3, 1, 1))
    positions = [
        bit
        for bit, vertex in enumerate(face_vertices(3, 1, 0))
        if vertex in shared_vertices
    ]
    indices = [
        state
        for state in range(16)
        if all(((state >> bit) & 1) == fixed_spin for bit in positions)
    ]
    return [[int(matrix[row][column]) for column in indices] for row in indices]


def radius_one_face_covectors() -> dict[str, list[int]]:
    """Two positive symmetry-invariant choices and their aligned-face direction."""
    free = [1] * 16
    aligned = [int(state in (0, 15)) for state in range(16)]
    biased = [free[state] + aligned[state] for state in range(16)]
    return {"free": free, "aligned": aligned, "biased": biased}


def permute_vertex(vertex: Vertex, permutation: Sequence[int]) -> Vertex:
    """Coordinate action new_vertex[i]=vertex[permutation[i]]."""
    if sorted(permutation) != list(range(len(vertex))):
        raise ValueError("not a coordinate permutation")
    return tuple(vertex[permutation[index]] for index in range(len(vertex)))


def permute_face_state(
    dimension: int,
    radius: int,
    old_axis: int,
    state: int,
    permutation: Sequence[int],
) -> tuple[int, int]:
    """Transport a face state under a coordinate permutation."""
    new_axis = list(permutation).index(old_axis)
    old_vertices = face_vertices(dimension, radius, old_axis)
    new_vertices = face_vertices(dimension, radius, new_axis)
    new_position = {vertex: bit for bit, vertex in enumerate(new_vertices)}
    transported = 0
    for bit, vertex in enumerate(old_vertices):
        spin = (state >> bit) & 1
        image = permute_vertex(vertex, permutation)
        transported |= spin << new_position[image]
    return new_axis, transported


def transform_square_face_state(
    state: int, transform: Callable[[int, int], tuple[int, int]]
) -> int:
    """Transport a four-spin radius-one face state under a square symmetry."""
    if not 0 <= state < 16:
        raise ValueError("face state out of range")
    positions = {(u, v): bit for bit, (u, v) in enumerate(product(range(2), repeat=2))}
    out = 0
    for coordinate, bit in positions.items():
        image = transform(*coordinate)
        out |= ((state >> bit) & 1) << positions[image]
    return out


def square_face_symmetry_maps() -> list[list[int]]:
    """The eight D4 permutations of the four radius-one face states."""
    maps: list[list[int]] = []
    for swap in (False, True):
        for flip_first in (False, True):
            for flip_second in (False, True):
                def transform(
                    first: int,
                    second: int,
                    swap: bool = swap,
                    flip_first: bool = flip_first,
                    flip_second: bool = flip_second,
                ) -> tuple[int, int]:
                    if swap:
                        first, second = second, first
                    if flip_first:
                        first = 1 - first
                    if flip_second:
                        second = 1 - second
                    return first, second

                maps.append([transform_square_face_state(state, transform) for state in range(16)])
    unique = {tuple(mapping) for mapping in maps}
    if len(unique) != 8:
        raise AssertionError("square symmetry construction did not produce D4")
    return [list(mapping) for mapping in sorted(unique)]


def tensor_nonzero_count(tensor: Iterable[Iterable[Iterable[int]]]) -> int:
    return sum(value != 0 for plane in tensor for row in plane for value in row)
