#!/usr/bin/env python3
"""Exact flat-Z2-twist audit of phase-free Hashimoto determinants.

The finite experiment deliberately stops at the open 2x2x2 cube after the
first mismatch is classified.  It also contains exact path/C4 controls that
separate an unphased nonbacktracking matrix from planar Kac--Ward turning
phases, plus one explicitly frame-dependent 3D phase ansatz.

Run from the repository root with

    .venv/bin/python experiments/e239_twist_hashimoto.py

The script is single-process and writes
``results/kac_ward/twist_hashimoto.json``.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import sys
import time
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "kac_ward" / "twist_hashimoto.json"
VERIFIER = ROOT / "tests" / "test_twist_hashimoto.py"
CPU_BUDGET_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3
V = sp.symbols("v")

Edge = tuple[int, int]
Arc = tuple[int, int]
Zeta8 = tuple[int, int, int, int]


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_check(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def canonical_edge(left: int, right: int) -> Edge:
    return (left, right) if left < right else (right, left)


def directed_arcs(edges: Sequence[Edge]) -> tuple[Arc, ...]:
    return tuple(arc for left, right in edges for arc in ((left, right), (right, left)))


def cube_graph() -> tuple[tuple[int, ...], tuple[Edge, ...]]:
    vertices = tuple(range(8))
    edges = tuple(
        sorted(
            (site, site ^ (1 << axis))
            for site in vertices
            for axis in range(3)
            if site < (site ^ (1 << axis))
        )
    )
    return vertices, edges


def cube_coordinates(site: int) -> tuple[int, int, int]:
    return (site & 1, (site >> 1) & 1, (site >> 2) & 1)


def adjacency(vertices: Sequence[int], edges: Sequence[Edge]) -> dict[int, list[int]]:
    rows = {vertex: [] for vertex in vertices}
    for left, right in edges:
        rows[left].append(right)
        rows[right].append(left)
    for row in rows.values():
        row.sort()
    return rows


def bfs_spanning_tree(
    vertices: Sequence[int], edges: Sequence[Edge], root: int
) -> tuple[Edge, ...]:
    rows = adjacency(vertices, edges)
    seen = {root}
    queue = deque([root])
    tree: list[Edge] = []
    while queue:
        vertex = queue.popleft()
        for neighbor in rows[vertex]:
            if neighbor in seen:
                continue
            seen.add(neighbor)
            queue.append(neighbor)
            tree.append(canonical_edge(vertex, neighbor))
    if len(seen) != len(vertices):
        raise ValueError("graph is disconnected")
    return tuple(tree)


def tree_path(tree: Sequence[Edge], start: int, finish: int) -> tuple[Edge, ...]:
    rows: dict[int, list[int]] = {}
    for left, right in tree:
        rows.setdefault(left, []).append(right)
        rows.setdefault(right, []).append(left)
    parent: dict[int, int | None] = {start: None}
    queue = deque([start])
    while queue and finish not in parent:
        vertex = queue.popleft()
        for neighbor in rows[vertex]:
            if neighbor not in parent:
                parent[neighbor] = vertex
                queue.append(neighbor)
    if finish not in parent:
        raise ValueError("tree path not found")
    path: list[Edge] = []
    vertex = finish
    while parent[vertex] is not None:
        previous = parent[vertex]
        assert previous is not None
        path.append(canonical_edge(previous, vertex))
        vertex = previous
    path.reverse()
    return tuple(path)


def even_subgraph_polynomial(
    vertices: Sequence[int], edges: Sequence[Edge]
) -> tuple[int, ...]:
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    coefficients = [0] * (len(edges) + 1)
    for mask in range(1 << len(edges)):
        parity = [0] * len(vertices)
        for index, (left, right) in enumerate(edges):
            if (mask >> index) & 1:
                parity[positions[left]] ^= 1
                parity[positions[right]] ^= 1
        if not any(parity):
            coefficients[mask.bit_count()] += 1
    return tuple(coefficients)


def convolve(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return tuple(result)


def pad(coefficients: Sequence[int], length: int) -> tuple[int, ...]:
    if len(coefficients) > length:
        raise ValueError("cannot shorten polynomial")
    return tuple(coefficients) + (0,) * (length - len(coefficients))


def first_mismatch(left: Sequence[int], right: Sequence[int]) -> int | None:
    if len(left) != len(right):
        raise ValueError("polynomial lengths differ")
    return next((degree for degree, pair in enumerate(zip(left, right)) if pair[0] != pair[1]), None)


def nonbacktracking_matrix(
    vertices: Sequence[int],
    edges: Sequence[Edge],
    edge_signs: dict[Edge, int],
    phase: Callable[[Arc, Arc], sp.Expr] | None = None,
) -> sp.Matrix:
    del vertices
    arcs = directed_arcs(edges)
    matrix = sp.zeros(len(arcs))
    for row, current in enumerate(arcs):
        for column, successor in enumerate(arcs):
            if current[1] != successor[0] or successor[1] == current[0]:
                continue
            weight: sp.Expr = sp.Integer(edge_signs[canonical_edge(*successor)])
            if phase is not None:
                weight *= phase(current, successor)
            matrix[row, column] = weight
    return matrix


def exact_matrix_determinant_coefficients(matrix: sp.Matrix) -> tuple[int, ...]:
    coefficients: list[int] = []
    for raw in matrix.charpoly().all_coeffs():
        coefficient = sp.simplify(raw)
        if coefficient.is_Integer is not True:
            raise ArithmeticError(f"nonintegral determinant coefficient: {coefficient}")
        coefficients.append(int(coefficient))
    return tuple(coefficients)


def planar_half_angle_phase(
    coordinates: dict[int, tuple[int, int]], current: Arc, successor: Arc
) -> sp.Expr:
    first = (
        coordinates[current[1]][0] - coordinates[current[0]][0],
        coordinates[current[1]][1] - coordinates[current[0]][1],
    )
    second = (
        coordinates[successor[1]][0] - coordinates[successor[0]][0],
        coordinates[successor[1]][1] - coordinates[successor[0]][1],
    )
    cross = first[0] * second[1] - first[1] * second[0]
    dot = first[0] * second[0] + first[1] * second[1]
    if cross == 1:
        return (1 + sp.I) / sp.sqrt(2)
    if cross == -1:
        return (1 - sp.I) / sp.sqrt(2)
    if dot == 1:
        return sp.Integer(1)
    raise ValueError(f"unsupported planar transition {current} -> {successor}")


def exact_average(polynomials: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if not polynomials:
        raise ValueError("empty polynomial family")
    numerators = tuple(sum(row[degree] for row in polynomials) for degree in range(len(polynomials[0])))
    denominator = len(polynomials)
    if not all(value % denominator == 0 for value in numerators):
        raise ArithmeticError("twist average is not integral")
    return tuple(value // denominator for value in numerators), numerators


def two_dimensional_controls() -> dict[str, object]:
    path_vertices = (0, 1, 2)
    path_edges = ((0, 1), (1, 2))
    path_coordinates = {0: (0, 0), 1: (1, 0), 2: (2, 0)}
    path_signs = {edge: 1 for edge in path_edges}
    path_p = even_subgraph_polynomial(path_vertices, path_edges)
    path_target = pad(convolve(path_p, path_p), 2 * len(path_edges) + 1)
    path_phase_free = exact_matrix_determinant_coefficients(
        nonbacktracking_matrix(path_vertices, path_edges, path_signs)
    )
    path_kac_ward = exact_matrix_determinant_coefficients(
        nonbacktracking_matrix(
            path_vertices,
            path_edges,
            path_signs,
            lambda current, successor: planar_half_angle_phase(
                path_coordinates, current, successor
            ),
        )
    )

    square_vertices = (0, 1, 2, 3)
    square_edges = ((0, 1), (0, 3), (1, 2), (2, 3))
    square_coordinates = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1)}
    square_p = even_subgraph_polynomial(square_vertices, square_edges)
    square_target = pad(convolve(square_p, square_p), 2 * len(square_edges) + 1)
    positive_signs = {edge: 1 for edge in square_edges}
    negative_holonomy_signs = dict(positive_signs)
    negative_holonomy_signs[(0, 3)] = -1

    phase_free = [
        exact_matrix_determinant_coefficients(
            nonbacktracking_matrix(square_vertices, square_edges, signs)
        )
        for signs in (positive_signs, negative_holonomy_signs)
    ]
    phase_free_average, phase_free_sum = exact_average(phase_free)
    kac_ward = [
        exact_matrix_determinant_coefficients(
            nonbacktracking_matrix(
                square_vertices,
                square_edges,
                signs,
                lambda current, successor: planar_half_angle_phase(
                    square_coordinates, current, successor
                ),
            )
        )
        for signs in (positive_signs, negative_holonomy_signs)
    ]
    kac_ward_average, kac_ward_sum = exact_average(kac_ward)

    return {
        "claim_tag": "[COMPUTATION][FINITE 2D CONTROLS]",
        "path_P": list(path_p),
        "path_target_P_squared": list(path_target),
        "path_phase_free_determinant": list(path_phase_free),
        "path_planar_kac_ward_determinant": list(path_kac_ward),
        "path_interpretation": (
            "The acyclic path has no closed nonbacktracking circuit, so both matrices give 1; "
            "this is a positive implementation control but cannot distinguish phase data."
        ),
        "square_P": list(square_p),
        "square_target_P_squared": list(square_target),
        "square_phase_free": {
            "trivial_holonomy": list(phase_free[0]),
            "nontrivial_holonomy": list(phase_free[1]),
            "twist_sum": list(phase_free_sum),
            "uniform_twist_average": list(phase_free_average),
        },
        "square_planar_kac_ward": {
            "trivial_holonomy": list(kac_ward[0]),
            "nontrivial_holonomy": list(kac_ward[1]),
            "twist_sum": list(kac_ward_sum),
            "uniform_twist_average": list(kac_ward_average),
        },
        "indispensable_local_data": (
            "On C4 the phase-free trivial sector is (1-v^4)^2.  The planar signed half-angle "
            "turns have loop product -1 and give (1+v^4)^2.  A fixed nontrivial Z2 holonomy "
            "can mimic that one loop sign, but uniform averaging gives 1+v^8 and deletes the "
            "degree-four term.  Thus averaging does not replace the local signed-turn/loop phase."
        ),
    }


def bass_hashimoto_coefficients(
    vertices: Sequence[int], edges: Sequence[Edge], edge_signs: dict[Edge, int]
) -> tuple[int, ...]:
    """Signed Bass--Ihara reduction of det(I-vB) to a vertex determinant."""

    size = len(vertices)
    positions = {vertex: index for index, vertex in enumerate(vertices)}
    signed_adjacency = sp.zeros(size)
    degrees = [0] * size
    for edge in edges:
        left, right = edge
        i, j = positions[left], positions[right]
        sign = edge_signs[edge]
        signed_adjacency[i, j] = sign
        signed_adjacency[j, i] = sign
        degrees[i] += 1
        degrees[j] += 1
    excess = len(edges) - len(vertices)
    if excess < 0:
        raise ValueError("Bass polynomial factor would have negative exponent")
    core = sp.eye(size) - V * signed_adjacency + V**2 * sp.diag(
        *(degree - 1 for degree in degrees)
    )
    expression = sp.expand((1 - V**2) ** excess * core.det(method="domain-ge"))
    polynomial = sp.Poly(expression, V, domain=sp.ZZ)
    return tuple(int(polynomial.nth(degree)) for degree in range(2 * len(edges) + 1))


def cube_faces(edges: Sequence[Edge]) -> tuple[dict[str, object], ...]:
    vertices = tuple(range(8))
    axis_names = ("x", "y", "z")
    faces: list[dict[str, object]] = []
    for fixed_axis in range(3):
        for fixed_value in (0, 1):
            face_vertices = {
                site for site in vertices if ((site >> fixed_axis) & 1) == fixed_value
            }
            face_edges = tuple(
                edge for edge in edges if edge[0] in face_vertices and edge[1] in face_vertices
            )
            faces.append(
                {
                    "label": f"{axis_names[fixed_axis]}={fixed_value}",
                    "edges": face_edges,
                }
            )
    return tuple(faces)


def enumerate_closed_nonbacktracking_walks(
    vertices: Sequence[int], edges: Sequence[Edge], length: int
) -> list[tuple[tuple[Arc, ...], int]]:
    if length <= 0:
        return []
    rows = adjacency(vertices, edges)
    edge_positions = {edge: index for index, edge in enumerate(edges)}
    walks: list[tuple[tuple[Arc, ...], int]] = []
    for start in directed_arcs(edges):
        path = [start]

        def extend(current: Arc) -> None:
            if len(path) == length:
                if current[1] == start[0] and start[1] != current[0]:
                    parity_mask = 0
                    for arc in path:
                        parity_mask ^= 1 << edge_positions[canonical_edge(*arc)]
                    walks.append((tuple(path), parity_mask))
                return
            for neighbor in rows[current[1]]:
                if neighbor == current[0]:
                    continue
                successor = (current[1], neighbor)
                path.append(successor)
                extend(successor)
                path.pop()

        extend(start)
    return walks


def circuit_certificate(
    vertices: Sequence[int], edges: Sequence[Edge]
) -> dict[str, object]:
    census: dict[str, object] = {}
    walks_by_length: dict[int, list[tuple[tuple[Arc, ...], int]]] = {}
    for length in (4, 6, 8):
        walks = enumerate_closed_nonbacktracking_walks(vertices, edges, length)
        walks_by_length[length] = walks
        classes = Counter(mask for _, mask in walks)
        census[str(length)] = {
            "rooted_closed_nonbacktracking_walks": len(walks),
            "homology_trivial_rooted_walks": classes[0],
            "distinct_mod2_edge_classes_including_zero_if_present": len(classes),
            "class_multiplicity_histogram": {
                str(multiplicity): count
                for multiplicity, count in sorted(Counter(classes.values()).items())
            },
        }
    trivial_eight = [walk for walk, mask in walks_by_length[8] if mask == 0]
    face_double = [
        walk
        for walk in trivial_eight
        if all(walk[index] == walk[index + 4] for index in range(4))
    ]
    return {
        "claim_tag": "[COMPUTATION][FINITE OPEN CUBE][THROUGH LENGTH 8]",
        "census": census,
        "length_four_classification": (
            "All 48 rooted walks are the 6 square faces times 2 orientations times 4 roots; "
            "all have nonzero mod-2 edge class, so the twist projection retains none."
        ),
        "length_eight_trivial_classification": {
            "rooted_count": len(trivial_eight),
            "face_double_rooted_count": len(face_double),
            "all_are_twice_around_one_face": len(face_double) == len(trivial_eight),
            "factorization_interpretation": (
                "The first nonconstant determinant survivors are the six products of the two "
                "oppositely oriented primitive face factors.  Their total class is zero and "
                "they contribute +6 v^8 to the uniform twist average."
            ),
        },
        "first_mismatch_interpretation": (
            "At v^4 character orthogonality deletes every oriented face circuit.  P(v)^2 instead "
            "has 12 ordered terms (empty,face) and (face,empty), so its coefficient is 12."
        ),
    }


def zeta8_add(left: Zeta8, right: Zeta8) -> Zeta8:
    return tuple(first + second for first, second in zip(left, right))  # type: ignore[return-value]


def zeta8_multiply(left: Zeta8, right: Zeta8) -> Zeta8:
    raw = [0] * 7
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            raw[i + j] += first * second
    for degree in range(6, 3, -1):
        raw[degree - 4] -= raw[degree]
    return tuple(raw[:4])  # type: ignore[return-value]


def zeta8_power(exponent: int) -> Zeta8:
    reduced = exponent % 8
    if reduced < 4:
        values = [0, 0, 0, 0]
        values[reduced] = 1
        return tuple(values)  # type: ignore[return-value]
    values = [0, 0, 0, 0]
    values[reduced - 4] = -1
    return tuple(values)  # type: ignore[return-value]


def coordinate_turn_exponent(current: Arc, successor: Arc) -> int:
    start = cube_coordinates(current[0])
    middle = cube_coordinates(current[1])
    finish = cube_coordinates(successor[1])
    first = tuple(middle[i] - start[i] for i in range(3))
    second = tuple(finish[i] - middle[i] for i in range(3))
    cross = (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )
    diagonal_pairing = sum(cross)
    if diagonal_pairing == 1:
        return 1
    if diagonal_pairing == -1:
        return -1
    if first == second:
        return 0
    raise ValueError(f"undefined coordinate phase for {current} -> {successor}")


def coordinate_phase_ansatz(
    vertices: Sequence[int], edges: Sequence[Edge], target: Sequence[int]
) -> dict[str, object]:
    zero: Zeta8 = (0, 0, 0, 0)
    one: Zeta8 = (1, 0, 0, 0)
    traces: list[Zeta8] = [zero]
    for length in range(1, 9):
        trace = zero
        for walk, _ in enumerate_closed_nonbacktracking_walks(vertices, edges, length):
            exponent = sum(
                coordinate_turn_exponent(walk[index], walk[(index + 1) % length])
                for index in range(length)
            )
            trace = zeta8_add(trace, zeta8_power(exponent))
        traces.append(trace)

    coefficients: list[Zeta8] = [one]
    for degree in range(1, 9):
        numerator = zero
        for index in range(1, degree + 1):
            numerator = zeta8_add(
                numerator, zeta8_multiply(coefficients[degree - index], traces[index])
            )
        if not all(component % degree == 0 for component in numerator):
            raise ArithmeticError("cyclotomic Newton division failed")
        coefficients.append(
            tuple(-component // degree for component in numerator)  # type: ignore[arg-type]
        )
    if not all(coefficient[1:] == (0, 0, 0) for coefficient in coefficients):
        raise ArithmeticError("coordinate ansatz prefix is unexpectedly nonrational")
    rational = tuple(coefficient[0] for coefficient in coefficients)
    mismatch = first_mismatch(rational, target[:9])
    return {
        "claim_tag": "[COMPUTATION][FINITE ANSATZ][OPEN CUBE THROUGH DEGREE 8]",
        "definition": (
            "Fix the oriented coordinate frame and n=(1,1,1).  A straight transition has "
            "phase 1.  An orthogonal transition d->d' has phase zeta_8^sigma, where "
            "sigma=sign((d cross d') dot n)."
        ),
        "coordinate_dependence": (
            "This rule uses the selected frame/diagonal, is not invariant under the full cubic "
            "group, and is tested only as an ansatz; it is not a 3D Kac--Ward formula."
        ),
        "twists_applied": False,
        "trace_zeta8_basis_1_zeta_zeta2_zeta3": [list(value) for value in traces],
        "determinant_prefix_zeta8_basis_1_zeta_zeta2_zeta3": [
            list(value) for value in coefficients
        ],
        "determinant_prefix": list(rational),
        "target_P_squared_prefix": list(target[:9]),
        "first_mismatch_degree": mismatch,
        "coefficient_at_first_mismatch": rational[mismatch] if mismatch is not None else None,
        "target_at_first_mismatch": target[mismatch] if mismatch is not None else None,
    }


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    controls = two_dimensional_controls()
    path_target = tuple(controls["path_target_P_squared"])
    square_target = tuple(controls["square_target_P_squared"])
    square_phase_free = controls["square_phase_free"]
    square_kac_ward = controls["square_planar_kac_ward"]
    add_check(
        checks,
        "acyclic_path_control",
        tuple(controls["path_phase_free_determinant"]) == path_target
        and tuple(controls["path_planar_kac_ward_determinant"]) == path_target,
        "path: phase-free = planar Kac--Ward = P^2 = 1",
    )
    add_check(
        checks,
        "square_signed_turn_positive_control",
        tuple(square_kac_ward["trivial_holonomy"]) == square_target
        and tuple(square_phase_free["nontrivial_holonomy"]) == square_target,
        "C4: planar half-angle phases, or one fixed negative loop holonomy, give (1+v^4)^2",
    )
    expected_square_average = (1, 0, 0, 0, 0, 0, 0, 0, 1)
    add_check(
        checks,
        "square_phase_free_and_twist_average_negative_controls",
        tuple(square_phase_free["trivial_holonomy"])
        == (1, 0, 0, 0, -2, 0, 0, 0, 1)
        and tuple(square_phase_free["uniform_twist_average"]) == expected_square_average
        and tuple(square_kac_ward["uniform_twist_average"]) == expected_square_average
        and expected_square_average != square_target,
        "C4: uniform H1 average is 1+v^8, not P^2; averaging cannot supply signed turning data",
    )

    vertices, edges = cube_graph()
    tree = bfs_spanning_tree(vertices, edges, root=0)
    tree_set = set(tree)
    chords = tuple(edge for edge in edges if edge not in tree_set)
    beta_one = len(edges) - len(vertices) + 1
    add_check(
        checks,
        "open_cube_topology_and_twist_count",
        len(vertices) == 8
        and len(edges) == 12
        and beta_one == 5
        and len(tree) == 7
        and len(chords) == 5,
        f"n={len(vertices)}, m={len(edges)}, b1={beta_one}, classes={1 << beta_one}",
    )

    fundamental_cycles = []
    for chord in chords:
        cycle = tree_path(tree, chord[0], chord[1]) + (chord,)
        fundamental_cycles.append(
            {
                "chord": list(chord),
                "cycle_edges": [list(edge) for edge in cycle],
                "edge_mask": sum(1 << edges.index(edge) for edge in cycle),
            }
        )

    faces = cube_faces(edges)
    twist_rows: list[dict[str, object]] = []
    determinant_family: list[tuple[int, ...]] = []
    for bits in itertools.product((0, 1), repeat=beta_one):
        signs = {edge: 1 for edge in edges}
        for bit, chord in zip(bits, chords):
            if bit:
                signs[chord] = -1
        coefficients = bass_hashimoto_coefficients(vertices, edges, signs)
        determinant_family.append(coefficients)
        face_holonomies = [
            int(sp.prod(signs[edge] for edge in face["edges"])) for face in faces
        ]
        row = {
            "id": "".join(map(str, bits)),
            "holonomy_bits_in_chord_order": list(bits),
            "chord_signs": [signs[chord] for chord in chords],
            "negative_chords": [list(chord) for bit, chord in zip(bits, chords) if bit],
            "face_holonomies_in_face_order": face_holonomies,
            "determinant_coefficients_ascending": list(coefficients),
            "coefficients_sha256": canonical_sha256(list(coefficients)),
        }
        twist_rows.append(row)

    average, sum_coefficients = exact_average(determinant_family)
    cube_p = even_subgraph_polynomial(vertices, edges)
    cube_target = pad(convolve(cube_p, cube_p), 2 * len(edges) + 1)
    mismatch = first_mismatch(average, cube_target)
    inventory_counter = Counter(determinant_family)
    inventory = [
        {
            "multiplicity": multiplicity,
            "determinant_coefficients_ascending": list(coefficients),
            "coefficients_sha256": canonical_sha256(list(coefficients)),
        }
        for coefficients, multiplicity in sorted(inventory_counter.items())
    ]
    add_check(
        checks,
        "complete_spanning_tree_gauge",
        len(twist_rows) == 32
        and len({row["id"] for row in twist_rows}) == 32
        and set(tree).isdisjoint(chords)
        and set(tree) | set(chords) == set(edges)
        and all(
            row["chord_signs"]
            == [1 if bit == 0 else -1 for bit in row["holonomy_bits_in_chord_order"]]
            for row in twist_rows
        )
        and len(fundamental_cycles) == 5
        and all(
            all(
                ((cycle["edge_mask"] >> edges.index(other_chord)) & 1)
                == int(chord == other_chord)
                for other_chord in chords
            )
            for chord, cycle in zip(chords, fundamental_cycles)
        ),
        "all 2^5 chord-sign representatives enumerated; tree signs fixed to +1",
    )
    add_check(
        checks,
        "cube_even_subgraph_polynomial",
        cube_p == (1, 0, 0, 0, 6, 0, 16, 0, 9, 0, 0, 0, 0)
        and sum(cube_p) == 32,
        f"P coefficients={cube_p}",
    )
    add_check(
        checks,
        "all_twisted_determinants_exact",
        len(determinant_family) == 32
        and all(len(row) == 25 and row[0] == 1 for row in determinant_family)
        and sorted(inventory_counter.values()) == [1, 1, 3, 3, 12, 12],
        f"32 integer polynomials; {len(inventory_counter)} distinct; multiplicities={sorted(inventory_counter.values())}",
    )
    add_check(
        checks,
        "uniform_twist_average_integral",
        all(value % 32 == 0 for value in sum_coefficients),
        f"average coefficients={average}",
    )
    add_check(
        checks,
        "first_cube_mismatch_is_degree_four",
        mismatch == 4 and average[4] == 0 and cube_target[4] == 12,
        f"first mismatch v^{mismatch}: average={average[mismatch]}, target={cube_target[mismatch]}",
    )

    all_face_minus = [
        row for row in twist_rows if row["face_holonomies_in_face_order"] == [-1] * 6
    ]
    if len(all_face_minus) != 1:
        raise AssertionError("expected a unique all-face-negative twist")
    face_minus_coefficients = tuple(all_face_minus[0]["determinant_coefficients_ascending"])
    face_minus_mismatch = first_mismatch(face_minus_coefficients, cube_target)
    add_check(
        checks,
        "fixed_all_face_negative_twist_still_fails",
        face_minus_mismatch == 6
        and face_minus_coefficients[4] == 12
        and face_minus_coefficients[6] == -16
        and cube_target[6] == 32,
        f"class={all_face_minus[0]['id']}; first mismatch v^6: -16 versus 32",
    )

    circuits = circuit_certificate(vertices, edges)
    circuit_census = circuits["census"]
    add_check(
        checks,
        "homology_trivial_circuit_classification",
        circuit_census["4"]["rooted_closed_nonbacktracking_walks"] == 48
        and circuit_census["4"]["homology_trivial_rooted_walks"] == 0
        and circuit_census["6"]["rooted_closed_nonbacktracking_walks"] == 192
        and circuit_census["6"]["homology_trivial_rooted_walks"] == 0
        and circuit_census["8"]["rooted_closed_nonbacktracking_walks"] == 336
        and circuit_census["8"]["homology_trivial_rooted_walks"] == 48
        and circuits["length_eight_trivial_classification"][
            "all_are_twice_around_one_face"
        ],
        "lengths 4,6,8: (total,trivial)=(48,0),(192,0),(336,48); all 48 are face doubles",
    )

    coordinate_ansatz = coordinate_phase_ansatz(vertices, edges, cube_target)
    add_check(
        checks,
        "coordinate_phase_ansatz_scoped_and_falsified",
        coordinate_ansatz["first_mismatch_degree"] == 6
        and coordinate_ansatz["coefficient_at_first_mismatch"] == 0
        and coordinate_ansatz["target_at_first_mismatch"] == 32
        and coordinate_ansatz["twists_applied"] is False,
        "frame-dependent untwisted ansatz matches v^4 and first fails at v^6: 0 versus 32",
    )

    comparison_digest_payload = {
        "P": list(cube_p),
        "P_squared": list(cube_target),
        "twist_sum": list(sum_coefficients),
        "twist_average": list(average),
        "first_mismatch_degree": mismatch,
    }
    content_sha256 = {
        "two_dimensional_controls": canonical_sha256(controls),
        "twist_class_determinants": canonical_sha256(twist_rows),
        "distinct_determinant_inventory": canonical_sha256(inventory),
        "cube_polynomial_comparison": canonical_sha256(comparison_digest_payload),
        "homology_circuit_certificate": canonical_sha256(circuits),
        "coordinate_phase_ansatz": canonical_sha256(coordinate_ansatz),
    }

    elapsed = time.process_time() - started
    peak_rss = max_rss_bytes()
    add_check(
        checks,
        "declared_resource_limits",
        elapsed <= CPU_BUDGET_SECONDS and peak_rss <= RSS_LIMIT_BYTES,
        f"process CPU={elapsed:.6f}s/{CPU_BUDGET_SECONDS}s; peak RSS={peak_rss}/{RSS_LIMIT_BYTES} bytes",
    )
    if not all(bool(check["passed"]) for check in checks):
        failures = [str(check["name"]) for check in checks if not check["passed"]]
        raise AssertionError(f"producer checks failed: {failures}")

    source_paths = (Path(__file__).resolve(), VERIFIER)
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e239_twist_hashimoto.py",
            "verifier": "tests/test_twist_hashimoto.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers and Z[zeta_8]; no floating point",
            "single_process": True,
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_Kc_used": False,
            "source_sha256": {
                str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
            },
            "content_sha256": content_sha256,
        },
        "data": {
            "claim_tag": "[COMPUTATION][FINITE NEGATIVE][OPEN 2x2x2 CUBE]",
            "question": (
                "Does the uniform average over all graph-H1 flat Z2 twists of the phase-free "
                "Hashimoto determinant equal the squared even-subgraph polynomial?"
            ),
            "answer": (
                "No on the finite controls and open cube.  On the cube the exact first mismatch "
                "is v^4: the average coefficient is 0 while [v^4]P(v)^2 is 12."
            ),
            "definition": {
                "hashimoto_entry": (
                    "B_s[(u,v),(v,w)]=s_{vw} for w!=u, and 0 otherwise; s_e is a symmetric "
                    "edge sign.  det(I-v B_s) is reported in ascending powers of v."
                ),
                "flat_twists": (
                    "The cube graph is treated as a 1-complex.  Every Z2 edge cochain is flat "
                    "there; vertex switching is gauge, and the 2^b1 classes are H^1(G;Z2)."
                ),
                "character_projection": (
                    "[LEMMA][ALL FINITE CONNECTED GRAPHS] Uniform character orthogonality keeps "
                    "exactly determinant circuit collections whose total mod-2 edge class is zero."
                ),
                "bass_ihara_method": (
                    "Each exact polynomial was computed as (1-v^2)^(m-n) det(I-v A_s+v^2(D-I)); "
                    "the independent verifier instead uses the 24x24 Hashimoto matrix, traces, "
                    "and Newton identities in a different spanning-tree gauge."
                ),
            },
            "all_graph_girth_obstruction": {
                "claim_tag": "[THEOREM]",
                "statement": "For every finite connected simple graph G containing a cycle, the uniform average over all graph-H1 flat Z2 twists of the phase-free Hashimoto determinant differs from P_G(v)^2. Its first mismatch is at degree g=girth(G): the twist average has coefficient 0 and P_G(v)^2 has coefficient 2*c_g>0, where c_g is the number of simple g-cycles.",
                "proof": [
                    "Every shortest closed nonbacktracking walk in a finite simple cyclic graph is an oriented simple cycle of length g; its nonempty mod-2 edge class is nonzero in H1 of the graph 1-complex.",
                    "Character orthogonality therefore kills every degree-g primitive-cycle contribution in the uniform twist average. No determinant collection of two positive-length circuits can occur before degree 2g.",
                    "The smallest nonempty even edge sets are exactly the c_g simple g-cycles, so [v^g]P_G=c_g and [v^g]P_G^2=2*c_g.",
                    "For a forest the Hashimoto matrix is nilpotent and P_G=1, so the uniform identity holds exactly only on the acyclic boundary.",
                ],
                "scope": "phase-free uniform averaging over every H1(G;Z2) character on finite connected simple graphs; fixed twists, nonuniform sums, turning-phase rules, surfaces with 2-cells, and other determinant constructions are outside the theorem",
            },
            "two_dimensional_controls": controls,
            "open_cube": {
                "vertices_binary_xyz": [list(cube_coordinates(site)) for site in vertices],
                "edges": [list(edge) for edge in edges],
                "n": len(vertices),
                "m": len(edges),
                "b1": beta_one,
                "twist_class_count": 1 << beta_one,
                "gauge": {
                    "tree_root": 0,
                    "tree_edges_sign_plus_one": [list(edge) for edge in tree],
                    "chords_in_bit_order": [list(edge) for edge in chords],
                    "bit_convention": "0 means chord sign +1; 1 means chord sign -1",
                    "fundamental_cycle_basis": fundamental_cycles,
                },
                "faces_in_holonomy_order": [
                    {"label": face["label"], "edges": [list(edge) for edge in face["edges"]]}
                    for face in faces
                ],
                "even_subgraph_P_coefficients_ascending": list(cube_p),
                "target_P_squared_coefficients_ascending": list(cube_target),
                "twist_classes": twist_rows,
                "distinct_determinant_inventory": inventory,
                "twist_sum_coefficients_ascending": list(sum_coefficients),
                "twist_average_denominator": 32,
                "uniform_twist_average_coefficients_ascending": list(average),
                "first_mismatch": {
                    "degree": mismatch,
                    "average_coefficient": average[mismatch],
                    "target_coefficient": cube_target[mismatch],
                    "average_minus_target": average[mismatch] - cube_target[mismatch],
                },
                "unique_all_face_negative_twist": {
                    "id": all_face_minus[0]["id"],
                    "determinant_coefficients_ascending": list(face_minus_coefficients),
                    "first_mismatch_degree": face_minus_mismatch,
                    "coefficient_at_first_mismatch": face_minus_coefficients[face_minus_mismatch],
                    "target_at_first_mismatch": cube_target[face_minus_mismatch],
                },
            },
            "homology_trivial_circuits": circuits,
            "coordinate_phase_ansatz": coordinate_ansatz,
            "scope": {
                "proved_all_finite_graphs": [
                    "uniform H1-character averaging projects onto total mod-2 homology class zero",
                    "the phase-free uniform twist average equals P(v)^2 for forests and fails at degree girth(G) for every connected simple cyclic graph",
                ],
                "computed_exactly": [
                    "path and C4 controls",
                    "all 32 gauge classes of the open 2x2x2 cube",
                    "closed nonbacktracking circuit census through length 8 on that cube",
                    "one untwisted frame-dependent coordinate phase ansatz through degree 8",
                ],
                "not_claimed": [
                    "a no-go for fixed twists, nonuniform twist sums, surface spin structures, or other local turning phases",
                    "a universal 3D Kac--Ward phase rule",
                    "a determinant representation of the 3D Ising model",
                    "a tractable thermodynamic limit, free energy, critical point, or critical exponent",
                ],
                "larger_graph_launched": False,
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for check in payload["checks"]:
        print(f"[{'PASS' if check['passed'] else 'FAIL'}] {check['name']}: {check['detail']}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} exact checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
