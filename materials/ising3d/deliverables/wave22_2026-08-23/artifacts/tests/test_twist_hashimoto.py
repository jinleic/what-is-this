#!/usr/bin/env python3
"""Independent exact verifier for the flat-Z2 Hashimoto twist audit.

This file does not import the producer.  It uses a different spanning tree and
reconstructs every 24x24 twisted Hashimoto determinant from exact traces and
Newton identities, whereas the producer uses the signed Bass--Ihara vertex
determinant.  Planar phase controls are rebuilt directly as directed-edge
matrices.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Callable, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "kac_ward" / "twist_hashimoto.json"
PRODUCER = ROOT / "experiments" / "e239_twist_hashimoto.py"
FAILURES: list[str] = []

Edge = tuple[int, int]
Arc = tuple[int, int]
Zeta8 = tuple[int, int, int, int]


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def edge(left: int, right: int) -> Edge:
    return (left, right) if left < right else (right, left)


def cube() -> tuple[tuple[int, ...], tuple[Edge, ...]]:
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


def arcs(edges: Sequence[Edge]) -> tuple[Arc, ...]:
    return tuple(arc for left, right in edges for arc in ((left, right), (right, left)))


def graph_adjacency(vertices: Sequence[int], edges: Sequence[Edge]) -> dict[int, list[int]]:
    result = {vertex: [] for vertex in vertices}
    for left, right in edges:
        result[left].append(right)
        result[right].append(left)
    for neighbors in result.values():
        neighbors.sort()
    return result


def even_polynomial(vertices: Sequence[int], edges: Sequence[Edge]) -> tuple[int, ...]:
    position = {vertex: index for index, vertex in enumerate(vertices)}
    coefficients = [0] * (len(edges) + 1)
    for mask in range(1 << len(edges)):
        odd = [0] * len(vertices)
        for index, (left, right) in enumerate(edges):
            if (mask >> index) & 1:
                odd[position[left]] ^= 1
                odd[position[right]] ^= 1
        if not any(odd):
            coefficients[mask.bit_count()] += 1
    return tuple(coefficients)


def multiply_polynomials(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return tuple(result)


def first_difference(left: Sequence[int], right: Sequence[int]) -> int | None:
    return next((index for index, values in enumerate(zip(left, right)) if values[0] != values[1]), None)


def direct_hashimoto_coefficients(edges: Sequence[Edge], signs: dict[Edge, int]) -> tuple[int, ...]:
    directed = arcs(edges)
    successors: list[list[tuple[int, int]]] = [[] for _ in directed]
    for row, current in enumerate(directed):
        for column, successor in enumerate(directed):
            if current[1] == successor[0] and successor[1] != current[0]:
                successors[row].append((column, signs[edge(*successor)]))

    size = len(directed)
    power = [[int(i == j) for j in range(size)] for i in range(size)]
    traces = [0]
    for _ in range(1, size + 1):
        following = [[0] * size for _ in range(size)]
        for i, row in enumerate(power):
            for middle, value in enumerate(row):
                if value == 0:
                    continue
                for column, weight in successors[middle]:
                    following[i][column] += value * weight
        power = following
        traces.append(sum(power[index][index] for index in range(size)))

    coefficients = [1]
    for degree in range(1, size + 1):
        numerator = sum(
            coefficients[degree - order] * traces[order]
            for order in range(1, degree + 1)
        )
        if numerator % degree:
            raise ArithmeticError(f"Newton division failed at degree {degree}")
        coefficients.append(-numerator // degree)
    return tuple(coefficients)


def direct_matrix(
    edges: Sequence[Edge],
    signs: dict[Edge, int],
    phase: Callable[[Arc, Arc], sp.Expr] | None = None,
) -> sp.Matrix:
    directed = arcs(edges)
    matrix = sp.zeros(len(directed))
    for row, current in enumerate(directed):
        for column, successor in enumerate(directed):
            if current[1] != successor[0] or successor[1] == current[0]:
                continue
            value: sp.Expr = sp.Integer(signs[edge(*successor)])
            if phase is not None:
                value *= phase(current, successor)
            matrix[row, column] = value
    return matrix


def integer_charpoly(matrix: sp.Matrix) -> tuple[int, ...]:
    result = []
    for value in matrix.charpoly().all_coeffs():
        exact = sp.simplify(value)
        if exact.is_Integer is not True:
            raise ArithmeticError(f"noninteger coefficient {exact}")
        result.append(int(exact))
    return tuple(result)


def planar_phase(
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
    raise ValueError("unexpected planar turn")


def rebuild_controls() -> dict[str, tuple[int, ...]]:
    path_edges = ((0, 1), (1, 2))
    path_signs = {item: 1 for item in path_edges}
    path_coordinates = {0: (0, 0), 1: (1, 0), 2: (2, 0)}
    path_phase_free = integer_charpoly(direct_matrix(path_edges, path_signs))
    path_kac_ward = integer_charpoly(
        direct_matrix(
            path_edges,
            path_signs,
            lambda current, successor: planar_phase(path_coordinates, current, successor),
        )
    )

    square_edges = ((0, 1), (0, 3), (1, 2), (2, 3))
    square_coordinates = {0: (0, 0), 1: (1, 0), 2: (1, 1), 3: (0, 1)}
    plus = {item: 1 for item in square_edges}
    minus = dict(plus)
    minus[(0, 3)] = -1
    phase_free = tuple(integer_charpoly(direct_matrix(square_edges, signs)) for signs in (plus, minus))
    kac_ward = tuple(
        integer_charpoly(
            direct_matrix(
                square_edges,
                signs,
                lambda current, successor: planar_phase(square_coordinates, current, successor),
            )
        )
        for signs in (plus, minus)
    )
    return {
        "path_phase_free": path_phase_free,
        "path_kac_ward": path_kac_ward,
        "square_phase_free_trivial": phase_free[0],
        "square_phase_free_nontrivial": phase_free[1],
        "square_phase_free_average": tuple((a + b) // 2 for a, b in zip(*phase_free)),
        "square_kac_ward_trivial": kac_ward[0],
        "square_kac_ward_nontrivial": kac_ward[1],
        "square_kac_ward_average": tuple((a + b) // 2 for a, b in zip(*kac_ward)),
    }


def closed_walks(
    vertices: Sequence[int], edges: Sequence[Edge], length: int
) -> list[tuple[tuple[Arc, ...], int]]:
    neighbors = graph_adjacency(vertices, edges)
    edge_index = {item: index for index, item in enumerate(edges)}
    result: list[tuple[tuple[Arc, ...], int]] = []
    for initial in arcs(edges):
        walk = [initial]

        def extend(current: Arc) -> None:
            if len(walk) == length:
                if current[1] == initial[0] and initial[1] != current[0]:
                    parity = 0
                    for step in walk:
                        parity ^= 1 << edge_index[edge(*step)]
                    result.append((tuple(walk), parity))
                return
            for finish in neighbors[current[1]]:
                if finish == current[0]:
                    continue
                successor = (current[1], finish)
                walk.append(successor)
                extend(successor)
                walk.pop()

        extend(initial)
    return result


def zadd(left: Zeta8, right: Zeta8) -> Zeta8:
    return tuple(first + second for first, second in zip(left, right))  # type: ignore[return-value]


def zmul(left: Zeta8, right: Zeta8) -> Zeta8:
    raw = [0] * 7
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            raw[i + j] += first * second
    for degree in range(6, 3, -1):
        raw[degree - 4] -= raw[degree]
    return tuple(raw[:4])  # type: ignore[return-value]


def zpower(exponent: int) -> Zeta8:
    exponent %= 8
    result = [0, 0, 0, 0]
    if exponent < 4:
        result[exponent] = 1
    else:
        result[exponent - 4] = -1
    return tuple(result)  # type: ignore[return-value]


def coordinates(site: int) -> tuple[int, int, int]:
    return (site & 1, (site >> 1) & 1, (site >> 2) & 1)


def phase_exponent(current: Arc, successor: Arc) -> int:
    start, middle, finish = coordinates(current[0]), coordinates(current[1]), coordinates(successor[1])
    first = tuple(middle[index] - start[index] for index in range(3))
    second = tuple(finish[index] - middle[index] for index in range(3))
    cross = (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )
    pairing = sum(cross)
    if pairing:
        return 1 if pairing > 0 else -1
    if first == second:
        return 0
    raise ValueError("coordinate ansatz transition undefined")


def phased_prefix_by_matrix(edges: Sequence[Edge], maximum: int = 8) -> tuple[list[Zeta8], list[Zeta8]]:
    directed = arcs(edges)
    transitions: list[list[tuple[int, Zeta8]]] = [[] for _ in directed]
    for row, current in enumerate(directed):
        for column, successor in enumerate(directed):
            if current[1] == successor[0] and successor[1] != current[0]:
                transitions[row].append((column, zpower(phase_exponent(current, successor))))

    zero: Zeta8 = (0, 0, 0, 0)
    one: Zeta8 = (1, 0, 0, 0)
    size = len(directed)
    power = [[one if i == j else zero for j in range(size)] for i in range(size)]
    traces = [zero]
    for _ in range(maximum):
        following = [[zero for _ in range(size)] for _ in range(size)]
        for i, row in enumerate(power):
            for middle, value in enumerate(row):
                if value == zero:
                    continue
                for column, weight in transitions[middle]:
                    following[i][column] = zadd(following[i][column], zmul(value, weight))
        power = following
        trace = zero
        for index in range(size):
            trace = zadd(trace, power[index][index])
        traces.append(trace)

    coefficients = [one]
    for degree in range(1, maximum + 1):
        numerator = zero
        for order in range(1, degree + 1):
            numerator = zadd(numerator, zmul(coefficients[degree - order], traces[order]))
        if any(component % degree for component in numerator):
            raise ArithmeticError("cyclotomic Newton division failed")
        coefficients.append(tuple(-component // degree for component in numerator))  # type: ignore[arg-type]
    return traces, coefficients


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    check("artifact shape", {"meta", "data", "checks"} <= set(artifact))
    producer_checks = artifact["checks"]
    names = [row["name"] for row in producer_checks]
    check(
        "producer checks unique and true",
        len(names) == len(set(names)) and all(row["passed"] for row in producer_checks),
    )
    meta = artifact["meta"]
    check(
        "source digests",
        meta["source_sha256"]
        == {
            "experiments/e239_twist_hashimoto.py": file_digest(PRODUCER),
            "tests/test_twist_hashimoto.py": file_digest(Path(__file__).resolve()),
        },
    )
    check(
        "resource and benchmark guards",
        meta["single_process"] is True
        and meta["process_cpu_seconds"] <= meta["process_cpu_budget_seconds"]
        and meta["peak_rss_bytes"] <= meta["rss_limit_bytes"]
        and meta["rss_limit_bytes"] == 2 * 1024**3
        and meta["benchmark_Kc_used"] is False,
    )

    rebuilt_controls = rebuild_controls()
    controls = artifact["data"]["two_dimensional_controls"]
    path_target = (1, 0, 0, 0, 0)
    square_target = (1, 0, 0, 0, 2, 0, 0, 0, 1)
    square_average = (1, 0, 0, 0, 0, 0, 0, 0, 1)
    check(
        "path positive control rebuilt",
        rebuilt_controls["path_phase_free"] == path_target
        and rebuilt_controls["path_kac_ward"] == path_target
        and tuple(controls["path_phase_free_determinant"]) == path_target
        and tuple(controls["path_planar_kac_ward_determinant"]) == path_target,
    )
    check(
        "C4 positive and negative controls rebuilt",
        rebuilt_controls["square_phase_free_trivial"]
        == (1, 0, 0, 0, -2, 0, 0, 0, 1)
        and rebuilt_controls["square_phase_free_nontrivial"] == square_target
        and rebuilt_controls["square_kac_ward_trivial"] == square_target
        and rebuilt_controls["square_kac_ward_nontrivial"]
        == (1, 0, 0, 0, -2, 0, 0, 0, 1)
        and rebuilt_controls["square_phase_free_average"] == square_average
        and rebuilt_controls["square_kac_ward_average"] == square_average
        and tuple(controls["square_phase_free"]["uniform_twist_average"]) == square_average,
        "the fixed -1 loop phase works; the uniform twist average deletes v^4",
    )
    triangle_vertices = (0, 1, 2)
    triangle_edges = ((0, 1), (1, 2), (0, 2))
    triangle_trivial = {item: 1 for item in triangle_edges}
    triangle_nontrivial = dict(triangle_trivial)
    triangle_nontrivial[(0, 2)] = -1
    triangle_determinants = (
        direct_hashimoto_coefficients(triangle_edges, triangle_trivial),
        direct_hashimoto_coefficients(triangle_edges, triangle_nontrivial),
    )
    triangle_average = tuple(
        (left + right) // 2 for left, right in zip(*triangle_determinants)
    )
    triangle_p = even_polynomial(triangle_vertices, triangle_edges)
    triangle_target = multiply_polynomials(triangle_p, triangle_p)
    theorem = artifact["data"]["all_graph_girth_obstruction"]
    check(
        "triangle independently witnesses all-graph girth theorem",
        triangle_p == (1, 0, 0, 1)
        and triangle_average == (1, 0, 0, 0, 0, 0, 1)
        and triangle_target == (1, 0, 0, 2, 0, 0, 1)
        and first_difference(triangle_average, triangle_target) == 3
        and theorem["claim_tag"] == "[THEOREM]"
        and "girth(G)" in theorem["statement"],
    )

    vertices, edges = cube()
    alternate_tree = {
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
        (4, 6),
        (5, 7),
        (6, 7),
    }
    alternate_chords = tuple(item for item in edges if item not in alternate_tree)
    producer_tree = {tuple(item) for item in artifact["data"]["open_cube"]["gauge"]["tree_edges_sign_plus_one"]}
    check(
        "independent gauge tree differs",
        len(alternate_tree) == 7
        and len(alternate_chords) == 5
        and alternate_tree != producer_tree,
        f"alternate chords={alternate_chords}",
    )

    alternate_determinants: list[tuple[int, ...]] = []
    representative_bits: set[tuple[int, ...]] = set()
    for bits in itertools.product((0, 1), repeat=5):
        signs = {item: 1 for item in edges}
        for bit, chord in zip(bits, alternate_chords):
            if bit:
                signs[chord] = -1
        representative_bits.add(bits)
        alternate_determinants.append(direct_hashimoto_coefficients(edges, signs))
    check(
        "32 alternate-tree representatives rebuilt",
        len(representative_bits) == 32 and len(alternate_determinants) == 32,
    )

    open_cube = artifact["data"]["open_cube"]
    stored_determinants = [
        tuple(row["determinant_coefficients_ascending"]) for row in open_cube["twist_classes"]
    ]
    check(
        "all 32 direct determinants match as gauge-invariant multiset",
        Counter(alternate_determinants) == Counter(stored_determinants)
        and sorted(Counter(alternate_determinants).values()) == [1, 1, 3, 3, 12, 12],
        f"distinct={len(Counter(alternate_determinants))}",
    )

    p = even_polynomial(vertices, edges)
    target = multiply_polynomials(p, p)
    summed = tuple(sum(row[degree] for row in alternate_determinants) for degree in range(25))
    check("twist sum divisible by 32", all(value % 32 == 0 for value in summed))
    average = tuple(value // 32 for value in summed)
    expected_average = (
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        6,
        0,
        0,
        0,
        16,
        0,
        -48,
        0,
        -87,
        0,
        -48,
        0,
        672,
        0,
        -768,
        0,
        256,
    )
    check(
        "cube P and uniform average rebuilt",
        p == (1, 0, 0, 0, 6, 0, 16, 0, 9, 0, 0, 0, 0)
        and target
        == (1, 0, 0, 0, 12, 0, 32, 0, 54, 0, 192, 0, 364, 0, 288, 0, 81, 0, 0, 0, 0, 0, 0, 0, 0)
        and average == expected_average
        and tuple(open_cube["twist_sum_coefficients_ascending"]) == summed
        and tuple(open_cube["uniform_twist_average_coefficients_ascending"]) == average,
    )
    mismatch = first_difference(average, target)
    check(
        "first mismatch exact",
        mismatch == 4
        and average[4] == 0
        and target[4] == 12
        and open_cube["first_mismatch"]
        == {
            "degree": 4,
            "average_coefficient": 0,
            "target_coefficient": 12,
            "average_minus_target": -12,
        },
    )

    circuit_rows: dict[int, list[tuple[tuple[Arc, ...], int]]] = {
        length: closed_walks(vertices, edges, length) for length in (4, 6, 8)
    }
    observed_census = {
        length: (len(rows), sum(mask == 0 for _, mask in rows))
        for length, rows in circuit_rows.items()
    }
    trivial_eight = [walk for walk, mask in circuit_rows[8] if mask == 0]
    doubles = sum(
        all(walk[index] == walk[index + 4] for index in range(4))
        for walk in trivial_eight
    )
    check(
        "homology-trivial circuit classification rebuilt",
        observed_census == {4: (48, 0), 6: (192, 0), 8: (336, 48)}
        and doubles == 48
        and open_cube["uniform_twist_average_coefficients_ascending"][8] == 6,
        f"census={observed_census}, face_doubles={doubles}",
    )

    traces, phased_coefficients = phased_prefix_by_matrix(edges)
    rational_prefix = tuple(value[0] for value in phased_coefficients)
    ansatz = artifact["data"]["coordinate_phase_ansatz"]
    check(
        "frame-dependent phase ansatz independently rebuilt",
        all(value[1:] == (0, 0, 0) for value in phased_coefficients)
        and rational_prefix == (1, 0, 0, 0, 12, 0, 0, 0, 30)
        and [list(value) for value in traces] == ansatz["trace_zeta8_basis_1_zeta_zeta2_zeta3"]
        and [list(value) for value in phased_coefficients]
        == ansatz["determinant_prefix_zeta8_basis_1_zeta_zeta2_zeta3"]
        and ansatz["first_mismatch_degree"] == 6
        and ansatz["twists_applied"] is False,
    )

    content = meta["content_sha256"]
    comparison_payload = {
        "P": list(p),
        "P_squared": list(target),
        "twist_sum": list(summed),
        "twist_average": list(average),
        "first_mismatch_degree": mismatch,
    }
    check(
        "content digests",
        content["two_dimensional_controls"] == digest(controls)
        and content["twist_class_determinants"] == digest(open_cube["twist_classes"])
        and content["distinct_determinant_inventory"]
        == digest(open_cube["distinct_determinant_inventory"])
        and content["cube_polynomial_comparison"] == digest(comparison_payload)
        and content["homology_circuit_certificate"]
        == digest(artifact["data"]["homology_trivial_circuits"])
        and content["coordinate_phase_ansatz"] == digest(ansatz),
    )
    check(
        "every stored row digest",
        all(
            row["coefficients_sha256"]
            == digest(row["determinant_coefficients_ascending"])
            for row in open_cube["twist_classes"]
        ),
    )

    scope = artifact["data"]["scope"]
    check(
        "finite scope and thermodynamic nonclaim",
        scope["larger_graph_launched"] is False
        and any("thermodynamic" in item for item in scope["not_claimed"])
        and any("universal" in item for item in scope["not_claimed"])
        and "[FINITE NEGATIVE]" in artifact["data"]["claim_tag"],
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} independent checks")
        return 1
    print("PASS: independent flat-Z2 Hashimoto twist verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
