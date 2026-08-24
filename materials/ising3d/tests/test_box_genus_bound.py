#!/usr/bin/env python3
"""Standalone exact verifier for the open-box genus-bound artifact.

This verifier intentionally does not import ``e241_box_genus_bound``.  It
rebuilds every graph from six-neighbor adjacency, enumerates xy/xz/yz unit
plaquettes in separate loops, validates every selected witness, reconstructs
the full 512-row table, and checks the theorem arithmetic and pinned rows.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import time
from collections import deque
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "topology" / "box_genus_bound.json"
PRODUCER = ROOT / "experiments" / "e241_box_genus_bound.py"
FAILURES: list[str] = []
CPU_BUDGET_SECONDS = 60.0
RSS_LIMIT_BYTES = 2 * 1024**3

Coordinate = tuple[int, int, int]
Edge = tuple[Coordinate, Coordinate]
Square = tuple[int, int, Coordinate]


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def edge(left: Coordinate, right: Coordinate) -> Edge:
    if left == right:
        raise AssertionError("loop in independently rebuilt box")
    return (left, right) if left < right else (right, left)


def rebuild_graph(a: int, b: int, c: int) -> tuple[tuple[Coordinate, ...], tuple[Edge, ...]]:
    sides = (a, b, c)
    vertices = tuple(
        (x, y, z) for x in range(a) for y in range(b) for z in range(c)
    )
    vertex_set = set(vertices)
    directions = (
        (-1, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
    )
    edges: set[Edge] = set()
    for vertex in vertices:
        for delta in directions:
            neighbor = tuple(vertex[index] + delta[index] for index in range(3))
            if neighbor in vertex_set:
                edges.add(edge(vertex, neighbor))
    expected_vertices = a * b * c
    if len(vertices) != expected_vertices or any(side < 1 for side in sides):
        raise AssertionError("invalid independent box dimensions")
    return vertices, tuple(sorted(edges))


def connected(vertices: Sequence[Coordinate], edges: Sequence[Edge]) -> bool:
    neighbors = {vertex: [] for vertex in vertices}
    for left, right in edges:
        neighbors[left].append(right)
        neighbors[right].append(left)
    reached = {vertices[-1]}
    queue = deque((vertices[-1],))
    while queue:
        current = queue.popleft()
        for neighbor in neighbors[current]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    return len(reached) == len(vertices)


def independent_squares(a: int, b: int, c: int) -> tuple[Square, ...]:
    squares: list[Square] = []
    for x in range(a - 1):
        for y in range(b - 1):
            for z in range(c):
                squares.append((0, 1, (x, y, z)))
    for x in range(a - 1):
        for y in range(b):
            for z in range(c - 1):
                squares.append((0, 2, (x, y, z)))
    for x in range(a):
        for y in range(b - 1):
            for z in range(c - 1):
                squares.append((1, 2, (x, y, z)))
    return tuple(squares)


def square_vertices(square: Square) -> tuple[Coordinate, Coordinate, Coordinate, Coordinate]:
    first_axis, second_axis, base = square
    first = list(base)
    first[first_axis] += 1
    opposite = list(first)
    opposite[second_axis] += 1
    second = list(base)
    second[second_axis] += 1
    return base, tuple(first), tuple(opposite), tuple(second)


def square_edges(square: Square) -> tuple[Edge, Edge, Edge, Edge]:
    vertices = square_vertices(square)
    return tuple(
        edge(vertices[index], vertices[(index + 1) % 4]) for index in range(4)
    )  # type: ignore[return-value]


def independent_witnesses(
    a: int, b: int, c: int, edges: Sequence[Edge]
) -> tuple[tuple[Square, ...], dict[Edge, tuple[Square, ...]], list[dict[str, object]]]:
    sides = (a, b, c)
    graph_edges = set(edges)
    squares = independent_squares(*sides)
    incidence: dict[Edge, list[Square]] = {current: [] for current in edges}
    for square in squares:
        vertices = square_vertices(square)
        current_edges = square_edges(square)
        valid = (
            len(set(vertices)) == 4
            and all(
                all(0 <= vertex[axis] < sides[axis] for axis in range(3))
                for vertex in vertices
            )
            and all(current in graph_edges for current in current_edges)
            and all(
                sum(abs(left[axis] - right[axis]) for axis in range(3)) == 1
                for left, right in current_edges
            )
        )
        if not valid:
            raise AssertionError(f"bad independently generated unit plaquette {square}")
        for current in current_edges:
            incidence[current].append(square)
    frozen = {current: tuple(sorted(rows)) for current, rows in incidence.items()}
    selected = []
    for current in sorted(edges):
        if not frozen[current]:
            continue
        witness = frozen[current][0]
        selected.append(
            {
                "edge": [list(current[0]), list(current[1])],
                "plaquette": {
                    "axes": [witness[0], witness[1]],
                    "base": list(witness[2]),
                },
            }
        )
    return squares, frozen, selected


def ceiling_quarter(value: int) -> int:
    quotient, remainder = divmod(value, 4)
    return quotient + int(remainder != 0)


def numerator(a: int, b: int, c: int) -> int:
    return a * b * c - a * b - b * c - c * a + 4


def planar_label(a: int, b: int, c: int) -> tuple[bool, str]:
    ordered = sorted((a, b, c))
    if ordered[0] == 1:
        return True, "planar_side_one_slab"
    if ordered[:2] == [2, 2]:
        return True, "planar_2x2xL"
    return False, "nonplanar_by_positive_genus_lower_bound"


def rebuild_row(a: int, b: int, c: int) -> dict[str, object]:
    vertices, edges = rebuild_graph(a, b, c)
    squares, incidence, selected = independent_witnesses(a, b, c, edges)
    even = sum((x + y + z) % 2 == 0 for x, y, z in vertices)
    bipartite = all((sum(left) + sum(right)) % 2 == 1 for left, right in edges)
    simple = len(edges) == len(set(edges)) and all(left != right for left, right in edges)
    axis_counts = {
        "xy": sum(square[:2] == (0, 1) for square in squares),
        "xz": sum(square[:2] == (0, 2) for square in squares),
        "yz": sum(square[:2] == (1, 2) for square in squares),
    }
    witnessed = sum(bool(rows) for rows in incidence.values())
    applies = min(a, b, c) >= 2
    lower: dict[str, object] | None = None
    if applies:
        value = numerator(a, b, c)
        quotient, remainder = divmod(value, 4)
        ceiling = ceiling_quarter(value)
        lower = {
            "numerator": value,
            "floor_quotient": quotient,
            "remainder": remainder,
            "ceiling": ceiling,
            "genus_lower_bound": max(0, ceiling),
            "positive": max(0, ceiling) > 0,
        }
    classified, label = planar_label(a, b, c)
    return {
        "dims": [a, b, c],
        "vertices": len(vertices),
        "edges": len(edges),
        "simple": simple,
        "connected": connected(vertices, edges),
        "bipartite": bipartite,
        "bipartition_sizes_even_odd": [even, len(vertices) - even],
        "unit_plaquette_axis_counts": axis_counts,
        "unit_plaquettes": len(squares),
        "plaquette_edge_incidences": sum(len(rows) for rows in incidence.values()),
        "edges_with_unit_plaquette_witness": witnessed,
        "all_edges_have_unit_plaquette": witnessed == len(edges),
        "selected_witnesses_sha256": digest(selected),
        "genus_bound_theorem_applies": applies,
        "lower_bound": lower,
        "classified_planar": classified,
        "classification": label,
    }


def rebuild_explicit_control(a: int, b: int, c: int) -> dict[str, object]:
    _, edges = rebuild_graph(a, b, c)
    squares, incidence, selected = independent_witnesses(a, b, c, edges)
    return {
        "dims": [a, b, c],
        "edges": len(edges),
        "unit_plaquettes": [
            {"axes": [first, second], "base": list(base)}
            for first, second, base in squares
        ],
        "selected_edge_witnesses": selected,
        "edge_incidence_counts": [
            {
                "edge": [list(current[0]), list(current[1])],
                "incident_unit_plaquettes": len(incidence[current]),
            }
            for current in edges
        ],
        "selected_witnesses_sha256": digest(selected),
    }


def pin_projection(row: dict[str, object]) -> dict[str, object]:
    lower = row["lower_bound"]
    return {
        "dims": row["dims"],
        "vertices": row["vertices"],
        "edges": row["edges"],
        "bipartition_sizes_even_odd": row["bipartition_sizes_even_odd"],
        "unit_plaquette_axis_counts": row["unit_plaquette_axis_counts"],
        "unit_plaquettes": row["unit_plaquettes"],
        "edges_with_unit_plaquette_witness": row["edges_with_unit_plaquette_witness"],
        "lower_bound_numerator": None if lower is None else lower["numerator"],
        "genus_lower_bound": None if lower is None else lower["genus_lower_bound"],
        "classification": row["classification"],
    }


PINNED_VALUES = {
    (1, 1, 1): (1, 0, [1, 0], {"xy": 0, "xz": 0, "yz": 0}, 0, 0, None, None, "planar_side_one_slab"),
    (1, 1, 8): (8, 7, [4, 4], {"xy": 0, "xz": 0, "yz": 0}, 0, 0, None, None, "planar_side_one_slab"),
    (1, 2, 3): (6, 7, [3, 3], {"xy": 0, "xz": 0, "yz": 2}, 2, 7, None, None, "planar_side_one_slab"),
    (2, 1, 2): (4, 4, [2, 2], {"xy": 0, "xz": 1, "yz": 0}, 1, 4, None, None, "planar_side_one_slab"),
    (2, 2, 2): (8, 12, [4, 4], {"xy": 2, "xz": 2, "yz": 2}, 6, 12, 0, 0, "planar_2x2xL"),
    (2, 2, 8): (32, 60, [16, 16], {"xy": 8, "xz": 14, "yz": 14}, 36, 60, 0, 0, "planar_2x2xL"),
    (2, 3, 3): (18, 33, [9, 9], {"xy": 6, "xz": 6, "yz": 8}, 20, 33, 1, 1, "nonplanar_by_positive_genus_lower_bound"),
    (2, 3, 8): (48, 98, [24, 24], {"xy": 16, "xz": 21, "yz": 28}, 65, 98, 6, 2, "nonplanar_by_positive_genus_lower_bound"),
    (3, 3, 3): (27, 54, [14, 13], {"xy": 12, "xz": 12, "yz": 12}, 36, 54, 4, 1, "nonplanar_by_positive_genus_lower_bound"),
    (3, 3, 8): (72, 159, [36, 36], {"xy": 32, "xz": 42, "yz": 42}, 116, 159, 19, 5, "nonplanar_by_positive_genus_lower_bound"),
    (4, 5, 6): (120, 286, [60, 60], {"xy": 72, "xz": 75, "yz": 80}, 227, 286, 50, 13, "nonplanar_by_positive_genus_lower_bound"),
    (8, 8, 8): (512, 1344, [256, 256], {"xy": 392, "xz": 392, "yz": 392}, 1176, 1344, 324, 81, "nonplanar_by_positive_genus_lower_bound"),
}


def pinned_expected() -> list[dict[str, object]]:
    rows = []
    for dims, values in PINNED_VALUES.items():
        (
            vertices,
            edges,
            parity,
            axis_counts,
            plaquettes,
            witnessed,
            lower_numerator,
            lower_bound,
            classification,
        ) = values
        rows.append(
            {
                "dims": list(dims),
                "vertices": vertices,
                "edges": edges,
                "bipartition_sizes_even_odd": parity,
                "unit_plaquette_axis_counts": axis_counts,
                "unit_plaquettes": plaquettes,
                "edges_with_unit_plaquette_witness": witnessed,
                "lower_bound_numerator": lower_numerator,
                "genus_lower_bound": lower_bound,
                "classification": classification,
            }
        )
    return rows


def rebuild_cubes() -> list[dict[str, int]]:
    rows = []
    for side in range(2, 13):
        direct = numerator(side, side, side)
        rows.append(
            {
                "side": side,
                "direct_numerator": direct,
                "factored_numerator": (side - 2) ** 2 * (side + 1),
                "genus_lower_bound": max(0, ceiling_quarter(direct)),
            }
        )
    return rows


def rebuild_elongated() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for a, b in ((2, 2), (2, 3), (2, 4), (3, 3), (3, 4), (4, 5)):
        slope = a * b - a - b
        length_rows = []
        for length in range(2, 17):
            direct = numerator(a, b, length)
            length_rows.append(
                {
                    "length": length,
                    "direct_numerator": direct,
                    "affine_numerator": slope * length - a * b + 4,
                    "genus_lower_bound": max(0, ceiling_quarter(direct)),
                }
            )
        rows.append(
            {
                "cross_section": [a, b],
                "numerator_slope": slope,
                "bound_grows_unbounded_with_length": slope > 0,
                "length_rows": length_rows,
            }
        )
    return rows


def validate_explicit_witnesses(control: dict[str, object]) -> bool:
    a, b, c = control["dims"]
    _, graph_edges = rebuild_graph(a, b, c)
    graph_edge_set = set(graph_edges)
    valid_squares = set(independent_squares(a, b, c))
    records = control["selected_edge_witnesses"]
    seen: set[Edge] = set()
    for record in records:
        endpoints = record["edge"]
        current = edge(tuple(endpoints[0]), tuple(endpoints[1]))
        descriptor = record["plaquette"]
        square = (
            descriptor["axes"][0],
            descriptor["axes"][1],
            tuple(descriptor["base"]),
        )
        if (
            current not in graph_edge_set
            or current in seen
            or square not in valid_squares
            or current not in square_edges(square)
        ):
            return False
        seen.add(current)
    return seen == graph_edge_set and digest(records) == control["selected_witnesses_sha256"]


def main() -> int:
    started = time.process_time()
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    meta = artifact["meta"]
    data = artifact["data"]

    check(
        "artifact_identity",
        meta["producer"] == "experiments/e241_box_genus_bound.py"
        and meta["verifier"] == "tests/test_box_genus_bound.py"
        and meta["arithmetic"] == "exact integers only; no floating point"
        and meta["single_process"] is True
        and meta["benchmark_used"] is False,
    )
    check(
        "source_digests",
        meta["source_sha256"]
        == {
            "experiments/e241_box_genus_bound.py": file_digest(PRODUCER),
            "tests/test_box_genus_bound.py": file_digest(Path(__file__).resolve()),
        },
        "producer and standalone verifier match the generated artifact",
    )
    check(
        "producer_checks",
        len(artifact["checks"]) == 12
        and all(row["passed"] is True for row in artifact["checks"]),
        f"{sum(row['passed'] is True for row in artifact['checks'])}/{len(artifact['checks'])} stored checks pass",
    )

    finite = data["finite_table"]
    rebuilt_rows = [
        rebuild_row(a, b, c)
        for a in range(1, 9)
        for b in range(1, 9)
        for c in range(1, 9)
    ]
    check(
        "independent_full_table",
        finite["side_range"] == [1, 8]
        and finite["ordered_row_count"] == 512
        and finite["rows"] == rebuilt_rows,
        "six-neighbor graphs and separate xy/xz/yz plaquette loops reproduce all 512 rows",
    )

    counts_and_parity = all(
        row["vertices"] == a * b * c
        and row["edges"] == (a - 1) * b * c + a * (b - 1) * c + a * b * (c - 1)
        and row["simple"]
        and row["connected"]
        and row["bipartite"]
        and sum(row["bipartition_sizes_even_odd"]) == row["vertices"]
        for row in rebuilt_rows
        for a, b, c in (tuple(row["dims"]),)
    )
    check(
        "independent_counts_connectivity_and_bipartition",
        counts_and_parity,
        "all ordered boxes satisfy direct direction-count and parity checks",
    )

    all_theorem_edges_witnessed = all(
        row["edges_with_unit_plaquette_witness"] == row["edges"]
        and row["all_edges_have_unit_plaquette"]
        and row["plaquette_edge_incidences"] == 4 * row["unit_plaquettes"]
        for row in rebuilt_rows
        if min(row["dims"]) >= 2
    )
    check(
        "independent_unit_plaquette_witnesses",
        all_theorem_edges_witnessed,
        "every theorem-domain edge has a validated selected 4-cycle witness",
    )

    ceiling_controls = [
        {
            "numerator": value,
            "ceiling_over_four": ceiling_quarter(value),
            "max_zero": max(0, ceiling_quarter(value)),
        }
        for value in range(-12, 21)
    ]
    theorem_ceiling_rows = all(
        row["lower_bound"]["ceiling"] == ceiling_quarter(numerator(*row["dims"]))
        and row["lower_bound"]["genus_lower_bound"]
        == max(0, ceiling_quarter(numerator(*row["dims"])))
        for row in rebuilt_rows
        if row["genus_bound_theorem_applies"]
    )
    check(
        "independent_ceiling_arithmetic",
        data["synthetic_ceiling_controls"] == ceiling_controls and theorem_ceiling_rows,
        "ceil(n/4) rebuilt exactly, including negative controls and every theorem row",
    )

    exact_classification = all(
        row["classified_planar"]
        == (min(row["dims"]) == 1 or sorted(row["dims"])[:2] == [2, 2])
        and (
            min(row["dims"]) == 1
            or row["lower_bound"]["positive"] == (not row["classified_planar"])
        )
        for row in rebuilt_rows
    )
    check(
        "independent_planarity_separation",
        exact_classification,
        "every nonplanar-class row has a positive bound and no planar-class row does",
    )

    expected_pins = pinned_expected()
    rebuilt_by_dims = {tuple(row["dims"]): row for row in rebuilt_rows}
    rebuilt_pins = [pin_projection(rebuilt_by_dims[tuple(row["dims"])]) for row in expected_pins]
    check(
        "hard_coded_pinned_rows",
        data["pinned_rows"] == expected_pins == rebuilt_pins,
        "12 pinned rows match independent literal values",
    )

    expected_cubes = rebuild_cubes()
    check(
        "cube_corollary_rows",
        data["cube_controls"] == expected_cubes
        and [row["genus_lower_bound"] for row in expected_cubes[:7]]
        == [0, 1, 5, 14, 28, 50, 81],
        "(L-2)^2(L+1) identity and L=2..8 bounds reproduced",
    )

    expected_elongated = rebuild_elongated()
    check(
        "elongated_corollary_rows",
        data["elongated_controls"] == expected_elongated
        and all(
            length_row["direct_numerator"] == length_row["affine_numerator"]
            for section in expected_elongated
            for length_row in section["length_rows"]
        ),
        "six fixed cross-sections through length 16 reproduced",
    )

    expected_witness_controls = [
        rebuild_explicit_control(2, 2, 2),
        rebuild_explicit_control(2, 3, 3),
    ]
    stored_witness_controls = data["explicit_unit_plaquette_witness_controls"]
    check(
        "explicit_witness_records",
        stored_witness_controls == expected_witness_controls
        and all(validate_explicit_witnesses(control) for control in stored_witness_controls),
        "all stored cube and 2x3x3 witness records are valid and exhaustive",
    )

    theorem = data["theorem"]
    cellularity_text = " ".join(theorem["cellularity_lemma"]["proof"])
    theorem_statement_text = " ".join(
        (
            theorem["domain"],
            theorem["face_and_euler_argument"]["bound"],
            theorem["planarity_classification"]["statement"],
        )
    )
    relevance = theorem["scoped_relevance_note"]
    check(
        "theorem_and_self_contained_cellularity",
        theorem["counts"]["vertices"] == "V=abc"
        and theorem["counts"]["edges"] == "E=3abc-(ab+bc+ca)"
        and "regular neighborhood" in cellularity_text
        and "Capping each boundary" in cellularity_text
        and "g=h+sum_i k_i+q-r" in cellularity_text
        and "without an external theorem" in cellularity_text
        and "Pfaffian" not in theorem_statement_text,
        "regular-neighborhood-and-cap proof and exact theorem formulas are present",
    )
    check(
        "scoped_surface_relevance",
        "4^g terms" in relevance
        and "not a lower bound" in relevance
        and "arbitrary formula" in relevance
        and "not used in this theorem" in relevance,
        "4^g is stated only as the classical construction's term count",
    )

    exclusions = data["scope"]["explicitly_excluded"]
    check(
        "explicit_exclusions",
        len(exclusions) == 4
        and any("periodic length-two multigraphs" in item for item in exclusions)
        and any("exact orientable genera beyond" in item for item in exclusions)
        and any("arbitrary formulas" in item for item in exclusions)
        and any("thermodynamic solution" in item for item in exclusions),
    )

    content_digests = {
        name: digest(value) for name, value in data.items()
    }
    check(
        "content_digests",
        meta["content_sha256"] == content_digests,
        "all theorem, table, pinned, corollary, witness, and scope sections are digest-bound",
    )
    check(
        "producer_resource_declarations",
        meta["process_cpu_seconds"] <= meta["process_cpu_budget_seconds"] <= CPU_BUDGET_SECONDS
        and meta["peak_rss_bytes"] <= meta["rss_limit_bytes"] == RSS_LIMIT_BYTES,
    )

    elapsed = time.process_time() - started
    check(
        "verifier_resource_limits",
        elapsed <= CPU_BUDGET_SECONDS and peak_rss_bytes() <= RSS_LIMIT_BYTES,
        f"process CPU={elapsed:.6f}s; peak RSS={peak_rss_bytes()} bytes",
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: standalone exact open-box genus-bound verification")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
