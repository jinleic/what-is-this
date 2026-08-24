#!/usr/bin/env python3
"""Exact finite controls for the open-box orientable-genus lower bound.

The mathematical payload proves an all-size theorem for
``P_a square P_b square P_c`` with open boundary conditions.  The computation
uses only exact integer graph data.  It checks the counting, parity,
unit-plaquette, ceiling, and classification arithmetic on every ordered box
with side lengths from 1 through 8, and records explicit witnesses for two
pinned boxes.

This file is the producer.  ``tests/test_box_genus_bound.py`` is deliberately
standalone and does not import it.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import sys
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "topology" / "box_genus_bound.json"
VERIFIER = ROOT / "tests" / "test_box_genus_bound.py"
CPU_BUDGET_SECONDS = 60.0
RSS_LIMIT_BYTES = 2 * 1024**3
TABLE_SIDE_MAXIMUM = 8

Coordinate = tuple[int, int, int]
Edge = tuple[Coordinate, Coordinate]
Plaquette = tuple[int, int, Coordinate]

AXIS_NAMES = ("x", "y", "z")
AXIS_PAIRS = ((0, 1), (0, 2), (1, 2))


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


def dimensions(a: int, b: int, c: int) -> tuple[int, int, int]:
    result = (a, b, c)
    if any(type(side) is not int or side < 1 for side in result):
        raise ValueError("open-box side lengths must be positive integers")
    return result


def canonical_edge(left: Coordinate, right: Coordinate) -> Edge:
    if left == right:
        raise ValueError("a simple-box edge cannot be a loop")
    return (left, right) if left < right else (right, left)


def box_graph(a: int, b: int, c: int) -> tuple[tuple[Coordinate, ...], tuple[Edge, ...]]:
    sides = dimensions(a, b, c)
    vertices = tuple(itertools.product(*(range(side) for side in sides)))
    edges: list[Edge] = []
    for vertex in vertices:
        for axis, side in enumerate(sides):
            if vertex[axis] + 1 >= side:
                continue
            neighbor = list(vertex)
            neighbor[axis] += 1
            edges.append(canonical_edge(vertex, tuple(neighbor)))
    return vertices, tuple(sorted(edges))


def graph_is_connected(vertices: Sequence[Coordinate], edges: Sequence[Edge]) -> bool:
    if not vertices:
        return False
    neighbors = {vertex: [] for vertex in vertices}
    for left, right in edges:
        neighbors[left].append(right)
        neighbors[right].append(left)
    reached = {vertices[0]}
    queue = deque((vertices[0],))
    while queue:
        vertex = queue.popleft()
        for neighbor in neighbors[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    return len(reached) == len(vertices)


def graph_is_simple(edges: Sequence[Edge]) -> bool:
    return len(edges) == len(set(edges)) and all(left != right for left, right in edges)


def parity_data(vertices: Sequence[Coordinate], edges: Sequence[Edge]) -> dict[str, object]:
    part_zero = sum(sum(vertex) % 2 == 0 for vertex in vertices)
    part_one = len(vertices) - part_zero
    crosses = all((sum(left) - sum(right)) % 2 == 1 for left, right in edges)
    return {
        "part_sizes_even_odd": [part_zero, part_one],
        "every_edge_crosses_parity": crosses,
    }


def enumerate_plaquettes(a: int, b: int, c: int) -> tuple[Plaquette, ...]:
    sides = dimensions(a, b, c)
    plaquettes: list[Plaquette] = []
    for first_axis, second_axis in AXIS_PAIRS:
        ranges = [
            range(side - 1) if axis in (first_axis, second_axis) else range(side)
            for axis, side in enumerate(sides)
        ]
        for base in itertools.product(*ranges):
            plaquettes.append((first_axis, second_axis, base))
    return tuple(plaquettes)


def plaquette_vertices(plaquette: Plaquette) -> tuple[Coordinate, Coordinate, Coordinate, Coordinate]:
    first_axis, second_axis, base = plaquette
    first = list(base)
    first[first_axis] += 1
    second = list(base)
    second[second_axis] += 1
    opposite = list(first)
    opposite[second_axis] += 1
    return base, tuple(first), tuple(opposite), tuple(second)


def plaquette_edges(plaquette: Plaquette) -> tuple[Edge, Edge, Edge, Edge]:
    vertices = plaquette_vertices(plaquette)
    return tuple(
        canonical_edge(vertices[index], vertices[(index + 1) % 4]) for index in range(4)
    )  # type: ignore[return-value]


def plaquette_is_valid(plaquette: Plaquette, sides: Sequence[int], graph_edges: set[Edge]) -> bool:
    first_axis, second_axis, base = plaquette
    vertices = plaquette_vertices(plaquette)
    return (
        0 <= first_axis < second_axis < 3
        and all(
            all(0 <= coordinate[axis] < sides[axis] for axis in range(3))
            for coordinate in vertices
        )
        and len(set(vertices)) == 4
        and all(edge in graph_edges for edge in plaquette_edges(plaquette))
        and all(
            sum(abs(left[axis] - right[axis]) for axis in range(3)) == 1
            for left, right in plaquette_edges(plaquette)
        )
        and all(base[axis] >= 0 for axis in range(3))
    )


def plaquette_witness_data(
    a: int, b: int, c: int, edges: Sequence[Edge]
) -> tuple[tuple[Plaquette, ...], dict[Edge, tuple[Plaquette, ...]], list[dict[str, object]]]:
    sides = dimensions(a, b, c)
    plaquettes = enumerate_plaquettes(*sides)
    edge_set = set(edges)
    incidence: dict[Edge, list[Plaquette]] = {edge: [] for edge in edges}
    for plaquette in plaquettes:
        if not plaquette_is_valid(plaquette, sides, edge_set):
            raise AssertionError(f"invalid unit plaquette {plaquette} in box {sides}")
        for edge in plaquette_edges(plaquette):
            incidence[edge].append(plaquette)
    frozen_incidence = {edge: tuple(sorted(rows)) for edge, rows in incidence.items()}
    selected = [
        {
            "edge": [list(edge[0]), list(edge[1])],
            "plaquette": {
                "axes": [witness[0], witness[1]],
                "base": list(witness[2]),
            },
        }
        for edge in sorted(edges)
        if (witness := frozen_incidence[edge][0] if frozen_incidence[edge] else None) is not None
    ]
    return plaquettes, frozen_incidence, selected


def ceil_div(numerator: int, denominator: int) -> int:
    if denominator <= 0:
        raise ValueError("ceil_div requires a positive denominator")
    return -((-numerator) // denominator)


def genus_numerator(a: int, b: int, c: int) -> int:
    dimensions(a, b, c)
    return a * b * c - a * b - b * c - c * a + 4


def classified_planar(a: int, b: int, c: int) -> bool:
    sides = dimensions(a, b, c)
    ordered = sorted(sides)
    return ordered[0] == 1 or ordered[:2] == [2, 2]


def classification_label(a: int, b: int, c: int) -> str:
    ordered = sorted(dimensions(a, b, c))
    if ordered[0] == 1:
        return "planar_side_one_slab"
    if ordered[:2] == [2, 2]:
        return "planar_2x2xL"
    return "nonplanar_by_positive_genus_lower_bound"


def box_control_row(a: int, b: int, c: int) -> dict[str, object]:
    sides = dimensions(a, b, c)
    vertices, edges = box_graph(*sides)
    plaquettes, incidence, selected_witnesses = plaquette_witness_data(*sides, edges)
    parity = parity_data(vertices, edges)

    axis_counts = {
        AXIS_NAMES[first] + AXIS_NAMES[second]: sum(
            plaquette[0:2] == (first, second) for plaquette in plaquettes
        )
        for first, second in AXIS_PAIRS
    }
    witnessed_edges = sum(bool(rows) for rows in incidence.values())
    theorem_applies = min(sides) >= 2
    lower_bound: dict[str, object] | None = None
    if theorem_applies:
        numerator = genus_numerator(*sides)
        quotient, remainder = divmod(numerator, 4)
        ceiling = quotient + int(remainder != 0)
        lower_bound = {
            "numerator": numerator,
            "floor_quotient": quotient,
            "remainder": remainder,
            "ceiling": ceiling,
            "genus_lower_bound": max(0, ceiling),
            "positive": max(0, ceiling) > 0,
        }

    return {
        "dims": list(sides),
        "vertices": len(vertices),
        "edges": len(edges),
        "simple": graph_is_simple(edges),
        "connected": graph_is_connected(vertices, edges),
        "bipartite": bool(parity["every_edge_crosses_parity"]),
        "bipartition_sizes_even_odd": parity["part_sizes_even_odd"],
        "unit_plaquette_axis_counts": axis_counts,
        "unit_plaquettes": len(plaquettes),
        "plaquette_edge_incidences": sum(len(rows) for rows in incidence.values()),
        "edges_with_unit_plaquette_witness": witnessed_edges,
        "all_edges_have_unit_plaquette": witnessed_edges == len(edges),
        "selected_witnesses_sha256": canonical_sha256(selected_witnesses),
        "genus_bound_theorem_applies": theorem_applies,
        "lower_bound": lower_bound,
        "classified_planar": classified_planar(*sides),
        "classification": classification_label(*sides),
    }


def explicit_witness_control(a: int, b: int, c: int) -> dict[str, object]:
    sides = dimensions(a, b, c)
    _, edges = box_graph(*sides)
    plaquettes, incidence, selected = plaquette_witness_data(*sides, edges)
    plaquette_rows = [
        {"axes": [first, second], "base": list(base)}
        for first, second, base in plaquettes
    ]
    incidence_rows = [
        {
            "edge": [list(edge[0]), list(edge[1])],
            "incident_unit_plaquettes": len(incidence[edge]),
        }
        for edge in edges
    ]
    return {
        "dims": list(sides),
        "edges": len(edges),
        "unit_plaquettes": plaquette_rows,
        "selected_edge_witnesses": selected,
        "edge_incidence_counts": incidence_rows,
        "selected_witnesses_sha256": canonical_sha256(selected),
    }


def pinned_projection(row: dict[str, object]) -> dict[str, object]:
    lower_bound = row["lower_bound"]
    return {
        "dims": row["dims"],
        "vertices": row["vertices"],
        "edges": row["edges"],
        "bipartition_sizes_even_odd": row["bipartition_sizes_even_odd"],
        "unit_plaquette_axis_counts": row["unit_plaquette_axis_counts"],
        "unit_plaquettes": row["unit_plaquettes"],
        "edges_with_unit_plaquette_witness": row["edges_with_unit_plaquette_witness"],
        "lower_bound_numerator": None if lower_bound is None else lower_bound["numerator"],
        "genus_lower_bound": None if lower_bound is None else lower_bound["genus_lower_bound"],
        "classification": row["classification"],
    }


PINNED_EXPECTED: list[dict[str, object]] = [
    {
        "dims": [1, 1, 1],
        "vertices": 1,
        "edges": 0,
        "bipartition_sizes_even_odd": [1, 0],
        "unit_plaquette_axis_counts": {"xy": 0, "xz": 0, "yz": 0},
        "unit_plaquettes": 0,
        "edges_with_unit_plaquette_witness": 0,
        "lower_bound_numerator": None,
        "genus_lower_bound": None,
        "classification": "planar_side_one_slab",
    },
    {
        "dims": [1, 1, 8],
        "vertices": 8,
        "edges": 7,
        "bipartition_sizes_even_odd": [4, 4],
        "unit_plaquette_axis_counts": {"xy": 0, "xz": 0, "yz": 0},
        "unit_plaquettes": 0,
        "edges_with_unit_plaquette_witness": 0,
        "lower_bound_numerator": None,
        "genus_lower_bound": None,
        "classification": "planar_side_one_slab",
    },
    {
        "dims": [1, 2, 3],
        "vertices": 6,
        "edges": 7,
        "bipartition_sizes_even_odd": [3, 3],
        "unit_plaquette_axis_counts": {"xy": 0, "xz": 0, "yz": 2},
        "unit_plaquettes": 2,
        "edges_with_unit_plaquette_witness": 7,
        "lower_bound_numerator": None,
        "genus_lower_bound": None,
        "classification": "planar_side_one_slab",
    },
    {
        "dims": [2, 1, 2],
        "vertices": 4,
        "edges": 4,
        "bipartition_sizes_even_odd": [2, 2],
        "unit_plaquette_axis_counts": {"xy": 0, "xz": 1, "yz": 0},
        "unit_plaquettes": 1,
        "edges_with_unit_plaquette_witness": 4,
        "lower_bound_numerator": None,
        "genus_lower_bound": None,
        "classification": "planar_side_one_slab",
    },
    {
        "dims": [2, 2, 2],
        "vertices": 8,
        "edges": 12,
        "bipartition_sizes_even_odd": [4, 4],
        "unit_plaquette_axis_counts": {"xy": 2, "xz": 2, "yz": 2},
        "unit_plaquettes": 6,
        "edges_with_unit_plaquette_witness": 12,
        "lower_bound_numerator": 0,
        "genus_lower_bound": 0,
        "classification": "planar_2x2xL",
    },
    {
        "dims": [2, 2, 8],
        "vertices": 32,
        "edges": 60,
        "bipartition_sizes_even_odd": [16, 16],
        "unit_plaquette_axis_counts": {"xy": 8, "xz": 14, "yz": 14},
        "unit_plaquettes": 36,
        "edges_with_unit_plaquette_witness": 60,
        "lower_bound_numerator": 0,
        "genus_lower_bound": 0,
        "classification": "planar_2x2xL",
    },
    {
        "dims": [2, 3, 3],
        "vertices": 18,
        "edges": 33,
        "bipartition_sizes_even_odd": [9, 9],
        "unit_plaquette_axis_counts": {"xy": 6, "xz": 6, "yz": 8},
        "unit_plaquettes": 20,
        "edges_with_unit_plaquette_witness": 33,
        "lower_bound_numerator": 1,
        "genus_lower_bound": 1,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
    {
        "dims": [2, 3, 8],
        "vertices": 48,
        "edges": 98,
        "bipartition_sizes_even_odd": [24, 24],
        "unit_plaquette_axis_counts": {"xy": 16, "xz": 21, "yz": 28},
        "unit_plaquettes": 65,
        "edges_with_unit_plaquette_witness": 98,
        "lower_bound_numerator": 6,
        "genus_lower_bound": 2,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
    {
        "dims": [3, 3, 3],
        "vertices": 27,
        "edges": 54,
        "bipartition_sizes_even_odd": [14, 13],
        "unit_plaquette_axis_counts": {"xy": 12, "xz": 12, "yz": 12},
        "unit_plaquettes": 36,
        "edges_with_unit_plaquette_witness": 54,
        "lower_bound_numerator": 4,
        "genus_lower_bound": 1,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
    {
        "dims": [3, 3, 8],
        "vertices": 72,
        "edges": 159,
        "bipartition_sizes_even_odd": [36, 36],
        "unit_plaquette_axis_counts": {"xy": 32, "xz": 42, "yz": 42},
        "unit_plaquettes": 116,
        "edges_with_unit_plaquette_witness": 159,
        "lower_bound_numerator": 19,
        "genus_lower_bound": 5,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
    {
        "dims": [4, 5, 6],
        "vertices": 120,
        "edges": 286,
        "bipartition_sizes_even_odd": [60, 60],
        "unit_plaquette_axis_counts": {"xy": 72, "xz": 75, "yz": 80},
        "unit_plaquettes": 227,
        "edges_with_unit_plaquette_witness": 286,
        "lower_bound_numerator": 50,
        "genus_lower_bound": 13,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
    {
        "dims": [8, 8, 8],
        "vertices": 512,
        "edges": 1344,
        "bipartition_sizes_even_odd": [256, 256],
        "unit_plaquette_axis_counts": {"xy": 392, "xz": 392, "yz": 392},
        "unit_plaquettes": 1176,
        "edges_with_unit_plaquette_witness": 1344,
        "lower_bound_numerator": 324,
        "genus_lower_bound": 81,
        "classification": "nonplanar_by_positive_genus_lower_bound",
    },
]


def theorem_payload() -> dict[str, object]:
    return {
        "name": "open simple-cubic box orientable-genus lower bound and planarity classification",
        "domain": (
            "G(a,b,c)=P_a square P_b square P_c for positive integers a,b,c, with open "
            "boundary conditions and no periodic identifications"
        ),
        "counts": {
            "vertices": "V=abc",
            "axis_edge_counts": ["(a-1)bc", "a(b-1)c", "ab(c-1)"],
            "edges": "E=3abc-(ab+bc+ca)",
            "derivation": (
                "Choose a starting coordinate for each positive coordinate-direction edge; "
                "the three disjoint direction classes have the displayed sizes."
            ),
        },
        "properties_for_minimum_side_at_least_two": {
            "simple": (
                "Every edge changes exactly one coordinate by one, so it has distinct endpoints; "
                "its endpoints and changed coordinate determine it uniquely."
            ),
            "connected": (
                "Move one coordinate at a time along its path to join any two vertices."
            ),
            "bipartite": (
                "The parity of x+y+z is a bipartition, since every edge flips it."
            ),
            "bridgeless": (
                "For an edge in one coordinate direction, choose either transverse coordinate. "
                "That side has length at least two, so the current transverse coordinate has an "
                "in-range neighbor. Shifting the edge to that neighbor supplies the other side of "
                "a unit plaquette. The remaining three plaquette edges join the endpoints after "
                "the original edge is removed, so the edge is not a bridge."
            ),
        },
        "cellularity_lemma": {
            "statement": (
                "Every finite connected graph has a cellular embedding in a closed orientable "
                "surface whose genus equals its minimum embedding genus."
            ),
            "proof": [
                "Start with any minimum-genus embedding in a closed orientable surface S_g and take a connected regular neighborhood N of the embedded graph.",
                "Build N from small vertex disks and narrow edge bands. Its boundary components are precisely the boundary walks of the induced ribbon graph. Capping each boundary component by a disk therefore gives a closed orientable surface S_h in which the same graph is cellular, with one disk face per cap.",
                "For completeness, let N have genus h and q boundary components. If the closures of the r components of S_g minus N have genera k_i and boundary counts q_i, then sum q_i=q and Euler-characteristic additivity along boundary circles gives g=h+sum_i k_i+q-r. Each complementary component has a boundary, so q>=r and h<=g.",
                "The capped surface S_h is another surface containing the graph, so minimality of g gives g<=h. Thus h=g, proving cellularity without an external theorem.",
            ],
        },
        "face_and_euler_argument": {
            "face_length": (
                "In the cellular minimum embedding, a face boundary of length one would be a "
                "loop. A length-two boundary would use parallel edges or traverse a bridge out and "
                "back. A length-three boundary is an odd closed walk. Simplicity, bridgelessness, "
                "and bipartiteness exclude these cases, so every face has length at least four."
            ),
            "incidence": "2E=sum_f length(f)>=4F, hence F<=E/2",
            "euler": "2-2 gamma=V-E+F<=V-E/2",
            "bound": "gamma(G)>=max(0,ceil((abc-ab-bc-ca+4)/4))",
        },
        "planarity_classification": {
            "statement": (
                "G(a,b,c) is planar exactly when at least one side is one or, up to permutation, "
                "the dimensions are 2x2xL."
            ),
            "planar_constructions": [
                "If a side is one, the graph is a rectangular path grid (or a path or point), embedded in the plane in its coordinate grid.",
                "For 2x2xL, each layer is a 4-cycle. Draw the L layer cycles as successively nested squares and join corresponding corners by four disjoint radial segments in each annulus.",
            ],
            "nonplanar_proof": [
                "Order 2<=a<=b<=c. If a=2, the lower-bound numerator factors as (b-2)(c-2), which is positive unless b=2.",
                "If a>=3, the numerator is coordinatewise nondecreasing for side lengths at least two because increasing a changes it by bc-b-c=(b-1)(c-1)-1>=0, and cyclically. Its value at 3x3x3 is 4.",
                "Thus every remaining box has genus at least one and is nonplanar, while the two displayed families have planar embeddings.",
            ],
        },
        "corollaries": {
            "cubes": (
                "For L>=2, gamma(P_L square P_L square P_L)>=ceil((L-2)^2(L+1)/4). "
                "The cubic boxes are planar exactly for L=1,2."
            ),
            "fixed_cross_section": (
                "For a,b,L>=2, gamma(P_a square P_b square P_L)>=max(0,ceil(((ab-a-b)L-ab+4)/4))."
            ),
            "two_by_b": (
                "For b,L>=2, gamma(P_2 square P_b square P_L)>=ceil((b-2)(L-2)/4)."
            ),
            "elongation": (
                "For fixed a,b>=2 other than a=b=2, the displayed lower bound grows linearly "
                "with L with positive numerator slope ab-a-b. This asserts unbounded graph genus "
                "for those open-box families, not an exact genus formula."
            ),
        },
        "scoped_relevance_note": (
            "When the classical genus-g spin-structure construction is invoked, that construction "
            "has 2^(2g)=4^g terms. This is only a term count internal to that standard surface "
            "construction; it is not used in this theorem and is not a lower bound on the number "
            "of Pfaffians in an arbitrary formula."
        ),
    }


def scope_payload() -> dict[str, object]:
    return {
        "proved": [
            "the graph properties, cellularity reduction, genus lower bound, and planarity classification for all positive open-box side lengths in the stated domains",
            "cube and fixed-cross-section lower-bound corollaries",
        ],
        "computed_exactly": [
            "all 8^3=512 ordered boxes with side lengths 1 through 8",
            "explicit selected unit-plaquette witnesses for 2x2x2 and 2x3x3",
            "cube lower-bound rows through side length 12",
            "fixed-cross-section elongated rows through length 16",
        ],
        "explicitly_excluded": [
            "periodic length-two multigraphs and every other periodic identification",
            "exact orientable genera beyond the cases proved planar",
            "lower bounds on the number of Pfaffians required by arbitrary formulas",
            "a thermodynamic solution, free energy, critical point, or critical exponent",
        ],
        "benchmark_used": False,
    }


def cube_controls() -> list[dict[str, int]]:
    rows = []
    for side in range(2, 13):
        direct = genus_numerator(side, side, side)
        factored = (side - 2) ** 2 * (side + 1)
        rows.append(
            {
                "side": side,
                "direct_numerator": direct,
                "factored_numerator": factored,
                "genus_lower_bound": max(0, ceil_div(direct, 4)),
            }
        )
    return rows


def elongated_controls() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for a, b in ((2, 2), (2, 3), (2, 4), (3, 3), (3, 4), (4, 5)):
        slope = a * b - a - b
        length_rows = []
        for length in range(2, 17):
            direct = genus_numerator(a, b, length)
            affine = slope * length - a * b + 4
            length_rows.append(
                {
                    "length": length,
                    "direct_numerator": direct,
                    "affine_numerator": affine,
                    "genus_lower_bound": max(0, ceil_div(direct, 4)),
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


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    table = [
        box_control_row(a, b, c)
        for a in range(1, TABLE_SIDE_MAXIMUM + 1)
        for b in range(1, TABLE_SIDE_MAXIMUM + 1)
        for c in range(1, TABLE_SIDE_MAXIMUM + 1)
    ]
    row_by_dims = {tuple(row["dims"]): row for row in table}
    add_check(
        checks,
        "broad_table_complete",
        len(table) == TABLE_SIDE_MAXIMUM**3 and len(row_by_dims) == len(table),
        f"{len(table)} ordered boxes, sides 1..{TABLE_SIDE_MAXIMUM}",
    )

    count_formulas_hold = all(
        row["vertices"] == a * b * c
        and row["edges"] == 3 * a * b * c - (a * b + b * c + c * a)
        for row in table
        for a, b, c in (tuple(row["dims"]),)
    )
    add_check(
        checks,
        "vertex_and_edge_formulas",
        count_formulas_hold,
        "V=abc and E=3abc-(ab+bc+ca) on all 512 rows",
    )

    structural_rows_hold = all(
        row["simple"] and row["connected"] and row["bipartite"] for row in table
    )
    parity_sizes_hold = all(
        row["bipartition_sizes_even_odd"]
        == [
            (a * b * c + int(a % 2 == b % 2 == c % 2 == 1)) // 2,
            (a * b * c - int(a % 2 == b % 2 == c % 2 == 1)) // 2,
        ]
        for row in table
        for a, b, c in (tuple(row["dims"]),)
    )
    add_check(
        checks,
        "simple_connected_bipartite",
        structural_rows_hold and parity_sizes_hold,
        "coordinate construction and parity class sizes verified on every table row",
    )

    plaquette_counts_hold = all(
        row["unit_plaquette_axis_counts"]
        == {
            "xy": (a - 1) * (b - 1) * c,
            "xz": (a - 1) * b * (c - 1),
            "yz": a * (b - 1) * (c - 1),
        }
        and row["unit_plaquettes"]
        == sum(row["unit_plaquette_axis_counts"].values())
        and row["plaquette_edge_incidences"] == 4 * row["unit_plaquettes"]
        for row in table
        for a, b, c in (tuple(row["dims"]),)
    )
    witness_theorem_rows_hold = all(
        row["all_edges_have_unit_plaquette"]
        and row["edges_with_unit_plaquette_witness"] == row["edges"]
        for row in table
        if row["genus_bound_theorem_applies"]
    )
    add_check(
        checks,
        "unit_plaquette_enumeration_and_witnesses",
        plaquette_counts_hold and witness_theorem_rows_hold,
        "all plaquette counts and incidences exact; every edge witnessed whenever min(a,b,c)>=2",
    )

    ceiling_rows_hold = all(
        lower is not None
        and lower["numerator"] == genus_numerator(a, b, c)
        and lower["numerator"] == 4 * lower["floor_quotient"] + lower["remainder"]
        and 0 <= lower["remainder"] < 4
        and lower["ceiling"] == ceil_div(lower["numerator"], 4)
        and lower["genus_lower_bound"] == max(0, lower["ceiling"])
        for row in table
        if row["genus_bound_theorem_applies"]
        for a, b, c in (tuple(row["dims"]),)
        for lower in (row["lower_bound"],)
    )
    synthetic_ceiling_controls = [
        {
            "numerator": numerator,
            "ceiling_over_four": ceil_div(numerator, 4),
            "max_zero": max(0, ceil_div(numerator, 4)),
        }
        for numerator in range(-12, 21)
    ]
    synthetic_ceiling_hold = all(
        4 * (row["ceiling_over_four"] - 1) < row["numerator"]
        <= 4 * row["ceiling_over_four"]
        for row in synthetic_ceiling_controls
    )
    add_check(
        checks,
        "exact_ceiling_arithmetic",
        ceiling_rows_hold and synthetic_ceiling_hold,
        "Euclidean quotient/remainder checked on theorem rows and ceil(n/4) checked for -12<=n<=20",
    )

    classification_holds = all(
        row["classified_planar"]
        == (min(a, b, c) == 1 or sorted((a, b, c))[:2] == [2, 2])
        and (
            min(a, b, c) == 1
            or (
                row["lower_bound"]["positive"]
                == (not row["classified_planar"])
            )
        )
        for row in table
        for a, b, c in (tuple(row["dims"]),)
    )
    add_check(
        checks,
        "planarity_classification_arithmetic",
        classification_holds,
        "positive lower bound separates every nonplanar class row from slabs and 2x2xL rows",
    )

    pinned = [pinned_projection(row_by_dims[tuple(expected["dims"])]) for expected in PINNED_EXPECTED]
    add_check(
        checks,
        "pinned_rows",
        pinned == PINNED_EXPECTED,
        "12 exact boundary, planar, first-nonplanar, cubic, elongated, and maximum rows",
    )

    cubes = cube_controls()
    add_check(
        checks,
        "cube_corollary",
        all(
            row["direct_numerator"] == row["factored_numerator"]
            and row["genus_lower_bound"]
            == max(0, ceil_div(row["factored_numerator"], 4))
            for row in cubes
        )
        and [row["genus_lower_bound"] for row in cubes[:7]] == [0, 1, 5, 14, 28, 50, 81],
        "L=2..8 lower bounds are 0,1,5,14,28,50,81",
    )

    elongated = elongated_controls()
    add_check(
        checks,
        "elongated_family_corollaries",
        all(
            length_row["direct_numerator"] == length_row["affine_numerator"]
            for cross_section in elongated
            for length_row in cross_section["length_rows"]
        )
        and all(
            length_row["direct_numerator"]
            == (cross_section["cross_section"][1] - 2) * (length_row["length"] - 2)
            for cross_section in elongated
            if cross_section["cross_section"][0] == 2
            for length_row in cross_section["length_rows"]
        ),
        "affine numerator and the 2xb factorization checked through L=16",
    )

    witness_controls = [explicit_witness_control(2, 2, 2), explicit_witness_control(2, 3, 3)]
    add_check(
        checks,
        "explicit_witness_controls",
        [len(control["selected_edge_witnesses"]) for control in witness_controls] == [12, 33]
        and all(
            all(row["incident_unit_plaquettes"] >= 1 for row in control["edge_incidence_counts"])
            for control in witness_controls
        ),
        "all 12 cube edges and all 33 first-nonplanar-box edges carry explicit plaquette witnesses",
    )

    theorem = theorem_payload()
    scope = scope_payload()
    exclusions = scope["explicitly_excluded"]
    add_check(
        checks,
        "scope_exclusions",
        len(exclusions) == 4
        and any("periodic length-two multigraphs" in item for item in exclusions)
        and any("exact orientable genera beyond" in item for item in exclusions)
        and any("arbitrary formulas" in item for item in exclusions)
        and any("thermodynamic solution" in item for item in exclusions)
        and "not a lower bound" in theorem["scoped_relevance_note"],
        "periodic multigraphs, unknown exact genera, arbitrary Pfaffian counts, and thermodynamics excluded",
    )

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

    content = {
        "theorem": theorem,
        "finite_table": {
            "side_range": [1, TABLE_SIDE_MAXIMUM],
            "ordered_row_count": len(table),
            "rows": table,
        },
        "pinned_rows": pinned,
        "synthetic_ceiling_controls": synthetic_ceiling_controls,
        "cube_controls": cubes,
        "elongated_controls": elongated,
        "explicit_unit_plaquette_witness_controls": witness_controls,
        "scope": scope,
    }
    content_sha256 = {
        name: canonical_sha256(value) for name, value in content.items()
    }
    source_paths = (Path(__file__).resolve(), VERIFIER)
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e241_box_genus_bound.py",
            "verifier": "tests/test_box_genus_bound.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers only; no floating point",
            "single_process": True,
            "process_cpu_seconds": round(elapsed, 6),
            "process_cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "peak_rss_bytes": peak_rss,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "source_sha256": {
                str(path.relative_to(ROOT)): file_sha256(path) for path in source_paths
            },
            "content_sha256": content_sha256,
        },
        "data": content,
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
