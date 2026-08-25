#!/usr/bin/env python3
"""Exact involution fixed-count frontier for ``srg(45,22,10,11)``.

Theorem
=======

Every nonidentity involution of a strongly regular graph with parameters
``(45,22,10,11)`` fixes exactly 1, 5, 9, or 13 vertices.

Let the involution have ``f`` fixed vertices and ``c=(45-f)/2`` transposition
orbits.  On the anti-invariant basis ``e_x-e_tau(x)``, adjacency restricts to
an integral symmetric matrix ``C`` with

``C^2 + C = 11 I``.

Its diagonal has ``c/2`` entries ``-1`` and ``c/2`` entries ``0``; its
off-diagonal entries are ``-1,0,1``.  Every row has exactly eleven nonzero
off-diagonal entries.  The spectral multiplicities first give ``f=1 mod 4``;
the row count gives ``c>=12`` and hence leaves ``f=1,5,9,13,17,21``.

For ``c=12`` the off-diagonal support is complete and no signing satisfies the
matrix equation.  For ``c=14`` the zero graph is 2-regular.  Reduction modulo
two forces every vertex to have either zero or two cross-part zero-neighbors,
so its cycles have exactly twelve partition-preserving forms.  A switching
gauge fixes a spanning tree positive; exhaustive row-prefix enumeration finds
no signing for any form.  Thus fixed counts 17 and 21 are impossible.

This theorem uses only the SRG equations.  It does not exclude the four
surviving fixed counts, does not constrain non-strongly-regular graphs, and
does not change a bound on ``R(5,5)``.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import sys
from pathlib import Path

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_srg_frontier"
DISPOSITION = "INVOLUTION_FIXED_COUNTS_REDUCED_TO_1_5_9_13_FOR_SRG"
PARAMETERS = (45, 22, 10, 11)
SURVIVING_FIXED_COUNTS = (1, 5, 9, 13)

# ``(cross vertices per part, cross-cycle half-lengths,
#    edge-part internal cycle lengths, nonedge-part internal cycle lengths)``.
# Parts contain seven vertices each.  A cycle part is at least three for an
# internal cycle and at least two for an alternating cross cycle.
_C14_SHAPE_SPECS = (
    (0, (), (7,), (7,)),
    (0, (), (7,), (4, 3)),
    (0, (), (4, 3), (7,)),
    (0, (), (4, 3), (4, 3)),
    (2, (2,), (5,), (5,)),
    (3, (3,), (4,), (4,)),
    (4, (4,), (3,), (3,)),
    (4, (2, 2), (3,), (3,)),
    (7, (7,), (), ()),
    (7, (5, 2), (), ()),
    (7, (4, 3), (), ()),
    (7, (3, 2, 2), (), ()),
)


class FrontierViolation(RuntimeError):
    """An exact involution invariant or finite exhaustion failed."""


def fixed_count_reduction() -> dict:
    """Derive the spectral and half-incidence fixed-count candidates."""
    vertices = PARAMETERS[0]
    spectral = [vertices - 2 * cycles for cycles in range(2, 23, 2)]
    spectral.sort()
    half_incidence = [fixed for fixed in spectral
                      if (vertices - fixed) // 2 >= 12]
    return {
        "anti_invariant_dimension": "c=(45-f)/2",
        "irreducible_polynomial": "x^2+x-11",
        "dimension_parity": "c is even",
        "fixed_count_congruence": "f=1 mod 4",
        "spectral_candidates": spectral,
        "half_incidence_degree": 11,
        "half_incidence_candidates": half_incidence,
        "surviving_fixed_counts": list(SURVIVING_FIXED_COUNTS),
    }


def anti_invariant_identity() -> dict:
    """Record the anti-invariant matrix equation and its graph meaning."""
    return {
        "matrix_equation": "C^2+C=11I",
        "diagonal_values": [-1, 0],
        "off_diagonal_values": [-1, 0, 1],
        "diagonal_minus_one_count": "c/2",
        "diagonal_zero_count": "c/2",
        "edge_pair_orbits": "c/2",
        "nonedge_pair_orbits": "c/2",
        "half_incidence_blocks_per_row": 11,
        "zero_graph_degree": "c-12",
        "switching_action": "C -> DCD for diagonal D in {+1,-1}",
    }


def _cycle_edges(vertices: list[int], parts: tuple[int, ...]) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    position = 0
    for size in parts:
        block = vertices[position:position + size]
        position += size
        if size < 3:
            raise FrontierViolation("an internal zero cycle has length below three")
        for index, vertex in enumerate(block):
            edges.add(tuple(sorted((vertex, block[(index + 1) % size]))))
    if position != len(vertices):
        raise FrontierViolation("internal zero cycles do not cover their vertices")
    return edges


def _cross_cycle_edges(edge_vertices: list[int], nonedge_vertices: list[int],
                       parts: tuple[int, ...]) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    position = 0
    for size in parts:
        left = edge_vertices[position:position + size]
        right = nonedge_vertices[position:position + size]
        position += size
        if size < 2:
            raise FrontierViolation("an alternating zero cycle has length below four")
        for index, vertex in enumerate(left):
            edges.add(tuple(sorted((vertex, right[index]))))
            edges.add(tuple(sorted((vertex, right[(index - 1) % size]))))
    if position != len(edge_vertices) or position != len(nonedge_vertices):
        raise FrontierViolation("alternating zero cycles do not cover their vertices")
    return edges


def _zero_edges_from_spec(spec: tuple) -> set[tuple[int, int]]:
    cross_vertices, cross_parts, edge_parts, nonedge_parts = spec
    edge_vertices = list(range(7))
    nonedge_vertices = list(range(7, 14))
    zero_edges: set[tuple[int, int]] = set()
    if cross_vertices:
        zero_edges.update(_cross_cycle_edges(
            edge_vertices[:cross_vertices],
            nonedge_vertices[:cross_vertices],
            cross_parts,
        ))
    remaining = 7 - cross_vertices
    if remaining:
        zero_edges.update(_cycle_edges(edge_vertices[cross_vertices:], edge_parts))
        zero_edges.update(_cycle_edges(
            nonedge_vertices[cross_vertices:], nonedge_parts))
    if len(zero_edges) != 14:
        raise FrontierViolation("a c=14 zero graph does not have fourteen edges")
    if any(sum(vertex in edge for edge in zero_edges) != 2
           for vertex in range(14)):
        raise FrontierViolation("a c=14 zero graph is not 2-regular")
    return zero_edges


def _shape_id(spec: tuple) -> str:
    cross_vertices, cross_parts, edge_parts, nonedge_parts = spec

    def part_text(parts):
        return "-".join(map(str, parts)) if parts else "none"

    return (
        f"x{cross_vertices}_cross_{part_text(cross_parts)}"
        f"_edge_{part_text(edge_parts)}"
        f"_nonedge_{part_text(nonedge_parts)}"
    )


def _positive_spanning_tree(matrix: list[list[int | None]],
                            support: list[set[int]]) -> None:
    """Gauge-fix the lexicographic breadth-first spanning tree positive."""
    seen = {0}
    queue = [0]
    while queue:
        vertex = queue.pop(0)
        for neighbor in sorted(support[vertex] - seen):
            matrix[vertex][neighbor] = matrix[neighbor][vertex] = 1
            seen.add(neighbor)
            queue.append(neighbor)
    if len(seen) != len(matrix):
        raise FrontierViolation("the half-incidence support graph is disconnected")


def _prefix_serialization(matrix: list[list[int | None]], row: int) -> bytes:
    pieces = []
    order = len(matrix)
    for left in range(row + 1):
        pieces.append(",".join(str(matrix[left][right])
                               for right in range(left, order)))
    return (";".join(pieces) + "\n").encode()


def _digest_records(records: list[bytes]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records):
        digest.update(record)
    return digest.hexdigest()


def _zero_neighborhood_parity_violations(
        order: int, zero_edges: set[tuple[int, int]]) -> list[dict]:
    """Find sign-independent parity failures in ``C^2+C=11I``."""
    zero_edges = {tuple(sorted(edge)) for edge in zero_edges}
    half = order // 2
    violations = []
    for left, right in itertools.combinations(range(order), 2):
        pair = (left, right)
        pair_is_half_incidence = pair not in zero_edges
        common_half_incidence = sum(
            tuple(sorted((left, vertex))) not in zero_edges
            and tuple(sorted((right, vertex))) not in zero_edges
            for vertex in range(order)
            if vertex not in pair
        )
        required_parity = int(
            pair_is_half_incidence
            and (left < half) == (right < half)
        )
        if common_half_incidence % 2 != required_parity:
            violations.append({
                "pair": [left, right],
                "pair_is_half_incidence": pair_is_half_incidence,
                "common_half_incidence_neighbors": common_half_incidence,
                "required_parity": required_parity,
            })
    return violations




def _signed_anti_search(order: int,
                        zero_edges: set[tuple[int, int]]) -> dict:
    """Exhaust signings after a canonical switching-tree gauge."""
    if order % 2:
        raise FrontierViolation("the anti-invariant dimension is odd")
    half = order // 2
    matrix: list[list[int | None]] = [[None] * order for _ in range(order)]
    for index in range(order):
        matrix[index][index] = -1 if index < half else 0
    for left, right in zero_edges:
        matrix[left][right] = matrix[right][left] = 0

    support = [
        {right for right in range(order)
         if right != left and tuple(sorted((left, right))) not in zero_edges}
        for left in range(order)
    ]
    if any(len(neighbors) != 11 for neighbors in support):
        raise FrontierViolation("a half-incidence row does not have degree eleven")
    _positive_spanning_tree(matrix, support)

    tested = [0] * order
    surviving_prefixes: list[list[bytes]] = [[] for _ in range(order)]
    solutions: list[bytes] = []

    def visit(row: int) -> None:
        if row == order:
            solutions.append(_prefix_serialization(matrix, order - 1))
            return
        tail = [column for column in range(row + 1, order)
                if matrix[row][column] is None]
        for signs in itertools.product((-1, 1), repeat=len(tail)):
            tested[row] += 1
            for column, sign in zip(tail, signs):
                matrix[row][column] = matrix[column][row] = sign
            if all(
                sum(matrix[previous][index] * matrix[index][row]
                    for index in range(order)) + matrix[previous][row] == 0
                for previous in range(row)
            ):
                surviving_prefixes[row].append(_prefix_serialization(matrix, row))
                visit(row + 1)
            for column in tail:
                matrix[row][column] = matrix[column][row] = None

    visit(0)
    tested_by_row = {str(row): count for row, count in enumerate(tested) if count}
    surviving_by_row = {
        str(row): len(records)
        for row, records in enumerate(surviving_prefixes)
        if tested[row]
    }
    survivor_sha256_by_row = {
        str(row): _digest_records(surviving_prefixes[row])
        for row, count in enumerate(tested)
        if count
    }
    trace_commitment = {
        "tested_by_row": tested_by_row,
        "surviving_by_row": surviving_by_row,
        "survivor_sha256_by_row": survivor_sha256_by_row,
        "solution_sha256": _digest_records(solutions),
    }
    search_domain = {
        "format_version": 1,
        "order": order,
        "diagonal_minus_one_indices": list(range(half)),
        "zero_edges": [list(edge) for edge in sorted(zero_edges)],
        "gauge": "lexicographic_bfs_spanning_tree_positive",
        "prefix_serialization": "upper_triangular_row_csv_v1",
    }
    search_sha256 = hashlib.sha256(json.dumps({
        "domain": search_domain,
        "trace": trace_commitment,
    }, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {
        **trace_commitment,
        "search_domain": search_domain,
        "search_sha256": search_sha256,
        "solutions": len(solutions),
    }


def _build_small_anti_obstructions() -> dict:
    """Close the c=12 and c=14 anti-invariant matrix cases exactly."""
    c12_parity = _zero_neighborhood_parity_violations(12, set())
    if not c12_parity:
        raise FrontierViolation("the c=12 parity obstruction disappeared")
    c12 = {
        "anti_invariant_dimension": 12,
        "fixed_points": 21,
        "zero_graph_degree": 0,
        "zero_neighborhood_parity_violations": len(c12_parity),
        "first_zero_neighborhood_parity_violation": c12_parity[0],
        "gauge": "lexicographic_bfs_spanning_tree_positive",
        **_signed_anti_search(12, set()),
    }
    shape_records = []
    for spec in _C14_SHAPE_SPECS:
        zero_edges = _zero_edges_from_spec(spec)
        cross_vertices, cross_parts, edge_parts, nonedge_parts = spec
        parity_violations = _zero_neighborhood_parity_violations(14, zero_edges)
        if not parity_violations:
            raise FrontierViolation("a c=14 zero shape passes the parity audit")
        shape_records.append({
            "shape_id": _shape_id(spec),
            "cross_vertices_per_part": cross_vertices,
            "cross_cycle_half_lengths": list(cross_parts),
            "edge_internal_cycle_lengths": list(edge_parts),
            "nonedge_internal_cycle_lengths": list(nonedge_parts),
            "zero_edges": [list(edge) for edge in sorted(zero_edges)],
            "zero_neighborhood_parity_violations": len(parity_violations),
            "first_zero_neighborhood_parity_violation": parity_violations[0],
            **_signed_anti_search(14, zero_edges),
        })
    c14 = {
        "anti_invariant_dimension": 14,
        "fixed_points": 17,
        "zero_graph_degree": 2,
        "partition_sizes": [7, 7],
        "cross_zero_degree_values": [0, 2],
        "shape_classification":
            "internal cycles or alternating even cycles",
        "gauge": "lexicographic_bfs_spanning_tree_positive",
        "zero_graph_shapes": shape_records,
        "solutions": sum(record["solutions"] for record in shape_records),
    }
    if c12["solutions"] or c14["solutions"]:
        raise FrontierViolation("a supposedly closed involution quotient survived")
    return {
        "closed_fixed_counts": [17, 21],
        "c12": c12,
        "c14": c14,
    }


@functools.cache
def _small_anti_obstructions_json() -> str:
    return json.dumps(
        _build_small_anti_obstructions(),
        sort_keys=True,
        separators=(",", ":"),
    )


def enumerate_small_anti_obstructions() -> dict:
    """Return a fresh copy of the cached exact anti-invariant exhaustion."""
    return json.loads(_small_anti_obstructions_json())




def _graph_from_upper_mask(order: int, mask: int) -> tuple[int, ...]:
    adjacency = [0] * order
    bit = 0
    for left in range(order):
        for right in range(left + 1, order):
            if mask >> bit & 1:
                adjacency[left] |= 1 << right
                adjacency[right] |= 1 << left
            bit += 1
    return tuple(adjacency)


def _fixed_graph_parity_holds(adjacency: tuple[int, ...]) -> bool:
    order = len(adjacency)
    if any(row.bit_count() % 2 for row in adjacency):
        return False
    for left in range(order):
        for right in range(left + 1, order):
            common_parity = (adjacency[left] & adjacency[right]).bit_count() % 2
            expected = 0 if adjacency[left] >> right & 1 else 1
            if common_parity != expected:
                return False
    return True


def _build_fixed_five_graph_enumeration() -> dict:
    """Enumerate the fixed graph when an involution has five fixed points."""
    solutions = []
    for mask in range(1 << 10):
        adjacency = _graph_from_upper_mask(5, mask)
        if mask in (0, (1 << 10) - 1):
            continue
        if _fixed_graph_parity_holds(adjacency):
            solutions.append(mask)
    if len(solutions) != 12:
        raise FrontierViolation("the fixed-five parity enumeration changed")
    if any(sorted(row.bit_count() for row in _graph_from_upper_mask(5, mask))
           != [2, 2, 2, 2, 2] for mask in solutions):
        raise FrontierViolation("a fixed-five graph is not a 5-cycle")
    digest = hashlib.sha256(
        ",".join(map(str, sorted(solutions))).encode()
    ).hexdigest()
    return {
        "tested_labelled_graphs": 1 << 10,
        "labelled_solutions": len(solutions),
        "isomorphism_types": 1,
        "unique_type": "C5",
        "solution_masks_sha256": digest,
    }


@functools.cache
def _fixed_five_graph_enumeration_json() -> str:
    return json.dumps(
        _build_fixed_five_graph_enumeration(),
        sort_keys=True,
        separators=(",", ":"),
    )


def fixed_five_graph_enumeration() -> dict:
    """Return a fresh fixed-five graph certificate."""
    return json.loads(_fixed_five_graph_enumeration_json())




def _rotate_mask(mask: int, shift: int, order: int = 5) -> int:
    return sum(1 << ((vertex + shift) % order)
               for vertex in range(order) if mask >> vertex & 1)


def _build_fixed_five_support_design() -> dict:
    """Give an exact support design showing that the f=5 packing survives."""
    edge_bases = (0b00011, 0b01111)
    nonedge_bases = (0b00001, 0b01011)
    edge_masks = sorted(_rotate_mask(mask, shift)
                        for mask in edge_bases for shift in range(5))
    nonedge_masks = sorted(_rotate_mask(mask, shift)
                           for mask in nonedge_bases for shift in range(5))
    masks = edge_masks + nonedge_masks
    replication = [sum(mask >> vertex & 1 for mask in masks)
                   for vertex in range(5)]
    pair_counts = [
        sum((mask >> left & 1) and (mask >> right & 1) for mask in masks)
        for left, right in itertools.combinations(range(5), 2)
    ]
    if replication != [10] * 5 or pair_counts != [5] * 10:
        raise FrontierViolation("the fixed-five support design is invalid")
    if any(mask.bit_count() % 2 for mask in edge_masks):
        raise FrontierViolation("an edge-pair support has odd size")
    if any(mask.bit_count() % 2 != 1 for mask in nonedge_masks):
        raise FrontierViolation("a nonedge-pair support has even size")
    return {
        "edge_pair_support_masks": edge_masks,
        "nonedge_pair_support_masks": nonedge_masks,
        "replication": replication,
        "pair_multiplicity": 5,
        "construction": "five rotations of support sizes 1,2,3,4",
        "disposition": "EXACT_SUPPORT_DESIGN_SURVIVES",
    }


@functools.cache
def _fixed_five_support_design_json() -> str:
    return json.dumps(
        _build_fixed_five_support_design(),
        sort_keys=True,
        separators=(",", ":"),
    )


def fixed_five_support_design() -> dict:
    """Return a fresh fixed-five support certificate."""
    return json.loads(_fixed_five_support_design_json())




def support_packing_identity() -> dict:
    """Record the exact binary Gram constraints on fixed supports."""
    return {
        "fixed_degree_parity": "d_F(v) is even",
        "fixed_common_neighbor_parity": {
            "adjacent": "even",
            "nonadjacent": "odd",
        },
        "matrix_parity_equation": "H^2+H+I=J mod 2",
        "replication": "r(v)=(22-d_F(v))/2",
        "adjacent_pair_multiplicity": "(10-c_F(u,v))/2",
        "nonadjacent_pair_multiplicity": "(11-c_F(u,v))/2",
        "edge_pair_support_size_parity": "even",
        "nonedge_pair_support_size_parity": "odd",
        "edge_pair_support_is_triangle_free_under_ramsey_hypothesis": True,
        "nonedge_pair_nonsupport_has_no_independent_triple_under_ramsey_hypothesis": True,
    }


def run_analysis() -> dict:
    fixed_counts = fixed_count_reduction()
    anti_search = enumerate_small_anti_obstructions()
    if fixed_counts["surviving_fixed_counts"] != list(SURVIVING_FIXED_COUNTS):
        raise FrontierViolation("the involution fixed-count frontier changed")
    fixed_five_graph = fixed_five_graph_enumeration()
    fixed_five_design = fixed_five_support_design()
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "claim": {
            "graph_family": "strongly_regular",
            "parameters": list(PARAMETERS),
            "involution_fixed_counts": list(SURVIVING_FIXED_COUNTS),
            "ramsey_hypothesis_used": False,
            "general_ramsey_bound_claimed": False,
        },
        "fixed_count_reduction": fixed_counts,
        "anti_invariant_identity": anti_invariant_identity(),
        "anti_invariant_search": anti_search,
        "fixed_support_identity": support_packing_identity(),
        "ramsey_frontier": {
            "disposition": "OPEN_RAMSEY_COMPLETION_FOR_FIXED_1_5_9_13",
            "fixed_1": {
                "fixed_graph": "K1",
                "edge_pair_supports": "eleven empty supports",
                "nonedge_pair_supports": "eleven singleton supports",
                "status": "OPEN_SIGNED_COMPLETION",
            },
            "fixed_5": {
                "fixed_graph_enumeration": fixed_five_graph,
                "support_design": fixed_five_design,
                "status": "PRIMARY_OPEN_SIGNED_COMPLETION_TARGET",
            },
            "fixed_9": {
                "status": "OPEN_FIXED_GRAPH_SUPPORT_AND_SIGNED_COMPLETION"
            },
            "fixed_13": {
                "status": "OPEN_FIXED_GRAPH_SUPPORT_AND_SIGNED_COMPLETION"
            },
            "non_strongly_regular_graphs_constrained": False,
            "general_ramsey_bound_claimed": False,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" /
        "involution_srg_frontier.json",
    )
    args = parser.parse_args(argv)
    try:
        document = run_analysis()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    except (OSError, FrontierViolation) as error:
        print(f"INVOLUTION FRONTIER FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "INVOLUTION FRONTIER:",
        "closed_fixed_counts=17,21",
        "surviving_fixed_counts=1,5,9,13",
        "c14_zero_graph_shapes=12",
    )
    print(
        "INVOLUTION RESULT:", DISPOSITION,
        "(all srg(45,22,10,11); surviving cases open; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
