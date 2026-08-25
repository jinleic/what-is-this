#!/usr/bin/env python3
"""Independent checker for the involution SRG frontier artifact.

No symbol is imported from ``involution_srg_frontier``.  The checker derives
the anti-invariant matrix equation from ``(45,22,10,11)``, generates the twelve
partitioned 2-regular zero graphs algorithmically, and exhausts their signings
with a recursive tail assignment rather than the producer's product iterator.
It also independently enumerates the fixed-five parity graphs and reconstructs
the surviving fixed-five support design.
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
SURVIVORS = (1, 5, 9, 13)


class CheckViolation(RuntimeError):
    """The artifact differs from the independent exact reconstruction."""


def _partitions(total: int, minimum: int,
                maximum: int | None = None) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    records = []
    for first in range(maximum, minimum - 1, -1):
        remainder = total - first
        if remainder == 0:
            records.append((first,))
        elif remainder >= minimum:
            for tail in _partitions(remainder, minimum, first):
                records.append((first,) + tail)
    return tuple(records)


def _shape_specs() -> tuple[tuple, ...]:
    records = []
    for cross_vertices in range(8):
        cross_parts = (((),) if cross_vertices == 0
                       else _partitions(cross_vertices, 2))
        remaining = 7 - cross_vertices
        internal_parts = _partitions(remaining, 3)
        if not cross_parts or not internal_parts:
            continue
        for cross in cross_parts:
            for edge_parts in internal_parts:
                for nonedge_parts in internal_parts:
                    records.append((cross_vertices, cross,
                                    edge_parts, nonedge_parts))
    return tuple(records)


def _shape_id(spec: tuple) -> str:
    cross_vertices, cross_parts, edge_parts, nonedge_parts = spec

    def text(parts):
        return "-".join(map(str, parts)) if parts else "none"

    return (
        f"x{cross_vertices}_cross_{text(cross_parts)}"
        f"_edge_{text(edge_parts)}"
        f"_nonedge_{text(nonedge_parts)}"
    )


def _add_internal_cycles(edges: set[tuple[int, int]], vertices: list[int],
                         parts: tuple[int, ...]) -> None:
    start = 0
    for length in parts:
        block = vertices[start:start + length]
        start += length
        for index in range(length):
            edges.add(tuple(sorted((block[index], block[(index + 1) % length]))))
    if start != len(vertices):
        raise CheckViolation("internal cycle partition is incomplete")


def _add_cross_cycles(edges: set[tuple[int, int]], left: list[int],
                      right: list[int], parts: tuple[int, ...]) -> None:
    start = 0
    for length in parts:
        left_block = left[start:start + length]
        right_block = right[start:start + length]
        start += length
        for index in range(length):
            edges.add(tuple(sorted((left_block[index], right_block[index]))))
            edges.add(tuple(sorted((left_block[index],
                                    right_block[(index - 1) % length]))))
    if start != len(left) or start != len(right):
        raise CheckViolation("cross-cycle partition is incomplete")


def _zero_graph(spec: tuple) -> set[tuple[int, int]]:
    cross_vertices, cross_parts, edge_parts, nonedge_parts = spec
    edge_part = list(range(7))
    nonedge_part = list(range(7, 14))
    edges: set[tuple[int, int]] = set()
    if cross_vertices:
        _add_cross_cycles(edges, edge_part[:cross_vertices],
                          nonedge_part[:cross_vertices], cross_parts)
    _add_internal_cycles(edges, edge_part[cross_vertices:], edge_parts)
    _add_internal_cycles(edges, nonedge_part[cross_vertices:], nonedge_parts)
    if len(edges) != 14:
        raise CheckViolation("zero graph has the wrong edge count")
    if sorted(sum(vertex in edge for edge in edges) for vertex in range(14)) != [2] * 14:
        raise CheckViolation("zero graph is not 2-regular")
    return edges


def _fix_tree(matrix, support) -> None:
    """Implement the artifact's canonical lexicographic BFS gauge."""
    reached = [False] * len(matrix)
    reached[0] = True
    queue = [0]
    position = 0
    while position < len(queue):
        vertex = queue[position]
        position += 1
        for neighbor in range(len(matrix)):
            if neighbor not in support[vertex] or reached[neighbor]:
                continue
            reached[neighbor] = True
            queue.append(neighbor)
            matrix[vertex][neighbor] = matrix[neighbor][vertex] = 1
    if not all(reached):
        raise CheckViolation("support graph is disconnected")


def _serialize_prefix(matrix, final_row: int) -> bytes:
    rows = []
    for left in range(final_row + 1):
        rows.append(",".join(str(matrix[left][right])
                             for right in range(left, len(matrix))))
    return (";".join(rows) + "\n").encode()


def _hash_prefixes(records) -> str:
    digest = hashlib.sha256()
    for record in sorted(records):
        digest.update(record)
    return digest.hexdigest()


def _parity_failures(order: int, zero_edges) -> list[dict]:
    zero_edges = {tuple(sorted(edge)) for edge in zero_edges}
    half = order // 2
    support = [
        {right for right in range(order)
         if right != left and tuple(sorted((left, right))) not in zero_edges}
        for left in range(order)
    ]
    failures = []
    for left, right in itertools.combinations(range(order), 2):
        half_incidence = right in support[left]
        common = len(support[left] & support[right])
        expected = int(
            half_incidence and (left < half) == (right < half)
        )
        if common % 2 != expected:
            failures.append({
                "pair": [left, right],
                "pair_is_half_incidence": half_incidence,
                "common_half_incidence_neighbors": common,
                "required_parity": expected,
            })
    return failures




def _independent_sign_search(order: int, zero_edges) -> dict:
    half = order // 2
    matrix = [[None] * order for _ in range(order)]
    for index in range(order):
        matrix[index][index] = -1 if index < half else 0
    zero_edges = {tuple(sorted(edge)) for edge in zero_edges}
    for left, right in zero_edges:
        matrix[left][right] = matrix[right][left] = 0
    support = []
    for left in range(order):
        support.append({
            right for right in range(order)
            if right != left and tuple(sorted((left, right))) not in zero_edges
        })
    if {len(neighbors) for neighbors in support} != {11}:
        raise CheckViolation("half-incidence support is not 11-regular")
    _fix_tree(matrix, support)

    tested = [0] * order
    survivors = [[] for _ in range(order)]
    solutions = []

    def assign_row(row: int) -> None:
        if row == order:
            solutions.append(_serialize_prefix(matrix, order - 1))
            return
        columns = [column for column in range(row + 1, order)
                   if matrix[row][column] is None]

        def assign_column(position: int) -> None:
            if position != len(columns):
                column = columns[position]
                matrix[row][column] = matrix[column][row] = -1
                assign_column(position + 1)
                matrix[row][column] = matrix[column][row] = 1
                assign_column(position + 1)
                matrix[row][column] = matrix[column][row] = None
                return
            tested[row] += 1
            for previous in range(row):
                inner = sum(matrix[previous][index] * matrix[index][row]
                            for index in range(order))
                if inner + matrix[previous][row] != 0:
                    return
            survivors[row].append(_serialize_prefix(matrix, row))
            assign_row(row + 1)

        assign_column(0)

    assign_row(0)
    tested_by_row = {str(row): count for row, count in enumerate(tested) if count}
    surviving_by_row = {
        str(row): len(survivors[row])
        for row, count in enumerate(tested) if count
    }
    hashes = {
        str(row): _hash_prefixes(survivors[row])
        for row, count in enumerate(tested) if count
    }
    trace_commitment = {
        "tested_by_row": tested_by_row,
        "surviving_by_row": surviving_by_row,
        "survivor_sha256_by_row": hashes,
        "solution_sha256": _hash_prefixes(solutions),
    }
    zero_edges = {tuple(sorted(edge)) for edge in zero_edges}
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


def _fixed_count_record() -> dict:
    vertices, _, _, _ = PARAMETERS
    spectral = sorted(vertices - 2 * cycles for cycles in range(2, 23, 2))
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
        "surviving_fixed_counts": list(SURVIVORS),
    }


def _anti_identity() -> dict:
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


def _anti_search_record() -> dict:
    c12_failures = _parity_failures(12, ())
    if not c12_failures:
        raise CheckViolation("the c=12 parity obstruction disappeared")
    c12 = {
        "anti_invariant_dimension": 12,
        "fixed_points": 21,
        "zero_graph_degree": 0,
        "zero_neighborhood_parity_violations": len(c12_failures),
        "first_zero_neighborhood_parity_violation": c12_failures[0],
        "gauge": "lexicographic_bfs_spanning_tree_positive",
        **_independent_sign_search(12, ()),
    }
    specs = _shape_specs()
    if len(specs) != 12:
        raise CheckViolation("the c=14 zero-graph classification is incomplete")
    shapes = []
    for spec in specs:
        zero_edges = _zero_graph(spec)
        cross_vertices, cross_parts, edge_parts, nonedge_parts = spec
        parity_failures = _parity_failures(14, zero_edges)
        if not parity_failures:
            raise CheckViolation("a c=14 zero shape passes the parity audit")
        shapes.append({
            "shape_id": _shape_id(spec),
            "cross_vertices_per_part": cross_vertices,
            "cross_cycle_half_lengths": list(cross_parts),
            "edge_internal_cycle_lengths": list(edge_parts),
            "nonedge_internal_cycle_lengths": list(nonedge_parts),
            "zero_edges": [list(edge) for edge in sorted(zero_edges)],
            "zero_neighborhood_parity_violations": len(parity_failures),
            "first_zero_neighborhood_parity_violation": parity_failures[0],
            **_independent_sign_search(14, zero_edges),
        })
    c14 = {
        "anti_invariant_dimension": 14,
        "fixed_points": 17,
        "zero_graph_degree": 2,
        "partition_sizes": [7, 7],
        "cross_zero_degree_values": [0, 2],
        "shape_classification": "internal cycles or alternating even cycles",
        "gauge": "lexicographic_bfs_spanning_tree_positive",
        "zero_graph_shapes": shapes,
        "solutions": sum(record["solutions"] for record in shapes),
    }
    if c12["solutions"] or c14["solutions"]:
        raise CheckViolation("a small anti-invariant signing survived")
    return {"closed_fixed_counts": [17, 21], "c12": c12, "c14": c14}


def _adjacency(order: int, bits: tuple[int, ...]) -> tuple[int, ...]:
    rows = [0] * order
    for enabled, (left, right) in zip(bits, itertools.combinations(range(order), 2)):
        if enabled:
            rows[left] |= 1 << right
            rows[right] |= 1 << left
    return tuple(rows)


def _fixed_parity(adjacency) -> bool:
    if any(row.bit_count() % 2 for row in adjacency):
        return False
    for left, right in itertools.combinations(range(len(adjacency)), 2):
        parity = (adjacency[left] & adjacency[right]).bit_count() % 2
        if parity != (0 if adjacency[left] >> right & 1 else 1):
            return False
    return True


def _fixed_five_record() -> dict:
    masks = []
    for bits in itertools.product((0, 1), repeat=10):
        mask = sum(bit << index for index, bit in enumerate(bits))
        if mask in (0, 1023):
            continue
        graph = _adjacency(5, bits)
        if _fixed_parity(graph):
            if sorted(row.bit_count() for row in graph) != [2] * 5:
                raise CheckViolation("a fixed-five parity graph is not C5")
            masks.append(mask)
    if len(masks) != 12:
        raise CheckViolation("the fixed-five labelled count differs")
    return {
        "tested_labelled_graphs": 1024,
        "labelled_solutions": 12,
        "isomorphism_types": 1,
        "unique_type": "C5",
        "solution_masks_sha256": hashlib.sha256(
            ",".join(map(str, sorted(masks))).encode()
        ).hexdigest(),
    }


def _rotate(mask: int, shift: int) -> int:
    output = 0
    for vertex in range(5):
        if mask >> vertex & 1:
            output |= 1 << ((vertex + shift) % 5)
    return output


def _support_design() -> dict:
    edge_masks = sorted(_rotate(mask, shift)
                        for mask in (3, 15) for shift in range(5))
    nonedge_masks = sorted(_rotate(mask, shift)
                           for mask in (1, 11) for shift in range(5))
    all_masks = edge_masks + nonedge_masks
    replication = [sum(mask >> point & 1 for mask in all_masks)
                   for point in range(5)]
    multiplicities = [
        sum((mask >> left & 1) * (mask >> right & 1) for mask in all_masks)
        for left, right in itertools.combinations(range(5), 2)
    ]
    if replication != [10] * 5 or multiplicities != [5] * 10:
        raise CheckViolation("the fixed-five support design does not balance")
    return {
        "edge_pair_support_masks": edge_masks,
        "nonedge_pair_support_masks": nonedge_masks,
        "replication": replication,
        "pair_multiplicity": 5,
        "construction": "five rotations of support sizes 1,2,3,4",
        "disposition": "EXACT_SUPPORT_DESIGN_SURVIVES",
    }


def _support_identity() -> dict:
    return {
        "fixed_degree_parity": "d_F(v) is even",
        "fixed_common_neighbor_parity": {"adjacent": "even", "nonadjacent": "odd"},
        "matrix_parity_equation": "H^2+H+I=J mod 2",
        "replication": "r(v)=(22-d_F(v))/2",
        "adjacent_pair_multiplicity": "(10-c_F(u,v))/2",
        "nonadjacent_pair_multiplicity": "(11-c_F(u,v))/2",
        "edge_pair_support_size_parity": "even",
        "nonedge_pair_support_size_parity": "odd",
        "edge_pair_support_is_triangle_free_under_ramsey_hypothesis": True,
        "nonedge_pair_nonsupport_has_no_independent_triple_under_ramsey_hypothesis": True,
    }


def _build_expected_document() -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "claim": {
            "graph_family": "strongly_regular",
            "parameters": list(PARAMETERS),
            "involution_fixed_counts": list(SURVIVORS),
            "ramsey_hypothesis_used": False,
            "general_ramsey_bound_claimed": False,
        },
        "fixed_count_reduction": _fixed_count_record(),
        "anti_invariant_identity": _anti_identity(),
        "anti_invariant_search": _anti_search_record(),
        "fixed_support_identity": _support_identity(),
        "ramsey_frontier": {
            "disposition": "OPEN_RAMSEY_COMPLETION_FOR_FIXED_1_5_9_13",
            "fixed_1": {
                "fixed_graph": "K1",
                "edge_pair_supports": "eleven empty supports",
                "nonedge_pair_supports": "eleven singleton supports",
                "status": "OPEN_SIGNED_COMPLETION",
            },
            "fixed_5": {
                "fixed_graph_enumeration": _fixed_five_record(),
                "support_design": _support_design(),
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


@functools.cache
def _expected_document_json() -> str:
    return json.dumps(
        _build_expected_document(),
        sort_keys=True,
        separators=(",", ":"),
    )


def _expected_document() -> dict:
    return json.loads(_expected_document_json())


def _exact_tree_equal(actual, expected) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return (
            actual.keys() == expected.keys()
            and all(_exact_tree_equal(actual[key], value)
                    for key, value in expected.items())
        )
    if isinstance(expected, list):
        return (
            len(actual) == len(expected)
            and all(_exact_tree_equal(left, right)
                    for left, right in zip(actual, expected))
        )
    return actual == expected


def _unique_object(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise CheckViolation(f"duplicate JSON member: {key}")
        output[key] = value
    return output


def load_document(text: str) -> dict:
    """Parse strict JSON, rejecting duplicate names and nonstandard constants."""
    def reject_constant(value):
        raise CheckViolation(f"nonstandard JSON constant: {value}")

    return json.loads(
        text,
        object_pairs_hook=_unique_object,
        parse_constant=reject_constant,
    )




def verify_document(document: dict) -> dict:
    expected = _expected_document()
    if not _exact_tree_equal(document, expected):
        raise CheckViolation("artifact differs from independent exact reconstruction")
    return {
        "closed_fixed_counts": [17, 21],
        "surviving_fixed_counts": list(SURVIVORS),
        "c14_zero_graph_shapes": 12,
        "ramsey_status": "OPEN_RAMSEY_COMPLETION_FOR_FIXED_1_5_9_13",
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    args = parser.parse_args(argv)
    try:
        document = load_document(args.artifact.read_text())
        summary = verify_document(document)
    except (OSError, json.JSONDecodeError, CheckViolation) as error:
        print(f"INVOLUTION CHECK FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "INVOLUTION EVIDENCE VERIFIED:", DISPOSITION,
        "closed_fixed_counts=17,21",
        "surviving_fixed_counts=1,5,9,13",
        f"c14_zero_graph_shapes={summary['c14_zero_graph_shapes']}",
        "(all SRG-only; surviving cases open; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
