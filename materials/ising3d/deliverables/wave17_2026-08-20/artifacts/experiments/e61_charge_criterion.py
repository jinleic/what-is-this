#!/usr/bin/env python3
"""Exact graph criteria and conceptual forms for the 3x3 support-seven charges.

The calculation has three parts.  It identifies the four exceptional kernels in
``flux_support.json`` with sparse rational permutation/Bell-projector formulas;
it tests the pointwise-stabilizer criterion on every stored support orbit and on
all proper connected support orbits of connected graphs through six vertices;
and it searches larger rectangles using graph automorphisms only.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
FLUX_RESULT = ROOT / "results" / "integrability" / "flux_support.json"
RESULT = ROOT / "results" / "integrability" / "charge_criterion.json"
SCRIPT = "experiments/e61_charge_criterion.py"
PRIMES = (2_147_483_647, 2_147_483_629)

Edge = tuple[int, int]
Permutation = tuple[int, ...]
Pauli = tuple[int, int]
Operator = dict[Pauli, Fraction]
Column = dict[tuple[int, int], int]


def popcount(mask: int) -> int:
    return bin(mask).count("1")


def grid_edges(rows: int, columns: int) -> tuple[Edge, ...]:
    return tuple(
        [(r * columns + c, (r + 1) * columns + c) for r in range(rows - 1) for c in range(columns)]
        + [(r * columns + c, r * columns + c + 1) for r in range(rows) for c in range(columns - 1)]
    )


def edge_set(edges: Iterable[Edge]) -> set[Edge]:
    return {tuple(sorted(edge)) for edge in edges}


def is_graph_automorphism(n: int, edges: tuple[Edge, ...], permutation: Permutation) -> bool:
    original = edge_set(edges)
    return {tuple(sorted((permutation[u], permutation[v]))) for u, v in edges} == original


def grid_automorphisms(rows: int, columns: int) -> tuple[Permutation, ...]:
    transforms = [
        lambda r, c: (r, c),
        lambda r, c: (rows - 1 - r, c),
        lambda r, c: (r, columns - 1 - c),
        lambda r, c: (rows - 1 - r, columns - 1 - c),
    ]
    if rows == columns:
        transforms.extend(
            [
                lambda r, c: (c, r),
                lambda r, c: (columns - 1 - c, rows - 1 - r),
                lambda r, c: (c, columns - 1 - r),
                lambda r, c: (rows - 1 - c, r),
            ]
        )
    edges = grid_edges(rows, columns)
    permutations = []
    for transform in transforms:
        permutation = tuple(
            transform(site // columns, site % columns)[0] * columns
            + transform(site // columns, site % columns)[1]
            for site in range(rows * columns)
        )
        if permutation not in permutations and is_graph_automorphism(rows * columns, edges, permutation):
            permutations.append(permutation)
    return tuple(permutations)


def transform_mask(mask: int, permutation: Permutation) -> int:
    return sum(1 << permutation[v] for v in range(len(permutation)) if (mask >> v) & 1)


def is_connected_support(n: int, edges: tuple[Edge, ...], support_mask: int) -> bool:
    if not support_mask:
        return False
    adjacency = [0] * n
    for u, v in edges:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    seen = support_mask & -support_mask
    frontier = seen
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        vertex = bit.bit_length() - 1
        new = adjacency[vertex] & support_mask & ~seen
        seen |= new
        frontier |= new
    return seen == support_mask


def pointwise_stabilizer(automorphisms: Iterable[Permutation], n: int, support_mask: int) -> tuple[Permutation, ...]:
    return tuple(
        permutation
        for permutation in automorphisms
        if all(permutation[v] == v for v in range(n) if not ((support_mask >> v) & 1))
    )


def moved_mask(permutation: Permutation) -> int:
    return sum(1 << v for v, image in enumerate(permutation) if image != v)


def op_clean(operator: Operator) -> Operator:
    return {key: Fraction(value) for key, value in operator.items() if value}


def op_add(*summands: tuple[Fraction | int, Operator]) -> Operator:
    result: defaultdict[Pauli, Fraction] = defaultdict(Fraction)
    for scale, operator in summands:
        scale_q = Fraction(scale)
        for key, value in operator.items():
            result[key] += scale_q * value
    return op_clean(dict(result))


def op_mul(left: Operator, right: Operator) -> Operator:
    result: defaultdict[Pauli, Fraction] = defaultdict(Fraction)
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            sign = -1 if popcount(b & c) % 2 else 1
            result[(a ^ c, b ^ d)] += sign * x * y
    return op_clean(dict(result))


def op_transpose(operator: Operator) -> Operator:
    return {
        (a, b): value * (-1 if popcount(a & b) % 2 else 1)
        for (a, b), value in operator.items()
    }


def identity_operator() -> Operator:
    return {(0, 0): Fraction(1)}


def pair_operator(u: int, v: int, kind: str) -> Operator:
    pair = (1 << u) | (1 << v)
    if kind == "swap":
        coefficients = (1, 1, 1, -1)
        denominator = 2
    elif kind == "singlet":
        coefficients = (1, -1, -1, 1)
        denominator = 4
    elif kind == "phi_minus":
        coefficients = (1, -1, 1, -1)
        denominator = 4
    else:
        raise ValueError(kind)
    keys = ((0, 0), (pair, 0), (0, pair), (pair, pair))
    return {key: Fraction(value, denominator) for key, value in zip(keys, coefficients)}


def product(operators: Iterable[Operator]) -> Operator:
    result = identity_operator()
    for operator in operators:
        result = op_mul(result, operator)
    return result


def x_times(vertex: int, operator: Operator) -> Operator:
    return op_mul({(1 << vertex, 0): Fraction(1)}, operator)


def decode_stored_operator(vector: dict[str, object]) -> Operator:
    return {
        (int(a), int(b)): Fraction(int(numerator), int(denominator))
        for a, b, numerator, denominator in vector["terms_a_b_num_den"]
    }


def encode_operator(operator: Operator) -> list[list[int]]:
    return [
        [a, b, coefficient.numerator, coefficient.denominator]
        for (a, b), coefficient in sorted(operator.items())
    ]


def operator_digest(operator: Operator) -> str:
    payload = json.dumps(encode_operator(operator), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def commutator_residual(operator: Operator, n: int, edges: tuple[Edge, ...]) -> tuple[Operator, Operator]:
    ad_a: defaultdict[Pauli, Fraction] = defaultdict(Fraction)
    ad_b: defaultdict[Pauli, Fraction] = defaultdict(Fraction)
    for (a_mask, b_mask), coefficient in operator.items():
        for site in range(n):
            if (b_mask >> site) & 1:
                ad_a[(a_mask ^ (1 << site), b_mask)] += coefficient
        for u, v in edges:
            if ((a_mask >> u) ^ (a_mask >> v)) & 1:
                ad_b[(a_mask, b_mask ^ (1 << u) ^ (1 << v))] -= coefficient
    return op_clean(dict(ad_a)), op_clean(dict(ad_b))


def projector_rank(operator: Operator, local_sites: int) -> int:
    trace = (1 << local_sites) * operator.get((0, 0), Fraction(0))
    if trace.denominator != 1:
        raise ArithmeticError(f"nonintegral projector trace {trace}")
    return trace.numerator


def projector_partition(projectors: list[Operator], local_sites: int) -> tuple[bool, list[int]]:
    identity = identity_operator()
    complete = op_add(*[(1, projector) for projector in projectors]) == identity
    idempotent = all(op_mul(projector, projector) == projector for projector in projectors)
    symmetric = all(op_transpose(projector) == projector for projector in projectors)
    orthogonal = all(
        not op_mul(projectors[i], projectors[j])
        for i in range(len(projectors))
        for j in range(len(projectors))
        if i != j
    )
    ranks = [projector_rank(projector, local_sites) for projector in projectors]
    return complete and idempotent and symmetric and orthogonal and sum(ranks) == 1 << local_sites, ranks


def cycles_of_involution(permutation: Permutation) -> tuple[list[int], list[tuple[int, int]]]:
    fixed = [v for v, image in enumerate(permutation) if image == v]
    pairs = sorted((v, permutation[v]) for v in range(len(permutation)) if v < permutation[v])
    if any(permutation[permutation[v]] != v for v in range(len(permutation))):
        raise ValueError("expected an involution")
    return fixed, pairs


def bell_assignments(
    edges: tuple[Edge, ...], fixed: list[int], pairs: list[tuple[int, int]]
) -> tuple[list[tuple[int, ...]], list[list[int]], list[list[int]]]:
    edges_q = edge_set(edges)
    anchors: list[list[int]] = []
    for u, v in pairs:
        anchors.append([w for w in fixed if tuple(sorted((u, w))) in edges_q or tuple(sorted((v, w))) in edges_q])
    pair_links: list[list[int]] = []
    for i, (u, v) in enumerate(pairs):
        linked = []
        for j, (x, y) in enumerate(pairs):
            if i < j and any(tuple(sorted(edge)) in edges_q for edge in ((u, x), (u, y), (v, x), (v, y))):
                linked.append(j)
        pair_links.append(linked)
    assignments = []
    for eta in itertools.product((-1, 1), repeat=len(pairs)):
        if any(anchors[i] and eta[i] != -1 for i in range(len(pairs))):
            continue
        if any(eta[i] * eta[j] != -1 for i, linked in enumerate(pair_links) for j in linked):
            continue
        assignments.append(tuple(eta))
    return assignments, anchors, pair_links


def individual_swap_permutation(n: int, pair: tuple[int, int]) -> Permutation:
    permutation = list(range(n))
    permutation[pair[0]], permutation[pair[1]] = permutation[pair[1]], permutation[pair[0]]
    return tuple(permutation)


def induced_c4s(sites: list[int], edges: tuple[Edge, ...]) -> list[list[int]]:
    original = edge_set(edges)
    cycles = []
    for vertices in itertools.combinations(sorted(sites), 4):
        internal = [edge for edge in original if edge[0] in vertices and edge[1] in vertices]
        degrees = {vertex: 0 for vertex in vertices}
        for u, v in internal:
            degrees[u] += 1
            degrees[v] += 1
        if len(internal) == 4 and set(degrees.values()) == {2}:
            cycles.append(list(vertices))
    return cycles


def rank_q(rows: list[list[Fraction]], columns: int) -> int:
    pivots: dict[int, dict[int, Fraction]] = {}
    for source in rows:
        row = {column: Fraction(value) for column, value in enumerate(source) if value}
        while row:
            lead = min(row)
            old = pivots.get(lead)
            if old is None:
                scale = row[lead]
                pivots[lead] = {column: value / scale for column, value in row.items()}
                break
            scale = row[lead]
            for column, value in old.items():
                reduced = row.get(column, Fraction(0)) - scale * value
                if reduced:
                    row[column] = reduced
                else:
                    row.pop(column, None)
    return len(pivots)


def span_intersection_dimension(basis: list[Operator], allowed_mask: int) -> int:
    keys = sorted(set().union(*(operator.keys() for operator in basis)))
    constraints = [
        [operator.get(key, Fraction(0)) for operator in basis]
        for key in keys
        if (key[0] | key[1]) & ~allowed_mask
    ]
    return len(basis) - rank_q(constraints, len(basis))


def analyze_exception(entry: dict[str, object], automorphisms: tuple[Permutation, ...], edges: tuple[Edge, ...]) -> dict[str, object]:
    support_mask = int(entry["support_mask"])
    sites = [int(site) for site in entry["sites"]]
    stabilizer = pointwise_stabilizer(automorphisms, 9, support_mask)
    nonidentity = [permutation for permutation in stabilizer if permutation != tuple(range(9))]
    if len(nonidentity) != 1:
        raise ArithmeticError(f"unexpected stabilizer for {support_mask}: {len(stabilizer)}")
    reflection = nonidentity[0]
    fixed, pairs = cycles_of_involution(reflection)
    swap = product(pair_operator(u, v, "swap") for u, v in pairs)
    stored_basis = [decode_stored_operator(vector) for vector in entry["exact_basis"]]
    assignments, anchors, pair_links = bell_assignments(edges, fixed, pairs)
    identity = identity_operator()
    conceptual_basis: list[Operator]
    projectors: list[Operator]
    projector_names: list[str]
    formulas: list[str]
    bell_data: dict[str, object]
    if assignments:
        if len(assignments) != 1:
            raise ArithmeticError(f"nonunique Bell assignment for {support_mask}: {assignments}")
        eta = assignments[0]
        bell = product(
            pair_operator(u, v, "singlet" if sign == -1 else "phi_minus")
            for (u, v), sign in zip(pairs, eta)
        )
        supported_fixed = [vertex for vertex in fixed if (support_mask >> vertex) & 1]
        eligible_fixed = [
            vertex
            for vertex in supported_fixed
            if not any(tuple(sorted((vertex, other))) in edge_set(edges) for other in fixed if other != vertex)
            and all(
                eta[index] == -1
                for index, pair in enumerate(pairs)
                if any(tuple(sorted((vertex, endpoint))) in edge_set(edges) for endpoint in pair)
            )
        ]
        if len(eligible_fixed) != 1:
            raise ArithmeticError(f"unexpected Bell fixed vertex for {support_mask}: {eligible_fixed}")
        fixed_vertex = eligible_fixed[0]
        x_bell = x_times(fixed_vertex, bell)
        conceptual_basis = [
            identity,
            op_add((4, swap), (-32, bell)),
            op_add((1, identity), (-4, swap), (-32, bell)),
            op_add((-64, x_bell)),
        ]
        reflection_minus = op_add((Fraction(1, 2), identity), (Fraction(-1, 2), swap))
        symmetric_bulk = op_add((Fraction(1, 2), identity), (Fraction(1, 2), swap), (-1, bell))
        code_plus = op_add((Fraction(1, 2), bell), (Fraction(1, 2), x_bell))
        code_minus = op_add((Fraction(1, 2), bell), (Fraction(-1, 2), x_bell))
        projectors = [reflection_minus, symmetric_bulk, code_plus, code_minus]
        projector_names = ["reflection_minus", "reflection_plus_outside_code", "bell_code_X_plus", "bell_code_X_minus"]
        formulas = ["I", "4 U - 32 P", "I - 4 U - 32 P", "-64 X_f P"]
        bell_data = {
            "claim_tag": "[THEOREM]",
            "assignment_exists": True,
            "eta_by_pair": [int(value) for value in eta],
            "fixed_vertex_f": fixed_vertex,
            "P_terms": len(bell),
            "X_f_P_terms": len(x_bell),
            "P_digest_sha256": operator_digest(bell),
            "X_f_P_digest_sha256": operator_digest(x_bell),
            "U_P_equals_P": op_mul(swap, bell) == bell,
            "basis_coefficients_in_projector_order": [
                [1, 1, 1, 1],
                [-4, 4, -28, -28],
                [5, -3, -35, -35],
                [0, 0, -64, 64],
            ],
        }
    else:
        conceptual_basis = [identity, op_add((1, identity), (-8, swap))]
        projectors = [
            op_add((Fraction(1, 2), identity), (Fraction(1, 2), swap)),
            op_add((Fraction(1, 2), identity), (Fraction(-1, 2), swap)),
        ]
        projector_names = ["reflection_plus", "reflection_minus"]
        formulas = ["I", "I - 8 U"]
        bell_data = {
            "claim_tag": "[COMPUTATION]",
            "assignment_exists": False,
            "eta_by_pair": None,
            "obstruction": "every reflection pair touches a fixed vertex, forcing eta=-1 on adjacent quotient cells that must have opposite signs",
            "basis_coefficients_in_projector_order": [[1, 1], [-7, 9]],
        }
    partition_ok, ranks = projector_partition(projectors, 7)
    c4s = induced_c4s(sites, edges)
    c4_intersections = [
        {
            "vertices": cycle,
            "kernel_span_intersection_dimension": span_intersection_dimension(
                stored_basis, sum(1 << vertex for vertex in cycle)
            ),
        }
        for cycle in c4s
    ]
    outside = [vertex for vertex in range(9) if not ((support_mask >> vertex) & 1)]
    outside_neighborhoods = []
    edges_q = edge_set(edges)
    for u, v in pairs:
        left = [w for w in outside if tuple(sorted((u, w))) in edges_q]
        right = [w for w in outside if tuple(sorted((v, w))) in edges_q]
        outside_neighborhoods.append({"pair": [u, v], "left": left, "right": right, "equal": left == right})
    return {
        "claim_tag": "[THEOREM]",
        "support_mask": support_mask,
        "sites": sites,
        "omitted_vertices": outside,
        "symmetry_orbit_size": int(entry["symmetry_orbit_size"]),
        "kernel_dimension_Q_from_flux_support": int(entry["kernel_dimension_Q"]),
        "pointwise_stabilizer_order": len(stabilizer),
        "reflection_permutation": list(reflection),
        "fixed_vertices": fixed,
        "reflection_pairs": [list(pair) for pair in pairs],
        "outside_neighborhood_cells": outside_neighborhoods,
        "individual_pair_swaps_are_automorphisms": [
            is_graph_automorphism(9, edges, individual_swap_permutation(9, pair)) for pair in pairs
        ],
        "simultaneous_pair_swap_is_automorphism": is_graph_automorphism(9, edges, reflection),
        "quotient_fixed_anchors": anchors,
        "quotient_pair_links": pair_links,
        "bell_mechanism": bell_data,
        "conceptual_basis_formulas": formulas,
        "conceptual_basis_term_counts": [len(operator) for operator in conceptual_basis],
        "stored_basis_term_counts": [len(operator) for operator in stored_basis],
        "conceptual_basis_exactly_equals_stored_basis": conceptual_basis == stored_basis,
        "conceptual_basis_digests_sha256": [operator_digest(operator) for operator in conceptual_basis],
        "all_conceptual_commutators_zero": all(
            commutator_residual(operator, 9, edges) == ({}, {}) for operator in conceptual_basis
        ),
        "global_commutant_projector_factorization": {
            "claim_tag": "[THEOREM]",
            "projector_order": projector_names,
            "local_ranks": ranks,
            "orthogonal_idempotent_partition_exact": partition_ok,
            "every_projector_commutes_with_A_and_B": all(
                commutator_residual(projector, 9, edges) == ({}, {}) for projector in projectors
            ),
            "interpretation": "after tensoring identity on the omitted vertices these are exact projectors in the full 3x3 joint commutant",
        },
        "induced_C4_subgraphs": c4s,
        "induced_C4_kernel_intersections": c4_intersections,
    }


def stored_grid_criterion(flux: dict[str, object]) -> dict[str, object]:
    rows_checked = 0
    mismatches = []
    positive = []
    per_layer = []
    for layer in flux["data"]["layers"]:
        rows, columns = map(int, layer["shape"])
        n = rows * columns
        automorphisms = grid_automorphisms(rows, columns)
        layer_positive = []
        for row in layer["table"]:
            rows_checked += 1
            support_mask = int(row["support_mask"])
            stabilizer = pointwise_stabilizer(automorphisms, n, support_mask)
            predicted = len(stabilizer) > 1
            observed = int(row["kernel_dimension_Q"]) > 1
            if predicted != observed:
                mismatches.append(
                    {
                        "graph": layer["graph"],
                        "support_mask": support_mask,
                        "stabilizer_order": len(stabilizer),
                        "kernel_dimension_Q": int(row["kernel_dimension_Q"]),
                    }
                )
            if predicted:
                record = {
                    "graph": layer["graph"],
                    "support_mask": support_mask,
                    "stabilizer_order": len(stabilizer),
                    "kernel_dimension_Q": int(row["kernel_dimension_Q"]),
                }
                positive.append(record)
                layer_positive.append(record)
        per_layer.append(
            {
                "graph": layer["graph"],
                "rows_checked": len(layer["table"]),
                "automorphism_group_order": len(automorphisms),
                "criterion_positive_rows": layer_positive,
            }
        )
    return {
        "claim_tag": "[COMPUTATION]",
        "rows_checked": rows_checked,
        "mismatches": mismatches,
        "criterion_positive_rows": positive,
        "per_layer": per_layer,
        "conclusion": "within all stored proper connected support orbits, kernel_dimension_Q>1 iff the complement has a nontrivial pointwise stabilizer",
    }


def complete_edges(n: int) -> tuple[Edge, ...]:
    return tuple((u, v) for u in range(n) for v in range(u + 1, n))


@lru_cache(maxsize=None)
def permutation_data(n: int) -> tuple[tuple[Permutation, tuple[int, ...]], ...]:
    edges = complete_edges(n)
    index = {edge: position for position, edge in enumerate(edges)}
    data = []
    for permutation in itertools.permutations(range(n)):
        edge_map = tuple(index[tuple(sorted((permutation[u], permutation[v])))] for u, v in edges)
        data.append((permutation, edge_map))
    return tuple(data)


def transform_graph_mask(graph_mask: int, edge_map: tuple[int, ...]) -> int:
    transformed = 0
    bits = graph_mask
    while bits:
        bit = bits & -bits
        source = bit.bit_length() - 1
        transformed |= 1 << edge_map[source]
        bits ^= bit
    return transformed


@lru_cache(maxsize=None)
def canonical_graph_mask(n: int, graph_mask: int) -> int:
    return min(transform_graph_mask(graph_mask, edge_map) for _, edge_map in permutation_data(n))


def graph_edges_from_mask(n: int, graph_mask: int) -> tuple[Edge, ...]:
    return tuple(edge for position, edge in enumerate(complete_edges(n)) if (graph_mask >> position) & 1)


def graph_mask_from_edges(n: int, edges: Iterable[Edge]) -> int:
    index = {edge: position for position, edge in enumerate(complete_edges(n))}
    return sum(1 << index[tuple(sorted(edge))] for edge in edges)


def extend_connected_graphs(previous: list[int], n: int) -> list[int]:
    representatives = set()
    for graph_mask in previous:
        old_edges = set(graph_edges_from_mask(n - 1, graph_mask))
        for neighborhood in range(1, 1 << (n - 1)):
            edges = old_edges | {(vertex, n - 1) for vertex in range(n - 1) if (neighborhood >> vertex) & 1}
            candidate = graph_mask_from_edges(n, edges)
            representatives.add(canonical_graph_mask(n, candidate))
    return sorted(representatives)


def graph_automorphisms(n: int, graph_mask: int) -> tuple[Permutation, ...]:
    return tuple(
        permutation
        for permutation, edge_map in permutation_data(n)
        if transform_graph_mask(graph_mask, edge_map) == graph_mask
    )


def local_columns(n: int, edges: tuple[Edge, ...], support_mask: int) -> list[Column]:
    sites = [vertex for vertex in range(n) if (support_mask >> vertex) & 1]
    incident = [edge for edge in edges if ((support_mask >> edge[0]) & 1) or ((support_mask >> edge[1]) & 1)]
    columns = []
    for column in range(1 << (2 * len(sites))):
        digits = column
        a_mask = 0
        b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a_mask |= 1 << site
            if label & 2:
                b_mask |= 1 << site
        pauli = a_mask | (b_mask << n)
        vector: Column = {}
        for site in sites:
            if (b_mask >> site) & 1:
                vector[(0, pauli ^ (1 << site))] = 1
        for u, v in incident:
            if ((a_mask >> u) ^ (a_mask >> v)) & 1:
                vector[(1, pauli ^ (1 << (n + u)) ^ (1 << (n + v)))] = -1
        columns.append(vector)
    return columns


def column_rank_mod(columns: list[Column], prime: int) -> int:
    frequencies = Counter(row for column in columns for row in column)
    order = {row: index for index, row in enumerate(sorted(frequencies, key=lambda row: (frequencies[row], row)))}
    pivots: dict[tuple[int, int], dict[tuple[int, int], int]] = {}
    rank = 0
    for source in sorted(columns, key=len):
        vector = {row: value % prime for row, value in source.items() if value % prime}
        while vector:
            lead = min(vector, key=order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = {row: value * inverse % prime for row, value in vector.items()}
                rank += 1
                break
            factor = vector[lead]
            for row, value in old.items():
                reduced = (vector.get(row, 0) - factor * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return rank


def finite_graph_scan() -> dict[str, object]:
    representatives: dict[int, list[int]] = {1: [0]}
    for n in range(2, 7):
        representatives[n] = extend_connected_graphs(representatives[n - 1], n)
    summaries = []
    case_lines = []
    mismatches = []
    prime_disagreements = []
    for n in range(2, 7):
        case_count = 0
        positive_count = 0
        trivial_count = 0
        nullity_histogram: Counter[int] = Counter()
        size_histogram: Counter[int] = Counter()
        for graph_mask in representatives[n]:
            edges = graph_edges_from_mask(n, graph_mask)
            automorphisms = graph_automorphisms(n, graph_mask)
            for support_mask in range(1, (1 << n) - 1):
                if not is_connected_support(n, edges, support_mask):
                    continue
                if support_mask != min(transform_mask(support_mask, permutation) for permutation in automorphisms):
                    continue
                stabilizer = pointwise_stabilizer(automorphisms, n, support_mask)
                columns = local_columns(n, edges, support_mask)
                dimension = len(columns)
                nullities = tuple(dimension - column_rank_mod(columns, prime) for prime in PRIMES)
                case_count += 1
                size_histogram[popcount(support_mask)] += 1
                nullity_histogram[nullities[0]] += 1
                if len(stabilizer) > 1:
                    positive_count += 1
                else:
                    trivial_count += 1
                if len(set(nullities)) != 1:
                    prime_disagreements.append([n, graph_mask, support_mask, *nullities])
                if (len(stabilizer) > 1) != (nullities[0] > 1):
                    mismatches.append([n, graph_mask, support_mask, len(stabilizer), *nullities])
                case_lines.append(
                    f"{n}:{graph_mask}:{support_mask}:{len(automorphisms)}:{len(stabilizer)}:{nullities[0]}:{nullities[1]}"
                )
        summaries.append(
            {
                "n": n,
                "connected_unlabeled_graphs": len(representatives[n]),
                "proper_connected_support_orbits": case_count,
                "nontrivial_pointwise_stabilizer_cases": positive_count,
                "trivial_pointwise_stabilizer_cases": trivial_count,
                "support_size_histogram": {str(key): size_histogram[key] for key in sorted(size_histogram)},
                "modular_nullity_histogram": {str(key): nullity_histogram[key] for key in sorted(nullity_histogram)},
            }
        )
    digest = hashlib.sha256(("\n".join(case_lines) + "\n").encode()).hexdigest()
    return {
        "claim_tag": "[COMPUTATION]",
        "scope": "all connected simple graphs on 2..6 vertices up to isomorphism and all nonempty proper connected supports up to each graph automorphism group",
        "primes": list(PRIMES),
        "per_n": summaries,
        "total_support_orbits": len(case_lines),
        "case_digest_sha256": digest,
        "prime_disagreements": prime_disagreements,
        "criterion_mismatches": mismatches,
        "exact_Q_logic": "trivial stabilizer cases have modular nullity one, so the nonzero good-prime minor plus explicit identity proves exact Q nullity one; positive cases have the non-scalar rational site-permutation supplied by the theorem",
        "finite_conclusion": "[COMPUTATION] On this finite budget, a proper connected support has a non-scalar exact-Q local joint commutant iff its complement pointwise stabilizer is nontrivial.",
        "all_size_conjecture": "[CONJECTURE] The same equivalence holds for every finite connected simple graph and every nonempty proper connected support.",
    }


def grid_mechanism_search(rows: int, columns: int) -> dict[str, object]:
    n = rows * columns
    edges = grid_edges(rows, columns)
    automorphisms = grid_automorphisms(rows, columns)
    connected_histogram: Counter[int] = Counter()
    positive_histogram: Counter[int] = Counter()
    for support_mask in range(1, 1 << n):
        if not is_connected_support(n, edges, support_mask):
            continue
        size = popcount(support_mask)
        connected_histogram[size] += 1
        if len(pointwise_stabilizer(automorphisms, n, support_mask)) > 1:
            positive_histogram[size] += 1
    nonidentity_moved_sizes = sorted(
        popcount(moved_mask(permutation))
        for permutation in automorphisms
        if permutation != tuple(range(n))
    )
    return {
        "graph": f"{rows}x{columns}_open_grid",
        "automorphism_group_order": len(automorphisms),
        "nonidentity_moved_set_sizes": nonidentity_moved_sizes,
        "connected_support_histogram": {str(key): connected_histogram[key] for key in sorted(connected_histogram)},
        "criterion_positive_histogram": {str(key): positive_histogram[key] for key in sorted(positive_histogram)},
        "minimum_connected_criterion_support": min(positive_histogram) if positive_histogram else None,
        "connected_size_7_supports": connected_histogram[7],
        "criterion_positive_size_7_supports": positive_histogram[7],
    }


def main() -> None:
    flux_bytes = FLUX_RESULT.read_bytes()
    flux = json.loads(flux_bytes)
    layers = {layer["graph"]: layer for layer in flux["data"]["layers"]}
    layer_3x3 = layers["3x3_open_grid"]
    automorphisms_3x3 = grid_automorphisms(3, 3)
    edges_3x3 = grid_edges(3, 3)
    exceptional_entries = [
        next(row for row in layer_3x3["table"] if int(row["support_mask"]) == support_mask)
        for support_mask in (239, 254, 367, 381)
    ]
    exceptional = [analyze_exception(entry, automorphisms_3x3, edges_3x3) for entry in exceptional_entries]
    stored_criterion = stored_grid_criterion(flux)
    graph_scan = finite_graph_scan()
    grid_searches = [grid_mechanism_search(*shape) for shape in ((2, 4), (3, 3), (3, 4), (4, 4), (3, 5))]
    search_by_graph = {entry["graph"]: entry for entry in grid_searches}
    scalar_c4_row = next(row for row in layer_3x3["table"] if int(row["support_mask"]) == 27)
    c4_rejection = {
        "claim_tag": "[THEOREM]",
        "not_necessary_counterexamples": [
            {"support_mask": entry["support_mask"], "kernel_dimension_Q": entry["kernel_dimension_Q_from_flux_support"]}
            for entry in exceptional
            if not entry["induced_C4_subgraphs"]
        ],
        "not_sufficient_counterexample": {
            "support_mask": 27,
            "sites": scalar_c4_row["sites"],
            "kernel_dimension_Q": scalar_c4_row["kernel_dimension_Q"],
            "induced_C4": [0, 1, 3, 4],
        },
        "only_exception_with_C4": {
            "support_mask": 254,
            "C4_intersections": next(
                entry["induced_C4_kernel_intersections"] for entry in exceptional if entry["support_mask"] == 254
            ),
            "conclusion": "each intersection is one-dimensional (identity only), so no nonidentity exceptional charge is supported on either induced C4",
        },
    }
    checks = [
        {
            "name": "four_flux_bases_rederived_exactly",
            "passed": all(entry["conceptual_basis_exactly_equals_stored_basis"] for entry in exceptional),
            "detail": "all rational Pauli coefficients agree term-for-term for masks 239, 254, 367, 381",
        },
        {
            "name": "conceptual_forms_commute_symbolically",
            "passed": all(entry["all_conceptual_commutators_zero"] for entry in exceptional),
            "detail": "ordered-Pauli ad_A and ad_B residual dictionaries are empty",
        },
        {
            "name": "global_projector_factorizations",
            "passed": all(
                entry["global_commutant_projector_factorization"]["orthogonal_idempotent_partition_exact"]
                and entry["global_commutant_projector_factorization"]["every_projector_commutes_with_A_and_B"]
                for entry in exceptional
            ),
            "detail": "local ranks are 56+70+1+1 for diagonal supports and 72+56 for axial supports",
        },
        {
            "name": "all_250_stored_orbit_rows_match_pointwise_stabilizer_criterion",
            "passed": stored_criterion["rows_checked"] == 250 and not stored_criterion["mismatches"],
            "detail": "30 + 38 + 182 stored support-orbit rows checked",
        },
        {
            "name": "finite_connected_graph_converse_budget",
            "passed": graph_scan["total_support_orbits"] == 2682
            and not graph_scan["criterion_mismatches"]
            and not graph_scan["prime_disagreements"],
            "detail": "all proper connected support orbits of all 142 connected unlabeled graphs on 2..6 vertices",
        },
        {
            "name": "no_3x4_size7_mechanism",
            "passed": search_by_graph["3x4_open_grid"]["connected_size_7_supports"] == 234
            and search_by_graph["3x4_open_grid"]["criterion_positive_size_7_supports"] == 0
            and search_by_graph["3x4_open_grid"]["minimum_connected_criterion_support"] == 9,
            "detail": "the smallest nonidentity moved set has size 8 and needs a ninth fixed-row connector",
        },
        {
            "name": "larger_grid_search",
            "passed": search_by_graph["4x4_open_grid"]["minimum_connected_criterion_support"] == 13
            and search_by_graph["3x5_open_grid"]["minimum_connected_criterion_support"] == 11
            and search_by_graph["4x4_open_grid"]["criterion_positive_size_7_supports"] == 0
            and search_by_graph["3x5_open_grid"]["criterion_positive_size_7_supports"] == 0,
            "detail": "automorphism-only search over every connected support; no 4^|S| kernels used",
        },
        {
            "name": "induced_C4_explanation_rejected",
            "passed": len(c4_rejection["not_necessary_counterexamples"]) == 3
            and scalar_c4_row["kernel_dimension_Q"] == 1
            and all(
                item["kernel_span_intersection_dimension"] == 1
                for item in c4_rejection["only_exception_with_C4"]["C4_intersections"]
            ),
            "detail": "C4 is neither necessary nor sufficient, and the mask-254 kernel has only scalar intersections with its two C4 operator spaces",
        },
    ]
    artifact = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "python": platform.python_version(),
            "source_flux_support": "results/integrability/flux_support.json",
            "source_flux_support_sha256": hashlib.sha256(flux_bytes).hexdigest(),
            "method": "exact Fraction sparse Pauli algebra; graph automorphism pointwise stabilizers; two-prime integer commutator ranks for the finite converse budget",
        },
        "data": {
            "claim_tags": ["[THEOREM]", "[LEMMA]", "[COMPUTATION]", "[CONJECTURE]", "[UNRESOLVED]"],
            "pointwise_stabilizer_theorem": {
                "claim_tag": "[THEOREM]",
                "statement": "If sigma is a graph automorphism fixing V\\S pointwise, its qubit-site permutation U_sigma is supported in S and commutes with A=sum_v X_v and B=sum_{uv in E} Z_u Z_v.",
                "boundary_neighborhood_criterion": "[LEMMA] A permutation pi of S extends by identity to a graph automorphism iff pi preserves G[S] and N(s) intersect (V\\S) = N(pi(s)) intersect (V\\S) for every s in S.",
                "mechanism_iff": "A support admits a nonidentity site-permutation charge iff the pointwise stabilizer Aut(G)_(V\\S) is nontrivial.",
            },
            "alternating_bell_chain_theorem": {
                "claim_tag": "[THEOREM]",
                "statement": "For an involutive graph automorphism with two-cycle cells C_i, assign eta_i in {+1,-1}; force eta_i=-1 when C_i touches a fixed vertex and eta_i eta_j=-1 when an edge orbit joins C_i,C_j. If consistent, the tensor product of Phi-minus Bell states on eta=+1 cells and singlets on eta=-1 cells defines a projector P commuting with A and B. If a fixed f has no fixed-fixed neighbor and touches only eta=-1 cells, X_f P also commutes.",
                "proof_identity": "(X_u+X_v)|beta_eta>=0 and Z_v|beta_eta>=eta Z_u|beta_eta>; fixed-cell bond sums vanish for eta=-1 and paired intercell bond orbits vanish for opposite eta.",
            },
            "exceptional_3x3_supports": exceptional,
            "stored_grid_criterion": stored_criterion,
            "finite_graph_converse_scan": graph_scan,
            "grid_mechanism_searches": grid_searches,
            "induced_C4_rejection": c4_rejection,
            "scope": {
                "claim_tag": "[UNRESOLVED]",
                "statement": "The pointwise-stabilizer and Bell constructions are general theorems, but necessity of the stabilizer criterion beyond the explicit n<=6 graph budget and the stored rectangular supports is conjectural. No thermodynamic integrability or non-integrability conclusion is claimed.",
            },
        },
        "checks": checks,
    }
    if not all(check["passed"] for check in checks):
        raise AssertionError([check for check in checks if not check["passed"]])
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"wrote {RESULT.relative_to(ROOT)}")
    print("PASS")


if __name__ == "__main__":
    main()
