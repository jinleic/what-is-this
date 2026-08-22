#!/usr/bin/env python3
"""Independent exact checks for the support-seven charge criterion artifact."""
from __future__ import annotations

import itertools
import json
import signal
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "charge_criterion.json"
FLUX = ROOT / "results" / "integrability" / "flux_support.json"
PRIME = 2_147_483_647
TIMEOUT_SECONDS = 300


def pc(mask: int) -> int:
    return bin(mask).count("1")


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        [(r * columns + c, (r + 1) * columns + c) for r in range(rows - 1) for c in range(columns)]
        + [(r * columns + c, r * columns + c + 1) for r in range(rows) for c in range(columns - 1)]
    )


def grid_autos(rows: int, columns: int) -> tuple[tuple[int, ...], ...]:
    maps = [
        lambda r, c: (r, c),
        lambda r, c: (rows - 1 - r, c),
        lambda r, c: (r, columns - 1 - c),
        lambda r, c: (rows - 1 - r, columns - 1 - c),
    ]
    if rows == columns:
        maps += [
            lambda r, c: (c, r),
            lambda r, c: (columns - 1 - c, rows - 1 - r),
            lambda r, c: (c, columns - 1 - r),
            lambda r, c: (rows - 1 - c, r),
        ]
    result = []
    for transform in maps:
        permutation = tuple(
            transform(v // columns, v % columns)[0] * columns
            + transform(v // columns, v % columns)[1]
            for v in range(rows * columns)
        )
        if permutation not in result:
            result.append(permutation)
    return tuple(result)


def transform(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[v] for v in range(len(permutation)) if (mask >> v) & 1)


def stabilizer(automorphisms: tuple[tuple[int, ...], ...], n: int, support: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        permutation
        for permutation in automorphisms
        if all(permutation[v] == v for v in range(n) if not ((support >> v) & 1))
    )


def connected(n: int, edges: tuple[tuple[int, int], ...], support: int) -> bool:
    if not support:
        return False
    adjacency = [0] * n
    for u, v in edges:
        adjacency[u] |= 1 << v
        adjacency[v] |= 1 << u
    reached = support & -support
    frontier = reached
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        vertex = bit.bit_length() - 1
        new = adjacency[vertex] & support & ~reached
        reached |= new
        frontier |= new
    return reached == support


def clean(operator: dict[tuple[int, int], Fraction]) -> dict[tuple[int, int], Fraction]:
    return {key: Fraction(value) for key, value in operator.items() if value}


def add(*terms: tuple[int | Fraction, dict[tuple[int, int], Fraction]]) -> dict[tuple[int, int], Fraction]:
    output: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for scale, operator in terms:
        for key, value in operator.items():
            output[key] += Fraction(scale) * value
    return clean(dict(output))


def multiply(
    left: dict[tuple[int, int], Fraction], right: dict[tuple[int, int], Fraction]
) -> dict[tuple[int, int], Fraction]:
    output: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            output[(a ^ c, b ^ d)] += (-1 if pc(b & c) % 2 else 1) * x * y
    return clean(dict(output))


def pair(u: int, v: int, kind: str) -> dict[tuple[int, int], Fraction]:
    mask = (1 << u) | (1 << v)
    if kind == "swap":
        values, denominator = (1, 1, 1, -1), 2
    elif kind == "singlet":
        values, denominator = (1, -1, -1, 1), 4
    elif kind == "phi":
        values, denominator = (1, -1, 1, -1), 4
    else:
        raise ValueError(kind)
    return {
        key: Fraction(value, denominator)
        for key, value in zip(((0, 0), (mask, 0), (0, mask), (mask, mask)), values)
    }


def product(operators: list[dict[tuple[int, int], Fraction]]) -> dict[tuple[int, int], Fraction]:
    result = {(0, 0): Fraction(1)}
    for operator in operators:
        result = multiply(result, operator)
    return result


def stored(vector: dict) -> dict[tuple[int, int], Fraction]:
    return {
        (a, b): Fraction(numerator, denominator)
        for a, b, numerator, denominator in vector["terms_a_b_num_den"]
    }


def commutators(
    operator: dict[tuple[int, int], Fraction], n: int, edges: tuple[tuple[int, int], ...]
) -> tuple[dict[tuple[int, int], Fraction], dict[tuple[int, int], Fraction]]:
    ad_a: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    ad_b: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (a, b), coefficient in operator.items():
        for vertex in range(n):
            if (b >> vertex) & 1:
                ad_a[(a ^ (1 << vertex), b)] += coefficient
        for u, v in edges:
            if ((a >> u) ^ (a >> v)) & 1:
                ad_b[(a, b ^ (1 << u) ^ (1 << v))] -= coefficient
    return clean(dict(ad_a)), clean(dict(ad_b))


def check_projectors(projectors: list[dict[tuple[int, int], Fraction]], ranks: list[int]) -> None:
    identity = {(0, 0): Fraction(1)}
    assert add(*[(1, projector) for projector in projectors]) == identity
    assert all(multiply(projector, projector) == projector for projector in projectors)
    assert all(
        not multiply(projectors[i], projectors[j])
        for i in range(len(projectors))
        for j in range(len(projectors))
        if i != j
    )
    assert [128 * projector.get((0, 0), 0) for projector in projectors] == ranks


def exact_exception_forms(flux: dict, result: dict) -> None:
    layer = next(layer for layer in flux["data"]["layers"] if layer["graph"] == "3x3_open_grid")
    rows = {row["support_mask"]: row for row in layer["table"]}
    identity = {(0, 0): Fraction(1)}
    cases = {
        239: ([(1, 3), (2, 6), (5, 7)], 0),
        254: ([(1, 3), (2, 6), (5, 7)], 4),
        367: ([(0, 2), (3, 5), (6, 8)], None),
        381: ([(0, 2), (3, 5), (6, 8)], None),
    }
    summaries = {entry["support_mask"]: entry for entry in result["data"]["exceptional_3x3_supports"]}
    for mask, (pairs, fixed_vertex) in cases.items():
        swap = product([pair(u, v, "swap") for u, v in pairs])
        if fixed_vertex is not None:
            bell = product([pair(*pairs[0], "singlet"), pair(*pairs[1], "phi"), pair(*pairs[2], "singlet")])
            x_bell = multiply({(1 << fixed_vertex, 0): Fraction(1)}, bell)
            conceptual = [
                identity,
                add((4, swap), (-32, bell)),
                add((1, identity), (-4, swap), (-32, bell)),
                add((-64, x_bell)),
            ]
            projectors = [
                add((Fraction(1, 2), identity), (Fraction(-1, 2), swap)),
                add((Fraction(1, 2), identity), (Fraction(1, 2), swap), (-1, bell)),
                add((Fraction(1, 2), bell), (Fraction(1, 2), x_bell)),
                add((Fraction(1, 2), bell), (Fraction(-1, 2), x_bell)),
            ]
            check_projectors(projectors, [56, 70, 1, 1])
            assert multiply(swap, bell) == bell
        else:
            conceptual = [identity, add((1, identity), (-8, swap))]
            projectors = [
                add((Fraction(1, 2), identity), (Fraction(1, 2), swap)),
                add((Fraction(1, 2), identity), (Fraction(-1, 2), swap)),
            ]
            check_projectors(projectors, [72, 56])
        exact_basis = [stored(vector) for vector in rows[mask]["exact_basis"]]
        assert conceptual == exact_basis
        assert all(commutators(operator, 9, grid_edges(3, 3)) == ({}, {}) for operator in conceptual + projectors)
        assert summaries[mask]["conceptual_basis_exactly_equals_stored_basis"]
        assert summaries[mask]["global_commutant_projector_factorization"]["orthogonal_idempotent_partition_exact"]


def verify_all_stored_rows(flux: dict, result: dict) -> None:
    checked = 0
    positives = []
    for layer in flux["data"]["layers"]:
        rows, columns = layer["shape"]
        autos = grid_autos(rows, columns)
        n = rows * columns
        for row in layer["table"]:
            checked += 1
            order = len(stabilizer(autos, n, row["support_mask"]))
            assert (order > 1) == (row["kernel_dimension_Q"] > 1)
            if order > 1:
                positives.append((layer["graph"], row["support_mask"], row["kernel_dimension_Q"], order))
    assert checked == 250
    assert positives == [
        ("3x3_open_grid", 239, 4, 2),
        ("3x3_open_grid", 367, 2, 2),
        ("3x3_open_grid", 381, 2, 2),
        ("3x3_open_grid", 254, 4, 2),
    ]
    stored_summary = result["data"]["stored_grid_criterion"]
    assert stored_summary["rows_checked"] == checked and not stored_summary["mismatches"]


def verify_grid_searches(result: dict) -> None:
    stored_searches = {entry["graph"]: entry for entry in result["data"]["grid_mechanism_searches"]}
    expected = {
        (2, 4): (8, 0),
        (3, 3): (7, 12),
        (3, 4): (9, 0),
        (4, 4): (13, 0),
        (3, 5): (11, 0),
    }
    for (rows, columns), (minimum, size_seven_count) in expected.items():
        n = rows * columns
        edges = grid_edges(rows, columns)
        autos = grid_autos(rows, columns)
        positive = Counter()
        connected_seven = 0
        for support in range(1, 1 << n):
            if not connected(n, edges, support):
                continue
            if pc(support) == 7:
                connected_seven += 1
            if len(stabilizer(autos, n, support)) > 1:
                positive[pc(support)] += 1
        graph = f"{rows}x{columns}_open_grid"
        assert min(positive) == minimum
        assert positive[7] == size_seven_count
        assert stored_searches[graph]["minimum_connected_criterion_support"] == minimum
        assert stored_searches[graph]["criterion_positive_size_7_supports"] == size_seven_count
        assert stored_searches[graph]["connected_size_7_supports"] == connected_seven


def complete_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((u, v) for u in range(n) for v in range(u + 1, n))


@lru_cache(maxsize=None)
def permutation_maps(n: int) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    edges = complete_edges(n)
    positions = {edge: index for index, edge in enumerate(edges)}
    return tuple(
        (
            permutation,
            tuple(positions[tuple(sorted((permutation[u], permutation[v])))] for u, v in edges),
        )
        for permutation in itertools.permutations(range(n))
    )


def map_graph(mask: int, edge_map: tuple[int, ...]) -> int:
    output = 0
    for index in range(len(edge_map)):
        if (mask >> index) & 1:
            output |= 1 << edge_map[index]
    return output


def graph_edges(n: int, mask: int) -> tuple[tuple[int, int], ...]:
    return tuple(edge for index, edge in enumerate(complete_edges(n)) if (mask >> index) & 1)


def graph_mask(n: int, edges: set[tuple[int, int]]) -> int:
    positions = {edge: index for index, edge in enumerate(complete_edges(n))}
    return sum(1 << positions[tuple(sorted(edge))] for edge in edges)


def canonical(n: int, mask: int) -> int:
    return min(map_graph(mask, edge_map) for _, edge_map in permutation_maps(n))


def extend(previous: list[int], n: int) -> list[int]:
    output = set()
    for old in previous:
        base = set(graph_edges(n - 1, old))
        for neighborhood in range(1, 1 << (n - 1)):
            edges = base | {(vertex, n - 1) for vertex in range(n - 1) if (neighborhood >> vertex) & 1}
            output.add(canonical(n, graph_mask(n, edges)))
    return sorted(output)


def graph_autos(n: int, mask: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutation for permutation, edge_map in permutation_maps(n) if map_graph(mask, edge_map) == mask)


def columns(n: int, edges: tuple[tuple[int, int], ...], support: int) -> list[dict[tuple[int, int], int]]:
    sites = [vertex for vertex in range(n) if (support >> vertex) & 1]
    incident = [edge for edge in edges if ((support >> edge[0]) & 1) or ((support >> edge[1]) & 1)]
    output = []
    for index in range(1 << (2 * len(sites))):
        digits = index
        a = b = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a |= 1 << site
            if label & 2:
                b |= 1 << site
        pauli = a | (b << n)
        column = {}
        for site in sites:
            if (b >> site) & 1:
                column[(0, pauli ^ (1 << site))] = 1
        for u, v in incident:
            if ((a >> u) ^ (a >> v)) & 1:
                column[(1, pauli ^ (1 << (n + u)) ^ (1 << (n + v)))] = -1
        output.append(column)
    return output


def rank_mod(vectors: list[dict[tuple[int, int], int]]) -> int:
    frequency = Counter(row for vector in vectors for row in vector)
    order = {row: index for index, row in enumerate(sorted(frequency, key=lambda row: (frequency[row], row)))}
    pivots = {}
    rank = 0
    for source in sorted(vectors, key=len):
        vector = {row: value % PRIME for row, value in source.items()}
        while vector:
            lead = min(vector, key=order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], PRIME - 2, PRIME)
                pivots[lead] = {row: value * inverse % PRIME for row, value in vector.items()}
                rank += 1
                break
            factor = vector[lead]
            for row, value in old.items():
                reduced = (vector.get(row, 0) - factor * value) % PRIME
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return rank


def verify_small_graph_budget(result: dict) -> None:
    representatives = {1: [0]}
    stored_per_n = {entry["n"]: entry for entry in result["data"]["finite_graph_converse_scan"]["per_n"]}
    for n in range(2, 6):
        representatives[n] = extend(representatives[n - 1], n)
        cases = positive = trivial = 0
        nullities = Counter()
        sizes = Counter()
        for mask in representatives[n]:
            edges = graph_edges(n, mask)
            autos = graph_autos(n, mask)
            for support in range(1, (1 << n) - 1):
                if not connected(n, edges, support):
                    continue
                if support != min(transform(support, permutation) for permutation in autos):
                    continue
                pointwise = len(stabilizer(autos, n, support)) > 1
                matrix = columns(n, edges, support)
                nullity = len(matrix) - rank_mod(matrix)
                cases += 1
                sizes[pc(support)] += 1
                nullities[nullity] += 1
                positive += int(pointwise)
                trivial += int(not pointwise)
                assert pointwise == (nullity > 1)
        stored_n = stored_per_n[n]
        assert stored_n["connected_unlabeled_graphs"] == len(representatives[n])
        assert stored_n["proper_connected_support_orbits"] == cases
        assert stored_n["nontrivial_pointwise_stabilizer_cases"] == positive
        assert stored_n["trivial_pointwise_stabilizer_cases"] == trivial
        assert stored_n["modular_nullity_histogram"] == {str(key): nullities[key] for key in sorted(nullities)}
        assert stored_n["support_size_histogram"] == {str(key): sizes[key] for key in sorted(sizes)}


def main() -> None:
    signal.alarm(TIMEOUT_SECONDS)
    result = json.loads(RESULT.read_text())
    flux = json.loads(FLUX.read_text())
    assert set(result) == {"provenance", "data", "checks"}
    assert result["provenance"]["script"] == "experiments/e61_charge_criterion.py"
    assert result["provenance"]["interpreter"] == ".venv/bin/python"
    assert result["checks"] and all(check["passed"] for check in result["checks"])
    exact_exception_forms(flux, result)
    verify_all_stored_rows(flux, result)
    verify_grid_searches(result)
    verify_small_graph_budget(result)
    rejection = result["data"]["induced_C4_rejection"]
    assert {entry["support_mask"] for entry in rejection["not_necessary_counterexamples"]} == {239, 367, 381}
    assert rejection["not_sufficient_counterexample"]["support_mask"] == 27
    assert all(
        entry["kernel_span_intersection_dimension"] == 1
        for entry in rejection["only_exception_with_C4"]["C4_intersections"]
    )
    assert result["data"]["finite_graph_converse_scan"]["total_support_orbits"] == 2682
    assert not result["data"]["finite_graph_converse_scan"]["criterion_mismatches"]
    signal.alarm(0)
    print("PASS")


if __name__ == "__main__":
    main()
