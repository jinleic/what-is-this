"""Exact bounded census of connected bipartite local-term Lie algebras.

The census is exhaustive up to graph isomorphism through six vertices.  All
arithmetic is integer or over GF(2); no numerical linear algebra is used.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import platform
import resource
import time
from collections import deque
from functools import lru_cache
from typing import Iterable, Sequence

CPU_BUDGET_SECONDS = 90.0
RSS_CAP_BYTES = 2_000_000_000
CENSUS_MAX_N = 6


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    used = time.process_time() - started
    if used > CPU_BUDGET_SECONDS:
        raise TimeoutError(
            f"process-time budget exceeded at {stage}: "
            f"{used:.6f}>{CPU_BUDGET_SECONDS:.6f}"
        )
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {rss}>={RSS_CAP_BYTES}")


@lru_cache(maxsize=None)
def vertex_pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(n), 2))


@lru_cache(maxsize=None)
def pair_indices(n: int) -> dict[tuple[int, int], int]:
    return {pair: index for index, pair in enumerate(vertex_pairs(n))}


def normalize_edges(
    n: int, edges: Iterable[tuple[int, int]]
) -> tuple[tuple[int, int], ...]:
    answer: set[tuple[int, int]] = set()
    for left, right in edges:
        if not 0 <= left < n or not 0 <= right < n or left == right:
            raise ValueError(f"invalid edge {(left, right)} for n={n}")
        answer.add((left, right) if left < right else (right, left))
    return tuple(sorted(answer))


def edge_code(n: int, edges: Iterable[tuple[int, int]]) -> int:
    indices = pair_indices(n)
    answer = 0
    for edge in normalize_edges(n, edges):
        answer |= 1 << indices[edge]
    return answer


def edges_from_code(n: int, code: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        pair
        for index, pair in enumerate(vertex_pairs(n))
        if (code >> index) & 1
    )


@lru_cache(maxsize=None)
def relabel_bit_maps(n: int) -> tuple[tuple[int, ...], ...]:
    indices = pair_indices(n)
    maps: list[tuple[int, ...]] = []
    for permutation in itertools.permutations(range(n)):
        mapped: list[int] = []
        for left, right in vertex_pairs(n):
            new_left, new_right = sorted((permutation[left], permutation[right]))
            mapped.append(1 << indices[(new_left, new_right)])
        maps.append(tuple(mapped))
    return tuple(maps)


def canonical_code(n: int, code: int) -> int:
    best = code
    for bit_map in relabel_bit_maps(n):
        candidate = 0
        live = code
        while live:
            bit = live & -live
            live -= bit
            candidate |= bit_map[bit.bit_length() - 1]
        if candidate < best:
            best = candidate
    return best


def adjacency_masks(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    answer = [0] * n
    for left, right in normalize_edges(n, edges):
        answer[left] |= 1 << right
        answer[right] |= 1 << left
    return tuple(answer)


def bipartition(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[bool, bool, tuple[int, ...]]:
    if n < 1:
        return False, True, ()
    adjacency = adjacency_masks(n, edges)
    colours = [-1] * n
    colours[0] = 0
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        live = adjacency[vertex]
        while live:
            bit = live & -live
            live -= bit
            neighbour = bit.bit_length() - 1
            if colours[neighbour] < 0:
                colours[neighbour] = colours[vertex] ^ 1
                queue.append(neighbour)
            elif colours[neighbour] == colours[vertex]:
                return False, False, tuple(colours)
    connected = all(colour >= 0 for colour in colours)
    return connected, connected, tuple(colours)


def maximum_degree(n: int, edges: Sequence[tuple[int, int]]) -> int:
    return max((row.bit_count() for row in adjacency_masks(n, edges)), default=0)


def path_endpoint_table(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    adjacency = adjacency_masks(n, edges)
    endpoints = [0] * (1 << n)
    for vertex in range(n):
        endpoints[1 << vertex] = 1 << vertex
    for mask in range(1 << n):
        live = endpoints[mask]
        while live:
            end_bit = live & -live
            live -= end_bit
            end = end_bit.bit_length() - 1
            choices = adjacency[end] & ~mask
            while choices:
                next_bit = choices & -choices
                choices -= next_bit
                endpoints[mask | next_bit] |= next_bit
    return tuple(endpoints)


def hamiltonian_path_exists(n: int, edges: Sequence[tuple[int, int]]) -> bool:
    return bool(path_endpoint_table(n, edges)[-1])


def minimum_path_cover(n: int, edges: Sequence[tuple[int, int]]) -> int:
    """Minimum number of vertex-disjoint paths covering every vertex."""
    endpoints = path_endpoint_table(n, edges)
    pathable = tuple(mask == 0 or bool(endpoints[mask]) for mask in range(1 << n))
    best = [n + 1] * (1 << n)
    best[0] = 0
    for mask in range(1, 1 << n):
        anchor = mask & -mask
        subset = mask
        while subset:
            if subset & anchor and pathable[subset]:
                best[mask] = min(best[mask], 1 + best[mask ^ subset])
            subset = (subset - 1) & mask
    return best[-1]


def maximum_matching_size(n: int, edges: Sequence[tuple[int, int]]) -> int:
    adjacency = adjacency_masks(n, edges)

    @lru_cache(maxsize=None)
    def solve(remaining: int) -> int:
        if not remaining:
            return 0
        vertex_bit = remaining & -remaining
        vertex = vertex_bit.bit_length() - 1
        without_vertex = remaining ^ vertex_bit
        answer = solve(without_vertex)
        choices = adjacency[vertex] & without_vertex
        while choices:
            neighbour_bit = choices & -choices
            choices -= neighbour_bit
            answer = max(answer, 1 + solve(without_vertex ^ neighbour_bit))
        return answer

    return solve((1 << n) - 1)


def raw_local_generators(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    fields = tuple(1 << vertex for vertex in range(n))
    bonds = tuple(
        (1 << (n + left)) | (1 << (n + right))
        for left, right in normalize_edges(n, edges)
    )
    return fields + bonds


def symplectic_form(left: int, right: int, n: int) -> int:
    mask = (1 << n) - 1
    x_left, z_left = left & mask, left >> n
    x_right, z_right = right & mask, right >> n
    return (
        (x_left & z_right).bit_count() + (z_left & x_right).bit_count()
    ) & 1


def exact_pauli_closure(
    n: int, edges: Sequence[tuple[int, int]]
) -> frozenset[int]:
    """Exact generator-adjoint closure of packed Pauli labels."""
    generators = raw_local_generators(n, edges)
    seen = set(generators)
    queue = deque(generators)
    while queue:
        value = queue.popleft()
        for generator in generators:
            if symplectic_form(value, generator, n):
                child = value ^ generator
                if child not in seen:
                    seen.add(child)
                    queue.append(child)
    return frozenset(seen)


def pairwise_closed(values: frozenset[int], n: int) -> bool:
    ordered = tuple(values)
    return all(
        (left ^ right) in values
        for index, left in enumerate(ordered)
        for right in ordered[index + 1 :]
        if symplectic_form(left, right, n)
    )


def label_digest(values: Iterable[int], n: int) -> str:
    width = max(1, (2 * n + 7) // 8)
    digest = hashlib.sha256()
    for value in sorted(values):
        digest.update(value.to_bytes(width, "little"))
    return digest.hexdigest()


def quadratic_value(
    label: int, n: int, colours: Sequence[int]
) -> int:
    """Q_c(a,b)=a.b+|a|+c.b over GF(2)."""
    mask = (1 << n) - 1
    x_part, z_part = label & mask, label >> n
    colour_dot = sum(
        colour for vertex, colour in enumerate(colours) if (z_part >> vertex) & 1
    )
    return (
        (x_part & z_part).bit_count() + x_part.bit_count() + colour_dot
    ) & 1


def quadratic_nonradical_roots(
    n: int, colours: Sequence[int]
) -> frozenset[int]:
    global_x = (1 << n) - 1
    answer: set[int] = set()
    for z_part in range(1 << n):
        if z_part.bit_count() & 1:
            continue
        for x_part in range(1 << n):
            label = x_part | (z_part << n)
            if label != global_x and quadratic_value(label, n, colours):
                answer.add(label)
    return frozenset(answer)


def quadratic_dimension(n: int, part_sizes: Sequence[int]) -> int:
    left, right = sorted(part_sizes)
    if left + right != n:
        raise ValueError((n, part_sizes))
    if n & 1:
        return (1 << (2 * n - 2)) - 1
    return (1 << (2 * n - 2)) - ((-1) ** left) * (1 << (n - 1))


def graph_branch(n: int, edges: Sequence[tuple[int, int]]) -> str:
    degree = maximum_degree(n, edges)
    if degree >= 3:
        return "branching_quadratic_root"
    if len(edges) == n - 1:
        return "path"
    if len(edges) == n:
        return "even_cycle"
    raise AssertionError(f"unexpected connected degree-two graph: n={n}, edges={edges}")


def branch_dimension(branch: str, n: int, part_sizes: Sequence[int]) -> int:
    if branch == "path":
        return n * (2 * n - 1)
    if branch == "even_cycle":
        return 2 * n * (2 * n - 1)
    if branch == "branching_quadratic_root":
        return quadratic_dimension(n, part_sizes)
    raise ValueError(branch)


def enumerate_representative_codes(
    max_n: int = CENSUS_MAX_N, started: float | None = None
) -> tuple[dict[int, int], dict[int, tuple[int, ...]]]:
    """Enumerate every connected bipartite graph, retaining canonical labels."""
    if max_n > CENSUS_MAX_N:
        raise ValueError(f"bounded census supports max_n<={CENSUS_MAX_N}")
    labeled_counts: dict[int, int] = {}
    representatives: dict[int, tuple[int, ...]] = {}
    for n in range(2, max_n + 1):
        pair_count = len(vertex_pairs(n))
        labeled = 0
        canonical: list[int] = []
        for code in range(1 << pair_count):
            if code.bit_count() < n - 1:
                continue
            edges = edges_from_code(n, code)
            connected, is_bipartite, _ = bipartition(n, edges)
            if not connected or not is_bipartite:
                continue
            labeled += 1
            if canonical_code(n, code) == code:
                canonical.append(code)
            if started is not None and code and code % 4096 == 0:
                budget_tick(started, f"graph enumeration n={n}, code={code}")
        labeled_counts[n] = labeled
        representatives[n] = tuple(canonical)
    return labeled_counts, representatives


def graph_row(n: int, code: int) -> dict[str, object]:
    edges = edges_from_code(n, code)
    connected, is_bipartite, colours = bipartition(n, edges)
    if not connected or not is_bipartite:
        raise AssertionError((n, code))
    part_sizes = tuple(sorted((colours.count(0), colours.count(1))))
    closure = exact_pauli_closure(n, edges)
    quadratic_roots = quadratic_nonradical_roots(n, colours)
    matching_size = maximum_matching_size(n, edges)
    path_cover = minimum_path_cover(n, edges)
    branch = graph_branch(n, edges)
    hamiltonian = hamiltonian_path_exists(n, edges)
    predicted = branch_dimension(branch, n, part_sizes)
    return {
        "tag": "[COMPUTATION]",
        "n": n,
        "canonical_code": code,
        "canonical_code_hex": hex(code),
        "edges": [list(edge) for edge in edges],
        "edge_count": len(edges),
        "connected": connected,
        "bipartite": is_bipartite,
        "bipartition_sizes": list(part_sizes),
        "colour_imbalance": part_sizes[1] - part_sizes[0],
        "smaller_colour_parity": part_sizes[0] & 1,
        "maximum_degree": maximum_degree(n, edges),
        "hamiltonian_path_exists": hamiltonian,
        "maximum_matching_size": matching_size,
        "matching_deficiency": n - 2 * matching_size,
        "minimum_path_cover": path_cover,
        "branch": branch,
        "dimension": len(closure),
        "predicted_dimension": predicted,
        "quadratic_nonradical_root_count": len(quadratic_roots),
        "closure_equals_quadratic_roots": closure == quadratic_roots,
        "pairwise_closed": pairwise_closed(closure, n),
        "raw_generator_sha256": label_digest(raw_local_generators(n, edges), n),
        "closure_sha256": label_digest(closure, n),
        "quadratic_root_sha256": label_digest(quadratic_roots, n),
    }


def run_census(max_n: int = CENSUS_MAX_N) -> dict[str, object]:
    started = time.process_time()
    labeled_counts, representative_codes = enumerate_representative_codes(
        max_n, started
    )
    rows: list[dict[str, object]] = []
    for n in range(2, max_n + 1):
        for code in representative_codes[n]:
            rows.append(graph_row(n, code))
        budget_tick(started, f"exact closures n={n}")

    branching = [row for row in rows if row["maximum_degree"] >= 3]
    nonhamiltonian = [
        row for row in rows if row["hamiltonian_path_exists"] is False
    ]
    checks = [
        {
            "name": "canonical census is nonempty and unique",
            "passed": bool(rows)
            and len({(row["n"], row["canonical_code"]) for row in rows})
            == len(rows)
            and all(
                canonical_code(int(row["n"]), int(row["canonical_code"]))
                == row["canonical_code"]
                for row in rows
            ),
            "detail": f"{len(rows)} unlabeled representatives",
        },
        {
            "name": "all raw closures are exact Lie-closed label sets",
            "passed": all(bool(row["pairwise_closed"]) for row in rows),
            "detail": f"{len(rows)} pairwise saturation checks",
        },
        {
            "name": "all branch dimensions match exact closures",
            "passed": all(
                row["dimension"] == row["predicted_dimension"] for row in rows
            ),
            "detail": f"{len(rows)} exact dimension comparisons",
        },
        {
            "name": "every branching closure is its quadratic root level",
            "passed": bool(branching)
            and all(
                row["closure_equals_quadratic_roots"]
                and row["closure_sha256"] == row["quadratic_root_sha256"]
                for row in branching
            ),
            "detail": f"{len(branching)} branching graphs",
        },
        {
            "name": "every non-Hamiltonian census graph is branching",
            "passed": bool(nonhamiltonian)
            and all(row["maximum_degree"] >= 3 for row in nonhamiltonian),
            "detail": f"{len(nonhamiltonian)} non-Hamiltonian graphs",
        },
        {
            "name": "resource cap",
            "passed": max_rss_bytes() < RSS_CAP_BYTES,
            "detail": f"peak_rss_bytes={max_rss_bytes()}",
        },
    ]
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)
    return {
        "tag": "[COMPUTATION]",
        "scope": {
            "minimum_n": 2,
            "maximum_n": max_n,
            "isomorphism_reduction": "minimum labelled edge bit-code over every vertex permutation",
            "closure_arithmetic": "exact GF(2) Pauli labels with integer counting",
        },
        "counts": {
            "connected_bipartite_labeled_by_n": {
                str(n): labeled_counts[n] for n in range(2, max_n + 1)
            },
            "connected_bipartite_unlabeled_by_n": {
                str(n): len(representative_codes[n])
                for n in range(2, max_n + 1)
            },
            "total_unlabeled": len(rows),
            "branching": len(branching),
            "nonhamiltonian": len(nonhamiltonian),
        },
        "graphs": rows,
        "checks": checks,
        "process_time_seconds": time.process_time() - started,
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_census()
    print(
        "PASS "
        f"(unlabeled={result['counts']['total_unlabeled']}, "
        f"nonhamiltonian={result['counts']['nonhamiltonian']}, "
        f"process_time={result['process_time_seconds']:.6f}s, "
        f"peak_rss={result['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
