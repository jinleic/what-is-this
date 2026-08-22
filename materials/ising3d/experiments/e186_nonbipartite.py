"""Exact support and raw-generator audit for the non-bipartite branch."""

from __future__ import annotations

import time
from collections import Counter
from math import comb
from typing import Iterable, Sequence

from e185_order_counterexample import (
    CPU_BUDGET_SECONDS,
    RSS_CAP_BYTES,
    bipartition,
    budget_tick,
    closure_digest,
    exact_pauli_closure,
    grade_histogram,
    hamiltonian_path_exists,
    max_rss_bytes,
    maximum_degree,
    normalize_edges,
    pairwise_closed,
)


def cycle_edges(n: int) -> tuple[tuple[int, int], ...]:
    if n < 3:
        raise ValueError("a simple cycle needs at least three vertices")
    return normalize_edges(
        n,
        tuple((vertex, vertex + 1) for vertex in range(n - 1)) + ((n - 1, 0),),
    )


def validate_hamiltonian_path(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> None:
    if len(path) != n or set(path) != set(range(n)):
        raise ValueError("path must visit every vertex once")
    edge_set = set(normalize_edges(n, edges))
    for left, right in zip(path, path[1:]):
        edge = (left, right) if left < right else (right, left)
        if edge not in edge_set:
            raise ValueError(f"missing path edge {edge}")


def chord_grades(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> tuple[int, ...]:
    validate_hamiltonian_path(n, edges, path)
    path_edges = {
        (left, right) if left < right else (right, left)
        for left, right in zip(path, path[1:])
    }
    position = {vertex: index for index, vertex in enumerate(path)}
    return tuple(
        2 * abs(position[left] - position[right])
        for left, right in normalize_edges(n, edges)
        if (left, right) not in path_edges
    )


def bracket_grade_outputs(N: int, left: int, right: int) -> tuple[int, ...]:
    answer: set[int] = set()
    for overlap in range(1, min(left, right) + 1, 2):
        if left + right - overlap <= N:
            grade = left + right - 2 * overlap
            if 0 < grade < N:
                answer.add(grade)
    return tuple(sorted(answer))


def grade_closure(N: int, seeds: Iterable[int]) -> tuple[int, ...]:
    grades = {int(grade) for grade in seeds}
    changed = True
    while changed:
        changed = False
        present = tuple(grades)
        for left in present:
            for right in present:
                for grade in bracket_grade_outputs(N, left, right):
                    if grade not in grades:
                        grades.add(grade)
                        changed = True
    return tuple(sorted(grades))


def all_noncentral_even_grades(N: int) -> tuple[int, ...]:
    return tuple(range(2, N, 2))


def full_even_dimension(n: int) -> int:
    return (1 << (2 * n - 1)) - 2


def full_even_histogram(n: int) -> dict[int, int]:
    return {grade: comb(2 * n, grade) for grade in range(2, 2 * n, 2)}


def endpoint_cycle_histogram(n: int) -> dict[int, int]:
    N = 2 * n
    return {2: comb(N, 2), N - 2: comb(N, N - 2)}


def support_generation_audit(max_N: int = 80) -> dict[str, object]:
    interior_rows: list[dict[str, int]] = []
    endpoint_rows: list[dict[str, int]] = []
    for N in range(8, max_N + 1, 2):
        target = all_noncentral_even_grades(N)
        for seed in range(4, N - 3, 4):
            reached = grade_closure(N, (2, seed))
            if reached != target:
                raise AssertionError(
                    f"interior 0-mod-4 seed failed N={N}, seed={seed}: {reached}"
                )
            interior_rows.append({"N": N, "seed": seed, "grade_count": len(reached)})
        if N % 4 == 2:
            for seed in range(6, N - 3, 4):
                reached = grade_closure(N, (2, seed, N - 2))
                if reached != target:
                    raise AssertionError(
                        f"endpoint-plus-interior seed failed N={N}, seed={seed}: {reached}"
                    )
                endpoint_rows.append(
                    {"N": N, "seed": seed, "grade_count": len(reached)}
                )
    return {
        "max_N": max_N,
        "interior_zero_mod_four_cases": len(interior_rows),
        "endpoint_plus_two_mod_four_cases": len(endpoint_rows),
        "interior_samples": interior_rows[:4] + interior_rows[-4:],
        "endpoint_samples": endpoint_rows[:4] + endpoint_rows[-4:],
    }


def exact_graph_row(
    *,
    label: str,
    n: int,
    edges: Sequence[tuple[int, int]],
    expected_branch: str,
    order: Sequence[int] | None,
) -> dict[str, object]:
    normalized = normalize_edges(n, edges)
    connected, is_bipartite, _ = bipartition(n, normalized)
    closure = exact_pauli_closure(n, normalized)
    representative_order = tuple(range(n)) if order is None else tuple(order)
    histogram = grade_histogram(closure, n, representative_order)
    expected_histogram = (
        endpoint_cycle_histogram(n)
        if expected_branch == "odd_cycle"
        else full_even_histogram(n)
    )
    expected_dimension = (
        2 * comb(2 * n, 2)
        if expected_branch == "odd_cycle"
        else full_even_dimension(n)
    )
    has_hamiltonian_path = hamiltonian_path_exists(n, normalized)
    chosen_chord_grades: tuple[int, ...] | None = None
    if order is not None:
        chosen_chord_grades = chord_grades(n, normalized, order)
    if not connected or is_bipartite:
        raise AssertionError(f"{label} is not connected non-bipartite")
    if len(closure) != expected_dimension or histogram != expected_histogram:
        raise AssertionError(
            f"{label}: dimension/histogram mismatch: {len(closure)}, {histogram}"
        )
    if not pairwise_closed(closure, n):
        raise AssertionError(f"{label}: generator closure is not pairwise saturated")
    return {
        "tag": "[COMPUTATION]",
        "label": label,
        "n": n,
        "edges": [list(edge) for edge in normalized],
        "connected": connected,
        "bipartite": is_bipartite,
        "maximum_degree": maximum_degree(n, normalized),
        "hamiltonian_path_exists": has_hamiltonian_path,
        "chosen_hamiltonian_order": list(order) if order is not None else None,
        "chord_grades": list(chosen_chord_grades)
        if chosen_chord_grades is not None
        else None,
        "branch": expected_branch,
        "dimension": len(closure),
        "expected_dimension": expected_dimension,
        "grade_histogram": {str(grade): count for grade, count in histogram.items()},
        "closure_sha256": closure_digest(closure),
        "pairwise_closed": True,
    }


def run_nonbipartite_audit() -> dict[str, object]:
    started = time.process_time()
    support_audit = support_generation_audit()
    budget_tick(started, "support generation")

    odd_cycle = exact_graph_row(
        label="odd_cycle_C5_exception",
        n=5,
        edges=cycle_edges(5),
        expected_branch="odd_cycle",
        order=tuple(range(5)),
    )
    paw = exact_graph_row(
        label="paw_C3_plus_leaf",
        n=4,
        edges=((0, 1), (1, 2), (2, 0), (0, 3)),
        expected_branch="full_even_noncentral",
        order=(3, 0, 1, 2),
    )
    c5_chord = exact_graph_row(
        label="C5_plus_chord_0_3",
        n=5,
        edges=(*cycle_edges(5), (0, 3)),
        expected_branch="full_even_noncentral",
        order=tuple(range(5)),
    )
    triangle_two_leaves = exact_graph_row(
        label="nonhamiltonian_triangle_two_leaves",
        n=5,
        edges=((0, 1), (1, 2), (2, 0), (0, 3), (0, 4)),
        expected_branch="full_even_noncentral",
        order=None,
    )
    c6_short_chord = exact_graph_row(
        label="C6_plus_even_distance_chord_0_2",
        n=6,
        edges=(*cycle_edges(6), (0, 2)),
        expected_branch="full_even_noncentral",
        order=tuple(range(6)),
    )
    graph_rows = [
        odd_cycle,
        paw,
        c5_chord,
        triangle_two_leaves,
        c6_short_chord,
    ]
    budget_tick(started, "raw graph closures")

    measured = {row["label"]: row["dimension"] for row in graph_rows}
    checks = [
        {
            "name": "all support-generation cases",
            "passed": support_audit["interior_zero_mod_four_cases"] > 0
            and support_audit["endpoint_plus_two_mod_four_cases"] > 0,
            "detail": str(support_audit),
        },
        {
            "name": "odd-cycle exception",
            "passed": measured["odd_cycle_C5_exception"] == 90,
            "detail": "C5 has grades 2 and 8 only",
        },
        {
            "name": "three parity structures",
            "passed": measured["paw_C3_plus_leaf"] == 126
            and measured["C5_plus_chord_0_3"] == 510
            and measured["C6_plus_even_distance_chord_0_2"] == 2046,
            "detail": "126, 510, 2046",
        },
        {
            "name": "non-Hamiltonian leaf extension control",
            "passed": measured["nonhamiltonian_triangle_two_leaves"] == 510
            and triangle_two_leaves["hamiltonian_path_exists"] is False,
            "detail": "triangle with two leaves: full dimension 510",
        },
        {
            "name": "all raw closures pairwise saturated",
            "passed": all(bool(row["pairwise_closed"]) for row in graph_rows),
            "detail": f"{len(graph_rows)} closures",
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
        "support_generation_audit": support_audit,
        "graph_rows": graph_rows,
        "checks": checks,
        "process_time_seconds": time.process_time() - started,
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_nonbipartite_audit()
    dimensions = {
        row["label"]: row["dimension"] for row in result["graph_rows"]
    }
    print(
        f"PASS (dimensions={dimensions}, "
        f"process_time={result['process_time_seconds']:.6f}s, "
        f"peak_rss={result['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
