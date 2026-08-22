"""Exact arbitrary-order Jordan--Wigner audit and non-Hamiltonian counterexamples.

This producer contains only finite exact integer/GF(2) computations.  The
all-size arguments audited here are written in proofs/trichotomy_extend.md.
"""

from __future__ import annotations

import hashlib
import itertools
import platform
import resource
import time
from collections import Counter, deque
from math import factorial
from typing import Iterable, Sequence

CPU_BUDGET_SECONDS = 120.0
RSS_CAP_BYTES = 2_000_000_000


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(started: float, stage: str) -> None:
    used = time.process_time() - started
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(
            f"process-time budget exceeded at {stage}: "
            f"{used:.6f}>{CPU_BUDGET_SECONDS:.6f}"
        )
    rss = max_rss_bytes()
    if rss >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {rss}>={RSS_CAP_BYTES}")


def normalize_edges(
    n: int, edges: Iterable[tuple[int, int]]
) -> tuple[tuple[int, int], ...]:
    answer: set[tuple[int, int]] = set()
    for left, right in edges:
        if not (0 <= left < n and 0 <= right < n) or left == right:
            raise ValueError(f"invalid simple edge {(left, right)} on {n} vertices")
        answer.add((left, right) if left < right else (right, left))
    return tuple(sorted(answer))


def adjacency_masks(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    answer = [0] * n
    for left, right in normalize_edges(n, edges):
        answer[left] |= 1 << right
        answer[right] |= 1 << left
    return tuple(answer)


def bipartition(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[bool, bool, tuple[int, ...]]:
    """Return connected, bipartite, and a 0/1 colouring."""
    adjacency = adjacency_masks(n, edges)
    colours = [-1] * n
    connected = True
    bipartite = True
    components = 0
    for root in range(n):
        if colours[root] != -1:
            continue
        components += 1
        colours[root] = 0
        queue = deque([root])
        while queue:
            vertex = queue.popleft()
            neighbours = adjacency[vertex]
            while neighbours:
                bit = neighbours & -neighbours
                neighbours -= bit
                other = bit.bit_length() - 1
                if colours[other] == -1:
                    colours[other] = 1 - colours[vertex]
                    queue.append(other)
                elif colours[other] == colours[vertex]:
                    bipartite = False
    connected = components == 1
    return connected, bipartite, tuple(colours)


def maximum_degree(n: int, edges: Sequence[tuple[int, int]]) -> int:
    return max((mask.bit_count() for mask in adjacency_masks(n, edges)), default=0)


def edge_grades_for_order(
    n: int, edges: Sequence[tuple[int, int]], order: Sequence[int]
) -> tuple[int, ...]:
    if len(order) != n or set(order) != set(range(n)):
        raise ValueError("order must be a permutation of all vertices")
    position = {vertex: index for index, vertex in enumerate(order)}
    return tuple(
        2 * abs(position[left] - position[right])
        for left, right in normalize_edges(n, edges)
    )


def grade_compatible_order(
    n: int, edges: Sequence[tuple[int, int]], order: Sequence[int]
) -> bool:
    return all(grade % 4 == 2 for grade in edge_grades_for_order(n, edges, order))


def compatible_order_count_formula(part_sizes: tuple[int, int]) -> int:
    left, right = sorted(part_sizes)
    if right - left > 1:
        return 0
    if left == right:
        return 2 * factorial(left) * factorial(right)
    return factorial(left) * factorial(right)


def hamiltonian_path_exists(n: int, edges: Sequence[tuple[int, int]]) -> bool:
    """Exact subset DP; used only at n<=6."""
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
    return bool(endpoints[-1])


def raw_local_generators(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    fields = [1 << vertex for vertex in range(n)]
    bonds = [
        (1 << (n + left)) | (1 << (n + right))
        for left, right in normalize_edges(n, edges)
    ]
    return tuple(fields + bonds)


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
    """Exact generator-adjoint closure of monomial Pauli labels."""
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


def closure_digest(values: Iterable[int]) -> str:
    digest = hashlib.sha256()
    for value in sorted(values):
        digest.update(value.to_bytes(16, "little"))
    return digest.hexdigest()


def path_majorana_masks(n: int, order: Sequence[int]) -> tuple[int, ...]:
    """GF(2) Pauli labels of the 2n Jordan--Wigner Majoranas."""
    if len(order) != n or set(order) != set(range(n)):
        raise ValueError("order must be a permutation")
    answer: list[int] = []
    prefix_x = 0
    for vertex in order:
        answer.append(prefix_x | (1 << (n + vertex)))
        answer.append(prefix_x | (1 << vertex) | (1 << (n + vertex)))
        prefix_x ^= 1 << vertex
    return tuple(answer)


def inverse_majorana_basis(
    n: int, order: Sequence[int]
) -> dict[int, tuple[int, int]]:
    rows: dict[int, tuple[int, int]] = {}
    for index, pauli in enumerate(path_majorana_masks(n, order)):
        coefficients = 1 << index
        while pauli:
            pivot = pauli.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = pauli, coefficients
                break
            row, row_coefficients = rows[pivot]
            pauli ^= row
            coefficients ^= row_coefficients
        else:
            raise AssertionError("Jordan--Wigner Majoranas are not independent")
    return rows


def majorana_grade(pauli: int, inverse: dict[int, tuple[int, int]]) -> int:
    coefficients = 0
    while pauli:
        pivot = pauli.bit_length() - 1
        row, row_coefficients = inverse[pivot]
        pauli ^= row
        coefficients ^= row_coefficients
    return coefficients.bit_count()


def grade_histogram(
    values: Iterable[int], n: int, order: Sequence[int]
) -> dict[int, int]:
    inverse = inverse_majorana_basis(n, order)
    return dict(
        sorted(Counter(majorana_grade(value, inverse) for value in values).items())
    )


def exhaustive_order_criterion(max_n: int = 5) -> dict[str, int]:
    """Audit the balanced-bipartition criterion on every connected graph."""
    graph_count = 0
    order_count = 0
    for n in range(2, max_n + 1):
        pairs = tuple(itertools.combinations(range(n), 2))
        permutations = tuple(itertools.permutations(range(n)))
        for edge_mask in range(1 << len(pairs)):
            edges = tuple(
                pairs[index]
                for index in range(len(pairs))
                if (edge_mask >> index) & 1
            )
            connected, is_bipartite, colours = bipartition(n, edges)
            if not connected or not is_bipartite:
                continue
            graph_count += 1
            part_sizes = (colours.count(0), colours.count(1))
            measured = sum(
                grade_compatible_order(n, edges, order) for order in permutations
            )
            expected = compatible_order_count_formula(part_sizes)
            if measured != expected:
                raise AssertionError(
                    f"order criterion mismatch n={n}, edges={edges}: "
                    f"{measured}!={expected}"
                )
            order_count += measured
    return {"graphs": graph_count, "compatible_orders": order_count}


def graph_counterexample_row(
    label: str, n: int, edges: Sequence[tuple[int, int]], predicted: int
) -> dict[str, object]:
    connected, is_bipartite, colours = bipartition(n, edges)
    closure = exact_pauli_closure(n, edges)
    all_orders = tuple(itertools.permutations(range(n)))
    compatible = tuple(
        order for order in all_orders if grade_compatible_order(n, edges, order)
    )
    representative = tuple(range(n))
    return {
        "tag": "[COMPUTATION]",
        "label": label,
        "n": n,
        "edges": [list(edge) for edge in normalize_edges(n, edges)],
        "connected": connected,
        "bipartite": is_bipartite,
        "bipartition_sizes": sorted((colours.count(0), colours.count(1))),
        "maximum_degree": maximum_degree(n, edges),
        "hamiltonian_path_exists": hamiltonian_path_exists(n, edges),
        "compatible_order_count": len(compatible),
        "compatible_order_count_formula": compatible_order_count_formula(
            (colours.count(0), colours.count(1))
        ),
        "representative_order": list(representative),
        "representative_generator_grades": list(
            edge_grades_for_order(n, edges, representative)
        ),
        "representative_closure_grade_histogram": {
            str(grade): count
            for grade, count in grade_histogram(closure, n, representative).items()
        },
        "dimension": len(closure),
        "hamiltonian_branching_formula_if_misapplied": predicted,
        "closure_sha256": closure_digest(closure),
        "pairwise_closed": pairwise_closed(closure, n),
    }


def run_order_audit() -> dict[str, object]:
    started = time.process_time()
    finite_criterion = exhaustive_order_criterion()
    budget_tick(started, "exhaustive order criterion")

    claw = graph_counterexample_row(
        "bipartite_claw_K1_3",
        4,
        ((0, 1), (0, 2), (0, 3)),
        56,
    )
    branch_tree = graph_counterexample_row(
        "bipartite_branch_tree_T6",
        6,
        ((0, 1), (1, 2), (1, 3), (3, 4), (4, 5)),
        1056,
    )
    budget_tick(started, "counterexample closures")

    checks = [
        {
            "name": "finite order criterion",
            "passed": finite_criterion["graphs"] > 0,
            "detail": str(finite_criterion),
        },
        {
            "name": "K1,3 is exact counterexample",
            "passed": claw["dimension"] == 72
            and claw["hamiltonian_branching_formula_if_misapplied"] == 56
            and claw["compatible_order_count"] == 0
            and claw["pairwise_closed"] is True,
            "detail": "72 != 56; no grade-compatible order",
        },
        {
            "name": "second non-Hamiltonian mismatch",
            "passed": branch_tree["dimension"] == 992
            and branch_tree["hamiltonian_branching_formula_if_misapplied"] == 1056
            and branch_tree["hamiltonian_path_exists"] is False,
            "detail": "T6: 992 != 1056",
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
        "order_criterion_finite_audit": finite_criterion,
        "counterexamples": [claw, branch_tree],
        "checks": checks,
        "process_time_seconds": time.process_time() - started,
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_order_audit()
    print(
        "PASS "
        f"(K1,3={result['counterexamples'][0]['dimension']}, "
        f"T6={result['counterexamples'][1]['dimension']}, "
        f"process_time={result['process_time_seconds']:.6f}s, "
        f"peak_rss={result['peak_rss_bytes']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
