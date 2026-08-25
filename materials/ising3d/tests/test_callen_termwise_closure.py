#!/usr/bin/env python3
"""Independent verifier for the Callen termwise-support closure theorem."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import resource
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "correlations" / "callen_termwise_closure.json"
PRODUCER = ROOT / "experiments" / "e243_callen_termwise_closure.py"
BASE_CALLEN_PRODUCER = ROOT / "experiments" / "e240_callen_identity_system.py"
RSS_LIMIT_BYTES = 2 * 1024**3
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def adjacency_from_edges(order: int, edges: set[tuple[int, int]]) -> tuple[tuple[int, ...], ...]:
    rows = [[] for _ in range(order)]
    for left, right in sorted(edges):
        rows[left].append(right)
        rows[right].append(left)
    return tuple(tuple(sorted(row)) for row in rows)


def k7() -> tuple[tuple[int, ...], ...]:
    return adjacency_from_edges(7, {(left, right) for left in range(7) for right in range(left + 1, 7)})


def k8_without_matching() -> tuple[tuple[int, ...], ...]:
    omitted = {(0, 1), (2, 3), (4, 5), (6, 7)}
    return adjacency_from_edges(
        8,
        {
            (left, right)
            for left in range(8)
            for right in range(left + 1, 8)
            if (left, right) not in omitted
        },
    )


def torus_3x3x3() -> tuple[tuple[int, ...], ...]:
    index = lambda x, y, z: (x * 3 + y) * 3 + z
    edges: set[tuple[int, int]] = set()
    for x, y, z in itertools.product(range(3), repeat=3):
        site = index(x, y, z)
        for neighbour in (
            index((x + 1) % 3, y, z),
            index(x, (y + 1) % 3, z),
            index(x, y, (z + 1) % 3),
        ):
            edges.add(tuple(sorted((site, neighbour))))
    return adjacency_from_edges(27, edges)


def mask(vertices: tuple[int, ...]) -> int:
    value = 0
    for vertex in vertices:
        value |= 1 << vertex
    return value


def fixed_point_closure(
    adjacency: tuple[tuple[int, ...], ...], sizes: tuple[int, ...]
) -> set[int]:
    family = {(1 << 0) | (1 << 1)}
    while True:
        enlarged = set(family)
        for support in family:
            for vertex, neighbours in enumerate(adjacency):
                if not (support >> vertex) & 1:
                    continue
                reduced = support ^ (1 << vertex)
                for size in sizes:
                    for subset in itertools.combinations(neighbours, size):
                        enlarged.add(reduced ^ mask(subset))
        if len(enlarged) == len(family):
            return family
        family = enlarged


def token_bfs(adjacency: tuple[tuple[int, ...], ...], tokens: int) -> int:
    start = mask(tuple(range(tokens)))
    seen = {start}
    queue = deque([start])
    while queue:
        support = queue.popleft()
        occupied = [site for site in range(len(adjacency)) if (support >> site) & 1]
        for site in occupied:
            for neighbour in adjacency[site]:
                if (support >> neighbour) & 1:
                    continue
                target = support ^ (1 << site) ^ (1 << neighbour)
                if target not in seen:
                    seen.add(target)
                    queue.append(target)
    return len(seen)

def valid_induction_witness(
    row: dict[str, object], adjacency: tuple[tuple[int, ...], ...]
) -> bool:
    pivot = int(row["pivot"])
    triple = tuple(int(vertex) for vertex in row["triple"])
    source_size = int(row["source_size"])
    if not (0 <= pivot < len(adjacency)):
        return False
    if len(triple) != 3 or len(set(triple)) != 3:
        return False
    if any(vertex not in adjacency[pivot] for vertex in triple):
        return False
    available = tuple(
        vertex
        for vertex in range(len(adjacency))
        if vertex not in {pivot, *triple}
    )
    if source_size - 1 > len(available):
        return False
    source = (pivot, *available[: source_size - 1])
    source_mask = mask(source)
    target_mask = (source_mask ^ (1 << pivot)) ^ mask(triple)
    target_size = target_mask.bit_count()
    return (
        int(row["target_size"]) == target_size == source_size + 2
        and not any((source_mask >> vertex) & 1 for vertex in triple)
        and not ((target_mask >> pivot) & 1)
        and bool(row["source_excludes_triple"])
        and bool(row["target_omits_pivot"])
    )


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    data = artifact["data"]
    check(
        "current producer hash",
        artifact["meta"]["source_sha256"]["experiments/e243_callen_termwise_closure.py"]
        == file_sha256(PRODUCER),
    )
    check(
        "current verifier hash",
        artifact["meta"]["source_sha256"]["tests/test_callen_termwise_closure.py"]
        == file_sha256(Path(__file__).resolve()),
    )
    check(
        "current e240 coefficient-source hash",
        artifact["meta"]["source_sha256"][
            "experiments/e240_callen_identity_system.py"
        ]
        == file_sha256(BASE_CALLEN_PRODUCER),
    )

    coefficients = data["Callen_coefficients"]
    check(
        "exact nonzero c1 and c3 formulas",
        coefficients["c1"]["numerator_coefficients_ascending"]
        == [0, 1, 0, 16, 0, 46, 0, 16, 0, 1]
        and coefficients["c1"]["denominator_coefficients_ascending"]
        == [1, 0, 21, 0, 106, 0, 106, 0, 21, 0, 1]
        and coefficients["c3"]["numerator_coefficients_ascending"]
        == [0, 0, 0, -2]
        and coefficients["c3"]["denominator_coefficients_ascending"]
        == [1, 0, 15, 0, 15, 0, 1],
        "positive denominator coefficients give c1>0 and c3<0 on 0<v<1",
    )

    for name, adjacency in (("K7", k7()), ("K8_minus_perfect_matching", k8_without_matching())):
        order = len(adjacency)
        closure = fixed_point_closure(adjacency, (1, 3))
        expected = {
            support
            for support in range(1 << order)
            if support.bit_count() % 2 == 0 and support != (1 << order) - 1
        }
        singleton = fixed_point_closure(adjacency, (1,))
        stored = data["finite_audits"][name]
        check(
            f"{name} independent closure",
            closure == expected
            and len(closure) == int(stored["closure_size"])
            and len(singleton) == math.comb(order, 2) + 1
            and len(singleton) == int(stored["singleton_only_size"]),
            f"full={len(closure)} singleton-only={len(singleton)}",
        )

    cubic = torus_3x3x3()
    check(
        "independent periodic graph",
        len(cubic) == 27
        and sum(map(len, cubic)) // 2 == 81
        and all(len(row) == 6 for row in cubic),
    )
    observed_tokens = {str(tokens): token_bfs(cubic, tokens) for tokens in (2, 4)}
    check(
        "independent token connectivity",
        observed_tokens == {"2": math.comb(27, 2), "4": math.comb(27, 4)}
        == data["periodic_3x3x3"]["token_component_sizes"],
        str(observed_tokens),
    )

    witnesses = data["periodic_3x3x3"]["induction_witnesses"]
    witness_sizes = [int(row["source_size"]) for row in witnesses]
    check(
        "independent induction witnesses",
        witness_sizes == list(range(2, 26, 2))
        and all(valid_induction_witness(row, cubic) for row in witnesses),
    )

    periodic_rows = data["periodic_cubic_rows"]
    formulas_hold = True
    for row in periodic_rows:
        order = int(row["order"])
        shape = tuple(int(length) for length in row["shape"])
        expected_count = (1 << (order - 1)) - (1 if order % 2 == 0 else 0)
        expected_nontrivial = expected_count - 1
        expected_support_orbits = (expected_count + order - 1) // order
        expected_nontrivial_orbits = (expected_nontrivial + order - 1) // order
        formulas_hold &= (
            len(shape) == 3
            and all(length >= 3 for length in shape)
            and order == math.prod(shape)
            and int(row["translation_group_order"]) == order
            and int(row["termwise_support_count_including_empty"]) == expected_count
            and int(row["nontrivial_correlator_support_count"])
            == expected_nontrivial
            and int(row["support_orbit_lower_bound_including_empty"])
            == expected_support_orbits
            and int(row["nontrivial_correlator_orbit_lower_bound"])
            == expected_nontrivial_orbits
            and int(row["largest_forced_support_size"])
            == order - (2 if order % 2 == 0 else 1)
        )
    check(
        "periodic orbit lower bounds",
        formulas_hold
        and int(periodic_rows[0]["termwise_support_count_including_empty"])
        == 1 << 26
        and int(periodic_rows[0]["nontrivial_correlator_orbit_lower_bound"])
        == 2485514,
    )

    theorem = data["theorem"]
    check(
        "scope preserves compression boundary",
        theorem["scope"]
        == "all-pivot {1,3}-termwise support closure using nonzero c1 and c3 Callen summands in the uniform zero-field nearest-neighbour model"
        and "no-go for schemes retaining only a selected subset of pivot rows"
        in theorem["not_proved"]
        and "no-go for arbitrary linear combinations of Callen rows"
        in theorem["not_proved"]
        and "no-go for transforms that aggregate distinct support orbits"
        in theorem["not_proved"],
    )
    check(
        "RSS wall",
        peak_rss_bytes() < RSS_LIMIT_BYTES,
        f"{peak_rss_bytes()}/{RSS_LIMIT_BYTES} bytes",
    )
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: Callen termwise-support closure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
