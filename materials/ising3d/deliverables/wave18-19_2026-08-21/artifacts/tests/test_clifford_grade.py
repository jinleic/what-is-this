"""Independent exact verifier for the Clifford-grade refutation/trichotomy.

No e178/e179/e180/e181 module and no project Clifford helper is imported.  The
five repository dimensions are rebuilt from the raw X_i and Z_u Z_v generator
bit-vectors by a separate left-normed bracket BFS.
"""

from __future__ import annotations

import hashlib
import json
import platform
import resource
import time
from collections import Counter, deque
from math import comb
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "clifford_grade.json"
CPU_BUDGET_SECONDS = 180.0
RSS_CAP_BYTES = 4_000_000_000
STARTED = time.process_time()
PASSED: list[str] = []
FAILED: list[str] = []


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(stage: str) -> None:
    used = time.process_time() - STARTED
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(
            f"NON-DECISIVE resource expiry: {stage} used {used:.3f}s > {CPU_BUDGET_SECONDS}s"
        )
    if max_rss_bytes() >= RSS_CAP_BYTES:
        raise RuntimeError(
            f"NON-DECISIVE resource expiry: {stage} reached {max_rss_bytes()} bytes"
        )


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(f"[{'PASS' if condition else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    (PASSED if condition else FAILED).append(name)


def normalize_edges(
    n: int, edges: Iterable[tuple[int, int]]
) -> tuple[tuple[int, int], ...]:
    values: set[tuple[int, int]] = set()
    for u, v in edges:
        if not 0 <= u < n or not 0 <= v < n or u == v:
            raise ValueError((n, u, v))
        values.add((u, v) if u < v else (v, u))
    return tuple(sorted(values))


def grid_edges(rows: int, cols: int) -> tuple[tuple[int, int], ...]:
    n = rows * cols
    edges: list[tuple[int, int]] = []
    for row in range(rows):
        for col in range(cols):
            site = row * cols + col
            if row + 1 < rows:
                edges.append((site, site + cols))
            if col + 1 < cols:
                edges.append((site, site + 1))
    return normalize_edges(n, edges)


def snake_path(rows: int, cols: int) -> tuple[int, ...]:
    values: list[int] = []
    for row in range(rows):
        columns = range(cols) if row % 2 == 0 else range(cols - 1, -1, -1)
        values.extend(row * cols + col for col in columns)
    return tuple(values)


def cycle_edges(n: int) -> tuple[tuple[int, int], ...]:
    return normalize_edges(
        n, [(site, site + 1) for site in range(n - 1)] + [(n - 1, 0)]
    )


def raw_generators(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    fields = [1 << site for site in range(n)]
    bonds = [(1 << (n + u)) | (1 << (n + v)) for u, v in normalize_edges(n, edges)]
    return tuple(fields + bonds)


def anticommutes(left: int, right: int, n: int) -> bool:
    mask = (1 << n) - 1
    x_left, z_left = left & mask, left >> n
    x_right, z_right = right & mask, right >> n
    return (
        (x_left & z_right).bit_count() + (z_left & x_right).bit_count()
    ) % 2 == 1


def raw_closure(n: int, edges: Sequence[tuple[int, int]]) -> frozenset[int]:
    """Exact closure under ad of the raw generators (left-normed Lie words)."""
    generators = raw_generators(n, edges)
    seen = set(generators)
    queue = deque(generators)
    visits = 0
    while queue:
        current = queue.popleft()
        for generator in generators:
            if anticommutes(current, generator, n):
                candidate = current ^ generator
                if candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
        visits += 1
        if visits % 16_384 == 0:
            budget_tick(f"raw closure n={n}, seen={len(seen)}")
    return frozenset(seen)


def digest(values: Iterable[int]) -> str:
    answer = hashlib.sha256()
    for value in sorted(values):
        answer.update(f"{value:x}\n".encode("ascii"))
    return answer.hexdigest()


def path_majorana_masks(n: int, path: Sequence[int]) -> tuple[int, ...]:
    prefix_x = 0
    values: list[int] = []
    for site in path:
        z = 1 << (n + site)
        values.extend((prefix_x | z, prefix_x | (1 << site) | z))
        prefix_x |= 1 << site
    return tuple(values)


def inverse_majorana_basis(
    n: int, path: Sequence[int]
) -> dict[int, tuple[int, int]]:
    rows: dict[int, tuple[int, int]] = {}
    for index, mask in enumerate(path_majorana_masks(n, path)):
        coefficients = 1 << index
        value = mask
        while value:
            pivot = value.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = (value, coefficients)
                break
            row, row_coefficients = rows[pivot]
            value ^= row
            coefficients ^= row_coefficients
        if not value:
            raise AssertionError("dependent Majorana basis")
    return rows


def majorana_grade(
    pauli: int, inverse: dict[int, tuple[int, int]]
) -> int:
    value = pauli
    coefficients = 0
    while value:
        pivot = value.bit_length() - 1
        row, row_coefficients = inverse[pivot]
        value ^= row
        coefficients ^= row_coefficients
    return coefficients.bit_count()


def grade_histogram(
    closure: Iterable[int], n: int, path: Sequence[int]
) -> dict[int, int]:
    inverse = inverse_majorana_basis(n, path)
    return dict(
        sorted(Counter(majorana_grade(value, inverse) for value in closure).items())
    )


def noncentral_two_mod_four_grades(n: int) -> tuple[int, ...]:
    N = 2 * n
    return tuple(grade for grade in range(2, N + 1, 4) if grade != N)


def two_mod_four_dimension(n: int) -> int:
    return sum(comb(2 * n, grade) for grade in noncentral_two_mod_four_grades(n))


def expected_histogram(n: int, branch: str) -> dict[int, int]:
    N = 2 * n
    if branch == "path":
        return {2: comb(N, 2)}
    if branch == "cycle":
        return {2: comb(N, 2), N - 2: comb(N, N - 2)}
    if branch == "branching":
        return {grade: comb(N, grade) for grade in noncentral_two_mod_four_grades(n)}
    if branch == "full_even_derived":
        return {grade: comb(N, grade) for grade in range(2, N, 2)}
    raise ValueError(branch)


def pairwise_closed(values: frozenset[int], n: int) -> bool:
    ordered = tuple(values)
    for index, left in enumerate(ordered):
        for right in ordered[index + 1 :]:
            if anticommutes(left, right, n) and (left ^ right) not in values:
                return False
    return True


def validate_path_and_chords(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> tuple[int, ...]:
    edge_set = set(normalize_edges(n, edges))
    if len(path) != n or set(path) != set(range(n)):
        raise AssertionError("not Hamiltonian")
    path_set = set()
    for u, v in zip(path, path[1:]):
        edge = tuple(sorted((u, v)))
        if edge not in edge_set:
            raise AssertionError(f"missing path edge {edge}")
        path_set.add(edge)
    position = {vertex: index for index, vertex in enumerate(path)}
    return tuple(
        2 * abs(position[u] - position[v])
        for u, v in sorted(edge_set - path_set)
    )


def is_bipartite(n: int, edges: Sequence[tuple[int, int]]) -> bool:
    adjacency = [[] for _ in range(n)]
    for u, v in normalize_edges(n, edges):
        adjacency[u].append(v)
        adjacency[v].append(u)
    colors = [-1] * n
    for root in range(n):
        if colors[root] >= 0:
            continue
        colors[root] = 0
        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if colors[v] < 0:
                    colors[v] = colors[u] ^ 1
                    queue.append(v)
                elif colors[v] == colors[u]:
                    return False
    return True


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    ok("artifact schema", set(artifact) == {"meta", "data", "checks"})
    ok(
        "producer checks decisive",
        bool(artifact["checks"]) and all(bool(row["passed"]) for row in artifact["checks"]),
        f"{len(artifact['checks'])} top-level checks",
    )
    ok(
        "artifact source digests",
        all(
            hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected
            for relative, expected in artifact["meta"]["source_sha256"].items()
        ),
    )

    stored_rows = {
        row["label"]: row
        for row in artifact["data"]["graph_verification"]["exact_rows"]
    }
    known_specs = (
        ("known_grid_2x2_C4", 2, 2, 56, "cycle"),
        ("known_grid_2x3", 2, 3, 1056, "branching"),
        ("known_grid_2x4", 2, 4, 16256, "branching"),
        ("known_grid_2x5", 2, 5, 262656, "branching"),
        ("known_grid_3x3", 3, 3, 65535, "branching"),
    )
    independently_measured: dict[str, int] = {}
    for label, rows, cols, expected_dimension, branch in known_specs:
        n = rows * cols
        edges = grid_edges(rows, cols)
        path = snake_path(rows, cols)
        closure = raw_closure(n, edges)
        histogram = grade_histogram(closure, n, path)
        independently_measured[label] = len(closure)
        ok(
            f"{label}: raw-generator dimension",
            len(closure) == expected_dimension,
            f"measured {len(closure)}",
        )
        ok(
            f"{label}: independent path-grade histogram",
            histogram == expected_histogram(n, branch),
            str(histogram),
        )
        stored = stored_rows[label]
        ok(
            f"{label}: artifact agreement",
            stored["dimension"] == len(closure)
            and stored["closure_sha256"] == digest(closure),
        )
        budget_tick(label)
    ok(
        "all five known values exactly reproduced",
        list(independently_measured.values()) == [56, 1056, 16256, 262656, 65535],
        str(independently_measured),
    )

    # Trichotomy controls: the path, two genuine cycles, and C4's coincidence.
    control_specs = (
        ("path_P6_control", 6, tuple((i, i + 1) for i in range(5)), 66, "path"),
        ("cycle_C6_counterexample", 6, cycle_edges(6), 132, "cycle"),
        ("cycle_C8_counterexample", 8, cycle_edges(8), 240, "cycle"),
    )
    for label, n, edges, expected_dimension, branch in control_specs:
        path = tuple(range(n))
        closure = raw_closure(n, edges)
        histogram = grade_histogram(closure, n, path)
        ok(
            f"{label}: exact trichotomy control",
            len(closure) == expected_dimension
            and histogram == expected_histogram(n, branch),
            f"dimension={len(closure)}, grades={histogram}",
        )
        ok(
            f"{label}: stored closure digest",
            stored_rows[label]["closure_sha256"] == digest(closure),
        )
        if branch == "cycle":
            ok(
                f"{label}: independently pairwise saturated",
                pairwise_closed(closure, n),
            )
    ok(
        "C6 refutes proposed classification",
        stored_rows["cycle_C6_counterexample"]["dimension"] == 132
        and two_mod_four_dimension(6) == 1056,
        "132=2*C(12,2), not 1056",
    )
    ok(
        "C4 is cycle branch but accidental full class",
        expected_histogram(4, "cycle") == expected_histogram(4, "branching")
        == {2: 28, 6: 28},
    )

    # New exact prediction: remove a non-path edge from the 3x3 grid.
    new_edges = tuple(edge for edge in grid_edges(3, 3) if edge != (0, 3))
    new_path = snake_path(3, 3)
    new_closure = raw_closure(9, new_edges)
    new_histogram = grade_histogram(new_closure, 9, new_path)
    ok(
        "new 3x3-minus-edge raw closure",
        len(new_closure) == 65535
        and new_histogram == expected_histogram(9, "branching"),
        f"dimension={len(new_closure)}, grades={new_histogram}",
    )
    ok(
        "new 3x3-minus-edge artifact agreement",
        stored_rows["new_grid_3x3_minus_edge_0_3"]["closure_sha256"]
        == digest(new_closure),
    )

    # Independent non-materialized predictions from only graph/path/binomial data.
    grid_3x4_edges = grid_edges(3, 4)
    grid_3x4_path = snake_path(3, 4)
    grid_3x4_grades = validate_path_and_chords(12, grid_3x4_edges, grid_3x4_path)
    grid_3x4_prediction = two_mod_four_dimension(12)
    ok(
        "new 3x4 prediction independently derived",
        is_bipartite(12, grid_3x4_edges)
        and 6 in grid_3x4_grades
        and all(grade % 4 == 2 for grade in grid_3x4_grades)
        and grid_3x4_prediction == 4_192_256,
        f"chord grades={grid_3x4_grades}, dimension={grid_3x4_prediction}",
    )
    characterized = {
        row["label"]: row
        for row in artifact["data"]["graph_verification"]["characterization_rows"]
    }
    ok(
        "new 3x4 artifact agreement",
        characterized["new_grid_3x4"]["dimension"] == grid_3x4_prediction,
    )

    c10_chord_edges = normalize_edges(10, (*cycle_edges(10), (0, 3)))
    c10_grades = validate_path_and_chords(10, c10_chord_edges, tuple(range(10)))
    c10_prediction = two_mod_four_dimension(10)
    ok(
        "new non-grid C10-plus-chord prediction",
        is_bipartite(10, c10_chord_edges)
        and set(c10_grades) == {6, 18}
        and c10_prediction == 262656
        and characterized["new_non_grid_C10_plus_chord_0_3"]["dimension"]
        == c10_prediction,
        f"chord grades={c10_grades}, dimension={c10_prediction}",
    )

    # Hypothesis-necessity falsification: odd cycle plus chord.
    nonbip_edges = normalize_edges(5, (*cycle_edges(5), (0, 3)))
    nonbip_path = tuple(range(5))
    nonbip_grades = validate_path_and_chords(5, nonbip_edges, nonbip_path)
    nonbip_closure = raw_closure(5, nonbip_edges)
    nonbip_histogram = grade_histogram(nonbip_closure, 5, nonbip_path)
    ok(
        "non-bipartite C5-plus-chord exact falsification",
        not is_bipartite(5, nonbip_edges)
        and set(nonbip_grades) == {6, 8}
        and len(nonbip_closure) == 510
        and nonbip_histogram == expected_histogram(5, "full_even_derived"),
        f"dimension={len(nonbip_closure)}, grades={nonbip_histogram}",
    )
    ok(
        "non-bipartite comparison values",
        two_mod_four_dimension(5) == 255
        and sum(comb(10, grade) for grade in range(2, 10, 2)) == 510,
        "measured 510; 2-mod-4 255; full noncentral even 510",
    )
    ok(
        "non-bipartite closure pairwise saturated",
        pairwise_closed(nonbip_closure, 5),
    )
    ok(
        "non-bipartite artifact digest",
        stored_rows["nonbipartite_C5_plus_chord_0_3"]["closure_sha256"]
        == digest(nonbip_closure),
    )

    # The central-volume correction is about generation, not matrix tracelessness.
    global_x_pauli = (1 << 9) - 1
    odd_volume_grade = 18
    ok(
        "odd-n volume is non-scalar on full physical space",
        global_x_pauli != 0
        and global_x_pauli not in raw_closure(9, grid_edges(3, 3))
        and odd_volume_grade not in expected_histogram(9, "branching"),
        "global X parity is a nonidentity Pauli string and is not generated",
    )
    ok(
        "artifact records corrected volume convention",
        artifact["data"]["centre_and_representation_convention"]["full_space_scalar"]
        is False
        and artifact["data"]["centre_and_representation_convention"]["half_spin_scalar"]
        is True,
    )

    budget_tick("final")
    ok("RSS cap", max_rss_bytes() < RSS_CAP_BYTES, f"peak={max_rss_bytes()} bytes")
    if FAILED:
        print(f"FAIL: {len(FAILED)} failed checks: {FAILED}")
        return 1
    print(
        f"PASS ({len(PASSED)} independent checks, "
        f"process_time={time.process_time()-STARTED:.3f}s, peak_rss={max_rss_bytes()})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
