"""Clean-room exact verifier for proofs/trichotomy_extend.md.

No experiment module or project Clifford helper is imported.  Every finite
closure starts from raw X_v and Z_u Z_v Pauli labels.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import time
from collections import Counter, deque
from math import comb, factorial
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "trichotomy_extend.json"
CPU_BUDGET_SECONDS = 180.0
RSS_CAP_BYTES = 2_000_000_000
STARTED = time.process_time()
PASSED: list[str] = []
FAILED: list[str] = []


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget_tick(stage: str) -> None:
    used = time.process_time() - STARTED
    if used > CPU_BUDGET_SECONDS:
        raise RuntimeError(f"process-time budget exceeded at {stage}: {used}")
    if max_rss_bytes() >= RSS_CAP_BYTES:
        raise MemoryError(f"RSS cap exceeded at {stage}: {max_rss_bytes()}")


def ok(name: str, condition: bool, detail: str = "") -> None:
    print(
        f"[{'PASS' if condition else 'FAIL'}] {name}"
        + (f": {detail}" if detail else ""),
        flush=True,
    )
    (PASSED if condition else FAILED).append(name)


def normalize_edges(
    n: int, edges: Iterable[tuple[int, int]]
) -> tuple[tuple[int, int], ...]:
    answer: set[tuple[int, int]] = set()
    for left, right in edges:
        if not (0 <= left < n and 0 <= right < n) or left == right:
            raise ValueError((n, left, right))
        answer.add((left, right) if left < right else (right, left))
    return tuple(sorted(answer))


def adjacency_masks(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    answer = [0] * n
    for left, right in normalize_edges(n, edges):
        answer[left] |= 1 << right
        answer[right] |= 1 << left
    return tuple(answer)


def graph_properties(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[bool, bool, tuple[int, ...], int]:
    adjacency = adjacency_masks(n, edges)
    colours = [-1] * n
    components = 0
    bipartite = True
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
    return (
        components == 1,
        bipartite,
        tuple(colours),
        max(mask.bit_count() for mask in adjacency),
    )


def raw_generators(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(
        [1 << vertex for vertex in range(n)]
        + [
            (1 << (n + left)) | (1 << (n + right))
            for left, right in normalize_edges(n, edges)
        ]
    )


def anticommutes(left: int, right: int, n: int) -> bool:
    mask = (1 << n) - 1
    x_left, z_left = left & mask, left >> n
    x_right, z_right = right & mask, right >> n
    return (
        (x_left & z_right).bit_count() + (z_left & x_right).bit_count()
    ) & 1 == 1


def raw_closure(n: int, edges: Sequence[tuple[int, int]]) -> frozenset[int]:
    generators = raw_generators(n, edges)
    seen = set(generators)
    queue = deque(generators)
    while queue:
        parent = queue.popleft()
        for generator in generators:
            if anticommutes(parent, generator, n):
                child = parent ^ generator
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
        if anticommutes(left, right, n)
    )


def digest(values: Iterable[int]) -> str:
    answer = hashlib.sha256()
    for value in sorted(values):
        answer.update(value.to_bytes(16, "little"))
    return answer.hexdigest()


def hamiltonian_path_exists(n: int, edges: Sequence[tuple[int, int]]) -> bool:
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


def edge_grades(
    n: int, edges: Sequence[tuple[int, int]], order: Sequence[int]
) -> tuple[int, ...]:
    if len(order) != n or set(order) != set(range(n)):
        raise ValueError("not an order")
    position = {vertex: index for index, vertex in enumerate(order)}
    return tuple(
        2 * abs(position[left] - position[right])
        for left, right in normalize_edges(n, edges)
    )


def compatible_order_count_formula(part_sizes: Sequence[int]) -> int:
    left, right = sorted(part_sizes)
    if right - left > 1:
        return 0
    if left == right:
        return 2 * factorial(left) * factorial(right)
    return factorial(left) * factorial(right)


def compatible_order_count(
    n: int, edges: Sequence[tuple[int, int]]
) -> int:
    return sum(
        all(grade % 4 == 2 for grade in edge_grades(n, edges, order))
        for order in itertools.permutations(range(n))
    )


def path_majoranas(n: int, order: Sequence[int]) -> tuple[int, ...]:
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
    for index, value in enumerate(path_majoranas(n, order)):
        coefficients = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in rows:
                rows[pivot] = value, coefficients
                break
            row, row_coefficients = rows[pivot]
            value ^= row
            coefficients ^= row_coefficients
        else:
            raise AssertionError("dependent Majoranas")
    return rows


def majorana_grade(value: int, inverse: dict[int, tuple[int, int]]) -> int:
    coefficients = 0
    while value:
        pivot = value.bit_length() - 1
        row, row_coefficients = inverse[pivot]
        value ^= row
        coefficients ^= row_coefficients
    return coefficients.bit_count()


def grade_histogram(
    values: Iterable[int], n: int, order: Sequence[int]
) -> dict[int, int]:
    inverse = inverse_majorana_basis(n, order)
    return dict(
        sorted(Counter(majorana_grade(value, inverse) for value in values).items())
    )


def full_even_histogram(n: int) -> dict[int, int]:
    return {grade: comb(2 * n, grade) for grade in range(2, 2 * n, 2)}


def full_even_dimension(n: int) -> int:
    return (1 << (2 * n - 1)) - 2


def bracket_outputs(N: int, left: int, right: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                left + right - 2 * overlap
                for overlap in range(1, min(left, right) + 1, 2)
                if left + right - overlap <= N
                and 0 < left + right - 2 * overlap < N
            }
        )
    )


def grade_closure(N: int, seeds: Iterable[int]) -> tuple[int, ...]:
    grades = set(seeds)
    changed = True
    while changed:
        changed = False
        for left in tuple(grades):
            for right in tuple(grades):
                for output in bracket_outputs(N, left, right):
                    if output not in grades:
                        grades.add(output)
                        changed = True
    return tuple(sorted(grades))


def grid_edges(rows: int, cols: int) -> tuple[tuple[int, int], ...]:
    edges: list[tuple[int, int]] = []
    for row in range(rows):
        for col in range(cols):
            vertex = row * cols + col
            if col + 1 < cols:
                edges.append((vertex, vertex + 1))
            if row + 1 < rows:
                edges.append((vertex, vertex + cols))
    return normalize_edges(rows * cols, edges)


def grid_path(rows: int, cols: int) -> tuple[int, ...]:
    answer: list[int] = []
    for row in range(rows):
        columns = range(cols) if row % 2 == 0 else range(cols - 1, -1, -1)
        answer.extend(row * cols + col for col in columns)
    return tuple(answer)


def box_vertex(a: int, b: int, x: int, y: int, z: int) -> int:
    return (z * b + y) * a + x


def box_edges(a: int, b: int, c: int) -> tuple[tuple[int, int], ...]:
    edges: list[tuple[int, int]] = []
    for z in range(c):
        for y in range(b):
            for x in range(a):
                vertex = box_vertex(a, b, x, y, z)
                if x + 1 < a:
                    edges.append((vertex, box_vertex(a, b, x + 1, y, z)))
                if y + 1 < b:
                    edges.append((vertex, box_vertex(a, b, x, y + 1, z)))
                if z + 1 < c:
                    edges.append((vertex, box_vertex(a, b, x, y, z + 1)))
    return normalize_edges(a * b * c, edges)


def box_path(a: int, b: int, c: int) -> tuple[int, ...]:
    first_layer: list[int] = []
    for y in range(b):
        xs = range(a) if y % 2 == 0 else range(a - 1, -1, -1)
        first_layer.extend(box_vertex(a, b, x, y, 0) for x in xs)
    answer: list[int] = []
    for z in range(c):
        layer = first_layer if z % 2 == 0 else list(reversed(first_layer))
        answer.extend(vertex + z * a * b for vertex in layer)
    return tuple(answer)


def validate_path(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> bool:
    edge_set = set(normalize_edges(n, edges))
    return (
        len(path) == n
        and set(path) == set(range(n))
        and all(
            ((left, right) if left < right else (right, left)) in edge_set
            for left, right in zip(path, path[1:])
        )
    )


def chord_grades(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> tuple[int, ...]:
    if not validate_path(n, edges, path):
        raise ValueError("invalid Hamiltonian path")
    used = {
        (left, right) if left < right else (right, left)
        for left, right in zip(path, path[1:])
    }
    position = {vertex: index for index, vertex in enumerate(path)}
    return tuple(
        2 * abs(position[left] - position[right])
        for left, right in normalize_edges(n, edges)
        if (left, right) not in used
    )


def branch_for_sides(sides: Sequence[int]) -> str:
    nontrivial = sorted(side for side in sides if side > 1)
    if len(nontrivial) <= 1:
        return "path"
    if nontrivial == [2, 2]:
        return "cycle"
    return "branching"


def branching_dimension(n: int) -> int:
    direct = sum(
        comb(2 * n, grade)
        for grade in range(2, 2 * n + 1, 4)
        if grade != 2 * n
    )
    closed = (
        (1 << (2 * n - 2)) - 1
        if n % 2
        else (1 << (2 * n - 2))
        - ((-1) ** (n // 2)) * (1 << (n - 1))
    )
    if direct != closed:
        raise AssertionError((n, direct, closed))
    return direct


def branch_dimension(branch: str, n: int) -> int:
    if branch == "path":
        return n * (2 * n - 1)
    if branch == "cycle":
        return 2 * n * (2 * n - 1)
    if branch == "branching":
        return branching_dimension(n)
    raise ValueError(branch)


def evidence_path_exists(reference: str) -> bool:
    relative = reference.split(":", 1)[0]
    return (ROOT / relative).exists()


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    ok("artifact schema", set(artifact) == {"meta", "data", "checks"})
    ok(
        "producer checks decisive",
        bool(artifact["checks"])
        and all(bool(check["passed"]) for check in artifact["checks"]),
        f"{len(artifact['checks'])} checks",
    )
    ok(
        "artifact source digests",
        all(
            hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected
            for path, expected in artifact["meta"]["source_sha256"].items()
        ),
    )

    # Arbitrary-order criterion and both exact bipartite counterexamples.
    stored_counterexamples = artifact["data"]["arbitrary_order"]["counterexamples"]
    counterexample_dimensions: dict[str, int] = {}
    for row in stored_counterexamples:
        n = int(row["n"])
        edges = tuple(tuple(edge) for edge in row["edges"])
        closure = raw_closure(n, edges)
        connected, bipartite, colours, degree = graph_properties(n, edges)
        measured_orders = compatible_order_count(n, edges)
        predicted_orders = compatible_order_count_formula(
            (colours.count(0), colours.count(1))
        )
        counterexample_dimensions[row["label"]] = len(closure)
        ok(
            f"{row['label']}: raw-generator counterexample",
            connected
            and bipartite
            and degree >= 3
            and len(closure) == row["dimension"]
            and digest(closure) == row["closure_sha256"]
            and measured_orders == predicted_orders == row["compatible_order_count"]
            and hamiltonian_path_exists(n, edges)
            == row["hamiltonian_path_exists"],
            f"dimension={len(closure)}, compatible_orders={measured_orders}",
        )
    claw_edges = ((0, 1), (0, 2), (0, 3))
    claw_closure = raw_closure(4, claw_edges)
    ok(
        "K1,3 refutes Hamiltonian-free branching formula",
        counterexample_dimensions["bipartite_claw_K1_3"] == 72
        and artifact["data"]["arbitrary_order"]["counterexamples"][0][
            "hamiltonian_branching_formula_if_misapplied"
        ]
        == 56
        and pairwise_closed(claw_closure, 4),
        "72 != 56, with exact pairwise saturation",
    )
    ok(
        "second mismatch is independently reproduced",
        counterexample_dimensions["bipartite_branch_tree_T6"] == 992,
        "T6 gives 992 rather than 1056",
    )

    # Exhaust the order criterion independently through n=4.
    order_graphs = 0
    for n in range(2, 5):
        pairs = tuple(itertools.combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            edges = tuple(
                pairs[index]
                for index in range(len(pairs))
                if (mask >> index) & 1
            )
            connected, bipartite, colours, _ = graph_properties(n, edges)
            if not connected or not bipartite:
                continue
            measured = compatible_order_count(n, edges)
            expected = compatible_order_count_formula(
                (colours.count(0), colours.count(1))
            )
            if measured != expected:
                FAILED.append(f"order criterion n={n}, edges={edges}")
            order_graphs += 1
    ok(
        "all connected bipartite graphs through n=4 obey order criterion",
        not any(name.startswith("order criterion n=") for name in FAILED),
        f"{order_graphs} graphs",
    )
    budget_tick("arbitrary orders")

    # Every stored non-bipartite row is rebuilt from raw local generators.
    nonbip_rows = artifact["data"]["nonbipartite"]["graph_rows"]
    nonbip_dimensions: dict[str, int] = {}
    for row in nonbip_rows:
        n = int(row["n"])
        edges = tuple(tuple(edge) for edge in row["edges"])
        closure = raw_closure(n, edges)
        order = (
            tuple(row["chosen_hamiltonian_order"])
            if row["chosen_hamiltonian_order"] is not None
            else tuple(range(n))
        )
        histogram = grade_histogram(closure, n, order)
        connected, bipartite, _, degree = graph_properties(n, edges)
        expected_dimension = (
            2 * comb(2 * n, 2)
            if row["branch"] == "odd_cycle"
            else full_even_dimension(n)
        )
        expected_histogram = (
            {2: comb(2 * n, 2), 2 * n - 2: comb(2 * n, 2)}
            if row["branch"] == "odd_cycle"
            else full_even_histogram(n)
        )
        nonbip_dimensions[row["label"]] = len(closure)
        ok(
            f"{row['label']}: independent non-bipartite closure",
            connected
            and not bipartite
            and (degree == 2 if row["branch"] == "odd_cycle" else degree >= 3)
            and len(closure) == expected_dimension == row["dimension"]
            and histogram == expected_histogram
            and digest(closure) == row["closure_sha256"],
            f"dimension={len(closure)}, grades={histogram}",
        )
    ok(
        "at least two non-bipartite dimensions re-derived",
        nonbip_dimensions["paw_C3_plus_leaf"] == 126
        and nonbip_dimensions["C5_plus_chord_0_3"] == 510
        and nonbip_dimensions["C6_plus_even_distance_chord_0_2"] == 2046,
        str(nonbip_dimensions),
    )
    ok(
        "non-Hamiltonian leaf-extension control",
        nonbip_dimensions["nonhamiltonian_triangle_two_leaves"] == 510
        and not hamiltonian_path_exists(
            5, ((0, 1), (0, 2), (1, 2), (0, 3), (0, 4))
        ),
    )
    ok(
        "bare odd cycle is the required exception",
        nonbip_dimensions["odd_cycle_C5_exception"] == 90
        and full_even_dimension(5) == 510,
        "C5=90, not 510",
    )
    budget_tick("non-bipartite closures")

    # Independently exhaust both grade-generation lemmas to the producer bound.
    interior_cases = 0
    endpoint_cases = 0
    for N in range(8, 81, 2):
        target = tuple(range(2, N, 2))
        for seed in range(4, N - 3, 4):
            if grade_closure(N, (2, seed)) != target:
                FAILED.append(f"interior grade seed N={N}, seed={seed}")
            interior_cases += 1
        if N % 4 == 2:
            for seed in range(6, N - 3, 4):
                if grade_closure(N, (2, seed, N - 2)) != target:
                    FAILED.append(f"endpoint grade seed N={N}, seed={seed}")
                endpoint_cases += 1
    stored_support = artifact["data"]["nonbipartite"]["support_generation_audit"]
    ok(
        "full-even support lemmas independently exhausted",
        not any("grade seed" in name for name in FAILED)
        and interior_cases == stored_support["interior_zero_mod_four_cases"]
        and endpoint_cases == stored_support["endpoint_plus_two_mod_four_cases"],
        f"interior={interior_cases}, endpoint={endpoint_cases}",
    )

    # Every concrete corollary row must have a raw graph/path/grade certificate.
    concrete = artifact["data"]["corollaries"]["concrete_rows"]
    expected_shapes = {
        "grid_2x6": (2, 6),
        "grid_3x4": (3, 4),
        "grid_3x5": (3, 5),
        "grid_4x4": (4, 4),
    }
    supported_labels: set[str] = set()
    derived_dimensions: dict[str, int] = {}
    for row in concrete:
        label = row["label"]
        if label not in expected_shapes or row["geometry"] != "open_rectangular_grid":
            FAILED.append(f"unsupported corollary row {label}")
            continue
        rows, cols = expected_shapes[label]
        if tuple(row["shape"]) != (rows, cols):
            FAILED.append(f"wrong shape for {label}")
            continue
        n = rows * cols
        edges = grid_edges(rows, cols)
        path = grid_path(rows, cols)
        generators = raw_generators(n, edges)
        grades = chord_grades(n, edges, path)
        connected, bipartite, _, degree = graph_properties(n, edges)
        branch = branch_for_sides((rows, cols))
        dimension = branch_dimension(branch, n)
        evidence_supported = all(
            evidence_path_exists(reference) for reference in row["previous_evidence"]
        )
        condition = (
            connected
            and bipartite
            and degree >= 3
            and branch == row["branch"] == "branching"
            and validate_path(n, edges, path)
            and tuple(row["hamiltonian_path"]) == path
            and tuple(row["chord_grades"]) == grades
            and all(grade % 4 == 2 for grade in grades)
            and any(6 <= grade <= 2 * n - 4 for grade in grades)
            and len(generators) == row["raw_generator_count"]
            and digest(generators) == row["raw_generator_sha256"]
            and dimension == row["dimension"] == row["dimension_by_binomial_sum"]
            and bool(row["previous_status"])
            and evidence_supported
            and row["two_generator_scope_separate"] is True
        )
        ok(
            f"{label}: raw-generator corollary",
            condition,
            f"dimension={dimension}, chord_grades={grades}",
        )
        if condition:
            supported_labels.add(label)
            derived_dimensions[label] = dimension
    ok(
        "no corollary row is unsupported",
        supported_labels == set(expected_shapes)
        and len(concrete) == len(expected_shapes)
        and not any(name.startswith("unsupported corollary") for name in FAILED),
        str(sorted(supported_labels)),
    )
    ok(
        "two new graph-specific dimensions re-derived",
        derived_dimensions["grid_3x5"] == 268_435_455
        and derived_dimensions["grid_4x4"] == 1_073_709_056,
        "3x5=268435455; 4x4=1073709056",
    )
    ok(
        "known n=12 corollaries remain correctly attributed",
        derived_dimensions["grid_2x6"] == 4_192_256
        and derived_dimensions["grid_3x4"] == 4_192_256
        and {
            row["label"]: row["previous_status_class"] for row in concrete
        }["grid_2x6"]
        == "previously_known"
        and {
            row["label"]: row["previous_status_class"] for row in concrete
        }["grid_3x4"]
        == "previously_known",
    )
    four_by_four = next(row for row in concrete if row["label"] == "grid_4x4")
    ok(
        "4x4 local-term/two-generator scopes stay separate",
        four_by_four["dimension"] == 1_073_709_056
        and "1,794" in four_by_four["previous_status"]
        and four_by_four["previous_status_class"]
        == "previously_implicit_two_generator_open",
    )

    # Audit the two all-size geometry corollaries on independently built graphs.
    families = artifact["data"]["corollaries"]
    rectangular = families["rectangular_grid_family"]
    cubic = families["simple_cubic_box_family"]
    ok(
        "family rows have complete branch support",
        {row["branch"] for row in rectangular["branch_cases"]}
        == {"path", "cycle", "branching"}
        and {row["branch"] for row in cubic["branch_cases"]}
        == {"path", "cycle", "branching"}
        and len(rectangular["branch_cases"]) == 3
        and len(cubic["branch_cases"]) == 3,
    )
    grid_cases = 0
    for rows in range(1, 6):
        for cols in range(1, 6):
            n = rows * cols
            edges = grid_edges(rows, cols)
            path = grid_path(rows, cols)
            branch = branch_for_sides((rows, cols))
            if not validate_path(n, edges, path):
                FAILED.append(f"grid family path {rows}x{cols}")
            if branch == "branching":
                grades = chord_grades(n, edges, path)
                if not all(grade % 4 == 2 for grade in grades):
                    FAILED.append(f"grid family grades {rows}x{cols}")
            branch_dimension(branch, n)
            grid_cases += 1
    box_cases = 0
    branch_box_cases = 0
    for a in range(1, 5):
        for b in range(1, 5):
            for c in range(1, 5):
                n = a * b * c
                if n > 24:
                    continue
                edges = box_edges(a, b, c)
                path = box_path(a, b, c)
                branch = branch_for_sides((a, b, c))
                connected, bipartite, _, _ = graph_properties(n, edges)
                if not connected or not bipartite or not validate_path(n, edges, path):
                    FAILED.append(f"box family graph {(a, b, c)}")
                if branch == "branching":
                    branch_box_cases += 1
                    grades = chord_grades(n, edges, path)
                    if not all(grade % 4 == 2 for grade in grades):
                        FAILED.append(f"box family grades {(a, b, c)}")
                branch_dimension(branch, n)
                box_cases += 1
    ok(
        "rectangular-grid family independently supported",
        not any(name.startswith("grid family") for name in FAILED),
        f"{grid_cases} shapes",
    )
    stored_box_audit = cubic["finite_path_audit"]
    ok(
        "simple-cubic box family independently supported",
        not any(name.startswith("box family") for name in FAILED)
        and box_cases == stored_box_audit["cases"]
        and branch_box_cases == stored_box_audit["branching_cases"],
        f"{box_cases} boxes, {branch_box_cases} branching",
    )
    ok(
        "three-dimensional scope exclusions are explicit",
        set(cubic["scope_exclusions"])
        == {
            "thermodynamic limit",
            "transfer spectrum",
            "two-generator algebra",
            "solvability of the three-dimensional Ising model",
        },
    )
    repository_rows = families["repository_audit"]
    ok(
        "repository novelty audit is complete and sourced",
        len(repository_rows) == 5
        and all(row["prior_status"] and row["evidence"] for row in repository_rows)
        and all(
            evidence_path_exists(reference)
            for row in repository_rows
            for reference in row["evidence"]
        ),
        f"{len(repository_rows)} prior-status rows",
    )

    budget_tick("final")
    ok("RSS cap", max_rss_bytes() < RSS_CAP_BYTES, f"peak={max_rss_bytes()}")
    if FAILED:
        print(f"FAIL: {len(FAILED)} failed checks: {FAILED}")
        return 1
    print(
        f"PASS ({len(PASSED)} independent checks, "
        f"process_time={time.process_time()-STARTED:.6f}s, "
        f"peak_rss={max_rss_bytes()})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
