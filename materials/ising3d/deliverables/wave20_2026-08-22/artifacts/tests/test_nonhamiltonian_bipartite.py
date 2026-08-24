"""Clean-room verifier for the non-Hamiltonian bipartite classification.

No experiment module is imported.  Graph isomorphism reduction, all-pairs
Pauli closure, quadratic-root enumeration, graph invariants, and Clifford
grades are independently reconstructed from raw fields and bonds.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import resource
import time
from collections import Counter, deque
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "nonhamiltonian_bipartite.json"
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


def pairs(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(itertools.combinations(range(n), 2))


def decoded_edges(n: int, code: int) -> tuple[tuple[int, int], ...]:
    return tuple(pair for i, pair in enumerate(pairs(n)) if (code >> i) & 1)


def packed_edges(n: int, edges: Iterable[tuple[int, int]]) -> int:
    index = {pair: i for i, pair in enumerate(pairs(n))}
    answer = 0
    for left, right in edges:
        edge = (left, right) if left < right else (right, left)
        answer |= 1 << index[edge]
    return answer


def adjacency(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    rows = [0] * n
    for left, right in edges:
        rows[left] |= 1 << right
        rows[right] |= 1 << left
    return tuple(rows)


def graph_colouring(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[bool, tuple[int, ...]] | None:
    rows = adjacency(n, edges)
    colour = [-1] * n
    colour[0] = 0
    queue = [0]
    for vertex in queue:
        live = rows[vertex]
        while live:
            bit = live & -live
            live ^= bit
            other = bit.bit_length() - 1
            if colour[other] < 0:
                colour[other] = colour[vertex] ^ 1
                queue.append(other)
            elif colour[other] == colour[vertex]:
                return None
    return all(value >= 0 for value in colour), tuple(colour)


@lru_cache(maxsize=None)
def vertex_permutations(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.permutations(range(n)))


def isomorphism_signature(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[int, ...]:
    edge_set = {tuple(sorted(edge)) for edge in edges}
    best: tuple[int, ...] | None = None
    for order in vertex_permutations(n):
        candidate = tuple(
            int(tuple(sorted((order[left], order[right]))) in edge_set)
            for left, right in pairs(n)
        )
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise AssertionError(n)
    return best


def enumerate_unlabeled(max_n: int) -> tuple[dict[int, int], dict[tuple[int, tuple[int, ...]], tuple[tuple[int, int], ...]]]:
    labeled: dict[int, int] = {}
    representatives: dict[
        tuple[int, tuple[int, ...]], tuple[tuple[int, int], ...]
    ] = {}
    for n in range(2, max_n + 1):
        count = 0
        for code in range(1 << len(pairs(n))):
            if code.bit_count() < n - 1:
                continue
            edges = decoded_edges(n, code)
            result = graph_colouring(n, edges)
            if result is None or not result[0]:
                continue
            count += 1
            signature = isomorphism_signature(n, edges)
            representatives.setdefault((n, signature), edges)
            if code and code % 4096 == 0:
                budget_tick(f"independent graph enumeration n={n}")
        labeled[n] = count
    return labeled, representatives


def raw_generators(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(1 << vertex for vertex in range(n)) + tuple(
        (1 << (n + left)) | (1 << (n + right)) for left, right in edges
    )


def bracket_nonzero(left: int, right: int, n: int) -> bool:
    low = (1 << n) - 1
    return bool(
        (
            ((left & low) & (right >> n)).bit_count()
            + ((left >> n) & (right & low)).bit_count()
        )
        & 1
    )


def all_pairs_closure(
    n: int, edges: Sequence[tuple[int, int]]
) -> frozenset[int]:
    """Saturate new labels against every label already reached."""
    reached = set(raw_generators(n, edges))
    queue = deque(reached)
    while queue:
        left = queue.popleft()
        for right in tuple(reached):
            if bracket_nonzero(left, right, n):
                child = left ^ right
                if child not in reached:
                    reached.add(child)
                    queue.append(child)
    return frozenset(reached)


def closure_digest(values: Iterable[int], n: int) -> str:
    width = max(1, (2 * n + 7) // 8)
    digest = hashlib.sha256()
    for value in sorted(values):
        digest.update(value.to_bytes(width, "little"))
    return digest.hexdigest()


def q_value(label: int, n: int, colour: Sequence[int]) -> int:
    low = (1 << n) - 1
    x_bits, z_bits = label & low, label >> n
    return (
        (x_bits & z_bits).bit_count()
        + x_bits.bit_count()
        + sum(colour[v] for v in range(n) if (z_bits >> v) & 1)
    ) & 1


def q_roots(n: int, colour: Sequence[int]) -> frozenset[int]:
    global_x = (1 << n) - 1
    return frozenset(
        x_bits | (z_bits << n)
        for z_bits in range(1 << n)
        if not z_bits.bit_count() & 1
        for x_bits in range(1 << n)
        if (x_bits | (z_bits << n)) != global_x
        and q_value(x_bits | (z_bits << n), n, colour)
    )


def endpoint_masks(n: int, edges: Sequence[tuple[int, int]]) -> tuple[int, ...]:
    rows = adjacency(n, edges)
    endpoints = [0] * (1 << n)
    for vertex in range(n):
        endpoints[1 << vertex] = 1 << vertex
    for used in range(1 << n):
        live = endpoints[used]
        while live:
            end_bit = live & -live
            live ^= end_bit
            end = end_bit.bit_length() - 1
            choices = rows[end] & ~used
            while choices:
                choice = choices & -choices
                choices ^= choice
                endpoints[used | choice] |= choice
    return tuple(endpoints)


def path_cover_number(n: int, edges: Sequence[tuple[int, int]]) -> tuple[bool, int]:
    endpoints = endpoint_masks(n, edges)
    best = [n + 1] * (1 << n)
    best[0] = 0
    for target in range(1, 1 << n):
        anchor = target & -target
        subset = target
        while subset:
            if subset & anchor and endpoints[subset]:
                best[target] = min(best[target], 1 + best[target ^ subset])
            subset = (subset - 1) & target
    return bool(endpoints[-1]), best[-1]


def matching_number(n: int, edges: Sequence[tuple[int, int]]) -> int:
    rows = adjacency(n, edges)

    @lru_cache(maxsize=None)
    def visit(vertices: int) -> int:
        if vertices == 0:
            return 0
        first = vertices & -vertices
        vertex = first.bit_length() - 1
        rest = vertices ^ first
        answer = visit(rest)
        choices = rows[vertex] & rest
        while choices:
            partner = choices & -choices
            choices ^= partner
            answer = max(answer, 1 + visit(rest ^ partner))
        return answer

    return visit((1 << n) - 1)


def direct_majorana_grade(label: int, n: int, order: Sequence[int]) -> int:
    """Invert Jordan--Wigner coefficients by a backward suffix recurrence."""
    low = (1 << n) - 1
    x_bits, z_bits = label & low, label >> n
    suffix_z = 0
    grade = 0
    for vertex in reversed(order):
        z_bit = (z_bits >> vertex) & 1
        odd_coefficient = ((x_bits >> vertex) & 1) ^ suffix_z
        even_coefficient = z_bit ^ odd_coefficient
        grade += odd_coefficient + even_coefficient
        suffix_z ^= z_bit
    return grade


def grade_profile(
    closure: frozenset[int], n: int
) -> Counter[tuple[tuple[int, int], ...]]:
    profiles: Counter[tuple[tuple[int, int], ...]] = Counter()
    for order in vertex_permutations(n):
        histogram = Counter(direct_majorana_grade(label, n, order) for label in closure)
        profiles[tuple(sorted(histogram.items()))] += 1
    return profiles


def rank_gf2(values: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for value in values:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def theorem_dimension(n: int, part_sizes: Sequence[int], branch: str) -> int:
    if branch == "path":
        return n * (2 * n - 1)
    if branch == "even_cycle":
        return 2 * n * (2 * n - 1)
    smaller = min(part_sizes)
    if n & 1:
        return (1 << (2 * n - 2)) - 1
    return (1 << (2 * n - 2)) - ((-1) ** smaller) * (1 << (n - 1))


def stored_profile_counter(row: dict[str, object]) -> Counter[tuple[tuple[int, int], ...]]:
    answer: Counter[tuple[tuple[int, int], ...]] = Counter()
    for profile in row["profiles"]:
        histogram = tuple(
            sorted((int(grade), int(count)) for grade, count in profile["histogram"].items())
        )
        answer[histogram] = int(profile["order_count"])
    return answer


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    ok("artifact top-level shape", set(artifact) == {"meta", "data", "checks"})
    ok(
        "stored top-level checks pass",
        bool(artifact["checks"])
        and all(check["passed"] is True for check in artifact["checks"]),
    )

    census = artifact["data"]["census"]
    max_n = int(census["scope"]["maximum_n"])
    labeled, representatives = enumerate_unlabeled(max_n)
    stored_rows = census["graphs"]
    stored_by_signature = {
        (
            int(row["n"]),
            isomorphism_signature(
                int(row["n"]), tuple(tuple(edge) for edge in row["edges"])
            ),
        ): row
        for row in stored_rows
    }
    ok(
        "independent unlabeled census matches",
        set(representatives) == set(stored_by_signature)
        and {str(n): count for n, count in labeled.items()}
        == census["counts"]["connected_bipartite_labeled_by_n"],
        f"{len(representatives)} representatives",
    )

    derived: dict[tuple[int, tuple[int, ...]], dict[str, object]] = {}
    closure_failures: list[str] = []
    for key in representatives:
        n, _ = key
        stored = stored_by_signature[key]
        edges = tuple(tuple(edge) for edge in stored["edges"])
        colouring = graph_colouring(n, edges)
        if colouring is None or not colouring[0]:
            closure_failures.append(f"colouring {key}")
            continue
        colour = colouring[1]
        parts = sorted((colour.count(0), colour.count(1)))
        closure = all_pairs_closure(n, edges)
        roots = q_roots(n, colour)
        hamiltonian, path_cover = path_cover_number(n, edges)
        matching = matching_number(n, edges)
        max_degree = max(mask.bit_count() for mask in adjacency(n, edges))
        branch = (
            "branching_quadratic_root"
            if max_degree >= 3
            else ("path" if len(edges) == n - 1 else "even_cycle")
        )
        formula = theorem_dimension(n, parts, branch)
        signature_data = {
            "closure": closure,
            "dimension": len(closure),
            "parts": parts,
            "hamiltonian": hamiltonian,
            "path_cover": path_cover,
            "matching": matching,
            "max_degree": max_degree,
            "branch": branch,
        }
        derived[key] = signature_data
        if not (
            len(closure) == stored["dimension"] == stored["predicted_dimension"] == formula
            and closure_digest(closure, n) == stored["closure_sha256"]
            and parts == stored["bipartition_sizes"]
            and hamiltonian == stored["hamiltonian_path_exists"]
            and path_cover == stored["minimum_path_cover"]
            and matching == stored["maximum_matching_size"]
            and n - 2 * matching == stored["matching_deficiency"]
            and branch == stored["branch"]
            and (max_degree < 3 or closure == roots)
        ):
            closure_failures.append(f"row n={n}, signature={key[1]}")
        budget_tick(f"independent closure n={n}")
    ok(
        "all exact dimensions and invariants independently reproduce",
        not closure_failures,
        f"{len(derived)} all-pairs closures",
    )

    stored_grades = artifact["data"]["candidate_invariants"]["grade_profiles"]
    grade_failures: list[str] = []
    for profile_row in stored_grades:
        n = int(profile_row["n"])
        edges = tuple(tuple(edge) for edge in profile_row["edges"])
        key = (n, isomorphism_signature(n, edges))
        closure = derived[key]["closure"]
        measured = grade_profile(closure, n)
        expected = stored_profile_counter(profile_row)
        if measured != expected or sum(measured.values()) != math.factorial(n):
            grade_failures.append(f"n={n}, code={profile_row['canonical_code_hex']}")
        for stored_profile in profile_row["profiles"]:
            order = tuple(stored_profile["example_order"])
            histogram = Counter(
                direct_majorana_grade(label, n, order) for label in closure
            )
            if {str(k): v for k, v in sorted(histogram.items())} != stored_profile[
                "histogram"
            ]:
                grade_failures.append(
                    f"example n={n}, code={profile_row['canonical_code_hex']}"
                )
        budget_tick(f"independent grades n={n}")
    ok(
        "all vertex-order grade profiles independently reproduce",
        not grade_failures,
        f"{len(stored_grades)} non-Hamiltonian graphs",
    )

    certificate_failures: list[str] = []
    for certificate in artifact["data"]["spanning_tree_e6_certificates"]:
        n = int(certificate["n"])
        tree = tuple(tuple(edge) for edge in certificate["spanning_tree_edges"])
        basis = raw_generators(n, tree)
        labels = tuple(int(label) for label in certificate["e6_labels"])
        pairing_edges = [
            [left, right]
            for left in range(6)
            for right in range(left + 1, 6)
            if bracket_nonzero(labels[left], labels[right], n)
        ]
        if not (
            len(tree) == n - 1
            and rank_gf2(basis) == 2 * n - 1
            and pairing_edges == certificate["expected_e6_pairing_edges"]
            and certificate["e6_induced"] is True
        ):
            certificate_failures.append(
                f"n={n}, code={certificate['canonical_code']}"
            )
    ok(
        "spanning-tree bases and induced E6 certificates reproduce",
        not certificate_failures,
        f"{len(artifact['data']['spanning_tree_e6_certificates'])} certificates",
    )

    counterexamples = artifact["data"]["candidate_invariants"]["counterexamples"]
    matching_pair = counterexamples["matching_deficiency_and_path_cover"]
    grade_pair = counterexamples["same_grade_set_different_dimension"]
    order_pair = counterexamples["order_dependent_grades"]
    colour_pair = counterexamples["colour_imbalance_magnitude"]
    ok(
        "candidate-invariant counterexamples are internally decisive",
        matching_pair["left"]["matching_deficiency"]
        == matching_pair["right"]["matching_deficiency"]
        and matching_pair["left"]["minimum_path_cover"]
        == matching_pair["right"]["minimum_path_cover"]
        and matching_pair["left"]["dimension"]
        != matching_pair["right"]["dimension"]
        and grade_pair["left"]["profile"]["grades"]
        == grade_pair["right"]["profile"]["grades"]
        and grade_pair["left"]["dimension"] != grade_pair["right"]["dimension"]
        and order_pair["first_profile"]["grades"]
        != order_pair["second_profile"]["grades"]
        and colour_pair["left"]["smaller_colour_parity"]
        == colour_pair["right"]["smaller_colour_parity"]
        and colour_pair["left"]["colour_imbalance"]
        != colour_pair["right"]["colour_imbalance"]
        and colour_pair["left"]["dimension"] == colour_pair["right"]["dimension"],
    )

    complete_controls = artifact["data"]["complete_bipartite_finite_controls"]
    complete_ok = True
    for row in complete_controls:
        left, right = row["part_sizes"]
        n = left + right
        edges = tuple(
            (u, left + v) for u in range(left) for v in range(right)
        )
        key = (n, isomorphism_signature(n, edges))
        complete_ok &= (
            derived[key]["hamiltonian"] is False
            and derived[key]["dimension"] == row["formula_dimension"]
            and closure_digest(derived[key]["closure"], n) == row["closure_sha256"]
        )
    ok(
        "complete-bipartite controls independently reproduce",
        complete_ok and bool(complete_controls),
        f"{len(complete_controls)} controls",
    )

    hash_failures = []
    for relative, expected in artifact["meta"]["source_sha256"].items():
        measured = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if measured != expected:
            hash_failures.append(relative)
    ok("artifact source hashes match", not hash_failures)
    ok(
        "scope exclusions remain explicit",
        artifact["data"]["limits"]["tag"] == "[UNRESOLVED]"
        and "3D" in artifact["data"]["verdict"]["ising_scope"]
        and "two-generator" in artifact["data"]["limits"]["excluded"],
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
