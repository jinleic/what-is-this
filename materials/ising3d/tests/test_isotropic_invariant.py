"""Clean-room verifier for the isotropic trace-invariant certificate.

This script imports none of e200--e202.  It independently rebuilds alignment
polynomials from spin tuples, Eulerian counts from edge-choice tuples, the
characteristic-coefficient invariant, physical-curve clearings, and every stored
finite audit row.  The finite rows test lemmas; they are not used to infer the
all-size theorem.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import platform
import resource
import time
from collections import deque
from fractions import Fraction
from math import comb, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "isotropic_invariant.json"
CPU_LIMIT_SECONDS = 30.0
RSS_LIMIT_BYTES = 2_000_000_000
STARTED = time.process_time()
FAILURES: list[str] = []
CHECK_COUNT = 0

GRAPHS: dict[str, tuple[int, tuple[tuple[int, int], ...]]] = {
    "path_4": (4, ((0, 1), (1, 2), (2, 3))),
    "square": (4, ((0, 1), (1, 2), (2, 3), (3, 0))),
    "triangle": (3, ((0, 1), (1, 2), (2, 0))),
    "paw": (4, ((0, 1), (1, 2), (2, 0), (0, 3))),
    "lollipop_2": (5, ((0, 1), (1, 2), (2, 0), (0, 3), (3, 4))),
    "cycle_5": (5, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))),
    "complete_4": (4, tuple((u, v) for u in range(4) for v in range(u + 1, 4))),
}


def rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def budget(stage: str) -> None:
    used = time.process_time() - STARTED
    if used >= CPU_LIMIT_SECONDS:
        raise RuntimeError(f"process-time budget exceeded at {stage}: {used:.3f}s")
    if rss_bytes() >= RSS_LIMIT_BYTES:
        raise RuntimeError(f"RSS cap exceeded at {stage}: {rss_bytes()}")


def check(name: str, passed: bool, detail: str = "") -> None:
    global CHECK_COUNT
    CHECK_COUNT += 1
    print(f"[{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    if not passed:
        FAILURES.append(name)


def clean(poly: list[int]) -> list[int]:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def plus(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return clean(answer)


def times(left: list[int], right: list[int]) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return clean(answer)


def repeated_power(poly: list[int], exponent: int) -> list[int]:
    """Deliberately use repeated multiplication, unlike the producer's binary power."""
    answer = [1]
    for _ in range(exponent):
        answer = times(answer, poly)
    return answer


def evaluate(poly: list[int], value: Fraction) -> Fraction:
    return sum((Fraction(coefficient) * value**degree for degree, coefficient in enumerate(poly)), Fraction(0))


def digest(poly: list[int]) -> str:
    return hashlib.sha256(",".join(str(value) for value in clean(poly)).encode()).hexdigest()


def alignment_by_spin_tuples(sites: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    counts = [0] * (len(edges) + 1)
    for spins in itertools.product((-1, 1), repeat=sites):
        aligned = sum(spins[left] == spins[right] for left, right in edges)
        counts[aligned] += 1
    return counts


def eulerian_by_edge_tuples(sites: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    counts = [0] * (len(edges) + 1)
    for chosen in itertools.product((0, 1), repeat=len(edges)):
        degrees = [0] * sites
        for include, (left, right) in zip(chosen, edges):
            if include:
                degrees[left] += 1
                degrees[right] += 1
        if all(degree % 2 == 0 for degree in degrees):
            counts[sum(chosen)] += 1
    return counts


def bipartite_coloring(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    color = [-1] * sites
    for root in range(sites):
        if color[root] >= 0:
            continue
        color[root] = 0
        queue = deque([root])
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if color[neighbor] < 0:
                    color[neighbor] = color[vertex] ^ 1
                    queue.append(neighbor)
                elif color[neighbor] == color[vertex]:
                    return False
    return True


def component_count(sites: int, edges: tuple[tuple[int, int], ...]) -> int:
    parent = list(range(sites))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right in edges:
        left_root, right_root = root(left), root(right)
        if left_root != right_root:
            parent[right_root] = left_root
    return len({root(vertex) for vertex in range(sites)})


def shifted_at_one(poly: list[int]) -> list[int]:
    answer = [0] * len(poly)
    for old_degree, coefficient in enumerate(poly):
        for new_degree in range(old_degree + 1):
            answer[new_degree] += coefficient * comb(old_degree, new_degree)
    return clean(answer)


def clear_physical_denominator(poly: list[int], power: int) -> list[int]:
    answer = [0] * (2 * power + 1)
    for y_degree, coefficient in enumerate(poly):
        for chosen_t_squared in range(y_degree + 1):
            degree = 2 * chosen_t_squared + power - y_degree
            answer[degree] += (
                coefficient
                * comb(y_degree, chosen_t_squared)
                * 2 ** (power - y_degree)
            )
    return clean(answer)


def fraction_matches(record: dict[str, object], value: Fraction) -> bool:
    return int(record["numerator"]) == value.numerator and int(record["denominator"]) == value.denominator


def verify_artifact_envelope(artifact: dict[str, object]) -> None:
    check("artifact has exactly meta/data/checks", set(artifact) == {"meta", "data", "checks"})
    checks = artifact["checks"]
    check(
        "all producer checks are stored passing",
        isinstance(checks, list) and len(checks) > 0 and all(bool(item["passed"]) for item in checks),
        f"stored={len(checks) if isinstance(checks, list) else 'invalid'}",
    )
    source_hashes = artifact["meta"]["source_sha256"]
    hashes_match = True
    for relative, expected in source_hashes.items():
        hashes_match &= hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected
    check("producer source hashes match artifact", hashes_match, f"files={len(source_hashes)}")
    check(
        "artifact declares exact non-modular arithmetic",
        "no floating point or modular rank" in str(artifact["meta"]["arithmetic"]),
    )


def verify_graph_expansions(artifact: dict[str, object]) -> None:
    stored_rows = {
        row["label"]: row
        for row in artifact["data"]["high_temperature_leading_terms"]["graph_rows"]
    }
    check("high-temperature graph labels match verifier fixtures", set(stored_rows) == set(GRAPHS))
    y0 = Fraction(5, 3)
    t0 = Fraction(1, 3)
    for label, (sites, edges) in GRAPHS.items():
        row = stored_rows[label]
        q = alignment_by_spin_tuples(sites, edges)
        q_star = list(reversed(q))
        skew = clean([left - right for left, right in zip(q, q_star)])
        eulerian = eulerian_by_edge_tuples(sites, edges)
        odd_sizes = [size for size, count in enumerate(eulerian) if size % 2 and count]
        odd_girth = min(odd_sizes) if odd_sizes else None
        shortest_count = 0 if odd_girth is None else eulerian[odd_girth]
        shifted = shifted_at_one(skew)
        nonzero_shifted = [(degree, coefficient) for degree, coefficient in enumerate(shifted) if coefficient]
        valuation = nonzero_shifted[0][0] if nonzero_shifted else None
        leading = nonzero_shifted[0][1] if nonzero_shifted else 0
        expected_leading = 0 if odd_girth is None else 2 ** (sites + 1 - odd_girth) * shortest_count
        physical = clear_physical_denominator(skew, len(edges))
        value_y0 = evaluate(skew, y0)
        row_matches = (
            row["sites"] == sites
            and row["edges"] == [list(edge) for edge in edges]
            and row["alignment_coefficients_ascending"] == q
            and row["reverse_alignment_coefficients_ascending"] == q_star
            and row["trace_skew_coefficients_ascending"] == skew
            and row["eulerian_counts_by_cardinality"] == eulerian
            and row["connected_components"] == component_count(sites, edges)
            and row["bipartite"] == bipartite_coloring(sites, edges)
            and row["odd_girth"] == odd_girth
            and row["shortest_odd_cycle_count"] == shortest_count
            and row["y_minus_one_valuation"] == valuation
            and row["leading_coefficient_at_y_one"] == leading == expected_leading
            and row["physical_cleared_coefficients_ascending"] == physical
            and fraction_matches(row["value_at_y_5_over_3"], value_y0)
            and evaluate(physical, t0) == (2 * t0) ** len(edges) * value_y0
            and sum(q) == 1 << sites
            and sum(eulerian) == 1 << (len(edges) - sites + component_count(sites, edges))
        )
        check(f"{label}: independent spin/Eulerian/Puiseux row", row_matches)
        budget(f"graph expansion {label}")


def verify_subset_invariant(artifact: dict[str, object]) -> None:
    stored = artifact["data"]["spectral_invariant"]["characteristic_coefficient_invariant"]
    audit_rows = stored["exact_slot_audits"]
    all_zero = True
    for modes in range(1, 5):
        base = Fraction(5, 7)
        factors = [Fraction(2 * index + 3, index + 2) for index in range(modes)]
        slots = []
        for bits in itertools.product((0, 1), repeat=modes):
            slots.append(base * prod((factor for bit, factor in zip(bits, factors) if bit), start=Fraction(1)))
        determinant = prod(slots, start=Fraction(1))
        e1 = sum(slots, start=Fraction(0))
        e_last = sum(
            (prod((slots[j] for j in range(len(slots)) if j != omitted), start=Fraction(1)) for omitted in range(len(slots))),
            start=Fraction(0),
        )
        half = len(slots) // 2
        phi = determinant ** (half - 1) * e1**half - e_last**half
        all_zero &= phi == 0
        artifact_row = audit_rows[modes - 1]
        all_zero &= artifact_row["modes"] == modes and artifact_row["slots"] == 1 << modes
    check("independent subset slots annihilate Phi_N for n=1..4", all_zero)

    sample = [Fraction(2), Fraction(3), Fraction(1, 5), Fraction(5, 6)]
    determinant = prod(sample, start=Fraction(1))
    e1 = sum(sample, start=Fraction(0))
    e3 = sum(
        (prod((sample[j] for j in range(4) if j != omitted), start=Fraction(1)) for omitted in range(4)),
        start=Fraction(0),
    )
    check("independent positive determinant-one sample has nonzero Phi_4", determinant == 1 and e1 * e1 != e3 * e3)


def verify_spectral_rows(artifact: dict[str, object]) -> None:
    rows = {
        row["label"]: row
        for row in artifact["data"]["spectral_invariant"]["physical_graph_reduction"]["graph_rows"]
    }
    check("spectral graph labels match verifier fixtures", set(rows) == set(GRAPHS))
    y0 = Fraction(5, 3)
    t0 = Fraction(1, 3)
    for label, (sites, edges) in GRAPHS.items():
        q = alignment_by_spin_tuples(sites, edges)
        q_star = list(reversed(q))
        half = 1 << (sites - 1)
        left_power = repeated_power(q, half)
        right_power = repeated_power(q_star, half)
        psi = plus(left_power, [-value for value in right_power])
        zero = psi == [0]
        expected_degree = len(edges) * half
        cleared = [0] if zero else clear_physical_denominator(psi, expected_degree)
        value = evaluate(psi, y0)
        row = rows[label]
        row_matches = (
            row["dimension"] == 1 << sites
            and row["half_dimension_M"] == half
            and row["psi_zero"] == zero
            and row["psi_degree"] == (None if zero else expected_degree)
            and row["psi_sha256"] == digest(psi)
            and row["psi_nonzero_coefficient_count"] == sum(coefficient != 0 for coefficient in psi)
            and row["psi_constant_coefficient"] == psi[0]
            and row["psi_leading_coefficient"] == psi[-1]
            and fraction_matches(row["psi_at_y_5_over_3"], value)
            and row["physical_cleared_sha256"] == digest(cleared)
            and row["physical_cleared_degree"] == (None if zero else 2 * expected_degree)
            and row["complex_exception_root_bound"] == (None if zero else 2 * expected_degree)
            and row["nonzero_modulo_physical_curve"] == (not zero)
            and (zero or evaluate(cleared, t0) == (2 * t0) ** expected_degree * value)
            and ((value == 0) == zero)
        )
        check(f"{label}: independent determinant-center-free invariant row", row_matches)
        budget(f"spectral row {label}")


def lollipop_edges(tail: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    edges = [(0, 1), (1, 2), (2, 0)]
    endpoint = 0
    for index in range(tail):
        new_vertex = index + 3
        edges.append((endpoint, new_vertex))
        endpoint = new_vertex
    return tail + 3, tuple(edges)


def verify_leaf_family(artifact: dict[str, object]) -> None:
    rows = artifact["data"]["spectral_invariant"]["leaf_extension"]["finite_lemma_tests"]
    rows_by_tail = {int(row["tail_edges_r"]): row for row in rows}
    check("stored lollipop audit range is r=0..5", set(rows_by_tail) == set(range(6)))
    triangle_skew = [-2, 6, -6, 2]
    for tail in range(6):
        sites, edges = lollipop_edges(tail)
        q = alignment_by_spin_tuples(sites, edges)
        skew = clean([q[index] - q[-1 - index] for index in range(len(q))])
        expected = times(triangle_skew, repeated_power([1, 1], tail))
        physical = clear_physical_denominator(skew, len(edges))
        expected_physical = times(
            [2 * comb(6, degree) * (-1) ** (6 - degree) for degree in range(7)],
            [comb(2 * tail, degree) for degree in range(2 * tail + 1)],
        )
        row = rows_by_tail[tail]
        check(
            f"lollipop r={tail}: independent leaf and physical recurrences",
            skew == expected
            and physical == expected_physical
            and row["sites"] == sites
            and row["edge_count"] == len(edges)
            and row["checks_passed"]
            and row["trace_skew_sha256"] == digest(skew)
            and row["physical_trace_skew_sha256"] == digest(physical),
        )


def verify_curve_and_scope(artifact: dict[str, object]) -> None:
    t = Fraction(1, 3)
    x = (1 + t) / (1 - t)
    y = (1 + t * t) / (2 * t)
    check("physical witness point is exact", x == 2 and y == Fraction(5, 3) and (y - 1) * x * x == y + 1)
    conclusion = artifact["data"]["conclusion"]
    boundary = artifact["data"]["route_boundary"]
    check(
        "artifact separates theorem from unresolved bipartite boundary",
        conclusion["tag"] == "[THEOREM]"
        and boundary["tag"] == "[UNRESOLVED]"
        and "non-bipartite" in conclusion["headline"]
        and "Bipartite" in boundary["statement"],
    )


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    verify_artifact_envelope(artifact)
    verify_graph_expansions(artifact)
    verify_subset_invariant(artifact)
    verify_spectral_rows(artifact)
    verify_leaf_family(artifact)
    verify_curve_and_scope(artifact)
    budget("final")
    check(
        "verifier process-time and RSS budgets",
        time.process_time() - STARTED < CPU_LIMIT_SECONDS and rss_bytes() < RSS_LIMIT_BYTES,
        f"CPU={time.process_time()-STARTED:.3f}s, RSS={rss_bytes()}",
    )
    if FAILURES:
        print(f"FAIL test_isotropic_invariant ({len(FAILURES)}/{CHECK_COUNT} failed)")
        return 1
    print(f"PASS test_isotropic_invariant ({CHECK_COUNT} checks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
