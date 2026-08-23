"""Exact graph-level audit of the Clifford-grade trichotomy.

Small and medium graphs are closed directly from the raw local Pauli generators.
The 12-site graph is not materialised: it is certified by the path/chord grades,
the corrected overlap lemma, and an explicit grade-2 orbit walk.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections import Counter
from math import comb
from typing import Iterable, Sequence

from e178_clifford_setup import (
    bipartite_coloring,
    budget_tick,
    cycle_edges,
    exact_pauli_closure,
    full_even_derived_dimension,
    grade_two_mod_four_dimension,
    grid_edges,
    max_rss_bytes,
    nonpath_edge_rows,
    normalize_edges,
    snake_path,
    symplectic_form,
    validate_hamiltonian_path,
)
from e179_grade_generation import (
    corrected_generation_witness,
    grade_closure,
    orbit_swap_witness,
    target_noncentral_grades,
)

CPU_BUDGET_SECONDS = 240.0


def path_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, site + 1) for site in range(n - 1))


def max_degree(n: int, edges: Sequence[tuple[int, int]]) -> int:
    degrees = [0] * n
    for u, v in normalize_edges(edges, n):
        degrees[u] += 1
        degrees[v] += 1
    return max(degrees, default=0)


def majorana_masks(n: int, path: Sequence[int]) -> tuple[int, ...]:
    values: list[int] = []
    prefix_x = 0
    for site in path:
        z = 1 << (n + site)
        values.append(prefix_x | z)
        values.append(prefix_x | (1 << site) | z)
        prefix_x |= 1 << site
    return tuple(values)


def majorana_echelon(
    n: int, path: Sequence[int]
) -> dict[int, tuple[int, int]]:
    echelon: dict[int, tuple[int, int]] = {}
    for index, vector in enumerate(majorana_masks(n, path)):
        coefficients = 1 << index
        value = vector
        while value:
            pivot = value.bit_length() - 1
            if pivot not in echelon:
                echelon[pivot] = (value, coefficients)
                break
            row, row_coefficients = echelon[pivot]
            value ^= row
            coefficients ^= row_coefficients
        if not value:
            raise AssertionError("Majorana masks are not independent")
    return echelon


def pauli_to_majorana_support(
    pauli: int, echelon: dict[int, tuple[int, int]]
) -> int:
    """Invert the exact GF(2) Jordan--Wigner basis, ignoring only scalar phase."""
    value = int(pauli)
    coefficients = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in echelon:
            raise AssertionError("Pauli vector is outside the Majorana span")
        row, row_coefficients = echelon[pivot]
        value ^= row
        coefficients ^= row_coefficients
    return coefficients


def grade_histogram(
    elements: Iterable[int], n: int, path: Sequence[int]
) -> dict[int, int]:
    echelon = majorana_echelon(n, path)
    histogram = Counter(
        pauli_to_majorana_support(pauli, echelon).bit_count() for pauli in elements
    )
    return dict(sorted(histogram.items()))


def closure_digest(elements: Iterable[int]) -> str:
    digest = hashlib.sha256()
    for value in sorted(elements):
        digest.update(f"{value:x}\n".encode("ascii"))
    return digest.hexdigest()


def pairwise_closed(elements: frozenset[int], n: int) -> bool:
    values = tuple(elements)
    for index, left in enumerate(values):
        for right in values[index + 1 :]:
            if symplectic_form(left, right, n) and (left ^ right) not in elements:
                return False
    return True


def expected_histogram(branch: str, n: int) -> dict[int, int]:
    N = 2 * n
    if branch == "path":
        return {2: comb(N, 2)}
    if branch == "even_cycle":
        return {2: comb(N, 2), N - 2: comb(N, N - 2)}
    if branch == "branching_bipartite":
        return {grade: comb(N, grade) for grade in target_noncentral_grades(N)}
    if branch == "nonbipartite_full_even_derived":
        return {grade: comb(N, grade) for grade in range(2, N, 2)}
    raise ValueError(branch)


def branch_dimension(branch: str, n: int) -> int:
    N = 2 * n
    if branch == "path":
        return n * (2 * n - 1)
    if branch == "even_cycle":
        return N * (N - 1)
    if branch == "branching_bipartite":
        return grade_two_mod_four_dimension(n)
    if branch == "nonbipartite_full_even_derived":
        return full_even_derived_dimension(n)
    raise ValueError(branch)


def exact_graph_row(
    label: str,
    n: int,
    edges: Sequence[tuple[int, int]],
    path: Sequence[int],
    branch: str,
    *,
    pairwise_audit: bool = False,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    validate_hamiltonian_path(n, edges, path)
    bipartite, _ = bipartite_coloring(n, edges)
    chords = nonpath_edge_rows(n, edges, path)
    closure = exact_pauli_closure(n, edges, budget=CPU_BUDGET_SECONDS)
    histogram = grade_histogram(closure.elements, n, path)
    expected_hist = expected_histogram(branch, n)
    expected_dimension = branch_dimension(branch, n)
    checks = [
        {
            "name": f"{label}: exact raw-generator dimension",
            "passed": closure.dimension == expected_dimension,
            "detail": f"measured={closure.dimension}, expected={expected_dimension}",
        },
        {
            "name": f"{label}: exact Clifford-grade histogram",
            "passed": histogram == expected_hist,
            "detail": f"measured={histogram}, expected={expected_hist}",
        },
    ]
    pairwise_passed: bool | None = None
    if pairwise_audit:
        pairwise_passed = pairwise_closed(closure.elements, n)
        checks.append(
            {
                "name": f"{label}: pairwise bracket saturation",
                "passed": pairwise_passed,
                "detail": f"all pairs checked in {closure.dimension}-element closure",
            }
        )
    row = {
        "tag": "[COMPUTATION]",
        "label": label,
        "n": n,
        "N_majoranas": 2 * n,
        "edges": [list(edge) for edge in normalize_edges(edges, n)],
        "hamiltonian_path": list(path),
        "bipartite": bipartite,
        "maximum_degree": max_degree(n, edges),
        "branch": branch,
        "nonpath_edges": list(chords),
        "method": "exact raw-local-generator Pauli BFS; grade histogram via independent GF(2) JW inversion",
        "dimension": closure.dimension,
        "expected_dimension": expected_dimension,
        "grade_histogram": {str(key): value for key, value in histogram.items()},
        "closure_sha256": closure_digest(closure.elements),
        "pairwise_saturation_audited": pairwise_passed,
        "process_time_seconds": round(closure.process_time_seconds, 6),
        "peak_rss_bytes": closure.peak_rss_bytes,
    }
    return row, checks


def characterization_row(
    label: str,
    n: int,
    edges: Sequence[tuple[int, int]],
    path: Sequence[int],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    validate_hamiltonian_path(n, edges, path)
    bipartite, _ = bipartite_coloring(n, edges)
    chords = nonpath_edge_rows(n, edges, path)
    N = 2 * n
    usable = [row for row in chords if 6 <= int(row["majorana_grade"]) <= N - 4]
    if not usable:
        raise AssertionError(f"{label} has no non-endpoint chord")
    seed_row = min(usable, key=lambda row: int(row["majorana_grade"]))
    seed_grade = int(seed_row["majorana_grade"])
    seed_support = tuple(int(value) for value in seed_row["support_zero_based"])
    target_support = tuple(range(N - seed_grade, N))
    if set(target_support) == set(seed_support):
        target_support = tuple(range(seed_grade))
    orbit = orbit_swap_witness(N, seed_support, target_support)
    generation = corrected_generation_witness(N, seed_grade)
    target = target_noncentral_grades(N)
    closure = grade_closure(N, (2, seed_grade))
    expected_dimension = grade_two_mod_four_dimension(n)
    checks = [
        {
            "name": f"{label}: bipartite Hamiltonian input",
            "passed": bipartite and bool(chords),
            "detail": f"{len(chords)} chords, Delta={max_degree(n, edges)}",
        },
        {
            "name": f"{label}: every chord has grade 2 mod 4",
            "passed": all(int(row["grade_mod_4"]) == 2 for row in chords),
            "detail": f"grades={[row['majorana_grade'] for row in chords]}",
        },
        {
            "name": f"{label}: corrected grade closure",
            "passed": closure == target,
            "detail": f"seed={seed_grade}, closure={closure}",
        },
        {
            "name": f"{label}: certified coordinate-orbit sample",
            "passed": bool(orbit["passed"]) and bool(orbit["steps"]),
            "detail": f"{len(orbit['steps'])} exact grade-2 swaps",
        },
    ]
    row = {
        "tag": "[THEOREM]",
        "label": label,
        "n": n,
        "N_majoranas": N,
        "edges": [list(edge) for edge in normalize_edges(edges, n)],
        "hamiltonian_path": list(path),
        "bipartite": bipartite,
        "maximum_degree": max_degree(n, edges),
        "branch": "branching_bipartite",
        "nonpath_edges": list(chords),
        "method": (
            "non-materialized exact grade characterization: path gives grade 2; one non-endpoint "
            "chord gives the displayed seed; overlap witnesses generate every allowed grade; "
            "the displayed Johnson walk certifies a sample coordinate orbit"
        ),
        "seed_chord": seed_row,
        "generation_witness": generation,
        "orbit_sample": orbit,
        "grade_set": list(target),
        "dimension": expected_dimension,
        "dimension_sum": sum(comb(N, grade) for grade in target),
        "materialized_closure": False,
    }
    return row, checks


def run_graph_audit() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    exact_rows: list[dict[str, object]] = []

    exact_specs: list[
        tuple[str, int, tuple[tuple[int, int], ...], tuple[int, ...], str, bool]
    ] = []
    exact_specs.append(("path_P6_control", 6, path_edges(6), tuple(range(6)), "path", True))
    exact_specs.append(("cycle_C6_counterexample", 6, cycle_edges(6), tuple(range(6)), "even_cycle", True))
    exact_specs.append(("cycle_C8_counterexample", 8, cycle_edges(8), tuple(range(8)), "even_cycle", True))

    for rows, cols, expected_label in (
        (2, 2, "known_grid_2x2_C4"),
        (2, 3, "known_grid_2x3"),
        (2, 4, "known_grid_2x4"),
        (2, 5, "known_grid_2x5"),
        (3, 3, "known_grid_3x3"),
    ):
        n = rows * cols
        branch = "even_cycle" if (rows, cols) == (2, 2) else "branching_bipartite"
        exact_specs.append(
            (expected_label, n, grid_edges(rows, cols), snake_path(rows, cols), branch, False)
        )

    removed_edge = (0, 3)
    grid_3x3_minus = tuple(edge for edge in grid_edges(3, 3) if edge != removed_edge)
    exact_specs.append(
        (
            "new_grid_3x3_minus_edge_0_3",
            9,
            grid_3x3_minus,
            snake_path(3, 3),
            "branching_bipartite",
            False,
        )
    )

    c5_plus_chord = normalize_edges((*cycle_edges(5), (0, 3)), 5)
    exact_specs.append(
        (
            "nonbipartite_C5_plus_chord_0_3",
            5,
            c5_plus_chord,
            tuple(range(5)),
            "nonbipartite_full_even_derived",
            True,
        )
    )

    for spec in exact_specs:
        row, row_checks = exact_graph_row(*spec[:-1], pairwise_audit=spec[-1])
        exact_rows.append(row)
        checks.extend(row_checks)
        budget_tick(started, f"exact graph {spec[0]}", CPU_BUDGET_SECONDS)

    characterization_rows: list[dict[str, object]] = []
    grid_3x4, grid_3x4_checks = characterization_row(
        "new_grid_3x4", 12, grid_edges(3, 4), snake_path(3, 4)
    )
    characterization_rows.append(grid_3x4)
    checks.extend(grid_3x4_checks)

    c10_chord_edges = normalize_edges((*cycle_edges(10), (0, 3)), 10)
    c10_chord, c10_checks = characterization_row(
        "new_non_grid_C10_plus_chord_0_3",
        10,
        c10_chord_edges,
        tuple(range(10)),
    )
    characterization_rows.append(c10_chord)
    checks.extend(c10_checks)

    by_label = {row["label"]: row for row in exact_rows}
    known_expected = {
        "known_grid_2x2_C4": 56,
        "known_grid_2x3": 1056,
        "known_grid_2x4": 16256,
        "known_grid_2x5": 262656,
        "known_grid_3x3": 65535,
    }
    checks.append(
        {
            "name": "five repository dimensions reproduced",
            "passed": all(by_label[label]["dimension"] == value for label, value in known_expected.items()),
            "detail": str({label: by_label[label]["dimension"] for label in known_expected}),
        }
    )
    checks.append(
        {
            "name": "C6 refutes the proposed universal dimension",
            "passed": (
                by_label["cycle_C6_counterexample"]["dimension"] == 132
                and grade_two_mod_four_dimension(6) == 1056
            ),
            "detail": "measured 132=2*C(12,2), proposed 1056",
        }
    )
    nonbip = by_label["nonbipartite_C5_plus_chord_0_3"]
    checks.append(
        {
            "name": "non-bipartite falsification reaches full even derived algebra",
            "passed": (
                nonbip["dimension"] == 510
                and grade_two_mod_four_dimension(5) == 255
                and full_even_derived_dimension(5) == 510
            ),
            "detail": "measured 510; 2-mod-4 prediction 255; full even derived 510",
        }
    )
    checks.append(
        {
            "name": "new n=12 prediction",
            "passed": grid_3x4["dimension"] == 4_192_256,
            "detail": f"3x4 grade characterization gives {grid_3x4['dimension']}",
        }
    )

    if not all(bool(row["passed"]) for row in checks):
        failed = [row for row in checks if not row["passed"]]
        raise AssertionError(f"graph audit failed: {failed}")
    return {
        "tag": "[COMPUTATION]",
        "exact_rows": exact_rows,
        "characterization_rows": characterization_rows,
        "checks": checks,
        "process_time_seconds": round(time.process_time() - started, 6),
        "peak_rss_bytes": max_rss_bytes(),
    }


def main() -> int:
    result = run_graph_audit()
    print(
        json.dumps(
            {
                "exact_rows": len(result["exact_rows"]),
                "characterization_rows": len(result["characterization_rows"]),
                "checks": len(result["checks"]),
                "process_time_seconds": result["process_time_seconds"],
                "peak_rss_bytes": result["peak_rss_bytes"],
            }
        )
    )
    print("PASS e180 graph verification")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
