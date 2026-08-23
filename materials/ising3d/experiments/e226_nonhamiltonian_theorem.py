"""Integrate the all-size non-Hamiltonian bipartite classification."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from e224_nonhamiltonian_census import (
    CPU_BUDGET_SECONDS,
    RSS_CAP_BYTES,
    adjacency_masks,
    bipartition,
    budget_tick,
    canonical_code,
    edge_code,
    max_rss_bytes,
    normalize_edges,
    quadratic_dimension,
    quadratic_value,
    raw_local_generators,
    run_census,
    symplectic_form,
)
from e225_matching_grade import run_matching_grade

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "nonhamiltonian_bipartite.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gf2_rank(values: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for original in values:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def branching_spanning_tree(
    n: int, edges: Sequence[tuple[int, int]]
) -> tuple[tuple[tuple[int, int], ...], int, tuple[int, int, int]]:
    normalized = normalize_edges(n, edges)
    adjacency = adjacency_masks(n, normalized)
    centre = next(
        vertex for vertex, mask in enumerate(adjacency) if mask.bit_count() >= 3
    )
    neighbours = tuple(
        vertex for vertex in range(n) if (adjacency[centre] >> vertex) & 1
    )[:3]
    selected = {
        (centre, neighbour) if centre < neighbour else (neighbour, centre)
        for neighbour in neighbours
    }

    parent = list(range(n))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> bool:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            return False
        parent[root_right] = root_left
        return True

    for left, right in sorted(selected):
        if not union(left, right):
            raise AssertionError("three incident edges cannot contain a cycle")
    for edge in normalized:
        if edge in selected:
            continue
        if union(*edge):
            selected.add(edge)
        if len(selected) == n - 1:
            break
    tree = tuple(sorted(selected))
    if len(tree) != n - 1 or len({find(vertex) for vertex in range(n)}) != 1:
        raise AssertionError(f"failed to extend branch star to a spanning tree: {tree}")
    return tree, centre, neighbours


def e6_certificate(row: dict[str, object]) -> dict[str, object]:
    n = int(row["n"])
    edges = tuple(tuple(edge) for edge in row["edges"])
    tree, centre, neighbours = branching_spanning_tree(n, edges)
    fields = tuple(1 << vertex for vertex in range(n))
    bonds = tuple(
        (1 << (n + left)) | (1 << (n + right)) for left, right in tree
    )
    basis = fields + bonds
    connected, is_bipartite, colours = bipartition(n, edges)
    if not connected or not is_bipartite:
        raise AssertionError((n, edges))

    incident_bonds = tuple(
        (1 << (n + centre)) | (1 << (n + neighbour))
        for neighbour in neighbours
    )
    labels = (1 << centre,) + incident_bonds + (
        1 << neighbours[0],
        1 << neighbours[1],
    )
    pairing_edges = tuple(
        (left, right)
        for left in range(len(labels))
        for right in range(left + 1, len(labels))
        if symplectic_form(labels[left], labels[right], n)
    )
    expected_pairing_edges = ((0, 1), (0, 2), (0, 3), (1, 4), (2, 5))
    return {
        "tag": "[COMPUTATION]",
        "n": n,
        "canonical_code": row["canonical_code"],
        "spanning_tree_edges": [list(edge) for edge in tree],
        "branch_centre": centre,
        "selected_neighbours": list(neighbours),
        "tree_generator_count": len(basis),
        "tree_generator_gf2_rank": gf2_rank(basis),
        "expected_basis_rank": 2 * n - 1,
        "all_tree_generators_have_Q_one": all(
            quadratic_value(label, n, colours) == 1 for label in basis
        ),
        "e6_labels": list(labels),
        "e6_pairing_edges": [list(edge) for edge in pairing_edges],
        "expected_e6_pairing_edges": [list(edge) for edge in expected_pairing_edges],
        "e6_induced": pairing_edges == expected_pairing_edges,
    }


def complete_bipartite_edges(
    left_size: int, right_size: int
) -> tuple[tuple[int, int], ...]:
    if left_size < 1 or right_size < 1:
        raise ValueError((left_size, right_size))
    return tuple(
        (left, left_size + right)
        for left in range(left_size)
        for right in range(right_size)
    )


def complete_bipartite_controls(
    graph_rows: Sequence[dict[str, object]], max_n: int = 6
) -> list[dict[str, object]]:
    by_key = {
        (int(row["n"]), int(row["canonical_code"])): row for row in graph_rows
    }
    answer: list[dict[str, object]] = []
    for left_size in range(1, max_n):
        for right_size in range(left_size, max_n - left_size + 1):
            n = left_size + right_size
            if n > max_n or right_size - left_size < 2:
                continue
            edges = complete_bipartite_edges(left_size, right_size)
            code = canonical_code(n, edge_code(n, edges))
            row = by_key[(n, code)]
            formula = quadratic_dimension(n, (left_size, right_size))
            answer.append(
                {
                    "tag": "[COMPUTATION]",
                    "family": f"K_{left_size},{right_size}",
                    "part_sizes": [left_size, right_size],
                    "n": n,
                    "canonical_code": code,
                    "hamiltonian_path_exists": row["hamiltonian_path_exists"],
                    "dimension": row["dimension"],
                    "formula_dimension": formula,
                    "closure_sha256": row["closure_sha256"],
                }
            )
    return answer


def run_all() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    census = run_census()
    matching_grade = run_matching_grade(census)
    graph_rows = census["graphs"]
    if not isinstance(graph_rows, list):
        raise TypeError("census graph rows must be a list")
    branching_rows = [
        row for row in graph_rows if int(row["maximum_degree"]) >= 3
    ]
    nonhamiltonian_rows = [
        row for row in graph_rows if row["hamiltonian_path_exists"] is False
    ]
    certificates = [e6_certificate(row) for row in branching_rows]
    complete_controls = complete_bipartite_controls(graph_rows)
    budget_tick(started, "integrated non-Hamiltonian theorem audit")

    matching_counterexample = matching_grade["counterexamples"][
        "matching_deficiency_and_path_cover"
    ]
    grade_counterexample = matching_grade["counterexamples"][
        "same_grade_set_different_dimension"
    ]
    checks = [
        {
            "name": "all bounded census checks pass",
            "passed": all(bool(check["passed"]) for check in census["checks"]),
            "detail": f"{len(census['checks'])} census checks",
        },
        {
            "name": "all candidate-invariant checks pass",
            "passed": all(
                bool(check["passed"]) for check in matching_grade["checks"]
            ),
            "detail": f"{len(matching_grade['checks'])} invariant checks",
        },
        {
            "name": "every branching spanning basis contains induced E6",
            "passed": bool(certificates)
            and all(
                certificate["tree_generator_gf2_rank"]
                == certificate["expected_basis_rank"]
                and certificate["all_tree_generators_have_Q_one"]
                and certificate["e6_induced"]
                for certificate in certificates
            ),
            "detail": f"{len(certificates)} exact spanning-tree certificates",
        },
        {
            "name": "every non-Hamiltonian graph reaches the complete Q=1 nonradical level",
            "passed": bool(nonhamiltonian_rows)
            and all(
                row["maximum_degree"] >= 3
                and row["closure_equals_quadratic_roots"]
                and row["dimension"] == row["predicted_dimension"]
                for row in nonhamiltonian_rows
            ),
            "detail": f"{len(nonhamiltonian_rows)} non-Hamiltonian representatives",
        },
        {
            "name": "matching deficiency and path cover are not complete invariants",
            "passed": matching_counterexample["left"]["matching_deficiency"]
            == matching_counterexample["right"]["matching_deficiency"]
            and matching_counterexample["left"]["minimum_path_cover"]
            == matching_counterexample["right"]["minimum_path_cover"]
            and matching_counterexample["left"]["dimension"]
            != matching_counterexample["right"]["dimension"],
            "detail": "exact same-invariants/different-dimension witness",
        },
        {
            "name": "occupied Clifford grades are not a complete invariant",
            "passed": grade_counterexample["left"]["profile"]["grades"]
            == grade_counterexample["right"]["profile"]["grades"]
            and grade_counterexample["left"]["dimension"]
            != grade_counterexample["right"]["dimension"],
            "detail": "same n and occupied grades, different exact dimensions",
        },
        {
            "name": "complete-bipartite non-Hamiltonian controls match the theorem",
            "passed": bool(complete_controls)
            and all(
                row["hamiltonian_path_exists"] is False
                and row["dimension"] == row["formula_dimension"]
                for row in complete_controls
            ),
            "detail": f"{len(complete_controls)} K_r,s controls through n=6",
        },
        {
            "name": "resource cap",
            "passed": max_rss_bytes() < RSS_CAP_BYTES,
            "detail": f"peak_rss_bytes={max_rss_bytes()}",
        },
    ]
    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError(checks)

    data = {
        "verdict": {
            "tag": "[THEOREM]+[EXTERNAL]",
            "headline": "Every connected bipartite graph with maximum degree at least three, Hamiltonian or not, generates exactly the nonradical Q_c=1 Pauli labels.",
            "nonhamiltonian_classification": "A connected bipartite graph without a Hamiltonian path necessarily has maximum degree at least three and therefore lies in the quadratic-root branch.",
            "dimension_formula": {
                "n_odd": "2^(2n-2)-1",
                "n_even": "2^(2n-2)-(-1)^r 2^(n-1), where r is either colour-class size",
            },
            "ising_scope": "finite-graph individual-local-term Lie algebra only; no partition-function, critical-point, transfer-spectrum, thermodynamic-limit, or 3D-solvability conclusion",
        },
        "classification": {
            "tag": "[THEOREM]",
            "algebra": "Lie_Q { i X_v, i Z_u Z_v : v in V, uv in E }",
            "packed_label_space": "H={(a,b) in F_2^n x F_2^n: |b| even}",
            "symplectic_form": "B((a,b),(a',b'))=a.b'+b.a'",
            "quadratic_form": "Q_c(a,b)=a.b+|a|+c.b mod 2, with c a bipartition indicator",
            "radical": "span{P_X=(1,...,1;0)}",
            "branches": [
                {
                    "condition": "connected path P_n",
                    "dimension": "n(2n-1)",
                },
                {
                    "condition": "connected even cycle C_n",
                    "dimension": "2n(2n-1)",
                },
                {
                    "condition": "connected bipartite and Delta>=3",
                    "label_set": "{p in H: Q_c(p)=1} minus {P_X when n is odd}",
                    "dimension": "2^(2n-2)-1 for odd n; 2^(2n-2)-(-1)^r 2^(n-1) for even n",
                },
            ],
            "complete_bipartite_corollary": "For K_r,s with |r-s|>=2, the graph is non-Hamiltonian and the branching formula applies for all r,s>=1.",
        },
        "external_orbit_input": {
            "tag": "[EXTERNAL]",
            "source": "Ahmet Seven, Orbits of groups generated by transvections over F_2, arXiv:math/0303098v1 (2003)",
            "url": "https://arxiv.org/abs/math/0303098",
            "result_used": "Theorem 3.2: for a basis equivalent to a tree containing induced E6, the transvection group has exactly the Q=0 and Q=1 orbits outside the radical.",
            "application_certificate": "A branching spanning tree gives a raw-generator basis whose commutation graph is its subdivision and contains induced E6 on one field, three incident bonds, and two neighbouring fields.",
        },
        "census": census,
        "candidate_invariants": matching_grade,
        "spanning_tree_e6_certificates": certificates,
        "complete_bipartite_finite_controls": complete_controls,
        "limits": {
            "tag": "[UNRESOLVED]",
            "finite_census_scope": "connected bipartite graphs through n=6 only",
            "proof_scope": "all finite connected bipartite graphs, conditional only on the cited external transvection-orbit theorem",
            "excluded": "two-generator algebras, thermodynamic limits, transfer spectra, and exact solution of the three-dimensional Ising model",
        },
    }
    return data, checks


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]]
) -> dict[str, object]:
    sources = (
        "experiments/e224_nonhamiltonian_census.py",
        "experiments/e225_matching_grade.py",
        "experiments/e226_nonhamiltonian_theorem.py",
        "tests/test_nonhamiltonian_bipartite.py",
        "proofs/nonhamiltonian_bipartite.md",
    )
    return {
        "meta": {
            "experiment": "e226_nonhamiltonian_theorem",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e226_nonhamiltonian_theorem.py",
            "working_directory": str(ROOT),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": "exact integer graph enumeration, GF(2) symplectic Pauli closure, exact matching/path-cover DP, and exact Jordan--Wigner basis inversion",
            "budgets": {
                "process_time_seconds": CPU_BUDGET_SECONDS,
                "rss_cap_bytes": RSS_CAP_BYTES,
            },
            "source_sha256": {
                source: sha256_file(ROOT / source) for source in sources
            },
        },
        "data": data,
        "checks": checks,
    }


def main() -> int:
    data, checks = run_all()
    artifact = make_artifact(data, checks)
    ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS ({len(checks)} top-level checks, artifact={ARTIFACT}, "
        f"peak_rss={max_rss_bytes()})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
