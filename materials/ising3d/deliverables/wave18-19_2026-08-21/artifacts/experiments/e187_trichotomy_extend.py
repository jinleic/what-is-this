"""Integrate the trichotomy extension, exact counterexample, and corollaries."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from math import comb
from pathlib import Path
from typing import Iterable, Sequence

from e185_order_counterexample import (
    CPU_BUDGET_SECONDS,
    RSS_CAP_BYTES,
    bipartition,
    budget_tick,
    closure_digest,
    max_rss_bytes,
    maximum_degree,
    normalize_edges,
    raw_local_generators,
    run_order_audit,
)
from e186_nonbipartite import chord_grades, run_nonbipartite_audit

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "algebra_growth" / "trichotomy_extend.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_edges(rows: int, cols: int) -> tuple[tuple[int, int], ...]:
    if rows < 1 or cols < 1:
        raise ValueError("grid dimensions must be positive")
    edges: list[tuple[int, int]] = []
    for row in range(rows):
        for col in range(cols):
            vertex = row * cols + col
            if col + 1 < cols:
                edges.append((vertex, vertex + 1))
            if row + 1 < rows:
                edges.append((vertex, vertex + cols))
    return normalize_edges(rows * cols, edges)


def grid_snake_path(rows: int, cols: int) -> tuple[int, ...]:
    path: list[int] = []
    for row in range(rows):
        columns: Iterable[int] = range(cols) if row % 2 == 0 else range(cols - 1, -1, -1)
        path.extend(row * cols + col for col in columns)
    return tuple(path)


def box_vertex(a: int, b: int, x: int, y: int, z: int) -> int:
    return (z * b + y) * a + x


def box_edges(a: int, b: int, c: int) -> tuple[tuple[int, int], ...]:
    if min(a, b, c) < 1:
        raise ValueError("box dimensions must be positive")
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


def box_snake_path(a: int, b: int, c: int) -> tuple[int, ...]:
    layer: list[int] = []
    for y in range(b):
        xs: Iterable[int] = range(a) if y % 2 == 0 else range(a - 1, -1, -1)
        layer.extend(box_vertex(a, b, x, y, 0) for x in xs)
    path: list[int] = []
    layer_size = a * b
    for z in range(c):
        offsets = layer if z % 2 == 0 else list(reversed(layer))
        path.extend(vertex + z * layer_size for vertex in offsets)
    return tuple(path)


def validate_path(
    n: int, edges: Sequence[tuple[int, int]], path: Sequence[int]
) -> bool:
    if len(path) != n or set(path) != set(range(n)):
        return False
    edge_set = set(normalize_edges(n, edges))
    return all(
        ((left, right) if left < right else (right, left)) in edge_set
        for left, right in zip(path, path[1:])
    )


def branch_for_sides(sides: Sequence[int]) -> str:
    nontrivial = sorted(side for side in sides if side > 1)
    if len(nontrivial) <= 1:
        return "path"
    if nontrivial == [2, 2]:
        return "cycle"
    return "branching"


def branching_dimension(n: int) -> int:
    if n % 2:
        return (1 << (2 * n - 2)) - 1
    return (1 << (2 * n - 2)) - ((-1) ** (n // 2)) * (1 << (n - 1))


def branch_dimension(branch: str, n: int) -> int:
    if branch == "path":
        return n * (2 * n - 1)
    if branch == "cycle":
        return 2 * n * (2 * n - 1)
    if branch == "branching":
        return branching_dimension(n)
    raise ValueError(branch)


def direct_branching_binomial_dimension(n: int) -> int:
    N = 2 * n
    return sum(comb(N, grade) for grade in range(2, N + 1, 4) if grade != N)


def concrete_grid_row(
    *,
    label: str,
    rows: int,
    cols: int,
    previous_status_class: str,
    previous_status: str,
    previous_evidence: Sequence[str],
) -> dict[str, object]:
    n = rows * cols
    edges = grid_edges(rows, cols)
    path = grid_snake_path(rows, cols)
    if not validate_path(n, edges, path):
        raise AssertionError(f"invalid grid Hamiltonian path for {rows}x{cols}")
    connected, is_bipartite, _ = bipartition(n, edges)
    branch = branch_for_sides((rows, cols))
    grades = chord_grades(n, edges, path)
    dimension = branch_dimension(branch, n)
    if branch == "branching":
        if not grades or not all(grade % 4 == 2 for grade in grades):
            raise AssertionError(f"bad branching chord grades for {rows}x{cols}: {grades}")
        if not any(6 <= grade <= 2 * n - 4 for grade in grades):
            raise AssertionError(f"no corrected grade seed for {rows}x{cols}: {grades}")
        if dimension != direct_branching_binomial_dimension(n):
            raise AssertionError(f"binomial formula mismatch for n={n}")
    generators = raw_local_generators(n, edges)
    return {
        "tag": "[THEOREM]",
        "label": label,
        "geometry": "open_rectangular_grid",
        "shape": [rows, cols],
        "n": n,
        "N_majoranas": 2 * n,
        "edge_count": len(edges),
        "raw_generator_count": len(generators),
        "raw_generator_sha256": closure_digest(generators),
        "hamiltonian_path": list(path),
        "hamiltonian_path_valid": True,
        "bipartite": is_bipartite,
        "connected": connected,
        "maximum_degree": maximum_degree(n, edges),
        "branch": branch,
        "chord_grades": list(grades),
        "grade_space": "all noncentral Clifford grades k=2 mod 4",
        "dimension": dimension,
        "dimension_by_binomial_sum": direct_branching_binomial_dimension(n)
        if branch == "branching"
        else dimension,
        "previous_status_class": previous_status_class,
        "previous_status": previous_status,
        "previous_evidence": list(previous_evidence),
        "two_generator_scope_separate": True,
    }


def finite_box_audit(max_side: int = 4, max_volume: int = 24) -> dict[str, int]:
    cases = 0
    path_cases = 0
    cycle_cases = 0
    branching_cases = 0
    for a in range(1, max_side + 1):
        for b in range(1, max_side + 1):
            for c in range(1, max_side + 1):
                n = a * b * c
                if n > max_volume:
                    continue
                edges = box_edges(a, b, c)
                path = box_snake_path(a, b, c)
                if not validate_path(n, edges, path):
                    raise AssertionError(f"box path failed at {(a, b, c)}")
                connected, is_bipartite, _ = bipartition(n, edges)
                if not connected or not is_bipartite:
                    raise AssertionError(f"box graph failed at {(a, b, c)}")
                branch = branch_for_sides((a, b, c))
                if branch == "path":
                    path_cases += 1
                elif branch == "cycle":
                    cycle_cases += 1
                else:
                    branching_cases += 1
                    grades = chord_grades(n, edges, path)
                    if not grades or not all(grade % 4 == 2 for grade in grades):
                        raise AssertionError(f"box grade failure at {(a, b, c)}")
                    if branch_dimension(branch, n) != direct_branching_binomial_dimension(n):
                        raise AssertionError(f"box dimension failure at {(a, b, c)}")
                cases += 1
    return {
        "cases": cases,
        "path_cases": path_cases,
        "cycle_cases": cycle_cases,
        "branching_cases": branching_cases,
    }


def repository_audit() -> list[dict[str, object]]:
    return [
        {
            "tag": "[EXTERNAL]",
            "object": "open 2x6 local-term algebra",
            "prior_status": "previously known exactly",
            "evidence": [
                "proofs/alll_saturation.md",
                "notes/hypotheses.csv:H431-H433",
            ],
            "detail": "L=6 was fully enumerated at dimension 4,192,256.",
        },
        {
            "tag": "[EXTERNAL]",
            "object": "open 3x4 local-term algebra",
            "prior_status": "previously known exactly",
            "evidence": [
                "proofs/clifford_grade_classification.md",
                "notes/hypotheses.csv:H486",
            ],
            "detail": "The wave-18 grade theorem explicitly recorded dimension 4,192,256 without materializing the closure.",
        },
        {
            "tag": "[EXTERNAL]",
            "object": "open 3x5 local-term algebra",
            "prior_status": "previously implicit, not graph-specifically tabled",
            "evidence": [
                "proofs/clifford_grade_classification.md",
                "results/algebra_growth/clifford_grade.json",
            ],
            "detail": "The wave-18 theorem and its n=15 binomial check imply the value; no prior 3x5 local-term row was found.",
        },
        {
            "tag": "[EXTERNAL]",
            "object": "open 4x4 local-term algebra",
            "prior_status": "previously implicit; the separately open object was the two-generator algebra",
            "evidence": [
                "proofs/clifford_grade_classification.md",
                "results/algebra_growth/clifford_grade.json",
                "proofs/char0_4x4.md",
            ],
            "detail": "The n=16 binomial value was already checked.  The 1,794 lower bound in char0_4x4 is for Lie<A=sum X,B=sum ZZ>, not for the individual local terms.",
        },
        {
            "tag": "[EXTERNAL]",
            "object": "simple-cubic a x b x c finite-box local-term algebra",
            "prior_status": "not previously stated as a box corollary",
            "evidence": ["proofs/clifford_grade_classification.md"],
            "detail": "It is an immediate finite-graph consequence once an explicit box Hamiltonian path is supplied.",
        },
    ]


def run_corollary_audit() -> dict[str, object]:
    started = time.process_time()
    rows = [
        concrete_grid_row(
            label="grid_2x6",
            rows=2,
            cols=6,
            previous_status_class="previously_known",
            previous_status="Previously known exactly as the L=6 local-term algebra: 4,192,256.",
            previous_evidence=(
                "proofs/alll_saturation.md",
                "notes/hypotheses.csv:H431-H433",
            ),
        ),
        concrete_grid_row(
            label="grid_3x4",
            rows=3,
            cols=4,
            previous_status_class="previously_known",
            previous_status="Previously known exactly from the wave-18 grade classification: 4,192,256.",
            previous_evidence=(
                "proofs/clifford_grade_classification.md",
                "notes/hypotheses.csv:H486",
            ),
        ),
        concrete_grid_row(
            label="grid_3x5",
            rows=3,
            cols=5,
            previous_status_class="previously_implicit_not_tabled",
            previous_status="The wave-18 theorem and stored n=15 binomial check already implied the value, but no graph-specific 3x5 local-term row was found.",
            previous_evidence=(
                "proofs/clifford_grade_classification.md",
                "results/algebra_growth/clifford_grade.json",
            ),
        ),
        concrete_grid_row(
            label="grid_4x4",
            rows=4,
            cols=4,
            previous_status_class="previously_implicit_two_generator_open",
            previous_status="The wave-18 theorem and stored n=16 check implied this local-term value.  The repository's 1,794 bound/open closure is for the different two-generator algebra Lie<A,B>.",
            previous_evidence=(
                "proofs/clifford_grade_classification.md",
                "results/algebra_growth/clifford_grade.json",
                "proofs/char0_4x4.md",
            ),
        ),
    ]
    expected = {
        "grid_2x6": 4_192_256,
        "grid_3x4": 4_192_256,
        "grid_3x5": 268_435_455,
        "grid_4x4": 1_073_709_056,
    }
    if {row["label"]: row["dimension"] for row in rows} != expected:
        raise AssertionError("concrete corollary values do not match")
    box_audit = finite_box_audit()
    budget_tick(started, "corollaries")
    return {
        "tag": "[THEOREM]",
        "concrete_rows": rows,
        "rectangular_grid_family": {
            "tag": "[THEOREM]",
            "scope": "open r x s grids with positive integer sides",
            "n": "r*s",
            "hamiltonian_path": "row snake",
            "branch_cases": [
                {
                    "condition": "min(r,s)=1",
                    "branch": "path",
                    "dimension": "n(2n-1)",
                },
                {
                    "condition": "r=s=2",
                    "branch": "cycle",
                    "dimension": "2n(2n-1)=56",
                },
                {
                    "condition": "r,s>=2 and (r,s)!=(2,2)",
                    "branch": "branching",
                    "dimension": "2^(2n-2)-(-1)^(n/2)2^(n-1) for even n; 2^(2n-2)-1 for odd n",
                },
            ],
            "previous_status": "The family statement was implicit in the wave-18 Hamiltonian-path theorem; this front records the corollary and audits named rows.",
        },
        "simple_cubic_box_family": {
            "tag": "[THEOREM]",
            "scope": "open a x b x c boxes with positive integer sides",
            "n": "a*b*c",
            "hamiltonian_path": "alternate a two-dimensional layer snake with its reversal between successive layers",
            "branch_cases": [
                {
                    "condition": "at most one side exceeds 1",
                    "branch": "path",
                    "dimension": "n(2n-1)",
                },
                {
                    "condition": "nontrivial side multiset is {2,2}",
                    "branch": "cycle",
                    "dimension": "2n(2n-1)=56",
                },
                {
                    "condition": "all other boxes; in particular every genuine box a,b,c>=2",
                    "branch": "branching",
                    "dimension": "2^(2n-2)-(-1)^(n/2)2^(n-1) for even n; 2^(2n-2)-1 for odd n",
                    "grade_space": "all noncentral Clifford grades k=2 mod 4",
                },
            ],
            "finite_path_audit": box_audit,
            "previous_status": "Not previously stated as a simple-cubic finite-box local-term corollary; it was implicit in the wave-18 theorem.",
            "scope_exclusions": [
                "thermodynamic limit",
                "transfer spectrum",
                "two-generator algebra",
                "solvability of the three-dimensional Ising model",
            ],
        },
        "repository_audit": repository_audit(),
        "process_time_seconds": time.process_time() - started,
        "peak_rss_bytes": max_rss_bytes(),
    }


def run_all() -> tuple[dict[str, object], list[dict[str, object]]]:
    started = time.process_time()
    order = run_order_audit()
    nonbipartite = run_nonbipartite_audit()
    corollaries = run_corollary_audit()
    budget_tick(started, "integrated audit")
    checks = [
        {
            "name": "all arbitrary-order checks",
            "passed": all(bool(row["passed"]) for row in order["checks"]),
            "detail": f"{len(order['checks'])} checks",
        },
        {
            "name": "all non-bipartite checks",
            "passed": all(bool(row["passed"]) for row in nonbipartite["checks"]),
            "detail": f"{len(nonbipartite['checks'])} checks",
        },
        {
            "name": "Hamiltonian-path hypothesis is necessary",
            "passed": order["counterexamples"][0]["dimension"] == 72
            and order["counterexamples"][0][
                "hamiltonian_branching_formula_if_misapplied"
            ]
            == 56,
            "detail": "K1,3: exact local-term dimension 72, not 56",
        },
        {
            "name": "general connected non-bipartite branch",
            "passed": {
                row["dimension"] for row in nonbipartite["graph_rows"]
            }
            >= {90, 126, 510, 2046},
            "detail": "odd-cycle exception and full-even branching controls",
        },
        {
            "name": "all concrete corollaries supported",
            "passed": all(
                row["hamiltonian_path_valid"]
                and row["bipartite"]
                and row["dimension"] == row["dimension_by_binomial_sum"]
                and bool(row["previous_status"])
                for row in corollaries["concrete_rows"]
            ),
            "detail": f"{len(corollaries['concrete_rows'])} rows",
        },
        {
            "name": "new graph-specific dimensions",
            "passed": {
                row["label"]: row["dimension"]
                for row in corollaries["concrete_rows"]
            }
            == {
                "grid_2x6": 4_192_256,
                "grid_3x4": 4_192_256,
                "grid_3x5": 268_435_455,
                "grid_4x4": 1_073_709_056,
            },
            "detail": "2x6, 3x4, 3x5, 4x4",
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
            "tag": "[THEOREM]",
            "hamiltonian_path_removal": "REFUTED for connected bipartite graphs by K1,3 (72 != 56)",
            "arbitrary_order_result": "Every edge has grade 2 mod 4 exactly when position parity is a proper two-colouring; for a connected bipartite graph such orders exist exactly when the colour-class sizes differ by at most one.",
            "nonbipartite_result": "For every connected simple non-bipartite graph, an odd cycle has grades {2,2n-2} and dimension 2n(2n-1); every other graph (equivalently Delta>=3) is the full noncentral even Pauli/Clifford algebra of dimension 2^(2n-1)-2.",
            "three_dimensional_scope": "finite-box local-term algebra only; not a thermodynamic, spectral, two-generator, or solvability theorem",
        },
        "arbitrary_order": order,
        "nonbipartite": nonbipartite,
        "corollaries": corollaries,
    }
    return data, checks


def make_artifact(
    data: dict[str, object], checks: list[dict[str, object]]
) -> dict[str, object]:
    sources = (
        "experiments/e185_order_counterexample.py",
        "experiments/e186_nonbipartite.py",
        "experiments/e187_trichotomy_extend.py",
    )
    return {
        "meta": {
            "experiment": "e187_trichotomy_extend",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "command": ".venv/bin/python experiments/e187_trichotomy_extend.py",
            "working_directory": str(ROOT),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "arithmetic": "exact integer binomial sums, GF(2) symplectic Pauli closure, finite graph/order enumeration, and exact Clifford-support overlap calculus",
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
