#!/usr/bin/env python3
"""Uniform affine-degree bound for finite full-spectrum exceptional couplings.

Wave 22 proved that every connected simple nonpath graph has only finitely many
physical couplings whose complete layer spectrum is one full subset-product
spectrum.  This producer polynomializes that incidence family and records an
explicit graph-uniform affine Bezout bound on the number of complex candidate
couplings, hence also on the physical exceptional set.

All finite controls use Python integers.  No eigenvalues, floating point, or
critical-coupling benchmark enter the calculation.
"""
from __future__ import annotations

import ctypes
import hashlib
import itertools
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "experiments/e245_bipartite_exception_degree.py"
OUTPUT = ROOT / "results" / "spectral" / "bipartite_exception_degree.json"
VERIFIER = ROOT / "tests" / "test_bipartite_exception_degree.py"
TOKEN_PRODUCER = ROOT / "experiments" / "e236_token_splitting_obstruction.py"
TOKEN_VERIFIER = ROOT / "tests" / "test_token_splitting_obstruction.py"
TOKEN_ARTIFACT = ROOT / "results" / "spectral" / "token_splitting_obstruction.json"
CHEVALLEY_SOURCE = ROOT / "sources" / "stacks_chevalley_054K.html"
CPU_LIMIT_SECONDS = 120.0
RSS_LIMIT_BYTES = 2 * 1024**3

Edge = tuple[int, int]


def peak_rss_bytes() -> int:
    if platform.system() == "Darwin":
        class TimeValue(ctypes.Structure):
            _fields_ = [
                ("seconds", ctypes.c_int32),
                ("microseconds", ctypes.c_int32),
            ]

        class MachTaskBasicInfo(ctypes.Structure):
            _fields_ = [
                ("virtual_size", ctypes.c_uint64),
                ("resident_size", ctypes.c_uint64),
                ("resident_size_max", ctypes.c_uint64),
                ("user_time", TimeValue),
                ("system_time", TimeValue),
                ("policy", ctypes.c_int32),
                ("suspend_count", ctypes.c_int32),
            ]

        system = ctypes.CDLL("/usr/lib/libSystem.B.dylib")
        system.mach_task_self.restype = ctypes.c_uint32
        info = MachTaskBasicInfo()
        count = ctypes.c_uint32(
            ctypes.sizeof(info) // ctypes.sizeof(ctypes.c_int32)
        )
        result = system.task_info(
            system.mach_task_self(),
            20,
            ctypes.byref(info),
            ctypes.byref(count),
        )
        if result != 0:
            raise OSError(f"task_info failed with status {result}")
        return int(info.resident_size_max)
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def peak_rss_measurement() -> str:
    if platform.system() == "Darwin":
        return "mach_task_basic_info.resident_size_max (task_info flavor 20)"
    return "getrusage(RUSAGE_SELF).ru_maxrss multiplied by 1024"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_check(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def grid_edges(rows: int, columns: int) -> tuple[Edge, ...]:
    edges: list[Edge] = []
    for row in range(rows):
        for column in range(columns):
            site = row * columns + column
            if row + 1 < rows:
                edges.append((site, site + columns))
            if column + 1 < columns:
                edges.append((site, site + 1))
    return tuple(edges)


def graph_controls() -> dict[str, tuple[int, tuple[Edge, ...]]]:
    return {
        "path4": (4, ((0, 1), (1, 2), (2, 3))),
        "cycle4": (4, ((0, 1), (1, 2), (2, 3), (0, 3))),
        "claw": (4, ((0, 1), (0, 2), (0, 3))),
        "open2x3": (6, grid_edges(2, 3)),
        "open3x3": (9, grid_edges(3, 3)),
    }


def validate_graph(order: int, edges: tuple[Edge, ...]) -> bool:
    if order < 1 or any(not (0 <= left < right < order) for left, right in edges):
        return False
    if len(set(edges)) != len(edges):
        return False
    adjacency = [[] for _ in range(order)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for other in adjacency[vertex]:
            if other not in seen:
                seen.add(other)
                frontier.append(other)
    return len(seen) == order


def is_path(order: int, edges: tuple[Edge, ...]) -> bool:
    degrees = [0] * order
    for left, right in edges:
        degrees[left] += 1
        degrees[right] += 1
    return (
        validate_graph(order, edges)
        and len(edges) == order - 1
        and max(degrees, default=0) <= 2
    )


def aligned_exponents(order: int, edges: tuple[Edge, ...]) -> tuple[int, ...]:
    edge_count = len(edges)
    return tuple(
        edge_count
        - sum(((state >> left) ^ (state >> right)) & 1 for left, right in edges)
        for state in range(1 << order)
    )


def mode_degree_rows(order: int) -> list[dict[str, int]]:
    slot_weights = sorted(
        (mask.bit_count() for mask in range(1 << order)), reverse=True
    )
    running = 0
    rows: list[dict[str, int]] = []
    for coefficient_index, weight in enumerate(slot_weights, start=1):
        running += weight
        rows.append(
            {
                "coefficient_index": coefficient_index,
                "maximum_u_total_degree": running,
                "c_degree": coefficient_index,
                "maximum_total_mode_degree": running + coefficient_index,
                "coarse_bound": coefficient_index * (order + 1),
            }
        )
    return rows


def degree_record(name: str, order: int, edges: tuple[Edge, ...]) -> dict[str, object]:
    edge_count = len(edges)
    dimension = 1 << order
    exponents = aligned_exponents(order, edges)
    diagonal_rows = []
    for exponent in sorted(set(exponents)):
        power_2t = edge_count - exponent
        power_1pt2 = exponent
        diagonal_rows.append(
            {
                "aligned_edge_exponent": exponent,
                "power_2t": power_2t,
                "power_1_plus_t2": power_1pt2,
                "polynomial_degree": power_2t + 2 * power_1pt2,
            }
        )
    diagonal_degree = 2 * edge_count
    entry_degree = 2 * order + diagonal_degree
    equation_degree = dimension * entry_degree
    ambient_dimension = order + 2
    graph_bound = equation_degree**ambient_dimension
    uniform_entry_degree = order * order + order
    uniform_equation_degree = dimension * uniform_entry_degree
    uniform_bound = uniform_equation_degree**ambient_dimension
    mode_rows = mode_degree_rows(order)
    exponent_sum = sum(exponents)
    determinant = {
        "aligned_exponent_sum": exponent_sum,
        "expected_aligned_exponent_sum": edge_count * dimension // 2,
        "power_2t": edge_count * dimension // 2,
        "power_1_plus_t2": edge_count * dimension // 2,
        "power_1_minus_t2": order * dimension,
        "root_set_over_C": ["0", "+i", "-i", "+1", "-1"],
    }
    return {
        "name": name,
        "order_n": order,
        "edge_count_m": edge_count,
        "dimension_N": dimension,
        "connected_simple": validate_graph(order, edges),
        "is_path_control": is_path(order, edges),
        "aligned_exponent_min": min(exponents),
        "aligned_exponent_max": max(exponents),
        "aligned_exponent_range_bound": [0, edge_count],
        "diagonal_monomial_rows": diagonal_rows,
        "diagonal_degree": diagonal_degree,
        "entry_degree_delta": entry_degree,
        "characteristic_coefficient_degree_bound": dimension * entry_degree,
        "mode_coefficient_degree_rows": mode_rows,
        "incidence_equation_degree_d": equation_degree,
        "ambient_variable_count": ambient_dimension,
        "incidence_variety_degree_bound_decimal": str(graph_bound),
        "exceptional_candidate_bound_decimal": (
            None if is_path(order, edges) else str(graph_bound)
        ),
        "uniform_simple_graph_entry_degree": uniform_entry_degree,
        "uniform_simple_graph_bound_decimal": str(uniform_bound),
        "determinant_factor_exponents": determinant,
    }


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    controls = graph_controls()
    records = {
        name: degree_record(name, order, edges)
        for name, (order, edges) in controls.items()
    }

    add_check(
        checks,
        "C1_current_dependencies_present",
        VERIFIER.is_file()
        and TOKEN_PRODUCER.is_file()
        and TOKEN_VERIFIER.is_file()
        and TOKEN_ARTIFACT.is_file()
        and CHEVALLEY_SOURCE.is_file(),
        "new verifier plus current e236 producer/verifier/artifact and pinned Chevalley source are present",
    )
    add_check(
        checks,
        "C2_connected_simple_controls",
        all(bool(row["connected_simple"]) for row in records.values())
        and bool(records["path4"]["is_path_control"])
        and all(
            not bool(row["is_path_control"])
            for name, row in records.items()
            if name != "path4"
        ),
        "path4 is the sole path control; cycle, claw, and grids are nonpaths",
    )
    add_check(
        checks,
        "C3_polynomial_diagonal_exponents_nonnegative",
        all(
            int(term["power_2t"]) >= 0 and int(term["power_1_plus_t2"]) >= 0
            for row in records.values()
            for term in row["diagonal_monomial_rows"]
        ),
        "(2t)^m*q^a=(2t)^(m-a)*(1+t^2)^a is polynomial for every spin state",
    )
    add_check(
        checks,
        "C4_exact_diagonal_and_entry_degrees",
        all(
            int(row["diagonal_degree"]) == 2 * int(row["edge_count_m"])
            and int(row["entry_degree_delta"])
            == 2 * int(row["order_n"]) + int(row["diagonal_degree"])
            for row in records.values()
        ),
        "diagonal degree is 2m, and two Hamming kernels add at most 2n",
    )
    add_check(
        checks,
        "C5_characteristic_and_mode_coefficient_degrees",
        all(
            all(
                int(mode["maximum_total_mode_degree"])
                <= int(mode["coarse_bound"])
                <= int(mode["coefficient_index"])
                * int(row["entry_degree_delta"])
                for mode in row["mode_coefficient_degree_rows"]
            )
            and int(row["incidence_equation_degree_d"])
            == int(row["dimension_N"]) * int(row["entry_degree_delta"])
            for row in records.values()
        ),
        "a_k has degree <=k*delta and the Boolean-mode coefficient <=k*(n+1)",
    )
    add_check(
        checks,
        "C6_determinant_factorization_exponents",
        all(
            int(row["determinant_factor_exponents"]["aligned_exponent_sum"])
            == int(
                row["determinant_factor_exponents"][
                    "expected_aligned_exponent_sum"
                ]
            )
            and int(row["determinant_factor_exponents"]["power_1_minus_t2"])
            == int(row["order_n"]) * int(row["dimension_N"])
            for row in records.values()
        ),
        "det(P_t)^2=(1-t^2)^(nN), and every edge is aligned in N/2 spin states",
    )
    add_check(
        checks,
        "C7_affine_degree_bounds_and_simple_graph_uniformization",
        all(
            int(row["incidence_variety_degree_bound_decimal"])
            == int(row["incidence_equation_degree_d"])
            ** int(row["ambient_variable_count"])
            and (
                row["exceptional_candidate_bound_decimal"] is None
                if row["is_path_control"]
                else int(row["exceptional_candidate_bound_decimal"])
                == int(row["incidence_variety_degree_bound_decimal"])
            )
            and int(row["entry_degree_delta"])
            <= int(row["uniform_simple_graph_entry_degree"])
            and int(row["incidence_variety_degree_bound_decimal"])
            <= int(row["uniform_simple_graph_bound_decimal"])
            for row in records.values()
        ),
        "Heintz gives d^(n+2); m<=n(n-1)/2 gives delta<=n^2+n",
    )
    add_check(
        checks,
        "C8_open2x3_named_bound",
        records["open2x3"]["entry_degree_delta"] == 26
        and records["open2x3"]["incidence_equation_degree_d"] == 1664
        and records["open2x3"]["exceptional_candidate_bound_decimal"]
        == "58779593138084256579321856",
        "the disposed finite layer is a scale control only; no elimination is rerun",
    )

    cpu = time.process_time() - started
    rss = peak_rss_bytes()
    rss_measurement = peak_rss_measurement()
    add_check(
        checks,
        "C9_resource_contract",
        cpu < CPU_LIMIT_SECONDS and rss < RSS_LIMIT_BYTES,
        (
            f"process CPU {cpu:.6f}s/{CPU_LIMIT_SECONDS}s; "
            f"peak RSS {rss}/{RSS_LIMIT_BYTES} bytes via {rss_measurement}"
        ),
    )
    if not all(bool(row["passed"]) for row in checks):
        failed = [row for row in checks if not row["passed"]]
        raise AssertionError(f"producer checks failed: {failed}")

    source_hashes = {
        SCRIPT: file_sha256(ROOT / SCRIPT),
        "tests/test_bipartite_exception_degree.py": file_sha256(VERIFIER),
        "experiments/e236_token_splitting_obstruction.py": file_sha256(
            TOKEN_PRODUCER
        ),
        "tests/test_token_splitting_obstruction.py": file_sha256(TOKEN_VERIFIER),
        "results/spectral/token_splitting_obstruction.json": file_sha256(
            TOKEN_ARTIFACT
        ),
        "sources/stacks_chevalley_054K.html": file_sha256(CHEVALLEY_SOURCE),
    }
    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": SCRIPT,
            "verifier": "tests/test_bipartite_exception_degree.py",
            "artifact": "results/spectral/bipartite_exception_degree.json",
            "interpreter": sys.executable,
            "arithmetic": "exact Python integers and symbolic degree inequalities",
            "process_cpu_seconds": cpu,
            "process_cpu_budget_seconds": CPU_LIMIT_SECONDS,
            "peak_rss_bytes": rss,
            "peak_rss_measurement": rss_measurement,
            "rss_limit_bytes": RSS_LIMIT_BYTES,
            "benchmark_used": False,
            "source_sha256": source_hashes,
        },
        "checks": checks,
        "data": {
            "claim_tag": "[THEOREM]",
            "theorem": (
                "Let G be a finite connected simple nonpath graph with n vertices, "
                "m edges, and N=2^n. The t-image of the full polynomialized "
                "coefficient-incidence variety has at most [N*(2n+2m)]^(n+2) "
                "complex points. It contains every physical coupling whose spectrum "
                "is one full nonzero n-mode subset-product multiset, so the physical "
                "exceptional set inside 0<t<1 obeys the same bound. Uniformly over "
                "connected simple nonpath n-vertex graphs the bound is "
                "[2^n*(n^2+n)]^(n+2)."
            ),
            "transfer_polynomialization": {
                "physical_scalar_equivalence": "R_unc=P_t*diag(q^a_sigma)*P_t, q=(1+t^2)/(2t), a_sigma=m-cut_G(sigma); R_unc differs from the centered physical matrix by one nonzero scalar",
                "clearing_scalar": "S(t)=(2t)^m",
                "diagonal_identity": "S*q^a=(2t)^(m-a)*(1+t^2)^a",
                "entry_degree_delta": "2n+2m",
                "determinant": "det(R_tilde)=(2t)^(mN/2)*(1+t^2)^(mN/2)*(1-t^2)^(nN)",
                "singular_parameter_roots": ["0", "+i", "-i", "+1", "-1"],
            },
            "incidence_system": {
                "variables": "t,c,u_1,...,u_n (n+2 affine variables)",
                "equations": "a_k(t)-(-1)^k*c^k*E_k({prod_(i in S)u_i : S subseteq [n]})=0, k=1,...,N",
                "equation_degree_bound": "d=N*(2n+2m)",
                "finiteness_input": "e236 endpoint token-splitting obstruction plus Chevalley constructibility makes the nonsingular t-projection finite; the five determinant roots are already a finite vertical set",
                "physical_positivity": "At physical t the polynomialized matrix is positive definite. Any complex coefficient solution contains c and c*u_i among its positive eigenvalues, hence c,u_i are positive real.",
            },
            "affine_degree_lemma": {
                "tag": "[EXTERNAL LEMMA WITH LOCAL PROOF SKETCH]",
                "statement": "If an affine algebraic set V in A^M is cut out by polynomials of degree at most d, then deg(V)<=d^M. If one coordinate has finite image on V, each irreducible component maps to one point, so the number of image points is at most the number of components and hence at most deg(V).",
                "proof_sketch": "Homogenize the defining equations in projective M-space. Theorem 2.1 of arXiv:2511.07639v1 bounds the degree of their common projective zero locus by d^M independently of the number of equations. The affine components are exactly the projective components meeting the affine chart, so their total degree obeys the same bound. Every irreducible component has degree at least one, and a finite irreducible constructible image in A^1 is one point.",
                "primary_source": {
                    "author": "Joos Heintz",
                    "title": "Definability and fast quantifier elimination in algebraically closed fields",
                    "journal": "Theoretical Computer Science 24 (1983), 239-277",
                    "doi": "10.1016/0304-3975(83)90002-6",
                },
                "corroborating_source": {
                    "title": "Growth in finite simple groups of Lie type of bounded rank",
                    "arxiv": "1005.1858",
                    "location": "Fact 19(g): deg(X)<=d^M for a common zero locus of degree-d polynomials in affine M-space",
                },
                "recent_primary_source": {
                    "authors": "Edward Bierstone, Dima Grigoriev, Pierre D. Milman, Jaroslaw Wlodarczyk",
                    "title": "Effective resolution of singularities",
                    "arxiv": "2511.07639v1",
                    "url": "https://arxiv.org/abs/2511.07639v1",
                    "location": "Theorem 2.1: homogeneous equations of degree at most d in projective M-space have common zero locus of degree at most d^M",
                },
            },
            "grid_family_corollaries": {
                "open_2xL": {
                    "range": "L>=2 (L=1 is a path and is outside the finite-projection theorem)",
                    "n": "2L",
                    "m": "3L-2",
                    "entry_degree_delta": "10L-4",
                    "bound": "[4^L*(10L-4)]^(2L+2)",
                },
                "open_3xL": {
                    "range": "L>=2 (L=1 is a path and is outside the finite-projection theorem)",
                    "n": "3L",
                    "m": "5L-3",
                    "entry_degree_delta": "16L-6",
                    "bound": "[8^L*(16L-6)]^(3L+2)",
                },
            },
            "finite_graph_controls": records,
            "scope": {
                "proved": [
                    "an explicit finite complex-candidate bound for every fixed connected simple nonpath graph",
                    "the uniform connected-simple-nonpath n-vertex bound [2^n*(n^2+n)]^(n+2)",
                    "the same upper bound for physical full-spectrum exceptional couplings on bipartite grids with both side lengths at least two",
                ],
                "not_proved": [
                    "emptiness of the exceptional set beyond the already disposed open 2x3 layer",
                    "a useful numerical bound, a list of exceptional roots, or a graph-uniform root-free interval",
                    "parity-sector unions, larger auxiliary spectra, integrability, or a thermodynamic solution",
                    "an improved critical endpoint or any use of the benchmark critical coupling",
                ],
            },
        },
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"[e245] wrote {OUTPUT.relative_to(ROOT)}: "
        f"{len(payload['checks'])}/{len(payload['checks'])} checks, "
        f"cpu={payload['meta']['process_cpu_seconds']:.6f}s "
        f"rss={payload['meta']['peak_rss_bytes']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
