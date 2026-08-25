#!/usr/bin/env python3
"""Clean-room verifier for the graph-uniform exceptional-coupling bound."""
from __future__ import annotations

import ctypes
import hashlib
import json
import math
import platform
import resource
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "bipartite_exception_degree.json"
PRODUCER = ROOT / "experiments" / "e245_bipartite_exception_degree.py"
TOKEN_PRODUCER = ROOT / "experiments" / "e236_token_splitting_obstruction.py"
TOKEN_VERIFIER = ROOT / "tests" / "test_token_splitting_obstruction.py"
TOKEN_ARTIFACT = ROOT / "results" / "spectral" / "token_splitting_obstruction.json"
CHEVALLEY_SOURCE = ROOT / "sources" / "stacks_chevalley_054K.html"
RSS_LIMIT_BYTES = 2 * 1024**3
FAILURES: list[str] = []

Edge = tuple[int, int]


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def rectangular_edges(height: int, width: int) -> tuple[Edge, ...]:
    horizontal = tuple(
        (row * width + column, row * width + column + 1)
        for row in range(height)
        for column in range(width - 1)
    )
    vertical = tuple(
        (row * width + column, (row + 1) * width + column)
        for row in range(height - 1)
        for column in range(width)
    )
    return tuple(sorted(horizontal + vertical))


def controls() -> dict[str, tuple[int, tuple[Edge, ...]]]:
    return {
        "path4": (4, ((0, 1), (1, 2), (2, 3))),
        "cycle4": (4, ((0, 1), (0, 3), (1, 2), (2, 3))),
        "claw": (4, ((0, 1), (0, 2), (0, 3))),
        "open2x3": (6, rectangular_edges(2, 3)),
        "open3x3": (9, rectangular_edges(3, 3)),
    }


def aligned_exponents(order: int, edges: tuple[Edge, ...]) -> list[int]:
    values: list[int] = []
    for bits in range(1 << order):
        cut = 0
        for left, right in edges:
            left_spin = -1 if bits & (1 << left) else 1
            right_spin = -1 if bits & (1 << right) else 1
            cut += int(left_spin != right_spin)
        values.append(len(edges) - cut)
    return values


def mode_degree_maxima(order: int) -> list[int]:
    weights = sorted(
        [sum((slot >> bit) & 1 for bit in range(order)) for slot in range(1 << order)],
        reverse=True,
    )
    output: list[int] = []
    running = 0
    for index, weight in enumerate(weights, start=1):
        running += weight
        output.append(index + running)
    return output


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    token_artifact = json.loads(TOKEN_ARTIFACT.read_text())
    meta = artifact["meta"]
    data = artifact["data"]
    hashes = meta["source_sha256"]
    check(
        "current producer hash",
        hashes["experiments/e245_bipartite_exception_degree.py"]
        == file_sha256(PRODUCER),
    )
    check(
        "current verifier hash",
        hashes["tests/test_bipartite_exception_degree.py"]
        == file_sha256(Path(__file__).resolve()),
    )
    check(
        "current e236 and Chevalley dependency hashes",
        hashes["experiments/e236_token_splitting_obstruction.py"]
        == file_sha256(TOKEN_PRODUCER)
        and hashes["tests/test_token_splitting_obstruction.py"]
        == file_sha256(TOKEN_VERIFIER)
        and hashes["results/spectral/token_splitting_obstruction.json"]
        == file_sha256(TOKEN_ARTIFACT)
        and hashes["sources/stacks_chevalley_054K.html"]
        == file_sha256(CHEVALLEY_SOURCE)
        and token_artifact["data"]["generic_finiteness"]["external_source"]
        == {
            "key": "stacks_chevalley_054K",
            "title": "Theorem 29.23.3 (054K): Chevalley's Theorem",
            "url": "https://stacks.math.columbia.edu/tag/054K",
            "local_path": "sources/stacks_chevalley_054K.html",
            "sha256": file_sha256(CHEVALLEY_SOURCE),
        },
    )

    stored = data["finite_graph_controls"]
    all_degrees = True
    all_modes = True
    all_determinants = True
    all_bounds = True
    for name, (order, edges) in controls().items():
        row = stored[name]
        edge_count = len(edges)
        dimension = 1 << order
        exponents = aligned_exponents(order, edges)
        diagonal_degree = max(
            edge_count - exponent + 2 * exponent for exponent in exponents
        )
        formula_diagonal = 2 * edge_count
        delta = 2 * order + formula_diagonal
        d = dimension * delta
        all_degrees &= (
            int(row["aligned_exponent_min"]) == min(exponents)
            and int(row["aligned_exponent_max"]) == max(exponents)
            and all(0 <= value <= edge_count for value in exponents)
            and diagonal_degree == formula_diagonal
            and int(row["diagonal_degree"]) == formula_diagonal
            and int(row["entry_degree_delta"]) == delta
            and int(row["incidence_equation_degree_d"]) == d
        )
        observed_mode_degrees = mode_degree_maxima(order)
        stored_mode_degrees = [
            int(item["maximum_total_mode_degree"])
            for item in row["mode_coefficient_degree_rows"]
        ]
        all_modes &= observed_mode_degrees == stored_mode_degrees and all(
            value <= index * (order + 1)
            for index, value in enumerate(observed_mode_degrees, start=1)
        )
        determinant = row["determinant_factor_exponents"]
        all_determinants &= (
            sum(exponents) == edge_count * dimension // 2
            == int(determinant["aligned_exponent_sum"])
            and int(determinant["power_2t"]) == edge_count * dimension // 2
            and int(determinant["power_1_plus_t2"])
            == edge_count * dimension // 2
            and int(determinant["power_1_minus_t2"]) == order * dimension
        )
        graph_bound = d ** (order + 2)
        uniform_delta = order * order + order
        uniform_bound = (dimension * uniform_delta) ** (order + 2)
        all_bounds &= (
            int(row["incidence_variety_degree_bound_decimal"]) == graph_bound
            and (
                row["exceptional_candidate_bound_decimal"] is None
                if name == "path4"
                else int(row["exceptional_candidate_bound_decimal"]) == graph_bound
            )
            and delta <= uniform_delta
            and int(row["uniform_simple_graph_bound_decimal"]) == uniform_bound
            and graph_bound <= uniform_bound
        )
    check("independent polynomial degree reconstruction", all_degrees)
    check("independent Boolean-mode degree reconstruction", all_modes)
    check("independent determinant factor exponents", all_determinants)
    check("independent affine and uniform bounds", all_bounds)

    open_grid = stored["open2x3"]
    check(
        "open 2x3 scale control",
        open_grid["order_n"] == 6
        and open_grid["edge_count_m"] == 7
        and open_grid["dimension_N"] == 64
        and open_grid["entry_degree_delta"] == 26
        and open_grid["incidence_equation_degree_d"] == 1664
        and open_grid["exceptional_candidate_bound_decimal"]
        == "58779593138084256579321856",
    )
    family = data["grid_family_corollaries"]
    check(
        "grid corollaries exclude path widths",
        family["open_2xL"]["range"]
        == "L>=2 (L=1 is a path and is outside the finite-projection theorem)"
        and family["open_3xL"]["range"]
        == "L>=2 (L=1 is a path and is outside the finite-projection theorem)",
    )

    transfer = data["transfer_polynomialization"]
    incidence = data["incidence_system"]
    lemma = data["affine_degree_lemma"]
    expected_finiteness = (
        "e236 endpoint token-splitting obstruction plus Chevalley constructibility "
        "makes the nonsingular t-projection finite; the five determinant roots "
        "are already a finite vertical set"
    )
    expected_positivity = (
        "At physical t the polynomialized matrix is positive definite. Any complex "
        "coefficient solution contains c and c*u_i among its positive eigenvalues, "
        "hence c,u_i are positive real."
    )
    expected_theorem = (
        "Let G be a finite connected simple nonpath graph with n vertices, m edges, "
        "and N=2^n. The t-image of the full polynomialized coefficient-incidence "
        "variety has at most [N*(2n+2m)]^(n+2) complex points. It contains every "
        "physical coupling whose spectrum is one full nonzero n-mode subset-product "
        "multiset, so the physical exceptional set inside 0<t<1 obeys the same "
        "bound. Uniformly over connected simple nonpath n-vertex graphs the bound "
        "is [2^n*(n^2+n)]^(n+2)."
    )
    check(
        "algebraic finiteness chain is explicit",
        data["theorem"] == expected_theorem
        and transfer["physical_scalar_equivalence"]
        == "R_unc=P_t*diag(q^a_sigma)*P_t, q=(1+t^2)/(2t), a_sigma=m-cut_G(sigma); R_unc differs from the centered physical matrix by one nonzero scalar"
        and transfer["clearing_scalar"] == "S(t)=(2t)^m"
        and transfer["diagonal_identity"]
        == "S*q^a=(2t)^(m-a)*(1+t^2)^a"
        and transfer["entry_degree_delta"] == "2n+2m"
        and transfer["determinant"]
        == "det(R_tilde)=(2t)^(mN/2)*(1+t^2)^(mN/2)*(1-t^2)^(nN)"
        and set(transfer["singular_parameter_roots"])
        == {"0", "+i", "-i", "+1", "-1"}
        and incidence["finiteness_input"] == expected_finiteness
        and incidence["physical_positivity"] == expected_positivity
        and incidence["equation_degree_bound"] == "d=N*(2n+2m)"
        and lemma["statement"]
        == "If an affine algebraic set V in A^M is cut out by polynomials of degree at most d, then deg(V)<=d^M. If one coordinate has finite image on V, each irreducible component maps to one point, so the number of image points is at most the number of components and hence at most deg(V)."
        and lemma["proof_sketch"]
        == "Homogenize the defining equations in projective M-space. Theorem 2.1 of arXiv:2511.07639v1 bounds the degree of their common projective zero locus by d^M independently of the number of equations. The affine components are exactly the projective components meeting the affine chart, so their total degree obeys the same bound. Every irreducible component has degree at least one, and a finite irreducible constructible image in A^1 is one point."
        and lemma["primary_source"]["doi"]
        == "10.1016/0304-3975(83)90002-6"
        and lemma["corroborating_source"]["arxiv"] == "1005.1858"
        and lemma["recent_primary_source"]["arxiv"] == "2511.07639v1"
        and lemma["recent_primary_source"]["url"]
        == "https://arxiv.org/abs/2511.07639v1",
    )
    scope = data["scope"]
    check(
        "scope preserves unresolved emptiness",
        scope["proved"]
        == [
            "an explicit finite complex-candidate bound for every fixed connected simple nonpath graph",
            "the uniform connected-simple-nonpath n-vertex bound [2^n*(n^2+n)]^(n+2)",
            "the same upper bound for physical full-spectrum exceptional couplings on bipartite grids with both side lengths at least two",
        ]
        and scope["not_proved"]
        == [
            "emptiness of the exceptional set beyond the already disposed open 2x3 layer",
            "a useful numerical bound, a list of exceptional roots, or a graph-uniform root-free interval",
            "parity-sector unions, larger auxiliary spectra, integrability, or a thermodynamic solution",
            "an improved critical endpoint or any use of the benchmark critical coupling",
        ],
    )
    check(
        "all producer checks passed",
        len(artifact["checks"]) == 9
        and all(bool(item["passed"]) for item in artifact["checks"]),
    )
    check(
        "producer resource metadata",
        meta["verifier"] == "tests/test_bipartite_exception_degree.py"
        and meta["artifact"] == "results/spectral/bipartite_exception_degree.json"
        and meta["process_cpu_budget_seconds"] == 120.0
        and meta["rss_limit_bytes"] == RSS_LIMIT_BYTES
        and meta["peak_rss_measurement"] == peak_rss_measurement()
        and meta["process_cpu_seconds"] < meta["process_cpu_budget_seconds"]
        and meta["peak_rss_bytes"] < meta["rss_limit_bytes"],
    )
    verifier_rss = peak_rss_bytes()
    check(
        "RSS wall",
        verifier_rss < RSS_LIMIT_BYTES,
        f"{verifier_rss}/{RSS_LIMIT_BYTES} bytes via {peak_rss_measurement()}",
    )
    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks failed: {', '.join(FAILURES)}")
        return 1
    print("PASS: graph-uniform full-spectrum exceptional-coupling degree bound")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
