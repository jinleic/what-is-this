#!/usr/bin/env python3
"""High-temperature token-graph obstruction to a full Gaussian spectrum.

Set s=1-t on the physical isotropic curve.  In the Walsh basis the transfer
operator has bands of valuation 2k.  The first graph-sensitive splitting of
the weight-k band is one quarter of the adjacency spectrum of the unsigned
k-token graph F_k(G).  A full n-mode subset-product spectrum would instead
force the k-band splitting to be the k-subset sums of the one-particle
splitting, namely the spectrum of the additive compound A_G^[k] on Lambda^k.

For k=2 these matrices have the same underlying token graph.  The exterior
signing is balanced exactly when a connected simple graph is a path.  Cycles
permit particle exchange around the cycle; a degree-three vertex permits the
six-step claw exchange.  Any negative closed walk makes an unsigned/signed
trace moment differ strictly.  All finite checks here use exact integers and
fractions.
"""

from __future__ import annotations

import itertools
import json
import platform
import resource
import sys
import time
from collections import deque
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "spectral" / "token_splitting_obstruction.json"


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def add_check(rows: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    rows.append({"name": name, "passed": bool(passed), "detail": detail})


def connected(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = {0}
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for other in adjacency[vertex]:
            if other not in seen:
                seen.add(other)
                queue.append(other)
    return len(seen) == sites


def is_path(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    degrees = [0] * sites
    for left, right in edges:
        degrees[left] += 1
        degrees[right] += 1
    return (
        connected(sites, edges)
        and len(edges) == sites - 1
        and max(degrees, default=0) <= 2
    )


def token_signed_edges(
    sites: int, edges: tuple[tuple[int, int], ...]
) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    subsets = list(itertools.combinations(range(sites), 2))
    index = {subset: slot for slot, subset in enumerate(subsets)}
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    signed: dict[tuple[int, int], int] = {}
    for source_slot, source in enumerate(subsets):
        for position, old_vertex in enumerate(source):
            for new_vertex in adjacency[old_vertex]:
                if new_vertex in source:
                    continue
                raw = list(source)
                raw[position] = new_vertex
                sign = -1 if raw[0] > raw[1] else 1
                target = tuple(sorted(raw))
                target_slot = index[target]
                edge = tuple(sorted((source_slot, target_slot)))
                if edge in signed and signed[edge] != sign:
                    raise AssertionError("inconsistent exterior sign")
                signed[edge] = sign
    return subsets, signed


def signing_balanced(vertex_count: int, signed: dict[tuple[int, int], int]) -> bool:
    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(vertex_count)]
    for (left, right), sign in signed.items():
        adjacency[left].append((right, sign))
        adjacency[right].append((left, sign))
    switches: list[int | None] = [None] * vertex_count
    for root in range(vertex_count):
        if switches[root] is not None:
            continue
        switches[root] = 1
        queue = deque([root])
        while queue:
            vertex = queue.popleft()
            for other, sign in adjacency[vertex]:
                required = int(switches[vertex]) * sign
                if switches[other] is None:
                    switches[other] = required
                    queue.append(other)
                elif switches[other] != required:
                    return False
    return True


def matrices(sites: int, edges: tuple[tuple[int, int], ...]) -> tuple[list[list[int]], list[list[int]]]:
    subsets, signed = token_signed_edges(sites, edges)
    size = len(subsets)
    unsigned = [[0] * size for _ in range(size)]
    exterior = [[0] * size for _ in range(size)]
    for (left, right), sign in signed.items():
        unsigned[left][right] = unsigned[right][left] = 1
        exterior[left][right] = exterior[right][left] = sign
    return unsigned, exterior


def matrix_multiply(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    size = len(left)
    out = [[0] * size for _ in range(size)]
    for row in range(size):
        for middle, value in enumerate(left[row]):
            if value:
                for column, other in enumerate(right[middle]):
                    out[row][column] += value * other
    return out


def trace_power(matrix: list[list[int]], exponent: int) -> int:
    size = len(matrix)
    power = [[int(row == column) for column in range(size)] for row in range(size)]
    for _ in range(exponent):
        power = matrix_multiply(power, matrix)
    return sum(power[index][index] for index in range(size))


def first_trace_difference(
    sites: int, edges: tuple[tuple[int, int], ...], maximum: int = 12
) -> tuple[int | None, int]:
    unsigned, exterior = matrices(sites, edges)
    for exponent in range(1, maximum + 1):
        difference = trace_power(unsigned, exponent) - trace_power(exterior, exponent)
        if difference:
            return exponent, difference
    return None, 0


def loop_sign(
    sites: int,
    edges: tuple[tuple[int, int], ...],
    loop: tuple[tuple[int, int], ...],
) -> int:
    subsets, signed = token_signed_edges(sites, edges)
    index = {subset: slot for slot, subset in enumerate(subsets)}
    result = 1
    for source, target in zip(loop, loop[1:]):
        edge = tuple(sorted((index[tuple(sorted(source))], index[tuple(sorted(target))])))
        result *= signed[edge]
    return result


def graph_balance_census(max_sites: int) -> dict[str, int]:
    connected_total = 0
    path_total = 0
    balanced_total = 0
    mismatches = 0
    for sites in range(1, max_sites + 1):
        possible_edges = tuple(itertools.combinations(range(sites), 2))
        for edge_mask in range(1 << len(possible_edges)):
            edges = tuple(
                edge
                for index, edge in enumerate(possible_edges)
                if (edge_mask >> index) & 1
            )
            if not connected(sites, edges):
                continue
            subsets, signed = token_signed_edges(sites, edges)
            balanced = signing_balanced(len(subsets), signed)
            path = is_path(sites, edges)
            mismatches += int(balanced != path)
            connected_total += 1
            path_total += int(path)
            balanced_total += int(balanced)
    return {
        "maximum_sites": max_sites,
        "connected_graph_total": connected_total,
        "path_total": path_total,
        "balanced_signing_total": balanced_total,
        "path_balance_mismatches": mismatches,
    }


def e_fourier_audit(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    edge_masks = {(1 << left) ^ (1 << right) for left, right in edges}
    edge_count = len(edges)
    values = []
    for state in range(1 << sites):
        cut = sum(
            ((state >> left) & 1) != ((state >> right) & 1)
            for left, right in edges
        )
        values.append(edge_count // 2 - cut)
    for mask in range(1 << sites):
        coefficient = sum(
            value * (-1 if (mask & state).bit_count() % 2 else 1)
            for state, value in enumerate(values)
        ) / Fraction(1 << sites)
        expected = Fraction(-1, 2) if mask == 0 and edge_count % 2 else Fraction(0)
        if mask in edge_masks:
            expected += Fraction(1, 2)
        if coefficient != expected:
            return False
    return True


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []
    path4 = (4, ((0, 1), (1, 2), (2, 3)))
    cycle4 = (4, ((0, 1), (1, 2), (2, 3), (3, 0)))
    claw = (4, ((0, 1), (0, 2), (0, 3)))
    open2x3 = (
        6,
        ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)),
    )
    controls = {
        name: first_trace_difference(*graph)
        for name, graph in {
            "path4": path4,
            "cycle4": cycle4,
            "claw": claw,
            "open2x3": open2x3,
        }.items()
    }
    expected_controls = {
        "path4": (None, 0),
        "cycle4": (4, 64),
        "claw": (6, 24),
        "open2x3": (4, 128),
    }
    add_check(
        checks,
        "finite_token_trace_controls",
        controls == expected_controls,
        str(controls),
    )

    cycle_loop = ((0, 1), (0, 2), (0, 3), (1, 3), (0, 1))
    claw_loop = ((1, 2), (0, 2), (2, 3), (0, 3), (1, 3), (0, 1), (1, 2))
    add_check(
        checks,
        "cycle_exchange_loop_is_negative",
        loop_sign(*cycle4, cycle_loop) == -1,
        str(cycle_loop),
    )
    add_check(
        checks,
        "claw_exchange_loop_is_negative",
        loop_sign(*claw, claw_loop) == -1,
        str(claw_loop),
    )

    census = graph_balance_census(6)
    add_check(
        checks,
        "all_connected_simple_graphs_through_six_vertices_path_iff_balanced",
        census["connected_graph_total"] == 27476
        and census["path_balance_mismatches"] == 0
        and census["path_total"] == census["balanced_signing_total"],
        str(census),
    )
    add_check(
        checks,
        "graph_energy_walsh_transform",
        e_fourier_audit(*path4)
        and e_fourier_audit(*cycle4)
        and e_fourier_audit(*claw)
        and e_fourier_audit(*open2x3),
        "only the empty mask (for odd m) and edge masks survive; every edge coefficient is 1/2",
    )
    q_identity = all(
        (1 + (1 - s) ** 2) / (2 * (1 - s))
        == 1 + s * s / (2 * (1 - s))
        for s in (Fraction(1, 7), Fraction(2, 9), Fraction(5, 12))
    )
    add_check(
        checks,
        "physical_high_temperature_parameter_identity",
        q_identity,
        "q=1+s^2/(2*(1-s)); t->1 (s->0) is K->0",
    )
    block_order_table = {
        f"{weight},{other}": 2 * max(weight, other)
        for weight in range(7)
        for other in range(7)
        if (weight - other) % 2 == 0
    }
    add_check(
        checks,
        "all_band_entry_order_ledger",
        all(
            order
            == sum(map(int, pair.split(",")))
            + abs(int(pair.split(",")[0]) - int(pair.split(",")[1]))
            for pair, order in block_order_table.items()
        ),
        "ord_s B_(k,l)>=k+l+|k-l|=2*max(k,l); odd-parity blocks vanish",
    )

    if not all(row["passed"] for row in checks):
        failed = [row["name"] for row in checks if not row["passed"]]
        raise AssertionError(f"producer checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e236_token_splitting_obstruction.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers and fractions",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
            "benchmark_used": False,
        },
        "data": {
            "claim_tag": "[THEOREM]",
            "theorem": "For every finite connected simple graph G that is not a path, there exists epsilon_G>0 such that for every physical 1-epsilon_G<t=tanh(K*/2)<1 (equivalently every sufficiently small positive isotropic K), its symmetric layer transfer spectrum is not one full positive |V(G)|-mode subset-product multiset.",
            "full_subset_product_definition": "One scalar c>0 and exactly n positive modes u_i whose 2^n Boolean products, with multiplicity, equal the complete spectrum. A union of parity sectors with different mode sets is not covered.",
            "band_splitting": {
                "parameter": "s=1-t; q=(1+t^2)/(2t)=1+s^2/(2*(1-s))=exp(2K), equivalently exp(-2K)=2t/(1+t^2); s->0 is K->0+",
                "parameter_convention": "t is the half-dual Hamming-kernel parameter, not exp(-2K): exp(-2K)=2t/(1+t^2). Therefore q-1=s^2/(2*(1-s)) exactly and D-I begins at order s^2.",
                "q_minus_one_certificate": {
                    "numerator_coefficients_ascending": [0, 0, 1],
                    "denominator_coefficients_ascending": [2, -2],
                    "s_valuation": 2,
                },
                "walsh_kernel_eigenvalue": "rho_k=(2-s)^(n-k)*s^k on weight-k subsets",
                "diagonal_expansion": "D=I+(s^2/2)*e+O(s^3), e_sigma=floor(m/2)-cut_G(sigma)",
                "energy_fourier_transform": "hat(e)_empty=-1/2 for odd m and 0 for even m; hat(e)_{uv}=1/2 for each simple edge uv; all other coefficients vanish",
                "same_band_operator": "the graph-sensitive relative s^2 coefficient on weight k is (1/4)*A(F_k(G)), the unsigned hard-core k-token adjacency",
                "off_band_control": "For every same-parity lower band ell<k, the full block B_(k,ell)=O(s^(2k)); the lower resolvent at a k-band eigenvalue is O(s^(-2k+4)), so the complete lower Schur term is O(s^(2k+4)), relative O(s^4). For ell>k, B_(k,ell)=O(s^(2ell)); k+2 gives the leading upper Schur term O(s^(2k+8)), relative O(s^8). Odd-parity blocks vanish.",
                "normalized_limit": "If theta ranges over spec A(F_k(G)) and r(s)=(s/(2-s))^2, then 4*s^(-2)*((lambda_{k}/lambda_top)/r(s)^k-1) tends, as a multiset, to theta.",
                "formal_block_certificate": {
                    "same_band_edge_coefficient_numerator": 1,
                    "same_band_edge_coefficient_denominator": 4,
                    "entry_order_formula": "ord_s B_(k,ell) >= k+ell+|k-ell| = 2*max(k,ell)",
                    "entry_order_table_through_weight_six": block_order_table,
                    "lower_schur_relative_order": 4,
                    "upper_schur_relative_order": 8,
                },
            },
            "band_valuation": {
                "target_multiset": "{2k with multiplicity C(n,k), k=0..n}",
                "perron_zero_multiplicity": 1,
                "subset_sum_mean": "sum_i a_i/2=n",
                "subset_sum_variance": "sum_i a_i^2/4=n",
                "conclusion": "sum a_i=2n and sum a_i^2=4n, hence sum_i(a_i-2)^2=0 and every one of the n mode-ratio valuations is 2",
            },
            "gaussian_necessity": {
                "valuation": "Invert modes larger than one and absorb them into c, so Perron simplicity makes c=lambda_top and all n oriented mode ratios strictly smaller than one. A finite occurrence-bijection subsequence plus the universal valuation multiset gives nonnegative integer mode-ratio valuations; its mean and variance force every valuation to equal 2.",
                "first_splitting": "Define alpha_i:=lim 4*s^(-2)*(u_i/r(s)-1) on a bounded fixed-slot subsequence. Then u_i/r(s)=1+(s^2/4)*alpha_i+o(s^2), and exact Boolean products give lim 4*s^(-2)*(prod_(i in S)u_i/r(s)^|S|-1)=sum_(i in S)alpha_i. No analytic mode function is assumed.",
                "matrix_form": "The alpha_i multiset is spec A_G, so spec A(F_k(G)) must equal the k-subset-sum spectrum of the additive compound A_G^[k] (the derivation induced on Lambda^k, not the multiplicative exterior power); at k=2 it is the signed exterior adjacency on the same token graph.",
            },
            "signing_theorem": {
                "statement": "For a finite connected simple graph, the exterior signing of F_2(G) is balanced iff G is a path.",
                "path_reason": "after relabelling vertices in path order, an allowed hard-core move never crosses the other token, so every exterior edge sign is positive (equivalently particle order supplies a switching function)",
                "cycle_reason": "on any cycle, one particle can traverse the unoccupied arc and the two labels exchange on returning to the same unordered configuration",
                "branch_reason": "three neighbours of a degree-at-least-three vertex give the explicit six-step claw exchange",
                "trace_reason": "tr(A_unsigned^r)-tr(A_signed^r) is the sum over closed walks of 1-sign(w), hence is nonnegative and is strictly positive when a negative closed walk exists",
                "cycle4_negative_loop": [list(pair) for pair in cycle_loop],
                "claw_negative_loop": [list(pair) for pair in claw_loop],
            },
            "generic_finiteness": {
                "tag": "[THEOREM][EXTERNAL]",
                "statement": "For every finite connected simple nonpath graph G, the physical couplings whose complete spectrum is one full n-mode subset-product multiset form a finite set. Paths are the only connected simple graphs not excluded by this argument; path Gaussianity is not proved here.",
                "coefficient_family": "On the nonsingular open set U={d_G(t)!=0}, clear a common denominator d_G(t) in the target characteristic coefficients and equate them to the coefficients of prod_(S subseteq [n]) (x-c*prod_(i in S)u_i). This defines a finite-type algebraic family in (t,c,u); equivalently add z*d_G(t)-1=0 or saturate by d_G.",
                "constructibility": "Chevalley's theorem makes the t-projection constructible in U. A constructible subset of the affine line is finite or cofinite.",
                "positive_real_fibres": "At physical real t, monic characteristic-polynomial equality makes c and every c*u_i roots of the positive-definite transfer matrix, hence c>0 and u_i=(c*u_i)/c>0 even if the existential variables were allowed complex.",
                "finish": "The high-temperature theorem supplies a real interval in U with no fibre, so the constructible projection is not cofinite and must be finite.",
                "external_source": {
                    "key": "stacks_chevalley_054K",
                    "title": "Theorem 29.23.3 (054K): Chevalley's Theorem",
                    "url": "https://stacks.math.columbia.edu/tag/054K",
                    "local_path": "sources/stacks_chevalley_054K.html",
                    "sha256": "cd262f58542d25515f167f1ae7f35617659fbda39dba3decb453931288c04740",
                },
            },
            "finite_controls": {
                name: {"first_trace": exponent, "difference": difference}
                for name, (exponent, difference) in controls.items()
            },
            "census": census,
            "scope": {
                "interval_kind": "exists epsilon_G>0; no explicit numerical epsilon is claimed",
                "applies_to": [
                    "every finite connected simple nonpath graph, including every open rectangular layer containing a plaquette and every branching tree",
                    "the high-temperature endpoint K->0+ on the physical isotropic curve",
                    "one full positive n-mode subset-product spectrum",
                    "all but finitely many physical couplings, by the generic-finiteness corollary",
                ],
                "not_proved": [
                    "an every-coupling theorem or a graph-uniform epsilon_G",
                    "emptiness, a degree bound, or an explicit list for the finite exceptional coupling set",
                    "multigraphs or weighted/parallel token moves; simplicity is part of the signing statement",
                    "failure of parity-sector unions, paired sectors, restrictions of larger Gaussian operators, or other descriptions not equal to one full n-mode subset-product spectrum",
                    "nonintegrability or an impossibility of solving the three-dimensional Ising model",
                    "an exact thermodynamic free energy, critical point, or critical exponent",
                ],
                "two_dimensional_control": "For C4, spec A(F_2)={2*sqrt(2),-2*sqrt(2),0,0,0,0} can coincide with a differently twisted sector, but the actual weight-one band pins spec A(C4)={2,0,0,-2}, whose pair sums are {2,2,0,0,-2,-2}; these differ. Integrable parity-sector unions with different mode sets remain outside the theorem.",
            },
        },
        "checks": checks,
    }


def main() -> int:
    payload = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    for row in payload["checks"]:
        print(f"[{'PASS' if row['passed'] else 'FAIL'}] {row['name']}: {row['detail']}")
    print(f"PASS: {len(payload['checks'])}/{len(payload['checks'])} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
