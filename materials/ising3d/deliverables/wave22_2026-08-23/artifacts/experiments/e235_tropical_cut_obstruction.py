#!/usr/bin/env python3
"""Tropical low-temperature obstruction for every cyclic simple layer graph.

For the symmetric rational layer representative R_G(t)=P_t D_G(t) P_t,
congruence min-max bounds show that its ordered eigenvalue valuations at t=0
are exactly the shifted cut sizes of G.  A full positive n-mode subset-product
spectrum would make that valuation multiset an affine Boolean subset-sum.

The cut distribution of a simple graph has mean |E|/2 and variance |E|/4.
The zero-cut multiplicity fixes the number of zero mode-valuations.  Positive
integer mode-valuations a_i would then satisfy sum a_i=sum a_i^2=|E|, forcing
all a_i=1 and |E|=|V|-components.  Thus only forests pass this necessary
tropical criterion.  No floating-point arithmetic is used.
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
OUTPUT = ROOT / "results" / "spectral" / "tropical_cut_obstruction.json"


def max_rss_bytes() -> int:
    raw = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return raw if platform.system() == "Darwin" else raw * 1024


def add_check(rows: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    rows.append({"name": name, "passed": bool(passed), "detail": detail})


def component_count(sites: int, edges: tuple[tuple[int, int], ...]) -> int:
    adjacency = [[] for _ in range(sites)]
    for left, right in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    seen = [False] * sites
    components = 0
    for root in range(sites):
        if seen[root]:
            continue
        components += 1
        seen[root] = True
        queue = deque([root])
        while queue:
            vertex = queue.popleft()
            for other in adjacency[vertex]:
                if not seen[other]:
                    seen[other] = True
                    queue.append(other)
    return components


def cut_counts(sites: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    counts = [0] * (len(edges) + 1)
    edge_masks = tuple((1 << left) ^ (1 << right) for left, right in edges)
    for state in range(1 << sites):
        cut_size = sum((state & mask).bit_count() == 1 for mask in edge_masks)
        counts[cut_size] += 1
    return counts


def cut_moments(counts: list[int]) -> tuple[Fraction, Fraction]:
    total = sum(counts)
    mean = sum(degree * count for degree, count in enumerate(counts)) / Fraction(total)
    second = sum(degree * degree * count for degree, count in enumerate(counts)) / Fraction(total)
    return mean, second - mean * mean


def forest_polynomial(sites: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    from math import comb

    components = component_count(sites, edges)
    return [(1 << components) * comb(len(edges), degree) for degree in range(len(edges) + 1)]


def exact_graph_census(max_sites: int) -> dict[str, int]:
    graph_total = 0
    forest_total = 0
    cyclic_total = 0
    moment_failures = 0
    factorization_mismatches = 0
    for sites in range(1, max_sites + 1):
        possible_edges = tuple(itertools.combinations(range(sites), 2))
        for edge_mask in range(1 << len(possible_edges)):
            edges = tuple(
                edge
                for index, edge in enumerate(possible_edges)
                if (edge_mask >> index) & 1
            )
            components = component_count(sites, edges)
            forest = len(edges) == sites - components
            counts = cut_counts(sites, edges)
            mean, variance = cut_moments(counts)
            moment_failures += int(
                mean != Fraction(len(edges), 2) or variance != Fraction(len(edges), 4)
            )
            factorization_mismatches += int(
                (counts == forest_polynomial(sites, edges)) != forest
            )
            graph_total += 1
            forest_total += int(forest)
            cyclic_total += int(not forest)
    return {
        "maximum_sites": max_sites,
        "graph_total": graph_total,
        "forest_total": forest_total,
        "cyclic_total": cyclic_total,
        "moment_failures": moment_failures,
        "forest_factorization_mismatches": factorization_mismatches,
    }


def hamming_kernel_audit(sites: int, t: Fraction) -> bool:
    dimension = 1 << sites
    matrix = [
        [t ** ((row ^ column).bit_count()) for column in range(dimension)]
        for row in range(dimension)
    ]
    for character in range(dimension):
        vector = [
            -1 if (character & state).bit_count() % 2 else 1
            for state in range(dimension)
        ]
        weight = character.bit_count()
        expected = (1 + t) ** (sites - weight) * (1 - t) ** weight
        image = [
            sum(
                (matrix[row][column] * vector[column] for column in range(dimension)),
                start=Fraction(0),
            )
            for row in range(dimension)
        ]
        if image != [expected * value for value in vector]:
            return False
    return True


def run() -> dict[str, object]:
    started = time.process_time()
    checks: list[dict[str, object]] = []

    path4 = (4, ((0, 1), (1, 2), (2, 3)))
    star4 = (4, ((0, 1), (0, 2), (0, 3)))
    cycle4 = (4, ((0, 1), (1, 2), (2, 3), (3, 0)))
    open2x3 = (
        6,
        ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)),
    )
    double_edge = (2, ((0, 1), (0, 1)))
    controls = {
        "path4": cut_counts(*path4),
        "cycle4": cut_counts(*cycle4),
        "open2x3": cut_counts(*open2x3),
    }
    add_check(
        checks,
        "nonisomorphic_tree_controls_share_binomial_cut_polynomial",
        cut_counts(*path4) == cut_counts(*star4) == [2, 6, 6, 2],
        "both four-vertex trees give 2*(1+z)^3",
    )
    add_check(
        checks,
        "cycle_and_grid_controls_fail_forest_polynomial",
        controls["cycle4"] == [2, 0, 12, 0, 2]
        and controls["open2x3"] == [2, 0, 12, 18, 18, 12, 0, 2],
        "each cyclic graph differs from 2*(1+z)^(n-1)",
    )
    double_counts = cut_counts(*double_edge)
    double_mean, double_variance = cut_moments(double_counts)
    add_check(
        checks,
        "parallel_edge_control_shows_simplicity_is_load_bearing",
        double_counts == [2, 0, 2]
        and double_mean == 1
        and double_variance == 1
        and double_variance != Fraction(len(double_edge[1]), 4),
        "a doubled edge has cut polynomial 2*(1+z^2), so the simple-edge variance formula fails",
    )

    for name, graph in {
        "path4": path4,
        "star4": star4,
        "cycle4": cycle4,
        "open2x3": open2x3,
    }.items():
        counts = cut_counts(*graph)
        mean, variance = cut_moments(counts)
        edges = len(graph[1])
        add_check(
            checks,
            f"{name}_cut_moments",
            mean == Fraction(edges, 2) and variance == Fraction(edges, 4),
            f"mean={mean}, variance={variance}",
        )

    census = exact_graph_census(6)
    add_check(
        checks,
        "all_simple_graphs_through_six_vertices_obey_sharp_dichotomy",
        census["graph_total"] == 33867
        and census["moment_failures"] == 0
        and census["forest_factorization_mismatches"] == 0,
        str(census),
    )
    add_check(
        checks,
        "hamming_kernel_walsh_eigenvalues_exact",
        hamming_kernel_audit(4, Fraction(3, 11)),
        "P_t has eigenvalues (1+t)^(n-k)*(1-t)^k on Walsh weight k",
    )

    if not all(row["passed"] for row in checks):
        failed = [row["name"] for row in checks if not row["passed"]]
        raise AssertionError(f"producer checks failed: {failed}")

    return {
        "meta": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "producer": "experiments/e235_tropical_cut_obstruction.py",
            "interpreter": sys.executable,
            "arithmetic": "exact integers and fractions",
            "process_cpu_seconds": time.process_time() - started,
            "peak_rss_bytes": max_rss_bytes(),
            "benchmark_used": False,
        },
        "data": {
            "claim_tag": "[THEOREM]",
            "theorem": "For every finite simple graph G containing a cycle, there exists delta_G>0 such that for every physical 0<t=tanh(K*/2)<delta_G (equivalently every sufficiently large isotropic K) the spectrum of its symmetric layer transfer operator is not a full positive |V(G)|-mode subset-product multiset.",
            "full_subset_product_definition": "One global scalar c>0 and exactly n positive modes u_i such that the complete 2^n eigenvalue multiset is {c*prod_(i in S) u_i : S subseteq [n]}. Sector-wise unions using different mode sets are outside the hypothesis.",
            "matrix_valuation": {
                "representative": "R_G(t)=P_t D_G(t) P_t, P_t[sigma,tau]=t^d_H(sigma,tau)",
                "physical_curve": "q=(1+t^2)/(2t)=exp(2K), equivalently exp(-2K)=2t/(1+t^2)=tanh(K*); therefore t->0 is the low-temperature/strong-coupling limit K->infinity",
                "hamming_kernel_spectrum": "spec(P_t)={(1+t)^(n-k)*(1-t)^k with multiplicity C(n,k)}",
                "congruence": "PDP and D^(1/2) P^2 D^(1/2) have the same eigenvalues",
                "minmax_bounds": "with both spectra ordered increasingly and multiplicity retained: (1-t)^(2n)*lambda_i(D) <= lambda_i(R) <= (1+t)^(2n)*lambda_i(D)",
                "minmax_proof": "PDP=XX^T for X=P*D^(1/2), while X^T X=D^(1/2)P^2D^(1/2). The Walsh bounds (1-t)^(2n) I <= P^2 <= (1+t)^(2n) I survive congruence by D^(1/2), and Courant-Fischer preserves the sorted-index inequalities.",
                "diagonal": "For m edges, energy_sigma=m-2*cut_G(sigma) and epsilon=m mod 2, so e_sigma=(energy_sigma-epsilon)/2=floor(m/2)-cut_G(sigma), and D_sigma=q^e_sigma.",
                "conclusion": "after a common shift, the ordered spectral valuations at t=0 are exactly the cut sizes of G with multiplicity",
            },
            "subsequence_argument": {
                "steps": [
                    "If Gaussian points accumulated at zero, choose a sequence t_j->0; invert modes as needed so u_i(t_j)>=1, label mode bits, and label repeated eigenvalue occurrences.",
                    "There are finitely many occurrence-bijections from labelled Boolean slots to the 2^n increasingly ordered eigenvalue occurrence slots, so one orientation, labelling, and bijection pi is constant on a subsequence.",
                    "On that subsequence u_i(t_j)=lambda_(pi({i}))(t_j)/lambda_(pi(empty))(t_j) exactly. Defining nu(x_j)=lim_j log(x_j)/log(t_j), the ordered-slot bounds give nu(u_i)=nu(lambda_pi({i}))-nu(lambda_pi(empty)), a nonpositive integer; no algebraic or Puiseux mode function is assumed.",
                    "Taking valuation limits in every exact product identity gives nu(lambda_pi(S))=nu(lambda_pi(empty))+sum_(i in S)nu(u_i). After shifting the minimum, the cut-size multiset is the subset-sum multiset of nonnegative integer weights a_i; its maximum is both maxcut(G) and sum_i a_i.",
                ],
                "logical_finish": "No accumulating sequence implies the existence of a punctured interval (0,delta_G) containing no Gaussian point; delta_G is not made effective here.",
            },
            "moment_obstruction": {
                "zero_cut_multiplicity": "2^components",
                "zero_valuation_weights": "components",
                "positive_valuation_weight_count": "n-components",
                "cut_mean": "m/2",
                "cut_variance": "m/4",
                "sum_weights": "m",
                "sum_weight_squares": "m",
                "conclusion": "all positive valuation weights are 1, hence m=n-components",
                "graph_conclusion": "a finite simple graph satisfies m=n-components iff it is a forest",
                "forest_realization": "For a forest, the cut map has rank n-components=m, every edge-cut pattern has 2^components preimages, and the cut polynomial is 2^components*(1+z)^m: components zero weights and m unit weights. Passing this tropical criterion does not prove the full transfer spectrum Gaussian.",
            },
            "finite_controls": controls,
            "multigraph_boundary_control": {
                "graph": "two vertices joined by two parallel edges",
                "cut_counts": double_counts,
                "mean": str(double_mean),
                "variance": str(double_variance),
                "meaning": "simplicity is essential: parallel cut indicators coincide, and the valuation multiset 2*(1+z^2) is itself Boolean-subset-sum",
            },
            "census": census,
            "scope": {
                "interval_kind": "exists delta_G>0; no explicit numerical delta is claimed",
                "applies_to": [
                    "every finite simple cyclic graph, including every open rectangular layer containing a plaquette",
                    "bipartite and non-bipartite graphs alike",
                    "the physical positive isotropic transfer spectrum up to an irrelevant scalar",
                ],
                "not_proved": [
                    "an every-coupling theorem or a graph-uniform value of delta_G",
                    "multigraphs or periodic size-two directions with parallel bonds; a doubled edge is an exact counterexample to extending the theorem beyond simple graphs",
                    "failure of parity-projected, paired-sector, or larger-auxiliary Gaussian descriptions not equal to one full n-mode subset-product spectrum",
                    "that forests are Gaussian; forests only pass this leading-valuation necessity",
                    "nonintegrability or an impossibility of solving the three-dimensional Ising model",
                    "an exact thermodynamic free energy, critical point, or critical exponent",
                ],
            },
            "two_dimensional_control": "C4 fails one global four-mode valuation cube but remains compatible with the standard integrable union of parity sectors using different mode sets; sector-wise descriptions are outside the theorem.",
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
