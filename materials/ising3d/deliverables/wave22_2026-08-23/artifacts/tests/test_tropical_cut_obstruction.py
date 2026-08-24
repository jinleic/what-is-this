#!/usr/bin/env python3
"""Independent exact verifier for the tropical cut-spectrum obstruction.

No producer imports.  The verifier rebuilds cut distributions, exhausts all
simple graphs through five vertices for the forest/cycle moment dichotomy, and
checks the exact Walsh diagonalization of the Hamming kernel P_t.
"""

from __future__ import annotations

import itertools
import json
from collections import deque
from fractions import Fraction
import sympy as sp
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "tropical_cut_obstruction.json"
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


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
    for state in range(1 << sites):
        size = sum(((state >> left) & 1) != ((state >> right) & 1) for left, right in edges)
        counts[size] += 1
    return counts


def binomial_cut_counts(sites: int, edges: tuple[tuple[int, int], ...]) -> list[int]:
    from math import comb

    components = component_count(sites, edges)
    return [(1 << components) * comb(len(edges), degree) for degree in range(len(edges) + 1)]


def moments(counts: list[int]) -> tuple[Fraction, Fraction]:
    total = sum(counts)
    mean = sum(index * value for index, value in enumerate(counts)) / Fraction(total)
    second = sum(index * index * value for index, value in enumerate(counts)) / Fraction(total)
    return mean, second - mean * mean


def matvec(matrix: list[list[Fraction]], vector: list[int]) -> list[Fraction]:
    return [sum((entry * value for entry, value in zip(row, vector, strict=True)), Fraction()) for row in matrix]


def hamming_kernel_walsh_audit(sites: int, t: Fraction) -> bool:
    dimension = 1 << sites
    matrix = [
        [t ** ((row ^ column).bit_count()) for column in range(dimension)]
        for row in range(dimension)
    ]
    for mask in range(dimension):
        vector = [(-1) ** ((mask & state).bit_count()) for state in range(dimension)]
        weight = mask.bit_count()
        eigenvalue = (1 + t) ** (sites - weight) * (1 - t) ** weight
        if matvec(matrix, vector) != [eigenvalue * value for value in vector]:
            return False
    return True
def roots_below_with_multiplicity(poly: sp.Poly, threshold: sp.Rational) -> int:
    return sum(
        int(sp.Poly(factor).count_roots(-sp.oo, threshold)) * multiplicity
        for factor, multiplicity in sp.factor_list(poly.as_expr())[1]
    )


def transfer_bridge_audit(
    sites: int, edges: tuple[tuple[int, int], ...], t: Fraction
) -> bool:
    tr = sp.Rational(t.numerator, t.denominator)
    q = (1 + tr * tr) / (2 * tr)
    dimension = 1 << sites
    p_matrix = sp.Matrix(
        dimension,
        dimension,
        lambda row, column: tr ** ((int(row) ^ int(column)).bit_count()),
    )
    edge_count = len(edges)
    parity = edge_count % 2
    exponents: list[int] = []
    for state in range(dimension):
        spins = [1 if (state >> site) & 1 else -1 for site in range(sites)]
        energy = sum(spins[left] * spins[right] for left, right in edges)
        cut = sum(spins[left] != spins[right] for left, right in edges)
        if energy != edge_count - 2 * cut:
            return False
        exponent = (energy - parity) // 2
        if exponent != edge_count // 2 - cut:
            return False
        exponents.append(exponent)
    diagonal_entries = [q**exponent for exponent in exponents]
    diagonal = sp.diag(*diagonal_entries)
    transfer = p_matrix * diagonal * p_matrix
    similar = diagonal * p_matrix * p_matrix
    if transfer.charpoly().as_expr() != similar.charpoly().as_expr():
        return False
    polynomial = sp.Poly(transfer.charpoly().as_expr())
    lower_scale = (1 - tr) ** (2 * sites)
    upper_scale = (1 + tr) ** (2 * sites)
    for index, diagonal_value in enumerate(sorted(diagonal_entries)):
        lower = lower_scale * diagonal_value
        upper = upper_scale * diagonal_value
        if polynomial.eval(lower) == 0 or polynomial.eval(upper) == 0:
            return False
        if roots_below_with_multiplicity(polynomial, lower) > index:
            return False
        if roots_below_with_multiplicity(polynomial, upper) < index + 1:
            return False
    return True


def subset_sum_counts(weights: tuple[int, ...]) -> list[int]:
    maximum = sum(weights)
    counts = [0] * (maximum + 1)
    for mask in range(1 << len(weights)):
        counts[sum(weight for index, weight in enumerate(weights) if (mask >> index) & 1)] += 1
    return counts




def exhaustive_graph_audit(max_sites: int) -> tuple[int, int]:
    graphs = 0
    forests = 0
    for sites in range(1, max_sites + 1):
        possible = tuple(itertools.combinations(range(sites), 2))
        for mask in range(1 << len(possible)):
            edges = tuple(edge for index, edge in enumerate(possible) if (mask >> index) & 1)
            components = component_count(sites, edges)
            is_forest = len(edges) == sites - components
            counts = cut_counts(sites, edges)
            mean, variance = moments(counts)
            if mean != Fraction(len(edges), 2) or variance != Fraction(len(edges), 4):
                return -1, -1
            if (counts == binomial_cut_counts(sites, edges)) != is_forest:
                return -1, -1
            graphs += 1
            forests += int(is_forest)
    return graphs, forests


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact has meta data checks", {"meta", "data", "checks"} <= set(artifact))
    names = [row["name"] for row in artifact["checks"]]
    check("producer checks unique and true", len(names) == len(set(names)) and all(row["passed"] for row in artifact["checks"]))

    path = (4, ((0, 1), (1, 2), (2, 3)))
    cycle = (4, ((0, 1), (1, 2), (2, 3), (3, 0)))
    grid = (6, ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)))
    controls = {
        "path4": cut_counts(*path),
        "cycle4": cut_counts(*cycle),
        "open2x3": cut_counts(*grid),
    }
    check("tree cut polynomial is Boolean-binomial", controls["path4"] == [2, 6, 6, 2])
    check("cycle cut polynomial differs", controls["cycle4"] == [2, 0, 12, 0, 2])
    check("open 2x3 cut polynomial rebuilt", controls["open2x3"] == [2, 0, 12, 18, 18, 12, 0, 2])
    check("stored cut controls match", artifact["data"]["finite_controls"] == controls)
    double_counts = cut_counts(2, ((0, 1), (0, 1)))
    double_mean, double_variance = moments(double_counts)
    check(
        "parallel-edge boundary control blocks multigraph overclaim",
        double_counts == [2, 0, 2]
        and double_mean == 1
        and double_variance == 1
        and artifact["data"]["multigraph_boundary_control"]["cut_counts"]
        == double_counts,
    )

    graph_count, forest_count = exhaustive_graph_audit(5)
    check("all simple graphs through n=5 satisfy sharp forest dichotomy", graph_count == 1099 and forest_count > 0, f"graphs={graph_count}, forests={forest_count}")
    check("Hamming kernel exact Walsh spectrum", hamming_kernel_walsh_audit(4, Fraction(2, 7)))
    check(
        "even-edge transfer exponent and ordered minmax bridge",
        transfer_bridge_audit(3, ((0, 1), (1, 2)), Fraction(1, 3)),
    )
    check(
        "odd-edge transfer exponent and ordered minmax bridge",
        transfer_bridge_audit(
            3, ((0, 1), (1, 2), (2, 0)), Fraction(1, 3)
        ),
    )
    t_values = (Fraction(1, 2), Fraction(1, 3), Fraction(1, 4))
    q_values = tuple((1 + t * t) / (2 * t) for t in t_values)
    check(
        "physical corner is low temperature",
        q_values[0] < q_values[1] < q_values[2]
        and all(1 / q == 2 * t / (1 + t * t) for t, q in zip(t_values, q_values, strict=True)),
        str(q_values),
    )
    check(
        "forest zero-unit valuation control",
        subset_sum_counts((0, 1, 1, 1)) == controls["path4"],
    )
    check(
        "C4 has no one-global-four-mode valuation cube",
        not any(
            subset_sum_counts(weights) == controls["cycle4"]
            for weights in itertools.combinations_with_replacement(range(5), 4)
        )
        and "parity sectors" in artifact["data"]["two_dimensional_control"],
    )
    check(
        "doubled edge realizes excluded multigraph cube",
        subset_sum_counts((0, 2)) == double_counts,
    )

    proof = artifact["data"]["moment_obstruction"]
    check(
        "stored moment proof has exact sum identities",
        proof["sum_weights"] == "m"
        and proof["sum_weight_squares"] == "m"
        and proof["conclusion"] == "all positive valuation weights are 1, hence m=n-components",
    )
    check(
        "full subset-product hypothesis is explicit",
        "one global scalar" in artifact["data"]["full_subset_product_definition"].lower()
        and "sector-wise" in artifact["data"]["full_subset_product_definition"].lower(),
    )
    scope = artifact["data"]["scope"]
    check(
        "scope is existential near zero and full-spectrum only",
        scope["interval_kind"] == "exists delta_G>0; no explicit numerical delta is claimed"
        and any("full" in item.lower() and "spectrum" in item.lower() for item in scope["not_proved"])
        and any("multigraph" in item.lower() for item in scope["not_proved"]),
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks")
        return 1
    print("PASS: independent tropical cut-spectrum verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
