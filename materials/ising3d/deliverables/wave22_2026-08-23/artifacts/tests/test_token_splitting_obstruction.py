#!/usr/bin/env python3
"""Independent exact verifier for the high-temperature token-splitting theorem.

No producer imports.  The verifier constructs unsigned two-token and signed
exterior-square graphs, checks their closed-walk traces, and exhausts every
connected simple graph through five vertices for the path/balance dichotomy.
"""

from __future__ import annotations

import itertools
import hashlib
import json
from collections import deque
import sympy as sp
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "spectral" / "token_splitting_obstruction.json"
FAILURES: list[str] = []


def check(name: str, passed: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not passed:
        FAILURES.append(name)


def connected(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    if sites == 0:
        return False
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
    return connected(sites, edges) and len(edges) == sites - 1 and max(degrees, default=0) <= 2


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
        for position, old in enumerate(source):
            for new_vertex in adjacency[old]:
                if new_vertex in source:
                    continue
                raw = list(source)
                raw[position] = new_vertex
                inversions = int(raw[0] > raw[1])
                target = tuple(sorted(raw))
                target_slot = index[target]
                key = tuple(sorted((source_slot, target_slot)))
                sign = -1 if inversions else 1
                if key in signed and signed[key] != sign:
                    raise AssertionError("exterior sign is orientation-inconsistent")
                signed[key] = sign
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


def token_matrices(sites: int, edges: tuple[tuple[int, int], ...]) -> tuple[list[list[int]], list[list[int]]]:
    subsets, signed = token_signed_edges(sites, edges)
    size = len(subsets)
    unsigned = [[0] * size for _ in range(size)]
    exterior = [[0] * size for _ in range(size)]
    for (left, right), sign in signed.items():
        unsigned[left][right] = unsigned[right][left] = 1
        exterior[left][right] = exterior[right][left] = sign
    return unsigned, exterior


def multiply(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
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
        power = multiply(power, matrix)
    return sum(power[index][index] for index in range(size))


def first_trace_difference(sites: int, edges: tuple[tuple[int, int], ...], maximum: int = 12) -> tuple[int | None, int]:
    unsigned, exterior = token_matrices(sites, edges)
    for exponent in range(1, maximum + 1):
        difference = trace_power(unsigned, exponent) - trace_power(exterior, exponent)
        if difference:
            return exponent, difference
    return None, 0
def subset_sum_counts(weights: tuple[int, ...]) -> list[int]:
    counts = [0] * (sum(weights) + 1)
    for mask in range(1 << len(weights)):
        value = sum(
            weight for index, weight in enumerate(weights) if (mask >> index) & 1
        )
        counts[value] += 1
    return counts




def graph_balance_census(max_sites: int) -> tuple[int, int]:
    connected_total = 0
    mismatches = 0
    for sites in range(1, max_sites + 1):
        possible = tuple(itertools.combinations(range(sites), 2))
        for mask in range(1 << len(possible)):
            edges = tuple(edge for index, edge in enumerate(possible) if (mask >> index) & 1)
            if not connected(sites, edges):
                continue
            subsets, signed = token_signed_edges(sites, edges)
            balanced = signing_balanced(len(subsets), signed)
            mismatches += int(balanced != is_path(sites, edges))
            connected_total += 1
    return connected_total, mismatches


def fourier_e_audit(sites: int, edges: tuple[tuple[int, int], ...]) -> bool:
    edge_masks = {(1 << left) ^ (1 << right) for left, right in edges}
    m = len(edges)
    values = []
    for state in range(1 << sites):
        cut = sum(
            ((state >> left) & 1) != ((state >> right) & 1)
            for left, right in edges
        )
        values.append(m // 2 - cut)
    for mask in range(1 << sites):
        coefficient = sum(
            value * (-1 if (mask & state).bit_count() % 2 else 1)
            for state, value in enumerate(values)
        ) / Fraction(1 << sites)
        expected = Fraction(-1, 2) if mask == 0 and m % 2 else Fraction(0)
        if mask in edge_masks:
            expected += Fraction(1, 2)
        if coefficient != expected:
            return False
    return True


def hamming_kernel_walsh_audit(sites: int, t: Fraction) -> bool:
    dimension = 1 << sites
    matrix = [
        [t ** ((row ^ column).bit_count()) for column in range(dimension)]
        for row in range(dimension)
    ]
    for mask in range(dimension):
        vector = [
            -1 if (mask & state).bit_count() % 2 else 1
            for state in range(dimension)
        ]
        expected = (1 + t) ** (sites - mask.bit_count()) * (
            1 - t
        ) ** mask.bit_count()
        image = [
            sum(
                (
                    matrix[row][column] * vector[column]
                    for column in range(dimension)
                ),
                start=Fraction(0),
            )
            for row in range(dimension)
        ]
        if image != [expected * value for value in vector]:
            return False
    return True


def main() -> int:
    artifact = json.loads(ARTIFACT.read_text())
    check("artifact has meta data checks", {"meta", "data", "checks"} <= set(artifact))
    names = [row["name"] for row in artifact["checks"]]
    check("producer checks unique and true", len(names) == len(set(names)) and all(row["passed"] for row in artifact["checks"]))

    path4 = (4, ((0, 1), (1, 2), (2, 3)))
    cycle4 = (4, ((0, 1), (1, 2), (2, 3), (3, 0)))
    claw = (4, ((0, 1), (0, 2), (0, 3)))
    grid = (6, ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)))
    differences = {
        "path4": first_trace_difference(*path4),
        "cycle4": first_trace_difference(*cycle4),
        "claw": first_trace_difference(*claw),
        "open2x3": first_trace_difference(*grid),
    }
    check("path control has balanced exterior signing", differences["path4"] == (None, 0))
    check("cycle first negative closed walks", differences["cycle4"] == (4, 64))
    check("claw exchange first appears at length six", differences["claw"] == (6, 24))
    check("open 2x3 splitting differs at trace four", differences["open2x3"] == (4, 128))
    stored = artifact["data"]["finite_controls"]
    check(
        "stored token trace controls match",
        stored == {name: {"first_trace": exponent, "difference": difference} for name, (exponent, difference) in differences.items()},
    )

    connected_total, mismatches = graph_balance_census(5)
    check(
        "all connected simple graphs through n=5 obey path iff balanced",
        connected_total == 772 and mismatches == 0,
        f"connected={connected_total}, mismatches={mismatches}",
    )
    check("graph-sensitive Walsh coefficient rebuilt", fourier_e_audit(*grid))
    check(
        "Hamming kernel Walsh bands rebuilt",
        hamming_kernel_walsh_audit(4, Fraction(3, 8)),
    )
    formal = artifact["data"]["band_splitting"]["formal_block_certificate"]
    expected_orders = {
        f"{weight},{other}": 2 * max(weight, other)
        for weight in range(7)
        for other in range(7)
        if (weight - other) % 2 == 0
    }
    check(
        "same-band coefficient and all-band order ledger rebuilt",
        Fraction(
            formal["same_band_edge_coefficient_numerator"],
            formal["same_band_edge_coefficient_denominator"],
        )
        == Fraction(1, 4)
        and formal["entry_order_table_through_weight_six"] == expected_orders
        and formal["lower_schur_relative_order"] == 4
        and formal["upper_schur_relative_order"] == 8,
    )
    unsigned_c4, exterior_c4 = token_matrices(*cycle4)
    variable = sp.symbols("lambda")
    unsigned_charpoly = sp.factor(
        sp.Matrix(unsigned_c4).charpoly(variable).as_expr()
    )
    exterior_charpoly = sp.factor(
        sp.Matrix(exterior_c4).charpoly(variable).as_expr()
    )
    check(
        "C4 actual-band additive compound mismatch exact",
        sp.expand(unsigned_charpoly - variable**4 * (variable**2 - 8)) == 0
        and sp.expand(
            exterior_charpoly - variable**2 * (variable**2 - 4) ** 2
        )
        == 0,
        f"unsigned={unsigned_charpoly}; additive={exterior_charpoly}",
    )
    z = sp.symbols("z")
    alphas = (sp.Integer(2), sp.Integer(-1), sp.Integer(3))
    check(
        "first correction of Boolean products is additive",
        all(
            sp.expand(
                sp.prod(1 + z * alphas[index] / 4 for index in subset)
            ).coeff(z, 1)
            == sum(
                (alphas[index] for index in subset), start=sp.Integer(0)
            )
            / 4
            for subset_size in range(4)
            for subset in itertools.combinations(range(3), subset_size)
        ),
    )
    from math import comb

    valuation_target = [0] * 9
    for weight in range(5):
        valuation_target[2 * weight] = comb(4, weight)
    matching_weights = [
        weights
        for weights in itertools.combinations_with_replacement(range(9), 4)
        if subset_sum_counts(weights) == valuation_target
    ]
    check(
        "universal band valuations force every mode valuation two",
        matching_weights == [(2, 2, 2, 2)]
        and "sum_i(a_i-2)^2=0"
        in artifact["data"]["band_valuation"]["conclusion"],
        str(matching_weights),
    )
    for numerator in (1, 2, 5):
        s = Fraction(numerator, 11)
        t = 1 - s
        q = (1 + t * t) / (2 * t)
        check(f"exact q expansion s={s}", q == 1 + s * s / (2 * (1 - s)))
    q_certificate = artifact["data"]["band_splitting"][
        "q_minus_one_certificate"
    ]
    check(
        "half-dual parameter makes graph splitting quadratic in s",
        q_certificate
        == {
            "numerator_coefficients_ascending": [0, 0, 1],
            "denominator_coefficients_ascending": [2, -2],
            "s_valuation": 2,
        }
        and "not exp(-2K)"
        in artifact["data"]["band_splitting"]["parameter_convention"],
    )
    check(
        "Gaussian necessity uses additive compound and actual band one",
        "additive compound"
        in artifact["data"]["gaussian_necessity"]["matrix_form"]
        and "actual weight-one band"
        in artifact["data"]["scope"]["two_dimensional_control"],
    )
    generic = artifact["data"]["generic_finiteness"]
    source = ROOT / generic["external_source"]["local_path"]
    source_bytes = source.read_bytes()
    check(
        "Chevalley source archived and hash verified",
        hashlib.sha256(source_bytes).hexdigest()
        == generic["external_source"]["sha256"]
        and b"Chevalley" in source_bytes
        and b"constructible" in source_bytes,
    )
    check(
        "generic finite-exception corollary is scoped",
        "finite set" in generic["statement"]
        and "not proved" in generic["statement"].lower()
        and "saturate" in generic["coefficient_family"]
        and "not cofinite" in generic["finish"],
    )

    scope = artifact["data"]["scope"]
    check(
        "scope is high-temperature endpoint and full-spectrum only",
        scope["interval_kind"] == "exists epsilon_G>0; no explicit numerical epsilon is claimed"
        and any("sector" in item.lower() for item in scope["not_proved"])
        and any("multigraph" in item.lower() for item in scope["not_proved"]),
    )

    if FAILURES:
        print(f"FAIL: {len(FAILURES)} checks")
        return 1
    print("PASS: independent token-splitting obstruction verifier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
