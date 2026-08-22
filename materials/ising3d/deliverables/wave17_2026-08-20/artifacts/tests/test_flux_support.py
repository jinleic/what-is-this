#!/usr/bin/env python3
"""Standalone exact checks for the support-resolved Ising commutant certificate.

The mathematical sample is rebuilt independently from the experiment module.
Stored exact witnesses are decoded and substituted directly in global ordered-
Pauli commutator equations.
"""
from __future__ import annotations

import json
import signal
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "integrability" / "flux_support.json"
TIMEOUT_SECONDS = 420
PRIMES = (2_147_483_647, 2_147_483_629)


def grid_edges(rows: int, columns: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        [(r * columns + c, (r + 1) * columns + c) for r in range(rows - 1) for c in range(columns)]
        + [(r * columns + c, r * columns + c + 1) for r in range(rows) for c in range(columns - 1)]
    )


def chain_edges(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((site, site + 1) for site in range(n - 1))


def columns(n: int, edges: tuple[tuple[int, int], ...], support_mask: int) -> list[dict[tuple[int, int], int]]:
    sites = [site for site in range(n) if (support_mask >> site) & 1]
    incident = [edge for edge in edges if ((support_mask >> edge[0]) & 1) or ((support_mask >> edge[1]) & 1)]
    output = []
    for column in range(4 ** len(sites)):
        digits = column
        a_mask = 0
        b_mask = 0
        for site in sites:
            label = digits & 3
            digits >>= 2
            if label & 1:
                a_mask |= 1 << site
            if label & 2:
                b_mask |= 1 << site
        pauli = a_mask | (b_mask << n)
        vector: dict[tuple[int, int], int] = {}
        for site in sites:
            if (b_mask >> site) & 1:
                vector[(0, pauli ^ (1 << site))] = 1
        for left, right in incident:
            if ((a_mask >> left) ^ ((a_mask >> right) & 1)) & 1:
                vector[(1, pauli ^ (1 << (n + left)) ^ (1 << (n + right)))] = -1
        output.append(vector)
    return output


def rank_mod(vectors: list[dict[tuple[int, int], int]], prime: int) -> int:
    frequencies = Counter(row for vector in vectors for row in vector)
    order = {row: index for index, row in enumerate(sorted(frequencies, key=lambda row: (frequencies[row], row)))}
    pivots: dict[tuple[int, int], dict[tuple[int, int], int]] = {}
    rank = 0
    for source in sorted(vectors, key=len):
        vector = {row: value % prime for row, value in source.items() if value % prime}
        while vector:
            lead = min(vector, key=order.__getitem__)
            old = pivots.get(lead)
            if old is None:
                inverse = pow(vector[lead], prime - 2, prime)
                pivots[lead] = {row: value * inverse % prime for row, value in vector.items()}
                rank += 1
                break
            factor = vector[lead]
            for row, value in old.items():
                reduced = (vector.get(row, 0) - factor * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return rank


def exact_commutator(
    terms: list[list[int]], n: int, edges: tuple[tuple[int, int], ...]
) -> tuple[dict[int, Fraction], dict[int, Fraction]]:
    ad_a: defaultdict[int, Fraction] = defaultdict(Fraction)
    ad_b: defaultdict[int, Fraction] = defaultdict(Fraction)
    for a_mask, b_mask, numerator, denominator in terms:
        coefficient = Fraction(numerator, denominator)
        pauli = a_mask | (b_mask << n)
        for site in range(n):
            if (b_mask >> site) & 1:
                ad_a[pauli ^ (1 << site)] += coefficient
        for left, right in edges:
            if ((a_mask >> left) ^ ((a_mask >> right) & 1)) & 1:
                ad_b[pauli ^ (1 << (n + left)) ^ (1 << (n + right))] -= coefficient
    return (
        {pauli: coefficient for pauli, coefficient in ad_a.items() if coefficient},
        {pauli: coefficient for pauli, coefficient in ad_b.items() if coefficient},
    )


def row_for_mask(layer: dict, support_mask: int) -> dict:
    return next(row for row in layer["table"] if row["support_mask"] == support_mask)


def main() -> None:
    signal.alarm(TIMEOUT_SECONDS)
    artifact = json.loads(RESULT.read_text())
    assert set(artifact) == {"provenance", "data", "checks"}
    provenance = artifact["provenance"]
    assert provenance["script"] == "experiments/e57_flux_support.py"
    assert provenance["interpreter"] == ".venv/bin/python"
    assert provenance["generated_utc"] and provenance["method"]
    assert artifact["checks"] and all(check["passed"] for check in artifact["checks"])

    layers = {layer["graph"]: layer for layer in artifact["data"]["layers"]}
    assert {name: layer["scalar_only_for_every_connected_support_through"] for name, layer in layers.items()} == {
        "2x4_open_grid": 7,
        "3x3_open_grid": 6,
        "3x4_open_grid": 7,
    }

    samples = [
        ("2x4_open_grid", 2, 4, 63),
        ("2x4_open_grid", 2, 4, 127),
        ("3x3_open_grid", 3, 3, 63),
        ("3x4_open_grid", 3, 4, 127),
    ]
    for graph, rows, columns_count, support_mask in samples:
        stored = row_for_mask(layers[graph], support_mask)
        matrix = columns(rows * columns_count, grid_edges(rows, columns_count), support_mask)
        nullities = [len(matrix) - rank_mod(matrix, prime) for prime in PRIMES]
        assert nullities == [1, 1]
        assert [entry["nullity"] for entry in stored["modular"]] == nullities
        assert stored["kernel_dimension_Q"] == 1

    anomaly_layer = layers["3x3_open_grid"]
    assert sorted(anomaly_layer["exceptional_support_masks"]) == [239, 254, 367, 381]
    for support_mask, expected_dimension in ((239, 4), (367, 2)):
        stored = row_for_mask(anomaly_layer, support_mask)
        assert stored["kernel_dimension_Q"] == expected_dimension
        assert stored["basis_exact_rank_Q"] == expected_dimension
        assert stored["exact_substitution_passed"]
        assert stored["parent_review_flag"].startswith("LOCAL CONSERVED CHARGE")
        for vector in stored["exact_basis"]:
            assert exact_commutator(vector["terms_a_b_num_den"], 9, grid_edges(3, 3)) == ({}, {})
        matrix = columns(9, grid_edges(3, 3), support_mask)
        assert [len(matrix) - rank_mod(matrix, prime) for prime in PRIMES] == [expected_dimension] * 2

    chain = next(control for control in artifact["data"]["controls"] if control["name"] == "open_chain_P4_full_support")
    chain_matrix = columns(4, chain_edges(4), 0b1111)
    assert [len(chain_matrix) - rank_mod(chain_matrix, prime) for prime in PRIMES] == [5, 5]
    assert chain["kernel_dimension_Q"] == 5
    for vector in chain["exact_basis"]:
        assert exact_commutator(vector["terms_a_b_num_den"], 4, chain_edges(4)) == ({}, {})

    c4 = {control["name"]: control for control in artifact["data"]["controls"]}
    assert c4["cycle_C4_three_site_support"]["kernel_dimension_Q"] == 3
    assert c4["cycle_C4_full_support"]["kernel_dimension_Q"] == 27
    signal.alarm(0)
    print("OK")


if __name__ == "__main__":
    main()
