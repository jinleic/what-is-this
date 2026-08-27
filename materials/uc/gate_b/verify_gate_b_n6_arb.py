#!/usr/bin/env python3
"""Arb certificate for the independent six-coordinate Gate B base.

The generic interval Bellman evaluator comes from the standalone n=7 verifier,
but this script independently reconstructs and checks the 25-row 3+3 family.
It provides an arithmetic cross-check of the all-720-order dyadic relaxation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys
from typing import Any

from flint import arb

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import verify_gate_b as verifier


DIMENSION = 6
LEFT_COORDINATES = (0, 1, 2)
RIGHT_COORDINATES = (3, 4, 5)
BLOCK_CELLS = frozenset(((0, 0), (0, 1), (1, 0), (1, 2), (2, 1)))
ALPHA = Fraction(356_069, 10_000_000)
FAMILY = (
    0, 1, 2, 4, 8, 11, 13, 14, 16, 19, 21, 22, 25,
    26, 28, 32, 35, 37, 38, 41, 42, 44, 49, 50, 52,
)


def _subsets_of_weight(
    coordinates: tuple[int, ...],
    weight: int,
) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in chosen)
        for chosen in combinations(coordinates, weight)
    )


def reconstruct_family() -> tuple[int, ...]:
    return tuple(sorted({
        left | right
        for left_weight, right_weight in BLOCK_CELLS
        for left in _subsets_of_weight(LEFT_COORDINATES, left_weight)
        for right in _subsets_of_weight(RIGHT_COORDINATES, right_weight)
    }))


def _exact_base() -> dict[str, Any]:
    family = reconstruct_family()
    assert family == FAMILY
    assert len(family) == len(set(family)) == 25
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    assert counts == (10,) * DIMENSION
    incidence = sum(counts)
    assert incidence == 60
    threshold = verifier.exact_reimer_threshold(len(family))
    assert threshold == 59
    assert incidence >= threshold
    assert max(counts) == 2 * len(family) // 5

    support = set(family)
    missing = sum(
        (left | right) not in support
        for left in family
        for right in family
    )
    defect = Fraction(missing, len(family) ** 2)
    assert missing == 444
    assert defect == Fraction(444, 625)

    columns = tuple(
        tuple((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    assert all(any(column) for column in columns)
    assert len(set(columns)) == DIMENSION
    assert 1 in support and 2 in support and 3 not in support
    return {
        "dimension": DIMENSION,
        "family_size": len(family),
        "block_cells": [list(cell) for cell in sorted(BLOCK_CELLS)],
        "family_rows": list(family),
        "coordinate_counts": list(counts),
        "total_incidence": incidence,
        "reimer_threshold": threshold,
        "missing_ordered_join_pairs": missing,
        "closure_defect": str(defect),
        "join_success_probability": str(1 - defect),
        "active": True,
        "separating": True,
        "ordered_pairs_with_replacement": True,
    }


def verify() -> dict[str, Any]:
    exact = _exact_base()
    family = FAMILY
    automorphisms, orders = verifier.automorphisms_and_order_representatives(
        family,
        DIMENSION,
    )
    assert len(automorphisms) == 72
    assert len(orders) == 10

    iid = verifier.certified_shapley_iid(family, DIMENSION)
    roots, state_count, ambiguous = verifier.certified_one_sided_costs(
        family,
        DIMENSION,
        orders,
    )
    one_sided = sum(roots, verifier.ZERO) / len(roots)
    entropy = arb(len(family)).log() / verifier.LOG2
    alpha = verifier.aq(ALPHA)
    a_plus = (
        (verifier.ONE - alpha) * iid
        + alpha * one_sided
        - entropy
    )
    negative_delta = Fraction(17, 1250)
    assert a_plus.upper() < verifier.aq(-negative_delta)
    assert state_count == 1434
    assert ambiguous == 305

    defect = Fraction(444, 625)
    ratio_lower = negative_delta / defect
    assert ratio_lower == Fraction(17, 888)
    assert 25**5 < 2**24

    return {
        "schema_version": 1,
        "verdict": "PROVED_GATE_B_UNBOUNDED_SECOND_BASE_ARB",
        "arithmetic": "python_flint_arb",
        "precision_bits": verifier.PRECISION_BITS,
        "alpha": {
            "numerator": ALPHA.numerator,
            "denominator": ALPHA.denominator,
            "dimension_independent": True,
        },
        "exact_base": exact,
        "symmetry": {
            "automorphism_count": len(automorphisms),
            "equal_order_orbit_count": len(orders),
            "order_orbit_size": len(automorphisms),
            "all_order_count": 720,
            "partition_verified_exactly": True,
        },
        "arb": {
            "shapley_iid": verifier.arb_bounds(iid),
            "one_sided_cost": verifier.arb_bounds(one_sided),
            "family_entropy": verifier.arb_bounds(entropy),
            "a_plus": verifier.arb_bounds(a_plus),
            "bellman_state_count": state_count,
            "ambiguous_clamp_state_count": ambiguous,
        },
        "rational_consequence": {
            "a_plus_upper_lt": "-17/1250",
            "base_ratio_lower_gt": str(ratio_lower),
        },
        "infinite_construction": {
            "family": "D^(boxtimes k)",
            "dimension": "6*k",
            "family_size": "25^k",
            "coordinate_count": "10*25^(k-1)",
            "closure_defect": "1-(181/625)^k",
            "a_plus_upper": "-17*k/1250",
            "repair_ratio_lower": (
                "(17*k/1250)/(1-(181/625)^k) > 17*k/1250"
            ),
            "reimer_integer_witness": "25^5 < 2^24",
            "zero_defect_convention_used": False,
        },
        "proof_boundary": (
            "The finite negative base is certified here with Arb. Product "
            "additivity and all-k admissibility are exact lemmas in PROOF.md."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", type=Path)
    args = parser.parse_args()
    certificate = verify()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write_certificate is not None:
        args.write_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.write_certificate.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
