#!/usr/bin/env python3
"""Third-party-free certificate for a second Gate B Cartesian-power base.

This verifier reconstructs the 25-row, six-coordinate family from its 3+3
weight cells and evaluates the clamp-free Bellman relaxation for all 720
coordinate orders.  It deliberately uses no symmetry quotient.  The outward-
rounded dyadic interval primitive and transcendental bounds come from the
independently audited ``verify_gate_b_dyadic`` module; the family, all-order
evaluation, exact filters, defect, and all-power consequence are separate from
the seven-coordinate certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import factorial
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_gate_b_dyadic import (
    ALPHA,
    EXP_TERMS,
    LN2,
    LOG_TERMS,
    PRECISION_BITS,
    Interval,
    certified_shapley_iid,
    exact_reimer_threshold,
    natural_log,
    relaxed_bellman_cost,
)


DIMENSION = 6
LEFT_COORDINATES = (0, 1, 2)
RIGHT_COORDINATES = (3, 4, 5)
BLOCK_CELLS = frozenset(((0, 0), (0, 1), (1, 0), (1, 2), (2, 1)))


def _subsets_of_weight(
    coordinates: tuple[int, ...],
    weight: int,
) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in chosen)
        for chosen in combinations(coordinates, weight)
    )


def reconstruct_family() -> tuple[int, ...]:
    rows = {
        left | right
        for left_weight, right_weight in BLOCK_CELLS
        for left in _subsets_of_weight(LEFT_COORDINATES, left_weight)
        for right in _subsets_of_weight(RIGHT_COORDINATES, right_weight)
    }
    return tuple(sorted(rows))



def _ordered_weight_union_count(
    block_size: int,
    left_weight: int,
    right_weight: int,
    union_weight: int,
) -> int:
    intersection = left_weight + right_weight - union_weight
    left_only = left_weight - intersection
    right_only = right_weight - intersection
    outside = block_size - union_weight
    if min(intersection, left_only, right_only, outside) < 0:
        return 0
    return factorial(block_size) // (
        factorial(intersection)
        * factorial(left_only)
        * factorial(right_only)
        * factorial(outside)
    )


def _cell_join_success_count() -> int:
    return sum(
        _ordered_weight_union_count(3, left[0], right[0], joined[0])
        * _ordered_weight_union_count(3, left[1], right[1], joined[1])
        for left in BLOCK_CELLS
        for right in BLOCK_CELLS
        for joined in BLOCK_CELLS
    )

def _exact_base(family: tuple[int, ...]) -> dict[str, Any]:
    expected = (
        0, 1, 2, 4, 8, 11, 13, 14, 16, 19, 21, 22, 25,
        26, 28, 32, 35, 37, 38, 41, 42, 44, 49, 50, 52,
    )
    assert family == expected
    assert len(family) == len(set(family)) == 25

    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    assert counts == (10,) * DIMENSION
    incidence = sum(counts)
    assert incidence == 60
    threshold = exact_reimer_threshold(len(family))
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
    cell_success = _cell_join_success_count()
    assert cell_success == len(family) ** 2 - missing == 181

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
        "block_partition": [list(LEFT_COORDINATES), list(RIGHT_COORDINATES)],
        "block_cells": [list(cell) for cell in sorted(BLOCK_CELLS)],
        "family_rows": list(family),
        "coordinate_counts": list(counts),
        "total_incidence": incidence,
        "reimer_threshold": threshold,
        "missing_ordered_join_pairs": missing,
        "cell_formula_join_success_pairs": cell_success,
        "closure_defect": str(defect),
        "join_success_probability": str(1 - defect),
        "active": True,
        "separating": True,
        "explicit_missing_join": [1, 2, 3],
        "ordered_pairs_with_replacement": True,
    }


def verify() -> dict[str, Any]:
    family = reconstruct_family()
    exact_base = _exact_base(family)

    iid = certified_shapley_iid(family, DIMENSION)
    roots: list[Interval] = []
    state_count = 0
    for order in permutations(range(DIMENSION)):
        root, states = relaxed_bellman_cost(family, order)
        roots.append(root)
        state_count += states
    assert len(roots) == factorial(DIMENSION) == 720

    relaxed = Interval.zero()
    for root in roots:
        relaxed = relaxed + root
    relaxed = relaxed.scale(Fraction(1, len(roots)))

    order_classes = Counter((root.lo, root.hi) for root in roots)
    assert len(order_classes) == 10
    assert set(order_classes.values()) == {72}

    family_entropy = natural_log(Interval.from_int(len(family))) / LN2
    a_plus_upper = (
        iid.scale(1 - ALPHA)
        + relaxed.scale(ALPHA)
        - family_entropy
    )
    negative_delta = Fraction(1, 80)
    assert a_plus_upper.upper_lt(-negative_delta)

    defect = Fraction(444, 625)
    success = 1 - defect
    ratio_lower = negative_delta / defect
    assert success == Fraction(181, 625)
    assert ratio_lower == Fraction(125, 7104)
    assert 25**5 < 2**24

    smallest_root = min(roots, key=lambda interval: interval.lo)
    largest_root = max(roots, key=lambda interval: interval.hi)
    return {
        "schema_version": 1,
        "verdict": "PROVED_GATE_B_UNBOUNDED_SECOND_BASE",
        "arithmetic": "pure_python_dyadic_intervals",
        "precision_bits": PRECISION_BITS,
        "transcendental_bounds": {
            "natural_log": (
                f"{LOG_TERMS}-term atanh series with exact rational tail"
            ),
            "exponential": (
                f"{EXP_TERMS}-term Taylor series after x/2^j <= 1/8, "
                "with exact geometric tail and repeated squaring"
            ),
            "third_party_dependencies": [],
            "arithmetic_module": "verify_gate_b_dyadic.py",
        },
        "alpha": {
            "numerator": ALPHA.numerator,
            "denominator": ALPHA.denominator,
            "dimension_independent": True,
        },
        "exact_base": exact_base,
        "all_order_evaluation": {
            "symmetry_quotient_used": False,
            "evaluated_order_count": len(roots),
            "distinct_root_interval_count": len(order_classes),
            "multiplicity_per_root_interval": 72,
            "total_bellman_state_count": state_count,
            "minimum_root_interval": smallest_root.to_json(),
            "maximum_root_interval": largest_root.to_json(),
        },
        "intervals": {
            "shapley_iid": iid.to_json(),
            "clamp_free_bellman_relaxation": relaxed.to_json(),
            "family_entropy": family_entropy.to_json(),
            "a_plus_upper_relaxation": a_plus_upper.to_json(),
        },
        "relaxation": {
            "clamp_free_upper_bound": True,
            "statement": (
                "At every nondegenerate Bellman state, maximize the affine "
                "entropy objective over [0,1], a superset of the feasible "
                "one-sided interval. Induction is monotone on feasible "
                "transitions because their four probabilities are nonnegative."
            ),
            "does_not_use_order_symmetry": True,
            "does_not_use_arb_clamp_classification": True,
        },
        "rational_consequence": {
            "a_plus_upper_lt": "-1/80",
            "base_ratio_lower_gt": str(ratio_lower),
        },
        "infinite_construction": {
            "family": "D^(boxtimes k)",
            "dimension": "6*k",
            "family_size": "25^k",
            "coordinate_count": "10*25^(k-1)",
            "closure_defect": "1-(181/625)^k",
            "a_plus_upper": "-k/80",
            "repair_ratio_lower": "(k/80)/(1-(181/625)^k) > k/80",
            "reimer_integer_witness": "25^5 < 2^24",
            "zero_defect_convention_used": False,
        },
        "proof_boundary": (
            "The second finite negative base and all 720 order values are "
            "certified here. Product additivity and all-k admissibility are "
            "exact lemmas proved in PROOF.md, not finite-n extrapolation."
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
