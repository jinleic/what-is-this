#!/usr/bin/env python3
"""Standalone verifier for the Gate B Cartesian-power counterexample.

Dependencies: Python >= 3.9 and python-flint.  This file deliberately imports
no module from ``math/uc``.  It independently checks the exact base-family
constraints, recomputes the Shapley iid term and one-sided Bellman upper
certificate with 256-bit Arb arithmetic, and directly checks the square of the
base family.  The accompanying proof establishes the product identities for
all powers.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path
from typing import Any

from flint import arb, ctx


PRECISION_BITS = 256
ctx.prec = PRECISION_BITS
ZERO = arb(0)
ONE = arb(1)
LOG2 = arb(2).log()


def aq(value: Fraction) -> arb:
    """Embed a rational exactly into an Arb ball."""
    return arb(value.numerator) / value.denominator


def binary_entropy(value: Fraction) -> arb:
    if value == 0 or value == 1:
        return ZERO
    point = aq(value)
    return -(point * point.log() + (ONE - point) * (ONE - point).log()) / LOG2


def exact_reimer_threshold(size: int) -> int:
    """Return ceil(size*log2(size)/2) using integer comparisons only."""
    if size <= 1:
        return 0
    target = size**size
    low = 0
    high = size * (size - 1).bit_length()
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def fraction_from_json(value: dict[str, int]) -> Fraction:
    return Fraction(value["numerator"], value["denominator"])


def family_bitset(family: tuple[int, ...]) -> int:
    return sum(1 << row for row in family)

def popcount(value: int) -> int:
    """Count set bits on Python 3.9 and newer."""
    if hasattr(int, "bit_count"):
        return value.bit_count()
    return bin(value).count("1")


def permute_subset(subset: int, permutation: tuple[int, ...]) -> int:
    image = 0
    for old_coordinate, new_coordinate in enumerate(permutation):
        if subset >> old_coordinate & 1:
            image |= 1 << new_coordinate
    return image


def predecessor_law(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Conditional-one probabilities and prefix masses, independently built."""
    groups: dict[int, list[int]] = {}
    for row in family:
        prefix = sum(
            ((row >> predecessor) & 1) << index
            for index, predecessor in enumerate(predecessors)
        )
        total, ones = groups.get(prefix, [0, 0])
        groups[prefix] = [total + 1, ones + ((row >> coordinate) & 1)]
    size = len(family)
    return tuple(
        (Fraction(ones, total), Fraction(total, size))
        for total, ones in groups.values()
    )


def certified_iid_local(law: tuple[tuple[Fraction, Fraction], ...]) -> arb:
    return sum(
        (
            aq(left_mass * right_mass)
            * binary_entropy(left + right - left * right)
            for left, left_mass in law
            for right, right_mass in law
        ),
        ZERO,
    )


def certified_shapley_iid(family: tuple[int, ...], dimension: int) -> arb:
    total = ZERO
    for coordinate in range(dimension):
        others = tuple(i for i in range(dimension) if i != coordinate)
        for predecessor_count in range(dimension):
            weight = Fraction(
                1, dimension * comb(dimension - 1, predecessor_count)
            )
            for predecessors in combinations(others, predecessor_count):
                total += aq(weight) * certified_iid_local(
                    predecessor_law(family, coordinate, predecessors)
                )
    return total


def sstar(left: Fraction, right: Fraction) -> Fraction:
    return max(left, right, min(left + right, Fraction(1, 2)))


def automorphisms_and_order_representatives(
    family: tuple[int, ...], dimension: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    orders = tuple(permutations(range(dimension)))
    mask = family_bitset(family)
    automorphisms = tuple(
        permutation
        for permutation in orders
        if family_bitset(
            tuple(permute_subset(row, permutation) for row in family)
        )
        == mask
    )

    unseen = set(orders)
    representatives: list[tuple[int, ...]] = []
    while unseen:
        order = min(unseen)
        orbit = {
            tuple(automorphism[coordinate] for coordinate in order)
            for automorphism in automorphisms
        }
        assert len(orbit) == len(automorphisms)
        assert orbit <= unseen
        representatives.append(order)
        unseen.difference_update(orbit)
    assert len(representatives) * len(automorphisms) == len(orders)
    return automorphisms, tuple(representatives)


def certified_one_sided_costs(
    family: tuple[int, ...],
    dimension: int,
    orders: tuple[tuple[int, ...], ...],
) -> tuple[tuple[arb, ...], int, int]:
    """Return rigorous upper balls for all supplied-order Bellman maxima."""
    row_ones = tuple(
        sum(
            1 << row
            for row in range(1 << dimension)
            if row >> coordinate & 1
        )
        for coordinate in range(dimension)
    )
    root = family_bitset(family)
    cache: dict[tuple[tuple[int, ...], int, int], arb] = {}
    ambiguous = 0

    def value(suffix: tuple[int, ...], left_rows: int, right_rows: int) -> arb:
        nonlocal ambiguous
        if not suffix:
            return ZERO
        if left_rows > right_rows:
            left_rows, right_rows = right_rows, left_rows
        key = (suffix, left_rows, right_rows)
        if key in cache:
            return cache[key]

        coordinate = suffix[0]
        rest = suffix[1:]
        ones = row_ones[coordinate]
        left_one = left_rows & ones
        right_one = right_rows & ones
        left_zero = left_rows ^ left_one
        right_zero = right_rows ^ right_one
        left = Fraction(popcount(left_one), popcount(left_rows))
        right = Fraction(popcount(right_one), popcount(right_rows))
        lower = sstar(left, right)
        upper = min(Fraction(1), left + right)

        children: dict[tuple[int, int], arb] = {}
        for left_bit, left_child in ((0, left_zero), (1, left_one)):
            if not left_child:
                continue
            for right_bit, right_child in ((0, right_zero), (1, right_one)):
                if not right_child:
                    continue
                children[left_bit, right_bit] = value(
                    rest, left_child, right_child
                )

        def boundary(action: Fraction) -> arb:
            probabilities = (
                (0, 0, 1 - action),
                (1, 0, action - right),
                (0, 1, action - left),
                (1, 1, left + right - action),
            )
            return binary_entropy(action) + sum(
                (
                    aq(probability) * children.get((left_bit, right_bit), ZERO)
                    for left_bit, right_bit, probability in probabilities
                ),
                ZERO,
            )

        if lower == upper:
            result = boundary(lower)
        else:
            assert len(children) == 4
            slope = (
                -children[0, 0]
                + children[1, 0]
                + children[0, 1]
                - children[1, 1]
            )
            constant = (
                children[0, 0]
                - aq(right) * children[1, 0]
                - aq(left) * children[0, 1]
                + aq(left + right) * children[1, 1]
            )
            power = (slope * LOG2).exp()
            logistic = power / (ONE + power)
            unconstrained = constant + (ONE + power).log() / LOG2
            lower_ball = aq(lower)
            upper_ball = aq(upper)
            if logistic.upper() < lower_ball:
                result = boundary(lower)
            elif logistic.lower() > upper_ball:
                result = boundary(upper)
            elif logistic.lower() > lower_ball and logistic.upper() < upper_ball:
                result = unconstrained
            else:
                ambiguous += 1
                candidates = [unconstrained]
                if logistic.lower() <= lower_ball:
                    candidates.append(boundary(lower))
                if logistic.upper() >= upper_ball:
                    candidates.append(boundary(upper))
                result = candidates[0]
                for candidate in candidates[1:]:
                    result = result.union(candidate)
        cache[key] = result
        return result

    roots = tuple(value(tuple(order), root, root) for order in orders)
    return roots, len(cache), ambiguous

def certified_square_order_cost(
    family: tuple[int, ...],
    dimension: int,
    global_order: tuple[int, ...],
) -> tuple[arb, int, int]:
    """Evaluate the actual product-family Bellman DP for one global order.

    A product fiber is stored as two base-family fiber masks.  This is an exact
    representation of the 2,025-row square, not an invocation of the
    tensorization identity.
    """
    assert len(global_order) == 2 * dimension
    assert set(global_order) == set(range(2 * dimension))
    row_ones = tuple(
        sum(
            1 << row
            for row in range(1 << dimension)
            if row >> coordinate & 1
        )
        for coordinate in range(dimension)
    )
    root = family_bitset(family)
    cache: dict[
        tuple[tuple[int, ...], tuple[int, int], tuple[int, int]], arb
    ] = {}
    ambiguous = 0

    def value(
        suffix: tuple[int, ...],
        left_fibers: tuple[int, int],
        right_fibers: tuple[int, int],
    ) -> arb:
        nonlocal ambiguous
        if not suffix:
            return ZERO

        active_blocks = {coordinate // dimension for coordinate in suffix}
        left_fibers = tuple(
            fiber if block in active_blocks else 0
            for block, fiber in enumerate(left_fibers)
        )
        right_fibers = tuple(
            fiber if block in active_blocks else 0
            for block, fiber in enumerate(right_fibers)
        )
        if left_fibers > right_fibers:
            left_fibers, right_fibers = right_fibers, left_fibers
        key = (suffix, left_fibers, right_fibers)
        if key in cache:
            return cache[key]

        global_coordinate = suffix[0]
        rest = suffix[1:]
        block, coordinate = divmod(global_coordinate, dimension)
        left_rows = left_fibers[block]
        right_rows = right_fibers[block]
        ones = row_ones[coordinate]
        left_one = left_rows & ones
        right_one = right_rows & ones
        left_zero = left_rows ^ left_one
        right_zero = right_rows ^ right_one
        left = Fraction(popcount(left_one), popcount(left_rows))
        right = Fraction(popcount(right_one), popcount(right_rows))
        lower = sstar(left, right)
        upper = min(Fraction(1), left + right)

        children: dict[tuple[int, int], arb] = {}
        for left_bit, left_child in ((0, left_zero), (1, left_one)):
            if not left_child:
                continue
            for right_bit, right_child in ((0, right_zero), (1, right_one)):
                if not right_child:
                    continue
                next_left = list(left_fibers)
                next_right = list(right_fibers)
                next_left[block] = left_child
                next_right[block] = right_child
                children[left_bit, right_bit] = value(
                    rest, tuple(next_left), tuple(next_right)
                )

        def boundary(action: Fraction) -> arb:
            probabilities = (
                (0, 0, 1 - action),
                (1, 0, action - right),
                (0, 1, action - left),
                (1, 1, left + right - action),
            )
            return binary_entropy(action) + sum(
                (
                    aq(probability) * children.get((left_bit, right_bit), ZERO)
                    for left_bit, right_bit, probability in probabilities
                ),
                ZERO,
            )

        if lower == upper:
            result = boundary(lower)
        else:
            assert len(children) == 4
            slope = (
                -children[0, 0]
                + children[1, 0]
                + children[0, 1]
                - children[1, 1]
            )
            constant = (
                children[0, 0]
                - aq(right) * children[1, 0]
                - aq(left) * children[0, 1]
                + aq(left + right) * children[1, 1]
            )
            power = (slope * LOG2).exp()
            logistic = power / (ONE + power)
            unconstrained = constant + (ONE + power).log() / LOG2
            lower_ball = aq(lower)
            upper_ball = aq(upper)
            if logistic.upper() < lower_ball:
                result = boundary(lower)
            elif logistic.lower() > upper_ball:
                result = boundary(upper)
            elif logistic.lower() > lower_ball and logistic.upper() < upper_ball:
                result = unconstrained
            else:
                ambiguous += 1
                candidates = [unconstrained]
                if logistic.lower() <= lower_ball:
                    candidates.append(boundary(lower))
                if logistic.upper() >= upper_ball:
                    candidates.append(boundary(upper))
                result = candidates[0]
                for candidate in candidates[1:]:
                    result = result.union(candidate)
        cache[key] = result
        return result

    result = value(tuple(global_order), (root, root), (root, root))
    return result, len(cache), ambiguous


def exact_family_checks(data: dict[str, Any]) -> dict[str, Any]:
    dimension = int(data["dimension"])
    family = tuple(int(row) for row in data["family_rows"])
    expected = data["expected"]
    assert family == tuple(sorted(set(family)))
    assert all(0 <= row < 1 << dimension for row in family)
    assert len(family) == expected["family_size"]

    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(dimension)
    )
    incidence = sum(counts)
    assert list(counts) == expected["coordinate_counts"]
    assert incidence == expected["total_incidence"]
    assert max(counts) <= 2 * len(family) // 5
    threshold = exact_reimer_threshold(len(family))
    assert incidence >= threshold
    assert 1 << (2 * incidence) >= len(family) ** len(family)

    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    defect = Fraction(missing, len(family) ** 2)
    assert missing == expected["missing_ordered_join_pairs"]
    assert defect == fraction_from_json(expected["closure_defect"])

    columns = tuple(
        tuple((row >> coordinate) & 1 for row in family)
        for coordinate in range(dimension)
    )
    active = all(any(column) for column in columns)
    separating = len(set(columns)) == dimension
    assert active is expected["active"]
    assert separating is expected["separating"]

    left_coordinates, right_coordinates = (
        tuple(coordinates) for coordinates in data["block_partition"]
    )
    assert sorted(left_coordinates + right_coordinates) == list(range(dimension))
    assert set(left_coordinates).isdisjoint(right_coordinates)
    declared_cells = {tuple(cell) for cell in data["block_cells"]}
    cell_union = tuple(
        row
        for row in range(1 << dimension)
        if (
            sum((row >> coordinate) & 1 for coordinate in left_coordinates),
            sum((row >> coordinate) & 1 for coordinate in right_coordinates),
        )
        in declared_cells
    )
    assert family == cell_union
    observed_cells = {
        (
            sum((row >> coordinate) & 1 for coordinate in left_coordinates),
            sum((row >> coordinate) & 1 for coordinate in right_coordinates),
        )
        for row in family
    }
    assert observed_cells == declared_cells
    left, right, joined = expected["explicit_missing_join"]
    assert left in support and right in support and joined == left | right
    assert joined not in support

    return {
        "family": family,
        "dimension": dimension,
        "counts": counts,
        "incidence": incidence,
        "reimer_threshold": threshold,
        "missing": missing,
        "defect": defect,
        "active": active,
        "separating": separating,
        "cell_union_verified": True,
    }


def cartesian_power_rows(
    family: tuple[int, ...], dimension: int, power: int
) -> tuple[int, ...]:
    rows = (0,)
    for block in range(power):
        rows = tuple(
            prefix | (row << (block * dimension))
            for prefix in rows
            for row in family
        )
    return tuple(sorted(rows))


def direct_square_checks(base: dict[str, Any]) -> dict[str, Any]:
    family = base["family"]
    dimension = base["dimension"]
    square = cartesian_power_rows(family, dimension, 2)
    assert len(square) == len(family) ** 2
    counts = tuple(
        sum((row >> coordinate) & 1 for row in square)
        for coordinate in range(2 * dimension)
    )
    expected_count = base["counts"][0] * len(family)
    assert counts == (expected_count,) * (2 * dimension)
    incidence = sum(counts)
    assert incidence == 2 * base["incidence"] * len(family)
    assert incidence >= exact_reimer_threshold(len(square))

    support = set(square)
    missing = sum(
        (left | right) not in support for left in square for right in square
    )
    expected_defect = 1 - (1 - base["defect"]) ** 2
    assert Fraction(missing, len(square) ** 2) == expected_defect

    columns = tuple(
        tuple((row >> coordinate) & 1 for row in square)
        for coordinate in range(2 * dimension)
    )
    assert all(any(column) for column in columns)
    assert len(set(columns)) == 2 * dimension
    return {
        "family_size": len(square),
        "dimension": 2 * dimension,
        "coordinate_count": expected_count,
        "total_incidence": incidence,
        "missing_ordered_join_pairs": missing,
        "closure_defect": expected_defect,
    }


def arb_bounds(value: arb) -> dict[str, str]:
    return {
        "ball": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
    }


def verify(candidate_path: Path) -> dict[str, Any]:
    data = json.loads(candidate_path.read_text())
    base = exact_family_checks(data)
    family = base["family"]
    dimension = base["dimension"]
    alpha = fraction_from_json(data["alpha"])
    assert alpha == Fraction(356069, 10000000)

    automorphisms, orders = automorphisms_and_order_representatives(
        family, dimension
    )
    expected = data["expected"]
    assert len(automorphisms) == expected["automorphism_count"]
    assert len(orders) == expected["order_representative_count"]

    iid = certified_shapley_iid(family, dimension)
    roots, states, ambiguous = certified_one_sided_costs(
        family, dimension, orders
    )
    one_sided = sum(roots, ZERO) / len(roots)
    entropy = arb(len(family)).log() / LOG2
    a_plus = (ONE - aq(alpha)) * iid + aq(alpha) * one_sided - entropy
    corrected = a_plus + aq(base["defect"] * Fraction(1, 50))
    assert base["defect"] > 0
    repair_ratio_lower = -a_plus / aq(base["defect"])

    a_plus_threshold = fraction_from_json(
        data["certified_rational_consequences"]["a_plus_upper_lt"]
    )
    corrected_threshold = fraction_from_json(
        data["certified_rational_consequences"][
            "corrected_one_fiftieth_upper_lt"
        ]
    )
    ratio_threshold = fraction_from_json(
        data["certified_rational_consequences"]["repair_ratio_lower_gt"]
    )
    assert a_plus.upper() < aq(a_plus_threshold)
    assert corrected.upper() < aq(corrected_threshold)
    assert repair_ratio_lower.lower() > aq(ratio_threshold)

    square = direct_square_checks(base)
    success = 1 - base["defect"]
    assert success == Fraction(17, 81)
    assert 45**5 < 2**28
    delta = -a_plus_threshold
    assert delta == Fraction(7, 250)

    identity_order = tuple(range(dimension))
    base_identity = roots[orders.index(identity_order)]
    square_orders = {
        "consecutive": tuple(range(2 * dimension)),
        "alternating": tuple(
            coordinate
            for local_coordinate in range(dimension)
            for coordinate in (local_coordinate, dimension + local_coordinate)
        ),
    }
    square_order_checks = {}
    for name, square_order in square_orders.items():
        square_cost, square_states, square_ambiguous = (
            certified_square_order_cost(family, dimension, square_order)
        )
        residual = square_cost - 2 * base_identity
        assert residual.contains(0)
        square_order_checks[name] = {
            "global_order": list(square_order),
            "actual_square_cost": arb_bounds(square_cost),
            "twice_base_order_cost": arb_bounds(2 * base_identity),
            "residual": arb_bounds(residual),
            "bellman_state_count": square_states,
            "ambiguous_clamp_state_count": square_ambiguous,
        }

    square_iid = 2 * iid
    square_one_sided = 2 * one_sided
    square_entropy = 2 * entropy
    square_a_plus = (
        (ONE - aq(alpha)) * square_iid
        + aq(alpha) * square_one_sided
        - square_entropy
    )

    checked_power_arithmetic = []
    for power in range(1, 9):
        power_size = len(family) ** power
        coordinate_count = base["counts"][0] * len(family) ** (power - 1)
        total_incidence = (
            base["incidence"] * power * len(family) ** (power - 1)
        )
        defect = 1 - success**power
        assert 5 * coordinate_count == 2 * power_size
        assert 5 * total_incidence == 14 * power * power_size
        assert 45 ** (5 * power) < 2 ** (28 * power)
        assert 0 < defect < 1
        ratio_lower_bound = delta * power / defect
        assert ratio_lower_bound > delta * power
        checked_power_arithmetic.append(
            {
                "power": power,
                "dimension": dimension * power,
                "family_size": power_size,
                "coordinate_count": coordinate_count,
                "total_incidence": total_incidence,
                "closure_defect": str(defect),
                "cap_holds_with_equality": True,
                "reimer_holds_strictly": True,
                "certified_ratio_lower_bound": str(ratio_lower_bound),
            }
        )

    return {
        "schema_version": 2,
        "verdict": "PROVED_GATE_B_UNBOUNDED",
        "candidate": str(candidate_path),
        "precision_bits": PRECISION_BITS,
        "alpha": {
            "numerator": alpha.numerator,
            "denominator": alpha.denominator,
            "dimension_independent": True,
        },
        "exact_base": {
            "dimension": dimension,
            "family_size": len(family),
            "coordinate_counts": list(base["counts"]),
            "total_incidence": base["incidence"],
            "reimer_threshold": base["reimer_threshold"],
            "missing_ordered_join_pairs": base["missing"],
            "closure_defect": str(base["defect"]),
            "join_success_probability": str(success),
            "active": base["active"],
            "separating": base["separating"],
            "automorphism_count": len(automorphisms),
            "order_representative_count": len(orders),
            "cell_union_verified": base["cell_union_verified"],
            "ordered_pairs_with_replacement": True,
        },
        "arb": {
            "shapley_iid": arb_bounds(iid),
            "one_sided_cost": arb_bounds(one_sided),
            "a_plus": arb_bounds(a_plus),
            "a_plus_plus_one_fiftieth_defect": arb_bounds(corrected),
            "repair_ratio_lower_certificate": {
                "lower": str(repair_ratio_lower.lower()),
                "proves_gt": str(ratio_threshold),
                "derived_from": "-upper(A_plus)/closure_defect",
            },
            "bellman_state_count": states,
            "ambiguous_clamp_state_count": ambiguous,
        },
        "direct_square": {
            **{key: value for key, value in square.items() if key != "closure_defect"},
            "closure_defect": str(square["closure_defect"]),
        },
        "tensorization_checks": {
            "actual_square_fixed_orders": square_order_checks,
        },
        "derived_from_proved_identities": {
            "evidence_origin": (
                "Computed from the product identities proved in PROOF.md; "
                "not independently recomputed on all square orders."
            ),
            "square_shapley_iid": arb_bounds(square_iid),
            "square_one_sided_cost": arb_bounds(square_one_sided),
            "square_a_plus": arb_bounds(square_a_plus),
        },
        "infinite_construction": {
            "base_dimension": dimension,
            "base_family_size": len(family),
            "base_coordinate_count": base["counts"][0],
            "base_total_incidence": base["incidence"],
            "base_join_success": {
                "numerator": success.numerator,
                "denominator": success.denominator,
            },
            "reimer_integer_witness": {
                "left": 45**5,
                "right": 2**28,
                "proves": "log2(45) < 28/5",
            },
            "base_columns_active_separating": (
                base["active"] and base["separating"]
            ),
            "base_columns_nonconstant": all(
                0 < count < len(family) for count in base["counts"]
            ),
            "certified_negative_delta": {
                "numerator": delta.numerator,
                "denominator": delta.denominator,
            },
            "base_repair_ratio_lower_gt": {
                "numerator": ratio_threshold.numerator,
                "denominator": ratio_threshold.denominator,
            },
            "zero_defect_extended_value_convention_used": False,
            "checked_power_arithmetic": checked_power_arithmetic,
            "formulas_proved_in_PROOF.md": {
                "dimension": "7*k",
                "family_size": "45^k",
                "coordinate_count": "18*45^(k-1)",
                "total_incidence": "126*k*45^(k-1)",
                "closure_defect": "1-(17/81)^k",
                "a_plus": "k*A_plus(F_1)",
                "ratio_lower_bound": "7*k/250",
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidate",
        type=Path,
        default=Path(__file__).parent / "candidates" / "n7_block_extremizer.json",
    )
    parser.add_argument("--write-certificate", type=Path)
    args = parser.parse_args()

    certificate = verify(args.candidate.resolve())
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write_certificate is not None:
        args.write_certificate.parent.mkdir(parents=True, exist_ok=True)
        args.write_certificate.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
