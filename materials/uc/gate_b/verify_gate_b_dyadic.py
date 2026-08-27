#!/usr/bin/env python3
"""Independent pure-Python certificate for the Gate B obstruction.

This checker deliberately does not import python-flint or any discovery module.
It reconstructs the seven-coordinate family from its weight cells, computes the
Shapley iid term with outward-rounded dyadic intervals, and upper-bounds the
one-sided Bellman term by relaxing every nondegenerate action interval to
``[0, 1]``.  The weaker relaxation still proves a strictly negative base.

All transcendental enclosures come from rational atanh/Taylor remainder bounds.
The exact Cartesian-product argument remains the mathematical theorem in
``PROOF.md``; this file supplies an arithmetically independent finite base.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path
from typing import Any


PRECISION_BITS = 160
SCALE = 1 << PRECISION_BITS
LOG_TERMS = 64
EXP_TERMS = 32

DIMENSION = 7
LEFT_COORDINATES = (5, 6)
RIGHT_COORDINATES = (0, 1, 2, 3, 4)
BLOCK_CELLS = frozenset(
    ((0, 1), (0, 3), (1, 0), (1, 2), (1, 5), (2, 0), (2, 1))
)
ALPHA = Fraction(356_069, 10_000_000)


def _ceil_div(numerator: int, denominator: int) -> int:
    assert denominator > 0
    return -((-numerator) // denominator)


@dataclass(frozen=True)
class Interval:
    """Closed dyadic interval with the fixed denominator ``2^PRECISION_BITS``."""

    lo: int
    hi: int

    def __post_init__(self) -> None:
        assert self.lo <= self.hi

    @staticmethod
    def zero() -> "Interval":
        return Interval(0, 0)

    @staticmethod
    def one() -> "Interval":
        return Interval(SCALE, SCALE)

    @staticmethod
    def from_int(value: int) -> "Interval":
        return Interval(value * SCALE, value * SCALE)

    @staticmethod
    def from_fraction(value: Fraction) -> "Interval":
        scaled = value.numerator * SCALE
        denominator = value.denominator
        return Interval(
            scaled // denominator,
            _ceil_div(scaled, denominator),
        )

    def contains_fraction(self, value: Fraction) -> bool:
        scaled = value.numerator * SCALE
        denominator = value.denominator
        return self.lo * denominator <= scaled <= self.hi * denominator

    def contains_interval(self, other: "Interval") -> bool:
        return self.lo <= other.lo and other.hi <= self.hi

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other: "Interval") -> "Interval":
        return Interval(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(
            min(products) // SCALE,
            _ceil_div(max(products), SCALE),
        )

    def scale(self, value: Fraction) -> "Interval":
        numerator = value.numerator
        denominator = value.denominator
        if numerator >= 0:
            return Interval(
                (self.lo * numerator) // denominator,
                _ceil_div(self.hi * numerator, denominator),
            )
        return Interval(
            (self.hi * numerator) // denominator,
            _ceil_div(self.lo * numerator, denominator),
        )

    def reciprocal(self) -> "Interval":
        assert self.lo > 0
        return Interval(
            (SCALE * SCALE) // self.hi,
            _ceil_div(SCALE * SCALE, self.lo),
        )

    def __truediv__(self, other: "Interval") -> "Interval":
        return self * other.reciprocal()

    def hull(self, other: "Interval") -> "Interval":
        return Interval(min(self.lo, other.lo), max(self.hi, other.hi))

    def upper_lt(self, value: Fraction) -> bool:
        return self.hi * value.denominator < value.numerator * SCALE

    def lower_gt(self, value: Fraction) -> bool:
        return self.lo * value.denominator > value.numerator * SCALE

    def to_json(self) -> dict[str, Any]:
        return {
            "denominator_power_of_two": PRECISION_BITS,
            "lower_numerator": self.lo,
            "upper_numerator": self.hi,
            "lower_decimal": _dyadic_decimal(self.lo),
            "upper_decimal": _dyadic_decimal(self.hi),
        }


def _dyadic_decimal(numerator: int, digits: int = 55) -> str:
    sign = "-" if numerator < 0 else ""
    numerator = abs(numerator)
    integer, remainder = divmod(numerator, SCALE)
    decimals = []
    for _ in range(digits):
        remainder *= 10
        digit, remainder = divmod(remainder, SCALE)
        decimals.append(str(digit))
    return f"{sign}{integer}." + "".join(decimals)


def _unit_log(value: Fraction) -> Interval:
    """Enclose ``ln(value)`` for ``1 <= value <= 2`` by an atanh series."""
    assert 1 <= value <= 2
    z = (value - 1) / (value + 1)
    if z == 0:
        return Interval.zero()

    z_interval = Interval.from_fraction(z)
    z_squared = z_interval * z_interval
    power = z_interval
    total = Interval.zero()
    for index in range(LOG_TERMS):
        total = total + power.scale(Fraction(1, 2 * index + 1))
        power = power * z_squared
    partial = total.scale(Fraction(2))

    exact_tail_upper = (
        2
        * z ** (2 * LOG_TERMS + 1)
        / ((2 * LOG_TERMS + 1) * (1 - z * z))
    )
    return Interval(
        partial.lo,
        (partial + Interval.from_fraction(exact_tail_upper)).hi,
    )


LN2 = _unit_log(Fraction(2))


def _log_point(value: Fraction) -> Interval:
    assert value > 0
    exponent = 0
    reduced = value
    while reduced >= 2:
        reduced /= 2
        exponent += 1
    while reduced < 1:
        reduced *= 2
        exponent -= 1
    return _unit_log(reduced) + LN2.scale(Fraction(exponent))


def natural_log(value: Interval) -> Interval:
    """Monotone enclosure of the natural logarithm."""
    assert value.lo > 0
    lower = _log_point(Fraction(value.lo, SCALE))
    upper = _log_point(Fraction(value.hi, SCALE))
    return Interval(lower.lo, upper.hi)


def _exp_nonnegative_point(value: Fraction) -> Interval:
    assert value >= 0
    reduced = value
    squarings = 0
    while reduced > Fraction(1, 8):
        reduced /= 2
        squarings += 1

    reduced_interval = Interval.from_fraction(reduced)
    term = Interval.one()
    total = Interval.one()
    for index in range(1, EXP_TERMS + 1):
        term = (term * reduced_interval).scale(Fraction(1, index))
        total = total + term

    next_term = (term * reduced_interval).scale(Fraction(1, EXP_TERMS + 1))
    ratio_denominator = 1 - reduced / (EXP_TERMS + 2)
    exact_tail_upper = Fraction(next_term.hi, SCALE) / ratio_denominator
    result = Interval(
        total.lo,
        (total + Interval.from_fraction(exact_tail_upper)).hi,
    )
    for _ in range(squarings):
        result = result * result
    return result


def _exp_point(value: Fraction) -> Interval:
    if value >= 0:
        return _exp_nonnegative_point(value)
    return _exp_nonnegative_point(-value).reciprocal()


def exponential(value: Interval) -> Interval:
    """Monotone enclosure of ``exp``."""
    lower = _exp_point(Fraction(value.lo, SCALE))
    upper = _exp_point(Fraction(value.hi, SCALE))
    return Interval(lower.lo, upper.hi)


@lru_cache(maxsize=None)
def binary_entropy(value: Fraction) -> Interval:
    """Binary entropy in bits at an exact rational probability."""
    assert 0 <= value <= 1
    if value == 0 or value == 1:
        return Interval.zero()
    numerator = -(
        natural_log(Interval.from_fraction(value)).scale(value)
        + natural_log(Interval.from_fraction(1 - value)).scale(1 - value)
    )
    return numerator / LN2


def _softplus2(value: Interval) -> Interval:
    """Enclose ``log2(1 + 2^value)``."""
    power = exponential(value * LN2)
    return natural_log(Interval.one() + power) / LN2


def reconstruct_family() -> tuple[int, ...]:
    return tuple(
        row
        for row in range(1 << DIMENSION)
        if (
            sum((row >> coordinate) & 1 for coordinate in LEFT_COORDINATES),
            sum((row >> coordinate) & 1 for coordinate in RIGHT_COORDINATES),
        )
        in BLOCK_CELLS
    )


def exact_reimer_threshold(size: int) -> int:
    low = 0
    high = size * size
    target = size**size
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def _permute_subset(subset: int, permutation: tuple[int, ...]) -> int:
    image = 0
    for coordinate, target in enumerate(permutation):
        image |= ((subset >> coordinate) & 1) << target
    return image


def order_representatives(
    family: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    """Build and exactly verify the 21 equal S2 x S5 order orbits."""
    all_orders = tuple(permutations(range(DIMENSION)))
    support = frozenset(family)
    automorphisms = tuple(
        permutation
        for permutation in all_orders
        if frozenset(_permute_subset(row, permutation) for row in family)
        == support
    )
    assert len(automorphisms) == 240

    representatives = []
    for left_positions in combinations(range(DIMENSION), len(LEFT_COORDINATES)):
        left_position_set = set(left_positions)
        left_index = 0
        right_index = 0
        order = []
        for position in range(DIMENSION):
            if position in left_position_set:
                order.append(LEFT_COORDINATES[left_index])
                left_index += 1
            else:
                order.append(RIGHT_COORDINATES[right_index])
                right_index += 1
        representatives.append(tuple(order))

    covered: set[tuple[int, ...]] = set()
    for representative in representatives:
        orbit = {
            tuple(automorphism[coordinate] for coordinate in representative)
            for automorphism in automorphisms
        }
        assert len(orbit) == len(automorphisms)
        assert covered.isdisjoint(orbit)
        covered.update(orbit)
    assert covered == set(all_orders)
    return tuple(representatives), len(automorphisms)


def _predecessor_law(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> tuple[tuple[Fraction, Fraction], ...]:
    groups: dict[tuple[int, ...], tuple[int, int]] = {}
    for row in family:
        prefix = tuple((row >> predecessor) & 1 for predecessor in predecessors)
        total, ones = groups.get(prefix, (0, 0))
        groups[prefix] = (total + 1, ones + ((row >> coordinate) & 1))
    size = len(family)
    return tuple(
        (Fraction(ones, total), Fraction(total, size))
        for total, ones in groups.values()
    )


def certified_shapley_iid(
    family: tuple[int, ...],
    dimension: int = DIMENSION,
) -> Interval:
    """Direct predecessor-subset Shapley average, independent of order orbits."""
    total = Interval.zero()
    for coordinate in range(dimension):
        others = tuple(index for index in range(dimension) if index != coordinate)
        for predecessor_count in range(dimension):
            weight = Fraction(
                1,
                dimension * comb(dimension - 1, predecessor_count),
            )
            for predecessors in combinations(others, predecessor_count):
                law = _predecessor_law(family, coordinate, predecessors)
                local = Interval.zero()
                for left, left_mass in law:
                    for right, right_mass in law:
                        probability = left + right - left * right
                        local = local + binary_entropy(probability).scale(
                            left_mass * right_mass
                        )
                total = total + local.scale(weight)
    return total


def _sstar(left: Fraction, right: Fraction) -> Fraction:
    return max(left, right, min(left + right, Fraction(1, 2)))


def relaxed_bellman_cost(
    family: tuple[int, ...],
    order: tuple[int, ...],
) -> tuple[Interval, int]:
    """Upper-bound one fixed-order Bellman optimum without clamp decisions.

    At a nondegenerate state, replace the feasible action interval by ``[0, 1]``.
    If the relaxed child values upper-bound the true children, nonnegative
    feasible transition weights preserve the inequality before this domain
    relaxation.  The unrestricted concave maximum is then ``softplus2(slope)``
    plus the affine constant.
    """
    dimension = len(order)

    @lru_cache(maxsize=None)
    def value(
        position: int,
        left_rows: tuple[int, ...],
        right_rows: tuple[int, ...],
    ) -> Interval:
        if position == dimension:
            return Interval.zero()
        if left_rows > right_rows:
            left_rows, right_rows = right_rows, left_rows

        coordinate = order[position]
        left_zero = tuple(row for row in left_rows if not (row >> coordinate) & 1)
        left_one = tuple(row for row in left_rows if (row >> coordinate) & 1)
        right_zero = tuple(row for row in right_rows if not (row >> coordinate) & 1)
        right_one = tuple(row for row in right_rows if (row >> coordinate) & 1)
        left = Fraction(len(left_one), len(left_rows))
        right = Fraction(len(right_one), len(right_rows))

        children = {
            (left_bit, right_bit): value(
                position + 1,
                left_child,
                right_child,
            )
            for left_bit, left_child in ((0, left_zero), (1, left_one))
            if left_child
            for right_bit, right_child in ((0, right_zero), (1, right_one))
            if right_child
        }

        if len(children) < 4:
            action = _sstar(left, right)
            probabilities = (
                (0, 0, 1 - action),
                (1, 0, action - right),
                (0, 1, action - left),
                (1, 1, left + right - action),
            )
            result = binary_entropy(action)
            for left_bit, right_bit, probability in probabilities:
                assert probability >= 0
                if probability:
                    result = result + children[left_bit, right_bit].scale(
                        probability
                    )
            return result

        slope = (
            -children[0, 0]
            + children[1, 0]
            + children[0, 1]
            - children[1, 1]
        )
        constant = (
            children[0, 0]
            - children[1, 0].scale(right)
            - children[0, 1].scale(left)
            + children[1, 1].scale(left + right)
        )
        return constant + _softplus2(slope)

    root = tuple(family)
    result = value(0, root, root)
    return result, value.cache_info().currsize


def _exact_base(family: tuple[int, ...]) -> dict[str, Any]:
    assert len(family) == 45
    assert family == tuple(sorted(set(family)))
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    assert counts == (18,) * DIMENSION
    incidence = sum(counts)
    assert incidence == 126
    threshold = exact_reimer_threshold(len(family))
    assert threshold == 124
    assert incidence >= threshold
    assert max(counts) == 2 * len(family) // 5

    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    defect = Fraction(missing, len(family) ** 2)
    assert missing == 1600
    assert defect == Fraction(64, 81)

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
        "coordinate_counts": list(counts),
        "total_incidence": incidence,
        "reimer_threshold": threshold,
        "missing_ordered_join_pairs": missing,
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
    representatives, automorphism_count = order_representatives(family)

    iid = certified_shapley_iid(family)
    relaxed_roots = []
    state_count = 0
    for order in representatives:
        root, states = relaxed_bellman_cost(family, order)
        relaxed_roots.append(root)
        state_count += states
    relaxed = Interval.zero()
    for root in relaxed_roots:
        relaxed = relaxed + root
    relaxed = relaxed.scale(Fraction(1, len(relaxed_roots)))

    family_entropy = natural_log(Interval.from_int(len(family))) / LN2
    a_plus_upper = (
        iid.scale(1 - ALPHA)
        + relaxed.scale(ALPHA)
        - family_entropy
    )
    negative_delta = Fraction(1, 40)
    assert a_plus_upper.upper_lt(-negative_delta)

    defect = Fraction(64, 81)
    ratio_lower = negative_delta / defect
    assert ratio_lower == Fraction(81, 2560)
    assert Fraction(17, 81) == 1 - defect
    assert 45**5 < 2**28

    return {
        "schema_version": 1,
        "verdict": "PROVED_GATE_B_UNBOUNDED",
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
        },
        "alpha": {
            "numerator": ALPHA.numerator,
            "denominator": ALPHA.denominator,
            "dimension_independent": True,
        },
        "exact_base": exact_base,
        "symmetry": {
            "automorphism_count": automorphism_count,
            "equal_order_orbit_count": len(representatives),
            "order_orbit_size": automorphism_count,
            "all_order_count": 5040,
            "partition_verified_exactly": True,
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
            "does_not_reuse_arb_clamp_classification": True,
        },
        "rational_consequence": {
            "a_plus_upper_lt": "-1/40",
            "base_ratio_lower_gt": str(ratio_lower),
        },
        "infinite_construction": {
            "family": "B^(boxtimes k)",
            "dimension": "7*k",
            "family_size": "45^k",
            "coordinate_count": "18*45^(k-1)",
            "closure_defect": "1-(17/81)^k",
            "a_plus_upper": "-k/40",
            "repair_ratio_lower": "(k/40)/(1-(17/81)^k) > k/40",
            "reimer_integer_witness": "45^5 < 2^28",
            "zero_defect_convention_used": False,
        },
        "proof_boundary": (
            "The finite negative base is certified here. Product additivity "
            "and all-k admissibility are exact lemmas proved in PROOF.md, not "
            "claims inferred from finite computation."
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
