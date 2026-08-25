"""Shared-prefix DAG evaluator for n=6 causal coupling costs.

The original evaluator rebuilds prefix tables and a private dynamic program for
each coordinate order.  A Bellman state actually depends only on

    (remaining coordinate suffix, left row fiber, right row fiber).

Representing row fibers as 64-bit masks over the absolute row universe lets all
720 orders share the same memoized DAG.  Canonicalizing the two fibers under
left/right exchange is exact because the coupling objective is symmetric.

This module provides both the one-sided Bellman optimum and the fixed s* policy.
Probabilities and action bounds are exact Fractions; entropy and continuation
values are float64, matching ``shapley_adaptive_coupling.py``.
"""

from fractions import Fraction

from shapley_global_coupling import binary_entropy, sstar


DIMENSION = 6
ROW_ONES = tuple(
    sum(1 << row for row in range(1 << DIMENSION)
        if (row >> coordinate) & 1)
    for coordinate in range(DIMENSION)
)


def family_mask(family):
    return sum(1 << row for row in family)


def _sigmoid2(value):
    if value >= 0.0:
        return 1.0 / (1.0 + 2.0 ** (-value))
    power = 2.0 ** value
    return power / (1.0 + power)


def coupling_costs(family, orders, optimize):
    """Return costs for all orders and the number of shared DAG states."""
    root = family_mask(family)
    cache = {}

    def value(suffix, left_rows, right_rows):
        if not suffix:
            return 0.0
        if left_rows > right_rows:
            left_rows, right_rows = right_rows, left_rows
        key = (suffix, left_rows, right_rows)
        cached = cache.get(key)
        if cached is not None:
            return cached

        coordinate = suffix[0]
        rest = suffix[1:]
        ones = ROW_ONES[coordinate]
        left_one = left_rows & ones
        right_one = right_rows & ones
        left_zero = left_rows ^ left_one
        right_zero = right_rows ^ right_one
        left = Fraction(left_one.bit_count(), left_rows.bit_count())
        right = Fraction(right_one.bit_count(), right_rows.bit_count())
        greedy = sstar(left, right)
        upper = min(Fraction(1), left + right)

        children = {}
        for left_bit, left_child in ((0, left_zero), (1, left_one)):
            if not left_child:
                continue
            for right_bit, right_child in ((0, right_zero), (1, right_one)):
                if not right_child:
                    continue
                children[left_bit, right_bit] = value(
                    rest, left_child, right_child)

        if optimize and greedy < upper:
            assert len(children) == 4
            slope = (
                -children[0, 0]
                + children[1, 0]
                + children[0, 1]
                - children[1, 1]
            )
            union = min(float(upper), max(float(greedy), _sigmoid2(slope)))
        else:
            union = float(greedy)

        left_float = float(left)
        right_float = float(right)
        continuation = sum(
            probability * children.get((left_bit, right_bit), 0.0)
            for left_bit, right_bit, probability in (
                (0, 0, 1.0 - union),
                (1, 0, union - right_float),
                (0, 1, union - left_float),
                (1, 1, left_float + right_float - union),
            )
        )
        result = binary_entropy(union) + continuation
        cache[key] = result
        return result

    costs = tuple(value(tuple(order), root, root) for order in orders)
    return costs, len(cache)


def one_sided_costs(family, orders):
    return coupling_costs(family, orders, optimize=True)


def sequential_costs(family, orders):
    return coupling_costs(family, orders, optimize=False)
