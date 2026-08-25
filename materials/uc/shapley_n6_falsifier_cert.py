"""Arb-certified n=6 falsifier of the arbitrary-family A_plus target.

The 25-row family is the union of five cells for the 3+3 coordinate partition:

    (|R intersect A|, |R intersect B|) in
    {(0,0), (0,1), (1,0), (1,2), (2,1)}.

It is simple, active, separating, has every coordinate count 10 (cap 2/5), and
has total incidence 60 > 25 log2(25)/2.  It is not union-closed: 1 and 2 are
members but 3 is not.

This script certifies Q, the fixed sequential cost, and the one-sided Bellman
upper envelope over all order orbits with 256-bit Arb.  At an interval-ambiguous
clamp it encloses both the boundary objective and the unconstrained softplus;
therefore the resulting A_plus upper bound remains valid.  The upper bound is
strictly negative, refuting the arbitrary-family cap-2/5 A_plus theorem.  This
does not refute a UC-restricted fiber theorem.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb

from flint import arb, ctx

from shapley_adaptive_coupling import replay_policy
from shapley_direct import predecessor_law
from shapley_direct_cert import ALPHA, LOG2, ONE, ZERO, aq, h
from shapley_global_coupling import sequential_coupling, sstar


ctx.prec = 256
DIMENSION = 6
FAMILY = (
    0, 1, 2, 4, 8, 11, 13, 14, 16, 19, 21, 22, 25,
    26, 28, 32, 35, 37, 38, 41, 42, 44, 49, 50, 52,
)
ORDERS = tuple(permutations(range(DIMENSION)))
ROW_ONES = tuple(
    sum(1 << row for row in range(64) if (row >> coordinate) & 1)
    for coordinate in range(DIMENSION)
)


def family_mask(family):
    return sum(1 << row for row in family)


def permute_subset(subset, permutation):
    image = 0
    for old_coordinate, new_coordinate in enumerate(permutation):
        if (subset >> old_coordinate) & 1:
            image |= 1 << new_coordinate
    return image


def permute_family(family, permutation):
    return family_mask(tuple(permute_subset(row, permutation) for row in family))


def representative_orders():
    mask = family_mask(FAMILY)
    automorphisms = tuple(
        permutation for permutation in ORDERS
        if permute_family(FAMILY, permutation) == mask
    )
    assert len(automorphisms) == 72
    unseen = set(ORDERS)
    result = []
    while unseen:
        order = min(unseen)
        orbit = {
            tuple(automorphism[coordinate] for coordinate in order)
            for automorphism in automorphisms
        }
        assert len(orbit) == len(automorphisms)
        assert orbit <= unseen
        result.append(order)
        unseen.difference_update(orbit)
    assert len(result) == 10
    return tuple(result)


def certified_iid(law):
    atoms = tuple(atom for atom, _weight in law)
    weights = tuple(weight for _atom, weight in law)
    return sum(
        (
            aq(weights[left] * weights[right])
            * h(atoms[left] + atoms[right] - atoms[left] * atoms[right])
            for left in range(len(law))
            for right in range(len(law))
        ),
        ZERO,
    )


def certified_shapley_iid():
    total = ZERO
    for coordinate in range(DIMENSION):
        others = tuple(
            other for other in range(DIMENSION) if other != coordinate)
        for size in range(DIMENSION):
            weight = Fraction(1, DIMENSION * comb(DIMENSION - 1, size))
            for predecessors in combinations(others, size):
                total += aq(weight) * certified_iid(
                    predecessor_law(FAMILY, coordinate, predecessors))
    return total


def entropy(value):
    if value == 0 or value == 1:
        return ZERO
    return h(value)


def causal_costs(orders, optimize):
    root = family_mask(FAMILY)
    cache = {}
    ambiguous = 0

    def value(suffix, left_rows, right_rows):
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

        slope = None
        constant = None
        if len(children) == 4:
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

        def boundary(action):
            probabilities = (
                (0, 0, 1 - action),
                (1, 0, action - right),
                (0, 1, action - left),
                (1, 1, left + right - action),
            )
            return entropy(action) + sum(
                (aq(probability) * children.get(
                    (left_bit, right_bit), ZERO)
                 for left_bit, right_bit, probability in probabilities),
                ZERO,
            )

        if not optimize or greedy == upper:
            result = boundary(greedy)
        else:
            assert len(children) == 4
            power = (slope * LOG2).exp()
            logistic = power / (ONE + power)
            softplus = constant + (ONE + power).log() / LOG2
            greedy_ball = aq(greedy)
            upper_ball = aq(upper)
            if logistic.upper() < greedy_ball:
                result = boundary(greedy)
            elif logistic.lower() > upper_ball:
                result = boundary(upper)
            elif (logistic.lower() > greedy_ball
                  and logistic.upper() < upper_ball):
                result = softplus
            else:
                ambiguous += 1
                candidates = [softplus]
                if logistic.lower() <= greedy_ball:
                    candidates.append(boundary(greedy))
                if logistic.upper() >= upper_ball:
                    candidates.append(boundary(upper))
                result = candidates[0]
                for candidate in candidates[1:]:
                    result = result.union(candidate)
        cache[key] = result
        return result

    roots = tuple(value(tuple(order), root, root) for order in orders)
    return roots, len(cache), ambiguous


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    assert len(FAMILY) == len(set(FAMILY)) == 25
    counts = tuple(
        sum((row >> coordinate) & 1 for row in FAMILY)
        for coordinate in range(DIMENSION)
    )
    assert counts == (10, 10, 10, 10, 10, 10)
    assert sum(row.bit_count() for row in FAMILY) == 60
    assert 1 in FAMILY and 2 in FAMILY and 3 not in FAMILY
    support = set(FAMILY)
    missing_join_pairs = sum(
        (left | right) not in support
        for left in FAMILY
        for right in FAMILY
    )
    assert missing_join_pairs == 444
    columns = tuple(
        sum(((row >> coordinate) & 1) << index
            for index, row in enumerate(FAMILY))
        for coordinate in range(DIMENSION)
    )
    assert all(columns) and len(set(columns)) == DIMENSION

    orders = representative_orders()
    iid = certified_shapley_iid()
    sequential_roots, sequential_states, sequential_ambiguous = causal_costs(
        orders, optimize=False)
    one_sided_roots, one_sided_states, one_sided_ambiguous = causal_costs(
        orders, optimize=True)
    sequential = sum(sequential_roots, ZERO) / len(sequential_roots)
    one_sided = sum(one_sided_roots, ZERO) / len(one_sided_roots)
    family_entropy = arb(len(FAMILY)).log() / LOG2
    reimer_slack = aq(Fraction(24, 5)) - family_entropy
    sequential_value = (
        (ONE - ALPHA) * iid + ALPHA * sequential - family_entropy)
    one_sided_value = (
        (ONE - ALPHA) * iid + ALPHA * one_sided - family_entropy)
    closure_defect = Fraction(missing_join_pairs, len(FAMILY) ** 2)
    corrected_value = (
        one_sided_value + aq(closure_defect * Fraction(1, 50)))

    assert reimer_slack.lower() > arb("0.15")
    assert sequential_value.upper() < arb("-0.0150")
    assert one_sided_value.upper() < arb("-0.0136")
    assert (one_sided - sequential).lower() > arb("0.0374")
    assert corrected_value.lower() > arb("0.00053")
    assert sequential_ambiguous == 0

    replayed_costs = tuple(
        replay_policy(FAMILY, order)[0] for order in ORDERS)
    replayed_sequential = tuple(
        sequential_coupling(FAMILY, order)[0] for order in ORDERS)
    replayed_average = sum(replayed_costs) / len(replayed_costs)
    replayed_sequential_average = (
        sum(replayed_sequential) / len(replayed_sequential))
    assert abs(replayed_average - 5.003955431250201) < 3e-12
    assert abs(replayed_sequential_average - 4.966516494328278) < 3e-12

    print("PROVED [exact filters/group + 256-bit Arb enclosure]")
    print("  family =", FAMILY)
    print("  block cells = ((0,0),(0,1),(1,0),(1,2),(2,1))")
    print("  coordinate counts =", counts)
    print("  total incidence = 60")
    print("  Reimer slack 24/5-log2(25) =", reimer_slack)
    print("  explicit missing join: 1 union 2 = 3")
    print("  missing join-pair probability =", closure_defect)
    print("  automorphisms / order representatives = 72 /", len(orders))
    print("  sequential / one-sided DAG states =",
          sequential_states, "/", one_sided_states)
    print("  one-sided ambiguous boundary states =", one_sided_ambiguous)
    print("  Q =", iid)
    print("  Cseq =", sequential)
    print("  Cplus upper enclosure =", one_sided)
    print("  Aseq =", sequential_value)
    print("  Aplus upper enclosure =", one_sided_value)
    print("  Cplus-Cseq =", one_sided - sequential)
    print("  Aplus + (1/50) closure defect =", corrected_value)
    print("  all 720 float policy replays: PASS")
    print()
    print("  independent 720-order Cseq / Cplus replay =",
          replayed_sequential_average, "/", replayed_average)
    print("FINITE VERDICT")
    print("  The nontrivial arbitrary-family cap-2/5 Aplus target is refuted")
    print("  in dimension six. UC-specific fiber structure is now mandatory.")
