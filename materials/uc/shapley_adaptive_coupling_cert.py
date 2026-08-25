"""Certified one-sided Bellman gain on the cap-2/5 boundary family.

For

    F5 = {000, 001, 010, 100, 111},

coordinate symmetry makes every one of the six orders identical.  The
one-sided action interval is ``[s*(p,r), min(1,p+r)]``; it excludes the
diagonal identity direction while retaining Cambie's fixed rule as a feasible
policy.  Let ``h=h_2(2/3)`` and ``d=1-h``.  Backward induction gives

    C_plus = h + 3/5 + log2(1 + 2**(1-h)).

Equivalently, relative to the exact fixed cost ``C_seq=151/60``,

    C_plus-C_seq
      = [log2(1+2**d) - 1 - d/2] + (h-5/6)/2.

The first bracket is the root Bernoulli KL advantage.  The second term is the
half-mass stage-two boundary advantage.  This script evaluates the closed form,
the exact-rational iid Shapley term, and the resulting sufficient functional
with 256-bit Arb intervals.

The certificate is a finite theorem about F5, which is not union-closed.  It
certifies the one-sided gain over the already certified sequential value; it
is not an arbitrary-dimensional cap-2/5 proof.

Run:
    ./.venv/bin/python uc/shapley_adaptive_coupling_cert.py
"""

from fractions import Fraction
from itertools import permutations

from flint import arb, ctx

from shapley_direct_cert import (
    ALPHA,
    LOG2,
    ONE,
    ZERO,
    aq,
    certified_direct_value,
)
from shapley_global_coupling_cert import (
    DIMENSION,
    FAMILY,
    certified_shapley_iid,
    h,
)
from shapley_global_coupling import fiber_tables, sstar


ctx.prec = 256


def permute_row(row, order):
    image = 0
    for new_coordinate, old_coordinate in enumerate(order):
        if (row >> old_coordinate) & 1:
            image |= 1 << new_coordinate
    return image


def exact_entropy(value):
    if value == Fraction(1, 2):
        return ONE
    return h(value)


def arb_entropy(value):
    assert value.lower() > 0
    assert value.upper() < 1
    return -(
        value * value.log()
        + (ONE - value) * (ONE - value).log()
    ) / LOG2


def certified_f5_bellman(order):
    """Derive the one-sided Bellman value from every F5 prefix state."""
    tables = fiber_tables(FAMILY, order)
    assert tables[0] == {(): Fraction(2, 5)}
    assert tables[1] == {
        (0,): Fraction(1, 3),
        (1,): Fraction(1, 2),
    }
    assert tables[2] == {
        (0, 0): Fraction(1, 2),
        (0, 1): Fraction(0),
        (1, 0): Fraction(0),
        (1, 1): Fraction(1),
    }

    terminal = {}
    for left_prefix, left in tables[2].items():
        for right_prefix, right in tables[2].items():
            terminal[left_prefix, right_prefix] = exact_entropy(
                sstar(left, right))

    expected_slopes = {
        ((0,), (0,)): ONE,
        ((0,), (1,)): -ONE,
        ((1,), (0,)): -ONE,
        ((1,), (1,)): ZERO,
    }
    middle = {}
    for (left_prefix, right_prefix), expected_slope in expected_slopes.items():
        left = tables[1][left_prefix]
        right = tables[1][right_prefix]
        children = {
            (left_bit, right_bit): terminal[
                left_prefix + (left_bit,),
                right_prefix + (right_bit,),
            ]
            for left_bit in (0, 1)
            for right_bit in (0, 1)
        }
        slope = (
            -children[0, 0]
            + children[1, 0]
            + children[0, 1]
            - children[1, 1]
        )
        assert slope == expected_slope
        greedy = sstar(left, right)
        upper = min(Fraction(1), left + right)
        if slope == ONE:
            assert (left_prefix, right_prefix) == ((0,), (0,))
            action = upper
            assert action == Fraction(2, 3)
        else:
            action = greedy
            assert action == Fraction(1, 2)
        probabilities = (
            (0, 0, 1 - action),
            (1, 0, action - right),
            (0, 1, action - left),
            (1, 1, left + right - action),
        )
        middle[left_prefix, right_prefix] = exact_entropy(action) + sum(
            (aq(probability) * children[left_bit, right_bit]
             for left_bit, right_bit, probability in probabilities),
            ZERO,
        )

    h_two_thirds = h(Fraction(2, 3))
    assert (middle[(0,), (0,)] - (ONE + h_two_thirds)).contains(0)
    assert (middle[(0,), (1,)] - aq(Fraction(3, 2))).contains(0)
    assert (middle[(1,), (0,)] - aq(Fraction(3, 2))).contains(0)
    assert middle[(1,), (1,)] == ONE

    slope = (
        -middle[(0,), (0,)]
        + middle[(1,), (0,)]
        + middle[(0,), (1,)]
        - middle[(1,), (1,)]
    )
    assert (slope - (ONE - h_two_thirds)).contains(0)
    power = (slope * LOG2).exp()
    action = power / (ONE + power)
    assert action.lower() > aq(Fraction(1, 2))
    assert action.upper() < aq(Fraction(4, 5))
    left = aq(Fraction(2, 5))
    value = (
        arb_entropy(action)
        + (ONE - action) * middle[(0,), (0,)]
        + (action - left) * middle[(1,), (0,)]
        + (action - left) * middle[(0,), (1,)]
        + (2 * left - action) * middle[(1,), (1,)]
    )
    return value, action, slope


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    traced = []
    for order in permutations(range(DIMENSION)):
        assert {permute_row(row, order) for row in FAMILY} == set(FAMILY)
        traced.append(certified_f5_bellman(order))

    h_two_thirds = h(Fraction(2, 3))
    slope = ONE - h_two_thirds
    power = (slope * LOG2).exp()
    root_action = power / (ONE + power)
    softplus = (ONE + power).log() / LOG2

    root_advantage = softplus - ONE - slope / 2
    stage_two_advantage = h_two_thirds - aq(Fraction(5, 6))
    sequential = aq(Fraction(151, 60))
    adaptive_from_gains = (
        sequential + root_advantage + stage_two_advantage / 2
    )
    adaptive_closed = (
        h_two_thirds + aq(Fraction(3, 5)) + softplus
    )
    assert (adaptive_from_gains - adaptive_closed).contains(0)
    for traced_value, traced_action, traced_slope in traced:
        assert (traced_value - adaptive_closed).contains(0)
        assert (traced_action - root_action).contains(0)
        assert (traced_slope - slope).contains(0)
    adaptive_recurrence = traced[0][0]

    iid = certified_shapley_iid()
    entropy = arb(len(FAMILY)).log() / LOG2
    adaptive_value = (
        (ONE - ALPHA) * iid + ALPHA * adaptive_recurrence - entropy
    )
    sequential_value = (
        (ONE - ALPHA) * iid + ALPHA * sequential - entropy
    )
    local_value, law_count, vertex_count = certified_direct_value(
        FAMILY, DIMENSION)

    assert root_action.lower() > arb("0.5")
    assert root_action.upper() < arb("0.8")
    assert root_advantage.lower() > arb("0.00057831")
    assert stage_two_advantage.lower() > arb("0.0849625")
    assert (adaptive_recurrence - sequential).lower() > arb("0.04305956")
    assert local_value.upper() < arb("-0.00633442")
    assert sequential_value.lower() > arb("0.00494109")
    assert adaptive_value.lower() > arb("0.00647431")
    assert (adaptive_value - local_value).lower() > arb("0.0128087")

    print("PROVED [exact prefix recurrence + closed form + 256-bit Arb]")
    print("  family =", FAMILY)
    print("  all six orders are equivalent under coordinate symmetry")
    print("  h2(2/3) =", h_two_thirds)
    print("  root continuation slope 1-h =", slope)
    print("  root one-sided Bellman action =", root_action)
    print("  root KL advantage =", root_advantage)
    print("  stage-two advantage before 1/2 occupancy =",
          stage_two_advantage)
    print("  Cseq = 151/60 =", sequential)
    print("  Cplus from prefix recurrence =", adaptive_recurrence)
    print("  Cplus-Cseq =", adaptive_recurrence - sequential)
    print("  Shapley iid Q =", iid)
    print("  H = log2(5) =", entropy)
    print("  local laws / transport vertices =", law_count, "/", vertex_count)
    print("  A_local =", local_value)
    print("  Aseq =", sequential_value)
    print("  Aplus =", adaptive_value)
    print("  Aplus-A_local =", adaptive_value - local_value)
    print()
    print("FINITE VERDICT")
    print("  One-sided Bellman lookahead strictly improves the certified")
    print("  sequential sign crossing at cap 2/5.  No arbitrary-dimensional")
    print("  theorem follows from this finite certificate.")
