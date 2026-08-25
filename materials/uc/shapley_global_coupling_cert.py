"""Certified global-coupling rescue of the first local-minimum failure.

The complete scan in ``shapley_global_coupling.py`` finds that the local-min
Shapley functional first becomes negative at frequency 2/5 on

    F5 = {000, 001, 010, 100, 111}.

The globally consistent sequential Sawin coupling reverses the sign.  This file
certifies both claims.  All six active-coordinate orders have exact sequential
cost

    Cseq = 1 + 1 + 31/60 = 151/60.

The iid term is evaluated from exact rational posterior laws with 256-bit Arb;
the independently minimized coupling value is certified by exhaustive rational
transport vertices through ``shapley_direct_cert.py``.

This is a finite theorem about F5, not a union-closed theorem.  F5 is not
union-closed.  Its role is to prove that global prefix-coupling consistency is
strictly stronger than minimizing every local coupling separately and can cross
a sign boundary that the previous frontier could not cross.

Run:
    ./.venv/bin/python uc/shapley_global_coupling_cert.py
"""

from fractions import Fraction
from itertools import combinations, permutations

from flint import arb, ctx

from shapley_direct import predecessor_law, predecessor_weight
from shapley_direct_cert import (
    ALPHA,
    LOG2,
    ONE,
    ZERO,
    aq,
    certified_direct_value,
    h,
)
from shapley_global_coupling import HALF, fiber_tables, sstar

ctx.prec = 256
FAMILY = (0, 1, 2, 4, 7)
DIMENSION = 3


def certified_iid(law):
    atoms = tuple(atom for atom, _ in law)
    weights = tuple(weight for _, weight in law)
    return sum(
        (
            aq(weights[i] * weights[j])
            * h(atoms[i] + atoms[j] - atoms[i] * atoms[j])
            for i in range(len(law))
            for j in range(len(law))
        ),
        ZERO,
    )


def certified_shapley_iid():
    total = ZERO
    for coordinate in range(DIMENSION):
        others = tuple(i for i in range(DIMENSION) if i != coordinate)
        for size in range(DIMENSION):
            coefficient = predecessor_weight(DIMENSION, size)
            for predecessors in combinations(others, size):
                law = predecessor_law(FAMILY, coordinate, predecessors)
                total += aq(coefficient) * certified_iid(law)
    return total


def exact_sequential_half_mass(order):
    """Return Cseq exactly when all positive-entropy s* values equal 1/2."""
    tables = fiber_tables(FAMILY, order)
    states = {((), ()): Fraction(1)}
    half_mass = Fraction()
    layers = []
    for position in range(DIMENSION):
        next_states = {}
        layer_half_mass = Fraction()
        for (left_prefix, right_prefix), mass in states.items():
            left = tables[position][left_prefix]
            right = tables[position][right_prefix]
            union = sstar(left, right)
            assert union in (Fraction(), HALF, Fraction(1))
            if union == HALF:
                half_mass += mass
                layer_half_mass += mass
            both = left + right - union
            for left_bit, right_bit, probability in (
                (0, 0, 1 - union),
                (1, 0, union - right),
                (0, 1, union - left),
                (1, 1, both),
            ):
                if probability:
                    key = (
                        left_prefix + (left_bit,),
                        right_prefix + (right_bit,),
                    )
                    next_states[key] = (
                        next_states.get(key, Fraction()) + mass * probability
                    )
        layers.append(layer_half_mass)
        states = next_states
    assert sum(states.values(), Fraction()) == 1
    return half_mass, tuple(layers)


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()

    sequential_values = set()
    layer_values = set()
    for order in permutations(range(DIMENSION)):
        value, layers = exact_sequential_half_mass(order)
        sequential_values.add(value)
        layer_values.add(layers)
    assert sequential_values == {Fraction(151, 60)}
    assert layer_values == {(Fraction(1), Fraction(1), Fraction(31, 60))}

    iid = certified_shapley_iid()
    entropy = arb(len(FAMILY)).log() / LOG2
    sequential = aq(Fraction(151, 60))
    sequential_value = (
        (ONE - ALPHA) * iid + ALPHA * sequential - entropy
    )
    local_value, law_count, vertex_count = certified_direct_value(
        FAMILY, DIMENSION)

    assert local_value.upper() < arb("-0.00633442")
    assert sequential_value.lower() > arb("0.00494109")
    assert (sequential_value - local_value).lower() > arb("0.0112755")

    print("PROVED [exact coupling propagation + rational transport vertices + Arb]")
    print("  family =", FAMILY)
    print("  maximum frequency = 2/5; Reimer slack = 12/5-log2(5) > 0")
    print("  all six sequential layer costs = (1, 1, 31/60)")
    print("  Cseq = 151/60")
    print("  Shapley iid Q =", iid)
    print("  H = log2(5) =", entropy)
    print("  local laws / transport vertices =", law_count, "/", vertex_count)
    print("  pessimistic local-min A =", local_value, "< 0")
    print("  globally consistent sequential A =", sequential_value, "> 0")
    print("  certified sign-crossing gain =", sequential_value - local_value)
    print()
    print("FINITE VERDICT")
    print("  Global sequential consistency repairs the first 2/5 local-minimum")
    print("  obstruction.  This validates the new mechanism on an exact boundary")
    print("  example; it does not prove the arbitrary-dimensional inequality.")
