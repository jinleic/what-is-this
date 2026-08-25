"""One-sided Bellman lookahead for the full-row sequential coupling.

For a fixed coordinate order and joint row prefixes ``u,v``, let

    p = Pr(A_i=1 | A_<i=u),  r = Pr(C_i=1 | C_<i=v).

Cambie's sequential rule uses

    g = s*(p,r) = clamp(1/2, max(p,r), min(1,p+r)).

Unrestricted Bellman optimization over the whole Bernoulli-coupling interval
is vacuous for the union-closed problem: the diagonal policy ``A=C`` is causal,
has cost ``H(A)=log2|F|`` for every family, and is optimal whenever the family
is union-closed.  This script therefore optimizes only in the nonidentity
direction

    g <= s <= min(1,p+r).

If ``V_ac`` is the continuation value after next bits ``a,c``, put
``D=-V_00+V_10+V_01-V_11``.  The one-sided Bellman action is

    s+ = clamp(1 / (1 + 2**(-D)), g, min(1,p+r)).

The fixed rule ``s=g`` remains feasible, so ``C+ >= Cseq`` orderwise.  Every
chosen kernel still preserves both uniform full-row marginals, and

    C+ <= H(A union C) <= log2|F|

for a union-closed family.  Thus ``C+`` is a genuine strengthening of the
fixed sequential sampler, unlike the unrestricted diagonal envelope.

The script performs complete numerical scans of nontrivial n=4 families,
replays selected policies to check both terminal marginals, and optionally
evaluates all 5,172 canonical nontrivial n=5 cap-2/5 representatives from the
existing exact C++ census.  Entropy and Bellman values are float64; family
filters and census membership are exact.

Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python \
        uc/shapley_adaptive_coupling.py \
        [/tmp/shapley_n5_global_coupling_enum]
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import log2
from pathlib import Path
import sys
import time

from shapley_entropy import ALPHA, union_closed
from shapley_global_coupling import (
    TARGET,
    THREE_SEVENTHS,
    TWO_FIFTHS,
    binary_entropy,
    fiber_tables,
    local_order_parts,
    members,
    n4_admissible,
    sequential_coupling,
    sstar,
)
import shapley_n5_global_coupling as n5


ORDERS4 = tuple(permutations(range(4)))


def _sigmoid2(value):
    """Stable base-two logistic function."""
    if value >= 0.0:
        return 1.0 / (1.0 + 2.0 ** (-value))
    power = 2.0 ** value
    return power / (1.0 + power)


def _binary_kl(value, reference):
    """Binary KL divergence in bits."""
    result = 0.0
    if value > 0.0:
        assert reference > 0.0
        result += value * log2(value / reference)
    if value < 1.0:
        assert reference < 1.0
        result += (1.0 - value) * log2(
            (1.0 - value) / (1.0 - reference))
    return result


def _transitions(left, right, union):
    return (
        (0, 0, 1.0 - union),
        (1, 0, union - right),
        (0, 1, union - left),
        (1, 1, left + right - union),
    )


def bellman_coupling(family, order, return_policy=False):
    """Return the optimal one-sided causal cost, optionally with its policy."""
    tables = fiber_tables(family, order)
    policy = {} if return_policy else None

    @lru_cache(maxsize=None)
    def value(position, left_prefix, right_prefix):
        if position == len(order):
            return 0.0

        left_exact = tables[position][left_prefix]
        right_exact = tables[position][right_prefix]
        lower_exact = sstar(left_exact, right_exact)
        upper_exact = min(Fraction(1), left_exact + right_exact)
        left = float(left_exact)
        right = float(right_exact)
        lower = float(lower_exact)
        upper = float(upper_exact)

        children = {}
        for left_bit in (0, 1):
            left_mass = left_exact if left_bit else 1 - left_exact
            if not left_mass:
                continue
            for right_bit in (0, 1):
                right_mass = right_exact if right_bit else 1 - right_exact
                if not right_mass:
                    continue
                children[left_bit, right_bit] = value(
                    position + 1,
                    left_prefix + (left_bit,),
                    right_prefix + (right_bit,),
                )

        if lower_exact < upper_exact:
            assert len(children) == 4
            slope = (
                -children[0, 0]
                + children[1, 0]
                + children[0, 1]
                - children[1, 1]
            )
            union = min(upper, max(lower, _sigmoid2(slope)))
        else:
            slope = 0.0
            union = lower

        continuation = sum(
            probability * children.get((left_bit, right_bit), 0.0)
            for left_bit, right_bit, probability in _transitions(
                left, right, union)
        )
        result = binary_entropy(union) + continuation
        if policy is not None:
            policy[position, left_prefix, right_prefix] = (
                union,
                lower,
                upper,
                slope,
                result,
            )
        return result

    result = value(0, (), ())
    if return_policy:
        return result, policy
    return result


def replay_policy(family, order):
    """Replay the Bellman policy and verify its full-row coupling contract."""
    optimum, policy = bellman_coupling(family, order, return_policy=True)
    tables = fiber_tables(family, order)
    states = {((), ()): 1.0}
    cost = 0.0
    closed = union_closed(family)
    causal_kl = 0.0

    for position in range(len(order)):
        next_states = {}
        for (left_prefix, right_prefix), mass in states.items():
            left = float(tables[position][left_prefix])
            right = float(tables[position][right_prefix])
            union = policy[position, left_prefix, right_prefix][0]
            cost += mass * binary_entropy(union)
            if closed:
                joined_prefix = tuple(
                    left_bit | right_bit
                    for left_bit, right_bit in zip(
                        left_prefix, right_prefix)
                )
                reference = float(tables[position][joined_prefix])
                causal_kl += mass * _binary_kl(union, reference)
            for left_bit, right_bit, probability in _transitions(
                    left, right, union):
                assert probability >= -2e-15
                if probability <= 0.0:
                    continue
                key = (
                    left_prefix + (left_bit,),
                    right_prefix + (right_bit,),
                )
                next_states[key] = next_states.get(key, 0.0) + mass * probability
        assert abs(sum(next_states.values()) - 1.0) < 3e-12
        states = next_states

    ordered_rows = {
        tuple((row >> coordinate) & 1 for coordinate in order)
        for row in family
    }
    uniform = 1.0 / len(family)
    for row in ordered_rows:
        assert abs(sum(
            mass for (left, _right), mass in states.items() if left == row
        ) - uniform) < 4e-12
        assert abs(sum(
            mass for (_left, right), mass in states.items() if right == row
        ) - uniform) < 4e-12

    join_mass = {}
    for (left, right), mass in states.items():
        joined = tuple(a | c for a, c in zip(left, right))
        join_mass[joined] = join_mass.get(joined, 0.0) + mass
    join_entropy = -sum(
        mass * log2(mass) for mass in join_mass.values() if mass > 0.0
    )
    assert abs(cost - optimum) < 4e-12
    assert cost <= join_entropy + 4e-12
    if closed:
        assert abs(log2(len(family)) - cost - causal_kl) < 5e-12
    return cost, join_entropy, states


def greedy_advantage_layers(family, order):
    """Exact performance-difference decomposition in float64 arithmetic."""
    optimum, policy = bellman_coupling(family, order, return_policy=True)
    tables = fiber_tables(family, order)
    states = {((), ()): 1.0}
    greedy_cost = 0.0
    layers = []

    for position in range(len(order)):
        next_states = {}
        layer = 0.0
        for (left_prefix, right_prefix), mass in states.items():
            left_exact = tables[position][left_prefix]
            right_exact = tables[position][right_prefix]
            left = float(left_exact)
            right = float(right_exact)
            greedy = float(sstar(left_exact, right_exact))
            adaptive, _lower, _upper, slope, _value = policy[
                position, left_prefix, right_prefix
            ]
            advantage = (
                binary_entropy(adaptive) + slope * adaptive
                - binary_entropy(greedy) - slope * greedy
            )
            assert advantage >= -2e-14
            layer += mass * advantage
            greedy_cost += mass * binary_entropy(greedy)
            for left_bit, right_bit, probability in _transitions(
                    left, right, greedy):
                if probability <= 0.0:
                    continue
                key = (
                    left_prefix + (left_bit,),
                    right_prefix + (right_bit,),
                )
                next_states[key] = next_states.get(key, 0.0) + mass * probability
        layers.append(layer)
        states = next_states

    assert abs(optimum - greedy_cost - sum(layers)) < 5e-12
    return optimum, greedy_cost, tuple(layers)




def all_order_metrics(family, orders):
    rows = []
    for order in orders:
        iid, _local_coupling, entropy, local_value = local_order_parts(
            family, order)
        sequential = sequential_coupling(family, order)[0]
        adaptive = bellman_coupling(family, order)
        assert adaptive >= sequential - 3e-12
        rows.append((iid, sequential, adaptive, entropy, local_value))

    averages = tuple(
        sum(row[index] for row in rows) / len(rows)
        for index in range(len(rows[0]))
    )
    iid, sequential, adaptive, entropy, local_value = averages
    sequential_value = (
        (1.0 - ALPHA) * iid + ALPHA * sequential - entropy
    )
    adaptive_value = (
        (1.0 - ALPHA) * iid + ALPHA * adaptive - entropy
    )
    return (
        local_value,
        sequential_value,
        adaptive_value,
        adaptive - sequential,
        iid,
        sequential,
        adaptive,
    )


def representative_controls():
    print("ONE-SIDED BELLMAN POLICY CONTROLS")
    controls = (
        ((0, 1, 2, 4, 7), 3),
        ((0, 3, 5, 6, 9, 10, 12), 4),
        ((0, 1, 3), 2),
    )
    for family, dimension in controls:
        orders = tuple(permutations(range(dimension)))
        metrics = all_order_metrics(family, orders)
        for order in orders:
            replay_policy(family, order)
        optimum, greedy, layers = greedy_advantage_layers(family, orders[0])
        print("  family", family)
        print("    A_local / A_seq / A_plus = %.15f / %.15f / %.15f" %
              metrics[:3])
        print("    C_plus-C_seq = %.15f" % metrics[3])
        print("    first-order advantage layers =", layers)
        assert abs(optimum - greedy - metrics[3]) < 5e-12

    assert abs(all_order_metrics(controls[0][0], tuple(
        permutations(range(3))))[2] - 0.006474313148642885) < 3e-12
    assert abs(all_order_metrics(controls[1][0], ORDERS4)[2]
               + 0.14687485799796152) < 3e-12
    assert abs(all_order_metrics(controls[2][0], tuple(
        permutations(range(2))))[3]) < 3e-12
    print("  terminal uniform marginals and cost<=join entropy: PASS")
    print()


def complete_n4_scan():
    records = []
    started = time.time()
    for family_mask in range(1, 1 << 16):
        family = members(family_mask)
        if family == (0,) or not n4_admissible(family, THREE_SEVENTHS):
            continue
        maximum_frequency = Fraction(max(
            sum((row >> coordinate) & 1 for row in family)
            for coordinate in range(4)
        ), len(family))
        records.append((maximum_frequency, all_order_metrics(
            family, ORDERS4), family))

    print("COMPLETE NONTRIVIAL n=4 ONE-SIDED CENSUS (float64 entropy)")
    print("  runtime: %.2f s" % (time.time() - started))

    expected = {TARGET: 110, TWO_FIFTHS: 358, THREE_SEVENTHS: 966}
    for cap, expected_count in expected.items():
        selected = [record for record in records if record[0] <= cap]
        assert len(selected) == expected_count
        worst = min(selected, key=lambda record: record[1][2])
        adaptive_negative = sum(record[1][2] < -1e-12 for record in selected)
        print("  cap", cap, "nontrivial families", len(selected))
        print("    min A_plus = %+.15f at %s" %
              (worst[1][2], worst[2]))
        print("    A_plus negative =", adaptive_negative)

    target = [record for record in records if record[0] <= TARGET]
    two_fifths = [record for record in records
                   if record[0] <= TWO_FIFTHS]
    worst_target = min(target, key=lambda record: record[1][2])
    worst_two_fifths = min(two_fifths, key=lambda record: record[1][2])
    worst_three_sevenths = min(records, key=lambda record: record[1][2])
    assert worst_target[2] == (0, 1, 2, 4, 7, 9, 10, 12)
    assert abs(worst_target[1][2] - 0.14506542577903236) < 3e-12
    assert worst_two_fifths[2] == (0, 1, 2, 4, 7)
    assert abs(worst_two_fifths[1][2] - 0.006474313148642885) < 3e-12
    assert worst_three_sevenths[2] == (0, 3, 5, 6, 9, 10, 12)
    assert abs(worst_three_sevenths[1][2] + 0.14687485799796152) < 3e-12
    assert sum(record[1][2] < -1e-12 for record in records) == 138
    print("  exact family/filter/order coverage: PASS")
    print()
    return records


def complete_n5_scan(binary):
    masks, enumerator_log = n5.read_enumerator(binary)
    counts = Counter(len(n5.members(mask)) for mask in masks)
    assert len(masks) == len(set(masks)) == n5.EXPECTED_CANONICAL_TOTAL
    assert dict(counts) == n5.EXPECTED_CANONICAL
    for mask in masks:
        assert n5.exact_admissible(n5.members(mask))
    assert "total=463343 canonical_total=5172" in enumerator_log

    records = []
    started = time.time()
    for index, mask in enumerate(masks, 1):
        family = n5.members(mask)
        iid = n5.shapley_iid(family)
        adaptive = sum(
            bellman_coupling(family, order) for order in n5.ORDERS
        ) / len(n5.ORDERS)
        entropy = log2(len(family))
        adaptive_value = (
            (1.0 - ALPHA) * iid + ALPHA * adaptive - entropy
        )
        records.append((adaptive_value, family, adaptive, iid))
        if index % 1000 == 0:
            print("  evaluated %d / %d" % (index, len(masks)))

    worst = min(records)
    assert worst[1] == (0, 3, 4, 8, 15)
    assert abs(worst[0] - 0.006474313148642441) < 4e-12
    assert all(record[0] > 1e-10 for record in records)
    assert not any(union_closed(record[1]) for record in records)
    for order in n5.ORDERS:
        replay_policy(worst[1], order)

    print("COMPLETE NONTRIVIAL n=5 CAP-2/5 ONE-SIDED CENSUS (float64)")
    print("  runtime: %.2f s" % (time.time() - started))
    print("  nontrivial canonical orbits:", len(records))
    print("  represented nontrivial labeled families:",
          n5.EXPECTED_LABELED_TOTAL)
    print("  nontrivial negative / zero / positive = 0 / 0 /", len(records))
    sequential_value, _iid, sequential = n5.sequential_value(worst[1])
    print("  nontrivial minimum")
    print("    family =", worst[1])
    print("    A_seq / A_plus = %.15f / %.15f" %
          (sequential_value, worst[0]))
    print("    C_plus-C_seq = %.15f" % (worst[2] - sequential))
    print("  exact census gates and worst-policy replay: PASS")
    print()
    return records


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    representative_controls()
    complete_n4_scan()
    if len(sys.argv) > 1:
        binary_path = Path(sys.argv[1])
        if not binary_path.is_file():
            raise SystemExit("missing n=5 enumerator binary")
        complete_n5_scan(binary_path)
    else:
        print("n=5 scan skipped; pass the exact enumerator binary to enable it.")
