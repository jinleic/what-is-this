"""Bounded structured n=6 falsification probe for the cap-2/5 frontier.

The complete nontrivial arbitrary-family n=5 census is positive, so n=6 is
the first unexamined dimension of the nontrivial entropy relaxation.  This
script tests two cheap, predeclared descendants of the first cap-3/7
sequential failure:

* all 127 unions of complete Hamming layers on six coordinates, with the
  trivial equality case reported separately; and
* the minimal one-new-coordinate dilution of
  ``{empty} union binom([4],2)`` from frequency 3/7 to 3/8.

Full Hamming-layer families are S_6-invariant, so one order is the exact
all-order average.  The dilution control evaluates all 5!=120 orders.  Filters
are exact except for the float64 evaluation of the Reimer logarithm; reported
entropy and Bellman values are float64.  This is a bounded falsification probe,
not a complete n=6 census or a union-closed theorem.

Run:
    ./.venv/bin/python uc/shapley_n6_structured_probe.py
"""

from itertools import permutations
from math import log2

from shapley_adaptive_coupling import bellman_coupling
from shapley_entropy import ALPHA
from shapley_global_coupling import sequential_coupling
from shapley_join_loss import shapley_iid


def reimer_admissible(family):
    incidence = sum(row.bit_count() for row in family)
    return 2.0 * incidence >= len(family) * log2(len(family)) - 1e-12


def value(family, dimension, orders):
    iid = shapley_iid(family, dimension)
    sequential = sum(
        sequential_coupling(family, order)[0] for order in orders
    ) / len(orders)
    one_sided = sum(
        bellman_coupling(family, order) for order in orders
    ) / len(orders)
    entropy = log2(len(family))
    return (
        (1.0 - ALPHA) * iid + ALPHA * sequential - entropy,
        (1.0 - ALPHA) * iid + ALPHA * one_sided - entropy,
        one_sided - sequential,
    )


def hamming_layer_scan():
    records = []
    trivial = 0
    order = tuple(range(6))
    for layer_mask in range(1, 1 << 7):
        family = tuple(
            row for row in range(1 << 6)
            if (layer_mask >> row.bit_count()) & 1
        )
        if family == (0,):
            assert reimer_admissible(family)
            trivial += 1
            continue
        counts = tuple(
            sum((row >> coordinate) & 1 for row in family)
            for coordinate in range(6)
        )
        if 5 * max(counts) > 2 * len(family):
            continue
        if not reimer_admissible(family):
            continue
        layers = tuple(
            weight for weight in range(7)
            if (layer_mask >> weight) & 1
        )
        records.append((value(family, 6, (order,)), layers, len(family)))

    assert trivial == 1
    assert len(records) == 5
    worst = min(records, key=lambda record: record[0][1])
    assert worst[1] == (2, 6)
    assert worst[2] == 16
    assert abs(worst[0][0] - 0.19395266706789815) < 3e-12
    assert abs(worst[0][1] - 0.19506114061458835) < 3e-12
    assert all(record[0][1] > 0.19 for record in records)
    return records, trivial


def diluted_first_failure():
    family = (0, 3, 5, 6, 9, 10, 12, 16)
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(5)
    )
    assert counts == (3, 3, 3, 3, 1)
    assert 5 * max(counts) <= 2 * len(family)
    assert reimer_admissible(family)
    result = value(family, 5, tuple(permutations(range(5))))
    assert abs(result[0] - 0.33713454898389517) < 3e-12
    assert abs(result[1] - 0.338994827344125) < 3e-12
    return family, result


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    records, trivial = hamming_layer_scan()
    print("COMPLETE 127-MASK S_6 HAMMING-LAYER SCREEN")
    print("  trivial cap/Reimer equality cases:", trivial)
    print("  nontrivial cap/Reimer admissible:", len(records))
    for metrics, layers, size in sorted(records, key=lambda record: record[1]):
        print("  layers=%s m=%d A_seq=%+.12f A_plus=%+.12f gain=%+.12f" %
              (layers, size, metrics[0], metrics[1], metrics[2]))
    family, metrics = diluted_first_failure()
    print()
    print("ALL-ORDER MINIMAL-DILUTION CONTROL")
    print("  family =", family)
    print("  A_seq / A_plus / gain = %.15f / %.15f / %.15f" % metrics)
    print()
    print("VERDICT")
    print("  No structured n=6 counterexample found; this does not cover")
    print("  arbitrary n=6 families.")
