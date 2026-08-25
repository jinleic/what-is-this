"""FULL-UNIFORM sequential-coupling frontier for the entropy sampler.

Cambie's dependent strategy samples two uniform rows A,C sequentially.  At a
joint prefix, write

    p = Pr(A_i=1 | A_prefix),   r = Pr(C_i=1 | C_prefix).

The Bernoulli coupling is chosen so that the conditional union probability is

    s*(p,r) = max(p, r, min(p+r, 1/2)).

Unlike ``shapley_direct.py``, which minimizes the coupling cost independently
for every coordinate and predecessor set, this file propagates the exact joint
prefix law.  The couplings at different coordinates therefore have to coexist
inside one coupling of two full uniform rows.  This is the global consistency
that projection-local invariants discard.

For an order pi define Q_pi using independent uniform rows and Cseq_pi using
the sequential coupling above.  Chain rule and data processing give

    Q_pi    <= H(A union B),
    Cseq_pi <= H(A union C).

If the family is union-closed, both output entropies are at most log2|F|.
Consequently a proof that

    E_pi[(1-alpha) Q_pi + alpha Cseq_pi] > log2|F|

under a low-frequency hypothesis would directly contradict union closure.  It
is a different sufficient route from nonnegativity of the pessimistic local-min
functional.  The balanced-label projection refutation does not touch it because
this process runs on the original full uniform rows.

This experiment performs:

* exact-Fraction propagation and exact uniform-marginal checks for every order;
* complete numerical scans of all nontrivial n=4 uniform families satisfying
  Reimer and caps 0.38261, 2/5, and 3/7;
* complete all-order controls on the two n=5 normalized-Shapley extremizers;
* a seeded 24-order stress test on the 120-row balanced-label family; and
* a negative control showing that simply adding the k=3 iid-union entropy term
  is worse at Cambie's two-atom equality obstruction.

Evidence labels:

* COMPLETE / EXACT: finite-family enumeration, cap/incidence gates, sequential
  transition probabilities, final uniform marginals, and all-order coverage.
* COMPLETE NUMERICAL: entropy values and local transport minima (float64
  logarithms and HiGHS).
* SEEDED NUMERICAL: the 120-row, 13-coordinate 24-order stress test.

The complete n=4 result is a genuine finite breakthrough, not a theorem for
arbitrary dimension: the globally consistent sampler rescues all 11 local-min
failures at cap 2/5.  It first fails in the next possible n=4 frequency band,
3/7, on {empty} plus all six two-subsets of [4].

Run:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python \
        uc/shapley_global_coupling.py
"""

from fractions import Fraction
from itertools import permutations
from math import log2
from random import Random

from shapley_contraction_refute import build_family
from shapley_entropy import ALPHA, functional, law_key, members

HALF = Fraction(1, 2)
TARGET = Fraction(38261, 100000)
TWO_FIFTHS = Fraction(2, 5)
THREE_SEVENTHS = Fraction(3, 7)
N4_INCIDENCE_NEED = {
    1: 0, 2: 1, 3: 3, 4: 4, 5: 6, 6: 8, 7: 10, 8: 12,
    9: 15, 10: 17, 11: 20, 12: 22, 13: 25, 14: 27, 15: 30, 16: 32,
}
N5_NONSEPARATING = (0, 1, 2, 4, 7, 25, 26, 28)
N5_SEPARATING = (0, 1, 6, 10, 13, 19, 20, 24)


def binary_entropy(value):
    value = float(value)
    if value <= 0.0 or value >= 1.0:
        return 0.0
    return -value * log2(value) - (1.0 - value) * log2(1.0 - value)


def sstar(left, right):
    return max(left, right, min(left + right, HALF))


def fiber_tables(family, order):
    """Exact next-bit probabilities for every occurring prefix."""
    tables = []
    for position, coordinate in enumerate(order):
        counts = {}
        for row in family:
            prefix = tuple((row >> order[j]) & 1 for j in range(position))
            total, ones = counts.get(prefix, (0, 0))
            counts[prefix] = (total + 1, ones + ((row >> coordinate) & 1))
        tables.append({
            prefix: Fraction(ones, total)
            for prefix, (total, ones) in counts.items()
        })
    return tuple(tables)


def entropy_of_mass(mass):
    return -sum(
        float(weight) * log2(float(weight))
        for weight in mass.values()
        if weight
    )


def sequential_coupling(family, order):
    """Return Cseq, join entropy, and the exact final row-pair coupling."""
    tables = fiber_tables(family, order)
    states = {((), ()): Fraction(1)}
    cost = 0.0

    for position, _coordinate in enumerate(order):
        next_states = {}
        for (left_prefix, right_prefix), mass in states.items():
            left = tables[position][left_prefix]
            right = tables[position][right_prefix]
            union = sstar(left, right)
            both = left + right - union
            transitions = (
                (0, 0, 1 - union),
                (1, 0, union - right),
                (0, 1, union - left),
                (1, 1, both),
            )
            assert all(probability >= 0 for _, _, probability in transitions)
            assert sum(probability for _, _, probability in transitions) == 1
            assert sum(probability for a, _, probability in transitions if a) == left
            assert sum(probability for _, c, probability in transitions if c) == right
            assert sum(probability for a, c, probability in transitions if a or c) \
                == union
            cost += float(mass) * binary_entropy(union)
            for left_bit, right_bit, probability in transitions:
                if not probability:
                    continue
                key = (
                    left_prefix + (left_bit,),
                    right_prefix + (right_bit,),
                )
                next_states[key] = (
                    next_states.get(key, Fraction()) + mass * probability
                )
        assert sum(next_states.values(), Fraction()) == 1
        states = next_states

    ordered_rows = {
        tuple((row >> coordinate) & 1 for coordinate in order)
        for row in family
    }
    uniform_mass = Fraction(1, len(family))
    assert len(ordered_rows) == len(family)
    for row in ordered_rows:
        assert sum(mass for (left, _), mass in states.items() if left == row) \
            == uniform_mass
        assert sum(mass for (_, right), mass in states.items() if right == row) \
            == uniform_mass

    join_mass = {}
    for (left, right), mass in states.items():
        joined = tuple(a | b for a, b in zip(left, right))
        join_mass[joined] = join_mass.get(joined, Fraction()) + mass
    join_entropy = entropy_of_mass(join_mass)
    assert cost <= join_entropy + 2e-12
    return cost, join_entropy, states


def local_order_parts(family, order):
    """Return Q, independently minimized C, L, and their local functional."""
    independent = 0.0
    coupled_minimum = 0.0
    entropy = 0.0
    local_value = 0.0
    for position, coordinate in enumerate(order):
        key = law_key(family, order, position)
        value, _mean, local_entropy, _atoms, _weights = functional(key)
        iid = sum(
            float(left_mass * right_mass)
            * binary_entropy(left + right - left * right)
            for left, left_mass in key
            for right, right_mass in key
        )
        coupled = (
            value + local_entropy - (1.0 - ALPHA) * iid
        ) / ALPHA
        independent += iid
        coupled_minimum += coupled
        entropy += local_entropy
        local_value += value
    assert abs(entropy - log2(len(family))) < 3e-12
    return independent, coupled_minimum, entropy, local_value


def all_order_metrics(family, dimension):
    rows = []
    for order in permutations(range(dimension)):
        independent, coupled_minimum, entropy, local_value = local_order_parts(
            family, order)
        sequential, join_entropy, _states = sequential_coupling(family, order)
        assert sequential >= coupled_minimum - 3e-10
        sequential_value = (
            (1.0 - ALPHA) * independent + ALPHA * sequential - entropy
        )
        rows.append((
            independent,
            coupled_minimum,
            sequential,
            entropy,
            local_value,
            sequential_value,
            join_entropy,
        ))
    return tuple(
        sum(row[index] for row in rows) / len(rows)
        for index in range(len(rows[0]))
    )


def n4_admissible(family, cap):
    size = len(family)
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(4)
    )
    return (
        max(counts) * cap.denominator <= cap.numerator * size
        and sum(counts) >= N4_INCIDENCE_NEED[size]
    )


def complete_n4_scan():
    records = []
    for family_mask in range(1, 1 << 16):
        family = members(family_mask)
        if family == (0,) or not n4_admissible(family, THREE_SEVENTHS):
            continue
        metrics = all_order_metrics(family, 4)
        maximum_frequency = Fraction(
            max(
                sum((row >> coordinate) & 1 for row in family)
                for coordinate in range(4)
            ),
            len(family),
        )
        records.append((maximum_frequency, metrics, family))

    def report(cap, expected_count):
        selected = [record for record in records if record[0] <= cap]
        assert len(selected) == expected_count
        worst_local = min(selected, key=lambda record: record[1][4])
        worst_sequential = min(selected, key=lambda record: record[1][5])
        local_negative = sum(record[1][4] < -1e-12 for record in selected)
        sequential_negative = sum(record[1][5] < -1e-12 for record in selected)
        rescued = sum(
            record[1][4] < -1e-12 and record[1][5] >= -1e-12
            for record in selected
        )
        print("  cap", cap, "nontrivial families", len(selected))
        print("    minimum local A     = %.15f at %s" %
              (worst_local[1][4], worst_local[2]))
        print("    minimum sequential A= %.15f at %s" %
              (worst_sequential[1][5], worst_sequential[2]))
        print("    local-negative / sequential-negative / rescued = %d / %d / %d" %
              (local_negative, sequential_negative, rescued))
        return worst_local, worst_sequential, local_negative, sequential_negative, rescued

    target = report(TARGET, 110)
    two_fifths = report(TWO_FIFTHS, 358)
    three_sevenths = report(THREE_SEVENTHS, 966)

    assert abs(target[1][1][5] - 0.14271122065723985) < 2e-12
    assert two_fifths[2:] == (11, 0, 11)
    assert abs(two_fifths[1][1][5] - 0.004941095444916943) < 2e-12
    assert three_sevenths[2:] == (173, 138, 35)
    assert three_sevenths[1][2] == (0, 3, 5, 6, 9, 10, 12)
    assert abs(three_sevenths[1][1][5] + 0.14797589694926927) < 2e-12
    return records


def representative_controls():
    print("COMPLETE ALL-ORDER CONTROLS")
    for name, family, dimension in (
        ("n4 target minimizer", (0, 1, 2, 4, 7, 9, 10, 12), 4),
        ("n5 nonseparating extremizer", N5_NONSEPARATING, 5),
        ("n5 separating extremizer", N5_SEPARATING, 5),
    ):
        metrics = all_order_metrics(family, dimension)
        print(" ", name)
        print("    Q / Cmin / Cseq / H = %.12f / %.12f / %.12f / %.12f" %
              (metrics[0], metrics[1], metrics[2], metrics[3]))
        print("    local A / sequential A / consistency gain = "
              "%.12f / %.12f / %.12f" %
              (metrics[4], metrics[5], metrics[5] - metrics[4]))
        assert metrics[5] > metrics[4] > 0


def balanced_label_stress():
    family = build_family()
    assert len(family) == len(set(family)) == 120
    rng = Random(20260824)
    orders = [tuple(rng.sample(range(13), 13)) for _ in range(24)]
    rows = []
    for order in orders:
        independent, coupled_minimum, entropy, local_value = local_order_parts(
            family, order)
        sequential, _join_entropy, _states = sequential_coupling(family, order)
        assert sequential >= coupled_minimum - 3e-10
        sequential_value = (
            (1.0 - ALPHA) * independent + ALPHA * sequential - entropy
        )
        rows.append((local_value, sequential_value, sequential - coupled_minimum))
    averages = tuple(sum(row[i] for row in rows) / len(rows) for i in range(3))
    print("SEEDED NUMERICAL BALANCED-LABEL STRESS")
    print("  family / coordinates / orders = 120 / 13 / 24")
    print("  mean local A / sequential A = %.12f / %.12f" % averages[:2])
    print("  mean C consistency gap = %.12f" % averages[2])
    print("  sequential-A range = [%.12f, %.12f]" %
          (min(row[1] for row in rows), max(row[1] for row in rows)))
    assert min(row[1] for row in rows) > 1.2
    assert averages[1] > averages[0]


def higher_iid_negative_control():
    # Cambie's numerical equality obstruction.  Q2=L and C=L by its defining
    # equation; the displayed k=3 deficit makes any positive replacement of
    # existing strategy weight by Q3 immediately worse at this obstruction.
    b = 0.32945473850303697239
    a = 1.0 - 1.0 / (2.0 - binary_entropy(b))
    entropy = (1.0 - a) * binary_entropy(b)
    values = {}
    for k in range(2, 9):
        qk = (1.0 - a) ** k * binary_entropy(1.0 - (1.0 - b) ** k)
        values[k] = qk - entropy
    print("HIGHER-IID-UNION NEGATIVE CONTROL [float64]")
    for k, difference in values.items():
        print("  k=%d: Q_k-L = %+.15f" % (k, difference))
    assert abs(values[2]) < 2e-15
    assert values[3] < -0.15


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    higher_iid_negative_control()
    print()
    representative_controls()
    print()
    print("COMPLETE n=4 REIMER + FREQUENCY SCAN [float64 entropy/HiGHS]")
    complete_n4_scan()
    print()
    balanced_label_stress()
    print()
    print("FRONTIER VERDICT")
    print("  PROMISING: globally consistent sequential coupling strictly improves")
    print("  every tested representative and rescues all 11 local-min failures")
    print("  through cap 2/5 in the complete n=4 relaxation.")
    print("  BOUNDARY: it fails at the next n=4 frequency band, 3/7; this is")
    print("  finite evidence, not an arbitrary-n theorem or a union-closed proof.")
    print("  REFUTED: simply adding the k=3 iid-union entropy term at the old")
    print("  equality obstruction.  The next proof target is a global consistency")
    print("  gap or telescoping identity for the sequential prefix coupling.")
