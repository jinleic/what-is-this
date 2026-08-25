"""COMPLETE n=5 cap-2/5 global sequential-coupling experiment.

The C++ companion enumerates exactly every nontrivial uniform simple family on
[5] satisfying

    max_i Pr(X_i=1) <= 2/5,
    2 E|X| >= log2|F|,

and emits one canonical representative of each S_5 coordinate-permutation
orbit.  There are 463,343 labeled families and 5,172 canonical orbits.

For every representative this script computes

    A_seq = E_pi[(1-alpha) Q_pi + alpha Cseq_pi] - log2|F|,

where Q uses independent uniform rows and Cseq propagates Cambie's dependent
Bernoulli coupling through one globally consistent joint-prefix process.  All
5!=120 orders are evaluated.  Coordinate-permutation invariance transfers the
comparison from canonical representatives to every labeled family.

Result: every one of the 463,343 families has A_seq>0 in float64.  The minimum
representative is

    (0,3,4,8,15),

a one-duplicate-coordinate plus one-dummy-coordinate lift of
{000,001,010,100,111}, with A_seq=0.0049410954449... .  The companion
``shapley_global_coupling_cert.py`` certifies the three-coordinate base
family's sign crossing with exact Fraction propagation, exhaustive rational
transport vertices, and 256-bit Arb.  It does not certify the five-coordinate
representative or its numerical minimality; those remain COMPLETE NUMERICAL.

No nontrivial enumerated family is union-closed.  The result is a complete
numerical n=5 entropy-relaxation census, not an arbitrary-dimensional
union-closed proof.

Build and run:
    c++ -O3 -std=c++17 uc/shapley_n5_global_coupling_enumerate.cpp \
        -o /tmp/shapley_n5_global_coupling_enum
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python \
        uc/shapley_n5_global_coupling.py \
        /tmp/shapley_n5_global_coupling_enum
"""

from collections import Counter
from itertools import combinations, permutations
from math import comb, log2
from pathlib import Path
import subprocess
import sys
import time

from shapley_direct import predecessor_law
from shapley_entropy import ALPHA, union_closed
from shapley_global_coupling import binary_entropy, fiber_tables

SPECS = {
    3: (1, 3),
    4: (1, 4),
    5: (2, 6),
    6: (2, 8),
    7: (2, 10),
    8: (3, 12),
    9: (3, 15),
    10: (4, 17),
    11: (4, 20),
    13: (5, 25),
    15: (6, 30),
}
EXPECTED_LABELED = {
    3: 145,
    4: 70,
    5: 6187,
    6: 2833,
    7: 296,
    8: 44120,
    9: 2910,
    10: 290488,
    11: 14160,
    13: 38770,
    15: 63364,
}
EXPECTED_CANONICAL = {
    3: 9,
    4: 5,
    5: 126,
    6: 62,
    7: 10,
    8: 573,
    9: 48,
    10: 3058,
    11: 179,
    13: 433,
    15: 669,
}
EXPECTED_LABELED_TOTAL = 463343
EXPECTED_CANONICAL_TOTAL = 5172
ORDERS = tuple(permutations(range(5)))


def members(mask):
    return tuple(subset for subset in range(32) if (mask >> subset) & 1)


def family_mask(family):
    return sum(1 << subset for subset in family)


def exact_admissible(family):
    coordinate_cap, incidence_need = SPECS[len(family)]
    counts = tuple(
        sum((subset >> coordinate) & 1 for subset in family)
        for coordinate in range(5)
    )
    return max(counts) <= coordinate_cap and sum(counts) >= incidence_need


def permute_subset(subset, permutation):
    image = 0
    for old, new in enumerate(permutation):
        if (subset >> old) & 1:
            image |= 1 << new
    return image


def permute_family(mask, permutation):
    image = 0
    for subset in range(32):
        if (mask >> subset) & 1:
            image |= 1 << permute_subset(subset, permutation)
    return image


def canonical(mask):
    return min(permute_family(mask, permutation) for permutation in ORDERS)


def brute_canonical(size):
    found = set()
    for family in combinations(range(32), size):
        if exact_admissible(family):
            mask = family_mask(family)
            if canonical(mask) == mask:
                found.add(mask)
    return found


def read_enumerator(binary):
    result = subprocess.run(
        [str(binary)], capture_output=True, text=True, check=True)
    masks = tuple(int(line, 16) for line in result.stdout.splitlines())
    return masks, result.stderr


def iid_local_cost(law):
    return sum(
        float(left_mass * right_mass)
        * binary_entropy(left + right - left * right)
        for left, left_mass in law
        for right, right_mass in law
    )


def shapley_iid(family):
    total = 0.0
    for coordinate in range(5):
        others = tuple(i for i in range(5) if i != coordinate)
        for size in range(5):
            weight = 1.0 / (5 * comb(4, size))
            for predecessors in combinations(others, size):
                total += weight * iid_local_cost(
                    predecessor_law(family, coordinate, predecessors)
                )
    return total


def sequential_cost_float(family, order):
    tables = fiber_tables(family, order)
    states = {((), ()): 1.0}
    cost = 0.0
    for position in range(5):
        next_states = {}
        for (left_prefix, right_prefix), mass in states.items():
            left = float(tables[position][left_prefix])
            right = float(tables[position][right_prefix])
            union = max(left, right, min(left + right, 0.5))
            both = left + right - union
            cost += mass * binary_entropy(union)
            for left_bit, right_bit, probability in (
                (0, 0, 1.0 - union),
                (1, 0, union - right),
                (0, 1, union - left),
                (1, 1, both),
            ):
                if probability > 1e-18:
                    key = (
                        left_prefix + (left_bit,),
                        right_prefix + (right_bit,),
                    )
                    next_states[key] = (
                        next_states.get(key, 0.0) + mass * probability
                    )
        assert abs(sum(next_states.values()) - 1.0) < 2e-12
        states = next_states
    return cost


def sequential_value(family):
    iid = shapley_iid(family)
    sequential = sum(
        sequential_cost_float(family, order) for order in ORDERS
    ) / len(ORDERS)
    value = (
        (1.0 - ALPHA) * iid + ALPHA * sequential - log2(len(family))
    )
    return value, iid, sequential


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    binary = Path(
        sys.argv[1] if len(sys.argv) > 1
        else "/tmp/shapley_n5_global_coupling_enum"
    )
    if not binary.is_file():
        raise SystemExit("missing enumerator binary; run the documented build")

    masks, enumerator_log = read_enumerator(binary)
    counts = Counter(len(members(mask)) for mask in masks)
    assert len(masks) == len(set(masks)) == EXPECTED_CANONICAL_TOTAL
    assert dict(counts) == EXPECTED_CANONICAL
    for mask in masks:
        family = members(mask)
        assert exact_admissible(family)
        assert canonical(mask) == mask
    for size in (3, 4):
        emitted = {mask for mask in masks if mask.bit_count() == size}
        assert emitted == brute_canonical(size)
    assert sum(EXPECTED_LABELED.values()) == EXPECTED_LABELED_TOTAL
    for size, labeled_count in EXPECTED_LABELED.items():
        expected_line = (
            "m=%d count=%d canonical=%d"
            % (size, labeled_count, EXPECTED_CANONICAL[size])
        )
        assert expected_line in enumerator_log
    assert "total=463343 canonical_total=5172" in enumerator_log

    print("COMPLETE EXACT ENUMERATION")
    print(enumerator_log.strip())
    print("  parsed canonical representatives:", len(masks))
    print("  Python exact constraint/canonical checks: PASS")
    print("  independent brute canonical checks m=3,4: PASS")
    print()

    records = []
    t0 = time.time()
    for index, mask in enumerate(masks, 1):
        family = members(mask)
        value, iid, sequential = sequential_value(family)
        records.append((value, family, iid, sequential))
        if index % 1000 == 0:
            print("  evaluated %d / %d" % (index, len(masks)))

    worst = min(records)
    assert worst[1] == (0, 3, 4, 8, 15)
    assert abs(worst[0] - 0.0049410954449165) < 2e-12
    assert all(record[0] > 1e-10 for record in records)
    assert not any(union_closed(record[1]) for record in records)

    print()
    print("COMPLETE NUMERICAL GLOBAL-COUPLING EVALUATION")
    print("  runtime: %.2f s" % (time.time() - t0))
    print("  canonical orbits evaluated:", len(records))
    print("  represented labeled families:", EXPECTED_LABELED_TOTAL)
    print("  negative / zero / positive = 0 / 0 /", len(records))
    print("  global minimum")
    print("    family =", worst[1])
    print("    Q / Cseq = %.15f / %.15f" % (worst[2], worst[3]))
    print("    Aseq = %.15f" % worst[0])
    print("  minimum by family size")
    for size in sorted(SPECS):
        size_worst = min(record for record in records
                         if len(record[1]) == size)
        print("    m=%2d orbits=%4d min=%+.12f family=%s" %
              (size, EXPECTED_CANONICAL[size], size_worst[0], size_worst[1]))
    print()
    print("FINITE VERDICT")
    print("  Every cap-2/5 + Reimer family on [5] has positive globally")
    print("  consistent sequential value in the complete float64 comparison.")
    print("  The numerical minimum is a duplicate-plus-dummy lift of the")
    print("  separately certified three-coordinate sign-crossing family.")
    print("  This is complete numerical n=5 evidence, not a theorem.")
