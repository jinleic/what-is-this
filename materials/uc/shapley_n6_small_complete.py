"""Complete numerical A_plus census for exact n=6 family sizes 3 through 7.

The C++ companion enumerates every labeled nontrivial simple family on six
coordinates with 3 <= |F| <= 7 satisfying the cap 2/5 and Reimer incidence
threshold, and emits one representative per S_6 coordinate orbit.  This script
independently checks every representative, canonicality, and orbit size; the
sum of orbit sizes must reproduce the C++ labeled counts.

For every orbit it computes the exact all-order average of the one-sided
Bellman functional in float64.  Orders are quotiented only by the family's
actual coordinate automorphism group; the action on coordinate orders is free,
so every retained order represents exactly |Aut(F)| original orders.

Evidence labels:

* COMPLETE EXACT: family filters, canonical S_6 coverage, orbit sizes, and
  order-orbit coverage.
* COMPLETE NUMERICAL: entropy, Bellman costs, and A_plus comparisons.

This is a finite n=6 relaxation census through family size seven, not a
union-closed theorem and not a census of larger admissible family sizes.

Build and run:
    c++ -O3 -std=c++20 uc/shapley_n6_small_enumerate.cpp \
        -o /tmp/shapley_n6_small_enum
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python \
        uc/shapley_n6_small_complete.py /tmp/shapley_n6_small_enum
"""

from collections import Counter
from itertools import permutations
from math import log2
from pathlib import Path
import re
import subprocess
import sys
import time

from shapley_adaptive_coupling import replay_policy
from shapley_entropy import ALPHA, union_closed
from shapley_join_loss import shapley_iid
from shapley_n6_shared_bellman import one_sided_costs, sequential_costs


DIMENSION = 6
ORDERS = tuple(permutations(range(DIMENSION)))
SPECS = {
    3: (1, 3),
    4: (1, 4),
    5: (2, 6),
    6: (2, 8),
    7: (2, 10),
}
EXPECTED_LABELED = {
    3: 636,
    4: 470,
    5: 118812,
    6: 101896,
    7: 40771,
}
EXPECTED_ORBITS = {
    3: 15,
    4: 10,
    5: 492,
    6: 390,
    7: 174,
}


def members(mask):
    return tuple(row for row in range(64) if (mask >> row) & 1)


def family_mask(family):
    return sum(1 << row for row in family)


def exact_admissible(family):
    cap, incidence_need = SPECS[len(family)]
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    return max(counts) <= cap and sum(counts) >= incidence_need


def permute_subset(subset, permutation):
    image = 0
    for old_coordinate, new_coordinate in enumerate(permutation):
        if (subset >> old_coordinate) & 1:
            image |= 1 << new_coordinate
    return image


PERMUTED_SUBSETS = tuple(
    tuple(permute_subset(subset, permutation) for subset in range(64))
    for permutation in ORDERS
)


def permute_family(family, permutation_index):
    images = PERMUTED_SUBSETS[permutation_index]
    return family_mask(tuple(images[row] for row in family))


def orbit_data(family):
    mask = family_mask(family)
    images = {
        permute_family(family, permutation_index)
        for permutation_index in range(len(ORDERS))
    }
    automorphisms = tuple(
        ORDERS[index]
        for index in range(len(ORDERS))
        if permute_family(family, index) == mask
    )
    assert len(images) * len(automorphisms) == len(ORDERS)
    return min(images), len(images), automorphisms


def representative_orders(automorphisms):
    unseen = set(ORDERS)
    representatives = []
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
    assert len(representatives) * len(automorphisms) == len(ORDERS)
    return tuple(representatives)


def normalized(family):
    columns = tuple(
        sum(((row >> coordinate) & 1) << index
            for index, row in enumerate(family))
        for coordinate in range(DIMENSION)
    )
    return all(columns) and len(set(columns)) == DIMENSION


def read_enumerator(binary):
    result = subprocess.run(
        [str(binary)], capture_output=True, text=True, check=True)
    masks = tuple(int(line, 16) for line in result.stdout.splitlines())
    counts = {}
    for size, labeled, orbits in re.findall(
            r"m=(\d+) labeled=(\d+) orbits=(\d+)", result.stderr):
        counts[int(size)] = (int(labeled), int(orbits))
    return masks, counts, result.stderr


def family_value(family, orders):
    iid = shapley_iid(family, DIMENSION)
    costs, _state_count = one_sided_costs(family, orders)
    one_sided = sum(costs) / len(costs)
    value = (
        (1.0 - ALPHA) * iid
        + ALPHA * one_sided
        - log2(len(family))
    )
    return value, iid, one_sided


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    binary = Path(
        sys.argv[1] if len(sys.argv) > 1 else "/tmp/shapley_n6_small_enum"
    )
    if not binary.is_file():
        raise SystemExit("missing exact n=6 enumerator binary")

    masks, enumerator_counts, enumerator_log = read_enumerator(binary)
    sizes = Counter(mask.bit_count() for mask in masks)
    assert len(masks) == len(set(masks))
    assert set(enumerator_counts) == set(SPECS)

    labeled_from_orbits = Counter()
    order_representative_counts = Counter()
    checked = []
    for mask in masks:
        family = members(mask)
        assert exact_admissible(family)
        canonical, orbit_size, automorphisms = orbit_data(family)
        assert canonical == mask
        labeled_from_orbits[len(family)] += orbit_size
        order_representatives = representative_orders(automorphisms)
        order_representative_counts[len(order_representatives)] += 1
        checked.append((family, order_representatives))

    for size, (labeled, orbits) in enumerator_counts.items():
        assert sizes[size] == orbits
        assert labeled_from_orbits[size] == labeled
    for size, expected in EXPECTED_LABELED.items():
        assert enumerator_counts[size][0] == expected
    for size, expected in EXPECTED_ORBITS.items():
        assert enumerator_counts[size][1] == expected

    print("COMPLETE EXACT ENUMERATION AND ORBIT CHECKS")
    print(enumerator_log.strip())
    print("  parsed canonical representatives:", len(masks))
    print("  order representatives per family:",
          dict(sorted(order_representative_counts.items())))
    print()

    started = time.time()
    records = []
    for index, (family, orders) in enumerate(checked, 1):
        value, iid, one_sided = family_value(family, orders)
        records.append((value, family, iid, one_sided, len(orders)))
        if index % 100 == 0:
            print("  evaluated %d / %d" % (index, len(checked)))

    worst = min(records)
    assert all(record[0] > 1e-10 for record in records)
    assert not any(union_closed(record[1]) for record in records)
    for order in ORDERS:
        replay_policy(worst[1], order)
    _canonical, _orbit_size, worst_automorphisms = orbit_data(worst[1])
    worst_orders = representative_orders(worst_automorphisms)
    sequential_roots, _state_count = sequential_costs(
        worst[1], worst_orders)
    sequential_cost = sum(sequential_roots) / len(sequential_roots)
    sequential = (
        (1.0 - ALPHA) * worst[2]
        + ALPHA * sequential_cost
        - log2(len(worst[1]))
    )

    print()
    print("COMPLETE NUMERICAL ONE-SIDED EVALUATION")
    print("  runtime: %.2f s" % (time.time() - started))
    print("  exact nontrivial labeled families:", sum(labeled_from_orbits.values()))
    print("  exact S6 orbits:", len(records))
    print("  negative / zero / positive = 0 / 0 /", len(records))
    print("  global small-size minimum")
    print("    family =", worst[1])
    print("    size =", len(worst[1]))
    print("    normalized =", normalized(worst[1]))
    print("    all-order representatives =", worst[4])
    print("    A_seq / A_plus = %.15f / %.15f" %
          (sequential, worst[0]))
    print("    C_plus-C_seq = %.15f" % (worst[3] - sequential_cost))
    print("  minimum by family size")
    for size in sorted(SPECS):
        size_worst = min(record for record in records
                         if len(record[1]) == size)
        print("    m=%d orbits=%d min=%+.12f family=%s" %
              (size, sizes[size], size_worst[0], size_worst[1]))
    normalized_records = [record for record in records if normalized(record[1])]
    if normalized_records:
        normalized_worst = min(normalized_records)
        print("  normalized minimum = %+.12f at %s" %
              (normalized_worst[0], normalized_worst[1]))
    else:
        print("  normalized families: none")
    print("  worst-policy full 720-order replay: PASS")
