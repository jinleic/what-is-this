"""Complete block-symmetric n=6 cap-2/5 relaxation census.

For k=1,2,3, fix a partition [6]=A union B with |A|=k and enumerate every
family invariant under S_k x S_{6-k}.  Such a family is a union of the
(k+1)(7-k)=12,15,16 weight cells indexed by (|R intersect A|,|R intersect B|),
so the three raw search spaces contain 2^12-1, 2^15-1, and 2^16-1 masks.

Every retained family is checked exactly for cap and Reimer incidence,
canonicalized under the full S_6 action, and evaluated on exact order-orbit
representatives under its full coordinate automorphism group.  Entropy and
coupling values are float64.

Evidence label: COMPLETE EXACT for the declared block-symmetric classes and
S_6/order orbit coverage; COMPLETE NUMERICAL for A_seq and A_plus.  This is not
an arbitrary n=6 census.
"""

from collections import Counter
from math import ceil, log2
import time

from shapley_entropy import ALPHA
from shapley_join_loss import shapley_iid
from shapley_n6_shared_bellman import one_sided_costs, sequential_costs
from shapley_n6_small_complete import (
    DIMENSION,
    members,
    orbit_data,
    representative_orders,
)


FEASIBLE_SPECS = {
    size: (
        2 * size // 5,
        ceil(size * log2(size) / 2.0 - 1e-12),
    )
    for size in range(3, 28)
    if ceil(size * log2(size) / 2.0 - 1e-12)
       <= DIMENSION * (2 * size // 5)
}


def exact_admissible(family):
    if family == (0,) or len(family) not in FEASIBLE_SPECS:
        return False
    cap, need = FEASIBLE_SPECS[len(family)]
    counts = tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    return max(counts) <= cap and sum(counts) >= need


def normalized(family):
    columns = tuple(
        sum(((row >> coordinate) & 1) << index
            for index, row in enumerate(family))
        for coordinate in range(DIMENSION)
    )
    return all(columns) and len(set(columns)) == DIMENSION


def cells(k):
    lower_mask = (1 << k) - 1
    return tuple(
        tuple(
            row for row in range(64)
            if (row & lower_mask).bit_count() == left_weight
            and (row >> k).bit_count() == right_weight
        )
        for left_weight in range(k + 1)
        for right_weight in range(DIMENSION - k + 1)
    )


def enumerate_class(k):
    class_cells = cells(k)
    accepted = []
    for cell_mask in range(1, 1 << len(class_cells)):
        family = tuple(sorted(
            row
            for index, cell in enumerate(class_cells)
            if (cell_mask >> index) & 1
            for row in cell
        ))
        if exact_admissible(family):
            accepted.append(family)
    return tuple(accepted), len(class_cells)


def values(family, orders):
    iid = shapley_iid(family, DIMENSION)
    sequential_roots, sequential_states = sequential_costs(family, orders)
    one_sided_roots, one_sided_states = one_sided_costs(family, orders)
    sequential = sum(sequential_roots) / len(sequential_roots)
    one_sided = sum(one_sided_roots) / len(one_sided_roots)
    entropy = log2(len(family))
    return (
        (1.0 - ALPHA) * iid + ALPHA * sequential - entropy,
        (1.0 - ALPHA) * iid + ALPHA * one_sided - entropy,
        one_sided - sequential,
        max(sequential_states, one_sided_states),
    )


def closure_defect(family):
    support = set(family)
    return sum(
        (left | right) not in support
        for left in family
        for right in family
    ) / (len(family) ** 2)


if __name__ == "__main__":
    print(__doc__.strip().splitlines()[0])
    print()
    started = time.time()
    accepted_by_k = {}
    canonical = {}
    for k in (1, 2, 3):
        accepted, cell_count = enumerate_class(k)
        accepted_by_k[k] = (cell_count, len(accepted))
        for family in accepted:
            canonical_mask, _orbit_size, _automorphisms = orbit_data(family)
            canonical[canonical_mask] = members(canonical_mask)

    records = []
    order_representative_counts = Counter()
    for index, family in enumerate(canonical.values(), 1):
        canonical_mask, _orbit_size, automorphisms = orbit_data(family)
        assert members(canonical_mask) == family
        orders = representative_orders(automorphisms)
        order_representative_counts[len(orders)] += 1
        metrics = values(family, orders)
        records.append((metrics[1], family, metrics))
        if index % 100 == 0:
            print("  evaluated %d / %d" % (index, len(canonical)))

    assert accepted_by_k == {1: (12, 33), 2: (15, 212), 3: (16, 391)}
    worst = min(records)
    normalized_records = [record for record in records if normalized(record[1])]
    normalized_worst = min(normalized_records) if normalized_records else None
    negative_records = sorted(
        record for record in records if record[0] < -1e-10)
    zero_records = sorted(
        record for record in records if abs(record[0]) <= 1e-10)
    corrected_records = sorted(
        (record[0] + 0.02 * closure_defect(record[1]), record)
        for record in records
    )
    corrected_worst = corrected_records[0]
    assert len(records) == 450
    assert len(negative_records) == 2
    assert worst[1] == (
        0, 1, 2, 4, 8, 11, 13, 14, 16, 19, 21, 22, 25,
        26, 28, 32, 35, 37, 38, 41, 42, 44, 49, 50, 52,
    )
    assert abs(worst[0] + 0.013672107732176642) < 3e-12
    assert corrected_worst[1][1] == worst[1]
    assert abs(corrected_worst[0] - 0.0005358922678233586) < 3e-12
    assert corrected_worst[0] > 0

    print()
    print("COMPLETE BLOCK-SYMMETRIC CENSUS")
    print("  raw accepted by k:", accepted_by_k)
    print("  distinct canonical S6 orbits:", len(records))
    print("  order representatives per family:",
          dict(sorted(order_representative_counts.items())))
    print("  negative / zero / positive = %d / %d / %d" % (
        len(negative_records),
        len(zero_records),
        len(records) - len(negative_records) - len(zero_records),
    ))
    print("  global minimum")
    print("    family =", worst[1])
    print("    size =", len(worst[1]))
    print("    A_seq / A_plus / gain = %.15f / %.15f / %.15f" %
          worst[2][:3])
    if negative_records:
        print("  negative families")
        for record in negative_records:
            print("    A_plus=%+.15f m=%d family=%s" %
                  (record[0], len(record[1]), record[1]))
    print("  minimum by family size")
    for size in sorted({len(record[1]) for record in records}):
        size_worst = min(
            record for record in records if len(record[1]) == size)
        print("    m=%d orbits=%d min=%+.12f family=%s" % (
            size,
            sum(len(record[1]) == size for record in records),
            size_worst[0],
            size_worst[1],
        ))
    print("  closure-corrected target")
    print("    min A_plus + (1/50) Pr[X union Y notin F] = %.15f" %
          corrected_worst[0])
    print("    closure defect at corrected minimum = %.15f" %
          closure_defect(corrected_worst[1][1]))
    if normalized_worst:
        print("  normalized minimum")
        print("    family =", normalized_worst[1])
        print("    size =", len(normalized_worst[1]))
        print("    A_seq / A_plus / gain = %.15f / %.15f / %.15f" %
              normalized_worst[2][:3])
    print("  runtime: %.2f s" % (time.time() - started))
