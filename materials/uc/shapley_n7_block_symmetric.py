"""Complete block-symmetric n=7 closure-defect census.

For k=1,2,3, fix a partition [7]=A union B with |A|=k and enumerate every
family invariant under S_k x S_{7-k}. Such a family is a union of the
(k+1)(8-k)=14,18,20 weight cells indexed by
(|R intersect A|, |R intersect B|), so the three raw search spaces contain
2^14-1, 2^18-1, and 2^20-1 masks.

Every retained family is checked exactly for the cap-2/5 and Reimer-incidence
conditions, canonicalized under the full S_7 action, and evaluated on exact
order-orbit representatives under its full coordinate automorphism group.
The closure defect is exact. Entropy and coupling values are float64.

Evidence label: COMPLETE EXACT for the declared block-symmetric classes,
S_7/order-orbit coverage, and closure defects; COMPLETE NUMERICAL for A_seq,
A_plus, and the closure-corrected comparison. This is not an arbitrary n=7
census or a proof of the union-closed conjecture.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import log2
import time

from shapley_adaptive_coupling import bellman_coupling
from shapley_entropy import ALPHA
from shapley_global_coupling import sequential_coupling
from shapley_join_loss import shapley_iid
from shapley_n6_shared_bellman import one_sided_costs, sequential_costs


DIMENSION = 7
SUBSET_COUNT = 1 << DIMENSION
ORDERS = tuple(permutations(range(DIMENSION)))
def exact_reimer_threshold(size):
    """Return ceil(size*log2(size)/2) by exact integer comparison."""
    target = size ** size
    low = 0
    high = size * (size - 1).bit_length()
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


FEASIBLE_SPECS = {
    size: (2 * size // 5, exact_reimer_threshold(size))
    for size in range(3, SUBSET_COUNT + 1)
    if exact_reimer_threshold(size) <= DIMENSION * (2 * size // 5)
}


def members(mask):
    return tuple(row for row in range(SUBSET_COUNT) if (mask >> row) & 1)


def family_mask(family):
    return sum(1 << row for row in family)


def exact_admissible(family):
    if len(family) not in FEASIBLE_SPECS:
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


def permute_subset(subset, permutation):
    image = 0
    for old_coordinate, new_coordinate in enumerate(permutation):
        if (subset >> old_coordinate) & 1:
            image |= 1 << new_coordinate
    return image


PERMUTED_SUBSETS = tuple(
    tuple(permute_subset(subset, permutation)
          for subset in range(SUBSET_COUNT))
    for permutation in ORDERS
)


def permute_family_mask(family, permutation_index):
    images = PERMUTED_SUBSETS[permutation_index]
    return sum(1 << images[row] for row in family)


def orbit_data(family):
    mask = family_mask(family)
    images = {
        permute_family_mask(family, permutation_index)
        for permutation_index in range(len(ORDERS))
    }
    automorphisms = tuple(
        ORDERS[index]
        for index in range(len(ORDERS))
        if permute_family_mask(family, index) == mask
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


def cells(k):
    lower_mask = (1 << k) - 1
    return tuple(
        tuple(
            row for row in range(SUBSET_COUNT)
            if (row & lower_mask).bit_count() == left_weight
            and (row >> k).bit_count() == right_weight
        )
        for left_weight in range(k + 1)
        for right_weight in range(DIMENSION - k + 1)
    )


def partition_representatives(k):
    """One coordinate permutation for each image of the k-coordinate block."""
    result = []
    coordinates = tuple(range(DIMENSION))
    for target_left in combinations(coordinates, k):
        target_left = tuple(target_left)
        target_right = tuple(
            coordinate for coordinate in coordinates
            if coordinate not in target_left
        )
        result.append(target_left + target_right)
    return tuple(result)


def canonical_block_family(family, k):
    """Canonical S_7 image, quotienting known S_k x S_(7-k) symmetry."""
    masks = []
    for permutation in partition_representatives(k):
        images = tuple(permute_subset(row, permutation) for row in family)
        masks.append(family_mask(images))
    return min(masks)


def enumerate_class(k):
    class_cells = cells(k)
    cell_masks = tuple(family_mask(cell) for cell in class_cells)
    cell_stats = tuple(
        (
            len(cell),
            sum((row >> 0) & 1 for row in cell),
            sum((row >> k) & 1 for row in cell),
            sum(row.bit_count() for row in cell),
        )
        for cell in class_cells
    )

    accepted = []
    current_family = 0
    size = left_count = right_count = incidence = 0
    previous_gray = 0
    for step in range(1, 1 << len(class_cells)):
        gray = step ^ (step >> 1)
        changed = gray ^ previous_gray
        index = changed.bit_length() - 1
        direction = 1 if gray & changed else -1
        cell_size, left, right, total = cell_stats[index]
        current_family ^= cell_masks[index]
        size += direction * cell_size
        left_count += direction * left
        right_count += direction * right
        incidence += direction * total
        previous_gray = gray

        spec = FEASIBLE_SPECS.get(size)
        if spec is None:
            continue
        cap, need = spec
        if max(left_count, right_count) > cap or incidence < need:
            continue
        family = members(current_family)
        assert len(family) == size
        assert exact_admissible(family)
        accepted.append(family)

    return tuple(accepted), len(class_cells)


def values(family, orders):
    iid = shapley_iid(family, DIMENSION)
    sequential_roots, sequential_states = sequential_costs(
        family, orders, dimension=DIMENSION)
    one_sided_roots, one_sided_states = one_sided_costs(
        family, orders, dimension=DIMENSION)
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
    missing = sum(
        (left | right) not in support
        for left in family
        for right in family
    )
    return Fraction(missing, len(family) ** 2)


def independent_replay(family, orders):
    shared_sequential, _ = sequential_costs(
        family, orders, dimension=DIMENSION)
    shared_one_sided, _ = one_sided_costs(
        family, orders, dimension=DIMENSION)
    direct_sequential = tuple(
        sequential_coupling(family, order)[0] for order in orders)
    direct_one_sided = tuple(
        bellman_coupling(family, order) for order in orders)
    assert max(abs(left - right) for left, right in zip(
        shared_sequential, direct_sequential)) < 3e-12
    assert max(abs(left - right) for left, right in zip(
        shared_one_sided, direct_one_sided)) < 3e-12


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
            canonical_mask = canonical_block_family(family, k)
            canonical[canonical_mask] = members(canonical_mask)

    assert accepted_by_k == {1: (14, 72), 2: (18, 680), 3: (20, 2238)}

    records = []
    order_representative_counts = Counter()
    for index, family in enumerate(canonical.values(), 1):
        canonical_mask, _orbit_size, automorphisms = orbit_data(family)
        assert members(canonical_mask) == family
        orders = representative_orders(automorphisms)
        order_representative_counts[len(orders)] += 1
        metrics = values(family, orders)
        defect = closure_defect(family)
        records.append((metrics[1], family, metrics, defect))
        if index % 100 == 0:
            print("  evaluated %d / %d" % (index, len(canonical)), flush=True)

    worst = min(records)
    corrected_records = sorted(
        (record[0] + float(record[3]) / 50.0, record)
        for record in records
    )
    corrected_worst = corrected_records[0]
    negative_records = sorted(
        record for record in records if record[0] < -1e-10)
    corrected_negative = sorted(
        record for record in corrected_records if record[0] < -1e-10)
    assert all(record[3] > 0 for record in negative_records)
    required_coefficient, required_family = max(
        (-record[0] / float(record[3]), record)
        for record in negative_records
    )
    one_twenty_fifth_records = sorted(
        (record[0] + float(record[3]) / 25.0, record)
        for record in records
    )
    normalized_records = [record for record in records if normalized(record[1])]
    normalized_worst = min(normalized_records) if normalized_records else None

    expected_worst = (
        1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 21, 22, 25, 26, 28,
        32, 35, 37, 38, 41, 42, 44, 49, 50, 52, 56, 63, 64, 67,
        69, 70, 73, 74, 76, 81, 82, 84, 88, 95, 96, 97, 98, 100,
        104, 112,
    )
    assert len(records) == 2980
    assert order_representative_counts == {
        1: 5, 7: 67, 21: 675, 35: 2233,
    }
    assert len(negative_records) == 16
    assert len(corrected_negative) == 4
    assert worst[1] == expected_worst
    assert worst[3] == Fraction(64, 81)
    assert abs(worst[0] + 0.028649186794459158) < 3e-12
    assert corrected_worst[1][1] == expected_worst
    assert abs(corrected_worst[0] + 0.012846717658656072) < 3e-12

    worst_mask, _orbit_size, worst_automorphisms = orbit_data(worst[1])
    assert members(worst_mask) == worst[1]
    worst_orders = representative_orders(worst_automorphisms)
    independent_replay(worst[1], worst_orders)

    print()
    print("COMPLETE BLOCK-SYMMETRIC n=7 CENSUS")
    print("  raw accepted by k:", accepted_by_k)
    print("  distinct canonical S7 orbits:", len(records))
    print("  order representatives per family:",
          dict(sorted(order_representative_counts.items())))
    print("  A_plus negative / nonnegative = %d / %d" % (
        len(negative_records), len(records) - len(negative_records)))
    print("  global A_plus minimum")
    print("    family =", worst[1])
    print("    size =", len(worst[1]))
    print("    closure defect =", worst[3])
    print("    A_seq / A_plus / gain = %.15f / %.15f / %.15f" %
          worst[2][:3])
    print("  closure-corrected target")
    print("    minimum A_plus + (1/50) Pr[X union Y notin F] = %.15f" %
          corrected_worst[0])
    print("    family =", corrected_worst[1][1])
    print("    size =", len(corrected_worst[1][1]))
    print("    closure defect =", corrected_worst[1][3])
    print("    corrected negative =", len(corrected_negative))
    print("    largest required scalar coefficient = %.15f" %
          required_coefficient)
    print("    coefficient extremizer =", required_family[1])
    print("    minimum with coefficient 1/25 = %.15f" %
          one_twenty_fifth_records[0][0])
    if normalized_worst:
        print("  normalized A_plus minimum")
        print("    family =", normalized_worst[1])
        print("    size =", len(normalized_worst[1]))
        print("    A_seq / A_plus / gain = %.15f / %.15f / %.15f" %
              normalized_worst[2][:3])
    print("  independent worst-family Bellman replay: PASS")
    print("  runtime: %.2f s" % (time.time() - started))
