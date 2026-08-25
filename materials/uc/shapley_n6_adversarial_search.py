"""Seeded adversarial search for an n=6 cap-2/5 A_plus counterexample.

This is the all-size heuristic complement to the exact size-3-through-7 census.
For each aggregate-feasible family size 8 through 25, it constructs distinct
row families from random bounded column-degree sequences, screens them on a
fixed eight-order sample, hill-climbs the best candidates, rescreens on 48
orders, and evaluates the finalists on all 6!=720 orders.

Every emitted family is checked exactly for cardinality, coordinate cap,
incidence threshold, activity, and distinct rows.  Entropy and Bellman values
are float64.  The search is deterministic at the recorded seed, but its
coverage is SEEDED HEURISTIC: failure to find a counterexample proves nothing.

Run with one numerical thread and low process priority:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./.venv/bin/python \
        uc/shapley_n6_adversarial_search.py \
        --output ../data/uc-n6-adversarial-search-2026-08-25.json
"""

import argparse
from dataclasses import dataclass
from itertools import permutations
import json
from math import ceil, log2
from pathlib import Path
from random import Random
import time

from shapley_adaptive_coupling import replay_policy
from shapley_entropy import ALPHA
from shapley_join_loss import shapley_iid
from shapley_n6_shared_bellman import one_sided_costs, sequential_costs


DIMENSION = 6
ORDERS = tuple(permutations(range(DIMENSION)))
SEED = 20260825
FEASIBLE_SPECS = {
    size: (
        2 * size // 5,
        ceil(size * log2(size) / 2.0 - 1e-12),
    )
    for size in range(3, 28)
    if ceil(size * log2(size) / 2.0 - 1e-12)
       <= DIMENSION * (2 * size // 5)
}
SEARCH_SIZES = tuple(size for size in FEASIBLE_SPECS if size >= 8)


@dataclass(frozen=True)
class Candidate:
    family: tuple[int, ...]
    iid: float
    sampled_cost: float
    sampled_value: float


def family_counts(family):
    return tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )


def exact_admissible(family, require_active=True):
    if len(family) not in FEASIBLE_SPECS:
        return False
    if len(set(family)) != len(family):
        return False
    cap, incidence_need = FEASIBLE_SPECS[len(family)]
    counts = family_counts(family)
    if require_active and not all(counts):
        return False
    return max(counts) <= cap and sum(counts) >= incidence_need


def normalized(family):
    columns = tuple(
        sum(((row >> coordinate) & 1) << index
            for index, row in enumerate(family))
        for coordinate in range(DIMENSION)
    )
    return all(columns) and len(set(columns)) == DIMENSION


def random_family(size, rng):
    cap, incidence_need = FEASIBLE_SPECS[size]
    for _attempt in range(2000):
        maximum_incidence = DIMENSION * cap
        if rng.random() < 0.8:
            incidence = incidence_need + rng.randrange(
                min(3, maximum_incidence - incidence_need) + 1)
        else:
            incidence = rng.randint(incidence_need, maximum_incidence)

        degrees = [1] * DIMENSION
        remaining = incidence - DIMENSION
        while remaining:
            choices = [
                coordinate for coordinate, degree in enumerate(degrees)
                if degree < cap
            ]
            if not choices:
                break
            degrees[rng.choice(choices)] += 1
            remaining -= 1
        if remaining:
            continue

        rows = [0] * size
        for coordinate, degree in enumerate(degrees):
            for index in rng.sample(range(size), degree):
                rows[index] |= 1 << coordinate
        family = tuple(sorted(rows))
        if exact_admissible(family):
            return family
    return None


def sample_value(family, orders):
    iid = shapley_iid(family, DIMENSION)
    costs, _state_count = one_sided_costs(family, orders)
    one_sided = sum(costs) / len(costs)
    value = (
        (1.0 - ALPHA) * iid
        + ALPHA * one_sided
        - log2(len(family))
    )
    return Candidate(family, iid, one_sided, value)


def mutate(candidate, orders, rng):
    family = list(candidate.family)
    remove_index = rng.randrange(len(family))
    existing = set(family)
    existing.remove(family[remove_index])
    replacement = rng.randrange(64)
    if replacement in existing:
        return candidate
    family[remove_index] = replacement
    proposed = tuple(sorted(family))
    if not exact_admissible(proposed):
        return candidate
    evaluated = sample_value(proposed, orders)
    return evaluated if evaluated.sampled_value < candidate.sampled_value else candidate


def full_value(candidate):
    costs, _state_count = one_sided_costs(candidate.family, ORDERS)
    one_sided = sum(costs) / len(costs)
    value = (
        (1.0 - ALPHA) * candidate.iid
        + ALPHA * one_sided
        - log2(len(candidate.family))
    )
    return value, one_sided


def record(family, value, iid, one_sided):
    counts = family_counts(family)
    cap, need = FEASIBLE_SPECS[len(family)]
    return {
        "family": list(family),
        "size": len(family),
        "coordinate_counts": list(counts),
        "coordinate_cap": cap,
        "incidence": sum(counts),
        "incidence_need": need,
        "normalized": normalized(family),
        "iid_Q": iid,
        "one_sided_C": one_sided,
        "A_plus": value,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--samples-per-size", type=int, default=100)
    parser.add_argument("--hill-starts", type=int, default=4)
    parser.add_argument("--hill-steps", type=int, default=150)
    parser.add_argument("--medium-finalists", type=int, default=64)
    parser.add_argument("--full-finalists", type=int, default=16)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    started = time.time()
    rng = Random(args.seed)
    order_rng = Random(args.seed ^ 0xA51CE)
    screening_orders = tuple(order_rng.sample(list(ORDERS), 8))
    medium_orders = tuple(order_rng.sample(list(ORDERS), 48))

    candidates = []
    generated_by_size = {}
    for size in SEARCH_SIZES:
        seen = set()
        attempts = 0
        while len(seen) < args.samples_per_size and attempts < 10000:
            attempts += 1
            family = random_family(size, rng)
            if family is None or family in seen:
                continue
            seen.add(family)
            candidates.append(sample_value(family, screening_orders))
        generated_by_size[size] = len(seen)
        print("  size %d generated %d" % (size, len(seen)))

    candidates.sort(key=lambda candidate: candidate.sampled_value)
    hill_candidates = []
    for start in candidates[:args.hill_starts]:
        current = start
        for _step in range(args.hill_steps):
            current = mutate(current, screening_orders, rng)
        hill_candidates.append(current)
        print("  hill size %d A_sample=%+.12f" %
              (len(current.family), current.sampled_value))

    unique = {
        candidate.family: candidate
        for candidate in candidates + hill_candidates
    }
    medium_pool = sorted(
        unique.values(), key=lambda candidate: candidate.sampled_value
    )[:args.medium_finalists]
    medium = [sample_value(candidate.family, medium_orders)
              for candidate in medium_pool]
    medium.sort(key=lambda candidate: candidate.sampled_value)

    full = []
    for index, candidate in enumerate(medium[:args.full_finalists], 1):
        value, one_sided = full_value(candidate)
        full.append((value, candidate, one_sided))
        print("  full %d/%d size=%d A_plus=%+.12f" %
              (index, min(args.full_finalists, len(medium)),
               len(candidate.family), value))
    full.sort(key=lambda item: item[0])
    best_value, best, best_cost = full[0]

    sequential_roots, _state_count = sequential_costs(best.family, ORDERS)
    sequential_cost = sum(sequential_roots) / len(sequential_roots)
    sequential_value = (
        (1.0 - ALPHA) * best.iid
        + ALPHA * sequential_cost
        - log2(len(best.family))
    )
    for order in ORDERS:
        replay_policy(best.family, order)

    result = {
        "schema": 1,
        "evidence_label": "SEEDED_HEURISTIC_ALL_SIZE_N6_SEARCH",
        "seed": args.seed,
        "config": {
            "samples_per_size": args.samples_per_size,
            "hill_starts": args.hill_starts,
            "hill_steps": args.hill_steps,
            "medium_finalists": args.medium_finalists,
            "full_finalists": args.full_finalists,
            "screening_orders": [list(order) for order in screening_orders],
            "medium_orders": [list(order) for order in medium_orders],
        },
        "feasible_specs": {
            str(size): {"cap": cap, "incidence_need": need}
            for size, (cap, need) in FEASIBLE_SPECS.items()
        },
        "searched_sizes": list(SEARCH_SIZES),
        "generated_by_size": {
            str(size): count for size, count in generated_by_size.items()
        },
        "total_generated": sum(generated_by_size.values()),
        "best": record(best.family, best_value, best.iid, best_cost),
        "best_A_seq": sequential_value,
        "best_C_seq": sequential_cost,
        "full_finalists": [
            record(candidate.family, value, candidate.iid, one_sided)
            for value, candidate, one_sided in full
        ],
        "runtime_seconds": time.time() - started,
        "verdict": (
            "COUNTEREXAMPLE_FOUND" if best_value < -1e-10
            else "NO_COUNTEREXAMPLE_IN_SEEDED_SEARCH"
        ),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print()
    print("SEEDED ALL-SIZE VERDICT")
    print("  generated:", result["total_generated"])
    print("  best family =", tuple(result["best"]["family"]))
    print("  best size =", result["best"]["size"])
    print("  best counts =", tuple(result["best"]["coordinate_counts"]))
    print("  best normalized =", result["best"]["normalized"])
    print("  A_seq / A_plus = %.15f / %.15f" %
          (sequential_value, best_value))
    print("  verdict =", result["verdict"])
    print("  runtime = %.2f s" % result["runtime_seconds"])
