#!/usr/bin/env python3
"""Resumable adversarial search for a counterexample to Gate B tensorization.

The Gate B theorem rests on exact additivity of the Shapley iid term and of the
one-sided Bellman term under disjoint-block Cartesian products, together with
multiplicativity of closure success.  Additivity is proved in ``PROOF.md``; this
script attacks it numerically on unequal-dimension, asymmetric, degenerate, and
interleaved instances that the earlier equal-dimension tests never touched.

For every sampled ordered pair of families and every global order of the
product, it compares

    Q_pi(F x G)     against Q_{pi_F}(F) + Q_{pi_G}(G),
    C_{+,pi}(F x G) against C_{+,pi_F}(F) + C_{+,pi_G}(G),

and checks the exact closure-defect product identity.  Records append to a
JSONL checkpoint, so an interrupted run resumes without repeating work.

Run from ``math/``:

    OMP_NUM_THREADS=1 nice -n 10 ./.venv/bin/python -B \\
        uc/gate_b/audit_tensorization.py [--limit N]
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
import json
from math import comb, log2
from pathlib import Path
from random import Random
import sys
import time
from typing import Any, Iterator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from shapley_global_coupling import binary_entropy  # noqa: E402
from shapley_n6_shared_bellman import one_sided_costs  # noqa: E402


CHECKPOINT = HERE / "experiments" / "tensorization_audit_checkpoint.jsonl"
SUMMARY = HERE / "candidates" / "tensorization_audit.json"
TOLERANCE = 1e-9
SEED = 20260826


def defect(family: tuple[int, ...]) -> Fraction:
    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    return Fraction(missing, len(family) ** 2)


def shapley_iid(family: tuple[int, ...], dimension: int) -> float:
    """Order-averaged iid OR-entropy, float64, computed independently here."""
    total = 0.0
    for coordinate in range(dimension):
        others = tuple(index for index in range(dimension) if index != coordinate)
        for count in range(dimension):
            weight = 1.0 / (dimension * comb(dimension - 1, count))
            for predecessors in _subsets(others, count):
                total += weight * _local_iid(family, coordinate, predecessors)
    return total


def _subsets(items: tuple[int, ...], size: int) -> Iterator[tuple[int, ...]]:
    if size == 0:
        yield ()
        return
    for index in range(len(items) - size + 1):
        for tail in _subsets(items[index + 1 :], size - 1):
            yield (items[index],) + tail


def _local_iid(
    family: tuple[int, ...], coordinate: int, predecessors: tuple[int, ...]
) -> float:
    groups: dict[int, list[int]] = {}
    for row in family:
        prefix = 0
        for index, predecessor in enumerate(predecessors):
            prefix |= ((row >> predecessor) & 1) << index
        entry = groups.get(prefix)
        if entry is None:
            groups[prefix] = [1, (row >> coordinate) & 1]
        else:
            entry[0] += 1
            entry[1] += (row >> coordinate) & 1
    size = len(family)
    law = [(ones / total, total / size) for total, ones in groups.values()]
    return sum(
        left_mass
        * right_mass
        * binary_entropy(left + right - left * right)
        for left, left_mass in law
        for right, right_mass in law
    )


def fixed_order_iid(
    family: tuple[int, ...], order: tuple[int, ...]
) -> float:
    """Fixed-order iid term, i.e. Q_pi rather than its order average."""
    return sum(
        _local_iid(family, coordinate, tuple(order[:position]))
        for position, coordinate in enumerate(order)
    )


@lru_cache(maxsize=None)
def fixed_order_bellman(
    family: tuple[int, ...], order: tuple[int, ...], dimension: int
) -> float:
    costs, _states = one_sided_costs(family, (order,), dimension=dimension)
    return costs[0]


def cartesian(
    left: tuple[int, ...], left_dimension: int, right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        sorted(
            left_row | (right_row << left_dimension)
            for left_row in left
            for right_row in right
        )
    )


def induced_orders(
    order: tuple[int, ...], left_dimension: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left_order = tuple(
        coordinate for coordinate in order if coordinate < left_dimension
    )
    right_order = tuple(
        coordinate - left_dimension
        for coordinate in order
        if coordinate >= left_dimension
    )
    return left_order, right_order


def _named_families() -> list[tuple[str, tuple[int, ...], int]]:
    """Structurally diverse probes, including degenerate and asymmetric ones."""
    return [
        ("singleton_empty", (0,), 1),
        ("singleton_full", (1,), 1),
        ("d1_both", (0, 1), 1),
        ("d2_union_closed", (0, 1, 2, 3), 2),
        ("d2_missing_join", (0, 1, 2), 2),
        ("d2_antichain", (1, 2), 2),
        ("d2_dead_column", (0, 1), 2),
        ("d3_odd_size", (0, 1, 2, 4, 7), 3),
        ("d3_no_empty", (1, 2, 4, 7), 3),
        ("d3_dead_and_duplicate", (0, 3), 3),
        ("d3_full_power_set", tuple(range(8)), 3),
        ("d3_middle_layer", (3, 5, 6), 3),
        ("d4_sparse", (0, 1, 6, 9, 15), 4),
        ("d4_biased", (1, 3, 7, 15), 4),
    ]


def _random_family(random: Random, dimension: int) -> tuple[int, ...]:
    universe = 1 << dimension
    size = random.randint(2, min(universe, 9))
    return tuple(sorted(random.sample(range(universe), size)))


def _cases(limit: int | None) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    named = _named_families()
    for left_name, left, left_dimension in named:
        for right_name, right, right_dimension in named:
            if left_dimension + right_dimension > 6:
                continue
            cases.append(
                {
                    "case_id": f"named:{left_name}|{right_name}",
                    "left": list(left),
                    "left_dimension": left_dimension,
                    "right": list(right),
                    "right_dimension": right_dimension,
                }
            )

    random = Random(SEED)
    index = 0
    while index < 120:
        left_dimension = random.randint(1, 3)
        right_dimension = random.randint(1, 6 - left_dimension)
        left = _random_family(random, left_dimension)
        right = _random_family(random, right_dimension)
        cases.append(
            {
                "case_id": f"random:{index}",
                "left": list(left),
                "left_dimension": left_dimension,
                "right": list(right),
                "right_dimension": right_dimension,
            }
        )
        index += 1

    if limit is not None:
        cases = cases[:limit]
    return cases


def _evaluate(case: dict[str, Any]) -> dict[str, Any]:
    left = tuple(case["left"])
    right = tuple(case["right"])
    left_dimension = case["left_dimension"]
    right_dimension = case["right_dimension"]
    product = cartesian(left, left_dimension, right)
    dimension = left_dimension + right_dimension

    worst_iid = 0.0
    worst_bellman = 0.0
    for order in permutations(range(dimension)):
        left_order, right_order = induced_orders(order, left_dimension)
        iid_gap = abs(
            fixed_order_iid(product, order)
            - fixed_order_iid(left, left_order)
            - fixed_order_iid(right, right_order)
        )
        bellman_gap = abs(
            fixed_order_bellman(product, order, dimension)
            - fixed_order_bellman(left, left_order, left_dimension)
            - fixed_order_bellman(right, right_order, right_dimension)
        )
        worst_iid = max(worst_iid, iid_gap)
        worst_bellman = max(worst_bellman, bellman_gap)

    averaged_iid_gap = abs(
        shapley_iid(product, dimension)
        - shapley_iid(left, left_dimension)
        - shapley_iid(right, right_dimension)
    )
    product_defect = defect(product)
    expected_defect = 1 - (1 - defect(left)) * (1 - defect(right))
    entropy_gap = abs(
        log2(len(product)) - log2(len(left)) - log2(len(right))
    )

    return {
        "case_id": case["case_id"],
        "left": list(left),
        "left_dimension": left_dimension,
        "right": list(right),
        "right_dimension": right_dimension,
        "product_size": len(product),
        "worst_fixed_order_iid_gap": worst_iid,
        "worst_fixed_order_bellman_gap": worst_bellman,
        "order_averaged_iid_gap": averaged_iid_gap,
        "log_size_gap": entropy_gap,
        "defect_identity_exact": product_defect == expected_defect,
        "product_defect": str(product_defect),
        "expected_defect": str(expected_defect),
        "counterexample": (
            worst_iid > TOLERANCE
            or worst_bellman > TOLERANCE
            or averaged_iid_gap > TOLERANCE
            or entropy_gap > TOLERANCE
            or product_defect != expected_defect
        ),
    }


def _load_checkpoint() -> dict[str, dict[str, Any]]:
    if not CHECKPOINT.exists():
        return {}
    records = {}
    for line in CHECKPOINT.read_text().splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        records[record["case_id"]] = record
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    cases = _cases(args.limit)
    done = _load_checkpoint()
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    started = time.time()
    new_records = 0
    with CHECKPOINT.open("a") as handle:
        for case in cases:
            if case["case_id"] in done:
                continue
            record = _evaluate(case)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            done[case["case_id"]] = record
            new_records += 1

    evaluated = [done[case["case_id"]] for case in cases]
    counterexamples = [
        record for record in evaluated if record["counterexample"]
    ]
    worst_iid = max(
        (record["worst_fixed_order_iid_gap"] for record in evaluated),
        default=0.0,
    )
    worst_bellman = max(
        (record["worst_fixed_order_bellman_gap"] for record in evaluated),
        default=0.0,
    )
    worst_average = max(
        (record["order_averaged_iid_gap"] for record in evaluated),
        default=0.0,
    )
    summary = {
        "schema_version": 1,
        "status": "complete" if not args.limit else "partial",
        "evidence_status": "adversarial_numerical_falsification_attempt",
        "case_count": len(cases),
        "float64_tolerance": TOLERANCE,
        "counterexample_count": len(counterexamples),
        "exact_defect_identity_failures": sum(
            0 if record["defect_identity_exact"] else 1
            for record in evaluated
        ),
        "worst_fixed_order_iid_gap": worst_iid,
        "worst_fixed_order_bellman_gap": worst_bellman,
        "worst_order_averaged_iid_gap": worst_average,
        "unequal_dimension_case_count": sum(
            1
            for record in evaluated
            if record["left_dimension"] != record["right_dimension"]
        ),
        "seed": SEED,
    }
    # Run-dependent counters stay out of the persisted artifact so that a
    # resume of a complete search reproduces it byte-identically.
    report = {
        **summary,
        "new_evaluations": new_records,
        "runtime_seconds": round(time.time() - started, 2),
    }
    if args.limit is None:
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
