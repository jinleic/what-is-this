#!/usr/bin/env python3
"""Exact rational falsification search for the Gate B product lemmas.

The three product lemmas of PROOF.md are the only human-proved step that the
all-\\(k\\) theorem consumes:

1.  the closure defect multiplies through join success,
2.  the order-averaged iid cost \\(Q\\) is additive under disjoint-block products,
3.  every fixed-order one-sided Bellman optimum \\(C_{+,\\pi}\\) is additive, with
    the two induced relative orders on the blocks.

Until now their machine support was float64: the 292-case search reported a
worst gap of \\(1.8\\times10^{-15}\\), i.e. rounding.  A float64 audit cannot
distinguish "additive" from "additive up to 1e-15", and it cannot see a defect
that only appears in exact arithmetic.

This audit repeats the attack in exact rational arithmetic.  For every global
order of every product it computes a certified rational enclosure of both sides
and demands that the two enclosures intersect.  Non-intersection is a
counterexample to the lemma, not a rounding artifact.  It also records the exact
Hausdorff distance between the two enclosures, which is bounded by the
certificate width rather than by machine epsilon.

The search is resumable at case granularity through an append-only checkpoint.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations
import json
from math import factorial
import os
from pathlib import Path
from random import Random
import sys
import time
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import verify_gate_b_rational as rational

DEFAULT_CHECKPOINT = HERE / "experiments" / "tensorization_exact_checkpoint.jsonl"
DEFAULT_RESULT = HERE / "candidates" / "tensorization_exact_audit.json"
RECORD_SCHEMA = 1
TOLERANCE = Fraction(1, 10**20)


def cartesian_product(
    left: tuple[int, ...], left_dimension: int, right: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        sorted(
            left_row | (right_row << left_dimension)
            for left_row in left
            for right_row in right
        )
    )


def closure_defect(family: tuple[int, ...]) -> Fraction:
    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    return Fraction(missing, len(family) ** 2)


def induced_orders(
    order: tuple[int, ...], left_dimension: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Relative orders that the global order induces on the two blocks."""
    left = tuple(
        coordinate for coordinate in order if coordinate < left_dimension
    )
    right = tuple(
        coordinate - left_dimension
        for coordinate in order
        if coordinate >= left_dimension
    )
    return left, right


def _relabel(order: tuple[int, ...]) -> tuple[int, ...]:
    """Rewrite a relative order as a permutation of its own index range."""
    rank = {coordinate: index for index, coordinate in enumerate(sorted(order))}
    positions = [0] * len(order)
    for position, coordinate in enumerate(order):
        positions[position] = rank[coordinate]
    return tuple(positions)


def _cases(seed: int) -> tuple[dict[str, Any], ...]:
    random = Random(seed)
    cases: list[dict[str, Any]] = []

    def add(
        name: str,
        left: tuple[int, ...],
        left_dimension: int,
        right: tuple[int, ...],
        right_dimension: int,
    ) -> None:
        cases.append(
            {
                "name": name,
                "left": tuple(sorted(set(left))),
                "left_dimension": left_dimension,
                "right": tuple(sorted(set(right))),
                "right_dimension": right_dimension,
            }
        )

    # Declared adversarial pairs: singletons, full power sets, dead columns,
    # duplicate columns, deterministic coordinates, and unequal dimensions.
    add("singleton_x_singleton", (0,), 1, (1,), 1)
    add("bit_x_bit", (0, 1), 1, (0, 1), 1)
    add("bit_x_pair", (0, 1), 1, (0, 3), 2)
    add("dead_column_x_bit", (0, 1), 2, (0, 1), 1)
    add("duplicate_columns_x_bit", (0, 3), 2, (0, 1), 1)
    add("power_set_x_bit", tuple(range(4)), 2, (0, 1), 1)
    add("power_set_x_power_set", tuple(range(4)), 2, tuple(range(4)), 2)
    add("antichain_x_pair", (1, 2), 2, (0, 3), 2)
    add("chain_x_triple", (0, 1, 3), 2, (1, 2, 4), 3)
    add("triple_x_triple", (1, 2, 4), 3, (3, 5, 6), 3)
    add("union_closed_x_not", (0, 1, 3), 2, (1, 2), 2)
    add("full_x_triple", tuple(range(8)), 3, (1, 2, 4), 3)
    add("weight_cells_x_bit", (0, 3, 5, 6), 3, (0, 1), 1)
    add("sparse_x_dense", (7,), 3, (0, 1, 2, 3), 2)

    # Full six-coordinate products, where the global order can interleave the
    # blocks in every one of 720 ways and the induced relative orders differ.
    add("bit_x_five_sparse", (0, 1), 1, (0, 1, 6, 25, 31), 5)
    add("five_sparse_x_bit", (0, 3, 12, 21, 31), 5, (0, 1), 1)
    add("pair_x_four_weights", (0, 3), 2, (0, 1, 6, 9, 15), 4)
    add("four_weights_x_pair", (0, 5, 10, 15), 4, (1, 2), 2)
    add("triple_deterministic_x_triple", (1, 3, 5, 7), 3, (0, 2, 4, 6), 3)
    add("triple_antichain_x_triple_chain", (1, 2, 4), 3, (0, 1, 3, 7), 3)
    add("triple_dead_x_triple_dup", (0, 1, 2, 3), 3, (0, 3, 4, 7), 3)
    add("triple_uc_x_triple_nonuc", (0, 1, 2, 3), 3, (1, 2, 4), 3)
    add("triple_weight_cells_x_triple", (0, 7, 3, 5, 6), 3, (1, 2, 4, 7), 3)

    # Seeded random pairs across the same dimension budget.
    for index in range(12):
        left_dimension = random.choice((1, 2, 3, 4, 5))
        right_dimension = 6 - left_dimension
        if index % 4 == 3:
            right_dimension = max(1, right_dimension - random.randint(0, 2))
        left_rows = random.sample(
            range(1 << left_dimension),
            random.randint(1, min(1 << left_dimension, 8)),
        )
        right_rows = random.sample(
            range(1 << right_dimension),
            random.randint(1, min(1 << right_dimension, 8)),
        )
        add(
            f"random_{index}",
            tuple(left_rows),
            left_dimension,
            tuple(right_rows),
            right_dimension,
        )
    return tuple(cases)


def _read_checkpoint(path: Path) -> tuple[dict[str, dict[str, Any]], int, int]:
    records: dict[str, dict[str, Any]] = {}
    truncated = 0
    foreign = 0
    if not path.exists():
        return records, truncated, foreign
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                truncated += 1
                continue
            if payload.get("kind") != "case":
                continue
            if payload.get("schema") != RECORD_SCHEMA:
                foreign += 1
                continue
            records[payload["name"]] = payload
    return records, truncated, foreign


def _append_record(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def evaluate_case(case: dict[str, Any]) -> dict[str, Any]:
    left = case["left"]
    right = case["right"]
    left_dimension = case["left_dimension"]
    right_dimension = case["right_dimension"]
    dimension = left_dimension + right_dimension
    product = cartesian_product(left, left_dimension, right)
    assert len(product) == len(left) * len(right)

    # Lemma 1 in exact integers.
    left_defect = closure_defect(left)
    right_defect = closure_defect(right)
    product_defect = closure_defect(product)
    defect_identity = (1 - product_defect) == (1 - left_defect) * (
        1 - right_defect
    )

    # Lemma 2 on the order averages, in exact rationals.
    left_q = rational.shapley_iid_bounds(left, left_dimension)
    right_q = rational.shapley_iid_bounds(right, right_dimension)
    product_q = rational.shapley_iid_bounds(product, dimension)
    q_sum = (left_q[0] + right_q[0], left_q[1] + right_q[1])
    q_intersects = (
        product_q[0] <= q_sum[1] + TOLERANCE
        and q_sum[0] <= product_q[1] + TOLERANCE
    )
    q_gap = max(
        abs(product_q[1] - q_sum[1]), abs(product_q[0] - q_sum[0])
    )

    # Lemma 3 for every global order, in exact rationals.
    worst_gap = Fraction(0)
    failures = 0
    orders = 0
    factor_cache: dict[tuple[str, tuple[int, ...]], tuple[Fraction, Fraction]] = {}

    def factor_bounds(
        side: str,
        family: tuple[int, ...],
        side_dimension: int,
        relative: tuple[int, ...],
    ) -> tuple[Fraction, Fraction]:
        key = (side, relative)
        cached = factor_cache.get(key)
        if cached is None:
            low, high, _states = rational.bellman_bounds(
                family, side_dimension, relative
            )
            cached = (low, high)
            factor_cache[key] = cached
        return cached

    for order in permutations(range(dimension)):
        orders += 1
        low, high, _states = rational.bellman_bounds(product, dimension, order)
        left_relative, right_relative = induced_orders(order, left_dimension)
        left_bounds = factor_bounds(
            "left", left, left_dimension, _relabel(left_relative)
        )
        right_bounds = factor_bounds(
            "right", right, right_dimension, _relabel(right_relative)
        )
        sum_low = left_bounds[0] + right_bounds[0]
        sum_high = left_bounds[1] + right_bounds[1]
        if not (low <= sum_high + TOLERANCE and sum_low <= high + TOLERANCE):
            failures += 1
        worst_gap = max(
            worst_gap, abs(high - sum_high), abs(low - sum_low)
        )

    assert orders == factorial(dimension)
    return {
        "kind": "case",
        "schema": RECORD_SCHEMA,
        "name": case["name"],
        "left_rows": list(left),
        "right_rows": list(right),
        "left_dimension": left_dimension,
        "right_dimension": right_dimension,
        "product_size": len(product),
        "orders_checked": orders,
        "left_defect": str(left_defect),
        "right_defect": str(right_defect),
        "product_defect": str(product_defect),
        "defect_identity_exact": defect_identity,
        "q_enclosures_intersect": q_intersects,
        "q_worst_gap_decimal": rational.decimal_upper(q_gap, 30),
        "bellman_intersection_failures": failures,
        "bellman_worst_gap_decimal": rational.decimal_upper(worst_gap, 30),
    }


def audit(
    checkpoint: Path | None,
    result_path: Path | None,
    seed: int,
    progress: bool,
) -> dict[str, Any]:
    cases = _cases(seed)
    stored, truncated, foreign = (
        ({}, 0, 0) if checkpoint is None else _read_checkpoint(checkpoint)
    )

    records: list[dict[str, Any]] = []
    evaluated = 0
    resumed = 0
    started = time.monotonic()
    for case in cases:
        cached = stored.get(case["name"])
        if cached is not None and cached["left_rows"] == list(case["left"]):
            records.append(cached)
            resumed += 1
            continue
        record = evaluate_case(case)
        records.append(record)
        evaluated += 1
        if checkpoint is not None:
            _append_record(checkpoint, record)
        if progress:
            print(
                f"{case['name']}: {record['orders_checked']} orders, "
                f"failures {record['bellman_intersection_failures']}, "
                f"{time.monotonic() - started:.1f} s",
                file=sys.stderr,
                flush=True,
            )

    counterexamples = sum(
        record["bellman_intersection_failures"] for record in records
    )
    defect_failures = sum(
        0 if record["defect_identity_exact"] else 1 for record in records
    )
    q_failures = sum(
        0 if record["q_enclosures_intersect"] else 1 for record in records
    )
    summary = {
        "schema_version": 1,
        "verdict": (
            "NO_EXACT_COUNTEREXAMPLE"
            if counterexamples == 0 and defect_failures == 0 and q_failures == 0
            else "COUNTEREXAMPLE_FOUND"
        ),
        "arithmetic": "exact_python_fractions",
        "third_party_dependencies": [],
        "tolerance": str(TOLERANCE),
        "tolerance_role": (
            "The two sides are certified enclosures, so the lemma is attacked "
            "by demanding that they intersect.  The tolerance only absorbs the "
            "certificate width, which is below 1e-24, and is far smaller than "
            "any additivity defect a wrong lemma would produce."
        ),
        "case_count": len(records),
        "total_orders_checked": sum(
            record["orders_checked"] for record in records
        ),
        "unequal_dimension_cases": sum(
            1
            for record in records
            if record["left_dimension"] != record["right_dimension"]
        ),
        "bellman_intersection_failures": counterexamples,
        "defect_identity_failures": defect_failures,
        "q_intersection_failures": q_failures,
        "worst_bellman_gap_decimal": max(
            record["bellman_worst_gap_decimal"] for record in records
        ),
        "worst_q_gap_decimal": max(
            record["q_worst_gap_decimal"] for record in records
        ),
        "resumed_cases": resumed,
        "new_evaluations": evaluated,
        "truncated_checkpoint_lines": truncated,
        "foreign_schema_checkpoint_lines": foreign,
        "seed": seed,
        "cases": records,
    }
    if result_path is not None:
        result_path.parent.mkdir(parents=True, exist_ok=True)
        result_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--seed", type=int, default=20260827)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    summary = audit(
        None if args.fresh else args.checkpoint,
        None if args.fresh else args.result,
        args.seed,
        not args.quiet,
    )
    printable = {
        key: value for key, value in summary.items() if key != "cases"
    }
    print(json.dumps(printable, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
