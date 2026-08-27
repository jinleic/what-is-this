#!/usr/bin/env python3
"""Resumable non-factored square audit for the second Gate B base.

The certificate proves negativity for the 25-row six-coordinate family.  This
script materializes its 625-row Cartesian square as ordinary 12-bit masks and
runs the generic float64 Bellman and iid evaluators on five hostile global
orders.  No product-aware representation is passed to either evaluator.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from random import Random
import sys
import time
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from shapley_n6_shared_bellman import one_sided_costs  # noqa: E402

from audit_square_direct import exact_reimer_threshold, fixed_order_iid  # noqa: E402
import verify_gate_b_n6_dyadic as n6  # noqa: E402


BASE_DIMENSION = 6
CHECKPOINT = HERE / "experiments" / "n6_square_direct_checkpoint.jsonl"
SUMMARY = HERE / "candidates" / "n6_square_direct_audit.json"
TOLERANCE = 1e-9


def build_square(base: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sorted(
            left | (right << BASE_DIMENSION)
            for left in base
            for right in base
        )
    )


def exact_square_checks(
    base: tuple[int, ...],
    square: tuple[int, ...],
) -> dict[str, Any]:
    dimension = 2 * BASE_DIMENSION
    assert len(square) == len(base) ** 2 == 625
    counts = tuple(
        sum((row >> coordinate) & 1 for row in square)
        for coordinate in range(dimension)
    )
    assert counts == (250,) * dimension
    incidence = sum(counts)
    assert incidence == 3000
    threshold = exact_reimer_threshold(len(square))

    support = set(square)
    missing = sum(
        (left | right) not in support
        for left in square
        for right in square
    )
    defect = Fraction(missing, len(square) ** 2)
    assert missing == 357_864
    assert defect == 1 - Fraction(181, 625) ** 2

    columns = tuple(
        tuple((row >> coordinate) & 1 for row in square)
        for coordinate in range(dimension)
    )
    return {
        "dimension": dimension,
        "family_size": len(square),
        "coordinate_counts": list(counts),
        "cap_bound": 2 * len(square) // 5,
        "cap_holds_with_equality": counts == (2 * len(square) // 5,) * dimension,
        "total_incidence": incidence,
        "reimer_threshold": threshold,
        "reimer_holds_strictly": incidence > threshold,
        "missing_ordered_join_pairs": missing,
        "closure_defect": str(defect),
        "closure_defect_matches_product": True,
        "active": all(any(column) for column in columns),
        "separating": len(set(columns)) == dimension,
    }


def declared_orders() -> dict[str, tuple[int, ...]]:
    consecutive = tuple(range(12))
    blockswap = tuple(range(6, 12)) + tuple(range(6))
    reversed_order = tuple(reversed(range(12)))
    alternating = tuple(
        coordinate
        for local in range(6)
        for coordinate in (local, 6 + local)
    )
    random = Random(20260826)
    shuffled = list(range(12))
    random.shuffle(shuffled)
    return {
        "consecutive": consecutive,
        "blockswap": blockswap,
        "reversed": reversed_order,
        "alternating": alternating,
        "seeded_interleaved": tuple(shuffled),
    }


def induced(order: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left = tuple(coordinate for coordinate in order if coordinate < BASE_DIMENSION)
    right = tuple(
        coordinate - BASE_DIMENSION
        for coordinate in order
        if coordinate >= BASE_DIMENSION
    )
    return left, right


def load_records() -> dict[str, dict[str, Any]]:
    if not CHECKPOINT.exists():
        return {}
    lines = CHECKPOINT.read_text().splitlines()
    records: dict[str, dict[str, Any]] = {}
    for index, line in enumerate(lines):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            if index == len(lines) - 1:
                break
            raise
        records[record["order_name"]] = record
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", nargs="*")
    args = parser.parse_args()

    base = n6.reconstruct_family()
    square = build_square(base)
    exact = exact_square_checks(base, square)
    assert exact["cap_holds_with_equality"]
    assert exact["reimer_holds_strictly"]
    assert exact["closure_defect_matches_product"]
    assert exact["active"] and exact["separating"]

    orders = declared_orders()
    selected = args.orders or list(orders)
    unknown = [name for name in selected if name not in orders]
    assert not unknown, unknown

    records = load_records()
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    new_records = 0
    run_started = time.time()
    with CHECKPOINT.open("a") as handle:
        for name in selected:
            if name in records:
                continue
            order = orders[name]
            left_order, right_order = induced(order)
            square_costs, states = one_sided_costs(
                square,
                (order,),
                dimension=2 * BASE_DIMENSION,
            )
            left_costs, _left_states = one_sided_costs(
                base,
                (left_order,),
                dimension=BASE_DIMENSION,
            )
            right_costs, _right_states = one_sided_costs(
                base,
                (right_order,),
                dimension=BASE_DIMENSION,
            )
            square_iid = fixed_order_iid(square, order)
            base_iid = fixed_order_iid(base, left_order) + fixed_order_iid(
                base,
                right_order,
            )
            record = {
                "order_name": name,
                "global_order": list(order),
                "left_block_order": list(left_order),
                "right_block_order": list(right_order),
                "direct_square_bellman": square_costs[0],
                "base_block_bellman_sum": left_costs[0] + right_costs[0],
                "bellman_gap": abs(
                    square_costs[0] - left_costs[0] - right_costs[0]
                ),
                "direct_square_iid": square_iid,
                "base_block_iid_sum": base_iid,
                "iid_gap": abs(square_iid - base_iid),
                "bellman_state_count": states,
                "representation": "plain_12_bit_row_masks_no_factorization",
            }
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            records[name] = record
            new_records += 1

    evaluated = [records[name] for name in selected if name in records]
    summary = {
        "schema_version": 1,
        "status": "complete" if len(evaluated) == len(orders) else "partial",
        "evidence_status": "direct_numerical_product_check_float64",
        "float64_tolerance": TOLERANCE,
        "exact_square": exact,
        "checked_order_count": len(evaluated),
        "declared_order_count": len(orders),
        "worst_bellman_gap": max(
            (record["bellman_gap"] for record in evaluated),
            default=0.0,
        ),
        "worst_iid_gap": max(
            (record["iid_gap"] for record in evaluated),
            default=0.0,
        ),
        "all_within_tolerance": all(
            record["bellman_gap"] <= TOLERANCE
            and record["iid_gap"] <= TOLERANCE
            for record in evaluated
        ),
        "orders": {
            record["order_name"]: {
                "direct_square_bellman": record["direct_square_bellman"],
                "base_block_bellman_sum": record["base_block_bellman_sum"],
                "bellman_state_count": record["bellman_state_count"],
            }
            for record in evaluated
        },
    }
    report = {
        **summary,
        "new_evaluations": new_records,
        "runtime_seconds": round(time.time() - run_started, 2),
    }
    if len(evaluated) == len(orders):
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
