#!/usr/bin/env python3
"""Direct, non-factored audit of the Gate B square family.

``verify_gate_b.py`` evaluates the 2,025-row square with an exact factored-fiber
representation, which already encodes the product structure it is testing.  This
script instead materializes the square as a plain list of 14-bit row masks and
runs the generic prefix-fiber Bellman evaluator and a generic fixed-order iid
evaluator on it, with no product-aware shortcut.  It then compares the results
against the sums of the corresponding base-order values.

Each declared global order is one append-only checkpoint record, so a long run
resumes exactly where it stopped.

Run from ``math/``:

    OMP_NUM_THREADS=1 nice -n 10 ./.venv/bin/python -B \\
        uc/gate_b/audit_square_direct.py [--orders NAME ...]
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

from shapley_global_coupling import binary_entropy  # noqa: E402
from shapley_n6_shared_bellman import one_sided_costs  # noqa: E402


CANDIDATE = HERE / "candidates" / "n7_block_extremizer.json"
CHECKPOINT = HERE / "experiments" / "square_direct_checkpoint.jsonl"
SUMMARY = HERE / "candidates" / "square_direct_audit.json"
BASE_DIMENSION = 7
TOLERANCE = 1e-9


def exact_reimer_threshold(size: int) -> int:
    low, high = 0, size * size
    target = size**size
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def fixed_order_iid(family: tuple[int, ...], order: tuple[int, ...]) -> float:
    """Fixed-order iid OR-entropy with probability-merged conditional laws."""
    total = 0.0
    for position, coordinate in enumerate(order):
        predecessors = order[:position]
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
        law: dict[Fraction, Fraction] = {}
        for count, ones in groups.values():
            probability = Fraction(ones, count)
            law[probability] = law.get(probability, Fraction(0)) + Fraction(
                count, size
            )
        assert sum(law.values()) == 1
        total += sum(
            float(left_mass * right_mass)
            * binary_entropy(left + right - left * right)
            for left, left_mass in law.items()
            for right, right_mass in law.items()
        )
    return total


def build_square(base: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sorted(
            left | (right << BASE_DIMENSION) for left in base for right in base
        )
    )


def exact_square_checks(
    base: tuple[int, ...], square: tuple[int, ...]
) -> dict[str, Any]:
    dimension = 2 * BASE_DIMENSION
    assert len(square) == len(base) ** 2 == 2025
    counts = tuple(
        sum((row >> coordinate) & 1 for row in square)
        for coordinate in range(dimension)
    )
    assert counts == (810,) * dimension
    incidence = sum(counts)
    threshold = exact_reimer_threshold(len(square))
    support = set(square)
    missing = sum(
        (left | right) not in support for left in square for right in square
    )
    defect = Fraction(missing, len(square) ** 2)
    columns = tuple(
        tuple((row >> coordinate) & 1 for row in square)
        for coordinate in range(dimension)
    )
    return {
        "dimension": dimension,
        "family_size": len(square),
        "coordinate_counts": list(counts),
        "cap_bound": 2 * len(square) // 5,
        "cap_holds": max(counts) <= 2 * len(square) // 5,
        "total_incidence": incidence,
        "reimer_threshold": threshold,
        "reimer_holds_strictly": incidence > threshold,
        "missing_ordered_join_pairs": missing,
        "closure_defect": str(defect),
        "closure_defect_matches_product": defect
        == 1 - Fraction(17, 81) ** 2,
        "active": all(any(column) for column in columns),
        "separating": len(set(columns)) == dimension,
    }


def declared_orders() -> dict[str, tuple[int, ...]]:
    consecutive = tuple(range(14))
    blockswap = tuple(range(7, 14)) + tuple(range(7))
    reversed_order = tuple(reversed(range(14)))
    alternating = tuple(
        coordinate
        for local in range(7)
        for coordinate in (local, 7 + local)
    )
    random = Random(20260826)
    shuffled = list(range(14))
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
    records = {}
    for line in CHECKPOINT.read_text().splitlines():
        if line.strip():
            record = json.loads(line)
            records[record["order_name"]] = record
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", nargs="*")
    args = parser.parse_args()

    base = tuple(json.loads(CANDIDATE.read_text())["family_rows"])
    square = build_square(base)
    exact = exact_square_checks(base, square)
    assert exact["cap_holds"]
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
    with CHECKPOINT.open("a") as handle:
        for name in selected:
            if name in records:
                continue
            order = orders[name]
            left_order, right_order = induced(order)
            started = time.time()
            square_costs, states = one_sided_costs(
                square, (order,), dimension=14
            )
            left_costs, _left_states = one_sided_costs(
                base, (left_order,), dimension=BASE_DIMENSION
            )
            right_costs, _right_states = one_sided_costs(
                base, (right_order,), dimension=BASE_DIMENSION
            )
            square_iid = fixed_order_iid(square, order)
            base_iid = fixed_order_iid(base, left_order) + fixed_order_iid(
                base, right_order
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
                "representation": "plain_14_bit_row_masks_no_factorization",
                "seconds": round(time.time() - started, 2),
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
            (record["bellman_gap"] for record in evaluated), default=0.0
        ),
        "worst_iid_gap": max(
            (record["iid_gap"] for record in evaluated), default=0.0
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
                "seconds": record["seconds"],
            }
            for record in evaluated
        },
    }
    # The persisted artifact excludes the run-dependent counter so a resume of
    # a complete audit reproduces it byte-identically.
    report = {**summary, "new_evaluations": new_records}
    if len(evaluated) == len(orders):
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
