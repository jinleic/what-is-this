#!/usr/bin/env python3
"""Checkpointed Gate B search in complete n=8 block-symmetric classes.

For a fixed k, enumerate every family invariant under S_k x S_(8-k), filter by
the authoritative cap-2/5 and Reimer constraints, quotient the full S_8 action,
and evaluate exact block-order representatives.  Enumeration, constraints,
canonicalization, order coverage, and closure defects are exact.  Entropy and
Bellman values are float64 and are therefore discovery evidence only.

The default k=1 class has 2^16-1 raw masks and eight order representatives per
family.  Checkpoints are append-only JSONL: every complete line is durable and
a truncated final line is ignored on resume.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import comb, log2
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from shapley_entropy import ALPHA
from shapley_join_loss import shapley_iid
from shapley_n6_shared_bellman import one_sided_costs, sequential_costs


DIMENSION = 8
SUBSET_COUNT = 1 << DIMENSION


def exact_reimer_threshold(size: int) -> int:
    if size <= 1:
        return 0
    target = size**size
    low = 0
    high = size * (size - 1).bit_length()
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def members(mask: int) -> tuple[int, ...]:
    return tuple(row for row in range(SUBSET_COUNT) if mask >> row & 1)


def family_mask(family: tuple[int, ...]) -> int:
    return sum(1 << row for row in family)


def coordinate_counts(family: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )


def exact_admissible(family: tuple[int, ...]) -> bool:
    size = len(family)
    if size < 3:
        return False
    counts = coordinate_counts(family)
    return (
        max(counts) <= 2 * size // 5
        and sum(counts) >= exact_reimer_threshold(size)
    )


def normalized(family: tuple[int, ...]) -> bool:
    columns = tuple(
        tuple((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    return all(any(column) for column in columns) and len(set(columns)) == DIMENSION


def permute_subset(subset: int, permutation: tuple[int, ...]) -> int:
    image = 0
    for old_coordinate, new_coordinate in enumerate(permutation):
        if subset >> old_coordinate & 1:
            image |= 1 << new_coordinate
    return image


def block_cells(k: int) -> tuple[tuple[int, ...], ...]:
    left_mask = (1 << k) - 1
    return tuple(
        tuple(
            row
            for row in range(SUBSET_COUNT)
            if (row & left_mask).bit_count() == left_weight
            and (row >> k).bit_count() == right_weight
        )
        for left_weight in range(k + 1)
        for right_weight in range(DIMENSION - k + 1)
    )


def canonical_block_family(family: tuple[int, ...], k: int) -> int:
    coordinates = tuple(range(DIMENSION))
    images = []
    for target_left in combinations(coordinates, k):
        target_right = tuple(
            coordinate for coordinate in coordinates if coordinate not in target_left
        )
        permutation = tuple(target_left) + target_right
        images.append(
            family_mask(
                tuple(permute_subset(row, permutation) for row in family)
            )
        )
    return min(images)


def enumerate_canonical(k: int) -> tuple[tuple[int, ...], int, int]:
    cells = block_cells(k)
    cell_masks = tuple(family_mask(cell) for cell in cells)
    cell_stats = tuple(
        (
            len(cell),
            sum((row >> 0) & 1 for row in cell),
            sum((row >> k) & 1 for row in cell),
            sum(row.bit_count() for row in cell),
        )
        for cell in cells
    )

    canonical: dict[int, tuple[int, ...]] = {}
    accepted_count = 0
    current_family = 0
    size = left_count = right_count = incidence = 0
    previous_gray = 0
    for step in range(1, 1 << len(cells)):
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

        if size < 3:
            continue
        cap = 2 * size // 5
        if (
            max(left_count, right_count) > cap
            or incidence < exact_reimer_threshold(size)
        ):
            continue
        family = members(current_family)
        assert exact_admissible(family)
        accepted_count += 1
        canonical_mask = canonical_block_family(family, k)
        canonical[canonical_mask] = members(canonical_mask)

    families = tuple(canonical[mask] for mask in sorted(canonical))
    return families, len(cells), accepted_count


def block_order_representatives(k: int) -> tuple[tuple[int, ...], ...]:
    """One order per S_k x S_(8-k) orbit, represented by its A/B pattern."""
    result = []
    positions = tuple(range(DIMENSION))
    for left_positions in combinations(positions, k):
        left_positions = set(left_positions)
        left_coordinate = 0
        right_coordinate = k
        order = []
        for position in positions:
            if position in left_positions:
                order.append(left_coordinate)
                left_coordinate += 1
            else:
                order.append(right_coordinate)
                right_coordinate += 1
        result.append(tuple(order))
    assert len(result) == comb(DIMENSION, k)
    return tuple(result)


def closure_defect(family: tuple[int, ...]) -> Fraction:
    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    return Fraction(missing, len(family) ** 2)


def evaluate(
    family: tuple[int, ...], orders: tuple[tuple[int, ...], ...]
) -> dict[str, Any]:
    iid = shapley_iid(family, DIMENSION)
    sequential_roots, sequential_states = sequential_costs(
        family, orders, dimension=DIMENSION
    )
    one_sided_roots, one_sided_states = one_sided_costs(
        family, orders, dimension=DIMENSION
    )
    sequential = sum(sequential_roots) / len(sequential_roots)
    one_sided = sum(one_sided_roots) / len(one_sided_roots)
    entropy = log2(len(family))
    a_seq = (1.0 - ALPHA) * iid + ALPHA * sequential - entropy
    a_plus = (1.0 - ALPHA) * iid + ALPHA * one_sided - entropy
    defect = closure_defect(family)
    ratio = -a_plus / float(defect) if a_plus < 0.0 and defect else None
    return {
        "family_mask_hex": hex(family_mask(family)),
        "family_rows": list(family),
        "family_size": len(family),
        "coordinate_counts": list(coordinate_counts(family)),
        "total_incidence": sum(row.bit_count() for row in family),
        "normalized": normalized(family),
        "closure_defect": {
            "numerator": defect.numerator,
            "denominator": defect.denominator,
        },
        "shapley_iid": iid,
        "sequential_cost": sequential,
        "one_sided_cost": one_sided,
        "a_seq": a_seq,
        "a_plus": a_plus,
        "repair_ratio": ratio,
        "shared_state_count": max(sequential_states, one_sided_states),
        "evidence_status": "exact_feasible_complete_numerical_objective",
    }


def load_checkpoint(path: Path, k: int) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    records: dict[str, dict[str, Any]] = {}
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.endswith("\n"):
                print(f"Ignoring truncated checkpoint line {line_number}")
                break
            entry = json.loads(line)
            assert entry["schema_version"] == 1
            assert entry["dimension"] == DIMENSION
            assert entry["k"] == k
            record = entry["record"]
            records[record["family_mask_hex"]] = record
    return records


def append_checkpoint(path: Path, k: int, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "schema_version": 1,
        "dimension": DIMENSION,
        "k": k,
        "record": record,
    }
    with path.open("a") as stream:
        stream.write(json.dumps(entry, sort_keys=True) + "\n")
        stream.flush()


def build_result(
    k: int,
    cell_count: int,
    raw_accepted_count: int,
    families: tuple[tuple[int, ...], ...],
    records: tuple[dict[str, Any], ...],
) -> dict[str, Any]:
    negative = tuple(record for record in records if record["a_plus"] < -1e-10)
    ratio_records = tuple(
        record for record in records if record["repair_ratio"] is not None
    )
    normalized_records = tuple(record for record in records if record["normalized"])
    worst = min(records, key=lambda record: record["a_plus"])
    ratio_worst = max(
        ratio_records, key=lambda record: record["repair_ratio"], default=None
    )
    normalized_worst = min(
        normalized_records, key=lambda record: record["a_plus"], default=None
    )
    return {
        "schema_version": 1,
        "status": "complete",
        "evidence_status": "complete_exact_class_complete_numerical_objective",
        "dimension": DIMENSION,
        "k": k,
        "cell_count": cell_count,
        "raw_mask_count": (1 << cell_count) - 1,
        "raw_accepted_count": raw_accepted_count,
        "canonical_family_count": len(families),
        "normalized_family_count": len(normalized_records),
        "order_representative_count": comb(DIMENSION, k),
        "negative_a_plus_count": len(negative),
        "worst_a_plus": worst,
        "largest_repair_ratio": ratio_worst,
        "normalized_worst_a_plus": normalized_worst,
    }


def main() -> None:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, choices=(1, 2), default=1)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=root / "experiments" / "n8_k1_checkpoint.jsonl",
    )
    parser.add_argument(
        "--result",
        type=Path,
        default=root / "candidates" / "n8_k1_census.json",
    )
    parser.add_argument(
        "--max-evaluations",
        type=int,
        help="Bound new evaluations in this invocation; omit for completion.",
    )
    args = parser.parse_args()

    families, cell_count, raw_accepted_count = enumerate_canonical(args.k)
    orders = block_order_representatives(args.k)
    checkpoint = load_checkpoint(args.checkpoint, args.k)
    expected_masks = {hex(family_mask(family)) for family in families}
    assert set(checkpoint) <= expected_masks

    new_evaluations = 0
    for index, family in enumerate(families, 1):
        key = hex(family_mask(family))
        if key in checkpoint:
            continue
        if (
            args.max_evaluations is not None
            and new_evaluations >= args.max_evaluations
        ):
            break
        record = evaluate(family, orders)
        append_checkpoint(args.checkpoint, args.k, record)
        checkpoint[key] = record
        new_evaluations += 1
        if new_evaluations % 10 == 0:
            print(
                f"evaluated {len(checkpoint)} / {len(families)} "
                f"({new_evaluations} new)",
                flush=True,
            )

    complete = len(checkpoint) == len(families)
    print(
        json.dumps(
            {
                "dimension": DIMENSION,
                "k": args.k,
                "cell_count": cell_count,
                "raw_accepted_count": raw_accepted_count,
                "canonical_family_count": len(families),
                "completed_evaluations": len(checkpoint),
                "new_evaluations": new_evaluations,
                "complete": complete,
            },
            sort_keys=True,
        )
    )
    if not complete:
        print("PARTIAL: rerun the same command to resume")
        return

    records = tuple(checkpoint[hex(family_mask(family))] for family in families)
    result = build_result(
        args.k, cell_count, raw_accepted_count, families, records
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.result.exists():
        assert args.result.read_text() == encoded
    else:
        args.result.parent.mkdir(parents=True, exist_ok=True)
        args.result.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
