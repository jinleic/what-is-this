#!/usr/bin/env python3
"""Resumable complete nontrivial n=5 cap/Reimer audit for Gate B.

The exact C++ producer supplies one representative of every coordinate-
permutation orbit. This independent numerical audit evaluates all 120 orders
with the legacy float64 Bellman implementation and checkpoints every family.
Exact claims: producer coverage, filters, defects, and orbit counts. Numerical
claims: Q, C_plus, A_plus, and repair-ratio comparisons.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import log2
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import shapley_n5_global_coupling as n5
from shapley_entropy import ALPHA
from shapley_n6_shared_bellman import one_sided_costs


DIMENSION = 5
SCHEMA_VERSION = 1


def popcount(value: int) -> int:
    if hasattr(int, "bit_count"):
        return value.bit_count()
    return bin(value).count("1")


def coordinate_counts(family: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )


def normalized(family: tuple[int, ...]) -> bool:
    columns = tuple(
        tuple((row >> coordinate) & 1 for row in family)
        for coordinate in range(DIMENSION)
    )
    return (
        all(any(column) for column in columns)
        and len(set(columns)) == DIMENSION
    )


def closure_defect(family: tuple[int, ...]) -> Fraction:
    support = set(family)
    missing = sum(
        (left | right) not in support for left in family for right in family
    )
    return Fraction(missing, len(family) ** 2)


def union_closed(family: tuple[int, ...]) -> bool:
    support = set(family)
    return all((left | right) in support for left in family for right in family)


def evaluate(mask: int) -> dict[str, Any]:
    family = n5.members(mask)
    assert n5.exact_admissible(family)
    iid = n5.shapley_iid(family)
    roots, state_count = one_sided_costs(family, n5.ORDERS, dimension=DIMENSION)
    one_sided = sum(roots) / len(roots)
    a_plus = (1.0 - ALPHA) * iid + ALPHA * one_sided - log2(len(family))
    defect = closure_defect(family)
    ratio = -a_plus / float(defect) if a_plus < 0.0 and defect else None
    return {
        "family_mask_hex": hex(mask),
        "family_rows": list(family),
        "family_size": len(family),
        "coordinate_counts": list(coordinate_counts(family)),
        "total_incidence": sum(popcount(row) for row in family),
        "normalized": normalized(family),
        "union_closed": union_closed(family),
        "closure_defect": {
            "numerator": defect.numerator,
            "denominator": defect.denominator,
        },
        "shapley_iid": iid,
        "one_sided_cost": one_sided,
        "a_plus": a_plus,
        "repair_ratio": ratio,
        "shared_state_count": state_count,
        "evidence_status": "exact_feasible_complete_numerical_objective",
    }


def load_checkpoint(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    records: dict[str, dict[str, Any]] = {}
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.endswith("\n"):
                print(f"Ignoring truncated checkpoint line {line_number}")
                break
            entry = json.loads(line)
            assert entry["schema_version"] == SCHEMA_VERSION
            assert entry["dimension"] == DIMENSION
            record = entry["record"]
            key = record["family_mask_hex"]
            assert key not in records
            records[key] = record
    return records


def append_checkpoint(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "schema_version": SCHEMA_VERSION,
        "dimension": DIMENSION,
        "record": record,
    }
    with path.open("a") as stream:
        stream.write(json.dumps(entry, sort_keys=True) + "\n")
        stream.flush()


def build_result(
    masks: tuple[int, ...], records: tuple[dict[str, Any], ...]
) -> dict[str, Any]:
    assert len(masks) == len(records) == n5.EXPECTED_CANONICAL_TOTAL
    negative = tuple(
        record for record in records if record["a_plus"] < -1e-10
    )
    near_zero = tuple(
        record for record in records if abs(record["a_plus"]) <= 1e-10
    )
    union_closed_records = tuple(
        record for record in records if record["union_closed"]
    )
    worst = min(records, key=lambda record: record["a_plus"])
    assert tuple(worst["family_rows"]) == (0, 3, 4, 8, 15)
    assert abs(worst["a_plus"] - 0.006474313148642441) < 4e-12
    assert not negative
    assert not near_zero
    assert not union_closed_records
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "evidence_status": "complete_exact_census_complete_numerical_objective",
        "dimension": DIMENSION,
        "nontrivial_only": True,
        "family_size_range": [min(n5.SPECS), max(n5.SPECS)],
        "coordinate_order_count": len(n5.ORDERS),
        "represented_labeled_family_count": n5.EXPECTED_LABELED_TOTAL,
        "canonical_family_count": len(masks),
        "family_size_counts": dict(
            sorted(
                Counter(record["family_size"] for record in records).items()
            )
        ),
        "negative_a_plus_count": len(negative),
        "near_zero_a_plus_count": len(near_zero),
        "union_closed_count": len(union_closed_records),
        "worst_a_plus": worst,
    }


def main() -> None:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser()
    parser.add_argument("enumerator", type=Path)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=root / "experiments" / "n5_complete_checkpoint.jsonl",
    )
    parser.add_argument(
        "--result",
        type=Path,
        default=root / "candidates" / "n5_complete_audit.json",
    )
    parser.add_argument(
        "--max-evaluations",
        type=int,
        help="Bound new evaluations in this invocation; omit for completion.",
    )
    args = parser.parse_args()
    if not args.enumerator.is_file():
        raise SystemExit(f"missing enumerator: {args.enumerator}")

    masks, enumerator_log = n5.read_enumerator(args.enumerator)
    assert len(masks) == len(set(masks)) == n5.EXPECTED_CANONICAL_TOTAL
    assert Counter(len(n5.members(mask)) for mask in masks) == Counter(
        n5.EXPECTED_CANONICAL
    )
    assert all(n5.exact_admissible(n5.members(mask)) for mask in masks)
    assert "total=463343 canonical_total=5172" in enumerator_log

    checkpoint = load_checkpoint(args.checkpoint)
    expected_keys = {hex(mask) for mask in masks}
    assert set(checkpoint) <= expected_keys

    new_evaluations = 0
    for mask in masks:
        key = hex(mask)
        if key in checkpoint:
            continue
        if (
            args.max_evaluations is not None
            and new_evaluations >= args.max_evaluations
        ):
            break
        record = evaluate(mask)
        append_checkpoint(args.checkpoint, record)
        checkpoint[key] = record
        new_evaluations += 1
        if new_evaluations % 100 == 0:
            print(
                f"evaluated {len(checkpoint)} / {len(masks)} "
                f"({new_evaluations} new)",
                flush=True,
            )

    complete = len(checkpoint) == len(masks)
    print(json.dumps({
        "dimension": DIMENSION,
        "canonical_family_count": len(masks),
        "completed_evaluations": len(checkpoint),
        "new_evaluations": new_evaluations,
        "complete": complete,
    }, sort_keys=True))
    if not complete:
        print("PARTIAL: rerun the same command to resume")
        return

    records = tuple(checkpoint[hex(mask)] for mask in masks)
    result = build_result(masks, records)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.result.exists():
        assert args.result.read_text() == encoded
    else:
        args.result.parent.mkdir(parents=True, exist_ok=True)
        args.result.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
