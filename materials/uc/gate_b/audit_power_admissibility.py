#!/usr/bin/env python3
"""Exact admissibility and defect audit of the Gate B powers by construction.

The theorem's admissibility and defect statements for \\(\\mathcal F_k\\) are
proved from the product lemmas, and the standalone Arb checker verifies the
resulting integer arithmetic for \\(k\\le8\\). This script closes the remaining
gap between "formula arithmetic" and "actual family" by building the powers
explicitly and recomputing everything from the constructed row set, including a
brute-force count of missing ordered joins with no product-aware shortcut.

Exact for \\(k\\le3\\): \\(45^3=91{,}125\\) rows give \\(8.3\\times10^9\\)
ordered pairs, counted with a dense membership table over the \\(2^{21}\\)
subsets. Records append to a JSONL checkpoint, one per power.

Run from ``math/``:

    OMP_NUM_THREADS=1 nice -n 10 ./.venv/bin/python -B \\
        uc/gate_b/audit_power_admissibility.py [--max-power K]
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
import time
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent

CANDIDATE = HERE / "candidates" / "n7_block_extremizer.json"
CHECKPOINT = HERE / "experiments" / "power_admissibility_checkpoint.jsonl"
SUMMARY = HERE / "candidates" / "power_admissibility_audit.json"
BASE_DIMENSION = 7
DEFAULT_MAX_POWER = 3


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


def power_rows(base: tuple[int, ...], power: int) -> np.ndarray:
    rows = np.zeros(1, dtype=np.int64)
    base_array = np.array(base, dtype=np.int64)
    for block in range(power):
        shifted = base_array << (block * BASE_DIMENSION)
        rows = (rows[:, None] | shifted[None, :]).reshape(-1)
    rows.sort()
    return rows


def missing_ordered_joins(rows: np.ndarray, dimension: int) -> int:
    """Brute-force count of ordered pairs whose union leaves the family."""
    present = np.zeros(1 << dimension, dtype=bool)
    present[rows] = True
    missing = 0
    for row in rows.tolist():
        missing += int(np.count_nonzero(~present[rows | row]))
    return missing


def audit_power(base: tuple[int, ...], power: int) -> dict[str, Any]:
    started = time.time()
    dimension = BASE_DIMENSION * power
    rows = power_rows(base, power)
    size = int(rows.size)
    assert size == len(base) ** power
    assert np.unique(rows).size == size

    bits = ((rows[:, None] >> np.arange(dimension)[None, :]) & 1).astype(np.int64)
    counts = bits.sum(axis=0)
    incidence = int(counts.sum())
    cap_bound = 2 * size // 5
    threshold = exact_reimer_threshold(size)
    columns = {tuple(column) for column in bits.T.tolist()}

    missing = missing_ordered_joins(rows, dimension)
    defect = Fraction(missing, size * size)
    predicted = 1 - Fraction(17, 81) ** power

    return {
        "power": power,
        "dimension": dimension,
        "family_size": size,
        "distinct_coordinate_counts": sorted({int(count) for count in counts}),
        "expected_coordinate_count": 18 * len(base) ** (power - 1),
        "cap_bound": cap_bound,
        "cap_holds_with_equality": int(counts.max()) == cap_bound,
        "total_incidence": incidence,
        "expected_total_incidence": 126 * power * len(base) ** (power - 1),
        "reimer_threshold": threshold,
        "reimer_holds_strictly": incidence > threshold,
        "active": bool(np.all(counts > 0)),
        "separating": len(columns) == dimension,
        "missing_ordered_join_pairs": missing,
        "closure_defect": str(defect),
        "predicted_closure_defect": str(predicted),
        "defect_matches_product_formula": defect == predicted,
        "brute_force_pair_count": size * size,
        "seconds": round(time.time() - started, 2),
    }


def load_records() -> dict[int, dict[str, Any]]:
    if not CHECKPOINT.exists():
        return {}
    records = {}
    for line in CHECKPOINT.read_text().splitlines():
        if line.strip():
            record = json.loads(line)
            records[record["power"]] = record
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-power", type=int, default=DEFAULT_MAX_POWER)
    args = parser.parse_args()

    base = tuple(json.loads(CANDIDATE.read_text())["family_rows"])
    records = load_records()
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    new_records = 0
    with CHECKPOINT.open("a") as handle:
        for power in range(1, args.max_power + 1):
            if power in records:
                continue
            record = audit_power(base, power)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            records[power] = record
            new_records += 1

    evaluated = [records[power] for power in range(1, args.max_power + 1)]
    summary = {
        "schema_version": 1,
        "status": "complete"
        if args.max_power >= DEFAULT_MAX_POWER
        else "partial",
        "evidence_status": "exact_constructed_family_audit",
        "max_power": args.max_power,
        "all_admissible": all(
            record["cap_holds_with_equality"]
            and record["reimer_holds_strictly"]
            and record["active"]
            and record["separating"]
            for record in evaluated
        ),
        "all_defects_match_product_formula": all(
            record["defect_matches_product_formula"] for record in evaluated
        ),
        "all_counts_uniform": all(
            len(record["distinct_coordinate_counts"]) == 1
            and record["distinct_coordinate_counts"][0]
            == record["expected_coordinate_count"]
            for record in evaluated
        ),
        "all_incidences_match": all(
            record["total_incidence"] == record["expected_total_incidence"]
            for record in evaluated
        ),
        "powers": {str(record["power"]): record for record in evaluated},
    }
    # The persisted artifact excludes the run-dependent counter so a resume of
    # a complete audit reproduces it byte-identically.
    report = {**summary, "new_evaluations": new_records}
    if args.max_power >= DEFAULT_MAX_POWER:
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
