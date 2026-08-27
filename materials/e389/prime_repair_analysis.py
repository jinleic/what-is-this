#!/usr/bin/env python3
"""Discriminate original-prime repair from new-prime interference."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import RatioSlackScanner, SPFFactorizer, slack

ROOT = Path(__file__).resolve().parent
DEFAULT_ATLAS = ROOT / "data" / "atlas_m1_50_k50000.json"
DEFAULT_OUTPUT = ROOT / "data" / "prime_repair_analysis_r64.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ordered_nonzero_deltas(radius: int) -> list[int]:
    return [delta for distance in range(1, radius + 1) for delta in (-distance, distance)]


def witness_bitmaps(min_m: int, max_m: int, max_k: int) -> list[bytearray | None]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    bitmaps: list[bytearray | None] = [None] * (max_m + 1)
    for m in range(min_m, max_m + 1):
        scanner = RatioSlackScanner(m, factorizer)
        bits = bytearray(max_k + 1)
        while True:
            if scanner.is_witness:
                bits[scanner.k] = 1
            if scanner.k == max_k:
                break
            scanner.advance()
        bitmaps[m] = bits
    return bitmaps


def analyze(
    *, min_m: int, max_m: int, max_k: int, radius: int, atlas_path: Path
) -> dict[str, Any]:
    if radius < 1:
        raise ValueError("radius must be positive")
    atlas_bytes = atlas_path.read_bytes()
    atlas = json.loads(atlas_bytes)
    if atlas["schema_version"] != 2:
        raise AssertionError("prime repair requires repaired atlas schema 2")
    for name, value in (("min_m", min_m), ("max_m", max_m), ("max_k", max_k)):
        if atlas["parameters"][name] != value:
            raise AssertionError(f"atlas parameter mismatch: {name}")

    factorizer = SPFFactorizer(max_m + 2 * max_k)
    bitmaps = witness_bitmaps(min_m, max_m, max_k)
    for row in atlas["m_rows"]:
        bits = bitmaps[row["m"]]
        assert bits is not None
        if sum(bits) != row["witness_count"]:
            raise AssertionError(f"witness bitmap mismatch at m={row['m']}")

    deltas = ordered_nonzero_deltas(radius)
    counts = Counter()
    by_m: dict[int, Counter[str]] = {}
    cases: list[dict[str, Any]] = []

    for m in range(min_m, max_m):
        source = bitmaps[m]
        target = bitmaps[m + 1]
        assert source is not None and target is not None
        row_counts: Counter[str] = Counter()
        by_m[m] = row_counts
        for k in range(2, max_k + 1):
            if not source[k] or target[k - 1]:
                continue
            counts["failed_natural_shifts"] += 1
            row_counts["failed_natural_shifts"] += 1
            base = k - 1
            if base - radius < 1 or base + radius > max_k:
                counts["boundary_excluded"] += 1
                row_counts["boundary_excluded"] += 1
                continue

            obstruction_primes = [
                p
                for p in sorted(factorizer.factor(m + 2 * k))
                if slack(m + 1, base, p) < 0
            ]
            if not obstruction_primes:
                raise AssertionError("failed natural shift had no exact obstruction prime")

            joint_delta = None
            actual_delta = None
            for delta in deltas:
                candidate_k = base + delta
                if joint_delta is None and all(
                    slack(m + 1, candidate_k, p) >= 0 for p in obstruction_primes
                ):
                    joint_delta = delta
                if actual_delta is None and target[candidate_k]:
                    actual_delta = delta
                if joint_delta is not None and actual_delta is not None:
                    break

            if actual_delta is not None and joint_delta is None:
                raise AssertionError("actual witness did not repair original obstructions")
            if actual_delta == joint_delta and actual_delta is not None:
                category = "exact_original_prime_alignment"
            elif actual_delta is not None:
                category = "new_prime_interference_delay"
            elif joint_delta is not None:
                category = "original_primes_repaired_but_no_witness"
            else:
                category = "original_primes_unrepaired_within_radius"
            counts[category] += 1
            row_counts[category] += 1

            cases.append(
                {
                    "m": m,
                    "k": k,
                    "target_base_k": base,
                    "obstruction_primes": obstruction_primes,
                    "joint_original_prime_repair_delta": joint_delta,
                    "actual_witness_repair_delta": actual_delta,
                    "category": category,
                }
            )

    eligible = counts["failed_natural_shifts"] - counts["boundary_excluded"]
    categorized = sum(
        counts[name]
        for name in (
            "exact_original_prime_alignment",
            "new_prime_interference_delay",
            "original_primes_repaired_but_no_witness",
            "original_primes_unrepaired_within_radius",
        )
    )
    if eligible != categorized or eligible != len(cases):
        raise AssertionError("prime-repair category partition mismatch")

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_PRIME_REPAIR_ANALYSIS",
        "claim_boundary": (
            "For each boundary-comparable failed natural shift, deltas are searched "
            "only within the declared radius. Repairing the original obstruction "
            "primes is necessary but may expose unrelated primes."
        ),
        "parameters": {
            "min_m": min_m,
            "max_m": max_m,
            "max_k": max_k,
            "radius": radius,
            "delta_order": "-1,+1,-2,+2,...",
        },
        "counts": dict(sorted(counts.items())),
        "by_m": [
            {"m": m, **dict(sorted(row_counts.items()))}
            for m, row_counts in sorted(by_m.items())
            if row_counts["failed_natural_shifts"]
        ],
        "cases": cases,
        "input": {
            "atlas_logical_name": atlas_path.name,
            "atlas_sha256": hashlib.sha256(atlas_bytes).hexdigest(),
        },
        "implementation_sha256": {
            "prime_repair_analysis.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-m", type=int, default=1)
    parser.add_argument("--max-m", type=int, default=50)
    parser.add_argument("--max-k", type=int, default=50_000)
    parser.add_argument("--radius", type=int, default=64)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        min_m=args.min_m,
        max_m=args.max_m,
        max_k=args.max_k,
        radius=args.radius,
        atlas_path=args.atlas.resolve(),
    )
    result["resources"] = {
        "wall_seconds": time.perf_counter() - started_wall,
        "process_cpu_seconds": time.process_time() - started_cpu,
        "processes": 1,
    }
    atomic_write_json(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "counts": result["counts"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
