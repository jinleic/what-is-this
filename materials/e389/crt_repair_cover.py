#!/usr/bin/env python3
"""Separate exact small-prime CRT repair from moving-window obstruction repair."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import SPFFactorizer, primes_up_to, slack
from residue_control_analysis import band_residue_slack
from shift_repair_analysis import scan_witness_bits

ROOT = Path(__file__).resolve().parent
DEFAULT_ATLAS = ROOT / "data" / "atlas_m1_50_k50000.json"
DEFAULT_OUTPUT = ROOT / "data" / "crt_repair_cover_m1_50_k50000_r64.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ordered_deltas(radius: int) -> list[int]:
    if radius < 0:
        raise ValueError("radius must be nonnegative")
    return [0] + [delta for distance in range(1, radius + 1) for delta in (-distance, distance)]


def _band_exponent(x: int, p: int) -> int:
    exponent = 0
    power = p
    while power <= x:
        exponent += 1
        power *= p
    if exponent == 0:
        raise AssertionError("controlled prime unexpectedly exceeds x")
    return exponent


def analyze(
    *,
    min_m: int,
    max_m: int,
    max_k: int,
    radius: int,
    atlas_path: Path,
) -> dict[str, Any]:
    if not 1 <= min_m < max_m or max_k < 2 or radius < 1:
        raise ValueError("require 1<=min_m<max_m, max_k>=2, and radius>=1")
    atlas_path = atlas_path.resolve()
    atlas_bytes = atlas_path.read_bytes()
    atlas = json.loads(atlas_bytes)
    if atlas["schema_version"] != 2:
        raise AssertionError("CRT repair cover requires repaired atlas schema 2")
    for name, expected in (("min_m", min_m), ("max_m", max_m), ("max_k", max_k)):
        if atlas["parameters"][name] != expected:
            raise AssertionError(f"atlas parameter mismatch: {name}")

    bitmaps = scan_witness_bits(min_m, max_m, max_k)
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    primes = primes_up_to(max_m)
    deltas = ordered_deltas(radius)
    counts = Counter()
    category_counts = Counter()
    tier_delta_counts = Counter()
    witness_delay_counts = Counter()
    blocker_counts = Counter()
    cases: list[dict[str, Any]] = []

    for source_m in range(min_m, max_m):
        source = bitmaps[source_m]
        target = bitmaps[source_m + 1]
        assert source is not None and target is not None
        target_m = source_m + 1
        controlled_primes = [p for p in primes if p <= target_m]
        for source_k in range(2, max_k + 1):
            if not source[source_k] or target[source_k - 1]:
                continue
            counts["failed_natural_shifts"] += 1
            base_k = source_k - 1
            if base_k - radius < 1 or base_k + radius > max_k:
                counts["boundary_excluded"] += 1
                continue

            first_tier_delta: int | None = None
            first_witness_delta: int | None = None
            for delta in deltas:
                candidate_k = base_k + delta
                x = target_m + 2 * candidate_k
                small_ok = True
                for p in controlled_primes:
                    exponent = _band_exponent(x, p)
                    actual = slack(target_m, candidate_k, p)
                    predicted = band_residue_slack(
                        target_m, p, exponent, candidate_k
                    )
                    if actual != predicted:
                        raise AssertionError("band-exact CRT slack mismatch")
                    counts["band_residue_checks"] += 1
                    if actual < 0:
                        small_ok = False
                        break
                if small_ok and first_tier_delta is None:
                    first_tier_delta = delta
                if target[candidate_k] and first_witness_delta is None:
                    if not small_ok:
                        raise AssertionError("witness failed the small-prime tier")
                    first_witness_delta = delta
                if first_tier_delta is not None and first_witness_delta is not None:
                    break

            if first_tier_delta is None:
                if first_witness_delta is not None:
                    raise AssertionError("witness found without a tier-one repair")
                category = "no_small_prime_repair_within_radius"
                category_counts[category] += 1
                cases.append(
                    {
                        "source_m": source_m,
                        "source_k": source_k,
                        "target_base_k": base_k,
                        "first_small_prime_tier_delta": None,
                        "first_witness_delta": None,
                        "window_blockers_at_first_tier": [],
                        "category": category,
                    }
                )
                continue

            tier_delta_counts[first_tier_delta] += 1
            tier_k = base_k + first_tier_delta
            x = target_m + 2 * tier_k
            seen: set[int] = set()
            blockers: list[dict[str, Any]] = []
            for position in range(target_m // 2 + 1, target_m + 1):
                term = tier_k + position
                for p in factorizer.factor(term):
                    if p <= target_m or p in seen:
                        continue
                    seen.add(p)
                    value = slack(target_m, tier_k, p)
                    if value >= 0:
                        continue
                    regime = "single_level" if p * p > x else "compensation_deficit"
                    blockers.append(
                        {
                            "prime": p,
                            "position": position,
                            "term": term,
                            "slack": value,
                            "regime": regime,
                        }
                    )
                    blocker_counts[regime] += 1

            tier_is_witness = bool(target[tier_k])
            if tier_is_witness:
                if blockers or first_witness_delta != first_tier_delta:
                    raise AssertionError("tier-one witness classification mismatch")
                category = "first_small_prime_tier_is_witness"
            else:
                if not blockers:
                    raise AssertionError("tier-one nonwitness had no window blocker")
                if first_witness_delta is None:
                    category = "window_blocked_no_witness_within_radius"
                else:
                    category = "window_blocked_before_later_witness"
                    witness_delay_counts[
                        deltas.index(first_witness_delta) - deltas.index(first_tier_delta)
                    ] += 1
            category_counts[category] += 1
            cases.append(
                {
                    "source_m": source_m,
                    "source_k": source_k,
                    "target_base_k": base_k,
                    "first_small_prime_tier_delta": first_tier_delta,
                    "first_witness_delta": first_witness_delta,
                    "window_blockers_at_first_tier": blockers,
                    "category": category,
                }
            )

    eligible = counts["failed_natural_shifts"] - counts["boundary_excluded"]
    if eligible != len(cases) or eligible != sum(category_counts.values()):
        raise AssertionError("CRT repair cover partition mismatch")

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_CRT_REPAIR_COVER",
        "claim_boundary": (
            "Every boundary-comparable failed natural shift is searched only within "
            "the declared radius. The small-prime tier is exact on each candidate's "
            "prime-power bands. Large blockers are exactly factored from that "
            "candidate's moving bad window. No global repair theorem is claimed."
        ),
        "parameters": {
            "min_m": min_m,
            "max_m": max_m,
            "max_k": max_k,
            "radius": radius,
            "delta_order": "0,-1,+1,-2,+2,...",
            "small_prime_bound": "target_m",
        },
        "counts": dict(sorted(counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "first_small_prime_tier_delta_counts": [
            {"delta": delta, "count": count}
            for delta, count in sorted(tier_delta_counts.items())
        ],
        "witness_order_delay_counts": [
            {"order_delay": delay, "count": count}
            for delay, count in sorted(witness_delay_counts.items())
        ],
        "window_blocker_counts": dict(sorted(blocker_counts.items())),
        "cases": cases,
        "input": {
            "atlas_logical_name": atlas_path.name,
            "atlas_sha256": hashlib.sha256(atlas_bytes).hexdigest(),
            "atlas_schema_version": atlas["schema_version"],
        },
        "implementation_sha256": {
            "crt_repair_cover.py": sha256_file(Path(__file__)),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "shift_repair_analysis.py": sha256_file(
                ROOT / "shift_repair_analysis.py"
            ),
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
        atlas_path=args.atlas,
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
                "category_counts": result["category_counts"],
                "window_blocker_counts": result["window_blocker_counts"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
