#!/usr/bin/env python3
"""Build an exact bounded near-witness atlas for Erdős Problem #389.

The output exhausts only the declared rectangle in (m, k).  Offset coverage is
an empirical finite screen, not a CRT construction and not evidence for the
unbounded conjecture.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import sys
import time
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import RatioSlackScanner, SPFFactorizer

ROOT = Path(__file__).resolve().parent
DEFAULT_KNOWN = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_OUTPUT = ROOT / "data" / "atlas_m1_50_k50000.json"
def _input_provenance(path: Path, content: bytes) -> dict[str, str]:
    resolved = path.resolve()
    try:
        logical_name = resolved.relative_to(ROOT).as_posix()
    except ValueError:
        logical_name = path.name
    source = (
        "Local transcription of OEIS A375071 b-file, retrieved 2026-08-25"
        if resolved == DEFAULT_KNOWN.resolve()
        else "Caller-supplied witness table; external source identity not asserted"
    )
    return {
        "logical_name": logical_name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "source": source,
    }


def base_digits(n: int, p: int) -> list[int]:
    if n == 0:
        return [0]
    digits: list[int] = []
    while n:
        digits.append(n % p)
        n //= p
    return digits


def carry_trace(a: int, b: int, p: int) -> dict[str, Any]:
    rows: list[dict[str, int]] = []
    incoming = 0
    position = 0
    while a or b or incoming:
        a_digit = a % p
        b_digit = b % p
        total = a_digit + b_digit + incoming
        outgoing = int(total >= p)
        rows.append(
            {
                "position": position,
                "a_digit": a_digit,
                "b_digit": b_digit,
                "carry_in": incoming,
                "carry_out": outgoing,
            }
        )
        incoming = outgoing
        a //= p
        b //= p
        position += 1
    return {
        "carry_positions": [row["position"] for row in rows if row["carry_out"]],
        "digits": rows,
    }


def _read_published(content: bytes) -> dict[int, int]:
    rows = list(csv.DictReader(io.StringIO(content.decode("utf-8"))))
    published: dict[int, int] = {}
    for row in rows:
        m = int(row["m"])
        original_n = int(row["original_n"])
        k = int(row["k"])
        if original_n != m + 1 or m in published:
            raise ValueError(f"invalid published witness row: {row}")
        published[m] = k
    return published


def _representative(m: int, k: int, p: int) -> dict[str, Any]:
    left = carry_trace(k, m, p)
    right = carry_trace(k, m + k, p)
    return {
        "prime": p,
        "m": m,
        "k": k,
        "base_p_digits_lsd_first": {
            "m": base_digits(m, p),
            "k": base_digits(k, p),
            "m_plus_k": base_digits(m + k, p),
        },
        "left_addition_k_plus_m": left,
        "right_addition_k_plus_m_plus_k": right,
    }


def _carry_positions(a: int, b: int, p: int) -> tuple[int, ...]:
    positions: list[int] = []
    incoming = 0
    position = 0
    while a or b or incoming:
        total = a % p + b % p + incoming
        incoming = int(total >= p)
        if incoming:
            positions.append(position)
        a //= p
        b //= p
        position += 1
    return tuple(positions)


def _near_pattern_key(
    p: int,
    deficit: int,
    left_positions: tuple[int, ...],
    right_positions: tuple[int, ...],
) -> tuple[Any, ...]:
    return p, deficit, left_positions, right_positions


def _shift_pattern_key(
    p: int,
    source_slack: int,
    added_valuation: int,
    removed_valuation: int,
    shifted_slack: int,
    left_positions: tuple[int, ...],
    right_positions: tuple[int, ...],
) -> tuple[Any, ...]:
    return (
        p,
        source_slack,
        added_valuation,
        removed_valuation,
        shifted_slack,
        left_positions,
        right_positions,
    )


def _record_pattern(
    groups: dict[tuple[Any, ...], dict[str, Any]],
    key: tuple[Any, ...],
    representative: tuple[int, int, int],
    representative_limit: int,
) -> None:
    if not key or representative[2] != key[0]:
        raise AssertionError("pattern key prime does not match representative prime")
    m, k, p = representative
    if (
        _carry_positions(k, m, p) != key[-2]
        or _carry_positions(k, m + k, p) != key[-1]
    ):
        raise AssertionError("pattern key carries do not match representative")
    group = groups.setdefault(key, {"count": 0, "representatives": []})
    group["count"] += 1
    if len(group["representatives"]) < representative_limit:
        group["representatives"].append(representative)


def _expand_representatives(
    coordinates: list[tuple[int, int, int]],
) -> list[dict[str, Any]]:
    return [_representative(m, k, p) for m, k, p in coordinates]


def _serialize_near_patterns(
    groups: dict[tuple[Any, ...], dict[str, Any]], limit: int
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for key, value in groups.items():
        p, deficit, left_positions, right_positions = key
        records.append(
            {
                "prime": p,
                "deficit": deficit,
                "left_carry_positions": list(left_positions),
                "right_carry_positions": list(right_positions),
                "count": value["count"],
                "_representatives": value["representatives"],
            }
        )
    records.sort(
        key=lambda row: (
            -row["count"],
            row["prime"],
            row["deficit"],
            row["left_carry_positions"],
            row["right_carry_positions"],
        )
    )
    reported = records[:limit]
    for record in reported:
        coordinates = record.pop("_representatives")
        record["representatives"] = _expand_representatives(coordinates)
    total_instances = sum(row["count"] for row in records)
    reported_instances = sum(row["count"] for row in reported)
    return {
        "key_fields": [
            "prime",
            "deficit",
            "left_carry_positions",
            "right_carry_positions",
        ],
        "distinct_patterns": len(records),
        "total_instances": total_instances,
        "reported_patterns": len(reported),
        "reported_instances": reported_instances,
        "reported_fraction": {
            "numerator": reported_instances,
            "denominator": total_instances,
        },
        "patterns": reported,
    }


def _serialize_shift_patterns(
    groups: dict[tuple[Any, ...], dict[str, Any]], limit: int
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for key, value in groups.items():
        (
            p,
            source_slack,
            added_valuation,
            removed_valuation,
            shifted_slack,
            left_positions,
            right_positions,
        ) = key
        records.append(
            {
                "prime": p,
                "source_slack": source_slack,
                "v_p_m_plus_one": added_valuation,
                "v_p_m_plus_two_k": removed_valuation,
                "shifted_slack": shifted_slack,
                "left_carry_positions": list(left_positions),
                "right_carry_positions": list(right_positions),
                "count": value["count"],
                "_representatives": value["representatives"],
            }
        )
    records.sort(
        key=lambda row: (
            -row["count"],
            row["prime"],
            row["source_slack"],
            row["v_p_m_plus_one"],
            row["v_p_m_plus_two_k"],
            row["shifted_slack"],
            row["left_carry_positions"],
            row["right_carry_positions"],
        )
    )
    reported = records[:limit]
    for record in reported:
        coordinates = record.pop("_representatives")
        record["representatives"] = _expand_representatives(coordinates)
    total_instances = sum(row["count"] for row in records)
    reported_instances = sum(row["count"] for row in reported)
    return {
        "key_fields": [
            "prime",
            "source_slack",
            "v_p_m_plus_one",
            "v_p_m_plus_two_k",
            "shifted_slack",
            "left_carry_positions",
            "right_carry_positions",
        ],
        "distinct_patterns": len(records),
        "total_instances": total_instances,
        "reported_patterns": len(reported),
        "reported_instances": reported_instances,
        "reported_fraction": {
            "numerator": reported_instances,
            "denominator": total_instances,
        },
        "patterns": reported,
    }


def _bounded_offset_screen(
    witness_bits: list[bytearray | None],
    min_m: int,
    max_m: int,
    max_k: int,
    radius: int,
) -> dict[str, Any]:
    deltas = list(range(-radius, radius + 1))
    coverage = Counter({delta: 0 for delta in deltas})
    mask_counts: Counter[int] = Counter()
    eligible_failures = 0
    boundary_failures = 0

    for m in range(min_m, max_m):
        source = witness_bits[m]
        target = witness_bits[m + 1]
        assert source is not None and target is not None
        for k in range(2, max_k + 1):
            if not source[k] or target[k - 1]:
                continue
            base = k - 1
            if base - radius < 1 or base + radius > max_k:
                boundary_failures += 1
                continue
            eligible_failures += 1
            mask = 0
            for bit, delta in enumerate(deltas):
                if target[base + delta]:
                    coverage[delta] += 1
                    mask |= 1 << bit
            mask_counts[mask] += 1

    unrepairable = mask_counts.pop(0, 0)
    remaining = dict(mask_counts)
    greedy_menu: list[dict[str, int]] = []
    while remaining:
        scores = [
            sum(count for mask, count in remaining.items() if mask & (1 << bit))
            for bit in range(len(deltas))
        ]
        best_bit = max(
            range(len(deltas)), key=lambda bit: (scores[bit], -abs(deltas[bit]))
        )
        newly_covered = scores[best_bit]
        if newly_covered == 0:
            break
        greedy_menu.append(
            {"delta": deltas[best_bit], "newly_covered": newly_covered}
        )
        remaining = {
            mask: count
            for mask, count in remaining.items()
            if not (mask & (1 << best_bit))
        }
    if remaining:
        raise AssertionError("greedy offset accounting left a nonzero repair mask")
    covered_by_menu = sum(item["newly_covered"] for item in greedy_menu)
    if covered_by_menu + unrepairable != eligible_failures:
        raise AssertionError("offset coverage partition does not match eligible failures")


    return {
        "status": "EXACT_FOR_DECLARED_BOUNDED_SCREEN_ONLY",
        "interpretation": (
            "Tests k' = k - 1 + delta. These are literal bounded offsets, not "
            "CRT-defined repairs and not a propagation theorem."
        ),
        "radius": radius,
        "eligible_natural_shift_failures": eligible_failures,
        "boundary_excluded_failures": boundary_failures,
        "coverage_by_offset": [
            {"delta": delta, "covered": coverage[delta]} for delta in deltas
        ],
        "greedy_offset_menu": greedy_menu,
        "uncovered_by_every_tested_offset": unrepairable,
        "covered_by_greedy_menu": covered_by_menu,
    }


def build_atlas(
    *,
    min_m: int,
    max_m: int,
    max_k: int,
    repair_radius: int,
    representative_limit: int,
    top_primes: int,
    pattern_limit: int,
    known_path: Path,
) -> dict[str, Any]:
    if not 0 <= min_m <= max_m:
        raise ValueError("require 0 <= min_m <= max_m")
    if min(max_k, representative_limit, top_primes, pattern_limit) < 1 or repair_radius < 0:
        raise ValueError(
            "max_k, representative_limit, top_primes, and pattern_limit must be "
            "positive; radius must be nonnegative"
        )

    known_path = known_path.resolve()
    known_bytes = known_path.read_bytes()
    published = _read_published(known_bytes)
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    witness_bits: list[bytearray | None] = [None] * (max_m + 1)
    rows: list[dict[str, Any]] = []
    near_patterns: dict[tuple[Any, ...], dict[str, Any]] = {}
    shift_patterns: dict[tuple[Any, ...], dict[str, Any]] = {}
    global_shift_prime_counts: Counter[int] = Counter()
    odd_large_prime_violations: list[dict[str, int]] = []

    for m in range(min_m, max_m + 1):
        scanner = RatioSlackScanner(m, factorizer)
        bits = bytearray(max_k + 1)
        witness_count = 0
        first: int | None = None
        near_count = 0
        near_obstructions: Counter[tuple[int, int]] = Counter()
        shift_eligible = 0
        shift_successes = 0
        shift_failures = 0
        shift_obstruction_count = 0
        shift_prime_counts: Counter[int] = Counter()
        m_plus_one_factors = factorizer.factor(m + 1)

        while True:
            k = scanner.k
            if scanner.is_witness:
                bits[k] = 1
                witness_count += 1
                if first is None:
                    first = k
                if m < max_m and k >= 2:
                    shift_eligible += 1
                    source_exponents = scanner.exponents
                    obstructions = []
                    for p, removed in factorizer.factor(m + 2 * k).items():
                        source_slack = source_exponents.get(p, 0)
                        added = m_plus_one_factors.get(p, 0)
                        shifted_slack = source_slack + added - removed
                        if shifted_slack >= 0:
                            continue
                        obstructions.append(p)
                        shift_obstruction_count += 1
                        shift_prime_counts[p] += 1
                        global_shift_prime_counts[p] += 1
                        left_positions = _carry_positions(k, m, p)
                        right_positions = _carry_positions(k, m + k, p)
                        key = _shift_pattern_key(
                            p,
                            source_slack,
                            added,
                            removed,
                            shifted_slack,
                            left_positions,
                            right_positions,
                        )
                        _record_pattern(
                            shift_patterns,
                            key,
                            (m, k, p),
                            representative_limit,
                        )
                        if m % 2 == 1 and p > m:
                            odd_large_prime_violations.append(
                                {"m": m, "k": k, "prime": p}
                            )
                    if obstructions:
                        shift_failures += 1
                    else:
                        shift_successes += 1
            elif scanner.negative_count == 1:
                near_count += 1
                p, exponent = next(iter(scanner.negative_exponents.items()))
                deficit = -exponent
                near_obstructions[(p, deficit)] += 1
                left_positions = _carry_positions(k, m, p)
                right_positions = _carry_positions(k, m + k, p)
                key = _near_pattern_key(
                    p, deficit, left_positions, right_positions
                )
                _record_pattern(
                    near_patterns, key, (m, k, p), representative_limit
                )

            if k == max_k:
                break
            scanner.advance()

        published_minimum = published.get(m)
        if published_minimum is not None:
            if published_minimum <= max_k and first != published_minimum:
                raise AssertionError(
                    f"bounded first witness m={m} is {first}, published {published_minimum}"
                )
            if published_minimum > max_k and first is not None:
                raise AssertionError(
                    f"found k={first} below published minimum {published_minimum} for m={m}"
                )

        sorted_near_obstructions = [
            {"prime": p, "deficit": deficit, "count": count}
            for (p, deficit), count in sorted(
                near_obstructions.items(),
                key=lambda item: (-item[1], item[0][0], item[0][1]),
            )
        ]
        reported_near_obstructions = sorted_near_obstructions[:top_primes]
        witness_bits[m] = bits
        rows.append(
            {
                "m": m,
                "original_n": m + 1,
                "witness_count": witness_count,
                "first_witness": first,
                "published_minimum": published_minimum,
                "published_minimum_reproduced": (
                    published_minimum is not None
                    and published_minimum <= max_k
                    and first == published_minimum
                ),
                "one_prime_near_witness_count": near_count,
                "one_prime_obstruction_classes": len(sorted_near_obstructions),
                "one_prime_obstruction_instances_reported": sum(
                    item["count"] for item in reported_near_obstructions
                ),
                "one_prime_obstruction_instances_omitted": near_count
                - sum(item["count"] for item in reported_near_obstructions),
                "one_prime_obstructions": reported_near_obstructions,
                "natural_shift": (
                    None
                    if m == max_m
                    else {
                        "eligible_source_witnesses": shift_eligible,
                        "successful": shift_successes,
                        "failed": shift_failures,
                        "obstruction_count": shift_obstruction_count,
                        "obstruction_prime_counts": [
                            {"prime": p, "count": count}
                            for p, count in sorted(
                                shift_prime_counts.items(),
                                key=lambda item: (-item[1], item[0]),
                            )
                        ],
                    }
                ),
            }
        )

    for m in range(min_m, max_m):
        source = witness_bits[m]
        target = witness_bits[m + 1]
        assert source is not None and target is not None
        row = rows[m - min_m]
        shift_summary = row["natural_shift"]
        assert shift_summary is not None
        direct_successes = sum(source[k] and target[k - 1] for k in range(2, max_k + 1))
        if direct_successes != shift_summary["successful"]:
            raise AssertionError(
                f"natural-shift identity mismatch at m={m}: "
                f"formula={shift_summary['successful']} direct={direct_successes}"
            )

    if odd_large_prime_violations:
        raise AssertionError(
            f"odd-m shift bound contradicted in bounded atlas: {odd_large_prime_violations[:3]}"
        )
    near_pattern_summary = _serialize_near_patterns(near_patterns, pattern_limit)
    shift_pattern_summary = _serialize_shift_patterns(shift_patterns, pattern_limit)
    expected_near_instances = sum(
        row["one_prime_near_witness_count"] for row in rows
    )
    expected_shift_instances = sum(
        (row["natural_shift"] or {}).get("obstruction_count", 0) for row in rows
    )
    if near_pattern_summary["total_instances"] != expected_near_instances:
        raise AssertionError("near-pattern grouping lost or duplicated scan events")
    if shift_pattern_summary["total_instances"] != expected_shift_instances:
        raise AssertionError("shift-pattern grouping lost or duplicated obstructions")


    return {
        "schema_version": 2,
        "status": "EXACT_BOUNDED_COMPUTATION",
        "claim_boundary": (
            "Every k in the declared rectangle was checked by an exact prime-exponent "
            "recurrence. Absence beyond max_k, general existence, and a finite CRT repair "
            "menu are not claimed."
        ),
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "parameters": {
            "min_m": min_m,
            "max_m": max_m,
            "max_k": max_k,
            "repair_radius": repair_radius,
            "representatives_per_pattern": representative_limit,
            "top_prime_classes_per_m": top_primes,
            "reported_pattern_limit": pattern_limit,
        },
        "classification_schema": {
            "base_sensitive": True,
            "prime_is_required_in_every_pattern_key": True,
        },
        "input": _input_provenance(known_path, known_bytes),
        "implementation_sha256": {
            "atlas.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "erdos389.py": hashlib.sha256(
                (ROOT / "erdos389.py").read_bytes()
            ).hexdigest(),
        },
        "m_rows": rows,
        "one_prime_near_witness_patterns": near_pattern_summary,
        "natural_shift_summary": {
            "identity": "s'_p = s_p + v_p(m+1) - v_p(m+2k)",
            "odd_m_large_prime_violations": 0,
            "obstruction_prime_counts": [
                {"prime": p, "count": count}
                for p, count in sorted(
                    global_shift_prime_counts.items(),
                    key=lambda item: (-item[1], item[0]),
                )
            ],
            "failure_patterns": shift_pattern_summary,
        },
        "bounded_offset_repair_screen": _bounded_offset_screen(
            witness_bits, min_m, max_m, max_k, repair_radius
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-m", type=int, default=1)
    parser.add_argument("--max-m", type=int, default=50)
    parser.add_argument("--max-k", type=int, default=50_000)
    parser.add_argument("--repair-radius", type=int, default=16)
    parser.add_argument("--representatives", type=int, default=2)
    parser.add_argument("--top-primes", type=int, default=32)
    parser.add_argument("--pattern-limit", type=int, default=256)
    parser.add_argument("--known", type=Path, default=DEFAULT_KNOWN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    atlas = build_atlas(
        min_m=args.min_m,
        max_m=args.max_m,
        max_k=args.max_k,
        repair_radius=args.repair_radius,
        representative_limit=args.representatives,
        top_primes=args.top_primes,
        pattern_limit=args.pattern_limit,
        known_path=args.known.resolve(),
    )
    atlas["resources"] = {
        "wall_seconds": time.perf_counter() - started_wall,
        "process_cpu_seconds": time.process_time() - started_cpu,
        "processes": 1,
    }
    atomic_write_json(args.output, atlas)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "m_rows": len(atlas["m_rows"]),
                "max_k": args.max_k,
                "wall_seconds": atlas["resources"]["wall_seconds"],
                "process_cpu_seconds": atlas["resources"]["process_cpu_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
