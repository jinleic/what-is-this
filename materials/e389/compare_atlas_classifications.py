#!/usr/bin/env python3
"""Recompute old and prime-sensitive atlas pattern classes and emit an exact diff.

This scanner intentionally does not import atlas.py's pattern-key or serialization
helpers.  It shares only the exact ratio-slack scanner from erdos389.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from artifact_io import atomic_text_writer, atomic_write_json
from erdos389 import RatioSlackScanner, SPFFactorizer

ROOT = Path(__file__).resolve().parent
DEFAULT_PRE = ROOT / "data" / "atlas_m1_50_k50000_pre_prime_key_fix.json"
DEFAULT_POST = ROOT / "data" / "atlas_m1_50_k50000.json"
DEFAULT_SUMMARY = ROOT / "data" / "atlas_prime_key_classification_diff.json"
DEFAULT_CASES = ROOT / "data" / "atlas_prime_key_changed_cases.jsonl"

PatternKey = tuple[Any, ...]
Event = tuple[int, int, int, PatternKey, PatternKey]


def carry_positions(a: int, b: int, p: int) -> tuple[int, ...]:
    positions: list[int] = []
    incoming = 0
    position = 0
    while a or b or incoming:
        incoming = int(a % p + b % p + incoming >= p)
        if incoming:
            positions.append(position)
        a //= p
        b //= p
        position += 1
    return tuple(positions)


def near_key_object(key: PatternKey, *, prime_sensitive: bool) -> dict[str, Any]:
    if prime_sensitive:
        p, deficit, left, right = key
    else:
        deficit, left, right = key
        p = None
    result = {
        "deficit": deficit,
        "left_carry_positions": list(left),
        "right_carry_positions": list(right),
    }
    if p is not None:
        result["prime"] = p
    return result


def shift_key_object(key: PatternKey, *, prime_sensitive: bool) -> dict[str, Any]:
    if prime_sensitive:
        p, source, added, removed, shifted, left, right = key
    else:
        source, added, removed, shifted, left, right = key
        p = None
    result = {
        "source_slack": source,
        "v_p_m_plus_one": added,
        "v_p_m_plus_two_k": removed,
        "shifted_slack": shifted,
        "left_carry_positions": list(left),
        "right_carry_positions": list(right),
    }
    if p is not None:
        result["prime"] = p
    return result


def pattern_id(kind: str, pattern: dict[str, Any]) -> str:
    payload = json.dumps(
        {"kind": kind, "pattern": pattern}, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def scan_events(min_m: int, max_m: int, max_k: int) -> dict[str, list[Event]]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    events: dict[str, list[Event]] = {"near": [], "shift": []}

    for m in range(min_m, max_m + 1):
        scanner = RatioSlackScanner(m, factorizer)
        m_plus_one_factors = factorizer.factor(m + 1)
        while True:
            k = scanner.k
            if scanner.is_witness and m < max_m and k >= 2:
                source_exponents = scanner.exponents
                for p, removed in sorted(factorizer.factor(m + 2 * k).items()):
                    source = source_exponents.get(p, 0)
                    added = m_plus_one_factors.get(p, 0)
                    shifted = source + added - removed
                    if shifted >= 0:
                        continue
                    left = carry_positions(k, m, p)
                    right = carry_positions(k, m + k, p)
                    old_key = (source, added, removed, shifted, left, right)
                    new_key = (p, source, added, removed, shifted, left, right)
                    events["shift"].append((m, k, p, old_key, new_key))
            elif scanner.negative_count == 1:
                p, exponent = next(iter(scanner.negative_exponents.items()))
                deficit = -exponent
                left = carry_positions(k, m, p)
                right = carry_positions(k, m + k, p)
                old_key = (deficit, left, right)
                new_key = (p, deficit, left, right)
                events["near"].append((m, k, p, old_key, new_key))

            if k == max_k:
                break
            scanner.advance()

    return events


def classify_kind(kind: str, events: list[Event]) -> tuple[dict[str, Any], set[PatternKey]]:
    old_counts: Counter[PatternKey] = Counter()
    new_counts: Counter[PatternKey] = Counter()
    old_to_new: dict[PatternKey, Counter[PatternKey]] = defaultdict(Counter)
    old_primes: dict[PatternKey, set[int]] = defaultdict(set)
    new_representatives: dict[PatternKey, tuple[int, int, int]] = {}
    prime_event_counts: Counter[int] = Counter()
    prime_patterns: dict[int, set[PatternKey]] = defaultdict(set)
    regime_event_counts: Counter[str] = Counter()

    for m, k, p, old_key, new_key in events:
        old_counts[old_key] += 1
        new_counts[new_key] += 1
        old_to_new[old_key][new_key] += 1
        old_primes[old_key].add(p)
        new_representatives.setdefault(new_key, (m, k, p))
        prime_event_counts[p] += 1
        prime_patterns[p].add(new_key)
        if kind == "near":
            regime = "p_le_m" if p <= m else "p_gt_m"
        else:
            parity = "odd_m" if m % 2 else "even_m"
            size = "p_le_m" if p <= m else "p_gt_m"
            regime = f"{parity}_{size}"
        regime_event_counts[regime] += 1

    split_keys = {key for key, primes in old_primes.items() if len(primes) > 1}
    split_groups: list[dict[str, Any]] = []
    old_object = near_key_object if kind == "near" else shift_key_object

    for old_key in sorted(split_keys, key=lambda key: pattern_id(kind, old_object(key, prime_sensitive=False))):
        old_pattern = old_object(old_key, prime_sensitive=False)
        children = []
        for new_key, count in sorted(
            old_to_new[old_key].items(),
            key=lambda item: pattern_id(kind, old_object(item[0], prime_sensitive=True)),
        ):
            new_pattern = old_object(new_key, prime_sensitive=True)
            children.append(
                {
                    "new_pattern_id": pattern_id(kind, new_pattern),
                    "count": count,
                    "pattern": new_pattern,
                }
            )
        split_groups.append(
            {
                "old_pattern_id": pattern_id(kind, old_pattern),
                "old_count": old_counts[old_key],
                "old_pattern": old_pattern,
                "prime_count": len(old_primes[old_key]),
                "new_patterns": children,
            }
        )

    top_new_patterns = []
    for new_key, count in sorted(
        new_counts.items(),
        key=lambda item: (
            -item[1],
            pattern_id(kind, old_object(item[0], prime_sensitive=True)),
        ),
    )[:50]:
        pattern = old_object(new_key, prime_sensitive=True)
        m, k, p = new_representatives[new_key]
        top_new_patterns.append(
            {
                "new_pattern_id": pattern_id(kind, pattern),
                "count": count,
                "pattern": pattern,
                "representative": {"m": m, "k": k, "prime": p},
            }
        )
    prime_distribution = [
        {
            "prime": p,
            "event_count": prime_event_counts[p],
            "distinct_pattern_count": len(prime_patterns[p]),
        }
        for p in sorted(prime_event_counts, key=lambda p: (-prime_event_counts[p], p))
    ]
    def describe_patterns(keys: list[PatternKey]) -> list[dict[str, Any]]:
        records = []
        for key in sorted(
            keys,
            key=lambda item: (
                -new_counts[item],
                pattern_id(kind, old_object(item, prime_sensitive=True)),
            ),
        )[:20]:
            pattern = old_object(key, prime_sensitive=True)
            m, k, p = new_representatives[key]
            records.append(
                {
                    "new_pattern_id": pattern_id(kind, pattern),
                    "count": new_counts[key],
                    "pattern": pattern,
                    "representative": {"m": m, "k": k, "prime": p},
                }
            )
        return records

    if kind == "near":
        extremal_values = {
            "maximum_deficit": max((key[1] for key in new_counts), default=0),
            "maximum_left_carry_count": max(
                (len(key[2]) for key in new_counts), default=0
            ),
            "maximum_right_carry_count": max(
                (len(key[3]) for key in new_counts), default=0
            ),
        }
        extremal_patterns = {
            "maximum_deficit": describe_patterns(
                [key for key in new_counts if key[1] == extremal_values["maximum_deficit"]]
            ),
            "maximum_left_carry_count": describe_patterns(
                [
                    key
                    for key in new_counts
                    if len(key[2]) == extremal_values["maximum_left_carry_count"]
                ]
            ),
            "maximum_right_carry_count": describe_patterns(
                [
                    key
                    for key in new_counts
                    if len(key[3]) == extremal_values["maximum_right_carry_count"]
                ]
            ),
        }
    else:
        extremal_values = {
            "minimum_shifted_slack": min(
                (key[4] for key in new_counts), default=0
            ),
            "maximum_source_slack": max(
                (key[1] for key in new_counts), default=0
            ),
            "maximum_removed_valuation": max(
                (key[3] for key in new_counts), default=0
            ),
        }
        extremal_patterns = {
            "minimum_shifted_slack": describe_patterns(
                [
                    key
                    for key in new_counts
                    if key[4] == extremal_values["minimum_shifted_slack"]
                ]
            ),
            "maximum_source_slack": describe_patterns(
                [
                    key
                    for key in new_counts
                    if key[1] == extremal_values["maximum_source_slack"]
                ]
            ),
            "maximum_removed_valuation": describe_patterns(
                [
                    key
                    for key in new_counts
                    if key[3] == extremal_values["maximum_removed_valuation"]
                ]
            ),
        }

    changed_case_count = sum(old_counts[key] for key in split_keys)
    return (
        {
            "event_count": len(events),
            "old_distinct_patterns": len(old_counts),
            "new_distinct_patterns": len(new_counts),
            "split_old_patterns": len(split_keys),
            "changed_case_count": changed_case_count,
            "unchanged_case_count": len(events) - changed_case_count,
            "split_groups": split_groups,
            "top_new_patterns": top_new_patterns,
            "prime_distribution": prime_distribution,
            "extremal_values": extremal_values,
            "extremal_patterns": extremal_patterns,
            "regime_event_counts": dict(sorted(regime_event_counts.items())),
        },
        split_keys,
    )


def recompute_reported_projection(
    kind: str,
    events: list[Event],
    *,
    prime_sensitive: bool,
    limit: int,
) -> dict[str, Any]:
    key_index = 4 if prime_sensitive else 3
    counts: Counter[PatternKey] = Counter(event[key_index] for event in events)
    if kind == "near":
        object_for = near_key_object
        key_fields = (
            ["prime"] if prime_sensitive else []
        ) + ["deficit", "left_carry_positions", "right_carry_positions"]
        if prime_sensitive:
            sort_key = lambda item: (
                -item[1],
                item[0][0],
                item[0][1],
                item[0][2],
                item[0][3],
            )
        else:
            sort_key = lambda item: (-item[1], item[0][0], item[0][1], item[0][2])
    else:
        object_for = shift_key_object
        key_fields = (
            ["prime"] if prime_sensitive else []
        ) + [
            "source_slack",
            "v_p_m_plus_one",
            "v_p_m_plus_two_k",
            "shifted_slack",
            "left_carry_positions",
            "right_carry_positions",
        ]
        if prime_sensitive:
            sort_key = lambda item: (
                -item[1],
                item[0][0],
                item[0][1],
                item[0][2],
                item[0][3],
                item[0][4],
                item[0][5],
                item[0][6],
            )
        else:
            # Schema 1 sorted only by shifted slack and carry positions; stable
            # sorting preserves first scan occurrence for remaining ties.
            sort_key = lambda item: (
                -item[1],
                item[0][3],
                item[0][4],
                item[0][5],
            )
    ordered = sorted(counts.items(), key=sort_key)
    reported = ordered[:limit]
    total_instances = sum(counts.values())
    reported_instances = sum(count for _, count in reported)
    return {
        "key_fields": key_fields,
        "distinct_patterns": len(counts),
        "total_instances": total_instances,
        "reported_patterns": len(reported),
        "reported_instances": reported_instances,
        "reported_fraction": {
            "numerator": reported_instances,
            "denominator": total_instances,
        },
        "patterns": [
            {
                **object_for(key, prime_sensitive=prime_sensitive),
                "count": count,
            }
            for key, count in reported
        ],
    }


def assert_artifact_projection(
    kind: str,
    artifact: dict[str, Any],
    expected: dict[str, Any],
    *,
    require_key_fields: bool,
) -> None:
    for field in (
        "distinct_patterns",
        "total_instances",
        "reported_patterns",
        "reported_instances",
    ):
        if artifact[field] != expected[field]:
            raise AssertionError(f"{kind} artifact field mismatch: {field}")
    if require_key_fields:
        if artifact.get("key_fields") != expected["key_fields"]:
            raise AssertionError(f"{kind} artifact key_fields mismatch")
        if artifact.get("reported_fraction") != expected["reported_fraction"]:
            raise AssertionError(f"{kind} artifact reported_fraction mismatch")
    fields = expected["key_fields"]
    normalized = [
        {**{field: row[field] for field in fields}, "count": row["count"]}
        for row in artifact["patterns"]
    ]
    if normalized != expected["patterns"]:
        raise AssertionError(f"{kind} reported pattern projection mismatch")


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_changed_cases(
    path: Path,
    events_by_kind: dict[str, list[Event]],
    split_keys_by_kind: dict[str, set[PatternKey]],
    parameters: dict[str, int],
) -> tuple[int, str]:
    digest = hashlib.sha256()
    count = 0
    with atomic_text_writer(path) as stream:
        header = json.dumps(
            {
                "type": "metadata",
                "schema_version": 1,
                "status": "EXACT_RECOMPUTATION",
                "parameters": parameters,
                "case_definition": (
                    "A scan event whose old prime-omitting equivalence class contains "
                    "more than one prime."
                ),
            },
            sort_keys=True,
            separators=(",", ":"),
        ) + "\n"
        stream.write(header)
        digest.update(header.encode("utf-8"))

        object_for = {"near": near_key_object, "shift": shift_key_object}
        for kind in ("near", "shift"):
            for m, k, p, old_key, new_key in events_by_kind[kind]:
                if old_key not in split_keys_by_kind[kind]:
                    continue
                old_pattern = object_for[kind](old_key, prime_sensitive=False)
                new_pattern = object_for[kind](new_key, prime_sensitive=True)
                line = json.dumps(
                    {
                        "type": "case",
                        "kind": kind,
                        "m": m,
                        "k": k,
                        "prime": p,
                        "old_pattern_id": pattern_id(kind, old_pattern),
                        "new_pattern_id": pattern_id(kind, new_pattern),
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                ) + "\n"
                stream.write(line)
                digest.update(line.encode("utf-8"))
                count += 1
    return count, digest.hexdigest()


def build_diff(
    pre_path: Path,
    post_path: Path,
    cases_path: Path,
    min_m: int,
    max_m: int,
    max_k: int,
) -> dict[str, Any]:
    pre_path = pre_path.resolve()
    post_path = post_path.resolve()
    cases_path = cases_path.resolve()
    if len({pre_path, post_path, cases_path}) != 3:
        raise ValueError("pre, post, and changed-cases paths must be distinct")
    pre_bytes = pre_path.read_bytes()
    post_bytes = post_path.read_bytes()
    pre = json.loads(pre_bytes)
    post = json.loads(post_bytes)
    if pre["schema_version"] != 1 or post["schema_version"] != 2:
        raise AssertionError("expected pre-fix schema 1 and repaired schema 2")
    if pre["m_rows"] != post["m_rows"]:
        raise AssertionError("prime-key repair changed per-m witness classifications")
    if pre["bounded_offset_repair_screen"] != post["bounded_offset_repair_screen"]:
        raise AssertionError("prime-key repair changed offset classifications")

    events = scan_events(min_m, max_m, max_k)
    summaries: dict[str, dict[str, Any]] = {}
    split_keys: dict[str, set[PatternKey]] = {}
    for kind in ("near", "shift"):
        summaries[kind], split_keys[kind] = classify_kind(kind, events[kind])

    pre_near = pre["one_prime_near_witness_patterns"]
    post_near = post["one_prime_near_witness_patterns"]
    pre_shift = pre["natural_shift_summary"]["failure_patterns"]
    post_shift = post["natural_shift_summary"]["failure_patterns"]
    pattern_limit = post["parameters"]["reported_pattern_limit"]
    if pre["parameters"]["reported_pattern_limit"] != pattern_limit:
        raise AssertionError("pre/post pattern limits differ")
    expected = {
        "near": (pre_near, post_near),
        "shift": (pre_shift, post_shift),
    }
    for kind, (old_artifact, new_artifact) in expected.items():
        old_projection = recompute_reported_projection(
            kind, events[kind], prime_sensitive=False, limit=pattern_limit
        )
        new_projection = recompute_reported_projection(
            kind, events[kind], prime_sensitive=True, limit=pattern_limit
        )
        assert_artifact_projection(
            f"old {kind}", old_artifact, old_projection, require_key_fields=False
        )
        assert_artifact_projection(
            f"new {kind}", new_artifact, new_projection, require_key_fields=True
        )
        summary = summaries[kind]
        if summary["old_distinct_patterns"] != old_projection["distinct_patterns"]:
            raise AssertionError(f"{kind} old grouping was not independently reproduced")
        if summary["new_distinct_patterns"] != new_projection["distinct_patterns"]:
            raise AssertionError(f"{kind} repaired grouping was not independently reproduced")

    parameters = {"min_m": min_m, "max_m": max_m, "max_k": max_k}
    changed_count, cases_sha256 = write_changed_cases(
        cases_path, events, split_keys, parameters
    )
    expected_changed = sum(summary["changed_case_count"] for summary in summaries.values())
    if changed_count != expected_changed:
        raise AssertionError("changed-case ledger count mismatch")

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_RECOMPUTATION",
        "claim_boundary": (
            "This diff covers every near-witness and natural-shift obstruction event "
            "in the declared rectangle. It compares only pattern equivalence classes; "
            "witness truth values and offset results are asserted unchanged."
        ),
        "parameters": parameters,
        "artifacts": {
            "pre_fix": {
                "logical_name": pre_path.name,
                "sha256": sha256_bytes(pre_bytes),
                "schema_version": pre["schema_version"],
            },
            "repaired": {
                "logical_name": post_path.name,
                "sha256": sha256_bytes(post_bytes),
                "schema_version": post["schema_version"],
            },
            "changed_cases": {
                "logical_name": cases_path.name,
                "sha256": cases_sha256,
                "case_count": changed_count,
            },
        },
        "unchanged": {
            "m_rows": True,
            "bounded_offset_repair_screen": True,
            "near_event_count": len(events["near"]),
            "shift_obstruction_event_count": len(events["shift"]),
        },
        "classifications": summaries,
        "implementation_sha256": {
            "compare_atlas_classifications.py": artifact_sha256(Path(__file__)),
            "atlas.py": artifact_sha256(ROOT / "atlas.py"),
            "erdos389.py": artifact_sha256(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre", type=Path, default=DEFAULT_PRE)
    parser.add_argument("--post", type=Path, default=DEFAULT_POST)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--min-m", type=int, default=1)
    parser.add_argument("--max-m", type=int, default=50)
    parser.add_argument("--max-k", type=int, default=50_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    resolved_paths = {
        args.pre.resolve(),
        args.post.resolve(),
        args.cases.resolve(),
        args.summary.resolve(),
    }
    if len(resolved_paths) != 4:
        raise ValueError("pre, post, cases, and summary paths must be distinct")
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = build_diff(
        args.pre.resolve(),
        args.post.resolve(),
        args.cases.resolve(),
        args.min_m,
        args.max_m,
        args.max_k,
    )
    result["resources"] = {
        "wall_seconds": time.perf_counter() - started_wall,
        "process_cpu_seconds": time.process_time() - started_cpu,
        "processes": 1,
    }
    atomic_write_json(args.summary, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "summary": str(args.summary),
                "cases": str(args.cases),
                "changed_cases": result["artifacts"]["changed_cases"]["case_count"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
