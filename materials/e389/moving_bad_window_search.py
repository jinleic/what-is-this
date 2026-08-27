#!/usr/bin/env python3
"""Exact moving bad-window sieve near a known natural-shift candidate.

Every skipped offset is certified by a prime larger than ``sqrt(m+2k)`` in its
bad window, hence by the proved single-level obstruction.  Only offsets that
survive this necessary sieve receive a full exact witness certificate.
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
from collections import Counter, deque
from math import isqrt
from pathlib import Path
from typing import Any, Iterable

from artifact_io import atomic_write_json
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    TrialFactorizer,
    is_witness_by_carries,
    natural_shift_obstructions,
    primes_up_to,
    witness_certificate,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_KNOWN = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_OUTPUT = ROOT / "data" / "moving_bad_window_m27_h1000000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_known(content: bytes) -> dict[int, int]:
    rows = csv.DictReader(io.StringIO(content.decode("utf-8")))
    return {int(row["m"]): int(row["k"]) for row in rows}


def segmented_large_cofactors(
    start: int,
    count: int,
    trial_primes: Iterable[int],
) -> list[int]:
    """Remove every supplied prime factor from a consecutive integer segment."""
    if start < 1 or count < 1:
        raise ValueError("require start>=1 and count>=1")
    remainders = list(range(start, start + count))
    end = start + count
    for p in trial_primes:
        first = (-start) % p
        if first >= count:
            continue
        for index in range(first, count, p):
            value = remainders[index]
            while value % p == 0:
                value //= p
            remainders[index] = value
    return remainders


def _canonical_line(record: dict[str, Any]) -> bytes:
    return (
        json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def _certificate_record(
    m: int,
    k: int,
    factorizer: TrialFactorizer,
) -> dict[str, Any]:
    certificate = witness_certificate(m, k, factorizer)
    carry_result = is_witness_by_carries(m, k, factorizer)
    if carry_result != certificate.is_witness:
        raise AssertionError("Legendre and Kummer witness checks disagree")
    obstructions = [
        {
            "prime": p,
            "left_valuation": certificate.left_valuations[p],
            "right_valuation": certificate.right_valuations[p],
            "slack": value,
            "regime": "small" if p <= m else "large",
        }
        for p, value in certificate.obstructions.items()
    ]
    return {
        "k": k,
        "is_witness": certificate.is_witness,
        "obstructions": obstructions,
        "tight_prime_count": len(certificate.tight_primes),
        "tight_prime_representatives": list(certificate.tight_primes[:20]),
    }


def search(
    *,
    target_m: int,
    start_k: int | None,
    offset_limit: int,
    chunk_size: int,
    known_path: Path,
) -> dict[str, Any]:
    if target_m < 1 or offset_limit < 1 or chunk_size < 1:
        raise ValueError("target_m, offset_limit, and chunk_size must be positive")

    known_path = known_path.resolve()
    known_bytes = known_path.read_bytes()
    known = _read_known(known_bytes)
    source_m = target_m - 1
    source_k = known.get(source_m)
    expected_start = source_k - 1 if source_k is not None else None
    if start_k is None:
        if expected_start is None:
            raise ValueError(
                "start_k is required when no published source row exists"
            )
        start_k = expected_start
    if start_k < 1:
        raise ValueError("start_k must be positive")
    max_k = start_k + offset_limit - 1
    max_x = target_m + 2 * max_k
    if max_x > MAX_TRIAL_FACTORIZER_VALUE:
        raise ValueError("search exceeds the exact retained factorizer bound")
    start_kind = (
        "natural_shift_of_published_source"
        if expected_start == start_k
        else "caller_supplied"
    )

    bad_first_position = target_m // 2 + 1
    bad_length = target_m - bad_first_position + 1
    first_term = start_k + bad_first_position
    integer_count = offset_limit + bad_length - 1
    last_term = first_term + integer_count - 1
    min_x = target_m + 2 * start_k
    if last_term >= min_x:
        raise ValueError(
            "offset interval is too wide for the residual-prime sieve; "
            "require the last bad-window term to be below the first x=m+2k"
        )
    trial_primes = primes_up_to(isqrt(last_term))

    fatal_digest = hashlib.sha256()
    survivor_digest = hashlib.sha256()
    fatal_examples: list[dict[str, int]] = []
    survivor_records: list[dict[str, Any]] = []
    witness_records: list[dict[str, Any]] = []
    counts = Counter()
    fatal_count_histogram: Counter[int] = Counter()
    selected_fatal_prime_counts: Counter[int] = Counter()
    window: deque[tuple[int, int]] = deque(maxlen=bad_length)
    direct_factorizer: TrialFactorizer | None = None

    processed = 0
    while processed < integer_count:
        current_count = min(chunk_size, integer_count - processed)
        segment_start = first_term + processed
        residuals = segmented_large_cofactors(
            segment_start, current_count, trial_primes
        )
        for local_index, residual in enumerate(residuals):
            term = segment_start + local_index
            window.append((term, residual))
            global_index = processed + local_index
            if global_index + 1 < bad_length:
                continue
            offset = global_index - bad_length + 1
            if offset >= offset_limit:
                continue
            k = start_k + offset
            x = target_m + 2 * k
            fatals = [
                (position, term_value, prime)
                for position, (term_value, prime) in enumerate(
                    window, start=bad_first_position
                )
                if prime > target_m and prime * prime > x
            ]
            fatal_count_histogram[len(fatals)] += 1
            fatal = fatals[0] if fatals else None
            if fatal is not None:
                position, term_value, prime = fatal
                if term_value % prime or not (
                    prime > target_m and prime * prime > x
                ):
                    raise AssertionError("invalid single-level obstruction certificate")
                record = {
                    "offset": offset,
                    "position": position,
                    "term": term_value,
                    "prime": prime,
                }
                fatal_digest.update(_canonical_line(record))
                counts["single_level_rejected_offsets"] += 1
                selected_fatal_prime_counts[prime] += 1
                if len(fatal_examples) < 40:
                    fatal_examples.append(record)
                continue

            counts["single_level_sieve_survivors"] += 1
            if direct_factorizer is None:
                direct_factorizer = TrialFactorizer(max_x)
            record = {"offset": offset, **_certificate_record(target_m, k, direct_factorizer)}
            survivor_digest.update(_canonical_line(record))
            if len(survivor_records) < 1_000:
                survivor_records.append(record)
            if record["is_witness"]:
                counts["witnesses"] += 1
                witness_records.append(record)
            else:
                counts["survivor_nonwitnesses"] += 1
                regimes = {item["regime"] for item in record["obstructions"]}
                if regimes == {"small"}:
                    counts["survivors_with_only_small_prime_obstructions"] += 1
                elif regimes == {"large"}:
                    counts["survivors_with_only_compensated_range_obstructions"] += 1
                else:
                    counts["survivors_with_mixed_obstructions"] += 1
        processed += current_count

    if counts["single_level_rejected_offsets"] + counts[
        "single_level_sieve_survivors"
    ] != offset_limit:
        raise AssertionError("moving-window sieve did not partition all offsets")

    source_record: dict[str, Any] | None = None
    if source_k is not None and source_k >= 2:
        source_factorizer = TrialFactorizer(max(source_m + 2 * source_k, max_x))
        source_certificate = witness_certificate(source_m, source_k, source_factorizer)
        if not source_certificate.is_witness:
            raise AssertionError("published source is not an exact witness")
        shift_obstructions = natural_shift_obstructions(
            source_m, source_k, source_factorizer
        )
        source_record = {
            "m": source_m,
            "k": source_k,
            "is_witness": True,
            "natural_shift_target_k": source_k - 1,
            "natural_shift_obstructions": [
                {
                    "prime": item.prime,
                    "source_slack": item.source_slack,
                    "removed_valuation": item.m_plus_two_k_valuation,
                    "shifted_slack": item.shifted_slack,
                }
                for item in shift_obstructions
            ],
        }

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_MOVING_BAD_WINDOW_SEARCH",
        "claim_boundary": (
            "Every offset in the declared forward interval is either rejected by an "
            "exact single-level large-prime certificate or checked by complete exact "
            "Legendre and Kummer witness tests. No claim is made outside the interval."
        ),
        "parameters": {
            "target_m": target_m,
            "start_k": start_k,
            "offset_min": 0,
            "offset_max": offset_limit - 1,
            "offset_count": offset_limit,
            "chunk_size": chunk_size,
            "bad_window_first_position": bad_first_position,
            "bad_window_length": bad_length,
            "factored_integer_interval": {"first": first_term, "last": last_term},
            "trial_prime_limit": isqrt(last_term),
            "residual_sieve_invariant": "last_bad_window_term < first_x",
        },
        "starting_point": {
            "kind": start_kind,
            "published_source": source_record,
        },
        "counts": dict(sorted(counts.items())),
        "single_level_rejections": {
            "certificate_sha256": fatal_digest.hexdigest(),
            "fatal_obstructions_per_window": [
                {"count": count, "windows": windows}
                for count, windows in sorted(fatal_count_histogram.items())
            ],
            "selected_certificate_prime_counts": [
                {"prime": prime, "count": count}
                for prime, count in sorted(
                    selected_fatal_prime_counts.items(),
                    key=lambda item: (-item[1], item[0]),
                )[:100]
            ],
            "representatives": fatal_examples,
        },
        "sieve_survivors": {
            "certificate_sha256": survivor_digest.hexdigest(),
            "reported_record_limit": 1_000,
            "reported_records": survivor_records,
            "witness_records": witness_records,
        },
        "input": {
            "known_logical_name": known_path.name,
            "known_sha256": hashlib.sha256(known_bytes).hexdigest(),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "moving_bad_window_search.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--start-k", type=int)
    parser.add_argument("--offset-limit", type=int, default=1_000_000)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--known", type=Path, default=DEFAULT_KNOWN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = search(
        target_m=args.target_m,
        start_k=args.start_k,
        offset_limit=args.offset_limit,
        chunk_size=args.chunk_size,
        known_path=args.known,
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
                "process_cpu_seconds": result["resources"]["process_cpu_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
