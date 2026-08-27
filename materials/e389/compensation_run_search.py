#!/usr/bin/env python3
"""Search exact runs of compensation-good bad-window integers.

Every rejected candidate carries a prime ``p > m`` whose exact local slack is
negative.  Every candidate whose entire bad window is compensation-good is
checked independently with complete Legendre and Kummer certificates, leaving
only the finite small-prime tier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import struct
import sys
import time
from collections import Counter, deque
from math import isqrt
from pathlib import Path
from typing import Any, Iterable

from artifact_io import atomic_write_json
from compensation_structure_analysis import term_compensation_formula
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    TrialFactorizer,
    is_witness_by_carries,
    primes_up_to,
    slack,
    witness_certificate,
)
from residue_control_analysis import upper_half_prefix_count

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "compensation_run_m27_k50001_h100000000.json"
_REJECTION_STRUCT = struct.Struct("<QIQQIi")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def segmented_compensation_obstructions(
    *,
    m: int,
    start: int,
    count: int,
    trial_primes: Iterable[int],
) -> tuple[list[int], bytearray, bytearray]:
    """Return the smallest exact obstructing ``p > m`` for consecutive terms.

    The supplied primes must be increasing and include every prime through
    ``isqrt(start + count - 1)``.  The three arrays contain the selected prime,
    its valuation, and the positive slack deficit.  Zero means that the term
    is compensation-good for every prime larger than ``m``.
    """
    _require_int("m", m, 1)
    _require_int("start", start, 1)
    _require_int("count", count, 1)
    primes = tuple(trial_primes)
    if any(left >= right for left, right in zip(primes, primes[1:])):
        raise ValueError("trial primes must be strictly increasing")

    remainders = list(range(start, start + count))
    obstruction_primes = [0] * count
    obstruction_exponents = bytearray(count)
    obstruction_deficits = bytearray(count)
    for p in primes:
        first = (-start) % p
        if first >= count:
            continue
        for index in range(first, count, p):
            value = remainders[index]
            exponent = 0
            prime_power = 1
            while value % p == 0:
                value //= p
                exponent += 1
                prime_power *= p
            remainders[index] = value
            if p <= m or obstruction_primes[index] or exponent == 0:
                continue
            term = start + index
            cofactor = term // prime_power
            predicted = upper_half_prefix_count(cofactor, p) - exponent
            if predicted < 0:
                obstruction_primes[index] = p
                obstruction_exponents[index] = exponent
                obstruction_deficits[index] = -predicted

    for index, residual in enumerate(remainders):
        if residual <= m or obstruction_primes[index] or residual == 1:
            continue
        term = start + index
        cofactor = term // residual
        predicted = upper_half_prefix_count(cofactor, residual) - 1
        if predicted < 0:
            obstruction_primes[index] = residual
            obstruction_exponents[index] = 1
            obstruction_deficits[index] = -predicted
    return obstruction_primes, obstruction_exponents, obstruction_deficits


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
    if is_witness_by_carries(m, k, factorizer) != certificate.is_witness:
        raise AssertionError("Legendre and Kummer witness checks disagree")
    obstructions = [
        {
            "prime": p,
            "left_valuation": certificate.left_valuations[p],
            "right_valuation": certificate.right_valuations[p],
            "slack": value,
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
    start_k: int,
    offset_limit: int,
    chunk_size: int,
) -> dict[str, Any]:
    """Classify every ``k`` in one interval by exact local compensation."""
    for name, value in (
        ("target_m", target_m),
        ("start_k", start_k),
        ("offset_limit", offset_limit),
        ("chunk_size", chunk_size),
    ):
        _require_int(name, value, 1)
    max_k = start_k + offset_limit - 1
    max_x = target_m + 2 * max_k
    if max_x > MAX_TRIAL_FACTORIZER_VALUE:
        raise ValueError("search exceeds the retained exact factorizer bound")

    first_position = target_m // 2 + 1
    bad_length = target_m - first_position + 1
    first_term = start_k + first_position
    integer_count = offset_limit + bad_length - 1
    last_term = first_term + integer_count - 1
    trial_prime_limit = isqrt(last_term)
    trial_primes = primes_up_to(trial_prime_limit)

    rejection_digest = hashlib.sha256()
    survivor_digest = hashlib.sha256()
    rejection_examples: list[dict[str, int]] = []
    survivor_records: list[dict[str, Any]] = []
    witness_records: list[dict[str, Any]] = []
    selected_prime_counts: Counter[int] = Counter()
    deficit_histogram: Counter[int] = Counter()
    counts = Counter()
    bad_queue: deque[tuple[int, int, int, int, int]] = deque()
    direct_factorizer: TrialFactorizer | None = None
    current_good_run = 0
    maximum_good_run = 0

    processed = 0
    while processed < integer_count:
        current_count = min(chunk_size, integer_count - processed)
        segment_start = first_term + processed
        primes, exponents, deficits = segmented_compensation_obstructions(
            m=target_m,
            start=segment_start,
            count=current_count,
            trial_primes=trial_primes,
        )
        for local_index, p in enumerate(primes):
            global_index = processed + local_index
            term = segment_start + local_index
            if p:
                exponent = exponents[local_index]
                deficit = deficits[local_index]
                bad_queue.append((global_index, term, p, exponent, deficit))
                counts["compensation_bad_terms"] += 1
                current_good_run = 0
            else:
                counts["compensation_good_terms"] += 1
                current_good_run += 1
                maximum_good_run = max(maximum_good_run, current_good_run)
            if global_index + 1 < bad_length:
                continue
            offset = global_index - bad_length + 1
            while bad_queue and bad_queue[0][0] < offset:
                bad_queue.popleft()
            if offset >= offset_limit:
                continue
            k = start_k + offset
            if bad_queue:
                bad_index, bad_term, p, exponent, deficit = bad_queue[0]
                position = first_position + bad_index - offset
                if bad_term != k + position:
                    raise AssertionError("rolling bad-window position drifted")
                predicted_slack = -deficit
                record = {
                    "offset": offset,
                    "position": position,
                    "term": bad_term,
                    "prime": p,
                    "exponent": exponent,
                    "slack": predicted_slack,
                }
                rejection_digest.update(
                    _REJECTION_STRUCT.pack(
                        offset,
                        position,
                        bad_term,
                        p,
                        exponent,
                        predicted_slack,
                    )
                )
                counts["large_prime_rejected_offsets"] += 1
                selected_prime_counts[p] += 1
                deficit_histogram[deficit] += 1
                if len(rejection_examples) < 40:
                    exact = term_compensation_formula(bad_term, p)
                    if exact["predicted_slack"] != predicted_slack:
                        raise AssertionError("representative local slack disagrees")
                    if slack(target_m, k, p) != predicted_slack:
                        raise AssertionError("representative exact slack disagrees")
                    counts["direct_rejection_checks"] += 1
                    rejection_examples.append(record)
                continue

            counts["large_prime_good_offsets"] += 1
            if direct_factorizer is None:
                direct_factorizer = TrialFactorizer(max_x)
            record = _certificate_record(target_m, k, direct_factorizer)
            if any(
                obstruction["prime"] > target_m
                for obstruction in record["obstructions"]
            ):
                raise AssertionError("compensation sieve missed a large prime")
            survivor_digest.update(_canonical_line(record))
            if len(survivor_records) < 1_000:
                survivor_records.append(record)
            if record["is_witness"]:
                counts["witnesses"] += 1
                witness_records.append(record)
            else:
                counts["small_prime_rejected_offsets"] += 1
        processed += current_count

    if (
        counts["large_prime_rejected_offsets"]
        + counts["large_prime_good_offsets"]
        != offset_limit
    ):
        raise AssertionError("compensation sieve did not partition all offsets")
    counts["maximum_consecutive_compensation_good_terms"] = maximum_good_run

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_COMPENSATION_RUN_SEARCH",
        "claim_boundary": (
            "Every k in the declared interval is either rejected by an exact "
            "negative local slack at a prime p>m or checked by complete Legendre "
            "and Kummer witness certificates. No claim is made outside the interval."
        ),
        "parameters": {
            "target_m": target_m,
            "start_k": start_k,
            "end_k": max_k,
            "offset_min": 0,
            "offset_max": offset_limit - 1,
            "offset_count": offset_limit,
            "chunk_size": chunk_size,
            "bad_window_first_position": first_position,
            "bad_window_length": bad_length,
            "factored_integer_interval": {
                "first": first_term,
                "last": last_term,
            },
            "trial_prime_limit": trial_prime_limit,
        },
        "counts": dict(sorted(counts.items())),
        "large_prime_rejections": {
            "certificate_sha256": rejection_digest.hexdigest(),
            "certificate_encoding": "little-endian <QIQQIi records in offset order",
            "record_fields": [
                "offset",
                "position",
                "term",
                "prime",
                "exponent",
                "slack",
            ],
            "selected_prime_counts": [
                {"prime": p, "count": count}
                for p, count in sorted(
                    selected_prime_counts.items(),
                    key=lambda item: (-item[1], item[0]),
                )[:100]
            ],
            "deficit_histogram": [
                {"deficit": deficit, "offsets": count}
                for deficit, count in sorted(deficit_histogram.items())
            ],
            "representatives": rejection_examples,
        },
        "survivors": {
            "certificate_sha256": survivor_digest.hexdigest(),
            "reported_record_limit": 1_000,
            "reported_records": survivor_records,
            "witness_records": witness_records,
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "compensation_run_search.py": sha256_file(Path(__file__)),
            "compensation_structure_analysis.py": sha256_file(
                ROOT / "compensation_structure_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--start-k", type=int, default=50_001)
    parser.add_argument("--offset-limit", type=int, default=100_000_000)
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
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
                "process_cpu_seconds": result["resources"][
                    "process_cpu_seconds"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
