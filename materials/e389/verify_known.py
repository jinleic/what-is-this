#!/usr/bin/env python3
"""Verify a 27-row witness table; the default is local OEIS A375071 data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import random
import time
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import (
    SPFFactorizer,
    TrialFactorizer,
    carry_count,
    first_witness,
    is_witness,
    is_witness_by_carries,
    is_witness_direct,
    natural_shift_obstructions,
    witness_certificate,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_OUTPUT = ROOT / "data" / "known_witness_verification.json"
def input_provenance(path: Path, content: bytes) -> dict[str, str]:
    resolved = path.resolve()
    try:
        logical_name = resolved.relative_to(ROOT).as_posix()
    except ValueError:
        logical_name = path.name
    source = (
        "Local transcription of OEIS A375071 b-file, retrieved 2026-08-25"
        if resolved == DEFAULT_INPUT.resolve()
        else "Caller-supplied witness table; external source identity not asserted"
    )
    return {
        "logical_name": logical_name,
        "sha256": hashlib.sha256(content).hexdigest(),
        "source": source,
    }


def read_rows(content: bytes) -> list[tuple[int, int, int]]:
    rows = [
        (int(row["m"]), int(row["original_n"]), int(row["k"]))
        for row in csv.DictReader(io.StringIO(content.decode("utf-8")))
    ]
    if [m for m, _, _ in rows] != list(range(27)):
        raise AssertionError("expected exactly one published row for every m=0,...,26")
    if any(original_n != m + 1 for m, original_n, _ in rows):
        raise AssertionError("original_n must equal m+1")
    return rows


def verify(path: Path, random_cases: int) -> dict[str, Any]:
    path = path.resolve()
    input_bytes = path.read_bytes()
    rows = read_rows(input_bytes)
    factorizer = TrialFactorizer(max(m + 2 * k for m, _, k in rows))
    records: list[dict[str, Any]] = []

    for m, original_n, k in rows:
        certificate = witness_certificate(m, k, factorizer)
        carry_agrees = all(
            left == carry_count(k, m, p)
            and certificate.right_valuations[p] == carry_count(k, m + k, p)
            for p, left in certificate.left_valuations.items()
        )
        if not certificate.is_witness or not carry_agrees:
            raise AssertionError(f"published row failed at (m,k)=({m},{k})")

        shift: dict[str, Any] | None = None
        if k >= 2:
            obstructions = natural_shift_obstructions(m, k, factorizer)
            shift = {
                "target_m": m + 1,
                "target_k": k - 1,
                "is_witness": not obstructions,
                "obstructions": [
                    {
                        "prime": item.prime,
                        "source_slack": item.source_slack,
                        "v_p_m_plus_one": item.m_plus_one_valuation,
                        "v_p_m_plus_two_k": item.m_plus_two_k_valuation,
                        "shifted_slack": item.shifted_slack,
                    }
                    for item in obstructions
                ],
            }
            if m % 2 == 1 and any(item.prime > m for item in obstructions):
                raise AssertionError(f"odd-m shift bound failed at (m,k)=({m},{k})")

        records.append(
            {
                "m": m,
                "original_n": original_n,
                "k": k,
                "left_prime_count": len(certificate.left_valuations),
                "tight_prime_count": len(certificate.tight_primes),
                "tight_primes": list(certificate.tight_primes),
                "minimum_slack": (
                    min(certificate.slacks.values()) if certificate.slacks else None
                ),
                "carry_check_agrees": carry_agrees,
                "natural_shift": shift,
            }
        )

    nontrivial = records[1:]
    if any(record["minimum_slack"] != 0 for record in nontrivial):
        raise AssertionError("a nontrivial published witness was not p-adically tight")

    rng = random.Random(389)
    direct_factorizer = TrialFactorizer(4_100)
    for _ in range(random_cases):
        m = rng.randrange(0, 31)
        k = rng.randrange(1, 2_001)
        direct = is_witness_direct(m, k, max_n=4_100)
        if is_witness(m, k, direct_factorizer) != direct:
            raise AssertionError(f"valuation/direct mismatch at (m,k)=({m},{k})")
        if is_witness_by_carries(m, k, direct_factorizer) != direct:
            raise AssertionError(f"carry/direct mismatch at (m,k)=({m},{k})")

    bounded_factorizer = SPFFactorizer(5_010)
    bounded_minima = {
        "m_3": first_witness(3, 207, bounded_factorizer),
        "m_5": first_witness(5, 2_475, bounded_factorizer),
    }
    if bounded_minima != {"m_3": 207, "m_5": 2_475}:
        raise AssertionError(f"bounded minima mismatch: {bounded_minima}")

    tight_counts = [record["tight_prime_count"] for record in nontrivial]
    return {
        "schema_version": 1,
        "status": "PASS",
        "claim_boundary": (
            "The listed k values are verified witnesses. Minimality is independently "
            "exhausted here only for m=3 and m=5."
        ),
        "input": input_provenance(path, input_bytes),
        "published_rows": len(records),
        "verified_witnesses": len(records),
        "nontrivial_tight_witnesses": len(nontrivial),
        "tight_prime_count_range": [min(tight_counts), max(tight_counts)],
        "random_cross_checks": random_cases,
        "bounded_minima": bounded_minima,
        "records": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--random-cases", type=int, default=2_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.random_cases < 0:
        raise ValueError("random-cases must be nonnegative")
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = verify(args.input.resolve(), args.random_cases)
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
                "verified_witnesses": result["verified_witnesses"],
                "random_cross_checks": result["random_cross_checks"],
                "tight_prime_count_range": result["tight_prime_count_range"],
                "output": str(args.output),
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
