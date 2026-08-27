#!/usr/bin/env python3
"""Exact bounded analysis of small-prime CRT control and shift spike barriers."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from bisect import bisect_right
from collections import Counter
from math import prod
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import (
    RatioSlackScanner,
    SPFFactorizer,
    carry_count,
    integer_valuation,
    primes_up_to,
    slack,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_ATLAS = ROOT / "data" / "atlas_m1_50_k50000.json"
DEFAULT_OUTPUT = ROOT / "data" / "shift_repair_analysis_m1_50_k50000.json"


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def zero_carry_modulus(m: int, primes: list[int]) -> int:
    modulus = 1
    for p in primes:
        power = p
        while power <= m:
            power *= p
        modulus *= power
    return modulus


def odd_shift_rough_modulus(m: int, primes: list[int]) -> tuple[int, int]:
    odd_primes = [p for p in primes if p != 2]
    return prod(odd_primes, start=1), prod((p - 1 for p in odd_primes), start=1)


def is_odd_shift_crt_safe(m: int, k: int, primes: list[int]) -> bool:
    x = m + 2 * k
    return all(x % p for p in primes)


def scan_witness_bits(min_m: int, max_m: int, max_k: int) -> list[bytearray | None]:
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    witness_bits: list[bytearray | None] = [None] * (max_m + 1)
    for m in range(min_m, max_m + 1):
        scanner = RatioSlackScanner(m, factorizer)
        bits = bytearray(max_k + 1)
        while True:
            if scanner.is_witness:
                bits[scanner.k] = 1
            if scanner.k == max_k:
                break
            scanner.advance()
        witness_bits[m] = bits
    return witness_bits


def analyze(
    *, min_m: int, max_m: int, max_k: int, atlas_path: Path
) -> dict[str, Any]:
    atlas_bytes = atlas_path.read_bytes()
    atlas = json.loads(atlas_bytes)
    if atlas["schema_version"] != 2:
        raise AssertionError("shift analysis requires repaired atlas schema 2")
    expected_parameters = {"min_m": min_m, "max_m": max_m, "max_k": max_k}
    for name, value in expected_parameters.items():
        if atlas["parameters"][name] != value:
            raise AssertionError(f"atlas parameter {name} does not match analysis")

    primes = primes_up_to(max_m)
    primes_by_m = {m: [p for p in primes if p <= m] for m in range(min_m, max_m + 1)}
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    witness_bits = scan_witness_bits(min_m, max_m, max_k)

    for row in atlas["m_rows"]:
        bits = witness_bits[row["m"]]
        assert bits is not None
        if sum(bits) != row["witness_count"]:
            raise AssertionError(f"witness count mismatch at m={row['m']}")

    zero_carry_rows: list[dict[str, Any]] = []
    zero_carry_checks = 0
    for m in range(min_m, max_m + 1):
        small_primes = primes_by_m[m]
        modulus = zero_carry_modulus(m, small_primes)
        bits = witness_bits[m]
        assert bits is not None
        candidate_count = max_k // modulus
        witness_multiples: list[int] = []
        for k in range(modulus, max_k + 1, modulus):
            for p in small_primes:
                if carry_count(k, m, p) != 0:
                    raise AssertionError(
                        f"zero-carry modulus failed at (m,k,p)=({m},{k},{p})"
                    )
                zero_carry_checks += 1
            if bits[k]:
                witness_multiples.append(k)
        zero_carry_rows.append(
            {
                "m": m,
                "modulus": modulus,
                "candidate_multiples_within_bound": candidate_count,
                "witness_multiple_count": len(witness_multiples),
                "first_witness_multiple": (
                    witness_multiples[0] if witness_multiples else None
                ),
                "witness_multiple_representatives": witness_multiples[:20],
            }
        )

    odd_rows: list[dict[str, Any]] = []
    odd_failure_records: list[dict[str, Any]] = []
    odd_safe_checks = 0
    for m in range(min_m, max_m):
        if m % 2 == 0:
            continue
        source = witness_bits[m]
        target = witness_bits[m + 1]
        assert source is not None and target is not None
        small_primes = primes_by_m[m]
        modulus, allowed_residues = odd_shift_rough_modulus(m, small_primes)
        safe_witnesses = [
            k
            for k in range(2, max_k + 1)
            if source[k] and is_odd_shift_crt_safe(m, k, small_primes)
        ]
        counts = Counter()
        repair_delta_counts: Counter[int] = Counter()
        repair_period_counts: Counter[int] = Counter()

        for k in range(2, max_k + 1):
            if not source[k]:
                continue
            success = bool(target[k - 1])
            safe = is_odd_shift_crt_safe(m, k, small_primes)
            counts["eligible"] += 1
            counts["successful" if success else "failed"] += 1
            counts["safe" if safe else "unsafe"] += 1
            if safe:
                odd_safe_checks += 1
                if not success:
                    raise AssertionError(f"CRT-safe odd witness failed to shift: {(m, k)}")
                counts["safe_successful"] += 1
            elif success:
                counts["unsafe_successful"] += 1
            else:
                counts["unsafe_failed"] += 1
                first_safe_delta = None
                for delta in range(1, modulus + 1):
                    if is_odd_shift_crt_safe(m, k + delta, small_primes):
                        first_safe_delta = delta
                        break
                index = bisect_right(safe_witnesses, k)
                safe_witness_delta = None
                safe_witness_delta_within_bound = None
                periods_to_safe_witness = None
                if index < len(safe_witnesses):
                    candidate = safe_witnesses[index]
                    safe_witness_delta_within_bound = candidate - k
                    periods_to_safe_witness = (
                        safe_witness_delta_within_bound + modulus - 1
                    ) // modulus
                    repair_period_counts[periods_to_safe_witness] += 1
                    if safe_witness_delta_within_bound <= modulus:
                        safe_witness_delta = safe_witness_delta_within_bound
                        repair_delta_counts[safe_witness_delta] += 1
                odd_failure_records.append(
                    {
                        "m": m,
                        "k": k,
                        "rough_modulus": modulus,
                        "first_safe_residue_delta": first_safe_delta,
                        "first_safe_residue_is_witness": bool(
                            first_safe_delta is not None
                            and k + first_safe_delta <= max_k
                            and source[k + first_safe_delta]
                        ),
                        "first_safe_witness_delta_within_one_period": safe_witness_delta,
                        "first_safe_witness_delta_within_bound": (
                            safe_witness_delta_within_bound
                        ),
                        "rough_periods_to_first_safe_witness": periods_to_safe_witness,
                    }
                )

        odd_rows.append(
            {
                "m": m,
                "odd_prime_modulus": modulus,
                "allowed_residue_count": allowed_residues,
                "eligible_source_witnesses": counts["eligible"],
                "successful_shifts": counts["successful"],
                "failed_shifts": counts["failed"],
                "crt_safe_source_witnesses": counts["safe"],
                "crt_safe_successful_shifts": counts["safe_successful"],
                "crt_unsafe_successful_shifts": counts["unsafe_successful"],
                "crt_unsafe_failed_shifts": counts["unsafe_failed"],
                "failed_sources_repaired_by_safe_witness_within_one_period": sum(
                    repair_delta_counts.values()
                ),
                "failed_sources_repaired_by_safe_witness_within_bound": sum(
                    repair_period_counts.values()
                ),
                "safe_witness_repair_period_counts": [
                    {"periods": periods, "count": count}
                    for periods, count in sorted(repair_period_counts.items())
                ],
                "safe_witness_repair_delta_counts": [
                    {"delta": delta, "count": count}
                    for delta, count in sorted(repair_delta_counts.items())
                ],
            }
        )

    even_rows: list[dict[str, Any]] = []
    fatal_examples: list[dict[str, Any]] = []
    fatal_prime_checks = 0
    for m in range(min_m, max_m):
        if m % 2:
            continue
        source = witness_bits[m]
        target = witness_bits[m + 1]
        assert source is not None and target is not None
        counts = Counter()
        for k in range(2, max_k + 1):
            if not source[k]:
                continue
            success = bool(target[k - 1])
            counts["eligible"] += 1
            counts["successful" if success else "failed"] += 1
            x = m + 2 * k
            fatal = []
            for p, exponent in sorted(factorizer.factor(x).items()):
                if p <= m + 1 or p ** (exponent + 1) <= x + m:
                    continue
                shifted_slack = (
                    slack(m, k, p) + integer_valuation(m + 1, p) - exponent
                )
                if shifted_slack >= 0:
                    raise AssertionError(
                        f"fatal spike lemma failed at (m,k,p)=({m},{k},{p})"
                    )
                fatal_prime_checks += 1
                fatal.append(
                    {
                        "prime": p,
                        "exponent": exponent,
                        "shifted_slack": shifted_slack,
                        "square_root_corollary": p * p > x + m,
                    }
                )
            if fatal:
                counts["with_fatal_spike"] += 1
                if success:
                    raise AssertionError(f"fatal spike survived natural shift: {(m, k)}")
                counts["failed_explained_by_fatal_spike"] += 1
                if len(fatal_examples) < 50:
                    fatal_examples.append({"m": m, "k": k, "x": x, "fatal": fatal})
            elif not success:
                counts["failed_without_fatal_spike"] += 1

        even_rows.append(
            {
                "m": m,
                "eligible_source_witnesses": counts["eligible"],
                "successful_shifts": counts["successful"],
                "failed_shifts": counts["failed"],
                "sources_with_fatal_spike": counts["with_fatal_spike"],
                "failed_shifts_explained_by_fatal_spike": counts[
                    "failed_explained_by_fatal_spike"
                ],
                "failed_shifts_without_fatal_spike": counts[
                    "failed_without_fatal_spike"
                ],
            }
        )

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_ANALYSIS",
        "claim_boundary": (
            "All counts are exact for the declared rectangle. The listed lemmas are "
            "elementary implications; bounded repair rates do not prove existence."
        ),
        "parameters": expected_parameters,
        "proved_lemmas": [
            {
                "name": "small_prime_zero_carry_modulus",
                "statement": (
                    "For each p<=m let a_p be minimal with p^a_p>m. If k is "
                    "divisible by product_p p^a_p, then v_p(C(m+k,m))=0 for "
                    "every p<=m."
                ),
            },
            {
                "name": "odd_shift_crt_roughness",
                "statement": (
                    "If m is odd, (m,k) is a witness, and m+2k is coprime to "
                    "the product of primes <=m, then (m+1,k-1) is a witness."
                ),
            },
            {
                "name": "even_shift_fatal_prime_power_spike",
                "statement": (
                    "Let m be even, (m,k) a witness, x=m+2k, and p>m+1 with "
                    "a=v_p(x)>=1. If p^(a+1)>x+m, then (m+1,k-1) fails at p."
                ),
            },
        ],
        "finite_checks": {
            "zero_carry_prime_checks": zero_carry_checks,
            "odd_crt_safe_shift_checks": odd_safe_checks,
            "fatal_prime_power_checks": fatal_prime_checks,
        },
        "zero_carry_corridor": zero_carry_rows,
        "odd_shift_crt_analysis": {
            "rows": odd_rows,
            "failure_records": odd_failure_records,
        },
        "even_shift_spike_analysis": {
            "rows": even_rows,
            "fatal_examples": fatal_examples,
        },
        "input": {
            "atlas_logical_name": atlas_path.name,
            "atlas_sha256": sha256_bytes(atlas_bytes),
            "atlas_schema_version": atlas["schema_version"],
        },
        "implementation_sha256": {
            "shift_repair_analysis.py": sha256_file(Path(__file__)),
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        min_m=args.min_m,
        max_m=args.max_m,
        max_k=args.max_k,
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
                "zero_carry_checks": result["finite_checks"][
                    "zero_carry_prime_checks"
                ],
                "odd_safe_checks": result["finite_checks"][
                    "odd_crt_safe_shift_checks"
                ],
                "fatal_spike_checks": result["finite_checks"][
                    "fatal_prime_power_checks"
                ],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
