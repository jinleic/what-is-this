#!/usr/bin/env python3
"""Search exact witnesses on the small-prime-zero corridor k=t*M_m."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
import time
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    TrialFactorizer,
    carry_count,
    primes_up_to,
    witness_certificate,
    natural_shift_obstructions,
)
from shift_repair_analysis import zero_carry_modulus

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "zero_carry_corridor_m5_18_t2000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def search(
    m_min: int, m_max: int, t_limit: int, prime_bound: int | None = None
) -> dict[str, Any]:
    if (
        not 1 <= m_min <= m_max
        or t_limit < 1
        or (prime_bound is not None and prime_bound < 2)
    ):
        raise ValueError(
            "require 1 <= m_min <= m_max, t_limit >= 1, and prime_bound >= 2"
        )

    rows_spec = []
    maximum_value = 1
    for m in range(m_min, m_max + 1):
        controlled_bound = m if prime_bound is None else min(m, prime_bound)
        primes = primes_up_to(controlled_bound)
        modulus = zero_carry_modulus(m, primes)
        effective_limit = min(
            t_limit, (MAX_TRIAL_FACTORIZER_VALUE - m) // (2 * modulus)
        )
        rows_spec.append((m, primes, modulus, effective_limit))
        if effective_limit:
            maximum_value = max(
                maximum_value, 2 * modulus * effective_limit + m
            )

    factorizer = TrialFactorizer(maximum_value)
    rows: list[dict[str, Any]] = []
    candidate_count = 0
    small_prime_checks = 0
    witness_count = 0

    for m, primes, modulus, effective_limit in rows_spec:
        witnesses = []
        row_witness_count = 0
        failure_categories: Counter[str] = Counter()
        obstruction_prime_counts: Counter[int] = Counter()
        failure_representatives: list[dict[str, Any]] = []
        for t in range(1, effective_limit + 1):
            candidate_count += 1
            k = modulus * t
            for p in primes:
                if carry_count(k, m, p) != 0:
                    raise AssertionError(
                        f"zero-carry corridor failed at (m,t,p)=({m},{t},{p})"
                    )
                small_prime_checks += 1
            certificate = witness_certificate(m, k, factorizer)
            if not certificate.is_witness:
                categories = set()
                x = m + 2 * k
                for p, deficit in certificate.obstructions.items():
                    obstruction_prime_counts[p] += 1
                    if p <= m:
                        if p in primes:
                            raise AssertionError(
                                "controlled small prime obstructed zero-carry corridor"
                            )
                        categories.add("uncontrolled_p_le_m")
                        continue
                    if not any(
                        (k + offset) % p == 0
                        for offset in range(m // 2 + 1, m + 1)
                    ):
                        raise AssertionError(
                            "large obstruction prime missed the bad window"
                        )
                    categories.add(
                        "single_level_large_prime"
                        if p * p > x
                        else "multi_level_large_prime"
                    )
                    if deficit >= 0:
                        raise AssertionError("recorded obstruction had nonnegative slack")
                for category in categories:
                    failure_categories[category] += 1
                if len(failure_representatives) < 20:
                    failure_representatives.append(
                        {
                            "t": t,
                            "k": k,
                            "categories": sorted(categories),
                            "obstructions": [
                                {"prime": p, "slack": value}
                                for p, value in sorted(
                                    certificate.obstructions.items()
                                )
                            ],
                        }
                    )
                continue
            witness_count += 1
            row_witness_count += 1
            if len(witnesses) < 20:
                shift_obstructions = natural_shift_obstructions(
                    m, k, factorizer
                )
                witnesses.append(
                    {
                        "t": t,
                        "k": k,
                        "left_prime_count": len(certificate.left_valuations),
                        "tight_prime_count": len(certificate.tight_primes),
                        "tight_primes": list(certificate.tight_primes),
                        "natural_shift_is_witness": not shift_obstructions,
                        "natural_shift_obstructions": [
                            {
                                "prime": item.prime,
                                "shifted_slack": item.shifted_slack,
                            }
                            for item in shift_obstructions
                        ],
                    }
                )
        rows.append(
            {
                "m": m,
                "modulus": modulus,
                "controlled_primes": primes,
                "requested_t_limit": t_limit,
                "effective_t_limit": effective_limit,
                "candidate_count": effective_limit,
                "witness_count": row_witness_count,
                "witness_representatives": witnesses,
                "truncated_representatives": row_witness_count > len(witnesses),
                "failure_count": effective_limit - row_witness_count,
                "failure_candidate_categories": dict(
                    sorted(failure_categories.items())
                ),
                "obstruction_prime_counts": [
                    {"prime": p, "count": count}
                    for p, count in sorted(
                        obstruction_prime_counts.items(),
                        key=lambda item: (-item[1], item[0]),
                    )
                ],
                "failure_representatives": failure_representatives,
            }
        )

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_CANDIDATE_FAMILY",
        "claim_boundary": (
            "Every declared k=t*M was checked, where M kills carries only for "
            "the recorded controlled_primes. Failure in this candidate family is "
            "not failure of Erdős #389."
        ),
        "parameters": {
            "m_min": m_min,
            "m_max": m_max,
            "t_limit": t_limit,
            "prime_bound": prime_bound,
        },
        "resource_bound": {
            "maximum_factored_value": maximum_value,
            "factorizer_cap": MAX_TRIAL_FACTORIZER_VALUE,
        },
        "finite_checks": {
            "candidate_count": candidate_count,
            "small_prime_zero_carry_checks": small_prime_checks,
            "witness_count": witness_count,
        },
        "rows": rows,
        "implementation_sha256": {
            "zero_carry_corridor_search.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "shift_repair_analysis.py": sha256_file(
                ROOT / "shift_repair_analysis.py"
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m-min", type=int, default=5)
    parser.add_argument("--m-max", type=int, default=18)
    parser.add_argument("--t-limit", type=int, default=2_000)
    parser.add_argument("--prime-bound", type=int)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = search(
        args.m_min, args.m_max, args.t_limit, prime_bound=args.prime_bound
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
                "candidate_count": result["finite_checks"]["candidate_count"],
                "witness_count": result["finite_checks"]["witness_count"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
