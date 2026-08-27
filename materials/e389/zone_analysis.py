#!/usr/bin/env python3
"""Finite exact certificate for prime-power zone formulas used in Erdős #389."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import (
    SPFFactorizer,
    left_valuations,
    local_floor_contribution,
    slack,
    zone_contribution,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "zone_theorem_m1_8_k2000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bad_window_offsets(m: int) -> range:
    return range(m // 2 + 1, m + 1)


def verify_zone_rectangle(max_m: int, max_k: int) -> dict[str, Any]:
    if max_m < 1 or max_k < 1:
        raise ValueError("max_m and max_k must be positive")
    factorizer = SPFFactorizer(max_m + 2 * max_k)
    counts = Counter()
    theorem_s_examples: list[dict[str, int]] = []

    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            counts["pairs"] += 1
            x = m + 2 * k
            n = m + k
            valuations = left_valuations(m, k, factorizer)
            pair_is_witness = True
            for p in valuations:
                counts["prime_checks"] += 1
                total = 0
                q = p
                while q <= x:
                    direct = local_floor_contribution(m, x, q)
                    zoned = zone_contribution(m, k, q)
                    if direct != zoned:
                        raise AssertionError(
                            f"zone mismatch at (m,k,p,q)=({m},{k},{p},{q})"
                        )
                    if direct not in (-1, 0, 1):
                        raise AssertionError("zone contribution escaped {-1,0,1}")
                    counts[{1: "good_levels", 0: "neutral_levels", -1: "bad_levels"}[direct]] += 1
                    counts["level_checks"] += 1
                    if n < q <= x:
                        if direct != 1:
                            raise AssertionError("free-good interval level was not +1")
                        counts["free_good_checks"] += 1
                    if k % q == 0:
                        if direct != 0:
                            raise AssertionError("q|k neutral congruence failed")
                        counts["neutral_congruence_checks"] += 1
                    if (k + m + 1) % q == 0:
                        expected = 0 if (m + 1) % q == 0 else 1
                        if direct != expected:
                            raise AssertionError("forcing congruence level failed")
                        counts["forcing_congruence_checks"] += 1
                    total += direct
                    q *= p

                exact_slack = slack(m, k, p)
                if total != exact_slack:
                    raise AssertionError(
                        f"zone sum mismatch at (m,k,p)=({m},{k},{p})"
                    )
                if exact_slack < 0:
                    pair_is_witness = False
                    if p > m and not any((k + s) % p == 0 for s in bad_window_offsets(m)):
                        raise AssertionError("large-prime obstruction missed bad window")
                    counts["obstruction_support_checks"] += 1

                if (
                    p > m
                    and p * p > x
                    and any((k + s) % p == 0 for s in bad_window_offsets(m))
                ):
                    if exact_slack != -1 or pair_is_witness:
                        raise AssertionError("single-level obstruction theorem failed")
                    counts["single_level_obstruction_checks"] += 1
                    if len(theorem_s_examples) < 50:
                        theorem_s_examples.append(
                            {"m": m, "k": k, "prime": p, "x": x}
                        )

            if pair_is_witness:
                counts["witnesses"] += 1

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_FINITE_CERTIFICATE",
        "claim_boundary": (
            "The zone, support, congruence, free-good, and single-level obstruction "
            "identities are exhaustively checked only on the declared rectangle. "
            "Their elementary proofs are recorded separately."
        ),
        "parameters": {"min_m": 1, "max_m": max_m, "min_k": 1, "max_k": max_k},
        "counts": dict(sorted(counts.items())),
        "single_level_obstruction_examples": theorem_s_examples,
        "implementation_sha256": {
            "zone_analysis.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=8)
    parser.add_argument("--max-k", type=int, default=2_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = verify_zone_rectangle(args.max_m, args.max_k)
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
                "pairs": result["counts"]["pairs"],
                "level_checks": result["counts"]["level_checks"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
