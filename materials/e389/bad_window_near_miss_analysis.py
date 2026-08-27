#!/usr/bin/env python3
"""Analyze the closest moving-bad-window misses and explicit CRT repairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
from collections import Counter, deque
from math import isqrt
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    TrialFactorizer,
    primes_up_to,
    slack,
    witness_certificate,
)
from moving_bad_window_search import segmented_large_cofactors
from residue_control_analysis import (
    crt_coprime,
    large_prime_bad_window_formula,
    local_safe_residues,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_KNOWN = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_OUTPUT = ROOT / "data" / "bad_window_near_miss_m27_h1000000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_known(content: bytes) -> dict[int, int]:
    rows = csv.DictReader(io.StringIO(content.decode("utf-8")))
    return {int(row["m"]): int(row["k"]) for row in rows}


def _obstruction_records(m: int, k: int, factorizer: TrialFactorizer) -> list[dict[str, int]]:
    certificate = witness_certificate(m, k, factorizer)
    return [
        {
            "prime": p,
            "left_valuation": certificate.left_valuations[p],
            "right_valuation": certificate.right_valuations[p],
            "slack": value,
        }
        for p, value in certificate.obstructions.items()
    ]


def _nearest_forward_residue(k: int, residues: tuple[int, ...], modulus: int) -> int:
    return min(residues, key=lambda residue: ((residue - k) % modulus, residue))


def _minimum_fatal_windows(
    *,
    m: int,
    start_k: int,
    offset_limit: int,
    chunk_size: int,
    representative_limit: int,
) -> tuple[int, int, list[dict[str, Any]]]:
    first_position = m // 2 + 1
    length = m - first_position + 1
    first_term = start_k + first_position
    integer_count = offset_limit + length - 1
    last_term = first_term + integer_count - 1
    first_x = m + 2 * start_k
    if last_term >= first_x:
        raise ValueError("require last bad-window term < first x for exact residual sieve")
    trial_primes = primes_up_to(isqrt(last_term))
    window: deque[tuple[int, int]] = deque(maxlen=length)
    minimum = length + 1
    minimum_count = 0
    representatives: list[dict[str, Any]] = []

    processed = 0
    while processed < integer_count:
        count = min(chunk_size, integer_count - processed)
        segment_start = first_term + processed
        residuals = segmented_large_cofactors(segment_start, count, trial_primes)
        for local_index, residual in enumerate(residuals):
            term = segment_start + local_index
            window.append((term, residual))
            global_index = processed + local_index
            if global_index + 1 < length:
                continue
            offset = global_index - length + 1
            if offset >= offset_limit:
                continue
            k = start_k + offset
            x = m + 2 * k
            fatal = [
                {
                    "position": position,
                    "term": term_value,
                    "prime": prime,
                }
                for position, (term_value, prime) in enumerate(
                    window, start=first_position
                )
                if prime > m and prime * prime > x
            ]
            fatal_count = len(fatal)
            if fatal_count < minimum:
                minimum = fatal_count
                minimum_count = 0
                representatives = []
            if fatal_count == minimum:
                minimum_count += 1
                if len(representatives) < representative_limit:
                    representatives.append(
                        {"offset": offset, "k": k, "fatal_obstructions": fatal}
                    )
        processed += count

    return minimum, minimum_count, representatives


def analyze(
    *,
    target_m: int,
    start_k: int | None,
    offset_limit: int,
    chunk_size: int,
    representative_limit: int,
    known_path: Path,
) -> dict[str, Any]:
    if min(target_m, offset_limit, chunk_size, representative_limit) < 1:
        raise ValueError("all numeric bounds must be positive")

    known_path = known_path.resolve()
    known_bytes = known_path.read_bytes()
    known = _read_known(known_bytes)
    expected_start = known.get(target_m - 1)
    expected_start = expected_start - 1 if expected_start is not None else None
    if start_k is None:
        if expected_start is None:
            raise ValueError(
                "start_k is required when no published source row exists"
            )
        start_k = expected_start
    if start_k < 1:
        raise ValueError("start_k must be positive")
    max_x = target_m + 2 * (start_k + offset_limit - 1)
    if max_x > MAX_TRIAL_FACTORIZER_VALUE:
        raise ValueError("analysis exceeds the retained exact factorizer bound")

    minimum, minimum_count, representatives = _minimum_fatal_windows(
        m=target_m,
        start_k=start_k,
        offset_limit=offset_limit,
        chunk_size=chunk_size,
        representative_limit=representative_limit,
    )

    base_factorizer = TrialFactorizer(max_x)
    controlled_systems: list[dict[str, Any]] = []
    counts = Counter()
    mirror_candidates: list[int] = []
    for representative in representatives:
        k = representative["k"]
        obstructions = _obstruction_records(target_m, k, base_factorizer)
        representative["exact_obstructions"] = obstructions
        counts["representatives_exactly_checked"] += 1
        counts["representative_obstruction_primes"] += len(obstructions)

        congruences: list[tuple[int, int]] = []
        small_controls: list[dict[str, int]] = []
        for p in primes_up_to(target_m):
            modulus, residues = local_safe_residues(target_m, p)
            residue = _nearest_forward_residue(k, residues, modulus)
            congruences.append((residue, modulus))
            small_controls.append(
                {"prime": p, "residue": residue, "modulus": modulus}
            )

        large_controls: list[dict[str, int]] = []
        for obstruction in obstructions:
            p = obstruction["prime"]
            if p <= target_m:
                continue
            formula = large_prime_bad_window_formula(target_m, k, p)
            exponent = formula["exponent"]
            power = p**exponent
            modulus = p ** (2 * exponent if p % 2 else 2 * exponent + 1)
            residue = (-formula["position"] - power) % modulus
            congruences.append((residue, modulus))
            large_controls.append(
                {
                    "prime": p,
                    "position": formula["position"],
                    "exponent": exponent,
                    "residue": residue,
                    "modulus": modulus,
                }
            )
            forward_delta = (residue - k) % modulus
            mirror_candidates.append(k + forward_delta)

        solution, modulus = crt_coprime(congruences)
        if solution == 0:
            solution = modulus
        for control in small_controls:
            if slack(target_m, solution, control["prime"]) < 0:
                raise AssertionError("combined CRT lost a small-prime local control")
            counts["combined_small_prime_controls_checked"] += 1
        for control in large_controls:
            if slack(target_m, solution, control["prime"]) < 0:
                raise AssertionError("combined CRT lost a large-prime mirror control")
            counts["combined_large_prime_controls_checked"] += 1
        controlled_systems.append(
            {
                "source_k": k,
                "combined_solution": solution,
                "combined_modulus": modulus,
                "within_factorizer_bound": (
                    target_m + 2 * solution <= MAX_TRIAL_FACTORIZER_VALUE
                ),
                "small_prime_controls": small_controls,
                "large_prime_controls": large_controls,
            }
        )

    evaluable_mirrors = sorted(
        {
            k
            for k in mirror_candidates
            if target_m + 2 * k <= MAX_TRIAL_FACTORIZER_VALUE
        }
    )
    mirror_evaluations: list[dict[str, Any]] = []
    if evaluable_mirrors:
        mirror_factorizer = TrialFactorizer(target_m + 2 * max(evaluable_mirrors))
        for k in evaluable_mirrors:
            obstructions = _obstruction_records(target_m, k, mirror_factorizer)
            mirror_evaluations.append({"k": k, "obstructions": obstructions})
            counts["mirror_candidates_exactly_checked"] += 1
            if obstructions:
                counts["mirror_candidates_with_new_prime_interference"] += 1
            else:
                counts["mirror_candidate_witnesses"] += 1
    counts["mirror_candidates_beyond_factorizer_bound"] = (
        len(set(mirror_candidates)) - len(evaluable_mirrors)
    )

    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_BAD_WINDOW_NEAR_MISS_ANALYSIS",
        "claim_boundary": (
            "The minimum is exact only over the declared offset interval. Full witness "
            "checks cover every reported minimum representative. CRT systems certify "
            "only their listed primes; unevaluated large solutions are not witnesses."
        ),
        "parameters": {
            "target_m": target_m,
            "start_k": start_k,
            "offset_count": offset_limit,
            "offset_min": 0,
            "offset_max": offset_limit - 1,
            "chunk_size": chunk_size,
            "representative_limit": representative_limit,
        },
        "starting_point_is_published_natural_shift": start_k == expected_start,
        "minimum_single_level_obstruction_count": minimum,
        "minimum_window_count": minimum_count,
        "minimum_window_representatives": representatives,
        "combined_crt_systems": controlled_systems,
        "mirror_candidate_evaluations": mirror_evaluations,
        "counts": dict(sorted(counts.items())),
        "input": {
            "known_logical_name": known_path.name,
            "known_sha256": hashlib.sha256(known_bytes).hexdigest(),
        },
        "implementation_sha256": {
            "bad_window_near_miss_analysis.py": sha256_file(Path(__file__)),
            "moving_bad_window_search.py": sha256_file(
                ROOT / "moving_bad_window_search.py"
            ),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--start-k", type=int)
    parser.add_argument("--offset-limit", type=int, default=1_000_000)
    parser.add_argument("--chunk-size", type=int, default=250_000)
    parser.add_argument("--representatives", type=int, default=100)
    parser.add_argument("--known", type=Path, default=DEFAULT_KNOWN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        target_m=args.target_m,
        start_k=args.start_k,
        offset_limit=args.offset_limit,
        chunk_size=args.chunk_size,
        representative_limit=args.representatives,
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
                "minimum_single_level_obstruction_count": result[
                    "minimum_single_level_obstruction_count"
                ],
                "minimum_window_count": result["minimum_window_count"],
                "counts": result["counts"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
