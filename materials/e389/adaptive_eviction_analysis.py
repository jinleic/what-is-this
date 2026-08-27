#!/usr/bin/env python3
"""Search one blocker-evicting progression for uniform small-prime-safe terms.

For each certified obstructed window, the progression K_t = k + L + tP keeps
all old blocker primes out of the translated bad window.  Because P is
coprime to the uniform small-prime modulus, t runs through every small-prime
residue class in one period.  This script finds the first locally safe t only
within the declared finite search bound; it does not factor the resulting new
window when the candidate exceeds the retained exact factorizer limit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter
from math import gcd, prod
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import MAX_TRIAL_FACTORIZER_VALUE, primes_up_to, slack
from residue_control_analysis import (
    large_prime_bad_window_formula,
    local_safe_residues,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_STRUCTURE = (
    ROOT / "data" / "compensation_structure_m1_20_k5000_m27_h200000000.json"
)
DEFAULT_OUTPUT = ROOT / "data" / "adaptive_eviction_m27_h200000000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(
    *,
    target_m: int,
    search_limit: int,
    trial_prime_limit: int,
    structure_path: Path,
) -> dict[str, Any]:
    if target_m < 1 or search_limit < 0:
        raise ValueError("require target_m>=1 and search_limit>=0")
    if target_m != 27:
        raise ValueError("the canonical structure sections are specific to m=27")
    if trial_prime_limit <= target_m:
        raise ValueError("trial_prime_limit must exceed target_m")
    structure_path = structure_path.resolve()
    structure_bytes = structure_path.read_bytes()
    structure = json.loads(structure_bytes)
    if structure.get("schema_version") != 3:
        raise ValueError("adaptive eviction requires structure schema 3")
    if structure.get("status") != "PASS_EXACT_COMPENSATION_STRUCTURE_ANALYSIS":
        raise ValueError("unexpected compensation-structure status")

    small_rows: list[dict[str, Any]] = []
    small_system: list[tuple[int, int, frozenset[int]]] = []
    small_modulus = 1
    safe_class_count = 1
    for p in primes_up_to(target_m):
        modulus, residues = local_safe_residues(target_m, p)
        safe = frozenset(residues)
        if not safe:
            raise AssertionError("uniform local-safe residue set is empty")
        small_system.append((p, modulus, safe))
        small_modulus *= modulus
        safe_class_count *= len(safe)
        small_rows.append(
            {
                "prime": p,
                "modulus": modulus,
                "safe_residue_count": len(safe),
            }
        )
    density_divisor = gcd(safe_class_count, small_modulus)

    section_names = (
        "m27_survivor_retained_prime_bounds",
        "m27_extended_survivor_retained_prime_bounds",
    )
    records: list[dict[str, Any]] = []
    counts = Counter()
    first_position = target_m // 2 + 1
    bad_window_length = target_m - first_position + 1
    trial_primes = tuple(
        p for p in primes_up_to(trial_prime_limit) if p > target_m
    )

    for section_name in section_names:
        section = structure.get(section_name)
        if section is None:
            raise ValueError(f"structure artifact lacks {section_name}")
        for source in section["rows"]:
            k = source["k"]
            blockers = tuple(source["obstructions"])
            if not blockers:
                raise AssertionError("adaptive source has no large-prime blocker")
            blocker_product = prod(record["prime"] for record in blockers)
            if blocker_product != source["retained_prime_product"]:
                raise AssertionError("recorded blocker product mismatch")
            if gcd(blocker_product, small_modulus) != 1:
                raise AssertionError("small and blocker moduli are not coprime")

            base_candidate = k + bad_window_length
            first_t: int | None = None
            candidate: int | None = None
            for t in range(search_limit + 1):
                value = base_candidate + t * blocker_product
                if all(
                    value % modulus in safe
                    for _, modulus, safe in small_system
                ):
                    first_t = t
                    candidate = value
                    break
            counts["source_systems"] += 1
            counts["old_blocker_prime_checks"] += len(blockers)
            if candidate is None or first_t is None:
                counts["systems_without_safe_t_in_bound"] += 1
                records.append(
                    {
                        "section": section_name,
                        "source_k": k,
                        "old_blocker_count": len(blockers),
                        "blocker_product": blocker_product,
                        "first_safe_t": None,
                    }
                )
                continue

            small_slacks = []
            for p, modulus, safe in small_system:
                if candidate % modulus not in safe:
                    raise AssertionError("candidate escaped its local-safe class")
                value = slack(target_m, candidate, p)
                if value < 0:
                    raise AssertionError(
                        "uniform local-safe candidate failed small tier"
                    )
                small_slacks.append({"prime": p, "slack": value})
                counts["small_prime_slack_checks"] += 1
            for blocker in blockers:
                p = blocker["prime"]
                position = blocker["position"]
                if (k + position) % p:
                    raise AssertionError("source blocker attachment mismatch")
                for new_position in range(first_position, target_m + 1):
                    if (candidate + new_position) % p == 0:
                        raise AssertionError("old blocker entered translated window")
                    counts["old_blocker_window_divisibility_checks"] += 1

            trial_obstruction: dict[str, int] | None = None
            for position in range(first_position, target_m + 1):
                term = candidate + position
                for p in trial_primes:
                    counts["trial_prime_divisibility_checks"] += 1
                    if term % p:
                        continue
                    counts["trial_prime_divisor_hits"] += 1
                    formula = large_prime_bad_window_formula(
                        target_m, candidate, p
                    )
                    if formula["position"] != position:
                        raise AssertionError("trial divisor position mismatch")
                    if formula["predicted_slack"] < 0:
                        trial_obstruction = {"prime": p, **formula}
                        break
                if trial_obstruction is not None:
                    break
            if trial_obstruction is None:
                counts["candidates_unclassified_after_trial_scan"] += 1
            else:
                counts["candidates_rejected_by_trial_prime"] += 1

            candidate_x = target_m + 2 * candidate
            exceeds_bound = candidate_x > MAX_TRIAL_FACTORIZER_VALUE
            counts["systems_with_safe_t_in_bound"] += 1
            counts["tested_t_values_through_first_safe"] += first_t + 1
            if exceeds_bound:
                counts["safe_candidates_beyond_factorizer_bound"] += 1
            records.append(
                {
                    "section": section_name,
                    "source_k": k,
                    "old_blocker_count": len(blockers),
                    "blocker_product": blocker_product,
                    "blocker_product_decimal_digits": len(str(blocker_product)),
                    "progression_base": base_candidate,
                    "first_safe_t": first_t,
                    "tested_t_count": first_t + 1,
                    "candidate": candidate,
                    "candidate_x": candidate_x,
                    "candidate_decimal_digits": len(str(candidate)),
                    "zero_carry_candidate_decimal_digits": source[
                        "zero_carry_eviction_candidate_decimal_digits"
                    ],
                    "candidate_exceeds_factorizer_bound": exceeds_bound,
                    "small_prime_slacks": small_slacks,
                    "trial_prime_obstruction": trial_obstruction,
                }
            )

    found_t = [
        record["first_safe_t"]
        for record in records
        if record["first_safe_t"] is not None
    ]
    found_digits = [
        record["candidate_decimal_digits"]
        for record in records
        if record["first_safe_t"] is not None
    ]
    return {
        "schema_version": 1,
        "status": "EXACT_BOUNDED_ADAPTIVE_EVICTION_SEARCH",
        "claim_boundary": (
            "For each recorded old-blocker system, every t from zero through the "
            "reported first safe t was checked exactly. The resulting candidate "
            "satisfies a proved uniform small-prime class and excludes every old "
            "blocker. The bounded trial-prime scan supplies exact rejections when "
            "reported. Candidates beyond the retained factorizer bound are not "
            "fully factored, so all other newly entering primes remain uncontrolled."
        ),
        "parameters": {
            "target_m": target_m,
            "search_t_min": 0,
            "search_t_max": search_limit,
            "bad_window_first_position": first_position,
            "bad_window_length": bad_window_length,
            "trial_prime_limit": trial_prime_limit,
            "progression": (
                "K_t = source_k + bad_window_length + t * blocker_product"
            ),
        },
        "uniform_small_prime_system": {
            "modulus": small_modulus,
            "safe_class_count": safe_class_count,
            "safe_t_count_per_full_period": safe_class_count,
            "density": {
                "numerator": safe_class_count // density_divisor,
                "denominator": small_modulus // density_divisor,
            },
            "prime_rows": small_rows,
        },
        "counts": dict(sorted(counts.items())),
        "extrema": {
            "minimum_first_safe_t": min(found_t) if found_t else None,
            "maximum_first_safe_t": max(found_t) if found_t else None,
            "minimum_candidate_decimal_digits": (
                min(found_digits) if found_digits else None
            ),
            "maximum_candidate_decimal_digits": (
                max(found_digits) if found_digits else None
            ),
        },
        "records": records,
        "input": {
            "structure_logical_name": structure_path.name,
            "structure_sha256": hashlib.sha256(structure_bytes).hexdigest(),
            "structure_schema_version": structure["schema_version"],
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "adaptive_eviction_analysis.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--search-limit", type=int, default=1_000_000)
    parser.add_argument("--trial-prime-limit", type=int, default=1_000_000)
    parser.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        target_m=args.target_m,
        search_limit=args.search_limit,
        trial_prime_limit=args.trial_prime_limit,
        structure_path=args.structure,
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
                "extrema": result["extrema"],
                "wall_seconds": result["resources"]["wall_seconds"],
                "process_cpu_seconds": result["resources"]["process_cpu_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
