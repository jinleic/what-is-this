#!/usr/bin/env python3
"""Certify the exact half-window translation eviction of large bad-window primes.

The producer checks three separate exact statements.  First, the uniform
translation ``k -> k + ceil(m/2)`` evicts every prime ``p > m`` attached to the
bad window, with no modulus and no size growth.  Second, that offset is the
least uniform one: every smaller positive offset has an explicit surviving
blocker.  Third, on the exact small-prime tier the evicted candidate is
classified completely, both across a bounded rectangle and at the recorded
``m = 27`` moving-window survivors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from compensation_structure_analysis import (
    large_prime_window_obstructions,
    large_prime_window_records,
)
from erdos389 import (
    MAX_TRIAL_FACTORIZER_VALUE,
    SPFFactorizer,
    TrialFactorizer,
    primes_up_to,
    slack,
    witness_certificate,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_STRUCTURE = (
    ROOT / "data" / "compensation_structure_m1_20_k5000_m27_h200000000.json"
)
DEFAULT_ADAPTIVE = ROOT / "data" / "adaptive_eviction_m27_h200000000.json"
DEFAULT_OUTPUT = ROOT / "data" / "exact_eviction_m1_20_k5000_m27.json"

Factorizer = TrialFactorizer | SPFFactorizer


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def translation_offset(m: int) -> int:
    """Return the proved uniform eviction offset ``ceil(m/2)``."""
    _require_int("m", m, 1)
    return (m + 1) // 2


def bad_window_positions(m: int) -> range:
    """Return the exact bad-window index set ``{floor(m/2)+1, ..., m}``."""
    _require_int("m", m, 1)
    return range(m // 2 + 1, m + 1)


def window_large_primes(m: int, k: int, factorizer: Factorizer) -> dict[int, int]:
    """Map each prime ``p > m`` dividing the bad window to its unique position."""
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    attached: dict[int, int] = {}
    for position in bad_window_positions(m):
        for p in factorizer.factor(k + position):
            if p <= m:
                continue
            if p in attached:
                raise AssertionError("a prime larger than m hit two bad-window terms")
            attached[p] = position
    return dict(sorted(attached.items()))


def surviving_blocker_offsets(m: int, source_position: int) -> tuple[int, ...]:
    """Return every offset that keeps a blocker at ``source_position`` attached.

    A prime ``p > m`` attached at ``source_position`` divides a term of the
    window translated by ``delta`` exactly when ``delta`` is congruent to
    ``source_position - position`` for some window ``position``.  Only offsets
    below ``p`` are enumerated here, so the result is the exact set of small
    positive offsets that fail to evict that blocker.
    """
    _require_int("m", m, 1)
    positions = bad_window_positions(m)
    if source_position not in positions:
        raise ValueError("source_position must be a bad-window position")
    return tuple(
        sorted(
            source_position - position
            for position in positions
            if source_position - position > 0
        )
    )


def verify_translation_identity(max_m: int) -> dict[str, Any]:
    """Check the residue core of the eviction proof for every window pair."""
    _require_int("max_m", max_m, 1)
    counts = Counter()
    extrema: dict[str, int] = {}
    for m in range(1, max_m + 1):
        offset = translation_offset(m)
        positions = tuple(bad_window_positions(m))
        if len(positions) != offset:
            raise AssertionError("bad-window length must equal ceil(m/2)")
        for source_position in positions:
            for position in positions:
                residue = offset + position - source_position
                if not 1 <= residue <= m:
                    raise AssertionError(
                        "translated window residue left the interval [1, m]"
                    )
                counts["residue_identity_checks"] += 1
                extrema["maximum_residue"] = max(
                    extrema.get("maximum_residue", 0), residue
                )
                extrema["minimum_residue"] = min(
                    extrema.get("minimum_residue", m), residue
                )
    return {"counts": dict(sorted(counts.items())), "extrema": extrema}


def verify_translation_eviction(max_m: int, max_k: int) -> dict[str, Any]:
    """Exhaustively evict real window primes across a bounded rectangle."""
    _require_int("max_m", max_m, 1)
    _require_int("max_k", max_k, 1)
    factorizer = SPFFactorizer(max_m + max_k + translation_offset(max_m))
    counts = Counter()
    for m in range(1, max_m + 1):
        offset = translation_offset(m)
        positions = tuple(bad_window_positions(m))
        for k in range(1, max_k + 1):
            attached = window_large_primes(m, k, factorizer)
            counts["scanned_pairs"] += 1
            counts["attached_large_primes"] += len(attached)
            for p in attached:
                for position in positions:
                    if (k + offset + position) % p == 0:
                        raise AssertionError(
                            "half-window translation failed to evict a blocker"
                        )
                    counts["eviction_divisibility_checks"] += 1
    return {"counts": dict(sorted(counts.items()))}


def verify_translation_sharpness(max_m: int) -> dict[str, Any]:
    """Exhibit an explicit surviving blocker for every smaller uniform offset."""
    _require_int("max_m", max_m, 1)
    counts = Counter()
    rows: list[dict[str, int]] = []
    for m in range(1, max_m + 1):
        offset = translation_offset(m)
        first_position = m // 2 + 1
        primes = primes_up_to(2 * m + 2)
        larger = [p for p in primes if p > m]
        if not larger:
            raise AssertionError("Bertrand's range produced no prime above m")
        p = larger[0]
        for delta in range(1, offset):
            source_position = first_position + delta
            if source_position > m:
                raise AssertionError("sharpness position left the bad window")
            k = p - source_position
            if k < 1:
                raise AssertionError("sharpness construction produced k < 1")
            if (k + source_position) % p:
                raise AssertionError("sharpness prime missed the source window")
            survivor_position = source_position - delta
            if survivor_position != first_position:
                raise AssertionError("sharpness survivor left the first position")
            if (k + delta + survivor_position) % p:
                raise AssertionError("sharpness prime failed to survive the offset")
            if delta not in surviving_blocker_offsets(m, source_position):
                raise AssertionError("survivor offset missing from the exact set")
            counts["sharpness_constructions"] += 1
            rows.append(
                {
                    "m": m,
                    "delta": delta,
                    "prime": p,
                    "k": k,
                    "source_position": source_position,
                    "survivor_position": survivor_position,
                }
            )
        counts["least_uniform_offsets_confirmed"] += 1
    return {"counts": dict(sorted(counts.items())), "rows": rows}


def classify_eviction_candidate(
    m: int,
    k: int,
    factorizer: Factorizer,
    small_primes: tuple[int, ...],
) -> dict[str, Any]:
    """Return the exact classification of the half-window eviction candidate."""
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    offset = translation_offset(m)
    candidate = k + offset
    source_primes = window_large_primes(m, k, factorizer)
    for p in source_primes:
        for position in bad_window_positions(m):
            if (candidate + position) % p == 0:
                raise AssertionError("eviction candidate retained an old blocker")
    small_slacks = [
        {"prime": p, "slack": slack(m, candidate, p)} for p in small_primes
    ]
    small_tier_ok = all(row["slack"] >= 0 for row in small_slacks)
    obstructions = large_prime_window_obstructions(m, candidate, factorizer)
    candidate_x = m + 2 * candidate
    records = [
        {
            **record,
            "regime": (
                "single_level"
                if record["prime"] * record["prime"] > candidate_x
                else "compensation_deficit"
            ),
        }
        for record in obstructions
    ]
    if small_tier_ok and not records:
        category = "witness"
    elif not small_tier_ok and not records:
        category = "small_prime_tier_only"
    elif small_tier_ok:
        category = "new_large_prime_only"
    else:
        category = "both_tiers_failed"
    return {
        "m": m,
        "source_k": k,
        "offset": offset,
        "candidate": candidate,
        "candidate_x": candidate_x,
        "candidate_decimal_digits": len(str(candidate)),
        "source_large_prime_count": len(source_primes),
        "small_prime_slacks": small_slacks,
        "small_prime_tier_ok": small_tier_ok,
        "new_large_prime_obstructions": records,
        "category": category,
    }


def bounded_family(max_m: int, max_k: int, cross_check_limit: int) -> dict[str, Any]:
    """Classify every eviction candidate of a large-prime-obstructed rectangle."""
    _require_int("max_m", max_m, 1)
    _require_int("max_k", max_k, 1)
    _require_int("cross_check_limit", cross_check_limit, 0)
    factorizer = SPFFactorizer(max_m + 2 * (max_k + translation_offset(max_m)))
    counts = Counter()
    category_counts = Counter()
    obstruction_counts = Counter()
    small_primes_by_m = {
        m: tuple(primes_up_to(m)) for m in range(1, max_m + 1)
    }
    witness_rows: list[dict[str, Any]] = []
    for m in range(1, max_m + 1):
        small_primes = small_primes_by_m[m]
        for k in range(1, max_k + 1):
            counts["scanned_pairs"] += 1
            if not large_prime_window_obstructions(m, k, factorizer):
                continue
            counts["large_prime_obstructed_sources"] += 1
            row = classify_eviction_candidate(m, k, factorizer, small_primes)
            category_counts[row["category"]] += 1
            obstruction_counts[len(row["new_large_prime_obstructions"])] += 1
            counts["eviction_candidates_classified"] += 1
            if row["category"] == "witness":
                witness_rows.append(
                    {
                        "m": m,
                        "source_k": k,
                        "candidate": row["candidate"],
                    }
                )
            if k <= cross_check_limit:
                certificate = witness_certificate(m, row["candidate"], factorizer)
                if certificate.is_witness != (row["category"] == "witness"):
                    raise AssertionError(
                        "separation classification disagreed with the full certificate"
                    )
                counts["full_certificate_cross_checks"] += 1
    return {
        "counts": dict(sorted(counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "new_obstruction_count_histogram": [
            {"new_large_prime_obstructions": key, "candidates": value}
            for key, value in sorted(obstruction_counts.items())
        ],
        "witness_rows": witness_rows,
    }


def survivor_classification(
    structure: dict[str, Any],
    adaptive: dict[str, Any],
    target_m: int,
) -> dict[str, Any]:
    """Classify the recorded moving-window survivors under exact eviction."""
    _require_int("target_m", target_m, 1)
    section_names = (
        "m27_survivor_retained_prime_bounds",
        "m27_extended_survivor_retained_prime_bounds",
    )
    adaptive_digits = {
        record["source_k"]: record.get("candidate_decimal_digits")
        for record in adaptive["records"]
    }
    small_primes = tuple(primes_up_to(target_m))
    counts = Counter()
    category_counts = Counter()
    rows: list[dict[str, Any]] = []
    factorizer: TrialFactorizer | None = None
    for section_name in section_names:
        section = structure.get(section_name)
        if section is None:
            raise ValueError(f"structure artifact lacks {section_name}")
        for source in section["rows"]:
            k = source["k"]
            candidate = k + translation_offset(target_m)
            candidate_x = target_m + 2 * candidate
            if candidate_x > MAX_TRIAL_FACTORIZER_VALUE:
                raise AssertionError("survivor eviction candidate left the bound")
            if factorizer is None:
                factorizer = TrialFactorizer(MAX_TRIAL_FACTORIZER_VALUE)
            recorded = tuple(record["prime"] for record in source["obstructions"])
            if not recorded:
                raise AssertionError("survivor row has no recorded blocker")
            for record in source["obstructions"]:
                if (k + record["position"]) % record["prime"]:
                    raise AssertionError("recorded survivor attachment mismatch")
                counts["recorded_blocker_attachment_checks"] += 1
            row = classify_eviction_candidate(
                target_m, k, factorizer, small_primes
            )
            for p in recorded:
                for position in bad_window_positions(target_m):
                    if (candidate + position) % p == 0:
                        raise AssertionError("recorded blocker survived the offset")
                    counts["recorded_blocker_eviction_checks"] += 1
            window = large_prime_window_records(target_m, candidate, factorizer)
            counts["new_window_large_prime_factors"] += len(window)
            counts["survivor_systems"] += 1
            category_counts[row["category"]] += 1
            rows.append(
                {
                    **row,
                    "section": section_name,
                    "recorded_blocker_count": len(recorded),
                    "retained_prime_product_decimal_digits": len(
                        str(source["retained_prime_product"])
                    ),
                    "zero_carry_candidate_decimal_digits": source[
                        "zero_carry_eviction_candidate_decimal_digits"
                    ],
                    "adaptive_candidate_decimal_digits": adaptive_digits.get(k),
                    "new_window_large_prime_factor_count": len(window),
                }
            )
    digits = [row["candidate_decimal_digits"] for row in rows]
    adaptive_reported = [
        row["adaptive_candidate_decimal_digits"]
        for row in rows
        if row["adaptive_candidate_decimal_digits"] is not None
    ]
    zero_carry_reported = [
        row["zero_carry_candidate_decimal_digits"] for row in rows
    ]
    obstruction_totals = [
        len(row["new_large_prime_obstructions"]) for row in rows
    ]
    return {
        "counts": dict(sorted(counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "extrema": {
            "minimum_candidate_decimal_digits": min(digits),
            "maximum_candidate_decimal_digits": max(digits),
            "minimum_adaptive_candidate_decimal_digits": (
                min(adaptive_reported) if adaptive_reported else None
            ),
            "maximum_adaptive_candidate_decimal_digits": (
                max(adaptive_reported) if adaptive_reported else None
            ),
            "minimum_zero_carry_candidate_decimal_digits": min(zero_carry_reported),
            "maximum_zero_carry_candidate_decimal_digits": max(zero_carry_reported),
            "minimum_new_large_prime_obstructions": min(obstruction_totals),
            "maximum_new_large_prime_obstructions": max(obstruction_totals),
        },
        "rows": rows,
    }


def analyze(
    *,
    identity_max_m: int,
    eviction_max_m: int,
    eviction_max_k: int,
    family_max_m: int,
    family_max_k: int,
    cross_check_limit: int,
    target_m: int,
    structure_path: Path,
    adaptive_path: Path,
) -> dict[str, Any]:
    structure_path = structure_path.resolve()
    adaptive_path = adaptive_path.resolve()
    structure_bytes = structure_path.read_bytes()
    adaptive_bytes = adaptive_path.read_bytes()
    structure = json.loads(structure_bytes)
    adaptive = json.loads(adaptive_bytes)
    if structure.get("schema_version") != 3:
        raise ValueError("exact eviction requires structure schema 3")
    if structure.get("status") != "PASS_EXACT_COMPENSATION_STRUCTURE_ANALYSIS":
        raise ValueError("unexpected compensation-structure status")
    if adaptive.get("status") != "EXACT_BOUNDED_ADAPTIVE_EVICTION_SEARCH":
        raise ValueError("unexpected adaptive-eviction status")

    identity = verify_translation_identity(identity_max_m)
    eviction = verify_translation_eviction(eviction_max_m, eviction_max_k)
    sharpness = verify_translation_sharpness(identity_max_m)
    family = bounded_family(family_max_m, family_max_k, cross_check_limit)
    survivors = survivor_classification(structure, adaptive, target_m)

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_HALF_WINDOW_EVICTION",
        "claim_boundary": (
            "The translation identity and its sharpness constructions are exact "
            "over their declared m ranges and mirror the proof in THEOREMS.md. "
            "The rectangle eviction scan factors every bad-window term exactly. "
            "Each classified eviction candidate is decided by the proved "
            "small/large separation, with every large-prime factor of its "
            "translated window exactly factored inside the retained bound. No "
            "statement is made about candidates outside these finite domains."
        ),
        "parameters": {
            "identity_max_m": identity_max_m,
            "eviction_max_m": eviction_max_m,
            "eviction_max_k": eviction_max_k,
            "family_max_m": family_max_m,
            "family_max_k": family_max_k,
            "family_cross_check_max_k": cross_check_limit,
            "target_m": target_m,
            "translation": "K = k + ceil(m/2)",
        },
        "translation_identity": identity,
        "translation_eviction": eviction,
        "translation_sharpness": sharpness,
        "bounded_family": family,
        "survivor_eviction": survivors,
        "input": {
            "structure_logical_name": structure_path.name,
            "structure_sha256": hashlib.sha256(structure_bytes).hexdigest(),
            "structure_schema_version": structure["schema_version"],
            "adaptive_logical_name": adaptive_path.name,
            "adaptive_sha256": hashlib.sha256(adaptive_bytes).hexdigest(),
            "adaptive_schema_version": adaptive["schema_version"],
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "exact_eviction_certificate.py": sha256_file(Path(__file__)),
            "compensation_structure_analysis.py": sha256_file(
                ROOT / "compensation_structure_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    parser.add_argument("--adaptive", type=Path, default=DEFAULT_ADAPTIVE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--identity-max-m", type=int, default=200)
    parser.add_argument("--eviction-max-m", type=int, default=20)
    parser.add_argument("--eviction-max-k", type=int, default=5_000)
    parser.add_argument("--family-max-m", type=int, default=20)
    parser.add_argument("--family-max-k", type=int, default=5_000)
    parser.add_argument("--family-cross-check-max-k", type=int, default=500)
    parser.add_argument("--target-m", type=int, default=27)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        identity_max_m=args.identity_max_m,
        eviction_max_m=args.eviction_max_m,
        eviction_max_k=args.eviction_max_k,
        family_max_m=args.family_max_m,
        family_max_k=args.family_max_k,
        cross_check_limit=args.family_cross_check_max_k,
        target_m=args.target_m,
        structure_path=args.structure,
        adaptive_path=args.adaptive,
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
                "identity_checks": result["translation_identity"]["counts"],
                "eviction_checks": result["translation_eviction"]["counts"],
                "sharpness_checks": result["translation_sharpness"]["counts"],
                "family_categories": result["bounded_family"]["category_counts"],
                "survivor_categories": result["survivor_eviction"][
                    "category_counts"
                ],
                "survivor_extrema": result["survivor_eviction"]["extrema"],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
