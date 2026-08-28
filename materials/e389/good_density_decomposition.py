#!/usr/bin/env python3
"""Decompose exactly why integers fail to be compensation-good.

Section 13 reduces Erdos #389 to Conjecture R, and the first open case is
R(1): do the compensation-good integers have positive lower density? The
elementary union bound over per-prime failure events exceeds 1, so it cannot
answer that. The natural repair is to claim the level failures overlap the
short-cofactor failures, since a short cofactor already accounts for density
log 2, and to keep only the non-overlapping part.

This producer refutes that repair with exact counts. Every integer in a
declared block is classified into four disjoint classes:

  good                    every prime p > m has C_p(w) >= v_p(w)
  short cofactor only     the only failures come from p with p^2 > 2w
  level failure only      some p with p^2 <= 2w fails, and no p^2 > 2w does
  both                    both kinds of failure occur

Under independence the level-failure-only class would have density
(1 - log 2) times the total level-failure mass. The measured value is far
larger, because a large prime factor leaves a small cofactor with few primes
available to fail a level test. The two failure modes are adversely
correlated, so the overlap that would rescue the union bound does not exist.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from erdos389 import SPFFactorizer, primes_up_to

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "data" / "good_density_decomposition.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def compensation_count(cofactor: int, p: int) -> int:
    """Return ``C_p = #{r>=1 : 2*(cofactor mod p^r) > p^r}`` exactly."""
    _require_int("cofactor", cofactor, 1)
    _require_int("p", p, 2)
    count = 0
    power = p
    while power <= 2 * cofactor:
        if 2 * (cofactor % power) > power:
            count += 1
        power *= p
    return count


def classify_term(
    value: int, target_m: int, factorizer: SPFFactorizer
) -> tuple[bool, bool]:
    """Return ``(short_cofactor_failure, level_failure)`` for one integer."""
    _require_int("value", value, 1)
    _require_int("target_m", target_m, 0)
    short = False
    level = False
    for p, exponent in factorizer.factor(value).items():
        if p <= target_m:
            continue
        cofactor = value // p**exponent
        if compensation_count(cofactor, p) >= exponent:
            continue
        if p * p > 2 * value:
            short = True
        else:
            level = True
    return short, level


def decompose_block(start: int, count: int, target_m: int) -> dict[str, Any]:
    """Classify every integer of ``[start, start+count)`` exactly."""
    _require_int("start", start, 2)
    _require_int("count", count, 1)
    _require_int("target_m", target_m, 0)
    factorizer = SPFFactorizer(start + count + 1)
    classes = {"good": 0, "short_only": 0, "level_only": 0, "both": 0}
    for value in range(start, start + count):
        short, level = classify_term(value, target_m, factorizer)
        if short and level:
            classes["both"] += 1
        elif short:
            classes["short_only"] += 1
        elif level:
            classes["level_only"] += 1
        else:
            classes["good"] += 1
    total = sum(classes.values())
    if total != count:
        raise AssertionError("classes do not partition the block")
    densities = {name: value / count for name, value in classes.items()}
    level_mass = densities["level_only"] + densities["both"]
    short_mass = densities["short_only"] + densities["both"]
    independent_level_only = level_mass * (1.0 - short_mass)
    return {
        "start": start,
        "count": count,
        "target_m": target_m,
        "counts": classes,
        "densities": densities,
        "total_level_failure_mass": level_mass,
        "total_short_cofactor_mass": short_mass,
        "independence_predicted_level_only": independent_level_only,
        "observed_over_predicted_level_only": (
            densities["level_only"] / independent_level_only
            if independent_level_only > 0
            else None
        ),
        "log_two_deviation": short_mass - math.log(2),
    }


MAX_UNION_LIMIT = 10**13


def union_bound_terms(limit: int, target_m: int) -> dict[str, Any]:
    """Model the per-prime level-failure sum that makes the union bound fail.

    Sieves to ``sqrt(2*limit)``, so the limit is capped; use
    :func:`asymptotic_union_bound` for the behaviour at large scales.
    """
    _require_int("limit", limit, 4)
    _require_int("target_m", target_m, 0)
    if limit > MAX_UNION_LIMIT:
        raise ValueError(
            "limit exceeds MAX_UNION_LIMIT; use asymptotic_union_bound instead"
        )
    total = 0.0
    bound = math.isqrt(2 * limit)
    for p in primes_up_to(bound):
        if p <= target_m:
            continue
        available = int(math.log(2 * limit) / math.log(p)) - 1
        if available < 1:
            continue
        total += (1.0 / p) * 2.0 ** (-available)
    return {
        "limit": limit,
        "modelled_level_failure_sum": total,
        "log_two": math.log(2),
        "union_bound_total": math.log(2) + total,
        "union_bound_exceeds_one": math.log(2) + total > 1.0,
    }


def asymptotic_union_bound(term_count: int = 400) -> dict[str, Any]:
    """Closed-form limit of the union bound, by Mertens grouping.

    A prime with exactly ``r`` usable levels satisfies
    ``X^{1/(r+2)} < p <= X^{1/(r+1)}``, and Mertens gives
    ``sum 1/p = log((r+2)/(r+1))`` over that range, so the modelled level mass
    tends to ``sum_r 2^{-r} log((r+2)/(r+1))``. Fixed small primes contribute
    nothing in the limit because their level counts grow with ``X``.
    """
    _require_int("term_count", term_count, 1)
    level_sum = sum(
        2.0 ** (-r) * math.log((r + 2) / (r + 1)) for r in range(1, term_count + 1)
    )
    total = math.log(2) + level_sum
    return {
        "term_count": term_count,
        "level_failure_sum": level_sum,
        "union_bound_total": total,
        "union_bound_exceeds_one": total > 1.0,
    }


def analyze(*, blocks: tuple[tuple[int, int], ...], target_m: int) -> dict[str, Any]:
    rows = [decompose_block(start, count, target_m) for start, count in blocks]
    for row in rows:
        if row["observed_over_predicted_level_only"] is None:
            raise AssertionError("degenerate block")
    asymptotic = asymptotic_union_bound()
    if not asymptotic["union_bound_exceeds_one"]:
        raise AssertionError(
            "the asymptotic union bound no longer exceeds 1; recheck the model"
        )
    if min(row["observed_over_predicted_level_only"] for row in rows) <= 1.0:
        raise AssertionError(
            "independence would not be refuted on some block; widen the blocks"
        )
    return {
        "schema_version": 1,
        "status": "PASS_EXACT_GOOD_DENSITY_DECOMPOSITION",
        "claim_boundary": (
            "Each block is classified exactly and exhaustively; the four "
            "classes partition it and the partition is asserted. Densities are "
            "exact counts over the declared blocks, not asymptotic claims and "
            "not samples. The independence prediction is a labelled heuristic "
            "shown here only to be refuted by the exact counts. Nothing here "
            "proves or disproves positive density for compensation-good "
            "integers; it removes one candidate proof strategy."
        ),
        "parameters": {"blocks": [list(block) for block in blocks], "target_m": target_m},
        "blocks": rows,
        "union_bound": [
            union_bound_terms(start + count, target_m) for start, count in blocks
        ],
        "asymptotic_union_bound": asymptotic,
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "good_density_decomposition.py": sha256_file(Path(__file__)),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--block-starts", type=int, nargs="+", default=[10**6, 10**7])
    parser.add_argument("--block-count", type=int, default=200_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        blocks=tuple((start, args.block_count) for start in args.block_starts),
        target_m=args.target_m,
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
                "blocks": [
                    {
                        "start": row["start"],
                        "good": row["densities"]["good"],
                        "short_only": row["densities"]["short_only"],
                        "level_only": row["densities"]["level_only"],
                        "both": row["densities"]["both"],
                        "independence_ratio": row[
                            "observed_over_predicted_level_only"
                        ],
                    }
                    for row in result["blocks"]
                ],
                "union_bound_totals": [
                    row["union_bound_total"] for row in result["union_bound"]
                ],
                "asymptotic_union_bound_total": result["asymptotic_union_bound"][
                    "union_bound_total"
                ],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
