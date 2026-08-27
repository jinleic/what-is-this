#!/usr/bin/env python3
"""Certify the powersmooth necessity bound and the witness density ceiling.

Three exact statements are checked.  First, every compensation-good integer
obeys ``p**e * (p**e + 1) <= 2w`` at each of its prime powers above ``m``, so a
prime factor with ``p*p > 2w`` is always fatal.  Second, the count of integers
carrying such a fatal short-cofactor prime is computed by two independent exact
methods and compared with the ``log 2`` density limit that bounds the witness
count.  Third, the exact compensation-good density and run-length spectrum are
measured at several scales, including a deliberately biased control block that
starts on a published witness window.

The closing growth model is a numerical extrapolation, not a proof, and is
labelled as such in the artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import time
from bisect import bisect_right
from collections import Counter
from math import isqrt
from pathlib import Path
from typing import Any

from artifact_io import atomic_write_json
from compensation_run_search import segmented_compensation_obstructions
from erdos389 import SPFFactorizer, primes_up_to, slack
from residue_control_analysis import upper_half_prefix_count

ROOT = Path(__file__).resolve().parent
DEFAULT_WITNESS_CSV = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_RUN_ARTIFACTS = (
    ROOT / "data" / "compensation_run_m27_k50001_h100000000.json",
    ROOT / "data" / "compensation_run_m27_k100050001_h900000000.json",
)
DEFAULT_OUTPUT = ROOT / "data" / "smooth_density_m1_20_k3000.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def compensation_levels(term: int, p: int, exponent: int) -> int:
    """Return ``C_p(term)`` for the exact prime power ``p**exponent``."""
    _require_int("term", term, 1)
    _require_int("p", p, 2)
    _require_int("exponent", exponent, 1)
    return upper_half_prefix_count(term // p**exponent, p)


def verify_powersmooth_necessity(max_m: int, max_k: int) -> dict[str, Any]:
    """Check the exact minimal-term bound at every large window prime power."""
    _require_int("max_m", max_m, 1)
    _require_int("max_k", max_k, 1)
    factorizer = SPFFactorizer(max_m + max_k)
    counts = Counter()
    tight_rows: list[dict[str, int]] = []
    for m in range(1, max_m + 1):
        for k in range(1, max_k + 1):
            for position in range(m // 2 + 1, m + 1):
                term = k + position
                for p, exponent in factorizer.factor(term).items():
                    if p <= m:
                        continue
                    counts["prime_power_checks"] += 1
                    power = p**exponent
                    levels = compensation_levels(term, p, exponent)
                    good = levels >= exponent
                    minimal_bound_holds = power * (power + 1) <= 2 * term
                    if good and not minimal_bound_holds:
                        raise AssertionError(
                            "compensated prime power violated the minimal-term bound"
                        )
                    if good:
                        counts["good_prime_powers"] += 1
                        if power * (power + 1) == 2 * term:
                            counts["tight_minimal_terms"] += 1
                            if len(tight_rows) < 40:
                                tight_rows.append(
                                    {
                                        "m": m,
                                        "k": k,
                                        "term": term,
                                        "prime": p,
                                        "exponent": exponent,
                                    }
                                )
                    if power * power > 2 * term:
                        counts["prime_powers_above_minimal_bound"] += 1
                        if good:
                            raise AssertionError(
                                "a prime power above the minimal bound compensated"
                            )
                    if p * p > 2 * term and good:
                        raise AssertionError(
                            "a prime above sqrt(2w) was compensated"
                        )
    return {
        "counts": dict(sorted(counts.items())),
        "tight_minimal_term_rows": tight_rows,
    }


def count_short_cofactor_integers(limit: int) -> dict[str, Any]:
    """Count ``n <= limit`` having a prime ``p | n`` with ``p*p > 2n``.

    The divisor sum and a direct sieve are independent exact methods; the
    producer asserts that they agree.
    """
    _require_int("limit", limit, 2)
    primes = primes_up_to(limit)

    def prime_count(bound: int) -> int:
        if bound < 0:
            return 0
        if bound > limit:
            raise ValueError("prime counting bound exceeded the sieve")
        return bisect_right(primes, bound)

    divisor_total = 0
    cofactor = 1
    while 2 * cofactor * cofactor < limit:
        upper = limit // cofactor
        divisor_total += max(0, prime_count(upper) - prime_count(2 * cofactor))
        cofactor += 1

    marked = bytearray(limit + 1)
    for p in primes:
        multiple = p
        while multiple <= limit:
            if p * p > 2 * multiple:
                marked[multiple] = 1
            multiple += p
    sieve_total = sum(marked)

    if divisor_total != sieve_total:
        raise AssertionError("short-cofactor counts disagree between methods")
    density = divisor_total / limit
    return {
        "limit": limit,
        "short_cofactor_integers": divisor_total,
        "density": density,
        "log_two_reference": math.log(2),
        "signed_deviation_from_log_two": density - math.log(2),
        "passing_density": 1.0 - density,
    }


def count_narrow_repair_integers(limit: int, direct_limit: int) -> dict[str, Any]:
    """Count ``n <= limit`` whose dominant prime is repaired by its cofactor.

    A prime power ``p**e || n`` with cofactor ``c = n / p**e < p`` compensates
    only when ``e == 1`` and ``p < 2c``, that is when ``n < p*p < 2n``. The
    divisor enumeration is exact; a direct factorization repeats the count for
    ``limit <= direct_limit``.
    """
    _require_int("limit", limit, 2)
    _require_int("direct_limit", direct_limit, 0)
    primes = primes_up_to(limit)
    divisor_total = 0
    cofactor = 1
    while cofactor * cofactor < limit:
        upper = min(2 * cofactor - 1, limit // cofactor)
        if upper > cofactor:
            divisor_total += bisect_right(primes, upper) - bisect_right(
                primes, cofactor
            )
        cofactor += 1

    direct_total: int | None = None
    if limit <= direct_limit:
        factorizer = SPFFactorizer(limit)
        direct_total = 0
        for value in range(2, limit + 1):
            for p, exponent in factorizer.factor(value).items():
                if exponent != 1:
                    continue
                quotient = value // p
                if quotient < p < 2 * quotient:
                    direct_total += 1
                    break
        if direct_total != divisor_total:
            raise AssertionError("narrow-repair counts disagree between methods")
    return {
        "limit": limit,
        "narrow_repair_integers": divisor_total,
        "density": divisor_total / limit,
        "direct_recount": direct_total,
        "limit_over_two_log_limit": limit / (2.0 * math.log(limit)),
    }


def measure_compensation_block(
    *,
    target_m: int,
    start: int,
    count: int,
    label: str,
) -> dict[str, Any]:
    """Measure the exact compensation-good density and run spectrum."""
    _require_int("target_m", target_m, 1)
    _require_int("start", start, 1)
    _require_int("count", count, 1)
    trial_primes = primes_up_to(isqrt(start + count - 1))
    obstructions, _, _ = segmented_compensation_obstructions(
        m=target_m,
        start=start,
        count=count,
        trial_primes=trial_primes,
    )
    run_histogram = Counter()
    good = 0
    current = 0
    for prime in obstructions:
        if prime:
            if current:
                run_histogram[current] += 1
            current = 0
            continue
        good += 1
        current += 1
    if current:
        run_histogram[current] += 1
    density = good / count
    model = {
        run: count * (density**run) * ((1.0 - density) ** 2)
        for run in sorted(run_histogram)
    }
    longest = max(run_histogram) if run_histogram else 0
    return {
        "label": label,
        "target_m": target_m,
        "start": start,
        "count": count,
        "compensation_good_terms": good,
        "density": density,
        "longest_run": longest,
        "run_histogram": [
            {
                "run_length": run,
                "observed": run_histogram[run],
                "independence_model": model[run],
            }
            for run in sorted(run_histogram)
        ],
    }


def measure_small_prime_tier(target_m: int, start: int, count: int) -> dict[str, Any]:
    """Measure the exact density of the small-prime tier over consecutive ``k``."""
    _require_int("target_m", target_m, 1)
    _require_int("start", start, 1)
    _require_int("count", count, 1)
    primes = primes_up_to(target_m)
    passing = 0
    for k in range(start, start + count):
        if all(slack(target_m, k, p) >= 0 for p in primes):
            passing += 1
    return {
        "target_m": target_m,
        "start": start,
        "count": count,
        "small_prime_tier_passing": passing,
        "density": passing / count,
        "primes": primes,
    }


def _linear_fit(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    size = len(points)
    if size < 3:
        raise ValueError("need at least three points for a fit")
    sum_x = sum(x for x, _ in points)
    sum_y = sum(y for _, y in points)
    sum_xx = sum(x * x for x, _ in points)
    sum_xy = sum(x * y for x, y in points)
    denominator = size * sum_xx - sum_x * sum_x
    if denominator == 0:
        raise ValueError("degenerate fit")
    slope = (size * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y - slope * sum_x) / size
    mean_y = sum_y / size
    total = sum((y - mean_y) ** 2 for _, y in points)
    residual = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    r_squared = 1.0 - residual / total if total else 1.0
    return slope, intercept, r_squared


def growth_model(
    witness_csv: Path,
    target_m: int,
    measured_density: float,
    small_tier_density: float,
) -> dict[str, Any]:
    """Extrapolate the published least witnesses; explicitly not a proof."""
    rows = [
        (int(row["m"]), int(row["k"]))
        for row in csv.DictReader(
            witness_csv.read_text(encoding="utf-8").splitlines()
        )
        if int(row["m"]) >= 1
    ]
    if not rows:
        raise ValueError("witness transcription produced no rows")
    points_all = [((m + 1) // 2, math.log(k)) for m, k in rows]
    points_odd = [((m + 1) // 2, math.log(k)) for m, k in rows if m % 2 == 1]
    slope_all, intercept_all, r2_all = _linear_fit(points_all)
    slope_odd, intercept_odd, r2_odd = _linear_fit(points_odd)
    window_length = (target_m + 1) // 2
    independence_scale = 1.0 / (
        (measured_density**window_length) * small_tier_density
    )
    return {
        "claim": (
            "NOT PROVED. A least-squares extrapolation of the published least "
            "witnesses plus an independence model for compensation-good runs. "
            "It predicts a search scale; it neither proves nor disproves the "
            "existence of a witness."
        ),
        "regression_all_m": {
            "slope_per_window_term": slope_all,
            "intercept": intercept_all,
            "r_squared": r2_all,
            "implied_per_term_density": math.exp(-slope_all),
            "predicted_least_k": math.exp(
                intercept_all + slope_all * window_length
            ),
        },
        "regression_odd_m": {
            "slope_per_window_term": slope_odd,
            "intercept": intercept_odd,
            "r_squared": r2_odd,
            "implied_per_term_density": math.exp(-slope_odd),
            "predicted_least_k": math.exp(
                intercept_odd + slope_odd * window_length
            ),
        },
        "independence_model": {
            "measured_per_term_density": measured_density,
            "small_prime_tier_density": small_tier_density,
            "window_length": window_length,
            "predicted_least_k": independence_scale,
        },
    }


def recorded_run_statistics(paths: tuple[Path, ...]) -> list[dict[str, Any]]:
    """Read the exact run-search artifacts and recompute their densities."""
    rows: list[dict[str, Any]] = []
    for path in paths:
        payload = json.loads(path.read_bytes())
        if payload.get("status") != "EXACT_BOUNDED_COMPENSATION_RUN_SEARCH":
            raise ValueError(f"unexpected run-search status in {path.name}")
        counts = payload["counts"]
        good = counts["compensation_good_terms"]
        bad = counts["compensation_bad_terms"]
        total = good + bad
        density = good / total
        longest = counts["maximum_consecutive_compensation_good_terms"]
        predicted = math.log(total * (1.0 - density) ** 2) / math.log(1.0 / density)
        rows.append(
            {
                "logical_name": path.name,
                "sha256": sha256_file(path),
                "factored_integer_interval": payload["parameters"][
                    "factored_integer_interval"
                ],
                "terms": total,
                "compensation_good_terms": good,
                "density": density,
                "longest_run": longest,
                "independence_model_longest_run": predicted,
            }
        )
    return rows


def analyze(
    *,
    max_m: int,
    max_k: int,
    count_limits: tuple[int, ...],
    target_m: int,
    block_count: int,
    small_tier_count: int,
    witness_csv: Path,
    run_paths: tuple[Path, ...],
) -> dict[str, Any]:
    witness_csv = witness_csv.resolve()
    necessity = verify_powersmooth_necessity(max_m, max_k)
    density_rows = [count_short_cofactor_integers(limit) for limit in count_limits]
    narrow_rows = [
        count_narrow_repair_integers(limit, min(count_limits))
        for limit in count_limits
    ]
    witness_rows = [
        (int(row["m"]), int(row["k"]))
        for row in csv.DictReader(
            witness_csv.read_text(encoding="utf-8").splitlines()
        )
    ]
    published = {m: k for m, k in witness_rows}
    control_m = target_m - 2
    control_k = published[control_m]
    blocks = [
        measure_compensation_block(
            target_m=target_m,
            start=10**6 + 1,
            count=block_count,
            label="scale_1e6",
        ),
        measure_compensation_block(
            target_m=target_m,
            start=10**9 + 1,
            count=block_count,
            label="scale_1e9",
        ),
        measure_compensation_block(
            target_m=target_m,
            start=control_k + target_m + 1 + block_count,
            count=block_count,
            label="scale_5e12_unbiased",
        ),
        measure_compensation_block(
            target_m=target_m,
            start=control_k + control_m // 2 + 1,
            count=block_count,
            label="scale_5e12_published_window_control",
        ),
    ]
    control = next(
        block
        for block in blocks
        if block["label"] == "scale_5e12_published_window_control"
    )
    unbiased = next(
        block for block in blocks if block["label"] == "scale_5e12_unbiased"
    )
    window_length = (target_m + 1) // 2
    control_window_length = (control_m + 1) // 2
    if control["longest_run"] < control_window_length:
        raise AssertionError(
            "published witness window did not appear as a compensation-good run"
        )
    if control_window_length >= window_length:
        raise AssertionError("control window must be shorter than the target window")
    if control["longest_run"] <= unbiased["longest_run"]:
        raise AssertionError("published-window control failed to show its bias")
    small_tier = measure_small_prime_tier(
        target_m, control_k + block_count, small_tier_count
    )
    model = growth_model(
        witness_csv,
        target_m,
        unbiased["density"],
        small_tier["density"],
    )
    return {
        "schema_version": 1,
        "status": "PASS_EXACT_POWERSMOOTH_DENSITY_CERTIFICATE",
        "claim_boundary": (
            "The powersmooth necessity checks and the short-cofactor counts are "
            "exact on their declared finite domains, and the two counting methods "
            "are compared for every limit. Every measured density and run "
            "spectrum is an exact count over its declared block, not an "
            "asymptotic claim. The published-window control block is reported to "
            "document a sampling bias, not as a random sample. The closing growth "
            "model is a labelled numerical extrapolation and proves nothing."
        ),
        "parameters": {
            "necessity_max_m": max_m,
            "necessity_max_k": max_k,
            "short_cofactor_limits": list(count_limits),
            "target_m": target_m,
            "block_count": block_count,
            "small_prime_tier_count": small_tier_count,
            "window_length": window_length,
            "published_control_source_m": control_m,
            "published_control_source_k": control_k,
            "published_control_window_length": control_window_length,
        },
        "powersmooth_necessity": necessity,
        "short_cofactor_density": density_rows,
        "narrow_repair_density": narrow_rows,
        "compensation_blocks": blocks,
        "small_prime_tier": small_tier,
        "recorded_run_statistics": recorded_run_statistics(run_paths),
        "heuristic_growth_model": model,
        "input": {
            "witness_csv_logical_name": witness_csv.name,
            "witness_csv_sha256": sha256_file(witness_csv),
            "published_witness_rows": len(witness_rows),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "smooth_density_analysis.py": sha256_file(Path(__file__)),
            "compensation_run_search.py": sha256_file(
                ROOT / "compensation_run_search.py"
            ),
            "residue_control_analysis.py": sha256_file(
                ROOT / "residue_control_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=20)
    parser.add_argument("--max-k", type=int, default=3_000)
    parser.add_argument("--target-m", type=int, default=27)
    parser.add_argument("--block-count", type=int, default=250_000)
    parser.add_argument("--small-tier-count", type=int, default=20_000)
    parser.add_argument("--witness-csv", type=Path, default=DEFAULT_WITNESS_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--short-cofactor-limits",
        type=int,
        nargs="+",
        default=[100_000, 1_000_000, 4_000_000],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        max_m=args.max_m,
        max_k=args.max_k,
        count_limits=tuple(args.short_cofactor_limits),
        target_m=args.target_m,
        block_count=args.block_count,
        small_tier_count=args.small_tier_count,
        witness_csv=args.witness_csv,
        run_paths=DEFAULT_RUN_ARTIFACTS,
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
                "necessity_counts": result["powersmooth_necessity"]["counts"],
                "short_cofactor_density": [
                    {
                        "limit": row["limit"],
                        "density": row["density"],
                        "deviation": row["signed_deviation_from_log_two"],
                    }
                    for row in result["short_cofactor_density"]
                ],
                "narrow_repair_density": [
                    {"limit": row["limit"], "density": row["density"]}
                    for row in result["narrow_repair_density"]
                ],
                "block_densities": [
                    {
                        "label": block["label"],
                        "density": block["density"],
                        "longest_run": block["longest_run"],
                    }
                    for block in result["compensation_blocks"]
                ],
                "small_prime_tier_density": result["small_prime_tier"]["density"],
                "predicted_least_k": {
                    "regression_all_m": result["heuristic_growth_model"][
                        "regression_all_m"
                    ]["predicted_least_k"],
                    "independence_model": result["heuristic_growth_model"][
                        "independence_model"
                    ]["predicted_least_k"],
                },
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
