#!/usr/bin/env python3
"""Reduce every natural shift of a witness to its adjoined bad-window term.

Goodness at ``m`` constrains only the primes above ``m``, so it is monotone in
``m``. Combined with the window identity of Section 10 this localizes the whole
large-prime cost of a natural shift: an odd source keeps its window unchanged
and can fail only at primes ``p <= m+1``, while an even source adjoins the
single integer ``x/2`` and fails at a large prime exactly when that one integer
is not ``(m+1)``-compensation-good.

The producer checks both statements exhaustively on a bounded rectangle and at
every published witness of OEIS A375071.
"""

from __future__ import annotations

import argparse
import csv
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
)

ROOT = Path(__file__).resolve().parent
DEFAULT_WITNESS_CSV = ROOT / "data" / "oeis_a375071.csv"
DEFAULT_OUTPUT = ROOT / "data" / "parity_shift_m1_20_k4000_published.json"

Factorizer = TrialFactorizer | SPFFactorizer


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_int(name: str, value: int, minimum: int) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def bad_window(m: int, k: int) -> tuple[int, ...]:
    """Return the exact bad window of ``(m, k)``."""
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    return tuple(range(k + m // 2 + 1, k + m + 1))


def shift_window_difference(m: int, k: int) -> tuple[int, ...]:
    """Return the terms the natural shift adjoins to the source bad window.

    The result is empty for odd ``m`` and the single term ``x/2`` for even
    ``m``. The function proves the inclusion it reports.
    """
    _require_int("m", m, 1)
    _require_int("k", k, 2)
    source = bad_window(m, k)
    target = bad_window(m + 1, k - 1)
    if not set(source) <= set(target):
        raise AssertionError("the natural shift dropped a bad-window term")
    adjoined = tuple(term for term in target if term not in set(source))
    if m % 2:
        if adjoined:
            raise AssertionError("an odd source changed its bad window")
    else:
        if adjoined != (k + m // 2,):
            raise AssertionError("an even source adjoined an unexpected term")
        if 2 * adjoined[0] != m + 2 * k:
            raise AssertionError("the adjoined term is not x/2")
    return adjoined


def small_prime_slacks(m: int, k: int) -> tuple[dict[str, int], ...]:
    """Return the exact slack of every prime ``p <= m``."""
    _require_int("m", m, 1)
    _require_int("k", k, 1)
    return tuple(
        {"prime": p, "slack": slack(m, k, p)} for p in primes_up_to(m)
    )


def classify_shift(
    m: int,
    k: int,
    factorizer: Factorizer,
) -> dict[str, Any]:
    """Classify one natural shift of a witness by its adjoined term only."""
    adjoined = shift_window_difference(m, k)
    source_obstructions = large_prime_window_obstructions(m, k, factorizer)
    if source_obstructions:
        raise ValueError("classify_shift requires a large-prime-good source")
    target_m, target_k = m + 1, k - 1
    inherited = bad_window(m, k)
    inherited_records = [
        record
        for record in large_prime_window_records(m, k, factorizer)
        if record["prime"] > target_m
    ]
    for record in inherited_records:
        if record["predicted_slack"] < 0:
            raise AssertionError("goodness failed to be monotone in m")
    target_obstructions = large_prime_window_obstructions(
        target_m, target_k, factorizer
    )
    adjoined_obstructions = [
        record
        for record in target_obstructions
        if record["term"] in set(adjoined)
    ]
    if len(adjoined_obstructions) != len(target_obstructions):
        raise AssertionError("a large-prime obstruction survived from the source")
    target_small = small_prime_slacks(target_m, target_k)
    small_ok = all(row["slack"] >= 0 for row in target_small)
    large_ok = not target_obstructions
    if m % 2 and not large_ok:
        raise AssertionError("an odd source produced a large-prime obstruction")
    return {
        "source_m": m,
        "source_k": k,
        "target_m": target_m,
        "target_k": target_k,
        "inherited_term_count": len(inherited),
        "inherited_large_prime_checks": len(inherited_records),
        "adjoined_terms": list(adjoined),
        "adjoined_obstructions": adjoined_obstructions,
        "target_small_prime_slacks": list(target_small),
        "small_prime_tier_ok": small_ok,
        "large_prime_tier_ok": large_ok,
        "shift_is_witness": small_ok and large_ok,
        "parity": "odd_source" if m % 2 else "even_source",
    }


def verify_rectangle(max_m: int, max_k: int) -> dict[str, Any]:
    """Check the reduction at every large-prime-good pair of a rectangle."""
    _require_int("max_m", max_m, 2)
    _require_int("max_k", max_k, 2)
    factorizer = SPFFactorizer(max_m + 1 + max_k)
    counts = Counter()
    category_counts = Counter()
    for m in range(1, max_m + 1):
        for k in range(2, max_k + 1):
            counts["scanned_pairs"] += 1
            if large_prime_window_obstructions(m, k, factorizer):
                continue
            if any(row["slack"] < 0 for row in small_prime_slacks(m, k)):
                continue
            counts["witness_sources"] += 1
            row = classify_shift(m, k, factorizer)
            counts["adjoined_terms"] += len(row["adjoined_terms"])
            counts["inherited_large_prime_checks"] += row[
                "inherited_large_prime_checks"
            ]
            if row["parity"] == "odd_source":
                if row["adjoined_terms"]:
                    raise AssertionError("odd source adjoined a term")
                category_counts[
                    "odd_source_shift_is_witness"
                    if row["shift_is_witness"]
                    else "odd_source_small_prime_failure"
                ] += 1
            else:
                if len(row["adjoined_terms"]) != 1:
                    raise AssertionError("even source adjoined the wrong count")
                if row["large_prime_tier_ok"] != (not row["adjoined_obstructions"]):
                    raise AssertionError("even-source large tier misclassified")
                if row["shift_is_witness"]:
                    category_counts["even_source_shift_is_witness"] += 1
                elif not row["adjoined_obstructions"]:
                    category_counts["even_source_small_prime_failure"] += 1
                elif row["small_prime_tier_ok"]:
                    category_counts["even_source_adjoined_term_failure"] += 1
                else:
                    category_counts["even_source_both_tiers_failure"] += 1
    return {
        "counts": dict(sorted(counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
    }


def verify_published(witness_csv: Path) -> dict[str, Any]:
    """Classify the natural shift of every published witness."""
    rows = [
        (int(row["m"]), int(row["k"]))
        for row in csv.DictReader(
            witness_csv.read_text(encoding="utf-8").splitlines()
        )
        if int(row["m"]) >= 1
    ]
    published = dict(rows)
    maximum_x = max(m + 1 + 2 * (k + m) for m, k in rows)
    if maximum_x > MAX_TRIAL_FACTORIZER_VALUE:
        raise ValueError("published shifts exceed the retained factorizer bound")
    factorizer = TrialFactorizer(maximum_x)
    counts = Counter()
    category_counts = Counter()
    records: list[dict[str, Any]] = []
    for m, k in rows:
        row = classify_shift(m, k, factorizer)
        counts["published_shifts"] += 1
        counts["inherited_large_prime_checks"] += row[
            "inherited_large_prime_checks"
        ]
        successor = published.get(m + 1)
        row["published_successor_k"] = successor
        row["successor_is_natural_shift"] = successor == k - 1
        if row["successor_is_natural_shift"] and not row["shift_is_witness"]:
            raise AssertionError(
                "the published successor equals a shift classified as failing"
            )
        if row["parity"] == "even_source":
            category_counts[
                "even_source_shift_is_witness"
                if row["shift_is_witness"]
                else "even_source_shift_fails"
            ] += 1
        else:
            category_counts[
                "odd_source_shift_is_witness"
                if row["shift_is_witness"]
                else "odd_source_shift_fails"
            ] += 1
        records.append(row)
    natural_steps = sum(
        1 for row in records if row["successor_is_natural_shift"]
    )
    return {
        "counts": dict(sorted(counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "published_successor_is_natural_shift": natural_steps,
        "records": records,
    }


def analyze(
    *,
    max_m: int,
    max_k: int,
    witness_csv: Path,
) -> dict[str, Any]:
    witness_csv = witness_csv.resolve()
    rectangle = verify_rectangle(max_m, max_k)
    published = verify_published(witness_csv)
    return {
        "schema_version": 1,
        "status": "PASS_EXACT_PARITY_SHIFT_REDUCTION",
        "claim_boundary": (
            "Every classification here is exact and applies only to natural "
            "shifts of sources that already pass both tiers. The rectangle is "
            "exhausted as declared, and every published witness is checked with "
            "complete factorizations of both bad windows inside the retained "
            "bound. No statement is made about non-natural repairs."
        ),
        "parameters": {
            "rectangle_max_m": max_m,
            "rectangle_max_k": max_k,
            "shift": "(m, k) -> (m + 1, k - 1)",
            "adjoined_term": "x / 2 = k + m / 2 for even m, none for odd m",
        },
        "rectangle": rectangle,
        "published": published,
        "input": {
            "witness_csv_logical_name": witness_csv.name,
            "witness_csv_sha256": sha256_file(witness_csv),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "processes": 1,
        },
        "implementation_sha256": {
            "parity_shift_certificate.py": sha256_file(Path(__file__)),
            "compensation_structure_analysis.py": sha256_file(
                ROOT / "compensation_structure_analysis.py"
            ),
            "erdos389.py": sha256_file(ROOT / "erdos389.py"),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-m", type=int, default=20)
    parser.add_argument("--max-k", type=int, default=4_000)
    parser.add_argument("--witness-csv", type=Path, default=DEFAULT_WITNESS_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    result = analyze(
        max_m=args.max_m,
        max_k=args.max_k,
        witness_csv=args.witness_csv,
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
                "rectangle_counts": result["rectangle"]["counts"],
                "rectangle_categories": result["rectangle"]["category_counts"],
                "published_counts": result["published"]["counts"],
                "published_categories": result["published"]["category_counts"],
                "published_natural_steps": result["published"][
                    "published_successor_is_natural_shift"
                ],
                "wall_seconds": result["resources"]["wall_seconds"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
