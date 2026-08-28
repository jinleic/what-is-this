#!/usr/bin/env python3
"""Exactly re-rank the negative families of a block-symmetric census.

The census searches (`search_n8_block_symmetric.py`) evaluate ``A_+`` in float64
and say so: their objective is discovery evidence, never a certificate.  This
module re-evaluates candidates in exact rational arithmetic using the same
primitives as `verify_gate_b_rational.py` -- the rational ``Q`` enclosure, the
one-sided Bellman recursion over the *exact* feasible action interval with no
relaxation and no clamp classification -- and records a two-sided enclosure per
family.

Why re-ranking is necessary, not merely tidy
--------------------------------------------
For the `n=8` `k=2` class the float and exact objectives disagree by far more
than float64 noise on families whose value is close to zero: the lowest-defect
witness reads ``-0.00285945`` in float64 and ``-0.00203238`` exactly, a 40%
relative gap.  A float ordering of near-zero values is therefore not a reliable
ordering of the exact ones, so the census's own "largest ratio" and "lowest
defect" picks have to be re-derived before any of them can carry a label.

Scope of each record
--------------------
Exact and complete: family size, distinct rows, per-coordinate degrees against
the cap, incidence against the integer Reimer threshold, the missing-join count,
activity, separation, the closure defect, the automorphism group, the order
orbits, and the ``A_+`` enclosure.  Nothing here is sampled and nothing is
float64.

Separation is reported, not required.  `DEFINITIONS.md` records normalization
(active and separating) as a search and reporting condition rather than a
condition in the displayed supremum, so a non-separating family is admissible
and counts; the flag is carried so that any claim can state which convention it
uses.

Resumability: one append-only JSONL record per family mask, keyed by the exact
row tuple.  A rerun skips finished families and reports zero new evaluations.
Records written under a different schema are counted and ignored, never trusted.

Run with one thread and low priority:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/certify_class_negatives.py \
        --census uc/gate_b/experiments/n8_k2_checkpoint.jsonl --threshold 0.0
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_gate_b_rational import (  # noqa: E402
    ACCUMULATOR_BITS,
    ALPHA,
    bellman_bounds,
    decimal_lower,
    decimal_upper,
    log2_lower,
    log2_upper,
    order_orbits,
    round_down,
    round_up,
    shapley_iid_bounds,
)
from bound_local_regime import (  # noqa: E402
    admissible,
    cap,
    defect,
    degrees,
    incidence,
    reimer_threshold,
    size_ceiling_holds,
    union_growth_floor,
)

RECORD_SCHEMA = 1
DEFAULT_CENSUS = HERE / "experiments" / "n8_k2_checkpoint.jsonl"
DEFAULT_CHECKPOINT = HERE / "experiments" / "n8_k2_exact_checkpoint.jsonl"
DEFAULT_REPORT = HERE / "candidates" / "n8_k2_exact_frontier.json"

# The two published comparison points this module is meant to move.
PUBLISHED_RATIO = Fraction(286_491, 10_000_000) / Fraction(64, 81)
PUBLISHED_FRONTIER_DEFECT = Fraction(1336, 2025)


def exact_a_plus(
    family: tuple[int, ...], dimension: int
) -> tuple[Fraction, Fraction, dict[str, object]]:
    """Two-sided rational enclosure of ``A_+`` over the exact feasible set."""
    representatives, group = order_orbits(family, dimension)
    values = [bellman_bounds(family, dimension, order) for order in representatives]
    count = len(values)
    average_low = round_down(
        sum((low for low, _high, _s in values), Fraction(0)) / count, ACCUMULATOR_BITS
    )
    average_high = round_up(
        sum((high for _low, high, _s in values), Fraction(0)) / count, ACCUMULATOR_BITS
    )
    iid_low, iid_high = shapley_iid_bounds(family, dimension)
    size = Fraction(len(family))
    low = round_down(
        (1 - ALPHA) * iid_low + ALPHA * average_low - log2_upper(size),
        ACCUMULATOR_BITS,
    )
    high = round_up(
        (1 - ALPHA) * iid_high + ALPHA * average_high - log2_lower(size),
        ACCUMULATOR_BITS,
    )
    assert low <= high
    distinct = Counter((a, b) for a, b, _s in values)
    detail = {
        "all_order_count": factorial(dimension),
        "order_orbit_count": len(representatives),
        "automorphism_count": len(group),
        "distinct_enclosure_count": len(distinct),
        "distinct_enclosure_multiplicities": sorted(distinct.values()),
        "bellman_state_count": sum(s for _a, _b, s in values),
        "shapley_iid_lower": str(iid_low),
        "shapley_iid_upper": str(iid_high),
        "bellman_average_lower": str(average_low),
        "bellman_average_upper": str(average_high),
    }
    return low, high, detail


def family_facts(family: tuple[int, ...], dimension: int) -> dict[str, object]:
    size = len(family)
    columns = [
        tuple((row >> coordinate) & 1 for row in family)
        for coordinate in range(dimension)
    ]
    observed = defect(family)
    return {
        "family_size": size,
        "coordinate_counts": list(degrees(family, dimension)),
        "cap": cap(size),
        "total_incidence": incidence(family),
        "reimer_threshold": reimer_threshold(size),
        "closure_defect": str(observed),
        "closure_defect_decimal": decimal_upper(observed, 12),
        "missing_ordered_join_pairs": int(observed * size * size),
        "largest_set": max(bin(row).count("1") for row in family),
        "dominant_set_threshold": str(Fraction(8, 5) * Fraction(incidence(family), size)),
        "union_growth_floor": (
            None
            if union_growth_floor(family, dimension) is None
            else str(union_growth_floor(family, dimension))
        ),
        "size_ceiling_holds": size_ceiling_holds(size, dimension),
        "active": all(any(column) for column in columns),
        "separating": len(set(columns)) == dimension,
        "admissible": admissible(family, dimension),
    }


def read_census(path: Path) -> list[dict]:
    records = []
    for line in path.open():
        line = line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            continue
        record = payload.get("record")
        if record is not None:
            records.append(record)
    return records


def read_checkpoint(path: Path) -> tuple[dict[tuple[int, ...], dict], int, int]:
    done: dict[tuple[int, ...], dict] = {}
    foreign = truncated = 0
    if not path.exists():
        return done, foreign, truncated
    for line in path.open():
        line = line.strip()
        if not line:
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError:
            truncated += 1
            continue
        if payload.get("schema") != RECORD_SCHEMA:
            foreign += 1
            continue
        done[tuple(payload["family_rows"])] = payload
    return done, foreign, truncated


def append_record(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--dimension", type=int, default=8)
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Re-evaluate families whose float A_+ is below this value.",
    )
    parser.add_argument("--max-evaluations", type=int, default=0)
    args = parser.parse_args()

    census = read_census(args.census)
    targets = []
    seen: set[tuple[int, ...]] = set()
    for record in census:
        if record.get("a_plus") is None or record["a_plus"] >= args.threshold:
            continue
        rows = tuple(sorted(record["family_rows"]))
        if rows in seen:
            continue
        seen.add(rows)
        targets.append((rows, record["a_plus"]))
    targets.sort(key=lambda item: item[1])

    done, foreign, truncated = read_checkpoint(args.checkpoint)
    evaluated = 0
    started = time.monotonic()
    for rows, float_value in targets:
        if rows in done:
            continue
        if args.max_evaluations and evaluated >= args.max_evaluations:
            break
        low, high, detail = exact_a_plus(rows, args.dimension)
        facts = family_facts(rows, args.dimension)
        observed = Fraction(facts["closure_defect"])
        payload = {
            "schema": RECORD_SCHEMA,
            "dimension": args.dimension,
            "family_rows": list(rows),
            "float_a_plus": float_value,
            "a_plus_lower": str(low),
            "a_plus_upper": str(high),
            "a_plus_upper_decimal": decimal_upper(high, 12),
            "a_plus_enclosure_width_decimal": decimal_upper(high - low, 12),
            "certified_negative": high < 0,
            "exact_ratio_lower": str(-high / observed) if observed else None,
            "exact_ratio_lower_decimal": decimal_lower(-high / observed, 12) if observed else None,
            "beats_published_ratio": bool(observed and -high / observed > PUBLISHED_RATIO),
            "beats_published_frontier_defect": observed < PUBLISHED_FRONTIER_DEFECT,
            "float_exact_gap_decimal": decimal_upper(abs(Fraction(high) - Fraction(float_value).limit_denominator(10**9)), 12),
            "order_evaluation": detail,
            "exact_facts": facts,
        }
        append_record(args.checkpoint, payload)
        done[rows] = payload
        evaluated += 1
        print(
            f"[{evaluated}/{len(targets)}] m={facts['family_size']} "
            f"defect={facts['closure_defect']} float={float_value:+.8f} "
            f"exact<={payload['a_plus_upper_decimal']} "
            f"neg={payload['certified_negative']} "
            f"ratio>={payload['exact_ratio_lower_decimal']} "
            f"({time.monotonic() - started:.0f}s)",
            file=sys.stderr,
            flush=True,
        )

    certified = [p for p in done.values() if p["certified_negative"]]
    frontier = min(
        certified, key=lambda p: Fraction(p["exact_facts"]["closure_defect"]), default=None
    )
    sharpest = max(
        certified, key=lambda p: Fraction(p["exact_ratio_lower"]), default=None
    )
    report = {
        "census": str(args.census),
        "dimension": args.dimension,
        "candidates_considered": len(targets),
        "records": len(done),
        # `new_evaluations` is deliberately NOT recorded here: it is 21 on a
        # first run and 0 on a resume, which would make the artifact fail its
        # own hash.  It is printed to stdout instead.
        "foreign_schema_lines": foreign,
        "truncated_lines": truncated,
        "certified_negative_count": len(certified),
        # A float-negative family that is exactly non-negative is a false
        # negative of the census objective.  A float-non-negative family that is
        # exactly negative would be a *missed* witness.  The two are recorded
        # separately because only the first can be produced by rounding of a
        # true negative, and only the second would move the frontier.
        "float_negative_but_exactly_nonnegative": sorted(
            p["exact_facts"]["closure_defect"]
            for p in done.values()
            if p["float_a_plus"] < 0 and not p["certified_negative"]
        ),
        "float_nonnegative_but_exactly_negative": sorted(
            p["exact_facts"]["closure_defect"]
            for p in done.values()
            if p["float_a_plus"] >= 0 and p["certified_negative"]
        ),
        "float_negative_count": sum(1 for p in done.values() if p["float_a_plus"] < 0),
        "largest_float_exact_gap_decimal": max(
            (p["float_exact_gap_decimal"] for p in done.values()), default=None
        ),
        "published_ratio": str(PUBLISHED_RATIO),
        "published_frontier_defect": str(PUBLISHED_FRONTIER_DEFECT),
        "lowest_defect_certified_negative": frontier,
        "highest_ratio_certified_negative": sharpest,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(
        json.dumps(
            {
                "records": len(done),
                "new_evaluations": evaluated,
                "certified_negative": len(certified),
                "foreign_schema_lines": foreign,
                "truncated_lines": truncated,
            }
        )
    )
    if frontier is not None:
        print(
            f"lowest certified defect: {frontier['exact_facts']['closure_defect']} "
            f"(A_+ <= {frontier['a_plus_upper_decimal']}, "
            f"separating={frontier['exact_facts']['separating']}, "
            f"beats published frontier={frontier['beats_published_frontier_defect']})"
        )
    if sharpest is not None:
        print(
            f"highest certified ratio: {sharpest['exact_ratio_lower_decimal']} "
            f"at defect {sharpest['exact_facts']['closure_defect']} "
            f"(separating={sharpest['exact_facts']['separating']}, "
            f"beats published ratio={sharpest['beats_published_ratio']})"
        )


if __name__ == "__main__":
    main()
