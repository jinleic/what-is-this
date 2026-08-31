#!/usr/bin/env python3
"""Chosen-coordinate cloning: exact finite ladder for the Gate B objective.

Duplicating a coordinate is defect-free for every coordinate, but its effect on
``A_+`` is **not sign-definite**.  This module measures the useful coordinate 2
of ``n8tiny``; ``audit_clone_limit.py`` gives the general order-reweighting
theorem and exact limit.

What cloning preserves, exactly
-------------------------------
Let ``F' `` be ``F`` with coordinate ``i`` duplicated.  Every row of ``F'`` is a
row of ``F`` with one bit repeated, so the map ``F -> F'`` is a bijection that
commutes with union.  Therefore

* the size ``m`` is unchanged;
* the closure defect ``eps_vee`` is unchanged, since ``X | Y`` lands in the
  family before cloning exactly when it does after;
* every coordinate degree is unchanged and the clone's degree equals its twin's,
  so the cap ``floor(2m/5)`` is still met;
* the incidence rises by ``deg_i`` while the Reimer threshold ``R_m`` depends
  only on ``m`` and so does not move.

Hence **cloning preserves admissibility unconditionally** and cannot change the
defect.  It is free in the only two quantities the local question measures.

What cloning changes
--------------------
``log2 m`` is fixed, but the direction of ``A_+`` depends on which fixed-order
objectives the clone multiplicity upweights.  On the ``n8tiny`` core, cloning
coordinate 2 lowers ``A_+``; cloning any of coordinates 0, 1, 3, or 4 instead
raises it from ``+0.010695694`` to ``+0.012500067``.  The earlier claim that
cloning always lowers both ``Q`` and ``C_+`` was false.  The exact statement is:
cloning only reweights the original fixed-order values.

That is exactly how the record low-defect witness exists.  The registered base
``n8tiny`` has four identical coordinates; collapsing them leaves an admissible
five-coordinate family with the same size 15 and the same defect ``4/9``, but

    A_+(collapsed, n=5) = +0.010695694     (positive)
    A_+(n8tiny,    n=8) = -0.000478465     (negative)

so the clones, not the ground set, are what make it a witness.  By contrast
``n8clone_lo`` and ``n8clone_hi`` collapse to *inadmissible* seven-coordinate
families (incidence 196 below ``R_70 = 215``), so for those two the clone buys
feasibility instead.  Both effects are real and they are different.

The ceiling
-----------
For the chosen coordinate, the finite deltas shrink and the objective converges
to the exact first-rank average ``-0.002338401508...`` certified by
``audit_clone_limit.py``.  The convergence is **not geometric**: with total
multiplicity ``r``, the probability that another coordinate precedes the first
clone is ``(n-1)/(n+r-1)=O(1/r)``.  The observed delta ratios therefore cannot
prove geometric decay; they quantify only the displayed finite ladder.  The
total available gain is bounded exactly by the convex hull of the original
fixed-order objectives.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/audit_clone_saturation.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from verify_gate_b_rational import (  # noqa: E402
    ACCUMULATOR_BITS,
    ALPHA,
    BASES,
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
)

DEFAULT_REPORT = HERE / "candidates" / "clone_saturation_audit.json"


def collapse(rows: tuple[int, ...], dimension: int) -> tuple[tuple[int, ...], list[int]]:
    """Drop duplicate columns; returns the family on the distinct coordinates."""
    columns = [
        tuple((row >> coordinate) & 1 for row in rows)
        for coordinate in range(dimension)
    ]
    keep: list[int] = []
    seen: set[tuple[int, ...]] = set()
    for coordinate, column in enumerate(columns):
        if column not in seen:
            seen.add(column)
            keep.append(coordinate)
    collapsed = tuple(
        sorted(
            sum(((row >> coordinate) & 1) << index for index, coordinate in enumerate(keep))
            for row in rows
        )
    )
    return collapsed, keep


def clone(
    rows: tuple[int, ...], dimension: int, coordinate: int, copies: int
) -> tuple[tuple[int, ...], int]:
    out = []
    for row in rows:
        bit = (row >> coordinate) & 1
        extra = 0
        for offset in range(copies):
            extra |= bit << (dimension + offset)
        out.append(row | extra)
    return tuple(sorted(out)), dimension + copies


def exact_a_plus(rows: tuple[int, ...], dimension: int) -> dict[str, object]:
    representatives, group = order_orbits(rows, dimension)
    values = [bellman_bounds(rows, dimension, order) for order in representatives]
    low = round_down(
        sum(value[0] for value in values) / len(values), ACCUMULATOR_BITS
    )
    high = round_up(
        sum(value[1] for value in values) / len(values), ACCUMULATOR_BITS
    )
    iid_low, iid_high = shapley_iid_bounds(rows, dimension)
    size = Fraction(len(rows))
    a_low = round_down(
        (1 - ALPHA) * iid_low + ALPHA * low - log2_upper(size), ACCUMULATOR_BITS
    )
    a_high = round_up(
        (1 - ALPHA) * iid_high + ALPHA * high - log2_lower(size), ACCUMULATOR_BITS
    )
    return {
        "a_plus_lower": str(a_low),
        "a_plus_upper": str(a_high),
        "a_plus_upper_decimal": decimal_upper(a_high, 12),
        "certified_negative": a_high < 0,
        "shapley_iid_upper_decimal": decimal_upper(iid_high, 12),
        "bellman_average_upper_decimal": decimal_upper(high, 12),
        "automorphism_count": len(group),
        "order_orbit_count": len(representatives),
        "_a_plus_upper": a_high,
    }


def structural(rows: tuple[int, ...], dimension: int) -> dict[str, object]:
    return {
        "dimension": dimension,
        "family_size": len(rows),
        "closure_defect": str(defect(rows)),
        "coordinate_counts": list(degrees(rows, dimension)),
        "cap": cap(len(rows)),
        "total_incidence": incidence(rows),
        "reimer_threshold": reimer_threshold(len(rows)),
        "admissible": admissible(rows, dimension),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=str, default="n8tiny")
    parser.add_argument("--coordinate", type=int, default=2)
    parser.add_argument("--max-clones", type=int, default=5)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    registered = BASES[args.base]
    full = registered.reconstruct()
    core, kept = collapse(full, registered.dimension)

    ladder = []
    previous = None
    for copies in range(args.max_clones + 1):
        rows, dimension = (
            (core, len(kept)) if copies == 0 else clone(core, len(kept), args.coordinate, copies)
        )
        record = structural(rows, dimension)
        record["clones"] = copies
        record.update(exact_a_plus(rows, dimension))
        current = record.pop("_a_plus_upper")
        record["delta_from_previous_decimal"] = (
            None if previous is None else decimal_lower(current - previous, 12)
        )
        previous = current
        ladder.append(record)

    defects = {row["closure_defect"] for row in ladder}
    sizes = {row["family_size"] for row in ladder}
    caps_met = all(row["admissible"] for row in ladder)
    signs = [row["certified_negative"] for row in ladder]
    flip = next(
        (row["clones"] for row in ladder if row["certified_negative"]), None
    )
    deltas = [
        Fraction(row["delta_from_previous_decimal"])
        for row in ladder
        if row["delta_from_previous_decimal"] is not None
    ]
    shrinking = all(
        abs(later) < abs(earlier) for earlier, later in zip(deltas, deltas[1:])
    )

    report = {
        "base": args.base,
        "cloned_coordinate": args.coordinate,
        "collapsed_dimension": len(kept),
        "collapsed_kept_coordinates": kept,
        "collapsed_rows": list(core),
        "ladder": ladder,
        "defect_is_invariant": len(defects) == 1,
        "size_is_invariant": len(sizes) == 1,
        "admissibility_preserved_throughout": caps_met,
        "sign_flips_at_clone_count": flip,
        "improvement_is_shrinking_over_reported_ladder": shrinking,
        "delta_ratios_decimal": [
            decimal_upper(abs(later / earlier), 6)
            for earlier, later in zip(deltas, deltas[1:])
        ],
        "convergence_status": (
            "FINITE LADDER ONLY; exact limit in clone_limit_audit.json; "
            "convergence is O(1/r), not geometric"
        ),
        "verdict": (
            "CHOSEN_CLONE_IS_DEFECT_FREE_AND_FINITE_LADDER_SHRINKS"
            if len(defects) == 1 and len(sizes) == 1 and caps_met and shrinking
            else "UNEXPECTED"
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(json.dumps({k: report[k] for k in (
        "base", "verdict", "collapsed_dimension", "defect_is_invariant",
        "size_is_invariant", "admissibility_preserved_throughout",
        "sign_flips_at_clone_count",
        "improvement_is_shrinking_over_reported_ladder",
        "delta_ratios_decimal")}, indent=2))
    print()
    print(f"{'clones':>6} {'n':>3} {'A_+ <=':>16} {'delta':>16} {'defect':>8} {'adm':>5}")
    for row in ladder:
        print(
            f"{row['clones']:6d} {row['dimension']:3d} "
            f"{row['a_plus_upper_decimal']:>16} "
            f"{str(row['delta_from_previous_decimal']):>16} "
            f"{row['closure_defect']:>8} {str(row['admissible']):>5}"
        )


if __name__ == "__main__":
    main()
