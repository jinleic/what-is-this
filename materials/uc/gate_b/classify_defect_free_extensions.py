#!/usr/bin/env python3
"""Every defect-free coordinate extension, classified -- and the record core audited.

The question this settles
-------------------------
Coordinate cloning (Proposition 18) improves ``A_+`` at exactly zero defect cost,
so the obvious next move is to ask what *else* does.  This module answers it
completely: it characterises **every** way to add a coordinate without moving the
closure defect, and then measures what each one does to the objective.

The classification
------------------
A coordinate extension of ``F`` is determined by the set of rows that receive a 1
in the new column.  Writing ``U`` for that subset and

    phi(A) = A u {new}  if A in U,      phi(A) = A  otherwise,

``phi`` is a bijection onto the extended family, so the size never moves.  For the
defect, ``phi(X) u phi(Y)`` has base part ``X u Y`` and flag ``[X in U] or
[Y in U]``, and it lands in the image exactly when ``X u Y`` is in ``F`` *and* the
flag agrees with ``[X u Y in U]``.  Hence:

    the defect never decreases, and it is preserved exactly when
        for all X, Y with X u Y in F:  [X u Y in U] = [X in U] or [Y in U].     (*)

Call such a ``U`` **join-consistent**.  Taking ``Y`` inside ``X`` in (*) shows a
join-consistent ``U`` is an up-set of ``(F, subset)``; the remaining content is
that a successful join landing in ``U`` must already have a part in ``U``.  Both
halves are audited here.

Three facts follow, all checked below:

* **Admissibility is automatic.**  Old degrees are untouched and the new degree is
  ``|U|``, so the cap needs only ``|U| <= floor(2m/5)``; the incidence rises by
  ``|U|`` while ``R_m`` depends only on ``m``, so Reimer cannot break.
* **The OR-columns are the nonconstant Boolean join-homomorphisms.**  A new
  column defined by a global Boolean function ``g`` preserves joins for every
  family iff ``g(x or y) = g(x) or g(y)``.  Apart from the two constant
  functions, these are precisely ``g(x) = OR_{i in S} x_i`` for nonempty ``S``;
  cloning is ``|S| = 1``.  Constant 0 is inert and constant 1 violates the cap.
* **The family-specific class is strictly larger.**  Because (*) only constrains
  *successful* joins, a defective family admits extensions that are not induced
  by a global join-homomorphism.  The collapsed core of ``n8tiny`` has exactly
  one such usable column, ``U* = {A in F : |A| >= 3}`` of size 5, while every
  nonconstant OR-column there has size at least 6.

And why it does not help
------------------------
Measured exactly, ``U*`` *raises* the objective from ``+0.010695694`` to
``+0.114769081`` -- ten times worse than doing nothing -- while a clone lowers it
to ``+0.003478204``.  A clone is determined by its twin at every prefix and so
contributes no entropy; a genuinely new column does, and the new entropy costs
more than the extra conditioning saves.

Concentration also matters: three clones of the *same* coordinate reach
``-0.000478465`` (the registered ``n8tiny``), whereas spreading one clone over
each of three distinct coordinates only reaches ``+0.006494611`` and never turns
negative.  So the amplifier is not merely "add redundant columns"; it is "make one
coordinate maximally redundant".

Together with the saturation ceiling of Proposition 18 this closes this route
**for the record core**: its defect-free extension class is fully enumerated,
only one direction improves the objective (a clone of coordinate 2), and
repeating that direction has a finite budget.  The classification theorem is
general; the objective comparison is deliberately not extrapolated beyond this
core.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/classify_defect_free_extensions.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
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
from audit_clone_saturation import collapse  # noqa: E402

DEFAULT_REPORT = HERE / "candidates" / "defect_free_extensions.json"


# --- the classification -----------------------------------------------------


def successful_joins(rows: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    """Ordered index triples (i, j, k) with rows[i] | rows[j] == rows[k]."""
    position = {row: index for index, row in enumerate(rows)}
    out = []
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            landing = position.get(left | right)
            if landing is not None:
                out.append((i, j, landing))
    return tuple(out)


def join_consistent(rows: tuple[int, ...]) -> tuple[int, ...]:
    """Every join-consistent subset, as a bitmask over row indices.

    ``U`` qualifies iff ``[X u Y in U] = [X in U] or [Y in U]`` for every
    *successful* join.  These are exactly the defect-free coordinate extensions.
    """
    joins = successful_joins(rows)
    out = []
    for candidate in range(1 << len(rows)):
        for i, j, k in joins:
            if ((candidate >> k) & 1) != (((candidate >> i) | (candidate >> j)) & 1):
                break
        else:
            out.append(candidate)
    return tuple(out)


def is_up_set(rows: tuple[int, ...], subset: int) -> bool:
    """Is the subset upward closed inside the family under inclusion?"""
    for i, lower in enumerate(rows):
        if not (subset >> i) & 1:
            continue
        for j, upper in enumerate(rows):
            if lower | upper == upper and not (subset >> j) & 1:
                return False
    return True


def is_join_prime(rows: tuple[int, ...], subset: int) -> bool:
    """Does every successful join landing in the subset have a part in it?"""
    for i, j, k in successful_joins(rows):
        if (subset >> k) & 1 and not (((subset >> i) | (subset >> j)) & 1):
            return False
    return True


def or_columns(rows: tuple[int, ...], dimension: int) -> dict[int, tuple[int, ...]]:
    """The OR-columns: subset masks of {A : A meets S} for nonempty S."""
    out: dict[int, tuple[int, ...]] = {}
    for width in range(1, dimension + 1):
        for support in combinations(range(dimension), width):
            mask = 0
            for index, row in enumerate(rows):
                if any((row >> coordinate) & 1 for coordinate in support):
                    mask |= 1 << index
            out.setdefault(mask, support)
    return out


def extend(
    rows: tuple[int, ...], dimension: int, subset: int
) -> tuple[tuple[int, ...], int]:
    """Append the column whose 1-rows are the given subset."""
    flag = 1 << dimension
    return (
        tuple(
            sorted(
                row | flag if (subset >> index) & 1 else row
                for index, row in enumerate(rows)
            )
        ),
        dimension + 1,
    )

def realize_columns(
    row_count: int, columns: tuple[int, ...]
) -> tuple[tuple[int, ...], int]:
    """Realize several columns at once without losing the original row identity."""
    rows = tuple(
        sorted(
            sum(
                ((column >> index) & 1) << coordinate
                for coordinate, column in enumerate(columns)
            )
            for index in range(row_count)
        )
    )
    assert len(set(rows)) == row_count
    return rows, len(columns)


# --- exact objective --------------------------------------------------------


def exact_a_plus(rows: tuple[int, ...], dimension: int) -> dict[str, object]:
    representatives, group = order_orbits(rows, dimension)
    values = [bellman_bounds(rows, dimension, order) for order in representatives]
    low = round_down(sum(value[0] for value in values) / len(values), ACCUMULATOR_BITS)
    high = round_up(sum(value[1] for value in values) / len(values), ACCUMULATOR_BITS)
    iid_low, iid_high = shapley_iid_bounds(rows, dimension)
    size = Fraction(len(rows))
    upper = round_up(
        (1 - ALPHA) * iid_high + ALPHA * high - log2_lower(size), ACCUMULATOR_BITS
    )
    lower = round_down(
        (1 - ALPHA) * iid_low + ALPHA * low - log2_upper(size), ACCUMULATOR_BITS
    )
    return {
        "a_plus_lower": str(lower),
        "a_plus_upper": str(upper),
        "a_plus_upper_decimal": decimal_upper(upper, 12),
        "certified_negative": upper < 0,
        "order_orbit_count": len(representatives),
        "automorphism_count": len(group),
        "_upper": upper,
    }


# --- audits -----------------------------------------------------------------


def classify_columns(
    rows: tuple[int, ...], dimension: int
) -> dict[str, object]:
    """Every defect-free column, split into clones, other OR-columns, and new ones."""
    size = len(rows)
    limit = cap(size)
    consistent = join_consistent(rows)
    ors = or_columns(rows, dimension)
    clones = {}
    for coordinate in range(dimension):
        mask = 0
        for index, row in enumerate(rows):
            if (row >> coordinate) & 1:
                mask |= 1 << index
        clones[mask] = coordinate
    usable = [
        candidate
        for candidate in consistent
        if 0 < bin(candidate).count("1") <= limit
    ]
    rows_of = lambda mask: [  # noqa: E731
        rows[index] for index in range(size) if (mask >> index) & 1
    ]
    entries = []
    for candidate in sorted(usable, key=lambda c: (bin(c).count("1"), c)):
        if candidate in clones:
            kind = f"clone of coordinate {clones[candidate]}"
        elif candidate in ors:
            kind = f"OR-column on {list(ors[candidate])}"
        else:
            kind = "new: not induced by a global join-homomorphism"
        entries.append(
            {
                "weight": bin(candidate).count("1"),
                "rows": rows_of(candidate),
                "kind": kind,
                "is_up_set": is_up_set(rows, candidate),
                "is_join_prime": is_join_prime(rows, candidate),
            }
        )
    return {
        "join_consistent_total": len(consistent),
        "cap": limit,
        "usable_total": len(usable),
        "or_column_weights": sorted(bin(mask).count("1") for mask in ors),
        "columns": entries,
        "_usable": usable,
        "_clones": clones,
        "_ors": ors,
    }


def equivalence_audit(rows: tuple[int, ...]) -> dict[str, object]:
    """join-consistent  <=>  up-set and join-prime, over every subset."""
    consistent = set(join_consistent(rows))
    mismatches = 0
    for candidate in range(1 << len(rows)):
        expected = is_up_set(rows, candidate) and is_join_prime(rows, candidate)
        if expected != (candidate in consistent):
            mismatches += 1
    return {
        "subsets_tested": 1 << len(rows),
        "join_consistent_total": len(consistent),
        "mismatches": mismatches,
    }


def invariance_audit(
    rows: tuple[int, ...], dimension: int, subset: int
) -> dict[str, object]:
    """A defect-free extension moves neither size nor defect, and stays admissible."""
    extended, width = extend(rows, dimension, subset)
    return {
        "size_before": len(rows),
        "size_after": len(extended),
        "defect_before": str(defect(rows)),
        "defect_after": str(defect(extended)),
        "size_preserved": len(rows) == len(extended),
        "defect_preserved": defect(rows) == defect(extended),
        "degrees_after": list(degrees(extended, width)),
        "new_degree": bin(subset).count("1"),
        "admissible_before": admissible(rows, dimension),
        "admissible_after": admissible(extended, width),
        "incidence_before": incidence(rows),
        "incidence_after": incidence(extended),
        "reimer_threshold": reimer_threshold(len(rows)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=str, default="n8tiny")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    arguments = parser.parse_args()

    base = BASES[arguments.base]
    registered = tuple(sorted(base.reconstruct()))
    core, kept = collapse(registered, base.dimension)
    width = len(kept)
    print(f"base {arguments.base}: n={base.dimension} m={len(registered)}")
    print(f"collapsed core: n={width} m={len(core)} defect={defect(core)} rows={list(core)}")

    equivalence = equivalence_audit(core)
    print(
        f"\nclassification audit: {equivalence['subsets_tested']} subsets, "
        f"{equivalence['join_consistent_total']} join-consistent, "
        f"{equivalence['mismatches']} mismatches against 'up-set and join-prime'"
    )

    classified = classify_columns(core, width)
    print(
        f"\ndefect-free columns respecting the cap {classified['cap']}: "
        f"{classified['usable_total']} of {classified['join_consistent_total']} "
        f"join-consistent subsets"
    )
    for entry in classified["columns"]:
        print(f"  weight {entry['weight']}  rows {entry['rows']}  -- {entry['kind']}")

    baseline = exact_a_plus(core, width)
    print(f"\ncore objective: A_+ <= {baseline['a_plus_upper_decimal']}")
    measured = []
    for candidate in classified["_usable"]:
        extended, extended_width = extend(core, width, candidate)
        facts = exact_a_plus(extended, extended_width)
        invariance = invariance_audit(core, width, candidate)
        kind = next(
            entry["kind"]
            for entry in classified["columns"]
            if entry["rows"] == [core[i] for i in range(len(core)) if (candidate >> i) & 1]
        )
        delta = facts["_upper"] - baseline["_upper"]
        measured.append(
            {
                "kind": kind,
                "weight": bin(candidate).count("1"),
                "a_plus_upper_decimal": facts["a_plus_upper_decimal"],
                "improves": delta < 0,
                "invariance": invariance,
            }
        )
        print(
            f"  + {kind:44s} -> A_+ <= {facts['a_plus_upper_decimal']:>15s}  "
            f"{'improves' if delta < 0 else 'WORSE'}  "
            f"defect preserved: {invariance['defect_preserved']}"
        )

    improving = [entry for entry in measured if entry["improves"]]
    only_clones = all(entry["kind"].startswith("clone") for entry in improving)
    all_preserved = all(
        entry["invariance"]["defect_preserved"] and entry["invariance"]["size_preserved"]
        and entry["invariance"]["admissible_after"]
        for entry in measured
    )
    print(
        f"\nimproving columns: {len(improving)} of {len(measured)}; "
        f"all of them clones: {only_clones}"
    )

    # Concentration: realize every appended column against the original row
    # indices.  Sequential ``extend`` calls would sort rows after each append and
    # invalidate the next index-bitmask.
    concentration = {}
    coordinate_columns = {
        coordinate: candidate for candidate, coordinate in classified["_clones"].items()
    }
    original_columns = tuple(coordinate_columns[i] for i in range(width))
    for label, multiset in (
        ("three clones of coordinate 2", (2, 2, 2)),
        ("one clone each of coordinates 0,1,2", (0, 1, 2)),
    ):
        columns = original_columns + tuple(
            coordinate_columns[coordinate] for coordinate in multiset
        )
        rows, dimension = realize_columns(len(core), columns)
        facts = exact_a_plus(rows, dimension)
        concentration[label] = {
            "dimension": dimension,
            "defect": str(defect(rows)),
            "a_plus_upper": facts["a_plus_upper"],
            "a_plus_upper_decimal": facts["a_plus_upper_decimal"],
            "certified_negative": facts["certified_negative"],
        }
        print(
            f"  {label:38s} n={dimension} defect={defect(rows)} "
            f"A_+ <= {facts['a_plus_upper_decimal']:>15s} "
            f"negative: {facts['certified_negative']}"
        )

    verdict = (
        "DEFECT_FREE_EXTENSIONS_CLASSIFIED_ONLY_CLONES_IMPROVE"
        if equivalence["mismatches"] == 0 and all_preserved and only_clones
        else "INCONCLUSIVE"
    )
    report = {
        "base": arguments.base,
        "core_rows": list(core),
        "core_dimension": width,
        "core_defect": str(defect(core)),
        "core_a_plus_upper_decimal": baseline["a_plus_upper_decimal"],
        "equivalence_audit": equivalence,
        "classification": {
            key: value
            for key, value in classified.items()
            if not key.startswith("_")
        },
        "measured": measured,
        "improving_are_only_clones": only_clones,
        "every_extension_preserved_size_defect_admissibility": all_preserved,
        "concentration": concentration,
        "verdict": verdict,
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"\n{verdict}")
    print(f"report: {arguments.report}")


if __name__ == "__main__":
    main()
