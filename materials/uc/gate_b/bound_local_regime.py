#!/usr/bin/env python3
"""Exact structural bounds for the Gate B local regime.

Gate B proved ``c_cl^star = sup -A_+/eps_vee = +infinity`` using Cartesian powers
whose closure defect tends to *one*.  The surviving question is the local one:
how does the repair ratio behave when the defect is *small*?  This module proves
the part of that question that is decidable by exact counting, and it isolates
precisely what is left.

Everything here is integer or ``Fraction`` arithmetic.  No floating point, no
interval library, no entropy evaluation, and no import from the search or
objective modules: the four statements below are re-derived from the definitions
in ``DEFINITIONS.md`` so that this file is an independent checker.

Statements
----------
``L7`` size ceiling.  Incidence is at most ``n * floor(2m/5) <= 2nm/5`` and
Reimer demands at least ``m log2 m / 2``, so every admissible family satisfies

    log2 m <= 4n/5,     equivalently the exact integer test  m^5 <= 2^(4n).

``C8`` numerator ceiling.  ``Q >= 0`` and ``C_+ >= 0`` give
``-A_+ = log2 m - (1-alpha)Q - alpha C_+ <= log2 m``, hence with ``L7``

    -A_+ <= 4n/5   for every admissible family on n coordinates.

The certified powers give ``-A_+(B^k) = k(-A_+(B))`` on ``7k`` coordinates, so
the numerator's growth is ``Theta(n)`` -- bounded on both sides.

``L9`` union growth.  Writing ``p_i`` for the coordinate frequencies, the cap
gives ``p_i <= 2/5`` and therefore

    E|X or Y| = sum_i p_i(2 - p_i) >= (8/5) * sbar,      sbar = (1/m) sum |A|.

A union that lands in the family has size at most ``M = max |A|``; otherwise it
has size at most ``n``.  So ``E|X or Y| <= M + eps_vee (n - M)`` and

    eps_vee >= ((8/5) sbar - M) / (n - M)        when M < n.

``L11`` downset capacity.  ``X or Y = A`` forces ``X, Y`` to be members inside
``A``, so at most ``N(A)^2`` ordered pairs have union ``A``, where
``N(A) = #{B in F : B subseteq A}``.  Summing over members,

    eps_vee >= 1 - (1/m^2) sum_{A in F} N(A)^2.

``P12`` maximality.  At ``n = 6`` and ``n = 7`` the largest admissible sizes are
25 and 45.  Both are attained by the certified bases with *every* coordinate
degree exactly at the cap, so no set at all can be added to them: a nonempty
addition breaks the cap and the empty set breaks Reimer.  The "add the missing
unions" amplification route is therefore blocked at every certified base, not
merely unpromising.

Consequence
-----------
``L9`` plus Reimer bounds the ratio on every family without a dominant set: if
``M <= theta * (8/5) * sbar`` with ``theta < 1`` then

    -A_+/eps_vee <= (4n/5)(n - M) / ((8/5) sbar - M) = O(n),

matching the certified ``Theta(n)`` lower bound.  Divergence faster than
``Theta(n)``, and any approach to the local regime, therefore *requires* a set of
size at least ``(8/5) sbar >= (4/5) log2 m``.  Both `n=7` bases have such a set
(sizes 6 and 7); this is the structural feature that makes them possible.

What is *not* decided.  A family with ``eps_vee = 0`` is union closed with every
frequency at most ``2/5``, so ``c_cl^star`` restricted to at most ``n``
coordinates is finite exactly when no such family has ``A_+ < 0``.  Frankl is
verified for ``n <= 12``, which settles those dimensions; a proof for all ``n``
would give a ``2/5`` frequency bound, beyond the published ``0.381966``.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/bound_local_regime.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys
from typing import Iterable

HERE = Path(__file__).resolve().parent
DEFAULT_REPORT = HERE / "candidates" / "local_regime_bounds.json"


# --- exact primitives, re-derived from DEFINITIONS.md ------------------------


def reimer_threshold(size: int) -> int:
    """Least integer r with 2^(2r) >= size^size, i.e. ceil(m log2 m / 2)."""
    if size <= 1:
        return 0
    target = size**size
    low, high = 0, size * max(1, (size - 1).bit_length())
    while low < high:
        middle = (low + high) // 2
        if 1 << (2 * middle) >= target:
            high = middle
        else:
            low = middle + 1
    return low


def cap(size: int) -> int:
    return 2 * size // 5


def degrees(rows: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    return tuple(
        sum((row >> coordinate) & 1 for row in rows)
        for coordinate in range(dimension)
    )


def incidence(rows: Iterable[int]) -> int:
    return sum(bin(row).count("1") for row in rows)


def admissible(rows: tuple[int, ...], dimension: int) -> bool:
    size = len(rows)
    if size < 3 or len(set(rows)) != size:
        return False
    if max(rows, default=0) >= 1 << dimension:
        return False
    if max(degrees(rows, dimension), default=0) > cap(size):
        return False
    return incidence(rows) >= reimer_threshold(size)


def defect_count(rows: tuple[int, ...]) -> int:
    members = set(rows)
    return sum(1 for a in rows for b in rows if (a | b) not in members)


def defect(rows: tuple[int, ...]) -> Fraction:
    return Fraction(defect_count(rows), len(rows) ** 2)


# --- L7: size ceiling -------------------------------------------------------


def size_ceiling_holds(size: int, dimension: int) -> bool:
    """Exact test of log2 m <= 4n/5, as the integer inequality m^5 <= 2^(4n)."""
    return size**5 <= 1 << (4 * dimension)


def feasible_sizes(dimension: int) -> tuple[int, ...]:
    """Sizes admitting the cap and Reimer constraints simultaneously."""
    return tuple(
        size
        for size in range(3, (1 << dimension) + 1)
        if reimer_threshold(size) <= dimension * cap(size)
    )


def max_admissible_size(dimension: int) -> int:
    return max(feasible_sizes(dimension))


# --- L9: union growth floor -------------------------------------------------


def union_growth_floor(rows: tuple[int, ...], dimension: int) -> Fraction | None:
    """Lower bound on eps_vee from the cap and membership.  None if vacuous."""
    size = len(rows)
    largest = max(bin(row).count("1") for row in rows)
    if largest >= dimension:
        return None
    mean = Fraction(incidence(rows), size)
    value = (Fraction(8, 5) * mean - largest) / (dimension - largest)
    return value if value > 0 else None


def reimer_union_growth_floor(size: int, dimension: int, largest: int) -> Fraction | None:
    """L9 with sbar replaced by its Reimer lower bound: depends only on shape."""
    if largest >= dimension:
        return None
    mean = Fraction(reimer_threshold(size), size)
    value = (Fraction(8, 5) * mean - largest) / (dimension - largest)
    return value if value > 0 else None


def ratio_ceiling(rows: tuple[int, ...], dimension: int) -> Fraction | None:
    """Fully rational upper bound on -A_+/eps_vee via C8 and L9."""
    floor_value = union_growth_floor(rows, dimension)
    if floor_value is None:
        return None
    numerator_ceiling = Fraction(4 * dimension, 5)
    return numerator_ceiling / floor_value


# --- L11: downset capacity --------------------------------------------------


def downset_counts(rows: tuple[int, ...]) -> dict[int, int]:
    return {a: sum(1 for b in rows if b | a == a) for a in rows}


def capacity_floor(rows: tuple[int, ...]) -> Fraction:
    size = len(rows)
    counts = downset_counts(rows)
    return 1 - Fraction(sum(v * v for v in counts.values()), size * size)


# --- P12: maximality --------------------------------------------------------


def admissible_additions(rows: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    """Every set whose addition keeps the family admissible."""
    present = set(rows)
    return tuple(
        candidate
        for candidate in range(1 << dimension)
        if candidate not in present
        and admissible(tuple(sorted(rows + (candidate,))), dimension)
    )


# --- audits -----------------------------------------------------------------


def audit_exhaustive(dimension: int) -> dict[str, object]:
    """Check every bound against every admissible family of this dimension."""
    tested = 0
    violations: list[str] = []
    growth_active = capacity_active = 0
    tightest_growth = tightest_capacity = None
    for size in feasible_sizes(dimension):
        for rows in combinations(range(1 << dimension), size):
            if not admissible(rows, dimension):
                continue
            tested += 1
            observed = defect(rows)
            if not size_ceiling_holds(size, dimension):
                violations.append(f"L7 size ceiling failed at m={size}")
            growth = union_growth_floor(rows, dimension)
            if growth is not None:
                growth_active += 1
                if observed < growth:
                    violations.append(f"L9 floor {growth} > defect {observed}")
                slack = observed - growth
                if tightest_growth is None or slack < tightest_growth:
                    tightest_growth = slack
            capacity = capacity_floor(rows)
            if capacity > 0:
                capacity_active += 1
                if observed < capacity:
                    violations.append(f"L11 floor {capacity} > defect {observed}")
                slack = observed - capacity
                if tightest_capacity is None or slack < tightest_capacity:
                    tightest_capacity = slack
    return {
        "dimension": dimension,
        "admissible_families": tested,
        "union_growth_active": growth_active,
        "capacity_active": capacity_active,
        "violations": violations,
        "tightest_union_growth_slack": None if tightest_growth is None else str(tightest_growth),
        "tightest_capacity_slack": None if tightest_capacity is None else str(tightest_capacity),
    }


def audit_shape_floor(dimensions: Iterable[int]) -> dict[str, object]:
    """Where the Reimer-only floor is active: it depends on (n, m, M) alone."""
    active = {}
    for dimension in dimensions:
        rows = []
        for size in feasible_sizes(dimension):
            best = None
            for largest in range(dimension):
                value = reimer_union_growth_floor(size, dimension, largest)
                if value is None:
                    continue
                if best is None or value > best[1]:
                    best = (largest, value)
            if best is not None:
                rows.append({"size": size, "largest_set": best[0], "floor": str(best[1])})
        threshold = {}
        for size in feasible_sizes(dimension):
            mean = Fraction(reimer_threshold(size), size)
            threshold[str(size)] = str(Fraction(8, 5) * mean)
        active[str(dimension)] = {
            "max_admissible_size": max_admissible_size(dimension),
            "size_ceiling_holds": size_ceiling_holds(max_admissible_size(dimension), dimension),
            "dominant_set_threshold": threshold[str(max_admissible_size(dimension))],
            "floor_rows": rows[-3:],
        }
    return active


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-dimension", type=int, default=4)
    parser.add_argument("--bases", type=str, default="")
    args = parser.parse_args()

    summary: dict[str, object] = {"statements": ["L7", "C8", "L9", "L11", "P12"]}

    ceiling = {}
    for dimension in range(3, 17):
        largest = max_admissible_size(dimension)
        ceiling[str(dimension)] = {
            "max_admissible_size": largest,
            "ceiling_holds": size_ceiling_holds(largest, dimension),
        }
    summary["size_ceiling"] = ceiling
    summary["size_ceiling_all_hold"] = all(v["ceiling_holds"] for v in ceiling.values())

    summary["exhaustive_audits"] = [
        audit_exhaustive(dimension) for dimension in range(3, args.max_dimension + 1)
    ]
    summary["shape_floor"] = audit_shape_floor(range(5, 9))

    if args.bases:
        sys.path.insert(0, str(HERE))
        import verify_gate_b_rational as rational

        rows_out = []
        for name in args.bases.split(","):
            base = rational.BASES[name.strip()]
            rows = tuple(sorted(base.reconstruct()))
            dimension = base.dimension
            additions = admissible_additions(rows, dimension)
            growth = union_growth_floor(rows, dimension)
            rows_out.append(
                {
                    "name": name.strip(),
                    "dimension": dimension,
                    "size": len(rows),
                    "is_max_admissible_size": len(rows) == max_admissible_size(dimension),
                    "degrees_all_at_cap": set(degrees(rows, dimension)) == {cap(len(rows))},
                    "largest_set": max(bin(row).count("1") for row in rows),
                    "dominant_set_threshold": str(Fraction(8, 5) * Fraction(incidence(rows), len(rows))),
                    "admissible_additions": len(additions),
                    "defect": str(defect(rows)),
                    "union_growth_floor": None if growth is None else str(growth),
                    "capacity_floor": str(capacity_floor(rows)),
                    "ratio_ceiling": (
                        None if ratio_ceiling(rows, dimension) is None
                        else str(ratio_ceiling(rows, dimension))
                    ),
                }
            )
        summary["certified_bases"] = rows_out
        summary["no_base_admits_an_addition"] = all(
            row["admissible_additions"] == 0 for row in rows_out
        )

    violations = [
        problem
        for audit in summary["exhaustive_audits"]
        for problem in audit["violations"]
    ]
    summary["verdict"] = (
        "ALL_BOUNDS_HOLD" if not violations and summary["size_ceiling_all_hold"]
        else "VIOLATION"
    )

    encoded = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(encoded)
    print(json.dumps({k: summary[k] for k in ("verdict", "size_ceiling_all_hold")}))
    for audit in summary["exhaustive_audits"]:
        print(
            f"n={audit['dimension']}: admissible={audit['admissible_families']} "
            f"L9_active={audit['union_growth_active']} "
            f"L11_active={audit['capacity_active']} "
            f"violations={len(audit['violations'])} "
            f"tightest_L9={audit['tightest_union_growth_slack']} "
            f"tightest_L11={audit['tightest_capacity_slack']}"
        )
    if "certified_bases" in summary:
        for row in summary["certified_bases"]:
            print(
                f"{row['name']}: m={row['size']} max_size={row['is_max_admissible_size']} "
                f"degrees_at_cap={row['degrees_all_at_cap']} additions={row['admissible_additions']} "
                f"largest_set={row['largest_set']} threshold={row['dominant_set_threshold']} "
                f"L9={row['union_growth_floor']} L11={row['capacity_floor']}"
            )


if __name__ == "__main__":
    main()
