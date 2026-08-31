#!/usr/bin/env python3
"""Exact cloning closure: clones realize precisely the fixed-order extrema.

Proposition
-----------
Fix a family ``F`` on distinct coordinates and replace coordinate ``i`` by
``r_i >= 1`` identical copies.  In a coordinate order of the cloned family,
retain the *first* copy of each clone class and delete all later copies.  Once
the first copy has been exposed, every later copy is deterministic in both
rows, contributes zero entropy, and leaves every row fiber unchanged.  Thus,
fixed order by fixed order,

    Q_pi(cloned F)   = Q_collapse(pi)(F),
    C+_pi(cloned F) = C+_collapse(pi)(F),
    A+_pi(cloned F) = A+_collapse(pi)(F).

For multiplicities ``r = (r_0,...,r_(n-1))``, the collapsed original order has
the exact Plackett--Luce law

    P_r(pi) = product_k r[pi[k]] / sum_{j >= k} r[pi[j]].

Therefore every finite cloned representation is a convex combination of the
original fixed-order values.  The converse is sharp at both extrema.  For any
target order ``pi`` and integer ``R >= 2``, set
``r[pi[k]] = R**(n-k-1)``.  Then

    P_r(pi) >= (1 - 1/R)**(n-1) >= 1 - (n-1)/R,

so the distribution converges to the point mass at ``pi``.  Consequently

    inf_r A+(cloned_r F) = min_pi A+_pi(F),
    sup_r A+(cloned_r F) = max_pi A+_pi(F).

The extrema are infimum/supremum statements; a nonconstant profile need not
attain them at finite multiplicity.  Sign reachability *is* finite: a negative
fixed-order value yields an explicit finite multiplicity vector with negative
average, while nonnegative fixed-order values exclude every clone multiset.

One-coordinate formula
----------------------
If coordinate ``i`` alone has total multiplicity ``r`` (``r=1`` means no added
clone), let ``M_j`` be the average fixed-order objective over original orders
in which ``i`` has rank ``j``.  The rank of its first copy has exact law

    P_r(j) = C(n+r-j-2, r-1) / C(n+r-1, r),    j=0,...,n-1.

Hence ``A_r = sum_j P_r(j) M_j`` and ``A_r -> M_0``.  This module encloses every
``M_j``, displayed ``A_r``, and finite sharp-reachability witness with exact
rational two-sided arithmetic.

Applications
------------
* ``n8tiny`` core (defect ``4/9``): negative fixed orders exist; cloning
  coordinate 2 crosses zero at total multiplicity 4.  The general law also
  constructs a finite multicoordinate witness from its best fixed order.
* ``n8m15_clone_reachable`` core (defect ``14/45 < 2/5``): its uniform
  objective is positive but one fixed order is exactly negative.  Ratio-41
  geometric multiplicities give a finite exact negative clone certificate at
  unchanged defect.
* ``n9m15`` low-defect core (defect ``14/45 < 2/5``): every one of its 5040
  fixed-order objectives is positive.  Its one genuinely new usable
  join-consistent extension likewise has all 40320 fixed orders positive.
  Join-consistency is invariant under a defect-free extension, so these two
  profiles cover every finite sequence of defect-free extensions.

Run:
    OMP_NUM_THREADS=1 nice -n 19 python3 -B uc/gate_b/audit_clone_limit.py
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations
from math import comb
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
    round_down,
    round_up,
)
from bound_local_regime import admissible, cap, defect, degrees, incidence, reimer_threshold  # noqa: E402
from audit_clone_saturation import collapse  # noqa: E402
from barrier_sawin_cambie import coordinate_term  # noqa: E402
from classify_defect_free_extensions import classify_columns, extend  # noqa: E402

DEFAULT_REPORT = HERE / "candidates" / "clone_limit_audit.json"
DEFAULT_CERTIFICATE = (
    HERE / "certificates" / "gate_b_subtwofifths_clone_rational_v1.json"
)

# Active, separating core produced by the deterministic n=9, m=15, restart-indexed
# search (seed 20260827).  The original nine-coordinate search row had one inert
# zero column and one duplicate; deleting both changes neither the family law nor
# the defect and leaves an admissible seven-coordinate representation.
N9_M15_ROWS = (0, 2, 4, 6, 8, 9, 16, 22, 31, 32, 96, 105, 112, 125, 127)

# Normalized eight-coordinate core of the only n=9, m=20 finalist that passed
# the Q-only necessary condition at defect 19/50.  The searched representation
# had one duplicate column; removing it makes the full-order Q-only lower bound
# positive again, while retaining the same defect and admissibility.
N9_M20_QGATE_ROWS = (
    0, 1, 2, 7, 12, 16, 20, 32, 43, 48,
    56, 63, 64, 192, 193, 199, 207, 250, 252, 255,
)

# Normalized core of the lowest-defect negative-fixed-order finalist from the
# strict n=9, m=15 row-changing search.  Its searched representation has one
# duplicate-column pair; first-occurrence collapse leaves this eight-coordinate
# family and the exact negative order (6,3,2,0,7,5,4,1).
N8_M15_CLONE_REACHABLE_ROWS = (
    0, 1, 2, 4, 5, 8, 10, 43, 64, 190, 192, 193, 245, 254, 255,
)

def n9_nonclone_extension() -> tuple[tuple[int, ...], int, int, dict[str, object]]:
    classification = classify_columns(N9_M15_ROWS, 7)
    new_columns = [
        candidate
        for candidate in classification["_usable"]
        if candidate not in classification["_clones"]
        and candidate not in classification["_ors"]
    ]
    assert len(new_columns) == 1
    candidate = new_columns[0]
    extended, dimension = extend(N9_M15_ROWS, 7, candidate)
    assert defect(extended) == defect(N9_M15_ROWS)
    return extended, dimension, candidate, classification


def candidates() -> dict[str, tuple[tuple[int, ...], int, str]]:
    tiny = tuple(sorted(BASES["n8tiny"].reconstruct()))
    tiny_core, kept = collapse(tiny, BASES["n8tiny"].dimension)
    assert len(kept) == 5
    extended, extended_dimension, _column, _classification = n9_nonclone_extension()
    return {
        "n8tiny_core": (
            tiny_core,
            5,
            "collapsed from registered n8tiny; Proposition 18 baseline",
        ),
        "n9m15_lowdefect": (
            N9_M15_ROWS,
            7,
            "search_local_defect schema 3, n=9 m=15 seed 20260827; inactive and duplicate columns removed",
        ),
        "n8m15_clone_reachable": (
            N8_M15_CLONE_REACHABLE_ROWS,
            8,
            "subtwofifths_min_a_n9_under_16_45 restart 11; duplicate column removed",
        ),
        "n9m15_plus_new": (
            extended,
            extended_dimension,
            "the unique usable non-clone join-consistent extension of n9m15_lowdefect",
        ),
        "n9m20_qgate": (
            N9_M20_QGATE_ROWS,
            8,
            "subtwofifths_q_n9 m=20 survivor; its one duplicate column removed",
        ),
    }


def normalized(rows: tuple[int, ...], dimension: int) -> bool:
    columns = tuple(
        tuple((row >> coordinate) & 1 for row in rows)
        for coordinate in range(dimension)
    )
    return all(any(column) for column in columns) and len(set(columns)) == dimension


def fixed_order_q_bounds(
    rows: tuple[int, ...], order: tuple[int, ...]
) -> tuple[Fraction, Fraction]:
    low = high = Fraction(0)
    predecessors: list[int] = []
    for coordinate in order:
        term_low, term_high = coordinate_term(
            rows, coordinate, tuple(predecessors)
        )
        low = round_down(low + term_low, ACCUMULATOR_BITS)
        high = round_up(high + term_high, ACCUMULATOR_BITS)
        predecessors.append(coordinate)
    return low, high


def fixed_order_a_bounds(
    rows: tuple[int, ...], dimension: int, order: tuple[int, ...]
) -> tuple[Fraction, Fraction]:
    q_low, q_high = fixed_order_q_bounds(rows, order)
    c_low, c_high, _states = bellman_bounds(rows, dimension, order)
    size = Fraction(len(rows))
    low = round_down(
        (1 - ALPHA) * q_low + ALPHA * c_low - log2_upper(size),
        ACCUMULATOR_BITS,
    )
    high = round_up(
        (1 - ALPHA) * q_high + ALPHA * c_high - log2_lower(size),
        ACCUMULATOR_BITS,
    )
    return low, high


def multiplicity_weights(
    dimension: int, multiplicity: int
) -> tuple[Fraction, ...]:
    """Rank law of the first clone among ``multiplicity`` identical copies."""
    assert multiplicity >= 1
    denominator = comb(dimension + multiplicity - 1, multiplicity)
    weights = tuple(
        Fraction(
            comb(dimension + multiplicity - rank - 2, multiplicity - 1),
            denominator,
        )
        for rank in range(dimension)
    )
    assert sum(weights) == 1
    return weights


def combine_rank_profile(
    profile: list[dict[str, object]], dimension: int, multiplicity: int
) -> tuple[Fraction, Fraction]:
    weights = multiplicity_weights(dimension, multiplicity)
    low = round_down(
        sum(weight * entry["a_plus_lower"] for weight, entry in zip(weights, profile)),
        ACCUMULATOR_BITS,
    )
    high = round_up(
        sum(weight * entry["a_plus_upper"] for weight, entry in zip(weights, profile)),
        ACCUMULATOR_BITS,
    )
    return low, high

def collapsed_order_probability(
    order: tuple[int, ...], multiplicities: tuple[int, ...]
) -> Fraction:
    """Exact first-appearance order law for uniformly ordered clone copies."""
    dimension = len(multiplicities)
    assert dimension >= 1
    assert tuple(sorted(order)) == tuple(range(dimension))
    assert all(isinstance(value, int) and value >= 1 for value in multiplicities)
    remaining = sum(multiplicities)
    probability = Fraction(1)
    for coordinate in order:
        weight = multiplicities[coordinate]
        probability *= Fraction(weight, remaining)
        remaining -= weight
    assert remaining == 0
    return probability


def collapsed_order_distribution(
    multiplicities: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], Fraction], ...]:
    """All collapsed orders and their exact Plackett--Luce probabilities."""
    distribution = tuple(
        (order, collapsed_order_probability(order, multiplicities))
        for order in permutations(range(len(multiplicities)))
    )
    assert sum(probability for _order, probability in distribution) == 1
    return distribution

def combine_order_bounds(
    order_bounds: tuple[
        tuple[tuple[int, ...], Fraction, Fraction], ...
    ],
    multiplicities: tuple[int, ...],
) -> tuple[Fraction, Fraction]:
    """Exact enclosure of a cloned average from every original order."""
    total_probability = Fraction(0)
    lower = upper = Fraction(0)
    for order, order_lower, order_upper in order_bounds:
        probability = collapsed_order_probability(order, multiplicities)
        total_probability += probability
        lower += probability * order_lower
        upper += probability * order_upper
    assert total_probability == 1
    return (
        round_down(lower, ACCUMULATOR_BITS),
        round_up(upper, ACCUMULATOR_BITS),
    )


def geometric_multiplicities(
    order: tuple[int, ...], ratio: int
) -> tuple[int, ...]:
    """Integer clone weights that concentrate first appearances on ``order``."""
    dimension = len(order)
    assert dimension >= 1
    assert tuple(sorted(order)) == tuple(range(dimension))
    assert ratio >= 2
    multiplicities = [0] * dimension
    for rank, coordinate in enumerate(order):
        multiplicities[coordinate] = ratio ** (dimension - rank - 1)
    return tuple(multiplicities)


def finite_negative_clone_witness(
    target_order: tuple[int, ...],
    target_lower: Fraction,
    target_upper: Fraction,
    minimum_lower: Fraction,
    maximum_upper: Fraction,
    order_bounds: tuple[
        tuple[tuple[int, ...], Fraction, Fraction], ...
    ] | None = None,
) -> dict[str, object] | None:
    """Construct an exact two-sided enclosure for a finite negative clone vector."""
    if target_upper >= 0:
        return None
    dimension = len(target_order)
    assert dimension >= 1
    assert minimum_lower <= target_lower <= target_upper <= maximum_upper
    if dimension == 1:
        return {
            "target_order": list(target_order),
            "geometric_ratio": 1,
            "multiplicities": [1],
            "cloned_dimension": 1,
            "target_order_probability": "1",
            "target_order_probability_lower": "1",
            "union_bound_probability_lower": "1",
            "a_plus_lower": str(target_lower),
            "a_plus_upper": str(target_upper),
            "a_plus_lower_decimal": decimal_lower(target_lower, 12),
            "a_plus_upper_decimal": decimal_upper(target_upper, 12),
            "certified_negative": True,
        }

    spread = maximum_upper - target_upper
    sufficient = Fraction(dimension - 1) * spread / (-target_upper)
    ratio = max(2, sufficient.numerator // sufficient.denominator + 1)
    while True:
        multiplicities = geometric_multiplicities(target_order, ratio)
        probability = collapsed_order_probability(target_order, multiplicities)
        probability_lower = Fraction(ratio - 1, ratio) ** (dimension - 1)
        union_lower = Fraction(1) - Fraction(dimension - 1, ratio)
        assert probability >= probability_lower >= union_lower
        lower = round_down(
            probability * target_lower
            + (1 - probability) * minimum_lower,
            ACCUMULATOR_BITS,
        )
        upper = round_up(
            probability * target_upper
            + (1 - probability) * maximum_upper,
            ACCUMULATOR_BITS,
        )
        if upper < 0:
            break
        ratio *= 2
    bound_method = "two-sided convex envelope over exact fixed-order bounds"
    if order_bounds is not None:
        exact_ratios = list(range(2, min(ratio, 64) + 1))
        probe = 128
        while probe < ratio:
            exact_ratios.append(probe)
            probe *= 2
        if not exact_ratios or exact_ratios[-1] != ratio:
            exact_ratios.append(ratio)
        for exact_ratio in exact_ratios:
            exact_multiplicities = geometric_multiplicities(
                target_order, exact_ratio
            )
            exact_lower, exact_upper = combine_order_bounds(
                order_bounds, exact_multiplicities
            )
            if exact_upper < 0:
                ratio = exact_ratio
                multiplicities = exact_multiplicities
                probability = collapsed_order_probability(
                    target_order, multiplicities
                )
                probability_lower = Fraction(ratio - 1, ratio) ** (
                    dimension - 1
                )
                union_lower = Fraction(1) - Fraction(dimension - 1, ratio)
                lower, upper = exact_lower, exact_upper
                bound_method = "exact reweighting of every fixed-order enclosure"
                break
        else:
            raise AssertionError("the convex certificate promised a finite witness")
    assert lower <= upper
    return {
        "target_order": list(target_order),
        "geometric_ratio": ratio,
        "multiplicities": list(multiplicities),
        "cloned_dimension": sum(multiplicities),
        "target_order_probability": str(probability),
        "target_order_probability_lower": str(probability_lower),
        "union_bound_probability_lower": str(union_lower),
        "a_plus_lower": str(lower),
        "a_plus_upper": str(upper),
        "a_plus_lower_decimal": decimal_lower(lower, 12),
        "a_plus_upper_decimal": decimal_upper(upper, 12),
        "bound_method": bound_method,
        "certified_negative": True,
    }


def audit_family(
    name: str, rows: tuple[int, ...], dimension: int, source: str, max_multiplicity: int
) -> dict[str, object]:
    count = 1
    sums = [
        [[Fraction(0), Fraction(0)] for _rank in range(dimension)]
        for _coordinate in range(dimension)
    ]
    fixed_min_low: Fraction | None = None
    fixed_min_high: Fraction | None = None
    fixed_max_low: Fraction | None = None
    fixed_max_high: Fraction | None = None
    min_order: tuple[int, ...] | None = None
    max_order: tuple[int, ...] | None = None
    total_low = total_high = Fraction(0)
    order_bounds: list[tuple[tuple[int, ...], Fraction, Fraction]] = []

    for order in permutations(range(dimension)):
        low, high = fixed_order_a_bounds(rows, dimension, order)
        total_low += low
        total_high += high
        order_bounds.append((order, low, high))
        if fixed_min_low is None or low < fixed_min_low:
            fixed_min_low = low
        if fixed_min_high is None or high < fixed_min_high:
            fixed_min_high = high
            min_order = order
        if fixed_max_low is None or low > fixed_max_low:
            fixed_max_low = low
            max_order = order
        if fixed_max_high is None or high > fixed_max_high:
            fixed_max_high = high
        for rank, coordinate in enumerate(order):
            sums[coordinate][rank][0] += low
            sums[coordinate][rank][1] += high
        count += 1
    order_count = count - 1
    per_rank_count = order_count // dimension
    assert per_rank_count * dimension == order_count
    assert fixed_min_low is not None
    assert fixed_min_high is not None
    assert fixed_max_low is not None
    assert fixed_max_high is not None
    assert min_order is not None
    assert max_order is not None

    profiles: dict[str, object] = {}
    best_coordinate = None
    best_limit_upper = None
    global_rank_min_lower = None
    for coordinate in range(dimension):
        profile: list[dict[str, object]] = []
        for rank in range(dimension):
            low = round_down(
                sums[coordinate][rank][0] / per_rank_count, ACCUMULATOR_BITS
            )
            high = round_up(
                sums[coordinate][rank][1] / per_rank_count, ACCUMULATOR_BITS
            )
            profile.append({
                "rank": rank,
                "a_plus_lower": low,
                "a_plus_upper": high,
                "a_plus_lower_decimal": decimal_lower(low, 12),
                "a_plus_upper_decimal": decimal_upper(high, 12),
            })
            if global_rank_min_lower is None or low < global_rank_min_lower:
                global_rank_min_lower = low
        limit_upper = profile[0]["a_plus_upper"]
        if best_limit_upper is None or limit_upper < best_limit_upper:
            best_limit_upper = limit_upper
            best_coordinate = coordinate
        finite = []
        for multiplicity in range(1, max_multiplicity + 1):
            low, high = combine_rank_profile(profile, dimension, multiplicity)
            finite.append({
                "total_multiplicity": multiplicity,
                "added_clones": multiplicity - 1,
                "a_plus_lower": str(low),
                "a_plus_upper": str(high),
                "a_plus_lower_decimal": decimal_lower(low, 12),
                "a_plus_upper_decimal": decimal_upper(high, 12),
                "certified_negative": high < 0,
                "certified_positive": low > 0,
            })
        profiles[str(coordinate)] = {
            "rank_profile": [
                {
                    key: str(value) if key in ("a_plus_lower", "a_plus_upper") else value
                    for key, value in entry.items()
                }
                for entry in profile
            ],
            "limit_lower": str(profile[0]["a_plus_lower"]),
            "limit_upper": str(profile[0]["a_plus_upper"]),
            "limit_lower_decimal": profile[0]["a_plus_lower_decimal"],
            "limit_upper_decimal": profile[0]["a_plus_upper_decimal"],
            "multiplicity_ladder": finite,
        }

    average_low = round_down(total_low / order_count, ACCUMULATOR_BITS)
    average_high = round_up(total_high / order_count, ACCUMULATOR_BITS)
    every_fixed_order_positive = fixed_min_low > 0
    negative_fixed_order = fixed_min_high < 0
    sharp_witness = finite_negative_clone_witness(
        min_order,
        fixed_min_low,
        fixed_min_high,
        fixed_min_low,
        fixed_max_high,
        tuple(order_bounds),
    )
    print(
        f"{name}: n={dimension} m={len(rows)} defect={defect(rows)} "
        f"A_+ in [{decimal_lower(average_low, 12)}, {decimal_upper(average_high, 12)}]"
    )
    print(
        f"  fixed-order min >= {decimal_lower(fixed_min_low, 12)}; "
        f"every order positive: {every_fixed_order_positive}"
    )
    if sharp_witness is not None:
        print(
            "  finite multicoordinate clone witness: "
            f"R={sharp_witness['geometric_ratio']}, "
            f"A_+ <= {sharp_witness['a_plus_upper_decimal']}"
        )
    print(
        f"  best one-coordinate clone limit: coordinate {best_coordinate}, "
        f"A_+ <= {decimal_upper(best_limit_upper, 12)}"
    )
    return {
        "name": name,
        "source": source,
        "dimension": dimension,
        "family_size": len(rows),
        "family_rows": list(rows),
        "normalized": normalized(rows, dimension),
        "admissible": admissible(rows, dimension),
        "coordinate_counts": list(degrees(rows, dimension)),
        "cap": cap(len(rows)),
        "incidence": incidence(rows),
        "reimer_threshold": reimer_threshold(len(rows)),
        "closure_defect": str(defect(rows)),
        "order_count": order_count,
        "a_plus_lower": str(average_low),
        "a_plus_upper": str(average_high),
        "a_plus_lower_decimal": decimal_lower(average_low, 12),
        "a_plus_upper_decimal": decimal_upper(average_high, 12),
        "fixed_order_minimum_lower": str(fixed_min_low),
        "fixed_order_minimum_upper": str(fixed_min_high),
        "fixed_order_minimum_lower_decimal": decimal_lower(fixed_min_low, 12),
        "fixed_order_minimum_upper_decimal": decimal_upper(fixed_min_high, 12),
        "fixed_order_maximum_lower": str(fixed_max_low),
        "fixed_order_maximum_upper": str(fixed_max_high),
        "fixed_order_maximum_lower_decimal": decimal_lower(fixed_max_low, 12),
        "fixed_order_maximum_upper_decimal": decimal_upper(fixed_max_high, 12),
        "minimum_upper_order": list(min_order),
        "maximum_lower_order": list(max_order),
        "sharp_clone_infimum_lower": str(fixed_min_low),
        "sharp_clone_infimum_upper": str(fixed_min_high),
        "sharp_clone_supremum_lower": str(fixed_max_low),
        "sharp_clone_supremum_upper": str(fixed_max_high),
        "negative_fixed_order_certified": negative_fixed_order,
        "finite_negative_clone_witness": sharp_witness,
        "every_fixed_order_certified_positive": every_fixed_order_positive,
        "every_clone_configuration_certified_positive": every_fixed_order_positive,
        "global_rank_average_minimum_lower": str(global_rank_min_lower),
        "global_rank_average_minimum_lower_decimal": decimal_lower(
            global_rank_min_lower, 12
        ),
        "best_limit_coordinate": best_coordinate,
        "profiles": profiles,
    }



def clone_reachability_certificate(entry: dict[str, object]) -> dict[str, object]:
    """Freeze the exact sub-two-fifths cloned witness in a compact certificate."""
    witness = entry["finite_negative_clone_witness"]
    assert witness is not None and witness["certified_negative"]
    lower = Fraction(witness["a_plus_lower"])
    upper = Fraction(witness["a_plus_upper"])
    closure_defect = Fraction(entry["closure_defect"])
    assert lower <= upper < 0
    assert closure_defect < Fraction(2, 5)
    multiplicities = tuple(witness["multiplicities"])
    counts = tuple(entry["coordinate_counts"])
    assert len(multiplicities) == len(counts) == entry["dimension"]
    cloned_incidence = sum(
        multiplicity * count
        for multiplicity, count in zip(multiplicities, counts)
    )
    return {
        "schema": 1,
        "verdict": "PROVED_SUB_TWO_FIFTHS_NEGATIVE_EXACT_RATIONAL",
        "proof_dependency": (
            "verify_gate_b_rational fixed-order two-sided enclosures plus "
            "PROOF.md Proposition 25 exact Plackett-Luce clone reweighting"
        ),
        "arithmetic": (
            "ordinary exact rational bounds; no float value or sampled order "
            "is used after selecting the certified fixed order"
        ),
        "displayed_supremum_convention": "cap and Reimer only",
        "normalization_convention": (
            "the core is active and separating; the cloned witness is active "
            "but non-separating because clone columns coincide"
        ),
        "frozen_targets": {
            "closure_defect_upper": "2/5",
            "a_plus_upper": "0",
        },
        "core": {
            "name": entry["name"],
            "dimension": entry["dimension"],
            "family_size": entry["family_size"],
            "family_rows": entry["family_rows"],
            "normalized": entry["normalized"],
            "admissible": entry["admissible"],
            "coordinate_counts": entry["coordinate_counts"],
            "cap": entry["cap"],
            "incidence": entry["incidence"],
            "reimer_threshold": entry["reimer_threshold"],
            "closure_defect": entry["closure_defect"],
            "negative_fixed_order": entry["minimum_upper_order"],
            "fixed_order_a_lower": entry["fixed_order_minimum_lower"],
            "fixed_order_a_upper": entry["fixed_order_minimum_upper"],
        },
        "cloned_family": {
            "symbolic_representation": (
                "replace core coordinate i by multiplicities[i] identical copies"
            ),
            "multiplicities": list(multiplicities),
            "dimension": witness["cloned_dimension"],
            "family_size": entry["family_size"],
            "active": True,
            "separating": False,
            "admissible": True,
            "cap": entry["cap"],
            "maximum_coordinate_count": max(counts),
            "incidence": cloned_incidence,
            "reimer_threshold": entry["reimer_threshold"],
            "closure_defect": entry["closure_defect"],
            "target_order_probability": witness["target_order_probability"],
            "a_plus_lower": str(lower),
            "a_plus_upper": str(upper),
            "a_plus_lower_decimal": witness["a_plus_lower_decimal"],
            "a_plus_upper_decimal": witness["a_plus_upper_decimal"],
            "repair_ratio_lower": str(-upper / closure_defect),
            "repair_ratio_upper": str(-lower / closure_defect),
            "bound_method": witness["bound_method"],
        },
    }

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidates",
        type=str,
        default=(
            "n8tiny_core,n8m15_clone_reachable,"
            "n9m15_lowdefect,n9m15_plus_new"
        ),
    )
    parser.add_argument("--max-multiplicity", type=int, default=8)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--certificate", type=Path, default=DEFAULT_CERTIFICATE
    )
    arguments = parser.parse_args()

    available = candidates()
    names = [name.strip() for name in arguments.candidates.split(",") if name.strip()]
    entries = []
    for name in names:
        rows, dimension, source = available[name]
        entries.append(
            audit_family(name, rows, dimension, source, arguments.max_multiplicity)
        )

    tiny = next((entry for entry in entries if entry["name"] == "n8tiny_core"), None)
    reachable = next(
        (
            entry
            for entry in entries
            if entry["name"] == "n8m15_clone_reachable"
        ),
        None,
    )
    low = next((entry for entry in entries if entry["name"] == "n9m15_lowdefect"), None)
    plus = next((entry for entry in entries if entry["name"] == "n9m15_plus_new"), None)
    tiny_crosses = bool(
        tiny and any(
            rung["certified_negative"]
            for rung in tiny["profiles"][str(tiny["best_limit_coordinate"])][
                "multiplicity_ladder"
            ]
        )
    )
    low_uncloneable = bool(low and low["every_clone_configuration_certified_positive"])
    plus_uncloneable = bool(
        plus and plus["every_clone_configuration_certified_positive"]
    )
    _extended, _dimension, new_column, classification = n9_nonclone_extension()
    nonclone_count = sum(
        1
        for candidate in classification["_usable"]
        if candidate not in classification["_clones"]
        and candidate not in classification["_ors"]
    )
    all_extensions_blocked = (
        low_uncloneable
        and plus_uncloneable
        and classification["usable_total"] == 8
        and nonclone_count == 1
    )
    reachable_witness = (
        reachable["finite_negative_clone_witness"] if reachable else None
    )
    sub_two_fifths_certified = bool(
        reachable
        and reachable["normalized"]
        and reachable["admissible"]
        and Fraction(reachable["closure_defect"]) < Fraction(2, 5)
        and reachable["negative_fixed_order_certified"]
        and reachable_witness
        and reachable_witness["certified_negative"]
        and Fraction(reachable_witness["a_plus_lower"])
        <= Fraction(reachable_witness["a_plus_upper"])
        < 0
    )
    verdict = (
        "SUB_TWO_FIFTHS_NEGATIVE_CLONE_WITNESS_CERTIFIED_AND_EXTENSION_CLOSURE_CLASSIFIED"
        if tiny_crosses and all_extensions_blocked and sub_two_fifths_certified
        else "INCONCLUSIVE"
    )
    report = {
        "theorem": (
            "Clone multiplicities induce the exact Plackett-Luce distribution "
            "on original first-appearance orders. The infimum and supremum over "
            "all finite clone multisets equal the minimum and maximum original "
            "fixed-order objectives; a negative fixed order has an explicit "
            "finite negative clone witness."
        ),
        "collapsed_order_probability": (
            "product_k r[pi[k]] / sum_{j>=k} r[pi[j]]"
        ),
        "extrema_status": "INFIMUM_AND_SUPREMUM; finite attainment is not claimed",
        "displayed_supremum_convention": "cap and Reimer only",
        "normalization_convention": (
            "active and separating is reported, not required by the displayed "
            "supremum; nontrivial cloning introduces duplicate columns"
        ),
        "candidates": entries,
        "sub_two_fifths_negative_clone_witness_certified": sub_two_fifths_certified,
        "sub_two_fifths_negative_clone_witness": {
            "core_name": reachable["name"],
            "core_dimension": reachable["dimension"],
            "core_rows": reachable["family_rows"],
            "core_normalized": reachable["normalized"],
            "closure_defect": reachable["closure_defect"],
            "negative_fixed_order": reachable["minimum_upper_order"],
            "fixed_order_a_lower": reachable["fixed_order_minimum_lower"],
            "fixed_order_a_upper": reachable["fixed_order_minimum_upper"],
            "clone_certificate": reachable_witness,
        } if reachable else None,
        "n8tiny_crosses_zero_by_cloning": tiny_crosses,
        "n9m15_every_clone_configuration_positive": low_uncloneable,
        "n9m15_plus_new_every_clone_configuration_positive": plus_uncloneable,
        "n9m15_extension_class": {
            "join_consistent_total": classification["join_consistent_total"],
            "usable_total": classification["usable_total"],
            "nonclone_usable_total": nonclone_count,
            "unique_nonclone_column_mask": new_column,
        },
        "n9m15_every_defect_free_extension_sequence_positive": all_extensions_blocked,
        "verdict": verdict,
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if sub_two_fifths_certified:
        certificate = clone_reachability_certificate(reachable)
        arguments.certificate.parent.mkdir(parents=True, exist_ok=True)
        arguments.certificate.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
        print(f"certificate: {arguments.certificate}")
    print(f"\n{verdict}")
    print(f"report: {arguments.report}")


if __name__ == "__main__":
    main()
