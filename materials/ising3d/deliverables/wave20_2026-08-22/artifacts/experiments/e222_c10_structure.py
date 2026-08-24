#!/usr/bin/env python3
"""Exact composition and connected-object structure at interlayer order ten.

Only finite partition algebra and exact rational formal series are used.  This
module launches no finite-lattice boxes and imports no earlier experiment.
"""

from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial
from typing import Iterable, Sequence


def record(
    checks: list[dict[str, object]], name: str, passed: bool, detail: str
) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


@lru_cache(maxsize=None)
def set_partitions(
    items: tuple[int, ...],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Independent Bell-partition generator for the ten-slot certificate."""

    if not items:
        return ((),)
    first, rest = items[0], items[1:]
    rows: list[tuple[tuple[int, ...], ...]] = []
    for partition in set_partitions(rest):
        rows.append(((first,),) + partition)
        for index in range(len(partition)):
            rows.append(
                partition[:index]
                + ((first,) + partition[index],)
                + partition[index + 1 :]
            )
    return tuple(rows)


@lru_cache(maxsize=None)
def compositions(total: int) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    return tuple(
        (first,) + tail
        for first in range(total, 0, -1)
        for tail in compositions(total - first)
    )


def denominator(composition: Sequence[int]) -> int:
    value = 1
    for part in composition:
        value *= factorial(2 * int(part))
    return value


def gap_map(composition: Sequence[int]) -> dict[int, int]:
    return {
        slot: gap
        for slot, gap in enumerate(
            (
                gap
                for gap, part in enumerate(composition)
                for _ in range(2 * int(part))
            ),
            start=1,
        )
    }


def block_gap_even(gaps: dict[int, int], block: Iterable[int]) -> bool:
    counts = Counter(gaps[slot] for slot in block)
    return all(count % 2 == 0 for count in counts.values())


def layer_degrees(composition: Sequence[int]) -> tuple[int, ...]:
    doubled = (0,) + tuple(2 * int(part) for part in composition) + (0,)
    return tuple(doubled[index] + doubled[index + 1] for index in range(len(doubled) - 1))


@lru_cache(maxsize=None)
def admissible_count_dp(group_sizes: tuple[int, ...]) -> int:
    """Count even-in-every-group partitions by a distinguished-label recurrence."""

    if not any(group_sizes):
        return 1
    first_group = next(index for index, size in enumerate(group_sizes) if size)
    choices = [tuple(range(0, size + 1, 2)) for size in group_sizes]
    choices[first_group] = tuple(
        value for value in choices[first_group] if value >= 2
    )
    total = 0
    for block_counts in product(*choices):
        ways = comb(group_sizes[first_group] - 1, block_counts[first_group] - 1)
        for index, count in enumerate(block_counts):
            if index != first_group:
                ways *= comb(group_sizes[index], count)
        remainder = tuple(
            size - count for size, count in zip(group_sizes, block_counts, strict=True)
        )
        total += ways * admissible_count_dp(remainder)
    return total


def explicit_partition_inventory(
    composition: tuple[int, ...],
    partitions: Sequence[tuple[tuple[int, ...], ...]],
) -> tuple[int, dict[tuple[int, ...], int]]:
    gaps = gap_map(composition)
    count = 0
    moebius: Counter[tuple[int, ...]] = Counter()
    for partition in partitions:
        if not all(block_gap_even(gaps, block) for block in partition):
            continue
        count += 1
        mu = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        block_class = tuple(sorted((len(block) for block in partition), reverse=True))
        moebius[block_class] += mu
    return count, dict(moebius)


def multiply_series(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    out = [Fraction(0) for _ in range(order + 1)]
    for i, first in enumerate(left):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                out[i + j] += first * second
    return tuple(out)


def divide_series(
    numerator: Sequence[Fraction], denominator_series: Sequence[Fraction], order: int
) -> tuple[Fraction, ...]:
    if denominator_series[0] == 0:
        raise ZeroDivisionError("formal denominator has zero constant term")
    out = [Fraction(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        value = numerator[degree] if degree < len(numerator) else Fraction(0)
        value -= sum(
            denominator_series[index] * out[degree - index]
            for index in range(1, degree + 1)
            if index < len(denominator_series)
        )
        out[degree] = value / denominator_series[0]
    return tuple(out)


def power_series(base: Sequence[Fraction], exponent: int, order: int) -> tuple[Fraction, ...]:
    out = (Fraction(1),) + (Fraction(0),) * order
    for _ in range(exponent):
        out = multiply_series(out, base, order)
    return out


def log_unit_series(base: Sequence[Fraction], order: int) -> tuple[Fraction, ...]:
    if base[0] != 1:
        raise ValueError("formal logarithm requires unit constant term")
    out = [Fraction(0) for _ in range(order + 1)]
    for degree in range(1, order + 1):
        correction = sum(
            (
                index * out[index] * base[degree - index]
                for index in range(1, degree)
            ),
            Fraction(0),
        )
        out[degree] = base[degree] - correction / degree
    return tuple(out)


def exact_low_order_c10() -> dict[str, object]:
    """Derive the v^0 and v^2 c10 data without an FLM computation."""

    order = 10
    exp_q = tuple(Fraction(1, factorial(index)) for index in range(order + 1))
    exp_minus_q = tuple(
        Fraction((-1) ** index, factorial(index)) for index in range(order + 1)
    )
    cosh_q = tuple((first + second) / 2 for first, second in zip(exp_q, exp_minus_q, strict=True))
    sinh_q = tuple((first - second) / 2 for first, second in zip(exp_q, exp_minus_q, strict=True))
    tanh_q = divide_series(sinh_q, cosh_q, order)
    log_cosh = log_unit_series(cosh_q, order)

    tanh_squared = multiply_series(tanh_q, tanh_q, order)
    one_minus_tanh_squared = list(-value for value in tanh_squared)
    one_minus_tanh_squared[0] += 1
    ladder_q = tuple(
        2 * value
        for value in divide_series(tanh_squared, one_minus_tanh_squared, order)
    )

    atanh_w = tuple(
        Fraction(1, degree) if degree % 2 else Fraction(0)
        for degree in range(order + 1)
    )
    conversion = {
        even: power_series(atanh_w, even, order)[order]
        for even in range(2, order + 1, 2)
    }
    converted_v0 = sum(
        conversion[even] * log_cosh[even] for even in range(2, order + 1, 2)
    )
    converted_v2 = sum(
        conversion[even] * ladder_q[even] for even in range(2, order + 1, 2)
    )
    return {
        "tag": "[THEOREM][COMPUTATION]",
        "q_to_w_weights": {str(key): str(value) for key, value in conversion.items()},
        "direct_q_v0": str(log_cosh[10]),
        "direct_q_v2": str(ladder_q[10]),
        "total_w_v0": str(converted_v0),
        "residual_w_v0": "0",
        "residual_w_v2": str(converted_v2),
        "derivation": (
            "[THEOREM] The v^0 column is log cosh(q). The all-order ladder theorem gives "
            "the residual v^2 column 2*w^2/(1-w^2); exact substitution w=tanh(q) "
            "and q=atanh(w) yields the displayed q and w coefficients."
        ),
    }


def build_c10_structure() -> tuple[dict[str, object], list[dict[str, object]]]:
    checks: list[dict[str, object]] = []
    partitions = set_partitions(tuple(range(1, 11)))
    rows: list[dict[str, object]] = []

    for composition in compositions(5):
        count, moebius = explicit_partition_inventory(composition, partitions)
        independent_count = admissible_count_dp(tuple(2 * part for part in composition))
        degrees = layer_degrees(composition)
        if len(composition) == 1:
            w10_count = admissible_count_dp((10,))
        elif len(composition) == 2:
            w10_count = admissible_count_dp((2 * composition[0],)) * admissible_count_dp(
                (2 * composition[1],)
            )
        else:
            w10_count = 0
        record(
            checks,
            f"c10 profile {composition} partition count two routes",
            count == independent_count,
            f"count={count}",
        )
        rows.append(
            {
                "composition": list(composition),
                "gap_multiplicities": [2 * part for part in composition],
                "layers": len(composition) + 1,
                "placement": str(Fraction(1, denominator(composition))),
                "ordered_gap_assignments": factorial(10) // denominator(composition),
                "layer_slot_degrees": list(degrees),
                "maximum_connected_rank": max(degrees),
                "connected_ranks_allowed_by_layer_moments": list(
                    range(2, max(degrees) + 1, 2)
                ),
                "admissible_bell_partitions": count,
                "w10_sector_monomials": w10_count,
                "moebius_by_block_class": {
                    "+".join(map(str, key)): str(value)
                    for key, value in sorted(moebius.items(), reverse=True)
                },
            }
        )

    by_composition = {tuple(row["composition"]): row for row in rows}
    reflection_ok = all(
        row["admissible_bell_partitions"]
        == by_composition[tuple(reversed(row["composition"]))][
            "admissible_bell_partitions"
        ]
        and row["w10_sector_monomials"]
        == by_composition[tuple(reversed(row["composition"]))][
            "w10_sector_monomials"
        ]
        and row["moebius_by_block_class"]
        == by_composition[tuple(reversed(row["composition"]))][
            "moebius_by_block_class"
        ]
        for row in rows
    )
    record(
        checks,
        "c10 Bell number and composition count",
        len(partitions) == 115975 and len(rows) == 16,
        f"Bell(10)={len(partitions)}, compositions={len(rows)}",
    )
    record(
        checks,
        "c10 every composition survives",
        all(row["admissible_bell_partitions"] > 0 for row in rows),
        "all sixteen exact counts are positive",
    )
    record(
        checks,
        "c10 reflection pairing",
        reflection_ok,
        "counts, W10 sectors, and Moebius block classes pair under reversal",
    )
    w10_profiles = {
        tuple(row["composition"])
        for row in rows
        if row["w10_sector_monomials"]
    }
    expected_w10_profiles = {(5,), (4, 1), (3, 2), (2, 3), (1, 4)}
    record(
        checks,
        "c10 exact W10-support profiles",
        w10_profiles == expected_w10_profiles,
        str(sorted(w10_profiles)),
    )
    record(
        checks,
        "c10 W10 square uniqueness",
        by_composition[(5,)]["placement"] == "1/3628800"
        and all(
            row["maximum_connected_rank"] < 10
            for row in rows
            if len(row["composition"]) >= 3
        ),
        "W10^2 occurs only in (5), with coefficient 1/10!",
    )

    low_order = exact_low_order_c10()
    record(
        checks,
        "c10 exact decoupled and ladder coefficients",
        low_order["direct_q_v0"] == "31/14175"
        and low_order["direct_q_v2"] == "4/14175"
        and low_order["total_w_v0"] == "1/10"
        and low_order["residual_w_v2"] == "2",
        (
            f"q(v0)={low_order['direct_q_v0']}, q(v2)={low_order['direct_q_v2']}, "
            f"w(v0)={low_order['total_w_v0']}, w(v2)={low_order['residual_w_v2']}"
        ),
    )

    if not all(bool(row["passed"]) for row in checks):
        failed = [str(row["name"]) for row in checks if not row["passed"]]
        raise AssertionError("c10 structure derivation failed: " + ", ".join(failed))

    data = {
        "tag": "[THEOREM][LEMMA][COMPUTATION][UNRESOLVED]",
        "connected_formula": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] For c_(10,q), sum over all sixteen compositions a of 5 with "
                "weight 1/prod_j(2a_j)!. In each profile sum Bell Moebius weights over "
                "partitions even in every gap and replace each induced layer moment by "
                "its even connected-partition expansion. This is an exact coefficientwise "
                "all-v formula using precisely C2=G, C4=U, C6=W, C8=W8, and C10=W10."
            ),
            "bell_number": len(partitions),
            "composition_count": len(rows),
            "maximum_layers": max(row["layers"] for row in rows),
            "connected_object_alphabet": {
                "2": "G",
                "4": "U",
                "6": "W",
                "8": "W8",
                "10": "W10",
            },
            "total_admissible_bell_partitions": sum(
                row["admissible_bell_partitions"] for row in rows
            ),
            "profiles": rows,
        },
        "w10_sector": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] Write E_10=W10+R_10. The W10 sector is "
                "(W10^2+2 W10 R_10)/10! in profile (5), and "
                "E(I_a) W10 E(J_a)/((2a)!(10-2a)!) in profiles (a,5-a), "
                "a=1,2,3,4. Every profile with at least three gaps has zero W10 sector."
            ),
            "supporting_compositions": [
                list(composition) for composition in sorted(expected_w10_profiles)
            ],
            "profile_monomial_counts": {
                str(tuple(row["composition"])): row["w10_sector_monomials"]
                for row in rows
                if row["w10_sector_monomials"]
            },
            "total_profile_indexed_monomials": sum(
                row["w10_sector_monomials"] for row in rows
            ),
            "w10_square_coefficient": "1/3628800",
        },
        "all_order_obstruction": {
            "tag": "[THEOREM]",
            "statement": (
                "[THEOREM] For every m>=1, the universal connected-correlation polynomial "
                "for c_(2m,q) contains C_[2m]^2 with coefficient 1/(2m)!. It comes uniquely "
                "from the one-gap composition (m), the one-block Bell partition, and the "
                "top connected term in each of its two layer moments. Hence no polynomial "
                "formula in connected objects of ranks below 2m can equal c_(2m,q) in the "
                "universal partition-cumulant algebra."
            ),
            "proof_scope": (
                "[LEMMA] A C_[2m] atom requires one Bell block containing all slots. Two "
                "such atoms require two layers each incident to all slots, which occurs only "
                "for the single-gap profile; its Taylor placement is 1/(2m)!."
            ),
            "limitation": (
                "[UNRESOLVED] The theorem is a universal polynomial-method obstruction. It "
                "does not prove algebraic independence after specializing to the square-lattice "
                "Ising state and does not rule out model-specific Pfaffian or non-polynomial "
                "representations of W8 or W10."
            ),
        },
        "exact_low_order_c10": low_order,
        "finite_series_scope": {
            "tag": "[UNRESOLVED]",
            "statement": (
                "[UNRESOLVED] No new c10 finite-lattice boxes were launched. Beyond the "
                "all-order v^0 and v^2 columns, coefficients v^4 through v^12 are not "
                "claimed here. The c8 coefficients through v^12 are audited separately."
            ),
        },
    }
    return data, checks


def main() -> int:
    data, checks = build_c10_structure()
    print(
        json.dumps(
            {
                "bell_number": data["connected_formula"]["bell_number"],
                "patterns": data["connected_formula"]["composition_count"],
                "checks": len(checks),
            },
            sort_keys=True,
        )
    )
    print("PASS e222 c10 exact composition structure and W10 obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
