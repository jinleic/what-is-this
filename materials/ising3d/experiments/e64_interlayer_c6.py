"""Exact sixth-order interlayer expansion of stacked square Ising layers.

Run from the repository root with
    .venv/bin/python experiments/e64_interlayer_c6.py

The cumulant variable is ``q = K_z``.  The repository high-temperature
variable is ``w = tanh(K_z)``.  This experiment keeps both conventions and
computes their exact change of variable rather than identifying them.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from sympy import Symbol, expand

from ising.interlayer import anisotropic_box_even_subgraph

SCRIPT = "experiments/e64_interlayer_c6.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "interlayer" / "c6_series.json"
V_ORDER = 12
DIRECT_ORDER = 6
STABILITY_ORDER = 8
MAX_CUMULANT = 6
MAX_LAYERS = 4

RationalSeries = tuple[Fraction, ...]
IntegerSeries = tuple[int, ...]
Series = tuple[int | Fraction, ...]


def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _integer_zero(order: int) -> list[int]:
    return [0 for _ in range(order + 1)]


def _fraction_strings(values: Iterable[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _add_scaled(
    target: list[int | Fraction],
    source: Sequence[int | Fraction],
    scale: int | Fraction = 1,
) -> None:
    exact_scale = Fraction(scale) if isinstance(scale, Fraction) else scale
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += exact_scale * coefficient


def _mul(
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
    order: int,
) -> Series:
    out: list[int | Fraction] = [0] * (order + 1)
    left_nonzero = [(degree, value) for degree, value in enumerate(left[: order + 1]) if value]
    right_nonzero = [(degree, value) for degree, value in enumerate(right[: order + 1]) if value]
    for i, first in left_nonzero:
        for j, second in right_nonzero:
            if i + j <= order:
                out[i + j] += first * second
    return tuple(out)


def _pow(base: Sequence[int | Fraction], exponent: int, order: int) -> Series:
    result: Series = (1,) + (0,) * order
    factor: Series = tuple(base[: order + 1]) + (0,) * max(0, order + 1 - len(base))
    power = int(exponent)
    while power:
        if power & 1:
            result = _mul(result, factor, order)
        power //= 2
        if power:
            factor = _mul(factor, factor, order)
    return result


def _divide(
    numerator: Sequence[int | Fraction],
    denominator: Sequence[int | Fraction],
    order: int,
) -> RationalSeries:
    den = [
        Fraction(denominator[degree]) if degree < len(denominator) else Fraction(0)
        for degree in range(order + 1)
    ]
    if not den[0]:
        raise ValueError("formal-series denominator has zero constant term")
    num = [
        Fraction(numerator[degree]) if degree < len(numerator) else Fraction(0)
        for degree in range(order + 1)
    ]
    quotient = _zero(order)
    for degree in range(order + 1):
        convolution = sum(
            (den[shift] * quotient[degree - shift] for shift in range(1, degree + 1)),
            Fraction(0),
        )
        quotient[degree] = (num[degree] - convolution) / den[0]
    return tuple(quotient)


def _square_bonds(width: int, height: int) -> tuple[tuple[int, int], ...]:
    bonds = []
    for y in range(height):
        for x in range(width):
            site = y * width + x
            if x + 1 < width:
                bonds.append((site, site + 1))
            if y + 1 < height:
                bonds.append((site, site + width))
    return tuple(bonds)


def parity_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """Return exact open-layer edge polynomials indexed by boundary mask."""

    site_count = width * height
    coefficients = [[0] * (order + 1) for _ in range(1 << site_count)]
    coefficients[0][0] = 1
    for left, right in _square_bonds(width, height):
        toggle = (1 << left) | (1 << right)
        for degree in range(order - 1, -1, -1):
            for mask, polynomial in enumerate(coefficients):
                value = polynomial[degree]
                if value:
                    coefficients[mask ^ toggle][degree + 1] += value
    return tuple(tuple(polynomial) for polynomial in coefficients)


def _set_partitions(items: tuple[int, ...]) -> Iterator[tuple[tuple[int, ...], ...]]:
    """Generate each unlabeled set partition once in canonical insertion order."""

    if not items:
        yield ()
        return
    first = items[0]
    for partition in _set_partitions(items[1:]):
        yield ((first,),) + partition
        for index, block in enumerate(partition):
            yield partition[:index] + ((first,) + block,) + partition[index + 1 :]


@lru_cache(maxsize=None)
def _partition_profile_coefficients(order_n: int) -> tuple[tuple[tuple[int, ...], int, int], ...]:
    """Compress the complete set-partition cumulant sum by block-size profile.

    Each returned row is ``(profile, partition_count, total_Moebius_weight)``.
    It is derived by enumerating set partitions, not from a hand-entered sixth
    cumulant formula.
    """

    counts: Counter[tuple[int, ...]] = Counter()
    weights: Counter[tuple[int, ...]] = Counter()
    for partition in _set_partitions(tuple(range(order_n))):
        profile = tuple(sorted((len(block) for block in partition), reverse=True))
        coefficient = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        counts[profile] += 1
        weights[profile] += coefficient
    return tuple(
        (profile, counts[profile], weights[profile])
        for profile in sorted(counts, reverse=True)
    )


def _cumulant_from_moments(
    moments: dict[int, RationalSeries], order_n: int, v_order: int
) -> RationalSeries:
    result = _zero(v_order)
    unit: RationalSeries = (Fraction(1),) + (Fraction(0),) * v_order
    for profile, _count, coefficient in _partition_profile_coefficients(order_n):
        term: Series = unit
        for block_size in profile:
            term = _mul(term, moments[block_size], v_order)
        _add_scaled(result, term, coefficient)
    return tuple(result)


def _symbolic_partition_certificate() -> dict[str, object]:
    """Certify the sixth-cumulant and vertical-gap bookkeeping symbolically."""

    partitions = tuple(_set_partitions(tuple(range(6))))
    profile_rows = []
    for profile, count, weight in _partition_profile_coefficients(6):
        profile_rows.append(
            {
                "block_sizes": list(profile),
                "partition_count": count,
                "total_mobius_weight": weight,
                "survives_zero_mean_symmetry": all(size % 2 == 0 for size in profile),
            }
        )

    gap_rows = []
    for gap_profile in ((6,), (4, 2), (2, 4), (2, 2, 2)):
        assignment = tuple(
            gap for gap, count in enumerate(gap_profile) for _ in range(count)
        )
        counts: Counter[tuple[int, ...]] = Counter()
        weights: Counter[tuple[int, ...]] = Counter()
        for partition in partitions:
            # A block moment vanishes if an outer-layer cut carries an odd
            # number of its spin insertions.  Testing every block and gap is
            # the symbolic zero-field specialization of the layer moments.
            admissible = all(
                all(
                    sum(assignment[index] == gap for index in block) % 2 == 0
                    for gap in range(len(gap_profile))
                )
                for block in partition
            )
            if not admissible:
                continue
            shape = tuple(sorted((len(block) for block in partition), reverse=True))
            coefficient = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
            counts[shape] += 1
            weights[shape] += coefficient
        gap_rows.append(
            {
                "ordered_gap_counts": list(gap_profile),
                "admissible_partition_counts": {
                    "+".join(map(str, shape)): counts[shape] for shape in sorted(counts, reverse=True)
                },
                "mobius_weight_sums": {
                    "+".join(map(str, shape)): weights[shape] for shape in sorted(weights, reverse=True)
                },
                "all_moments_one_cancellation": sum(weights.values()),
            }
        )

    # A mixed cumulant of two independent nonempty gap groups must vanish.
    # Use a distinct symbol for every subset moment, so this check does not
    # accidentally identify inequivalent spin correlations.
    left_labels = frozenset((0, 1))
    right_labels = frozenset((2, 3, 4, 5))
    subset_symbols: dict[tuple[str, tuple[int, ...]], object] = {}

    def subset_moment(family: str, labels: tuple[int, ...]):
        if not labels:
            return 1
        key = (family, labels)
        if key not in subset_symbols:
            subset_symbols[key] = Symbol(f"{family}_{''.join(map(str, labels))}")
        return subset_symbols[key]

    disconnected = 0
    all_labels = left_labels | right_labels
    for partition in partitions:
        term = 1
        for block in partition:
            left = tuple(sorted(set(block) & left_labels))
            right = tuple(sorted(set(block) & right_labels))
            term *= subset_moment("A", left) * subset_moment("B", right)
        coefficient = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        disconnected += coefficient * term
    disconnected_zero = expand(disconnected) == 0 and all_labels == frozenset(range(6))

    c2 = Symbol("c2")
    c4 = Symbol("c4")
    c6 = Symbol("c6")
    w = Symbol("w")
    q_truncated = w + w**3 / 3 + w**5 / 5
    transformed = expand(c2 * q_truncated**2 + c4 * q_truncated**4 + c6 * q_truncated**6)
    conversion_residual = expand(
        transformed.coeff(w, 6) - (c6 + Fraction(4, 3) * c4 + Fraction(23, 45) * c2)
    )

    passed = (
        len(partitions) == 203
        and sum(row["partition_count"] for row in profile_rows) == 203
        and disconnected_zero
        and conversion_residual == 0
        and [row["ordered_gap_counts"] for row in gap_rows]
        == [[6], [4, 2], [2, 4], [2, 2, 2]]
    )
    return {
        "tag": "[LEMMA]",
        "passed": bool(passed),
        "bell_number_B6": len(partitions),
        "set_partition_profiles": profile_rows,
        "zero_mean_surviving_formula": "kappa6=mu6-15*mu4*mu2+30*mu2^3",
        "connected_adjacent_gap_profiles": gap_rows,
        "gap_profile_scope": (
            "after anchoring the lowest used gap, the only nonzero connected profiles are "
            "6, 4+2, 2+4, and 2+2+2; skipped gaps factor into independent families"
        ),
        "disconnected_gap_symbolic_cancellation": bool(disconnected_zero),
        "change_of_variable": "c6_w=c6_q+(4/3)c4_q+(23/45)c2_q",
    }


@lru_cache(maxsize=None)
def _q_transform_weights(vertical_bonds: int) -> tuple[tuple[Fraction, ...], ...]:
    """Return n! [q^n] cosh(q)^B tanh(q)^d for 0 <= n,d <= 6."""

    q_order = MAX_CUMULANT
    cosh_q: RationalSeries = (
        Fraction(1),
        Fraction(0),
        Fraction(1, 2),
        Fraction(0),
        Fraction(1, 24),
        Fraction(0),
        Fraction(1, 720),
    )
    tanh_q: RationalSeries = (
        Fraction(0),
        Fraction(1),
        Fraction(0),
        Fraction(-1, 3),
        Fraction(0),
        Fraction(2, 15),
        Fraction(0),
    )
    prefactor = _pow(cosh_q, vertical_bonds, q_order)
    rows: list[tuple[Fraction, ...]] = []
    for order_n in range(q_order + 1):
        row = []
        for vertical_degree in range(q_order + 1):
            factor = _mul(prefactor, _pow(tanh_q, vertical_degree, q_order), q_order)
            row.append(Fraction(factor[order_n]) * factorial(order_n))
        rows.append(tuple(row))
    return tuple(rows)


def _active_boundary_masks(
    parity: Sequence[IntegerSeries], max_weight: int
) -> tuple[tuple[tuple[int, ...], ...], bool]:
    by_weight: list[list[int]] = [[] for _ in range(max_weight + 1)]
    odd_zero = True
    for mask, polynomial in enumerate(parity):
        weight = mask.bit_count()
        if weight <= max_weight and any(polynomial):
            by_weight[weight].append(mask)
            if weight % 2:
                odd_zero = False
    return tuple(tuple(masks) for masks in by_weight), odd_zero


def slab_subset_polynomials(
    parity: Sequence[IntegerSeries], layers: int, order: int, max_weight: int = 6
) -> tuple[IntegerSeries, ...]:
    """Build finite-slab vertical-subset columns from independent layer moments.

    For gap masks ``S_1,...,S_(h-1)``, layer independence gives the numerator
    ``P[S_1] P[S_1 xor S_2] ... P[S_(h-2) xor S_(h-1)] P[S_(h-1)]``.
    A transfer over the last mask sums all such sequences without enumerating
    site tuples or omitting an adjacent-gap pattern.
    """

    if layers < 1 or layers > MAX_LAYERS:
        raise ValueError("layers must lie between one and four")
    if layers == 1:
        columns = [tuple(_integer_zero(order)) for _ in range(max_weight + 1)]
        columns[0] = tuple(parity[0])
        return tuple(columns)

    by_weight, odd_zero = _active_boundary_masks(parity, max_weight)
    if not odd_zero:
        raise AssertionError("an odd boundary polynomial is nonzero")
    active = tuple(
        (mask, weight)
        for weight, masks in enumerate(by_weight)
        for mask in masks
    )

    # The first outer layer contributes P[S_1].
    state: dict[tuple[int, int], IntegerSeries] = {
        (mask, weight): tuple(parity[mask]) for mask, weight in active
    }
    for _gap in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int]] = {}
        for (previous_mask, degree_w), accumulated in state.items():
            for current_mask, current_weight in active:
                total_weight = degree_w + current_weight
                if total_weight > max_weight:
                    continue
                middle = parity[previous_mask ^ current_mask]
                if not any(middle):
                    continue
                contribution = _mul(accumulated, middle, order)
                if not any(contribution):
                    continue
                key = (current_mask, total_weight)
                target = updated.get(key)
                if target is None:
                    target = _integer_zero(order)
                    updated[key] = target
                _add_scaled(target, contribution)
        state = {key: tuple(values) for key, values in updated.items()}

    columns = [_integer_zero(order) for _ in range(max_weight + 1)]
    for (last_mask, degree_w), accumulated in state.items():
        _add_scaled(columns[degree_w], _mul(accumulated, parity[last_mask], order))
    return tuple(tuple(column) for column in columns)


def _raw_moments_from_layer_moments(
    columns: Sequence[IntegerSeries],
    plane_partition: IntegerSeries,
    layers: int,
    area: int,
    order: int,
) -> tuple[dict[int, RationalSeries], bool]:
    """Recover raw moments of the vertical-bond sum from 2D layer moments."""

    denominator = _pow(plane_partition, layers, order)
    normalized = tuple(_divide(column, denominator, order) for column in columns)
    odd_columns_zero = all(not any(normalized[degree]) for degree in (1, 3, 5))
    weights = _q_transform_weights(area * (layers - 1))
    moments: dict[int, RationalSeries] = {}
    for order_n in range(1, MAX_CUMULANT + 1):
        values = _zero(order)
        for vertical_degree in range(MAX_CUMULANT + 1):
            coefficient = weights[order_n][vertical_degree]
            if coefficient:
                _add_scaled(values, normalized[vertical_degree], coefficient)
        moments[order_n] = tuple(values)
    return moments, odd_columns_zero


@lru_cache(maxsize=None)
def _rectangle_cumulant_components_canonical(
    width: int, height: int, order: int
) -> tuple[tuple[str, RationalSeries], ...]:
    parity = parity_polynomials(width, height, order)
    area = width * height
    finite: dict[int, dict[int, RationalSeries]] = {1: {n: tuple(_zero(order)) for n in range(1, 7)}}
    odd_columns_zero = True
    for layers in range(2, MAX_LAYERS + 1):
        columns = slab_subset_polynomials(parity, layers, order, MAX_CUMULANT)
        moments, box_odd_zero = _raw_moments_from_layer_moments(
            columns, parity[0], layers, area, order
        )
        odd_columns_zero = odd_columns_zero and box_odd_zero
        finite[layers] = {
            order_n: tuple(
                value / factorial(order_n)
                for value in _cumulant_from_moments(moments, order_n, order)
            )
            for order_n in range(1, MAX_CUMULANT + 1)
        }

    vertical_weights: dict[tuple[int, int], RationalSeries] = {}
    for order_n in range(1, MAX_CUMULANT + 1):
        for layers in range(1, MAX_LAYERS + 1):
            values = list(finite[layers][order_n])
            for sublayers in range(1, layers):
                placements = layers - sublayers + 1
                _add_scaled(values, vertical_weights[(order_n, sublayers)], -placements)
            vertical_weights[(order_n, layers)] = tuple(values)

    result = {
        f"c{order_n}_h{layers}": vertical_weights[(order_n, layers)]
        for order_n in range(1, MAX_CUMULANT + 1)
        for layers in range(1, MAX_LAYERS + 1)
    }
    result["odd_columns_zero"] = (Fraction(int(odd_columns_zero)),) + (Fraction(0),) * order
    return tuple(sorted(result.items()))


def rectangle_cumulant_components(
    width: int, height: int, order: int
) -> dict[str, RationalSeries]:
    first, second = sorted((int(width), int(height)))
    return dict(_rectangle_cumulant_components_canonical(first, second, int(order)))


def _plane_shapes(order: int, bound_slack: int = 0) -> tuple[tuple[int, int], ...]:
    span_budget = order // 2 + bound_slack
    return tuple(
        sorted(
            (
                (width, height)
                for width in range(1, span_budget + 2)
                for height in range(1, span_budget + 2)
                if (width - 1) + (height - 1) <= span_budget
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )


def cumulant_flm(
    order: int, bound_slack: int = 0
) -> tuple[dict[str, RationalSeries], dict[tuple[int, int], dict[str, RationalSeries]], tuple[tuple[int, int], ...]]:
    """Planar rectangular FLM of the layer-moment cumulant coefficients."""

    shapes = _plane_shapes(order, bound_slack)
    quantity_names = tuple(
        f"c{order_n}_h{layers}"
        for order_n in range(1, MAX_CUMULANT + 1)
        for layers in range(1, MAX_LAYERS + 1)
    )
    weights: dict[tuple[int, int], dict[str, RationalSeries]] = {}
    bulk = {name: _zero(order) for name in quantity_names}
    for width, height in shapes:
        finite = rectangle_cumulant_components(width, height, order)
        values = {name: list(finite[name]) for name in quantity_names}
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for name in quantity_names:
                    _add_scaled(values[name], subweight[name], -placements)
        exact = {name: tuple(series_values) for name, series_values in values.items()}
        weights[(width, height)] = exact
        for name in quantity_names:
            _add_scaled(bulk[name], exact[name])
    return {name: tuple(values) for name, values in bulk.items()}, weights, shapes


def _sum_components(
    bulk: dict[str, RationalSeries], order_n: int, layers: Iterable[int], order: int
) -> RationalSeries:
    total = _zero(order)
    for layer_count in layers:
        _add_scaled(total, bulk[f"c{order_n}_h{layer_count}"])
    return tuple(total)


def _log_even_columns(
    polynomial: Sequence[Sequence[int | Fraction]], order: int
) -> tuple[RationalSeries, RationalSeries, RationalSeries, bool]:
    columns = tuple(
        tuple(Fraction(polynomial[v_degree][w_degree]) for v_degree in range(order + 1))
        for w_degree in range(MAX_CUMULANT + 1)
    )
    odd_zero = all(not any(columns[degree]) for degree in (1, 3, 5))
    ratio2 = _divide(columns[2], columns[0], order)
    ratio4 = _divide(columns[4], columns[0], order)
    ratio6 = _divide(columns[6], columns[0], order)
    ratio2_square = _mul(ratio2, ratio2, order)
    log4 = tuple(
        ratio4[degree] - Fraction(ratio2_square[degree], 2)
        for degree in range(order + 1)
    )
    ratio2_ratio4 = _mul(ratio2, ratio4, order)
    ratio2_cube = _mul(ratio2_square, ratio2, order)
    log6 = tuple(
        ratio6[degree]
        - ratio2_ratio4[degree]
        + Fraction(ratio2_cube[degree], 3)
        for degree in range(order + 1)
    )
    return ratio2, log4, log6, odd_zero


@lru_cache(maxsize=None)
def _direct_box_log(
    width: int, height: int, layers: int, order: int
) -> tuple[RationalSeries, RationalSeries, RationalSeries, bool]:
    width, height = sorted((width, height))
    if layers == 1:
        zero = tuple(_zero(order))
        return zero, zero, zero, True
    polynomial = anisotropic_box_even_subgraph(
        (width, height, layers), order, MAX_CUMULANT
    )
    return _log_even_columns(polynomial, order)


def direct_spin_dos_flm(
    order: int,
) -> tuple[
    tuple[RationalSeries, RationalSeries, RationalSeries],
    dict[tuple[int, int, int], tuple[RationalSeries, RationalSeries, RationalSeries]],
    bool,
    tuple[tuple[int, int, int], ...],
]:
    """Independent exact anisotropic FLM from three-dimensional spin DOS."""

    plane_shapes = _plane_shapes(order)
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in plane_shapes
                for layers in range(1, MAX_LAYERS + 1)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    weights: dict[
        tuple[int, int, int], tuple[RationalSeries, RationalSeries, RationalSeries]
    ] = {}
    bulk = tuple(_zero(order) for _ in range(3))
    odd_zero = True
    for width, height, layers in shapes:
        c2, c4, c6, box_odd_zero = _direct_box_log(width, height, layers, order)
        odd_zero = odd_zero and box_odd_zero
        values = [list(c2), list(c4), list(c6)]
        for (sub_width, sub_height, sublayers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sublayers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sublayers + 1)
                )
                for quantity in range(3):
                    _add_scaled(values[quantity], subweight[quantity], -placements)
        exact = tuple(tuple(series_values) for series_values in values)
        weights[(width, height, layers)] = exact
        for quantity in range(3):
            _add_scaled(bulk[quantity], exact[quantity])
    return tuple(tuple(values) for values in bulk), weights, odd_zero, shapes


def _sum_direct_height(
    weights: dict[tuple[int, int, int], tuple[RationalSeries, ...]],
    layers: int,
    quantity: int,
    order: int,
) -> RationalSeries:
    total = _zero(order)
    for shape, values in weights.items():
        if shape[2] == layers:
            _add_scaled(total, values[quantity])
    return tuple(total)


def _subtract_constant(values: RationalSeries, constant: Fraction) -> RationalSeries:
    result = list(values)
    result[0] -= constant
    return tuple(result)


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    checks: list[dict[str, object]] = []

    symbolic = _symbolic_partition_certificate()
    _record(
        checks,
        "complete sixth-cumulant set partitions",
        bool(symbolic["passed"]),
        "[LEMMA] all 203 partitions, adjacent profiles 6/4+2/2+4/2+2+2, and q-to-w mixing were generated exactly",
    )

    bulk, _weights, shapes = cumulant_flm(V_ORDER)
    c2_q = _sum_components(bulk, 2, (2,), V_ORDER)
    c4_q_h2 = bulk["c4_h2"]
    c4_q_h3 = bulk["c4_h3"]
    c4_q = _sum_components(bulk, 4, (2, 3), V_ORDER)
    c6_q_h2 = bulk["c6_h2"]
    c6_q_h3 = bulk["c6_h3"]
    c6_q_h4 = bulk["c6_h4"]
    c6_q = _sum_components(bulk, 6, (2, 3, 4), V_ORDER)

    c4_w_total = tuple(
        c4_q[degree] + Fraction(2, 3) * c2_q[degree]
        for degree in range(V_ORDER + 1)
    )
    c6_w_total = tuple(
        c6_q[degree]
        + Fraction(4, 3) * c4_q[degree]
        + Fraction(23, 45) * c2_q[degree]
        for degree in range(V_ORDER + 1)
    )
    c2_w_residual = _subtract_constant(c2_q, Fraction(1, 2))
    c4_w_residual = _subtract_constant(c4_w_total, Fraction(1, 4))
    c6_w_residual = _subtract_constant(c6_w_total, Fraction(1, 6))

    c6_w_height2_residual = _subtract_constant(
        tuple(
            c6_q_h2[degree]
            + Fraction(4, 3) * c4_q_h2[degree]
            + Fraction(23, 45) * c2_q[degree]
            for degree in range(V_ORDER + 1)
        ),
        Fraction(1, 6),
    )
    c6_w_height3 = tuple(
        c6_q_h3[degree] + Fraction(4, 3) * c4_q_h3[degree]
        for degree in range(V_ORDER + 1)
    )
    c6_w_height4 = c6_q_h4

    known_c2 = (
        Fraction(1, 2), Fraction(0), Fraction(2), Fraction(0), Fraction(18),
        Fraction(0), Fraction(118), Fraction(0), Fraction(778), Fraction(0),
        Fraction(4978), Fraction(0), Fraction(31398),
    )
    known_c4_q = (
        Fraction(-1, 12), Fraction(0), Fraction(2, 3), Fraction(0), Fraction(51),
        Fraction(0), Fraction(2914, 3), Fraction(0), Fraction(41113, 3), Fraction(0),
        Fraction(481894, 3), Fraction(0), Fraction(1679961),
    )
    _record(
        checks,
        "lower-order reproduction",
        c2_q == known_c2 and c4_q == known_c4_q,
        "[COMPUTATION] the generic partition engine reproduces the exact c2 and c4 prefixes through v^12",
    )

    odd_series_zero = all(
        not any(bulk[f"c{order_n}_h{layers}"])
        for order_n in (1, 3, 5)
        for layers in range(1, MAX_LAYERS + 1)
    )
    forbidden_extents_zero = (
        not any(bulk["c2_h3"])
        and not any(bulk["c2_h4"])
        and not any(bulk["c4_h4"])
    )
    _record(
        checks,
        "odd cumulants and vertical support",
        odd_series_zero and forbidden_extents_zero,
        "[LEMMA] c3=c5=0 and orders 2/4 have no linked weights beyond heights 2/3 through v^12",
    )

    direct_bulk, direct_weights, direct_odd_zero, direct_shapes = direct_spin_dos_flm(
        DIRECT_ORDER
    )
    direct_c2_residual, direct_c4_residual, direct_c6_residual = direct_bulk
    _record(
        checks,
        "independent spin-DOS agreement",
        direct_odd_zero
        and direct_c2_residual == c2_w_residual[: DIRECT_ORDER + 1]
        and direct_c4_residual == c4_w_residual[: DIRECT_ORDER + 1]
        and direct_c6_residual == c6_w_residual[: DIRECT_ORDER + 1],
        "[COMPUTATION] direct anisotropic 3D spin-DOS FLM agrees for c2/c4/c6 through v^6 and has w^1,w^3,w^5=0",
    )

    direct_height2 = _sum_direct_height(
        direct_weights, 2, 2, DIRECT_ORDER
    )
    direct_height3 = _sum_direct_height(
        direct_weights, 3, 2, DIRECT_ORDER
    )
    direct_height4 = _sum_direct_height(
        direct_weights, 4, 2, DIRECT_ORDER
    )
    _record(
        checks,
        "c6 exact vertical-extent agreement",
        direct_height2 == c6_w_height2_residual[: DIRECT_ORDER + 1]
        and direct_height3 == c6_w_height3[: DIRECT_ORDER + 1]
        and direct_height4 == c6_w_height4[: DIRECT_ORDER + 1],
        "[COMPUTATION] height-2/3/4 slab weights match same-gap/two-gap/three-gap cumulant components through v^6",
    )

    _record(
        checks,
        "zero-inplane control",
        c2_q[0] == Fraction(1, 2)
        and c4_q[0] == Fraction(-1, 12)
        and c6_q[0] == Fraction(1, 45)
        and c6_w_total[0] == Fraction(1, 6)
        and c6_w_residual[0] == 0,
        "[COMPUTATION] K=0 exactly reproduces log cosh(q) and -log(1-w^2)/2 through sixth order",
    )
    _record(
        checks,
        "bipartite v parity",
        all(
            series_values[degree] == 0
            for series_values in (c2_q, c4_q, c6_q, c6_w_residual)
            for degree in range(1, V_ORDER + 1, 2)
        ),
        "[COMPUTATION] every odd v coefficient vanishes through v^12",
    )

    # A one-unit enlargement is inexpensive at order eight and directly checks
    # the projected-span cutoff used for the reported order-twelve calculation.
    minimal_stability, _, minimal_shapes = cumulant_flm(STABILITY_ORDER)
    enlarged_stability, _, enlarged_shapes = cumulant_flm(
        STABILITY_ORDER, bound_slack=1
    )
    _record(
        checks,
        "finite-lattice bound stability",
        minimal_stability == enlarged_stability,
        "[COMPUTATION] adding all rectangles with one extra span unit changes no cumulant component through v^8",
    )

    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError("one or more exact c6 checks failed")

    order_provenance = [
        {
            "v_degree": degree,
            "c6_direct_q": str(c6_q[degree]),
            "c6_total_w": str(c6_w_total[degree]),
            "c6_residual_w": str(c6_w_residual[degree]),
            "same_gap_q_height2": str(c6_q_h2[degree]),
            "two_adjacent_gaps_q_height3": str(c6_q_h3[degree]),
            "three_adjacent_gaps_q_height4": str(c6_q_h4[degree]),
            "direct_spin_dos_residual_w": (
                str(direct_c6_residual[degree]) if degree <= DIRECT_ORDER else None
            ),
            "two_route_agreement_available": degree <= DIRECT_ORDER,
            "maximum_required_inplane_span": degree // 2,
        }
        for degree in range(V_ORDER + 1)
    ]

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "method": (
                "exact 2D boundary-polynomial layer moments plus complete set-partition cumulants "
                "and planar FLM; independent exact anisotropic 3D spin-DOS FLM through v^6"
            ),
        },
        "data": {
            "claim": (
                "[COMPUTATION] Exact finite high-temperature c6 prefix through v^12, "
                "independently checked through v^6; no all-orders coefficient claim is made."
            ),
            "variables": {
                "v": "tanh(K)",
                "q": "K_z, direct Boltzmann coupling for cumulants",
                "w": "tanh(K_z), interlayer high-temperature variable",
                "conversion": (
                    "q=atanh(w); c4_w=c4_q+(2/3)c2_q; "
                    "c6_w=c6_q+(4/3)c4_q+(23/45)c2_q"
                ),
            },
            "v_order": V_ORDER,
            "independent_check_v_order": DIRECT_ORDER,
            "symbolic_sixth_cumulant": symbolic,
            "reproduced_lower_orders": {
                "tag": "[COMPUTATION]",
                "c2_total_q_or_w": _fraction_strings(c2_q),
                "c2_residual_w": _fraction_strings(c2_w_residual),
                "c4_direct_q": _fraction_strings(c4_q),
                "c4_total_w": _fraction_strings(c4_w_total),
                "c4_residual_w": _fraction_strings(c4_w_residual),
            },
            "c6_direct_coupling_q": {
                "tag": "[COMPUTATION]",
                "full_infinite_stack": _fraction_strings(c6_q),
                "same_gap_height2": _fraction_strings(c6_q_h2),
                "two_adjacent_gaps_height3": _fraction_strings(c6_q_h3),
                "three_adjacent_gaps_height4": _fraction_strings(c6_q_h4),
            },
            "c6_wave_variable_w_tanh_Kz": {
                "tag": "[COMPUTATION]",
                "total": _fraction_strings(c6_w_total),
                "residual_after_log_cosh": _fraction_strings(c6_w_residual),
                "residual_by_exact_vertical_extent": {
                    "2": _fraction_strings(c6_w_height2_residual),
                    "3": _fraction_strings(c6_w_height3),
                    "4": _fraction_strings(c6_w_height4),
                },
                "nonzero_residual_coefficients": {
                    str(degree): str(value)
                    for degree, value in enumerate(c6_w_residual)
                    if value
                },
            },
            "odd_interlayer_orders": {
                "tag": "[LEMMA]",
                "c3_through_v12": ["0"] * (V_ORDER + 1),
                "c5_through_v12": ["0"] * (V_ORDER + 1),
                "reason": (
                    "every layer boundary has even cardinality; equivalently the open stack is "
                    "mapped q to -q by flipping alternate layers"
                ),
            },
            "independent_spin_dos_route": {
                "tag": "[COMPUTATION]",
                "method": (
                    "exact joint spin density of states for each open 3D box, transformed directly "
                    "to P(v,w), followed by [w^2],[w^4],[w^6] log extraction and 3D Moebius inversion"
                ),
                "box_count": len(direct_shapes),
                "boxes": [list(shape) for shape in direct_shapes],
                "c2_residual_w_through_v6": _fraction_strings(direct_c2_residual),
                "c4_residual_w_through_v6": _fraction_strings(direct_c4_residual),
                "c6_residual_w_through_v6": _fraction_strings(direct_c6_residual),
                "c6_residual_by_exact_vertical_extent": {
                    "2": _fraction_strings(direct_height2),
                    "3": _fraction_strings(direct_height3),
                    "4": _fraction_strings(direct_height4),
                },
                "odd_w_columns_zero": direct_odd_zero,
                "agreement": True,
            },
            "finite_lattice_certificate": {
                "tag": "[LEMMA]",
                "2d_rectangle_count": len(shapes),
                "2d_rectangles": [list(shape) for shape in shapes],
                "maximum_vertical_extent": MAX_LAYERS,
                "completeness_bound": (
                    "a linked sixth-order term crosses each used layer cut a positive even number "
                    "of times, hence at most three consecutive gaps (height four); a projected "
                    "connected even multigraph spanning a-by-b uses at least "
                    "2*((a-1)+(b-1)) in-plane edges, hence span at most six through v^12"
                ),
                "stability_check_order": STABILITY_ORDER,
                "minimal_stability_rectangle_count": len(minimal_shapes),
                "enlarged_stability_rectangle_count": len(enlarged_shapes),
            },
            "order_by_order_provenance": order_provenance,
            "scope": (
                "[COMPUTATION] These exact coefficients are a finite expansion around K=K_z=0. "
                "They do not determine the convergence radius, justify evaluation at K_c, or solve the 3D Ising model."
            ),
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e64_interlayer_c6: {error}")
        raise
