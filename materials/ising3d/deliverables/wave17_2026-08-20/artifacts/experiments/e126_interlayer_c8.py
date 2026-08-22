"""Exact eighth-order interlayer expansion of stacked square Ising layers.

This experiment computes the next interlayer coefficient

    f_3D(K, K_z) = f_2D(K) + c2 w^2 + c4 w^4 + c6 w^6 + c8 w^8 + O(w^10),
    w = tanh(K_z),  v = tanh(K),

as an exact high-temperature series in ``v``, generalizing the certified
sixth-order pipeline (``experiments/e64_interlayer_c6.py``,
``experiments/e99_c6_closed.py``) with the same conventions:

* primary route: exact open-rectangle boundary-mask layer polynomials,
  complete Bell(8) = 4140 set-partition cumulants of the vertical-bond sum,
  vertical Moebius inversion through six layers, and planar rectangular
  finite-lattice inversion;
* the vertical-subset transfer is accelerated by the rectangle symmetry
  group and validated against the plain transfer;
* second route (composition machinery): exact ``[w^d] log`` columns of the
  slab even-subgraph polynomial and three-dimensional box Moebius inversion;
* third route (engine independent): the repository anisotropic spin density
  of states through ``v^6``;
* exact ``q -> w`` conversion with every ``atanh`` weight derived from the
  exponential series (nothing hand entered);
* complete Bell(8) connected-moment term inventory for ``c8``: which of the
  2D connected objects G / U4 / W6 / W8 appear and in how many terms,
  anchored by regenerating the certified sixth-order 13-term formula.

Everything is ``int`` / ``fractions.Fraction``; no floating point enters any
result.
"""

from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial
from pathlib import Path
from typing import Iterable, Sequence

from ising.interlayer import anisotropic_box_even_subgraph

SCRIPT = "experiments/e126_interlayer_c8.py"
ROOT = Path(__file__).resolve().parents[1]
C6_SERIES = ROOT / "results" / "interlayer" / "c6_series.json"
C6_CLOSED = ROOT / "results" / "interlayer" / "c6_closed.json"
RESULT = ROOT / "results" / "interlayer" / "c8_series.json"

V_ORDER = 12        # in-plane order of every exact series below
Q_ORDER = 8         # eighth vertical cumulant
MAX_WEIGHT = 8      # vertical weight columns of every slab
MAX_LAYERS = 5      # maximal contributing height of c8 (pattern (1,1,1,1))
EXTENT_LAYERS = 6   # one extra height for the vertical extent control
SPIN_DOS_ORDER = 6  # depth of the engine-independent spin-DOS cross-check
IDENTITY_ORDER = 8  # v-order of the sampled connected-form verification
STABILITY_SHAPES = ((8, 1), (7, 2))  # partial span-7 stability control
WALL_BUDGET_SECONDS = 1800


RationalSeries = tuple[Fraction, ...]
IntegerSeries = tuple[int, ...]
Series = tuple[int | Fraction, ...]

# ---------------------------------------------------------------------------
# exact truncated power series (int and Fraction)
# ---------------------------------------------------------------------------
def _int_zero(order: int) -> list[int]:
    return [0] * (order + 1)


def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _add_scaled(
    target: list[int | Fraction], source: Sequence[int | Fraction], scale: int | Fraction = 1
) -> None:
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += scale * coefficient


def _mul_int(left: Sequence[int], right: Sequence[int], order: int) -> list[int]:
    out = _int_zero(order)
    for first, value in enumerate(left[: order + 1]):
        if not value:
            continue
        for second, other in enumerate(right[: order + 1 - first]):
            if other:
                out[first + second] += value * other
    return out


def _mul(
    left: Sequence[int | Fraction], right: Sequence[int | Fraction], order: int
) -> list[Fraction]:
    out = _zero(order)
    for first, value in enumerate(left[: order + 1]):
        if not value:
            continue
        for second, other in enumerate(right[: order + 1 - first]):
            if other:
                out[first + second] += Fraction(value) * Fraction(other)
    return out


def _pow_series(base: Sequence[int | Fraction], exponent: int, order: int) -> list[Fraction]:
    result = [Fraction(1)] + [Fraction(0)] * order
    for _ in range(exponent):
        result = _mul(result, base, order)
    return result


def _divide_unit(numerator: Sequence[int], denominator: Sequence[int], order: int) -> list[int]:
    if denominator[0] != 1:
        raise ValueError("formal division requires unit constant denominator")
    quotient = _int_zero(order)
    for degree in range(order + 1):
        value = numerator[degree] if degree < len(numerator) else 0
        for lower in range(degree):
            value -= quotient[lower] * denominator[degree - lower]
        quotient[degree] = value
    return quotient


def _divide(
    numerator: Sequence[int | Fraction], denominator: Sequence[int | Fraction], order: int
) -> RationalSeries:
    if not denominator or not denominator[0]:
        raise ValueError("formal division requires nonzero constant denominator")
    quotient = _zero(order)
    for degree in range(order + 1):
        value = (
            Fraction(numerator[degree])
            if degree < len(numerator) and numerator[degree]
            else Fraction(0)
        )
        for lower in range(degree):
            if quotient[lower] and degree - lower < len(denominator):
                value -= quotient[lower] * Fraction(denominator[degree - lower])
        quotient[degree] = value / Fraction(denominator[0])
    return quotient


# ---------------------------------------------------------------------------
# exact q-series machinery: cosh, sinh, tanh, atanh from the exponential
# ---------------------------------------------------------------------------
def _exp_series(order: int, sign: int = 1) -> tuple[Fraction, ...]:
    return tuple(Fraction(sign**k, factorial(k)) for k in range(order + 1))


COSH_Q = tuple(
    (Fraction(1, factorial(k)) + Fraction((-1) ** k, factorial(k))) / 2
    for k in range(Q_ORDER + 1)
)
SINH_Q = tuple(
    (Fraction(1, factorial(k)) - Fraction((-1) ** k, factorial(k))) / 2
    for k in range(Q_ORDER + 1)
)
TANH_Q = tuple(_divide(SINH_Q, COSH_Q, Q_ORDER))
# q(w)=atanh(w) is the unique formal primitive of 1/(1-w^2).  The
# denominator is obtained from the exact identity (tanh q)'=1-tanh(q)^2,
# so no conversion coefficient is hand entered.
_ATANH_DENOMINATOR = [Fraction(0)] * (Q_ORDER + 1)
_ATANH_DENOMINATOR[0] = Fraction(1)
if Q_ORDER >= 2:
    _ATANH_DENOMINATOR[2] = Fraction(-1)
_ATANH_DERIVATIVE = _divide(
    [Fraction(1)] + [Fraction(0)] * Q_ORDER, _ATANH_DENOMINATOR, Q_ORDER
)
ATANH_W = [Fraction(0)] * (Q_ORDER + 1)
for _k in range(1, Q_ORDER + 1):
    ATANH_W[_k] = _ATANH_DERIVATIVE[_k - 1] / _k


@lru_cache(maxsize=None)
def _q_transform_weights(vertical_bonds: int) -> tuple[tuple[Fraction, ...], ...]:
    """Rows ``n`` of ``n! [q^n] cosh(q)^B tanh(q)^d`` for ``0<=n,d<=Q_ORDER``."""

    prefactor = _pow_series(COSH_Q, vertical_bonds, Q_ORDER)
    rows: list[tuple[Fraction, ...]] = []
    for order_n in range(Q_ORDER + 1):
        row: list[Fraction] = []
        for vertical_degree in range(Q_ORDER + 1):
            factor = _mul(prefactor, _pow_series(TANH_Q, vertical_degree, Q_ORDER), Q_ORDER)
            row.append(Fraction(factor[order_n]) * factorial(order_n))
        rows.append(tuple(row))
    return tuple(rows)


def _conversion_coefficients(order_n: int = Q_ORDER) -> dict[int, Fraction]:
    """``[w^order_n] atanh(w)^n`` for every even ``n <= order_n``, exactly."""

    out: dict[int, Fraction] = {}
    for power in range(2, order_n + 1, 2):
        series = _pow_series(ATANH_W, power, order_n)
        out[power] = Fraction(series[order_n])
    return out


# ---------------------------------------------------------------------------
# open-rectangle boundary-mask polynomials P_S(v)  (c4/c6 convention)
# ---------------------------------------------------------------------------
def _square_bonds(width: int, height: int) -> tuple[tuple[int, int], ...]:
    bonds: list[tuple[int, int]] = []
    for y in range(height):
        for x in range(width):
            site = y * width + x
            if x + 1 < width:
                bonds.append((site, site + 1))
            if y + 1 < height:
                bonds.append((site, site + width))
    return tuple(bonds)


@lru_cache(maxsize=None)
def boundary_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """``P_S(v)`` for every odd-degree mask ``S`` of one open 2D rectangle."""

    sites = width * height
    values = [[0] * (order + 1) for _ in range(1 << sites)]
    values[0][0] = 1
    for first, second in _square_bonds(width, height):
        toggle = (1 << first) | (1 << second)
        for degree in range(order - 1, -1, -1):
            for mask in range(1 << sites):
                coefficient = values[mask][degree]
                if coefficient:
                    values[mask ^ toggle][degree + 1] += coefficient
    return tuple(tuple(row) for row in values)


# ---------------------------------------------------------------------------
# rectangle symmetry group and orbit-accelerated slab transfer
# ---------------------------------------------------------------------------
def rectangle_group(width: int, height: int) -> tuple[dict[int, int], ...]:
    """Site permutations of the open rectangle that preserve its bond set."""

    def idx(x: int, y: int) -> int:
        return y * width + x

    coordinates = [
        lambda x, y: (x, y),
        lambda x, y: (x, height - 1 - y),
        lambda x, y: (width - 1 - x, y),
        lambda x, y: (width - 1 - x, height - 1 - y),
    ]
    if width == height:
        coordinates += [
            lambda x, y: (y, x),
            lambda x, y: (y, width - 1 - x),
            lambda x, y: (height - 1 - y, x),
            lambda x, y: (height - 1 - y, width - 1 - x),
        ]
    bond_set = set(_square_bonds(width, height))
    permutations_out: list[dict[int, int]] = []
    for transform in coordinates:
        mapping: dict[int, int] = {}
        for y in range(height):
            for x in range(width):
                u, v = transform(x, y)
                mapping[idx(x, y)] = idx(u, v)
        if len(set(mapping.values())) != width * height:
            continue
        mapped = set(
            (min(mapping[a], mapping[b]), max(mapping[a], mapping[b])) for a, b in bond_set
        )
        if mapped == bond_set:
            permutations_out.append(mapping)
    if not permutations_out:
        raise AssertionError("rectangle group is never empty")
    return tuple(permutations_out)


class Orbits:
    """Canonical representatives and orbit sizes of subsets under the group."""

    def __init__(self, width: int, height: int) -> None:
        self.group = rectangle_group(width, height)
        sites = width * height
        self.rep = [0] * (1 << sites)
        self.size = [0] * (1 << sites)
        for mask in range(1 << sites):
            orbit = {self._permute(mask, perm) for perm in self.group}
            representative = min(orbit)
            self.rep[mask] = representative
            self.size[mask] = len(orbit)

    def _permute(self, mask: int, perm: dict[int, int]) -> int:
        out = 0
        for source, target in perm.items():
            if (mask >> source) & 1:
                out |= 1 << target
        return out


def _active_masks(
    parity: Sequence[IntegerSeries], max_weight: int
) -> tuple[tuple[int, int], ...]:
    """Even masks of weight <= max_weight with a nonzero truncated series.

    Dropping heavier-boundary masks is exact through the working order: a
    slab path through an intermediate mask S costs at least 2*delta(S)
    in-plane edges, where delta(S) is the minimal T-join size of S.
    """
    active: list[tuple[int, int]] = []
    for mask, series in enumerate(parity):
        weight = mask.bit_count()
        if weight <= max_weight and weight % 2 == 0 and any(series):
            active.append((mask, weight))
    return tuple(active)


def slab_columns_plain(
    parity: Sequence[IntegerSeries], layers: int, order: int, max_weight: int
) -> tuple[IntegerSeries, ...]:
    """Reference vertical-subset column transfer without orbit compression."""

    active = _active_masks(parity, max_weight)
    if layers == 1:
        columns = [_int_zero(order) for _ in range(max_weight + 1)]
        columns[0] = list(parity[0])
        return tuple(tuple(column) for column in columns)
    state: dict[tuple[int, int], list[int]] = {
        (mask, weight): list(parity[mask]) for mask, weight in active
    }
    for _gap in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int]] = {}
        for (previous, previous_weight), accumulated in state.items():
            for current, current_weight in active:
                total = previous_weight + current_weight
                if total > max_weight:
                    continue
                middle = parity[previous ^ current]
                if not any(middle):
                    continue
                contribution = _mul_int(accumulated, middle, order)
                if not any(contribution):
                    continue
                key = (current, total)
                if key not in updated:
                    updated[key] = contribution
                else:
                    _add_scaled(updated[key], contribution)
        state = updated
    columns = [_int_zero(order) for _ in range(max_weight + 1)]
    for (last, vertical_degree), accumulated in state.items():
        _add_scaled(columns[vertical_degree], _mul_int(accumulated, parity[last], order))
    return tuple(tuple(column) for column in columns)


def slab_columns(
    parity: Sequence[IntegerSeries], layers: int, order: int, max_weight: int, orbits: Orbits
) -> tuple[IntegerSeries, ...]:
    """Orbit-accelerated vertical-subset columns; equals ``slab_columns_plain``.

    The accumulated state is invariant under the rectangle group because the
    in-plane coupling is isotropic, so targets are accumulated at orbit
    representatives while sources are expanded over their orbits; the final
    assembly multiplies each representative by its orbit size.
    """

    active = _active_masks(parity, max_weight)
    if layers == 1:
        columns = [_int_zero(order) for _ in range(max_weight + 1)]
        columns[0] = list(parity[0])
        return tuple(tuple(column) for column in columns)
    rep, size = orbits.rep, orbits.size
    targets_by_weight: dict[int, list[int]] = defaultdict(list)
    state: dict[tuple[int, int], list[int]] = {}
    members: dict[int, list[int]] = defaultdict(list)
    for mask, weight in active:
        representative = rep[mask]
        targets_by_weight[weight].append(representative)
        members[representative].append(mask)
        key = (representative, weight)
        if key not in state:
            state[key] = list(parity[mask])
    for weight in targets_by_weight:
        targets_by_weight[weight] = sorted(set(targets_by_weight[weight]))
    for _gap in range(1, layers - 1):
        updated: dict[tuple[int, int], list[int]] = {}
        for (source, used), accumulated in state.items():
            for mask in members[source]:
                for target_weight, targets in targets_by_weight.items():
                    if used + target_weight > max_weight:
                        continue
                    for target in targets:
                        middle = parity[mask ^ target]
                        if not any(middle):
                            continue
                        contribution = _mul_int(accumulated, middle, order)
                        if not any(contribution):
                            continue
                        key = (target, used + target_weight)
                        if key not in updated:
                            updated[key] = contribution
                        else:
                            _add_scaled(updated[key], contribution)
        state = updated
    columns = [_int_zero(order) for _ in range(max_weight + 1)]
    for (last, vertical_degree), accumulated in state.items():
        tail = _mul_int(accumulated, parity[last], order)
        _add_scaled(columns[vertical_degree], tail, size[last])
    return tuple(tuple(column) for column in columns)


# ---------------------------------------------------------------------------
# primary route: raw moments, Bell cumulants, vertical + planar inversion
# ---------------------------------------------------------------------------
def set_partitions(items: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Every unlabeled set partition once, in canonical insertion order."""

    if not items:
        return ((),)

    def build(remaining: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
        if not remaining:
            return ((),)
        head, tail = remaining[0], remaining[1:]
        rows: list[tuple[tuple[int, ...], ...]] = []
        for partition in build(tail):
            rows.append(((head,),) + partition)
            for index in range(len(partition)):
                rows.append(
                    partition[:index] + ((head,) + partition[index],) + partition[index + 1 :]
                )
        return tuple(rows)

    return build(items)


@lru_cache(maxsize=None)
def profile_coefficients(order_n: int) -> tuple[tuple[tuple[int, ...], int], ...]:
    """Total Moebius weight of every block-size profile of ``[order_n]``."""

    weights: Counter[tuple[int, ...]] = Counter()
    for partition in set_partitions(tuple(range(order_n))):
        profile = tuple(sorted((len(block) for block in partition), reverse=True))
        weights[profile] += (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
    return tuple(sorted(weights.items(), reverse=True))


def cumulant_from_moments(
    moments: dict[int, RationalSeries], order_n: int, order: int
) -> RationalSeries:
    values = _zero(order)
    unit: list[Fraction] = [Fraction(1)] + [Fraction(0)] * order
    for profile, coefficient in profile_coefficients(order_n):
        term: list[Fraction] = unit
        for size in profile:
            term = _mul(term, moments[size], order)
        _add_scaled(values, term, coefficient)
    return tuple(values)


def raw_moments(
    columns: Sequence[IntegerSeries],
    partition: IntegerSeries,
    layers: int,
    area: int,
    order: int,
) -> dict[int, RationalSeries]:
    """Raw moments ``mu_n = n![q^n] cosh(q)^B Z_h(tanh q)`` for ``n<=8``."""

    denominator = [Fraction(value) for value in _pow_series(partition, layers, order)]
    normalized = [
        _divide([Fraction(value) for value in columns[d]], denominator, order)
        for d in range(MAX_WEIGHT + 1)
    ]
    weights = _q_transform_weights(area * (layers - 1))
    moments: dict[int, RationalSeries] = {}
    for order_n in range(1, Q_ORDER + 1):
        values = _zero(order)
        for vertical_degree in range(MAX_WEIGHT + 1):
            scale = weights[order_n][vertical_degree]
            if scale:
                _add_scaled(values, normalized[vertical_degree], scale)
        moments[order_n] = tuple(values)
    return moments


@lru_cache(maxsize=None)
def rectangle_components(
    width: int, height: int, order: int, max_layers: int
) -> dict[tuple[int, int], RationalSeries]:
    """Exact-height cumulant weights ``c_n(h)/n!`` of one open rectangle."""

    first, second = sorted((width, height))
    parity = boundary_polynomials(first, second, order)
    orbits = Orbits(first, second)
    area = first * second
    finite: dict[int, dict[int, RationalSeries]] = {
        1: {order_n: tuple(_zero(order)) for order_n in range(1, Q_ORDER + 1)}
    }
    for layers in range(2, max_layers + 1):
        columns = slab_columns(parity, layers, order, MAX_WEIGHT, orbits)
        moments = raw_moments(columns, parity[0], layers, area, order)
        finite[layers] = {
            order_n: tuple(
                value / factorial(order_n)
                for value in cumulant_from_moments(moments, order_n, order)
            )
            for order_n in range(1, Q_ORDER + 1)
        }
    vertical: dict[tuple[int, int], RationalSeries] = {}
    for order_n in range(1, Q_ORDER + 1):
        for layers in range(1, max_layers + 1):
            values = list(finite[layers][order_n])
            for lower in range(1, layers):
                _add_scaled(values, vertical[(order_n, lower)], -(layers - lower + 1))
            vertical[(order_n, layers)] = tuple(values)
    return vertical


def plane_shapes(order: int, slack: int = 0) -> tuple[tuple[int, int], ...]:
    span = order // 2 + slack
    return tuple(
        sorted(
            (
                (width, height)
                for width in range(1, span + 2)
                for height in range(1, span + 2)
                if (width - 1) + (height - 1) <= span
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )


def cumulant_flm(
    order: int, slack: int = 0, max_layers: int = EXTENT_LAYERS, shapes=None
) -> tuple[dict[str, RationalSeries], dict[tuple[int, int], dict[str, RationalSeries]], tuple]:
    """Planar rectangular finite-lattice inversion of every cumulant weight."""

    shape_tuple = tuple(shapes) if shapes else plane_shapes(order, slack)
    names = tuple(
        f"c{order_n}_h{layers}"
        for order_n in range(1, Q_ORDER + 1)
        for layers in range(1, max_layers + 1)
    )
    weights: dict[tuple[int, int], dict[str, RationalSeries]] = {}
    bulk = {name: _zero(order) for name in names}
    for width, height in shape_tuple:
        vertical = rectangle_components(width, height, order, max_layers)
        local = {
            f"c{order_n}_h{layers}": list(vertical[(order_n, layers)])
            for order_n in range(1, Q_ORDER + 1)
            for layers in range(1, max_layers + 1)
        }
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for name in names:
                    _add_scaled(local[name], subweight[name], -placements)
        exact = {name: tuple(values) for name, values in local.items()}
        weights[(width, height)] = exact
        for name in names:
            _add_scaled(bulk[name], exact[name])
    return {name: tuple(values) for name, values in bulk.items()}, weights, shape_tuple


def _sum_series(
    bulk: dict[str, RationalSeries], order_n: int, heights: Iterable[int], order: int
) -> RationalSeries:
    total = _zero(order)
    for layers in heights:
        _add_scaled(total, bulk[f"c{order_n}_h{layers}"])
    return tuple(total)

def _subtract_constant(values: Sequence[Fraction], constant: Fraction) -> RationalSeries:
    out = list(values)
    out[0] -= constant
    return tuple(out)


# ---------------------------------------------------------------------------
# second route: [w^d] log columns and three-dimensional box inversion
# ---------------------------------------------------------------------------
def log_columns(
    columns: Sequence[IntegerSeries], order: int, w_order: int = Q_ORDER
) -> dict[int, RationalSeries]:
    """Exact ``[w^d] log(P/P_0)`` by the Newton recursion."""

    ratios = {
        degree: _divide(
            [Fraction(value) for value in columns[degree]],
            [Fraction(value) for value in columns[0]],
            order,
        )
        for degree in range(1, w_order + 1)
    }
    out: dict[int, RationalSeries] = {}
    for degree in range(1, w_order + 1):
        accumulated = _zero(order)
        _add_scaled(accumulated, ratios[degree], degree)
        for lower in range(1, degree):
            _add_scaled(accumulated, _mul(out[lower], ratios[degree - lower], order), -lower)
        out[degree] = tuple(value / degree for value in accumulated)
    return out


def wlog_bulk(
    order: int,
    max_layers: int = EXTENT_LAYERS,
    w_order: int = Q_ORDER,
    shapes=None,
    slab_fn=None,
) -> tuple[dict[int, RationalSeries], dict[int, dict[int, RationalSeries]], dict]:
    """Three-dimensional box Moebius inversion of the ``w``-log columns."""

    shape_tuple = tuple(shapes) if shapes else plane_shapes(order)
    if slab_fn is None:
        cache: dict[tuple[int, int], tuple[IntegerSeries, ...]] = {}

        def slab_fn(width: int, height: int, layers: int) -> tuple[IntegerSeries, ...]:
            key = (width, height)
            if key not in cache:
                parity = boundary_polynomials(width, height, order)
                orbits = Orbits(width, height)
                cache[key] = (parity, orbits)
            parity, orbits = cache[key]
            return slab_columns(parity, layers, order, w_order, orbits)

    boxes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in shape_tuple
                for layers in range(1, max_layers + 1)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    degrees = tuple(range(2, w_order + 1, 2))
    weights: dict[tuple[int, int, int], dict[int, RationalSeries]] = {}
    bulk = {degree: _zero(order) for degree in degrees}
    by_layers = {
        layers: {degree: _zero(order) for degree in degrees}
        for layers in range(1, max_layers + 1)
    }
    for width, height, layers in boxes:
        columns = slab_fn(width, height, layers)
        logs = log_columns(columns, order, w_order)
        values = {degree: list(logs[degree]) for degree in degrees}
        for (sub_w, sub_h, sub_l), subweight in weights.items():
            if sub_w <= width and sub_h <= height and sub_l <= layers:
                placements = (
                    (width - sub_w + 1) * (height - sub_h + 1) * (layers - sub_l + 1)
                )
                for degree in degrees:
                    _add_scaled(values[degree], subweight[degree], -placements)
        exact = {degree: tuple(value) for degree, value in values.items()}
        weights[(width, height, layers)] = exact
        for degree in degrees:
            _add_scaled(bulk[degree], exact[degree])
            _add_scaled(by_layers[layers][degree], exact[degree])
    return (
        {degree: tuple(values) for degree, values in bulk.items()},
        {
            layers: {degree: tuple(values) for degree, values in rows.items()}
            for layers, rows in by_layers.items()
        },
        weights,
    )


# ---------------------------------------------------------------------------
# third route: anisotropic spin-DOS three-dimensional inversion
# ---------------------------------------------------------------------------
def spin_dos_box_columns(
    shape: tuple[int, int, int], order: int, w_order: int = Q_ORDER
) -> tuple[RationalSeries, ...]:
    polynomial = anisotropic_box_even_subgraph(shape, order, w_order)
    return tuple(
        tuple(Fraction(polynomial[v_degree][w_degree]) for v_degree in range(order + 1))
        for w_degree in range(w_order + 1)
    )


def spin_dos_bulk(
    order: int, heights: range, shapes=None
) -> tuple[dict[int, RationalSeries], dict[int, dict[int, RationalSeries]], bool, dict]:
    shape_tuple = tuple(shapes) if shapes else plane_shapes(order)
    boxes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in shape_tuple
                for layers in heights
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    degrees = tuple(range(2, Q_ORDER + 1, 2))
    weights: dict[tuple[int, int, int], dict[int, RationalSeries]] = {}
    bulk = {degree: _zero(order) for degree in degrees}
    by_layers = {layers: {degree: _zero(order) for degree in degrees} for layers in heights}
    odd_zero = True
    for width, height, layers in boxes:
        canonical_width, canonical_height = sorted((width, height))
        columns = spin_dos_box_columns((canonical_width, canonical_height, layers), order)
        odd_zero = odd_zero and all(not any(columns[d]) for d in (1, 3, 5, 7))
        logs = log_columns(columns, order)
        values = {degree: list(logs[degree]) for degree in degrees}
        for (sub_w, sub_h, sub_l), subweight in weights.items():
            if sub_w <= width and sub_h <= height and sub_l <= layers:
                placements = (
                    (width - sub_w + 1) * (height - sub_h + 1) * (layers - sub_l + 1)
                )
                for degree in degrees:
                    _add_scaled(values[degree], subweight[degree], -placements)
        exact = {degree: tuple(value) for degree, value in values.items()}
        weights[(width, height, layers)] = exact
        for degree in degrees:
            _add_scaled(bulk[degree], exact[degree])
            _add_scaled(by_layers[layers][degree], exact[degree])
    return (
        {degree: tuple(values) for degree, values in bulk.items()},
        {
            layers: {degree: tuple(values) for degree, values in rows.items()}
            for layers, rows in by_layers.items()
        },
        odd_zero,
        weights,
    )


# ---------------------------------------------------------------------------
# complete Bell(8) connected-moment inventory (G / U / W6 / W8 alphabet)
# ---------------------------------------------------------------------------
Atom = tuple[str, tuple[int, ...]]
Monomial = tuple[Atom, ...]
Polynomial = dict[Monomial, Fraction]

_ATOM_KIND = {2: "G", 4: "U", 6: "W", 8: "W8"}


def _atom(kind: str, slots: Iterable[int]) -> Atom:
    return kind, tuple(sorted(slots))


def _monomial(*atoms: Atom) -> Monomial:
    return tuple(sorted(atoms))


def _poly_add(target: Polynomial, source: Polynomial, scale: Fraction = Fraction(1)) -> None:
    for monomial, coefficient in source.items():
        target[monomial] = target.get(monomial, Fraction(0)) + scale * coefficient


def _poly_mul_into(target: Polynomial, left: Polynomial, right: Polynomial) -> None:
    for first, coefficient in left.items():
        for second, other in right.items():
            monomial = _monomial(*first, *second)
            target[monomial] = target.get(monomial, Fraction(0)) + coefficient * other


def _even_block_partitions(slots: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Every partition of ``slots`` into even blocks."""

    if not slots:
        return ((),)
    first, rest = slots[0], slots[1:]
    rows: list[tuple[tuple[int, ...], ...]] = []
    for size in (2, 4, 6, 8):
        if size > len(slots):
            break
        for partners in combinations(rest, size - 1):
            block = (first,) + partners
            remaining = tuple(slot for slot in rest if slot not in partners)
            for tail in _even_block_partitions(remaining):
                rows.append((block,) + tail)
    return tuple(rows)


@lru_cache(maxsize=None)
def connected_expansion(slots: tuple[int, ...]) -> Polynomial:
    """``M_I`` as a sum over even partitions of products of connected atoms."""

    if not slots:
        return {(): Fraction(1)}
    if len(slots) % 2:
        return {}
    out: Polynomial = {}
    for partition in _even_block_partitions(slots):
        monomial = _monomial(*(_atom(_ATOM_KIND[len(block)], block) for block in partition))
        out[monomial] = out.get(monomial, Fraction(0)) + Fraction(1)
    return out


def compositions(total: int) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    out: list[tuple[int, ...]] = []
    for first in range(total, 0, -1):
        for tail in compositions(total - first):
            out.append((first,) + tail)
    return tuple(out)


def composition_gap_map(composition: Sequence[int]) -> dict[int, int]:
    gaps: dict[int, int] = {}
    slot = 1
    for gap_index, part in enumerate(composition):
        for _ in range(2 * int(part)):
            gaps[slot] = gap_index
            slot += 1
    return gaps


def block_layer_slot_sets(gaps: dict[int, int], block: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Slot tuples carried by one block, one per layer it touches."""

    per_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        per_gap[gaps[slot]].append(slot)
    top_gap = max(gaps.values())
    sets: list[tuple[int, ...]] = []
    for layer in range(top_gap + 2):
        members = sorted(per_gap.get(layer - 1, []) + per_gap.get(layer, []))
        if members:
            sets.append(tuple(members))
    return tuple(sets)

def indexed_block_layer_slot_sets(
    gaps: dict[int, int], block: Sequence[int]
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Layer-labelled moment factors for a literal formal Bell polynomial."""

    per_gap: dict[int, list[int]] = defaultdict(list)
    for slot in block:
        per_gap[gaps[slot]].append(slot)
    out: list[tuple[int, tuple[int, ...]]] = []
    for layer in range(max(gaps.values()) + 2):
        slots = tuple(sorted(per_gap.get(layer - 1, []) + per_gap.get(layer, [])))
        if slots:
            out.append((layer, slots))
    return tuple(out)


def normalized_gap_multisets(count: int) -> tuple[tuple[int, ...], ...]:
    """All nondecreasing normalized gap multisets, with no selection predicate."""

    rows: list[tuple[int, ...]] = []

    def walk(prefix: tuple[int, ...]) -> None:
        if len(prefix) == count:
            rows.append(prefix)
            return
        if not prefix:
            walk((0,))
            return
        for value in range(prefix[-1], count):
            walk(prefix + (value,))

    walk(())
    return tuple(rows)


def literal_profile_moment_polynomial(
    gap_values: Sequence[int],
) -> dict[tuple[tuple[int, tuple[int, ...]], ...], int]:
    """The literal layer-labelled Bell/moment polynomial of one gap multiset."""

    gaps = {slot + 1: value for slot, value in enumerate(gap_values)}
    raw: dict[tuple[tuple[int, tuple[int, ...]], ...], int] = {}
    for partition in set_partitions(tuple(range(1, len(gap_values) + 1))):
        if not all(block_gap_even(gaps, block) for block in partition):
            continue
        coefficient = (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        monomial = tuple(
            sorted(
                factor
                for block in partition
                for factor in indexed_block_layer_slot_sets(gaps, block)
            )
        )
        raw[monomial] = raw.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in raw.items() if coefficient}


def placement_from_ordered_assignments(composition: Sequence[int]) -> tuple[Fraction, int]:
    """Derive the bulk placement factor by explicit labelled-gap assignments."""

    gap_multiset = tuple(
        gap for gap, part in enumerate(composition) for _ in range(2 * int(part))
    )
    ordered_assignments = len(set(permutations(gap_multiset)))
    return Fraction(ordered_assignments, factorial(len(gap_multiset))), ordered_assignments


def block_gap_even(gaps: dict[int, int], block: Sequence[int]) -> bool:
    return all(value % 2 == 0 for value in Counter(gaps[slot] for slot in block).values())


def block_layer_even(gaps: dict[int, int], block: Sequence[int]) -> bool:
    return all(len(slots) % 2 == 0 for slots in block_layer_slot_sets(gaps, block))


def profile_connected_polynomial(
    gaps: dict[int, int], placement: Fraction, substitute_during_build: bool = True
) -> tuple[Polynomial, dict[str, object]]:
    """Collected connected-cumulant polynomial of one gap pattern."""

    n = len(gaps)
    out: Polynomial = {}
    raw: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {}
    admissible = 0
    moebius_by_class: Counter[tuple[int, ...]] = Counter()
    partition_list = set_partitions(tuple(range(1, n + 1)))
    for partition in partition_list:
        if not all(block_gap_even(gaps, block) for block in partition):
            continue
        weight = placement * Fraction(
            (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        )
        admissible += 1
        moebius_by_class[
            tuple(sorted((len(block) for block in partition), reverse=True))
        ] += int(weight / placement)
        term_raw: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {(): Fraction(1)}
        if substitute_during_build:
            term: Polynomial = {(): Fraction(1)}
        for block in partition:
            layer_sets = block_layer_slot_sets(gaps, block)
            block_raw: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {(): Fraction(1)}
            for slots in layer_sets:
                fresh: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {}
                for monomial, coefficient in block_raw.items():
                    key = monomial + ((0, slots),)
                    fresh[key] = fresh.get(key, Fraction(0)) + coefficient
                block_raw = fresh
            merged: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction] = {}
            for left, coefficient in term_raw.items():
                for right, other in block_raw.items():
                    key = tuple(sorted(left + right))
                    merged[key] = merged.get(key, Fraction(0)) + coefficient * other
            term_raw = merged
            if substitute_during_build:
                block_term: Polynomial = {(): Fraction(1)}
                for slots in layer_sets:
                    product: Polynomial = {}
                    _poly_mul_into(product, block_term, connected_expansion(slots))
                    block_term = {k: v for k, v in product.items() if v}
                    if not block_term:
                        break
                product = {}
                _poly_mul_into(product, term, block_term)
                term = {k: v for k, v in product.items() if v}
                if not term:
                    break
        if substitute_during_build and not term:
            continue
        if not term_raw:
            continue
        for monomial, coefficient in term_raw.items():
            raw[monomial] = raw.get(monomial, Fraction(0)) + weight * coefficient
        if substitute_during_build:
            for monomial, coefficient in term.items():
                out[monomial] = out.get(monomial, Fraction(0)) + weight * coefficient
    return {k: v for k, v in out.items() if v}, {
        "admissible_partitions": admissible,
        "bell_number": len(partition_list),
        "moebius_by_class": dict(moebius_by_class),
        "raw_monomials": {k: v for k, v in raw.items() if v},
    }


def substitute_raw(raw: dict[tuple[tuple[int, tuple[int, ...]], ...], Fraction]) -> Polynomial:
    """Independent second path: substitute after collecting the raw form."""

    out: Polynomial = {}
    for monomial, coefficient in raw.items():
        term: Polynomial = {(): Fraction(1)}
        for _layer, slots in monomial:
            product: Polynomial = {}
            _poly_mul_into(product, term, connected_expansion(slots))
            term = {k: v for k, v in product.items() if v}
            if not term:
                break
        if not term:
            continue
        _poly_add(out, term, coefficient)
    return {k: v for k, v in out.items() if v}


def atom_kinds(monomial: Monomial) -> tuple[str, ...]:
    return tuple(sorted(kind for kind, _slots in monomial))


def monomial_even_incidence(monomial: Monomial, n: int) -> bool:
    incidence: Counter[int] = Counter(slot for _kind, slots in monomial for slot in slots)
    return all(incidence[slot] % 2 == 0 for slot in range(1, n + 1))


def monomial_anchor_connected(monomial: Monomial, n: int) -> bool:
    reached = {1}
    changed = True
    while changed:
        changed = False
        for _kind, slots in monomial:
            if reached.intersection(slots) and not set(slots).issubset(reached):
                reached.update(slots)
                changed = True
    return reached == set(range(1, n + 1))


_PERMUTATION_CACHE: dict[int, tuple[dict[int, int], ...]] = {}


def _slot_permutations(n: int) -> tuple[dict[int, int], ...]:
    if n not in _PERMUTATION_CACHE:
        _PERMUTATION_CACHE[n] = tuple(
            {1: 1, **dict(zip(range(2, n + 1), ordering))}
            for ordering in permutations(range(2, n + 1))
        )
    return _PERMUTATION_CACHE[n]


def canonical_monomial(monomial: Monomial, n: int) -> Monomial:
    return min(
        _monomial(*(_atom(kind, (mapping[slot] for slot in slots)) for kind, slots in monomial))
        for mapping in _slot_permutations(n)
    )


def canonical_form(polynomial: Polynomial, n: int) -> Polynomial:
    out: Polynomial = {}
    for monomial, coefficient in polynomial.items():
        key = canonical_monomial(monomial, n)
        out[key] = out.get(key, Fraction(0)) + coefficient
    return {k: v for k, v in out.items() if v}


# ---------------------------------------------------------------------------
# sampled numeric verification of the connected inventory
# ---------------------------------------------------------------------------
def _deterministic_tuples(sites: int, count: int, anchor: int = 0) -> tuple[tuple[int, ...], ...]:
    """Reproducible pseudo-random anchored site tuples (LCG, fixed seed)."""

    state = 0x1234_5678
    out: list[tuple[int, ...]] = []
    while len(out) < count:
        state = (6364136223846793005 * state + 1442695040888963407) % (1 << 63)
        tuple_out = [anchor]
        for _slot in range(7):
            state = (6364136223846793005 * state + 1442695040888963407) % (1 << 63)
            tuple_out.append(state % sites)
        out.append(tuple(tuple_out))
    return tuple(out)


def layer_moment_values(width: int, height: int, order: int) -> tuple[RationalSeries, ...]:
    parity = boundary_polynomials(width, height, order)
    denominator = [Fraction(v) for v in parity[0]]
    return tuple(_divide([Fraction(v) for v in row], denominator, order) for row in parity)


def joint_cumulant_by_recursion(
    gaps: dict[int, int],
    positions: Sequence[int],
    moments: Sequence[RationalSeries],
    order: int,
) -> RationalSeries:
    """Joint cumulant of the vertical bonds via the subset recursion."""

    n = len(positions)
    site_of = [1 << positions[i] for i in range(n)]

    def layer_masks(subset: int) -> dict[int, int]:
        masks: dict[int, int] = defaultdict(int)
        for slot in range(n):
            if (subset >> slot) & 1:
                gap = gaps[slot + 1]
                masks[gap] ^= site_of[slot]
                masks[gap + 1] ^= site_of[slot]
        return masks

    joint: dict[int, RationalSeries] = {}
    for subset in range(1, 1 << n):
        value = [Fraction(1)] + [Fraction(0)] * order
        for _gap, mask in layer_masks(subset).items():
            value = _mul(value, moments[mask], order)
            if not any(value):
                break
        joint[subset] = value
    joint[0] = [Fraction(1)] + [Fraction(0)] * order

    cumulant = {0: joint[0]}

    def kappa(subset: int) -> RationalSeries:
        if subset in cumulant:
            return cumulant[subset]
        lowest_bit = subset & -subset
        rest = subset ^ lowest_bit
        value = list(joint[subset])
        chosen = rest
        while True:
            block = chosen | lowest_bit
            if block != subset:
                _add_scaled(value, _mul(kappa(block), joint[subset ^ block], order), -1)
            if chosen == 0:
                break
            chosen = (chosen - 1) & rest
        cumulant[subset] = value
        return value

    return tuple(kappa((1 << n) - 1))


def evaluate_collected(
    polynomial: Polynomial,
    gaps: dict[int, int],
    positions: Sequence[int],
    moments: Sequence[RationalSeries],
    order: int,
) -> RationalSeries:
    """Evaluate a collected connected polynomial on one anchored tuple."""

    n = len(positions)
    site_of = [1 << positions[i] for i in range(n)]
    masks: dict[tuple[int, ...], int] = {}

    def mask_of(slots: tuple[int, ...]) -> int:
        if slots not in masks:
            mask = 0
            for slot in slots:
                mask ^= site_of[slot - 1]
            masks[slots] = mask
        return masks[slots]

    connected: dict[tuple[int, ...], RationalSeries] = {}

    def connected_value(slots: tuple[int, ...]) -> RationalSeries:
        if slots in connected:
            return connected[slots]
        if not slots:
            return tuple([Fraction(1)] + [Fraction(0)] * order)
        value = list(moments[mask_of(slots)])
        first, rest = slots[0], slots[1:]
        for size in range(2, len(slots), 2):
            for partners in combinations(rest, size - 1):
                block = (first,) + partners
                remaining = tuple(slot for slot in rest if slot not in partners)
                _add_scaled(
                    value,
                    _mul(connected_value(block), moments[mask_of(remaining)], order),
                    -1,
                )
        connected[slots] = tuple(value)
        return connected[slots]

    total = _zero(order)
    for monomial, coefficient in polynomial.items():
        product: list[Fraction] = [Fraction(1)] + [Fraction(0)] * order
        for kind, slots in monomial:
            product = _mul(product, connected_value(slots), order)
            if not any(product):
                break
        _add_scaled(total, product, coefficient)
    return tuple(total)


def evaluate_raw_partitions(
    gaps: dict[int, int],
    placement: Fraction,
    positions: Sequence[int],
    moments: Sequence[RationalSeries],
    order: int,
) -> RationalSeries:
    """Literal partition-sum evaluation of one pattern's cumulant."""

    n = len(positions)
    site_of = [1 << positions[i] for i in range(n)]
    total = _zero(order)
    for partition in set_partitions(tuple(range(1, n + 1))):
        if not all(block_gap_even(gaps, block) for block in partition):
            continue
        weight = placement * Fraction(
            (-1) ** (len(partition) - 1) * factorial(len(partition) - 1)
        )
        product: list[Fraction] = [Fraction(1)] + [Fraction(0)] * order
        for block in partition:
            for slots in block_layer_slot_sets(gaps, block):
                mask = 0
                for slot in slots:
                    mask ^= site_of[slot - 1]
                product = _mul(product, moments[mask], order)
                if not any(product):
                    break
            if not any(product):
                break
        _add_scaled(total, product, weight)
    return tuple(total)


# ---------------------------------------------------------------------------
# artifact helpers
# ---------------------------------------------------------------------------
def _strings(values: Iterable[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _fractions(values: Sequence[str]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def _polynomial_digest(polynomials: dict[str, Polynomial]) -> str:
    digest = hashlib.sha256()
    for name in sorted(polynomials):
        digest.update(name.encode("utf-8"))
        for monomial, coefficient in sorted(polynomials[name].items()):
            digest.update(
                repr((monomial, coefficient)).encode("utf-8")
            )
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    started = time.monotonic()
    checks: list[dict[str, object]] = []

    # -- sanity of the layer engine ------------------------------------
    shapes = plane_shapes(V_ORDER)
    engine_ok = True
    for width, height in shapes:
        parity = boundary_polynomials(width, height, V_ORDER)
        n_bonds = len(_square_bonds(width, height))
        engine_ok &= all(
            not any(parity[mask])
            for mask in range(1 << (width * height))
            if mask.bit_count() % 2
        )
        engine_ok &= all(
            sum(parity[mask][degree] for mask in range(1 << (width * height)))
            == comb(n_bonds, degree)
            for degree in range(min(n_bonds, V_ORDER) + 1)
        )
    _record(
        checks,
        "boundary polynomial sanity",
        engine_ok,
        "[LEMMA] odd boundaries vanish and every retained degree sums to binomial(n_bonds, degree)",
    )

    # -- orbit transfer equals plain transfer ---------------------------
    transfer_ok = True
    for width, height in ((1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (2, 2), (3, 2), (4, 2)):
        parity = boundary_polynomials(width, height, V_ORDER)
        orbits = Orbits(width, height)
        for layers in range(2, EXTENT_LAYERS + 1):
            for weight in (4, MAX_WEIGHT):
                transfer_ok &= slab_columns_plain(
                    parity, layers, V_ORDER, weight
                ) == slab_columns(parity, layers, V_ORDER, weight, orbits)
    _record(
        checks,
        "orbit-accelerated transfer equals plain transfer",
        transfer_ok,
        "[COMPUTATION] exact equality on all rectangles of area <= 8, heights 2..6, weights 4 and 8",
    )

    # -- primary cumulant FLM ------------------------------------------
    bulk, rectangle_weights, used_shapes = cumulant_flm(V_ORDER)
    c2_q = _sum_series(bulk, 2, (2,), V_ORDER)
    c4_q = _sum_series(bulk, 4, (2, 3), V_ORDER)
    c6_q = _sum_series(bulk, 6, (2, 3, 4), V_ORDER)
    c6_heights = {layers: bulk[f"c6_h{layers}"] for layers in (2, 3, 4)}
    c8_heights = {layers: bulk[f"c8_h{layers}"] for layers in (2, 3, 4, 5)}
    c8_q = _sum_series(bulk, 8, (2, 3, 4, 5), V_ORDER)
    odd_q = {
        order_n: _sum_series(bulk, order_n, range(1, EXTENT_LAYERS + 1), V_ORDER)
        for order_n in (1, 3, 5, 7)
    }

    stored = json.loads(C6_SERIES.read_text(encoding="utf-8"))["data"]
    lower = stored["reproduced_lower_orders"]
    stored_c6 = stored["c6_direct_coupling_q"]
    _record(
        checks,
        "lower orders reproduced through v^12",
        c2_q == _fractions(lower["c2_total_q_or_w"])
        and c4_q == _fractions(lower["c4_direct_q"])
        and c6_q == _fractions(stored_c6["full_infinite_stack"])
        and all(
            c6_heights[layers] == _fractions(
                stored_c6[
                    {
                        2: "same_gap_height2",
                        3: "two_adjacent_gaps_height3",
                        4: "three_adjacent_gaps_height4",
                    }[layers]
                ]
            )
            for layers in (2, 3, 4)
        ),
        "[COMPUTATION] the eighth-order engine reproduces every stored c2/c4/c6 coefficient and height split",
    )

    extent_ok = (
        all(not any(bulk[f"c2_h{layers}"]) for layers in (3, 4, 5, 6))
        and all(not any(bulk[f"c4_h{layers}"]) for layers in (4, 5, 6))
        and all(not any(bulk[f"c6_h{layers}"]) for layers in (5, 6))
        and all(not any(bulk[f"c8_h{layers}"]) for layers in (6,))
    )
    _record(
        checks,
        "vertical extent control",
        extent_ok,
        "[THEOREM] By the all-order vertical-extent theorem, cumulant heights above m+1 vanish "
        "identically through v^12 (c8: height 6)",
    )

    odd_vanishes = (
        all(not any(values) for values in odd_q.values())
        and all(c8_q[degree] == 0 for degree in range(1, V_ORDER + 1, 2))
        and c8_q[0] == Fraction(-17, 2520)
    )
    _record(
        checks,
        "odd interlayer orders, c8 parity, and normalization",
        odd_vanishes,
        "[COMPUTATION] primary Bell cumulants c1,c3,c5,c7 vanish through v^12; odd v-degrees of c8 vanish and "
        "[v^0] c8_q = -17/2520 = [q^8] log cosh q",
    )

    # -- exact q -> w conversion ----------------------------------------
    conversion = _conversion_coefficients()
    c8_w_total = tuple(
        c8_q[degree]
        + conversion[6] * c6_q[degree]
        + conversion[4] * c4_q[degree]
        + conversion[2] * c2_q[degree]
        for degree in range(V_ORDER + 1)
    )
    c8_w_residual = list(c8_w_total)
    c8_w_residual[0] -= Fraction(1, 8)
    c2_w_residual = _subtract_constant(c2_q, Fraction(1, 2))
    c4_w_total = tuple(c4_q[d] + Fraction(2, 3) * c2_q[d] for d in range(V_ORDER + 1))
    c4_w_residual = _subtract_constant(c4_w_total, Fraction(1, 4))
    c6_w_total = tuple(
        c6_q[d] + Fraction(4, 3) * c4_q[d] + Fraction(23, 45) * c2_q[d]
        for d in range(V_ORDER + 1)
    )
    c6_w_residual = _subtract_constant(c6_w_total, Fraction(1, 6))
    conversion_ok = (
        conversion == {2: Fraction(44, 105), 4: Fraction(22, 15), 6: Fraction(2), 8: Fraction(1)}
        and c8_w_total[0] == Fraction(1, 8)
        and c8_w_residual[2] == 2
    )
    _record(
        checks,
        "exact q-to-w conversion and ladder identity",
        conversion_ok,
        "[COMPUTATION] [w^8]atanh^2=44/105, atanh^4=22/15, atanh^6=2, atanh^8=1 from the exp series; "
        "[v^0] c8_w = 1/8 and the emergent [v^2 w^8] = 2",
    )

    per_height_w = {}
    for layers in (2, 3, 4, 5):
        values = [
            bulk[f"c8_h{layers}"][degree]
            + conversion[6] * bulk[f"c6_h{layers}"][degree]
            + conversion[4] * bulk[f"c4_h{layers}"][degree]
            + conversion[2] * bulk[f"c2_h{layers}"][degree]
            for degree in range(V_ORDER + 1)
        ]
        if layers == 2:
            values[0] -= Fraction(1, 8)
        per_height_w[layers] = tuple(values)
    height_sum_ok = all(
        sum(per_height_w[layers][degree] for layers in (2, 3, 4, 5)) == c8_w_residual[degree]
        for degree in range(V_ORDER + 1)
    )
    _record(
        checks,
        "per-height w decomposition is additive",
        height_sum_ok,
        "[COMPUTATION] converted height components sum coefficientwise to the residual",
    )

    # -- second route: w-log three-dimensional inversion -----------------
    log_bulk, log_by_layers, log_weights = wlog_bulk(V_ORDER)
    second_ok = (
        log_bulk[2] == c2_w_residual
        and log_bulk[4] == c4_w_residual
        and log_bulk[6] == c6_w_residual
        and log_bulk[8] == tuple(c8_w_residual)
    )
    height_route_ok = all(
        log_by_layers[layers][8] == per_height_w[layers] for layers in (2, 3, 4, 5)
    )
    height6_zero = all(
        not any(log_by_layers[layers][degree])
        for layers in (1, 6)
        for degree in (2, 4, 6, 8)
    )
    _record(
        checks,
        "second route: w-log box inversion",
        second_ok and height_route_ok and height6_zero,
        "[COMPUTATION] independent [w^d] log + 3D Moebius route matches c2/c4/c6/c8 residuals and every height "
        "component through v^12; heights 1 and 6 carry no weight",
    )

    # -- third route: spin DOS through v^6 -------------------------------
    dos_bulk, dos_by_layers, dos_odd_zero, dos_weights = spin_dos_bulk(
        SPIN_DOS_ORDER, range(1, MAX_LAYERS + 1)
    )
    dos_ok = dos_odd_zero and all(
        dos_bulk[degree] == tuple(c8_w_residual)[: SPIN_DOS_ORDER + 1]
        if degree == 8
        else dos_bulk[degree] == _lower_residual(degree, c2_q, c4_q, c6_q, conversion)[
            : SPIN_DOS_ORDER + 1
        ]
        for degree in (2, 4, 6, 8)
    )
    dos_heights_ok = all(
        dos_by_layers[layers][8] == per_height_w[layers][: SPIN_DOS_ORDER + 1]
        for layers in (2, 3, 4, 5)
    )
    extent_shapes = ((1, 1), (2, 1), (3, 1), (4, 1), (2, 2))
    _extent_bulk, _extent_by_layers, extent_odd, extent_weights = spin_dos_bulk(
        SPIN_DOS_ORDER,
        range(1, EXTENT_LAYERS + 1),
        shapes=extent_shapes,
    )
    extent_boxes = tuple((width, height, EXTENT_LAYERS) for width, height in extent_shapes)
    extent_per_box = {
        f"{width}x{height}x{layers}": {
            str(degree): not any(extent_weights[(width, height, layers)][degree])
            for degree in (2, 4, 6, 8)
        }
        for width, height, layers in extent_boxes
    }
    extent_ok_dos = extent_odd and all(
        all(columns_zero.values()) for columns_zero in extent_per_box.values()
    )
    _record(
        checks,
        "third route: anisotropic spin DOS through v^6",
        dos_ok and dos_heights_ok and extent_ok_dos,
        "[COMPUTATION] full spin enumeration of open boxes reproduces the residual and its height split "
        "through v^6; every odd column and every height-6 box in the listed five-shape control family vanishes",
    )

    extended_shapes = tuple(
        sorted(set(used_shapes).union(STABILITY_SHAPES), key=lambda shape: (sum(shape), shape))
    )
    _extended_bulk, extended_weights, _ = cumulant_flm(V_ORDER, shapes=extended_shapes)
    stability_ok = all(
        not any(series)
        for shape in STABILITY_SHAPES
        for series in extended_weights[shape].values()
    )
    _record(
        checks,
        "partial bound stability",
        stability_ok,
        "[COMPUTATION] span-7 rectangles (8x1) and (7x2) contribute nothing at any height through v^12",
    )

    # -- Bell(8) connected inventory -------------------------------------
    inventory, inventory_checks_ok, inventory_details = derive_inventory(checks)
    _record(
        checks,
        "eighth-order connected inventory",
        inventory_checks_ok,
        "[COMPUTATION] " + inventory_details,
    )

    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError("one or more e126 checks failed")
    elapsed = time.monotonic() - started
    if elapsed > WALL_BUDGET_SECONDS:
        raise TimeoutError(
            f"e126 exceeded its declared {WALL_BUDGET_SECONDS}s wall budget: {elapsed:.3f}s"
        )
    rss_raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


    # -- artifact ---------------------------------------------------------
    provenance_rows = []
    for degree in range(V_ORDER + 1):
        row = {
            "v_degree": degree,
            "c8_direct_q": str(c8_q[degree]),
            "c8_total_w": str(c8_w_total[degree]),
            "c8_residual_w": str(c8_w_residual[degree]),
            "spin_dos_residual_w": (
                str(dos_bulk[8][degree]) if degree <= SPIN_DOS_ORDER else None
            ),
            "two_route_agreement_available": degree <= SPIN_DOS_ORDER,
        }
        provenance_rows.append(row)

    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "method": (
                "exact open-rectangle boundary-mask layer transfer with rectangle-symmetry orbit "
                "acceleration, complete Bell(8)=4140 set-partition cumulants, vertical Moebius "
                "inversion through six layers, planar rectangular FLM; second route [w^d] log "
                "columns with 3D box Moebius inversion; third route exact anisotropic spin-DOS FLM "
                "through v^6"
            ),
            "resources": {
                "tag": "[COMPUTATION]",
                "wall_budget_seconds": WALL_BUDGET_SECONDS,
                "observed_wall_seconds": round(elapsed, 3),
                "observed_wall_is_measured": True,
                "observed_ru_maxrss_raw": rss_raw,
                "ru_maxrss_units": "bytes on Darwin; KiB on Linux",
                "peak_rss_budget_bytes": 8 * 1024**3,
                "peak_rss_budget_is_preflight": True,
            },
        },
        "data": {
            "claim": (
                "[COMPUTATION] Exact finite high-temperature c8 prefix through v^12 by the generated-"
                "cumulant route, cross-checked by two further exact routes through v^12 (w-log box "
                "inversion) and v^6 (spin DOS). No all-orders coefficient claim is made."
            ),
            "variables": {
                "v": "tanh(K), the in-plane high-temperature variable",
                "q": "K_z, direct interlayer coupling of the cumulant expansion",
                "w": "tanh(K_z)",
            },
            "v_order": V_ORDER,
            "independent_check_v_order": SPIN_DOS_ORDER,
            "second_route_v_order": V_ORDER,
            "reproduced_lower_orders": {
                "c2_total_q_or_w": _strings(c2_q),
                "c4_direct_q": _strings(c4_q),
                "c6_direct_q": _strings(c6_q),
            },
            "c8_direct_coupling_q": {
                "tag": "[COMPUTATION]",
                "full_infinite_stack": _strings(c8_q),
                "by_exact_vertical_extent": {
                    str(layers): _strings(c8_heights[layers]) for layers in (2, 3, 4, 5)
                },
                "extent_note": (
                    "height 2 is the single gap pattern (4); height 5 is the single pattern "
                    "(1,1,1,1); heights 3 and 4 carry three patterns each"
                ),
                "height6_identically_zero": _strings(bulk["c8_h6"]),
            },
            "q_to_w_conversion": {
                "tag": "[LEMMA]",
                "atanh_w8_powers": {
                    str(power): str(coefficient) for power, coefficient in conversion.items()
                },
                "formula": "c8_w_total = c8_q + 2*c6_q + (22/15)*c4_q + (44/105)*c2_q",
                "derivation": "powers of the exact atanh series from the formal identity d(atanh w)/dw = 1/(1-w^2)",
            },
            "c8_wave_variable_w_tanh_Kz": {
                "tag": "[COMPUTATION]",
                "total": _strings(c8_w_total),
                "residual_after_log_cosh": _strings(c8_w_residual),
                "residual_by_exact_vertical_extent": {
                    str(layers): _strings(per_height_w[layers]) for layers in (2, 3, 4, 5)
                },
                "nonzero_residual_coefficients": {
                    str(degree): str(c8_w_residual[degree])
                    for degree in range(2, V_ORDER + 1, 2)
                },
                "ladder_identity_v2": str(c8_w_residual[2]),
            },
            "odd_interlayer_orders": {
                "tag": "[COMPUTATION]",
                "series_q": {str(order_n): _strings(values) for order_n, values in odd_q.items()},
                "reason": (
                    "the open stack is invariant under q -> -q (flip alternate layers), and the "
                    "primary complete Bell-cumulant computation independently returns zero"
                ),
            },
            "second_route_w_log": {
                "tag": "[COMPUTATION]",
                "method": (
                    "exact [w^d] log columns of the slab even-subgraph polynomial via the Newton "
                    "recursion, then three-dimensional box Moebius inversion over 28 rectangles x 6 heights"
                ),
                "agreement_through_v12": True,
                "residual_w8": _strings(log_bulk[8]),
                "height6_and_height1_weights_zero": True,
            },
            "independent_spin_dos_route": {
                "tag": "[COMPUTATION]",
                "method": (
                    "exact joint spin density of states of open 3D boxes transformed to P(v,w), "
                    "[w^d] log extraction, 3D Moebius inversion through v^6"
                ),
                "v_order": SPIN_DOS_ORDER,
                "residual_w8": _strings(dos_bulk[8]),
                "residual_w8_by_height": {
                    str(layers): _strings(dos_by_layers[layers][8]) for layers in (2, 3, 4, 5)
                },
                "box_count": len(dos_weights),
                "odd_columns_zero": bool(dos_odd_zero),
                "height6_control": {
                    "scope": "five-shape control family only; not a claim about every open box shape",
                    "shapes": [list(shape) for shape in extent_boxes],
                    "per_box_w_columns_zero": extent_per_box,
                    "all_per_box_w_columns_zero": bool(extent_ok_dos),
                },
            },
            "bell_eighth_order_inventory": inventory,
            "finite_lattice_certificate": {
                "tag": "[COMPUTATION]",
                "rectangles": [list(shape) for shape in used_shapes],
                "rectangle_count": len(used_shapes),
                "span_bound": (
                    "a connected even subgraph spanning an a x b rectangle uses at least "
                    "2((a-1)+(b-1)) in-plane edges, so span 6 exhausts v^12"
                ),
                "stability_control_shapes": [list(shape) for shape in STABILITY_SHAPES],
                "rectangle_mobius_weights_c8_h5": {
                    f"{width}x{height}": _strings(values["c8_h5"])
                    for (width, height), values in sorted(rectangle_weights.items())
                },
            },
            "order_by_order_provenance": provenance_rows,
            "scope": (
                "[UNRESOLVED] These exact coefficients are a finite expansion around K=K_z=0. They do "
                "not determine a convergence radius, justify evaluation at the benchmark critical "
                "coupling, or solve the 3D Ising model. The connected-moment inventory is a "
                "coefficientwise formal-series identity, not an evaluated closed form of the "
                "infinite 2D lattice sums."
            ),
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"PASS e126_interlayer_c8 ({elapsed:.1f}s)")


def _lower_residual(degree, c2_q, c4_q, c6_q, conversion):
    if degree == 2:
        return _subtract_constant(c2_q, Fraction(1, 2))
    if degree == 4:
        total = tuple(c4_q[d] + Fraction(2, 3) * c2_q[d] for d in range(V_ORDER + 1))
        return _subtract_constant(total, Fraction(1, 4))
    total = tuple(
        c6_q[d] + Fraction(4, 3) * c4_q[d] + Fraction(23, 45) * c2_q[d]
        for d in range(V_ORDER + 1)
    )
    return _subtract_constant(total, Fraction(1, 6))


# ---------------------------------------------------------------------------
# eighth-order connected inventory derivation
# ---------------------------------------------------------------------------
def derive_inventory(checks: list[dict[str, object]]) -> tuple[dict[str, object], bool, str]:
    n = 8
    patterns = compositions(4)
    placement: dict[tuple[int, ...], Fraction] = {}
    assignment_counts: dict[tuple[int, ...], int] = {}
    placement_ok = True
    for composition in patterns:
        weight, count = placement_from_ordered_assignments(composition)
        placement[composition] = weight
        assignment_counts[composition] = count
        placement_ok &= weight == Fraction(1, prod_factorials(composition))
    collected: dict[str, Polynomial] = {}
    rows: list[dict[str, object]] = []
    structural_ok = True
    substitution_ok = True
    reflection_ok = True
    total_substituted = 0
    total_admissible = 0

    for composition in patterns:
        gaps = composition_gap_map(composition)
        polynomial, info = profile_connected_polynomial(gaps, placement[composition])
        cross = substitute_raw(info["raw_monomials"])
        substitution_ok &= cross == polynomial
        even_ok = all(monomial_even_incidence(m, n) for m in polynomial)
        anchor_ok = all(monomial_anchor_connected(m, n) for m in polynomial)
        structural_ok &= even_ok and anchor_ok
        histogram = Counter(atom_kinds(m) for m in polynomial)
        total_substituted += len(polynomial)
        total_admissible += int(info["admissible_partitions"])
        collected[str(composition)] = polynomial
        rows.append(
            {
                "pattern": list(composition),
                "gaps": str(dict(sorted(gaps.items()))),
                "placement_weight": str(placement[composition]),
                "ordered_gap_assignments": assignment_counts[composition],
                "admissible_partitions": int(info["admissible_partitions"]),
                "raw_layer_moment_monomials": len(info["raw_monomials"]),
                "connected_monomials": len(polynomial),
                "connected_polynomial_sha256": _polynomial_digest({str(composition): polynomial}),
                "atom_kind_histogram": {
                    ",".join(kinds): count for kinds, count in sorted(histogram.items())
                },
                "moebius_by_block_class": {
                    "+".join(map(str, cls)): str(value)
                    for cls, value in sorted(info["moebius_by_class"].items())
                },
                "all_even_incidence": bool(even_ok),
                "all_anchor_connected": bool(anchor_ok),
            }
        )

    # reflection pairing: reversing the composition relabels slots.
    reversal = {slot: n + 1 - slot for slot in range(1, n + 1)}
    for composition in patterns:
        mirrored = tuple(reversed(composition))
        if mirrored == composition:
            continue
        mapped: Polynomial = {}
        for monomial, coefficient in collected[str(composition)].items():
            key = _monomial(
                *(_atom(kind, (reversal[slot] for slot in slots)) for kind, slots in monomial)
            )
            mapped[key] = mapped.get(key, Fraction(0)) + coefficient
        reflection_ok &= mapped == collected[str(mirrored)]

    # Selection rule: literal formal Bell/moment evaluation of every normalized
    # gap multiset.  No contiguity or multiplicity predicate is applied here.
    every_multiset = normalized_gap_multisets(n)
    literal_nonzero = {
        gap_values
        for gap_values in every_multiset
        if literal_profile_moment_polynomial(gap_values)
    }
    expected = {
        tuple(sorted(composition_gap_map(composition).values())) for composition in patterns
    }
    selection_ok = literal_nonzero == expected and len(literal_nonzero) == len(patterns)
    vanish_ok = selection_ok and len(every_multiset) == 3432

    # block parity lemma: per-gap even <=> per-layer even, over all partitions
    parity_ok = True
    for partition in set_partitions(tuple(range(1, n + 1))):
        for composition in patterns:
            gaps = composition_gap_map(composition)
            for block in partition:
                parity_ok &= block_gap_even(gaps, block) == block_layer_even(gaps, block)

    # regression: sixth-order canonical form reproduces the certified theorem
    regression_ok = True
    if C6_CLOSED.exists():
        stored_formula = json.loads(C6_CLOSED.read_text(encoding="utf-8"))["data"][
            "connected_correlation_identity"
        ]["formula"]
        stored_rows = {
            row["latex"]: Fraction(row["coefficient"]) for row in stored_formula
        }
        sixth = (
            ((3,), Fraction(1, factorial(6)), Fraction(1)),
            ((2, 1), Fraction(1, factorial(4) * factorial(2)), Fraction(2)),
            ((1, 1, 1), Fraction(1, factorial(2) ** 3), Fraction(1)),
        )
        rebuilt: dict[str, Fraction] = {}
        for composition, weight, scale in sixth:
            polynomial, _info = profile_connected_polynomial(composition_gap_map(composition), weight)
            for monomial, coefficient in canonical_form(polynomial, 6).items():
                latex = " \\cdot ".join(
                    f"{kind}({','.join(map(str, slots))})" for kind, slots in monomial
                )
                rebuilt[latex] = rebuilt.get(latex, Fraction(0)) + coefficient * scale
        regression_ok &= rebuilt == stored_rows
    else:
        regression_ok = False

    # second-order anchor: c2 = (1/2) sum' G(1,2)^2
    second_gaps = composition_gap_map((1,))
    second_poly, _ = profile_connected_polynomial(second_gaps, Fraction(1, 2))
    regression_ok &= second_poly == {_monomial(_atom("G", (1, 2)), _atom("G", (1, 2))): Fraction(1, 2)}

    # sampled numeric verification of raw and collected evaluations
    numeric_ok = True
    for width, height in ((2, 2), (1, 3)):
        moments = layer_moment_values(width, height, IDENTITY_ORDER)
        tuples = _deterministic_tuples(width * height, 4)
        for composition in patterns:
            gaps = composition_gap_map(composition)
            for positions in tuples:
                recursion = joint_cumulant_by_recursion(gaps, positions, moments, IDENTITY_ORDER)
                literal = evaluate_raw_partitions(
                    gaps, placement[composition], positions, moments, IDENTITY_ORDER
                )
                numeric_ok &= literal == tuple(
                    placement[composition] * value for value in recursion
                )
        tuples_small = _deterministic_tuples(width * height, 2)
        for composition in patterns:
            gaps = composition_gap_map(composition)
            for positions in tuples_small:
                recursion = joint_cumulant_by_recursion(gaps, positions, moments, IDENTITY_ORDER)
                collected_value = evaluate_collected(
                    collected[str(composition)], gaps, positions, moments, IDENTITY_ORDER
                )
                numeric_ok &= collected_value == tuple(
                    placement[composition] * value for value in recursion
                )

    inventory = {
        "tag": "[COMPUTATION]",
        "order": n,
        "bell_number_B8": bell_number(8),
        "patterns": rows,
        "pattern_count": len(patterns),
        "total_admissible_partitions": total_admissible,
        "total_connected_monomials": total_substituted,
        "atom_kinds_present": sorted(
            {kind for poly in collected.values() for monomial in poly for kind, _ in monomial}
        ),
        "collected_digest_sha256": _polynomial_digest(collected),
        "placement_factors_independently_derived": bool(placement_ok),
        "normalized_gap_multisets_tested": len(every_multiset),
        "literal_nonzero_gap_multisets": [list(values) for values in sorted(literal_nonzero)],
        "selection_rule_verified": bool(selection_ok),
        "inadmissible_patterns_vanish": bool(vanish_ok),
        "block_parity_equivalence": bool(parity_ok),
        "reflection_pairing": bool(reflection_ok),
        "substitution_two_paths_equal": bool(substitution_ok),
        "sixth_order_regression": bool(regression_ok),
        "sampled_numeric_verification": bool(numeric_ok),
        "counts_note": (
            "connected_monomials counts distinct slot-indexed monomials before S_7 orbit "
            "collapsing; the same pipeline collapses to the certified 13-row sixth-order form"
        ),
    }
    details = (
        f"Bell(8)={bell_number(8)}; {len(patterns)} gap patterns; "
        f"{len(every_multiset)} literal normalized gap multisets; "
        f"{total_substituted} distinct connected monomials over G/U/W/W8; "
        f"assignment-derived placement factors, exhaustive literal selection, reflection "
        f"pairing, block parity, two-path substitution, sixth-order regression and sampled "
        f"numeric identity all verified"
    )
    ok = (
        structural_ok
        and substitution_ok
        and reflection_ok
        and placement_ok
        and selection_ok
        and vanish_ok
        and parity_ok
        and regression_ok
        and numeric_ok
    )
    return inventory, ok, details


def prod_factorials(composition: Sequence[int]) -> int:
    out = 1
    for part in composition:
        out *= factorial(2 * int(part))
    return out


def bell_number(n: int) -> int:
    return len(set_partitions(tuple(range(n))))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:  # pragma: no cover
        print(f"FAIL e126_interlayer_c8: {error}")
        raise
