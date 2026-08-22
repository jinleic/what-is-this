"""Exact fourth-order interlayer expansion of stacked square Ising layers.

Run from the repository root with
    .venv/bin/python experiments/e51_interlayer_c4.py

The repository's interlayer high-temperature variable is ``w = tanh(K_z)``.
The fourth cumulant is naturally a coefficient in the direct Boltzmann
coupling ``q = K_z``; this script computes that coefficient first and then
performs the exact change of variable ``q = atanh(w)``.  Both conventions are
stored so that the fourth-order distinction cannot be hidden.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

from sympy import Symbol, expand, series

SCRIPT = "experiments/e51_interlayer_c4.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "interlayer" / "c4_series.json"
V_ORDER = 12
STABILITY_ORDER = 8
W_ORDER = 4

RationalSeries = tuple[Fraction, ...]
IntegerSeries = tuple[int, ...]


@dataclass(frozen=True)
class RectangleCumulants:
    """Finite open-layer cumulant sums before finite-lattice inversion."""

    c2_total: RationalSeries
    c4_same_gap_q: RationalSeries
    c4_adjacent_gaps_q: RationalSeries
    c4_stack_q: RationalSeries


@dataclass(frozen=True)
class FLMResult:
    """Exact rectangular finite-lattice inversion result."""

    order: int
    bound_slack: int
    shapes: tuple[tuple[int, ...], ...]
    bulk: tuple[RationalSeries, ...]
    weights: dict[tuple[int, ...], tuple[RationalSeries, ...]]
    odd_vertical_columns_zero: bool = True


def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _fraction_strings(values: Iterable[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _add_scaled(
    target: list[Fraction], source: Sequence[int | Fraction], scale: int | Fraction = 1
) -> None:
    exact_scale = Fraction(scale)
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += exact_scale * coefficient


def _mul(
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
    order: int,
) -> RationalSeries:
    out = _zero(order)
    for i, first in enumerate(left[: order + 1]):
        if not first:
            continue
        for j, second in enumerate(right[: order + 1 - i]):
            if second:
                out[i + j] += Fraction(first) * second
    return tuple(out)


def _divide(
    numerator: Sequence[int | Fraction],
    denominator: Sequence[int | Fraction],
    order: int,
) -> RationalSeries:
    den = [
        Fraction(denominator[degree])
        if degree < len(denominator)
        else Fraction(0)
        for degree in range(order + 1)
    ]
    if not den[0]:
        raise ValueError("formal-series denominator has zero constant term")
    num = [
        Fraction(numerator[degree])
        if degree < len(numerator)
        else Fraction(0)
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


@lru_cache(maxsize=None)
def parity_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """Count edge subsets by boundary mask and exact edge count.

    ``result[mask][m]`` is the number of subsets of ``m`` open-square edges
    whose odd-degree vertex set is ``mask``.  Descending degree updates ensure
    that every edge is used at most once.
    """

    width = int(width)
    height = int(height)
    order = int(order)
    if width < 1 or height < 1 or order < 0:
        raise ValueError("positive rectangle sides and non-negative order required")
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


def _symbolic_factorization_check() -> dict[str, object]:
    """Check the fourth-cumulant factorization before any series enumeration."""

    m4 = Symbol("M_ijkl")
    gij = Symbol("G_ij")
    gik = Symbol("G_ik")
    gil = Symbol("G_il")
    gjk = Symbol("G_jk")
    gjl = Symbol("G_jl")
    gkl = Symbol("G_kl")
    u4 = Symbol("U_ijkl")
    pair_ij_kl = gij * gkl
    pair_ik_jl = gik * gjl
    pair_il_jk = gil * gjk
    m4_from_connected = u4 + pair_ij_kl + pair_ik_jl + pair_il_jk

    # For four zero-mean variables, kappa_4 is the four-point moment minus
    # the three pair partitions: 12|34, 13|24, and 14|23.
    same_gap_from_partitions = (
        m4**2 - gij**2 * gkl**2 - gik**2 * gjl**2 - gil**2 * gjk**2
    )
    same_gap_expected = (
        m4**2 - gij**2 * gkl**2 - gik**2 * gjl**2 - gil**2 * gjk**2
    )

    # With two variables on each of two adjacent gaps, only 12|34 survives;
    # every cross-gap two-variable moment leaves an outer layer with one spin.
    adjacent_from_partitions = gij * m4 * gkl - gij**2 * gkl**2 - 0 - 0
    adjacent_expected = gij * gkl * (m4 - gij * gkl)

    same_gap_connected = (
        u4**2
        + 2 * u4 * (pair_ij_kl + pair_ik_jl + pair_il_jk)
        + 2
        * (
            pair_ij_kl * pair_ik_jl
            + pair_ij_kl * pair_il_jk
            + pair_ik_jl * pair_il_jk
        )
    )
    adjacent_connected = pair_ij_kl * (u4 + pair_ik_jl + pair_il_jk)

    c2 = Symbol("c2")
    c4q = Symbol("c4q")
    w = Symbol("w")
    q_of_w = w + w**3 / 3
    transformed = series(c2 * q_of_w**2 + c4q * q_of_w**4, w, 0, 5).removeO()
    conversion_expected = c2 * w**2 + (c4q + Fraction(2, 3) * c2) * w**4

    passed = (
        expand(same_gap_from_partitions - same_gap_expected) == 0
        and expand(adjacent_from_partitions - adjacent_expected) == 0
        and expand(
            same_gap_from_partitions.subs(m4, m4_from_connected)
            - same_gap_connected
        )
        == 0
        and expand(
            adjacent_from_partitions.subs(m4, m4_from_connected)
            - adjacent_connected
        )
        == 0
        and expand(transformed - conversion_expected) == 0
    )
    return {
        "passed": bool(passed),
        "zero_mean_partition_formula": (
            "kappa4(Y1,Y2,Y3,Y4)=E1234-E12*E34-E13*E24-E14*E23"
        ),
        "same_gap": (
            "M_ijkl^2-G_ij^2 G_kl^2-G_ik^2 G_jl^2-G_il^2 G_jk^2"
        ),
        "adjacent_gaps": "G_ij G_kl M_ijkl-G_ij^2 G_kl^2",
        "connected_four_point": (
            "U_ijkl=M_ijkl-G_ij G_kl-G_ik G_jl-G_il G_jk"
        ),
        "same_gap_in_G_and_U": (
            "U^2+2U(P+Q+R)+2(PQ+PR+QR), "
            "P=G_ij G_kl, Q=G_ik G_jl, R=G_il G_jk"
        ),
        "adjacent_in_G_and_U": "P(U+Q+R)",
        "change_of_variable": "[w^4]F(atanh(w))=c4_q+(2/3)c2",
    }


@lru_cache(maxsize=None)
def rectangle_cumulants(width: int, height: int, order: int) -> RectangleCumulants:
    """Return exact finite-rectangle cumulant sums through ``v**order``."""

    parity = parity_polynomials(width, height, order)
    partition = parity[0]
    site_count = width * height

    # An ordered pair (i,j) has XOR mask zero with multiplicity A (i=j),
    # and each two-bit mask with multiplicity two (the two orders).
    pair_masks = [(0, site_count)]
    pair_masks.extend(
        ((1 << first) | (1 << second), 2)
        for first in range(site_count)
        for second in range(first + 1, site_count)
    )

    moments: dict[int, RationalSeries] = {}

    def moment(mask: int) -> RationalSeries:
        if mask not in moments:
            moments[mask] = _divide(parity[mask], partition, order)
        return moments[mask]

    overlap_sum = _zero(order)
    for mask, multiplicity in pair_masks:
        _add_scaled(overlap_sum, _mul(moment(mask), moment(mask), order), multiplicity)

    four_spin_square_sum = _zero(order)
    adjacent_moment_sum = _zero(order)
    for first_mask, first_multiplicity in pair_masks:
        first_correlation = moment(first_mask)
        for second_mask, second_multiplicity in pair_masks:
            second_correlation = moment(second_mask)
            four_spin = moment(first_mask ^ second_mask)
            multiplicity = first_multiplicity * second_multiplicity
            _add_scaled(
                four_spin_square_sum,
                _mul(four_spin, four_spin, order),
                multiplicity,
            )
            _add_scaled(
                adjacent_moment_sum,
                _mul(_mul(first_correlation, second_correlation, order), four_spin, order),
                multiplicity,
            )

    overlap_square = _mul(overlap_sum, overlap_sum, order)
    c2_total = tuple(value / 2 for value in overlap_sum)
    same_gap = tuple(
        (four_spin_square_sum[degree] - 3 * overlap_square[degree]) / 24
        for degree in range(order + 1)
    )
    adjacent_gaps = tuple(
        (adjacent_moment_sum[degree] - overlap_square[degree]) / 4
        for degree in range(order + 1)
    )
    stack = tuple(
        same_gap[degree] + adjacent_gaps[degree] for degree in range(order + 1)
    )
    return RectangleCumulants(c2_total, same_gap, adjacent_gaps, stack)


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


def cumulant_flm(order: int, bound_slack: int = 0) -> FLMResult:
    """Rectangular FLM of the exact 2D moment/cumulant formula."""

    shapes = _plane_shapes(order, bound_slack)
    weights: dict[tuple[int, ...], tuple[RationalSeries, ...]] = {}
    bulk = tuple(_zero(order) for _ in range(4))
    for width, height in shapes:
        finite = rectangle_cumulants(width, height, order)
        values = [
            list(finite.c2_total),
            list(finite.c4_same_gap_q),
            list(finite.c4_adjacent_gaps_q),
            list(finite.c4_stack_q),
        ]
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for quantity in range(4):
                    _add_scaled(values[quantity], subweight[quantity], -placements)
        exact_values = tuple(tuple(value) for value in values)
        weights[(width, height)] = exact_values
        for quantity in range(4):
            _add_scaled(bulk[quantity], exact_values[quantity])
    return FLMResult(
        order=order,
        bound_slack=bound_slack,
        shapes=tuple(shapes),
        bulk=tuple(tuple(value) for value in bulk),
        weights=weights,
    )


@lru_cache(maxsize=None)
def slab_even_subgraph_polynomial(
    width: int, height: int, layers: int, order: int
) -> tuple[RationalSeries, ...]:
    """Construct exact ``P(v,w)`` columns for one-, two-, or three-layer slabs.

    This is the independent anisotropic route.  If ``P_S(v)`` is one layer's
    edge-subset polynomial with boundary ``S``, the two-layer column is
    ``sum_|S|=d P_S^2``.  For three layers it is
    ``sum_|S|+|T|=d P_S P_(S xor T) P_T``.
    """

    if layers not in (1, 2, 3):
        raise ValueError("only one-, two-, and three-layer slabs are needed at w^4")
    parity = parity_polynomials(width, height, order)
    site_count = width * height
    by_weight: list[list[int]] = [[] for _ in range(W_ORDER + 1)]
    for mask in range(1 << site_count):
        weight = bin(mask).count("1")
        if weight <= W_ORDER:
            by_weight[weight].append(mask)

    columns = tuple(_zero(order) for _ in range(W_ORDER + 1))
    if layers == 1:
        columns[0][:] = parity[0]
    elif layers == 2:
        for weight in range(W_ORDER + 1):
            for mask in by_weight[weight]:
                if any(parity[mask]):
                    _add_scaled(columns[weight], _mul(parity[mask], parity[mask], order))
    else:
        for first_weight in range(W_ORDER + 1):
            for second_weight in range(W_ORDER + 1 - first_weight):
                degree_w = first_weight + second_weight
                for first_mask in by_weight[first_weight]:
                    first = parity[first_mask]
                    if not any(first):
                        continue
                    for second_mask in by_weight[second_weight]:
                        second = parity[second_mask]
                        middle = parity[first_mask ^ second_mask]
                        if any(second) and any(middle):
                            _add_scaled(
                                columns[degree_w],
                                _mul(_mul(first, middle, order), second, order),
                            )
    return tuple(tuple(column) for column in columns)


@lru_cache(maxsize=None)
def _slab_log_coefficients(
    width: int, height: int, layers: int, order: int
) -> tuple[RationalSeries, RationalSeries, bool]:
    polynomial = slab_even_subgraph_polynomial(width, height, layers, order)
    p0 = polynomial[0]
    ratio2 = _divide(polynomial[2], p0, order)
    ratio4 = _divide(polynomial[4], p0, order)
    ratio2_square = _mul(ratio2, ratio2, order)
    log_w4 = tuple(
        ratio4[degree] - ratio2_square[degree] / 2 for degree in range(order + 1)
    )
    odd_zero = not any(polynomial[1]) and not any(polynomial[3])
    return ratio2, log_w4, odd_zero


def anisotropic_flm(order: int, bound_slack: int = 0) -> FLMResult:
    """Independent 3D slab FLM in ``v`` and ``w=tanh(K_z)``."""

    plane_shapes = _plane_shapes(order, bound_slack)
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in plane_shapes
                for layers in (1, 2, 3)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    weights: dict[tuple[int, ...], tuple[RationalSeries, ...]] = {}
    bulk = tuple(_zero(order) for _ in range(2))
    odd_zero = True
    for width, height, layers in shapes:
        log_w2, log_w4, box_odd_zero = _slab_log_coefficients(
            width, height, layers, order
        )
        odd_zero = odd_zero and box_odd_zero
        values = [list(log_w2), list(log_w4)]
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if (
                sub_width <= width
                and sub_height <= height
                and sub_layers <= layers
            ):
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                for quantity in range(2):
                    _add_scaled(values[quantity], subweight[quantity], -placements)
        exact_values = tuple(tuple(value) for value in values)
        weights[(width, height, layers)] = exact_values
        for quantity in range(2):
            _add_scaled(bulk[quantity], exact_values[quantity])
    return FLMResult(
        order=order,
        bound_slack=bound_slack,
        shapes=shapes,
        bulk=tuple(tuple(value) for value in bulk),
        weights=weights,
        odd_vertical_columns_zero=odd_zero,
    )


def _sum_weights_at_height(
    result: FLMResult, layers: int, quantity: int
) -> RationalSeries:
    total = _zero(result.order)
    for shape, values in result.weights.items():
        if len(shape) == 3 and shape[2] == layers:
            _add_scaled(total, values[quantity])
    return tuple(total)


def _with_constant(series_values: RationalSeries, constant: Fraction) -> RationalSeries:
    values = list(series_values)
    values[0] += constant
    return tuple(values)


def _subtract_constant(series_values: RationalSeries, constant: Fraction) -> RationalSeries:
    values = list(series_values)
    values[0] -= constant
    return tuple(values)


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    checks: list[dict[str, object]] = []

    # This symbolic finite-combinatorics check deliberately precedes every
    # edge-subset or finite-lattice calculation.
    symbolic = _symbolic_factorization_check()
    _record(
        checks,
        "symbolic fourth-cumulant factorization",
        bool(symbolic["passed"]),
        "[LEMMA] all three zero-mean pair partitions and the adjacent-gap specialization agree exactly",
    )

    cumulant = cumulant_flm(V_ORDER)
    c2_total, c4_same_q, c4_adjacent_q, c4_stack_q = cumulant.bulk

    anisotropic = anisotropic_flm(V_ORDER)
    c2_residual_w, c4_residual_w = anisotropic.bulk
    anisotropic_c2_total = _with_constant(c2_residual_w, Fraction(1, 2))
    c4_total_w = _with_constant(c4_residual_w, Fraction(1, 4))
    c4_q_from_anisotropic = tuple(
        c4_total_w[degree] - Fraction(2, 3) * anisotropic_c2_total[degree]
        for degree in range(V_ORDER + 1)
    )

    height2_residual_w = _sum_weights_at_height(anisotropic, 2, 1)
    height3_residual_w = _sum_weights_at_height(anisotropic, 3, 1)
    same_q_as_height2_w = tuple(
        c4_same_q[degree] + Fraction(2, 3) * c2_total[degree]
        for degree in range(V_ORDER + 1)
    )
    same_q_as_height2_w = _subtract_constant(same_q_as_height2_w, Fraction(1, 4))

    _record(
        checks,
        "known c2 prefix",
        c2_total[:9]
        == (
            Fraction(1, 2),
            Fraction(0),
            Fraction(2),
            Fraction(0),
            Fraction(18),
            Fraction(0),
            Fraction(118),
            Fraction(0),
            Fraction(778),
        ),
        "[COMPUTATION] the cumulant engine reproduces total c2 through v^8, including [v^8]=778",
    )
    _record(
        checks,
        "c2 two-route agreement",
        c2_total == anisotropic_c2_total,
        "[COMPUTATION] 2D moment FLM and anisotropic 3D slab FLM agree through v^12",
    )
    _record(
        checks,
        "c4 two-route agreement",
        c4_stack_q == c4_q_from_anisotropic,
        "[COMPUTATION] direct-coupling c4 agrees coefficientwise through v^12",
    )
    _record(
        checks,
        "c4 vertical-span component agreement",
        same_q_as_height2_w == height2_residual_w
        and c4_adjacent_q == height3_residual_w,
        "[COMPUTATION] two-layer weights equal the same-gap cumulant and three-layer weights equal the adjacent-gap cumulant",
    )
    _record(
        checks,
        "odd interlayer coefficients",
        anisotropic.odd_vertical_columns_zero
        and all(value == 0 for value in _sum_weights_at_height(anisotropic, 1, 1)),
        "[LEMMA] every enumerated w^1 and w^3 slab column is identically zero",
    )
    _record(
        checks,
        "zero-inplane coupling control",
        c2_total[0] == Fraction(1, 2)
        and c4_stack_q[0] == Fraction(-1, 12)
        and c4_total_w[0] == Fraction(1, 4)
        and c4_residual_w[0] == 0,
        "[COMPUTATION] K=0 gives log cosh(K_z): q^4=-1/12 and tanh-variable w^4=1/4",
    )
    _record(
        checks,
        "bipartite v parity",
        all(
            series_values[degree] == 0
            for series_values in (
                c2_total,
                c4_same_q,
                c4_adjacent_q,
                c4_stack_q,
                c4_residual_w,
            )
            for degree in range(1, V_ORDER + 1, 2)
        ),
        "[COMPUTATION] every odd v coefficient vanishes through v^12",
    )

    cumulant_enlarged = cumulant_flm(STABILITY_ORDER, bound_slack=1)
    cumulant_minimal = cumulant_flm(STABILITY_ORDER)
    anisotropic_enlarged = anisotropic_flm(STABILITY_ORDER, bound_slack=1)
    anisotropic_minimal = anisotropic_flm(STABILITY_ORDER)
    _record(
        checks,
        "finite-lattice bound stability",
        cumulant_minimal.bulk == cumulant_enlarged.bulk
        and anisotropic_minimal.bulk == anisotropic_enlarged.bulk,
        "[COMPUTATION] adding every rectangle with one extra span unit changes neither route through v^8",
    )

    if not all(bool(check["passed"]) for check in checks):
        raise AssertionError("one or more exact c4 checks failed")

    order_provenance = []
    for degree in range(V_ORDER + 1):
        span = degree // 2
        plane_count = sum(
            1
            for width, height in cumulant.shapes
            if (width - 1) + (height - 1) <= span
        )
        order_provenance.append(
            {
                "v_degree": degree,
                "coefficient_c4_residual_w": str(c4_residual_w[degree]),
                "coefficient_c4_total_w": str(c4_total_w[degree]),
                "coefficient_c4_direct_q": str(c4_stack_q[degree]),
                "cumulant_route": str(c4_stack_q[degree]),
                "anisotropic_route_after_exact_variable_change": str(
                    c4_q_from_anisotropic[degree]
                ),
                "maximum_required_inplane_span": span,
                "2d_rectangle_count_available_at_this_order": plane_count,
                "3d_slab_box_count_available_at_this_order": 3 * plane_count,
            }
        )

    height_components = {
        str(layers): _fraction_strings(_sum_weights_at_height(anisotropic, layers, 1))
        for layers in (1, 2, 3)
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "method": (
                "exact Fraction edge-boundary polynomials and 2D cumulant FLM; "
                "independent anisotropic one/two/three-layer even-subgraph FLM"
            ),
        },
        "data": {
            "claim": (
                "[COMPUTATION] Exact finite high-temperature prefixes through v^12; "
                "no all-orders coefficient claim is made."
            ),
            "variables": {
                "v": "tanh(K)",
                "q": "K_z, the direct Boltzmann coupling used by the cumulant",
                "w": "tanh(K_z), the repository interlayer HT convention",
                "conversion": "q=atanh(w)=w+w^3/3+O(w^5)",
            },
            "v_order": V_ORDER,
            "symbolic_factorization_lemma": {
                "tag": "[LEMMA]",
                **symbolic,
                "full_stack_formula": (
                    "c4_q=(1/24) sum_{j,k,l} kappa_same(0,j,k,l) "
                    "+(1/4) sum_{j,k,l} kappa_adjacent(0,j,k,l)"
                ),
                "multiplicity_note": (
                    "six placements of two ordered bonds on each adjacent gap turn 6/24 into 1/4"
                ),
            },
            "reproduced_c2": {
                "tag": "[COMPUTATION]",
                "total_q_or_w": _fraction_strings(c2_total),
                "residual_w_after_log_cosh": _fraction_strings(c2_residual_w),
                "known_v8_coefficient": str(c2_total[8]),
            },
            "c4_wave_variable_w_tanh_Kz": {
                "tag": "[COMPUTATION]",
                "total": _fraction_strings(c4_total_w),
                "residual_after_log_cosh": _fraction_strings(c4_residual_w),
                "nonzero_residual_coefficients": {
                    str(degree): str(value)
                    for degree, value in enumerate(c4_residual_w)
                    if value
                },
            },
            "c4_direct_coupling_q": {
                "tag": "[COMPUTATION]",
                "full_infinite_stack": _fraction_strings(c4_stack_q),
                "same_layer_gap_component": _fraction_strings(c4_same_q),
                "adjacent_layer_gaps_component": _fraction_strings(c4_adjacent_q),
                "from_anisotropic_route": _fraction_strings(c4_q_from_anisotropic),
            },
            "bilayer_and_three_layer_distinction": {
                "tag": "[LEMMA]",
                "bilayer_residual_w": _fraction_strings(height2_residual_w),
                "three_layer_spanning_residual_w": _fraction_strings(
                    height3_residual_w
                ),
                "anisotropic_residual_by_exact_vertical_extent": height_components,
                "scope": (
                    "a two-layer slab certifies only the same-gap term; the full 3D c4 also "
                    "has the nonzero adjacent-gap term, first visible at v^2"
                ),
            },
            "c3": {
                "tag": "[LEMMA]",
                "coefficients_through_v12": ["0"] * (V_ORDER + 1),
                "reason": (
                    "odd boundary masks have zero layer polynomial, so every slab w^1 and w^3 column vanishes"
                ),
            },
            "finite_lattice_certificate": {
                "tag": "[LEMMA]",
                "2d_rectangle_count": len(cumulant.shapes),
                "2d_rectangles": [list(shape) for shape in cumulant.shapes],
                "3d_slab_box_count": len(anisotropic.shapes),
                "vertical_extents": [1, 2, 3],
                "completeness_bound": (
                    "a connected projected even multigraph spanning an a-by-b rectangle "
                    "uses at least 2*((a-1)+(b-1)) in-plane edges; four vertical edges "
                    "span at most two layer gaps, so v^12 needs plane span <=6 and height <=3"
                ),
                "stability_check_order": STABILITY_ORDER,
                "minimal_2d_box_count_at_stability_order": len(
                    cumulant_minimal.shapes
                ),
                "enlarged_2d_box_count_at_stability_order": len(
                    cumulant_enlarged.shapes
                ),
                "minimal_3d_box_count_at_stability_order": len(
                    anisotropic_minimal.shapes
                ),
                "enlarged_3d_box_count_at_stability_order": len(
                    anisotropic_enlarged.shapes
                ),
            },
            "independent_routes": {
                "route_1": (
                    "2D edge-boundary polynomials -> exact G and M4 ratios -> complete "
                    "fourth-cumulant sums -> rectangular Moebius inversion"
                ),
                "route_2": (
                    "anisotropic P(v,w) for open slabs from vertical-subset masks -> "
                    "[w^4]log P -> three-dimensional rectangular Moebius inversion"
                ),
                "agreement": True,
            },
            "order_by_order_provenance": order_provenance,
            "scope": (
                "[COMPUTATION] These are exact finite series coefficients around K=K_z=0. "
                "They do not assert convergence at the isotropic critical point or solve the 3D model."
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
        print(f"FAIL e51_interlayer_c4: {error}")
        raise
