"""Clean-room exact certificate for the interlayer c6/c4 closed-form result.

This verifier deliberately does not import the producer.  It independently rebuilds
both the anisotropic even-subgraph finite-lattice calculation and the c4 contraction
from integer boundary-polynomial transfer data.
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

from ising.interlayer import anisotropic_box_even_subgraph


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "closed_form.json"
PRIMARY_C6 = ROOT / "results" / "interlayer" / "c6_series.json"
PRIMARY_C4 = ROOT / "results" / "interlayer" / "c4_series.json"
ORDER = 12
MAX_LAYERS = 4

IntegerSeries = tuple[int, ...]
RationalSeries = tuple[Fraction, ...]


def _zero(order: int = ORDER) -> list[int]:
    return [0] * (order + 1)


def _fraction_zero(order: int = ORDER) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _add_scaled(target: list[int], source: Sequence[int], scale: int = 1) -> None:
    for degree, value in enumerate(source):
        if value:
            target[degree] += scale * value


def _mul(left: Sequence[int], right: Sequence[int], order: int = ORDER) -> IntegerSeries:
    out = _zero(order)
    for left_degree, left_value in enumerate(left[: order + 1]):
        if not left_value:
            continue
        for right_degree, right_value in enumerate(right[: order + 1 - left_degree]):
            if right_value:
                out[left_degree + right_degree] += left_value * right_value
    return tuple(out)


def _square(values: Sequence[int], order: int = ORDER) -> IntegerSeries:
    return _mul(values, values, order)


def _divide_unit(numerator: Sequence[int], denominator: Sequence[int], order: int = ORDER) -> IntegerSeries:
    """Divide exact formal series with unit constant term without rationals."""
    if denominator[0] != 1:
        raise AssertionError("boundary partition polynomial is not monic")
    quotient = _zero(order)
    for degree in range(order + 1):
        quotient[degree] = numerator[degree] - sum(
            denominator[lower] * quotient[degree - lower]
            for lower in range(1, degree + 1)
        )
    return tuple(quotient)


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
def _boundary_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """Independent edge-by-edge transfer indexed by the odd boundary mask."""
    sites = width * height
    values = [[0] * (order + 1) for _ in range(1 << sites)]
    values[0][0] = 1
    for left, right in _square_bonds(width, height):
        toggle = (1 << left) | (1 << right)
        for degree in range(order - 1, -1, -1):
            for mask in range(1 << sites):
                coefficient = values[mask][degree]
                if coefficient:
                    values[mask ^ toggle][degree + 1] += coefficient
    return tuple(tuple(row) for row in values)


@lru_cache(maxsize=None)
def _moment(width: int, height: int, mask: int, order: int = ORDER) -> IntegerSeries:
    boundary = _boundary_polynomials(width, height, order)
    return _divide_unit(boundary[mask], boundary[0], order)


def _masks_by_weight(width: int, height: int) -> dict[int, tuple[int, ...]]:
    sites = width * height
    return {
        weight: tuple(mask for mask in range(1 << sites) if mask.bit_count() == weight)
        for weight in (0, 2, 4, 6)
    }


def _sum_squares(width: int, height: int, masks: Iterable[int]) -> IntegerSeries:
    total = _zero()
    for mask in masks:
        _add_scaled(total, _square(_moment(width, height, mask)))
    return tuple(total)


def _three_layer_pair_sum(
    width: int,
    height: int,
    first_masks: Iterable[int],
    second_masks: Iterable[int],
) -> IntegerSeries:
    """Sum m_A m_(A xor B) m_B for the indicated vertical-mask weights."""
    total = _zero()
    first = tuple(first_masks)
    second = tuple(second_masks)
    for first_mask in first:
        first_moment = _moment(width, height, first_mask)
        for second_mask in second:
            middle = _moment(width, height, first_mask ^ second_mask)
            term = _mul(_mul(first_moment, middle), _moment(width, height, second_mask))
            _add_scaled(total, term)
    return tuple(total)


def _three_gap_two_two_two_sum(width: int, height: int, masks_two: Iterable[int]) -> IntegerSeries:
    """Factor the three-gap w^6 sum through its middle boundary mask.

    For A,B,C of weight two,
      sum m_A m_(A xor B) m_(B xor C) m_C
        = sum_B (sum_A m_A m_(A xor B))^2.
    """
    total = _zero()
    masks = tuple(masks_two)
    for middle_mask in masks:
        half = _zero()
        for outside_mask in masks:
            _add_scaled(
                half,
                _mul(
                    _moment(width, height, outside_mask),
                    _moment(width, height, outside_mask ^ middle_mask),
                ),
            )
        _add_scaled(total, _square(tuple(half)))
    return tuple(total)


def _normalized_slab_columns(width: int, height: int, layers: int) -> tuple[IntegerSeries, ...]:
    """Build P(v,w)/P(v,0) through w^6 for one open slab.

    This is a closed-form boundary-mask transfer rather than a call into the
    producer.  Parity makes every odd w column identically zero.
    """
    one = (1,) + (0,) * ORDER
    zero = (0,) * (ORDER + 1)
    if layers == 1:
        return one, zero, zero, zero

    masks = _masks_by_weight(width, height)
    square_two = _sum_squares(width, height, masks[2])
    square_four = _sum_squares(width, height, masks[4])
    square_six = _sum_squares(width, height, masks[6])
    if layers == 2:
        return one, square_two, square_four, square_six

    pair_two_two = _three_layer_pair_sum(width, height, masks[2], masks[2])
    pair_two_four = _three_layer_pair_sum(width, height, masks[2], masks[4])
    if layers == 3:
        return (
            one,
            tuple(2 * value for value in square_two),
            tuple(2 * square_four[degree] + pair_two_two[degree] for degree in range(ORDER + 1)),
            tuple(2 * square_six[degree] + 2 * pair_two_four[degree] for degree in range(ORDER + 1)),
        )

    if layers != 4:
        raise AssertionError("the w^6 certificate uses at most four layers")
    triple_two = _three_gap_two_two_two_sum(width, height, masks[2])
    square_two_square = _square(square_two)
    two_four_product = _mul(square_two, square_four)
    return (
        one,
        tuple(3 * value for value in square_two),
        tuple(
            3 * square_four[degree] + 2 * pair_two_two[degree] + square_two_square[degree]
            for degree in range(ORDER + 1)
        ),
        tuple(
            3 * square_six[degree]
            + 4 * pair_two_four[degree]
            + 2 * two_four_product[degree]
            + triple_two[degree]
            for degree in range(ORDER + 1)
        ),
    )


def _log_w6(columns: tuple[IntegerSeries, ...]) -> RationalSeries:
    _, r2, r4, r6 = columns
    r2r4 = _mul(r2, r4)
    r2cubed = _mul(_square(r2), r2)
    return tuple(
        Fraction(r6[degree] - r2r4[degree]) + Fraction(r2cubed[degree], 3)
        for degree in range(ORDER + 1)
    )


def _box_shapes(order: int = ORDER) -> tuple[tuple[int, int, int], ...]:
    span = order // 2
    return tuple(
        sorted(
            (
                (width, height, layers)
                for width in range(1, span + 2)
                for height in range(1, span + 2)
                if (width - 1) + (height - 1) <= span
                for layers in range(1, MAX_LAYERS + 1)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )


def _independent_c6_flm() -> tuple[RationalSeries, dict[tuple[int, int, int], RationalSeries]]:
    """Clean-room 3D rectangular finite-lattice inversion of [w^6] log P."""
    weights: dict[tuple[int, int, int], RationalSeries] = {}
    bulk = _fraction_zero()
    for width, height, layers in _box_shapes():
        if layers == 1:
            local = _fraction_zero()
        else:
            local = list(_log_w6(_normalized_slab_columns(width, height, layers)))
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                for degree, coefficient in enumerate(subweight):
                    local[degree] -= placements * coefficient
        exact = tuple(local)
        weights[(width, height, layers)] = exact
        for degree, coefficient in enumerate(exact):
            bulk[degree] += coefficient
    return tuple(bulk), weights


def _sum_height(
    weights: dict[tuple[int, int, int], RationalSeries], layers: int
) -> RationalSeries:
    return tuple(
        sum(
            (values[degree] for shape, values in weights.items() if shape[2] == layers),
            Fraction(0),
        )
        for degree in range(ORDER + 1)
    )


def _series_subtract(*terms: tuple[int, Sequence[int]]) -> IntegerSeries:
    result = _zero()
    for scale, values in terms:
        _add_scaled(result, values, scale)
    return tuple(result)


def _closed_c4_rectangle(width: int, height: int) -> tuple[IntegerSeries, IntegerSeries, IntegerSeries, IntegerSeries, IntegerSeries]:
    """Evaluate c2 and the three G/U contractions on one exact open rectangle."""
    sites = width * height
    pairs = tuple((i, j, (1 << i) ^ (1 << j)) for i, j in product(range(sites), repeat=2))
    overlap = _zero()
    u_square = _zero()
    g_u = _zero()
    g_cross = _zero()
    direct = _zero()
    for _, _, mask in pairs:
        _add_scaled(overlap, _square(_moment(width, height, mask)))
    for i, j, ij in pairs:
        gij = _moment(width, height, ij)
        for k, ell, kell in pairs:
            gkl = _moment(width, height, kell)
            pair_a = _mul(gij, gkl)
            pair_b = _mul(
                _moment(width, height, (1 << i) ^ (1 << k)),
                _moment(width, height, (1 << j) ^ (1 << ell)),
            )
            pair_c = _mul(
                _moment(width, height, (1 << i) ^ (1 << ell)),
                _moment(width, height, (1 << j) ^ (1 << k)),
            )
            m4 = _moment(width, height, ij ^ kell)
            u4 = _series_subtract((1, m4), (-1, pair_a), (-1, pair_b), (-1, pair_c))
            _add_scaled(u_square, _square(u4))
            _add_scaled(g_u, _mul(pair_a, u4))
            _add_scaled(g_cross, _mul(pair_a, pair_b))
            # The direct same-gap plus adjacent-gap cumulant formula, used here
            # as an independently summed algebraic comparison.
            _add_scaled(direct, _square(m4), 1)
            _add_scaled(direct, _mul(pair_a, m4), 6)
            _add_scaled(direct, _square(pair_a), -9)
    return tuple(overlap), tuple(u_square), tuple(g_u), tuple(g_cross), tuple(direct)


def _closed_c4_flm() -> tuple[tuple[RationalSeries, ...], dict[tuple[int, int], tuple[RationalSeries, ...]]]:
    """Independent 2D FLM for c2 and the closed c4 correlation contractions."""
    span = ORDER // 2
    shapes = tuple(
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
    weights: dict[tuple[int, int], tuple[RationalSeries, ...]] = {}
    bulk = [[Fraction(0) for _ in range(ORDER + 1)] for _ in range(5)]
    for width, height in shapes:
        overlap, u_square, g_u, g_cross, direct_scaled = _closed_c4_rectangle(width, height)
        local = [
            [Fraction(value, 2) for value in overlap],
            [Fraction(value, 24) for value in u_square],
            [Fraction(value, 2) for value in g_u],
            [Fraction(3 * value, 4) for value in g_cross],
            [Fraction(value, 24) for value in direct_scaled],
        ]
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for quantity in range(len(local)):
                    for degree, coefficient in enumerate(subweight[quantity]):
                        local[quantity][degree] -= placements * coefficient
        exact = tuple(tuple(values) for values in local)
        weights[(width, height)] = exact
        for quantity, values in enumerate(exact):
            for degree, coefficient in enumerate(values):
                bulk[quantity][degree] += coefficient
    return tuple(tuple(values) for values in bulk), weights


def _fractions(values: Sequence[str]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def _strings(values: Sequence[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def main() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    primary_c6 = json.loads(PRIMARY_C6.read_text(encoding="utf-8"))
    primary_c4 = json.loads(PRIMARY_C4.read_text(encoding="utf-8"))

    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e84_interlayer_closed.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert all(check["passed"] for check in payload["checks"])

    data = payload["data"]
    assert data["v_order"] == ORDER
    independent = data["independent_c6_even_subgraph_flm"]
    assert independent["tag"] == "[COMPUTATION]"
    assert independent["odd_vertical_columns_zero"]
    assert independent["box_count"] == len(_box_shapes())
    assert independent["plane_span_bound"] == 6
    assert independent["maximum_layers"] == MAX_LAYERS

    primary_residual = _fractions(
        primary_c6["data"]["c6_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]
    )
    artifact_c6 = _fractions(independent["residual_c6_w"])
    assert artifact_c6 == primary_residual

    rebuilt_c6, c6_weights = _independent_c6_flm()
    assert rebuilt_c6 == artifact_c6
    assert rebuilt_c6 == primary_residual
    assert rebuilt_c6 == (
        Fraction(0), Fraction(0), Fraction(2), Fraction(0), Fraction(138), Fraction(0),
        Fraction(13682, 3), Fraction(0), Fraction(109718), Fraction(0), Fraction(2061698),
        Fraction(0), Fraction(32654478),
    )
    for layers in range(1, MAX_LAYERS + 1):
        assert _sum_height(c6_weights, layers) == _fractions(independent["exact_vertical_extent"][str(layers)])
    assert tuple(
        sum((_sum_height(c6_weights, layers)[degree] for layers in range(1, MAX_LAYERS + 1)), Fraction(0))
        for degree in range(ORDER + 1)
    ) == rebuilt_c6
    spin_dos = anisotropic_box_even_subgraph((2, 2, 4), ORDER, 6)
    spin_dos_columns = tuple(
        tuple(int(spin_dos[v_degree][w_degree]) for v_degree in range(ORDER + 1))
        for w_degree in range(7)
    )
    assert all(not any(spin_dos_columns[w_degree]) for w_degree in (1, 3, 5))
    normalized_spin_dos = tuple(
        _divide_unit(spin_dos_columns[w_degree], spin_dos_columns[0])
        for w_degree in (0, 2, 4, 6)
    )
    assert normalized_spin_dos == _normalized_slab_columns(2, 2, 4)


    closed = data["c4_closed_form"]
    assert closed["tag"] == "[THEOREM]"
    assert closed["finite_volume_derivation_verified"]
    assert closed["formal_series_scope"]
    expected_c4_q = _fractions(primary_c4["data"]["c4_direct_coupling_q"]["full_infinite_stack"])
    expected_c4_w = _fractions(primary_c4["data"]["c4_wave_variable_w_tanh_Kz"]["total"])
    artifact_c2 = _fractions(closed["c2_direct_q"])
    artifact_c4_q = _fractions(closed["c4_direct_q"])
    artifact_c4_w = _fractions(closed["c4_total_w"])
    artifact_terms = tuple(
        _fractions(closed["direct_q_terms"][key])
        for key in ("U4_squared_over_24", "G_pair_U4_over_2", "G_cross_over_4_times_3")
    )
    assert artifact_c4_q == expected_c4_q
    assert artifact_c4_w == expected_c4_w
    assert artifact_c4_w == tuple(
        artifact_c4_q[degree] + Fraction(2, 3) * artifact_c2[degree]
        for degree in range(ORDER + 1)
    )
    assert artifact_c4_q == tuple(sum(values[degree] for values in artifact_terms) for degree in range(ORDER + 1))

    rebuilt_closed, _ = _closed_c4_flm()
    rebuilt_c2, rebuilt_u2, rebuilt_gu, rebuilt_cross, rebuilt_direct = rebuilt_closed
    assert rebuilt_c2 == artifact_c2
    assert rebuilt_u2 == artifact_terms[0]
    assert rebuilt_gu == artifact_terms[1]
    assert rebuilt_cross == artifact_terms[2]
    assert rebuilt_direct == artifact_c4_q
    assert rebuilt_direct == tuple(
        rebuilt_u2[degree] + rebuilt_gu[degree] + rebuilt_cross[degree]
        for degree in range(ORDER + 1)
    )
    assert tuple(
        rebuilt_direct[degree] + Fraction(2, 3) * rebuilt_c2[degree]
        for degree in range(ORDER + 1)
    ) == artifact_c4_w

    witnesses = closed["exact_2d_transfer"]["finite_rectangle_witnesses"]
    assert witnesses["rectangle"] == [3, 3]
    assert _fractions(witnesses["G_0_1"]) == _moment(3, 3, (1 << 0) ^ (1 << 1))
    assert _fractions(witnesses["G_0_4"]) == _moment(3, 3, (1 << 0) ^ (1 << 4))
    witness_u4 = _series_subtract(
        (1, _moment(3, 3, (1 << 0) ^ (1 << 1) ^ (1 << 3) ^ (1 << 4))),
        (-1, _mul(_moment(3, 3, 3), _moment(3, 3, (1 << 3) ^ (1 << 4)))),
        (-1, _mul(_moment(3, 3, (1 << 0) ^ (1 << 3)), _moment(3, 3, (1 << 1) ^ (1 << 4)))),
        (-1, _mul(_moment(3, 3, (1 << 0) ^ (1 << 4)), _moment(3, 3, (1 << 1) ^ (1 << 3)))),
    )
    assert _fractions(witnesses["U4_0_1_3_4"]) == witness_u4

    assert _strings(rebuilt_c6) == independent["residual_c6_w"]
    print("PASS test_interlayer_closed")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_interlayer_closed: {error}")
        raise
