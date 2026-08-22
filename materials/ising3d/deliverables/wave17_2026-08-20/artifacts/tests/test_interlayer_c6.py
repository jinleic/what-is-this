"""Independent exact controls for the sixth-order interlayer coefficient.

This test deliberately does not import the experiment.  It transforms exact
spin densities of states of open anisotropic 3D boxes into ``P(v,w)``, extracts
``[w^2]``, ``[w^4]``, and ``[w^6] log P``, and performs its own three-dimensional
finite-lattice inversion through ``v^6``.
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Sequence

from ising.interlayer import anisotropic_box_even_subgraph

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "interlayer" / "c6_series.json"
TEST_ORDER = 6
MAX_W_ORDER = 6
MAX_LAYERS = 4

RationalSeries = tuple[Fraction, ...]


def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


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
    den = [Fraction(value) for value in denominator]
    num = [Fraction(value) for value in numerator]
    if not den or not den[0]:
        raise ValueError("zero formal denominator")
    quotient = _zero(order)
    for degree in range(order + 1):
        convolution = sum(
            (
                den[shift] * quotient[degree - shift]
                for shift in range(1, degree + 1)
            ),
            Fraction(0),
        )
        quotient[degree] = (num[degree] - convolution) / den[0]
    return tuple(quotient)


def _add_scaled(
    target: list[Fraction], source: Sequence[int | Fraction], scale: int
) -> None:
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += scale * coefficient


@lru_cache(maxsize=None)
def _box_log_coefficients(
    width: int, height: int, layers: int, order: int
) -> tuple[RationalSeries, RationalSeries, RationalSeries, bool]:
    """Extract even vertical log columns from an exact spin-DOS transform."""

    width, height = sorted((width, height))
    if layers == 1:
        zero = tuple(_zero(order))
        return zero, zero, zero, True
    polynomial = anisotropic_box_even_subgraph(
        (width, height, layers), order, MAX_W_ORDER
    )
    columns = tuple(
        tuple(Fraction(polynomial[v_degree][w_degree]) for v_degree in range(order + 1))
        for w_degree in range(MAX_W_ORDER + 1)
    )
    odd_zero = all(not any(columns[degree]) for degree in (1, 3, 5))
    ratio2 = _divide(columns[2], columns[0], order)
    ratio4 = _divide(columns[4], columns[0], order)
    ratio6 = _divide(columns[6], columns[0], order)
    ratio2_square = _mul(ratio2, ratio2, order)
    log4 = tuple(
        ratio4[degree] - ratio2_square[degree] / 2
        for degree in range(order + 1)
    )
    ratio2_ratio4 = _mul(ratio2, ratio4, order)
    ratio2_cube = _mul(ratio2_square, ratio2, order)
    log6 = tuple(
        ratio6[degree]
        - ratio2_ratio4[degree]
        + ratio2_cube[degree] / 3
        for degree in range(order + 1)
    )
    return ratio2, log4, log6, odd_zero


def _independent_spin_dos_flm(
    order: int,
) -> tuple[
    tuple[RationalSeries, RationalSeries, RationalSeries],
    dict[tuple[int, int, int], tuple[RationalSeries, RationalSeries, RationalSeries]],
    bool,
]:
    """Recompute the residual bulk coefficients from open 3D spin DOS."""

    span_budget = order // 2
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width in range(1, span_budget + 2)
                for height in range(1, span_budget + 2)
                if (width - 1) + (height - 1) <= span_budget
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
        c2, c4, c6, box_odd_zero = _box_log_coefficients(
            width, height, layers, order
        )
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
            _add_scaled(bulk[quantity], exact[quantity], 1)
    return tuple(tuple(values) for values in bulk), weights, odd_zero


def _sum_height(
    weights: dict[tuple[int, int, int], tuple[RationalSeries, ...]],
    layers: int,
    quantity: int,
    order: int,
) -> RationalSeries:
    total = _zero(order)
    for shape, values in weights.items():
        if shape[2] == layers:
            _add_scaled(total, values[quantity], 1)
    return tuple(total)


def _fractions(values: Sequence[str]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def main() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e64_interlayer_c6.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert all(check["passed"] for check in payload["checks"])

    data = payload["data"]
    stored_order = int(data["v_order"])
    assert stored_order == 12
    assert int(data["independent_check_v_order"]) >= TEST_ORDER

    symbolic = data["symbolic_sixth_cumulant"]
    assert symbolic["passed"]
    assert symbolic["bell_number_B6"] == 203
    assert sum(row["partition_count"] for row in symbolic["set_partition_profiles"]) == 203
    assert symbolic["zero_mean_surviving_formula"] == "kappa6=mu6-15*mu4*mu2+30*mu2^3"
    assert [
        row["ordered_gap_counts"]
        for row in symbolic["connected_adjacent_gap_profiles"]
    ] == [[6], [4, 2], [2, 4], [2, 2, 2]]
    assert symbolic["disconnected_gap_symbolic_cancellation"]

    lower = data["reproduced_lower_orders"]
    c2 = _fractions(lower["c2_total_q_or_w"])
    c4_q = _fractions(lower["c4_direct_q"])
    c4_w_residual = _fractions(lower["c4_residual_w"])
    c6_q_data = data["c6_direct_coupling_q"]
    c6_q = _fractions(c6_q_data["full_infinite_stack"])
    c6_q_h2 = _fractions(c6_q_data["same_gap_height2"])
    c6_q_h3 = _fractions(c6_q_data["two_adjacent_gaps_height3"])
    c6_q_h4 = _fractions(c6_q_data["three_adjacent_gaps_height4"])
    c6_w_data = data["c6_wave_variable_w_tanh_Kz"]
    c6_w_total = _fractions(c6_w_data["total"])
    c6_w_residual = _fractions(c6_w_data["residual_after_log_cosh"])

    assert all(len(values) == stored_order + 1 for values in (c2, c4_q, c6_q, c6_w_total))
    assert c6_q == tuple(
        c6_q_h2[degree] + c6_q_h3[degree] + c6_q_h4[degree]
        for degree in range(stored_order + 1)
    )
    assert c6_w_total == tuple(
        c6_q[degree]
        + Fraction(4, 3) * c4_q[degree]
        + Fraction(23, 45) * c2[degree]
        for degree in range(stored_order + 1)
    )
    assert c6_w_total[0] == Fraction(1, 6)
    assert c6_w_residual[0] == 0
    assert c6_w_total[1:] == c6_w_residual[1:]

    assert c2[: TEST_ORDER + 1] == (
        Fraction(1, 2), Fraction(0), Fraction(2), Fraction(0),
        Fraction(18), Fraction(0), Fraction(118),
    )
    assert c4_w_residual[: TEST_ORDER + 1] == (
        Fraction(0), Fraction(0), Fraction(2), Fraction(0),
        Fraction(63), Fraction(0), Fraction(1050),
    )
    assert c6_q[: TEST_ORDER + 1] == (
        Fraction(1, 45), Fraction(0), Fraction(4, 45), Fraction(0),
        Fraction(304, 5), Fraction(0), Fraction(144236, 45),
    )
    assert c6_w_residual[: TEST_ORDER + 1] == (
        Fraction(0), Fraction(0), Fraction(2), Fraction(0),
        Fraction(138), Fraction(0), Fraction(13682, 3),
    )
    assert c6_w_residual[8] == 109718
    assert c6_w_residual[12] == 32654478

    (direct_c2, direct_c4, direct_c6), weights, odd_zero = _independent_spin_dos_flm(
        TEST_ORDER
    )
    assert odd_zero
    assert direct_c2 == _fractions(lower["c2_residual_w"])[: TEST_ORDER + 1]
    assert direct_c4 == c4_w_residual[: TEST_ORDER + 1]
    assert direct_c6 == c6_w_residual[: TEST_ORDER + 1]

    stored_heights = c6_w_data["residual_by_exact_vertical_extent"]
    for layers in (2, 3, 4):
        direct_height = _sum_height(weights, layers, 2, TEST_ORDER)
        assert direct_height == _fractions(stored_heights[str(layers)])[: TEST_ORDER + 1]
    assert tuple(
        sum((_sum_height(weights, layers, 2, TEST_ORDER)[degree] for layers in (2, 3, 4)), Fraction(0))
        for degree in range(TEST_ORDER + 1)
    ) == direct_c6

    odd = data["odd_interlayer_orders"]
    assert _fractions(odd["c3_through_v12"]) == (Fraction(0),) * (stored_order + 1)
    assert _fractions(odd["c5_through_v12"]) == (Fraction(0),) * (stored_order + 1)

    provenance_rows = data["order_by_order_provenance"]
    assert len(provenance_rows) == stored_order + 1
    for degree, row in enumerate(provenance_rows):
        assert row["v_degree"] == degree
        assert Fraction(row["c6_direct_q"]) == c6_q[degree]
        assert Fraction(row["c6_total_w"]) == c6_w_total[degree]
        assert Fraction(row["c6_residual_w"]) == c6_w_residual[degree]
        if degree <= TEST_ORDER:
            assert Fraction(row["direct_spin_dos_residual_w"]) == direct_c6[degree]
            assert row["two_route_agreement_available"]
        else:
            assert row["direct_spin_dos_residual_w"] is None
            assert not row["two_route_agreement_available"]

    print("PASS test_interlayer_c6")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_interlayer_c6: {error}")
        raise
