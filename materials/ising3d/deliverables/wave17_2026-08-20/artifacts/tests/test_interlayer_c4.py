"""Independent exact controls for the fourth-order interlayer expansion.

The lower-order recomputation deliberately does not import the experiment.  It
uses the repository's spin-density-of-states transform for open 3D boxes,
extracts ``[w^2]log P`` and ``[w^4]log P``, and performs its own three-dimensional
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
ARTIFACT = ROOT / "results" / "interlayer" / "c4_series.json"
TEST_ORDER = 6

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
        quotient[degree] = (
            num[degree]
            - sum(
                (
                    den[shift] * quotient[degree - shift]
                    for shift in range(1, degree + 1)
                ),
                Fraction(0),
            )
        ) / den[0]
    return tuple(quotient)


def _add_scaled(
    target: list[Fraction], source: Sequence[Fraction], scale: int
) -> None:
    for degree, coefficient in enumerate(source):
        if coefficient:
            target[degree] += scale * coefficient


@lru_cache(maxsize=None)
def _box_log_coefficients(
    width: int, height: int, layers: int, order: int
) -> tuple[RationalSeries, RationalSeries, bool]:
    """Independent spin-enumeration transform for one open 3D box."""

    # The open square directions are isomorphic after swapping their lengths.
    width, height = sorted((width, height))
    if layers == 1:
        zero = tuple(_zero(order))
        return zero, zero, True
    polynomial = anisotropic_box_even_subgraph(
        (width, height, layers), order, 4
    )
    columns = tuple(
        tuple(Fraction(polynomial[v_degree][w_degree]) for v_degree in range(order + 1))
        for w_degree in range(5)
    )
    odd_zero = not any(columns[1]) and not any(columns[3])
    ratio2 = _divide(columns[2], columns[0], order)
    ratio4 = _divide(columns[4], columns[0], order)
    ratio2_square = _mul(ratio2, ratio2, order)
    log_w4 = tuple(
        ratio4[degree] - ratio2_square[degree] / 2 for degree in range(order + 1)
    )
    return ratio2, log_w4, odd_zero


def _brute_force_slab_flm(order: int) -> tuple[RationalSeries, RationalSeries, bool]:
    """Return residual c2 and c4 from exact spin-enumerated 3D boxes."""

    span_budget = order // 2
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width in range(1, span_budget + 2)
                for height in range(1, span_budget + 2)
                if (width - 1) + (height - 1) <= span_budget
                for layers in (1, 2, 3)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    weights: dict[tuple[int, int, int], tuple[RationalSeries, RationalSeries]] = {}
    bulk = (_zero(order), _zero(order))
    odd_zero = True
    for width, height, layers in shapes:
        log_w2, log_w4, box_odd_zero = _box_log_coefficients(
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
        exact = (tuple(values[0]), tuple(values[1]))
        weights[(width, height, layers)] = exact
        _add_scaled(bulk[0], exact[0], 1)
        _add_scaled(bulk[1], exact[1], 1)
    return tuple(bulk[0]), tuple(bulk[1]), odd_zero


def _fractions(values: Sequence[str]) -> RationalSeries:
    return tuple(Fraction(value) for value in values)


def main() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["provenance"]["script"] == "experiments/e51_interlayer_c4.py"
    assert payload["provenance"]["interpreter"] == ".venv/bin/python"
    assert all(check["passed"] for check in payload["checks"])

    data = payload["data"]
    stored_order = int(data["v_order"])
    assert stored_order >= 8
    stored_c2 = _fractions(data["reproduced_c2"]["total_q_or_w"])
    stored_c4_residual_w = _fractions(
        data["c4_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]
    )
    stored_c4_total_w = _fractions(data["c4_wave_variable_w_tanh_Kz"]["total"])
    stored_c4_q = _fractions(
        data["c4_direct_coupling_q"]["full_infinite_stack"]
    )
    stored_c4_q_second_route = _fractions(
        data["c4_direct_coupling_q"]["from_anisotropic_route"]
    )

    assert len(stored_c2) == stored_order + 1
    assert len(stored_c4_residual_w) == stored_order + 1
    assert len(stored_c4_total_w) == stored_order + 1
    assert len(stored_c4_q) == stored_order + 1
    assert stored_c4_q == stored_c4_q_second_route
    assert stored_c2[8] == 778
    assert stored_c4_residual_w[8] == 14223

    residual_c2, residual_c4_w, odd_zero = _brute_force_slab_flm(TEST_ORDER)
    total_c2 = list(residual_c2)
    total_c2[0] += Fraction(1, 2)
    total_c4_w = list(residual_c4_w)
    total_c4_w[0] += Fraction(1, 4)
    direct_c4_q = tuple(
        total_c4_w[degree] - Fraction(2, 3) * total_c2[degree]
        for degree in range(TEST_ORDER + 1)
    )

    assert odd_zero
    assert tuple(total_c2) == (
        Fraction(1, 2),
        Fraction(0),
        Fraction(2),
        Fraction(0),
        Fraction(18),
        Fraction(0),
        Fraction(118),
    )
    assert residual_c4_w == (
        Fraction(0),
        Fraction(0),
        Fraction(2),
        Fraction(0),
        Fraction(63),
        Fraction(0),
        Fraction(1050),
    )
    assert direct_c4_q == (
        Fraction(-1, 12),
        Fraction(0),
        Fraction(2, 3),
        Fraction(0),
        Fraction(51),
        Fraction(0),
        Fraction(2914, 3),
    )
    assert tuple(total_c2) == stored_c2[: TEST_ORDER + 1]
    assert residual_c4_w == stored_c4_residual_w[: TEST_ORDER + 1]
    assert tuple(total_c4_w) == stored_c4_total_w[: TEST_ORDER + 1]
    assert direct_c4_q == stored_c4_q[: TEST_ORDER + 1]

    assert _fractions(data["c3"]["coefficients_through_v12"]) == (
        Fraction(0),
    ) * (stored_order + 1)
    provenance_rows = data["order_by_order_provenance"]
    assert len(provenance_rows) == stored_order + 1
    for degree, row in enumerate(provenance_rows):
        assert row["v_degree"] == degree
        assert Fraction(row["coefficient_c4_residual_w"]) == stored_c4_residual_w[degree]
        assert Fraction(row["coefficient_c4_total_w"]) == stored_c4_total_w[degree]
        assert Fraction(row["coefficient_c4_direct_q"]) == stored_c4_q[degree]
        assert Fraction(row["cumulant_route"]) == Fraction(
            row["anisotropic_route_after_exact_variable_change"]
        )

    print("PASS test_interlayer_c4")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_interlayer_c4: {error}")
        raise
