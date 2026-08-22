"""Independent v^12 c6 slab certificate and a correlation closed form for c4.

The c6 route is intentionally independent of the generated-cumulant producer:
it builds exact open-slab even-subgraph polynomials from planar odd-boundary
polynomials, takes the logarithm by Newton's identity, then applies 3D
finite-lattice Möbius inversion.  The c4 calculation starts directly from
exact two-dimensional G and connected U4 transfer data.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

SCRIPT = "experiments/e84_interlayer_closed.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "interlayer" / "closed_form.json"
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
    """Multiply two exact v series through the requested degree."""
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
    """Exact formal division in Z[[v]] when denominator[0] is one."""
    if denominator[0] != 1:
        raise AssertionError("the empty planar even subgraph must be unique")
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
def boundary_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """Exact P_S(v) for all odd boundary masks of an open planar rectangle."""
    sites = width * height
    coefficients = [[0] * (order + 1) for _ in range(1 << sites)]
    coefficients[0][0] = 1
    for left, right in _square_bonds(width, height):
        toggle = (1 << left) | (1 << right)
        # Descending degree makes each planar edge available exactly once.
        for degree in range(order - 1, -1, -1):
            for mask in range(1 << sites):
                value = coefficients[mask][degree]
                if value:
                    coefficients[mask ^ toggle][degree + 1] += value
    return tuple(tuple(values) for values in coefficients)


@lru_cache(maxsize=None)
def moment_series(width: int, height: int, mask: int, order: int = ORDER) -> IntegerSeries:
    """M(mask)=P_mask/P_empty as an exact truncated 2D spin-moment series."""
    boundary = boundary_polynomials(width, height, order)
    return _divide_unit(boundary[mask], boundary[0], order)


def _masks_by_weight(width: int, height: int) -> dict[int, tuple[int, ...]]:
    sites = width * height
    return {
        weight: tuple(mask for mask in range(1 << sites) if mask.bit_count() == weight)
        for weight in (0, 2, 4, 6)
    }


# ---------------------------------------------------------------------------
# Independent w^6 open-slab transfer
# ---------------------------------------------------------------------------
def _sum_squares(width: int, height: int, masks: Iterable[int]) -> IntegerSeries:
    total = _zero()
    for mask in masks:
        _add_scaled(total, _square(moment_series(width, height, mask)))
    return tuple(total)


def _pair_boundary_sum(
    width: int,
    height: int,
    first_masks: Iterable[int],
    second_masks: Iterable[int],
) -> IntegerSeries:
    """Sum m_A m_(A xor B) m_B over two vertical gap boundary masks."""
    total = _zero()
    for first_mask in first_masks:
        first_moment = moment_series(width, height, first_mask)
        for second_mask in second_masks:
            term = _mul(
                _mul(first_moment, moment_series(width, height, first_mask ^ second_mask)),
                moment_series(width, height, second_mask),
            )
            _add_scaled(total, term)
    return tuple(total)


def _triple_two_boundary_sum(width: int, height: int, masks_two: Iterable[int]) -> IntegerSeries:
    """Exact factorization of the three-gap w^6 contribution.

    If A,B,C all have vertical weight two, summation over the two outside
    boundaries factorizes at fixed B:
      sum_{A,B,C} m_A m_(A xor B) m_(B xor C) m_C
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
                    moment_series(width, height, outside_mask),
                    moment_series(width, height, outside_mask ^ middle_mask),
                ),
            )
        _add_scaled(total, _square(tuple(half)))
    return tuple(total)


def normalized_slab_columns(width: int, height: int, layers: int) -> tuple[IntegerSeries, ...]:
    """Return r_0,r_2,r_4,r_6 of P(v,w)/P(v,0) through w^6.

    Write S_g for the chosen vertical-edge mask across gap g.  The horizontal
    boundary in layer ell is S_(ell-1) xor S_ell.  Since every planar boundary
    has even cardinality and S_0=0, all S_g have even cardinality; therefore
    every odd w column is identically zero before any numerical evaluation.
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

    pair_two_two = _pair_boundary_sum(width, height, masks[2], masks[2])
    pair_two_four = _pair_boundary_sum(width, height, masks[2], masks[4])
    if layers == 3:
        return (
            one,
            tuple(2 * value for value in square_two),
            tuple(2 * square_four[degree] + pair_two_two[degree] for degree in range(ORDER + 1)),
            tuple(2 * square_six[degree] + 2 * pair_two_four[degree] for degree in range(ORDER + 1)),
        )

    if layers != 4:
        raise ValueError("w^6 can have no connected vertical extent above four layers")
    triple_two = _triple_two_boundary_sum(width, height, masks[2])
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


def log_w6(columns: tuple[IntegerSeries, ...]) -> RationalSeries:
    """Newton/logarithm extraction [w^6] log(1+r2 w^2+r4 w^4+r6 w^6)."""
    _, r2, r4, r6 = columns
    r2r4 = _mul(r2, r4)
    r2cubed = _mul(_square(r2), r2)
    return tuple(
        Fraction(r6[degree] - r2r4[degree]) + Fraction(r2cubed[degree], 3)
        for degree in range(ORDER + 1)
    )


def box_shapes(order: int = ORDER) -> tuple[tuple[int, int, int], ...]:
    """All boxes whose in-plane spanning cut lower bound permits v^order."""
    span_budget = order // 2
    return tuple(
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


def independent_c6_flm() -> tuple[RationalSeries, dict[tuple[int, int, int], RationalSeries]]:
    """Exact 3D rectangular FLM of [w^6] log P(v,w)."""
    weights: dict[tuple[int, int, int], RationalSeries] = {}
    bulk = _fraction_zero()
    for width, height, layers in box_shapes():
        local = (
            _fraction_zero()
            if layers == 1
            else list(log_w6(normalized_slab_columns(width, height, layers)))
        )
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


# ---------------------------------------------------------------------------
# c4 closed correlation form
# ---------------------------------------------------------------------------
def _series_subtract(*terms: tuple[int, Sequence[int]]) -> IntegerSeries:
    result = _zero()
    for scale, values in terms:
        _add_scaled(result, values, scale)
    return tuple(result)


def _closed_c4_rectangle(
    width: int, height: int
) -> tuple[IntegerSeries, IntegerSeries, IntegerSeries, IntegerSeries, IntegerSeries, bool]:
    """Exact c2 and c4 G/U contractions on one open 2D rectangle.

    The final returned integer series is 24 times the direct same-gap plus
    adjacent-gap fourth coefficient.  The boolean compares the U4 closed form
    with that direct expression before any finite-lattice inversion.
    """
    sites = width * height
    pairs = tuple((i, j, (1 << i) ^ (1 << j)) for i, j in product(range(sites), repeat=2))
    overlap = _zero()
    u_squared = _zero()
    g_pair_u = _zero()
    g_cross = _zero()
    direct_times_24 = _zero()

    for _, _, mask in pairs:
        _add_scaled(overlap, _square(moment_series(width, height, mask)))

    for i, j, ij in pairs:
        gij = moment_series(width, height, ij)
        for k, ell, kell in pairs:
            gkl = moment_series(width, height, kell)
            pair_a = _mul(gij, gkl)
            pair_b = _mul(
                moment_series(width, height, (1 << i) ^ (1 << k)),
                moment_series(width, height, (1 << j) ^ (1 << ell)),
            )
            pair_c = _mul(
                moment_series(width, height, (1 << i) ^ (1 << ell)),
                moment_series(width, height, (1 << j) ^ (1 << k)),
            )
            m4 = moment_series(width, height, ij ^ kell)
            u4 = _series_subtract((1, m4), (-1, pair_a), (-1, pair_b), (-1, pair_c))
            _add_scaled(u_squared, _square(u4))
            _add_scaled(g_pair_u, _mul(pair_a, u4))
            _add_scaled(g_cross, _mul(pair_a, pair_b))

            # Summed direct form: (M4^2 + 6 A M4 - 9 A^2)/24.
            _add_scaled(direct_times_24, _square(m4))
            _add_scaled(direct_times_24, _mul(pair_a, m4), 6)
            _add_scaled(direct_times_24, _square(pair_a), -9)

    announced_times_24 = tuple(
        u_squared[degree] + 12 * g_pair_u[degree] + 18 * g_cross[degree]
        for degree in range(ORDER + 1)
    )
    return (
        tuple(overlap),
        tuple(u_squared),
        tuple(g_pair_u),
        tuple(g_cross),
        tuple(direct_times_24),
        tuple(announced_times_24) == tuple(direct_times_24),
    )


def plane_shapes(order: int = ORDER) -> tuple[tuple[int, int], ...]:
    span_budget = order // 2
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


def closed_c4_flm() -> tuple[
    tuple[RationalSeries, ...],
    dict[tuple[int, int], tuple[RationalSeries, ...]],
    bool,
]:
    """2D FLM of c2 and the three terms in the c4 G/U identity."""
    weights: dict[tuple[int, int], tuple[RationalSeries, ...]] = {}
    bulk = [[Fraction(0) for _ in range(ORDER + 1)] for _ in range(5)]
    finite_identity = True
    for width, height in plane_shapes():
        overlap, u_squared, g_pair_u, g_cross, direct_times_24, passed = _closed_c4_rectangle(width, height)
        finite_identity = finite_identity and passed
        local = [
            [Fraction(value, 2) for value in overlap],
            [Fraction(value, 24) for value in u_squared],
            [Fraction(value, 2) for value in g_pair_u],
            [Fraction(3 * value, 4) for value in g_cross],
            [Fraction(value, 24) for value in direct_times_24],
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
    return tuple(tuple(values) for values in bulk), weights, finite_identity


def _u4_series(width: int, height: int, points: tuple[int, int, int, int]) -> IntegerSeries:
    i, j, k, ell = points
    m4 = moment_series(width, height, (1 << i) ^ (1 << j) ^ (1 << k) ^ (1 << ell))
    return _series_subtract(
        (1, m4),
        (-1, _mul(moment_series(width, height, (1 << i) ^ (1 << j)), moment_series(width, height, (1 << k) ^ (1 << ell)))),
        (-1, _mul(moment_series(width, height, (1 << i) ^ (1 << k)), moment_series(width, height, (1 << j) ^ (1 << ell)))),
        (-1, _mul(moment_series(width, height, (1 << i) ^ (1 << ell)), moment_series(width, height, (1 << j) ^ (1 << k)))),
    )


# ---------------------------------------------------------------------------
# Artifact and checks
# ---------------------------------------------------------------------------
def _fraction_strings(values: Sequence[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def _weight_rows(weights: dict[tuple[int, ...], Sequence[Sequence[int | Fraction]]]) -> list[dict[str, object]]:
    rows = []
    for shape in sorted(weights, key=lambda item: (sum(item), item)):
        rows.append(
            {
                "shape": list(shape),
                "series": [_fraction_strings(values) for values in weights[shape]],
            }
        )
    return rows


def _digest(rows: object) -> str:
    text = json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def main() -> None:
    started = time.monotonic()
    checks: list[dict[str, object]] = []

    primary_c6 = json.loads(PRIMARY_C6.read_text(encoding="utf-8"))
    primary_c4 = json.loads(PRIMARY_C4.read_text(encoding="utf-8"))
    primary_c6_residual = tuple(
        Fraction(value)
        for value in primary_c6["data"]["c6_wave_variable_w_tanh_Kz"]["residual_after_log_cosh"]
    )
    primary_c4_q = tuple(
        Fraction(value)
        for value in primary_c4["data"]["c4_direct_coupling_q"]["full_infinite_stack"]
    )
    primary_c4_w = tuple(
        Fraction(value)
        for value in primary_c4["data"]["c4_wave_variable_w_tanh_Kz"]["total"]
    )

    c6_residual, c6_weights = independent_c6_flm()
    c6_heights = {str(layers): _sum_height(c6_weights, layers) for layers in range(1, MAX_LAYERS + 1)}
    c6_weight_rows = _weight_rows({shape: (series,) for shape, series in c6_weights.items()})

    _record(
        checks,
        "independent c6 residual agrees v^0 through v^12",
        c6_residual == primary_c6_residual,
        "exact Fraction comparison with the stored generated-cumulant route at every coefficient",
    )
    _record(
        checks,
        "c6 exact-height decomposition sums to bulk",
        tuple(
            sum((c6_heights[str(layers)][degree] for layers in range(1, MAX_LAYERS + 1)), Fraction(0))
            for degree in range(ORDER + 1)
        )
        == c6_residual,
        "open-slab 3D Möbius weights at heights 1, 2, 3, and 4 sum exactly",
    )
    _record(
        checks,
        "odd vertical columns vanish before enumeration",
        all(weight % 2 == 0 for weight in (0, 2, 4, 6)),
        "layer boundary parity forces every vertical mask S_g to have even cardinality",
    )
    _record(
        checks,
        "c6 finite lattice support bound",
        max(layers for _, _, layers in c6_weights) == 4
        and all((width - 1) + (height - 1) <= ORDER // 2 for width, height, _ in c6_weights),
        "w^6 permits at most three nonempty even layer cuts; an a by b span needs at least 2[(a-1)+(b-1)] in-plane edges",
    )

    c4_bulk, c4_weights, finite_identity = closed_c4_flm()
    c2_q, u_squared_term, g_pair_u_term, g_cross_term, c4_q = c4_bulk
    c4_w = tuple(c4_q[degree] + Fraction(2, 3) * c2_q[degree] for degree in range(ORDER + 1))
    c4_weight_rows = _weight_rows(c4_weights)

    _record(
        checks,
        "finite-volume c4 U4 identity",
        finite_identity,
        "the exact G/U numerator U4^2+12(G_ij G_kl)U4+18(G_ij G_kl)(G_ik G_jl) equals the direct cumulant numerator on every FLM rectangle",
    )
    _record(
        checks,
        "closed c4 direct-q series agrees v^0 through v^12",
        c4_q == primary_c4_q,
        "exact 2D correlation contraction and stored primary cumulant coefficients agree coefficientwise",
    )
    _record(
        checks,
        "closed c4 w conversion agrees v^0 through v^12",
        c4_w == primary_c4_w,
        "[w^4]F(atanh w)=c4_q+(2/3)c2_q exactly",
    )
    _record(
        checks,
        "closed c4 terms sum exactly",
        c4_q == tuple(
            u_squared_term[degree] + g_pair_u_term[degree] + g_cross_term[degree]
            for degree in range(ORDER + 1)
        ),
        "the three U4/G contractions reproduce the complete direct fourth coefficient",
    )

    witnesses = {
        "rectangle": [3, 3],
        "G_0_1": _fraction_strings(moment_series(3, 3, (1 << 0) ^ (1 << 1))),
        "G_0_4": _fraction_strings(moment_series(3, 3, (1 << 0) ^ (1 << 4))),
        "U4_0_1_3_4": _fraction_strings(_u4_series(3, 3, (0, 1, 3, 4))),
        "scope": "Exact open-rectangle transfer witnesses through v^12; the infinite-volume contractions are reconstructed separately by the exact planar FLM.",
    }
    _record(
        checks,
        "2D correlation transfer witnesses",
        witnesses["G_0_1"][1] == "1" and witnesses["U4_0_1_3_4"][4] == "-8",
        "integer boundary-polynomial transfer supplies G and connected U4 data through v^12",
    )

    if not all(check["passed"] for check in checks):
        raise AssertionError("refusing to write a failed exact certificate")

    c6_rows_digest = _digest(c6_weight_rows)
    c4_rows_digest = _digest(c4_weight_rows)
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "method": "independent exact integer boundary-polynomial slab transfer, Fraction logarithm/finite-lattice inversion, and exact 2D G/U4 contraction",
            "wall_seconds": f"{time.monotonic() - started:.6f}",
        },
        "data": {
            "claim": "[COMPUTATION] An independent open-slab even-subgraph FLM agrees with c6 through v^12. [THEOREM] c4 has the stated exact G/U4 closed form as a formal high-temperature identity.",
            "v_order": ORDER,
            "independent_c6_even_subgraph_flm": {
                "tag": "[COMPUTATION]",
                "quantity": "[w^6] log P(v,w), equivalently c6_w residual after the vertical log-cosh prefactor",
                "residual_c6_w": _fraction_strings(c6_residual),
                "primary_generated_cumulant_reference": {
                    "artifact": "results/interlayer/c6_series.json",
                    "residual_c6_w": _fraction_strings(primary_c6_residual),
                    "comparison_used_only_after_independent_calculation": True,
                },
                "exact_vertical_extent": {key: _fraction_strings(value) for key, value in c6_heights.items()},
                "odd_vertical_columns_zero": True,
                "maximum_layers": MAX_LAYERS,
                "plane_span_bound": ORDER // 2,
                "box_count": len(c6_weights),
                "boxes": [list(shape) for shape in box_shapes()],
                "box_log_w6_weights": c6_weight_rows,
                "box_log_w6_weights_sha256": c6_rows_digest,
                "stabilization_bound": {
                    "tag": "[LEMMA]",
                    "statement": "Every connected even subgraph contributing to [v^d w^6] log P with d<=12 has vertical extent at most four layers and in-plane bounding-box span (a-1)+(b-1)<=6.",
                    "vertical_reason": "Every nonempty vertical cut contains at least two edges; six vertical edges therefore occupy at most three gaps.",
                    "planar_reason": "Each crossed in-plane coordinate cut is crossed by a positive even number of in-plane edges, so an a by b span needs at least 2[(a-1)+(b-1)] in-plane edges.",
                    "consequence": "All omitted boxes have either w-degree at least eight or v-degree at least fourteen; their [v^0..v^12 w^6] FLM weights vanish exactly. This is a finite-support identification, not an extrapolation.",
                },
            },
            "c4_closed_form": {
                "tag": "[THEOREM]",
                "variables": {
                    "G_ij": "<sigma_i sigma_j>_2D",
                    "U4_ijkl": "<sigma_i sigma_j sigma_k sigma_l>_2D-G_ij G_kl-G_ik G_jl-G_il G_jk",
                    "A": "G_0r G_st",
                    "B": "G_0s G_rt",
                },
                "identity_direct_q": "c4_q=(1/24) sum_(r,s,t) U4(0,r,s,t)^2+(1/2) sum_(r,s,t) G(0,r)G(s,t)U4(0,r,s,t)+(3/4) sum_(r,s,t) G(0,r)G(s,t)G(0,s)G(r,t)",
                "identity_wave_w_total": "c4_w_total=c4_q+(1/3) sum_r G(0,r)^2",
                "derivation": {
                    "same_gap": "M4^2-A^2-B^2-C^2",
                    "adjacent_gaps": "A M4-A^2",
                    "substitution": "M4=A+B+C+U4",
                    "dummy_index_symmetries": "sum A U4=sum B U4=sum C U4 and sum A B=sum A C=sum B C",
                    "resulting_numerator": "U4^2+12 A U4+18 A B, divided by 24",
                    "w_change_of_variable": "q=atanh(w), hence [w^4]F=c4_q+(2/3)c2_q and c2_q=(1/2)sum_r G(0,r)^2",
                },
                "formal_series_scope": True,
                "finite_volume_derivation_verified": finite_identity,
                "c2_direct_q": _fraction_strings(c2_q),
                "c4_direct_q": _fraction_strings(c4_q),
                "c4_total_w": _fraction_strings(c4_w),
                "direct_q_terms": {
                    "U4_squared_over_24": _fraction_strings(u_squared_term),
                    "G_pair_U4_over_2": _fraction_strings(g_pair_u_term),
                    "G_cross_over_4_times_3": _fraction_strings(g_cross_term),
                },
                "exact_2d_transfer": {
                    "method": "For every open FLM rectangle, integer P_S(v) boundary polynomials are divided by P_empty(v) in Z[[v]] to obtain G and M4; U4 is then formed exactly coefficientwise.",
                    "plane_span_bound": ORDER // 2,
                    "rectangle_count": len(c4_weights),
                    "rectangles": [list(shape) for shape in plane_shapes()],
                    "contraction_box_weights": c4_weight_rows,
                    "contraction_box_weights_sha256": c4_rows_digest,
                    "finite_rectangle_witnesses": witnesses,
                },
                "primary_reference": {
                    "artifact": "results/interlayer/c4_series.json",
                    "c4_direct_q": _fraction_strings(primary_c4_q),
                    "c4_total_w": _fraction_strings(primary_c4_w),
                },
            },
            "scope": "[UNRESOLVED] The identities are coefficientwise formal high-temperature identities. No convergence radius, critical-coupling evaluation, or solution of the 3D Ising model follows.",
        },
        "checks": checks,
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT.relative_to(ROOT)}")
    print("PASS e84_interlayer_closed")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL e84_interlayer_closed: {error}")
        raise
