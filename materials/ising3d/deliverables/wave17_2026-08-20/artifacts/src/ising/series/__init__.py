"""Exact finite-lattice high- and low-temperature Ising series.

The implementation uses Enting's finite-lattice decomposition.  If ``L(A)`` is
the logarithm of the normalized partition polynomial for a free box ``A``,

    L(A) = sum_{R <= A} prod_i (A_i - R_i + 1) W(R),

then ``W(R)`` is the contribution of connected clusters with bounding box
exactly ``R``.  The bulk interaction series is ``sum_R W(R)``.  Equivalently,
``W = Delta_1^2 ... Delta_d^2 L`` after extending ``L`` by zero at non-positive
side lengths.  One mixed *first* difference is only a cumulative sum of
weights; it must not itself be summed over all boxes.

Every polynomial and formal-series operation below uses Python integers and
``fractions.Fraction``.  Floating-point arithmetic is not used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, prod
from typing import Callable, Mapping, Sequence

from ising.transfer_matrix import box_broken_bond_poly

__all__ = [
    "FLMSeries",
    "broken_bond_to_even_subgraph",
    "high_temperature_free_energy",
    "log_series",
    "low_temperature_free_energy",
    "onsager_square_ht_series",
]

RationalSeries = tuple[Fraction, ...]
Shape = tuple[int, ...]


@dataclass(frozen=True)
class FLMSeries:
    """A truncated exact bulk free-energy series.

    ``coefficients`` excludes the non-polynomial constant ``log(2)`` in the
    high-temperature expansion and excludes the ground-state term ``d*K`` in
    the low-temperature expansion.  Thus

    * HT: ``phi = log(2) + sum_n coefficients[n] v**n``;
    * LT: ``phi = d*K + sum_n coefficients[n] x**n``.

    For HT, ``interaction_coefficients`` is the finite-lattice contribution
    before adding ``d*log(cosh K)``.  For LT it equals ``coefficients``.
    ``box_weights`` contains one exact ``W`` series for every ordered box used.
    """

    dimension: int
    order: int
    variable: str
    coefficients: RationalSeries
    interaction_coefficients: RationalSeries
    boxes: tuple[Shape, ...]
    box_weights: Mapping[Shape, RationalSeries]
    bound_budget: int
    bound_slack: int
    order_bound: str


def _check_order(order: int) -> int:
    order = int(order)
    if order < 0:
        raise ValueError("order must be non-negative")
    return order


def _poly_mul_int(
    left: Sequence[int], right: Sequence[int], max_degree: int | None = None
) -> list[int]:
    if max_degree is None:
        max_degree = len(left) + len(right) - 2
    out = [0] * (min(max_degree, len(left) + len(right) - 2) + 1)
    for i, a in enumerate(left):
        if not a or i > max_degree:
            continue
        stop = min(len(right), max_degree - i + 1)
        for j in range(stop):
            b = right[j]
            if b:
                out[i + j] += int(a) * int(b)
    return out


def _box_bond_count(shape: Shape) -> int:
    sites = prod(shape)
    return sum((side - 1) * (sites // side) for side in shape)


def broken_bond_to_even_subgraph(
    broken_coefficients: Sequence[int],
    shape: Sequence[int],
    order: int | None = None,
) -> tuple[int, ...]:
    """Convert a free-box broken-bond polynomial to ``P(v)`` exactly.

    If ``c[q]`` counts configurations with ``q`` unsatisfied bonds, comparison
    of

    ``Z = exp(K*n_b) sum_q c[q] x^q`` and
    ``Z = 2^N cosh(K)^n_b P(v)``

    with ``x=(1-v)/(1+v)`` gives

    ``P(v) = 2^-N sum_q c[q] (1-v)^q (1+v)^(n_b-q)``.

    Every returned numerator is asserted divisible by ``2^N``.  Passing
    ``order`` computes and checks only the requested truncation.
    """

    box = tuple(int(side) for side in shape)
    if len(box) not in (2, 3) or any(side < 1 for side in box):
        raise ValueError("shape must contain two or three positive side lengths")
    coefficients = tuple(int(value) for value in broken_coefficients)
    if any(value < 0 for value in coefficients):
        raise ValueError("broken-bond coefficients must be non-negative")

    sites = prod(box)
    bonds = _box_bond_count(box)
    if len(coefficients) > bonds + 1:
        raise ValueError("broken-bond polynomial degree exceeds the free-box bond count")
    if sum(coefficients) != 1 << sites:
        raise ValueError("broken-bond coefficients do not sum to 2^N")

    limit = bonds if order is None else min(_check_order(order), bonds)
    positive_powers = [[1]]
    negative_powers = [[1]]
    for _ in range(bonds):
        positive_powers.append(_poly_mul_int(positive_powers[-1], [1, 1], limit))
        negative_powers.append(_poly_mul_int(negative_powers[-1], [1, -1], limit))

    accumulated = [0] * (limit + 1)
    for q, count in enumerate(coefficients):
        if not count:
            continue
        term = _poly_mul_int(negative_powers[q], positive_powers[bonds - q], limit)
        for degree, value in enumerate(term):
            accumulated[degree] += count * value

    denominator = 1 << sites
    converted: list[int] = []
    for numerator in accumulated:
        assert numerator % denominator == 0, (
            "non-integral even-subgraph coefficient: broken-bond convention error"
        )
        converted.append(numerator // denominator)
    assert converted[0] == 1, "P(0) must count exactly the empty even subgraph"
    assert all(value >= 0 for value in converted), "even-subgraph counts must be non-negative"
    return tuple(converted)


def log_series(polynomial: Sequence[int | Fraction], order: int) -> RationalSeries:
    """Return the exact truncation of ``log(polynomial)`` through ``order``.

    The input must have constant coefficient one.  From ``A' = (log A)' A``,

    ``ell[n] = a[n] - (1/n) sum_{k=1}^{n-1} k*ell[k]*a[n-k]``.
    """

    order = _check_order(order)
    if not polynomial or Fraction(polynomial[0]) != 1:
        raise ValueError("log_series requires constant coefficient one")
    a = [Fraction(0) for _ in range(order + 1)]
    for degree, value in enumerate(polynomial[: order + 1]):
        a[degree] = Fraction(value)
    logarithm = [Fraction(0) for _ in range(order + 1)]
    for degree in range(1, order + 1):
        convolution = sum(
            (k * logarithm[k] * a[degree - k] for k in range(1, degree)),
            Fraction(0),
        )
        logarithm[degree] = a[degree] - convolution / degree
    return tuple(logarithm)


def _canonical_shape(shape: Shape) -> Shape:
    # Put the largest side in the open transfer direction, minimizing the
    # cross-section state space.  Isotropy makes permutations equivalent.
    return tuple(sorted(shape))


@lru_cache(maxsize=None)
def _free_broken_polynomial(shape: Shape) -> tuple[int, ...]:
    return tuple(box_broken_bond_poly(_canonical_shape(shape)))


@lru_cache(maxsize=None)
def _plus_broken_polynomial(shape: Shape) -> tuple[int, ...]:
    canonical = _canonical_shape(shape)
    coefficients = tuple(box_broken_bond_poly(canonical, plus_boundary=True))
    # NOTE: a local repair for the length-one box used to live here, because the transfer engine
    # applied the two end ghosts on the first and last propagation steps and those steps coincide
    # when c == 1.  That was a genuine bug in ising.transfer_matrix and it has now been fixed at
    # the source (see tests/test_plus_boundary.py, which compares against an independent brute
    # force on 22 boxes including every c == 1 case).  The workaround is therefore removed: one
    # source of truth, not two.  Canonicalisation sorts the shape ascending, so c == 1 occurred
    # only for (1,1,1); the published low-temperature coefficients were computed with the repaired
    # value and are unchanged, which experiments/e26 confirms by an independent route.
    if canonical == (1, 1, 1):
        assert coefficients == (1, 0, 0, 0, 0, 0, 1), (
            f"one-site plus-boundary box must be 1 + x^6, got {coefficients}"
        )
    assert coefficients[0] == 1, (
        "plus-boundary Xi(0) must be one: only the all-up state has no broken bonds"
    )
    return coefficients


@lru_cache(maxsize=None)
def _ht_box_log(shape: Shape, order: int) -> RationalSeries:
    polynomial = broken_bond_to_even_subgraph(
        _free_broken_polynomial(shape), _canonical_shape(shape), order
    )
    return log_series(polynomial, order)


@lru_cache(maxsize=None)
def _lt_box_log(shape: Shape, order: int) -> RationalSeries:
    return log_series(_plus_broken_polynomial(shape), order)


def _ht_shapes(dimension: int, order: int, bound_slack: int) -> tuple[Shape, ...]:
    extent_budget = order // 2 + bound_slack
    shapes = [
        tuple(sides)
        for sides in product(range(1, extent_budget + 2), repeat=dimension)
        if sum(side - 1 for side in sides) <= extent_budget
    ]
    return tuple(sorted(shapes, key=lambda shape: (sum(shape), shape)))


def _lt_shapes(order: int, bound_slack: int) -> tuple[Shape, ...]:
    side_sum_budget = (order + 6) // 4 + bound_slack
    shapes = [
        tuple(sides)
        for sides in product(range(1, side_sum_budget + 1), repeat=3)
        if sum(sides) <= side_sum_budget
    ]
    return tuple(sorted(shapes, key=lambda shape: (sum(shape), shape)))


def _finite_lattice_weights(
    shapes: tuple[Shape, ...],
    order: int,
    box_logarithm: Callable[[Shape, int], RationalSeries],
    proven_zero: Callable[[Shape], bool] | None = None,
) -> tuple[dict[Shape, RationalSeries], RationalSeries]:
    """Möbius-invert box logs, using a proved zero only for unreachable boxes."""

    weights: dict[Shape, RationalSeries] = {}
    bulk = [Fraction(0) for _ in range(order + 1)]
    for shape in shapes:
        if proven_zero is not None and proven_zero(shape):
            exact_weight = tuple(Fraction(0) for _ in range(order + 1))
        else:
            weight = list(box_logarithm(_canonical_shape(shape), order))
            for subshape, subweight in weights.items():
                if all(sub <= side for sub, side in zip(subshape, shape)):
                    placements = prod(
                        side - sub + 1 for sub, side in zip(subshape, shape)
                    )
                    for degree in range(order + 1):
                        weight[degree] -= placements * subweight[degree]
            exact_weight = tuple(weight)
        weights[shape] = exact_weight
        for degree in range(order + 1):
            bulk[degree] += exact_weight[degree]
    return weights, tuple(bulk)


def high_temperature_free_energy(
    dimension: int,
    order: int,
    *,
    bound_slack: int = 0,
) -> FLMSeries:
    """Generate the exact free-energy HT series for square or simple-cubic Ising.

    A connected even subgraph with bounding box ``s`` must traverse each
    coordinate span in both directions, so it has at least
    ``2*sum_i(s_i-1)`` edges.  ``bound_slack`` includes that many additional
    units of coordinate span and is useful for self-consistency checks.
    """

    dimension = int(dimension)
    order = _check_order(order)
    bound_slack = int(bound_slack)
    if dimension not in (2, 3):
        raise ValueError("high-temperature FLM supports dimensions 2 and 3")
    if bound_slack < 0:
        raise ValueError("bound_slack must be non-negative")

    shapes = _ht_shapes(dimension, order, bound_slack)

    def unreachable_but_proven_zero(shape: Shape) -> bool:
        canonical = _canonical_shape(shape)
        transfer_unreachable = (
            prod(canonical) > 62 or prod(canonical[:-1]) > 22
        )
        minimum_order = 2 * sum(side - 1 for side in shape)
        return transfer_unreachable and minimum_order > order

    weights, interaction = _finite_lattice_weights(
        shapes, order, _ht_box_log, unreachable_but_proven_zero
    )
    coefficients = list(interaction)
    # cosh(K) = (1-v^2)^(-1/2), hence
    # d log cosh(K) = sum_{m>=1} d/(2m) v^(2m).
    for degree in range(2, order + 1, 2):
        coefficients[degree] += Fraction(dimension, degree)
    return FLMSeries(
        dimension=dimension,
        order=order,
        variable="v",
        coefficients=tuple(coefficients),
        interaction_coefficients=interaction,
        boxes=shapes,
        box_weights=weights,
        bound_budget=order // 2 + bound_slack,
        bound_slack=bound_slack,
        order_bound="2*sum_i(side_i-1) <= truncation_order",
    )


def low_temperature_free_energy(
    dimension: int,
    order: int,
    *,
    bound_slack: int = 0,
) -> FLMSeries:
    """Generate the exact simple-cubic low-temperature free-energy series.

    For a connected flipped-spin set spanning ``a*b*c``, each coordinate-line
    projection contributes two boundary faces.  Its three projections are
    connected sets spanning rectangles ``a*b``, ``a*c``, and ``b*c`` and
    therefore contain at least ``a+b-1``, ``a+c-1``, and ``b+c-1`` cells.
    Thus the broken-bond surface is at least ``4(a+b+c)-6``.  The same bound
    applies to a connected Mayer cluster of droplets because its union is
    connected and the sum of droplet surfaces is no smaller than the union's
    exterior surface.
    """

    dimension = int(dimension)
    order = _check_order(order)
    bound_slack = int(bound_slack)
    if dimension != 3:
        raise ValueError("low-temperature FLM is implemented for dimension 3")
    if bound_slack < 0:
        raise ValueError("bound_slack must be non-negative")

    shapes = _lt_shapes(order, bound_slack)
    weights, interaction = _finite_lattice_weights(shapes, order, _lt_box_log)
    return FLMSeries(
        dimension=dimension,
        order=order,
        variable="x",
        coefficients=interaction,
        interaction_coefficients=interaction,
        boxes=shapes,
        box_weights=weights,
        bound_budget=(order + 6) // 4 + bound_slack,
        bound_slack=bound_slack,
        order_bound="4*(a+b+c)-6 <= truncation_order",
    )


def _bivariate_product(
    left: Mapping[tuple[int, int], Fraction],
    right: Mapping[tuple[int, int], Fraction],
    order: int,
) -> dict[tuple[int, int], Fraction]:
    result: dict[tuple[int, int], Fraction] = {}
    for (v_left, s_left), a in left.items():
        for (v_right, s_right), b in right.items():
            v_degree = v_left + v_right
            if v_degree > order:
                continue
            key = (v_degree, s_left + s_right)
            result[key] = result.get(key, Fraction(0)) + a * b
    return {key: value for key, value in result.items() if value}


@lru_cache(maxsize=None)
def _cosine_sum_moment(power: int) -> Fraction:
    """Exact normalized double integral of ``(cos p + cos q)**power``."""

    def cosine_moment(exponent: int) -> Fraction:
        if exponent % 2:
            return Fraction(0)
        return Fraction(comb(exponent, exponent // 2), 1 << exponent)

    return sum(
        Fraction(comb(power, j))
        * cosine_moment(j)
        * cosine_moment(power - j)
        for j in range(power + 1)
    )


def onsager_square_ht_series(order: int) -> RationalSeries:
    """Expand Onsager's square-lattice double integral exactly in ``v``.

    Independently of finite boxes,

    ``phi = log(2) - log(1-v^2) + (1/2)<log Q>``,

    where angular brackets are the normalized double integral over
    ``p,q in [0,2*pi]`` and

    ``Q = (1+v^2)^2 - 2*v*(1-v^2)*(cos p + cos q)``.

    Powers of ``cos p + cos q`` are integrated using exact central-binomial
    moments, so no numerical recognition is involved.
    """

    order = _check_order(order)
    # t = Q - 1, represented by (v degree, S degree) -> rational,
    # with S = cos(p) + cos(q).
    t = {
        (1, 1): Fraction(-2),
        (2, 0): Fraction(2),
        (3, 1): Fraction(2),
        (4, 0): Fraction(1),
    }
    power: dict[tuple[int, int], Fraction] = {(0, 0): Fraction(1)}
    logarithm: dict[tuple[int, int], Fraction] = {}
    for exponent in range(1, order + 1):
        power = _bivariate_product(power, t, order)
        sign = 1 if exponent % 2 else -1
        scale = Fraction(sign, exponent)
        for key, value in power.items():
            logarithm[key] = logarithm.get(key, Fraction(0)) + scale * value

    coefficients = [Fraction(0) for _ in range(order + 1)]
    for (v_degree, s_degree), value in logarithm.items():
        coefficients[v_degree] += Fraction(1, 2) * value * _cosine_sum_moment(
            s_degree
        )
    for degree in range(2, order + 1, 2):
        coefficients[degree] += Fraction(2, degree)  # -log(1-v^2)
    return tuple(coefficients)