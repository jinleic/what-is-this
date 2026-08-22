"""All-order structure of the interlayer expansion of a stack of Ising layers.

Run from the repository root with
    .venv/bin/python experiments/e75_interlayer_allorders.py

The structural statements are proved in ``proofs/interlayer_allorders.md``.
Everything computed here is exact integer/rational arithmetic: the generic
cumulant generator, the high-temperature moment series, the finite-lattice
inversion, and the independent even-subgraph route.  No floating point is used.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial
from pathlib import Path
from typing import Callable, Iterable, Sequence

SCRIPT = "experiments/e75_interlayer_allorders.py"
ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "interlayer" / "allorders.json"

# Orders chosen so that every exact route below finishes in seconds.
LOW_ORDER = 8       # v-order for the m=1,2 patterns (c2 and c4)
HIGH_ORDER = 6      # v-order for the m=3 patterns (c6) and the two-route check
STABILITY_ORDER = 4  # prefix re-verified with one extra in-plane span unit

GapPattern = tuple[int, ...]
MomentAtom = tuple[int, tuple[int, ...]]
Monomial = tuple[MomentAtom, ...]
Polynomial = dict[Monomial, int]
RationalSeries = tuple[Fraction, ...]
IntegerSeries = tuple[int, ...]

# Certified prefixes from the existing artifact results/interlayer/c4_series.json
# (experiments/e51_interlayer_c4.py).  They are targets here, never inputs.
STORED_C2_TOTAL = ("1/2", "0", "2", "0", "18", "0", "118", "0", "778")
STORED_C4_TOTAL_W = ("1/4", "0", "2", "0", "63", "0", "1050", "0", "14223")
STORED_C4_DIRECT_Q = ("-1/12", "0", "2/3", "0", "51", "0", "2914/3", "0", "41113/3")
STORED_C4_SAME_GAP = ("-1/12", "0", "-4/3", "0", "-17", "0", "-464/3", "0", "-3527/3")
STORED_C4_ADJACENT = ("0", "0", "2", "0", "68", "0", "1126", "0", "14880")


# --------------------------------------------------------------------------
# exact power-series helpers
# --------------------------------------------------------------------------
def _zero(order: int) -> list[Fraction]:
    return [Fraction(0) for _ in range(order + 1)]


def _one(order: int) -> RationalSeries:
    values = _zero(order)
    values[0] = Fraction(1)
    return tuple(values)


def _add_scaled(
    target: list[Fraction], source: Sequence[int | Fraction], scale: int | Fraction = 1
) -> None:
    factor = Fraction(scale)
    for degree, coefficient in enumerate(source):
        if coefficient and degree < len(target):
            target[degree] += factor * coefficient


def _mul(
    left: Sequence[int | Fraction], right: Sequence[int | Fraction], order: int
) -> RationalSeries:
    out = _zero(order)
    for left_degree, left_value in enumerate(left):
        if not left_value or left_degree > order:
            continue
        for right_degree, right_value in enumerate(right):
            total = left_degree + right_degree
            if total > order:
                break
            if right_value:
                out[total] += left_value * right_value
    return tuple(out)


def _divide(
    numerator: Sequence[int | Fraction], denominator: Sequence[int | Fraction], order: int
) -> RationalSeries:
    den = [Fraction(value) for value in denominator]
    while len(den) <= order:
        den.append(Fraction(0))
    if den[0] == 0:
        raise ZeroDivisionError("high-temperature denominators have constant term one")
    quotient = _zero(order)
    for degree in range(order + 1):
        value = Fraction(numerator[degree]) if degree < len(numerator) else Fraction(0)
        for lower in range(degree):
            value -= quotient[lower] * den[degree - lower]
        quotient[degree] = value / den[0]
    return tuple(quotient)


def _compose_q_of_w(order: int) -> RationalSeries:
    """Exact series of q=atanh(w) through ``w**order``."""

    values = _zero(order)
    for degree in range(1, order + 1, 2):
        values[degree] = Fraction(1, degree)
    return tuple(values)


def _powers_of_q(order: int) -> dict[int, RationalSeries]:
    """Exact ``[w^*] q**n`` tables for even ``n`` up to ``order``."""

    q_series = _compose_q_of_w(order)
    powers = {1: q_series}
    current = q_series
    for exponent in range(2, order + 1):
        current = _mul(current, q_series, order)
        powers[exponent] = current
    return powers


def _log_cosh_in_w(order: int) -> RationalSeries:
    """log cosh(atanh w) = -(1/2) log(1-w^2), exactly."""

    values = _zero(order)
    for half in range(1, order // 2 + 1):
        values[2 * half] = Fraction(1, 2 * half)
    return tuple(values)


# --------------------------------------------------------------------------
# generic joint-cumulant generator
# --------------------------------------------------------------------------
@lru_cache(maxsize=None)
def set_partitions(size: int) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Every set partition of ``range(size)`` exactly once."""

    if size < 0:
        raise ValueError("partition size must be non-negative")
    partitions: list[tuple[tuple[int, ...], ...]] = [()]
    for item in range(size):
        extended: list[tuple[tuple[int, ...], ...]] = []
        for partition in partitions:
            extended.append(partition + ((item,),))
            for index in range(len(partition)):
                block = partition[index] + (item,)
                extended.append(partition[:index] + (block,) + partition[index + 1 :])
        partitions = extended
    return tuple(partitions)


def mobius_coefficient(block_count: int) -> int:
    """Möbius coefficient of the partition lattice: (-1)^(k-1) (k-1)!."""

    if block_count < 1:
        raise ValueError("a nonempty cumulant partition has at least one block")
    return (-1) ** (block_count - 1) * factorial(block_count - 1)


def odd_multiplicity_tuple(values: Iterable[int]) -> tuple[int, ...]:
    """Ising parity reduction of a spin product using sigma_x^2=1."""

    counts = Counter(values)
    return tuple(sorted(value for value, count in counts.items() if count % 2))


def normalized_gap_pattern(gaps: Sequence[int]) -> GapPattern:
    """Sort and translate a gap multiset so that its minimum is zero."""

    if not gaps:
        return ()
    minimum = min(gaps)
    return tuple(sorted(int(gap) - minimum for gap in gaps))


def gap_support_connected(gaps: Sequence[int]) -> bool:
    """Whether the occupied gaps form one interval of the gap path."""

    if not gaps:
        return False
    occupied = set(int(gap) for gap in gaps)
    return occupied == set(range(min(occupied), max(occupied) + 1))


def even_gap_multiplicity(gaps: Sequence[int]) -> bool:
    """Whether every occupied gap carries an even number of bonds."""

    counts = Counter(int(gap) for gap in gaps)
    return all(multiplicity % 2 == 0 for multiplicity in counts.values())


def even_layer_incidence(gaps: Sequence[int]) -> bool:
    """Whether every layer meets an even number of the selected bonds."""

    incidence: Counter[int] = Counter()
    for gap in gaps:
        incidence[int(gap)] += 1
        incidence[int(gap) + 1] += 1
    return all(degree % 2 == 0 for degree in incidence.values())


def admissible_gap_pattern(gaps: Sequence[int]) -> bool:
    """Selection rule: interval support with every gap multiplicity even."""

    return bool(gaps) and gap_support_connected(gaps) and even_gap_multiplicity(gaps)


def compositions(total: int) -> tuple[tuple[int, ...], ...]:
    """Every ordered composition of ``total`` into positive parts."""

    if total < 0:
        raise ValueError("compositions need a non-negative total")
    if total == 0:
        return ((),)
    out: list[tuple[int, ...]] = []
    for first in range(1, total + 1):
        for rest in compositions(total - first):
            out.append((first,) + rest)
    return tuple(out)


@lru_cache(maxsize=None)
def admissible_gap_patterns(order: int) -> tuple[GapPattern, ...]:
    """All-order construction: admissible patterns are compositions of m.

    An order ``2m`` pattern occupies gaps ``0..k-1`` with even multiplicities
    ``2a_1,...,2a_k`` where ``(a_1,...,a_k)`` is a composition of ``m``.  Hence
    there are exactly ``2**(m-1)`` patterns and at most ``m+1`` layers.
    """

    if order < 1 or order % 2:
        raise ValueError("only even positive cumulant orders can be admissible")
    patterns: list[GapPattern] = []
    for parts in compositions(order // 2):
        pattern: list[int] = []
        for index, part in enumerate(parts):
            pattern.extend([index] * (2 * part))
        patterns.append(tuple(pattern))
    return tuple(sorted(patterns))


def brute_force_gap_patterns(order: int) -> tuple[GapPattern, ...]:
    """Independent exhaustive search over normalized gap multisets."""

    if order < 1:
        raise ValueError("order must be positive")
    found: set[GapPattern] = set()
    for width in range(1, order // 2 + 1):
        for gaps in product(range(width), repeat=order):
            if normalized_gap_pattern(gaps) != gaps:
                continue
            if admissible_gap_pattern(gaps):
                found.add(gaps)
    return tuple(sorted(found))


def pattern_weight(pattern: Sequence[int]) -> Fraction:
    """Ordered-placement weight ``1/prod (2 a_g)!`` of one gap pattern."""

    counts = Counter(int(gap) for gap in pattern)
    weight = Fraction(1)
    for multiplicity in counts.values():
        weight /= factorial(multiplicity)
    return weight


def block_layer_parity_ok(
    block: Sequence[int], gaps: Sequence[int], positions: Sequence[int]
) -> bool:
    """Whether every layer sees an even parity-reduced spin set in ``block``."""

    per_layer: dict[int, list[int]] = defaultdict(list)
    for bond in block:
        gap = int(gaps[bond])
        per_layer[gap].append(int(positions[bond]))
        per_layer[gap + 1].append(int(positions[bond]))
    return all(
        len(odd_multiplicity_tuple(layer_positions)) % 2 == 0
        for layer_positions in per_layer.values()
    )


def block_gap_parity_ok(block: Sequence[int], gaps: Sequence[int]) -> bool:
    """Whether ``block`` has an even number of bonds in every single gap."""

    counts = Counter(int(gaps[bond]) for bond in block)
    return all(multiplicity % 2 == 0 for multiplicity in counts.values())


def block_monomial(
    block: Sequence[int], gaps: Sequence[int], positions: Sequence[int]
) -> Monomial | None:
    """Factor one block moment of the decoupled measure into 2D moments.

    ``None`` means an odd spin count on some layer, hence exactly zero in the
    zero-field 2D Gibbs state.  Empty parity masks are identity factors.
    """

    per_layer: dict[int, list[int]] = defaultdict(list)
    for bond in block:
        gap = int(gaps[bond])
        position = int(positions[bond])
        per_layer[gap].append(position)
        per_layer[gap + 1].append(position)
    atoms: list[MomentAtom] = []
    for layer, layer_positions in per_layer.items():
        reduced = odd_multiplicity_tuple(layer_positions)
        if len(reduced) % 2:
            return None
        if reduced:
            atoms.append((layer, reduced))
    return tuple(sorted(atoms))


def cumulant_polynomial(gaps: Sequence[int], positions: Sequence[int]) -> Polynomial:
    r"""Exact moment-partition polynomial of one vertical-bond tuple.

    Returns the integer combination representing

        kappa_0(Y_1,...,Y_n)
            = sum_(pi) (-1)^(|pi|-1) (|pi|-1)! prod_(B in pi) prod_l M_l(S_(l,B)),

    with Ising parity reduction applied and equal monomials combined.
    """

    gaps = tuple(int(gap) for gap in gaps)
    positions = tuple(int(position) for position in positions)
    if len(gaps) != len(positions) or not gaps:
        raise ValueError("gaps and positions must have equal positive length")
    polynomial: defaultdict[Monomial, int] = defaultdict(int)
    for partition in set_partitions(len(gaps)):
        factors: list[MomentAtom] = []
        survives = True
        for block in partition:
            monomial = block_monomial(block, gaps, positions)
            if monomial is None:
                survives = False
                break
            factors.extend(monomial)
        if survives:
            polynomial[tuple(sorted(factors))] += mobius_coefficient(len(partition))
    return {monomial: value for monomial, value in polynomial.items() if value}


def coincidence_signature(positions: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Split a position tuple into its equality pattern and distinct sites."""

    order: dict[int, int] = {}
    canonical: list[int] = []
    sites: list[int] = []
    for position in positions:
        if position not in order:
            order[position] = len(sites)
            sites.append(int(position))
        canonical.append(order[position])
    return tuple(canonical), tuple(sites)


@lru_cache(maxsize=None)
def canonical_mask_polynomial(
    pattern: GapPattern, canonical: tuple[int, ...]
) -> tuple[tuple[tuple[tuple[int, ...], ...], int], ...]:
    """Layer-free encoding of the cumulant polynomial for one equality pattern.

    Every layer of the stack carries the same 2D Gibbs state, so a monomial is
    determined by the multiset of parity-reduced spin sets it multiplies.
    """

    polynomial = cumulant_polynomial(pattern, canonical)
    encoded: list[tuple[tuple[tuple[int, ...], ...], int]] = []
    for monomial, coefficient in polynomial.items():
        encoded.append((tuple(sorted(spins for _, spins in monomial)), coefficient))
    return tuple(sorted(encoded))


def evaluate_polynomial(polynomial: Polynomial, moment: Callable) -> Fraction:
    """Evaluate a generated polynomial with an exact moment callback."""

    total = Fraction(0)
    for monomial, coefficient in polynomial.items():
        term = Fraction(coefficient)
        for layer, positions in monomial:
            term *= Fraction(moment(layer, positions))
        total += term
    return total


def direct_joint_cumulant(
    gaps: Sequence[int], positions: Sequence[int], moment: Callable
) -> Fraction:
    """Independent literal evaluation of the defining partition formula."""

    total = Fraction(0)
    for partition in set_partitions(len(gaps)):
        term = Fraction(mobius_coefficient(len(partition)))
        for block in partition:
            per_layer: dict[int, list[int]] = defaultdict(list)
            for bond in block:
                per_layer[gaps[bond]].append(positions[bond])
                per_layer[gaps[bond] + 1].append(positions[bond])
            for layer, layer_positions in per_layer.items():
                term *= Fraction(moment(layer, odd_multiplicity_tuple(layer_positions)))
        total += term
    return total


def symbolic_key(polynomial: Polynomial) -> tuple[tuple[Monomial, int], ...]:
    return tuple(sorted(polynomial.items()))


def same_gap_c4_expected() -> Polynomial:
    """M_ijkl^2 minus the three squared pair products."""

    output: defaultdict[Monomial, int] = defaultdict(int)
    output[((0, (0, 1, 2, 3)), (1, (0, 1, 2, 3)))] += 1
    for first, second, third, fourth in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
        left = tuple(sorted((first, second)))
        right = tuple(sorted((third, fourth)))
        output[tuple(sorted(((0, left), (0, right), (1, left), (1, right))))] -= 1
    return {monomial: value for monomial, value in output.items() if value}


def adjacent_c4_expected() -> Polynomial:
    """G_ij M_ijkl G_kl minus its single surviving pair partition."""

    positive = tuple(sorted(((0, (0, 1)), (1, (0, 1, 2, 3)), (2, (2, 3)))))
    negative = tuple(sorted(((0, (0, 1)), (1, (0, 1)), (1, (2, 3)), (2, (2, 3)))))
    return {positive: 1, negative: -1}


# --------------------------------------------------------------------------
# exact two-dimensional high-temperature moments
# --------------------------------------------------------------------------
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
def parity_polynomials(width: int, height: int, order: int) -> tuple[IntegerSeries, ...]:
    """``result[mask][m]``: number of ``m``-edge subsets with odd-degree ``mask``."""

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


@lru_cache(maxsize=None)
def moment_series(width: int, height: int, order: int) -> tuple[RationalSeries, ...]:
    """Exact open-rectangle spin moments ``<prod_(x in S) sigma_x>`` in ``v``."""

    parity = parity_polynomials(width, height, order)
    partition = parity[0]
    return tuple(_divide(polynomial, partition, order) for polynomial in parity)


@lru_cache(maxsize=None)
def pattern_rectangle_sum(
    width: int, height: int, pattern: GapPattern, order: int
) -> RationalSeries:
    """Sum of one gap pattern's cumulant over all in-plane position tuples."""

    site_count = width * height
    bonds = len(pattern)
    counts: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for positions in product(range(site_count), repeat=bonds):
        canonical, sites = coincidence_signature(positions)
        for spin_sets, coefficient in canonical_mask_polynomial(pattern, canonical):
            key = []
            for spins in spin_sets:
                mask = 0
                for index in spins:
                    mask |= 1 << sites[index]
                key.append(mask)
            counts[tuple(sorted(key))] += coefficient
    moments = moment_series(width, height, order)
    total = _zero(order)
    for key, coefficient in counts.items():
        if not coefficient:
            continue
        product_series = _one(order)
        for mask in key:
            product_series = _mul(product_series, moments[mask], order)
            if not any(product_series):
                break
        _add_scaled(total, product_series, coefficient)
    return tuple(total)


def plane_shapes(order: int, bound_slack: int = 0) -> tuple[tuple[int, int], ...]:
    """Rectangles needed for exact finite-lattice inversion through ``v**order``."""

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


def pattern_bulk_series(
    patterns: Sequence[GapPattern], order: int, bound_slack: int = 0
) -> dict[GapPattern, RationalSeries]:
    """Rectangular finite-lattice inversion of every pattern's position sum."""

    shapes = plane_shapes(order, bound_slack)
    weights: dict[tuple[int, int], tuple[RationalSeries, ...]] = {}
    bulk = [_zero(order) for _ in patterns]
    for width, height in shapes:
        values = [
            list(pattern_rectangle_sum(width, height, pattern, order))
            for pattern in patterns
        ]
        for (sub_width, sub_height), subweight in weights.items():
            if sub_width <= width and sub_height <= height:
                placements = (width - sub_width + 1) * (height - sub_height + 1)
                for index in range(len(patterns)):
                    _add_scaled(values[index], subweight[index], -placements)
        exact = tuple(tuple(value) for value in values)
        weights[(width, height)] = exact
        for index in range(len(patterns)):
            _add_scaled(bulk[index], exact[index])
    return {pattern: tuple(bulk[index]) for index, pattern in enumerate(patterns)}


def cumulant_coefficient(
    order_2m: int, v_order: int, bound_slack: int = 0
) -> tuple[RationalSeries, dict[GapPattern, RationalSeries]]:
    """Exact ``c_(2m)(K)`` in the direct coupling ``q`` and its pattern parts."""

    patterns = admissible_gap_patterns(order_2m)
    bulk = pattern_bulk_series(patterns, v_order, bound_slack)
    total = _zero(v_order)
    contributions: dict[GapPattern, RationalSeries] = {}
    for pattern in patterns:
        weight = pattern_weight(pattern)
        scaled = tuple(weight * value for value in bulk[pattern])
        contributions[pattern] = scaled
        _add_scaled(total, scaled)
    return tuple(total), contributions


def to_tanh_variable(
    coefficients: dict[int, RationalSeries], order: int, v_order: int
) -> RationalSeries:
    """Convert ``sum_n c_n q^n`` into the ``[w^order]`` coefficient, w=tanh(q)."""

    powers = _powers_of_q(order)
    out = _zero(v_order)
    for exponent, series_values in coefficients.items():
        factor = powers[exponent][order]
        if factor:
            _add_scaled(out, series_values, factor)
    return tuple(out)


# --------------------------------------------------------------------------
# independent even-subgraph route
# --------------------------------------------------------------------------
@lru_cache(maxsize=None)
def slab_even_subgraph_columns(
    width: int, height: int, layers: int, v_order: int, w_order: int
) -> tuple[RationalSeries, ...]:
    """Exact ``[w^d] P(v,w)`` columns of an open ``width x height x layers`` slab.

    Vertical bond subsets are labelled by one parity mask per layer gap; layer
    ``i`` then carries the mask ``S_i xor S_(i+1)`` with ``S_0=S_layers=0``.
    """

    if layers < 1:
        raise ValueError("a slab needs at least one layer")
    parity = parity_polynomials(width, height, v_order)
    site_count = width * height
    masks_by_weight: list[list[int]] = [[] for _ in range(w_order + 1)]
    for mask in range(1 << site_count):
        weight = bin(mask).count("1")
        if weight <= w_order:
            masks_by_weight[weight].append(mask)
    columns = [_zero(v_order) for _ in range(w_order + 1)]
    gaps = layers - 1
    if gaps == 0:
        columns[0][:] = [Fraction(value) for value in parity[0]]
        return tuple(tuple(column) for column in columns)

    def walk(index: int, used: int, chosen: tuple[int, ...]) -> None:
        if index == gaps:
            series = _one(v_order)
            boundaries = (0,) + chosen + (0,)
            for layer in range(layers):
                mask = boundaries[layer] ^ boundaries[layer + 1]
                series = _mul(series, parity[mask], v_order)
                if not any(series):
                    return
            _add_scaled(columns[used], series)
            return
        for weight in range(w_order + 1 - used):
            for mask in masks_by_weight[weight]:
                walk(index + 1, used + weight, chosen + (mask,))

    walk(0, 0, ())
    return tuple(tuple(column) for column in columns)


def slab_log_columns(
    width: int, height: int, layers: int, v_order: int, w_order: int
) -> dict[int, RationalSeries]:
    """Exact ``[w^d] log P`` columns of one slab for every ``d<=w_order``.

    Uses the exact Newton recursion ``d l_d = d r_d - sum_(j<d) j l_j r_(d-j)``
    for ``P/P_0=exp(L)``, so no order-by-order hand expansion is needed.
    """

    columns = slab_even_subgraph_columns(width, height, layers, v_order, w_order)
    ratios = {
        degree: _divide(columns[degree], columns[0], v_order)
        for degree in range(1, w_order + 1)
    }
    out: dict[int, RationalSeries] = {}
    for degree in range(1, w_order + 1):
        accumulated = _zero(v_order)
        _add_scaled(accumulated, ratios[degree], degree)
        for lower in range(1, degree):
            _add_scaled(
                accumulated,
                _mul(out[lower], ratios[degree - lower], v_order),
                -lower,
            )
        out[degree] = tuple(value / degree for value in accumulated)
    out[-1] = tuple(
        Fraction(int(any(columns[degree]))) for degree in range(1, w_order + 1, 2)
    )
    return out


def anisotropic_bulk(
    v_order: int, w_order: int, bound_slack: int = 0, max_layers: int | None = None
) -> tuple[dict[int, RationalSeries], dict[int, dict[int, RationalSeries]], bool]:
    """Independent 3D finite-lattice inversion in ``v`` and ``w=tanh(K_z)``."""

    layer_budget = max_layers if max_layers is not None else w_order // 2 + 1
    shapes = tuple(
        sorted(
            (
                (width, height, layers)
                for width, height in plane_shapes(v_order, bound_slack)
                for layers in range(1, layer_budget + 1)
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )
    degrees = tuple(range(2, w_order + 1, 2))
    weights: dict[tuple[int, int, int], dict[int, RationalSeries]] = {}
    bulk = {degree: _zero(v_order) for degree in degrees}
    by_layers = {
        layers: {degree: _zero(v_order) for degree in degrees}
        for layers in range(1, layer_budget + 1)
    }
    odd_zero = True
    for width, height, layers in shapes:
        columns = slab_log_columns(width, height, layers, v_order, w_order)
        odd_zero = odd_zero and not any(columns[-1])
        values = {degree: list(columns[degree]) for degree in degrees}
        for (sub_width, sub_height, sub_layers), subweight in weights.items():
            if sub_width <= width and sub_height <= height and sub_layers <= layers:
                placements = (
                    (width - sub_width + 1)
                    * (height - sub_height + 1)
                    * (layers - sub_layers + 1)
                )
                for degree in degrees:
                    _add_scaled(values[degree], subweight[degree], -placements)
        exact = {degree: tuple(values[degree]) for degree in degrees}
        weights[(width, height, layers)] = exact
        for degree in degrees:
            _add_scaled(bulk[degree], exact[degree])
            _add_scaled(by_layers[layers][degree], exact[degree])
    return (
        {degree: tuple(bulk[degree]) for degree in degrees},
        {
            layers: {degree: tuple(values[degree]) for degree in values}
            for layers, values in by_layers.items()
        },
        odd_zero,
    )


# --------------------------------------------------------------------------
# finite-volume symmetry controls
# --------------------------------------------------------------------------
def slab_partition_in_q(
    horizontal_edges: Sequence[tuple[int, int]],
    sites_per_layer: int,
    layers: int,
    horizontal_weight: Fraction,
    max_degree: int,
    vertical_periodic: bool = False,
) -> dict[int, Fraction]:
    """Exact ``q``-Taylor coefficients of a finite open-layer slab ``Z``.

    Horizontal factors are evaluated exactly at the rational weight
    ``a=exp(K)`` raised to integer energies; the vertical coupling stays
    formal through ``exp(q E)=sum q^n E^n/n!``.
    """

    if layers < 1 or sites_per_layer < 1:
        raise ValueError("positive slab dimensions required")
    spin_count = layers * sites_per_layer
    vertical_pairs = [
        (layer * sites_per_layer + site, ((layer + 1) % layers) * sites_per_layer + site)
        for layer in range(layers if vertical_periodic else layers - 1)
        for site in range(sites_per_layer)
    ]
    coefficients = {degree: Fraction(0) for degree in range(max_degree + 1)}
    for state in range(1 << spin_count):
        spins = tuple(1 if (state >> site) & 1 else -1 for site in range(spin_count))
        horizontal_energy = 0
        for layer in range(layers):
            offset = layer * sites_per_layer
            horizontal_energy += sum(
                spins[offset + left] * spins[offset + right]
                for left, right in horizontal_edges
            )
        weight = horizontal_weight**horizontal_energy
        vertical_energy = sum(spins[left] * spins[right] for left, right in vertical_pairs)
        term = Fraction(1)
        for degree in range(max_degree + 1):
            if degree:
                term *= Fraction(vertical_energy, degree)
            coefficients[degree] += weight * term
    return coefficients


def _record(checks: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    checks.append({"name": name, "passed": bool(passed), "detail": detail})
    print(f"{'PASS' if passed else 'FAIL'} {name}: {detail}")


def _strings(values: Iterable[int | Fraction]) -> list[str]:
    return [str(Fraction(value)) for value in values]


def main() -> None:
    checks: list[dict[str, object]] = []
    started = time.time()

    # ---------------- structural combinatorics ----------------
    bells = {size: len(set_partitions(size)) for size in (2, 4, 6)}
    _record(
        checks,
        "set-partition counts m=1,2,3",
        bells == {2: 2, 4: 15, 6: 203},
        f"exact Bell numbers B(2m) are {bells}",
    )

    pattern_table = {m: admissible_gap_patterns(2 * m) for m in range(1, 8)}
    brute_force_match = all(
        pattern_table[m] == brute_force_gap_patterns(2 * m) for m in (1, 2, 3, 4)
    )
    count_match = all(len(pattern_table[m]) == 2 ** (m - 1) for m in range(1, 8))
    layer_match = all(
        max(max(pattern) + 2 for pattern in admissible_gap_patterns(2 * m)) == m + 1
        for m in range(1, 8)
    )
    _record(
        checks,
        "layer-gap selection rule m=1..4",
        brute_force_match and count_match and layer_match,
        "composition construction equals exhaustive search through m=4; counts are 2^(m-1) and the maximal vertical extent is m+1 layers through m=7",
    )

    parity_equivalence = True
    partition_survivors: dict[str, int] = {}
    for m in (1, 2, 3):
        positions = tuple(range(2 * m))
        for pattern in pattern_table[m]:
            layer_ok = 0
            for partition in set_partitions(2 * m):
                by_layer = all(
                    block_layer_parity_ok(block, pattern, positions)
                    for block in partition
                )
                by_gap = all(block_gap_parity_ok(block, pattern) for block in partition)
                if by_layer != by_gap:
                    parity_equivalence = False
                layer_ok += int(by_layer)
            partition_survivors[",".join(map(str, pattern))] = layer_ok
    _record(
        checks,
        "block parity lemma m=1,2,3",
        parity_equivalence,
        f"even layer parity per block equals even per-gap multiplicity per block on all patterns; surviving partition counts {partition_survivors}",
    )

    # ---------------- generic generator against the certified c4 lemma -------
    c2_polynomial = cumulant_polynomial((0, 0), (0, 1))
    same_gap = cumulant_polynomial((0, 0, 0, 0), (0, 1, 2, 3))
    adjacent = cumulant_polynomial((0, 0, 1, 1), (0, 1, 2, 3))
    formulas_pass = (
        symbolic_key(c2_polynomial) == symbolic_key({((0, (0, 1)), (1, (0, 1))): 1})
        and symbolic_key(same_gap) == symbolic_key(same_gap_c4_expected())
        and symbolic_key(adjacent) == symbolic_key(adjacent_c4_expected())
    )
    _record(
        checks,
        "generic generator reproduces c2 and c4 formulas",
        formulas_pass,
        "kappa_2=G_ij^2, same-gap kappa_4=M_ijkl^2-sum of three squared pair products, adjacent kappa_4=G_ij G_kl M_ijkl-G_ij^2 G_kl^2",
    )

    def oracle(layer: int, positions: tuple[int, ...]) -> Fraction:
        if len(positions) % 2:
            return Fraction(0)
        if not positions:
            return Fraction(1)
        code = sum((index + 3) * (position + 2) for index, position in enumerate(positions))
        return Fraction((layer + 2) * (code + 5), code + 7)

    generator_pass = True
    monomial_counts: dict[str, int] = {}
    for m in (1, 2, 3):
        positions = tuple(range(2 * m))
        for pattern in pattern_table[m]:
            polynomial = cumulant_polynomial(pattern, positions)
            monomial_counts[",".join(map(str, pattern))] = len(polynomial)
            if evaluate_polynomial(polynomial, oracle) != direct_joint_cumulant(
                pattern, positions, oracle
            ):
                generator_pass = False
    _record(
        checks,
        "generic generator versus literal definition m=1,2,3",
        generator_pass,
        f"all seven patterns agree exactly on a non-Gaussian rational moment oracle; simplified monomial counts {monomial_counts}",
    )

    inadmissible = (
        (0, 0, 0, 1),
        (0, 1, 1, 2),
        (0, 0, 2, 2),
        (0, 2, 2, 4),
        (0, 0, 0, 0, 2, 2),
        (0, 0, 1, 1, 3, 3),
        (0, 0, 0, 1, 1, 1),
    )
    vanishing_pass = all(
        cumulant_polynomial(pattern, tuple(range(len(pattern)))) == {}
        and direct_joint_cumulant(pattern, tuple(range(len(pattern))), oracle) == 0
        for pattern in inadmissible
    )
    _record(
        checks,
        "inadmissible patterns vanish formally",
        vanishing_pass,
        "odd gap multiplicity and disconnected gap support both give the empty moment polynomial, so the cancellation is formal and not numerical",
    )

    # ---------------- exact series: decoupled control ----------------
    log_cosh_q = {2: Fraction(1, 2), 4: Fraction(-1, 12), 6: Fraction(1, 45)}
    zero_coupling = {}
    for m in (1, 2, 3):
        total, _ = cumulant_coefficient(2 * m, 0)
        zero_coupling[2 * m] = total[0]
    decoupled_pass = zero_coupling == log_cosh_q
    _record(
        checks,
        "K=0 log cosh control",
        decoupled_pass,
        f"c2,c4,c6 at v=0 are {[str(zero_coupling[degree]) for degree in (2, 4, 6)]}, exactly [q^n] log cosh q",
    )

    # ---------------- exact series: c2 and c4 ----------------
    c2_total, c2_parts = cumulant_coefficient(2, LOW_ORDER)
    c4_total, c4_parts = cumulant_coefficient(4, LOW_ORDER)
    c2_pass = _strings(c2_total) == list(STORED_C2_TOTAL)
    c4_same = c4_parts[(0, 0, 0, 0)]
    c4_adjacent = c4_parts[(0, 0, 1, 1)]
    c4_pass = (
        _strings(c4_total) == list(STORED_C4_DIRECT_Q)
        and _strings(c4_same) == list(STORED_C4_SAME_GAP)
        and _strings(c4_adjacent) == list(STORED_C4_ADJACENT)
    )
    _record(
        checks,
        "stored c2 prefix reproduced",
        c2_pass,
        f"generic route gives c2={_strings(c2_total)} through v^{LOW_ORDER}, including [v^8]=778",
    )
    _record(
        checks,
        "stored c4 prefix and gap split reproduced",
        c4_pass,
        f"generic route reproduces c4, its same-gap part and the adjacent-gap correction through v^{LOW_ORDER}",
    )

    c4_total_w = to_tanh_variable({2: c2_total, 4: c4_total}, 4, LOW_ORDER)
    c2_total_w = to_tanh_variable({2: c2_total}, 2, LOW_ORDER)
    w_pass = _strings(c4_total_w) == list(STORED_C4_TOTAL_W) and _strings(
        c2_total_w
    ) == list(STORED_C2_TOTAL)
    _record(
        checks,
        "exact q to w change of variable",
        w_pass,
        f"[w^4]=c4+(2/3)c2 reproduces the stored tanh-variable prefix {list(STORED_C4_TOTAL_W)}",
    )

    # ---------------- exact series: c6, derived independently ----------------
    c6_total, c6_parts = cumulant_coefficient(6, HIGH_ORDER)
    c2_high, _ = cumulant_coefficient(2, HIGH_ORDER)
    c4_high, _ = cumulant_coefficient(4, HIGH_ORDER)
    c6_total_w = to_tanh_variable({2: c2_high, 4: c4_high, 6: c6_total}, 6, HIGH_ORDER)
    log_cosh_w = _log_cosh_in_w(6)
    c6_residual_w = tuple(
        c6_total_w[degree] - (log_cosh_w[6] if degree == 0 else 0)
        for degree in range(HIGH_ORDER + 1)
    )
    c6_v0_pass = c6_total[0] == Fraction(1, 45) and c6_total_w[0] == Fraction(1, 6)
    _record(
        checks,
        "c6 decoupled normalization",
        c6_v0_pass,
        "[v^0] c6 is 1/45 in q and 1/6 in w, matching log cosh q and -(1/2)log(1-w^2)",
    )

    nonvanishing = {
        ",".join(map(str, pattern)): any(values)
        for pattern, values in c6_parts.items()
    }
    _record(
        checks,
        "all four c6 gap patterns contribute",
        all(nonvanishing.values()),
        f"every admissible m=3 pattern has a nonzero exact series through v^{HIGH_ORDER}",
    )

    # ---------------- independent even-subgraph route ----------------
    aniso_bulk, aniso_by_layers, aniso_odd_zero = anisotropic_bulk(HIGH_ORDER, 6)
    route_c2 = tuple(
        aniso_bulk[2][degree] + (Fraction(1, 2) if degree == 0 else 0)
        for degree in range(HIGH_ORDER + 1)
    )
    route_c4_w = tuple(
        aniso_bulk[4][degree] + (Fraction(1, 4) if degree == 0 else 0)
        for degree in range(HIGH_ORDER + 1)
    )
    route_c6_w = tuple(
        aniso_bulk[6][degree] + (Fraction(1, 6) if degree == 0 else 0)
        for degree in range(HIGH_ORDER + 1)
    )
    c2_w_high = to_tanh_variable({2: c2_high}, 2, HIGH_ORDER)
    c4_w_high = to_tanh_variable({2: c2_high, 4: c4_high}, 4, HIGH_ORDER)
    two_route_pass = (
        route_c2 == c2_w_high
        and route_c4_w == c4_w_high
        and route_c6_w == c6_total_w
    )
    _record(
        checks,
        "two-route agreement through w^6",
        two_route_pass,
        f"independent even-subgraph slab inversion reproduces [w^2], [w^4] and [w^6] coefficientwise through v^{HIGH_ORDER}",
    )
    _record(
        checks,
        "odd vertical columns vanish",
        aniso_odd_zero,
        "every enumerated w^1, w^3 and w^5 slab column of every box is identically zero",
    )

    max_extent_pattern = (0, 0, 1, 1, 2, 2)
    height4_pass = aniso_by_layers[4][6] == c6_parts[max_extent_pattern]
    height1_pass = all(
        not any(aniso_by_layers[1][degree]) for degree in (2, 4, 6)
    )
    _record(
        checks,
        "vertical extent four isolates the maximal c6 pattern",
        height4_pass and height1_pass,
        "the four-layer finite-lattice weight at w^6 equals the (0,0,1,1,2,2) cumulant contribution exactly, and single-layer weights vanish",
    )

    reflection_pass = True
    for m in (1, 2, 3):
        _, parts = (
            (c2_total, c2_parts)
            if m == 1
            else (c4_total, c4_parts)
            if m == 2
            else (c6_total, c6_parts)
        )
        for pattern, values in parts.items():
            counts = Counter(pattern)
            reversed_parts = tuple(
                counts[gap] for gap in reversed(range(max(counts) + 1))
            )
            mirrored: list[int] = []
            for index, multiplicity in enumerate(reversed_parts):
                mirrored.extend([index] * multiplicity)
            if parts[tuple(mirrored)] != values:
                reflection_pass = False
    _record(
        checks,
        "layer-reflection pairing of gap patterns",
        reflection_pass,
        "every pattern contribution equals the contribution of the reversed composition, exactly, for m=1,2,3",
    )

    # ---------------- all-order verifications at second in-plane order ------
    v2_bulk, v2_by_layers, v2_odd_zero = anisotropic_bulk(2, 12)
    v2_identity = {
        m: v2_bulk[2 * m][2] for m in range(1, 7)
    }
    v2_pass = all(value == 2 for value in v2_identity.values()) and v2_odd_zero
    _record(
        checks,
        "all-order [v^2 w^(2m)] ladder identity",
        v2_pass,
        f"the exact residual second-order in-plane coefficient is 2 for every m=1..6: {[str(v2_identity[m]) for m in range(1, 7)]}",
    )

    extent_pass = True
    extent_witness: dict[str, str] = {}
    for m in range(1, 7):
        for layers, columns in v2_by_layers.items():
            value = columns[2 * m][2]
            if layers > m + 1 and value != 0:
                extent_pass = False
            if layers == m + 1:
                extent_witness[str(m)] = str(value)
    _record(
        checks,
        "vertical extent bound m+1 verified to m=6",
        extent_pass,
        f"independent even-subgraph weights at w^(2m) vanish for every slab taller than m+1 layers; extremal weights {extent_witness}",
    )

    v2_slack_bulk, _, _ = anisotropic_bulk(2, 8, 1)
    v2_slack_pass = all(
        v2_slack_bulk[2 * m][2] == v2_bulk[2 * m][2] for m in range(1, 5)
    )
    _record(
        checks,
        "all-order identity bound stability",
        v2_slack_pass,
        "one extra in-plane span unit leaves the [v^2 w^(2m)] coefficients unchanged for m=1..4",
    )

    # ---------------- finite-lattice bound stability ----------------
    stability_pass = True
    for order_2m in (2, 4, 6):
        tight, _ = cumulant_coefficient(order_2m, STABILITY_ORDER, 0)
        slack, _ = cumulant_coefficient(order_2m, STABILITY_ORDER, 1)
        if tight != slack:
            stability_pass = False
    _record(
        checks,
        "finite-lattice bound stability",
        stability_pass,
        f"adding every rectangle with one extra in-plane span unit changes no c2, c4 or c6 coefficient through v^{STABILITY_ORDER}",
    )

    # ---------------- finite-volume evenness controls ----------------
    bipartite_slab = slab_partition_in_q(
        horizontal_edges=((0, 1), (1, 2), (0, 2)),
        sites_per_layer=3,
        layers=4,
        horizontal_weight=Fraction(3, 2),
        max_degree=9,
    )
    even_pass = all(bipartite_slab[degree] == 0 for degree in range(1, 10, 2))
    positive_pass = all(bipartite_slab[degree] > 0 for degree in (0, 2))
    _record(
        checks,
        "finite slab interlayer evenness",
        even_pass and positive_pass,
        "a four-layer slab over a frustrated in-plane triangle has exactly zero q^1,q^3,q^5,q^7,q^9 coefficients, so in-plane bipartiteness is not needed",
    )

    odd_ring = slab_partition_in_q(
        horizontal_edges=(),
        sites_per_layer=1,
        layers=3,
        horizontal_weight=Fraction(1),
        max_degree=3,
        vertical_periodic=True,
    )
    scope_pass = odd_ring[1] == 0 and odd_ring[3] != 0
    _record(
        checks,
        "odd periodic layer ring breaks evenness",
        scope_pass,
        f"three layers glued periodically give [q^3]Z={odd_ring[3]} != 0, so the alternating flip needs a layer-parity two-colouring",
    )

    elapsed = round(time.time() - started, 3)
    payload = {
        "provenance": {
            "script": SCRIPT,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "interpreter": ".venv/bin/python",
            "actual_executable": sys.executable,
            "python": platform.python_version(),
            "arithmetic": "Python int and Fraction only",
            "elapsed_seconds": elapsed,
            "method": "generic joint-cumulant set-partition generator over exact 2D high-temperature moments, rectangular finite-lattice inversion, independent even-subgraph slab inversion, and exact finite-slab spin enumeration",
            "independence": "c6 was derived here from the generic cumulant generator and, independently, from the even-subgraph slab route; a third full spin-enumeration route runs in tests/test_interlayer_allorders.py. No external c6 artifact was read. A blind coefficient exchange with the concurrent agent InterlayerC6 reported exact agreement on every c6 coefficient and on the pattern split.",
        },
        "data": {
            "claim": "[THEOREM] For any zero-field Ising stack whose interlayer bonds respect the layer-parity two-colouring, every finite-volume pressure is even in the interlayer coupling, so all odd coefficients vanish; the even coefficient c_(2m) is a finite sum over exactly 2^(m-1) layer-gap patterns of set-partition cumulants of 2D moments, involving at most m+1 layers.",
            "variables": {
                "v": "tanh(K), the in-plane high-temperature variable",
                "q": "the direct interlayer coupling K_z",
                "w": "tanh(K_z), the repository interlayer convention",
                "Y": "Y_(g,x)=sigma_(x,g) sigma_(x,g+1), the vertical bond variable of gap g at in-plane position x",
                "pressure": "f(K,q)=f_2D(K)+sum_(n>=1) c_n(K) q^n",
            },
            "theorem_all_odd_vanishing": {
                "tag": "[THEOREM]",
                "hypotheses": "zero external field; every interlayer bond joins layers of opposite parity (nearest-neighbour stacking, free vertical boundary, or periodic stacking with an even number of layers); arbitrary in-plane graph and arbitrary layer-dependent in-plane couplings",
                "transformation": "sigma_(x,l) -> (-1)^l sigma_(x,l)",
                "effect": "every in-plane product sigma_(x,l) sigma_(y,l) is fixed and every interlayer product sigma_(x,l) sigma_(x,l+1) changes sign",
                "finite_volume": "Z_Lambda(K,q)=Z_Lambda(K,-q) identically, hence f_Lambda(K,q)=(1/|Lambda|) log Z_Lambda is even and entire in q, so every odd Taylor coefficient of f_Lambda vanishes exactly",
                "w_variable": "w=tanh(q) is odd, so evenness and all-odd vanishing hold verbatim in w",
                "thermodynamic_limit": "any pointwise limit of the even functions f_Lambda is even in q; where the limit is differentiable all odd derivatives vanish, and where it is analytic every odd Taylor coefficient vanishes. The evenness statement needs no analyticity; the coefficient statement is exactly as strong as the regularity assumed for the limit.",
                "necessity_of_parity_colouring": "for three layers glued periodically the flip is unavailable and the exact q^3 coefficient of Z is nonzero, so the hypothesis cannot be dropped",
                "in_plane_bipartiteness_not_needed": "verified exactly on a four-layer slab over an in-plane triangle",
            },
            "theorem_even_cumulant_formula": {
                "tag": "[THEOREM]",
                "finite_volume": "c_n=(1/(n! |Lambda|)) sum over all n-tuples of vertical bonds of kappa_0(Y_(b_1),...,Y_(b_n)), where kappa_0 is the joint cumulant in the decoupled product of independent 2D layers",
                "partition_formula": "kappa_0(Y_1,...,Y_n)=sum_(pi in Pi_n) (-1)^(|pi|-1)(|pi|-1)! prod_(B in pi) prod_l M_l(S_(l,B))",
                "moment_factor": "S_(l,B) is the parity-reduced set of in-plane positions that bonds of B place on layer l; M_l(S)=<prod_(x in S) sigma_x>_2D and M_l(S)=0 whenever |S| is odd",
                "bulk_form": "c_(2m)=sum over admissible normalized gap patterns P of (1/prod_g m_g(P)!) sum_(x_2,...,x_(2m)) kappa_0(Y_(g_1,0),Y_(g_2,x_2),...,Y_(g_(2m),x_(2m))), with the pattern's gaps in non-decreasing order; the in-plane sum is per site by translation invariance",
                "convergence_hypothesis": "the bulk form requires absolute convergence of the in-plane cumulant sum; each individual coefficient computed here is an exact finite-lattice-inverted series prefix and needs no such hypothesis",
                "algorithm": [
                    "enumerate the 2^(m-1) admissible gap patterns as compositions of m",
                    "for each pattern enumerate set partitions of the 2m bond indices",
                    "reject a partition as soon as one block has odd bond multiplicity in some gap",
                    "factor every surviving block into one 2D moment per layer, parity-reducing repeated positions",
                    "weight the product by (-1)^(|pi|-1)(|pi|-1)! and combine equal moment monomials over the integers",
                    "sum over in-plane position tuples and apply rectangular finite-lattice inversion",
                ],
            },
            "theorem_layer_gap_selection": {
                "tag": "[THEOREM]",
                "statement": "kappa_0(Y_(g_1,x_1),...,Y_(g_n,x_n)) vanishes identically unless every occupied gap carries an even number of bonds and the occupied gaps form one interval. Consequently the admissible order-2m patterns are exactly the compositions (a_1,...,a_k) of m with 2a_i bonds in gap i-1, there are 2^(m-1) of them, and at most m+1 consecutive layers participate.",
                "proof_step_block_parity": "inside a surviving partition, the topmost occupied gap forces each block to have an even count there; layer by layer downward, equal parity of the two gaps meeting a layer forces every block to have an even count in every gap",
                "proof_step_interval": "an unoccupied gap splits the bonds into two families supported on disjoint layer sets, which are independent under the decoupled measure, and mixed cumulants of independent families vanish; the cancellation is formal, so the generator returns the empty polynomial",
                "direction": "the two conditions are proved necessary; they are not proved sufficient in general, but every admissible pattern with m<=3 is verified here to have a nonzero exact series",
                "pattern_counts_m1_to_m7": {
                    str(m): 2 ** (m - 1) for m in range(1, 8)
                },
                "patterns_m1_m2_m3": {
                    str(m): [list(pattern) for pattern in pattern_table[m]]
                    for m in (1, 2, 3)
                },
                "surviving_partition_counts_m1_m2_m3": partition_survivors,
                "maximal_vertical_extent": "m+1 layers, attained only by the pattern with all a_i=1",
            },
            "theorem_adjacent_gap_correction": {
                "tag": "[THEOREM]",
                "origin": "the adjacent-gap term of c4 is the composition (1,1) of m=2, the unique admissible pattern beyond the single-gap composition (2)",
                "same_gap_kappa4": "M_ijkl^2-G_ij^2 G_kl^2-G_ik^2 G_jl^2-G_il^2 G_jk^2",
                "adjacent_gap_kappa4": "G_ij G_kl M_ijkl-G_ij^2 G_kl^2",
                "weights": "1/4! for (2) and 1/(2!2!)=1/4 for (1,1)",
                "generalization": "at order 2m the multi-gap patterns are the 2^(m-1)-1 compositions with k>=2 parts, weight 1/prod (2a_i)!, and the k=m composition needs exactly m+1 layers; the c4 adjacent term is the first instance, not a boundary artefact",
                "certified_identity": "the four-layer finite-lattice weight at w^6 equals the (0,0,1,1,2,2) cumulant contribution exactly, confirming the extent bound at m=3",
            },
            "theorem_reflection_pairing": {
                "tag": "[THEOREM]",
                "statement": "The layer reflection l -> L-1-l maps the gap pattern of composition (a_1,...,a_k) to the reversed composition (a_k,...,a_1) and preserves the decoupled layer measure, so the two pattern contributions to c_(2m) are identical.",
                "consequence": "admissible patterns come in reflection pairs, self-paired exactly when the composition is a palindrome; the c6 contributions of (0,0,0,0,1,1) and (0,0,1,1,1,1) coincide term by term",
                "verified": "exactly for m=1,2,3 on the computed series",
            },
            "theorem_ladder_identity": {
                "tag": "[THEOREM]",
                "statement": "For the simple-cubic stack with in-plane variables v_x,v_y and interlayer variable w, the second in-plane order of the residual pressure is [v^2 w^(2m)](f-f_2D-log cosh K_z)=2 for every m>=1, and in the anisotropic case the coefficient of v_x^2 w^(2m) and of v_y^2 w^(2m) is 1 each.",
                "proof": "an even subgraph with exactly two in-plane edges is a single cycle alternating two vertical straight paths with two copies of one in-plane edge; the two vertical paths have equal length m, so per site the cycles are indexed by the in-plane edge direction only, and the logarithm has no multi-cycle correction at this in-plane order because every cycle in an open stack uses at least two in-plane edges",
                "extent_corollary": "each such cycle spans exactly m+1 layers, so the whole coefficient sits in the (m+1)-layer finite-lattice weight, saturating the extent bound",
                "verified_m": [str(v2_identity[m]) for m in range(1, 7)],
                "extremal_layer_weights": extent_witness,
            },
            "theorem_extent_bound_check": {
                "tag": "[COMPUTATION]",
                "statement": "Independent even-subgraph finite-lattice weights at w^(2m) vanish for every slab of more than m+1 layers, for all m=1..6 at second in-plane order.",
                "role": "route-independent confirmation of the m+1 vertical-extent consequence of the selection theorem",
            },
            "exact_series": {
                "tag": "[COMPUTATION]",
                "v_order_low": LOW_ORDER,
                "v_order_high": HIGH_ORDER,
                "c2_direct_q": _strings(c2_total),
                "c4_direct_q": _strings(c4_total),
                "c4_same_gap_component": _strings(c4_same),
                "c4_adjacent_gap_component": _strings(c4_adjacent),
                "c4_total_w": _strings(c4_total_w),
                "c6_direct_q": _strings(c6_total),
                "c6_total_w": _strings(c6_total_w),
                "c6_residual_w_after_log_cosh": _strings(c6_residual_w),
                "c6_pattern_components_direct_q": {
                    ",".join(map(str, pattern)): _strings(values)
                    for pattern, values in sorted(c6_parts.items())
                },
                "independent_even_subgraph_route": {
                    "c2_w": _strings(route_c2),
                    "c4_w": _strings(route_c4_w),
                    "c6_w": _strings(route_c6_w),
                },
                "note": "c6 is new here and was obtained twice independently: once from the generic cumulant generator over 2D moments, once from four-layer even-subgraph slab inversion.",
            },
            "structural_counts": {
                "bell_numbers": {str(size): bells[size] for size in (2, 4, 6)},
                "simplified_monomial_counts": monomial_counts,
                "finite_slab_odd_coefficients": _strings(
                    [bipartite_slab[degree] for degree in range(1, 10, 2)]
                ),
                "odd_periodic_ring_q3": str(odd_ring[3]),
                "ladder_identity_v2_coefficients": {
                    str(m): str(v2_identity[m]) for m in range(1, 7)
                },
                "residual_w_columns_at_v_order_two": {
                    str(2 * m): _strings(v2_bulk[2 * m]) for m in range(1, 7)
                },
            },
            "limitations": {
                "tag": "[UNRESOLVED]",
                "sign": "No sign or positivity theorem is claimed for c_(2m) or for its pattern components; the Mobius weights alternate and the computed c4 same-gap component is negative while the adjacent component is positive.",
                "radius": "No radius of convergence is claimed, in q, in w, or in v. Nothing here bounds the interlayer expansion at or near the isotropic critical point.",
                "moments": "The formula reduces c_(2m) to 2D multipoint moments; it does not evaluate those moments in closed form, so the 3D model is not solved.",
                "resummation": "The number of admissible patterns grows as 2^(m-1) and the number of surviving partitions grows superexponentially, so the theorem is structural rather than a closed-form solution.",
            },
        },
        "checks": checks,
    }
    if not all(check["passed"] for check in checks):
        raise AssertionError("at least one exact all-orders check failed")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {RESULT.relative_to(ROOT)} in {elapsed}s")
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # one clear failure line for artifact scripts
        print(f"FAIL {exc}")
        raise
