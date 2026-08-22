"""Exact interlayer-coupling expansion tools for stacked square-lattice Ising models.

The expansion variable is ``w = tanh(K_z)``.  Two coefficient conventions are
kept distinct throughout this module:

``residual_c2``
    ``[w**2] (phi_3D - phi_2D - log(cosh(K_z)))``;

``total_c2``
    ``[w**2] (phi_3D - phi_2D) = residual_c2 + 1/2``.

At zero field, for a translation-invariant square layer,

``residual_c2 = 1/2 * sum_{r != 0} G(r)**2`` and
``total_c2    = 1/2 * sum_r G(r)**2``.

All formal-series routines below use exact Python integers and
:class:`fractions.Fraction`.  The finite-torus numerical routines evaluate an
algebraically exact row transfer matrix either with mpmath or, when explicitly
requested, IEEE-754 binary64 arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb
from typing import Sequence

import mpmath as mp
import numpy as np

from ising.exact_enumeration import joint_dos
from ising.lattices import cubic

__all__ = [
    "InterlayerC2Series",
    "anisotropic_box_even_subgraph",
    "anisotropic_flm_c2_series",
    "critical_coupling_2d",
    "decoupled_chain_phi",
    "exact_2d_overlap_series",
    "log_cosh_w_series",
    "torus_overlap_sum_float",
    "torus_overlap_sum_mp",
]

RationalSeries = tuple[Fraction, ...]
IntegerSeries = tuple[int, ...]


@dataclass(frozen=True)
class InterlayerC2Series:
    """Exact ``w**2`` series obtained by anisotropic finite-lattice inversion.

    ``residual`` excludes the separate ``log(cosh(K_z))`` bond prefactor.
    ``total`` includes its exact ``w**2/2`` contribution.
    """

    order: int
    residual: RationalSeries
    total: RationalSeries
    box_weights: dict[tuple[int, int, int], RationalSeries]
    boxes: tuple[tuple[int, int, int], ...]
    bound_slack: int


def _check_order(order: int) -> int:
    value = int(order)
    if value < 0:
        raise ValueError("order must be non-negative")
    return value


def _truncated_mul(
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
    order: int,
) -> tuple[int | Fraction, ...]:
    out: list[int | Fraction] = [0] * (order + 1)
    for i, a in enumerate(left[: order + 1]):
        if not a:
            continue
        for j, b in enumerate(right[: order + 1 - i]):
            if b:
                out[i + j] += a * b
    return tuple(out)


def _series_divide(
    numerator: Sequence[int | Fraction],
    denominator: Sequence[int | Fraction],
    order: int,
) -> RationalSeries:
    if not denominator or not denominator[0]:
        raise ValueError("series denominator must have nonzero constant coefficient")
    den = [Fraction(denominator[k]) if k < len(denominator) else Fraction(0) for k in range(order + 1)]
    num = [Fraction(numerator[k]) if k < len(numerator) else Fraction(0) for k in range(order + 1)]
    quotient = [Fraction(0) for _ in range(order + 1)]
    for degree in range(order + 1):
        convolution = sum(
            (den[k] * quotient[degree - k] for k in range(1, degree + 1)),
            Fraction(0),
        )
        quotient[degree] = (num[degree] - convolution) / den[0]
    return tuple(quotient)


def log_cosh_w_series(order: int) -> RationalSeries:
    """Return ``log(cosh(atanh(w))) = -log(1-w**2)/2`` exactly."""

    order = _check_order(order)
    coefficients = [Fraction(0) for _ in range(order + 1)]
    for degree in range(2, order + 1, 2):
        coefficients[degree] = Fraction(1, degree)
    return tuple(coefficients)


@lru_cache(maxsize=None)
def _signed_binomial(satisfied: int, unsatisfied: int, order: int) -> IntegerSeries:
    """Coefficients of ``(1+x)^satisfied (1-x)^unsatisfied``."""

    out = [0] * (order + 1)
    for left_degree in range(min(satisfied, order) + 1):
        left = comb(satisfied, left_degree)
        for right_degree in range(min(unsatisfied, order - left_degree) + 1):
            out[left_degree + right_degree] += (
                left * comb(unsatisfied, right_degree) * (-1) ** right_degree
            )
    return tuple(out)


@lru_cache(maxsize=None)
def anisotropic_box_even_subgraph(
    shape: tuple[int, int, int],
    v_order: int,
    w_order: int = 2,
) -> tuple[tuple[int, ...], ...]:
    """Return the exact open-box ``P(v,w)`` coefficient table.

    Directions zero and one share ``v`` and direction two uses ``w``.  The
    returned table is indexed ``[v_degree][w_degree]``.  It is obtained by an
    exact transform of the integer anisotropic density of states, independently
    of the finite-lattice Möbius inversion.
    """

    box = tuple(int(side) for side in shape)
    if len(box) != 3 or any(side < 1 for side in box):
        raise ValueError("shape must contain three positive side lengths")
    v_order = _check_order(v_order)
    w_order = _check_order(w_order)

    lattice = cubic(*box, periodic=False)
    density = joint_dos(lattice, with_field=False)
    inplane_bonds = len(lattice.bonds_by_dir[0]) + len(lattice.bonds_by_dir[1])
    vertical_bonds = len(lattice.bonds_by_dir[2])
    inplane_factors = tuple(
        _signed_binomial(satisfied, inplane_bonds - satisfied, v_order)
        for satisfied in range(inplane_bonds + 1)
    )
    vertical_factors = tuple(
        _signed_binomial(satisfied, vertical_bonds - satisfied, w_order)
        for satisfied in range(vertical_bonds + 1)
    )

    accumulated = [[0 for _ in range(w_order + 1)] for _ in range(v_order + 1)]
    for index, raw_count in np.ndenumerate(density):
        count = int(raw_count)
        if not count:
            continue
        bx, by, bz = index
        pv = inplane_factors[bx + by]
        pw = vertical_factors[bz]
        for iv, coefficient_v in enumerate(pv):
            if not coefficient_v:
                continue
            for iw, coefficient_w in enumerate(pw):
                if coefficient_w:
                    accumulated[iv][iw] += count * coefficient_v * coefficient_w

    denominator = 1 << lattice.n_sites
    result: list[tuple[int, ...]] = []
    for row in accumulated:
        exact_row: list[int] = []
        for numerator in row:
            if numerator % denominator:
                raise AssertionError("anisotropic high-temperature transform is non-integral")
            exact_row.append(numerator // denominator)
        result.append(tuple(exact_row))
    if result[0][0] != 1:
        raise AssertionError("the empty even subgraph must be unique")
    if any(result[degree][1] for degree in range(v_order + 1)):
        raise AssertionError("an even subgraph crosses every layer cut an even number of times")
    return tuple(result)


def _box_log_w2(shape: tuple[int, int, int], order: int) -> RationalSeries:
    polynomial = anisotropic_box_even_subgraph(shape, order, 2)
    p0 = tuple(row[0] for row in polynomial)
    p2 = tuple(row[2] for row in polynomial)
    # Since [w]P=0, [w^2]log(P)=P_2/P_0 exactly.
    return _series_divide(p2, p0, order)


def anisotropic_flm_c2_series(
    order: int = 6,
    *,
    bound_slack: int = 0,
) -> InterlayerC2Series:
    """Compute ``[w**2]`` of the anisotropic bulk free energy exactly.

    A connected even graph with exactly two vertical edges has vertical extent
    one and needs at least ``2*((a-1)+(b-1))`` in-plane edges to span an
    ``a x b`` rectangle.  Hence the listed boxes are complete through the
    requested ``v`` order.  ``bound_slack`` adds provably irrelevant boxes as a
    stability check.
    """

    order = _check_order(order)
    bound_slack = int(bound_slack)
    if bound_slack < 0:
        raise ValueError("bound_slack must be non-negative")
    span_budget = order // 2 + bound_slack
    boxes = tuple(
        sorted(
            (
                (a, b, 2)
                for a in range(1, span_budget + 2)
                for b in range(1, span_budget + 2)
                if (a - 1) + (b - 1) <= span_budget
            ),
            key=lambda shape: (sum(shape), shape),
        )
    )

    weights: dict[tuple[int, int, int], RationalSeries] = {}
    bulk = [Fraction(0) for _ in range(order + 1)]
    for shape in boxes:
        a, b, _ = shape
        weight = list(_box_log_w2(shape, order))
        for subshape, subweight in weights.items():
            sa, sb, _ = subshape
            if sa <= a and sb <= b:
                placements = (a - sa + 1) * (b - sb + 1)
                for degree in range(order + 1):
                    weight[degree] -= placements * subweight[degree]
        exact_weight = tuple(weight)
        weights[shape] = exact_weight
        for degree in range(order + 1):
            bulk[degree] += exact_weight[degree]

    residual = tuple(bulk)
    total = list(residual)
    total[0] += Fraction(1, 2)
    return InterlayerC2Series(
        order=order,
        residual=residual,
        total=tuple(total),
        box_weights=weights,
        boxes=boxes,
        bound_slack=bound_slack,
    )


def _row_spin_table(width: int) -> np.ndarray:
    states = np.arange(1 << width, dtype=np.uint64)
    bits = ((states[:, None] >> np.arange(width, dtype=np.uint64)) & 1).astype(np.int8)
    return 2 * bits - 1


def _poly_forward_step(
    values: np.ndarray,
    transition: tuple[np.ndarray, ...],
    horizontal: np.ndarray,
    order: int,
) -> np.ndarray:
    inner = np.zeros(values.shape, dtype=object)
    for left_degree in range(order + 1):
        for edge_degree in range(order + 1 - left_degree):
            inner[left_degree + edge_degree] += values[left_degree] @ transition[edge_degree]
    out = np.zeros(values.shape, dtype=object)
    for inner_degree in range(order + 1):
        for row_degree in range(order + 1 - inner_degree):
            out[inner_degree + row_degree] += inner[inner_degree] * horizontal[row_degree]
    return out


def _poly_backward_step(
    values: np.ndarray,
    transition: tuple[np.ndarray, ...],
    horizontal: np.ndarray,
    order: int,
) -> np.ndarray:
    weighted = np.zeros(values.shape, dtype=object)
    for value_degree in range(order + 1):
        for row_degree in range(order + 1 - value_degree):
            weighted[value_degree + row_degree] += values[value_degree] * horizontal[row_degree]
    out = np.zeros(values.shape, dtype=object)
    for edge_degree in range(order + 1):
        for value_degree in range(order + 1 - edge_degree):
            out[edge_degree + value_degree] += transition[edge_degree] @ weighted[value_degree]
    return out


def _formal_ratio(numerator: Sequence[int], denominator: Sequence[int], order: int) -> RationalSeries:
    return _series_divide(numerator, denominator, order)


@lru_cache(maxsize=None)
def exact_2d_overlap_series(order: int = 6) -> RationalSeries:
    """Return ``sum_r G(r)^2`` on the infinite square lattice through ``v^order``.

    This is an independent exact two-dimensional calculation.  It uses an
    open strip of height ``2*order+1``, a periodic row of width ``order+1``, and
    exact integer polynomial row transfer.  The source is ``order`` rows from
    either open boundary.  A horizontal wrapping contribution first affects a
    squared correlation at degree ``order+1``; an open-boundary contribution
    must travel farther still.  Thus the returned truncation is the infinite-
    lattice series.  Runtime is intentionally bounded to the validation order.
    """

    order = _check_order(order)
    if order > 6:
        raise ValueError("the exact correlation transfer is certified and budgeted only through v^6")
    if order % 2:
        raise ValueError("overlap-susceptibility truncation order must be even")
    if order == 0:
        return (Fraction(1),)

    width = order + 1
    height = 2 * order + 1
    center = order
    spins = _row_spin_table(width)
    state_count = len(spins)

    horizontal_satisfied = np.sum(spins == np.roll(spins, -1, axis=1), axis=1)
    horizontal = np.empty((order + 1, state_count), dtype=object)
    for state, satisfied in enumerate(horizontal_satisfied):
        polynomial = _signed_binomial(int(satisfied), width - int(satisfied), order)
        for degree in range(order + 1):
            horizontal[degree, state] = polynomial[degree]

    state_indices = np.arange(state_count, dtype=np.uint64)
    xor = state_indices[:, None] ^ state_indices[None, :]
    different = np.bitwise_count(xor).astype(np.int16)
    transition_arrays = [np.empty((state_count, state_count), dtype=np.int64) for _ in range(order + 1)]
    transition_by_distance = [
        _signed_binomial(width - distance, distance, order) for distance in range(width + 1)
    ]
    for degree in range(order + 1):
        lookup = np.array(
            [transition_by_distance[distance][degree] for distance in range(width + 1)],
            dtype=np.int64,
        )
        transition_arrays[degree][:] = lookup[different]
    transition = tuple(transition_arrays)

    forward: list[np.ndarray] = [horizontal.copy()]
    for _ in range(1, height):
        forward.append(_poly_forward_step(forward[-1], transition, horizontal, order))

    backward: list[np.ndarray] = [np.empty((0, 0), dtype=object) for _ in range(height)]
    top = np.zeros((order + 1, state_count), dtype=object)
    top[0] = 1
    backward[-1] = top
    for row in range(height - 2, -1, -1):
        backward[row] = _poly_backward_step(backward[row + 1], transition, horizontal, order)

    partition = [0] * (order + 1)
    for degree in range(order + 1):
        partition[degree] = sum(
            int(np.sum(forward[center][left_degree] * backward[center][degree - left_degree]))
            for left_degree in range(degree + 1)
        )
    if partition[0] != 1 << (width * height):
        raise AssertionError("zero-coupling strip partition count is not 2^N")

    source = forward[center] * spins[:, 0]
    overlap = [Fraction(0) for _ in range(order + 1)]
    max_distance = order // 2
    current = source
    for vertical_distance in range(max_distance + 1):
        row = center + vertical_distance
        for x in range(width):
            horizontal_distance = min(x, width - x)
            if horizontal_distance + vertical_distance > max_distance:
                continue
            numerator = [0] * (order + 1)
            target_sign = spins[:, x]
            for degree in range(order + 1):
                numerator[degree] = sum(
                    int(
                        np.sum(
                            target_sign
                            * current[left_degree]
                            * backward[row][degree - left_degree]
                        )
                    )
                    for left_degree in range(degree + 1)
                )
            correlation = _formal_ratio(numerator, partition, order)
            square = tuple(Fraction(value) for value in _truncated_mul(correlation, correlation, order))
            multiplicity = 1 if vertical_distance == 0 else 2
            for degree in range(order + 1):
                overlap[degree] += multiplicity * square[degree]
        if vertical_distance < max_distance:
            current = _poly_forward_step(current, transition, horizontal, order)

    if overlap[0] != 1:
        raise AssertionError("only G(0)^2 contributes at v=0")
    if any(overlap[degree] for degree in range(1, order + 1, 2)):
        raise AssertionError("the square-lattice overlap sum must be even in v")
    return tuple(overlap)


def critical_coupling_2d(dps: int = 80) -> mp.mpf:
    with mp.workdps(int(dps)):
        return +mp.log(1 + mp.sqrt(2)) / 2


def decoupled_chain_phi(kz: mp.mpf | str | float, dps: int = 80) -> mp.mpf:
    """Exact thermodynamic reduced free energy at ``K_x=K_y=0``."""

    with mp.workdps(int(dps)):
        value = mp.mpf(kz)
        return +mp.log(2 * mp.cosh(value))


def _symmetric_row_transfer_float(side: int, coupling: float) -> tuple[np.ndarray, np.ndarray]:
    spins = _row_spin_table(side).astype(np.float64)
    horizontal_energy = np.sum(spins * np.roll(spins, -1, axis=1), axis=1)
    vertical_energy = spins @ spins.T
    transfer = np.exp(
        coupling
        * (
            vertical_energy
            + horizontal_energy[:, None] / 2
            + horizontal_energy[None, :] / 2
        )
    )
    return transfer, spins


def torus_overlap_sum_float(side: int, coupling: float) -> dict[str, float | int]:
    """Evaluate ``sum_r G_L(r)^2`` on an ``L x L`` torus in binary64.

    This fast path is used only for finite-size extrapolation.  The returned
    matrix residual and the independent mpmath comparison recorded by the
    experiment quantify numerical, separately from finite-size, error.
    """

    side = int(side)
    if side < 2:
        raise ValueError("side must be at least two under the repository bond convention")
    transfer, spins = _symmetric_row_transfer_float(side, float(coupling))
    powers = [np.eye(len(transfer), dtype=np.float64)]
    for _ in range(side):
        powers.append(powers[-1] @ transfer)
    partition = float(np.trace(powers[side]))
    source_sign = spins[:, 0]
    overlap = 0.0
    maximum_correlation_asymmetry = 0.0
    correlations: dict[tuple[int, int], float] = {}
    for y in range(side):
        left = powers[y]
        right = powers[side - y].T
        for x in range(side):
            target_sign = spins[:, x]
            numerator = np.sum(
                source_sign[:, None] * left * target_sign[None, :] * right
            )
            correlation = float(numerator / partition)
            correlations[(x, y)] = correlation
            overlap += correlation * correlation
    for (x, y), value in correlations.items():
        reflected = correlations[((-x) % side, (-y) % side)]
        maximum_correlation_asymmetry = max(
            maximum_correlation_asymmetry, abs(value - reflected)
        )
    identity_residual = abs(correlations[(0, 0)] - 1.0)
    return {
        "side": side,
        "overlap_sum": overlap,
        "total_c2": overlap / 2,
        "residual_c2": (overlap - 1) / 2,
        "partition": partition,
        "identity_residual": identity_residual,
        "maximum_correlation_asymmetry": maximum_correlation_asymmetry,
    }


def torus_overlap_sum_mp(
    side: int,
    coupling: mp.mpf | str | float,
    *,
    dps: int = 60,
) -> dict[str, mp.mpf | int]:
    """Arbitrary-precision finite-torus counterpart of :func:`torus_overlap_sum_float`."""

    side = int(side)
    if side < 2:
        raise ValueError("side must be at least two under the repository bond convention")
    if side > 6:
        raise ValueError("the arbitrary-precision dense transfer is budgeted for side <= 6")
    with mp.workdps(int(dps)):
        spins_array = _row_spin_table(side)
        states = [[int(value) for value in row] for row in spins_array]
        count = len(states)
        horizontal = [
            sum(row[i] * row[(i + 1) % side] for i in range(side)) for row in states
        ]
        k = mp.mpf(coupling)
        transfer = mp.matrix(count, count)
        for left in range(count):
            for right in range(count):
                vertical = sum(
                    states[left][i] * states[right][i] for i in range(side)
                )
                transfer[left, right] = mp.exp(
                    k * (vertical + mp.mpf(horizontal[left] + horizontal[right]) / 2)
                )
        powers = [mp.eye(count)]
        for _ in range(side):
            powers.append(powers[-1] * transfer)
        partition = mp.fsum(powers[side][state, state] for state in range(count))
        overlap = mp.mpf(0)
        identity = None
        for y in range(side):
            left_power = powers[y]
            right_power = powers[side - y]
            for x in range(side):
                numerator = mp.fsum(
                    states[left][0]
                    * left_power[left, right]
                    * states[right][x]
                    * right_power[right, left]
                    for left in range(count)
                    for right in range(count)
                )
                correlation = numerator / partition
                if x == 0 and y == 0:
                    identity = correlation
                overlap += correlation * correlation
        if identity is None:
            raise AssertionError("origin correlation was not evaluated")
        return {
            "side": side,
            "overlap_sum": +overlap,
            "total_c2": +(overlap / 2),
            "residual_c2": +((overlap - 1) / 2),
            "partition": +partition,
            "identity_residual": +abs(identity - 1),
        }
