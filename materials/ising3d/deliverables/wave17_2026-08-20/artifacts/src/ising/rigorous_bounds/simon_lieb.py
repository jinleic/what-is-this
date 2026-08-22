"""Computer-assisted Simon--Lieb lower bounds for the simple-cubic Ising model.

For a finite open box ``B`` and ``v = tanh(K)``, write

``P(v) = sum_{F: partial F = empty} v**|F|`` and
``Q(v) = sum_x q_B(x) sum_{F: partial F = {origin, x}} v**|F|``,

where ``q_B(x)`` is the number of nearest-neighbour edges from ``x`` to the
complement of ``B``.  The Simon--Lieb coefficient is exactly

``kappa_B(K) = v * Q(v) / P(v)``.

The routines below provide three independent arithmetic paths:

* exhaustive spin enumeration (with one spin fixed by global spin flip),
* exact integer parity transfer at rational ``v``, and
* positive-term float/mpmath parity transfer for root finding and auditing.

All boxes have free boundary conditions.  Bits within a transfer layer are
ordered with the first coordinate fastest; the third coordinate is the
propagation direction.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from typing import Iterable

import mpmath as mp
import numpy as np


Shape3D = tuple[int, int, int]
Origin3D = tuple[int, int, int]

__all__ = [
    "atanh_rational_interval",
    "boundary_multiplicities",
    "enumerate_pq_polynomials",
    "exact_criterion_residual",
    "exact_pq_at_rational",
    "kappa_float_from_v",
    "kappa_mpf_from_v",
    "mpmath_relative_roundoff_bound",
    "solve_box_root_float",
]


@dataclass(frozen=True)
class _Geometry:
    shape: Shape3D
    origin: Origin3D
    layer_sites: int
    state_count: int
    indices: np.ndarray
    popcounts: np.ndarray
    intralayer_edge_masks: tuple[int, ...]
    boundary_by_layer: tuple[tuple[int, ...], ...]
    origin_layer: int
    origin_bit: int
    site_count: int
    edge_count: int


def _validate_shape(shape: Iterable[int]) -> Shape3D:
    values = tuple(shape)
    if len(values) != 3:
        raise ValueError("shape must contain exactly three side lengths")
    if any(isinstance(side, bool) or not isinstance(side, int) or side <= 0 for side in values):
        raise ValueError("all side lengths must be positive integers")
    return values  # type: ignore[return-value]


def _validate_origin(shape: Shape3D, origin: Iterable[int] | None) -> Origin3D:
    if origin is None:
        return tuple((side - 1) // 2 for side in shape)  # type: ignore[return-value]
    values = tuple(origin)
    if len(values) != 3:
        raise ValueError("origin must contain exactly three coordinates")
    if any(
        isinstance(coordinate, bool)
        or not isinstance(coordinate, int)
        or coordinate < 0
        or coordinate >= side
        for coordinate, side in zip(values, shape, strict=True)
    ):
        raise ValueError("origin must be an integer site inside the box")
    return values  # type: ignore[return-value]


def _edge_count(shape: Shape3D) -> int:
    nx, ny, nz = shape
    return (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1)


def boundary_multiplicities(shape: Iterable[int]) -> tuple[int, ...]:
    """Return ``q_B(x)=6-deg_B(x)`` in ``z,y,x`` (``x`` fastest) order."""

    nx, ny, nz = _validate_shape(shape)
    values: list[int] = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                internal_degree = (
                    int(x > 0)
                    + int(x + 1 < nx)
                    + int(y > 0)
                    + int(y + 1 < ny)
                    + int(z > 0)
                    + int(z + 1 < nz)
                )
                values.append(6 - internal_degree)
    return tuple(values)


@lru_cache(maxsize=16)
def _geometry(shape: Shape3D, origin: Origin3D) -> _Geometry:
    nx, ny, nz = shape
    layer_sites = nx * ny
    if layer_sites >= 63:
        raise ValueError("a transfer cross-section must contain fewer than 63 sites")
    state_count = 1 << layer_sites
    index_dtype = np.uint32 if layer_sites <= 32 else np.uint64
    indices = np.arange(state_count, dtype=index_dtype)
    popcounts = np.bitwise_count(indices).astype(np.uint8, copy=False)

    edge_masks: list[int] = []
    for y in range(ny):
        for x in range(nx):
            site = x + nx * y
            if x + 1 < nx:
                edge_masks.append((1 << site) | (1 << (site + 1)))
            if y + 1 < ny:
                edge_masks.append((1 << site) | (1 << (site + nx)))

    flat_boundary = boundary_multiplicities(shape)
    boundary_by_layer = tuple(
        tuple(flat_boundary[z * layer_sites : (z + 1) * layer_sites])
        for z in range(nz)
    )
    ox, oy, oz = origin
    return _Geometry(
        shape=shape,
        origin=origin,
        layer_sites=layer_sites,
        state_count=state_count,
        indices=indices,
        popcounts=popcounts,
        intralayer_edge_masks=tuple(edge_masks),
        boundary_by_layer=boundary_by_layer,
        origin_layer=oz,
        origin_bit=1 << (ox + nx * oy),
        site_count=nx * ny * nz,
        edge_count=_edge_count(shape),
    )


def _get_geometry(
    shape: Iterable[int], origin: Iterable[int] | None = None
) -> _Geometry:
    checked_shape = _validate_shape(shape)
    checked_origin = _validate_origin(checked_shape, origin)
    return _geometry(checked_shape, checked_origin)


def _transfer_pq_float(geometry: _Geometry, v: float) -> tuple[float, float]:
    n = geometry.state_count
    plain = np.zeros(n, dtype=np.float64)
    fixed = np.zeros(n, dtype=np.float64)
    marked = np.zeros(n, dtype=np.float64)
    plain[0] = 1.0
    fixed[0] = 1.0
    vertical_weights = np.power(v, geometry.popcounts, dtype=np.float64)

    nz = geometry.shape[2]
    for z in range(nz):
        for edge_mask in geometry.intralayer_edge_masks:
            shifted = geometry.indices ^ edge_mask
            plain = plain + v * plain[shifted]
            fixed = fixed + v * fixed[shifted]
            marked = marked + v * marked[shifted]

        fixed_target = geometry.origin_bit if z == geometry.origin_layer else 0
        boundary = geometry.boundary_by_layer[z]
        if z + 1 == nz:
            p_value = float(plain[0])
            q_value = float(marked[fixed_target])
            for site, multiplicity in enumerate(boundary):
                q_value += multiplicity * float(fixed[fixed_target ^ (1 << site)])
            return p_value, q_value

        shifted_target = geometry.indices ^ fixed_target
        next_marked = marked[shifted_target].copy()
        for site, multiplicity in enumerate(boundary):
            if multiplicity:
                next_marked += multiplicity * fixed[shifted_target ^ (1 << site)]
        marked = next_marked * vertical_weights
        plain = plain * vertical_weights
        fixed = fixed[shifted_target] * vertical_weights

    raise AssertionError("unreachable transfer termination")


def kappa_float_from_v(
    shape: Iterable[int], v: float, origin: Iterable[int] | None = None
) -> float:
    """Evaluate ``kappa_B=v*Q/P`` with positive-term binary64 transfer."""

    if not isinstance(v, (float, int)) or not math.isfinite(float(v)):
        raise ValueError("v must be a finite real number")
    v_float = float(v)
    if not 0.0 <= v_float < 1.0:
        raise ValueError("v must satisfy 0 <= v < 1")
    geometry = _get_geometry(shape, origin)
    p_value, q_value = _transfer_pq_float(geometry, v_float)
    return v_float * q_value / p_value


def solve_box_root_float(
    shape: Iterable[int],
    origin: Iterable[int] | None = None,
    *,
    iterations: int = 58,
) -> tuple[float, float]:
    """Return a binary64 bisection estimate ``(v_B, atanh(v_B))``.

    This point estimate is only a locator.  A rigorous bound must evaluate a
    downward-rounded rational ``v`` with ``exact_criterion_residual`` or the
    audited multiprecision evaluator.
    """

    if not isinstance(iterations, int) or iterations < 20:
        raise ValueError("iterations must be an integer of at least 20")
    geometry = _get_geometry(shape, origin)
    low = 0.0
    high = 0.5
    if high * _transfer_pq_float(geometry, high)[1] <= _transfer_pq_float(
        geometry, high
    )[0]:
        raise ValueError("failed to bracket the finite-box root below v=1/2")
    for _ in range(iterations):
        midpoint = (low + high) / 2.0
        p_value, q_value = _transfer_pq_float(geometry, midpoint)
        if midpoint * q_value < p_value:
            low = midpoint
        else:
            high = midpoint
    v_root = (low + high) / 2.0
    return v_root, math.atanh(v_root)


def _integer_vertical_weights(
    geometry: _Geometry, numerator: int, denominator: int
) -> np.ndarray:
    weights_by_count = np.asarray(
        [
            numerator**count * denominator ** (geometry.layer_sites - count)
            for count in range(geometry.layer_sites + 1)
        ],
        dtype=object,
    )
    return weights_by_count[geometry.popcounts]


def _validate_rational(numerator: int, denominator: int) -> tuple[int, int]:
    if (
        isinstance(numerator, bool)
        or not isinstance(numerator, int)
        or isinstance(denominator, bool)
        or not isinstance(denominator, int)
        or denominator <= 0
        or numerator < 0
        or numerator >= denominator
    ):
        raise ValueError("the rational must satisfy 0 <= numerator < denominator")
    divisor = math.gcd(numerator, denominator)
    return numerator // divisor, denominator // divisor


def exact_pq_at_rational(
    shape: Iterable[int],
    numerator: int,
    denominator: int,
    origin: Iterable[int] | None = None,
) -> tuple[int, int]:
    """Return exact scaled values ``(q**E P(p/q), q**E Q(p/q))``.

    Here ``p/q`` is reduced internally and ``E`` is the number of internal
    box edges.  Python integers and object arrays are used throughout.
    """

    p, q = _validate_rational(numerator, denominator)
    geometry = _get_geometry(shape, origin)
    n = geometry.state_count
    plain = np.zeros(n, dtype=object)
    fixed = np.zeros(n, dtype=object)
    marked = np.zeros(n, dtype=object)
    plain[0] = 1
    fixed[0] = 1
    vertical_weights = _integer_vertical_weights(geometry, p, q)

    nz = geometry.shape[2]
    for z in range(nz):
        for edge_mask in geometry.intralayer_edge_masks:
            shifted = geometry.indices ^ edge_mask
            plain = q * plain + p * plain[shifted]
            fixed = q * fixed + p * fixed[shifted]
            marked = q * marked + p * marked[shifted]

        fixed_target = geometry.origin_bit if z == geometry.origin_layer else 0
        boundary = geometry.boundary_by_layer[z]
        if z + 1 == nz:
            p_value = int(plain[0])
            q_value = int(marked[fixed_target])
            for site, multiplicity in enumerate(boundary):
                q_value += multiplicity * int(fixed[fixed_target ^ (1 << site)])
            return p_value, q_value

        shifted_target = geometry.indices ^ fixed_target
        next_marked = marked[shifted_target].copy()
        for site, multiplicity in enumerate(boundary):
            if multiplicity:
                next_marked += multiplicity * fixed[shifted_target ^ (1 << site)]
        marked = next_marked * vertical_weights
        plain = plain * vertical_weights
        fixed = fixed[shifted_target] * vertical_weights

    raise AssertionError("unreachable transfer termination")


def exact_criterion_residual(
    shape: Iterable[int],
    numerator: int,
    denominator: int,
    origin: Iterable[int] | None = None,
) -> int:
    """Return an exact integer with the sign of ``v*Q(v)-P(v)``.

    This two-stream transfer directly accumulates
    ``p*sum_x q_B(x) P_{0x}(p/q) - q*P(p/q)``.  It avoids constructing the
    two enormous positive integers separately when only the safe side of the
    criterion is required.
    """

    p, q = _validate_rational(numerator, denominator)
    geometry = _get_geometry(shape, origin)
    n = geometry.state_count
    fixed = np.zeros(n, dtype=object)
    combined = np.zeros(n, dtype=object)
    fixed[0] = 1
    vertical_weights = _integer_vertical_weights(geometry, p, q)
    ox, oy, oz = geometry.origin
    origin_site = ox + geometry.shape[0] * oy

    nz = geometry.shape[2]
    for z in range(nz):
        for edge_mask in geometry.intralayer_edge_masks:
            shifted = geometry.indices ^ edge_mask
            fixed = q * fixed + p * fixed[shifted]
            combined = q * combined + p * combined[shifted]

        fixed_target = geometry.origin_bit if z == geometry.origin_layer else 0
        boundary = geometry.boundary_by_layer[z]
        if z + 1 == nz:
            residual = int(combined[fixed_target])
            for site, multiplicity in enumerate(boundary):
                branch_weight = p * multiplicity - q * int(z == oz and site == origin_site)
                if branch_weight:
                    residual += branch_weight * int(
                        fixed[fixed_target ^ (1 << site)]
                    )
            return residual

        shifted_target = geometry.indices ^ fixed_target
        next_combined = combined[shifted_target].copy()
        for site, multiplicity in enumerate(boundary):
            branch_weight = p * multiplicity - q * int(z == oz and site == origin_site)
            if branch_weight:
                next_combined += branch_weight * fixed[
                    shifted_target ^ (1 << site)
                ]
        combined = next_combined * vertical_weights
        fixed = fixed[shifted_target] * vertical_weights

    raise AssertionError("unreachable transfer termination")


def _poly_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_power, left_value in enumerate(left):
        if left_value:
            for right_power, right_value in enumerate(right):
                if right_value:
                    result[left_power + right_power] += left_value * right_value
    return result


def _histogram_to_ht_polynomial(
    histogram: np.ndarray, edge_count: int, normalization: int
) -> tuple[int, ...]:
    positive_powers: list[list[int]] = [[1]]
    negative_powers: list[list[int]] = [[1]]
    for _ in range(edge_count):
        positive_powers.append(_poly_multiply(positive_powers[-1], [1, 1]))
        negative_powers.append(_poly_multiply(negative_powers[-1], [1, -1]))

    accumulator = [0] * (edge_count + 1)
    for satisfied, raw_count in enumerate(histogram):
        count = int(raw_count)
        if not count:
            continue
        term = _poly_multiply(
            positive_powers[satisfied], negative_powers[edge_count - satisfied]
        )
        for power, coefficient in enumerate(term):
            accumulator[power] += count * coefficient

    coefficients: list[int] = []
    for value in accumulator:
        quotient, remainder = divmod(value, normalization)
        if remainder:
            raise AssertionError("spin histogram did not produce an integer HT coefficient")
        coefficients.append(quotient)
    return tuple(coefficients)


def enumerate_pq_polynomials(
    shape: Iterable[int],
    origin: Iterable[int] | None = None,
    *,
    chunk_size: int = 1 << 22,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Exhaustively enumerate the exact integer polynomials ``P`` and ``Q``.

    One spin is fixed to ``+1``; global spin flip accounts for the omitted
    half of the configurations.  The implementation refuses more than 30
    sites, matching the repository's brute-force policy.  Integer histograms
    are accumulated in chunks.  ``numpy.bincount``'s weighted binary64 path
    is exact here because every partial integer sum is explicitly bounded by
    ``2**53`` and checked to be integral before conversion back to int64.
    """

    geometry = _get_geometry(shape, origin)
    if geometry.site_count > 30:
        raise ValueError(
            f"brute-force enumeration refuses N={geometry.site_count} > 30"
        )
    if not isinstance(chunk_size, int) or chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    nx, ny, nz = geometry.shape
    site_count = geometry.site_count
    edge_count = geometry.edge_count
    origin_index = (
        geometry.origin[0]
        + nx * geometry.origin[1]
        + nx * ny * geometry.origin[2]
    )

    directional_masks: list[tuple[int, int]] = []
    for stride, coordinates in (
        (1, ((x, y, z) for z in range(nz) for y in range(ny) for x in range(nx - 1))),
        (nx, ((x, y, z) for z in range(nz) for y in range(ny - 1) for x in range(nx))),
        (
            nx * ny,
            ((x, y, z) for z in range(nz - 1) for y in range(ny) for x in range(nx)),
        ),
    ):
        mask = 0
        for x, y, z in coordinates:
            mask |= 1 << (x + nx * y + nx * ny * z)
        directional_masks.append((stride, mask))

    boundary = boundary_multiplicities(geometry.shape)
    boundary_masks: dict[int, int] = {}
    for site, multiplicity in enumerate(boundary):
        boundary_masks[multiplicity] = boundary_masks.get(multiplicity, 0) | (1 << site)
    total_boundary = sum(boundary)
    fixed_state_count = 1 << (site_count - 1)
    if fixed_state_count * total_boundary >= 1 << 53:
        raise ValueError("weighted bincount would exceed exact binary64 integer range")

    energy_histogram = np.zeros(edge_count + 1, dtype=np.int64)
    correlation_histogram = np.zeros(edge_count + 1, dtype=np.int64)
    lower_mask = (1 << origin_index) - 1
    origin_mask = np.uint32(1 << origin_index)

    for start in range(0, fixed_state_count, chunk_size):
        stop = min(start + chunk_size, fixed_state_count)
        compressed = np.arange(start, stop, dtype=np.uint32)
        states = (
            (compressed & np.uint32(lower_mask))
            | ((compressed >> np.uint32(origin_index)) << np.uint32(origin_index + 1))
            | origin_mask
        )

        unsatisfied = np.zeros(stop - start, dtype=np.int16)
        for stride, mask in directional_masks:
            differences = (states ^ (states >> np.uint32(stride))) & np.uint32(mask)
            unsatisfied += np.bitwise_count(differences).astype(np.int16, copy=False)
        satisfied = edge_count - unsatisfied
        energy_histogram += np.bincount(satisfied, minlength=edge_count + 1).astype(
            np.int64, copy=False
        )

        weighted_spin_sum = np.full(stop - start, -total_boundary, dtype=np.int16)
        for multiplicity, mask in boundary_masks.items():
            weighted_spin_sum += (
                2
                * multiplicity
                * np.bitwise_count(states & np.uint32(mask)).astype(np.int16, copy=False)
            )
        weighted_counts = np.bincount(
            satisfied,
            weights=weighted_spin_sum.astype(np.float64),
            minlength=edge_count + 1,
        )
        rounded = np.rint(weighted_counts)
        if not np.array_equal(weighted_counts, rounded):
            raise AssertionError("weighted integer histogram lost exactness")
        correlation_histogram += rounded.astype(np.int64)

    if int(energy_histogram.sum()) != fixed_state_count:
        raise AssertionError("enumeration histogram does not sum to 2**(N-1)")
    normalization = fixed_state_count
    p_polynomial = _histogram_to_ht_polynomial(
        energy_histogram, edge_count, normalization
    )
    q_polynomial = _histogram_to_ht_polynomial(
        correlation_histogram, edge_count, normalization
    )
    if p_polynomial[0] != 1 or any(value < 0 for value in p_polynomial):
        raise AssertionError("invalid even-subgraph polynomial")
    if any(value < 0 for value in q_polynomial):
        raise AssertionError("invalid two-defect correlation polynomial")
    return p_polynomial, q_polynomial


def kappa_mpf_from_v(
    shape: Iterable[int],
    v: object,
    origin: Iterable[int] | None = None,
    *,
    dps: int = 50,
) -> mp.mpf:
    """Evaluate ``kappa_B`` with positive-term mpmath transfer at ``dps``."""

    if not isinstance(dps, int) or dps < 30:
        raise ValueError("dps must be an integer of at least 30")
    geometry = _get_geometry(shape, origin)
    with mp.workdps(dps):
        value = mp.mpf(v)
        if not mp.isfinite(value) or value < 0 or value >= 1:
            raise ValueError("v must satisfy 0 <= v < 1")

        powers = [mp.mpf(1)]
        for _ in range(geometry.layer_sites):
            powers.append(powers[-1] * value)
        vertical_weights = np.asarray(powers, dtype=object)[geometry.popcounts]

        n = geometry.state_count
        zero = mp.mpf(0)
        plain = np.full(n, zero, dtype=object)
        fixed = np.full(n, zero, dtype=object)
        marked = np.full(n, zero, dtype=object)
        plain[0] = mp.mpf(1)
        fixed[0] = mp.mpf(1)

        nz = geometry.shape[2]
        for z in range(nz):
            for edge_mask in geometry.intralayer_edge_masks:
                shifted = geometry.indices ^ edge_mask
                plain = plain + value * plain[shifted]
                fixed = fixed + value * fixed[shifted]
                marked = marked + value * marked[shifted]

            fixed_target = geometry.origin_bit if z == geometry.origin_layer else 0
            boundary = geometry.boundary_by_layer[z]
            if z + 1 == nz:
                p_value = plain[0]
                q_value = marked[fixed_target]
                for site, multiplicity in enumerate(boundary):
                    q_value += multiplicity * fixed[fixed_target ^ (1 << site)]
                return +(value * q_value / p_value)

            shifted_target = geometry.indices ^ fixed_target
            next_marked = marked[shifted_target].copy()
            for site, multiplicity in enumerate(boundary):
                if multiplicity:
                    next_marked += multiplicity * fixed[
                        shifted_target ^ (1 << site)
                    ]
            marked = next_marked * vertical_weights
            plain = plain * vertical_weights
            fixed = fixed[shifted_target] * vertical_weights

    raise AssertionError("unreachable transfer termination")


def mpmath_relative_roundoff_bound(
    shape: Iterable[int], dps: int, origin: Iterable[int] | None = None
) -> mp.mpf:
    """Conservative relative-error allowance for ``kappa_mpf_from_v``.

    Every transfer operation is on non-negative numbers, so there is no
    cancellation in ``P`` or ``Q``.  We overcount the arithmetic depth by
    ``32*(E+N+1)`` and use a unit roundoff ``10**(-(dps-5))``, five decimal
    orders larger than mpmath's working unit.  The returned ``8*gamma_n``
    also covers forming ``v*Q/P`` and the two independently rounded streams.
    """

    if not isinstance(dps, int) or dps < 30:
        raise ValueError("dps must be an integer of at least 30")
    geometry = _get_geometry(shape, origin)
    operation_depth = 32 * (geometry.edge_count + geometry.site_count + 1)
    with mp.workdps(max(30, dps)):
        unit = mp.power(10, -(dps - 5))
        product = operation_depth * unit
        if product >= mp.mpf("0.5"):
            raise ValueError("precision too low for the roundoff estimate")
        gamma = product / (1 - product)
        return +(8 * gamma)


def _interval_endpoints(value: object) -> tuple[str, str]:
    text = str(value).strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError(f"unexpected mpmath interval representation: {text}")
    lower, upper = text[1:-1].split(",", maxsplit=1)
    return lower.strip(), upper.strip()


def atanh_rational_interval(
    numerator: int, denominator: int, *, dps: int = 80
) -> tuple[str, str]:
    """Directed-rounding enclosure of ``atanh(numerator/denominator)``."""

    p, q = _validate_rational(numerator, denominator)
    if not isinstance(dps, int) or dps < 30:
        raise ValueError("dps must be an integer of at least 30")
    previous_dps = mp.iv.dps
    try:
        mp.iv.dps = dps
        value = mp.iv.mpf(p) / q
        interval = mp.iv.log((1 + value) / (1 - value)) / 2
        return _interval_endpoints(interval)
    finally:
        mp.iv.dps = previous_dps
