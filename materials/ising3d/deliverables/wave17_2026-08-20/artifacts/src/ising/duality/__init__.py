r"""Exact finite-cell Ising--:math:`\mathbb Z_2` gauge duality.

The implementation uses the cellular chain complex of a three-dimensional rectangular cubical
lattice over :math:`\mathbb F_2`.  Links carry gauge variables, plaquettes carry the gauge action,
and cubes carry the dual Ising spins.  Every combinatorial calculation is performed with Python
integers (or with overflow-checked ``numpy.int64`` counters during brute-force enumeration).

For a complex with boundary maps ``d3: C3 -> C2`` and ``d2: C2 -> C1``, the exact polynomial
identity checked here is

``sum_[h] Q_h(x) = 2**b3 * S(x)``,

where ``S`` enumerates closed plaquette surfaces ``ker(d2)``, ``[h]`` runs over
``H2 = ker(d2) / im(d3)``, and ``Q_h`` enumerates the broken bonds of the dual signed Ising model.
The corresponding partition-function theorem, including its normalization, is derived in
``proofs/duality_3d.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import prod
from typing import Iterable, Sequence, cast

import numpy as np

__all__ = [
    "CubicalComplex",
    "closed_surface_polynomial",
    "cubical_complex",
    "dual_ising_sector_polynomials",
    "flat_gauge_orbit_count",
    "gauge_orbit_size",
    "gauge_satisfied_plaquette_dos",
    "homology_dimensions",
    "ising_torus_sector_polynomials_2d",
    "kw_character",
    "signed_high_temperature_polynomial",
    "surface_polynomial_from_gauge_dos",
    "surface_sector_basis",
    "surface_sector_representatives",
]

_ENUMERATION_CHUNK = 1 << 20

Link = tuple[tuple[int, int, int], int]
Plaquette = tuple[tuple[int, int, int], int, int]


@dataclass(frozen=True)
class CubicalComplex:
    """A finite three-dimensional rectangular cubical cell complex.

    Cells are distinct even when a side has length two and periodicity makes two cells share the
    same endpoints.  This is the same doubled-cell convention as the repository's doubled-bond
    convention.  Boundary maps are packed as Python integers: bit ``j`` of a plaquette boundary
    names link ``j`` and bit ``p`` of a cube boundary names plaquette ``p``.
    """

    shape: tuple[int, int, int]
    periodic: tuple[bool, bool, bool]
    links: tuple[Link, ...]
    plaquettes: tuple[Plaquette, ...]
    cubes: tuple[tuple[int, int, int], ...]
    link_boundaries: tuple[int, ...]
    plaquette_boundaries: tuple[int, ...]
    cube_boundaries: tuple[int, ...]

    @property
    def n_vertices(self) -> int:
        return prod(self.shape)

    @property
    def n_sites(self) -> int:
        """Alias for ``n_vertices`` (the sites on which gauge transformations live)."""

        return self.n_vertices

    @property
    def n_links(self) -> int:
        return len(self.links)

    @property
    def n_plaquettes(self) -> int:
        return len(self.plaquettes)

    @property
    def n_cubes(self) -> int:
        return len(self.cubes)

    def describe(self) -> str:
        bc = "".join("P" if value else "O" for value in self.periodic)
        size = "x".join(str(length) for length in self.shape)
        return (
            f"{size}[{bc}] V={self.n_vertices} L={self.n_links} "
            f"P={self.n_plaquettes} C={self.n_cubes}"
        )


def _normalise_periodic(periodic: bool | Sequence[bool]) -> tuple[bool, bool, bool]:
    if isinstance(periodic, bool):
        return (periodic, periodic, periodic)
    result = tuple(bool(value) for value in periodic)
    if len(result) != 3:
        raise ValueError("periodic must be a bool or a length-three sequence")
    return cast(tuple[bool, bool, bool], result)


def cubical_complex(
    shape: Sequence[int], periodic: bool | Sequence[bool] = True
) -> CubicalComplex:
    """Construct the link--plaquette--cube complex for a rectangular three-dimensional lattice.

    All side lengths must be at least two.  In a periodic direction every coordinate roots a cell;
    in an open direction only coordinates ``0, ..., L-2`` do.  In particular, a periodic side of
    length two contains two distinct parallel links, two distinct layers of plaquettes, and two
    distinct layers of cubes.
    """

    shape_tuple = tuple(int(length) for length in shape)
    if len(shape_tuple) != 3:
        raise ValueError("a gauge cubical complex must have exactly three dimensions")
    if any(length < 2 for length in shape_tuple):
        raise ValueError("all gauge-lattice side lengths must be at least two")
    shape3 = cast(tuple[int, int, int], shape_tuple)
    periodic3 = _normalise_periodic(periodic)

    vertex_coords = tuple(product(*(range(length) for length in shape3)))
    vertex_index = {coord: index for index, coord in enumerate(vertex_coords)}

    def advance(coord: tuple[int, int, int], axis: int) -> tuple[int, int, int]:
        result = list(coord)
        result[axis] += 1
        if result[axis] == shape3[axis]:
            if not periodic3[axis]:
                raise ValueError("attempted to advance through an open boundary")
            result[axis] = 0
        return cast(tuple[int, int, int], tuple(result))

    span_counts = tuple(
        shape3[axis] if periodic3[axis] else shape3[axis] - 1 for axis in range(3)
    )

    links: list[Link] = []
    link_boundaries: list[int] = []
    for axis in range(3):
        ranges = [range(length) for length in shape3]
        ranges[axis] = range(span_counts[axis])
        for coord in product(*ranges):
            endpoint = advance(coord, axis)
            links.append((coord, axis))
            link_boundaries.append((1 << vertex_index[coord]) ^ (1 << vertex_index[endpoint]))
    link_index = {link: index for index, link in enumerate(links)}

    plaquettes: list[Plaquette] = []
    plaquette_boundaries: list[int] = []
    for axis_a, axis_b in combinations(range(3), 2):
        ranges = [range(length) for length in shape3]
        ranges[axis_a] = range(span_counts[axis_a])
        ranges[axis_b] = range(span_counts[axis_b])
        for coord in product(*ranges):
            shifted_a = advance(coord, axis_a)
            shifted_b = advance(coord, axis_b)
            edge_keys = (
                (coord, axis_a),
                (shifted_a, axis_b),
                (shifted_b, axis_a),
                (coord, axis_b),
            )
            boundary = 0
            for edge_key in edge_keys:
                boundary ^= 1 << link_index[edge_key]
            plaquettes.append((coord, axis_a, axis_b))
            plaquette_boundaries.append(boundary)
    plaquette_index = {face: index for index, face in enumerate(plaquettes)}

    cubes: list[tuple[int, int, int]] = []
    cube_boundaries: list[int] = []
    cube_ranges = [range(count) for count in span_counts]
    for coord in product(*cube_ranges):
        boundary = 0
        for normal_axis in range(3):
            tangent_axes = tuple(axis for axis in range(3) if axis != normal_axis)
            upper_coord = advance(coord, normal_axis)
            lower_face = (coord, tangent_axes[0], tangent_axes[1])
            upper_face = (upper_coord, tangent_axes[0], tangent_axes[1])
            boundary ^= 1 << plaquette_index[lower_face]
            boundary ^= 1 << plaquette_index[upper_face]
        cubes.append(coord)
        cube_boundaries.append(boundary)

    for cube_boundary in cube_boundaries:
        link_boundary = 0
        remaining = cube_boundary
        while remaining:
            least_bit = remaining & -remaining
            face_index = least_bit.bit_length() - 1
            link_boundary ^= plaquette_boundaries[face_index]
            remaining ^= least_bit
        if link_boundary:
            raise AssertionError("constructed boundary maps violate d2 o d3 = 0")

    return CubicalComplex(
        shape=shape3,
        periodic=periodic3,
        links=tuple(links),
        plaquettes=tuple(plaquettes),
        cubes=tuple(cubes),
        link_boundaries=tuple(link_boundaries),
        plaquette_boundaries=tuple(plaquette_boundaries),
        cube_boundaries=tuple(cube_boundaries),
    )


def _gf2_reduce(vector: int, basis: dict[int, int]) -> int:
    result = vector
    while result:
        pivot = result.bit_length() - 1
        basis_vector = basis.get(pivot)
        if basis_vector is None:
            break
        result ^= basis_vector
    return result


def _gf2_add(vector: int, basis: dict[int, int]) -> bool:
    remainder = _gf2_reduce(vector, basis)
    if not remainder:
        return False
    basis[remainder.bit_length() - 1] = remainder
    return True


def _gf2_rank(vectors: Iterable[int]) -> int:
    basis: dict[int, int] = {}
    for vector in vectors:
        _gf2_add(int(vector), basis)
    return len(basis)


def _kernel_basis_of_columns(columns: Sequence[int]) -> tuple[int, ...]:
    """Return a basis of relations among packed GF(2) column vectors."""

    pivots: dict[int, tuple[int, int]] = {}
    kernel: list[int] = []
    for column_index, column in enumerate(columns):
        image = int(column)
        combination = 1 << column_index
        while image:
            pivot = image.bit_length() - 1
            existing = pivots.get(pivot)
            if existing is None:
                pivots[pivot] = (image, combination)
                break
            image ^= existing[0]
            combination ^= existing[1]
        else:
            kernel.append(combination)
    return tuple(kernel)


def homology_dimensions(complex_: CubicalComplex) -> dict[str, int]:
    """Return exact GF(2) boundary ranks and Betti numbers ``b0, ..., b3``."""

    rank_d1 = _gf2_rank(complex_.link_boundaries)
    rank_d2 = _gf2_rank(complex_.plaquette_boundaries)
    rank_d3 = _gf2_rank(complex_.cube_boundaries)
    result = {
        "rank_d1": rank_d1,
        "rank_d2": rank_d2,
        "rank_d3": rank_d3,
        "b0": complex_.n_vertices - rank_d1,
        "b1": complex_.n_links - rank_d1 - rank_d2,
        "b2": complex_.n_plaquettes - rank_d2 - rank_d3,
        "b3": complex_.n_cubes - rank_d3,
    }
    if any(result[f"b{degree}"] < 0 for degree in range(4)):
        raise AssertionError("negative Betti number: invalid chain complex")
    euler_cells = (
        complex_.n_vertices - complex_.n_links + complex_.n_plaquettes - complex_.n_cubes
    )
    euler_homology = result["b0"] - result["b1"] + result["b2"] - result["b3"]
    if euler_cells != euler_homology:
        raise AssertionError("Euler characteristic does not match homology")
    return result


def gauge_orbit_size(complex_: CubicalComplex) -> int:
    """Number of distinct link configurations in every local gauge orbit."""

    return 1 << homology_dimensions(complex_)["rank_d1"]


def flat_gauge_orbit_count(complex_: CubicalComplex) -> int:
    """Number of gauge-inequivalent configurations with every plaquette satisfied."""

    return 1 << homology_dimensions(complex_)["b1"]


def gauge_satisfied_plaquette_dos(complex_: CubicalComplex) -> tuple[int, ...]:
    """Brute-force the raw gauge density of states by enumerating all link variables.

    Entry ``b`` counts configurations with exactly ``b`` positive plaquette products.  A set bit
    represents ``U_l = -1``.  The implementation deliberately enumerates all ``2**n_links`` raw
    configurations, including gauge-equivalent ones, and refuses sizes beyond 30 links.
    """

    n_links = complex_.n_links
    if n_links > 30:
        raise ValueError(f"brute-force gauge enumeration refuses {n_links} links > 30")
    n_plaquettes = complex_.n_plaquettes
    counts = np.zeros(n_plaquettes + 1, dtype=np.int64)
    n_states = 1 << n_links
    state_dtype = np.uint32

    for start in range(0, n_states, _ENUMERATION_CHUNK):
        stop = min(start + _ENUMERATION_CHUNK, n_states)
        states = np.arange(start, stop, dtype=state_dtype)
        frustrated = np.zeros(stop - start, dtype=np.int32)
        for boundary in complex_.plaquette_boundaries:
            parity = np.bitwise_count(states & state_dtype(boundary)) & np.uint8(1)
            frustrated += parity.astype(np.int32)
        satisfied = n_plaquettes - frustrated
        counts += np.bincount(satisfied, minlength=n_plaquettes + 1)

    if int(counts.sum()) != n_states:
        raise AssertionError("gauge density of states does not sum to 2**n_links")
    return tuple(int(value) for value in counts)


def _poly_mul(left: Sequence[int], right: Sequence[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        if left_coefficient:
            for right_degree, right_coefficient in enumerate(right):
                if right_coefficient:
                    result[left_degree + right_degree] += (
                        int(left_coefficient) * int(right_coefficient)
                    )
    return result


def _signed_binomial_transform(
    counts: Sequence[int], normalisation_power: int, positive_power_is_reversed: bool
) -> tuple[int, ...]:
    degree = len(counts) - 1
    positive_powers: list[list[int]] = [[1]]
    negative_powers: list[list[int]] = [[1]]
    for _ in range(degree):
        positive_powers.append(_poly_mul(positive_powers[-1], (1, 1)))
        negative_powers.append(_poly_mul(negative_powers[-1], (1, -1)))

    accumulator = [0] * (degree + 1)
    for index, count in enumerate(counts):
        if not count:
            continue
        positive_exponent = degree - index if positive_power_is_reversed else index
        negative_exponent = index if positive_power_is_reversed else degree - index
        term = _poly_mul(
            positive_powers[positive_exponent], negative_powers[negative_exponent]
        )
        for term_degree, coefficient in enumerate(term):
            accumulator[term_degree] += int(count) * coefficient

    denominator = 1 << normalisation_power
    result: list[int] = []
    for coefficient in accumulator:
        if coefficient % denominator:
            raise AssertionError("binomial transform has a non-integral coefficient")
        result.append(coefficient // denominator)
    return tuple(result)


def surface_polynomial_from_gauge_dos(
    complex_: CubicalComplex, dos: Sequence[int] | None = None
) -> tuple[int, ...]:
    """Recover ``S(v) = sum_{closed surfaces} v**area`` from the raw gauge DOS."""

    counts = tuple(dos) if dos is not None else gauge_satisfied_plaquette_dos(complex_)
    if len(counts) != complex_.n_plaquettes + 1:
        raise ValueError("gauge DOS length does not match the number of plaquettes")
    if sum(counts) != 1 << complex_.n_links:
        raise ValueError("gauge DOS does not count all link configurations")
    result = _signed_binomial_transform(
        counts, normalisation_power=complex_.n_links, positive_power_is_reversed=False
    )
    if result[0] != 1:
        raise AssertionError("the empty closed surface must have coefficient one")
    return result


def closed_surface_polynomial(complex_: CubicalComplex) -> tuple[int, ...]:
    """Independently enumerate ``ker(d2)`` from an exact GF(2) nullspace basis."""

    cycle_basis = _kernel_basis_of_columns(complex_.plaquette_boundaries)
    if len(cycle_basis) > 30:
        raise ValueError(f"closed-surface enumeration refuses dimension {len(cycle_basis)} > 30")
    counts = [0] * (complex_.n_plaquettes + 1)
    surface = 0
    previous_gray = 0
    for index in range(1 << len(cycle_basis)):
        gray = index ^ (index >> 1)
        if index:
            changed = gray ^ previous_gray
            basis_index = changed.bit_length() - 1
            surface ^= cycle_basis[basis_index]
        counts[surface.bit_count()] += 1
        previous_gray = gray
    return tuple(counts)


def _algebraic_surface_sector_basis(complex_: CubicalComplex) -> tuple[int, ...]:
    image_span: dict[int, int] = {}
    for boundary in complex_.cube_boundaries:
        _gf2_add(boundary, image_span)

    quotient_basis: list[int] = []
    for cycle in _kernel_basis_of_columns(complex_.plaquette_boundaries):
        remainder = _gf2_reduce(cycle, image_span)
        if remainder:
            quotient_basis.append(remainder)
            _gf2_add(remainder, image_span)
    return tuple(quotient_basis)


def surface_sector_basis(complex_: CubicalComplex) -> tuple[int, ...]:
    """Return closed-surface representatives generating ``H2``.

    For rectangular product boundary conditions the preferred representatives are coordinate
    two-tori: one plaquette sheet for each pair of periodic axes.  An exact algebraic quotient
    computation validates that choice and is retained as a general fallback.
    """

    expected_dimension = homology_dimensions(complex_)["b2"]
    image_span: dict[int, int] = {}
    for boundary in complex_.cube_boundaries:
        _gf2_add(boundary, image_span)

    preferred: list[int] = []
    for axis_a, axis_b in combinations(range(3), 2):
        if not (complex_.periodic[axis_a] and complex_.periodic[axis_b]):
            continue
        normal_axis = next(axis for axis in range(3) if axis not in (axis_a, axis_b))
        sheet = 0
        for face_index, (coord, face_axis_a, face_axis_b) in enumerate(complex_.plaquettes):
            if (
                (face_axis_a, face_axis_b) == (axis_a, axis_b)
                and coord[normal_axis] == 0
            ):
                sheet ^= 1 << face_index
        boundary = 0
        remaining = sheet
        while remaining:
            least_bit = remaining & -remaining
            boundary ^= complex_.plaquette_boundaries[least_bit.bit_length() - 1]
            remaining ^= least_bit
        if boundary:
            raise AssertionError("coordinate sheet is not a closed surface")
        if _gf2_reduce(sheet, image_span):
            preferred.append(sheet)
            _gf2_add(sheet, image_span)

    if len(preferred) == expected_dimension:
        return tuple(preferred)

    algebraic = _algebraic_surface_sector_basis(complex_)
    if len(algebraic) != expected_dimension:
        raise AssertionError("failed to construct a basis of H2")
    return algebraic


def surface_sector_representatives(complex_: CubicalComplex) -> tuple[int, ...]:
    """Return one packed closed surface from every class in ``H2``."""

    basis = surface_sector_basis(complex_)
    representatives: list[int] = []
    for sector in range(1 << len(basis)):
        representative = 0
        for basis_index, basis_surface in enumerate(basis):
            if (sector >> basis_index) & 1:
                representative ^= basis_surface
        representatives.append(representative)
    return tuple(representatives)


def dual_ising_sector_polynomials(
    complex_: CubicalComplex, representatives: Sequence[int] | None = None
) -> tuple[tuple[int, ...], ...]:
    """Enumerate broken-bond polynomials of every signed dual Ising sector.

    A dual spin lives on each cube.  Flipping cube ``c`` toggles the six plaquette interactions in
    ``d3(c)``.  On an open boundary, a plaquette incident on one cube is therefore a coupling to a
    fixed exterior ``+`` spin.  ``representative`` supplies the negative bonds for a topological
    sector.  The returned coefficient at degree ``q`` counts spin configurations with ``q``
    unsatisfied signed interactions.
    """

    if complex_.n_cubes > 30:
        raise ValueError(f"dual Ising enumeration refuses {complex_.n_cubes} spins > 30")
    sectors = (
        tuple(int(value) for value in representatives)
        if representatives is not None
        else surface_sector_representatives(complex_)
    )

    polynomials: list[tuple[int, ...]] = []
    for representative in sectors:
        if representative < 0 or representative >> complex_.n_plaquettes:
            raise ValueError("surface representative has an out-of-range plaquette bit")
        link_boundary = 0
        remaining = representative
        while remaining:
            least_bit = remaining & -remaining
            link_boundary ^= complex_.plaquette_boundaries[least_bit.bit_length() - 1]
            remaining ^= least_bit
        if link_boundary:
            raise ValueError("dual Ising sector representative is not a closed surface")

        counts = [0] * (complex_.n_plaquettes + 1)
        unsatisfied_surface = representative
        previous_gray = 0
        for spin_index in range(1 << complex_.n_cubes):
            gray = spin_index ^ (spin_index >> 1)
            if spin_index:
                changed = gray ^ previous_gray
                cube_index = changed.bit_length() - 1
                unsatisfied_surface ^= complex_.cube_boundaries[cube_index]
            counts[unsatisfied_surface.bit_count()] += 1
            previous_gray = gray
        if sum(counts) != 1 << complex_.n_cubes:
            raise AssertionError("dual Ising DOS does not sum to 2**n_cubes")
        polynomials.append(tuple(counts))
    return tuple(polynomials)


def _row_horizontal_cost(state: int, width: int, twist_x: int) -> int:
    cost = 0
    for x_coord in range(width):
        next_coord = (x_coord + 1) % width
        broken = ((state >> x_coord) ^ (state >> next_coord)) & 1
        if x_coord == width - 1:
            broken ^= twist_x
        cost += broken
    return cost


def ising_torus_sector_polynomials_2d(shape: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Return exact broken-bond polynomials for all four 2D torus twist sectors.

    Sector order is ``(tx, ty) = (0,0), (1,0), (0,1), (1,1)``.  ``tx`` changes the sign of the
    horizontal wrap bonds and ``ty`` changes the sign of the vertical wrap bonds.  A row transfer
    calculation uses exact ``numpy.int64`` counts; ``N < 63`` is enforced so overflow is impossible.
    """

    shape_tuple = tuple(int(length) for length in shape)
    if len(shape_tuple) != 2 or any(length < 2 for length in shape_tuple):
        raise ValueError("a 2D torus shape must contain two side lengths >= 2")
    width, height = shape_tuple
    n_sites = width * height
    if n_sites >= 63:
        raise ValueError("exact int64 transfer enumeration requires fewer than 63 spins")
    n_edges = 2 * n_sites
    n_row_states = 1 << width
    full_row = n_row_states - 1
    vertical_cost = np.empty((n_row_states, n_row_states), dtype=np.int16)
    for previous in range(n_row_states):
        for current in range(n_row_states):
            vertical_cost[previous, current] = (previous ^ current).bit_count()

    sector_polynomials: list[tuple[int, ...]] = []
    for sector_index in range(4):
        twist_x = sector_index & 1
        twist_y = (sector_index >> 1) & 1
        horizontal_cost = [
            _row_horizontal_cost(state, width, twist_x) for state in range(n_row_states)
        ]

        paths = np.zeros((n_row_states, n_row_states, n_edges + 1), dtype=np.int64)
        for first_state in range(n_row_states):
            paths[first_state, first_state, horizontal_cost[first_state]] = 1

        for _ in range(1, height):
            extended = np.zeros_like(paths)
            for previous_state in range(n_row_states):
                source = paths[:, previous_state, :]
                for current_state in range(n_row_states):
                    shift = int(vertical_cost[previous_state, current_state]) + horizontal_cost[
                        current_state
                    ]
                    if shift <= n_edges:
                        extended[:, current_state, shift:] += source[:, : n_edges + 1 - shift]
            paths = extended

        counts = np.zeros(n_edges + 1, dtype=np.int64)
        seam_mask = full_row if twist_y else 0
        for first_state in range(n_row_states):
            for last_state in range(n_row_states):
                shift = (last_state ^ first_state ^ seam_mask).bit_count()
                counts[shift:] += paths[first_state, last_state, : n_edges + 1 - shift]
        if int(counts.sum()) != 1 << n_sites:
            raise AssertionError("2D sector polynomial does not count all spin configurations")
        sector_polynomials.append(tuple(int(value) for value in counts))
    return tuple(sector_polynomials)


def signed_high_temperature_polynomial(
    broken_bond_dos: Sequence[int], n_sites: int
) -> tuple[int, ...]:
    """Convert a signed Ising broken-bond DOS to its high-temperature polynomial.

    If ``q`` counts negative signed interactions, this computes exactly
    ``2**(-N) sum_q g[q] (1+v)**(E-q) (1-v)**q``.  Coefficients may be negative in a twisted
    sector, but are always integers.
    """

    counts = tuple(int(value) for value in broken_bond_dos)
    if n_sites < 0 or sum(counts) != 1 << n_sites:
        raise ValueError("broken-bond DOS does not count 2**n_sites configurations")
    return _signed_binomial_transform(
        counts, normalisation_power=n_sites, positive_power_is_reversed=True
    )


def kw_character(alpha: int, beta: int) -> int:
    """Character in the full four-sector 2D torus Kramers--Wannier transform.

    Sector bits use ``tx + 2*ty``.  Primal and dual seams intersect crosswise, hence the exponent
    ``alpha_x*beta_y + alpha_y*beta_x`` rather than the ordinary bit dot product.
    """

    if alpha not in range(4) or beta not in range(4):
        raise ValueError("Kramers--Wannier sector labels must lie in range(4)")
    alpha_x, alpha_y = alpha & 1, (alpha >> 1) & 1
    beta_x, beta_y = beta & 1, (beta >> 1) & 1
    parity = (alpha_x * beta_y + alpha_y * beta_x) & 1
    return -1 if parity else 1
