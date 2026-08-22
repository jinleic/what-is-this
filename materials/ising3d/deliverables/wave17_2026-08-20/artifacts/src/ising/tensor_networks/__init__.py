"""Exact finite-field and tensor-network tests of Ising-layer integrability.

The local-charge calculation uses translation-orbit sums of Pauli strings.  A Pauli density is
admitted at range ``r`` precisely when its non-identity sites are contained in some connected
lattice animal with at most ``r`` vertices.  Translation is the only quotient: rotations and
reflections remain independent coefficients.  Commutators are assembled with the packed
symplectic convention from :mod:`ising.clifford` and their ranks are computed exactly in prime
fields.

The transfer matrix is built without transcendental numbers.  If ``x = exp(-2 K) = p/q``, then
``tanh(K*) = x`` and ``tanh(K) = (q-p)/(q+p)``.  Dropping positive scalar factors gives an integer
matrix with the same eigenspaces as ``exp(K* sum X) exp(K sum ZZ)``.

The rank-six Ising tensor and the constant tetrahedron equation are treated separately.  This is
a test of one particular reshaping of the site tensor as an ``8 x 8`` operator; it is not a
classification of tetrahedron-equation solutions.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from math import gcd, prod
from typing import Iterable, Sequence

import numpy as np

from ising.clifford import pauli_x, zz

__all__ = [
    "exact_transfer_commutant",
    "ising_tetrahedron_analysis",
    "scaled_ising_transfer_matrix",
    "scan_local_conserved_charges",
    "stable_transfer_commutant",
]


@dataclass
class _PeriodicGeometry:
    shape: tuple[int, ...]
    coordinates: tuple[tuple[int, ...], ...]
    translations: tuple[tuple[int, ...], ...]
    translation_maps: tuple[tuple[int, ...], ...]
    neighbours: tuple[tuple[int, ...], ...]
    bonds: tuple[tuple[int, int], ...]

    @property
    def n_sites(self) -> int:
        return len(self.coordinates)


def _site_index(coord: Sequence[int], shape: Sequence[int]) -> int:
    index = 0
    for value, length in zip(coord, shape):
        index = index * length + value
    return index


def _periodic_geometry(shape: Sequence[int]) -> _PeriodicGeometry:
    shape = tuple(int(length) for length in shape)
    if not shape or any(length < 1 for length in shape):
        raise ValueError("shape must contain positive side lengths")
    coordinates = tuple(product(*(range(length) for length in shape)))
    translations = coordinates
    maps = []
    for delta in translations:
        maps.append(
            tuple(
                _site_index(tuple((coord[axis] + delta[axis]) % shape[axis] for axis in range(len(shape))), shape)
                for coord in coordinates
            )
        )

    neighbour_sets = [set() for _ in coordinates]
    bonds = []
    for coord in coordinates:
        i = _site_index(coord, shape)
        for axis, length in enumerate(shape):
            if length == 1:
                continue
            forward = list(coord)
            forward[axis] = (forward[axis] + 1) % length
            j = _site_index(forward, shape)
            bonds.append((i, j))
            neighbour_sets[i].add(j)
            neighbour_sets[j].add(i)
    return _PeriodicGeometry(
        shape=shape,
        coordinates=coordinates,
        translations=translations,
        translation_maps=tuple(maps),
        neighbours=tuple(tuple(sorted(values)) for values in neighbour_sets),
        bonds=tuple(bonds),
    )


def _iter_set_bits(value: int) -> Iterable[int]:
    while value:
        low = value & -value
        yield low.bit_length() - 1
        value ^= low


def _translated_site_mask(mask: int, mapping: Sequence[int]) -> int:
    translated = 0
    for site in _iter_set_bits(mask):
        translated |= 1 << mapping[site]
    return translated


def _canonical_site_mask(mask: int, geometry: _PeriodicGeometry) -> int:
    return min(_translated_site_mask(mask, mapping) for mapping in geometry.translation_maps)


def _translated_pauli(value: int, mapping: Sequence[int], n_sites: int) -> int:
    site_mask = (1 << n_sites) - 1
    x_part = value & site_mask
    z_part = value >> n_sites
    translated_x = _translated_site_mask(x_part, mapping)
    translated_z = _translated_site_mask(z_part, mapping)
    return translated_x | (translated_z << n_sites)


def _canonical_pauli(value: int, geometry: _PeriodicGeometry, cache: dict[int, int]) -> int:
    cached = cache.get(value)
    if cached is not None:
        return cached
    canonical = min(
        _translated_pauli(value, mapping, geometry.n_sites) for mapping in geometry.translation_maps
    )
    cache[value] = canonical
    return canonical


def _connected_animals(geometry: _PeriodicGeometry, max_sites: int) -> dict[int, set[int]]:
    animals: dict[int, set[int]] = {1: {_canonical_site_mask(1, geometry)}}
    for size in range(1, max_sites):
        enlarged: set[int] = set()
        for animal in animals[size]:
            frontier: set[int] = set()
            for site in _iter_set_bits(animal):
                frontier.update(geometry.neighbours[site])
            for site in frontier:
                if not (animal >> site) & 1:
                    enlarged.add(_canonical_site_mask(animal | (1 << site), geometry))
        animals[size + 1] = enlarged
    return animals


def _local_pauli_orbits(geometry: _PeriodicGeometry, max_sites: int) -> tuple[list[int], dict[int, int], dict[int, int]]:
    """Return canonical Pauli representatives and their minimum connected-cover sizes."""

    animals = _connected_animals(geometry, max_sites)
    canonical_cache: dict[int, int] = {0: 0}
    minimum_cover: dict[int, int] = {0: 0}
    n_sites = geometry.n_sites
    for size in range(1, max_sites + 1):
        for animal in animals[size]:
            sites = tuple(_iter_set_bits(animal))
            for labels in product(range(4), repeat=size):
                if not any(labels):
                    continue
                value = 0
                for site, label in zip(sites, labels):
                    if label & 1:
                        value |= 1 << site
                    if label & 2:
                        value |= 1 << (n_sites + site)
                canonical = _canonical_pauli(value, geometry, canonical_cache)
                minimum_cover.setdefault(canonical, size)
    representatives = sorted(minimum_cover, key=lambda value: (minimum_cover[value], value))
    return representatives, minimum_cover, canonical_cache


def _hamiltonian_terms(
    geometry: _PeriodicGeometry,
    field_numerator: int,
    field_denominator: int,
) -> dict[int, int]:
    """Integer coefficients of ``field_denominator * H``; an overall sign is immaterial."""

    if field_denominator <= 0 or field_numerator == 0:
        raise ValueError("the transverse field must be a nonzero rational with positive denominator")
    terms: dict[int, int] = defaultdict(int)
    n_sites = geometry.n_sites
    for site in range(n_sites):
        terms[pauli_x(n_sites, site)] -= field_numerator
    for i, j in geometry.bonds:
        terms[zz(n_sites, i, j)] -= field_denominator
    return dict(terms)


def _commutator_columns(
    representatives: Sequence[int],
    geometry: _PeriodicGeometry,
    hamiltonian_terms: dict[int, int],
    canonical_cache: dict[int, int],
) -> list[dict[int, int]]:
    """Columns of ``Q -> [Q,H]`` in full translation-orbit-sum bases."""

    n_sites = geometry.n_sites
    site_mask = (1 << n_sites) - 1
    columns = []
    for value in representatives:
        x_value, z_value = value & site_mask, value >> n_sites
        column: dict[int, int] = defaultdict(int)
        for term, coefficient in hamiltonian_terms.items():
            x_term, z_term = term & site_mask, term >> n_sites
            left_phase = (z_value & x_term).bit_count() & 1
            right_phase = (z_term & x_value).bit_count() & 1
            if left_phase == right_phase:
                continue
            structure_constant = 2 if left_phase == 0 else -2
            output = _canonical_pauli(value ^ term, geometry, canonical_cache)
            column[output] += structure_constant * coefficient
        columns.append({row: coefficient for row, coefficient in column.items() if coefficient})
    return columns


def _sparse_rank_mod_prime(columns: Sequence[dict[int, int]], prime: int) -> int:
    """Exact column rank over ``F_prime`` with a sparse rare-row pivot order."""

    if prime <= 2:
        raise ValueError("prime must be an odd prime")
    frequencies = Counter(row for column in columns for row in column)
    row_order = {row: order for order, row in enumerate(sorted(frequencies, key=lambda row: (frequencies[row], row)))}
    basis: dict[int, dict[int, int]] = {}
    rank = 0
    for original in sorted(columns, key=len):
        vector = {row: coefficient % prime for row, coefficient in original.items() if coefficient % prime}
        while vector:
            lead = min(vector, key=row_order.__getitem__)
            prior = basis.get(lead)
            if prior is None:
                inverse = pow(vector[lead], prime - 2, prime)
                vector = {row: (coefficient * inverse) % prime for row, coefficient in vector.items()}
                basis[lead] = vector
                rank += 1
                break
            factor = vector[lead]
            for row, coefficient in prior.items():
                updated = (vector.get(row, 0) - factor * coefficient) % prime
                if updated:
                    vector[row] = updated
                else:
                    vector.pop(row, None)
    return rank


def scan_local_conserved_charges(
    *,
    dimension: int,
    max_support: int,
    linear_size: int | None = None,
    field_numerator: int = 2,
    field_denominator: int = 1,
    primes: Sequence[int] = (2_147_483_647, 2_147_483_629),
) -> dict:
    """Search translation-invariant local charges on a large periodic chain or square torus.

    The finite torus is only a faithful container for the local equations when it is wider than
    twice the largest connected cover.  This function enforces ``linear_size >= 2*r+1``.  Returned
    dimensions are exact dimensions over each listed finite field; agreement of two fields guards
    against, but is not presented as a proof excluding, an unlucky characteristic.
    """

    if dimension not in (1, 2):
        raise ValueError("only one- and two-dimensional layers are supported")
    if max_support < 2:
        raise ValueError("max_support must be at least 2")
    minimum_size = 2 * max_support + 1
    if linear_size is None:
        linear_size = minimum_size
    if linear_size < minimum_size:
        raise ValueError(f"linear_size must be at least 2*max_support+1 = {minimum_size}")
    geometry = _periodic_geometry((linear_size,) * dimension)
    representatives, minimum_cover, canonical_cache = _local_pauli_orbits(geometry, max_support)
    hamiltonian = _hamiltonian_terms(geometry, field_numerator, field_denominator)
    all_columns = _commutator_columns(representatives, geometry, hamiltonian, canonical_cache)

    rows = []
    for support in range(2, max_support + 1):
        indices = [index for index, value in enumerate(representatives) if minimum_cover[value] <= support]
        columns = [all_columns[index] for index in indices]
        ranks = [_sparse_rank_mod_prime(columns, int(prime)) for prime in primes]
        nullities = [len(columns) - rank for rank in ranks]
        if len(set(nullities)) != 1:
            raise ArithmeticError(f"finite-field nullities disagree at support {support}: {nullities}")
        kernel_dimension = nullities[0]
        trivial_dimension = 2  # full translation sum of I, and H itself
        if kernel_dimension < trivial_dimension:
            raise ArithmeticError("the computed kernel omitted identity or the Hamiltonian")
        output_rows = set().union(*(column.keys() for column in columns)) if columns else set()
        rows.append(
            {
                "max_support": support,
                "ansatz_dimension": len(columns),
                "commutator_target_dimension": len(output_rows),
                "rank_by_prime": ranks,
                "kernel_dimension_by_prime": nullities,
                "kernel_dimension": kernel_dimension,
                "trivial_dimension": trivial_dimension,
                "nontrivial_dimension": kernel_dimension - trivial_dimension,
            }
        )
    return {
        "dimension": dimension,
        "periodic_shape": list(geometry.shape),
        "field_g": f"{field_numerator}/{field_denominator}",
        "primes": [int(prime) for prime in primes],
        "ansatz_definition": (
            "one full translation-orbit sum for every Pauli string whose non-identity support is "
            "contained in a connected lattice animal of at most r sites; rotations and reflections "
            "are not quotiented"
        ),
        "trivial_definition": (
            "the two-dimensional span of the full-orbit identity density and H; global spin flip "
            "and lattice permutation operators are nonlocal and absent from the r-local ansatz"
        ),
        "rows": rows,
    }


def _open_bonds(shape: Sequence[int]) -> tuple[tuple[int, int], ...]:
    shape = tuple(shape)
    coordinates = tuple(product(*(range(length) for length in shape)))
    bonds = []
    for coord in coordinates:
        i = _site_index(coord, shape)
        for axis, length in enumerate(shape):
            if coord[axis] + 1 < length:
                neighbour = list(coord)
                neighbour[axis] += 1
                bonds.append((i, _site_index(neighbour, shape)))
    return tuple(bonds)


def scaled_ising_transfer_matrix(
    shape: Sequence[int],
    *,
    x_numerator: int,
    x_denominator: int,
) -> tuple[np.ndarray, dict]:
    """Return an exactly integer-scaled open-layer Ising transfer matrix.

    The omitted scalar is ``cosh(K*)**n * cosh(K)**n_b`` together with the integer denominators.
    It cannot affect eigenvalue multiplicities or the commutant dimension.
    """

    shape = tuple(int(length) for length in shape)
    if not shape or len(shape) > 2 or any(length < 1 for length in shape):
        raise ValueError("shape must describe a nonempty one- or two-dimensional layer")
    if not (0 < x_numerator < x_denominator):
        raise ValueError("x_numerator/x_denominator must lie strictly between zero and one")
    common = gcd(x_numerator, x_denominator)
    p, q = x_numerator // common, x_denominator // common
    v_numerator, v_denominator = q - p, q + p
    common_v = gcd(v_numerator, v_denominator)
    v_numerator //= common_v
    v_denominator //= common_v

    n_sites = prod(shape)
    dimension = 1 << n_sites
    bonds = _open_bonds(shape)
    diagonal = []
    for state in range(dimension):
        weight = 1
        for i, j in bonds:
            equal = ((state >> i) & 1) == ((state >> j) & 1)
            weight *= v_denominator + (v_numerator if equal else -v_numerator)
        diagonal.append(weight)

    matrix = np.empty((dimension, dimension), dtype=object)
    for row in range(dimension):
        for column in range(dimension):
            distance = (row ^ column).bit_count()
            matrix[row, column] = q ** (n_sites - distance) * p**distance * diagonal[column]
    metadata = {
        "shape": list(shape),
        "matrix_dimension": dimension,
        "sites": n_sites,
        "open_bonds": len(bonds),
        "x": f"{p}/{q}",
        "v": f"{v_numerator}/{v_denominator}",
        "integer_scaling": {
            "transverse_factor_denominator": q**n_sites,
            "interaction_factor_denominator": v_denominator ** len(bonds),
        },
    }
    return matrix, metadata


def exact_transfer_commutant(
    shape: Sequence[int],
    *,
    x_numerator: int,
    x_denominator: int,
) -> dict:
    """Compute the full matrix commutant dimension from an exact characteristic polynomial.

    This routine is intended for matrices through dimension 64.  Since the transfer matrix is
    similar to a real symmetric positive-definite matrix, it is diagonalizable.  If its
    characteristic polynomial has square-free decomposition ``prod f_e**e``, every root of
    ``f_e`` has multiplicity ``e`` and the commutant dimension is ``sum deg(f_e)*e**2``.
    """

    import sympy as sp

    matrix, metadata = scaled_ising_transfer_matrix(
        shape, x_numerator=x_numerator, x_denominator=x_denominator
    )
    if metadata["matrix_dimension"] > 64:
        raise ValueError("exact characteristic-polynomial mode is capped at matrix dimension 64")
    variable = sp.Symbol("lambda")
    characteristic = sp.Matrix(matrix.tolist()).charpoly(variable).as_poly()
    _, square_free = sp.sqf_list(characteristic)
    factors = [
        {"degree": int(factor.degree()), "eigenvalue_multiplicity": int(multiplicity)}
        for factor, multiplicity in square_free
    ]
    commutant_dimension = sum(row["degree"] * row["eigenvalue_multiplicity"] ** 2 for row in factors)
    return {
        **metadata,
        "method": "exact_characteristic_polynomial",
        "arithmetic": "integer matrix and Q[lambda] square-free factorization",
        "diagonalizable": True,
        "diagonalizability_reason": (
            "A D is similar to sqrt(D) A sqrt(D), which is real symmetric positive definite for 0<x<1"
        ),
        "characteristic_polynomial_degree": int(characteristic.degree()),
        "square_free_factors": factors,
        "distinct_eigenvalues": sum(row["degree"] for row in factors),
        "commutant_dimension": int(commutant_dimension),
    }


def _state_permutation(shape: tuple[int, ...], axis: int | None) -> tuple[int, ...]:
    """Permutation of computational states: spin flip if axis is None, else coordinate reflection."""

    n_sites = prod(shape)
    dimension = 1 << n_sites
    if axis is None:
        return tuple(state ^ (dimension - 1) for state in range(dimension))
    coordinates = tuple(product(*(range(length) for length in shape)))
    site_map = []
    for coord in coordinates:
        reflected = list(coord)
        reflected[axis] = shape[axis] - 1 - reflected[axis]
        site_map.append(_site_index(reflected, shape))
    permutation = []
    for state in range(dimension):
        mapped = 0
        for site in _iter_set_bits(state):
            mapped |= 1 << site_map[site]
        permutation.append(mapped)
    return tuple(permutation)


def _symmetry_sector_basis(
    dimension: int,
    generators: Sequence[Sequence[int]],
    signs: Sequence[int],
) -> list[dict[int, int]]:
    """Orthonormal-after-normalization signed orbit vectors for commuting involutions."""

    group_elements = []
    for powers in product((0, 1), repeat=len(generators)):
        permutation = tuple(range(dimension))
        character = 1
        for power, generator, sign in zip(powers, generators, signs):
            if power:
                permutation = tuple(generator[permutation[state]] for state in range(dimension))
                character *= sign
        group_elements.append((permutation, character))

    seen: set[int] = set()
    basis = []
    for state in range(dimension):
        orbit = {permutation[state] for permutation, _ in group_elements}
        if state in seen:
            continue
        seen.update(orbit)
        coefficients: dict[int, int] = defaultdict(int)
        for permutation, character in group_elements:
            coefficients[permutation[state]] += character
        coefficients = {index: coefficient for index, coefficient in coefficients.items() if coefficient}
        if coefficients:
            basis.append(coefficients)
    return basis


def _cluster_eigenvalues(values, tolerance):
    values = sorted(values)
    clusters = []
    for value in values:
        if not clusters or abs(value - clusters[-1][-1]) > tolerance * max(1, abs(value), abs(clusters[-1][-1])):
            clusters.append([value])
        else:
            clusters[-1].append(value)
    return clusters


def _multiprecision_commutant_once(
    shape: tuple[int, ...],
    x_numerator: int,
    x_denominator: int,
    dps: int,
    tolerance_exponent: int,
) -> dict:
    import mpmath as mp

    matrix, metadata = scaled_ising_transfer_matrix(
        shape, x_numerator=x_numerator, x_denominator=x_denominator
    )
    mp.mp.dps = dps
    dimension = metadata["matrix_dimension"]
    # Reconstructing from metadata-independent integer factors avoids taking square roots of the
    # nonsymmetric matrix itself.  A D is similar to S=sqrt(D) A sqrt(D).
    n_sites = prod(shape)
    p = x_numerator // gcd(x_numerator, x_denominator)
    q = x_denominator // gcd(x_numerator, x_denominator)
    a_entries = [
        [
            q ** (n_sites - (row ^ column).bit_count()) * p ** (row ^ column).bit_count()
            for column in range(dimension)
        ]
        for row in range(dimension)
    ]
    # matrix[0,j] = A[0,j] D[j], so the division is exact.
    diagonal = [int(matrix[0, column]) // a_entries[0][column] for column in range(dimension)]
    scale = mp.mpf(q**n_sites * max(diagonal))
    sqrt_diagonal = [mp.sqrt(value) for value in diagonal]

    generators = [_state_permutation(shape, None)]
    generators.extend(_state_permutation(shape, axis) for axis in range(len(shape)))
    sector_rows = []
    all_values = []
    for signs in product((-1, 1), repeat=len(generators)):
        basis = _symmetry_sector_basis(dimension, generators, signs)
        block_dimension = len(basis)
        if not basis:
            continue
        block = mp.matrix(block_dimension)
        norms = [mp.sqrt(sum(coefficient * coefficient for coefficient in vector.values())) for vector in basis]
        for row_index, row_vector in enumerate(basis):
            for column_index in range(row_index + 1):
                column_vector = basis[column_index]
                total = mp.mpf("0")
                for row_state, row_coefficient in row_vector.items():
                    sqrt_row = sqrt_diagonal[row_state]
                    for column_state, column_coefficient in column_vector.items():
                        total += (
                            row_coefficient
                            * column_coefficient
                            * sqrt_row
                            * a_entries[row_state][column_state]
                            * sqrt_diagonal[column_state]
                        )
                total /= norms[row_index] * norms[column_index] * scale
                block[row_index, column_index] = total
                block[column_index, row_index] = total
        eigenvalues = list(mp.eigsy(block, eigvals_only=True))
        all_values.extend(eigenvalues)
        sector_rows.append({"signs": list(signs), "dimension": block_dimension})

    tolerance = mp.mpf(10) ** (-tolerance_exponent)
    clusters = _cluster_eigenvalues(all_values, tolerance)
    multiplicities = Counter(len(cluster) for cluster in clusters)
    return {
        **metadata,
        "working_precision_decimal_digits": dps,
        "relative_cluster_tolerance": mp.nstr(tolerance, 8),
        "sector_generators": ["global_spin_flip"] + [f"reflection_axis_{axis}" for axis in range(len(shape))],
        "sector_dimensions": sector_rows,
        "distinct_eigenvalues": len(clusters),
        "multiplicity_histogram": {str(key): value for key, value in sorted(multiplicities.items())},
        "commutant_dimension": sum(len(cluster) ** 2 for cluster in clusters),
    }


def stable_transfer_commutant(
    shape: Sequence[int],
    *,
    x_numerator: int,
    x_denominator: int,
    precisions: Sequence[int] = (50, 80),
    tolerance_exponents: Sequence[int] = (30, 45),
) -> dict:
    """Determine spectral multiplicities twice with symmetry-blocked ``mpmath`` diagonalization."""

    shape = tuple(int(length) for length in shape)
    if len(precisions) != 2 or len(tolerance_exponents) != 2:
        raise ValueError("exactly two precision and tolerance settings are required")
    runs = [
        _multiprecision_commutant_once(
            shape,
            x_numerator,
            x_denominator,
            int(dps),
            int(tolerance_exponent),
        )
        for dps, tolerance_exponent in zip(precisions, tolerance_exponents)
    ]
    stable = all(
        run["commutant_dimension"] == runs[0]["commutant_dimension"]
        and run["multiplicity_histogram"] == runs[0]["multiplicity_histogram"]
        for run in runs[1:]
    )
    if not stable:
        raise ArithmeticError("spectral multiplicities were not stable across precision settings")
    return {
        **{key: value for key, value in runs[-1].items() if key not in {"working_precision_decimal_digits", "relative_cluster_tolerance"}},
        "method": "mpmath_symmetry_block_eigenspectra",
        "diagonalizable": True,
        "precision_stability_passed": stable,
        "precision_runs": [
            {
                "working_precision_decimal_digits": run["working_precision_decimal_digits"],
                "relative_cluster_tolerance": run["relative_cluster_tolerance"],
                "distinct_eigenvalues": run["distinct_eigenvalues"],
                "multiplicity_histogram": run["multiplicity_histogram"],
                "commutant_dimension": run["commutant_dimension"],
            }
            for run in runs
        ],
        "numerical_scope": (
            "The integer transfer matrix is exact; only the decision that two algebraic eigenvalues "
            "coincide is numerical.  Agreement at both recorded precisions is a stability check, not "
            "an interval-arithmetic proof of every nondegeneracy."
        ),
    }


def _apply_local_r(
    amplitudes: dict[tuple[int, int], int],
    positions: tuple[int, int, int],
) -> dict[tuple[int, int], int]:
    """Apply the parity Ising R; keys are ``(six_bit_state, power_of_w)``."""

    output: dict[tuple[int, int], int] = defaultdict(int)
    position_mask = sum(1 << position for position in positions)
    for (state, exponent), coefficient in amplitudes.items():
        local_input = sum(((state >> position) & 1) << local for local, position in enumerate(positions))
        input_weight = local_input.bit_count()
        for local_output in range(8):
            if (input_weight + local_output.bit_count()) & 1:
                continue
            new_state = state & ~position_mask
            for local, position in enumerate(positions):
                if (local_output >> local) & 1:
                    new_state |= 1 << position
            output[(new_state, exponent + input_weight + local_output.bit_count())] += coefficient
    return dict(output)


def _tetrahedron_residual_polynomials() -> dict[tuple[int, int], tuple[int, ...]]:
    positions = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))
    residuals: dict[tuple[int, int], tuple[int, ...]] = {}
    for input_state in range(64):
        left = {(input_state, 0): 1}
        for local_positions in reversed(positions):
            left = _apply_local_r(left, local_positions)
        right = {(input_state, 0): 1}
        for local_positions in positions:
            right = _apply_local_r(right, local_positions)
        by_output: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        for (output_state, exponent), coefficient in left.items():
            by_output[output_state][exponent] += coefficient
        for (output_state, exponent), coefficient in right.items():
            by_output[output_state][exponent] -= coefficient
        for output_state, polynomial in by_output.items():
            if not any(polynomial.values()):
                continue
            degree = max(polynomial)
            coefficients = tuple(polynomial.get(power, 0) for power in range(degree + 1))
            if any(coefficients):
                residuals[(output_state, input_state)] = coefficients
    return residuals


def _ghz_tetrahedron_residual_is_zero() -> bool:
    positions = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))

    def apply(amplitudes, local_positions):
        output = {}
        mask = sum(1 << position for position in local_positions)
        for state, coefficient in amplitudes.items():
            local = tuple((state >> position) & 1 for position in local_positions)
            if local == (0, 0, 0) or local == (1, 1, 1):
                output[state & ~mask | sum(bit << position for bit, position in zip(local, local_positions))] = coefficient
        return output

    for input_state in range(64):
        left = {input_state: 1}
        right = {input_state: 1}
        for local_positions in reversed(positions):
            left = apply(left, local_positions)
        for local_positions in positions:
            right = apply(right, local_positions)
        if left != right:
            return False
    return True


@lru_cache(maxsize=1)
def ising_tetrahedron_analysis() -> dict:
    """Specialize the constant tetrahedron equation to the isotropic bond-dimension-two tensor."""

    import sympy as sp

    w = sp.Symbol("w")
    residual_coefficients = _tetrahedron_residual_polynomials()
    polynomials = [
        sp.Poly(sum(coefficient * w**power for power, coefficient in enumerate(coefficients)), w, domain=sp.ZZ)
        for coefficients in residual_coefficients.values()
    ]
    common = polynomials[0]
    for polynomial in polynomials[1:]:
        common = sp.gcd(common, polynomial)
        if common.degree() == 0:
            break
    common = sp.Poly(common.monic(), w)
    factorization = sp.factor(common.as_expr())
    roots = sp.solve(common.as_expr(), w)
    physical_roots = [root for root in roots if root.is_real is True and root.is_positive is True and (1 - root).is_positive is True]
    unique_polynomials = {
        coefficients for coefficients in residual_coefficients.values()
    }
    unique_expressions = [
        sum(coefficient * w**power for power, coefficient in enumerate(coefficients))
        for coefficients in sorted(unique_polynomials)
    ]
    groebner = sp.groebner(unique_expressions, w, order="lex", domain=sp.QQ)
    groebner_generators = [sp.Poly(polynomial, w).monic() for polynomial in groebner.polys]
    if groebner_generators != [common]:
        raise ArithmeticError("the univariate Groebner basis disagreed with the polynomial gcd")

    ghz_passed = _ghz_tetrahedron_residual_is_zero()
    if not ghz_passed:
        raise ArithmeticError("the explicitly gauge-related GHZ tensor failed the tetrahedron equation")
    return {
        "local_tensor": {
            "bond_dimension": 2,
            "definition": "T[a,b,c,d,e,f] = sum_{sigma=+-1} product_l u_sigma[index_l], u_sigma=(1,sigma*w)",
            "entries": "T = 2*w^(a+b+c+d+e+f) for even index sum and 0 for odd index sum",
            "parameter_relation": "w^2 = v = tanh(K)",
            "partition_function_prefactor": "(cosh K)^number_of_bonds",
            "operator_reshaping": "R[out_x,out_y,out_z ; in_x,in_y,in_z] = T[in_x,in_y,in_z,out_x,out_y,out_z]",
        },
        "general_system": {
            "variables": 64,
            "equations": 4096,
            "degree": 4,
            "terms_per_side_before_collection": 64,
            "equation": "R_123 R_145 R_246 R_356 = R_356 R_246 R_145 R_123 on (C^2)^tensor_6",
        },
        "isotropic_substitution": {
            "nonzero_equations": len(residual_coefficients),
            "distinct_nonzero_polynomials": len(unique_polynomials),
            "maximum_terms_in_residual": max(sum(coefficient != 0 for coefficient in values) for values in residual_coefficients.values()),
            "common_gcd": str(factorization),
            "reduced_groebner_basis": [str(polynomial.as_expr()) for polynomial in groebner_generators],
            "common_roots": [str(root) for root in roots],
            "physical_parameter_interval": "0 < w < 1 (finite ferromagnetic K)",
            "physical_interval_solutions": [str(root) for root in physical_roots],
            "method": (
                "exact integer coefficient propagation; in Q[w] the polynomial gcd is the "
                "reduced univariate Groebner-basis generator"
            ),
        },
        "arbitrary_leg_gl2_gauge": {
            "matrix": "G = [[1/2, 1/(2*w)], [1/2, -1/(2*w)]]",
            "action": "apply G covariantly to each of the six tensor legs",
            "normal_form": "G^tensor_6 T = |000000> + |111111> (the six-leg GHZ tensor)",
            "satisfies_tetrahedron_equation": ghz_passed,
            "scope_warning": (
                "This non-orthogonal transformation changes the standard delta bond metric if used "
                "identically at both ends of every edge.  It preserves the partition function only "
                "when the inverse metric, or alternating inverse gauges, is retained; it therefore "
                "does not by itself furnish a commuting 3D-Ising transfer family."
            ),
        },
        "delta_preserving_identical_tensor_gauge": {
            "gauge_class": "the same real G on every half-edge with G^T G = I (up to an overall scalar)",
            "possible_for_finite_positive_K": False,
            "reason": (
                "such a gauge is an operator similarity and tetrahedron-equation satisfaction is "
                "invariant; the raw tensor fails for every 0<w<1 by the exact common-gcd test"
            ),
        },
        "not_established": (
            "No classification was made for alternating bipartite tensors, changed bond metrics, "
            "higher bond dimension, spectral-parameter-dependent R matrices, or other vertex/IRF "
            "representations of the 3D Ising model."
        ),
    }
