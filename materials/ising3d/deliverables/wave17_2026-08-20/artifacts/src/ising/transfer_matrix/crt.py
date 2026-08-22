"""Exact transfer-matrix propagation by modular arithmetic and CRT.

The public :mod:`ising.transfer_matrix` engine stores polynomial coefficients in
``numpy.int64`` and therefore rejects boxes with more than 62 sites.  The
propagation itself only shifts and adds coefficient arrays.  This module runs
that same propagation modulo distinct 31-bit primes and reconstructs the
non-negative integer coefficients with the Chinese remainder theorem.

For a box with ``N`` sites, every coefficient is at most ``2**N``.  The prime
product is therefore chosen strictly larger than ``2**N``; reconstruction in
``[0, product)`` is then the unique exact answer.  A single inter-layer
butterfly can increase a residue by at most ``2**ns``, where ``ns <= 22`` is the
cross-section guard inherited from the original engine.  Thus a complete
butterfly remains below ``2**53`` for primes below ``2**31`` and is safe in
signed int64 before the next modular reduction.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import prod
from typing import Iterator, Sequence

import numpy as np

from . import (
    _broken_in_layer,
    _shift_rows,
    _weighted_downcount,
    box_broken_bond_poly as _int64_box_broken_bond_poly,
    layer_bonds,
    torus_broken_bond_poly as _int64_torus_broken_bond_poly,
)

__all__ = [
    "CRTModuli",
    "box_broken_bond_poly",
    "estimate_box_peak_bytes",
    "hybrid_box_broken_bond_poly",
    "hybrid_torus_broken_bond_poly",
    "layer_bonds",
    "primes_for_sites",
    "torus_broken_bond_poly",
]

# Audited lift 22 -> 25 (wave 16, notes/ht_v24_x54_feasibility.md section 0): with primes
# below 2^31 a complete butterfly pass stays below 2^(31+25) = 2^56 < 2^63 - 1, so the lift
# cannot change any output value; it only admits cross-sections 23..25 (needed for the
# LT x^54 box (5,5,5)).  The dense int64 engine keeps its own ns <= 22 guard.
_MAX_CROSS_SECTION = 25
_FIRST_31_BIT_CANDIDATE = (1 << 31) - 1


@dataclass(frozen=True)
class CRTModuli(Sequence[int]):
    """The deterministic prime set used for one reconstruction."""

    primes: tuple[int, ...]
    product: int
    required_sites: int

    def __getitem__(self, index):
        return self.primes[index]

    def __len__(self) -> int:
        return len(self.primes)

    def __iter__(self) -> Iterator[int]:
        return iter(self.primes)

    @property
    def product_bits(self) -> int:
        return self.product.bit_length()


def _is_prime_32(value: int) -> bool:
    """Deterministic Miller--Rabin primality test on unsigned 32-bit inputs."""

    if value < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small_primes:
        if value == prime:
            return True
        if value % prime == 0:
            return False

    odd_part = value - 1
    powers_of_two = 0
    while odd_part % 2 == 0:
        powers_of_two += 1
        odd_part //= 2
    for base in (2, 3, 5, 7, 11, 13, 17):
        if base >= value:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(powers_of_two - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


@lru_cache(maxsize=None)
def _descending_primes(count: int) -> tuple[int, ...]:
    if count < 1:
        return ()
    primes: list[int] = []
    candidate = _FIRST_31_BIT_CANDIDATE
    while len(primes) < count:
        if _is_prime_32(candidate):
            primes.append(candidate)
        candidate -= 2
    return tuple(primes)


@lru_cache(maxsize=None)
def primes_for_sites(sites: int) -> CRTModuli:
    """Return enough distinct 31-bit primes that their product exceeds ``2**sites``."""

    sites = int(sites)
    if sites < 0:
        raise ValueError("sites must be non-negative")
    target = 1 << sites
    count = max(1, sites // 30)
    while True:
        primes = _descending_primes(count)
        modulus_product = prod(primes)
        if modulus_product > target:
            return CRTModuli(primes, modulus_product, sites)
        count += 1


def _shape_data(shape, *, periodic=None):
    box = tuple(int(side) for side in shape)
    if len(box) not in (2, 3):
        raise ValueError("dim must be 2 or 3")
    if any(side < 1 for side in box):
        raise ValueError("shape sides must be positive")
    cross = box[:-1]
    layers = box[-1]
    if periodic is None:
        per = (False,) * len(cross)
    else:
        supplied = tuple(bool(value) for value in periodic)
        if len(supplied) < len(cross):
            raise ValueError("periodic must specify every cross-section direction")
        per = supplied[: len(cross)]
    cross_sites = prod(cross)
    if cross_sites > _MAX_CROSS_SECTION:
        raise ValueError(
            f"cross-section {cross} has {cross_sites} sites: too large "
            f"(maximum {_MAX_CROSS_SECTION})"
        )
    return box, cross, layers, per, cross_sites, prod(box)


def _degree_limit(full_degree: int, max_degree: int | None) -> int:
    if max_degree is None:
        return full_degree
    limit = int(max_degree)
    if limit < 0:
        raise ValueError("max_degree must be non-negative")
    return min(full_degree, limit)


def _apply_interlayer_mod(
    vector: np.ndarray, cross_sites: int, modulus: int
) -> np.ndarray:
    """Apply one inter-layer operator and reduce once after its safe butterfly."""

    tail = vector.shape[1:]
    degree = vector.shape[-1] - 1
    transformed = vector
    for bit in range(cross_sites):
        high = 1 << (cross_sites - bit - 1)
        low = 1 << bit
        reshaped = transformed.reshape((high, 2, low) + tail)
        added = reshaped.copy()
        added[..., 1:] += reshaped[:, ::-1][..., :degree]
        transformed = added.reshape((1 << cross_sites,) + tail)
    np.remainder(transformed, modulus, out=transformed)
    return transformed


def _accumulate_crt(
    reconstructed: list[int] | None,
    current_product: int,
    residues: np.ndarray,
    prime: int,
) -> tuple[list[int], int]:
    flat_residues = residues.reshape(-1)
    if reconstructed is None:
        return [int(value) for value in flat_residues], prime
    if len(reconstructed) != len(flat_residues):
        raise ValueError("inconsistent residue-vector lengths")
    inverse = pow(current_product % prime, -1, prime)
    for index, residue in enumerate(flat_residues):
        correction = ((int(residue) - reconstructed[index]) % prime) * inverse % prime
        reconstructed[index] += current_product * correction
    return reconstructed, current_product * prime


def _trim(coefficients: list[int]) -> list[int]:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def _box_residues(
    shape,
    periodic,
    plus_boundary: bool,
    max_degree: int | None,
    modulus: int,
) -> tuple[np.ndarray, int, int]:
    box, cross, layers, per, cross_sites, total_sites = _shape_data(
        shape, periodic=periodic
    )
    in_layer_bonds = layer_bonds(cross, per)
    physical_bonds = len(in_layer_bonds) * layers + cross_sites * (layers - 1)

    inplane_coordination = 2 * (len(box) - 1)
    inplane_degree = [0] * cross_sites
    for left, right in in_layer_bonds:
        inplane_degree[left] += 1
        inplane_degree[right] += 1
    ghost_inplane = [inplane_coordination - degree for degree in inplane_degree]
    if not all(missing >= 0 for missing in ghost_inplane):
        raise AssertionError("negative ghost-bond count")
    ghost_bonds = (
        sum(ghost_inplane) * layers + 2 * cross_sites if plus_boundary else 0
    )
    full_degree = physical_bonds + ghost_bonds
    degree = _degree_limit(full_degree, max_degree)

    broken_layer = _broken_in_layer(cross_sites, in_layer_bonds)
    if plus_boundary:
        broken_ghost_inplane = _weighted_downcount(cross_sites, ghost_inplane)
        broken_ghost_transfer = _weighted_downcount(
            cross_sites, [1] * cross_sites
        )
    else:
        broken_ghost_inplane = np.zeros(1 << cross_sites, dtype=np.int64)
        broken_ghost_transfer = np.zeros(1 << cross_sites, dtype=np.int64)

    vector = np.zeros((1 << cross_sites, degree + 1), dtype=np.int64)
    vector[:, 0] = 1
    # A one-layer box touches both frozen transfer-direction faces.
    first_shift = (
        broken_layer
        + broken_ghost_inplane
        + broken_ghost_transfer
        + (broken_ghost_transfer if layers == 1 else 0)
    )
    vector = _shift_rows(vector, first_shift)
    for layer in range(1, layers):
        vector = _apply_interlayer_mod(vector, cross_sites, modulus)
        shift = broken_layer + broken_ghost_inplane
        if layer == layers - 1:
            shift = shift + broken_ghost_transfer
        vector = _shift_rows(vector, shift)
    total = vector.sum(axis=0, dtype=np.int64)
    np.remainder(total, modulus, out=total)
    return total, total_sites, full_degree


def box_broken_bond_poly(
    shape,
    periodic=None,
    plus_boundary: bool = False,
    *,
    max_degree: int | None = None,
) -> list[int]:
    """Return the exact broken-bond polynomial for an open-transfer box.

    The first three arguments have the same meaning as
    :func:`ising.transfer_matrix.box_broken_bond_poly`.  ``max_degree`` is an
    optional exact low-degree truncation; omitting it returns the full
    polynomial and checks that its coefficients sum to ``2**N``.
    """

    total_sites = prod(tuple(int(side) for side in shape))
    moduli = primes_for_sites(total_sites)
    reconstructed: list[int] | None = None
    current_product = 1
    full_degree = 0
    for prime in moduli:
        residues, checked_sites, full_degree = _box_residues(
            shape, periodic, bool(plus_boundary), max_degree, prime
        )
        if checked_sites != total_sites:
            raise AssertionError("site-count mismatch")
        reconstructed, current_product = _accumulate_crt(
            reconstructed, current_product, residues, prime
        )
    if reconstructed is None:
        raise AssertionError("CRT requires at least one modulus")
    if current_product != moduli.product:
        raise AssertionError("CRT modulus-product mismatch")
    if max_degree is None or int(max_degree) >= full_degree:
        expected = 1 << total_sites
        if sum(reconstructed) != expected:
            raise AssertionError(
                f"coefficient sum {sum(reconstructed)} != 2^{total_sites}"
            )
    return _trim(reconstructed)


def _torus_residues(
    shape,
    block: int,
    max_degree: int | None,
    modulus: int,
) -> tuple[np.ndarray, int, int]:
    box, cross, layers, _, cross_sites, total_sites = _shape_data(
        shape, periodic=(True,) * (len(tuple(shape)) - 1)
    )
    in_layer_bonds = layer_bonds(cross, (True,) * len(cross))
    full_degree = len(in_layer_bonds) * layers + cross_sites * layers
    degree = _degree_limit(full_degree, max_degree)
    broken_layer = _broken_in_layer(cross_sites, in_layer_bonds)
    state_count = 1 << cross_sites
    if block <= 0:
        block = state_count
    trace = np.zeros(degree + 1, dtype=np.int64)
    for start in range(0, state_count, block):
        stop = min(start + block, state_count)
        vector = np.zeros(
            (state_count, stop - start, degree + 1), dtype=np.int64
        )
        for state in range(start, stop):
            vector[state, state - start, 0] = 1
        for _ in range(layers):
            vector = _shift_rows(vector, broken_layer)
            vector = _apply_interlayer_mod(vector, cross_sites, modulus)
        for state in range(start, stop):
            trace += vector[state, state - start]
        np.remainder(trace, modulus, out=trace)
    return trace, total_sites, full_degree


def torus_broken_bond_poly(
    shape,
    block: int = 0,
    *,
    max_degree: int | None = None,
) -> list[int]:
    """Return the exact fully-periodic broken-bond polynomial by modular CRT."""

    total_sites = prod(tuple(int(side) for side in shape))
    moduli = primes_for_sites(total_sites)
    reconstructed: list[int] | None = None
    current_product = 1
    full_degree = 0
    for prime in moduli:
        residues, checked_sites, full_degree = _torus_residues(
            shape, int(block), max_degree, prime
        )
        if checked_sites != total_sites:
            raise AssertionError("site-count mismatch")
        reconstructed, current_product = _accumulate_crt(
            reconstructed, current_product, residues, prime
        )
    if reconstructed is None:
        raise AssertionError("CRT requires at least one modulus")
    if max_degree is None or int(max_degree) >= full_degree:
        expected = 1 << total_sites
        if sum(reconstructed) != expected:
            raise AssertionError(
                f"coefficient sum {sum(reconstructed)} != 2^{total_sites}"
            )
    return _trim(reconstructed)


def hybrid_box_broken_bond_poly(
    shape,
    periodic=None,
    plus_boundary: bool = False,
) -> list[int]:
    """Use the original int64 engine when safe and CRT only beyond 62 sites."""

    total_sites = prod(tuple(int(side) for side in shape))
    if total_sites <= 62:
        return _int64_box_broken_bond_poly(
            shape, periodic=periodic, plus_boundary=plus_boundary
        )
    return box_broken_bond_poly(
        shape, periodic=periodic, plus_boundary=plus_boundary
    )


def hybrid_torus_broken_bond_poly(shape, block: int = 0) -> list[int]:
    """Use the original torus engine when safe and CRT only beyond 62 sites."""

    total_sites = prod(tuple(int(side) for side in shape))
    if total_sites <= 62:
        return _int64_torus_broken_bond_poly(shape, block=block)
    return torus_broken_bond_poly(shape, block=block)


def estimate_box_peak_bytes(
    shape,
    periodic=None,
    plus_boundary: bool = False,
    *,
    max_degree: int | None = None,
) -> int:
    """Conservative peak for three simultaneous int64 propagation arrays."""

    box, cross, layers, per, cross_sites, _ = _shape_data(
        shape, periodic=periodic
    )
    bonds = layer_bonds(cross, per)
    physical = len(bonds) * layers + cross_sites * (layers - 1)
    if plus_boundary:
        target = 2 * (len(box) - 1)
        degree = [0] * cross_sites
        for left, right in bonds:
            degree[left] += 1
            degree[right] += 1
        physical += sum(target - value for value in degree) * layers + 2 * cross_sites
    retained = _degree_limit(physical, max_degree)
    return 3 * (1 << cross_sites) * (retained + 1) * np.dtype(np.int64).itemsize
