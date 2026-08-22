"""Standalone acceptance checks for the exact Lee--Yang machinery."""

from __future__ import annotations

import mpmath as mp
import numpy as np

from ising.exact_enumeration import dos_bonds
from ising.lattices import cubic, square
from ising.lee_yang import (
    circle_residual,
    field_dos_from_enumeration,
    field_polynomial_at_K,
    lee_yang_roots,
    normalized_root_residual,
    rational_field_polynomial,
    rational_field_polynomial_transfer,
    transfer_field_dos,
)
from ising.transfer_matrix import box_broken_bond_poly, torus_broken_bond_poly


DPS = 100
CIRCLE_TOLERANCE = mp.mpf("1e-40")
ROOT_TOLERANCE = mp.mpf("1e-70")


def _trim(values: list[int]) -> list[int]:
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def _check_enumerated_lattice(label, lat) -> tuple[int, str]:
    field_dos = field_dos_from_enumeration(lat)
    n_sites = lat.n_sites
    n_bonds = lat.n_bonds

    # Spin reversal leaves the broken-bond count q fixed and sends
    # k=#down spins to N-k, hence c[k,q] = c[N-k,q] exactly.
    assert np.array_equal(field_dos, field_dos[::-1]), label

    # At z=1 the field grading disappears.  dos_bonds is indexed by the
    # satisfied-bond count b=N_b-q, whereas field_dos is indexed by q.
    expected = [int(value) for value in dos_bonds(lat)][::-1]
    actual = [int(value) for value in field_dos.sum(axis=0)]
    assert actual == expected, label

    coefficients = rational_field_polynomial(field_dos, 2, 3)
    assert coefficients == coefficients[::-1], label
    roots = lee_yang_roots(coefficients, dps=DPS)
    assert len(roots) == n_sites, label
    radial_error = circle_residual(roots)
    root_error = normalized_root_residual(coefficients, roots)
    assert radial_error < CIRCLE_TOLERANCE, (label, mp.nstr(radial_error, 8))
    assert root_error < ROOT_TOLERANCE, (label, mp.nstr(root_error, 8))
    if n_sites % 2:
        numerical = field_polynomial_at_K(field_dos, "0.221654626", dps=DPS)
        numerical_roots = lee_yang_roots(numerical, dps=DPS)
        assert len(numerical_roots) == n_sites, label
    return n_sites, f"{label}: | |z|-1 |={mp.nstr(radial_error, 5)}"


def _check_transfer_square(side: int, periodic: bool) -> tuple[int, str]:
    label = f"square_{side}x{side}_{'periodic' if periodic else 'open'}"
    field_dos = transfer_field_dos((side, side), periodic=periodic)
    n_sites = side * side
    assert np.array_equal(field_dos, field_dos[::-1]), label

    expected = (
        torus_broken_bond_poly((side, side))
        if periodic
        else box_broken_bond_poly((side, side))
    )
    actual = _trim([int(value) for value in field_dos.sum(axis=0)])
    assert actual == expected, label

    coefficients = rational_field_polynomial(field_dos, 2, 3)
    roots = lee_yang_roots(coefficients, dps=DPS)
    assert len(roots) == n_sites, label
    radial_error = circle_residual(roots)
    root_error = normalized_root_residual(coefficients, roots)
    assert radial_error < CIRCLE_TOLERANCE, (label, mp.nstr(radial_error, 8))
    assert root_error < ROOT_TOLERANCE, (label, mp.nstr(root_error, 8))
    return n_sites, f"{label}: | |z|-1 |={mp.nstr(radial_error, 5)}"


def test_lee_yang_acceptance() -> None:
    checked: list[str] = []
    largest = 0

    lattices = []
    for shape in ((2, 2, 2), (3, 3, 2), (3, 3, 3)):
        for periodic in (False, True):
            label = f"cubic_{'x'.join(map(str, shape))}_{'periodic' if periodic else 'open'}"
            lattices.append((label, cubic(*shape, periodic=periodic)))
    for side in (4, 5):
        for periodic in (False, True):
            label = f"square_{side}x{side}_{'periodic' if periodic else 'open'}"
            lattices.append((label, square(side, side, periodic=periodic)))

    for label, lat in lattices:
        n_sites, detail = _check_enumerated_lattice(label, lat)
        largest = max(largest, n_sites)
        checked.append(detail)

    # Validate both new transfer axes against full enumeration before using
    # the transfer engine beyond joint_dos's size limit.
    for periodic in (False, True):
        enumerated = field_dos_from_enumeration(square(4, 4, periodic=periodic))
        transferred = transfer_field_dos((4, 4), periodic=periodic)
        assert np.array_equal(transferred, enumerated)

    for periodic in (False, True):
        n_sites, detail = _check_transfer_square(6, periodic)
        largest = max(largest, n_sites)
        checked.append(detail)

    enumerated_3d = field_dos_from_enumeration(cubic(3, 3, 2, periodic=False))
    expected_rational = rational_field_polynomial(enumerated_3d, 2, 3)
    transferred_rational = rational_field_polynomial_transfer(
        (3, 3, 2), numerator=2, denominator=3, periodic=False
    )
    assert transferred_rational == expected_rational

    # The 64-site open cube is beyond joint_dos's N<=30 limit.  This is an
    # exact integer polynomial at x=2/3: a common factor 3^N_b clears all
    # bond-weight denominators and does not change the roots.
    coefficients_64 = rational_field_polynomial_transfer(
        (4, 4, 4), numerator=2, denominator=3, periodic=False
    )
    assert coefficients_64 == coefficients_64[::-1]
    roots_64 = lee_yang_roots(coefficients_64, dps=DPS)
    assert len(roots_64) == 64
    radial_error_64 = circle_residual(roots_64)
    root_error_64 = normalized_root_residual(coefficients_64, roots_64)
    assert radial_error_64 < CIRCLE_TOLERANCE, mp.nstr(radial_error_64, 8)
    assert root_error_64 < ROOT_TOLERANCE, mp.nstr(root_error_64, 8)
    largest = max(largest, 64)
    checked.append(f"cubic_4x4x4_open: | |z|-1 |={mp.nstr(radial_error_64, 5)}")

    assert largest == 64
    print("PASS: Lee-Yang circle, spin-reversal palindrome, Z(1), and root count")
    print(f"PASS: {len(checked)} lattices; largest lattice reached = 4x4x4 ({largest} sites)")
    for detail in checked:
        print("  ", detail)


if __name__ == "__main__":
    test_lee_yang_acceptance()
