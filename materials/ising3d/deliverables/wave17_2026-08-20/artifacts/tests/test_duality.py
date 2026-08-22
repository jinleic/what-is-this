"""Standalone exact-integer checks for finite-volume Ising--gauge duality."""

from __future__ import annotations

from ising.duality import (
    closed_surface_polynomial,
    cubical_complex,
    dual_ising_sector_polynomials,
    gauge_satisfied_plaquette_dos,
    homology_dimensions,
    ising_torus_sector_polynomials_2d,
    kw_character,
    signed_high_temperature_polynomial,
    surface_polynomial_from_gauge_dos,
    surface_sector_representatives,
)
from ising.exact_enumeration import dos_bonds
from ising.lattices import cubic
from ising.transfer_matrix import torus_broken_bond_poly


CASES_3D = (
    ((2, 2, 2), (False, False, False)),
    ((3, 2, 2), (False, False, False)),
    ((2, 2, 2), (True, False, False)),
    ((2, 2, 2), (True, True, False)),
    ((2, 2, 2), (True, True, True)),
)

KNOWN_4X4_EVEN_SUBGRAPH = [
    1,
    0,
    0,
    0,
    24,
    0,
    128,
    0,
    876,
    0,
    3584,
    0,
    13160,
    0,
    28032,
    0,
    39462,
    0,
    28032,
    0,
    13160,
    0,
    3584,
    0,
    876,
    0,
    128,
    0,
    24,
    0,
    0,
    0,
    1,
]


def _sum_polynomials(polynomials: tuple[tuple[int, ...], ...]) -> list[int]:
    degree = max(len(poly) for poly in polynomials)
    return [sum(poly[k] if k < len(poly) else 0 for poly in polynomials) for k in range(degree)]


def _check_3d_case(shape: tuple[int, int, int], periodic: tuple[bool, bool, bool]) -> None:
    complex_ = cubical_complex(shape, periodic)
    homology = homology_dimensions(complex_)

    gauge_dos = gauge_satisfied_plaquette_dos(complex_)
    assert sum(gauge_dos) == 1 << complex_.n_links

    from_gauge = surface_polynomial_from_gauge_dos(complex_, gauge_dos)
    direct = closed_surface_polynomial(complex_)
    assert from_gauge == direct
    assert direct[0] == 1
    assert sum(direct) == 1 << (complex_.n_plaquettes - homology["rank_d2"])

    representatives = surface_sector_representatives(complex_)
    assert len(representatives) == 1 << homology["b2"]
    sector_polynomials = dual_ising_sector_polynomials(complex_, representatives)
    assert all(sum(poly) == 1 << complex_.n_cubes for poly in sector_polynomials)
    dual_sum = _sum_polynomials(sector_polynomials)
    kernel_factor = 1 << homology["b3"]
    assert dual_sum == [kernel_factor * coefficient for coefficient in direct]
    if periodic == (True, True, True) and shape == (2, 2, 2):
        independent_ising_dos = dos_bonds(cubic(*shape, periodic=True))
        independent_broken_bond_polynomial = tuple(
            int(value) for value in reversed(independent_ising_dos)
        )
        assert sector_polynomials[0] == independent_broken_bond_polynomial

    expected_flat = 1 << (complex_.n_links - homology["rank_d2"])
    assert gauge_dos[complex_.n_plaquettes] == expected_flat
    assert expected_flat == (1 << homology["rank_d1"]) * (1 << homology["b1"])
    print(f"3D {complex_.describe()}: PASS")


def _check_2d_control(shape: tuple[int, int]) -> None:
    sector_polynomials = ising_torus_sector_polynomials_2d(shape)
    n_sites = shape[0] * shape[1]
    n_edges = 2 * n_sites
    assert len(sector_polynomials) == 4
    assert all(sum(poly) == 1 << n_sites for poly in sector_polynomials)
    independent_polynomial = torus_broken_bond_poly(shape)
    assert list(sector_polynomials[0][: len(independent_polynomial)]) == independent_polynomial
    assert not any(sector_polynomials[0][len(independent_polynomial) :])

    high_temperature = tuple(
        signed_high_temperature_polynomial(poly, n_sites) for poly in sector_polynomials
    )
    for alpha in range(4):
        rhs = [
            sum(kw_character(alpha, beta) * sector_polynomials[beta][degree] for beta in range(4))
            for degree in range(n_edges + 1)
        ]
        lhs = [2 * coefficient for coefficient in high_temperature[alpha]]
        assert lhs == rhs, (shape, alpha)

    if shape == (4, 4):
        assert list(high_temperature[0]) == KNOWN_4X4_EVEN_SUBGRAPH
    assert high_temperature[0][0] == 1
    assert sum(high_temperature[0]) == 1 << (n_edges - n_sites + 1)
    print(f"2D {shape[0]}x{shape[1]} torus, four sectors: PASS")


def main() -> None:
    open_box = cubical_complex((2, 2, 2), False)
    assert (open_box.n_vertices, open_box.n_links, open_box.n_plaquettes, open_box.n_cubes) == (
        8,
        12,
        6,
        1,
    )
    assert tuple(homology_dimensions(open_box)[f"b{k}"] for k in range(4)) == (1, 0, 0, 0)

    torus = cubical_complex((2, 2, 2), True)
    assert (torus.n_vertices, torus.n_links, torus.n_plaquettes, torus.n_cubes) == (8, 24, 24, 8)
    assert tuple(homology_dimensions(torus)[f"b{k}"] for k in range(4)) == (1, 3, 3, 1)
    print("cell counts and homology: PASS")

    for shape, periodic in CASES_3D:
        _check_3d_case(shape, periodic)

    _check_2d_control((4, 4))
    _check_2d_control((5, 5))
    print("PASS")


if __name__ == "__main__":
    main()
