"""Exact validation gates for the finite-lattice series engine.

Run standalone from the repository root:
    .venv/bin/python tests/test_series.py
"""

from fractions import Fraction
from functools import lru_cache
from math import comb

from ising.lattices import cubic

from ising.series import (
    broken_bond_to_even_subgraph,
    high_temperature_free_energy,
    log_series,
    low_temperature_free_energy,
    onsager_square_ht_series,
)
from ising.transfer_matrix import box_broken_bond_poly


@lru_cache(maxsize=None)
def _square_flm():
    return high_temperature_free_energy(2, 24)


@lru_cache(maxsize=None)
def _sc_ht_validation(slack=0):
    return high_temperature_free_energy(3, 16, bound_slack=slack)


@lru_cache(maxsize=None)
def _sc_ht_result():
    return _sc_ht_validation(0)


@lru_cache(maxsize=None)
def _sc_lt(slack=0):
    return low_temperature_free_energy(3, 20, bound_slack=slack)


def _assert_exact(coefficients):
    assert all(isinstance(value, Fraction) for value in coefficients)


def test_exact_log_composition():
    # log(1 + 2t + t^2) = 2 log(1+t), including nonintegral coefficients.
    got = log_series([1, 2, 1], 6)
    expected = tuple([Fraction(0)] + [Fraction(2 * (-1) ** (n + 1), n) for n in range(1, 7)])
    assert got == expected
    _assert_exact(got)


def test_broken_bond_conversion_is_integral():
    # The free 2x2 square has exactly the empty even subgraph and its four-edge plaquette.
    broken = box_broken_bond_poly((2, 2))
    got = broken_bond_to_even_subgraph(broken, (2, 2))
    assert got == (1, 0, 0, 0, 1)
    assert all(isinstance(value, int) for value in got)


def test_square_flm_matches_twelve_onsager_orders():
    flm = _square_flm()
    onsager = onsager_square_ht_series(24)
    assert flm.coefficients == onsager
    nonzero_orders = [n for n in range(1, 25) if onsager[n] != 0]
    assert nonzero_orders == list(range(2, 25, 2))
    assert len(nonzero_orders) == 12
    _assert_exact(flm.coefficients)
    _assert_exact(flm.interaction_coefficients)


def test_sc_ht_bound_increase_is_inert():
    minimal = _sc_ht_validation(0)
    enlarged = _sc_ht_validation(1)
    assert minimal.coefficients == enlarged.coefficients
    assert minimal.interaction_coefficients == enlarged.interaction_coefficients
    extra = set(enlarged.boxes) - set(minimal.boxes)
    assert extra
    for shape in extra:
        assert all(value == 0 for value in enlarged.box_weights[shape])
    _assert_exact(minimal.coefficients)
    for weight in enlarged.box_weights.values():
        _assert_exact(weight)


def _count_four_cycles(lattice):
    adjacency = [set() for _ in range(lattice.n_sites)]
    for left, right in lattice.bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    opposite_pair_count = sum(
        comb(len(adjacency[left] & adjacency[right]), 2)
        for left in range(lattice.n_sites)
        for right in range(left + 1, lattice.n_sites)
    )
    assert opposite_pair_count % 2 == 0
    return opposite_pair_count // 2


def test_sc_ht_plaquette_coefficient():
    series = _sc_ht_result()
    assert series is _sc_ht_validation(0)
    # A 5^3 torus has no length-four winding loops.  Direct graph enumeration
    # counts each four-cycle through its two pairs of opposite vertices.
    lattice = cubic(5, 5, 5, periodic=True)
    cycles = _count_four_cycles(lattice)
    assert cycles == 3 * lattice.n_sites
    # At order v^4 log P has no products of lower graphs, so every square has
    # cluster weight +1.
    assert series.interaction_coefficients[4] == Fraction(cycles, lattice.n_sites)
    # 3 log(cosh K) = -(3/2) log(1-v^2) contributes 3/4 at v^4.
    assert series.coefficients[4] == Fraction(15, 4)
    assert series.order == 16
    _assert_exact(series.coefficients)


def test_sc_lt_bound_increase_and_single_spin_term():
    minimal = _sc_lt(0)
    enlarged = _sc_lt(1)
    assert minimal.coefficients == enlarged.coefficients
    extra = set(enlarged.boxes) - set(minimal.boxes)
    assert extra
    for shape in extra:
        assert all(value == 0 for value in enlarged.box_weights[shape])
    # A single flipped spin has exactly six broken bonds and one translation per site.
    assert minimal.coefficients[6] == 1
    assert minimal.coefficients[10] == 3
    assert minimal.coefficients[12] == Fraction(-7, 2)
    _assert_exact(minimal.coefficients)


def main():
    tests = [
        test_exact_log_composition,
        test_broken_bond_conversion_is_integral,
        test_square_flm_matches_twelve_onsager_orders,
        test_sc_ht_bound_increase_is_inert,
        test_sc_ht_plaquette_coefficient,
        test_sc_lt_bound_increase_and_single_spin_term,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(
        "PASS: 2D HT through v^24 (12 nonzero orders); "
        "3D HT through v^16 (enlarged-box check through v^16); "
        "3D LT through x^20"
    )


if __name__ == "__main__":
    main()
