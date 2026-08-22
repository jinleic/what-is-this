"""Standalone acceptance checks for the integrability-obstruction experiments."""

from __future__ import annotations

from ising.tensor_networks import (
    exact_transfer_commutant,
    ising_tetrahedron_analysis,
    scan_local_conserved_charges,
)


PRIMES = (2_147_483_647, 2_147_483_629)


def main() -> None:
    chain = scan_local_conserved_charges(
        dimension=1,
        max_support=4,
        linear_size=9,
        field_numerator=2,
        primes=PRIMES,
    )
    chain_rows = chain["rows"]
    assert [row["max_support"] for row in chain_rows] == [2, 3, 4]
    assert all(row["kernel_dimension_by_prime"][0] == row["kernel_dimension_by_prime"][1] for row in chain_rows)
    assert [row["kernel_dimension"] for row in chain_rows] == [3, 5, 7]
    assert [row["nontrivial_dimension"] for row in chain_rows] == [1, 3, 5]
    print(
        "1D free-fermion local-charge tower: PASS",
        [(row["max_support"], row["kernel_dimension"], row["nontrivial_dimension"]) for row in chain_rows],
    )

    grid = scan_local_conserved_charges(
        dimension=2,
        max_support=4,
        linear_size=9,
        field_numerator=2,
        primes=PRIMES,
    )
    grid_rows = grid["rows"]
    assert [row["max_support"] for row in grid_rows] == [2, 3, 4]
    assert all(row["kernel_dimension_by_prime"][0] == row["kernel_dimension_by_prime"][1] for row in grid_rows)
    assert all(row["kernel_dimension"] == row["trivial_dimension"] == 2 for row in grid_rows)
    print(
        "2D layer local-charge search: PASS",
        [(row["max_support"], row["ansatz_dimension"], row["kernel_dimension"]) for row in grid_rows],
    )

    commutants = []
    for shape in ((4,), (2, 2)):
        for x in ((1, 2), (1, 3)):
            row = exact_transfer_commutant(shape, x_numerator=x[0], x_denominator=x[1])
            assert row["method"] == "exact_characteristic_polynomial"
            assert row["diagonalizable"] is True
            assert row["commutant_dimension"] >= row["matrix_dimension"]
            commutants.append((shape, x, row["commutant_dimension"]))
    assert [row[2] for row in commutants] == [16, 16, 32, 32]
    print("exact transfer-matrix commutants: PASS", commutants)

    tetrahedron = ising_tetrahedron_analysis()
    assert tetrahedron["general_system"]["variables"] == 64
    assert tetrahedron["general_system"]["equations"] == 4096
    assert tetrahedron["isotropic_substitution"]["nonzero_equations"] == 832
    assert tetrahedron["isotropic_substitution"]["common_gcd"] == "w**2*(w - 1)**3*(w + 1)**3*(w**2 + 1)**2"
    assert len(tetrahedron["isotropic_substitution"]["reduced_groebner_basis"]) == 1
    assert tetrahedron["isotropic_substitution"]["physical_interval_solutions"] == []
    assert tetrahedron["arbitrary_leg_gl2_gauge"]["satisfies_tetrahedron_equation"] is True
    assert tetrahedron["delta_preserving_identical_tensor_gauge"]["possible_for_finite_positive_K"] is False
    print("tetrahedron-equation specialization and gauge scope: PASS")
    print("PASS")


if __name__ == "__main__":
    main()
