"""Standalone controls for the exact interlayer-coupling expansion."""

from __future__ import annotations

from fractions import Fraction

import mpmath as mp

from ising.interlayer import (
    anisotropic_flm_c2_series,
    anisotropic_box_even_subgraph,
    decoupled_chain_phi,
    exact_2d_overlap_series,
    log_cosh_w_series,
    torus_overlap_sum_mp,
)
from ising.lattices import cubic, square


DPS = 60


def _direct_torus_overlap(side: int, coupling: mp.mpf) -> mp.mpf:
    """Independent full-spin enumeration of sum_r <s_0 s_r>^2."""

    lattice = square(side, side, periodic=True)
    weights = []
    spin_rows = []
    for state in range(1 << lattice.n_sites):
        spins = [1 if (state >> site) & 1 else -1 for site in range(lattice.n_sites)]
        energy = sum(spins[left] * spins[right] for left, right in lattice.bonds)
        weights.append(mp.exp(coupling * energy))
        spin_rows.append(spins)
    partition = mp.fsum(weights)
    correlations = []
    for target in range(lattice.n_sites):
        numerator = mp.fsum(
            weight * spins[0] * spins[target]
            for weight, spins in zip(weights, spin_rows, strict=True)
        )
        correlations.append(numerator / partition)
    return mp.fsum(value * value for value in correlations)


def _finite_k0_chain_control() -> None:
    """Check a finite 3D torus against independent periodic 1D chains."""

    lx, ly, lz = 2, 2, 3
    lattice = cubic(lx, ly, lz, periodic=True)
    kz = mp.mpf("0.271")
    partition = mp.mpf(0)
    for state in range(1 << lattice.n_sites):
        spins = [1 if (state >> site) & 1 else -1 for site in range(lattice.n_sites)]
        vertical_energy = sum(
            spins[left] * spins[right] for left, right in lattice.bonds_by_dir[2]
        )
        partition += mp.exp(kz * vertical_energy)
    lambda_plus = 2 * mp.cosh(kz)
    lambda_minus = 2 * mp.sinh(kz)
    one_chain = lambda_plus**lz + lambda_minus**lz
    expected = one_chain ** (lx * ly)
    assert abs(partition - expected) / expected < mp.mpf("1e-55")


def main() -> None:
    overlap = exact_2d_overlap_series(6)
    assert overlap == (
        Fraction(1),
        Fraction(0),
        Fraction(4),
        Fraction(0),
        Fraction(36),
        Fraction(0),
        Fraction(236),
    )

    flm = anisotropic_flm_c2_series(6)
    expected_residual = tuple(value / 2 for value in overlap)
    expected_residual = (
        expected_residual[0] - Fraction(1, 2),
        *expected_residual[1:],
    )
    assert flm.residual == expected_residual
    assert flm.residual == anisotropic_flm_c2_series(6, bound_slack=1).residual
    assert flm.total == tuple(value / 2 for value in overlap)
    for shape in flm.boxes:
        assert all(
            row[1] == 0 for row in anisotropic_box_even_subgraph(shape, 6, 2)
        )

    prefactor = log_cosh_w_series(12)
    for degree, coefficient in enumerate(prefactor):
        expected = Fraction(1, degree) if degree and degree % 2 == 0 else Fraction(0)
        assert coefficient == expected
    # K=0 leaves no finite connected even graph in the z-chain, so every
    # computed residual coefficient vanishes and only the exact prefactor remains.
    assert flm.residual[0] == 0

    with mp.workdps(DPS):
        for kz in (mp.mpf("0"), mp.mpf("0.19"), mp.mpf("0.71")):
            exact_chain = mp.log(mp.e**kz + mp.e ** (-kz))
            assert abs(decoupled_chain_phi(kz, DPS) - exact_chain) < mp.mpf("1e-55")
        _finite_k0_chain_control()

        coupling = mp.mpf("0.23")
        transfer = torus_overlap_sum_mp(3, coupling, dps=DPS)
        direct = _direct_torus_overlap(3, coupling)
        assert abs(transfer["overlap_sum"] - direct) < mp.mpf("1e-50")
        assert transfer["identity_residual"] < mp.mpf("1e-55")

    print("PASS test_interlayer")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL test_interlayer: {error}")
        raise
