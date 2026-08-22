"""Standalone acceptance checks for the Onsager--Kaufman control experiment."""

from __future__ import annotations

import mpmath as mp

from ising.onsager import (
    anisotropic_torus_broken_bond_counts,
    collapse_bivariate_counts,
    evaluate_bivariate_partition,
    free_fermion_diagnostics,
    kaufman_torus_Z,
    onsager_free_energy,
    reconstruct_torus_Z_from_generators,
)
from ising.transfer_matrix import torus_broken_bond_poly


PRECISION = 60
EXACT_CASES = (
    (3, 3, "0.17", "0.17"),
    (4, 4, "0.29", "0.29"),
    (5, 5, "0.4406867935097715126163046624898961545", "0.4406867935097715126163046624898961545"),
    (6, 6, "0.61", "0.61"),
    (3, 4, "0.19", "0.37"),
    (4, 5, "0.31", "0.47"),
    (5, 6, "0.58", "0.21"),
    (6, 3, "0.26", "0.53"),
)


def _relative_error(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    return abs(a - b) / abs(b)


def main() -> None:
    mp.mp.dps = PRECISION

    # (a) Kaufman's four-sector formula against exact integer transfer matrices.
    exact_rows = []
    for m, n, sx, sy in EXACT_CASES:
        kx, ky = mp.mpf(sx), mp.mpf(sy)
        counts = anisotropic_torus_broken_bond_counts(m, n)
        collapsed = collapse_bivariate_counts(counts)
        assert collapsed == torus_broken_bond_poly((n, m)), (m, n)
        exact = evaluate_bivariate_partition(counts, kx, ky)
        kaufman = kaufman_torus_Z(m, n, kx, ky)
        rel = _relative_error(kaufman, exact)
        assert rel < mp.mpf("1e-40"), (m, n, mp.nstr(rel, 12))
        exact_rows.append((m, n, rel))
    print("Kaufman vs exact integer TM: PASS", len(exact_rows), "cases; max rel =", mp.nstr(max(r[2] for r in exact_rows), 8))

    # (b) The finite torus approaches the normalized double-integral free energy.
    kc = mp.asinh(1) / 2
    phi = onsager_free_energy(kc, kc)
    # At criticality the integral also has this independent closed evaluation.
    phi_critical = mp.log(2) / 2 + 2 * mp.catalan / mp.pi
    assert abs(phi - phi_critical) < mp.mpf("1e-50")
    convergence = []
    for length in (64, 128, 256):
        finite_phi = mp.log(kaufman_torus_Z(length, length, kc, kc)) / (length * length)
        residual = finite_phi - phi
        assert residual > 0
        convergence.append((length, residual))
        print(f"L={length:3d} free-energy residual = {mp.nstr(residual, 18)}")
    assert convergence[1][1] < convergence[0][1] / mp.mpf("3.9")
    assert convergence[2][1] < convergence[1][1] / mp.mpf("3.9")

    # (c) The zero mode and logarithmically growing curvature locate the singular line.
    assert abs(mp.sinh(2 * kc) ** 2 - 1) < mp.mpf("1e-55")
    hs = (mp.mpf("0.002"), mp.mpf("0.001"), mp.mpf("0.0005"))
    curvatures = []
    for h in hs:
        curvature = (
            onsager_free_energy(kc + h, kc + h)
            - 2 * phi
            + onsager_free_energy(kc - h, kc - h)
        ) / h**2
        curvatures.append(curvature)
    assert curvatures[0] < curvatures[1] < curvatures[2]
    print("critical sinh(2Kx)sinh(2Ky)=1 and growing curvature: PASS", [mp.nstr(v, 10) for v in curvatures])

    # (d) log(V) has only scalar and Majorana-bilinear Pauli support after parity resolution.
    off_support = []
    reconstruction_errors = []
    energy_errors = []
    for n in (3, 4, 5, 6):
        for parity in (+1, -1):
            diagnostic = free_fermion_diagnostics(n, 0.31, 0.47, parity_sector=parity)
            off_support.append(diagnostic["max_non_bilinear_coefficient"])
            energy_errors.append(diagnostic["max_energy_residual"])
            assert diagnostic["max_non_bilinear_coefficient"] < 1e-11, (n, parity, diagnostic)
            assert diagnostic["max_energy_residual"] < 1e-11, (n, parity, diagnostic)
            assert diagnostic["max_generator_identity_residual"] < 1e-13, (n, parity, diagnostic)
            assert diagnostic["boundary_parity_identity_residual"] < 1e-13, (n, parity, diagnostic)
        reconstructed = reconstruct_torus_Z_from_generators(n, n, 0.31, 0.47)
        counts = anisotropic_torus_broken_bond_counts(n, n)
        exact = evaluate_bivariate_partition(counts, mp.mpf("0.31"), mp.mpf("0.47"))
        rel = abs(mp.mpf(reconstructed["Z"]) - exact) / exact
        reconstruction_errors.append(float(rel))
        assert rel < mp.mpf("1e-10"), (n, rel)
    print("Majorana-bilinear support: PASS; max off =", f"{max(off_support):.3e}")
    print("generator energies and torus reconstruction: PASS; max rel =", f"{max(reconstruction_errors):.3e}")
    print("PASS")


if __name__ == "__main__":
    main()
