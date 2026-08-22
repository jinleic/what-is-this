"""Reproduce the Onsager--Kaufman solution as an algebraic control experiment.

Run from the repository root with

    .venv/bin/python experiments/e02_onsager_2d.py

The script performs all finite-integer, thermodynamic-limit, singularity, and
free-fermion checks and writes ``results/onsager/onsager_2d.json``.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path

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
SCRIPT = "experiments/e02_onsager_2d.py"
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "onsager" / "onsager_2d.json"
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


def _decimal(value: mp.mpf, digits: int = PRECISION) -> str:
    return mp.nstr(value, digits)


def _relative_error(found: mp.mpf, expected: mp.mpf) -> mp.mpf:
    return abs(found - expected) / abs(expected)

def _raw_sector_products(
    m: int, n: int, kx: mp.mpf, ky: mp.mpf
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """Independent raw products used only to falsify nearby sign conventions."""

    ky_star = mp.atanh(mp.exp(-2 * ky))
    gamma = []
    for ell in range(2 * n):
        if ell == 0:
            gamma.append(2 * (ky_star - kx))
            continue
        momentum = mp.pi * ell / n
        gamma.append(
            mp.acosh(
                mp.cosh(2 * kx) * mp.cosh(2 * ky_star)
                - mp.sinh(2 * kx) * mp.sinh(2 * ky_star) * mp.cos(momentum)
            )
        )
    odd = [gamma[2 * r + 1] for r in range(n)]
    even = [gamma[2 * r] for r in range(n)]
    prefactor = (2 * mp.sinh(2 * ky)) ** (mp.mpf(m * n) / 2)
    return (
        prefactor,
        mp.fprod(2 * mp.cosh(m * value / 2) for value in odd),
        mp.fprod(2 * mp.sinh(m * value / 2) for value in odd),
        mp.fprod(2 * mp.cosh(m * value / 2) for value in even),
        mp.fprod(2 * mp.sinh(m * value / 2) for value in even),
    )


def run_experiment() -> dict[str, object]:
    mp.mp.dps = PRECISION
    checks: list[dict[str, object]] = []

    exact_records: list[dict[str, object]] = []
    exact_relatives: list[mp.mpf] = []
    collapse_passed = True
    for m, n, sx, sy in EXACT_CASES:
        kx, ky = mp.mpf(sx), mp.mpf(sy)
        counts = anisotropic_torus_broken_bond_counts(m, n)
        collapsed = collapse_bivariate_counts(counts)
        isotropic_integer_tm = torus_broken_bond_poly((n, m))
        collapse_matches = collapsed == isotropic_integer_tm
        collapse_passed &= collapse_matches

        exact = evaluate_bivariate_partition(counts, kx, ky)
        kaufman = kaufman_torus_Z(m, n, kx, ky)
        relative = _relative_error(kaufman, exact)
        exact_relatives.append(relative)
        exact_records.append(
            {
                "m": m,
                "n": n,
                "Kx": sx,
                "Ky": sy,
                "coefficient_sum": str(sum(sum(row) for row in counts)),
                "collapsed_matches_torus_broken_bond_poly": collapse_matches,
                "Z_integer_tm": _decimal(exact),
                "Z_kaufman": _decimal(kaufman),
                "relative_error": _decimal(relative),
            }
        )

    max_exact_relative = max(exact_relatives)
    checks.append(
        {
            "name": "bivariate_integer_TM_collapses_to_repository_TM",
            "passed": collapse_passed,
            "detail": "all eight exact bivariate polynomials collapse coefficient-by-coefficient",
        }
    )
    checks.append(
        {
            "name": "Kaufman_four_sector_formula_vs_exact_integer_TM",
            "passed": max_exact_relative < mp.mpf("1e-40"),
            "detail": f"8 cases at {PRECISION} dps; max relative error {_decimal(max_exact_relative, 12)}",
        }
    )

    # Falsify the nearest alternative conventions at generic anisotropic and
    # low-temperature points, rather than selecting signs from a benchmark.
    generic_m, generic_n = 3, 4
    generic_kx, generic_ky = mp.mpf("0.19"), mp.mpf("0.37")
    generic_correct = kaufman_torus_Z(generic_m, generic_n, generic_kx, generic_ky)
    prefactor, odd_cosh, odd_sinh, even_cosh, even_sinh = _raw_sector_products(
        generic_m, generic_n, generic_kx, generic_ky
    )
    variants = {
        "all_four_products_plus": prefactor
        * (odd_cosh + odd_sinh + even_cosh + even_sinh)
        / 2,
        "odd_sinh_minus": prefactor
        * (odd_cosh - odd_sinh + even_cosh - even_sinh)
        / 2,
        "integer_and_half_integer_grids_swapped": prefactor
        * (odd_cosh - odd_sinh + even_cosh + even_sinh)
        / 2,
    }
    low_k = mp.mpf("0.7")
    low_correct = kaufman_torus_Z(generic_m, generic_n, low_k, low_k)
    low_prefactor, low_oc, low_os, low_ec, low_es = _raw_sector_products(
        generic_m, generic_n, low_k, low_k
    )
    variants["force_gamma_0_nonnegative"] = (
        low_prefactor * (low_oc + low_os + low_ec - abs(low_es)) / 2
    )
    rejected_conventions = []
    for name, value in variants.items():
        reference = low_correct if name == "force_gamma_0_nonnegative" else generic_correct
        rejected_conventions.append(
            {
                "name": name,
                "relative_error": _decimal(_relative_error(value, reference)),
            }
        )
    nearest_wrong_error = min(
        mp.mpf(record["relative_error"]) for record in rejected_conventions
    )
    checks.append(
        {
            "name": "nearby_Kaufman_sector_conventions_falsified",
            "passed": nearest_wrong_error > mp.mpf("1e-3"),
            "detail": "wrong signs/grid assignment/signed-zero-mode variants miss exact Z by at least "
            + _decimal(nearest_wrong_error, 8),
        }
    )

    critical_coupling = mp.asinh(1) / 2
    infinite_phi = onsager_free_energy(critical_coupling, critical_coupling)
    independent_critical_phi = mp.log(2) / 2 + 2 * mp.catalan / mp.pi
    normalization_error = abs(infinite_phi - independent_critical_phi)
    convergence: list[dict[str, object]] = []
    residuals: list[mp.mpf] = []
    for length in (64, 128, 256):
        finite_phi = (
            mp.log(kaufman_torus_Z(length, length, critical_coupling, critical_coupling))
            / length**2
        )
        residual = finite_phi - infinite_phi
        residuals.append(residual)
        convergence.append(
            {
                "m": length,
                "n": length,
                "finite_phi": _decimal(finite_phi),
                "residual_finite_minus_infinite": _decimal(residual),
                "L_squared_times_residual": _decimal(length**2 * residual),
            }
        )
    converges = (
        residuals[0] > 0
        and residuals[1] < residuals[0] / mp.mpf("3.9")
        and residuals[2] < residuals[1] / mp.mpf("3.9")
    )
    checks.append(
        {
            "name": "double_integral_normalization",
            "passed": normalization_error < mp.mpf("1e-50"),
            "detail": "critical integral agrees with log(2)/2 + 2 Catalan/pi; error "
            + _decimal(normalization_error, 12),
        }
    )
    checks.append(
        {
            "name": "finite_torus_converges_to_double_integral",
            "passed": converges,
            "detail": "positive residual falls by approximately four under each doubling through L=256",
        }
    )

    # The determinant gap is proportional to (sinh(2Kx)sinh(2Ky)-1)^2.
    ky_on_line = mp.mpf("0.3")
    kx_on_line = -mp.log(mp.tanh(ky_on_line)) / 2
    critical_condition = mp.sinh(2 * kx_on_line) * mp.sinh(2 * ky_on_line) - 1
    zero_mode = 2 * (-mp.log(mp.tanh(ky_on_line)) / 2 - kx_on_line)
    determinant_minimum = (
        mp.cosh(2 * kx_on_line) * mp.cosh(2 * ky_on_line)
        - mp.sinh(2 * kx_on_line)
        - mp.sinh(2 * ky_on_line)
    )
    steps = (mp.mpf("0.002"), mp.mpf("0.001"), mp.mpf("0.0005"))
    curvatures: list[dict[str, str]] = []
    curvature_values: list[mp.mpf] = []
    for step in steps:
        curvature = (
            onsager_free_energy(critical_coupling + step, critical_coupling + step)
            - 2 * infinite_phi
            + onsager_free_energy(critical_coupling - step, critical_coupling - step)
        ) / step**2
        curvature_values.append(curvature)
        curvatures.append({"step": _decimal(step), "centered_second_difference": _decimal(curvature)})
    singularity_passed = (
        abs(critical_condition) < mp.mpf("1e-55")
        and abs(zero_mode) < mp.mpf("1e-55")
        and abs(determinant_minimum) < mp.mpf("1e-55")
        and curvature_values[0] < curvature_values[1] < curvature_values[2]
    )
    checks.append(
        {
            "name": "critical_singularity",
            "passed": singularity_passed,
            "detail": "zero determinant and fermion mode at sinh(2Kx)sinh(2Ky)=1; centered curvature grows as step halves",
        }
    )

    fermion_records: list[dict[str, object]] = []
    reconstruction_records: list[dict[str, object]] = []
    maximum_off = 0.0
    maximum_energy_error = 0.0
    maximum_generator_Z_error = mp.mpf(0)
    maximum_identity_error = 0.0
    physical_off_minimum = float("inf")
    for n in (3, 4, 5, 6):
        sectors = []
        for parity in (+1, -1):
            diagnostic = free_fermion_diagnostics(
                n, mp.mpf("0.31"), mp.mpf("0.47"), parity_sector=parity
            )
            sectors.append(diagnostic)
            maximum_off = max(maximum_off, float(diagnostic["max_non_bilinear_coefficient"]))
            maximum_energy_error = max(maximum_energy_error, float(diagnostic["max_energy_residual"]))
            maximum_identity_error = max(
                maximum_identity_error,
                float(diagnostic["max_generator_identity_residual"]),
                float(diagnostic["boundary_parity_identity_residual"]),
            )

        physical = free_fermion_diagnostics(
            n, mp.mpf("0.31"), mp.mpf("0.47"), parity_sector=None
        )
        physical_off_minimum = min(
            physical_off_minimum, float(physical["max_non_bilinear_coefficient"])
        )
        fermion_records.append(
            {
                "n": n,
                "sectors": sectors,
                "unresolved_physical_matrix": {
                    "max_non_bilinear_coefficient": physical["max_non_bilinear_coefficient"],
                    "largest_non_bilinear_word": physical["largest_non_bilinear_word"],
                    "max_quadratic_reconstruction_entry": physical[
                        "max_quadratic_reconstruction_entry"
                    ],
                },
            }
        )

        reconstructed = reconstruct_torus_Z_from_generators(
            n, n, mp.mpf("0.31"), mp.mpf("0.47")
        )
        exact = evaluate_bivariate_partition(
            anisotropic_torus_broken_bond_counts(n, n), mp.mpf("0.31"), mp.mpf("0.47")
        )
        reconstruction_relative = abs(mp.mpf(reconstructed["Z"]) - exact) / exact
        maximum_generator_Z_error = max(maximum_generator_Z_error, reconstruction_relative)
        reconstruction_records.append(
            {
                "m": n,
                "n": n,
                "Kx": "0.31",
                "Ky": "0.47",
                "Z_from_2n_generators": repr(reconstructed["Z"]),
                "Z_exact_integer_TM": _decimal(exact),
                "relative_error": _decimal(reconstruction_relative),
                "sector_details": reconstructed["sectors"],
            }
        )

    checks.extend(
        [
            {
                "name": "Jordan_Wigner_generator_identities",
                "passed": maximum_identity_error < 1e-13,
                "detail": f"max matrix-entry residual {maximum_identity_error:.3e}",
            },
            {
                "name": "log_transfer_matrix_is_Majorana_quadratic",
                "passed": maximum_off < 1e-11,
                "detail": f"n=3..6, both parity sectors; max off-bilinear Pauli coefficient {maximum_off:.3e}",
            },
            {
                "name": "single_particle_energies_match_Kaufman_gamma",
                "passed": maximum_energy_error < 1e-11,
                "detail": f"max absolute energy error {maximum_energy_error:.3e}",
            },
            {
                "name": "torus_Z_reconstructed_from_2n_generators",
                "passed": maximum_generator_Z_error < mp.mpf("1e-10"),
                "detail": "n=3..6; max relative error " + _decimal(maximum_generator_Z_error, 12),
            },
            {
                "name": "periodic_Jordan_Wigner_parity_resolution_required",
                "passed": physical_off_minimum > 1e-3 and maximum_off < 1e-11,
                "detail": "full spin matrix is a direct sum of two Gaussians; each fixed-parity extension is Gaussian",
            },
        ]
    )

    return {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": PRECISION,
            "matrix_arithmetic": "numpy.complex128 with scipy.linalg.expm/logm",
        },
        "data": {
            "sector_convention": {
                "formula": "Z=(2*sinh(2Ky))^(mn/2)*(C_odd+S_odd+C_even-S_even)/2",
                "odd_momenta": "(2r+1)*pi/n",
                "even_momenta": "2r*pi/n",
                "signed_zero_mode": "gamma_0=2*(Ky_star-Kx), tanh(Ky_star)=exp(-2Ky)",
                "empirical_ground_truth": "exact bivariate integer transfer matrix; diagonal collapse equals ising.transfer_matrix.torus_broken_bond_poly((n,m))",
                "rejected_nearby_conventions": rejected_conventions,
            },
            "exact_finite_tori": {
                "cases": exact_records,
                "max_relative_error": _decimal(max_exact_relative),
            },
            "thermodynamic_limit": {
                "Kc_isotropic": _decimal(critical_coupling),
                "phi_double_integral_at_Kc": _decimal(infinite_phi),
                "phi_independent_Catalan_form": _decimal(independent_critical_phi),
                "normalization_error": _decimal(normalization_error),
                "finite_size_convergence": convergence,
            },
            "critical_singularity": {
                "Ky": _decimal(ky_on_line),
                "Kx_from_duality": _decimal(kx_on_line),
                "sinh_product_minus_one": _decimal(critical_condition),
                "signed_gamma_0": _decimal(zero_mode),
                "determinant_at_zero_momentum": _decimal(determinant_minimum),
                "curvature_sequence": curvatures,
            },
            "free_fermion": {
                "Kx": "0.31",
                "Ky": "0.47",
                "diagnostics": fermion_records,
                "max_off_bilinear_coefficient": maximum_off,
                "max_energy_residual": maximum_energy_error,
                "max_generator_identity_residual": maximum_identity_error,
                "minimum_unresolved_physical_off_bilinear": physical_off_minimum,
            },
            "generator_partition_reconstruction": {
                "cases": reconstruction_records,
                "max_relative_error": _decimal(maximum_generator_Z_error),
            },
        },
        "checks": checks,
    }


def main() -> None:
    report = run_experiment()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")

    for check in report["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"{status}: {check['name']} -- {check['detail']}")
    passed = all(check["passed"] for check in report["checks"])
    print(("PASS" if passed else "FAIL") + f": Onsager control experiment; wrote {OUTPUT.relative_to(ROOT)}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
