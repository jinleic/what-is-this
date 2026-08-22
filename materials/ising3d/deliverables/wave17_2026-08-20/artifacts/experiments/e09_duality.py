"""Exact finite-volume verification of 3D Ising--Z2 gauge and 2D KW duality.

Run from the repository root with::

    .venv/bin/python experiments/e09_duality.py

The script uses no floating point.  It writes every density of states and every polynomial needed
for independent inspection to ``results/duality/duality.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from ising.duality import (
    closed_surface_polynomial,
    cubical_complex,
    dual_ising_sector_polynomials,
    flat_gauge_orbit_count,
    gauge_orbit_size,
    gauge_satisfied_plaquette_dos,
    homology_dimensions,
    ising_torus_sector_polynomials_2d,
    kw_character,
    signed_high_temperature_polynomial,
    surface_polynomial_from_gauge_dos,
    surface_sector_basis,
    surface_sector_representatives,
)
from ising.exact_enumeration import dos_bonds
from ising.lattices import cubic
from ising.transfer_matrix import torus_broken_bond_poly


THREE_DIMENSIONAL_CASES = (
    ((2, 2, 2), (False, False, False)),
    ((3, 2, 2), (False, False, False)),
    ((2, 2, 2), (True, False, False)),
    ((2, 2, 2), (True, True, False)),
    ((2, 2, 2), (True, True, True)),
)
TWO_DIMENSIONAL_CASES = ((4, 4), (5, 5))
SECTOR_LABELS_2D = ("PP", "AP_x", "AP_y", "AP_x_AP_y")


def _sum_polynomials(polynomials: tuple[tuple[int, ...], ...]) -> list[int]:
    degree = max(len(polynomial) for polynomial in polynomials)
    return [
        sum(
            polynomial[power] if power < len(polynomial) else 0
            for polynomial in polynomials
        )
        for power in range(degree)
    ]


def _signed_sum_polynomials(
    polynomials: tuple[tuple[int, ...], ...], signs: tuple[int, ...]
) -> list[int]:
    degree = max(len(polynomial) for polynomial in polynomials)
    return [
        sum(
            signs[index] * (polynomial[power] if power < len(polynomial) else 0)
            for index, polynomial in enumerate(polynomials)
        )
        for power in range(degree)
    ]


def main() -> None:
    checks: list[dict[str, object]] = []
    data_3d: list[dict[str, object]] = []
    data_2d: list[dict[str, object]] = []

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        print(f"{'PASS' if passed else 'FAIL'}: {name} -- {detail}", flush=True)

    for shape, periodic in THREE_DIMENSIONAL_CASES:
        complex_ = cubical_complex(shape, periodic)
        case_name = f"3d_{'x'.join(map(str, shape))}_{''.join('P' if p else 'O' for p in periodic)}"
        homology = homology_dimensions(complex_)
        gauge_dos = gauge_satisfied_plaquette_dos(complex_)
        gauge_surface = surface_polynomial_from_gauge_dos(complex_, gauge_dos)
        direct_surface = closed_surface_polynomial(complex_)
        sector_basis = surface_sector_basis(complex_)
        representatives = surface_sector_representatives(complex_)
        sector_polynomials = dual_ising_sector_polynomials(complex_, representatives)
        dual_sum = _sum_polynomials(sector_polynomials)
        kernel_factor = 1 << homology["b3"]
        expected_dual_sum = [kernel_factor * coefficient for coefficient in direct_surface]
        expected_flat_configurations = 1 << (complex_.n_links - homology["rank_d2"])
        expected_surface_count = 1 << (
            complex_.n_plaquettes - homology["rank_d2"]
        )

        record(
            f"{case_name}.gauge_enumeration",
            sum(gauge_dos) == 1 << complex_.n_links,
            f"sum DOS={sum(gauge_dos)}=2^{complex_.n_links}",
        )
        record(
            f"{case_name}.gauge_HT_surfaces",
            gauge_surface == direct_surface,
            "DOS binomial transform equals independent GF(2) closed-surface enumerator",
        )
        record(
            f"{case_name}.dual_sector_identity",
            dual_sum == expected_dual_sum,
            f"sum_h Q_h = 2^{homology['b3']} S coefficient by coefficient",
        )
        record(
            f"{case_name}.topological_multiplicities",
            (
                len(representatives) == 1 << homology["b2"]
                and sum(direct_surface) == expected_surface_count
                and gauge_dos[complex_.n_plaquettes] == expected_flat_configurations
                and expected_flat_configurations
                == gauge_orbit_size(complex_) * flat_gauge_orbit_count(complex_)
            ),
            (
                f"dual sectors=2^{homology['b2']}; closed surfaces={expected_surface_count}; "
                f"flat links=2^{homology['rank_d1']}*2^{homology['b1']}"
            ),
        )
        if shape == (2, 2, 2) and periodic == (True, True, True):
            independent_ising_dos = dos_bonds(cubic(*shape, periodic=True))
            independent_broken_bond_polynomial = tuple(
                int(value) for value in reversed(independent_ising_dos)
            )
            record(
                f"{case_name}.independent_ising_enumerator",
                sector_polynomials[0] == independent_broken_bond_polynomial,
                "periodic dual sector agrees with ising.exact_enumeration",
            )

        data_3d.append(
            {
                "shape": list(shape),
                "periodic": list(periodic),
                "description": complex_.describe(),
                "cell_counts": {
                    "n0_vertices": complex_.n_vertices,
                    "n1_links": complex_.n_links,
                    "n2_plaquettes": complex_.n_plaquettes,
                    "n3_cubes": complex_.n_cubes,
                },
                "homology_F2": homology,
                "gauge_orbit_size": gauge_orbit_size(complex_),
                "flat_gauge_orbit_count": flat_gauge_orbit_count(complex_),
                "raw_gauge_satisfied_plaquette_dos": list(gauge_dos),
                "closed_surface_polynomial": list(direct_surface),
                "surface_sector_basis_hex": [hex(value) for value in sector_basis],
                "surface_sector_representatives_hex": [
                    hex(value) for value in representatives
                ],
                "dual_ising_broken_bond_polynomials": [
                    list(polynomial) for polynomial in sector_polynomials
                ],
                "integer_identity": {
                    "left_sum_dual_sector_polynomials": dual_sum,
                    "right_2_pow_b3_times_surface_polynomial": expected_dual_sum,
                },
                "partition_function_prefactor": {
                    "coupling_map": "exp(-2*K_star) = tanh(K_g)",
                    "A": "sqrt(sinh(2*K_g)/2)",
                    "raw_identity": (
                        "Z_gauge = 2^(n1-b3) * A^n2 * "
                        "sum_[h in H2] Z_Ising_dual^[h]"
                    ),
                    "power_of_two_n1_minus_b3": complex_.n_links - homology["b3"],
                    "power_of_A_n2": complex_.n_plaquettes,
                    "gauge_orbit_normalized_power_of_two": (
                        complex_.n_links - homology["b3"] - homology["rank_d1"]
                    ),
                },
            }
        )

    for shape in TWO_DIMENSIONAL_CASES:
        n_sites = shape[0] * shape[1]
        n_edges = 2 * n_sites
        sector_polynomials = ising_torus_sector_polynomials_2d(shape)
        high_temperature_polynomials = tuple(
            signed_high_temperature_polynomial(polynomial, n_sites)
            for polynomial in sector_polynomials
        )
        row_results: list[dict[str, object]] = []
        all_rows_pass = True
        for alpha in range(4):
            signs = tuple(kw_character(alpha, beta) for beta in range(4))
            right = _signed_sum_polynomials(sector_polynomials, signs)
            left = [2 * coefficient for coefficient in high_temperature_polynomials[alpha]]
            row_pass = left == right
            all_rows_pass = all_rows_pass and row_pass
            row_results.append(
                {
                    "primal_sector": SECTOR_LABELS_2D[alpha],
                    "character_signs": list(signs),
                    "left_2_times_high_temperature_polynomial": left,
                    "right_signed_sum_dual_low_temperature_polynomials": right,
                    "passed": row_pass,
                }
            )
        cycle_space_dimension = n_edges - n_sites + 1
        invariants_pass = (
            high_temperature_polynomials[0][0] == 1
            and sum(high_temperature_polynomials[0]) == 1 << cycle_space_dimension
            and all(sum(polynomial) == 1 << n_sites for polynomial in sector_polynomials)
        )
        case_name = f"2d_{shape[0]}x{shape[1]}_torus"
        record(
            f"{case_name}.four_sector_KW",
            all_rows_pass,
            "2 R_alpha(v) = sum_beta chi(alpha,beta) Q_beta(v), all four rows",
        )
        record(
            f"{case_name}.enumeration_invariants",
            invariants_pass,
            f"each sector has 2^{n_sites} spins; P_PP(1)=2^{cycle_space_dimension}",
        )
        independent_tm_polynomial = torus_broken_bond_poly(shape)
        independent_tm_pass = (
            list(sector_polynomials[0][: len(independent_tm_polynomial)])
            == independent_tm_polynomial
            and not any(sector_polynomials[0][len(independent_tm_polynomial) :])
        )
        record(
            f"{case_name}.independent_transfer_matrix",
            independent_tm_pass,
            "untwisted row-sector polynomial agrees with ising.transfer_matrix",
        )
        data_2d.append(
            {
                "shape": list(shape),
                "sector_order": list(SECTOR_LABELS_2D),
                "n_sites": n_sites,
                "n_edges": n_edges,
                "broken_bond_polynomials": [
                    list(polynomial) for polynomial in sector_polynomials
                ],
                "signed_high_temperature_polynomials": [
                    list(polynomial) for polynomial in high_temperature_polynomials
                ],
                "independent_untwisted_transfer_matrix_polynomial": independent_tm_polynomial,
                "four_sector_integer_transform": row_results,
                "partition_function_identity_PP": (
                    "Z_PP(K) = (sinh(2*K))^N / 2 * "
                    "sum_beta Z_beta(K_star), exp(-2*K_star)=tanh(K)"
                ),
            }
        )

    output = {
        "provenance": {
            "script": "experiments/e09_duality.py",
            "interpreter": sys.executable,
            "python_version": sys.version.split()[0],
            "numpy_version": np.__version__,
            "arithmetic": (
                "exact Python integers; overflow-bounded numpy.int64 enumeration; no floating point"
            ),
            "method": (
                "raw link brute force -> exact binomial transform, independent GF(2) surface "
                "enumeration, and independent dual-spin/row-transfer enumeration"
            ),
        },
        "data": {
            "three_dimensional": data_3d,
            "two_dimensional_control": data_2d,
        },
        "checks": checks,
    }
    output_path = Path(__file__).resolve().parents[1] / "results" / "duality" / "duality.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")

    passed_count = sum(bool(check["passed"]) for check in checks)
    print(f"wrote {output_path.relative_to(output_path.parents[2])}")
    print(f"SUMMARY: {passed_count}/{len(checks)} checks passed")
    if passed_count != len(checks):
        raise AssertionError("one or more exact duality checks failed")
    print("PASS")


if __name__ == "__main__":
    main()
