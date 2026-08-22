"""Exact finite-lattice Lee--Yang zeros, edge diagnostics, and scaling tests.

Run from the repository root with

    .venv/bin/python experiments/e04_lee_yang.py

The script writes ``results/lee_yang/lee_yang_analysis.json``.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import least_squares

from ising.exact_enumeration import dos_bonds
from ising.lattices import cubic, square
from ising.lee_yang import (
    circle_residual,
    field_dos_from_enumeration,
    field_polynomial_at_K,
    lee_yang_roots,
    near_edge_density,
    normalized_root_residual,
    positive_zero_angles,
    rational_field_polynomial,
    rational_field_polynomial_transfer,
    transfer_field_dos,
)
from ising.transfer_matrix import box_broken_bond_poly, torus_broken_bond_poly


DPS = 100
mp.mp.dps = DPS
REFERENCE_K = mp.mpf("0.221654626")
CIRCLE_TOLERANCE = mp.mpf("1e-40")
ROOT_TOLERANCE = mp.mpf("1e-70")
SCRIPT = "experiments/e04_lee_yang.py"
OUTPUT = Path(__file__).resolve().parents[1] / "results" / "lee_yang" / "lee_yang_analysis.json"

# These exact rational x=e^{-2K} points straddle the quoted 3D critical
# coupling.  The 9/14 point is near it but is not identified with it.  The
# published K_c is evaluated separately and is never used as a fit target.
RATIONAL_GRID = (
    ("x_7_over_10", 7, 10),
    ("x_2_over_3", 2, 3),
    ("x_9_over_14", 9, 14),
    ("x_5_over_8", 5, 8),
    ("x_3_over_5", 3, 5),
    ("x_1_over_2", 1, 2),
)


def mp_text(value, digits: int = 70) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def trim(values: list[int]) -> list[int]:
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def roots_record(coefficients: list, n_sites: int) -> dict:
    roots = lee_yang_roots(coefficients, dps=DPS)
    angles = positive_zero_angles(roots)
    densities = near_edge_density(roots, n_sites, intervals=4)
    return {
        "root_count": len(roots),
        "theta_1": mp_text(angles[0]),
        "positive_zero_angles": [mp_text(angle) for angle in angles],
        "near_edge_density": [
            {key: mp_text(value) for key, value in interval.items()}
            for interval in densities
        ],
        "max_abs_modulus_minus_one": mp_text(circle_residual(roots)),
        "max_normalized_polynomial_residual": mp_text(
            normalized_root_residual(coefficients, roots)
        ),
    }


def log_slope(points: list[tuple[int, float]]) -> dict:
    lengths = np.array([point[0] for point in points], dtype=np.float64)
    theta = np.array([point[1] for point in points], dtype=np.float64)
    slope, intercept = np.polyfit(np.log(lengths), np.log(theta), 1)
    prediction = intercept + slope * np.log(lengths)
    rms = float(np.sqrt(np.mean((np.log(theta) - prediction) ** 2)))
    return {
        "points": [{"L": int(L), "theta_1": repr(float(value))} for L, value in points],
        "effective_power_y_from_theta1_proportional_L_minus_y": repr(float(-slope)),
        "log_fit_rms": repr(rms),
        "warning": "Effective finite-size diagnostic only; an edge offset is not removed.",
    }


def _edge_fit(angle_sets: list[list[float]], sizes: list[int], j_max: int) -> dict:
    quantiles: list[float] = []
    angles: list[float] = []
    for zeros, n_sites in zip(angle_sets, sizes):
        for index, theta in enumerate(zeros[: min(j_max, len(zeros))], start=1):
            quantiles.append((index - 0.5) / n_sites)
            angles.append(theta)
    x = np.array(quantiles, dtype=np.float64)
    y = np.array(angles, dtype=np.float64)

    # Integrated density model:
    #   (j-1/2)/N = C' (theta_j-theta_edge)^(sigma+1),
    # or theta_j = theta_edge + C x^p with p=1/(sigma+1).
    def residual(parameters):
        edge, log_amplitude, power = parameters
        return edge + np.exp(log_amplitude) * x**power - y

    candidates = []
    for initial_power in (0.5, 0.9, 1.2):
        fit = least_squares(
            residual,
            (max(0.0, 0.5 * min(y)), math.log(max(y)), initial_power),
            bounds=((0.0, -20.0, 0.2), (math.pi, 20.0, 3.0)),
        )
        candidates.append(fit)
    fit = min(candidates, key=lambda result: float(np.dot(result.fun, result.fun)))
    edge, log_amplitude, power = map(float, fit.x)
    sigma = 1.0 / power - 1.0
    dof = max(1, len(y) - 3)
    residual_variance = float(np.dot(fit.fun, fit.fun) / dof)
    try:
        covariance = np.linalg.inv(fit.jac.T @ fit.jac) * residual_variance
        standard_error = np.sqrt(np.diag(covariance))
        edge_error = float(standard_error[0])
        sigma_error = float(standard_error[2] / power**2)
    except np.linalg.LinAlgError:
        edge_error = float("nan")
        sigma_error = float("nan")
    return {
        "theta_edge": edge,
        "theta_edge_formal_standard_error": edge_error,
        "sigma": sigma,
        "sigma_formal_standard_error": sigma_error,
        "power_1_over_sigma_plus_1": power,
        "rms_angle_residual": float(np.sqrt(np.mean(fit.fun**2))),
        "point_count": len(y),
        "j_max": j_max,
    }


def edge_fit_with_sensitivity(
    angle_sets: list[list[float]], sizes: list[int], dimension: int
) -> dict:
    central = _edge_fit(angle_sets, sizes, j_max=4)
    alternatives = []
    for j_max in (3, 4, 5):
        alternatives.append(_edge_fit(angle_sets, sizes, j_max))
        if len(angle_sets) > 2:
            alternatives.append(_edge_fit(angle_sets[1:], sizes[1:], j_max))
    edge_values = [fit["theta_edge"] for fit in alternatives]
    sigma_values = [fit["sigma"] for fit in alternatives]
    return {
        "model": "theta_j=theta_edge+C*((j-1/2)/N)^(1/(sigma+1))",
        "dimension": dimension,
        "central_fit": {key: repr(value) if isinstance(value, float) else value for key, value in central.items()},
        "sensitivity_range": {
            "theta_edge": [repr(min(edge_values)), repr(max(edge_values))],
            "sigma": [repr(min(sigma_values)), repr(max(sigma_values))],
            "variants": "j_max=3,4,5, with and without the smallest lattice",
        },
        "precision_note": "Roots use 100-digit mpmath; nonlinear regression and covariance use float64 scipy.",
    }


def main() -> None:
    mp.mp.dps = DPS
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    rational_grid = []
    for grid_id, numerator, denominator in RATIONAL_GRID:
        K = -mp.log(mp.mpf(numerator) / denominator) / 2
        rational_grid.append(
            {
                "id": grid_id,
                "numerator": numerator,
                "denominator": denominator,
                "x": f"{numerator}/{denominator}",
                "K": mp_text(K),
                "relative_to_reference_K": "below" if K < REFERENCE_K else "above",
            }
        )

    lattice_specs = []
    for shape in ((2, 2, 2), (3, 3, 2), (3, 3, 3)):
        for periodic in (False, True):
            label = f"cubic_{'x'.join(map(str, shape))}_{'periodic' if periodic else 'open'}"
            lattice_specs.append(
                (label, cubic(*shape, periodic=periodic), "joint_dos")
            )
    for side in (4, 5):
        for periodic in (False, True):
            label = f"square_{side}x{side}_{'periodic' if periodic else 'open'}"
            lattice_specs.append(
                (label, square(side, side, periodic=periodic), "joint_dos")
            )

    exact_polynomials: dict[str, dict] = {}
    zeros: dict[str, dict] = {}
    field_dos_by_label: dict[str, np.ndarray] = {}
    for label, lattice, engine in lattice_specs:
        field_dos = field_dos_from_enumeration(lattice)
        field_dos_by_label[label] = field_dos
        palindrome = np.array_equal(field_dos, field_dos[::-1])
        broken_from_field = [int(value) for value in field_dos.sum(axis=0)]
        broken_from_dos = [int(value) for value in dos_bonds(lattice)][::-1]
        z_one_ok = broken_from_field == broken_from_dos
        check(
            f"{label}_exact_symmetry_and_z1",
            palindrome and z_one_ok,
            f"c[k,q]=c[N-k,q]: {palindrome}; sum_k c[k,q]=dos_bonds[N_b-q]: {z_one_ok}",
        )
        exact_polynomials[label] = {
            "shape": list(lattice.shape),
            "periodic": list(lattice.periodic),
            "n_sites": lattice.n_sites,
            "n_bonds": lattice.n_bonds,
            "engine": engine,
            "axes": ["k_down", "q_unsatisfied_bonds"],
            "coefficients_c_k_q": [[int(value) for value in row] for row in field_dos],
        }

    for side in (6,):
        for periodic in (False, True):
            label = f"square_{side}x{side}_{'periodic' if periodic else 'open'}"
            field_dos = transfer_field_dos((side, side), periodic=periodic)
            field_dos_by_label[label] = field_dos
            zero_field = (
                torus_broken_bond_poly((side, side))
                if periodic
                else box_broken_bond_poly((side, side))
            )
            z_one = trim([int(value) for value in field_dos.sum(axis=0)])
            palindrome = np.array_equal(field_dos, field_dos[::-1])
            z_one_ok = z_one == zero_field
            check(
                f"{label}_exact_symmetry_and_z1",
                palindrome and z_one_ok,
                f"c[k,q]=c[N-k,q]: {palindrome}; agrees with independent zero-field transfer: {z_one_ok}",
            )
            n_bonds = field_dos.shape[1] - 1
            exact_polynomials[label] = {
                "shape": [side, side],
                "periodic": [periodic, periodic],
                "n_sites": side * side,
                "n_bonds": n_bonds,
                "engine": "bivariate_transfer_matrix",
                "axes": ["k_down", "q_unsatisfied_bonds"],
                "coefficients_c_k_q": [[int(value) for value in row] for row in field_dos],
            }

    # All exact c[k,q] polynomials are evaluated on the rational grid, and
    # separately at the quoted critical-coupling benchmark.
    for label, polynomial in exact_polynomials.items():
        n_sites = polynomial["n_sites"]
        field_dos = field_dos_by_label[label]
        records: dict[str, dict] = {}
        for grid_point in rational_grid:
            coefficients = rational_field_polynomial(
                field_dos, grid_point["numerator"], grid_point["denominator"]
            )
            record = roots_record(coefficients, n_sites)
            record.update(
                {
                    "K": grid_point["K"],
                    "x": grid_point["x"],
                    "coefficient_arithmetic": "exact integers after a common denominator^n_b scaling",
                }
            )
            records[grid_point["id"]] = record
        reference_coefficients = field_polynomial_at_K(field_dos, REFERENCE_K, dps=DPS)
        reference_record = roots_record(reference_coefficients, n_sites)
        reference_record.update(
            {
                "K": mp_text(REFERENCE_K),
                "x": mp_text(mp.exp(-2 * REFERENCE_K)),
                "coefficient_arithmetic": "100-digit evaluation of exact integer c[k,q]",
                "role": "published coupling used only as a falsification/evaluation point, never fitted",
            }
        )
        records["reference_K_0.221654626"] = reference_record
        zeros[label] = records

        worst_circle = max(
            mp.mpf(record["max_abs_modulus_minus_one"]) for record in records.values()
        )
        correct_counts = all(record["root_count"] == n_sites for record in records.values())
        good_residuals = all(
            mp.mpf(record["max_normalized_polynomial_residual"]) < ROOT_TOLERANCE
            for record in records.values()
        )
        check(
            f"{label}_roots",
            worst_circle < CIRCLE_TOLERANCE and correct_counts and good_residuals,
            f"7 couplings; max radial error={mp_text(worst_circle, 12)}; root count={n_sites}",
        )

    # Exact 64-site integer polynomials at every rational grid point.  A full
    # q axis would be too large; modular transfer evaluates x exactly while
    # retaining all 65 field coefficients.
    large_cube_polynomials: dict[str, dict] = {}
    large_cube_zeros: dict[str, dict] = {}
    for grid_point in rational_grid:
        numerator = grid_point["numerator"]
        denominator = grid_point["denominator"]
        coefficients = rational_field_polynomial_transfer(
            (4, 4, 4), numerator, denominator, periodic=False
        )
        record = roots_record(coefficients, 64)
        record.update(
            {
                "K": grid_point["K"],
                "x": grid_point["x"],
                "coefficient_arithmetic": "exact modular transfer plus CRT",
            }
        )
        large_cube_zeros[grid_point["id"]] = record
        large_cube_polynomials[grid_point["id"]] = {
            "x": grid_point["x"],
            "K": grid_point["K"],
            "common_scale": f"{denominator}^144",
            "A_k": coefficients,
        }
        radial = mp.mpf(record["max_abs_modulus_minus_one"])
        residual = mp.mpf(record["max_normalized_polynomial_residual"])
        check(
            f"cubic_4x4x4_open_{grid_point['id']}_roots",
            coefficients == coefficients[::-1]
            and record["root_count"] == 64
            and radial < CIRCLE_TOLERANCE
            and residual < ROOT_TOLERANCE,
            f"exact degree 64; radial error={mp_text(radial, 12)}",
        )

    # theta_1 finite-size diagnostics.  At rational points the open cubic
    # family reaches L=4; at the benchmark only L=2,3 are available.
    scaling: dict[str, dict] = {}
    family_specs = {
        "cubic_open": (("cubic_2x2x2_open", 2), ("cubic_3x3x3_open", 3)),
        "cubic_periodic": (
            ("cubic_2x2x2_periodic", 2),
            ("cubic_3x3x3_periodic", 3),
        ),
        "square_open_control": (
            ("square_4x4_open", 4),
            ("square_5x5_open", 5),
            ("square_6x6_open", 6),
        ),
        "square_periodic_control": (
            ("square_4x4_periodic", 4),
            ("square_5x5_periodic", 5),
            ("square_6x6_periodic", 6),
        ),
    }
    grid_ids = [point["id"] for point in rational_grid] + ["reference_K_0.221654626"]
    for family, members in family_specs.items():
        family_records = {}
        for grid_id in grid_ids:
            points = [
                (L, float(mp.mpf(zeros[label][grid_id]["theta_1"])))
                for label, L in members
            ]
            if family == "cubic_open" and grid_id in large_cube_zeros:
                points.append(
                    (4, float(mp.mpf(large_cube_zeros[grid_id]["theta_1"])))
                )
            family_records[grid_id] = log_slope(points)
        scaling[family] = family_records

    # High-temperature edge fits must not take ANY literature number as an input.  The previous
    # version selected the disordered-phase grid points by `K < REFERENCE_K` with
    # REFERENCE_K = 0.221654626 (Ferrenberg-Xu-Landau).  Even though it was not an adjustable
    # parameter, it made a published benchmark an input to the fit selection -- flagged by
    # reports/adversarial_audit.md and fixed here.
    #
    # Replacement: this repository's OWN certified rigorous lower bound on K_c, proved in
    # proofs/kc_bounds.md (SAW / connective-constant bound, K_c >= atanh(c_14^{-1/14}) with the
    # exact SAW count c_14 = 4468911678 computed in-repo).  `K < CERTIFIED_KC_LOWER` therefore
    # RIGOROUSLY implies `K < K_c`, which is strictly stronger than the old heuristic and uses
    # no external number at all.
    CERTIFIED_KC_LOWER = mp.mpf("0.2074277114992039908436804465100887760027")
    edge_fits_3d: dict[str, dict] = {}
    edge_fits_2d_control: dict[str, dict] = {}
    high_temperature_ids = [
        point["id"]
        for point in rational_grid
        if mp.mpf(point["K"]) < CERTIFIED_KC_LOWER
    ]
    for grid_id in high_temperature_ids:
        cubic_sets = [
            [float(mp.mpf(value)) for value in zeros[label][grid_id]["positive_zero_angles"]]
            for label in ("cubic_2x2x2_open", "cubic_3x3x3_open")
        ]
        cubic_sets.append(
            [
                float(mp.mpf(value))
                for value in large_cube_zeros[grid_id]["positive_zero_angles"]
            ]
        )
        edge_fits_3d[grid_id] = edge_fit_with_sensitivity(
            cubic_sets, [8, 27, 64], dimension=3
        )

        square_sets = [
            [float(mp.mpf(value)) for value in zeros[label][grid_id]["positive_zero_angles"]]
            for label in (
                "square_4x4_periodic",
                "square_5x5_periodic",
                "square_6x6_periodic",
            )
        ]
        edge_fits_2d_control[grid_id] = edge_fit_with_sensitivity(
            square_sets, [16, 25, 36], dimension=2
        )

    data = {
        "convention": {
            "k": "number of down spins",
            "q": "number of unsatisfied bonds",
            "x": "exp(-2K)",
            "z": "exp(-2H_f)",
            "identity": "Z=exp(K*n_b+H_f*N)*sum_{k,q} c[k,q] x^q z^k",
            "spin_reversal": "c[k,q]=c[N-k,q], hence A_k=A_{N-k}",
        },
        "coupling_grid": rational_grid
        + [
            {
                "id": "reference_K_0.221654626",
                "K": mp_text(REFERENCE_K),
                "x": mp_text(mp.exp(-2 * REFERENCE_K)),
                "role": "falsification/evaluation benchmark only; not fitted",
            }
        ],
        "exact_field_dos_polynomials": exact_polynomials,
        "zeros": zeros,
        "largest_lattice": {
            "shape": [4, 4, 4],
            "periodic": [False, False, False],
            "n_sites": 64,
            "n_bonds": 144,
            "engine": "exact rational modular transfer with a field-polynomial axis and CRT",
            "polynomials": large_cube_polynomials,
            "zeros": large_cube_zeros,
        },
        "theta_1_finite_size_scaling": scaling,
        "edge_analysis": {
            "definition": "g(theta) proportional to (theta-theta_edge)^sigma on the positive arc",
            "three_dimensional_open_cubes": edge_fits_3d,
            "two_dimensional_periodic_control": edge_fits_2d_control,
            "literature_falsification_targets_not_used_in_fit": {
                "sigma_3d_representative_range": ["0.07", "0.09"],
                "sigma_2d_exact": "-0.1666666666666666666666666666666666667",
            },
            "assessment": (
                "The 3D L=2,3,4 fits give positive edges for all sampled K<K_c, "
                "but sigma is strongly fit-window and boundary sensitive.  The same small-size "
                "procedure does not recover the exact 2D sigma=-1/6 control.  Therefore these "
                "lattices locate only provisional edges and do not support a defensible "
                "asymptotic 3D sigma estimate."
            ),
        },
    }
    payload = {
        "provenance": {
            "script": SCRIPT,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "precision": DPS,
            "root_method": "Chebyshev half-polynomial isolation plus 100-digit real and complex Newton refinement",
            "fit_method": "float64 scipy least_squares on 100-digit zero angles",
        },
        "data": data,
        "checks": checks,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    passed = sum(item["passed"] for item in checks)
    total = len(checks)
    worst = max(
        mp.mpf(record["max_abs_modulus_minus_one"])
        for lattice_records in list(zeros.values()) + [large_cube_zeros]
        for record in lattice_records.values()
    )
    if passed == total:
        print(f"PASS: {passed}/{total} exact and numerical Lee-Yang checks")
        print(f"PASS: largest lattice = 4x4x4 open (64 sites)")
        print(f"PASS: worst max ||z|-1| = {mp_text(worst, 12)} < 1e-40")
        print(f"WROTE: {OUTPUT}")
    else:
        print(f"FAIL: {passed}/{total} checks passed; see {OUTPUT}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
