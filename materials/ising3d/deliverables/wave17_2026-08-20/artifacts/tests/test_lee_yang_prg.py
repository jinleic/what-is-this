"""Standalone checks for the calibration-frozen Lee--Yang PRG barrier."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from ising.transfer_matrix import layer_bonds


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "lee_yang" / "prg.json"
SOURCE = ROOT / "results" / "lee_yang" / "scaling.json"
PROOF = ROOT / "proofs" / "lee_yang_prg.md"
DPS = 120
CENTERS = tuple(range(5, 15))
THETA_RADIUS = Fraction(1, 10**60)


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ratio_value(power: mp.mpf, sizes: tuple[int, int, int]) -> mp.mpf:
    left, middle, right = map(mp.mpf, sizes)
    return (left ** (-power) - middle ** (-power)) / (
        middle ** (-power) - right ** (-power)
    )


def independent_prg_power(sizes: tuple[int, int, int], values: tuple[mp.mpf, ...]) -> mp.mpf:
    observed = (values[0] - values[1]) / (values[1] - values[2])
    lower = mp.mpf("1e-12")
    upper = mp.mpf(16)
    f_lower = ratio_value(lower, sizes) - observed
    f_upper = ratio_value(upper, sizes) - observed
    assert f_lower * f_upper < 0
    for _ in range(600):
        midpoint = (lower + upper) / 2
        f_midpoint = ratio_value(midpoint, sizes) - observed
        if f_lower * f_midpoint <= 0:
            upper = midpoint
            f_upper = f_midpoint
        else:
            lower = midpoint
            f_lower = f_midpoint
        if upper - lower < mp.mpf("1e-100"):
            break
    return (lower + upper) / 2


def fit_sigma(centers: tuple[int, ...], powers: dict[int, mp.mpf]) -> mp.mpf:
    rows = [[mp.mpf(1), mp.mpf(1) / center, mp.mpf(1) / center**2] for center in centers]
    normal = mp.matrix(
        [[sum(row[i] * row[j] for row in rows) for j in range(3)] for i in range(3)]
    )
    rhs = mp.matrix(
        [sum(row[i] * powers[center] for row, center in zip(rows, centers)) for i in range(3)]
    )
    coefficients = mp.lu_solve(normal, rhs)
    return mp.mpf(2) / coefficients[0] - 1


def invert_matrix(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    work = [row[:] + [Fraction(int(i == j)) for j in range(size)] for i, row in enumerate(matrix)]
    for pivot_index in range(size):
        pivot_row = next(row for row in range(pivot_index, size) if work[row][pivot_index])
        work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
        pivot = work[pivot_index][pivot_index]
        work[pivot_index] = [value / pivot for value in work[pivot_index]]
        for row in range(size):
            if row != pivot_index:
                factor = work[row][pivot_index]
                work[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(work[row], work[pivot_index])
                ]
    return [row[size:] for row in work]


def independent_barrier_coefficients(theta: list[Fraction], power: int) -> tuple[list[Fraction], list[Fraction]]:
    sizes = (2, 3, 4)
    edge = Fraction(1, 50)
    matrix = [[Fraction(1), Fraction(1, size), Fraction(1, size**2)] for size in sizes]
    inverse = invert_matrix(matrix)
    rhs = [Fraction(size**power) * (value - edge) for size, value in zip(sizes, theta)]
    coefficients = [
        sum((inverse[row][column] * rhs[column] for column in range(3)), Fraction(0))
        for row in range(3)
    ]
    radii = [
        THETA_RADIUS
        * sum(
            (abs(inverse[row][column]) * sizes[column] ** power for column in range(3)),
            Fraction(0),
        )
        for row in range(3)
    ]
    return coefficients, radii


def fixed_field_value(theta: mp.mpf, length: int, periodic_cross: tuple[bool, bool]) -> mp.mpf:
    bonds = layer_bonds((2, 2), periodic_cross)
    broken = [
        sum(((state >> left) ^ (state >> right)) & 1 for left, right in bonds)
        for state in range(16)
    ]
    down = [4 - state.bit_count() for state in range(16)]
    x = mp.mpf(2) / 3
    z = mp.exp(mp.j * theta)
    weights = [x ** broken[state] * z ** down[state] for state in range(16)]
    vector = weights[:]
    for _ in range(1, length):
        transformed = vector
        for bit in range(4):
            transformed = [
                transformed[state] + x * transformed[state ^ (1 << bit)]
                for state in range(16)
            ]
        vector = [weights[state] * transformed[state] for state in range(16)]
        scale = max(abs(value) for value in vector)
        vector = [value / scale for value in vector]
    phased = mp.exp(-mp.j * 2 * length * theta) * sum(vector)
    assert abs(mp.im(phased)) < mp.mpf("1e-100")
    return mp.re(phased)


def test_lee_yang_prg_acceptance() -> None:
    mp.mp.dps = DPS
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e73_lee_yang_prg.py"
    assert payload["provenance"]["three_dimensional_benchmark_exponent_used"] is False
    assert payload["checks"] and all(check["passed"] for check in payload["checks"])

    data = payload["data"]
    protocol = data["frozen_protocol"]
    assert protocol["selected_using"] == "two-dimensional control only"
    assert protocol["three_dimensional_benchmark_used"] is False
    assert canonical_sha256(protocol) == data["frozen_protocol_sha256"]

    # Independently reconstruct every local ratio exponent and every member of
    # the predeclared robust calibration hull.
    theta = {
        int(size): mp.mpf(record["theta_1"])
        for size, record in source["data"]["two_dimensional_calibration"]["zero_data"].items()
    }
    powers = {
        center: independent_prg_power(
            (center - 1, center, center + 1),
            (theta[center - 1], theta[center], theta[center + 1]),
        )
        for center in CENTERS
    }
    calibration = data["two_dimensional_calibration"]
    for center in CENTERS:
        stored = mp.mpf(calibration["local_prg"][str(center)]["power_p"])
        assert abs(stored - powers[center]) < mp.mpf("1e-80")

    central = fit_sigma(CENTERS, powers)
    leave_one_out = [
        fit_sigma(tuple(center for center in CENTERS if center != omitted), powers)
        for omitted in CENTERS
    ]
    suffixes = [
        fit_sigma(tuple(range(start, 15)), powers) for start in range(5, 11)
    ]
    estimates = [central, *leave_one_out, *suffixes]
    independent_interval = [min(estimates), max(estimates)]
    stored_interval = list(map(mp.mpf, calibration["empirical_sigma_interval"]))
    assert abs(central - mp.mpf(calibration["central_overdetermined_fit"]["sigma"])) < mp.mpf("1e-80")
    assert max(abs(left - right) for left, right in zip(independent_interval, stored_interval)) < mp.mpf("1e-80")
    target = -mp.mpf(1) / 6
    assert independent_interval[0] <= target <= independent_interval[1]
    assert calibration["passed"] and calibration["target_inside_empirical_interval"]

    application = data["three_dimensional_application"]
    assert application["protocol_sha256"] == data["frozen_protocol_sha256"]
    assert application["available_isotropic_cube_sizes"] == [2, 3, 4]
    assert application["required_isotropic_cube_sizes"] == list(range(4, 16))
    assert application["required_consecutive_isotropic_size_count"] == 12
    assert not application["protocol_instantiable"]
    assert application["edge_exponent_estimate"] is None
    assert application["external_3d_exponent_used_for_selection_or_comparison"] is None
    assert all(
        item["tag"] == "UNCONTROLLED_DIAGNOSTIC_NOT_AN_ESTIMATE"
        for item in application["three_size_diagnostics"].values()
    )

    # Reconstruct both exact rational interpolation certificates, including the
    # angle-error coefficient boxes and the global positivity discriminants.
    barrier = data["exact_identifiability_barrier"]
    assert barrier["claim_tag"] == "THEOREM"
    source_cubes = source["data"]["three_dimensional"]["exact_open_cube_zero_data"]
    for grid_id, coupling in barrier["couplings"].items():
        assert set(coupling["theta_angle_certificates"]) == {"2", "3", "4"}
        assert all(item["passed"] for item in coupling["theta_angle_certificates"].values())
        theta_centers = [
            Fraction(item["theta_1"])
            for item in source_cubes[grid_id]["theta_1_sequence"]
        ]
        for model in coupling["models"]:
            power = model["p"]
            coefficients, radii = independent_barrier_coefficients(theta_centers, power)
            stored_coefficients = [
                Fraction(model["nominal_coefficients_exact"][name])
                for name in ("c0", "c1", "c2")
            ]
            assert coefficients == stored_coefficients
            for index, name in enumerate(("c0", "c1", "c2")):
                stored_bounds = list(
                    map(Fraction, model["coefficient_intervals_covering_exact_angles"][name])
                )
                assert stored_bounds == [
                    coefficients[index] - radii[index],
                    coefficients[index] + radii[index],
                ]
            c0, c1, c2 = coefficients
            assert c0 > 0 and c2 > 0
            assert c1 * c1 - 4 * c0 * c2 < 0
            assert (power + 1) ** 2 * c1 * c1 - 4 * power * (power + 2) * c0 * c2 < 0
            assert model["nominal_exact_fit_residuals"] == ["0", "0", "0"]
            assert model[
                "all_coefficient_box_models_above_edge_and_strictly_decreasing_for_L_positive"
            ]
        assert [model["sigma_exact_d3"] for model in coupling["models"]] == ["1/2", "-1/4"]
        assert coupling["sigma_separation_exact"] == "3/4"
        assert coupling["both_models_certified"]

    tube = data["fixed_field_strip_cylinder_diagnostic"]
    assert tube["claim_tag"] == "EMPIRICAL_DIAGNOSTIC"
    boundary_map = {
        "open_2x2_bar": (False, False),
        "one_transverse_periodic_cylinder": (True, False),
    }
    for name, periodic_cross in boundary_map.items():
        sequence = tube["sequences"][name]
        assert [record["length"] for record in sequence] == [2, 4, 8, 16, 32, 64]
        angles = [mp.mpf(record["theta_1"]) for record in sequence]
        assert all(left > right > 0 for left, right in zip(angles, angles[1:]))
        assert all(mp.mpf(record["bracket_width_160_dps"]) <= mp.mpf("1e-75") for record in sequence)
        assert all(
            mp.mpf(record["120_vs_160_dps_midpoint_difference"]) < mp.mpf("1e-60")
            for record in sequence
        )
        lower, upper = map(mp.mpf, sequence[0]["theta_sign_bracket_160_dps"])
        assert fixed_field_value(lower, 2, periodic_cross) * fixed_field_value(upper, 2, periodic_cross) < 0
    open_theta_2 = mp.mpf(tube["sequences"]["open_2x2_bar"][0]["theta_1"])
    cube_theta_2 = mp.mpf(source_cubes["x_2_over_3"]["theta_1_sequence"][0]["theta_1"])
    assert abs(open_theta_2 - cube_theta_2) < mp.mpf("1e-60")
    assert "cannot be inserted" in tube["geometric_warning"]

    proof = PROOF.read_text(encoding="utf-8")
    assert "[EMPIRICAL CALIBRATION]" in proof
    assert "[THEOREM]" in proof
    assert "No thermodynamic 3D edge exponent" in proof
    assert proof.index("Frozen 2D-only protocol") < proof.index("3D application")

    print("PASS")


if __name__ == "__main__":
    try:
        test_lee_yang_prg_acceptance()
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise SystemExit(1) from error
