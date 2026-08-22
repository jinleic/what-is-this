"""Standalone acceptance checks for calibration-first Lee--Yang scaling."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from ising.lee_yang import (
    lee_yang_roots,
    positive_zero_angles,
    rational_field_polynomial_transfer,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results" / "lee_yang" / "scaling.json"
NOTE = ROOT / "notes" / "lee_yang_scaling.md"
DPS = 100


def coefficient_sha256(coefficients: list[int]) -> str:
    encoded = json.dumps(coefficients, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def test_lee_yang_scaling_acceptance() -> None:
    mp.mp.dps = DPS
    payload = json.loads(RESULT.read_text())
    assert set(payload) == {"provenance", "data", "checks"}
    assert payload["provenance"]["script"] == "experiments/e34_lee_yang_scaling.py"
    assert payload["provenance"]["precision"]["mpmath_dps"] == DPS
    assert payload["checks"] and all(item["passed"] for item in payload["checks"])

    data = payload["data"]
    calibration = data["two_dimensional_calibration"]
    assert calibration["external_control"]["sigma_2d_exact"] == "-1/6"
    assert calibration["external_control"]["role"] == "post-fit calibration target only"
    assert calibration["size_sequence"] == list(range(4, 16))
    assert calibration["largest_exact_lattice"]["L"] == 15
    assert calibration["largest_exact_lattice"]["n_sites"] == 225
    assert calibration["estimator"]["fit_sizes"] == list(range(8, 16))
    assert calibration["estimator"]["ranks_per_size"] == 6
    assert calibration["coupling"]["phase"] == "high_temperature (K<Kc)"
    assert mp.mpf(calibration["coupling"]["K"]) < mp.mpf(
        calibration["coupling"]["square_lattice_Kc_exact"]
    )

    # Recompute a complete exact calibration polynomial and its first root
    # independently through the public Lee--Yang API.
    coefficients_4 = rational_field_polynomial_transfer(
        (4, 4), numerator=2, denominator=3, periodic=False
    )
    stored_4 = calibration["exact_field_polynomials"]["4"]["A_k"]
    assert coefficients_4 == stored_4
    roots_4 = lee_yang_roots(coefficients_4, dps=DPS)
    theta_4 = positive_zero_angles(roots_4)[0]
    assert abs(theta_4 - mp.mpf(calibration["zero_data"]["4"]["theta_1"])) < mp.mpf(
        "1e-75"
    )

    # The largest exact polynomial is present literally, is palindromic, and
    # matches its independent content fingerprint.
    largest = calibration["exact_field_polynomials"]["15"]
    coefficients_15 = largest["A_k"]
    assert len(coefficients_15) == 226
    assert coefficients_15 == coefficients_15[::-1]
    assert coefficient_sha256(coefficients_15) == largest["coefficient_sha256"]
    roots_15 = calibration["zero_data"]["15"]
    assert roots_15["exact_reduced_real_root_count"] == 112
    assert roots_15["exact_reduced_degree"] == 112
    assert roots_15["positive_zero_count"] == 113
    assert mp.mpf(roots_15["theta_1"]) > 0
    for record in roots_15["positive_zero_angles"][:-1]:
        left, right = map(Fraction, record["t_isolating_interval"])
        assert -1 <= left < right <= 1

    # The scientific control is deliberately allowed to fail, but a failed
    # control must close the 3D-exponent gate without softening the decision.
    fit = calibration["fit_result"]
    target = -mp.mpf(1) / 6
    central_sigma = mp.mpf(fit["central"]["sigma"])
    assert not fit["passed"]
    assert abs(central_sigma - target) > mp.mpf("0.03")
    assert not fit["conditions"]["central_within_0.03_of_minus_one_sixth"]
    assert not fit["conditions"]["every_leave_one_size_out_within_0.03"]
    assert fit["conditions"]["leave_one_size_out_span_at_most_0.03"]
    assert fit["conditions"]["all_fits_converged_away_from_bounds"]
    loo = [mp.mpf(item["sigma"]) for item in fit["leave_one_size_out"]]
    assert len(loo) == 8
    assert max(loo) - min(loo) <= mp.mpf("0.03")
    drift = calibration["observed_drift_and_size_diagnostic"]
    assert drift["conditional_combined_projected_endpoint"] == 72
    resource = drift["conditional_projection_resource_scale"]
    assert resource["current_L15_final_layer_entries"] == 7_405_568
    assert resource["projected_endpoint_final_layer_entries"] > 10**24
    assert "unreliable on the available L<=15" in calibration["failure_scope"]

    three_dimensional = data["three_dimensional"]
    assert three_dimensional["edge_exponent"] is None
    assert three_dimensional["external_3d_exponent_used_or_reported"] is None
    assert three_dimensional["edge_exponent_status"].startswith("REFUSAL:")
    assert three_dimensional["certified_high_temperature_grid_ids"] == [
        "x_7_over_10",
        "x_2_over_3",
    ]

    cube_data = three_dimensional["exact_open_cube_zero_data"]
    assert set(cube_data) == {
        "x_7_over_10",
        "x_2_over_3",
        "x_9_over_14",
        "x_5_over_8",
        "x_3_over_5",
        "x_1_over_2",
    }
    for coupling in cube_data.values():
        sequence = coupling["theta_1_sequence"]
        assert [item["L"] for item in sequence] == [2, 3, 4]
        theta = [mp.mpf(item["theta_1"]) for item in sequence]
        assert theta[0] > theta[1] > theta[2] > 0
        assert coupling["finite_size_trend"]["theta_1_strictly_decreases_L2_to_L4"]
        assert coupling["finite_volume_edge_away_from_zero"][
            "certified_by_exact_P_at_z_1"
        ]
        for side in (2, 3, 4):
            record = coupling["sizes"][str(side)]
            coefficients = record["A_k"]
            assert coefficients == coefficients[::-1]
            assert record["exact_field_polynomial_at_z_1"] == sum(coefficients) > 0
            assert record["agrees_with_frozen_angles_to_1e-68"]
            assert mp.mpf(record["max_abs_angle_difference_from_frozen"]) < mp.mpf(
                "1e-68"
            )
            expected_intervals = min(8, record["positive_zero_count"] - 1)
            assert len(record["near_edge_density"]) == expected_intervals
            assert all(
                mp.mpf(item["density_1_over_N_delta_theta"]) > 0
                for item in record["near_edge_density"]
            )

    # The document contract requires the criterion to appear before any 3D
    # numbers and requires an explicit refusal when the control fails.
    note = NOTE.read_text()
    criterion_position = note.index("Estimator and stability criterion")
    result_position = note.index("Independent 3D exact-zero and density data")
    assert criterion_position < result_position
    assert "**The 2D control fails.**" in note
    assert "no 3D edge exponent was fitted" in note

    print("PASS: calibration-first Lee-Yang scaling artifact and exact sequences")
    print("PASS: 2D control failed and the 3D exponent gate remained closed")
    print("PASS: exact finite-volume theta_1 and near-edge density data through 4x4x4")


if __name__ == "__main__":
    try:
        test_lee_yang_scaling_acceptance()
    except Exception as error:
        print(f"FAIL: {type(error).__name__}: {error}")
        raise SystemExit(1) from error
