from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.e06_falsify_claims import (  # noqa: E402
    PRECISION,
    box_broken_bond_magnetization_poly,
    claimed_degang_free_energy_m,
    degang_critical_data,
    exact_ht_bulk_v_series,
    falsify_free_energy,
    low_temperature_magnetization_series,
    run_analysis,
    winding_projected_torus_polynomial,
    zhang_critical_data,
)
from ising.transfer_matrix import box_broken_bond_poly  # noqa: E402

def _brute_plus_joint(shape: tuple[int, int, int]) -> list[list[int]]:
    coordinates = list(product(*(range(side) for side in shape)))
    site_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    counts: dict[tuple[int, int], int] = {}
    for mask in range(1 << len(coordinates)):
        down = [(mask >> index) & 1 for index in range(len(coordinates))]
        broken = 0
        for coordinate, index in site_index.items():
            for axis, side in enumerate(shape):
                neighbor = list(coordinate)
                neighbor[axis] += 1
                if neighbor[axis] < side:
                    broken += down[index] != down[site_index[tuple(neighbor)]]
                else:
                    broken += down[index]
                if coordinate[axis] == 0:
                    broken += down[index]
        key = (broken, sum(down))
        counts[key] = counts.get(key, 0) + 1
    result = [
        [counts.get((broken, n_down), 0) for n_down in range(len(coordinates) + 1)]
        for broken in range(max(key[0] for key in counts) + 1)
    ]
    return result


def test_zhang_critical_condition() -> None:
    data = zhang_critical_data()
    assert mp.mpf(data["root_closed_form_error"]) < mp.mpf("1e-70")
    assert mp.mpf(data["absolute_discrepancy"]) > mp.mpf("0.018")
    assert mp.mpf(data["relative_discrepancy"]) > mp.mpf("0.085")
    assert mp.mpf(data["standard_deviations"]) > mp.mpf("3.7e6")


def test_magnetization_axis_and_bulk_series() -> None:
    for shape in ((1, 1, 2), (1, 2, 2), (2, 2, 2)):
        joint = box_broken_bond_magnetization_poly(shape)
        marginal = [sum(row) for row in joint]
        while len(marginal) > 1 and marginal[-1] == 0:
            marginal.pop()
        assert marginal == box_broken_bond_poly(shape, plus_boundary=True)

    assert box_broken_bond_magnetization_poly((2, 2, 2)) == _brute_plus_joint((2, 2, 2))
    try:
        box_broken_bond_magnetization_poly((1, 1, 1))
    except ValueError as error:
        assert "transfer depth" in str(error)
    else:
        raise AssertionError("depth-one plus boundary must be rejected")

    series = low_temperature_magnetization_series(max_u_order=6)
    assert series == [
        Fraction(1),
        Fraction(0),
        Fraction(0),
        Fraction(-2),
        Fraction(0),
        Fraction(-12),
        Fraction(14),
    ]
    assert series[6] != Fraction(-18)


def test_degang_critical_condition_breaks_rotational_symmetry() -> None:
    data = degang_critical_data()
    roots = [mp.mpf(item["beta_c"]) for item in data["cyclic_permutations"]]
    assert max(roots) - min(roots) > mp.mpf("0.05")
    assert data["rotational_invariance_passed"] is False
    assert abs(mp.mpf(data["isotropic_beta_c"]) - mp.mpf("0.304688931718003")) < mp.mpf(
        "1e-15"
    )


def test_exact_ht_reference_and_general_falsifier() -> None:
    p0 = winding_projected_torus_polynomial(3)
    assert p0[:7] == (1, 0, 0, 0, 81, 0, 702)

    bulk = exact_ht_bulk_v_series(12)
    assert bulk[4] == Fraction(3)
    assert bulk[6] == Fraction(22)
    assert bulk[8] == Fraction(375, 2)
    assert bulk[10] == Fraction(1980)
    assert bulk[12] == Fraction(24044)

    def deliberately_bad_phi(k: mp.mpf) -> mp.mpf:
        return mp.log(2) + mp.mpf(3) * k**2 / 2 + mp.mpf(3) * k**4 / 4

    report = falsify_free_energy(deliberately_bad_phi, "controlled-K4-error")
    assert report["first_failure_order"] == 4
    assert abs(mp.mpf(report["candidate_coefficients"]["4"]) - mp.mpf("0.75")) < mp.mpf(
        "1e-20"
    )
    assert abs(mp.mpf(report["reference_coefficients"]["4"]) - mp.mpf("2.75")) < mp.mpf(
        "1e-60"
    )


def test_recorded_free_energy_discrepancy_is_reproducible() -> None:
    result_path = ROOT / "results" / "falsification" / "falsification.json"
    if not result_path.exists():
        run_analysis(write_result=True)
    recorded = json.loads(result_path.read_text())
    assert recorded["provenance"]["precision"] == PRECISION
    free = recorded["data"]["degang_2021"]["free_energy"]
    k = mp.mpf(free["K"])
    recomputed = claimed_degang_free_energy_m(k, k, k, 32)
    assert abs(recomputed - mp.mpf(free["phi_by_m"]["32"])) < mp.mpf("1e-65")
    assert mp.mpf(free["absolute_discrepancy_to_finite_lattice"]) > mp.mpf(
        free["finite_size_error_bound"]
    )
    assert "discrepancy_minus_finite_size_bound" in free
    assert "discrepancy_lower_bound_to_thermodynamic_limit" not in free
    assert recorded["data"]["general_falsifier"]["degang_2021"]["first_failure_order"] == 4
    assert all(check["passed"] for check in recorded["checks"])


def main() -> None:
    mp.mp.dps = PRECISION
    test_zhang_critical_condition()
    test_magnetization_axis_and_bulk_series()
    test_degang_critical_condition_breaks_rotational_symmetry()
    test_exact_ht_reference_and_general_falsifier()
    test_recorded_free_energy_discrepancy_is_reproducible()
    print("PASS: all falsification tests")


if __name__ == "__main__":
    main()
