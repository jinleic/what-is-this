"""Standalone acceptance checks for the exact Dolan--Grady deformation search."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from experiments.e17_dg_deformation import run_all


EXPECTED_COUNTS = {
    "grid_2x3": {
        "F1": (27, 9, 8, 18, 5),
        "F2": (81, 27, 26, 286, 78),
        "F3": (621, 186, 185, 1044, 282),
        "F4": (2, 1, 1, 22, 8),
    },
    "grid_3x3": {
        "F1": (49, 11, 10, 60, 10),
        "F2": (135, 27, 26, 884, 126),
        "F3": (1800, 318, 317, 5113, 711),
        "F4": (8, 2, 2, 80, 13),
    },
    "torus_3x3": {
        "F1": (90, 5, 4, 252, 5),
        "F2": (189, 9, 8, 2304, 40),
        "F3": (6102, 136, 135, 9909, 188),
        "F4": (36, 1, 1, 432, 10),
    },
}


def main() -> None:
    report = run_all(write_results=False, include_expansion_counts=False)
    data = report["data"]

    ring = data["controls"]["cycle_C6"]
    assert ring["dg2_floating_lambda"]["solution_exists"] is True
    assert ring["dg2_floating_lambda"]["lambda"] == "1"
    assert ring["residual_at_solution_terms"] == 0
    assert ring["tridiagonal_relation"]["solution"] == {
        "beta": "2",
        "gamma": "0",
        "rho": "16",
    }

    extra = data["controls"]["path_P6_plus_bond_1_4"]
    assert extra["dg2_floating_lambda"]["solution_exists"] is False
    assert extra["dg2_floating_lambda"]["quartic_obstruction_terms"] > 0
    assert extra["tridiagonal_relation"]["solution_exists"] is False
    print("1D controls: PASS")

    for graph_name, family_counts in EXPECTED_COUNTS.items():
        graph = data["graphs"][graph_name]
        for family_name, expected in family_counts.items():
            row = graph["families"][family_name]
            linear = row["b_linearized_quartic"]
            observed = (
                row["literal_word_count"],
                row["orbit_parameter_count"],
                linear["scale_normalized_unknowns"],
                linear["pauli_equations"],
                linear["symmetry_reduced_equations"],
            )
            assert observed == expected, (graph_name, family_name, observed, expected)
            print(f"{graph_name} {family_name} counts: PASS {observed}")

    torus = data["graphs"]["torus_3x3"]
    for family_name in ("F1", "F2", "F3", "F4", "F_ALL"):
        linear = torus["families"][family_name]["b_linearized_quartic"]
        assert linear["solution_exists"] is False
        assert linear["augmented_rank"] == linear["rank"] + 1

    f3_exact = torus["families"]["F3"]["a_exact_fixed_b"]
    assert f3_exact["commuting_solution_exists"] is False
    assert [row["lambda"] for row in f3_exact["noncommuting_solution_spaces"]] == ["1"]
    solution = f3_exact["noncommuting_solution_spaces"][0]
    assert solution["affine_dimension"] == 33
    assert solution["particular_nonzero_terms"] == [
        {
            "coefficient": "-1/2",
            "orbit_index": 69,
            "orbit_size": 18,
            "representative": "X0 Z1 Z2",
        }
    ]

    explicit = data["explicit_torus_F3_solution"]
    assert explicit["dg2_residual_terms"] == 0
    assert explicit["formula"] == (
        "sum_v X_v [1 - (Z_{v+x}Z_{v-x} + Z_{v+y}Z_{v-y})/2]"
    )
    assert explicit["dg1_floating_constant"]["solution_exists"] is False
    assert explicit["dg1_floating_constant"]["witness"]["coefficient_in_ad_A_cubed_B"] == "3"
    print("3x3 torus F3 exact cancellation and DG1 follow-up: PASS")

    nonlinear = data["nonlinear_exact"]
    assert nonlinear["torus_F1"]["only_decoupled_solutions"] is True
    assert nonlinear["torus_F1"]["groebner_basis_is_unit"] is False
    assert all(row["groebner_basis_is_unit"] for row in nonlinear["F4_b_deformation"].values())
    assert all(
        not row["nontrivial_connected_solution_exists"]
        for row in nonlinear["weighted_nearest_neighbour_ZZ"].values()
    )
    print("nonlinear classifications: PASS")

    assert all(check["passed"] for check in report["checks"]), report["checks"]
    print("PASS")


if __name__ == "__main__":
    main()
