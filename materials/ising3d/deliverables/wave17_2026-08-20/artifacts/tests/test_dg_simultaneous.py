"""Standalone acceptance checks for simultaneous DG/tridiagonal deformations."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from experiments.e36_dg_simultaneous import run_all


def main() -> None:
    try:
        report = run_all(write_results=False)
        data = report["data"]

        controls = data["controls"]
        cycle = controls["cycle_C6"]
        assert cycle["dg2_floating_lambda"]["pauli_equations"] == 12
        assert cycle["dg2_floating_lambda"]["lambda"] == "1"
        assert cycle["residual_at_solution_terms"] == 0
        assert cycle["tridiagonal_relation"]["solution"] == {
            "beta": "2",
            "gamma": "0",
            "rho": "16",
        }
        print("C6 zero-residual control: PASS")

        obstruction = controls["path_P6_plus_bond_1_4"]
        assert obstruction["dg2_floating_lambda"]["pauli_equations"] == 14
        assert obstruction["dg2_floating_lambda"]["solution_exists"] is False
        assert obstruction["dg2_floating_lambda"]["quartic_obstruction_terms"] == 2
        assert obstruction["tridiagonal_relation"]["solution_exists"] is False
        print("P6-plus-bond obstruction control: PASS")

        f3 = controls["torus_3x3_F3_DG2_cancellation"]
        assert f3["orbit_representative"] == "X0 Z1 Z2"
        assert f3["orbit_size"] == 18
        assert f3["lambda_direct"] == "16"
        assert f3["dg2_residual_terms"] == 0
        print("F3 DG2 cancellation control: PASS")

        witness = controls["torus_3x3_F3_DG1_witness"]
        assert witness == {
            "pauli_word": "Y0 Z1 X2 X3",
            "coefficient_in_ad_A_cubed_B": "3",
            "coefficient_in_ad_A_B": "0",
        }
        print("F3 DG1 (3,0) witness control: PASS")

        scope = data["scope"]
        assert scope["common_orbit_basis"] == ["O_X", "O_E", "O_C"]
        assert scope["common_orbit_basis_dimension"] == 3
        assert scope["raw_pair_coefficients"] == 6
        assert scope["projective_chart"]["free_deformation_coefficients"] == 4
        enumerations = scope["enumeration_by_size"]
        assert [row["group_order"] for row in enumerations] == [72, 128]
        assert [
            [orbit["literal_pauli_terms"] for orbit in row["orbits"]]
            for row in enumerations
        ] == [[9, 18, 18], [16, 32, 32]]
        print("minimal three-orbit ansatz enumeration: PASS")

        dg = data["dolan_grady"]
        assert dg["orbit_reduction"]["DG1"]["full_pauli_equations"] == 3384
        assert dg["orbit_reduction"]["DG1"]["symmetry_orbit_equations"] == 57
        assert dg["orbit_reduction"]["DG2"]["full_pauli_equations"] == 3384
        assert dg["orbit_reduction"]["DG2"]["symmetry_orbit_equations"] == 57
        assert dg["orbit_reduction"]["distinct_polynomials_both_relations"] == 60
        assert dg["linear_elimination"]["combined_distinct_polynomials"] == 160
        assert dg["groebner_after_linear_elimination"]["basis"] == [
            "(alpha*delta - 1)**2",
            "(alpha*delta - 1)*(chi*delta - eta)",
            "(chi*delta - eta)**2",
            "(alpha*delta - 1)*(alpha*eta - chi)",
            "(alpha*eta - chi)*(chi*delta - eta)",
            "(alpha*eta - chi)**2",
        ]
        assert all(
            basis == ["1"]
            for basis in dg["groebner_after_linear_elimination"][
                "nonproportional_saturation_charts"
            ].values()
        )
        assert dg["exact_variety"]["equivalent_operator_statement"] == "B_prime=delta*A_prime"
        print("simultaneous DG empty noncommuting variety: PASS")

        td = data["tridiagonal"]
        assert td["orbit_reduction"]["TD1"]["full_pauli_equations"] == 16344
        assert td["orbit_reduction"]["TD1"]["symmetry_orbit_equations"] == 263
        assert td["orbit_reduction"]["TD2"]["full_pauli_equations"] == 16344
        assert td["orbit_reduction"]["TD2"]["symmetry_orbit_equations"] == 263
        assert td["orbit_reduction"]["distinct_polynomials_both_relations"] == 266
        assert td["linear_elimination_then_groebner"] == {
            "eliminated_variable": "gamma_star=0 from the support-2 terminal row",
            "remaining_exact_constant": "48",
            "groebner_basis": ["1"],
        }
        assert td["exact_variety"]["equivalent_operator_statement"] == "B_prime=delta*A_prime"
        print("simultaneous tridiagonal empty noncommuting variety: PASS")

        size_rows = data["multiple_size_checks"]["sizes"]
        assert [row["layer"] for row in size_rows] == ["3x3_torus", "4x4_torus"]
        for row in size_rows:
            assert row["fixed_coefficients"] == {
                "alpha": "1/2",
                "chi": "-1/2",
                "delta": "2",
                "eta": "-1",
            }
            assert row["commutator_terms"] == 0
            assert row["DG1_residual_terms_at_mu_7"] == 0
            assert row["DG2_residual_terms_at_lambda_11"] == 0
            assert row["TD1_residual_terms_at_beta_3_gamma_5_rho_7"] == 0
            assert row["TD2_residual_terms_at_beta_3_gamma_star_11_rho_star_13"] == 0
            assert row["F3_DG2_residual_terms_at_lambda_16"] == 0
            assert row["F3_DG1_floating_constant_solution_exists"] is False
        assert data["decision"]["integrability"] is False
        assert all(check["passed"] for check in report["checks"])
        print("two-size fixed-coefficient checks: PASS")
        print("PASS")
    except Exception as error:
        print(f"failure: {type(error).__name__}: {error}")
        print("FAIL")
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
