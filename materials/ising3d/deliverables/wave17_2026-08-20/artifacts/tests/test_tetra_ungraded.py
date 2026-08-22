"""Standalone exact checks for the unrestricted auxiliary-R tetrahedron no-go."""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "e70_tetra_ungraded.py"
ARTIFACT = ROOT / "results" / "integrability" / "tetra_ungraded.json"


def _load_experiment():
    spec = importlib.util.spec_from_file_location("e70_tetra_ungraded_test", EXPERIMENT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {EXPERIMENT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_check(name: str, function) -> None:
    try:
        detail = function()
    except Exception as exc:
        print(f"{name}: FAIL ({exc})")
        raise
    print(f"{name}: PASS ({detail})")


def _basis_matrix_from_record(rows: list[dict], dimension: int) -> sp.Matrix:
    vectors = []
    for row in rows:
        vector = sp.zeros(dimension, 1)
        for entry in row["entries"]:
            vector[int(entry["coordinate_index"])] = sp.sympify(entry["coefficient"])
        vectors.append(vector)
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(dimension, 0)


def main() -> None:
    experiment = _load_experiment()
    cache: dict[str, object] = {}

    def check_full_system() -> str:
        matrix, coordinates = experiment.full_equation_matrix(Fraction(1, 4), scale=1)
        assert matrix.shape == (4096, 64)
        assert len(coordinates) == 64
        assert len(set(coordinates)) == 64
        assert set(coordinates) == set(product(product((0, 1), repeat=3), repeat=2))
        assert all(
            (sum(output_bits) - sum(input_bits)) % 2 == 0
            for output_bits, input_bits in coordinates[:32]
        )
        assert all(
            (sum(output_bits) - sum(input_bits)) % 2 == 1
            for output_bits, input_bits in coordinates[32:]
        )
        assert sum(
            any(matrix[row, column] != 0 for column in range(matrix.cols))
            for row in range(matrix.rows)
        ) == 4096
        direct = experiment.GRADED.direct_embedding_cross_check(
            matrix, coordinates, q=Fraction(1, 4), scale=1
        )
        assert direct["passed"]
        assert direct["compared_coefficients"] == 4096 * 64
        cache["fixed_matrix"] = matrix
        return "4096x64, all R entries, 262144 direct coefficients"

    def check_symbolic_minors() -> str:
        q = sp.Symbol("q")
        matrix, _ = experiment.full_equation_matrix(q, scale=1)
        even = matrix.extract(experiment.EVEN_MINOR_ROWS, range(32))
        odd_primary = matrix.extract(
            experiment.ODD_PRIMARY_MINOR_ROWS, range(32, 64)
        )
        odd_alternate = matrix.extract(
            experiment.ODD_ALTERNATE_MINOR_ROWS, range(32, 64)
        )
        assert matrix.extract(
            experiment.EVEN_MINOR_ROWS, range(32, 64)
        ).is_zero_matrix
        assert matrix.extract(
            experiment.ODD_PRIMARY_MINOR_ROWS, range(32)
        ).is_zero_matrix
        assert matrix.extract(
            experiment.ODD_ALTERNATE_MINOR_ROWS, range(32)
        ).is_zero_matrix

        even_det = sp.factor(even.det(method="domain-ge"))
        odd_primary_det = sp.factor(odd_primary.det(method="domain-ge"))
        odd_alternate_det = sp.factor(odd_alternate.det(method="domain-ge"))
        assert even_det == -2 * q**56 * (q - 1) ** 45 * (q + 1) ** 39
        assert odd_primary_det == (
            -q**53
            * (q - 1) ** 47
            * (q + 1) ** 38
            * (q**3 + 2 * q - 1)
        )
        assert odd_alternate_det == (
            -q**54 * (q - 1) ** 47 * (q + 1) ** 38 * (q**2 + 1)
        )

        full_primary_det = sp.factor(even_det * odd_primary_det)
        full_alternate_det = sp.factor(even_det * odd_alternate_det)
        common = q**109 * (q - 1) ** 92 * (q + 1) ** 77
        assert sp.cancel(
            full_primary_det - 2 * common * (q**3 + 2 * q - 1)
        ) == 0
        assert sp.cancel(
            full_alternate_det - 2 * common * q * (q**2 + 1)
        ) == 0
        determinant_gcd = sp.factor(
            sp.gcd(
                sp.Poly(full_primary_det, q, domain=sp.QQ),
                sp.Poly(full_alternate_det, q, domain=sp.QQ),
            ).as_expr()
        )
        assert determinant_gcd == common
        s = -(q**2 + q + 2) / 2
        t = (q**2 + q + 3) / 2
        assert sp.expand(s * (q**3 + 2 * q - 1) + t * q * (q**2 + 1)) == 1

        fixed = cache["fixed_matrix"]
        fixed_even = fixed.extract(experiment.EVEN_MINOR_ROWS, range(32)).det(
            method="domain-ge"
        )
        fixed_odd_primary = fixed.extract(
            experiment.ODD_PRIMARY_MINOR_ROWS, range(32, 64)
        ).det(method="domain-ge")
        assert fixed_even != 0 and fixed_odd_primary != 0
        cache["symbolic_determinants"] = (
            even_det,
            odd_primary_det,
            odd_alternate_det,
        )
        return "two reconstructed 64x64 block minors; gcd roots exactly 0,+1,-1"

    def check_algebraic_candidates() -> str:
        q = sp.Symbol("q")
        cubic = q**3 + 2 * q - 1
        quadratic = q**2 + 1
        inverse_alternate = sp.invert(q * quadratic, cubic, domain=sp.QQ)
        inverse_primary = sp.invert(cubic, quadratic, domain=sp.QQ)
        assert sp.rem(inverse_alternate * q * quadratic, cubic, domain=sp.QQ) == 1
        assert sp.rem(inverse_primary * cubic, quadratic, domain=sp.QQ) == 1
        analysis = experiment.algebraic_candidate_analysis()
        assert analysis["passed"]
        assert analysis["primary_only_candidate"]["exact_rank_at_every_root"] == 64
        assert analysis["primary_only_candidate"]["nullity_at_every_root"] == 0
        assert analysis["alternate_only_candidate"]["exact_rank_at_every_root"] == 64
        assert analysis["alternate_only_candidate"]["nullity_at_every_root"] == 0
        assert not analysis["primary_only_candidate"]["invertible_R_solution_exists"]
        assert not analysis["alternate_only_candidate"]["invertible_R_solution_exists"]
        return "cubic and q^2+1 roots have exact nullspace {0} by complementary minors"

    def check_exceptional_specializations() -> str:
        expected = {
            0: (14, 50, 26, 24),
            1: (26, 38, 19, 19),
            -1: (0, 64, 32, 32),
        }
        details = []
        for q_value, expected_row in expected.items():
            result = experiment.exceptional_specialization(q_value)
            assert (
                result["rank"],
                result["nullity"],
                result["sector_nullities"]["parity_preserving"],
                result["sector_nullities"]["parity_reversing"],
            ) == expected_row
            assert result["complete_exact_nullspace"]
            assert result["basis_residual_zero"]
            assert result["basis_independent"]
            assert result["rank_plus_nullity"] == 64

            matrix, coordinates = experiment.full_equation_matrix(q_value, scale=1)
            basis = _basis_matrix_from_record(result["nullspace_basis"], 64)
            assert basis.shape == (64, result["nullity"])
            assert matrix * basis == sp.zeros(4096, result["nullity"])
            assert basis.rank() == result["nullity"]
            _, pivots = matrix.rref()
            assert len(pivots) == result["rank"]

            identity = sp.zeros(64, 1)
            for index, (output_bits, input_bits) in enumerate(coordinates):
                if output_bits == input_bits:
                    identity[index] = 1
            assert matrix * identity == sp.zeros(4096, 1)
            identity_r = experiment._vector_to_r(identity, coordinates)
            assert identity_r == sp.eye(8)
            assert identity_r.det() == 1
            assert result["invertibility"] == {
                "invertible_solution_exists": True,
                "witness": "R=I_8",
                "witness_entries": result["invertibility"]["witness_entries"],
                "RLLL_residual_zero": True,
                "determinant": "1",
            }
            details.append(f"q={q_value}:({result['rank']},{result['nullity']})")
        return ", ".join(details)

    def check_artifact_envelope() -> str:
        artifact = json.loads(ARTIFACT.read_text())
        assert artifact["meta"]["provenance"] == "experiments/e70_tetra_ungraded.py"
        assert artifact["meta"]["predecessor"] == "experiments/e44_tetra_graded.py"
        assert set(artifact) == {"meta", "data", "checks"}
        assert all(check["passed"] for check in artifact["checks"])
        outcome = artifact["data"]["outcome"]
        assert outcome["class"] == "UNRESTRICTED_8X8_RLLL_NOGO_Q"
        assert artifact["data"]["generic_rank_certificate"]["generic_nullity"] == 0
        exceptional = artifact["data"]["true_exceptional_specializations"]
        assert {key: row["nullity"] for key, row in exceptional.items()} == {
            "0": 50,
            "1": 38,
            "-1": 64,
        }
        return "provenance/data/checks present and all stored checks pass"

    _run_check("full arbitrary-R component reconstruction", check_full_system)
    _run_check("complementary exact minor reconstruction", check_symbolic_minors)
    _run_check("candidate algebraic specialization clearance", check_algebraic_candidates)
    _run_check("true exceptional exact nullspaces", check_exceptional_specializations)
    _run_check("artifact envelope", check_artifact_envelope)
    print("PASS")


if __name__ == "__main__":
    main()
