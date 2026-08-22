"""Standalone exact checks for the graded auxiliary-R tetrahedron no-go."""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "e44_tetra_graded.py"


def _load_experiment():
    spec = importlib.util.spec_from_file_location("e44_tetra_graded_test", EXPERIMENT)
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


def main() -> None:
    experiment = _load_experiment()

    def check_partitions() -> str:
        controls = experiment.physical_tensor_controls()
        expected = {
            "2x2x2": 36450,
            "2x2x3": 16394562,
            "2x3x3": 246853161090,
        }
        assert controls["actual_scaled_partition_integers"] == expected
        assert controls["scaled_partition_integer_matches"]
        assert controls["coefficientwise_matches"]
        assert controls["local_entry_reconstruction_mismatches"] == []
        return str(expected)

    def check_grading() -> str:
        violations = []
        support = 0
        for input_bits in product((0, 1), repeat=3):
            for output_bits in product((0, 1), repeat=3):
                value = experiment.physical_l_entry(output_bits, input_bits, Fraction(1, 4))
                if value != 0:
                    support += 1
                    if (sum(input_bits) - sum(output_bits)) % 2:
                        violations.append((output_bits, input_bits, value))
        assert support == 32
        assert violations == []
        analysis = experiment.physical_grade_analysis(Fraction(1, 4))
        assert analysis["grade_preserving"]
        assert analysis["grade_block_ranks"] == {"even": 1, "odd": 1}
        return "32/64 supported entries, zero grade violations"

    cache: dict[str, object] = {}

    def check_rank_and_symbolic_minor() -> str:
        matrix, coordinates = experiment.graded_equation_matrix(
            Fraction(1, 4), scale=experiment.INTEGER_L_SCALE
        )
        assert matrix.shape == (4096, 32)
        assert len(coordinates) == 32
        assert sum(
            any(matrix[row, column] != 0 for column in range(matrix.cols))
            for row in range(matrix.rows)
        ) == 2048
        certificate = experiment.exact_rank_certificate(matrix)
        assert certificate["rank"] == 32
        assert certificate["determinant"] != 0
        assert certificate["determinant_factorization"] == {"2": 297, "3": 45, "5": 39}
        symbolic = experiment.symbolic_minor_certificate(
            certificate["pivot_rows"], certificate["determinant"]
        )
        q = sp.Symbol("q")
        assert sp.expand(
            sp.sympify(symbolic["determinant"])
            + 2 * q**56 * (q - 1) ** 45 * (q + 1) ** 39
        ) == 0
        assert symbolic["certificate_valid"]
        assert symbolic["fixed_q_scaled_determinant_cross_check"]
        cache["certificate"] = certificate
        cache["coordinates"] = coordinates
        return "rank 32; det=-2*q^56*(q-1)^45*(q+1)^39"

    def check_normalized_groebner() -> str:
        certificate = cache["certificate"]
        coordinates = cache["coordinates"]
        variables = experiment.coordinate_symbols(coordinates)
        minor = certificate["minor"]
        equations = [
            sp.Add(*(minor[row, column] * variables[column] for column in range(32)))
            for row in range(32)
        ]
        result = experiment.run_groebner(
            equations + [variables[0] - 1],
            variables,
            timeout_seconds=120,
        )
        assert result["status"] == "completed"
        assert result["basis"] == ["1"]
        assert result["unit_ideal"]
        return f"raw normalization patch basis {result['basis']} over {result['field']}"

    _run_check("physical-tensor partition controls", check_partitions)
    _run_check("physical-tensor Z2 grading", check_grading)
    _run_check("exact graded-system rank certificate", check_rank_and_symbolic_minor)
    _run_check("representative saturated Groebner certificate", check_normalized_groebner)
    print("PASS")


if __name__ == "__main__":
    main()
