"""Exact arbitrary-8x8 auxiliary-R tetrahedron/RLLL classification.

This removes the parity-preserving ansatz of ``e44_tetra_graded.py``.  All
64 entries of R are independent.  Two complementary exact 64x64 component
minors prove full column rank for every characteristic-zero parameter
q outside 0, +1, and -1; those three specializations are classified by exact
RREFs.

Run from the repository root with::

    .venv/bin/python experiments/e70_tetra_ungraded.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "tetra_ungraded.json"
GRADED_EXPERIMENT = ROOT / "experiments" / "e44_tetra_graded.py"
PHYSICAL_Q = Fraction(1, 4)
INTEGER_L_SCALE = 64

# The even-sector rows are the e44 full-rank certificate.  The first odd
# certificate has one additional cubic factor.  Replacing row 262 by row 274
# gives a complementary certificate whose extra factor is q^2+1.
EVEN_MINOR_ROWS = (
    0,
    3,
    5,
    6,
    9,
    10,
    12,
    15,
    17,
    18,
    33,
    34,
    36,
    39,
    40,
    43,
    45,
    46,
    65,
    66,
    68,
    71,
    192,
    195,
    257,
    263,
    264,
    267,
    269,
    270,
    272,
    275,
)
ODD_PRIMARY_MINOR_ROWS = (
    1,
    2,
    4,
    7,
    8,
    11,
    13,
    14,
    16,
    19,
    32,
    35,
    37,
    38,
    41,
    42,
    44,
    47,
    64,
    67,
    69,
    70,
    193,
    194,
    256,
    259,
    261,
    262,
    265,
    266,
    268,
    271,
)
ODD_ALTERNATE_MINOR_ROWS = (
    1,
    2,
    4,
    7,
    8,
    11,
    13,
    14,
    16,
    19,
    32,
    35,
    37,
    38,
    41,
    42,
    44,
    47,
    64,
    67,
    69,
    70,
    193,
    194,
    256,
    259,
    261,
    265,
    266,
    268,
    271,
    274,
)
TRUE_EXCEPTIONAL_Q = (0, 1, -1)


def _load_graded_experiment():
    spec = importlib.util.spec_from_file_location(
        "e44_tetra_graded_for_e70", GRADED_EXPERIMENT
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {GRADED_EXPERIMENT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GRADED = _load_graded_experiment()


def odd_coordinates() -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """Return the 32 parity-reversing entries of an auxiliary 8x8 R."""

    coordinates = []
    for input_bits in product((0, 1), repeat=3):
        for output_bits in product((0, 1), repeat=3):
            if (sum(output_bits) - sum(input_bits)) % 2:
                coordinates.append((output_bits, input_bits))
    return tuple(coordinates)


def ungraded_coordinates() -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """Return all 64 R entries, block-ordered only to expose the certificate."""

    coordinates = tuple(GRADED.graded_coordinates()) + odd_coordinates()
    if len(coordinates) != 64 or len(set(coordinates)) != 64:
        raise ArithmeticError("the ungraded coordinate list is not all 64 R entries")
    return coordinates


def coordinate_name(coordinate: tuple[tuple[int, ...], tuple[int, ...]]) -> str:
    return GRADED.coordinate_name(coordinate)


def coordinate_symbols(
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]] | None = None,
) -> tuple[sp.Symbol, ...]:
    if coordinates is None:
        coordinates = ungraded_coordinates()
    return tuple(sp.Symbol(coordinate_name(coordinate)) for coordinate in coordinates)


def full_equation_matrix(
    q: int | Fraction | sp.Expr = PHYSICAL_Q, scale: int = INTEGER_L_SCALE
) -> tuple[sp.Matrix, tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]]:
    """Build all 4096 equations B(q) vec(R)=0 for arbitrary 8x8 R."""

    coordinates = ungraded_coordinates()
    forward, reverse = GRADED._lll_products(q, scale=scale)
    rows = [
        GRADED._component_row(output_state, input_state, forward, reverse, coordinates)
        for output_state in range(64)
        for input_state in range(64)
    ]
    return sp.Matrix(rows), coordinates


def _row_parity(row_index: int) -> int:
    output_state, input_state = divmod(int(row_index), 64)
    return (output_state.bit_count() - input_state.bit_count()) % 2


def _sector_minor(
    q: int | Fraction | sp.Expr,
    rows: Sequence[int],
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
    scale: int = 1,
) -> sp.Matrix:
    forward, reverse = GRADED._lll_products(q, scale=scale)
    return sp.Matrix(
        [
            GRADED._component_row(
                int(row) // 64,
                int(row) % 64,
                forward,
                reverse,
                coordinates,
            )
            for row in rows
        ]
    )


def _factor_list(expression: sp.Expr, variable: sp.Symbol) -> dict:
    coefficient, factors = sp.factor_list(expression, variable)
    return {
        "coefficient": str(coefficient),
        "factors": [[str(factor), int(exponent)] for factor, exponent in factors],
    }


def _integer_factorization(value: int) -> dict[str, int]:
    return {
        str(prime): int(exponent)
        for prime, exponent in sp.factorint(abs(int(value))).items()
    }


def generic_rank_certificate() -> dict:
    """Construct and factor two complementary exact full-column-rank minors."""

    started = time.perf_counter()
    q = sp.Symbol("q")
    even_coordinates = ungraded_coordinates()[:32]
    odd = ungraded_coordinates()[32:]

    even_minor = _sector_minor(q, EVEN_MINOR_ROWS, even_coordinates)
    odd_primary_minor = _sector_minor(q, ODD_PRIMARY_MINOR_ROWS, odd)
    odd_alternate_minor = _sector_minor(q, ODD_ALTERNATE_MINOR_ROWS, odd)

    even_determinant = sp.factor(even_minor.det(method="domain-ge"))
    odd_primary_determinant = sp.factor(
        odd_primary_minor.det(method="domain-ge")
    )
    odd_alternate_determinant = sp.factor(
        odd_alternate_minor.det(method="domain-ge")
    )
    full_primary_determinant = sp.factor(
        even_determinant * odd_primary_determinant
    )
    full_alternate_determinant = sp.factor(
        even_determinant * odd_alternate_determinant
    )

    expected_even = -2 * q**56 * (q - 1) ** 45 * (q + 1) ** 39
    expected_odd_primary = (
        -q**53 * (q - 1) ** 47 * (q + 1) ** 38 * (q**3 + 2 * q - 1)
    )
    expected_odd_alternate = (
        -q**54 * (q - 1) ** 47 * (q + 1) ** 38 * (q**2 + 1)
    )
    common_root_factor = q**109 * (q - 1) ** 92 * (q + 1) ** 77
    primary_extra = q**3 + 2 * q - 1
    alternate_extra = q * (q**2 + 1)
    bezout_primary = -(q**2 + q + 2) / 2
    bezout_alternate = (q**2 + q + 3) / 2
    bezout_identity = sp.expand(
        bezout_primary * primary_extra + bezout_alternate * alternate_extra
    )
    determinant_gcd = sp.factor(
        sp.gcd(
            sp.Poly(full_primary_determinant, q, domain=sp.QQ),
            sp.Poly(full_alternate_determinant, q, domain=sp.QQ),
        ).as_expr()
    )

    fixed_substitution = sp.Rational(PHYSICAL_Q.numerator, PHYSICAL_Q.denominator)
    scaled_even = sp.cancel(
        even_determinant.subs(q, fixed_substitution) * INTEGER_L_SCALE ** (3 * 32)
    )
    scaled_odd_primary = sp.cancel(
        odd_primary_determinant.subs(q, fixed_substitution)
        * INTEGER_L_SCALE ** (3 * 32)
    )
    scaled_odd_alternate = sp.cancel(
        odd_alternate_determinant.subs(q, fixed_substitution)
        * INTEGER_L_SCALE ** (3 * 32)
    )
    scaled_full_primary = sp.cancel(scaled_even * scaled_odd_primary)
    scaled_full_alternate = sp.cancel(scaled_even * scaled_odd_alternate)

    expected_match = (
        sp.cancel(even_determinant - expected_even) == 0
        and sp.cancel(odd_primary_determinant - expected_odd_primary) == 0
        and sp.cancel(odd_alternate_determinant - expected_odd_alternate) == 0
    )
    certificate_valid = (
        expected_match
        and bezout_identity == 1
        and sp.cancel(full_primary_determinant - 2 * common_root_factor * primary_extra)
        == 0
        and sp.cancel(
            full_alternate_determinant
            - 2 * common_root_factor * alternate_extra
        )
        == 0
        and sp.cancel(determinant_gcd - common_root_factor) == 0
    )
    if not certificate_valid:
        raise ArithmeticError("the complementary symbolic-minor certificate failed")

    return {
        "variable": "q=w^2=tanh(K)",
        "matrix_shape": [4096, 64],
        "exact_upper_bound": 64,
        "certified_rank_outside_exceptional_set": 64,
        "generic_nullity": 0,
        "true_exceptional_q": [str(value) for value in TRUE_EXCEPTIONAL_Q],
        "even_sector": {
            "rows": list(EVEN_MINOR_ROWS),
            "shape": [32, 32],
            "determinant": str(even_determinant),
            "factorization": _factor_list(even_determinant, q),
        },
        "odd_sector_primary": {
            "rows": list(ODD_PRIMARY_MINOR_ROWS),
            "shape": [32, 32],
            "determinant": str(odd_primary_determinant),
            "factorization": _factor_list(odd_primary_determinant, q),
        },
        "odd_sector_alternate": {
            "rows": list(ODD_ALTERNATE_MINOR_ROWS),
            "shape": [32, 32],
            "determinant": str(odd_alternate_determinant),
            "factorization": _factor_list(odd_alternate_determinant, q),
        },
        "full_primary_minor": {
            "row_order": list(EVEN_MINOR_ROWS + ODD_PRIMARY_MINOR_ROWS),
            "column_order": "32 parity-preserving coordinates, then 32 parity-reversing coordinates",
            "shape": [64, 64],
            "determinant": str(full_primary_determinant),
            "factorization": _factor_list(full_primary_determinant, q),
        },
        "full_alternate_minor": {
            "row_order": list(EVEN_MINOR_ROWS + ODD_ALTERNATE_MINOR_ROWS),
            "column_order": "32 parity-preserving coordinates, then 32 parity-reversing coordinates",
            "shape": [64, 64],
            "determinant": str(full_alternate_determinant),
            "factorization": _factor_list(full_alternate_determinant, q),
        },
        "coverage": {
            "common_root_factor": str(common_root_factor),
            "monic_polynomial_gcd_of_full_minors": str(determinant_gcd),
            "primary_extra_factor": str(primary_extra),
            "alternate_extra_factor": str(alternate_extra),
            "bezout_primary_coefficient": str(bezout_primary),
            "bezout_alternate_coefficient": str(bezout_alternate),
            "bezout_identity": str(bezout_identity),
            "conclusion": (
                "At q not in {0,+1,-1}, the common root factor is nonzero and the "
                "Bezout identity prevents both 64x64 minors from vanishing."
            ),
        },
        "fixed_q_scaled_cross_check": {
            "q": "1/4",
            "uniform_L_scale": INTEGER_L_SCALE,
            "minor_scaling_rule": "an n-column minor is multiplied by 64^(3*n)",
            "even_determinant": str(scaled_even),
            "even_factorization": _integer_factorization(int(scaled_even)),
            "odd_primary_determinant": str(scaled_odd_primary),
            "odd_primary_factorization": _integer_factorization(
                int(scaled_odd_primary)
            ),
            "odd_alternate_determinant": str(scaled_odd_alternate),
            "odd_alternate_factorization": _integer_factorization(
                int(scaled_odd_alternate)
            ),
            "full_primary_determinant": str(scaled_full_primary),
            "full_primary_factorization": _integer_factorization(
                int(scaled_full_primary)
            ),
            "full_alternate_determinant": str(scaled_full_alternate),
            "full_alternate_factorization": _integer_factorization(
                int(scaled_full_alternate)
            ),
        },
        "certificate_valid": certificate_valid,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def algebraic_candidate_analysis() -> dict:
    """Clear the non-common cubic and quadratic roots by the other minor."""

    q = sp.Symbol("q")
    cubic = q**3 + 2 * q - 1
    quadratic = q**2 + 1
    alternate_extra = q * quadratic
    cubic_inverse_of_alternate = sp.invert(alternate_extra, cubic, domain=sp.QQ)
    quadratic_inverse_of_primary = sp.invert(cubic, quadratic, domain=sp.QQ)
    cubic_check = sp.rem(
        sp.expand(cubic_inverse_of_alternate * alternate_extra),
        cubic,
        domain=sp.QQ,
    )
    quadratic_check = sp.rem(
        sp.expand(quadratic_inverse_of_primary * cubic),
        quadratic,
        domain=sp.QQ,
    )
    exact_gcd = sp.gcd(sp.Poly(cubic, q), sp.Poly(alternate_extra, q)).as_expr()
    passed = cubic_check == 1 and quadratic_check == 1 and exact_gcd == 1
    if not passed:
        raise ArithmeticError("failed to clear a candidate algebraic minor root")
    return {
        "candidate_sets_are_not_true_exceptions": True,
        "exact_gcd": str(exact_gcd),
        "primary_only_candidate": {
            "minimal_polynomial": str(cubic),
            "irreducible_over_Q": bool(sp.Poly(cubic, q, domain=sp.QQ).is_irreducible),
            "number_of_conjugate_roots": 3,
            "cleared_by": "full_alternate_minor",
            "inverse_of_alternate_extra_mod_minimal_polynomial": str(
                cubic_inverse_of_alternate
            ),
            "inverse_check_remainder": str(cubic_check),
            "exact_rank_at_every_root": 64,
            "nullity_at_every_root": 0,
            "nullspace_basis_at_every_root": [],
            "invertible_R_solution_exists": False,
        },
        "alternate_only_candidate": {
            "minimal_polynomial": str(quadratic),
            "irreducible_over_Q": bool(
                sp.Poly(quadratic, q, domain=sp.QQ).is_irreducible
            ),
            "number_of_conjugate_roots": 2,
            "cleared_by": "full_primary_minor",
            "inverse_of_primary_extra_mod_minimal_polynomial": str(
                quadratic_inverse_of_primary
            ),
            "inverse_check_remainder": str(quadratic_check),
            "exact_rank_at_every_root": 64,
            "nullity_at_every_root": 0,
            "nullspace_basis_at_every_root": [],
            "invertible_R_solution_exists": False,
        },
        "passed": passed,
    }


def _nullspace_from_rref(
    reduced: sp.Matrix, pivot_columns: Sequence[int], column_count: int
) -> tuple[sp.Matrix, ...]:
    pivot_columns = tuple(int(column) for column in pivot_columns)
    free_columns = tuple(
        column for column in range(column_count) if column not in set(pivot_columns)
    )
    basis = []
    for free_column in free_columns:
        vector = sp.zeros(column_count, 1)
        vector[free_column] = 1
        for row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = -reduced[row, free_column]
        basis.append(vector)
    return tuple(basis)


def _vector_to_r(
    vector: sp.Matrix,
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> sp.Matrix:
    matrix = sp.zeros(8, 8)
    for index, (output_bits, input_bits) in enumerate(coordinates):
        matrix[GRADED._state(output_bits), GRADED._state(input_bits)] = vector[index]
    return matrix


def _sparse_basis_record(
    basis: Sequence[sp.Matrix],
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[dict]:
    records = []
    for basis_index, vector in enumerate(basis):
        entries = []
        for coordinate_index, coefficient in enumerate(vector):
            if coefficient != 0:
                output_bits, input_bits = coordinates[coordinate_index]
                entries.append(
                    {
                        "coordinate_index": coordinate_index,
                        "coordinate": coordinate_name(coordinates[coordinate_index]),
                        "output": list(output_bits),
                        "input": list(input_bits),
                        "coefficient": str(coefficient),
                    }
                )
        records.append({"basis_index": basis_index, "entries": entries})
    return records


def exceptional_specialization(q_value: int) -> dict:
    """Compute the full exact RREF/nullspace and an invertible witness."""

    if q_value not in TRUE_EXCEPTIONAL_Q:
        raise ValueError(f"q={q_value} is not a true exceptional specialization")
    started = time.perf_counter()
    matrix, coordinates = full_equation_matrix(q_value, scale=1)
    reduced, pivot_columns = matrix.rref()
    basis = _nullspace_from_rref(reduced, pivot_columns, matrix.cols)
    basis_matrix = sp.Matrix.hstack(*basis) if basis else sp.zeros(matrix.cols, 0)
    residual_zero = matrix * basis_matrix == sp.zeros(matrix.rows, len(basis))
    basis_independent = basis_matrix.rank() == len(basis)

    identity_vector = sp.zeros(64, 1)
    for index, (output_bits, input_bits) in enumerate(coordinates):
        if output_bits == input_bits:
            identity_vector[index] = 1
    identity_r = _vector_to_r(identity_vector, coordinates)
    identity_residual_zero = matrix * identity_vector == sp.zeros(matrix.rows, 1)
    identity_determinant = sp.factor(identity_r.det())

    even_pivots = sum(int(column) < 32 for column in pivot_columns)
    odd_pivots = len(pivot_columns) - even_pivots
    completeness = (
        residual_zero
        and basis_independent
        and len(pivot_columns) + len(basis) == 64
    )
    if not completeness or not identity_residual_zero or identity_determinant == 0:
        raise ArithmeticError(f"exceptional classification failed at q={q_value}")

    return {
        "q": str(q_value),
        "field": "QQ",
        "matrix_shape": [matrix.rows, matrix.cols],
        "nonzero_component_equations": sum(
            any(matrix[row, column] != 0 for column in range(matrix.cols))
            for row in range(matrix.rows)
        ),
        "rank": len(pivot_columns),
        "nullity": len(basis),
        "sector_ranks": {
            "parity_preserving": even_pivots,
            "parity_reversing": odd_pivots,
        },
        "sector_nullities": {
            "parity_preserving": 32 - even_pivots,
            "parity_reversing": 32 - odd_pivots,
        },
        "pivot_columns": list(pivot_columns),
        "free_columns": [
            column for column in range(64) if column not in set(pivot_columns)
        ],
        "nullspace_basis": _sparse_basis_record(basis, coordinates),
        "basis_residual_zero": residual_zero,
        "basis_independent": basis_independent,
        "rank_plus_nullity": len(pivot_columns) + len(basis),
        "invertibility": {
            "invertible_solution_exists": True,
            "witness": "R=I_8",
            "witness_entries": [
                {
                    "coordinate_index": index,
                    "coordinate": coordinate_name(coordinates[index]),
                    "coefficient": "1",
                }
                for index, (output_bits, input_bits) in enumerate(coordinates)
                if output_bits == input_bits
            ],
            "RLLL_residual_zero": identity_residual_zero,
            "determinant": str(identity_determinant),
        },
        "complete_exact_nullspace": completeness,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def _coordinate_table(
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]]
) -> list[dict]:
    return [
        {
            "index": index,
            "name": coordinate_name(coordinate),
            "output": list(coordinate[0]),
            "input": list(coordinate[1]),
            "parity_change": (sum(coordinate[0]) - sum(coordinate[1])) % 2,
        }
        for index, coordinate in enumerate(coordinates)
    ]


def _fixed_system_analysis() -> dict:
    started = time.perf_counter()
    matrix, coordinates = full_equation_matrix(PHYSICAL_Q, scale=INTEGER_L_SCALE)
    direct = GRADED.direct_embedding_cross_check(
        matrix, coordinates, q=PHYSICAL_Q, scale=INTEGER_L_SCALE
    )
    even_columns = range(32)
    odd_columns = range(32, 64)
    even_rows = [row for row in range(4096) if _row_parity(row) == 0]
    odd_rows = [row for row in range(4096) if _row_parity(row) == 1]
    cross_even_rows_odd_columns = sum(
        matrix[row, column] != 0 for row in even_rows for column in odd_columns
    )
    cross_odd_rows_even_columns = sum(
        matrix[row, column] != 0 for row in odd_rows for column in even_columns
    )

    even_minor = matrix.extract(EVEN_MINOR_ROWS, range(32))
    odd_primary_minor = matrix.extract(ODD_PRIMARY_MINOR_ROWS, range(32, 64))
    odd_alternate_minor = matrix.extract(ODD_ALTERNATE_MINOR_ROWS, range(32, 64))
    fixed_even = int(even_minor.det(method="domain-ge"))
    fixed_odd_primary = int(odd_primary_minor.det(method="domain-ge"))
    fixed_odd_alternate = int(odd_alternate_minor.det(method="domain-ge"))

    nonzero_rows = sum(
        any(matrix[row, column] != 0 for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    return {
        "q": "1/4",
        "uniform_L_scale": INTEGER_L_SCALE,
        "shape": [matrix.rows, matrix.cols],
        "coordinates": _coordinate_table(coordinates),
        "component_equations": matrix.rows,
        "unknowns": matrix.cols,
        "nonzero_component_equations": nonzero_rows,
        "parity_block_decomposition": {
            "even_component_rows": len(even_rows),
            "odd_component_rows": len(odd_rows),
            "parity_preserving_columns": 32,
            "parity_reversing_columns": 32,
            "even_row_to_odd_column_nonzeros": cross_even_rows_odd_columns,
            "odd_row_to_even_column_nonzeros": cross_odd_rows_even_columns,
            "interpretation": (
                "This is a direct sum decomposition of the full 64-variable system, "
                "not a restriction on R."
            ),
        },
        "direct_operator_cross_check": direct,
        "fixed_minor_determinants": {
            "even": str(fixed_even),
            "odd_primary": str(fixed_odd_primary),
            "odd_alternate": str(fixed_odd_alternate),
            "full_primary_product": str(fixed_even * fixed_odd_primary),
            "full_alternate_product": str(fixed_even * fixed_odd_alternate),
        },
        "fixed_minor_factorizations": {
            "even": _integer_factorization(fixed_even),
            "odd_primary": _integer_factorization(fixed_odd_primary),
            "odd_alternate": _integer_factorization(fixed_odd_alternate),
            "full_primary_product": _integer_factorization(
                fixed_even * fixed_odd_primary
            ),
            "full_alternate_product": _integer_factorization(
                fixed_even * fixed_odd_alternate
            ),
        },
        "exact_rank": 64,
        "nullity": 0,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def build_artifact() -> dict:
    controls_started = time.perf_counter()
    controls = GRADED.physical_tensor_controls()
    controls_elapsed = time.perf_counter() - controls_started
    fixed = _fixed_system_analysis()
    generic = generic_rank_certificate()
    algebraic = algebraic_candidate_analysis()
    exceptional = {
        str(q_value): exceptional_specialization(q_value)
        for q_value in TRUE_EXCEPTIONAL_Q
    }

    fixed_symbolic = generic["fixed_q_scaled_cross_check"]
    fixed_direct = fixed["fixed_minor_determinants"]
    fixed_cross_match = (
        fixed_direct["even"] == fixed_symbolic["even_determinant"]
        and fixed_direct["odd_primary"]
        == fixed_symbolic["odd_primary_determinant"]
        and fixed_direct["odd_alternate"]
        == fixed_symbolic["odd_alternate_determinant"]
        and fixed_direct["full_primary_product"]
        == fixed_symbolic["full_primary_determinant"]
        and fixed_direct["full_alternate_product"]
        == fixed_symbolic["full_alternate_determinant"]
    )
    expected_exceptional = {
        "0": (14, 50, 26, 24),
        "1": (26, 38, 19, 19),
        "-1": (0, 64, 32, 32),
    }
    exceptional_match = all(
        (
            row["rank"],
            row["nullity"],
            row["sector_nullities"]["parity_preserving"],
            row["sector_nullities"]["parity_reversing"],
        )
        == expected_exceptional[key]
        and row["complete_exact_nullspace"]
        and row["invertibility"]["invertible_solution_exists"]
        and row["invertibility"]["RLLL_residual_zero"]
        and row["invertibility"]["determinant"] == "1"
        for key, row in exceptional.items()
    )

    checks = [
        {
            "name": "physical Ising tensor provenance controls pass",
            "passed": (
                controls["coefficientwise_matches"]
                and controls["scaled_partition_integer_matches"]
                and not controls["local_entry_reconstruction_mismatches"]
            ),
            "detail": str(controls["actual_scaled_partition_integers"]),
        },
        {
            "name": "full arbitrary-8x8 component system matches direct embedding",
            "passed": (
                fixed["shape"] == [4096, 64]
                and len(fixed["coordinates"]) == 64
                and fixed["direct_operator_cross_check"]["passed"]
                and fixed["nonzero_component_equations"] == 4096
            ),
            "detail": (
                f"{fixed['component_equations']} equations, {fixed['unknowns']} independent R entries, "
                f"{fixed['direct_operator_cross_check']['compared_coefficients']} coefficients compared"
            ),
        },
        {
            "name": "parity organization retains both 32-variable sectors without leakage",
            "passed": (
                fixed["parity_block_decomposition"][
                    "even_row_to_odd_column_nonzeros"
                ]
                == 0
                and fixed["parity_block_decomposition"][
                    "odd_row_to_even_column_nonzeros"
                ]
                == 0
                and fixed["parity_block_decomposition"][
                    "parity_preserving_columns"
                ]
                == 32
                and fixed["parity_block_decomposition"][
                    "parity_reversing_columns"
                ]
                == 32
            ),
            "detail": "all 64 R coordinates are present as a direct sum of two equation sectors",
        },
        {
            "name": "two symbolic 64x64 minors certify rank 64 away from 0,+1,-1",
            "passed": generic["certificate_valid"] and fixed_cross_match,
            "detail": (
                f"gcd={generic['coverage']['monic_polynomial_gcd_of_full_minors']}; "
                f"Bezout={generic['coverage']['bezout_identity']}"
            ),
        },
        {
            "name": "all additional algebraic roots from individual minors are cleared",
            "passed": algebraic["passed"],
            "detail": (
                "cubic roots have alternate-minor rank 64; q^2+1 roots have primary-minor rank 64"
            ),
        },
        {
            "name": "true exceptional specializations have complete nullspaces and invertible witnesses",
            "passed": exceptional_match,
            "detail": "; ".join(
                f"q={key}: rank {row['rank']}, nullity {row['nullity']}, det(I_8)={row['invertibility']['determinant']}"
                for key, row in exceptional.items()
            ),
        },
    ]

    return {
        "meta": {
            "provenance": "experiments/e70_tetra_ungraded.py",
            "predecessor": "experiments/e44_tetra_graded.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "arithmetic": (
                "exact Python integers/Fraction and SymPy QQ polynomial/RREF algebra; "
                "floating point is used only for elapsed wall times"
            ),
            "sympy_version": sp.__version__,
            "physical_controls_elapsed_seconds": f"{controls_elapsed:.6f}",
        },
        "data": {
            "ansatz": {
                "equation": (
                    "R_123 L_145(q) L_246(q) L_356(q) = "
                    "L_356(q) L_246(q) L_145(q) R_123"
                ),
                "operator_convention": "L[output,input]=T[input,output]/2",
                "global_operator_shape": [64, 64],
                "auxiliary_R_shape": [8, 8],
                "R_entries": 64,
                "independent_unknowns": 64,
                "grading_assumption_on_R": "none",
                "other_R_symmetries": "none",
                "identical_local_factors": True,
            },
            "physical_tensor_controls": controls,
            "full_system_at_q_1_4": fixed,
            "generic_rank_certificate": generic,
            "algebraic_candidate_specializations": algebraic,
            "true_exceptional_specializations": exceptional,
            "outcome": {
                "class": "UNRESTRICTED_8X8_RLLL_NOGO_Q",
                "theorem": (
                    "Over every characteristic-zero field and for every q not in {0,+1,-1}, "
                    "the only arbitrary 8x8 solution of the displayed identical-L Ising RLLL "
                    "equation is R=0. Thus no nonzero or invertible R exists there."
                ),
                "physical_corollary": (
                    "In particular, no nonzero arbitrary 8x8 auxiliary R exists for finite "
                    "ferromagnetic weights 0<q<1."
                ),
                "exceptional_statement": (
                    "At q=0,+1,-1 the exact nullities are 50,38,64 respectively, and R=I_8 "
                    "is an invertible solution at each specialization."
                ),
                "scope_limit": (
                    "Higher auxiliary dimension, nonidentical or spectral L factors, IRF/dynamical "
                    "relations, and dimension-changing projections are not classified."
                ),
            },
            "status": (
                "[THEOREM] arbitrary-8x8 identical-L RLLL no-go for every characteristic-zero "
                "q outside {0,+1,-1}; [EXACT EXCEPTIONS] complete nullspaces at 0,+1,-1."
            ),
        },
        "checks": checks,
    }


def main() -> None:
    artifact = build_artifact()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    for check in artifact["checks"]:
        print(check["name"], "PASS" if check["passed"] else "FAIL", check["detail"])
    if not all(check["passed"] for check in artifact["checks"]):
        raise AssertionError("one or more ungraded tetrahedron checks failed")
    print("written", OUTPUT.relative_to(ROOT))
    print("PASS")


if __name__ == "__main__":
    main()
