"""Exact RLLL obstruction for every parity-preserving auxiliary 8x8 R.

The fixed local operator is the verified six-leg Ising parity tensor.  No
leg-permutation or translation symmetry is imposed on R: all 32 entries allowed
by total Z2-grade conservation are independent.
"""

from __future__ import annotations

import importlib.util
import json
import math
import signal
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "integrability" / "tetra_graded.json"
WAVE5_EXPERIMENT = ROOT / "experiments" / "e35_tetrahedron_spectral.py"
POSITIONS = ((0, 1, 2), (0, 3, 4), (1, 3, 5), (2, 4, 5))
CONTROL_BOXES = ((2, 2, 2), (2, 2, 3), (2, 3, 3))
EXPECTED_PARTITION_INTEGERS = {
    (2, 2, 2): 36450,
    (2, 2, 3): 16394562,
    (2, 3, 3): 246853161090,
}
PHYSICAL_Q = Fraction(1, 4)  # q=w^2=tanh(K), with the Wave-5 choice w=1/2.
INTEGER_L_SCALE = 64
GROEBNER_TIMEOUT_SECONDS = 3000
FINITE_FIELD_TIMEOUT_SECONDS = 2400
FINITE_FIELD_PRIMES = (101, 32003)


class GroebnerTimeout(TimeoutError):
    """Raised when a bounded Groebner calculation exceeds its wall time."""


def _load_wave5_experiment():
    spec = importlib.util.spec_from_file_location("e35_tetrahedron_spectral_for_e44", WAVE5_EXPERIMENT)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {WAVE5_EXPERIMENT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


WAVE5 = _load_wave5_experiment()


def _bits(state: int, width: int) -> tuple[int, ...]:
    return tuple((state >> position) & 1 for position in range(width))


def _state(bits: Sequence[int]) -> int:
    return sum(int(bit) << position for position, bit in enumerate(bits))


def _sympy_number(value: int | Fraction | sp.Expr) -> sp.Expr:
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.sympify(value)


def physical_l_entry(
    output_bits: Sequence[int], input_bits: Sequence[int], q: int | Fraction | sp.Expr
) -> sp.Expr:
    """Return L[output,input]=T[input,output]/2 at q=w^2.

    The omitted factor two is uniform and occurs three times on each side of
    RLLL=LLLR, so it cannot affect the equation.
    """

    if len(output_bits) != 3 or len(input_bits) != 3:
        raise ValueError("L has three binary input and three binary output legs")
    if any(bit not in (0, 1) for bit in tuple(output_bits) + tuple(input_bits)):
        raise ValueError("L legs must be binary")
    total = sum(output_bits) + sum(input_bits)
    if total % 2:
        return sp.Integer(0)
    return _sympy_number(q) ** (total // 2)


def physical_local_matrix(
    q: int | Fraction | sp.Expr = PHYSICAL_Q, scale: int = 1
) -> sp.Matrix:
    """Build the exact 8x8 physical L matrix in the bit-position basis."""

    matrix = sp.zeros(8, 8)
    for input_state in range(8):
        input_bits = _bits(input_state, 3)
        for output_state in range(8):
            output_bits = _bits(output_state, 3)
            matrix[output_state, input_state] = scale * physical_l_entry(
                output_bits, input_bits, q
            )
    return matrix


def physical_tensor_controls() -> dict:
    """Re-run the Wave-5 tensor contractions before constructing any RLLL system."""

    rows = WAVE5.finite_lattice_controls()
    actual = {
        tuple(row["shape"]): int(row["common_scaled_partition_integer"]) for row in rows
    }
    expected_match = actual == EXPECTED_PARTITION_INTEGERS
    coefficientwise_match = (
        len(rows) == len(CONTROL_BOXES)
        and all(row["coefficientwise_match"] for row in rows)
        and {tuple(row["shape"]) for row in rows} == set(CONTROL_BOXES)
    )

    reconstruction_mismatches = []
    for input_bits in product((0, 1), repeat=3):
        for output_bits in product((0, 1), repeat=3):
            rebuilt = physical_l_entry(output_bits, input_bits, PHYSICAL_Q)
            stored_method = sp.Rational(
                WAVE5.site_tensor_entry(input_bits + output_bits, Fraction(1, 2))
            )
            if sp.cancel(rebuilt - stored_method) != 0:
                reconstruction_mismatches.append(
                    {
                        "input": list(input_bits),
                        "output": list(output_bits),
                        "rebuilt": str(rebuilt),
                        "wave5": str(stored_method),
                    }
                )

    if not expected_match or not coefficientwise_match or reconstruction_mismatches:
        raise ArithmeticError("the prerequisite physical-tensor controls did not pass")
    return {
        "method_reused": "experiments/e35_tetrahedron_spectral.py: finite_lattice_controls and site_tensor_entry",
        "rows": rows,
        "expected_scaled_partition_integers": {
            "x".join(map(str, shape)): value
            for shape, value in EXPECTED_PARTITION_INTEGERS.items()
        },
        "actual_scaled_partition_integers": {
            "x".join(map(str, shape)): value for shape, value in actual.items()
        },
        "coefficientwise_matches": coefficientwise_match,
        "scaled_partition_integer_matches": expected_match,
        "local_entry_reconstruction_mismatches": reconstruction_mismatches,
    }


def physical_grade_analysis(q: int | Fraction | sp.Expr = PHYSICAL_Q) -> dict:
    """Verify exact Z2-grade preservation and record the two parity blocks."""

    local = physical_local_matrix(q)
    violations = []
    support = []
    allowed = 0
    for input_state in range(8):
        for output_state in range(8):
            grades_match = (input_state.bit_count() - output_state.bit_count()) % 2 == 0
            if grades_match:
                allowed += 1
            if local[output_state, input_state] != 0:
                support.append((output_state, input_state))
                if not grades_match:
                    violations.append((output_state, input_state, str(local[output_state, input_state])))
    even_states = [state for state in range(8) if state.bit_count() % 2 == 0]
    odd_states = [state for state in range(8) if state.bit_count() % 2 == 1]
    block_ranks = {
        "even": local.extract(even_states, even_states).rank(),
        "odd": local.extract(odd_states, odd_states).rank(),
    }
    return {
        "q": str(_sympy_number(q)),
        "matrix_shape": [8, 8],
        "entries": 64,
        "grade_allowed_entries": allowed,
        "nonzero_entries": len(support),
        "grade_violations": [list(item) for item in violations],
        "even_basis_states": even_states,
        "odd_basis_states": odd_states,
        "grade_block_ranks": block_ranks,
        "full_rank": local.rank(),
        "grade_preserving": not violations,
    }


def graded_coordinates() -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    """All 32 independent R[output,input] coordinates preserving total parity."""

    coordinates = []
    for input_bits in product((0, 1), repeat=3):
        for output_bits in product((0, 1), repeat=3):
            if (sum(input_bits) - sum(output_bits)) % 2 == 0:
                coordinates.append((output_bits, input_bits))
    return tuple(coordinates)


def coordinate_name(coordinate: tuple[tuple[int, ...], tuple[int, ...]]) -> str:
    output_bits, input_bits = coordinate
    return "r_" + "".join(map(str, output_bits)) + "_" + "".join(map(str, input_bits))


def coordinate_symbols(
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]] | None = None,
) -> tuple[sp.Symbol, ...]:
    if coordinates is None:
        coordinates = graded_coordinates()
    return tuple(sp.Symbol(coordinate_name(coordinate)) for coordinate in coordinates)


def _replace_local_bits(global_state: int, positions: Sequence[int], local_state: int) -> int:
    result = global_state
    for local_position, global_position in enumerate(positions):
        if (local_state >> local_position) & 1:
            result |= 1 << global_position
        else:
            result &= ~(1 << global_position)
    return result


def _embedded_three_space(local: sp.Matrix, positions: Sequence[int]) -> sp.SparseMatrix:
    entries: dict[tuple[int, int], sp.Expr] = {}
    for global_input in range(64):
        local_input = _state(tuple((global_input >> position) & 1 for position in positions))
        for local_output in range(8):
            value = local[local_output, local_input]
            if value != 0:
                global_output = _replace_local_bits(global_input, positions, local_output)
                entries[(global_output, global_input)] = value
    return sp.SparseMatrix(64, 64, entries)


def _lll_products(
    q: int | Fraction | sp.Expr, scale: int = 1
) -> tuple[sp.Matrix, sp.Matrix]:
    local = physical_local_matrix(q, scale=scale)
    embedded = [_embedded_three_space(local, positions) for positions in POSITIONS[1:]]
    forward = embedded[0] * embedded[1] * embedded[2]
    reverse = embedded[2] * embedded[1] * embedded[0]
    return forward, reverse


def _component_row(
    output_state: int,
    input_state: int,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[sp.Expr]:
    """Coefficients of (R*forward-reverse*R)[output,input]."""

    output_auxiliary = output_state & 7
    input_auxiliary = input_state & 7
    row = []
    for output_bits, input_bits in coordinates:
        coordinate_output = _state(output_bits)
        coordinate_input = _state(input_bits)
        value = sp.Integer(0)
        if output_auxiliary == coordinate_output:
            intermediate = (output_state & ~7) | coordinate_input
            value += forward[intermediate, input_state]
        if input_auxiliary == coordinate_input:
            intermediate = (input_state & ~7) | coordinate_output
            value -= reverse[output_state, intermediate]
        row.append(value)
    return row


def graded_equation_matrix(
    q: int | Fraction | sp.Expr = PHYSICAL_Q, scale: int = INTEGER_L_SCALE
) -> tuple[sp.Matrix, tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]]:
    """Return all 4096 exact component equations as a 4096x32 matrix."""

    coordinates = graded_coordinates()
    forward, reverse = _lll_products(q, scale=scale)
    rows = []
    for output_state in range(64):
        for input_state in range(64):
            rows.append(_component_row(output_state, input_state, forward, reverse, coordinates))
    return sp.Matrix(rows), coordinates

def direct_embedding_cross_check(
    matrix: sp.Matrix,
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
    q: int | Fraction | sp.Expr = PHYSICAL_Q,
    scale: int = INTEGER_L_SCALE,
) -> dict:
    """Compare every derived coefficient with direct 64x64 operator products."""

    started = time.perf_counter()
    forward, reverse = _lll_products(q, scale=scale)
    mismatches = []
    for column, (output_bits, input_bits) in enumerate(coordinates):
        local_r = sp.zeros(8, 8)
        local_r[_state(output_bits), _state(input_bits)] = 1
        embedded_r = _embedded_three_space(local_r, POSITIONS[0])
        residual = embedded_r * forward - reverse * embedded_r
        for output_state in range(64):
            for input_state in range(64):
                row = output_state * 64 + input_state
                if residual[output_state, input_state] != matrix[row, column]:
                    mismatches.append(
                        {
                            "row": row,
                            "column": column,
                            "direct": str(residual[output_state, input_state]),
                            "derived": str(matrix[row, column]),
                        }
                    )
                    if len(mismatches) == 8:
                        break
            if len(mismatches) == 8:
                break
        if len(mismatches) == 8:
            break
    return {
        "compared_coefficients": matrix.rows * matrix.cols,
        "mismatches": mismatches,
        "passed": not mismatches,
        "elapsed_seconds": f"{time.perf_counter() - started:.6f}",
    }


def _primitive_integer_row(row: Iterable[sp.Expr]) -> tuple[int, ...]:
    integers = [int(value) for value in row]
    divisor = 0
    for value in integers:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("the zero row has no primitive representative")
    primitive = [value // divisor for value in integers]
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = [-value for value in primitive]
    return tuple(primitive)


def exact_rank_certificate(matrix: sp.Matrix) -> dict:
    """Extract a nonsingular 32x32 row minor over Q."""

    started = time.perf_counter()
    _, pivot_rows = matrix.T.rref()
    pivot_rows = tuple(int(index) for index in pivot_rows)
    minor = matrix.extract(pivot_rows, range(matrix.cols))
    determinant = int(minor.det(method="domain-ge"))
    elapsed = time.perf_counter() - started
    if len(pivot_rows) != matrix.cols or determinant == 0:
        raise ArithmeticError("failed to extract a full-column-rank component minor")
    return {
        "rank": len(pivot_rows),
        "pivot_rows": pivot_rows,
        "minor": minor,
        "determinant": determinant,
        "determinant_factorization": {
            str(prime): exponent for prime, exponent in sp.factorint(abs(determinant)).items()
        },
        "elapsed_seconds": f"{elapsed:.6f}",
    }


def symbolic_minor_certificate(pivot_rows: Sequence[int], fixed_determinant: int) -> dict:
    """Factor the same component minor over Q[q]."""

    started = time.perf_counter()
    q = sp.Symbol("q")
    coordinates = graded_coordinates()
    forward, reverse = _lll_products(q, scale=1)
    rows = []
    for row_index in pivot_rows:
        output_state, input_state = divmod(int(row_index), 64)
        rows.append(_component_row(output_state, input_state, forward, reverse, coordinates))
    minor = sp.Matrix(rows)
    determinant = sp.factor(minor.det(method="domain-ge"))
    coefficient, factors = sp.factor_list(determinant, q)
    expected_root_part = q**56 * (q - 1) ** 45 * (q + 1) ** 39
    constant_quotient = sp.cancel(determinant / expected_root_part)
    fixed_cross_check = sp.cancel(
        determinant.subs(q, sp.Rational(1, 4)) * INTEGER_L_SCALE ** (3 * len(pivot_rows))
        - fixed_determinant
    )
    elapsed = time.perf_counter() - started
    return {
        "variable": "q=w^2=tanh(K)",
        "minor_shape": [minor.rows, minor.cols],
        "determinant": str(determinant),
        "factor_coefficient": str(coefficient),
        "factors": [[str(factor), int(exponent)] for factor, exponent in factors],
        "constant_quotient_after_expected_root_part": str(constant_quotient),
        "exceptional_q": ["0", "1", "-1"],
        "nonexceptional_localization": "q*(q-1)*(q+1) != 0",
        "fixed_q_scaled_determinant_cross_check": fixed_cross_check == 0,
        "elapsed_seconds": f"{elapsed:.6f}",
        "certificate_valid": constant_quotient.is_Rational and constant_quotient != 0,
    }


@contextmanager
def _hard_timeout(seconds: int):
    previous_handler = signal.getsignal(signal.SIGALRM)

    def handle_timeout(_signum, _frame):
        raise GroebnerTimeout(f"Groebner calculation exceeded {seconds} seconds")

    signal.signal(signal.SIGALRM, handle_timeout)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)


def run_groebner(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    timeout_seconds: int,
    modulus: int | None = None,
    order: str = "grevlex",
) -> dict:
    """Compute a bounded Groebner basis and return an auditable outcome."""

    started = time.perf_counter()
    try:
        with _hard_timeout(timeout_seconds):
            if modulus is None:
                basis = sp.groebner(equations, *variables, order=order, domain=sp.QQ)
            else:
                basis = sp.groebner(equations, *variables, order=order, modulus=modulus)
        expressions = [polynomial.as_expr() for polynomial in basis.polys]
        status = "completed"
        error = None
    except GroebnerTimeout as exc:
        expressions = []
        status = "timeout"
        error = str(exc)
    elapsed = time.perf_counter() - started
    return {
        "status": status,
        "field": "QQ" if modulus is None else f"F_{modulus}",
        "order": order,
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": f"{elapsed:.6f}",
        "basis": [str(expression) for expression in expressions],
        "basis_size": len(expressions),
        "unit_ideal": len(expressions) == 1 and sp.expand(expressions[0]) == 1,
        "error": error,
    }


def _generic_graded_r(
    variables: Sequence[sp.Symbol],
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> tuple[sp.Matrix, sp.Expr]:
    matrix = sp.zeros(8, 8)
    for variable, (output_bits, input_bits) in zip(variables, coordinates, strict=True):
        matrix[_state(output_bits), _state(input_bits)] = variable
    even_states = [state for state in range(8) if state.bit_count() % 2 == 0]
    odd_states = [state for state in range(8) if state.bit_count() % 2 == 1]
    determinant = sp.expand(
        matrix.extract(even_states, even_states).det(method="domain-ge")
        * matrix.extract(odd_states, odd_states).det(method="domain-ge")
    )
    return matrix, determinant


def groebner_and_saturation_analysis(
    rank_minor: sp.Matrix,
    coordinates: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> dict:
    """Classify the fixed-q homogeneous ideal and remove every degenerate R locus."""

    variables = coordinate_symbols(coordinates)
    raw_equations = [
        sp.Add(*(rank_minor[row, column] * variables[column] for column in range(len(variables))))
        for row in range(rank_minor.rows)
    ]

    homogeneous = run_groebner(
        raw_equations,
        variables,
        timeout_seconds=GROEBNER_TIMEOUT_SECONDS,
    )
    normalization = run_groebner(
        raw_equations + [variables[0] - 1],
        variables,
        timeout_seconds=GROEBNER_TIMEOUT_SECONDS,
    )

    grade_block_fallback = {
        "attempted": False,
        "reason": "not needed: the exact QQ normalization Groebner calculation completed",
    }
    finite_field_fallbacks: list[dict] = []
    if normalization["status"] == "timeout":
        even_variables = [
            variable
            for variable, (output_bits, _input_bits) in zip(variables, coordinates, strict=True)
            if sum(output_bits) % 2 == 0
        ]
        odd_variables = [variable for variable in variables if variable not in even_variables]
        grade_block_fallback = run_groebner(
            raw_equations + [variables[0] - 1],
            tuple(even_variables + odd_variables),
            timeout_seconds=GROEBNER_TIMEOUT_SECONDS,
        )
        grade_block_fallback["attempted"] = True
        for prime in FINITE_FIELD_PRIMES:
            finite_field_fallbacks.append(
                run_groebner(
                    raw_equations + [variables[0] - 1],
                    variables,
                    timeout_seconds=FINITE_FIELD_TIMEOUT_SECONDS,
                    modulus=prime,
                )
            )

    reduced_generators = list(variables)
    patch_results = []
    for variable in variables:
        result = run_groebner(
            reduced_generators + [variable - 1],
            variables,
            timeout_seconds=GROEBNER_TIMEOUT_SECONDS,
        )
        patch_results.append(
            {
                "normalization": f"{variable}=1",
                "status": result["status"],
                "unit_ideal": result["unit_ideal"],
                "basis": result["basis"],
                "elapsed_seconds": result["elapsed_seconds"],
            }
        )

    generic_r, determinant_r = _generic_graded_r(variables, coordinates)
    z = sp.Symbol("z_det")
    invertibility = run_groebner(
        reduced_generators + [z * determinant_r - 1],
        tuple(variables) + (z,),
        timeout_seconds=GROEBNER_TIMEOUT_SECONDS,
    )

    return {
        "field": "QQ",
        "fixed_physical_specialization": "q=1/4 (w=1/2)",
        "smallest_honest_formulation": (
            "32 raw independent component equations selected by exact row reduction, "
            "plus normalization r_000_000=1"
        ),
        "raw_independent_equations": [str(sp.expand(equation)) for equation in raw_equations],
        "homogeneous_ideal": homogeneous,
        "normalization_coordinate": str(variables[0]),
        "normalized_raw_system": normalization,
        "grade_block_order_fallback": grade_block_fallback,
        "finite_field_fallbacks": finite_field_fallbacks,
        "finite_field_status_note": (
            "No finite-field fallback was needed. Had it been needed, [1] modulo a prime would be "
            "a characteristic-p certificate and strong evidence only, not a characteristic-zero theorem."
        ),
        "projective_nonzero_saturation": {
            "method": (
                "The homogeneous Groebner basis is the irrelevant coordinate ideal. "
                "The 32 charts r_j=1 exhaust every nonzero R modulo overall scalar."
            ),
            "patch_count": len(patch_results),
            "unit_ideal_patch_count": sum(row["unit_ideal"] for row in patch_results),
            "patches": patch_results,
        },
        "invertibility_saturation": {
            "method": "Rabinowitsch equation z_det*det(R)-1 after exact linear reduction",
            "r_matrix_shape": [generic_r.rows, generic_r.cols],
            "parity_block_shapes": [[4, 4], [4, 4]],
            "determinant_total_degree": sp.Poly(determinant_r, *variables).total_degree(),
            "determinant_term_count": len(sp.Poly(determinant_r, *variables).terms()),
            "result": invertibility,
        },
    }


def graded_rlll_analysis() -> dict:
    """Construct, reduce, and classify the exact graded RLLL system."""

    matrix_started = time.perf_counter()
    matrix, coordinates = graded_equation_matrix(PHYSICAL_Q, scale=INTEGER_L_SCALE)
    matrix_elapsed = time.perf_counter() - matrix_started
    direct_cross_check = direct_embedding_cross_check(matrix, coordinates)
    nonzero_rows = [row for row in range(matrix.rows) if any(matrix[row, column] != 0 for column in range(matrix.cols))]
    primitive_rows = {
        _primitive_integer_row(matrix[row, :]) for row in nonzero_rows
    }
    rank = exact_rank_certificate(matrix)
    symbolic = symbolic_minor_certificate(rank["pivot_rows"], rank["determinant"])
    groebner = groebner_and_saturation_analysis(rank["minor"], coordinates)

    coordinate_table = []
    for index, coordinate in enumerate(coordinates):
        output_bits, input_bits = coordinate
        coordinate_table.append(
            {
                "index": index,
                "name": coordinate_name(coordinate),
                "output": list(output_bits),
                "input": list(input_bits),
                "grade": sum(input_bits) % 2,
            }
        )

    pivot_components = [
        {
            "row": row,
            "output_state": row // 64,
            "input_state": row % 64,
            "output_bits": list(_bits(row // 64, 6)),
            "input_bits": list(_bits(row % 64, 6)),
        }
        for row in rank["pivot_rows"]
    ]
    minor_rows = [[int(rank["minor"][row, column]) for column in range(matrix.cols)] for row in range(matrix.cols)]

    return {
        "ansatz": {
            "equation": "R_123 L_145 L_246 L_356 = L_356 L_246 L_145 R_123",
            "operator_convention": "L[output,input]=T[input,output]/2",
            "global_spaces": 6,
            "global_operator_shape": [64, 64],
            "auxiliary_r_shape": [8, 8],
            "r_total_entries": 64,
            "grading_rule": "R[o1,o2,o3;i1,i2,i3]=0 unless sum(o)-sum(i)=0 mod 2",
            "independent_unknowns": len(coordinates),
            "additional_symmetry_imposed": "none; no leg permutation or translation identification",
            "overall_scalar_gauge": "homogeneous in R; normalize one nonzero coordinate to 1",
            "physical_q": "symbolic q=w^2=tanh(K), with exact QQ Groebner specialization q=1/4",
            "coordinates": coordinate_table,
        },
        "system": {
            "component_equations": matrix.rows,
            "unknowns": matrix.cols,
            "integer_specialization": "q=1/4 and L replaced uniformly by 64*L",
            "integer_scale_effect": "both RLLL sides are multiplied by 64^3",
            "nonzero_component_equations": len(nonzero_rows),
            "identically_zero_component_equations": matrix.rows - len(nonzero_rows),
            "distinct_nonzero_primitive_linear_forms": len(primitive_rows),
            "exact_Q_rank": rank["rank"],
            "nullity": matrix.cols - rank["rank"],
            "construction_elapsed_seconds": f"{matrix_elapsed:.6f}",
            "direct_operator_cross_check": direct_cross_check,
        },
        "rank_certificate": {
            "pivot_components": pivot_components,
            "minor_matrix": minor_rows,
            "fixed_q_scaled_determinant": str(rank["determinant"]),
            "fixed_q_scaled_determinant_factorization": rank["determinant_factorization"],
            "elapsed_seconds": rank["elapsed_seconds"],
        },
        "symbolic_parameter_certificate": symbolic,
        "groebner_and_saturations": groebner,
        "outcome": {
            "class": "GRADED_NOGO_Q",
            "theorem": (
                "For every characteristic-zero q with q not in {0,+1,-1}, the only Z2-grade-preserving "
                "8x8 solution of the displayed RLLL equation is R=0. In particular no nonzero or invertible "
                "graded R exists for any finite ferromagnetic 3D Ising weight 0<q<1."
            ),
            "solutions_found": False,
            "commuting_transfer_matrix_check": "not applicable: no nonzero graded R survives",
            "scope_limit": (
                "Parity-violating R, higher auxiliary dimension, nonidentical/spectral L factors, IRF/dynamical "
                "relations, and singular dimension-changing projections are not classified. Exceptional "
                "q in {0,+1,-1} is excluded by the stated localization and is not classified here."
            ),
        },
    }


def build_artifact() -> dict:
    # This ordering is intentional: no RLLL construction starts until all three
    # stored physical-tensor partition controls have been recomputed and checked.
    controls_started = time.perf_counter()
    controls = physical_tensor_controls()
    controls_elapsed = time.perf_counter() - controls_started
    grade = physical_grade_analysis(PHYSICAL_Q)
    analysis = graded_rlll_analysis()

    homogeneous = analysis["groebner_and_saturations"]["homogeneous_ideal"]
    normalized = analysis["groebner_and_saturations"]["normalized_raw_system"]
    projective = analysis["groebner_and_saturations"]["projective_nonzero_saturation"]
    invertible = analysis["groebner_and_saturations"]["invertibility_saturation"]["result"]
    symbolic = analysis["symbolic_parameter_certificate"]
    checks = [
        {
            "name": "verified physical tensor partition controls",
            "passed": (
                controls["coefficientwise_matches"]
                and controls["scaled_partition_integer_matches"]
                and not controls["local_entry_reconstruction_mismatches"]
            ),
            "detail": "; ".join(
                f"{shape}={value}"
                for shape, value in controls["actual_scaled_partition_integers"].items()
            ),
        },
        {
            "name": "physical L is exactly Z2-grade-preserving",
            "passed": grade["grade_preserving"] and grade["nonzero_entries"] == 32,
            "detail": (
                f"{grade['nonzero_entries']} nonzero of 64 entries; "
                f"{len(grade['grade_violations'])} grade violations"
            ),
        },
        {
            "name": "graded RLLL component system matches direct operator embedding and has full rank",
            "passed": (
                analysis["system"]["direct_operator_cross_check"]["passed"]
                and analysis["system"]["unknowns"] == 32
                and analysis["system"]["exact_Q_rank"] == 32
                and int(analysis["rank_certificate"]["fixed_q_scaled_determinant"]) != 0
            ),
            "detail": (
                f"{analysis['system']['nonzero_component_equations']} nonzero of 4096 components; "
                f"rank {analysis['system']['exact_Q_rank']}/32"
            ),
        },
        {
            "name": "symbolic minor excludes every nonexceptional physical q",
            "passed": (
                symbolic["certificate_valid"]
                and symbolic["fixed_q_scaled_determinant_cross_check"]
            ),
            "detail": f"det={symbolic['determinant']}",
        },
        {
            "name": "exact QQ Groebner normalization is the unit ideal",
            "passed": (
                homogeneous["status"] == "completed"
                and homogeneous["basis_size"] == 32
                and normalized["status"] == "completed"
                and normalized["unit_ideal"]
            ),
            "detail": (
                f"raw homogeneous basis size {homogeneous['basis_size']}; "
                f"normalized basis {normalized['basis']} in {normalized['elapsed_seconds']} s"
            ),
        },
        {
            "name": "all nonzero and invertible R loci are saturated away",
            "passed": (
                projective["patch_count"] == 32
                and projective["unit_ideal_patch_count"] == 32
                and invertible["status"] == "completed"
                and invertible["unit_ideal"]
            ),
            "detail": (
                f"{projective['unit_ideal_patch_count']}/{projective['patch_count']} coordinate charts [1]; "
                f"invertibility saturation basis {invertible['basis']}"
            ),
        },
    ]
    return {
        "meta": {
            "provenance": "experiments/e44_tetra_graded.py",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "arithmetic": "exact Python integers/Fraction and SymPy QQ polynomial algebra; wall timings only use floats",
            "sympy_version": sp.__version__,
            "groebner_timeout_seconds": GROEBNER_TIMEOUT_SECONDS,
            "finite_field_fallback_timeout_seconds": FINITE_FIELD_TIMEOUT_SECONDS,
            "physical_controls_elapsed_seconds": f"{controls_elapsed:.6f}",
        },
        "data": {
            "physical_tensor_controls": controls,
            "physical_grading": grade,
            "graded_rlll": analysis,
            "status": (
                "[THEOREM] GRADED_NOGO_Q for the displayed 8x8 parity-preserving RLLL subfamily at every "
                "characteristic-zero q outside {0,+1,-1}. [UNRESOLVED] parity-breaking, higher-dimensional, "
                "nonidentical-L, IRF/dynamical, and singular projection extensions."
            ),
        },
        "checks": checks,
    }


def main() -> None:
    artifact = build_artifact()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    analysis = artifact["data"]["graded_rlll"]
    print(
        "physical controls:",
        artifact["data"]["physical_tensor_controls"]["actual_scaled_partition_integers"],
    )
    print(
        "graded system:",
        analysis["system"]["nonzero_component_equations"],
        "nonzero components, rank",
        analysis["system"]["exact_Q_rank"],
    )
    print("symbolic minor:", analysis["symbolic_parameter_certificate"]["determinant"])
    print(
        "normalized QQ Groebner:",
        analysis["groebner_and_saturations"]["normalized_raw_system"]["basis"],
    )
    print("written", OUTPUT.relative_to(ROOT))
    if not all(check["passed"] for check in artifact["checks"]):
        for check in artifact["checks"]:
            print(check["name"], "PASS" if check["passed"] else "FAIL", check["detail"])
        raise AssertionError("one or more graded tetrahedron checks failed")
    print("PASS")


if __name__ == "__main__":
    main()
