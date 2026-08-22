"""Standalone exact verifier for the cabled tetra16 specialization-locus certificate.

The producer is deliberately not imported.  This verifier rebuilds raw seven-bit
operators, each stored 128x128 minor, fresh determinant values, inherited
full-rank controls, the complementary-minor gcd, and every listed rational
exceptional specialization.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcdex
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "integrability" / "tetra16_locus.json"
SOURCE = ROOT / "results" / "integrability" / "tetra16.json"
Q = sp.Symbol("q")
DIMENSION = 128
AUXILIARY_MASK = (1 << 4) - 1
SECTOR_SIZE = 128
SCALE_EXPONENT = 12 * SECTOR_SIZE
PRIME_CEILING = 2_000_000_000


def _state(bits: tuple[int, ...]) -> int:
    return sum(int(bit) << position for position, bit in enumerate(bits))


def _replace_local_bits(global_state: int, positions: tuple[int, ...], local_state: int) -> int:
    result = global_state
    for local_position, global_position in enumerate(positions):
        if (local_state >> local_position) & 1:
            result |= 1 << global_position
        else:
            result &= ~(1 << global_position)
    return result


def _integer_local_matrix(q_value: Fraction) -> sp.SparseMatrix:
    """Return denominator^3 times the exact binary Ising L(q)."""

    entries: dict[tuple[int, int], int] = {}
    numerator = q_value.numerator
    denominator = q_value.denominator
    for input_state in range(8):
        for output_state in range(8):
            total = input_state.bit_count() + output_state.bit_count()
            if total % 2:
                continue
            exponent = total // 2
            entries[(output_state, input_state)] = numerator**exponent * denominator ** (3 - exponent)
    return sp.SparseMatrix(8, 8, entries)


def _embedded_three_space(local: sp.SparseMatrix, positions: tuple[int, ...]) -> sp.SparseMatrix:
    entries: dict[tuple[int, int], int] = {}
    for global_input in range(DIMENSION):
        local_input = _state(tuple((global_input >> position) & 1 for position in positions))
        for local_output in range(8):
            value = local[local_output, local_input]
            if value:
                entries[(_replace_local_bits(global_input, positions, local_output), global_input)] = int(value)
    return sp.SparseMatrix(DIMENSION, DIMENSION, entries)


def _products(q_value: Fraction) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix]:
    """Build scaled forward/reverse cabled products independently of the producer."""

    local = _integer_local_matrix(q_value)
    first = _embedded_three_space(local, (0, 4, 5))
    second = _embedded_three_space(local, (1, 4, 5))
    cable = first * second
    reverse_cable = second * first
    middle = _embedded_three_space(local, (2, 4, 6))
    last = _embedded_three_space(local, (3, 5, 6))
    return cable * middle * last, last * middle * cable, cable, reverse_cable

def _polynomial_local_matrix() -> sp.SparseMatrix:
    entries: dict[tuple[int, int], sp.Expr] = {}
    for input_state in range(8):
        for output_state in range(8):
            total = input_state.bit_count() + output_state.bit_count()
            if total % 2 == 0:
                entries[(output_state, input_state)] = Q ** (total // 2)
    return sp.SparseMatrix(8, 8, entries)


def _embedded_three_space_polynomial(
    local: sp.SparseMatrix, positions: tuple[int, ...]
) -> sp.SparseMatrix:
    entries: dict[tuple[int, int], sp.Expr] = {}
    for global_input in range(DIMENSION):
        local_input = _state(tuple((global_input >> position) & 1 for position in positions))
        for local_output in range(8):
            value = local[local_output, local_input]
            if value:
                entries[(_replace_local_bits(global_input, positions, local_output), global_input)] = value
    return sp.SparseMatrix(DIMENSION, DIMENSION, entries)


def _symbolic_products() -> tuple[sp.Matrix, sp.Matrix]:
    """Rebuild the unscaled Z[q] forward/reverse products without the producer."""

    local = _polynomial_local_matrix()
    first = _embedded_three_space_polynomial(local, (0, 4, 5))
    second = _embedded_three_space_polynomial(local, (1, 4, 5))
    cable = first * second
    middle = _embedded_three_space_polynomial(local, (2, 4, 6))
    last = _embedded_three_space_polynomial(local, (3, 5, 6))
    return cable * middle * last, last * middle * cable


def _minor_entries(
    rows: list[int],
    sector: str,
    forward: sp.Matrix,
    reverse: sp.Matrix,
) -> tuple[tuple[tuple[int, int, tuple[int, ...]], ...], int, int, int]:
    """Return expanded nonzero entries and independent tropical bounds."""

    coordinates = SECTOR_COORDINATES[sector]
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    infinite = 10**6
    low = np.full((SECTOR_SIZE, SECTOR_SIZE), infinite, dtype=np.int64)
    negative_high = np.full((SECTOR_SIZE, SECTOR_SIZE), infinite, dtype=np.int64)
    entries: list[tuple[int, int, tuple[int, ...]]] = []
    row_squares = [0] * SECTOR_SIZE
    for matrix_row, row_index in enumerate(rows):
        output_state, input_state = divmod(int(row_index), DIMENSION)
        output_auxiliary = output_state & AUXILIARY_MASK
        input_auxiliary = input_state & AUXILIARY_MASK
        values: list[sp.Expr] = [sp.Integer(0)] * SECTOR_SIZE
        for coordinate_input in range(16):
            column = coordinate_index.get((output_auxiliary, coordinate_input))
            if column is not None:
                intermediate = (output_state & ~AUXILIARY_MASK) | coordinate_input
                values[column] += forward[intermediate, input_state]
        for coordinate_output in range(16):
            column = coordinate_index.get((coordinate_output, input_auxiliary))
            if column is not None:
                intermediate = (input_state & ~AUXILIARY_MASK) | coordinate_output
                values[column] -= reverse[output_state, intermediate]
        for column, value in enumerate(values):
            polynomial = sp.Poly(sp.expand(value), Q, domain=sp.ZZ)
            if polynomial.is_zero:
                continue
            degree = int(polynomial.degree())
            coefficients = [0] * (degree + 1)
            for (power,), coefficient in polynomial.as_dict().items():
                coefficients[int(power)] = int(coefficient)
            valuation = next(power for power, coefficient in enumerate(coefficients) if coefficient)
            entries.append((matrix_row, column, tuple(coefficients)))
            low[matrix_row, column] = valuation
            negative_high[matrix_row, column] = -degree
            row_squares[matrix_row] += sum(abs(coefficient) for coefficient in coefficients) ** 2
    assert all(value > 0 for value in row_squares)
    minimum_rows, minimum_columns = linear_sum_assignment(low)
    maximum_rows, maximum_columns = linear_sum_assignment(negative_high)
    assert len(minimum_rows) == SECTOR_SIZE == len(maximum_rows)
    minimum = sum(int(low[row, column]) for row, column in zip(minimum_rows, minimum_columns, strict=True))
    maximum = -sum(
        int(negative_high[row, column])
        for row, column in zip(maximum_rows, maximum_columns, strict=True)
    )
    assert all(int(low[row, column]) < infinite for row, column in zip(minimum_rows, minimum_columns, strict=True))
    assert all(
        int(negative_high[row, column]) < infinite
        for row, column in zip(maximum_rows, maximum_columns, strict=True)
    )
    product = math.prod(row_squares)
    height = math.isqrt(product)
    if height * height < product:
        height += 1
    return tuple(entries), minimum, maximum, height


def _modular_determinant(
    entries: tuple[tuple[int, int, tuple[int, ...]], ...],
    node: int,
    prime: int,
) -> int:
    """Exact F_p determinant, using only range-checked int64 lane operations."""

    assert prime < PRIME_CEILING and prime * prime < 2**63
    maximum_degree = max(len(coefficients) - 1 for _, _, coefficients in entries)
    powers = [1]
    for _ in range(maximum_degree):
        powers.append(powers[-1] * (node % prime) % prime)
    matrix = np.zeros((SECTOR_SIZE, SECTOR_SIZE), dtype=np.int64)
    for row, column, coefficients in entries:
        matrix[row, column] = sum(
            coefficient * powers[power] for power, coefficient in enumerate(coefficients)
        ) % prime
    determinant = 1
    for column in range(SECTOR_SIZE):
        pivot = next((row for row in range(column, SECTOR_SIZE) if matrix[row, column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[[column, pivot]] = matrix[[pivot, column]]
            determinant = -determinant
        pivot_value = int(matrix[column, column])
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, prime - 2, prime)
        multipliers = (matrix[column + 1 :, column] * inverse) % prime
        matrix[column + 1 :, column + 1 :] = (
            matrix[column + 1 :, column + 1 :]
            - multipliers[:, None] * matrix[column, column + 1 :]
        ) % prime
        matrix[column + 1 :, column] = 0
    return int(determinant % prime)


def _verify_modular_interpolation(record: dict, entries: tuple[tuple[int, int, tuple[int, ...]], ...], valuation: int, degree: int, height: int) -> int:
    """Independently prove the stored shifted coefficients equal the determinant."""

    point_count = degree - valuation + 1
    coefficients = [int(value) for value in record["shifted_coefficients_low_to_high"]]
    assert record["valuation_lower_bound"] == valuation
    assert record["determinant_degree_upper_bound"] == degree
    assert record["shifted_degree_upper_bound"] == point_count - 1
    assert len(coefficients) == point_count
    assert max(abs(value) for value in coefficients) <= height
    assert int(record["height_bound"]) >= height
    all_primes = [int(value) for value in record["crt_primes"]]
    assert len(all_primes) == int(record["crt_prime_count"])
    assert len(all_primes) == len(set(all_primes))
    full_modulus = math.prod(all_primes)
    assert str(full_modulus) == record["crt_modulus"]
    modulus = 1
    selected_primes = []
    for prime in all_primes:
        assert sp.isprime(prime)
        assert prime < PRIME_CEILING and prime * prime < 2**63
        assert prime > point_count
        selected_primes.append(prime)
        modulus *= prime
        if modulus > 2 * height:
            break
    assert modulus > 2 * height
    for prime in selected_primes:
        for node in range(1, point_count + 1):
            determinant = _modular_determinant(entries, node, prime)
            actual = determinant * pow(pow(node, valuation, prime), prime - 2, prime) % prime
            expected = 0
            for coefficient in reversed(coefficients):
                expected = (expected * node + coefficient) % prime
            assert actual == expected
    return len(selected_primes)


def _coordinates() -> tuple[tuple[int, int], ...]:
    preserving: list[tuple[int, int]] = []
    reversing: list[tuple[int, int]] = []
    for input_state in range(16):
        for output_state in range(16):
            target = preserving if (output_state.bit_count() - input_state.bit_count()) % 2 == 0 else reversing
            target.append((output_state, input_state))
    coordinates = tuple(preserving + reversing)
    assert len(coordinates) == 256 and len(set(coordinates)) == 256
    return coordinates


COORDINATES = _coordinates()
SECTOR_COORDINATES = {
    "parity_preserving": COORDINATES[:SECTOR_SIZE],
    "parity_reversing": COORDINATES[SECTOR_SIZE:],
}


def _component_row(
    output_state: int,
    input_state: int,
    forward: sp.Matrix,
    reverse: sp.Matrix,
    coordinates: tuple[tuple[int, int], ...],
) -> list[int]:
    output_auxiliary = output_state & AUXILIARY_MASK
    input_auxiliary = input_state & AUXILIARY_MASK
    values: list[int] = []
    for coordinate_output, coordinate_input in coordinates:
        value = 0
        if output_auxiliary == coordinate_output:
            intermediate = (output_state & ~AUXILIARY_MASK) | coordinate_input
            value += int(forward[intermediate, input_state])
        if input_auxiliary == coordinate_input:
            intermediate = (input_state & ~AUXILIARY_MASK) | coordinate_output
            value -= int(reverse[output_state, intermediate])
        values.append(value)
    return values


def _minor(
    rows: list[int],
    sector: str,
    forward: sp.Matrix,
    reverse: sp.Matrix,
) -> sp.Matrix:
    coordinates = SECTOR_COORDINATES[sector]
    return sp.Matrix(
        [
            _component_row(row // DIMENSION, row % DIMENSION, forward, reverse, coordinates)
            for row in rows
        ]
    )


def _q_value(text: str) -> Fraction:
    if "/" in text:
        numerator, denominator = map(int, text.split("/"))
        return Fraction(numerator, denominator)
    return Fraction(int(text), 1)


def _poly(record: dict) -> sp.Poly:
    coefficients = [int(value) for value in record["shifted_coefficients_low_to_high"]]
    shifted = sum(coefficient * Q**power for power, coefficient in enumerate(coefficients))
    return sp.Poly(Q ** int(record["valuation_lower_bound"]) * shifted, Q, domain=sp.ZZ)


def _linear_multiplicity(polynomial: sp.Poly, root: int) -> int:
    divisor = sp.Poly(Q - root, Q, domain=sp.ZZ)
    work = polynomial
    multiplicity = 0
    while True:
        quotient, remainder = sp.div(work, divisor, domain=sp.ZZ)
        if not remainder.is_zero:
            return multiplicity
        work = quotient
        multiplicity += 1


def _common_linear_factor(primary: sp.Poly, alternate: sp.Poly) -> sp.Poly:
    expression = sp.Integer(1)
    for root in (0, 1, -1):
        multiplicity = min(
            _linear_multiplicity(primary, root),
            _linear_multiplicity(alternate, root),
        )
        expression *= (Q - root) ** multiplicity
    return sp.Poly(sp.expand(expression), Q, domain=sp.ZZ)


def _modular_residual_gcd(
    primary: sp.Poly, alternate: sp.Poly, common: sp.Poly, prime: int
) -> list[int]:
    reduced_primary, remainder_primary = sp.div(primary, common, domain=sp.ZZ)
    reduced_alternate, remainder_alternate = sp.div(alternate, common, domain=sp.ZZ)
    assert remainder_primary.is_zero and remainder_alternate.is_zero
    primary_mod = [int(value) % prime for value in reduced_primary.all_coeffs()]
    alternate_mod = [int(value) % prime for value in reduced_alternate.all_coeffs()]
    _, _, gcd = gf_gcdex(primary_mod, alternate_mod, prime, ZZ)
    return [int(value) for value in gcd]


def _verify_modular_bezout(primary: sp.Poly, alternate: sp.Poly, common: sp.Poly, record: dict) -> None:
    reduced_primary, remainder_primary = sp.div(primary, common, domain=sp.ZZ)
    reduced_alternate, remainder_alternate = sp.div(alternate, common, domain=sp.ZZ)
    assert remainder_primary.is_zero and remainder_alternate.is_zero
    assert record["status"] == "completed"
    prime = int(record["prime"])
    assert sp.isprime(prime)
    assert int(reduced_primary.LC()) % prime == int(
        record["reduced_primary_leading_coefficient_mod_prime"]
    )
    primary_mod = sp.Poly(reduced_primary.as_expr(), Q, modulus=prime)
    alternate_mod = sp.Poly(reduced_alternate.as_expr(), Q, modulus=prime)
    coefficient_a = sp.Poly(
        sum(
            int(coefficient) * Q**power
            for power, coefficient in enumerate(
                record["bezout_a_coefficients_low_to_high_mod_prime"]
            )
        ),
        Q,
        modulus=prime,
    )
    coefficient_b = sp.Poly(
        sum(
            int(coefficient) * Q**power
            for power, coefficient in enumerate(
                record["bezout_b_coefficients_low_to_high_mod_prime"]
            )
        ),
        Q,
        modulus=prime,
    )
    identity = sp.Poly(
        coefficient_a.as_expr() * primary_mod.as_expr()
        + coefficient_b.as_expr() * alternate_mod.as_expr(),
        Q,
        modulus=prime,
    )
    assert identity == sp.Poly(1, Q, modulus=prime)
    assert _modular_residual_gcd(primary, alternate, common, prime) == [1]


def _verify_factor_record(polynomial: sp.Poly, record: dict) -> None:
    if record["status"] != "completed":
        assert record["status"] in {"not_attempted", "timed_out"}
        return
    product = sp.Integer(int(record["coefficient"]))
    for factor in record["factors"]:
        coefficients = [int(value) for value in factor["coefficients_low_to_high"]]
        expression = sum(coefficient * Q**power for power, coefficient in enumerate(coefficients))
        assert str(sp.expand(expression)) == factor["expression"]
        product *= expression ** int(factor["exponent"])
    assert sp.Poly(sp.expand(product), Q, domain=sp.ZZ) == polynomial


def _scaled_polynomial_value(polynomial: sp.Poly, q_value: Fraction) -> int:
    value = polynomial.eval(sp.Rational(q_value.numerator, q_value.denominator))
    scaled = sp.cancel(value * q_value.denominator**SCALE_EXPONENT)
    assert scaled.is_Integer
    return int(scaled)


def _full_sector_ranks(q_value: Fraction) -> dict[str, int]:
    forward, reverse, _, _ = _products(q_value)
    ranks = {}
    for sector, coordinates in SECTOR_COORDINATES.items():
        rows = [
            _component_row(output_state, input_state, forward, reverse, coordinates)
            for output_state in range(DIMENSION)
            for input_state in range(DIMENSION)
            if (output_state.bit_count() - input_state.bit_count()) % 2
            == (0 if sector == "parity_preserving" else 1)
        ]
        ranks[sector] = int(sp.Matrix(rows).rank())
    return ranks


def _verify_quadratic_full_rank_certificates(certificates: list[dict]) -> None:
    expected_residues = {
        ("i", 5): 2,
        ("i", 13): 5,
        ("-i", 5): 3,
        ("-i", 13): 8,
    }
    assert {certificate["q"] for certificate in certificates} == {"i", "-i"}
    symbolic_forward, symbolic_reverse = _symbolic_products()
    seen = set()
    for certificate in certificates:
        assert certificate["claim_tag"] == "[THEOREM]"
        assert certificate["minimal_polynomial"] == "q**2 + 1"
        assert int(certificate["full_system_rank"]) == 256
        assert int(certificate["sector_rank"]) == 128
        for reduction in certificate["reductions"]:
            prime = int(reduction["prime"])
            residue = int(reduction["q_residue"])
            assert expected_residues[(certificate["q"], prime)] == residue
            assert (residue * residue + 1) % prime == 0
            assert sp.isprime(prime)
            seen.add((certificate["q"], prime))
            for sector, sector_record in reduction["sectors"].items():
                rows = [int(row) for row in sector_record["rows"]]
                assert len(rows) == 128 and len(set(rows)) == 128
                entries, _, _, _ = _minor_entries(
                    rows, sector, symbolic_forward, symbolic_reverse
                )
                determinant = _modular_determinant(entries, residue, prime)
                assert determinant != 0
                assert determinant == int(sector_record["determinant_mod_prime"])
                assert int(sector_record["rank_lower_bound"]) == 128
                assert int(sector_record["rank_upper_bound"]) == 128
    assert seen == set(expected_residues)


def _run_check(name: str, check) -> None:
    try:
        detail = check()
    except Exception as exc:
        print(f"{name}: FAIL ({exc})")
        raise
    print(f"{name}: PASS ({detail})")


def main() -> None:
    artifact = json.loads(ARTIFACT.read_text())
    source = json.loads(SOURCE.read_text())
    data = artifact["data"]
    records = data["minor_determinants"]
    polynomials = {label: _poly(record) for label, record in records.items()}

    def check_envelope_and_interpolation() -> str:
        assert artifact["meta"]["provenance"] == "experiments/e95_tetra16_locus.py"
        assert data["source_artifact"] == "results/integrability/tetra16.json"
        interpolation = data["interpolation"]
        assert interpolation["status"] == "completed"
        assert interpolation["per_point_exact_method"] == "modular Gaussian elimination of the evaluated 128x128 minor"
        assert set(records) == {
            "primary_preserving",
            "primary_reversing",
            "alternate_preserving",
            "alternate_reversing",
        }
        symbolic_forward, symbolic_reverse = _symbolic_products()
        verified_primes = []
        for label, record in records.items():
            assert record["shape"] == [128, 128]
            rows = [int(row) for row in record["rows"]]
            assert len(rows) == 128 and len(set(rows)) == 128
            entries, valuation, degree, height = _minor_entries(
                rows, record["sector"], symbolic_forward, symbolic_reverse
            )
            assert len(entries) == int(record["entry_count"])
            assert record["actual_degree"] <= degree
            assert record["actual_shifted_degree"] <= degree - valuation
            verified_primes.append(
                _verify_modular_interpolation(record, entries, valuation, degree, height)
            )
            _verify_factor_record(polynomials[label], record["factorization"])
        return (
            "four independently rebuilt Z[q] minors agree at every interpolation "
            f"node modulo enough CRT primes ({', '.join(map(str, verified_primes))})"
        )

    def check_fresh_exact_determinants() -> str:
        # Neither value is an interpolation node (all nodes are positive integers),
        # nor is either used by the producer's own fresh-evaluation step.
        fresh_points = (Fraction(-5, 1), Fraction(5, 7))
        for q_value in fresh_points:
            forward, reverse, _, _ = _products(q_value)
            for label, record in records.items():
                determinant = int(
                    _minor([int(row) for row in record["rows"]], record["sector"], forward, reverse).det(
                        method="domain-ge"
                    )
                )
                assert determinant == _scaled_polynomial_value(polynomials[label], q_value)
        return "direct Bareiss determinants agree at q=-5 and q=5/7"

    def check_stored_full_rank_controls() -> str:
        controls = source["data"]["fixed_q_certificates"]
        for control in controls:
            q_value = _q_value(control["q"])
            forward, reverse, _, _ = _products(q_value)
            for sector, label in (
                ("parity_preserving", "primary_preserving"),
                ("parity_reversing", "primary_reversing"),
            ):
                determinant = int(
                    _minor([int(row) for row in records[label]["rows"]], sector, forward, reverse).det(
                        method="domain-ge"
                    )
                )
                assert determinant != 0
                assert determinant == int(control["sectors"][sector]["determinant"])
                assert determinant == _scaled_polynomial_value(polynomials[label], q_value)
        return "q=1/4, 1/3, and 2/5 retain both nonzero exact primary sector minors"

    def check_complementary_gcd_and_rank_drops() -> str:
        primary_full = sp.Poly(
            polynomials["primary_preserving"].as_expr()
            * polynomials["primary_reversing"].as_expr(),
            Q,
            domain=sp.ZZ,
        )
        alternate_full = sp.Poly(
            polynomials["alternate_preserving"].as_expr()
            * polynomials["alternate_reversing"].as_expr(),
            Q,
            domain=sp.ZZ,
        )
        locus = data["complementary_minor_locus"]
        candidate = _common_linear_factor(primary_full, alternate_full)
        stored_candidate = sp.Poly(
            sum(
                int(coefficient) * Q**power
                for power, coefficient in enumerate(
                    locus["candidate_linear_common_factor_coefficients_low_to_high"]
                )
            ),
            Q,
            domain=sp.ZZ,
        )
        assert candidate == stored_candidate
        candidate_factorization = locus["candidate_linear_common_factorization"]
        _verify_factor_record(candidate, candidate_factorization)
        initial_bezout = locus["initial_candidate_bezout"]
        if initial_bezout["status"] == "completed":
            _verify_modular_bezout(primary_full, alternate_full, candidate, initial_bezout)
        else:
            assert initial_bezout["status"] == "unresolved_nonconstant_modular_residual_gcd"
            rejected = initial_bezout["rejected_modular_primes"]
            assert rejected
            for row in rejected:
                prime = int(row["prime"])
                assert sp.isprime(prime)
                if row["reason"] == "nonconstant residual gcd modulo prime":
                    gcd = _modular_residual_gcd(
                        primary_full, alternate_full, candidate, prime
                    )
                    assert len(gcd) - 1 == int(row["residual_gcd_degree"])
                    assert len(gcd) > 1
        exact_record = locus["exact_common_gcd"]
        exact_coefficients = locus["exact_common_factor_coefficients_low_to_high"]
        if exact_coefficients is None:
            assert exact_record["status"] == "timed_out"
            assert not artifact["checks"][3]["passed"]
            factor_record = candidate_factorization
            rank_roots = {
                Fraction(
                    -int(factor["coefficients_low_to_high"][0]),
                    int(factor["coefficients_low_to_high"][1]),
                )
                for factor in factor_record["factors"]
            }
            assert data["unresolved_algebraic_remainders"]
            proof_detail = "candidate linear roots and an independently reproduced exact-gcd wall record"
        else:
            common = sp.Poly(
                sum(
                    int(coefficient) * Q**power
                    for power, coefficient in enumerate(exact_coefficients)
                ),
                Q,
                domain=sp.ZZ,
            )
            assert exact_record["status"] in {"certified_by_safe_prime", "completed"}
            _, remainder = sp.div(common, candidate, domain=sp.ZZ)
            assert remainder.is_zero
            _verify_modular_bezout(
                primary_full, alternate_full, common, locus["bezout"]
            )
            assert artifact["checks"][3]["passed"]
            factor_record = locus["exact_common_factorization"]
            assert factor_record is not None
            _verify_factor_record(common, factor_record)
            if factor_record["status"] == "completed":
                rank_roots = {
                    Fraction(
                        -int(factor["coefficients_low_to_high"][0]),
                        int(factor["coefficients_low_to_high"][1]),
                    )
                    for factor in factor_record["factors"]
                    if len(factor["coefficients_low_to_high"]) == 2
                }
                quadratic_factors = [
                    factor
                    for factor in factor_record["factors"]
                    if factor["expression"] == "q**2 + 1"
                ]
                if quadratic_factors:
                    assert len(quadratic_factors) == 1
                    assert int(quadratic_factors[0]["exponent"]) == 4
                    _verify_quadratic_full_rank_certificates(
                        data["algebraic_full_rank_specializations"]
                    )
                    assert artifact["checks"][4]["passed"]
                    assert not data["unresolved_algebraic_remainders"]
                    proof_detail = (
                        "exact common gcd, safe-prime Bezout identity, q=±i "
                        "two-prime full-rank squeezes, and rational rank drops"
                    )
                else:
                    proof_detail = (
                        "exact common gcd, safe-prime Bezout identity, and listed "
                        "rational rank drops"
                    )
            else:
                rank_roots = {
                    Fraction(
                        -int(factor["coefficients_low_to_high"][0]),
                        int(factor["coefficients_low_to_high"][1]),
                    )
                    for factor in candidate_factorization["factors"]
                }
                assert data["unresolved_algebraic_remainders"]
                proof_detail = "exact common gcd with factorization wall and candidate rational rank drops"
        claims = {
            _q_value(record["q"]): record
            for record in data["exceptional_specializations"]
        }
        assert rank_roots == set(claims)
        for root, claim in claims.items():
            ranks = _full_sector_ranks(root)
            assert ranks == {
                key: int(value) for key, value in claim["sector_ranks"].items()
            }
            total = sum(ranks.values())
            assert total == int(claim["full_system_rank"])
            assert 256 - total == int(claim["nullity"])
            assert total < 256
        return proof_detail

    _run_check("interpolation envelope and factor products", check_envelope_and_interpolation)
    _run_check("fresh exact determinant evaluations", check_fresh_exact_determinants)
    _run_check("stored full-rank controls", check_stored_full_rank_controls)
    _run_check("complementary gcd and exceptional rank drops", check_complementary_gcd_and_rank_drops)
    print("PASS")


if __name__ == "__main__":
    main()
