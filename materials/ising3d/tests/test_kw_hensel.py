"""Standalone exact checks for the Kac--Ward p-adic lifting experiment."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import sympy as sp

from ising.fermions.kac_ward import DIRECTION_GAUGE_TREE, full_weight_finite_system

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/kac_ward/full_family.json"
RESULT = ROOT / "results/kac_ward/hensel_lift.json"
SHAPES = ((3, 3, 2), (2, 2, 3))
ORDERS = (4, 6, 8)


def rref_pivots_mod(matrix: list[list[int]], prime: int) -> tuple[int, list[int]]:
    rows = [[value % prime for value in row] for row in matrix]
    rank = 0
    pivots = []
    for column in range(len(rows[0])):
        pivot = next(
            (row for row in range(rank, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(value * inverse) % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (left - factor * right) % prime
                for left, right in zip(rows[row], rows[rank], strict=True)
            ]
        pivots.append(column)
        rank += 1
        if rank == len(rows):
            break
    return rank, pivots


def solve_mod_prime(
    matrix: list[list[int]], right_hand_side: list[int], prime: int
) -> list[int]:
    size = len(matrix)
    rows = [
        [value % prime for value in row] + [right_hand_side[index] % prime]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = pow(rows[column][column], -1, prime)
        rows[column] = [(value * inverse) % prime for value in rows[column]]
        for row in range(size):
            if row == column or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (left - factor * right) % prime
                for left, right in zip(rows[row], rows[column], strict=True)
            ]
    return [rows[index][-1] for index in range(size)]


def rational_reconstruct_independent(
    residue: int, modulus: int, bound: int
) -> Fraction | None:
    residue %= modulus
    if residue == 0:
        return Fraction(0, 1)
    old_r, r = modulus, residue
    old_t, t = 0, 1
    while r > bound:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_t, t = t, old_t - quotient * t
    numerator, denominator = r, t
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if (
        abs(numerator) > bound
        or not 1 <= denominator <= bound
        or math.gcd(numerator, denominator) != 1
        or math.gcd(denominator, modulus) != 1
        or (numerator - residue * denominator) % modulus
    ):
        return None
    return Fraction(numerator, denominator)


def load_witnesses() -> dict[int, list[int]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    return {
        int(row["prime"]): [int(value) for value in row["values_in_variable_order"]]
        for row in source["data"]["exact_modular_solution_witnesses"]
    }


def construction_system():
    system = full_weight_finite_system(SHAPES, ORDERS, gauge_fix=True)
    assert len(system.variables) == 25
    assert len(system.equations) == 6
    return system


def check_one_hensel_step() -> None:
    system = construction_system()
    witness = load_witnesses()[5]
    substitutions = dict(zip(system.variables, witness, strict=True))
    exact_residuals = [int(equation.subs(substitutions)) for equation in system.equations]
    assert not any(residual % 5 for residual in exact_residuals)

    jacobian = sp.Matrix(system.equations).jacobian(system.variables)
    rows = [
        [int(entry.subs(substitutions)) % 5 for entry in row]
        for row in jacobian.tolist()
    ]
    rank, pivots = rref_pivots_mod(rows, 5)
    assert rank == 6
    pivots = pivots[:6]
    square = [[row[column] for column in pivots] for row in rows]
    right_hand_side = [-(residual // 5) % 5 for residual in exact_residuals]
    correction = solve_mod_prime(square, right_hand_side, 5)
    lifted = witness.copy()
    for column, value in zip(pivots, correction, strict=True):
        lifted[column] += 5 * value
    lifted_substitutions = dict(zip(system.variables, lifted, strict=True))
    lifted_residuals = [
        int(equation.subs(lifted_substitutions)) % 25
        for equation in system.equations
    ]
    assert lifted_residuals == [0] * 6


def check_F7_nonlift_certificate() -> None:
    system = construction_system()
    witness = load_witnesses()[7]
    substitutions = dict(zip(system.variables, witness, strict=True))
    residuals = [int(equation.subs(substitutions)) for equation in system.equations]
    right_hand_side = [-(residual // 7) % 7 for residual in residuals]
    jacobian = sp.Matrix(system.equations).jacobian(system.variables)
    rows = [
        [int(entry.subs(substitutions)) % 7 for entry in row]
        for row in jacobian.tolist()
    ]
    rank, _ = rref_pivots_mod(rows, 7)
    assert rank == 5
    left_kernel = [4, 0, 0, 1, 0, 0]
    assert all(
        sum(left_kernel[row] * rows[row][column] for row in range(6)) % 7 == 0
        for column in range(25)
    )
    assert sum(
        left_kernel[row] * right_hand_side[row] for row in range(6)
    ) % 7 == 5


def check_reconstruction_bookkeeping() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    lifts = {row["prime"]: row for row in payload["data"]["lifts"]}
    assert {prime for prime, row in lifts.items() if row["status"] == "lifted"} == {5, 11}
    for prime in (5, 11):
        row = lifts[prime]
        modulus = int(row["modulus"])
        exponent = int(row["precision_exponent"])
        assert modulus == prime**exponent
        assert len(str(modulus)) - 1 >= 60
        bound = int(row["rational_reconstruction"]["height_bound"])
        assert 2 * bound * bound < modulus
        assert 2 * (bound + 1) * (bound + 1) >= modulus
        reconstructed = []
        for coordinate, residue in zip(
            row["rational_reconstruction"]["coordinates"],
            row["residues"],
            strict=True,
        ):
            value = rational_reconstruct_independent(int(residue), modulus, bound)
            reconstructed.append(value)
            if coordinate["status"] == "failure":
                assert value is None
            else:
                assert value == Fraction(
                    int(coordinate["numerator"]), int(coordinate["denominator"])
                )
        assert sum(value is not None for value in reconstructed) == row[
            "rational_reconstruction"
        ]["successful_coordinate_count"]
        assert not all(value is not None for value in reconstructed)

    modulus = 5**12
    bound = math.isqrt((modulus - 1) // 2)
    target = Fraction(3, 4)
    residue = target.numerator * pow(target.denominator, -1, modulus) % modulus
    assert rational_reconstruct_independent(residue, modulus, bound) == target


def check_exact_branch_certificate() -> None:
    ungauged = full_weight_finite_system(SHAPES, ORDERS, gauge_fix=False)
    symbol_map = dict(ungauged.all_symbols)
    substitutions = {symbol_map[pair]: 0 for pair in DIRECTION_GAUGE_TREE}
    equations = [sp.expand(equation.subs(substitutions)) for equation in ungauged.equations]
    assert sp.expand(2 * equations[0] - 3 * equations[3]) == 56

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    certificate = payload["data"]["exact_smallest_branch_certificate"]
    assert certificate["branch"] == "00000"
    assert certificate["combination_coefficients"] == [2, 0, 0, -3, 0, 0]
    assert certificate["constant"] == 56
    assert certificate["identity_verified"]
    attempts = payload["data"]["groebner_retries"]
    assert {row["prime"] for row in attempts} == {101, 32003}
    assert all(row["status"] == "completed" and row["contains_one"] for row in attempts)


def check_result_scope() -> None:
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    assert payload["data"]["headline"].startswith("UNRESOLVED")
    assert payload["data"]["exact_reconstruction_candidates"] == []
    assert all(check["passed"] for check in payload["checks"])
    cross = payload["data"]["cross_prime_reconstruction"]
    assert {
        row["variable"] for row in cross["consistently_reconstructed_coordinates"]
    } == {"u_pz_mx", "u_pz_my"}
    held = {
        name
        for row in payload["data"]["lifts"]
        if row["status"] == "lifted"
        for name in row["local_section"]["held_coordinates"]
    }
    assert all(
        row["variable"] in held
        for row in cross["consistently_reconstructed_coordinates"]
    )


def run_check(name: str, function) -> None:
    try:
        function()
    except Exception as error:
        print(f"FAIL {name}: {type(error).__name__}: {error}")
        raise
    print(f"PASS {name}")


def main() -> None:
    run_check("one Hensel step modulo 25", check_one_hensel_step)
    run_check("F_7 point nonlift certificate", check_F7_nonlift_certificate)
    run_check("height bounds and rational reconstruction", check_reconstruction_bookkeeping)
    run_check("exact branch-00000 constant certificate", check_exact_branch_certificate)
    run_check("artifact scope", check_result_scope)
    print("PASS")


if __name__ == "__main__":
    main()
