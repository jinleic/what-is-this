#!/usr/bin/env python3
"""Fast exact F_5-unit census and holdout classifier for branch ``11111``.

The 13 unit coordinates are written as powers of the generator 2 of F_5^x.
Every Laurent monomial is then a character of (Z/4Z)^13.  Chunked uint8
character evaluation, with int32 coefficient accumulation, enumerates the full
4^13 grid exactly while keeping memory well below the project limit.

This module is imported by e171.  Its direct entry point performs the complete
census under a hard process-time cap and prints PASS only after all points have
an exact per-point disposition.
"""

from __future__ import annotations

import argparse
import itertools
import time
from collections import Counter, defaultdict
from typing import Any

import numpy as np

import e169_kw_section_group as group_engine

UNIT_VALUES = (1, 2, 3, 4)
GENERATOR_VALUES = (1, 2, 4, 3)  # 2**log mod 5 for log=0,1,2,3.
TOTAL_POINTS = 4**13
TOTAL_SECTIONS = 4**7
SOLVE_GRID_SIZE = 4**6
DEFAULT_CAP_SECONDS = 120.0
DEFAULT_CHUNK_SIZE = 1 << 20
HOLDOUT_LOG_P = {4: 36, 6: 164, 8: 663, 10: 2280, 12: 972}


def cap_expired(deadline: float) -> bool:
    return time.process_time() >= deadline


def polynomial_character_data(polynomials) -> tuple[list[np.ndarray], list[np.ndarray]]:
    exponents: list[np.ndarray] = []
    coefficients: list[np.ndarray] = []
    for polynomial in polynomials:
        terms = polynomial.terms()
        exponents.append(
            np.array([[int(entry) % 4 for entry in monomial] for monomial, _ in terms], dtype=np.uint8)
        )
        coefficients.append(np.array([int(value) % 5 for _, value in terms], dtype=np.int16))
    return exponents, coefficients


def enumerate_unit_solutions(
    polynomials,
    deadline: float,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, Any]:
    """Enumerate the entire F_5-unit grid by exact character evaluation."""
    exponents, coefficients = polynomial_character_data(polynomials)
    polynomial_order = sorted(range(len(polynomials)), key=lambda index: len(coefficients[index]))
    powers = np.array(GENERATOR_VALUES, dtype=np.int16)
    solutions: list[int] = []
    stage_survivors = [0] * len(polynomials)
    chunks_completed = 0
    next_index = 0
    started = time.process_time()

    for start in range(0, TOTAL_POINTS, chunk_size):
        if cap_expired(deadline):
            break
        count = min(chunk_size, TOTAL_POINTS - start)
        indices = np.arange(start, start + count, dtype=np.uint32)
        logs = np.empty((count, 13), dtype=np.uint8)
        for coordinate in range(13):
            logs[:, coordinate] = (indices >> (2 * coordinate)) & 3
        live = np.arange(count, dtype=np.int32)
        chunk_complete = True
        for stage, polynomial_index in enumerate(polynomial_order):
            if cap_expired(deadline):
                chunk_complete = False
                break
            # The largest dot product is 13*3*3=117, so uint8 matmul is exact.
            phases = (logs[live] @ exponents[polynomial_index].T) & 3
            residues = (
                powers[phases] * coefficients[polynomial_index]
            ).sum(axis=1, dtype=np.int32) % 5
            live = live[residues == 0]
            stage_survivors[stage] += int(len(live))
            if not len(live):
                break
        if not chunk_complete:
            break
        solutions.extend(int(value) for value in indices[live])
        chunks_completed += 1
        next_index = start + count

    return {
        "complete": next_index == TOTAL_POINTS,
        "solution_indices": solutions,
        "next_global_index": next_index,
        "chunks_completed": chunks_completed,
        "chunk_size": chunk_size,
        "polynomial_order": polynomial_order,
        "polynomial_term_counts": [len(coefficients[index]) for index in range(len(polynomials))],
        "stage_survivor_totals": stage_survivors,
        "process_seconds": time.process_time() - started,
        "method": (
            "write each unit as 2^a, reduce monomial exponents modulo 4, and evaluate the six "
            "exact F_5 character sums in ascending term-count order on uint8 chunks"
        ),
        "integer_safety": (
            "uint8 phase dot products <=117; int32 coefficient sums <=277*4*4=4432; "
            "all reductions are exact modulo 5"
        ),
    }


def decode_point(global_index: int) -> tuple[int, ...]:
    return tuple(GENERATOR_VALUES[(global_index >> (2 * coordinate)) & 3] for coordinate in range(13))


def terms_from_polynomial(polynomial) -> list[tuple[tuple[int, ...], int]]:
    return [(tuple(map(int, monomial)), int(coefficient)) for monomial, coefficient in polynomial.terms()]


def derivative_terms(
    terms: list[tuple[tuple[int, ...], int]], coordinate: int
) -> list[tuple[tuple[int, ...], int]]:
    result: list[tuple[tuple[int, ...], int]] = []
    for exponents, coefficient in terms:
        exponent = exponents[coordinate]
        if not exponent:
            continue
        differentiated = list(exponents)
        differentiated[coordinate] -= 1
        result.append((tuple(differentiated), coefficient * exponent))
    return result


def eval_terms_mod(
    terms: list[tuple[tuple[int, ...], int]], values: tuple[int, ...] | list[int], modulus: int
) -> int:
    total = 0
    for exponents, coefficient in terms:
        term = coefficient % modulus
        for value, exponent in zip(values, exponents, strict=True):
            if exponent:
                term = term * pow(int(value) % modulus, exponent, modulus) % modulus
        total = (total + term) % modulus
    return total


def determinant_mod5(matrix: list[list[int]]) -> int:
    work = [[int(value) % 5 for value in row] for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % 5
        inverse = pow(pivot_value, -1, 5)
        for row in range(column + 1, len(work)):
            factor = work[row][column] * inverse % 5
            if factor:
                work[row] = [
                    (left - factor * right) % 5
                    for left, right in zip(work[row], work[column], strict=True)
                ]
    return determinant % 5


def rank_mod5(matrix: list[list[int]]) -> int:
    work = [[int(value) % 5 for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((candidate for candidate in range(row, len(work)) if work[candidate][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, 5)
        work[row] = [(value * inverse) % 5 for value in work[row]]
        for other in range(len(work)):
            if other != row and work[other][column]:
                factor = work[other][column]
                work[other] = [
                    (left - factor * right) % 5
                    for left, right in zip(work[other], work[row], strict=True)
                ]
        row += 1
    return row


def solve_nonsingular_mod5(matrix: list[list[int]], rhs: list[int]) -> list[int]:
    solutions = affine_solutions_mod5(matrix, rhs)
    if len(solutions) != 1:
        raise AssertionError("expected one correction for a nonsingular Jacobian")
    return solutions[0]


def affine_solutions_mod5(matrix: list[list[int]], rhs: list[int]) -> list[list[int]]:
    rows = len(matrix)
    columns = len(matrix[0])
    augmented = [
        [int(matrix[row][column]) % 5 for column in range(columns)] + [int(rhs[row]) % 5]
        for row in range(rows)
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if augmented[row][column]),
            None,
        )
        if pivot is None:
            continue
        augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
        inverse = pow(augmented[pivot_row][column], -1, 5)
        augmented[pivot_row] = [(value * inverse) % 5 for value in augmented[pivot_row]]
        for row in range(rows):
            if row != pivot_row and augmented[row][column]:
                factor = augmented[row][column]
                augmented[row] = [
                    (left - factor * right) % 5
                    for left, right in zip(augmented[row], augmented[pivot_row], strict=True)
                ]
        pivot_columns.append(column)
        pivot_row += 1
    for row in range(pivot_row, rows):
        if all(augmented[row][column] == 0 for column in range(columns)) and augmented[row][-1]:
            return []
    free_columns = [column for column in range(columns) if column not in pivot_columns]
    particular = [0] * columns
    for row, column in enumerate(pivot_columns):
        particular[column] = augmented[row][-1]
    basis: list[list[int]] = []
    for free in free_columns:
        vector = [0] * columns
        vector[free] = 1
        for row, column in enumerate(pivot_columns):
            vector[column] = -augmented[row][free] % 5
        basis.append(vector)
    result: list[list[int]] = []
    for coefficients in itertools.product(range(5), repeat=len(basis)):
        vector = particular.copy()
        for coefficient, basis_vector in zip(coefficients, basis, strict=True):
            if coefficient:
                vector = [
                    (left + coefficient * right) % 5
                    for left, right in zip(vector, basis_vector, strict=True)
                ]
        result.append(vector)
    return result


def jacobian_at_point(
    derivative_table: list[list[list[tuple[tuple[int, ...], int]]]],
    point: tuple[int, ...],
) -> list[list[int]]:
    return [
        [eval_terms_mod(derivative_table[row][column], point, 5) for column in range(6)]
        for row in range(6)
    ]


def raw_weights(e139, point: tuple[int, ...], modulus: int) -> dict[tuple[int, int], int]:
    value = dict(zip(e139.NONDIAGONAL_NAMES, point, strict=True))
    result: dict[tuple[int, int], int] = {}
    gauge = set(e139.DIRECTION_GAUGE_TREE)
    for pair in e139.ALLOWED_DIRECTION_PAIRS:
        name = group_engine.pair_name(pair)
        if pair in gauge or name in e139.DIAGONAL_NAMES or name in {"u_mx_mx", "u_my_my", "u_mz_mz"}:
            result[pair] = 1
        elif name in value:
            result[pair] = int(value[name]) % modulus
    inverse_formulas = {
        "u_mx_my": (-1, ("u_my_px",)),
        "u_mx_py": (-1, ("u_my_mx", "u_py_px")),
        "u_mx_mz": (-1, ("u_mz_px", "u_pz_mx")),
        "u_mx_pz": (-1, ("u_mz_mx", "u_pz_px")),
        "u_my_mz": (-1, ("u_mz_py", "u_py_pz", "u_pz_my")),
        "u_my_pz": (-1, ("u_mz_my", "u_py_mz", "u_pz_py")),
    }
    for pair in e139.ALLOWED_DIRECTION_PAIRS:
        name = group_engine.pair_name(pair)
        if name not in inverse_formulas:
            continue
        numerator, denominator_names = inverse_formulas[name]
        denominator = 1
        for denominator_name in denominator_names:
            denominator = denominator * int(value[denominator_name]) % modulus
        result[pair] = numerator * pow(denominator, -1, modulus) % modulus
    if set(result) != set(e139.ALLOWED_DIRECTION_PAIRS):
        raise AssertionError("numeric thin chart did not reconstruct all weights")
    return result


def point_values25(e139, construction: dict[str, Any], point: tuple[int, ...], modulus: int) -> list[int]:
    weights = raw_weights(e139, point, modulus)
    by_name = {group_engine.pair_name(pair): value for pair, value in weights.items()}
    return [by_name[str(variable)] for variable in construction["variables"]]


def trace_power_mod(
    e139,
    structure: dict[str, Any],
    point: tuple[int, ...],
    order: int,
    modulus: int,
) -> int:
    weights = raw_weights(e139, point, modulus)
    size = int(structure["n"])
    matrix = np.zeros((size, size), dtype=np.int64)
    directions = structure["directions"]
    for edge, next_edges in enumerate(structure["transitions"]):
        for next_edge in next_edges:
            matrix[edge, next_edge] = weights[(directions[edge], directions[next_edge])]
    result = np.eye(size, dtype=np.int64)
    base = matrix
    exponent = order
    while exponent:
        if exponent & 1:
            result = (result @ base) % modulus
        exponent //= 2
        if exponent:
            base = (base @ base) % modulus
    return int(np.trace(result) % modulus)


def holdout_residue(
    e139,
    structure: dict[str, Any],
    point: tuple[int, ...],
    order: int,
    modulus: int,
) -> int:
    return (trace_power_mod(e139, structure, point, order, modulus) + 2 * order * HOLDOUT_LOG_P[order]) % modulus


def all_construction_zero(
    e139,
    construction: dict[str, Any],
    point: tuple[int, ...],
    modulus: int,
) -> bool:
    values25 = point_values25(e139, construction, point, modulus)
    return all(
        e139.eval_sparse_mod(polynomial, values25, modulus) == 0
        for polynomial in construction["sparse42"]
    )


def hensel_lift_nonsingular(
    primitive_terms: list[list[tuple[tuple[int, ...], int]]],
    point: tuple[int, ...],
    jacobian: list[list[int]],
) -> tuple[int, ...]:
    solve = list(point[:6])
    held = tuple(point[6:])
    modulus = 5
    for _exponent in range(2, 5):
        next_modulus = modulus * 5
        residues = [eval_terms_mod(terms, tuple(solve) + held, next_modulus) for terms in primitive_terms]
        if any(residue % modulus for residue in residues):
            raise AssertionError("Hensel input does not solve the previous modulus")
        rhs = [-(residue // modulus) % 5 for residue in residues]
        correction = solve_nonsingular_mod5(jacobian, rhs)
        solve = [
            (value + modulus * digit) % next_modulus
            for value, digit in zip(solve, correction, strict=True)
        ]
        modulus = next_modulus
    lifted = tuple(solve) + held
    if any(eval_terms_mod(terms, lifted, 625) for terms in primitive_terms):
        raise AssertionError("nonsingular Hensel lift failed")
    return lifted


def singular_obstruction_tree(
    e139,
    construction: dict[str, Any],
    structure: dict[str, Any],
    primitive_terms: list[list[tuple[tuple[int, ...], int]]],
    point: tuple[int, ...],
    jacobian: list[list[int]],
    deadline: float,
) -> dict[str, Any]:
    current = [tuple(point[:6])]
    held = tuple(point[6:])
    modulus = 5
    levels: list[dict[str, Any]] = []
    for _exponent in range(2, 5):
        if cap_expired(deadline):
            return {"complete": False, "levels": levels, "survivors": current, "last_modulus": modulus}
        next_modulus = modulus * 5
        next_points: list[tuple[int, ...]] = []
        branch_records: list[dict[str, Any]] = []
        for solve in current:
            if cap_expired(deadline):
                return {"complete": False, "levels": levels, "survivors": current, "last_modulus": modulus}
            residues = [eval_terms_mod(terms, solve + held, next_modulus) for terms in primitive_terms]
            if any(residue % modulus for residue in residues):
                raise AssertionError("singular lift input failed its previous modulus")
            rhs = [-(residue // modulus) % 5 for residue in residues]
            corrections = affine_solutions_mod5(jacobian, rhs)
            candidates: list[dict[str, Any]] = []
            for correction in corrections:
                lifted_solve = tuple(
                    (value + modulus * digit) % next_modulus
                    for value, digit in zip(solve, correction, strict=True)
                )
                lifted_point = lifted_solve + held
                primitive_zero = all(
                    eval_terms_mod(terms, lifted_point, next_modulus) == 0
                    for terms in primitive_terms
                )
                full42_zero = all_construction_zero(e139, construction, lifted_point, next_modulus)
                if not primitive_zero or not full42_zero:
                    raise AssertionError("linearized singular lift did not satisfy the full construction system")
                residue8 = holdout_residue(e139, structure, lifted_point, 8, next_modulus)
                candidates.append({
                    "correction_mod5": correction,
                    "solve6": list(lifted_solve),
                    "construction_residues_all_zero": True,
                    "holdout_order8_residue": int(residue8),
                })
                if residue8 == 0:
                    next_points.append(lifted_solve)
            branch_records.append({
                "input_solve6": list(solve),
                "linear_rhs_mod5": rhs,
                "correction_count": len(corrections),
                "candidates": candidates,
            })
        levels.append({
            "modulus": next_modulus,
            "input_count": len(current),
            "construction_lift_count": sum(row["correction_count"] for row in branch_records),
            "holdout_zero_lift_count": len(next_points),
            "branches": branch_records,
        })
        current = next_points
        modulus = next_modulus
        if not current:
            break
    return {"complete": True, "levels": levels, "survivors": [list(value) for value in current], "last_modulus": modulus}


def classify_solutions(
    e139,
    construction: dict[str, Any],
    solution_indices: list[int],
    deadline: float,
) -> dict[str, Any]:
    started = time.process_time()
    points = sorted(decode_point(index) for index in solution_indices)
    primitive_terms = [terms_from_polynomial(polynomial) for polynomial in construction["primitives"]]
    derivative_table = [
        [derivative_terms(primitive_terms[row], column) for column in range(6)]
        for row in range(6)
    ]
    structure = e139.transfer_matrix_structure(e139.HOLDOUT_SHAPE)
    preliminary: list[dict[str, Any]] = []

    for point in points:
        if cap_expired(deadline):
            return {"complete": False, "point_records": preliminary, "process_seconds": time.process_time() - started}
        residue8 = holdout_residue(e139, structure, point, 8, 5)
        if residue8:
            preliminary.append({
                "point13_mod5": list(point),
                "disposition": "KILLED_H8_MOD5",
                "holdout_residues": {"8_mod_5": int(residue8)},
                "certificate_modulus": 5,
            })
            continue
        residue12 = holdout_residue(e139, structure, point, 12, 5)
        if residue12:
            preliminary.append({
                "point13_mod5": list(point),
                "disposition": "KILLED_H12_MOD5",
                "holdout_residues": {"8_mod_5": 0, "12_mod_5": int(residue12)},
                "certificate_modulus": 5,
            })
            continue
        residue10 = holdout_residue(e139, structure, point, 10, 5)
        jacobian = jacobian_at_point(derivative_table, point)
        determinant = determinant_mod5(jacobian)
        preliminary.append({
            "point13_mod5": list(point),
            "disposition": "PENDING_LIFT",
            "holdout_residues": {"8_mod_5": 0, "10_mod_5": int(residue10), "12_mod_5": 0},
            "jacobian_mod5": jacobian,
            "jacobian_determinant_mod5": determinant,
            "jacobian_rank_mod5": rank_mod5(jacobian),
        })

    final: list[dict[str, Any]] = []
    for record in preliminary:
        if record["disposition"] != "PENDING_LIFT":
            final.append(record)
            continue
        if cap_expired(deadline):
            final.append(record)
            continue
        point = tuple(record["point13_mod5"])
        determinant = int(record["jacobian_determinant_mod5"])
        if determinant:
            lift = hensel_lift_nonsingular(primitive_terms, point, record["jacobian_mod5"])
            if not all_construction_zero(e139, construction, lift, 625):
                raise AssertionError("nonsingular lift did not satisfy all 42 construction equations")
            residue8 = holdout_residue(e139, structure, lift, 8, 625)
            record["lift_solve6_mod625"] = list(lift[:6])
            record["full42_construction_residues_mod625_all_zero"] = True
            record["holdout_residues"]["8_mod_625"] = int(residue8)
            record["certificate_modulus"] = 625
            record["disposition"] = (
                "KILLED_H8_MOD625_NONSINGULAR" if residue8 else "UNDECIDED_NONSINGULAR_H8_SURVIVOR"
            )
        else:
            tree = singular_obstruction_tree(
                e139,
                construction,
                structure,
                primitive_terms,
                point,
                record["jacobian_mod5"],
                deadline,
            )
            record["singular_lift_tree"] = tree
            if tree["complete"] and not tree["survivors"]:
                record["disposition"] = "KILLED_BY_SINGULAR_LIFT_OBSTRUCTION"
                record["certificate_modulus"] = int(tree["last_modulus"])
            else:
                record["disposition"] = "UNDECIDED_SINGULAR_LIFT_SURVIVOR"
                record["certificate_modulus"] = int(tree["last_modulus"])
        final.append(record)

    killed_prefix = "KILLED_"
    complete = len(final) == len(points) and all(
        record["disposition"].startswith(killed_prefix) for record in final
    )
    for index, record in enumerate(final):
        record["point_id"] = f"P{index:04d}"
    return {
        "complete": complete,
        "point_records": final,
        "process_seconds": time.process_time() - started,
        "disposition_counts": dict(sorted(Counter(record["disposition"] for record in final).items())),
    }


def lex_section_index(section: tuple[int, ...]) -> int:
    index = 0
    for value in section:
        index = 4 * index + UNIT_VALUES.index(value)
    return index


def make_status_bitmap(statuses: list[bool]) -> str:
    packed = bytearray((len(statuses) + 7) // 8)
    for index, status in enumerate(statuses):
        if status:
            packed[index // 8] |= 1 << (index % 8)
    return packed.hex()


def classify_sections(point_records: list[dict[str, Any]]) -> dict[str, Any]:
    by_section: dict[tuple[int, ...], list[dict[str, Any]]] = defaultdict(list)
    for record in point_records:
        point = tuple(record["point13_mod5"])
        by_section[point[6:]].append(record)
    rows: list[dict[str, Any]] = []
    status_bits: list[bool] = []
    undecided_indices: list[int] = []
    for index, section in enumerate(itertools.product(UNIT_VALUES, repeat=7)):
        records = by_section.get(section, [])
        if records and all(record["disposition"].startswith("KILLED_") for record in records):
            status = "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT"
            reason = (
                "Every F_5-unit construction point in this section has a recorded exact holdout "
                "residue or finite lift obstruction."
            )
            bit = True
        else:
            status = "UNDECIDED"
            if not records:
                reason = (
                    "The exhaustive F_5-unit solve grid has no construction point.  This excludes "
                    "only integral 5-adic unit points; it is not an EMPTY_OVER_Q certificate."
                )
            else:
                reason = "At least one recorded local construction point remains undecided."
            bit = False
            undecided_indices.append(index)
        status_bits.append(bit)
        rows.append({
            "lex_index": index,
            "section": list(section),
            "orbit_size": 1,
            "status": status,
            "solution_count_mod5": len(records),
            "point_ids": [record["point_id"] for record in records],
            "certificate": {
                "method": "complete_F5_unit_solve_grid" if not records else "per_point_exact_holdout_or_lift",
                "detail": reason,
            },
        })
    counts = Counter(row["status"] for row in rows)
    counts["EMPTY_OVER_Q"] = 0
    return {
        "classification_rows": rows,
        "status_counts": {
            "EMPTY_OVER_Q": 0,
            "HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT": counts["HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT"],
            "UNDECIDED": counts["UNDECIDED"],
        },
        "nonempty_section_count": len(by_section),
        "undecided_lex_indices": undecided_indices,
        "status_bitmap": {
            "encoding": (
                "2048 bytes in hex; lex_index i is HAS_LOCAL_POINT_BUT_FAILS_HOLDOUT iff "
                "bit (i mod 8) of byte floor(i/8) is 1; zero means UNDECIDED"
            ),
            "hex": make_status_bitmap(status_bits),
        },
    }


def run_complete_census(
    e139,
    construction: dict[str, Any],
    deadline: float,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, Any]:
    enumeration = enumerate_unit_solutions(construction["primitives"], deadline, chunk_size)
    if not enumeration["complete"]:
        return {"complete": False, "enumeration": enumeration, "holdout": None, "sections": None}
    holdout = classify_solutions(e139, construction, enumeration["solution_indices"], deadline)
    sections = classify_sections(holdout["point_records"])
    complete = holdout["complete"] and len(sections["classification_rows"]) == TOTAL_SECTIONS
    return {"complete": complete, "enumeration": enumeration, "holdout": holdout, "sections": sections}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap-seconds", type=float, default=DEFAULT_CAP_SECONDS)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    args = parser.parse_args()
    if args.cap_seconds <= 0 or args.chunk_size <= 0:
        raise ValueError("cap and chunk size must be positive")
    started = time.process_time()
    deadline = started + args.cap_seconds
    e139 = group_engine.load_e139()
    construction = group_engine.build_construction(e139)
    result = run_complete_census(e139, construction, deadline, args.chunk_size)
    if not result["complete"]:
        print(
            "PARTIAL cap reached; next global index=",
            result["enumeration"]["next_global_index"],
        )
        return
    print(
        f"solutions={len(result['holdout']['point_records'])}; "
        f"nonempty sections={result['sections']['nonempty_section_count']}; "
        f"statuses={result['sections']['status_counts']}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
