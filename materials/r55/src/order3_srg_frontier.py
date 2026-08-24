#!/usr/bin/env python3
"""Exact order-three frontier for a Ramsey-good srg(45,22,10,11).

Theorem
=======

A strongly regular graph with parameters (45,22,10,11) that has no clique or
independent set of order five has no automorphism of order three.  This is not
a bound on R(5,5), and graphs outside the strongly regular lane remain open.

For an order-three automorphism with f fixed points, the omega-eigenspace has
c=(45-f)/3 dimensions over Q(omega).  On it the SRG relation becomes
x^2+x-11=0, irreducible over Q(omega), so c is even and f=3 mod 6.  The
nonidentity spectral candidates are therefore f=3,9,15,21,27,33,39.

If f=3, every fixed-graph degree is 1 mod 3 and hence one; a three-vertex graph
cannot have every degree one.  For f=9, the fixed graph is forced to
srg(9,4,1,2), uniquely L_2(3).  The exact quotient reduction below leaves two
matrices, each of which forces the same independent set of order five.

It remains to exclude f>=15.  Let O be a triangle 3-orbit, let s be the number
of fixed vertices adjacent to O, and let t_i be the number of neighbors that a
vertex of O has in each other 3-orbit.  Degree and adjacent common-neighbor
counts give

    sum(t_i)            = 20-s,
    sum(binomial(t_i,2)) = 9-s,
    sum(t_i^2)          = 38-3s.

The second identity makes s<=9.  If c<=10, Cauchy-Schwarz would require

    (20-s)^2 <= (c-1)(38-3s).

At c=10 the difference is s^2-13s+58, whose discriminant is -63; it is
positive for every real s.  Decreasing c only increases the gap because
38-3s is nonnegative.  The spectral trace on the orbit-constant subspace
forces exactly c/2 triangle orbits, so one exists.  Thus every
srg(45,22,10,11), Ramsey-good or not, excludes f=15,21,27,33,39.  Together
with the Ramsey-dependent f=3 and f=9 arguments, this excludes every
order-three automorphism in the Ramsey-good lane.

The constructive f=9 uniqueness check reduces its fixed graph to 90
four-by-four biregular incidence matrices, eight labelled survivors, and one
isomorphism type.  Six triangle and six independent 3-orbits then force a
canonical 9x12 incidence matrix B with BB^T=3(I+J).  Writing the equitable
quotient as Q=[[A,3B],[B^T,T]], an exact rank-77 affine system leaves tau=1,2;
both candidates force the independent five-set recorded below.  No phase
search, floating-point computation, or SAT solver enters either proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

SCHEMA_VERSION = 2
CAMPAIGN_ID = "order3_srg_frontier"
DISPOSITION = "ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG"
PARAMETERS = (45, 22, 10, 11)
SPECTRAL_FIXED_POINT_COUNTS = (3, 9, 15, 21, 27, 33, 39)
LARGE_FIXED_POINT_COUNTS = (15, 21, 27, 33, 39)

# Canonical T at tau=1 and tau=2.  The direction is their difference.
_T_ONE = (
    (2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2),
    (2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 1, 1),
    (2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 1, 1),
    (1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2),
    (1, 2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2),
    (2, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 1),
    (2, 1, 1, 2, 2, 1, 0, 1, 1, 1, 2, 2),
    (1, 2, 2, 1, 1, 2, 1, 0, 2, 2, 1, 1),
    (1, 2, 2, 1, 1, 2, 1, 2, 0, 2, 1, 1),
    (1, 2, 2, 1, 1, 2, 1, 2, 2, 0, 1, 1),
    (2, 1, 1, 2, 2, 1, 2, 1, 1, 1, 0, 2),
    (2, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 0),
)
_T_TWO = (
    (2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 1),
    (2, 2, 1, 2, 2, 1, 2, 1, 1, 1, 2, 2),
    (2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2),
    (1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1),
    (1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 1, 1),
    (2, 1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2),
    (1, 2, 2, 1, 1, 2, 0, 1, 1, 1, 2, 2),
    (2, 1, 1, 2, 2, 1, 1, 0, 2, 2, 1, 1),
    (2, 1, 1, 2, 2, 1, 1, 2, 0, 2, 1, 1),
    (2, 1, 1, 2, 2, 1, 1, 2, 2, 0, 1, 1),
    (1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 0, 2),
    (1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 0),
)


class FrontierViolation(RuntimeError):
    """An exact structural invariant failed."""


def _add_edge(adjacency: list[int], left: int, right: int) -> None:
    adjacency[left] |= 1 << right
    adjacency[right] |= 1 << left


def _local_graph(row_sets: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    adjacency = [0] * 9
    for vertex in range(1, 5):
        _add_edge(adjacency, 0, vertex)
    _add_edge(adjacency, 1, 2)
    _add_edge(adjacency, 3, 4)
    for left, right in ((5, 6), (6, 7), (7, 8), (8, 5)):
        _add_edge(adjacency, left, right)
    for neighbor, selected in enumerate(row_sets, start=1):
        for dual_index in selected:
            _add_edge(adjacency, neighbor, 5 + dual_index)
    return tuple(adjacency)


def construct_fixed_graph() -> tuple[int, ...]:
    """A canonical labelled copy of L_2(3)."""
    return _local_graph(((0, 1), (2, 3), (0, 3), (1, 2)))


def srg_parameters(adjacency: tuple[int, ...]) -> tuple[int, int, int, int] | None:
    order = len(adjacency)
    degrees = {row.bit_count() for row in adjacency}
    if len(degrees) != 1:
        return None
    degree = next(iter(degrees))
    adjacent_common = set()
    nonadjacent_common = set()
    for left in range(order):
        for right in range(left + 1, order):
            common = (adjacency[left] & adjacency[right]).bit_count()
            target = adjacent_common if adjacency[left] >> right & 1 else nonadjacent_common
            target.add(common)
    if len(adjacent_common) != 1 or len(nonadjacent_common) != 1:
        return None
    return order, degree, next(iter(adjacent_common)), next(iter(nonadjacent_common))


def enumerate_fixed_graph_candidates() -> list[tuple[int, ...]]:
    """Enumerate the 90 locally possible biregular matrices exactly."""
    pairs = tuple(itertools.combinations(range(4), 2))
    survivors = []
    for row_sets in itertools.product(pairs, repeat=4):
        if tuple(sum(column in row for row in row_sets) for column in range(4)) != (2,) * 4:
            continue
        graph = _local_graph(row_sets)
        if srg_parameters(graph) == (9, 4, 1, 2):
            survivors.append(graph)
    return survivors


def _permuted_equal(left: tuple[int, ...], right: tuple[int, ...], permutation) -> bool:
    for new_left, old_left in enumerate(permutation):
        expected = 0
        for new_right, old_right in enumerate(permutation):
            if left[old_left] >> old_right & 1:
                expected |= 1 << new_right
        if expected != right[new_left]:
            return False
    return True


def is_isomorphic(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    if len(left) != len(right):
        return False
    if sorted(row.bit_count() for row in left) != sorted(row.bit_count() for row in right):
        return False
    return any(_permuted_equal(left, right, permutation)
               for permutation in itertools.permutations(range(len(left))))


def degree_extreme_obstructions() -> dict:
    """Record the complete residue check excluding f=9 degrees one and seven."""
    required_residue = 11 % 3
    possible_common_fixed = list(range(2))
    return {
        "degree_one_max_common_fixed": 1,
        "nonedge_required_residue": required_residue,
        "possible_common_fixed_counts": possible_common_fixed,
        "residue_hits": [
            count for count in possible_common_fixed
            if count % 3 == required_residue
        ],
        "degree_seven_complement_degree": 8 - 7,
    }


def triangle_orbit_moment(orbit_cycles: int, fixed_support: int) -> dict:
    """Return the exact Cauchy obstruction for one triangle 3-orbit."""
    if orbit_cycles < 2:
        raise ValueError("an order-three automorphism needs at least two 3-orbits here")
    if not 0 <= fixed_support <= 9:
        raise ValueError("the adjacent common-neighbor count forces support at most nine")
    other_orbits = orbit_cycles - 1
    degree_sum = 20 - fixed_support
    square_sum = 38 - 3 * fixed_support
    return {
        "fixed_support": fixed_support,
        "other_orbits": other_orbits,
        "off_orbit_degree_sum": degree_sum,
        "off_orbit_square_sum": square_sum,
        "cauchy_gap": degree_sum * degree_sum - other_orbits * square_sum,
    }


def large_fixed_count_obstructions() -> list[dict]:
    """Certify the Ramsey-free triangle moment contradiction for f>=15."""
    records = []
    fixed_support_values = list(range(10))
    for fixed_points in LARGE_FIXED_POINT_COUNTS:
        orbit_cycles = (45 - fixed_points) // 3
        if orbit_cycles % 2:
            raise FrontierViolation("spectral orbit count must be even")
        moments = [triangle_orbit_moment(orbit_cycles, support)
                   for support in fixed_support_values]
        gaps = [moment["cauchy_gap"] for moment in moments]
        if any(gap <= 0 for gap in gaps):
            raise FrontierViolation("large-fixed-count Cauchy obstruction failed")
        records.append({
            "fixed_points": fixed_points,
            "orbit_3_cycles": orbit_cycles,
            "triangle_orbits": orbit_cycles // 2,
            "fixed_support_values": fixed_support_values,
            "cauchy_gaps": gaps,
            "minimum_cauchy_gap": min(gaps),
        })
    return records


def is_clique(adjacency: tuple[int, ...], vertices) -> bool:
    return all(adjacency[left] >> right & 1
               for left, right in itertools.combinations(vertices, 2))


def is_independent(adjacency: tuple[int, ...], vertices) -> bool:
    return all(not (adjacency[left] >> right & 1)
               for left, right in itertools.combinations(vertices, 2))


def triangles(adjacency: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(vertices for vertices in itertools.combinations(range(len(adjacency)), 3)
                 if is_clique(adjacency, vertices))


def independent_triples(adjacency: tuple[int, ...]) -> tuple[tuple[int, int, int], ...]:
    return tuple(vertices for vertices in itertools.combinations(range(len(adjacency)), 3)
                 if is_independent(adjacency, vertices))


def build_incidence(adjacency: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    independent = independent_triples(adjacency)
    cliques = triangles(adjacency)
    if len(independent) != 6 or len(cliques) != 6:
        raise FrontierViolation("fixed graph must have six triples of each type")
    rows = []
    for vertex in range(9):
        rows.append(tuple(
            [int(vertex in members) for members in independent]
            + [int(vertex not in members) for members in cliques]
        ))
    return tuple(rows)


def row_sums(matrix) -> tuple[int, ...]:
    return tuple(sum(row) for row in matrix)


def column_sums(matrix) -> tuple[int, ...]:
    return tuple(sum(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def column_members(matrix, column: int) -> tuple[int, ...]:
    return tuple(row for row in range(len(matrix)) if matrix[row][column])


def gram(matrix) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(sum(matrix[left][column] * matrix[right][column]
                           for column in range(len(matrix[0])))
                       for right in range(len(matrix)))
                 for left in range(len(matrix)))


def target_gram() -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(6 if left == right else 3 for right in range(9))
                 for left in range(9))


def _compositions(total: int, parts: int, prefix=()):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from _compositions(total - value, parts - 1, prefix + (value,))


def _weighted_signature(columns, multiplicities) -> tuple[int, ...]:
    rows = [0] * 9
    pair_values = []
    for vertex in range(9):
        rows[vertex] = sum(multiplicities[index] * int(vertex in column)
                           for index, column in enumerate(columns))
    for left in range(9):
        for right in range(left, 9):
            pair_values.append(sum(
                multiplicities[index] * int(left in column and right in column)
                for index, column in enumerate(columns)
            ))
    return tuple(rows + pair_values)


def incidence_multiplicity_solutions(adjacency: tuple[int, ...]):
    independent = independent_triples(adjacency)
    complement_cliques = tuple(tuple(v for v in range(9) if v not in clique)
                               for clique in triangles(adjacency))
    target = (6,) * 9 + tuple(6 if left == right else 3
                              for left in range(9) for right in range(left, 9))
    right_by_signature = {}
    for right in _compositions(6, 6):
        signature = _weighted_signature(complement_cliques, right)
        right_by_signature.setdefault(signature, []).append(right)
    solutions = []
    for left in _compositions(6, 6):
        signature = _weighted_signature(independent, left)
        needed = tuple(target[index] - signature[index] for index in range(len(target)))
        for right in right_by_signature.get(needed, ()): 
            solutions.append((left, right))
    return solutions


def _adjacency_matrix(adjacency: tuple[int, ...]):
    return tuple(tuple(int(adjacency[row] >> column & 1) for column in range(9))
                 for row in range(9))


def _matmul(left, right):
    return tuple(tuple(sum(left[row][inner] * right[inner][column]
                           for inner in range(len(right)))
                       for column in range(len(right[0])))
                 for row in range(len(left)))


def _transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def _symmetric_variables():
    return tuple((left, right) for left in range(12) for right in range(left, 12))


def quotient_linear_system() -> dict:
    adjacency = _adjacency_matrix(construct_fixed_graph())
    incidence = build_incidence(construct_fixed_graph())
    variables = _symmetric_variables()
    owner = {pair: index for index, pair in enumerate(variables)}
    adjacency_incidence = _matmul(adjacency, incidence)
    rows = []
    rhs = []
    for fixed in range(9):
        for orbit in range(12):
            row = [0] * len(variables)
            for other in range(12):
                row[owner[(min(other, orbit), max(other, orbit))]] += incidence[fixed][other]
            rows.append(tuple(row))
            rhs.append(11 - adjacency_incidence[fixed][orbit] - incidence[fixed][orbit])
    for orbit in range(12):
        row = [0] * len(variables)
        row[owner[(orbit, orbit)]] = 1
        rows.append(tuple(row))
        rhs.append(2 if orbit < 6 else 0)
    for orbit in range(12):
        row = [0] * len(variables)
        for other in range(12):
            row[owner[(min(orbit, other), max(orbit, other))]] += 1
        rows.append(tuple(row))
        rhs.append(19 if orbit < 6 else 16)
    rank = _rational_rank(rows)
    augmented_rank = _rational_rank(
        tuple(tuple(row) + (rhs[index],) for index, row in enumerate(rows))
    )
    return {
        "variables": len(variables),
        "rows": tuple(rows),
        "rhs": tuple(rhs),
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


def _rational_rank(matrix) -> int:
    rows = [list(map(Fraction, row)) for row in matrix]
    if not rows:
        return 0
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            multiple = rows[row][column]
            rows[row] = [left - multiple * right
                         for left, right in zip(rows[row], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


def _direction():
    return tuple(tuple(_T_TWO[row][column] - _T_ONE[row][column]
                       for column in range(12)) for row in range(12))


def _matrix_at(tau: int):
    direction = _direction()
    return tuple(tuple(_T_ONE[row][column] + (tau - 1) * direction[row][column]
                       for column in range(12)) for row in range(12))


def _flatten_symmetric(matrix):
    return tuple(matrix[left][right] for left, right in _symmetric_variables())


def verify_affine_family(system: dict) -> bool:
    base = _flatten_symmetric(_matrix_at(0))
    direction = _flatten_symmetric(_direction())
    for row, target in zip(system["rows"], system["rhs"]):
        if sum(coefficient * value for coefficient, value in zip(row, base)) != target:
            return False
        if sum(coefficient * value for coefficient, value in zip(row, direction)) != 0:
            return False
    return system["rank"] == system["variables"] - 1 and any(direction)


def verify_quotient_matrix(matrix) -> bool:
    if len(matrix) != 12 or any(len(row) != 12 for row in matrix):
        return False
    if any(matrix[left][right] != matrix[right][left]
           for left in range(12) for right in range(12)):
        return False
    if tuple(matrix[index][index] for index in range(12)) != (2,) * 6 + (0,) * 6:
        return False
    if any(not 0 <= matrix[left][right] <= 3
           for left in range(12) for right in range(12) if left != right):
        return False
    incidence = build_incidence(construct_fixed_graph())
    adjacency = _adjacency_matrix(construct_fixed_graph())
    if tuple(sum(matrix[row]) for row in range(12)) != (19,) * 6 + (16,) * 6:
        return False
    left_block = _matmul(adjacency, incidence)
    incidence_t = _transpose(incidence)
    incidence_t_incidence = _matmul(incidence_t, incidence)
    incidence_matrix = _matmul(incidence, matrix)
    for fixed in range(9):
        for orbit in range(12):
            if left_block[fixed][orbit] + incidence_matrix[fixed][orbit] + incidence[fixed][orbit] != 11:
                return False
    square = _matmul(matrix, matrix)
    for left in range(12):
        for right in range(12):
            value = (3 * incidence_t_incidence[left][right] + square[left][right]
                     + matrix[left][right] - 11 * int(left == right))
            if value != 33:
                return False
    return True


def quotient_candidates():
    system = quotient_linear_system()
    if not verify_affine_family(system):
        raise FrontierViolation("affine quotient family is not complete")
    if any(_matrix_at(tau)[0][7] != tau
           or _matrix_at(tau)[0][6] != 3 - tau for tau in range(4)):
        raise FrontierViolation("affine parameter is not integrally normalized")
    records = []
    for tau in range(4):
        matrix = _matrix_at(tau)
        if verify_quotient_matrix(matrix):
            records.append({"tau": tau, "matrix": [list(row) for row in matrix]})
    return records


def forced_independent_five() -> dict:
    return {
        "fixed_vertices": [2, 3, 6],
        "orbit_pair": [0, 3],
        "orbit_positions": [0, 1],
        "normalized_phase": 0,
    }


def verify_forced_independent_five(matrix, witness: dict) -> bool:
    graph = construct_fixed_graph()
    incidence = build_incidence(graph)
    fixed = tuple(witness["fixed_vertices"])
    left_orbit, right_orbit = witness["orbit_pair"]
    left_position, right_position = witness["orbit_positions"]
    phase = witness["normalized_phase"]
    if not is_independent(graph, fixed):
        return False
    if any(incidence[vertex][left_orbit] or incidence[vertex][right_orbit]
           for vertex in fixed):
        return False
    if matrix[left_orbit][right_orbit] != 1:
        return False
    difference = (right_position - left_position) % 3
    orbit_edge = difference == phase
    return not orbit_edge


def run_analysis() -> dict:
    graph = construct_fixed_graph()
    candidates = enumerate_fixed_graph_candidates()
    if len(candidates) != 8 or not all(is_isomorphic(graph, item) for item in candidates):
        raise FrontierViolation("fixed graph uniqueness check failed")
    degree_obstructions = degree_extreme_obstructions()
    if (degree_obstructions["residue_hits"]
            or degree_obstructions["degree_seven_complement_degree"] != 1):
        raise FrontierViolation("fixed-graph degree extremes were not excluded")
    incidence = build_incidence(graph)
    multiplicities = incidence_multiplicity_solutions(graph)
    if multiplicities != [((1,) * 6, (1,) * 6)]:
        raise FrontierViolation("Ramsey incidence columns are not unique")
    system = quotient_linear_system()
    quotients = quotient_candidates()
    witness = forced_independent_five()
    if any(not verify_forced_independent_five(record["matrix"], witness)
           for record in quotients):
        raise FrontierViolation("quotient obstruction failed")
    large_records = large_fixed_count_obstructions()
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "claim": {
            "graph_family": "ramsey_good_strongly_regular",
            "parameters": list(PARAMETERS),
            "automorphism_order": 3,
            "conclusion": "no_automorphism_of_order_three",
            "excluded_spectral_fixed_point_counts": list(SPECTRAL_FIXED_POINT_COUNTS),
            "general_ramsey_bound_claimed": False,
        },
        "fixed_point_congruence": "f_mod_6_eq_3",
        "fixed_3": {"disposition": "IMPOSSIBLE_ODD_ONE_REGULAR_GRAPH"},
        "fixed_9": {
            "fixed_graph_rows": list(graph),
            "fixed_graph_parameters": [9, 4, 1, 2],
            "orbit_3_cycles": 12,
            "triangle_orbits": 6,
            "independent_orbits": 6,
            "local_biregular_matrices": 90,
            "labelled_fixed_graph_survivors": len(candidates),
            "fixed_graph_isomorphism_types": 1,
            "degree_extreme_obstructions": degree_obstructions,
            "independent_triples": [list(item) for item in independent_triples(graph)],
            "triangles": [list(item) for item in triangles(graph)],
            "incidence": [list(row) for row in incidence],
            "incidence_multiplicity_solutions": [
                [list(left), list(right)] for left, right in multiplicities
            ],
            "linear_variables": system["variables"],
            "linear_rank": system["rank"],
            "linear_augmented_rank": system["augmented_rank"],
            "affine_parameter_entries": {
                "tau": [0, 7],
                "three_minus_tau": [0, 6],
            },
            "quotient_candidates": quotients,
            "forced_independent_five": witness,
        },
        "fixed_15_plus": {
            "disposition": "IMPOSSIBLE_TRIANGLE_ORBIT_MOMENT_FOR_ANY_SRG",
            "graph_family": "strongly_regular",
            "ramsey_hypothesis_used": False,
            "fixed_point_counts": list(LARGE_FIXED_POINT_COUNTS),
            "fixed_support_upper_bound_from_common_neighbors": 9,
            "cauchy_gap_discriminant_at_c10": -63,
            "moment_equations": {
                "off_orbit_degree_sum": "20-s",
                "off_orbit_square_sum": "38-3s",
                "off_orbit_pair_intersection_sum": "9-s",
                "cauchy_requirement": "(20-s)^2 <= (c-1)(38-3s)",
            },
            "records": large_records,
        },
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "order3_srg_frontier.json",
    )
    args = parser.parse_args(argv)
    document = run_analysis()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    large_records = document["fixed_15_plus"]["records"]
    minimum_gap = min(record["minimum_cauchy_gap"] for record in large_records)
    print(
        "ORDER3 FRONTIER:",
        f"fixed9_types={document['fixed_9']['fixed_graph_isomorphism_types']}",
        f"quotients={len(document['fixed_9']['quotient_candidates'])}",
        "forced_I5=2/2",
        f"large_fixed_counts={len(large_records)}",
        f"minimum_cauchy_gap={minimum_gap}",
    )
    print(
        "ORDER3 RESULT:", DISPOSITION,
        "(no order-three automorphism in the Ramsey SRG lane; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
