"""Independent checker for the order-three SRG frontier artifact.

This module imports nothing from ``order3_srg_frontier``.  It constructs
L_2(3) as a rook graph rather than from a vertex neighborhood, independently
re-enumerates the f=9 fixed graph and quotient, and checks its forced
independent five-set.  It also re-derives the trace and triangle-orbit moment
obstruction for f=15,21,27,33,39 using exact integer arithmetic.  Together
these disjoint checks certify that a Ramsey-good srg(45,22,10,11) has no
automorphism of order three.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

SCHEMA_VERSION = 2
CAMPAIGN_ID = "order3_srg_frontier"
DISPOSITION = "ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG"
SPECTRAL_FIXED_POINT_COUNTS = (3, 9, 15, 21, 27, 33, 39)
PARAMETERS = (45, 22, 10, 11)
CLAIM = {
    "graph_family": "ramsey_good_strongly_regular",
    "parameters": [45, 22, 10, 11],
    "automorphism_order": 3,
    "conclusion": "no_automorphism_of_order_three",
    "excluded_spectral_fixed_point_counts": list(SPECTRAL_FIXED_POINT_COUNTS),
    "general_ramsey_bound_claimed": False,
}


class CheckViolation(RuntimeError):
    """The artifact differs from the independent exact reconstruction."""


def _edge(rows, left, right):
    rows[left] |= 1 << right
    rows[right] |= 1 << left


def _rook_graph() -> tuple[int, ...]:
    rows = [0] * 9
    for left in range(9):
        lr, lc = divmod(left, 3)
        for right in range(left + 1, 9):
            rr, rc = divmod(right, 3)
            if lr == rr or lc == rc:
                _edge(rows, left, right)
    return tuple(rows)


def _parameters(graph):
    degrees = {row.bit_count() for row in graph}
    if len(degrees) != 1:
        return None
    adjacent = set()
    nonadjacent = set()
    for left in range(len(graph)):
        for right in range(left + 1, len(graph)):
            target = adjacent if graph[left] >> right & 1 else nonadjacent
            target.add((graph[left] & graph[right]).bit_count())
    if len(adjacent) != 1 or len(nonadjacent) != 1:
        return None
    return len(graph), next(iter(degrees)), next(iter(adjacent)), next(iter(nonadjacent))


def _isomorphic(left, right):
    if _parameters(left) != _parameters(right):
        return False
    for permutation in itertools.permutations(range(9)):
        valid = True
        for new_left, old_left in enumerate(permutation):
            row = 0
            for new_right, old_right in enumerate(permutation):
                if left[old_left] >> old_right & 1:
                    row |= 1 << new_right
            if row != right[new_left]:
                valid = False
                break
        if valid:
            return True
    return False


def _local_candidate(row_sets):
    rows = [0] * 9
    for vertex in range(1, 5):
        _edge(rows, 0, vertex)
    _edge(rows, 1, 2)
    _edge(rows, 3, 4)
    for left, right in ((5, 6), (6, 7), (7, 8), (8, 5)):
        _edge(rows, left, right)
    for neighbor, selected in enumerate(row_sets, 1):
        for column in selected:
            _edge(rows, neighbor, 5 + column)
    return tuple(rows)


def _local_survivors():
    pairs = tuple(itertools.combinations(range(4), 2))
    candidates = []
    local_matrices = 0
    for rows in itertools.product(pairs, repeat=4):
        if tuple(sum(column in row for row in rows) for column in range(4)) != (2,) * 4:
            continue
        local_matrices += 1
        graph = _local_candidate(rows)
        if _parameters(graph) == (9, 4, 1, 2):
            candidates.append(graph)
    return local_matrices, candidates


def _clique(graph, vertices):
    return all(graph[left] >> right & 1
               for left, right in itertools.combinations(vertices, 2))


def _independent(graph, vertices):
    return all(not (graph[left] >> right & 1)
               for left, right in itertools.combinations(vertices, 2))


def _triple_families(graph):
    triples = tuple(itertools.combinations(range(9), 3))
    independent = tuple(item for item in triples if _independent(graph, item))
    cliques = tuple(item for item in triples if _clique(graph, item))
    return independent, cliques


def _incidence(graph):
    independent, cliques = _triple_families(graph)
    return tuple(tuple(
        [int(vertex in triple) for triple in independent]
        + [int(vertex not in triangle) for triangle in cliques]
    ) for vertex in range(9))


def _gram(matrix):
    return tuple(tuple(sum(matrix[left][column] * matrix[right][column]
                           for column in range(12)) for right in range(9))
                 for left in range(9))


def _compositions(total, parts, prefix=()):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from _compositions(total - value, parts - 1, prefix + (value,))


def _signature(columns, counts):
    values = []
    for vertex in range(9):
        values.append(sum(counts[index] * int(vertex in column)
                          for index, column in enumerate(columns)))
    for left in range(9):
        for right in range(left, 9):
            values.append(sum(counts[index] * int(left in column and right in column)
                              for index, column in enumerate(columns)))
    return tuple(values)


def _multiplicity_solutions(graph):
    independent, cliques = _triple_families(graph)
    complements = tuple(tuple(v for v in range(9) if v not in triangle)
                        for triangle in cliques)
    target = (6,) * 9 + tuple(6 if left == right else 3
                              for left in range(9) for right in range(left, 9))
    rights = {}
    for counts in _compositions(6, 6):
        rights.setdefault(_signature(complements, counts), []).append(counts)
    result = []
    for counts in _compositions(6, 6):
        signature = _signature(independent, counts)
        needed = tuple(target[index] - signature[index] for index in range(len(target)))
        result.extend((counts, right) for right in rights.get(needed, ()))
    return result


def _matmul(left, right):
    return tuple(tuple(sum(left[row][inner] * right[inner][column]
                           for inner in range(len(right)))
                       for column in range(len(right[0])))
                 for row in range(len(left)))


def _transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def _rank(matrix):
    rows = [list(map(Fraction, row)) for row in matrix]
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next((row for row in range(pivot_row, len(rows))
                      if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [entry / divisor for entry in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            multiple = rows[row][column]
            rows[row] = [left - multiple * right
                         for left, right in zip(rows[row], rows[pivot_row])]
        pivot_row += 1
    return pivot_row


def _linear_system(graph, incidence):
    variables = tuple((left, right) for left in range(12) for right in range(left, 12))
    owner = {pair: index for index, pair in enumerate(variables)}
    adjacency = tuple(tuple(int(graph[row] >> column & 1) for column in range(9))
                      for row in range(9))
    ab = _matmul(adjacency, incidence)
    rows = []
    rhs = []
    for fixed in range(9):
        for orbit in range(12):
            row = [0] * 78
            for other in range(12):
                row[owner[min(other, orbit), max(other, orbit)]] += incidence[fixed][other]
            rows.append(tuple(row))
            rhs.append(11 - ab[fixed][orbit] - incidence[fixed][orbit])
    for orbit in range(12):
        row = [0] * 78
        row[owner[orbit, orbit]] = 1
        rows.append(tuple(row))
        rhs.append(2 if orbit < 6 else 0)
    for orbit in range(12):
        row = [0] * 78
        for other in range(12):
            row[owner[min(orbit, other), max(orbit, other)]] += 1
        rows.append(tuple(row))
        rhs.append(19 if orbit < 6 else 16)
    return variables, tuple(rows), tuple(rhs)


def _flatten(matrix, variables):
    return tuple(matrix[left][right] for left, right in variables)


def _linear_solution(rows, rhs, values):
    return all(sum(coefficient * value for coefficient, value in zip(row, values)) == target
               for row, target in zip(rows, rhs))


def _quotient_valid(graph, incidence, matrix):
    if len(matrix) != 12 or any(len(row) != 12 for row in matrix):
        return False
    if any(matrix[left][right] != matrix[right][left]
           for left in range(12) for right in range(12)):
        return False
    if tuple(matrix[index][index] for index in range(12)) != (2,) * 6 + (0,) * 6:
        return False
    if any(not 0 <= matrix[left][right] <= 3 for left in range(12)
           for right in range(12) if left != right):
        return False
    adjacency = tuple(tuple(int(graph[row] >> column & 1) for column in range(9))
                      for row in range(9))
    if tuple(sum(row) for row in matrix) != (19,) * 6 + (16,) * 6:
        return False
    if any(_matmul(adjacency, incidence)[fixed][orbit]
           + _matmul(incidence, matrix)[fixed][orbit]
           + incidence[fixed][orbit] != 11
           for fixed in range(9) for orbit in range(12)):
        return False
    bt_b = _matmul(_transpose(incidence), incidence)
    square = _matmul(matrix, matrix)
    return all(3 * bt_b[left][right] + square[left][right] + matrix[left][right]
               - 11 * int(left == right) == 33
               for left in range(12) for right in range(12))


def _witness_valid(graph, incidence, matrix, witness):
    fixed = tuple(witness["fixed_vertices"])
    first, second = witness["orbit_pair"]
    left_position, right_position = witness["orbit_positions"]
    phase = witness["normalized_phase"]
    if not _independent(graph, fixed):
        return False
    if any(incidence[vertex][first] or incidence[vertex][second] for vertex in fixed):
        return False
    return (matrix[first][second] == 1
            and (right_position - left_position) % 3 != phase)


def _degree_extreme_record():
    fixed_order = 9
    degree_one = 1
    degree_seven = 7
    required_residue = PARAMETERS[3] % 3
    possible = list(range(degree_one + 1))
    return {
        "degree_one_max_common_fixed": degree_one,
        "nonedge_required_residue": required_residue,
        "possible_common_fixed_counts": possible,
        "residue_hits": [value for value in possible
                         if value % 3 == required_residue],
        "degree_seven_complement_degree": fixed_order - 1 - degree_seven,
    }


def _trace_triangle_orbits(fixed_points):
    vertices, degree, adjacent_common, nonadjacent_common = PARAMETERS
    if (vertices - fixed_points) % 3:
        raise CheckViolation("fixed-point count is incompatible with 3-orbits")
    orbit_cycles = (vertices - fixed_points) // 3
    polynomial_linear = nonadjacent_common - adjacent_common
    polynomial_constant = nonadjacent_common - degree
    discriminant = polynomial_linear**2 - 4 * polynomial_constant
    # A positive rational square in Q(sqrt(-3)) must already be a rational
    # square: the cross term in (a+b*sqrt(-3))^2 forces a*b=0.
    square_root = math.isqrt(discriminant)
    if discriminant != 45 or square_root * square_root == discriminant:
        raise CheckViolation("nonprincipal polynomial is not the required irreducible one")
    if orbit_cycles % 2:
        raise CheckViolation("irreducibility requires even nontrivial eigenspaces")
    full_root_multiplicity = (vertices - 1) // 2
    fixed_root_multiplicity = full_root_multiplicity - orbit_cycles
    root_sum = -polynomial_linear
    quotient_trace = degree + fixed_root_multiplicity * root_sum
    if quotient_trace != orbit_cycles or quotient_trace % 2:
        raise CheckViolation("orbit-constant spectral trace differs")
    return orbit_cycles, quotient_trace // 2


def _triangle_moment(orbit_cycles, fixed_support):
    _, degree, adjacent_common, _ = PARAMETERS
    internal_degree = 2
    internal_common = 1
    other_orbits = orbit_cycles - 1
    degree_sum = degree - internal_degree - fixed_support
    # For two adjacent vertices in the triangle orbit, lambda consists of the
    # fixed support, their third orbitmate, and C(t_i,2) in each cyclic
    # bipartite 3x3 block.  Hence sum C(t_i,2)=lambda-1-s.
    pair_intersection_sum = adjacent_common - internal_common - fixed_support
    square_sum = degree_sum + 2 * pair_intersection_sum
    return {
        "fixed_support": fixed_support,
        "other_orbits": other_orbits,
        "off_orbit_degree_sum": degree_sum,
        "off_orbit_square_sum": square_sum,
        "cauchy_gap": degree_sum * degree_sum - other_orbits * square_sum,
    }


def _large_fixed_records():
    _, degree, adjacent_common, _ = PARAMETERS
    support_values = list(range(adjacent_common))
    records = []
    for fixed_points in SPECTRAL_FIXED_POINT_COUNTS[2:]:
        orbit_cycles, triangle_orbits = _trace_triangle_orbits(fixed_points)
        moments = [_triangle_moment(orbit_cycles, support)
                   for support in support_values]
        gaps = [moment["cauchy_gap"] for moment in moments]
        if any(moment["off_orbit_square_sum"] < 0 for moment in moments):
            raise CheckViolation("derived two-walk square sum is negative")
        if any(gap <= 0 for gap in gaps):
            raise CheckViolation("triangle-orbit Cauchy gap is not positive")
        records.append({
            "fixed_points": fixed_points,
            "orbit_3_cycles": orbit_cycles,
            "triangle_orbits": triangle_orbits,
            "fixed_support_values": support_values,
            "cauchy_gaps": gaps,
            "minimum_cauchy_gap": min(gaps),
        })
    leading_cycles = records[0]["orbit_3_cycles"]
    linear = -2 * (degree - 2) + 3 * (leading_cycles - 1)
    constant = ((degree - 2)**2
                - (leading_cycles - 1) * (degree + 2 * adjacent_common - 4))
    if linear * linear - 4 * constant != -63:
        raise CheckViolation("independent Cauchy discriminant differs")
    return records


def _require_keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise CheckViolation(f"{label}: key set differs")


def verify_document(document: dict) -> dict:
    _require_keys(document, {
        "schema_version", "campaign_id", "disposition", "claim",
        "fixed_point_congruence", "fixed_3", "fixed_9", "fixed_15_plus",
    }, "document")
    if document["schema_version"] != SCHEMA_VERSION:
        raise CheckViolation("schema version differs")
    if document["campaign_id"] != CAMPAIGN_ID or document["disposition"] != DISPOSITION:
        raise CheckViolation("campaign identity or disposition differs")
    if document["claim"] != CLAIM:
        raise CheckViolation("structured claim differs")
    if document["fixed_point_congruence"] != "f_mod_6_eq_3":
        raise CheckViolation("fixed-point congruence differs")
    if document["fixed_3"] != {"disposition": "IMPOSSIBLE_ODD_ONE_REGULAR_GRAPH"}:
        raise CheckViolation("fixed-three proof record differs")

    fixed = document["fixed_9"]
    _require_keys(fixed, {
        "fixed_graph_rows", "fixed_graph_parameters", "orbit_3_cycles",
        "triangle_orbits", "independent_orbits", "local_biregular_matrices",
        "labelled_fixed_graph_survivors", "fixed_graph_isomorphism_types",
        "degree_extreme_obstructions", "independent_triples", "triangles",
        "incidence", "incidence_multiplicity_solutions", "linear_variables",
        "linear_rank", "linear_augmented_rank", "affine_parameter_entries",
        "quotient_candidates", "forced_independent_five",
    }, "fixed_9")
    graph = tuple(fixed["fixed_graph_rows"])
    if _parameters(graph) != (9, 4, 1, 2) or not _isomorphic(graph, _rook_graph()):
        raise CheckViolation("fixed graph is not L_2(3)")
    fixed_nine_cycles, fixed_nine_triangles = _trace_triangle_orbits(9)
    if (fixed["orbit_3_cycles"], fixed["triangle_orbits"],
            fixed["independent_orbits"]) != (
                fixed_nine_cycles,
                fixed_nine_triangles,
                fixed_nine_cycles - fixed_nine_triangles,
            ):
        raise CheckViolation("order-three orbit counts differ")
    degree_extremes = _degree_extreme_record()
    if fixed["degree_extreme_obstructions"] != degree_extremes:
        raise CheckViolation("fixed-degree residue obstruction differs")
    if degree_extremes["residue_hits"] or degree_extremes[
            "degree_seven_complement_degree"] != 1:
        raise CheckViolation("fixed-degree extremes were not excluded")
    local_count, survivors = _local_survivors()
    if local_count != 90 or len(survivors) != 8:
        raise CheckViolation("independent local census differs")
    if not all(_isomorphic(_rook_graph(), candidate) for candidate in survivors):
        raise CheckViolation("local census has another isomorphism type")
    if fixed["local_biregular_matrices"] != 90 or fixed["labelled_fixed_graph_survivors"] != 8:
        raise CheckViolation("fixed-graph census counts differ")
    if fixed["fixed_graph_isomorphism_types"] != 1:
        raise CheckViolation("fixed-graph type count differs")
    independent, cliques = _triple_families(graph)
    if fixed["independent_triples"] != [list(item) for item in independent]:
        raise CheckViolation("independent triples differ")
    if fixed["triangles"] != [list(item) for item in cliques]:
        raise CheckViolation("triangles differ")
    incidence = _incidence(graph)
    if fixed["incidence"] != [list(row) for row in incidence]:
        raise CheckViolation("incidence matrix differs")
    target_gram = tuple(tuple(6 if left == right else 3 for right in range(9))
                        for left in range(9))
    if _gram(incidence) != target_gram:
        raise CheckViolation("incidence Gram matrix differs")
    multiplicities = _multiplicity_solutions(graph)
    expected_multiplicities = [[[1] * 6, [1] * 6]]
    if multiplicities != [((1,) * 6, (1,) * 6)] or fixed["incidence_multiplicity_solutions"] != expected_multiplicities:
        raise CheckViolation("incidence multiplicities differ")

    variables, rows, rhs = _linear_system(graph, incidence)
    rank = _rank(rows)
    augmented = _rank(tuple(tuple(row) + (rhs[index],)
                            for index, row in enumerate(rows)))
    if (fixed["linear_variables"], fixed["linear_rank"],
            fixed["linear_augmented_rank"]) != (78, rank, augmented):
        raise CheckViolation("linear ranks differ")
    if rank != 77 or augmented != 77:
        raise CheckViolation("linear family does not have dimension one")
    if fixed["affine_parameter_entries"] != {
        "tau": [0, 7], "three_minus_tau": [0, 6],
    }:
        raise CheckViolation("affine parameter coordinates differ")
    records = fixed["quotient_candidates"]
    if not isinstance(records, list) or [record.get("tau") for record in records] != [1, 2]:
        raise CheckViolation("quotient parameter set differs")
    matrices = [tuple(tuple(row) for row in record["matrix"]) for record in records]
    if any(not _quotient_valid(graph, incidence, matrix) for matrix in matrices):
        raise CheckViolation("recorded quotient matrix is invalid")
    direction = tuple(tuple(matrices[1][row][column] - matrices[0][row][column]
                            for column in range(12)) for row in range(12))
    if not any(value for row in direction for value in row):
        raise CheckViolation("quotient direction is zero")
    flat_direction = _flatten(direction, variables)
    if not _linear_solution(rows, (0,) * len(rhs), flat_direction):
        raise CheckViolation("quotient direction is not homogeneous")
    if not _linear_solution(rows, rhs, _flatten(matrices[0], variables)):
        raise CheckViolation("quotient base is not a linear solution")
    # Rank 77 proves this affine line is complete.  A direction entry of unit
    # magnitude makes its parameter integral for every integer T; entry bounds
    # leave exactly tau=0,1,2,3, and the quadratic block accepts 1,2.
    if not any(abs(value) == 1 for row in direction for value in row):
        raise CheckViolation("affine parameter is not integrally normalized")
    if not (matrices[0][0][7] == 1 and matrices[1][0][7] == 2
            and matrices[0][0][6] == 2 and matrices[1][0][6] == 1):
        raise CheckViolation("affine parameter bounds are not tau and 3-tau")
    tested = []
    for tau in range(4):
        matrix = tuple(tuple(matrices[0][row][column] + (tau - 1) * direction[row][column]
                             for column in range(12)) for row in range(12))
        if _quotient_valid(graph, incidence, matrix):
            tested.append(tau)
    if tested != [1, 2]:
        raise CheckViolation("independent quotient enumeration differs")
    witness = fixed["forced_independent_five"]
    if witness != {
        "fixed_vertices": [2, 3, 6],
        "orbit_pair": [0, 3],
        "orbit_positions": [0, 1],
        "normalized_phase": 0,
    }:
        raise CheckViolation("forced independent-set record differs")
    if any(not _witness_valid(graph, incidence, matrix, witness) for matrix in matrices):
        raise CheckViolation("forced independent set is invalid")

    large = document["fixed_15_plus"]
    _require_keys(large, {
        "disposition", "graph_family", "ramsey_hypothesis_used",
        "fixed_point_counts", "fixed_support_upper_bound_from_common_neighbors",
        "cauchy_gap_discriminant_at_c10", "moment_equations", "records",
    }, "fixed_15_plus")
    expected_records = _large_fixed_records()
    expected_large = {
        "disposition": "IMPOSSIBLE_TRIANGLE_ORBIT_MOMENT_FOR_ANY_SRG",
        "graph_family": "strongly_regular",
        "ramsey_hypothesis_used": False,
        "fixed_point_counts": [record["fixed_points"] for record in expected_records],
        "fixed_support_upper_bound_from_common_neighbors": PARAMETERS[2] - 1,
        "cauchy_gap_discriminant_at_c10": -63,
        "moment_equations": {
            "off_orbit_degree_sum": "20-s",
            "off_orbit_square_sum": "38-3s",
            "off_orbit_pair_intersection_sum": "9-s",
            "cauchy_requirement": "(20-s)^2 <= (c-1)(38-3s)",
        },
        "records": expected_records,
    }
    if large != expected_large:
        raise CheckViolation("large-fixed-count moment record differs")
    if _triangle_moment(12, 4)["cauchy_gap"] != -30:
        raise CheckViolation("moment boundary was extended into the fixed-nine case")
    spectral_candidates = tuple(
        fixed_points for fixed_points in range(0, PARAMETERS[0], 3)
        if ((PARAMETERS[0] - fixed_points) // 3) % 2 == 0
    )
    if spectral_candidates != SPECTRAL_FIXED_POINT_COUNTS:
        raise CheckViolation("spectral fixed-point coverage is incomplete")
    return {
        "quotient_candidates": 2,
        "forced_independent_sets": 2,
        "large_fixed_counts": len(expected_records),
        "minimum_cauchy_gap": min(record["minimum_cauchy_gap"]
                                  for record in expected_records),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    args = parser.parse_args(argv)
    try:
        document = json.loads(args.artifact.read_text())
        summary = verify_document(document)
    except (OSError, json.JSONDecodeError, CheckViolation) as error:
        print(f"ORDER3 CHECK FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "ORDER3 EVIDENCE VERIFIED:", DISPOSITION,
        f"quotients={summary['quotient_candidates']}",
        f"forced_I5={summary['forced_independent_sets']}",
        f"large_fixed_counts={summary['large_fixed_counts']}",
        f"minimum_cauchy_gap={summary['minimum_cauchy_gap']}",
        "(no order-three automorphism in the Ramsey SRG lane; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
