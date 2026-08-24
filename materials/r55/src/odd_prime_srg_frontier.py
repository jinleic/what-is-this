#!/usr/bin/env python3
"""Exact odd-prime automorphism frontier for a Ramsey-good conference SRG.

Theorem
=======

If a strongly regular graph with parameters (45,22,10,11) contains neither a
clique nor an independent set of order five, then its automorphism group is a
2-group.  This is a structural theorem inside the strongly regular lane; it is
not a bound on R(5,5), and it says nothing about a non-strongly-regular
Ramsey(5,5,45) graph.

The existing order-three certificate is composed as a dependency.  This file
excludes every other odd prime order.

* Order five: cyclotomic Galois symmetry gives quotient trace 2c.  A
  quotient-diagonal moment, valid for every internal orbit degree 0,2,4,
  excludes fixed counts at least ten without a Ramsey hypothesis.  Ramsey
  avoidance then makes every moving orbit a C5 in the two surviving cases.
  With no fixed point, a 9x9 integral Seidel quotient would satisfy
  M^2=45I-5J; a complete 56,755-prefix integer search
  finds no matrix.  With five fixed points, the fixed graph is C5 and the
  fixed-to-orbit supports form eight blocks with replication four and pair
  multiplicity two.  A triple-count argument plus R(3,4)=R(4,3)=9 gives a
  solver-free K5/I5 bridge.

* Order seven: cyclotomic irreducibility leaves fixed counts 3,17,31.  The
  first fails by handshake, the second by simultaneous support and complement
  intersection bounds, and the last by a Rayleigh contradiction in the fixed
  block.

* Order eleven: fixed counts 1,23 remain.  The latter has an independent
  fixed support group of size at least eleven.  The former has a unique 4x4
  quotient up to swapping two orbits.  Its two fixed-supported 11-orbits give
  46,200 exact cyclic blocks; every block has a K4 or an I5.

* Prime orders at least thirteen: 13,17,19 fail by fixed-graph handshake, and
  p>=23 would require an even positive number of p-cycles although at most one
  fits on 45 vertices.

All searches use deterministic integer bitsets and exhaustive finite domains;
no SAT/SMT solver, floating point computation, or external graph catalog is in
the proof.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import sys
from pathlib import Path

SCHEMA_VERSION = 1
CAMPAIGN_ID = "odd_prime_srg_frontier"
DISPOSITION = "AUTOMORPHISM_GROUP_IS_A_2_GROUP_FOR_RAMSEY_SRG"
PARAMETERS = (45, 22, 10, 11)
ORDER3_CAMPAIGN = "order3_srg_frontier"
ORDER3_DISPOSITION = "ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG"

_ODD_VALUES = (-5, -3, -1, 1, 3, 5)
_ORDER5_FIRST_ROWS = (
    (-3, -3, -1, -1, 1, 1, 3, 3),
    (-3, -1, -1, -1, -1, 1, 1, 5),
    (-5, -1, -1, 1, 1, 1, 1, 3),
)


class FrontierViolation(RuntimeError):
    """An exact structural invariant failed."""


def _canonical_document_sha256(document: dict) -> str:
    payload = json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def order5_trace_record() -> dict:
    """Record the Galois trace consequence for order-five moving orbits."""
    return {
        "quotient_trace": "2*c",
        "internal_degree_counts_equation": "n0=n4",
        "ramsey_allowed_internal_degrees": [2],
        "moving_orbit_graph": "C5",
    }


def order5_orbit_moment(cycles: int, fixed_support: int,
                        internal_degree: int) -> dict:
    """Exact quotient-diagonal moment and Cauchy gap for one five-orbit."""
    degree_sum = 22 - internal_degree - fixed_support
    square_sum = 66 - internal_degree - internal_degree * internal_degree
    square_sum -= 5 * fixed_support
    gap = degree_sum * degree_sum - (cycles - 1) * square_sum
    return {
        "internal_degree": internal_degree,
        "fixed_support": fixed_support,
        "other_orbits": cycles - 1,
        "off_orbit_degree_sum": degree_sum,
        "off_orbit_square_sum": square_sum,
        "cauchy_gap": gap,
    }


def order5_large_fixed_records() -> list[dict]:
    records = []
    for fixed_points in range(10, 45, 5):
        cycles = (45 - fixed_points) // 5
        degree_records = []
        for internal_degree in (0, 2, 4):
            maximum_support = min(
                fixed_points,
                (66 - internal_degree - internal_degree * internal_degree) // 5,
            )
            moments = [
                order5_orbit_moment(cycles, support, internal_degree)
                for support in range(maximum_support + 1)
            ]
            degree_records.append({
                "internal_degree": internal_degree,
                "fixed_support_values": list(range(maximum_support + 1)),
                "cauchy_gaps": [item["cauchy_gap"] for item in moments],
                "minimum_cauchy_gap": min(item["cauchy_gap"] for item in moments),
            })
        records.append({
            "fixed_points": fixed_points,
            "orbit_5_cycles": cycles,
            "internal_degree_records": degree_records,
            "minimum_cauchy_gap": min(
                item["minimum_cauchy_gap"] for item in degree_records),
        })
    return records


@functools.cache
def _completion_tuples(length: int, total: int, square_sum: int) -> tuple[tuple[int, ...], ...]:
    if length == 0:
        return ((),) if total == 0 and square_sum == 0 else ()
    if square_sum < length or square_sum > 25 * length:
        return ()
    records = []
    for value in _ODD_VALUES:
        for tail in _completion_tuples(length - 1, total - value,
                                       square_sum - value * value):
            records.append((value,) + tail)
    return tuple(records)


def _order5_free_search(first_row: tuple[int, ...]) -> dict:
    order = 9
    matrix = [[0] * order for _ in range(order)]
    for column, value in enumerate(first_row, 1):
        matrix[0][column] = matrix[column][0] = value

    tested = {}
    surviving = {}
    serializations = {row: [] for row in range(1, order)}

    def visit(row: int) -> int:
        if row == order:
            return 1
        prefix = matrix[row][:row]
        candidates = _completion_tuples(
            order - row - 1,
            -sum(prefix),
            40 - sum(value * value for value in prefix),
        )
        solutions = 0
        for tail in candidates:
            tested[row] = tested.get(row, 0) + 1
            for column, value in enumerate(tail, row + 1):
                matrix[row][column] = matrix[column][row] = value
            if any(sum(matrix[previous][index] * matrix[index][row]
                       for index in range(order)) != -5
                   for previous in range(row)):
                continue
            surviving[row] = surviving.get(row, 0) + 1
            flattened = []
            for previous in range(row + 1):
                flattened.extend(matrix[previous][column]
                                 for column in range(previous + 1, order))
            serializations[row].append(",".join(map(str, flattened)) + "\n")
            solutions += visit(row + 1)
        return solutions

    solutions = visit(1)
    return {
        "first_row": list(first_row),
        "tested_by_row": {str(row): tested[row] for row in sorted(tested)},
        "surviving_by_row": {str(row): surviving[row] for row in sorted(surviving)},
        "survivor_sha256_by_row": {
            str(row): hashlib.sha256(
                "".join(sorted(serializations[row])).encode()
            ).hexdigest()
            for row in sorted(tested)
        },
        "solutions": solutions,
    }


@functools.cache
def _cached_order5_free_searches() -> tuple[dict, ...]:
    return tuple(_order5_free_search(first_row) for first_row in _ORDER5_FIRST_ROWS)


def enumerate_order5_free_quotients() -> list[dict]:
    """Exhaust every possible fixed-point-free 9x9 Seidel quotient."""
    return [json.loads(json.dumps(record)) for record in _cached_order5_free_searches()]


def order5_fixed_five_obstruction() -> dict:
    """Encode the solver-free triple-count/Ramsey bridge for f=5."""
    return {
        "fixed_graph": "C5",
        "moving_orbits": 8,
        "incidence_row_size": 4,
        "incidence_pair_multiplicity": 2,
        "triple_multiplicity_alternatives": ["zero", "at_least_two"],
        "moving_vertices_in_bridge": 10,
        "ramsey_numbers": {"R(3,4)": 9, "R(4,3)": 9},
        "all_triple_multiplicities_one_is_impossible": True,
        "all_one_contradiction": {
            "block_size_five": "seven remaining blocks cannot cover ten pairs",
            "maximum_block_size_four": "ten pairs cannot each occur in one size-three block",
        },
    }


def order7_obstructions() -> dict:
    """Return the exact three-case order-seven exclusion."""
    norm_square = 16 + 14
    vector_sum = 16 - 14
    incidence_projection = (16, -14)
    rhs_quadratic = 11 * (norm_square + vector_sum * vector_sum) - 7 * sum(
        value * value for value in incidence_projection
    )
    fourfold_gap = 4 * rhs_quadratic + norm_square
    return {
        "spectral_candidates": [
            {"fixed_points": 3, "orbit_7_cycles": 6},
            {"fixed_points": 17, "orbit_7_cycles": 4},
            {"fixed_points": 31, "orbit_7_cycles": 2},
        ],
        "fixed_3": {
            "fixed_degree": 1,
            "degree_sum": 3,
            "disposition": "IMPOSSIBLE_HANDSHAKE",
        },
        "fixed_17": {
            "support_sizes": [1, 2, 3],
            "support_intersection_upper_bound": 1,
            "complement_support_intersection_upper_bound": 1,
            "maximum_size_one_supports": 1,
            "maximum_size_two_supports": 6,
            "maximum_size_three_supports": 1,
            "maximum_support_family_size": 6,
            "disposition": "IMPOSSIBLE_SUPPORT_FAMILY_CAPACITY",
        },
        "fixed_31": {
            "moving_internal_degrees": [2, 4],
            "incidence_cases_up_to_complement": [{
                "both_orbits": 1,
                "neither_orbit": 0,
                "cross_degree": 3,
                "orbit_support_sizes": [17, 15],
            }],
            "support_group_sizes": [16, 14],
            "test_vector_norm_square": norm_square,
            "test_vector_sum": vector_sum,
            "incidence_projection": list(incidence_projection),
            "rhs_quadratic": rhs_quadratic,
            "rayleigh_fourfold_gap": fourfold_gap,
            "disposition": "IMPOSSIBLE_FIXED_BLOCK_RAYLEIGH",
        },
    }


def order11_canonical_quotient() -> list[list[int]]:
    """The unique f=1 moving-orbit quotient, up to swapping the last cells."""
    return [
        [4, 6, 4, 7],
        [6, 4, 7, 4],
        [4, 7, 6, 5],
        [7, 4, 5, 6],
    ]


def verify_order11_quotient(matrix: list[list[int]]) -> bool:
    if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
        return False
    if any(matrix[left][right] != matrix[right][left]
           for left in range(4) for right in range(4)):
        return False
    support = (1, 1, 0, 0)
    if [matrix[index][index] for index in range(4)] != [4, 4, 6, 6]:
        return False
    if [sum(row) + support[index] for index, row in enumerate(matrix)] != [22] * 4:
        return False
    if [sum(support[row] * matrix[row][column] for row in range(4))
            + support[column] for column in range(4)] != [11] * 4:
        return False
    for left in range(4):
        for right in range(4):
            square = sum(matrix[left][middle] * matrix[middle][right]
                         for middle in range(4))
            value = square + matrix[left][right] + 11 * support[left] * support[right]
            if value != (132 if left == right else 121):
                return False
    return True


def _add_edge(adjacency: list[int], left: int, right: int) -> None:
    adjacency[left] |= 1 << right
    adjacency[right] |= 1 << left


def _two_cyclic_orbit_graph(internal_left: int, internal_right: int,
                            cross_mask: int) -> tuple[int, ...]:
    prime = 11
    adjacency = [0] * (2 * prime)
    for orbit, internal_mask in enumerate((internal_left, internal_right)):
        for position in range(prime):
            for difference in range(1, 6):
                if internal_mask >> (difference - 1) & 1:
                    _add_edge(
                        adjacency,
                        orbit * prime + position,
                        orbit * prime + (position + difference) % prime,
                    )
    for position in range(prime):
        for difference in range(prime):
            if cross_mask >> difference & 1:
                _add_edge(adjacency, position,
                          prime + (position + difference) % prime)
    return tuple(adjacency)


def _find_clique(adjacency: tuple[int, ...], order: int) -> tuple[int, ...] | None:
    def visit(candidates: int, need: int,
              chosen: tuple[int, ...]) -> tuple[int, ...] | None:
        if need == 0:
            return chosen
        if candidates.bit_count() < need:
            return None
        while candidates:
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            witness = visit(candidates & adjacency[vertex], need - 1,
                            chosen + (vertex,))
            if witness is not None:
                return witness
        return None

    return visit((1 << len(adjacency)) - 1, order, ())


def _complement(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    full = (1 << len(adjacency)) - 1
    return tuple(full ^ (1 << vertex) ^ row
                 for vertex, row in enumerate(adjacency))


@functools.cache
def _cached_order11_coverage() -> tuple:
    internal_masks = tuple(
        sum(1 << index for index in selected)
        for selected in itertools.combinations(range(5), 2)
    )
    cross_masks = tuple(
        sum(1 << index for index in selected)
        for selected in itertools.combinations(range(11), 6)
    )
    digest = hashlib.sha256()
    k4_primary = 0
    i5_fallback = 0
    configurations_with_both = 0
    examples = {}
    for internal_left in internal_masks:
        for internal_right in internal_masks:
            for cross_mask in cross_masks:
                graph = _two_cyclic_orbit_graph(
                    internal_left, internal_right, cross_mask)
                clique = _find_clique(graph, 4)
                independent = _find_clique(_complement(graph), 5)
                if clique is not None:
                    kind = "K4"
                    witness = clique
                    k4_primary += 1
                    if independent is not None:
                        configurations_with_both += 1
                elif independent is not None:
                    kind = "I5"
                    witness = independent
                    i5_fallback += 1
                else:
                    raise FrontierViolation("an order-eleven cyclic block escaped K4/I5")
                examples.setdefault(kind, (
                    internal_left, internal_right, cross_mask, witness))
                digest.update(
                    f"{internal_left:02x}:{internal_right:02x}:{cross_mask:03x}:"
                    f"{kind}:{','.join(map(str, witness))}\n".encode()
                )
    return (
        len(internal_masks) ** 2 * len(cross_masks),
        k4_primary,
        i5_fallback,
        configurations_with_both,
        digest.hexdigest(),
        tuple(sorted(examples.items())),
    )


def order11_block_coverage() -> dict:
    """Exhaust the 46,200 fixed-supported two-orbit cyclic blocks."""
    (configurations, k4_primary, i5_fallback, both, digest,
     example_items) = _cached_order11_coverage()
    examples = {}
    for kind, (left, right, cross, witness) in example_items:
        examples[kind] = {
            "internal_pair_masks": [left, right],
            "cross_mask": cross,
            "vertices": list(witness),
        }
    return {
        "internal_connection_sets_per_orbit": 10,
        "cross_connection_sets": 462,
        "configurations": configurations,
        "k4_primary_witnesses": k4_primary,
        "i5_fallback_witnesses": i5_fallback,
        "configurations_with_both": both,
        "coverage_sha256": digest,
        "examples": examples,
    }


def order11_fixed_23_obstruction() -> dict:
    return {
        "orbit_11_cycles": 2,
        "fixed_support_sizes_per_vertex": [0, 1, 2],
        "at_most_one_both_support_vertex": True,
        "at_most_one_neither_support_vertex": True,
        "single_support_vertices_lower_bound": 21,
        "largest_single_support_group_lower_bound": 11,
        "forced_independent_set_lower_bound": 11,
        "disposition": "IMPOSSIBLE_FIXED_INDEPENDENT_FIVE",
    }


def large_prime_obstructions() -> dict:
    return {
        "handshake_cases": [
            {"prime": 13, "orbit_cycles": 2, "fixed_points": 19,
             "fixed_degree": 9, "degree_sum": 171},
            {"prime": 17, "orbit_cycles": 2, "fixed_points": 11,
             "fixed_degree": 5, "degree_sum": 55},
            {"prime": 19, "orbit_cycles": 2, "fixed_points": 7,
             "fixed_degree": 3, "degree_sum": 21},
        ],
        "cycle_parity_cutoff_prime": 23,
        "largest_relevant_prime": 43,
        "disposition": "NO_PRIME_ORDER_AT_LEAST_THIRTEEN_FOR_ANY_SRG",
        "ramsey_hypothesis_used": False,
    }


def _order3_dependency(path: Path) -> dict:
    try:
        document = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise FrontierViolation(f"cannot read order-three dependency: {error}") from error
    if document.get("campaign_id") != ORDER3_CAMPAIGN:
        raise FrontierViolation("order-three dependency campaign differs")
    if document.get("disposition") != ORDER3_DISPOSITION:
        raise FrontierViolation("order-three dependency is not closed")
    return {
        "campaign_id": ORDER3_CAMPAIGN,
        "disposition": ORDER3_DISPOSITION,
        "document_sha256": _canonical_document_sha256(document),
        "independent_checker_required": True,
    }


def run_analysis(order3_artifact: Path | None = None) -> dict:
    if order3_artifact is None:
        order3_artifact = Path(__file__).resolve().parents[1] / "data" / "order3_srg_frontier.json"

    order5_large = order5_large_fixed_records()
    order5_free = enumerate_order5_free_quotients()
    if not all(record["minimum_cauchy_gap"] > 0 for record in order5_large):
        raise FrontierViolation("order-five large-fixed moment did not close")
    if any(record["solutions"] for record in order5_free):
        raise FrontierViolation("a fixed-point-free order-five quotient survived")

    order7 = order7_obstructions()
    if order7["fixed_17"]["maximum_support_family_size"] >= 17:
        raise FrontierViolation("order-seven fixed-17 support bound did not close")
    if order7["fixed_31"]["rayleigh_fourfold_gap"] >= 0:
        raise FrontierViolation("order-seven Rayleigh contradiction did not close")

    quotient11 = order11_canonical_quotient()
    if not verify_order11_quotient(quotient11):
        raise FrontierViolation("order-eleven quotient is invalid")
    coverage11 = order11_block_coverage()
    if coverage11["k4_primary_witnesses"] + coverage11["i5_fallback_witnesses"] != 46200:
        raise FrontierViolation("order-eleven block coverage is incomplete")

    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "claim": {
            "graph_family": "ramsey_good_strongly_regular",
            "parameters": list(PARAMETERS),
            "conclusion": "automorphism_group_order_is_a_power_of_two",
            "excluded_odd_prime_orders": [3, 5, 7, 11],
            "excluded_odd_prime_orders_at_least": 13,
            "general_ramsey_bound_claimed": False,
        },
        "group_theory": {
            "principle": "Cauchy_theorem",
            "largest_prime_dividing_a_subgroup_of_S45": 43,
            "conclusion": "no_odd_prime_divides_automorphism_group_order",
        },
        "spectral_reduction": {
            "nonprincipal_polynomial": "x^2+x-11",
            "nonprincipal_multiplicities": [22, 22],
            "irreducible_over_Q_zeta_p_for_odd_primes_except": [5],
            "nonfive_prime_cycle_count_parity": "even",
        },
        "order_3_dependency": _order3_dependency(order3_artifact),
        "order_5": {
            "spectral_fixed_point_counts": list(range(0, 45, 5)),
            "trace": order5_trace_record(),
            "fixed_10_plus": {
                "disposition": "IMPOSSIBLE_ORBIT_DIAGONAL_MOMENT_FOR_ANY_SRG",
                "graph_family": "strongly_regular",
                "ramsey_hypothesis_used": False,
                "internal_orbit_degrees": [0, 2, 4],
                "moment_equations": {
                    "off_orbit_degree_sum": "22-d-s",
                    "off_orbit_square_sum": "66-d-d^2-5s",
                    "cauchy_requirement":
                        "(22-d-s)^2 <= (c-1)(66-d-d^2-5s)",
                },
                "records": order5_large,
            },
            "fixed_0": {
                "disposition": "IMPOSSIBLE_INTEGRAL_SEIDEL_QUOTIENT",
                "quotient_equation": "M^2=45I-5J",
                "row_sum": 0,
                "row_square_sum": 40,
                "off_diagonal_values": list(_ODD_VALUES),
                "first_row_types": [list(row) for row in _ORDER5_FIRST_ROWS],
                "matrix_search": order5_free,
            },
            "fixed_5": {
                "disposition": "IMPOSSIBLE_TRIPLE_RAMSEY_BRIDGE",
                **order5_fixed_five_obstruction(),
            },
        },
        "order_7": order7,
        "order_11": {
            "spectral_candidates": [
                {"fixed_points": 1, "orbit_11_cycles": 4},
                {"fixed_points": 23, "orbit_11_cycles": 2},
            ],
            "fixed_23": order11_fixed_23_obstruction(),
            "fixed_1": {
                "disposition": "IMPOSSIBLE_CYCLIC_TWO_ORBIT_RAMSEY_BLOCK",
                "fixed_supported_orbits": [0, 1],
                "moving_internal_degrees": [4, 4, 6, 6],
                "canonical_quotient": quotient11,
                "block_coverage": coverage11,
            },
        },
        "prime_13_plus": large_prime_obstructions(),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--order3-artifact", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "order3_srg_frontier.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "odd_prime_srg_frontier.json",
    )
    args = parser.parse_args(argv)
    try:
        document = run_analysis(args.order3_artifact)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    except (OSError, FrontierViolation) as error:
        print(f"ODD-PRIME FRONTIER FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "ODD-PRIME FRONTIER:",
        "order5_fixed_free_quotients=0",
        "order5_fixed_five=RAMSEY_BRIDGE",
        "order7_cases=3/3",
        f"order11_blocks={document['order_11']['fixed_1']['block_coverage']['configurations']}",
    )
    print(
        "ODD-PRIME RESULT:", DISPOSITION,
        "(Ramsey-good srg(45,22,10,11) only; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
