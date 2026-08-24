#!/usr/bin/env python3
"""Independent checker for the odd-prime SRG frontier artifact.

This checker imports no symbol from ``odd_prime_srg_frontier``.  It composes the
separate order-three checker, reconstructs the order-five matrix search with a
row-multiset DFS rather than the producer's moment-completion recursion, derives
the order-seven arithmetic afresh, and resweeps all 46,200 order-eleven cyclic
blocks with separate graph construction and clique kernels.
"""

from __future__ import annotations

import argparse
import collections
import functools
import hashlib
import itertools
import json
import sys
from pathlib import Path

import check_order3_srg_frontier as order3_check

SCHEMA_VERSION = 1
CAMPAIGN_ID = "odd_prime_srg_frontier"
DISPOSITION = "AUTOMORPHISM_GROUP_IS_A_2_GROUP_FOR_RAMSEY_SRG"
ORDER3_CAMPAIGN = "order3_srg_frontier"
ORDER3_DISPOSITION = "ORDER3_AUTOMORPHISMS_EXCLUDED_FOR_RAMSEY_SRG"
PARAMETERS = (45, 22, 10, 11)
VALUES = (-5, -3, -1, 1, 3, 5)


class CheckViolation(RuntimeError):
    """The artifact differs from the independent exact reconstruction."""


def _canonical_sha256(document: dict) -> str:
    payload = json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _order5_row_types() -> tuple[tuple[int, ...], ...]:
    rows = []
    for row in itertools.combinations_with_replacement(VALUES, 8):
        if sum(row) == 0 and sum(value * value for value in row) == 40:
            rows.append(row)
    return tuple(sorted(rows, key=lambda row: (sum(abs(value) == 5 for value in row),
                                               -max(row))))


def _multiset_permutations(counter: collections.Counter, length: int):
    values = tuple(sorted(counter))
    row = [0] * length

    def visit(index):
        if index == length:
            yield tuple(row)
            return
        for value in values:
            if counter[value] == 0:
                continue
            counter[value] -= 1
            row[index] = value
            yield from visit(index + 1)
            counter[value] += 1

    yield from visit(0)


def _independent_free_quotient_search(first_row, row_types):
    order = 9
    matrix = [[0] * order for _ in range(order)]
    for column, value in enumerate(first_row, 1):
        matrix[0][column] = matrix[column][0] = value
    tested = collections.Counter()
    surviving = collections.Counter()
    serializations = collections.defaultdict(list)

    def visit(row):
        if row == order:
            return 1
        prefix = matrix[row][:row]
        candidates = []
        prefix_counts = collections.Counter(prefix)
        for row_type in row_types:
            remaining = collections.Counter(row_type)
            if any(prefix_counts[value] > remaining[value] for value in prefix_counts):
                continue
            remaining.subtract(prefix_counts)
            if any(count < 0 for count in remaining.values()):
                continue
            remaining += collections.Counter()
            candidates.extend(_multiset_permutations(remaining, order - row - 1))
        candidates.sort()
        solutions = 0
        for tail in candidates:
            tested[row] += 1
            for column, value in enumerate(tail, row + 1):
                matrix[row][column] = matrix[column][row] = value
            valid = True
            for previous in range(row):
                dot = sum(matrix[previous][index] * matrix[index][row]
                          for index in range(order))
                if dot != -5:
                    valid = False
                    break
            if not valid:
                continue
            surviving[row] += 1
            flattened = []
            for previous in range(row + 1):
                flattened.extend(matrix[previous][column]
                                 for column in range(previous + 1, order))
            serializations[row].append(",".join(map(str, flattened)) + "\n")
            solutions += visit(row + 1)
        return solutions

    solutions = visit(1)
    digests = {}
    for row in sorted(tested):
        digest = hashlib.sha256()
        for serialization in sorted(serializations[row]):
            digest.update(serialization.encode())
        digests[str(row)] = digest.hexdigest()
    return {
        "first_row": list(first_row),
        "tested_by_row": {str(row): tested[row] for row in sorted(tested)},
        "surviving_by_row": {str(row): surviving[row] for row in sorted(surviving)},
        "survivor_sha256_by_row": digests,
        "solutions": solutions,
    }


def _order5_free_searches():
    row_types = _order5_row_types()
    if len(row_types) != 3:
        raise CheckViolation("order-five row classification differs")
    return [_independent_free_quotient_search(first, row_types) for first in row_types]


def _order5_large_records():
    records = []
    for fixed_points in range(10, 45, 5):
        cycles = (45 - fixed_points) // 5
        degree_records = []
        for internal_degree in (0, 2, 4):
            maximum_support = min(
                fixed_points,
                (66 - internal_degree - internal_degree * internal_degree) // 5,
            )
            gaps = []
            for support in range(maximum_support + 1):
                degree_sum = 22 - internal_degree - support
                square_sum = 66 - internal_degree - internal_degree ** 2
                square_sum -= 5 * support
                gaps.append(degree_sum ** 2 - (cycles - 1) * square_sum)
            degree_records.append({
                "internal_degree": internal_degree,
                "fixed_support_values": list(range(maximum_support + 1)),
                "cauchy_gaps": gaps,
                "minimum_cauchy_gap": min(gaps),
            })
        records.append({
            "fixed_points": fixed_points,
            "orbit_5_cycles": cycles,
            "internal_degree_records": degree_records,
            "minimum_cauchy_gap": min(
                item["minimum_cauchy_gap"] for item in degree_records),
        })
    return records


def _fixed_five_record():
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

def _order5_large_minimum_gap_from_parameters():
    vertices, degree, adjacent_common, nonadjacent_common = PARAMETERS
    if adjacent_common - nonadjacent_common != -1:
        raise CheckViolation("SRG polynomial is not x^2+x-11")
    identity_coefficient = degree - nonadjacent_common
    quotient_diagonal = identity_coefficient + 5 * nonadjacent_common
    minima = []
    for fixed_points in range(10, vertices, 5):
        cycles = (vertices - fixed_points) // 5
        case_gaps = []
        for internal_degree in range(0, 5, 2):
            maximum_support = min(
                fixed_points,
                (quotient_diagonal - internal_degree ** 2
                 - internal_degree) // 5,
            )
            for support in range(maximum_support + 1):
                degree_sum = degree - internal_degree - support
                square_sum = quotient_diagonal - internal_degree ** 2
                square_sum -= internal_degree + 5 * support
                case_gaps.append(
                    degree_sum ** 2 - (cycles - 1) * square_sum)
        minima.append(min(case_gaps))
    if any(gap <= 0 for gap in minima):
        raise CheckViolation("parameter-derived order-five moment has a survivor")
    return min(minima)


def _fixed5_block_design_audit():
    prime = 5
    fixed_points = 5
    moving_orbits = (PARAMETERS[0] - fixed_points) // prime
    points = range(fixed_points)
    fixed_degree = PARAMETERS[1] % prime
    replication = (PARAMETERS[1] - fixed_degree) // prime
    adjacent_pair_multiplicity = PARAMETERS[2] // prime
    nonadjacent_pair_multiplicity = (PARAMETERS[3] - 1) // prime
    if (fixed_degree != 2
            or PARAMETERS[2] % prime
            or (PARAMETERS[3] - 1) % prime
            or adjacent_pair_multiplicity != nonadjacent_pair_multiplicity):
        raise CheckViolation("fixed-five incidence parameters do not reduce")
    pair_multiplicity = adjacent_pair_multiplicity
    total_memberships = fixed_points * replication
    total_pair_incidences = fixed_points * (fixed_points - 1) // 2
    total_pair_incidences *= pair_multiplicity
    masks_by_size = {
        size: tuple(mask for mask in range(1 << fixed_points)
                    if mask.bit_count() == size)
        for size in range(fixed_points + 1)
    }
    size_patterns = tuple(
        sizes for sizes in itertools.combinations_with_replacement(
            range(fixed_points + 1), moving_orbits)
        if sum(sizes) == total_memberships
        and sum(size * (size - 1) // 2 for size in sizes)
        == total_pair_incidences
    )
    designs = set()
    for sizes in size_patterns:
        multiplicities = collections.Counter(sizes)
        choices = [
            tuple(itertools.combinations_with_replacement(
                masks_by_size[size], count))
            for size, count in sorted(multiplicities.items())
        ]
        for groups in itertools.product(*choices):
            blocks = tuple(sorted(itertools.chain.from_iterable(groups)))
            if any(sum(mask >> point & 1 for mask in blocks) != replication
                   for point in points):
                continue
            if any(sum((mask >> left & 1) and (mask >> right & 1)
                       for mask in blocks) != pair_multiplicity
                   for left, right in itertools.combinations(points, 2)):
                continue
            designs.add(blocks)

    def fixed_edge(left, right):
        return (left - right) % 5 in (1, 4)

    def has_ramsey_bridge(blocks):
        for triple in itertools.combinations(points, 3):
            triple_mask = sum(1 << point for point in triple)
            multiplicity = sum(
                mask & triple_mask == triple_mask for mask in blocks)
            if multiplicity == 0:
                for left, right in itertools.combinations(triple, 2):
                    if not fixed_edge(left, right):
                        continue
                    omitted = next(point for point in triple
                                   if point not in (left, right))
                    eligible = sum(
                        (mask >> left & 1) and (mask >> right & 1)
                        and not (mask >> omitted & 1)
                        for mask in blocks
                    )
                    if eligible >= 2:
                        return True
            elif multiplicity >= 2:
                for left, right in itertools.combinations(triple, 2):
                    if fixed_edge(left, right):
                        continue
                    included = next(point for point in triple
                                    if point not in (left, right))
                    eligible = sum(
                        (mask >> included & 1)
                        and not (mask >> left & 1)
                        and not (mask >> right & 1)
                        for mask in blocks
                    )
                    if eligible >= 2:
                        return True
        return False

    bridged = sum(has_ramsey_bridge(blocks) for blocks in designs)
    if not designs or bridged != len(designs):
        raise CheckViolation("fixed-five incidence audit found an unbridged design")
    return len(designs), bridged


def _order7_fixed17_support_capacity():
    universe = (1 << 4) - 1
    supports = []
    for mask in range(1 << 4):
        fixed_degree = PARAMETERS[1] - 7 * mask.bit_count()
        if 0 <= fixed_degree <= 16:
            supports.append(mask)

    def compatible(left, right):
        return ((left & right).bit_count() <= 1
                and ((universe ^ left) & (universe ^ right)).bit_count() <= 1)

    for size in range(len(supports), -1, -1):
        for family in itertools.combinations(supports, size):
            if all(compatible(left, right)
                   for left, right in itertools.combinations(family, 2)):
                return size
    raise CheckViolation("order-seven support family enumeration failed")


def _order7_fixed31_rayleigh_audit():
    prime = 7
    fixed_points = 31
    cycles = 2
    trace = (prime - 1) * cycles // 2
    allowed_internal = tuple(
        value for value in range(0, prime, 2)
        if value not in (0, prime - 1)
    )
    degree_pairs = tuple(
        pair for pair in itertools.combinations_with_replacement(
            allowed_internal, cycles)
        if sum(pair) == trace
    )
    if degree_pairs != ((2, 4),):
        raise CheckViolation("order-seven moving degrees differ")
    incidence_cases = []
    for both in range(2):
        for neither in range(2):
            numerator = 7 + neither - both
            if numerator % 2:
                continue
            cross_degree = numerator // 2
            supports = (
                PARAMETERS[1] - degree_pairs[0][0] - cross_degree,
                PARAMETERS[1] - degree_pairs[0][1] - cross_degree,
            )
            if sum(supports) == fixed_points - neither + both:
                incidence_cases.append((both, neither, cross_degree, supports))
    if incidence_cases != [(0, 1, 4, (16, 14)), (1, 0, 3, (17, 15))]:
        raise CheckViolation("order-seven incidence cases differ")
    group_sizes = (16, 14)
    norm = sum(group_sizes)
    vector_sum = group_sizes[0] - group_sizes[1]
    projection = (group_sizes[0], -group_sizes[1])
    identity_coefficient = PARAMETERS[1] - PARAMETERS[3]
    rhs = identity_coefficient * norm + PARAMETERS[3] * vector_sum ** 2
    rhs -= prime * sum(value ** 2 for value in projection)
    gap = 4 * rhs + norm
    if gap >= 0:
        raise CheckViolation("order-seven Rayleigh audit did not contradict")
    return gap


def _order11_fixed23_group_lower_bound():
    prime = 11
    fixed_points = 23
    cycles = 2
    support_sizes = [
        support for support in range(cycles + 1)
        if 0 <= PARAMETERS[1] - prime * support <= fixed_points - 1
    ]
    if support_sizes != [0, 1, 2]:
        raise CheckViolation("order-eleven fixed support sizes differ")
    maximum_shared_orbits = max(PARAMETERS[2], PARAMETERS[3]) // prime
    if maximum_shared_orbits != 1 or cycles <= maximum_shared_orbits:
        raise CheckViolation("order-eleven support intersection bound differs")
    maximum_both = maximum_neither = 1
    single_support_vertices = fixed_points - maximum_both - maximum_neither
    lower_bound = (single_support_vertices + cycles - 1) // cycles
    if prime <= PARAMETERS[2] or lower_bound < 5:
        raise CheckViolation("order-eleven support pigeonhole is insufficient")
    return lower_bound


def _is_prime(value):
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return value >= 2


def _prime13_plus_audit():
    excluded = []
    for prime in range(13, 44):
        if not _is_prime(prime):
            continue
        cycle_counts = [
            cycles for cycles in range(2, 45 // prime + 1, 2)
        ]
        if not cycle_counts:
            excluded.append(prime)
            continue
        if cycle_counts != [2]:
            raise CheckViolation("unexpected large-prime cycle count")
        fixed_points = 45 - 2 * prime
        fixed_degrees = [
            PARAMETERS[1] - prime * support
            for support in range(3)
            if 0 <= PARAMETERS[1] - prime * support <= fixed_points - 1
        ]
        if len(fixed_degrees) != 1:
            raise CheckViolation("large-prime fixed degree is not forced")
        if fixed_points * fixed_degrees[0] % 2 == 0:
            raise CheckViolation("large-prime handshake did not contradict")
        excluded.append(prime)
    return excluded


@functools.cache
def independent_proof_summary():
    designs, bridged = _fixed5_block_design_audit()
    if 3 * (PARAMETERS[1] % 7) % 2 == 0:
        raise CheckViolation("order-seven fixed-three handshake did not contradict")
    fixed17_capacity = _order7_fixed17_support_capacity()
    if fixed17_capacity >= 17:
        raise CheckViolation("order-seven fixed support family fits")
    _order7_fixed31_rayleigh_audit()
    large_primes = _prime13_plus_audit()
    relevant_odd_primes = {
        prime for prime in range(3, 44) if _is_prime(prime)
    }
    if {3, 5, 7, 11, *large_primes} != relevant_odd_primes:
        raise CheckViolation("odd-prime coverage is incomplete")
    return {
        "order5_large_minimum_gap":
            _order5_large_minimum_gap_from_parameters(),
        "order5_fixed5_block_designs": designs,
        "order5_fixed5_bridged_designs": bridged,
        "order7_fixed17_maximum_support_family": fixed17_capacity,
        "order11_fixed23_support_group_lower_bound":
            _order11_fixed23_group_lower_bound(),
        "prime13_plus_excluded": large_primes,
    }


def _order7_record():
    norm = 30
    vector_sum = 2
    projection = (16, -14)
    rhs = 11 * (norm + vector_sum * vector_sum) - 7 * sum(x * x for x in projection)
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
            "test_vector_norm_square": norm,
            "test_vector_sum": vector_sum,
            "incidence_projection": list(projection),
            "rhs_quadratic": rhs,
            "rayleigh_fourfold_gap": 4 * rhs + norm,
            "disposition": "IMPOSSIBLE_FIXED_BLOCK_RAYLEIGH",
        },
    }


def _quotient_valid(matrix):
    support = (1, 1, 0, 0)
    if any(matrix[i][j] != matrix[j][i] for i in range(4) for j in range(4)):
        return False
    if [matrix[i][i] for i in range(4)] not in ([4, 4, 6, 6], [6, 6, 4, 4]):
        return False
    if any(sum(matrix[i]) + support[i] != 22 for i in range(4)):
        return False
    if any(sum(support[i] * matrix[i][j] for i in range(4)) + support[j] != 11
           for j in range(4)):
        return False
    for i in range(4):
        for j in range(4):
            lhs = sum(matrix[i][k] * matrix[k][j] for k in range(4))
            lhs += matrix[i][j] + 11 * support[i] * support[j]
            if lhs != (132 if i == j else 121):
                return False
    return True


def _derive_order11_quotients():
    records = []
    for supported_degree in (4, 6):
        unsupported_degree = 10 - supported_degree
        supported_cross = 10 - supported_degree
        for cross_02 in range(12):
            cross_03 = 11 - cross_02
            matrix = [
                [supported_degree, supported_cross, cross_02, cross_03],
                [supported_cross, supported_degree, cross_03, cross_02],
                [cross_02, cross_03, unsupported_degree, 11 - unsupported_degree],
                [cross_03, cross_02, 11 - unsupported_degree, unsupported_degree],
            ]
            if _quotient_valid(matrix):
                records.append(matrix)
    return records


def _cyclic_graph(left_mask, right_mask, cross_mask):
    prime = 11
    rows = [0] * 22
    for left in range(22):
        orbit_left, position_left = divmod(left, prime)
        for right in range(left + 1, 22):
            orbit_right, position_right = divmod(right, prime)
            if orbit_left == orbit_right:
                difference = (position_right - position_left) % prime
                pair = min(difference, prime - difference) - 1
                internal = left_mask if orbit_left == 0 else right_mask
                edge = bool(internal >> pair & 1)
            else:
                difference = (position_right - position_left) % prime
                edge = bool(cross_mask >> difference & 1)
            if edge:
                rows[left] |= 1 << right
                rows[right] |= 1 << left
    return tuple(rows)


def _lex_clique_by_loops(rows, order):
    size = len(rows)
    masks_after = [((1 << size) - 1) ^ ((1 << (index + 1)) - 1)
                   for index in range(size)]
    if order == 4:
        for a in range(size):
            abits = rows[a] & masks_after[a]
            while abits:
                bbit = abits & -abits
                b = bbit.bit_length() - 1
                abits ^= bbit
                common = rows[a] & rows[b] & masks_after[b]
                cbits = common
                while cbits:
                    cbit = cbits & -cbits
                    c = cbit.bit_length() - 1
                    cbits ^= cbit
                    dbits = common & rows[c] & masks_after[c]
                    if dbits:
                        return a, b, c, (dbits & -dbits).bit_length() - 1
        return None
    if order == 5:
        for a in range(size):
            abits = rows[a] & masks_after[a]
            while abits:
                bbit = abits & -abits
                b = bbit.bit_length() - 1
                abits ^= bbit
                common2 = rows[a] & rows[b] & masks_after[b]
                cbits = common2
                while cbits:
                    cbit = cbits & -cbits
                    c = cbit.bit_length() - 1
                    cbits ^= cbit
                    common3 = common2 & rows[c] & masks_after[c]
                    dbits = common3
                    while dbits:
                        dbit = dbits & -dbits
                        d = dbit.bit_length() - 1
                        dbits ^= dbit
                        ebits = common3 & rows[d] & masks_after[d]
                        if ebits:
                            return a, b, c, d, (ebits & -ebits).bit_length() - 1
        return None
    raise CheckViolation("unsupported independent clique order")


def _complement(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ (1 << vertex) ^ row for vertex, row in enumerate(rows))


def _order11_coverage():
    internal_masks = tuple(sum(1 << index for index in pair)
                           for pair in itertools.combinations(range(5), 2))
    cross_masks = tuple(sum(1 << index for index in selected)
                        for selected in itertools.combinations(range(11), 6))
    digest = hashlib.sha256()
    primary = fallback = both = 0
    examples = {}
    for left in internal_masks:
        for right in internal_masks:
            for cross in cross_masks:
                graph = _cyclic_graph(left, right, cross)
                clique = _lex_clique_by_loops(graph, 4)
                independent = _lex_clique_by_loops(_complement(graph), 5)
                if clique is not None:
                    kind, witness = "K4", clique
                    primary += 1
                    if independent is not None:
                        both += 1
                elif independent is not None:
                    kind, witness = "I5", independent
                    fallback += 1
                else:
                    raise CheckViolation("independent order-eleven resweep found an escape")
                examples.setdefault(kind, {
                    "internal_pair_masks": [left, right],
                    "cross_mask": cross,
                    "vertices": list(witness),
                })
                digest.update(
                    f"{left:02x}:{right:02x}:{cross:03x}:{kind}:"
                    f"{','.join(map(str, witness))}\n".encode()
                )
    return {
        "internal_connection_sets_per_orbit": 10,
        "cross_connection_sets": 462,
        "configurations": len(internal_masks) ** 2 * len(cross_masks),
        "k4_primary_witnesses": primary,
        "i5_fallback_witnesses": fallback,
        "configurations_with_both": both,
        "coverage_sha256": digest.hexdigest(),
        "examples": examples,
    }


def _expected_document(order3_document):
    order3_check.verify_document(order3_document)
    if order3_document.get("campaign_id") != ORDER3_CAMPAIGN:
        raise CheckViolation("order-three dependency campaign differs")
    if order3_document.get("disposition") != ORDER3_DISPOSITION:
        raise CheckViolation("order-three dependency disposition differs")

    row_types = _order5_row_types()
    free_searches = _order5_free_searches()
    if any(record["solutions"] for record in free_searches):
        raise CheckViolation("independent order-five quotient search found a solution")
    large5 = _order5_large_records()
    if not all(record["minimum_cauchy_gap"] > 0 for record in large5):
        raise CheckViolation("independent order-five moment did not close")

    quotient_candidates = _derive_order11_quotients()
    if len(quotient_candidates) != 2:
        raise CheckViolation("order-eleven quotient count differs")
    quotient11 = min(quotient_candidates)
    coverage11 = _order11_coverage()

    order7 = _order7_record()
    fixed5 = _fixed_five_record()
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
        "order_3_dependency": {
            "campaign_id": ORDER3_CAMPAIGN,
            "disposition": ORDER3_DISPOSITION,
            "document_sha256": _canonical_sha256(order3_document),
            "independent_checker_required": True,
        },
        "order_5": {
            "spectral_fixed_point_counts": list(range(0, 45, 5)),
            "trace": {
                "quotient_trace": "2*c",
                "internal_degree_counts_equation": "n0=n4",
                "ramsey_allowed_internal_degrees": [2],
                "moving_orbit_graph": "C5",
            },
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
                "records": large5,
            },
            "fixed_0": {
                "disposition": "IMPOSSIBLE_INTEGRAL_SEIDEL_QUOTIENT",
                "quotient_equation": "M^2=45I-5J",
                "row_sum": 0,
                "row_square_sum": 40,
                "off_diagonal_values": list(VALUES),
                "first_row_types": [list(row) for row in row_types],
                "matrix_search": free_searches,
            },
            "fixed_5": {
                "disposition": "IMPOSSIBLE_TRIPLE_RAMSEY_BRIDGE",
                **fixed5,
            },
        },
        "order_7": order7,
        "order_11": {
            "spectral_candidates": [
                {"fixed_points": 1, "orbit_11_cycles": 4},
                {"fixed_points": 23, "orbit_11_cycles": 2},
            ],
            "fixed_23": {
                "orbit_11_cycles": 2,
                "fixed_support_sizes_per_vertex": [0, 1, 2],
                "at_most_one_both_support_vertex": True,
                "at_most_one_neither_support_vertex": True,
                "single_support_vertices_lower_bound": 21,
                "largest_single_support_group_lower_bound": 11,
                "forced_independent_set_lower_bound": 11,
                "disposition": "IMPOSSIBLE_FIXED_INDEPENDENT_FIVE",
            },
            "fixed_1": {
                "disposition": "IMPOSSIBLE_CYCLIC_TWO_ORBIT_RAMSEY_BLOCK",
                "fixed_supported_orbits": [0, 1],
                "moving_internal_degrees": [4, 4, 6, 6],
                "canonical_quotient": quotient11,
                "block_coverage": coverage11,
            },
        },
        "prime_13_plus": {
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
        },
    }


def verify_document(document: dict, order3_document: dict) -> dict:
    independent_proof_summary()
    expected = _expected_document(order3_document)
    if document != expected:
        raise CheckViolation("artifact differs from independent exact reconstruction")
    return {
        "excluded_odd_primes": "all",
        "order5_free_quotient_solutions": 0,
        "order11_configurations": expected["order_11"]["fixed_1"]["block_coverage"]["configurations"],
        "automorphism_group": "2-group",
    }


def main(argv=None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument(
        "--order3-artifact", type=Path,
        default=root / "data" / "order3_srg_frontier.json",
    )
    args = parser.parse_args(argv)
    try:
        document = json.loads(args.artifact.read_text())
        order3_document = json.loads(args.order3_artifact.read_text())
        summary = verify_document(document, order3_document)
    except (OSError, json.JSONDecodeError, CheckViolation,
            order3_check.CheckViolation) as error:
        print(f"ODD-PRIME CHECK FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "ODD-PRIME EVIDENCE VERIFIED:", DISPOSITION,
        f"order11_blocks={summary['order11_configurations']}",
        "odd_primes=ALL",
        "(Ramsey SRG lane only; no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
