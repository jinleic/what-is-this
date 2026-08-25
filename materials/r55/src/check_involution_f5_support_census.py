#!/usr/bin/env python3
"""Independent checker for the fixed-five balanced-support census.

The producer uses a Walsh closed-neighborhood enumeration with ternary
multiplicities.  This checker does not import it.  Instead it constructs the
two 16-by-16 incidence-moment matrices, inverts the odd matrix over the
rationals, and recursively enumerates every weak composition of ten even
blocks.  Nonnegative integral odd partners are recovered from the solved
linear system.  The C5 automorphism group is generated independently by
filtering all 120 point permutations.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_support_census"
DISPOSITION = "F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS"
PARAMETERS = (45, 22, 10, 11)
POINTS = tuple(range(5))
FULL_MASK = 31
EVEN_MASKS = tuple(mask for mask in range(32) if mask.bit_count() % 2 == 0)
ODD_MASKS = tuple(mask for mask in range(32) if mask.bit_count() % 2 == 1)
EVEN_INDEX = {mask: index for index, mask in enumerate(EVEN_MASKS)}
ODD_INDEX = {mask: index for index, mask in enumerate(ODD_MASKS)}
PAIRS = tuple(itertools.combinations(POINTS, 2))
ROTATION_EDGE_SUPPORTS = (3, 6, 12, 15, 17, 23, 24, 27, 29, 30)
ROTATION_NONEDGE_SUPPORTS = (1, 2, 4, 8, 11, 13, 16, 21, 22, 26)


class CheckViolation(RuntimeError):
    """The artifact differs from the independent exact reconstruction."""


def _features(mask: int) -> tuple[int, ...]:
    return (
        1,
        *(int(mask >> vertex & 1) for vertex in POINTS),
        *(int((mask >> left & 1) and (mask >> right & 1))
          for left, right in PAIRS),
    )


def _feature_matrix(mask_order: tuple[int, ...]) -> list[list[int]]:
    columns = [_features(mask) for mask in mask_order]
    return [[columns[column][row] for column in range(len(mask_order))]
            for row in range(16)]


def _inverse(matrix: list[list[int]]) -> list[list[Fraction]]:
    order = len(matrix)
    augmented = [
        [Fraction(value) for value in row]
        + [Fraction(int(left == right)) for right in range(order)]
        for left, row in enumerate(matrix)
    ]
    for column in range(order):
        pivot = next((row for row in range(column, order)
                      if augmented[row][column]), None)
        if pivot is None:
            raise CheckViolation("the odd incidence-moment matrix is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(order):
            if row == column or not augmented[row][column]:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                left - factor * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [row[order:] for row in augmented]


def _matvec(matrix, vector):
    return [sum(coefficient * value
                for coefficient, value in zip(row, vector))
            for row in matrix]


def _partner_affine_map() -> tuple[int, list[int], list[list[int]]]:
    even_matrix = _feature_matrix(EVEN_MASKS)
    odd_inverse = _inverse(_feature_matrix(ODD_MASKS))
    target = [20] + [10] * 5 + [5] * 10
    constant = _matvec(odd_inverse, target)
    coefficients = []
    for row in odd_inverse:
        coefficients.append([
            -sum(row[moment] * even_matrix[moment][column]
                 for moment in range(16))
            for column in range(16)
        ])
    denominator = 1
    for value in itertools.chain(constant, *coefficients):
        denominator = math.lcm(denominator, value.denominator)
    constant_integer = [int(value * denominator) for value in constant]
    coefficient_integer = [
        [int(value * denominator) for value in row]
        for row in coefficients
    ]
    return denominator, constant_integer, coefficient_integer


def _moments(edge_counts: tuple[int, ...],
             odd_counts: tuple[int, ...]) -> tuple[int, ...]:
    even_matrix = _feature_matrix(EVEN_MASKS)
    odd_matrix = _feature_matrix(ODD_MASKS)
    return tuple(
        sum(even_matrix[row][column] * edge_counts[column]
            for column in range(16))
        + sum(odd_matrix[row][column] * odd_counts[column]
              for column in range(16))
        for row in range(16)
    )


def _enumerate_designs() -> tuple[
        list[tuple[tuple[int, ...], tuple[int, ...]]], dict]:
    denominator, constant, coefficients = _partner_affine_map()
    edge_counts = [0] * 16
    partner_numerators = constant.copy()
    designs = []
    stats = {"nodes": 0, "leaves": 0}
    suffix_max = [
        [max(coefficients[row][column:])
         for row in range(16)]
        for column in range(16)
    ]

    def visit(column: int, remaining: int) -> None:
        stats["nodes"] += 1
        if column == 15:
            edge_counts[column] = remaining
            for row in range(16):
                partner_numerators[row] += coefficients[row][column] * remaining
            stats["leaves"] += 1
            if all(value >= 0 and value % denominator == 0
                   for value in partner_numerators):
                odd_counts = tuple(value // denominator
                                   for value in partner_numerators)
                edge_tuple = tuple(edge_counts)
                if _moments(edge_tuple, odd_counts) != (
                        20, *([10] * 5), *([5] * 10)):
                    raise CheckViolation("linear partner failed direct moments")
                designs.append((edge_tuple, odd_counts))
            for row in range(16):
                partner_numerators[row] -= coefficients[row][column] * remaining
            return

        for value in range(remaining + 1):
            edge_counts[column] = value
            if value:
                for row in range(16):
                    partner_numerators[row] += coefficients[row][column] * value
            rest = remaining - value
            possible = all(
                partner_numerators[row] + suffix_max[column + 1][row] * rest >= 0
                for row in range(16)
            )
            if possible:
                visit(column + 1, rest)
            if value:
                for row in range(16):
                    partner_numerators[row] -= coefficients[row][column] * value

    visit(0, 10)
    designs.sort()
    if len(designs) != len(set(designs)):
        raise CheckViolation("the matrix census produced duplicate designs")
    return designs, stats


def _c5_edges() -> set[tuple[int, int]]:
    return {tuple(sorted((vertex, (vertex + 1) % 5))) for vertex in POINTS}


def _c5_automorphisms() -> tuple[tuple[int, ...], ...]:
    edges = _c5_edges()
    automorphisms = []
    for permutation in itertools.permutations(POINTS):
        image = {tuple(sorted((permutation[left], permutation[right])))
                 for left, right in edges}
        if image == edges:
            automorphisms.append(permutation)
    if len(automorphisms) != 10:
        raise CheckViolation("the independently generated C5 group is not D5")
    return tuple(sorted(automorphisms))


AUTOMORPHISMS = _c5_automorphisms()


def _permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    image = 0
    for source, target in enumerate(permutation):
        if mask >> source & 1:
            image |= 1 << target
    return image


def _permute_counts(counts: tuple[int, ...], mask_order: tuple[int, ...],
                    index: dict[int, int],
                    permutation: tuple[int, ...]) -> tuple[int, ...]:
    image = [0] * len(mask_order)
    for mask, count in zip(mask_order, counts):
        image[index[_permute_mask(mask, permutation)]] = count
    return tuple(image)


def _edge_orbit(edge_counts: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted({
        _permute_counts(edge_counts, EVEN_MASKS, EVEN_INDEX, permutation)
        for permutation in AUTOMORPHISMS
    }))


def _support_counts(supports: tuple[int, ...], order: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(supports.count(mask) for mask in order)


def _design_line(edge_counts: tuple[int, ...], odd_counts: tuple[int, ...],
                 orbit_size: int | None = None) -> bytes:
    edge_text = "".join(map(str, edge_counts))
    odd_text = "".join(map(str, odd_counts))
    suffix = "" if orbit_size is None else f"|orbit={orbit_size}"
    return f"f5-support-v1|edge={edge_text}|odd={odd_text}{suffix}\n".encode()


def _digest(lines: list[bytes]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line)
    return digest.hexdigest()


def check_representative(record: dict) -> None:
    expected_keys = {
        "index", "edge_even_counts", "nonedge_odd_counts",
        "orbit_size", "stabilizer_size",
    }
    if type(record) is not dict or record.keys() != expected_keys:
        raise CheckViolation("representative fields differ")
    if any(type(record[key]) is not int
           for key in ("index", "orbit_size", "stabilizer_size")):
        raise CheckViolation("representative scalar types differ")
    edge = record["edge_even_counts"]
    odd = record["nonedge_odd_counts"]
    if (type(edge) is not list or type(odd) is not list
            or len(edge) != 16 or len(odd) != 16
            or any(type(value) is not int or value < 0
                   for value in edge + odd)):
        raise CheckViolation("representative multiplicities differ")
    edge_tuple, odd_tuple = tuple(edge), tuple(odd)
    if sum(edge_tuple) != 10 or sum(odd_tuple) != 10:
        raise CheckViolation("representative class sizes differ")
    if _moments(edge_tuple, odd_tuple) != (20, *([10] * 5), *([5] * 10)):
        raise CheckViolation("representative is not balanced")
    orbit = _edge_orbit(edge_tuple)
    if edge_tuple != orbit[0]:
        raise CheckViolation("representative is not D5-canonical")
    if record["orbit_size"] != len(orbit):
        raise CheckViolation("representative orbit size differs")
    if record["stabilizer_size"] * record["orbit_size"] != 10:
        raise CheckViolation("representative stabilizer size differs")


def _build_support_census() -> dict:
    designs, _stats = _enumerate_designs()
    if len(designs) != 7872:
        raise CheckViolation("the independently reconstructed support count differs")
    design_map = {edge: odd for edge, odd in designs}
    canonical = {}
    for edge, _odd in designs:
        orbit = _edge_orbit(edge)
        if any(image not in design_map for image in orbit):
            raise CheckViolation("the matrix census is not D5-invariant")
        canonical.setdefault(orbit[0], len(orbit))
    if len(canonical) != 844:
        raise CheckViolation("the independently reconstructed orbit count differs")

    representatives = []
    for index, edge in enumerate(sorted(canonical)):
        record = {
            "index": index,
            "edge_even_counts": list(edge),
            "nonedge_odd_counts": list(design_map[edge]),
            "orbit_size": canonical[edge],
            "stabilizer_size": 10 // canonical[edge],
        }
        check_representative(record)
        representatives.append(record)

    fixed_by = [
        sum(_permute_counts(edge, EVEN_MASKS, EVEN_INDEX, permutation) == edge
            for edge, _odd in designs)
        for permutation in AUTOMORPHISMS
    ]
    identity = next(index for index, permutation in enumerate(AUTOMORPHISMS)
                    if permutation == POINTS)
    nonidentity_rotation_counts = sorted(
        count for index, (permutation, count) in enumerate(zip(AUTOMORPHISMS, fixed_by))
        if index != identity and all(
            (permutation[(vertex + 1) % 5] - permutation[vertex]) % 5 == 1
            for vertex in POINTS
        )
    )
    reflection_counts = sorted(
        count for permutation, count in zip(AUTOMORPHISMS, fixed_by)
        if all((permutation[(vertex + 1) % 5] - permutation[vertex]) % 5 == 4
               for vertex in POINTS)
    )
    if nonidentity_rotation_counts != [2] * 4 or reflection_counts != [112] * 5:
        raise CheckViolation("the independent Burnside classes differ")
    if sum(fixed_by) // 10 != len(representatives):
        raise CheckViolation("independent Burnside averaging differs")

    control_edge = _support_counts(ROTATION_EDGE_SUPPORTS, EVEN_MASKS)
    control_odd = _support_counts(ROTATION_NONEDGE_SUPPORTS, ODD_MASKS)
    if (control_edge, control_odd) not in set(designs):
        raise CheckViolation("the rotation support control disappeared")
    control_canonical = min(_edge_orbit(control_edge))
    canonical_order = sorted(canonical)

    raw_lines = [_design_line(edge, odd) for edge, odd in designs]
    representative_lines = [
        _design_line(tuple(record["edge_even_counts"]),
                     tuple(record["nonedge_odd_counts"]),
                     record["orbit_size"])
        for record in representatives
    ]
    permutation_lines = [
        ("f5-support-d5-v1|" + ",".join(map(str, permutation)) + "\n").encode()
        for permutation in AUTOMORPHISMS
    ]
    histogram = {str(size): sum(value == size for value in canonical.values())
                 for size in (1, 5, 10)}
    return {
        "fixed_graph": "C5",
        "edge_mask_order": list(EVEN_MASKS),
        "nonedge_mask_order": list(ODD_MASKS),
        "all_even_multisets_before_walsh_bound": math.comb(25, 15),
        "walsh_ternary_candidates": sum(
            math.comb(16, doubled)
            * math.comb(16 - doubled, 10 - 2 * doubled)
            for doubled in range(6)
        ),
        "balanced_labelled_designs": len(designs),
        "dihedral_group_order": len(AUTOMORPHISMS),
        "dihedral_orbits": len(representatives),
        "orbit_size_histogram": histogram,
        "burnside_fixed_designs": {
            "identity": fixed_by[identity],
            "nonidentity_rotations_each": nonidentity_rotation_counts[0],
            "reflections_each": reflection_counts[0],
        },
        "raw_designs_sha256": _digest(raw_lines),
        "representatives_sha256": _digest(representative_lines),
        "dihedral_action_sha256": _digest(permutation_lines),
        "representatives": representatives,
        "positive_control": {
            "construction": "five rotations of one support of each size 1,2,3,4",
            "edge_support_masks": list(ROTATION_EDGE_SUPPORTS),
            "nonedge_support_masks": list(ROTATION_NONEDGE_SUPPORTS),
            "canonical_representative_index": canonical_order.index(control_canonical),
            "dihedral_orbit_size": len(_edge_orbit(control_edge)),
        },
        "disposition": "EXACT_BALANCED_SUPPORT_MULTISET_CENSUS_COMPLETE",
    }


def _signed_completion_model() -> dict:
    return {
        "fixed_incidence_sign_matrix": "R_vi=1-2*N_vi",
        "fixed_seidel_matrix": "Q_F=J-I-2*A(C5)",
        "invariant_diagonal": "W_ii=-1 on edge pairs and +1 on nonedge pairs",
        "invariant_off_diagonal_values": [-2, 0, 2],
        "invariant_cross_equation": "Q_F R+R W=-J",
        "invariant_square": "W^2+2R^T R=45I-2J",
        "invariant_row_nonzeros": 8,
        "anti_invariant_diagonal": "T_ii=+1 on edge pairs and -1 on nonedge pairs",
        "anti_invariant_off_diagonal_values": [-2, 0, 2],
        "anti_invariant_square": "T^2=45I",
        "anti_invariant_row_nonzeros": 11,
        "off_diagonal_compatibility": "W_ij^2+T_ij^2=4",
        "relation_recovery": "b_ij=1-W_ij/2, c_ij=-T_ij/2, p=(b+c)/2, q=(b-c)/2",
        "status": "OPEN_SIGNED_COMPLETION_FOR_844_BALANCED_SUPPORT_ORBITS",
        "ramsey_constraints_included": False,
        "ramsey_bound_claimed": False,
    }


def _build_expected_document() -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "claim": {
            "graph_family": "strongly_regular",
            "parameters": list(PARAMETERS),
            "involution_fixed_vertices": 5,
            "fixed_graph": "C5",
            "balanced_support_multiset_census_complete": True,
            "signed_completion_proved": False,
            "ramsey_hypothesis_used": False,
            "general_ramsey_bound_claimed": False,
        },
        "walsh_reduction": {
            "even_support_sign_vectors": "chi_E in {+1,-1}^5 with product +1",
            "odd_support_complement": "z_T counts odd support V\\T",
            "singleton_characters": "hat(x)=hat(z)",
            "pair_characters": "hat(x)=-hat(z)",
            "partner_formula": "z_T=(5-sum_{|E xor T| in {0,4}} x_E)/2",
            "valid_closed_neighborhood_sums": [1, 3, 5],
            "valid_multiplicities": [0, 1, 2],
            "transform_is_involutive": True,
        },
        "support_census": _build_support_census(),
        "signed_completion": _signed_completion_model(),
        "disposition": DISPOSITION,
    }


@functools.cache
def _expected_document_json() -> str:
    return json.dumps(_build_expected_document(), sort_keys=True,
                      separators=(",", ":"))


def _expected_document() -> dict:
    return json.loads(_expected_document_json())


def _exact_tree_equal(actual, expected) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return (actual.keys() == expected.keys()
                and all(_exact_tree_equal(actual[key], value)
                        for key, value in expected.items()))
    if isinstance(expected, list):
        return (len(actual) == len(expected)
                and all(_exact_tree_equal(left, right)
                        for left, right in zip(actual, expected)))
    return actual == expected


def _unique_object(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise CheckViolation(f"duplicate JSON member: {key}")
        output[key] = value
    return output


def load_document(text: str) -> dict:
    """Parse strict JSON, rejecting duplicate names and nonstandard constants."""
    def reject_constant(value):
        raise CheckViolation(f"nonstandard JSON constant: {value}")

    return json.loads(text, object_pairs_hook=_unique_object,
                      parse_constant=reject_constant)


def verify_document(document: dict) -> dict:
    expected = _expected_document()
    if not _exact_tree_equal(document, expected):
        raise CheckViolation("artifact differs from independent matrix reconstruction")
    return {
        "balanced_labelled_designs": 7872,
        "dihedral_orbits": 844,
        "signed_completion_status": "OPEN_SIGNED_COMPLETION_FOR_844_BALANCED_SUPPORT_ORBITS",
        "ramsey_bound_claimed": False,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    args = parser.parse_args(argv)
    try:
        document = load_document(args.artifact.read_text())
        summary = verify_document(document)
    except (OSError, json.JSONDecodeError, CheckViolation) as error:
        print(f"F5 SUPPORT CHECK FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "F5 SUPPORT EVIDENCE VERIFIED:", DISPOSITION,
        f"labelled={summary['balanced_labelled_designs']}",
        f"D5_orbits={summary['dihedral_orbits']}",
        ("balanced-support upper bound only; signed/Ramsey completion open; "
         "no R(5,5) bound claimed"),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
