#!/usr/bin/env python3
"""Exact balanced-support census for the fixed-five involution branch.

Let an involution of an ``srg(45,22,10,11)`` fix five vertices.  The fixed
subgraph is ``C5`` and the other vertices form ten internally adjacent and
ten internally nonadjacent transposition orbits.  Their fixed-neighbor
supports are even and odd subsets of the five fixed points respectively.
Across all twenty supports, every point has replication ten and every pair
has multiplicity five.

This module enumerates the support multisets satisfying those necessary
incidence moments exactly.  Signed adjacency completion can eliminate some or
all of them.  No Ramsey hypothesis is used and no Ramsey-number bound changes.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_support_census"
DISPOSITION = "F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS"
PARAMETERS = (45, 22, 10, 11)
POINTS = tuple(range(5))
FULL_MASK = (1 << len(POINTS)) - 1
EVEN_MASKS = tuple(mask for mask in range(1 << len(POINTS))
                   if mask.bit_count() % 2 == 0)
ODD_MASKS = tuple(mask for mask in range(1 << len(POINTS))
                  if mask.bit_count() % 2 == 1)
EVEN_INDEX = {mask: index for index, mask in enumerate(EVEN_MASKS)}
ODD_INDEX = {mask: index for index, mask in enumerate(ODD_MASKS)}
ROTATION_EDGE_SUPPORTS = (3, 6, 12, 15, 17, 23, 24, 27, 29, 30)
ROTATION_NONEDGE_SUPPORTS = (1, 2, 4, 8, 11, 13, 16, 21, 22, 26)


class CensusViolation(RuntimeError):
    """An exact support invariant or finite enumeration failed."""


def _dihedral_permutations() -> tuple[tuple[int, ...], ...]:
    rotations = [tuple((vertex + shift) % 5 for vertex in POINTS)
                 for shift in range(5)]
    reflections = [tuple((shift - vertex) % 5 for vertex in POINTS)
                   for shift in range(5)]
    permutations = tuple(rotations + reflections)
    if len(set(permutations)) != 10:
        raise CensusViolation("the C5 dihedral action is not faithful")
    return permutations


DIHEDRAL_PERMUTATIONS = _dihedral_permutations()
_CLOSED_WALSH_NEIGHBORHOODS = tuple(
    sum(1 << EVEN_INDEX[source] for source in EVEN_MASKS
        if (source ^ target).bit_count() in (0, 4))
    for target in EVEN_MASKS
)


def counts_from_supports(supports: tuple[int, ...],
                         mask_order: tuple[int, ...]) -> tuple[int, ...]:
    """Convert a sorted or unsorted support multiset to multiplicities."""
    index = {mask: position for position, mask in enumerate(mask_order)}
    counts = [0] * len(mask_order)
    for support in supports:
        if support not in index:
            raise CensusViolation(f"support mask {support} has the wrong parity")
        counts[index[support]] += 1
    return tuple(counts)


def supports_from_counts(counts: tuple[int, ...] | list[int],
                         mask_order: tuple[int, ...]) -> tuple[int, ...]:
    """Expand multiplicities in the declared mask order."""
    if len(counts) != len(mask_order):
        raise CensusViolation("support-count vector has the wrong length")
    supports = []
    for mask, count in zip(mask_order, counts, strict=True):
        if type(count) is not int or count < 0:
            raise CensusViolation("support multiplicities must be nonnegative integers")
        supports.extend([mask] * count)
    return tuple(supports)


def _point_pair_moments(edge_counts: tuple[int, ...],
                        odd_counts: tuple[int, ...]) -> tuple[list[int], list[int]]:
    supports = (supports_from_counts(edge_counts, EVEN_MASKS)
                + supports_from_counts(odd_counts, ODD_MASKS))
    replications = [sum(mask >> vertex & 1 for mask in supports)
                    for vertex in POINTS]
    pair_counts = [
        sum((mask >> left & 1) and (mask >> right & 1)
            for mask in supports)
        for left, right in itertools.combinations(POINTS, 2)
    ]
    return replications, pair_counts


def walsh_partner_odd_counts(
        edge_even_counts: tuple[int, ...]) -> tuple[int, ...] | None:
    """Recover the unique odd-support multiset by the even Walsh transform.

    Complement every odd support, obtaining counts ``z_T`` on even subsets.
    The replication and pair equations give, for every even ``T``,

    ``z_T = (5 - sum_{|E xor T| in {0,4}} x_E) / 2``.

    A nonnegative integral partner therefore exists exactly when all sixteen
    closed-neighborhood sums are in ``{1,3,5}``.
    """
    if len(edge_even_counts) != len(EVEN_MASKS):
        raise CensusViolation("edge-count vector has the wrong length")
    if any(type(count) is not int or count < 0 for count in edge_even_counts):
        raise CensusViolation("edge multiplicities must be nonnegative integers")
    if sum(edge_even_counts) != 10:
        return None

    complement_even_counts = []
    for target in EVEN_MASKS:
        neighborhood_sum = sum(
            edge_even_counts[EVEN_INDEX[source]]
            for source in EVEN_MASKS
            if (source ^ target).bit_count() in (0, 4)
        )
        if neighborhood_sum not in (1, 3, 5):
            return None
        complement_even_counts.append((5 - neighborhood_sum) // 2)

    odd_counts = [0] * len(ODD_MASKS)
    for even_mask, count in zip(
            EVEN_MASKS, complement_even_counts, strict=True):
        odd_counts[ODD_INDEX[FULL_MASK ^ even_mask]] = count
    if sum(odd_counts) != 10:
        raise CensusViolation("the Walsh partner has the wrong block count")
    return tuple(odd_counts)


def _counts_from_bitsets(twice: int, once: int) -> tuple[int, ...]:
    return tuple(2 if twice >> index & 1 else 1 if once >> index & 1 else 0
                 for index in range(len(EVEN_MASKS)))


def _enumerate_balanced_designs() -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """Enumerate all labelled-C5 support multisets via the Walsh criterion."""
    designs = []
    for doubled_count in range(6):
        singleton_count = 10 - 2 * doubled_count
        for doubled_indices in itertools.combinations(
                range(len(EVEN_MASKS)), doubled_count):
            doubled = sum(1 << index for index in doubled_indices)
            remaining = [index for index in range(len(EVEN_MASKS))
                         if not doubled >> index & 1]
            for singleton_indices in itertools.combinations(
                    remaining, singleton_count):
                singletons = sum(1 << index for index in singleton_indices)
                neighborhood_sums = tuple(
                    (singletons & neighborhood).bit_count()
                    + 2 * (doubled & neighborhood).bit_count()
                    for neighborhood in _CLOSED_WALSH_NEIGHBORHOODS
                )
                if any(value not in (1, 3, 5)
                       for value in neighborhood_sums):
                    continue
                edge_counts = _counts_from_bitsets(doubled, singletons)
                complement_even_counts = tuple(
                    (5 - value) // 2 for value in neighborhood_sums
                )
                odd_counts = [0] * len(ODD_MASKS)
                for even_mask, count in zip(
                        EVEN_MASKS, complement_even_counts, strict=True):
                    odd_counts[ODD_INDEX[FULL_MASK ^ even_mask]] = count
                designs.append((edge_counts, tuple(odd_counts)))
    designs.sort()
    if len(designs) != len(set(designs)):
        raise CensusViolation("the Walsh enumeration produced duplicate designs")
    return designs


def _permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[vertex] for vertex in POINTS
               if mask >> vertex & 1)


def _permute_even_counts(counts: tuple[int, ...],
                         permutation: tuple[int, ...]) -> tuple[int, ...]:
    image = [0] * len(EVEN_MASKS)
    for source, count in zip(EVEN_MASKS, counts, strict=True):
        target = _permute_mask(source, permutation)
        image[EVEN_INDEX[target]] = count
    return tuple(image)


def _dihedral_orbit(edge_counts: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted({_permute_even_counts(edge_counts, permutation)
                         for permutation in DIHEDRAL_PERMUTATIONS}))


def _design_line(edge_counts: tuple[int, ...], odd_counts: tuple[int, ...],
                 orbit_size: int | None = None) -> bytes:
    edge_text = "".join(str(count) for count in edge_counts)
    odd_text = "".join(str(count) for count in odd_counts)
    suffix = "" if orbit_size is None else f"|orbit={orbit_size}"
    return f"f5-support-v1|edge={edge_text}|odd={odd_text}{suffix}\n".encode()


def _digest(lines: list[bytes]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line)
    return digest.hexdigest()


def _build_support_census() -> dict:
    designs = _enumerate_balanced_designs()
    design_map = {edge: odd for edge, odd in designs}
    if len(designs) != 7872:
        raise CensusViolation("the fixed-five support count changed")

    canonical = {}
    for edge_counts, odd_counts in designs:
        orbit = _dihedral_orbit(edge_counts)
        if any(image not in design_map for image in orbit):
            raise CensusViolation("the balanced domain is not D5-invariant")
        canonical.setdefault(orbit[0], len(orbit))
    if len(canonical) != 844:
        raise CensusViolation("the fixed-five dihedral orbit count changed")

    representatives = []
    for index, edge_counts in enumerate(sorted(canonical)):
        odd_counts = walsh_partner_odd_counts(edge_counts)
        if odd_counts is None:
            raise CensusViolation("a canonical edge multiset lost its partner")
        orbit_size = canonical[edge_counts]
        representatives.append({
            "index": index,
            "edge_even_counts": list(edge_counts),
            "nonedge_odd_counts": list(odd_counts),
            "orbit_size": orbit_size,
            "stabilizer_size": 10 // orbit_size,
        })

    histogram = {str(size): sum(value == size for value in canonical.values())
                 for size in (1, 5, 10)}
    fixed_by = [
        sum(_permute_even_counts(edge_counts, permutation) == edge_counts
            for edge_counts, _ in designs)
        for permutation in DIHEDRAL_PERMUTATIONS
    ]
    if fixed_by != [7872, 2, 2, 2, 2, 112, 112, 112, 112, 112]:
        raise CensusViolation("the D5 Burnside audit changed")
    if sum(fixed_by) // 10 != len(representatives):
        raise CensusViolation("Burnside does not reproduce the orbit count")

    control_edge = counts_from_supports(ROTATION_EDGE_SUPPORTS, EVEN_MASKS)
    control_odd = counts_from_supports(ROTATION_NONEDGE_SUPPORTS, ODD_MASKS)
    if walsh_partner_odd_counts(control_edge) != control_odd:
        raise CensusViolation("the rotation support control is not balanced")
    control_canonical = min(_dihedral_orbit(control_edge))
    representative_index = sorted(canonical).index(control_canonical)

    raw_lines = [_design_line(edge, odd) for edge, odd in designs]
    representative_lines = [
        _design_line(tuple(record["edge_even_counts"]),
                     tuple(record["nonedge_odd_counts"]),
                     record["orbit_size"])
        for record in representatives
    ]
    permutation_lines = [
        ("f5-support-d5-v1|" + ",".join(map(str, permutation)) + "\n").encode()
        for permutation in sorted(DIHEDRAL_PERMUTATIONS)
    ]
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
        "dihedral_group_order": len(DIHEDRAL_PERMUTATIONS),
        "dihedral_orbits": len(representatives),
        "orbit_size_histogram": histogram,
        "burnside_fixed_designs": {
            "identity": fixed_by[0],
            "nonidentity_rotations_each": fixed_by[1],
            "reflections_each": fixed_by[5],
        },
        "raw_designs_sha256": _digest(raw_lines),
        "representatives_sha256": _digest(representative_lines),
        "dihedral_action_sha256": _digest(permutation_lines),
        "representatives": representatives,
        "positive_control": {
            "construction": "five rotations of one support of each size 1,2,3,4",
            "edge_support_masks": list(ROTATION_EDGE_SUPPORTS),
            "nonedge_support_masks": list(ROTATION_NONEDGE_SUPPORTS),
            "canonical_representative_index": representative_index,
            "dihedral_orbit_size": len(_dihedral_orbit(control_edge)),
        },
        "disposition": "EXACT_BALANCED_SUPPORT_MULTISET_CENSUS_COMPLETE",
    }


@functools.cache
def _support_census_json() -> str:
    return json.dumps(_build_support_census(), sort_keys=True,
                      separators=(",", ":"))


def support_census() -> dict:
    """Return a fresh copy of the exact balanced-support census."""
    return json.loads(_support_census_json())


def signed_completion_model() -> dict:
    """Record the exact next-stage Seidel equations without solving them."""
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


def _build_analysis() -> dict:
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
        "support_census": support_census(),
        "signed_completion": signed_completion_model(),
        "disposition": DISPOSITION,
    }


@functools.cache
def _analysis_json() -> str:
    return json.dumps(_build_analysis(), sort_keys=True, separators=(",", ":"))


def run_analysis() -> dict:
    """Return a fresh copy of the complete support-census artifact."""
    return json.loads(_analysis_json())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "data" / "involution_f5_support_census.json",
    )
    args = parser.parse_args(argv)
    document = run_analysis()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(DISPOSITION)
    print(
        "balanced moment designs:",
        document["support_census"]["balanced_labelled_designs"],
        "labelled C5 /",
        document["support_census"]["dihedral_orbits"],
        "candidate D5 orbits; signed completion remains open",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
