#!/usr/bin/env python3
"""Replay Maksimovic's 288 published S3-invariant SRG(45) records.

The source is a GAP/GRAPE assignment, not graph6.  This validator parses the
adjacency records without executing GAP, verifies every graph has parameters
``(45,22,10,11)``, and checks both published ``S3`` generators: the order-three
permutation with nine fixed points and the involution with five fixed points.
For the involution it reconstructs the balanced support multiset and both
signed Seidel square equations.  It also runs the existing bitset clique
kernel on each graph and complement.  This is source-verified positive-control
evidence, not catalog completeness or an ``R(5,5)`` bound.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import check_ramsey

SCHEMA_VERSION = 1
CAMPAIGN_ID = "published_s3_validation"
DEFAULT_URL = "https://www.math.uniri.hr/~mmaksimovic/srg45.txt"


class ValidationViolation(RuntimeError):
    """The source payload or one of its graph records is invalid."""


def _unique_object(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise ValidationViolation(f"duplicate JSON member: {key}")
        output[key] = value
    return output


def load_document(text: str):
    """Parse strict JSON, rejecting duplicate names and nonstandard constants."""
    def reject_constant(value):
        raise ValidationViolation(f"nonstandard JSON constant: {value}")

    return json.loads(
        text,
        object_pairs_hook=_unique_object,
        parse_constant=reject_constant,
    )


def _exact_tree_equal(actual, expected) -> bool:
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return (
            actual.keys() == expected.keys()
            and all(_exact_tree_equal(actual[key], value)
                    for key, value in expected.items())
        )
    if isinstance(expected, list):
        return (
            len(actual) == len(expected)
            and all(_exact_tree_equal(left, right)
                    for left, right in zip(actual, expected))
        )
    return actual == expected

def _load_support_census(
        path: Path) -> tuple[set[tuple[tuple[int, ...], tuple[int, ...]]], dict]:
    try:
        payload = path.read_bytes()
        document = load_document(payload.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValidationViolation(f"cannot read support census: {error}") from error
    if type(document) is not dict:
        raise ValidationViolation("support census root is not an object")
    if document.get("campaign_id") != "involution_f5_support_census":
        raise ValidationViolation("support census campaign id differs")
    if document.get("disposition") != (
            "F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS"):
        raise ValidationViolation("support census disposition differs")
    census = document.get("support_census")
    if type(census) is not dict or type(census.get("representatives")) is not list:
        raise ValidationViolation("support census representatives are missing")
    if (type(census.get("dihedral_orbits")) is not int
            or census["dihedral_orbits"] != 844
            or len(census["representatives"]) != 844):
        raise ValidationViolation("support census representative count differs")
    representatives = set()
    for record in census["representatives"]:
        if type(record) is not dict:
            raise ValidationViolation("support census representative is not an object")
        edge = record.get("edge_even_counts")
        odd = record.get("nonedge_odd_counts")
        if (type(edge) is not list or type(odd) is not list
                or len(edge) != 16 or len(odd) != 16
                or any(type(value) is not int or value < 0
                       for value in edge + odd)
                or sum(edge) != 10 or sum(odd) != 10):
            raise ValidationViolation("support census representative counts differ")
        representatives.add((tuple(edge), tuple(odd)))
    if len(representatives) != 844:
        raise ValidationViolation("support census has duplicate representatives")
    return representatives, {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "disposition": document["disposition"],
        "representatives": len(representatives),
    }



def _balanced_block(text: str, start: int) -> tuple[str, int]:
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "[":
            depth += 1
        elif text[index] == "]":
            depth -= 1
            if depth == 0:
                return text[start:index + 1], index + 1
            if depth < 0:
                break
    raise ValidationViolation("unbalanced GAP adjacency list")


def _parse_rows(block: str) -> tuple[int, ...]:
    rows = []
    row_start = None
    depth = 0
    for index, character in enumerate(block):
        if character == "[":
            depth += 1
            if depth == 2:
                row_start = index + 1
        elif character == "]":
            if depth == 2:
                if row_start is None:
                    raise ValidationViolation("malformed GAP adjacency row")
                body = block[row_start:index]
                values = []
                for token in body.replace(",", " ").split():
                    try:
                        values.append(int(token))
                    except ValueError as error:
                        raise ValidationViolation("noninteger GAP adjacency entry") from error
                if any(value <= 0 for value in values):
                    raise ValidationViolation("GAP vertices must be positive")
                rows.append(sum(1 << (value - 1) for value in values))
                row_start = None
            depth -= 1
    if depth != 0:
        raise ValidationViolation("unbalanced GAP adjacency rows")
    return tuple(rows)


def parse_gap_adjacencies(payload: bytes) -> list[tuple[int, ...]]:
    """Parse every ``adjacencies := [...]`` record from a GAP text payload."""
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as error:
        raise ValidationViolation("published payload is not ASCII") from error
    marker = "adjacencies :="
    position = 0
    graphs = []
    while True:
        record = text.find(marker, position)
        if record < 0:
            break
        start = text.find("[", record + len(marker))
        if start < 0:
            raise ValidationViolation("adjacency marker has no list")
        block, position = _balanced_block(text, start)
        rows = _parse_rows(block)
        if not rows:
            raise ValidationViolation("empty GAP adjacency record")
        graphs.append(rows)
    if not graphs:
        raise ValidationViolation("payload contains no adjacency records")
    return graphs


def _specified_order3_permutation() -> tuple[int, ...]:
    permutation = list(range(45))
    for start in range(9, 45, 3):
        permutation[start] = start + 1
        permutation[start + 1] = start + 2
        permutation[start + 2] = start
    return tuple(permutation)

def _specified_involution() -> tuple[int, ...]:
    permutation = list(range(45))
    cycles = (
        (2, 3), (4, 5), (6, 7), (8, 9),
        (11, 12), (14, 15), (17, 18), (20, 21),
        (22, 25), (23, 27), (24, 26),
        (28, 31), (29, 33), (30, 32),
        (34, 37), (35, 39), (36, 38),
        (40, 43), (41, 45), (42, 44),
    )
    for left, right in cycles:
        left -= 1
        right -= 1
        permutation[left] = right
        permutation[right] = left
    return tuple(permutation)


def _fixed_cycle_order(adjacency: tuple[int, ...],
                       fixed: tuple[int, ...]) -> tuple[int, ...]:
    if len(fixed) != 5:
        raise ValidationViolation("the specified involution does not fix five points")
    start = min(fixed)
    neighbors = sorted(vertex for vertex in fixed
                       if adjacency[start] >> vertex & 1)
    if len(neighbors) != 2:
        raise ValidationViolation("the involution fixed graph is not C5")
    order = [start, neighbors[0]]
    while len(order) < 5:
        previous, current = order[-2:]
        following = [vertex for vertex in fixed
                     if adjacency[current] >> vertex & 1
                     and vertex != previous]
        if len(following) != 1 or following[0] == start:
            raise ValidationViolation("the fixed C5 traversal failed")
        order.append(following[0])
    if not adjacency[order[-1]] >> start & 1:
        raise ValidationViolation("the fixed C5 does not close")
    return tuple(order)


def _d5_permutations() -> tuple[tuple[int, ...], ...]:
    return tuple(
        [tuple((vertex + shift) % 5 for vertex in range(5))
         for shift in range(5)]
        + [tuple((shift - vertex) % 5 for vertex in range(5))
           for shift in range(5)]
    )


def _permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    return sum(1 << permutation[vertex] for vertex in range(5)
               if mask >> vertex & 1)


def _support_key(supports: tuple[int, ...],
                 internal: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    even_masks = tuple(mask for mask in range(32) if mask.bit_count() % 2 == 0)
    odd_masks = tuple(mask for mask in range(32) if mask.bit_count() % 2 == 1)
    edge = tuple(sum(a and support == mask
                     for support, a in zip(supports, internal))
                 for mask in even_masks)
    nonedge = tuple(sum(not a and support == mask
                        for support, a in zip(supports, internal))
                    for mask in odd_masks)
    return edge, nonedge


def _canonical_support_key(
        supports: tuple[int, ...],
        internal: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    images = []
    for permutation in _d5_permutations():
        images.append(_support_key(
            tuple(_permute_mask(mask, permutation) for mask in supports),
            internal,
        ))
    return min(images)


def _matrix_product_entry(left: list[list[int]], right: list[list[int]],
                          row: int, column: int) -> int:
    return sum(left[row][index] * right[index][column]
               for index in range(len(right)))


def _involution_completion_record(
        adjacency: tuple[int, ...], record_index: int,
        involution: tuple[int, ...]) -> tuple[tuple, tuple, bytes]:
    fixed = tuple(vertex for vertex in range(45)
                  if involution[vertex] == vertex)
    fixed_order = _fixed_cycle_order(adjacency, fixed)
    pairs = tuple((vertex, involution[vertex]) for vertex in range(45)
                  if vertex < involution[vertex])
    if len(pairs) != 20:
        raise ValidationViolation("the specified involution lacks twenty pairs")

    supports = []
    internal = []
    for left, right in pairs:
        if any(((adjacency[left] >> vertex) & 1)
               != ((adjacency[right] >> vertex) & 1)
               for vertex in fixed_order):
            raise ValidationViolation("an involution pair has unequal fixed support")
        supports.append(sum(((adjacency[left] >> vertex) & 1) << position
                            for position, vertex in enumerate(fixed_order)))
        internal.append((adjacency[left] >> right) & 1)
    supports = tuple(supports)
    internal = tuple(internal)
    if internal.count(1) != 10 or internal.count(0) != 10:
        raise ValidationViolation("the involution pair diagonals are not balanced")
    if any(mask.bit_count() % 2 != 1 - edge
           for mask, edge in zip(supports, internal)):
        raise ValidationViolation("an involution support has the wrong parity")
    replications = [sum(mask >> vertex & 1 for mask in supports)
                    for vertex in range(5)]
    pair_counts = [
        sum((mask >> left & 1) and (mask >> right & 1)
            for mask in supports)
        for left, right in itertools.combinations(range(5), 2)
    ]
    if replications != [10] * 5 or pair_counts != [5] * 10:
        raise ValidationViolation("the involution support design is unbalanced")

    p = [[0] * 20 for _ in range(20)]
    q = [[0] * 20 for _ in range(20)]
    for orbit, edge in enumerate(internal):
        q[orbit][orbit] = edge
    for left in range(20):
        x, y = pairs[left]
        for right in range(left + 1, 20):
            u, v = pairs[right]
            parallel = (adjacency[x] >> u) & 1
            crossed = (adjacency[x] >> v) & 1
            if ((adjacency[y] >> v) & 1) != parallel:
                raise ValidationViolation("parallel orbit edges violate the involution")
            if ((adjacency[y] >> u) & 1) != crossed:
                raise ValidationViolation("crossed orbit edges violate the involution")
            p[left][right] = p[right][left] = parallel
            q[left][right] = q[right][left] = crossed

    incidence_sign = [
        [1 - 2 * ((supports[orbit] >> vertex) & 1)
         for orbit in range(20)]
        for vertex in range(5)
    ]
    fixed_seidel = [
        [0 if left == right
         else -1 if (left - right) % 5 in (1, 4)
         else 1
         for right in range(5)]
        for left in range(5)
    ]
    invariant = [[0] * 20 for _ in range(20)]
    anti = [[0] * 20 for _ in range(20)]
    for left in range(20):
        invariant[left][left] = 1 - 2 * internal[left]
        anti[left][left] = 2 * internal[left] - 1
        for right in range(left + 1, 20):
            invariant[left][right] = invariant[right][left] = (
                2 - 2 * (p[left][right] + q[left][right])
            )
            anti[left][right] = anti[right][left] = (
                2 * (q[left][right] - p[left][right])
            )
            if invariant[left][right] ** 2 + anti[left][right] ** 2 != 4:
                raise ValidationViolation("the signed relation blocks are incompatible")
    if any(sum(value != 0 for column, value in enumerate(row)
               if column != index) != 8
           for index, row in enumerate(invariant)):
        raise ValidationViolation("an invariant Seidel row does not have eight blocks")
    if any(sum(value != 0 for column, value in enumerate(row)
               if column != index) != 11
           for index, row in enumerate(anti)):
        raise ValidationViolation("an anti-invariant Seidel row does not have eleven blocks")

    for vertex in range(5):
        for orbit in range(20):
            value = sum(fixed_seidel[vertex][other]
                        * incidence_sign[other][orbit]
                        for other in range(5))
            value += sum(incidence_sign[vertex][other]
                         * invariant[other][orbit]
                         for other in range(20))
            if value != -1:
                raise ValidationViolation("the invariant cross equation failed")
    for left in range(20):
        for right in range(20):
            invariant_square = _matrix_product_entry(
                invariant, invariant, left, right
            )
            gram = sum(incidence_sign[vertex][left]
                       * incidence_sign[vertex][right]
                       for vertex in range(5))
            expected = 45 * (left == right) - 2
            if invariant_square + 2 * gram != expected:
                raise ValidationViolation("the invariant square equation failed")
            anti_square = _matrix_product_entry(anti, anti, left, right)
            if anti_square != 45 * (left == right):
                raise ValidationViolation("the anti-invariant square equation failed")

    raw_key = _support_key(supports, internal)
    canonical_key = _canonical_support_key(supports, internal)
    relation_bits = [
        str(internal[index]) for index in range(20)
    ] + [
        f"{p[left][right]}{q[left][right]}"
        for left in range(20) for right in range(left + 1, 20)
    ]
    digest_line = (
        f"{record_index}|fixed={','.join(str(vertex + 1) for vertex in fixed)}"
        f"|edge={','.join(map(str, canonical_key[0]))}"
        f"|odd={','.join(map(str, canonical_key[1]))}"
        f"|relations={','.join(relation_bits)}\n"
    ).encode()
    return canonical_key, raw_key, digest_line


def _is_automorphism(adjacency: tuple[int, ...], permutation: tuple[int, ...]) -> bool:
    for vertex, row in enumerate(adjacency):
        mapped = 0
        remaining = row
        while remaining:
            bit = remaining & -remaining
            neighbor = bit.bit_length() - 1
            remaining ^= bit
            mapped |= 1 << permutation[neighbor]
        if mapped != adjacency[permutation[vertex]]:
            return False
    return True


def _srg_parameters_are_valid(adjacency: tuple[int, ...]) -> bool:
    if len(adjacency) != 45:
        return False
    for vertex, row in enumerate(adjacency):
        if row >> vertex & 1 or row.bit_count() != 22:
            return False
    for left in range(45):
        for right in range(left + 1, 45):
            if (adjacency[left] >> right & 1) != (adjacency[right] >> left & 1):
                return False
            expected = 10 if adjacency[left] >> right & 1 else 11
            if (adjacency[left] & adjacency[right]).bit_count() != expected:
                return False
    return True


def _complement(adjacency: tuple[int, ...]) -> tuple[int, ...]:
    full = (1 << len(adjacency)) - 1
    return tuple(full ^ (1 << vertex) ^ row
                 for vertex, row in enumerate(adjacency))


def _first_clique(adjacency: tuple[int, ...], order: int) -> tuple[int, ...] | None:
    def visit(candidates: int, need: int, chosen: tuple[int, ...]):
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


def run_validation(
        payload: bytes, source_url: str = DEFAULT_URL,
        support_representatives: set[
            tuple[tuple[int, ...], tuple[int, ...]]
        ] | None = None,
        support_dependency: dict | None = None) -> dict:
    graphs = parse_gap_adjacencies(payload)
    order_three = _specified_order3_permutation()
    involution = _specified_involution()
    coverage = hashlib.sha256()
    involution_coverage = hashlib.sha256()
    srg_valid = order3_valid = involution_valid = 0
    with_k5 = with_i5 = with_both = 0
    raw_supports = set()
    canonical_supports = set()
    for index, graph in enumerate(graphs):
        if not _srg_parameters_are_valid(graph):
            raise ValidationViolation(f"record {index} is not srg(45,22,10,11)")
        srg_valid += 1
        if not _is_automorphism(graph, order_three):
            raise ValidationViolation(
                f"record {index} lacks the specified order-three action"
            )
        order3_valid += 1
        if not _is_automorphism(graph, involution):
            raise ValidationViolation(
                f"record {index} lacks the specified involution"
            )
        involution_valid += 1
        canonical_support, raw_support, involution_line = (
            _involution_completion_record(graph, index, involution)
        )
        canonical_supports.add(canonical_support)
        raw_supports.add(raw_support)
        involution_coverage.update(involution_line)

        complement = _complement(graph)
        clique = _first_clique(graph, 5)
        independent = _first_clique(complement, 5)
        full = (1 << 45) - 1
        if bool(clique) != check_ramsey.has_clique(graph, full, 5):
            raise ValidationViolation(
                f"record {index} disagrees with the existing K5 kernel"
            )
        if bool(independent) != check_ramsey.has_clique(complement, full, 5):
            raise ValidationViolation(
                f"record {index} disagrees with the existing I5 kernel"
            )
        with_k5 += clique is not None
        with_i5 += independent is not None
        with_both += clique is not None and independent is not None
        if clique is None or independent is None:
            raise ValidationViolation(f"record {index} did not contain both K5 and I5")
        coverage.update(
            f"{index}:{','.join(map(str, clique))}:"
            f"{','.join(map(str, independent))}\n".encode()
        )
    membership_verified = support_representatives is not None
    if membership_verified:
        missing = canonical_supports - support_representatives
        if missing:
            raise ValidationViolation(
                f"{len(missing)} published support types are absent from the census"
            )
        if type(support_dependency) is not dict:
            raise ValidationViolation("support census dependency metadata is missing")
    elif support_dependency is not None:
        raise ValidationViolation("support census metadata lacks representatives")

    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "source": {
            "url": source_url,
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "format": "GAP_grape_records",
        },
        "support_census_dependency": support_dependency,
        "records": len(graphs),
        "validation": {
            "srg_45_22_10_11": srg_valid,
            "distinct_labelled_adjacency_records": len(set(graphs)),
            "specified_order3_automorphism": order3_valid,
            "order3_fixed_points": sum(
                order_three[index] == index for index in range(45)
            ),
            "order3_cycles": 12,
            "specified_involution_automorphism": involution_valid,
            "involution_fixed_points": sum(
                involution[index] == index for index in range(45)
            ),
            "involution_transpositions": 20,
            "fixed_graph_C5": involution_valid,
            "balanced_involution_support": involution_valid,
            "signed_seidel_completion_equations": involution_valid,
            "labelled_involution_support_signatures": len(raw_supports),
            "dihedral_involution_support_types": len(canonical_supports),
            "support_types_in_balanced_census": (
                len(canonical_supports) if membership_verified else 0
            ),
            "support_census_membership_verified": membership_verified,
            "involution_coverage_sha256": involution_coverage.hexdigest(),
            "with_K5": with_k5,
            "with_I5": with_i5,
            "with_both": with_both,
            "witness_coverage_sha256": coverage.hexdigest(),
        },
        "scope": {
            "role": (
                "source_verified_positive_controls_for_order3_and_"
                "fixed5_involution_signed_completion"
            ),
            "catalog_completeness_claimed": False,
            "novelty_claimed": False,
            "general_ramsey_bound_claimed": False,
        },
    }


def _read_source(source: str) -> bytes:
    parsed = urllib.parse.urlparse(source)
    if parsed.scheme in ("http", "https"):
        with urllib.request.urlopen(source, timeout=60) as response:
            return response.read()
    try:
        return Path(source).read_bytes()
    except OSError as error:
        raise ValidationViolation(f"cannot read source: {error}") from error


def main(argv=None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=DEFAULT_URL)
    parser.add_argument("--output", type=Path,
                        default=root / "data" / "published_s3_validation.json")
    parser.add_argument(
        "--support-census",
        type=Path,
        default=root / "data" / "involution_f5_support_census.json",
    )
    parser.add_argument("--verify-existing", action="store_true")
    args = parser.parse_args(argv)
    try:
        support_representatives, support_dependency = _load_support_census(
            args.support_census
        )
        try:
            relative_support_path = args.support_census.resolve().relative_to(
                root.parent.resolve()
            )
        except ValueError:
            relative_support_path = args.support_census
        support_dependency["relative_path"] = str(relative_support_path)
        document = run_validation(
            _read_source(args.source),
            args.source,
            support_representatives,
            support_dependency,
        )
        if args.verify_existing:
            expected = load_document(args.output.read_text())
            if not _exact_tree_equal(document, expected):
                raise ValidationViolation("committed artifact differs from live replay")
        else:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    except (OSError, json.JSONDecodeError, ValidationViolation) as error:
        print(f"PUBLISHED S3 VALIDATION FAILED: {error}", file=sys.stderr)
        return 1
    print(
        "PUBLISHED S3 DATA VERIFIED:",
        f"records={document['records']}",
        f"srg={document['validation']['srg_45_22_10_11']}",
        f"order3={document['validation']['specified_order3_automorphism']}",
        f"involution={document['validation']['specified_involution_automorphism']}",
        f"f5_support_types={document['validation']['dihedral_involution_support_types']}",
        f"f5_types_in_census={document['validation']['support_types_in_balanced_census']}",
        f"K5_and_I5={document['validation']['with_both']}",
        "(empirical corroboration only)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
