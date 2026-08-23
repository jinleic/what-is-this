#!/usr/bin/env python3
"""Disjoint exact checker for the order-45 Cayley Ramsey exhaustion.

This file imports nothing from ``cayley_r55_search``.  Its default proof path
rebuilds both groups, normalizes every possible translated K5 to a clause over
the 22 inverse-pair variables, and checks the resulting CNF on all 2^22 truth
assignments using exact Python-integer bitsets.  This takes seconds and shares
neither the producer's candidate enumeration nor its clique kernel.

With ``--resweep`` it additionally rebuilds all 2^21 complement classes per
group, every degree-class count, every monochromatic K5 witness, both coverage
digests, and the (45,22,10,11) partial-difference-set count.  That defense-in-
depth path constructs connection-set adjacency from subtraction membership and
finds K4 by edge-common-neighbour intersections, rather than producer recursion.

The accepted conclusion is only ``NO_CAYLEY_RAMSEY_5_5_45``.  It is not a new
bound on R(5,5).
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

SCHEMA_VERSION = 2
CAMPAIGN_ID = "cayley_r55_order45"
DISPOSITION = "NO_CAYLEY_RAMSEY_5_5_45"
ORDER = 45
CLASSES_PER_GROUP = 1 << 21
TOTAL_CLASSES = 1 << 22
TOTAL_CONNECTION_SETS = 1 << 23


class CheckViolation(RuntimeError):
    """The artifact disagrees with the independent exact reconstruction."""


@dataclass(frozen=True)
class AuditGroup:
    name: str
    moduli: tuple[int, ...]
    minus: tuple[tuple[int, ...], ...]
    inverse_pairs: tuple[tuple[int, int], ...]
    nonzero_mask: int


def _coordinates(value: int, moduli: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * len(moduli)
    for index in range(len(moduli) - 1, -1, -1):
        result[index] = value % moduli[index]
        value //= moduli[index]
    if value:
        raise CheckViolation("group element is outside its mixed-radix range")
    return tuple(result)


def _number(coordinates: tuple[int, ...], moduli: tuple[int, ...]) -> int:
    result = 0
    for value, modulus in zip(coordinates, moduli, strict=True):
        result = modulus * result + value
    return result


def _audit_group(name: str, moduli: tuple[int, ...]) -> AuditGroup:
    order = math.prod(moduli)
    rows = tuple(_coordinates(value, moduli) for value in range(order))
    minus = tuple(
        tuple(
            _number(
                tuple((a - b) % modulus
                      for a, b, modulus in zip(left, right, moduli, strict=True)),
                moduli,
            )
            for right in rows
        )
        for left in rows
    )
    inverse = tuple(minus[0][value] for value in range(order))
    unseen = set(range(1, order))
    pairs = []
    while unseen:
        first = min(unseen)
        second = inverse[first]
        if second == first or second not in unseen:
            raise CheckViolation("order-45 inverse pairing is malformed")
        unseen.remove(first)
        unseen.remove(second)
        pairs.append((first, second))
    return AuditGroup(
        name=name,
        moduli=moduli,
        minus=minus,
        inverse_pairs=tuple(pairs),
        nonzero_mask=((1 << order) - 1) - 1,
    )


def _groups() -> tuple[AuditGroup, AuditGroup]:
    return _audit_group("Z45", (45,)), _audit_group("Z15xZ3", (15, 3))


def _representatives(group: AuditGroup):
    """Independent reconstruction of the producer's canonical Gray traversal."""
    pair_masks = tuple((1 << left) | (1 << right)
                       for left, right in group.inverse_pairs)
    connection = 0
    old_code = 0
    for number in range(CLASSES_PER_GROUP):
        code = number ^ (number // 2)
        if number:
            changed = code ^ old_code
            connection ^= pair_masks[changed.bit_length() - 1]
        old_code = code
        other = group.nonzero_mask ^ connection
        left_degree = connection.bit_count()
        if left_degree < 44 - left_degree:
            yield connection, other
        elif left_degree > 44 - left_degree:
            yield other, connection
        else:
            yield min(connection, other), max(connection, other)


def _members(mask: int) -> tuple[int, ...]:
    result = []
    while mask:
        bit = mask & -mask
        result.append(bit.bit_length() - 1)
        mask ^= bit
    return tuple(result)


def _first_k5(group: AuditGroup, connection: int) -> tuple[int, ...] | None:
    """Pair-intersection K4 search; independent of producer recursion."""
    vertices = _members(connection)
    neighbor_cache: dict[int, int] = {}

    def neighbors(vertex: int) -> int:
        answer = neighbor_cache.get(vertex)
        if answer is None:
            answer = 0
            for other in vertices:
                difference = group.minus[other][vertex]
                if connection >> difference & 1:
                    answer |= 1 << other
            neighbor_cache[vertex] = answer
        return answer

    for first_index, first in enumerate(vertices):
        first_neighbors = neighbors(first)
        for second in vertices[first_index + 1:]:
            if not (first_neighbors >> second & 1):
                continue
            common = first_neighbors & neighbors(second)
            common &= ~((1 << (second + 1)) - 1)
            while common:
                third_bit = common & -common
                third = third_bit.bit_length() - 1
                common ^= third_bit
                fourths = common & neighbors(third)
                if fourths:
                    fourth = (fourths & -fourths).bit_length() - 1
                    return 0, first, second, third, fourth
    return None

def _cnf_patterns(group: AuditGroup) -> tuple[tuple[int, ...], ...]:
    variable = [0] * ORDER
    for index, pair in enumerate(group.inverse_pairs, start=1):
        for element in pair:
            variable[element] = index
    patterns = set()
    for four in itertools.combinations(range(1, ORDER), 4):
        vertices = (0,) + four
        patterns.add(tuple(sorted({
            variable[group.minus[left][right]]
            for left, right in itertools.combinations(vertices, 2)
        })))
    return tuple(sorted(patterns))


def _cnf_clauses(
    patterns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    clauses = []
    for pattern in patterns:
        clauses.extend((pattern, tuple(-value for value in pattern)))
    return tuple(clauses)


def _cnf_sha256(clauses: tuple[tuple[int, ...], ...]) -> str:
    digest = hashlib.sha256()
    for clause in clauses:
        digest.update((" ".join(str(value) for value in clause)
                       + " 0\n").encode("ascii"))
    return digest.hexdigest()


def _assignment_masks(variable_count: int) -> tuple[int, tuple[int, ...]]:
    assignment_count = 1 << variable_count
    all_assignments = (1 << assignment_count) - 1
    masks = []
    for index in range(variable_count):
        block = 1 << index
        period = 2 * block
        repeated = all_assignments // ((1 << period) - 1)
        masks.append((((1 << block) - 1) << block) * repeated)
    return all_assignments, tuple(masks)


def _truth_table_survivors(
    clauses: tuple[tuple[int, ...], ...],
    variable_count: int,
) -> tuple[int, int]:
    """Return surviving assignment count and clauses consumed."""
    all_assignments, true_masks = _assignment_masks(variable_count)
    survivors = all_assignments
    used = 0
    for clause in sorted(clauses, key=lambda row: (len(row), row)):
        falsifies = all_assignments
        for literal in clause:
            truth = true_masks[abs(literal) - 1]
            falsifies &= ((all_assignments ^ truth) if literal > 0 else truth)
        survivors &= all_assignments ^ falsifies
        used += 1
        if not survivors:
            break
    return survivors.bit_count(), used


def _cnf_evidence(group: AuditGroup) -> dict:
    patterns = _cnf_patterns(group)
    clauses = _cnf_clauses(patterns)
    survivors, used = _truth_table_survivors(clauses, len(group.inverse_pairs))
    return {
        "normalized_k5_patterns": len(patterns),
        "ramsey_cnf_clauses": len(clauses),
        "ramsey_cnf_sha256": _cnf_sha256(clauses),
        "truth_table_survivors": survivors,
        "truth_table_clauses_used": used,
    }


def _is_conference_pds(group: AuditGroup, connection: int) -> bool:
    if connection.bit_count() != 22:
        return False
    vertices = _members(connection)
    for delta in range(1, ORDER):
        common = sum(
            1 for vertex in vertices
            if connection >> group.minus[vertex][delta] & 1
        )
        if common != (10 if connection >> delta & 1 else 11):
            return False
    return True


def _record(index: int, low: int, high: int, side: str,
            witness: tuple[int, ...]) -> bytes:
    values = ",".join(str(value) for value in witness)
    return f"{index}|{low:012x}|{high:012x}|{side}|{values}\n".encode("ascii")


def _resweep(group: AuditGroup) -> dict:
    rows = {
        degree: {
            "minimum_degree": degree,
            "complement_classes": 0,
            "low_side_k5_free": 0,
            "ramsey_graphs": 0,
        }
        for degree in range(0, 23, 2)
    }
    digest = hashlib.sha256()
    witnessed = 0
    ramsey = 0
    pds = 0

    for index, (low, high) in enumerate(_representatives(group)):
        row = rows[low.bit_count()]
        row["complement_classes"] += 1
        witness = _first_k5(group, low)
        if witness is not None:
            digest.update(_record(index, low, high, "low", witness))
            witnessed += 1
        else:
            row["low_side_k5_free"] += 1
            witness = _first_k5(group, high)
            if witness is None:
                digest.update(_record(index, low, high, "ramsey", ()))
                row["ramsey_graphs"] += 1
                ramsey += 1
            else:
                digest.update(_record(index, low, high, "high", witness))
                witnessed += 1
        if low.bit_count() == 22 and _is_conference_pds(group, low):
            pds += 1

    patterns = _cnf_patterns(group)
    clauses = _cnf_clauses(patterns)

    return {
        "name": group.name,
        "moduli": list(group.moduli),
        "order": ORDER,
        "inverse_pairs": 22,
        "complement_classes": CLASSES_PER_GROUP,
        "connection_sets": 2 * CLASSES_PER_GROUP,
        "witnessed_classes": witnessed,
        "ramsey_graphs": ramsey,
        "conference_pds": pds,
        "coverage_sha256": digest.hexdigest(),
        "normalized_k5_patterns": len(patterns),
        "ramsey_cnf_clauses": len(clauses),
        "ramsey_cnf_sha256": _cnf_sha256(clauses),
        "degree_classes": list(rows.values()),
    }


def _require_exact_keys(value, expected: set[str], label: str) -> None:
    if not isinstance(value, dict):
        raise CheckViolation(f"{label}: expected object")
    if set(value) != expected:
        raise CheckViolation(
            f"{label}: key set differs; got {sorted(value)}, expected {sorted(expected)}"
        )


def _validate_static(document: dict) -> None:
    _require_exact_keys(document, {
        "schema_version", "campaign_id", "disposition", "claim",
        "total_complement_classes", "total_connection_sets", "ramsey_graphs",
        "groups",
    }, "document")
    if document["schema_version"] != SCHEMA_VERSION:
        raise CheckViolation("schema version differs")
    if document["campaign_id"] != CAMPAIGN_ID:
        raise CheckViolation("campaign id differs")
    if document["disposition"] != DISPOSITION:
        raise CheckViolation("disposition differs")
    if document["total_complement_classes"] != TOTAL_CLASSES:
        raise CheckViolation("total complement-class count differs")
    if document["total_connection_sets"] != TOTAL_CONNECTION_SETS:
        raise CheckViolation("connection-set count differs")
    if document["ramsey_graphs"] != 0:
        raise CheckViolation("negative disposition records a surviving graph")
    expected_claim = {
        "graph_family": "undirected_cayley",
        "order": ORDER,
        "forbidden_clique_order": 5,
        "forbidden_independent_set_order": 5,
        "general_ramsey_bound_claimed": False,
    }
    if document["claim"] != expected_claim:
        raise CheckViolation("structured claim differs")
    groups = document["groups"]
    if not isinstance(groups, list) or len(groups) != 2:
        raise CheckViolation("expected exactly two group results")

    expected_groups = (("Z45", [45]), ("Z15xZ3", [15, 3]))
    expected_row_keys = {
        "minimum_degree", "complement_classes", "low_side_k5_free",
        "ramsey_graphs",
    }
    for record, (name, moduli) in zip(groups, expected_groups, strict=True):
        _require_exact_keys(record, {
            "name", "moduli", "order", "inverse_pairs", "complement_classes",
            "connection_sets", "witnessed_classes", "ramsey_graphs",
            "conference_pds", "coverage_sha256", "normalized_k5_patterns",
            "ramsey_cnf_clauses", "ramsey_cnf_sha256", "degree_classes",
        }, f"group {name}")
        if record["name"] != name or record["moduli"] != moduli:
            raise CheckViolation(f"group {name}: identity differs")
        if record["order"] != ORDER or record["inverse_pairs"] != 22:
            raise CheckViolation(f"group {name}: order or inverse-pair count differs")
        if record["complement_classes"] != CLASSES_PER_GROUP:
            raise CheckViolation(f"group {name}: complement-class count differs")
        if record["connection_sets"] != 2 * CLASSES_PER_GROUP:
            raise CheckViolation(f"group {name}: connection-set count differs")
        if record["witnessed_classes"] != CLASSES_PER_GROUP:
            raise CheckViolation(f"group {name}: not every class is witnessed")
        if record["ramsey_graphs"] != 0 or record["conference_pds"] != 0:
            raise CheckViolation(f"group {name}: negative count differs")
        digest = record["coverage_sha256"]
        if not isinstance(digest, str) or len(digest) != 64:
            raise CheckViolation(f"group {name}: malformed coverage digest")
        if type(record["normalized_k5_patterns"]) is not int:
            raise CheckViolation(f"group {name}: malformed K5-pattern count")
        if type(record["ramsey_cnf_clauses"]) is not int:
            raise CheckViolation(f"group {name}: malformed CNF-clause count")
        cnf_digest = record["ramsey_cnf_sha256"]
        if not isinstance(cnf_digest, str) or len(cnf_digest) != 64:
            raise CheckViolation(f"group {name}: malformed CNF digest")
        rows = record["degree_classes"]
        if not isinstance(rows, list) or len(rows) != 12:
            raise CheckViolation(f"group {name}: expected 12 degree classes")
        for index, row in enumerate(rows):
            _require_exact_keys(row, expected_row_keys,
                                f"group {name} degree row {index}")
            degree = 2 * index
            expected = (math.comb(22, index) if index < 11
                        else math.comb(22, 11) // 2)
            if row["minimum_degree"] != degree:
                raise CheckViolation(f"group {name}: degree order differs")
            if row["complement_classes"] != expected:
                raise CheckViolation(f"group {name}: binomial class count differs")
            if not 0 <= row["low_side_k5_free"] <= expected:
                raise CheckViolation(f"group {name}: K5-free count is out of range")
            if row["ramsey_graphs"] != 0:
                raise CheckViolation(f"group {name}: degree row has a survivor")


@lru_cache(maxsize=1)
def _independent_results() -> tuple[dict, dict]:
    """Cache one process-local resweep; CLI invocations still start from zero."""
    return tuple(_resweep(group) for group in _groups())


def verify_document(document: dict, *, resweep: bool = False) -> dict:
    """Prove CNF unsatisfiability and optionally repeat every witness search."""
    _validate_static(document)
    proof_rows = []
    for record, group in zip(document["groups"], _groups(), strict=True):
        evidence = _cnf_evidence(group)
        for field in (
            "normalized_k5_patterns", "ramsey_cnf_clauses", "ramsey_cnf_sha256",
        ):
            if record[field] != evidence[field]:
                raise CheckViolation(
                    f"group {record['name']}: {field} differs from reconstruction"
                )
        if evidence["truth_table_survivors"] != 0:
            raise CheckViolation(
                f"group {record['name']}: Ramsey CNF has surviving assignments"
            )
        proof_rows.append(evidence)
    if resweep:
        observed = list(_independent_results())
        if observed != document["groups"]:
            for expected, actual in zip(document["groups"], observed, strict=True):
                if expected != actual:
                    raise CheckViolation(
                        f"group {expected['name']}: independent resweep differs"
                    )
            raise CheckViolation("independent resweep differs")
    return {
        "complement_classes": TOTAL_CLASSES,
        "connection_sets": TOTAL_CONNECTION_SETS,
        "ramsey_graphs": 0,
        "truth_table_survivors": sum(
            row["truth_table_survivors"] for row in proof_rows
        ),
        "truth_table_clauses_used": [
            row["truth_table_clauses_used"] for row in proof_rows
        ],
        "witness_replay_verified": resweep,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument(
        "--resweep",
        action="store_true",
        help="also repeat all 4,194,304 witness searches and coverage digests",
    )
    args = parser.parse_args(argv)
    try:
        document = json.loads(args.artifact.read_text())
        summary = verify_document(document, resweep=args.resweep)
    except (OSError, json.JSONDecodeError, CheckViolation) as error:
        print(f"CAYLEY CHECK FAILED: {error}", file=sys.stderr)
        return 1
    label = "CAYLEY EVIDENCE VERIFIED:" if args.resweep else "CAYLEY THEOREM VERIFIED:"
    print(
        label, DISPOSITION,
        f"classes={summary['complement_classes']}",
        f"connection_sets={summary['connection_sets']}",
        f"truth_survivors={summary['truth_table_survivors']}",
        f"witness_replay={'VERIFIED' if args.resweep else 'NOT_RUN'}",
        "(no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
