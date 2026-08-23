#!/usr/bin/env python3
"""Exhaust every undirected Cayley graph of order 45 for Ramsey(5,5).

Result scope
============

This computation proves only the following finite statement:

    No undirected Cayley graph on 45 vertices has both clique number and
    independence number at most four.

It does NOT prove a new bound on R(5,5): a hypothetical Ramsey(5,5,45) graph
need not be Cayley or even vertex-transitive.

Why the search is complete
==========================

There are exactly two groups of order 45.  Sylow's theorem makes the subgroup
of order five normal.  Its conjugation action has image with order dividing
both 9 and |Aut(C5)| = 4, so the action is trivial and the subgroup is central.
A Sylow 3-subgroup has order nine and is either C9 or C3 x C3.  Consequently
G is C45 or C15 x C3; both groups are abelian.

Neither group has an involution.  Its 44 nonidentity elements therefore split
into 22 inverse pairs.  An undirected simple Cayley graph Cay(G,S) chooses any
subset of these pairs.  Complementation replaces S by G \\ ({0} union S), so
fixing one selected inverse pair to be absent chooses exactly one member of
each complementary pair: 2^21 cases per group, 2^22 in total.  No degree-window
or Ramsey-number assumption enters this exhaustion.

The conference subcase has a shorter obstruction
=================================================

Either group has a character of order three.  Because S is inverse-closed,
each inverse pair contributes either 2 or omega + omega^2 = -1 to that
character's Cayley eigenvalue, so the eigenvalue is an integer.  In contrast,
a strongly regular graph with parameters (45,22,10,11) requires every
nonprincipal eigenvalue to solve x^2 + x - 11 = 0.  Both roots are irrational.
Thus no Cayley graph has those conference parameters; the exhaustive
``conference_pds = 0`` count is an independent mechanical replay of this
character-theoretic obstruction.

Why the K5 test is exact
========================

Translation is an automorphism of a Cayley graph.  Any K5 can therefore be
translated so that it contains zero.  The other four vertices lie in S and are
pairwise adjacent exactly when they form a K4 in the subgraph induced by S.
The same test on the complementary connection set detects an independent set
of order five.  Every rejected complement class is committed to an explicit,
deterministically chosen monochromatic K5 witness in ``coverage_sha256``.

The discovery module uses a recursive translated-bitset clique kernel.  The
disjoint checker reconstructs the equivalent 22-variable Ramsey CNF and checks
all 2^22 assignments at once with exact Python-integer bitsets; an optional
full resweep uses a separately written pair-intersection kernel.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

SCHEMA_VERSION = 2
CAMPAIGN_ID = "cayley_r55_order45"
DISPOSITION = "NO_CAYLEY_RAMSEY_5_5_45"
ORDER = 45
COMPLEMENT_CLASSES_PER_GROUP = 1 << 21
TOTAL_COMPLEMENT_CLASSES = 1 << 22
TOTAL_CONNECTION_SETS = 1 << 23

comb = math.comb


class SearchViolation(RuntimeError):
    """A supposedly impossible state occurred during the exact search."""


@dataclass(frozen=True)
class Group:
    """A small abelian group with elements encoded by mixed-radix integers."""

    name: str
    moduli: tuple[int, ...]
    order: int
    add: tuple[tuple[int, ...], ...]
    sub: tuple[tuple[int, ...], ...]
    neg: tuple[int, ...]
    inverse_pairs: tuple[tuple[int, int], ...]
    nonzero_mask: int


def _decode(value: int, moduli: tuple[int, ...]) -> tuple[int, ...]:
    coordinates = []
    for modulus in reversed(moduli):
        coordinates.append(value % modulus)
        value //= modulus
    if value:
        raise ValueError("mixed-radix value is out of range")
    return tuple(reversed(coordinates))


def _encode(coordinates: tuple[int, ...], moduli: tuple[int, ...]) -> int:
    value = 0
    for coordinate, modulus in zip(coordinates, moduli, strict=True):
        if not 0 <= coordinate < modulus:
            raise ValueError("mixed-radix coordinate is out of range")
        value = value * modulus + coordinate
    return value


def make_group(name: str, moduli: tuple[int, ...]) -> Group:
    """Build exact addition/subtraction tables for a product of cyclic groups."""
    if not moduli or any(type(m) is not int or m < 2 for m in moduli):
        raise ValueError("moduli must be integers >= 2")
    order = math.prod(moduli)
    coordinates = tuple(_decode(x, moduli) for x in range(order))
    neg = tuple(
        _encode(tuple((-a) % m for a, m in zip(row, moduli, strict=True)), moduli)
        for row in coordinates
    )
    add = tuple(
        tuple(
            _encode(
                tuple((a + b) % m for a, b, m in zip(left, right, moduli,
                                                       strict=True)),
                moduli,
            )
            for right in coordinates
        )
        for left in coordinates
    )
    sub = tuple(tuple(add[x][neg[y]] for y in range(order))
                for x in range(order))

    seen = {0}
    inverse_pairs = []
    for x in range(1, order):
        if x in seen:
            continue
        y = neg[x]
        if y == x:
            raise ValueError("this search requires a group with no involution")
        inverse_pairs.append((x, y))
        seen.update((x, y))
    if len(seen) != order:
        raise SearchViolation("inverse pairs do not cover the group")

    return Group(
        name=name,
        moduli=moduli,
        order=order,
        add=add,
        sub=sub,
        neg=neg,
        inverse_pairs=tuple(inverse_pairs),
        nonzero_mask=((1 << order) - 1) ^ 1,
    )


def groups_of_order_45() -> tuple[Group, Group]:
    """The two group isomorphism types of order 45."""
    return make_group("Z45", (45,)), make_group("Z15xZ3", (15, 3))


def complement_representatives(group: Group) -> Iterator[tuple[int, int]]:
    """Yield one canonically oriented pair {S, nonzero\\S} for every class.

    A Gray-code traversal toggles one inverse pair per step.  The last of the
    22 pairs stays absent, selecting exactly one of S and its complement.  The
    yielded first member has smaller degree; equal-degree pairs use the smaller
    integer bitmask.  This orientation makes all summary counts deterministic.
    """
    if group.order != ORDER or len(group.inverse_pairs) != 22:
        raise ValueError("the order-45 search requires exactly 22 inverse pairs")
    pair_masks = tuple((1 << x) | (1 << y) for x, y in group.inverse_pairs)
    connection = 0
    previous_gray = 0
    for index in range(COMPLEMENT_CLASSES_PER_GROUP):
        gray = index ^ (index >> 1)
        if index:
            changed = gray ^ previous_gray
            connection ^= pair_masks[changed.bit_length() - 1]
        previous_gray = gray
        complement = group.nonzero_mask ^ connection
        degree = connection.bit_count()
        if degree < ORDER - 1 - degree:
            yield connection, complement
        elif degree > ORDER - 1 - degree:
            yield complement, connection
        elif connection < complement:
            yield connection, complement
        else:
            yield complement, connection


def _translate_mask(group: Group, connection: int, vertex: int) -> int:
    translated = 0
    remaining = connection
    while remaining:
        bit = remaining & -remaining
        element = bit.bit_length() - 1
        remaining ^= bit
        translated |= 1 << group.add[vertex][element]
    return translated


def _first_clique(
    group: Group,
    connection: int,
    candidates: int,
    needed: int,
    prefix: tuple[int, ...],
    neighbor_cache: dict[int, int],
) -> tuple[int, ...] | None:
    """Return the lexicographically first clique extension, if one exists."""
    if needed == 0:
        return prefix
    while candidates.bit_count() >= needed:
        bit = candidates & -candidates
        vertex = bit.bit_length() - 1
        candidates ^= bit
        neighbors = neighbor_cache.get(vertex)
        if neighbors is None:
            neighbors = _translate_mask(group, connection, vertex)
            neighbor_cache[vertex] = neighbors
        found = _first_clique(
            group,
            connection,
            candidates & neighbors,
            needed - 1,
            prefix + (vertex,),
            neighbor_cache,
        )
        if found is not None:
            return found
    return None


def first_cayley_k5(group: Group, connection: int) -> tuple[int, ...] | None:
    """Return a K5 translated to contain identity 0, or ``None``."""
    if connection & ~group.nonzero_mask:
        raise ValueError("connection set contains an invalid element")
    clique4 = _first_clique(
        group, connection, connection, 4, (), {},
    )
    return None if clique4 is None else (0,) + clique4

def ramsey_cnf_patterns(group: Group) -> tuple[tuple[int, ...], ...]:
    """Return all distinct normalized K5 edge-orbit patterns.

    Translation lets every K5 contain identity zero.  For each remaining
    four-set, the ten pair differences name inverse-pair variables.  Repeated
    variables collapse because an edge orbit is selected only once.
    """
    owner = {
        element: index + 1
        for index, pair in enumerate(group.inverse_pairs)
        for element in pair
    }
    patterns = set()
    for four in itertools.combinations(range(1, group.order), 4):
        vertices = (0,) + four
        pattern = tuple(sorted({
            owner[group.sub[left][right]]
            for left, right in itertools.combinations(vertices, 2)
        }))
        patterns.add(pattern)
    return tuple(sorted(patterns))


def ramsey_cnf_clauses(
    patterns: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    """Forbid every pattern in S and in its complementary connection set."""
    clauses = []
    for pattern in patterns:
        clauses.append(pattern)
        clauses.append(tuple(-variable for variable in pattern))
    return tuple(clauses)


def _ramsey_cnf_sha256(patterns: tuple[tuple[int, ...], ...]) -> str:
    digest = hashlib.sha256()
    for clause in ramsey_cnf_clauses(patterns):
        digest.update((" ".join(str(literal) for literal in clause)
                       + " 0\n").encode("ascii"))
    return digest.hexdigest()


def is_conference_pds(group: Group, connection: int) -> bool:
    """Test the (45,22,10,11) partial-difference-set equations exactly."""
    if connection.bit_count() != 22:
        return False
    for delta in range(1, ORDER):
        common = (connection & _translate_mask(group, connection, delta)).bit_count()
        expected = 10 if connection >> delta & 1 else 11
        if common != expected:
            return False
    return True


def _witness_record(
    index: int,
    low: int,
    high: int,
    side: str,
    witness: tuple[int, ...],
) -> bytes:
    values = ",".join(str(v) for v in witness)
    return f"{index}|{low:012x}|{high:012x}|{side}|{values}\n".encode("ascii")


def search_group(group: Group) -> dict:
    """Exhaust one group and return a deterministic evidence summary."""
    degree_rows = {
        degree: {
            "minimum_degree": degree,
            "complement_classes": 0,
            "low_side_k5_free": 0,
            "ramsey_graphs": 0,
        }
        for degree in range(0, 23, 2)
    }
    coverage = hashlib.sha256()
    witnessed = 0
    ramsey_graphs = 0
    conference_pds = 0

    for index, (low, high) in enumerate(complement_representatives(group)):
        degree = low.bit_count()
        row = degree_rows[degree]
        row["complement_classes"] += 1

        low_witness = first_cayley_k5(group, low)
        if low_witness is not None:
            coverage.update(_witness_record(
                index, low, high, "low", low_witness,
            ))
            witnessed += 1
        else:
            row["low_side_k5_free"] += 1
            high_witness = first_cayley_k5(group, high)
            if high_witness is None:
                row["ramsey_graphs"] += 1
                ramsey_graphs += 1
                coverage.update(_witness_record(
                    index, low, high, "ramsey", (),
                ))
            else:
                coverage.update(_witness_record(
                    index, low, high, "high", high_witness,
                ))
                witnessed += 1

        if degree == 22 and is_conference_pds(group, low):
            conference_pds += 1

    if witnessed + ramsey_graphs != COMPLEMENT_CLASSES_PER_GROUP:
        raise SearchViolation("search did not account for every complement class")
    expected_counts = {
        degree: (comb(22, degree // 2) if degree < 22 else comb(22, 11) // 2)
        for degree in degree_rows
    }
    observed_counts = {
        degree: row["complement_classes"] for degree, row in degree_rows.items()
    }
    if observed_counts != expected_counts:
        raise SearchViolation("complement-class degree counts are incomplete")

    cnf_patterns = ramsey_cnf_patterns(group)

    return {
        "name": group.name,
        "moduli": list(group.moduli),
        "order": group.order,
        "inverse_pairs": len(group.inverse_pairs),
        "complement_classes": COMPLEMENT_CLASSES_PER_GROUP,
        "connection_sets": 2 * COMPLEMENT_CLASSES_PER_GROUP,
        "witnessed_classes": witnessed,
        "ramsey_graphs": ramsey_graphs,
        "conference_pds": conference_pds,
        "coverage_sha256": coverage.hexdigest(),
        "normalized_k5_patterns": len(cnf_patterns),
        "ramsey_cnf_clauses": 2 * len(cnf_patterns),
        "ramsey_cnf_sha256": _ramsey_cnf_sha256(cnf_patterns),
        "degree_classes": list(degree_rows.values()),
    }


def run_search() -> dict:
    """Run both complete group searches and return the canonical document."""
    groups = [search_group(group) for group in groups_of_order_45()]
    ramsey_graphs = sum(group["ramsey_graphs"] for group in groups)
    disposition = DISPOSITION if ramsey_graphs == 0 else "CAYLEY_RAMSEY_GRAPH_FOUND"
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": disposition,
        "claim": {
            "graph_family": "undirected_cayley",
            "order": ORDER,
            "forbidden_clique_order": 5,
            "forbidden_independent_set_order": 5,
            "general_ramsey_bound_claimed": False,
        },
        "total_complement_classes": TOTAL_COMPLEMENT_CLASSES,
        "total_connection_sets": TOTAL_CONNECTION_SETS,
        "ramsey_graphs": ramsey_graphs,
        "groups": groups,
    }


def _validate_document(document: dict) -> None:
    if document.get("schema_version") != SCHEMA_VERSION:
        raise SearchViolation("wrong schema version")
    if document.get("campaign_id") != CAMPAIGN_ID:
        raise SearchViolation("wrong campaign id")
    expected_claim = {
        "graph_family": "undirected_cayley",
        "order": ORDER,
        "forbidden_clique_order": 5,
        "forbidden_independent_set_order": 5,
        "general_ramsey_bound_claimed": False,
    }
    if document.get("claim") != expected_claim:
        raise SearchViolation("wrong structured claim")
    if document.get("total_complement_classes") != TOTAL_COMPLEMENT_CLASSES:
        raise SearchViolation("wrong total complement-class count")
    if document.get("total_connection_sets") != TOTAL_CONNECTION_SETS:
        raise SearchViolation("wrong connection-set count")
    if document.get("ramsey_graphs") != 0:
        raise SearchViolation("the claimed negative has a surviving graph")
    if document.get("disposition") != DISPOSITION:
        raise SearchViolation("wrong negative disposition")


def write_artifact(path: Path) -> dict:
    document = run_search()
    _validate_document(document)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    return document


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data" / "cayley_r55_45.json",
    )
    args = parser.parse_args(argv)
    document = write_artifact(args.output)
    for group in document["groups"]:
        print(
            "CAYLEY GROUP:", group["name"],
            f"classes={group['complement_classes']}",
            f"witnessed={group['witnessed_classes']}",
            f"ramsey={group['ramsey_graphs']}",
            f"conference_pds={group['conference_pds']}",
        )
    print(
        "CAYLEY EXHAUSTION:",
        f"classes={document['total_complement_classes']}",
        f"connection_sets={document['total_connection_sets']}",
    )
    print(
        "CAYLEY RESULT:", document["disposition"],
        "(no R(5,5) bound claimed)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
