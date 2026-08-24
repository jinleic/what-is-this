#!/usr/bin/env python3
"""Replay Maksimovic's 288 published S3-invariant SRG(45) records.

The source is a GAP/GRAPE assignment, not graph6.  This validator parses the
adjacency records without executing GAP, verifies every graph has parameters
(45,22,10,11), verifies the published order-three permutation with nine fixed
points and twelve 3-cycles, and runs the existing bitset clique kernel on both
the graph and its complement.  The replay is independent empirical
corroboration; it is not a catalog-completeness premise or an R(5,5) bound.
"""

from __future__ import annotations

import argparse
import hashlib
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


def run_validation(payload: bytes, source_url: str = DEFAULT_URL) -> dict:
    graphs = parse_gap_adjacencies(payload)
    permutation = _specified_order3_permutation()
    coverage = hashlib.sha256()
    srg_valid = automorphism_valid = with_k5 = with_i5 = with_both = 0
    for index, graph in enumerate(graphs):
        if not _srg_parameters_are_valid(graph):
            raise ValidationViolation(f"record {index} is not srg(45,22,10,11)")
        srg_valid += 1
        if not _is_automorphism(graph, permutation):
            raise ValidationViolation(f"record {index} lacks the specified order-three action")
        automorphism_valid += 1
        complement = _complement(graph)
        clique = _first_clique(graph, 5)
        independent = _first_clique(complement, 5)
        full = (1 << 45) - 1
        if bool(clique) != check_ramsey.has_clique(graph, full, 5):
            raise ValidationViolation(f"record {index} disagrees with the existing K5 kernel")
        if bool(independent) != check_ramsey.has_clique(complement, full, 5):
            raise ValidationViolation(f"record {index} disagrees with the existing I5 kernel")
        with_k5 += clique is not None
        with_i5 += independent is not None
        with_both += clique is not None and independent is not None
        if clique is None or independent is None:
            raise ValidationViolation(f"record {index} did not contain both K5 and I5")
        coverage.update(
            f"{index}:{','.join(map(str, clique))}:"
            f"{','.join(map(str, independent))}\n".encode()
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "source": {
            "url": source_url,
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "format": "GAP_grape_records",
        },
        "records": len(graphs),
        "validation": {
            "srg_45_22_10_11": srg_valid,
            "distinct_labelled_adjacency_records": len(set(graphs)),
            "specified_order3_automorphism": automorphism_valid,
            "order3_fixed_points": sum(permutation[index] == index for index in range(45)),
            "order3_cycles": 12,
            "with_K5": with_k5,
            "with_I5": with_i5,
            "with_both": with_both,
            "witness_coverage_sha256": coverage.hexdigest(),
        },
        "scope": {
            "role": "independent_empirical_corroboration_of_order3_fixed9_exclusion",
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
    parser.add_argument("--verify-existing", action="store_true")
    args = parser.parse_args(argv)
    try:
        document = run_validation(_read_source(args.source), args.source)
        if args.verify_existing:
            expected = json.loads(args.output.read_text())
            if document != expected:
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
        f"K5_and_I5={document['validation']['with_both']}",
        "(empirical corroboration only)",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
