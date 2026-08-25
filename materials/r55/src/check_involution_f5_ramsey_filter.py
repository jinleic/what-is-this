#!/usr/bin/env python3
"""Independent checker for the fixed-five exact R(3,3) support filter.

This checker does not import the producer.  It hash-binds the prior exact
balanced-support artifact, exhaustively verifies R(3,3)=6 on labelled graphs,
re-derives the four-vertex relation domains, and reconstructs all 844 filter
records before exact recursive comparison.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

SCHEMA_VERSION = 1
CAMPAIGN_ID = "involution_f5_ramsey_filter"
DISPOSITION = "F5_RAMSEY_R33_SUPPORT_FILTER_EXACT_6627_LABELLED_705_D5_ORBITS"
INPUT_DISPOSITION = "F5_BALANCED_SUPPORT_CENSUS_EXACT_7872_LABELLED_844_D5_ORBITS"
SIGNED_STATUS = "OPEN_SIGNED_COMPLETION_FOR_705_R33_SURVIVING_SUPPORT_ORBITS"
INPUT_BYTES = 503_066
INPUT_SHA256 = "520e2cb453cc2efee7925636af018316acca3e375a8107886c5e9c90c3727aaa"
INPUT_REPRESENTATIVES_SHA256 = (
    "8d53fbe601442c9a481ce729cd7d0baaf4f3727caebcc94c2a14afeb1ae8b2f1"
)
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "data" / "involution_f5_support_census.json"
DEFAULT_ARTIFACT = ROOT / "data" / "involution_f5_ramsey_filter.json"
WINDOWS = (
    ((0, 1), (2, 4)),
    ((1, 2), (0, 3)),
    ((2, 3), (1, 4)),
    ((3, 4), (0, 2)),
    ((0, 4), (1, 3)),
)


class CheckViolation(RuntimeError):
    """The artifact, dependency, or finite proof failed verification."""


def _reject_constant(token: str):
    raise CheckViolation(f"nonstandard JSON constant {token}")


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise CheckViolation(f"duplicate JSON member {key!r}")
        result[key] = value
    return result


def load_document(payload: str | bytes) -> dict:
    """Parse strict JSON, rejecting duplicates and nonstandard constants."""
    try:
        if isinstance(payload, bytes):
            payload = payload.decode("utf-8")
        if type(payload) is not str:
            raise CheckViolation("JSON payload must be str or bytes")
        document = json.loads(
            payload,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise CheckViolation(f"invalid JSON: {exc}") from exc
    if type(document) is not dict:
        raise CheckViolation("JSON root must be an object")
    return document


def _read_strict(path: Path) -> tuple[dict, bytes]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise CheckViolation(f"cannot read {path}: {exc}") from exc
    return load_document(raw), raw


def _require(value, expected_type: type, label: str) -> None:
    if type(value) is not expected_type:
        raise CheckViolation(
            f"{label} must have exact type {expected_type.__name__}")


def _homogeneous_triangle(edge_bits: int, order: int, independent: bool) -> bool:
    bit = 0
    adjacency = [[False] * order for _ in range(order)]
    for left in range(order):
        for right in range(left + 1, order):
            present = bool(edge_bits >> bit & 1)
            adjacency[left][right] = adjacency[right][left] = present
            bit += 1
    for left in range(order):
        for middle in range(left + 1, order):
            for right in range(middle + 1, order):
                values = (
                    adjacency[left][middle],
                    adjacency[left][right],
                    adjacency[middle][right],
                )
                if independent and not any(values):
                    return True
                if not independent and all(values):
                    return True
    return False


@functools.cache
def verify_r33_six_vertex_lemma() -> bool:
    """Exhaust all 2^15 labelled six-vertex graphs and verify R(3,3)=6."""
    for edge_bits in range(1 << 15):
        if not (_homogeneous_triangle(edge_bits, 6, independent=False)
                or _homogeneous_triangle(edge_bits, 6, independent=True)):
            raise CheckViolation("found a six-vertex R(3,3) counterexample")
    # C5 verifies that the bound is genuinely six rather than five.
    cycle_bits = 0
    bit = 0
    cycle_edges = {tuple(sorted((vertex, (vertex + 1) % 5)))
                   for vertex in range(5)}
    for left in range(5):
        for right in range(left + 1, 5):
            if (left, right) in cycle_edges:
                cycle_bits |= 1 << bit
            bit += 1
    if (_homogeneous_triangle(cycle_bits, 5, independent=False)
            or _homogeneous_triangle(cycle_bits, 5, independent=True)):
        raise CheckViolation("C5 failed the R(3,3) sharpness control")
    return True


def _four_vertex_graph(left_internal: int, right_internal: int,
                       parallel: int, crossed: int) -> int:
    edges = set()
    if left_internal:
        edges.add((0, 1))
    if right_internal:
        edges.add((2, 3))
    if parallel:
        edges.update(((0, 2), (1, 3)))
    if crossed:
        edges.update(((0, 3), (1, 2)))
    bits = 0
    bit = 0
    for left in range(4):
        for right in range(left + 1, 4):
            if (left, right) in edges:
                bits |= 1 << bit
            bit += 1
    return bits


def enumerate_four_vertex_allowed_w_values(
        left_internal: int, right_internal: int) -> tuple[int, ...]:
    """Independently enumerate p/q states with neither K3 nor I3."""
    if (type(left_internal) is not int or left_internal not in (0, 1)
            or type(right_internal) is not int
            or right_internal not in (0, 1)):
        raise CheckViolation("internal bits must be exact integers zero or one")
    allowed = set()
    for parallel in (0, 1):
        for crossed in (0, 1):
            bits = _four_vertex_graph(
                left_internal, right_internal, parallel, crossed)
            if (_homogeneous_triangle(bits, 4, independent=False)
                    or _homogeneous_triangle(bits, 4, independent=True)):
                continue
            allowed.add(1 - parallel - crossed)
    return tuple(sorted(allowed))


def _count_vector(value, label: str) -> tuple[int, ...]:
    _require(value, list, label)
    if len(value) != 16:
        raise CheckViolation(f"{label} must have length sixteen")
    result = []
    for index, count in enumerate(value):
        _require(count, int, f"{label}[{index}]")
        if count < 0:
            raise CheckViolation(f"{label}[{index}] is negative")
        result.append(count)
    if sum(result) != 10:
        raise CheckViolation(f"{label} must sum to ten")
    return tuple(result)


def _expand(census: dict, representative: dict
            ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    _require(census, dict, "support_census")
    _require(representative, dict, "representative")
    even_order = census.get("edge_mask_order")
    odd_order = census.get("nonedge_mask_order")
    _require(even_order, list, "edge_mask_order")
    _require(odd_order, list, "nonedge_mask_order")
    if len(even_order) != 16 or len(odd_order) != 16:
        raise CheckViolation("mask orders must have length sixteen")
    for index, mask in enumerate(even_order):
        _require(mask, int, f"edge_mask_order[{index}]")
        if not 0 <= mask < 32 or mask.bit_count() % 2:
            raise CheckViolation("invalid even mask order")
    for index, mask in enumerate(odd_order):
        _require(mask, int, f"nonedge_mask_order[{index}]")
        if not 0 <= mask < 32 or mask.bit_count() % 2 != 1:
            raise CheckViolation("invalid odd mask order")
    even_counts = _count_vector(
        representative.get("edge_even_counts"), "edge_even_counts")
    odd_counts = _count_vector(
        representative.get("nonedge_odd_counts"), "nonedge_odd_counts")
    masks = []
    internal = []
    for mask, count in zip(even_order, even_counts, strict=True):
        masks.extend([mask] * count)
        internal.extend([1] * count)
    for mask, count in zip(odd_order, odd_counts, strict=True):
        masks.extend([mask] * count)
        internal.extend([0] * count)
    if len(masks) != 20:
        raise CheckViolation("support did not expand to twenty orbits")
    return tuple(masks), tuple(internal)


def _selected(masks, edge, nonedge):
    return tuple(index for index, mask in enumerate(masks)
                 if all(mask >> vertex & 1 for vertex in edge)
                 and all(not (mask >> vertex & 1) for vertex in nonedge))


def _multiplicities(masks):
    return tuple(len(_selected(masks, edge, nonedge))
                 for edge, nonedge in WINDOWS)


def _restrictions(masks, internal):
    records = []
    seen = set()
    for window_index, (edge, nonedge) in enumerate(WINDOWS):
        shared = _selected(masks, edge, nonedge)
        if len(shared) != 2:
            continue
        pair = tuple(shared)
        if pair in seen:
            raise CheckViolation("duplicate relation restriction")
        seen.add(pair)
        left, right = pair
        records.append({
            "left_orbit": left,
            "right_orbit": right,
            "source_window_index": window_index,
            "w_ij_over_2": list(enumerate_four_vertex_allowed_w_values(
                internal[left], internal[right])),
        })
    return records


def _histogram(values):
    counts = Counter(values)
    return {str(key): counts[key] for key in (1, 5, 10)}


def _digest(records):
    digest = hashlib.sha256()
    for record in records:
        digest.update(json.dumps(
            record, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _signed_frontier():
    degrees = {}
    for support_size in range(6):
        plus = 3 + support_size // 2
        degrees[str(support_size)] = {
            "W=+2": plus,
            "W=-2": 8 - plus,
        }
    return {
        "deterministic_orbit_order": (
            "ten even supports in edge_mask_order, then ten odd supports in "
            "nonedge_mask_order"
        ),
        "w_action_on_one": "W 1=-R^T 1",
        "w_action_on_support_sum": "W R^T 1=-5 1",
        "w_sign_degrees_by_support_size": degrees,
        "global_unordered_relation_counts": {
            "W=+2": 40,
            "W=-2": 40,
            "T-supported": 110,
        },
        "w_characteristic_polynomial": "(t^2-5)^3(t^2-45)^7",
        "t_characteristic_polynomial": "(t^2-45)^10",
        "four_vertex_domains": {
            "both_internal_nonedges": "W_ij/2 in {-1,0}",
            "mixed_internal_type": "W_ij=0",
            "both_internal_edges": "W_ij/2 in {0,+1}",
        },
        "status": SIGNED_STATUS,
    }


def _build_expected(input_path: Path) -> dict:
    if input_path.resolve() != DEFAULT_INPUT.resolve():
        raise CheckViolation(
            "only the canonical balanced-support artifact is accepted")
    support, raw = _read_strict(input_path)
    digest = hashlib.sha256(raw).hexdigest()
    if len(raw) != INPUT_BYTES or digest != INPUT_SHA256:
        raise CheckViolation("balanced-support dependency commitment changed")
    _require(support.get("schema_version"), int, "support schema_version")
    if support["schema_version"] != 1:
        raise CheckViolation("unexpected support schema version")
    if support.get("campaign_id") != "involution_f5_support_census":
        raise CheckViolation("unexpected support campaign")
    if support.get("disposition") != INPUT_DISPOSITION:
        raise CheckViolation("unexpected support disposition")
    census = support.get("support_census")
    _require(census, dict, "support_census")
    if census.get("balanced_labelled_designs") != 7872:
        raise CheckViolation("unexpected labelled support count")
    if census.get("dihedral_orbits") != 844:
        raise CheckViolation("unexpected support orbit count")
    if census.get("orbit_size_histogram") != {"1": 2, "5": 110, "10": 732}:
        raise CheckViolation("unexpected support orbit-size histogram")
    if census.get("representatives_sha256") != INPUT_REPRESENTATIVES_SHA256:
        raise CheckViolation("unexpected support representative commitment")
    representatives = census.get("representatives")
    _require(representatives, list, "representatives")
    if len(representatives) != 844:
        raise CheckViolation("expected 844 representatives")
    verify_r33_six_vertex_lemma()

    candidates = []
    rejected_sizes = []
    survivor_sizes = []
    occurrence = Counter()
    labelled_occurrence = Counter()
    domain_histogram = Counter()
    total_restrictions = 0
    for index, representative in enumerate(representatives):
        _require(representative.get("index"), int,
                 f"representative[{index}].index")
        if representative["index"] != index:
            raise CheckViolation("noncontiguous representative index")
        _require(representative.get("orbit_size"), int,
                 f"representative[{index}].orbit_size")
        orbit_size = representative["orbit_size"]
        if orbit_size not in (1, 5, 10):
            raise CheckViolation("invalid representative orbit size")
        masks, internal = _expand(census, representative)
        multiplicities = _multiplicities(masks)
        for value in multiplicities:
            occurrence[value] += 1
            labelled_occurrence[value] += orbit_size
        violating = [window for window, value in enumerate(multiplicities)
                     if value >= 3]
        if violating:
            status = "REJECTED_FORCED_R33_SIX_SET"
            restrictions = []
            rejected_sizes.append(orbit_size)
        else:
            status = "SURVIVES_R33_SUPPORT_FILTER"
            restrictions = _restrictions(masks, internal)
            survivor_sizes.append(orbit_size)
            total_restrictions += len(restrictions)
            for restriction in restrictions:
                domain_histogram[",".join(map(
                    str, restriction["w_ij_over_2"]))] += 1
        candidates.append({
            "source_index": index,
            "orbit_size": orbit_size,
            "opposite_window_multiplicities": list(multiplicities),
            "status": status,
            "violating_window_indices": violating,
            "relation_restrictions": restrictions,
        })

    if (len(rejected_sizes), sum(rejected_sizes),
            len(survivor_sizes), sum(survivor_sizes)) != (139, 1245, 705, 6627):
        raise CheckViolation("reconstructed frontier counts changed")
    if total_restrictions != 1227:
        raise CheckViolation("reconstructed restriction count changed")

    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_support_census.json",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "disposition": support["disposition"],
            "representatives_sha256": census["representatives_sha256"],
        },
        "theorem": {
            "fixed_graph": "C5",
            "ramsey_fact": "R(3,3)=6",
            "opposite_windows": [
                {"fixed_edge": list(edge), "fixed_nonedge": list(nonedge)}
                for edge, nonedge in WINDOWS
            ],
            "maximum_shared_transposition_orbits": 2,
            "proof": [
                "Each selected orbit contributes both of its vertices.",
                "A triangle among six selected vertices joins the fixed edge to form K5.",
                "An independent triple joins the fixed nonedge to form I5.",
                "R(3,3)=6 therefore forbids three selected transposition orbits.",
            ],
        },
        "census": {
            "input_d5_orbits": 844,
            "input_labelled_supports": 7872,
            "input_orbit_size_histogram": census["orbit_size_histogram"],
            "rejected_d5_orbits": len(rejected_sizes),
            "rejected_labelled_supports": sum(rejected_sizes),
            "rejected_orbit_size_histogram": _histogram(rejected_sizes),
            "surviving_d5_orbits": len(survivor_sizes),
            "surviving_labelled_supports": sum(survivor_sizes),
            "surviving_orbit_size_histogram": _histogram(survivor_sizes),
            "window_multiplicity_histogram_over_representatives": {
                str(key): occurrence[key] for key in sorted(occurrence)
            },
            "window_multiplicity_histogram_over_labelled_supports": {
                str(key): labelled_occurrence[key]
                for key in sorted(labelled_occurrence)
            },
            "relation_restrictions": total_restrictions,
            "relation_domain_histogram": dict(sorted(domain_histogram.items())),
            "candidates_sha256": _digest(candidates),
            "candidates": candidates,
        },
        "signed_frontier": _signed_frontier(),
        "claim": {
            "exact_scope": (
                "necessary R(3,3) support filter for Ramsey-good signed "
                "completions of the f=5 srg(45,22,10,11) involution branch"
            ),
            "signed_completion_status": SIGNED_STATUS,
            "strongly_regular_fixed_five_branch_only": True,
            "general_ramsey_bound_claimed": False,
        },
    }


def _compare_exact(actual, expected, path="root") -> None:
    if type(actual) is not type(expected):
        raise CheckViolation(
            f"{path}: type {type(actual).__name__} != "
            f"{type(expected).__name__}")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise CheckViolation(f"{path}: object keys differ")
        for key in expected:
            _compare_exact(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise CheckViolation(f"{path}: list lengths differ")
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            _compare_exact(left, right, f"{path}[{index}]")
    elif actual != expected:
        raise CheckViolation(f"{path}: value differs")


def verify_document(document: dict, input_path: Path = DEFAULT_INPUT) -> dict:
    """Reconstruct and exactly compare the complete frontier artifact."""
    _require(document, dict, "artifact root")
    expected = _build_expected(input_path)
    _compare_exact(document, expected)
    claim = document["claim"]
    census = document["census"]
    return {
        "rejected_d5_orbits": census["rejected_d5_orbits"],
        "surviving_d5_orbits": census["surviving_d5_orbits"],
        "surviving_labelled_supports": census["surviving_labelled_supports"],
        "signed_completion_status": claim["signed_completion_status"],
        "general_ramsey_bound_claimed": claim[
            "general_ramsey_bound_claimed"],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", nargs="?", type=Path,
                        default=DEFAULT_ARTIFACT)
    args = parser.parse_args(argv)
    document, _ = _read_strict(args.artifact)
    summary = verify_document(document)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
