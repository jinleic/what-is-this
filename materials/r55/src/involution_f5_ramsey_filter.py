#!/usr/bin/env python3
"""Exact R(3,3) support filter for the fixed-five involution branch.

The input is the exact balanced-support census.  If a fixed C5 edge and its
opposite fixed nonedge select at least three common transposition orbits, the
resulting six vertices contain a triangle or an independent triple.  Together
with the fixed pair this forces a K5 or I5.  The filter is necessary for a
Ramsey-good completion; it does not solve signed SRG completion.
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
DEFAULT_OUTPUT = ROOT / "data" / "involution_f5_ramsey_filter.json"
POINTS = tuple(range(5))
OPPOSITE_WINDOWS = (
    ((0, 1), (2, 4)),
    ((1, 2), (0, 3)),
    ((2, 3), (1, 4)),
    ((3, 4), (0, 2)),
    ((0, 4), (1, 3)),
)


class FrontierViolation(RuntimeError):
    """An input invariant or exact frontier count failed."""


def _dihedral_permutations() -> tuple[tuple[int, ...], ...]:
    rotations = [tuple((vertex + shift) % 5 for vertex in POINTS)
                 for shift in POINTS]
    reflections = [tuple((shift - vertex) % 5 for vertex in POINTS)
                   for shift in POINTS]
    permutations = tuple(dict.fromkeys(rotations + reflections))
    if len(permutations) != 10:
        raise FrontierViolation("failed to construct D5")
    return permutations


DIHEDRAL_PERMUTATIONS = _dihedral_permutations()


def permute_mask(mask: int, permutation: tuple[int, ...]) -> int:
    """Apply a point permutation to one five-bit support mask."""
    if type(mask) is not int or not 0 <= mask < 32:
        raise FrontierViolation("support mask must be an integer in 0..31")
    if (type(permutation) is not tuple or len(permutation) != 5
            or set(permutation) != set(POINTS)):
        raise FrontierViolation("invalid point permutation")
    return sum(1 << permutation[vertex] for vertex in POINTS
               if mask >> vertex & 1)


def _reject_constant(token: str):
    raise FrontierViolation(f"nonstandard JSON constant {token}")


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise FrontierViolation(f"duplicate JSON member {key!r}")
        result[key] = value
    return result


def _load_json_bytes(raw: bytes) -> dict:
    try:
        text = raw.decode("utf-8")
        document = json.loads(
            text,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise FrontierViolation(f"invalid support artifact JSON: {exc}") from exc
    if type(document) is not dict:
        raise FrontierViolation("support artifact root must be an object")
    return document


def _require_exact(value, expected_type: type, label: str) -> None:
    if type(value) is not expected_type:
        raise FrontierViolation(
            f"{label} must have exact type {expected_type.__name__}")


def load_support_artifact(path: Path = DEFAULT_INPUT) -> tuple[dict, bytes]:
    """Strictly load the byte-pinned exact balanced-support dependency."""
    if path.resolve() != DEFAULT_INPUT.resolve():
        raise FrontierViolation(
            "only the canonical balanced-support artifact is accepted")
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise FrontierViolation(f"cannot read support artifact: {exc}") from exc
    digest = hashlib.sha256(raw).hexdigest()
    if len(raw) != INPUT_BYTES or digest != INPUT_SHA256:
        raise FrontierViolation("balanced-support dependency commitment changed")
    document = _load_json_bytes(raw)
    if type(document.get("schema_version")) is not int:
        raise FrontierViolation("schema_version must have exact type int")
    if document["schema_version"] != 1:
        raise FrontierViolation("unexpected support schema version")
    if document.get("campaign_id") != "involution_f5_support_census":
        raise FrontierViolation("unexpected support campaign")
    if document.get("disposition") != INPUT_DISPOSITION:
        raise FrontierViolation("unexpected support disposition")
    census = document.get("support_census")
    _require_exact(census, dict, "support_census")
    if census.get("balanced_labelled_designs") != 7872:
        raise FrontierViolation("unexpected labelled support count")
    if census.get("dihedral_orbits") != 844:
        raise FrontierViolation("unexpected D5 support count")
    if census.get("orbit_size_histogram") != {"1": 2, "5": 110, "10": 732}:
        raise FrontierViolation("unexpected support orbit-size histogram")
    if census.get("representatives_sha256") != INPUT_REPRESENTATIVES_SHA256:
        raise FrontierViolation("unexpected support representative commitment")
    representatives = census.get("representatives")
    _require_exact(representatives, list, "support representatives")
    if len(representatives) != 844:
        raise FrontierViolation("expected 844 support representatives")
    return document, raw


def _count_vector(value, label: str) -> tuple[int, ...]:
    _require_exact(value, list, label)
    if len(value) != 16:
        raise FrontierViolation(f"{label} must have length sixteen")
    counts = []
    for index, count in enumerate(value):
        _require_exact(count, int, f"{label}[{index}]")
        if count < 0:
            raise FrontierViolation(f"{label}[{index}] is negative")
        counts.append(count)
    if sum(counts) != 10:
        raise FrontierViolation(f"{label} must sum to ten")
    return tuple(counts)


def expand_representative(census: dict, representative: dict
                          ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Expand one canonical count record in the declared deterministic order."""
    _require_exact(census, dict, "support census")
    _require_exact(representative, dict, "representative")
    even_order = census.get("edge_mask_order")
    odd_order = census.get("nonedge_mask_order")
    _require_exact(even_order, list, "edge_mask_order")
    _require_exact(odd_order, list, "nonedge_mask_order")
    if (len(even_order) != 16 or len(odd_order) != 16
            or any(type(mask) is not int for mask in even_order + odd_order)):
        raise FrontierViolation("mask orders must contain sixteen integers each")
    if any(mask.bit_count() % 2 for mask in even_order):
        raise FrontierViolation("edge mask order contains an odd support")
    if any(mask.bit_count() % 2 != 1 for mask in odd_order):
        raise FrontierViolation("nonedge mask order contains an even support")
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
        raise FrontierViolation("representative did not expand to twenty orbits")
    return tuple(masks), tuple(internal)


def _window_orbits(masks: tuple[int, ...], edge: tuple[int, int],
                    nonedge: tuple[int, int]) -> tuple[int, ...]:
    return tuple(
        orbit for orbit, mask in enumerate(masks)
        if all(mask >> vertex & 1 for vertex in edge)
        and all(not (mask >> vertex & 1) for vertex in nonedge)
    )


def opposite_window_multiplicities(
        masks: tuple[int, ...]) -> tuple[int, ...]:
    """Count transposition orbits in the five opposite edge/nonedge windows."""
    if type(masks) is not tuple or len(masks) != 20:
        raise FrontierViolation("expected a tuple of twenty support masks")
    if any(type(mask) is not int or not 0 <= mask < 32 for mask in masks):
        raise FrontierViolation("invalid support mask")
    return tuple(len(_window_orbits(masks, edge, nonedge))
                 for edge, nonedge in OPPOSITE_WINDOWS)


def four_vertex_allowed_w_values(left_internal: int,
                                 right_internal: int) -> tuple[int, ...]:
    """Allowed W_ij/2 values when four shared vertices avoid K3 and I3."""
    if (type(left_internal) is not int or left_internal not in (0, 1)
            or type(right_internal) is not int
            or right_internal not in (0, 1)):
        raise FrontierViolation("internal relation bits must be zero or one")
    if left_internal != right_internal:
        return (0,)
    return (-1, 0) if left_internal == 0 else (0, 1)


def r33_relation_restrictions(
        masks: tuple[int, ...], internal: tuple[int, ...]) -> list[dict]:
    """Derive all four-vertex relation domains from size-two intersections."""
    if type(internal) is not tuple or len(internal) != 20:
        raise FrontierViolation("expected twenty internal relation bits")
    if any(type(bit) is not int or bit not in (0, 1) for bit in internal):
        raise FrontierViolation("invalid internal relation bit")
    restrictions = []
    seen = set()
    for window_index, (edge, nonedge) in enumerate(OPPOSITE_WINDOWS):
        shared = _window_orbits(masks, edge, nonedge)
        if len(shared) != 2:
            continue
        left, right = shared
        if (left, right) in seen:
            raise FrontierViolation("duplicate R(3,3) relation restriction")
        seen.add((left, right))
        restrictions.append({
            "left_orbit": left,
            "right_orbit": right,
            "source_window_index": window_index,
            "w_ij_over_2": list(four_vertex_allowed_w_values(
                internal[left], internal[right])),
        })
    return restrictions


def signed_frontier_model() -> dict:
    """Record exact signed consequences exposed by the new support frontier."""
    sign_degrees = {}
    for support_size in range(6):
        plus = 3 + support_size // 2
        sign_degrees[str(support_size)] = {
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
        "w_sign_degrees_by_support_size": sign_degrees,
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


def _histogram(values, keys=(1, 5, 10)) -> dict[str, int]:
    counts = Counter(values)
    return {str(key): counts[key] for key in keys}


def _digest_records(records: list[dict]) -> str:
    digest = hashlib.sha256()
    for record in records:
        digest.update(json.dumps(
            record, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def build_analysis(input_path: Path = DEFAULT_INPUT) -> dict:
    """Filter all 844 exact support representatives and build the artifact."""
    support_document, raw = load_support_artifact(input_path)
    support_census = support_document["support_census"]
    candidates = []
    rejected_orbit_sizes = []
    surviving_orbit_sizes = []
    occurrence_histogram = Counter()
    labelled_occurrence_histogram = Counter()
    restriction_domain_histogram = Counter()
    total_restrictions = 0

    for expected_index, representative in enumerate(
            support_census["representatives"]):
        _require_exact(representative.get("index"), int,
                       f"representative[{expected_index}].index")
        if representative["index"] != expected_index:
            raise FrontierViolation("representative indices are not contiguous")
        _require_exact(representative.get("orbit_size"), int,
                       f"representative[{expected_index}].orbit_size")
        orbit_size = representative["orbit_size"]
        if orbit_size not in (1, 5, 10):
            raise FrontierViolation("invalid D5 orbit size")
        masks, internal = expand_representative(support_census, representative)
        multiplicities = opposite_window_multiplicities(masks)
        for multiplicity in multiplicities:
            occurrence_histogram[multiplicity] += 1
            labelled_occurrence_histogram[multiplicity] += orbit_size
        violating = [index for index, value in enumerate(multiplicities)
                     if value >= 3]
        if violating:
            status = "REJECTED_FORCED_R33_SIX_SET"
            restrictions = []
            rejected_orbit_sizes.append(orbit_size)
        else:
            status = "SURVIVES_R33_SUPPORT_FILTER"
            restrictions = r33_relation_restrictions(masks, internal)
            surviving_orbit_sizes.append(orbit_size)
            total_restrictions += len(restrictions)
            for restriction in restrictions:
                key = ",".join(map(str, restriction["w_ij_over_2"]))
                restriction_domain_histogram[key] += 1
        candidates.append({
            "source_index": expected_index,
            "orbit_size": orbit_size,
            "opposite_window_multiplicities": list(multiplicities),
            "status": status,
            "violating_window_indices": violating,
            "relation_restrictions": restrictions,
        })

    rejected_labelled = sum(rejected_orbit_sizes)
    surviving_labelled = sum(surviving_orbit_sizes)
    if (len(rejected_orbit_sizes), rejected_labelled,
            len(surviving_orbit_sizes), surviving_labelled) != (
                139, 1245, 705, 6627):
        raise FrontierViolation("unexpected exact R(3,3) frontier counts")
    if total_restrictions != 1227:
        raise FrontierViolation("unexpected relation-restriction count")

    return {
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "disposition": DISPOSITION,
        "input_dependency": {
            "relative_path": "r55/data/involution_f5_support_census.json",
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "disposition": support_document["disposition"],
            "representatives_sha256": support_census[
                "representatives_sha256"],
        },
        "theorem": {
            "fixed_graph": "C5",
            "ramsey_fact": "R(3,3)=6",
            "opposite_windows": [
                {"fixed_edge": list(edge), "fixed_nonedge": list(nonedge)}
                for edge, nonedge in OPPOSITE_WINDOWS
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
            "input_orbit_size_histogram": support_census[
                "orbit_size_histogram"],
            "rejected_d5_orbits": len(rejected_orbit_sizes),
            "rejected_labelled_supports": rejected_labelled,
            "rejected_orbit_size_histogram": _histogram(
                rejected_orbit_sizes),
            "surviving_d5_orbits": len(surviving_orbit_sizes),
            "surviving_labelled_supports": surviving_labelled,
            "surviving_orbit_size_histogram": _histogram(
                surviving_orbit_sizes),
            "window_multiplicity_histogram_over_representatives": {
                str(key): occurrence_histogram[key]
                for key in sorted(occurrence_histogram)
            },
            "window_multiplicity_histogram_over_labelled_supports": {
                str(key): labelled_occurrence_histogram[key]
                for key in sorted(labelled_occurrence_histogram)
            },
            "relation_restrictions": total_restrictions,
            "relation_domain_histogram": dict(sorted(
                restriction_domain_histogram.items())),
            "candidates_sha256": _digest_records(candidates),
            "candidates": candidates,
        },
        "signed_frontier": signed_frontier_model(),
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


@functools.cache
def _analysis_json() -> str:
    return json.dumps(
        build_analysis(), sort_keys=True, separators=(",", ":"))


def run_analysis() -> dict:
    """Return a fresh copy of the exact default frontier artifact."""
    return json.loads(_analysis_json())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    document = build_analysis()
    payload = json.dumps(document, indent=2, sort_keys=True) + "\n"
    args.output.write_text(payload)
    census = document["census"]
    print(
        f"rejected={census['rejected_d5_orbits']} "
        f"surviving={census['surviving_d5_orbits']} "
        f"labelled_surviving={census['surviving_labelled_supports']} "
        f"disposition={document['disposition']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
