#!/usr/bin/env python3
"""Disjoint checker for the mixed-tier Engstrom identity artifact.

Role and trust boundary
-----------------------
This module re-verifies ``data/engstrom_identity.json`` (schema version 3)
without importing any cone, identity, search, or producer module.  Its only
non-standard-library import is the frozen graph6 parser and bit counter from
``check_ramsey``.  It independently validates the artifact schema and trust
roots, hashes every recorded input, reconstructs all 3,215 cone states,
checks the eight-coefficient dual certificates and primal witnesses exactly,
and derives each route status and the terminal disposition.  A negative
terminal disposition says only that the recorded frozen continuous basis did
not prove one of its acceptance routes; this checker makes no Ramsey-number
bound claim.

Independent motif kernel
------------------------
The fifteen streamed coordinates are counted by a pair/triangle-incidence
kernel, not by the producer's closed forms.  For every vertex pair, the
kernel counts edges and non-edges inside its common-neighbour set.  Adjacent
pairs account for K4s; non-adjacent pairs account for diamonds and opposite
pairs of C4s.  The same incidence pass on the complement gives i4, 2K2, and
K2+2K1.  A separately ordered triangle pass counts the unique triangle of
K3+K1 and the paw; applying it to the complement gives claws and P3+K1.
P4 is the remaining member of the eleven induced four-vertex classes.
Triple incidences give t, P3, complement-P3, and q.  ``brute_motif_vector``
classifies induced triples and quadruples directly with
``itertools.combinations``; ``cross_validate_motif_kernel(6)`` compares every
coordinate on all 33,867 labeled graphs of orders 1 through 6.

Row reconstruction
------------------
For g, integer degree histograms are reconstructed independently.  A dynamic
program carries the R1--R3 triangle-incidence lower/upper sums modulo three,
the capped interval gap, and the extrema of the two linear q objectives.  It
is exactly equivalent to enumerating every histogram but does not spend an
hour materializing them.  The g endpoints then follow the signed m=3 row.
The h endpoints are independently evaluated from

    3*(47-2*d)*t(X) + 4*diamond(X) - 6*(K3+K1)(X) - 12*i4(Z),

and F is independently interval-evaluated from the full mixed row printed in
the campaign plan.  No g or h value is copied from or compared with the
frozen m=4 state table.

Scoped gap (verbatim checker disclosure)
----------------------------------------
For classes represented by a hash-pinned catalog, ``--sweep`` recomputes all
fifteen motif windows with the independent kernel.  For missing classes, the
checker independently recomputes the m=3 q envelope from integer degree
histograms and R(3,5) edge windows, but does not rerun the producer's
degree-histogram p3/pc3 derivation or exact R1-R5 envelope LP for the
remaining four-motif windows.  It validates those recorded outer windows for
exact schema, combinatorial caps, and containment of every applicable catalog
class; h and F endpoints in missing strata
therefore remain conditional on those hash-pinned windows.  This gap does not
turn the recorded negative disposition into a bound: the witnesses only
certify non-provability over the recorded frozen cone.  Any future
MIXED_ACCEPTED_CUT requires an independent R1-R5 LP recomputation before an
evidence-grade positive verdict.

CLI
---
``check_mixed_certificate.py [--root DIR] [--no-input-hashes] [--sweep]
ANALYSIS`` exits 0 and ends stdout with
``MIXED EVIDENCE VERIFIED: <disposition>``.  A mismatch exits 1 and prints
``MIXED EVIDENCE FAILED: <field>: <counterexample>`` to stderr.  The optional
catalog sweep streams 8,500,211 graphs and prints progress at every 500,000
completed graphs.  ``--no-input-hashes`` is an explicit scope reduction for
machines without the large catalogs and always prints a loud
``MIXED PROVENANCE UNPROVED`` notice.
"""

import argparse
import hashlib
import itertools
import json
import multiprocessing
import os
import sys
from fractions import Fraction
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_ramsey import parse_graph6_line, popcount  # noqa: E402


# ---------------------------------------------------------------- constants

N_TARGET = 45
VERTICES = Fraction(N_TARGET)
SCHEMA_VERSION = 3
CAMPAIGN_ID = "higher_order_identity_positive_deficiency_mixed"
EXPECTED_STATES = 3215
EXPECTED_CATALOG_GRAPHS = 8_500_211

STREAM_KEYS = (
    "q", "t", "p3", "pc3", "i4", "k4", "diamond", "k3k1", "c4",
    "claw", "paw", "p4", "two_k2", "k2_2k1", "p3_k1",
)
FOUR_KEYS = (
    "i4", "k2_2k1", "two_k2", "p3_k1", "p4", "claw", "k3k1",
    "c4", "paw", "diamond", "k4",
)
STATE_DEGREES = tuple(range(20, 25))
FROZEN_EDGE_BOUNDS = {
    17: (41, 79), 18: (50, 85), 19: (57, 92), 20: (68, 100),
    21: (77, 107), 22: (88, 114), 23: (101, 122), 24: (116, 132),
}
E45 = {order: FROZEN_EDGE_BOUNDS[order][1] for order in STATE_DEGREES}
R35_COUNTS = (1, 2, 3, 7, 13, 32, 71, 179, 290, 313, 105, 12, 1)

EXTREME_CLASSES = (
    (17, 41), (17, 42), (17, 77), (17, 78), (17, 79),
    (18, 50), (18, 51), (18, 83), (18, 84), (18, 85),
    (19, 57), (19, 58), (19, 90), (19, 91), (19, 92),
    (20, 68), (20, 69), (20, 98), (20, 99), (20, 100),
    (21, 77), (21, 78), (21, 106), (21, 107),
    (22, 88), (22, 89), (22, 113), (22, 114),
    (23, 101), (23, 102), (23, 103), (23, 104),
    (23, 119), (23, 120), (23, 121), (23, 122),
)
CATALOG_CLASSES = EXTREME_CLASSES + tuple((24, edge) for edge in range(116, 133))
OUTER_CLASSES = tuple(
    (order, edge)
    for order in STATE_DEGREES
    for edge in range(FROZEN_EDGE_BOUNDS[order][0],
                      FROZEN_EDGE_BOUNDS[order][1] + 1)
)

PRODUCER_SOURCES = (
    "check_ramsey.py",
    "m3_deficiency_cone.py",
    "m4_deficiency_cone.py",
    "m4_subgraph_identities.py",
    "mixed_deficiency_cone.py",
    "mixed_subgraph_identities.py",
    "search_mixed_cuts.py",
    "subgraph_identities.py",
)
STATIC_INPUTS = (
    "r55/data/VALIDATION.json",
    "r55/data/VALIDATION_extreme.json",
    "r55/data/higher_identity_m3.json",
    "r55/data/higher_identity_m4.json",
    "r55/data/structural_tables.json",
)

# Each statement is pinned independently as its UTF-8 SHA-256; URLs and order
# are pinned literally.  This catches any label wording or source drift while
# keeping this checker independent of producer constants.
TRUST_ROOT_PINS = (
    ("94331103c4fe6bf493ad4eb4f52ebacc7427ab43a4262484eb9633c929c475ae",
     "https://users.cecs.anu.edu.au/~bdm/data/ramsey.html"),
    ("bb81b903bc9abaeb84ae52f96d6bdbaa48c1c6eca734100ebf59052878a0ccfb",
     "https://arxiv.org/abs/2409.15709"),
    ("bed3bb471f4dcaa170ca1368d767c6a2e9cd2867a43057dc161b5bdee88c4fb1",
     "https://arxiv.org/abs/1703.08768"),
    ("4d5f1f5a60b15e7d876d22a34658a5a365cc8f4ba84f79d9f36a23a25e5f07d9",
     "https://users.cecs.anu.edu.au/~bdm/data/ramsey.html"),
    ("85c6082a2fe7cd074213cc75267586280165e15bab81309661ff6b79a79242da",
     "https://users.cecs.anu.edu.au/~bdm/papers/r55.pdf"),
    ("6e1a2626f115638f172242b6abc9bfd4203f68ce7bb31d9f821e2b79a7f2af50",
     "https://arxiv.org/abs/2409.15709"),
    ("c48718b9eebae4bb36275d1eb3707c87a7f6d03ef68d25b1eb936a265e30a776",
     "https://arxiv.org/abs/1002.4304"),
    ("07336f960e71ca2024445425ab4e309c6bcc98c5418ef07c4ba41afba78f16e2",
     "https://arxiv.org/abs/2409.15709"),
)

ACCEPTED = "ACCEPTED_EXACT_CUT"
REJECTED = "REJECTED_BY_EXACT_WITNESS"
UNRESOLVED = "CERTIFICATION_UNRESOLVED"
UNAVAILABLE = "UNAVAILABLE_NO_COVER_CERTIFICATE"
ROUTE_IDS = ("total_deficiency", "degree20_count", "deficiency_ge8_count")
UNAVAILABLE_ROUTE = "required_local_family"
ROUTE_EDGES = {
    "total_deficiency": (Fraction(315), True),
    "degree20_count": (Fraction(1), False),
    "deficiency_ge8_count": (Fraction(1), False),
}
COEFFICIENTS = (
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
)
DISPOSITION_ACCEPTED = "MIXED_ACCEPTED_CUT"
DISPOSITION_NO_CUT = "MIXED_NO_CUT_IN_FROZEN_BASIS"
DISPOSITION_UNRESOLVED = "MIXED_CERTIFICATION_UNRESOLVED"
DISPOSITIONS = (
    DISPOSITION_ACCEPTED, DISPOSITION_NO_CUT, DISPOSITION_UNRESOLVED,
)

REPLAY_RECORD = {"published_graphs": 328, "complements": 328, "residual_zero": 656}
N49_RECORD = {
    "extremal_graphs": 2,
    "combos": 4,
    "row_values": [0, 144, 288, 432],
    "zero_corners": 1,
}
UNAVAILABLE_RECORD = {
    "objective_id": UNAVAILABLE_ROUTE,
    "certificate": None,
    "exact_upper_bound": None,
    "primal_witness": None,
    "exact_witness_value": None,
    "route_status": UNAVAILABLE,
}

_TOP_LEVEL_KEYS = {
    "campaign_id", "catalog_motif_histograms", "catalog_motif_windows",
    "disposition", "inputs", "n49_replay", "outer_motif_windows",
    "r35_edge_windows", "replay", "schema_version", "searches", "states",
    "trust_roots",
}
_INPUT_KEYS = {"bytes", "graph_count", "relative_path", "sha256"}
_TRUST_KEYS = {"statement", "url"}
_R35_KEYS = {"edge_max", "edge_min", "order"}
_HISTOGRAM_KEYS = {"edges", "graph_count", "histograms", "order", "source"}
_WINDOW_KEYS = {"edges", "order", "source", "windows"}
_OUTER_KEYS = {"edges", "method", "order", "windows"}
_HIST_ENTRY_KEYS = {"count", "value"}
_STATE_KEYS = {
    "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
    "h_lo", "h_hi", "f_lo", "f_hi", "x_interval_source",
    "y_interval_source", "x_motif_source", "y_motif_source",
}
_STATE_INT_KEYS = (
    "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
    "h_lo", "h_hi", "f_lo", "f_hi",
)
_SEARCH_KEYS = {
    "certificate", "exact_upper_bound", "exact_witness_value", "objective_id",
    "primal_witness", "route_status",
}
_CERTIFICATE_KEYS = set(COEFFICIENTS) | {"objective_id"}
_WITNESS_KEYS = {"objective_id", "weights"}
_WEIGHT_KEYS = {"state_index", "weight"}

_DEFAULT_ROOT = Path(__file__).resolve().parents[2]
_STREAM_BLOCK = 1 << 20
_SWEEP_CHUNK = 10_000
_SWEEP_PROGRESS = 500_000
_STATE_CACHE = {}


class CheckViolation(ValueError):
    """One artifact field disagrees with independently derived exact data."""


class InfeasibleEdgeClass(ValueError):
    """No integer degree histogram survives for an edge class."""


# -------------------------------------------------------------- exact guards

def _c2(value):
    return value * (value - 1) // 2


def _c3(value):
    return value * (value - 1) * (value - 2) // 6


def _c4(value):
    return value * (value - 1) * (value - 2) * (value - 3) // 24


def _exact_int(value, label):
    if type(value) is not int:
        raise CheckViolation(f"{label}: expected an exact int, got {value!r}")
    return value


def _exact_string(value, label):
    if not isinstance(value, str) or not value:
        raise CheckViolation(f"{label}: expected a nonempty string, got {value!r}")
    return value


def _require_keys(record, keys, label):
    if not isinstance(record, dict):
        raise CheckViolation(f"{label}: expected an object, got {type(record).__name__}")
    actual = set(record)
    if actual != set(keys):
        raise CheckViolation(
            f"{label}: key set is {sorted(actual)!r}, expected {sorted(keys)!r}"
        )
    return record


def _require_list(value, label):
    if not isinstance(value, list):
        raise CheckViolation(f"{label}: expected a list, got {type(value).__name__}")
    return value


def _canonical_fraction(text, label):
    if not isinstance(text, str):
        raise CheckViolation(f"{label}: expected a fraction string, got {text!r}")
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError):
        raise CheckViolation(f"{label}: invalid rational {text!r}") from None
    if str(value) != text:
        raise CheckViolation(f"{label}: noncanonical rational {text!r}")
    return value


def _is_sha256(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _relative_path(value, label):
    text = _exact_string(value, label)
    pure = PurePosixPath(text)
    if pure.is_absolute() or str(pure) != text or ".." in pure.parts or "." in pure.parts:
        raise CheckViolation(f"{label}: noncanonical repository path {text!r}")
    if not pure.parts or pure.parts[0] != "r55":
        raise CheckViolation(f"{label}: path must live under r55/: {text!r}")
    return text


def _load_json(path, label):
    try:
        return json.loads(Path(path).read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise CheckViolation(f"{label}: unreadable JSON ({exc})") from None


def load_analysis(path):
    """Load a candidate artifact with the checker's fail-closed error type."""
    return _load_json(path, f"analysis: {path}")


def _sha256_of(path):
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(_STREAM_BLOCK), b""):
                digest.update(block)
    except OSError as exc:
        raise CheckViolation(f"input {path}: unreadable ({exc})") from None
    return digest.hexdigest()


# ------------------------------------------------------- independent motifs

def complement(adj):
    """Return the loop-free complement of one bit-mask adjacency list."""
    order = len(adj)
    full = (1 << order) - 1
    return [full & ~row & ~(1 << vertex) for vertex, row in enumerate(adj)]


def _edges_inside(adj, mask):
    """Count graph edges with both endpoints in ``mask`` exactly once."""
    edges = 0
    remaining = mask
    while remaining:
        bit = remaining & -remaining
        vertex = bit.bit_length() - 1
        remaining ^= bit
        edges += popcount(adj[vertex] & remaining)
    return edges


def _pair_incidence(adj):
    """Return (diamond, C4, K4) from common-neighbour pair incidences."""
    diamonds = c4_twice = k4_sixfold = 0
    order = len(adj)
    for left in range(order):
        for right in range(left + 1, order):
            common = adj[left] & adj[right]
            common_size = popcount(common)
            linked = _edges_inside(adj, common)
            if (adj[left] >> right) & 1:
                k4_sixfold += linked
            else:
                diamonds += linked
                c4_twice += _c2(common_size) - linked
    if k4_sixfold % 6 or c4_twice % 2:
        raise ValueError("pair-incidence divisibility invariant failed")
    return diamonds, c4_twice // 2, k4_sixfold // 6


def _triangle_incidence(adj):
    """Return (triangles, K3+K1, paw) by ordered triangle extensions."""
    order = len(adj)
    full = (1 << order) - 1
    triangles = k3k1 = paws = 0
    for first in range(order):
        later_first = adj[first] & ~((1 << (first + 1)) - 1) & full
        while later_first:
            second_bit = later_first & -later_first
            second = second_bit.bit_length() - 1
            later_first ^= second_bit
            thirds = (
                adj[first] & adj[second]
                & ~((1 << (second + 1)) - 1) & full
            )
            while thirds:
                third_bit = thirds & -thirds
                third = third_bit.bit_length() - 1
                thirds ^= third_bit
                triangles += 1
                first_row = adj[first]
                second_row = adj[second]
                third_row = adj[third]
                isolated = full & ~first_row & ~second_row & ~third_row
                k3k1 += popcount(isolated)
                exactly_one = (
                    (first_row & ~second_row & ~third_row)
                    | (second_row & ~first_row & ~third_row)
                    | (third_row & ~first_row & ~second_row)
                ) & full
                paws += popcount(exactly_one)
    return triangles, k3k1, paws


def motif_vector(adj):
    """Count all fifteen streamed coordinates with the independent kernel."""
    order = len(adj)
    degrees = [popcount(row) for row in adj]
    wedges = sum(_c2(degree) for degree in degrees)
    dual = complement(adj)
    dual_wedges = sum(_c2(popcount(row)) for row in dual)

    triangles, k3k1, paw = _triangle_incidence(adj)
    dual_triangles, claw, p3_k1 = _triangle_incidence(dual)
    diamond, c4, k4 = _pair_incidence(adj)
    k2_2k1, two_k2, i4 = _pair_incidence(dual)
    p3 = wedges - 3 * triangles
    pc3 = dual_wedges - 3 * dual_triangles

    four_known = {
        "i4": i4,
        "k2_2k1": k2_2k1,
        "two_k2": two_k2,
        "p3_k1": p3_k1,
        "claw": claw,
        "k3k1": k3k1,
        "c4": c4,
        "paw": paw,
        "diamond": diamond,
        "k4": k4,
    }
    p4 = _c4(order) - sum(four_known.values())
    vector = {
        "q": wedges - triangles,
        "t": triangles,
        "p3": p3,
        "pc3": pc3,
        **four_known,
        "p4": p4,
    }
    if set(vector) != set(STREAM_KEYS) or any(value < 0 for value in vector.values()):
        raise ValueError(f"motif partition invariant failed: {vector}")
    if sum(vector[key] for key in FOUR_KEYS) != _c4(order):
        raise ValueError("four-vertex motif partition does not sum to C(n,4)")
    return {key: vector[key] for key in STREAM_KEYS}


_FOUR_CLASS_BY_DEGREES = {
    (0, 0, 0, 0): "i4",
    (0, 0, 1, 1): "k2_2k1",
    (1, 1, 1, 1): "two_k2",
    (0, 1, 1, 2): "p3_k1",
    (1, 1, 2, 2): "p4",
    (1, 1, 1, 3): "claw",
    (0, 2, 2, 2): "k3k1",
    (2, 2, 2, 2): "c4",
    (1, 2, 2, 3): "paw",
    (2, 2, 3, 3): "diamond",
    (3, 3, 3, 3): "k4",
}


def brute_motif_vector(adj):
    """Slow combinations-only reference, intentionally unrelated to the kernel."""
    order = len(adj)
    degrees = [popcount(row) for row in adj]
    triangles = p3 = pc3 = 0
    for vertices in itertools.combinations(range(order), 3):
        edges = sum(
            (adj[left] >> right) & 1
            for left, right in itertools.combinations(vertices, 2)
        )
        triangles += int(edges == 3)
        p3 += int(edges == 2)
        pc3 += int(edges == 1)

    counts = {key: 0 for key in FOUR_KEYS}
    for vertices in itertools.combinations(range(order), 4):
        local_degrees = tuple(sorted(
            sum((adj[vertex] >> other) & 1 for other in vertices if other != vertex)
            for vertex in vertices
        ))
        try:
            key = _FOUR_CLASS_BY_DEGREES[local_degrees]
        except KeyError:
            raise AssertionError(f"unknown four-vertex degree type {local_degrees}") from None
        counts[key] += 1
    result = {
        "q": sum(_c2(degree) for degree in degrees) - triangles,
        "t": triangles,
        "p3": p3,
        "pc3": pc3,
        **counts,
    }
    return {key: result[key] for key in STREAM_KEYS}


def _adjacency_from_mask(order, mask):
    adjacency = [0] * order
    position = 0
    for right in range(1, order):
        for left in range(right):
            if (mask >> position) & 1:
                adjacency[left] |= 1 << right
                adjacency[right] |= 1 << left
            position += 1
    return adjacency


def cross_validate_motif_kernel(max_order=6):
    """Exhaust all labeled graphs of orders 1..``max_order``; return count."""
    _exact_int(max_order, "max_order")
    if not 1 <= max_order <= 6:
        raise ValueError("cross-validation max_order must lie in 1..6")
    checked = 0
    for order in range(1, max_order + 1):
        for mask in range(1 << _c2(order)):
            adjacency = _adjacency_from_mask(order, mask)
            fast = motif_vector(adjacency)
            reference = brute_motif_vector(adjacency)
            if fast != reference:
                differing = next(key for key in STREAM_KEYS if fast[key] != reference[key])
                raise AssertionError(
                    f"n={order} mask={mask} key={differing}: "
                    f"kernel={fast[differing]} brute={reference[differing]}"
                )
            checked += 1
    return checked


# ---------------------------------------------------------- schema and tables

def _catalog_source(cls):
    order, edges = cls
    if order == 24:
        return "r55/data/r45_24.g6"
    return f"r55/data/r45extreme/r45{order}.{edges}.g6"


def _expected_inputs():
    paths = set(STATIC_INPUTS)
    paths.update(f"r55/data/r35_{order}.g6" for order in range(1, 14))
    paths.update(_catalog_source(cls) for cls in CATALOG_CLASSES)
    paths.add("r55/data/r55_42some.g6")
    paths.update(f"r55/src/{name}" for name in PRODUCER_SOURCES)
    return paths


def _motif_cap(order, key):
    if key == "q":
        return 2 * _c3(order)
    if key in ("t", "p3", "pc3"):
        return _c3(order)
    return _c4(order)


def _window_pair(value, label, cap):
    pair = _require_list(value, label)
    if len(pair) != 2:
        raise CheckViolation(f"{label}: expected [low, high], got {pair!r}")
    low = _exact_int(pair[0], f"{label}[0]")
    high = _exact_int(pair[1], f"{label}[1]")
    # Envelope relaxations may extend below the nonnegative motif orthant;
    # that only widens them.  A negative upper endpoint would exclude every
    # realizable nonnegative value and is therefore not sound.
    if low > high or high < 0 or high > cap:
        raise CheckViolation(f"{label}: unsound or inverted window [{low}, {high}]")
    return low, high


def _check_trust_roots(document):
    roots = _require_list(document["trust_roots"], "trust_roots")
    if len(roots) != len(TRUST_ROOT_PINS):
        raise CheckViolation(
            f"trust_roots: expected {len(TRUST_ROOT_PINS)}, got {len(roots)}"
        )
    observed = []
    for index, record in enumerate(roots):
        _require_keys(record, _TRUST_KEYS, f"trust_roots[{index}]")
        statement = _exact_string(record["statement"], f"trust_roots[{index}].statement")
        url = _exact_string(record["url"], f"trust_roots[{index}].url")
        observed.append((hashlib.sha256(statement.encode("utf-8")).hexdigest(), url))
    if tuple(observed) != TRUST_ROOT_PINS:
        mismatch = next(
            index for index, pair in enumerate(observed)
            if index >= len(TRUST_ROOT_PINS) or pair != TRUST_ROOT_PINS[index]
        )
        raise CheckViolation(
            f"trust_roots[{mismatch}]: statement label or URL drifted"
        )


def _catalog_tables(document):
    histogram_records = _require_list(
        document["catalog_motif_histograms"], "catalog_motif_histograms"
    )
    if len(histogram_records) != len(CATALOG_CLASSES):
        raise CheckViolation(
            "catalog_motif_histograms: expected "
            f"{len(CATALOG_CLASSES)} classes, got {len(histogram_records)}"
        )
    class_counts = {}
    owners = {}
    derived_windows = {}
    observed_classes = []
    total_graphs = 0
    for index, record in enumerate(histogram_records):
        label = f"catalog_motif_histograms[{index}]"
        _require_keys(record, _HISTOGRAM_KEYS, label)
        order = _exact_int(record["order"], f"{label}.order")
        edges = _exact_int(record["edges"], f"{label}.edges")
        graph_count = _exact_int(record["graph_count"], f"{label}.graph_count")
        if graph_count <= 0:
            raise CheckViolation(f"{label}.graph_count: must be positive")
        cls = (order, edges)
        observed_classes.append(cls)
        source = _relative_path(record["source"], f"{label}.source")
        if source != _catalog_source(cls):
            raise CheckViolation(f"{label}.source: wrong owner for class {cls}")
        histograms = _require_keys(record["histograms"], STREAM_KEYS,
                                   f"{label}.histograms")
        windows = {}
        for key in STREAM_KEYS:
            entries = _require_list(histograms[key], f"{label}.histograms.{key}")
            if not entries:
                raise CheckViolation(f"{label}.histograms.{key}: empty histogram")
            previous = None
            subtotal = 0
            for entry_index, entry in enumerate(entries):
                entry_label = f"{label}.histograms.{key}[{entry_index}]"
                _require_keys(entry, _HIST_ENTRY_KEYS, entry_label)
                value = _exact_int(entry["value"], f"{entry_label}.value")
                count = _exact_int(entry["count"], f"{entry_label}.count")
                if value < 0 or value > _motif_cap(order, key):
                    raise CheckViolation(f"{entry_label}.value: outside combinatorial cap")
                if count <= 0:
                    raise CheckViolation(f"{entry_label}.count: must be positive")
                if previous is not None and value <= previous:
                    raise CheckViolation(f"{entry_label}.value: histogram is not increasing")
                previous = value
                subtotal += count
            if subtotal != graph_count:
                raise CheckViolation(
                    f"{label}.histograms.{key}: counts sum to {subtotal}, "
                    f"not graph_count {graph_count}"
                )
            windows[key] = (entries[0]["value"], entries[-1]["value"])
        if windows["k4"] != (0, 0):
            raise CheckViolation(f"{label}.histograms.k4: catalog is not K4-free")
        class_counts[cls] = graph_count
        owners[cls] = source
        derived_windows[cls] = windows
        total_graphs += graph_count

    if tuple(observed_classes) != CATALOG_CLASSES:
        raise CheckViolation(
            "catalog_motif_histograms: class set/order differs from the frozen 53 classes"
        )
    if total_graphs != EXPECTED_CATALOG_GRAPHS:
        raise CheckViolation(
            f"catalog_motif_histograms: total {total_graphs}, "
            f"expected {EXPECTED_CATALOG_GRAPHS}"
        )

    window_records = _require_list(
        document["catalog_motif_windows"], "catalog_motif_windows"
    )
    if len(window_records) != len(CATALOG_CLASSES):
        raise CheckViolation(
            f"catalog_motif_windows: expected {len(CATALOG_CLASSES)} classes"
        )
    window_classes = []
    for index, record in enumerate(window_records):
        label = f"catalog_motif_windows[{index}]"
        _require_keys(record, _WINDOW_KEYS, label)
        order = _exact_int(record["order"], f"{label}.order")
        edges = _exact_int(record["edges"], f"{label}.edges")
        cls = (order, edges)
        window_classes.append(cls)
        source = _relative_path(record["source"], f"{label}.source")
        if source != owners.get(cls):
            raise CheckViolation(f"{label}.source: differs from histogram owner")
        windows = _require_keys(record["windows"], STREAM_KEYS, f"{label}.windows")
        claimed = {
            key: _window_pair(windows[key], f"{label}.windows.{key}",
                              _motif_cap(order, key))
            for key in STREAM_KEYS
        }
        if claimed != derived_windows.get(cls):
            raise CheckViolation(f"{label}.windows: differs from its histograms")
    if tuple(window_classes) != CATALOG_CLASSES:
        raise CheckViolation("catalog_motif_windows: class set/order drifted")
    return class_counts, owners, derived_windows


def _outer_tables(document, catalog_windows):
    records = _require_list(document["outer_motif_windows"], "outer_motif_windows")
    if len(records) != len(OUTER_CLASSES):
        raise CheckViolation(
            f"outer_motif_windows: expected {len(OUTER_CLASSES)} classes, "
            f"got {len(records)}"
        )
    outer = {}
    methods = {}
    classes = []
    for index, record in enumerate(records):
        label = f"outer_motif_windows[{index}]"
        _require_keys(record, _OUTER_KEYS, label)
        order = _exact_int(record["order"], f"{label}.order")
        edges = _exact_int(record["edges"], f"{label}.edges")
        cls = (order, edges)
        classes.append(cls)
        method = _exact_string(record["method"], f"{label}.method")
        expected_method = "catalog" if cls in catalog_windows else "envelope_lp"
        if method != expected_method:
            raise CheckViolation(
                f"{label}.method: {method!r}, expected {expected_method!r}"
            )
        windows = _require_keys(record["windows"], STREAM_KEYS, f"{label}.windows")
        resolved = {
            key: _window_pair(windows[key], f"{label}.windows.{key}",
                              _motif_cap(order, key))
            for key in STREAM_KEYS
        }
        if cls in catalog_windows and resolved != catalog_windows[cls]:
            raise CheckViolation(f"{label}.windows: catalog window drifted")
        outer[cls] = resolved
        methods[cls] = method
    if tuple(classes) != OUTER_CLASSES:
        raise CheckViolation("outer_motif_windows: class set/order drifted")
    return outer, methods


def _check_input_manifest(document, class_counts, owners):
    records = _require_list(document["inputs"], "inputs")
    expected_paths = _expected_inputs()
    if len(records) != len(expected_paths):
        raise CheckViolation(
            f"inputs: expected {len(expected_paths)} records, got {len(records)}"
        )
    by_path = {}
    paths = []
    for index, record in enumerate(records):
        label = f"inputs[{index}]"
        _require_keys(record, _INPUT_KEYS, label)
        relative = _relative_path(record["relative_path"], f"{label}.relative_path")
        paths.append(relative)
        size = _exact_int(record["bytes"], f"{label}.bytes")
        if size < 0:
            raise CheckViolation(f"{label}.bytes: negative")
        if not _is_sha256(record["sha256"]):
            raise CheckViolation(f"{label}.sha256: not a lowercase SHA-256")
        graph_count = record["graph_count"]
        if graph_count is not None:
            graph_count = _exact_int(graph_count, f"{label}.graph_count")
            if graph_count <= 0:
                raise CheckViolation(f"{label}.graph_count: must be positive")
        if relative in by_path:
            raise CheckViolation(f"{label}.relative_path: duplicate {relative!r}")
        by_path[relative] = record
    if paths != sorted(paths):
        raise CheckViolation("inputs: records are not sorted by relative_path")
    if set(paths) != expected_paths:
        missing = sorted(expected_paths - set(paths))
        extra = sorted(set(paths) - expected_paths)
        raise CheckViolation(f"inputs: path set drifted; missing={missing}, extra={extra}")

    expected_graph_counts = {
        f"r55/data/r35_{order}.g6": R35_COUNTS[order - 1]
        for order in range(1, 14)
    }
    expected_graph_counts["r55/data/r55_42some.g6"] = 328
    for cls, count in class_counts.items():
        source = owners[cls]
        expected_graph_counts[source] = expected_graph_counts.get(source, 0) + count
    for relative, record in by_path.items():
        expected = expected_graph_counts.get(relative)
        if record["graph_count"] != expected:
            raise CheckViolation(
                f"inputs[{relative}].graph_count: {record['graph_count']!r}, "
                f"expected {expected!r}"
            )
    return by_path


def _check_replay_shapes(document):
    replay = _require_keys(document["replay"], REPLAY_RECORD, "replay")
    for key in REPLAY_RECORD:
        _exact_int(replay[key], f"replay.{key}")
    if replay != REPLAY_RECORD:
        raise CheckViolation(f"replay: expected {REPLAY_RECORD}, got {replay}")

    n49 = _require_keys(document["n49_replay"], N49_RECORD, "n49_replay")
    for key in ("extremal_graphs", "combos", "zero_corners"):
        _exact_int(n49[key], f"n49_replay.{key}")
    values = _require_list(n49["row_values"], "n49_replay.row_values")
    for index, value in enumerate(values):
        _exact_int(value, f"n49_replay.row_values[{index}]")
    if n49 != N49_RECORD:
        raise CheckViolation(f"n49_replay: expected {N49_RECORD}, got {n49}")


def _check_state_shapes(document):
    states = _require_list(document["states"], "states")
    if len(states) != EXPECTED_STATES:
        raise CheckViolation(f"states: expected {EXPECTED_STATES}, got {len(states)}")
    previous = None
    for index, state in enumerate(states):
        label = f"states[{index}]"
        _require_keys(state, _STATE_KEYS, label)
        for field in _STATE_INT_KEYS:
            _exact_int(state[field], f"{label}.{field}")
        if state["d"] not in STATE_DEGREES:
            raise CheckViolation(f"{label}.d: {state['d']} is outside 20..24")
        if state["a"] < 0 or state["b"] < 0:
            raise CheckViolation(f"{label}: negative deficiency coordinate")
        if state["deficiency"] != state["a"] + state["b"]:
            raise CheckViolation(f"{label}.deficiency: not a+b")
        for low, high in (("g_lo", "g_hi"), ("h_lo", "h_hi"),
                          ("f_lo", "f_hi")):
            if state[low] > state[high]:
                raise CheckViolation(f"{label}: inverted {low}/{high} interval")
        for field in (
            "x_interval_source", "y_interval_source", "x_motif_source",
            "y_motif_source",
        ):
            _exact_string(state[field], f"{label}.{field}")
        coordinate = (state["d"], state["a"], state["b"])
        if previous is not None and coordinate <= previous:
            raise CheckViolation(
                f"{label}: state coordinates are duplicate or out of order: {coordinate}"
            )
        previous = coordinate


def _check_search_shapes(document):
    searches = _require_list(document["searches"], "searches")
    expected = ROUTE_IDS + (UNAVAILABLE_ROUTE,)
    if len(searches) != len(expected):
        raise CheckViolation(f"searches: expected four records, got {len(searches)}")
    observed = []
    for index, record in enumerate(searches):
        label = f"searches[{index}]"
        _require_keys(record, _SEARCH_KEYS, label)
        objective = _exact_string(record["objective_id"], f"{label}.objective_id")
        observed.append(objective)
        _exact_string(record["route_status"], f"{label}.route_status")
        for field in ("exact_upper_bound", "exact_witness_value"):
            if record[field] is not None and not isinstance(record[field], str):
                raise CheckViolation(f"{label}.{field}: expected string or null")
        certificate = record["certificate"]
        if certificate is not None:
            _require_keys(certificate, _CERTIFICATE_KEYS, f"{label}.certificate")
            if certificate["objective_id"] != objective:
                raise CheckViolation(f"{label}.certificate.objective_id: drifted")
            for coefficient in COEFFICIENTS:
                if not isinstance(certificate[coefficient], str):
                    raise CheckViolation(
                        f"{label}.certificate.{coefficient}: expected string"
                    )
        witness = record["primal_witness"]
        if witness is not None:
            _require_keys(witness, _WITNESS_KEYS, f"{label}.primal_witness")
            if witness["objective_id"] != objective:
                raise CheckViolation(f"{label}.primal_witness.objective_id: drifted")
            weights = _require_list(witness["weights"], f"{label}.primal_witness.weights")
            if not weights:
                raise CheckViolation(f"{label}.primal_witness.weights: empty")
            for pair_index, pair in enumerate(weights):
                pair_label = f"{label}.primal_witness.weights[{pair_index}]"
                _require_keys(pair, _WEIGHT_KEYS, pair_label)
                _exact_int(pair["state_index"], f"{pair_label}.state_index")
                if not isinstance(pair["weight"], str):
                    raise CheckViolation(f"{pair_label}.weight: expected string")
    if tuple(observed) != expected:
        raise CheckViolation(f"searches: objective order/set drifted: {observed}")


def _check_schema(document):
    _require_keys(document, _TOP_LEVEL_KEYS, "analysis")
    version = _exact_int(document["schema_version"], "schema_version")
    if version != SCHEMA_VERSION:
        raise CheckViolation(f"schema_version: expected {SCHEMA_VERSION}, got {version}")
    if document["campaign_id"] != CAMPAIGN_ID:
        raise CheckViolation(
            f"campaign_id: expected {CAMPAIGN_ID!r}, got {document['campaign_id']!r}"
        )
    if document["disposition"] not in DISPOSITIONS:
        raise CheckViolation(f"disposition: unknown value {document['disposition']!r}")
    _check_trust_roots(document)
    _check_replay_shapes(document)

    r35_records = _require_list(document["r35_edge_windows"], "r35_edge_windows")
    if len(r35_records) != 14:
        raise CheckViolation(f"r35_edge_windows: expected 14 records, got {len(r35_records)}")
    r35 = {}
    for index, record in enumerate(r35_records):
        label = f"r35_edge_windows[{index}]"
        _require_keys(record, _R35_KEYS, label)
        order = _exact_int(record["order"], f"{label}.order")
        low = _exact_int(record["edge_min"], f"{label}.edge_min")
        high = _exact_int(record["edge_max"], f"{label}.edge_max")
        if order != index or low < 0 or low > high:
            raise CheckViolation(f"{label}: malformed order/window ({order}, {low}, {high})")
        r35[order] = (low, high)

    class_counts, owners, catalog_windows = _catalog_tables(document)
    outer_windows, outer_methods = _outer_tables(document, catalog_windows)
    input_records = _check_input_manifest(document, class_counts, owners)
    _check_state_shapes(document)
    _check_search_shapes(document)
    return {
        "r35": r35,
        "class_counts": class_counts,
        "owners": owners,
        "catalog_windows": catalog_windows,
        "outer_windows": outer_windows,
        "outer_methods": outer_methods,
        "inputs": input_records,
    }


# ---------------------------------------------------------- provenance input

def _check_hashes(context, root):
    root = Path(root).resolve()
    for relative in sorted(context["inputs"]):
        record = context["inputs"][relative]
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError:
            raise CheckViolation(f"inputs[{relative}]: resolves outside repository root")
        if not path.is_file():
            raise CheckViolation(f"inputs[{relative}]: missing input file")
        try:
            size = path.stat().st_size
        except OSError as exc:
            raise CheckViolation(f"inputs[{relative}]: cannot stat ({exc})") from None
        if size != record["bytes"]:
            raise CheckViolation(
                f"inputs[{relative}].bytes: observed {size}, recorded {record['bytes']}"
            )
        observed = _sha256_of(path)
        if observed != record["sha256"]:
            raise CheckViolation(
                f"inputs[{relative}].sha256: observed {observed}, "
                f"recorded {record['sha256']}"
            )


def _load_edge_bounds(root):
    path = Path(root) / "r55/data/structural_tables.json"
    document = _load_json(path, "structural_tables.json")
    if not isinstance(document, dict):
        raise CheckViolation("structural_tables.json: top level is not an object")
    table = document.get("extremal_edge_tables_R45")
    if not isinstance(table, dict):
        raise CheckViolation("structural_tables.json: missing extremal_edge_tables_R45")
    observed = {}
    for order in range(17, 25):
        record = table.get(str(order))
        _require_keys(record, {"min", "max"},
                      f"structural_tables.extremal_edge_tables_R45[{order}]")
        low = _exact_int(record["min"], f"R45[{order}].min")
        high = _exact_int(record["max"], f"R45[{order}].max")
        if low > high:
            raise CheckViolation(f"R45[{order}]: inverted edge bounds")
        observed[order] = (low, high)
    if observed != FROZEN_EDGE_BOUNDS:
        raise CheckViolation(
            f"structural_tables.extremal_edge_tables_R45: observed {observed}, "
            f"expected {FROZEN_EDGE_BOUNDS}"
        )
    return observed


def _derive_r35_windows(root):
    windows = {0: (0, 0)}
    for order in range(1, 14):
        path = Path(root) / f"r55/data/r35_{order}.g6"
        low = high = None
        count = 0
        try:
            stream = path.open("r", encoding="ascii")
        except (OSError, UnicodeError) as exc:
            raise CheckViolation(f"r35_{order}.g6: unreadable ({exc})") from None
        with stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    parsed_order, adjacency = parse_graph6_line(line)
                except (TypeError, ValueError) as exc:
                    raise CheckViolation(
                        f"r35_{order}.g6:{line_number}: invalid graph6 ({exc})"
                    ) from None
                if parsed_order != order:
                    raise CheckViolation(
                        f"r35_{order}.g6:{line_number}: order {parsed_order}, expected {order}"
                    )
                edges = sum(popcount(row) for row in adjacency) // 2
                low = edges if low is None else min(low, edges)
                high = edges if high is None else max(high, edges)
                count += 1
        if count != R35_COUNTS[order - 1]:
            raise CheckViolation(
                f"r35_{order}.g6: observed {count} graphs, "
                f"expected {R35_COUNTS[order - 1]}"
            )
        windows[order] = (low, high)
    return windows


# --------------------------------------------------------- catalog re-sweep

def _sweep_chunk_worker(task):
    source, raw_lines = task
    stats = {}
    for offset, raw in enumerate(raw_lines):
        try:
            line = raw.decode("ascii")
            order, adjacency = parse_graph6_line(line)
            edges = sum(popcount(row) for row in adjacency) // 2
            vector = motif_vector(adjacency)
        except (UnicodeError, TypeError, ValueError) as exc:
            return source, None, f"chunk line {offset + 1}: {exc}"
        if vector["k4"] != 0:
            return source, None, f"class {(order, edges)} contains K4={vector['k4']}"
        cls = (order, edges)
        entry = stats.get(cls)
        values = [vector[key] for key in STREAM_KEYS]
        if entry is None:
            stats[cls] = [1, values[:], values[:]]
        else:
            entry[0] += 1
            for index, value in enumerate(values):
                entry[1][index] = min(entry[1][index], value)
                entry[2][index] = max(entry[2][index], value)
    return source, stats, None


def _catalog_chunks(root, sources):
    for source in sources:
        path = Path(root) / source
        try:
            stream = path.open("rb")
        except OSError as exc:
            raise CheckViolation(f"sweep {source}: unreadable ({exc})") from None
        with stream:
            chunk = []
            for raw in stream:
                if not raw.strip():
                    continue
                chunk.append(raw)
                if len(chunk) == _SWEEP_CHUNK:
                    yield source, tuple(chunk)
                    chunk = []
            if chunk:
                yield source, tuple(chunk)


def _merge_sweep_stats(total_stats, partial):
    for cls, entry in partial.items():
        target = total_stats.get(cls)
        if target is None:
            total_stats[cls] = [entry[0], entry[1][:], entry[2][:]]
            continue
        target[0] += entry[0]
        for index in range(len(STREAM_KEYS)):
            target[1][index] = min(target[1][index], entry[1][index])
            target[2][index] = max(target[2][index], entry[2][index])


def sweep_catalogs(context, root, echo=None, jobs=None):
    """Re-stream all 8,500,211 graphs and require every class window."""
    progress = echo if callable(echo) else (lambda line: None)
    sources = sorted(set(context["owners"].values()))
    tasks = _catalog_chunks(root, sources)
    if jobs is None:
        jobs = min(32, os.cpu_count() or 1)
    _exact_int(jobs, "sweep jobs")
    if jobs < 1:
        raise CheckViolation(f"sweep jobs: expected positive int, got {jobs}")

    total_stats = {}
    completed = 0
    next_notice = _SWEEP_PROGRESS

    def consume(results):
        nonlocal completed, next_notice
        for source, partial, error in results:
            if error is not None:
                raise CheckViolation(f"sweep {source}: {error}")
            partial_count = sum(entry[0] for entry in partial.values())
            completed += partial_count
            _merge_sweep_stats(total_stats, partial)
            while completed >= next_notice:
                progress(f"MIXED SWEEP PROGRESS: {next_notice:,} graphs completed")
                next_notice += _SWEEP_PROGRESS

    if jobs == 1:
        consume(map(_sweep_chunk_worker, tasks))
    else:
        try:
            with multiprocessing.Pool(processes=jobs) as pool:
                consume(pool.imap_unordered(_sweep_chunk_worker, tasks, chunksize=1))
        except CheckViolation:
            raise
        except (OSError, RuntimeError, ValueError) as exc:
            raise CheckViolation(f"catalog sweep worker failure: {exc}") from None

    expected_classes = set(CATALOG_CLASSES)
    if set(total_stats) != expected_classes:
        missing = sorted(expected_classes - set(total_stats))
        extra = sorted(set(total_stats) - expected_classes)
        raise CheckViolation(f"catalog sweep class set: missing={missing}, extra={extra}")
    for cls in CATALOG_CLASSES:
        count, lows, highs = total_stats[cls]
        expected_count = context["class_counts"][cls]
        if count != expected_count:
            raise CheckViolation(
                f"catalog sweep class {cls}: observed {count} graphs, "
                f"recorded {expected_count}"
            )
        for index, key in enumerate(STREAM_KEYS):
            observed = (lows[index], highs[index])
            expected = context["catalog_windows"][cls][key]
            if observed != expected:
                raise CheckViolation(
                    f"catalog sweep class {cls} {key}: observed {observed}, "
                    f"recorded {expected}"
                )
    if completed != EXPECTED_CATALOG_GRAPHS:
        raise CheckViolation(
            f"catalog sweep total: observed {completed}, expected {EXPECTED_CATALOG_GRAPHS}"
        )
    return completed, len(total_stats)


# ----------------------------------------------------------- row algebra

def _degree_bounds(order):
    return max(0, order - 18), min(13, order - 1)


def _q_envelopes(order, r35_windows):
    """All exact q envelopes for one order via histogram dynamic programming.

    For a histogram, write W=sum C(deg,2), L=sum r35_min(deg), and
    H=sum r35_max(deg).  Triangle incidence gives
    ceil(L/3) <= t <= floor(H/3), hence
    q_lo=W-floor(H/3), q_hi=W-ceil(L/3).  The DP tracks H-L capped at two,
    L/H modulo three, and the extrema of 3W-H and 3W-L.  Those statistics
    are sufficient for the floors and for deciding whether the integer
    triangle interval is empty.
    """
    low_degree, high_degree = _degree_bounds(order)
    for degree in range(low_degree, high_degree + 1):
        if degree not in r35_windows:
            raise CheckViolation(
                f"r35_edge_windows: degree {degree} missing for order {order}"
            )

    # (degree_sum, L mod 3, H mod 3, min(2,H-L)) -> (min(3W-H), max(3W-L))
    dynamic = {(0, 0, 0, 0): (0, 0)}
    for _vertex in range(order):
        following = {}
        for (degree_sum, low_mod, high_mod, gap), extrema in dynamic.items():
            for degree in range(low_degree, high_degree + 1):
                local_low, local_high = r35_windows[degree]
                wedge = _c2(degree)
                key = (
                    degree_sum + degree,
                    (low_mod + local_low) % 3,
                    (high_mod + local_high) % 3,
                    min(2, gap + local_high - local_low),
                )
                low_objective = extrema[0] + 3 * wedge - local_high
                high_objective = extrema[1] + 3 * wedge - local_low
                old = following.get(key)
                if old is None:
                    following[key] = (low_objective, high_objective)
                else:
                    following[key] = (
                        min(old[0], low_objective),
                        max(old[1], high_objective),
                    )
        dynamic = following

    envelopes = {}
    for edges in range(_c2(order) + 1):
        lower = upper = None
        target_sum = 2 * edges
        for (degree_sum, low_mod, high_mod, gap), extrema in dynamic.items():
            if degree_sum != target_sum:
                continue
            feasible = (
                gap >= 2
                or (gap == 0 and low_mod == 0)
                or (gap == 1 and low_mod != 1)
            )
            if not feasible:
                continue
            local_lower = (extrema[0] + high_mod) // 3
            low_correction = (-low_mod) % 3
            local_upper = (extrema[1] - low_correction) // 3
            lower = local_lower if lower is None else min(lower, local_lower)
            upper = local_upper if upper is None else max(upper, local_upper)
        if lower is not None:
            envelopes[edges] = (lower, upper)
    return envelopes


def _budget2(degree, edge_bounds):
    dual_order = N_TARGET - 1 - degree
    h2_min = (
        2 * (_c2(dual_order) - edge_bounds[dual_order][1]
             - edge_bounds[degree][1])
        - degree * (N_TARGET - 2 * degree)
    )
    return -h2_min


def _h_interval(degree, x_windows, z_windows):
    coefficient = 3 * (47 - 2 * degree)
    t_products = (
        coefficient * x_windows["t"][0],
        coefficient * x_windows["t"][1],
    )
    return (
        min(t_products)
        + 4 * x_windows["diamond"][0]
        - 6 * x_windows["k3k1"][1]
        - 12 * z_windows["i4"][1],
        max(t_products)
        + 4 * x_windows["diamond"][1]
        - 6 * x_windows["k3k1"][0]
        - 12 * z_windows["i4"][0],
    )


def _affine_interval(coefficient, window):
    values = (coefficient * window[0], coefficient * window[1])
    return min(values), max(values)


def _f_interval(degree, a_value, b_value, windows):
    dual_order = N_TARGET - 1 - degree
    x_edges = E45[degree] - a_value
    z_edges = E45[dual_order] - b_value
    y_edges = _c2(dual_order) - z_edges
    x_class = (degree, x_edges)
    z_class = (dual_order, z_edges)
    exact = (
        1890 * degree
        - 2109 * degree ** 2
        + 135 * degree ** 3
        - 2 * degree ** 4
        + 4124 * x_edges
        - 12 * x_edges ** 2
        - 528 * degree * x_edges
        + 12 * degree ** 2 * x_edges
        + 4 * y_edges ** 2
        - 86 * degree * y_edges
        + 4 * degree ** 2 * y_edges
    )
    terms = (
        (516, "t", x_class),
        (72, "c4", x_class),
        (24, "claw", x_class),
        (24, "p4", x_class),
        (24, "paw", x_class),
        (564 - 24 * degree, "p3", x_class),
        (32, "diamond", x_class),
        (74 + 4 * degree, "pc3", z_class),
        (-12, "k3k1", z_class),
        (-8, "two_k2", z_class),
        (-8, "p3_k1", z_class),
        (-24, "k2_2k1", z_class),
    )
    low = high = exact
    for coefficient, key, cls in terms:
        try:
            interval = windows[cls][key]
        except KeyError:
            raise CheckViolation(f"F row: no {key} window for class {cls}") from None
        term_low, term_high = _affine_interval(coefficient, interval)
        low += term_low
        high += term_high
    return low, high


def _state_cache_key(context, edge_bounds, r35_windows):
    return (
        tuple(sorted(edge_bounds.items())),
        tuple(sorted(r35_windows.items())),
        tuple(
            (cls, context["owners"].get(cls),
             tuple((key, context["catalog_windows"].get(cls, {}).get(key))
                   for key in STREAM_KEYS))
            for cls in sorted(context["catalog_windows"])
        ),
        tuple(
            (cls, context["outer_methods"][cls],
             tuple((key, context["outer_windows"][cls][key]) for key in STREAM_KEYS))
            for cls in sorted(context["outer_windows"])
        ),
    )


def rebuild_states(context, edge_bounds, r35_windows):
    """Rebuild all state coordinates and g/h/F endpoints independently."""
    cache_key = _state_cache_key(context, edge_bounds, r35_windows)
    cached = _STATE_CACHE.get(cache_key)
    if cached is not None:
        return cached

    q_outer = {
        order: _q_envelopes(order, r35_windows)
        for order in STATE_DEGREES
    }

    def q_window(cls):
        if cls in context["catalog_windows"]:
            return context["catalog_windows"][cls]["q"]
        return q_outer[cls[0]].get(cls[1])

    states = []
    for degree in STATE_DEGREES:
        dual_order = N_TARGET - 1 - degree
        degree_low, degree_high = edge_bounds[degree]
        dual_low, dual_high = edge_bounds[dual_order]
        budget = _budget2(degree, edge_bounds)
        for a_value in range(degree_high - degree_low + 1):
            x_edges = degree_high - a_value
            x_class = (degree, x_edges)
            x_q = q_window(x_class)
            if x_q is None:
                continue
            x_windows = context["outer_windows"][x_class]
            for b_value in range(dual_high - dual_low + 1):
                z_edges = dual_high - b_value
                z_class = (dual_order, z_edges)
                z_q = q_window(z_class)
                if z_q is None:
                    continue
                z_windows = context["outer_windows"][z_class]
                deficiency = a_value + b_value
                excess_balance = 2 * deficiency - budget
                dual_base = 3 * (_c3(dual_order) - z_edges * (dual_order - 2))
                x_base = (48 - 3 * degree) * x_edges
                g_low = dual_base - x_base + 3 * z_q[0] - 3 * x_q[1]
                g_high = dual_base - x_base + 3 * z_q[1] - 3 * x_q[0]
                h_low, h_high = _h_interval(degree, x_windows, z_windows)
                f_low, f_high = _f_interval(
                    degree, a_value, b_value, context["outer_windows"]
                )
                x_owner = context["owners"].get(x_class)
                z_owner = context["owners"].get(z_class)
                states.append({
                    "d": degree,
                    "a": a_value,
                    "b": b_value,
                    "deficiency": deficiency,
                    "excess_balance": excess_balance,
                    "g_lo": g_low,
                    "g_hi": g_high,
                    "h_lo": h_low,
                    "h_hi": h_high,
                    "f_lo": f_low,
                    "f_hi": f_high,
                    "x_interval_source": x_owner or "outer_degree_histogram",
                    "y_interval_source": z_owner or "outer_degree_histogram",
                    "x_motif_source": (
                        f"catalog:{x_owner}" if x_owner
                        else context["outer_methods"][x_class]
                    ),
                    "y_motif_source": (
                        f"catalog:{z_owner}" if z_owner
                        else context["outer_methods"][z_class]
                    ),
                })
    states.sort(key=lambda state: (state["d"], state["a"], state["b"]))
    if len(states) != EXPECTED_STATES:
        raise CheckViolation(
            f"state rebuild: derived {len(states)} states, expected {EXPECTED_STATES}"
        )
    result = tuple(states)
    _STATE_CACHE[cache_key] = result
    return result


def _compare_states(recorded, rebuilt):
    if len(recorded) != len(rebuilt):
        raise CheckViolation(
            f"states: recorded {len(recorded)}, rebuilt {len(rebuilt)}"
        )
    fields = (
        "d", "a", "b", "deficiency", "excess_balance", "g_lo", "g_hi",
        "h_lo", "h_hi", "f_lo", "f_hi", "x_interval_source",
        "y_interval_source", "x_motif_source", "y_motif_source",
    )
    for index, (record, truth) in enumerate(zip(recorded, rebuilt)):
        for field in fields:
            if record[field] != truth[field]:
                coordinate = (truth["d"], truth["a"], truth["b"])
                raise CheckViolation(
                    f"states[{index}] {coordinate}.{field}: recorded "
                    f"{record[field]!r}, rebuilt {truth[field]!r}"
                )


# ------------------------------------------------------ exact search evidence

def _objective_value(objective, state):
    if objective == "total_deficiency":
        return Fraction(state["deficiency"])
    if objective == "degree20_count":
        return Fraction(int(state["d"] == 20))
    if objective == "deficiency_ge8_count":
        return Fraction(int(state["deficiency"] >= 8))
    raise CheckViolation(f"objective_id: unregistered objective {objective!r}")


def _verify_certificate(states, record, objective):
    certificate = record["certificate"]
    if certificate is None:
        if record["exact_upper_bound"] is not None:
            raise CheckViolation(f"search {objective}: bound without certificate")
        return None
    values = {
        key: _canonical_fraction(certificate[key], f"search {objective}.{key}")
        for key in COEFFICIENTS
    }
    for key in ("gamma", "delta", "epsilon", "zeta", "eta", "theta"):
        if values[key] < 0:
            raise CheckViolation(f"search {objective}.{key}: sign gate is negative")

    required_alpha = None
    required_index = None
    for index, state in enumerate(states):
        needed = (
            _objective_value(objective, state)
            - values["beta"] * state["excess_balance"]
            - values["gamma"] * state["g_lo"]
            + values["delta"] * state["g_hi"]
            - values["epsilon"] * state["h_lo"]
            + values["zeta"] * state["h_hi"]
            - values["eta"] * state["f_lo"]
            + values["theta"] * state["f_hi"]
        )
        if required_alpha is None or needed > required_alpha:
            required_alpha = needed
            required_index = index
    if values["alpha"] < required_alpha:
        raise CheckViolation(
            f"search {objective}.alpha: {values['alpha']} violates state "
            f"{required_index}, which requires {required_alpha}"
        )
    if values["alpha"] != required_alpha:
        raise CheckViolation(
            f"search {objective}.alpha: {values['alpha']} is not the least "
            f"admissible {required_alpha}"
        )
    bound = VERTICES * values["alpha"]
    recorded_bound = _canonical_fraction(
        record["exact_upper_bound"], f"search {objective}.exact_upper_bound"
    )
    if recorded_bound != bound:
        raise CheckViolation(
            f"search {objective}.exact_upper_bound: recorded {recorded_bound}, "
            f"recomputed {bound}"
        )
    return bound


def _verify_witness(states, record, objective):
    witness = record["primal_witness"]
    if witness is None:
        if record["exact_witness_value"] is not None:
            raise CheckViolation(f"search {objective}: value without primal witness")
        return None
    total = Fraction(0)
    balance = Fraction(0)
    aggregates = {
        key: Fraction(0)
        for key in ("g_lo", "g_hi", "h_lo", "h_hi", "f_lo", "f_hi")
    }
    objective_total = Fraction(0)
    previous = None
    for pair_index, pair in enumerate(witness["weights"]):
        index = pair["state_index"]
        if not 0 <= index < len(states):
            raise CheckViolation(
                f"search {objective}.weights[{pair_index}].state_index: out of range"
            )
        if previous is not None and index <= previous:
            raise CheckViolation(
                f"search {objective}.weights: state indices do not strictly increase"
            )
        previous = index
        weight = _canonical_fraction(
            pair["weight"], f"search {objective}.weights[{pair_index}].weight"
        )
        if weight <= 0:
            raise CheckViolation(
                f"search {objective}.weights[{pair_index}].weight: not positive"
            )
        state = states[index]
        total += weight
        balance += weight * state["excess_balance"]
        for field in aggregates:
            aggregates[field] += weight * state[field]
        objective_total += weight * _objective_value(objective, state)
    if total != VERTICES:
        raise CheckViolation(f"search {objective}.weights: sum {total}, expected 45")
    if balance != 0:
        raise CheckViolation(f"search {objective}.witness balance: {balance}, expected 0")
    for low_field, high_field in (
        ("g_lo", "g_hi"), ("h_lo", "h_hi"), ("f_lo", "f_hi")
    ):
        if aggregates[low_field] > 0:
            raise CheckViolation(
                f"search {objective}.witness {low_field}: "
                f"{aggregates[low_field]} > 0"
            )
        if aggregates[high_field] < 0:
            raise CheckViolation(
                f"search {objective}.witness {high_field}: "
                f"{aggregates[high_field]} < 0"
            )
    recorded = _canonical_fraction(
        record["exact_witness_value"], f"search {objective}.exact_witness_value"
    )
    if recorded != objective_total:
        raise CheckViolation(
            f"search {objective}.exact_witness_value: recorded {recorded}, "
            f"recomputed {objective_total}"
        )
    return objective_total


def _clears(value, limit, inclusive):
    return value <= limit if inclusive else value < limit


def _verify_searches(states, searches, echo=None):
    progress = echo if callable(echo) else (lambda line: None)
    statuses = []
    for record in searches:
        objective = record["objective_id"]
        if objective == UNAVAILABLE_ROUTE:
            if record != UNAVAILABLE_RECORD:
                raise CheckViolation(
                    "search required_local_family: unavailable record carries drift/evidence"
                )
            continue
        bound = _verify_certificate(states, record, objective)
        witness_value = _verify_witness(states, record, objective)
        limit, inclusive = ROUTE_EDGES[objective]
        accepts = bound is not None and _clears(bound, limit, inclusive)
        rejects = (
            witness_value is not None
            and not _clears(witness_value, limit, inclusive)
        )
        if accepts and rejects:
            raise CheckViolation(f"search {objective}: contradictory upper/lower evidence")
        derived = ACCEPTED if accepts else REJECTED if rejects else UNRESOLVED
        if record["route_status"] != derived:
            raise CheckViolation(
                f"search {objective}.route_status: recorded "
                f"{record['route_status']!r}, derived {derived!r}"
            )
        statuses.append(derived)
        progress(
            f"MIXED ROUTE {objective}: bound={bound} witness={witness_value} "
            f"status={derived}"
        )
    if len(statuses) != len(ROUTE_IDS):
        raise CheckViolation(f"searches: derived {len(statuses)} instantiable routes")
    return tuple(statuses)


def _terminal_disposition(statuses):
    if any(status == ACCEPTED for status in statuses):
        return DISPOSITION_ACCEPTED
    if all(status == REJECTED for status in statuses):
        return DISPOSITION_NO_CUT
    return DISPOSITION_UNRESOLVED


# ------------------------------------------------------------------ verdict

def verify_analysis(document, root=None, echo=None, *, check_input_hashes=True,
                    sweep=False, sweep_jobs=None):
    """Re-verify one mixed-tier artifact; return its derived disposition."""
    progress = echo if callable(echo) else (lambda line: None)
    repository = Path(root if root is not None else _DEFAULT_ROOT).resolve()

    context = _check_schema(document)
    progress(
        f"MIXED SCHEMA: v{SCHEMA_VERSION}, {EXPECTED_STATES} states, "
        f"{len(CATALOG_CLASSES)} catalog classes"
    )

    if check_input_hashes:
        _check_hashes(context, repository)
        progress(f"MIXED HASHES: {len(context['inputs'])} inputs byte-identical")
    else:
        progress(
            "MIXED PROVENANCE UNPROVED: --no-input-hashes skipped SHA-256 and "
            f"byte checks for {len(context['inputs'])} recorded inputs; catalog "
            "completeness and input provenance are UNVERIFIED in this run"
        )

    edge_bounds = _load_edge_bounds(repository)
    if check_input_hashes:
        observed_r35 = _derive_r35_windows(repository)
        if observed_r35 != context["r35"]:
            raise CheckViolation(
                f"r35_edge_windows: recorded {context['r35']}, observed {observed_r35}"
            )
        r35_windows = observed_r35
        progress("MIXED R35: 14 edge windows independently re-derived")
    else:
        r35_windows = context["r35"]
        progress("MIXED R35: using recorded windows under unproved-provenance mode")

    if sweep:
        total, classes = sweep_catalogs(
            context, repository, echo=progress, jobs=sweep_jobs
        )
        progress(
            f"MIXED SWEEP: {total:,} graphs across {classes} classes; "
            "every per-class motif window agreed"
        )
    else:
        progress(
            f"MIXED CATALOG TABLES: {EXPECTED_CATALOG_GRAPHS:,} recorded graphs "
            "internally reconciled; byte re-sweep not requested"
        )

    rebuilt = rebuild_states(context, edge_bounds, r35_windows)
    _compare_states(document["states"], rebuilt)
    progress(
        f"MIXED STATES: {len(rebuilt)} balance/g/h/F rows independently rebuilt"
    )

    statuses = _verify_searches(rebuilt, document["searches"], echo=progress)
    disposition = _terminal_disposition(statuses)
    if document["disposition"] != disposition:
        raise CheckViolation(
            f"disposition: recorded {document['disposition']!r}, derived {disposition!r}"
        )
    progress(
        "MIXED SEARCHES: eight-coefficient certificates, least-alpha gates, "
        "and all g/h/F primal interval relations verified"
    )
    return disposition


def _echo(line):
    print(line, flush=True)


def main(argv=None):
    """CLI entry point with m4-checker-compatible exit discipline."""
    parser = argparse.ArgumentParser(
        description=(
            "Disjointly re-verify the mixed Engstrom identity certificate "
            "artifact with exact arithmetic."
        )
    )
    parser.add_argument(
        "--root",
        default=None,
        help=(
            "math repository root for recorded r55/... inputs "
            "(default: the checker's own repository root)"
        ),
    )
    parser.add_argument(
        "--no-input-hashes",
        action="store_true",
        help="loudly skip recorded input byte/SHA-256 verification",
    )
    parser.add_argument(
        "--sweep",
        action="store_true",
        help="opt in to the 8,500,211-graph independent catalog sweep",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="catalog sweep workers (default: up to 32 local CPUs)",
    )
    parser.add_argument("analysis", help="path to engstrom_identity.json")
    args = parser.parse_args(argv)
    if args.jobs is not None and not args.sweep:
        parser.error("--jobs requires --sweep")
    repository = Path(args.root) if args.root is not None else _DEFAULT_ROOT
    try:
        document = load_analysis(args.analysis)
        disposition = verify_analysis(
            document,
            repository,
            echo=_echo,
            check_input_hashes=not args.no_input_hashes,
            sweep=args.sweep,
            sweep_jobs=args.jobs,
        )
    except CheckViolation as exc:
        print(f"MIXED EVIDENCE FAILED: {exc}", file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print(f"MIXED EVIDENCE FAILED: analysis: {exc}", file=sys.stderr)
        return 1
    print(f"MIXED EVIDENCE VERIFIED: {disposition}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
